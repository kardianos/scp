/* chain.c — v103 lane c1: voids in a chain medium (../PREREG.md).
 *
 * Cells on a periodic lattice hold Q_i field quanta (q = 2 normally; a void
 * is a cell missing quanta).  Links carry n_l bonds; a bond pairs one quantum
 * at each end, so deg_i = sum_l n_l <= Q_i.  The uniform measure over valid
 * (bonds, quanta) states is the Rokhsar-Kivelson ground state of the quantum
 * medium that re-forms chains; diagonal ground-state correlations are its
 * averages.
 *
 * Candidates: a  exclusion per link end (n_l <= 1; off: n_l <= q)
 *             b  chain conservation (deg_i == Q_i; off: deg_i <= Q_i)
 *             c  bipartite fabric (ring/square/cubic; off: kagome/pyrochlore)
 *
 *   ./chain mode=run lat=cubic L=12 rules=bc D=2 seeds=8 sweeps=... warm=... out=f
 *   ./chain mode=t0 lat=square L1=2 L2=4 rules=ab D=2 seeds=32 sweeps=...
 *   ./chain mode=lat lat=pyro L=5            (lattice + Green's-function gates)
 */
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define QQ 2 /* quanta per normal cell: chains */

/* ---------------- lattice ---------------- */
static int N, NL, Z, DIM;
static double *pos;   /* N x 3 cartesian */
static double box[3][3], boxinv[3][3];
static int *lnk;      /* NL x 2 */
static double *ldisp; /* NL x 3: displacement from lnk[0] to lnk[1] */
static int *inc, *nbr; /* N x Z */
static int *color;    /* 2-colouring (c on) or all 0 */
static int bipartite;
static int *upbond;   /* NL: 1 if in the initial full chain configuration */
static int *cellid, *basisid; /* for initial configurations */

static void inv3(double m[3][3], double o[3][3])
{
    double d = m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) +
               m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
    o[0][0] = (m[1][1] * m[2][2] - m[1][2] * m[2][1]) / d;
    o[0][1] = (m[0][2] * m[2][1] - m[0][1] * m[2][2]) / d;
    o[0][2] = (m[0][1] * m[1][2] - m[0][2] * m[1][1]) / d;
    o[1][0] = (m[1][2] * m[2][0] - m[1][0] * m[2][2]) / d;
    o[1][1] = (m[0][0] * m[2][2] - m[0][2] * m[2][0]) / d;
    o[1][2] = (m[0][2] * m[1][0] - m[0][0] * m[1][2]) / d;
    o[2][0] = (m[1][0] * m[2][1] - m[1][1] * m[2][0]) / d;
    o[2][1] = (m[0][1] * m[2][0] - m[0][0] * m[2][1]) / d;
    o[2][2] = (m[0][0] * m[1][1] - m[0][1] * m[1][0]) / d;
}
/* box rows are the box vectors; frac = cart * boxinv (row vector convention) */
static void tofrac(const double *c, double *f)
{
    for (int k = 0; k < 3; k++) f[k] = c[0] * boxinv[0][k] + c[1] * boxinv[1][k] + c[2] * boxinv[2][k];
}
static void tocart(const double *f, double *c)
{
    for (int k = 0; k < 3; k++) c[k] = f[0] * box[0][k] + f[1] * box[1][k] + f[2] * box[2][k];
}
/* minimum-image distance between cells i and j */
static double mindist(int i, int j)
{
    double d[3], f[3], best = 1e300;
    for (int k = 0; k < 3; k++) d[k] = pos[3 * j + k] - pos[3 * i + k];
    tofrac(d, f);
    for (int k = 0; k < DIM; k++) f[k] -= floor(f[k] + 0.5);
    int r0 = DIM > 0 ? 1 : 0, r1 = DIM > 1 ? 1 : 0, r2 = DIM > 2 ? 1 : 0;
    for (int a = -r0; a <= r0; a++)
        for (int b = -r1; b <= r1; b++)
            for (int c = -r2; c <= r2; c++) {
                double g[3] = {f[0] + a, f[1] + b, f[2] + c}, e[3];
                tocart(g, e);
                double s = e[0] * e[0] + e[1] * e[1] + e[2] * e[2];
                if (s < best) best = s;
            }
    return sqrt(best);
}

static void add_sites(int n)
{
    N = n;
    pos = calloc(3 * (size_t)N, sizeof(double));
    cellid = calloc(N, sizeof(int));
    basisid = calloc(N, sizeof(int));
}

/* find every nearest-neighbour link (distance 1) over images; multi-links allowed */
static void find_links(void)
{
    int cap = 64, nl = 0;
    lnk = malloc(sizeof(int) * 2 * cap);
    ldisp = malloc(sizeof(double) * 3 * cap);
    int R = 2;
    for (int i = 0; i < N; i++)
        for (int j = i; j < N; j++) {
            double d[3], f[3];
            for (int k = 0; k < 3; k++) d[k] = pos[3 * j + k] - pos[3 * i + k];
            tofrac(d, f);
            int r0 = DIM > 0 ? R : 0, r1 = DIM > 1 ? R : 0, r2 = DIM > 2 ? R : 0;
            for (int a = -r0; a <= r0; a++)
                for (int b = -r1; b <= r1; b++)
                    for (int c = -r2; c <= r2; c++) {
                        double g[3] = {f[0] + a, f[1] + b, f[2] + c}, e[3];
                        tocart(g, e);
                        double s = sqrt(e[0] * e[0] + e[1] * e[1] + e[2] * e[2]);
                        if (fabs(s - 1.0) > 1e-6) continue;
                        if (i == j) {
                            /* a self-link: count it once per +/- pair by keeping one orientation */
                            fprintf(stderr, "self-link at cell %d\n", i);
                            continue;
                        }
                        if (nl == cap) {
                            cap *= 2;
                            lnk = realloc(lnk, sizeof(int) * 2 * cap);
                            ldisp = realloc(ldisp, sizeof(double) * 3 * cap);
                        }
                        lnk[2 * nl] = i;
                        lnk[2 * nl + 1] = j;
                        for (int k = 0; k < 3; k++) ldisp[3 * nl + k] = e[k];
                        nl++;
                    }
        }
    NL = nl;
    int *cnt = calloc((size_t)(N > 0 ? N : 1), sizeof(int));
    for (int l = 0; l < NL; l++) {
        cnt[lnk[2 * l]]++;
        cnt[lnk[2 * l + 1]]++;
    }
    Z = cnt[0];
    for (int i = 0; i < N; i++)
        if (cnt[i] != Z) {
            fprintf(stderr, "degree mismatch: cell %d has %d, cell 0 has %d\n", i, cnt[i], Z);
            exit(1);
        }
    inc = malloc(sizeof(int) * N * Z);
    nbr = malloc(sizeof(int) * N * Z);
    memset(cnt, 0, sizeof(int) * N);
    for (int l = 0; l < NL; l++) {
        int i = lnk[2 * l], j = lnk[2 * l + 1];
        inc[i * Z + cnt[i]] = l;
        nbr[i * Z + cnt[i]++] = j;
        inc[j * Z + cnt[j]] = l;
        nbr[j * Z + cnt[j]++] = i;
    }
    free(cnt);
    /* 2-colouring */
    color = malloc(sizeof(int) * N);
    for (int i = 0; i < N; i++) color[i] = -1;
    int *q = malloc(sizeof(int) * N), qh = 0, qt = 0, ok = 1;
    color[0] = 0;
    q[qt++] = 0;
    while (qh < qt) {
        int i = q[qh++];
        for (int k = 0; k < Z; k++) {
            int j = nbr[i * Z + k];
            if (color[j] < 0) {
                color[j] = 1 - color[i];
                q[qt++] = j;
            } else if (color[j] == color[i])
                ok = 0;
        }
    }
    free(q);
    bipartite = ok;
    if (!ok)
        for (int i = 0; i < N; i++) color[i] = 0;
}

static int same_disp(const double *a, const double *b)
{
    return fabs(a[0] - b[0]) < 1e-6 && fabs(a[1] - b[1]) < 1e-6 && fabs(a[2] - b[2]) < 1e-6;
}

static char latname[16];
static int LL[3];

static void build(const char *lat, int L1, int L2, int L3)
{
    strncpy(latname, lat, 15);
    memset(box, 0, sizeof(box));
    LL[0] = L1;
    LL[1] = L2;
    LL[2] = L3;
    if (!strcmp(lat, "ring")) {
        DIM = 1;
        add_sites(L1);
        for (int i = 0; i < N; i++) pos[3 * i] = i;
        box[0][0] = L1;
        box[1][1] = box[2][2] = 1;
    } else if (!strcmp(lat, "square")) {
        DIM = 2;
        add_sites(L1 * L2);
        for (int y = 0; y < L2; y++)
            for (int x = 0; x < L1; x++) {
                int i = y * L1 + x;
                pos[3 * i] = x;
                pos[3 * i + 1] = y;
            }
        box[0][0] = L1;
        box[1][1] = L2;
        box[2][2] = 1;
    } else if (!strcmp(lat, "cubic")) {
        DIM = 3;
        add_sites(L1 * L2 * L3);
        for (int z = 0; z < L3; z++)
            for (int y = 0; y < L2; y++)
                for (int x = 0; x < L1; x++) {
                    int i = (z * L2 + y) * L1 + x;
                    pos[3 * i] = x;
                    pos[3 * i + 1] = y;
                    pos[3 * i + 2] = z;
                }
        box[0][0] = L1;
        box[1][1] = L2;
        box[2][2] = L3;
    } else if (!strcmp(lat, "kagome")) {
        DIM = 2;
        add_sites(3 * L1 * L2);
        double a1[2] = {2, 0}, a2[2] = {1, sqrt(3.0)};
        double bs[3][2] = {{0, 0}, {1, 0}, {0.5, sqrt(3.0) / 2}};
        for (int y = 0; y < L2; y++)
            for (int x = 0; x < L1; x++)
                for (int s = 0; s < 3; s++) {
                    int i = 3 * (y * L1 + x) + s;
                    pos[3 * i] = x * a1[0] + y * a2[0] + bs[s][0];
                    pos[3 * i + 1] = x * a1[1] + y * a2[1] + bs[s][1];
                    cellid[i] = y * L1 + x;
                    basisid[i] = s;
                }
        box[0][0] = L1 * a1[0];
        box[0][1] = L1 * a1[1];
        box[1][0] = L2 * a2[0];
        box[1][1] = L2 * a2[1];
        box[2][2] = 1;
    } else if (!strcmp(lat, "pyro")) {
        DIM = 3;
        int Lc = L1;
        add_sites(16 * Lc * Lc * Lc);
        double a = 2 * sqrt(2.0);
        double f[4][3] = {{0, 0, 0}, {0, .5, .5}, {.5, 0, .5}, {.5, .5, 0}};
        double b[4][3] = {{0, 0, 0}, {0, .25, .25}, {.25, 0, .25}, {.25, .25, 0}};
        int i = 0;
        for (int z = 0; z < Lc; z++)
            for (int y = 0; y < Lc; y++)
                for (int x = 0; x < Lc; x++)
                    for (int fi = 0; fi < 4; fi++)
                        for (int bi = 0; bi < 4; bi++) {
                            pos[3 * i] = a * (x + f[fi][0] + b[bi][0]);
                            pos[3 * i + 1] = a * (y + f[fi][1] + b[bi][1]);
                            pos[3 * i + 2] = a * (z + f[fi][2] + b[bi][2]);
                            cellid[i] = ((z * Lc + y) * Lc + x) * 4 + fi;
                            basisid[i] = bi;
                            i++;
                        }
        box[0][0] = box[1][1] = box[2][2] = Lc * a;
    } else {
        fprintf(stderr, "unknown lattice %s\n", lat);
        exit(1);
    }
    inv3(box, boxinv);
    find_links();
    /* the initial full chain configuration (every cell deg = 2) */
    upbond = calloc(NL, sizeof(int));
    int *dg = calloc(N, sizeof(int));
    for (int l = 0; l < NL; l++) {
        int i = lnk[2 * l], j = lnk[2 * l + 1];
        int take = 0;
        if (!strcmp(lat, "ring")) take = 1;
        else if (!strcmp(lat, "square") || !strcmp(lat, "cubic")) take = fabs(fabs(ldisp[3 * l]) - 1) < 1e-6;
        else if (!strcmp(lat, "kagome")) {
            /* up triangle: same unit cell, displacement = basis difference */
            if (cellid[i] == cellid[j]) {
                double d[3] = {pos[3 * j] - pos[3 * i], pos[3 * j + 1] - pos[3 * i + 1], 0};
                take = same_disp(d, ldisp + 3 * l);
            }
        } else if (!strcmp(lat, "pyro")) {
            /* up tetrahedron (same cube and fcc point), 4-cycle 0-1-2-3-0 */
            if (cellid[i] == cellid[j]) {
                double d[3] = {pos[3 * j] - pos[3 * i], pos[3 * j + 1] - pos[3 * i + 1], pos[3 * j + 2] - pos[3 * i + 2]};
                int bi = basisid[i], bj = basisid[j];
                int cyc = (abs(bi - bj) == 1) || (bi + bj == 3 && abs(bi - bj) == 3);
                take = cyc && same_disp(d, ldisp + 3 * l);
            }
        }
        if (take && dg[i] < QQ && dg[j] < QQ) {
            upbond[l] = 1;
            dg[i]++;
            dg[j]++;
        }
    }
    for (int i = 0; i < N; i++)
        if (dg[i] != QQ) {
            fprintf(stderr, "initial chain configuration failed at cell %d (deg %d)\n", i, dg[i]);
            exit(1);
        }
    free(dg);
}

/* ---------------- Green's functions ---------------- */
static void lap(const double *x, double *y, double m2)
{
    for (int i = 0; i < N; i++) {
        double s = (Z + m2) * x[i];
        for (int k = 0; k < Z; k++) s -= x[nbr[i * Z + k]];
        y[i] = s;
    }
}
/* solve (L + m2) g = rhs by CG; returns relative residual */
static double cg(double m2, const double *rhs, double *g)
{
    double *r = malloc(sizeof(double) * N), *p = malloc(sizeof(double) * N), *Ap = malloc(sizeof(double) * N);
    memset(g, 0, sizeof(double) * N);
    memcpy(r, rhs, sizeof(double) * N);
    memcpy(p, r, sizeof(double) * N);
    double rr = 0, b2 = 0;
    for (int i = 0; i < N; i++) rr += r[i] * r[i], b2 += rhs[i] * rhs[i];
    for (int it = 0; it < 20 * N && rr > 1e-28 * b2; it++) {
        lap(p, Ap, m2);
        double pAp = 0;
        for (int i = 0; i < N; i++) pAp += p[i] * Ap[i];
        double al = rr / pAp, rn = 0;
        for (int i = 0; i < N; i++) {
            g[i] += al * p[i];
            r[i] -= al * Ap[i];
            rn += r[i] * r[i];
        }
        for (int i = 0; i < N; i++) p[i] = r[i] + rn / rr * p[i];
        rr = rn;
    }
    lap(g, Ap, m2);
    double res = 0;
    for (int i = 0; i < N; i++) res += (Ap[i] - rhs[i]) * (Ap[i] - rhs[i]);
    free(r);
    free(p);
    free(Ap);
    return sqrt(res / b2);
}

/* ---------------- shells ---------------- */
static int NS;
static double *shr;       /* shell radius */
static uint16_t *shell;   /* N x N */
static double *ideal[2];  /* per shell: unordered pairs (same, opp colour); same-cell pairs in shell 0 */
static void build_shells(void)
{
    double *d = malloc(sizeof(double) * (size_t)N * N);
    int cap = 256;
    shr = malloc(sizeof(double) * cap);
    NS = 0;
    shr[NS++] = 0;
    for (int i = 0; i < N; i++)
        for (int j = i + 1; j < N; j++) {
            double r = mindist(i, j);
            d[(size_t)i * N + j] = r;
        }
    /* distinct radii */
    for (int i = 0; i < N; i++)
        for (int j = i + 1; j < N; j++) {
            double r = d[(size_t)i * N + j];
            int s;
            for (s = 0; s < NS; s++)
                if (fabs(shr[s] - r) < 1e-6) break;
            if (s == NS) {
                if (NS == cap) cap *= 2, shr = realloc(shr, sizeof(double) * cap);
                shr[NS++] = r;
            }
        }
    /* sort radii */
    for (int a = 1; a < NS; a++)
        for (int b = a; b > 0 && shr[b] < shr[b - 1]; b--) {
            double t = shr[b];
            shr[b] = shr[b - 1];
            shr[b - 1] = t;
        }
    shell = malloc(sizeof(uint16_t) * (size_t)N * N);
    ideal[0] = calloc(NS, sizeof(double));
    ideal[1] = calloc(NS, sizeof(double));
    for (int i = 0; i < N; i++) {
        shell[(size_t)i * N + i] = 0;
        ideal[0][0] += 1;
        for (int j = i + 1; j < N; j++) {
            double r = d[(size_t)i * N + j];
            int lo = 0, hi = NS - 1;
            while (lo < hi) {
                int mid = (lo + hi) / 2;
                if (shr[mid] < r - 1e-6) lo = mid + 1;
                else hi = mid;
            }
            shell[(size_t)i * N + j] = shell[(size_t)j * N + i] = (uint16_t)lo;
            ideal[color[i] != color[j]][lo] += 1;
        }
    }
    free(d);
}

/* Green's functions averaged over each shell, per colour type, on the m grid */
#define NM 13
static const double MGRID[NM] = {0, 0.02, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.7, 1.0, 1.5, 2.0};
static double *Gsh[NM][2];
static double gres_max = 0;
static void build_green(void)
{
    /* representatives: one per (cellid-free) translation class.  Translation
       invariance holds within a basis/colour class; use every cell with a
       distinct basisid among the first unit cell's, or cell 0 for Bravais
       lattices.  Simplest exact route: average over a set of sources that
       covers each basis class once, weighted by class size. */
    int nb = 1;
    if (!strcmp(latname, "kagome")) nb = 3;
    if (!strcmp(latname, "pyro")) nb = 16; /* the 16 sites of the first cube */
    double *rhs = malloc(sizeof(double) * N), *g = malloc(sizeof(double) * N);
    for (int m = 0; m < NM; m++)
        for (int t = 0; t < 2; t++) Gsh[m][t] = calloc(NS, sizeof(double));
    double *w[2] = {calloc(NS, sizeof(double)), calloc(NS, sizeof(double))};
    for (int s = 0; s < nb; s++) {
        int src = s;
        for (int m = 0; m < NM; m++) {
            double m2 = MGRID[m] * MGRID[m];
            for (int i = 0; i < N; i++) rhs[i] = (i == src) - (m == 0 ? 1.0 / N : 0.0);
            double res = cg(m2, rhs, g);
            if (res > gres_max) gres_max = res;
            if (m == 0) {
                double mean = 0;
                for (int i = 0; i < N; i++) mean += g[i];
                mean /= N;
                for (int i = 0; i < N; i++) g[i] -= mean;
            }
            for (int j = 0; j < N; j++) {
                int sh = shell[(size_t)src * N + j], t = color[src] != color[j];
                Gsh[m][t][sh] += g[j];
                if (m == 0) w[t][sh] += 1;
            }
        }
    }
    for (int m = 0; m < NM; m++)
        for (int t = 0; t < 2; t++)
            for (int s = 0; s < NS; s++)
                if (w[t][s] > 0) Gsh[m][t][s] /= w[t][s];
    free(rhs);
    free(g);
    free(w[0]);
    free(w[1]);
}

/* ---------------- RNG ---------------- */
static uint64_t sm(uint64_t *x)
{
    uint64_t z = (*x += 0x9e3779b97f4a7c15ULL);
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}
static uint64_t rs[4];
static inline uint64_t rotl(uint64_t x, int k) { return (x << k) | (x >> (64 - k)); }
static inline uint64_t xnext(void)
{
    uint64_t r = rotl(rs[1] * 5, 7) * 9, t = rs[1] << 17;
    rs[2] ^= rs[0];
    rs[3] ^= rs[1];
    rs[1] ^= rs[2];
    rs[0] ^= rs[3];
    rs[2] ^= t;
    rs[3] = rotl(rs[3], 45);
    return r;
}
static inline double urand(void) { return ((xnext() >> 11) + 0.5) * (1.0 / 9007199254740992.0); }
static inline int irand(int n) { return (int)(urand() * n); }
static void rseed(uint64_t s)
{
    uint64_t x = s * 0x2545F4914F6CDD1DULL + 777;
    for (int i = 0; i < 4; i++) rs[i] = sm(&x);
}

/* ---------------- state and moves ---------------- */
static int ra, rb; /* rules a, b (c is the lattice) */
static int CAP;    /* link capacity */
static int D;      /* missing quanta */
static int *Q, *deg, *nl;
static int dunits[16]; /* cells of the D deficit units */
static long acc[3], att[3];

static int deficit(int i) { return QQ - Q[i]; }

static void init_state(void)
{
    Q = malloc(sizeof(int) * N);
    deg = calloc(N, sizeof(int));
    nl = calloc(NL, sizeof(int));
    for (int i = 0; i < N; i++) Q[i] = QQ;
    if (rb) {
        for (int l = 0; l < NL; l++)
            if (upbond[l]) {
                nl[l] = 1;
                deg[lnk[2 * l]]++;
                deg[lnk[2 * l + 1]]++;
            }
        /* remove D/2 distinct bonds: each makes two deficits */
        int k = 0;
        while (k < D / 2) {
            int l = irand(NL);
            if (nl[l] < 1) continue;
            int i = lnk[2 * l], j = lnk[2 * l + 1];
            if (Q[i] < QQ || Q[j] < QQ) continue;
            nl[l]--;
            deg[i]--;
            deg[j]--;
            Q[i]--;
            Q[j]--;
            dunits[2 * k] = i;
            dunits[2 * k + 1] = j;
            k++;
        }
    } else {
        int k = 0;
        while (k < D) {
            int i = irand(N);
            if (Q[i] < QQ) continue; /* start distinct; they may merge later */
            Q[i]--;
            dunits[k++] = i;
        }
    }
}

/* b on: void hop.  Returns 1 if accepted. */
static int move_hop(void)
{
    att[0]++;
    int u = irand(D), v = dunits[u];
    int kk = irand(Z), l1 = inc[v * Z + kk], k = nbr[v * Z + kk];
    if (nl[l1] >= CAP) return 0;
    int dk = deg[k]; /* units at k excluding the one to be added */
    if (dk == 0) return 0;
    int pick = irand(dk), l2 = -1, lcell = -1;
    for (int t = 0; t < Z; t++) {
        int l = inc[k * Z + t];
        int nu = nl[l];
        if (pick < nu) {
            l2 = l;
            lcell = nbr[k * Z + t];
            break;
        }
        pick -= nu;
    }
    if (l2 == l1) return 0; /* null move */
    /* acceptance for the uniform measure: proposal counts in each direction */
    int defv = deficit(v), defl_after = deficit(lcell) + 1 - (lcell == v ? 1 : 0);
    double pf = (double)defv * nl[l2];             /* x 1/(D z dk) */
    double pr = (double)defl_after * (nl[l1] + 1); /* reverse: unit at l, link l2, remove on l1 */
    if (lcell == v) {
        /* the hop returns to v through a parallel link: deficits unchanged */
        defl_after = defv;
        pr = (double)defv * (nl[l1] + 1);
    }
    if (pr < pf && urand() * pf >= pr) return 0;
    nl[l1]++;
    nl[l2]--;
    if (lcell != v) {
        Q[v]++;
        deg[v]++;
        Q[lcell]--;
        deg[lcell]--;
        dunits[u] = lcell;
    }
    acc[0]++;
    return 1;
}
/* b off: link flip */
static int move_flip(void)
{
    att[1]++;
    int l = irand(NL), dlt = (xnext() >> 63) ? 1 : -1;
    int n2 = nl[l] + dlt;
    if (n2 < 0 || n2 > CAP) return 0;
    int i = lnk[2 * l], j = lnk[2 * l + 1];
    if (deg[i] + dlt > Q[i] || deg[j] + dlt > Q[j]) return 0;
    nl[l] = n2;
    deg[i] += dlt;
    deg[j] += dlt;
    acc[1]++;
    return 1;
}
/* b off: a quantum moves from neighbour k into void cell v */
static int move_qhop(void)
{
    att[2]++;
    int u = irand(D), v = dunits[u];
    int kk = irand(Z), k = nbr[v * Z + kk];
    if (Q[k] < 1 || deg[k] > Q[k] - 1 || Q[v] + 1 > QQ) return 0;
    double pf = deficit(v), pr = deficit(k) + 1;
    if (pr < pf && urand() * pf >= pr) return 0;
    Q[k]--;
    Q[v]++;
    dunits[u] = k;
    acc[2]++;
    return 1;
}
static inline void step(void)
{
    if (rb) move_hop();
    else if (xnext() >> 63) move_flip();
    else move_qhop();
}
static long ledger_viol = 0;
static void check_ledger(void)
{
    long sq = 0;
    int *dd = calloc(N, sizeof(int));
    for (int l = 0; l < NL; l++) {
        if (nl[l] < 0 || nl[l] > CAP) ledger_viol++;
        dd[lnk[2 * l]] += nl[l];
        dd[lnk[2 * l + 1]] += nl[l];
    }
    for (int i = 0; i < N; i++) {
        sq += Q[i];
        if (dd[i] != deg[i]) ledger_viol++;
        if (rb ? deg[i] != Q[i] : deg[i] > Q[i]) ledger_viol++;
        if (Q[i] < 0 || Q[i] > QQ) ledger_viol++;
    }
    if (sq != (long)QQ * N - D) ledger_viol++;
    int *cnt = calloc(N, sizeof(int));
    for (int u = 0; u < D; u++) cnt[dunits[u]]++;
    for (int i = 0; i < N; i++)
        if (cnt[i] != QQ - Q[i]) ledger_viol++;
    free(cnt);
    free(dd);
}

/* ---------------- exact enumeration (T0) ---------------- */
static double *exact_bins; /* multiset {i<=j} -> weight */

/* ---------------- exact transition matrix (T0 gate 1b) ---------------- */
typedef struct {
    uint64_t a, b, c;
} Key;
static int dbmode = 0;
static Key *skeys;
static long nst = 0, capst = 0;
static long *htab;
static long hcap = 0;
static uint64_t hkey(Key k)
{
    uint64_t h = k.a * 0x9E3779B97F4A7C15ULL ^ (k.b + 0x632BE59BD9B4E019ULL) * 0xBF58476D1CE4E5B9ULL ^ (k.c + 0x1234567ULL) * 0x94D049BB133111EBULL;
    return h ^ (h >> 29);
}
static Key mkkey(const int *n, const int *q)
{
    Key k = {0, 0, 0};
    int h = NL / 2;
    for (int l = h - 1; l >= 0; l--) k.a = k.a * (CAP + 1) + n[l];
    for (int l = NL - 1; l >= h; l--) k.b = k.b * (CAP + 1) + n[l];
    for (int i = N - 1; i >= 0; i--) k.c = k.c * 3 + q[i];
    return k;
}
static void unkey(Key k, int *n, int *q)
{
    int h = NL / 2;
    for (int l = 0; l < h; l++) n[l] = k.a % (CAP + 1), k.a /= (CAP + 1);
    for (int l = h; l < NL; l++) n[l] = k.b % (CAP + 1), k.b /= (CAP + 1);
    for (int i = 0; i < N; i++) q[i] = k.c % 3, k.c /= 3;
}
static long hfind(Key k, int insert)
{
    if (insert && 2 * (nst + 1) > hcap) {
        long nc = hcap ? 2 * hcap : 1 << 16;
        long *nt = malloc(sizeof(long) * nc);
        for (long i = 0; i < nc; i++) nt[i] = -1;
        for (long s2 = 0; s2 < nst; s2++) {
            uint64_t h = hkey(skeys[s2]) & (nc - 1);
            while (nt[h] >= 0) h = (h + 1) & (nc - 1);
            nt[h] = s2;
        }
        free(htab);
        htab = nt;
        hcap = nc;
    }
    uint64_t h = hkey(k) & (hcap - 1);
    while (htab[h] >= 0) {
        Key o = skeys[htab[h]];
        if (o.a == k.a && o.b == k.b && o.c == k.c) return htab[h];
        h = (h + 1) & (hcap - 1);
    }
    if (!insert) return -1;
    if (nst == capst) capst = capst ? 2 * capst : 1 << 16, skeys = realloc(skeys, sizeof(Key) * capst);
    skeys[nst] = k;
    htab[h] = nst;
    return nst++;
}
static int *ecur, *edeg, *erem;
static void db_record(const int *def)
{
    int *q = malloc(sizeof(int) * N);
    for (int i = 0; i < N; i++) q[i] = QQ - def[i];
    hfind(mkkey(ecur, q), 1);
    free(q);
}

static int bin_of(int i, int j)
{
    if (i > j) {
        int t = i;
        i = j;
        j = t;
    }
    return i * N + j;
}
static void enum_leaf(void)
{
    if (rb) {
        int s = 0, c[2], nc = 0;
        for (int i = 0; i < N; i++) {
            int d = QQ - edeg[i];
            if (d < 0) return;
            s += d;
            for (int t = 0; t < d && nc < 2; t++) c[nc++] = i;
        }
        if (s != 2) return;
        if (dbmode) {
            int *def = calloc((size_t)(N > 0 ? N : 1), sizeof(int));
            for (int i = 0; i < N; i++) def[i] = QQ - edeg[i];
            db_record(def);
            free(def);
            return;
        }
        exact_bins[bin_of(c[0], c[1])] += 1;
    } else {
        for (int i = 0; i < N; i++) {
            int si = QQ - edeg[i];
            if (si < 0) return;
        }
        if (dbmode) {
            int *def = calloc((size_t)(N > 0 ? N : 1), sizeof(int));
            for (int i = 0; i < N; i++) {
                int si = QQ - edeg[i];
                if (si >= 2) def[i] = 2, db_record(def), def[i] = 0;
                if (si >= 1)
                    for (int j = i + 1; j < N; j++)
                        if (QQ - edeg[j] >= 1) def[i] = def[j] = 1, db_record(def), def[i] = def[j] = 0;
            }
            free(def);
            return;
        }
        for (int i = 0; i < N; i++) {
            int si = QQ - edeg[i];
            if (si >= 2) exact_bins[bin_of(i, i)] += 1;
            if (si >= 1)
                for (int j = i + 1; j < N; j++)
                    if (QQ - edeg[j] >= 1) exact_bins[bin_of(i, j)] += 1;
        }
    }
}
static void enum_rec(int l)
{
    if (l == NL) {
        enum_leaf();
        return;
    }
    int i = lnk[2 * l], j = lnk[2 * l + 1];
    erem[i] -= CAP;
    erem[j] -= CAP;
    for (int v = 0; v <= CAP; v++) {
        if (edeg[i] + v > QQ || edeg[j] + v > QQ) break;
        edeg[i] += v;
        edeg[j] += v;
        ecur[l] = v;
        int prune = 0;
        if (rb) {
            /* minimal achievable total deficit must not exceed 2 */
            int md = 0;
            for (int c = 0; c < N && !prune; c++) {
                int best = edeg[c] + erem[c];
                if (best > QQ) best = QQ;
                md += QQ - best;
                if (md > 2) prune = 1;
            }
        }
        if (!prune) enum_rec(l + 1);
        edeg[i] -= v;
        edeg[j] -= v;
    }
    erem[i] += CAP;
    erem[j] += CAP;
}

/* ---------------- main ---------------- */
static const char *arg(int argc, char **argv, const char *k, const char *def)
{
    size_t n = strlen(k);
    for (int i = 1; i < argc; i++)
        if (!strncmp(argv[i], k, n) && argv[i][n] == '=') return argv[i] + n + 1;
    return def;
}

int main(int argc, char **argv)
{
    const char *mode = arg(argc, argv, "mode", "run");
    const char *lat = arg(argc, argv, "lat", "cubic");
    int L = atoi(arg(argc, argv, "L", "4"));
    int L1 = atoi(arg(argc, argv, "L1", "0")), L2 = atoi(arg(argc, argv, "L2", "0")), L3 = atoi(arg(argc, argv, "L3", "0"));
    if (!L1) L1 = L;
    if (!L2) L2 = L;
    if (!L3) L3 = L;
    const char *rules = arg(argc, argv, "rules", "none");
    ra = strchr(rules, 'a') != NULL;
    rb = strchr(rules, 'b') != NULL;
    int rc_req = strchr(rules, 'c') != NULL;
    CAP = ra ? 1 : QQ;
    D = atoi(arg(argc, argv, "D", "2"));
    int seeds = atoi(arg(argc, argv, "seeds", "8"));
    long sweeps = atol(arg(argc, argv, "sweeps", "1000"));
    long warm = atol(arg(argc, argv, "warm", "100"));
    uint64_t seed0 = strtoull(arg(argc, argv, "seed", "1"), NULL, 10);

    build(lat, L1, L2, L3);
    printf("# chain lat=%s L=%d,%d,%d N=%d NL=%d Z=%d dim=%d bipartite=%d rules=%s (a=%d b=%d c=%d) CAP=%d q=%d D=%d\n", lat, L1, L2,
           L3, N, NL, Z, DIM, bipartite, rules, ra, rb, rc_req, CAP, QQ, D);
    if (rc_req != bipartite) {
        printf("# WARNING: rule c=%d but lattice bipartite=%d\n", rc_req, bipartite);
    }
    printf("LATTICE links=%d expect=%d selfcheck=%s\n", NL, N * Z / 2, NL == N * Z / 2 ? "ok" : "FAIL");

    if (!strcmp(mode, "lat")) {
        build_shells();
        build_green();
        printf("GREEN NS=%d residual_max=%.3e\n", NS, gres_max);
        for (int s = 0; s < NS && s < 12; s++)
            printf("  shell %d r=%.4f ideal_same=%.0f ideal_opp=%.0f G0=%.6f/%.6f\n", s, shr[s], ideal[0][s], ideal[1][s], Gsh[0][0][s],
                   Gsh[0][1][s]);
        return 0;
    }


    if (!strcmp(mode, "db")) {
        /* every valid state, then the sampler's exact transition probabilities into each:
           uniform is stationary iff every column sums to 1 */
        dbmode = 1;
        exact_bins = calloc(1, sizeof(double));
        ecur = calloc(NL, sizeof(int));
        edeg = calloc(N, sizeof(int));
        erem = calloc(N, sizeof(int));
        for (int l = 0; l < NL; l++) erem[lnk[2 * l]] += CAP, erem[lnk[2 * l + 1]] += CAP;
        enum_rec(0);
        double *in = calloc(nst, sizeof(double));
        int *n0 = malloc(sizeof(int) * NL), *q0 = malloc(sizeof(int) * N), *n1 = malloc(sizeof(int) * NL), *q1 = malloc(sizeof(int) * N);
        int *dg = malloc(sizeof(int) * N), du[16];
        long missing = 0;
        for (long st = 0; st < nst; st++) {
            unkey(skeys[st], n0, q0);
            memset(dg, 0, sizeof(int) * N);
            for (int l = 0; l < NL; l++) dg[lnk[2 * l]] += n0[l], dg[lnk[2 * l + 1]] += n0[l];
            int nd = 0;
            for (int i = 0; i < N; i++)
                for (int t = 0; t < QQ - q0[i]; t++) du[nd++] = i;
            double stay = 0;
#define GO(prob)                                                     \
    do {                                                             \
        long id = hfind(mkkey(n1, q1), 0);                           \
        if (id < 0) missing++;                                       \
        else in[id] += (prob);                                       \
    } while (0)
            if (rb) {
                for (int u = 0; u < D; u++) {
                    int v = du[u];
                    for (int kk = 0; kk < Z; kk++) {
                        double p1 = 1.0 / D / Z;
                        int l1 = inc[v * Z + kk], k = nbr[v * Z + kk];
                        if (n0[l1] >= CAP || dg[k] == 0) { stay += p1; continue; }
                        int dk = dg[k];
                        for (int t = 0; t < Z; t++) {
                            int l2 = inc[k * Z + t], lc = nbr[k * Z + t];
                            if (n0[l2] == 0) continue;
                            double p2 = p1 * n0[l2] / dk;
                            if (l2 == l1) { stay += p2; continue; }
                            int defv = QQ - q0[v];
                            double pf = (double)defv * n0[l2], pr;
                            if (lc == v) pr = (double)defv * (n0[l1] + 1);
                            else pr = (double)(QQ - q0[lc] + 1) * (n0[l1] + 1);
                            double A = pr < pf ? pr / pf : 1;
                            memcpy(n1, n0, sizeof(int) * NL);
                            memcpy(q1, q0, sizeof(int) * N);
                            n1[l1]++;
                            n1[l2]--;
                            if (lc != v) q1[v]++, q1[lc]--;
                            GO(p2 * A);
                            stay += p2 * (1 - A);
                        }
                    }
                }
            } else {
                for (int l = 0; l < NL; l++)
                    for (int dl = -1; dl <= 1; dl += 2) {
                        double p = 0.5 / NL * 0.5;
                        int i = lnk[2 * l], j = lnk[2 * l + 1], m2 = n0[l] + dl;
                        if (m2 < 0 || m2 > CAP || dg[i] + dl > q0[i] || dg[j] + dl > q0[j]) { stay += p; continue; }
                        memcpy(n1, n0, sizeof(int) * NL);
                        memcpy(q1, q0, sizeof(int) * N);
                        n1[l] = m2;
                        GO(p);
                    }
                for (int u = 0; u < D; u++) {
                    int v = du[u];
                    for (int kk = 0; kk < Z; kk++) {
                        double p = 0.5 / D / Z;
                        int k = nbr[v * Z + kk];
                        if (q0[k] < 1 || dg[k] > q0[k] - 1 || q0[v] + 1 > QQ) { stay += p; continue; }
                        double pf = QQ - q0[v], pr = QQ - q0[k] + 1, A = pr < pf ? pr / pf : 1;
                        memcpy(n1, n0, sizeof(int) * NL);
                        memcpy(q1, q0, sizeof(int) * N);
                        q1[k]--;
                        q1[v]++;
                        GO(p * A);
                        stay += p * (1 - A);
                    }
                }
            }
            in[st] += stay;
        }
        double dmax = 0;
        for (long st = 0; st < nst; st++)
            if (fabs(in[st] - 1) > dmax) dmax = fabs(in[st] - 1);
        printf("DB states=%ld max|colsum-1|=%.3e missing_targets=%ld  %s\n", nst, dmax, missing, dmax < 1e-12 && missing == 0 ? "PASS" : "FAIL");
        return 0;
    }

    if (!strcmp(mode, "t0")) {
        if (D != 2) {
            fprintf(stderr, "t0 needs D=2\n");
            return 1;
        }
        long nb = (long)N * N;
        exact_bins = calloc(nb, sizeof(double));
        ecur = calloc(NL, sizeof(int));
        edeg = calloc(N, sizeof(int));
        erem = calloc(N, sizeof(int));
        for (int l = 0; l < NL; l++) erem[lnk[2 * l]] += CAP, erem[lnk[2 * l + 1]] += CAP;
        enum_rec(0);
        double tot = 0;
        for (long b = 0; b < nb; b++) tot += exact_bins[b];
        printf("EXACT states=%.0f\n", tot);
        /* Monte Carlo, many seeds */
        double *f = calloc((size_t)seeds * nb, sizeof(double));
        for (int s = 0; s < seeds; s++) {
            rseed(seed0 + 1000003ULL * s);
            init_state();
            for (long w = 0; w < warm * N; w++) step();
            for (long t = 0; t < sweeps; t++) {
                for (int m = 0; m < N; m++) step();
                f[(size_t)s * nb + bin_of(dunits[0], dunits[1])] += 1;
            }
            check_ledger();
            for (long b = 0; b < nb; b++) f[(size_t)s * nb + b] /= sweeps;
            free(Q);
            free(deg);
            free(nl);
        }
        double chi2 = 0, zmax = 0;
        int dof = 0, bad0 = 0;
        for (long b = 0; b < nb; b++) {
            double p = exact_bins[b] / tot, m = 0, v = 0;
            for (int s = 0; s < seeds; s++) m += f[(size_t)s * nb + b];
            m /= seeds;
            for (int s = 0; s < seeds; s++) v += (f[(size_t)s * nb + b] - m) * (f[(size_t)s * nb + b] - m);
            v /= (seeds - 1);
            double se = sqrt(v / seeds);
            if (p == 0) {
                if (m > 0) bad0++;
                continue;
            }
            if (se == 0) continue;
            double z = (m - p) / se;
            chi2 += z * z;
            dof++;
            if (fabs(z) > zmax) zmax = fabs(z);
        }
        int pass = dof > 0 && chi2 / dof <= 1.5 && zmax <= 4.5 && bad0 == 0 && ledger_viol == 0;
        printf("T0 bins=%d chi2/dof=%.3f zmax=%.2f forbidden_hits=%d ledger_viol=%ld acc_hop=%.3f acc_flip=%.3f acc_qhop=%.3f  %s\n", dof,
               dof ? chi2 / dof : 0, zmax, bad0, ledger_viol, att[0] ? (double)acc[0] / att[0] : 0, att[1] ? (double)acc[1] / att[1] : 0,
               att[2] ? (double)acc[2] / att[2] : 0, pass ? "PASS" : "FAIL");
        return 0;
    }

    /* production run */
    build_shells();
    build_green();
    printf("GREEN NS=%d residual_max=%.3e\n", NS, gres_max);
    int npair = D * (D - 1) / 2;
    double *cnt = calloc((size_t)seeds * 2 * NS, sizeof(double)); /* [seed][type][shell] */
    /* anemia: unbonded quanta per cell vs shell to nearest void */
    double *an = calloc((size_t)seeds * NS, sizeof(double)), *anw = calloc((size_t)seeds * NS, sizeof(double));
    for (int s = 0; s < seeds; s++) {
        rseed(seed0 + 1000003ULL * s);
        init_state();
        for (long w = 0; w < warm * N; w++) step();
        for (long t = 0; t < sweeps; t++) {
            for (int m = 0; m < N; m++) step();
            for (int u = 0; u < D; u++)
                for (int w = u + 1; w < D; w++) {
                    int i = dunits[u], j = dunits[w];
                    cnt[((size_t)s * 2 + (color[i] != color[j])) * NS + shell[(size_t)i * N + j]] += 1;
                }
            if (!rb && (t % 16 == 0)) {
                for (int i = 0; i < N; i++) {
                    int best = 1 << 30;
                    for (int u = 0; u < D; u++) {
                        int sh = shell[(size_t)i * N + dunits[u]];
                        if (sh < best) best = sh;
                    }
                    an[(size_t)s * NS + best] += (Q[i] - deg[i]);
                    anw[(size_t)s * NS + best] += 1;
                }
            }
            if (t % 1024 == 0) check_ledger();
        }
        check_ledger();
        free(Q);
        free(deg);
        free(nl);
        fprintf(stderr, "seed %d done\n", s);
    }
    printf("RUN sweeps=%ld warm=%ld seeds=%d pairs_per_sample=%d\n", sweeps, warm, seeds, npair);
    printf("ACC hop=%.4f flip=%.4f qhop=%.4f\n", att[0] ? (double)acc[0] / att[0] : 0, att[1] ? (double)acc[1] / att[1] : 0,
           att[2] ? (double)acc[2] / att[2] : 0);
    printf("LEDGER violations=%ld\n", ledger_viol);
    printf("MGRID");
    for (int m = 0; m < NM; m++) printf(" %g", MGRID[m]);
    printf("\n");
    for (int sh = 0; sh < NS; sh++) {
        for (int t = 0; t < 2; t++) {
            if (ideal[t][sh] == 0) continue;
            printf("SHELL %d r=%.6f type=%s ideal=%.0f counts", sh, shr[sh], t ? "opp" : "same", ideal[t][sh]);
            for (int s = 0; s < seeds; s++) printf(" %.0f", cnt[((size_t)s * 2 + t) * NS + sh]);
            printf(" G");
            for (int m = 0; m < NM; m++) printf(" %.9e", Gsh[m][t][sh]);
            printf("\n");
        }
    }
    if (!rb) {
        for (int sh = 0; sh < NS; sh++) {
            double tw = 0;
            for (int s = 0; s < seeds; s++) tw += anw[(size_t)s * NS + sh];
            if (tw == 0) continue;
            printf("ANEMIA %d r=%.6f", sh, shr[sh]);
            for (int s = 0; s < seeds; s++) printf(" %.6f", anw[(size_t)s * NS + sh] > 0 ? an[(size_t)s * NS + sh] / anw[(size_t)s * NS + sh] : -1);
            printf("\n");
        }
    }
    printf("CHECK ledger_violations=%ld green_residual=%.3e\n", ledger_viol, gres_max);
    return 0;
}
