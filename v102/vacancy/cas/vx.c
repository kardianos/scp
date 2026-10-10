/* vx.c — v102 vacancy step (../PREREG.md): exact diagonalization of the
 * lock in one charge sector.
 *
 *   model=lo  the leading-order effective lock of ../../perm/cas/lo.py
 *             (reading V): exchange -(g^2/Delta_type) A^+A on every link,
 *             the T-floor hop t on mobile sites, R12 (singlet core) on
 *             (1,1) links as lo.py's 'score' (P11_CORE = 1.0 rate units);
 *   model=ao  the microscopic S-law of ../../perm/cas/graphm.py: field
 *             energies E_f(N), link oscillators at Omega, conversion
 *             g sum P (c^+ b b + h.c.), hop t (T2: refused into a cell at
 *             its cap; the T-floor: refused if it would empty a cell);
 *             basis = closure of the seeds under H, sum c <= cmax.
 *
 * M = 2 (3+2, P = eps) or M = 4 (W, P = eps x eps, modes (s,tau) =
 * uu, ud, du, dd).  geom=ring or open (R6 end cells, cap 1).  vac=none,
 * static (cell vsite holds one quantum less), mobile (one vacancy, any
 * bulk cell; needs t > 0 to move).  Energies: rate units internally,
 * printed also in J0 = 2 g^2/Delta.
 *
 * Lanczos with full reorthogonalization (LAPACK dstev on the
 * tridiagonal).  pen=1 adds lambda (S^2 + I^2), lambda = 5 J0, so the
 * lowest states are S = 0 (and I = 0).
 *
 * Build: gcc -O3 -march=native -fopenmp -o vx vx.c -llapack -lm
 * Usage: vx key=value ...   (see parse() for keys)                     */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <omp.h>

extern void dgemv_(const char *trans, const int *m, const int *n, const double *alpha, const double *a, const int *lda,
                   const double *x, const int *incx, const double *beta, double *y, const int *incy);
extern void dstev_(const char *jobz, const int *n, double *d, double *e,
                   double *z, const int *ldz, double *work, int *info);

/* ---------------- constants (../../options/cas/graph.py) ---------------- */
static double EF[8], OMEGA = 2.9, DELTA, G, J0;
static void constants(double gr)
{
    const double W1 = 1.65, QD = 1.2, CAPL = 2.5, A0 = 1.15;
    double hb = A0 / (2 * M_PI), a = QD * hb / CAPL;
    EF[0] = 0;
    for (int N = 1; N < 8; N++) EF[N] = (-1 + sqrt(1 + 4 * a * N * W1)) / (2 * a);
    DELTA = OMEGA - 2 * (EF[2] - EF[1]);
    G = gr * DELTA;
    J0 = 2 * G * G / DELTA;
}

/* ---------------- parameters ---------------- */
static int M = 4, n = 6, ring = 1, aomodel = 0, vac = 0, vsite = -1, pen = 0;
static int sz2 = 0, iz2 = 0, kwant = 4, nvec = 1, mmax = 300, cmax = -1, seed = 12345;
static double tJ0 = 0, gr = 0.1, lam_J0 = 5.0, tol = 1e-10;
static int obs_prof = 0, obs_corr = 0, obs_cs = 0, nrestart = 8;

static double P[4][4];
static int ms2[4], mt2[4];       /* (2 s_z, 2 tau_z) of each mode */

/* ---------------- local cell spaces ---------------- */
#define MAXD 16
typedef struct {
    int d, cap;
    int occ[MAXD][4];
    int N[MAXD], s2[MAXD], t2[MAXD];
    int key2idx[81];
} Loc;

static int okey(const int *o) { int k = 0, p = 1; for (int a = 0; a < M; a++) { k += o[a] * p; p *= 3; } return k; }

static void gen_N(Loc *L, int N, int *o, int a, int left)
{
    if (a == M - 1) {
        o[a] = left;
        int d = L->d++;
        int s = 0, t = 0;
        for (int b = 0; b < M; b++) { L->occ[d][b] = o[b]; s += o[b] * ms2[b]; t += o[b] * mt2[b]; }
        L->N[d] = N; L->s2[d] = s; L->t2[d] = t;
        L->key2idx[okey(o)] = d;
        return;
    }
    for (int x = left; x >= 0; x--) { o[a] = x; gen_N(L, N, o, a + 1, left - x); }
}

/* the cell space holding N in [Nlo, Nhi], N descending (lo.py order) */
static void loc_build(Loc *L, int Nlo, int Nhi, int cap)
{
    int o[4] = {0, 0, 0, 0};
    L->d = 0; L->cap = cap;
    for (int k = 0; k < 81; k++) L->key2idx[k] = -1;
    for (int N = Nhi; N >= Nlo; N--) gen_N(L, N, o, 0, N);
}

/* b^+_r b_m on state k: new index or -1; amplitude */
static int Eop(const Loc *L, int k, int r, int m, double *amp)
{
    if (L->occ[k][m] == 0) return -1;
    int o[4] = {0, 0, 0, 0};
    for (int a = 0; a < M; a++) o[a] = L->occ[k][a];
    o[m]--; o[r]++;
    if (o[r] > 2) return -1;
    int k2 = L->key2idx[okey(o)];
    if (k2 < 0) return -1;
    *amp = (r != m) ? sqrt((double)L->occ[k][m] * o[r]) : (double)L->occ[k][m];
    return k2;
}

/* b_m (dir -1) or b^+_m (dir +1) */
static int bop(const Loc *L, int k, int m, int dir, double *amp)
{
    int o[4] = {0, 0, 0, 0};
    for (int a = 0; a < M; a++) o[a] = L->occ[k][a];
    if (dir < 0) { if (o[m] == 0) return -1; *amp = sqrt((double)o[m]); o[m]--; }
    else { if (o[m] >= 2) return -1; *amp = sqrt((double)o[m] + 1); o[m]++; }
    return L->key2idx[okey(o)];
}

/* ---------------- the chain ---------------- */
#define NMAX 20
static Loc locs[6];              /* distinct cell spaces */
static Loc *site[NMAX];
static int nlink, li[NMAX], lj[NMAX];
static int nsl;                  /* digits: n cells (+ nlink links in AO) */
static int64_t stride[2 * NMAX];
static int radix[2 * NMAX];
static int clmax = 2;            /* AO: per-link c radix - 1 */
static int kind[NMAX];           /* 0 full/cap-2 bulk, 1 end (cap 1), 2 static vacancy, 3 mobile bulk */

/* LO link tables, column-sparse over the pair index a*dB + b */
typedef struct { int dA, dB; int *cp, *ri; double *v; } Two;
static Two two[NMAX];

static double dtype(int ni, int nj)
{
    if (ni == 2 && nj == 2) return DELTA;
    return OMEGA - (EF[2] - EF[1]) - EF[1];
}

static void two_build(Two *T, const Loc *A, const Loc *B, int Ahop, int Bhop, double t)
{
    int dA = A->d, dB = B->d, dd = dA * dB;
    double *X = calloc((size_t)dd * dd, sizeof(double));
    double P2 = 0; for (int a = 0; a < M; a++) for (int b = 0; b < M; b++) P2 += P[a][b] * P[a][b];
    for (int a = 0; a < dA; a++)
        for (int b = 0; b < dB; b++) {
            int ni = A->N[a], nj = B->N[b], col = a * dB + b;
            double c = (ni == 1 && nj == 1) ? 1.0 / P2 : -G * G / dtype(ni, nj);   /* R12: 'score' */
            for (int r = 0; r < M; r++) for (int m = 0; m < M; m++)
                for (int s = 0; s < M; s++) for (int q = 0; q < M; q++) {
                    double w = P[r][s] * P[m][q];
                    if (w == 0) continue;
                    double ea, eb;
                    int a2 = Eop(A, a, r, m, &ea); if (a2 < 0) continue;
                    int b2 = Eop(B, b, s, q, &eb); if (b2 < 0) continue;
                    X[(size_t)(a2 * dB + b2) * dd + col] += c * w * ea * eb;
                }
            if (t != 0 && Ahop && Bhop)
                for (int m = 0; m < M; m++) {
                    double ea, eb;
                    /* quantum A -> B: allowed N_A 2 -> 1, N_B 1 -> 2 */
                    if (ni == 2 && nj == 1) {
                        int a2 = bop(A, a, m, -1, &ea), b2 = a2 >= 0 ? bop(B, b, m, +1, &eb) : -1;
                        if (a2 >= 0 && b2 >= 0) X[(size_t)(a2 * dB + b2) * dd + col] += t * ea * eb;
                    }
                    if (ni == 1 && nj == 2) {
                        int a2 = bop(A, a, m, +1, &ea), b2 = a2 >= 0 ? bop(B, b, m, -1, &eb) : -1;
                        if (a2 >= 0 && b2 >= 0) X[(size_t)(a2 * dB + b2) * dd + col] += t * ea * eb;
                    }
                }
        }
    T->dA = dA; T->dB = dB;
    T->cp = calloc(dd + 1, sizeof(int));
    int nnz = 0;
    for (int col = 0; col < dd; col++) for (int row = 0; row < dd; row++) if (fabs(X[(size_t)row * dd + col]) > 1e-15) nnz++;
    T->ri = malloc((nnz + 1) * sizeof(int)); T->v = malloc((nnz + 1) * sizeof(double));
    nnz = 0;
    for (int col = 0; col < dd; col++) {
        T->cp[col] = nnz;
        for (int row = 0; row < dd; row++) {
            double x = X[(size_t)row * dd + col];
            if (fabs(x) > 1e-15) { T->ri[nnz] = row; T->v[nnz] = x; nnz++; }
        }
    }
    T->cp[dd] = nnz;
    free(X);
}

/* ---------------- basis ---------------- */
static int64_t *codes;
static long D;
static double Eref = 0;

static inline void decode(int64_t c, int *dg) { for (int i = 0; i < nsl; i++) dg[i] = (int)((c / stride[i]) % radix[i]); }

static inline long findc(int64_t c)
{
    long lo = 0, hi = D - 1;
    while (lo <= hi) {
        long mid = (lo + hi) >> 1;
        if (codes[mid] == c) return mid;
        if (codes[mid] < c) lo = mid + 1; else hi = mid - 1;
    }
    return -1;
}

static int sector_ok(const int *dg)
{
    int s = 0, t = 0, nh = 0;
    for (int i = 0; i < n; i++) {
        const Loc *L = site[i];
        s += L->s2[dg[i]]; t += L->t2[dg[i]];
        if (kind[i] == 3 && L->N[dg[i]] == 1) nh++;
    }
    if (s != sz2) return 0;
    if (M == 4 && t != iz2) return 0;
    if (!aomodel && vac == 2 && nh != 1) return 0;
    return 1;
}

/* LO: every product state in the sector */
static void basis_lo(void)
{
    int64_t tot = 1;
    for (int i = 0; i < n; i++) tot *= radix[i];
    long cap = 1 << 20; codes = malloc(cap * sizeof(int64_t)); D = 0;
    int dg[2 * NMAX];
    for (int64_t c = 0; c < tot; c++) {
        decode(c, dg);
        if (!sector_ok(dg)) continue;
        if (D == cap) { cap *= 2; codes = realloc(codes, cap * sizeof(int64_t)); }
        codes[D++] = c;
    }
}

/* AO: the closure of the seeds under H, by breadth-first search */
typedef struct { int64_t *k; long cap, cnt; } HSet;
static inline uint64_t hmix(uint64_t x) { x ^= x >> 33; x *= 0xff51afd7ed558ccdULL; x ^= x >> 33; x *= 0xc4ceb9fe1a85ec53ULL; x ^= x >> 33; return x; }
static void hs_init(HSet *h, long cap) { h->cap = cap; h->cnt = 0; h->k = malloc(cap * sizeof(int64_t)); for (long i = 0; i < cap; i++) h->k[i] = -1; }
static int hs_add(HSet *h, int64_t key);
static void hs_grow(HSet *h)
{
    HSet g; hs_init(&g, h->cap * 2);
    for (long i = 0; i < h->cap; i++) if (h->k[i] >= 0) hs_add(&g, h->k[i]);
    free(h->k); *h = g;
}
static int hs_add(HSet *h, int64_t key)
{
    if (2 * (h->cnt + 1) > h->cap) hs_grow(h);
    uint64_t m = h->cap - 1, p = hmix((uint64_t)key) & m;
    while (h->k[p] >= 0) { if (h->k[p] == key) return 0; p = (p + 1) & m; }
    h->k[p] = key; h->cnt++;
    return 1;
}

/* AO off-diagonal moves from a digit vector: callback(newcode, amp) */
typedef void (*MoveFn)(int64_t c2, double amp, void *ctx);
static double t_rate = 0;
static int sumc_of(const int *dg) { int s = 0; for (int l = 0; l < nlink; l++) s += dg[n + l]; return s; }

static void ao_moves(int64_t c, const int *dg, MoveFn f, void *ctx)
{
    int sc = sumc_of(dg);
    for (int l = 0; l < nlink; l++) {
        int i = li[l], j = lj[l], cl = dg[n + l];
        const Loc *A = site[i], *B = site[j];
        int a = dg[i], b = dg[j];
        for (int mu = 0; mu < M; mu++) for (int nu = 0; nu < M; nu++) {
            double p = P[mu][nu];
            if (p == 0) continue;
            double ea, eb;
            /* c^+ b_i mu b_j nu */
            if ((cmax < 0 || sc + 1 <= cmax) && cl + 1 <= clmax) {
                int a2 = bop(A, a, mu, -1, &ea), b2 = a2 >= 0 ? bop(B, b, nu, -1, &eb) : -1;
                if (a2 >= 0 && b2 >= 0) {
                    int64_t c2 = c + (int64_t)(a2 - a) * stride[i] + (int64_t)(b2 - b) * stride[j] + stride[n + l];
                    f(c2, G * p * ea * eb * sqrt((double)cl + 1), ctx);
                }
            }
            /* b^+ b^+ c */
            if (cl > 0) {
                int a2 = bop(A, a, mu, +1, &ea), b2 = a2 >= 0 ? bop(B, b, nu, +1, &eb) : -1;
                if (a2 >= 0 && b2 >= 0 && A->N[a2] <= A->cap && B->N[b2] <= B->cap) {
                    int64_t c2 = c + (int64_t)(a2 - a) * stride[i] + (int64_t)(b2 - b) * stride[j] - stride[n + l];
                    f(c2, G * p * ea * eb * sqrt((double)cl), ctx);
                }
            }
        }
        if (t_rate != 0 && kind[i] != 1 && kind[j] != 1)
            for (int dir = 0; dir < 2; dir++) {
                int src = dir ? j : i, dst = dir ? i : j;
                int as = dg[src], ad = dg[dst];
                const Loc *S = site[src], *Dd = site[dst];
                if (S->N[as] < 2) continue;               /* T-floor: may not empty a cell */
                if (Dd->N[ad] >= Dd->cap) continue;       /* T2: refused into a cell at its cap */
                for (int m = 0; m < M; m++) {
                    double e1, e2;
                    int s2 = bop(S, as, m, -1, &e1), d2 = s2 >= 0 ? bop(Dd, ad, m, +1, &e2) : -1;
                    if (s2 < 0 || d2 < 0) continue;
                    int64_t c2 = c + (int64_t)(s2 - as) * stride[src] + (int64_t)(d2 - ad) * stride[dst];
                    f(c2, t_rate * e1 * e2, ctx);
                }
            }
    }
}

typedef struct { HSet *hs; int64_t *q; long *qn, *qcap; } BfsCtx;
static void bfs_cb(int64_t c2, double amp, void *ctx)
{
    (void)amp;
    BfsCtx *b = ctx;
    if (hs_add(b->hs, c2)) {
        if (*b->qn == *b->qcap) { *b->qcap *= 2; b->q = realloc(b->q, *b->qcap * sizeof(int64_t)); }
        b->q[(*b->qn)++] = c2;
    }
}
static int cmp64(const void *a, const void *b) { int64_t x = *(const int64_t *)a, y = *(const int64_t *)b; return x < y ? -1 : x > y; }

static void basis_ao(void)
{
    /* seeds: every field-only product state in the sector with each cell at its target N */
    HSet hs; hs_init(&hs, 1 << 16);
    long qcap = 1 << 16, qn = 0, qh = 0;
    int64_t *q = malloc(qcap * sizeof(int64_t));
    int dg[2 * NMAX];
    int nvpos = 0, vpos[NMAX];
    if (vac == 2) { for (int i = 0; i < n; i++) if (kind[i] == 3) vpos[nvpos++] = i; }
    else { vpos[0] = vac == 1 ? vsite : -1; nvpos = 1; }
    for (int vp = 0; vp < nvpos; vp++) {
        int tgt[NMAX];
        for (int i = 0; i < n; i++) tgt[i] = site[i]->cap - (i == vpos[vp] ? 1 : 0);
        /* enumerate the cells' states with N = tgt by an odometer over each cell's list */
        int idx[NMAX], lst[NMAX][MAXD], nl[NMAX];
        for (int i = 0; i < n; i++) {
            nl[i] = 0;
            for (int k = 0; k < site[i]->d; k++) if (site[i]->N[k] == tgt[i]) lst[i][nl[i]++] = k;
            idx[i] = 0;
        }
        for (;;) {
            for (int i = 0; i < n; i++) dg[i] = lst[i][idx[i]];
            for (int l = 0; l < nlink; l++) dg[n + l] = 0;
            int s = 0, t = 0;
            for (int i = 0; i < n; i++) { s += site[i]->s2[dg[i]]; t += site[i]->t2[dg[i]]; }
            if (s == sz2 && (M == 2 || t == iz2)) {
                int64_t c = 0;
                for (int i = 0; i < nsl; i++) c += (int64_t)dg[i] * stride[i];
                if (hs_add(&hs, c)) { if (qn == qcap) { qcap *= 2; q = realloc(q, qcap * sizeof(int64_t)); } q[qn++] = c; }
            }
            int i = n - 1;
            while (i >= 0 && ++idx[i] == nl[i]) { idx[i] = 0; i--; }
            if (i < 0) break;
        }
    }
    BfsCtx ctx = {&hs, q, &qn, &qcap};
    while (qh < qn) {
        int64_t c = ctx.q[qh++];
        decode(c, dg);
        ao_moves(c, dg, bfs_cb, &ctx);
    }
    q = ctx.q;
    D = hs.cnt;
    codes = malloc(D * sizeof(int64_t));
    long k = 0;
    for (long i = 0; i < hs.cap; i++) if (hs.k[i] >= 0) codes[k++] = hs.k[i];
    qsort(codes, D, sizeof(int64_t), cmp64);
    free(hs.k); free(q);
}

/* ---------------- H x and S^2 x ---------------- */
static double *diagv;            /* diagonal of H (AO: field + space energies; LO: from the tables) */

typedef struct { double acc; const double *x; } GCtx;
static void gather_cb(int64_t c2, double amp, void *ctx)
{
    GCtx *g = ctx;
    long k = findc(c2);
    if (k >= 0) g->acc += amp * g->x[k];
}

/* spin (which=0) or isospin (which=1) raising / lowering on a cell */
static int rl_modes(int which, int raise, int *from, int *to)
{
    if (M == 2) { if (which) return 0; from[0] = raise ? 1 : 0; to[0] = raise ? 0 : 1; return 1; }
    if (which == 0) { if (raise) { from[0] = 2; to[0] = 0; from[1] = 3; to[1] = 1; } else { from[0] = 0; to[0] = 2; from[1] = 1; to[1] = 3; } }
    else { if (raise) { from[0] = 1; to[0] = 0; from[1] = 3; to[1] = 2; } else { from[0] = 0; to[0] = 1; from[1] = 2; to[1] = 3; } }
    return 2;
}


/* ---- S^2 (I^2) through the neighbouring sector: S^2 = S^- S^+ + Sz^2 + Sz ---- */
static int64_t *auxc[2];         /* sector (sz2+2, iz2) [0] and (sz2, iz2+2) [1] */
static long auxD[2];
static double *phi[2];

static inline long findc_in(const int64_t *cs, long nn, int64_t c)
{
    long lo = 0, hi = nn - 1;
    while (lo <= hi) {
        long mid = (lo + hi) >> 1;
        if (cs[mid] == c) return mid;
        if (cs[mid] < c) lo = mid + 1; else hi = mid - 1;
    }
    return -1;
}

/* phi = S^+ x (which 0) or I^+ x (which 1), gathered over the aux basis by lowering */
static void raise_into(int which, const double *x)
{
    int fl[2], tl[2];
    int nlw = rl_modes(which, 0, fl, tl);
    #pragma omp parallel for schedule(dynamic, 256)
    for (long u = 0; u < auxD[which]; u++) {
        int dg[2 * NMAX];
        int64_t c = auxc[which][u];
        decode(c, dg);
        double acc = 0;
        for (int i = 0; i < n; i++)
            for (int w = 0; w < nlw; w++) {
                double e;
                int b = Eop(site[i], dg[i], tl[w], fl[w], &e);
                if (b < 0) continue;
                long q = findc(c + (int64_t)(b - dg[i]) * stride[i]);
                if (q >= 0) acc += e * x[q];
            }
        phi[which][u] = acc;
    }
}

/* (S^- phi)_s for one main-basis row, gathered by raising */
static double lower_row(int which, const int *dg, int64_t c)
{
    int fr[2], to[2];
    int nr = rl_modes(which, 1, fr, to);
    double acc = 0;
    for (int i = 0; i < n; i++)
        for (int u = 0; u < nr; u++) {
            double e;
            int a = Eop(site[i], dg[i], to[u], fr[u], &e);
            if (a < 0) continue;
            long q = findc_in(auxc[which], auxD[which], c + (int64_t)(a - dg[i]) * stride[i]);
            if (q >= 0) acc += e * phi[which][q];
        }
    return acc;
}

static double zz_of(int which, const int *dg)
{
    int z2 = 0;
    for (int i = 0; i < n; i++) z2 += which ? site[i]->t2[dg[i]] : site[i]->s2[dg[i]];
    double z = z2 / 2.0;
    return z * z + z;
}

/* <x|S^2|x> = |S^+ x|^2 + (Sz^2 + Sz) */
static double expect_sq(int which, const double *x)
{
    if (which == 1 && M == 2) return 0;
    raise_into(which, x);
    double s = 0;
    for (long u = 0; u < auxD[which]; u++) s += phi[which][u] * phi[which][u];
    int dg[2 * NMAX];
    decode(codes[0], dg);
    return s + zz_of(which, dg);
}

static double lam = 0;
static void Hx(const double *x, double *y)
{
    if (lam != 0) { raise_into(0, x); if (M == 4) raise_into(1, x); }
    #pragma omp parallel for schedule(dynamic, 256)
    for (long s = 0; s < D; s++) {
        int dg[2 * NMAX];
        int64_t c = codes[s];
        decode(c, dg);
        double acc = diagv[s] * x[s];
        if (!aomodel) {
            for (int l = 0; l < nlink; l++) {
                int i = li[l], j = lj[l];
                const Two *T = &two[l];
                int a = dg[i], b = dg[j], col = a * T->dB + b;
                for (int k = T->cp[col]; k < T->cp[col + 1]; k++) {
                    int row = T->ri[k];
                    if (row == col) continue;
                    int a2 = row / T->dB, b2 = row % T->dB;
                    int64_t c2 = c + (int64_t)(a2 - a) * stride[i] + (int64_t)(b2 - b) * stride[j];
                    long q = findc(c2);
                    if (q >= 0) acc += T->v[k] * x[q];
                }
            }
        } else {
            GCtx g = {0, x};
            ao_moves(c, dg, gather_cb, &g);
            acc += g.acc;
        }
        if (lam != 0) {
            acc += lam * (lower_row(0, dg, c) + zz_of(0, dg) * x[s]);
            if (M == 4) acc += lam * (lower_row(1, dg, c) + zz_of(1, dg) * x[s]);
        }
        y[s] = acc;
    }
}

static void build_diag(void)
{
    diagv = malloc(D * sizeof(double));
    #pragma omp parallel for
    for (long s = 0; s < D; s++) {
        int dg[2 * NMAX];
        decode(codes[s], dg);
        double e = 0;
        if (!aomodel) {
            for (int l = 0; l < nlink; l++) {
                const Two *T = &two[l];
                int col = dg[li[l]] * T->dB + dg[lj[l]];
                for (int k = T->cp[col]; k < T->cp[col + 1]; k++) if (T->ri[k] == col) e += T->v[k];
            }
        } else {
            for (int i = 0; i < n; i++) e += EF[site[i]->N[dg[i]]];
            for (int l = 0; l < nlink; l++) e += OMEGA * dg[n + l];
            e -= Eref;
        }
        diagv[s] = e;
    }
}

/* ---------------- Lanczos ---------------- */
static double dot(const double *a, const double *b)
{
    double s = 0;
    #pragma omp parallel for reduction(+:s)
    for (long i = 0; i < D; i++) s += a[i] * b[i];
    return s;
}

/* one Lanczos pass from start vector v0 (unit norm): lowest nwant Ritz pairs
 * into ev/res/ritz (ritz may be NULL); returns steps; *convp = converged */
static int lanczos_pass(const double *v0, int nwant, double *ev, double *resid_out, double **ritz, int *convp)
{
    int m = mmax < D ? mmax : (int)D;
    double **V = calloc(m + 2, sizeof(double *));
    double *Vb = malloc((size_t)(m + 1) * D * sizeof(double));   /* contiguous: column j = V[j] */
    double *al = calloc(m + 1, sizeof(double)), *be = calloc(m + 1, sizeof(double));
    double *w = malloc(D * sizeof(double)), *hh = malloc((m + 1) * sizeof(double));
    V[0] = Vb;
    memcpy(V[0], v0, D * sizeof(double));
    int k, conv = 0, used = 0;
    double *d = malloc((m + 1) * sizeof(double)), *e = malloc((m + 1) * sizeof(double));
    double *Z = NULL, *work = malloc(2 * (m + 1) * sizeof(double));
    for (k = 0; k < m; k++) {
        Hx(V[k], w);
        al[k] = dot(w, V[k]);
        #pragma omp parallel for
        for (long i = 0; i < D; i++) w[i] -= al[k] * V[k][i] + (k ? be[k] * V[k - 1][i] : 0);
        /* full reorthogonalization, classical Gram-Schmidt twice, by BLAS dgemv */
        for (int pass = 0; pass < 2; pass++) {
            int Dn = (int)D, k1 = k + 1, inc = 1;
            double one = 1, zero = 0, mone = -1;
            dgemv_("T", &Dn, &k1, &one, Vb, &Dn, w, &inc, &zero, hh, &inc);
            dgemv_("N", &Dn, &k1, &mone, Vb, &Dn, hh, &inc, &one, w, &inc);
        }
        be[k + 1] = sqrt(dot(w, w));
        used = k + 1;
        int check = ((k + 1) % 10 == 0) || be[k + 1] < 1e-12 || k + 1 == m;
        if (check) {
            int nn = used, ldz = nn, info;
            for (int i = 0; i < nn; i++) { d[i] = al[i]; e[i] = be[i + 1]; }
            free(Z); Z = malloc((size_t)nn * nn * sizeof(double));
            dstev_("V", &nn, d, e, Z, &ldz, work, &info);
            int nw = nwant < nn ? nwant : nn;
            conv = 1;
            for (int i = 0; i < nw; i++) if (fabs(be[k + 1] * Z[(size_t)i * nn + nn - 1]) > tol) conv = 0;
            if (nw < nwant && be[k + 1] >= 1e-12) conv = 0;
            if (conv || be[k + 1] < 1e-12 || k + 1 == m) {
                for (int i = 0; i < nw; i++) { ev[i] = d[i]; resid_out[i] = fabs(be[k + 1] * Z[(size_t)i * nn + nn - 1]); }
                for (int v = 0; ritz && v < nw; v++) {
                    double *x = ritz[v];
                    memset(x, 0, D * sizeof(double));
                    for (int j = 0; j < nn; j++) {
                        double c = Z[(size_t)v * nn + j];
                        #pragma omp parallel for
                        for (long i = 0; i < D; i++) x[i] += c * V[j][i];
                    }
                    double xn = sqrt(dot(x, x));
                    for (long i = 0; i < D; i++) x[i] /= xn;
                }
                if (be[k + 1] < 1e-12) conv = 1;
                break;
            }
        }
        V[k + 1] = Vb + (size_t)(k + 1) * D;
        #pragma omp parallel for
        for (long i = 0; i < D; i++) V[k + 1][i] = w[i] / be[k + 1];
    }
    free(Vb); free(hh);
    free(V); free(al); free(be); free(w); free(d); free(e); free(Z); free(work);
    *convp = conv;
    return used;
}

/* restarted Lanczos: up to nrestart passes; each restart starts from the sum
 * of the previous pass's lowest nwant Ritz vectors */

static int lanczos(int nwant, double *ev, int nv, double **vecs, double *resid_out)
{
    double *v0 = malloc(D * sizeof(double));
    uint64_t r = (uint64_t)seed * 0x9E3779B97F4A7C15ULL + 1;
    for (long i = 0; i < D; i++) { r ^= r << 13; r ^= r >> 7; r ^= r << 17; v0[i] = (double)(r >> 11) / 9007199254740992.0 - 0.5; }
    double nrm = sqrt(dot(v0, v0));
    for (long i = 0; i < D; i++) v0[i] /= nrm;
    int nw = nwant < D ? nwant : (int)D;
    double **ritz = malloc(nw * sizeof(double *));
    for (int v = 0; v < nw; v++) ritz[v] = malloc(D * sizeof(double));
    int total = 0, conv = 0, pass;
    for (pass = 0; pass < nrestart; pass++) {
        total += lanczos_pass(v0, nwant, ev, resid_out, ritz, &conv);
        if (conv) break;
        memset(v0, 0, D * sizeof(double));
        for (int v = 0; v < nw; v++) for (long i = 0; i < D; i++) v0[i] += ritz[v][i];
        nrm = sqrt(dot(v0, v0));
        for (long i = 0; i < D; i++) v0[i] /= nrm;
    }
    for (int v = 0; v < nv && v < nw; v++) memcpy(vecs[v], ritz[v], D * sizeof(double));
    for (int v = 0; v < nw; v++) free(ritz[v]);
    free(ritz); free(v0);
    fprintf(stderr, "lanczos: %d steps in %d passes, converged %d\n", total, pass < nrestart ? pass + 1 : nrestart, conv);
    return total;
}

/* ---------------- observables ---------------- */

/* vacancy indicator of cell h in a digit vector */
static int is_vac(const int *dg, int h)
{
    if (kind[h] == 1) return 0;
    int N = site[h]->N[dg[h]], u = 0;
    if (aomodel) for (int l = 0; l < nlink; l++) if (li[l] == h || lj[l] == h) u += dg[n + l];
    return N + u == site[h]->cap - 1;
}

static void observe(const double *x, int which_state)
{
    double *sz = calloc(n, sizeof(double)), *tz = calloc(n, sizeof(double)), *nv = calloc(n, sizeof(double));
    double *cz = calloc(n * n, sizeof(double)), *cl = calloc(nlink + 1, sizeof(double)), *Nc = calloc(n, sizeof(double));
    double *cond = calloc(n * n, sizeof(double));   /* <n_h Sz_j> */
    for (long s = 0; s < D; s++) {
        double p = x[s] * x[s];
        if (p == 0) continue;
        int dg[2 * NMAX];
        decode(codes[s], dg);
        double zi[NMAX];
        int vh[NMAX];
        for (int i = 0; i < n; i++) {
            zi[i] = site[i]->s2[dg[i]] / 2.0;
            sz[i] += p * zi[i]; tz[i] += p * site[i]->t2[dg[i]] / 2.0;
            Nc[i] += p * site[i]->N[dg[i]];
            vh[i] = is_vac(dg, i); nv[i] += p * vh[i];
        }
        for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) {
            cz[i * n + j] += p * zi[i] * zi[j];
            if (vh[i]) cond[i * n + j] += p * zi[j];
        }
        if (aomodel) for (int l = 0; l < nlink; l++) cl[l] += p * dg[n + l];
    }
    printf("state %d profile Sz", which_state); for (int i = 0; i < n; i++) printf(" %+.6f", sz[i]); printf("\n");
    if (M == 4) { printf("state %d profile Iz", which_state); for (int i = 0; i < n; i++) printf(" %+.6f", tz[i]); printf("\n"); }
    printf("state %d profile vac", which_state); for (int i = 0; i < n; i++) printf(" %.6f", nv[i]); printf("\n");
    printf("state %d profile N", which_state); for (int i = 0; i < n; i++) printf(" %.6f", Nc[i]); printf("\n");
    if (aomodel) { printf("state %d link c", which_state); for (int l = 0; l < nlink; l++) printf(" %.6f", cl[l]); printf("\n"); }
    if (obs_corr) {
        /* SzSz correlation, ring: averaged over i at distance r; open: from the middle cell */
        printf("state %d corr SzSz", which_state);
        if (ring) for (int r = 0; r <= n / 2; r++) { double a = 0; for (int i = 0; i < n; i++) a += cz[i * n + (i + r) % n]; printf(" %+.8f", a / n); }
        else { int c0 = n / 2; for (int r = 0; c0 + r < n; r++) printf(" %+.8f", cz[c0 * n + c0 + r]); }
        printf("\n");
        if (vac == 2) {
            /* conditional profile C(d) = sum_h <n_h Sz_{h+d}> / sum_h <n_h> (ring) */
            printf("state %d cond C(d)", which_state);
            double nt = 0; for (int h = 0; h < n; h++) nt += nv[h];
            if (ring) for (int dd = 0; dd <= n / 2; dd++) { double a = 0; for (int h = 0; h < n; h++) a += cond[h * n + (h + dd) % n]; printf(" %+.8f", a / nt); }
            printf("\n");
        }
    }
    free(sz); free(tz); free(nv); free(cz); free(cl); free(Nc); free(cond);
}

/* |CS> = prod_links A^+_l |0> (vac == 1: then b_{vsite, mu} on it), overlap with x */
static int cs_links[NMAX], cs_nl;

/* accumulate the |CS> vector exactly in a hash map, then overlap (handles coherent sums) */
typedef struct { int64_t *k; double *v; long cap, cnt; } HMap;
static HMap csmap;
static void hm_add(HMap *h, int64_t key, double val)
{
    if (2 * (h->cnt + 1) > h->cap) {
        HMap g = {malloc(2 * h->cap * sizeof(int64_t)), malloc(2 * h->cap * sizeof(double)), 2 * h->cap, 0};
        for (long i = 0; i < g.cap; i++) g.k[i] = -1;
        for (long i = 0; i < h->cap; i++) if (h->k[i] >= 0) hm_add(&g, h->k[i], h->v[i]);
        free(h->k); free(h->v); *h = g;
    }
    uint64_t m = h->cap - 1, p = hmix((uint64_t)key) & m;
    while (h->k[p] >= 0 && h->k[p] != key) p = (p + 1) & m;
    if (h->k[p] < 0) { h->k[p] = key; h->v[p] = 0; h->cnt++; }
    h->v[p] += val;
}
static void cs_rec2(int l, int occ[][4], double amp, int mu_v)
{
    if (l == cs_nl) {
        int o2[NMAX][4];
        memcpy(o2, occ, sizeof(int) * NMAX * 4);
        double a = amp;
        for (int i = 0; i < n; i++) for (int m = 0; m < M; m++) if (o2[i][m] == 2) a *= sqrt(2.0);
        if (mu_v >= 0) { int h = vsite; if (o2[h][mu_v] == 0) return; a *= sqrt((double)o2[h][mu_v]); o2[h][mu_v]--; }
        int64_t c = 0;
        for (int i = 0; i < n; i++) {
            int k = site[i]->key2idx[okey(o2[i])];
            if (k < 0) return;
            c += (int64_t)k * stride[i];
        }
        hm_add(&csmap, c, a);
        return;
    }
    int i = li[cs_links[l]], j = lj[cs_links[l]];
    for (int mu = 0; mu < M; mu++) for (int nu = 0; nu < M; nu++) {
        if (P[mu][nu] == 0) continue;
        if (occ[i][mu] >= 2 || occ[j][nu] >= 2) continue;
        occ[i][mu]++; occ[j][nu]++;
        cs_rec2(l + 1, occ, amp * P[mu][nu], mu_v);
        occ[i][mu]--; occ[j][nu]--;
    }
}
static void cs_overlap(const double *x, int which_state)
{
    cs_nl = nlink;
    for (int l = 0; l < nlink; l++) cs_links[l] = l;
    int mus[4], nmu = 0;
    if (vac == 1) { for (int m = 0; m < M; m++) if (-ms2[m] == sz2 && (M == 2 || -mt2[m] == iz2)) mus[nmu++] = m; }
    else mus[nmu++] = -1;
    if (vac == 2) { printf("state %d cs n/a (mobile)\n", which_state); return; }
    for (int u = 0; u < nmu; u++) {
        csmap.cap = 1 << 16; csmap.cnt = 0;
        csmap.k = malloc(csmap.cap * sizeof(int64_t)); csmap.v = malloc(csmap.cap * sizeof(double));
        for (long i = 0; i < csmap.cap; i++) csmap.k[i] = -1;
        int occ[NMAX][4]; memset(occ, 0, sizeof(occ));
        cs_rec2(0, occ, 1.0, mus[u]);
        double nr = 0, ov = 0, inb = 0;
        for (long i = 0; i < csmap.cap; i++) if (csmap.k[i] >= 0) {
            double v = csmap.v[i];
            nr += v * v;
            long q = findc(csmap.k[i]);
            if (q >= 0) { ov += v * x[q]; inb += v * v; }
        }
        double f = ov * ov / nr;
        printf("state %d cs mode %d overlap^2 %.8f per-cell %.6f (weight of |CS> inside the basis %.6f, terms %ld)\n",
               which_state, mus[u], f, pow(f, 1.0 / n), inb / nr, csmap.cnt);
        free(csmap.k); free(csmap.v);
    }
}

/* ---------------- main ---------------- */
static void parse(int argc, char **argv)
{
    for (int a = 1; a < argc; a++) {
        char k[64]; char v[64];
        if (sscanf(argv[a], "%63[^=]=%63s", k, v) != 2) { fprintf(stderr, "bad arg %s\n", argv[a]); exit(2); }
        if (!strcmp(k, "model")) aomodel = !strcmp(v, "ao");
        else if (!strcmp(k, "M")) M = atoi(v);
        else if (!strcmp(k, "n")) n = atoi(v);
        else if (!strcmp(k, "geom")) ring = !strcmp(v, "ring");
        else if (!strcmp(k, "vac")) vac = !strcmp(v, "static") ? 1 : !strcmp(v, "mobile") ? 2 : 0;
        else if (!strcmp(k, "vsite")) vsite = atoi(v);
        else if (!strcmp(k, "t")) tJ0 = atof(v);
        else if (!strcmp(k, "sz")) sz2 = atoi(v);
        else if (!strcmp(k, "iz")) iz2 = atoi(v);
        else if (!strcmp(k, "k")) kwant = atoi(v);
        else if (!strcmp(k, "nvec")) nvec = atoi(v);
        else if (!strcmp(k, "m")) mmax = atoi(v);
        else if (!strcmp(k, "cmax")) cmax = atoi(v);
        else if (!strcmp(k, "pen")) pen = atoi(v);
        else if (!strcmp(k, "lam")) lam_J0 = atof(v);
        else if (!strcmp(k, "gr")) gr = atof(v);
        else if (!strcmp(k, "tol")) tol = atof(v);
        else if (!strcmp(k, "seed")) seed = atoi(v);
        else if (!strcmp(k, "rs")) nrestart = atoi(v);
        else if (!strcmp(k, "prof")) obs_prof = atoi(v);
        else if (!strcmp(k, "corr")) obs_corr = atoi(v);
        else if (!strcmp(k, "cs")) obs_cs = atoi(v);
        else { fprintf(stderr, "unknown key %s\n", k); exit(2); }
    }
}

int main(int argc, char **argv)
{
    parse(argc, argv);
    constants(gr);
    t_rate = tJ0 * J0;
    memset(P, 0, sizeof(P));
    if (M == 2) { P[0][1] = 1; P[1][0] = -1; ms2[0] = 1; ms2[1] = -1; mt2[0] = mt2[1] = 0; }
    else {
        double eps[2][2] = {{0, 1}, {-1, 0}};
        for (int s = 0; s < 2; s++) for (int t = 0; t < 2; t++) for (int s2 = 0; s2 < 2; s2++) for (int t2 = 0; t2 < 2; t2++)
            P[2 * s + t][2 * s2 + t2] = eps[s][s2] * eps[t][t2];
        int sv[4] = {1, 1, -1, -1}, tv[4] = {1, -1, 1, -1};
        for (int a = 0; a < 4; a++) { ms2[a] = sv[a]; mt2[a] = tv[a]; }
    }
    /* cells */
    if (!aomodel) {
        loc_build(&locs[0], 2, 2, 2);        /* full */
        loc_build(&locs[1], 1, 1, 1);        /* end / static vacancy (N = 1) */
        loc_build(&locs[2], 1, 2, 2);        /* mobile */
    } else {
        loc_build(&locs[0], 0, 2, 2);        /* bulk, cap 2 */
        loc_build(&locs[1], 0, 1, 1);        /* end, cap 1 */
    }
    for (int i = 0; i < n; i++) {
        int end = !ring && (i == 0 || i == n - 1);
        kind[i] = end ? 1 : 0;
        if (!end && vac == 1 && i == vsite) kind[i] = 2;
        if (!end && vac == 2) kind[i] = 3;
        if (!aomodel) site[i] = end ? &locs[1] : kind[i] == 2 ? &locs[1] : kind[i] == 3 ? &locs[2] : &locs[0];
        else site[i] = end ? &locs[1] : &locs[0];
    }
    nlink = ring ? n : n - 1;
    for (int l = 0; l < nlink; l++) { li[l] = l; lj[l] = (l + 1) % n; }
    nsl = aomodel ? n + nlink : n;
    clmax = 2;
    if (cmax >= 0 && cmax < clmax) clmax = cmax > 0 ? cmax : 0;
    for (int i = 0; i < n; i++) radix[i] = site[i]->d;
    for (int l = 0; l < nlink && aomodel; l++) radix[n + l] = clmax + 1;
    stride[nsl - 1] = 1;
    for (int i = nsl - 2; i >= 0; i--) stride[i] = stride[i + 1] * radix[i + 1];
    if (!aomodel)
        for (int l = 0; l < nlink; l++)
            two_build(&two[l], site[li[l]], site[lj[l]], kind[li[l]] == 3, kind[lj[l]] == 3, t_rate);
    if (aomodel) {
        Eref = 0;
        for (int i = 0; i < n; i++) Eref += EF[site[i]->cap - (vac == 1 && i == vsite ? 1 : 0)];
        if (vac == 2) Eref -= EF[2] - EF[1];
    }
    double t0 = omp_get_wtime();
    for (int w = 0; w < (M == 4 ? 2 : 1); w++) {
        if (w == 0) sz2 += 2; else iz2 += 2;
        if (aomodel) basis_ao(); else basis_lo();
        auxc[w] = codes; auxD[w] = D;
        phi[w] = malloc((D + 1) * sizeof(double));
        if (w == 0) sz2 -= 2; else iz2 -= 2;
    }
    if (aomodel) basis_ao(); else basis_lo();
    build_diag();
    lam = pen ? lam_J0 * J0 : 0;
    printf("# vx model=%s M=%d n=%d geom=%s vac=%s vsite=%d t=%g J0 sector=(%d,%d) cmax=%d pen=%d gr=%g  D=%ld  J0=%.10f  [basis %.1fs]\n",
           aomodel ? "ao" : "lo", M, n, ring ? "ring" : "open", vac == 0 ? "none" : vac == 1 ? "static" : "mobile", vsite, tJ0,
           sz2, iz2, cmax, pen, gr, D, J0, omp_get_wtime() - t0);
    if (D == 0) { printf("empty sector\n"); return 0; }
    double *ev = calloc(kwant + 1, sizeof(double)), *res = calloc(kwant + 1, sizeof(double));
    double **vecs = malloc((nvec + 1) * sizeof(double *));
    for (int v = 0; v < nvec; v++) vecs[v] = malloc(D * sizeof(double));
    int steps = lanczos(kwant, ev, nvec, vecs, res);
    for (int i = 0; i < kwant && i < D; i++)
        printf("E %d %.12f rate  %.8f J0  rel %.8f J0  resid %.1e\n", i, ev[i], ev[i] / J0, (ev[i] - ev[0]) / J0, res[i]);
    for (int v = 0; v < nvec && v < D; v++) {
        double S2 = expect_sq(0, vecs[v]), I2 = M == 4 ? expect_sq(1, vecs[v]) : 0;
        printf("state %d E %.8f J0  <S2> %.6f  <I2> %.6f\n", v, ev[v] / J0, S2, I2);
        if (obs_prof || obs_corr) observe(vecs[v], v);
        if (obs_cs) cs_overlap(vecs[v], v);
    }
    printf("# steps %d  [%.1fs]\n", steps, omp_get_wtime() - t0);
    return 0;
}
