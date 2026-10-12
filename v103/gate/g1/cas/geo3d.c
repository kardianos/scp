/* geo3d.c — v103 gate g1: can a medium with no container hold three
 * dimensions, and what does it do with a gap? (../PREREG.md)
 *
 * The medium is N quanta and their contacts, nothing else: no coordinates,
 * no cells, no slots.  Every quantum has six contacts, all of them paired
 * (an unpaired contact would be an edge inside the medium), and two quanta
 * touch at most once.  "Where" is counted through the medium: the distance
 * between two quanta is the number of contacts on the shortest path.  The
 * medium re-forms by local flips: along a path i-j-k-l, the contacts ij
 * and kl become ik and jl.  A flip keeps six contacts per quantum and keeps
 * the medium in one piece.
 *
 * candidates (all 8 combinations); energy in units of t, 2t per unit
 * deviation squared (c2's matched scale):
 *   a  flat squares:     every contact lies in exactly four squares
 *   b  octahedral link:  at every quantum, each contact shares a square
 *                        with exactly four of the other five
 *   c  no triangles:     2t per triangle
 * The state is sampled in the random process at temperature t (Gibbs).
 *
 * sectors per run (PREREG §4): M the medium from the flat start; G one gap
 * (a quantum removed, its six contacts re-paired by Gibbs weight); X one
 * excess (a quantum added into three contacts by Gibbs weight); A assembly
 * (the flat start melted, then cooled).
 *
 * The flat start is built from bookkeeping coordinates that are dropped at
 * once; no move and no measurement uses a coordinate.
 *
 * build: gcc -O3 -march=native -Wall -Wextra -o geo3d geo3d.c \
 *            -l:liblapack.so.3 -l:libopenblas.so.0 -lm
 * usage: ./geo3d mode=run L=12 a=1 b=0 c=1     (one registered run)
 *        ./geo3d mode=t0                       (tool gates, PREREG §6)
 */
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

extern void dsyev_(const char *, const char *, const int *, double *, const int *, double *, double *, const int *, int *);
extern void dstev_(const char *, const int *, double *, double *, double *, const int *, double *, int *);

#define DEG 6
#define RMAX 256
static int OPT_A, OPT_B, OPT_C;
static const double TREG = 1.0; /* the registered temperature, t = 1 */

/* the registered schedule (PREREG §4) */
static const int EQ_SWEEPS = 2000, MEAS_SWEEPS = 2000, SAMPLE_EVERY = 100, MELT_SWEEPS = 2000, ANNEAL_SWEEPS = 300;
static const double ANNEAL_T[] = {64, 32, 16, 8, 4, 2, 1};
#define NANNEAL 7

/* ---------------- rng: xoshiro256** seeded by splitmix64 ---------------- */
typedef struct {
    uint64_t s[4];
} Rng;
static uint64_t splitmix64(uint64_t *x)
{
    uint64_t z = (*x += 0x9e3779b97f4a7c15ULL);
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}
static void rng_seed(Rng *r, uint64_t seed)
{
    uint64_t x = seed;
    for (int i = 0; i < 4; i++) r->s[i] = splitmix64(&x);
}
static inline uint64_t rotl(uint64_t x, int k) { return (x << k) | (x >> (64 - k)); }
static inline uint64_t rnext(Rng *r)
{
    uint64_t *s = r->s, res = rotl(s[1] * 5, 7) * 9, t = s[1] << 17;
    s[2] ^= s[0];
    s[3] ^= s[1];
    s[1] ^= s[2];
    s[0] ^= s[3];
    s[2] ^= t;
    s[3] = rotl(s[3], 45);
    return res;
}
static inline double runif(Rng *r) { return (double)(rnext(r) >> 11) * 0x1.0p-53; }
static inline int rint_(Rng *r, int n) { return (int)(((rnext(r) >> 32) * (uint64_t)n) >> 32); }

/* ---------------- the medium: quanta and contacts ---------------- */
typedef struct {
    int n, cap, W;
    int *deg, *nb; /* nb[i*DEG + k] */
    uint64_t *adj; /* cap x W bits */
} G;

static G *g_new(int cap)
{
    G *g = calloc(1, sizeof *g);
    g->cap = cap;
    g->W = (cap + 63) / 64;
    g->deg = calloc((size_t)cap, sizeof(int));
    g->nb = calloc((size_t)cap * DEG, sizeof(int));
    g->adj = calloc((size_t)cap * g->W, sizeof(uint64_t));
    if (!g->deg || !g->nb || !g->adj) {
        fprintf(stderr, "out of memory\n");
        exit(2);
    }
    return g;
}
static void g_free(G *g)
{
    free(g->deg);
    free(g->nb);
    free(g->adj);
    free(g);
}
static void g_reset(G *g, int n)
{
    g->n = n;
    memset(g->deg, 0, sizeof(int) * (size_t)g->cap);
    memset(g->adj, 0, sizeof(uint64_t) * (size_t)g->cap * g->W);
}
static void g_copy(G *d, const G *s)
{
    if (d->cap != s->cap) {
        fprintf(stderr, "g_copy: cap mismatch\n");
        exit(2);
    }
    d->n = s->n;
    memcpy(d->deg, s->deg, sizeof(int) * (size_t)s->cap);
    memcpy(d->nb, s->nb, sizeof(int) * (size_t)s->cap * DEG);
    memcpy(d->adj, s->adj, sizeof(uint64_t) * (size_t)s->cap * s->W);
}
static inline int ADJ(const G *g, int i, int j) { return (int)((g->adj[(size_t)i * g->W + (j >> 6)] >> (j & 63)) & 1u); }
static inline void setb(G *g, int i, int j) { g->adj[(size_t)i * g->W + (j >> 6)] |= 1ULL << (j & 63); }
static inline void clrb(G *g, int i, int j) { g->adj[(size_t)i * g->W + (j >> 6)] &= ~(1ULL << (j & 63)); }
static void add_edge(G *g, int i, int j)
{
    if (g->deg[i] >= DEG || g->deg[j] >= DEG || i == j || ADJ(g, i, j)) {
        fprintf(stderr, "add_edge: invalid %d-%d\n", i, j);
        exit(2);
    }
    g->nb[i * DEG + g->deg[i]++] = j;
    g->nb[j * DEG + g->deg[j]++] = i;
    setb(g, i, j);
    setb(g, j, i);
}
static void rm_nb(G *g, int i, int j)
{
    int *p = g->nb + i * DEG;
    for (int k = 0; k < g->deg[i]; k++)
        if (p[k] == j) {
            p[k] = p[--g->deg[i]];
            return;
        }
    fprintf(stderr, "rm_nb: %d not a contact of %d\n", j, i);
    exit(2);
}
static void del_edge(G *g, int i, int j)
{
    rm_nb(g, i, j);
    rm_nb(g, j, i);
    clrb(g, i, j);
    clrb(g, j, i);
}
/* drop quantum q (no contacts left): the last quantum's bookkeeping index becomes q */
static void g_drop(G *g, int q)
{
    int last = g->n - 1;
    if (g->deg[q]) {
        fprintf(stderr, "g_drop: quantum still has contacts\n");
        exit(2);
    }
    if (q != last) {
        for (int k = 0; k < g->deg[last]; k++) {
            int u = g->nb[last * DEG + k];
            int *p = g->nb + u * DEG;
            for (int m = 0; m < g->deg[u]; m++)
                if (p[m] == last) p[m] = q;
            clrb(g, u, last);
            setb(g, u, q);
            g->nb[q * DEG + k] = u;
            setb(g, q, u);
        }
        g->deg[q] = g->deg[last];
        g->deg[last] = 0;
        memset(g->adj + (size_t)last * g->W, 0, sizeof(uint64_t) * g->W);
    }
    g->n--;
}

/* ---------------- initial patterns (bookkeeping coordinates, dropped at once) ---------------- */
static void build_cubic(G *g, int L)
{
    g_reset(g, L * L * L);
    for (int z = 0; z < L; z++)
        for (int y = 0; y < L; y++)
            for (int x = 0; x < L; x++) {
                int i = x + L * (y + L * z);
                add_edge(g, i, (x + 1) % L + L * (y + L * z));
                add_edge(g, i, x + L * ((y + 1) % L + L * z));
                add_edge(g, i, x + L * (y + L * ((z + 1) % L)));
            }
}
static void build_tri(G *g, int L) /* the 2-D triangular arrangement, closed */
{
    g_reset(g, L * L);
    for (int y = 0; y < L; y++)
        for (int x = 0; x < L; x++) {
            int i = x + L * y;
            add_edge(g, i, (x + 1) % L + L * y);
            add_edge(g, i, x + L * ((y + 1) % L));
            add_edge(g, i, (x + 1) % L + L * ((y - 1 + L) % L));
        }
}
static void build_circ(G *g, int n, int a, int b) /* contacts i~i+1, i~i+a, i~i+b mod n */
{
    g_reset(g, n);
    for (int i = 0; i < n; i++) {
        add_edge(g, i, (i + 1) % n);
        add_edge(g, i, (i + a) % n);
        add_edge(g, i, (i + b) % n);
    }
}

/* ---------------- local structure ---------------- */
static int sq(const G *g, int u, int v) /* squares through contact u-v */
{
    int c = 0;
    const int *Nu = g->nb + u * DEG, *Nv = g->nb + v * DEG;
    for (int a = 0; a < g->deg[u]; a++) {
        int w = Nu[a];
        if (w == v) continue;
        for (int b = 0; b < g->deg[v]; b++) {
            int x = Nv[b];
            if (x == u || x == w) continue;
            if (ADJ(g, w, x)) c++;
        }
    }
    return c;
}
static int lk(const G *g, int i, int j) /* contacts k of i (k != j) that share a square with i and j */
{
    int c = 0;
    const int *Ni = g->nb + i * DEG, *Nj = g->nb + j * DEG;
    for (int a = 0; a < g->deg[i]; a++) {
        int k = Ni[a];
        if (k == j) continue;
        for (int b = 0; b < g->deg[j]; b++) {
            int m = Nj[b];
            if (m == i || m == k) continue;
            if (ADJ(g, m, k)) {
                c++;
                break;
            }
        }
    }
    return c;
}
static int cn(const G *g, int u, int v) /* common contacts of u and v (triangles through u-v) */
{
    int c = 0;
    const int *Nu = g->nb + u * DEG;
    for (int a = 0; a < g->deg[u]; a++)
        if (Nu[a] != v && ADJ(g, Nu[a], v)) c++;
    return c;
}

/* ---------------- energy (in units of 2t) ---------------- */
static char *MARK;
static int *VL, SCAP;
static void scratch(int cap)
{
    if (cap <= SCAP) return;
    free(MARK);
    free(VL);
    MARK = calloc((size_t)cap, 1);
    VL = malloc(sizeof(int) * (size_t)cap);
    SCAP = cap;
}
static inline int anyopt(void) { return OPT_A || OPT_B || OPT_C; }
/* energy carried by the structures inside the marked set V (all = every quantum) */
static long eloc_impl(const G *g, const int *V, int nv, int all)
{
    long e = 0;
    for (int a = 0; a < nv; a++) {
        int u = all ? a : V[a];
        const int *Nu = g->nb + u * DEG;
        for (int k = 0; k < g->deg[u]; k++) {
            int v = Nu[k];
            if (v <= u || !(all || MARK[v])) continue;
            if (OPT_A) {
                int s = sq(g, u, v) - 4;
                e += (long)s * s;
            }
            if (OPT_C)
                for (int m = 0; m < g->deg[u]; m++) {
                    int w = Nu[m];
                    if (w > v && (all || MARK[w]) && ADJ(g, v, w)) e++;
                }
        }
        if (OPT_B)
            for (int k = 0; k < g->deg[u]; k++) {
                int q = lk(g, u, Nu[k]) - 4;
                e += (long)q * q;
            }
    }
    return e;
}
static long eloc(const G *g, const int *V, int nv) { return eloc_impl(g, V, nv, 0); }
static long efull(const G *g) { return eloc_impl(g, NULL, g->n, 1); }
/* V1 = V0 and their contacts, marked; returns its size (in VL) */
static int build_V1(const G *g, const int *V0, int n0)
{
    int nv = 0;
    for (int a = 0; a < n0; a++)
        if (!MARK[V0[a]]) {
            MARK[V0[a]] = 1;
            VL[nv++] = V0[a];
        }
    for (int a = 0; a < n0; a++)
        for (int k = 0; k < g->deg[V0[a]]; k++) {
            int y = g->nb[V0[a] * DEG + k];
            if (!MARK[y]) {
                MARK[y] = 1;
                VL[nv++] = y;
            }
        }
    return nv;
}
static void unmark(int nv)
{
    for (int a = 0; a < nv; a++) MARK[VL[a]] = 0;
}
/* the three term counts, independent of the options: sum (s-4)^2, sum (q-4)^2, triangles */
static void term_counts(const G *g, long *A2, long *B2, long *T3)
{
    *A2 = *B2 = *T3 = 0;
    for (int u = 0; u < g->n; u++)
        for (int k = 0; k < g->deg[u]; k++) {
            int v = g->nb[u * DEG + k], q = lk(g, u, v) - 4;
            *B2 += (long)q * q;
            if (v <= u) continue;
            int s = sq(g, u, v) - 4;
            *A2 += (long)s * s;
            *T3 += cn(g, u, v);
        }
    *T3 /= 3;
}

/* ---------------- the re-forming: local flips ---------------- */
typedef struct {
    int i, j, k, l;
} Flip;
static long NPROP, NACC;
static int pick_other(const G *g, Rng *r, int j, int excl)
{
    int a = rint_(r, g->deg[j] - 1), c = 0;
    const int *Nj = g->nb + j * DEG;
    for (int m = 0; m < g->deg[j]; m++) {
        if (Nj[m] == excl) continue;
        if (c++ == a) return Nj[m];
    }
    return -1;
}
/* uniform directed contact j->k, uniform i in N(j)\k, uniform l in N(k)\j: symmetric (PREREG §2) */
static int propose(const G *g, Rng *r, Flip *f)
{
    int j = rint_(r, g->n), k = g->nb[j * DEG + rint_(r, g->deg[j])];
    int i = pick_other(g, r, j, k), l = pick_other(g, r, k, j);
    if (i < 0 || l < 0 || i == l || ADJ(g, i, k) || ADJ(g, j, l)) return 0;
    f->i = i;
    f->j = j;
    f->k = k;
    f->l = l;
    return 1;
}
static void flip_apply(G *g, const Flip *f)
{
    del_edge(g, f->i, f->j);
    del_edge(g, f->k, f->l);
    add_edge(g, f->i, f->k);
    add_edge(g, f->j, f->l);
}
static void flip_revert(G *g, const Flip *f)
{
    del_edge(g, f->i, f->k);
    del_edge(g, f->j, f->l);
    add_edge(g, f->i, f->j);
    add_edge(g, f->k, f->l);
}
/* the energy change of a flip, in units of 2t; the flip is left undone */
static long flip_de(G *g, const Flip *f)
{
    int V0[4] = {f->i, f->j, f->k, f->l};
    int nv = build_V1(g, V0, 4);
    long e0 = eloc(g, VL, nv);
    flip_apply(g, f);
    long e1 = eloc(g, VL, nv);
    flip_revert(g, f);
    unmark(nv);
    return e1 - e0;
}
static int mc_step(G *g, Rng *r, double T, int free_)
{
    Flip f;
    NPROP++;
    if (!propose(g, r, &f)) return 0;
    if (free_ || !anyopt()) {
        flip_apply(g, &f);
        NACC++;
        return 1;
    }
    int V0[4] = {f.i, f.j, f.k, f.l};
    int nv = build_V1(g, V0, 4);
    long e0 = eloc(g, VL, nv);
    flip_apply(g, &f);
    long e1 = eloc(g, VL, nv);
    double dE = 2.0 * (double)(e1 - e0); /* in t */
    int acc = dE <= 0 || runif(r) < exp(-dE / T);
    if (acc)
        NACC++;
    else
        flip_revert(g, &f);
    unmark(nv);
    return acc;
}
static void sweeps(G *g, Rng *r, int ns, double T, int free_)
{
    for (int s = 0; s < ns; s++) {
        long m = 3L * g->n; /* one sweep = one proposal per contact */
        for (long p = 0; p < m; p++) mc_step(g, r, T, free_);
    }
}

/* ---------------- graph instruments ---------------- */
static int bfs(const G *g, const int *src, int ns, int *dist, int *que)
{
    for (int i = 0; i < g->n; i++) dist[i] = -1;
    int h = 0, t = 0;
    for (int a = 0; a < ns; a++)
        if (dist[src[a]] < 0) {
            dist[src[a]] = 0;
            que[t++] = src[a];
        }
    while (h < t) {
        int u = que[h++];
        for (int k = 0; k < g->deg[u]; k++) {
            int v = g->nb[u * DEG + k];
            if (dist[v] < 0) {
                dist[v] = dist[u] + 1;
                que[t++] = v;
            }
        }
    }
    return t; /* number reached */
}
static int pieces(const G *g)
{
    int *dist = malloc(sizeof(int) * (size_t)g->n), *que = malloc(sizeof(int) * (size_t)g->n), *seen = calloc((size_t)g->n, sizeof(int));
    int np = 0;
    for (int s = 0; s < g->n; s++) {
        if (seen[s]) continue;
        np++;
        int nr = bfs(g, &s, 1, dist, que);
        for (int a = 0; a < nr; a++) seen[que[a]] = 1;
    }
    free(dist);
    free(que);
    free(seen);
    return np;
}
static void lap_mv(const G *g, const double *x, double *y)
{
    for (int i = 0; i < g->n; i++) {
        double s = g->deg[i] * x[i];
        for (int k = 0; k < g->deg[i]; k++) s -= x[g->nb[i * DEG + k]];
        y[i] = s;
    }
}
static double dot(const double *a, const double *b, int n)
{
    double s = 0;
    for (int i = 0; i < n; i++) s += a[i] * b[i];
    return s;
}
static void center(double *a, int n)
{
    double s = 0;
    for (int i = 0; i < n; i++) s += a[i];
    s /= n;
    for (int i = 0; i < n; i++) a[i] -= s;
}
/* lowest nonzero eigenvalue of the contact Laplacian: Lanczos with full reorthogonalisation */
static double lambda2(const G *g, Rng *r, int *ok)
{
    int n = g->n, mmax = n - 1 < 400 ? n - 1 : 400;
    double *V = malloc(sizeof(double) * (size_t)(mmax + 1) * n), *w = malloc(sizeof(double) * n);
    double *al = malloc(sizeof(double) * mmax), *be = malloc(sizeof(double) * mmax);
    double *d = malloc(sizeof(double) * mmax), *e = malloc(sizeof(double) * mmax);
    double *Z = malloc(sizeof(double) * (size_t)mmax * mmax), *work = malloc(sizeof(double) * (2 * (size_t)mmax + 2));
    for (int i = 0; i < n; i++) V[i] = runif(r) - 0.5;
    center(V, n);
    double nr = sqrt(dot(V, V, n));
    for (int i = 0; i < n; i++) V[i] /= nr;
    double theta = NAN;
    *ok = 0;
    for (int m = 0; m < mmax; m++) {
        double *vj = V + (size_t)m * n;
        lap_mv(g, vj, w);
        al[m] = dot(w, vj, n);
        for (int pass = 0; pass < 2; pass++) {
            for (int p = 0; p <= m; p++) {
                double *vp = V + (size_t)p * n, c = dot(w, vp, n);
                for (int i = 0; i < n; i++) w[i] -= c * vp[i];
            }
            center(w, n);
        }
        be[m] = sqrt(dot(w, w, n));
        int mm = m + 1;
        if (mm % 20 == 0 || mm == mmax || be[m] < 1e-10) {
            for (int i = 0; i < mm; i++) d[i] = al[i];
            for (int i = 0; i < mm - 1; i++) e[i] = be[i];
            int info;
            dstev_("V", &mm, d, e, Z, &mm, work, &info);
            if (info) break;
            theta = d[0];
            double res = fabs(be[m] * Z[mm - 1]);
            if (res < 1e-10 * (1 + fabs(theta)) || be[m] < 1e-10) {
                *ok = 1;
                break;
            }
        }
        double *vn = V + (size_t)(m + 1) * n;
        for (int i = 0; i < n; i++) vn[i] = w[i] / be[m];
    }
    free(V);
    free(w);
    free(al);
    free(be);
    free(d);
    free(e);
    free(Z);
    free(work);
    return theta;
}
static double lambda2_dense(const G *g)
{
    int n = g->n, lw = 3 * n, info;
    double *A = calloc((size_t)n * n, sizeof(double)), *wv = malloc(sizeof(double) * n), *wk = malloc(sizeof(double) * lw);
    for (int i = 0; i < n; i++) {
        A[(size_t)i * n + i] = g->deg[i];
        for (int k = 0; k < g->deg[i]; k++) A[(size_t)i * n + g->nb[i * DEG + k]] -= 1;
    }
    dsyev_("N", "U", &n, A, &n, wv, wk, &lw, &info);
    double l2 = info ? NAN : wv[1];
    free(A);
    free(wv);
    free(wk);
    return l2;
}
/* the medium's Green's function from a unit source spread over C (zero mean), by CG */
static int green(const G *g, const int *C, int nc, double *x)
{
    int n = g->n;
    double *b = malloc(sizeof(double) * n), *r = malloc(sizeof(double) * n), *p = malloc(sizeof(double) * n), *Ap = malloc(sizeof(double) * n);
    for (int i = 0; i < n; i++) b[i] = -1.0 / n;
    for (int c = 0; c < nc; c++) b[C[c]] += 1.0 / nc;
    for (int i = 0; i < n; i++) {
        x[i] = 0;
        r[i] = b[i];
        p[i] = b[i];
    }
    double rr = dot(r, r, n), bn = sqrt(rr);
    int it, conv = -1;
    for (it = 0; it < 50000; it++) {
        lap_mv(g, p, Ap);
        double alpha = rr / dot(p, Ap, n);
        for (int i = 0; i < n; i++) {
            x[i] += alpha * p[i];
            r[i] -= alpha * Ap[i];
        }
        center(r, n);
        double rn = dot(r, r, n);
        if (sqrt(rn) < 1e-11 * bn) {
            conv = it + 1;
            break;
        }
        double beta = rn / rr;
        for (int i = 0; i < n; i++) p[i] = r[i] + beta * p[i];
        rr = rn;
    }
    center(x, n);
    free(b);
    free(r);
    free(p);
    free(Ap);
    return conv;
}
/* mean distance over pairs of distinct quanta, diameter (one piece assumed) */
static void apsp_basic(const G *g, double *meand, int *diam)
{
    int n = g->n, *dist = malloc(sizeof(int) * n), *que = malloc(sizeof(int) * n);
    double tot = 0;
    int dm = 0;
    for (int s = 0; s < n; s++) {
        int nr = bfs(g, &s, 1, dist, que);
        for (int a = 0; a < nr; a++) {
            int d = dist[que[a]];
            tot += d;
            if (d > dm) dm = d;
        }
    }
    *meand = tot / ((double)n * (n - 1));
    *diam = dm;
    free(dist);
    free(que);
}
static int check(const G *g, int nexp)
{
    int ok = g->n == nexp;
    for (int i = 0; i < g->n && ok; i++) {
        if (g->deg[i] != DEG) ok = 0;
        int pc = 0;
        for (int w = 0; w < g->W; w++) pc += __builtin_popcountll(g->adj[(size_t)i * g->W + w]);
        if (pc != g->deg[i]) ok = 0;
        for (int k = 0; k < g->deg[i]; k++) {
            int j = g->nb[i * DEG + k];
            if (j == i || j < 0 || j >= g->n || !ADJ(g, i, j) || !ADJ(g, j, i)) ok = 0;
            for (int m = k + 1; m < g->deg[i]; m++)
                if (g->nb[i * DEG + m] == j) ok = 0;
        }
    }
    if (ok && pieces(g) != 1) ok = 0;
    return ok;
}

/* ---------------- union-find ---------------- */
static int uf_find(int *p, int x)
{
    while (p[x] != x) {
        p[x] = p[p[x]];
        x = p[x];
    }
    return x;
}
static void uf_union(int *p, int a, int b)
{
    a = uf_find(p, a);
    b = uf_find(p, b);
    if (a != b) p[a] = b;
}

/* ---------------- measurement (PREREG §5) ---------------- */
/* mode 0: medium readings only; 1: reference profiles from quantum 0; 2: footprint around the mismatch centre */
static void measure(G *g, Rng *r, const char *sec, int k, int mode, double *Sacc, int *Smax)
{
    int n = g->n, B = 0;
    int *eu = malloc(sizeof(int) * 3 * (size_t)n), *ev = malloc(sizeof(int) * 3 * (size_t)n), *es = malloc(sizeof(int) * 3 * (size_t)n);
    char *emis = malloc(3 * (size_t)n), *okv = malloc((size_t)n), *isD = calloc((size_t)n, 1);
    long A2 = 0, B2 = 0, T3 = 0, sqdef = 0;
    int ns4 = 0, nmis = 0;
    memset(okv, 1, (size_t)n);
    for (int u = 0; u < n; u++)
        for (int a = 0; a < g->deg[u]; a++) {
            int v = g->nb[u * DEG + a];
            if (v <= u) continue;
            int s = sq(g, u, v), t = cn(g, u, v), l1 = lk(g, u, v), l2 = lk(g, v, u);
            eu[B] = u;
            ev[B] = v;
            es[B] = s;
            emis[B] = (s != 4 || t > 0 || l1 != 4 || l2 != 4);
            if (s == 4) ns4++;
            sqdef += 4 - s;
            A2 += (long)(s - 4) * (s - 4);
            B2 += (long)(l1 - 4) * (l1 - 4) + (long)(l2 - 4) * (l2 - 4);
            T3 += t;
            if (l1 != 4) okv[u] = 0;
            if (l2 != 4) okv[v] = 0;
            if (emis[B]) {
                nmis++;
                isD[u] = isD[v] = 1;
            }
            B++;
        }
    T3 /= 3;
    int nO = 0, nD = 0;
    for (int i = 0; i < n; i++) {
        nO += okv[i];
        nD += isD[i];
    }
    double E = 2.0 * ((OPT_A ? A2 : 0) + (OPT_B ? B2 : 0) + (OPT_C ? T3 : 0));
    int np = pieces(g);
    /* all-pairs distances */
    int *dist = malloc(sizeof(int) * n), *que = malloc(sizeof(int) * n), *eccD = malloc(sizeof(int) * n), *uf = malloc(sizeof(int) * n);
    double *mx = malloc(sizeof(double) * n);
    long *cnt = calloc(RMAX, sizeof(long));
    for (int i = 0; i < n; i++) uf[i] = i;
    double tot = 0;
    int diam = 0;
    for (int s = 0; s < n; s++) {
        int nr = bfs(g, &s, 1, dist, que), ecd = 0;
        double sm = 0;
        for (int a = 0; a < nr; a++) {
            int y = que[a], d = dist[y];
            sm += d;
            if (d > diam) diam = d;
            if (d < RMAX) cnt[d]++;
            if (isD[y] && d > ecd) ecd = d;
            if (isD[s] && isD[y] && d <= 2) uf_union(uf, s, y);
        }
        mx[s] = nr > 1 ? sm / (nr - 1) : 0;
        tot += mx[s];
        eccD[s] = ecd;
    }
    for (int d = 0; d < RMAX; d++) {
        Sacc[d] += (double)cnt[d];
        if (cnt[d] && d > *Smax) *Smax = d;
    }
    double meand = tot / n;
    int ncl = 0, dspread = -1, minecc = 1 << 30, nC = 0;
    for (int i = 0; i < n; i++)
        if (isD[i]) {
            if (uf_find(uf, i) == i) ncl++;
            if (eccD[i] > dspread) dspread = eccD[i];
        }
    int lamok = 0;
    double lam2 = np == 1 ? lambda2(g, r, &lamok) : NAN;
    printf("SAMPLE sec=%s k=%d n=%d E=%.6f psiT=%.8f psiO=%.8f tri=%.8f sqdef=%ld meand=%.9f diam=%d pieces=%d lam2=%.12g lamok=%d nmis=%d nD=%d ncl=%d dspread=%d\n",
           sec, k, n, E, (double)ns4 / B, (double)nO / n, (double)T3 / n, sqdef, meand, diam, np, lam2, lamok, nmis, nD, ncl, dspread);
    if (mode == 1 || (mode == 2 && nD > 0)) {
        int *C = malloc(sizeof(int) * n);
        if (mode == 1) {
            C[0] = 0;
            nC = 1;
        } else {
            for (int i = 0; i < n; i++)
                if (eccD[i] < minecc) minecc = eccD[i];
            for (int i = 0; i < n; i++)
                if (eccD[i] == minecc) C[nC++] = i;
        }
        bfs(g, C, nC, dist, que);
        double *phi = malloc(sizeof(double) * n);
        int gconv = green(g, C, nC, phi);
        int rmax = 0;
        for (int i = 0; i < n; i++)
            if (dist[i] > rmax) rmax = dist[i];
        if (rmax >= RMAX) rmax = RMAX - 1;
        long *sn = calloc(RMAX, sizeof(long)), *ssq = calloc(RMAX, sizeof(long)), *smis = calloc(RMAX, sizeof(long));
        double *sdm = calloc(RMAX, sizeof(double)), *sgf = calloc(RMAX, sizeof(double));
        for (int i = 0; i < n; i++) {
            int d = dist[i] < RMAX ? dist[i] : RMAX - 1;
            sn[d]++;
            sdm[d] += mx[i];
            sgf[d] += phi[i];
        }
        for (int e = 0; e < B; e++) {
            int d = dist[eu[e]] < dist[ev[e]] ? dist[eu[e]] : dist[ev[e]];
            if (d >= RMAX) d = RMAX - 1;
            ssq[d] += 4 - es[e];
            smis[d] += emis[e];
        }
        printf("CENTRE sec=%s k=%d nC=%d ecc=%d greenit=%d mmean=%.9f\n", sec, k, nC, mode == 2 ? minecc : 0, gconv, meand);
        for (int d = 0; d <= rmax; d++)
            if (sn[d] || ssq[d] || smis[d])
                printf("PROF sec=%s k=%d r=%d n=%ld sqd=%ld mis=%ld dm=%.9f gf=%.12g\n", sec, k, d, sn[d], ssq[d], smis[d], sn[d] ? sdm[d] / sn[d] : 0.0,
                       sn[d] ? sgf[d] / sn[d] : 0.0);
        free(sn);
        free(ssq);
        free(smis);
        free(sdm);
        free(sgf);
        free(phi);
        free(C);
    }
    fflush(stdout);
    free(eu);
    free(ev);
    free(es);
    free(emis);
    free(okv);
    free(isD);
    free(dist);
    free(que);
    free(eccD);
    free(uf);
    free(mx);
    free(cnt);
}

/* ---------------- making a gap and an excess (PREREG §4) ---------------- */
static int PM[15][6];
static void init_pm(void)
{
    int c = 0;
    for (int x = 1; x < 6; x++) {
        int rest[4], nr = 0;
        for (int y = 1; y < 6; y++)
            if (y != x) rest[nr++] = y;
        for (int z = 1; z < 4; z++) {
            int o[2], no = 0;
            for (int w = 1; w < 4; w++)
                if (w != z) o[no++] = rest[w];
            int *P = PM[c++];
            P[0] = 0;
            P[1] = x;
            P[2] = rest[0];
            P[3] = rest[z];
            P[4] = o[0];
            P[5] = o[1];
        }
    }
}
/* energy change (2t units) of each re-pairing when q is removed; valid[p] = 0 if it would double a contact.
 * If dfull != NULL, also the full recomputation (gate T0.6). Leaves the medium unchanged. */
static int gap_candidates(G *g, int q, int *valid, long *dloc, long *dfull)
{
    int J[DEG];
    memcpy(J, g->nb + q * DEG, sizeof J);
    int V0[DEG + 1];
    V0[0] = q;
    for (int m = 0; m < DEG; m++) V0[m + 1] = J[m];
    int nv = build_V1(g, V0, DEG + 1);
    long e0 = eloc(g, VL, nv), E0 = dfull ? efull(g) : 0;
    int nvalid = 0;
    for (int p = 0; p < 15; p++) {
        const int *P = PM[p];
        valid[p] = !ADJ(g, J[P[0]], J[P[1]]) && !ADJ(g, J[P[2]], J[P[3]]) && !ADJ(g, J[P[4]], J[P[5]]);
        if (!valid[p]) continue;
        nvalid++;
        for (int m = 0; m < DEG; m++) del_edge(g, q, J[m]);
        for (int m = 0; m < 6; m += 2) add_edge(g, J[P[m]], J[P[m + 1]]);
        dloc[p] = eloc(g, VL, nv) - e0;
        if (dfull) dfull[p] = efull(g) - E0;
        for (int m = 0; m < 6; m += 2) del_edge(g, J[P[m]], J[P[m + 1]]);
        for (int m = 0; m < DEG; m++) add_edge(g, q, J[m]);
    }
    unmark(nv);
    return nvalid;
}
static int gibbs_pick(Rng *r, const long *dE, const int *valid, int n, double T, double *wsel)
{
    long emin = 0;
    int first = 1;
    for (int p = 0; p < n; p++)
        if (!valid || valid[p])
            if (first || dE[p] < emin) {
                emin = dE[p];
                first = 0;
            }
    double Z = 0;
    for (int p = 0; p < n; p++)
        if (!valid || valid[p]) Z += exp(-2.0 * (double)(dE[p] - emin) / T);
    double u = runif(r) * Z, c = 0;
    int pick = -1;
    for (int p = 0; p < n; p++)
        if (!valid || valid[p]) {
            pick = p;
            c += exp(-2.0 * (double)(dE[p] - emin) / T);
            if (u < c) break;
        }
    *wsel = exp(-2.0 * (double)(dE[pick] - emin) / T) / Z;
    return pick;
}
static int make_gap(G *g, Rng *r)
{
    for (int att = 0; att < 1000; att++) {
        int q = rint_(r, g->n), valid[15];
        long dl[15];
        if (g->deg[q] != DEG) continue;
        int nvalid = gap_candidates(g, q, valid, dl, NULL);
        if (!nvalid) continue;
        double w;
        int p = gibbs_pick(r, dl, valid, 15, TREG, &w);
        long dmin = dl[p];
        for (int m = 0; m < 15; m++)
            if (valid[m] && dl[m] < dmin) dmin = dl[m];
        int J[DEG];
        memcpy(J, g->nb + q * DEG, sizeof J);
        for (int m = 0; m < DEG; m++) del_edge(g, q, J[m]);
        for (int m = 0; m < 6; m += 2) add_edge(g, J[PM[p][m]], J[PM[p][m + 1]]);
        g_drop(g, q);
        printf("GAP attempts=%d valid=%d dE_chosen=%.6g dE_min=%.6g w_chosen=%.6g\n", att + 1, nvalid, 2.0 * dl[p], 2.0 * dmin, w);
        return 1;
    }
    printf("GAP failed\n");
    return 0;
}
/* candidate triples of pairwise disjoint contacts with both ends within distance 2 of p */
typedef struct {
    int e[3];
} Tri;
static int excess_candidates(G *g, int p, int **eu_o, int **ev_o, Tri **tr_o, long **dl_o, long **df_o, int want_full)
{
    int n = g->n, *dist = malloc(sizeof(int) * n), *que = malloc(sizeof(int) * n);
    bfs(g, &p, 1, dist, que);
    int ne = 0, *eu = malloc(sizeof(int) * 3 * (size_t)n), *ev = malloc(sizeof(int) * 3 * (size_t)n);
    for (int u = 0; u < n; u++)
        if (dist[u] >= 0 && dist[u] <= 2)
            for (int k = 0; k < g->deg[u]; k++) {
                int v = g->nb[u * DEG + k];
                if (v > u && dist[v] >= 0 && dist[v] <= 2) {
                    eu[ne] = u;
                    ev[ne] = v;
                    ne++;
                }
            }
    int cap = 1024, nt = 0;
    Tri *tr = malloc(sizeof(Tri) * cap);
    for (int a = 0; a < ne; a++)
        for (int b = a + 1; b < ne; b++) {
            if (eu[a] == eu[b] || eu[a] == ev[b] || ev[a] == eu[b] || ev[a] == ev[b]) continue;
            for (int c = b + 1; c < ne; c++) {
                int x = eu[c], y = ev[c];
                if (x == eu[a] || x == ev[a] || y == eu[a] || y == ev[a] || x == eu[b] || x == ev[b] || y == eu[b] || y == ev[b]) continue;
                if (nt == cap) {
                    cap *= 2;
                    tr = realloc(tr, sizeof(Tri) * cap);
                }
                tr[nt].e[0] = a;
                tr[nt].e[1] = b;
                tr[nt].e[2] = c;
                nt++;
            }
        }
    long *dl = malloc(sizeof(long) * (nt ? nt : 1)), *df = want_full ? malloc(sizeof(long) * (nt ? nt : 1)) : NULL;
    int q = g->n;
    g->n++; /* the new quantum, no contacts yet */
    long E0 = want_full ? efull(g) : 0;
    for (int t = 0; t < nt; t++) {
        int V0[7], nv;
        V0[0] = q;
        for (int m = 0; m < 3; m++) {
            V0[1 + 2 * m] = eu[tr[t].e[m]];
            V0[2 + 2 * m] = ev[tr[t].e[m]];
        }
        nv = build_V1(g, V0, 7);
        long e0 = eloc(g, VL, nv);
        for (int m = 0; m < 3; m++) {
            int u = eu[tr[t].e[m]], v = ev[tr[t].e[m]];
            del_edge(g, u, v);
            add_edge(g, u, q);
            add_edge(g, q, v);
        }
        dl[t] = eloc(g, VL, nv) - e0;
        if (want_full) df[t] = efull(g) - E0;
        for (int m = 0; m < 3; m++) {
            int u = eu[tr[t].e[m]], v = ev[tr[t].e[m]];
            del_edge(g, u, q);
            del_edge(g, q, v);
            add_edge(g, u, v);
        }
        unmark(nv);
    }
    g->n--;
    free(dist);
    free(que);
    *eu_o = eu;
    *ev_o = ev;
    *tr_o = tr;
    *dl_o = dl;
    if (df_o) *df_o = df;
    return nt;
}
static int make_excess(G *g, Rng *r)
{
    for (int att = 0; att < 1000; att++) {
        int p = rint_(r, g->n), *eu, *ev;
        Tri *tr;
        long *dl;
        int nt = excess_candidates(g, p, &eu, &ev, &tr, &dl, NULL, 0);
        if (!nt) {
            free(eu);
            free(ev);
            free(tr);
            free(dl);
            continue;
        }
        double w;
        int t = gibbs_pick(r, dl, NULL, nt, TREG, &w);
        long dmin = dl[0];
        for (int m = 0; m < nt; m++)
            if (dl[m] < dmin) dmin = dl[m];
        int q = g->n;
        g->n++;
        for (int m = 0; m < 3; m++) {
            int u = eu[tr[t].e[m]], v = ev[tr[t].e[m]];
            del_edge(g, u, v);
            add_edge(g, u, q);
            add_edge(g, q, v);
        }
        printf("EXCESS attempts=%d triples=%d dE_chosen=%.6g dE_min=%.6g w_chosen=%.6g\n", att + 1, nt, 2.0 * dl[t], 2.0 * dmin, w);
        free(eu);
        free(ev);
        free(tr);
        free(dl);
        return 1;
    }
    printf("EXCESS failed\n");
    return 0;
}

/* ---------------- the registered run ---------------- */
static const char *setname(void)
{
    static char s[8];
    int c = 0;
    if (OPT_A) s[c++] = 'a';
    if (OPT_B) s[c++] = 'b';
    if (OPT_C) s[c++] = 'c';
    s[c] = 0;
    return c ? s : "none";
}
/* rigidity of the flat medium: every flip, its energy change (in t) */
static void rigid(G *g)
{
    long nflip = 0, nle0 = 0, nmin = 0;
    double dmin = INFINITY;
    for (int j = 0; j < g->n; j++) {
        int Nj[DEG], dj = g->deg[j];
        memcpy(Nj, g->nb + j * DEG, sizeof Nj);
        for (int a = 0; a < dj; a++) {
            int k = Nj[a], Nk[DEG], dk = g->deg[k];
            memcpy(Nk, g->nb + k * DEG, sizeof Nk);
            for (int b = 0; b < dj; b++) {
                int i = Nj[b];
                if (i == k) continue;
                for (int c = 0; c < dk; c++) {
                    int l = Nk[c];
                    if (l == j || i == l || ADJ(g, i, k) || ADJ(g, j, l)) continue;
                    Flip f = {i, j, k, l};
                    double dE = anyopt() ? 2.0 * (double)flip_de(g, &f) : 0.0;
                    nflip++;
                    if (dE <= 0) nle0++;
                    if (dE < dmin) {
                        dmin = dE;
                        nmin = 0;
                    }
                    if (dE == dmin) nmin++;
                }
            }
        }
    }
    printf("RIGID flips=%ld min_dE=%.6g at_min=%ld dE_le0=%ld\n", nflip, dmin, nmin, nle0);
}
static void stage(G *g, Rng *r, const char *sec, const char *name, int ns, double T, int free_)
{
    NPROP = NACC = 0;
    sweeps(g, r, ns, T, free_);
    printf("STAGE sec=%s stage=%s sweeps=%d T=%s%g acc=%.8g\n", sec, name, ns, free_ ? "inf:" : "", T, NPROP ? (double)NACC / NPROP : 0.0);
    fflush(stdout);
}
static void window(G *g, Rng *r, const char *sec, int mode)
{
    double *S = calloc(RMAX, sizeof(double));
    int smax = 0, ns = MEAS_SWEEPS / SAMPLE_EVERY;
    long np = 0, na = 0;
    for (int k = 1; k <= ns; k++) {
        NPROP = NACC = 0;
        sweeps(g, r, SAMPLE_EVERY, TREG, 0);
        np += NPROP;
        na += NACC;
        measure(g, r, sec, k, mode, S, &smax);
    }
    printf("STAGE sec=%s stage=meas sweeps=%d T=%g acc=%.8g\n", sec, MEAS_SWEEPS, TREG, np ? (double)na / np : 0.0);
    for (int d = 0; d <= smax; d++) printf("SPHERE sec=%s r=%d S=%.6f\n", sec, d, S[d] / ((double)g->n * ns));
    free(S);
}
static void run(int L, uint64_t seed)
{
    int N = L * L * L, cap = N + 2;
    scratch(cap);
    G *g = g_new(cap), *s0 = g_new(cap), *h = g_new(cap);
    Rng r;
    rng_seed(&r, seed);
    printf("RUN set=%s L=%d N=%d seed=%llu T=%g eq=%d meas=%d every=%d melt=%d anneal=%d\n", setname(), L, N, (unsigned long long)seed, TREG, EQ_SWEEPS, MEAS_SWEEPS,
           SAMPLE_EVERY, MELT_SWEEPS, ANNEAL_SWEEPS);
    build_cubic(g, L);
    printf("FLAT E=%.6g\n", 2.0 * (double)efull(g));
    rigid(g);
    fflush(stdout);
    /* M: the medium */
    stage(g, &r, "M", "eq", EQ_SWEEPS, TREG, 0);
    g_copy(s0, g);
    window(g, &r, "M", 1);
    printf("CHECK sec=M ok=%d\n", check(g, N));
    /* G: one gap */
    g_copy(h, s0);
    if (make_gap(h, &r)) {
        stage(h, &r, "G", "eq", EQ_SWEEPS, TREG, 0);
        window(h, &r, "G", 2);
        printf("CHECK sec=G ok=%d\n", check(h, N - 1));
    }
    /* X: one excess */
    g_copy(h, s0);
    if (make_excess(h, &r)) {
        stage(h, &r, "X", "eq", EQ_SWEEPS, TREG, 0);
        window(h, &r, "X", 2);
        printf("CHECK sec=X ok=%d\n", check(h, N + 1));
    }
    /* A: assembly from the melt */
    build_cubic(h, L);
    stage(h, &r, "A", "melt", MELT_SWEEPS, 1, 1);
    {
        double S0[RMAX] = {0};
        int sm0 = 0;
        measure(h, &r, "A0", 0, 0, S0, &sm0); /* the melt itself, one reading */
    }
    for (int s = 0; s < NANNEAL; s++) stage(h, &r, "A", "anneal", ANNEAL_SWEEPS, ANNEAL_T[s], 0);
    window(h, &r, "A", 0);
    printf("CHECK sec=A ok=%d\n", check(h, N));
    printf("DONE\n");
    g_free(g);
    g_free(s0);
    g_free(h);
}

/* ---------------- tool gates (PREREG §6) ---------------- */
static void setopts(int m)
{
    OPT_A = (m >> 2) & 1;
    OPT_B = (m >> 1) & 1;
    OPT_C = m & 1;
}
static int ctype9(const G *g) /* cycle type of the complement of a 6-regular graph on 9 quanta */
{
    int nn[9][2], c[9] = {0};
    for (int i = 0; i < 9; i++) {
        int m = 0;
        for (int j = 0; j < 9; j++)
            if (j != i && !ADJ(g, i, j)) {
                if (m == 2) return -1;
                nn[i][m++] = j;
            }
        if (m != 2) return -1;
    }
    int seen[9] = {0}, len[9], nl = 0;
    for (int s = 0; s < 9; s++) {
        if (seen[s]) continue;
        int prev = -1, cur = s, l = 0;
        do {
            seen[cur] = 1;
            l++;
            int nx = nn[cur][0] != prev ? nn[cur][0] : nn[cur][1];
            prev = cur;
            cur = nx;
        } while (cur != s && l < 10);
        len[nl++] = l;
    }
    for (int i = 0; i < nl; i++) c[len[i] - 1]++;
    if (nl == 1 && c[8] == 1) return 0;
    if (nl == 2 && c[2] == 1 && c[5] == 1) return 1;
    if (nl == 2 && c[3] == 1 && c[4] == 1) return 2;
    if (nl == 3 && c[2] == 3) return 3;
    return -1;
}
static void build_k9_minus(G *g, const int *cyc, int ncyc, const int *lens)
{
    int F[9][9] = {{0}}, off = 0;
    for (int c = 0; c < ncyc; c++) {
        for (int m = 0; m < lens[c]; m++) {
            int u = cyc[off + m], v = cyc[off + (m + 1) % lens[c]];
            F[u][v] = F[v][u] = 1;
        }
        off += lens[c];
    }
    g_reset(g, 9);
    for (int i = 0; i < 9; i++)
        for (int j = i + 1; j < 9; j++)
            if (!F[i][j]) add_edge(g, i, j);
}
static int cmp_int(const void *a, const void *b) { return (*(const int *)a > *(const int *)b) - (*(const int *)a < *(const int *)b); }
static int code8(const G *g)
{
    int code = 0, pw = 1;
    for (int i = 0; i < 8; i++) {
        int p = -1;
        for (int j = 0; j < 8; j++)
            if (j != i && !ADJ(g, i, j)) p = j;
        code += p * pw;
        pw *= 8;
    }
    return code;
}
static int t0(void)
{
    int pass = 1;
    init_pm();
    scratch(4200);
    /* T0.1 invariants and T0.2 local = full energy changes, L = 6, every set, flat and melted starts */
    {
        long bad = 0, tested = 0;
        int inv = 0;
        for (int m = 0; m < 8; m++) {
            setopts(m);
            for (int start = 0; start < 2; start++) {
                G *g = g_new(218);
                build_cubic(g, 6);
                Rng r;
                rng_seed(&r, 101 + 2 * m + start);
                if (start) sweeps(g, &r, 200, 1, 1);
                for (int p = 0; p < 4000; p++) {
                    Flip f;
                    if (!propose(g, &r, &f)) continue;
                    long E0 = efull(g), dl = flip_de(g, &f);
                    flip_apply(g, &f);
                    long E1 = efull(g);
                    tested++;
                    if (E1 - E0 != dl) bad++;
                    if (runif(&r) < 0.5) flip_revert(g, &f);
                }
                for (int p = 0; p < 10; p++) {
                    NPROP = NACC = 0;
                    for (int s = 0; s < 10000; s++) mc_step(g, &r, TREG, 0);
                    if (!check(g, 216)) inv++;
                }
                g_free(g);
            }
        }
        int ok = bad == 0 && inv == 0;
        printf("T0.1 invariants (six distinct contacts, no self contact, adjacency consistent, one piece) violations=%d %s\n", inv, inv ? "FAIL" : "PASS");
        printf("T0.2 local dE = full dE: tested=%ld mismatches=%ld %s\n", tested, bad, bad ? "FAIL" : "PASS");
        pass &= ok;
    }
    /* T0.3 detailed balance */
    {
        /* N = 8: every state is K8 minus a perfect matching (105, all alike): uniform target */
        int codes[105], nc = 0;
        {
            /* enumerate the perfect matchings of 8: partner of each quantum, coded base 8 */
            int part[8], used[8];
            for (int x1 = 1; x1 < 8; x1++) {
                memset(used, 0, sizeof used);
                used[0] = used[x1] = 1;
                int s1 = 0;
                while (used[s1]) s1++;
                for (int x2 = s1 + 1; x2 < 8; x2++) {
                    if (used[x2]) continue;
                    used[s1] = used[x2] = 1;
                    int s2 = 0;
                    while (used[s2]) s2++;
                    for (int x3 = s2 + 1; x3 < 8; x3++) {
                        if (used[x3]) continue;
                        used[s2] = used[x3] = 1;
                        int s3 = 0;
                        while (used[s3]) s3++;
                        int x4 = s3 + 1;
                        while (used[x4]) x4++;
                        part[0] = x1;
                        part[x1] = 0;
                        part[s1] = x2;
                        part[x2] = s1;
                        part[s2] = x3;
                        part[x3] = s2;
                        part[s3] = x4;
                        part[x4] = s3;
                        int code = 0, pw = 1;
                        for (int i = 0; i < 8; i++) {
                            code += part[i] * pw;
                            pw *= 8;
                        }
                        codes[nc++] = code;
                        used[s2] = used[x3] = 0;
                    }
                    used[s1] = used[x2] = 0;
                }
            }
        }
        qsort(codes, (size_t)nc, sizeof(int), cmp_int);
        setopts(7);
        G *g = g_new(10);
        g_reset(g, 8);
        for (int i = 0; i < 8; i++)
            for (int j = i + 1; j < 8; j++)
                if (!(j == i + 1 && i % 2 == 0)) add_edge(g, i, j);
        Rng r;
        rng_seed(&r, 7);
        long *h = calloc(105, sizeof(long)), ns = 0;
        for (long p = 0; p < 10000000; p++) {
            mc_step(g, &r, TREG, 0);
            if (p % 10 == 9) {
                int c = code8(g), *f = bsearch(&c, codes, (size_t)nc, sizeof(int), cmp_int);
                if (f) h[f - codes]++;
                ns++;
            }
        }
        double tv = 0;
        for (int i = 0; i < 105; i++) tv += fabs((double)h[i] / ns - 1.0 / 105);
        tv *= 0.5;
        int ok8 = nc == 105 && tv < 0.01;
        printf("T0.3a detailed balance N=8 (105 alike states, uniform target): states=%d TV=%.5f %s\n", nc, tv, ok8 ? "PASS" : "FAIL");
        free(h);
        g_free(g);
        pass &= ok8;
        /* N = 9: four kinds of state (complement C9, C6+C3, C5+C4, C3+C3+C3) with exact counts */
        static const int cyc[4][9] = {{0, 1, 2, 3, 4, 5, 6, 7, 8}, {0, 1, 2, 3, 4, 5, 6, 7, 8}, {0, 1, 2, 3, 4, 5, 6, 7, 8}, {0, 1, 2, 3, 4, 5, 6, 7, 8}};
        static const int lens[4][3] = {{9, 0, 0}, {6, 3, 0}, {5, 4, 0}, {3, 3, 3}}, ncy[4] = {1, 2, 2, 3};
        static const double count[4] = {20160, 5040, 4536, 280};
        int ok9 = 1;
        double tvmax = 0;
        for (int m = 0; m < 8; m++) {
            setopts(m);
            G *g9 = g_new(10);
            double Et[4];
            for (int t = 0; t < 4; t++) {
                build_k9_minus(g9, cyc[t], ncy[t], lens[t]);
                Et[t] = 2.0 * (double)efull(g9);
                if (ctype9(g9) != t) ok9 = 0;
            }
            double emin = Et[0], emax = Et[0];
            for (int t = 1; t < 4; t++) {
                if (Et[t] < emin) emin = Et[t];
                if (Et[t] > emax) emax = Et[t];
            }
            double Tt = (emax - emin) / 4 > 1 ? (emax - emin) / 4 : 1, pex[4], Z = 0;
            for (int t = 0; t < 4; t++) Z += pex[t] = count[t] * exp(-(Et[t] - emin) / Tt);
            build_k9_minus(g9, cyc[0], ncy[0], lens[0]);
            Rng r9;
            rng_seed(&r9, 900 + m);
            long hh[4] = {0}, n9 = 0;
            for (long p = 0; p < 10000000; p++) {
                mc_step(g9, &r9, Tt, 0);
                if (p % 10 == 9) {
                    int t = ctype9(g9);
                    if (t < 0) ok9 = 0;
                    else hh[t]++;
                    n9++;
                }
            }
            double tv9 = 0;
            for (int t = 0; t < 4; t++) tv9 += fabs((double)hh[t] / n9 - pex[t] / Z);
            tv9 *= 0.5;
            if (tv9 > tvmax) tvmax = tv9;
            if (tv9 >= 0.01) ok9 = 0;
            printf("T0.3b N=9 set=%s E(types)=%.0f %.0f %.0f %.0f T_test=%.4g exact=%.4f %.4f %.4f %.4f mc=%.4f %.4f %.4f %.4f TV=%.5f\n", setname(), Et[0], Et[1], Et[2],
                   Et[3], Tt, pex[0] / Z, pex[1] / Z, pex[2] / Z, pex[3] / Z, (double)hh[0] / n9, (double)hh[1] / n9, (double)hh[2] / n9, (double)hh[3] / n9, tv9);
            g_free(g9);
        }
        printf("T0.3b detailed balance N=9, all sets: max TV=%.5f %s\n", tvmax, ok9 ? "PASS" : "FAIL");
        pass &= ok9;
    }
    /* T0.4 the flat medium: zero energy, exact distance and lambda2 */
    {
        int Ls[5] = {6, 8, 10, 12, 16}, ok = 1;
        double md[5], l2[5];
        for (int a = 0; a < 5; a++) {
            int L = Ls[a], N = L * L * L, diam;
            G *g = g_new(N + 2);
            build_cubic(g, L);
            long A2, B2, T3;
            term_counts(g, &A2, &B2, &T3);
            apsp_basic(g, &md[a], &diam);
            Rng r;
            rng_seed(&r, 400 + L);
            int lok;
            l2[a] = lambda2(g, &r, &lok);
            double mex = 0.75 * L * N / (N - 1.0), lex = 2 - 2 * cos(2 * M_PI / L);
            int oka = A2 == 0 && B2 == 0 && T3 == 0 && fabs(md[a] - mex) < 1e-9 && diam == 3 * L / 2 && lok && fabs(l2[a] - lex) < 1e-8 && check(g, N);
            ok &= oka;
            printf("T0.4 flat L=%d: sum(s-4)^2=%ld sum(q-4)^2=%ld triangles=%ld meand=%.9f (exact %.9f) diam=%d (exact %d) lam2=%.12f (exact %.12f) %s\n", L, A2, B2, T3,
                   md[a], mex, diam, 3 * L / 2, l2[a], lex, oka ? "PASS" : "FAIL");
            g_free(g);
        }
        double r12 = log(4096.0 / 1728.0);
        printf("T0.4 flat readings 12->16: d_H=%.4f d_s=%.4f\n", r12 / log(md[4] / md[3]), -2 * r12 / log(l2[4] / l2[3]));
        pass &= ok;
    }
    /* T0.5 the 2-D competitor, and the instruments on known media */
    {
        int ok = 1;
        for (int L = 6; L <= 8; L += 2) {
            G *g = g_new(L * L + 2);
            build_tri(g, L);
            long A2, B2, T3;
            term_counts(g, &A2, &B2, &T3);
            int oka = A2 == 0 && B2 == 0 && T3 == 2L * L * L;
            ok &= oka;
            printf("T0.5 triangular 2-D L=%d: sum(s-4)^2=%ld sum(q-4)^2=%ld triangles=%ld (2L^2=%d) %s\n", L, A2, B2, T3, 2 * L * L, oka ? "PASS" : "FAIL");
            g_free(g);
        }
        int Lt[2] = {42, 64};
        double md[2], l2[2], nn[2];
        for (int a = 0; a < 2; a++) {
            int diam, lok;
            G *g = g_new(Lt[a] * Lt[a] + 2);
            build_tri(g, Lt[a]);
            apsp_basic(g, &md[a], &diam);
            Rng r;
            rng_seed(&r, 500 + a);
            l2[a] = lambda2(g, &r, &lok);
            nn[a] = (double)g->n;
            g_free(g);
        }
        double dH = log(nn[1] / nn[0]) / log(md[1] / md[0]), ds = -2 * log(nn[1] / nn[0]) / log(l2[1] / l2[0]);
        int ok2 = dH > 1.8 && dH < 2.2 && ds > 1.8 && ds < 2.2;
        printf("T0.5 instruments on the 2-D medium (N=%d, %d): d_H=%.4f d_s=%.4f %s\n", (int)nn[0], (int)nn[1], dH, ds, ok2 ? "PASS" : "FAIL");
        ok &= ok2;
        int Lr[2] = {12, 16};
        setopts(0);
        for (int a = 0; a < 2; a++) {
            int diam, lok, N = Lr[a] * Lr[a] * Lr[a];
            G *g = g_new(N + 2);
            build_cubic(g, Lr[a]);
            Rng r;
            rng_seed(&r, 600 + a);
            sweeps(g, &r, 2000, 1, 1);
            apsp_basic(g, &md[a], &diam);
            l2[a] = lambda2(g, &r, &lok);
            nn[a] = N;
            printf("T0.5 random six-contact medium N=%d: meand=%.4f diam=%d lam2=%.6f\n", N, md[a], diam, l2[a]);
            g_free(g);
        }
        dH = log(nn[1] / nn[0]) / log(md[1] / md[0]);
        ds = -2 * log(nn[1] / nn[0]) / log(l2[1] / l2[0]);
        int ok3 = dH > 4 && (ds > 4 || ds < 0 || !isfinite(ds));
        printf("T0.5 instruments on the random medium: d_H=%.4f d_s=%.4f (must read far from 3) %s\n", dH, ds, ok3 ? "PASS" : "FAIL");
        ok &= ok3;
        pass &= ok;
    }
    /* T0.6 making a gap and an excess: candidate counts, local = full, invariants after */
    {
        int ok = 1;
        for (int m = 0; m < 8; m++) {
            setopts(m);
            G *g = g_new(218);
            build_cubic(g, 6);
            int valid[15], nvalid;
            long dl[15], dfu[15];
            nvalid = gap_candidates(g, 0, valid, dl, dfu);
            int bad = 0;
            for (int p = 0; p < 15; p++)
                if (valid[p] && dl[p] != dfu[p]) bad++;
            int *eu, *ev;
            Tri *tr;
            long *xl, *xf;
            int nt = excess_candidates(g, 0, &eu, &ev, &tr, &xl, &xf, 1), xbad = 0;
            for (int t = 0; t < nt; t++)
                if (xl[t] != xf[t]) xbad++;
            free(eu);
            free(ev);
            free(tr);
            free(xl);
            free(xf);
            int intact = check(g, 216);
            Rng r;
            rng_seed(&r, 700 + m);
            G *h = g_new(218);
            g_copy(h, g);
            int okg = make_gap(h, &r) && check(h, 215);
            g_copy(h, g);
            int okx = make_excess(h, &r) && check(h, 217);
            int oka = nvalid == 15 && bad == 0 && xbad == 0 && intact && okg && okx;
            ok &= oka;
            printf("T0.6 set=%s gap: valid pairings=%d (expect 15) local!=full=%d; excess: triples=%d local!=full=%d; medium intact=%d; after gap ok=%d after excess ok=%d %s\n",
                   setname(), nvalid, bad, nt, xbad, intact, okg, okx, oka ? "PASS" : "FAIL");
            g_free(g);
            g_free(h);
        }
        pass &= ok;
    }
    /* T0.7 Lanczos lambda2 against dense */
    {
        int ok = 1;
        setopts(0);
        int Ls[3] = {6, 6, 8}, melt[3] = {0, 1, 1};
        for (int a = 0; a < 3; a++) {
            int N = Ls[a] * Ls[a] * Ls[a], lok;
            G *g = g_new(N + 2);
            build_cubic(g, Ls[a]);
            Rng r;
            rng_seed(&r, 800 + a);
            if (melt[a]) sweeps(g, &r, 200, 1, 1);
            double lz = lambda2(g, &r, &lok), ld = lambda2_dense(g);
            int oka = lok && fabs(lz - ld) < 1e-8;
            ok &= oka;
            printf("T0.7 N=%d %s: lanczos=%.12f dense=%.12f %s\n", N, melt[a] ? "melted" : "flat", lz, ld, oka ? "PASS" : "FAIL");
            g_free(g);
        }
        pass &= ok;
    }
    /* T0.8 a perfect flat medium of N-1 quanta exists (a sheared closure, C_n(1,a,b)) */
    {
        int ok = 1;
        setopts(7);
        for (int L = 6; L <= 8; L += 2) {
            int n = L * L * L - 1, found = 0;
            G *g = g_new(n + 2);
            for (int a = 2; a < n / 2 && !found; a++)
                for (int b = a + 1; b < (n + 1) / 2 && !found; b++) {
                    if (2 * b == n || 2 * a == n) continue;
                    build_circ(g, n, a, b);
                    if (efull(g) == 0 && check(g, n)) {
                        double md;
                        int diam;
                        apsp_basic(g, &md, &diam);
                        printf("T0.8 L=%d: n=L^3-1=%d flat with zero energy at (1,%d,%d): meand=%.4f diam=%d\n", L, n, a, b, md, diam);
                        found = 1;
                    }
                }
            if (!found) printf("T0.8 L=%d: no zero-energy circulant found FAIL\n", L);
            ok &= found;
            g_free(g);
        }
        pass &= ok;
        printf("T0.8 %s\n", ok ? "PASS" : "FAIL");
    }
    printf("T0 %s\n", pass ? "PASS" : "FAIL");
    return pass ? 0 : 1;
}

int main(int argc, char **argv)
{
    const char *mode = "run";
    int L = 6, a = 0, b = 0, c = 0;
    long long seed = -1;
    for (int i = 1; i < argc; i++) {
        if (!strncmp(argv[i], "mode=", 5)) mode = argv[i] + 5;
        else if (!strncmp(argv[i], "L=", 2)) L = atoi(argv[i] + 2);
        else if (!strncmp(argv[i], "a=", 2)) a = atoi(argv[i] + 2);
        else if (!strncmp(argv[i], "b=", 2)) b = atoi(argv[i] + 2);
        else if (!strncmp(argv[i], "c=", 2)) c = atoi(argv[i] + 2);
        else if (!strncmp(argv[i], "seed=", 5)) seed = atoll(argv[i] + 5);
        else {
            fprintf(stderr, "unknown argument %s\n", argv[i]);
            return 2;
        }
    }
    if (!strcmp(mode, "t0")) return t0();
    if (strcmp(mode, "run") || L < 5) {
        fprintf(stderr, "usage: ./geo3d mode=run L=12 a=1 b=0 c=1 | mode=t0\n");
        return 2;
    }
    OPT_A = a != 0;
    OPT_B = b != 0;
    OPT_C = c != 0;
    init_pm();
    /* the registered seed (PREREG §4) */
    uint64_t s = seed >= 0 ? (uint64_t)seed : (uint64_t)(1000003LL * L + 7LL * (4 * OPT_A + 2 * OPT_B + OPT_C) + 1);
    run(L, s);
    return 0;
}
