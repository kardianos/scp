/* gap1d.c — v103 lane c2: gaps in a one-dimensional medium (../PREREG.md).
 *
 * A periodic ring of L cells; each cell is a piece of the medium holding
 * n_i quanta, 1 <= n_i <= cap (never empty, R14), full n = 2.  A gap is
 * n_i = 1, an excess n_i = 3.  Distance between two cells is counted in
 * quanta (the quanta strictly between them, the shorter way round).
 *
 * candidates (all 16 combinations):
 *   a  lack recurses: cap 2 (off: cap 3)
 *   b  tweaked links: first-moment-conserving moves T3, T4 (off: single hop)
 *   c  inertia: quantum lowest manifold (off: Gibbs e^{-E_d/t} of the sector)
 *   d  flat level: t sum (n_i - n_{i+1})^2 (off: stiff 2t sum (n_i - 2)^2)
 *
 * build: gcc -O3 -march=native -fopenmp -Wall -Wextra -o gap1d gap1d.c \
 *            -l:liblapack.so.3 -l:libopenblas.so.0 -lm
 * usage: ./gap1d mode=run L=12 a=1 b=0 c=1 d=0      (one registered run)
 *        ./gap1d mode=t0                            (tool gates)
 */
#include <math.h>
#include <omp.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

extern void dsyev_(const char *, const char *, const int *, double *, const int *, double *, double *, const int *, int *);
extern void dstev_(const char *, const int *, double *, double *, double *, const int *, double *, int *);

static int L, CAP, OPT_A, OPT_B, OPT_C, OPT_D, FORCE_LANCZOS = 0;
static const double T = 1.0;
#define DENSE_MAX 3000
#define DEGTOL 1e-8
#define MAXMAN 48

/* ---------------- sectors: configurations with sum n = N (and moment P mod L if P >= 0) ---------------- */
typedef struct {
    int N, P, PL; /* PL = L if P constrained, else 1 */
    long D;
    uint32_t *cfg; /* 2 bits per cell: n - 1 */
    long *cnt;     /* (L+1) x S x PL */
    int S;
} Sector;

static inline int cell(uint32_t c, int i) { return (int)((c >> (2 * i)) & 3u) + 1; }
static inline long CNT(const Sector *s, int i, int sum, int p) { return s->cnt[((long)i * s->S + sum) * s->PL + p]; }
static inline int pmod(long x, int m) { long r = x % m; return (int)(r < 0 ? r + m : r); }

static long rank_of(const Sector *s, const int *n)
{
    long r = 0;
    int srem = s->N, prem = s->PL > 1 ? s->P : 0;
    for (int i = 0; i < L; i++) {
        for (int u = 1; u < n[i]; u++) {
            int s2 = srem - u;
            if (s2 < 0) break;
            int p2 = s->PL > 1 ? pmod(prem - (long)i * u, L) : 0;
            r += CNT(s, i + 1, s2, p2);
        }
        srem -= n[i];
        if (s->PL > 1) prem = pmod(prem - (long)i * n[i], L);
        if (srem < 0) return -1;
    }
    return (srem == 0 && prem == 0) ? r : -1;
}

static void gen(Sector *s, int i, int srem, int prem, uint32_t c, long *k)
{
    if (i == L) {
        if (srem == 0 && prem == 0) s->cfg[(*k)++] = c;
        return;
    }
    for (int u = 1; u <= CAP; u++) {
        int s2 = srem - u;
        if (s2 < 0) break;
        int p2 = s->PL > 1 ? pmod(prem - (long)i * u, L) : 0;
        if (CNT(s, i + 1, s2, p2) == 0) continue;
        gen(s, i + 1, s2, p2, c | ((uint32_t)(u - 1) << (2 * i)), k);
    }
}

static Sector *make_sector(int N, int P)
{
    Sector *s = calloc(1, sizeof(Sector));
    s->N = N;
    s->P = P;
    s->PL = P >= 0 ? L : 1;
    s->S = CAP * L + 1;
    s->cnt = calloc((size_t)(L + 1) * s->S * s->PL, sizeof(long));
    s->cnt[((long)L * s->S + 0) * s->PL + 0] = 1;
    for (int i = L - 1; i >= 0; i--)
        for (int sum = 0; sum < s->S; sum++)
            for (int p = 0; p < s->PL; p++) {
                long v = 0;
                for (int u = 1; u <= CAP && u <= sum; u++) v += CNT(s, i + 1, sum - u, s->PL > 1 ? pmod(p - (long)i * u, L) : 0);
                s->cnt[((long)i * s->S + sum) * s->PL + p] = v;
            }
    s->D = (N >= 0 && N < s->S) ? CNT(s, 0, N, s->PL > 1 ? pmod(P, L) : 0) : 0;
    s->cfg = malloc(sizeof(uint32_t) * (s->D ? s->D : 1));
    long k = 0;
    if (s->D) gen(s, 0, N, s->PL > 1 ? pmod(P, L) : 0, 0, &k);
    if (k != s->D) {
        fprintf(stderr, "sector enumeration mismatch %ld vs %ld\n", k, s->D);
        exit(1);
    }
    return s;
}
static void free_sector(Sector *s)
{
    free(s->cfg);
    free(s->cnt);
    free(s);
}
static void decode(uint32_t c, int *n)
{
    for (int i = 0; i < L; i++) n[i] = cell(c, i);
}

/* ---------------- energies and moves ---------------- */
static double ediag(const int *n)
{
    double e = 0;
    if (OPT_D)
        for (int i = 0; i < L; i++) {
            int d = n[i] - n[(i + 1) % L];
            e += T * d * d;
        }
    else
        for (int i = 0; i < L; i++) e += 2 * T * (n[i] - 2) * (n[i] - 2);
    return e;
}

/* enumerate moves from n: calls cb(newn, amp); amp is the (negative) matrix element */
typedef void (*movecb)(const int *nn, double amp, void *ctx);
static void moves(const int *n, movecb cb, void *ctx)
{
    int nn[32];
    if (!OPT_B) {
        for (int i = 0; i < L; i++) {
            int j = (i + 1) % L;
            /* i -> j */
            if (n[i] >= 2 && n[j] < CAP) {
                memcpy(nn, n, sizeof(int) * L);
                nn[i]--, nn[j]++;
                cb(nn, -T * sqrt((double)n[i] * (n[j] + 1)), ctx);
            }
            if (n[j] >= 2 && n[i] < CAP) {
                memcpy(nn, n, sizeof(int) * L);
                nn[j]--, nn[i]++;
                cb(nn, -T * sqrt((double)n[j] * (n[i] + 1)), ctx);
            }
        }
        return;
    }
    /* T3 centred at i: (+1,-2,+1) and (-1,+2,-1) */
    for (int i = 0; i < L; i++) {
        int l = (i + L - 1) % L, r = (i + 1) % L;
        if (n[i] - 2 >= 1 && n[l] + 1 <= CAP && n[r] + 1 <= CAP) {
            memcpy(nn, n, sizeof(int) * L);
            nn[l]++, nn[r]++, nn[i] -= 2;
            cb(nn, -T * sqrt((double)(n[l] + 1) * (n[r] + 1) * n[i] * (n[i] - 1)), ctx);
        }
        if (n[l] >= 2 && n[r] >= 2 && n[i] + 2 <= CAP) {
            memcpy(nn, n, sizeof(int) * L);
            nn[l]--, nn[r]--, nn[i] += 2;
            cb(nn, -T * sqrt((double)n[l] * n[r] * (n[i] + 1) * (n[i] + 2)), ctx);
        }
    }
    /* T4 starting at i: (+1,-1,-1,+1) and (-1,+1,+1,-1) */
    for (int i = 0; i < L; i++) {
        int c0 = i, c1 = (i + 1) % L, c2 = (i + 2) % L, c3 = (i + 3) % L;
        if (n[c0] + 1 <= CAP && n[c1] >= 2 && n[c2] >= 2 && n[c3] + 1 <= CAP) {
            memcpy(nn, n, sizeof(int) * L);
            nn[c0]++, nn[c1]--, nn[c2]--, nn[c3]++;
            cb(nn, -T * sqrt((double)(n[c0] + 1) * n[c1] * n[c2] * (n[c3] + 1)), ctx);
        }
        if (n[c0] >= 2 && n[c1] + 1 <= CAP && n[c2] + 1 <= CAP && n[c3] >= 2) {
            memcpy(nn, n, sizeof(int) * L);
            nn[c0]--, nn[c1]++, nn[c2]++, nn[c3]--;
            cb(nn, -T * sqrt((double)n[c0] * (n[c1] + 1) * (n[c2] + 1) * n[c3]), ctx);
        }
    }
}

static long moment(const int *n)
{
    long p = 0;
    for (int i = 0; i < L; i++) p += (long)i * n[i];
    return pmod(p, L);
}

/* ---------------- sparse Hamiltonian ---------------- */
typedef struct {
    long D, nnz;
    long *rp;
    int32_t *ci;
    double *v, *dg;
} CSR;

typedef struct {
    const Sector *s;
    long cnt;
    int32_t *ci;
    double *v;
    long p0;
    long viol; /* moment violations */
    long src_moment;
} BCtx;
static void cb_count(const int *nn, double amp, void *vctx)
{
    (void)amp;
    BCtx *c = vctx;
    long r = rank_of(c->s, nn);
    if (r >= 0) c->cnt++;
    else if (c->s->PL > 1) c->viol++;
}
static void cb_fill(const int *nn, double amp, void *vctx)
{
    BCtx *c = vctx;
    long r = rank_of(c->s, nn);
    if (r < 0) return;
    c->ci[c->p0 + c->cnt] = (int32_t)r;
    c->v[c->p0 + c->cnt] = amp;
    c->cnt++;
}
static long g_moment_viol;
static CSR *build(const Sector *s)
{
    CSR *H = calloc(1, sizeof(CSR));
    H->D = s->D;
    H->rp = calloc(s->D + 1, sizeof(long));
    H->dg = malloc(sizeof(double) * s->D);
    long viol = 0;
#pragma omp parallel for schedule(dynamic, 4096) reduction(+ : viol)
    for (long r = 0; r < s->D; r++) {
        int n[32];
        decode(s->cfg[r], n);
        H->dg[r] = ediag(n);
        BCtx c = {s, 0, NULL, NULL, 0, 0, 0};
        moves(n, cb_count, &c);
        H->rp[r + 1] = c.cnt;
        viol += c.viol;
    }
    g_moment_viol += viol;
    for (long r = 0; r < s->D; r++) H->rp[r + 1] += H->rp[r];
    H->nnz = H->rp[s->D];
    H->ci = malloc(sizeof(int32_t) * (H->nnz ? H->nnz : 1));
    H->v = malloc(sizeof(double) * (H->nnz ? H->nnz : 1));
#pragma omp parallel for schedule(dynamic, 4096)
    for (long r = 0; r < s->D; r++) {
        int n[32];
        decode(s->cfg[r], n);
        BCtx c = {s, 0, H->ci, H->v, H->rp[r], 0, 0};
        moves(n, cb_fill, &c);
    }
    return H;
}
static void free_csr(CSR *H)
{
    free(H->rp);
    free(H->ci);
    free(H->v);
    free(H->dg);
    free(H);
}
static void matvec(const CSR *H, const double *x, double *y)
{
#pragma omp parallel for schedule(static)
    for (long r = 0; r < H->D; r++) {
        double s = H->dg[r] * x[r];
        for (long k = H->rp[r]; k < H->rp[r + 1]; k++) s += H->v[k] * x[H->ci[k]];
        y[r] = s;
    }
}
static double dot(const double *a, const double *b, long D)
{
    double s = 0;
#pragma omp parallel for reduction(+ : s) schedule(static)
    for (long k = 0; k < D; k++) s += a[k] * b[k];
    return s;
}
static void axpy(double a, const double *x, double *y, long D)
{
#pragma omp parallel for schedule(static)
    for (long k = 0; k < D; k++) y[k] += a * x[k];
}
static void scal(double a, double *x, long D)
{
#pragma omp parallel for schedule(static)
    for (long k = 0; k < D; k++) x[k] *= a;
}

static uint64_t rng_s = 0x9e3779b97f4a7c15ULL;
static double urand(void)
{
    uint64_t z = (rng_s += 0x9e3779b97f4a7c15ULL);
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    z ^= z >> 31;
    return (double)(z >> 11) / 9007199254740992.0;
}

/* ---------------- lowest manifold ---------------- */
typedef struct {
    int dim;
    double E0, Enext, resid;
    double **vec; /* dim vectors */
    int trunc, dense;
} Manifold;

static void dense_eigs(const CSR *H, double **evec_out, double **eval_out)
{
    int n = (int)H->D, lw = -1, info;
    double *A = calloc((size_t)n * n, sizeof(double)), *w = malloc(sizeof(double) * n), q;
    for (long r = 0; r < H->D; r++) {
        A[(size_t)r * n + r] += H->dg[r];
        for (long k = H->rp[r]; k < H->rp[r + 1]; k++) A[(size_t)r * n + H->ci[k]] += H->v[k];
    }
    dsyev_("V", "U", &n, A, &n, w, &q, &lw, &info);
    lw = (int)q;
    double *work = malloc(sizeof(double) * lw);
    dsyev_("V", "U", &n, A, &n, w, work, &lw, &info);
    if (info) {
        fprintf(stderr, "dsyev info %d\n", info);
        exit(1);
    }
    free(work);
    *evec_out = A;
    *eval_out = w;
}

/* restarted Lanczos with full reorthogonalisation, deflated against defl[0..nd-1] */
static double lanczos_lowest(const CSR *H, double **defl, int nd, double *psi, double *resid)
{
    long D = H->D;
    int m = 100;
    if (m > D - nd) m = (int)(D - nd);
    if (m < 1) m = 1;
    double **V = malloc(sizeof(double *) * m);
    for (int j = 0; j < m; j++) V[j] = malloc(sizeof(double) * D);
    double *w = malloc(sizeof(double) * D), *al = malloc(sizeof(double) * m), *be = malloc(sizeof(double) * m);
    double *d = malloc(sizeof(double) * m), *e = malloc(sizeof(double) * m), *Z = malloc(sizeof(double) * m * m),
           *work = malloc(sizeof(double) * (2 * m + 2));
    double E = 0;
    for (long k = 0; k < D; k++) psi[k] = urand() - 0.5;
    for (int it = 0; it < 400; it++) {
        for (int q = 0; q < nd; q++) axpy(-dot(defl[q], psi, D), defl[q], psi, D);
        scal(1 / sqrt(dot(psi, psi, D)), psi, D);
        memcpy(V[0], psi, sizeof(double) * D);
        int steps = 0;
        for (int j = 0; j < m; j++) {
            matvec(H, V[j], w);
            al[j] = dot(V[j], w, D);
            steps = j + 1;
            /* full reorthogonalisation, twice */
            for (int pass = 0; pass < 2; pass++) {
                for (int q = 0; q < nd; q++) axpy(-dot(defl[q], w, D), defl[q], w, D);
                for (int k = 0; k <= j; k++) axpy(-dot(V[k], w, D), V[k], w, D);
            }
            double b = sqrt(dot(w, w, D));
            if (j == m - 1 || b < 1e-10 * (1 + fabs(al[j]))) break;
            be[j] = b;
            for (long k = 0; k < D; k++) V[j + 1][k] = w[k] / b;
        }
        int ms = steps, info;
        for (int j = 0; j < ms; j++) d[j] = al[j], e[j] = j < ms - 1 ? be[j] : 0;
        dstev_("V", &ms, d, e, Z, &ms, work, &info);
        if (info) {
            fprintf(stderr, "dstev info %d\n", info);
            exit(1);
        }
        memset(psi, 0, sizeof(double) * D);
        for (int j = 0; j < ms; j++) axpy(Z[j], V[j], psi, D);
        scal(1 / sqrt(dot(psi, psi, D)), psi, D);
        matvec(H, psi, w);
        E = dot(psi, w, D);
        axpy(-E, psi, w, D);
        *resid = sqrt(dot(w, w, D));
        if (getenv("GDBG")) fprintf(stderr, "  it %d steps %d E %.14f resid %.3e\n", it, ms, E, *resid);
        if (*resid < 1e-10) break;
    }
    for (int j = 0; j < m; j++) free(V[j]);
    free(V);
    free(w);
    free(al);
    free(be);
    free(d);
    free(e);
    free(Z);
    free(work);
    return E;
}

static Manifold lowest_manifold(const CSR *H)
{
    Manifold M = {0};
    M.vec = malloc(sizeof(double *) * MAXMAN);
    long D = H->D;
    if (D <= DENSE_MAX && !FORCE_LANCZOS) {
        double *A, *w;
        dense_eigs(H, &A, &w);
        M.dense = 1;
        M.E0 = w[0];
        int k = 0;
        while (k < D && w[k] - w[0] < DEGTOL * (1 + fabs(w[0]))) {
            if (k < MAXMAN) {
                M.vec[k] = malloc(sizeof(double) * D);
                memcpy(M.vec[k], A + (size_t)k * D, sizeof(double) * D);
            }
            k++;
        }
        M.trunc = k > MAXMAN;
        M.dim = k > MAXMAN ? MAXMAN : k;
        M.Enext = k < D ? w[k] : NAN;
        M.resid = 0;
        free(A);
        free(w);
        return M;
    }
    double r;
    M.vec[0] = malloc(sizeof(double) * D);
    M.E0 = lanczos_lowest(H, NULL, 0, M.vec[0], &r);
    M.resid = r;
    M.dim = 1;
    while (M.dim < D) {
        double *v = malloc(sizeof(double) * D);
        double E = lanczos_lowest(H, M.vec, M.dim, v, &r);
        if (E - M.E0 < DEGTOL * (1 + fabs(M.E0))) {
            if (M.dim == MAXMAN) {
                M.trunc = 1;
                free(v);
                M.Enext = E;
                break;
            }
            M.vec[M.dim++] = v;
            if (r > M.resid) M.resid = r;
        } else {
            M.Enext = E;
            free(v);
            break;
        }
    }
    return M;
}
static void free_manifold(Manifold *M)
{
    for (int k = 0; k < M->dim; k++) free(M->vec[k]);
    free(M->vec);
}

/* ---------------- states as diagonal weights over a full (unconstrained-P) sector ---------------- */
typedef struct {
    Sector *full; /* N, P = -1 */
    double *w;    /* probability per configuration */
    double E;     /* lowest energy (c on) or free energy (c off) */
    int mdim, trunc, nsec;
    double resid, gapnext;
} State;

/* c on: lowest manifold over the needed P sectors, mapped onto the full sector */
static State state_quantum(int N)
{
    State st = {0};
    st.full = make_sector(N, -1);
    st.w = calloc(st.full->D ? st.full->D : 1, sizeof(double));
    st.gapnext = INFINITY;
    if (st.full->D == 0) {
        st.E = NAN;
        return st;
    }
    int Ps[64], nP = 0;
    if (!OPT_B) Ps[nP++] = -1;
    else {
        int shift = pmod(N, L);
        /* representatives of the orbits P -> P + N (mod L) */
        int seen[64] = {0};
        for (int p = 0; p < L; p++) {
            if (seen[p]) continue;
            Ps[nP++] = p;
            int q = p;
            do {
                seen[q] = 1;
                q = pmod(q + shift, L);
            } while (!seen[q]);
        }
    }
    double Emin = INFINITY;
    Manifold Ms[64];
    Sector *Ss[64];
    for (int k = 0; k < nP; k++) {
        Ss[k] = make_sector(N, Ps[k]);
        if (Ss[k]->D == 0) {
            Ms[k].dim = 0;
            Ms[k].E0 = INFINITY;
            Ms[k].vec = NULL;
            continue;
        }
        CSR *H = build(Ss[k]);
        Ms[k] = lowest_manifold(H);
        free_csr(H);
        if (Ms[k].E0 < Emin) Emin = Ms[k].E0;
    }
    /* orbit sizes: each representative stands for its orbit; translation maps observables onto
       each other, so an orbit's members contribute equally. Weight each orbit by its size. */
    int shift = pmod(N, L);
    double totw = 0;
    for (int k = 0; k < nP; k++) {
        if (Ms[k].dim == 0) continue;
        if (Ms[k].E0 - Emin < DEGTOL * (1 + fabs(Emin))) {
            int osz = 1;
            if (OPT_B) {
                int q = pmod(Ps[k] + shift, L);
                while (q != Ps[k]) osz++, q = pmod(q + shift, L);
            }
            st.mdim += Ms[k].dim * osz;
            st.trunc |= Ms[k].trunc;
            st.nsec += osz;
            if (Ms[k].resid > st.resid) st.resid = Ms[k].resid;
            for (int v = 0; v < Ms[k].dim; v++)
                for (long r = 0; r < Ss[k]->D; r++) {
                    int n[32];
                    decode(Ss[k]->cfg[r], n);
                    long R = rank_of(st.full, n);
                    double x = Ms[k].vec[v][r];
                    /* translation averaging is done by the observables (all translation invariant);
                       here weight the representative by the orbit size */
                    st.w[R] += osz * x * x;
                    totw += osz * x * x;
                }
            if (Ms[k].Enext - Emin < st.gapnext) st.gapnext = Ms[k].Enext - Emin;
        } else if (Ms[k].E0 - Emin < st.gapnext)
            st.gapnext = Ms[k].E0 - Emin;
    }
    for (long r = 0; r < st.full->D; r++) st.w[r] /= totw;
    st.E = Emin;
    for (int k = 0; k < nP; k++) {
        if (Ms[k].vec) free_manifold(&Ms[k]);
        free_sector(Ss[k]);
    }
    return st;
}

/* c off: Gibbs weight over the full sector */
static State state_gibbs(int N)
{
    State st = {0};
    st.full = make_sector(N, -1);
    st.w = calloc(st.full->D ? st.full->D : 1, sizeof(double));
    if (st.full->D == 0) {
        st.E = NAN;
        return st;
    }
    double emin = INFINITY;
    double *e = malloc(sizeof(double) * st.full->D);
    for (long r = 0; r < st.full->D; r++) {
        int n[32];
        decode(st.full->cfg[r], n);
        e[r] = ediag(n);
        if (e[r] < emin) emin = e[r];
    }
    double Z = 0;
    for (long r = 0; r < st.full->D; r++) Z += (st.w[r] = exp(-(e[r] - emin) / T));
    for (long r = 0; r < st.full->D; r++) st.w[r] /= Z;
    st.E = emin - T * log(Z);
    st.mdim = 0;
    free(e);
    return st;
}
static State get_state(int N) { return OPT_C ? state_quantum(N) : state_gibbs(N); }
static void free_state(State *s)
{
    free(s->w);
    free_sector(s->full);
}

/* ---------------- pair observables ---------------- */
/* minimal composition: exactly ng gaps (n=1) and nx excesses (n=3), the rest full */
typedef struct {
    double Wmin, meanD, meanDref, rho;
    int nD;
    double p[64], pref[64];
} Pair;
static Pair pair_obs(const State *st, int ng, int nx, int sameKind)
{
    Pair P = {0};
    int maxD = 2 * (L / 2 - 1);
    P.nD = maxD / 2 + 1;
    double wtot = 0;
    for (long r = 0; r < st->full->D; r++) {
        int n[32], g = 0, x = 0, pos[2], np = 0, ok = 1;
        decode(st->full->cfg[r], n);
        for (int i = 0; i < L; i++) {
            if (n[i] == 1) g++;
            else if (n[i] == 3) x++;
            if (n[i] != 2) {
                if (np < 2) pos[np] = i;
                np++;
            }
        }
        if (g != ng || x != nx || np != 2) ok = 0;
        if (!ok) continue;
        int dl = pos[1] - pos[0], dc = dl < L - dl ? dl : L - dl;
        int Dq = 2 * (dc - 1); /* quanta strictly between, the shorter way: every cell between is full */
        P.p[Dq / 2] += st->w[r];
        wtot += st->w[r];
        P.pref[Dq / 2] += 1;
    }
    (void)sameKind;
    P.Wmin = wtot;
    double rt = 0;
    for (int k = 0; k < P.nD; k++) rt += P.pref[k];
    P.meanD = P.meanDref = 0;
    for (int k = 0; k < P.nD; k++) {
        if (wtot > 0) P.p[k] /= wtot;
        if (rt > 0) P.pref[k] /= rt;
        P.meanD += 2 * k * P.p[k];
        P.meanDref += 2 * k * P.pref[k];
    }
    P.rho = (wtot > 0 && P.meanDref > 0) ? P.meanD / P.meanDref : NAN;
    return P;
}
static void print_pair(const char *tag, const Pair *P)
{
    printf("PAIR %s Wmin=%.9f meanD=%.6f meanDref=%.6f rho=%.6f\n", tag, P->Wmin, P->meanD, P->meanDref, P->rho);
    printf("PROFILE %s", tag);
    for (int k = 0; k < P->nD; k++) {
        double V = (P->p[k] > 0 && P->pref[k] > 0) ? -log(P->p[k] / P->pref[k]) : INFINITY;
        printf(" D=%d:p=%.6e,pref=%.6e,V=%.6f", 2 * k, P->p[k], P->pref[k], V);
    }
    printf("\n");
}

/* ---------------- Mach: reachability in the move graph (structure only) ---------------- */
typedef struct {
    const Sector *s;
    uint8_t *seen;
    long *queue, qt;
} RCtx;
static void cb_reach(const int *nn, double amp, void *vctx)
{
    (void)amp;
    RCtx *c = vctx;
    long r = rank_of(c->s, nn);
    if (r >= 0 && !c->seen[r]) c->seen[r] = 1, c->queue[c->qt++] = r;
}
static long fragment(const Sector *s, const int *start, uint8_t *seen)
{
    memset(seen, 0, s->D);
    long *q = malloc(sizeof(long) * s->D), qh = 0;
    RCtx c = {s, seen, q, 0};
    long r0 = rank_of(s, start);
    seen[r0] = 1;
    q[c.qt++] = r0;
    while (qh < c.qt) {
        int n[32];
        decode(s->cfg[q[qh++]], n);
        moves(n, cb_reach, &c);
    }
    free(q);
    return c.qt;
}
static int set_gaps(int *n, int g1, int g2)
{
    for (int i = 0; i < L; i++) n[i] = 2;
    if (g1 >= 0) n[g1] = 1;
    if (g2 >= 0) n[g2] = 1;
    return 0;
}
static void mach(void)
{
    int n[32];
    /* isolated gap */
    Sector *s1 = make_sector(2 * L - 1, -1);
    uint8_t *seen = malloc(s1->D);
    set_gaps(n, 0, -1);
    long f1 = fragment(s1, n, seen);
    set_gaps(n, 1, -1);
    int iso = seen[rank_of(s1, n)];
    free(seen);
    free_sector(s1);
    /* pairs */
    Sector *s2 = make_sector(2 * L - 2, -1);
    seen = malloc(s2->D);
    set_gaps(n, 0, 1);
    long fa = fragment(s2, n, seen);
    set_gaps(n, 1, 2);
    int pairmove = seen[rank_of(s2, n)];
    set_gaps(n, 0, L / 2);
    long ff = fragment(s2, n, seen);
    int join = 0;
    for (int i = 0; i < L && !join; i++) {
        set_gaps(n, i, (i + 1) % L);
        if (seen[rank_of(s2, n)]) join = 1;
    }
    free(seen);
    free_sector(s2);
    printf("MACH iso=%d iso_frag=%ld pairmove=%d adjpair_frag=%ld join=%d farpair_frag=%ld\n", iso, f1, pairmove, fa, join, ff);
}

/* ---------------- one registered run ---------------- */
static void run(void)
{
    printf("C2 L=%d a=%d b=%d c=%d d=%d cap=%d\n", L, OPT_A, OPT_B, OPT_C, OPT_D, CAP);
    int Ns[5] = {2 * L - 2, 2 * L - 1, 2 * L, 2 * L + 1, 2 * L + 2};
    State S[5];
    double maxres = 0;
    int trunc = 0;
    for (int k = 0; k < 5; k++) {
        S[k] = get_state(Ns[k]);
        printf("SECTOR N=%d D=%ld E=%.12f mdim=%d nsec=%d gapnext=%.3e resid=%.1e trunc=%d\n", Ns[k], S[k].full->D, S[k].E, S[k].mdim,
               S[k].nsec, S[k].gapnext, S[k].resid, S[k].trunc);
        if (S[k].resid > maxres) maxres = S[k].resid;
        trunc |= S[k].trunc;
    }
    double Bgg = S[0].E - 2 * S[1].E + S[2].E, Bee = S[4].E - 2 * S[3].E + S[2].E;
    if (!OPT_C) {
        /* the random process: measured against counting, F_count = -t ln D, so that a medium
           without interaction gives B = 0 (PREREG 3.2, amended before any run) */
        double Fc[5];
        for (int k = 0; k < 5; k++) Fc[k] = S[k].full->D ? -T * log((double)S[k].full->D) : NAN;
        Bgg -= Fc[0] - 2 * Fc[1] + Fc[2];
        Bee -= Fc[4] - 2 * Fc[3] + Fc[2];
    }
    printf("BIND gg=%.12f ee=%.12f\n", Bgg, Bee);
    Pair gg = pair_obs(&S[0], 2, 0, 1), ee = pair_obs(&S[4], 0, 2, 1), ge = pair_obs(&S[2], 1, 1, 0);
    print_pair("gg", &gg);
    print_pair("ee", &ee);
    print_pair("ge", &ge);
    mach();
    printf("CHECK maxresid=%.1e trunc=%d moment_viol=%ld\n", maxres, trunc, g_moment_viol);
    for (int k = 0; k < 5; k++) free_state(&S[k]);
}

/* ---------------- gates ---------------- */
static void setopts(int a, int b, int c, int d)
{
    OPT_A = a, OPT_B = b, OPT_C = c, OPT_D = d;
    CAP = a ? 2 : 3;
}
static long brute_count(int N)
{
    long total = 1, k = 0;
    for (int i = 0; i < L; i++) total *= CAP;
    for (long code = 0; code < total; code++) {
        long x = code;
        int s = 0;
        for (int i = 0; i < L; i++) s += (int)(x % CAP) + 1, x /= CAP;
        if (s == N) k++;
    }
    return k;
}
static int t0(void)
{
    int fails = 0;
    /* gate 5: sector sizes */
    for (int Lv = 8; Lv <= 12; Lv += 4)
        for (int a = 0; a <= 1; a++) {
            L = Lv;
            setopts(a, 0, 1, 0);
            for (int N = 2 * L - 2; N <= 2 * L + 2; N++) {
                Sector *s = make_sector(N, -1);
                long bc = brute_count(N);
                long sumP = 0;
                for (int p = 0; p < L; p++) {
                    Sector *sp = make_sector(N, p);
                    sumP += sp->D;
                    free_sector(sp);
                }
                int ok = bc == s->D && sumP == s->D;
                if (!ok) fails++;
                printf("T0.5 L=%d cap=%d N=%d D=%ld brute=%ld sumP=%ld %s\n", L, CAP, N, s->D, bc, sumP, ok ? "PASS" : "FAIL");
                free_sector(s);
            }
        }
    /* gates 1, 2 and 4 at L = 8, every set */
    L = 8;
    for (int m = 0; m < 16; m++) {
        int a = m >> 3 & 1, b = m >> 2 & 1, c = m >> 1 & 1, d = m & 1;
        setopts(a, b, c, d);
        double herm = 0;
        g_moment_viol = 0;
        long moved = 0;
        for (int N = 2 * L - 2; N <= 2 * L + 2; N++) {
            Sector *s = make_sector(N, -1);
            if (!s->D) {
                free_sector(s);
                continue;
            }
            CSR *H = build(s);
            moved += H->nnz;
            double *x = malloc(sizeof(double) * s->D), *y = malloc(sizeof(double) * s->D), *Hx = malloc(sizeof(double) * s->D),
                   *Hy = malloc(sizeof(double) * s->D);
            for (long k = 0; k < s->D; k++) x[k] = urand() - 0.5, y[k] = urand() - 0.5;
            matvec(H, x, Hx);
            matvec(H, y, Hy);
            double h = fabs(dot(x, Hy, s->D) - dot(Hx, y, s->D)) / (1 + fabs(dot(x, Hy, s->D)));
            if (h > herm) herm = h;
            /* moment conservation: every move of a b-on set must keep P; count violations on the full sector */
            if (b)
                for (long r = 0; r < s->D; r++)
                    for (long k = H->rp[r]; k < H->rp[r + 1]; k++) {
                        int n1[32], n2[32];
                        decode(s->cfg[r], n1);
                        decode(s->cfg[H->ci[k]], n2);
                        if (moment(n1) != moment(n2)) g_moment_viol++;
                    }
            free(x);
            free(y);
            free(Hx);
            free(Hy);
            free_csr(H);
            free_sector(s);
        }
        int ok = herm < 1e-12 && g_moment_viol == 0;
        if (!ok) fails++;
        printf("T0.1/2 set a=%d b=%d c=%d d=%d herm=%.1e moment_viol=%ld elements=%ld %s\n", a, b, c, d, herm, g_moment_viol, moved, ok ? "PASS" : "FAIL");
        /* gate 4: iterative against dense, c on sets */
        if (c) {
            double worst = 0, worstobs = 0;
            for (int N = 2 * L - 2; N <= 2 * L + 2; N++) {
                FORCE_LANCZOS = 0;
                State s1 = get_state(N);
                if (!s1.full->D) {
                    free_state(&s1);
                    continue;
                }
                FORCE_LANCZOS = 1;
                State s2 = get_state(N);
                FORCE_LANCZOS = 0;
                double de = fabs(s1.E - s2.E);
                if (de > worst) worst = de;
                if (s1.mdim != s2.mdim) worst = INFINITY;
                for (long r = 0; r < s1.full->D; r++) {
                    double dw = fabs(s1.w[r] - s2.w[r]);
                    if (dw > worstobs) worstobs = dw;
                }
                free_state(&s1);
                free_state(&s2);
            }
            int ok4 = worst < 1e-9 && worstobs < 1e-8;
            if (!ok4) fails++;
            printf("T0.4 set a=%d b=%d c=%d d=%d dense_vs_lanczos dE=%.1e dw=%.1e %s\n", a, b, c, d, worst, worstobs, ok4 ? "PASS" : "FAIL");
        }
    }
    /* gate 3: exact limit (a on, b off, d off, c on) */
    for (int Lv = 8; Lv <= 16; Lv += 4) {
        L = Lv;
        setopts(1, 0, 1, 0);
        State s1 = get_state(2 * L - 1), s2 = get_state(2 * L - 2);
        double e1 = 2 * T - 4 * T, e2 = 4 * T - 8 * T * cos(M_PI / L);
        double d1 = fabs(s1.E - e1), d2 = fabs(s2.E - e2);
        int ok = d1 < 1e-10 && d2 < 1e-10;
        if (!ok) fails++;
        printf("T0.3 L=%d one gap E=%.12f exact=%.12f  two gaps E=%.12f exact=%.12f %s\n", L, s1.E, e1, s2.E, e2, ok ? "PASS" : "FAIL");
        free_state(&s1);
        free_state(&s2);
    }
    printf("T0 %s (%d fails)\n", fails ? "FAIL" : "PASS", fails);
    return fails;
}

int main(int argc, char **argv)
{
    const char *mode = "run";
    L = 12;
    int a = 1, b = 0, c = 1, d = 0;
    for (int i = 1; i < argc; i++) {
        if (!strncmp(argv[i], "mode=", 5)) mode = argv[i] + 5;
        else if (!strncmp(argv[i], "L=", 2)) L = atoi(argv[i] + 2);
        else if (!strncmp(argv[i], "a=", 2)) a = atoi(argv[i] + 2);
        else if (!strncmp(argv[i], "b=", 2)) b = atoi(argv[i] + 2);
        else if (!strncmp(argv[i], "c=", 2)) c = atoi(argv[i] + 2);
        else if (!strncmp(argv[i], "d=", 2)) d = atoi(argv[i] + 2);
        else if (!strncmp(argv[i], "seed=", 5)) rng_s = strtoull(argv[i] + 5, NULL, 10) * 0x9e3779b97f4a7c15ULL + 1;
    }
    if (L < 4 || L > 16) {
        fprintf(stderr, "L must be 4..16\n");
        return 1;
    }
    if (!strcmp(mode, "t0")) return t0() ? 1 : 0;
    setopts(a, b, c, d);
    run();
    return 0;
}
