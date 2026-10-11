/* bx.c — v103 shared tool (lanes a and b): bosonic medium on a periodic ring,
 * exact diagonalization by restarted Lanczos (complex, OpenMP).
 *
 * Cells i = 0..L-1 (a ring: no edge).  Each holds bosonic quanta in M modes
 * (M = 1, or M = 2: the two-plane field, SU(2)), at most C per cell.
 *   H = -t sum_i sum_mu ( e^{i theta} b+_{i+1,mu} b_{i,mu} + h.c. )      theta = phi/L
 *       + (U/2) sum_i N_i (N_i - 1)
 *       - Jx sum_i A+_{i,i+1} A_{i,i+1},  A_ij = sum_{mu nu} P_{mu nu} b_{i mu} b_{j nu}
 *   P = 1 (M = 1) or the antisymmetric eps (M = 2).  t = 1 sets the unit.
 *
 *   ./bx mode=t0
 *   ./bx mode=a1 L=8 M=2 C=2 nb=2 Jx=5 [U=0]        lane a1 ledger
 *   ./bx mode=b1 L=4 M=1 C=0 nb=4 U=4               lane b1 (C=0: no cap)
 */
#include <complex.h>
#include <math.h>
#include <omp.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef double complex cplx;

/* ---------------- model ---------------- */
static int L = 6, M = 1, C = 2;
static double tt = 1, U = 0, Jx = 0, phi = 0;

/* local states: tuples (n_0..n_{M-1}) with sum <= C */
static int NLOC;
static int ln[1024][2], lN[1024], lS[1024]; /* occupations, total, 2Sz */
static int ldec[1024][2], linc[1024][2];    /* index after removing / adding a quantum of mode mu, or -1 */
static int lidx(int n0, int n1)
{
    for (int s = 0; s < NLOC; s++)
        if (ln[s][0] == n0 && (M == 1 || ln[s][1] == n1)) return s;
    return -1;
}
static void build_local(void)
{
    NLOC = 0;
    for (int tot = 0; tot <= C; tot++) {
        if (M == 1) {
            ln[NLOC][0] = tot;
            ln[NLOC][1] = 0;
            lN[NLOC] = tot;
            lS[NLOC] = 0;
            NLOC++;
        } else
            for (int u = tot; u >= 0; u--) {
                ln[NLOC][0] = u;
                ln[NLOC][1] = tot - u;
                lN[NLOC] = tot;
                lS[NLOC] = u - (tot - u);
                NLOC++;
            }
        if (NLOC > 1000) {
            fprintf(stderr, "local dimension too large\n");
            exit(1);
        }
    }
    for (int s = 0; s < NLOC; s++)
        for (int mu = 0; mu < M; mu++) {
            int a = ln[s][0], b = ln[s][1];
            if (mu == 0) {
                ldec[s][0] = a > 0 ? lidx(a - 1, b) : -1;
                linc[s][0] = lN[s] < C ? lidx(a + 1, b) : -1;
            } else {
                ldec[s][1] = b > 0 ? lidx(a, b - 1) : -1;
                linc[s][1] = lN[s] < C ? lidx(a, b + 1) : -1;
            }
        }
}

/* ---------------- sectors ---------------- */
typedef struct {
    int N, S2;          /* total quanta, 2Sz */
    long D;
    uint64_t *code;     /* sorted */
    long *hash;
    long hmask;
} Sector;

static uint64_t powl_[24];
static inline int cell_of(uint64_t code, int i) { return (int)((code / powl_[i]) % (uint64_t)NLOC); }
static inline uint64_t set_cell(uint64_t code, int i, int olds, int news)
{
    return code + ((uint64_t)news - (uint64_t)olds) * powl_[i];
}
static uint64_t hmix(uint64_t x)
{
    x ^= x >> 33;
    x *= 0xff51afd7ed558ccdULL;
    x ^= x >> 33;
    x *= 0xc4ceb9fe1a85ec53ULL;
    return x ^ (x >> 33);
}
static long sfind(const Sector *S, uint64_t c)
{
    uint64_t h = hmix(c) & S->hmask;
    while (S->hash[h] >= 0) {
        if (S->code[S->hash[h]] == c) return S->hash[h];
        h = (h + 1) & S->hmask;
    }
    return -1;
}
static long cap_codes;
static void enum_cells(Sector *S, int i, int Nrem, int Srem, uint64_t code)
{
    if (i == L) {
        if (Nrem == 0 && (M == 1 || Srem == 0)) {
            if (S->D == cap_codes) cap_codes *= 2, S->code = realloc(S->code, sizeof(uint64_t) * cap_codes);
            S->code[S->D++] = code;
        }
        return;
    }
    int cellsleft = L - i - 1;
    for (int s = 0; s < NLOC; s++) {
        if (lN[s] > Nrem) continue;
        if (Nrem - lN[s] > (long)C * cellsleft) continue;
        enum_cells(S, i + 1, Nrem - lN[s], Srem - lS[s], code + (uint64_t)s * powl_[i]);
    }
}
static int cmpu(const void *a, const void *b)
{
    uint64_t x = *(const uint64_t *)a, y = *(const uint64_t *)b;
    return x < y ? -1 : x > y;
}
static Sector *make_sector(int N, int S2)
{
    Sector *S = calloc(1, sizeof(Sector));
    S->N = N;
    S->S2 = S2;
    cap_codes = 1024;
    S->code = malloc(sizeof(uint64_t) * cap_codes);
    if (N >= 0 && N <= (long)C * L) enum_cells(S, 0, N, S2, 0);
    qsort(S->code, S->D, sizeof(uint64_t), cmpu);
    long hc = 1;
    while (hc < 2 * S->D + 2) hc <<= 1;
    S->hash = malloc(sizeof(long) * hc);
    for (long h = 0; h < hc; h++) S->hash[h] = -1;
    S->hmask = hc - 1;
    for (long k = 0; k < S->D; k++) {
        uint64_t h = hmix(S->code[k]) & S->hmask;
        while (S->hash[h] >= 0) h = (h + 1) & S->hmask;
        S->hash[h] = k;
    }
    return S;
}
static void free_sector(Sector *S)
{
    free(S->code);
    free(S->hash);
    free(S);
}

/* ---------------- H ---------------- */
static double Pm[2][2];
static void Happly(const Sector *S, const cplx *x, cplx *y)
{
    cplx eth = cexp(I * phi / L);
#pragma omp parallel for schedule(dynamic, 1024)
    for (long r = 0; r < S->D; r++) {
        /* gather: y[r] = sum_c H[r][c] x[c].  H is Hermitian: H[r][c] = conj(H[c][r]);
           we generate the terms H[c][r] (from config r) and use their conjugates. */
        uint64_t code = S->code[r];
        int st[24];
        for (int i = 0; i < L; i++) st[i] = cell_of(code, i);
        cplx acc = 0;
        double diag = 0;
        for (int i = 0; i < L; i++) diag += 0.5 * U * lN[st[i]] * (lN[st[i]] - 1);
        for (int i = 0; i < L; i++) {
            int j = (i + 1) % L;
            /* hops: b+_j b_i (i -> j, amplitude -t e^{i th}) and b+_i b_j (amplitude -t e^{-i th}) */
            for (int mu = 0; mu < M; mu++) {
                int si = st[i], sj = st[j];
                if (ldec[si][mu] >= 0 && linc[sj][mu] >= 0) {
                    double amp = sqrt((double)ln[si][mu]) * sqrt((double)ln[sj][mu] + 1);
                    uint64_t c2 = set_cell(set_cell(code, i, si, ldec[si][mu]), j, sj, linc[sj][mu]);
                    long cidx = sfind(S, c2);
                    /* <c2| H |r> = -t e^{i th} amp ; so H[r][c2] = conj = -t e^{-i th} amp */
                    if (cidx >= 0) acc += -tt * conj(eth) * amp * x[cidx];
                }
                if (ldec[sj][mu] >= 0 && linc[si][mu] >= 0) {
                    double amp = sqrt((double)ln[sj][mu]) * sqrt((double)ln[si][mu] + 1);
                    uint64_t c2 = set_cell(set_cell(code, j, sj, ldec[sj][mu]), i, si, linc[si][mu]);
                    long cidx = sfind(S, c2);
                    if (cidx >= 0) acc += -tt * eth * amp * x[cidx];
                }
            }
            /* exchange -Jx A+A on link (i,j) */
            if (Jx != 0) {
                for (int mu = 0; mu < M; mu++)
                    for (int nu = 0; nu < M; nu++) {
                        double p1 = Pm[mu][nu];
                        if (p1 == 0) continue;
                        int si = st[i], sj = st[j];
                        if (ldec[si][mu] < 0 || ldec[sj][nu] < 0) continue;
                        double a1 = sqrt((double)ln[si][mu]) * sqrt((double)ln[sj][nu]);
                        int ti = ldec[si][mu], tj = ldec[sj][nu];
                        for (int mp = 0; mp < M; mp++)
                            for (int np = 0; np < M; np++) {
                                double p2 = Pm[mp][np];
                                if (p2 == 0) continue;
                                if (linc[ti][mp] < 0 || linc[tj][np] < 0) continue;
                                double a2 = sqrt((double)ln[ti][mp] + 1) * sqrt((double)ln[tj][np] + 1);
                                uint64_t c2 = set_cell(set_cell(code, i, si, linc[ti][mp]), j, sj, linc[tj][np]);
                                long cidx = sfind(S, c2);
                                if (cidx >= 0) acc += -Jx * p1 * p2 * a1 * a2 * x[cidx];
                            }
                    }
            }
        }
        y[r] = acc + diag * x[r];
    }
}

/* ---------------- vectors ---------------- */
static double vnorm(const cplx *x, long D)
{
    double s = 0;
#pragma omp parallel for reduction(+ : s)
    for (long k = 0; k < D; k++) s += creal(x[k] * conj(x[k]));
    return sqrt(s);
}
static cplx vdot(const cplx *x, const cplx *y, long D) /* <x|y> */
{
    double re = 0, im = 0;
#pragma omp parallel for reduction(+ : re, im)
    for (long k = 0; k < D; k++) {
        cplx z = conj(x[k]) * y[k];
        re += creal(z);
        im += cimag(z);
    }
    return re + I * im;
}
static uint64_t rng = 88172645463325252ULL;
static double urand(void)
{
    rng ^= rng << 13;
    rng ^= rng >> 7;
    rng ^= rng << 17;
    return (rng >> 11) * (1.0 / 9007199254740992.0);
}

/* LAPACK tridiagonal eigensolver */
extern void dstev_(char *jobz, int *n, double *d, double *e, double *z, int *ldz, double *work, int *info);
static void tri_eig(int m, const double *al, const double *be, double *ev, double *Zm)
{
    double *d = malloc(sizeof(double) * m), *e = malloc(sizeof(double) * m), *w = malloc(sizeof(double) * 2 * m);
    memcpy(d, al, sizeof(double) * m);
    for (int k = 0; k < m - 1; k++) e[k] = be[k];
    char jz = 'V';
    int info, ld = m;
    dstev_(&jz, &m, d, e, Zm, &ld, w, &info);
    memcpy(ev, d, sizeof(double) * m);
    free(d);
    free(e);
    free(w);
}

/* Lanczos with full reorthogonalization: V[0] = v0 (normalized, orthogonal to defl).
   Fills al, be and the basis V[0..steps-1]; returns steps. */
static int lanczos_full(const Sector *S, const cplx *v0, int m, double *al, double *be, cplx **V, cplx **defl, int ndefl)
{
    long D = S->D;
    cplx *w = malloc(sizeof(cplx) * D);
    memcpy(V[0], v0, sizeof(cplx) * D);
    int steps = 0;
    for (int j = 0; j < m; j++) {
        Happly(S, V[j], w);
        double a = creal(vdot(V[j], w, D));
        al[j] = a;
        steps = j + 1;
        for (int pass = 0; pass < 2; pass++) {
            for (int q = 0; q <= j; q++) {
                cplx c = vdot(V[q], w, D);
                cplx *vq = V[q];
#pragma omp parallel for
                for (long k = 0; k < D; k++) w[k] -= c * vq[k];
            }
            for (int q = 0; q < ndefl; q++) {
                cplx c = vdot(defl[q], w, D);
                cplx *dq = defl[q];
#pragma omp parallel for
                for (long k = 0; k < D; k++) w[k] -= c * dq[k];
            }
        }
        double b = vnorm(w, D);
        /* breakdown, relative to the scale of H v: the Krylov space is exhausted and what is
           left is roundoff, which must not be followed into other symmetry sectors */
        if (j == m - 1 || b < 1e-9 * (1 + fabs(a))) break;
        be[j] = b;
        cplx *vn = V[j + 1];
#pragma omp parallel for
        for (long k = 0; k < D; k++) vn[k] = w[k] / b;
    }
    free(w);
    return steps;
}
static cplx **alloc_basis(int m, long D)
{
    cplx **V = malloc(sizeof(cplx *) * m);
    for (int j = 0; j < m; j++) V[j] = malloc(sizeof(cplx) * D);
    return V;
}
static void free_basis(cplx **V, int m)
{
    for (int j = 0; j < m; j++) free(V[j]);
    free(V);
}
extern void zheev_(char *jobz, char *uplo, int *n, cplx *a, int *lda, double *w, cplx *work, int *lwork, double *rwork, int *info);
/* dense H of a small sector: eigenvalues ascending, eigenvectors in columns */
static double *dense_eig(const Sector *S, cplx **vecs)
{
    int D = (int)S->D;
    cplx *A = calloc((size_t)D * D, sizeof(cplx)), *e = calloc(D, sizeof(cplx)), *col = malloc(sizeof(cplx) * D);
    for (int c = 0; c < D; c++) {
        memset(e, 0, sizeof(cplx) * D);
        e[c] = 1;
        Happly(S, e, col);
        for (int r = 0; r < D; r++) A[(size_t)c * D + r] = col[r];
    }
    double *w = malloc(sizeof(double) * D), *rw = malloc(sizeof(double) * 3 * D);
    int lw = 4 * D, info;
    cplx *work = malloc(sizeof(cplx) * lw);
    char jz = vecs ? 'V' : 'N', up = 'U';
    zheev_(&jz, &up, &D, A, &D, w, work, &lw, rw, &info);
    if (vecs) *vecs = A;
    else free(A);
    free(e);
    free(col);
    free(rw);
    free(work);
    return w;
}
#define DENSE_MAX 2500

/* Lowest eigenpair in S, orthogonal to defl.  Small sectors: dense; otherwise restarted Lanczos
   with full reorthogonalization.  Returns E; psi normalized; resid = ||H psi - E psi||. */
static double lowest(const Sector *S, const cplx *start, cplx *psi, double *resid, cplx **defl, int ndefl)
{
    long D = S->D;
    cplx *w = malloc(sizeof(cplx) * D);
    double E;
    if (D <= DENSE_MAX) {
        cplx *A;
        double *ev = dense_eig(S, &A);
        /* the lowest eigenvector orthogonal to defl: project the dense eigenvectors in order */
        int pick = 0;
        for (int c = 0; c < D; c++) {
            double ov = 0;
            for (int q = 0; q < ndefl; q++) {
                cplx z = 0;
                for (long k = 0; k < D; k++) z += conj(defl[q][k]) * A[(size_t)c * D + k];
                ov += creal(z * conj(z));
            }
            if (ov < 0.5) {
                pick = c;
                break;
            }
        }
        for (long k = 0; k < D; k++) psi[k] = A[(size_t)pick * D + k];
        free(A);
        free(ev);
    } else {
        cplx *v0 = malloc(sizeof(cplx) * D);
        if (start) memcpy(v0, start, sizeof(cplx) * D);
        else
            for (long k = 0; k < D; k++) v0[k] = (urand() - 0.5) + I * (urand() - 0.5);
        int m = 100;
        cplx **V = alloc_basis(m, D);
        double *al = malloc(sizeof(double) * m), *be = calloc(m, sizeof(double)), *ev = malloc(sizeof(double) * m),
               *Zm = malloc(sizeof(double) * m * m);
        for (int it = 0; it < 200; it++) {
            for (int q = 0; q < ndefl; q++) {
                cplx c = vdot(defl[q], v0, D);
                for (long k = 0; k < D; k++) v0[k] -= c * defl[q][k];
            }
            double nn = vnorm(v0, D);
            for (long k = 0; k < D; k++) v0[k] /= nn;
            int mm = lanczos_full(S, v0, m, al, be, V, defl, ndefl);
            tri_eig(mm, al, be, ev, Zm);
            memset(psi, 0, sizeof(cplx) * D);
            for (int j = 0; j < mm; j++) {
                double cj = Zm[j];
                cplx *vj = V[j];
#pragma omp parallel for
                for (long k = 0; k < D; k++) psi[k] += cj * vj[k];
            }
            double n2 = vnorm(psi, D);
            for (long k = 0; k < D; k++) psi[k] /= n2;
            Happly(S, psi, w);
            E = creal(vdot(psi, w, D));
            double r = 0;
            for (long k = 0; k < D; k++) r += creal((w[k] - E * psi[k]) * conj(w[k] - E * psi[k]));
            *resid = sqrt(r);
            if (getenv("BXDBG")) fprintf(stderr, "  restart %d: D=%ld mm=%d E=%.12f resid=%.3e\n", it, D, mm, E, *resid);
            if (*resid < 1e-9) break;
            memcpy(v0, psi, sizeof(cplx) * D);
        }
        free_basis(V, m);
        free(v0);
        free(al);
        free(be);
        free(ev);
        free(Zm);
    }
    Happly(S, psi, w);
    E = creal(vdot(psi, w, D));
    double r = 0;
    for (long k = 0; k < D; k++) r += creal((w[k] - E * psi[k]) * conj(w[k] - E * psi[k]));
    *resid = sqrt(r);
    free(w);
    return E;
}

/* Spectral Lanczos from vector x (not normalized): lowest pole energy and its weight
   (fraction of |x|^2), with full reorthogonalization; conv = change of the lowest pole between
   m/2 and m steps.  If psi != NULL, the lowest pole's state is returned. */
static double spectral(const Sector *S, const cplx *x, double *weight, double *conv, cplx *psi)
{
    long D = S->D;
    double nx = vnorm(x, D);
    if (nx < 1e-14) {
        *weight = 0;
        *conv = 0;
        return NAN;
    }
    int m = 160;
    if (m > D) m = (int)D;
    cplx *v0 = malloc(sizeof(cplx) * D);
    for (long k = 0; k < D; k++) v0[k] = x[k] / nx;
    cplx **V = alloc_basis(m, D);
    double *al = malloc(sizeof(double) * m), *be = calloc(m, sizeof(double)), *ev = malloc(sizeof(double) * m),
           *Zm = malloc(sizeof(double) * m * m);
    int mm = lanczos_full(S, v0, m, al, be, V, NULL, 0);
    tri_eig(mm, al, be, ev, Zm);
    /* the lowest pole: the lowest Ritz level that carries weight */
    int j0 = 0;
    while (j0 < mm - 1 && Zm[(size_t)j0 * mm] * Zm[(size_t)j0 * mm] < 1e-10) j0++;
    double E = ev[j0], w0 = 0;
    for (int j = 0; j < mm; j++)
        if (fabs(ev[j] - E) < 1e-7 * (1 + fabs(E))) w0 += Zm[(size_t)j * mm] * Zm[(size_t)j * mm];
    int mh = mm / 2 > 2 ? mm / 2 : mm;
    double *ev2 = malloc(sizeof(double) * mh), *Z2 = malloc(sizeof(double) * mh * mh);
    tri_eig(mh, al, be, ev2, Z2);
    int j2 = 0;
    while (j2 < mh - 1 && Z2[(size_t)j2 * mh] * Z2[(size_t)j2 * mh] < 1e-10) j2++;
    *conv = mm < m ? 0 : fabs(ev2[j2] - E); /* an exhausted Krylov space is exact */
    *weight = w0;
    if (psi) {
        memset(psi, 0, sizeof(cplx) * D);
        for (int j = 0; j < mm; j++)
            for (long k = 0; k < D; k++) psi[k] += Zm[(size_t)j0 * mm + j] * V[j][k];
        double nn = vnorm(psi, D);
        for (long k = 0; k < D; k++) psi[k] /= nn;
    }
    free_basis(V, m);
    free(v0);
    free(al);
    free(be);
    free(ev);
    free(Zm);
    free(ev2);
    free(Z2);
    return E;
}

/* ---------------- operators between sectors ---------------- */
/* y (in T) = O x (in S), O = sum_j f_j b_{j mu} (dir=-1) or sum_j f_j b+_{j mu} (dir=+1) */
static void apply_b(const Sector *S, const cplx *x, const Sector *T, cplx *y, int mu, int dir, const cplx *f)
{
    memset(y, 0, sizeof(cplx) * T->D);
    for (long r = 0; r < S->D; r++) {
        if (x[r] == 0) continue;
        uint64_t code = S->code[r];
        for (int j = 0; j < L; j++) {
            int s = cell_of(code, j);
            int s2 = dir < 0 ? ldec[s][mu] : linc[s][mu];
            if (s2 < 0) continue;
            double amp = dir < 0 ? sqrt((double)ln[s][mu]) : sqrt((double)ln[s][mu] + 1);
            long t = sfind(T, set_cell(code, j, s, s2));
            if (t >= 0) y[t] += f[j] * amp * x[r];
        }
    }
}
/* S+ = sum_j b+_{j0} b_{j1}: S (2Sz) -> T (2Sz+2) */
static void apply_Splus(const Sector *S, const cplx *x, const Sector *T, cplx *y)
{
    memset(y, 0, sizeof(cplx) * T->D);
    for (long r = 0; r < S->D; r++) {
        if (x[r] == 0) continue;
        uint64_t code = S->code[r];
        for (int j = 0; j < L; j++) {
            int s = cell_of(code, j);
            int a = ln[s][0], b = ln[s][1];
            if (b == 0) continue;
            int s2 = lidx(a + 1, b - 1);
            double amp = sqrt((double)b) * sqrt((double)a + 1);
            long t = sfind(T, set_cell(code, j, s, s2));
            if (t >= 0) y[t] += amp * x[r];
        }
    }
}
static double S2expect(const Sector *S, const cplx *psi)
{
    if (M != 2) return 0;
    Sector *T = make_sector(S->N, S->S2 + 2);
    double sz = S->S2 / 2.0, nrm2 = 0;
    if (T->D > 0) {
        cplx *y = malloc(sizeof(cplx) * T->D);
        apply_Splus(S, psi, T, y);
        nrm2 = vnorm(y, T->D);
        nrm2 *= nrm2;
        free(y);
    }
    free_sector(T);
    return nrm2 + sz * (sz + 1);
}

/* ---------------- diagonal observables ---------------- */
static double expect_diag(const Sector *S, const cplx *psi, double (*f)(const int *))
{
    double s = 0;
#pragma omp parallel for reduction(+ : s)
    for (long r = 0; r < S->D; r++) {
        int st[24];
        for (int i = 0; i < L; i++) st[i] = lN[cell_of(S->code[r], i)];
        s += creal(psi[r] * conj(psi[r])) * f(st);
    }
    return s;
}
static int g_nb; /* medium filling */
static double f_double(const int *n)
{
    int c = 0;
    for (int i = 0; i < L; i++) c += n[i] <= g_nb - 2;
    return c;
}
static double f_N0sq(const int *n)
{
    double s = 0;
    for (int i = 0; i < L; i++) s += (double)n[i] * n[i];
    return s / L;
}
static double f_half(const int *n) /* translation-averaged N_half */
{
    double s = 0;
    for (int o = 0; o < L; o++) {
        int h = 0;
        for (int k = 0; k < L / 2; k++) h += n[(o + k) % L];
        s += h;
    }
    return s / L;
}
static double f_half2(const int *n)
{
    double s = 0;
    for (int o = 0; o < L; o++) {
        int h = 0;
        for (int k = 0; k < L / 2; k++) h += n[(o + k) % L];
        s += (double)h * h;
    }
    return s / L;
}
static double f_fluct(const int *n)
{
    double s = 0;
    for (int i = 0; i < L; i++) s += (double)(n[i] - g_nb) * (n[i] - g_nb);
    return s;
}
/* g1(r) = sum_mu <b+_{i mu} b_{i+r mu}>, translation-averaged */
static double g1(const Sector *S, const cplx *psi, int rr)
{
    double tot = 0;
    for (int mu = 0; mu < M; mu++) {
        double re = 0;
#pragma omp parallel for reduction(+ : re)
        for (long r = 0; r < S->D; r++) {
            uint64_t code = S->code[r];
            for (int i = 0; i < L; i++) {
                int j = (i + rr) % L;
                /* <psi| b+_i b_j |psi>: term from config r: b+_i b_j |r> */
                int si = cell_of(code, i), sj = cell_of(code, j);
                if (i == j) {
                    re += creal(psi[r] * conj(psi[r])) * ln[si][mu];
                    continue;
                }
                if (ldec[sj][mu] < 0 || linc[si][mu] < 0) continue;
                double amp = sqrt((double)ln[sj][mu]) * sqrt((double)ln[si][mu] + 1);
                long t = sfind(S, set_cell(set_cell(code, j, sj, ldec[sj][mu]), i, si, linc[si][mu]));
                if (t >= 0) re += creal(conj(psi[t]) * amp * psi[r]);
            }
        }
        tot += re / L;
    }
    return tot;
}

/* ---------------- full ED for gates ---------------- */
static double full_lowest(const Sector *S)
{
    double *w = dense_eig(S, NULL), E = w[0];
    free(w);
    return E;
}
/* Lanczos-only lowest (forces the iterative path, for the gate) */
static double lanczos_lowest(const Sector *S, double *resid)
{
    long D = S->D;
    cplx *v0 = malloc(sizeof(cplx) * D), *psi = malloc(sizeof(cplx) * D), *w = malloc(sizeof(cplx) * D);
    for (long k = 0; k < D; k++) v0[k] = (urand() - 0.5) + I * (urand() - 0.5);
    int m = 100 < D ? 100 : (int)D;
    cplx **V = alloc_basis(m, D);
    double *al = malloc(sizeof(double) * m), *be = calloc(m, sizeof(double)), *ev = malloc(sizeof(double) * m), *Zm = malloc(sizeof(double) * m * m), E = 0;
    for (int it = 0; it < 200; it++) {
        double nn = vnorm(v0, D);
        for (long k = 0; k < D; k++) v0[k] /= nn;
        int mm = lanczos_full(S, v0, m, al, be, V, NULL, 0);
        tri_eig(mm, al, be, ev, Zm);
        memset(psi, 0, sizeof(cplx) * D);
        for (int j = 0; j < mm; j++)
            for (long k = 0; k < D; k++) psi[k] += Zm[j] * V[j][k];
        nn = vnorm(psi, D);
        for (long k = 0; k < D; k++) psi[k] /= nn;
        Happly(S, psi, w);
        E = creal(vdot(psi, w, D));
        double r = 0;
        for (long k = 0; k < D; k++) r += creal((w[k] - E * psi[k]) * conj(w[k] - E * psi[k]));
        *resid = sqrt(r);
        if (*resid < 1e-9) break;
        memcpy(v0, psi, sizeof(cplx) * D);
    }
    free_basis(V, m);
    free(v0);
    free(psi);
    free(w);
    free(al);
    free(be);
    free(ev);
    free(Zm);
    return E;
}
static double herm_check(const Sector *S)
{
    long D = S->D;
    cplx *x = malloc(sizeof(cplx) * D), *y = malloc(sizeof(cplx) * D), *Hx = malloc(sizeof(cplx) * D), *Hy = malloc(sizeof(cplx) * D);
    for (long k = 0; k < D; k++) x[k] = urand() - 0.5 + I * (urand() - 0.5), y[k] = urand() - 0.5 + I * (urand() - 0.5);
    Happly(S, x, Hx);
    Happly(S, y, Hy);
    cplx a = vdot(x, Hy, D), b = conj(vdot(y, Hx, D));
    double r = cabs(a - b) / (cabs(a) + 1e-300);
    free(x);
    free(y);
    free(Hx);
    free(Hy);
    return r;
}
/* SU(2): || H S+ x - S+ H x || / ||S+ H x|| */
static double su2_check(int N, int S2)
{
    Sector *S = make_sector(N, S2), *T = make_sector(N, S2 + 2);
    if (S->D == 0 || T->D == 0) return 0;
    cplx *x = malloc(sizeof(cplx) * S->D), *Hx = malloc(sizeof(cplx) * S->D), *a = malloc(sizeof(cplx) * T->D), *b = malloc(sizeof(cplx) * T->D),
         *c = malloc(sizeof(cplx) * T->D);
    for (long k = 0; k < S->D; k++) x[k] = urand() - 0.5 + I * (urand() - 0.5);
    Happly(S, x, Hx);
    apply_Splus(S, Hx, T, a);
    apply_Splus(S, x, T, c);
    Happly(T, c, b);
    double num = 0;
    for (long k = 0; k < T->D; k++) num += creal((a[k] - b[k]) * conj(a[k] - b[k]));
    double r = sqrt(num) / (vnorm(a, T->D) + 1e-300);
    free(x);
    free(Hx);
    free(a);
    free(b);
    free(c);
    free_sector(S);
    free_sector(T);
    return r;
}

static void setup(void)
{
    build_local();
    powl_[0] = 1;
    for (int i = 1; i < 24; i++) powl_[i] = powl_[i - 1] * (uint64_t)NLOC;
    if (M == 1) Pm[0][0] = 1, Pm[0][1] = Pm[1][0] = Pm[1][1] = 0;
    else Pm[0][0] = Pm[1][1] = 0, Pm[0][1] = 1, Pm[1][0] = -1;
}

static const char *arg(int argc, char **argv, const char *k, const char *def)
{
    size_t n = strlen(k);
    for (int i = 1; i < argc; i++)
        if (!strncmp(argv[i], k, n) && argv[i][n] == '=') return argv[i] + n + 1;
    return def;
}

/* ---------------- gates ---------------- */
static int t0(void)
{
    int fails = 0;
    /* 1. hard-core bosons (M=1, C=1) on a ring = free fermions (boundary twist by parity) */
    for (int l = 4; l <= 10; l += 2)
        for (int n = 1; n < l; n++)
            for (int pk = 0; pk < 2; pk++) {
                L = l;
                M = 1;
                C = 1;
                U = 0;
                Jx = 0;
                phi = pk ? 1.3 : 0;
                setup();
                Sector *S = make_sector(n, 0);
                cplx *psi = malloc(sizeof(cplx) * S->D);
                double res, E = lowest(S, NULL, psi, &res, NULL, 0);
                /* free fermions: k = 2pi(m + delta + phi/2pi)/L, delta = 0 (n odd), 1/2 (n even) */
                double delta = (n % 2) ? 0 : 0.5, ek[64];
                for (int m = 0; m < l; m++) ek[m] = -2 * cos(2 * M_PI * (m + delta + phi / (2 * M_PI)) / l);
                for (int a = 0; a < l; a++)
                    for (int b = a + 1; b < l; b++)
                        if (ek[b] < ek[a]) {
                            double t = ek[a];
                            ek[a] = ek[b];
                            ek[b] = t;
                        }
                double Ef = 0;
                for (int m = 0; m < n; m++) Ef += ek[m];
                if (fabs(E - Ef) > 1e-8 || res > 1e-8) {
                    printf("T0 FAIL hard-core L=%d n=%d phi=%.2f E=%.12f free-fermion=%.12f resid=%.2e\n", l, n, phi, E, Ef, res);
                    fails++;
                }
                free(psi);
                free_sector(S);
            }
    printf("T0 1 hard-core bosons = free fermions (L=4..10, all n, phi=0 and 1.3): %s\n", fails ? "FAIL" : "PASS");
    /* 2. free bosons, no cap: E0 = -2 t N cos(phi/L) */
    int f2 = 0;
    for (int n = 1; n <= 8; n++) {
        L = 5;
        M = 1;
        C = n;
        U = 0;
        Jx = 0;
        phi = 0.7;
        setup();
        Sector *S = make_sector(n, 0);
        cplx *psi = malloc(sizeof(cplx) * S->D);
        double res, E = lowest(S, NULL, psi, &res, NULL, 0);
        double Ex = -2 * n * cos(phi / L);
        if (fabs(E - Ex) > 1e-8) {
            printf("T0 FAIL free bosons n=%d E=%.12f exact=%.12f\n", n, E, Ex);
            f2++;
        }
        free(psi);
        free_sector(S);
    }
    printf("T0 2 free bosons E0 = -2tN cos(phi/L) (L=5, N=1..8): %s\n", f2 ? "FAIL" : "PASS");
    fails += f2;
    /* 3, 4, 5: Hermiticity, SU(2), Lanczos vs full ED for every rule-set shape */
    int f3 = 0;
    double hmax = 0, smax = 0, emax = 0;
    int Ms[2] = {1, 2};
    for (int im = 0; im < 2; im++)
        for (int cap = 2; cap <= 3; cap++)
            for (int jx = 0; jx <= 1; jx++)
                for (int uu = 0; uu <= 1; uu++) {
                    L = 4;
                    M = Ms[im];
                    C = cap;
                    Jx = jx ? 5 : 0;
                    U = uu ? 4 : 0;
                    phi = 0.9;
                    setup();
                    for (int n = 5; n <= 9; n++) {
                        Sector *S = make_sector(n, M == 2 ? (n % 2) : 0);
                        if (S->D == 0) {
                            free_sector(S);
                            continue;
                        }
                        double h = herm_check(S);
                        if (h > hmax) hmax = h;
                        if (M == 2) {
                            double s = su2_check(n, M == 2 ? (n % 2) - 2 : 0);
                            if (s > smax) smax = s;
                        }
                        if (S->D <= 1500) {
                            double res, E = lanczos_lowest(S, &res), Ef = full_lowest(S);
                            if (fabs(E - Ef) > emax) emax = fabs(E - Ef);
                        }
                        free_sector(S);
                    }
                }
    /* 6. spectral(): lowest pole and weight against dense diagonalization, for b_k|G> */
    double spmax = 0, swmax = 0;
    for (int im = 0; im < 2; im++)
        for (int cap = 2; cap <= 3; cap++) {
            L = 4;
            M = Ms[im];
            C = cap;
            Jx = 5;
            U = 0;
            phi = 0;
            setup();
            int N = 8;
            Sector *G = make_sector(N, 0), *V = make_sector(N - 1, M == 2 ? -1 : 0);
            cplx *g = malloc(sizeof(cplx) * G->D), *x = malloc(sizeof(cplx) * V->D), f[16];
            double res;
            lowest(G, NULL, g, &res, NULL, 0);
            for (int mk = 0; mk < L; mk++) {
                for (int j = 0; j < L; j++) f[j] = cexp(-I * 2 * M_PI * mk * j / L) / sqrt(L);
                apply_b(G, g, V, x, 0, -1, f);
                double nx = vnorm(x, V->D);
                if (nx < 1e-10) continue;
                double w, c, E = spectral(V, x, &w, &c, NULL);
                cplx *A;
                double *ev = dense_eig(V, &A);
                /* exact: lowest eigenvalue level with overlap > 1e-12, and its summed weight */
                double Ex = NAN, Wx = 0;
                for (int col = 0; col < V->D; col++) {
                    cplx z = 0;
                    for (long k = 0; k < V->D; k++) z += conj(A[(size_t)col * V->D + k]) * x[k];
                    double ww = creal(z * conj(z)) / (nx * nx);
                    if (isnan(Ex) && ww > 1e-12) Ex = ev[col];
                    if (!isnan(Ex) && fabs(ev[col] - Ex) < 1e-9) Wx += ww;
                }
                if (getenv("BXDBG")) fprintf(stderr, "spec M=%d C=%d k=%d D=%ld E=%.10f Ex=%.10f w=%.6f Wx=%.6f conv=%.1e\n", M, C, mk, V->D, E, Ex, w, Wx, c);
                if (fabs(E - Ex) > spmax) spmax = fabs(E - Ex);
                if (fabs(w - Wx) > swmax) swmax = fabs(w - Wx);
                free(A);
                free(ev);
            }
            free(g);
            free(x);
            free_sector(G);
            free_sector(V);
        }
    printf("T0 6 spectral lowest pole vs dense: max |dE| %.2e, max |dweight| %.2e: %s\n", spmax, swmax, spmax < 1e-8 && swmax < 1e-6 ? "PASS" : "FAIL");
    if (spmax >= 1e-8 || swmax >= 1e-6) f3 = 1;
    if (hmax > 1e-12 || smax > 1e-12 || emax > 1e-9) f3 = 1;
    printf("T0 3-5 Hermiticity max %.2e, SU(2) [H,S+] max %.2e, Lanczos vs full ED max |dE| %.2e: %s\n", hmax, smax, emax, f3 ? "FAIL" : "PASS");
    fails += f3;
    printf("T0 %s\n", fails ? "FAIL" : "PASS");
    return fails;
}

/* ---------------- lane a1 ---------------- */
static void run_a1(int nb)
{
    g_nb = nb;
    int N = nb * L;
    int s2v = M == 2 ? -1 : 0, s2e = M == 2 ? 1 : 0; /* remove an up quantum: 2Sz -1; add an up quantum: +1 */
    printf("A1 model L=%d M=%d C=%d nb=%d Jx=%g U=%g\n", L, M, C, nb, Jx, U);
    /* medium */
    phi = 0;
    Sector *G = make_sector(N, 0);
    cplx *g = malloc(sizeof(cplx) * G->D), *g2 = malloc(sizeof(cplx) * G->D);
    double res, EG = lowest(G, NULL, g, &res, NULL, 0);
    cplx *dfl[1] = {g};
    double res2, EG2 = G->D > 1 ? lowest(G, NULL, g2, &res2, dfl, 1) : INFINITY;
    double S2G = S2expect(G, g);
    printf("MEDIUM D=%ld E=%.12f resid=%.2e gap=%.6e S2=%.6f fluct=%.6f\n", G->D, EG, res, EG2 - EG, S2G, expect_diag(G, g, f_fluct));
    /* void band: b_{k, up}|G> for every k */
    Sector *V = make_sector(N - 1, s2v);
    Sector *X = make_sector(N + 1, s2e);
    cplx *x = malloc(sizeof(cplx) * (V->D > X->D ? V->D : X->D) + 16);
    cplx *f = malloc(sizeof(cplx) * L);
    double Ev[64], Zv[64], Ee[64], Ze[64], cv = 0, ce = 0;
    double Emin = INFINITY;
    int kmin = 0;
    for (int m = 0; m < L; m++) {
        double k = 2 * M_PI * m / L, w, c;
        for (int j = 0; j < L; j++) f[j] = cexp(-I * k * j) / sqrt(L);
        apply_b(G, g, V, x, 0, -1, f);
        double nn = vnorm(x, V->D);
        Ev[m] = spectral(V, x, &w, &c, NULL);
        Zv[m] = w;
        if (c > cv) cv = c;
        (void)nn;
        if (Ev[m] < Emin) Emin = Ev[m], kmin = m;
        if (X->D > 0) {
            apply_b(G, g, X, x, 0, +1, f);
            if (vnorm(x, X->D) > 1e-12) {
                Ee[m] = spectral(X, x, &w, &c, NULL);
                Ze[m] = w;
                if (c > ce) ce = c;
            } else
                Ee[m] = NAN, Ze[m] = 0;
        } else
            Ee[m] = NAN, Ze[m] = 0;
    }
    double Emax = -INFINITY;
    for (int m = 0; m < L; m++)
        if (Ev[m] > Emax) Emax = Ev[m];
    printf("VOID band:");
    for (int m = 0; m < L; m++) printf(" k%d E=%.9f Z=%.6f;", m, Ev[m] - EG, Zv[m]);
    printf("\nVOID W=%.9f Zmin=%.6f kmin=%d conv=%.2e\n", Emax - Emin, Zv[kmin], kmin, cv);
    /* void spin: the lowest pole state at kmin */
    {
        double k = 2 * M_PI * kmin / L, w, c;
        for (int j = 0; j < L; j++) f[j] = cexp(-I * k * j) / sqrt(L);
        apply_b(G, g, V, x, 0, -1, f);
        cplx *ps = malloc(sizeof(cplx) * V->D);
        spectral(V, x, &w, &c, ps);
        printf("VOID S2=%.6f (medium S2=%.6f)\n", S2expect(V, ps), S2G);
        free(ps);
    }
    /* global lowest of the one-void sector, and Aharonov-Bohm */
    cplx *pv = malloc(sizeof(cplx) * V->D), *pg = malloc(sizeof(cplx) * G->D);
    double r1, r0, Ev0 = lowest(V, NULL, pv, &r1, NULL, 0);
    phi = M_PI;
    double Gp = lowest(G, NULL, pg, &r0, NULL, 0), r3, Vp = lowest(V, NULL, pv, &r3, NULL, 0);
    phi = 0;
    printf("AB E1(0)-E0(0)=%.12f E1(pi)-E0(pi)=%.12f dAB=%.6e resid=%.1e,%.1e,%.1e\n", Ev0 - EG, Vp - Gp, (Vp - Gp) - (Ev0 - EG), r1, r0, r3);
    /* excess */
    if (X->D > 0) {
        double Exmin = INFINITY;
        int km = 0;
        printf("EXCESS band:");
        for (int m = 0; m < L; m++) {
            printf(" k%d E=%.9f Z=%.6f;", m, Ee[m] - EG, Ze[m]);
            if (Ee[m] < Exmin) Exmin = Ee[m], km = m;
        }
        printf("\nEXCESS Zmin=%.6f kmin=%d conv=%.2e\n", Ze[km], km, ce);
        double thr = INFINITY;
        for (int m = 0; m < L; m++)
            if (!isnan(Ee[m]) && Ev[m] - EG + Ee[m] - EG < thr) thr = Ev[m] - EG + Ee[m] - EG;
        /* breathing: O|G> with O = sum_i (N_i - nb)^2 - <.>, k = 0, spin 0 */
        double fl = expect_diag(G, g, f_fluct);
        cplx *ob = malloc(sizeof(cplx) * G->D);
        for (long r = 0; r < G->D; r++) {
            int st[24];
            for (int i = 0; i < L; i++) st[i] = lN[cell_of(G->code[r], i)];
            ob[r] = (f_fluct(st) - fl) * g[r];
        }
        double w, c;
        /* project out the ground state */
        cplx ov = vdot(g, ob, G->D);
        for (long r = 0; r < G->D; r++) ob[r] -= ov * g[r];
        double Eb = spectral(G, ob, &w, &c, NULL);
        printf("BREATH E=%.9f weight=%.6f conv=%.2e threshold=%.9f R=%.6f\n", Eb - EG, w, c, thr, (Eb - EG) / thr);
        free(ob);
    } else {
        printf("EXCESS none (cap = filling)\n");
        printf("BREATH none (the saturated medium cannot breathe: O|G> = 0)\n");
    }
    /* two voids */
    Sector *W2 = make_sector(N - 2, 0);
    cplx *p2 = malloc(sizeof(cplx) * W2->D);
    double r2, E2 = lowest(W2, NULL, p2, &r2, NULL, 0);
    double pd = expect_diag(W2, p2, f_double);
    printf("TWOVOID D=%ld E=%.12f resid=%.1e bind=%.9f Pdouble=%.6f\n", W2->D, E2 - EG, r2, (E2 - EG) - 2 * (Ev0 - EG), pd);
    printf("CHECK medium_resid=%.1e void_resid=%.1e spectral_conv=%.1e,%.1e\n", res, r1, cv, ce);
}

/* ---------------- lane b1 ---------------- */
static void run_b1(int nb)
{
    g_nb = nb;
    int N = nb * L;
    phi = 0;
    printf("B1 model L=%d M=%d C=%d nb=%d U=%g\n", L, M, C, nb, U);
    Sector *G = make_sector(N, 0), *V = make_sector(N - 1, M == 2 ? -1 : 0);
    cplx *g = malloc(sizeof(cplx) * G->D), *v = malloc(sizeof(cplx) * V->D);
    double r0, r1, EG = lowest(G, NULL, g, &r0, NULL, 0);
    double EV = lowest(V, NULL, v, &r1, NULL, 0);
    double gr = g1(G, g, L / 2) / nb;
    double sig = sqrt(expect_diag(G, g, f_N0sq) - (double)nb * nb) / nb;
    double h0 = expect_diag(G, g, f_half), h0s = expect_diag(G, g, f_half2), h1 = expect_diag(V, v, f_half), h1s = expect_diag(V, v, f_half2);
    double var0 = h0s - h0 * h0, var1 = h1s - h1 * h1, p = h0 - h1;
    double beta = (var1 - var0) / (p * (1 - p));
    /* quasiparticle weight at the band minimum */
    cplx *f = malloc(sizeof(cplx) * L), *x = malloc(sizeof(cplx) * V->D);
    double Zmin = 0, Emin = INFINITY;
    for (int m = 0; m < L; m++) {
        double k = 2 * M_PI * m / L, w, c;
        for (int j = 0; j < L; j++) f[j] = cexp(-I * k * j) / sqrt(L);
        apply_b(G, g, V, x, 0, -1, f);
        double E = spectral(V, x, &w, &c, NULL);
        if (E < Emin) Emin = E, Zmin = w;
    }
    printf("MEDIUM D=%ld E=%.12f resid=%.1e g1(L/2)=%.6f sigmaN/nb=%.6f\n", G->D, EG, r0, gr, sig);
    printf("VOID D=%ld E=%.12f resid=%.1e p=%.6f var0=%.6f var1=%.6f beta=%.6f Z=%.6f\n", V->D, EV - EG, r1, p, var0, var1, beta, Zmin);
    printf("CHECK resid=%.1e,%.1e\n", r0, r1);
}

int main(int argc, char **argv)
{
    const char *mode = arg(argc, argv, "mode", "t0");
    L = atoi(arg(argc, argv, "L", "6"));
    M = atoi(arg(argc, argv, "M", "1"));
    C = atoi(arg(argc, argv, "C", "2"));
    U = atof(arg(argc, argv, "U", "0"));
    Jx = atof(arg(argc, argv, "Jx", "0"));
    int nb = atoi(arg(argc, argv, "nb", "2"));
    rng = strtoull(arg(argc, argv, "seed", "88172645463325252"), NULL, 10);
    if (C == 0) C = nb * L; /* no cap */
    if (!strcmp(mode, "t0")) return t0() ? 1 : 0;
    setup();
    if (!strcmp(mode, "dims")) {
        int N = nb * L;
        for (int dn = -2; dn <= 1; dn++) {
            Sector *S = make_sector(N + dn, M == 2 ? ((N + dn) % 2 ? (dn > 0 ? 1 : -1) : 0) : 0);
            printf("DIMS L=%d M=%d C=%d nb=%d N=%d D=%ld\n", L, M, C, nb, N + dn, S->D);
            free_sector(S);
        }
        return 0;
    }
    if (!strcmp(mode, "a1")) run_a1(nb);
    else if (!strcmp(mode, "b1")) run_b1(nb);
    return 0;
}
