/* v102 kernel step K2a — lockq: an exact quantum lock inside the
 * classical kernel (v102/kernel/K2a_PREREG.md, K2a_DECIDE.md).  Included
 * by freecell.c after the kernel globals.  Every entry point returns at
 * once unless P.lockq != 0, so lockq = 0 (the default) is byte-inert.
 *
 * The lock is the leading-order effective model of the S-law fiber
 * (v102/perm/cas/lo.py, reading V):
 *   - R6 end cells (N = 1) at both ends of an open chain;
 *   - mobile bulk cells (N = 2, or N = 1 = a hole moved by the T-floor
 *     hop t, R9);
 *   - on every chain link the exchange -(g^2/Delta_type) A^+A, with
 *     A = sum P[mu,nu] b_i mu b_j nu, P = eps (3+2, M = 2) or
 *     eps (x) eps (W, M = 4; R8).
 * Energies inside the module are in the S-law's rate units; the kernel
 * carries hbar * E (hbar = A0/2pi).
 *
 * The lock is fenced: every slot touching a lock cell has lens area 0,
 * so no classical field, space or dense transport crosses it (R9, R10).
 * Its one door is the formation event at a boundary link (R10, K1-1):
 * a field atom of a neighbour converts with a lock quantum, the jump
 * J = sum_mu P[mu,0] b_i mu acts on the lock (photon mode 0: the kernel
 * field is one complex mode until K2b), and the space quantum is
 * credited to the neighbour's Es with the exact energy balance.  The
 * golden-rule rate and the band window use the exact energy change of
 * each candidate jump (K2a_DECIDE.md §G5); events are drawn as a
 * Poisson process inside each lock step.
 *
 * K2b (v102/kernel/K2b_PREREG.md): with fmodes = F > 1 the kernel field
 * carries F complex modes and lockq_M = F.  A photon in mode nu makes
 * the jump J_nu = sum_mu P[mu,nu] b_i mu; atoms are counted per mode,
 * and the event takes its atom from mode nu.  At F = 1 only nu = 0
 * exists and the module is K2a's. */
#include <complex.h>

#define LQ_NMAX 16
#define LQ_DMAX 16      /* site dimension: W mobile = 10 + 4 = 14 */
#define LQ_KMAX 64      /* Lanczos / Krylov dimension cap */

typedef struct {
    int d, kind;                 /* kind 0 = end (N = 1), 1 = mobile (N = 2 first, then N = 1) */
    int occ[LQ_DMAX][4];
    int N[LQ_DMAX];
    int lab[LQ_DMAX][3];         /* (2Sz, 2Iz or n_hole, n_hole) */
} LqLoc;

typedef struct { int dA, dB; int *cp, *ri; double *v; } LqTwo;   /* column-sparse, pair index a*dB + b */

typedef struct {                 /* a charge sector: basis + column-sparse H */
    long D; int Nh; int tgt[3];
    int64_t *basis; long *cp; int *ri; double *hv;
} LqSec;

static int      lq_M, lq_n, lq_nlab, lq_nm, lq_muof[4];
static double   lq_P[4][4];
static int      lq_mlab[4][2];
static double   lq_Ef1, lq_Ef2, lq_hb, lq_Om, lq_Dl, lq_g, lq_J0, lq_t;
/* v102/bind: the (1,1) link's pricing (lockq_p11: 0 none, 1 core, 2 singlet
 * core, 3 field branch, 4 space branch), as lo.py P11; lq_p11c its coefficient */
static double   lq_p11c;
/* v102 K2e: lockq_proj projects each jump onto the allowed subspace of the
 * (1,1) constraint (p11 = 2: no conversion singlet on any (1,1) link; p11 = 1:
 * no (1,1) link at all), by alternating projections; lq_projmax = most sweeps */
static long     lq_projmax;
static LqLoc    lq_loc[2];
static LqTwo    lq_two[LQ_NMAX];
static int     *lq_cell;         /* chain site -> kernel cell */
static int     *lq_site;         /* kernel cell -> chain site, or -1 */
static double  *lq_rawA;         /* lens area of fenced slots, before the fence */
static int      lq_dim[LQ_NMAX];
static int64_t  lq_stride[LQ_NMAX];
static LqSec    lq_cur, lq_nxt[4];   /* the state's sector, and the sector one jump (photon mode nu) away */
static double complex *lq_psi, *lq_V, *lq_w, *lq_phi[4], *lq_phw;
static double   lq_Egs, lq_time, lq_kerr;
static long     lq_events, lq_winviol, lq_steps, lq_maxper, lq_chkviol, lq_evmode[4];
static double   lq_Rmax, lq_Ein, lq_Eout, lq_dxmax;
static uint64_t lq_rng;

static double lq_rand(void)
{
    uint64_t x = lq_rng; x ^= x << 13; x ^= x >> 7; x ^= x << 17; lq_rng = x;
    return (double)(x >> 11) * (1.0 / 9007199254740992.0);
}

/* ---- site spaces (cell_states order: descending occupation tuples) ---- */
static void lq_gen(LqLoc *L, int Ntot, int *o, int k, int left)
{
    if (k == lq_M - 1) {
        o[k] = left;
        int d = L->d++;
        int s = 0, t = 0;
        for (int a = 0; a < lq_M; a++) {
            L->occ[d][a] = o[a];
            if (lq_M == 4) { s += o[a] * lq_mlab[a][0]; t += o[a] * lq_mlab[a][1]; }
        }
        if (lq_M == 2) s = o[0] - o[1];
        int nh = (L->kind == 1 && Ntot == 1);
        L->N[d] = Ntot;
        L->lab[d][0] = s;
        L->lab[d][1] = lq_M == 4 ? t : nh;
        L->lab[d][2] = nh;
        return;
    }
    for (int x = left; x >= 0; x--) { o[k] = x; lq_gen(L, Ntot, o, k + 1, left - x); }
}

static void lq_loc_build(LqLoc *L, int kind)
{
    int o[4] = {0, 0, 0, 0};
    L->d = 0; L->kind = kind;
    if (kind == 1) lq_gen(L, 2, o, 0, 2);
    lq_gen(L, 1, o, 0, 1);
}

static int lq_find(const LqLoc *L, const int *o)
{
    for (int k = 0; k < L->d; k++) {
        int same = 1;
        for (int a = 0; a < lq_M; a++) if (L->occ[k][a] != o[a]) { same = 0; break; }
        if (same) return k;
    }
    return -1;
}

/* E_rm = b^+_r b_m on state k: new index (or -1) and amplitude */
static int lq_Eop(const LqLoc *L, int k, int r, int m, double *amp)
{
    if (L->occ[k][m] == 0) return -1;
    int o[4];
    for (int a = 0; a < lq_M; a++) o[a] = L->occ[k][a];
    o[m]--; o[r]++;
    int k2 = lq_find(L, o);
    if (k2 < 0) return -1;
    *amp = sqrt((double)L->occ[k][m] * (r != m ? o[r] : L->occ[k][m]));
    return k2;
}

/* b_m (dir = -1) or b^+_m (dir = +1) on state k, within the site space */
static int lq_bop(const LqLoc *L, int k, int m, int dir, double *amp)
{
    int o[4];
    for (int a = 0; a < lq_M; a++) o[a] = L->occ[k][a];
    if (dir < 0) { if (o[m] == 0) return -1; *amp = sqrt((double)o[m]); o[m]--; }
    else { *amp = sqrt((double)o[m] + 1); o[m]++; }
    return lq_find(L, o);
}

static double lq_dtype(int ni, int nj)
{
    if (ni == 2 && nj == 2) return lq_Dl;
    return lq_Om - (lq_Ef2 - lq_Ef1) - lq_Ef1;    /* (2,1) and (1,2) */
}

/* the two-site link matrix, as lo.two_site_matrix(link_terms + hop_terms) */
static void lq_two_build(LqTwo *T, const LqLoc *A, const LqLoc *B)
{
    int dA = A->d, dB = B->d, dd = dA * dB;
    double *X = calloc((size_t)dd * dd, sizeof(double));
    for (int a = 0; a < dA; a++)
        for (int b = 0; b < dB; b++) {
            int ni = A->N[a], nj = B->N[b], col = a * dB + b;
            int n11 = (ni == 1 && nj == 1);
            if (n11 && P.lockq_p11 == 1) X[(size_t)col * dd + col] += lq_p11c;
            if (!n11 || P.lockq_p11 >= 2) {
                double c = n11 ? lq_p11c : -lq_g * lq_g / lq_dtype(ni, nj);
                for (int r = 0; r < lq_M; r++) for (int m = 0; m < lq_M; m++)
                    for (int s = 0; s < lq_M; s++) for (int n = 0; n < lq_M; n++) {
                        double w = lq_P[r][s] * lq_P[m][n];
                        if (w == 0) continue;
                        double ea, eb;
                        int a2 = lq_Eop(A, a, r, m, &ea);
                        if (a2 < 0) continue;
                        int b2 = lq_Eop(B, b, s, n, &eb);
                        if (b2 < 0) continue;
                        X[(size_t)(a2 * dB + b2) * dd + col] += c * w * ea * eb;
                    }
            }
            if (lq_t != 0 && A->kind == 1 && B->kind == 1)
                for (int m = 0; m < lq_M; m++) {
                    double ea, eb;
                    int a2 = lq_bop(A, a, m, -1, &ea), b2 = a2 >= 0 ? lq_bop(B, b, m, +1, &eb) : -1;
                    if (a2 >= 0 && b2 >= 0) X[(size_t)(a2 * dB + b2) * dd + col] += lq_t * ea * eb;
                    a2 = lq_bop(A, a, m, +1, &ea); b2 = a2 >= 0 ? lq_bop(B, b, m, -1, &eb) : -1;
                    if (a2 >= 0 && b2 >= 0) X[(size_t)(a2 * dB + b2) * dd + col] += lq_t * ea * eb;
                }
        }
    T->dA = dA; T->dB = dB;
    T->cp = calloc(dd + 1, sizeof(int));
    int nnz = 0;
    for (int col = 0; col < dd; col++)
        for (int row = 0; row < dd; row++) if (fabs(X[(size_t)row * dd + col]) > 1e-14) nnz++;
    T->ri = malloc((nnz + 1) * sizeof(int)); T->v = malloc((nnz + 1) * sizeof(double));
    nnz = 0;
    for (int col = 0; col < dd; col++) {
        T->cp[col] = nnz;
        for (int row = 0; row < dd; row++) {
            double x = X[(size_t)row * dd + col];
            if (fabs(x) > 1e-14) { T->ri[nnz] = row; T->v[nnz] = x; nnz++; }
        }
    }
    T->cp[dd] = nnz;
    free(X);
}

static const LqLoc *lq_locof(int site) { return &lq_loc[(site == 0 || site == lq_n - 1) ? 0 : 1]; }
static int lq_digit(int64_t code, int site) { return (int)((code / lq_stride[site]) % lq_dim[site]); }

static long lq_bsearch(const LqSec *S, int64_t c)
{
    long lo = 0, hi = S->D - 1;
    while (lo <= hi) {
        long mid = (lo + hi) >> 1;
        if (S->basis[mid] == c) return mid;
        if (S->basis[mid] < c) lo = mid + 1; else hi = mid - 1;
    }
    return -1;
}

static void lq_sec_free(LqSec *S)
{
    free(S->basis); free(S->cp); free(S->ri); free(S->hv);
    S->basis = NULL; S->cp = NULL; S->ri = NULL; S->hv = NULL; S->D = 0;
}

/* the sector (2Sz[, 2Iz], n_hole) = tgt: basis + column-sparse H.
 * An empty sector is left with D = 0. */
static void lq_sec_build(LqSec *S, const int *tgt, int Nh)
{
    lq_sec_free(S);
    S->Nh = Nh;
    for (int q = 0; q < 3; q++) S->tgt[q] = tgt[q];
    int64_t all = 1;
    for (int s = lq_n - 1; s >= 0; s--) { lq_stride[s] = all; lq_dim[s] = lq_locof(s)->d; all *= lq_dim[s]; }
    long cnt = 0;
    for (int pass = 0; pass < 2; pass++) {
        if (pass == 1) {
            if (cnt == 0) return;
            S->basis = malloc(cnt * sizeof(int64_t));
        }
        cnt = 0;
        for (int64_t c = 0; c < all; c++) {
            int tot[3] = {0, 0, 0};
            int64_t x = c;
            for (int s = lq_n - 1; s >= 0; s--) {
                int dg = (int)(x % lq_dim[s]); x /= lq_dim[s];
                const LqLoc *L = lq_locof(s);
                for (int q = 0; q < lq_nlab; q++) tot[q] += L->lab[dg][q];
            }
            int ok = 1;
            for (int q = 0; q < lq_nlab; q++) if (tot[q] != tgt[q]) { ok = 0; break; }
            if (!ok) continue;
            if (pass == 1) S->basis[cnt] = c;
            cnt++;
        }
    }
    S->D = cnt;
    long nnz = 0;
    for (long k = 0; k < S->D; k++)
        for (int l = 0; l + 1 < lq_n; l++) {
            int p = lq_digit(S->basis[k], l) * lq_two[l].dB + lq_digit(S->basis[k], l + 1);
            nnz += lq_two[l].cp[p + 1] - lq_two[l].cp[p];
        }
    S->cp = malloc((S->D + 1) * sizeof(long));
    S->ri = malloc((nnz + 1) * sizeof(int));
    S->hv = malloc((nnz + 1) * sizeof(double));
    nnz = 0;
    for (long k = 0; k < S->D; k++) {
        S->cp[k] = nnz;
        int64_t c = S->basis[k];
        for (int l = 0; l + 1 < lq_n; l++) {
            const LqTwo *T = &lq_two[l];
            int da = lq_digit(c, l), db = lq_digit(c, l + 1), p = da * T->dB + db;
            for (int e = T->cp[p]; e < T->cp[p + 1]; e++) {
                int q = T->ri[e], a2 = q / T->dB, b2 = q % T->dB;
                int64_t c2 = c + (int64_t)(a2 - da) * lq_stride[l] + (int64_t)(b2 - db) * lq_stride[l + 1];
                long row = lq_bsearch(S, c2);
                if (row < 0) { fprintf(stderr, "FATAL lockq: a term leaves the sector\n"); exit(3); }
                S->ri[nnz] = (int)row; S->hv[nnz] = T->v[e]; nnz++;
            }
        }
    }
    S->cp[S->D] = nnz;
}

static void lq_matvec(const LqSec *S, const double complex *x, double complex *y)
{
    for (long k = 0; k < S->D; k++) y[k] = 0;
    for (long k = 0; k < S->D; k++) {
        double complex xk = x[k];
        if (xk == 0) continue;
        for (long e = S->cp[k]; e < S->cp[k + 1]; e++) y[S->ri[e]] += S->hv[e] * xk;
    }
}

static double complex lq_dot(long D, const double complex *a, const double complex *b)
{
    double complex s = 0;
    for (long k = 0; k < D; k++) s += conj(a[k]) * b[k];
    return s;
}

static double lq_nrm(long D, const double complex *a) { return sqrt(creal(lq_dot(D, a, a))); }

/* <x|H|x>/<x|x> in sector S (y: scratch of length S->D) */
static double lq_expect(const LqSec *S, const double complex *x, double complex *y)
{
    lq_matvec(S, x, y);
    return creal(lq_dot(S->D, x, y)) / creal(lq_dot(S->D, x, x));
}

/* implicit QL on a symmetric tridiagonal (d diagonal, off[i] couples i,i+1);
 * z (n x n, row-major) starts as the identity and ends with eigenvectors in columns */
static void lq_tqli(double *d, const double *off, int n, double *z)
{
    double e[LQ_KMAX];
    for (int i = 0; i < n - 1; i++) e[i] = off[i];
    e[n - 1] = 0.0;
    for (int l = 0; l < n; l++) {
        int iter = 0, m;
        do {
            for (m = l; m < n - 1; m++) {
                double dd = fabs(d[m]) + fabs(d[m + 1]);
                if (fabs(e[m]) <= 2.2e-16 * dd) break;
            }
            if (m != l) {
                if (iter++ == 60) break;
                double g = (d[l + 1] - d[l]) / (2.0 * e[l]);
                double r = hypot(g, 1.0);
                g = d[m] - d[l] + e[l] / (g + (g >= 0 ? r : -r));
                double s = 1.0, c = 1.0, p = 0.0;
                int i;
                for (i = m - 1; i >= l; i--) {
                    double f = s * e[i], b = c * e[i];
                    e[i + 1] = (r = hypot(f, g));
                    if (r == 0.0) { d[i + 1] -= p; e[m] = 0.0; break; }
                    s = f / r; c = g / r;
                    g = d[i + 1] - p;
                    r = (d[i] - g) * s + 2.0 * c * b;
                    d[i + 1] = g + (p = s * r);
                    g = c * r - b;
                    for (int k = 0; k < n; k++) {
                        f = z[k * n + i + 1];
                        z[k * n + i + 1] = s * z[k * n + i] + c * f;
                        z[k * n + i] = c * z[k * n + i] - s * f;
                    }
                }
                if (r == 0.0 && i >= l) continue;
                d[l] -= p; e[l] = g; e[m] = 0.0;
            }
        } while (m != l);
    }
}

/* Lanczos in the current sector with full (twice-repeated)
 * reorthogonalization from v0; returns the number of basis vectors built
 * (stops early on an invariant subspace) */
static int lq_lanczos(const double complex *v0, int m, double *a, double *b)
{
    long D = lq_cur.D;
    double nr = lq_nrm(D, v0);
    for (long k = 0; k < D; k++) lq_V[k] = v0[k] / nr;
    if (m > D) m = (int)D;
    for (int j = 0; j < m; j++) {
        double complex *vj = lq_V + (size_t)j * D;
        lq_matvec(&lq_cur, vj, lq_w);
        a[j] = creal(lq_dot(D, vj, lq_w));
        for (int pass = 0; pass < 2; pass++)
            for (int i = 0; i <= j; i++) {
                double complex *vi = lq_V + (size_t)i * D;
                double complex c = lq_dot(D, vi, lq_w);
                for (long k = 0; k < D; k++) lq_w[k] -= c * vi[k];
            }
        b[j] = lq_nrm(D, lq_w);
        if (j + 1 == m || b[j] < 1e-13 * (fabs(a[j]) + 1e-30)) return j + 1;
        double complex *vn = lq_V + (size_t)(j + 1) * D;
        for (long k = 0; k < D; k++) vn[k] = lq_w[k] / b[j];
    }
    return m;
}

/* scratch sized to the current and next sectors */
static void lq_scratch(void)
{
    free(lq_V); free(lq_w); free(lq_phw);
    lq_V = malloc((size_t)LQ_KMAX * lq_cur.D * sizeof(double complex));
    lq_w = malloc((size_t)lq_cur.D * sizeof(double complex));
    long Dmax = 1;
    for (int nu = 0; nu < lq_nm; nu++) {
        free(lq_phi[nu]);
        long Dn = lq_nxt[nu].D > 0 ? lq_nxt[nu].D : 1;
        lq_phi[nu] = malloc((size_t)Dn * sizeof(double complex));
        if (Dn > Dmax) Dmax = Dn;
    }
    lq_phw = malloc((size_t)Dmax * sizeof(double complex));
}

/* the sectors one jump (photon mode nu) away from the current one */
static void lq_next_build(void)
{
    for (int nu = 0; nu < lq_nm; nu++) {
        if (lq_cur.Nh >= lq_n - 2) { lq_sec_free(&lq_nxt[nu]); continue; }
        int tgt[3];
        for (int q = 0; q < 3; q++) tgt[q] = lq_cur.tgt[q];
        for (int q = 0; q < lq_nlab - 1; q++) tgt[q] -= lq_mlab[lq_muof[nu]][q];
        tgt[lq_nlab - 1] += 1;
        lq_sec_build(&lq_nxt[nu], tgt, lq_cur.Nh + 1);
    }
}

/* the current sector's ground state into lq_psi (restarted Lanczos) */
static double lq_ground(void)
{
    double a[LQ_KMAX], b[LQ_KMAX], d[LQ_KMAX], z[LQ_KMAX * LQ_KMAX];
    long D = lq_cur.D;
    uint64_t r = 0x9E3779B97F4A7C15ULL;
    for (long k = 0; k < D; k++) {
        r ^= r << 13; r ^= r >> 7; r ^= r << 17;
        lq_psi[k] = (double)(r >> 11) / 9007199254740992.0 - 0.5;
    }
    double E = 0;
    for (int it = 0; it < 200; it++) {
        int m = lq_lanczos(lq_psi, LQ_KMAX, a, b);
        for (int i = 0; i < m; i++) { d[i] = a[i]; for (int j = 0; j < m; j++) z[i * m + j] = i == j; }
        lq_tqli(d, b, m, z);
        int lo = 0;
        for (int i = 1; i < m; i++) if (d[i] < d[lo]) lo = i;
        for (long k = 0; k < D; k++) lq_psi[k] = 0;
        for (int j = 0; j < m; j++) {
            double c = z[j * m + lo];
            const double complex *vj = lq_V + (size_t)j * D;
            for (long k = 0; k < D; k++) lq_psi[k] += c * vj[k];
        }
        double nr = lq_nrm(D, lq_psi);
        for (long k = 0; k < D; k++) lq_psi[k] /= nr;
        lq_matvec(&lq_cur, lq_psi, lq_w);
        E = creal(lq_dot(D, lq_psi, lq_w));
        double res = 0;
        for (long k = 0; k < D; k++) { double complex q = lq_w[k] - E * lq_psi[k]; res += creal(q * conj(q)); }
        if (sqrt(res) < 1e-12 || m == D) break;
    }
    return E;
}

/* psi <- exp(-i H tau) psi (Krylov; split tau until the a-posteriori bound holds) */
static void lq_evolve(double tau, int depth)
{
    double a[LQ_KMAX], b[LQ_KMAX], d[LQ_KMAX], z[LQ_KMAX * LQ_KMAX];
    double complex cf[LQ_KMAX];
    long D = lq_cur.D;
    double nr = lq_nrm(D, lq_psi);
    int m = lq_lanczos(lq_psi, LQ_KMAX < 24 ? LQ_KMAX : 24, a, b);
    for (int i = 0; i < m; i++) { d[i] = a[i]; for (int j = 0; j < m; j++) z[i * m + j] = i == j; }
    lq_tqli(d, b, m, z);
    for (int j = 0; j < m; j++) {
        double complex s = 0;
        for (int i = 0; i < m; i++) s += z[j * m + i] * cexp(-I * d[i] * tau) * z[0 * m + i];
        cf[j] = s;
    }
    int full = (m == D) || (b[m - 1] < 1e-13 * (fabs(a[m - 1]) + 1e-30));
    double err = full ? 0 : b[m - 1] * cabs(cf[m - 1]);
    if (err > 1e-14 && depth < 12) { lq_evolve(0.5 * tau, depth + 1); lq_evolve(0.5 * tau, depth + 1); return; }
    if (err > lq_kerr) lq_kerr = err;
    for (long k = 0; k < D; k++) lq_psi[k] = 0;
    for (int j = 0; j < m; j++) {
        const double complex *vj = lq_V + (size_t)j * D;
        for (long k = 0; k < D; k++) lq_psi[k] += nr * cf[j] * vj[k];
    }
}

/* sum of E_f(N_i): fixed by the hole count */
static double lq_efsum(int Nh) { return 2 * lq_Ef1 + (lq_n - 2 - Nh) * lq_Ef2 + Nh * lq_Ef1; }

/* the lock's energy in kernel units: hbar * (sum E_f + <H>) */
static double lq_energy_kernel(void)
{
    if (!P.lockq || !lq_psi) return 0;
    return lq_hb * (lq_efsum(lq_cur.Nh) + lq_expect(&lq_cur, lq_psi, lq_w));
}

/* lq_phi[nu] <- J_nu,i psi in the next sector (unnormalized); returns |J psi|^2 */
/* apply (1 - Pi_l) to phi on link l; returns the squared norm removed */
static double lq_proj_link(const LqSec *X, double complex *phi, int l)
{
    const LqLoc *A = lq_locof(l), *B = lq_locof(l + 1);
    double c2 = 0, rem = 0;
    int cmp[4], d1A[4], d1B[4];
    for (int r = 0; r < lq_M; r++) for (int q = 0; q < lq_M; q++) c2 += lq_P[r][q] * lq_P[r][q];
    for (int r = 0; r < lq_M; r++) {
        int o[4] = {0, 0, 0, 0};
        cmp[r] = -1;
        for (int q = 0; q < lq_M; q++) if (lq_P[r][q] != 0) cmp[r] = q;
        o[r] = 1; d1A[r] = lq_find(A, o); d1B[r] = lq_find(B, o);
    }
    double cn = 1.0 / sqrt(c2);
    for (long k = 0; k < X->D; k++) {
        int64_t c = X->basis[k];
        int da = lq_digit(c, l), db = lq_digit(c, l + 1);
        if (A->N[da] != 1 || B->N[db] != 1) continue;
        if (P.lockq_p11 == 1) { rem += creal(phi[k] * conj(phi[k])); phi[k] = 0; continue; }
        if (da != d1A[0] || db != d1B[cmp[0]]) continue;          /* the group's representative: (0, cmp(0)) */
        long row[4];
        double complex a = 0;
        for (int r = 0; r < lq_M; r++) {
            row[r] = lq_bsearch(X, c + (int64_t)(d1A[r] - da) * lq_stride[l] + (int64_t)(d1B[cmp[r]] - db) * lq_stride[l + 1]);
            if (row[r] < 0) { fprintf(stderr, "FATAL lockq: singlet partner outside the sector\n"); exit(3); }
            a += lq_P[r][cmp[r]] * cn * phi[row[r]];
        }
        for (int r = 0; r < lq_M; r++) phi[row[r]] -= lq_P[r][cmp[r]] * cn * a;
        rem += creal(a * conj(a));
    }
    return rem;
}

static void lq_proj(const LqSec *X, double complex *phi)
{
    double n2 = creal(lq_dot(X->D, phi, phi));
    long sw;
    for (sw = 1; sw <= 200; sw++) {
        double rem = 0;
        for (int l = 0; l + 1 < lq_n; l++) rem += lq_proj_link(X, phi, l);
        if (rem <= 1e-30 * n2) break;
    }
    if (sw > lq_projmax) lq_projmax = sw;
}

static double lq_map(int i, int nu)
{
    const LqLoc *L = lq_locof(i);
    const LqSec *X = &lq_nxt[nu];
    double complex *phi = lq_phi[nu];
    for (long k = 0; k < X->D; k++) phi[k] = 0;
    for (long k = 0; k < lq_cur.D; k++) {
        if (lq_psi[k] == 0) continue;
        int dg = lq_digit(lq_cur.basis[k], i);
        if (L->N[dg] != 2) continue;
        for (int mu = 0; mu < lq_M; mu++) {
            if (lq_P[mu][nu] == 0) continue;
            double amp;
            int d2 = lq_bop(L, dg, mu, -1, &amp);
            if (d2 < 0) continue;
            long row = lq_bsearch(X, lq_cur.basis[k] + (int64_t)(d2 - dg) * lq_stride[i]);
            if (row < 0) { fprintf(stderr, "FATAL lockq: jump leaves the target sector\n"); exit(3); }
            phi[row] += lq_P[mu][nu] * amp * lq_psi[k];
        }
    }
    if (P.lockq_proj) lq_proj(X, phi);
    return creal(lq_dot(X->D, phi, phi)) / creal(lq_dot(lq_cur.D, lq_psi, lq_psi));
}

/* the energy of photon mode nu in kernel cell o (F = 1: the whole field) */
static double lq_mode_e(int nu, int o) { return P.fmodes > 1 ? fx_mode_e(nu, o) : Ee[o]; }

/* the jump at bulk site i: the state moves to the next sector.  o is the
 * kernel neighbour that gave the atom (or -1 for a seeded jump) */
static void lq_fire(int i, int o, int nu)
{
    double Eb = lq_energy_kernel();
    double jj = lq_map(i, nu);
    if (jj <= 0) { fprintf(stderr, "FATAL lockq: zero jump\n"); exit(3); }
    if (o < 0 && P.lockq_proj) printf("# LOCKQ seeded jump weight (projected) %.15f\n", jj);
    int otgt[3] = {lq_cur.tgt[0], lq_cur.tgt[1], lq_cur.tgt[2]}, oNh = lq_cur.Nh;
    LqSec old = lq_cur;
    lq_cur = lq_nxt[nu];
    lq_nxt[nu].basis = NULL; lq_nxt[nu].cp = NULL; lq_nxt[nu].ri = NULL; lq_nxt[nu].hv = NULL; lq_nxt[nu].D = 0;
    lq_sec_free(&old);
    free(lq_psi);
    lq_psi = lq_phi[nu]; lq_phi[nu] = NULL;
    double nr = lq_nrm(lq_cur.D, lq_psi);
    for (long k = 0; k < lq_cur.D; k++) lq_psi[k] /= nr;
    lq_next_build();
    lq_scratch();
    double Ea = lq_energy_kernel();
    /* in-kernel bookkeeping check: the lock moves by mode nu's charges */
    int bad = lq_cur.Nh != oNh + 1;
    for (int q = 0; q < lq_nlab - 1; q++) if (lq_cur.tgt[q] - otgt[q] != -lq_mlab[lq_muof[nu]][q]) bad = 1;
    if (o >= 0) {
        double eps = A0eff * w1e[o] / TWO_PI;
        double em0 = lq_mode_e(nu, o);
        if (P.fmodes > 1) {
            double fac = sqrt((em0 - eps) / em0);
            if (nu == 0) { fa1[o] *= fac; fa2[o] *= fac; }
            else { fxr_(nu)[o] *= fac; fxi_(nu)[o] *= fac; }
            Ee[o] -= eps;
            if (fabs(em0 - lq_mode_e(nu, o) - eps) > 1e-12 * eps) bad = 1;
        } else {
            double ee = Ee[o], fac = sqrt((ee - eps) / ee);
            fa1[o] *= fac; fa2[o] *= fac;
            Ee[o] = ee - eps;
        }
        Es[o] += eps - (Ea - Eb);
        double dx = w1e[o] - (Ea - Eb) / lq_hb - lq_Om;
        if (fabs(dx) >= 2 * P.s_k) lq_winviol++;
        if (fabs(dx) > lq_dxmax) lq_dxmax = fabs(dx);
        lq_Ein += eps; lq_Eout += eps - (Ea - Eb);
        lq_events++; lq_evmode[nu]++;
        if (P.fmodes > 1)
            printf("# LQ event t=%.4f site=%d cell=%d o=%d nu=%d w1e_o=%.6f eps=%.9f dE_lock=%.9f space_out=%.9f delta_exact=%+.6f Nh=%d charges=(%d,%d)\n",
                   sim_t, i, lq_cell[i], o, nu, w1e[o], eps, Ea - Eb, eps - (Ea - Eb), dx, lq_cur.Nh,
                   lq_cur.tgt[0], lq_nlab == 3 ? lq_cur.tgt[1] : 0);
        else
            printf("# LQ event t=%.4f site=%d cell=%d o=%d w1e_o=%.6f eps=%.9f dE_lock=%.9f space_out=%.9f delta_exact=%+.6f Nh=%d\n",
                   sim_t, i, lq_cell[i], o, w1e[o], eps, Ea - Eb, eps - (Ea - Eb), dx, lq_cur.Nh);
    }
    if (bad) lq_chkviol++;
}

/* fence: called wherever the kernel sets a slot's lens area */
static double lq_fence(int s, double A)
{
    if (!P.lockq || !lq_site) return A;
    if (lq_site[sli[s]] >= 0 || lq_site[slj[s]] >= 0) { lq_rawA[s] = A; return 0; }
    return A;
}

static int lq_skip(int i) { return P.lockq && lq_site && lq_site[i] >= 0; }

/* one lock step: exact evolution over dt_L, then formation events as a
 * Poisson process over the step (rates recomputed after each event) */
#define LQ_CMAX 4096
static void lq_step(void)
{
    if (!P.lockq) return;
    static long kst = 0;
    if (++kst % P.lockq_every) return;
    double dtL = P.lockq_every * P.dt;
    lq_evolve(dtL, 0);
    lq_time += dtL; lq_steps++;
    if (P.lockq_rate <= 0) return;
    static int ci[LQ_CMAX], co[LQ_CMAX], cn[LQ_CMAX];
    static double cr_[LQ_CMAX];
    double band = 2 * P.s_k, rem = dtL;
    long nper = 0;
    for (;;) {
        if (lq_cur.Nh >= lq_n - 2) break;
        double Ex = lq_expect(&lq_cur, lq_psi, lq_w);
        int nc = 0;
        double R = 0;
        for (int i = 1; i < lq_n - 1; i++) {
            int c = lq_cell[i];
            for (int nu = 0; nu < lq_nm; nu++) {
                if (lq_nxt[nu].D == 0) continue;
                int lit = 0;
                for (int q = cls_[c]; q < cls_[c + 1] && !lit; q++) {
                    int s = clidx[q];
                    if (sst[s] == S_FREE || lq_rawA[s] <= 0) continue;
                    int o = sli[s] == c ? slj[s] : sli[s];
                    if (lq_site[o] >= 0) continue;
                    double eps = A0eff * w1e[o] / TWO_PI;
                    if (eps > 0 && lq_mode_e(nu, o) >= eps) lit = 1;
                }
                if (!lit) continue;
                double jj = lq_map(i, nu);
                if (jj <= 0) continue;
                double dE = (lq_Ef1 - lq_Ef2) + lq_expect(&lq_nxt[nu], lq_phi[nu], lq_phw) - Ex;   /* rate units */
                for (int q = cls_[c]; q < cls_[c + 1]; q++) {
                    int s = clidx[q];
                    if (sst[s] == S_FREE || lq_rawA[s] <= 0) continue;
                    int o = sli[s] == c ? slj[s] : sli[s];
                    if (lq_site[o] >= 0) continue;
                    double eps = A0eff * w1e[o] / TWO_PI;
                    if (eps <= 0) continue;
                    double na = floor(lq_mode_e(nu, o) / eps);
                    if (na < 1) continue;
                    double del = w1e[o] - dE - lq_Om;
                    if (fabs(del) >= band) continue;
                    double G = lq_g * lq_g * sqrt(band * band - del * del) / (P.s_k * P.s_k);
                    double r = P.lockq_rate * G * jj * na;
                    if (nc < LQ_CMAX) { ci[nc] = i; co[nc] = o; cn[nc] = nu; cr_[nc] = r; nc++; R += r; }
                }
            }
        }
        if (nper == 0 && R * dtL > lq_Rmax) lq_Rmax = R * dtL;
        if (R <= 0) break;
        double tau = -log(1.0 - lq_rand()) / R;
        if (tau >= rem) break;
        rem -= tau;
        double u = lq_rand() * R, acc = 0;
        int k = 0;
        for (k = 0; k < nc - 1; k++) { acc += cr_[k]; if (u < acc) break; }
        lq_fire(ci[k], co[k], cn[k]);
        nper++;
    }
    if (nper > lq_maxper) lq_maxper = nper;
}

static void lq_profile(double *nh, double *sz, double *iz)
{
    for (int s = 0; s < lq_n; s++) nh[s] = sz[s] = iz[s] = 0;
    double nn = 0;
    for (long k = 0; k < lq_cur.D; k++) {
        double p = creal(lq_psi[k] * conj(lq_psi[k]));
        nn += p;
        for (int s = 0; s < lq_n; s++) {
            const LqLoc *L = lq_locof(s);
            int dg = lq_digit(lq_cur.basis[k], s);
            nh[s] += p * L->lab[dg][2];
            sz[s] += p * 0.5 * L->lab[dg][0];
            if (lq_M == 4) iz[s] += p * 0.5 * L->lab[dg][1];
        }
    }
    for (int s = 0; s < lq_n; s++) { nh[s] /= nn; sz[s] /= nn; iz[s] /= nn; }
}

static void lq_diag(void)
{
    if (!P.lockq) return;
    double nh[LQ_NMAX], sz[LQ_NMAX], iz[LQ_NMAX];
    lq_profile(nh, sz, iz);
    printf("# LQ t=%.4f tL=%.6f E_lock=%.12f Ex/J0=%.10f Nh=%d dim=%ld events=%ld norm-1=%.2e kerr=%.1e | nh:",
           sim_t, lq_time, lq_energy_kernel(), lq_expect(&lq_cur, lq_psi, lq_w) / lq_J0, lq_cur.Nh,
           lq_cur.D, lq_events, lq_nrm(lq_cur.D, lq_psi) - 1, lq_kerr);
    for (int s = 0; s < lq_n; s++) printf(" %.10f", nh[s]);
    printf(" | sz:");
    for (int s = 0; s < lq_n; s++) printf(" %+.10f", sz[s]);
    if (lq_M == 4) { printf(" | iz:"); for (int s = 0; s < lq_n; s++) printf(" %+.10f", iz[s]); }
    printf("\n");
}

static void lq_result(void)
{
    if (!P.lockq) return;
    printf("# RESULT lockq events=%ld Nh=%d window_viol=%ld max|delta_exact|=%.6f max_expected_per_step=%.4f max_events_per_step=%ld atoms_in=%.9f space_out=%.9f lock_gain=%.9f kerr_max=%.1e tL=%.6f\n",
           lq_events, lq_cur.Nh, lq_winviol, lq_dxmax, lq_Rmax, lq_maxper, lq_Ein, lq_Eout, lq_Ein - lq_Eout, lq_kerr, lq_time);
    if (P.fmodes > 1)
        printf("# RESULT lockq_modes events_by_mode=%ld,%ld,%ld,%ld check_violations=%ld final_charges=(%d,%d) Nh=%d\n",
               lq_evmode[0], lq_evmode[1], lq_evmode[2], lq_evmode[3], lq_chkviol,
               lq_cur.tgt[0], lq_nlab == 3 ? lq_cur.tgt[1] : 0, lq_cur.Nh);
    if (P.lockq_proj) printf("# RESULT lockq_proj max_sweeps=%ld\n", lq_projmax);
    lq_diag();
}

/* installation at t = 0, after the initial topology and before E0 */
static void lq_install(double cx, double cy, double cz)
{
    if (!P.lockq) return;
    if (P.bath_xb > 0 || P.door_g > 0 || P.bh_r > 0 || P.k_cant > 0 || P.amp_nat > 0 || P.pattern_only) {
        printf("# CONFIG ERROR: lockq is fenced only for the V3a surface (bath_xb, door_g, bh_r, k_cant, amp_nat, pattern_only must be off)\n");
        exit(1);
    }
    lq_M = P.lockq_M; lq_n = P.lockq_n;
    if ((lq_M != 2 && lq_M != 4) || lq_n < 3 || lq_n > (lq_M == 2 ? 12 : 7) || P.lockq_every < 1) {
        printf("# CONFIG ERROR: lockq_M must be 2 or 4, 3 <= lockq_n <= %d, lockq_every >= 1\n", lq_M == 2 ? 12 : 7);
        exit(1);
    }
    lq_rng = P.lockq_seed ? P.lockq_seed : 1;
    /* S-law constants from the law table */
    lq_hb = A0eff / TWO_PI;
    double aa = P.q_detune * lq_hb / P.cap;
    lq_Ef1 = (-1.0 + sqrt(1.0 + 4.0 * aa * 1.0 * P.w1)) / (2.0 * aa);
    lq_Ef2 = (-1.0 + sqrt(1.0 + 4.0 * aa * 2.0 * P.w1)) / (2.0 * aa);
    lq_Om = P.w2;
    lq_Dl = lq_Om - 2.0 * (lq_Ef2 - lq_Ef1);
    lq_g = P.lockq_g * lq_Dl;
    lq_J0 = 2.0 * lq_g * lq_g / lq_Dl;
    lq_t = P.lockq_t * lq_J0;
    double eps2[2][2] = {{0, 1}, {-1, 0}};
    static const int mw[4][2] = {{+1, +1}, {+1, -1}, {-1, +1}, {-1, -1}};
    for (int a = 0; a < 4; a++) for (int b = 0; b < 4; b++) lq_P[a][b] = 0;
    if (lq_M == 2) {
        for (int a = 0; a < 2; a++) for (int b = 0; b < 2; b++) lq_P[a][b] = eps2[a][b];
        lq_mlab[0][0] = +1; lq_mlab[1][0] = -1;
        lq_nlab = 2;
    } else {
        for (int a = 0; a < 4; a++) for (int b = 0; b < 4; b++)
            lq_P[a][b] = eps2[a / 2][b / 2] * eps2[a % 2][b % 2];
        for (int a = 0; a < 4; a++) { lq_mlab[a][0] = mw[a][0]; lq_mlab[a][1] = mw[a][1]; }
        lq_nlab = 3;
    }
    if (P.lockq_p11 < 0 || P.lockq_p11 > 4) { printf("# CONFIG ERROR: lockq_p11 must be 0 .. 4\n"); exit(1); }
    if (P.lockq_proj && P.lockq_p11 != 1 && P.lockq_p11 != 2) { printf("# CONFIG ERROR: lockq_proj needs lockq_p11 = 1 or 2\n"); exit(1); }
    if (P.lockq_proj) printf("# LOCKQ proj=1: each jump is projected onto the allowed subspace of the (1,1) constraint\n");
    if (P.lockq_p11) {
        double c2 = 0, d11 = lq_Om - 2.0 * lq_Ef1;
        for (int a = 0; a < lq_M; a++) for (int b = 0; b < lq_M; b++) c2 += lq_P[a][b] * lq_P[a][b];
        double root = sqrt(d11 * d11 / 4 + lq_g * lq_g * c2);
        lq_p11c = P.lockq_p11 == 1 ? 1.0 : P.lockq_p11 == 2 ? 1.0 / c2
                : P.lockq_p11 == 3 ? (d11 / 2 + root) / c2 : (d11 / 2 - root) / c2;
        printf("# LOCKQ p11=%d: the (1,1) link priced, coefficient %.15f (Delta11=%.15f |P|^2=%g)\n",
               P.lockq_p11, lq_p11c, d11, c2);
    }
    if (P.fmodes > 1 && P.fmodes != lq_M) {
        printf("# CONFIG ERROR: with fmodes > 1, lockq_M must equal fmodes\n");
        exit(1);
    }
    lq_nm = P.fmodes > 1 ? P.fmodes : 1;
    for (int nu = 0; nu < lq_nm; nu++) {
        lq_muof[nu] = -1;
        for (int mu = 0; mu < lq_M; mu++) if (lq_P[mu][nu] != 0) {
            if (lq_muof[nu] >= 0) { printf("# CONFIG ERROR: lockq photon mode couples to two lock modes\n"); exit(1); }
            lq_muof[nu] = mu;
        }
    }
    lq_loc_build(&lq_loc[0], 0);
    lq_loc_build(&lq_loc[1], 1);
    for (int l = 0; l + 1 < lq_n; l++) lq_two_build(&lq_two[l], lq_locof(l), lq_locof(l + 1));

    /* the chain: nearest cell to (cx + dx, cy, cz), then greedily the
     * lens-contact neighbour farthest along +x not yet used */
    lq_cell = malloc(lq_n * sizeof(int));
    lq_site = malloc(NC * sizeof(int));
    for (int i = 0; i < NC; i++) lq_site[i] = -1;
    double tx = cx + P.lockq_dx, best = 1e300;
    int c0 = -1;
    for (int i = 0; i < NC; i++) {
        if (pin[i]) continue;
        double dx = wr(px_[i] - tx), dy = wr(py_[i] - cy), dz = wr(pz_[i] - cz);
        double r2 = dx * dx + dy * dy + dz * dz;
        if (r2 < best) { best = r2; c0 = i; }
    }
    lq_cell[0] = c0; lq_site[c0] = 0;
    for (int k = 1; k < lq_n; k++) {
        int cur = lq_cell[k - 1], nx = -1;
        double bx = -1e300;
        for (int q = cls_[cur]; q < cls_[cur + 1]; q++) {
            int s = clidx[q];
            if (sst[s] == S_FREE || sA[s] <= 0) continue;
            int o = sli[s] == cur ? slj[s] : sli[s];
            if (lq_site[o] >= 0 || pin[o]) continue;
            double ax = wr(px_[o] - px_[cur]);
            if (ax > bx) { bx = ax; nx = o; }
        }
        if (nx < 0) { printf("# CONFIG ERROR: lockq chain stalls at site %d\n", k); exit(1); }
        lq_cell[k] = nx; lq_site[nx] = k;
    }
    /* seed: remove the lock cells' classical contents; in-flight dense
     * energy on fenced slots goes to the non-lock end's Em */
    double removed = 0, moved = 0;
    for (int k = 0; k < lq_n; k++) {
        int c = lq_cell[k];
        removed += Ee[c] + Em[c];
        Ee[c] = Em[c] = 0; fa1[c] = fa2[c] = 0;
        if (P.fmodes > 1) for (int nu = 1; nu < P.fmodes; nu++) { fxr_(nu)[c] = 0; fxi_(nu)[c] = 0; }   /* K2b: every mode */
        pin[c] = 1;
    }
    lq_rawA = calloc(NLMAX, sizeof(double));   /* the slot table's capacity: NSLOT is only its high-water mark */
    int nfence = 0;
    for (int s = 0; s < NSLOT; s++) {
        if (sst[s] == S_FREE) continue;
        int li = lq_site[sli[s]] >= 0, lj = lq_site[slj[s]] >= 0;
        if (!li && !lj) continue;
        double fl = slem[2 * s] + slem[2 * s + 1];
        if (fl != 0) {
            if (!li) Em[sli[s]] += fl; else if (!lj) Em[slj[s]] += fl; else removed += fl;
            moved += fl;
            slem[2 * s] = slem[2 * s + 1] = 0;
        }
        sA[s] = lq_fence(s, sA[s]);
        nfence++;
    }
    /* the lock state: the no-hole, zero-charge sector's ground state */
    int z0[3] = {0, 0, 0};
    lq_sec_build(&lq_cur, z0, 0);
    lq_next_build();
    lq_scratch();
    lq_psi = calloc(lq_cur.D, sizeof(double complex));
    lq_Egs = lq_ground();
    printf("# LOCKQ v102 K2a: M=%d n=%d g/Delta=%g t/J0=%g rate=%g every=%d seed=%lu dx=%g\n",
           lq_M, lq_n, P.lockq_g, P.lockq_t, P.lockq_rate, P.lockq_every, P.lockq_seed, P.lockq_dx);
    printf("# LOCKQ S-law: Ef1=%.15f Ef2=%.15f Omega=%g Delta=%.15f g=%.15f J0=%.15f t=%.6e hbar=%.15f\n",
           lq_Ef1, lq_Ef2, lq_Om, lq_Dl, lq_g, lq_J0, lq_t, lq_hb);
    printf("# LOCKQ chain cells:");
    for (int k = 0; k < lq_n; k++) printf(" %d", lq_cell[k]);
    printf(" | fenced slots %d | removed classical %.9f (in-flight moved %.9f)\n", nfence, removed, moved);
    printf("# LOCKQ ground: sector dim %ld E_gs=%.15f (rate) = %.12f J0 ; lock energy %.12f (kernel)\n",
           lq_cur.D, lq_Egs, lq_Egs / lq_J0, lq_energy_kernel());
    printf("# LOCKQ resonance: photon pitch %.6f +- %.6f on the E_f part (vacuum w1=%g); the exact window is per jump\n",
           lq_Om - (lq_Ef2 - lq_Ef1), 2 * P.s_k, P.w1);
    for (int sj = 0; sj < 2; sj++) {
    int fsite = sj == 0 ? P.lockq_force : P.lockq_force2;
    if (fsite >= 1 && fsite <= lq_n - 2) {
        int fm = sj == 0 ? P.lockq_fmode : P.lockq_fmode2;
        if (fm < 0 || fm >= lq_M) { printf("# CONFIG ERROR: lockq_fmode must be 0 .. lockq_M-1\n"); exit(1); }
        if (fm >= lq_nm) {   /* a seeded jump may use any lock mode, even with a one-mode field */
            int nm0 = lq_nm;
            lq_nm = fm + 1;
            for (int nu = nm0; nu < lq_nm; nu++)
                for (int mu = 0; mu < lq_M; mu++) if (lq_P[mu][nu] != 0) lq_muof[nu] = mu;
            lq_next_build(); lq_scratch();
            lq_fire(fsite, -1, fm);
            for (int nu = nm0; nu < lq_nm; nu++) { lq_sec_free(&lq_nxt[nu]); free(lq_phi[nu]); lq_phi[nu] = NULL; }
            lq_nm = nm0;
            lq_scratch();
        } else lq_fire(fsite, -1, fm);
        printf("# LOCKQ seeded jump at site %d, photon mode %d: sector dim %ld, E_lock %.12f (kernel), <H>/J0 %.12f\n",
               fsite, fm, lq_cur.D, lq_energy_kernel(), lq_expect(&lq_cur, lq_psi, lq_w) / lq_J0);
    }
    }
    lq_diag();
}
