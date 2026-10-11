/* gx.c — v102 gather step (../PREREG.md): the leading-order fiber on a 1-D
 * chain with golden-rule events, under the eight rule sets of candidates
 * a (reverse boundary event), b (conversion in an unlabelled fiber) and
 * c (membership by saturation).
 *
 * Chain: n cells, cell 0 at a wall (z = 1), cell n-1 on a light lead.
 * Cells hold bosons in M modes (M = 2: 3+2, P = eps; M = 4: W, P = eps x eps).
 * H = sum_i E_f(N_i) + exchange (LO: -(g^2/Delta_type) A^+A, R12 core on
 * (1,1)) + hops (t_in / t_L by membership).  Every sector (N, 2Sz, 2Iz) is
 * diagonalized in full (LAPACK dsyevd).  Events (light out/in at the lead,
 * space out/in on links, side formation / side a) are golden-rule jumps
 * between eigenvalue clusters; kinetic Monte Carlo carries the state
 * vector coherently within its cluster.
 *
 * Build: gcc -O3 -march=native -fopenmp -Wall -Wextra -o gx gx.c \
 *        -l:liblapack.so.3 -l:libopenblas.so.0 -lm
 * Usage: gx key=value ...  (see parse())                                   */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <time.h>
#include <complex.h>

extern void dsyevd_(const char *jobz, const char *uplo, const int *n, double *a, const int *lda, double *w,
                    double *work, const int *lwork, int *iwork, const int *liwork, int *info);
extern void dgemv_(const char *trans, const int *m, const int *n, const double *alpha, const double *a, const int *lda,
                   const double *x, const int *incx, const double *beta, double *y, const int *incy);
extern void dgemm_(const char *ta, const char *tb, const int *m, const int *n, const int *k, const double *alpha,
                   const double *a, const int *lda, const double *b, const int *ldb, const double *beta, double *c, const int *ldc);

/* ---------------- constants ---------------- */
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
static int M = 2, n = 6, kL = 4, ra = 0, rb = 0, rc = 0, side = 0;
static double gr = 0.1, tin_J0 = 0.1, tL = 0.9, sk = 0.06, tin;
static char exp_[16] = "static";
static double nu = 0, Tmax = 5000, Rmin = 1e-12, Tdark = 1e9;
static int ntraj = 200, Amax = 100, ntraj_dense = 0;
static double Tdense = 0;
static uint64_t seed0 = 1;
static double ns_env = 0, nL_env = 0;

static double P[4][4];
static int ms2[4], mt2[4];
static double P2;

/* ---------------- cells ---------------- */
#define MAXL 16
typedef struct { int d, cap; int occ[MAXL][4]; int N[MAXL], s2[MAXL], t2[MAXL]; int key2idx[81]; } Loc;
static Loc loc[3];
static Loc *site[16];
static int cap[16], z[16], inL[16];
static int radix[16];


static int okey(const int *o) { int k = 0, p = 1; for (int a = 0; a < M; a++) { k += o[a] * p; p *= 3; } return k; }

static void loc_add(Loc *L, const int *o)
{
    int d = L->d++, s = 0, t = 0, N = 0;
    for (int b = 0; b < M; b++) { L->occ[d][b] = o[b]; N += o[b]; s += o[b] * ms2[b]; t += o[b] * mt2[b]; }
    for (int b = M; b < 4; b++) L->occ[d][b] = 0;
    L->N[d] = N; L->s2[d] = s; L->t2[d] = t;
    L->key2idx[okey(o)] = d;
}
static void loc_build(Loc *L, int cp)
{
    L->d = 0; L->cap = cp;
    for (int k = 0; k < 81; k++) L->key2idx[k] = -1;
    int o[4];
    for (int N = 0; N <= cp; N++) {
        /* all occupations of M modes summing to N */
        int tot = 1; for (int a = 0; a < M; a++) tot *= (N + 1);
        for (int c = 0; c < tot; c++) {
            int x = c, s = 0;
            for (int a = 0; a < M; a++) { o[a] = x % (N + 1); x /= (N + 1); s += o[a]; }
            if (s == N) loc_add(L, o);
        }
    }
}
/* b^+_r b_m on local state k */
static int Eop(const Loc *L, int k, int r, int m, double *amp)
{
    if (L->occ[k][m] == 0) return -1;
    int o[4]; for (int a = 0; a < M; a++) o[a] = L->occ[k][a];
    double am = sqrt((double)o[m]); o[m]--;
    am *= sqrt((double)o[r] + 1); o[r]++;
    int k2 = L->key2idx[okey(o)];
    if (k2 < 0) return -1;
    *amp = am; return k2;
}
/* b_m (dir -1) or b^+_m (dir +1) */
static int bop(const Loc *L, int k, int m, int dir, double *amp)
{
    int o[4]; for (int a = 0; a < M; a++) o[a] = L->occ[k][a];
    if (dir < 0) { if (o[m] == 0) return -1; *amp = sqrt((double)o[m]); o[m]--; }
    else { if (L->N[k] >= L->cap) return -1; *amp = sqrt((double)o[m] + 1); o[m]++; }
    return L->key2idx[okey(o)];
}

/* ---------------- configurations and sectors ---------------- */
static inline void decode(int64_t c, int *dg) { for (int i = 0; i < n; i++) { dg[i] = (int)(c % radix[i]); c /= radix[i]; } }
static inline int64_t encode(const int *dg) { int64_t c = 0; for (int i = n - 1; i >= 0; i--) c = c * radix[i] + dg[i]; return c; }

typedef struct Ent Ent;
typedef struct {
    int N, s2, t2; long D; int64_t *codes;
    int ready; double *E, *U;
    int ncl; int *cl0, *clof;
    Ent **cache; int *ncache;     /* per eigenstate (single-state clusters) */
    double *obs;                  /* per eigenstate: cached observables (NOBS each), NaN if not yet */
    double *Gam; int gready;      /* per eigenstate: decay width into the dark continua (sharp golden rule) */
} Sec;
static Sec *secs; static int nsec;
static int *secmap; static int SW;

static int sec_find(int N, int s2, int t2)
{
    if (N < 0 || N > 2 * n) return -1;
    int off = 2 * n;
    if (s2 + off < 0 || s2 + off >= SW || t2 + off < 0 || t2 + off >= SW) return -1;
    return secmap[(N * SW + s2 + off) * SW + t2 + off];
}
static long sec_idx(const Sec *S, int64_t c)
{
    long lo = 0, hi = S->D - 1;
    while (lo <= hi) { long mid = (lo + hi) >> 1; if (S->codes[mid] == c) return mid; if (S->codes[mid] < c) lo = mid + 1; else hi = mid - 1; }
    return -1;
}

static void build_configs(void)
{
    int64_t tot = 1;
    for (int i = 0; i < n; i++) { site[i] = &loc[cap[i]]; radix[i] = site[i]->d; tot *= radix[i]; }
    SW = 4 * n + 1;
    secmap = malloc(sizeof(int) * (2 * n + 1) * SW * SW);
    for (int i = 0; i < (2 * n + 1) * SW * SW; i++) secmap[i] = -1;
    int cap_s = 64; secs = calloc(cap_s, sizeof(Sec)); nsec = 0;
    long *cnt = NULL; long *capv = NULL;
    int dg[16];
    /* two passes: count, then fill (codes ascending) */
    for (int pass = 0; pass < 2; pass++) {
        for (int64_t c = 0; c < tot; c++) {
            decode(c, dg);
            int N = 0, s = 0, t = 0;
            for (int i = 0; i < n; i++) { N += site[i]->N[dg[i]]; s += site[i]->s2[dg[i]]; t += site[i]->t2[dg[i]]; }
            int key = (N * SW + s + 2 * n) * SW + t + 2 * n;
            int id = secmap[key];
            if (id < 0) {
                if (nsec == cap_s) { cap_s *= 2; secs = realloc(secs, cap_s * sizeof(Sec)); memset(secs + nsec, 0, (cap_s - nsec) * sizeof(Sec)); }
                id = nsec++; secmap[key] = id; secs[id].N = N; secs[id].s2 = s; secs[id].t2 = t; secs[id].D = 0;
            }
            if (pass == 0) secs[id].D++;
            else secs[id].codes[cnt[id]++] = c;
        }
        if (pass == 0) {
            cnt = calloc(nsec, sizeof(long)); capv = calloc(nsec, sizeof(long));
            for (int s = 0; s < nsec; s++) { secs[s].codes = malloc(secs[s].D * sizeof(int64_t)); capv[s] = secs[s].D; }
        }
    }
    free(cnt); free(capv);
}

/* membership of every cell of a configuration */
static void members(const int *dg, int *Nc, int *mu)
{
    for (int i = 0; i < n; i++) Nc[i] = site[i]->N[dg[i]];
    for (int i = 0; i < n; i++) {
        if (!rc) { mu[i] = inL[i]; continue; }
        int m = (Nc[i] == z[i]);
        if (i > 0 && Nc[i - 1] < 1) m = 0;
        if (i == n - 1) m = 0;                 /* the lead is an empty neighbour */
        else if (Nc[i + 1] < 1) m = 0;
        mu[i] = m;
    }
}
/* saturation in c's sense, for the uniform readout */
static int nsat_of(const int *Nc)
{
    int s = 0;
    for (int i = 0; i < n - 1; i++) {
        if (Nc[i] != z[i]) continue;
        if (i > 0 && Nc[i - 1] < 1) continue;
        if (Nc[i + 1] < 1) continue;
        s++;
    }
    return s;
}

static double dtype(int ni, int nj) { if (ni == 2 && nj == 2) return DELTA; return OMEGA - EF[2]; }

/* ---------------- Hamiltonian of a sector (dense, column-major) ---------------- */
static double maxasym = 0;
static void build_H(const Sec *S, double *H)
{
    long D = S->D;
    memset(H, 0, sizeof(double) * D * D);
    int dg[16], dy[16], Nx[16], mux[16], Ny[16], muy[16];
    for (long x = 0; x < D; x++) {
        decode(S->codes[x], dg);
        members(dg, Nx, mux);
        double e = 0; for (int i = 0; i < n; i++) e += EF[Nx[i]];
        H[x + x * D] += e;
        for (int i = 0; i + 1 < n; i++) {
            int j = i + 1;
            /* exchange */
            int act = rb ? 1 : (rc ? (mux[i] || mux[j]) : (inL[i] && inL[j]));
            if (act && Nx[i] >= 1 && Nx[j] >= 1) {
                double cc = (Nx[i] == 1 && Nx[j] == 1) ? 1.0 / P2 : -G * G / dtype(Nx[i], Nx[j]);
                for (int r = 0; r < M; r++) for (int s = 0; s < M; s++) {
                    if (P[r][s] == 0) continue;
                    for (int m = 0; m < M; m++) for (int q = 0; q < M; q++) {
                        if (P[m][q] == 0) continue;
                        double ea, eb;
                        int a2 = Eop(site[i], dg[i], r, m, &ea); if (a2 < 0) continue;
                        int b2 = Eop(site[j], dg[j], s, q, &eb); if (b2 < 0) continue;
                        memcpy(dy, dg, sizeof(int) * n); dy[i] = a2; dy[j] = b2;
                        long y = sec_idx(S, encode(dy));
                        if (y < 0) { fprintf(stderr, "exchange left the sector\n"); exit(1); }
                        H[y + x * D] += cc * P[r][s] * P[m][q] * ea * eb;
                    }
                }
            }
            /* hops */
            for (int dir = 0; dir < 2; dir++) {
                int src = dir ? j : i, dst = dir ? i : j;
                for (int m = 0; m < M; m++) {
                    double ea, eb;
                    int s2_ = bop(site[src], dg[src], m, -1, &ea); if (s2_ < 0) continue;
                    int d2_ = bop(site[dst], dg[dst], m, +1, &eb); if (d2_ < 0) continue;
                    memcpy(dy, dg, sizeof(int) * n); dy[src] = s2_; dy[dst] = d2_;
                    members(dy, Ny, muy);
                    double t;
                    if (!rc) {
                        if (inL[i] && inL[j]) {
                            t = tin;
                            if (Nx[i] == 0 || Nx[j] == 0 || Ny[i] == 0 || Ny[j] == 0) t = 0;
                        } else if (!inL[i] && !inL[j]) t = tL;
                        else t = 0;
                    } else {
                        int any = mux[i] || mux[j] || muy[i] || muy[j];
                        if (any) {
                            t = tin;
                            if ((mux[i] && Ny[i] == 0) || (mux[j] && Ny[j] == 0) || (muy[i] && Nx[i] == 0) || (muy[j] && Nx[j] == 0)) t = 0;
                        } else t = tL;
                    }
                    if (t == 0) continue;
                    long y = sec_idx(S, encode(dy));
                    if (y < 0) { fprintf(stderr, "hop left the sector\n"); exit(1); }
                    H[y + x * D] += t * ea * eb;
                }
            }
        }
    }
    for (long x = 0; x < D; x++) for (long y = 0; y < x; y++) {
        double d = fabs(H[y + x * D] - H[x + y * D]);
        if (d > maxasym) maxasym = d;
    }
}

#define NOBS 4   /* w_lock, F, nsat, Nlight */
static int lockN, Lsec = -1;  /* lock sector (2k-2, 0, 0) */
static double *L0 = NULL;     /* |L0> in Lsec's occupation basis */
static int is_lockpat(const int *dg)
{
    for (int i = 0; i < n; i++) {
        int N = site[i]->N[dg[i]];
        int want = (i < kL) ? ((i == 0 || i == kL - 1) ? 1 : 2) : 0;
        if (N != want) return 0;
    }
    return 1;
}

static double resid_max = 0;
static int check_resid = 0;
static void sec_ready(int s)
{
    Sec *S = &secs[s];
    if (S->ready) return;
    long D = S->D; int Di = (int)D;
    double *H = malloc(sizeof(double) * D * D);
    build_H(S, H);
    double *Hc = NULL;
    if (check_resid) { Hc = malloc(sizeof(double) * D * D); memcpy(Hc, H, sizeof(double) * D * D); }
    S->E = malloc(sizeof(double) * D);
    int lwork = -1, liwork = -1, info, iwq; double wq;
    dsyevd_("V", "U", &Di, H, &Di, S->E, &wq, &lwork, &iwq, &liwork, &info);
    lwork = (int)wq + 1; liwork = iwq + 1;
    double *work = malloc(sizeof(double) * lwork); int *iwork = malloc(sizeof(int) * liwork);
    dsyevd_("V", "U", &Di, H, &Di, S->E, work, &lwork, iwork, &liwork, &info);
    if (info) { fprintf(stderr, "dsyevd info %d\n", info); exit(1); }
    free(work); free(iwork);
    S->U = H;
    if (Hc) {
        double *R = malloc(sizeof(double) * D * D);
        double one = 1, zero = 0;
        dgemm_("N", "N", &Di, &Di, &Di, &one, Hc, &Di, S->U, &Di, &zero, R, &Di);
        for (long k = 0; k < D; k++) for (long x = 0; x < D; x++) {
            double r = fabs(R[x + k * D] - S->E[k] * S->U[x + k * D]);
            if (r > resid_max) resid_max = r;
        }
        free(R); free(Hc);
    }
    /* clusters */
    S->cl0 = malloc(sizeof(int) * (D + 1)); S->clof = malloc(sizeof(int) * D);
    S->ncl = 0;
    for (long k = 0; k < D; k++) {
        if (k == 0 || S->E[k] - S->E[k - 1] > 1e-9) S->cl0[S->ncl++] = (int)k;
        S->clof[k] = S->ncl - 1;
    }
    S->cl0[S->ncl] = (int)D;
    S->cache = calloc(D, sizeof(Ent *)); S->ncache = calloc(D, sizeof(int));
    for (long k = 0; k < D; k++) S->ncache[k] = -1;
    S->obs = malloc(sizeof(double) * D * NOBS);
    for (long k = 0; k < D * NOBS; k++) S->obs[k] = NAN;
    S->ready = 1;
}

/* ---------------- observables of an occupation-basis vector ---------------- */
static void observe(int s, const double *psi, double *o)
{
    Sec *S = &secs[s];
    int dg[16], Nc[16];
    double wl = 0, ns_ = 0, nl = 0, F = 0;
    for (long x = 0; x < S->D; x++) {
        double p = psi[x] * psi[x];
        if (p == 0) continue;
        decode(S->codes[x], dg);
        for (int i = 0; i < n; i++) Nc[i] = site[i]->N[dg[i]];
        if (S->N == lockN && is_lockpat(dg)) wl += p;
        ns_ += p * nsat_of(Nc);
        int L = 0; for (int i = kL; i < n; i++) L += Nc[i];
        nl += p * L;
    }
    if (s == Lsec && L0) { double a = 0; for (long x = 0; x < S->D; x++) a += L0[x] * psi[x]; F = a * a; }
    o[0] = wl; o[1] = F; o[2] = ns_; o[3] = nl;
}

/* ---------------- events ---------------- */
enum { LOUT, LIN, SOUT, SIN, SFORM, SA, NTYPE };
static const char *tname[NTYPE] = {"light_out", "light_in", "space_out", "space_in", "side_form", "side_a"};
static int is_arrival(int t) { return t == LIN || t == SIN || t == SFORM || t == SA; }

struct Ent { int type, idx, tsec, tcl; double rate; double *coef; };

/* apply event operator (type, idx) to psi in sector s; returns target sector, fills out (target occupation basis) */

static int nu_partner(int nuv, double *p) { for (int a = 0; a < M; a++) if (P[a][nuv] != 0) { *p = P[a][nuv]; return a; } return -1; }

static int target_sector(int s, int type, int idx)
{
    Sec *S = &secs[s];
    int N = S->N, s2 = S->s2, t2 = S->t2;
    switch (type) {
    case LOUT: return sec_find(N - 1, s2 - ms2[idx], t2 - mt2[idx]);
    case LIN:  return sec_find(N + 1, s2 + ms2[idx], t2 + mt2[idx]);
    case SOUT: return sec_find(N - 2, s2, t2);
    case SIN:  return sec_find(N + 2, s2, t2);
    case SFORM: { double p; int mu = nu_partner(idx % M, &p); return sec_find(N - 1, s2 - ms2[mu], t2 - mt2[mu]); }
    case SA:    { double p; int mu = nu_partner(idx % M, &p); return sec_find(N + 1, s2 + ms2[mu], t2 + mt2[mu]); }
    }
    return -1;
}

static int apply_op(int s, int type, int idx, const double *psi, double *out, int ts)
{
    Sec *S = &secs[s], *T = &secs[ts];
    memset(out, 0, sizeof(double) * T->D);
    int dg[16], dy[16], Nx[16], mux[16], Ny[16], muy[16];
    int any = 0;
    for (long x = 0; x < S->D; x++) {
        double px = psi[x];
        if (px == 0) continue;
        decode(S->codes[x], dg);
        if (type == LOUT || type == LIN) {
            int c = n - 1; double a;
            int c2 = bop(site[c], dg[c], idx, type == LOUT ? -1 : +1, &a);
            if (c2 < 0) continue;
            memcpy(dy, dg, sizeof(int) * n); dy[c] = c2;
            long y = sec_idx(T, encode(dy)); if (y < 0) { fprintf(stderr, "op sector\n"); exit(1); }
            out[y] += a * px; any = 1;
        } else if (type == SOUT || type == SIN) {
            int i = idx, j = idx + 1;
            members(dg, Nx, mux);
            for (int a = 0; a < M; a++) for (int b = 0; b < M; b++) {
                if (P[a][b] == 0) continue;
                double ea, eb;
                int dir = type == SOUT ? -1 : +1;
                int a2 = bop(site[i], dg[i], a, dir, &ea); if (a2 < 0) continue;
                int b2 = bop(site[j], dg[j], b, dir, &eb); if (b2 < 0) continue;
                memcpy(dy, dg, sizeof(int) * n); dy[i] = a2; dy[j] = b2;
                members(dy, Ny, muy);
                int memlink = (mux[i] && mux[j]) || (muy[i] && muy[j]);
                int anym = mux[i] || mux[j] || muy[i] || muy[j];
                int cls = memlink ? 0 : (anym ? 1 : 2);
                int ok;
                if (type == SOUT) ok = (cls == 1) || (cls == 2 && rb);
                else ok = (cls == 1 && ra) || (cls == 2 && rb);
                if (!ok) continue;
                long y = sec_idx(T, encode(dy)); if (y < 0) { fprintf(stderr, "op sector A\n"); exit(1); }
                out[y] += P[a][b] * ea * eb * px; any = 1;
            }
        } else {   /* side: idx = i*M + nu; J_{i nu} = P[mu][nu] b_{i mu} */
            int i = idx / M, nuv = idx % M; double p = 0;
            int mu = nu_partner(nuv, &p);
            double a;
            int c2 = bop(site[i], dg[i], mu, type == SFORM ? -1 : +1, &a);
            if (c2 < 0) continue;
            memcpy(dy, dg, sizeof(int) * n); dy[i] = c2;
            long y = sec_idx(T, encode(dy)); if (y < 0) { fprintf(stderr, "op sector side\n"); exit(1); }
            out[y] += p * a * px; any = 1;
        }
    }
    return any;
}

/* two-continuum kernels, tabulated */
#define KN 40001
static double Kf[KN], Ka[KN], Klo, Khi, Kdx;
static double rhoS(double d) { double x = 4 * sk * sk - d * d; return x > 0 ? sqrt(x) / (2 * M_PI * sk * sk) : 0; }
static double rhoL(double e) { double x = 4 * tL * tL - e * e; return x > 0 ? sqrt(x) / (2 * M_PI * tL * tL) : 0; }
static void kernels(void)
{
    Klo = OMEGA - EF[1] - 2 * tL - 2 * sk - 0.01; Khi = OMEGA - EF[1] + 2 * tL + 2 * sk + 0.01;
    Kdx = (Khi - Klo) / (KN - 1);
    const int Q = 2000;
    for (int k = 0; k < KN; k++) {
        double dE = Klo + k * Kdx, sf = 0, sa = 0;
        for (int q = 0; q < Q; q++) {
            double th = (q + 0.5) * M_PI / Q, w = (2.0 / M_PI) * sin(th) * sin(th) * (M_PI / Q);
            double eps = 2 * tL * cos(th);              /* side formation: light in at eps */
            sf += w * rhoS(dE + EF[1] + eps - OMEGA);
            double del = 2 * sk * cos(th);              /* side a: space in at del */
            sa += w * rhoL(OMEGA + del - dE - EF[1]);   /* table variable = E_tgt - E_src */
        }
        Kf[k] = 2 * M_PI * G * G * sf; Ka[k] = 2 * M_PI * G * G * sa;
    }
}
static double kint(const double *K, double dE)
{
    double u = (dE - Klo) / Kdx; if (u <= 0 || u >= KN - 1) return 0;
    int k = (int)u; double f = u - k; return K[k] * (1 - f) + K[k + 1] * f;
}
#define EDGE 1e-9
/* the semicircle band sqrt(a^2 - w^2) convolved with a Lorentzian of half-width g:
   int sqrt(a^2-w^2) (g/pi)/((x-w)^2+g^2) dw = Im sqrt((x+ig)^2 - a^2) - g  (branch Im >= 0).
   g = 0: the sharp band; a transition within EDGE of a band edge has zero density of states (exact edge
   coincidences of the two-cell light region otherwise give round-off rates ~1e-8). */
static double semiL(double x, double a, double g)
{
    if (g <= 0) { if (fabs(x) >= a - EDGE) return 0; return sqrt(a * a - x * x); }
    double complex zz = x + I * g;
    double complex s = csqrt(zz * zz - a * a);
    if (cimag(s) < 0) s = -s;
    double r = cimag(s) - g;
    return r > 0 ? r : 0;
}
/* rate factor; dE = E_src - E_tgt (chain energy released); gam = (Gamma_src + Gamma_tgt)/2 */
static double fac(int type, double dE, double gam)
{
    switch (type) {
    case LOUT: return (1 + nL_env) * semiL(dE - EF[1], 2 * tL, gam);
    case LIN:  return nL_env * semiL(-dE - EF[1], 2 * tL, gam);
    case SOUT: return (1 + ns_env) * G * G / (sk * sk) * semiL(dE - OMEGA, 2 * sk, gam);
    case SIN:  return ns_env * G * G / (sk * sk) * semiL(-dE - OMEGA, 2 * sk, gam);
    case SFORM: return nL_env * (1 + ns_env) * kint(Kf, dE);
    case SA:    return ns_env * (1 + nL_env) * kint(Ka, -dE);
    }
    return 0;
}
/* the dark, sharp rate of a spontaneous event (for the widths) */
static double fsharp(int type, double dE)
{
    if (type == LOUT) return semiL(dE - EF[1], 2 * tL, 0);
    if (type == SOUT) return G * G / (sk * sk) * semiL(dE - OMEGA, 2 * sk, 0);
    return 0;
}

/* event list of the current run */
static int nev; static int evt[256], evi[256];
static void events_init(void)
{
    nev = 0;
    for (int m = 0; m < M; m++) { evt[nev] = LOUT; evi[nev++] = m; }
    if (nL_env > 0) for (int m = 0; m < M; m++) { evt[nev] = LIN; evi[nev++] = m; }
    for (int l = 0; l + 1 < n; l++) { evt[nev] = SOUT; evi[nev++] = l; }
    if (ns_env > 0 && (ra || rb)) for (int l = 0; l + 1 < n; l++) { evt[nev] = SIN; evi[nev++] = l; }
    if (side) for (int i = 1; i < kL - 1; i++) for (int v = 0; v < M; v++) {
        if (nL_env > 0) { evt[nev] = SFORM; evi[nev++] = i * M + v; }
        if (ns_env > 0 && ra) { evt[nev] = SA; evi[nev++] = i * M + v; }
    }
}

static double compl_max = 0, clus_max = 0;
static double *vbuf, *cbuf, *fbuf; static long vcap = 0;
static void ensure_buf(long D) { if (D > vcap) { vcap = D; vbuf = realloc(vbuf, sizeof(double) * D); cbuf = realloc(cbuf, sizeof(double) * D); fbuf = realloc(fbuf, sizeof(double) * D); } }

/* widths: Gamma_k = sum over spontaneous events (light out at the lead, space out where allowed) of the
   sharp dark golden-rule rates; computed per sector by dgemm, once per rule set */
static void sec_gamma(int s)
{
    Sec *S = &secs[s];
    if (S->gready) return;
    sec_ready(s);
    S->Gam = calloc(S->D, sizeof(double));
    int Ds = (int)S->D;
    for (int e = 0; e < M + n - 1; e++) {
        int type = e < M ? LOUT : SOUT, idx = e < M ? e : e - M;
        int ts = target_sector(s, type, idx);
        if (ts < 0) continue;
        sec_ready(ts);
        Sec *T = &secs[ts];
        int Dt = (int)T->D;
        double *V = malloc(sizeof(double) * (size_t)Dt * Ds), *C = malloc(sizeof(double) * (size_t)Dt * Ds);
        int any = 0;
        for (int k = 0; k < Ds; k++) any |= apply_op(s, type, idx, S->U + (size_t)k * Ds, V + (size_t)k * Dt, ts);
        if (any) {
            double one = 1, zero = 0;
            dgemm_("T", "N", &Dt, &Ds, &Dt, &one, T->U, &Dt, V, &Dt, &zero, C, &Dt);
            for (int k = 0; k < Ds; k++) {
                double g = 0;
                for (int q = 0; q < Dt; q++) { double c = C[q + (size_t)k * Dt]; if (c != 0) g += c * c * fsharp(type, S->E[k] - T->E[q]); }
                S->Gam[k] += g;
            }
        }
        free(V); free(C);
    }
    S->gready = 1;
}

/* compute the entries from (s, psi) with energy E; returns count, list malloc'd */
static void sec_gamma(int s);
static int compute_entries(int s, const double *psi, double E, double gsrc, Ent **outp)
{
    int cap_e = 64, ne = 0; Ent *L = malloc(sizeof(Ent) * cap_e);
    for (int e = 0; e < nev; e++) {
        int type = evt[e], idx = evi[e];
        int ts = target_sector(s, type, idx);
        if (ts < 0) continue;
        sec_ready(ts); sec_gamma(ts);
        Sec *T = &secs[ts];
        ensure_buf(T->D);
        if (!apply_op(s, type, idx, psi, vbuf, ts)) continue;
        int Di = (int)T->D, one_i = 1; double one = 1, zero = 0;
        dgemv_("T", &Di, &Di, &one, T->U, &Di, vbuf, &one_i, &zero, cbuf, &one_i);
        double vn = 0, cn = 0;
        for (long k = 0; k < T->D; k++) { vn += vbuf[k] * vbuf[k]; cn += cbuf[k] * cbuf[k]; }
        if (vn > 0) { double d = fabs(cn - vn) / vn; if (d > compl_max) compl_max = d; }
        for (int c = 0; c < T->ncl; c++) {
            int k0 = T->cl0[c], k1 = T->cl0[c + 1];
            double w = 0; for (int k = k0; k < k1; k++) w += cbuf[k] * cbuf[k];
            if (w < 1e-20 * vn) continue;          /* numerical zero */
            double r = 0, rsum = 0;
            for (int k = k0; k < k1; k++) {
                double f = fac(type, E - T->E[k], 0.5 * (gsrc + T->Gam[k]));
                fbuf[k - k0] = f;
                r += cbuf[k] * cbuf[k] * f;
                rsum += cbuf[k] * cbuf[k] * fac(type, E - T->E[k0], 0.5 * (gsrc + T->Gam[k]));
            }
            if (r <= 0) continue;
            /* rsum: the same with every member at the cluster's first energy (tests the clustering) */
            if (w > 1e-12 * vn) {
                double d = fabs(rsum - r) / r;
                if (d > clus_max) clus_max = d;
                if (d > 1e-6 && getenv("GXDBG"))
                    fprintf(stderr, "clus type=%s idx=%d size=%d dE0=%.15f dE1=%.15f w=%.3e r=%.3e rsum=%.3e\n",
                            tname[type], idx, k1 - k0, E - T->E[k0], E - T->E[k1 - 1], w, r, rsum);
            }
            if (ne == cap_e) { cap_e *= 2; L = realloc(L, sizeof(Ent) * cap_e); }
            Ent *en = &L[ne++];
            en->type = type; en->idx = idx; en->tsec = ts; en->tcl = c; en->rate = r; en->coef = NULL;
            if (k1 - k0 > 1) { en->coef = malloc(sizeof(double) * (k1 - k0)); for (int k = k0; k < k1; k++) en->coef[k - k0] = cbuf[k] * sqrt(fbuf[k - k0]); }
        }
    }
    *outp = L; return ne;
}

/* ---------------- RNG ---------------- */
static uint64_t splitmix64(uint64_t *x) { uint64_t z_ = (*x += 0x9e3779b97f4a7c15ULL); z_ = (z_ ^ (z_ >> 30)) * 0xbf58476d1ce4e5b9ULL; z_ = (z_ ^ (z_ >> 27)) * 0x94d049bb133111ebULL; return z_ ^ (z_ >> 31); }
static uint64_t rs[4];
static inline uint64_t rotl(uint64_t x, int k) { return (x << k) | (x >> (64 - k)); }
static uint64_t xnext(void) { uint64_t r = rotl(rs[1] * 5, 7) * 9, t = rs[1] << 17; rs[2] ^= rs[0]; rs[3] ^= rs[1]; rs[1] ^= rs[2]; rs[0] ^= rs[3]; rs[2] ^= t; rs[3] = rotl(rs[3], 45); return r; }
static double urand(void) { return ((xnext() >> 11) + 0.5) * (1.0 / 9007199254740992.0); }
static void rseed(uint64_t s) { uint64_t x = s * 0x2545F4914F6CDD1DULL + 12345; for (int i = 0; i < 4; i++) rs[i] = splitmix64(&x); }

/* ---------------- KMC ---------------- */
typedef struct { int sec, cl, ncoef; double *coef; } State;   /* coef: capacity gDmax */
static long gDmax = 0;

static void state_vec(const State *st, double *psi)
{
    Sec *S = &secs[st->sec];
    int k0 = S->cl0[st->cl];
    memset(psi, 0, sizeof(double) * S->D);
    for (int a = 0; a < st->ncoef; a++) { const double *u = S->U + (long)(k0 + a) * S->D; double c = st->coef[a]; for (long x = 0; x < S->D; x++) psi[x] += c * u[x]; }
}
static double state_E(const State *st) { Sec *S = &secs[st->sec]; return S->E[S->cl0[st->cl]]; }

static double *psibuf; static long psicap = 0;
static void ensure_psi(long D) { if (D > psicap) { psicap = D; psibuf = realloc(psibuf, sizeof(double) * D); } }

static void state_obs(const State *st, double *o)
{
    Sec *S = &secs[st->sec];
    int k0 = S->cl0[st->cl];
    if (st->ncoef == 1 && !isnan(S->obs[(long)k0 * NOBS])) { memcpy(o, S->obs + (long)k0 * NOBS, sizeof(double) * NOBS); return; }
    ensure_psi(S->D); state_vec(st, psibuf); observe(st->sec, psibuf, o);
    if (st->ncoef == 1) memcpy(S->obs + (long)k0 * NOBS, o, sizeof(double) * NOBS);
}

static long CACHE_MAX = 20000000L;   /* cached entries per process (key cache=); past it, rates are recomputed (speed only) */
static long cached_entries = 0;
static Ent *get_entries(const State *st, int *ne, int *owned)
{
    Sec *S = &secs[st->sec];
    int k0 = S->cl0[st->cl];
    if (st->ncoef == 1 && S->ncache[k0] >= 0) { *ne = S->ncache[k0]; *owned = 0; return S->cache[k0]; }
    ensure_psi(S->D); state_vec(st, psibuf);
    sec_gamma(st->sec);
    double gsrc = 0; for (int a = 0; a < st->ncoef; a++) gsrc += st->coef[a] * st->coef[a] * S->Gam[k0 + a];
    Ent *L; int c = compute_entries(st->sec, psibuf, state_E(st), gsrc, &L);
    if (st->ncoef == 1 && cached_entries + c <= CACHE_MAX) { S->cache[k0] = L; S->ncache[k0] = c; cached_entries += c; *owned = 0; }
    else *owned = 1;
    *ne = c; return L;
}
static void free_entries(Ent *L, int ne) { for (int i = 0; i < ne; i++) free(L[i].coef); free(L); }

/* choose a start cluster from an occupation-basis vector */
static void start_from(int s, const double *psi0, State *st)
{
    sec_ready(s);
    Sec *S = &secs[s];
    ensure_buf(S->D);
    int Di = (int)S->D, one_i = 1; double one = 1, zero = 0;
    dgemv_("T", &Di, &Di, &one, S->U, &Di, psi0, &one_i, &zero, cbuf, &one_i);
    double u = urand(), acc = 0, tot = 0;
    for (long k = 0; k < S->D; k++) tot += cbuf[k] * cbuf[k];
    int pick = S->ncl - 1;
    for (int c = 0; c < S->ncl; c++) {
        double w = 0; for (int k = S->cl0[c]; k < S->cl0[c + 1]; k++) w += cbuf[k] * cbuf[k];
        acc += w / tot;
        if (u < acc) { pick = c; break; }
    }
    st->sec = s; st->cl = pick; st->ncoef = S->cl0[pick + 1] - S->cl0[pick];
    double nr = 0; for (int a = 0; a < st->ncoef; a++) { st->coef[a] = cbuf[S->cl0[pick] + a]; nr += st->coef[a] * st->coef[a]; }
    nr = sqrt(nr); for (int a = 0; a < st->ncoef; a++) st->coef[a] /= nr;
}

typedef struct {
    double obs_end[NOBS], obs_avg[NOBS];
    long cnt[NTYPE]; double t; long steps; int arrivals; int stuck; double Nend;
    double t_refill;
} TrajOut;

static double ledger_max = 0;
static void run_traj(const State *st0, int dense, TrajOut *R)
{
    double *bufA = malloc(sizeof(double) * gDmax), *bufB = malloc(sizeof(double) * gDmax);
    State st = *st0; st.coef = bufA; memcpy(bufA, st0->coef, sizeof(double) * st0->ncoef);
    memset(R, 0, sizeof(*R)); R->t_refill = -1;
    double t = 0, Tend = dense ? Tmax : 1e300;
    int relax = 0;
    long N0 = secs[st.sec].N, Lout = 0, Lin = 0, Sout = 0, Sin = 0;
    double o[NOBS];
    for (long step = 0; step < 20000000; step++) {
        int ne, owned;
        Ent *L = get_entries(&st, &ne, &owned);
        double Rt = 0;
        for (int i = 0; i < ne; i++) if (!relax || !is_arrival(L[i].type)) Rt += L[i].rate;
        state_obs(&st, o);
        if (R->t_refill < 0 && o[0] >= 0.9) R->t_refill = t;
        if (Rt < Rmin) {                    /* nothing can fire */
            if (dense) { double a = fmax(t, Tmax / 2), b = Tmax; if (b > a) for (int q = 0; q < NOBS; q++) R->obs_avg[q] += o[q] * (b - a); t = Tmax; }
            if (owned) free_entries(L, ne);
            R->stuck = 1; break;
        }
        double dt = -log(urand()) / Rt;
        if (dense) {
            double a = fmax(t, Tmax / 2), b = fmin(t + dt, Tmax);
            if (b > a) for (int q = 0; q < NOBS; q++) R->obs_avg[q] += o[q] * (b - a);
            if (t + dt >= Tend) { t = Tend; if (owned) free_entries(L, ne); break; }
        } else if (t + dt > Tdark && !relax && nu == 0) { t = Tdark; if (owned) free_entries(L, ne); break; }
        t += dt;
        double u = urand() * Rt, acc = 0; int pick = -1;
        for (int i = 0; i < ne; i++) { if (relax && is_arrival(L[i].type)) continue; acc += L[i].rate; if (u < acc) { pick = i; break; } }
        if (pick < 0) { for (int i = ne - 1; i >= 0; i--) if (!relax || !is_arrival(L[i].type)) { pick = i; break; } }
        Ent *en = &L[pick];
        R->cnt[en->type]++;
        switch (en->type) {
        case LOUT: Lout++; break; case LIN: Lin++; break;
        case SOUT: Sout++; break; case SIN: Sin++; break;
        case SFORM: Lin++; Sout++; break; case SA: Sin++; Lout++; break;
        }
        if (is_arrival(en->type)) R->arrivals++;
        /* jump */
        Sec *T = &secs[en->tsec];
        State ns; ns.sec = en->tsec; ns.cl = en->tcl; ns.ncoef = T->cl0[en->tcl + 1] - T->cl0[en->tcl];
        ns.coef = (st.coef == bufA) ? bufB : bufA;
        if (ns.ncoef == 1) ns.coef[0] = 1;
        else { double nr = 0; for (int a = 0; a < ns.ncoef; a++) { ns.coef[a] = en->coef[a]; nr += ns.coef[a] * ns.coef[a]; } nr = sqrt(nr); for (int a = 0; a < ns.ncoef; a++) ns.coef[a] /= nr; }
        if (owned) free_entries(L, ne);
        st = ns;
        long Nc = secs[st.sec].N;
        double led = fabs((double)(Nc + Lout - Lin + 2 * Sout - 2 * Sin - N0));
        if (led > ledger_max) ledger_max = led;
        R->steps++;
        if (!dense && !relax && nu > 0 && R->arrivals >= Amax) relax = 1;
    }
    state_obs(&st, o);
    memcpy(R->obs_end, o, sizeof(o));
    if (dense) for (int q = 0; q < NOBS; q++) R->obs_avg[q] /= (Tmax / 2);
    else memcpy(R->obs_avg, o, sizeof(o));
    R->t = t; R->Nend = secs[st.sec].N;
    free(bufA); free(bufB);
}

/* ---------------- start vectors ---------------- */
/* lowest eigenstate of the c-off, b-off H restricted to configs passing pred, in sector s */
static double last_E = 0;
static double *lowest_restricted(int s, int (*pred)(const int *))
{
    int sra = ra, srb = rb, src_ = rc; ra = rb = rc = 0;
    Sec *S = &secs[s];
    double *H = malloc(sizeof(double) * S->D * S->D);
    build_H(S, H);
    ra = sra; rb = srb; rc = src_;
    long *ix = malloc(sizeof(long) * S->D); int m = 0;
    int dg[16];
    for (long x = 0; x < S->D; x++) { decode(S->codes[x], dg); if (pred(dg)) ix[m++] = x; }
    /* invariance check: no H element from the subspace to outside */
    double leak = 0;
    for (int a = 0; a < m; a++) for (long y = 0; y < S->D; y++) {
        int in = 0; for (int b = 0; b < m; b++) if (ix[b] == y) { in = 1; break; }
        if (!in) leak = fmax(leak, fabs(H[y + ix[a] * S->D]));
    }
    double *h = malloc(sizeof(double) * m * m);
    for (int a = 0; a < m; a++) for (int b = 0; b < m; b++) h[a + b * m] = H[ix[a] + ix[b] * S->D];
    double *w = malloc(sizeof(double) * m);
    int lwork = -1, liwork = -1, info, iwq; double wq;
    dsyevd_("V", "U", &m, h, &m, w, &wq, &lwork, &iwq, &liwork, &info);
    lwork = (int)wq + 1; liwork = iwq + 1;
    double *work = malloc(sizeof(double) * lwork); int *iwork = malloc(sizeof(int) * liwork);
    dsyevd_("V", "U", &m, h, &m, w, work, &lwork, iwork, &liwork, &info);
    double *v = calloc(S->D, sizeof(double));
    for (int a = 0; a < m; a++) v[ix[a]] = h[a];
    printf("start_restricted sector=(%d,%d,%d) dim=%d E=%.12f leak=%.3e gap=%.6e\n", S->N, S->s2, S->t2, m, w[0], leak, m > 1 ? w[1] - w[0] : 0.0);
    last_E = w[0];
    free(H); free(ix); free(h); free(w); free(work); free(iwork);
    return v;
}
static int pred_lock(const int *dg) { return is_lockpat(dg); }
static int pred_endvac(const int *dg)
{
    for (int i = 0; i < n; i++) {
        int N = site[i]->N[dg[i]];
        int want = (i < kL) ? ((i == 0) ? 1 : (i == kL - 1 ? 0 : 2)) : 0;
        if (N != want) return 0;
    }
    return 1;
}
static int pred_bulkvac(const int *dg)
{
    int ones = 0;
    for (int i = 0; i < n; i++) {
        int N = site[i]->N[dg[i]];
        if (i >= kL) { if (N) return 0; continue; }
        if (i == 0 || i == kL - 1) { if (N != 1) return 0; continue; }
        if (N == 1) ones++; else if (N != 2) return 0;
    }
    return ones == 1;
}

/* ---------------- runs ---------------- */
static void set_env(const char *e, double v)
{
    ns_env = nL_env = 0;
    if (!strcmp(e, "S")) ns_env = v;
    else if (!strcmp(e, "L")) nL_env = v;
    else if (!strcmp(e, "SL")) ns_env = nL_env = v;
}

static void summarize(const char *tag, TrajOut *R, int nt, int dense)
{
    double m[NOBS] = {0}, m2[NOBS] = {0}, me[NOBS] = {0};
    double cnt[NTYPE] = {0}, tt = 0, arr = 0, st = 0, fl = 0, rel = 0, both = 0, anys = 0, stuckq = 0, tref = 0; int nref = 0;
    double frel = 0; int nfrel = 0;
    for (int r = 0; r < nt; r++) {
        for (int q = 0; q < NOBS; q++) { m[q] += R[r].obs_avg[q]; m2[q] += R[r].obs_avg[q] * R[r].obs_avg[q]; me[q] += R[r].obs_end[q]; }
        for (int y = 0; y < NTYPE; y++) cnt[y] += R[r].cnt[y];
        tt += R[r].t; arr += R[r].arrivals; st += R[r].stuck;
        if (R[r].obs_avg[0] >= 0.9) fl += 1;
        if (R[r].obs_avg[0] >= 0.5) { frel += R[r].obs_avg[1] / R[r].obs_avg[0]; nfrel++; }
        if (R[r].obs_end[0] >= 0.9) rel += 1;
        if (R[r].cnt[LOUT] == 2 && R[r].cnt[SOUT] == 0 && R[r].Nend == 0) both += 1;
        if (R[r].cnt[SOUT] > 0) anys += 1;
        if (R[r].Nend > 0) stuckq += 1;
        if (R[r].t_refill >= 0) { tref += R[r].t_refill; nref++; }
    }
    printf("RUN %s ntraj=%d dense=%d\n", tag, nt, dense);
    const char *on[NOBS] = {"w_lock", "F", "nsat", "Nlight"};
    for (int q = 0; q < NOBS; q++) {
        double mu = m[q] / nt, var = m2[q] / nt - mu * mu; if (var < 0) var = 0;
        printf("  %-7s mean=%.6f se=%.6f final=%.6f\n", on[q], mu, sqrt(var / nt), me[q] / nt);
    }
    printf("  P(w_lock>=0.9)=%.4f P(final w_lock>=0.9)=%.4f F/w_lock(w>=0.5)=%.6f (n=%d)\n", fl / nt, rel / nt, nfrel ? frel / nfrel : -1.0, nfrel);
    printf("  P(both light out)=%.4f P(any space out)=%.4f P(quanta left in chain)=%.4f P(ended dark-stable)=%.4f\n", both / nt, anys / nt, stuckq / nt, st / nt);
    printf("  mean t=%.6g arrivals=%.3f  t_first(w_lock>=0.9)=%.6g (n=%d)\n", tt / nt, arr / nt, nref ? tref / nref : -1.0, nref);
    printf("  events/traj:");
    for (int y = 0; y < NTYPE; y++) printf(" %s=%.3f", tname[y], cnt[y] / nt);
    printf("\n");
    fflush(stdout);
}

/* rates depend on the environment: drop every cached entry list */
static void cache_clear(void)
{
    for (int s = 0; s < nsec; s++) {
        Sec *S = &secs[s];
        if (!S->ready) continue;
        for (long k = 0; k < S->D; k++) if (S->ncache[k] >= 0) { free_entries(S->cache[k], S->ncache[k]); S->cache[k] = NULL; S->ncache[k] = -1; }
    }
    cached_entries = 0;
}

static void run_set(const char *tag, int s0, const double *psi0, const char *e, double v, int dense)
{
    cache_clear();
    set_env(e, v); nu = v; events_init();
    TrajOut *R = calloc(ntraj > ntraj_dense ? ntraj : ntraj_dense, sizeof(TrajOut));
    double *stbuf = malloc(sizeof(double) * gDmax);
    clock_t c0 = clock();
    int nt = (dense && ntraj_dense > 0) ? ntraj_dense : ntraj;
    double Tsave = Tmax; if (dense && Tdense > 0) Tmax = Tdense;
    for (int r = 0; r < nt; r++) {
        rseed(seed0 * 1000003ULL + r);
        State st; st.coef = stbuf; start_from(s0, psi0, &st);
        run_traj(&st, dense, &R[r]);
    }
    summarize(tag, R, nt, dense);
    if (dense) printf("  T=%g\n", Tmax);
    Tmax = Tsave;
    printf("  cpu=%.1fs\n", (double)(clock() - c0) / CLOCKS_PER_SEC);
    free(R); free(stbuf);
}

static void static_report(void)
{
    sec_ready(Lsec);
    Sec *S = &secs[Lsec];
    set_env("dark", 0); nu = 0; events_init();
    ensure_psi(S->D);
    printf("STATIC lock sector (%d,%d,%d) D=%ld\n", S->N, S->s2, S->t2, S->D);
    double best = -1; long kb = -1;
    for (long k = 0; k < S->D; k++) { double a = 0; for (long x = 0; x < S->D; x++) a += L0[x] * S->U[x + k * S->D]; if (a * a > best) { best = a * a; kb = k; } }
    for (long k = 0; k < S->D && k < 8; k++) {
        double o[NOBS];
        observe(Lsec, S->U + k * S->D, o);
        sec_gamma(Lsec);
        Ent *L; int ne = compute_entries(Lsec, S->U + k * S->D, S->E[k], S->Gam[k], &L);
        double rt = 0; for (int i = 0; i < ne; i++) rt += L[i].rate;
        free_entries(L, ne);
        printf("  k=%ld E=%.10f w_lock=%.6f F=%.6f nsat=%.4f Nlight=%.4f rate_out=%.6e\n", k, S->E[k], o[0], o[1], o[2], o[3], rt);
    }
    {
        double o[NOBS];
        observe(Lsec, S->U + kb * S->D, o);
        Ent *L; int ne = compute_entries(Lsec, S->U + kb * S->D, S->E[kb], S->Gam[kb], &L);
        double rt = 0; for (int i = 0; i < ne; i++) rt += L[i].rate;
        free_entries(L, ne);
        printf("  maxF k=%ld E=%.10f w_lock=%.6f F=%.6f rate_out=%.6e  (E-E0=%.6e J0)\n", kb, S->E[kb], o[0], o[1], rt, (S->E[kb] - S->E[0]) / J0);
    }
    /* is L0 an eigenstate of this rule set? */
    double *H = malloc(sizeof(double) * S->D * S->D);
    build_H(S, H);
    double e = 0, r = 0;
    for (long x = 0; x < S->D; x++) { double hv = 0; for (long y = 0; y < S->D; y++) hv += H[x + y * S->D] * L0[y]; e += L0[x] * hv; }
    for (long x = 0; x < S->D; x++) { double hv = 0; for (long y = 0; y < S->D; y++) hv += H[x + y * S->D] * L0[y]; r += (hv - e * L0[x]) * (hv - e * L0[x]); }
    printf("  L0: <H>=%.12f residual=%.3e\n", e, sqrt(r));
    sec_gamma(Lsec);
    { double gl = 0, gmax = 0; Sec *T = &secs[Lsec];
      for (long k = 0; k < T->D; k++) { double c = 0; for (long x = 0; x < T->D; x++) c += L0[x] * T->U[x + k * T->D]; gl += c * c * T->Gam[k]; if (T->Gam[k] > gmax) gmax = T->Gam[k]; }
      printf("  widths: <L0|Gamma|L0>=%.6e  max Gamma in sector=%.6e  Gamma(k=0)=%.6e\n", gl, gmax, T->Gam[0]); }
    free(H);
}

/* ---------------- main ---------------- */
static void parse(int argc, char **argv)
{
    for (int i = 1; i < argc; i++) {
        char *eq = strchr(argv[i], '='); if (!eq) continue;
        *eq = 0; const char *k = argv[i], *v = eq + 1;
        if (!strcmp(k, "M")) M = atoi(v);
        else if (!strcmp(k, "n")) n = atoi(v);
        else if (!strcmp(k, "k")) kL = atoi(v);
        else if (!strcmp(k, "rules")) { ra = strchr(v, 'a') != NULL; rb = strchr(v, 'b') != NULL; rc = strchr(v, 'c') != NULL; }
        else if (!strcmp(k, "side")) side = atoi(v);
        else if (!strcmp(k, "gr")) gr = atof(v);
        else if (!strcmp(k, "tin")) tin_J0 = atof(v);
        else if (!strcmp(k, "tL")) tL = atof(v);
        else if (!strcmp(k, "sk")) sk = atof(v);
        else if (!strcmp(k, "exp")) snprintf(exp_, sizeof(exp_), "%s", v);
        else if (!strcmp(k, "ntraj")) ntraj = atoi(v);
        else if (!strcmp(k, "Amax")) Amax = atoi(v);
        else if (!strcmp(k, "T")) Tmax = atof(v);
        else if (!strcmp(k, "ntd")) ntraj_dense = atoi(v);
        else if (!strcmp(k, "cache")) CACHE_MAX = atol(v);
        else if (!strcmp(k, "Td")) Tdense = atof(v);
        else if (!strcmp(k, "seed")) seed0 = strtoull(v, NULL, 10);
        else if (!strcmp(k, "resid")) check_resid = atoi(v);
        else { fprintf(stderr, "unknown key %s\n", k); exit(1); }
    }
}

int main(int argc, char **argv)
{
    parse(argc, argv);
    constants(gr);
    tin = tin_J0 * J0;
    if (M == 2) { memset(P, 0, sizeof(P)); P[0][1] = 1; P[1][0] = -1; ms2[0] = 1; ms2[1] = -1; mt2[0] = mt2[1] = 0; }
    else {
        double e[2][2] = {{0, 1}, {-1, 0}};
        int S_[4] = {1, 1, -1, -1}, T_[4] = {1, -1, 1, -1};
        for (int a = 0; a < 4; a++) for (int b = 0; b < 4; b++) P[a][b] = e[a / 2][b / 2] * e[a % 2][b % 2];
        for (int a = 0; a < 4; a++) { ms2[a] = S_[a]; mt2[a] = T_[a]; }
    }
    P2 = 0; for (int a = 0; a < M; a++) for (int b = 0; b < M; b++) P2 += P[a][b] * P[a][b];
    loc_build(&loc[1], 1); loc_build(&loc[2], 2);
    for (int i = 0; i < n; i++) {
        z[i] = (i == 0) ? 1 : 2;
        inL[i] = i < kL;
        if (!rc) cap[i] = (i < kL) ? ((i == 0 || i == kL - 1) ? 1 : 2) : 2;
        else cap[i] = z[i];
    }
    /* the configuration space must hold the lock pattern under every rule set:
       c on has cap 2 at cell kL-1 (z = 2); c off has cap 1.  Use the same
       per-cell caps as the rule set; the lock pattern is in both.  */
    build_configs();
    kernels();
    lockN = 2 * kL - 2;
    Lsec = sec_find(lockN, 0, 0);
    printf("# gx M=%d n=%d k=%d rules=%s%s%s%s side=%d g/D=%g tin=%g J0 tL=%g sk=%g exp=%s nsec=%d J0=%.10f\n",
           M, n, kL, ra ? "a" : "", rb ? "b" : "", rc ? "c" : "", (ra || rb || rc) ? "" : "none", side, gr, tin_J0, tL, sk, exp_, nsec, J0);
    long Dmax = 0; for (int s = 0; s < nsec; s++) if (secs[s].D > Dmax) Dmax = secs[s].D;
    gDmax = Dmax;
    printf("# caps:"); for (int i = 0; i < n; i++) printf(" %d", cap[i]); printf("  Dmax=%ld\n", Dmax);
    fflush(stdout);

    /* |L0>: built in the c-off cap layout; under c on the cap layout differs, so build L0 there by the same predicate */
    if (side && rc) { fprintf(stderr, "side coupling is c-off only (PREREG X5)\n"); return 1; }
    L0 = lowest_restricted(Lsec, pred_lock);
    double EL0 = last_E;

    if (!strcmp(exp_, "t0")) {
        check_resid = 1;
        for (int s = 0; s < nsec; s++) sec_ready(s);
        printf("T0 maxasym=%.3e resid_max=%.3e\n", maxasym, resid_max);
        static_report();
        /* completeness over a sweep of entries from every eigenstate of the lock sector, in env SL dense with side */
        set_env("SL", 1.0); nu = 1; events_init();
        Sec *S = &secs[Lsec];
        sec_gamma(Lsec);
        for (long k = 0; k < S->D; k++) { Ent *L; int ne = compute_entries(Lsec, S->U + k * S->D, S->E[k], S->Gam[k], &L); free_entries(L, ne); }
        printf("T0 completeness_max=%.3e cluster_rate_max=%.3e\n", compl_max, clus_max);
        /* the broadened band against direct quadrature */
        { double worst = 0; double xs[6] = {-2.5, -1.7, -0.3, 0.0, 1.1, 3.0}, gs[4] = {1e-3, 0.05, 0.5, 2.0};
          for (int i = 0; i < 6; i++) for (int j = 0; j < 4; j++) {
              double a = 1.8, x = xs[i], g = gs[j], sum = 0; int Q = 2000000;
              for (int q = 0; q < Q; q++) { double th = (q + 0.5) * M_PI / Q, w = a * cos(th), sn = a * sin(th);
                  sum += sn * sn * (g / M_PI) / ((x - w) * (x - w) + g * g) * (M_PI / Q); }
              double d = fabs(sum - semiL(x, a, g)) / fmax(sum, 1e-300); if (d > worst) worst = d; }
          printf("T0 semicircle*Lorentzian closed form vs quadrature: max rel diff %.3e\n", worst); }
        /* side kernels: formation from f and side a into f are one double integral, Kf(x) = Ka(x) */
        { double mx = 0, dk = 0; for (int q = 0; q < KN; q++) { mx = fmax(mx, Kf[q]); dk = fmax(dk, fabs(Kf[q] - Ka[q])); }
          printf("T0 side kernels max|Kf-Ka|/max Kf = %.3e (max Kf = %.6e at g^2 = %.6e)\n", dk / mx, mx, G * G); }
        double EfL = 0; for (int i = 0; i < kL; i++) EfL += EF[(i == 0 || i == kL - 1) ? 1 : 2];
        printf("T0 lock E0(c-off,b-off restricted) - sum E_f = %.12f rate  (%.8f J0)\n", EL0 - EfL, (EL0 - EfL) / J0);
        return 0;
    }
    if (!strcmp(exp_, "static")) { static_report(); return 0; }

    int all = !strcmp(exp_, "all");
    if (all || !strcmp(exp_, "static")) static_report();
    if (all || !strcmp(exp_, "gather")) {
        char tag[64];
        int s0 = sec_find(0, 0, 0); double *psi0 = calloc(1, sizeof(double)); psi0[0] = 1;
        const char *envs[3] = {"S", "L", "SL"};
        for (int e = 0; e < 3; e++) {
            snprintf(tag, sizeof(tag), "X1_%s_dilute", envs[e]); run_set(tag, s0, psi0, envs[e], 1e-6, 0);
            snprintf(tag, sizeof(tag), "X1_%s_dense", envs[e]);  run_set(tag, s0, psi0, envs[e], 1.0, 1);
        }
    }
    if (all || !strcmp(exp_, "hold")) {
        run_set("X2_dark", Lsec, L0, "dark", 0, 0);
        run_set("X2_S_dense", Lsec, L0, "S", 1.0, 1);
        run_set("X2_SL_dense", Lsec, L0, "SL", 1.0, 1);
    }
    if (all || !strcmp(exp_, "light")) {
        int s = sec_find(2, ms2[0] + ms2[M - 1], mt2[0] + mt2[M - 1]);
        sec_ready(s);
        int dg[16] = {0};
        int o[4] = {0, 0, 0, 0};
        o[0] = 1; dg[n - 2] = site[n - 2]->key2idx[okey(o)];
        o[0] = 0; o[M - 1] = 1; dg[n - 1] = site[n - 1]->key2idx[okey(o)];
        long x = sec_idx(&secs[s], encode(dg));
        double *v = calloc(secs[s].D, sizeof(double)); v[x] = 1;
        run_set("X3_dark", s, v, "dark", 0, 0);
    }
    if (all || !strcmp(exp_, "vac")) {
        int s = sec_find(2 * kL - 3, 1, M == 4 ? 1 : 0);
        double *ve = lowest_restricted(s, pred_endvac);
        double *vb = lowest_restricted(s, pred_bulkvac);
        if (!side) {
            run_set("X5e_S_dilute", s, ve, "S", 1e-6, 0);
            run_set("X5e_S_dense", s, ve, "S", 1.0, 1);
        }
        run_set(side ? "X5b_side_S_dilute" : "X5b_S_dilute", s, vb, "S", 1e-6, 0);
        run_set(side ? "X5b_side_S_dense" : "X5b_S_dense", s, vb, "S", 1.0, 1);
        if (side) run_set("X5b_side_SL_dense_supplement", s, vb, "SL", 1.0, 1);
    }
    printf("CHECK ledger_max=%.3e completeness_max=%.3e cluster_rate_max=%.3e maxasym=%.3e\n", ledger_max, compl_max, clus_max, maxasym);
    return 0;
}
