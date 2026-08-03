/* harmfrust.c — v89 — does harmonic size-quantisation survive FRUSTRATION?
 *
 * Standalone. Build: gcc -O2 -Wall -Wextra -o harmfrust harmfrust.c -lm
 *
 * THE QUESTION THAT CAN KILL IDEA (1)
 * -----------------------------------
 * `harmgrow.c` measured a clean staircase — a cell's size locks to
 * 1/1, 2/1, 3/1, 4/1, 5/1, 6/1 as its preferred size sweeps — on a RING.
 * On a ring the driven cell has TWO neighbours and everything else is
 * uniform, so satisfying the comb is trivially possible.
 *
 * In a foam a cell has degree ~7 and must hold intervals with all of
 * them AT ONCE. Rungs are a discrete set, so generically you cannot put
 * seven pair-ratios on low-height rungs simultaneously — the sizes would
 * all have to lie in a common small-integer lattice. If the staircase is
 * an artifact of a degenerate topology it will die as degree rises, and
 * idea (1) is rejected as a bounding mechanism.
 *
 * But if it dies in a particular WAY — if consonance can be satisfied
 * only inside small groups — that failure is itself the bound being
 * looked for: consonant domains would have a maximum size fixed by
 * combinatorics, not by any length. That is "harmonics within,
 * dissonance between" as a mechanism rather than a slogan, and it would
 * be a better answer than the comb-span bound harmgrow found.
 *
 * So this measures BOTH: does quantisation survive, and if not, is there
 * a characteristic consonant-domain size?
 *
 * Same potential as harmgrow:
 *   U = sum_edges D(u_ij) + sum_i (kE/2)(x_i - x*_i)^2,  u = log2(r_j/r_i)
 *   D(u) = min over rungs [ (u-u_k)^2/(2 sigma^2) + lambda H_k ]
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAXDEG 24

typedef struct { int p, q; double u, H; } Rung;
static Rung RUNGS[8192];
static int NRUNG = 0;

static int igcd(int a, int b) { while (b) { int t = a % b; a = b; b = t; } return a; }

static void build_rungs(double climb)
{
    NRUNG = 0;
    /* climb < 0 means "unison only" — the negative control. The height
     * test below rejects even 1/1 (0 > -1), which left NRUNG = 0: the
     * energy became 3e300, the backtracking comparison degenerated to
     * equality in floating point, and the relaxation exited after one
     * iteration. The control then read a perfectly flat response and
     * looked like it had passed. It had not run. */
    if (climb < 0) {
        RUNGS[0].p = 1; RUNGS[0].q = 1; RUNGS[0].u = 0.0; RUNGS[0].H = 0.0;
        NRUNG = 1;
        return;
    }
    int lim = (int)(pow(2.0, climb > 0 ? climb : 0.0) + 0.5);
    if (lim < 1) lim = 1;
    for (int p = 1; p <= lim; p++)
        for (int q = 1; q <= lim; q++) {
            if (p * q > lim || igcd(p, q) != 1) continue;
            double H = log2((double)p * (double)q);
            if (H > climb + 1e-12 || NRUNG >= 8192) continue;
            RUNGS[NRUNG].p = p; RUNGS[NRUNG].q = q;
            RUNGS[NRUNG].u = log2((double)p / (double)q);
            RUNGS[NRUNG].H = H; NRUNG++;
        }
}

/* D(u) and dD/du tabulated once per comb: the inner loop was O(NRUNG)
 * per gradient evaluation, which is ~211 rung comparisons per edge per
 * iteration. Table lookup makes it O(1). The exact routine is kept for
 * diagnostics, where accuracy matters and cost does not. */
#define TBN 400001
#define TBLO (-10.0)
#define TBHI (10.0)
static double TBD[TBN], TBG[TBN];
static double diss_exact(double u, double sg, double lm, int *kb);
static void build_table(double sg, double lm)
{
    for (int i = 0; i < TBN; i++) {
        double u = TBLO + (TBHI - TBLO) * i / (TBN - 1);
        int k; TBD[i] = diss_exact(u, sg, lm, &k);
        TBG[i] = (u - RUNGS[k].u) / (sg * sg);
    }
}
static inline int tbidx(double u)
{
    double f = (u - TBLO) / (TBHI - TBLO) * (TBN - 1);
    int i = (int)(f + 0.5);
    return i < 0 ? 0 : (i >= TBN ? TBN - 1 : i);
}
static double diss_exact(double u, double sg, double lm, int *kb)
{
    double best = 1e300; int b = 0;
    for (int k = 0; k < NRUNG; k++) {
        double d = u - RUNGS[k].u;
        double v = d * d / (2 * sg * sg) + lm * RUNGS[k].H;
        if (v < best) { best = v; b = k; }
    }
    if (kb) *kb = b;
    return best;
}
static double diss(double u, double sg, double lm, int *kb)
{ if (kb) return diss_exact(u, sg, lm, kb); return TBD[tbidx(u)]; }
static double dgrad(double u, double sg, double lm)
{ (void)sg; (void)lm; return TBG[tbidx(u)]; }

typedef struct {
    int N, ne, *ei, *ej, *nd;
    double *x, *xt, sigma, lambda, kE;
} Sys;

static unsigned long long RS;
static double rnd(void)
{ RS ^= RS << 13; RS ^= RS >> 7; RS ^= RS << 17;
  return (double)((RS >> 11) & 0xFFFFFFFFFFFFFULL) / 9007199254740992.0; }

static Sys *sys_new(int N)
{
    Sys *s = calloc(1, sizeof(Sys));
    s->N = N; s->sigma = 0.05; s->lambda = 0.10; s->kE = 1.0;
    s->x = calloc(N, sizeof(double));
    s->xt = calloc(N, sizeof(double));
    s->nd = calloc(N, sizeof(int));
    s->ei = calloc((size_t)N * MAXDEG, sizeof(int));
    s->ej = calloc((size_t)N * MAXDEG, sizeof(int));
    return s;
}
static void sys_free(Sys *s)
{ free(s->x); free(s->xt); free(s->nd); free(s->ei); free(s->ej); free(s); }

/* random regular-ish graph of the given degree; ring backbone for
 * connectivity, then random chords. deg=2 reproduces harmgrow's ring. */
static void build_graph(Sys *s, int deg, unsigned long long seed)
{
    int N = s->N;
    RS = seed;
    s->ne = 0;
    memset(s->nd, 0, N * sizeof(int));
    for (int i = 0; i < N; i++) {
        int j = (i + 1) % N;
        s->ei[s->ne] = i; s->ej[s->ne] = j; s->ne++;
        s->nd[i]++; s->nd[j]++;
    }
    int want = N * (deg - 2) / 2;
    for (int a = 0; a < want * 20 && s->ne < N * MAXDEG / 2; a++) {
        int cnt = 0;
        for (int i = 0; i < N; i++) cnt += s->nd[i];
        if (cnt >= N * deg) break;
        int i = (int)(rnd() * N), j = (int)(rnd() * N);
        if (i == j || s->nd[i] >= deg || s->nd[j] >= deg) continue;
        int dup = 0;
        for (int e = 0; e < s->ne; e++)
            if ((s->ei[e] == i && s->ej[e] == j) || (s->ei[e] == j && s->ej[e] == i))
                { dup = 1; break; }
        if (dup) continue;
        s->ei[s->ne] = i; s->ej[s->ne] = j; s->ne++;
        s->nd[i]++; s->nd[j]++;
    }
}

static double uij(Sys *s, int e)
{ return (s->x[s->ej[e]] - s->x[s->ei[e]]) / M_LN2; }

static double energy(Sys *s)
{
    double U = 0;
    for (int e = 0; e < s->ne; e++)
        U += diss(uij(s, e), s->sigma, s->lambda, NULL);
    for (int i = 0; i < s->N; i++) {
        double d = s->x[i] - s->xt[i];
        U += 0.5 * s->kE * d * d;
    }
    return U;
}

static void grad(Sys *s, double *G)
{
    int N = s->N;
    memset(G, 0, N * sizeof(double));
    for (int e = 0; e < s->ne; e++) {
        double dd = dgrad(uij(s, e), s->sigma, s->lambda);
        G[s->ei[e]] += dd * (-1.0 / M_LN2);
        G[s->ej[e]] += dd * (+1.0 / M_LN2);
    }
    for (int i = 0; i < N; i++) G[i] += s->kE * (s->x[i] - s->xt[i]);
    double m = 0;
    for (int i = 0; i < N; i++) m += G[i];
    m /= N;
    for (int i = 0; i < N; i++) G[i] -= m;
}

static double relax(Sys *s, int iters, double lr, double *gm)
{
    int N = s->N;
    double *G = malloc(N * sizeof(double)), *xb = malloc(N * sizeof(double));
    double U = energy(s), step = lr;
    for (int it = 0; it < iters; it++) {
        grad(s, G);
        memcpy(xb, s->x, N * sizeof(double));
        double Un = 0; int ok = 0;
        for (int tr = 0; tr < 40; tr++) {
            for (int i = 0; i < N; i++) s->x[i] = xb[i] - step * G[i];
            double m = 0;
            for (int i = 0; i < N; i++) m += s->x[i];
            m /= N;
            for (int i = 0; i < N; i++) s->x[i] -= m;
            Un = energy(s);
            if (isfinite(Un) && Un <= U) { ok = 1; break; }
            step *= 0.5;
        }
        if (!ok) { memcpy(s->x, xb, N * sizeof(double)); break; }
        if (U - Un < 1e-15 * (1 + fabs(U))) { U = Un; break; }
        U = Un; step *= 1.1; if (step > 1.0) step = 1.0;
    }
    grad(s, G);
    *gm = 0;
    for (int i = 0; i < N; i++) if (fabs(G[i]) > *gm) *gm = fabs(G[i]);
    U = energy(s);
    free(G); free(xb);
    return U;
}

/* rung-seeded search on cell 0, as harmgrow requires */
static double solve(Sys *s, double *gm)
{
    int N = s->N;
    double *best = malloc(N * sizeof(double));
    double Ub = 0; int have = 0; *gm = 0;
    for (int k = 0; k < NRUNG; k++) {
        for (int i = 0; i < N; i++) s->x[i] = (i == 0) ? RUNGS[k].u * M_LN2 : 0.0;
        double m = 0;
        for (int i = 0; i < N; i++) m += s->x[i];
        m /= N;
        for (int i = 0; i < N; i++) s->x[i] -= m;
        double g2, U = relax(s, 800, 0.002, &g2);
        if (isfinite(U) && (!have || U < Ub))
            { Ub = U; *gm = g2; have = 1; memcpy(best, s->x, N * sizeof(double)); }
    }
    memcpy(s->x, best, N * sizeof(double));
    free(best);
    return Ub;
}

/* fraction of edges sitting on a rung, and the largest connected set of
 * cells joined only by locked edges — the consonant domain */
static double locked_and_domain(Sys *s, double tol, int *maxdom, double *mH)
{
    int N = s->N, lock = 0;
    double sH = 0;
    int *par = malloc(N * sizeof(int));
    for (int i = 0; i < N; i++) par[i] = i;
    for (int e = 0; e < s->ne; e++) {
        int k; double u = uij(s, e);
        diss(u, s->sigma, s->lambda, &k);
        if (fabs(u - RUNGS[k].u) < tol) {
            lock++; sH += RUNGS[k].H;
            int a = s->ei[e], b = s->ej[e];
            while (par[a] != a) a = par[a];
            while (par[b] != b) b = par[b];
            if (a != b) par[a] = b;
        }
    }
    int *cnt = calloc(N, sizeof(int));
    for (int i = 0; i < N; i++) {
        int a = i; while (par[a] != a) a = par[a];
        cnt[a]++;
    }
    *maxdom = 0;
    for (int i = 0; i < N; i++) if (cnt[i] > *maxdom) *maxdom = cnt[i];
    free(par); free(cnt);
    if (mH) *mH = lock ? sH / lock : -1;
    return s->ne ? (double)lock / s->ne : 0;
}

int main(void)
{
    const int N = 48;
    printf("# harmfrust — does harmonic size-quantisation survive degree?\n");
    printf("# N=%d, comb_limit=6, sigma=0.05, lambda=0.10, kE=1\n", N);
    printf("# deg=2 is harmgrow's ring (the control that must reproduce it)\n\n");

    build_rungs(6.0);
    build_table(0.05, 0.10);

    printf("== F1. FRUSTRATION AT FIXED COMPETITION RATIO ==\n");
    printf("The comb term grows with degree (one D per edge) while a fixed\n");
    printf("kE does not, so comparing degrees at fixed kE confounds\n");
    printf("frustration with total coupling: at kE=1 every degree collapses\n");
    printf("to the unison crystal (locked=1.000, mean_H=0) and nothing is\n");
    printf("measured. kE is scaled with degree to hold the ratio fixed, and\n");
    printf("EVERY cell is given its own preferred size so the uniform\n");
    printf("solution is not trivially optimal.\n");
    printf(" degree     kE   locked_frac  mean_H  max_domain  frac_at_unison\n");
    {
        int degs[] = {2, 3, 4, 6, 8, 10, 12};
        for (unsigned d = 0; d < 7; d++) {
            double kE = 12.0 * degs[d];
            Sys *s = sys_new(N);
            build_graph(s, degs[d], 20260802);
            s->kE = kE;
            RS = 424242;
            for (int i = 0; i < N; i++) s->xt[i] = 3.0 * M_LN2 * (rnd() - 0.5);
            double gm;
            solve(s, &gm);
            double mH; int md;
            double lf = locked_and_domain(s, 1e-3, &md, &mH);
            int uni = 0;
            for (int e = 0; e < s->ne; e++) {
                int k; diss(uij(s, e), s->sigma, s->lambda, &k);
                if (RUNGS[k].p == 1 && RUNGS[k].q == 1) uni++;
            }
            printf("%7d  %6.0f  %11.3f  %6.2f  %10d  %14.3f\n",
                   degs[d], kE, lf, mH, md, (double)uni / s->ne);
            sys_free(s);
        }
    }
    printf("\nreading: if consonance can be satisfied at every degree,\n");
    printf("locked_frac stays ~1 and mean_H stays low — frustration is not\n");
    printf("a barrier and idea (1) survives. If locked_frac falls and mean_H\n");
    printf("climbs with degree, generic sizes CANNOT all sit on low rungs at\n");
    printf("once, and max_domain is then the consonant-group size — the\n");
    printf("bound worth having. frac_at_unison near 1 means the run collapsed\n");
    printf("to the crystal and measured nothing either way.\n\n");

    printf("== F2. DOMAIN SIZE vs comb_limit (is the bound countable?) ==\n");
    printf("comb_limit  nrungs  locked_frac  max_domain  mean_H\n");
    {
        for (double cl = 1; cl <= 8.001; cl += 1.0) {
            build_rungs(cl);
            build_table(0.05, 0.10);
            Sys *s = sys_new(N);
            build_graph(s, 7, 20260802);
            s->kE = 12.0 * 7;
            RS = 424242;
            for (int i = 0; i < N; i++) s->xt[i] = 3.0 * M_LN2 * (rnd() - 0.5);
            double gm;
            solve(s, &gm);
            double mH; int md;
            double lf = locked_and_domain(s, 1e-3, &md, &mH);
            printf("%10.0f  %6d  %11.3f  %10d  %6.2f\n", cl, NRUNG, lf, md, mH);
            sys_free(s);
        }
    }
    printf("\nreading: a max_domain that grows with comb_limit and saturates\n");
    printf("is a countable particle scale. One that is always N (or always 1)\n");
    printf("means consonance is not selecting anything.\n");
    return 0;
}
