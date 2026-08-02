/* harmgrow.c — v89 — "harmonics within, dissonance between":
 * can a cell's MORPHOLOGY BOUND come out of the consonance comb instead
 * of being imposed as a constant?
 *
 * Standalone. Does not touch cellfab.c. Build:
 *   gcc -O2 -Wall -Wextra -o harmgrow harmgrow.c -lm
 *
 * THE QUESTION
 * ------------
 * Idea (1) of the 2026-08-02 design thread: cells are not held in place,
 * they GROW, with a bound on morphology. The objection to that idea is
 * that a fixed maximum extent is a background — a length that persists,
 * is never converted, and that nothing can act on. This file tests the
 * proposed escape: let the bound come from HARMONY. A cell may grow as
 * long as growing keeps its intervals with its neighbours on the comb;
 * the bound is where the next increment would cost more dissonance than
 * the growth is worth.
 *
 * THE MODEL (deliberately minimal — a mechanism test, not a substrate)
 * -------------------------------------------------------------------
 * N cells with extents r_i > 0. Work in x_i = ln r_i, with sum x_i = 0.
 * That constraint is scale fixing, and it is REQUIRED: the potential
 * depends on ratios plus a drive that is linear in log-size, so without
 * it the whole configuration inflates without bound. Fixing the
 * geometric mean says "space is redistributed, not created" — nothing
 * here privileges any particular absolute size.
 *
 * Pitch is set by extent, w_i = 1/r_i: a bigger cavity rings lower. So
 * the interval between two cells IS their size ratio, and "consonant
 * neighbours" and "commensurate sizes" become the same statement. That
 * identification is the whole reason this idea might work.
 *
 * Pair term — dissonance of the interval, over the comb:
 *   u_ij   = log2(r_j/r_i) = (x_j - x_i)/ln2
 *   D(u)   = min over rungs k of [ (u-u_k)^2/(2 sigma^2) + lambda*H_k ]
 *   rung k = p/q reduced, Tenney height H = log2(p*q) <= comb_limit
 * D is low ON a rung, and low-height rungs are cheaper: the program's
 * own comb, ported from frequency ratios to size ratios.
 *
 * Drive — each cell claims space in proportion to its energy:  -g_i x_i.
 * NOTHING in U contains a maximum size. If a bound appears it is made
 * of intervals, not of lengths.
 *
 *   U = sum_{i<j} w_ij D(u_ij)  -  sum_i g_i x_i     subject to sum x = 0
 *
 * WHAT WOULD COUNT AS THE IDEA WORKING
 * ------------------------------------
 * If the bound is harmonic then r(g) must be a STAIRCASE: sweep one
 * cell's energy continuously and its size should sit on plateaux (locked
 * to a rung against its neighbours) and jump between them. Smooth r(g)
 * means harmony bounds nothing and idea (1) needs an imposed constant.
 *
 * The negative control is the same sweep with the comb collapsed to the
 * unison. A staircase THERE would mean the steps come from the
 * constraint or the descent rather than from the comb.
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ------------------------------------------------------------ rungs */

typedef struct { int p, q; double u, H; } Rung;

static Rung RUNGS[8192];
static int NRUNG = 0;

static int igcd(int a, int b) { while (b) { int t = a % b; a = b; b = t; } return a; }

/* every reduced p/q with Tenney height log2(p*q) <= climb, both
 * orientations. climb < 0 collapses the comb to the unison (control). */
static void build_rungs(double climb)
{
    NRUNG = 0;
    int lim = (int)(pow(2.0, climb > 0 ? climb : 0.0) + 0.5);
    if (lim < 1) lim = 1;
    for (int p = 1; p <= lim; p++)
        for (int q = 1; q <= lim; q++) {
            if (p * q > lim) continue;
            if (igcd(p, q) != 1) continue;
            double H = log2((double)p * (double)q);
            if (H > climb + 1e-12) continue;
            if (NRUNG >= 8192) continue;
            RUNGS[NRUNG].p = p; RUNGS[NRUNG].q = q;
            RUNGS[NRUNG].u = log2((double)p / (double)q);
            RUNGS[NRUNG].H = H;
            NRUNG++;
        }
}

static double diss(double u, double sigma, double lambda, int *kbest)
{
    double best = 1e300; int kb = 0;
    for (int k = 0; k < NRUNG; k++) {
        double d = u - RUNGS[k].u;
        double v = d * d / (2 * sigma * sigma) + lambda * RUNGS[k].H;
        if (v < best) { best = v; kb = k; }
    }
    if (kbest) *kbest = kb;
    return best;
}

static double diss_grad(double u, double sigma, double lambda)
{
    int k; diss(u, sigma, lambda, &k);
    return (u - RUNGS[k].u) / (sigma * sigma);
}

/* ------------------------------------------------------- the system */

typedef struct {
    int N;
    double sigma, lambda;
    double kE;           /* stiffness of a cell's own preferred size */
    double *x, *g, *xt;  /* x = ln(extent); g unused; xt = preferred ln size */
    double *cw;          /* coupling weights, N*N */
} Sys;

static Sys *sys_new(int N)
{
    Sys *s = calloc(1, sizeof(Sys));
    s->N = N; s->sigma = 0.05; s->lambda = 0.10; s->kE = 1.0;
    s->x = calloc(N, sizeof(double));
    s->xt = calloc(N, sizeof(double));
    s->g = malloc(N * sizeof(double));
    s->cw = calloc((size_t)N * N, sizeof(double));
    for (int i = 0; i < N; i++) s->g[i] = 1.0;
    return s;
}
static void sys_free(Sys *s)
{ free(s->x); free(s->xt); free(s->g); free(s->cw); free(s); }

static void couple_ring(Sys *s)
{
    int N = s->N;
    memset(s->cw, 0, (size_t)N * N * sizeof(double));
    for (int i = 0; i < N; i++) {
        int j = (i + 1) % N;
        s->cw[i * N + j] = s->cw[j * N + i] = 1.0;
    }
}

static void couple_decay(Sys *s, double alpha)
{
    int N = s->N;
    memset(s->cw, 0, (size_t)N * N * sizeof(double));
    for (int i = 0; i < N; i++)
        for (int j = i + 1; j < N; j++) {
            int d = j - i; if (d > N - d) d = N - d;
            double w = exp(-alpha * (d - 1));
            s->cw[i * N + j] = s->cw[j * N + i] = w;
        }
}

static double uij(Sys *s, int i, int j) { return (s->x[j] - s->x[i]) / M_LN2; }

static double energy(Sys *s)
{
    int N = s->N; double U = 0;
    for (int i = 0; i < N; i++) {
        for (int j = i + 1; j < N; j++) {
            double w = s->cw[i * N + j];
            if (w == 0) continue;
            U += w * diss(uij(s, i, j), s->sigma, s->lambda, NULL);
        }
        double d = s->x[i] - s->xt[i];
        U += 0.5 * s->kE * d * d;
    }
    return U;
}

static void grad(Sys *s, double *G)
{
    int N = s->N;
    memset(G, 0, N * sizeof(double));
    for (int i = 0; i < N; i++) {
        for (int j = i + 1; j < N; j++) {
            double w = s->cw[i * N + j];
            if (w == 0) continue;
            double dd = w * diss_grad(uij(s, i, j), s->sigma, s->lambda);
            G[i] += dd * (-1.0 / M_LN2);
            G[j] += dd * (+1.0 / M_LN2);
        }
        G[i] += s->kE * (s->x[i] - s->xt[i]);
    }
    double m = 0;                 /* project onto sum x = 0 */
    for (int i = 0; i < N; i++) m += G[i];
    m /= N;
    for (int i = 0; i < N; i++) G[i] -= m;
}

static double relax(Sys *s, int iters, double lr, double *gmax)
{
    /* Adaptive-step descent. The comb curvature is 1/sigma^2 per pair
     * (1/0.0025 = 400 at sigma=0.05), so the stable fixed step is
     * ~1e-3; a fixed lr chosen by eye diverges silently, U goes NaN,
     * and — before the guard below — an all-NaN restart set left the
     * caller reading uninitialised memory. Backtracking removes the
     * guesswork; the NaN guard makes a failure loud instead of subtle. */
    int N = s->N;
    double *G = malloc(N * sizeof(double));
    double *xb = malloc(N * sizeof(double));
    double U = energy(s), step = lr;
    for (int it = 0; it < iters; it++) {
        grad(s, G);
        memcpy(xb, s->x, N * sizeof(double));
        double Un = 0;
        int ok = 0;
        for (int tries = 0; tries < 40; tries++) {
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
        U = Un;
        step *= 1.1;
        if (step > 1.0) step = 1.0;
    }
    grad(s, G);
    *gmax = 0;
    for (int i = 0; i < N; i++)
        if (fabs(G[i]) > *gmax) *gmax = fabs(G[i]);
    U = energy(s);
    free(G); free(xb);
    return U;
}

/* deterministic xorshift — runs reproduce byte-for-byte */
static unsigned long long RS = 88172645463325252ULL;
static void rseed(unsigned long long s) { RS = s ? s : 1; }
static double rnd(void)
{
    RS ^= RS << 13; RS ^= RS >> 7; RS ^= RS << 17;
    return (double)((RS >> 11) & 0xFFFFFFFFFFFFFULL) / 9007199254740992.0;
}

/* The landscape is a min-over-rungs and therefore rugged: one seed
 * selects one basin. That is exactly how the FK saddle got mistaken for
 * the minimum in the charge work. Always restart. */
static double relax_best(Sys *s, int nrest, int iters, double lr,
                         double spread, double *gmax)
{
    int N = s->N;
    double *best = malloc(N * sizeof(double));
    double Ub = 1e300, gb = 0;
    int have = 0;              /* never let an all-failed restart set
                                * leak uninitialised memory into results */
    for (int t = 0; t < nrest; t++) {
        for (int i = 0; i < N; i++)
            s->x[i] = (t == 0) ? 0.0 : spread * (rnd() - 0.5);
        double m = 0;
        for (int i = 0; i < N; i++) m += s->x[i];
        m /= N;
        for (int i = 0; i < N; i++) s->x[i] -= m;
        double gm, U = relax(s, iters, lr, &gm);
        if (isfinite(U) && (!have || U < Ub)) {
            Ub = U; gb = gm; have = 1;
            memcpy(best, s->x, N * sizeof(double));
        }
    }
    if (!have) { fprintf(stderr, "FATAL: every restart failed\n"); exit(2); }
    memcpy(s->x, best, N * sizeof(double));
    free(best);
    *gmax = gb;
    return Ub;
}

/* ------------------------------------------------------- diagnostics */

static int nsizes(Sys *s, double reltol)
{
    int N = s->N, n = 0; double v[4096];
    for (int i = 0; i < N; i++) {
        int f = 0;
        for (int k = 0; k < n; k++)
            if (fabs(s->x[i] - v[k]) < reltol) { f = 1; break; }
        if (!f && n < 4096) v[n++] = s->x[i];
    }
    return n;
}

/* ------------------------------------------------------------ main */

/* the ratio cell0 settles at, against its ring neighbour, at drive g0 */
static double probe2(int N, double g0, double sigma, double climb,
                     double kE, int *pk, double *gm)
{
    /* RUNG-SEEDED SEARCH. The landscape is a min-over-rungs: with 211
     * rungs at comb_limit=6, random restarts find a different basin at
     * adjacent drives and the sweep fills with jumps that are search
     * noise, not physics. Seeding one relaxation AT each candidate rung
     * makes the search exhaustive over exactly the basins that exist. */
    build_rungs(climb);
    Sys *s = sys_new(N);
    couple_ring(s);
    s->sigma = sigma;
    s->kE = kE;
    s->xt[0] = g0 * M_LN2;   /* g0 is now the PREFERRED log2-size of cell 0 */
    double *best = malloc(N * sizeof(double));
    double Ub = 0, gb = 0;
    int have = 0;
    for (int k = 0; k < NRUNG; k++) {
        double u0 = RUNGS[k].u * M_LN2;       /* put cell0 on rung k */
        for (int i = 0; i < N; i++) s->x[i] = (i == 0) ? u0 : 0.0;
        double m = 0;
        for (int i = 0; i < N; i++) m += s->x[i];
        m /= N;
        for (int i = 0; i < N; i++) s->x[i] -= m;
        double gm2, U = relax(s, 3000, 0.002, &gm2);
        if (isfinite(U) && (!have || U < Ub)) {
            Ub = U; gb = gm2; have = 1;
            memcpy(best, s->x, N * sizeof(double));
        }
    }
    if (!have) { fprintf(stderr, "FATAL: no basin converged\n"); exit(2); }
    memcpy(s->x, best, N * sizeof(double));
    free(best);
    *gm = gb;
    double u = uij(s, 1, 0);
    diss(u, s->sigma, s->lambda, pk);
    sys_free(s);
    return u;
}

static double probe(int N, double g0, double sigma, double climb,
                    int *pk, double *gm)
{ return probe2(N, g0, sigma, climb, 1.0, pk, gm); }

int main(void)
{
    printf("# harmgrow — is the morphology bound harmonic?\n");
    printf("# x=ln r, sum x = 0; pitch w=1/r; pair cost = comb dissonance\n");
    printf("# of log2(r_j/r_i); drive -g_i x_i. NO maximum size in U.\n\n");

    printf("== 1. the rung set (what 'consonant' means here) ==\n");
    for (double cl = 0; cl <= 6.0; cl += 2.0) {
        build_rungs(cl);
        double hi = -1e9;
        for (int k = 0; k < NRUNG; k++) if (RUNGS[k].u > hi) hi = RUNGS[k].u;
        printf("comb_limit=%.0f  nrungs=%4d  u span=[%+.2f,%+.2f] "
               "=> ratio span up to %.0fx\n", cl, NRUNG, -hi, hi, pow(2, hi));
    }
    printf("\nNOTE: the comb has a FINITE SPAN. Past the outermost rung the\n");
    printf("cost grows quadratically with no rung left to relax to. That is\n");
    printf("the only thing here that could bound a size, and it is made of\n");
    printf("intervals (dimensionless), not of lengths.\n\n");

    /* ---------------- 2. THE STAIRCASE ---------------- */
    printf("== 2. STAIRCASE TEST — sweep one cell's PREFERRED size ==\n");
    printf("ring N=12, neighbours only, comb_limit=6, sigma=0.05, kE=1\n");
    printf("Cell 0's own energy prefers log2-size t (the x axis). The comb\n");
    printf("prefers a rung against its neighbours. If harmony quantises,\n");
    printf("the ACHIEVED size lags the preferred one in steps.\n");
    printf("  t_pref  r0/r_nb    log2     rung    |u-u_k|   max|grad|\n");
    {
        int pk = 0; double gm = 0, prevu = -1e9;
        for (double g0 = 0; g0 <= 3.001; g0 += 0.02) {
            double u = probe(12, g0, 0.05, 6.0, &pk, &gm);
            const char *mark = (fabs(u - prevu) > 1e-3) ? " <-- JUMP" : "";
            if (prevu < -1e8) mark = "";
            printf("%7.3f  %8.4f  %+7.4f  %3d/%-3d  %8.2e  %8.1e%s\n",
                   g0, pow(2, u), u, RUNGS[pk].p, RUNGS[pk].q,
                   fabs(u - RUNGS[pk].u), gm, mark);
            prevu = u;
        }
    }
    printf("\nreading: flat runs = the size is LOCKED to a rung and further\n");
    printf("drive is absorbed by the lock; jumps = the lock breaking to the\n");
    printf("next rung. Smooth = no harmonic quantisation.\n\n");

    /* ---------------- 3. NEGATIVE CONTROL ---------------- */
    printf("== 3. NEGATIVE CONTROL — same sweep, comb collapsed to unison ==\n");
    printf("  t_pref  r0/r_nb    log2\n");
    {
        int pk = 0; double gm = 0;
        for (double g0 = 0; g0 <= 3.001; g0 += 0.15) {
            double u = probe(12, g0, 0.05, -1.0, &pk, &gm);
            printf("%7.3f  %8.4f  %+7.4f\n", g0, pow(2, u), u);
        }
    }
    printf("\nreading: MUST be smooth. A staircase here would falsify §2.\n\n");

    /* ---------------- 4. LOCK WIDTH ---------------- */
    printf("== 4. HOW SHARP IS THE QUANTISATION? vs lock width sigma ==\n");
    printf("sweep g0 = 4..80 at each sigma; count distinct rungs visited and\n");
    printf("how much of the sweep is spent locked (|u-u_k| < sigma/10).\n");
    printf(" sigma  rungs_visited  jumps  frac_locked  max_ratio\n");
    {
        double sig[] = {0.02, 0.05, 0.10, 0.15, 0.20, 0.30, 0.50};
        for (unsigned t = 0; t < sizeof(sig) / sizeof(sig[0]); t++) {
            int pk = 0, seen[512], ns = 0, jumps = 0, nlock = 0, n = 0;
            double gm = 0, prevu = -1e9, mx = 0;
            for (double g0 = 0; g0 <= 3.001; g0 += 0.05) {
                double u = probe(12, g0, sig[t], 6.0, &pk, &gm);
                int f = 0;
                for (int k = 0; k < ns; k++) if (seen[k] == pk) { f = 1; break; }
                if (!f && ns < 512) seen[ns++] = pk;
                if (prevu > -1e8 && fabs(u - prevu) > 1e-3) jumps++;
                if (fabs(u - RUNGS[pk].u) < sig[t] / 10) nlock++;
                if (pow(2, u) > mx) mx = pow(2, u);
                prevu = u; n++;
            }
            printf("%6.2f  %13d  %5d  %11.3f  %9.2f\n",
                   sig[t], ns, jumps, (double)nlock / n, mx);
        }
    }
    printf("\nreading: a stiff comb (small sigma) locks hard but reaches few\n");
    printf("rungs; a loose comb reaches many but holds none. Quantisation\n");
    printf("needs a window in sigma, not a limit.\n\n");

    /* ---------------- 4b. THE COMPETITION RATIO ---------------- */
    printf("== 4b. WHEN DOES QUANTISATION EXIST? comb vs the cell's own\n");
    printf("energetic preference. Comb curvature at sigma=0.05 is 400.\n");
    printf("     kE  rungs_visited  jumps/151  max_log2_achieved (pref 3.000)\n");
    {
        double kes[] = {1, 10, 50, 200, 400, 2000, 8000};
        for (unsigned t = 0; t < sizeof(kes) / sizeof(kes[0]); t++) {
            int pk, ns = 0, jumps = 0, seen[512];
            double gm, prevu = -1e9, mx = 0;
            for (double tp = 0; tp <= 3.001; tp += 0.02) {
                double u = probe2(12, tp, 0.05, 6.0, kes[t], &pk, &gm);
                int f = 0;
                for (int k = 0; k < ns; k++) if (seen[k] == pk) { f = 1; break; }
                if (!f && ns < 512) seen[ns++] = pk;
                if (prevu > -1e8 && fabs(u - prevu) > 1e-3) jumps++;
                if (u > mx) mx = u;
                prevu = u;
            }
            printf("%7.0f  %13d  %9d  %.3f\n", kes[t], ns, jumps, mx);
        }
    }
    printf("\nreading: few rungs + few jumps = long plateaux = QUANTISED.\n");
    printf("jumps ~ every sample = the elastic term wins and size follows\n");
    printf("preference continuously. Quantisation is a WINDOW in kE, and\n");
    printf("max_log2 < 3.000 means the comb held the cell BELOW the size\n");
    printf("its own energy wanted — a bound, made of intervals.\n\n");

    /* ---------------- 5. comb_limit window ---------------- */
    printf("== 5. IS THERE A WINDOW IN comb_limit? ==\n");
    printf("(too few rungs -> nothing to lock to; too many -> rungs dense\n");
    printf(" everywhere -> no selection). g0 = 48, sigma = 0.15.\n");
    printf("comb_limit  nrungs   ratio    rung    |u-u_k|\n");
    {
        for (double cl = 1; cl <= 8.001; cl += 1.0) {
            int pk = 0; double gm = 0;
            double u = probe(12, 1.7, 0.15, cl, &pk, &gm);
            printf("%8.0f  %6d  %8.3f  %3d/%-3d  %9.2e\n",
                   cl, NRUNG, pow(2, u), RUNGS[pk].p, RUNGS[pk].q,
                   fabs(u - RUNGS[pk].u));
        }
    }
    printf("\n");

    /* ---------------- 6. within vs between ---------------- */
    printf("== 6. HARMONICS WITHIN, DISSONANCE BETWEEN ==\n");
    printf("ring of 16, all-pairs weight exp(-alpha*(d-1)), sigma=0.15,\n");
    printf("drives ramped 1..64 so the comb is actually engaged.\n");
    printf(" alpha  n_sizes  medoff_near  medoff_far  H_near  H_far\n");
    {
        build_rungs(6.0);
        for (double alpha = 0.25; alpha <= 2.001; alpha += 0.25) {
            rseed(20260802);
            Sys *s = sys_new(16);
            couple_decay(s, alpha);
            s->sigma = 0.15;
            for (int i = 0; i < 16; i++) s->xt[i] = (i * 4.0 / 15.0) * M_LN2;
            double gm;
            relax_best(s, 8, 3000, 0.002, 2.0, &gm);
            int N = 16, nn = 0, nf = 0;
            double on[256], of[256], Hn = 0, Hf = 0;
            for (int i = 0; i < N; i++)
                for (int j = i + 1; j < N; j++) {
                    int d = j - i; if (d > N - d) d = N - d;
                    int k; double u = uij(s, i, j);
                    diss(u, s->sigma, s->lambda, &k);
                    double off = fabs(u - RUNGS[k].u);
                    if (d <= 2 && nn < 256) { on[nn++] = off; Hn += RUNGS[k].H; }
                    else if (d >= 4 && nf < 256) { of[nf++] = off; Hf += RUNGS[k].H; }
                }
            /* medians */
            for (int a = 0; a < nn; a++) for (int b = a + 1; b < nn; b++)
                if (on[b] < on[a]) { double t2 = on[a]; on[a] = on[b]; on[b] = t2; }
            for (int a = 0; a < nf; a++) for (int b = a + 1; b < nf; b++)
                if (of[b] < of[a]) { double t2 = of[a]; of[a] = of[b]; of[b] = t2; }
            printf("%6.2f  %7d  %11.4f  %10.4f  %6.2f  %6.2f\n",
                   alpha, nsizes(s, 1e-2),
                   nn ? on[nn / 2] : -1, nf ? of[nf / 2] : -1,
                   nn ? Hn / nn : -1, nf ? Hf / nf : -1);
            sys_free(s);
        }
    }
    printf("\nreading: medoff_near << medoff_far, with H_near < H_far, IS\n");
    printf("'harmonics within, dissonance between' as a measurement.\n\n");

    /* ---------------- 7. the bound ---------------- */
    printf("== 7. THE BOUND — push one cell arbitrarily hard ==\n");
    printf("   g0      r0/r_nb    log2     rung    max|grad|\n");
    {
        double gs[] = {0.5, 1, 2, 3, 4, 5, 6, 7, 8, 10};
        for (unsigned t = 0; t < sizeof(gs) / sizeof(gs[0]); t++) {
            int pk = 0; double gm = 0;
            double u = probe(12, gs[t], 0.05, 6.0, &pk, &gm);
            printf("%6.1f  %10.3f  %+7.3f  %3d/%-3d  %8.1e\n",
                   gs[t], pow(2, u), u, RUNGS[pk].p, RUNGS[pk].q, gm);
        }
    }
    printf("\nreading: the outermost rung is 64/1 at u=6. Past it the wall is\n");
    printf("quadratic and the drive is linear, so equilibrium sits at\n");
    printf("u = u_max + g*sigma^2*ln2/2 — finite for any g, but UNBOUNDED as\n");
    printf("g grows. The comb QUANTISES size; it does not CAP it.\n");
    return 0;
}
