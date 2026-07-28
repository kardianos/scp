/*  bell_grid.c — v87 B1: exact-integer, streaming, fine-grid CHSH for the
 *  fabric phase construction.
 *
 *  Three things make the large run possible:
 *
 *  1. STREAMING SAMPLES. lambda is never stored. Each thread generates and
 *     histograms on the fly, so N is limited by time, not by memory. (The
 *     previous version malloc'd N doubles; at N = 1e11 that is 800 GB.)
 *
 *  2. EXACT INTEGER ARITHMETIC. Angles are uint32 "turns" (2^32 = one
 *     revolution) so modular angle arithmetic is exact -- no fmod, no pi. The
 *     response square waves are bit-range tests, not cos() calls:
 *         sgn cos(u) = +1  iff  turn in [0,1/4) U [3/4,1)
 *     Histogram bins and the correlation sum are __int128. For arm 1 the
 *     result E = k/N is an exact rational; for arm 2 the weights are Q30
 *     fixed point.
 *
 *  3. O(grid^2) SEARCH, not O(grid^3). With a = 0,
 *         S = [E(-b) + E(a'-b)] + [E(-b') - E(a'-b')] = f_{a'}(b) + g_{a'}(b')
 *     which is SEPARABLE in b and b'. For each a' the two maxima (and minima,
 *     for |S|) are found independently in O(grid). A grid of 65536 settings --
 *     0.0055 degree resolution -- costs 4.3e9 table reads instead of 2.8e14.
 *
 *  Also computes the ANALYTIC table (no sampling at all) as a control: the
 *  search over an exact E(theta) must return exactly 2 (arm 1) and 2sqrt2
 *  (arm 2). Any excess in the sampled run is then unambiguously statistical.
 *
 *  Build: gcc -O3 -march=native -fopenmp -o bell_grid bell_grid.c -lm
 *  Usage: bell_grid [N] [grid] [nb] [nth]
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <omp.h>

#define TAU 6.283185307179586476925287

static inline uint64_t rotl(const uint64_t x, int k) { return (x << k) | (x >> (64 - k)); }
typedef struct { uint64_t s[4]; } rng_t;
static inline uint64_t rnext(rng_t *r) {
    const uint64_t res = rotl(r->s[0] + r->s[3], 23) + r->s[0];
    const uint64_t t = r->s[1] << 17;
    r->s[2] ^= r->s[0]; r->s[3] ^= r->s[1]; r->s[1] ^= r->s[2]; r->s[0] ^= r->s[3];
    r->s[2] ^= t; r->s[3] = rotl(r->s[3], 45);
    return res;
}
static void rseed(rng_t *r, uint64_t seed) {
    for (int i = 0; i < 4; i++) {
        seed += 0x9E3779B97F4A7C15ULL;
        uint64_t z = seed;
        z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
        z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
        r->s[i] = z ^ (z >> 31);
    }
}
static double now(void) { return omp_get_wtime(); }

/* +1 on the first and last quarter turn -- sgn cos, with no transcendental */
static inline int sq(uint32_t turn) { return (turn < 0x40000000u || turn >= 0xC0000000u) ? 1 : -1; }

/* ------------------------------------------------------- streaming histogram */
static void build_table(long double *E, long long N, int nb, int nth, int arm,
                        int logb, __int128 *out_den, uint64_t seed,
                        double *out_meanw) {
    __int128 *h = calloc(nb, sizeof(__int128));
    const int SH = 30;
    long long *wq = malloc(sizeof(long long) * nb);
    for (int j = 0; j < nb; j++) {
        double u = TAU * (j + 0.5) / nb;
        wq[j] = (arm == 2) ? (long long)llround(fabs(cos(u)) * (1LL << SH)) : 1LL;
    }
    #pragma omp parallel
    {
        int tid = omp_get_thread_num(), nt = omp_get_num_threads();
        rng_t r; rseed(&r, seed + 7919ULL * tid);
        __int128 *hl = calloc(nb, sizeof(__int128));
        long long lo = (long long)((__int128)N * tid / nt);
        long long hi = (long long)((__int128)N * (tid + 1) / nt);
        for (long long i = lo; i < hi; i++) {
            uint32_t turn = (uint32_t)(rnext(&r) >> 32);   /* uniform on the circle */
            hl[turn >> (32 - logb)] += wq[turn >> (32 - logb)];
        }
        #pragma omp critical
        for (int j = 0; j < nb; j++) h[j] += hl[j];
        free(hl);
    }
    __int128 den = 0;
    for (int j = 0; j < nb; j++) den += h[j];
    *out_den = den;
    /* mean sampling weight = detector firing probability under Reading 2.
     * arm 1: w == 1 identically, so the mean is 1 and nothing is discarded.
     * arm 2: w = |cos u| in Q30, so <w> = den / (N * 2^SH). */
    if (out_meanw)
        *out_meanw = (arm == 2)
            ? (double)((long double)den / ((long double)N * (long double)(1LL << SH)))
            : 1.0;
    signed char *As = malloc(nb);
    for (int j = 0; j < nb; j++) As[j] = (j < nb / 4 || j >= 3 * nb / 4) ? 1 : -1;
    #pragma omp parallel for schedule(static)
    for (int k = 0; k < nth; k++) {
        int shift = (int)((long long)k * nb / nth);
        __int128 num = 0;
        for (int j = 0; j < nb; j++) {
            int jj = j + shift; if (jj >= nb) jj -= nb;
            num += h[j] * (__int128)(-(int)As[jj] * (int)As[j]);
        }
        E[k] = (long double)num / (long double)den;
    }
    free(h); free(As); free(wq);
}

/* ------------------------------------------------------------ analytic table */
static void build_table_exact(long double *E, int nth, int arm) {
    for (int k = 0; k < nth; k++) {
        long double th = (long double)TAU * k / nth;
        long double t = th > (long double)M_PI ? (long double)TAU - th : th;
        E[k] = (arm == 1) ? (-1.0L + 2.0L * t / (long double)M_PI) : -cosl(th);
    }
}

/* -------------------------------------------------- O(grid^2) separable search
 * Returns max|S| and, in `st`, the settings that attained it. Keeping the
 * argmax is what makes the search-bias null possible: the settings found on
 * one table can be replayed, WITHOUT re-searching, on an independent one. */
typedef struct { int iap, ib, ibp; } settings_t;

static long double chsh_sep_arg(const long double *E, int nth, int grid,
                                settings_t *st) {
    settings_t *cand = malloc(sizeof(settings_t) * grid);
    long double *cval = malloc(sizeof(long double) * grid);
    #pragma omp parallel for schedule(static)
    for (int ia = 0; ia < grid; ia++) {
        long long apk = (long long)ia * nth / grid;
        long double fmax = -1e30L, fmin = 1e30L, gmax = -1e30L, gmin = 1e30L;
        int bfx = 0, bfn = 0, bgx = 0, bgn = 0;
        for (int ib = 0; ib < grid; ib++) {
            long long bk = (long long)ib * nth / grid;
            long long i1 = ((-bk) % nth + nth) % nth;
            long long i2 = ((apk - bk) % nth + nth) % nth;
            long double f = E[i1] + E[i2];
            long double g = E[i1] - E[i2];
            if (f > fmax) { fmax = f; bfx = ib; }
            if (f < fmin) { fmin = f; bfn = ib; }
            if (g > gmax) { gmax = g; bgx = ib; }
            if (g < gmin) { gmin = g; bgn = ib; }
        }
        long double s1 = fabsl(fmax + gmax), s2 = fabsl(fmin + gmin);
        if (s1 >= s2) { cval[ia] = s1; cand[ia] = (settings_t){ia, bfx, bgx}; }
        else          { cval[ia] = s2; cand[ia] = (settings_t){ia, bfn, bgn}; }
    }
    long double best = -1.0L; settings_t bs = {0, 0, 0};
    for (int ia = 0; ia < grid; ia++)
        if (cval[ia] > best) { best = cval[ia]; bs = cand[ia]; }
    free(cand); free(cval);
    if (st) *st = bs;
    return best;
}

static long double chsh_sep(const long double *E, int nth, int grid) {
    return chsh_sep_arg(E, nth, grid, NULL);
}

/* |S| at FIXED, pre-registered settings — no maximisation, hence unbiased. */
static long double chsh_at(const long double *E, int nth, int grid,
                           settings_t st) {
    long long apk = (long long)st.iap * nth / grid;
    long long bk  = (long long)st.ib  * nth / grid;
    long long bpk = (long long)st.ibp * nth / grid;
    #define IDX(x) (long double)E[(int)(((x) % nth + nth) % nth)]
    long double S = (IDX(-bk) + IDX(apk - bk)) + (IDX(-bpk) - IDX(apk - bpk));
    #undef IDX
    return fabsl(S);
}

/* --------------------------------------------------------- discard fraction
 * Reading 2 of rho_p (CRANK2 §4): the source is uniform and the DETECTOR fires
 * with probability proportional to the weight. Then <w> is a coincidence
 * efficiency and 1-<w> is thrown away. Report it against the standard
 * thresholds so a "violation" bought with discards is labelled as what it is. */
static void report_discard(int arm, double meanw, long double S) {
    const double thr_chsh = 2.0 * (sqrt(2.0) - 1.0);   /* 0.828427 */
    const double thr_eb   = 2.0 / 3.0;                 /* Eberhard/CH, ~0.6667 */
    printf("  DISCARD  <w> = eta = %.6f   discarded = %.6f\n",
           meanw, 1.0 - meanw);
    if (arm == 1) {
        printf("            arm 1 is unweighted: nothing is discarded, so no\n"
               "            detection loophole is available to it. |S| = %.6Lf\n", S);
        return;
    }
    printf("            CHSH detection threshold 2(sqrt2-1) = %.6f  -> %s\n",
           thr_chsh, meanw >= thr_chsh ? "ABOVE (loophole closed)"
                                       : "BELOW (loophole OPEN)");
    printf("            Eberhard/CH threshold      2/3      = %.6f  -> %s\n",
           thr_eb, meanw >= thr_eb ? "ABOVE" : "BELOW");
    printf("            VERDICT: read detector-side, this arm's |S| = %.6Lf is\n"
           "            obtained at eta = %.4f, i.e. inside the detection\n"
           "            loophole. It is NOT a loophole-free violation.\n",
           S, meanw);
}

int main(int argc, char **argv) {
    long long N = (argc > 1) ? atoll(argv[1]) : 100000000000LL;
    int grid = (argc > 2) ? atoi(argv[2]) : 65536;
    int nb   = (argc > 3) ? atoi(argv[3]) : 262144;
    int nth  = (argc > 4) ? atoi(argv[4]) : 262144;
    int reps = 0; long long Nnull = 0;
    for (int i = 1; i < argc; i++)
        if (!strcmp(argv[i], "--null") && i + 2 < argc) {
            reps = atoi(argv[i + 1]); Nnull = atoll(argv[i + 2]);
        }
    int logb = (int)llround(log2((double)nb));
    if ((1 << logb) != nb) { fprintf(stderr, "nb must be a power of two\n"); return 1; }

    printf("bell_grid   N=%lld  grid=%d (%.5f deg)  nb=%d  nth=%d  threads=%d\n",
           N, grid, 360.0 / grid, nb, nth, omp_get_max_threads());
    printf("  exact targets: arm1 = 2 (triangle) ; arm2 = 2sqrt2 = %.12f\n\n",
           2.0 * sqrt(2.0));

    long double *E = malloc(sizeof(long double) * nth);

    for (int arm = 1; arm <= 2; arm++) {
        long double tgt = (arm == 1) ? 2.0L : 2.0L * sqrtl(2.0L);
        printf("--- ARM %d ---\n", arm);

        /* ---- REQUIREMENT 1: analytic control. The search is run on an exact
         * E(theta) whose max|S| is known in closed form. Any deviation here is
         * grid discretisation, not statistics, and it bounds how much of a
         * sampled excess the search itself can manufacture from a clean table. */
        double t0 = now();
        build_table_exact(E, nth, arm);
        settings_t amx;
        long double Sx = chsh_sep_arg(E, nth, grid, &amx);
        /* Pre-registration must not be read off the analytic argmax: for the
         * triangle the maximum sits on a plateau, and the search happens to
         * return the DEGENERATE point a'=b'=0, where S = 2E(pi) with the two
         * noisy terms cancelling identically. That estimator has no variance
         * and tests nothing. Pin the canonical CHSH quadruple instead --
         * (a,a',b,b') = (0, 90, 225, 135) degrees -- which attains the maximum
         * for BOTH arms (2 for the triangle, 2sqrt2 for the cosine) with four
         * distinct, independently noisy table entries. */
        settings_t pre = { grid / 4, 5 * grid / 8, 3 * grid / 8 };
        printf("  ANALYTIC (no sampling)  max|S| = %.12Lf   excess %+.3e   %6.2f s\n",
               Sx, (double)(Sx - tgt), now() - t0);
        printf("            analytic argmax:  a'=%.4f  b=%.4f  b'=%.4f deg\n",
               360.0 * amx.iap / grid, 360.0 * amx.ib / grid, 360.0 * amx.ibp / grid);
        printf("            pre-registered:   a'=%.4f  b=%.4f  b'=%.4f deg"
               "   -> analytic |S| = %.12Lf\n",
               360.0 * pre.iap / grid, 360.0 * pre.ib / grid, 360.0 * pre.ibp / grid,
               chsh_at(E, nth, grid, pre));

        t0 = now();
        __int128 den = 0; double meanw = 1.0;
        build_table(E, N, nb, nth, arm, logb, &den, 20260726ULL, &meanw);
        double tb = now() - t0;
        t0 = now();
        settings_t got;
        long double Ss = chsh_sep_arg(E, nth, grid, &got);
        double ts = now() - t0;
        long double Sfix = chsh_at(E, nth, grid, pre);
        printf("  SAMPLED  (int128 exact) max|S| = %.12Lf   excess %+.3e   "
               "%6.2f s build + %5.2f s search\n", Ss, (double)(Ss - tgt), tb, ts);
        printf("            1/sqrt(N) = %.3e   ->  excess/noise = %.2f\n",
               1.0 / sqrt((double)N), (double)(Ss - tgt) * sqrt((double)N));
        printf("  FIXED    (pre-registered, no maximisation, UNBIASED)\n");
        printf("            |S| = %.12Lf   excess %+.3e  -> %.2f sigma\n",
               Sfix, (double)(Sfix - tgt), (double)(Sfix - tgt) * sqrt((double)N));
        printf("            search bias = max|S| - |S|_fixed = %+.3e\n",
               (double)(Ss - Sfix));
        printf("            E(0)=%.12Lf E(60)=%.12Lf E(90)=%.12Lf E(180)=%.12Lf\n",
               E[0], E[nth / 6], E[nth / 4], E[nth / 2]);

        /* ---- REQUIREMENT 3: discard fraction. */
        report_discard(arm, meanw, Ss);
        printf("\n");
    }

    /* ---- REQUIREMENT 2: search-bias null.
     * The searched max|S| is a maximum over grid^2 noisy settings, so it is
     * biased UP even when the true S is exactly the target. Measure that bias
     * directly: replicate the null ensemble at N_null with independent seeds,
     * and report the excess in units of 1/sqrt(N). If c = excess*sqrt(N) is
     * stable in N, the null band at production N is c/sqrt(N), and no claimed
     * violation below that band counts. */
    if (reps > 0) {
        printf("=== SEARCH-BIAS NULL  (%d replicas at N=%lld) ===\n", reps, Nnull);
        printf("    the true |S| is exactly the target by construction; every\n"
               "    excess below is manufactured by maximising over the grid.\n\n");
        for (int arm = 1; arm <= 2; arm++) {
            long double tgt = (arm == 1) ? 2.0L : 2.0L * sqrtl(2.0L);
            settings_t pre = { grid / 4, 5 * grid / 8, 3 * grid / 8 };
            double cs = 0, cs2 = 0, cmax = -1e30, cmin = 1e30;
            double fs = 0, fs2 = 0;
            for (int r = 0; r < reps; r++) {
                __int128 den = 0; double mw = 1.0;
                build_table(E, Nnull, nb, nth, arm, logb, &den,
                            1000003ULL * (r + 1) + 17ULL * arm, &mw);
                long double Sr = chsh_sep(E, nth, grid);
                long double Sf = chsh_at(E, nth, grid, pre);
                double c = (double)(Sr - tgt) * sqrt((double)Nnull);
                double f = (double)(Sf - tgt) * sqrt((double)Nnull);
                cs += c; cs2 += c * c; fs += f; fs2 += f * f;
                if (c > cmax) cmax = c;  if (c < cmin) cmin = c;
                printf("    arm %d rep %2d:  searched c = %+7.3f   fixed c = %+7.3f\n",
                       arm, r, c, f);
            }
            double cm = cs / reps, csd = sqrt(fmax(0.0, cs2 / reps - cm * cm));
            double fm = fs / reps, fsd = sqrt(fmax(0.0, fs2 / reps - fm * fm));
            printf("    arm %d SEARCHED  c = %+.3f +/- %.3f   [%.3f, %.3f]\n",
                   arm, cm, csd, cmin, cmax);
            printf("    arm %d FIXED     c = %+.3f +/- %.3f   (should straddle 0)\n",
                   arm, fm, fsd);
            printf("    arm %d null band at production N=%lld: excess < %.3e\n\n",
                   arm, N, (cm + 3.0 * csd) / sqrt((double)N));
        }
    }
    free(E);
    return 0;
}
