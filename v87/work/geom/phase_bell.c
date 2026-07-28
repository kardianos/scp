/*  phase_bell.c — v87 B1: in-kernel Bell construction from the fabric's phase
 *  structure. C reference implementation, and a race against phase_bell.py.
 *
 *  PROTOCOL (identical physics to phase_bell.py):
 *    object 1 phase  lambda + omega t          -> wing A, detector offset a
 *    object 2 phase  lambda + pi + omega t     -> wing B, detector offset b
 *    detectors carry the SAME omega, so omega*t cancels in both relative phases
 *    readout (A6 fuse/repel dichotomy):  A = sgn cos(lambda - a)
 *                                        B = sgn cos(lambda + pi - b)
 *    arm 1: lambda ~ Uniform                       -- (MI) respected
 *    arm 2: lambda ~ |cos(lambda - b)|             -- (MI) violated, p = 1 tilt
 *
 *  THREE ALGORITHMS, so the race separates language from algorithm:
 *    NAIVE  O(grid^3 * N)   — what the Python does: re-sum all N samples for
 *                             every one of grid^3 setting triples.
 *    TABLE  O(grid*N + grid^3) — note E(a,b) depends only on (a-b): substituting
 *                             u = lambda - b makes the weight and B functions of
 *                             u alone and A a function of u + b - a. So one
 *                             correlation table E(theta) serves every triple.
 *    HIST   O(N + nb^2 + grid^3) — histogram lambda once, then the whole E(theta)
 *                             table is a circular correlation of two square waves
 *                             against that histogram. The N-dependence collapses
 *                             to a single pass.
 *
 *  Build: gcc -O3 -march=native -fopenmp -o phase_bell phase_bell.c -lm
 *  Usage: phase_bell [N] [grid] [--race]
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <omp.h>

#define TAU 6.283185307179586476925287

/* xoshiro256++ — fast, good, and reproducible across runs */
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
static inline double runif(rng_t *r) { return (rnext(r) >> 11) * 0x1.0p-53; }

static double now(void) { return omp_get_wtime(); }

/* ---------------------------------------------------------------- NAIVE */
/* exactly the Python algorithm: for every (a',b,b') resum all N samples */
static double chsh_naive(const double *lam, long N, int grid, int arm) {
    double best = 0.0;
    #pragma omp parallel for schedule(dynamic) reduction(max:best)
    for (int ia = 0; ia < grid; ia++) {
        double ap = TAU * ia / grid;
        for (int ib = 0; ib < grid; ib++) {
            double b = TAU * ib / grid;
            for (int ic = 0; ic < grid; ic++) {
                double bp = TAU * ic / grid;
                double e[4] = {0, 0, 0, 0};
                double wsum[2] = {0, 0};
                for (long i = 0; i < N; i++) {
                    double L = lam[i];
                    double sb = (cos(L - b) >= 0) ? -1.0 : 1.0;    /* B at b (with the pi) */
                    double sbp = (cos(L - bp) >= 0) ? -1.0 : 1.0;
                    double s0 = (cos(L) >= 0) ? 1.0 : -1.0;        /* A at a = 0 */
                    double sa = (cos(L - ap) >= 0) ? 1.0 : -1.0;   /* A at a' */
                    double wb = (arm == 2) ? fabs(cos(L - b)) : 1.0;
                    double wbp = (arm == 2) ? fabs(cos(L - bp)) : 1.0;
                    e[0] += wb * s0 * sb;  e[1] += wbp * s0 * sbp;
                    e[2] += wb * sa * sb;  e[3] += wbp * sa * sbp;
                    wsum[0] += wb; wsum[1] += wbp;
                }
                double S = e[0] / wsum[0] + e[1] / wsum[1]
                         + e[2] / wsum[0] - e[3] / wsum[1];
                if (fabs(S) > best) best = fabs(S);
            }
        }
    }
    return best;
}

/* ---------------------------------------------------------------- TABLE */
/* E depends only on theta = a - b: build one table, then search by lookup */
static void build_table_direct(const double *lam, long N, int nth, int arm, double *E) {
    #pragma omp parallel for schedule(static)
    for (int k = 0; k < nth; k++) {
        double th = TAU * k / nth;
        double num = 0.0, den = 0.0;
        for (long i = 0; i < N; i++) {
            double u = lam[i];                       /* u = lambda - b */
            double w = (arm == 2) ? fabs(cos(u)) : 1.0;
            double B = (cos(u) >= 0) ? -1.0 : 1.0;   /* includes the pi offset */
            double A = (cos(u + th) >= 0) ? 1.0 : -1.0;
            num += w * A * B; den += w;
        }
        E[k] = num / den;
    }
}

/* ---------------------------------------------------------------- HIST */
/* one O(N) pass to histogram lambda; then the table is a circular correlation */
static void build_table_hist(const double *lam, long N, int nb, int nth, int arm,
                             double *E) {
    double *h = calloc(nb, sizeof(double));
    #pragma omp parallel
    {
        double *hl = calloc(nb, sizeof(double));
        #pragma omp for schedule(static)
        for (long i = 0; i < N; i++) {
            double u = fmod(lam[i], TAU); if (u < 0) u += TAU;
            int j = (int)(u * nb / TAU); if (j >= nb) j = nb - 1;
            hl[j] += (arm == 2) ? fabs(cos(u)) : 1.0;
        }
        #pragma omp critical
        for (int j = 0; j < nb; j++) h[j] += hl[j];
        free(hl);
    }
    double den = 0.0;
    for (int j = 0; j < nb; j++) den += h[j];
    /* B(u) and A(u+theta) are square waves; correlate against h */
    signed char *Bs = malloc(nb), *As = malloc(nb);
    for (int j = 0; j < nb; j++) {
        double u = TAU * (j + 0.5) / nb;
        Bs[j] = (cos(u) >= 0) ? -1 : 1;
        As[j] = (cos(u) >= 0) ? 1 : -1;
    }
    #pragma omp parallel for schedule(static)
    for (int k = 0; k < nth; k++) {
        int shift = (int)llround((double)k * nb / nth);
        double num = 0.0;
        for (int j = 0; j < nb; j++) {
            int jj = j + shift; if (jj >= nb) jj -= nb;
            num += h[j] * (double)As[jj] * (double)Bs[j];
        }
        E[k] = num / den;
    }
    free(h); free(Bs); free(As);
}

/* ---------------------------------------------------------------- EXACT INT
 * Angles as uint32 "turns": 2^32 == one full revolution. Modular angle
 * arithmetic is then EXACT -- no fmod, no representation of pi, and the
 * response square waves need no cos() call at all:
 *     sgn cos(u) = +1  iff  turn in [0, 1/4) U [3/4, 1)
 * For arm 1 the histogram bins are integer COUNTS and the correlation is an
 * integer sum, so E = (integer)/N is computed with ZERO rounding.
 * For arm 2 the weights are irrational; they are taken to Q30 fixed point and
 * accumulated in __int128, which removes overflow as a concern entirely.   */
static void build_table_int(const double *lam, long N, int nb, int nth, int arm,
                            double *E, long *exact_num, long *exact_den) {
    __int128 *h = calloc(nb, sizeof(__int128));
    const int SH = 30;                       /* Q30 fixed point for arm 2 */
    long long *wq = malloc(sizeof(long long) * nb);
    for (int j = 0; j < nb; j++) {
        double u = TAU * (j + 0.5) / nb;
        wq[j] = (arm == 2) ? (long long)llround(fabs(cos(u)) * (1LL << SH))
                           : 1LL;
    }
    #pragma omp parallel
    {
        __int128 *hl = calloc(nb, sizeof(__int128));
        #pragma omp for schedule(static)
        for (long i = 0; i < N; i++) {
            /* exact modular reduction into a uint32 turn */
            double frac = lam[i] / TAU;
            frac -= floor(frac);
            uint32_t turn = (uint32_t)(frac * 4294967296.0);
            int j = (int)(turn >> (32 - (int)llround(log2((double)nb))));
            hl[j] += wq[j];
        }
        #pragma omp critical
        for (int j = 0; j < nb; j++) h[j] += hl[j];
        free(hl);
    }
    __int128 den = 0;
    for (int j = 0; j < nb; j++) den += h[j];
    /* square waves without cos(): +1 on the first and last quarter */
    signed char *As = malloc(nb);
    for (int j = 0; j < nb; j++)
        As[j] = (j < nb / 4 || j >= 3 * nb / 4) ? 1 : -1;
    #pragma omp parallel for schedule(static)
    for (int k = 0; k < nth; k++) {
        int shift = (int)llround((double)k * nb / nth);
        __int128 num = 0;
        for (int j = 0; j < nb; j++) {
            int jj = j + shift; if (jj >= nb) jj -= nb;
            /* B = -As[j] (the pi offset), A = As[jj] */
            num += (__int128)h[j] * (__int128)(-(int)As[jj] * (int)As[j]);
        }
        E[k] = (double)num / (double)den;
        if (k == 0 && exact_num) { *exact_num = (long)num; *exact_den = (long)den; }
    }
    free(h); free(As); free(wq);
}

static double chsh_from_table(const double *E, int nth, int grid) {
    double best = 0.0;
    #pragma omp parallel for schedule(static) reduction(max:best)
    for (int ia = 0; ia < grid; ia++) {
        double ap = TAU * ia / grid;
        for (int ib = 0; ib < grid; ib++) {
            double b = TAU * ib / grid;
            for (int ic = 0; ic < grid; ic++) {
                double bp = TAU * ic / grid;
                #define ELOOK(th) E[((int)llround(fmod((th) + TAU * 4, TAU) * nth / TAU)) % nth]
                double S = ELOOK(0.0 - b) + ELOOK(0.0 - bp)
                         + ELOOK(ap - b) - ELOOK(ap - bp);
                if (fabs(S) > best) best = fabs(S);
            }
        }
    }
    return best;
}

int main(int argc, char **argv) {
    long N = (argc > 1) ? atol(argv[1]) : 4000000L;
    int grid = (argc > 2) ? atoi(argv[2]) : 49;
    int race = 0;
    for (int i = 1; i < argc; i++) if (!strcmp(argv[i], "--race")) race = 1;

    printf("phase_bell (C)  N=%ld  grid=%d  threads=%d\n", N, grid, omp_get_max_threads());

    double *lam = malloc(sizeof(double) * N);
    rng_t r; rseed(&r, 20260726ULL);
    for (long i = 0; i < N; i++) lam[i] = runif(&r) * TAU;

    int nth = 2880, nb = 4096;
    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "--nb") && i + 1 < argc) nb = atoi(argv[i + 1]);
        if (!strcmp(argv[i], "--nth") && i + 1 < argc) nth = atoi(argv[i + 1]);
    }
    double *E = malloc(sizeof(double) * nth);

    for (int arm = 1; arm <= 2; arm++) {
        printf("\n--- ARM %d (%s) ---\n", arm,
               arm == 1 ? "MI respected, uniform creation phase"
                        : "MI violated by design, p=1 tilt");
        double t0, S;

        long en = 0, ed = 0;
        t0 = now(); build_table_int(lam, N, nb, nth, arm, E, &en, &ed);
        double Si = chsh_from_table(E, nth, grid); double th_i = now() - t0;
        printf("  INT    max|S| = %.9f  %8.3f s   [exact E(0) = %ld/%ld = %.15f]\n",
               Si, th_i, en, ed, (double)en / (double)ed);

        t0 = now(); build_table_hist(lam, N, nb, nth, arm, E);
        double Sh = chsh_from_table(E, nth, grid); double th_h = now() - t0;
        printf("  HIST   max|S| = %.9f  %8.3f s   [nb=%d nth=%d]   f64-vs-int diff %.2e\n",
               Sh, th_h, nb, nth, fabs(Sh - Si));

        if (race) {
            t0 = now(); build_table_direct(lam, N, nth, arm, E);
            double Sd = chsh_from_table(E, nth, grid); double th_d = now() - t0;
            printf("  TABLE  max|S| = %.5f      %8.3f s\n", Sd, th_d);
        }

        if (race) {
            t0 = now(); S = chsh_naive(lam, N, grid, arm); double th_n = now() - t0;
            printf("  NAIVE  max|S| = %.5f      %8.3f s   (%.0fx slower than HIST)\n",
                   S, th_n, th_n / th_h);
        }
        printf("  E(theta) samples: ");
        for (int d = 0; d <= 180; d += 30)
            printf("%+.4f ", E[(int)llround((double)d / 360.0 * nth)]);
        printf("\n");
    }
    printf("\n  targets: arm1 triangle -> 2.00000 ; arm2 cosine -> 2sqrt2 = %.5f\n",
           2.0 * sqrt(2.0));
    free(lam); free(E);
    return 0;
}
