/*
 * critical_density_1d.c -- 1D critical density scanner
 *
 * For each point in (alpha, beta) parameter space, binary-searches the
 * initial blob amplitude A to find the critical value A_crit separating
 * dispersal from collapse.
 *
 * Three-field equation (symmetric ansatz phi_1=phi_2=phi_3=phi, P=phi^3):
 *   d^2 phi/dt^2 = d^2 phi/dx^2 - m_eff^2 * phi - V'(P)
 *   V(P) = (mu/2) P^2 / (1 + kappa P^2)
 *   V'(phi) = dV/dphi = 3 mu phi^5 / (1 + kappa phi^6)^2
 *
 * Mode 1 (pure inverse): m_eff^2 = alpha / (1 + beta * 3*phi^2)
 * Mode 2 (hybrid):       m_eff^2 = m_const^2 + alpha / (1 + beta * 3*phi^2)
 * Mode 3 (density kappa): m_eff^2 = alpha / (1 + beta * 3*phi^2)
 *                          kappa_eff = kappa / (1 + gamma * 3*phi^2)
 *
 * Collapse criterion: |E_pot(T)| > |E_pot(0)| * threshold  (binding grew)
 *
 * Build: gcc -O3 -fopenmp -o critical_density_1d src/critical_density_1d.c -lm
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

/* Standard physics parameters */
static double MU       = -41.345;
static double KAPPA    = 50.0;
static double M_CONST2 = 2.25;     /* constant mass^2 for mode 2 (hybrid) */
static double GAMMA    = 0.0;      /* density-dependent kappa strength (mode 3) */

/* Phase offsets for the three-field symmetric ansatz (reference) */
/* delta = {0.0, 3.0005, 4.4325} */

/* Scan parameters */
static double ALPHA_MIN   = 0.5;
static double ALPHA_MAX   = 5.0;
static int    ALPHA_STEPS = 20;
static double BETA_MIN    = 0.5;
static double BETA_MAX    = 20.0;
static int    BETA_STEPS  = 20;

/* Gamma scan for mode 3 */
static double GAMMA_MIN   = 0.1;
static double GAMMA_MAX   = 10.0;
static int    GAMMA_STEPS = 20;

/* Grid and evolution */
static int    NX     = 1024;
static double L      = 30.0;
static double T_EVOL = 50.0;
static double W_BLOB = 5.0;        /* blob half-width */

/* Binary search */
static double A_LOW_INIT  = 0.01;
static double A_HIGH_INIT = 2.0;
static int    BISECT_ITER = 20;
static double COLLAPSE_THRESHOLD = 1.5; /* |E_pot(T)/E_pot(0)| > this => collapse */

/* Mode */
static int MODE = 1;

/* Output */
static char OUTPATH[512] = "phase_diagram.tsv";

static void parse_args(int argc, char **argv)
{
    for (int i = 1; i < argc; i++) {
        if      (!strcmp(argv[i], "-mode"))        MODE          = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-alpha_min"))   ALPHA_MIN     = atof(argv[++i]);
        else if (!strcmp(argv[i], "-alpha_max"))   ALPHA_MAX     = atof(argv[++i]);
        else if (!strcmp(argv[i], "-alpha_steps")) ALPHA_STEPS   = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-beta_min"))    BETA_MIN      = atof(argv[++i]);
        else if (!strcmp(argv[i], "-beta_max"))    BETA_MAX      = atof(argv[++i]);
        else if (!strcmp(argv[i], "-beta_steps"))  BETA_STEPS    = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-gamma_min"))   GAMMA_MIN     = atof(argv[++i]);
        else if (!strcmp(argv[i], "-gamma_max"))   GAMMA_MAX     = atof(argv[++i]);
        else if (!strcmp(argv[i], "-gamma_steps")) GAMMA_STEPS   = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-W"))           W_BLOB        = atof(argv[++i]);
        else if (!strcmp(argv[i], "-T"))           T_EVOL        = atof(argv[++i]);
        else if (!strcmp(argv[i], "-N"))           NX            = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-L"))           L             = atof(argv[++i]);
        else if (!strcmp(argv[i], "-A_low"))       A_LOW_INIT    = atof(argv[++i]);
        else if (!strcmp(argv[i], "-A_high"))      A_HIGH_INIT   = atof(argv[++i]);
        else if (!strcmp(argv[i], "-bisect"))      BISECT_ITER   = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-threshold"))   COLLAPSE_THRESHOLD = atof(argv[++i]);
        else if (!strcmp(argv[i], "-mu"))          MU            = atof(argv[++i]);
        else if (!strcmp(argv[i], "-kappa"))       KAPPA         = atof(argv[++i]);
        else if (!strcmp(argv[i], "-m_const"))     M_CONST2      = atof(argv[++i]);
        else if (!strcmp(argv[i], "-gamma"))       GAMMA         = atof(argv[++i]);
        else if (!strcmp(argv[i], "-o"))           strncpy(OUTPATH, argv[++i], sizeof(OUTPATH)-1);
        else { fprintf(stderr, "Unknown flag: %s\n", argv[i]); exit(1); }
    }
}

/*
 * Run a single evolution at given (alpha, beta, gamma, A) and return
 * the ratio |E_pot(T)| / |E_pot(0)|.
 *
 * Returns the E_pot ratio (>1 means collapse grew the binding).
 * phi, vel, acc are pre-allocated work arrays of size NX.
 */
static double run_evolution(double alpha, double beta, double gamma_val,
                            double A, int mode,
                            double *phi, double *vel, double *acc, double *damp)
{
    const double dx  = 2.0 * L / (NX - 1);
    const double dx2 = dx * dx;

    /* Initialize blob: phi_a = A * cos(delta_a) for |x| < W, smooth edges */
    memset(phi, 0, NX * sizeof(double));
    memset(vel, 0, NX * sizeof(double));
    memset(acc, 0, NX * sizeof(double));

    /*
     * Symmetric ansatz: phi_1 = phi_2 = phi_3 = phi.
     * The effective single-field amplitude is the "common field" value.
     * In three-field product P = phi1*phi2*phi3, with phi_a = A*cos(delta_a)*env,
     * we work with the single effective field whose amplitude encodes
     * the three-field structure via the potential.
     *
     * For the purpose of collapse detection we use a single field phi
     * with the EOM:
     *   phi'' = phi_xx - m_eff^2 phi - 3 mu phi^5 / (1 + kappa phi^6)^2
     * where the factor 3 comes from dP/dphi_a with P=phi^3.
     *
     * Sigma_phi2 = 3 phi^2  (three equal fields).
     */
    for (int i = 0; i < NX; i++) {
        double x = -L + i * dx;
        if (fabs(x) < W_BLOB) {
            /* Smooth (cosine) taper at edges */
            double env = 1.0;
            double edge = W_BLOB - 1.0; /* taper starts 1 unit from edge */
            if (fabs(x) > edge && edge > 0) {
                double s = (fabs(x) - edge) / (W_BLOB - edge);
                env = 0.5 * (1.0 + cos(M_PI * s));
            }
            phi[i] = A * env;
        }
    }

    /* CFL: for standard wave eq + mass, dt < dx/c where c=1.
     * With the potential nonlinearity, use safety factor. */
    double dt = 0.4 * dx;

    /* Estimate max effective mass for CFL tightening */
    double m2_max = (mode == 2) ? M_CONST2 + alpha : alpha;
    if (m2_max > 0) {
        double kmax = M_PI / dx;
        double omega_max = sqrt(kmax * kmax + m2_max);
        double dt_cfl = 1.8 / omega_max; /* Verlet CFL */
        if (dt_cfl < dt) dt = dt_cfl;
    }

    int Nt = (int)(T_EVOL / dt) + 1;

    /* Absorbing layer: outer 15% of domain */
    double x_abs = L * 0.85;
    for (int i = 0; i < NX; i++) {
        double x = -L + i * dx;
        double ax = fabs(x);
        if (ax > x_abs) {
            double f = (ax - x_abs) / (L - x_abs);
            damp[i] = 1.0 - 0.98 * f * f;
        } else {
            damp[i] = 1.0;
        }
    }

    /* Compute initial E_pot */
    double E_pot_init = 0.0;
    for (int i = 0; i < NX; i++) {
        double p = phi[i];
        double p6 = p * p * p * p * p * p;
        double kap = KAPPA;
        if (mode == 3 && gamma_val > 0) {
            double sigma_phi2 = 3.0 * p * p;
            kap = KAPPA / (1.0 + gamma_val * sigma_phi2);
        }
        double V = 0.5 * MU * p6 / (1.0 + kap * p6);
        E_pot_init += V * dx;
    }

    /* Macro: compute acceleration for all interior points */
    #define COMPUTE_ACC_1D() do { \
        for (int ii = 1; ii < NX - 1; ii++) { \
            double lapl = (phi[ii+1] - 2.0*phi[ii] + phi[ii-1]) / dx2; \
            double p = phi[ii]; \
            double p2 = p * p; \
            double sigma_phi2 = 3.0 * p2; /* three equal fields */ \
            \
            /* Effective mass */ \
            double m2eff; \
            if (mode == 1) { \
                m2eff = alpha / (1.0 + beta * sigma_phi2); \
            } else if (mode == 2) { \
                m2eff = M_CONST2 + alpha / (1.0 + beta * sigma_phi2); \
            } else { /* mode 3 */ \
                m2eff = alpha / (1.0 + beta * sigma_phi2); \
            } \
            \
            /* Potential force: -dV/dphi */ \
            /* V = (mu/2) phi^6 / (1 + kappa_eff phi^6) */ \
            /* dV/dphi = 3 mu phi^5 / (1 + kappa_eff phi^6)^2 */ \
            double p5 = p2 * p2 * p; \
            double p6 = p5 * p; \
            double kap = KAPPA; \
            if (mode == 3 && gamma_val > 0) { \
                kap = KAPPA / (1.0 + gamma_val * sigma_phi2); \
            } \
            double denom = 1.0 + kap * p6; \
            double dVdp = 3.0 * MU * p5 / (denom * denom); \
            \
            acc[ii] = lapl - m2eff * p - dVdp; \
        } \
        acc[0] = 0.0; \
        acc[NX-1] = 0.0; \
    } while (0)

    COMPUTE_ACC_1D();

    /* Velocity Verlet evolution */
    for (int n = 0; n < Nt; n++) {
        /* Half-kick */
        for (int i = 1; i < NX - 1; i++)
            vel[i] += 0.5 * dt * acc[i];
        /* Drift */
        for (int i = 1; i < NX - 1; i++)
            phi[i] += dt * vel[i];
        /* Recompute forces */
        COMPUTE_ACC_1D();
        /* Half-kick */
        for (int i = 1; i < NX - 1; i++)
            vel[i] += 0.5 * dt * acc[i];
        /* Absorbing boundary */
        for (int i = 0; i < NX; i++) {
            vel[i] *= damp[i];
            phi[i] *= damp[i];
        }
    }

    #undef COMPUTE_ACC_1D

    /* Compute final E_pot */
    double E_pot_final = 0.0;
    for (int i = 0; i < NX; i++) {
        double p = phi[i];
        double p6 = p * p * p * p * p * p;
        double kap = KAPPA;
        if (mode == 3 && gamma_val > 0) {
            double sigma_phi2 = 3.0 * p * p;
            kap = KAPPA / (1.0 + gamma_val * sigma_phi2);
        }
        double V = 0.5 * MU * p6 / (1.0 + kap * p6);
        E_pot_final += V * dx;
    }

    /* Return ratio of |E_pot|. mu < 0 so E_pot < 0 (binding). */
    double ratio;
    if (fabs(E_pot_init) < 1e-30)
        ratio = (fabs(E_pot_final) > 1e-20) ? 1e6 : 1.0;
    else
        ratio = fabs(E_pot_final) / fabs(E_pot_init);

    return ratio;
}

/*
 * Binary search for A_crit at a given (alpha, beta, gamma).
 * Returns A_crit and the E_pot ratio at A_crit via *ratio_out.
 */
static double find_A_crit(double alpha, double beta, double gamma_val, int mode,
                           double *ratio_out)
{
    /* Per-thread work arrays */
    double *phi  = malloc(NX * sizeof(double));
    double *vel  = malloc(NX * sizeof(double));
    double *acc  = malloc(NX * sizeof(double));
    double *damp = malloc(NX * sizeof(double));

    double A_low  = A_LOW_INIT;
    double A_high = A_HIGH_INIT;

    /* First check: does A_high actually collapse? If not, report A_crit = A_high */
    double r_high = run_evolution(alpha, beta, gamma_val, A_high, mode,
                                  phi, vel, acc, damp);
    if (r_high < COLLAPSE_THRESHOLD) {
        /* Even max amplitude doesn't collapse */
        *ratio_out = r_high;
        free(phi); free(vel); free(acc); free(damp);
        return A_high; /* saturated: no collapse found */
    }

    /* Check: does A_low disperse? If it collapses, report A_crit = A_low */
    double r_low = run_evolution(alpha, beta, gamma_val, A_low, mode,
                                 phi, vel, acc, damp);
    if (r_low >= COLLAPSE_THRESHOLD) {
        /* Even min amplitude collapses */
        *ratio_out = r_low;
        free(phi); free(vel); free(acc); free(damp);
        return A_low;
    }

    /* Binary search */
    double ratio_mid = 1.0;
    for (int iter = 0; iter < BISECT_ITER; iter++) {
        double A_mid = 0.5 * (A_low + A_high);
        ratio_mid = run_evolution(alpha, beta, gamma_val, A_mid, mode,
                                  phi, vel, acc, damp);

        if (ratio_mid > COLLAPSE_THRESHOLD) {
            A_high = A_mid; /* collapse: lower the ceiling */
        } else {
            A_low = A_mid;  /* dispersal: raise the floor */
        }
    }

    *ratio_out = ratio_mid;
    free(phi); free(vel); free(acc); free(damp);
    return 0.5 * (A_low + A_high);
}

int main(int argc, char **argv)
{
    parse_args(argc, argv);

    printf("=== Critical Density 1D Scanner ===\n");
    printf("Mode %d: ", MODE);
    switch (MODE) {
        case 1: printf("pure inverse coupling m_eff^2 = alpha/(1+beta*3*phi^2)\n"); break;
        case 2: printf("hybrid m_eff^2 = %.4f + alpha/(1+beta*3*phi^2)\n", M_CONST2); break;
        case 3: printf("density-dependent kappa: kappa_eff = %.1f/(1+gamma*3*phi^2)\n", KAPPA); break;
        default: fprintf(stderr, "Unknown mode %d\n", MODE); return 1;
    }
    printf("mu=%.3f kappa=%.1f W=%.1f T=%.1f N=%d L=%.1f\n", MU, KAPPA, W_BLOB, T_EVOL, NX, L);
    printf("A search: [%.3f, %.3f] x %d bisections, threshold=%.2f\n",
           A_LOW_INIT, A_HIGH_INIT, BISECT_ITER, COLLAPSE_THRESHOLD);

    double dx = 2.0 * L / (NX - 1);
    double dt_est = 0.4 * dx;
    int Nt_est = (int)(T_EVOL / dt_est) + 1;
    printf("dx=%.5f dt~%.5f Nt~%d per evolution\n", dx, dt_est, Nt_est);

    FILE *fout = fopen(OUTPATH, "w");
    if (!fout) {
        fprintf(stderr, "Cannot open %s for writing\n", OUTPATH);
        return 1;
    }

    if (MODE == 3) {
        /* Mode 3: scan (alpha, beta) grid at each gamma, plus gamma sweep at fixed (alpha, beta) */
        fprintf(fout, "alpha\tbeta\tgamma\tA_crit\tE_pot_ratio\n");

        /* For mode 3, also scan gamma at several fixed (alpha, beta) points */
        int done = 0;

        printf("\n--- Phase A: (alpha, beta) grid at gamma=%.2f ---\n", GAMMA);
        printf("Grid: alpha [%.2f, %.2f] x %d, beta [%.2f, %.2f] x %d\n",
               ALPHA_MIN, ALPHA_MAX, ALPHA_STEPS, BETA_MIN, BETA_MAX, BETA_STEPS);

        double wall0 = omp_get_wtime();

        /* Flatten the 2D grid for OpenMP */
        int n_ab = ALPHA_STEPS * BETA_STEPS;

        #pragma omp parallel for schedule(dynamic, 1)
        for (int idx = 0; idx < n_ab; idx++) {
            int ia = idx / BETA_STEPS;
            int ib = idx % BETA_STEPS;

            double alpha = (ALPHA_STEPS > 1)
                ? ALPHA_MIN + ia * (ALPHA_MAX - ALPHA_MIN) / (ALPHA_STEPS - 1)
                : ALPHA_MIN;
            double beta = (BETA_STEPS > 1)
                ? BETA_MIN + ib * (BETA_MAX - BETA_MIN) / (BETA_STEPS - 1)
                : BETA_MIN;

            double ratio;
            double Ac = find_A_crit(alpha, beta, GAMMA, MODE, &ratio);

            #pragma omp critical
            {
                fprintf(fout, "%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n",
                        alpha, beta, GAMMA, Ac, ratio);
                fflush(fout);
                done++;
                if (done % 10 == 0 || done == n_ab) {
                    double wall = omp_get_wtime() - wall0;
                    printf("  [%d/%d] alpha=%.2f beta=%.2f -> A_crit=%.4f ratio=%.3f (%.0fs)\n",
                           done, n_ab, alpha, beta, Ac, ratio, wall);
                }
            }
        }

        /* Gamma sweep at a few representative (alpha, beta) points */
        printf("\n--- Phase B: gamma sweep ---\n");
        double ab_points[][2] = {
            {1.0, 2.0}, {2.0, 5.0}, {3.0, 10.0}, {5.0, 15.0}
        };
        int n_ab_pts = 4;

        for (int p = 0; p < n_ab_pts; p++) {
            double alpha = ab_points[p][0];
            double beta  = ab_points[p][1];
            printf("  Sweeping gamma at alpha=%.1f beta=%.1f\n", alpha, beta);

            #pragma omp parallel for schedule(dynamic, 1)
            for (int ig = 0; ig < GAMMA_STEPS; ig++) {
                double gv = (GAMMA_STEPS > 1)
                    ? GAMMA_MIN + ig * (GAMMA_MAX - GAMMA_MIN) / (GAMMA_STEPS - 1)
                    : GAMMA_MIN;

                double ratio;
                double Ac = find_A_crit(alpha, beta, gv, MODE, &ratio);

                #pragma omp critical
                {
                    fprintf(fout, "%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n",
                            alpha, beta, gv, Ac, ratio);
                    fflush(fout);
                    printf("    gamma=%.3f -> A_crit=%.4f ratio=%.3f\n", gv, Ac, ratio);
                }
            }
        }

    } else {
        /* Mode 1 or 2: scan (alpha, beta) grid */
        fprintf(fout, "alpha\tbeta\tA_crit\tE_pot_ratio\n");

        int n_total = ALPHA_STEPS * BETA_STEPS;
        int done = 0;

        printf("Grid: alpha [%.2f, %.2f] x %d, beta [%.2f, %.2f] x %d = %d points\n",
               ALPHA_MIN, ALPHA_MAX, ALPHA_STEPS,
               BETA_MIN, BETA_MAX, BETA_STEPS, n_total);

        double wall0 = omp_get_wtime();

        #pragma omp parallel for schedule(dynamic, 1)
        for (int idx = 0; idx < n_total; idx++) {
            int ia = idx / BETA_STEPS;
            int ib = idx % BETA_STEPS;

            double alpha = (ALPHA_STEPS > 1)
                ? ALPHA_MIN + ia * (ALPHA_MAX - ALPHA_MIN) / (ALPHA_STEPS - 1)
                : ALPHA_MIN;
            double beta = (BETA_STEPS > 1)
                ? BETA_MIN + ib * (BETA_MAX - BETA_MIN) / (BETA_STEPS - 1)
                : BETA_MIN;

            double ratio;
            double Ac = find_A_crit(alpha, beta, 0.0, MODE, &ratio);

            #pragma omp critical
            {
                fprintf(fout, "%.6f\t%.6f\t%.6f\t%.6f\n", alpha, beta, Ac, ratio);
                fflush(fout);
                done++;
                if (done % 10 == 0 || done == n_total) {
                    double wall = omp_get_wtime() - wall0;
                    double eta = (done > 0) ? wall * (n_total - done) / done : 0;
                    printf("  [%d/%d] alpha=%.2f beta=%.2f -> A_crit=%.4f ratio=%.3f "
                           "(%.0fs elapsed, ~%.0fs remaining)\n",
                           done, n_total, alpha, beta, Ac, ratio, wall, eta);
                }
            }
        }

        double wall = omp_get_wtime() - wall0;
        printf("\nComplete: %d points in %.1f s (%.3f s/point)\n",
               n_total, wall, wall / n_total);
    }

    fclose(fout);
    printf("\nOutput: %s\n", OUTPATH);

    /* Print summary statistics */
    printf("\n=== SUMMARY ===\n");
    printf("Mode %d, mu=%.3f, kappa=%.1f, W=%.1f, T=%.1f\n",
           MODE, MU, KAPPA, W_BLOB, T_EVOL);
    printf("Grid: N=%d, L=%.1f, dx=%.5f\n", NX, L, dx);
    printf("Bisection: %d iterations, threshold=%.2f\n",
           BISECT_ITER, COLLAPSE_THRESHOLD);
    printf("Phase diagram written to %s\n", OUTPATH);

    return 0;
}
