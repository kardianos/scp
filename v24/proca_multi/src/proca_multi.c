/*
 * proca_multi.c — Asymmetric pairwise coupling: two-scale force measurement
 *
 * Three massive scalars with triple-product AND asymmetric pairwise coupling:
 *   L = sum_a [ (1/2)(dt phi_a)^2 - (1/2)(dx phi_a)^2 - (m^2/2)phi_a^2 ]
 *     - (mu/2) P^2 / (1 + kappa P^2)
 *     - lambda_12 phi_1 phi_2 - lambda_23 phi_2 phi_3 - lambda_31 phi_3 phi_1
 *
 * Mass matrix:
 *   M^2 = | m^2      lambda_12   lambda_31  |
 *         | lambda_12 m^2        lambda_23  |
 *         | lambda_31 lambda_23  m^2        |
 *
 * Three distinct eigenvalues -> three mediator masses -> up to three ranges.
 *
 * Phase 1: Single oscillon equilibration with asymmetric coupling
 * Phase 2: Two-oscillon force measurement at multiple separations
 * Phase 3: Coupling ratio scan
 *
 * Compile: gcc -O3 -Wall -o proca_multi src/proca_multi.c -lm
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ---- Parameters ---- */
static double mu       = -20.0;
static double kappa    = 20.0;
static double mass     = 1.0;
static double lam12    = 0.85;
static double lam23    = 0.95;
static double lam31    = 0.95;
static double A_init   = 0.8;
static double sigma_w  = 3.0;
static int    Nx       = 8000;
static double xmax     = 200.0;
static double tfinal   = 10000.0;
static int    phase    = 1;
static double sep      = 20.0;
static char   outdir[512] = "data";

static void parse_args(int argc, char **argv)
{
    for (int i = 1; i < argc - 1; i += 2) {
        if      (!strcmp(argv[i], "-mu"))     mu      = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-kappa"))  kappa   = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-mass"))   mass    = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-lam12"))  lam12   = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-lam23"))  lam23   = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-lam31"))  lam31   = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-A"))      A_init  = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-sigma"))  sigma_w = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-Nx"))     Nx      = atoi(argv[i+1]);
        else if (!strcmp(argv[i], "-xmax"))   xmax    = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-tfinal")) tfinal  = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-phase"))  phase   = atoi(argv[i+1]);
        else if (!strcmp(argv[i], "-sep"))    sep     = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-o"))      strncpy(outdir, argv[i+1], sizeof(outdir)-1);
        else { fprintf(stderr, "Unknown flag: %s\n", argv[i]); exit(1); }
    }
}

/* -dV_triple/dphi_a where V = (mu/2)P^2/(1+kappa P^2), P = phi1 phi2 phi3 */
static inline double force_triple(double p1, double p2, double p3, int a)
{
    double P  = p1 * p2 * p3;
    double P2 = P * P;
    double denom2 = (1.0 + kappa * P2) * (1.0 + kappa * P2);
    double dP;
    switch (a) {
        case 0: dP = p2 * p3; break;
        case 1: dP = p1 * p3; break;
        case 2: dP = p1 * p2; break;
        default: dP = 0.0;
    }
    return -mu * P * dP / denom2;
}

/* Pairwise coupling constants for field a:
 *   -dV_pw/dphi_1 = -lambda_12 phi_2 - lambda_31 phi_3
 *   -dV_pw/dphi_2 = -lambda_12 phi_1 - lambda_23 phi_3
 *   -dV_pw/dphi_3 = -lambda_23 phi_2 - lambda_31 phi_1
 */
static inline double force_pair(double p1, double p2, double p3, int a,
                                double l12, double l23, double l31)
{
    switch (a) {
        case 0: return -l12 * p2 - l31 * p3;
        case 1: return -l12 * p1 - l23 * p3;
        case 2: return -l23 * p2 - l31 * p1;
        default: return 0.0;
    }
}

/* ---- Analytic 3x3 eigenvalues ----
 * M^2 = m^2*I + L  where L = {{0,l12,l31},{l12,0,l23},{l31,l23,0}}
 * Eigenvalues of L via cubic formula (trace=0 simplifies):
 *   det(L - eI) = -e^3 + (l12^2+l23^2+l31^2)e + 2*l12*l23*l31 = 0
 *   i.e. e^3 - p*e - q = 0 with p = l12^2+l23^2+l31^2, q = 2*l12*l23*l31
 */
static void compute_eigenvalues(double l12, double l23, double l31,
                                double m2, double evals[3])
{
    double p = l12*l12 + l23*l23 + l31*l31;
    double q = 2.0 * l12 * l23 * l31;

    /* Cardano with trigonometric form (all real roots when discriminant >= 0) */
    double disc = 4.0*p*p*p - 27.0*q*q;
    if (disc >= 0 && p > 1e-30) {
        /* Three real roots */
        double theta = acos(3.0*q*sqrt(3.0/p) / (2.0*p)) / 3.0;
        double r = 2.0 * sqrt(p / 3.0);
        evals[0] = m2 + r * cos(theta);
        evals[1] = m2 + r * cos(theta - 2.0*M_PI/3.0);
        evals[2] = m2 + r * cos(theta + 2.0*M_PI/3.0);
    } else {
        /* Fallback: direct Cardano */
        double D = q*q/4.0 - p*p*p/27.0;
        if (D >= 0) {
            double sqD = sqrt(D);
            double u = cbrt(-q/2.0 + sqD);
            double v = cbrt(-q/2.0 - sqD);
            evals[0] = m2 + u + v;
            /* Complex conjugate pair - take real part */
            evals[1] = m2 - (u+v)/2.0;
            evals[2] = evals[1]; /* degenerate */
        } else {
            /* All real (should be caught above) */
            double theta = acos(3.0*q*sqrt(3.0/p) / (2.0*p)) / 3.0;
            double r = 2.0 * sqrt(p / 3.0);
            evals[0] = m2 + r * cos(theta);
            evals[1] = m2 + r * cos(theta - 2.0*M_PI/3.0);
            evals[2] = m2 + r * cos(theta + 2.0*M_PI/3.0);
        }
    }

    /* Sort ascending */
    for (int i = 0; i < 2; i++)
        for (int j = i+1; j < 3; j++)
            if (evals[j] < evals[i]) {
                double tmp = evals[i]; evals[i] = evals[j]; evals[j] = tmp;
            }
}

/* ---- Print eigenvalue analysis ---- */
static void print_eigenvalues(double l12, double l23, double l31)
{
    double m2 = mass * mass;
    double evals[3];
    compute_eigenvalues(l12, l23, l31, m2, evals);

    printf("\n  Mass matrix eigenvalues (m^2 + lambda_i):\n");
    for (int i = 0; i < 3; i++) {
        double mi = (evals[i] > 0) ? sqrt(evals[i]) : 0.0;
        double range = (mi > 1e-10) ? 1.0 / mi : 1e30;
        int tach = (evals[i] < 0);
        printf("    mode %d: m^2_eff = %.6f  m_eff = %.6f  range = %.4f%s\n",
               i, evals[i], mi, range, tach ? "  TACHYONIC!" : "");
    }
    printf("  Sum of eigenvalues: %.6f (should be 3*m^2 = %.6f)\n",
           evals[0]+evals[1]+evals[2], 3.0*m2);
}

/* ======================================================================
 * Core evolution engine (shared by all phases)
 * init_mode: 0 = single oscillon at origin
 *            1 = two oscillons at +/- sep/2
 * Returns: arrays phi[3], vel[3] (caller must free if non-NULL)
 * ====================================================================== */
typedef struct {
    double fc_final;
    double omega_peak;
    double E_final;
    double peak_final[3];
} EvResult;

static EvResult evolve(double l12, double l23, double l31,
                       int init_mode, double separation, double t_end,
                       int quiet, const char *ts_prefix,
                       double **phi_ret, double **vel_ret)
{
    EvResult res = {0};
    double dx  = 2.0 * xmax / (Nx - 1);
    double dx2 = dx * dx;
    double m2  = mass * mass;

    /* CFL: based on heaviest mode eigenvalue */
    double evals[3];
    compute_eigenvalues(l12, l23, l31, m2, evals);
    double m2_max = evals[2]; /* largest eigenvalue */
    if (m2_max < m2) m2_max = m2;

    double kmax = M_PI / dx;
    double dt = 0.8 * 2.0 / sqrt(kmax * kmax + m2_max);
    int Nt = (int)(t_end / dt) + 1;

    if (!quiet)
        printf("  Nx=%d xmax=%.1f dx=%.5f dt=%.6f Nt=%d t_end=%.0f\n",
               Nx, xmax, dx, dt, Nt, t_end);

    double *phi[3], *vel[3], *acc[3];
    for (int a = 0; a < 3; a++) {
        phi[a] = calloc(Nx, sizeof(double));
        vel[a] = calloc(Nx, sizeof(double));
        acc[a] = calloc(Nx, sizeof(double));
    }

    /* Absorbing boundary: outer 25% */
    double *damp = malloc(Nx * sizeof(double));
    double x_abs = xmax * 0.75;
    for (int i = 0; i < Nx; i++) {
        double x = -xmax + i * dx;
        double ax = fabs(x);
        if (ax > x_abs) {
            double f = (ax - x_abs) / (xmax - x_abs);
            damp[i] = 1.0 - 0.98 * f * f;
        } else {
            damp[i] = 1.0;
        }
    }

    /* Initialize */
    if (init_mode == 0) {
        /* Single oscillon at origin */
        for (int a = 0; a < 3; a++)
            for (int i = 0; i < Nx; i++) {
                double x = -xmax + i * dx;
                phi[a][i] = A_init * exp(-x * x / (2.0 * sigma_w * sigma_w));
            }
    } else {
        /* Two oscillons at +/- sep/2 */
        double x1 = -separation / 2.0;
        double x2 = +separation / 2.0;
        for (int a = 0; a < 3; a++)
            for (int i = 0; i < Nx; i++) {
                double x = -xmax + i * dx;
                double g1 = A_init * exp(-(x-x1)*(x-x1) / (2.0 * sigma_w * sigma_w));
                double g2 = A_init * exp(-(x-x2)*(x-x2) / (2.0 * sigma_w * sigma_w));
                phi[a][i] = g1 + g2;
            }
    }

    /* Macro: acceleration with asymmetric pairwise */
    #define COMPUTE_ACC_M() do { \
        for (int a = 0; a < 3; a++) { \
            acc[a][0] = acc[a][1] = acc[a][Nx-2] = acc[a][Nx-1] = 0; \
            for (int i = 1; i < Nx - 1; i++) { \
                double lapl = (phi[a][i+1] - 2.0*phi[a][i] + phi[a][i-1]) / dx2; \
                double fp = force_triple(phi[0][i], phi[1][i], phi[2][i], a); \
                double fpw = force_pair(phi[0][i], phi[1][i], phi[2][i], a, \
                                        l12, l23, l31); \
                acc[a][i] = lapl - m2*phi[a][i] + fp + fpw; \
            } \
        } \
    } while(0)

    COMPUTE_ACC_M();

    /* Timeseries output */
    FILE *fts = NULL;
    if (ts_prefix) {
        char tspath[600];
        snprintf(tspath, sizeof(tspath), "%s/%s.tsv", outdir, ts_prefix);
        fts = fopen(tspath, "w");
        if (!fts) { fprintf(stderr, "Cannot open %s\n", tspath); ts_prefix = NULL; }
        if (fts && init_mode == 0)
            fprintf(fts, "time\tphi1_0\tphi2_0\tphi3_0\tpeak1\tpeak2\tpeak3\t"
                         "E_kin\tE_grad\tE_mass\tE_pair\tE_pot\tE_total\tfc\n");
        else if (fts)
            fprintf(fts, "time\tx_left\tx_right\tsep\tE_left\tE_right\tE_mid\n");
    }

    /* DFT storage (for single-oscillon omega measurement) */
    int max_dft = 50000;
    double *phi0_hist = malloc(max_dft * sizeof(double));
    double *t_hist = malloc(max_dft * sizeof(double));
    int n_dft = 0;

    int rec_every = Nt / 20000;
    if (rec_every < 1) rec_every = 1;
    int print_every = Nt / 40;
    if (print_every < 1) print_every = 1;
    int dft_every = Nt / max_dft;
    if (dft_every < 1) dft_every = 1;

    double core_r = 3.0 * sigma_w;
    int ic = Nx / 2;

    for (int n = 0; n <= Nt; n++) {
        double t = n * dt;

        /* DFT sampling */
        if (init_mode == 0 && n % dft_every == 0 && n_dft < max_dft) {
            phi0_hist[n_dft] = phi[0][ic];
            t_hist[n_dft] = t;
            n_dft++;
        }

        int do_rec = (n % rec_every == 0);
        int do_print = (n % print_every == 0);

        if (do_rec || do_print) {
            if (init_mode == 0) {
                /* Single oscillon diagnostics */
                double Ek = 0, Eg = 0, Em = 0, Epw = 0, Ep = 0;
                double Ecore = 0, Eall = 0;
                double peak[3] = {0};

                for (int i = 1; i < Nx - 1; i++) {
                    double x = -xmax + i * dx;
                    for (int a = 0; a < 3; a++) {
                        Ek += 0.5 * vel[a][i] * vel[a][i] * dx;
                        double dp = (phi[a][i+1] - phi[a][i-1]) / (2.0*dx);
                        Eg += 0.5 * dp * dp * dx;
                        Em += 0.5 * m2 * phi[a][i] * phi[a][i] * dx;
                        if (fabs(phi[a][i]) > peak[a]) peak[a] = fabs(phi[a][i]);
                    }

                    double pw = l12 * phi[0][i]*phi[1][i]
                              + l23 * phi[1][i]*phi[2][i]
                              + l31 * phi[2][i]*phi[0][i];
                    Epw += pw * dx;

                    double P = phi[0][i] * phi[1][i] * phi[2][i];
                    double P2 = P * P;
                    double V = 0.5 * mu * P2 / (1.0 + kappa * P2);
                    Ep += V * dx;

                    double e = 0;
                    for (int a = 0; a < 3; a++) {
                        e += 0.5*vel[a][i]*vel[a][i];
                        double dp = (phi[a][i+1] - phi[a][i-1]) / (2.0*dx);
                        e += 0.5*dp*dp + 0.5*m2*phi[a][i]*phi[a][i];
                    }
                    e += pw + V;
                    Eall += e * dx;
                    if (fabs(x) < core_r) Ecore += e * dx;
                }

                double Et = Ek + Eg + Em + Epw + Ep;
                double fc = (Eall > 1e-20) ? Ecore / Eall : 0.0;
                res.fc_final = fc;
                res.E_final = Et;
                for (int a = 0; a < 3; a++) res.peak_final[a] = peak[a];

                if (fts && do_rec)
                    fprintf(fts, "%.6f\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t"
                                 "%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6f\n",
                            t, phi[0][ic], phi[1][ic], phi[2][ic],
                            peak[0], peak[1], peak[2],
                            Ek, Eg, Em, Epw, Ep, Et, fc);

                if (!quiet && do_print)
                    printf("  t=%7.1f  p0=(%+.3f,%+.3f,%+.3f)  pk=(%.3f,%.3f,%.3f)  "
                           "E=%+.4f  fc=%.3f\n",
                           t, phi[0][ic], phi[1][ic], phi[2][ic],
                           peak[0], peak[1], peak[2], Et, fc);
            } else {
                /* Two-oscillon diagnostics: track center of energy */
                double E_L = 0, E_R = 0, E_M = 0;
                double xE_L = 0, xE_R = 0;

                for (int i = 1; i < Nx - 1; i++) {
                    double x = -xmax + i * dx;
                    double e = 0;
                    for (int a = 0; a < 3; a++) {
                        e += 0.5 * vel[a][i] * vel[a][i];
                        double dp = (phi[a][i+1] - phi[a][i-1]) / (2.0 * dx);
                        e += 0.5 * dp * dp + 0.5 * m2 * phi[a][i] * phi[a][i];
                    }
                    e += l12 * phi[0][i]*phi[1][i]
                       + l23 * phi[1][i]*phi[2][i]
                       + l31 * phi[2][i]*phi[0][i];
                    double P = phi[0][i] * phi[1][i] * phi[2][i];
                    e += 0.5 * mu * P * P / (1.0 + kappa * P * P);
                    double ed = e * dx;

                    if (x < -2.0) {
                        E_L += ed;
                        xE_L += x * ed;
                    } else if (x > 2.0) {
                        E_R += ed;
                        xE_R += x * ed;
                    } else {
                        E_M += ed;
                    }
                }

                double xcL = (E_L > 1e-20) ? xE_L / E_L : -separation/2.0;
                double xcR = (E_R > 1e-20) ? xE_R / E_R : +separation/2.0;
                double d = xcR - xcL;

                if (fts && do_rec)
                    fprintf(fts, "%.4f\t%.4f\t%.4f\t%.4f\t%.6e\t%.6e\t%.6e\n",
                            t, xcL, xcR, d, E_L, E_R, E_M);

                if (!quiet && do_print)
                    printf("  t=%7.1f  xL=%.2f  xR=%.2f  sep=%.2f  "
                           "EL=%.4f  ER=%.4f  EM=%.4e\n",
                           t, xcL, xcR, d, E_L, E_R, E_M);
            }
        }

        if (n == Nt) break;

        /* Velocity Verlet */
        for (int a = 0; a < 3; a++)
            for (int i = 1; i < Nx - 1; i++)
                vel[a][i] += 0.5 * dt * acc[a][i];
        for (int a = 0; a < 3; a++)
            for (int i = 1; i < Nx - 1; i++)
                phi[a][i] += dt * vel[a][i];
        COMPUTE_ACC_M();
        for (int a = 0; a < 3; a++)
            for (int i = 1; i < Nx - 1; i++)
                vel[a][i] += 0.5 * dt * acc[a][i];

        /* absorbing boundary */
        for (int a = 0; a < 3; a++)
            for (int i = 0; i < Nx; i++) {
                vel[a][i] *= damp[i];
                phi[a][i] *= damp[i];
            }
    }
    #undef COMPUTE_ACC_M

    if (fts) fclose(fts);

    /* DFT: find peak omega (single oscillon only) */
    if (init_mode == 0) {
        int dft_start = n_dft / 2;
        if (n_dft - dft_start > 100) {
            double T = t_hist[n_dft-1] - t_hist[dft_start];
            int nf = 500;
            double peak_pow = 0, peak_om = 0;
            for (int k = 0; k < nf; k++) {
                double omega = 3.0 * mass * k / nf;
                double re = 0, im = 0;
                for (int j = dft_start; j < n_dft; j++) {
                    double dtj = (j > dft_start) ?
                        (t_hist[j]-t_hist[j-1]) : (t_hist[dft_start+1]-t_hist[dft_start]);
                    re += phi0_hist[j] * cos(omega * t_hist[j]) * dtj;
                    im += phi0_hist[j] * sin(omega * t_hist[j]) * dtj;
                }
                double pw = (re*re + im*im) / (T*T);
                if (pw > peak_pow) { peak_pow = pw; peak_om = omega; }
            }
            res.omega_peak = peak_om;
            if (!quiet)
                printf("\n  Spectrum: peak omega = %.4f (mass = %.4f)\n", peak_om, mass);
        }
    }

    free(phi0_hist);
    free(t_hist);
    free(damp);

    /* Return field state or free */
    if (phi_ret && vel_ret) {
        for (int a = 0; a < 3; a++) {
            phi_ret[a] = phi[a];
            vel_ret[a] = vel[a];
        }
        for (int a = 0; a < 3; a++) free(acc[a]);
    } else {
        for (int a = 0; a < 3; a++) { free(phi[a]); free(vel[a]); free(acc[a]); }
    }

    return res;
}

/* ======================================================================
 * PHASE 1: Single oscillon with asymmetric coupling
 * ====================================================================== */
static void phase1_single(void)
{
    printf("=== PHASE 1: Asymmetric pairwise oscillon ===\n");
    printf("  mu=%.1f kappa=%.1f mass=%.4f A=%.3f sigma=%.3f\n",
           mu, kappa, mass, A_init, sigma_w);
    printf("  lam12=%.4f lam23=%.4f lam31=%.4f\n", lam12, lam23, lam31);
    printf("  Nx=%d xmax=%.1f tfinal=%.0f\n", Nx, xmax, tfinal);

    print_eigenvalues(lam12, lam23, lam31);

    /* Verify eigenvalues: also compute eigenvectors */
    double m2 = mass * mass;
    double evals[3];
    compute_eigenvalues(lam12, lam23, lam31, m2, evals);

    printf("\n  Predicted mediator properties:\n");
    for (int i = 0; i < 3; i++) {
        double mi = (evals[i] > 0) ? sqrt(evals[i]) : 0.0;
        double range = (mi > 1e-10) ? 1.0 / mi : 1e30;
        printf("    mode %d: m_eff=%.6f  range=%.4f\n", i, mi, range);
    }

    /* Equilibrate oscillon */
    printf("\n  Equilibrating oscillon (t=%.0f)...\n", tfinal);
    EvResult res = evolve(lam12, lam23, lam31, 0, 0, tfinal, 0,
                          "phase1_ts", NULL, NULL);

    printf("\n=== PHASE 1 RESULTS ===\n");
    printf("  fc_final = %.4f\n", res.fc_final);
    printf("  omega    = %.4f\n", res.omega_peak);
    printf("  E_final  = %.4f\n", res.E_final);
    printf("  peak     = (%.4f, %.4f, %.4f)\n",
           res.peak_final[0], res.peak_final[1], res.peak_final[2]);
    printf("  Stable?  %s (fc > 0.9)\n", (res.fc_final > 0.9) ? "YES" : "NO");

    /* Check if omega < min(sqrt(eigenvalue)) */
    double m_min = (evals[0] > 0) ? sqrt(evals[0]) : 0.0;
    printf("  omega < m_min(%.4f)?  %s (sub-threshold => oscillon)\n",
           m_min, (res.omega_peak < m_min) ? "YES" : "NO");
}

/* ======================================================================
 * PHASE 2: Force measurement at multiple separations
 * ====================================================================== */
static void phase2_force(void)
{
    printf("=== PHASE 2: Two-oscillon force measurement ===\n");
    printf("  lam12=%.4f lam23=%.4f lam31=%.4f\n", lam12, lam23, lam31);

    print_eigenvalues(lam12, lam23, lam31);

    double m2 = mass * mass;
    double evals[3];
    compute_eigenvalues(lam12, lam23, lam31, m2, evals);

    /* First equilibrate a single oscillon briefly to check stability */
    printf("\n  Pre-check: equilibrating single oscillon (t=5000)...\n");
    EvResult pre = evolve(lam12, lam23, lam31, 0, 0, 5000.0, 1, NULL, NULL, NULL);
    printf("  Pre-check: fc=%.4f, stable=%s\n", pre.fc_final,
           (pre.fc_final > 0.8) ? "YES" : "NO");
    if (pre.fc_final < 0.5) {
        printf("  ERROR: oscillon not stable, aborting force measurement\n");
        return;
    }

    /* Separations to test */
    double seps[] = {10.0, 15.0, 20.0, 30.0, 40.0};
    int n_sep = sizeof(seps) / sizeof(seps[0]);

    char fpath[600];
    snprintf(fpath, sizeof(fpath), "%s/phase2_force.tsv", outdir);
    FILE *ff = fopen(fpath, "w");
    fprintf(ff, "D\tvel_early\taccel_early\n");

    printf("\n  Running two-oscillon at separations:\n");

    double D_arr[10], F_arr[10];
    int n_good = 0;

    for (int k = 0; k < n_sep; k++) {
        double D = seps[k];
        printf("\n  --- D = %.1f ---\n", D);

        char ts_name[128];
        snprintf(ts_name, sizeof(ts_name), "phase2_D%.0f", D);

        /* Short run: enough to measure initial force, not so long
         * that oscillons move far from initial position.
         * Use t=2000 — enough for several oscillon periods (~4 each). */
        double t_total = 2000.0;

        evolve(lam12, lam23, lam31, 1, D, t_total, 1,
               ts_name, NULL, NULL);

        /* Read back the timeseries to extract force */
        char tspath[600];
        snprintf(tspath, sizeof(tspath), "%s/%s.tsv", outdir, ts_name);
        FILE *fread = fopen(tspath, "r");
        if (!fread) {
            printf("  ERROR: cannot read %s\n", tspath);
            continue;
        }

        char line[1024];
        if (!fgets(line, sizeof(line), fread)) { fclose(fread); continue; }

        double *t_data = malloc(30000 * sizeof(double));
        double *s_data = malloc(30000 * sizeof(double));
        int npts = 0;
        while (fgets(line, sizeof(line), fread) && npts < 30000) {
            double t, xl, xr, s, el, er, em;
            if (sscanf(line, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf",
                       &t, &xl, &xr, &s, &el, &er, &em) == 7) {
                t_data[npts] = t;
                s_data[npts] = s;
                npts++;
            }
        }
        fclose(fread);

        if (npts < 100) {
            printf("  Too few data points (%d)\n", npts);
            free(t_data); free(s_data);
            continue;
        }

        /* Measure force in EARLY window (first 50% of run)
         * when oscillons are still near initial positions.
         * Fit sep(t) = a0 + a1*t + a2*t^2 => accel = 2*a2 */
        int i_end = npts / 2;
        int N = i_end;
        double t0 = t_data[0];
        double S1=0, S2=0, S3=0, S4=0, Sy=0, Sty=0, St2y=0;
        for (int i = 0; i < i_end; i++) {
            double ti = t_data[i] - t0;
            double yi = s_data[i];
            S1 += ti; S2 += ti*ti; S3 += ti*ti*ti; S4 += ti*ti*ti*ti;
            Sy += yi; Sty += ti*yi; St2y += ti*ti*yi;
        }
        /* Solve 3x3: [N S1 S2; S1 S2 S3; S2 S3 S4] [a0;a1;a2] = [Sy;Sty;St2y] */
        double MV[3][4] = {
            {(double)N, S1, S2, Sy},
            {S1, S2, S3, Sty},
            {S2, S3, S4, St2y}
        };
        for (int col = 0; col < 3; col++) {
            int piv = col;
            for (int row = col+1; row < 3; row++)
                if (fabs(MV[row][col]) > fabs(MV[piv][col])) piv = row;
            for (int j = 0; j < 4; j++) {
                double tmp = MV[col][j]; MV[col][j] = MV[piv][j]; MV[piv][j] = tmp;
            }
            for (int row = col+1; row < 3; row++) {
                double fac = MV[row][col] / MV[col][col];
                for (int j = col; j < 4; j++)
                    MV[row][j] -= fac * MV[col][j];
            }
        }
        double a[3];
        for (int i = 2; i >= 0; i--) {
            a[i] = MV[i][3];
            for (int j = i+1; j < 3; j++)
                a[i] -= MV[i][j] * a[j];
            a[i] /= MV[i][i];
        }
        double vel = a[1];     /* initial velocity dsep/dt */
        double accel = 2.0 * a[2]; /* d^2sep/dt^2 = force/mass */

        printf("  D=%.1f: sep(0)=%.4f vel=%.4e accel=%.4e\n",
               D, s_data[0], vel, accel);

        fprintf(ff, "%.1f\t%.6e\t%.6e\n", D, vel, accel);

        D_arr[n_good] = D;
        F_arr[n_good] = fabs(accel); /* use |accel| for fitting */
        n_good++;

        free(t_data);
        free(s_data);
    }

    fclose(ff);

    /* Fit to two-Yukawa: F = F_s*exp(-D/R_s) + F_w*exp(-D/R_w) */
    printf("\n  --- Two-Yukawa fit ---\n");
    if (n_good >= 4) {
        /* Grid search over R_s, R_w */
        double best_err = 1e30;
        double best_Fs = 0, best_Rs = 0, best_Fw = 0, best_Rw = 0;

        for (int is = 1; is <= 50; is++) {
            double Rs = 0.5 + is * 0.5; /* 1.0 to 25.5 */
            for (int iw = is+1; iw <= 80; iw++) {
                double Rw = 0.5 + iw * 0.5; /* must be > Rs */

                /* Linear least-squares for Fs, Fw given Rs, Rw */
                /* F_k = Fs*exp(-Dk/Rs) + Fw*exp(-Dk/Rw) */
                double A11=0, A12=0, A22=0, b1=0, b2=0;
                for (int k = 0; k < n_good; k++) {
                    double es = exp(-D_arr[k] / Rs);
                    double ew = exp(-D_arr[k] / Rw);
                    A11 += es*es;
                    A12 += es*ew;
                    A22 += ew*ew;
                    b1  += F_arr[k]*es;
                    b2  += F_arr[k]*ew;
                }
                double det = A11*A22 - A12*A12;
                if (fabs(det) < 1e-30) continue;
                double Fs = (A22*b1 - A12*b2) / det;
                double Fw = (A11*b2 - A12*b1) / det;

                double err = 0;
                for (int k = 0; k < n_good; k++) {
                    double fit = Fs*exp(-D_arr[k]/Rs) + Fw*exp(-D_arr[k]/Rw);
                    double resid = fit - F_arr[k];
                    err += resid * resid;
                }
                if (err < best_err) {
                    best_err = err;
                    best_Fs = Fs; best_Rs = Rs;
                    best_Fw = Fw; best_Rw = Rw;
                }
            }
        }

        printf("  Two-Yukawa: F = %.4e * exp(-D/%.2f) + %.4e * exp(-D/%.2f)\n",
               best_Fs, best_Rs, best_Fw, best_Rw);
        printf("  Residual: %.4e\n", best_err);

        /* Compare ranges to eigenvalue predictions */
        printf("\n  Predicted ranges from eigenvalues:\n");
        for (int i = 0; i < 3; i++) {
            double mi = (evals[i] > 0) ? sqrt(evals[i]) : 0.0;
            double range = (mi > 1e-10) ? 1.0 / mi : 1e30;
            printf("    mode %d: range = %.4f\n", i, range);
        }
        printf("  Fitted ranges: R_strong = %.2f, R_weak = %.2f\n",
               best_Rs, best_Rw);
        printf("  TWO DISTINCT RANGES? %s (ratio = %.2f)\n",
               (best_Rw / best_Rs > 1.5) ? "YES" : "NO",
               best_Rw / best_Rs);

        /* Save fit */
        char fitpath[600];
        snprintf(fitpath, sizeof(fitpath), "%s/phase2_fit.tsv", outdir);
        FILE *ffit = fopen(fitpath, "w");
        fprintf(ffit, "D\tF_num\tF_fit\n");
        for (int k = 0; k < n_good; k++) {
            double fit = best_Fs*exp(-D_arr[k]/best_Rs) + best_Fw*exp(-D_arr[k]/best_Rw);
            fprintf(ffit, "%.1f\t%.6e\t%.6e\n", D_arr[k], F_arr[k], fit);
        }
        fclose(ffit);

        /* Also single-Yukawa for comparison */
        double best1_err = 1e30, best1_F = 0, best1_R = 0;
        for (int ir = 1; ir <= 100; ir++) {
            double R = 0.5 + ir * 0.5;
            double Asum = 0, Bsum = 0;
            for (int k = 0; k < n_good; k++) {
                double e = exp(-D_arr[k] / R);
                Asum += e * e;
                Bsum += F_arr[k] * e;
            }
            if (Asum < 1e-30) continue;
            double F = Bsum / Asum;
            double err = 0;
            for (int k = 0; k < n_good; k++) {
                double resid = F*exp(-D_arr[k]/R) - F_arr[k];
                err += resid * resid;
            }
            if (err < best1_err) {
                best1_err = err; best1_F = F; best1_R = R;
            }
        }
        printf("\n  Single-Yukawa: F = %.4e * exp(-D/%.2f), residual = %.4e\n",
               best1_F, best1_R, best1_err);
        printf("  Two-Yukawa improvement: %.2fx better residual\n",
               (best_err > 1e-30) ? best1_err / best_err : 0.0);
    } else {
        printf("  Not enough data points for two-Yukawa fit (need >= 4, got %d)\n", n_good);
    }

    printf("\n  Output: %s\n", fpath);
}

/* ======================================================================
 * PHASE 3: Coupling ratio scan
 * ====================================================================== */
static void phase3_scan(void)
{
    printf("=== PHASE 3: Coupling ratio scan ===\n");
    printf("  mu=%.1f kappa=%.1f mass=%.4f\n", mu, kappa, mass);

    /* (lam12, lam23=lam31) pairs — chosen for stability (no tachyonic modes)
     * Sweep from "strongly split" to "weakly split" coupling
     * Each pair satisfies: min eigenvalue of mass matrix > 0 */
    double ratios[][2] = {
        {0.90, 0.95},   /* small split: ranges 5.49, 3.16, 0.59 */
        {0.85, 0.95},   /* moderate split: ranges 7.94, 2.58, 0.59 */
        {0.80, 0.90},   /* wider split: ranges 3.90, 2.24, 0.60 */
        {0.70, 0.85},   /* large split: ranges 3.19, 1.83, 0.62 */
        {0.50, 0.75},   /* very large split: ranges 2.50, 1.41, 0.65 */
    };
    int n_ratio = sizeof(ratios) / sizeof(ratios[0]);

    char path[600];
    snprintf(path, sizeof(path), "%s/phase3_scan.tsv", outdir);
    FILE *fp = fopen(path, "w");
    fprintf(fp, "lam12\tlam23\tlam31\teval0\teval1\teval2\t"
                "m0\tm1\tm2\trange0\trange1\trange2\t"
                "fc\tomega\tstable\ttachyonic\n");

    double m2 = mass * mass;
    double t_equil = 5000.0;

    for (int k = 0; k < n_ratio; k++) {
        double l12 = ratios[k][0];
        double l23 = ratios[k][1];
        double l31 = l23; /* symmetric: lam23 = lam31 */

        printf("\n=== Ratio: lam12=%.2f, lam23=lam31=%.2f ===\n", l12, l23);

        double evals[3];
        compute_eigenvalues(l12, l23, l31, m2, evals);

        int tachyonic = 0;
        for (int i = 0; i < 3; i++)
            if (evals[i] < 0) tachyonic = 1;

        print_eigenvalues(l12, l23, l31);

        if (tachyonic) {
            printf("  TACHYONIC — skipping evolution\n");
            fprintf(fp, "%.4f\t%.4f\t%.4f\t%.6f\t%.6f\t%.6f\t"
                        "0\t0\t0\t0\t0\t0\t0\t0\t0\t1\n",
                    l12, l23, l31, evals[0], evals[1], evals[2]);
            continue;
        }

        char ts_name[128];
        snprintf(ts_name, sizeof(ts_name), "phase3_l12_%.2f_l23_%.2f", l12, l23);

        EvResult res = evolve(l12, l23, l31, 0, 0, t_equil, 0,
                              ts_name, NULL, NULL);

        int stable = (res.fc_final > 0.9);
        printf("  => fc=%.4f omega=%.4f stable=%s\n",
               res.fc_final, res.omega_peak, stable ? "YES" : "NO");

        fprintf(fp, "%.4f\t%.4f\t%.4f", l12, l23, l31);
        for (int i = 0; i < 3; i++) fprintf(fp, "\t%.6f", evals[i]);
        for (int i = 0; i < 3; i++) {
            double mi = (evals[i] > 0) ? sqrt(evals[i]) : 0.0;
            fprintf(fp, "\t%.6f", mi);
        }
        for (int i = 0; i < 3; i++) {
            double mi = (evals[i] > 0) ? sqrt(evals[i]) : 0.0;
            double range = (mi > 1e-10) ? 1.0 / mi : 1e30;
            fprintf(fp, "\t%.6f", range);
        }
        fprintf(fp, "\t%.4f\t%.4f\t%d\t%d\n",
                res.fc_final, res.omega_peak, stable, tachyonic);
    }

    fclose(fp);
    printf("\n  Output: %s\n", path);
}

/* ====================================================================== */
int main(int argc, char **argv)
{
    parse_args(argc, argv);

    switch (phase) {
        case 1: phase1_single();  break;
        case 2: phase2_force();   break;
        case 3: phase3_scan();    break;
        default:
            fprintf(stderr, "Unknown phase %d (use 1, 2, or 3)\n", phase);
            return 1;
    }
    return 0;
}
