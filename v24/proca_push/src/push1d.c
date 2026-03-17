/*
 * push1d.c — Push lambda to extreme: how close to m^2 can we get?
 *
 * Three massive scalars with pairwise + triple-product coupling:
 *   L = sum_a [ (1/2)(dt phi_a)^2 - (1/2)(dx phi_a)^2 - (m^2/2)phi_a^2 ]
 *     - lambda(phi_1 phi_2 + phi_2 phi_3 + phi_3 phi_1)
 *     - (mu/2) P^2 / (1 + kappa P^2)
 *
 * EOM: d^2 phi_a/dt^2 = d^2 phi_a/dx^2 - m^2 phi_a - lambda(phi_b+phi_c)
 *      - mu P dP/dphi_a / (1+kappa P^2)^2
 *
 * Mass spectrum:
 *   antisymmetric: m^2_A = m^2 - lambda  -> 0 as lambda -> m^2
 *   symmetric:     m^2_S = m^2 + 2*lambda
 *
 * Phase 1: Scan lambda = {0.999, 0.9995, 0.9999, 0.99995, 0.99999}
 *          For each: evolve t=10000, check oscillon survival + vacuum stability
 * Phase 2: Two oscillons at D=50,100 at most extreme stable lambda
 *
 * Compile: gcc -O3 -Wall -o push1d src/push1d.c -lm
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ---- Parameters ---- */
static double mu      = -20.0;
static double kappa   = 20.0;
static double mass    = 1.0;
static double lambda  = 0.0;
static double A_init  = 0.8;
static double sigma_w = 3.0;
static int    Nx      = 16000;
static double xmax    = 500.0;
static double tfinal  = 10000.0;
static int    phase   = 1;
static double sep     = 50.0;
static char   outdir[512] = "v24/proca_push/data";

static void parse_args(int argc, char **argv)
{
    for (int i = 1; i < argc - 1; i += 2) {
        if      (!strcmp(argv[i], "-mu"))     mu      = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-kappa"))  kappa   = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-mass"))   mass    = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-lambda")) lambda  = atof(argv[i+1]);
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

/* Compute energy density at grid point i (needs phi, vel arrays and grid info) */
static double energy_density(double *phi[3], double *vel[3], int i,
                             double m2, double lam, double dx)
{
    double e = 0;
    for (int a = 0; a < 3; a++) {
        e += 0.5 * vel[a][i] * vel[a][i];
        double dp = (phi[a][i+1] - phi[a][i-1]) / (2.0 * dx);
        e += 0.5 * dp * dp + 0.5 * m2 * phi[a][i] * phi[a][i];
    }
    e += lam * (phi[0][i]*phi[1][i] + phi[1][i]*phi[2][i]
               + phi[2][i]*phi[0][i]);
    double P = phi[0][i] * phi[1][i] * phi[2][i];
    e += 0.5 * mu * P * P / (1.0 + kappa * P * P);
    return e;
}

/* ======================================================================
 * PHASE 1: Extreme lambda scan
 * ====================================================================== */
static void phase1_scan(void)
{
    double lam_vals[] = {0.999, 0.9995, 0.9999, 0.99995, 0.99999};
    int n_lam = sizeof(lam_vals) / sizeof(lam_vals[0]);
    double m2 = mass * mass;

    printf("=== PHASE 1: Extreme Lambda Scan ===\n");
    printf("  mu=%.1f kappa=%.1f mass=%.4f A=%.3f sigma=%.3f\n",
           mu, kappa, mass, A_init, sigma_w);
    printf("  Nx=%d xmax=%.1f tfinal=%.0f\n", Nx, xmax, tfinal);

    char sumpath[600];
    snprintf(sumpath, sizeof(sumpath), "%s/phase1_summary.tsv", outdir);
    FILE *fsum = fopen(sumpath, "w");
    if (!fsum) { fprintf(stderr, "Cannot open %s\n", sumpath); return; }
    fprintf(fsum, "lambda\tm2_anti\tm_A\tpred_range\tfc_final\tpeak_amp\tvac_max\t"
                  "E_final\tE_init\tdE_frac\tomega\tstable\n");

    double lam_best = 0;
    double fc_best  = 0;

    for (int k = 0; k < n_lam; k++) {
        double lam = lam_vals[k];
        double m2_anti = m2 - lam;
        double m_A = sqrt(m2_anti);
        double pred_range = 1.0 / m_A;

        printf("\n========================================\n");
        printf("  lambda=%.6f  m2_anti=%.6e  m_A=%.6f  range=%.1f\n",
               lam, m2_anti, m_A, pred_range);
        printf("========================================\n");

        double dx  = 2.0 * xmax / (Nx - 1);
        double dx2 = dx * dx;

        /* CFL based on heaviest mode */
        double m2_sym = m2 + 2.0 * lam;
        double kmax_grid = M_PI / dx;
        double dt = 0.8 * 2.0 / sqrt(kmax_grid * kmax_grid + m2_sym);
        int Nt = (int)(tfinal / dt) + 1;

        printf("  dx=%.6f dt=%.6f Nt=%d\n", dx, dt, Nt);

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

        /* Initialize: symmetric oscillon (all fields equal) */
        for (int a = 0; a < 3; a++)
            for (int i = 0; i < Nx; i++) {
                double x = -xmax + i * dx;
                phi[a][i] = A_init * exp(-x * x / (2.0 * sigma_w * sigma_w));
                vel[a][i] = 0.0;
            }

        /* Compute acceleration macro */
        #define COMPUTE_ACC() do { \
            for (int a = 0; a < 3; a++) { \
                acc[a][0] = acc[a][1] = acc[a][Nx-2] = acc[a][Nx-1] = 0; \
                int b = (a+1)%3, c_idx = (a+2)%3; \
                for (int i = 1; i < Nx - 1; i++) { \
                    double lapl = (phi[a][i+1] - 2.0*phi[a][i] + phi[a][i-1]) / dx2; \
                    double fp = force_triple(phi[0][i], phi[1][i], phi[2][i], a); \
                    acc[a][i] = lapl - m2*phi[a][i] - lam*(phi[b][i]+phi[c_idx][i]) + fp; \
                } \
            } \
        } while(0)

        COMPUTE_ACC();

        /* DFT storage for frequency measurement */
        int max_dft = 50000;
        double *phi0_hist = malloc(max_dft * sizeof(double));
        double *t_hist = malloc(max_dft * sizeof(double));
        int n_dft = 0;
        int dft_every = Nt / max_dft;
        if (dft_every < 1) dft_every = 1;

        double core_r = 3.0 * sigma_w;
        int ic = Nx / 2;
        int print_every = Nt / 20;
        if (print_every < 1) print_every = 1;

        double E_init = 0, E_final = 0;
        double fc_final = 0;
        double peak_amp = 0;
        double vac_max = 0; /* max |phi| at |x| > 50 */

        /* Time series file per lambda */
        char tspath[600];
        snprintf(tspath, sizeof(tspath), "%s/phase1_lam%.6f_ts.tsv", outdir, lam);
        FILE *fts = fopen(tspath, "w");
        fprintf(fts, "time\tphi0_center\tpeak\tfc\tE_total\tvac_max\n");

        int rec_every = Nt / 5000;
        if (rec_every < 1) rec_every = 1;

        for (int n = 0; n <= Nt; n++) {
            double t = n * dt;

            /* DFT sampling */
            if (n % dft_every == 0 && n_dft < max_dft) {
                phi0_hist[n_dft] = phi[0][ic];
                t_hist[n_dft] = t;
                n_dft++;
            }

            int do_rec = (n % rec_every == 0);
            int do_print = (n % print_every == 0);

            if (do_rec || do_print || n == 0 || n == Nt) {
                double Eall = 0, Ecore = 0;
                double pk = 0, vm = 0;

                for (int i = 1; i < Nx - 1; i++) {
                    double x = -xmax + i * dx;
                    double e = energy_density(phi, vel, i, m2, lam, dx);
                    Eall += e * dx;
                    if (fabs(x) < core_r) Ecore += e * dx;

                    /* Peak amplitude */
                    for (int a = 0; a < 3; a++)
                        if (fabs(phi[a][i]) > pk) pk = fabs(phi[a][i]);

                    /* Vacuum stability: max |phi| at |x| > 50 */
                    if (fabs(x) > 50.0)
                        for (int a = 0; a < 3; a++)
                            if (fabs(phi[a][i]) > vm) vm = fabs(phi[a][i]);
                }

                fc_final = (Eall > 1e-20) ? Ecore / Eall : 0.0;
                peak_amp = pk;
                vac_max = vm;
                E_final = Eall;
                if (n == 0) E_init = Eall;

                if (do_rec)
                    fprintf(fts, "%.4f\t%.6e\t%.6e\t%.4f\t%.6e\t%.6e\n",
                            t, phi[0][ic], pk, fc_final, Eall, vm);

                if (do_print)
                    printf("  t=%7.1f  phi0=%.4f  pk=%.4f  fc=%.4f  E=%.4f  vac=%.2e\n",
                           t, phi[0][ic], pk, fc_final, Eall, vm);
            }

            if (n == Nt) break;

            /* Velocity Verlet */
            for (int a = 0; a < 3; a++)
                for (int i = 1; i < Nx - 1; i++)
                    vel[a][i] += 0.5 * dt * acc[a][i];
            for (int a = 0; a < 3; a++)
                for (int i = 1; i < Nx - 1; i++)
                    phi[a][i] += dt * vel[a][i];
            COMPUTE_ACC();
            for (int a = 0; a < 3; a++)
                for (int i = 1; i < Nx - 1; i++)
                    vel[a][i] += 0.5 * dt * acc[a][i];

            /* Absorbing boundary */
            for (int a = 0; a < 3; a++)
                for (int i = 0; i < Nx; i++) {
                    vel[a][i] *= damp[i];
                    phi[a][i] *= damp[i];
                }
        }
        #undef COMPUTE_ACC

        fclose(fts);

        /* DFT: measure oscillon frequency from second half of run */
        double omega_meas = 0;
        int dft_start = n_dft / 2;
        if (n_dft - dft_start > 100) {
            double T = t_hist[n_dft-1] - t_hist[dft_start];
            int nf = 1000;
            double peak_pow = 0;
            for (int kk = 1; kk < nf; kk++) {
                double omega = 3.0 * mass * kk / nf;
                double re = 0, im = 0;
                for (int j = dft_start; j < n_dft; j++) {
                    double dtj = (j > dft_start) ?
                        (t_hist[j]-t_hist[j-1]) : (t_hist[dft_start+1]-t_hist[dft_start]);
                    re += phi0_hist[j] * cos(omega * t_hist[j]) * dtj;
                    im += phi0_hist[j] * sin(omega * t_hist[j]) * dtj;
                }
                double pw = (re*re + im*im) / (T*T);
                if (pw > peak_pow) { peak_pow = pw; omega_meas = omega; }
            }
        }

        double dE_frac = (E_init > 1e-20) ? (E_final - E_init) / E_init : 0.0;
        int stable = (fc_final > 0.99) && (vac_max < 0.01);

        printf("\n  RESULT: lam=%.6f  fc=%.4f  pk=%.4f  vac=%.2e  dE=%.4f  omega=%.4f  %s\n",
               lam, fc_final, peak_amp, vac_max, dE_frac, omega_meas,
               stable ? "STABLE" : "UNSTABLE");

        fprintf(fsum, "%.6f\t%.6e\t%.6f\t%.1f\t%.4f\t%.4f\t%.2e\t%.6e\t%.6e\t%.6f\t%.4f\t%d\n",
                lam, m2_anti, m_A, pred_range, fc_final, peak_amp, vac_max,
                E_final, E_init, dE_frac, omega_meas, stable);
        fflush(fsum);

        if (stable && lam > lam_best) {
            lam_best = lam;
            fc_best = fc_final;
        }

        for (int a = 0; a < 3; a++) { free(phi[a]); free(vel[a]); free(acc[a]); }
        free(damp); free(phi0_hist); free(t_hist);
    }

    fclose(fsum);

    printf("\n=== PHASE 1 SUMMARY ===\n");
    if (lam_best > 0) {
        double m_A_best = sqrt(m2 - lam_best);
        printf("  Best stable lambda = %.6f\n", lam_best);
        printf("  m_A = %.6f\n", m_A_best);
        printf("  Predicted range = %.1f\n", 1.0 / m_A_best);
        printf("  fc = %.4f\n", fc_best);
    } else {
        printf("  No stable lambda found!\n");
    }
    printf("  Output: %s\n", sumpath);
}

/* ======================================================================
 * PHASE 2: Two-oscillon interaction at given lambda and separation
 * ====================================================================== */
static void phase2_two_oscillon(void)
{
    double m2 = mass * mass;
    double m2_anti = m2 - lambda;
    if (m2_anti <= 0) {
        printf("ERROR: tachyonic at lambda=%.6f\n", lambda);
        return;
    }
    double m_A = sqrt(m2_anti);
    double pred_range = 1.0 / m_A;

    printf("=== PHASE 2: Two-Oscillon Interaction ===\n");
    printf("  lambda=%.6f  m_A=%.6f  pred_range=%.1f  separation=%.1f\n",
           lambda, m_A, pred_range, sep);

    double dx  = 2.0 * xmax / (Nx - 1);
    double dx2 = dx * dx;
    double m2_sym = m2 + 2.0 * lambda;
    double kmax_grid = M_PI / dx;
    double dt = 0.8 * 2.0 / sqrt(kmax_grid * kmax_grid + m2_sym);
    int Nt = (int)(tfinal / dt) + 1;

    printf("  Nx=%d xmax=%.1f dx=%.6f dt=%.6f Nt=%d\n", Nx, xmax, dx, dt, Nt);

    double *phi[3], *vel[3], *acc[3];
    for (int a = 0; a < 3; a++) {
        phi[a] = calloc(Nx, sizeof(double));
        vel[a] = calloc(Nx, sizeof(double));
        acc[a] = calloc(Nx, sizeof(double));
    }

    /* Absorbing boundary */
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

    /* Initialize: two oscillons at x = +/- sep/2 */
    double x1 = -sep / 2.0;
    double x2 = +sep / 2.0;
    for (int a = 0; a < 3; a++)
        for (int i = 0; i < Nx; i++) {
            double x = -xmax + i * dx;
            double g1 = A_init * exp(-(x-x1)*(x-x1) / (2.0 * sigma_w * sigma_w));
            double g2 = A_init * exp(-(x-x2)*(x-x2) / (2.0 * sigma_w * sigma_w));
            phi[a][i] = g1 + g2;
            vel[a][i] = 0.0;
        }

    #define COMPUTE_ACC2() do { \
        for (int a = 0; a < 3; a++) { \
            acc[a][0] = acc[a][1] = acc[a][Nx-2] = acc[a][Nx-1] = 0; \
            int b = (a+1)%3, c_idx = (a+2)%3; \
            for (int i = 1; i < Nx - 1; i++) { \
                double lapl = (phi[a][i+1] - 2.0*phi[a][i] + phi[a][i-1]) / dx2; \
                double fp = force_triple(phi[0][i], phi[1][i], phi[2][i], a); \
                acc[a][i] = lapl - m2*phi[a][i] - lambda*(phi[b][i]+phi[c_idx][i]) + fp; \
            } \
        } \
    } while(0)

    COMPUTE_ACC2();

    /* Output */
    char tspath[600];
    snprintf(tspath, sizeof(tspath), "%s/phase2_sep%.0f.tsv", outdir, sep);
    FILE *fts = fopen(tspath, "w");
    fprintf(fts, "time\tx_left\tx_right\tsep\tE_left\tE_right\tE_mid\tfc_left\tfc_right\n");

    int rec_every = Nt / 10000;
    if (rec_every < 1) rec_every = 1;
    int print_every = Nt / 20;
    if (print_every < 1) print_every = 1;

    double core_r = 3.0 * sigma_w;
    double sep_init = sep;

    for (int n = 0; n <= Nt; n++) {
        double t = n * dt;

        int do_rec = (n % rec_every == 0);
        int do_print = (n % print_every == 0);

        if (do_rec || do_print || n == 0 || n == Nt) {
            /* Find center of energy for left (x<0) and right (x>0) halves */
            double E_L = 0, E_R = 0, E_M = 0;
            double xE_L = 0, xE_R = 0;
            (void)0;

            for (int i = 1; i < Nx - 1; i++) {
                double x = -xmax + i * dx;
                double e = energy_density(phi, vel, i, m2, lambda, dx);
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

            double xcL = (E_L > 1e-20) ? xE_L / E_L : x1;
            double xcR = (E_R > 1e-20) ? xE_R / E_R : x2;
            double d = xcR - xcL;

            /* fc for each oscillon: fraction within core_r of their center */
            /* (approximate: use initial centers) */
            double fcL = 0, fcR = 0, EcL = 0, EcR = 0;
            for (int i = 1; i < Nx - 1; i++) {
                double x = -xmax + i * dx;
                double e = energy_density(phi, vel, i, m2, lambda, dx) * dx;
                if (fabs(x - xcL) < core_r) { EcL += e; }
                if (fabs(x - xcR) < core_r) { EcR += e; }
            }
            fcL = (E_L > 1e-20) ? EcL / E_L : 0;
            fcR = (E_R > 1e-20) ? EcR / E_R : 0;

            if (do_rec)
                fprintf(fts, "%.4f\t%.4f\t%.4f\t%.4f\t%.6e\t%.6e\t%.6e\t%.4f\t%.4f\n",
                        t, xcL, xcR, d, E_L, E_R, E_M, fcL, fcR);

            if (do_print)
                printf("  t=%7.1f  xL=%.2f  xR=%.2f  sep=%.2f  dsep=%.3f  EL=%.3f  ER=%.3f  EM=%.2e\n",
                       t, xcL, xcR, d, d - sep_init, E_L, E_R, E_M);
        }

        if (n == Nt) break;

        /* Velocity Verlet */
        for (int a = 0; a < 3; a++)
            for (int i = 1; i < Nx - 1; i++)
                vel[a][i] += 0.5 * dt * acc[a][i];
        for (int a = 0; a < 3; a++)
            for (int i = 1; i < Nx - 1; i++)
                phi[a][i] += dt * vel[a][i];
        COMPUTE_ACC2();
        for (int a = 0; a < 3; a++)
            for (int i = 1; i < Nx - 1; i++)
                vel[a][i] += 0.5 * dt * acc[a][i];

        /* Absorbing boundary */
        for (int a = 0; a < 3; a++)
            for (int i = 0; i < Nx; i++) {
                vel[a][i] *= damp[i];
                phi[a][i] *= damp[i];
            }
    }
    #undef COMPUTE_ACC2

    fclose(fts);
    printf("\n  Output: %s\n", tspath);

    for (int a = 0; a < 3; a++) { free(phi[a]); free(vel[a]); free(acc[a]); }
    free(damp);
}

/* ====================================================================== */
int main(int argc, char **argv)
{
    parse_args(argc, argv);

    switch (phase) {
        case 1: phase1_scan();        break;
        case 2: phase2_two_oscillon(); break;
        default:
            fprintf(stderr, "Unknown phase %d (use 1 or 2)\n", phase);
            return 1;
    }
    return 0;
}
