/*
 * maxwell_e.c — Proca field from mass-split sector
 *
 * Three massive scalars with triple-product AND pairwise coupling:
 *   L = sum_a [ (1/2)(dt phi_a)^2 - (1/2)(dx phi_a)^2 - (m^2/2)phi_a^2 ]
 *     - (mu/2) P^2 / (1 + kappa P^2)
 *     - lambda(phi_1 phi_2 + phi_2 phi_3 + phi_3 phi_1)
 *
 * Mass spectrum:
 *   antisymmetric mode: m^2_A = m^2 - lambda  (lighter)
 *   symmetric mode:     m^2_S = m^2 + 2*lambda (heavier)
 *
 * Phases:
 *   1: Scan lambda to find lambda_max where oscillon survives (fc > 0.9)
 *   2: At lambda_max, excite antisymmetric perturbation, measure propagation range
 *   3: Two-oscillon interaction via Proca at lambda_max
 *
 * Compile: gcc -O3 -Wall -o maxwell_e src/maxwell_e.c -lm
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
static int    Nx      = 4000;
static double xmax    = 100.0;
static double tfinal  = 5000.0;
static int    phase   = 1;     /* 1=scan, 2=propagation, 3=two-oscillon */
static double sep     = 30.0;  /* separation for phase 3 */
static char   outdir[512] = "v24/maxwell_e/data";

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

/* ---- Single-oscillon evolution returning fc_final ---- */
static double run_single(double lam, double t_end, int quiet,
                         double **phi_out, double **vel_out, int *Nx_out)
{
    double dx  = 2.0 * xmax / (Nx - 1);
    double dx2 = dx * dx;
    double m2  = mass * mass;

    double m2_anti = m2 - lam;
    if (m2_anti <= 0) return -1.0; /* tachyonic */

    /* CFL: based on heaviest mode */
    double m2_sym = m2 + 2.0 * lam;
    double kmax = M_PI / dx;
    double dt = 0.8 * 2.0 / sqrt(kmax * kmax + m2_sym);
    int Nt = (int)(t_end / dt) + 1;

    if (!quiet)
        printf("  lam=%.4f m2_anti=%.4f Nt=%d dt=%.6f\n", lam, m2_anti, Nt, dt);

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

    /* Initialize: symmetric 0-degree oscillon */
    double omega_init = sqrt(m2 + 2.0 * lam) * 0.9; /* use symmetric mass */
    for (int a = 0; a < 3; a++)
        for (int i = 0; i < Nx; i++) {
            double x = -xmax + i * dx;
            phi[a][i] = A_init * exp(-x * x / (2.0 * sigma_w * sigma_w));
            vel[a][i] = 0.0; /* start at max amplitude, zero velocity */
        }

    /* Macro for acceleration */
    #define COMPUTE_ACC_S() do { \
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

    COMPUTE_ACC_S();

    double core_r = 3.0 * sigma_w;
    int print_every = Nt / 20;
    if (print_every < 1) print_every = 1;

    double E_init = 0;
    double fc_final = 0;

    for (int n = 0; n <= Nt; n++) {
        if (n % print_every == 0 || n == Nt) {
            double Eall = 0, Ecore = 0;
            for (int i = 1; i < Nx - 1; i++) {
                double x = -xmax + i * dx;
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
                Eall += e * dx;
                if (fabs(x) < core_r) Ecore += e * dx;
            }
            fc_final = (Eall > 1e-20) ? Ecore / Eall : 0.0;
            if (n == 0) E_init = Eall;

            if (!quiet && n % print_every == 0) {
                double t = n * dt;
                printf("    t=%7.1f  fc=%.4f  E=%.4f  pk=(%.4f,%.4f,%.4f)\n",
                       t, fc_final, Eall,
                       fabs(phi[0][Nx/2]), fabs(phi[1][Nx/2]), fabs(phi[2][Nx/2]));
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
        COMPUTE_ACC_S();
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
    #undef COMPUTE_ACC_S

    /* Return field state if requested */
    if (phi_out && vel_out && Nx_out) {
        *Nx_out = Nx;
        for (int a = 0; a < 3; a++) {
            phi_out[a] = phi[a];
            vel_out[a] = vel[a];
        }
        /* don't free phi/vel */
        for (int a = 0; a < 3; a++) free(acc[a]);
    } else {
        for (int a = 0; a < 3; a++) { free(phi[a]); free(vel[a]); free(acc[a]); }
    }
    free(damp);
    return fc_final;
}

/* ======================================================================
 * PHASE 1: Lambda scan — find max stable lambda
 * ====================================================================== */
static void phase1_scan(void)
{
    printf("=== PHASE 1: Lambda scan for maximum stable coupling ===\n");
    printf("  mu=%.1f kappa=%.1f mass=%.4f A=%.3f sigma=%.3f\n",
           mu, kappa, mass, A_init, sigma_w);
    printf("  Nx=%d xmax=%.1f tfinal=%.0f\n", Nx, xmax, tfinal);

    char path[600];
    snprintf(path, sizeof(path), "%s/phase1_scan.tsv", outdir);
    FILE *fp = fopen(path, "w");
    if (!fp) { fprintf(stderr, "Cannot open %s\n", path); return; }
    fprintf(fp, "lambda\tm2_anti\tm_A\tfc_final\tstable\n");

    double lam_vals[] = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95};
    int n_lam = sizeof(lam_vals) / sizeof(lam_vals[0]);

    double lam_max = 0.0;
    double fc_at_max = 0.0;

    for (int k = 0; k < n_lam; k++) {
        double lam = lam_vals[k];
        double m2 = mass * mass;
        double m2a = m2 - lam;
        if (m2a <= 0) {
            printf("\nlambda=%.4f: TACHYONIC (m2_anti=%.4f), skipping\n", lam, m2a);
            fprintf(fp, "%.4f\t%.4f\t%.4f\t%.4f\t%d\n", lam, m2a, 0.0, 0.0, 0);
            continue;
        }

        printf("\n--- lambda=%.4f  m2_anti=%.4f  m_A=%.4f ---\n", lam, m2a, sqrt(m2a));
        double fc = run_single(lam, tfinal, 0, NULL, NULL, NULL);
        int stable = (fc > 0.9);
        printf("  => fc_final=%.4f  stable=%s\n", fc, stable ? "YES" : "NO");

        fprintf(fp, "%.4f\t%.4f\t%.4f\t%.4f\t%d\n", lam, m2a, sqrt(m2a), fc, stable);
        fflush(fp);

        if (stable && lam > lam_max) {
            lam_max = lam;
            fc_at_max = fc;
        }
    }

    fclose(fp);

    double m2 = mass * mass;
    double m_A_min = (lam_max < m2) ? sqrt(m2 - lam_max) : 0.0;
    double range = (m_A_min > 0) ? 1.0 / m_A_min : 1e30;

    printf("\n=== PHASE 1 RESULTS ===\n");
    printf("  lambda_max = %.4f (fc=%.4f at t=%.0f)\n", lam_max, fc_at_max, tfinal);
    printf("  m_A_min    = %.6f (lightest Proca mass)\n", m_A_min);
    printf("  range_pred = %.4f (1/m_A)\n", range);
    printf("  Output: %s\n", path);
}

/* ======================================================================
 * PHASE 2: Antisymmetric mode propagation at given lambda
 * ====================================================================== */
static void phase2_propagation(void)
{
    double m2 = mass * mass;
    double m2_anti = m2 - lambda;
    double m_A = sqrt(fabs(m2_anti));
    double pred_range = (m_A > 1e-10) ? 1.0 / m_A : 1e30;

    printf("=== PHASE 2: Antisymmetric mode propagation ===\n");
    printf("  lambda=%.4f m_A=%.6f predicted_range=%.4f\n", lambda, m_A, pred_range);

    if (m2_anti <= 0) {
        printf("  ERROR: tachyonic!\n");
        return;
    }

    double dx  = 2.0 * xmax / (Nx - 1);
    double dx2 = dx * dx;
    double m2_sym = m2 + 2.0 * lambda;

    /* CFL */
    double kmax = M_PI / dx;
    double dt = 0.8 * 2.0 / sqrt(kmax * kmax + m2_sym);
    int Nt = (int)(tfinal / dt) + 1;

    printf("  Nx=%d xmax=%.1f dt=%.6f Nt=%d tfinal=%.0f\n", Nx, xmax, dt, Nt, tfinal);

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

    /* Initialize: symmetric oscillon + antisymmetric perturbation */
    /* Base: all three equal (symmetric oscillon) */
    /* Perturbation: A_1 = (phi_1 - phi_2)/sqrt(2) excited as a localized pulse */
    double eps = 0.05;  /* perturbation amplitude (small) */
    double sig_pert = 2.0; /* perturbation width */
    for (int a = 0; a < 3; a++)
        for (int i = 0; i < Nx; i++) {
            double x = -xmax + i * dx;
            double g = A_init * exp(-x * x / (2.0 * sigma_w * sigma_w));
            phi[a][i] = g;
        }
    /* Add antisymmetric perturbation: phi_1 += eps*f, phi_2 -= eps*f, phi_3 unchanged */
    for (int i = 0; i < Nx; i++) {
        double x = -xmax + i * dx;
        double pert = eps * exp(-x * x / (2.0 * sig_pert * sig_pert));
        phi[0][i] += pert;
        phi[1][i] -= pert;
        /* phi[2] unchanged — this excites A_1 = (phi_1-phi_2)/sqrt(2) */
    }

    #define COMPUTE_ACC_P() do { \
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

    COMPUTE_ACC_P();

    /* Output: measure antisymmetric amplitude A(x) = (phi_1 - phi_2)/sqrt(2) at snapshots */
    char path[600];
    snprintf(path, sizeof(path), "%s/phase2_propagation.tsv", outdir);
    FILE *fp = fopen(path, "w");
    fprintf(fp, "time\tx\tA_anti\tphi_sym\n");

    char tspath[600];
    snprintf(tspath, sizeof(tspath), "%s/phase2_ts.tsv", outdir);
    FILE *fts = fopen(tspath, "w");
    fprintf(fts, "time\tA_anti_peak\tA_anti_rms\tprop_dist\tfc\n");

    int snap_every = Nt / 20;
    if (snap_every < 1) snap_every = 1;
    int rec_every = Nt / 10000;
    if (rec_every < 1) rec_every = 1;
    int print_every = Nt / 20;
    if (print_every < 1) print_every = 1;
    double core_r = 3.0 * sigma_w;

    for (int n = 0; n <= Nt; n++) {
        double t = n * dt;

        int do_snap = (n % snap_every == 0);
        int do_rec  = (n % rec_every == 0);
        int do_print = (n % print_every == 0);

        if (do_snap || do_rec || do_print) {
            /* Compute antisymmetric amplitude profile */
            double A_peak = 0, A_rms = 0, prop_dist = 0;
            double Eall = 0, Ecore = 0;
            int n_rms = 0;

            for (int i = 1; i < Nx - 1; i++) {
                double x = -xmax + i * dx;
                double A = (phi[0][i] - phi[1][i]) / sqrt(2.0);
                double sym = (phi[0][i] + phi[1][i] + phi[2][i]) / sqrt(3.0);

                if (do_snap)
                    fprintf(fp, "%.2f\t%.4f\t%.6e\t%.6e\n", t, x, A, sym);

                double absA = fabs(A);
                if (absA > A_peak) A_peak = absA;
                A_rms += A * A * dx;
                /* propagation distance: furthest x where |A| > threshold */
                double thresh = eps * 0.01; /* 1% of initial perturbation */
                if (absA > thresh && fabs(x) > prop_dist)
                    prop_dist = fabs(x);

                /* energy for fc */
                double e = 0;
                for (int a = 0; a < 3; a++) {
                    e += 0.5 * vel[a][i] * vel[a][i];
                    double dp = (phi[a][i+1] - phi[a][i-1]) / (2.0 * dx);
                    e += 0.5 * dp * dp + 0.5 * m2 * phi[a][i] * phi[a][i];
                }
                e += lambda * (phi[0][i]*phi[1][i] + phi[1][i]*phi[2][i]
                              + phi[2][i]*phi[0][i]);
                double P = phi[0][i] * phi[1][i] * phi[2][i];
                e += 0.5 * mu * P * P / (1.0 + kappa * P * P);
                Eall += e * dx;
                if (fabs(x) < core_r) Ecore += e * dx;
            }
            A_rms = sqrt(A_rms / (2.0 * xmax));
            double fc = (Eall > 1e-20) ? Ecore / Eall : 0.0;

            if (do_rec)
                fprintf(fts, "%.4f\t%.6e\t%.6e\t%.4f\t%.4f\n",
                        t, A_peak, A_rms, prop_dist, fc);

            if (do_print)
                printf("  t=%7.1f  A_peak=%.4e  A_rms=%.4e  prop=%.2f  fc=%.4f\n",
                       t, A_peak, A_rms, prop_dist, fc);
        }

        if (n == Nt) break;

        /* Velocity Verlet */
        for (int a = 0; a < 3; a++)
            for (int i = 1; i < Nx - 1; i++)
                vel[a][i] += 0.5 * dt * acc[a][i];
        for (int a = 0; a < 3; a++)
            for (int i = 1; i < Nx - 1; i++)
                phi[a][i] += dt * vel[a][i];
        COMPUTE_ACC_P();
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
    #undef COMPUTE_ACC_P

    fclose(fp);
    fclose(fts);

    printf("\n  Output: %s\n", path);
    printf("  Output: %s\n", tspath);

    for (int a = 0; a < 3; a++) { free(phi[a]); free(vel[a]); free(acc[a]); }
    free(damp);
}

/* ======================================================================
 * PHASE 3: Two-oscillon interaction via Proca exchange
 * ====================================================================== */
static void phase3_two_oscillon(void)
{
    double m2 = mass * mass;
    double m2_anti = m2 - lambda;
    double m_A = sqrt(fabs(m2_anti));
    double pred_range = (m_A > 1e-10) ? 1.0 / m_A : 1e30;

    printf("=== PHASE 3: Two-oscillon interaction at lambda=%.4f ===\n", lambda);
    printf("  m_A=%.6f predicted_range=%.4f separation=%.1f\n", m_A, pred_range, sep);

    if (m2_anti <= 0) {
        printf("  ERROR: tachyonic!\n");
        return;
    }

    double dx  = 2.0 * xmax / (Nx - 1);
    double dx2 = dx * dx;
    double m2_sym = m2 + 2.0 * lambda;

    double kmax = M_PI / dx;
    double dt = 0.8 * 2.0 / sqrt(kmax * kmax + m2_sym);
    int Nt = (int)(tfinal / dt) + 1;

    printf("  Nx=%d xmax=%.1f dt=%.6f Nt=%d tfinal=%.0f\n", Nx, xmax, dt, Nt, tfinal);

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
        }

    #define COMPUTE_ACC_T() do { \
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

    COMPUTE_ACC_T();

    /* Track center of energy of each oscillon */
    char tspath[600];
    snprintf(tspath, sizeof(tspath), "%s/phase3_two_osc.tsv", outdir);
    FILE *fts = fopen(tspath, "w");
    fprintf(fts, "time\tx_left\tx_right\tsep\tE_left\tE_right\tE_mid\tfc_left\tfc_right\n");

    int rec_every = Nt / 10000;
    if (rec_every < 1) rec_every = 1;
    int print_every = Nt / 20;
    if (print_every < 1) print_every = 1;

    for (int n = 0; n <= Nt; n++) {
        double t = n * dt;

        int do_rec  = (n % rec_every == 0);
        int do_print = (n % print_every == 0);

        if (do_rec || do_print) {
            /* Find center of energy for left (x<0) and right (x>0) halves */
            double E_L = 0, E_R = 0, E_M = 0;
            double xE_L = 0, xE_R = 0;
            double core_L = 0, core_R = 0, tot_L = 0, tot_R = 0;
            double core_r = 3.0 * sigma_w;

            for (int i = 1; i < Nx - 1; i++) {
                double x = -xmax + i * dx;
                double e = 0;
                for (int a = 0; a < 3; a++) {
                    e += 0.5 * vel[a][i] * vel[a][i];
                    double dp = (phi[a][i+1] - phi[a][i-1]) / (2.0 * dx);
                    e += 0.5 * dp * dp + 0.5 * m2 * phi[a][i] * phi[a][i];
                }
                e += lambda * (phi[0][i]*phi[1][i] + phi[1][i]*phi[2][i]
                              + phi[2][i]*phi[0][i]);
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

            double xcL = (E_L > 1e-20) ? xE_L / E_L : x1;
            double xcR = (E_R > 1e-20) ? xE_R / E_R : x2;
            double d = xcR - xcL;

            if (do_rec)
                fprintf(fts, "%.4f\t%.4f\t%.4f\t%.4f\t%.6e\t%.6e\t%.6e\t%.4f\t%.4f\n",
                        t, xcL, xcR, d, E_L, E_R, E_M, 0.0, 0.0);

            if (do_print)
                printf("  t=%7.1f  xL=%.2f  xR=%.2f  sep=%.2f  EL=%.4f  ER=%.4f  EM=%.4e\n",
                       t, xcL, xcR, d, E_L, E_R, E_M);
        }

        if (n == Nt) break;

        /* Velocity Verlet */
        for (int a = 0; a < 3; a++)
            for (int i = 1; i < Nx - 1; i++)
                vel[a][i] += 0.5 * dt * acc[a][i];
        for (int a = 0; a < 3; a++)
            for (int i = 1; i < Nx - 1; i++)
                phi[a][i] += dt * vel[a][i];
        COMPUTE_ACC_T();
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
    #undef COMPUTE_ACC_T

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
        case 2: phase2_propagation(); break;
        case 3: phase3_two_oscillon(); break;
        default:
            fprintf(stderr, "Unknown phase %d (use 1, 2, or 3)\n", phase);
            return 1;
    }
    return 0;
}
