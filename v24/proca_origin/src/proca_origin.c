/*
 * proca_origin.c — V24-S5: What Determines lambda?
 *
 * Three tests for the origin of the pairwise coupling lambda:
 *   A: Gap margin universality — is delta = 1 - omega/m_S constant across lambda?
 *   B: Energy minimum — is there a minimum in E_osc(lambda)?
 *   C: Dynamical modulus — promote lambda to a field Lambda(x,t)
 *
 * Lagrangian (Tests A,B):
 *   L = sum_a [ (1/2)(dt phi_a)^2 - (1/2)(dx phi_a)^2 - (m^2/2)phi_a^2 ]
 *     - lambda(phi_1 phi_2 + phi_2 phi_3 + phi_3 phi_1)
 *     - (mu/2) P^2 / (1 + kappa P^2),   P = phi_1 phi_2 phi_3
 *
 * Test C adds: L_Lambda = (1/2)(dt Lambda)^2 - (1/2)(dx Lambda)^2
 *              - (m_Lambda^2/2)(Lambda - Lambda_0)^2
 *   with pairwise coupling: Lambda(x)(phi_1 phi_2 + phi_2 phi_3 + phi_3 phi_1)
 *   and source: -g_Lambda * (phi_1 phi_2 + phi_2 phi_3 + phi_3 phi_1)
 *
 * Compile: gcc -O3 -Wall -o proca_origin src/proca_origin.c -lm
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ---- Parameters ---- */
static double mu      = -20.0;
static double kappa   = 20.0;
static double mass    = 1.0;
static double A_init  = 0.8;
static double sigma_w = 3.0;
static int    Nx      = 4000;
static double xmax    = 100.0;
static double tfinal  = 10000.0;
static int    test_sel = 0;  /* 0=all, 1=A, 2=B, 3=C */
static char   outdir[512] = "v24/proca_origin/data";

/* Test C parameters */
static double Lambda_0 = 0.5;
static double m_Lambda = 1.0;
static double g_Lambda = 0.01;
static double gamma_L  = 0.01;  /* friction on Lambda for settling */

static void parse_args(int argc, char **argv)
{
    for (int i = 1; i < argc - 1; i += 2) {
        if      (!strcmp(argv[i], "-mu"))       mu       = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-kappa"))    kappa    = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-mass"))     mass     = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-A"))        A_init   = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-sigma"))    sigma_w  = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-Nx"))       Nx       = atoi(argv[i+1]);
        else if (!strcmp(argv[i], "-xmax"))     xmax     = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-tfinal"))   tfinal   = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-test"))     test_sel = atoi(argv[i+1]);
        else if (!strcmp(argv[i], "-Lambda0"))  Lambda_0 = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-mLambda"))  m_Lambda = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-gLambda"))  g_Lambda = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-gammaL"))  gamma_L  = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-o"))        strncpy(outdir, argv[i+1], sizeof(outdir)-1);
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

/* ======================================================================
 * Core oscillon solver for a given lambda.  Returns measured omega and
 * total oscillon energy (core fraction > 0.5).
 * ====================================================================== */
static void run_oscillon(double lam, double t_run,
                         double *omega_out, double *energy_out,
                         double *peak_out, double *fc_out,
                         const char *tsfile)
{
    double m2  = mass * mass;
    double dx  = 2.0 * xmax / (Nx - 1);
    double dx2 = dx * dx;

    /* CFL based on heaviest mode */
    double m2_eff = m2 + 2.0 * fabs(lam);
    double kmax_grid = M_PI / dx;
    double dt = 0.8 * 2.0 / sqrt(kmax_grid * kmax_grid + m2_eff);
    int Nt = (int)(t_run / dt) + 1;

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

    /* Initialize: symmetric Gaussians */
    for (int a = 0; a < 3; a++)
        for (int i = 0; i < Nx; i++) {
            double x = -xmax + i * dx;
            phi[a][i] = A_init * exp(-x * x / (2.0 * sigma_w * sigma_w));
        }

    int ic = Nx / 2;
    double core_r = 3.0 * sigma_w;

    /* DFT storage */
    int max_dft = 50000;
    double *phi0_hist = malloc(max_dft * sizeof(double));
    double *t_hist    = malloc(max_dft * sizeof(double));
    int n_dft = 0;
    int dft_every = Nt / max_dft;
    if (dft_every < 1) dft_every = 1;

    /* Time series file (optional) */
    FILE *fts = NULL;
    if (tsfile) {
        fts = fopen(tsfile, "w");
        if (fts) fprintf(fts, "time\tphi0\tpeak\tE_total\tfc\n");
    }
    int rec_every = Nt / 2000;
    if (rec_every < 1) rec_every = 1;
    int print_every = Nt / 10;
    if (print_every < 1) print_every = 1;

    /* Acceleration macro */
    #define COMPUTE_ACC_AB() do { \
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

    COMPUTE_ACC_AB();

    double E_final = 0, fc_final = 0, pk_final = 0;

    for (int n = 0; n <= Nt; n++) {
        double t = n * dt;

        /* DFT sampling */
        if (n % dft_every == 0 && n_dft < max_dft) {
            phi0_hist[n_dft] = phi[0][ic];
            t_hist[n_dft] = t;
            n_dft++;
        }

        int do_rec   = (fts && n % rec_every == 0);
        int do_print = (n % print_every == 0);

        if (do_rec || do_print || n == Nt) {
            double Eall = 0, Ecore = 0, pk = 0;
            for (int i = 1; i < Nx - 1; i++) {
                double x = -xmax + i * dx;
                double e = 0;
                for (int a = 0; a < 3; a++) {
                    e += 0.5 * vel[a][i] * vel[a][i];
                    double dp = (phi[a][i+1] - phi[a][i-1]) / (2.0*dx);
                    e += 0.5 * dp * dp + 0.5 * m2 * phi[a][i] * phi[a][i];
                    if (fabs(phi[a][i]) > pk) pk = fabs(phi[a][i]);
                }
                e += lam * (phi[0][i]*phi[1][i] + phi[1][i]*phi[2][i]
                           + phi[2][i]*phi[0][i]);
                double P = phi[0][i] * phi[1][i] * phi[2][i];
                e += 0.5 * mu * P * P / (1.0 + kappa * P * P);
                Eall += e * dx;
                if (fabs(x) < core_r) Ecore += e * dx;
            }
            double fc = (Eall > 1e-20) ? Ecore / Eall : 0.0;
            E_final  = Eall;
            fc_final = fc;
            pk_final = pk;

            if (do_rec)
                fprintf(fts, "%.4f\t%.6e\t%.6e\t%.6e\t%.4f\n",
                        t, phi[0][ic], pk, Eall, fc);
            if (do_print)
                printf("    t=%7.1f  phi0=%+.4f  pk=%.4f  E=%.4f  fc=%.4f\n",
                       t, phi[0][ic], pk, Eall, fc);
        }

        if (n == Nt) break;

        /* Velocity Verlet */
        for (int a = 0; a < 3; a++)
            for (int i = 1; i < Nx - 1; i++)
                vel[a][i] += 0.5 * dt * acc[a][i];
        for (int a = 0; a < 3; a++)
            for (int i = 1; i < Nx - 1; i++)
                phi[a][i] += dt * vel[a][i];
        COMPUTE_ACC_AB();
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
    #undef COMPUTE_ACC_AB

    if (fts) fclose(fts);

    /* DFT: measure frequency from second half */
    double omega_meas = 0;
    int dft_start = n_dft / 2;
    if (n_dft - dft_start > 100) {
        double T = t_hist[n_dft-1] - t_hist[dft_start];
        int nf = 2000;
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

    *omega_out  = omega_meas;
    *energy_out = E_final;
    *peak_out   = pk_final;
    *fc_out     = fc_final;

    for (int a = 0; a < 3; a++) { free(phi[a]); free(vel[a]); free(acc[a]); }
    free(damp); free(phi0_hist); free(t_hist);
}

/* ======================================================================
 * TEST A: Gap Margin Universality
 * ====================================================================== */
static void test_A(void)
{
    double lam_vals[] = {0.0, 0.1, 0.3, 0.5, 0.7, 0.9, 0.95, 0.99};
    int n_lam = sizeof(lam_vals) / sizeof(lam_vals[0]);
    double m2 = mass * mass;

    printf("\n");
    printf("============================================================\n");
    printf("  TEST A: Gap Margin Universality\n");
    printf("============================================================\n");
    printf("  mu=%.1f kappa=%.1f mass=%.4f A=%.3f sigma=%.3f\n",
           mu, kappa, mass, A_init, sigma_w);
    printf("  Nx=%d xmax=%.1f tfinal=%.0f\n", Nx, xmax, tfinal);

    char sumpath[600];
    snprintf(sumpath, sizeof(sumpath), "%s/testA_gap_margin.tsv", outdir);
    FILE *fsum = fopen(sumpath, "w");
    if (!fsum) { fprintf(stderr, "Cannot open %s\n", sumpath); return; }
    fprintf(fsum, "lambda\tm_S\tomega\tgap_margin\tE_osc\tE_over_mS\tpeak\tfc\n");

    printf("\n  %-8s  %-8s  %-8s  %-12s  %-10s  %-10s  %-8s  %-6s\n",
           "lambda", "m_S", "omega", "gap_margin", "E_osc", "E/m_S", "peak", "fc");
    printf("  %-8s  %-8s  %-8s  %-12s  %-10s  %-10s  %-8s  %-6s\n",
           "------", "----", "-----", "----------", "-----", "-----", "----", "--");

    for (int k = 0; k < n_lam; k++) {
        double lam = lam_vals[k];
        double m_S = sqrt(m2 + 2.0 * lam);

        printf("\n  --- lambda = %.2f, m_S = %.4f ---\n", lam, m_S);

        char tspath[600];
        snprintf(tspath, sizeof(tspath), "%s/testA_lam%.2f_ts.tsv", outdir, lam);

        double omega, energy, peak, fc;
        run_oscillon(lam, tfinal, &omega, &energy, &peak, &fc, tspath);

        double gap_margin = (m_S > 1e-10) ? 1.0 - omega / m_S : -1.0;
        double E_over_mS  = (m_S > 1e-10) ? energy / m_S : 0.0;

        printf("  RESULT: lam=%.2f  m_S=%.4f  omega=%.4f  delta=%.6f  E=%.4f  E/mS=%.4f  pk=%.4f  fc=%.4f\n",
               lam, m_S, omega, gap_margin, energy, E_over_mS, peak, fc);

        fprintf(fsum, "%.4f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.4f\n",
                lam, m_S, omega, gap_margin, energy, E_over_mS, peak, fc);
        fflush(fsum);
    }

    fclose(fsum);
    printf("\n  Test A output: %s\n", sumpath);
}

/* ======================================================================
 * TEST B: Energy Minimum (uses same data as Test A)
 * ====================================================================== */
static void test_B(void)
{
    printf("\n");
    printf("============================================================\n");
    printf("  TEST B: Energy Minimum (combined with Test A data)\n");
    printf("============================================================\n");
    printf("  (Test B uses the same runs as Test A.  The energy values\n");
    printf("   are already in testA_gap_margin.tsv.)\n");
    printf("  Columns E_osc and E_over_mS in that file contain the data.\n");
    printf("  See RESULTS.md for analysis.\n");
}

/* ======================================================================
 * TEST C: Dynamical Modulus
 * ====================================================================== */
static void test_C(void)
{
    double m2 = mass * mass;
    double dx  = 2.0 * xmax / (Nx - 1);
    double dx2 = dx * dx;

    /* CFL: include modulus mass and maximum effective pairwise coupling */
    double m2_eff = m2 + 2.0 * fabs(Lambda_0) + 2.0;  /* conservative */
    double kmax_grid = M_PI / dx;
    double dt = 0.8 * 2.0 / sqrt(kmax_grid * kmax_grid + m2_eff);
    int Nt = (int)(tfinal / dt) + 1;

    printf("\n");
    printf("============================================================\n");
    printf("  TEST C: Dynamical Modulus\n");
    printf("============================================================\n");
    printf("  mu=%.1f kappa=%.1f mass=%.4f A=%.3f sigma=%.3f\n",
           mu, kappa, mass, A_init, sigma_w);
    printf("  Lambda_0=%.3f  m_Lambda=%.4f  g_Lambda=%.4f  gamma_L=%.4f\n",
           Lambda_0, m_Lambda, g_Lambda, gamma_L);
    printf("  Nx=%d xmax=%.1f dx=%.6f dt=%.6f Nt=%d tfinal=%.0f\n",
           Nx, xmax, dx, dt, Nt, tfinal);

    /* Allocate fields: phi[3], vel[3], acc[3] + Lambda, vel_L, acc_L */
    double *phi[3], *vel[3], *acc[3];
    for (int a = 0; a < 3; a++) {
        phi[a] = calloc(Nx, sizeof(double));
        vel[a] = calloc(Nx, sizeof(double));
        acc[a] = calloc(Nx, sizeof(double));
    }
    double *Lam     = calloc(Nx, sizeof(double));
    double *vel_L   = calloc(Nx, sizeof(double));
    double *acc_L   = calloc(Nx, sizeof(double));

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

    /* Initialize: symmetric Gaussians for phi, uniform Lambda_0 */
    for (int a = 0; a < 3; a++)
        for (int i = 0; i < Nx; i++) {
            double x = -xmax + i * dx;
            phi[a][i] = A_init * exp(-x * x / (2.0 * sigma_w * sigma_w));
        }
    for (int i = 0; i < Nx; i++)
        Lam[i] = Lambda_0;

    int ic = Nx / 2;
    double core_r = 3.0 * sigma_w;

    /* Acceleration macro for Test C */
    #define COMPUTE_ACC_C() do { \
        for (int a = 0; a < 3; a++) { \
            acc[a][0] = acc[a][1] = acc[a][Nx-2] = acc[a][Nx-1] = 0; \
            int b = (a+1)%3, c_idx = (a+2)%3; \
            for (int i = 1; i < Nx - 1; i++) { \
                double lapl = (phi[a][i+1] - 2.0*phi[a][i] + phi[a][i-1]) / dx2; \
                double fp = force_triple(phi[0][i], phi[1][i], phi[2][i], a); \
                acc[a][i] = lapl - m2*phi[a][i] \
                    - Lam[i]*(phi[b][i]+phi[c_idx][i]) + fp; \
            } \
        } \
        /* Lambda acceleration */ \
        acc_L[0] = acc_L[1] = acc_L[Nx-2] = acc_L[Nx-1] = 0; \
        for (int i = 1; i < Nx - 1; i++) { \
            double lapl_L = (Lam[i+1] - 2.0*Lam[i] + Lam[i-1]) / dx2; \
            double pw = phi[0][i]*phi[1][i] + phi[1][i]*phi[2][i] \
                       + phi[2][i]*phi[0][i]; \
            acc_L[i] = lapl_L - m_Lambda*m_Lambda*(Lam[i] - Lambda_0) \
                       - g_Lambda * pw; \
        } \
    } while(0)

    COMPUTE_ACC_C();

    /* Output */
    char tspath[600];
    snprintf(tspath, sizeof(tspath), "%s/testC_modulus_ts.tsv", outdir);
    FILE *fts = fopen(tspath, "w");
    fprintf(fts, "time\tLam_center\tLam_edge\tphi0\tpeak\tE_total\tfc\tLam_min\tLam_max\n");

    char profpath[600];
    snprintf(profpath, sizeof(profpath), "%s/testC_modulus_profile.tsv", outdir);

    int rec_every = Nt / 5000;
    if (rec_every < 1) rec_every = 1;
    int print_every = Nt / 20;
    if (print_every < 1) print_every = 1;

    /* Profile snapshots */
    int prof_snap = 0;
    double prof_times[] = {0, 100, 500, 1000, 2000, 5000, 10000};
    int n_prof_times = sizeof(prof_times) / sizeof(prof_times[0]);

    FILE *fprof = fopen(profpath, "w");
    if (fprof) fprintf(fprof, "time\tx\tphi0\tLambda\n");

    for (int n = 0; n <= Nt; n++) {
        double t = n * dt;

        /* Profile snapshot */
        if (prof_snap < n_prof_times && t >= prof_times[prof_snap] - 0.5*dt) {
            if (fprof) {
                int stride = Nx / 500;
                if (stride < 1) stride = 1;
                for (int i = 0; i < Nx; i += stride)
                    fprintf(fprof, "%.2f\t%.4f\t%.6e\t%.6e\n",
                            t, -xmax + i*dx, phi[0][i], Lam[i]);
            }
            prof_snap++;
        }

        int do_rec   = (n % rec_every == 0);
        int do_print = (n % print_every == 0);

        if (do_rec || do_print || n == Nt) {
            double Eall = 0, Ecore = 0, pk = 0;
            double Lam_min = 1e30, Lam_max = -1e30;
            int i_edge = (int)((xmax * 0.5 + xmax) / dx);  /* x = xmax/2 */

            for (int i = 1; i < Nx - 1; i++) {
                /* phi energy */
                double e = 0;
                for (int a = 0; a < 3; a++) {
                    e += 0.5 * vel[a][i] * vel[a][i];
                    double dp = (phi[a][i+1] - phi[a][i-1]) / (2.0*dx);
                    e += 0.5 * dp * dp + 0.5 * m2 * phi[a][i] * phi[a][i];
                    if (fabs(phi[a][i]) > pk) pk = fabs(phi[a][i]);
                }
                double pw = phi[0][i]*phi[1][i] + phi[1][i]*phi[2][i]
                           + phi[2][i]*phi[0][i];
                e += Lam[i] * pw;
                double P = phi[0][i] * phi[1][i] * phi[2][i];
                e += 0.5 * mu * P * P / (1.0 + kappa * P * P);

                /* Lambda energy */
                e += 0.5 * vel_L[i] * vel_L[i];
                double dL = (Lam[i+1] - Lam[i-1]) / (2.0*dx);
                e += 0.5 * dL * dL;
                e += 0.5 * m_Lambda * m_Lambda * (Lam[i] - Lambda_0) * (Lam[i] - Lambda_0);

                double x = -xmax + i * dx;
                Eall += e * dx;
                if (fabs(x) < core_r) Ecore += e * dx;

                if (Lam[i] < Lam_min) Lam_min = Lam[i];
                if (Lam[i] > Lam_max) Lam_max = Lam[i];
            }

            double fc = (Eall > 1e-20) ? Ecore / Eall : 0.0;

            if (do_rec)
                fprintf(fts, "%.4f\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.4f\t%.6e\t%.6e\n",
                        t, Lam[ic], Lam[i_edge < Nx ? i_edge : Nx-2],
                        phi[0][ic], pk, Eall, fc, Lam_min, Lam_max);

            if (do_print)
                printf("  t=%7.1f  Lam(0)=%.5f  Lam(edge)=%.5f  phi0=%+.4f  pk=%.4f  E=%.4f  fc=%.4f\n",
                       t, Lam[ic], Lam[i_edge < Nx ? i_edge : Nx-2],
                       phi[0][ic], pk, Eall, fc);
        }

        if (n == Nt) break;

        /* Velocity Verlet: half-step velocity */
        for (int a = 0; a < 3; a++)
            for (int i = 1; i < Nx - 1; i++)
                vel[a][i] += 0.5 * dt * acc[a][i];
        for (int i = 1; i < Nx - 1; i++)
            vel_L[i] += 0.5 * dt * acc_L[i];

        /* Full-step position */
        for (int a = 0; a < 3; a++)
            for (int i = 1; i < Nx - 1; i++)
                phi[a][i] += dt * vel[a][i];
        for (int i = 1; i < Nx - 1; i++)
            Lam[i] += dt * vel_L[i];

        /* Recompute acceleration */
        COMPUTE_ACC_C();

        /* Half-step velocity */
        for (int a = 0; a < 3; a++)
            for (int i = 1; i < Nx - 1; i++)
                vel[a][i] += 0.5 * dt * acc[a][i];
        for (int i = 1; i < Nx - 1; i++)
            vel_L[i] += 0.5 * dt * acc_L[i];

        /* Friction on Lambda to help it settle */
        double fric = exp(-gamma_L * dt);
        for (int i = 0; i < Nx; i++)
            vel_L[i] *= fric;

        /* Absorbing boundary */
        for (int a = 0; a < 3; a++)
            for (int i = 0; i < Nx; i++) {
                vel[a][i] *= damp[i];
                phi[a][i] *= damp[i];
            }
        for (int i = 0; i < Nx; i++) {
            vel_L[i] *= damp[i];
            /* Lambda should relax to Lambda_0 at boundaries, not zero */
            Lam[i] = Lambda_0 + (Lam[i] - Lambda_0) * damp[i];
        }
    }
    #undef COMPUTE_ACC_C

    fclose(fts);
    if (fprof) fclose(fprof);

    printf("\n  Lambda(center) = %.6f  (vacuum = %.6f, shift = %.6f)\n",
           Lam[ic], Lambda_0, Lam[ic] - Lambda_0);
    printf("  Test C output: %s\n", tspath);
    printf("  Profiles: %s\n", profpath);

    for (int a = 0; a < 3; a++) { free(phi[a]); free(vel[a]); free(acc[a]); }
    free(Lam); free(vel_L); free(acc_L); free(damp);
}

/* ====================================================================== */
int main(int argc, char **argv)
{
    parse_args(argc, argv);

    printf("proca_origin: V24-S5 — What Determines Lambda?\n");
    printf("  mu=%.1f kappa=%.1f mass=%.4f\n", mu, kappa, mass);

    if (test_sel == 0 || test_sel == 1) test_A();
    if (test_sel == 0 || test_sel == 2) test_B();
    if (test_sel == 0 || test_sel == 3) test_C();

    printf("\nDone.\n");
    return 0;
}
