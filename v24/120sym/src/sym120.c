/*
 * sym120.c — V24-SYM: Symmetric penalty (φ₁+φ₂+φ₃)² forcing 120° phase separation
 *
 * Lagrangian:
 *   L = sum_a [ (1/2)(∂t φ_a)² - (1/2)(∂x φ_a)² - (m²/2)φ_a² ]
 *     - (μ/2) P² / (1 + κ P²)        [triple product, P = φ₁φ₂φ₃]
 *     - λ_Q (φ₁ + φ₂ + φ₃)²          [symmetric penalty]
 *
 * EOM: ∂²φ_a/∂t² = ∂²φ_a/∂x² - m²φ_a - 2λ_Q(φ₁+φ₂+φ₃) + f_triple_a
 *
 * Test modes:
 *   1: 0° initial (all three identical Gaussians) — should transition to 120° at large λ_Q
 *   2: 120° initial (Gaussians × cos(ωt + 2πa/3) seeded via velocity)
 *
 * Compile: gcc -O3 -Wall -o sym120 src/sym120.c -lm
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static double mu       = -20.0;
static double kappa    = 20.0;
static double mass     = 1.0;
static double lambda_Q = 0.0;
static double A_init   = 0.8;
static double sigma_w  = 3.0;
static int    Nx       = 4000;
static double xmax     = 100.0;
static double tfinal   = 10000.0;
static int    test     = 1;
static char   outdir[512] = "v24/120sym/data";

static void parse_args(int argc, char **argv)
{
    for (int i = 1; i < argc - 1; i += 2) {
        if      (!strcmp(argv[i], "-mu"))      mu       = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-kappa"))   kappa    = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-mass"))    mass     = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-lQ"))      lambda_Q = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-A"))       A_init   = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-sigma"))   sigma_w  = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-Nx"))      Nx       = atoi(argv[i+1]);
        else if (!strcmp(argv[i], "-xmax"))    xmax     = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-tfinal"))  tfinal   = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-test"))    test     = atoi(argv[i+1]);
        else if (!strcmp(argv[i], "-o"))       strncpy(outdir, argv[i+1], sizeof(outdir)-1);
        else { fprintf(stderr, "Unknown flag: %s\n", argv[i]); exit(1); }
    }
}

/* -dV_triple/dφ_a where V = (μ/2)P²/(1+κP²), P = φ₁φ₂φ₃ */
static double force_triple(double p1, double p2, double p3, int a)
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

int main(int argc, char **argv)
{
    parse_args(argc, argv);

    double dx  = 2.0 * xmax / (Nx - 1);
    double dx2 = dx * dx;
    double m2  = mass * mass;

    /* CFL: max frequency is sqrt(k_max^2 + m^2 + 6*lambda_Q) */
    double kmax = M_PI / dx;
    double omega_max2 = kmax * kmax + m2 + 6.0 * lambda_Q;
    double dt = 0.8 * 2.0 / sqrt(omega_max2);
    int Nt = (int)(tfinal / dt) + 1;

    const char *test_desc;
    switch (test) {
        case 1: test_desc = "0-degree init (all identical)"; break;
        case 2: test_desc = "120-degree init (phase-shifted)"; break;
        default: fprintf(stderr, "Unknown test %d\n", test); return 1;
    }

    printf("sym120 test %d: %s\n", test, test_desc);
    printf("  mu=%.1f kappa=%.1f mass=%.4f lambda_Q=%.4f\n", mu, kappa, mass, lambda_Q);
    printf("  A=%.3f sigma=%.3f\n", A_init, sigma_w);
    printf("  Nx=%d xmax=%.1f dx=%.5f dt=%.6f Nt=%d tfinal=%.0f\n",
           Nx, xmax, dx, dt, Nt, tfinal);

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
    double omega_init = 0.9 * mass;  /* sub-mass frequency for oscillon */
    for (int i = 0; i < Nx; i++) {
        double x = -xmax + i * dx;
        double env = A_init * exp(-x * x / (2.0 * sigma_w * sigma_w));
        if (test == 1) {
            /* 0°: all three identical */
            for (int a = 0; a < 3; a++) {
                phi[a][i] = env;
                vel[a][i] = 0.0;
            }
        } else {
            /* 120°: at t=0, φ_a = env * cos(2πa/3), dφ_a/dt = -ω*env*sin(2πa/3) */
            for (int a = 0; a < 3; a++) {
                double phase = 2.0 * M_PI * a / 3.0;
                phi[a][i] = env * cos(phase);
                vel[a][i] = -omega_init * env * sin(phase);
            }
        }
    }

    /* Compute acceleration */
    #define COMPUTE_ACC() do { \
        for (int a = 0; a < 3; a++) { \
            acc[a][0] = acc[a][1] = acc[a][Nx-2] = acc[a][Nx-1] = 0; \
        } \
        for (int i = 1; i < Nx - 1; i++) { \
            double S = phi[0][i] + phi[1][i] + phi[2][i]; \
            double pen = -2.0 * lambda_Q * S; \
            for (int a = 0; a < 3; a++) { \
                double lapl = (phi[a][i+1] - 2.0*phi[a][i] + phi[a][i-1]) / dx2; \
                double ft = force_triple(phi[0][i], phi[1][i], phi[2][i], a); \
                acc[a][i] = lapl - m2 * phi[a][i] + pen + ft; \
            } \
        } \
    } while(0)

    COMPUTE_ACC();

    /* Time series output */
    char tspath[600];
    snprintf(tspath, sizeof(tspath), "%s/sym120_test%d_lQ%.2f_ts.tsv", outdir, test, lambda_Q);
    FILE *fts = fopen(tspath, "w");
    if (!fts) { perror(tspath); return 1; }
    fprintf(fts, "time\tphi1_0\tphi2_0\tphi3_0\tpeak1\tpeak2\tpeak3\t"
                 "E_kin\tE_grad\tE_mass\tE_pot\tE_pen\tE_total\tf_core\tsum_center\n");

    /* DFT storage — per-field */
    int max_dft = 50000;
    double *phi_hist[3];
    for (int a = 0; a < 3; a++)
        phi_hist[a] = malloc(max_dft * sizeof(double));
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

        if (n % dft_every == 0 && n_dft < max_dft) {
            for (int a = 0; a < 3; a++)
                phi_hist[a][n_dft] = phi[a][ic];
            t_hist[n_dft] = t;
            n_dft++;
        }

        int do_rec = (n % rec_every == 0);
        int do_print = (n % print_every == 0);

        if (do_rec || do_print) {
            double Ek = 0, Eg = 0, Em = 0, Ep = 0, Epen = 0;
            double Ecore = 0, Eall = 0;
            double peak[3] = {0};

            for (int i = 1; i < Nx - 1; i++) {
                double x = -xmax + i * dx;
                for (int a = 0; a < 3; a++) {
                    Ek += 0.5 * vel[a][i] * vel[a][i] * dx;
                    double dp = (phi[a][i+1] - phi[a][i-1]) / (2.0 * dx);
                    Eg += 0.5 * dp * dp * dx;
                    Em += 0.5 * m2 * phi[a][i] * phi[a][i] * dx;
                    if (fabs(phi[a][i]) > peak[a]) peak[a] = fabs(phi[a][i]);
                }

                double P  = phi[0][i] * phi[1][i] * phi[2][i];
                double P2 = P * P;
                double V  = 0.5 * mu * P2 / (1.0 + kappa * P2);
                Ep += V * dx;

                double S  = phi[0][i] + phi[1][i] + phi[2][i];
                double Vp = lambda_Q * S * S;
                Epen += Vp * dx;

                double e = V + Vp;
                for (int a = 0; a < 3; a++) {
                    e += 0.5 * vel[a][i] * vel[a][i];
                    double dp = (phi[a][i+1] - phi[a][i-1]) / (2.0 * dx);
                    e += 0.5 * dp * dp + 0.5 * m2 * phi[a][i] * phi[a][i];
                }
                Eall += e * dx;
                if (fabs(x) < core_r) Ecore += e * dx;
            }

            double Et = Ek + Eg + Em + Ep + Epen;
            double fc = (Eall > 1e-20) ? Ecore / Eall : 0.0;
            double sum_c = phi[0][ic] + phi[1][ic] + phi[2][ic];

            if (do_rec)
                fprintf(fts, "%.6f\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t"
                             "%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6f\t%.6e\n",
                        t, phi[0][ic], phi[1][ic], phi[2][ic],
                        peak[0], peak[1], peak[2],
                        Ek, Eg, Em, Ep, Epen, Et, fc, sum_c);
            if (do_print)
                printf("  t=%7.1f  p0=(%+.3f,%+.3f,%+.3f) S=%+.4f  pk=(%.3f,%.3f,%.3f)  "
                       "E=%+.4f  Ep=%+.4f  Epen=%+.4f  fc=%.3f\n",
                       t, phi[0][ic], phi[1][ic], phi[2][ic], sum_c,
                       peak[0], peak[1], peak[2], Et, Ep, Epen, fc);
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

        /* absorbing boundary */
        for (int a = 0; a < 3; a++)
            for (int i = 0; i < Nx; i++) {
                vel[a][i] *= damp[i];
                phi[a][i] *= damp[i];
            }
    }

    fclose(fts);

    /* DFT of each field at center — second half only */
    int dft_start = n_dft / 2;
    if (n_dft - dft_start > 100) {
        char dftpath[600];
        snprintf(dftpath, sizeof(dftpath), "%s/sym120_test%d_lQ%.2f_spectrum.tsv",
                 outdir, test, lambda_Q);
        FILE *fdft = fopen(dftpath, "w");
        fprintf(fdft, "omega\tpower1\tpower2\tpower3\tpower_sum\n");
        double T = t_hist[n_dft-1] - t_hist[dft_start];
        int nf = 500;
        double peak_pow[3] = {0}, peak_om[3] = {0};
        double peak_sum_pow = 0, peak_sum_om = 0;
        for (int k = 0; k < nf; k++) {
            double omega = 3.0 * mass * k / nf;
            double re[3] = {0}, im[3] = {0};
            double re_s = 0, im_s = 0;
            for (int j = dft_start; j < n_dft; j++) {
                double dtj = (j > dft_start) ?
                    (t_hist[j] - t_hist[j-1]) : (t_hist[dft_start+1] - t_hist[dft_start]);
                double c = cos(omega * t_hist[j]) * dtj;
                double s = sin(omega * t_hist[j]) * dtj;
                double sum_phi = 0;
                for (int a = 0; a < 3; a++) {
                    re[a] += phi_hist[a][j] * c;
                    im[a] += phi_hist[a][j] * s;
                    sum_phi += phi_hist[a][j];
                }
                re_s += sum_phi * c;
                im_s += sum_phi * s;
            }
            double pw_s = (re_s * re_s + im_s * im_s) / (T * T);
            fprintf(fdft, "%.6f", omega);
            for (int a = 0; a < 3; a++) {
                double pw = (re[a]*re[a] + im[a]*im[a]) / (T*T);
                fprintf(fdft, "\t%.6e", pw);
                if (pw > peak_pow[a]) { peak_pow[a] = pw; peak_om[a] = omega; }
            }
            fprintf(fdft, "\t%.6e\n", pw_s);
            if (pw_s > peak_sum_pow) { peak_sum_pow = pw_s; peak_sum_om = omega; }
        }
        fclose(fdft);

        printf("\nSpectrum:\n");
        for (int a = 0; a < 3; a++)
            printf("  phi%d: peak omega = %.4f (mass gap = %.4f)  %s\n",
                   a+1, peak_om[a], mass,
                   (peak_om[a] > 0.01 && peak_om[a] < mass) ? "OSCILLON" : "no");
        printf("  sum:  peak omega = %.4f  power = %.2e\n", peak_sum_om, peak_sum_pow);
        printf("  sum suppressed vs phi1? %s (ratio=%.2e)\n",
               (peak_sum_pow < 0.1 * peak_pow[0]) ? "YES (120-deg signature)" : "NO",
               (peak_pow[0] > 1e-30) ? peak_sum_pow / peak_pow[0] : 0.0);
    }

    /* Phase analysis: cross-correlate pairs at center over second half */
    if (n_dft - dft_start > 100) {
        printf("\nPhase analysis (second half):\n");
        for (int a = 0; a < 3; a++) {
            int b = (a + 1) % 3;
            double corr = 0, norm_a = 0, norm_b = 0;
            for (int j = dft_start; j < n_dft; j++) {
                corr   += phi_hist[a][j] * phi_hist[b][j];
                norm_a += phi_hist[a][j] * phi_hist[a][j];
                norm_b += phi_hist[b][j] * phi_hist[b][j];
            }
            double cos_phase = (norm_a > 1e-30 && norm_b > 1e-30) ?
                corr / sqrt(norm_a * norm_b) : 0.0;
            /* 120° → cos = -0.5; 0° → cos = +1.0 */
            printf("  cos(phi%d,phi%d) = %+.4f  [120deg: -0.50, 0deg: +1.00]\n",
                   a+1, b+1, cos_phase);
        }
    }

    printf("\nOutput: %s\n", tspath);

    for (int a = 0; a < 3; a++) {
        free(phi[a]); free(vel[a]); free(acc[a]);
        free(phi_hist[a]);
    }
    free(damp); free(t_hist);
    return 0;
}
