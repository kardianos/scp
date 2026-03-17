/*
 * pw120.c — Pairwise coupling 120-degree phase binding
 *
 * Three massive scalars with triple-product AND pairwise linear coupling:
 *   L = sum_a [ (1/2)(dt phi_a)^2 - (1/2)(dx phi_a)^2 - (m^2/2)phi_a^2 ]
 *     - (mu/2) P^2 / (1 + kappa P^2)
 *     - (lambda/2)(phi_1 phi_2 + phi_2 phi_3 + phi_3 phi_1)
 *
 * EOM: d^2 phi_a/dt^2 = d^2 phi_a/dx^2 - m^2 phi_a - lambda(phi_b+phi_c) + F_triple
 *
 * Mass spectrum: antisymmetric mode m^2_anti = m^2 - lambda (lighter)
 *                symmetric mode      m^2_sym  = m^2 + 2*lambda (heavier)
 *
 * Modes:
 *   -mode 120: Initialize with 120-degree phases
 *   -mode 0:   Initialize all in-phase (0-degree)
 *
 * Compile: gcc -O3 -Wall -o pw120 src/pw120.c -lm
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static double mu     = -20.0;
static double kappa  = 20.0;
static double mass   = 1.0;
static double lambda = 0.0;
static double A_init = 0.8;
static double sigma_w= 3.0;
static int    Nx     = 4000;
static double xmax   = 100.0;
static double tfinal = 10000.0;
static int    mode   = 120;   /* 120 or 0 */
static char   outdir[512] = "data";

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
        else if (!strcmp(argv[i], "-mode"))   mode    = atoi(argv[i+1]);
        else if (!strcmp(argv[i], "-o"))      strncpy(outdir, argv[i+1], sizeof(outdir)-1);
        else { fprintf(stderr, "Unknown flag: %s\n", argv[i]); exit(1); }
    }
}

/* -dV_triple/dphi_a where V = (mu/2)P^2/(1+kappa P^2), P = phi1 phi2 phi3 */
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

/* Measure relative phase between fields at center using Hilbert-like approach:
 * Track sign changes of phi_a(0,t) to measure period, use zero-crossing times
 * to extract phase. Simpler: measure instantaneous phase from (phi, vel) */
static double measure_phase_diff(double p1, double v1, double p2, double v2, double omega)
{
    /* theta_a = atan2(-v_a/omega, p_a) */
    if (omega < 1e-8) return 0.0;
    double th1 = atan2(-v1/omega, p1);
    double th2 = atan2(-v2/omega, p2);
    double dp = th2 - th1;
    /* wrap to [-pi, pi] */
    while (dp > M_PI)  dp -= 2*M_PI;
    while (dp < -M_PI) dp += 2*M_PI;
    return dp;
}

int main(int argc, char **argv)
{
    parse_args(argc, argv);

    double dx = 2.0 * xmax / (Nx - 1);
    double dx2 = dx * dx;
    double m2 = mass * mass;

    /* Effective mass for antisymmetric mode */
    double m2_anti = m2 - lambda;
    double m2_sym  = m2 + 2.0*lambda;
    if (m2_anti <= 0) {
        fprintf(stderr, "ERROR: lambda=%.4f >= m^2=%.4f, antisymmetric mode tachyonic!\n",
                lambda, m2);
        return 1;
    }

    /* Base frequency for initialization */
    double omega_anti = sqrt(m2_anti) * 0.9; /* slightly below mass gap */
    double omega_sym  = sqrt(m2_sym)  * 0.9;

    /* CFL: fastest mode is symmetric with m^2_sym */
    double kmax = M_PI / dx;
    double dt = 0.8 * 2.0 / sqrt(kmax*kmax + m2_sym);
    int Nt = (int)(tfinal / dt) + 1;

    printf("pw120: mode=%d lambda=%.4f\n", mode, lambda);
    printf("  mu=%.1f kappa=%.1f mass=%.4f A=%.3f sigma=%.3f\n",
           mu, kappa, mass, A_init, sigma_w);
    printf("  m2_anti=%.4f m2_sym=%.4f omega_anti=%.4f omega_sym=%.4f\n",
           m2_anti, m2_sym, omega_anti, omega_sym);
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
    double phases[3], omega_init;
    if (mode == 120) {
        phases[0] = 0.0;
        phases[1] = 2.0*M_PI/3.0;
        phases[2] = 4.0*M_PI/3.0;
        omega_init = omega_anti;
        printf("  Init: 120-degree phases, omega=%.4f\n", omega_init);
    } else {
        phases[0] = phases[1] = phases[2] = 0.0;
        omega_init = omega_sym;
        printf("  Init: 0-degree phases (symmetric), omega=%.4f\n", omega_init);
    }

    for (int a = 0; a < 3; a++)
        for (int i = 0; i < Nx; i++) {
            double x = -xmax + i * dx;
            double g = exp(-x * x / (2.0 * sigma_w * sigma_w));
            phi[a][i] =  A_init * g * cos(phases[a]);
            vel[a][i] = -omega_init * A_init * g * sin(phases[a]);
        }

    /* Compute acceleration */
    #define COMPUTE_ACC() do { \
        for (int a = 0; a < 3; a++) { \
            acc[a][0] = acc[a][1] = acc[a][Nx-2] = acc[a][Nx-1] = 0; \
            int b = (a+1)%3, c = (a+2)%3; \
            for (int i = 1; i < Nx - 1; i++) { \
                double lapl = (phi[a][i+1] - 2.0*phi[a][i] + phi[a][i-1]) / dx2; \
                double fp = force_triple(phi[0][i], phi[1][i], phi[2][i], a); \
                acc[a][i] = lapl - m2*phi[a][i] - lambda*(phi[b][i]+phi[c][i]) + fp; \
            } \
        } \
    } while(0)

    COMPUTE_ACC();

    /* Output file */
    char tspath[600];
    snprintf(tspath, sizeof(tspath), "%s/pw120_mode%d_lam%.3f_ts.tsv", outdir, mode, lambda);
    FILE *fts = fopen(tspath, "w");
    if (!fts) { fprintf(stderr, "Cannot open %s\n", tspath); return 1; }
    fprintf(fts, "time\tphi1_0\tphi2_0\tphi3_0\tpeak1\tpeak2\tpeak3\t"
                 "E_kin\tE_grad\tE_mass\tE_pw\tE_pot\tE_total\tf_core\t"
                 "phase12\tphase13\n");

    /* DFT storage — all three fields */
    int max_dft = 50000;
    double *phi_hist[3], *t_hist;
    for (int a = 0; a < 3; a++)
        phi_hist[a] = malloc(max_dft * sizeof(double));
    t_hist = malloc(max_dft * sizeof(double));
    int n_dft = 0;

    int rec_every = Nt / 20000;
    if (rec_every < 1) rec_every = 1;
    int print_every = Nt / 40;
    if (print_every < 1) print_every = 1;
    int dft_every = Nt / max_dft;
    if (dft_every < 1) dft_every = 1;

    double core_r = 3.0 * sigma_w;
    int ic = Nx / 2;

    double E_init = 0.0;

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

                /* Pairwise coupling energy */
                Epw += lambda * (phi[0][i]*phi[1][i] + phi[1][i]*phi[2][i]
                                + phi[2][i]*phi[0][i]) * dx;

                /* Triple product potential */
                double P = phi[0][i] * phi[1][i] * phi[2][i];
                double P2 = P * P;
                double V = 0.5 * mu * P2 / (1.0 + kappa * P2);
                Ep += V * dx;

                /* total density for core fraction */
                double e = V + lambda*(phi[0][i]*phi[1][i]+phi[1][i]*phi[2][i]
                                       +phi[2][i]*phi[0][i]);
                for (int a = 0; a < 3; a++) {
                    e += 0.5*vel[a][i]*vel[a][i];
                    double dp2 = (phi[a][i+1] - phi[a][i-1]) / (2.0*dx);
                    e += 0.5*dp2*dp2 + 0.5*m2*phi[a][i]*phi[a][i];
                }
                Eall += e * dx;
                if (fabs(x) < core_r) Ecore += e * dx;
            }

            double Et = Ek + Eg + Em + Epw + Ep;
            double fc = (Eall > 1e-20) ? Ecore / Eall : 0.0;
            if (n == 0) E_init = Et;

            /* Phase measurement using omega_anti */
            double omega_meas = (mode == 120) ? sqrt(fabs(m2_anti)) : sqrt(m2_sym);
            double ph12 = measure_phase_diff(phi[0][ic], vel[0][ic],
                                             phi[1][ic], vel[1][ic], omega_meas);
            double ph13 = measure_phase_diff(phi[0][ic], vel[0][ic],
                                             phi[2][ic], vel[2][ic], omega_meas);

            if (do_rec)
                fprintf(fts, "%.6f\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t"
                             "%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6f\t"
                             "%.4f\t%.4f\n",
                        t, phi[0][ic], phi[1][ic], phi[2][ic],
                        peak[0], peak[1], peak[2],
                        Ek, Eg, Em, Epw, Ep, Et, fc, ph12, ph13);
            if (do_print)
                printf("  t=%7.1f  p0=(%+.3f,%+.3f,%+.3f)  pk=(%.3f,%.3f,%.3f)  "
                       "E=%+.4f  dE/E0=%.2e  fc=%.3f  ph=%.2f,%.2f\n",
                       t, phi[0][ic], phi[1][ic], phi[2][ic],
                       peak[0], peak[1], peak[2], Et,
                       (E_init != 0) ? (Et-E_init)/fabs(E_init) : 0.0,
                       fc, ph12, ph13);
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

    /* DFT of all three fields — second half only */
    int dft_start = n_dft / 2;
    double peak_om[3] = {0}, peak_pow[3] = {0};
    if (n_dft - dft_start > 100) {
        char dftpath[600];
        snprintf(dftpath, sizeof(dftpath), "%s/pw120_mode%d_lam%.3f_spectrum.tsv",
                 outdir, mode, lambda);
        FILE *fdft = fopen(dftpath, "w");
        fprintf(fdft, "omega\tpower1\tpower2\tpower3\n");
        double T = t_hist[n_dft-1] - t_hist[dft_start];
        int nf = 500;
        for (int k = 0; k < nf; k++) {
            double omega = 3.0 * mass * k / nf;
            double pw[3] = {0};
            for (int a = 0; a < 3; a++) {
                double re = 0, im = 0;
                for (int j = dft_start; j < n_dft; j++) {
                    double dtj = (j > dft_start) ?
                        (t_hist[j]-t_hist[j-1]) : (t_hist[dft_start+1]-t_hist[dft_start]);
                    re += phi_hist[a][j] * cos(omega * t_hist[j]) * dtj;
                    im += phi_hist[a][j] * sin(omega * t_hist[j]) * dtj;
                }
                pw[a] = (re*re + im*im) / (T*T);
                if (pw[a] > peak_pow[a]) { peak_pow[a] = pw[a]; peak_om[a] = omega; }
            }
            fprintf(fdft, "%.6f\t%.6e\t%.6e\t%.6e\n", omega, pw[0], pw[1], pw[2]);
        }
        fclose(fdft);
    }

    /* Summary line for scanning */
    /* Read final-quarter diagnostics from time series */
    double fc_final = 0, dE_rel = 0;
    {
        /* Re-read final state from arrays */
        double Eall = 0, Ecore = 0;
        double Ek = 0, Eg = 0, Em = 0, Epw = 0, Ep = 0;
        for (int i = 1; i < Nx - 1; i++) {
            double x = -xmax + i * dx;
            for (int a = 0; a < 3; a++) {
                Ek += 0.5 * vel[a][i] * vel[a][i] * dx;
                double dp = (phi[a][i+1] - phi[a][i-1]) / (2.0*dx);
                Eg += 0.5 * dp * dp * dx;
                Em += 0.5 * m2 * phi[a][i] * phi[a][i] * dx;
            }
            Epw += lambda * (phi[0][i]*phi[1][i] + phi[1][i]*phi[2][i]
                            + phi[2][i]*phi[0][i]) * dx;
            double P = phi[0][i] * phi[1][i] * phi[2][i];
            double P2 = P * P;
            Ep += 0.5 * mu * P2 / (1.0 + kappa * P2) * dx;

            double e = 0.5*mu*P*P/(1.0+kappa*P*P)
                     + lambda*(phi[0][i]*phi[1][i]+phi[1][i]*phi[2][i]+phi[2][i]*phi[0][i]);
            for (int a = 0; a < 3; a++) {
                e += 0.5*vel[a][i]*vel[a][i];
                double dp2 = (phi[a][i+1] - phi[a][i-1]) / (2.0*dx);
                e += 0.5*dp2*dp2 + 0.5*m2*phi[a][i]*phi[a][i];
            }
            Eall += e * dx;
            if (fabs(x) < core_r) Ecore += e * dx;
        }
        double Et = Ek + Eg + Em + Epw + Ep;
        fc_final = (Eall > 1e-20) ? Ecore / Eall : 0.0;
        dE_rel = (E_init != 0) ? (Et - E_init) / fabs(E_init) : 0.0;
    }

    double omega_meas_f = (mode==120) ? sqrt(fabs(m2_anti)) : sqrt(m2_sym);
    double ph12_f = measure_phase_diff(phi[0][ic], vel[0][ic],
                                       phi[1][ic], vel[1][ic], omega_meas_f);
    double ph13_f = measure_phase_diff(phi[0][ic], vel[0][ic],
                                       phi[2][ic], vel[2][ic], omega_meas_f);

    printf("\n=== SUMMARY mode=%d lambda=%.4f ===\n", mode, lambda);
    printf("  omega_peak: %.4f %.4f %.4f\n", peak_om[0], peak_om[1], peak_om[2]);
    printf("  mass_gap_anti=%.4f mass_gap_sym=%.4f\n", sqrt(m2_anti), sqrt(m2_sym));
    printf("  fc_final=%.4f dE/E0=%.3e\n", fc_final, dE_rel);
    printf("  phase12=%.4f phase13=%.4f (target: %.4f)\n",
           ph12_f, ph13_f, (mode==120) ? 2*M_PI/3 : 0.0);
    printf("  oscillon? %s\n",
           (peak_om[0] > 0.01 && peak_om[0] < sqrt(m2_anti) && fc_final > 0.3) ? "YES" : "NO");
    printf("Output: %s\n", tspath);

    for (int a = 0; a < 3; a++) { free(phi[a]); free(vel[a]); free(acc[a]); free(phi_hist[a]); }
    free(damp); free(t_hist);
    return 0;
}
