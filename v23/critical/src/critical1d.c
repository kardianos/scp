/*
 * critical1d.c — V23-C: Correlation length near the gap edge
 *
 * Based on v21/src/triad1d.c. Three massive scalars with saturating
 * triple-product coupling:
 *
 *   L = sum_a [ (1/2)(dt phi_a)^2 - (1/2)(dx phi_a)^2 - (m^2/2)phi_a^2 ]
 *     - (mu/2) P^2 / (1 + kappa P^2)
 *   P = phi_1 phi_2 phi_3
 *
 * Phase 1: Scan mu from -20 to -5, measure oscillation frequency omega,
 *   tail decay rate kappa_tail, correlation length xi = 1/kappa_tail.
 *   Uses time-averaged RMS profile for robust tail fitting.
 *
 * Phase 2: Two-oscillon interaction at selected mu values.
 *   Initialize two oscillons at separation D, measure force F(D).
 *   Fit F(D) ~ F0*exp(-D/lambda), compare lambda with xi from Phase 1.
 *
 * Compile: gcc -O3 -Wall -o critical1d src/critical1d.c -lm
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* --- Parameters --- */
static double kappa    = 20.0;
static double mass     = 1.0;
static double A_init   = 1.0;
static double sig      = 3.0;
static int    Nx       = 4000;
static double xmax     = 100.0;
static double tfinal   = 10000.0;
static int    phase    = 0;      /* 0 = both, 1 = Phase 1 only, 2 = Phase 2 only */
static char   outdir[512] = "data";

/* Phase 2 parameters */
static double sep_list[] = {20.0, 30.0, 40.0, 50.0};
static int    n_sep = 4;

/* mu scan list */
static double mu_scan[] = {-20, -18, -16, -14, -12, -10, -9, -8, -7, -6, -5};
static int    n_mu = 11;

/* Phase 2: selected mu values (filled after Phase 1) */
static double mu_phase2[3];
static int    n_mu_p2 = 0;

static void parse_args(int argc, char **argv)
{
    for (int i = 1; i < argc - 1; i += 2) {
        if      (!strcmp(argv[i], "-kappa"))  kappa   = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-mass"))   mass    = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-A"))      A_init  = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-sigma"))  sig     = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-Nx"))     Nx      = atoi(argv[i+1]);
        else if (!strcmp(argv[i], "-xmax"))   xmax    = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-tfinal")) tfinal  = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-phase"))  phase   = atoi(argv[i+1]);
        else if (!strcmp(argv[i], "-o"))      strncpy(outdir, argv[i+1], sizeof(outdir)-1);
        else { fprintf(stderr, "Unknown flag: %s\n", argv[i]); exit(1); }
    }
}

/* -dV/dphi_a where V = (mu/2)P^2/(1+kappa P^2), P = phi1 phi2 phi3 */
static inline double force_pot(double p1, double p2, double p3, int a, double mu)
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

/* ===================================================================
 *  Single-oscillon evolution
 * =================================================================== */
typedef struct {
    double omega;
    double kappa_tail;
    double delta;
    double xi;
    double A_peak;
    double E_total;
    double fc_avg;       /* time-averaged f_core in second half */
    int    alive;
} oscillon_result_t;

static oscillon_result_t run_single(double mu, int save_files)
{
    oscillon_result_t res = {0};

    double dx  = 2.0 * xmax / (Nx - 1);
    double dx2 = dx * dx;
    double m2  = mass * mass;

    /* CFL */
    double kmax = M_PI / dx;
    double dt   = 0.8 * 2.0 / sqrt(kmax*kmax + m2);
    int    Nt   = (int)(tfinal / dt) + 1;

    /* Allocate fields */
    double *phi[3], *vel[3], *acc[3];
    for (int a = 0; a < 3; a++) {
        phi[a] = calloc(Nx, sizeof(double));
        vel[a] = calloc(Nx, sizeof(double));
        acc[a] = calloc(Nx, sizeof(double));
    }

    /* Time-averaged RMS profile: accumulate phi^2 during second half */
    double *phi_rms = calloc(Nx, sizeof(double));
    int rms_count = 0;

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
    int ic = Nx / 2;
    for (int a = 0; a < 3; a++)
        for (int i = 0; i < Nx; i++) {
            double x = -xmax + i * dx;
            phi[a][i] = A_init * exp(-x * x / (2.0 * sig * sig));
        }

    /* Compute acceleration */
    #define COMPUTE_ACC_S() do { \
        for (int a = 0; a < 3; a++) { \
            acc[a][0] = acc[a][1] = acc[a][Nx-2] = acc[a][Nx-1] = 0; \
            for (int i = 1; i < Nx - 1; i++) { \
                double lapl = (phi[a][i+1] - 2.0*phi[a][i] + phi[a][i-1]) / dx2; \
                double fp = force_pot(phi[0][i], phi[1][i], phi[2][i], a, mu); \
                acc[a][i] = lapl - m2*phi[a][i] + fp; \
            } \
        } \
    } while(0)

    COMPUTE_ACC_S();

    /* DFT storage: record phi1(x=0) */
    int max_dft = 100000;
    double *phi0_hist = malloc(max_dft * sizeof(double));
    double *t_hist    = malloc(max_dft * sizeof(double));
    int n_dft = 0;
    int dft_every = Nt / max_dft;
    if (dft_every < 1) dft_every = 1;

    /* Output file */
    FILE *fts = NULL;
    if (save_files) {
        char path[600];
        snprintf(path, sizeof(path), "%s/critical_mu%.0f_ts.tsv", outdir, mu);
        fts = fopen(path, "w");
        if (fts) fprintf(fts, "time\tphi1_0\tpeak\tE_kin\tE_grad\tE_mass\tE_pot\tE_total\tf_core\n");
    }

    int rec_every  = Nt / 10000;
    if (rec_every < 1) rec_every = 1;
    int print_every = Nt / 20;
    if (print_every < 1) print_every = 1;
    /* RMS accumulation: every ~100 time steps in second half */
    int rms_every = Nt / 50000;
    if (rms_every < 1) rms_every = 1;

    double core_r = 3.0 * sig;
    double last_E = 0;
    double fc_sum = 0;
    int fc_count = 0;
    double max_peak_second_half = 0;

    for (int n = 0; n <= Nt; n++) {
        double t = n * dt;
        int second_half = (t > tfinal * 0.5);

        /* Record for DFT */
        if (n % dft_every == 0 && n_dft < max_dft) {
            phi0_hist[n_dft] = phi[0][ic];
            t_hist[n_dft] = t;
            n_dft++;
        }

        /* Accumulate RMS profile in second half */
        if (second_half && (n % rms_every == 0)) {
            for (int i = 0; i < Nx; i++)
                phi_rms[i] += phi[0][i] * phi[0][i];
            rms_count++;
        }

        int do_rec   = (n % rec_every == 0);
        int do_print = (n % print_every == 0);

        if (do_rec || do_print) {
            double Ek = 0, Eg = 0, Em = 0, Ep = 0;
            double Ecore = 0, Eall = 0;
            double peak = 0;

            for (int i = 1; i < Nx - 1; i++) {
                double x = -xmax + i * dx;
                for (int a = 0; a < 3; a++) {
                    Ek += 0.5 * vel[a][i] * vel[a][i] * dx;
                    double dp = (phi[a][i+1] - phi[a][i-1]) / (2.0*dx);
                    Eg += 0.5 * dp * dp * dx;
                    Em += 0.5 * m2 * phi[a][i] * phi[a][i] * dx;
                    if (fabs(phi[a][i]) > peak) peak = fabs(phi[a][i]);
                }

                double P  = phi[0][i] * phi[1][i] * phi[2][i];
                double P2 = P * P;
                double V  = 0.5 * mu * P2 / (1.0 + kappa * P2);
                Ep += V * dx;

                double e = V;
                for (int a = 0; a < 3; a++) {
                    e += 0.5*vel[a][i]*vel[a][i];
                    double dp = (phi[a][i+1] - phi[a][i-1]) / (2.0*dx);
                    e += 0.5*dp*dp + 0.5*m2*phi[a][i]*phi[a][i];
                }
                Eall += e * dx;
                if (fabs(x) < core_r) Ecore += e * dx;
            }

            double Et = Ek + Eg + Em + Ep;
            double fc = (Eall > 1e-20) ? Ecore / Eall : 0.0;
            last_E = Et;

            if (second_half) {
                fc_sum += fc;
                fc_count++;
                if (peak > max_peak_second_half) max_peak_second_half = peak;
            }

            if (do_rec && fts)
                fprintf(fts, "%.6f\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6f\n",
                        t, phi[0][ic], peak, Ek, Eg, Em, Ep, Et, fc);

            if (do_print)
                printf("  mu=%.0f t=%7.1f  phi0=%+.4f  pk=%.4f  E=%+.4f  fc=%.3f\n",
                       mu, t, phi[0][ic], peak, Et, fc);
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

        /* Absorbing boundary */
        for (int a = 0; a < 3; a++)
            for (int i = 0; i < Nx; i++) {
                vel[a][i] *= damp[i];
                phi[a][i] *= damp[i];
            }
    }

    if (fts) fclose(fts);

    /* --- Compute time-averaged RMS envelope --- */
    if (rms_count > 0) {
        for (int i = 0; i < Nx; i++)
            phi_rms[i] = sqrt(phi_rms[i] / rms_count);
    }

    /* --- Save RMS profile --- */
    if (save_files) {
        char path[600];
        snprintf(path, sizeof(path), "%s/critical_mu%.0f_profile.tsv", outdir, mu);
        FILE *fp = fopen(path, "w");
        if (fp) {
            fprintf(fp, "x\tphi1_rms\tphi1_final\n");
            for (int i = 0; i < Nx; i++) {
                double x = -xmax + i * dx;
                fprintf(fp, "%.6f\t%.6e\t%.6e\n", x, phi_rms[i], phi[0][i]);
            }
            fclose(fp);
        }
    }

    /* --- Measure oscillation frequency from DFT (second half) --- */
    int dft_start = n_dft / 2;
    double peak_omega = 0;
    double peak_pow_val = 0;
    if (n_dft - dft_start > 100) {
        double T = t_hist[n_dft-1] - t_hist[dft_start];
        int nf = 2000;
        for (int k = 1; k < nf; k++) {
            double omega = 2.0 * mass * k / nf;
            double re = 0, im = 0;
            for (int j = dft_start; j < n_dft; j++) {
                double dtj = (j > dft_start) ?
                    (t_hist[j]-t_hist[j-1]) : (t_hist[dft_start+1]-t_hist[dft_start]);
                re += phi0_hist[j] * cos(omega * t_hist[j]) * dtj;
                im += phi0_hist[j] * sin(omega * t_hist[j]) * dtj;
            }
            double pw = (re*re + im*im) / (T*T);
            if (pw > peak_pow_val) { peak_pow_val = pw; peak_omega = omega; }
        }

        /* Save spectrum */
        if (save_files) {
            char path[600];
            snprintf(path, sizeof(path), "%s/critical_mu%.0f_spectrum.tsv", outdir, mu);
            FILE *fs = fopen(path, "w");
            if (fs) {
                fprintf(fs, "omega\tpower\n");
                for (int k = 1; k < nf; k++) {
                    double omega = 2.0 * mass * k / nf;
                    double re = 0, im = 0;
                    for (int j = dft_start; j < n_dft; j++) {
                        double dtj = (j > dft_start) ?
                            (t_hist[j]-t_hist[j-1]) : (t_hist[dft_start+1]-t_hist[dft_start]);
                        re += phi0_hist[j] * cos(omega * t_hist[j]) * dtj;
                        im += phi0_hist[j] * sin(omega * t_hist[j]) * dtj;
                    }
                    double pw = (re*re + im*im) / (T*T);
                    fprintf(fs, "%.6f\t%.6e\n", omega, pw);
                }
                fclose(fs);
            }
        }
    }

    /* --- Measure tail decay rate kappa_tail from RMS envelope --- */
    double kappa_tail = 0;
    {
        /* Fit log(phi_rms) vs x for x in [3, 20] on right side.
         * The predicted decay length is xi = 1/sqrt(m^2-w^2) ~ 2-3,
         * so the exponential tail dominates in the range [3, 20].
         * Beyond x~20 the signal is dominated by radiation/noise. */
        double sum_x = 0, sum_y = 0, sum_xx = 0, sum_xy = 0;
        int nfit = 0;
        /* First find the peak of the RMS envelope to start fitting after it */
        double rms_peak = 0;
        for (int i = ic; i < Nx; i++) {
            if (phi_rms[i] > rms_peak) rms_peak = phi_rms[i];
        }
        double fit_floor = rms_peak * 1e-6;
        for (int i = ic; i < Nx; i++) {
            double x = (i - ic) * dx;
            double amp = phi_rms[i];
            if (x > 3.0 && x < 25.0 && amp > fit_floor) {
                double lna = log(amp);
                sum_x  += x;
                sum_y  += lna;
                sum_xx += x * x;
                sum_xy += x * lna;
                nfit++;
            }
        }
        if (nfit > 10) {
            double denom = nfit * sum_xx - sum_x * sum_x;
            if (fabs(denom) > 1e-30) {
                double slope = (nfit * sum_xy - sum_x * sum_y) / denom;
                kappa_tail = -slope;
            }
        }
    }

    /* --- Fill result --- */
    res.omega      = peak_omega;
    res.kappa_tail = kappa_tail;
    res.delta      = (mass - peak_omega) / mass;
    res.xi         = (kappa_tail > 0.001) ? 1.0 / kappa_tail : 9999.0;
    res.A_peak     = max_peak_second_half;
    res.E_total    = last_E;
    res.fc_avg     = (fc_count > 0) ? fc_sum / fc_count : 0.0;
    /* Alive: need significant peak amplitude AND f_core > 0.8 AND omega below gap */
    res.alive      = (max_peak_second_half > 0.1 && res.fc_avg > 0.8 &&
                      peak_omega > 0.01 && peak_omega < mass);

    /* Cleanup */
    for (int a = 0; a < 3; a++) { free(phi[a]); free(vel[a]); free(acc[a]); }
    free(damp); free(phi0_hist); free(t_hist); free(phi_rms);

    return res;
}

/* ===================================================================
 *  Phase 2: Two-oscillon interaction
 * =================================================================== */
typedef struct {
    double mu;
    double D_init;
    double D_final;
    double F_avg;
    double lambda;
} interaction_result_t;

static double find_centroid_right(double *phi_arr, int Nx_l, double dx_l, double xmax_l)
{
    double num = 0, den = 0;
    int ic_l = Nx_l / 2;
    for (int i = ic_l; i < Nx_l; i++) {
        double x = -xmax_l + i * dx_l;
        double w = phi_arr[i] * phi_arr[i];
        num += x * w;
        den += w;
    }
    return (den > 1e-30) ? num / den : 0;
}

static double find_centroid_left(double *phi_arr, int Nx_l, double dx_l, double xmax_l)
{
    double num = 0, den = 0;
    int ic_l = Nx_l / 2;
    for (int i = 0; i < ic_l; i++) {
        double x = -xmax_l + i * dx_l;
        double w = phi_arr[i] * phi_arr[i];
        num += x * w;
        den += w;
    }
    return (den > 1e-30) ? num / den : 0;
}

static interaction_result_t run_two_oscillon(double mu, double D,
                                              double tfinal_2, int save)
{
    interaction_result_t res = {0};
    res.mu = mu;
    res.D_init = D;

    int Nx2 = 8000;
    double xmax2 = 200.0;
    double dx  = 2.0 * xmax2 / (Nx2 - 1);
    double dx2 = dx * dx;
    double m2  = mass * mass;

    double kmax_l = M_PI / dx;
    double dt   = 0.8 * 2.0 / sqrt(kmax_l*kmax_l + m2);
    int    Nt   = (int)(tfinal_2 / dt) + 1;

    double *phi[3], *vel[3], *acc[3];
    for (int a = 0; a < 3; a++) {
        phi[a] = calloc(Nx2, sizeof(double));
        vel[a] = calloc(Nx2, sizeof(double));
        acc[a] = calloc(Nx2, sizeof(double));
    }

    double *damp = malloc(Nx2 * sizeof(double));
    double x_abs = xmax2 * 0.75;
    for (int i = 0; i < Nx2; i++) {
        double x = -xmax2 + i * dx;
        double ax = fabs(x);
        if (ax > x_abs) {
            double f = (ax - x_abs) / (xmax2 - x_abs);
            damp[i] = 1.0 - 0.98 * f * f;
        } else {
            damp[i] = 1.0;
        }
    }

    /* Initialize: two Gaussians at x = +/- D/2 */
    for (int a = 0; a < 3; a++)
        for (int i = 0; i < Nx2; i++) {
            double x = -xmax2 + i * dx;
            double g1 = A_init * exp(-(x - D/2)*(x - D/2) / (2.0 * sig * sig));
            double g2 = A_init * exp(-(x + D/2)*(x + D/2) / (2.0 * sig * sig));
            phi[a][i] = g1 + g2;
        }

    /* Compute acceleration - use local macro to avoid Nx conflict */
    int Nx_2 = Nx2;  /* local copy for macro */

    /* Track separation vs time */
    int max_track = 10000;
    double *sep_hist = malloc(max_track * sizeof(double));
    double *t_track  = malloc(max_track * sizeof(double));
    int n_track = 0;
    int track_every = Nt / max_track;
    if (track_every < 1) track_every = 1;

    FILE *fts = NULL;
    if (save) {
        char path[600];
        snprintf(path, sizeof(path), "%s/interact_mu%.0f_D%.0f_ts.tsv", outdir, mu, D);
        fts = fopen(path, "w");
        if (fts) fprintf(fts, "time\tseparation\tx_right\tx_left\n");
    }

    int print_every = Nt / 20;
    if (print_every < 1) print_every = 1;

    /* Initial acceleration */
    for (int a = 0; a < 3; a++) {
        acc[a][0] = acc[a][1] = acc[a][Nx_2-2] = acc[a][Nx_2-1] = 0;
        for (int i = 1; i < Nx_2 - 1; i++) {
            double lapl = (phi[a][i+1] - 2.0*phi[a][i] + phi[a][i-1]) / dx2;
            double fp = force_pot(phi[0][i], phi[1][i], phi[2][i], a, mu);
            acc[a][i] = lapl - m2*phi[a][i] + fp;
        }
    }

    for (int n = 0; n <= Nt; n++) {
        double t = n * dt;

        if (n % track_every == 0) {
            double xr = find_centroid_right(phi[0], Nx_2, dx, xmax2);
            double xl = find_centroid_left(phi[0], Nx_2, dx, xmax2);
            double sep = xr - xl;

            if (n_track < max_track) {
                sep_hist[n_track] = sep;
                t_track[n_track]  = t;
                n_track++;
            }

            if (fts)
                fprintf(fts, "%.6f\t%.6f\t%.6f\t%.6f\n", t, sep, xr, xl);

            if (n % print_every == 0)
                printf("  mu=%.0f D=%.0f t=%7.1f  sep=%.3f  xR=%.3f  xL=%.3f\n",
                       mu, D, t, sep, xr, xl);
        }

        if (n == Nt) break;

        /* Velocity Verlet */
        for (int a = 0; a < 3; a++)
            for (int i = 1; i < Nx_2 - 1; i++)
                vel[a][i] += 0.5 * dt * acc[a][i];
        for (int a = 0; a < 3; a++)
            for (int i = 1; i < Nx_2 - 1; i++)
                phi[a][i] += dt * vel[a][i];

        /* Recompute acc */
        for (int a = 0; a < 3; a++) {
            acc[a][0] = acc[a][1] = acc[a][Nx_2-2] = acc[a][Nx_2-1] = 0;
            for (int i = 1; i < Nx_2 - 1; i++) {
                double lapl = (phi[a][i+1] - 2.0*phi[a][i] + phi[a][i-1]) / dx2;
                double fp = force_pot(phi[0][i], phi[1][i], phi[2][i], a, mu);
                acc[a][i] = lapl - m2*phi[a][i] + fp;
            }
        }

        for (int a = 0; a < 3; a++)
            for (int i = 1; i < Nx_2 - 1; i++)
                vel[a][i] += 0.5 * dt * acc[a][i];

        for (int a = 0; a < 3; a++)
            for (int i = 0; i < Nx_2; i++) {
                vel[a][i] *= damp[i];
                phi[a][i] *= damp[i];
            }
    }

    if (fts) fclose(fts);

    /* Measure average force from second half */
    if (n_track > 20) {
        int half = n_track / 2;
        double t0 = t_track[half];

        /* Accumulate sums for quadratic fit s(t) = a0 + a1*t + a2*t^2 */
        double S[5] = {0};  /* S[k] = sum(tj^k) */
        double Sy[3] = {0}; /* Sy[k] = sum(sj * tj^k) */
        int nn = n_track - half;

        for (int j = half; j < n_track; j++) {
            double tj = t_track[j] - t0;
            double sj = sep_hist[j];
            double tk = 1.0;
            for (int k = 0; k < 5; k++) { S[k] += tk; tk *= tj; }
            Sy[0] += sj;
            Sy[1] += sj * tj;
            Sy[2] += sj * tj * tj;
        }

        /* 3x3 normal equations: [S0 S1 S2; S1 S2 S3; S2 S3 S4] [a0;a1;a2] = [Sy0;Sy1;Sy2] */
        double M[3][4] = {
            {S[0], S[1], S[2], Sy[0]},
            {S[1], S[2], S[3], Sy[1]},
            {S[2], S[3], S[4], Sy[2]}
        };

        /* Gaussian elimination with partial pivoting */
        for (int col = 0; col < 3; col++) {
            int piv = col;
            for (int r = col+1; r < 3; r++)
                if (fabs(M[r][col]) > fabs(M[piv][col])) piv = r;
            if (piv != col)
                for (int c = 0; c < 4; c++) {
                    double tmp = M[col][c]; M[col][c] = M[piv][c]; M[piv][c] = tmp;
                }
            if (fabs(M[col][col]) < 1e-30) continue;
            for (int r = col+1; r < 3; r++) {
                double fac = M[r][col] / M[col][col];
                for (int c = col; c < 4; c++)
                    M[r][c] -= fac * M[col][c];
            }
        }
        double coeff[3] = {0};
        for (int r = 2; r >= 0; r--) {
            if (fabs(M[r][r]) < 1e-30) continue;
            coeff[r] = M[r][3];
            for (int c = r+1; c < 3; c++)
                coeff[r] -= M[r][c] * coeff[c];
            coeff[r] /= M[r][r];
        }

        res.F_avg = 2.0 * coeff[2];
        res.D_final = sep_hist[n_track - 1];
        (void)nn;
    }

    for (int a = 0; a < 3; a++) { free(phi[a]); free(vel[a]); free(acc[a]); }
    free(damp); free(sep_hist); free(t_track);

    return res;
}

/* =================================================================== */
int main(int argc, char **argv)
{
    parse_args(argc, argv);

    printf("=== V23-C: Critical Gravity -- Correlation Length Near Gap Edge ===\n\n");
    printf("Parameters: kappa=%.1f  mass=%.3f  A=%.2f  sigma=%.2f\n", kappa, mass, A_init, sig);
    printf("Grid: Nx=%d  xmax=%.1f  tfinal=%.0f\n\n", Nx, xmax, tfinal);

    /* ===== PHASE 1: Correlation length scan ===== */
    if (phase == 0 || phase == 1) {
        printf("===== PHASE 1: Correlation Length vs Gap Margin =====\n\n");
        printf("%-6s  %-8s  %-10s  %-10s  %-10s  %-8s  %-10s  %-6s  %-5s  %-10s  %-10s\n",
               "mu", "omega", "kappa_tail", "delta", "xi", "A_peak", "E_total",
               "fc_avg", "alive", "xi_pred", "xi/xi_p");
        printf("------  --------  ----------  ----------  ----------  --------  "
               "----------  ------  -----  ----------  ----------\n");

        oscillon_result_t results[20];
        int n_alive = 0;

        for (int mi = 0; mi < n_mu; mi++) {
            double mu_i = mu_scan[mi];
            printf("\n--- Running mu = %.0f ---\n", mu_i);

            results[mi] = run_single(mu_i, 1);
            oscillon_result_t *r = &results[mi];

            /* Predicted xi = 1/sqrt(m^2 - omega^2) */
            double xi_pred = 0;
            if (r->delta > 0 && r->delta < 1.0) {
                double k2 = mass*mass - r->omega*r->omega;
                if (k2 > 0) xi_pred = 1.0 / sqrt(k2);
            }

            printf("\n%-6.0f  %-8.4f  %-10.4f  %-10.6f  %-10.4f  %-8.4f  %-10.4f  %-6.3f  %-5s  %-10.4f  %-10.4f\n",
                   mu_i, r->omega, r->kappa_tail, r->delta, r->xi,
                   r->A_peak, r->E_total, r->fc_avg, r->alive ? "YES" : "NO",
                   xi_pred, (xi_pred > 0.01) ? r->xi / xi_pred : 0.0);

            if (r->alive) n_alive++;
        }

        /* Print summary table */
        printf("\n\n===== PHASE 1 SUMMARY =====\n\n");
        printf("%-6s  %-8s  %-10s  %-10s  %-10s  %-8s  %-10s  %-6s  %-5s  %-10s  %-10s\n",
               "mu", "omega", "kappa_tail", "delta", "xi", "A_peak", "E_total",
               "fc_avg", "alive", "xi_pred", "xi/xi_p");
        printf("------  --------  ----------  ----------  ----------  --------  "
               "----------  ------  -----  ----------  ----------\n");

        int idx_standard = -1, idx_interm = -1, idx_critical = -1;

        for (int mi = 0; mi < n_mu; mi++) {
            oscillon_result_t *r = &results[mi];
            double xi_pred = 0;
            if (r->delta > 0 && r->delta < 1.0) {
                double k2 = mass*mass - r->omega*r->omega;
                if (k2 > 0) xi_pred = 1.0 / sqrt(k2);
            }
            printf("%-6.0f  %-8.4f  %-10.4f  %-10.6f  %-10.4f  %-8.4f  %-10.4f  %-6.3f  %-5s  %-10.4f  %-10.4f\n",
                   mu_scan[mi], r->omega, r->kappa_tail, r->delta, r->xi,
                   r->A_peak, r->E_total, r->fc_avg, r->alive ? "YES" : "NO",
                   xi_pred, (xi_pred > 0.01) ? r->xi / xi_pred : 0.0);

            if (r->alive) {
                if (idx_standard < 0) idx_standard = mi;
                idx_critical = mi;
            }
        }

        /* Intermediate */
        if (idx_standard >= 0 && idx_critical >= 0) {
            idx_interm = (idx_standard + idx_critical) / 2;
            if (idx_interm == idx_standard && idx_critical > idx_standard)
                idx_interm = idx_standard + 1;
        }

        /* Find critical mu_c */
        printf("\n");
        int found_critical = 0;
        for (int mi = 0; mi < n_mu - 1; mi++) {
            if (results[mi].alive && !results[mi+1].alive) {
                printf("CRITICAL mu_c: between %.0f and %.0f\n",
                       mu_scan[mi], mu_scan[mi+1]);
                found_critical = 1;
                break;
            }
        }
        if (!found_critical) {
            /* Check if last alive borders end of scan */
            if (idx_critical == n_mu - 1)
                printf("All oscillons survived -- mu_c is below %.0f\n", mu_scan[n_mu-1]);
            else if (idx_standard < 0)
                printf("No oscillons survived at any mu!\n");
        }

        printf("Alive oscillons: %d / %d\n\n", n_alive, n_mu);

        /* Scaling check: does xi ~ 1/sqrt(2*m*delta)? */
        printf("Scaling check: xi vs 1/sqrt(m^2 - omega^2)\n");
        printf("%-6s  %-10s  %-10s  %-10s\n", "mu", "xi_meas", "xi_pred", "ratio");
        printf("------  ----------  ----------  ----------\n");
        for (int mi = 0; mi < n_mu; mi++) {
            oscillon_result_t *r = &results[mi];
            if (!r->alive) continue;
            double k2 = mass*mass - r->omega*r->omega;
            double xi_pred = (k2 > 0) ? 1.0/sqrt(k2) : 0;
            printf("%-6.0f  %-10.4f  %-10.4f  %-10.4f\n",
                   mu_scan[mi], r->xi, xi_pred,
                   (xi_pred > 0.01) ? r->xi / xi_pred : 0.0);
        }
        printf("\n");

        /* Set up Phase 2 candidates */
        if (idx_standard >= 0) {
            mu_phase2[n_mu_p2++] = mu_scan[idx_standard];
            if (idx_interm >= 0 && idx_interm != idx_standard && results[idx_interm].alive)
                mu_phase2[n_mu_p2++] = mu_scan[idx_interm];
            if (idx_critical >= 0 && idx_critical != idx_standard &&
                idx_critical != idx_interm && results[idx_critical].alive)
                mu_phase2[n_mu_p2++] = mu_scan[idx_critical];
        }
    }

    /* ===== PHASE 2: Two-oscillon interaction ===== */
    if ((phase == 0 && n_mu_p2 > 0) || phase == 2) {
        printf("===== PHASE 2: Two-Oscillon Interaction =====\n\n");

        if (phase == 2) {
            mu_phase2[0] = -20.0;
            mu_phase2[1] = -12.0;
            mu_phase2[2] = -8.0;
            n_mu_p2 = 3;
        }

        double tfinal_2 = 5000.0;

        printf("Phase 2 mu values:");
        for (int k = 0; k < n_mu_p2; k++) printf(" %.0f", mu_phase2[k]);
        printf("\nSeparations:");
        for (int s = 0; s < n_sep; s++) printf(" %.0f", sep_list[s]);
        printf("\n\n");

        printf("%-6s  %-6s  %-10s  %-10s  %-12s\n",
               "mu", "D", "D_final", "F_avg", "dD");
        printf("------  ------  ----------  ----------  ------------\n");

        for (int k = 0; k < n_mu_p2; k++) {
            double mu_k = mu_phase2[k];
            printf("\n--- mu = %.0f ---\n", mu_k);

            double F_vals[10], D_vals[10];
            int nf = 0;

            for (int s = 0; s < n_sep; s++) {
                double D = sep_list[s];
                printf("\n  Running D=%.0f ...\n", D);
                interaction_result_t ir = run_two_oscillon(mu_k, D, tfinal_2, 1);

                printf("  %-6.0f  %-6.0f  %-10.3f  %-10.3e  %-12.3f\n",
                       mu_k, D, ir.D_final, ir.F_avg, ir.D_final - D);

                if (fabs(ir.F_avg) > 1e-15) {
                    F_vals[nf] = ir.F_avg;
                    D_vals[nf] = D;
                    nf++;
                }
            }

            /* Fit F(D) = F0 * exp(-D/lambda) */
            if (nf >= 2) {
                double sx = 0, sy = 0, sxx = 0, sxy = 0;
                int nfit = 0;
                for (int j = 0; j < nf; j++) {
                    if (fabs(F_vals[j]) > 1e-20) {
                        double lnF = log(fabs(F_vals[j]));
                        sx += D_vals[j];
                        sy += lnF;
                        sxx += D_vals[j] * D_vals[j];
                        sxy += D_vals[j] * lnF;
                        nfit++;
                    }
                }
                if (nfit >= 2) {
                    double det = nfit * sxx - sx * sx;
                    if (fabs(det) > 1e-30) {
                        double slope = (nfit * sxy - sx * sy) / det;
                        double lam = -1.0 / slope;
                        printf("\n  Fit: lambda = %.3f (interaction range)\n", lam);
                    }
                }
            }
            printf("\n");
        }
    }

    printf("\n=== Done ===\n");
    return 0;
}
