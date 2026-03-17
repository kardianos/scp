/*
 * maxwell_f.c — V24-MF: Z2 -> U(1) Promotion via Goldstone Mode
 *
 * Three massive scalar fields phi_a with triple-product coupling (v21 model)
 * plus a massless scalar theta coupled via:
 *
 *   □theta = g * d_t(P^2)   where P = phi_1 * phi_2 * phi_3
 *
 * The oscillating oscillon (P ~ cos^3(omega*t)) sources theta at frequencies
 * omega, 3*omega. Since theta is massless, it propagates at c and creates
 * a long-range (linear in 1D) profile around the oscillon.
 *
 * Modes:
 *   1: Single oscillon + theta field (Goldstone radiation)
 *   2: Two oscillons separated by D, theta mediates interaction
 *   3: Control: single oscillon, g=0 (no theta coupling)
 *
 * Compile: gcc -O3 -Wall -o maxwell_f src/maxwell_f.c -lm
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Physics parameters */
static double mu     = -20.0;
static double kappa  = 20.0;
static double mass   = 1.0;
static double A_init = 0.8;
static double sigma  = 3.0;
static double g_coup = 0.1;    /* theta coupling */

/* Grid */
static int    Nx     = 4000;
static double xmax   = 100.0;
static double tfinal = 10000.0;

/* Two-oscillon */
static double D_sep  = 40.0;   /* separation for mode 2 */
static double phase2 = 0.0;    /* phase offset of second oscillon */

static int    mode   = 1;
static char   outdir[512] = "v24/maxwell_f/data";
static char   tag[128]    = "";

static void parse_args(int argc, char **argv)
{
    for (int i = 1; i < argc - 1; i += 2) {
        if      (!strcmp(argv[i], "-mu"))     mu     = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-kappa"))  kappa  = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-mass"))   mass   = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-A"))      A_init = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-sigma"))  sigma  = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-g"))      g_coup = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-Nx"))     Nx     = atoi(argv[i+1]);
        else if (!strcmp(argv[i], "-xmax"))   xmax   = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-tfinal")) tfinal = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-D"))      D_sep  = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-phase"))  phase2 = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-mode"))   mode   = atoi(argv[i+1]);
        else if (!strcmp(argv[i], "-tag"))    strncpy(tag, argv[i+1], sizeof(tag)-1);
        else if (!strcmp(argv[i], "-o"))      strncpy(outdir, argv[i+1], sizeof(outdir)-1);
        else { fprintf(stderr, "Unknown flag: %s\n", argv[i]); exit(1); }
    }
}

/* -dV/dphi_a where V = (mu/2)P^2/(1+kappa*P^2), P = phi1*phi2*phi3 */
static inline double force_pot(double p1, double p2, double p3, int a)
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

    double dx = 2.0 * xmax / (Nx - 1);
    double dx2 = dx * dx;
    double m2 = mass * mass;

    /* CFL: massless theta has c=1, so dt < dx. Also need dt for massive fields. */
    double dt_mass = 0.8 * 2.0 / sqrt((M_PI/dx)*(M_PI/dx) + m2);
    double dt_wave = 0.8 * dx;  /* CFL for massless wave */
    double dt = (dt_mass < dt_wave) ? dt_mass : dt_wave;
    int Nt = (int)(tfinal / dt) + 1;

    const char *mode_desc;
    switch (mode) {
        case 1: mode_desc = "Single oscillon + theta"; break;
        case 2: mode_desc = "Two oscillons + theta"; break;
        case 3: mode_desc = "Control (g=0)"; g_coup = 0.0; break;
        default: fprintf(stderr, "Unknown mode %d\n", mode); return 1;
    }

    printf("maxwell_f mode %d: %s\n", mode, mode_desc);
    printf("  mu=%.3f kappa=%.4f mass=%.4f A=%.3f sigma=%.3f g=%.4f\n",
           mu, kappa, mass, A_init, sigma, g_coup);
    printf("  Nx=%d xmax=%.1f dx=%.5f dt=%.6f Nt=%d tfinal=%.0f\n",
           Nx, xmax, dx, dt, Nt, tfinal);
    if (mode == 2)
        printf("  D=%.1f phase=%.4f\n", D_sep, phase2);
    fflush(stdout);

    /* Allocate: 3 phi fields + 1 theta field */
    double *phi[3], *vel[3], *acc[3];
    double *theta, *vth, *ath;
    double *P2_cur, *P2_prev;  /* P^2 at current and previous step for d_t(P^2) */

    for (int a = 0; a < 3; a++) {
        phi[a] = calloc(Nx, sizeof(double));
        vel[a] = calloc(Nx, sizeof(double));
        acc[a] = calloc(Nx, sizeof(double));
    }
    theta  = calloc(Nx, sizeof(double));
    vth    = calloc(Nx, sizeof(double));
    ath    = calloc(Nx, sizeof(double));
    P2_cur  = calloc(Nx, sizeof(double));
    P2_prev = calloc(Nx, sizeof(double));

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

    /* Initialize phi fields */
    if (mode == 2) {
        /* Two oscillons at x = +D/2 and x = -D/2 */
        double x1 = +D_sep / 2.0;
        double x2 = -D_sep / 2.0;
        double cos_ph = cos(phase2);
        for (int a = 0; a < 3; a++)
            for (int i = 0; i < Nx; i++) {
                double x = -xmax + i * dx;
                double g1 = A_init * exp(-(x-x1)*(x-x1) / (2.0*sigma*sigma));
                double g2 = A_init * exp(-(x-x2)*(x-x2) / (2.0*sigma*sigma));
                phi[a][i] = g1 + cos_ph * g2;
            }
    } else {
        /* Single oscillon at x=0 */
        for (int a = 0; a < 3; a++)
            for (int i = 0; i < Nx; i++) {
                double x = -xmax + i * dx;
                phi[a][i] = A_init * exp(-x*x / (2.0*sigma*sigma));
            }
    }

    /* Initialize P^2 */
    for (int i = 0; i < Nx; i++) {
        double P = phi[0][i] * phi[1][i] * phi[2][i];
        P2_cur[i] = P * P;
        P2_prev[i] = P2_cur[i];  /* initially stationary */
    }

    /* Compute acceleration macro */
    /* Lagrangian coupling: L_int = g * (d_t theta) * P^2
     * phi_a EOM: acc[a] += 2*g*(d_t theta)*P*(dP/dphi_a)  [backreaction]
     * theta EOM: d_tt theta = d_xx theta + g * d_t(P^2)
     * d_t(P^2) approximated as (P2_cur - P2_prev) / dt */
    #define COMPUTE_ACC() do { \
        for (int a = 0; a < 3; a++) { \
            acc[a][0] = acc[a][1] = acc[a][Nx-2] = acc[a][Nx-1] = 0; \
            for (int i = 1; i < Nx - 1; i++) { \
                double lapl = (phi[a][i+1] - 2.0*phi[a][i] + phi[a][i-1]) / dx2; \
                double fp = force_pot(phi[0][i], phi[1][i], phi[2][i], a); \
                /* Backreaction from theta: +2*g*(d_t theta)*P*(dP/dphi_a) */ \
                double P_loc = phi[0][i] * phi[1][i] * phi[2][i]; \
                double dPda; \
                switch (a) { \
                    case 0: dPda = phi[1][i] * phi[2][i]; break; \
                    case 1: dPda = phi[0][i] * phi[2][i]; break; \
                    case 2: dPda = phi[0][i] * phi[1][i]; break; \
                    default: dPda = 0.0; \
                } \
                double f_back = 2.0 * g_coup * vth[i] * P_loc * dPda; \
                acc[a][i] = lapl - m2*phi[a][i] + fp + f_back; \
            } \
        } \
        /* theta acceleration: □theta = g * d_t(P^2) */ \
        /* => d_tt theta = d_xx theta + g * d_t(P^2) */ \
        ath[0] = ath[1] = ath[Nx-2] = ath[Nx-1] = 0; \
        for (int i = 1; i < Nx - 1; i++) { \
            double lapl_th = (theta[i+1] - 2.0*theta[i] + theta[i-1]) / dx2; \
            double dP2dt = (P2_cur[i] - P2_prev[i]) / dt; \
            ath[i] = lapl_th + g_coup * dP2dt; \
        } \
    } while(0)

    COMPUTE_ACC();

    /* Output files */
    char tspath[600], profpath[600];
    if (tag[0])
        snprintf(tspath, sizeof(tspath), "%s/mf_mode%d_%s_ts.tsv", outdir, mode, tag);
    else
        snprintf(tspath, sizeof(tspath), "%s/mf_mode%d_g%.4f_ts.tsv", outdir, mode, g_coup);
    FILE *fts = fopen(tspath, "w");
    if (!fts) { fprintf(stderr, "Cannot open %s\n", tspath); return 1; }

    if (mode == 2) {
        fprintf(fts, "time\tphi1_0\tpeak_phi\tE_phi\tE_theta\tE_total\tf_core\t"
                     "theta_0\ttheta_mid\ttheta_far\t"
                     "z_upper\tz_lower\tseparation\tpk_upper\tpk_lower\n");
    } else {
        fprintf(fts, "time\tphi1_0\tpeak_phi\tE_phi\tE_theta\tE_total\tf_core\t"
                     "theta_0\ttheta_10\ttheta_30\ttheta_50\tmax_theta\n");
    }

    /* DFT storage */
    int max_dft = 50000;
    double *phi0_hist = malloc(max_dft * sizeof(double));
    double *th0_hist  = malloc(max_dft * sizeof(double));
    double *t_hist    = malloc(max_dft * sizeof(double));
    int n_dft = 0;

    int rec_every = Nt / 20000;
    if (rec_every < 1) rec_every = 1;
    int print_every = Nt / 100;
    if (print_every < 1) print_every = 1;
    int dft_every = Nt / max_dft;
    if (dft_every < 1) dft_every = 1;

    /* Profile snapshots */
    int n_profiles = 20;
    double profile_interval = tfinal / n_profiles;
    double next_profile_t = profile_interval;

    double core_r = 3.0 * sigma;
    int ic = Nx / 2;

    /* Two-oscillon centroid indices */
    int ic_upper = -1, ic_lower = -1;
    if (mode == 2) {
        ic_upper = ic + (int)(D_sep / (2.0 * dx));
        ic_lower = ic - (int)(D_sep / (2.0 * dx));
    }

    printf("  Starting time evolution...\n");
    fflush(stdout);

    for (int n = 0; n <= Nt; n++) {
        double t = n * dt;

        /* DFT sampling */
        if (n % dft_every == 0 && n_dft < max_dft) {
            phi0_hist[n_dft] = phi[0][ic];
            th0_hist[n_dft]  = theta[ic];
            t_hist[n_dft] = t;
            n_dft++;
        }

        int do_rec   = (n % rec_every == 0);
        int do_print = (n % print_every == 0);

        /* Profile snapshots */
        if (t >= next_profile_t - dt/2 && t > 0.5) {
            if (tag[0])
                snprintf(profpath, sizeof(profpath), "%s/mf_mode%d_%s_prof_t%06d.tsv",
                         outdir, mode, tag, (int)t);
            else
                snprintf(profpath, sizeof(profpath), "%s/mf_mode%d_g%.4f_prof_t%06d.tsv",
                         outdir, mode, g_coup, (int)t);
            FILE *fp = fopen(profpath, "w");
            if (fp) {
                fprintf(fp, "x\tphi1\tphi2\tphi3\ttheta\tP2\trho_phi\n");
                for (int i = 0; i < Nx; i += 2) {  /* every other point for size */
                    double x = -xmax + i * dx;
                    double P = phi[0][i] * phi[1][i] * phi[2][i];
                    double rho = 0;
                    for (int a = 0; a < 3; a++) {
                        rho += 0.5 * vel[a][i] * vel[a][i];
                        if (i > 0 && i < Nx-1) {
                            double dp = (phi[a][i+1] - phi[a][i-1]) / (2.0*dx);
                            rho += 0.5 * dp * dp;
                        }
                        rho += 0.5 * m2 * phi[a][i] * phi[a][i];
                    }
                    rho += 0.5 * mu * P*P / (1.0 + kappa * P*P);
                    fprintf(fp, "%.4f\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\n",
                            x, phi[0][i], phi[1][i], phi[2][i], theta[i], P*P, rho);
                }
                fclose(fp);
                if (do_print)
                    printf("  [profile saved: t=%d]\n", (int)t);
            }
            next_profile_t += profile_interval;
        }

        if (do_rec || do_print) {
            double Ek = 0, Eg = 0, Em = 0, Ep = 0;
            double Eth_kin = 0, Eth_grad = 0;
            double Ecore = 0, Eall = 0;
            double peak_phi = 0, max_theta = 0;

            for (int i = 1; i < Nx - 1; i++) {
                double x = -xmax + i * dx;
                for (int a = 0; a < 3; a++) {
                    Ek += 0.5 * vel[a][i] * vel[a][i] * dx;
                    double dp = (phi[a][i+1] - phi[a][i-1]) / (2.0*dx);
                    Eg += 0.5 * dp * dp * dx;
                    Em += 0.5 * m2 * phi[a][i] * phi[a][i] * dx;
                    if (fabs(phi[a][i]) > peak_phi) peak_phi = fabs(phi[a][i]);
                }
                double P = phi[0][i] * phi[1][i] * phi[2][i];
                double P2 = P * P;
                double V = 0.5 * mu * P2 / (1.0 + kappa * P2);
                Ep += V * dx;

                /* theta energy */
                Eth_kin += 0.5 * vth[i] * vth[i] * dx;
                double dth = (theta[i+1] - theta[i-1]) / (2.0*dx);
                Eth_grad += 0.5 * dth * dth * dx;

                if (fabs(theta[i]) > max_theta) max_theta = fabs(theta[i]);

                /* core fraction */
                double e = V;
                for (int a = 0; a < 3; a++) {
                    e += 0.5*vel[a][i]*vel[a][i];
                    double dp2 = (phi[a][i+1] - phi[a][i-1]) / (2.0*dx);
                    e += 0.5*dp2*dp2 + 0.5*m2*phi[a][i]*phi[a][i];
                }
                Eall += e * dx;
                if (fabs(x) < core_r) Ecore += e * dx;
            }

            double Ephi = Ek + Eg + Em + Ep;
            double Etheta = Eth_kin + Eth_grad;
            double Etotal = Ephi + Etheta;
            double fc = (Eall > 1e-20) ? Ecore / Eall : 0.0;

            /* theta at specific distances */
            int i10 = ic + (int)(10.0/dx);
            int i30 = ic + (int)(30.0/dx);
            int i50 = ic + (int)(50.0/dx);
            double th0  = theta[ic];
            double th10 = (i10 < Nx) ? theta[i10] : 0;
            double th30 = (i30 < Nx) ? theta[i30] : 0;
            double th50 = (i50 < Nx) ? theta[i50] : 0;

            if (mode == 2) {
                /* Two-oscillon diagnostics */
                double wz_up = 0, w_up = 0, wz_lo = 0, w_lo = 0;
                double pk_up = 0, pk_lo = 0;
                for (int i = 1; i < Nx - 1; i++) {
                    double x = -xmax + i * dx;
                    double e = 0;
                    for (int a = 0; a < 3; a++) {
                        e += 0.5 * vel[a][i] * vel[a][i];
                        double dp = (phi[a][i+1] - phi[a][i-1]) / (2.0*dx);
                        e += 0.5 * dp * dp;
                        e += 0.5 * m2 * phi[a][i] * phi[a][i];
                    }
                    double P = phi[0][i] * phi[1][i] * phi[2][i];
                    e += 0.5 * mu * P*P / (1.0 + kappa * P*P);
                    double edx = e * dx;

                    if (x > 0) {
                        wz_up += x * edx;
                        w_up += edx;
                        if (fabs(phi[0][i]) > pk_up) pk_up = fabs(phi[0][i]);
                    } else {
                        wz_lo += x * edx;
                        w_lo += edx;
                        if (fabs(phi[0][i]) > pk_lo) pk_lo = fabs(phi[0][i]);
                    }
                }
                double z_up = (w_up > 1e-20) ? wz_up / w_up : D_sep/2;
                double z_lo = (w_lo > 1e-20) ? wz_lo / w_lo : -D_sep/2;
                double sep = z_up - z_lo;
                double th_mid = theta[ic];

                if (do_rec)
                    fprintf(fts, "%.6f\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6f\t"
                                 "%.6e\t%.6e\t%.6e\t"
                                 "%.6f\t%.6f\t%.6f\t%.6e\t%.6e\n",
                            t, phi[0][ic], peak_phi, Ephi, Etheta, Etotal, fc,
                            th0, th_mid, th50,
                            z_up, z_lo, sep, pk_up, pk_lo);
                if (do_print)
                    printf("  t=%7.0f  sep=%.3f  pk=(%.3f,%.3f)  E=%.3f(phi)+%.3e(th)  "
                           "fc=%.3f  th_mid=%.3e\n",
                           t, sep, pk_up, pk_lo, Ephi, Etheta, fc, th_mid);
            } else {
                if (do_rec)
                    fprintf(fts, "%.6f\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6f\t"
                                 "%.6e\t%.6e\t%.6e\t%.6e\t%.6e\n",
                            t, phi[0][ic], peak_phi, Ephi, Etheta, Etotal, fc,
                            th0, th10, th30, th50, max_theta);
                if (do_print)
                    printf("  t=%7.0f  p0=%+.4f  pk=%.3f  E=%.3f(phi)+%.3e(th)  "
                           "fc=%.3f  th(0,10,30,50)=(%.2e,%.2e,%.2e,%.2e)\n",
                           t, phi[0][ic], peak_phi, Ephi, Etheta, fc,
                           th0, th10, th30, th50);
            }
        }

        if (n == Nt) break;

        /* --- Velocity Verlet --- */
        /* 1. Half-kick velocities */
        for (int a = 0; a < 3; a++)
            for (int i = 1; i < Nx - 1; i++)
                vel[a][i] += 0.5 * dt * acc[a][i];
        for (int i = 1; i < Nx - 1; i++)
            vth[i] += 0.5 * dt * ath[i];

        /* 2. Drift positions */
        for (int a = 0; a < 3; a++)
            for (int i = 1; i < Nx - 1; i++)
                phi[a][i] += dt * vel[a][i];
        for (int i = 1; i < Nx - 1; i++)
            theta[i] += dt * vth[i];

        /* 3. Update P^2 history */
        for (int i = 0; i < Nx; i++)
            P2_prev[i] = P2_cur[i];
        for (int i = 0; i < Nx; i++) {
            double P = phi[0][i] * phi[1][i] * phi[2][i];
            P2_cur[i] = P * P;
        }

        /* 4. Recompute accelerations */
        COMPUTE_ACC();

        /* 5. Half-kick velocities */
        for (int a = 0; a < 3; a++)
            for (int i = 1; i < Nx - 1; i++)
                vel[a][i] += 0.5 * dt * acc[a][i];
        for (int i = 1; i < Nx - 1; i++)
            vth[i] += 0.5 * dt * ath[i];

        /* 6. Absorbing boundary */
        for (int a = 0; a < 3; a++)
            for (int i = 0; i < Nx; i++) {
                vel[a][i] *= damp[i];
                phi[a][i] *= damp[i];
            }
        for (int i = 0; i < Nx; i++) {
            vth[i]   *= damp[i];
            theta[i] *= damp[i];
        }
    }

    fclose(fts);

    /* DFT of phi_1(0,t) and theta(0,t) — second half */
    int dft_start = n_dft / 2;
    if (n_dft - dft_start > 100) {
        char dftpath[600];
        if (tag[0])
            snprintf(dftpath, sizeof(dftpath), "%s/mf_mode%d_%s_spectrum.tsv", outdir, mode, tag);
        else
            snprintf(dftpath, sizeof(dftpath), "%s/mf_mode%d_g%.4f_spectrum.tsv", outdir, mode, g_coup);
        FILE *fdft = fopen(dftpath, "w");
        fprintf(fdft, "omega\tpower_phi\tpower_theta\n");
        double T = t_hist[n_dft-1] - t_hist[dft_start];
        int nf = 500;
        double peak_pow = 0, peak_om = 0;
        double peak_th_pow = 0, peak_th_om = 0;
        for (int k = 0; k < nf; k++) {
            double omega = 3.0 * mass * k / nf;
            double re_p = 0, im_p = 0, re_t = 0, im_t = 0;
            for (int j = dft_start; j < n_dft; j++) {
                double dtj = (j > dft_start) ?
                    (t_hist[j]-t_hist[j-1]) : (t_hist[dft_start+1]-t_hist[dft_start]);
                double c = cos(omega * t_hist[j]);
                double s = sin(omega * t_hist[j]);
                re_p += phi0_hist[j] * c * dtj;
                im_p += phi0_hist[j] * s * dtj;
                re_t += th0_hist[j] * c * dtj;
                im_t += th0_hist[j] * s * dtj;
            }
            double pw_p = (re_p*re_p + im_p*im_p) / (T*T);
            double pw_t = (re_t*re_t + im_t*im_t) / (T*T);
            fprintf(fdft, "%.6f\t%.6e\t%.6e\n", omega, pw_p, pw_t);
            if (pw_p > peak_pow) { peak_pow = pw_p; peak_om = omega; }
            if (pw_t > peak_th_pow) { peak_th_pow = pw_t; peak_th_om = omega; }
        }
        fclose(fdft);
        printf("\nPhi spectrum: peak omega = %.4f (mass gap = %.4f)\n", peak_om, mass);
        printf("Theta spectrum: peak omega = %.4f\n", peak_th_om);
        printf("Oscillon (omega < m)? %s\n",
               (peak_om > 0.01 && peak_om < mass) ? "YES" : "NO");
    }

    /* Final theta profile analysis (single oscillon only) */
    if (mode == 1 || mode == 3) {
        printf("\n=== THETA PROFILE ANALYSIS ===\n");
        printf("  x\t\ttheta\t\t|theta|*x (1D: should be const for 'linear')\n");

        /* In 1D, a massless field sourced by a localized source gives
         * theta ~ const * |x| (growing linearly, NOT decaying).
         * Actually: d^2 theta/dx^2 = -S(x) gives theta ~ -Q/2 * |x| + const
         * where Q = integral of source. But we have a WAVE equation, not Poisson.
         * The retarded solution in 1D: theta(x,t) = (g/2) integral from 0 to t of
         * integral S(x',t-|x-x'|) dx' dt'. For oscillating source at frequency 2*omega:
         * theta ~ amplitude * cos(2*omega*(t - |x|/c)) / 1 (1D: NO geometric falloff).
         * So theta does NOT decay in 1D — it's a traveling wave with constant amplitude.
         */
        /* Measure theta envelope at various distances */
        int test_dists[] = {5, 10, 15, 20, 30, 40, 50, 60, 70};
        int n_test = sizeof(test_dists) / sizeof(test_dists[0]);
        for (int k = 0; k < n_test; k++) {
            int idx = ic + (int)(test_dists[k] / dx);
            if (idx >= 0 && idx < Nx) {
                printf("  x=%3d\t\ttheta=%.6e\n", test_dists[k], theta[idx]);
            }
        }

        /* Time-averaged theta envelope from late-time profiles */
        printf("\n  Note: In 1D, massless radiation has CONSTANT amplitude (no geometric decay).\n");
        printf("  The theta 'profile' is a traveling wave, not a static potential.\n");
        printf("  For long-range force, need the DC (time-averaged) component.\n");
    }

    printf("\nOutput: %s\n", tspath);

    for (int a = 0; a < 3; a++) { free(phi[a]); free(vel[a]); free(acc[a]); }
    free(theta); free(vth); free(ath);
    free(P2_cur); free(P2_prev);
    free(damp); free(phi0_hist); free(th0_hist); free(t_hist);
    return 0;
}
