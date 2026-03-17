/*
 * locality_ema.c — Causal gravity via EMA-sourced wave equation
 *
 * Three-field oscillon (triad) with gravitational potential Phi:
 *   □Phi = alpha * <rho>   (massless wave equation, NOT Poisson)
 *   <rho> = EMA of energy density with time constant tau
 *
 * Backreaction: m²_eff = m²(1+2Phi), c²_eff = 1+4Phi
 *
 * Phase 1: tau scan {2,4,7,10,15,20,30}, measure Phi oscillation
 * Phase 2: boost oscillon, measure Phi response delay at x=50
 *
 * Compile: gcc -O3 -Wall -o locality_ema src/locality_ema.c -lm
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
static double alpha  = -1e-4;   /* gravitational coupling (negative = attractive) */
static double tau    = 7.0;     /* EMA time constant */
static double c_phi  = 1.0;    /* Phi propagation speed */

/* Grid */
static int    Nx     = 8000;
static double xmax   = 200.0;
static double t_equil = 5000.0;
static double t_test  = 5000.0;

/* Mode: 1 = tau scan, 2 = causality test */
static int    mode   = 1;
static double boost_v = 0.01;   /* boost velocity for causality test */
static double x_probe = 50.0;   /* probe distance for causality test */

static char   outdir[512] = "data";

static void parse_args(int argc, char **argv)
{
    for (int i = 1; i < argc - 1; i += 2) {
        if      (!strcmp(argv[i], "-mu"))      mu       = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-kappa"))   kappa    = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-mass"))    mass     = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-A"))       A_init   = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-sigma"))   sigma    = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-alpha"))   alpha    = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-tau"))     tau      = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-cphi"))    c_phi    = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-Nx"))      Nx       = atoi(argv[i+1]);
        else if (!strcmp(argv[i], "-xmax"))    xmax     = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-tequil"))  t_equil  = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-ttest"))   t_test   = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-mode"))    mode     = atoi(argv[i+1]);
        else if (!strcmp(argv[i], "-boost"))   boost_v  = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-xprobe")) x_probe  = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-o"))       strncpy(outdir, argv[i+1], sizeof(outdir)-1);
        else { fprintf(stderr, "Unknown flag: %s\n", argv[i]); exit(1); }
    }
}

/* -dV/dphi_a where V = (mu/2)P^2/(1+kappa P^2), P = phi1 phi2 phi3 */
static double force_pot(double p1, double p2, double p3, int a)
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

/* Compute total energy density at point i */
static double energy_density(double *phi[3], double *vel[3], double dx, double m2,
                             int i, double Phi_i)
{
    double e = 0.0;
    double c2_eff = 1.0 + 4.0 * Phi_i;
    double m2_eff = m2 * (1.0 + 2.0 * Phi_i);
    if (c2_eff < 0.1) c2_eff = 0.1;  /* safety floor */
    if (m2_eff < 0.0) m2_eff = 0.0;

    for (int a = 0; a < 3; a++) {
        e += 0.5 * vel[a][i] * vel[a][i];
        double dp = (phi[a][i+1] - phi[a][i-1]) / (2.0 * dx);
        e += 0.5 * c2_eff * dp * dp;
        e += 0.5 * m2_eff * phi[a][i] * phi[a][i];
    }
    double P = phi[0][i] * phi[1][i] * phi[2][i];
    double P2 = P * P;
    e += 0.5 * mu * P2 / (1.0 + kappa * P2);
    return e;
}

/* Find centroid of energy density (center of mass) */
static double find_centroid(double *phi[3], double *vel[3], double *Phi_field,
                            double dx, double m2, int Nx_loc, double xmax_loc)
{
    double num = 0.0, den = 0.0;
    for (int i = 1; i < Nx_loc - 1; i++) {
        double x = -xmax_loc + i * dx;
        double e = energy_density(phi, vel, dx, m2, i, Phi_field[i]);
        if (e > 0) {
            num += x * e * dx;
            den += e * dx;
        }
    }
    return (den > 1e-20) ? num / den : 0.0;
}

int main(int argc, char **argv)
{
    parse_args(argc, argv);

    double dx  = 2.0 * xmax / (Nx - 1);
    double dx2 = dx * dx;
    double m2  = mass * mass;

    /* CFL: must satisfy for both scalar fields and Phi */
    double kmax = M_PI / dx;
    double dt = 0.4 * 2.0 / sqrt(kmax * kmax + m2);
    /* Also check CFL for Phi wave: dt < dx / c_phi */
    double dt_phi = 0.4 * dx / c_phi;
    if (dt_phi < dt) dt = dt_phi;

    int Nt_equil = (int)(t_equil / dt) + 1;
    int Nt_test  = (int)(t_test / dt) + 1;
    (void)(Nt_equil + Nt_test);  /* total steps for reference */

    int ic = Nx / 2;  /* center index */

    printf("locality_ema: mode=%d, alpha=%.2e, tau=%.1f, c_phi=%.2f\n",
           mode, alpha, tau, c_phi);
    printf("  mu=%.1f kappa=%.1f mass=%.3f A=%.3f sigma=%.3f\n",
           mu, kappa, mass, A_init, sigma);
    printf("  Nx=%d xmax=%.1f dx=%.5f dt=%.6f\n", Nx, xmax, dx, dt);
    printf("  t_equil=%.0f (%d steps), t_test=%.0f (%d steps)\n",
           t_equil, Nt_equil, t_test, Nt_test);

    /* Allocate scalar fields */
    double *phi[3], *vel[3], *acc[3];
    for (int a = 0; a < 3; a++) {
        phi[a] = calloc(Nx, sizeof(double));
        vel[a] = calloc(Nx, sizeof(double));
        acc[a] = calloc(Nx, sizeof(double));
    }

    /* Allocate Phi field (gravitational potential) and its velocity/acceleration */
    double *Phi      = calloc(Nx, sizeof(double));
    double *Phi_vel  = calloc(Nx, sizeof(double));
    double *Phi_acc  = calloc(Nx, sizeof(double));

    /* EMA of energy density */
    double *rho_ema  = calloc(Nx, sizeof(double));

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

    /* Initialize: symmetric triad Gaussians at origin */
    for (int a = 0; a < 3; a++)
        for (int i = 0; i < Nx; i++) {
            double x = -xmax + i * dx;
            phi[a][i] = A_init * exp(-x * x / (2.0 * sigma * sigma));
        }

    /* Macro: compute acceleration for scalar fields with backreaction from Phi */
    #define COMPUTE_ACC_SCALAR() do { \
        for (int a = 0; a < 3; a++) { \
            acc[a][0] = acc[a][1] = acc[a][Nx-2] = acc[a][Nx-1] = 0; \
            for (int i = 1; i < Nx - 1; i++) { \
                double c2_eff = 1.0 + 4.0 * Phi[i]; \
                double m2_eff = m2 * (1.0 + 2.0 * Phi[i]); \
                if (c2_eff < 0.1) c2_eff = 0.1; \
                double lapl = (phi[a][i+1] - 2.0*phi[a][i] + phi[a][i-1]) / dx2; \
                double fp = force_pot(phi[0][i], phi[1][i], phi[2][i], a); \
                acc[a][i] = c2_eff * lapl - m2_eff * phi[a][i] + fp; \
            } \
        } \
    } while(0)

    /* Macro: compute acceleration for Phi (wave equation) */
    #define COMPUTE_ACC_PHI() do { \
        Phi_acc[0] = Phi_acc[Nx-1] = 0; \
        for (int i = 1; i < Nx - 1; i++) { \
            double lapl = (Phi[i+1] - 2.0*Phi[i] + Phi[i-1]) / dx2; \
            Phi_acc[i] = c_phi * c_phi * lapl + alpha * rho_ema[i]; \
        } \
    } while(0)

    /* Initial acceleration (no Phi yet during equilibration) */
    COMPUTE_ACC_SCALAR();

    /* Phase 1: equilibrate oscillon without gravity */
    printf("\n--- Equilibrating (no gravity) for t=%.0f ---\n", t_equil);
    int print_every_eq = Nt_equil / 20;
    if (print_every_eq < 1) print_every_eq = 1;

    for (int n = 0; n < Nt_equil; n++) {
        double t = n * dt;

        /* Initialize EMA with current density during equilibration */
        for (int i = 1; i < Nx - 1; i++) {
            double rho = energy_density(phi, vel, dx, m2, i, 0.0);
            double w = dt / tau;
            if (w > 1.0) w = 1.0;
            rho_ema[i] = (1.0 - w) * rho_ema[i] + w * rho;
        }

        if (n % print_every_eq == 0) {
            double fc_num = 0, fc_den = 0;
            for (int i = 1; i < Nx - 1; i++) {
                double x = -xmax + i * dx;
                double e = energy_density(phi, vel, dx, m2, i, 0.0);
                fc_den += e * dx;
                if (fabs(x) < 3.0 * sigma) fc_num += e * dx;
            }
            printf("  t=%7.1f  phi0=%.4f  fc=%.4f  rho_ema(0)=%.4e\n",
                   t, phi[0][ic], fc_num/(fc_den+1e-30), rho_ema[ic]);
        }

        /* Velocity Verlet for scalars only */
        for (int a = 0; a < 3; a++)
            for (int i = 1; i < Nx - 1; i++)
                vel[a][i] += 0.5 * dt * acc[a][i];
        for (int a = 0; a < 3; a++)
            for (int i = 1; i < Nx - 1; i++)
                phi[a][i] += dt * vel[a][i];
        COMPUTE_ACC_SCALAR();
        for (int a = 0; a < 3; a++)
            for (int i = 1; i < Nx - 1; i++)
                vel[a][i] += 0.5 * dt * acc[a][i];

        /* Absorbing BC for scalars */
        for (int a = 0; a < 3; a++)
            for (int i = 0; i < Nx; i++) {
                vel[a][i] *= damp[i];
                phi[a][i] *= damp[i];
            }
    }

    printf("Equilibration done. phi0(center)=%.6f, rho_ema(center)=%.6e\n",
           phi[0][ic], rho_ema[ic]);

    /* For causality test mode 2: Phi starts at zero, source turns on at t=0.
     * The Phi wave should arrive at x_probe after delay = x_probe/c_phi.
     * Also apply boost to oscillon to create a moving-source signal. */
    if (mode == 2) {
        printf("\n--- Causality test: Phi=0, source ON, boost v=%.4f ---\n", boost_v);
        /* Reset Phi to zero */
        memset(Phi, 0, Nx * sizeof(double));
        memset(Phi_vel, 0, Nx * sizeof(double));
        memset(Phi_acc, 0, Nx * sizeof(double));
        /* Apply boost to oscillon */
        for (int a = 0; a < 3; a++)
            for (int i = 1; i < Nx - 1; i++) {
                double dphi_dx = (phi[a][i+1] - phi[a][i-1]) / (2.0 * dx);
                vel[a][i] += -boost_v * dphi_dx;
            }
    }

    /* Phase 2: test phase with gravity ON */
    printf("\n--- Test phase (gravity ON) for t=%.0f ---\n", t_test);

    /* Initialize Phi acceleration */
    COMPUTE_ACC_PHI();

    /* Output file */
    char tspath[600];
    if (mode == 1)
        snprintf(tspath, sizeof(tspath), "%s/ema_tau%.0f_ts.tsv", outdir, tau);
    else
        snprintf(tspath, sizeof(tspath), "%s/ema_causality_tau%.0f_ts.tsv", outdir, tau);

    FILE *fts = fopen(tspath, "w");
    if (!fts) { fprintf(stderr, "Cannot open %s\n", tspath); return 1; }

    if (mode == 1)
        fprintf(fts, "time\tphi0_center\tPhi_center\tPhi_min\tPhi_max\t"
                     "rho_ema_center\tfc\tE_total\tomega_inst\n");
    else
        fprintf(fts, "time\tphi0_center\tPhi_center\tPhi_at_probe\t"
                     "centroid\trho_ema_center\tfc\n");

    int rec_every = Nt_test / 20000;
    if (rec_every < 1) rec_every = 1;
    int print_every = Nt_test / 40;
    if (print_every < 1) print_every = 1;

    /* DFT storage for Phi oscillation analysis */
    int max_dft = 50000;
    double *Phi0_hist = malloc(max_dft * sizeof(double));
    double *phi0_hist = malloc(max_dft * sizeof(double));
    double *t_hist    = malloc(max_dft * sizeof(double));
    int n_dft = 0;
    int dft_every = Nt_test / max_dft;
    if (dft_every < 1) dft_every = 1;

    /* Probe index for causality test */
    int i_probe = ic + (int)(x_probe / dx);
    if (i_probe >= Nx - 1) i_probe = Nx - 2;

    /* Track Phi statistics */
    double Phi_center_sum = 0, Phi_center_sq_sum = 0;
    int Phi_count = 0;
    /* Track initial Phi at probe for causality */
    double Phi_probe_init = 0.0;
    int probe_init_set = 0;
    double delay_time = -1.0;

    double prev_phi0 = phi[0][ic];
    double prev_t = 0;

    for (int n = 0; n <= Nt_test; n++) {
        double t = n * dt;

        /* Update EMA */
        for (int i = 1; i < Nx - 1; i++) {
            double rho = energy_density(phi, vel, dx, m2, i, Phi[i]);
            double w = dt / tau;
            if (w > 1.0) w = 1.0;
            rho_ema[i] = (1.0 - w) * rho_ema[i] + w * rho;
        }

        /* DFT recording */
        if (n % dft_every == 0 && n_dft < max_dft) {
            Phi0_hist[n_dft] = Phi[ic];
            phi0_hist[n_dft] = phi[0][ic];
            t_hist[n_dft] = t;
            n_dft++;
        }

        /* Phi statistics (skip first 10% for settling) */
        if (n > Nt_test / 10) {
            Phi_center_sum += Phi[ic];
            Phi_center_sq_sum += Phi[ic] * Phi[ic];
            Phi_count++;
        }

        /* Causality: detect when Phi at probe changes */
        if (mode == 2) {
            if (!probe_init_set && n > 10) {
                Phi_probe_init = Phi[i_probe];
                probe_init_set = 1;
            }
            if (probe_init_set && delay_time < 0) {
                double delta = fabs(Phi[i_probe] - Phi_probe_init);
                if (delta > 1e-10) {
                    delay_time = t;
                }
            }
        }

        int do_rec = (n % rec_every == 0);
        int do_print = (n % print_every == 0);

        if (do_rec || do_print) {
            /* Compute diagnostics */
            double Etot = 0, fc_num = 0, fc_den = 0;
            double Phi_min_val = 1e30, Phi_max_val = -1e30;
            double core_r = 3.0 * sigma;

            for (int i = 1; i < Nx - 1; i++) {
                double x = -xmax + i * dx;
                double e = energy_density(phi, vel, dx, m2, i, Phi[i]);
                /* Add Phi field energy */
                e += 0.5 * Phi_vel[i] * Phi_vel[i];
                double dPhi = (Phi[i+1] - Phi[i-1]) / (2.0 * dx);
                e += 0.5 * c_phi * c_phi * dPhi * dPhi;

                Etot += e * dx;
                fc_den += e * dx;
                if (fabs(x) < core_r) fc_num += e * dx;

                if (Phi[i] < Phi_min_val) Phi_min_val = Phi[i];
                if (Phi[i] > Phi_max_val) Phi_max_val = Phi[i];
            }

            double fc = (fc_den > 1e-20) ? fc_num / fc_den : 0.0;

            /* Instantaneous frequency from zero crossings */
            double omega_inst = 0;
            if (phi[0][ic] * prev_phi0 < 0 && t > 0) {
                omega_inst = M_PI / (t - prev_t + 1e-30);
                prev_t = t;
            }
            if (phi[0][ic] * prev_phi0 < 0) prev_t = t;
            prev_phi0 = phi[0][ic];

            if (do_rec) {
                if (mode == 1) {
                    fprintf(fts, "%.6f\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6f\t%.6e\t%.4f\n",
                            t, phi[0][ic], Phi[ic], Phi_min_val, Phi_max_val,
                            rho_ema[ic], fc, Etot, omega_inst);
                } else {
                    double centroid = find_centroid(phi, vel, Phi, dx, m2, Nx, xmax);
                    fprintf(fts, "%.6f\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6f\n",
                            t, phi[0][ic], Phi[ic], Phi[i_probe], centroid,
                            rho_ema[ic], fc);
                }
            }

            if (do_print) {
                if (mode == 1)
                    printf("  t=%7.1f  phi0=%+.4f  Phi(0)=%+.4e  Phi_range=[%+.3e,%+.3e]  "
                           "rho_ema=%.3e  fc=%.3f\n",
                           t, phi[0][ic], Phi[ic], Phi_min_val, Phi_max_val,
                           rho_ema[ic], fc);
                else {
                    double centroid = find_centroid(phi, vel, Phi, dx, m2, Nx, xmax);
                    printf("  t=%7.1f  phi0=%+.4f  Phi(0)=%+.4e  Phi(probe)=%+.4e  "
                           "x_cm=%.3f  fc=%.3f\n",
                           t, phi[0][ic], Phi[ic], Phi[i_probe], centroid, fc);
                }
            }
        }

        if (n == Nt_test) break;

        /* Velocity Verlet: half-kick */
        for (int a = 0; a < 3; a++)
            for (int i = 1; i < Nx - 1; i++)
                vel[a][i] += 0.5 * dt * acc[a][i];
        for (int i = 1; i < Nx - 1; i++)
            Phi_vel[i] += 0.5 * dt * Phi_acc[i];

        /* Drift */
        for (int a = 0; a < 3; a++)
            for (int i = 1; i < Nx - 1; i++)
                phi[a][i] += dt * vel[a][i];
        for (int i = 1; i < Nx - 1; i++)
            Phi[i] += dt * Phi_vel[i];

        /* New accelerations */
        COMPUTE_ACC_SCALAR();
        COMPUTE_ACC_PHI();

        /* Half-kick */
        for (int a = 0; a < 3; a++)
            for (int i = 1; i < Nx - 1; i++)
                vel[a][i] += 0.5 * dt * acc[a][i];
        for (int i = 1; i < Nx - 1; i++)
            Phi_vel[i] += 0.5 * dt * Phi_acc[i];

        /* Absorbing BC for all fields */
        for (int a = 0; a < 3; a++)
            for (int i = 0; i < Nx; i++) {
                vel[a][i] *= damp[i];
                phi[a][i] *= damp[i];
            }
        for (int i = 0; i < Nx; i++) {
            Phi_vel[i] *= damp[i];
            Phi[i]     *= damp[i];
        }
    }

    fclose(fts);

    /* Analysis */
    double Phi_mean = Phi_center_sum / (Phi_count + 1);
    double Phi_rms  = sqrt(Phi_center_sq_sum / (Phi_count + 1) - Phi_mean * Phi_mean);

    printf("\n=== RESULTS (tau=%.1f) ===\n", tau);
    printf("  Phi(0) mean  = %+.6e\n", Phi_mean);
    printf("  Phi(0) rms   = %.6e (oscillation amplitude)\n", Phi_rms);
    printf("  Phi osc/mean = %.4f\n", fabs(Phi_rms / (Phi_mean + 1e-30)));

    /* DFT of Phi(0) to check oscillation frequency */
    if (n_dft > 200) {
        int dft_start = n_dft / 2;
        double T_span = t_hist[n_dft-1] - t_hist[dft_start];
        int nf = 500;
        double peak_pow = 0, peak_om = 0;
        double dc_pow = 0;

        char dftpath[600];
        if (mode == 1)
            snprintf(dftpath, sizeof(dftpath), "%s/ema_tau%.0f_Phi_spectrum.tsv", outdir, tau);
        else
            snprintf(dftpath, sizeof(dftpath), "%s/ema_causality_Phi_spectrum.tsv", outdir);
        FILE *fdft = fopen(dftpath, "w");
        fprintf(fdft, "omega\tPhi_power\tphi_power\n");

        for (int k = 0; k < nf; k++) {
            double omega = 3.0 * mass * k / nf;
            double re_P = 0, im_P = 0, re_p = 0, im_p = 0;
            for (int j = dft_start; j < n_dft; j++) {
                double dtj = (j > dft_start) ?
                    (t_hist[j]-t_hist[j-1]) : (t_hist[dft_start+1]-t_hist[dft_start]);
                re_P += Phi0_hist[j] * cos(omega * t_hist[j]) * dtj;
                im_P += Phi0_hist[j] * sin(omega * t_hist[j]) * dtj;
                re_p += phi0_hist[j] * cos(omega * t_hist[j]) * dtj;
                im_p += phi0_hist[j] * sin(omega * t_hist[j]) * dtj;
            }
            double pw_P = (re_P*re_P + im_P*im_P) / (T_span*T_span);
            double pw_p = (re_p*re_p + im_p*im_p) / (T_span*T_span);
            fprintf(fdft, "%.6f\t%.6e\t%.6e\n", omega, pw_P, pw_p);
            if (k == 0) dc_pow = pw_P;
            if (k > 0 && pw_P > peak_pow) { peak_pow = pw_P; peak_om = omega; }
        }
        fclose(fdft);

        printf("  Phi spectrum: DC power = %.4e, peak AC at omega=%.3f power=%.4e\n",
               dc_pow, peak_om, peak_pow);
        printf("  AC/DC ratio = %.4e\n", peak_pow / (dc_pow + 1e-30));

        /* Also find phi oscillation frequency */
        double phi_peak_pow = 0, phi_peak_om = 0;
        for (int k = 1; k < nf; k++) {
            double omega = 3.0 * mass * k / nf;
            double re = 0, im = 0;
            for (int j = dft_start; j < n_dft; j++) {
                double dtj = (j > dft_start) ?
                    (t_hist[j]-t_hist[j-1]) : (t_hist[dft_start+1]-t_hist[dft_start]);
                re += phi0_hist[j] * cos(omega * t_hist[j]) * dtj;
                im += phi0_hist[j] * sin(omega * t_hist[j]) * dtj;
            }
            double pw = (re*re + im*im) / (T_span*T_span);
            if (pw > phi_peak_pow) { phi_peak_pow = pw; phi_peak_om = omega; }
        }
        printf("  Oscillon omega = %.4f (mass gap = %.4f)\n", phi_peak_om, mass);
    }

    if (mode == 2) {
        printf("\n=== CAUSALITY TEST ===\n");
        printf("  Probe at x = %.1f\n", x_probe);
        printf("  Expected delay (causal) = %.1f\n", x_probe / c_phi);
        if (delay_time >= 0)
            printf("  Measured delay = %.2f\n", delay_time);
        else
            printf("  No signal detected at probe (threshold too high?)\n");
    }

    printf("\nOutput: %s\n", tspath);

    /* Cleanup */
    for (int a = 0; a < 3; a++) { free(phi[a]); free(vel[a]); free(acc[a]); }
    free(Phi); free(Phi_vel); free(Phi_acc); free(rho_ema);
    free(damp); free(Phi0_hist); free(phi0_hist); free(t_hist);
    return 0;
}
