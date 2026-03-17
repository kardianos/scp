/*
 * locality_retarded.c — Causal gravity via wave equation for oscillon
 *
 * Instead of instantaneous Poisson (d²Phi/dx² = alpha*rho), solve the
 * RETARDED wave equation:
 *
 *   (1/c²) d²Phi/dt² - d²Phi/dx² = -alpha*rho(x,t)
 *
 * This is exactly causal: signals propagate at speed c.
 * No free parameters beyond alpha and c.
 *
 * The wave equation is solved directly as a PDE using leapfrog, alongside
 * the field evolution. For a static source, Phi converges locally to the
 * Poisson solution (plus outgoing transient waves).
 *
 * Phase 1: Single oscillon, verify Phi(0) matches Poisson steady-state.
 * Phase 2: Boost oscillon at t_boost, verify delay at x=50 is exactly 50/c.
 *
 * Compile: gcc -O3 -Wall -o locality_retarded src/locality_retarded.c -lm
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Parameters */
static double mu       = -20.0;
static double kappa    = 20.0;
static double mass     = 1.0;
static double A_init   = 0.8;
static double sigma_init = 3.0;
static int    Nx       = 4000;
static double xmax     = 100.0;
static double t_equil  = 5000.0;
static double t_grav   = 5000.0;
static double alpha    = -1e-4;
static double c_light  = 1.0;
static double t_ramp   = 200.0;  /* adiabatic ramp-up time */
static char   outdir[512] = "v24/fundamental/locality_retarded/data";

/* Phase 2 */
static double t_boost  = 2500.0;
static double v_boost  = 0.01;
static double x_probe  = 50.0;

static void parse_args(int argc, char **argv)
{
    for (int i = 1; i < argc - 1; i += 2) {
        if      (!strcmp(argv[i], "-mu"))      mu       = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-kappa"))   kappa    = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-mass"))    mass     = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-A"))       A_init   = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-sigma"))   sigma_init = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-Nx"))      Nx       = atoi(argv[i+1]);
        else if (!strcmp(argv[i], "-xmax"))    xmax     = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-equil"))   t_equil  = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-tgrav"))   t_grav   = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-alpha"))   alpha    = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-c"))       c_light  = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-tramp"))   t_ramp   = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-tboost"))  t_boost  = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-vboost"))  v_boost  = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-xprobe"))  x_probe  = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-o"))       strncpy(outdir, argv[i+1], sizeof(outdir)-1);
        else { fprintf(stderr, "Unknown flag: %s\n", argv[i]); exit(1); }
    }
}

/* -dV/dphi_a */
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

/* Compute energy density rho(x) */
static void compute_rho(double *rho, double **phi, double **vel, int N, double dx, double m2)
{
    for (int i = 0; i < N; i++) rho[i] = 0.0;
    for (int i = 1; i < N - 1; i++) {
        double e = 0.0;
        for (int a = 0; a < 3; a++) {
            double dp = (phi[a][i+1] - phi[a][i-1]) / (2.0 * dx);
            e += 0.5 * vel[a][i] * vel[a][i];
            e += 0.5 * dp * dp;
            e += 0.5 * m2 * phi[a][i] * phi[a][i];
        }
        double P = phi[0][i] * phi[1][i] * phi[2][i];
        double P2 = P * P;
        e += 0.5 * mu * P2 / (1.0 + kappa * P2);
        rho[i] = e;
    }
}

/* Compute total energy and core fraction */
static void compute_energy(double *Et, double *Ec, double *fc,
                          double **phi, double **vel, int N, double dx,
                          double m2, double core_r, double xmin)
{
    *Et = 0; *Ec = 0;
    for (int i = 1; i < N - 1; i++) {
        double x = xmin + i * dx;
        double e = 0;
        for (int a = 0; a < 3; a++) {
            double dp = (phi[a][i+1] - phi[a][i-1]) / (2.0*dx);
            e += 0.5*vel[a][i]*vel[a][i] + 0.5*dp*dp + 0.5*m2*phi[a][i]*phi[a][i];
        }
        double P = phi[0][i]*phi[1][i]*phi[2][i];
        double P2 = P*P;
        e += 0.5*mu*P2/(1.0+kappa*P2);
        *Et += e * dx;
        if (fabs(x) < core_r) *Ec += e * dx;
    }
    *fc = (*Et > 1e-20) ? *Ec / *Et : 0.0;
}

/* Solve Poisson for comparison */
static void solve_poisson(double *Phi, const double *rho, int N, double dx, double alp)
{
    double *F = malloc(N * sizeof(double));
    F[0] = 0.0;
    for (int i = 1; i < N; i++)
        F[i] = F[i-1] + alp * rho[i] * dx;
    Phi[0] = 0.0;
    for (int i = 1; i < N; i++)
        Phi[i] = Phi[i-1] + F[i] * dx;
    double L = (N - 1) * dx;
    if (L > 0) {
        double slope = Phi[N-1] / L;
        for (int i = 0; i < N; i++)
            Phi[i] -= slope * (i * dx);
    }
    free(F);
}


int main(int argc, char **argv)
{
    parse_args(argc, argv);

    double dx  = 2.0 * xmax / (Nx - 1);
    double dx2 = dx * dx;
    double m2  = mass * mass;

    /* CFL: need stability for both KG fields AND wave equation for Phi.
     * KG: dt < 2/sqrt(kmax^2 + m^2)
     * Wave: dt < dx/c
     * Use the more restrictive. */
    double kmax = M_PI / dx;
    double dt_kg   = 0.4 * 2.0 / sqrt(kmax * kmax + m2);
    double dt_wave = 0.4 * dx / c_light;
    double dt = (dt_kg < dt_wave) ? dt_kg : dt_wave;

    int Nt_equil = (int)(t_equil / dt) + 1;
    int Nt_grav  = (int)(t_grav / dt) + 1;
    int Nt_ramp  = (int)(t_ramp / dt) + 1;

    printf("locality_retarded: Causal gravity via wave equation\n");
    printf("  mu=%.1f kappa=%.1f mass=%.4f A=%.3f sigma=%.3f\n",
           mu, kappa, mass, A_init, sigma_init);
    printf("  Nx=%d xmax=%.1f dx=%.5f dt=%.6f\n", Nx, xmax, dx, dt);
    printf("  alpha=%.6e c=%.4f t_ramp=%.1f\n", alpha, c_light, t_ramp);
    printf("  t_equil=%.0f t_grav=%.0f\n", t_equil, t_grav);
    printf("  t_boost=%.0f v_boost=%.4f x_probe=%.1f\n", t_boost, v_boost, x_probe);
    printf("  Nt_equil=%d Nt_grav=%d\n", Nt_equil, Nt_grav);

    /* Allocate fields */
    double *phi[3], *vel[3], *acc[3];
    for (int a = 0; a < 3; a++) {
        phi[a] = calloc(Nx, sizeof(double));
        vel[a] = calloc(Nx, sizeof(double));
        acc[a] = calloc(Nx, sizeof(double));
    }

    /* Gravitational potential: wave equation fields */
    double *Phi     = calloc(Nx, sizeof(double));  /* potential */
    double *Phi_vel = calloc(Nx, sizeof(double));  /* dPhi/dt */
    double *Phi_acc = calloc(Nx, sizeof(double));  /* d²Phi/dt² */
    double *Phi_poi = calloc(Nx, sizeof(double));  /* Poisson comparison */
    double *rho     = calloc(Nx, sizeof(double));
    double *rho_avg = calloc(Nx, sizeof(double));
    double *dPhi    = calloc(Nx, sizeof(double));  /* dPhi/dx for field evolution */

    /* Absorbing boundary for Phi: damp outgoing waves */
    double *damp_phi = malloc(Nx * sizeof(double));

    /* Absorbing boundary for fields: outer 25% */
    double *damp_arr = malloc(Nx * sizeof(double));
    double x_abs = xmax * 0.75;
    for (int i = 0; i < Nx; i++) {
        double x = -xmax + i * dx;
        double ax = fabs(x);
        if (ax > x_abs) {
            double f = (ax - x_abs) / (xmax - x_abs);
            damp_arr[i] = 1.0 - 0.98 * f * f;
            damp_phi[i] = 1.0 - 0.999 * f * f;  /* stronger damping for Phi */
        } else {
            damp_arr[i] = 1.0;
            damp_phi[i] = 1.0;
        }
    }

    /* Initialize: Gaussians (symmetric triad) */
    for (int a = 0; a < 3; a++)
        for (int i = 0; i < Nx; i++) {
            double x = -xmax + i * dx;
            phi[a][i] = A_init * exp(-x * x / (2.0 * sigma_init * sigma_init));
        }

    int ic = Nx / 2;
    double core_r = 3.0 * sigma_init;

    /* Flat-metric acceleration */
    #define COMPUTE_ACC_FLAT() do { \
        for (int a = 0; a < 3; a++) { \
            acc[a][0] = acc[a][1] = acc[a][Nx-2] = acc[a][Nx-1] = 0; \
            for (int i = 1; i < Nx - 1; i++) { \
                double lapl = (phi[a][i+1] - 2.0*phi[a][i] + phi[a][i-1]) / dx2; \
                double fp = force_pot(phi[0][i], phi[1][i], phi[2][i], a); \
                acc[a][i] = lapl - m2*phi[a][i] + fp; \
            } \
        } \
    } while(0)

    /* Curved-metric acceleration */
    #define COMPUTE_ACC_GRAV() do { \
        for (int a = 0; a < 3; a++) { \
            acc[a][0] = acc[a][1] = acc[a][Nx-2] = acc[a][Nx-1] = 0; \
            for (int i = 1; i < Nx - 1; i++) { \
                double lapl = (phi[a][i+1] - 2.0*phi[a][i] + phi[a][i-1]) / dx2; \
                double fp = force_pot(phi[0][i], phi[1][i], phi[2][i], a); \
                double dphi_dx = (phi[a][i+1] - phi[a][i-1]) / (2.0 * dx); \
                acc[a][i] = (1.0 + 4.0*Phi[i]) * lapl \
                          + 2.0 * dPhi[i] * dphi_dx \
                          - m2 * (1.0 + 2.0*Phi[i]) * phi[a][i] \
                          + fp; \
            } \
        } \
    } while(0)

    /* Wave equation acceleration for Phi:
     * (1/c²) Phi_tt = Phi_xx - alpha_eff * rho
     * => Phi_tt = c² * Phi_xx - alpha_eff * c² * rho
     */
    #define COMPUTE_PHI_ACC(alpha_eff) do { \
        double c2 = c_light * c_light; \
        Phi_acc[0] = Phi_acc[Nx-1] = 0; \
        for (int i = 1; i < Nx - 1; i++) { \
            double lapl_phi = (Phi[i+1] - 2.0*Phi[i] + Phi[i-1]) / dx2; \
            Phi_acc[i] = c2 * lapl_phi - (alpha_eff) * c2 * rho[i]; \
        } \
    } while(0)

    /* === Phase 0: Equilibrate on flat metric === */
    printf("\n=== Phase 0: Equilibrate (flat metric, t=0..%.0f) ===\n", t_equil);

    COMPUTE_ACC_FLAT();

    int print_every = Nt_equil / 20;
    if (print_every < 1) print_every = 1;

    for (int n = 0; n < Nt_equil; n++) {
        if (n % print_every == 0) {
            double Et, Ec, fc;
            compute_energy(&Et, &Ec, &fc, phi, vel, Nx, dx, m2, core_r, -xmax);
            printf("  t=%7.1f  phi0=(%+.4f,%+.4f,%+.4f)  E=%.4f  fc=%.3f\n",
                   n * dt, phi[0][ic], phi[1][ic], phi[2][ic], Et, fc);
        }

        /* Velocity Verlet for fields */
        for (int a = 0; a < 3; a++)
            for (int i = 1; i < Nx - 1; i++)
                vel[a][i] += 0.5 * dt * acc[a][i];
        for (int a = 0; a < 3; a++)
            for (int i = 1; i < Nx - 1; i++)
                phi[a][i] += dt * vel[a][i];
        COMPUTE_ACC_FLAT();
        for (int a = 0; a < 3; a++)
            for (int i = 1; i < Nx - 1; i++)
                vel[a][i] += 0.5 * dt * acc[a][i];

        for (int a = 0; a < 3; a++)
            for (int i = 0; i < Nx; i++) {
                vel[a][i] *= damp_arr[i];
                phi[a][i] *= damp_arr[i];
            }
    }

    double E_equil, Ec_equil, fc_equil;
    compute_energy(&E_equil, &Ec_equil, &fc_equil, phi, vel, Nx, dx, m2, core_r, -xmax);
    printf("\n  Equilibrated: E=%.4f  fc=%.4f  phi0=%.4f\n", E_equil, fc_equil, phi[0][ic]);

    /* Compute time-averaged rho for Poisson comparison */
    printf("\n  Computing time-averaged rho...\n");
    double T_osc = 8.0;
    int N_avg = (int)(T_osc / dt);
    for (int i = 0; i < Nx; i++) rho_avg[i] = 0.0;

    double *phi_tmp[3], *vel_tmp[3], *acc_tmp[3];
    for (int a = 0; a < 3; a++) {
        phi_tmp[a] = malloc(Nx * sizeof(double));
        vel_tmp[a] = malloc(Nx * sizeof(double));
        acc_tmp[a] = malloc(Nx * sizeof(double));
        memcpy(phi_tmp[a], phi[a], Nx * sizeof(double));
        memcpy(vel_tmp[a], vel[a], Nx * sizeof(double));
    }
    for (int a = 0; a < 3; a++) {
        acc_tmp[a][0] = acc_tmp[a][1] = acc_tmp[a][Nx-2] = acc_tmp[a][Nx-1] = 0;
        for (int i = 1; i < Nx - 1; i++) {
            double lapl = (phi_tmp[a][i+1] - 2.0*phi_tmp[a][i] + phi_tmp[a][i-1]) / dx2;
            double fp = force_pot(phi_tmp[0][i], phi_tmp[1][i], phi_tmp[2][i], a);
            acc_tmp[a][i] = lapl - m2*phi_tmp[a][i] + fp;
        }
    }
    for (int n = 0; n < N_avg; n++) {
        compute_rho(rho, phi_tmp, vel_tmp, Nx, dx, m2);
        for (int i = 0; i < Nx; i++) rho_avg[i] += rho[i] / N_avg;

        for (int a = 0; a < 3; a++)
            for (int i = 1; i < Nx - 1; i++)
                vel_tmp[a][i] += 0.5 * dt * acc_tmp[a][i];
        for (int a = 0; a < 3; a++)
            for (int i = 1; i < Nx - 1; i++)
                phi_tmp[a][i] += dt * vel_tmp[a][i];
        for (int a = 0; a < 3; a++) {
            acc_tmp[a][0] = acc_tmp[a][1] = acc_tmp[a][Nx-2] = acc_tmp[a][Nx-1] = 0;
            for (int i = 1; i < Nx - 1; i++) {
                double lapl = (phi_tmp[a][i+1] - 2.0*phi_tmp[a][i] + phi_tmp[a][i-1]) / dx2;
                double fp = force_pot(phi_tmp[0][i], phi_tmp[1][i], phi_tmp[2][i], a);
                acc_tmp[a][i] = lapl - m2*phi_tmp[a][i] + fp;
            }
        }
        for (int a = 0; a < 3; a++)
            for (int i = 1; i < Nx - 1; i++)
                vel_tmp[a][i] += 0.5 * dt * acc_tmp[a][i];
        for (int a = 0; a < 3; a++)
            for (int i = 0; i < Nx; i++) {
                vel_tmp[a][i] *= damp_arr[i];
                phi_tmp[a][i] *= damp_arr[i];
            }
    }
    for (int a = 0; a < 3; a++) { free(phi_tmp[a]); free(vel_tmp[a]); free(acc_tmp[a]); }

    double rho_peak_avg = 0;
    for (int i = 0; i < Nx; i++)
        if (rho_avg[i] > rho_peak_avg) rho_peak_avg = rho_avg[i];
    printf("  <rho> peak = %.6f\n", rho_peak_avg);

    solve_poisson(Phi_poi, rho_avg, Nx, dx, alpha);
    printf("  Poisson Phi(0) = %.6e\n", Phi_poi[ic]);

    /* === Phase 1+2: Wave equation gravity with boost === */
    printf("\n=== Phase 1+2: Retarded gravity via wave eq (t=%.0f..%.0f, boost at +%.0f) ===\n",
           t_equil, t_equil + t_grav, t_boost);

    /* Initialize Phi to Poisson solution (skip transient ringing) */
    /* NO -- start from zero to test causality properly. Use adiabatic ramp. */
    memset(Phi, 0, Nx * sizeof(double));
    memset(Phi_vel, 0, Nx * sizeof(double));
    memset(Phi_acc, 0, Nx * sizeof(double));

    COMPUTE_ACC_FLAT();
    compute_rho(rho, phi, vel, Nx, dx, m2);
    COMPUTE_PHI_ACC(0.0);

    /* Time series */
    char tspath[600];
    snprintf(tspath, sizeof(tspath), "%s/retarded_ts.tsv", outdir);
    FILE *fts = fopen(tspath, "w");
    if (!fts) { fprintf(stderr, "Cannot open %s\n", tspath); return 1; }
    fprintf(fts, "time\tphi1_0\tE_total\tfc\tPhi_wave_0\tPhi_poi_0\talpha_eff\n");

    /* Probe */
    char probepath[600];
    snprintf(probepath, sizeof(probepath), "%s/retarded_probe.tsv", outdir);
    FILE *fprobe = fopen(probepath, "w");
    if (!fprobe) { fprintf(stderr, "Cannot open %s\n", probepath); return 1; }
    fprintf(fprobe, "time\tPhi_probe\tPhi_origin\tPhi_probe_ema\trho_origin\n");

    int print_every2 = Nt_grav / 60;
    if (print_every2 < 1) print_every2 = 1;
    int rec_every = Nt_grav / 20000;
    if (rec_every < 1) rec_every = 1;
    int probe_every = (int)(1.0 / dt);
    if (probe_every < 1) probe_every = 1;

    int boosted = 0;
    int i_probe = (int)((x_probe + xmax) / dx);
    if (i_probe < 0) i_probe = 0;
    if (i_probe >= Nx) i_probe = Nx - 1;

    /* Track pre-boost Phi at probe for causality measurement.
     * Use a running average to filter out oscillations from the static oscillon.
     * Detect when the running average deviates from its pre-boost value. */
    double Phi_probe_avg = 0.0;       /* EMA of Phi at probe */
    double Phi_probe_pre_boost = 0.0; /* baseline EMA before boost */
    double Phi_probe_max_osc = 0.0;   /* max oscillation amplitude (pre-boost) */
    int measuring_delay = 0;
    int calibrating = 0;       /* 1 = recording pre-boost oscillation stats */
    double delay_arrival = -1.0;
    double t_calib_start = 0.0;
    double ema_tau = 30.0;     /* EMA timescale: several oscillation periods */

    /* DFT storage */
    int max_dft = 50000;
    double *phi0_hist = malloc(max_dft * sizeof(double));
    double *t_hist = malloc(max_dft * sizeof(double));
    int n_dft = 0;
    int dft_every = Nt_grav / max_dft;
    if (dft_every < 1) dft_every = 1;

    for (int n = 0; n <= Nt_grav; n++) {
        double t = t_equil + n * dt;
        double T_elapsed = n * dt;

        /* Adiabatic ramp */
        double ramp = (n < Nt_ramp) ? (double)n / Nt_ramp : 1.0;
        ramp = ramp * ramp * (3.0 - 2.0 * ramp);  /* Hermite smoothstep */
        double alpha_eff = alpha * ramp;

        /* DFT */
        if (n % dft_every == 0 && n_dft < max_dft) {
            phi0_hist[n_dft] = phi[0][ic];
            t_hist[n_dft] = t;
            n_dft++;
        }

        /* Time series */
        int do_rec = (n % rec_every == 0);
        int do_print = (n % print_every2 == 0);

        if (do_rec || do_print) {
            double Et, Ec, fc;
            compute_energy(&Et, &Ec, &fc, phi, vel, Nx, dx, m2, core_r, -xmax);

            if (do_print) {
                /* Update Poisson for comparison */
                compute_rho(rho, phi, vel, Nx, dx, m2);
                double beta_avg = 0.002;
                for (int i = 0; i < Nx; i++)
                    rho_avg[i] = (1.0 - beta_avg) * rho_avg[i] + beta_avg * rho[i];
                solve_poisson(Phi_poi, rho_avg, Nx, dx, alpha_eff);
            }

            if (do_rec)
                fprintf(fts, "%.2f\t%.6e\t%.6e\t%.4f\t%.6e\t%.6e\t%.6e\n",
                        t, phi[0][ic], Et, fc, Phi[ic], Phi_poi[ic], alpha_eff);

            if (do_print)
                printf("  t=%7.1f  phi0=%+.4f  E=%.4f  fc=%.3f  "
                       "Phi_w(0)=%.4e  Phi_p(0)=%.4e  a_eff=%.4e\n",
                       t, phi[0][ic], Et, fc, Phi[ic], Phi_poi[ic], alpha_eff);
        }

        /* Probe: update running average */
        {
            double beta_ema = dt / ema_tau;
            if (beta_ema > 1.0) beta_ema = 1.0;
            Phi_probe_avg = (1.0 - beta_ema) * Phi_probe_avg + beta_ema * Phi[i_probe];
        }

        if (n % probe_every == 0) {
            fprintf(fprobe, "%.4f\t%.6e\t%.6e\t%.6e\t%.6e\n",
                    t, Phi[i_probe], Phi[ic], Phi_probe_avg,
                    rho[ic > 0 ? ic : 0]);
        }

        /* Calibrate pre-boost oscillation amplitude at probe */
        if (calibrating && T_elapsed < t_boost) {
            double dev = fabs(Phi_probe_avg - Phi_probe_pre_boost);
            if (dev > Phi_probe_max_osc) Phi_probe_max_osc = dev;
            /* Update baseline (slow drift) */
            double beta_base = dt / (5.0 * ema_tau);
            Phi_probe_pre_boost = (1.0 - beta_base) * Phi_probe_pre_boost
                                + beta_base * Phi_probe_avg;
        }

        /* Causality measurement: detect when EMA at probe deviates */
        if (measuring_delay && delay_arrival < 0) {
            double delta = fabs(Phi_probe_avg - Phi_probe_pre_boost);
            /* Threshold: 3x the pre-boost oscillation amplitude */
            double threshold = 3.0 * Phi_probe_max_osc + 1e-12;
            if (delta > threshold) {
                delay_arrival = t;
                printf("  >>> SIGNAL ARRIVED at probe x=%.1f: t=%.2f "
                       "(delay=%.2f, expected=%.2f) <<<\n",
                       x_probe, t, t - (t_equil + t_boost),
                       x_probe / c_light);
                printf("  EMA delta=%.4e  threshold=%.4e  max_osc=%.4e\n",
                       delta, threshold, Phi_probe_max_osc);
            }
        }

        /* Start calibrating 500 time units before boost */
        if (!calibrating && T_elapsed >= t_boost - 500.0 && T_elapsed < t_boost) {
            calibrating = 1;
            t_calib_start = t; (void)t_calib_start;
            Phi_probe_pre_boost = Phi_probe_avg;
            Phi_probe_max_osc = 0.0;
        }

        /* Phase 2: Perturbation for causality test.
         * DOUBLE alpha instantaneously. The change in source at x=0 should
         * reach x_probe at t_boost + x_probe/c. */
        if (!boosted && T_elapsed >= t_boost) {
            printf("\n  >>> PERTURBATION at t=%.1f: doubling alpha <<<\n", t);
            alpha *= 2.0;  /* double the coupling */
            printf("  New alpha = %.6e\n", alpha);
            boosted = 1;
            calibrating = 0;
            measuring_delay = 1;
            printf("  Phi_probe EMA baseline = %.6e\n", Phi_probe_pre_boost);
            printf("  Pre-boost oscillation amplitude = %.6e\n", Phi_probe_max_osc);
            printf("  Detection threshold (3x) = %.6e\n", 3.0 * Phi_probe_max_osc);
            printf("  Expected signal at x=%.0f: t=%.1f (delay=%.1f)\n",
                   x_probe, t + x_probe / c_light, x_probe / c_light);
        }

        if (n == Nt_grav) break;

        /* === Co-evolve fields and Phi using Velocity Verlet === */

        /* Half-step velocity update for fields */
        for (int a = 0; a < 3; a++)
            for (int i = 1; i < Nx - 1; i++)
                vel[a][i] += 0.5 * dt * acc[a][i];

        /* Half-step velocity update for Phi */
        for (int i = 1; i < Nx - 1; i++)
            Phi_vel[i] += 0.5 * dt * Phi_acc[i];

        /* Full-step position update */
        for (int a = 0; a < 3; a++)
            for (int i = 1; i < Nx - 1; i++)
                phi[a][i] += dt * vel[a][i];
        for (int i = 1; i < Nx - 1; i++)
            Phi[i] += dt * Phi_vel[i];

        /* Compute dPhi/dx for gravity coupling */
        for (int i = 1; i < Nx - 1; i++)
            dPhi[i] = (Phi[i+1] - Phi[i-1]) / (2.0 * dx);
        dPhi[0] = dPhi[Nx-1] = 0.0;

        /* Recompute accelerations */
        compute_rho(rho, phi, vel, Nx, dx, m2);
        COMPUTE_PHI_ACC(alpha_eff);

        if (fabs(alpha_eff) > 1e-15) {
            COMPUTE_ACC_GRAV();
        } else {
            COMPUTE_ACC_FLAT();
        }

        /* Second half-step velocity update */
        for (int a = 0; a < 3; a++)
            for (int i = 1; i < Nx - 1; i++)
                vel[a][i] += 0.5 * dt * acc[a][i];
        for (int i = 1; i < Nx - 1; i++)
            Phi_vel[i] += 0.5 * dt * Phi_acc[i];

        /* Absorbing boundaries */
        for (int a = 0; a < 3; a++)
            for (int i = 0; i < Nx; i++) {
                vel[a][i] *= damp_arr[i];
                phi[a][i] *= damp_arr[i];
            }
        for (int i = 0; i < Nx; i++) {
            Phi_vel[i] *= damp_phi[i];
            Phi[i] *= damp_phi[i];
        }
    }

    fclose(fts);
    fclose(fprobe);

    /* === Analysis === */
    printf("\n=== Analysis ===\n");

    /* Final Phi profile comparison */
    char phipath[600];
    snprintf(phipath, sizeof(phipath), "%s/retarded_phi_profile.tsv", outdir);
    FILE *fphi = fopen(phipath, "w");
    fprintf(fphi, "x\tPhi_wave\tPhi_poisson\trho_avg\n");
    /* Recompute final Poisson */
    compute_rho(rho, phi, vel, Nx, dx, m2);
    solve_poisson(Phi_poi, rho_avg, Nx, dx, alpha);
    for (int i = 0; i < Nx; i++) {
        double x = -xmax + i * dx;
        fprintf(fphi, "%.6f\t%.6e\t%.6e\t%.6e\n", x, Phi[i], Phi_poi[i], rho_avg[i]);
    }
    fclose(fphi);

    printf("  Phi_wave(0) = %.6e\n", Phi[ic]);
    printf("  Phi_poi(0)  = %.6e\n", Phi_poi[ic]);
    if (fabs(Phi_poi[ic]) > 1e-30)
        printf("  Ratio Phi_wave/Phi_poi = %.6f\n", Phi[ic] / Phi_poi[ic]);

    printf("\n  Causality test: probe at x=%.1f\n", x_probe);
    printf("  Boost time: t=%.1f\n", t_equil + t_boost);
    if (delay_arrival > 0) {
        double measured_delay = delay_arrival - (t_equil + t_boost);
        printf("  Signal arrival: t=%.2f\n", delay_arrival);
        printf("  Measured delay: %.2f  (expected: %.2f)\n",
               measured_delay, x_probe / c_light);
        printf("  Delay error: %.2f%%\n",
               100.0 * fabs(measured_delay - x_probe / c_light) / (x_probe / c_light));
    } else {
        printf("  Signal NOT detected at probe (threshold too high or time too short)\n");
    }

    /* DFT */
    int dft_start = n_dft / 2;
    double omega_peak = 0, power_peak = 0;
    if (n_dft - dft_start > 100) {
        char dftpath[600];
        snprintf(dftpath, sizeof(dftpath), "%s/retarded_spectrum.tsv", outdir);
        FILE *fdft = fopen(dftpath, "w");
        fprintf(fdft, "omega\tpower\n");
        double T = t_hist[n_dft-1] - t_hist[dft_start];
        int nf = 500;
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
            fprintf(fdft, "%.6f\t%.6e\n", omega, pw);
            if (pw > power_peak) { power_peak = pw; omega_peak = omega; }
        }
        fclose(fdft);
        printf("  omega_peak = %.4f (mass gap = %.4f)\n", omega_peak, mass);
    }

    printf("  Oscillon survived? %s\n",
           (omega_peak > 0.01 && omega_peak < mass) ? "YES" : "NO");

    double E_final, Ec_final, fc_final;
    compute_energy(&E_final, &Ec_final, &fc_final, phi, vel, Nx, dx, m2, core_r, -xmax);
    printf("  E_equil=%.4f  E_final=%.4f  dE/E=%.4e\n",
           E_equil, E_final, (E_final - E_equil) / fabs(E_equil));

    printf("\nOutput:\n");
    printf("  %s/retarded_ts.tsv\n", outdir);
    printf("  %s/retarded_probe.tsv\n", outdir);
    printf("  %s/retarded_phi_profile.tsv\n", outdir);
    printf("  %s/retarded_spectrum.tsv\n", outdir);

    for (int a = 0; a < 3; a++) { free(phi[a]); free(vel[a]); free(acc[a]); }
    free(Phi); free(Phi_vel); free(Phi_acc); free(Phi_poi);
    free(rho); free(rho_avg); free(dPhi);
    free(damp_arr); free(damp_phi);
    free(phi0_hist); free(t_hist);

    return 0;
}
