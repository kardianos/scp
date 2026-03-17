/*
 * selfref.c — Self-consistent metric from oscillon stress-energy
 *
 * Test F: Solve Poisson for gravitational potential Phi from the oscillon's
 * TIME-AVERAGED energy density rho, then evolve the fields on the curved metric.
 *
 * Physics:
 *   Poisson:    d²Phi/dx² = alpha * <rho(x)>_t
 *   KG on curved metric (weak field):
 *     d²phi/dt² = (1+4Phi) d²phi/dx² + 2(dPhi/dx)(dphi/dx) - m²(1+2Phi)phi + F_pot
 *
 * Key improvements over naive approach:
 *   - Time-averaged rho prevents oscillation-driven pumping instability
 *   - Adiabatic ramp-up of alpha over t_ramp period
 *   - Scans multiple alpha values in one run
 *
 * Compile: gcc -O3 -Wall -o selfref src/selfref.c -lm
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Parameters */
static double mu      = -20.0;
static double kappa   = 20.0;
static double mass    = 1.0;
static double A_init  = 0.8;
static double sigma_init = 3.0;
static int    Nx      = 4000;
static double xmax    = 100.0;
static double t_equil = 5000.0;
static double t_grav  = 5000.0;   /* gravity-on evolution time per alpha */
static double t_ramp  = 500.0;    /* adiabatic ramp-up time */
static char   outdir[512] = "v24/fundamental/testF_selfref/data";

/* Alpha scan values */
static double alpha_scan[] = {0.0001, 0.001, 0.01, 0.1};
static int n_alpha = 4;

static void parse_args(int argc, char **argv)
{
    for (int i = 1; i < argc - 1; i += 2) {
        if      (!strcmp(argv[i], "-mu"))     mu      = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-kappa"))  kappa   = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-mass"))   mass    = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-A"))      A_init  = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-sigma"))  sigma_init = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-Nx"))     Nx      = atoi(argv[i+1]);
        else if (!strcmp(argv[i], "-xmax"))   xmax    = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-equil"))  t_equil = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-tgrav"))  t_grav  = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-tramp"))  t_ramp  = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-o"))      strncpy(outdir, argv[i+1], sizeof(outdir)-1);
        else { fprintf(stderr, "Unknown flag: %s\n", argv[i]); exit(1); }
    }
}

/* -dV/dphi_a where V = (mu/2)P^2/(1+kappa P^2), P = phi1*phi2*phi3 */
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

/*
 * Solve Poisson: d²Phi/dx² = alp_eff * rho(x)
 * Direct double integration with Phi=0 at boundaries.
 */
static void solve_poisson(double *Phi, const double *rho, int N, double dx, double alp_eff)
{
    /* First integration: F(x) = integral_{0}^{x} rho dx' */
    double *F = malloc(N * sizeof(double));
    F[0] = 0.0;
    for (int i = 1; i < N; i++)
        F[i] = F[i-1] + alp_eff * rho[i] * dx;

    /* Second integration: G(x) = integral_{0}^{x} F dx' */
    Phi[0] = 0.0;
    for (int i = 1; i < N; i++)
        Phi[i] = Phi[i-1] + F[i] * dx;

    /* Subtract linear trend: Phi -> Phi - Phi(N-1)*x/L so Phi(N-1)=0 */
    double L = (N - 1) * dx;
    if (L > 0) {
        double slope = Phi[N-1] / L;
        for (int i = 0; i < N; i++)
            Phi[i] -= slope * (i * dx);
    }

    free(F);
}

/*
 * Compute energy density rho(x)
 */
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

/*
 * Compute total energy and core fraction
 */
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

int main(int argc, char **argv)
{
    parse_args(argc, argv);

    double dx  = 2.0 * xmax / (Nx - 1);
    double dx2 = dx * dx;
    double m2  = mass * mass;

    /* CFL */
    double kmax = M_PI / dx;
    double dt = 0.4 * 2.0 / sqrt(kmax * kmax + m2);

    int Nt_equil = (int)(t_equil / dt) + 1;

    printf("selfref: Self-consistent metric test\n");
    printf("  mu=%.1f kappa=%.1f mass=%.4f A=%.3f sigma=%.3f\n",
           mu, kappa, mass, A_init, sigma_init);
    printf("  Nx=%d xmax=%.1f dx=%.5f dt=%.6f\n", Nx, xmax, dx, dt);
    printf("  t_equil=%.0f t_grav=%.0f t_ramp=%.0f\n", t_equil, t_grav, t_ramp);
    printf("  Alpha scan:");
    for (int ia = 0; ia < n_alpha; ia++) printf(" %.4f", alpha_scan[ia]);
    printf("\n");

    /* Allocate fields */
    double *phi[3], *vel[3], *acc[3];
    double *phi_save[3], *vel_save[3];  /* save state for resetting between alphas */
    for (int a = 0; a < 3; a++) {
        phi[a] = calloc(Nx, sizeof(double));
        vel[a] = calloc(Nx, sizeof(double));
        acc[a] = calloc(Nx, sizeof(double));
        phi_save[a] = calloc(Nx, sizeof(double));
        vel_save[a] = calloc(Nx, sizeof(double));
    }
    double *Phi       = calloc(Nx, sizeof(double));
    double *rho       = calloc(Nx, sizeof(double));
    double *rho_avg   = calloc(Nx, sizeof(double));  /* time-averaged rho */
    double *dPhi      = calloc(Nx, sizeof(double));

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
                double dphi = (phi[a][i+1] - phi[a][i-1]) / (2.0 * dx); \
                acc[a][i] = (1.0 + 4.0*Phi[i]) * lapl \
                          + 2.0 * dPhi[i] * dphi \
                          - m2 * (1.0 + 2.0*Phi[i]) * phi[a][i] \
                          + fp; \
            } \
        } \
    } while(0)

    /* === Phase 1: Equilibrate on flat metric === */
    printf("\n=== Phase 1: Equilibrate (flat metric, t=0..%.0f) ===\n", t_equil);

    COMPUTE_ACC_FLAT();

    int print_every = Nt_equil / 20;
    if (print_every < 1) print_every = 1;

    for (int n = 0; n < Nt_equil; n++) {
        double t = n * dt;

        if (n % print_every == 0) {
            double Et, Ec, fc;
            compute_energy(&Et, &Ec, &fc, phi, vel, Nx, dx, m2, core_r, -xmax);
            printf("  t=%7.1f  phi0=(%+.4f,%+.4f,%+.4f)  E=%.4f  fc=%.3f\n",
                   t, phi[0][ic], phi[1][ic], phi[2][ic], Et, fc);
        }

        /* Velocity Verlet */
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

        /* Absorbing boundary */
        for (int a = 0; a < 3; a++)
            for (int i = 0; i < Nx; i++) {
                vel[a][i] *= damp[i];
                phi[a][i] *= damp[i];
            }
    }

    /* Record equilibrated properties */
    double E_equil, Ec_equil, fc_equil;
    compute_energy(&E_equil, &Ec_equil, &fc_equil, phi, vel, Nx, dx, m2, core_r, -xmax);
    printf("\n  Equilibrated: E=%.4f  fc=%.4f  phi0=%.4f\n", E_equil, fc_equil, phi[0][ic]);

    /* Save equilibrated state */
    for (int a = 0; a < 3; a++) {
        memcpy(phi_save[a], phi[a], Nx * sizeof(double));
        memcpy(vel_save[a], vel[a], Nx * sizeof(double));
    }

    /* Compute initial time-averaged rho (from the equilibrated oscillon).
     * Average over one full oscillation period (~2pi/omega ~ 8). */
    printf("\n  Computing time-averaged rho over one period...\n");
    double T_osc = 8.0;  /* approximate period */
    int N_avg = (int)(T_osc / dt);
    for (int i = 0; i < Nx; i++) rho_avg[i] = 0.0;

    /* Restore state and evolve flat for one period to compute average */
    for (int a = 0; a < 3; a++) {
        memcpy(phi[a], phi_save[a], Nx * sizeof(double));
        memcpy(vel[a], vel_save[a], Nx * sizeof(double));
    }
    COMPUTE_ACC_FLAT();
    for (int n = 0; n < N_avg; n++) {
        compute_rho(rho, phi, vel, Nx, dx, m2);
        for (int i = 0; i < Nx; i++) rho_avg[i] += rho[i] / N_avg;

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
                vel[a][i] *= damp[i];
                phi[a][i] *= damp[i];
            }
    }
    double rho_peak_avg = 0;
    for (int i = 0; i < Nx; i++)
        if (rho_avg[i] > rho_peak_avg) rho_peak_avg = rho_avg[i];
    printf("  <rho> peak = %.6f\n", rho_peak_avg);

    /* === Summary file === */
    char sumpath[600];
    snprintf(sumpath, sizeof(sumpath), "%s/selfref_summary.tsv", outdir);
    FILE *fsum = fopen(sumpath, "w");
    fprintf(fsum, "alpha\tPhi_0\tPhi_min\tE_final\tE_equil\tdE_pct\tfc\t"
                  "omega\tsurvived\tcollapse_t\n");

    /* === Phase 2: Scan alpha values === */
    for (int ia = 0; ia < n_alpha; ia++) {
        double alpha = alpha_scan[ia];
        int Nt_grav = (int)(t_grav / dt) + 1;
        int Nt_ramp = (int)(t_ramp / dt) + 1;

        printf("\n========================================================\n");
        printf("=== Alpha = %.6f (scan %d/%d) ===\n", alpha, ia+1, n_alpha);
        printf("========================================================\n");

        /* Restore equilibrated state */
        for (int a = 0; a < 3; a++) {
            memcpy(phi[a], phi_save[a], Nx * sizeof(double));
            memcpy(vel[a], vel_save[a], Nx * sizeof(double));
        }
        memset(Phi, 0, Nx * sizeof(double));
        memset(dPhi, 0, Nx * sizeof(double));

        /* Solve Poisson with time-averaged rho at full alpha to show initial Phi */
        solve_poisson(Phi, rho_avg, Nx, dx, alpha);
        printf("  Phi(0) from <rho>: %.6e  (c_eff^2 = %.4f, m_eff^2/m^2 = %.4f)\n",
               Phi[ic], 1.0 + 4.0*Phi[ic], 1.0 + 2.0*Phi[ic]);

        /* Check weak-field validity */
        double Phi0_full = Phi[ic];
        int weak_field = (fabs(Phi0_full) < 0.5);
        if (!weak_field)
            printf("  WARNING: |Phi(0)| = %.3f > 0.5 — weak-field approximation breaks down!\n",
                   fabs(Phi0_full));

        /* Reset Phi for adiabatic ramp */
        memset(Phi, 0, Nx * sizeof(double));
        memset(dPhi, 0, Nx * sizeof(double));

        COMPUTE_ACC_FLAT();

        /* Time series output */
        char tspath[600];
        snprintf(tspath, sizeof(tspath), "%s/selfref_alpha%.4f_ts.tsv", outdir, alpha);
        FILE *fts = fopen(tspath, "w");
        fprintf(fts, "time\tphi1_0\tE_total\tfc\tPhi_0\trho_peak\talpha_eff\n");

        /* DFT storage */
        int max_dft = 50000;
        double *phi0_hist = malloc(max_dft * sizeof(double));
        double *t_hist    = malloc(max_dft * sizeof(double));
        int n_dft = 0;
        int dft_every = Nt_grav / max_dft;
        if (dft_every < 1) dft_every = 1;

        int print_every2 = Nt_grav / 40;
        if (print_every2 < 1) print_every2 = 1;
        int rec_every = Nt_grav / 10000;
        if (rec_every < 1) rec_every = 1;

        /* Running average of rho (exponential moving average) */
        double *rho_run = calloc(Nx, sizeof(double));
        memcpy(rho_run, rho_avg, Nx * sizeof(double));
        double tau_avg = 20.0;  /* averaging timescale in code units (~2.5 periods) */

        int collapsed = 0;
        double collapse_t = -1.0;

        for (int n = 0; n <= Nt_grav; n++) {
            double t = n * dt;
            double t_abs = t_equil + t;

            /* Adiabatic ramp: alpha_eff = alpha * ramp(t) */
            double ramp = (n < Nt_ramp) ? (double)n / Nt_ramp : 1.0;
            ramp = ramp * ramp * (3.0 - 2.0*ramp);  /* smooth Hermite step */
            double alpha_eff = alpha * ramp;

            /* Update running average of rho every 10 steps */
            if (n % 10 == 0 && !collapsed) {
                compute_rho(rho, phi, vel, Nx, dx, m2);
                double beta = dt * 10.0 / tau_avg;
                if (beta > 1.0) beta = 1.0;
                for (int i = 0; i < Nx; i++)
                    rho_run[i] = (1.0 - beta) * rho_run[i] + beta * rho[i];

                /* Solve Poisson with running-average rho */
                solve_poisson(Phi, rho_run, Nx, dx, alpha_eff);
                for (int i = 1; i < Nx - 1; i++)
                    dPhi[i] = (Phi[i+1] - Phi[i-1]) / (2.0 * dx);
                dPhi[0] = dPhi[Nx-1] = 0.0;
            }

            /* Record DFT samples */
            if (n % dft_every == 0 && n_dft < max_dft) {
                phi0_hist[n_dft] = phi[0][ic];
                t_hist[n_dft] = t_abs;
                n_dft++;
            }

            int do_print = (n % print_every2 == 0);
            int do_rec = (n % rec_every == 0);

            if ((do_rec || do_print) && !collapsed) {
                double Et, Ec, fc;
                compute_energy(&Et, &Ec, &fc, phi, vel, Nx, dx, m2, core_r, -xmax);
                double rp = 0;
                for (int i = 0; i < Nx; i++) if (rho_run[i] > rp) rp = rho_run[i];

                /* Check for NaN or collapse */
                if (isnan(Et) || isnan(phi[0][ic]) || Et > 100.0 * E_equil) {
                    printf("  t=%7.1f  *** COLLAPSED (E=%.4e) ***\n", t, Et);
                    collapsed = 1;
                    collapse_t = t;
                }

                if (do_rec && !collapsed)
                    fprintf(fts, "%.2f\t%.6e\t%.6e\t%.4f\t%.6e\t%.6e\t%.6e\n",
                            t_abs, phi[0][ic], Et, fc, Phi[ic], rp, alpha_eff);

                if (do_print && !collapsed)
                    printf("  t=%7.1f  phi0=%+.4f  E=%.4f  fc=%.3f  Phi(0)=%.4e  a_eff=%.4e\n",
                           t, phi[0][ic], Et, fc, Phi[ic], alpha_eff);
            }

            if (n == Nt_grav || collapsed) break;

            /* Velocity Verlet with gravity */
            for (int a = 0; a < 3; a++)
                for (int i = 1; i < Nx - 1; i++)
                    vel[a][i] += 0.5 * dt * acc[a][i];
            for (int a = 0; a < 3; a++)
                for (int i = 1; i < Nx - 1; i++)
                    phi[a][i] += dt * vel[a][i];

            /* Use gravity-corrected acceleration when Phi is nonzero */
            if (alpha_eff > 1e-15) {
                COMPUTE_ACC_GRAV();
            } else {
                COMPUTE_ACC_FLAT();
            }

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
        fclose(fts);

        /* Final measurements */
        double E_final, Ec_final, fc_final;
        double omega_peak = 0, power_peak = 0;

        if (!collapsed) {
            compute_energy(&E_final, &Ec_final, &fc_final, phi, vel, Nx, dx, m2, core_r, -xmax);

            printf("\n  --- Results for alpha=%.6f ---\n", alpha);
            printf("  E_equil=%.4f  E_final=%.4f  dE/E=%.4e\n",
                   E_equil, E_final, (E_final - E_equil) / fabs(E_equil));
            printf("  fc=%.4f  Phi(0)=%.6e\n", fc_final, Phi[ic]);
        } else {
            E_final = 0; fc_final = 0;
            printf("\n  --- Results for alpha=%.6f: COLLAPSED at t=%.1f ---\n",
                   alpha, collapse_t);
        }

        /* Write Phi profile */
        char phipath[600];
        snprintf(phipath, sizeof(phipath), "%s/selfref_alpha%.4f_phi.tsv", outdir, alpha);
        FILE *fphi = fopen(phipath, "w");
        fprintf(fphi, "x\tPhi\trho_avg\n");
        for (int i = 0; i < Nx; i++) {
            double x = -xmax + i * dx;
            fprintf(fphi, "%.6f\t%.6e\t%.6e\n", x, Phi[i], rho_run[i]);
        }
        fclose(fphi);

        /* DFT of phi_1(0,t) — second half */
        int dft_start = n_dft / 2;
        if (!collapsed && n_dft - dft_start > 100) {
            char dftpath[600];
            snprintf(dftpath, sizeof(dftpath), "%s/selfref_alpha%.4f_spectrum.tsv", outdir, alpha);
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

        int survived = (!collapsed && omega_peak > 0.01 && omega_peak < mass);
        printf("  Oscillon survived? %s\n", survived ? "YES" : "NO");

        /* Summary line */
        fprintf(fsum, "%.6f\t%.6e\t%.6e\t%.6e\t%.6e\t%.4f\t%.4f\t%.4f\t%d\t%.1f\n",
                alpha, Phi[ic], Phi0_full, E_final, E_equil,
                collapsed ? 999.0 : 100.0*(E_final-E_equil)/fabs(E_equil),
                fc_final, omega_peak, survived, collapse_t);

        free(phi0_hist);
        free(t_hist);
        free(rho_run);
    }

    fclose(fsum);
    printf("\nSummary: %s\n", sumpath);

    /* Cleanup */
    for (int a = 0; a < 3; a++) {
        free(phi[a]); free(vel[a]); free(acc[a]);
        free(phi_save[a]); free(vel_save[a]);
    }
    free(Phi); free(rho); free(rho_avg); free(dPhi); free(damp);
    return 0;
}
