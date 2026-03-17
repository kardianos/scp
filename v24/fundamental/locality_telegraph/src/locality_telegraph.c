/*
 * locality_telegraph.c — Telegraph equation for gravitational potential
 *
 * Option 3: Interpolate between wave (causal) and Poisson (instantaneous)
 * using the telegraph equation:
 *
 *   d²Phi/dt² + gamma * dPhi/dt = c² d²Phi/dx² + alpha * rho
 *
 * gamma controls damping:
 *   gamma=0:    pure wave (oscillating Phi, radiation problem)
 *   gamma->inf: overdamped -> Poisson (instantaneous)
 *   gamma~c/L:  causal propagation with rapid relaxation
 *
 * Phase 1: gamma scan with single oscillon
 * Phase 2: causality test at gamma* (boost oscillon, measure delay at x=50)
 *
 * Compile: gcc -O3 -Wall -o locality_telegraph src/locality_telegraph.c -lm
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
static int    Nx      = 8000;
static double xmax    = 200.0;
static double t_equil = 5000.0;
static double t_test  = 5000.0;
static double alpha   = -1e-4;
static char   outdir[512] = "v24/fundamental/locality_telegraph/data";

/* Gamma scan values */
static double gamma_scan[] = {0.01, 0.05, 0.1, 0.2, 0.5, 1.0, 5.0, 20.0};
static int n_gamma = 8;

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
        else if (!strcmp(argv[i], "-ttest"))  t_test  = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-alpha"))  alpha   = atof(argv[i+1]);
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

/* Compute energy density rho(x) */
static void compute_rho(double *rho, double **phi, double **vel,
                        int N, double dx, double m2)
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

/* Solve Poisson: d²Phi/dx² = alp * rho(x), for reference comparison */
static void solve_poisson(double *Phi, const double *rho, int N,
                         double dx, double alp)
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

    /* CFL: need dt < dx/c for wave equation, also dt < 2/gamma for stability
     * with large gamma. Use conservative factor. */
    double kmax = M_PI / dx;
    double dt = 0.4 * 2.0 / sqrt(kmax * kmax + m2);

    int Nt_equil = (int)(t_equil / dt) + 1;

    printf("locality_telegraph: Telegraph equation for gravity\n");
    printf("  mu=%.1f kappa=%.1f mass=%.4f A=%.3f sigma=%.3f\n",
           mu, kappa, mass, A_init, sigma_init);
    printf("  Nx=%d xmax=%.1f dx=%.5f dt=%.6f\n", Nx, xmax, dx, dt);
    printf("  alpha=%.6e  t_equil=%.0f  t_test=%.0f\n", alpha, t_equil, t_test);
    printf("  Gamma scan:");
    for (int ig = 0; ig < n_gamma; ig++) printf(" %.2f", gamma_scan[ig]);
    printf("\n");

    /* Allocate fields */
    double *phi[3], *vel[3], *acc[3];
    double *phi_save[3], *vel_save[3];
    for (int a = 0; a < 3; a++) {
        phi[a] = calloc(Nx, sizeof(double));
        vel[a] = calloc(Nx, sizeof(double));
        acc[a] = calloc(Nx, sizeof(double));
        phi_save[a] = calloc(Nx, sizeof(double));
        vel_save[a] = calloc(Nx, sizeof(double));
    }
    double *Phi       = calloc(Nx, sizeof(double));  /* telegraph potential */
    double *Phi_vel   = calloc(Nx, sizeof(double));  /* dPhi/dt */
    double *Phi_acc   = calloc(Nx, sizeof(double));  /* d²Phi/dt² */
    double *rho       = calloc(Nx, sizeof(double));
    double *rho_avg   = calloc(Nx, sizeof(double));
    double *Phi_poisson = calloc(Nx, sizeof(double)); /* reference Poisson solution */

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

    int ic = Nx / 2;  /* center index */
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

    /* ============================================================ */
    /* Phase 0: Equilibrate on flat metric                          */
    /* ============================================================ */
    printf("\n=== Phase 0: Equilibrate (flat metric, t=0..%.0f) ===\n", t_equil);

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
    printf("\n  Equilibrated: E=%.4f  fc=%.4f  phi0=%.4f\n",
           E_equil, fc_equil, phi[0][ic]);

    /* Save equilibrated state */
    for (int a = 0; a < 3; a++) {
        memcpy(phi_save[a], phi[a], Nx * sizeof(double));
        memcpy(vel_save[a], vel[a], Nx * sizeof(double));
    }

    /* Compute time-averaged rho over one oscillation period */
    printf("\n  Computing time-averaged rho...\n");
    double T_osc = 8.0;
    int N_avg = (int)(T_osc / dt);
    for (int i = 0; i < Nx; i++) rho_avg[i] = 0.0;

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

    /* Solve Poisson for reference */
    solve_poisson(Phi_poisson, rho_avg, Nx, dx, alpha);
    printf("  Poisson Phi(0) = %.6e\n", Phi_poisson[ic]);

    /* ============================================================ */
    /* Phase 1: Gamma scan                                          */
    /* ============================================================ */

    /* Summary file */
    char sumpath[600];
    snprintf(sumpath, sizeof(sumpath), "%s/telegraph_summary.tsv", outdir);
    FILE *fsum = fopen(sumpath, "w");
    fprintf(fsum, "gamma\tPhi_static\tPhi_osc_amp\tPhi_osc_ratio\tfc\t"
                  "E_final\tomega\tsurvived\tPhi_poisson\n");

    /* Best gamma tracking */
    double best_gamma = -1.0;
    double best_ratio = 1e30;

    for (int ig = 0; ig < n_gamma; ig++) {
        double gam = gamma_scan[ig];
        int Nt_test = (int)(t_test / dt) + 1;

        printf("\n========================================================\n");
        printf("=== Gamma = %.4f (scan %d/%d) ===\n", gam, ig+1, n_gamma);
        printf("========================================================\n");

        /* Restore equilibrated state */
        for (int a = 0; a < 3; a++) {
            memcpy(phi[a], phi_save[a], Nx * sizeof(double));
            memcpy(vel[a], vel_save[a], Nx * sizeof(double));
        }
        memset(Phi, 0, Nx * sizeof(double));
        memset(Phi_vel, 0, Nx * sizeof(double));
        memset(Phi_acc, 0, Nx * sizeof(double));

        COMPUTE_ACC_FLAT();

        /* Time series output */
        char tspath[600];
        snprintf(tspath, sizeof(tspath), "%s/telegraph_gamma%.4f_ts.tsv", outdir, gam);
        FILE *fts = fopen(tspath, "w");
        fprintf(fts, "time\tphi1_0\tE_total\tfc\tPhi_0\tPhi_vel_0\trho_0\n");

        /* DFT storage for phi_1(0,t) */
        int max_dft = 50000;
        double *phi0_hist = malloc(max_dft * sizeof(double));
        double *t_hist    = malloc(max_dft * sizeof(double));
        int n_dft = 0;
        int dft_every = Nt_test / max_dft;
        if (dft_every < 1) dft_every = 1;

        /* Track Phi(0) for oscillation measurement (last 20% of run) */
        int phi_track_start = (int)(0.8 * Nt_test);
        double Phi0_max = -1e30, Phi0_min = 1e30;
        double Phi0_sum = 0.0;
        int Phi0_count = 0;

        int print_every2 = Nt_test / 40;
        if (print_every2 < 1) print_every2 = 1;
        int rec_every = Nt_test / 10000;
        if (rec_every < 1) rec_every = 1;

        /* Running average of rho for source term */
        double *rho_run = calloc(Nx, sizeof(double));
        memcpy(rho_run, rho_avg, Nx * sizeof(double));
        double tau_avg = 20.0;

        /* Adiabatic ramp for alpha */
        double t_ramp = 500.0;
        int Nt_ramp = (int)(t_ramp / dt) + 1;

        int collapsed = 0;

        for (int n = 0; n <= Nt_test; n++) {
            double t = n * dt;

            /* Adiabatic ramp */
            double ramp = (n < Nt_ramp) ? (double)n / Nt_ramp : 1.0;
            ramp = ramp * ramp * (3.0 - 2.0 * ramp);  /* Hermite smoothstep */
            double alpha_eff = alpha * ramp;

            /* Update running average of rho every 10 steps */
            if (n % 10 == 0 && !collapsed) {
                compute_rho(rho, phi, vel, Nx, dx, m2);
                double beta = dt * 10.0 / tau_avg;
                if (beta > 1.0) beta = 1.0;
                for (int i = 0; i < Nx; i++)
                    rho_run[i] = (1.0 - beta) * rho_run[i] + beta * rho[i];
            }

            /* Record DFT samples */
            if (n % dft_every == 0 && n_dft < max_dft) {
                phi0_hist[n_dft] = phi[0][ic];
                t_hist[n_dft] = t_equil + t;
                n_dft++;
            }

            /* Track Phi(0) oscillations in last 20% */
            if (n >= phi_track_start && !collapsed) {
                if (Phi[ic] > Phi0_max) Phi0_max = Phi[ic];
                if (Phi[ic] < Phi0_min) Phi0_min = Phi[ic];
                Phi0_sum += Phi[ic];
                Phi0_count++;
            }

            int do_print = (n % print_every2 == 0);
            int do_rec = (n % rec_every == 0);

            if ((do_rec || do_print) && !collapsed) {
                double Et, Ec, fc;
                compute_energy(&Et, &Ec, &fc, phi, vel, Nx, dx, m2, core_r, -xmax);

                if (isnan(Et) || isnan(phi[0][ic]) || isnan(Phi[ic]) ||
                    Et > 100.0 * fabs(E_equil)) {
                    printf("  t=%7.1f  *** COLLAPSED ***\n", t);
                    collapsed = 1;
                }

                if (do_rec && !collapsed)
                    fprintf(fts, "%.2f\t%.6e\t%.6e\t%.4f\t%.6e\t%.6e\t%.6e\n",
                            t, phi[0][ic], Et, fc, Phi[ic], Phi_vel[ic], rho_run[ic]);

                if (do_print && !collapsed)
                    printf("  t=%7.1f  phi0=%+.4f  E=%.4f  fc=%.3f  Phi(0)=%.4e  vPhi=%.4e\n",
                           t, phi[0][ic], Et, fc, Phi[ic], Phi_vel[ic]);
            }

            if (n == Nt_test || collapsed) break;

            /* ===== Velocity Verlet for scalar fields (no damping on phi_a) ===== */
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

            /* Absorbing boundary on scalar fields */
            for (int a = 0; a < 3; a++)
                for (int i = 0; i < Nx; i++) {
                    vel[a][i] *= damp[i];
                    phi[a][i] *= damp[i];
                }

            /* ===== Velocity Verlet for Phi (telegraph equation with damping) ===== */
            /*
             * d²Phi/dt² + gamma * dPhi/dt = c² d²Phi/dx² + alpha * rho
             *
             * Rewrite: d²Phi/dt² = c² d²Phi/dx² + alpha*rho - gamma*dPhi/dt
             *
             * Velocity Verlet with velocity-dependent damping:
             *   Phi_vel += 0.5*dt*Phi_acc                      (half-kick)
             *   Phi += dt*Phi_vel                               (drift)
             *   Phi_acc_new = c²*lapl(Phi) + alpha*rho - gamma*Phi_vel
             *   Phi_vel += 0.5*dt*Phi_acc_new                  (half-kick)
             *
             * For stability with large gamma, use implicit treatment of damping:
             *   After full step: Phi_vel_new = Phi_vel_explicit / (1 + 0.5*gamma*dt)
             *   This is the standard semi-implicit velocity Verlet for friction.
             */

            /* Half-kick Phi velocity */
            for (int i = 1; i < Nx - 1; i++)
                Phi_vel[i] += 0.5 * dt * Phi_acc[i];

            /* Drift Phi */
            for (int i = 1; i < Nx - 1; i++)
                Phi[i] += dt * Phi_vel[i];

            /* Compute new Phi acceleration: c²*lapl(Phi) + alpha*rho */
            for (int i = 1; i < Nx - 1; i++) {
                double lapl_Phi = (Phi[i+1] - 2.0*Phi[i] + Phi[i-1]) / dx2;
                /* c=1 in code units */
                Phi_acc[i] = lapl_Phi + alpha_eff * rho_run[i];
            }
            Phi_acc[0] = Phi_acc[Nx-1] = 0;

            /* Half-kick + semi-implicit damping */
            for (int i = 1; i < Nx - 1; i++) {
                Phi_vel[i] += 0.5 * dt * Phi_acc[i];
                /* Semi-implicit damping: divide by (1 + gamma*dt)
                 * This approximates the exponential decay e^{-gamma*dt} */
                Phi_vel[i] /= (1.0 + gam * dt);
            }

            /* Absorbing boundary on Phi */
            for (int i = 0; i < Nx; i++) {
                Phi[i]     *= damp[i];
                Phi_vel[i] *= damp[i];
            }

            /* Boundary conditions: Phi=0 at edges */
            Phi[0] = Phi[Nx-1] = 0;
            Phi_vel[0] = Phi_vel[Nx-1] = 0;
        }
        fclose(fts);

        /* Final measurements */
        double E_final, Ec_final, fc_final;
        double omega_peak = 0, power_peak = 0;
        double Phi_static = 0, Phi_osc_amp = 0, Phi_osc_ratio = 0;

        if (!collapsed) {
            compute_energy(&E_final, &Ec_final, &fc_final, phi, vel, Nx, dx, m2, core_r, -xmax);

            /* Phi oscillation analysis */
            if (Phi0_count > 0) {
                Phi_static = Phi0_sum / Phi0_count;
                Phi_osc_amp = 0.5 * (Phi0_max - Phi0_min);
                if (fabs(Phi_static) > 1e-20)
                    Phi_osc_ratio = Phi_osc_amp / fabs(Phi_static);
            }

            printf("\n  --- Results for gamma=%.4f ---\n", gam);
            printf("  E_equil=%.4f  E_final=%.4f  dE/E=%.4e\n",
                   E_equil, E_final, (E_final - E_equil) / fabs(E_equil));
            printf("  fc=%.4f\n", fc_final);
            printf("  Phi_static=%.6e  Phi_osc_amp=%.6e  ratio=%.4f\n",
                   Phi_static, Phi_osc_amp, Phi_osc_ratio);
            printf("  Phi_poisson=%.6e\n", Phi_poisson[ic]);
        } else {
            E_final = 0; fc_final = 0;
            printf("\n  --- gamma=%.4f: COLLAPSED ---\n", gam);
        }

        /* DFT of phi_1(0,t) */
        int dft_start = n_dft / 2;
        if (!collapsed && n_dft - dft_start > 100) {
            char dftpath[600];
            snprintf(dftpath, sizeof(dftpath),
                     "%s/telegraph_gamma%.4f_spectrum.tsv", outdir, gam);
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

        /* Track best gamma: causal (small ratio) with survival */
        if (survived && Phi_osc_ratio < 0.1 && Phi_osc_ratio < best_ratio) {
            best_ratio = Phi_osc_ratio;
            best_gamma = gam;
        }
        /* If none has ratio<0.1, track smallest ratio among survivors */
        if (survived && best_gamma < 0 && Phi_osc_ratio < best_ratio) {
            best_ratio = Phi_osc_ratio;
            best_gamma = gam;
        }

        /* Write Phi profile */
        char phipath[600];
        snprintf(phipath, sizeof(phipath),
                 "%s/telegraph_gamma%.4f_phi_profile.tsv", outdir, gam);
        FILE *fphi = fopen(phipath, "w");
        fprintf(fphi, "x\tPhi_telegraph\tPhi_poisson\trho_avg\n");
        for (int i = 0; i < Nx; i++) {
            double x = -xmax + i * dx;
            fprintf(fphi, "%.6f\t%.6e\t%.6e\t%.6e\n",
                    x, Phi[i], Phi_poisson[i], rho_run[i]);
        }
        fclose(fphi);

        /* Summary line */
        fprintf(fsum, "%.4f\t%.6e\t%.6e\t%.6f\t%.4f\t%.6e\t%.4f\t%d\t%.6e\n",
                gam, Phi_static, Phi_osc_amp, Phi_osc_ratio, fc_final,
                E_final, omega_peak, survived, Phi_poisson[ic]);

        free(phi0_hist);
        free(t_hist);
        free(rho_run);
    }
    fclose(fsum);

    printf("\n\n========================================================\n");
    printf("=== Phase 1 Summary ===\n");
    printf("========================================================\n");
    printf("  Poisson Phi(0) = %.6e\n", Phi_poisson[ic]);
    if (best_gamma > 0)
        printf("  Best gamma* = %.4f (osc ratio = %.4f)\n", best_gamma, best_ratio);
    else
        printf("  No gamma found with osc_ratio < 0.1 AND survival\n");

    /* ============================================================ */
    /* Phase 2: Causality test at gamma*                            */
    /* ============================================================ */

    /* Pick gamma*: use best from scan, or fall back to smallest-ratio survivor */
    double gamma_star = best_gamma;
    if (gamma_star < 0) {
        /* Fall back: pick gamma=0.2 (physically motivated) */
        gamma_star = 0.2;
        printf("  Using fallback gamma* = %.4f for causality test\n", gamma_star);
    }

    printf("\n========================================================\n");
    printf("=== Phase 2: Causality test at gamma* = %.4f ===\n", gamma_star);
    printf("========================================================\n");

    /* Restore equilibrated state */
    for (int a = 0; a < 3; a++) {
        memcpy(phi[a], phi_save[a], Nx * sizeof(double));
        memcpy(vel[a], vel_save[a], Nx * sizeof(double));
    }
    memset(Phi, 0, Nx * sizeof(double));
    memset(Phi_vel, 0, Nx * sizeof(double));
    memset(Phi_acc, 0, Nx * sizeof(double));

    COMPUTE_ACC_FLAT();

    /* First: let telegraph Phi settle with static oscillon for t_settle */
    double t_settle = 2000.0;
    int Nt_settle = (int)(t_settle / dt) + 1;
    double t_ramp2 = 500.0;
    int Nt_ramp2 = (int)(t_ramp2 / dt) + 1;

    printf("  Settling Phi for t=%.0f...\n", t_settle);

    double *rho_run2 = calloc(Nx, sizeof(double));
    memcpy(rho_run2, rho_avg, Nx * sizeof(double));
    double tau_avg2 = 20.0;

    for (int n = 0; n < Nt_settle; n++) {
        double ramp = (n < Nt_ramp2) ? (double)n / Nt_ramp2 : 1.0;
        ramp = ramp * ramp * (3.0 - 2.0 * ramp);
        double alpha_eff = alpha * ramp;

        if (n % 10 == 0) {
            compute_rho(rho, phi, vel, Nx, dx, m2);
            double beta = dt * 10.0 / tau_avg2;
            if (beta > 1.0) beta = 1.0;
            for (int i = 0; i < Nx; i++)
                rho_run2[i] = (1.0 - beta) * rho_run2[i] + beta * rho[i];
        }

        /* Velocity Verlet: scalars */
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

        /* Velocity Verlet: telegraph Phi */
        for (int i = 1; i < Nx - 1; i++)
            Phi_vel[i] += 0.5 * dt * Phi_acc[i];
        for (int i = 1; i < Nx - 1; i++)
            Phi[i] += dt * Phi_vel[i];
        for (int i = 1; i < Nx - 1; i++) {
            double lapl_Phi = (Phi[i+1] - 2.0*Phi[i] + Phi[i-1]) / dx2;
            Phi_acc[i] = lapl_Phi + alpha_eff * rho_run2[i];
        }
        Phi_acc[0] = Phi_acc[Nx-1] = 0;
        for (int i = 1; i < Nx - 1; i++) {
            Phi_vel[i] += 0.5 * dt * Phi_acc[i];
            Phi_vel[i] /= (1.0 + gamma_star * dt);
        }
        for (int i = 0; i < Nx; i++) {
            Phi[i]     *= damp[i];
            Phi_vel[i] *= damp[i];
        }
        Phi[0] = Phi[Nx-1] = 0;
        Phi_vel[0] = Phi_vel[Nx-1] = 0;
    }

    printf("  Phi(0) after settle = %.6e\n", Phi[ic]);

    /* Record Phi at x=+50 before boost */
    double x_probe = 50.0;
    int i_probe = ic + (int)(x_probe / dx);
    if (i_probe >= Nx) i_probe = Nx - 2;
    double Phi_probe_before = Phi[i_probe];
    printf("  Phi(x=%.0f) before boost = %.6e\n", x_probe, Phi_probe_before);

    /* Now boost the oscillon: shift phi_a and vel_a by +v_boost */
    double v_boost = 0.3;  /* moderate speed */
    printf("  Boosting oscillon with v = %.2f c...\n", v_boost);

    /* Apply Lorentz boost: phi(x) -> phi(gamma_L*(x - v*0)), v -> gamma_L*(v - v_boost*dphi/dx) */
    /* Simpler: just add velocity kick v_boost * dphi/dx to each field */
    for (int a = 0; a < 3; a++)
        for (int i = 1; i < Nx - 1; i++) {
            double dphi_dx = (phi[a][i+1] - phi[a][i-1]) / (2.0 * dx);
            vel[a][i] += -v_boost * dphi_dx;
        }

    /* Evolve after boost and track Phi at x=+50 */
    double t_causal = 200.0;  /* enough for signal at c=1 to reach x=50 */
    int Nt_causal = (int)(t_causal / dt) + 1;

    char causepath[600];
    snprintf(causepath, sizeof(causepath), "%s/telegraph_causality.tsv", outdir);
    FILE *fcause = fopen(causepath, "w");
    fprintf(fcause, "time\tPhi_probe\tPhi_0\tphi1_0\tphi1_probe\n");

    int rec_cause = Nt_causal / 5000;
    if (rec_cause < 1) rec_cause = 1;
    int print_cause = Nt_causal / 20;
    if (print_cause < 1) print_cause = 1;

    /* Track when Phi at probe deviates significantly from pre-boost value */
    double Phi_threshold = fabs(Phi_probe_before) * 0.1;  /* 10% deviation */
    if (Phi_threshold < 1e-15) Phi_threshold = 1e-10;
    double t_arrival = -1.0;

    for (int n = 0; n <= Nt_causal; n++) {
        double t = n * dt;

        /* Update rho_run */
        if (n % 10 == 0) {
            compute_rho(rho, phi, vel, Nx, dx, m2);
            double beta = dt * 10.0 / tau_avg2;
            if (beta > 1.0) beta = 1.0;
            for (int i = 0; i < Nx; i++)
                rho_run2[i] = (1.0 - beta) * rho_run2[i] + beta * rho[i];
        }

        double Phi_probe_now = Phi[i_probe];
        double delta_Phi = fabs(Phi_probe_now - Phi_probe_before);

        if (t_arrival < 0 && delta_Phi > Phi_threshold && t > 1.0) {
            t_arrival = t;
            printf("  ** Signal arrival at x=%.0f: t=%.2f (expected t=%.0f/c=%.0f) **\n",
                   x_probe, t_arrival, x_probe, x_probe);
        }

        if (n % rec_cause == 0)
            fprintf(fcause, "%.4f\t%.6e\t%.6e\t%.6e\t%.6e\n",
                    t, Phi[i_probe], Phi[ic], phi[0][ic], phi[0][i_probe]);

        if (n % print_cause == 0)
            printf("  t=%7.2f  Phi(0)=%.4e  Phi(%.0f)=%.4e  dPhi=%.4e\n",
                   t, Phi[ic], x_probe, Phi[i_probe], delta_Phi);

        if (n == Nt_causal) break;

        /* Velocity Verlet: scalars */
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

        /* Velocity Verlet: telegraph Phi */
        for (int i = 1; i < Nx - 1; i++)
            Phi_vel[i] += 0.5 * dt * Phi_acc[i];
        for (int i = 1; i < Nx - 1; i++)
            Phi[i] += dt * Phi_vel[i];
        for (int i = 1; i < Nx - 1; i++) {
            double lapl_Phi = (Phi[i+1] - 2.0*Phi[i] + Phi[i-1]) / dx2;
            Phi_acc[i] = lapl_Phi + alpha * rho_run2[i];
        }
        Phi_acc[0] = Phi_acc[Nx-1] = 0;
        for (int i = 1; i < Nx - 1; i++) {
            Phi_vel[i] += 0.5 * dt * Phi_acc[i];
            Phi_vel[i] /= (1.0 + gamma_star * dt);
        }
        for (int i = 0; i < Nx; i++) {
            Phi[i]     *= damp[i];
            Phi_vel[i] *= damp[i];
        }
        Phi[0] = Phi[Nx-1] = 0;
        Phi_vel[0] = Phi_vel[Nx-1] = 0;
    }
    fclose(fcause);

    if (t_arrival > 0) {
        printf("\n  Propagation delay: %.2f code units (expected %.0f for c=1)\n",
               t_arrival, x_probe);
        printf("  Effective speed: c_eff = %.0f / %.2f = %.4f\n",
               x_probe, t_arrival, x_probe / t_arrival);
    } else {
        printf("\n  No clear signal arrival detected at x=%.0f within t=%.0f\n",
               x_probe, t_causal);
    }

    printf("\n  Phase 2 complete. Output: %s\n", causepath);
    printf("  Summary: %s\n", sumpath);

    /* Cleanup */
    for (int a = 0; a < 3; a++) {
        free(phi[a]); free(vel[a]); free(acc[a]);
        free(phi_save[a]); free(vel_save[a]);
    }
    free(Phi); free(Phi_vel); free(Phi_acc);
    free(rho); free(rho_avg); free(Phi_poisson);
    free(damp); free(rho_run2);
    return 0;
}
