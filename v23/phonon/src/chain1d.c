/*
 * chain1d.c — V23-D Phase 2: N-oscillon chain phonon spectrum
 *
 * Step 1: Equilibrate single oscillon (t=5000), save at breathing peak.
 * Step 2: Place N=8 oscillons at spacing d on PERIODIC grid (ring topology).
 *         Periodic BC eliminates edge effects / radiation pressure asymmetry.
 * Step 3: Evolve for t_chain=10000.
 * Step 4: Track oscillon positions via energy-weighted centroids.
 * Step 5: Compute phonon dispersion from normal mode decomposition + DFT.
 *
 * The periodic grid has length L = N_osc * d. Grid index wraps: i+1 at i=Nx-1
 * maps to i=0. This eliminates absorbing boundaries for the chain run (still
 * uses absorbing BC for the single-oscillon equilibration).
 *
 * Compile: gcc -O3 -Wall -o chain1d src/chain1d.c -lm
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* --- Parameters --- */
static double mu    = -20.0;
static double kappa = 20.0;
static double mass  = 1.0;
static double A_init = 1.0;
static double sig   = 3.0;
static char   outdir[512] = "data";

/* Phase 1a: single oscillon equilibration */
static int    Nx_eq   = 4000;
static double xmax_eq = 100.0;
static double t_equil = 5000.0;

/* Phase 2: chain */
static int    N_osc    = 8;
static double d_space  = 16.0;
static double t_chain  = 10000.0;
static double t_discard = 2000.0;
static double diag_dt  = 10.0;
static double delta_max = 0.5;
static unsigned int rng_seed = 42;
static int    pts_per_d = 320;  /* grid points per inter-oscillon spacing */
static double t_pre_equil = 0;  /* chain pre-equilibration time (0 = skip) */
static double gamma_pre = 0.001; /* pre-equilibration damping rate */
static double A_kick = 0.01;     /* phonon kick amplitude */

/* Saved equilibrated profile */
static double *eq_phi[3], *eq_vel[3];
static int    eq_Nx = 0;
static double eq_dx = 0;
static double eq_xmax = 0;
static double eq_mass_osc = 0;

/* -dV/dphi_a where V = (mu/2)P^2/(1+kappa P^2), P = phi1 phi2 phi3 */
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

/* Energy density at grid point i (non-periodic, for equilibration) */
static double energy_density_np(double *phi0, double *phi1, double *phi2,
                                double *v0, double *v1, double *v2,
                                int i, int Nx, double dx)
{
    double m2 = mass * mass;
    double e = 0;
    double *phi_a[3] = {phi0, phi1, phi2};
    double *vel_a[3] = {v0, v1, v2};

    for (int a = 0; a < 3; a++) {
        e += 0.5 * vel_a[a][i] * vel_a[a][i];
        if (i > 0 && i < Nx - 1) {
            double dp = (phi_a[a][i+1] - phi_a[a][i-1]) / (2.0 * dx);
            e += 0.5 * dp * dp;
        }
        e += 0.5 * m2 * phi_a[a][i] * phi_a[a][i];
    }

    double P = phi0[i] * phi1[i] * phi2[i];
    double P2 = P * P;
    e += 0.5 * mu * P2 / (1.0 + kappa * P2);

    return e;
}

/* Energy density at grid point i (periodic) */
static double energy_density_p(double *phi0, double *phi1, double *phi2,
                               double *v0, double *v1, double *v2,
                               int i, int Nx, double dx)
{
    double m2 = mass * mass;
    double e = 0;
    double *phi_a[3] = {phi0, phi1, phi2};
    double *vel_a[3] = {v0, v1, v2};

    int ip = (i + 1) % Nx;
    int im = (i - 1 + Nx) % Nx;

    for (int a = 0; a < 3; a++) {
        e += 0.5 * vel_a[a][i] * vel_a[a][i];
        double dp = (phi_a[a][ip] - phi_a[a][im]) / (2.0 * dx);
        e += 0.5 * dp * dp;
        e += 0.5 * m2 * phi_a[a][i] * phi_a[a][i];
    }

    double P = phi0[i] * phi1[i] * phi2[i];
    double P2 = P * P;
    e += 0.5 * mu * P2 / (1.0 + kappa * P2);

    return e;
}

/* Simple LCG random in [0,1) */
static double rand01(unsigned int *state)
{
    *state = *state * 1103515245u + 12345u;
    return (double)((*state >> 16) & 0x7fff) / 32768.0;
}

/* ===================================================================
 *  Phase 1a: Equilibrate single oscillon (identical to phonon1d.c)
 * =================================================================== */
static void run_equilibrate(void)
{
    printf("===== Phase 1a: Single oscillon equilibration =====\n");

    int Nx = Nx_eq;
    double xmax = xmax_eq;
    double dx  = 2.0 * xmax / (Nx - 1);
    double dx2 = dx * dx;
    double m2  = mass * mass;

    double kmax = M_PI / dx;
    double dt   = 0.8 * 2.0 / sqrt(kmax * kmax + m2);
    int    Nt   = (int)(t_equil / dt) + 1;

    printf("  Nx=%d xmax=%.1f dx=%.5f dt=%.6f Nt=%d tfinal=%.0f\n",
           Nx, xmax, dx, dt, Nt, t_equil);

    double *phi[3], *vel[3], *acc[3];
    for (int a = 0; a < 3; a++) {
        phi[a] = calloc(Nx, sizeof(double));
        vel[a] = calloc(Nx, sizeof(double));
        acc[a] = calloc(Nx, sizeof(double));
    }

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

    int ic = Nx / 2;
    for (int a = 0; a < 3; a++)
        for (int i = 0; i < Nx; i++) {
            double x = -xmax + i * dx;
            phi[a][i] = A_init * exp(-x * x / (2.0 * sig * sig));
        }

    #define COMPUTE_ACC_EQ() do { \
        for (int a = 0; a < 3; a++) { \
            acc[a][0] = acc[a][1] = acc[a][Nx-2] = acc[a][Nx-1] = 0; \
            for (int ii = 1; ii < Nx - 1; ii++) { \
                double lapl = (phi[a][ii+1] - 2.0*phi[a][ii] + phi[a][ii-1]) / dx2; \
                double fp = force_pot(phi[0][ii], phi[1][ii], phi[2][ii], a); \
                acc[a][ii] = lapl - m2*phi[a][ii] + fp; \
            } \
        } \
    } while(0)

    COMPUTE_ACC_EQ();

    double prev_phi0 = phi[0][ic];
    double prev_prev_phi0 = phi[0][ic];
    int print_every = Nt / 20;
    if (print_every < 1) print_every = 1;

    double *save_phi[3], *save_vel[3];
    for (int a = 0; a < 3; a++) {
        save_phi[a] = calloc(Nx, sizeof(double));
        save_vel[a] = calloc(Nx, sizeof(double));
    }
    double save_t = -1;
    int n_maxima = 0;
    double t_start_save = 0.5 * t_equil;

    for (int n = 0; n <= Nt; n++) {
        double t = n * dt;

        if (n > 2 && t > t_start_save) {
            double cur = phi[0][ic];
            if (prev_phi0 > prev_prev_phi0 && prev_phi0 > cur && prev_phi0 > 0.1) {
                for (int a = 0; a < 3; a++) {
                    memcpy(save_phi[a], phi[a], Nx * sizeof(double));
                    memcpy(save_vel[a], vel[a], Nx * sizeof(double));
                }
                save_t = t - dt;
                n_maxima++;
            }
        }

        if (n % print_every == 0) {
            double Etot = 0;
            for (int i = 1; i < Nx - 1; i++)
                Etot += energy_density_np(phi[0], phi[1], phi[2],
                                          vel[0], vel[1], vel[2], i, Nx, dx) * dx;
            printf("  t=%7.1f  phi0=%+.4f  E=%.4f  maxima=%d\n",
                   t, phi[0][ic], Etot, n_maxima);
        }

        prev_prev_phi0 = prev_phi0;
        prev_phi0 = phi[0][ic];

        if (n == Nt) break;

        for (int a = 0; a < 3; a++)
            for (int i = 1; i < Nx - 1; i++)
                vel[a][i] += 0.5 * dt * acc[a][i];
        for (int a = 0; a < 3; a++)
            for (int i = 1; i < Nx - 1; i++)
                phi[a][i] += dt * vel[a][i];
        COMPUTE_ACC_EQ();
        for (int a = 0; a < 3; a++)
            for (int i = 1; i < Nx - 1; i++)
                vel[a][i] += 0.5 * dt * acc[a][i];

        for (int a = 0; a < 3; a++)
            for (int i = 0; i < Nx; i++) {
                vel[a][i] *= damp[i];
                phi[a][i] *= damp[i];
            }
    }

    if (n_maxima == 0) {
        printf("ERROR: No breathing maxima found! Using final profile.\n");
        for (int a = 0; a < 3; a++) {
            memcpy(save_phi[a], phi[a], Nx * sizeof(double));
            memcpy(save_vel[a], vel[a], Nx * sizeof(double));
        }
        save_t = t_equil;
    }

    printf("  Saved profile at t=%.1f (%d maxima)\n", save_t, n_maxima);

    eq_mass_osc = 0;
    for (int i = 1; i < Nx - 1; i++)
        eq_mass_osc += energy_density_np(save_phi[0], save_phi[1], save_phi[2],
                                          save_vel[0], save_vel[1], save_vel[2],
                                          i, Nx, dx) * dx;
    printf("  Oscillon mass = %.4f\n", eq_mass_osc);

    eq_Nx = Nx;
    eq_dx = dx;
    eq_xmax = xmax;
    for (int a = 0; a < 3; a++) {
        eq_phi[a] = save_phi[a];
        eq_vel[a] = save_vel[a];
    }

    for (int a = 0; a < 3; a++) { free(phi[a]); free(vel[a]); free(acc[a]); }
    free(damp);
    printf("  Phase 1a complete.\n\n");
}

/* ===================================================================
 *  Initialize chain on periodic domain [0, L)
 * =================================================================== */
static void init_chain(double *phi[3], double *vel[3],
                       int Nx, double dx, double L,
                       double *x_osc, int N)
{
    for (int a = 0; a < 3; a++) {
        memset(phi[a], 0, Nx * sizeof(double));
        memset(vel[a], 0, Nx * sizeof(double));
    }

    for (int n = 0; n < N; n++) {
        double xc = x_osc[n];
        for (int a = 0; a < 3; a++) {
            for (int i = 0; i < Nx; i++) {
                double x = i * dx;
                /* Distance on periodic domain */
                double x_rel = x - xc;
                /* Wrap to [-L/2, L/2) */
                while (x_rel > L / 2) x_rel -= L;
                while (x_rel < -L / 2) x_rel += L;

                /* Interpolate from equilibrated profile (centered at 0) */
                double j_f = (x_rel + eq_xmax) / eq_dx;
                if (j_f >= 0 && j_f < eq_Nx - 1) {
                    int j0 = (int)j_f;
                    double frac = j_f - j0;
                    phi[a][i] += (1 - frac) * eq_phi[a][j0] + frac * eq_phi[a][j0 + 1];
                    vel[a][i] += (1 - frac) * eq_vel[a][j0] + frac * eq_vel[a][j0 + 1];
                }
            }
        }
    }
}

/* ===================================================================
 *  Track oscillon positions on periodic domain
 *  Energy-weighted centroid within each Voronoi cell
 * =================================================================== */
static void find_positions_periodic(double *phi0, double *phi1, double *phi2,
                                     double *v0, double *v1, double *v2,
                                     int Nx, double dx, double L,
                                     double *x_osc, int N)
{
    for (int n = 0; n < N; n++) {
        /* Voronoi cell: points closer to x_osc[n] than to any neighbor.
         * On a 1D ring, this is [midpoint_left, midpoint_right).
         * We compute energy-weighted centroid using angular method to handle
         * periodic wrapping: theta = 2*pi*x/L, centroid via atan2. */

        double x_left = x_osc[(n - 1 + N) % N];
        double x_right = x_osc[(n + 1) % N];
        double xc = x_osc[n];

        /* Half-distances to neighbors (on periodic ring) */
        double dl = xc - x_left;
        while (dl < 0) dl += L;
        while (dl > L) dl -= L;
        double dr = x_right - xc;
        while (dr < 0) dr += L;
        while (dr > L) dr -= L;

        /* Cell boundaries relative to xc */
        double cell_lo = -dl / 2.0;
        double cell_hi = dr / 2.0;

        /* Accumulate centroid using displacement from xc */
        double num = 0, den = 0;
        for (int i = 0; i < Nx; i++) {
            double x = i * dx;
            double rel = x - xc;
            while (rel > L / 2) rel -= L;
            while (rel < -L / 2) rel += L;

            if (rel >= cell_lo && rel < cell_hi) {
                double e = energy_density_p(phi0, phi1, phi2, v0, v1, v2, i, Nx, dx);
                if (e > 0) {
                    num += rel * e;
                    den += e;
                }
            }
        }

        if (den > 1e-30) {
            double shift = num / den;
            double new_x = xc + shift;
            /* Wrap to [0, L) */
            while (new_x < 0) new_x += L;
            while (new_x >= L) new_x -= L;
            x_osc[n] = new_x;
        }
    }
}

/* ===================================================================
 *  Main
 * =================================================================== */
int main(int argc, char **argv)
{
    for (int i = 1; i < argc - 1; i += 2) {
        if      (!strcmp(argv[i], "-mu"))       mu       = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-kappa"))    kappa    = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-mass"))     mass     = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-A"))        A_init   = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-sigma"))    sig      = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-o"))        strncpy(outdir, argv[i+1], sizeof(outdir)-1);
        else if (!strcmp(argv[i], "-N"))        N_osc    = atoi(argv[i+1]);
        else if (!strcmp(argv[i], "-d"))        d_space  = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-t_chain"))  t_chain  = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-t_discard")) t_discard = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-seed"))     rng_seed = (unsigned int)atoi(argv[i+1]);
        else if (!strcmp(argv[i], "-diag_dt"))  diag_dt  = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-pts_per_d")) pts_per_d = atoi(argv[i+1]);
        else if (!strcmp(argv[i], "-t_pre"))    t_pre_equil = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-gamma"))    gamma_pre = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-kick"))     A_kick = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-delta"))    delta_max = atof(argv[i+1]);
        else { fprintf(stderr, "Unknown flag: %s\n", argv[i]); exit(1); }
    }

    printf("=== V23-D Phase 2: %d-Oscillon Chain (Periodic) ===\n\n", N_osc);

    /* ===== Step 1: Equilibrate ===== */
    run_equilibrate();

    /* ===== Step 2: Initialize periodic chain ===== */
    printf("===== Step 2: Initialize %d-oscillon periodic chain =====\n", N_osc);

    double L = N_osc * d_space;  /* periodic box length */
    int Nx = N_osc * pts_per_d;  /* grid points */
    double dx = L / Nx;
    double dx2 = dx * dx;
    double m2 = mass * mass;

    double kmax_c = M_PI / dx;
    double dt = 0.8 * 2.0 / sqrt(kmax_c * kmax_c + m2);
    int Nt = (int)(t_chain / dt) + 1;

    printf("  L=%.1f Nx=%d dx=%.5f dt=%.6f Nt=%d\n", L, Nx, dx, dt, Nt);
    printf("  d=%.1f delta_max=%.2f seed=%u\n", d_space, delta_max, rng_seed);

    /* Oscillon positions on [0, L) */
    double *x_osc = malloc(N_osc * sizeof(double));
    double *x_osc_ref = malloc(N_osc * sizeof(double));  /* reference positions */
    unsigned int rng = rng_seed;

    printf("  Oscillon positions:\n");
    for (int n = 0; n < N_osc; n++) {
        double xc = n * d_space + d_space / 2.0;  /* centered in each cell */
        double delta = delta_max * (2.0 * rand01(&rng) - 1.0);
        x_osc[n] = xc + delta;
        x_osc_ref[n] = x_osc[n];  /* will be updated at t_discard */
        printf("    n=%d: x=%.3f (nominal=%.1f, delta=%+.3f)\n",
               n, x_osc[n], xc, delta);
    }

    /* Allocate fields */
    double *phi[3], *vel[3], *acc[3];
    for (int a = 0; a < 3; a++) {
        phi[a] = calloc(Nx, sizeof(double));
        vel[a] = calloc(Nx, sizeof(double));
        acc[a] = calloc(Nx, sizeof(double));
    }

    init_chain(phi, vel, Nx, dx, L, x_osc, N_osc);

    /* Compute initial energy */
    double E_init = 0;
    for (int i = 0; i < Nx; i++)
        E_init += energy_density_p(phi[0], phi[1], phi[2],
                                    vel[0], vel[1], vel[2], i, Nx, dx) * dx;
    printf("  Initial energy: %.4f (expected: %d * %.4f = %.4f)\n",
           E_init, N_osc, eq_mass_osc, N_osc * eq_mass_osc);

    /* ===== Step 3: Evolve with periodic BC ===== */
    printf("\n===== Step 3: Evolving chain (periodic BC) for t=%.0f =====\n", t_chain);

    /* Compute acceleration with periodic Laplacian */
    #define COMPUTE_ACC_P() do { \
        for (int a = 0; a < 3; a++) { \
            for (int ii = 0; ii < Nx; ii++) { \
                int ip1 = (ii + 1) % Nx; \
                int im1 = (ii - 1 + Nx) % Nx; \
                double lapl = (phi[a][ip1] - 2.0*phi[a][ii] + phi[a][im1]) / dx2; \
                double fp = force_pot(phi[0][ii], phi[1][ii], phi[2][ii], a); \
                acc[a][ii] = lapl - m2*phi[a][ii] + fp; \
            } \
        } \
    } while(0)

    COMPUTE_ACC_P();

    /* ===== Chain pre-equilibration with damping ===== */
    {
        double t_pre = t_pre_equil;
        double gamma_damp = gamma_pre;
        int Nt_pre = (int)(t_pre / dt) + 1;
        int print_pre = Nt_pre / 10;
        if (print_pre < 1) print_pre = 1;

        printf("\n  Pre-equilibrating chain with damping (gamma=%.4f, t=%.0f)...\n",
               gamma_damp, t_pre);



        for (int step = 0; step < Nt_pre; step++) {
            for (int a = 0; a < 3; a++)
                for (int i = 0; i < Nx; i++)
                    vel[a][i] += 0.5 * dt * acc[a][i];
            for (int a = 0; a < 3; a++)
                for (int i = 0; i < Nx; i++)
                    phi[a][i] += dt * vel[a][i];
            COMPUTE_ACC_P();
            for (int a = 0; a < 3; a++)
                for (int i = 0; i < Nx; i++)
                    vel[a][i] += 0.5 * dt * acc[a][i];

            /* Apply position-dependent damping: damp strongly between
             * oscillons (where field is small), weakly inside cores.
             * This removes inter-oscillon radiation while preserving breathing. */
            for (int a = 0; a < 3; a++)
                for (int i = 0; i < Nx; i++) {
                    /* Field amplitude squared at this point */
                    double p2 = phi[0][i]*phi[0][i] + phi[1][i]*phi[1][i]
                              + phi[2][i]*phi[2][i];
                    /* Protection factor: 1 inside core (large p2), 0 in vacuum */
                    double protect = p2 / (p2 + 0.01);
                    /* Effective damping: strong in vacuum, weak in core */
                    double local_damp = 1.0 - gamma_damp * dt * (1.0 - protect);
                    vel[a][i] *= local_damp;
                }

            if (step % print_pre == 0) {
                double Etot = 0;
                for (int i = 0; i < Nx; i++)
                    Etot += energy_density_p(phi[0], phi[1], phi[2],
                                              vel[0], vel[1], vel[2], i, Nx, dx) * dx;
                printf("    pre t=%7.1f  E=%.4f\n", step * dt, Etot);
            }
        }

        /* Compute post-damping energy */
        double E_post = 0;
        for (int i = 0; i < Nx; i++)
            E_post += energy_density_p(phi[0], phi[1], phi[2],
                                        vel[0], vel[1], vel[2], i, Nx, dx) * dx;
        printf("  Pre-equilibration done. E_post=%.4f\n", E_post);

        /* Track positions after pre-equilibration */
        find_positions_periodic(phi[0], phi[1], phi[2],
                                 vel[0], vel[1], vel[2],
                                 Nx, dx, L, x_osc, N_osc);
        printf("  Post-damping positions:\n");
        for (int n = 0; n < N_osc; n++)
            printf("    n=%d: x=%.3f\n", n, x_osc[n]);

        /* Now apply small velocity kicks to seed phonon modes.
         * Add a q=1 mode perturbation: v_kick(x) ~ sin(2*pi*x/L) * A_kick
         * This kicks each oscillon proportional to sin(2*pi*x_n/L). */
        double kick_amp = A_kick;
        printf("  Applying phonon kick (A=%.4f)...\n", A_kick);
        for (int i = 0; i < Nx; i++) {
            double x = i * dx;
            double kick = kick_amp * sin(2.0 * M_PI * x / L);
            for (int a = 0; a < 3; a++)
                vel[a][i] += kick * phi[a][i];
            /* Scale kick by phi so only oscillon cores are kicked, not vacuum */
        }

        /* Re-compute acceleration */
        COMPUTE_ACC_P();
    }

    printf("\n===== Step 3: Evolving chain (conservative) for t=%.0f =====\n", t_chain);

    /* Diagnostic storage */
    int max_diag = (int)(t_chain / diag_dt) + 10;
    double *diag_t = malloc(max_diag * sizeof(double));
    double **diag_x = malloc(N_osc * sizeof(double*));
    for (int n = 0; n < N_osc; n++)
        diag_x[n] = malloc(max_diag * sizeof(double));
    double *diag_E = malloc(max_diag * sizeof(double));
    double *diag_maxdisp = malloc(max_diag * sizeof(double));
    double *diag_mean_spacing = malloc(max_diag * sizeof(double));
    int n_diag = 0;

    int diag_every = (int)(diag_dt / dt);
    if (diag_every < 1) diag_every = 1;
    int print_every = Nt / 40;
    if (print_every < 1) print_every = 1;

    /* Reference positions updated at t_discard */
    int ref_set = 0;

    for (int step = 0; step <= Nt; step++) {
        double t = step * dt;

        /* Diagnostic */
        if (step % diag_every == 0 && n_diag < max_diag) {
            /* Track positions (only every diag_every to save CPU) */
            if (step > 0)
                find_positions_periodic(phi[0], phi[1], phi[2],
                                         vel[0], vel[1], vel[2],
                                         Nx, dx, L, x_osc, N_osc);

            /* Set reference at t_discard */
            if (!ref_set && t >= t_discard) {
                for (int n = 0; n < N_osc; n++)
                    x_osc_ref[n] = x_osc[n];
                ref_set = 1;
                printf("  ** Reference positions set at t=%.1f **\n", t);
            }

            /* Energy */
            double Etot = 0;
            for (int i = 0; i < Nx; i++)
                Etot += energy_density_p(phi[0], phi[1], phi[2],
                                          vel[0], vel[1], vel[2], i, Nx, dx) * dx;

            /* Mean spacing and max displacement */
            double ms = 0, md = 0;
            for (int n = 0; n < N_osc; n++) {
                int nn = (n + 1) % N_osc;
                double sp = x_osc[nn] - x_osc[n];
                while (sp < 0) sp += L;
                while (sp > L) sp -= L;
                ms += sp;

                double disp = x_osc[n] - x_osc_ref[n];
                while (disp > L / 2) disp -= L;
                while (disp < -L / 2) disp += L;
                if (fabs(disp) > md) md = fabs(disp);
            }
            ms /= N_osc;

            diag_t[n_diag] = t;
            for (int n = 0; n < N_osc; n++)
                diag_x[n][n_diag] = x_osc[n];
            diag_E[n_diag] = Etot;
            diag_mean_spacing[n_diag] = ms;
            diag_maxdisp[n_diag] = md;
            n_diag++;

            if (step % print_every == 0) {
                printf("  t=%7.1f  E=%.4f  <d>=%.3f  max_u=%.4f\n",
                       t, Etot, ms, md);
            }
        }

        if (step == Nt) break;

        /* Velocity Verlet with periodic BC */
        for (int a = 0; a < 3; a++)
            for (int i = 0; i < Nx; i++)
                vel[a][i] += 0.5 * dt * acc[a][i];
        for (int a = 0; a < 3; a++)
            for (int i = 0; i < Nx; i++)
                phi[a][i] += dt * vel[a][i];
        COMPUTE_ACC_P();
        for (int a = 0; a < 3; a++)
            for (int i = 0; i < Nx; i++)
                vel[a][i] += 0.5 * dt * acc[a][i];
    }

    printf("\n  Evolution complete. %d diagnostic samples.\n", n_diag);

    /* ===== Step 4: Output data ===== */
    printf("\n===== Step 4: Writing output files =====\n");

    /* chain_positions.tsv */
    {
        char path[600];
        snprintf(path, sizeof(path), "%s/chain_positions.tsv", outdir);
        FILE *fp = fopen(path, "w");
        if (fp) {
            fprintf(fp, "time");
            for (int n = 0; n < N_osc; n++) fprintf(fp, "\tX_%d", n+1);
            fprintf(fp, "\n");
            for (int j = 0; j < n_diag; j++) {
                fprintf(fp, "%.4f", diag_t[j]);
                for (int n = 0; n < N_osc; n++)
                    fprintf(fp, "\t%.6f", diag_x[n][j]);
                fprintf(fp, "\n");
            }
            fclose(fp);
            printf("  Wrote %s (%d rows)\n", path, n_diag);
        }
    }

    /* chain_ts.tsv */
    {
        char path[600];
        snprintf(path, sizeof(path), "%s/chain_ts.tsv", outdir);
        FILE *fp = fopen(path, "w");
        if (fp) {
            fprintf(fp, "time\tE_total\tmean_spacing\tmax_displacement\n");
            for (int j = 0; j < n_diag; j++)
                fprintf(fp, "%.4f\t%.6e\t%.6f\t%.6e\n",
                        diag_t[j], diag_E[j], diag_mean_spacing[j], diag_maxdisp[j]);
            fclose(fp);
            printf("  Wrote %s\n", path);
        }
    }

    /* ===== Steps 5-6: Phonon spectrum ===== */
    printf("\n===== Steps 5-6: Phonon spectrum =====\n");

    int j_start = 0;
    for (int j = 0; j < n_diag; j++) {
        if (diag_t[j] >= t_discard) { j_start = j; break; }
    }
    int n_analysis = n_diag - j_start;
    printf("  Analysis: t=[%.0f, %.0f], %d samples\n",
           diag_t[j_start], diag_t[n_diag-1], n_analysis);

    if (n_analysis < 20) {
        printf("ERROR: Too few samples!\n");
        goto cleanup;
    }

    /* Compute displacements u_n(t) on the periodic ring.
     * u_n(t) = x_n(t) - x_n_ref, unwrapped for continuity. */
    double **u_n = malloc(N_osc * sizeof(double*));
    for (int n = 0; n < N_osc; n++) {
        u_n[n] = calloc(n_analysis, sizeof(double));
        double prev_x = diag_x[n][j_start];
        u_n[n][0] = 0;
        for (int j = 1; j < n_analysis; j++) {
            int jj = j + j_start;
            double dx_step = diag_x[n][jj] - prev_x;
            /* Unwrap periodic jumps */
            while (dx_step > L / 2) dx_step -= L;
            while (dx_step < -L / 2) dx_step += L;
            u_n[n][j] = u_n[n][j - 1] + dx_step;
            prev_x = diag_x[n][jj];
        }
    }

    /* Write displacements */
    {
        char path[600];
        snprintf(path, sizeof(path), "%s/chain_displacements.tsv", outdir);
        FILE *fp = fopen(path, "w");
        if (fp) {
            fprintf(fp, "time");
            for (int n = 0; n < N_osc; n++) fprintf(fp, "\tu_%d", n+1);
            fprintf(fp, "\n");
            for (int j = 0; j < n_analysis; j++) {
                fprintf(fp, "%.4f", diag_t[j + j_start]);
                for (int n = 0; n < N_osc; n++)
                    fprintf(fp, "\t%.8e", u_n[n][j]);
                fprintf(fp, "\n");
            }
            fclose(fp);
            printf("  Wrote %s\n", path);
        }
    }

    /* Remove linear drift from each u_n (subtract center-of-mass motion) */
    for (int j = 0; j < n_analysis; j++) {
        double com = 0;
        for (int n = 0; n < N_osc; n++) com += u_n[n][j];
        com /= N_osc;
        for (int n = 0; n < N_osc; n++) u_n[n][j] -= com;
    }

    /* Normal mode decomposition:
     * Q_q(t) = (1/sqrt(N)) sum_n u_n(t) exp(-i 2pi q n / N) */
    double T_analysis = diag_t[n_diag-1] - diag_t[j_start];
    double dt_diag = T_analysis / (n_analysis - 1);

    double **Q_re = malloc(N_osc * sizeof(double*));
    double **Q_im = malloc(N_osc * sizeof(double*));
    for (int q = 0; q < N_osc; q++) {
        Q_re[q] = calloc(n_analysis, sizeof(double));
        Q_im[q] = calloc(n_analysis, sizeof(double));
    }

    for (int j = 0; j < n_analysis; j++) {
        for (int q = 0; q < N_osc; q++) {
            double re = 0, im = 0;
            for (int n = 0; n < N_osc; n++) {
                double phase = -2.0 * M_PI * q * n / N_osc;
                re += u_n[n][j] * cos(phase);
                im += u_n[n][j] * sin(phase);
            }
            Q_re[q][j] = re / sqrt((double)N_osc);
            Q_im[q][j] = im / sqrt((double)N_osc);
        }
    }

    /* Power spectrum for each mode q */
    int n_omega = n_analysis / 2;
    if (n_omega > 4000) n_omega = 4000;
    double omega_nyquist = M_PI / dt_diag;
    /* Focus on low frequencies where phonons live */
    double omega_max = 1.0;  /* well below mass gap */
    if (omega_max > omega_nyquist) omega_max = omega_nyquist;

    double *omega_peak = calloc(N_osc, sizeof(double));
    double *power_peak = calloc(N_osc, sizeof(double));

    /* 2D spectrum output */
    char path2d[600];
    snprintf(path2d, sizeof(path2d), "%s/chain_spectrum_2d.tsv", outdir);
    FILE *f2d = fopen(path2d, "w");
    if (f2d) fprintf(f2d, "k\tomega\tpower\n");

    for (int q = 0; q <= N_osc / 2; q++) {
        double k_q = 2.0 * M_PI * q / (N_osc * d_space);
        double best_pow = 0, best_omega = 0;

        for (int l = 1; l <= n_omega; l++) {
            double omega = omega_max * l / n_omega;

            double dft_re = 0, dft_im = 0;
            for (int j = 0; j < n_analysis; j++) {
                double t_j = j * dt_diag;
                double c = cos(omega * t_j);
                double s = sin(omega * t_j);
                dft_re += (Q_re[q][j] * c + Q_im[q][j] * s) * dt_diag;
                dft_im += (-Q_re[q][j] * s + Q_im[q][j] * c) * dt_diag;
            }
            double power = (dft_re * dft_re + dft_im * dft_im) / (T_analysis * T_analysis);

            if (f2d) fprintf(f2d, "%.6f\t%.6f\t%.8e\n", k_q, omega, power);

            if (power > best_pow) {
                best_pow = power;
                best_omega = omega;
            }
        }

        omega_peak[q] = best_omega;
        power_peak[q] = best_pow;
        printf("  Mode q=%d: k=%.5f  omega_peak=%.5f  power=%.4e\n",
               q, k_q, best_omega, best_pow);
    }

    if (f2d) { fclose(f2d); printf("  Wrote %s\n", path2d); }

    /* chain_modes.tsv */
    {
        char path[600];
        snprintf(path, sizeof(path), "%s/chain_modes.tsv", outdir);
        FILE *fp = fopen(path, "w");
        if (fp) {
            fprintf(fp, "q\tk\tomega_peak\tpower_peak\n");
            for (int q = 0; q <= N_osc / 2; q++) {
                double k_q = 2.0 * M_PI * q / (N_osc * d_space);
                fprintf(fp, "%d\t%.6f\t%.6f\t%.8e\n",
                        q, k_q, omega_peak[q], power_peak[q]);
            }
            fclose(fp);
            printf("  Wrote %s\n", path);
        }
    }

    /* ===== Summary ===== */
    printf("\n===== SUMMARY =====\n");
    printf("  Chain: %d oscillons, spacing d=%.1f, periodic L=%.1f\n",
           N_osc, d_space, L);
    printf("  Oscillon mass M = %.4f\n", eq_mass_osc);
    printf("  Initial energy: %.4f, Final: %.4f (%.1f%%)\n",
           E_init, diag_E[n_diag-1],
           100.0 * diag_E[n_diag-1] / E_init);
    printf("  Mean spacing: initial=%.3f, final=%.3f\n",
           diag_mean_spacing[0], diag_mean_spacing[n_diag-1]);
    printf("  Max displacement at end: %.4f\n", diag_maxdisp[n_diag-1]);

    printf("\n  Dispersion relation:\n");
    printf("  %-4s  %-10s  %-10s  %-10s\n", "q", "k", "omega", "c_s=w/k");
    printf("  ----  ----------  ----------  ----------\n");
    for (int q = 0; q <= N_osc / 2; q++) {
        double k_q = 2.0 * M_PI * q / (N_osc * d_space);
        double cs = (k_q > 1e-10) ? omega_peak[q] / k_q : 0;
        printf("  %-4d  %-10.5f  %-10.5f  %-10.4f\n",
               q, k_q, omega_peak[q], cs);
    }

    if (N_osc >= 2) {
        double k1 = 2.0 * M_PI / (N_osc * d_space);
        double cs_est = omega_peak[1] / k1;
        printf("\n  Sound speed (q=1): c_s = %.4f\n", cs_est);
        printf("  Phase 1 prediction: c_s ~ 0.44\n");

        /* Check stability: is max displacement < d/4 ? */
        double max_u_final = diag_maxdisp[n_diag-1];
        printf("\n  Stability: max_u/d = %.3f",
               max_u_final / d_space);
        if (max_u_final < d_space / 4)
            printf(" (STABLE)\n");
        else if (max_u_final < d_space / 2)
            printf(" (MARGINAL)\n");
        else
            printf(" (UNSTABLE)\n");
    }

cleanup:
    for (int a = 0; a < 3; a++) {
        free(phi[a]); free(vel[a]); free(acc[a]);
    }
    free(x_osc); free(x_osc_ref);
    for (int n = 0; n < N_osc; n++) free(diag_x[n]);
    free(diag_x); free(diag_t); free(diag_E);
    free(diag_mean_spacing); free(diag_maxdisp);
    for (int q = 0; q < N_osc; q++) { free(Q_re[q]); free(Q_im[q]); }
    free(Q_re); free(Q_im);
    free(omega_peak); free(power_peak);
    for (int n = 0; n < N_osc; n++) free(u_n[n]);
    free(u_n);
    for (int a = 0; a < 3; a++) {
        if (eq_phi[a]) free(eq_phi[a]);
        if (eq_vel[a]) free(eq_vel[a]);
    }

    printf("\n=== V23-D Phase 2 Complete ===\n");
    return 0;
}
