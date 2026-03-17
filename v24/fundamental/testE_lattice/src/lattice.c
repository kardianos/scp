/*
 * lattice.c — Test E: Pairwise-Coupled Oscillon Lattice
 *
 * 8-oscillon periodic chain WITH pairwise coupling lambda.
 * Combines V23-D chain1d.c (periodic chain + phonon analysis)
 * with V24-ME maxwell_e.c (pairwise coupling).
 *
 * Lagrangian:
 *   L = sum_a [ (1/2)(dt phi_a)^2 - (1/2)(dx phi_a)^2 - (m^2/2)phi_a^2 ]
 *     - (mu/2) P^2 / (1 + kappa P^2)         [triple product]
 *     - lambda(phi_1 phi_2 + phi_2 phi_3 + phi_3 phi_1)  [pairwise]
 *
 * Steps:
 *   1. Equilibrate single oscillon with pairwise coupling (t=5000)
 *   2. Build 8-oscillon periodic chain at spacing d=16
 *   3. Add small random displacements delta in [-0.5, 0.5]
 *   4. Evolve t=10000
 *   5. Track all 8 positions
 *   6. Compute phonon spectrum via normal mode decomposition
 *   7. Measure phase coherence (time crystal test)
 *
 * Compile: gcc -O3 -Wall -o lattice src/lattice.c -lm
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* --- Parameters --- */
static double mu       = -20.0;
static double kappa    = 20.0;
static double mass     = 1.0;
static double lambda   = 0.0;    /* pairwise coupling */
static double A_init   = 1.0;
static double sig      = 3.0;
static char   outdir[512] = "data";

/* Phase 1: single oscillon equilibration */
static int    Nx_eq    = 4000;
static double xmax_eq  = 100.0;
static double t_equil  = 5000.0;

/* Phase 2: chain */
static int    N_osc    = 8;
static double d_space  = 16.0;
static double t_chain  = 10000.0;
static double t_discard = 2000.0;
static double diag_dt  = 10.0;
static double delta_max = 0.5;
static unsigned int rng_seed = 42;
static int    pts_per_d = 320;
static double t_pre_equil = 500;   /* chain pre-equilibration with damping */
static double gamma_pre = 0.001;
static double A_kick = 0.01;

/* Saved equilibrated profile */
static double *eq_phi[3], *eq_vel[3];
static int    eq_Nx = 0;
static double eq_dx = 0;
static double eq_xmax = 0;
static double eq_mass_osc = 0;
static double eq_omega_osc = 0;  /* breathing frequency */

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

    /* Pairwise coupling energy */
    e += lambda * (phi0[i]*phi1[i] + phi1[i]*phi2[i] + phi2[i]*phi0[i]);

    /* Triple product energy */
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

    /* Pairwise coupling energy */
    e += lambda * (phi0[i]*phi1[i] + phi1[i]*phi2[i] + phi2[i]*phi0[i]);

    /* Triple product energy */
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
 *  Phase 1: Equilibrate single oscillon WITH pairwise coupling
 * =================================================================== */
static void run_equilibrate(void)
{
    printf("===== Phase 1: Single oscillon equilibration (lambda=%.4f) =====\n", lambda);

    int Nx = Nx_eq;
    double xmax = xmax_eq;
    double dx  = 2.0 * xmax / (Nx - 1);
    double dx2 = dx * dx;
    double m2  = mass * mass;

    /* CFL: use heaviest mode mass for stability */
    double m2_eff = m2 + 2.0 * lambda;  /* symmetric mode is heaviest */
    if (m2_eff < m2) m2_eff = m2;
    double kmax = M_PI / dx;
    double dt   = 0.8 * 2.0 / sqrt(kmax * kmax + m2_eff);
    int    Nt   = (int)(t_equil / dt) + 1;

    /* Check for tachyonic antisymmetric mode */
    double m2_anti = m2 - lambda;
    if (m2_anti <= 0) {
        printf("  WARNING: m2_anti = %.4f <= 0 (tachyonic!)\n", m2_anti);
    }

    printf("  Nx=%d xmax=%.1f dx=%.5f dt=%.6f Nt=%d tfinal=%.0f\n",
           Nx, xmax, dx, dt, Nt, t_equil);
    printf("  m2_anti=%.4f m_A=%.4f  m2_sym=%.4f m_S=%.4f\n",
           m2_anti, m2_anti > 0 ? sqrt(m2_anti) : 0.0,
           m2 + 2.0*lambda, sqrt(m2 + 2.0*lambda));

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

    /* Initialize: symmetric Gaussians (all three fields identical) */
    for (int a = 0; a < 3; a++)
        for (int i = 0; i < Nx; i++) {
            double x = -xmax + i * dx;
            phi[a][i] = A_init * exp(-x * x / (2.0 * sig * sig));
        }

    /* Acceleration macro with pairwise coupling */
    #define COMPUTE_ACC_EQ() do { \
        for (int a = 0; a < 3; a++) { \
            int b = (a+1)%3, c_idx = (a+2)%3; \
            acc[a][0] = acc[a][1] = acc[a][Nx-2] = acc[a][Nx-1] = 0; \
            for (int ii = 1; ii < Nx - 1; ii++) { \
                double lapl = (phi[a][ii+1] - 2.0*phi[a][ii] + phi[a][ii-1]) / dx2; \
                double fp = force_triple(phi[0][ii], phi[1][ii], phi[2][ii], a); \
                acc[a][ii] = lapl - m2*phi[a][ii] \
                           - lambda*(phi[b][ii] + phi[c_idx][ii]) + fp; \
            } \
        } \
    } while(0)

    COMPUTE_ACC_EQ();

    int ic = Nx / 2;
    double prev_phi0 = phi[0][ic];
    double prev_prev_phi0 = phi[0][ic];
    int print_every = Nt / 20;
    if (print_every < 1) print_every = 1;

    /* Save best profile at breathing peak */
    double *save_phi[3], *save_vel[3];
    for (int a = 0; a < 3; a++) {
        save_phi[a] = calloc(Nx, sizeof(double));
        save_vel[a] = calloc(Nx, sizeof(double));
    }
    double save_t = -1;
    int n_maxima = 0;
    double t_start_save = 0.5 * t_equil;

    /* DFT storage for breathing frequency */
    int max_dft = 50000;
    double *phi0_hist = malloc(max_dft * sizeof(double));
    double *t_hist = malloc(max_dft * sizeof(double));
    int n_dft = 0;
    int dft_every = Nt / max_dft;
    if (dft_every < 1) dft_every = 1;

    for (int n = 0; n <= Nt; n++) {
        double t = n * dt;

        /* Record for DFT */
        if (n % dft_every == 0 && n_dft < max_dft) {
            phi0_hist[n_dft] = phi[0][ic];
            t_hist[n_dft] = t;
            n_dft++;
        }

        /* Save at breathing peak */
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

        /* Velocity Verlet */
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

        /* Absorbing boundary */
        for (int a = 0; a < 3; a++)
            for (int i = 0; i < Nx; i++) {
                vel[a][i] *= damp[i];
                phi[a][i] *= damp[i];
            }
    }
    #undef COMPUTE_ACC_EQ

    if (n_maxima == 0) {
        printf("  WARNING: No breathing maxima found! Using final profile.\n");
        for (int a = 0; a < 3; a++) {
            memcpy(save_phi[a], phi[a], Nx * sizeof(double));
            memcpy(save_vel[a], vel[a], Nx * sizeof(double));
        }
        save_t = t_equil;
    }
    printf("  Saved profile at t=%.1f (%d maxima)\n", save_t, n_maxima);

    /* Compute oscillon mass */
    eq_mass_osc = 0;
    for (int i = 1; i < Nx - 1; i++)
        eq_mass_osc += energy_density_np(save_phi[0], save_phi[1], save_phi[2],
                                          save_vel[0], save_vel[1], save_vel[2],
                                          i, Nx, dx) * dx;
    printf("  Oscillon mass = %.4f\n", eq_mass_osc);

    /* Compute breathing frequency from DFT of second half */
    eq_omega_osc = 0;
    {
        int dft_start = n_dft / 2;
        if (n_dft - dft_start > 100) {
            double T = t_hist[n_dft-1] - t_hist[dft_start];
            int nf = 500;
            double peak_pow = 0, peak_om = 0;
            for (int k = 1; k < nf; k++) {
                double omega = 3.0 * mass * k / nf;
                double re = 0, im = 0;
                for (int j = dft_start; j < n_dft; j++) {
                    double dtj = (j > dft_start) ?
                        (t_hist[j]-t_hist[j-1]) : (t_hist[dft_start+1]-t_hist[dft_start]);
                    re += phi0_hist[j] * cos(omega * t_hist[j]) * dtj;
                    im += phi0_hist[j] * sin(omega * t_hist[j]) * dtj;
                }
                double pw = re*re + im*im;
                if (pw > peak_pow) { peak_pow = pw; peak_om = omega; }
            }
            eq_omega_osc = peak_om;
            printf("  Breathing frequency = %.4f (mass gap = %.4f)\n",
                   eq_omega_osc, mass);
        }
    }

    /* Store equilibrated profile */
    eq_Nx = Nx;
    eq_dx = dx;
    eq_xmax = xmax;
    for (int a = 0; a < 3; a++) {
        eq_phi[a] = save_phi[a];
        eq_vel[a] = save_vel[a];
    }

    for (int a = 0; a < 3; a++) { free(phi[a]); free(vel[a]); free(acc[a]); }
    free(damp);
    free(phi0_hist);
    free(t_hist);
    printf("  Phase 1 complete.\n\n");
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
                double x_rel = x - xc;
                while (x_rel > L / 2) x_rel -= L;
                while (x_rel < -L / 2) x_rel += L;

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
 *  Track oscillon positions: energy-weighted centroid in Voronoi cell
 * =================================================================== */
static void find_positions_periodic(double *phi0, double *phi1, double *phi2,
                                     double *v0, double *v1, double *v2,
                                     int Nx, double dx, double L,
                                     double *x_osc, int N)
{
    for (int n = 0; n < N; n++) {
        double x_left = x_osc[(n - 1 + N) % N];
        double x_right = x_osc[(n + 1) % N];
        double xc = x_osc[n];

        double dl = xc - x_left;
        while (dl < 0) dl += L;
        while (dl > L) dl -= L;
        double dr = x_right - xc;
        while (dr < 0) dr += L;
        while (dr > L) dr -= L;

        double cell_lo = -dl / 2.0;
        double cell_hi = dr / 2.0;

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
        else if (!strcmp(argv[i], "-lambda"))   lambda   = atof(argv[i+1]);
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

    printf("=== Test E: %d-Oscillon Pairwise-Coupled Lattice ===\n", N_osc);
    printf("  lambda=%.4f  mu=%.1f  kappa=%.1f  mass=%.4f\n", lambda, mu, kappa, mass);

    /* ===== Step 1: Equilibrate single oscillon with pairwise coupling ===== */
    run_equilibrate();

    /* ===== Step 2: Initialize periodic chain ===== */
    printf("===== Step 2: Initialize %d-oscillon periodic chain =====\n", N_osc);

    double L = N_osc * d_space;
    int Nx = N_osc * pts_per_d;
    double dx = L / Nx;
    double dx2 = dx * dx;
    double m2 = mass * mass;

    /* CFL with pairwise coupling */
    double m2_eff = m2 + 2.0 * lambda;
    if (m2_eff < m2) m2_eff = m2;
    double kmax_c = M_PI / dx;
    double dt = 0.8 * 2.0 / sqrt(kmax_c * kmax_c + m2_eff);
    int Nt = (int)(t_chain / dt) + 1;

    printf("  L=%.1f Nx=%d dx=%.5f dt=%.6f Nt=%d\n", L, Nx, dx, dt, Nt);
    printf("  d=%.1f delta_max=%.2f seed=%u\n", d_space, delta_max, rng_seed);

    /* Oscillon positions on [0, L) */
    double *x_osc = malloc(N_osc * sizeof(double));
    double *x_osc_ref = malloc(N_osc * sizeof(double));
    unsigned int rng = rng_seed;

    printf("  Oscillon positions:\n");
    for (int n = 0; n < N_osc; n++) {
        double xc = n * d_space + d_space / 2.0;
        double delta = delta_max * (2.0 * rand01(&rng) - 1.0);
        x_osc[n] = xc + delta;
        x_osc_ref[n] = x_osc[n];
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

    /* Compute acceleration with periodic Laplacian + pairwise coupling */
    #define COMPUTE_ACC_P() do { \
        for (int a = 0; a < 3; a++) { \
            int b = (a+1)%3, c_idx = (a+2)%3; \
            for (int ii = 0; ii < Nx; ii++) { \
                int ip1 = (ii + 1) % Nx; \
                int im1 = (ii - 1 + Nx) % Nx; \
                double lapl = (phi[a][ip1] - 2.0*phi[a][ii] + phi[a][im1]) / dx2; \
                double fp = force_triple(phi[0][ii], phi[1][ii], phi[2][ii], a); \
                acc[a][ii] = lapl - m2*phi[a][ii] \
                           - lambda*(phi[b][ii] + phi[c_idx][ii]) + fp; \
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

            /* Position-dependent damping: strong between oscillons, weak inside cores */
            for (int a = 0; a < 3; a++)
                for (int i = 0; i < Nx; i++) {
                    double p2 = phi[0][i]*phi[0][i] + phi[1][i]*phi[1][i]
                              + phi[2][i]*phi[2][i];
                    double protect = p2 / (p2 + 0.01);
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

        double E_post = 0;
        for (int i = 0; i < Nx; i++)
            E_post += energy_density_p(phi[0], phi[1], phi[2],
                                        vel[0], vel[1], vel[2], i, Nx, dx) * dx;
        printf("  Pre-equilibration done. E_post=%.4f\n", E_post);

        find_positions_periodic(phi[0], phi[1], phi[2],
                                 vel[0], vel[1], vel[2],
                                 Nx, dx, L, x_osc, N_osc);
        printf("  Post-damping positions:\n");
        for (int n = 0; n < N_osc; n++)
            printf("    n=%d: x=%.3f\n", n, x_osc[n]);

        /* Apply phonon kick */
        printf("  Applying phonon kick (A=%.4f)...\n", A_kick);
        for (int i = 0; i < Nx; i++) {
            double x = i * dx;
            double kick = A_kick * sin(2.0 * M_PI * x / L);
            for (int a = 0; a < 3; a++)
                vel[a][i] += kick * phi[a][i];
        }

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

    /* Phase coherence: record phi_a(x_n, t) at each oscillon center */
    double **diag_phi0_center = malloc(N_osc * sizeof(double*));
    for (int n = 0; n < N_osc; n++)
        diag_phi0_center[n] = malloc(max_diag * sizeof(double));

    int diag_every = (int)(diag_dt / dt);
    if (diag_every < 1) diag_every = 1;
    int print_every = Nt / 40;
    if (print_every < 1) print_every = 1;

    int ref_set = 0;

    for (int step = 0; step <= Nt; step++) {
        double t = step * dt;

        if (step % diag_every == 0 && n_diag < max_diag) {
            if (step > 0)
                find_positions_periodic(phi[0], phi[1], phi[2],
                                         vel[0], vel[1], vel[2],
                                         Nx, dx, L, x_osc, N_osc);

            if (!ref_set && t >= t_discard) {
                for (int n = 0; n < N_osc; n++)
                    x_osc_ref[n] = x_osc[n];
                ref_set = 1;
                printf("  ** Reference positions set at t=%.1f **\n", t);
            }

            double Etot = 0;
            for (int i = 0; i < Nx; i++)
                Etot += energy_density_p(phi[0], phi[1], phi[2],
                                          vel[0], vel[1], vel[2], i, Nx, dx) * dx;

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

            /* Record phi_0 at each oscillon center for phase coherence */
            for (int n = 0; n < N_osc; n++) {
                double xc = x_osc[n];
                int idx = (int)(xc / dx);
                if (idx < 0) idx = 0;
                if (idx >= Nx) idx = Nx - 1;
                diag_phi0_center[n][n_diag] = phi[0][idx];
            }

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

        /* Velocity Verlet with periodic BC + pairwise coupling */
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

    /* Compute displacements u_n(t) */
    double **u_n = malloc(N_osc * sizeof(double*));
    for (int n = 0; n < N_osc; n++) {
        u_n[n] = calloc(n_analysis, sizeof(double));
        double prev_x = diag_x[n][j_start];
        u_n[n][0] = 0;
        for (int j = 1; j < n_analysis; j++) {
            int jj = j + j_start;
            double dx_step = diag_x[n][jj] - prev_x;
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

    /* Remove center-of-mass drift */
    for (int j = 0; j < n_analysis; j++) {
        double com = 0;
        for (int n = 0; n < N_osc; n++) com += u_n[n][j];
        com /= N_osc;
        for (int n = 0; n < N_osc; n++) u_n[n][j] -= com;
    }

    /* Normal mode decomposition */
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
    double omega_max = 1.0;
    if (omega_max > omega_nyquist) omega_max = omega_nyquist;

    double *omega_peak = calloc(N_osc, sizeof(double));
    double *power_peak = calloc(N_osc, sizeof(double));

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

    /* ===== Step 7: Phase coherence analysis ===== */
    printf("\n===== Step 7: Phase coherence (time crystal test) =====\n");

    /* For each oscillon, compute DFT of phi_0(x_n, t) to find breathing frequency.
     * Phase-locked = all oscillons have same frequency AND small phase spread. */
    double *osc_omega = calloc(N_osc, sizeof(double));
    double *osc_phase = calloc(N_osc, sizeof(double));
    double *osc_amp   = calloc(N_osc, sizeof(double));

    for (int n = 0; n < N_osc; n++) {
        /* DFT of phi0 at oscillon n center, analysis window only */
        double best_pow_n = 0, best_om_n = 0;
        double best_re = 0, best_im = 0;

        /* Scan around expected breathing frequency */
        double omega_lo = 0.3;
        double omega_hi = 2.0;
        int n_scan = 2000;

        for (int k = 0; k < n_scan; k++) {
            double omega = omega_lo + (omega_hi - omega_lo) * k / n_scan;
            double re = 0, im = 0;
            for (int j = j_start; j < n_diag; j++) {
                double tj = diag_t[j];
                double val = diag_phi0_center[n][j];
                re += val * cos(omega * tj) * diag_dt;
                im += val * sin(omega * tj) * diag_dt;
            }
            double pw = re * re + im * im;
            if (pw > best_pow_n) {
                best_pow_n = pw;
                best_om_n = omega;
                best_re = re;
                best_im = im;
            }
        }

        osc_omega[n] = best_om_n;
        osc_phase[n] = atan2(best_im, best_re);
        osc_amp[n] = sqrt(best_pow_n);
    }

    /* Compute frequency spread and phase spread */
    double omega_mean = 0;
    for (int n = 0; n < N_osc; n++) omega_mean += osc_omega[n];
    omega_mean /= N_osc;
    double omega_std = 0;
    for (int n = 0; n < N_osc; n++) {
        double d = osc_omega[n] - omega_mean;
        omega_std += d * d;
    }
    omega_std = sqrt(omega_std / N_osc);

    /* Phase spread: compute circular variance */
    double cos_sum = 0, sin_sum = 0;
    for (int n = 0; n < N_osc; n++) {
        cos_sum += cos(osc_phase[n]);
        sin_sum += sin(osc_phase[n]);
    }
    double R_phase = sqrt(cos_sum*cos_sum + sin_sum*sin_sum) / N_osc;
    /* R_phase = 1 means perfect phase coherence, 0 means random */
    double phase_mean = atan2(sin_sum, cos_sum);

    printf("  Breathing frequencies:\n");
    for (int n = 0; n < N_osc; n++)
        printf("    osc %d: omega=%.4f  phase=%.3f  amp=%.4e\n",
               n, osc_omega[n], osc_phase[n], osc_amp[n]);
    printf("  Mean omega = %.4f +/- %.4f (%.1f%% spread)\n",
           omega_mean, omega_std, 100.0 * omega_std / (omega_mean + 1e-30));
    printf("  Phase coherence R = %.4f (1=locked, 0=random)\n", R_phase);
    printf("  Mean phase = %.3f rad\n", phase_mean);

    int phase_locked = (omega_std / (omega_mean + 1e-30) < 0.01 && R_phase > 0.9);
    printf("  Phase-locked (time crystal)? %s\n", phase_locked ? "YES" : "NO");

    /* Write phase coherence data */
    {
        char path[600];
        snprintf(path, sizeof(path), "%s/phase_coherence.tsv", outdir);
        FILE *fp = fopen(path, "w");
        if (fp) {
            fprintf(fp, "osc\tomega\tphase\tamp\n");
            for (int n = 0; n < N_osc; n++)
                fprintf(fp, "%d\t%.6f\t%.6f\t%.8e\n",
                        n, osc_omega[n], osc_phase[n], osc_amp[n]);
            fclose(fp);
            printf("  Wrote %s\n", path);
        }
    }

    /* ===== Summary ===== */
    printf("\n===== SUMMARY =====\n");
    printf("  lambda = %.4f\n", lambda);
    printf("  Chain: %d oscillons, spacing d=%.1f, periodic L=%.1f\n",
           N_osc, d_space, L);
    printf("  Oscillon mass M = %.4f\n", eq_mass_osc);
    printf("  Breathing freq = %.4f\n", eq_omega_osc);
    printf("  Initial energy: %.4f, Final: %.4f (%.2f%%)\n",
           E_init, diag_E[n_diag-1],
           100.0 * (diag_E[n_diag-1] - E_init) / (fabs(E_init) + 1e-30));
    printf("  Mean spacing: initial=%.3f, final=%.3f\n",
           diag_mean_spacing[0], diag_mean_spacing[n_diag-1]);
    printf("  Max displacement at end: %.4f\n", diag_maxdisp[n_diag-1]);

    /* Max displacement over all time */
    double max_u_ever = 0;
    for (int j = j_start; j < n_diag; j++)
        if (diag_maxdisp[j] > max_u_ever) max_u_ever = diag_maxdisp[j];
    printf("  Max displacement ever (t>%.0f): %.4f\n", t_discard, max_u_ever);

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

        double max_u_final = diag_maxdisp[n_diag-1];
        printf("  Stability: max_u/d = %.3f", max_u_final / d_space);
        if (max_u_final < d_space / 4)
            printf(" (STABLE)\n");
        else if (max_u_final < d_space / 2)
            printf(" (MARGINAL)\n");
        else
            printf(" (UNSTABLE)\n");
    }

    printf("\n  Phase coherence: R=%.4f, omega_spread=%.4f/%.4f = %.2f%%\n",
           R_phase, omega_std, omega_mean,
           100.0 * omega_std / (omega_mean + 1e-30));
    printf("  TIME CRYSTAL? %s\n", phase_locked ? "YES" : "NO");

    /* Pairwise coupling diagnostics */
    if (lambda > 0) {
        double m2_anti = mass * mass - lambda;
        double m_A = m2_anti > 0 ? sqrt(m2_anti) : 0.0;
        double range_A = m_A > 0 ? 1.0 / m_A : 1e30;
        printf("\n  Pairwise coupling:\n");
        printf("    m_A (Proca mass) = %.4f\n", m_A);
        printf("    Proca range = %.4f\n", range_A);
        printf("    d/range = %.2f\n", d_space / range_A);
    }

cleanup:
    for (int a = 0; a < 3; a++) {
        free(phi[a]); free(vel[a]); free(acc[a]);
    }
    free(x_osc); free(x_osc_ref);
    for (int n = 0; n < N_osc; n++) { free(diag_x[n]); free(diag_phi0_center[n]); }
    free(diag_x); free(diag_phi0_center);
    free(diag_t); free(diag_E);
    free(diag_mean_spacing); free(diag_maxdisp);
    for (int q = 0; q < N_osc; q++) { free(Q_re[q]); free(Q_im[q]); }
    free(Q_re); free(Q_im);
    free(omega_peak); free(power_peak);
    for (int n = 0; n < N_osc; n++) free(u_n[n]);
    free(u_n);
    free(osc_omega); free(osc_phase); free(osc_amp);
    for (int a = 0; a < 3; a++) {
        if (eq_phi[a]) free(eq_phi[a]);
        if (eq_vel[a]) free(eq_vel[a]);
    }

    printf("\n=== Test E Complete ===\n");
    return 0;
}
