/*
 * combo25.c — Combo 2+5: Inertia + Condensate — Lattice Deformation
 *
 * Build a stable 8-oscillon periodic chain (lambda=0.5, from Test E).
 * After equilibration: apply a LOCALIZED antisymmetric perturbation at
 * oscillon #4: boost phi_1 by +eps, phi_2 by -eps (Gaussian envelope).
 * Evolve and track how the antisymmetric perturbation propagates.
 *
 * The antisymmetric mode breaks the phi_1=phi_2=phi_3 symmetry.
 * Symmetric phonon = compression wave (scalar, spin-0 analog).
 * Antisymmetric phonon = shear-like wave (candidate spin-2 in higher dim).
 *
 * Lagrangian:
 *   L = sum_a [ (1/2)(dt phi_a)^2 - (1/2)(dx phi_a)^2 - (m^2/2)phi_a^2 ]
 *     - (mu/2) P^2 / (1 + kappa P^2)         [triple product]
 *     - lambda(phi_1 phi_2 + phi_2 phi_3 + phi_3 phi_1)  [pairwise]
 *
 * Compile: gcc -O3 -Wall -o combo25 src/combo25.c -lm
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* --- Parameters --- */
static double mu       = -20.0;
static double kappa    = 20.0;
static double mass     = 1.0;
static double lambda   = 0.5;    /* pairwise coupling */
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
static double t_pre_equil = 500.0;
static double gamma_pre = 0.001;
static double t_settle = 5000.0;    /* conservative settle before perturbation */
static double t_evolve = 10000.0;   /* post-perturbation evolution */
static double diag_dt  = 5.0;
static int    pts_per_d = 320;

/* Perturbation */
static double eps_pert = 0.01;
static int    pert_osc = 4;        /* which oscillon to perturb (0-indexed: #4) */
static double pert_sigma = 3.0;    /* Gaussian envelope width */

/* Saved equilibrated profile */
static double *eq_phi[3], *eq_vel[3];
static int    eq_Nx = 0;
static double eq_dx = 0;
static double eq_xmax = 0;
static double eq_mass_osc = 0;
static double eq_omega_osc = 0;

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

/* Simple LCG random in [0,1) */
static double rand01(unsigned int *state)
{
    *state = *state * 1103515245u + 12345u;
    return (double)((*state >> 16) & 0x7fff) / 32768.0;
}

/* Energy density (periodic) */
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
    e += lambda * (phi0[i]*phi1[i] + phi1[i]*phi2[i] + phi2[i]*phi0[i]);
    double P = phi0[i] * phi1[i] * phi2[i];
    double P2 = P * P;
    e += 0.5 * mu * P2 / (1.0 + kappa * P2);
    return e;
}

/* Energy density (non-periodic, for equilibration) */
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
    e += lambda * (phi0[i]*phi1[i] + phi1[i]*phi2[i] + phi2[i]*phi0[i]);
    double P = phi0[i] * phi1[i] * phi2[i];
    double P2 = P * P;
    e += 0.5 * mu * P2 / (1.0 + kappa * P2);
    return e;
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

    double m2_eff = m2 + 2.0 * lambda;
    if (m2_eff < m2) m2_eff = m2;
    double kmax = M_PI / dx;
    double dt   = 0.8 * 2.0 / sqrt(kmax * kmax + m2_eff);
    int    Nt   = (int)(t_equil / dt) + 1;

    printf("  Nx=%d xmax=%.1f dx=%.5f dt=%.6f Nt=%d tfinal=%.0f\n",
           Nx, xmax, dx, dt, Nt, t_equil);

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

    /* Initialize: symmetric Gaussians */
    for (int a = 0; a < 3; a++)
        for (int i = 0; i < Nx; i++) {
            double x = -xmax + i * dx;
            phi[a][i] = A_init * exp(-x * x / (2.0 * sig * sig));
        }

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
    int print_every = Nt / 10;
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

        if (n % dft_every == 0 && n_dft < max_dft) {
            phi0_hist[n_dft] = phi[0][ic];
            t_hist[n_dft] = t;
            n_dft++;
        }

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

    /* Breathing frequency from DFT */
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
                       int Nx, double dx, double L)
{
    for (int a = 0; a < 3; a++) {
        memset(phi[a], 0, Nx * sizeof(double));
        memset(vel[a], 0, Nx * sizeof(double));
    }

    for (int n = 0; n < N_osc; n++) {
        double xc = n * d_space + d_space / 2.0;
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
static void find_positions(double *phi0, double *phi1, double *phi2,
                           double *v0, double *v1, double *v2,
                           int Nx, double dx, double L,
                           double *x_osc)
{
    for (int n = 0; n < N_osc; n++) {
        double x_left = x_osc[(n - 1 + N_osc) % N_osc];
        double x_right = x_osc[(n + 1) % N_osc];
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
            double new_x = xc + num / den;
            while (new_x < 0) new_x += L;
            while (new_x >= L) new_x -= L;
            x_osc[n] = new_x;
        }
    }
}

/* ===================================================================
 *  Measure antisymmetric content at each oscillon
 *  S_n = (phi1 + phi2 + phi3) / 3   (symmetric)
 *  A12_n = phi1 - phi2               (antisymmetric 12)
 *  A13_n = phi1 - phi3               (antisymmetric 13)
 *  Measured as amplitude-weighted integral near oscillon center
 * =================================================================== */
static void measure_asym(double *phi[3], int Nx, double dx, double L,
                         double *x_osc,
                         double *S_amp, double *A12_amp, double *A13_amp)
{
    double core_r = 2.0 * sig;  /* integrate within 2*sigma of center */

    for (int n = 0; n < N_osc; n++) {
        double xc = x_osc[n];
        double s_num = 0, a12_num = 0, a13_num = 0;
        double w_sum = 0;

        for (int i = 0; i < Nx; i++) {
            double x = i * dx;
            double rel = x - xc;
            while (rel > L / 2) rel -= L;
            while (rel < -L / 2) rel += L;

            if (fabs(rel) < core_r) {
                double w = exp(-rel * rel / (2.0 * sig * sig));
                double p1 = phi[0][i], p2 = phi[1][i], p3 = phi[2][i];
                s_num   += w * (p1 + p2 + p3) / 3.0;
                a12_num += w * (p1 - p2);
                a13_num += w * (p1 - p3);
                w_sum   += w;
            }
        }

        if (w_sum > 1e-30) {
            S_amp[n]   = s_num / w_sum;
            A12_amp[n] = a12_num / w_sum;
            A13_amp[n] = a13_num / w_sum;
        } else {
            S_amp[n] = A12_amp[n] = A13_amp[n] = 0;
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
        else if (!strcmp(argv[i], "-eps"))      eps_pert = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-pert_osc")) pert_osc = atoi(argv[i+1]);
        else if (!strcmp(argv[i], "-t_settle")) t_settle = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-t_evolve")) t_evolve = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-diag_dt"))  diag_dt  = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-pts_per_d")) pts_per_d = atoi(argv[i+1]);
        else if (!strcmp(argv[i], "-t_pre"))    t_pre_equil = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-gamma"))    gamma_pre = atof(argv[i+1]);
        else { fprintf(stderr, "Unknown flag: %s\n", argv[i]); exit(1); }
    }

    printf("=== Combo 2+5: Antisymmetric Lattice Phonon (Shear Mode) ===\n");
    printf("  lambda=%.4f  mu=%.1f  kappa=%.1f  mass=%.4f\n", lambda, mu, kappa, mass);
    printf("  N_osc=%d  d=%.1f  eps=%.6f  pert_osc=%d\n", N_osc, d_space, eps_pert, pert_osc);
    printf("  t_settle=%.0f  t_evolve=%.0f\n", t_settle, t_evolve);

    /* ===== Step 1: Equilibrate single oscillon ===== */
    run_equilibrate();

    /* ===== Step 2: Build periodic chain ===== */
    printf("===== Step 2: Build %d-oscillon periodic chain =====\n", N_osc);

    double L = N_osc * d_space;
    int Nx = N_osc * pts_per_d;
    double dx = L / Nx;
    double dx2 = dx * dx;
    double m2 = mass * mass;

    double m2_eff = m2 + 2.0 * lambda;
    if (m2_eff < m2) m2_eff = m2;
    double kmax_c = M_PI / dx;
    double dt = 0.8 * 2.0 / sqrt(kmax_c * kmax_c + m2_eff);

    printf("  L=%.1f Nx=%d dx=%.5f dt=%.6f\n", L, Nx, dx, dt);

    double *x_osc = malloc(N_osc * sizeof(double));
    for (int n = 0; n < N_osc; n++)
        x_osc[n] = n * d_space + d_space / 2.0;

    double *phi[3], *vel[3], *acc[3];
    for (int a = 0; a < 3; a++) {
        phi[a] = calloc(Nx, sizeof(double));
        vel[a] = calloc(Nx, sizeof(double));
        acc[a] = calloc(Nx, sizeof(double));
    }

    init_chain(phi, vel, Nx, dx, L);

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

    /* ===== Pre-equilibrate with damping ===== */
    {
        int Nt_pre = (int)(t_pre_equil / dt) + 1;
        int print_pre = Nt_pre / 5;
        if (print_pre < 1) print_pre = 1;

        printf("\n  Pre-equilibrating with damping (gamma=%.4f, t=%.0f)...\n",
               gamma_pre, t_pre_equil);

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

            /* Position-dependent damping: strong between oscillons, weak in cores */
            for (int a = 0; a < 3; a++)
                for (int i = 0; i < Nx; i++) {
                    double p2 = phi[0][i]*phi[0][i] + phi[1][i]*phi[1][i]
                              + phi[2][i]*phi[2][i];
                    double protect = p2 / (p2 + 0.01);
                    double local_damp = 1.0 - gamma_pre * dt * (1.0 - protect);
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
        COMPUTE_ACC_P();
    }

    /* ===== Conservative settle (no perturbation yet) ===== */
    printf("\n===== Step 3: Conservative settle for t=%.0f =====\n", t_settle);
    {
        int Nt_settle = (int)(t_settle / dt) + 1;
        int print_every = Nt_settle / 10;
        if (print_every < 1) print_every = 1;

        for (int step = 0; step <= Nt_settle; step++) {
            if (step % print_every == 0) {
                double Etot = 0;
                for (int i = 0; i < Nx; i++)
                    Etot += energy_density_p(phi[0], phi[1], phi[2],
                                              vel[0], vel[1], vel[2], i, Nx, dx) * dx;

                /* Check max antisymmetric content */
                find_positions(phi[0], phi[1], phi[2],
                               vel[0], vel[1], vel[2], Nx, dx, L, x_osc);
                double S_amp[8], A12_amp[8], A13_amp[8];
                measure_asym(phi, Nx, dx, L, x_osc, S_amp, A12_amp, A13_amp);
                double max_asym = 0;
                for (int n = 0; n < N_osc; n++) {
                    double a = fabs(A12_amp[n]);
                    if (a > max_asym) max_asym = a;
                }
                printf("  settle t=%7.1f  E=%.4f  max|A12|=%.6e\n",
                       step * dt, Etot, max_asym);
            }

            if (step == Nt_settle) break;

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
    }

    /* ===== Step 4: Apply antisymmetric perturbation ===== */
    printf("\n===== Step 4: Apply antisymmetric perturbation at oscillon #%d =====\n", pert_osc);
    find_positions(phi[0], phi[1], phi[2],
                   vel[0], vel[1], vel[2], Nx, dx, L, x_osc);

    /* Measure pre-perturbation baseline */
    double S_pre[8], A12_pre[8], A13_pre[8];
    measure_asym(phi, Nx, dx, L, x_osc, S_pre, A12_pre, A13_pre);
    printf("  Pre-perturbation A12 at each oscillon:\n");
    for (int n = 0; n < N_osc; n++)
        printf("    osc %d: S=%.6f  A12=%.6e  A13=%.6e\n",
               n, S_pre[n], A12_pre[n], A13_pre[n]);

    double E_before = 0;
    for (int i = 0; i < Nx; i++)
        E_before += energy_density_p(phi[0], phi[1], phi[2],
                                      vel[0], vel[1], vel[2], i, Nx, dx) * dx;
    printf("  E_before = %.6f\n", E_before);

    /* Apply: phi_1 += +eps * Gaussian,  phi_2 += -eps * Gaussian at osc #pert_osc */
    {
        double xc = x_osc[pert_osc];
        for (int i = 0; i < Nx; i++) {
            double x = i * dx;
            double rel = x - xc;
            while (rel > L / 2) rel -= L;
            while (rel < -L / 2) rel += L;
            double g = eps_pert * exp(-rel * rel / (2.0 * pert_sigma * pert_sigma));
            phi[0][i] += g;    /* phi_1 += eps */
            phi[1][i] -= g;    /* phi_2 -= eps */
        }
    }

    COMPUTE_ACC_P();

    double E_after = 0;
    for (int i = 0; i < Nx; i++)
        E_after += energy_density_p(phi[0], phi[1], phi[2],
                                     vel[0], vel[1], vel[2], i, Nx, dx) * dx;
    printf("  E_after = %.6f  (Delta E = %.6e)\n", E_after, E_after - E_before);

    /* Measure post-perturbation */
    double S_post[8], A12_post[8], A13_post[8];
    measure_asym(phi, Nx, dx, L, x_osc, S_post, A12_post, A13_post);
    printf("  Post-perturbation A12:\n");
    for (int n = 0; n < N_osc; n++)
        printf("    osc %d: A12=%.6e (delta=%.6e)\n",
               n, A12_post[n], A12_post[n] - A12_pre[n]);

    /* ===== Step 5: Post-perturbation evolution ===== */
    printf("\n===== Step 5: Post-perturbation evolution for t=%.0f =====\n", t_evolve);

    int Nt_evolve = (int)(t_evolve / dt) + 1;
    int diag_every = (int)(diag_dt / dt);
    if (diag_every < 1) diag_every = 1;
    int max_diag = Nt_evolve / diag_every + 10;
    int print_every = Nt_evolve / 40;
    if (print_every < 1) print_every = 1;

    /* Diagnostic arrays */
    double *diag_t      = malloc(max_diag * sizeof(double));
    double *diag_E      = malloc(max_diag * sizeof(double));
    double **diag_S     = malloc(N_osc * sizeof(double*));
    double **diag_A12   = malloc(N_osc * sizeof(double*));
    double **diag_A13   = malloc(N_osc * sizeof(double*));
    double **diag_xpos  = malloc(N_osc * sizeof(double*));
    for (int n = 0; n < N_osc; n++) {
        diag_S[n]    = malloc(max_diag * sizeof(double));
        diag_A12[n]  = malloc(max_diag * sizeof(double));
        diag_A13[n]  = malloc(max_diag * sizeof(double));
        diag_xpos[n] = malloc(max_diag * sizeof(double));
    }
    int n_diag = 0;

    for (int step = 0; step <= Nt_evolve; step++) {
        double t = step * dt;

        if (step % diag_every == 0 && n_diag < max_diag) {
            if (step > 0)
                find_positions(phi[0], phi[1], phi[2],
                               vel[0], vel[1], vel[2], Nx, dx, L, x_osc);

            double Etot = 0;
            for (int i = 0; i < Nx; i++)
                Etot += energy_density_p(phi[0], phi[1], phi[2],
                                          vel[0], vel[1], vel[2], i, Nx, dx) * dx;

            double S_amp[8], A12_amp[8], A13_amp[8];
            measure_asym(phi, Nx, dx, L, x_osc, S_amp, A12_amp, A13_amp);

            diag_t[n_diag] = t;
            diag_E[n_diag] = Etot;
            for (int n = 0; n < N_osc; n++) {
                diag_S[n][n_diag]    = S_amp[n];
                diag_A12[n][n_diag]  = A12_amp[n];
                diag_A13[n][n_diag]  = A13_amp[n];
                diag_xpos[n][n_diag] = x_osc[n];
            }

            if (step % print_every == 0) {
                double max_a12 = 0;
                int max_n = 0;
                for (int n = 0; n < N_osc; n++) {
                    if (fabs(A12_amp[n]) > max_a12) {
                        max_a12 = fabs(A12_amp[n]);
                        max_n = n;
                    }
                }
                printf("  t=%7.1f  E=%.4f  max|A12|=%.6e at osc %d  A12[%d]=%.6e\n",
                       t, Etot, max_a12, max_n, pert_osc, A12_amp[pert_osc]);
            }

            n_diag++;
        }

        if (step == Nt_evolve) break;

        /* Velocity Verlet */
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

    /* ===== Step 6: Output data ===== */
    printf("\n===== Step 6: Writing output files =====\n");

    /* combo25_asym_ts.tsv — antisymmetric content vs time at each oscillon */
    {
        char path[600];
        snprintf(path, sizeof(path), "%s/combo25_asym_ts.tsv", outdir);
        FILE *fp = fopen(path, "w");
        if (fp) {
            fprintf(fp, "time");
            for (int n = 0; n < N_osc; n++) fprintf(fp, "\tA12_%d", n);
            for (int n = 0; n < N_osc; n++) fprintf(fp, "\tA13_%d", n);
            for (int n = 0; n < N_osc; n++) fprintf(fp, "\tS_%d", n);
            fprintf(fp, "\n");
            for (int j = 0; j < n_diag; j++) {
                fprintf(fp, "%.4f", diag_t[j]);
                for (int n = 0; n < N_osc; n++)
                    fprintf(fp, "\t%.8e", diag_A12[n][j]);
                for (int n = 0; n < N_osc; n++)
                    fprintf(fp, "\t%.8e", diag_A13[n][j]);
                for (int n = 0; n < N_osc; n++)
                    fprintf(fp, "\t%.8e", diag_S[n][j]);
                fprintf(fp, "\n");
            }
            fclose(fp);
            printf("  Wrote %s (%d rows)\n", path, n_diag);
        }
    }

    /* combo25_energy_ts.tsv */
    {
        char path[600];
        snprintf(path, sizeof(path), "%s/combo25_energy_ts.tsv", outdir);
        FILE *fp = fopen(path, "w");
        if (fp) {
            fprintf(fp, "time\tE_total\n");
            for (int j = 0; j < n_diag; j++)
                fprintf(fp, "%.4f\t%.8e\n", diag_t[j], diag_E[j]);
            fclose(fp);
            printf("  Wrote %s\n", path);
        }
    }

    /* combo25_positions.tsv */
    {
        char path[600];
        snprintf(path, sizeof(path), "%s/combo25_positions.tsv", outdir);
        FILE *fp = fopen(path, "w");
        if (fp) {
            fprintf(fp, "time");
            for (int n = 0; n < N_osc; n++) fprintf(fp, "\tX_%d", n);
            fprintf(fp, "\n");
            for (int j = 0; j < n_diag; j++) {
                fprintf(fp, "%.4f", diag_t[j]);
                for (int n = 0; n < N_osc; n++)
                    fprintf(fp, "\t%.6f", diag_xpos[n][j]);
                fprintf(fp, "\n");
            }
            fclose(fp);
            printf("  Wrote %s\n", path);
        }
    }

    /* ===== Step 7: Normal mode decomposition (symmetric + antisymmetric) ===== */
    printf("\n===== Step 7: Normal mode decomposition =====\n");

    /* For the antisymmetric mode A12_n(t), decompose into spatial Fourier modes
     * Q_q(t) = (1/sqrt(N)) sum_n A12_n(t) * exp(-2pi i q n / N)
     * Then DFT in time to get omega for each q => dispersion relation.
     *
     * Similarly for symmetric displacement u_n(t). */

    double T_analysis = diag_t[n_diag-1] - diag_t[0];
    double dt_diag_eff = T_analysis / (n_diag - 1);

    /* --- Antisymmetric phonon spectrum --- */
    printf("\n  --- Antisymmetric (shear) mode spectrum ---\n");

    /* Spatial Fourier of A12 */
    double **AQ_re = malloc(N_osc * sizeof(double*));
    double **AQ_im = malloc(N_osc * sizeof(double*));
    for (int q = 0; q < N_osc; q++) {
        AQ_re[q] = calloc(n_diag, sizeof(double));
        AQ_im[q] = calloc(n_diag, sizeof(double));
    }

    for (int j = 0; j < n_diag; j++) {
        for (int q = 0; q < N_osc; q++) {
            double re = 0, im = 0;
            for (int n = 0; n < N_osc; n++) {
                double phase = -2.0 * M_PI * q * n / N_osc;
                re += diag_A12[n][j] * cos(phase);
                im += diag_A12[n][j] * sin(phase);
            }
            AQ_re[q][j] = re / sqrt((double)N_osc);
            AQ_im[q][j] = im / sqrt((double)N_osc);
        }
    }

    /* Power spectrum in time for each spatial mode */
    int n_omega = n_diag / 2;
    if (n_omega > 4000) n_omega = 4000;
    double omega_nyquist = M_PI / dt_diag_eff;
    double omega_max = 2.0;
    if (omega_max > omega_nyquist) omega_max = omega_nyquist;

    double *asym_omega_peak = calloc(N_osc, sizeof(double));
    double *asym_power_peak = calloc(N_osc, sizeof(double));

    char aspec_path[600];
    snprintf(aspec_path, sizeof(aspec_path), "%s/combo25_asym_spectrum_2d.tsv", outdir);
    FILE *faspec = fopen(aspec_path, "w");
    if (faspec) fprintf(faspec, "q\tk\tomega\tpower\n");

    for (int q = 0; q <= N_osc / 2; q++) {
        double k_q = 2.0 * M_PI * q / (N_osc * d_space);
        double best_pow = 0, best_omega = 0;

        for (int l = 1; l <= n_omega; l++) {
            double omega = omega_max * l / n_omega;

            double dft_re = 0, dft_im = 0;
            for (int j = 0; j < n_diag; j++) {
                double t_j = j * dt_diag_eff;
                double c = cos(omega * t_j);
                double s = sin(omega * t_j);
                dft_re += (AQ_re[q][j] * c + AQ_im[q][j] * s) * dt_diag_eff;
                dft_im += (-AQ_re[q][j] * s + AQ_im[q][j] * c) * dt_diag_eff;
            }
            double power = (dft_re * dft_re + dft_im * dft_im) / (T_analysis * T_analysis);

            if (faspec) fprintf(faspec, "%d\t%.6f\t%.6f\t%.8e\n", q, k_q, omega, power);

            if (power > best_pow) {
                best_pow = power;
                best_omega = omega;
            }
        }

        asym_omega_peak[q] = best_omega;
        asym_power_peak[q] = best_pow;
        double cs = (k_q > 1e-10) ? best_omega / k_q : 0;
        printf("  Asym mode q=%d: k=%.5f  omega=%.5f  power=%.4e  c_s=%.4f\n",
               q, k_q, best_omega, best_pow, cs);
    }
    if (faspec) { fclose(faspec); printf("  Wrote %s\n", aspec_path); }

    /* --- Symmetric (positional) phonon spectrum --- */
    printf("\n  --- Symmetric (compression) mode spectrum ---\n");

    /* Compute displacements u_n(t) from positions */
    double **u_n = malloc(N_osc * sizeof(double*));
    for (int n = 0; n < N_osc; n++) {
        u_n[n] = calloc(n_diag, sizeof(double));
        double prev_x = diag_xpos[n][0];
        u_n[n][0] = 0;
        for (int j = 1; j < n_diag; j++) {
            double dx_step = diag_xpos[n][j] - prev_x;
            while (dx_step > L / 2) dx_step -= L;
            while (dx_step < -L / 2) dx_step += L;
            u_n[n][j] = u_n[n][j - 1] + dx_step;
            prev_x = diag_xpos[n][j];
        }
    }

    /* Remove COM drift */
    for (int j = 0; j < n_diag; j++) {
        double com = 0;
        for (int n = 0; n < N_osc; n++) com += u_n[n][j];
        com /= N_osc;
        for (int n = 0; n < N_osc; n++) u_n[n][j] -= com;
    }

    /* Spatial Fourier of u_n */
    double **SQ_re = malloc(N_osc * sizeof(double*));
    double **SQ_im = malloc(N_osc * sizeof(double*));
    for (int q = 0; q < N_osc; q++) {
        SQ_re[q] = calloc(n_diag, sizeof(double));
        SQ_im[q] = calloc(n_diag, sizeof(double));
    }

    for (int j = 0; j < n_diag; j++) {
        for (int q = 0; q < N_osc; q++) {
            double re = 0, im = 0;
            for (int n = 0; n < N_osc; n++) {
                double phase = -2.0 * M_PI * q * n / N_osc;
                re += u_n[n][j] * cos(phase);
                im += u_n[n][j] * sin(phase);
            }
            SQ_re[q][j] = re / sqrt((double)N_osc);
            SQ_im[q][j] = im / sqrt((double)N_osc);
        }
    }

    double *sym_omega_peak = calloc(N_osc, sizeof(double));
    double *sym_power_peak = calloc(N_osc, sizeof(double));

    char sspec_path[600];
    snprintf(sspec_path, sizeof(sspec_path), "%s/combo25_sym_spectrum_2d.tsv", outdir);
    FILE *fsspec = fopen(sspec_path, "w");
    if (fsspec) fprintf(fsspec, "q\tk\tomega\tpower\n");

    for (int q = 0; q <= N_osc / 2; q++) {
        double k_q = 2.0 * M_PI * q / (N_osc * d_space);
        double best_pow = 0, best_omega = 0;

        for (int l = 1; l <= n_omega; l++) {
            double omega = omega_max * l / n_omega;
            double dft_re = 0, dft_im = 0;
            for (int j = 0; j < n_diag; j++) {
                double t_j = j * dt_diag_eff;
                double c = cos(omega * t_j);
                double s = sin(omega * t_j);
                dft_re += (SQ_re[q][j] * c + SQ_im[q][j] * s) * dt_diag_eff;
                dft_im += (-SQ_re[q][j] * s + SQ_im[q][j] * c) * dt_diag_eff;
            }
            double power = (dft_re * dft_re + dft_im * dft_im) / (T_analysis * T_analysis);

            if (fsspec) fprintf(fsspec, "%d\t%.6f\t%.6f\t%.8e\n", q, k_q, omega, power);

            if (power > best_pow) {
                best_pow = power;
                best_omega = omega;
            }
        }

        sym_omega_peak[q] = best_omega;
        sym_power_peak[q] = best_pow;
        double cs = (k_q > 1e-10) ? best_omega / k_q : 0;
        printf("  Sym  mode q=%d: k=%.5f  omega=%.5f  power=%.4e  c_s=%.4f\n",
               q, k_q, best_omega, best_pow, cs);
    }
    if (fsspec) { fclose(fsspec); printf("  Wrote %s\n", sspec_path); }

    /* ===== Step 8: Propagation speed measurement ===== */
    printf("\n===== Step 8: Propagation analysis =====\n");

    /* For each oscillon, find the TIME at which |A12| first exceeds a threshold.
     * This gives the arrival time, and distance/time = speed. */
    double threshold = 0.1 * fabs(diag_A12[pert_osc][1]);  /* 10% of initial perturbation */
    printf("  Detection threshold: %.6e (10%% of initial A12 at osc %d)\n",
           threshold, pert_osc);

    double *arrival_t = malloc(N_osc * sizeof(double));
    for (int n = 0; n < N_osc; n++) arrival_t[n] = -1;

    for (int n = 0; n < N_osc; n++) {
        /* Skip the perturbed oscillon itself */
        if (n == pert_osc) { arrival_t[n] = 0; continue; }

        /* Subtract baseline: average of first few samples */
        double baseline = 0;
        int n_base = n_diag < 20 ? n_diag : 20;
        for (int j = 0; j < n_base; j++) baseline += diag_A12[n][j];
        baseline /= n_base;

        for (int j = 1; j < n_diag; j++) {
            if (fabs(diag_A12[n][j] - baseline) > threshold) {
                arrival_t[n] = diag_t[j];
                break;
            }
        }
    }

    printf("  Arrival times:\n");
    for (int n = 0; n < N_osc; n++) {
        int dist = abs(n - pert_osc);
        if (dist > N_osc / 2) dist = N_osc - dist;
        double d_travel = dist * d_space;
        double speed = (arrival_t[n] > 0 && dist > 0) ? d_travel / arrival_t[n] : 0;
        printf("    osc %d: dist=%d (%.1f)  t_arr=%.1f  c_shear=%.4f\n",
               n, dist, d_travel, arrival_t[n], speed);
    }

    /* ===== Step 9: Amplitude decay measurement ===== */
    printf("\n===== Step 9: Amplitude decay =====\n");

    /* Peak |A12| at each oscillon over all time */
    double *peak_A12 = calloc(N_osc, sizeof(double));
    for (int n = 0; n < N_osc; n++) {
        for (int j = 0; j < n_diag; j++) {
            double a = fabs(diag_A12[n][j]);
            if (a > peak_A12[n]) peak_A12[n] = a;
        }
    }

    printf("  Peak |A12| at each oscillon:\n");
    for (int n = 0; n < N_osc; n++) {
        int dist = abs(n - pert_osc);
        if (dist > N_osc / 2) dist = N_osc - dist;
        printf("    osc %d: dist=%d  peak|A12|=%.6e  ratio=%.4f\n",
               n, dist, peak_A12[n],
               peak_A12[pert_osc] > 0 ? peak_A12[n] / peak_A12[pert_osc] : 0);
    }

    /* Measure damping at source oscillon */
    double A12_initial = fabs(diag_A12[pert_osc][1]);
    double A12_final = 0;
    /* Average over last 10% of data */
    {
        int j_start = (int)(0.9 * n_diag);
        double sum = 0;
        int count = 0;
        for (int j = j_start; j < n_diag; j++) {
            sum += fabs(diag_A12[pert_osc][j]);
            count++;
        }
        if (count > 0) A12_final = sum / count;
    }
    double lifetime_ratio = (A12_initial > 0) ? A12_final / A12_initial : 0;
    printf("\n  Source oscillon (#%d) A12 decay:\n", pert_osc);
    printf("    Initial |A12| = %.6e\n", A12_initial);
    printf("    Final   |A12| = %.6e (last 10%%)\n", A12_final);
    printf("    Ratio final/initial = %.4f\n", lifetime_ratio);

    /* ===== Step 10: Summary ===== */
    printf("\n===== SUMMARY =====\n");
    printf("  Chain: %d oscillons, d=%.1f, lambda=%.4f, periodic L=%.1f\n",
           N_osc, d_space, lambda, L);
    printf("  Oscillon mass = %.4f, breathing freq = %.4f\n", eq_mass_osc, eq_omega_osc);
    printf("  Perturbation: eps=%.6f at oscillon #%d (phi1 += eps, phi2 -= eps)\n",
           eps_pert, pert_osc);
    printf("\n  Energy: initial=%.4f, final=%.4f (conservation: %.2e)\n",
           diag_E[0], diag_E[n_diag-1],
           (diag_E[n_diag-1] - diag_E[0]) / (fabs(diag_E[0]) + 1e-30));

    printf("\n  Dispersion comparison (symmetric vs antisymmetric):\n");
    printf("  %-4s  %-10s  %-12s  %-12s  %-12s  %-12s\n",
           "q", "k", "omega_sym", "c_s_sym", "omega_asym", "c_s_asym");
    printf("  ----  ----------  ------------  ------------  ------------  ------------\n");
    for (int q = 0; q <= N_osc / 2; q++) {
        double k_q = 2.0 * M_PI * q / (N_osc * d_space);
        double cs_s = (k_q > 1e-10) ? sym_omega_peak[q] / k_q : 0;
        double cs_a = (k_q > 1e-10) ? asym_omega_peak[q] / k_q : 0;
        printf("  %-4d  %-10.5f  %-12.5f  %-12.4f  %-12.5f  %-12.4f\n",
               q, k_q, sym_omega_peak[q], cs_s, asym_omega_peak[q], cs_a);
    }

    /* Overall speed estimates */
    if (N_osc >= 4) {
        double k1 = 2.0 * M_PI / (N_osc * d_space);
        double cs_sym = (k1 > 0) ? sym_omega_peak[1] / k1 : 0;
        double cs_asym = (k1 > 0) ? asym_omega_peak[1] / k1 : 0;
        printf("\n  Sound speeds (q=1 long-wavelength):\n");
        printf("    Symmetric (compression): c_s = %.4f\n", cs_sym);
        printf("    Antisymmetric (shear):   c_a = %.4f\n", cs_asym);
        if (cs_sym > 0)
            printf("    Ratio c_a/c_s = %.4f\n", cs_asym / cs_sym);
    }

    /* Coherence assessment */
    int propagates = 0;
    for (int n = 0; n < N_osc; n++) {
        if (n == pert_osc) continue;
        int dist = abs(n - pert_osc);
        if (dist > N_osc / 2) dist = N_osc - dist;
        if (dist >= 2 && peak_A12[n] > 0.05 * peak_A12[pert_osc])
            propagates = 1;
    }
    printf("\n  Propagating antisymmetric phonon (dist >= 2)? %s\n",
           propagates ? "YES" : "NO");

    int coherent = propagates && (lifetime_ratio > 0.1);
    printf("  Coherent shear mode? %s\n", coherent ? "YES" : "NO");
    printf("  (Criterion: propagates to dist>=2 at >5%% amplitude AND source retains >10%%)\n");

    /* Write summary modes file */
    {
        char path[600];
        snprintf(path, sizeof(path), "%s/combo25_modes.tsv", outdir);
        FILE *fp = fopen(path, "w");
        if (fp) {
            fprintf(fp, "q\tk\tomega_sym\tpower_sym\tomega_asym\tpower_asym\n");
            for (int q = 0; q <= N_osc / 2; q++) {
                double k_q = 2.0 * M_PI * q / (N_osc * d_space);
                fprintf(fp, "%d\t%.6f\t%.6f\t%.8e\t%.6f\t%.8e\n",
                        q, k_q, sym_omega_peak[q], sym_power_peak[q],
                        asym_omega_peak[q], asym_power_peak[q]);
            }
            fclose(fp);
            printf("  Wrote %s\n", path);
        }
    }

    /* Write propagation data */
    {
        char path[600];
        snprintf(path, sizeof(path), "%s/combo25_propagation.tsv", outdir);
        FILE *fp = fopen(path, "w");
        if (fp) {
            fprintf(fp, "osc\tdist\tpeak_A12\tarrival_t\tspeed\n");
            for (int n = 0; n < N_osc; n++) {
                int dist = abs(n - pert_osc);
                if (dist > N_osc / 2) dist = N_osc - dist;
                double d_travel = dist * d_space;
                double speed = (arrival_t[n] > 0 && dist > 0) ? d_travel / arrival_t[n] : 0;
                fprintf(fp, "%d\t%d\t%.8e\t%.4f\t%.6f\n",
                        n, dist, peak_A12[n], arrival_t[n], speed);
            }
            fclose(fp);
            printf("  Wrote %s\n", path);
        }
    }

    printf("\n=== Combo 2+5 Complete ===\n");

    /* Cleanup */
    for (int a = 0; a < 3; a++) { free(phi[a]); free(vel[a]); free(acc[a]); }
    free(x_osc);
    free(diag_t); free(diag_E);
    for (int n = 0; n < N_osc; n++) {
        free(diag_S[n]); free(diag_A12[n]); free(diag_A13[n]); free(diag_xpos[n]);
        free(u_n[n]);
    }
    free(diag_S); free(diag_A12); free(diag_A13); free(diag_xpos);
    free(u_n);
    for (int q = 0; q < N_osc; q++) {
        free(AQ_re[q]); free(AQ_im[q]);
        free(SQ_re[q]); free(SQ_im[q]);
    }
    free(AQ_re); free(AQ_im);
    free(SQ_re); free(SQ_im);
    free(asym_omega_peak); free(asym_power_peak);
    free(sym_omega_peak); free(sym_power_peak);
    free(arrival_t); free(peak_A12);
    for (int a = 0; a < 3; a++) {
        if (eq_phi[a]) free(eq_phi[a]);
        if (eq_vel[a]) free(eq_vel[a]);
    }

    return 0;
}
