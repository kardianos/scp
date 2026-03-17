/*
 * combo1257.c — Combo 1+2+5+7: EOS + Inertia + Condensate + Self-Ref = Full GR Analog
 *
 * 8-oscillon periodic chain with pairwise coupling lambda=0.5,
 * PLUS self-consistent gravitational potential Phi from Poisson equation.
 *
 * At each timestep:
 *   1. Compute energy density rho(x)
 *   2. Solve Poisson: d^2 Phi/dx^2 = alpha * rho(x)  (Thomas algorithm, periodic BC)
 *   3. Modify dynamics: m^2_eff = m^2*(1+2*Phi), Laplacian coeff = (1+4*Phi)
 *   4. Evolve fields with modified dynamics
 *
 * Phases:
 *   1. Equilibrate self-gravitating lattice (t=5000)
 *   2. Add random position perturbations, phonon spectrum (t=5000)
 *   3. Perturb one oscillon, track Phi(x,t) propagation
 *   4. Perturb two oscillons differently, equivalence principle test
 *
 * Compile: gcc -O3 -Wall -o combo1257 src/combo1257.c -lm
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
static double alpha    = -0.0001; /* gravitational coupling (negative = attractive) */
static double A_init   = 1.0;
static double sig      = 3.0;
static char   outdir[512] = "data";

/* Phase 1: single oscillon equilibration */
static int    Nx_eq    = 4000;
static double xmax_eq  = 100.0;
static double t_equil  = 5000.0;

/* Chain parameters */
static int    N_osc    = 8;
static double d_space  = 16.0;
static int    pts_per_d = 320;
static double t_pre_equil = 500.0;
static double gamma_pre = 0.001;
static double delta_max = 0.5;
static unsigned int rng_seed = 42;

/* Phase timing */
static double t_phase1  = 5000.0;   /* self-gravitating equilibration */
static double t_phase2  = 5000.0;   /* phonon spectrum measurement */
static double t_phase3  = 3000.0;   /* Phi propagation test */
static double t_phase4  = 3000.0;   /* equivalence principle test */
static double diag_dt   = 5.0;

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
 *  Solve Poisson equation: d^2 Phi/dx^2 = alpha * rho
 *  Periodic BC: use Sherman-Morrison on tridiagonal
 *  The periodic tridiagonal system: A*Phi = b where
 *    A_ii = -2/dx^2, A_{i,i+1} = A_{i,i-1} = 1/dx^2
 *    A_{0,N-1} = A_{N-1,0} = 1/dx^2 (periodic wrap)
 *    b_i = alpha * rho_i
 *  Subtract mean of rho to ensure solvability (integral = 0).
 *  Fix Phi mean = 0.
 * =================================================================== */
static void solve_poisson_periodic(double *Phi, double *rho, int Nx, double dx, double alp)
{
    double dx2 = dx * dx;

    /* Subtract mean from source for solvability */
    double rho_mean = 0;
    for (int i = 0; i < Nx; i++) rho_mean += rho[i];
    rho_mean /= Nx;

    double *b = malloc(Nx * sizeof(double));
    for (int i = 0; i < Nx; i++)
        b[i] = alp * (rho[i] - rho_mean) * dx2;

    /* Sherman-Morrison: decompose A = B + u*v^T
     * B is tridiagonal with B_{0,N-1} = B_{N-1,0} = 0 instead of 1/dx^2
     * Correction: u = (gamma, 0, ..., 0, c/gamma) where gamma = -A_{00}
     *             v = (1, 0, ..., 0, a_{N-1,0}/gamma)
     * Actually simpler: use the standard periodic tridiagonal solver.
     *
     * Standard form: a_i x_{i-1} + b_i x_i + c_i x_{i+1} = d_i
     * with a_0 = c_{N-1} = 1 (the periodic connections).
     * All a_i = c_i = 1, b_i = -2, d_i = b[i] (our RHS).
     */
    double gam = 1.0; /* arbitrary nonzero, use gamma = -b_0 = 2 */
    gam = 2.0;

    /* Modified diagonal for Sherman-Morrison */
    double *bb = malloc(Nx * sizeof(double));
    double *dd = malloc(Nx * sizeof(double));
    double *uu = calloc(Nx, sizeof(double));

    for (int i = 0; i < Nx; i++) {
        bb[i] = -2.0;
        dd[i] = b[i];
    }
    bb[0]    -= gam;
    bb[Nx-1] -= 1.0 / gam;  /* c_{N-1} * a_0 / gamma */

    uu[0]    = gam;
    uu[Nx-1] = 1.0;

    double *vv = calloc(Nx, sizeof(double));
    vv[0]    = 1.0;
    vv[Nx-1] = 1.0 / gam;

    /* Solve B*y = d and B*z = u via Thomas algorithm */
    double *y = malloc(Nx * sizeof(double));
    double *z = malloc(Nx * sizeof(double));

    /* Thomas forward sweep for both RHS */
    double *cc = malloc(Nx * sizeof(double)); /* modified upper diagonal */
    double *dd2 = malloc(Nx * sizeof(double));
    double *uu2 = malloc(Nx * sizeof(double));

    /* Lower diag = 1, main = bb[i], upper = 1 */
    cc[0] = 1.0 / bb[0];
    dd2[0] = dd[0] / bb[0];
    uu2[0] = uu[0] / bb[0];

    for (int i = 1; i < Nx; i++) {
        double w = bb[i] - 1.0 * cc[i-1];  /* bb[i] - a_i * cc[i-1] */
        cc[i]  = 1.0 / w;
        dd2[i] = (dd[i] - 1.0 * dd2[i-1]) / w;
        uu2[i] = (uu[i] - 1.0 * uu2[i-1]) / w;
    }

    /* Back substitution */
    y[Nx-1] = dd2[Nx-1];
    z[Nx-1] = uu2[Nx-1];
    for (int i = Nx - 2; i >= 0; i--) {
        y[i] = dd2[i] - cc[i] * y[i+1];
        z[i] = uu2[i] - cc[i] * z[i+1];
    }

    /* Sherman-Morrison combination: x = y - (v^T y)/(1 + v^T z) * z */
    double vy = 0, vz = 0;
    for (int i = 0; i < Nx; i++) {
        vy += vv[i] * y[i];
        vz += vv[i] * z[i];
    }

    double factor = vy / (1.0 + vz);
    for (int i = 0; i < Nx; i++)
        Phi[i] = y[i] - factor * z[i];

    /* Subtract mean to fix gauge */
    double Phi_mean = 0;
    for (int i = 0; i < Nx; i++) Phi_mean += Phi[i];
    Phi_mean /= Nx;
    for (int i = 0; i < Nx; i++) Phi[i] -= Phi_mean;

    free(b); free(bb); free(dd); free(uu); free(vv);
    free(y); free(z); free(cc); free(dd2); free(uu2);
}

/* ===================================================================
 *  Phase 0: Equilibrate single oscillon (no gravity)
 * =================================================================== */
static void run_equilibrate(void)
{
    printf("===== Phase 0: Single oscillon equilibration (lambda=%.4f) =====\n", lambda);

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

    printf("  Nx=%d xmax=%.1f dx=%.5f dt=%.6f Nt=%d\n", Nx, xmax, dx, dt, Nt);

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
    double prev_phi0 = phi[0][ic], prev_prev_phi0 = phi[0][ic];
    int print_every = Nt / 10;
    if (print_every < 1) print_every = 1;

    double *save_phi[3], *save_vel[3];
    for (int a = 0; a < 3; a++) {
        save_phi[a] = calloc(Nx, sizeof(double));
        save_vel[a] = calloc(Nx, sizeof(double));
    }
    double save_t = -1;
    int n_maxima = 0;
    double t_start_save = 0.5 * t_equil;

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

    eq_mass_osc = 0;
    for (int i = 1; i < Nx - 1; i++)
        eq_mass_osc += energy_density_np(save_phi[0], save_phi[1], save_phi[2],
                                          save_vel[0], save_vel[1], save_vel[2],
                                          i, Nx, dx) * dx;
    printf("  Oscillon mass = %.4f\n", eq_mass_osc);

    eq_omega_osc = 0;
    {
        int dft_start = n_dft / 2;
        if (n_dft - dft_start > 100) {
            (void)(t_hist[n_dft-1] - t_hist[dft_start]); /* T unused */
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

    eq_Nx = Nx;
    eq_dx = dx;
    eq_xmax = xmax;
    for (int a = 0; a < 3; a++) {
        eq_phi[a] = save_phi[a];
        eq_vel[a] = save_vel[a];
    }

    for (int a = 0; a < 3; a++) { free(phi[a]); free(vel[a]); free(acc[a]); }
    free(damp); free(phi0_hist); free(t_hist);
    printf("  Phase 0 complete.\n\n");
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
 *  Compute acceleration with gravity: periodic Laplacian + pairwise
 *  coupling + gravitational modifications.
 *
 *  EOM with gravity:
 *    ddot phi_a = (1+4*Phi) * Laplacian(phi_a) - m^2*(1+2*Phi)*phi_a
 *                 - lambda*(1+2*Phi)*(phi_b + phi_c) + (1+2*Phi)*F_triple
 *
 *  The (1+2*Phi) on mass and potential terms comes from the effective
 *  metric modifying the potential sector; (1+4*Phi) on the Laplacian
 *  from the spatial metric correction.
 * =================================================================== */
static void compute_acc_gravity(double *phi[3], double *acc[3],
                                double *Phi, int Nx, double dx)
{
    double dx2 = dx * dx;
    double m2 = mass * mass;

    for (int a = 0; a < 3; a++) {
        int b = (a+1)%3, c_idx = (a+2)%3;
        for (int i = 0; i < Nx; i++) {
            int ip = (i + 1) % Nx;
            int im = (i - 1 + Nx) % Nx;
            double lapl = (phi[a][ip] - 2.0*phi[a][i] + phi[a][im]) / dx2;
            double fp = force_triple(phi[0][i], phi[1][i], phi[2][i], a);

            double p = Phi[i];
            acc[a][i] = (1.0 + 4.0*p) * lapl
                      - m2 * (1.0 + 2.0*p) * phi[a][i]
                      - lambda * (1.0 + 2.0*p) * (phi[b][i] + phi[c_idx][i])
                      + (1.0 + 2.0*p) * fp;
        }
    }
}

/* Acceleration without gravity (for comparison) */
static void compute_acc_nograv(double *phi[3], double *acc[3], int Nx, double dx)
{
    double dx2 = dx * dx;
    double m2 = mass * mass;

    for (int a = 0; a < 3; a++) {
        int b = (a+1)%3, c_idx = (a+2)%3;
        for (int i = 0; i < Nx; i++) {
            int ip = (i + 1) % Nx;
            int im = (i - 1 + Nx) % Nx;
            double lapl = (phi[a][ip] - 2.0*phi[a][i] + phi[a][im]) / dx2;
            double fp = force_triple(phi[0][i], phi[1][i], phi[2][i], a);
            acc[a][i] = lapl - m2*phi[a][i]
                       - lambda*(phi[b][i] + phi[c_idx][i]) + fp;
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
        else if (!strcmp(argv[i], "-alpha"))    alpha    = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-A"))        A_init   = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-sigma"))    sig      = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-o"))        strncpy(outdir, argv[i+1], sizeof(outdir)-1);
        else if (!strcmp(argv[i], "-N"))        N_osc    = atoi(argv[i+1]);
        else if (!strcmp(argv[i], "-d"))        d_space  = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-seed"))     rng_seed = (unsigned int)atoi(argv[i+1]);
        else if (!strcmp(argv[i], "-delta"))    delta_max = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-t_phase1")) t_phase1 = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-t_phase2")) t_phase2 = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-t_phase3")) t_phase3 = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-t_phase4")) t_phase4 = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-diag_dt"))  diag_dt  = atof(argv[i+1]);
        else { fprintf(stderr, "Unknown flag: %s\n", argv[i]); exit(1); }
    }

    printf("=== Combo 1+2+5+7: Full GR Analog ===\n");
    printf("  lambda=%.4f  alpha=%.6f  mu=%.1f  kappa=%.1f  mass=%.4f\n",
           lambda, alpha, mu, kappa, mass);
    printf("  N_osc=%d  d=%.1f  seed=%u\n", N_osc, d_space, rng_seed);

    /* Step 0: Equilibrate single oscillon */
    run_equilibrate();

    /* Chain setup */
    double L = N_osc * d_space;
    int Nx = N_osc * pts_per_d;
    double dx = L / Nx;
    double m2 = mass * mass;

    double m2_eff = m2 + 2.0 * lambda;
    if (m2_eff < m2) m2_eff = m2;
    double kmax_c = M_PI / dx;
    double dt = 0.8 * 2.0 / sqrt(kmax_c * kmax_c + m2_eff);

    printf("  Chain: L=%.1f Nx=%d dx=%.5f dt=%.6f\n", L, Nx, dx, dt);

    /* Oscillon positions */
    double *x_osc = malloc(N_osc * sizeof(double));
    unsigned int rng = rng_seed;
    for (int n = 0; n < N_osc; n++) {
        x_osc[n] = n * d_space + d_space / 2.0;
        /* No random displacements yet; those come in Phase 2 */
    }

    /* Allocate fields */
    double *phi[3], *vel[3], *acc[3];
    for (int a = 0; a < 3; a++) {
        phi[a] = calloc(Nx, sizeof(double));
        vel[a] = calloc(Nx, sizeof(double));
        acc[a] = calloc(Nx, sizeof(double));
    }
    double *Phi_grav = calloc(Nx, sizeof(double));
    double *rho_dens = calloc(Nx, sizeof(double));

    init_chain(phi, vel, Nx, dx, L, x_osc, N_osc);

    /* Pre-equilibration with damping (no gravity yet) */
    {
        int Nt_pre = (int)(t_pre_equil / dt) + 1;
        printf("\n  Pre-equilibrating chain (gamma=%.4f, t=%.0f, no gravity)...\n",
               gamma_pre, t_pre_equil);

        compute_acc_nograv(phi, acc, Nx, dx);

        for (int step = 0; step < Nt_pre; step++) {
            for (int a = 0; a < 3; a++)
                for (int i = 0; i < Nx; i++)
                    vel[a][i] += 0.5 * dt * acc[a][i];
            for (int a = 0; a < 3; a++)
                for (int i = 0; i < Nx; i++)
                    phi[a][i] += dt * vel[a][i];
            compute_acc_nograv(phi, acc, Nx, dx);
            for (int a = 0; a < 3; a++)
                for (int i = 0; i < Nx; i++)
                    vel[a][i] += 0.5 * dt * acc[a][i];

            /* Position-dependent damping */
            for (int a = 0; a < 3; a++)
                for (int i = 0; i < Nx; i++) {
                    double p2 = phi[0][i]*phi[0][i] + phi[1][i]*phi[1][i]
                              + phi[2][i]*phi[2][i];
                    double protect = p2 / (p2 + 0.01);
                    double local_damp = 1.0 - gamma_pre * dt * (1.0 - protect);
                    vel[a][i] *= local_damp;
                }
        }

        double E_post = 0;
        for (int i = 0; i < Nx; i++)
            E_post += energy_density_p(phi[0], phi[1], phi[2],
                                        vel[0], vel[1], vel[2], i, Nx, dx) * dx;
        printf("  Pre-equilibration done. E=%.4f\n", E_post);
    }

    /* ===================================================================
     *  PHASE 1: Self-gravitating equilibration
     * =================================================================== */
    printf("\n===== PHASE 1: Self-gravitating lattice equilibration (t=%.0f) =====\n", t_phase1);

    int Nt1 = (int)(t_phase1 / dt) + 1;
    int print_every = Nt1 / 20;
    if (print_every < 1) print_every = 1;

    /* Ramp up alpha gradually to avoid shock; use damping for first half */
    double gamma_grav = 0.0005; /* gentle damping during gravity equilibration */
    double t_ramp = t_phase1 * 0.3; /* ramp alpha over first 30% */

    /* Compute initial Phi */
    for (int i = 0; i < Nx; i++)
        rho_dens[i] = energy_density_p(phi[0], phi[1], phi[2],
                                         vel[0], vel[1], vel[2], i, Nx, dx);
    double alpha_eff = alpha * 0.01; /* start small */
    solve_poisson_periodic(Phi_grav, rho_dens, Nx, dx, alpha_eff);
    compute_acc_gravity(phi, acc, Phi_grav, Nx, dx);

    double Phi_min_init = 1e30, Phi_max_init = -1e30;
    for (int i = 0; i < Nx; i++) {
        if (Phi_grav[i] < Phi_min_init) Phi_min_init = Phi_grav[i];
        if (Phi_grav[i] > Phi_max_init) Phi_max_init = Phi_grav[i];
    }
    printf("  Initial Phi: min=%.6e  max=%.6e  range=%.6e\n",
           Phi_min_init, Phi_max_init, Phi_max_init - Phi_min_init);

    /* Time series file for Phase 1 */
    char tspath[600];
    snprintf(tspath, sizeof(tspath), "%s/phase1_ts.tsv", outdir);
    FILE *fts1 = fopen(tspath, "w");
    if (fts1) fprintf(fts1, "time\tE_total\tPhi_min\tPhi_max\tPhi_range\talpha_eff\n");

    for (int step = 0; step <= Nt1; step++) {
        double t = step * dt;

        /* Ramp alpha gradually */
        double ramp = (t < t_ramp) ? t / t_ramp : 1.0;
        alpha_eff = alpha * ramp;

        /* Damping: strong early, zero for last 40% */
        double damp_frac = (t < t_phase1 * 0.6) ? gamma_grav : 0.0;

        if (step % print_every == 0) {
            double Etot = 0;
            for (int i = 0; i < Nx; i++)
                Etot += energy_density_p(phi[0], phi[1], phi[2],
                                          vel[0], vel[1], vel[2], i, Nx, dx) * dx;

            double pmin = 1e30, pmax = -1e30;
            for (int i = 0; i < Nx; i++) {
                if (Phi_grav[i] < pmin) pmin = Phi_grav[i];
                if (Phi_grav[i] > pmax) pmax = Phi_grav[i];
            }

            printf("  t=%7.1f  E=%.4f  Phi=[%.6e, %.6e]  range=%.6e  alpha_eff=%.6e\n",
                   t, Etot, pmin, pmax, pmax - pmin, alpha_eff);

            if (fts1)
                fprintf(fts1, "%.4f\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\n",
                        t, Etot, pmin, pmax, pmax - pmin, alpha_eff);
        }

        if (step == Nt1) break;

        /* Velocity Verlet with Poisson solve at each step */
        for (int a = 0; a < 3; a++)
            for (int i = 0; i < Nx; i++)
                vel[a][i] += 0.5 * dt * acc[a][i];
        for (int a = 0; a < 3; a++)
            for (int i = 0; i < Nx; i++)
                phi[a][i] += dt * vel[a][i];

        /* Recompute rho and Phi */
        for (int i = 0; i < Nx; i++)
            rho_dens[i] = energy_density_p(phi[0], phi[1], phi[2],
                                             vel[0], vel[1], vel[2], i, Nx, dx);
        solve_poisson_periodic(Phi_grav, rho_dens, Nx, dx, alpha_eff);
        compute_acc_gravity(phi, acc, Phi_grav, Nx, dx);

        for (int a = 0; a < 3; a++)
            for (int i = 0; i < Nx; i++)
                vel[a][i] += 0.5 * dt * acc[a][i];

        /* Apply gentle damping (position-dependent: protect cores) */
        if (damp_frac > 0) {
            for (int a = 0; a < 3; a++)
                for (int i = 0; i < Nx; i++) {
                    double p2 = phi[0][i]*phi[0][i] + phi[1][i]*phi[1][i]
                              + phi[2][i]*phi[2][i];
                    double protect = p2 / (p2 + 0.01);
                    double local_damp = 1.0 - damp_frac * dt * (1.0 - protect);
                    vel[a][i] *= local_damp;
                }
        }
    }

    if (fts1) fclose(fts1);

    /* Save Phi profile at end of Phase 1 */
    {
        char path[600];
        snprintf(path, sizeof(path), "%s/phase1_phi_profile.tsv", outdir);
        FILE *fp = fopen(path, "w");
        if (fp) {
            fprintf(fp, "x\tPhi\trho\n");
            for (int i = 0; i < Nx; i += 4)
                fprintf(fp, "%.4f\t%.8e\t%.8e\n", i*dx, Phi_grav[i], rho_dens[i]);
            fclose(fp);
            printf("  Wrote %s\n", path);
        }
    }

    /* Track positions at end of Phase 1 */
    find_positions_periodic(phi[0], phi[1], phi[2], vel[0], vel[1], vel[2],
                             Nx, dx, L, x_osc, N_osc);
    printf("  Phase 1 positions:\n");
    for (int n = 0; n < N_osc; n++)
        printf("    n=%d: x=%.3f\n", n, x_osc[n]);

    double E_phase1 = 0;
    for (int i = 0; i < Nx; i++)
        E_phase1 += energy_density_p(phi[0], phi[1], phi[2],
                                      vel[0], vel[1], vel[2], i, Nx, dx) * dx;
    printf("  Phase 1 final energy: %.4f\n", E_phase1);

    /* ===================================================================
     *  PHASE 2: Phonon spectrum with gravity
     *  Save state, then run WITH gravity. Also run a copy WITHOUT gravity.
     * =================================================================== */
    printf("\n===== PHASE 2: Phonon spectrum (t=%.0f) =====\n", t_phase2);

    /* Save state for no-gravity comparison */
    double *phi_ng[3], *vel_ng[3], *acc_ng[3];
    for (int a = 0; a < 3; a++) {
        phi_ng[a] = malloc(Nx * sizeof(double));
        vel_ng[a] = malloc(Nx * sizeof(double));
        acc_ng[a] = malloc(Nx * sizeof(double));
        memcpy(phi_ng[a], phi[a], Nx * sizeof(double));
        memcpy(vel_ng[a], vel[a], Nx * sizeof(double));
    }
    double *x_osc_ng = malloc(N_osc * sizeof(double));
    memcpy(x_osc_ng, x_osc, N_osc * sizeof(double));

    /* Add random position perturbations by velocity kicks */
    printf("  Adding random phonon perturbations (delta=%.3f)...\n", delta_max);
    double A_kick = 0.01;
    rng = rng_seed + 12345;
    for (int i = 0; i < Nx; i++) {
        double kick = A_kick * (2.0 * rand01(&rng) - 1.0);
        for (int a = 0; a < 3; a++) {
            vel[a][i]    += kick * phi[a][i];
            vel_ng[a][i] += kick * phi_ng[a][i];
        }
    }

    int Nt2 = (int)(t_phase2 / dt) + 1;
    int diag_every = (int)(diag_dt / dt);
    if (diag_every < 1) diag_every = 1;
    int max_diag = Nt2 / diag_every + 10;

    /* Position tracking */
    double *diag_t2 = malloc(max_diag * sizeof(double));
    double **diag_x_grav = malloc(N_osc * sizeof(double*));
    double **diag_x_nograv = malloc(N_osc * sizeof(double*));
    for (int n = 0; n < N_osc; n++) {
        diag_x_grav[n] = malloc(max_diag * sizeof(double));
        diag_x_nograv[n] = malloc(max_diag * sizeof(double));
    }
    double *diag_E_grav = malloc(max_diag * sizeof(double));
    double *diag_E_nograv = malloc(max_diag * sizeof(double));
    int n_diag2 = 0;

    /* Initial acceleration */
    for (int i = 0; i < Nx; i++)
        rho_dens[i] = energy_density_p(phi[0], phi[1], phi[2],
                                         vel[0], vel[1], vel[2], i, Nx, dx);
    solve_poisson_periodic(Phi_grav, rho_dens, Nx, dx, alpha);
    compute_acc_gravity(phi, acc, Phi_grav, Nx, dx);
    compute_acc_nograv(phi_ng, acc_ng, Nx, dx);

    print_every = Nt2 / 20;
    if (print_every < 1) print_every = 1;

    for (int step = 0; step <= Nt2; step++) {
        double t = step * dt;

        if (step % diag_every == 0 && n_diag2 < max_diag) {
            if (step > 0) {
                find_positions_periodic(phi[0], phi[1], phi[2],
                                         vel[0], vel[1], vel[2],
                                         Nx, dx, L, x_osc, N_osc);
                find_positions_periodic(phi_ng[0], phi_ng[1], phi_ng[2],
                                         vel_ng[0], vel_ng[1], vel_ng[2],
                                         Nx, dx, L, x_osc_ng, N_osc);
            }

            double Eg = 0, Eng = 0;
            for (int i = 0; i < Nx; i++) {
                Eg += energy_density_p(phi[0], phi[1], phi[2],
                                        vel[0], vel[1], vel[2], i, Nx, dx) * dx;
                Eng += energy_density_p(phi_ng[0], phi_ng[1], phi_ng[2],
                                         vel_ng[0], vel_ng[1], vel_ng[2], i, Nx, dx) * dx;
            }

            diag_t2[n_diag2] = t;
            for (int n = 0; n < N_osc; n++) {
                diag_x_grav[n][n_diag2] = x_osc[n];
                diag_x_nograv[n][n_diag2] = x_osc_ng[n];
            }
            diag_E_grav[n_diag2] = Eg;
            diag_E_nograv[n_diag2] = Eng;
            n_diag2++;
        }

        if (step % print_every == 0) {
            double Eg = 0;
            for (int i = 0; i < Nx; i++)
                Eg += energy_density_p(phi[0], phi[1], phi[2],
                                        vel[0], vel[1], vel[2], i, Nx, dx) * dx;
            printf("  t=%7.1f  E_grav=%.4f\n", t, Eg);
        }

        if (step == Nt2) break;

        /* Evolve WITH gravity */
        for (int a = 0; a < 3; a++)
            for (int i = 0; i < Nx; i++)
                vel[a][i] += 0.5 * dt * acc[a][i];
        for (int a = 0; a < 3; a++)
            for (int i = 0; i < Nx; i++)
                phi[a][i] += dt * vel[a][i];
        for (int i = 0; i < Nx; i++)
            rho_dens[i] = energy_density_p(phi[0], phi[1], phi[2],
                                             vel[0], vel[1], vel[2], i, Nx, dx);
        solve_poisson_periodic(Phi_grav, rho_dens, Nx, dx, alpha);
        compute_acc_gravity(phi, acc, Phi_grav, Nx, dx);
        for (int a = 0; a < 3; a++)
            for (int i = 0; i < Nx; i++)
                vel[a][i] += 0.5 * dt * acc[a][i];

        /* Evolve WITHOUT gravity */
        for (int a = 0; a < 3; a++)
            for (int i = 0; i < Nx; i++)
                vel_ng[a][i] += 0.5 * dt * acc_ng[a][i];
        for (int a = 0; a < 3; a++)
            for (int i = 0; i < Nx; i++)
                phi_ng[a][i] += dt * vel_ng[a][i];
        compute_acc_nograv(phi_ng, acc_ng, Nx, dx);
        for (int a = 0; a < 3; a++)
            for (int i = 0; i < Nx; i++)
                vel_ng[a][i] += 0.5 * dt * acc_ng[a][i];
    }

    printf("  Phase 2 complete. %d diagnostic samples.\n", n_diag2);

    /* Phonon spectrum analysis: compute displacements and normal modes */
    int j_start2 = n_diag2 / 4; /* discard first quarter */
    int n_an2 = n_diag2 - j_start2;

    if (n_an2 > 20) {
        /* Displacements u_n(t) for both grav and nograv */
        double **u_grav = malloc(N_osc * sizeof(double*));
        double **u_nograv = malloc(N_osc * sizeof(double*));

        for (int n = 0; n < N_osc; n++) {
            u_grav[n] = calloc(n_an2, sizeof(double));
            u_nograv[n] = calloc(n_an2, sizeof(double));

            double prev_g = diag_x_grav[n][j_start2];
            double prev_ng = diag_x_nograv[n][j_start2];
            u_grav[n][0] = 0;
            u_nograv[n][0] = 0;

            for (int j = 1; j < n_an2; j++) {
                int jj = j + j_start2;
                double dxg = diag_x_grav[n][jj] - prev_g;
                while (dxg > L/2) dxg -= L;
                while (dxg < -L/2) dxg += L;
                u_grav[n][j] = u_grav[n][j-1] + dxg;
                prev_g = diag_x_grav[n][jj];

                double dxng = diag_x_nograv[n][jj] - prev_ng;
                while (dxng > L/2) dxng -= L;
                while (dxng < -L/2) dxng += L;
                u_nograv[n][j] = u_nograv[n][j-1] + dxng;
                prev_ng = diag_x_nograv[n][jj];
            }
        }

        /* Remove COM drift */
        for (int j = 0; j < n_an2; j++) {
            double com_g = 0, com_ng = 0;
            for (int n = 0; n < N_osc; n++) {
                com_g += u_grav[n][j];
                com_ng += u_nograv[n][j];
            }
            com_g /= N_osc; com_ng /= N_osc;
            for (int n = 0; n < N_osc; n++) {
                u_grav[n][j] -= com_g;
                u_nograv[n][j] -= com_ng;
            }
        }

        /* Normal mode decomposition Q_q(t) */
        double T_an = diag_t2[n_diag2-1] - diag_t2[j_start2];
        double dt_diag = T_an / (n_an2 - 1);

        /* DFT for each mode q to get omega_peak */
        int n_omega = n_an2 / 2;
        if (n_omega > 4000) n_omega = 4000;
        double omega_max = 1.0;
        double omega_nyq = M_PI / dt_diag;
        if (omega_max > omega_nyq) omega_max = omega_nyq;

        double *omega_grav = calloc(N_osc, sizeof(double));
        double *omega_nograv = calloc(N_osc, sizeof(double));
        double *power_grav = calloc(N_osc, sizeof(double));
        double *power_nograv = calloc(N_osc, sizeof(double));

        for (int q = 0; q <= N_osc / 2; q++) {
            double best_pg = 0, best_og = 0;
            double best_png = 0, best_ong = 0;

            for (int l = 1; l <= n_omega; l++) {
                double omega = omega_max * l / n_omega;

                /* Q_q for grav */
                double re_g = 0, im_g = 0;
                double re_ng = 0, im_ng = 0;

                for (int j = 0; j < n_an2; j++) {
                    double Q_re_g = 0, Q_im_g = 0;
                    double Q_re_ng = 0, Q_im_ng = 0;
                    for (int n = 0; n < N_osc; n++) {
                        double phase = -2.0 * M_PI * q * n / N_osc;
                        double cp = cos(phase), sp = sin(phase);
                        Q_re_g += u_grav[n][j] * cp;
                        Q_im_g += u_grav[n][j] * sp;
                        Q_re_ng += u_nograv[n][j] * cp;
                        Q_im_ng += u_nograv[n][j] * sp;
                    }
                    Q_re_g /= sqrt((double)N_osc);
                    Q_im_g /= sqrt((double)N_osc);
                    Q_re_ng /= sqrt((double)N_osc);
                    Q_im_ng /= sqrt((double)N_osc);

                    double tj = j * dt_diag;
                    double c = cos(omega * tj), s = sin(omega * tj);
                    re_g += (Q_re_g * c + Q_im_g * s) * dt_diag;
                    im_g += (-Q_re_g * s + Q_im_g * c) * dt_diag;
                    re_ng += (Q_re_ng * c + Q_im_ng * s) * dt_diag;
                    im_ng += (-Q_re_ng * s + Q_im_ng * c) * dt_diag;
                }

                double pg = (re_g*re_g + im_g*im_g) / (T_an*T_an);
                double png = (re_ng*re_ng + im_ng*im_ng) / (T_an*T_an);

                if (pg > best_pg) { best_pg = pg; best_og = omega; }
                if (png > best_png) { best_png = png; best_ong = omega; }
            }

            omega_grav[q] = best_og;
            omega_nograv[q] = best_ong;
            power_grav[q] = best_pg;
            power_nograv[q] = best_png;
        }

        /* Write phonon comparison */
        {
            char path[600];
            snprintf(path, sizeof(path), "%s/phase2_phonons.tsv", outdir);
            FILE *fp = fopen(path, "w");
            if (fp) {
                fprintf(fp, "q\tk\tomega_grav\tomega_nograv\tdomega_pct\tcs_grav\tcs_nograv\n");
                for (int q = 0; q <= N_osc / 2; q++) {
                    double k_q = 2.0 * M_PI * q / (N_osc * d_space);
                    double cs_g = (k_q > 1e-10) ? omega_grav[q] / k_q : 0;
                    double cs_ng = (k_q > 1e-10) ? omega_nograv[q] / k_q : 0;
                    double domega = (omega_nograv[q] > 1e-10) ?
                        100.0*(omega_grav[q] - omega_nograv[q])/omega_nograv[q] : 0;
                    fprintf(fp, "%d\t%.6f\t%.6f\t%.6f\t%.4f\t%.4f\t%.4f\n",
                            q, k_q, omega_grav[q], omega_nograv[q], domega, cs_g, cs_ng);
                }
                fclose(fp);
                printf("  Wrote %s\n", path);
            }
        }

        printf("\n  Phonon dispersion comparison:\n");
        printf("  %-4s  %-10s  %-10s  %-10s  %-10s  %-10s  %-10s\n",
               "q", "k", "w_grav", "w_nograv", "dw(%)", "cs_grav", "cs_nograv");
        printf("  ----  ----------  ----------  ----------  ----------  ----------  ----------\n");
        for (int q = 0; q <= N_osc / 2; q++) {
            double k_q = 2.0 * M_PI * q / (N_osc * d_space);
            double cs_g = (k_q > 1e-10) ? omega_grav[q] / k_q : 0;
            double cs_ng = (k_q > 1e-10) ? omega_nograv[q] / k_q : 0;
            double domega = (omega_nograv[q] > 1e-10) ?
                100.0*(omega_grav[q] - omega_nograv[q])/omega_nograv[q] : 0;
            printf("  %-4d  %-10.5f  %-10.5f  %-10.5f  %-+10.4f  %-10.4f  %-10.4f\n",
                   q, k_q, omega_grav[q], omega_nograv[q], domega, cs_g, cs_ng);
        }

        if (N_osc >= 2) {
            double k1 = 2.0 * M_PI / (N_osc * d_space);
            double cs_g1 = omega_grav[1] / k1;
            double cs_ng1 = omega_nograv[1] / k1;
            printf("\n  Sound speed (q=1): c_s_grav=%.4f  c_s_nograv=%.4f  delta=%.4f%%\n",
                   cs_g1, cs_ng1, 100.0*(cs_g1-cs_ng1)/(cs_ng1+1e-30));
        }

        /* Clean up Phase 2 analysis arrays */
        for (int n = 0; n < N_osc; n++) { free(u_grav[n]); free(u_nograv[n]); }
        free(u_grav); free(u_nograv);
        free(omega_grav); free(omega_nograv);
        free(power_grav); free(power_nograv);
    }

    /* Free no-grav arrays */
    for (int a = 0; a < 3; a++) {
        free(phi_ng[a]); free(vel_ng[a]); free(acc_ng[a]);
    }
    free(x_osc_ng);

    /* Free Phase 2 diagnostics */
    free(diag_t2); free(diag_E_grav); free(diag_E_nograv);
    for (int n = 0; n < N_osc; n++) { free(diag_x_grav[n]); free(diag_x_nograv[n]); }
    free(diag_x_grav); free(diag_x_nograv);

    /* ===================================================================
     *  PHASE 3: Phi propagation test
     *  Perturb oscillon #0, track Phi(x,t)
     * =================================================================== */
    printf("\n===== PHASE 3: Phi propagation test (t=%.0f) =====\n", t_phase3);

    /* Save Phi profile before perturbation */
    double *Phi_before = malloc(Nx * sizeof(double));
    memcpy(Phi_before, Phi_grav, Nx * sizeof(double));

    /* Perturb oscillon #0: kick velocity */
    double pert_amp = 0.05;
    int pert_osc = 0;
    find_positions_periodic(phi[0], phi[1], phi[2], vel[0], vel[1], vel[2],
                             Nx, dx, L, x_osc, N_osc);
    double xc_pert = x_osc[pert_osc];
    printf("  Perturbing oscillon #%d at x=%.2f with amplitude %.4f\n",
           pert_osc, xc_pert, pert_amp);

    for (int i = 0; i < Nx; i++) {
        double x = i * dx;
        double rel = x - xc_pert;
        while (rel > L/2) rel -= L;
        while (rel < -L/2) rel += L;
        double env = exp(-rel*rel / (2.0 * sig * sig));
        for (int a = 0; a < 3; a++)
            vel[a][i] += pert_amp * env * phi[a][i];
    }

    int Nt3 = (int)(t_phase3 / dt) + 1;
    int snap_every = Nt3 / 50;
    if (snap_every < 1) snap_every = 1;

    /* Store Phi snapshots at oscillon positions */
    int max_snap = 60;
    double *snap_t = malloc(max_snap * sizeof(double));
    double **snap_Phi_osc = malloc(N_osc * sizeof(double*));
    for (int n = 0; n < N_osc; n++)
        snap_Phi_osc[n] = malloc(max_snap * sizeof(double));
    int n_snap = 0;

    /* Also store full Phi(x) at a few times for wavefront detection */
    int n_full_snap = 0;
    int max_full = 10;
    double *full_snap_t = malloc(max_full * sizeof(double));
    double **full_snap_Phi = malloc(max_full * sizeof(double*));
    for (int k = 0; k < max_full; k++)
        full_snap_Phi[k] = malloc(Nx * sizeof(double));
    int full_snap_interval = Nt3 / max_full;
    if (full_snap_interval < 1) full_snap_interval = 1;

    /* Recompute */
    for (int i = 0; i < Nx; i++)
        rho_dens[i] = energy_density_p(phi[0], phi[1], phi[2],
                                         vel[0], vel[1], vel[2], i, Nx, dx);
    solve_poisson_periodic(Phi_grav, rho_dens, Nx, dx, alpha);
    compute_acc_gravity(phi, acc, Phi_grav, Nx, dx);

    print_every = Nt3 / 10;
    if (print_every < 1) print_every = 1;

    for (int step = 0; step <= Nt3; step++) {
        double t = step * dt;

        if (step % snap_every == 0 && n_snap < max_snap) {
            find_positions_periodic(phi[0], phi[1], phi[2], vel[0], vel[1], vel[2],
                                     Nx, dx, L, x_osc, N_osc);
            snap_t[n_snap] = t;
            for (int n = 0; n < N_osc; n++) {
                int idx = (int)(x_osc[n] / dx);
                if (idx < 0) idx = 0;
                if (idx >= Nx) idx = Nx - 1;
                snap_Phi_osc[n][n_snap] = Phi_grav[idx];
            }
            n_snap++;
        }

        if (step % full_snap_interval == 0 && n_full_snap < max_full) {
            full_snap_t[n_full_snap] = t;
            memcpy(full_snap_Phi[n_full_snap], Phi_grav, Nx * sizeof(double));
            n_full_snap++;
        }

        if (step % print_every == 0) {
            double pmin = 1e30, pmax = -1e30;
            for (int i = 0; i < Nx; i++) {
                if (Phi_grav[i] < pmin) pmin = Phi_grav[i];
                if (Phi_grav[i] > pmax) pmax = Phi_grav[i];
            }
            printf("  t=%7.1f  Phi=[%.6e, %.6e]\n", t, pmin, pmax);
        }

        if (step == Nt3) break;

        for (int a = 0; a < 3; a++)
            for (int i = 0; i < Nx; i++)
                vel[a][i] += 0.5 * dt * acc[a][i];
        for (int a = 0; a < 3; a++)
            for (int i = 0; i < Nx; i++)
                phi[a][i] += dt * vel[a][i];
        for (int i = 0; i < Nx; i++)
            rho_dens[i] = energy_density_p(phi[0], phi[1], phi[2],
                                             vel[0], vel[1], vel[2], i, Nx, dx);
        solve_poisson_periodic(Phi_grav, rho_dens, Nx, dx, alpha);
        compute_acc_gravity(phi, acc, Phi_grav, Nx, dx);
        for (int a = 0; a < 3; a++)
            for (int i = 0; i < Nx; i++)
                vel[a][i] += 0.5 * dt * acc[a][i];
    }

    /* Write Phi at oscillon positions over time */
    {
        char path[600];
        snprintf(path, sizeof(path), "%s/phase3_phi_propagation.tsv", outdir);
        FILE *fp = fopen(path, "w");
        if (fp) {
            fprintf(fp, "time");
            for (int n = 0; n < N_osc; n++) fprintf(fp, "\tPhi_%d", n);
            fprintf(fp, "\n");
            for (int j = 0; j < n_snap; j++) {
                fprintf(fp, "%.4f", snap_t[j]);
                for (int n = 0; n < N_osc; n++)
                    fprintf(fp, "\t%.8e", snap_Phi_osc[n][j]);
                fprintf(fp, "\n");
            }
            fclose(fp);
            printf("  Wrote %s\n", path);
        }
    }

    /* Write full Phi snapshots */
    {
        char path[600];
        snprintf(path, sizeof(path), "%s/phase3_phi_snapshots.tsv", outdir);
        FILE *fp = fopen(path, "w");
        if (fp) {
            fprintf(fp, "x");
            for (int k = 0; k < n_full_snap; k++) fprintf(fp, "\tPhi_t%.0f", full_snap_t[k]);
            fprintf(fp, "\n");
            for (int i = 0; i < Nx; i += 8) {
                fprintf(fp, "%.4f", i * dx);
                for (int k = 0; k < n_full_snap; k++)
                    fprintf(fp, "\t%.8e", full_snap_Phi[k][i]);
                fprintf(fp, "\n");
            }
            fclose(fp);
            printf("  Wrote %s\n", path);
        }
    }

    /* Analyze Phi propagation: compute delta_Phi = Phi(t) - Phi(t=0) at each oscillon
     * and check if perturbation arrives at distant oscillons with a delay */
    printf("\n  Phi propagation analysis:\n");
    printf("  %-6s  %-12s  %-12s  %-12s  %-12s\n",
           "osc", "dist_from_0", "dPhi_peak", "t_peak", "v_apparent");
    printf("  ------  ------------  ------------  ------------  ------------\n");

    for (int n = 0; n < N_osc; n++) {
        double dist = fabs(n - pert_osc) * d_space;
        if (dist > L/2) dist = L - dist;

        /* Find time of max |delta Phi| */
        double dphi_peak = 0, t_peak_n = 0;
        for (int j = 1; j < n_snap; j++) {
            double dphi = fabs(snap_Phi_osc[n][j] - snap_Phi_osc[n][0]);
            if (dphi > dphi_peak) {
                dphi_peak = dphi;
                t_peak_n = snap_t[j];
            }
        }
        double v_app = (t_peak_n > 1e-10 && dist > 1e-10) ? dist / t_peak_n : 0;
        printf("  %-6d  %-12.2f  %-12.6e  %-12.2f  %-12.4f\n",
               n, dist, dphi_peak, t_peak_n, v_app);
    }

    free(Phi_before); free(snap_t);
    for (int n = 0; n < N_osc; n++) free(snap_Phi_osc[n]);
    free(snap_Phi_osc);
    free(full_snap_t);
    for (int k = 0; k < max_full; k++) free(full_snap_Phi[k]);
    free(full_snap_Phi);

    /* ===================================================================
     *  PHASE 4: Equivalence principle test
     *  Perturb oscillons #2 and #5 with DIFFERENT amplitudes.
     *  Track if they experience same acceleration in same Phi gradient.
     * =================================================================== */
    printf("\n===== PHASE 4: Equivalence principle test (t=%.0f) =====\n", t_phase4);

    /* Perturb oscillon #2 (small) and #5 (large) */
    int osc_A = 2, osc_B = 5;
    double amp_A = 0.02, amp_B = 0.10;

    find_positions_periodic(phi[0], phi[1], phi[2], vel[0], vel[1], vel[2],
                             Nx, dx, L, x_osc, N_osc);

    printf("  Perturbing osc #%d (amp=%.4f) and osc #%d (amp=%.4f)\n",
           osc_A, amp_A, osc_B, amp_B);

    /* Apply perturbations */
    for (int i = 0; i < Nx; i++) {
        double x = i * dx;
        for (int p = 0; p < 2; p++) {
            int osc_n = (p == 0) ? osc_A : osc_B;
            double amp = (p == 0) ? amp_A : amp_B;
            double xc = x_osc[osc_n];
            double rel = x - xc;
            while (rel > L/2) rel -= L;
            while (rel < -L/2) rel += L;
            double env = exp(-rel*rel / (2.0 * sig * sig));
            for (int a = 0; a < 3; a++)
                vel[a][i] += amp * env * phi[a][i];
        }
    }

    /* Recompute */
    for (int i = 0; i < Nx; i++)
        rho_dens[i] = energy_density_p(phi[0], phi[1], phi[2],
                                         vel[0], vel[1], vel[2], i, Nx, dx);
    solve_poisson_periodic(Phi_grav, rho_dens, Nx, dx, alpha);
    compute_acc_gravity(phi, acc, Phi_grav, Nx, dx);

    int Nt4 = (int)(t_phase4 / dt) + 1;
    diag_every = (int)(diag_dt / dt);
    if (diag_every < 1) diag_every = 1;
    int max_diag4 = Nt4 / diag_every + 10;

    double *diag_t4 = malloc(max_diag4 * sizeof(double));
    double *diag_xA = malloc(max_diag4 * sizeof(double));
    double *diag_xB = malloc(max_diag4 * sizeof(double));
    double *diag_PhiA = malloc(max_diag4 * sizeof(double));
    double *diag_PhiB = malloc(max_diag4 * sizeof(double));
    double *diag_EA = malloc(max_diag4 * sizeof(double));
    double *diag_EB = malloc(max_diag4 * sizeof(double));
    int n_diag4 = 0;

    double x_osc_init_A = x_osc[osc_A];
    double x_osc_init_B = x_osc[osc_B];

    print_every = Nt4 / 10;
    if (print_every < 1) print_every = 1;

    for (int step = 0; step <= Nt4; step++) {
        double t = step * dt;

        if (step % diag_every == 0 && n_diag4 < max_diag4) {
            if (step > 0)
                find_positions_periodic(phi[0], phi[1], phi[2], vel[0], vel[1], vel[2],
                                         Nx, dx, L, x_osc, N_osc);

            diag_t4[n_diag4] = t;
            diag_xA[n_diag4] = x_osc[osc_A];
            diag_xB[n_diag4] = x_osc[osc_B];

            int idxA = (int)(x_osc[osc_A] / dx);
            int idxB = (int)(x_osc[osc_B] / dx);
            if (idxA < 0) idxA = 0;
            if (idxA >= Nx) idxA = Nx - 1;
            if (idxB < 0) idxB = 0;
            if (idxB >= Nx) idxB = Nx - 1;
            diag_PhiA[n_diag4] = Phi_grav[idxA];
            diag_PhiB[n_diag4] = Phi_grav[idxB];

            /* Local energy of each oscillon */
            double EA = 0, EB = 0;
            for (int i = 0; i < Nx; i++) {
                double x = i * dx;
                double relA = x - x_osc[osc_A];
                while (relA > L/2) relA -= L;
                while (relA < -L/2) relA += L;
                double relB = x - x_osc[osc_B];
                while (relB > L/2) relB -= L;
                while (relB < -L/2) relB += L;

                double e = energy_density_p(phi[0], phi[1], phi[2],
                                             vel[0], vel[1], vel[2], i, Nx, dx);
                if (fabs(relA) < d_space/2) EA += e * dx;
                if (fabs(relB) < d_space/2) EB += e * dx;
            }
            diag_EA[n_diag4] = EA;
            diag_EB[n_diag4] = EB;
            n_diag4++;
        }

        if (step % print_every == 0) {
            printf("  t=%7.1f  xA=%.3f  xB=%.3f\n", t, x_osc[osc_A], x_osc[osc_B]);
        }

        if (step == Nt4) break;

        for (int a = 0; a < 3; a++)
            for (int i = 0; i < Nx; i++)
                vel[a][i] += 0.5 * dt * acc[a][i];
        for (int a = 0; a < 3; a++)
            for (int i = 0; i < Nx; i++)
                phi[a][i] += dt * vel[a][i];
        for (int i = 0; i < Nx; i++)
            rho_dens[i] = energy_density_p(phi[0], phi[1], phi[2],
                                             vel[0], vel[1], vel[2], i, Nx, dx);
        solve_poisson_periodic(Phi_grav, rho_dens, Nx, dx, alpha);
        compute_acc_gravity(phi, acc, Phi_grav, Nx, dx);
        for (int a = 0; a < 3; a++)
            for (int i = 0; i < Nx; i++)
                vel[a][i] += 0.5 * dt * acc[a][i];
    }

    /* Write equivalence data */
    {
        char path[600];
        snprintf(path, sizeof(path), "%s/phase4_equivalence.tsv", outdir);
        FILE *fp = fopen(path, "w");
        if (fp) {
            fprintf(fp, "time\txA\txB\tPhiA\tPhiB\tEA\tEB\tuA\tuB\n");
            for (int j = 0; j < n_diag4; j++) {
                double uA = diag_xA[j] - x_osc_init_A;
                while (uA > L/2) uA -= L;
                while (uA < -L/2) uA += L;
                double uB = diag_xB[j] - x_osc_init_B;
                while (uB > L/2) uB -= L;
                while (uB < -L/2) uB += L;
                fprintf(fp, "%.4f\t%.6f\t%.6f\t%.8e\t%.8e\t%.6e\t%.6e\t%.8e\t%.8e\n",
                        diag_t4[j], diag_xA[j], diag_xB[j],
                        diag_PhiA[j], diag_PhiB[j],
                        diag_EA[j], diag_EB[j], uA, uB);
            }
            fclose(fp);
            printf("  Wrote %s\n", path);
        }
    }

    /* Analyze equivalence: compute effective acceleration from u(t) ≈ (1/2)a*t^2 */
    printf("\n  Equivalence principle analysis:\n");
    if (n_diag4 > 10) {
        /* Use first 20% of data for quadratic fit: u = a0 + a1*t + (1/2)*a2*t^2 */
        int n_fit = n_diag4 / 5;
        if (n_fit < 5) n_fit = 5;
        if (n_fit > n_diag4) n_fit = n_diag4;

        /* Simple: compute average acceleration = 2*u_final / t_final^2 */
        double uA_end = diag_xA[n_fit-1] - x_osc_init_A;
        while (uA_end > L/2) uA_end -= L;
        while (uA_end < -L/2) uA_end += L;
        double uB_end = diag_xB[n_fit-1] - x_osc_init_B;
        while (uB_end > L/2) uB_end -= L;
        while (uB_end < -L/2) uB_end += L;

        double t_fit = diag_t4[n_fit-1] - diag_t4[0];
        double a_eff_A = 2.0 * uA_end / (t_fit * t_fit);
        double a_eff_B = 2.0 * uB_end / (t_fit * t_fit);

        printf("  Oscillon A (#%d): amp=%.4f, E=%.4f, u_end=%.6e, a_eff=%.6e\n",
               osc_A, amp_A, diag_EA[0], uA_end, a_eff_A);
        printf("  Oscillon B (#%d): amp=%.4f, E=%.4f, u_end=%.6e, a_eff=%.6e\n",
               osc_B, amp_B, diag_EB[0], uB_end, a_eff_B);

        /* Compare Phi gradients at the two positions */
        printf("  Phi at A: %.6e, Phi at B: %.6e\n", diag_PhiA[0], diag_PhiB[0]);

        if (fabs(a_eff_A) > 1e-20 && fabs(a_eff_B) > 1e-20) {
            double ratio = a_eff_A / a_eff_B;
            printf("  a_A / a_B = %.4f\n", ratio);
            printf("  E_A / E_B = %.4f\n", diag_EA[0] / (diag_EB[0] + 1e-30));
            printf("  Equivalence principle: %s\n",
                   (fabs(ratio - 1.0) < 0.2) ? "CONSISTENT (ratio near 1)" :
                   "VIOLATED (ratio far from 1)");
        } else {
            printf("  Accelerations too small to compare meaningfully.\n");
            printf("  |a_A|=%.2e, |a_B|=%.2e\n", fabs(a_eff_A), fabs(a_eff_B));
        }

        /* Also check: do both oscillons see the same Phi gradient? */
        printf("\n  Phi gradient comparison (averaged over evolution):\n");
        double dPhi_A_sum = 0, dPhi_B_sum = 0;
        int n_grad = 0;
        for (int j = 0; j < n_diag4; j++) {
            int idxA = (int)(diag_xA[j] / dx);
            int idxB = (int)(diag_xB[j] / dx);
            if (idxA >= 1 && idxA < Nx-1 && idxB >= 1 && idxB < Nx-1) {
                /* This is a rough estimate; Phi is only stored at snap times */
                dPhi_A_sum += diag_PhiA[j];
                dPhi_B_sum += diag_PhiB[j];
                n_grad++;
            }
        }
        if (n_grad > 0) {
            printf("  <Phi_A> = %.6e\n", dPhi_A_sum / n_grad);
            printf("  <Phi_B> = %.6e\n", dPhi_B_sum / n_grad);
        }
    }

    free(diag_t4); free(diag_xA); free(diag_xB);
    free(diag_PhiA); free(diag_PhiB);
    free(diag_EA); free(diag_EB);

    /* ===================================================================
     *  SUMMARY
     * =================================================================== */
    printf("\n========================================\n");
    printf("  SUMMARY: Combo 1+2+5+7 (Full GR Analog)\n");
    printf("========================================\n");
    printf("  Parameters: lambda=%.4f, alpha=%.6f, N_osc=%d, d=%.1f\n",
           lambda, alpha, N_osc, d_space);
    printf("  Oscillon mass = %.4f\n", eq_mass_osc);
    printf("  Breathing freq = %.4f\n", eq_omega_osc);
    printf("\n  KEY QUESTIONS:\n");
    printf("  1. Self-gravitating lattice stable? (Phase 1 survived)\n");
    printf("  2. Phonon spectrum modified by gravity? (Phase 2 comparison)\n");
    printf("  3. Does Phi propagate? (Phase 3 wavefront analysis)\n");
    printf("  4. Equivalence principle? (Phase 4 acceleration comparison)\n");

    /* Cleanup */
    for (int a = 0; a < 3; a++) {
        free(phi[a]); free(vel[a]); free(acc[a]);
    }
    free(Phi_grav); free(rho_dens);
    free(x_osc);
    for (int a = 0; a < 3; a++) {
        if (eq_phi[a]) free(eq_phi[a]);
        if (eq_vel[a]) free(eq_vel[a]);
    }

    printf("\n=== Combo 1+2+5+7 Complete ===\n");
    return 0;
}
