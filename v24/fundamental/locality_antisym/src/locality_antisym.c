/*
 * locality_antisym.c -- Locality Option 2: Antisymmetric Mode as Causal Gravity
 *
 * Does the near-gapless antisymmetric mode from the confining potential
 * mediate a FORCE between oscillons?  Key test: apply antisymmetric
 * perturbation at oscillon #1 only, track whether oscillons #2-8 feel
 * a force (center displacement) and with what time delay.
 *
 * Potential: V = -sigma*sqrt(P^2+eps^2) + (kc/2)*P^4  (confining + quartic)
 * Pairwise:  lambda*(phi1*phi2 + phi2*phi3 + phi3*phi1)
 * P = phi1*phi2*phi3
 *
 * Phases:
 *   1: Equilibrate single oscillon (t=5000, absorbing BC)
 *   2: Build 8-oscillon periodic chain, pre-equilibrate, settle (t=5000)
 *   3: Antisymmetric perturbation at oscillon #1, evolve t=10000
 *      Track center displacements + antisymmetric content at all oscillons
 *   4: Two isolated oscillons at D=40, same test
 *
 * Compile: gcc -O3 -Wall -o locality_antisym src/locality_antisym.c -lm
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* --- Parameters --- */
static double mass     = 1.0;
static double A_init   = 1.0;
static double sig      = 3.0;     /* Gaussian width for initialization */
static char   outdir[512] = "data";

/* Confining potential parameters */
static double sigma_c  = 1.0;     /* string tension */
static double eps_reg  = 1e-6;    /* regularization */
static double kc       = 0.01;    /* quartic stabilizer */

/* Pairwise coupling */
static double lambda   = 0.5;

/* Phase 1: equilibration */
static int    Nx_eq    = 4000;
static double xmax_eq  = 100.0;
static double t_equil  = 5000.0;

/* Phase 2+3: chain */
static int    N_osc    = 8;
static double d_space  = 16.0;
static int    pts_per_d = 320;
static double t_pre_equil = 500.0;
static double gamma_pre = 0.001;
static double t_settle  = 5000.0;
static double t_evolve  = 10000.0;

/* Perturbation */
static double eps_pert  = 0.01;
static int    pert_osc  = 0;      /* oscillon #1 (0-indexed) */
static double pert_sigma = 3.0;
static double diag_dt   = 5.0;

/* Phase 4: two isolated oscillons */
static int    Nx_iso    = 8000;
static double xmax_iso  = 200.0;
static double D_iso     = 40.0;
static double t_settle_iso = 2000.0;
static double t_evolve_iso = 10000.0;

/* Saved equilibrated profile */
static double *eq_phi[3], *eq_vel[3];
static int    eq_Nx = 0;
static double eq_dx = 0;
static double eq_xmax = 0;
static double eq_mass_osc = 0;
static double eq_omega_osc = 0;

/* ===================================================================
 *  Potential V(P) and force -dV/dphi_a
 * =================================================================== */

/* Confining + quartic: V = -sigma*sqrt(P^2+eps^2) + (kc/2)*P^4 */
static double V_pot(double P)
{
    double P2 = P * P;
    return -sigma_c * sqrt(P2 + eps_reg * eps_reg) + 0.5 * kc * P2 * P2;
}

/* Force: -dV/dphi_a */
static inline double force_pot(double p1, double p2, double p3, int a)
{
    double P = p1 * p2 * p3;
    double dP;
    switch (a) {
        case 0: dP = p2 * p3; break;
        case 1: dP = p1 * p3; break;
        case 2: dP = p1 * p2; break;
        default: dP = 0.0;
    }
    double P2 = P * P;
    double sq = sqrt(P2 + eps_reg * eps_reg);
    return sigma_c * P * dP / sq - 2.0 * kc * P2 * P * dP;
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
    e += V_pot(P);
    return e;
}

/* Energy density (periodic, for chain) */
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
    e += V_pot(P);
    return e;
}

/* ===================================================================
 *  Phase 1: Equilibrate single oscillon
 * =================================================================== */
static void run_equilibrate(void)
{
    printf("===== Phase 1: Single oscillon equilibration =====\n");

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

    printf("  sigma_c=%.3f kc=%.4f lambda=%.4f mass=%.4f\n",
           sigma_c, kc, lambda, mass);
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

    /* Compute acceleration with pairwise coupling */
    #define COMPUTE_ACC_EQ() do { \
        for (int a = 0; a < 3; a++) { \
            int b = (a+1)%3, c_idx = (a+2)%3; \
            acc[a][0] = acc[a][1] = acc[a][Nx-2] = acc[a][Nx-1] = 0; \
            for (int ii = 1; ii < Nx - 1; ii++) { \
                double lapl = (phi[a][ii+1] - 2.0*phi[a][ii] + phi[a][ii-1]) / dx2; \
                double fp = force_pot(phi[0][ii], phi[1][ii], phi[2][ii], a); \
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

    /* Breathing frequency from DFT of second half */
    eq_omega_osc = 0;
    {
        int dft_start = n_dft / 2;
        if (n_dft - dft_start > 100) {
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

    /* Write equilibrated profile */
    {
        char ppath[600];
        snprintf(ppath, sizeof(ppath), "%s/locality_profile.tsv", outdir);
        FILE *fp = fopen(ppath, "w");
        if (fp) {
            fprintf(fp, "x\tphi\tvel\n");
            for (int i = 0; i < Nx; i++) {
                double x = -xmax + i * dx;
                fprintf(fp, "%.6f\t%.8e\t%.8e\n", x, eq_phi[0][i], eq_vel[0][i]);
            }
            fclose(fp);
            printf("  Wrote %s\n", ppath);
        }
    }

    for (int a = 0; a < 3; a++) { free(phi[a]); free(vel[a]); free(acc[a]); }
    free(damp); free(phi0_hist); free(t_hist);
    printf("  Phase 1 complete.\n\n");
}

/* ===================================================================
 *  Initialize chain on periodic domain [0, L)
 * =================================================================== */
static void init_chain_periodic(double *phi[3], double *vel[3],
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

/* Initialize two oscillons on non-periodic domain [-xmax, xmax] */
static void init_two_oscillons(double *phi[3], double *vel[3],
                               int Nx, double dx, double xmax,
                               double x1, double x2)
{
    for (int a = 0; a < 3; a++) {
        memset(phi[a], 0, Nx * sizeof(double));
        memset(vel[a], 0, Nx * sizeof(double));
    }
    double centers[2] = {x1, x2};
    for (int c = 0; c < 2; c++) {
        for (int a = 0; a < 3; a++) {
            for (int i = 0; i < Nx; i++) {
                double x = -xmax + i * dx;
                double x_rel = x - centers[c];
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
 *  (periodic version)
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
            double new_x = xc + num / den;
            while (new_x < 0) new_x += L;
            while (new_x >= L) new_x -= L;
            x_osc[n] = new_x;
        }
    }
}

/* Track positions for two isolated oscillons (non-periodic) */
static void find_positions_isolated(double *phi0, double *phi1, double *phi2,
                                    double *v0, double *v1, double *v2,
                                    int Nx, double dx, double xmax,
                                    double *x_osc, double midpoint)
{
    /* oscillon 0: left half, oscillon 1: right half */
    for (int osc = 0; osc < 2; osc++) {
        double xc = x_osc[osc];
        double num = 0, den = 0;
        for (int i = 1; i < Nx - 1; i++) {
            double x = -xmax + i * dx;
            /* Only look in half corresponding to this oscillon */
            if (osc == 0 && x > midpoint) continue;
            if (osc == 1 && x < midpoint) continue;
            double rel = x - xc;
            if (fabs(rel) < 20.0) {  /* within 20 code units */
                double e = energy_density_np(phi0, phi1, phi2, v0, v1, v2, i, Nx, dx);
                if (e > 0) {
                    num += rel * e;
                    den += e;
                }
            }
        }
        if (den > 1e-30)
            x_osc[osc] = xc + num / den;
    }
}

/* ===================================================================
 *  Measure antisymmetric content at each oscillon
 * =================================================================== */
static void measure_asym_periodic(double *phi[3], int Nx, double dx, double L,
                                  double *x_osc, int N,
                                  double *S_amp, double *A12_amp, double *A13_amp)
{
    double core_r = 2.0 * sig;
    for (int n = 0; n < N; n++) {
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

static void measure_asym_isolated(double *phi[3], int Nx, double dx, double xmax,
                                  double *x_osc, int N,
                                  double *S_amp, double *A12_amp, double *A13_amp)
{
    double core_r = 2.0 * sig;
    for (int n = 0; n < N; n++) {
        double xc = x_osc[n];
        double s_num = 0, a12_num = 0, a13_num = 0;
        double w_sum = 0;
        for (int i = 0; i < Nx; i++) {
            double x = -xmax + i * dx;
            double rel = x - xc;
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
 *  Phase 2+3: Chain test
 * =================================================================== */
static void run_chain_test(void)
{
    printf("\n========== Phase 2+3: 8-Oscillon Chain Test ==========\n");

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
    printf("  N_osc=%d d=%.1f pert_osc=%d eps_pert=%.6f\n",
           N_osc, d_space, pert_osc, eps_pert);

    /* Oscillon positions */
    double *x_osc = malloc(N_osc * sizeof(double));
    double *x_osc_init = malloc(N_osc * sizeof(double));
    for (int n = 0; n < N_osc; n++) {
        x_osc[n] = n * d_space + d_space / 2.0;
        x_osc_init[n] = x_osc[n];
    }

    double *phi[3], *vel[3], *acc[3];
    for (int a = 0; a < 3; a++) {
        phi[a] = calloc(Nx, sizeof(double));
        vel[a] = calloc(Nx, sizeof(double));
        acc[a] = calloc(Nx, sizeof(double));
    }

    init_chain_periodic(phi, vel, Nx, dx, L, x_osc, N_osc);

    double E_init = 0;
    for (int i = 0; i < Nx; i++)
        E_init += energy_density_p(phi[0], phi[1], phi[2],
                                    vel[0], vel[1], vel[2], i, Nx, dx) * dx;
    printf("  Initial energy: %.4f (expected: %d * %.4f = %.4f)\n",
           E_init, N_osc, eq_mass_osc, N_osc * eq_mass_osc);

    /* Periodic Laplacian + pairwise coupling */
    #define COMPUTE_ACC_P() do { \
        for (int a = 0; a < 3; a++) { \
            int b = (a+1)%3, c_idx = (a+2)%3; \
            for (int ii = 0; ii < Nx; ii++) { \
                int ip1 = (ii + 1) % Nx; \
                int im1 = (ii - 1 + Nx) % Nx; \
                double lapl = (phi[a][ip1] - 2.0*phi[a][ii] + phi[a][im1]) / dx2; \
                double fp = force_pot(phi[0][ii], phi[1][ii], phi[2][ii], a); \
                acc[a][ii] = lapl - m2*phi[a][ii] \
                           - lambda*(phi[b][ii] + phi[c_idx][ii]) + fp; \
            } \
        } \
    } while(0)

    COMPUTE_ACC_P();

    /* ===== Pre-equilibrate chain with damping ===== */
    {
        int Nt_pre = (int)(t_pre_equil / dt) + 1;
        int print_pre = Nt_pre / 5;
        if (print_pre < 1) print_pre = 1;

        printf("  Pre-equilibrating (gamma=%.4f, t=%.0f)...\n",
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

            /* Position-dependent damping */
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

    /* ===== Conservative settle ===== */
    printf("  Conservative settle for t=%.0f...\n", t_settle);
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
                find_positions_periodic(phi[0], phi[1], phi[2],
                                        vel[0], vel[1], vel[2],
                                        Nx, dx, L, x_osc, N_osc);
                double max_asym = 0;
                double S_a[8], A12_a[8], A13_a[8];
                measure_asym_periodic(phi, Nx, dx, L, x_osc, N_osc,
                                      S_a, A12_a, A13_a);
                for (int n = 0; n < N_osc; n++) {
                    double av = fabs(A12_a[n]);
                    if (av > max_asym) max_asym = av;
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

    /* ===== Apply antisymmetric perturbation at oscillon #pert_osc ===== */
    printf("\n  Applying antisymmetric perturbation at oscillon #%d (eps=%.6f)\n",
           pert_osc, eps_pert);
    find_positions_periodic(phi[0], phi[1], phi[2],
                            vel[0], vel[1], vel[2],
                            Nx, dx, L, x_osc, N_osc);

    /* Record reference positions */
    double *x_osc_ref = malloc(N_osc * sizeof(double));
    for (int n = 0; n < N_osc; n++)
        x_osc_ref[n] = x_osc[n];

    /* Measure pre-perturbation baseline */
    double S_pre[8], A12_pre[8], A13_pre[8];
    measure_asym_periodic(phi, Nx, dx, L, x_osc, N_osc, S_pre, A12_pre, A13_pre);
    printf("  Pre-perturbation A12 at each oscillon:\n");
    for (int n = 0; n < N_osc; n++)
        printf("    osc %d: S=%.6f  A12=%.6e  A13=%.6e\n",
               n, S_pre[n], A12_pre[n], A13_pre[n]);

    double E_before = 0;
    for (int i = 0; i < Nx; i++)
        E_before += energy_density_p(phi[0], phi[1], phi[2],
                                      vel[0], vel[1], vel[2], i, Nx, dx) * dx;
    printf("  E_before = %.6f\n", E_before);

    /* Apply: phi_1 += +eps, phi_2 -= eps at oscillon #pert_osc */
    {
        double xc = x_osc[pert_osc];
        for (int i = 0; i < Nx; i++) {
            double x = i * dx;
            double rel = x - xc;
            while (rel > L / 2) rel -= L;
            while (rel < -L / 2) rel += L;
            double g = eps_pert * exp(-rel * rel / (2.0 * pert_sigma * pert_sigma));
            phi[0][i] += g;
            phi[1][i] -= g;
        }
    }
    COMPUTE_ACC_P();

    double E_after = 0;
    for (int i = 0; i < Nx; i++)
        E_after += energy_density_p(phi[0], phi[1], phi[2],
                                     vel[0], vel[1], vel[2], i, Nx, dx) * dx;
    printf("  E_after = %.6f  (Delta E = %.6e)\n", E_after, E_after - E_before);

    /* ===== Post-perturbation evolution ===== */
    printf("\n  Post-perturbation evolution for t=%.0f...\n", t_evolve);

    int Nt_evolve = (int)(t_evolve / dt) + 1;
    int diag_every_step = (int)(diag_dt / dt);
    if (diag_every_step < 1) diag_every_step = 1;
    int max_diag = Nt_evolve / diag_every_step + 10;
    int print_every = Nt_evolve / 40;
    if (print_every < 1) print_every = 1;

    /* Diagnostic arrays */
    double *d_t      = malloc(max_diag * sizeof(double));
    double *d_E      = malloc(max_diag * sizeof(double));
    double **d_A12   = malloc(N_osc * sizeof(double*));
    double **d_A13   = malloc(N_osc * sizeof(double*));
    double **d_S     = malloc(N_osc * sizeof(double*));
    double **d_xpos  = malloc(N_osc * sizeof(double*));
    double **d_disp  = malloc(N_osc * sizeof(double*));  /* displacement from ref */
    for (int n = 0; n < N_osc; n++) {
        d_A12[n]  = malloc(max_diag * sizeof(double));
        d_A13[n]  = malloc(max_diag * sizeof(double));
        d_S[n]    = malloc(max_diag * sizeof(double));
        d_xpos[n] = malloc(max_diag * sizeof(double));
        d_disp[n] = malloc(max_diag * sizeof(double));
    }
    int n_diag = 0;

    for (int step = 0; step <= Nt_evolve; step++) {
        double t = step * dt;

        if (step % diag_every_step == 0 && n_diag < max_diag) {
            if (step > 0)
                find_positions_periodic(phi[0], phi[1], phi[2],
                                        vel[0], vel[1], vel[2],
                                        Nx, dx, L, x_osc, N_osc);

            double Etot = 0;
            for (int i = 0; i < Nx; i++)
                Etot += energy_density_p(phi[0], phi[1], phi[2],
                                          vel[0], vel[1], vel[2], i, Nx, dx) * dx;

            double SA[8], A12A[8], A13A[8];
            measure_asym_periodic(phi, Nx, dx, L, x_osc, N_osc, SA, A12A, A13A);

            d_t[n_diag] = t;
            d_E[n_diag] = Etot;
            for (int n = 0; n < N_osc; n++) {
                d_S[n][n_diag]    = SA[n];
                d_A12[n][n_diag]  = A12A[n];
                d_A13[n][n_diag]  = A13A[n];
                d_xpos[n][n_diag] = x_osc[n];

                /* Displacement relative to reference (periodic) */
                double disp = x_osc[n] - x_osc_ref[n];
                while (disp > L / 2) disp -= L;
                while (disp < -L / 2) disp += L;
                d_disp[n][n_diag] = disp;
            }

            if (step % print_every == 0) {
                double max_a12 = 0;
                int max_n = 0;
                for (int n = 0; n < N_osc; n++) {
                    if (fabs(A12A[n]) > max_a12) {
                        max_a12 = fabs(A12A[n]);
                        max_n = n;
                    }
                }
                double max_disp_v = 0;
                for (int n = 0; n < N_osc; n++) {
                    double dd = fabs(d_disp[n][n_diag]);
                    if (dd > max_disp_v) max_disp_v = dd;
                }
                printf("  t=%7.1f  E=%.4f  max|A12|=%.6e at osc %d  max|disp|=%.6e\n",
                       t, Etot, max_a12, max_n, max_disp_v);
            }
            n_diag++;
        }

        if (step == Nt_evolve) break;

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
    #undef COMPUTE_ACC_P

    printf("\n  Evolution complete. %d diagnostic samples.\n", n_diag);

    /* ===== Write output files ===== */
    printf("\n  Writing chain output files...\n");

    /* Antisymmetric content time series */
    {
        char path[600];
        snprintf(path, sizeof(path), "%s/chain_asym_ts.tsv", outdir);
        FILE *fp = fopen(path, "w");
        if (fp) {
            fprintf(fp, "time");
            for (int n = 0; n < N_osc; n++) fprintf(fp, "\tA12_%d", n);
            for (int n = 0; n < N_osc; n++) fprintf(fp, "\tA13_%d", n);
            for (int n = 0; n < N_osc; n++) fprintf(fp, "\tS_%d", n);
            fprintf(fp, "\n");
            for (int j = 0; j < n_diag; j++) {
                fprintf(fp, "%.4f", d_t[j]);
                for (int n = 0; n < N_osc; n++)
                    fprintf(fp, "\t%.8e", d_A12[n][j]);
                for (int n = 0; n < N_osc; n++)
                    fprintf(fp, "\t%.8e", d_A13[n][j]);
                for (int n = 0; n < N_osc; n++)
                    fprintf(fp, "\t%.8e", d_S[n][j]);
                fprintf(fp, "\n");
            }
            fclose(fp);
            printf("  Wrote %s (%d rows)\n", path, n_diag);
        }
    }

    /* Displacement time series (KEY output) */
    {
        char path[600];
        snprintf(path, sizeof(path), "%s/chain_displacement_ts.tsv", outdir);
        FILE *fp = fopen(path, "w");
        if (fp) {
            fprintf(fp, "time");
            for (int n = 0; n < N_osc; n++) fprintf(fp, "\tdisp_%d", n);
            fprintf(fp, "\n");
            for (int j = 0; j < n_diag; j++) {
                fprintf(fp, "%.4f", d_t[j]);
                for (int n = 0; n < N_osc; n++)
                    fprintf(fp, "\t%.8e", d_disp[n][j]);
                fprintf(fp, "\n");
            }
            fclose(fp);
            printf("  Wrote %s\n", path);
        }
    }

    /* Energy time series */
    {
        char path[600];
        snprintf(path, sizeof(path), "%s/chain_energy_ts.tsv", outdir);
        FILE *fp = fopen(path, "w");
        if (fp) {
            fprintf(fp, "time\tE_total\n");
            for (int j = 0; j < n_diag; j++)
                fprintf(fp, "%.4f\t%.8e\n", d_t[j], d_E[j]);
            fclose(fp);
            printf("  Wrote %s\n", path);
        }
    }

    /* ===== Analysis: force detection ===== */
    printf("\n  ===== FORCE ANALYSIS =====\n");

    /* For each oscillon, find the first time |displacement| exceeds threshold */
    double disp_threshold = 1e-6;
    printf("  Displacement threshold: %.1e\n", disp_threshold);
    printf("  Reference: perturbation at oscillon #%d\n\n", pert_osc);

    /* Compute distance from perturbed oscillon (periodic, in lattice spacings) */
    for (int n = 0; n < N_osc; n++) {
        int dist = n - pert_osc;
        if (dist < 0) dist += N_osc;
        if (dist > N_osc / 2) dist = N_osc - dist;
        double phys_dist = dist * d_space;

        /* Find first time displacement exceeds threshold */
        double t_first = -1;
        double max_disp_this = 0;
        double t_at_max = 0;
        double avg_disp_late = 0;
        int n_late = 0;

        for (int j = 0; j < n_diag; j++) {
            double dd = fabs(d_disp[n][j]);
            if (dd > max_disp_this) {
                max_disp_this = dd;
                t_at_max = d_t[j];
            }
            if (t_first < 0 && dd > disp_threshold)
                t_first = d_t[j];
            /* Late-time average (last 20%) */
            if (j > 0.8 * n_diag) {
                avg_disp_late += d_disp[n][j];
                n_late++;
            }
        }
        if (n_late > 0) avg_disp_late /= n_late;

        /* Also find first time A12 exceeds threshold */
        double a12_threshold = 1e-6;
        double t_a12_first = -1;
        double max_a12_this = 0;
        for (int j = 0; j < n_diag; j++) {
            double aa = fabs(d_A12[n][j]);
            if (aa > max_a12_this) max_a12_this = aa;
            if (t_a12_first < 0 && aa > a12_threshold)
                t_a12_first = d_t[j];
        }

        /* Direction: is avg_disp toward or away from pert_osc? */
        /* On periodic ring, "toward" means displacement decreases the distance */
        const char *direction;
        if (fabs(avg_disp_late) < 1e-10)
            direction = "NONE";
        else {
            /* Compute which direction is toward pert_osc */
            double dx_toward = x_osc_ref[pert_osc] - x_osc_ref[n];
            while (dx_toward > L / 2) dx_toward -= L;
            while (dx_toward < -L / 2) dx_toward += L;
            if (dx_toward * avg_disp_late > 0)
                direction = "ATTRACTIVE";
            else
                direction = "REPULSIVE";
        }

        printf("  Osc #%d: dist=%d (%.0f)  t_A12_arrive=%.1f  t_force=%.1f  "
               "max|disp|=%.3e  avg_disp=%.3e  %s\n",
               n, dist, phys_dist, t_a12_first, t_first,
               max_disp_this, avg_disp_late, direction);

        if (dist > 0 && t_first > 0 && phys_dist > 0) {
            double c_signal = phys_dist / t_first;
            printf("         signal speed ~ %.4f (c=1)\n", c_signal);
        }
    }

    /* Asymmetric mode speed from arrival times */
    printf("\n  ===== ARRIVAL TIME vs DISTANCE =====\n");
    {
        char path[600];
        snprintf(path, sizeof(path), "%s/chain_arrival.tsv", outdir);
        FILE *fp = fopen(path, "w");
        if (fp) {
            fprintf(fp, "osc\tdist_lattice\tdist_phys\tt_A12_arrive\tt_force_arrive\t"
                        "max_A12\tmax_disp\tavg_disp_late\n");
            for (int n = 0; n < N_osc; n++) {
                int dist = n - pert_osc;
                if (dist < 0) dist += N_osc;
                if (dist > N_osc / 2) dist = N_osc - dist;
                double phys_dist = dist * d_space;

                double t_a12 = -1, t_disp = -1;
                double max_a12_this = 0, max_disp_this = 0;
                double avg_late = 0;
                int n_late = 0;

                for (int j = 0; j < n_diag; j++) {
                    double aa = fabs(d_A12[n][j]);
                    double dd = fabs(d_disp[n][j]);
                    if (aa > max_a12_this) max_a12_this = aa;
                    if (dd > max_disp_this) max_disp_this = dd;
                    if (t_a12 < 0 && aa > 1e-6) t_a12 = d_t[j];
                    if (t_disp < 0 && dd > 1e-6) t_disp = d_t[j];
                    if (j > 0.8 * n_diag) { avg_late += d_disp[n][j]; n_late++; }
                }
                if (n_late > 0) avg_late /= n_late;

                fprintf(fp, "%d\t%d\t%.1f\t%.2f\t%.2f\t%.6e\t%.6e\t%.6e\n",
                        n, dist, phys_dist, t_a12, t_disp,
                        max_a12_this, max_disp_this, avg_late);
            }
            fclose(fp);
            printf("  Wrote %s\n", path);
        }
    }

    /* Normal mode decomposition (antisymmetric spatial modes) */
    printf("\n  ===== ANTISYMMETRIC MODE SPECTRUM =====\n");
    {
        double T_an = d_t[n_diag-1] - d_t[0];
        double dt_eff = T_an / (n_diag - 1);

        /* Spatial Fourier of A12 */
        int n_omega = n_diag / 2;
        if (n_omega > 4000) n_omega = 4000;
        double omega_max = 2.0;

        char spath[600];
        snprintf(spath, sizeof(spath), "%s/chain_asym_spectrum.tsv", outdir);
        FILE *fsp = fopen(spath, "w");
        if (fsp) fprintf(fsp, "q\tk\tomega_peak\tpower\tc_phase\n");

        for (int q = 0; q <= N_osc / 2; q++) {
            double k_q = 2.0 * M_PI * q / (N_osc * d_space);

            /* Spatial DFT of A12 at each time step */
            double *AQ_re = calloc(n_diag, sizeof(double));
            double *AQ_im = calloc(n_diag, sizeof(double));
            for (int j = 0; j < n_diag; j++) {
                double re = 0, im = 0;
                for (int n = 0; n < N_osc; n++) {
                    double phase = -2.0 * M_PI * q * n / N_osc;
                    re += d_A12[n][j] * cos(phase);
                    im += d_A12[n][j] * sin(phase);
                }
                AQ_re[j] = re / sqrt((double)N_osc);
                AQ_im[j] = im / sqrt((double)N_osc);
            }

            /* Temporal DFT */
            double best_pow = 0, best_omega = 0;
            for (int l = 1; l <= n_omega; l++) {
                double omega = omega_max * l / n_omega;
                double dft_re = 0, dft_im = 0;
                for (int j = 0; j < n_diag; j++) {
                    double tj = j * dt_eff;
                    dft_re += (AQ_re[j] * cos(omega * tj) + AQ_im[j] * sin(omega * tj)) * dt_eff;
                    dft_im += (-AQ_re[j] * sin(omega * tj) + AQ_im[j] * cos(omega * tj)) * dt_eff;
                }
                double power = (dft_re * dft_re + dft_im * dft_im) / (T_an * T_an);
                if (power > best_pow) { best_pow = power; best_omega = omega; }
            }

            double cs = (k_q > 1e-10) ? best_omega / k_q : 0;
            printf("  q=%d: k=%.5f  omega=%.5f  power=%.4e  c_phase=%.4f\n",
                   q, k_q, best_omega, best_pow, cs);
            if (fsp)
                fprintf(fsp, "%d\t%.6f\t%.6f\t%.8e\t%.6f\n",
                        q, k_q, best_omega, best_pow, cs);

            free(AQ_re); free(AQ_im);
        }
        if (fsp) { fclose(fsp); printf("  Wrote %s\n", spath); }
    }

    /* Energy conservation check */
    {
        double dE = (d_E[n_diag-1] - d_E[0]) / (fabs(d_E[0]) + 1e-30);
        printf("\n  Energy conservation: dE/E = %.3e\n", dE);
    }

    /* Cleanup */
    free(d_t); free(d_E); free(x_osc); free(x_osc_init); free(x_osc_ref);
    for (int n = 0; n < N_osc; n++) {
        free(d_A12[n]); free(d_A13[n]); free(d_S[n]);
        free(d_xpos[n]); free(d_disp[n]);
    }
    free(d_A12); free(d_A13); free(d_S); free(d_xpos); free(d_disp);
    for (int a = 0; a < 3; a++) { free(phi[a]); free(vel[a]); free(acc[a]); }

    printf("  Phase 2+3 complete.\n\n");
}

/* ===================================================================
 *  Phase 4: Two isolated oscillons at D=40
 * =================================================================== */
static void run_isolated_test(void)
{
    printf("\n========== Phase 4: Two Isolated Oscillons (D=%.0f) ==========\n", D_iso);

    int Nx = Nx_iso;
    double xmax = xmax_iso;
    double dx  = 2.0 * xmax / (Nx - 1);
    double dx2 = dx * dx;
    double m2  = mass * mass;

    double m2_eff = m2 + 2.0 * lambda;
    if (m2_eff < m2) m2_eff = m2;
    double kmax_c = M_PI / dx;
    double dt = 0.8 * 2.0 / sqrt(kmax_c * kmax_c + m2_eff);

    printf("  Nx=%d xmax=%.1f dx=%.5f dt=%.6f\n", Nx, xmax, dx, dt);

    /* Place oscillons at -D/2 and +D/2 */
    double x1 = -D_iso / 2.0;
    double x2 = +D_iso / 2.0;
    printf("  Oscillon centers: x1=%.1f  x2=%.1f  D=%.1f\n", x1, x2, D_iso);

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

    init_two_oscillons(phi, vel, Nx, dx, xmax, x1, x2);

    /* Non-periodic acceleration */
    #define COMPUTE_ACC_NP() do { \
        for (int a = 0; a < 3; a++) { \
            int b = (a+1)%3, c_idx = (a+2)%3; \
            acc[a][0] = acc[a][Nx-1] = 0; \
            for (int ii = 1; ii < Nx - 1; ii++) { \
                double lapl = (phi[a][ii+1] - 2.0*phi[a][ii] + phi[a][ii-1]) / dx2; \
                double fp = force_pot(phi[0][ii], phi[1][ii], phi[2][ii], a); \
                acc[a][ii] = lapl - m2*phi[a][ii] \
                           - lambda*(phi[b][ii] + phi[c_idx][ii]) + fp; \
            } \
        } \
    } while(0)

    COMPUTE_ACC_NP();

    double E_init = 0;
    for (int i = 1; i < Nx - 1; i++)
        E_init += energy_density_np(phi[0], phi[1], phi[2],
                                     vel[0], vel[1], vel[2], i, Nx, dx) * dx;
    printf("  Initial energy: %.4f (expected: 2 * %.4f = %.4f)\n",
           E_init, eq_mass_osc, 2 * eq_mass_osc);

    /* Pre-equilibrate with damping between oscillons */
    {
        int Nt_pre = (int)(t_pre_equil / dt) + 1;
        printf("  Pre-equilibrating (t=%.0f)...\n", t_pre_equil);

        for (int step = 0; step < Nt_pre; step++) {
            for (int a = 0; a < 3; a++)
                for (int i = 1; i < Nx - 1; i++)
                    vel[a][i] += 0.5 * dt * acc[a][i];
            for (int a = 0; a < 3; a++)
                for (int i = 1; i < Nx - 1; i++)
                    phi[a][i] += dt * vel[a][i];
            COMPUTE_ACC_NP();
            for (int a = 0; a < 3; a++)
                for (int i = 1; i < Nx - 1; i++)
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

            /* Absorbing boundary */
            for (int a = 0; a < 3; a++)
                for (int i = 0; i < Nx; i++) {
                    vel[a][i] *= damp[i];
                    phi[a][i] *= damp[i];
                }
        }
        COMPUTE_ACC_NP();
    }

    /* Conservative settle */
    printf("  Conservative settle for t=%.0f...\n", t_settle_iso);
    {
        int Nt_settle = (int)(t_settle_iso / dt) + 1;
        int print_every = Nt_settle / 5;
        if (print_every < 1) print_every = 1;

        for (int step = 0; step <= Nt_settle; step++) {
            if (step % print_every == 0) {
                double Etot = 0;
                for (int i = 1; i < Nx - 1; i++)
                    Etot += energy_density_np(phi[0], phi[1], phi[2],
                                              vel[0], vel[1], vel[2], i, Nx, dx) * dx;
                printf("  settle t=%7.1f  E=%.4f\n", step * dt, Etot);
            }
            if (step == Nt_settle) break;

            for (int a = 0; a < 3; a++)
                for (int i = 1; i < Nx - 1; i++)
                    vel[a][i] += 0.5 * dt * acc[a][i];
            for (int a = 0; a < 3; a++)
                for (int i = 1; i < Nx - 1; i++)
                    phi[a][i] += dt * vel[a][i];
            COMPUTE_ACC_NP();
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
    }

    /* Apply antisymmetric perturbation at oscillon #1 (x1) */
    double x_osc[2] = {x1, x2};
    find_positions_isolated(phi[0], phi[1], phi[2],
                            vel[0], vel[1], vel[2],
                            Nx, dx, xmax, x_osc, 0.0);

    double x_osc_ref[2] = {x_osc[0], x_osc[1]};

    printf("  Positions before perturbation: x1=%.4f  x2=%.4f  D=%.4f\n",
           x_osc[0], x_osc[1], x_osc[1] - x_osc[0]);

    /* Measure pre-perturbation */
    double S_pre[2], A12_pre[2], A13_pre[2];
    measure_asym_isolated(phi, Nx, dx, xmax, x_osc, 2, S_pre, A12_pre, A13_pre);
    printf("  Pre-perturbation: A12[0]=%.6e  A12[1]=%.6e\n", A12_pre[0], A12_pre[1]);

    double E_before = 0;
    for (int i = 1; i < Nx - 1; i++)
        E_before += energy_density_np(phi[0], phi[1], phi[2],
                                       vel[0], vel[1], vel[2], i, Nx, dx) * dx;

    /* Apply: phi_1 += +eps, phi_2 -= eps at oscillon #1 (x1) */
    {
        double xc = x_osc[0];
        for (int i = 0; i < Nx; i++) {
            double x = -xmax + i * dx;
            double rel = x - xc;
            double g = eps_pert * exp(-rel * rel / (2.0 * pert_sigma * pert_sigma));
            phi[0][i] += g;
            phi[1][i] -= g;
        }
    }
    COMPUTE_ACC_NP();

    double E_after = 0;
    for (int i = 1; i < Nx - 1; i++)
        E_after += energy_density_np(phi[0], phi[1], phi[2],
                                      vel[0], vel[1], vel[2], i, Nx, dx) * dx;
    printf("  E_before=%.6f  E_after=%.6f  Delta=%.6e\n",
           E_before, E_after, E_after - E_before);

    /* Evolve */
    printf("  Evolving for t=%.0f...\n", t_evolve_iso);

    int Nt_ev = (int)(t_evolve_iso / dt) + 1;
    int diag_ev_step = (int)(diag_dt / dt);
    if (diag_ev_step < 1) diag_ev_step = 1;
    int max_diag = Nt_ev / diag_ev_step + 10;
    int print_every = Nt_ev / 40;
    if (print_every < 1) print_every = 1;

    double *d_t    = malloc(max_diag * sizeof(double));
    double *d_E    = malloc(max_diag * sizeof(double));
    double *d_x1   = malloc(max_diag * sizeof(double));
    double *d_x2   = malloc(max_diag * sizeof(double));
    double *d_d1   = malloc(max_diag * sizeof(double));  /* disp of osc 1 */
    double *d_d2   = malloc(max_diag * sizeof(double));  /* disp of osc 2 */
    double *d_a12_1 = malloc(max_diag * sizeof(double)); /* A12 at osc 1 */
    double *d_a12_2 = malloc(max_diag * sizeof(double)); /* A12 at osc 2 */
    double *d_sep  = malloc(max_diag * sizeof(double));  /* separation */
    int n_diag = 0;

    for (int step = 0; step <= Nt_ev; step++) {
        double t = step * dt;

        if (step % diag_ev_step == 0 && n_diag < max_diag) {
            if (step > 0)
                find_positions_isolated(phi[0], phi[1], phi[2],
                                        vel[0], vel[1], vel[2],
                                        Nx, dx, xmax, x_osc, 0.0);

            double Etot = 0;
            for (int i = 1; i < Nx - 1; i++)
                Etot += energy_density_np(phi[0], phi[1], phi[2],
                                          vel[0], vel[1], vel[2], i, Nx, dx) * dx;

            double SA[2], A12A[2], A13A[2];
            measure_asym_isolated(phi, Nx, dx, xmax, x_osc, 2, SA, A12A, A13A);

            d_t[n_diag] = t;
            d_E[n_diag] = Etot;
            d_x1[n_diag] = x_osc[0];
            d_x2[n_diag] = x_osc[1];
            d_d1[n_diag] = x_osc[0] - x_osc_ref[0];
            d_d2[n_diag] = x_osc[1] - x_osc_ref[1];
            d_a12_1[n_diag] = A12A[0];
            d_a12_2[n_diag] = A12A[1];
            d_sep[n_diag] = x_osc[1] - x_osc[0];

            if (step % print_every == 0) {
                printf("  t=%7.1f  E=%.4f  sep=%.4f  d1=%+.3e  d2=%+.3e  "
                       "A12_1=%.3e  A12_2=%.3e\n",
                       t, Etot, d_sep[n_diag],
                       d_d1[n_diag], d_d2[n_diag],
                       d_a12_1[n_diag], d_a12_2[n_diag]);
            }
            n_diag++;
        }

        if (step == Nt_ev) break;

        for (int a = 0; a < 3; a++)
            for (int i = 1; i < Nx - 1; i++)
                vel[a][i] += 0.5 * dt * acc[a][i];
        for (int a = 0; a < 3; a++)
            for (int i = 1; i < Nx - 1; i++)
                phi[a][i] += dt * vel[a][i];
        COMPUTE_ACC_NP();
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
    #undef COMPUTE_ACC_NP

    printf("\n  Evolution complete. %d diagnostic samples.\n", n_diag);

    /* Write output */
    {
        char path[600];
        snprintf(path, sizeof(path), "%s/isolated_ts.tsv", outdir);
        FILE *fp = fopen(path, "w");
        if (fp) {
            fprintf(fp, "time\tE\tx1\tx2\tdisp1\tdisp2\tA12_1\tA12_2\tsep\n");
            for (int j = 0; j < n_diag; j++) {
                fprintf(fp, "%.4f\t%.8e\t%.6f\t%.6f\t%.8e\t%.8e\t%.8e\t%.8e\t%.6f\n",
                        d_t[j], d_E[j], d_x1[j], d_x2[j],
                        d_d1[j], d_d2[j], d_a12_1[j], d_a12_2[j], d_sep[j]);
            }
            fclose(fp);
            printf("  Wrote %s\n", path);
        }
    }

    /* Analysis */
    printf("\n  ===== ISOLATED PAIR ANALYSIS =====\n");
    {
        /* Find when A12 arrives at oscillon 2 */
        double t_a12_arrive = -1;
        for (int j = 0; j < n_diag; j++) {
            if (fabs(d_a12_2[j]) > 1e-6) {
                t_a12_arrive = d_t[j];
                break;
            }
        }

        /* Find when displacement of oscillon 2 exceeds threshold */
        double t_force_arrive = -1;
        for (int j = 0; j < n_diag; j++) {
            if (fabs(d_d2[j]) > 1e-6) {
                t_force_arrive = d_t[j];
                break;
            }
        }

        /* Late-time displacement */
        double avg_d2 = 0;
        int n_late = 0;
        for (int j = (int)(0.8 * n_diag); j < n_diag; j++) {
            avg_d2 += d_d2[j];
            n_late++;
        }
        if (n_late > 0) avg_d2 /= n_late;

        /* Maximum displacement */
        double max_d2 = 0;
        for (int j = 0; j < n_diag; j++) {
            double dd = fabs(d_d2[j]);
            if (dd > max_d2) max_d2 = dd;
        }

        printf("  Distance: D = %.1f\n", D_iso);
        printf("  A12 arrival at osc 2: t = %.1f", t_a12_arrive);
        if (t_a12_arrive > 0)
            printf("  (signal speed = %.4f)", D_iso / t_a12_arrive);
        printf("\n");
        printf("  Force arrival at osc 2: t = %.1f", t_force_arrive);
        if (t_force_arrive > 0)
            printf("  (signal speed = %.4f)", D_iso / t_force_arrive);
        printf("\n");
        printf("  Max |displacement| of osc 2: %.6e\n", max_d2);
        printf("  Late-time avg displacement of osc 2: %+.6e\n", avg_d2);
        if (fabs(avg_d2) > 1e-10) {
            if (avg_d2 < 0) /* osc 2 moves left, toward osc 1 */
                printf("  Direction: ATTRACTIVE (osc 2 moves toward osc 1)\n");
            else
                printf("  Direction: REPULSIVE (osc 2 moves away from osc 1)\n");
        } else {
            printf("  Direction: NEGLIGIBLE\n");
        }

        /* Separation change */
        double sep_init = d_sep[0];
        double sep_final = d_sep[n_diag - 1];
        printf("  Separation: initial=%.4f  final=%.4f  change=%.6e\n",
               sep_init, sep_final, sep_final - sep_init);
    }

    /* Cleanup */
    free(d_t); free(d_E); free(d_x1); free(d_x2);
    free(d_d1); free(d_d2); free(d_a12_1); free(d_a12_2); free(d_sep);
    free(damp);
    for (int a = 0; a < 3; a++) { free(phi[a]); free(vel[a]); free(acc[a]); }

    printf("  Phase 4 complete.\n\n");
}

/* ===================================================================
 *  Main
 * =================================================================== */
int main(int argc, char **argv)
{
    for (int i = 1; i < argc - 1; i += 2) {
        if      (!strcmp(argv[i], "-mass"))     mass     = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-sigma_c"))  sigma_c  = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-kc"))       kc       = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-lambda"))   lambda   = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-A"))        A_init   = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-sig"))      sig      = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-o"))        strncpy(outdir, argv[i+1], sizeof(outdir)-1);
        else if (!strcmp(argv[i], "-N"))        N_osc    = atoi(argv[i+1]);
        else if (!strcmp(argv[i], "-d"))        d_space  = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-eps"))      eps_pert = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-pert_osc")) pert_osc = atoi(argv[i+1]);
        else if (!strcmp(argv[i], "-D"))        D_iso    = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-t_equil"))  t_equil  = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-t_settle")) t_settle = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-t_evolve")) t_evolve = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-diag_dt"))  diag_dt  = atof(argv[i+1]);
        else { fprintf(stderr, "Unknown flag: %s\n", argv[i]); exit(1); }
    }

    printf("=== Locality Option 2: Antisymmetric Mode as Causal Gravity ===\n");
    printf("  sigma_c=%.3f  kc=%.4f  lambda=%.4f  mass=%.4f\n",
           sigma_c, kc, lambda, mass);
    printf("  N_osc=%d  d=%.1f  eps_pert=%.6f  pert_osc=%d\n",
           N_osc, d_space, eps_pert, pert_osc);
    printf("  D_iso=%.1f  t_equil=%.0f  t_settle=%.0f  t_evolve=%.0f\n",
           D_iso, t_equil, t_settle, t_evolve);

    /* Phase 1: Equilibrate */
    run_equilibrate();

    /* Phase 2+3: Chain test */
    run_chain_test();

    /* Phase 4: Two isolated oscillons */
    run_isolated_test();

    /* Cleanup */
    for (int a = 0; a < 3; a++) {
        free(eq_phi[a]);
        free(eq_vel[a]);
    }

    printf("=== All phases complete. ===\n");
    return 0;
}
