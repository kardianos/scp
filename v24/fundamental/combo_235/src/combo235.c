/*
 * combo235.c -- Combo 2+3+5: Confined Lattice with Inertial Deformation
 *
 * Combines:
 *   Test C (Confine): V = -sigma*sqrt(P^2+eps^2) + (kc/2)*P^4
 *   Test E (Lattice): pairwise coupling lambda*(phi1*phi2 + phi2*phi3 + phi3*phi1)
 *   Combo 2+5: antisymmetric perturbation tracking
 *
 * Force on phi_a:
 *   acc = laplacian(phi_a) - m^2 phi_a
 *       + sigma * P * (dP/dphi_a) / sqrt(P^2 + eps^2)
 *       - 2*kc * P^3 * (dP/dphi_a)
 *       - lambda * (phi_b + phi_c)
 *
 * Phases:
 *   1: Equilibrate single oscillon (t=5000)
 *   2: Build 8-oscillon chain, evolve t=10000, phonon spectrum
 *   3: Antisymmetric perturbation at oscillon #4, track propagation
 *
 * Also runs standard potential (mu/2)P^2/(1+kappa*P^2) for comparison.
 *
 * Compile: gcc -O3 -Wall -o combo235 src/combo235.c -lm
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* --- Parameters --- */
static double mass     = 1.0;
static double A_init   = 1.0;
static double sig      = 3.0;
static char   outdir[512] = "data";

/* Confining potential parameters */
static double sigma_c  = 1.0;    /* string tension */
static double eps_reg   = 1e-6;   /* regularization */
static double kc       = 0.01;   /* quartic stabilizer */

/* Standard potential parameters (for comparison run) */
static double mu_std   = -20.0;
static double kappa_std = 20.0;

/* Pairwise coupling */
static double lambda   = 0.5;

/* Phase 1: equilibration */
static int    Nx_eq    = 4000;
static double xmax_eq  = 100.0;
static double t_equil  = 5000.0;

/* Phase 2: chain */
static int    N_osc    = 8;
static double d_space  = 16.0;
static double t_chain  = 10000.0;
static double t_discard = 2000.0;
static double diag_dt  = 10.0;
static int    pts_per_d = 320;
static double t_pre_equil = 500.0;
static double gamma_pre = 0.001;
static double A_kick   = 0.01;

/* Phase 3: antisymmetric perturbation */
static double eps_pert  = 0.01;
static int    pert_osc  = 4;
static double pert_sigma = 3.0;
static double t_settle  = 5000.0;
static double t_evolve  = 10000.0;

/* Potential mode: 0=confining+quartic, 1=standard saturating */
static int pot_mode = 0;

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
/* Standard saturating: V = (mu/2)*P^2/(1+kappa*P^2) */
static double V_pot(double P)
{
    double P2 = P * P;
    if (pot_mode == 0) {
        return -sigma_c * sqrt(P2 + eps_reg * eps_reg) + 0.5 * kc * P2 * P2;
    } else {
        return 0.5 * mu_std * P2 / (1.0 + kappa_std * P2);
    }
}

/*
 * Force: -dV/dphi_a (does NOT include pairwise coupling)
 * dP/dphi_a: for a=0, dP=phi2*phi3, etc.
 */
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
    if (pot_mode == 0) {
        /* Confining + quartic:
         * -dV/dphi_a = sigma * P * dP / sqrt(P^2 + eps^2) - 2*kc*P^3*dP */
        double sq = sqrt(P2 + eps_reg * eps_reg);
        return sigma_c * P * dP / sq - 2.0 * kc * P2 * P * dP;
    } else {
        /* Standard: -dV/dphi_a = -mu*P*dP/(1+kappa*P^2)^2 */
        double denom2 = (1.0 + kappa_std * P2) * (1.0 + kappa_std * P2);
        return -mu_std * P * dP / denom2;
    }
}

/* Simple LCG random in [0,1) */
static double rand01(unsigned int *state)
{
    *state = *state * 1103515245u + 12345u;
    return (double)((*state >> 16) & 0x7fff) / 32768.0;
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

    /* Pairwise coupling energy */
    e += lambda * (phi0[i]*phi1[i] + phi1[i]*phi2[i] + phi2[i]*phi0[i]);

    /* Potential energy */
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
static void run_equilibrate(const char *label)
{
    printf("===== Phase 1: Single oscillon equilibration [%s] (lambda=%.4f) =====\n",
           label, lambda);

    int Nx = Nx_eq;
    double xmax = xmax_eq;
    double dx  = 2.0 * xmax / (Nx - 1);
    double dx2 = dx * dx;
    double m2  = mass * mass;

    /* CFL: use heaviest mode mass for stability */
    double m2_eff = m2 + 2.0 * lambda;
    if (m2_eff < m2) m2_eff = m2;
    double kmax = M_PI / dx;
    double dt   = 0.8 * 2.0 / sqrt(kmax * kmax + m2_eff);
    int    Nt   = (int)(t_equil / dt) + 1;

    printf("  Nx=%d xmax=%.1f dx=%.5f dt=%.6f Nt=%d tfinal=%.0f\n",
           Nx, xmax, dx, dt, Nt, t_equil);
    if (pot_mode == 0)
        printf("  sigma=%.3f kc=%.4f eps=%.2e lambda=%.4f\n",
               sigma_c, kc, eps_reg, lambda);
    else
        printf("  mu=%.1f kappa=%.1f lambda=%.4f\n", mu_std, kappa_std, lambda);

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

            /* Write spectrum */
            char spath[600];
            snprintf(spath, sizeof(spath), "%s/combo235_%s_spectrum.tsv",
                     outdir, label);
            FILE *fspec = fopen(spath, "w");
            if (fspec) {
                fprintf(fspec, "omega\tpower\n");
                for (int k = 0; k < nf; k++) {
                    double omega = 3.0 * mass * k / nf;
                    double re2 = 0, im2 = 0;
                    for (int j = dft_start; j < n_dft; j++) {
                        double dtj2 = (j > dft_start) ?
                            (t_hist[j]-t_hist[j-1]) : (t_hist[dft_start+1]-t_hist[dft_start]);
                        re2 += phi0_hist[j] * cos(omega * t_hist[j]) * dtj2;
                        im2 += phi0_hist[j] * sin(omega * t_hist[j]) * dtj2;
                    }
                    double pw2 = (re2*re2 + im2*im2) / (T*T);
                    fprintf(fspec, "%.6f\t%.6e\n", omega, pw2);
                }
                fclose(fspec);
                printf("  Wrote %s\n", spath);
            }
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
        snprintf(ppath, sizeof(ppath), "%s/combo235_%s_profile.tsv", outdir, label);
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
    printf("  Phase 1 [%s] complete.\n\n", label);
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
static void find_positions(double *phi0, double *phi1, double *phi2,
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

/* ===================================================================
 *  Measure antisymmetric content at each oscillon
 * =================================================================== */
static void measure_asym(double *phi[3], int Nx, double dx, double L,
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

/* ===================================================================
 *  Run chain evolution (Phase 2 or 3)
 *  phase=2: build chain, phonon kick, evolve, spectrum
 *  phase=3: antisymmetric perturbation, evolve, track
 * =================================================================== */
static void run_chain(const char *label, int do_antisym,
                      /* output arrays for comparison */
                      double *out_omega_sym, double *out_omega_asym,
                      double *out_max_disp, double *out_E_consv,
                      double *out_asym_gap)
{
    printf("\n===== Phase %d: %s chain [%s] =====\n",
           do_antisym ? 3 : 2,
           do_antisym ? "Antisymmetric perturbation" : "Phonon spectrum",
           label);

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

    /* Oscillon positions */
    double *x_osc = malloc(N_osc * sizeof(double));
    double *x_osc_ref = malloc(N_osc * sizeof(double));

    for (int n = 0; n < N_osc; n++) {
        x_osc[n] = n * d_space + d_space / 2.0;
        x_osc_ref[n] = x_osc[n];
    }

    /* Allocate fields */
    double *phi[3], *vel[3], *acc[3];
    for (int a = 0; a < 3; a++) {
        phi[a] = calloc(Nx, sizeof(double));
        vel[a] = calloc(Nx, sizeof(double));
        acc[a] = calloc(Nx, sizeof(double));
    }

    init_chain(phi, vel, Nx, dx, L, x_osc, N_osc);

    double E_init = 0;
    for (int i = 0; i < Nx; i++)
        E_init += energy_density_p(phi[0], phi[1], phi[2],
                                    vel[0], vel[1], vel[2], i, Nx, dx) * dx;
    printf("  Initial energy: %.4f (expected: %d * %.4f = %.4f)\n",
           E_init, N_osc, eq_mass_osc, N_osc * eq_mass_osc);

    /* Compute acceleration with periodic BC + pairwise coupling */
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

    if (!do_antisym) {
        /* Phase 2: Apply phonon kick and evolve */
        printf("  Applying phonon kick (A=%.4f)...\n", A_kick);
        for (int i = 0; i < Nx; i++) {
            double x = i * dx;
            double kick = A_kick * sin(2.0 * M_PI * x / L);
            for (int a = 0; a < 3; a++)
                vel[a][i] += kick * phi[a][i];
        }
        COMPUTE_ACC_P();

        printf("  Evolving chain for t=%.0f...\n", t_chain);

        int Nt = (int)(t_chain / dt) + 1;
        int max_diag = (int)(t_chain / diag_dt) + 10;
        double *diag_t = malloc(max_diag * sizeof(double));
        double **diag_x = malloc(N_osc * sizeof(double*));
        for (int n = 0; n < N_osc; n++)
            diag_x[n] = malloc(max_diag * sizeof(double));
        double *diag_E = malloc(max_diag * sizeof(double));
        double *diag_maxdisp = malloc(max_diag * sizeof(double));
        int n_diag = 0;

        int diag_every = (int)(diag_dt / dt);
        if (diag_every < 1) diag_every = 1;
        int print_every = Nt / 20;
        if (print_every < 1) print_every = 1;
        int ref_set = 0;

        for (int step = 0; step <= Nt; step++) {
            double t = step * dt;

            if (step % diag_every == 0 && n_diag < max_diag) {
                if (step > 0)
                    find_positions(phi[0], phi[1], phi[2],
                                   vel[0], vel[1], vel[2],
                                   Nx, dx, L, x_osc, N_osc);

                if (!ref_set && t >= t_discard) {
                    for (int n = 0; n < N_osc; n++)
                        x_osc_ref[n] = x_osc[n];
                    ref_set = 1;
                }

                double Etot = 0;
                for (int i = 0; i < Nx; i++)
                    Etot += energy_density_p(phi[0], phi[1], phi[2],
                                              vel[0], vel[1], vel[2], i, Nx, dx) * dx;

                double md = 0;
                for (int n = 0; n < N_osc; n++) {
                    double disp = x_osc[n] - x_osc_ref[n];
                    while (disp > L / 2) disp -= L;
                    while (disp < -L / 2) disp += L;
                    if (fabs(disp) > md) md = fabs(disp);
                }

                diag_t[n_diag] = t;
                for (int n = 0; n < N_osc; n++)
                    diag_x[n][n_diag] = x_osc[n];
                diag_E[n_diag] = Etot;
                diag_maxdisp[n_diag] = md;
                n_diag++;

                if (step % print_every == 0)
                    printf("  t=%7.1f  E=%.4f  max_u=%.4f\n", t, Etot, md);
            }

            if (step == Nt) break;

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

        /* Write chain data */
        {
            char path[600];
            snprintf(path, sizeof(path), "%s/combo235_%s_chain_ts.tsv", outdir, label);
            FILE *fp = fopen(path, "w");
            if (fp) {
                fprintf(fp, "time\tE_total\tmax_displacement\n");
                for (int j = 0; j < n_diag; j++)
                    fprintf(fp, "%.4f\t%.6e\t%.6e\n",
                            diag_t[j], diag_E[j], diag_maxdisp[j]);
                fclose(fp);
                printf("  Wrote %s\n", path);
            }
        }

        /* ===== Phonon spectrum (normal mode decomposition) ===== */
        int j_start = 0;
        for (int j = 0; j < n_diag; j++) {
            if (diag_t[j] >= t_discard) { j_start = j; break; }
        }
        int n_analysis = n_diag - j_start;
        printf("  Phonon analysis: t=[%.0f, %.0f], %d samples\n",
               diag_t[j_start], diag_t[n_diag-1], n_analysis);

        if (n_analysis >= 20) {
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

            /* Remove COM drift */
            for (int j = 0; j < n_analysis; j++) {
                double com = 0;
                for (int n = 0; n < N_osc; n++) com += u_n[n][j];
                com /= N_osc;
                for (int n = 0; n < N_osc; n++) u_n[n][j] -= com;
            }

            /* Normal mode decomposition */
            double T_an = diag_t[n_diag-1] - diag_t[j_start];
            double dt_diag_eff = T_an / (n_analysis - 1);

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

            int n_omega = n_analysis / 2;
            if (n_omega > 4000) n_omega = 4000;
            double omega_nyquist = M_PI / dt_diag_eff;
            double omega_max = 1.0;
            if (omega_max > omega_nyquist) omega_max = omega_nyquist;

            double *omega_peak = calloc(N_osc, sizeof(double));
            double *power_peak = calloc(N_osc, sizeof(double));

            char spath[600];
            snprintf(spath, sizeof(spath), "%s/combo235_%s_phonon_2d.tsv", outdir, label);
            FILE *f2d = fopen(spath, "w");
            if (f2d) fprintf(f2d, "k\tomega\tpower\n");

            for (int q = 0; q <= N_osc / 2; q++) {
                double k_q = 2.0 * M_PI * q / (N_osc * d_space);
                double best_pow = 0, best_omega = 0;

                for (int l = 1; l <= n_omega; l++) {
                    double omega = omega_max * l / n_omega;
                    double dft_re = 0, dft_im = 0;
                    for (int j = 0; j < n_analysis; j++) {
                        double t_j = j * dt_diag_eff;
                        double c = cos(omega * t_j);
                        double s = sin(omega * t_j);
                        dft_re += (Q_re[q][j] * c + Q_im[q][j] * s) * dt_diag_eff;
                        dft_im += (-Q_re[q][j] * s + Q_im[q][j] * c) * dt_diag_eff;
                    }
                    double power = (dft_re * dft_re + dft_im * dft_im) / (T_an * T_an);

                    if (f2d) fprintf(f2d, "%.6f\t%.6f\t%.8e\n", k_q, omega, power);

                    if (power > best_pow) {
                        best_pow = power;
                        best_omega = omega;
                    }
                }

                omega_peak[q] = best_omega;
                power_peak[q] = best_pow;
                double cs = (k_q > 1e-10) ? best_omega / k_q : 0;
                printf("  Phonon q=%d: k=%.5f  omega=%.5f  c_s=%.4f\n",
                       q, k_q, best_omega, cs);
            }

            if (f2d) { fclose(f2d); printf("  Wrote %s\n", spath); }

            /* Write modes */
            {
                char mpath[600];
                snprintf(mpath, sizeof(mpath), "%s/combo235_%s_modes.tsv", outdir, label);
                FILE *fp = fopen(mpath, "w");
                if (fp) {
                    fprintf(fp, "q\tk\tomega_peak\tpower_peak\n");
                    for (int q = 0; q <= N_osc / 2; q++) {
                        double k_q = 2.0 * M_PI * q / (N_osc * d_space);
                        fprintf(fp, "%d\t%.6f\t%.6f\t%.8e\n",
                                q, k_q, omega_peak[q], power_peak[q]);
                    }
                    fclose(fp);
                    printf("  Wrote %s\n", mpath);
                }
            }

            /* Store output for comparison */
            for (int q = 0; q <= N_osc / 2; q++)
                out_omega_sym[q] = omega_peak[q];

            /* Max displacement */
            double max_u = 0;
            for (int j = j_start; j < n_diag; j++)
                if (diag_maxdisp[j] > max_u) max_u = diag_maxdisp[j];
            *out_max_disp = max_u;

            /* Energy conservation */
            *out_E_consv = (diag_E[n_diag-1] - diag_E[0]) / (fabs(diag_E[0]) + 1e-30);

            printf("  Max displacement: %.4f (%.2f%% of d)\n",
                   max_u, 100.0 * max_u / d_space);
            printf("  Energy conservation: %.2e\n", *out_E_consv);

            if (N_osc >= 2) {
                double k1 = 2.0 * M_PI / (N_osc * d_space);
                double cs_est = omega_peak[1] / k1;
                printf("  Sound speed (q=1): c_s = %.4f\n", cs_est);
            }

            for (int q = 0; q < N_osc; q++) { free(Q_re[q]); free(Q_im[q]); }
            free(Q_re); free(Q_im);
            for (int n = 0; n < N_osc; n++) free(u_n[n]);
            free(u_n);
            free(omega_peak); free(power_peak);
        }

        free(diag_t); free(diag_E); free(diag_maxdisp);
        for (int n = 0; n < N_osc; n++) free(diag_x[n]);
        free(diag_x);

    } else {
        /* ===== Phase 3: Conservative settle + antisymmetric perturbation ===== */
        printf("  Conservative settle for t=%.0f...\n", t_settle);
        {
            int Nt_settle = (int)(t_settle / dt) + 1;
            int print_ev = Nt_settle / 5;
            if (print_ev < 1) print_ev = 1;

            for (int step = 0; step <= Nt_settle; step++) {
                if (step % print_ev == 0) {
                    double Etot = 0;
                    for (int i = 0; i < Nx; i++)
                        Etot += energy_density_p(phi[0], phi[1], phi[2],
                                                  vel[0], vel[1], vel[2], i, Nx, dx) * dx;
                    printf("    settle t=%7.1f  E=%.4f\n", step * dt, Etot);
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

        /* Apply antisymmetric perturbation */
        printf("  Applying antisymmetric perturbation at oscillon #%d (eps=%.6f)\n",
               pert_osc, eps_pert);
        find_positions(phi[0], phi[1], phi[2],
                       vel[0], vel[1], vel[2], Nx, dx, L, x_osc, N_osc);

        double S_pre[8], A12_pre[8], A13_pre[8];
        measure_asym(phi, Nx, dx, L, x_osc, N_osc, S_pre, A12_pre, A13_pre);

        double E_before = 0;
        for (int i = 0; i < Nx; i++)
            E_before += energy_density_p(phi[0], phi[1], phi[2],
                                          vel[0], vel[1], vel[2], i, Nx, dx) * dx;

        /* phi_1 += +eps, phi_2 -= eps at oscillon #pert_osc */
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
        printf("  E_before=%.6f  E_after=%.6f  Delta=%.6e\n",
               E_before, E_after, E_after - E_before);

        /* Evolve post-perturbation */
        printf("  Evolving post-perturbation for t=%.0f...\n", t_evolve);

        int Nt_ev = (int)(t_evolve / dt) + 1;
        double diag_dt_asym = 5.0;
        int diag_ev = (int)(diag_dt_asym / dt);
        if (diag_ev < 1) diag_ev = 1;
        int max_diag = Nt_ev / diag_ev + 10;
        int print_every = Nt_ev / 20;
        if (print_every < 1) print_every = 1;

        double *diag_t = malloc(max_diag * sizeof(double));
        double *diag_E = malloc(max_diag * sizeof(double));
        double **diag_A12 = malloc(N_osc * sizeof(double*));
        double **diag_A13 = malloc(N_osc * sizeof(double*));
        double **diag_S   = malloc(N_osc * sizeof(double*));
        double **diag_xpos = malloc(N_osc * sizeof(double*));
        for (int n = 0; n < N_osc; n++) {
            diag_A12[n]  = malloc(max_diag * sizeof(double));
            diag_A13[n]  = malloc(max_diag * sizeof(double));
            diag_S[n]    = malloc(max_diag * sizeof(double));
            diag_xpos[n] = malloc(max_diag * sizeof(double));
        }
        int n_diag = 0;

        for (int step = 0; step <= Nt_ev; step++) {
            double t = step * dt;

            if (step % diag_ev == 0 && n_diag < max_diag) {
                if (step > 0)
                    find_positions(phi[0], phi[1], phi[2],
                                   vel[0], vel[1], vel[2],
                                   Nx, dx, L, x_osc, N_osc);

                double Etot = 0;
                for (int i = 0; i < Nx; i++)
                    Etot += energy_density_p(phi[0], phi[1], phi[2],
                                              vel[0], vel[1], vel[2], i, Nx, dx) * dx;

                double SA[8], A12A[8], A13A[8];
                measure_asym(phi, Nx, dx, L, x_osc, N_osc, SA, A12A, A13A);

                diag_t[n_diag] = t;
                diag_E[n_diag] = Etot;
                for (int n = 0; n < N_osc; n++) {
                    diag_S[n][n_diag]    = SA[n];
                    diag_A12[n][n_diag]  = A12A[n];
                    diag_A13[n][n_diag]  = A13A[n];
                    diag_xpos[n][n_diag] = x_osc[n];
                }

                if (step % print_every == 0) {
                    double max_a12 = 0;
                    int max_n = 0;
                    for (int n = 0; n < N_osc; n++) {
                        double a = fabs(A12A[n]);
                        if (a > max_a12) { max_a12 = a; max_n = n; }
                    }
                    printf("  t=%7.1f  E=%.4f  max|A12|=%.6e at osc %d\n",
                           t, Etot, max_a12, max_n);
                }
                n_diag++;
            }

            if (step == Nt_ev) break;

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

        printf("  Evolution complete. %d diagnostic samples.\n", n_diag);

        /* Write antisymmetric time series */
        {
            char path[600];
            snprintf(path, sizeof(path), "%s/combo235_%s_asym_ts.tsv", outdir, label);
            FILE *fp = fopen(path, "w");
            if (fp) {
                fprintf(fp, "time");
                for (int n = 0; n < N_osc; n++) fprintf(fp, "\tA12_%d", n);
                for (int n = 0; n < N_osc; n++) fprintf(fp, "\tS_%d", n);
                fprintf(fp, "\n");
                for (int j = 0; j < n_diag; j++) {
                    fprintf(fp, "%.4f", diag_t[j]);
                    for (int n = 0; n < N_osc; n++)
                        fprintf(fp, "\t%.8e", diag_A12[n][j]);
                    for (int n = 0; n < N_osc; n++)
                        fprintf(fp, "\t%.8e", diag_S[n][j]);
                    fprintf(fp, "\n");
                }
                fclose(fp);
                printf("  Wrote %s\n", path);
            }
        }

        /* ===== Antisymmetric mode spectrum (spatial + temporal DFT) ===== */
        printf("\n  --- Antisymmetric mode spectrum ---\n");

        double T_an = diag_t[n_diag-1] - diag_t[0];
        double dt_diag_eff = T_an / (n_diag - 1);

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

        int n_omega = n_diag / 2;
        if (n_omega > 4000) n_omega = 4000;
        double omega_nyquist = M_PI / dt_diag_eff;
        double omega_max = 2.0;
        if (omega_max > omega_nyquist) omega_max = omega_nyquist;

        double *asym_omega_peak = calloc(N_osc, sizeof(double));
        double *asym_power_peak = calloc(N_osc, sizeof(double));

        char aspec_path[600];
        snprintf(aspec_path, sizeof(aspec_path), "%s/combo235_%s_asym_spectrum_2d.tsv",
                 outdir, label);
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
                double power = (dft_re * dft_re + dft_im * dft_im) / (T_an * T_an);

                if (faspec) fprintf(faspec, "%d\t%.6f\t%.6f\t%.8e\n", q, k_q, omega, power);

                if (power > best_pow) {
                    best_pow = power;
                    best_omega = omega;
                }
            }

            asym_omega_peak[q] = best_omega;
            asym_power_peak[q] = best_pow;
            double cs = (k_q > 1e-10) ? best_omega / k_q : 0;
            printf("  Asym q=%d: k=%.5f  omega=%.5f  power=%.4e  c_s=%.4f\n",
                   q, k_q, best_omega, best_pow, cs);
        }
        if (faspec) { fclose(faspec); printf("  Wrote %s\n", aspec_path); }

        /* Store asym omega for comparison */
        for (int q = 0; q <= N_osc / 2; q++)
            out_omega_asym[q] = asym_omega_peak[q];
        *out_asym_gap = asym_omega_peak[0];  /* q=0 gap */

        /* Also do symmetric displacement spectrum */
        printf("\n  --- Symmetric displacement spectrum ---\n");

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

        for (int j = 0; j < n_diag; j++) {
            double com = 0;
            for (int n = 0; n < N_osc; n++) com += u_n[n][j];
            com /= N_osc;
            for (int n = 0; n < N_osc; n++) u_n[n][j] -= com;
        }

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
                double power = (dft_re * dft_re + dft_im * dft_im) / (T_an * T_an);
                if (power > best_pow) {
                    best_pow = power;
                    best_omega = omega;
                }
            }

            sym_omega_peak[q] = best_omega;
            double cs = (k_q > 1e-10) ? best_omega / k_q : 0;
            printf("  Sym  q=%d: k=%.5f  omega=%.5f  c_s=%.4f\n",
                   q, k_q, best_omega, cs);
        }

        for (int q = 0; q <= N_osc / 2; q++)
            out_omega_sym[q] = sym_omega_peak[q];

        /* ===== Propagation analysis ===== */
        printf("\n  --- Propagation analysis ---\n");

        double threshold = 0.1 * fabs(diag_A12[pert_osc][1]);
        double *arrival_t = malloc(N_osc * sizeof(double));
        for (int n = 0; n < N_osc; n++) arrival_t[n] = -1;

        for (int n = 0; n < N_osc; n++) {
            if (n == pert_osc) { arrival_t[n] = 0; continue; }
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

        double *peak_A12 = calloc(N_osc, sizeof(double));
        for (int n = 0; n < N_osc; n++) {
            for (int j = 0; j < n_diag; j++) {
                double a = fabs(diag_A12[n][j]);
                if (a > peak_A12[n]) peak_A12[n] = a;
            }
        }

        printf("  Propagation results:\n");
        for (int n = 0; n < N_osc; n++) {
            int dist = abs(n - pert_osc);
            if (dist > N_osc / 2) dist = N_osc - dist;
            double d_travel = dist * d_space;
            double speed = (arrival_t[n] > 0 && dist > 0) ? d_travel / arrival_t[n] : 0;
            printf("    osc %d: dist=%d  peak|A12|=%.6e  t_arr=%.1f  c=%.4f\n",
                   n, dist, peak_A12[n], arrival_t[n], speed);
        }

        /* Decay at source */
        double A12_initial = fabs(diag_A12[pert_osc][1]);
        double A12_final = 0;
        {
            int j_start2 = (int)(0.9 * n_diag);
            double sum = 0;
            int count = 0;
            for (int j = j_start2; j < n_diag; j++) {
                sum += fabs(diag_A12[pert_osc][j]);
                count++;
            }
            if (count > 0) A12_final = sum / count;
        }
        double lifetime_ratio = (A12_initial > 0) ? A12_final / A12_initial : 0;
        printf("  Source decay: initial=%.6e  final=%.6e  ratio=%.4f\n",
               A12_initial, A12_final, lifetime_ratio);

        /* Write propagation data */
        {
            char path[600];
            snprintf(path, sizeof(path), "%s/combo235_%s_propagation.tsv", outdir, label);
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

        /* Write modes */
        {
            char mpath[600];
            snprintf(mpath, sizeof(mpath), "%s/combo235_%s_asym_modes.tsv", outdir, label);
            FILE *fp = fopen(mpath, "w");
            if (fp) {
                fprintf(fp, "q\tk\tomega_asym\tpower_asym\tomega_sym\n");
                for (int q = 0; q <= N_osc / 2; q++) {
                    double k_q = 2.0 * M_PI * q / (N_osc * d_space);
                    fprintf(fp, "%d\t%.6f\t%.6f\t%.8e\t%.6f\n",
                            q, k_q, asym_omega_peak[q], asym_power_peak[q],
                            sym_omega_peak[q]);
                }
                fclose(fp);
                printf("  Wrote %s\n", mpath);
            }
        }

        /* Energy conservation */
        *out_E_consv = (diag_E[n_diag-1] - diag_E[0]) / (fabs(diag_E[0]) + 1e-30);

        /* Max displacement */
        {
            double max_u = 0;
            for (int j = 0; j < n_diag; j++) {
                for (int n = 0; n < N_osc; n++) {
                    double disp = diag_xpos[n][j] - (n * d_space + d_space / 2.0);
                    while (disp > L / 2) disp -= L;
                    while (disp < -L / 2) disp += L;
                    if (fabs(disp) > max_u) max_u = fabs(disp);
                }
            }
            *out_max_disp = max_u;
        }

        /* Cleanup */
        free(diag_t); free(diag_E);
        for (int n = 0; n < N_osc; n++) {
            free(diag_A12[n]); free(diag_A13[n]);
            free(diag_S[n]); free(diag_xpos[n]);
            free(u_n[n]);
        }
        free(diag_A12); free(diag_A13); free(diag_S); free(diag_xpos); free(u_n);
        for (int q = 0; q < N_osc; q++) {
            free(AQ_re[q]); free(AQ_im[q]);
            free(SQ_re[q]); free(SQ_im[q]);
        }
        free(AQ_re); free(AQ_im); free(SQ_re); free(SQ_im);
        free(asym_omega_peak); free(asym_power_peak);
        free(sym_omega_peak);
        free(arrival_t); free(peak_A12);
    }

    #undef COMPUTE_ACC_P

    for (int a = 0; a < 3; a++) { free(phi[a]); free(vel[a]); free(acc[a]); }
    free(x_osc); free(x_osc_ref);
}

/* ===================================================================
 *  Main
 * =================================================================== */
int main(int argc, char **argv)
{
    for (int i = 1; i < argc - 1; i += 2) {
        if      (!strcmp(argv[i], "-sigma_c"))  sigma_c  = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-kc"))       kc       = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-eps_reg"))  eps_reg  = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-mu"))       mu_std   = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-kappa"))    kappa_std = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-mass"))     mass     = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-lambda"))   lambda   = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-A"))        A_init   = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-sigma"))    sig      = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-o"))        strncpy(outdir, argv[i+1], sizeof(outdir)-1);
        else if (!strcmp(argv[i], "-N"))        N_osc    = atoi(argv[i+1]);
        else if (!strcmp(argv[i], "-d"))        d_space  = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-t_equil"))  t_equil  = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-t_chain"))  t_chain  = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-t_settle")) t_settle = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-t_evolve")) t_evolve = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-eps_pert")) eps_pert = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-pert_osc")) pert_osc = atoi(argv[i+1]);
        else if (!strcmp(argv[i], "-diag_dt"))  diag_dt  = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-pts_per_d")) pts_per_d = atoi(argv[i+1]);
        else { fprintf(stderr, "Unknown flag: %s\n", argv[i]); exit(1); }
    }

    printf("=== Combo 2+3+5: Confined Lattice with Inertial Deformation ===\n");
    printf("  Confining: sigma=%.3f kc=%.4f eps=%.2e\n", sigma_c, kc, eps_reg);
    printf("  Standard:  mu=%.1f kappa=%.1f\n", mu_std, kappa_std);
    printf("  Pairwise:  lambda=%.4f\n", lambda);
    printf("  Chain: N=%d d=%.1f pts_per_d=%d\n", N_osc, d_space, pts_per_d);
    printf("  Perturbation: eps=%.6f at osc #%d\n", eps_pert, pert_osc);

    /* Storage for comparison results */
    double conf_omega_sym[5], conf_omega_asym[5];
    double conf_max_disp_p2, conf_E_consv_p2;
    double conf_max_disp_p3, conf_E_consv_p3;
    double conf_asym_gap = 0;
    double conf_mass, conf_omega;

    double std_omega_sym[5], std_omega_asym[5];
    double std_max_disp_p2, std_E_consv_p2;
    double std_max_disp_p3, std_E_consv_p3;
    double std_asym_gap = 0;
    double std_mass, std_omega;

    memset(conf_omega_sym, 0, sizeof(conf_omega_sym));
    memset(conf_omega_asym, 0, sizeof(conf_omega_asym));
    memset(std_omega_sym, 0, sizeof(std_omega_sym));
    memset(std_omega_asym, 0, sizeof(std_omega_asym));

    /* ============================================================
     *  Run 1: CONFINING potential
     * ============================================================ */
    printf("\n\n########################################\n");
    printf("### RUN 1: CONFINING POTENTIAL ###\n");
    printf("########################################\n\n");

    pot_mode = 0;

    /* Phase 1: equilibrate */
    run_equilibrate("confine");
    conf_mass = eq_mass_osc;
    conf_omega = eq_omega_osc;

    /* Phase 2: chain + phonon spectrum */
    run_chain("confine_p2", 0,
              conf_omega_sym, conf_omega_asym,
              &conf_max_disp_p2, &conf_E_consv_p2,
              &conf_asym_gap);

    /* Free eq profile, re-equilibrate for Phase 3 */
    for (int a = 0; a < 3; a++) { free(eq_phi[a]); free(eq_vel[a]); }
    run_equilibrate("confine_p3");

    /* Phase 3: antisymmetric perturbation */
    run_chain("confine_p3", 1,
              conf_omega_sym, conf_omega_asym,
              &conf_max_disp_p3, &conf_E_consv_p3,
              &conf_asym_gap);

    for (int a = 0; a < 3; a++) { free(eq_phi[a]); free(eq_vel[a]); }

    /* ============================================================
     *  Run 2: STANDARD potential (comparison)
     * ============================================================ */
    printf("\n\n########################################\n");
    printf("### RUN 2: STANDARD POTENTIAL ###\n");
    printf("########################################\n\n");

    pot_mode = 1;

    /* Phase 1 */
    run_equilibrate("standard");
    std_mass = eq_mass_osc;
    std_omega = eq_omega_osc;

    /* Phase 2 */
    run_chain("standard_p2", 0,
              std_omega_sym, std_omega_asym,
              &std_max_disp_p2, &std_E_consv_p2,
              &std_asym_gap);

    for (int a = 0; a < 3; a++) { free(eq_phi[a]); free(eq_vel[a]); }
    run_equilibrate("standard_p3");

    /* Phase 3 */
    run_chain("standard_p3", 1,
              std_omega_sym, std_omega_asym,
              &std_max_disp_p3, &std_E_consv_p3,
              &std_asym_gap);

    for (int a = 0; a < 3; a++) { free(eq_phi[a]); free(eq_vel[a]); }

    /* ============================================================
     *  COMPARISON SUMMARY
     * ============================================================ */
    printf("\n\n========================================\n");
    printf("=== COMPARISON: CONFINING vs STANDARD ===\n");
    printf("========================================\n\n");

    printf("  Single oscillon:\n");
    printf("    %-20s  %-12s  %-12s\n", "", "Confining", "Standard");
    printf("    %-20s  %-12.4f  %-12.4f\n", "Mass M:", conf_mass, std_mass);
    printf("    %-20s  %-12.4f  %-12.4f\n", "Breathing omega:", conf_omega, std_omega);

    printf("\n  Phase 2 (phonon chain):\n");
    printf("    %-20s  %-12.4f  %-12.4f\n", "Max displacement:", conf_max_disp_p2, std_max_disp_p2);
    printf("    %-20s  %-12.2e  %-12.2e\n", "E conservation:", conf_E_consv_p2, std_E_consv_p2);

    printf("\n  Phonon dispersion (q, omega_confine, omega_standard):\n");
    for (int q = 0; q <= N_osc / 2; q++) {
        double k_q = 2.0 * M_PI * q / (N_osc * d_space);
        printf("    q=%d k=%.5f: conf=%.5f  std=%.5f\n",
               q, k_q, conf_omega_sym[q], std_omega_sym[q]);
    }

    if (N_osc >= 2) {
        double k1 = 2.0 * M_PI / (N_osc * d_space);
        double cs_conf = (k1 > 0) ? conf_omega_sym[1] / k1 : 0;
        double cs_std  = (k1 > 0) ? std_omega_sym[1] / k1 : 0;
        printf("    Sound speed: conf=%.4f  std=%.4f\n", cs_conf, cs_std);
    }

    printf("\n  Phase 3 (antisymmetric perturbation):\n");
    printf("    %-20s  %-12.4f  %-12.4f\n", "Max displacement:", conf_max_disp_p3, std_max_disp_p3);
    printf("    %-20s  %-12.2e  %-12.2e\n", "E conservation:", conf_E_consv_p3, std_E_consv_p3);
    printf("    %-20s  %-12.5f  %-12.5f\n", "Asym gap (q=0):", conf_asym_gap, std_asym_gap);

    printf("\n  Antisymmetric dispersion (q, omega_confine, omega_standard):\n");
    for (int q = 0; q <= N_osc / 2; q++) {
        double k_q = 2.0 * M_PI * q / (N_osc * d_space);
        printf("    q=%d k=%.5f: conf=%.5f  std=%.5f\n",
               q, k_q, conf_omega_asym[q], std_omega_asym[q]);
    }

    printf("\n  KEY QUESTIONS:\n");
    printf("    1. Confined lattice MORE stable? max_disp: %.4f vs %.4f => %s\n",
           conf_max_disp_p2, std_max_disp_p2,
           conf_max_disp_p2 < std_max_disp_p2 ? "YES (confined tighter)" : "NO");
    printf("    2. Antisymmetric gap lower? conf=%.5f vs std=%.5f => %s\n",
           conf_asym_gap, std_asym_gap,
           conf_asym_gap < std_asym_gap ? "YES (lower gap)" : "NO (higher or same)");
    printf("    3. Breathing omega: conf=%.4f vs std=%.4f (deeper binding => lower omega)\n",
           conf_omega, std_omega);

    printf("\n=== Combo 2+3+5 Complete ===\n");
    return 0;
}
