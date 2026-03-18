/*  t12_char_m0.c — Control: standard EOM, no self-limiting mechanism
 *  N=128, L=30, T=500, mass2=2.25, bimodal braid + background
 *  Runs braid+bg and bg-only simultaneously for differential measurement.
 *
 *  Build: gcc -O3 -fopenmp -o t12_char_m0 src/t12_char_m0.c -lm
 */

#include "../../src/braid_core.h"

#define NBINS   100
#define A_BG    0.1
#define T_TOTAL 500.0
#define SNAP_DT 10.0
#define N_SNAPS 50
#define N_DUMPS 6

static double rho0_bg, mass2_phys;

/* ================================================================
   Shared utilities
   ================================================================ */

static void add_background(Grid *g) {
    int N = g->N, NN = N*N;
    double dx = g->dx, L = g->L;
    double k_bg = PI / L, omega_bg = sqrt(k_bg*k_bg + mass2_phys);
    for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
    for (int k = 0; k < N; k++) {
        int idx = i*NN + j*N + k;
        double z = -L + k*dx;
        for (int a = 0; a < NFIELDS; a++) {
            double phase = k_bg*z + 2.0*PI*a/3.0;
            g->phi[a][idx] += A_BG * cos(phase);
            g->vel[a][idx] += omega_bg * A_BG * sin(phase);
        }
    }
}

static void init_bg_only(Grid *g) {
    int N = g->N, NN = N*N;
    double dx = g->dx, L = g->L;
    double k_bg = PI / L, omega_bg = sqrt(k_bg*k_bg + mass2_phys);
    grid_zero(g);
    for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
    for (int k = 0; k < N; k++) {
        int idx = i*NN + j*N + k;
        double z = -L + k*dx;
        for (int a = 0; a < NFIELDS; a++) {
            double phase = k_bg*z + 2.0*PI*a/3.0;
            g->phi[a][idx] = A_BG * cos(phase);
            g->vel[a][idx] = omega_bg * A_BG * sin(phase);
        }
    }
}

static void apply_edge_damping(Grid *g) {
    int N = g->N, NN = N*N;
    double dx = g->dx, L = g->L;
    double r_start = 0.85*L, r_end = 0.98*L;
    double inv_dr = 1.0 / (r_end - r_start + 1e-30);
    for (int i = 0; i < N; i++) {
        double x = -L + i*dx;
        for (int j = 0; j < N; j++) {
            double y = -L + j*dx;
            double rp = sqrt(x*x + y*y);
            if (rp <= r_start) continue;
            double f = (rp - r_start) * inv_dr;
            if (f > 1.0) f = 1.0;
            double damp = 1.0 - 0.95*f*f;
            for (int kk = 0; kk < N; kk++) {
                int idx = i*NN + j*N + kk;
                for (int a = 0; a < NFIELDS; a++)
                    g->vel[a][idx] *= damp;
            }
        }
    }
}

/* Energy density at a single point (kinetic + gradient + mass) */
static double point_energy(Grid *g, int i, int j, int k) {
    int N = g->N, NN = N*N;
    double dx = g->dx;
    int idx = i*NN + j*N + k;
    int kp = (k+1)%N, km = (k-1+N)%N;
    double e = 0;
    for (int a = 0; a < NFIELDS; a++) {
        e += 0.5 * g->vel[a][idx] * g->vel[a][idx];
        double gx = (g->phi[a][idx+NN] - g->phi[a][idx-NN]) / (2*dx);
        double gy = (g->phi[a][idx+N]  - g->phi[a][idx-N])  / (2*dx);
        double gz = (g->phi[a][i*NN+j*N+kp] - g->phi[a][i*NN+j*N+km]) / (2*dx);
        e += 0.5*(gx*gx + gy*gy + gz*gz);
        e += 0.5*mass2_phys * g->phi[a][idx] * g->phi[a][idx];
    }
    return e;
}

/* Radial profile of energy density */
static void compute_radial_rho(Grid *g, double *r_bins, double *rho_bins,
                                int *counts, double dr) {
    int N = g->N;
    double dx = g->dx, L = g->L;
    for (int b = 0; b < NBINS; b++) { r_bins[b] = (b+0.5)*dr; rho_bins[b] = 0; counts[b] = 0; }
    for (int i = 1; i < N-1; i++) {
        double x = -L + i*dx;
        for (int j = 1; j < N-1; j++) {
            double y = -L + j*dx;
            double rp = sqrt(x*x + y*y);
            int b = (int)(rp / dr);
            if (b >= NBINS) continue;
            for (int k = 0; k < N; k++) {
                rho_bins[b] += point_energy(g, i, j, k);
                counts[b]++;
            }
        }
    }
    for (int b = 0; b < NBINS; b++)
        if (counts[b] > 0) rho_bins[b] /= counts[b];
}

/* Energy in radial shells */
static void compute_shell_energy(Grid *g, double *E_shells) {
    /* shells: 0-5, 5-10, 10-15, 15-20, 20-25 */
    double bounds[6] = {0, 5, 10, 15, 20, 25};
    int N = g->N;
    double dx = g->dx, L = g->L, dV = dx*dx*dx;
    for (int s = 0; s < 5; s++) E_shells[s] = 0;
    for (int i = 1; i < N-1; i++) {
        double x = -L + i*dx;
        for (int j = 1; j < N-1; j++) {
            double y = -L + j*dx;
            double rp = sqrt(x*x + y*y);
            for (int s = 0; s < 5; s++) {
                if (rp >= bounds[s] && rp < bounds[s+1]) {
                    for (int k = 0; k < N; k++)
                        E_shells[s] += point_energy(g, i, j, k) * dV;
                    break;
                }
            }
        }
    }
}

static double compute_fc(Grid *g) {
    int N = g->N, NN = N*N;
    double dx = g->dx, L = g->L;
    double Rc2 = 64.0, total = 0, core = 0;
    for (int i = 1; i < N-1; i++) {
        double x = -L + i*dx;
        for (int j = 1; j < N-1; j++) {
            double y = -L + j*dx;
            double rp2 = x*x + y*y;
            for (int k = 0; k < N; k++) {
                int idx = i*NN + j*N + k;
                double rho = 0;
                for (int a = 0; a < NFIELDS; a++)
                    rho += g->phi[a][idx] * g->phi[a][idx];
                total += rho;
                if (rp2 < Rc2) core += rho;
            }
        }
    }
    return core / (total + 1e-30);
}

static double compute_energy(Grid *g, const double *phys) {
    int N = g->N, NN = N*N;
    double dx = g->dx, dV = dx*dx*dx;
    double mu = phys[12], kappa = phys[13], lpw = phys[15];
    double E = 0;
    for (int i = 1; i < N-1; i++)
    for (int j = 1; j < N-1; j++)
    for (int k = 0; k < N; k++) {
        int idx = i*NN + j*N + k;
        int kp = (k+1)%N, km = (k-1+N)%N;
        double ek = 0, eg = 0, em = 0;
        for (int a = 0; a < NFIELDS; a++) {
            ek += 0.5*g->vel[a][idx]*g->vel[a][idx];
            double gx = (g->phi[a][idx+NN] - g->phi[a][idx-NN]) / (2*dx);
            double gy = (g->phi[a][idx+N]  - g->phi[a][idx-N])  / (2*dx);
            double gz = (g->phi[a][i*NN+j*N+kp] - g->phi[a][i*NN+j*N+km]) / (2*dx);
            eg += 0.5*(gx*gx + gy*gy + gz*gz);
            em += 0.5*mass2_phys*g->phi[a][idx]*g->phi[a][idx];
        }
        double p0 = g->phi[0][idx], p1 = g->phi[1][idx], p2 = g->phi[2][idx];
        double P = p0*p1*p2;
        double ep = (mu/2.0)*P*P / (1.0 + kappa*P*P);
        double epw = lpw*(p0*p1 + p1*p2 + p2*p0);
        E += (ek + eg + em + ep + epw) * dV;
    }
    return E;
}

static double compute_peak_P(Grid *g) {
    int N = g->N, NN = N*N;
    double peak = 0;
    for (int i = 1; i < N-1; i++)
    for (int j = 1; j < N-1; j++)
    for (int k = 0; k < N; k++) {
        int idx = i*NN + j*N + k;
        double P = fabs(g->phi[0][idx] * g->phi[1][idx] * g->phi[2][idx]);
        if (P > peak) peak = P;
    }
    return peak;
}

/* Write xy-slice of energy density at z=z_mid */
static void dump_rho_xy(Grid *g, double t, const char *prefix) {
    int N = g->N, NN = N*N;
    double dx = g->dx, L = g->L;
    int k_mid = N/2;
    char fname[256];
    snprintf(fname, sizeof(fname), "data/%s_rho_xy_t%03d.tsv", prefix, (int)t);
    FILE *f = fopen(fname, "w");
    fprintf(f, "x\ty\trho\n");
    for (int i = 1; i < N-1; i++) {
        double x = -L + i*dx;
        for (int j = 1; j < N-1; j++) {
            double y = -L + j*dx;
            double e = point_energy(g, i, j, k_mid);
            fprintf(f, "%.4f\t%.4f\t%.6e\n", x, y, e);
        }
    }
    fclose(f);
}

/* Write z-axis field profile at (x_mid, y_mid) */
static void dump_zaxis(Grid *g, double t, const char *prefix) {
    int N = g->N, NN = N*N;
    double dx = g->dx, L = g->L;
    int i_mid = N/2, j_mid = N/2;
    char fname[256];
    snprintf(fname, sizeof(fname), "data/%s_zaxis_t%03d.tsv", prefix, (int)t);
    FILE *f = fopen(fname, "w");
    fprintf(f, "z\tphi0\tphi1\tphi2\n");
    for (int k = 0; k < N; k++) {
        double z = -L + k*dx;
        int idx = i_mid*NN + j_mid*N + k;
        fprintf(f, "%.4f\t%.6e\t%.6e\t%.6e\n", z, g->phi[0][idx], g->phi[1][idx], g->phi[2][idx]);
    }
    fclose(f);
}

/* Fit delta_rho(r) to power law and exponential */
static void fit_profiles(double *r, double *drho, int nb, const char *prefix) {
    /* Power law: log|drho| = log(A) - alpha*log(r) -> linear regression */
    /* Exponential: log|drho| = log(A) - r/lambda -> linear regression */
    double sx = 0, sy_p = 0, sxx = 0, sxy_p = 0;
    double sx_e = 0, sy_e = 0, sxx_e = 0, sxy_e = 0;
    int n_p = 0, n_e = 0;

    for (int b = 0; b < nb; b++) {
        if (r[b] < 2.0 || r[b] > 22.0) continue;
        if (fabs(drho[b]) < 1e-15) continue;
        double ldr = log(fabs(drho[b]));

        /* Power law: use log(r) as x */
        double lr = log(r[b]);
        sx += lr; sy_p += ldr; sxx += lr*lr; sxy_p += lr*ldr; n_p++;

        /* Exponential: use r as x */
        sx_e += r[b]; sy_e += ldr; sxx_e += r[b]*r[b]; sxy_e += r[b]*ldr; n_e++;
    }

    double alpha = 0, A_pw = 0, res_pw = 0;
    double lambda = 0, A_ex = 0, res_ex = 0;

    if (n_p > 2) {
        double det = n_p*sxx - sx*sx;
        double slope = (n_p*sxy_p - sx*sy_p) / (det + 1e-30);
        double intercept = (sy_p - slope*sx) / n_p;
        alpha = -slope;
        A_pw = exp(intercept);

        /* Compute residuals */
        for (int b = 0; b < nb; b++) {
            if (r[b] < 2.0 || r[b] > 22.0 || fabs(drho[b]) < 1e-15) continue;
            double pred = A_pw / pow(r[b], alpha);
            double diff = fabs(drho[b]) - pred;
            res_pw += diff*diff;
        }
    }
    if (n_e > 2) {
        double det = n_e*sxx_e - sx_e*sx_e;
        double slope = (n_e*sxy_e - sx_e*sy_e) / (det + 1e-30);
        double intercept = (sy_e - slope*sx_e) / n_e;
        lambda = -1.0 / (slope + 1e-30);
        A_ex = exp(intercept);

        for (int b = 0; b < nb; b++) {
            if (r[b] < 2.0 || r[b] > 22.0 || fabs(drho[b]) < 1e-15) continue;
            double pred = A_ex * exp(-r[b] / lambda);
            double diff = fabs(drho[b]) - pred;
            res_ex += diff*diff;
        }
    }

    char fname[256];
    snprintf(fname, sizeof(fname), "data/%s_fit.tsv", prefix);
    FILE *f = fopen(fname, "w");
    fprintf(f, "fit_type\tA\texponent\tresidual\n");
    fprintf(f, "power_law\t%.6e\t%.4f\t%.6e\n", A_pw, alpha, res_pw);
    fprintf(f, "exponential\t%.6e\t%.4f\t%.6e\n", A_ex, lambda, res_ex);
    fprintf(f, "# best_fit\t%s\n", (res_pw < res_ex) ? "power_law" : "exponential");
    fclose(f);
    printf("  Fit: power_law alpha=%.3f (res=%.2e), exponential lambda=%.3f (res=%.2e) -> %s\n",
           alpha, res_pw, lambda, res_ex, (res_pw < res_ex) ? "POWER_LAW" : "EXPONENTIAL");
}

/* ================================================================
   Main
   ================================================================ */

int main(void) {
    bimodal_init_params();
    mass2_phys = 2.25;
    int N = 128;
    double L = 30.0;
    double k_bg = PI / L, omega_bg = sqrt(k_bg*k_bg + mass2_phys);
    rho0_bg = 3.0 * A_BG*A_BG * omega_bg*omega_bg;

    printf("=== T12 Char M0: Control (no mechanism) ===\n");
    printf("N=%d, L=%.1f, T=%.1f, mass2=%.4f, rho0_bg=%.6f\n", N, L, T_TOTAL, mass2_phys, rho0_bg);

    /* Braid + background grid */
    Grid *gb = grid_alloc(N, L);
    init_braid(gb, BIMODAL, -1);
    add_background(gb);
    compute_forces(gb, BIMODAL, mass2_phys);

    /* Control: background only */
    Grid *gc = grid_alloc(N, L);
    init_bg_only(gc);
    compute_forces(gc, BIMODAL, mass2_phys);

    double dr = L / NBINS;  /* 0.3 */
    double r_bins[NBINS], rho_braid[NBINS], rho_ctrl[NBINS];
    int counts_b[NBINS], counts_c[NBINS];

    /* Previous profiles for ddrho/dt */
    double prev_drho[NBINS];
    memset(prev_drho, 0, sizeof(prev_drho));
    double prev_t = -1;

    /* Dump timepoints */
    double dump_times[N_DUMPS] = {0, 100, 200, 300, 400, 500};
    int dump_idx = 0;

    /* Late-time profiles for fitting (accumulate T>300) */
    double late_drho[NBINS];
    int late_count = 0;
    memset(late_drho, 0, sizeof(late_drho));

    FILE *fts = fopen("data/charM0_timeseries.tsv", "w");
    fprintf(fts, "t\tfc\tpeak_P\twinding\tenergy\tE_s0\tE_s1\tE_s2\tE_s3\tE_s4\n");

    FILE *fpr = fopen("data/charM0_profiles.tsv", "w");
    fprintf(fpr, "t\tr\trho_braid\trho_ctrl\tdrho\n");

    int n_total = (int)(T_TOTAL / gb->dt);
    int snap_every = (int)(SNAP_DT / gb->dt);
    int snap_count = 0;

    printf("dt=%.5f, steps=%d, snap_every=%d\n", gb->dt, n_total, snap_every);
    double wt0 = omp_get_wtime();

    for (int step = 0; step <= n_total; step++) {
        if (step > 0) {
            /* Advance braid grid */
            verlet_full_step(gb, BIMODAL, mass2_phys);
            apply_edge_damping(gb);

            /* Advance control grid */
            verlet_full_step(gc, BIMODAL, mass2_phys);
            apply_edge_damping(gc);
        }

        if (step % snap_every == 0 && snap_count < N_SNAPS + 1) {
            double t = step * gb->dt;

            /* Diagnostics */
            double fc = compute_fc(gb);
            double peak_P = compute_peak_P(gb);
            double winding = compute_winding(gb);
            double energy = compute_energy(gb, BIMODAL);
            double E_shells[5];
            compute_shell_energy(gb, E_shells);

            /* Radial profiles */
            compute_radial_rho(gb, r_bins, rho_braid, counts_b, dr);
            compute_radial_rho(gc, r_bins, rho_ctrl, counts_c, dr);

            /* Write time series */
            fprintf(fts, "%.2f\t%.6f\t%.6e\t%.4f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n",
                    t, fc, peak_P, winding, energy,
                    E_shells[0], E_shells[1], E_shells[2], E_shells[3], E_shells[4]);

            /* Write profiles */
            for (int b = 0; b < NBINS; b++) {
                double drho = rho_braid[b] - rho_ctrl[b];
                fprintf(fpr, "%.2f\t%.4f\t%.6e\t%.6e\t%.6e\n",
                        t, r_bins[b], rho_braid[b], rho_ctrl[b], drho);
            }

            /* Depletion rate */
            double max_drate = 0;
            if (prev_t >= 0 && t > prev_t) {
                for (int b = 0; b < NBINS; b++) {
                    double drho = rho_braid[b] - rho_ctrl[b];
                    double rate = (drho - prev_drho[b]) / (t - prev_t);
                    if (fabs(rate) > max_drate) max_drate = fabs(rate);
                }
            }

            /* Accumulate late-time profiles */
            if (t > 300.0) {
                for (int b = 0; b < NBINS; b++)
                    late_drho[b] += rho_braid[b] - rho_ctrl[b];
                late_count++;
            }

            /* Save for next rate computation */
            for (int b = 0; b < NBINS; b++)
                prev_drho[b] = rho_braid[b] - rho_ctrl[b];
            prev_t = t;

            printf("t=%6.1f  fc=%.4f  |P|=%.4e  w=%.3f  E=%.1f  max|drate|=%.2e  [%.0fs]\n",
                   t, fc, peak_P, winding, energy, max_drate, omp_get_wtime()-wt0);

            /* Full field dumps */
            if (dump_idx < N_DUMPS && fabs(t - dump_times[dump_idx]) < SNAP_DT*0.5) {
                dump_rho_xy(gb, t, "charM0");
                dump_zaxis(gb, t, "charM0");
                dump_idx++;
            }

            snap_count++;
            if (check_blowup(gb)) { printf("BLOWUP at t=%.1f\n", t); break; }
        }
    }

    fclose(fts);
    fclose(fpr);

    /* Fit late-time average profile */
    if (late_count > 0) {
        for (int b = 0; b < NBINS; b++) late_drho[b] /= late_count;
        printf("\nLate-time (T>300) average depletion profile fit:\n");
        fit_profiles(r_bins, late_drho, NBINS, "charM0");
    }

    grid_free(gb);
    grid_free(gc);
    printf("\n=== M0 Control Complete (%.0f s) ===\n", omp_get_wtime()-wt0);
    return 0;
}
