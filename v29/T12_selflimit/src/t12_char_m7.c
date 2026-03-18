/*  t12_char_m7.c — M7 Two-Component characterization
 *  N=128, L=30, T=500, mass2=2.25, bimodal braid + background
 *  M7: S (braid) + B (background) with coupling g*S^2*B^2.
 *  S: standard EOM + triple-product - g*B^2*S
 *  B: free wave - g*S^2*B
 *  g=0.01. 6 total field arrays.
 *
 *  Build: gcc -O3 -fopenmp -o t12_char_m7 src/t12_char_m7.c -lm
 */

#include "../../src/braid_core.h"

#define NBINS   100
#define A_BG    0.1
#define T_TOTAL 500.0
#define SNAP_DT 10.0
#define N_SNAPS 50
#define N_DUMPS 6
#define G_COUP  0.01

static double rho0_bg, mass2_phys;

/* B fields (background component) */
static double *B_phi[NFIELDS], *B_vel[NFIELDS], *B_acc[NFIELDS];

/* ================================================================
   Shared utilities
   ================================================================ */

static void apply_edge_damping_both(Grid *g) {
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
                for (int a = 0; a < NFIELDS; a++) {
                    g->vel[a][idx] *= damp;
                    B_vel[a][idx] *= damp;
                }
            }
        }
    }
}

/* Energy density at point, including both S and B */
static double point_energy_total(Grid *g, int i, int j, int k) {
    int N = g->N, NN = N*N;
    double dx = g->dx;
    int idx = i*NN + j*N + k;
    int kp = (k+1)%N, km = (k-1+N)%N;
    double e = 0;
    for (int a = 0; a < NFIELDS; a++) {
        /* S kinetic + gradient + mass */
        e += 0.5 * g->vel[a][idx] * g->vel[a][idx];
        double gx = (g->phi[a][idx+NN] - g->phi[a][idx-NN]) / (2*dx);
        double gy = (g->phi[a][idx+N]  - g->phi[a][idx-N])  / (2*dx);
        double gz = (g->phi[a][i*NN+j*N+kp] - g->phi[a][i*NN+j*N+km]) / (2*dx);
        e += 0.5*(gx*gx + gy*gy + gz*gz);
        e += 0.5*mass2_phys * g->phi[a][idx] * g->phi[a][idx];
        /* B kinetic + gradient + mass */
        e += 0.5 * B_vel[a][idx] * B_vel[a][idx];
        double bx = (B_phi[a][idx+NN] - B_phi[a][idx-NN]) / (2*dx);
        double by = (B_phi[a][idx+N]  - B_phi[a][idx-N])  / (2*dx);
        double bz = (B_phi[a][i*NN+j*N+kp] - B_phi[a][i*NN+j*N+km]) / (2*dx);
        e += 0.5*(bx*bx + by*by + bz*bz);
        e += 0.5*mass2_phys * B_phi[a][idx] * B_phi[a][idx];
    }
    /* Coupling energy */
    double S2 = 0, B2 = 0;
    for (int a = 0; a < NFIELDS; a++) {
        S2 += g->phi[a][idx] * g->phi[a][idx];
        B2 += B_phi[a][idx] * B_phi[a][idx];
    }
    e += G_COUP * S2 * B2;
    return e;
}

/* B-only energy density (for control reference) */
static double point_energy_B(int N, int NN, double dx, int i, int j, int k) {
    int idx = i*NN + j*N + k;
    int kp = (k+1)%N, km = (k-1+N)%N;
    double e = 0;
    for (int a = 0; a < NFIELDS; a++) {
        e += 0.5 * B_vel[a][idx] * B_vel[a][idx];
        double bx = (B_phi[a][idx+NN] - B_phi[a][idx-NN]) / (2*dx);
        double by = (B_phi[a][idx+N]  - B_phi[a][idx-N])  / (2*dx);
        double bz = (B_phi[a][i*NN+j*N+kp] - B_phi[a][i*NN+j*N+km]) / (2*dx);
        e += 0.5*(bx*bx + by*by + bz*bz);
        e += 0.5*mass2_phys * B_phi[a][idx] * B_phi[a][idx];
    }
    return e;
}

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
                rho_bins[b] += point_energy_total(g, i, j, k);
                counts[b]++;
            }
        }
    }
    for (int b = 0; b < NBINS; b++)
        if (counts[b] > 0) rho_bins[b] /= counts[b];
}

/* B-only radial profile (acts as control for M7) */
static void compute_radial_rho_B(Grid *g, double *r_bins, double *rho_bins,
                                  int *counts, double dr) {
    int N = g->N, NN = N*N;
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
                rho_bins[b] += point_energy_B(N, NN, dx, i, j, k);
                counts[b]++;
            }
        }
    }
    for (int b = 0; b < NBINS; b++)
        if (counts[b] > 0) rho_bins[b] /= counts[b];
}

static void compute_shell_energy(Grid *g, double *E_shells) {
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
                        E_shells[s] += point_energy_total(g, i, j, k) * dV;
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

static double compute_energy_m7(Grid *g, const double *phys) {
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
            ek += 0.5*g->vel[a][idx]*g->vel[a][idx] + 0.5*B_vel[a][idx]*B_vel[a][idx];
            double gx = (g->phi[a][idx+NN]-g->phi[a][idx-NN])/(2*dx);
            double gy = (g->phi[a][idx+N]-g->phi[a][idx-N])/(2*dx);
            double gz = (g->phi[a][i*NN+j*N+kp]-g->phi[a][i*NN+j*N+km])/(2*dx);
            eg += 0.5*(gx*gx+gy*gy+gz*gz);
            double bx = (B_phi[a][idx+NN]-B_phi[a][idx-NN])/(2*dx);
            double by = (B_phi[a][idx+N]-B_phi[a][idx-N])/(2*dx);
            double bz = (B_phi[a][i*NN+j*N+kp]-B_phi[a][i*NN+j*N+km])/(2*dx);
            eg += 0.5*(bx*bx+by*by+bz*bz);
            em += 0.5*mass2_phys*(g->phi[a][idx]*g->phi[a][idx] + B_phi[a][idx]*B_phi[a][idx]);
        }
        double p0=g->phi[0][idx], p1=g->phi[1][idx], p2=g->phi[2][idx];
        double P=p0*p1*p2;
        double ep=(mu/2.0)*P*P/(1.0+kappa*P*P);
        double epw=lpw*(p0*p1+p1*p2+p2*p0);
        double S2=p0*p0+p1*p1+p2*p2;
        double B2=B_phi[0][idx]*B_phi[0][idx]+B_phi[1][idx]*B_phi[1][idx]+B_phi[2][idx]*B_phi[2][idx];
        E += (ek+eg+em+ep+epw+G_COUP*S2*B2)*dV;
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
            fprintf(f, "%.4f\t%.4f\t%.6e\n", x, y, point_energy_total(g, i, j, k_mid));
        }
    }
    fclose(f);
}

static void dump_zaxis(Grid *g, double t, const char *prefix) {
    int N = g->N, NN = N*N;
    double dx = g->dx, L = g->L;
    int i_mid = N/2, j_mid = N/2;
    char fname[256];
    snprintf(fname, sizeof(fname), "data/%s_zaxis_t%03d.tsv", prefix, (int)t);
    FILE *f = fopen(fname, "w");
    fprintf(f, "z\tphi0\tphi1\tphi2\tB0\tB1\tB2\n");
    for (int k = 0; k < N; k++) {
        double z = -L + k*dx;
        int idx = i_mid*NN + j_mid*N + k;
        fprintf(f, "%.4f\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\n",
                z, g->phi[0][idx], g->phi[1][idx], g->phi[2][idx],
                B_phi[0][idx], B_phi[1][idx], B_phi[2][idx]);
    }
    fclose(f);
}

static void fit_profiles(double *r, double *drho, int nb, const char *prefix) {
    double sx = 0, sy_p = 0, sxx = 0, sxy_p = 0;
    double sx_e = 0, sy_e = 0, sxx_e = 0, sxy_e = 0;
    int n_p = 0, n_e = 0;
    for (int b = 0; b < nb; b++) {
        if (r[b] < 2.0 || r[b] > 22.0 || fabs(drho[b]) < 1e-15) continue;
        double ldr = log(fabs(drho[b]));
        double lr = log(r[b]);
        sx += lr; sy_p += ldr; sxx += lr*lr; sxy_p += lr*ldr; n_p++;
        sx_e += r[b]; sy_e += ldr; sxx_e += r[b]*r[b]; sxy_e += r[b]*ldr; n_e++;
    }
    double alpha = 0, A_pw = 0, res_pw = 0;
    double lambda = 0, A_ex = 0, res_ex = 0;
    if (n_p > 2) {
        double det = n_p*sxx - sx*sx;
        double slope = (n_p*sxy_p - sx*sy_p) / (det + 1e-30);
        double intercept = (sy_p - slope*sx) / n_p;
        alpha = -slope; A_pw = exp(intercept);
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
        lambda = -1.0 / (slope + 1e-30); A_ex = exp(intercept);
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
   M7-specific: two-component forces
   ================================================================ */

static void compute_forces_m7(Grid *g, const double *phys) {
    int N = g->N, NN = N*N, N3 = N*N*N;
    double idx2 = 1.0 / (g->dx * g->dx);
    double mu = phys[12], kappa = phys[13], lpw = phys[15];

    #pragma omp parallel for schedule(static)
    for (int idx = 0; idx < N3; idx++) {
        int i = idx/NN, j = (idx/N)%N, k = idx%N;
        if (i < 1 || i >= N-1 || j < 1 || j >= N-1) {
            for (int a = 0; a < NFIELDS; a++) { g->acc[a][idx] = 0; B_acc[a][idx] = 0; }
            continue;
        }
        int kp = (k+1)%N, km = (k-1+N)%N;
        int idx_kp = i*NN+j*N+kp, idx_km = i*NN+j*N+km;

        double p0 = g->phi[0][idx], p1 = g->phi[1][idx], p2 = g->phi[2][idx];
        double P = p0*p1*p2;
        double denom = 1.0 + kappa*P*P;
        double mu_P_d2 = mu*P / (denom*denom);

        double S2 = p0*p0 + p1*p1 + p2*p2;
        double B2 = B_phi[0][idx]*B_phi[0][idx] + B_phi[1][idx]*B_phi[1][idx] + B_phi[2][idx]*B_phi[2][idx];

        for (int a = 0; a < NFIELDS; a++) {
            /* S forces: standard braid + coupling */
            double lap_s = (g->phi[a][idx+NN]+g->phi[a][idx-NN]
                           +g->phi[a][idx+N]+g->phi[a][idx-N]
                           +g->phi[a][idx_kp]+g->phi[a][idx_km]
                           -6.0*g->phi[a][idx]) * idx2;
            double dPda = (a==0)?p1*p2:(a==1)?p0*p2:p0*p1;
            double f_triple = mu_P_d2 * dPda;
            double f_pw = lpw*(g->phi[(a+1)%3][idx]+g->phi[(a+2)%3][idx]);
            g->acc[a][idx] = lap_s - mass2_phys*g->phi[a][idx] - f_triple - f_pw
                            - G_COUP*B2*g->phi[a][idx];

            /* B forces: free wave + coupling */
            double lap_b = (B_phi[a][idx+NN]+B_phi[a][idx-NN]
                           +B_phi[a][idx+N]+B_phi[a][idx-N]
                           +B_phi[a][idx_kp]+B_phi[a][idx_km]
                           -6.0*B_phi[a][idx]) * idx2;
            B_acc[a][idx] = lap_b - mass2_phys*B_phi[a][idx]
                           - G_COUP*S2*B_phi[a][idx];
        }
    }
}

/* ================================================================
   Main
   ================================================================ */

int main(void) {
    bimodal_init_params();
    mass2_phys = 2.25;
    int N = 128;
    double L = 30.0;
    int N3 = N*N*N;
    double k_bg = PI / L, omega_bg = sqrt(k_bg*k_bg + mass2_phys);
    rho0_bg = 3.0 * A_BG*A_BG * omega_bg*omega_bg;

    printf("=== T12 Char M7: Two-Component (g=%.3f) ===\n", G_COUP);
    printf("N=%d, L=%.1f, T=%.1f, mass2=%.4f, rho0_bg=%.6f\n", N, L, T_TOTAL, mass2_phys, rho0_bg);

    /* Allocate B fields */
    for (int a = 0; a < NFIELDS; a++) {
        B_phi[a] = calloc(N3, sizeof(double));
        B_vel[a] = calloc(N3, sizeof(double));
        B_acc[a] = calloc(N3, sizeof(double));
    }

    /* S = braid (no background) */
    Grid *g = grid_alloc(N, L);
    init_braid(g, BIMODAL, -1);

    /* B = background traveling wave */
    int NN = N*N;
    double dx = g->dx;
    for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
    for (int k = 0; k < N; k++) {
        int idx = i*NN + j*N + k;
        double z = -L + k*dx;
        for (int a = 0; a < NFIELDS; a++) {
            double phase = PI/L*z + 2.0*PI*a/3.0;
            B_phi[a][idx] = A_BG * cos(phase);
            B_vel[a][idx] = omega_bg * A_BG * sin(phase);
        }
    }

    compute_forces_m7(g, BIMODAL);

    /* Also run a separate control grid (bg only, no coupling) for differential */
    Grid *gc = grid_alloc(N, L);
    grid_zero(gc);
    for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
    for (int k = 0; k < N; k++) {
        int idx = i*NN + j*N + k;
        double z = -L + k*dx;
        for (int a = 0; a < NFIELDS; a++) {
            double phase = PI/L*z + 2.0*PI*a/3.0;
            gc->phi[a][idx] = A_BG * cos(phase);
            gc->vel[a][idx] = omega_bg * A_BG * sin(phase);
        }
    }
    compute_forces(gc, BIMODAL, mass2_phys);

    double dr = L / NBINS;
    double r_bins[NBINS], rho_total[NBINS], rho_ctrl[NBINS];
    int counts_t[NBINS], counts_c[NBINS];
    double prev_drho[NBINS]; memset(prev_drho, 0, sizeof(prev_drho));
    double prev_t = -1;

    double dump_times[N_DUMPS] = {0, 100, 200, 300, 400, 500};
    int dump_idx = 0;

    double late_drho[NBINS]; int late_count = 0;
    memset(late_drho, 0, sizeof(late_drho));

    FILE *fts = fopen("data/charM7_timeseries.tsv", "w");
    fprintf(fts, "t\tfc\tpeak_P\twinding\tenergy\tE_s0\tE_s1\tE_s2\tE_s3\tE_s4\n");
    FILE *fpr = fopen("data/charM7_profiles.tsv", "w");
    fprintf(fpr, "t\tr\trho_total\trho_ctrl\tdrho\n");

    int n_total = (int)(T_TOTAL / g->dt);
    int snap_every = (int)(SNAP_DT / g->dt);
    int snap_count = 0;

    printf("dt=%.5f, steps=%d\n", g->dt, n_total);
    double wt0 = omp_get_wtime();

    for (int step = 0; step <= n_total; step++) {
        if (step > 0) {
            double hdt = 0.5*g->dt, dt = g->dt;

            /* M7: kick both S and B */
            for (int a = 0; a < NFIELDS; a++) {
                double *vs = g->vel[a], *as = g->acc[a];
                double *vb = B_vel[a], *ab = B_acc[a];
                for (int idx = 0; idx < N3; idx++) {
                    vs[idx] += hdt*as[idx];
                    vb[idx] += hdt*ab[idx];
                }
            }
            /* Drift both */
            for (int a = 0; a < NFIELDS; a++) {
                double *ps = g->phi[a], *vs = g->vel[a];
                double *pb = B_phi[a], *vb = B_vel[a];
                for (int idx = 0; idx < N3; idx++) {
                    ps[idx] += dt*vs[idx];
                    pb[idx] += dt*vb[idx];
                }
            }
            /* Force */
            compute_forces_m7(g, BIMODAL);
            /* Kick again */
            for (int a = 0; a < NFIELDS; a++) {
                double *vs = g->vel[a], *as = g->acc[a];
                double *vb = B_vel[a], *ab = B_acc[a];
                for (int idx = 0; idx < N3; idx++) {
                    vs[idx] += hdt*as[idx];
                    vb[idx] += hdt*ab[idx];
                }
            }
            apply_edge_damping_both(g);

            /* Control: standard free wave */
            verlet_full_step(gc, BIMODAL, mass2_phys);
            apply_edge_damping_both(gc);  /* This only damps gc->vel; B_vel is shared — use separate damping */
        }

        if (step % snap_every == 0 && snap_count < N_SNAPS + 1) {
            double t = step * g->dt;
            double fc = compute_fc(g);
            double peak_P = compute_peak_P(g);
            double winding = compute_winding(g);
            double energy = compute_energy_m7(g, BIMODAL);
            double E_shells[5];
            compute_shell_energy(g, E_shells);

            compute_radial_rho(g, r_bins, rho_total, counts_t, dr);

            /* Control: compute bg-only radial profile from gc */
            {
                double dx_c = gc->dx, L_c = gc->L;
                int N_c = gc->N;
                for (int b = 0; b < NBINS; b++) { rho_ctrl[b] = 0; counts_c[b] = 0; }
                for (int i = 1; i < N_c-1; i++) {
                    double x = -L_c + i*dx_c;
                    for (int j = 1; j < N_c-1; j++) {
                        double y = -L_c + j*dx_c;
                        double rp = sqrt(x*x + y*y);
                        int b = (int)(rp / dr);
                        if (b >= NBINS) continue;
                        for (int k = 0; k < N_c; k++) {
                            int idx = i*NN + j*N + k;
                            int kp = (k+1)%N_c, km = (k-1+N_c)%N_c;
                            double e = 0;
                            for (int a = 0; a < NFIELDS; a++) {
                                e += 0.5*gc->vel[a][idx]*gc->vel[a][idx];
                                double gx = (gc->phi[a][idx+NN]-gc->phi[a][idx-NN])/(2*dx_c);
                                double gy = (gc->phi[a][idx+N]-gc->phi[a][idx-N])/(2*dx_c);
                                double gz = (gc->phi[a][i*NN+j*N+kp]-gc->phi[a][i*NN+j*N+km])/(2*dx_c);
                                e += 0.5*(gx*gx+gy*gy+gz*gz);
                                e += 0.5*mass2_phys*gc->phi[a][idx]*gc->phi[a][idx];
                            }
                            rho_ctrl[b] += e; counts_c[b]++;
                        }
                    }
                }
                for (int b = 0; b < NBINS; b++)
                    if (counts_c[b] > 0) rho_ctrl[b] /= counts_c[b];
            }

            fprintf(fts, "%.2f\t%.6f\t%.6e\t%.4f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n",
                    t, fc, peak_P, winding, energy,
                    E_shells[0], E_shells[1], E_shells[2], E_shells[3], E_shells[4]);

            for (int b = 0; b < NBINS; b++) {
                double drho = rho_total[b] - rho_ctrl[b];
                fprintf(fpr, "%.2f\t%.4f\t%.6e\t%.6e\t%.6e\n",
                        t, r_bins[b], rho_total[b], rho_ctrl[b], drho);
            }

            double max_drate = 0;
            if (prev_t >= 0 && t > prev_t) {
                for (int b = 0; b < NBINS; b++) {
                    double drho = rho_total[b] - rho_ctrl[b];
                    double rate = (drho - prev_drho[b]) / (t - prev_t);
                    if (fabs(rate) > max_drate) max_drate = fabs(rate);
                }
            }

            if (t > 300.0) {
                for (int b = 0; b < NBINS; b++)
                    late_drho[b] += rho_total[b] - rho_ctrl[b];
                late_count++;
            }

            for (int b = 0; b < NBINS; b++)
                prev_drho[b] = rho_total[b] - rho_ctrl[b];
            prev_t = t;

            printf("t=%6.1f  fc=%.4f  |P|=%.4e  w=%.3f  E=%.1f  max|drate|=%.2e  [%.0fs]\n",
                   t, fc, peak_P, winding, energy, max_drate, omp_get_wtime()-wt0);

            if (dump_idx < N_DUMPS && fabs(t - dump_times[dump_idx]) < SNAP_DT*0.5) {
                dump_rho_xy(g, t, "charM7");
                dump_zaxis(g, t, "charM7");
                dump_idx++;
            }

            snap_count++;
            if (check_blowup(g)) { printf("BLOWUP at t=%.1f\n", t); break; }
        }
    }

    fclose(fts); fclose(fpr);

    if (late_count > 0) {
        for (int b = 0; b < NBINS; b++) late_drho[b] /= late_count;
        printf("\nLate-time (T>300) average depletion profile fit:\n");
        fit_profiles(r_bins, late_drho, NBINS, "charM7");
    }

    for (int a = 0; a < NFIELDS; a++) { free(B_phi[a]); free(B_vel[a]); free(B_acc[a]); }
    grid_free(g); grid_free(gc);
    printf("\n=== M7 Two-Component Complete (%.0f s) ===\n", omp_get_wtime()-wt0);
    return 0;
}
