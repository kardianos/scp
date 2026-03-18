/*  t12_m1.c — Mechanism 1: Density-dependent sink
 *  After each Verlet step, scale velocities in braid core (r<5)
 *  by factor (1 - alpha*dt*rho_local/rho0). Extracts KE proportional to local density.
 *
 *  Build: gcc -O3 -fopenmp -o t12_m1 src/t12_m1.c -lm
 */

#include "../../src/braid_core.h"

/* ================================================================
   Shared measurement infrastructure
   ================================================================ */

#define NBINS 60
#define A_BG  0.1
#define T_TOTAL 300.0
#define SNAP_DT 50.0

static double rho0_bg;   /* expected background energy density */
static double mass2_phys;

static void compute_radial_rho(Grid *g, double *r_bins, double *rho_bins,
                                int *counts, double dr) {
    int N = g->N, NN = N*N;
    double dx = g->dx, L = g->L;

    for (int b = 0; b < NBINS; b++) {
        r_bins[b] = (b + 0.5) * dr;
        rho_bins[b] = 0;
        counts[b] = 0;
    }

    for (int i = 1; i < N-1; i++) {
        double x = -L + i * dx;
        for (int j = 1; j < N-1; j++) {
            double y = -L + j * dx;
            double rp = sqrt(x*x + y*y);
            int b = (int)(rp / dr);
            if (b >= NBINS) continue;

            for (int k = 0; k < N; k++) {
                int idx = i*NN + j*N + k;
                int kp = (k+1)%N, km = (k-1+N)%N;

                double ek = 0, eg = 0, em = 0;
                for (int a = 0; a < NFIELDS; a++) {
                    ek += 0.5 * g->vel[a][idx] * g->vel[a][idx];
                    double gx = (g->phi[a][idx+NN] - g->phi[a][idx-NN]) / (2*dx);
                    double gy = (g->phi[a][idx+N]  - g->phi[a][idx-N])  / (2*dx);
                    double gz = (g->phi[a][i*NN+j*N+kp] - g->phi[a][i*NN+j*N+km]) / (2*dx);
                    eg += 0.5 * (gx*gx + gy*gy + gz*gz);
                    em += 0.5 * mass2_phys * g->phi[a][idx] * g->phi[a][idx];
                }
                rho_bins[b] += ek + eg + em;
                counts[b]++;
            }
        }
    }

    for (int b = 0; b < NBINS; b++)
        if (counts[b] > 0) rho_bins[b] /= counts[b];
}

static double get_rho_at_r(double *r_bins, double *rho_bins, int *counts, double r_target) {
    double best = 0;
    double best_dist = 1e30;
    for (int b = 0; b < NBINS; b++) {
        if (counts[b] < 5) continue;
        double d = fabs(r_bins[b] - r_target);
        if (d < best_dist) { best_dist = d; best = rho_bins[b]; }
    }
    return best;
}

static void add_background(Grid *g) {
    int N = g->N, NN = N*N;
    double dx = g->dx, L = g->L;
    double k_bg = PI / L;
    double omega_bg = sqrt(k_bg * k_bg + mass2_phys);

    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            for (int k = 0; k < N; k++) {
                int idx = i*NN + j*N + k;
                double z = -L + k * dx;
                for (int a = 0; a < NFIELDS; a++) {
                    double phase = k_bg * z + 2.0*PI*a/3.0;
                    g->phi[a][idx] += A_BG * cos(phase);
                    g->vel[a][idx] += omega_bg * A_BG * sin(phase);
                }
            }
}

static void apply_edge_damping(Grid *g) {
    int N = g->N, NN = N * N;
    double dx = g->dx, L = g->L;
    double r_start = 0.85 * L, r_end = 0.98 * L;
    double inv_dr = 1.0 / (r_end - r_start + 1e-30);

    for (int i = 0; i < N; i++) {
        double x = -L + i * dx;
        for (int j = 0; j < N; j++) {
            double y = -L + j * dx;
            double rp = sqrt(x*x + y*y);
            if (rp <= r_start) continue;
            double f = (rp - r_start) * inv_dr;
            if (f > 1.0) f = 1.0;
            double damp = 1.0 - 0.95 * f * f;
            for (int kk = 0; kk < N; kk++) {
                int idx = i * NN + j * N + kk;
                for (int a = 0; a < NFIELDS; a++) {
                    g->vel[a][idx] *= damp;
                }
            }
        }
    }
}

static double compute_fc(Grid *g) {
    int N = g->N, NN = N*N;
    double dx = g->dx, L = g->L;
    double R_core = 8.0, Rc2 = R_core * R_core;
    double total = 0, core = 0;

    for (int i = 1; i < N-1; i++) {
        double x = -L + i * dx;
        for (int j = 1; j < N-1; j++) {
            double y = -L + j * dx;
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
    double dx = g->dx, L = g->L, dV = dx*dx*dx;
    double mu = phys[12], kappa = phys[13], lpw = phys[15];
    double E = 0;

    for (int i = 1; i < N-1; i++)
        for (int j = 1; j < N-1; j++)
            for (int k = 0; k < N; k++) {
                int idx = i*NN + j*N + k;
                int kp = (k+1)%N, km = (k-1+N)%N;
                double ek = 0, eg = 0, em = 0;
                for (int a = 0; a < NFIELDS; a++) {
                    ek += 0.5 * g->vel[a][idx] * g->vel[a][idx];
                    double gx = (g->phi[a][idx+NN] - g->phi[a][idx-NN]) / (2*dx);
                    double gy = (g->phi[a][idx+N]  - g->phi[a][idx-N])  / (2*dx);
                    double gz = (g->phi[a][i*NN+j*N+kp] - g->phi[a][i*NN+j*N+km]) / (2*dx);
                    eg += 0.5*(gx*gx + gy*gy + gz*gz);
                    em += 0.5 * mass2_phys * g->phi[a][idx] * g->phi[a][idx];
                }
                double p0 = g->phi[0][idx], p1 = g->phi[1][idx], p2 = g->phi[2][idx];
                double P = p0*p1*p2;
                double ep = (mu/2.0)*P*P/(1.0+kappa*P*P);
                double epw = lpw*(p0*p1 + p1*p2 + p2*p0);
                E += (ek+eg+em+ep+epw) * dV;
            }
    return E;
}

/* ================================================================
   M1-specific: density-dependent sink
   ================================================================ */

static void apply_m1_sink(Grid *g, double alpha) {
    int N = g->N, NN = N*N;
    double dx = g->dx, L = g->L, dt = g->dt;
    double R_sink = 5.0, R_sink2 = R_sink * R_sink;

    for (int i = 1; i < N-1; i++) {
        double x = -L + i * dx;
        for (int j = 1; j < N-1; j++) {
            double y = -L + j * dx;
            double rp2 = x*x + y*y;
            if (rp2 > R_sink2) continue;

            for (int k = 0; k < N; k++) {
                int idx = i*NN + j*N + k;
                /* local rho ~ sum(v^2 + m^2 phi^2)/2 */
                double rho_local = 0;
                for (int a = 0; a < NFIELDS; a++)
                    rho_local += 0.5*g->vel[a][idx]*g->vel[a][idx]
                               + 0.5*mass2_phys*g->phi[a][idx]*g->phi[a][idx];

                double factor = 1.0 - alpha * dt * rho_local / (rho0_bg + 1e-30);
                if (factor < 0.9) factor = 0.9;  /* clamp */
                for (int a = 0; a < NFIELDS; a++)
                    g->vel[a][idx] *= factor;
            }
        }
    }
}

int main(void) {
    bimodal_init_params();
    mass2_phys = 2.25;  /* from BIMODAL mass=1.5 */
    double k_bg = PI / 30.0;
    double omega_bg = sqrt(k_bg * k_bg + mass2_phys);
    rho0_bg = 3.0 * A_BG * A_BG * omega_bg * omega_bg;

    printf("=== T12 M1: Density-Dependent Sink ===\n");
    printf("mass2=%.4f, rho0_bg=%.6f, omega_bg=%.4f\n", mass2_phys, rho0_bg, omega_bg);

    int N = 96;
    double L = 30.0;
    Grid *g = grid_alloc(N, L);
    init_braid(g, BIMODAL, -1);
    add_background(g);
    compute_forces(g, BIMODAL, mass2_phys);

    double dr = L / NBINS;
    double r_bins[NBINS], rho_bins[NBINS];
    int counts[NBINS];

    FILE *fout = fopen("data/t12_m1_timeseries.dat", "w");
    fprintf(fout, "# t  fc  energy  drho_r10  drho_r15  drho_r20\n");

    /* Store previous depletion for rate */
    double prev_drho15 = 0;
    double prev_t = 0;

    int n_total = (int)(T_TOTAL / g->dt);
    int snap_every = (int)(SNAP_DT / g->dt);
    double alpha = 0.01;

    printf("dt=%.5f, n_total=%d, snap_every=%d\n", g->dt, n_total, snap_every);

    for (int step = 0; step <= n_total; step++) {
        if (step > 0) {
            verlet_full_step(g, BIMODAL, mass2_phys);
            apply_m1_sink(g, alpha);
            apply_edge_damping(g);
        }

        if (step % snap_every == 0) {
            double t = step * g->dt;
            compute_radial_rho(g, r_bins, rho_bins, counts, dr);
            double rho10 = get_rho_at_r(r_bins, rho_bins, counts, 10.0);
            double rho15 = get_rho_at_r(r_bins, rho_bins, counts, 15.0);
            double rho20 = get_rho_at_r(r_bins, rho_bins, counts, 20.0);
            double drho10 = rho10 - rho0_bg;
            double drho15 = rho15 - rho0_bg;
            double drho20 = rho20 - rho0_bg;
            double fc = compute_fc(g);
            double E = compute_energy(g, BIMODAL);
            double rate15 = (t > prev_t + 1.0) ? (drho15 - prev_drho15)/(t - prev_t) : 0;

            printf("t=%6.1f  fc=%.4f  E=%.1f  drho(10)=%+.4e  drho(15)=%+.4e  drho(20)=%+.4e  d/dt=%.2e\n",
                   t, fc, E, drho10, drho15, drho20, rate15);
            fprintf(fout, "%.2f  %.6f  %.2f  %.8e  %.8e  %.8e\n",
                    t, fc, E, drho10, drho15, drho20);

            prev_drho15 = drho15;
            prev_t = t;

            if (check_blowup(g)) { printf("BLOWUP at t=%.1f\n", t); break; }
        }
    }

    fclose(fout);
    grid_free(g);
    printf("=== M1 Complete ===\n");
    return 0;
}
