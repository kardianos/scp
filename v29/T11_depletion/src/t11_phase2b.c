/*  t11_phase2b.c — Braid depletion of background: differential measurement
 *
 *  Problem with Phase 2: absorbing BC kills the background wave too.
 *
 *  Better approach: run TWO simulations in parallel:
 *    Run A: background only (control)
 *    Run B: braid + background
 *  Then delta_rho(r) = rho_B(r) - rho_A(r) is the braid's effect.
 *  If delta_rho < 0 somewhere outside core: DEPLETION exists.
 *
 *  Build: gcc -O3 -fopenmp -o t11p2b src/t11_phase2b.c -lm
 */

#include "../../src/braid_core.h"
#include <float.h>

/* Radial profile */
static void compute_radial_profile(Grid *g, const double *phys, double mass2,
                                    double *rho_bins, int *counts, int nbins, double dr) {
    int N = g->N, NN = N*N;
    double dx = g->dx, L = g->L;
    double mu = phys[12], kappa = phys[13];

    for (int b = 0; b < nbins; b++) { rho_bins[b] = 0; counts[b] = 0; }

    for (int i = 1; i < N-1; i++) {
        double x = -L + i * dx;
        for (int j = 1; j < N-1; j++) {
            double y = -L + j * dx;
            double rp = sqrt(x*x + y*y);
            int b = (int)(rp / dr);
            if (b >= nbins) continue;

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
                    em += 0.5 * mass2 * g->phi[a][idx] * g->phi[a][idx];
                }
                double p0 = g->phi[0][idx], p1 = g->phi[1][idx], p2 = g->phi[2][idx];
                double P = p0 * p1 * p2;
                double ep = (mu/2.0) * P*P / (1.0 + kappa*P*P);
                rho_bins[b] += ek + eg + em + ep;
                counts[b]++;
            }
        }
    }
    for (int b = 0; b < nbins; b++)
        if (counts[b] > 0) rho_bins[b] /= counts[b];
}

/* Mild damping — only far corners */
static void apply_mild_damping(Grid *g) {
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
            double damp = 1.0 - 0.5 * f * f;
            for (int kk = 0; kk < N; kk++) {
                int idx = i * NN + j * N + kk;
                for (int a = 0; a < NFIELDS; a++) {
                    g->phi[a][idx] *= damp;
                    g->vel[a][idx] *= damp;
                }
            }
        }
    }
}

int main(int argc, char **argv) {
    bimodal_init_params();
    double mass2 = BIMODAL[14] * BIMODAL[14];

    int N = 128;
    double L = 30.0;
    double A_bg = 0.1;
    double T_total = 400.0;

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-N") == 0 && i+1 < argc) N = atoi(argv[++i]);
        else if (strcmp(argv[i], "-L") == 0 && i+1 < argc) L = atof(argv[++i]);
        else if (strcmp(argv[i], "-bg") == 0 && i+1 < argc) A_bg = atof(argv[++i]);
        else if (strcmp(argv[i], "-T") == 0 && i+1 < argc) T_total = atof(argv[++i]);
    }

    printf("=== T11 Phase 2b: Differential Depletion ===\n");
    printf("N=%d  L=%.1f  A_bg=%.3f  T=%.1f\n", N, L, A_bg, T_total);

    double dx = 2.0 * L / (N - 1);
    double dt = 0.20 * dx;
    printf("dx=%.4f  dt=%.5f  mass²=%.4f\n", dx, dt, mass2);

    double k_bg = PI / L;
    double omega_bg = sqrt(k_bg * k_bg + mass2);
    double rho_bg_theory = NFIELDS * A_bg * A_bg * omega_bg * omega_bg;
    printf("Background: k=%.4f omega=%.4f rho_theory=%.6e\n", k_bg, omega_bg, rho_bg_theory);

    int nbins = 60;
    double dr = L / nbins;

    /* Allocate both grids */
    Grid *g_ctrl  = grid_alloc(N, L);
    Grid *g_braid = grid_alloc(N, L);

    /* Init control: background only */
    grid_zero(g_ctrl);
    int NN = N*N;
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            for (int k = 0; k < N; k++) {
                int idx = i*NN + j*N + k;
                double z = -L + k * dx;
                for (int a = 0; a < NFIELDS; a++) {
                    g_ctrl->phi[a][idx] = A_bg * cos(k_bg * z);
                    g_ctrl->vel[a][idx] = A_bg * omega_bg * sin(k_bg * z);
                }
            }
    compute_forces(g_ctrl, BIMODAL, -1);

    /* Init braid: braid + background */
    init_braid(g_braid, BIMODAL, -1);
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            for (int k = 0; k < N; k++) {
                int idx = i*NN + j*N + k;
                double z = -L + k * dx;
                for (int a = 0; a < NFIELDS; a++) {
                    g_braid->phi[a][idx] += A_bg * cos(k_bg * z);
                    g_braid->vel[a][idx] += A_bg * omega_bg * sin(k_bg * z);
                }
            }
    compute_forces(g_braid, BIMODAL, -1);

    /* Profiles */
    double *rho_ctrl  = calloc(nbins, sizeof(double));
    double *rho_braid = calloc(nbins, sizeof(double));
    int    *counts    = calloc(nbins, sizeof(int));

    /* Accumulators for time averaging */
    double *avg_ctrl  = calloc(nbins, sizeof(double));
    double *avg_braid = calloc(nbins, sizeof(double));
    int n_avg = 0;

    int n_steps = (int)(T_total / dt);
    int snap_every = n_steps / 20;
    int avg_start = n_steps / 2;

    /* Time series output */
    FILE *fts = fopen("data/phase2b_timeseries.tsv", "w");
    fprintf(fts, "t\trho_braid_r0\trho_ctrl_r0\trho_braid_r5\trho_ctrl_r5\trho_braid_r10\trho_ctrl_r10\trho_braid_r15\trho_ctrl_r15\n");

    printf("\nRunning both simulations (%d steps)...\n", n_steps);

    for (int step = 0; step <= n_steps; step++) {
        if (step > 0) {
            verlet_full_step(g_ctrl,  BIMODAL, -1);
            verlet_full_step(g_braid, BIMODAL, -1);
            apply_mild_damping(g_ctrl);
            apply_mild_damping(g_braid);
        }

        if (step % snap_every == 0) {
            double t = step * dt;
            compute_radial_profile(g_ctrl,  BIMODAL, mass2, rho_ctrl,  counts, nbins, dr);
            compute_radial_profile(g_braid, BIMODAL, mass2, rho_braid, counts, nbins, dr);

            int b0 = 0, b5 = (int)(5.0/dr), b10 = (int)(10.0/dr), b15 = (int)(15.0/dr);
            if (b15 >= nbins) b15 = nbins - 1;

            printf("  t=%6.1f  B-C at r=5: %+.4e  r=10: %+.4e  r=15: %+.4e\n",
                   t,
                   rho_braid[b5] - rho_ctrl[b5],
                   rho_braid[b10] - rho_ctrl[b10],
                   rho_braid[b15] - rho_ctrl[b15]);

            fprintf(fts, "%.2f\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\n",
                    t, rho_braid[b0], rho_ctrl[b0],
                    rho_braid[b5], rho_ctrl[b5],
                    rho_braid[b10], rho_ctrl[b10],
                    rho_braid[b15], rho_ctrl[b15]);
            fflush(fts);

            if (check_blowup(g_braid)) { printf("BRAID BLOWUP!\n"); break; }
            if (check_blowup(g_ctrl))  { printf("CTRL BLOWUP!\n"); break; }
        }

        /* Accumulate for time averaging */
        if (step >= avg_start && step % (snap_every/2 + 1) == 0) {
            compute_radial_profile(g_ctrl,  BIMODAL, mass2, rho_ctrl,  counts, nbins, dr);
            compute_radial_profile(g_braid, BIMODAL, mass2, rho_braid, counts, nbins, dr);
            for (int b = 0; b < nbins; b++) {
                avg_ctrl[b]  += rho_ctrl[b];
                avg_braid[b] += rho_braid[b];
            }
            n_avg++;
        }
    }
    fclose(fts);

    /* Finalize averages */
    if (n_avg > 0) {
        for (int b = 0; b < nbins; b++) {
            avg_ctrl[b]  /= n_avg;
            avg_braid[b] /= n_avg;
        }
    }

    /* Save differential profile */
    FILE *fp = fopen("data/phase2b_profiles.tsv", "w");
    fprintf(fp, "r\trho_braid\trho_ctrl\tdelta_rho\tdelta_frac\n");

    printf("\n=== DIFFERENTIAL PROFILE (time-averaged, T=%.0f..%.0f) ===\n", T_total/2, T_total);
    printf("  r       rho_braid      rho_ctrl       delta          delta/ctrl\n");
    printf("  -----  -------------  -------------  -------------  ----------\n");

    int found_depletion = 0;
    double max_depl_frac = 0, r_max_depl = 0;

    for (int b = 0; b < nbins; b++) {
        double r = (b + 0.5) * dr;
        if (r > 0.80 * L) continue;
        double delta = avg_braid[b] - avg_ctrl[b];
        double frac = delta / (fabs(avg_ctrl[b]) + 1e-30);

        fprintf(fp, "%.4f\t%.8e\t%.8e\t%+.8e\t%+.8f\n",
                r, avg_braid[b], avg_ctrl[b], delta, frac);

        if (b % 3 == 0) {
            printf("  %5.2f  %13.6e  %13.6e  %+13.6e  %+8.4f\n",
                   r, avg_braid[b], avg_ctrl[b], delta, frac);
        }

        /* Look for depletion outside core (r > 5) */
        if (r > 5.0 && r < 0.7*L && delta < 0) {
            if (frac < max_depl_frac) {
                max_depl_frac = frac;
                r_max_depl = r;
            }
            found_depletion = 1;
        }
    }
    fclose(fp);

    printf("\n=== VERDICT ===\n");
    if (found_depletion) {
        printf("DEPLETION DETECTED at r=%.2f: delta/rho_ctrl = %.4f (%.2f%% reduction)\n",
               r_max_depl, max_depl_frac, -100.0 * max_depl_frac);
        printf("The braid REDUCES the background field energy density.\n");
    } else {
        printf("NO DEPLETION: braid adds energy (delta > 0) everywhere outside core.\n");
        printf("The braid is a pure energy source, not a sink.\n");
    }

    /* Energy analysis */
    printf("\n=== ENERGY BALANCE ===\n");
    Result res;
    snprintf(res.label, sizeof(res.label), "final");
    compute_diagnostics(g_braid, BIMODAL, &res);
    double E_braid = res.energy;
    compute_diagnostics(g_ctrl, BIMODAL, &res);
    double E_ctrl = res.energy;
    printf("E_braid = %.2f  E_ctrl = %.2f  delta_E = %.2f\n",
           E_braid, E_ctrl, E_braid - E_ctrl);
    printf("If delta_E > 0: braid adds net energy (no global depletion)\n");
    printf("If delta_E < 0: braid consumes net energy (global depletion)\n");

    /* Potential energy component */
    printf("\nmu = %.2f (negative → V < 0 where P ≠ 0)\n", BIMODAL[12]);
    printf("The triple-product potential is ATTRACTIVE (mu < 0).\n");
    printf("V(P) < 0 at the braid core, but this doesn't deplete the BACKGROUND.\n");
    printf("It just lowers the braid's own energy below the kinetic+gradient+mass.\n");

    free(rho_ctrl); free(rho_braid); free(counts);
    free(avg_ctrl); free(avg_braid);
    grid_free(g_ctrl); grid_free(g_braid);

    printf("\nPhase 2b complete.\n");
    return 0;
}
