/*  t11_phase2c.c — Depletion test: NO DAMPING, smaller box, short time
 *
 *  To eliminate damping artifacts: periodic BC in all directions.
 *  Use short run time T < L so reflections don't contaminate.
 *  Focus on the question: does the braid SUBTRACT from background?
 *
 *  Key insight: if the braid depletes the field, we should see:
 *  rho_braid(r) < rho_ctrl(r) for some r > core_radius
 *  even at early times, before reflections arrive.
 *
 *  Also test: is the braid energy LESS than braid_alone + background_alone?
 *  If E(braid+bg) < E(braid) + E(bg): nonlinear coupling is ATTRACTIVE → depletion
 *  If E(braid+bg) > E(braid) + E(bg): nonlinear coupling REPELS → no depletion
 *
 *  Build: gcc -O3 -fopenmp -o t11p2c src/t11_phase2c.c -lm
 */

#include "../../src/braid_core.h"

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

int main(void) {
    bimodal_init_params();
    double mass2 = BIMODAL[14] * BIMODAL[14];
    int N = 128;
    double L = 20.0;
    
    printf("=== T11 Phase 2c: No-Damping Depletion Test ===\n");
    printf("N=%d L=%.1f mass²=%.4f\n", N, L, mass2);

    double dx = 2.0 * L / (N - 1);
    double dt = 0.20 * dx;
    int NN = N*N;

    /* Test multiple A_bg values */
    double A_vals[] = {0.05, 0.1, 0.3, 0.5};
    int nA = 4;

    for (int ia = 0; ia < nA; ia++) {
        double A_bg = A_vals[ia];
        double k_bg = PI / L;
        double omega_bg = sqrt(k_bg * k_bg + mass2);
        double rho_bg_theory = NFIELDS * A_bg * A_bg * omega_bg * omega_bg;

        printf("\n===== A_bg = %.3f (rho_bg_theory = %.4e) =====\n", A_bg, rho_bg_theory);

        /* Three runs: braid only, bg only, braid+bg */
        Grid *g_braid_only = grid_alloc(N, L);
        Grid *g_bg_only    = grid_alloc(N, L);
        Grid *g_both       = grid_alloc(N, L);

        /* Braid only */
        init_braid(g_braid_only, BIMODAL, -1);

        /* BG only */
        grid_zero(g_bg_only);
        for (int i = 0; i < N; i++)
            for (int j = 0; j < N; j++)
                for (int k = 0; k < N; k++) {
                    int idx = i*NN + j*N + k;
                    double z = -L + k * dx;
                    for (int a = 0; a < NFIELDS; a++) {
                        g_bg_only->phi[a][idx] = A_bg * cos(k_bg * z);
                        g_bg_only->vel[a][idx] = A_bg * omega_bg * sin(k_bg * z);
                    }
                }

        /* Both */
        init_braid(g_both, BIMODAL, -1);
        for (int i = 0; i < N; i++)
            for (int j = 0; j < N; j++)
                for (int k = 0; k < N; k++) {
                    int idx = i*NN + j*N + k;
                    double z = -L + k * dx;
                    for (int a = 0; a < NFIELDS; a++) {
                        g_both->phi[a][idx] += A_bg * cos(k_bg * z);
                        g_both->vel[a][idx] += A_bg * omega_bg * sin(k_bg * z);
                    }
                }

        compute_forces(g_braid_only, BIMODAL, -1);
        compute_forces(g_bg_only,    BIMODAL, -1);
        compute_forces(g_both,       BIMODAL, -1);

        /* Measure initial energies */
        Result res;
        compute_diagnostics(g_braid_only, BIMODAL, &res); double E_braid_0 = res.energy;
        compute_diagnostics(g_bg_only,    BIMODAL, &res); double E_bg_0    = res.energy;
        compute_diagnostics(g_both,       BIMODAL, &res); double E_both_0  = res.energy;

        double E_interaction_0 = E_both_0 - E_braid_0 - E_bg_0;
        printf("  t=0: E_braid=%.1f  E_bg=%.1f  E_both=%.1f\n", E_braid_0, E_bg_0, E_both_0);
        printf("  t=0: E_interaction = E_both - E_braid - E_bg = %.4f\n", E_interaction_0);
        printf("  Sign: %s (interaction %s braid energy)\n",
               E_interaction_0 < 0 ? "NEGATIVE" : "POSITIVE",
               E_interaction_0 < 0 ? "REDUCES" : "INCREASES");

        /* Evolve for short time (T=20, well before box reflections at T=L=20) */
        double T = 15.0;
        int n_steps = (int)(T / dt);
        printf("  Evolving %d steps (T=%.1f, no damping)...\n", n_steps, T);

        for (int step = 0; step < n_steps; step++) {
            verlet_full_step(g_braid_only, BIMODAL, -1);
            verlet_full_step(g_bg_only,    BIMODAL, -1);
            verlet_full_step(g_both,       BIMODAL, -1);
            /* NO damping — all periodic */
        }

        /* Measure final energies (should be conserved for no-damping) */
        compute_diagnostics(g_braid_only, BIMODAL, &res); double E_braid_f = res.energy;
        compute_diagnostics(g_bg_only,    BIMODAL, &res); double E_bg_f    = res.energy;
        compute_diagnostics(g_both,       BIMODAL, &res); double E_both_f  = res.energy;

        double E_interaction_f = E_both_f - E_braid_f - E_bg_f;
        printf("  t=%.0f: E_braid=%.1f  E_bg=%.1f  E_both=%.1f\n", T, E_braid_f, E_bg_f, E_both_f);
        printf("  t=%.0f: E_interaction = %.4f (was %.4f at t=0)\n", T, E_interaction_f, E_interaction_0);

        /* Radial differential profile */
        int nbins = 40;
        double dr = L / nbins;
        double *rho_both  = calloc(nbins, sizeof(double));
        double *rho_bg    = calloc(nbins, sizeof(double));
        double *rho_braid = calloc(nbins, sizeof(double));
        int    *counts    = calloc(nbins, sizeof(int));

        compute_radial_profile(g_both,       BIMODAL, mass2, rho_both,  counts, nbins, dr);
        compute_radial_profile(g_bg_only,    BIMODAL, mass2, rho_bg,    counts, nbins, dr);
        compute_radial_profile(g_braid_only, BIMODAL, mass2, rho_braid, counts, nbins, dr);

        printf("\n  Radial profile at t=%.0f:\n", T);
        printf("  r     rho_both       rho_sum(b+bg)  delta          delta/rho_bg\n");
        int found_depl = 0;
        for (int b = 0; b < nbins; b++) {
            double r = (b + 0.5) * dr;
            if (r > 0.85 * L) continue;
            double rho_sum = rho_braid[b] + rho_bg[b];  /* linear superposition */
            double delta = rho_both[b] - rho_sum;
            double frac = delta / (fabs(rho_bg[b]) + 1e-30);
            if (b % 2 == 0) {
                printf("  %4.1f  %13.6e  %13.6e  %+13.6e  %+8.4f\n",
                       r, rho_both[b], rho_sum, delta, frac);
            }
            if (r > 4.0 && delta < -1e-6 * fabs(rho_bg[b])) found_depl = 1;
        }

        if (found_depl)
            printf("  >>> NONLINEAR DEPLETION: rho_both < rho_braid + rho_bg somewhere\n");
        else
            printf("  >>> NO NONLINEAR DEPLETION: rho_both >= rho_braid + rho_bg\n");

        free(rho_both); free(rho_bg); free(rho_braid); free(counts);
        grid_free(g_braid_only); grid_free(g_bg_only); grid_free(g_both);
    }

    printf("\n=== Phase 2c Complete ===\n");
    return 0;
}
