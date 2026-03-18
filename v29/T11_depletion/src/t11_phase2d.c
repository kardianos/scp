/*  t11_phase2d.c — Decompose depletion into field-level effects
 *
 *  The question: does the braid CONSUME the background field (amplitude reduction)
 *  or does it just add negative potential energy (V < 0 from mu < 0)?
 *
 *  Method: compute phi² contributions separately.
 *  For the background-only run: phi_bg² at each r
 *  For the braid+bg run: phi_both² at each r
 *  If the braid consumes background: phi_both² should be LESS than phi_braid² + phi_bg²
 *  in the kinetic energy and mass terms (not just potential).
 *
 *  Also: measure the BACKGROUND AMPLITUDE specifically.
 *  The bg is cos(k*z), so project phi onto this mode: 
 *    A_eff(r) = (2/Nz) * sum_z phi(r,z) * cos(k*z)
 *  If A_eff(r) < A_bg at the braid location: field IS consumed.
 *
 *  Build: gcc -O3 -fopenmp -o t11p2d src/t11_phase2d.c -lm
 */

#include "../../src/braid_core.h"

int main(void) {
    bimodal_init_params();
    double mass2 = BIMODAL[14] * BIMODAL[14];
    int N = 128;
    double L = 20.0;
    double A_bg = 0.1;
    double k_bg = PI / L;
    double omega_bg = sqrt(k_bg * k_bg + mass2);

    printf("=== T11 Phase 2d: Field Amplitude Decomposition ===\n");
    printf("N=%d L=%.1f A_bg=%.3f k=%.4f omega=%.4f\n", N, L, A_bg, k_bg, omega_bg);

    double dx = 2.0 * L / (N - 1);
    int NN = N*N;

    /* Three grids */
    Grid *g_braid = grid_alloc(N, L);
    Grid *g_bg    = grid_alloc(N, L);
    Grid *g_both  = grid_alloc(N, L);

    init_braid(g_braid, BIMODAL, -1);
    grid_zero(g_bg);
    init_braid(g_both, BIMODAL, -1);

    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            for (int k = 0; k < N; k++) {
                int idx = i*NN + j*N + k;
                double z = -L + k * dx;
                for (int a = 0; a < NFIELDS; a++) {
                    g_bg->phi[a][idx]    = A_bg * cos(k_bg * z);
                    g_bg->vel[a][idx]    = A_bg * omega_bg * sin(k_bg * z);
                    g_both->phi[a][idx] += A_bg * cos(k_bg * z);
                    g_both->vel[a][idx] += A_bg * omega_bg * sin(k_bg * z);
                }
            }

    compute_forces(g_braid, BIMODAL, -1);
    compute_forces(g_bg,    BIMODAL, -1);
    compute_forces(g_both,  BIMODAL, -1);

    /* Evolve for T=15 (no damping) */
    double T = 10.0;
    int n_steps = (int)(T / g_braid->dt);
    printf("Evolving %d steps (T=%.1f)...\n", n_steps, T);
    for (int step = 0; step < n_steps; step++) {
        verlet_full_step(g_braid, BIMODAL, -1);
        verlet_full_step(g_bg,    BIMODAL, -1);
        verlet_full_step(g_both,  BIMODAL, -1);
    }

    /* Project onto background mode: A_eff(r) */
    int nbins = 40;
    double dr = L / nbins;
    double *A_eff_both  = calloc(nbins, sizeof(double));
    double *A_eff_bg    = calloc(nbins, sizeof(double));
    double *A_eff_braid = calloc(nbins, sizeof(double));
    int    *counts      = calloc(nbins, sizeof(int));

    /* Also compute per-component energy densities */
    double *rho_ke_both  = calloc(nbins, sizeof(double));
    double *rho_ke_sum   = calloc(nbins, sizeof(double));
    double *rho_mass_both = calloc(nbins, sizeof(double));
    double *rho_mass_sum  = calloc(nbins, sizeof(double));
    double *rho_pot_both  = calloc(nbins, sizeof(double));
    double *rho_pot_sum   = calloc(nbins, sizeof(double));

    double mu = BIMODAL[12], kappa = BIMODAL[13];

    for (int i = 1; i < N-1; i++) {
        double x = -L + i * dx;
        for (int j = 1; j < N-1; j++) {
            double y = -L + j * dx;
            double rp = sqrt(x*x + y*y);
            int b = (int)(rp / dr);
            if (b >= nbins) continue;

            /* Fourier projection: A_eff = (2/Nz) * sum_k phi(k) * cos(k_bg * z) */
            double proj_both = 0, proj_bg = 0, proj_braid = 0;
            double ek_both = 0, ek_braid = 0, ek_bg = 0;
            double em_both = 0, em_braid = 0, em_bg = 0;
            double ep_both = 0, ep_braid = 0, ep_bg = 0;

            for (int k = 0; k < N; k++) {
                int idx = i*NN + j*N + k;
                double z = -L + k * dx;
                double ck = cos(k_bg * z);

                /* Project field 0 onto cos(k*z) */
                proj_both  += g_both->phi[0][idx]  * ck;
                proj_bg    += g_bg->phi[0][idx]    * ck;
                proj_braid += g_braid->phi[0][idx] * ck;

                /* KE */
                for (int a = 0; a < NFIELDS; a++) {
                    ek_both  += 0.5 * g_both->vel[a][idx]  * g_both->vel[a][idx];
                    ek_braid += 0.5 * g_braid->vel[a][idx] * g_braid->vel[a][idx];
                    ek_bg    += 0.5 * g_bg->vel[a][idx]    * g_bg->vel[a][idx];
                    em_both  += 0.5 * mass2 * g_both->phi[a][idx]  * g_both->phi[a][idx];
                    em_braid += 0.5 * mass2 * g_braid->phi[a][idx] * g_braid->phi[a][idx];
                    em_bg    += 0.5 * mass2 * g_bg->phi[a][idx]    * g_bg->phi[a][idx];
                }

                /* Potential */
                double P_b  = g_both->phi[0][idx]  * g_both->phi[1][idx]  * g_both->phi[2][idx];
                double P_br = g_braid->phi[0][idx] * g_braid->phi[1][idx] * g_braid->phi[2][idx];
                double P_bg = g_bg->phi[0][idx]    * g_bg->phi[1][idx]    * g_bg->phi[2][idx];
                ep_both  += (mu/2.0) * P_b*P_b   / (1.0 + kappa*P_b*P_b);
                ep_braid += (mu/2.0) * P_br*P_br / (1.0 + kappa*P_br*P_br);
                ep_bg    += (mu/2.0) * P_bg*P_bg / (1.0 + kappa*P_bg*P_bg);
            }

            A_eff_both[b]  += fabs(proj_both)  * 2.0 / N;
            A_eff_bg[b]    += fabs(proj_bg)    * 2.0 / N;
            A_eff_braid[b] += fabs(proj_braid) * 2.0 / N;

            rho_ke_both[b]  += ek_both / N;
            rho_ke_sum[b]   += (ek_braid + ek_bg) / N;
            rho_mass_both[b] += em_both / N;
            rho_mass_sum[b]  += (em_braid + em_bg) / N;
            rho_pot_both[b]  += ep_both / N;
            rho_pot_sum[b]   += (ep_braid + ep_bg) / N;
            counts[b]++;
        }
    }

    /* Normalize */
    for (int b = 0; b < nbins; b++) {
        if (counts[b] > 0) {
            A_eff_both[b]  /= counts[b];
            A_eff_bg[b]    /= counts[b];
            A_eff_braid[b] /= counts[b];
            rho_ke_both[b]  /= counts[b];
            rho_ke_sum[b]   /= counts[b];
            rho_mass_both[b] /= counts[b];
            rho_mass_sum[b]  /= counts[b];
            rho_pot_both[b]  /= counts[b];
            rho_pot_sum[b]   /= counts[b];
        }
    }

    /* Output */
    printf("\n=== BACKGROUND MODE AMPLITUDE ===\n");
    printf("  If A_eff_both < A_eff_bg at the braid core: FIELD IS CONSUMED\n");
    printf("  r     A_both   A_bg     A_braid  A_both-A_bg-A_braid\n");
    for (int b = 0; b < nbins && (b+0.5)*dr < 0.8*L; b++) {
        if (b % 2 != 0) continue;
        double r = (b + 0.5) * dr;
        double excess = A_eff_both[b] - A_eff_bg[b] - A_eff_braid[b];
        printf("  %4.1f  %.5f  %.5f  %.5f  %+.6f\n",
               r, A_eff_both[b], A_eff_bg[b], A_eff_braid[b], excess);
    }

    printf("\n=== ENERGY COMPONENT DECOMPOSITION ===\n");
    printf("  delta_KE = KE_both - (KE_braid + KE_bg)\n");
    printf("  delta_mass = mass_both - (mass_braid + mass_bg)\n");
    printf("  delta_pot = pot_both - (pot_braid + pot_bg)\n");
    printf("  r     delta_KE       delta_mass     delta_pot\n");
    for (int b = 0; b < nbins && (b+0.5)*dr < 0.8*L; b++) {
        if (b % 2 != 0) continue;
        double r = (b + 0.5) * dr;
        printf("  %4.1f  %+13.6e  %+13.6e  %+13.6e\n",
               r, rho_ke_both[b] - rho_ke_sum[b],
               rho_mass_both[b] - rho_mass_sum[b],
               rho_pot_both[b] - rho_pot_sum[b]);
    }

    printf("\n=== CONCLUSION ===\n");
    /* Check which component dominates the interaction */
    double total_dke = 0, total_dmass = 0, total_dpot = 0;
    for (int b = 0; b < nbins; b++) {
        total_dke   += (rho_ke_both[b] - rho_ke_sum[b]) * counts[b];
        total_dmass += (rho_mass_both[b] - rho_mass_sum[b]) * counts[b];
        total_dpot  += (rho_pot_both[b] - rho_pot_sum[b]) * counts[b];
    }
    printf("Integrated interaction energy components:\n");
    printf("  delta_KE   = %+.4f\n", total_dke);
    printf("  delta_mass = %+.4f\n", total_dmass);
    printf("  delta_pot  = %+.4f\n", total_dpot);
    printf("  TOTAL      = %+.4f\n", total_dke + total_dmass + total_dpot);

    if (total_dpot < 0 && fabs(total_dpot) > fabs(total_dke) + fabs(total_dmass)) {
        printf("\nDominated by POTENTIAL: mu < 0 → V < 0 → nonlinear binding.\n");
        printf("This is NOT field consumption. It's the triple-product attraction.\n");
        printf("The background field is not consumed — the potential is just lowered.\n");
    } else if (total_dke + total_dmass < 0) {
        printf("\nKE + mass contribution is NEGATIVE → background field IS reduced!\n");
        printf("This would be genuine field consumption (depletion).\n");
    } else {
        printf("\nKE + mass contribution is positive → no field consumption.\n");
        printf("The negative interaction energy is purely from V(P) < 0.\n");
    }

    /* Cleanup */
    free(A_eff_both); free(A_eff_bg); free(A_eff_braid); free(counts);
    free(rho_ke_both); free(rho_ke_sum);
    free(rho_mass_both); free(rho_mass_sum);
    free(rho_pot_both); free(rho_pot_sum);
    grid_free(g_braid); grid_free(g_bg); grid_free(g_both);

    return 0;
}
