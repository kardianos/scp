/*  t11_phase1.c — Measure radial energy density profile around a braid
 *
 *  THE KEY QUESTION: does the braid create a region where the energy
 *  density rho(r) < rho_vac OUTSIDE the core?
 *
 *  If yes: the braid DEPLETES the field → depletion gravity is possible.
 *  If no: the braid ADDS energy everywhere → depletion mechanism fails.
 *
 *  Method:
 *  1. Initialize bimodal braid (BIMODAL params, mass2_override=-1)
 *  2. Equilibrate for T=300 with absorbing BC
 *  3. Compute radial energy density profile rho(r_perp):
 *     - Average over z (periodic direction) and azimuthal angle
 *     - Energy density = KE + grad + mass + potential
 *  4. Compare to vacuum energy density at large r
 *  5. Output radial profile to data file
 *
 *  Build: gcc -O3 -fopenmp -o t11p1 src/t11_phase1.c -lm
 */

#include "../../src/braid_core.h"
#include <float.h>

/* ================================================================
   Radial energy density profile
   Average over z and azimuthal angle → rho(r_perp)
   ================================================================ */

static void compute_radial_profile(Grid *g, const double *phys, double mass2,
                                    double *r_bins, double *rho_bins,
                                    double *rho_ke, double *rho_grad,
                                    double *rho_mass, double *rho_pot,
                                    int *counts, int nbins, double dr) {
    int N = g->N, NN = N*N;
    double dx = g->dx, L = g->L;
    double mu = phys[12], kappa = phys[13];

    /* Zero bins */
    for (int b = 0; b < nbins; b++) {
        r_bins[b] = (b + 0.5) * dr;
        rho_bins[b] = 0;
        rho_ke[b] = 0;
        rho_grad[b] = 0;
        rho_mass[b] = 0;
        rho_pot[b] = 0;
        counts[b] = 0;
    }

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

                double ek = 0, eg = 0, em_val = 0;
                for (int a = 0; a < NFIELDS; a++) {
                    /* Kinetic energy density */
                    ek += 0.5 * g->vel[a][idx] * g->vel[a][idx];

                    /* Gradient energy density */
                    double gx = (g->phi[a][idx+NN] - g->phi[a][idx-NN]) / (2*dx);
                    double gy = (g->phi[a][idx+N]  - g->phi[a][idx-N])  / (2*dx);
                    double gz = (g->phi[a][i*NN+j*N+kp] - g->phi[a][i*NN+j*N+km]) / (2*dx);
                    eg += 0.5 * (gx*gx + gy*gy + gz*gz);

                    /* Mass energy density */
                    em_val += 0.5 * mass2 * g->phi[a][idx] * g->phi[a][idx];
                }

                /* Potential energy density */
                double p0 = g->phi[0][idx], p1 = g->phi[1][idx], p2 = g->phi[2][idx];
                double P = p0 * p1 * p2;
                double ep = (mu/2.0) * P*P / (1.0 + kappa*P*P);

                rho_ke[b]   += ek;
                rho_grad[b] += eg;
                rho_mass[b] += em_val;
                rho_pot[b]  += ep;
                rho_bins[b] += ek + eg + em_val + ep;
                counts[b]++;
            }
        }
    }

    /* Normalize by count (average over z and azimuth) */
    for (int b = 0; b < nbins; b++) {
        if (counts[b] > 0) {
            rho_bins[b] /= counts[b];
            rho_ke[b]   /= counts[b];
            rho_grad[b] /= counts[b];
            rho_mass[b] /= counts[b];
            rho_pot[b]  /= counts[b];
        }
    }
}

/* Compute the "vacuum" energy density from the field far from the braid.
   In the vacuum, phi_a oscillate as plane waves: phi_a = A cos(kz - wt + delta).
   The time-averaged energy density is: rho_vac = (omega² + k² + m²) * A² / 2
   = omega² * A² (using dispersion).
   We measure this at the outermost reliable bin. */

int main(int argc, char **argv) {
    bimodal_init_params();

    int N = 128;
    double L = 20.0;
    double T_equil = 300.0;

    printf("=== T11 Phase 1: Field Depletion Profile ===\n");
    printf("N=%d  L=%.1f  T_equil=%.1f\n", N, L, T_equil);
    printf("mass²=%.4f (from BIMODAL[14]=%.4f)\n",
           BIMODAL[14]*BIMODAL[14], BIMODAL[14]);

    /* Allocate and initialize braid */
    Grid *g = grid_alloc(N, L);
    init_braid(g, BIMODAL, -1);  /* mass_init_override=-1 → use BIMODAL[14]=1.5 */
    compute_forces(g, BIMODAL, -1); /* mass2_override=-1 → use BIMODAL[14]²=2.25 */

    /* Measure initial profile before equilibration */
    int nbins = 60;
    double dr = L / nbins;
    double *r_bins   = calloc(nbins, sizeof(double));
    double *rho_bins = calloc(nbins, sizeof(double));
    double *rho_ke   = calloc(nbins, sizeof(double));
    double *rho_grad = calloc(nbins, sizeof(double));
    double *rho_mass = calloc(nbins, sizeof(double));
    double *rho_pot  = calloc(nbins, sizeof(double));
    int    *counts   = calloc(nbins, sizeof(int));

    double mass2 = BIMODAL[14] * BIMODAL[14];

    compute_radial_profile(g, BIMODAL, mass2, r_bins, rho_bins,
                           rho_ke, rho_grad, rho_mass, rho_pot, counts, nbins, dr);

    /* Save initial profile */
    FILE *fi = fopen("data/phase1_initial.dat", "w");
    fprintf(fi, "# r  rho_total  rho_ke  rho_grad  rho_mass  rho_pot  count\n");
    for (int b = 0; b < nbins; b++) {
        fprintf(fi, "%.4f  %.8e  %.8e  %.8e  %.8e  %.8e  %d\n",
                r_bins[b], rho_bins[b], rho_ke[b], rho_grad[b],
                rho_mass[b], rho_pot[b], counts[b]);
    }
    fclose(fi);
    printf("Initial profile saved.\n");

    /* Equilibrate with absorbing BC */
    int n_steps = (int)(T_equil / g->dt);
    int snap_every = n_steps / 10;
    printf("Equilibrating: %d steps, dt=%.5f\n", n_steps, g->dt);

    /* Track energy over time */
    FILE *fe = fopen("data/phase1_energy_vs_t.dat", "w");
    fprintf(fe, "# t  E_total  fc  winding\n");

    for (int step = 0; step <= n_steps; step++) {
        if (step > 0) {
            verlet_full_step(g, BIMODAL, -1);
            apply_damping_xy(g);
        }

        if (step % snap_every == 0) {
            Result res;
            snprintf(res.label, sizeof(res.label), "step_%d", step);
            compute_diagnostics(g, BIMODAL, &res);
            double t = step * g->dt;
            fprintf(fe, "%.4f  %.4f  %.4f  %.4f\n", t, res.energy, res.fc, res.winding);
            printf("  t=%6.1f  E=%10.1f  fc=%.4f  w=%+.3f\n",
                   t, res.energy, res.fc, res.winding);
            fflush(stdout);

            if (check_blowup(g)) {
                printf("BLOWUP at step %d!\n", step);
                break;
            }
        }
    }
    fclose(fe);

    /* Measure equilibrated profile */
    compute_radial_profile(g, BIMODAL, mass2, r_bins, rho_bins,
                           rho_ke, rho_grad, rho_mass, rho_pot, counts, nbins, dr);

    /* Save equilibrated profile */
    FILE *fp = fopen("data/phase1_equilibrated.dat", "w");
    fprintf(fp, "# r  rho_total  rho_ke  rho_grad  rho_mass  rho_pot  count\n");
    for (int b = 0; b < nbins; b++) {
        fprintf(fp, "%.4f  %.8e  %.8e  %.8e  %.8e  %.8e  %d\n",
                r_bins[b], rho_bins[b], rho_ke[b], rho_grad[b],
                rho_mass[b], rho_pot[b], counts[b]);
    }
    fclose(fp);

    /* Analysis: find vacuum reference and check for depletion */
    /* Use the outermost bins (r > 0.5*L) as vacuum reference,
       but BEFORE the damping layer (r < 0.7*L) */
    double rho_vac = 0;
    int n_vac = 0;
    for (int b = 0; b < nbins; b++) {
        if (r_bins[b] > 0.5*L && r_bins[b] < 0.65*L && counts[b] > 0) {
            rho_vac += rho_bins[b];
            n_vac++;
        }
    }
    if (n_vac > 0) rho_vac /= n_vac;

    printf("\n=== DEPLETION ANALYSIS ===\n");
    printf("Vacuum reference rho_vac = %.8e (avg of r=%.1f..%.1f)\n",
           rho_vac, 0.5*L, 0.65*L);

    /* Also compute time-averaged profile from multiple snapshots */
    /* For now, use the equilibrated snapshot */
    int found_depletion = 0;
    double min_deficit = 0;
    double r_min_deficit = 0;

    printf("\n  r       rho(r)         rho-rho_vac    delta/rho_vac\n");
    printf("  -----  -------------  -------------  -------------\n");
    for (int b = 0; b < nbins; b++) {
        if (counts[b] < 10) continue;
        double deficit = rho_bins[b] - rho_vac;
        double frac = deficit / (fabs(rho_vac) + 1e-30);
        if (r_bins[b] > 1.0 && r_bins[b] < 0.65*L) {
            printf("  %5.2f  %13.6e  %+13.6e  %+10.4f\n",
                   r_bins[b], rho_bins[b], deficit, frac);
        }
        /* Check for depletion outside core (r > 4) */
        if (r_bins[b] > 4.0 && r_bins[b] < 0.65*L && deficit < min_deficit) {
            min_deficit = deficit;
            r_min_deficit = r_bins[b];
            found_depletion = 1;
        }
    }

    printf("\n=== VERDICT ===\n");
    if (found_depletion && min_deficit < -1e-6 * fabs(rho_vac)) {
        printf("DEPLETION FOUND at r=%.2f: delta_rho/rho_vac = %.6f\n",
               r_min_deficit, min_deficit / (fabs(rho_vac) + 1e-30));
        printf("The braid creates a region where rho < rho_vac outside the core.\n");
        printf("→ Proceed to Phase 2.\n");
    } else if (rho_vac < 1e-10) {
        printf("Vacuum energy density is essentially ZERO (rho_vac=%.2e).\n", rho_vac);
        printf("With no background field, there is nothing to deplete.\n");
        printf("The braid IS the field — it adds energy, doesn't subtract.\n");
        printf("→ Phase 2 with explicit background needed to test depletion.\n");
    } else {
        printf("NO DEPLETION: rho(r) >= rho_vac everywhere outside the core.\n");
        printf("The braid ADDS energy everywhere. Minimum deficit outside core:\n");
        printf("  r=%.2f, delta_rho/rho_vac = %+.6f\n",
               r_min_deficit, min_deficit / (fabs(rho_vac) + 1e-30));
        printf("→ The braid doesn't consume background field.\n");
        printf("→ Phase 2 with explicit background still worth testing.\n");
    }

    /* Also check: does the braid have negative energy ANYWHERE?
       This would indicate field consumption */
    printf("\n=== COMPONENT ANALYSIS ===\n");
    printf("At the core (r≈0):\n");
    if (counts[0] > 0) {
        printf("  rho_ke   = %.6e\n", rho_ke[0]);
        printf("  rho_grad = %.6e\n", rho_grad[0]);
        printf("  rho_mass = %.6e\n", rho_mass[0]);
        printf("  rho_pot  = %.6e (sign: %s)\n", rho_pot[0],
               rho_pot[0] < 0 ? "NEGATIVE — potential energy consumed!" : "positive");
        printf("  rho_tot  = %.6e\n", rho_bins[0]);
    }

    printf("\nAt intermediate r (r≈5):\n");
    int b5 = (int)(5.0 / dr);
    if (b5 < nbins && counts[b5] > 0) {
        printf("  rho_ke   = %.6e\n", rho_ke[b5]);
        printf("  rho_grad = %.6e\n", rho_grad[b5]);
        printf("  rho_mass = %.6e\n", rho_mass[b5]);
        printf("  rho_pot  = %.6e\n", rho_pot[b5]);
        printf("  rho_tot  = %.6e\n", rho_bins[b5]);
    }

    /* Cleanup */
    free(r_bins); free(rho_bins);
    free(rho_ke); free(rho_grad); free(rho_mass); free(rho_pot);
    free(counts);
    grid_free(g);

    printf("\nPhase 1 complete. Data saved to data/phase1_*.dat\n");
    return 0;
}
