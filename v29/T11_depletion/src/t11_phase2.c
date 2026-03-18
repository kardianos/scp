/*  t11_phase2.c — Braid in explicit background field: measure depletion
 *
 *  Phase 1 showed: with no background, the braid IS the field. rho→0 at large r.
 *  There is nothing to deplete.
 *
 *  Phase 2: Add a UNIFORM BACKGROUND field and see if the braid depletes it.
 *
 *  Setup:
 *  - Initialize braid as usual
 *  - ADD a spatially uniform standing wave in each field:
 *    phi_a → phi_a_braid + A_bg * cos(k_bg * z)  (static initial condition)
 *    vel_a → vel_a_braid + A_bg * omega_bg * sin(k_bg * z)
 *    where k_bg = 2*pi/Lz, omega_bg = sqrt(k_bg² + m²)
 *  - Run dynamics with absorbing BC in xy, periodic in z
 *  - Measure radial energy density profile at multiple times
 *  - Compare: does rho decrease NEAR the braid relative to FAR from it?
 *
 *  The background k_bg mode is independent of x,y → uniform rho_bg.
 *  If the braid consumes this background, we'll see rho_bg decrease near the core.
 *
 *  Test multiple A_bg values: 0.01, 0.05, 0.1, 0.3
 *
 *  Build: gcc -O3 -fopenmp -o t11p2 src/t11_phase2.c -lm
 */

#include "../../src/braid_core.h"
#include <float.h>

/* ================================================================
   Radial profile (same as Phase 1)
   ================================================================ */
static void compute_radial_profile(Grid *g, const double *phys, double mass2,
                                    double *r_bins, double *rho_bins,
                                    int *counts, int nbins, double dr) {
    int N = g->N, NN = N*N;
    double dx = g->dx, L = g->L;
    double mu = phys[12], kappa = phys[13];

    for (int b = 0; b < nbins; b++) {
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
            if (b >= nbins) continue;

            for (int k = 0; k < N; k++) {
                int idx = i*NN + j*N + k;
                int kp = (k+1)%N, km = (k-1+N)%N;

                double ek = 0, eg = 0, em_val = 0;
                for (int a = 0; a < NFIELDS; a++) {
                    ek += 0.5 * g->vel[a][idx] * g->vel[a][idx];
                    double gx = (g->phi[a][idx+NN] - g->phi[a][idx-NN]) / (2*dx);
                    double gy = (g->phi[a][idx+N]  - g->phi[a][idx-N])  / (2*dx);
                    double gz = (g->phi[a][i*NN+j*N+kp] - g->phi[a][i*NN+j*N+km]) / (2*dx);
                    eg += 0.5 * (gx*gx + gy*gy + gz*gz);
                    em_val += 0.5 * mass2 * g->phi[a][idx] * g->phi[a][idx];
                }

                double p0 = g->phi[0][idx], p1 = g->phi[1][idx], p2 = g->phi[2][idx];
                double P = p0 * p1 * p2;
                double ep = (mu/2.0) * P*P / (1.0 + kappa*P*P);

                rho_bins[b] += ek + eg + em_val + ep;
                counts[b]++;
            }
        }
    }

    for (int b = 0; b < nbins; b++) {
        if (counts[b] > 0) rho_bins[b] /= counts[b];
    }
}

/* ================================================================
   Run one A_bg test
   ================================================================ */
static void run_test(double A_bg, const double *phys, double mass2) {
    int N = 128;
    double L = 20.0;
    double T_equil = 200.0;  /* pre-equilibrate braid */
    double T_measure = 200.0; /* measurement period with background */

    printf("\n===== A_bg = %.3f =====\n", A_bg);

    Grid *g = grid_alloc(N, L);
    init_braid(g, phys, -1);
    compute_forces(g, phys, -1);

    /* Pre-equilibrate braid without background */
    printf("Pre-equilibrating braid for T=%.0f...\n", T_equil);
    int n_equil = (int)(T_equil / g->dt);
    for (int step = 0; step < n_equil; step++) {
        verlet_full_step(g, phys, -1);
        apply_damping_xy(g);
    }
    printf("Pre-equilibration done.\n");

    /* Now add background standing wave:
       phi_a += A_bg * cos(k_bg * z)
       vel_a += A_bg * omega_bg * sin(k_bg * z)
       Use k_bg = 2*pi/(2*L) = pi/L (fundamental mode of periodic box) */
    double k_bg = PI / L;
    double omega_bg = sqrt(k_bg * k_bg + mass2);
    double rho_bg_theory = 0.5 * NFIELDS * A_bg * A_bg * (omega_bg*omega_bg + k_bg*k_bg + mass2);
    /* For standing wave: <rho> = (omega² + k² + m²)/2 * A² per field
       = omega² * A² per field (using dispersion) ... actually:
       phi = A cos(kz) cos(wt) → <v²>=A²w²/2, <(dz phi)²>=A²k²/2, <m²phi²>=A²m²/2
       <rho> = A²(w²+k²+m²)/4 per field ... let me be careful.
       phi = A cos(kz), vel = 0 at t=0. This evolves as A cos(kz) cos(wt).
       <v²> = A² w² sin²(wt) cos²(kz), time avg = A²w²/4
       <(∂z phi)²> = A²k² sin²(kz) cos²(wt), time avg = A²k²/4
       <m² phi²> = A²m² cos²(kz) cos²(wt), time avg = A²m²/4
       So <rho> = A²(w²+k²+m²)/4 = A²(2w²)/4 = A²w²/2 per field.
       NO — at t=0 we add both phi and vel:
       phi += A cos(kz), vel += A w sin(kz)
       This is a traveling wave: A cos(kz - wt).
       <v²> = A²w² sin²(kz-wt), time avg = A²w²/2
       <(∂z)²> = A²k² sin²(kz-wt), time avg = A²k²/2
       <m² phi²> = A²m² cos²(kz-wt), time avg = A²m²/2
       <rho> = A²(w²+k²+m²)/2 = A²·2w²/2 = A²w² per field.
       Actually w² = k²+m² so w²+k²+m² = 2w². So <rho> = A²w² per field.
       Total: 3 * A² * w². */
    rho_bg_theory = NFIELDS * A_bg * A_bg * omega_bg * omega_bg;

    printf("Background: k=%.4f, omega=%.4f, A_bg=%.3f\n", k_bg, omega_bg, A_bg);
    printf("Expected background rho_bg = %.6e (time-averaged)\n", rho_bg_theory);

    int NN = N*N, N3 = N*N*N;
    double dx = g->dx;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            for (int k = 0; k < N; k++) {
                int idx = i*NN + j*N + k;
                double z = -L + k * dx;
                for (int a = 0; a < NFIELDS; a++) {
                    /* Traveling wave in +z direction */
                    g->phi[a][idx] += A_bg * cos(k_bg * z);
                    g->vel[a][idx] += A_bg * omega_bg * sin(k_bg * z);
                }
            }
        }
    }

    /* Recompute forces after adding background */
    compute_forces(g, phys, -1);

    /* Measure profile at multiple times */
    int nbins = 50;
    double dr = L / nbins;
    double *r_bins   = calloc(nbins, sizeof(double));
    double *rho_bins = calloc(nbins, sizeof(double));
    int    *counts   = calloc(nbins, sizeof(int));

    /* Accumulator for time-averaged profile */
    double *rho_avg  = calloc(nbins, sizeof(double));
    int n_avg = 0;

    int n_measure = (int)(T_measure / g->dt);
    int snap_every = n_measure / 20;
    int avg_start = n_measure / 4;  /* start averaging after initial transients */

    char fname[256];
    snprintf(fname, sizeof(fname), "data/phase2_Abg%.3f_timeseries.dat", A_bg);
    FILE *ft = fopen(fname, "w");
    fprintf(ft, "# t  rho_core(r<3)  rho_mid(r=6-8)  rho_far(r=10-13)  E_total  fc\n");

    for (int step = 0; step <= n_measure; step++) {
        if (step > 0) {
            verlet_full_step(g, phys, -1);
            apply_damping_xy(g);
        }

        if (step % snap_every == 0) {
            compute_radial_profile(g, phys, mass2, r_bins, rho_bins, counts, nbins, dr);

            /* Compute regional averages */
            double rho_core = 0, rho_mid = 0, rho_far = 0;
            int nc = 0, nm = 0, nf = 0;
            for (int b = 0; b < nbins; b++) {
                if (counts[b] == 0) continue;
                if (r_bins[b] < 3.0) { rho_core += rho_bins[b]; nc++; }
                if (r_bins[b] > 6.0 && r_bins[b] < 8.0) { rho_mid += rho_bins[b]; nm++; }
                if (r_bins[b] > 10.0 && r_bins[b] < 13.0) { rho_far += rho_bins[b]; nf++; }
            }
            if (nc > 0) rho_core /= nc;
            if (nm > 0) rho_mid /= nm;
            if (nf > 0) rho_far /= nf;

            double t = (T_equil + step * g->dt);
            Result res;
            snprintf(res.label, sizeof(res.label), "s%d", step);
            compute_diagnostics(g, phys, &res);

            fprintf(ft, "%.2f  %.8e  %.8e  %.8e  %.2f  %.4f\n",
                    t, rho_core, rho_mid, rho_far, res.energy, res.fc);

            if (step % (snap_every * 5) == 0) {
                printf("  t=%6.1f  rho_core=%.4e  rho_mid=%.4e  rho_far=%.4e  fc=%.4f\n",
                       t, rho_core, rho_mid, rho_far, res.fc);
            }

            /* Accumulate for time-average */
            if (step >= avg_start) {
                for (int b = 0; b < nbins; b++)
                    rho_avg[b] += rho_bins[b];
                n_avg++;
            }

            if (check_blowup(g)) {
                printf("BLOWUP!\n");
                break;
            }
        }
    }
    fclose(ft);

    /* Time-averaged profile */
    if (n_avg > 0) {
        for (int b = 0; b < nbins; b++)
            rho_avg[b] /= n_avg;
    }

    /* Save time-averaged profile */
    snprintf(fname, sizeof(fname), "data/phase2_Abg%.3f_profile.dat", A_bg);
    FILE *fp = fopen(fname, "w");
    fprintf(fp, "# r  rho_avg  count  rho_bg_theory=%.6e\n", rho_bg_theory);
    for (int b = 0; b < nbins; b++) {
        fprintf(fp, "%.4f  %.8e  %d\n", r_bins[b], rho_avg[b], counts[b]);
    }
    fclose(fp);

    /* Analysis */
    double rho_far_avg = 0;
    int nf = 0;
    for (int b = 0; b < nbins; b++) {
        if (r_bins[b] > 10.0 && r_bins[b] < 13.0 && counts[b] > 0) {
            rho_far_avg += rho_avg[b];
            nf++;
        }
    }
    if (nf > 0) rho_far_avg /= nf;

    printf("\n  RESULT for A_bg=%.3f:\n", A_bg);
    printf("  Background rho_bg_theory = %.6e\n", rho_bg_theory);
    printf("  Far-field rho_far_avg    = %.6e\n", rho_far_avg);
    printf("  Radial profile (time-averaged):\n");
    printf("  r       rho_avg        delta/rho_far\n");
    for (int b = 0; b < nbins; b++) {
        if (counts[b] < 10) continue;
        if (r_bins[b] > 0.65 * L) continue;  /* skip damping layer */
        double frac = (rho_avg[b] - rho_far_avg) / (fabs(rho_far_avg) + 1e-30);
        if (b % 2 == 0 || r_bins[b] < 5.0) {
            printf("  %5.2f  %.6e  %+.4f\n", r_bins[b], rho_avg[b], frac);
        }
    }

    /* Check: did the background survive? */
    printf("  Background survival: rho_far/rho_theory = %.4f\n",
           rho_far_avg / (rho_bg_theory + 1e-30));

    /* Check for depletion between core and far field */
    double min_rho = 1e30;
    double r_min = 0;
    for (int b = 0; b < nbins; b++) {
        if (r_bins[b] > 4.0 && r_bins[b] < 0.6*L && counts[b] > 0) {
            if (rho_avg[b] < min_rho) {
                min_rho = rho_avg[b];
                r_min = r_bins[b];
            }
        }
    }
    if (min_rho < rho_far_avg * 0.99) {
        printf("  >>> DEPLETION: min rho = %.6e at r=%.2f (%.2f%% below far-field)\n",
               min_rho, r_min, 100.0*(1.0 - min_rho/rho_far_avg));
    } else {
        printf("  >>> NO DEPLETION: rho >= rho_far everywhere outside core\n");
    }

    free(r_bins); free(rho_bins); free(counts); free(rho_avg);
    grid_free(g);
}

/* Also run a CONTROL: same background, NO braid */
static void run_control(double A_bg, const double *phys, double mass2) {
    int N = 128;
    double L = 20.0;
    double T_measure = 200.0;

    printf("\n===== CONTROL (no braid), A_bg = %.3f =====\n", A_bg);

    Grid *g = grid_alloc(N, L);
    grid_zero(g);

    double k_bg = PI / L;
    double omega_bg = sqrt(k_bg * k_bg + mass2);

    int NN = N*N;
    double dx = g->dx;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            for (int k = 0; k < N; k++) {
                int idx = i*NN + j*N + k;
                double z = -L + k * dx;
                for (int a = 0; a < NFIELDS; a++) {
                    g->phi[a][idx] = A_bg * cos(k_bg * z);
                    g->vel[a][idx] = A_bg * omega_bg * sin(k_bg * z);
                }
            }
        }
    }

    compute_forces(g, phys, -1);

    int nbins = 50;
    double dr = L / nbins;
    double *r_bins   = calloc(nbins, sizeof(double));
    double *rho_bins = calloc(nbins, sizeof(double));
    int    *counts   = calloc(nbins, sizeof(int));

    /* Run and measure at end */
    int n_measure = (int)(T_measure / g->dt);
    for (int step = 0; step < n_measure; step++) {
        verlet_full_step(g, phys, -1);
        apply_damping_xy(g);
    }

    compute_radial_profile(g, phys, mass2, r_bins, rho_bins, counts, nbins, dr);

    printf("  Control profile (after T=%.0f):\n", T_measure);
    printf("  r       rho\n");
    double rho_center = 0, rho_far = 0;
    int nc0 = 0, nf0 = 0;
    for (int b = 0; b < nbins; b++) {
        if (counts[b] < 10) continue;
        if (r_bins[b] < 3.0) { rho_center += rho_bins[b]; nc0++; }
        if (r_bins[b] > 10.0 && r_bins[b] < 13.0) { rho_far += rho_bins[b]; nf0++; }
        if (b % 5 == 0) printf("  %5.2f  %.6e\n", r_bins[b], rho_bins[b]);
    }
    if (nc0 > 0) rho_center /= nc0;
    if (nf0 > 0) rho_far /= nf0;
    printf("  Control: rho_center/rho_far = %.6f (should be ~1 if uniform)\n",
           rho_center / (rho_far + 1e-30));

    free(r_bins); free(rho_bins); free(counts);
    grid_free(g);
}

int main(int argc, char **argv) {
    bimodal_init_params();
    double mass2 = BIMODAL[14] * BIMODAL[14];

    printf("=== T11 Phase 2: Background Field Depletion ===\n");
    printf("mass²=%.4f\n", mass2);

    /* Control run first (no braid, just background) */
    run_control(0.1, BIMODAL, mass2);

    /* Braid + background at different amplitudes */
    double A_values[] = {0.01, 0.05, 0.1, 0.3};
    int nA = 4;
    for (int i = 0; i < nA; i++) {
        run_test(A_values[i], BIMODAL, mass2);
    }

    printf("\n=== Phase 2 Complete ===\n");
    return 0;
}
