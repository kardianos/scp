/*  t10d.c — Kibble-Zurek Defect Formation
 *
 *  Start with random high-energy field, cool via damping.
 *  Do braided structures form spontaneously?
 *
 *  If YES: braids are natural topological defects of this field theory
 *  If NO: braids are artifacts of our specific initialization
 *
 *  Method:
 *  1. Initialize random field with controlled amplitude and spectrum
 *  2. Evolve with viscous damping (γ·v term) that slowly extracts energy
 *  3. Monitor for localized structures, winding, fc
 *
 *  Build: gcc -O3 -fopenmp -o t10d src/t10d.c -lm
 */

#include "../../src/braid_core.h"
#include <float.h>

/* Apply uniform viscous damping */
static void apply_viscous_damping(Grid *g, double gamma_damp) {
    int N3 = g->N * g->N * g->N;
    double factor = 1.0 - gamma_damp * g->dt;
    if (factor < 0) factor = 0;
    for (int a = 0; a < NFIELDS; a++)
        for (int idx = 0; idx < N3; idx++)
            g->vel[a][idx] *= factor;
}

/* Initialize random field with specified amplitude and maximum k */
static void init_random_field(Grid *g, double amplitude, int k_max, unsigned int seed) {
    int N = g->N, NN = N*N;
    double L = g->L, dx = g->dx;
    srand(seed);

    grid_zero(g);

    /* Direct-space initialization: random amplitudes with spatial smoothing */
    /* First set random values */
    for (int a = 0; a < NFIELDS; a++) {
        for (int i = 0; i < N; i++)
            for (int j = 0; j < N; j++)
                for (int k = 0; k < N; k++) {
                    int idx = i*NN + j*N + k;
                    g->phi[a][idx] = amplitude * (2.0 * rand() / RAND_MAX - 1.0);
                    g->vel[a][idx] = amplitude * (2.0 * rand() / RAND_MAX - 1.0);
                }

        /* Smooth: apply diffusion (Gaussian blur) to cut off high-k modes */
        /* sigma_smooth = L/(pi*k_max) → smooths out modes with k > k_max */
        int n_smooth = N / (2 * k_max + 1);
        if (n_smooth < 1) n_smooth = 1;

        double *tmp = calloc(N*N*N, sizeof(double));
        for (int iter = 0; iter < n_smooth; iter++) {
            /* Jacobi smoothing step */
            for (int i = 1; i < N-1; i++)
                for (int j = 1; j < N-1; j++)
                    for (int k = 0; k < N; k++) {
                        int idx = i*NN + j*N + k;
                        int kp = (k+1)%N, km = (k-1+N)%N;
                        tmp[idx] = (g->phi[a][idx+NN] + g->phi[a][idx-NN]
                                  + g->phi[a][idx+N] + g->phi[a][idx-N]
                                  + g->phi[a][i*NN+j*N+kp] + g->phi[a][i*NN+j*N+km]
                                  + 6.0*g->phi[a][idx]) / 12.0;
                    }
            memcpy(g->phi[a], tmp, N*N*N*sizeof(double));
        }
        free(tmp);
    }
}

/* Measure local energy density and find peaks */
static int find_structures(Grid *g, const double *phys, double mass2,
                           double *max_rho, double *avg_rho,
                           int *n_peaks, double peak_threshold) {
    int N = g->N, NN = N*N;
    double dx = g->dx, L = g->L;
    double rho_sum = 0, rho_max = 0;
    int count = 0;
    int peaks = 0;

    for (int i = 2; i < N-2; i++) {
        for (int j = 2; j < N-2; j++) {
            for (int k = 0; k < N; k++) {
                int idx = i*NN + j*N + k;
                double rho = 0;
                for (int a = 0; a < NFIELDS; a++)
                    rho += g->phi[a][idx] * g->phi[a][idx];
                rho_sum += rho;
                count++;
                if (rho > rho_max) rho_max = rho;
            }
        }
    }

    *avg_rho = rho_sum / count;
    *max_rho = rho_max;

    /* Count peaks: points where rho > threshold * avg_rho and is local max */
    double thresh = peak_threshold * (*avg_rho);
    for (int i = 3; i < N-3; i++) {
        for (int j = 3; j < N-3; j++) {
            for (int k = 0; k < N; k++) {
                int idx = i*NN + j*N + k;
                double rho = 0;
                for (int a = 0; a < NFIELDS; a++)
                    rho += g->phi[a][idx] * g->phi[a][idx];
                if (rho < thresh) continue;

                /* Check if local max in all 6 neighbors */
                int kp = (k+1)%N, km = (k-1+N)%N;
                int is_max = 1;
                int nb[6] = {idx+NN, idx-NN, idx+N, idx-N,
                             i*NN+j*N+kp, i*NN+j*N+km};
                for (int n = 0; n < 6 && is_max; n++) {
                    double rho_n = 0;
                    for (int a = 0; a < NFIELDS; a++)
                        rho_n += g->phi[a][nb[n]] * g->phi[a][nb[n]];
                    if (rho_n >= rho) is_max = 0;
                }
                if (is_max) peaks++;
            }
        }
    }
    *n_peaks = peaks;
    return peaks > 0;
}

int main(int argc, char **argv) {
    int N = 96;
    double L = 20.0;

    printf("=== T10D: Kibble-Zurek Defect Formation ===\n\n");
    printf("Grid: N=%d, L=%.1f\n\n", N, L);

    bimodal_init_params();
    double mass2 = 2.25;

    FILE *fout = fopen("data/t10d_cooling.tsv", "w");
    fprintf(fout, "seed\tgamma\tt\tenergy\tmax_rho\tavg_rho\tn_peaks\tfc\twinding\n");

    /* Test different cooling rates and initial amplitudes */
    double gammas[] = {0.01, 0.05, 0.1, 0.5};
    int ngamma = 4;
    double init_amps[] = {0.5, 1.0, 2.0};
    int namps = 3;

    printf("=== Cooling sweep ===\n");
    printf("gamma\tamp\tseed\tT_final\tE_final\tfc\tmax_rho\tn_peaks\twinding\n");

    for (int ig = 0; ig < ngamma; ig++) {
        double gamma = gammas[ig];
        for (int ia = 0; ia < namps; ia++) {
            double amp0 = init_amps[ia];
            for (int seed = 42; seed <= 44; seed++) {  /* 3 seeds per config */
                Grid *g = grid_alloc(N, L);
                init_random_field(g, amp0, 8, seed);
                compute_forces(g, BIMODAL, mass2);

                double T_cool = 500.0;
                int cool_steps = (int)(T_cool / g->dt);
                int meas_every = cool_steps / 50;
                if (meas_every < 1) meas_every = 1;

                double final_fc = 0, final_E = 0, final_wind = 0;
                double final_max_rho = 0, final_avg_rho = 0;
                int final_peaks = 0;

                for (int s = 0; s < cool_steps; s++) {
                    verlet_full_step(g, BIMODAL, mass2);
                    apply_viscous_damping(g, gamma);

                    if (s % meas_every == 0) {
                        double max_rho, avg_rho;
                        int n_peaks;
                        find_structures(g, BIMODAL, mass2, &max_rho, &avg_rho,
                                       &n_peaks, 5.0);

                        Result res;
                        compute_diagnostics(g, BIMODAL, &res);

                        fprintf(fout, "%d\t%.2f\t%.1f\t%.1f\t%.4f\t%.4f\t%d\t%.4f\t%.3f\n",
                                seed, gamma, s*g->dt, res.energy, max_rho, avg_rho,
                                n_peaks, res.fc, res.winding);

                        final_fc = res.fc;
                        final_E = res.energy;
                        final_wind = res.winding;
                        final_max_rho = max_rho;
                        final_avg_rho = avg_rho;
                        final_peaks = n_peaks;

                        /* Check for blowup */
                        if (check_blowup(g)) {
                            printf("%.2f\t%.1f\t%d\tBLOWUP at t=%.0f\n",
                                   gamma, amp0, seed, s*g->dt);
                            break;
                        }
                    }
                }

                printf("%.2f\t%.1f\t%d\t%.0f\t%.1f\t%.4f\t%.4f\t%d\t%.3f\n",
                       gamma, amp0, seed, T_cool, final_E, final_fc,
                       final_max_rho, final_peaks, final_wind);

                grid_free(g);
            }
        }
    }

    fclose(fout);

    /* ============================================================
       Part 2: Check if any cooled configuration has braid-like structure
       by measuring fc, winding, and energy localization
       ============================================================ */
    printf("\n=== Part 2: Detailed analysis of best cooling ===\n");
    {
        /* Use moderate cooling (gamma=0.05) with largest amplitude */
        Grid *g = grid_alloc(N, L);
        init_random_field(g, 2.0, 8, 42);
        compute_forces(g, BIMODAL, mass2);

        double T_cool = 1000.0;
        int cool_steps = (int)(T_cool / g->dt);

        printf("  Cooling gamma=0.05, amp=2.0, T=%.0f...\n", T_cool);
        fflush(stdout);

        for (int s = 0; s < cool_steps; s++) {
            verlet_full_step(g, BIMODAL, mass2);
            apply_viscous_damping(g, 0.05);

            if (s == cool_steps/4 || s == cool_steps/2 ||
                s == 3*cool_steps/4 || s == cool_steps-1) {
                Result res;
                compute_diagnostics(g, BIMODAL, &res);
                printf("  t=%6.0f: E=%8.1f fc=%.4f |P|=%.4f w=%.3f\n",
                       s*g->dt, res.energy, res.fc, res.peak_P, res.winding);
            }
        }

        /* Now remove damping and check stability */
        printf("\n  Removing damping, evolving for T=200...\n");
        double T_free = 200.0;
        int free_steps = (int)(T_free / g->dt);
        for (int s = 0; s < free_steps; s++) {
            verlet_full_step(g, BIMODAL, mass2);
            apply_damping_xy(g); /* still absorb edges */
        }

        Result res_free;
        compute_diagnostics(g, BIMODAL, &res_free);
        printf("  After free evolution: E=%.1f fc=%.4f |P|=%.4f w=%.3f\n",
               res_free.energy, res_free.fc, res_free.peak_P, res_free.winding);

        if (res_free.fc > 0.5)
            printf("  → LOCALIZED STRUCTURE PERSISTS\n");
        else
            printf("  → Field dispersed (no stable structure formed)\n");

        grid_free(g);
    }

    printf("\n=== T10D Complete ===\n");
    printf("Data in data/t10d_cooling.tsv\n");
    return 0;
}
