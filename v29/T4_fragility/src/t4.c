/*  T4: Topological Fragility Under Perturbation
 *
 *  Phase 1: Settle bimodal braid for T=200
 *  Phase 2: For each epsilon, restore settled state, apply Gaussian
 *           perturbation at core, evolve T=300 more, track winding(t).
 *
 *  Build: gcc -O3 -fopenmp -o t4 src/t4.c -lm
 */

#include "../../src/braid_core.h"

#define N_GRID   96
#define L_BOX    20.0
#define T_SETTLE 200.0
#define T_PERTURB 300.0
#define N_EPS    11
#define N_TRACK  10   /* intermediate winding samples during perturbation phase */

static const double EPS_LIST[N_EPS] = {
    0.05, 0.1, 0.2, 0.3, 0.5, 0.8, 1.0, 1.5, 2.0, 3.0, 5.0
};

/* Apply Gaussian perturbation centered at origin:
 *   phi_a += eps * exp(-r^2/(2*sigma^2)) * n_a
 * with n = (1/sqrt(3), 1/sqrt(3), 1/sqrt(3))
 */
static void apply_perturbation(Grid *g, double eps, double sigma) {
    int N = g->N, NN = N * N;
    double dx = g->dx, L = g->L;
    double inv_2s2 = 1.0 / (2.0 * sigma * sigma);
    double na = 1.0 / sqrt(3.0);

    for (int i = 0; i < N; i++) {
        double x = -L + i * dx;
        for (int j = 0; j < N; j++) {
            double y = -L + j * dx;
            double rp2 = x * x + y * y;
            /* Skip if far from core (>4 sigma in xy) */
            if (rp2 > 16.0 * sigma * sigma) continue;
            for (int k = 0; k < N; k++) {
                double z = -L + k * dx;
                double r2 = rp2 + z * z;
                double bump = eps * na * exp(-r2 * inv_2s2);
                int idx = i * NN + j * N + k;
                for (int a = 0; a < NFIELDS; a++)
                    g->phi[a][idx] += bump;
            }
        }
    }
}

/* Save/restore grid state (phi + vel only) */
static void save_state(Grid *g, double **phi_bak, double **vel_bak) {
    int N3 = g->N * g->N * g->N;
    for (int a = 0; a < NFIELDS; a++) {
        memcpy(phi_bak[a], g->phi[a], N3 * sizeof(double));
        memcpy(vel_bak[a], g->vel[a], N3 * sizeof(double));
    }
}

static void restore_state(Grid *g, double **phi_bak, double **vel_bak) {
    int N3 = g->N * g->N * g->N;
    for (int a = 0; a < NFIELDS; a++) {
        memcpy(g->phi[a], phi_bak[a], N3 * sizeof(double));
        memcpy(g->vel[a], vel_bak[a], N3 * sizeof(double));
    }
}

int main(void) {
    printf("=== T4: Topological Fragility ===\n");
    printf("N=%d, L=%.1f, T_settle=%.0f, T_perturb=%.0f\n\n",
           N_GRID, L_BOX, T_SETTLE, T_PERTURB);

    bimodal_init_params();

    /* Allocate grid */
    Grid *g = grid_alloc(N_GRID, L_BOX);
    int N3 = g->N * g->N * g->N;

    /* Allocate backup arrays */
    double *phi_bak[NFIELDS], *vel_bak[NFIELDS];
    for (int a = 0; a < NFIELDS; a++) {
        phi_bak[a] = malloc(N3 * sizeof(double));
        vel_bak[a] = malloc(N3 * sizeof(double));
    }

    /* ============================================================
       Phase 1: Settle
       ============================================================ */
    printf("Phase 1: Settling for T=%.0f ...\n", T_SETTLE);
    fflush(stdout);

    init_braid(g, BIMODAL, -1);
    compute_forces(g, BIMODAL, 0.0);  /* mass2=0 */

    int n_settle = (int)(T_SETTLE / g->dt);
    int damp_every = 5;
    double t_wall0 = omp_get_wtime();

    for (int step = 0; step < n_settle; step++) {
        verlet_full_step(g, BIMODAL, 0.0);
        if (step % damp_every == 0)
            apply_damping_xy(g);
        if (check_blowup(g)) {
            printf("  BLOWUP during settle at step %d!\n", step);
            return 1;
        }
    }

    double t_settle_wall = omp_get_wtime() - t_wall0;
    printf("  Settled in %.1f sec wall time\n", t_settle_wall);

    /* Diagnostics after settle */
    Result res_settle;
    compute_diagnostics(g, BIMODAL, &res_settle);
    snprintf(res_settle.label, sizeof(res_settle.label), "settled");
    res_settle.stable = 1;
    printf("  After settle: w=%+.3f, fc=%.4f, |P|=%.4f, E=%.1f\n",
           res_settle.winding, res_settle.fc, res_settle.peak_P, res_settle.energy);

    double w_init = res_settle.winding;
    printf("  Initial winding (reference): %.3f\n\n", w_init);
    fflush(stdout);

    /* Save settled state */
    save_state(g, phi_bak, vel_bak);

    /* ============================================================
       Phase 2: Perturbation scan
       ============================================================ */
    printf("Phase 2: Scanning %d perturbation amplitudes ...\n\n", N_EPS);
    printf("%-8s  %-8s  %-8s  %-8s  %-12s  %-6s  winding trajectory\n",
           "eps", "w_final", "fc", "|P|", "energy", "stable");
    printf("------   ------   ------   ------   ----------   ------  ------------------\n");
    fflush(stdout);

    /* Results storage for TSV */
    double eps_arr[N_EPS], w_final_arr[N_EPS], fc_arr[N_EPS];
    double peak_P_arr[N_EPS], energy_arr[N_EPS];
    int stable_arr[N_EPS];
    double w_track[N_EPS][N_TRACK]; /* intermediate winding values */

    int n_perturb = (int)(T_PERTURB / g->dt);
    int track_interval = n_perturb / N_TRACK;

    double sigma_pert = 3.0;  /* = R_tube */

    for (int ie = 0; ie < N_EPS; ie++) {
        double eps = EPS_LIST[ie];
        double t_w0 = omp_get_wtime();

        /* Restore settled state */
        restore_state(g, phi_bak, vel_bak);

        /* Apply perturbation */
        apply_perturbation(g, eps, sigma_pert);

        /* Recompute forces after perturbation */
        compute_forces(g, BIMODAL, 0.0);

        /* Evolve */
        int blown = 0;
        int track_idx = 0;

        for (int step = 0; step < n_perturb; step++) {
            verlet_full_step(g, BIMODAL, 0.0);
            if (step % damp_every == 0)
                apply_damping_xy(g);

            /* Track winding at intermediate times */
            if (track_interval > 0 && step > 0 && (step % track_interval == 0) && track_idx < N_TRACK) {
                w_track[ie][track_idx] = compute_winding(g);
                track_idx++;
            }

            if (check_blowup(g)) {
                blown = 1;
                /* Fill remaining track slots */
                for (; track_idx < N_TRACK; track_idx++)
                    w_track[ie][track_idx] = 999.0;
                break;
            }
        }

        /* Fill any remaining track slots */
        for (; track_idx < N_TRACK; track_idx++)
            w_track[ie][track_idx] = compute_winding(g);

        /* Final diagnostics */
        Result res;
        compute_diagnostics(g, BIMODAL, &res);
        res.stable = !blown;

        eps_arr[ie] = eps;
        w_final_arr[ie] = res.winding;
        fc_arr[ie] = res.fc;
        peak_P_arr[ie] = res.peak_P;
        energy_arr[ie] = res.energy;
        stable_arr[ie] = res.stable;

        double dt_wall = omp_get_wtime() - t_w0;

        /* Print row with winding trajectory */
        printf("%-8.3f  %+7.3f  %7.4f  %7.4f  %11.1f  %-6s  [",
               eps, res.winding, res.fc, res.peak_P, res.energy,
               res.stable ? "OK" : "BLOW");
        for (int t = 0; t < N_TRACK; t++) {
            if (t > 0) printf(",");
            if (w_track[ie][t] > 900.0)
                printf("X");
            else
                printf("%+.2f", w_track[ie][t]);
        }
        printf("]  (%.0fs)\n", dt_wall);
        fflush(stdout);
    }

    /* ============================================================
       Find epsilon_crit
       ============================================================ */
    printf("\n--- Critical Perturbation Analysis ---\n");
    double eps_crit = -1.0;
    for (int ie = 0; ie < N_EPS; ie++) {
        if (!stable_arr[ie] || fabs(w_final_arr[ie] - w_init) > 0.3) {
            eps_crit = eps_arr[ie];
            printf("eps_crit = %.3f (first eps where |dw| > 0.3 or blowup)\n", eps_crit);
            printf("  w_init = %.3f, w_final = %.3f, |dw| = %.3f\n",
                   w_init, w_final_arr[ie], fabs(w_final_arr[ie] - w_init));
            break;
        }
    }
    if (eps_crit < 0) {
        printf("No topological failure detected up to eps=%.1f\n", EPS_LIST[N_EPS - 1]);
        eps_crit = EPS_LIST[N_EPS - 1]; /* lower bound */
    }

    /* A0 is the field amplitude */
    double A0 = BIMODAL[0];
    printf("eps_crit / A0 = %.2f / %.2f = %.2f\n", eps_crit, A0, eps_crit / A0);
    if (eps_crit / A0 > 5.0)
        printf("Assessment: EFFECTIVELY TOPOLOGICAL (ratio > 5)\n");
    else if (eps_crit / A0 > 1.0)
        printf("Assessment: MODERATELY ROBUST (ratio 1-5)\n");
    else
        printf("Assessment: FRAGILE (ratio < 1)\n");

    /* ============================================================
       Write TSV
       ============================================================ */
    FILE *fp = fopen("data/t4_results.tsv", "w");
    if (fp) {
        fprintf(fp, "eps\tw_init\tw_final\tdw\tfc\tpeak_P\tenergy\tstable");
        for (int t = 0; t < N_TRACK; t++)
            fprintf(fp, "\tw_%d", t + 1);
        fprintf(fp, "\n");
        for (int ie = 0; ie < N_EPS; ie++) {
            fprintf(fp, "%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.1f\t%d",
                    eps_arr[ie], w_init, w_final_arr[ie],
                    w_final_arr[ie] - w_init,
                    fc_arr[ie], peak_P_arr[ie], energy_arr[ie], stable_arr[ie]);
            for (int t = 0; t < N_TRACK; t++)
                fprintf(fp, "\t%.4f", w_track[ie][t]);
            fprintf(fp, "\n");
        }
        fclose(fp);
        printf("\nResults written to data/t4_results.tsv\n");
    } else {
        printf("\nWARNING: Could not open data/t4_results.tsv for writing\n");
    }

    /* Cleanup */
    for (int a = 0; a < NFIELDS; a++) {
        free(phi_bak[a]);
        free(vel_bak[a]);
    }
    grid_free(g);

    printf("\nDone.\n");
    return 0;
}
