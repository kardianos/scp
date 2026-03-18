/*  T2: Initialization mass independence test
 *
 *  Scans m_init in {0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0}
 *  with dynamics mass = 0 (massless EOM) throughout.
 *
 *  Also runs two extras:
 *    - m_init=0, dynamics mass=0   (pure massless)
 *    - m_init=0, dynamics mass²=2.25 (mass in dynamics only)
 *
 *  Build: gcc -O3 -fopenmp -o t2 src/t2.c -lm
 */

#include "../../src/braid_core.h"

#define N_GRID  128
#define L_BOX   20.0
#define T_FINAL 500.0

typedef struct {
    double m_init;
    double mass2_dyn;
    const char *tag;
} RunConfig;

int main(void) {
    bimodal_init_params();

    printf("=== T2: Initialization Mass Independence ===\n");
    printf("N=%d  L=%.1f  T=%.1f  bimodal t=0.85\n", N_GRID, L_BOX, T_FINAL);
    printf("Params:");
    for (int d = 0; d < NDIM; d++) printf(" %s=%.3f", PNAME[d], BIMODAL[d]);
    printf("\n\n");

    /* Define runs: 8 m_init scan + 2 extras */
    RunConfig runs[] = {
        {0.00, 0.0, "m_init=0.00 dyn=0"},
        {0.25, 0.0, "m_init=0.25 dyn=0"},
        {0.50, 0.0, "m_init=0.50 dyn=0"},
        {0.75, 0.0, "m_init=0.75 dyn=0"},
        {1.00, 0.0, "m_init=1.00 dyn=0"},
        {1.25, 0.0, "m_init=1.25 dyn=0"},
        {1.50, 0.0, "m_init=1.50 dyn=0"},
        {2.00, 0.0, "m_init=2.00 dyn=0"},
        /* extras */
        {0.00, 0.0,  "EXTRA: pure massless"},
        {0.00, 2.25, "EXTRA: m_init=0 dyn=2.25"},
    };
    int nruns = sizeof(runs) / sizeof(runs[0]);

    Result results[16];

    for (int r = 0; r < nruns; r++) {
        printf("--- Run %d/%d: %s ---\n", r+1, nruns, runs[r].tag);
        fflush(stdout);

        Grid *g = grid_alloc(N_GRID, L_BOX);
        init_braid(g, BIMODAL, runs[r].m_init);
        compute_forces(g, BIMODAL, runs[r].mass2_dyn);

        int nsteps = (int)(T_FINAL / g->dt + 0.5);
        int print_every = nsteps / 10;
        if (print_every < 1) print_every = 1;
        int blown = 0;

        double t0 = omp_get_wtime();
        for (int s = 0; s < nsteps; s++) {
            verlet_full_step(g, BIMODAL, runs[r].mass2_dyn);
            apply_damping_xy(g);
            if (s % print_every == 0) {
                double pct = 100.0 * s / nsteps;
                double elapsed = omp_get_wtime() - t0;
                printf("  step %d/%d (%.0f%%)  elapsed %.1fs\n", s, nsteps, pct, elapsed);
                fflush(stdout);
            }
            if (check_blowup(g)) { blown = 1; break; }
        }
        double wall = omp_get_wtime() - t0;

        Result *res = &results[r];
        snprintf(res->label, sizeof(res->label), "%s", runs[r].tag);
        res->stable = !blown;

        if (!blown) {
            compute_diagnostics(g, BIMODAL, res);
        } else {
            res->energy = 0;
            res->fc = 0;
            res->l2_frac = 0;
            res->transverse_l2 = 0;
            res->torsion_flux = 0;
            res->peak_P = 0;
            res->winding = 0;
        }
        print_result(res);
        printf("  wall time: %.1fs\n\n", wall);
        fflush(stdout);

        grid_free(g);
    }

    /* Summary table */
    printf("\n========== SUMMARY ==========\n");
    printf("%-28s %10s %7s %7s %7s %7s %7s %7s %s\n",
           "Config", "Energy", "fc", "l2", "trans", "torsion", "|P|", "wind", "status");
    for (int r = 0; r < nruns; r++) {
        Result *res = &results[r];
        printf("%-28s %10.1f %7.4f %7.4f %7.4f %7.4f %7.4f %+7.3f %s\n",
               res->label, res->energy, res->fc, res->l2_frac,
               res->transverse_l2, res->torsion_flux, res->peak_P,
               res->winding, res->stable ? "OK" : "BLOWUP");
    }

    /* Write TSV */
    FILE *fp = fopen("data/t2_results.tsv", "w");
    if (fp) {
        fprintf(fp, "config\tm_init\tmass2_dyn\tenergy\tfc\tl2\ttrans_l2\ttorsion\tpeak_P\twinding\tstable\n");
        for (int r = 0; r < nruns; r++) {
            Result *res = &results[r];
            fprintf(fp, "%s\t%.2f\t%.2f\t%.1f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.4f\t%d\n",
                    res->label, runs[r].m_init, runs[r].mass2_dyn,
                    res->energy, res->fc, res->l2_frac, res->transverse_l2,
                    res->torsion_flux, res->peak_P, res->winding, res->stable);
        }
        fclose(fp);
        printf("\nResults written to data/t2_results.tsv\n");
    } else {
        printf("\nWARNING: could not open data/t2_results.tsv for writing\n");
    }

    return 0;
}
