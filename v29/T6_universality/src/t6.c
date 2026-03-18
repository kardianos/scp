/*  T6: Universality — Bimodal Synergy Across Parameter Space
 *
 *  For each (mu, kappa) pair, run three configs:
 *    A) Pure A geometry with overridden (mu, kappa)
 *    B) Pure B geometry with overridden (mu, kappa)
 *    C) Bimodal (t=0.85) geometry with overridden (mu, kappa)
 *  All use mass=1.50 (m^2=2.25 via mass2_override=-1).
 *
 *  Synergy = min(bimodal_trans / max(A_trans, eps),
 *                bimodal_tor  / max(B_tor,  eps))
 *
 *  Build: cd v29/T6_universality && gcc -O3 -fopenmp -o t6 src/t6.c -lm
 */

#include "../../src/braid_core.h"

#define N_GRID  80
#define L_BOX   20.0
#define T_EVOL  200.0
#define DAMP_EVERY 5

#define N_PAIRS 10

static const double MU_LIST[N_PAIRS] = {
    -20.0, -30.0, -41.3, -60.0, -80.0,
    -100.0, -41.3, -41.3, -20.0, -80.0
};
static const double KAPPA_LIST[N_PAIRS] = {
    20.0, 30.0, 50.0, 60.0, 80.0,
    100.0, 20.0, 100.0, 50.0, 50.0
};

typedef struct {
    double mu, kappa;
    double A_trans, A_tor;
    double B_trans, B_tor;
    double bm_trans, bm_tor;
    double synergy;
    int A_ok, B_ok, bm_ok;
} PairResult;

/* Run one configuration: init, evolve T_EVOL, return diagnostics */
static void run_one(const double *params, const char *label,
                    Result *res, Grid *g)
{
    init_braid(g, params, -1);           /* m_init from params[14]=1.50 */
    compute_forces(g, params, -1);       /* m^2 from params[14]^2=2.25 */

    int n_steps = (int)(T_EVOL / g->dt);
    int blown = 0;

    for (int step = 0; step < n_steps; step++) {
        verlet_full_step(g, params, -1);
        if (step % DAMP_EVERY == 0)
            apply_damping_xy(g);
        if (check_blowup(g)) { blown = 1; break; }
    }

    compute_diagnostics(g, params, res);
    snprintf(res->label, sizeof(res->label), "%s", label);
    res->stable = !blown;
}

int main(void) {
    printf("=== T6: Universality — Bimodal Synergy Across Parameter Space ===\n");
    printf("N=%d, L=%.1f, T=%.0f, mass=1.50 (m^2=2.25)\n", N_GRID, L_BOX, T_EVOL);
    printf("Scanning %d (mu,kappa) pairs x 3 configs = %d runs\n\n",
           N_PAIRS, N_PAIRS * 3);

    bimodal_init_params();

    Grid *g = grid_alloc(N_GRID, L_BOX);

    PairResult pr[N_PAIRS];
    double t_wall0 = omp_get_wtime();

    /* Header */
    printf("%-6s %-6s | %-8s %-8s | %-8s %-8s | %-8s %-8s | %-7s | status\n",
           "mu", "kappa",
           "A_trans", "A_tor",
           "B_trans", "B_tor",
           "bm_trans", "bm_tor",
           "synergy");
    printf("------ ------ | -------- -------- | -------- -------- | -------- -------- | ------- | ------\n");
    fflush(stdout);

    for (int ip = 0; ip < N_PAIRS; ip++) {
        double mu = MU_LIST[ip];
        double kap = KAPPA_LIST[ip];
        pr[ip].mu = mu;
        pr[ip].kappa = kap;

        double pA[NDIM], pB[NDIM], pBM[NDIM];
        Result rA, rB, rBM;

        /* --- Pure A --- */
        memcpy(pA, PATH_A, sizeof(pA));
        pA[12] = mu;  pA[13] = kap;  pA[14] = 1.50;
        run_one(pA, "A", &rA, g);

        /* --- Pure B --- */
        memcpy(pB, PATH_B, sizeof(pB));
        pB[12] = mu;  pB[13] = kap;  pB[14] = 1.50;
        run_one(pB, "B", &rB, g);

        /* --- Bimodal --- */
        memcpy(pBM, BIMODAL, sizeof(pBM));
        pBM[12] = mu;  pBM[13] = kap;  pBM[14] = 1.50;
        run_one(pBM, "BM", &rBM, g);

        pr[ip].A_trans = rA.transverse_l2;
        pr[ip].A_tor   = rA.torsion_flux;
        pr[ip].B_trans = rB.transverse_l2;
        pr[ip].B_tor   = rB.torsion_flux;
        pr[ip].bm_trans = rBM.transverse_l2;
        pr[ip].bm_tor   = rBM.torsion_flux;
        pr[ip].A_ok  = rA.stable;
        pr[ip].B_ok  = rB.stable;
        pr[ip].bm_ok = rBM.stable;

        /* Synergy: bimodal must beat both controls */
        double eps_floor = 0.001;
        double s_trans = pr[ip].bm_trans / fmax(pr[ip].A_trans, eps_floor);
        double s_tor   = pr[ip].bm_tor   / fmax(pr[ip].B_tor,  eps_floor);
        pr[ip].synergy = fmin(s_trans, s_tor);

        /* If any config blew up, mark synergy as 0 */
        if (!rBM.stable) pr[ip].synergy = 0.0;

        char status[32];
        snprintf(status, sizeof(status), "%s%s%s",
                 rA.stable ? "A" : "a",
                 rB.stable ? "B" : "b",
                 rBM.stable ? "M" : "m");

        printf("%6.1f %6.1f | %8.4f %8.4f | %8.4f %8.4f | %8.4f %8.4f | %7.3f | %s\n",
               mu, kap,
               pr[ip].A_trans, pr[ip].A_tor,
               pr[ip].B_trans, pr[ip].B_tor,
               pr[ip].bm_trans, pr[ip].bm_tor,
               pr[ip].synergy, status);
        fflush(stdout);
    }

    double t_total = omp_get_wtime() - t_wall0;

    /* === Summary === */
    printf("\n=== Summary ===\n");
    int n_synergy = 0, n_total_ok = 0;
    for (int ip = 0; ip < N_PAIRS; ip++) {
        if (pr[ip].bm_ok) n_total_ok++;
        if (pr[ip].synergy > 1.0) n_synergy++;
    }
    printf("Bimodal stable: %d / %d\n", n_total_ok, N_PAIRS);
    printf("Synergy > 1.0:  %d / %d\n", n_synergy, N_PAIRS);
    printf("Total wall time: %.1f sec\n", t_total);

    /* Best and worst synergy */
    int best = 0, worst = 0;
    for (int ip = 1; ip < N_PAIRS; ip++) {
        if (pr[ip].synergy > pr[best].synergy) best = ip;
        if (pr[ip].synergy < pr[worst].synergy) worst = ip;
    }
    printf("Best synergy:  %.3f at mu=%.1f, kappa=%.1f\n",
           pr[best].synergy, pr[best].mu, pr[best].kappa);
    printf("Worst synergy: %.3f at mu=%.1f, kappa=%.1f\n",
           pr[worst].synergy, pr[worst].mu, pr[worst].kappa);

    if (n_synergy >= 7)
        printf("\nVERDICT: UNIVERSAL — synergy persists across >= 70%% of parameter space\n");
    else if (n_synergy >= 4)
        printf("\nVERDICT: BROAD — synergy in %d/10 points (moderate region)\n", n_synergy);
    else if (n_synergy >= 1)
        printf("\nVERDICT: NARROW — synergy only at %d/10 points\n", n_synergy);
    else
        printf("\nVERDICT: NO SYNERGY — bimodal does not consistently beat controls\n");

    /* === Write TSV === */
    FILE *fp = fopen("data/t6_results.tsv", "w");
    if (fp) {
        fprintf(fp, "mu\tkappa\tA_trans\tA_tor\tB_trans\tB_tor\tbm_trans\tbm_tor\tsynergy\tA_ok\tB_ok\tbm_ok\n");
        for (int ip = 0; ip < N_PAIRS; ip++) {
            fprintf(fp, "%.1f\t%.1f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.4f\t%d\t%d\t%d\n",
                    pr[ip].mu, pr[ip].kappa,
                    pr[ip].A_trans, pr[ip].A_tor,
                    pr[ip].B_trans, pr[ip].B_tor,
                    pr[ip].bm_trans, pr[ip].bm_tor,
                    pr[ip].synergy,
                    pr[ip].A_ok, pr[ip].B_ok, pr[ip].bm_ok);
        }
        fclose(fp);
        printf("\nResults written to data/t6_results.tsv\n");
    } else {
        printf("\nWARNING: Could not write data/t6_results.tsv\n");
    }

    grid_free(g);
    printf("\nDone.\n");
    return 0;
}
