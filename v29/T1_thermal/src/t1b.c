/*  t1b.c — Thermal bath equilibrium test WITH correct mass AND pre-equilibration
 *
 *  Phase 1: Pre-equilibrate T=300 with absorbing BC to shed init radiation
 *  Phase 2: Switch to fully periodic BC + noise, evolve T=700 more
 *
 *  Key difference from t1: mass²=2.25 in dynamics, pre-equilibration phase.
 *
 *  Build: gcc -O3 -fopenmp -o t1b src/t1b.c -lm
 */

#include "../../src/braid_core.h"

/* ================================================================
   Box-Muller Gaussian RNG (deterministic, seeded)
   ================================================================ */

static unsigned long rng_state;

static void rng_seed(unsigned long s) { rng_state = s; }

static double rng_uniform(void) {
    rng_state = rng_state * 6364136223846793005ULL + 1442695040888963407ULL;
    return (double)(rng_state >> 11) / (double)(1ULL << 53);
}

static double rng_gaussian(void) {
    double u1 = rng_uniform() + 1e-30;
    double u2 = rng_uniform();
    return sqrt(-2.0 * log(u1)) * cos(2.0 * PI * u2);
}

/* ================================================================
   Fully periodic force computation (all 3 directions wrap)
   ================================================================ */

static void compute_forces_allperiodic(Grid *g, const double *phys, double mass2_override) {
    int N = g->N, NN = N*N, N3 = N*N*N;
    double idx2 = 1.0 / (g->dx * g->dx);
    double mu    = phys[12];
    double kappa = phys[13];
    double mass2 = (mass2_override >= 0) ? mass2_override : phys[14]*phys[14];
    double lpw   = phys[15];

    #pragma omp parallel for schedule(static)
    for (int idx = 0; idx < N3; idx++) {
        int i = idx / NN, j = (idx / N) % N, k = idx % N;

        /* Periodic neighbors in all 3 directions */
        int ip = (i+1)%N, im = (i-1+N)%N;
        int jp = (j+1)%N, jm = (j-1+N)%N;
        int kp = (k+1)%N, km = (k-1+N)%N;

        int idx_ip = ip*NN + j*N + k;
        int idx_im = im*NN + j*N + k;
        int idx_jp = i*NN + jp*N + k;
        int idx_jm = i*NN + jm*N + k;
        int idx_kp = i*NN + j*N + kp;
        int idx_km = i*NN + j*N + km;

        double p0 = g->phi[0][idx], p1 = g->phi[1][idx], p2 = g->phi[2][idx];
        double P  = p0 * p1 * p2;
        double denom = 1.0 + kappa * P * P;
        double mu_P_d2 = mu * P / (denom * denom);

        for (int a = 0; a < NFIELDS; a++) {
            double lap = (g->phi[a][idx_ip] + g->phi[a][idx_im]
                        + g->phi[a][idx_jp] + g->phi[a][idx_jm]
                        + g->phi[a][idx_kp] + g->phi[a][idx_km]
                        - 6.0 * g->phi[a][idx]) * idx2;
            double dPda = (a==0) ? p1*p2 : (a==1) ? p0*p2 : p0*p1;
            double f_triple = mu_P_d2 * dPda;
            double f_pw = lpw * (g->phi[(a+1)%3][idx] + g->phi[(a+2)%3][idx]);
            g->acc[a][idx] = lap - mass2 * g->phi[a][idx] - f_triple - f_pw;
        }
    }
}

/* Full Verlet step with all-periodic forces */
static void verlet_full_step_periodic(Grid *g, const double *phys, double mass2_ov) {
    double hdt = 0.5 * g->dt;
    verlet_kick(g, hdt);
    verlet_drift(g);
    compute_forces_allperiodic(g, phys, mass2_ov);
    verlet_kick(g, hdt);
}

/* ================================================================
   Grid copy (save settled state for reuse)
   ================================================================ */

static void grid_copy(Grid *dst, const Grid *src) {
    int N3 = src->N * src->N * src->N;
    dst->N = src->N; dst->L = src->L; dst->dx = src->dx; dst->dt = src->dt;
    for (int a = 0; a < NFIELDS; a++) {
        memcpy(dst->phi[a], src->phi[a], N3 * sizeof(double));
        memcpy(dst->vel[a], src->vel[a], N3 * sizeof(double));
        memcpy(dst->acc[a], src->acc[a], N3 * sizeof(double));
    }
}

/* ================================================================
   Main
   ================================================================ */

int main(void) {
    int N = 96;
    double L = 20.0;

    bimodal_init_params();

    printf("=== T1b: Thermal Bath with mass²=2.25 + Pre-equilibration ===\n");
    printf("N=%d  L=%.1f  dx=%.4f  dt=%.5f\n", N, L, 2.0*L/(N-1), 0.20*2.0*L/(N-1));
    printf("Bimodal params: mu=%.1f  kappa=%.1f  mass=%.2f  mass²=%.2f\n",
           BIMODAL[12], BIMODAL[13], BIMODAL[14], BIMODAL[14]*BIMODAL[14]);
    printf("\n");

    /* ---- Phase 1: Pre-equilibrate with absorbing BC ---- */
    printf("=== PHASE 1: Pre-equilibration with absorbing BC (T=0..300) ===\n");
    fflush(stdout);

    Grid *g_base = grid_alloc(N, L);
    init_braid(g_base, BIMODAL, -1);   /* mass from BIMODAL = 1.50 */
    compute_forces(g_base, BIMODAL, -1);  /* initial forces (absorbing: Dirichlet x,y) */

    double T_pre = 300.0;
    double dt = g_base->dt;
    int steps_pre = (int)(T_pre / dt + 0.5);
    int diag_interval_pre = (int)(50.0 / dt + 0.5);

    Result res;
    for (int s = 0; s <= steps_pre; s++) {
        if (s % diag_interval_pre == 0) {
            double T = s * dt;
            compute_diagnostics(g_base, BIMODAL, &res);
            printf("  T=%6.1f  E=%10.1f  fc=%.4f  |P|=%.4f  w=%+.3f\n",
                   T, res.energy, res.fc, res.peak_P, res.winding);
            fflush(stdout);
        }
        if (s < steps_pre) {
            verlet_full_step(g_base, BIMODAL, -1);  /* mass²=2.25 */
            apply_damping_xy(g_base);
        }
        if (check_blowup(g_base)) {
            printf("  BLOWUP during pre-equilibration at step %d!\n", s);
            return 1;
        }
    }

    /* Record settled state */
    compute_diagnostics(g_base, BIMODAL, &res);
    double E_settled = res.energy, fc_settled = res.fc, P_settled = res.peak_P;
    double w_settled = res.winding;
    printf("\nSettled state: E=%.1f  fc=%.4f  |P|=%.4f  w=%+.3f\n\n",
           E_settled, fc_settled, P_settled, w_settled);
    fflush(stdout);

    /* Save settled state */
    Grid *g_saved = grid_alloc(N, L);
    grid_copy(g_saved, g_base);

    /* ---- Phase 2: Periodic BC + noise for each noise level ---- */
    double noise_levels[] = {0.0, 0.001, 0.005, 0.01, 0.05, 0.1};
    int n_noise = 6;

    double T_phase2 = 700.0;
    int steps_p2 = (int)(T_phase2 / dt + 0.5);
    int diag_interval = (int)(50.0 / dt + 0.5);
    int n_diag = (int)(T_phase2 / 50.0) + 1;

    /* Summary arrays */
    double late_fc[6], late_P[6], late_dEdt[6];

    printf("=== PHASE 2: Periodic BC + noise (T=300..1000) ===\n\n");
    fflush(stdout);

    for (int ni = 0; ni < n_noise; ni++) {
        double A_noise = noise_levels[ni];
        printf("--- Noise amplitude = %.4f ---\n", A_noise);
        fflush(stdout);

        /* Copy settled state */
        grid_copy(g_base, g_saved);

        /* Add noise */
        int N3 = N * N * N;
        rng_seed(42 + ni * 12345);
        for (int a = 0; a < NFIELDS; a++) {
            for (int idx = 0; idx < N3; idx++) {
                g_base->phi[a][idx] += A_noise * rng_gaussian();
                g_base->vel[a][idx] += A_noise * 0.5 * rng_gaussian();
            }
        }

        /* Initial forces with all-periodic */
        compute_forces_allperiodic(g_base, BIMODAL, -1);

        /* Open data file */
        char fname[256];
        snprintf(fname, sizeof(fname), "data/t1b_noise_%.4f.tsv", A_noise);
        FILE *fp = fopen(fname, "w");
        if (!fp) { printf("Cannot open %s\n", fname); return 1; }
        fprintf(fp, "T\tE\tfc\tpeak_P\twinding\n");

        /* Time series storage for late-time averaging */
        double sum_fc = 0, sum_P = 0, n_late = 0;
        double E_first_late = 0, E_last_late = 0;
        double T_first_late = 0, T_last_late = 0;

        for (int s = 0; s <= steps_p2; s++) {
            if (s % diag_interval == 0) {
                double T = T_pre + s * dt;
                compute_diagnostics(g_base, BIMODAL, &res);
                printf("  T=%6.1f  E=%10.1f  fc=%.4f  |P|=%.4f  w=%+.3f\n",
                       T, res.energy, res.fc, res.peak_P, res.winding);
                fprintf(fp, "%.1f\t%.4f\t%.6f\t%.6f\t%.4f\n",
                        T, res.energy, res.fc, res.peak_P, res.winding);
                fflush(stdout);

                /* Late-time (T > 600 total, i.e. s*dt > 300 in phase 2) */
                if (T > 600.0) {
                    sum_fc += res.fc;
                    sum_P += res.peak_P;
                    n_late += 1.0;
                    if (n_late == 1.0) {
                        E_first_late = res.energy;
                        T_first_late = T;
                    }
                    E_last_late = res.energy;
                    T_last_late = T;
                }
            }
            if (s < steps_p2) {
                verlet_full_step_periodic(g_base, BIMODAL, -1);  /* mass²=2.25, all periodic */
            }
            if (check_blowup(g_base)) {
                printf("  BLOWUP at step %d (T=%.1f)!\n", s, T_pre + s*dt);
                late_fc[ni] = -1; late_P[ni] = -1; late_dEdt[ni] = 0;
                fprintf(fp, "# BLOWUP\n");
                fclose(fp);
                goto next_noise;
            }
        }

        fclose(fp);

        /* Summary */
        if (n_late > 0) {
            late_fc[ni] = sum_fc / n_late;
            late_P[ni] = sum_P / n_late;
            late_dEdt[ni] = (T_last_late > T_first_late) ?
                (E_last_late - E_first_late) / (T_last_late - T_first_late) : 0.0;
        } else {
            late_fc[ni] = res.fc;
            late_P[ni] = res.peak_P;
            late_dEdt[ni] = 0;
        }
        printf("  Late-time avg: fc=%.4f  |P|=%.4f  dE/dt=%.2f\n\n",
               late_fc[ni], late_P[ni], late_dEdt[ni]);
        fflush(stdout);

        next_noise: ;
    }

    /* ---- Final summary ---- */
    printf("\n=== SUMMARY (late-time T>600 averages) ===\n");
    printf("  %-10s  %-10s  %-10s  %-10s\n", "A_noise", "avg_fc", "avg_|P|", "dE/dt");
    printf("  %-10s  %-10s  %-10s  %-10s\n", "-------", "------", "-------", "-----");
    for (int ni = 0; ni < n_noise; ni++) {
        if (late_fc[ni] < 0)
            printf("  %-10.4f  BLOWUP\n", noise_levels[ni]);
        else
            printf("  %-10.4f  %-10.4f  %-10.4f  %-10.2f\n",
                   noise_levels[ni], late_fc[ni], late_P[ni], late_dEdt[ni]);
    }
    printf("\nSettled (end of Phase 1): E=%.1f  fc=%.4f  |P|=%.4f  w=%+.3f\n",
           E_settled, fc_settled, P_settled, w_settled);

    grid_free(g_base);
    grid_free(g_saved);

    return 0;
}
