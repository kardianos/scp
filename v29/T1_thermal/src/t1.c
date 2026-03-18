/*  T1: Thermal bath equilibrium test
 *
 *  Tests whether the bimodal braid reaches dynamic equilibrium in a
 *  radiation bath (periodic BC, no absorbing layer).
 *
 *  Scans noise levels A_noise in {0, 0.001, 0.005, 0.01, 0.05, 0.1}.
 *  Runs to T=1000 with diagnostics every T=50.
 */

#include <stdint.h>
#include "../../src/braid_core.h"

/* ================================================================
   Fully periodic force computation (replaces compute_forces)
   ================================================================ */

static void compute_forces_periodic(Grid *g, const double *phys, double mass2_override) {
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

/* Full Verlet step with periodic forces */
static void verlet_full_step_periodic(Grid *g, const double *phys, double mass2_ov) {
    double hdt = 0.5 * g->dt;
    verlet_kick(g, hdt);
    verlet_drift(g);
    compute_forces_periodic(g, phys, mass2_ov);
    verlet_kick(g, hdt);
}

/* ================================================================
   RNG: splitmix64-based Gaussian
   ================================================================ */

static uint64_t splitmix64(uint64_t *state) {
    uint64_t z = (*state += 0x9e3779b97f4a7c15ULL);
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}

static double uniform01(uint64_t *state) {
    return (splitmix64(state) >> 11) * (1.0 / 9007199254740992.0);
}

/* Box-Muller */
static double gaussian_random(uint64_t *state) {
    double u1 = uniform01(state);
    double u2 = uniform01(state);
    if (u1 < 1e-30) u1 = 1e-30;
    return sqrt(-2.0 * log(u1)) * cos(2.0 * PI * u2);
}

/* ================================================================
   Far-field probe: average |phi| at 4 corner points
   ================================================================ */

static double probe_farfield(Grid *g) {
    int N = g->N, NN = N*N;
    double dx = g->dx, L = g->L;
    /* Probe points: (+-8, +-8, 0) */
    double probes[4][2] = {{8,8},{-8,8},{8,-8},{-8,-8}};
    double sum = 0;
    for (int p = 0; p < 4; p++) {
        int i = (int)((probes[p][0] + L) / dx + 0.5);
        int j = (int)((probes[p][1] + L) / dx + 0.5);
        int k = N / 2;  /* z=0 */
        if (i < 0) i = 0; if (i >= N) i = N-1;
        if (j < 0) j = 0; if (j >= N) j = N-1;
        int idx = i*NN + j*N + k;
        double amp2 = 0;
        for (int a = 0; a < NFIELDS; a++)
            amp2 += g->phi[a][idx] * g->phi[a][idx];
        sum += sqrt(amp2);
    }
    return sum / 4.0;
}

/* ================================================================
   Main
   ================================================================ */

int main(void) {
    bimodal_init_params();

    int N = 96;
    double L = 20.0;
    double T_end = 1000.0;
    double T_diag = 50.0;
    double omega_avg = 0.5;

    double noise_levels[] = {0.0, 0.001, 0.005, 0.01, 0.05, 0.1};
    int n_noise = 6;
    int n_snap = (int)(T_end / T_diag) + 1;  /* 21 snapshots (t=0..1000) */

    printf("T1: Thermal Bath Equilibrium\n");
    printf("N=%d, L=%.1f, T_end=%.0f, T_diag=%.0f\n", N, L, T_end, T_diag);

    /* Temporary grid to get dt */
    Grid *tmp = grid_alloc(N, L);
    double dt = tmp->dt;
    double dx = tmp->dx;
    grid_free(tmp);
    printf("dx=%.4f, dt=%.5f\n\n", dx, dt);

    /* Storage for late-time averages */
    double late_fc[6], late_P[6], late_wind[6], late_E[6];
    int late_count[6];
    double E_first_late[6], E_last_late[6];
    double T_first_late[6], T_last_late[6];

    for (int ni = 0; ni < n_noise; ni++) {
        double A_noise = noise_levels[ni];
        printf("=== Noise level A = %.4f ===\n", A_noise);

        Grid *g = grid_alloc(N, L);
        init_braid(g, BIMODAL, -1);

        /* Add noise */
        int N3 = N*N*N;
        for (int a = 0; a < NFIELDS; a++) {
            for (int idx = 0; idx < N3; idx++) {
                uint64_t seed = (uint64_t)ni * 1000003ULL + (uint64_t)a * 100003ULL + (uint64_t)idx;
                g->phi[a][idx] += A_noise * gaussian_random(&seed);
                seed += 7777777ULL;
                g->vel[a][idx] += A_noise * omega_avg * gaussian_random(&seed);
            }
        }

        /* Initial force */
        compute_forces_periodic(g, BIMODAL, 0.0);

        /* Open TSV file */
        char fname[128];
        snprintf(fname, sizeof(fname), "data/t1_noise_%.4f.tsv", A_noise);
        FILE *fp = fopen(fname, "w");
        fprintf(fp, "time\tE_total\tfc\tpeak_P\twinding\tff_amp\n");

        /* Time loop */
        double t = 0;
        int step = 0;
        int steps_per_diag = (int)(T_diag / dt + 0.5);
        int total_steps = (int)(T_end / dt + 0.5);

        late_fc[ni] = late_P[ni] = late_wind[ni] = late_E[ni] = 0;
        late_count[ni] = 0;
        E_first_late[ni] = E_last_late[ni] = 0;
        T_first_late[ni] = T_last_late[ni] = 0;

        int snap = 0;
        double t_wall_start = omp_get_wtime();

        /* Snapshot at t=0 */
        {
            Result res;
            compute_diagnostics(g, BIMODAL, &res);
            double ff = probe_farfield(g);
            printf("  t=%7.1f  E=%10.1f  fc=%.4f  |P|=%.4f  w=%+.3f  ff=%.4f\n",
                   0.0, res.energy, res.fc, res.peak_P, res.winding, ff);
            fprintf(fp, "%.1f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n",
                    0.0, res.energy, res.fc, res.peak_P, res.winding, ff);
            fflush(fp);
        }

        for (step = 1; step <= total_steps; step++) {
            verlet_full_step_periodic(g, BIMODAL, 0.0);
            t = step * dt;

            if (step % steps_per_diag == 0) {
                snap++;
                if (check_blowup(g)) {
                    printf("  BLOWUP at t=%.1f\n", t);
                    break;
                }
                Result res;
                compute_diagnostics(g, BIMODAL, &res);
                double ff = probe_farfield(g);
                printf("  t=%7.1f  E=%10.1f  fc=%.4f  |P|=%.4f  w=%+.3f  ff=%.4f\n",
                       t, res.energy, res.fc, res.peak_P, res.winding, ff);
                fprintf(fp, "%.1f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n",
                        t, res.energy, res.fc, res.peak_P, res.winding, ff);
                fflush(fp);

                /* Accumulate late-time stats (T > 500) */
                if (t > 500.0) {
                    late_fc[ni] += res.fc;
                    late_P[ni] += res.peak_P;
                    late_wind[ni] += res.winding;
                    late_E[ni] += res.energy;
                    late_count[ni]++;
                    if (E_first_late[ni] == 0) {
                        E_first_late[ni] = res.energy;
                        T_first_late[ni] = t;
                    }
                    E_last_late[ni] = res.energy;
                    T_last_late[ni] = t;
                }
            }
        }
        fclose(fp);

        double t_wall = omp_get_wtime() - t_wall_start;
        printf("  Wall time: %.1f s\n\n", t_wall);

        grid_free(g);
    }

    /* ================================================================
       Summary
       ================================================================ */
    printf("\n========================================\n");
    printf("SUMMARY: Late-time averages (T > 500)\n");
    printf("========================================\n");
    printf("%-10s %10s %8s %8s %8s %12s %10s\n",
           "A_noise", "fc", "|P|", "winding", "dE/E/T", "E_avg", "stable?");
    printf("%-10s %10s %8s %8s %8s %12s %10s\n",
           "-------", "------", "------", "-------", "------", "-------", "-------");

    for (int ni = 0; ni < n_noise; ni++) {
        double A_noise = noise_levels[ni];
        if (late_count[ni] == 0) {
            printf("%-10.4f  %s\n", A_noise, "BLOWUP before T=500");
            continue;
        }
        double avg_fc = late_fc[ni] / late_count[ni];
        double avg_P  = late_P[ni]  / late_count[ni];
        double avg_w  = late_wind[ni] / late_count[ni];
        double avg_E  = late_E[ni]  / late_count[ni];
        double dT = T_last_late[ni] - T_first_late[ni];
        double dE_rel = (dT > 0) ? (E_last_late[ni] - E_first_late[ni]) / (fabs(avg_E) * dT + 1e-30) : 0;
        int stable = (avg_fc > 0.3) && (fabs(dE_rel) < 0.01);
        printf("%-10.4f %10.4f %8.4f %8.3f %8.5f %12.1f %10s\n",
               A_noise, avg_fc, avg_P, avg_w, dE_rel, avg_E,
               stable ? "YES" : "NO");
    }

    printf("\nSuccess criteria: fc > 0.3, |dE/E/T| < 0.01\n");
    printf("Done.\n");
    return 0;
}
