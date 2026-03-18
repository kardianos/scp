/*  T5: Self-Consistent Metric Coupling
 *
 *  The fields propagate on the metric they create:
 *    g^{ij} = delta^{ij} - alpha_g * h_{ij}
 *    h_{ij} = d_i phi_j + d_j phi_i   (symmetrized strain)
 *
 *  Modified EOM:
 *    d^2 phi_a / dt^2 = g^{ij} d_i d_j phi_a - m^2 phi_a - dV/dphi_a
 *
 *  Scan alpha_g in {0, 0.0001, 0.001, 0.005, 0.01, 0.05}.
 *  N=96, L=20, T=300. Mass=1.50, m^2=2.25.
 *
 *  Build: cd v29/T5_metric && gcc -O3 -fopenmp -o t5 src/t5.c -lm
 */

#include "../../src/braid_core.h"

/* ================================================================
   Modified force computation with metric backreaction
   ================================================================ */

static void compute_forces_metric(Grid *g, const double *phys,
                                  double mass2_override, double alpha_g)
{
    int N = g->N, NN = N*N, N3 = N*N*N;
    double dx = g->dx;
    double idx2 = 1.0 / (dx * dx);
    double idx4 = 1.0 / (4.0 * dx * dx);  /* for mixed partials */
    double idx1 = 1.0 / (2.0 * dx);       /* for first derivatives */
    double mu    = phys[12];
    double kappa = phys[13];
    double mass2 = (mass2_override >= 0) ? mass2_override : phys[14]*phys[14];
    double lpw   = phys[15];

    #pragma omp parallel for schedule(static)
    for (int idx = 0; idx < N3; idx++) {
        int i = idx / NN, j = (idx / N) % N, k = idx % N;

        /* Need 2 cells of interior for mixed partials */
        if (i < 2 || i >= N-2 || j < 2 || j >= N-2) {
            g->acc[0][idx] = g->acc[1][idx] = g->acc[2][idx] = 0;
            continue;
        }

        int kp  = (k+1) % N, km  = (k-1+N) % N;
        int kp2 = (k+2) % N, km2 = (k-2+N) % N;

        int idx_kp = i*NN + j*N + kp;
        int idx_km = i*NN + j*N + km;

        /* Potential terms */
        double p0 = g->phi[0][idx], p1 = g->phi[1][idx], p2 = g->phi[2][idx];
        double P  = p0 * p1 * p2;
        double denom = 1.0 + kappa * P * P;
        double mu_P_d2 = mu * P / (denom * denom);

        /* If alpha_g == 0, skip metric computation entirely */
        if (alpha_g == 0.0) {
            for (int a = 0; a < NFIELDS; a++) {
                double lap = (g->phi[a][idx+NN] + g->phi[a][idx-NN]
                            + g->phi[a][idx+N]  + g->phi[a][idx-N]
                            + g->phi[a][idx_kp] + g->phi[a][idx_km]
                            - 6.0 * g->phi[a][idx]) * idx2;
                double dPda = (a==0) ? p1*p2 : (a==1) ? p0*p2 : p0*p1;
                double f_triple = mu_P_d2 * dPda;
                double f_pw = lpw * (g->phi[(a+1)%3][idx] + g->phi[(a+2)%3][idx]);
                g->acc[a][idx] = lap - mass2 * g->phi[a][idx] - f_triple - f_pw;
            }
            continue;
        }

        /* ---- Compute gradient tensor G_{ia} = d_i phi_a ---- */
        /* i=0 -> x, i=1 -> y, i=2 -> z; a=field index */
        double G[3][3];
        for (int a = 0; a < 3; a++) {
            G[0][a] = (g->phi[a][idx+NN] - g->phi[a][idx-NN]) * idx1;
            G[1][a] = (g->phi[a][idx+N]  - g->phi[a][idx-N])  * idx1;
            G[2][a] = (g->phi[a][idx_kp] - g->phi[a][idx_km])  * idx1;
        }

        /* ---- Strain tensor h_{ij} = G_{ij} + G_{ji} ---- */
        /* Note: h_{ij} is 3x3 symmetric, where G_{ij} means d_i phi_j
           For NFIELDS=3, field index a maps to spatial index a */
        double h[3][3];
        for (int ii = 0; ii < 3; ii++)
            for (int jj = ii; jj < 3; jj++) {
                h[ii][jj] = G[ii][jj] + G[jj][ii];
                h[jj][ii] = h[ii][jj];
            }

        /* ---- For each field a, compute modified Laplacian ---- */
        for (int a = 0; a < NFIELDS; a++) {
            /* Standard Laplacian */
            double lap = (g->phi[a][idx+NN] + g->phi[a][idx-NN]
                        + g->phi[a][idx+N]  + g->phi[a][idx-N]
                        + g->phi[a][idx_kp] + g->phi[a][idx_km]
                        - 6.0 * g->phi[a][idx]) * idx2;

            /* Metric correction: delta_lap = -alpha_g * h_{ij} * d^2 phi_a / dx_i dx_j */

            /* Diagonal second derivatives d^2 phi_a / dx_i^2 */
            double d2xx = (g->phi[a][idx+NN] - 2.0*g->phi[a][idx] + g->phi[a][idx-NN]) * idx2;
            double d2yy = (g->phi[a][idx+N]  - 2.0*g->phi[a][idx] + g->phi[a][idx-N])  * idx2;
            double d2zz = (g->phi[a][idx_kp] - 2.0*g->phi[a][idx] + g->phi[a][idx_km]) * idx2;

            /* Off-diagonal mixed second derivatives */
            /* d^2 phi_a / dxdy */
            double d2xy = (g->phi[a][(i+1)*NN + (j+1)*N + k]
                         - g->phi[a][(i+1)*NN + (j-1)*N + k]
                         - g->phi[a][(i-1)*NN + (j+1)*N + k]
                         + g->phi[a][(i-1)*NN + (j-1)*N + k]) * idx4;

            /* d^2 phi_a / dxdz */
            double d2xz = (g->phi[a][(i+1)*NN + j*N + kp]
                         - g->phi[a][(i+1)*NN + j*N + km]
                         - g->phi[a][(i-1)*NN + j*N + kp]
                         + g->phi[a][(i-1)*NN + j*N + km]) * idx4;

            /* d^2 phi_a / dydz */
            double d2yz = (g->phi[a][i*NN + (j+1)*N + kp]
                         - g->phi[a][i*NN + (j+1)*N + km]
                         - g->phi[a][i*NN + (j-1)*N + kp]
                         + g->phi[a][i*NN + (j-1)*N + km]) * idx4;

            /* Metric correction: -alpha_g * h_{ij} * d2_{ij} */
            double delta_lap = -alpha_g * (
                h[0][0] * d2xx + h[1][1] * d2yy + h[2][2] * d2zz
              + 2.0 * (h[0][1] * d2xy + h[0][2] * d2xz + h[1][2] * d2yz)
            );

            /* Potential forces */
            double dPda = (a==0) ? p1*p2 : (a==1) ? p0*p2 : p0*p1;
            double f_triple = mu_P_d2 * dPda;
            double f_pw = lpw * (g->phi[(a+1)%3][idx] + g->phi[(a+2)%3][idx]);

            g->acc[a][idx] = lap + delta_lap - mass2 * g->phi[a][idx] - f_triple - f_pw;
        }
    }
}

/* ================================================================
   Verlet step using metric force
   ================================================================ */

static void verlet_metric_step(Grid *g, const double *phys,
                               double mass2_ov, double alpha_g)
{
    double hdt = 0.5 * g->dt;
    verlet_kick(g, hdt);
    verlet_drift(g);
    compute_forces_metric(g, phys, mass2_ov, alpha_g);
    verlet_kick(g, hdt);
}

/* ================================================================
   Far-field h_{ij} decomposition on spherical shell
   Golden-ratio spiral sampling, l=0 (trace) vs l=2 (deviatoric)
   ================================================================ */

static double compute_l2_metric(Grid *g, double R_shell) {
    int N = g->N, NN = N*N;
    double dx = g->dx, L = g->L;
    double idx1 = 1.0 / (2.0 * dx);

    int N_pts = 50;
    double golden_angle = PI * (3.0 - sqrt(5.0));

    double trace2_sum = 0.0, dev2_sum = 0.0;

    for (int n = 0; n < N_pts; n++) {
        /* Golden-ratio spiral point on unit sphere */
        double theta = acos(1.0 - 2.0 * (n + 0.5) / N_pts);
        double phi_ang = golden_angle * n;

        double x = R_shell * sin(theta) * cos(phi_ang);
        double y = R_shell * sin(theta) * sin(phi_ang);
        double z = R_shell * cos(theta);

        /* Map to grid indices */
        double fi = (x + L) / dx;
        double fj = (y + L) / dx;
        double fk = (z + L) / dx;

        int i = (int)fi, j = (int)fj, k = (int)fk;

        /* Bounds check (need +-1 for gradients) */
        if (i < 2 || i >= N-2 || j < 2 || j >= N-2) continue;
        /* z is periodic, but skip if out of main domain */
        if (k < 0 || k >= N) continue;

        int kp = (k+1)%N, km = (k-1+N)%N;
        int idx_c = i*NN + j*N + k;
        int idx_kp = i*NN + j*N + kp;
        int idx_km = i*NN + j*N + km;

        /* Compute gradient tensor G_{ia} = d_i phi_a at this point */
        double G[3][3];
        for (int a = 0; a < 3; a++) {
            G[0][a] = (g->phi[a][idx_c+NN] - g->phi[a][idx_c-NN]) * idx1;
            G[1][a] = (g->phi[a][idx_c+N]  - g->phi[a][idx_c-N])  * idx1;
            G[2][a] = (g->phi[a][idx_kp]   - g->phi[a][idx_km])    * idx1;
        }

        /* Strain h_{ij} = G_{ij} + G_{ji} */
        double h[3][3];
        for (int ii = 0; ii < 3; ii++)
            for (int jj = ii; jj < 3; jj++) {
                h[ii][jj] = G[ii][jj] + G[jj][ii];
                h[jj][ii] = h[ii][jj];
            }

        /* Trace (l=0 part) */
        double tr = h[0][0] + h[1][1] + h[2][2];
        trace2_sum += tr * tr;

        /* Deviatoric (traceless = l>=2 part) */
        double dev2 = 0.0;
        for (int ii = 0; ii < 3; ii++)
            for (int jj = 0; jj < 3; jj++) {
                double dij = (ii == jj) ? 1.0 : 0.0;
                double dev_ij = h[ii][jj] - (tr / 3.0) * dij;
                dev2 += dev_ij * dev_ij;
            }
        dev2_sum += dev2;
    }

    /* l2_h = |deviatoric|^2 / (|trace|^2 + |deviatoric|^2) */
    double total = trace2_sum + dev2_sum;
    if (total < 1e-30) return 0.0;
    return dev2_sum / total;
}

/* ================================================================
   MAIN
   ================================================================ */

int main(void) {
    printf("=== T5: Self-Consistent Metric Coupling ===\n\n");

    bimodal_init_params();
    printf("BIMODAL params:\n");
    for (int d = 0; d < NDIM; d++)
        printf("  %8s = %8.4f\n", PNAME[d], BIMODAL[d]);

    double mass = 1.50;
    double mass2 = mass * mass;  /* 2.25 */
    printf("\nmass = %.4f => mass^2 = %.4f\n", mass, mass2);

    int N_grid = 96;
    double L = 20.0;
    double dx = 2.0*L/(N_grid-1);
    printf("Grid: N=%d, L=%.1f, dx=%.4f\n", N_grid, L, dx);

    double alpha_g_vals[] = {0.0, 0.0001, 0.001, 0.005, 0.01, 0.05};
    int n_alpha = 6;
    double T_max = 300.0;
    double R_shell = 10.0;

    /* Results storage */
    double res_fc[6], res_trans[6], res_torsion[6], res_l2h[6];
    int    res_stable[6];
    double res_l2[6], res_energy[6];

    printf("T_max = %.0f, R_shell = %.1f for l=2 decomposition\n", T_max, R_shell);
    printf("alpha_g values: ");
    for (int s = 0; s < n_alpha; s++) printf("%.4f ", alpha_g_vals[s]);
    printf("\n\n");

    double wall_total_start = omp_get_wtime();

    for (int s = 0; s < n_alpha; s++) {
        double alpha_g = alpha_g_vals[s];
        printf("--- alpha_g = %.4f ---\n", alpha_g);

        Grid *g = grid_alloc(N_grid, L);
        init_braid(g, BIMODAL, -1);
        compute_forces_metric(g, BIMODAL, mass2, alpha_g);

        int total_steps = (int)(T_max / g->dt + 0.5);
        int diag_interval = total_steps / 6;  /* ~6 checkpoints */
        if (diag_interval < 1) diag_interval = 1;

        double wall0 = omp_get_wtime();
        int blown = 0;

        for (int step = 0; step < total_steps; step++) {
            verlet_metric_step(g, BIMODAL, mass2, alpha_g);
            apply_damping_xy(g);

            if ((step+1) % diag_interval == 0 || step == total_steps-1) {
                if (check_blowup(g)) {
                    blown = 1;
                    double T = (step+1) * g->dt;
                    printf("  BLOWUP at T=%.1f (step %d)\n", T, step+1);
                    break;
                }
                double T = (step+1) * g->dt;
                double frac = (double)(step+1) / total_steps;
                double wall_el = omp_get_wtime() - wall0;
                printf("  T=%.1f  (%.0f%%, wall=%.1fs)\n", T, 100*frac, wall_el);
                fflush(stdout);
            }
        }

        /* Diagnostics */
        Result res;
        memset(&res, 0, sizeof(res));
        if (!blown) {
            compute_diagnostics(g, BIMODAL, &res);
            res.stable = !check_blowup(g);
        }

        double l2h = 0.0;
        if (!blown) {
            l2h = compute_l2_metric(g, R_shell);
        }

        res_fc[s]      = blown ? 0.0 : res.fc;
        res_trans[s]    = blown ? 0.0 : res.transverse_l2;
        res_torsion[s]  = blown ? 0.0 : res.torsion_flux;
        res_l2h[s]      = l2h;
        res_stable[s]   = blown ? 0 : res.stable;
        res_l2[s]       = blown ? 0.0 : res.l2_frac;
        res_energy[s]   = blown ? 0.0 : res.energy;

        double wall_el = omp_get_wtime() - wall0;
        printf("  Done: fc=%.4f trans_l2=%.4f torsion=%.4f l2_h=%.6f E=%.1f %s  [%.1fs]\n\n",
               res_fc[s], res_trans[s], res_torsion[s], res_l2h[s], res_energy[s],
               res_stable[s] ? "OK" : "BLOWUP", wall_el);
        fflush(stdout);

        grid_free(g);
    }

    double wall_total = omp_get_wtime() - wall_total_start;

    /* ================================================================
       Summary
       ================================================================ */
    printf("\n=== SUMMARY ===\n\n");
    printf("%10s %8s %10s %10s %10s %8s\n",
           "alpha_g", "fc", "trans_l2", "torsion", "l2_h", "stable");
    printf("---------- -------- ---------- ---------- ---------- --------\n");
    for (int s = 0; s < n_alpha; s++) {
        printf("%10.4f %8.4f %10.4f %10.4f %10.6f %8s\n",
               alpha_g_vals[s], res_fc[s], res_trans[s], res_torsion[s],
               res_l2h[s], res_stable[s] ? "OK" : "BLOWUP");
    }
    printf("\nTotal wall time: %.1f s\n", wall_total);

    /* ================================================================
       Write TSV
       ================================================================ */
    FILE *fp = fopen("data/t5_results.tsv", "w");
    if (fp) {
        fprintf(fp, "alpha_g\tfc\ttrans_l2\ttorsion\tl2_h\tenergy\tstable\n");
        for (int s = 0; s < n_alpha; s++) {
            fprintf(fp, "%.4f\t%.6f\t%.6f\t%.6f\t%.6f\t%.2f\t%d\n",
                    alpha_g_vals[s], res_fc[s], res_trans[s], res_torsion[s],
                    res_l2h[s], res_energy[s], res_stable[s]);
        }
        fclose(fp);
        printf("\nWrote data/t5_results.tsv\n");
    }

    printf("\n=== T5 Complete ===\n");
    return 0;
}
