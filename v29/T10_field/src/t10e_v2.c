/*  t10e_v2.c — Field Response (differential method)
 *
 *  Run TWO copies of the braid simulation:
 *  - Copy A: braid evolves normally (control)
 *  - Copy B: braid + impulse at t=0
 *  The response = B - A at each point and time.
 *
 *  This cleanly isolates the impulse response from the braid dynamics.
 *
 *  Build: gcc -O3 -fopenmp -o t10e_v2 src/t10e_v2.c -lm
 */

#include "../../src/braid_core.h"
#include <float.h>

int main(int argc, char **argv) {
    int N = 96;
    double L = 20.0;
    double T_equil = 200.0;
    double T_resp = 30.0;
    double impulse_amp = 0.01; /* small enough for linearity */

    printf("=== T10E v2: Differential Response Function ===\n\n");
    printf("Grid: N=%d, L=%.1f, dx=%.4f\n", N, L, 2.0*L/(N-1));
    printf("Impulse: %.4f (velocity kick)\n", impulse_amp);
    printf("Response time: T=%.0f\n\n", T_resp);

    bimodal_init_params();
    double dx = 2.0*L/(N-1);
    double dt_grid = 0.20 * dx;
    int resp_steps = (int)(T_resp / dt_grid);
    int N3 = N*N*N;
    int NN = N*N;

    /* Equilibrate */
    printf("Equilibrating braid...\n");
    fflush(stdout);
    Grid *g_ref = grid_alloc(N, L);
    init_braid(g_ref, BIMODAL, -1);
    compute_forces(g_ref, BIMODAL, 2.25);
    int equil_steps = (int)(T_equil / g_ref->dt);
    for (int s = 0; s < equil_steps; s++) {
        verlet_full_step(g_ref, BIMODAL, 2.25);
        apply_damping_xy(g_ref);
    }
    Result res;
    compute_diagnostics(g_ref, BIMODAL, &res);
    printf("  fc=%.4f, E=%.1f\n\n", res.fc, res.energy);

    /* Create copy for impulse run */
    Grid *g_imp = grid_alloc(N, L);
    for (int a = 0; a < NFIELDS; a++) {
        memcpy(g_imp->phi[a], g_ref->phi[a], N3*sizeof(double));
        memcpy(g_imp->vel[a], g_ref->vel[a], N3*sizeof(double));
        memcpy(g_imp->acc[a], g_ref->acc[a], N3*sizeof(double));
    }

    FILE *fout = fopen("data/t10e_v2_response.tsv", "w");
    fprintf(fout, "location\tcomp_exc\tt\tdir\tr\tdelta_0\tdelta_1\tdelta_2\n");

    /* Test impulse at center and off-center */
    double imp_locs[] = {0.0, 5.0, 10.0};
    int n_locs = 3;
    int imp_comps[] = {0, 1}; /* test exciting φ₀ and φ₁ */
    int n_comps = 2;

    for (int il = 0; il < n_locs; il++) {
        for (int ic = 0; ic < n_comps; ic++) {
            double x_imp = imp_locs[il];
            int comp_exc = imp_comps[ic];

            printf("=== Impulse at x=%.0f, exciting φ_%d ===\n", x_imp, comp_exc);
            fflush(stdout);

            /* Reset both grids to equilibrated state */
            /* Re-equilibrate from scratch for clean state */
            /* Actually, just re-copy from the stored equilibrated state */
            /* We need to store the equilibrated state... */

            /* For the first iteration, g_ref is at T_equil.
             * For subsequent ones, it has evolved further.
             * Better: store the equilibrated state separately. */
            /* Quick fix: re-init and re-equilibrate each time.
             * This is slow but correct. For speed, just do first loc/comp. */

            if (il > 0 || ic > 0) {
                /* Re-equilibrate */
                init_braid(g_ref, BIMODAL, -1);
                compute_forces(g_ref, BIMODAL, 2.25);
                for (int s = 0; s < equil_steps; s++) {
                    verlet_full_step(g_ref, BIMODAL, 2.25);
                    apply_damping_xy(g_ref);
                }
            }

            /* Copy to impulse grid */
            for (int a = 0; a < NFIELDS; a++) {
                memcpy(g_imp->phi[a], g_ref->phi[a], N3*sizeof(double));
                memcpy(g_imp->vel[a], g_ref->vel[a], N3*sizeof(double));
                memcpy(g_imp->acc[a], g_ref->acc[a], N3*sizeof(double));
            }

            /* Apply impulse */
            int i_imp = (int)((x_imp + L) / dx + 0.5);
            int jc = N/2, kc = N/2;
            int idx_imp = i_imp * NN + jc * N + kc;
            g_imp->vel[comp_exc][idx_imp] += impulse_amp / dt_grid;

            /* Recompute forces for the impulse grid (force depends on phi, not vel,
             * so forces are actually the same. But to be safe:) */
            compute_forces(g_imp, BIMODAL, 2.25);

            int meas_every = resp_steps / 100;
            if (meas_every < 1) meas_every = 1;

            double max_delta[3] = {0,0,0};
            double arrival_time_5 = -1, arrival_time_10 = -1;
            double threshold = 1e-8;

            /* Time arrays for speed measurement */
            double speed_x_5 = -1, speed_x_10 = -1;
            double speed_z_5 = -1;

            for (int s = 0; s < resp_steps; s++) {
                verlet_full_step(g_ref, BIMODAL, 2.25);
                apply_damping_xy(g_ref);
                verlet_full_step(g_imp, BIMODAL, 2.25);
                apply_damping_xy(g_imp);

                double t = (s+1) * dt_grid;

                /* Check arrival at r=5 and r=10 along x */
                {
                    int i5 = i_imp + (int)(5.0/dx);
                    if (i5 < N-1 && arrival_time_5 < 0) {
                        int idx5 = i5*NN + jc*N + kc;
                        double d = fabs(g_imp->phi[comp_exc][idx5] - g_ref->phi[comp_exc][idx5]);
                        if (d > threshold) {
                            arrival_time_5 = t;
                            speed_x_5 = 5.0 / t;
                        }
                    }
                    int i10 = i_imp + (int)(10.0/dx);
                    if (i10 < N-1 && arrival_time_10 < 0) {
                        int idx10 = i10*NN + jc*N + kc;
                        double d = fabs(g_imp->phi[comp_exc][idx10] - g_ref->phi[comp_exc][idx10]);
                        if (d > threshold) {
                            arrival_time_10 = t;
                            speed_x_10 = 10.0 / t;
                        }
                    }
                }

                /* Check arrival at r=5 along z */
                if (speed_z_5 < 0) {
                    int kz5 = kc + (int)(5.0/dx);
                    if (kz5 < N) {
                        int idxz5 = i_imp*NN + jc*N + (kz5 % N);
                        double d = fabs(g_imp->phi[comp_exc][idxz5] - g_ref->phi[comp_exc][idxz5]);
                        if (d > threshold) {
                            speed_z_5 = 5.0 / t;
                        }
                    }
                }

                /* Record profiles */
                if (s % meas_every == 0) {
                    /* Along x from impulse */
                    for (int i = i_imp; i < N-1; i++) {
                        double r = (i - i_imp) * dx;
                        if (r > 15.0) break;
                        int idx = i*NN + jc*N + kc;
                        double d0 = g_imp->phi[0][idx] - g_ref->phi[0][idx];
                        double d1 = g_imp->phi[1][idx] - g_ref->phi[1][idx];
                        double d2 = g_imp->phi[2][idx] - g_ref->phi[2][idx];
                        fprintf(fout, "%.0f\t%d\t%.4f\tx\t%.4f\t%.6e\t%.6e\t%.6e\n",
                                x_imp, comp_exc, t, r, d0, d1, d2);
                        for (int a = 0; a < 3; a++) {
                            double da = (a==0)?fabs(d0):(a==1)?fabs(d1):fabs(d2);
                            if (da > max_delta[a]) max_delta[a] = da;
                        }
                    }

                    /* Along z from impulse (at same x) */
                    for (int k = kc; k < N; k++) {
                        double r = (k - kc) * dx;
                        if (r > 15.0) break;
                        int idx = i_imp*NN + jc*N + (k % N);
                        double d0 = g_imp->phi[0][idx] - g_ref->phi[0][idx];
                        double d1 = g_imp->phi[1][idx] - g_ref->phi[1][idx];
                        double d2 = g_imp->phi[2][idx] - g_ref->phi[2][idx];
                        fprintf(fout, "%.0f\t%d\t%.4f\tz\t%.4f\t%.6e\t%.6e\t%.6e\n",
                                x_imp, comp_exc, t, r, d0, d1, d2);
                    }
                }
            }

            printf("  Max deviations: δφ₀=%.4e, δφ₁=%.4e, δφ₂=%.4e\n",
                   max_delta[0], max_delta[1], max_delta[2]);
            printf("  Channel mixing: max(other)/excited = %.4f\n",
                   (comp_exc==0) ?
                   fmax(max_delta[1],max_delta[2])/(max_delta[0]+1e-30) :
                   fmax(max_delta[0],max_delta[2])/(max_delta[1]+1e-30));
            printf("  Signal speed (x, to r=5): %.4f", speed_x_5);
            if (speed_x_5 > 0) printf(" (t_arr=%.2f)", 5.0/speed_x_5);
            printf("\n");
            printf("  Signal speed (x, to r=10): %.4f\n", speed_x_10 > 0 ? speed_x_10 : 0);
            printf("  Signal speed (z, to r=5): %.4f\n", speed_z_5 > 0 ? speed_z_5 : 0);
            printf("  Anisotropy (vz/vx): %.4f\n",
                   (speed_x_5 > 0 && speed_z_5 > 0) ? speed_z_5/speed_x_5 : 0);
            printf("\n");
        }
    }

    fclose(fout);
    grid_free(g_ref);
    grid_free(g_imp);

    printf("=== T10E v2 Complete ===\n");
    printf("Data in data/t10e_v2_response.tsv\n");
    return 0;
}
