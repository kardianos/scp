/*  t10e.c — Field Response Function (Green's Function)
 *
 *  Apply a sharp impulse to the field at a point and measure the response
 *  at all later times. This gives the retarded Green's function G_R(r,t).
 *
 *  Test in two backgrounds:
 *  (a) Vacuum (φ=0): should give standard massive Klein-Gordon Green's function
 *  (b) Braid background: the braid modifies the Green's function
 *
 *  From G_R, extract:
 *  - Propagation speed (light cone)
 *  - Effective mass (oscillation frequency at late times)
 *  - Tensor structure (does the response differ by direction?)
 *  - Decay law (1/r? 1/r²? exponential?)
 *  - Channel mixing (does impulse in φ₀ excite φ₁, φ₂?)
 *
 *  Build: gcc -O3 -fopenmp -o t10e src/t10e.c -lm
 */

#include "../../src/braid_core.h"
#include <float.h>

/* ================================================================
   Measure response along a line from the impulse point
   ================================================================ */

typedef struct {
    double t;
    double r[40];       /* distances */
    double resp[3][40]; /* response for each component */
    int nr;
} Snapshot;

int main(int argc, char **argv) {
    int N = 96;
    double L = 20.0;
    double T_equil = 200.0;
    double T_resp = 40.0;  /* response measurement time */
    double impulse_amp = 0.1; /* impulse amplitude */

    printf("=== T10E: Field Response Function ===\n\n");
    printf("Grid: N=%d, L=%.1f, dx=%.4f\n", N, L, 2.0*L/(N-1));
    printf("Impulse amplitude: %.3f\n", impulse_amp);
    printf("Response time: T=%.0f\n\n", T_resp);

    bimodal_init_params();

    double dx = 2.0*L/(N-1);
    double dt = 0.20 * dx;
    int resp_steps = (int)(T_resp / dt);
    int meas_every = resp_steps / 200;
    if (meas_every < 1) meas_every = 1;

    FILE *fvac = fopen("data/t10e_vacuum.tsv", "w");
    FILE *fbr  = fopen("data/t10e_braid.tsv", "w");
    fprintf(fvac, "comp_excited\tt\tdirection\tr\tresp_0\tresp_1\tresp_2\n");
    fprintf(fbr,  "comp_excited\tt\tdirection\tr\tresp_0\tresp_1\tresp_2\n");

    /* ============================================================
       Part A: Vacuum Green's function (validation)
       ============================================================ */
    printf("=== Part A: Vacuum response ===\n");
    {
        Grid *g = grid_alloc(N, L);
        grid_zero(g);

        /* Apply impulse at center: δφ₀(0,0,0) = ε */
        int ic = N/2, jc = N/2, kc = N/2;
        int idx_c = ic*N*N + jc*N + kc;
        g->vel[0][idx_c] = impulse_amp / dt; /* impulse = Δvel */
        compute_forces(g, BIMODAL, 2.25);

        printf("  Exciting φ₀ at center...\n");

        /* Measure response along x-axis, y-axis, z-axis */
        for (int s = 0; s < resp_steps; s++) {
            verlet_full_step(g, BIMODAL, 2.25);
            /* No damping — want to see the full response */

            if (s % meas_every == 0) {
                double t = s * dt;

                /* Along x-axis (j=jc, k=kc) */
                for (int i = ic; i < N-1; i++) {
                    double r = (i - ic) * dx;
                    if (r > L * 0.8) break;
                    int idx = i*N*N + jc*N + kc;
                    fprintf(fvac, "0\t%.4f\tx\t%.4f\t%.6e\t%.6e\t%.6e\n",
                            t, r, g->phi[0][idx], g->phi[1][idx], g->phi[2][idx]);
                }

                /* Along z-axis (i=ic, j=jc) */
                for (int k = kc; k < N; k++) {
                    double r = (k - kc) * dx;
                    if (r > L * 0.8) break;
                    int idx = ic*N*N + jc*N + k;
                    fprintf(fvac, "0\t%.4f\tz\t%.4f\t%.6e\t%.6e\t%.6e\n",
                            t, r, g->phi[0][idx], g->phi[1][idx], g->phi[2][idx]);
                }
            }
        }

        /* Final snapshot: radial decay */
        printf("  Final snapshot (t=%.1f), radial decay along x:\n", T_resp);
        printf("    r\tφ₀\tφ₁\tφ₂\n");
        for (int i = ic; i < N-1; i += 4) {
            double r = (i - ic) * dx;
            if (r > L * 0.7) break;
            int idx = i*N*N + jc*N + kc;
            printf("    %.2f\t%.4e\t%.4e\t%.4e\n",
                   r, g->phi[0][idx], g->phi[1][idx], g->phi[2][idx]);
        }

        /* Measure: did the impulse in φ₀ excite φ₁ or φ₂? */
        double max_0 = 0, max_1 = 0, max_2 = 0;
        for (int i = 0; i < N*N*N; i++) {
            if (fabs(g->phi[0][i]) > max_0) max_0 = fabs(g->phi[0][i]);
            if (fabs(g->phi[1][i]) > max_1) max_1 = fabs(g->phi[1][i]);
            if (fabs(g->phi[2][i]) > max_2) max_2 = fabs(g->phi[2][i]);
        }
        printf("  Max amplitudes: φ₀=%.4e, φ₁=%.4e, φ₂=%.4e\n", max_0, max_1, max_2);
        printf("  Channel mixing (φ₁/φ₀): %.4e\n", max_1 / (max_0 + 1e-30));
        printf("  → Should be ~0 in vacuum (no mixing)\n\n");

        grid_free(g);
    }

    /* ============================================================
       Part B: Braid background response
       ============================================================ */
    printf("=== Part B: Braid background response ===\n");
    {
        /* Equilibrate braid */
        printf("  Equilibrating braid...\n");
        fflush(stdout);
        Grid *g = grid_alloc(N, L);
        init_braid(g, BIMODAL, -1);
        compute_forces(g, BIMODAL, 2.25);
        int equil_steps = (int)(T_equil / g->dt);
        for (int s = 0; s < equil_steps; s++) {
            verlet_full_step(g, BIMODAL, 2.25);
            apply_damping_xy(g);
        }
        Result res;
        compute_diagnostics(g, BIMODAL, &res);
        printf("  fc=%.4f\n", res.fc);

        /* Save background state */
        int N3 = N*N*N;
        double *bg[NFIELDS];
        for (int a = 0; a < NFIELDS; a++) {
            bg[a] = malloc(N3 * sizeof(double));
            memcpy(bg[a], g->phi[a], N3 * sizeof(double));
        }

        /* Test impulse at different locations:
         * (a) At braid center
         * (b) At offset r=5 from center
         * (c) At offset r=10 from center (near vacuum)
         */

        double offsets[] = {0, 5, 10};
        int noff = 3;

        for (int io = 0; io < noff; io++) {
            double x_imp = offsets[io];
            int i_imp = (int)((x_imp + L) / dx + 0.5);
            if (i_imp >= N-1) i_imp = N-2;
            int jc = N/2, kc = N/2;
            int idx_imp = i_imp * N*N + jc * N + kc;

            printf("\n  Impulse at x=%.0f (r=%.0f from center):\n", x_imp, x_imp);

            /* Reset to background + impulse */
            for (int a = 0; a < NFIELDS; a++) {
                memcpy(g->phi[a], bg[a], N3 * sizeof(double));
                memset(g->vel[a], 0, N3 * sizeof(double));
            }
            g->vel[0][idx_imp] += impulse_amp / dt;
            compute_forces(g, BIMODAL, 2.25);

            /* Track response */
            double max_mix = 0;
            for (int s = 0; s < resp_steps; s++) {
                verlet_full_step(g, BIMODAL, 2.25);

                if (s % meas_every == 0) {
                    double t = s * dt;

                    /* Along x-axis from impulse point */
                    for (int i = i_imp; i < N-1 && i < i_imp + 40; i++) {
                        double r = (i - i_imp) * dx;
                        int idx = i*N*N + jc*N + kc;
                        /* Subtract background */
                        double d0 = g->phi[0][idx] - bg[0][idx];
                        double d1 = g->phi[1][idx] - bg[1][idx];
                        double d2 = g->phi[2][idx] - bg[2][idx];
                        fprintf(fbr, "0_at%.0f\t%.4f\tx\t%.4f\t%.6e\t%.6e\t%.6e\n",
                                x_imp, t, r, d0, d1, d2);
                    }

                    /* Check mixing at this timestep */
                    for (int idx = 0; idx < N3; idx++) {
                        double d1 = fabs(g->phi[1][idx] - bg[1][idx]);
                        double d2 = fabs(g->phi[2][idx] - bg[2][idx]);
                        double mx = (d1 > d2) ? d1 : d2;
                        if (mx > max_mix) max_mix = mx;
                    }
                }
            }

            /* Final: measure amplitude of each component's deviation from bg */
            double max_d[3] = {0,0,0};
            for (int idx = 0; idx < N3; idx++) {
                for (int a = 0; a < 3; a++) {
                    double d = fabs(g->phi[a][idx] - bg[a][idx]);
                    if (d > max_d[a]) max_d[a] = d;
                }
            }
            printf("    Max deviations: δφ₀=%.4e, δφ₁=%.4e, δφ₂=%.4e\n",
                   max_d[0], max_d[1], max_d[2]);
            printf("    Channel mixing (max(δφ₁,δφ₂)/δφ₀): %.4f\n",
                   fmax(max_d[1], max_d[2]) / (max_d[0] + 1e-30));

            /* Measure effective speed: when does response reach r=5? */
            {
                int i_detect = i_imp + (int)(5.0/dx);
                if (i_detect < N-1) {
                    double first_arrival = -1;
                    double threshold = 1e-6;

                    /* Re-run to find arrival time */
                    for (int a = 0; a < NFIELDS; a++) {
                        memcpy(g->phi[a], bg[a], N3 * sizeof(double));
                        memset(g->vel[a], 0, N3 * sizeof(double));
                    }
                    g->vel[0][idx_imp] += impulse_amp / dt;
                    compute_forces(g, BIMODAL, 2.25);

                    for (int s = 0; s < resp_steps; s++) {
                        verlet_full_step(g, BIMODAL, 2.25);
                        int idx_det = i_detect*N*N + jc*N + kc;
                        double d0 = fabs(g->phi[0][idx_det] - bg[0][idx_det]);
                        if (d0 > threshold && first_arrival < 0) {
                            first_arrival = s * dt;
                            break;
                        }
                    }
                    if (first_arrival > 0) {
                        double v_eff = 5.0 / first_arrival;
                        printf("    Signal speed to r=5: v=%.4f (arrival at t=%.2f)\n",
                               v_eff, first_arrival);
                        printf("    Light speed = 1.0, massive vg(k→inf)=1.0\n");
                    }
                }
            }
        }

        for (int a = 0; a < NFIELDS; a++) free(bg[a]);
        grid_free(g);
    }

    fclose(fvac);
    fclose(fbr);

    printf("\n=== SUMMARY ===\n");
    printf("Key findings:\n");
    printf("  1. Does the braid background cause channel mixing?\n");
    printf("     (impulse in φ₀ excites φ₁,φ₂?) → indicator of mode coupling\n");
    printf("  2. Is the signal speed modified by the background?\n");
    printf("     (slower near braid = gravity, faster = anti-gravity)\n");
    printf("  3. Does the response differ by direction?\n");
    printf("     (anisotropy = tensor structure)\n");

    return 0;
}
