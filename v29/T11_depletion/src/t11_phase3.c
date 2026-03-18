/*  t11_phase3.c — Two braids in background: do they attract via depletion?
 *
 *  Phase 2 showed: the braid DOES reduce background amplitude at its core.
 *  Now: does this create an attractive force between two braids?
 *
 *  Setup:
 *  - Two braids at (x,y) = (-D/2, 0) and (+D/2, 0), separated by D
 *  - Background field: A_bg * cos(k_bg * z) as before
 *  - Measure: does the separation D(t) decrease?
 *  - Compare: same two braids WITHOUT background → is attraction stronger with bg?
 *
 *  The depletion hypothesis predicts: attraction should be STRONGER with background
 *  because the bg creates an additional depletion-mediated force.
 *
 *  Build: gcc -O3 -fopenmp -o t11p3 src/t11_phase3.c -lm
 */

#include "../../src/braid_core.h"

/* Initialize braid centered at (cx, cy) */
static void init_braid_at(Grid *g, const double *phys, double cx, double cy, double m_init_ov) {
    double A[3]     = {phys[0], phys[1], phys[2]};
    double delta[3] = {0.0, phys[3], phys[4]};
    double R_tube   = phys[5];
    double ellip    = phys[6];
    double ell_ang  = phys[7];
    double k_fac    = phys[8];
    double A_bg_p   = phys[9];
    double R_disp   = phys[10];
    double ell_rot  = phys[11];
    double m_init   = (m_init_ov >= 0) ? m_init_ov : phys[14];

    int N = g->N, NN = N * N;
    double dx = g->dx, L = g->L;
    double k = k_fac * PI / L;
    double omega = sqrt(k * k + m_init * m_init);
    double inv_2R2 = 1.0 / (2.0 * R_tube * R_tube);
    double sx = 1.0 + ellip, sy = 1.0 - ellip;

    double strand_cx[3], strand_cy[3], ea[3];
    for (int a = 0; a < 3; a++) {
        double ang = 2.0 * PI * a / 3.0;
        strand_cx[a] = cx + R_disp * cos(ang);
        strand_cy[a] = cy + R_disp * sin(ang);
        ea[a] = (ell_rot > 0.5) ? ell_ang + 2.0*PI*a/3.0 : ell_ang;
    }

    for (int i = 0; i < N; i++) {
        double x = -L + i * dx;
        for (int j = 0; j < N; j++) {
            double y = -L + j * dx;
            for (int kk = 0; kk < N; kk++) {
                double z = -L + kk * dx;
                int idx = i * NN + j * N + kk;
                for (int a = 0; a < NFIELDS; a++) {
                    double xc = x - strand_cx[a], yc = y - strand_cy[a];
                    double ca = cos(ea[a]), sa = sin(ea[a]);
                    double xr = xc*ca + yc*sa;
                    double yr = -xc*sa + yc*ca;
                    double r2e = xr*xr/(sx*sx) + yr*yr/(sy*sy);
                    double env = exp(-r2e * inv_2R2);
                    double ph = k * z + delta[a];
                    double amp = A[a] * env + A_bg_p;
                    g->phi[a][idx] += amp * cos(ph);
                    g->vel[a][idx] += omega * amp * sin(ph);
                }
            }
        }
    }
}

/* Find centroid of phi² in a region */
static double find_centroid_x(Grid *g, double xmin, double xmax) {
    int N = g->N, NN = N*N;
    double dx = g->dx, L = g->L;
    double M = 0, Mx = 0;
    for (int i = 1; i < N-1; i++) {
        double x = -L + i * dx;
        if (x < xmin || x > xmax) continue;
        for (int j = 1; j < N-1; j++) {
            double y = -L + j * dx;
            for (int k = 0; k < N; k++) {
                int idx = i*NN + j*N + k;
                double rho = 0;
                for (int a = 0; a < NFIELDS; a++)
                    rho += g->phi[a][idx] * g->phi[a][idx];
                M += rho;
                Mx += rho * x;
            }
        }
    }
    return Mx / (M + 1e-30);
}

int main(void) {
    bimodal_init_params();
    double mass2 = BIMODAL[14] * BIMODAL[14];
    int N = 192;
    double L = 40.0;
    double D = 20.0;  /* initial separation */
    double T = 100.0;

    printf("=== T11 Phase 3: Two-Braid Depletion Attraction ===\n");
    printf("N=%d L=%.1f D=%.1f T=%.1f\n", N, L, D, T);

    double dx = 2.0 * L / (N - 1);
    double dt = 0.20 * dx;
    int NN = N*N;

    double A_bg_vals[] = {0.0, 0.1, 0.3};
    int nA = 3;

    for (int ia = 0; ia < nA; ia++) {
        double A_bg = A_bg_vals[ia];
        double k_bg = PI / L;
        double omega_bg = sqrt(k_bg * k_bg + mass2);

        printf("\n===== A_bg = %.3f =====\n", A_bg);

        Grid *g = grid_alloc(N, L);
        grid_zero(g);

        /* Place two braids */
        init_braid_at(g, BIMODAL, -D/2, 0.0, -1);
        init_braid_at(g, BIMODAL, +D/2, 0.0, -1);

        /* Add background */
        if (A_bg > 0) {
            for (int i = 0; i < N; i++)
                for (int j = 0; j < N; j++)
                    for (int k = 0; k < N; k++) {
                        int idx = i*NN + j*N + k;
                        double z = -L + k * dx;
                        for (int a = 0; a < NFIELDS; a++) {
                            g->phi[a][idx] += A_bg * cos(k_bg * z);
                            g->vel[a][idx] += A_bg * omega_bg * sin(k_bg * z);
                        }
                    }
        }

        compute_forces(g, BIMODAL, -1);

        int n_steps = (int)(T / dt);
        int snap_every = n_steps / 40;

        char fname[256];
        snprintf(fname, sizeof(fname), "data/phase3_Abg%.3f_separation.dat", A_bg);
        FILE *fp = fopen(fname, "w");
        fprintf(fp, "# t  x1  x2  D  E_total\n");

        printf("  Running %d steps...\n", n_steps);

        for (int step = 0; step <= n_steps; step++) {
            if (step > 0) {
                verlet_full_step(g, BIMODAL, -1);
                apply_damping_xy(g);
            }

            if (step % snap_every == 0) {
                double t = step * dt;
                double x1 = find_centroid_x(g, -L, 0);
                double x2 = find_centroid_x(g, 0, L);
                double sep = x2 - x1;

                Result res;
                snprintf(res.label, sizeof(res.label), "s%d", step);
                compute_diagnostics(g, BIMODAL, &res);

                fprintf(fp, "%.4f  %.4f  %.4f  %.4f  %.2f\n",
                        t, x1, x2, sep, res.energy);

                if (step % (snap_every * 10) == 0) {
                    printf("  t=%6.1f  x1=%+6.2f  x2=%+6.2f  D=%.2f  E=%.0f\n",
                           t, x1, x2, sep, res.energy);
                }

                if (check_blowup(g)) { printf("BLOWUP!\n"); break; }
            }
        }
        fclose(fp);

        grid_free(g);
    }

    printf("\n=== Phase 3 Complete ===\n");
    printf("Compare separation D(t) between A_bg=0 (no bg) and A_bg>0.\n");
    printf("If background makes braids attract MORE: depletion-mediated gravity.\n");
    printf("If same or less: no depletion gravity.\n");

    return 0;
}
