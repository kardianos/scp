/*  t10b.c — Effective Metric Extraction
 *
 *  Test whether a frozen braid background creates an effective metric
 *  for test waves. The background is set up and equilibrated, then frozen.
 *  Small-amplitude test wave packets are launched from different directions.
 *  We measure:
 *    1. Speed change (compared to vacuum propagation)
 *    2. Deflection angle (bending toward/away from braid)
 *    3. Time delay (Shapiro-like delay)
 *    4. Amplitude change (focusing/defocusing)
 *
 *  If waves slow down and bend TOWARD the braid, that is gravity.
 *
 *  Method:
 *  - Phase 1: Evolve braid for T_equil time with absorbing BC to equilibrate
 *  - Phase 2: Freeze the background (copy phi to bg_phi, set bg_vel=0)
 *  - Phase 3: For each test direction:
 *      (a) Initialize small Gaussian wave packet at offset from braid
 *      (b) Evolve ONLY the test perturbation on the frozen background
 *      (c) Track the wave packet centroid to measure speed + deflection
 *      (d) Compare with identical wave in vacuum (no background)
 *
 *  Build: gcc -O3 -fopenmp -o t10b src/t10b.c -lm
 */

#include "../../src/braid_core.h"
#include <float.h>

/* ================================================================
   Test wave on frozen background
   ================================================================ */

typedef struct {
    double *bg[NFIELDS];     /* frozen background fields */
    double *dph[NFIELDS];    /* test perturbation phi */
    double *dvel[NFIELDS];   /* test perturbation velocity */
    double *dacc[NFIELDS];   /* test perturbation acceleration */
    int N;
    double L, dx, dt;
} TestWave;

static TestWave *tw_alloc(int N, double L) {
    TestWave *tw = calloc(1, sizeof(TestWave));
    int N3 = N * N * N;
    for (int a = 0; a < NFIELDS; a++) {
        tw->bg[a]   = calloc(N3, sizeof(double));
        tw->dph[a]  = calloc(N3, sizeof(double));
        tw->dvel[a] = calloc(N3, sizeof(double));
        tw->dacc[a] = calloc(N3, sizeof(double));
    }
    tw->N = N; tw->L = L;
    tw->dx = 2.0 * L / (N - 1);
    tw->dt = 0.20 * tw->dx;
    return tw;
}

static void tw_free(TestWave *tw) {
    for (int a = 0; a < NFIELDS; a++) {
        free(tw->bg[a]); free(tw->dph[a]);
        free(tw->dvel[a]); free(tw->dacc[a]);
    }
    free(tw);
}

/* Copy braid background from Grid */
static void tw_set_background(TestWave *tw, Grid *g) {
    int N3 = g->N * g->N * g->N;
    for (int a = 0; a < NFIELDS; a++)
        memcpy(tw->bg[a], g->phi[a], N3 * sizeof(double));
}

/* Initialize Gaussian wave packet in perturbation field
 * Component 'comp', centered at (x0,y0,z0), wave vector (kx,ky,kz), width sigma */
static void tw_init_packet(TestWave *tw, int comp,
                           double x0, double y0, double z0,
                           double kx, double ky, double kz,
                           double sigma, double amplitude) {
    int N = tw->N, NN = N * N;
    double dx = tw->dx, L = tw->L;
    double inv_2s2 = 1.0 / (2.0 * sigma * sigma);
    double k2 = kx*kx + ky*ky + kz*kz;
    double mass2 = 2.25; /* m=1.5 */
    double omega = sqrt(k2 + mass2);

    for (int a = 0; a < NFIELDS; a++) {
        int N3 = N * N * N;
        memset(tw->dph[a], 0, N3 * sizeof(double));
        memset(tw->dvel[a], 0, N3 * sizeof(double));
        memset(tw->dacc[a], 0, N3 * sizeof(double));
    }

    for (int i = 0; i < N; i++) {
        double x = -L + i * dx;
        for (int j = 0; j < N; j++) {
            double y = -L + j * dx;
            for (int kk = 0; kk < N; kk++) {
                double z = -L + kk * dx;
                int idx = i * NN + j * N + kk;
                double rx = x - x0, ry = y - y0, rz = z - z0;
                double r2 = rx*rx + ry*ry + rz*rz;
                double env = amplitude * exp(-r2 * inv_2s2);
                double phase = kx * x + ky * y + kz * z;
                tw->dph[comp][idx] = env * cos(phase);
                tw->dvel[comp][idx] = -omega * env * cos(phase);
                /* Note: vel = -omega*cos for a packet moving in +k direction */
            }
        }
    }
}

/* Compute linearized force on the perturbation around the frozen background.
 *
 * Full EOM: d²φ/dt² = ∇²φ - m²φ - (μP/denom²) * ∂P/∂φ_a
 * where P = φ₀φ₁φ₂, denom = 1 + κP²
 *
 * Linearize: φ_a = bg_a + δφ_a
 * Keep terms linear in δφ. The key is the derivative of the triple-product force.
 */
static void tw_compute_forces(TestWave *tw, const double *phys) {
    int N = tw->N, NN = N*N, N3 = N*N*N;
    double idx2 = 1.0 / (tw->dx * tw->dx);
    double mu    = phys[12];
    double kappa = phys[13];
    double mass2 = 2.25;

    #pragma omp parallel for schedule(static)
    for (int idx = 0; idx < N3; idx++) {
        int i = idx / NN, j = (idx / N) % N, k = idx % N;
        if (i < 1 || i >= N-1 || j < 1 || j >= N-1) {
            tw->dacc[0][idx] = tw->dacc[1][idx] = tw->dacc[2][idx] = 0;
            continue;
        }
        int kp = (k+1)%N, km = (k-1+N)%N;

        /* Background at this point */
        double b0 = tw->bg[0][idx], b1 = tw->bg[1][idx], b2 = tw->bg[2][idx];
        double P_bg = b0 * b1 * b2;
        double denom_bg = 1.0 + kappa * P_bg * P_bg;
        double denom2 = denom_bg * denom_bg;

        /* Perturbations */
        double d0 = tw->dph[0][idx], d1 = tw->dph[1][idx], d2 = tw->dph[2][idx];

        for (int a = 0; a < NFIELDS; a++) {
            /* Laplacian of perturbation */
            int idx_kp = i*NN + j*N + kp;
            int idx_km = i*NN + j*N + km;
            double lap = (tw->dph[a][idx+NN] + tw->dph[a][idx-NN]
                        + tw->dph[a][idx+N]  + tw->dph[a][idx-N]
                        + tw->dph[a][idx_kp] + tw->dph[a][idx_km]
                        - 6.0 * tw->dph[a][idx]) * idx2;

            /* Linearized triple-product force:
             * F_a = -d/dφ_a [ (μ/2)P²/(1+κP²) ]
             *      = -μ P / (1+κP²)² * ∂P/∂φ_a
             *
             * Linearizing F_a around bg: dF_a/dφ_b evaluated at bg
             *
             * For the full nonlinear potential V = (μ/2)P²/(1+κP²):
             * ∂²V/∂φ_a∂φ_b = μ * [∂P/∂φ_a * ∂P/∂φ_b * (1-3κP²) / (1+κP²)³
             *                      + P/(1+κP²)² * ∂²P/∂φ_a∂φ_b]
             *
             * where ∂P/∂φ_0 = φ_1 φ_2, etc.
             * ∂²P/∂φ_0∂φ_1 = φ_2, ∂²P/∂φ_0∂φ_0 = 0, etc.
             */
            double dP_da[3];
            dP_da[0] = b1 * b2;
            dP_da[1] = b0 * b2;
            dP_da[2] = b0 * b1;

            double factor1 = mu * (1.0 - 3.0*kappa*P_bg*P_bg) / (denom2 * denom_bg);
            double factor2 = mu * P_bg / denom2;

            double f_lin = 0;
            for (int b = 0; b < NFIELDS; b++) {
                double db = (b==0) ? d0 : (b==1) ? d1 : d2;
                /* Term 1: dP/da * dP/db * factor1 * db */
                f_lin += dP_da[a] * dP_da[b] * factor1 * db;
                /* Term 2: P/(1+kP²)² * d²P/(da db) * db */
                if (a != b) {
                    /* d²P/(da db) = the third field */
                    int c = 3 - a - b; /* 0+1+2=3 */
                    double d2P = (c==0) ? b0 : (c==1) ? b1 : b2;
                    f_lin += factor2 * d2P * db;
                }
                /* d²P/da² = 0 for triple product */
            }

            tw->dacc[a][idx] = lap - mass2 * tw->dph[a][idx] - f_lin;
        }
    }
}

/* Verlet step for test wave */
static void tw_verlet_step(TestWave *tw, const double *phys) {
    int N3 = tw->N * tw->N * tw->N;
    double hdt = 0.5 * tw->dt;

    /* kick */
    for (int a = 0; a < NFIELDS; a++)
        for (int idx = 0; idx < N3; idx++)
            tw->dvel[a][idx] += hdt * tw->dacc[a][idx];
    /* drift */
    for (int a = 0; a < NFIELDS; a++)
        for (int idx = 0; idx < N3; idx++)
            tw->dph[a][idx] += tw->dt * tw->dvel[a][idx];
    /* force */
    tw_compute_forces(tw, phys);
    /* kick */
    for (int a = 0; a < NFIELDS; a++)
        for (int idx = 0; idx < N3; idx++)
            tw->dvel[a][idx] += hdt * tw->dacc[a][idx];
}

/* Absorbing damping in xy (same as braid_core) */
static void tw_apply_damping_xy(TestWave *tw) {
    int N = tw->N, NN = N * N;
    double dx = tw->dx, L = tw->L;
    double r_start = 0.70 * L, r_end = 0.95 * L;
    double inv_dr = 1.0 / (r_end - r_start + 1e-30);

    for (int i = 0; i < N; i++) {
        double x = -L + i * dx;
        for (int j = 0; j < N; j++) {
            double y = -L + j * dx;
            double rp = sqrt(x*x + y*y);
            if (rp <= r_start) continue;
            double f = (rp - r_start) * inv_dr;
            if (f > 1.0) f = 1.0;
            double damp = 1.0 - 0.98 * f * f;
            for (int kk = 0; kk < N; kk++) {
                int idx = i * NN + j * N + kk;
                for (int a = 0; a < NFIELDS; a++) {
                    tw->dph[a][idx] *= damp;
                    tw->dvel[a][idx] *= damp;
                }
            }
        }
    }
}

/* Track wave packet centroid (energy-weighted) for a specific component */
static void tw_centroid(TestWave *tw, int comp, double *cx, double *cy, double *cz, double *total_E) {
    int N = tw->N, NN = N*N;
    double dx = tw->dx, L = tw->L;
    double wx = 0, wy = 0, wz = 0, wtot = 0;

    for (int i = 1; i < N-1; i++) {
        double x = -L + i * dx;
        for (int j = 1; j < N-1; j++) {
            double y = -L + j * dx;
            for (int kk = 0; kk < N; kk++) {
                double z = -L + kk * dx;
                int idx = i * NN + j * N + kk;
                /* energy density of perturbation */
                double v2 = tw->dvel[comp][idx] * tw->dvel[comp][idx];
                double p2 = tw->dph[comp][idx] * tw->dph[comp][idx];
                /* Simple proxy: kinetic + field^2 */
                double e = v2 + p2;
                wx += e * x;
                wy += e * y;
                wz += e * z;
                wtot += e;
            }
        }
    }
    *cx = wx / (wtot + 1e-30);
    *cy = wy / (wtot + 1e-30);
    *cz = wz / (wtot + 1e-30);
    *total_E = wtot;
}

/* Track centroid using ALL components */
static void tw_centroid_all(TestWave *tw, double *cx, double *cy, double *cz, double *total_E) {
    int N = tw->N, NN = N*N;
    double dx = tw->dx, L = tw->L;
    double wx = 0, wy = 0, wz = 0, wtot = 0;

    for (int i = 1; i < N-1; i++) {
        double x = -L + i * dx;
        for (int j = 1; j < N-1; j++) {
            double y = -L + j * dx;
            for (int kk = 0; kk < N; kk++) {
                double z = -L + kk * dx;
                int idx = i * NN + j * N + kk;
                double e = 0;
                for (int a = 0; a < NFIELDS; a++) {
                    e += tw->dvel[a][idx] * tw->dvel[a][idx];
                    e += tw->dph[a][idx] * tw->dph[a][idx];
                }
                wx += e * x;
                wy += e * y;
                wz += e * z;
                wtot += e;
            }
        }
    }
    *cx = wx / (wtot + 1e-30);
    *cy = wy / (wtot + 1e-30);
    *cz = wz / (wtot + 1e-30);
    *total_E = wtot;
}

/* ================================================================
   Main experiment
   ================================================================ */

int main(int argc, char **argv) {
    /* Parameters */
    int N = 96;
    double L = 20.0;
    double T_equil = 200.0;   /* equilibration time */
    double T_prop = 60.0;     /* propagation time for test wave */
    double amp = 0.01;        /* test wave amplitude (small!) */
    double sigma = 2.0;       /* wave packet width */
    double k_test = 2.0;      /* test wave number */

    /* Impact parameters to test */
    double b_values[] = {0.0, 3.0, 6.0, 9.0, 12.0, 15.0};
    int n_b = 6;

    printf("=== T10B: Effective Metric Extraction ===\n\n");
    printf("Grid: N=%d, L=%.1f, dx=%.4f\n", N, L, 2.0*L/(N-1));
    printf("T_equil=%.0f, T_prop=%.0f\n", T_equil, T_prop);
    printf("Test wave: amp=%.3f, sigma=%.1f, k=%.1f\n", amp, sigma, k_test);
    printf("Mass dispersion: omega = sqrt(k^2 + m^2) = %.4f\n",
           sqrt(k_test*k_test + 2.25));
    printf("Group velocity: v_g = k/omega = %.4f\n",
           k_test / sqrt(k_test*k_test + 2.25));
    printf("\n");

    /* Phase 1: Create and equilibrate the braid */
    printf("Phase 1: Equilibrating braid (T=%.0f)...\n", T_equil);
    fflush(stdout);

    bimodal_init_params();
    Grid *g = grid_alloc(N, L);
    init_braid(g, BIMODAL, -1);
    compute_forces(g, BIMODAL, 2.25);

    double dx_grid = g->dx, dt = g->dt;
    int equil_steps = (int)(T_equil / dt);
    for (int s = 0; s < equil_steps; s++) {
        verlet_full_step(g, BIMODAL, 2.25);
        apply_damping_xy(g);
    }

    Result res;
    snprintf(res.label, sizeof(res.label), "equilibrated");
    compute_diagnostics(g, BIMODAL, &res);
    printf("  After equilibration: ");
    print_result(&res);
    printf("\n");

    /* Phase 2: Freeze background */
    printf("Phase 2: Freezing background...\n");
    TestWave *tw = tw_alloc(N, L);
    tw_set_background(tw, g);
    printf("  Background frozen. Braid fc=%.4f\n\n", res.fc);

    /* Open data file */
    FILE *fdata = fopen("data/t10b_trajectories.tsv", "w");
    if (!fdata) { fdata = fopen("T10_field/data/t10b_trajectories.tsv", "w"); }
    if (!fdata) {
        printf("ERROR: cannot open data file\n");
        return 1;
    }
    fprintf(fdata, "mode\tb\tt\tcx\tcy\tcz\tenergy\n");

    /* Summary data file */
    FILE *fsum = fopen("data/t10b_summary.tsv", "w");
    if (!fsum) { fsum = fopen("T10_field/data/t10b_summary.tsv", "w"); }
    if (!fsum) {
        printf("ERROR: cannot open summary file\n");
        return 1;
    }
    fprintf(fsum, "b\tvx_vac\tvy_vac\tvx_braid\tvy_braid\tdeflection_rad\tspeed_ratio\ttime_delay\n");

    double vg_theory = k_test / sqrt(k_test * k_test + 2.25);
    int prop_steps = (int)(T_prop / dt);
    int output_every = prop_steps / 100;
    if (output_every < 1) output_every = 1;

    /* Phase 3: For each impact parameter, run vacuum + braid */
    for (int ib = 0; ib < n_b; ib++) {
        double b = b_values[ib];
        printf("=== Impact parameter b = %.1f ===\n", b);

        /* Wave packet starts at x=-12, moving in +x direction
         * at transverse offset y = b */
        double x0 = -12.0, y0 = b, z0 = 0.0;
        double kx = k_test, ky = 0.0, kz = 0.0;

        /* --- Vacuum run --- */
        printf("  Vacuum run...\n");
        fflush(stdout);

        /* Zero background for vacuum */
        TestWave *tw_vac = tw_alloc(N, L);
        /* bg stays zero = vacuum */

        tw_init_packet(tw_vac, 0, x0, y0, z0, kx, ky, kz, sigma, amp);
        tw_compute_forces(tw_vac, BIMODAL);

        double vac_cx0, vac_cy0, vac_cz0, vac_E0;
        tw_centroid_all(tw_vac, &vac_cx0, &vac_cy0, &vac_cz0, &vac_E0);

        double vac_cx_final = vac_cx0, vac_cy_final = vac_cy0;
        double vac_cx_prev = vac_cx0, vac_cy_prev = vac_cy0;

        for (int s = 0; s < prop_steps; s++) {
            tw_verlet_step(tw_vac, BIMODAL);
            tw_apply_damping_xy(tw_vac);
            if (s % output_every == 0) {
                double cx, cy, cz, E;
                tw_centroid_all(tw_vac, &cx, &cy, &cz, &E);
                fprintf(fdata, "vacuum\t%.1f\t%.3f\t%.6f\t%.6f\t%.6f\t%.6e\n",
                        b, s*dt, cx, cy, cz, E);
                vac_cx_prev = vac_cx_final;
                vac_cy_prev = vac_cy_final;
                vac_cx_final = cx;
                vac_cy_final = cy;
            }
        }

        double vac_vx = (vac_cx_final - vac_cx0) / T_prop;
        double vac_vy = (vac_cy_final - vac_cy0) / T_prop;
        printf("    Vacuum: cx %.4f -> %.4f, vx=%.4f, vy=%.6f\n",
               vac_cx0, vac_cx_final, vac_vx, vac_vy);

        tw_free(tw_vac);

        /* --- Braid run --- */
        printf("  Braid background run...\n");
        fflush(stdout);

        /* Reset perturbation on braid background */
        tw_init_packet(tw, 0, x0, y0, z0, kx, ky, kz, sigma, amp);
        tw_compute_forces(tw, BIMODAL);

        double br_cx0, br_cy0, br_cz0, br_E0;
        tw_centroid_all(tw, &br_cx0, &br_cy0, &br_cz0, &br_E0);

        double br_cx_final = br_cx0, br_cy_final = br_cy0;

        for (int s = 0; s < prop_steps; s++) {
            tw_verlet_step(tw, BIMODAL);
            tw_apply_damping_xy(tw);
            if (s % output_every == 0) {
                double cx, cy, cz, E;
                tw_centroid_all(tw, &cx, &cy, &cz, &E);
                fprintf(fdata, "braid\t%.1f\t%.3f\t%.6f\t%.6f\t%.6f\t%.6e\n",
                        b, s*dt, cx, cy, cz, E);
                br_cx_final = cx;
                br_cy_final = cy;
            }
        }

        double br_vx = (br_cx_final - br_cx0) / T_prop;
        double br_vy = (br_cy_final - br_cy0) / T_prop;
        printf("    Braid:  cx %.4f -> %.4f, vx=%.4f, vy=%.6f\n",
               br_cx0, br_cx_final, br_vx, br_vy);

        /* Deflection angle */
        double dvx = br_vx - vac_vx;
        double dvy = br_vy - vac_vy;
        double deflection = atan2(dvy, vac_vx); /* angle of deflection */
        double speed_vac = sqrt(vac_vx*vac_vx + vac_vy*vac_vy);
        double speed_br = sqrt(br_vx*br_vx + br_vy*br_vy);
        double speed_ratio = speed_br / (speed_vac + 1e-30);

        /* Time delay: how much later does the braid-run packet arrive at same x? */
        double time_delay = (vac_cx_final - br_cx_final) / (vac_vx + 1e-30);

        printf("    Deflection: %.6f rad (%.3f deg)\n", deflection, deflection*180/PI);
        printf("    Speed ratio (braid/vacuum): %.6f\n", speed_ratio);
        printf("    Time delay: %.4f\n", time_delay);
        printf("    dvy: %.6e (attractive if negative for b>0)\n\n", dvy);

        fprintf(fsum, "%.1f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.4f\n",
                b, vac_vx, vac_vy, br_vx, br_vy, deflection, speed_ratio, time_delay);
    }

    fclose(fdata);
    fclose(fsum);

    /* Phase 4: Also test wave along z-axis (through the braid) for speed measurement */
    printf("\n=== Z-axis propagation (speed of sound through braid) ===\n");
    {
        double x0 = 0.0, y0 = 0.0, z0 = -10.0;
        double kx = 0.0, ky = 0.0, kz = k_test;

        /* Vacuum */
        TestWave *tw_vac = tw_alloc(N, L);
        tw_init_packet(tw_vac, 0, x0, y0, z0, kx, ky, kz, sigma, amp);
        tw_compute_forces(tw_vac, BIMODAL);
        double vac_cz0, dum1, dum2, dum3;
        tw_centroid_all(tw_vac, &dum1, &dum2, &vac_cz0, &dum3);
        for (int s = 0; s < prop_steps; s++) {
            tw_verlet_step(tw_vac, BIMODAL);
            tw_apply_damping_xy(tw_vac);
        }
        double vac_czf;
        tw_centroid_all(tw_vac, &dum1, &dum2, &vac_czf, &dum3);
        double vz_vac = (vac_czf - vac_cz0) / T_prop;
        tw_free(tw_vac);

        /* Braid */
        tw_init_packet(tw, 0, x0, y0, z0, kx, ky, kz, sigma, amp);
        tw_compute_forces(tw, BIMODAL);
        double br_cz0;
        tw_centroid_all(tw, &dum1, &dum2, &br_cz0, &dum3);
        for (int s = 0; s < prop_steps; s++) {
            tw_verlet_step(tw, BIMODAL);
            tw_apply_damping_xy(tw);
        }
        double br_czf;
        tw_centroid_all(tw, &dum1, &dum2, &br_czf, &dum3);
        double vz_br = (br_czf - br_cz0) / T_prop;

        printf("  Z vacuum speed: %.6f\n", vz_vac);
        printf("  Z braid speed:  %.6f\n", vz_br);
        printf("  Z speed ratio:  %.6f\n", vz_br / (vz_vac + 1e-30));
        printf("  Theory vg = k/omega = %.6f\n", vg_theory);
    }

    /* Clean up */
    tw_free(tw);
    grid_free(g);

    printf("\n=== T10B Complete ===\n");
    printf("Data written to data/t10b_trajectories.tsv and data/t10b_summary.tsv\n");
    return 0;
}
