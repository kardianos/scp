/*  t10b_v2.c — Effective Metric Extraction (v2: improved tracking)
 *
 *  Key improvements over v1:
 *  - Higher k (k=4) for less dispersive, more relativistic wave packets
 *  - Track peak position (max |δφ|²) instead of energy centroid
 *  - Shorter propagation through smaller region of interest
 *  - Also measure speed along slices (time-of-flight at fixed x)
 *  - Explicit check: does the wave BEND TOWARD the braid (attractive) or away?
 *
 *  Build: gcc -O3 -fopenmp -o t10b_v2 src/t10b_v2.c -lm
 */

#include "../../src/braid_core.h"
#include <float.h>

/* ================================================================
   Test wave on frozen background (same as v1)
   ================================================================ */

typedef struct {
    double *bg[NFIELDS];
    double *dph[NFIELDS];
    double *dvel[NFIELDS];
    double *dacc[NFIELDS];
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

static void tw_set_background(TestWave *tw, Grid *g) {
    int N3 = g->N * g->N * g->N;
    for (int a = 0; a < NFIELDS; a++)
        memcpy(tw->bg[a], g->phi[a], N3 * sizeof(double));
}

static void tw_zero_perturbation(TestWave *tw) {
    int N3 = tw->N * tw->N * tw->N;
    for (int a = 0; a < NFIELDS; a++) {
        memset(tw->dph[a], 0, N3 * sizeof(double));
        memset(tw->dvel[a], 0, N3 * sizeof(double));
        memset(tw->dacc[a], 0, N3 * sizeof(double));
    }
}

/* Initialize wave packet: all 3 components, same polarization (scalar wave) */
static void tw_init_packet(TestWave *tw,
                           double x0, double y0, double z0,
                           double kx, double ky, double kz,
                           double sigma, double amplitude) {
    int N = tw->N, NN = N * N;
    double dx = tw->dx, L = tw->L;
    double inv_2s2 = 1.0 / (2.0 * sigma * sigma);
    double k2 = kx*kx + ky*ky + kz*kz;
    double mass2 = 2.25;
    double omega = sqrt(k2 + mass2);

    tw_zero_perturbation(tw);

    /* Put the wave in component 0 only (clean single-field test) */
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
                double phase = kx * rx + ky * ry + kz * rz;
                tw->dph[0][idx] = env * cos(phase);
                tw->dvel[0][idx] = -omega * env * cos(phase)
                                   + env * sin(phase) * (kx*rx + ky*ry + kz*rz)
                                     * 0; /* just use -omega*cos for simplicity */
            }
        }
    }
}

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

        double b0 = tw->bg[0][idx], b1 = tw->bg[1][idx], b2 = tw->bg[2][idx];
        double P_bg = b0 * b1 * b2;
        double denom_bg = 1.0 + kappa * P_bg * P_bg;
        double denom2 = denom_bg * denom_bg;

        double d0 = tw->dph[0][idx], d1 = tw->dph[1][idx], d2 = tw->dph[2][idx];

        double dP_da[3] = {b1*b2, b0*b2, b0*b1};
        double factor1 = mu * (1.0 - 3.0*kappa*P_bg*P_bg) / (denom2 * denom_bg);
        double factor2 = mu * P_bg / denom2;
        double dd[3] = {d0, d1, d2};

        for (int a = 0; a < NFIELDS; a++) {
            int idx_kp = i*NN + j*N + kp;
            int idx_km = i*NN + j*N + km;
            double lap = (tw->dph[a][idx+NN] + tw->dph[a][idx-NN]
                        + tw->dph[a][idx+N]  + tw->dph[a][idx-N]
                        + tw->dph[a][idx_kp] + tw->dph[a][idx_km]
                        - 6.0 * tw->dph[a][idx]) * idx2;

            double f_lin = 0;
            for (int bb = 0; bb < NFIELDS; bb++) {
                f_lin += dP_da[a] * dP_da[bb] * factor1 * dd[bb];
                if (a != bb) {
                    int c = 3 - a - bb;
                    double bg_c = (c==0)?b0:(c==1)?b1:b2;
                    f_lin += factor2 * bg_c * dd[bb];
                }
            }

            tw->dacc[a][idx] = lap - mass2 * tw->dph[a][idx] - f_lin;
        }
    }
}

static void tw_verlet_step(TestWave *tw, const double *phys) {
    int N3 = tw->N * tw->N * tw->N;
    double hdt = 0.5 * tw->dt;
    for (int a = 0; a < NFIELDS; a++)
        for (int idx = 0; idx < N3; idx++)
            tw->dvel[a][idx] += hdt * tw->dacc[a][idx];
    for (int a = 0; a < NFIELDS; a++)
        for (int idx = 0; idx < N3; idx++)
            tw->dph[a][idx] += tw->dt * tw->dvel[a][idx];
    tw_compute_forces(tw, phys);
    for (int a = 0; a < NFIELDS; a++)
        for (int idx = 0; idx < N3; idx++)
            tw->dvel[a][idx] += hdt * tw->dacc[a][idx];
}

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

/* ================================================================
   Measurement: 1D profile along x at (y=b, z=0) slice
   Returns the x-position of the peak |δφ|²
   ================================================================ */
static double find_peak_x(TestWave *tw, double y_target, double z_target) {
    int N = tw->N, NN = N*N;
    double dx = tw->dx, L = tw->L;

    /* Find nearest j, k */
    int j = (int)((y_target + L) / dx + 0.5);
    int k = (int)((z_target + L) / dx + 0.5);
    if (j < 0) j = 0; if (j >= N) j = N-1;
    if (k < 0) k = 0; if (k >= N) k = N-1;

    double max_amp2 = 0;
    double peak_x = 0;
    for (int i = 1; i < N-1; i++) {
        int idx = i * NN + j * N + k;
        double amp2 = 0;
        for (int a = 0; a < NFIELDS; a++)
            amp2 += tw->dph[a][idx] * tw->dph[a][idx];
        if (amp2 > max_amp2) {
            max_amp2 = amp2;
            peak_x = -L + i * dx;
        }
    }
    return peak_x;
}

/* Measure transverse deflection: find the y-centroid of the wave near expected x position */
static double find_y_centroid(TestWave *tw, double x_expect, double z_target, double window) {
    int N = tw->N, NN = N*N;
    double dx = tw->dx, L = tw->L;
    int k = (int)((z_target + L) / dx + 0.5);
    if (k < 0) k = 0; if (k >= N) k = N-1;

    double wy = 0, wtot = 0;
    for (int i = 1; i < N-1; i++) {
        double x = -L + i * dx;
        if (fabs(x - x_expect) > window) continue;
        for (int j = 1; j < N-1; j++) {
            double y = -L + j * dx;
            int idx = i * NN + j * N + k;
            double amp2 = 0;
            for (int a = 0; a < NFIELDS; a++)
                amp2 += tw->dph[a][idx] * tw->dph[a][idx];
            wy += amp2 * y;
            wtot += amp2;
        }
    }
    return wy / (wtot + 1e-30);
}

/* ================================================================
   Measure: x-profile of wave at z=0 mid-plane for a given y-slice
   This gives us the 1D snapshot to track peak position
   ================================================================ */
static void measure_x_profile(TestWave *tw, double y_target,
                               double *x_arr, double *amp_arr, int *n_out) {
    int N = tw->N, NN = N*N;
    double dx = tw->dx, L = tw->L;
    int j = (int)((y_target + L) / dx + 0.5);
    int k = N / 2; /* z=0 */
    if (j < 0) j = 0; if (j >= N) j = N-1;

    *n_out = N - 2;
    for (int i = 1; i < N-1; i++) {
        int idx = i * NN + j * N + k;
        x_arr[i-1] = -L + i * dx;
        double a2 = 0;
        for (int a = 0; a < NFIELDS; a++)
            a2 += tw->dph[a][idx] * tw->dph[a][idx];
        amp_arr[i-1] = a2;
    }
}

/* ================================================================
   Main: cleaner experiment with multiple approaches
   ================================================================ */

int main(int argc, char **argv) {
    int N = 96;
    double L = 20.0;
    double T_equil = 200.0;
    double T_prop = 30.0;      /* shorter propagation */
    double amp = 0.01;
    double sigma = 1.5;        /* tighter packet */
    double k_test = 4.0;       /* higher k for more relativistic packet */

    double mass2 = 2.25;
    double omega = sqrt(k_test*k_test + mass2);
    double vg = k_test / omega;  /* group velocity */

    printf("=== T10B v2: Effective Metric Extraction ===\n\n");
    printf("Grid: N=%d, L=%.1f, dx=%.4f\n", N, L, 2.0*L/(N-1));
    printf("Test wave: k=%.1f, omega=%.4f, v_group=%.4f\n", k_test, omega, vg);
    printf("Expected transit from x=-12 to x=+12: T=%.1f\n", 24.0/vg);
    printf("Packet sigma=%.1f, amplitude=%.4f\n\n", sigma, amp);

    /* Phase 1: Equilibrate braid */
    printf("Phase 1: Equilibrating braid...\n");
    fflush(stdout);
    bimodal_init_params();
    Grid *g = grid_alloc(N, L);
    init_braid(g, BIMODAL, -1);
    compute_forces(g, BIMODAL, 2.25);

    int equil_steps = (int)(T_equil / g->dt);
    for (int s = 0; s < equil_steps; s++) {
        verlet_full_step(g, BIMODAL, 2.25);
        apply_damping_xy(g);
    }
    Result res;
    snprintf(res.label, sizeof(res.label), "equilibrated");
    compute_diagnostics(g, BIMODAL, &res);
    printf("  fc=%.4f, energy=%.1f\n\n", res.fc, res.energy);

    /* Phase 2: Freeze background */
    printf("Phase 2: Setting up test wave system...\n\n");
    TestWave *tw_br = tw_alloc(N, L);
    tw_set_background(tw_br, g);

    TestWave *tw_vac = tw_alloc(N, L);
    /* tw_vac bg stays zero = vacuum */

    double dt = tw_br->dt;
    int prop_steps = (int)(T_prop / dt);

    /* Open output files */
    FILE *fdata = fopen("data/t10b_v2_timeseries.tsv", "w");
    FILE *fsum = fopen("data/t10b_v2_summary.tsv", "w");
    fprintf(fdata, "b\tmode\tt\tpeak_x\ty_centroid\tpeak_amp2\n");
    fprintf(fsum, "direction\tb\tvac_speed\tbraid_speed\tspeed_ratio\tdeflection_y\tnotes\n");

    /* ============================================================
       Experiment A: X-directed waves at various impact parameters
       ============================================================ */
    printf("=== Experiment A: X-directed test waves ===\n");
    double b_vals[] = {0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 14.0};
    int nb = 7;

    for (int ib = 0; ib < nb; ib++) {
        double b = b_vals[ib];
        double x0 = -12.0;
        double kx = k_test, ky = 0, kz = 0;

        printf("  b=%.1f: ", b);
        fflush(stdout);

        /* Vacuum run */
        tw_init_packet(tw_vac, x0, b, 0, kx, ky, kz, sigma, amp);
        tw_compute_forces(tw_vac, BIMODAL);

        double vac_peak0 = find_peak_x(tw_vac, b, 0);
        double vac_peak_final = vac_peak0;
        double vac_y0 = find_y_centroid(tw_vac, x0, 0, 4.0);

        int meas_count = 0;
        int meas_interval = prop_steps / 50;
        if (meas_interval < 1) meas_interval = 1;

        double vac_peaks[100], vac_times[100];
        int n_vac_meas = 0;

        for (int s = 0; s < prop_steps; s++) {
            tw_verlet_step(tw_vac, BIMODAL);
            tw_apply_damping_xy(tw_vac);
            if (s % meas_interval == 0 && n_vac_meas < 100) {
                double px = find_peak_x(tw_vac, b, 0);
                double yc = find_y_centroid(tw_vac, px, 0, 4.0);
                double t = s * dt;
                fprintf(fdata, "%.1f\tvacuum\t%.4f\t%.4f\t%.6f\t0\n", b, t, px, yc);
                vac_peaks[n_vac_meas] = px;
                vac_times[n_vac_meas] = t;
                n_vac_meas++;
                vac_peak_final = px;
            }
        }
        double vac_yfinal = find_y_centroid(tw_vac, vac_peak_final, 0, 4.0);

        /* Vacuum speed from linear fit of peak positions */
        double vac_speed = 0;
        if (n_vac_meas > 2) {
            /* Simple: use first third and last third */
            int n3 = n_vac_meas / 3;
            double t1 = 0, x1 = 0, t2 = 0, x2 = 0;
            for (int i = 0; i < n3; i++) {
                t1 += vac_times[i]; x1 += vac_peaks[i];
            }
            for (int i = n_vac_meas - n3; i < n_vac_meas; i++) {
                t2 += vac_times[i]; x2 += vac_peaks[i];
            }
            t1 /= n3; x1 /= n3; t2 /= n3; x2 /= n3;
            vac_speed = (x2 - x1) / (t2 - t1 + 1e-30);
        }

        /* Braid run */
        tw_init_packet(tw_br, x0, b, 0, kx, ky, kz, sigma, amp);
        tw_compute_forces(tw_br, BIMODAL);

        double br_peaks[100], br_times[100];
        int n_br_meas = 0;

        for (int s = 0; s < prop_steps; s++) {
            tw_verlet_step(tw_br, BIMODAL);
            tw_apply_damping_xy(tw_br);
            if (s % meas_interval == 0 && n_br_meas < 100) {
                double px = find_peak_x(tw_br, b, 0);
                double yc = find_y_centroid(tw_br, px, 0, 4.0);
                double t = s * dt;
                fprintf(fdata, "%.1f\tbraid\t%.4f\t%.4f\t%.6f\t0\n", b, t, px, yc);
                br_peaks[n_br_meas] = px;
                br_times[n_br_meas] = t;
                n_br_meas++;
            }
        }
        double br_yfinal = find_y_centroid(tw_br, br_peaks[n_br_meas-1], 0, 4.0);

        /* Braid speed */
        double br_speed = 0;
        if (n_br_meas > 2) {
            int n3 = n_br_meas / 3;
            double t1 = 0, x1 = 0, t2 = 0, x2 = 0;
            for (int i = 0; i < n3; i++) {
                t1 += br_times[i]; x1 += br_peaks[i];
            }
            for (int i = n_br_meas - n3; i < n_br_meas; i++) {
                t2 += br_times[i]; x2 += br_peaks[i];
            }
            t1 /= n3; x1 /= n3; t2 /= n3; x2 /= n3;
            br_speed = (x2 - x1) / (t2 - t1 + 1e-30);
        }

        double speed_ratio = br_speed / (vac_speed + 1e-30);
        double defl_y = br_yfinal - vac_yfinal;

        printf("vac_speed=%.4f, braid_speed=%.4f, ratio=%.4f, defl_y=%.4f",
               vac_speed, br_speed, speed_ratio, defl_y);
        if (b > 0 && defl_y < 0)
            printf(" [ATTRACTIVE!]");
        else if (b > 0 && defl_y > 0)
            printf(" [repulsive]");
        printf("\n");

        fprintf(fsum, "x\t%.1f\t%.6f\t%.6f\t%.6f\t%.6f\t%s\n",
                b, vac_speed, br_speed, speed_ratio, defl_y,
                (b > 0 && defl_y < 0) ? "attractive" : "repulsive_or_zero");
    }

    /* ============================================================
       Experiment B: Z-directed wave (through braid core)
       ============================================================ */
    printf("\n=== Experiment B: Z-directed test waves ===\n");
    {
        /* Z-direction: use periodic BC naturally */
        double z0 = -10.0;
        double kz = k_test;

        printf("  Through center (r=0): ");
        fflush(stdout);

        /* Vacuum */
        tw_init_packet(tw_vac, 0, 0, z0, 0, 0, kz, sigma, amp);
        tw_compute_forces(tw_vac, BIMODAL);

        /* For z-directed, track peak in z */
        /* Find peak z by scanning along z at x=0, y=0 */
        int ic = tw_vac->N / 2, jc = tw_vac->N / 2;
        int NN = tw_vac->N * tw_vac->N;
        double dx_g = tw_vac->dx, L_g = tw_vac->L;

        double vac_zpk[100], vac_zt[100];
        int nvz = 0;
        int meas_int = prop_steps / 50;
        if (meas_int < 1) meas_int = 1;

        for (int s = 0; s < prop_steps; s++) {
            tw_verlet_step(tw_vac, BIMODAL);
            tw_apply_damping_xy(tw_vac);
            if (s % meas_int == 0 && nvz < 100) {
                double maxamp = 0, zpk = z0;
                for (int kk = 0; kk < tw_vac->N; kk++) {
                    int idx = ic * NN + jc * tw_vac->N + kk;
                    double a2 = 0;
                    for (int a = 0; a < NFIELDS; a++)
                        a2 += tw_vac->dph[a][idx] * tw_vac->dph[a][idx];
                    if (a2 > maxamp) {
                        maxamp = a2;
                        zpk = -L_g + kk * dx_g;
                    }
                }
                vac_zpk[nvz] = zpk;
                vac_zt[nvz] = s * dt;
                nvz++;
            }
        }

        /* Braid */
        tw_init_packet(tw_br, 0, 0, z0, 0, 0, kz, sigma, amp);
        tw_compute_forces(tw_br, BIMODAL);

        double br_zpk[100], br_zt[100];
        int nbz = 0;

        for (int s = 0; s < prop_steps; s++) {
            tw_verlet_step(tw_br, BIMODAL);
            tw_apply_damping_xy(tw_br);
            if (s % meas_int == 0 && nbz < 100) {
                double maxamp = 0, zpk = z0;
                for (int kk = 0; kk < tw_br->N; kk++) {
                    int idx = ic * NN + jc * tw_br->N + kk;
                    double a2 = 0;
                    for (int a = 0; a < NFIELDS; a++)
                        a2 += tw_br->dph[a][idx] * tw_br->dph[a][idx];
                    if (a2 > maxamp) {
                        maxamp = a2;
                        zpk = -L_g + kk * dx_g;
                    }
                }
                br_zpk[nbz] = zpk;
                br_zt[nbz] = s * dt;
                nbz++;
            }
        }

        /* Speeds */
        double vac_vz = 0, br_vz = 0;
        if (nvz > 4) {
            int n3 = nvz / 3;
            double t1=0,z1=0,t2=0,z2=0;
            for(int i=0;i<n3;i++){t1+=vac_zt[i];z1+=vac_zpk[i];}
            for(int i=nvz-n3;i<nvz;i++){t2+=vac_zt[i];z2+=vac_zpk[i];}
            t1/=n3;z1/=n3;t2/=n3;z2/=n3;
            vac_vz=(z2-z1)/(t2-t1+1e-30);
        }
        if (nbz > 4) {
            int n3 = nbz / 3;
            double t1=0,z1=0,t2=0,z2=0;
            for(int i=0;i<n3;i++){t1+=br_zt[i];z1+=br_zpk[i];}
            for(int i=nbz-n3;i<nbz;i++){t2+=br_zt[i];z2+=br_zpk[i];}
            t1/=n3;z1/=n3;t2/=n3;z2/=n3;
            br_vz=(z2-z1)/(t2-t1+1e-30);
        }

        printf("vac_vz=%.4f, braid_vz=%.4f, ratio=%.4f\n", vac_vz, br_vz,
               br_vz/(vac_vz+1e-30));
        fprintf(fsum, "z_core\t0.0\t%.6f\t%.6f\t%.6f\t0\tthrough_core\n",
                vac_vz, br_vz, br_vz/(vac_vz+1e-30));
    }

    /* ============================================================
       Experiment C: Measure effective potential V_eff(r)
       Launch wave at many b, measure speed change -> extract n(r) = c/v(r)
       ============================================================ */
    printf("\n=== Experiment C: Effective refractive index n(r) ===\n");
    printf("  (Extracted from speed ratio at each impact parameter)\n");
    printf("  b\tn_eff\t(n>1 = slower = attractive lensing)\n");

    /* Use the speed ratios already computed above */
    /* Re-run quickly with focus on speed measurement at more b values */
    double b_fine[] = {0, 1, 2, 3, 4, 5, 6, 8, 10, 12, 14, 16};
    int nb_fine = 12;
    double n_eff[12];

    for (int ib = 0; ib < nb_fine; ib++) {
        double b = b_fine[ib];
        double x0 = -12.0;

        /* Just do a quick speed measurement: track peak at t=0 and t=T */
        tw_init_packet(tw_vac, x0, b, 0, k_test, 0, 0, sigma, amp);
        tw_compute_forces(tw_vac, BIMODAL);
        double vac_x0 = find_peak_x(tw_vac, b, 0);
        for (int s = 0; s < prop_steps; s++) {
            tw_verlet_step(tw_vac, BIMODAL);
            tw_apply_damping_xy(tw_vac);
        }
        double vac_xf = find_peak_x(tw_vac, b, 0);
        double vac_v = (vac_xf - vac_x0) / T_prop;

        tw_init_packet(tw_br, x0, b, 0, k_test, 0, 0, sigma, amp);
        tw_compute_forces(tw_br, BIMODAL);
        double br_x0 = find_peak_x(tw_br, b, 0);
        for (int s = 0; s < prop_steps; s++) {
            tw_verlet_step(tw_br, BIMODAL);
            tw_apply_damping_xy(tw_br);
        }
        double br_xf = find_peak_x(tw_br, b, 0);
        double br_v = (br_xf - br_x0) / T_prop;

        n_eff[ib] = vac_v / (br_v + 1e-30);
        printf("  %.1f\t%.6f\n", b, n_eff[ib]);
    }

    fclose(fdata);
    fclose(fsum);

    /* Summary */
    printf("\n=== SUMMARY ===\n");
    printf("Theory group velocity: vg = k/omega = %.4f\n", vg);
    printf("Key question: does n(r) > 1 near the braid? (slower = attractive = gravity-like)\n");
    printf("Key question: does wave deflect TOWARD braid? (dvy < 0 for b > 0)\n\n");

    int any_attractive = 0;
    for (int ib = 0; ib < nb_fine; ib++) {
        if (n_eff[ib] > 1.001) { any_attractive = 1; break; }
    }
    if (any_attractive) {
        printf("RESULT: Waves ARE slower near the braid (n > 1) → ATTRACTIVE effective metric\n");
        printf("This is the GRAVITY signal at the field level.\n");
    } else {
        int any_repulsive = 0;
        for (int ib = 0; ib < nb_fine; ib++) {
            if (n_eff[ib] < 0.999) { any_repulsive = 1; break; }
        }
        if (any_repulsive) {
            printf("RESULT: Waves are FASTER near the braid (n < 1) → REPULSIVE effective metric\n");
            printf("This is anti-gravity. The braid defocuses waves.\n");
        } else {
            printf("RESULT: No significant speed change detected. The braid does NOT curve the field.\n");
            printf("Gravity is NOT emergent from this Lagrangian at the linear level.\n");
        }
    }

    tw_free(tw_br);
    tw_free(tw_vac);
    grid_free(g);

    return 0;
}
