/*  t10f.c — Massless Mediator Test
 *
 *  From T10B: no 1/r because all modes are massive (m=1.5).
 *  Question: what if we make one field massless while keeping others massive?
 *
 *  The massless field would propagate as 1/r (in 3D), potentially mediating
 *  a long-range force. The triple-product coupling sources the massless field
 *  through the braid.
 *
 *  Test configurations:
 *  1. m₀=0, m₁=m₂=1.5: φ₀ is massless mediator
 *  2. m₀=m₁=m₂=0 with separate confinement: all massless, triple product confines
 *  3. m₀=0.1, m₁=m₂=1.5: slightly massive (long Yukawa range)
 *
 *  For each: init braid, equilibrate, measure:
 *  - Does braid survive with mixed masses?
 *  - What is the far-field profile of the massless field?
 *  - Does it decay as 1/r or faster?
 *  - Hessian eigenvalues: is one mode truly massless on the background?
 *
 *  Build: gcc -O3 -fopenmp -o t10f src/t10f.c -lm
 */

#include "../../src/braid_core.h"
#include <float.h>

/* Modified force computation with per-field mass */
static void compute_forces_mixed_mass(Grid *g, const double *phys,
                                       double m0, double m1, double m2) {
    int N = g->N, NN = N*N, N3 = N*N*N;
    double idx2 = 1.0 / (g->dx * g->dx);
    double mu    = phys[12];
    double kappa = phys[13];
    double mass2[3] = {m0*m0, m1*m1, m2*m2};

    #pragma omp parallel for schedule(static)
    for (int idx = 0; idx < N3; idx++) {
        int i = idx / NN, j = (idx / N) % N, k = idx % N;
        if (i < 1 || i >= N-1 || j < 1 || j >= N-1) {
            g->acc[0][idx] = g->acc[1][idx] = g->acc[2][idx] = 0;
            continue;
        }
        int kp = (k+1)%N, km = (k-1+N)%N;
        int idx_kp = i*NN + j*N + kp;
        int idx_km = i*NN + j*N + km;

        double p0 = g->phi[0][idx], p1 = g->phi[1][idx], p2 = g->phi[2][idx];
        double P  = p0 * p1 * p2;
        double denom = 1.0 + kappa * P * P;
        double mu_P_d2 = mu * P / (denom * denom);

        for (int a = 0; a < NFIELDS; a++) {
            double lap = (g->phi[a][idx+NN] + g->phi[a][idx-NN]
                        + g->phi[a][idx+N]  + g->phi[a][idx-N]
                        + g->phi[a][idx_kp] + g->phi[a][idx_km]
                        - 6.0 * g->phi[a][idx]) * idx2;
            double dPda = (a==0) ? p1*p2 : (a==1) ? p0*p2 : p0*p1;
            double f_triple = mu_P_d2 * dPda;
            g->acc[a][idx] = lap - mass2[a] * g->phi[a][idx] - f_triple;
        }
    }
}

/* Verlet with mixed mass */
static void verlet_mixed(Grid *g, const double *phys,
                         double m0, double m1, double m2) {
    int N3 = g->N * g->N * g->N;
    double hdt = 0.5 * g->dt;
    for (int a = 0; a < NFIELDS; a++)
        for (int idx = 0; idx < N3; idx++)
            g->vel[a][idx] += hdt * g->acc[a][idx];
    for (int a = 0; a < NFIELDS; a++)
        for (int idx = 0; idx < N3; idx++)
            g->phi[a][idx] += g->dt * g->vel[a][idx];
    compute_forces_mixed_mass(g, phys, m0, m1, m2);
    for (int a = 0; a < NFIELDS; a++)
        for (int idx = 0; idx < N3; idx++)
            g->vel[a][idx] += hdt * g->acc[a][idx];
}

/* Measure radial profile of each field (azimuthal average at z=0) */
static void measure_radial_profile(Grid *g, int nr, double dr,
                                    double *r_arr, double *rms_arr[3]) {
    int N = g->N, NN = N*N;
    double dx = g->dx, L = g->L;
    int k_mid = N/2;

    int *count = calloc(nr, sizeof(int));
    for (int e = 0; e < 3; e++)
        memset(rms_arr[e], 0, nr * sizeof(double));

    for (int i = 1; i < N-1; i++) {
        double x = -L + i * dx;
        for (int j = 1; j < N-1; j++) {
            double y = -L + j * dx;
            double r = sqrt(x*x + y*y);
            int ir = (int)(r / dr);
            if (ir >= nr) continue;
            int idx = i*NN + j*N + k_mid;
            for (int a = 0; a < 3; a++)
                rms_arr[a][ir] += g->phi[a][idx] * g->phi[a][idx];
            count[ir]++;
        }
    }

    for (int ir = 0; ir < nr; ir++) {
        r_arr[ir] = (ir + 0.5) * dr;
        if (count[ir] > 0) {
            for (int a = 0; a < 3; a++)
                rms_arr[a][ir] = sqrt(rms_arr[a][ir] / count[ir]);
        }
    }
    free(count);
}

int main(int argc, char **argv) {
    int N = 96;
    double L = 30.0; /* larger box to see 1/r decay */
    double T_equil = 300.0;

    printf("=== T10F: Massless Mediator Test ===\n\n");
    printf("Grid: N=%d, L=%.1f\n\n", N, L);

    bimodal_init_params();

    FILE *fout = fopen("data/t10f_profiles.tsv", "w");
    fprintf(fout, "config\tr\trms_0\trms_1\trms_2\n");

    FILE *fsum = fopen("data/t10f_summary.tsv", "w");
    fprintf(fsum, "config\tm0\tm1\tm2\tfc\tenergy\tfit_power\tfit_range\tnotes\n");

    /* Test configurations */
    struct {
        const char *name;
        double m0, m1, m2;
    } configs[] = {
        {"standard",  1.5, 1.5, 1.5},
        {"m0_zero",   0.0, 1.5, 1.5},
        {"m0_small",  0.1, 1.5, 1.5},
        {"m0_half",   0.75, 1.5, 1.5},
        {"all_small", 0.3, 0.3, 0.3},
    };
    int nconfigs = 5;

    int nr = 60;
    double dr = L / nr;
    double *r_arr = malloc(nr * sizeof(double));
    double *rms[3];
    for (int a = 0; a < 3; a++)
        rms[a] = malloc(nr * sizeof(double));

    for (int ic = 0; ic < nconfigs; ic++) {
        printf("=== Config: %s (m = %.2f, %.2f, %.2f) ===\n",
               configs[ic].name, configs[ic].m0, configs[ic].m1, configs[ic].m2);
        fflush(stdout);

        Grid *g = grid_alloc(N, L);

        /* Initialize with bimodal params but modified mass for init */
        /* Use m_init = max(m0,m1,m2) to set the initial oscillation frequency */
        double m_init = configs[ic].m0;
        if (configs[ic].m1 > m_init) m_init = configs[ic].m1;
        if (configs[ic].m2 > m_init) m_init = configs[ic].m2;
        init_braid(g, BIMODAL, m_init);

        /* Equilibrate with mixed masses */
        compute_forces_mixed_mass(g, BIMODAL,
                                   configs[ic].m0, configs[ic].m1, configs[ic].m2);
        int equil_steps = (int)(T_equil / g->dt);
        for (int s = 0; s < equil_steps; s++) {
            verlet_mixed(g, BIMODAL,
                        configs[ic].m0, configs[ic].m1, configs[ic].m2);
            apply_damping_xy(g);

            if (check_blowup(g)) {
                printf("  BLOWUP at t=%.0f\n", s * g->dt);
                break;
            }
        }

        /* Diagnostics */
        Result res;
        compute_diagnostics(g, BIMODAL, &res);
        printf("  fc=%.4f, E=%.1f, |P|=%.4f, w=%.3f\n",
               res.fc, res.energy, res.peak_P, res.winding);

        /* Radial profile */
        measure_radial_profile(g, nr, dr, r_arr, rms);
        for (int ir = 0; ir < nr; ir++) {
            fprintf(fout, "%s\t%.4f\t%.6e\t%.6e\t%.6e\n",
                    configs[ic].name, r_arr[ir], rms[0][ir], rms[1][ir], rms[2][ir]);
        }

        /* Fit power law to the massless field's far-field decay */
        /* For φ₀ (the potentially massless one), fit log(rms) = A - n*log(r) */
        int fit_start = nr/3, fit_end = 2*nr/3;
        double sum_logr = 0, sum_logrms = 0, sum_logr2 = 0, sum_lr_lrms = 0;
        int nfit = 0;
        for (int ir = fit_start; ir < fit_end; ir++) {
            if (rms[0][ir] > 1e-10) {
                double lr = log(r_arr[ir]);
                double la = log(rms[0][ir]);
                sum_logr += lr;
                sum_logrms += la;
                sum_logr2 += lr * lr;
                sum_lr_lrms += lr * la;
                nfit++;
            }
        }
        double power = 0, fit_range = 0;
        if (nfit > 2) {
            power = -(nfit * sum_lr_lrms - sum_logr * sum_logrms) /
                     (nfit * sum_logr2 - sum_logr * sum_logr2 + 1e-30);
            /* Actually use proper formula */
            double denom = nfit * sum_logr2 - sum_logr * sum_logr;
            if (fabs(denom) > 1e-30)
                power = -(nfit * sum_lr_lrms - sum_logr * sum_logrms) / denom;
            fit_range = r_arr[fit_end] - r_arr[fit_start];
        }

        printf("  Far-field power law fit (φ₀): rms ~ r^(-%.2f) over r=[%.1f,%.1f]\n",
               power, r_arr[fit_start], r_arr[fit_end-1]);

        /* Also print explicit decay at a few radii */
        printf("  φ₀ profile: ");
        for (int ir = 0; ir < nr; ir += nr/6) {
            if (rms[0][ir] > 1e-15)
                printf("r=%.1f:%.3e  ", r_arr[ir], rms[0][ir]);
        }
        printf("\n");

        /* Check: does the far field of the massless component show 1/r? */
        if (configs[ic].m0 < 0.01 && power > 0.8 && power < 1.5) {
            printf("  >>> 1/r DECAY DETECTED for massless field! <<<\n");
        } else if (configs[ic].m0 < 0.01) {
            printf("  Power law: %.2f (1/r would be 1.0)\n", power);
        }

        fprintf(fsum, "%s\t%.2f\t%.2f\t%.2f\t%.4f\t%.1f\t%.2f\t%.1f\t%s\n",
                configs[ic].name, configs[ic].m0, configs[ic].m1, configs[ic].m2,
                res.fc, res.energy, power, fit_range,
                res.fc > 0.5 ? "braid_survives" : "braid_lost");

        printf("\n");
        grid_free(g);
    }

    fclose(fout);
    fclose(fsum);

    free(r_arr);
    for (int a = 0; a < 3; a++) free(rms[a]);

    printf("=== T10F Complete ===\n");
    printf("Data in data/t10f_*.tsv\n");
    return 0;
}
