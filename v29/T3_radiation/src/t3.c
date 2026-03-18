/*  T3: Late-time radiation rate for bimodal braid
 *
 *  Large domain (L=60, N=192), damping only at r_perp > 54.
 *  Tracks E_core (R<8), E_shell (8<R<30), E_total over T=2000.
 *  Uses BIMODAL params with mass^2=2.25 (mass2_override=-1).
 *
 *  Build: gcc -O3 -fopenmp -o t3 src/t3.c -lm
 */

#include "../../src/braid_core.h"

/* Custom damping: only at r_perp > 0.90*L (=54 for L=60) */
static void apply_damping_far(Grid *g) {
    int N = g->N, NN = N * N;
    double dx = g->dx, L = g->L;
    double r_start = 0.90 * L;   /* 54.0 */
    double r_end   = 0.98 * L;   /* 58.8 */
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
                    g->phi[a][idx] *= damp;
                    g->vel[a][idx] *= damp;
                }
            }
        }
    }
}

/* Compute energy in a cylindrical shell r_perp in [r_min, r_max) */
static double energy_in_shell(Grid *g, const double *phys, double r_min, double r_max) {
    int N = g->N, NN = N*N;
    double dx = g->dx, L = g->L, dV = dx*dx*dx;
    double mass2 = phys[14]*phys[14];  /* 2.25 */
    double mu = phys[12], kappa = phys[13], lpw = phys[15];
    double rmin2 = r_min*r_min, rmax2 = r_max*r_max;

    double E = 0;
    #pragma omp parallel for reduction(+:E) schedule(static)
    for (int i = 1; i < N-1; i++) {
        double x = -L + i*dx;
        for (int j = 1; j < N-1; j++) {
            double y = -L + j*dx;
            double rp2 = x*x + y*y;
            if (rp2 < rmin2 || rp2 >= rmax2) continue;
            for (int k = 0; k < N; k++) {
                int idx = i*NN + j*N + k;
                double p0 = g->phi[0][idx], p1 = g->phi[1][idx], p2 = g->phi[2][idx];
                double rho = p0*p0 + p1*p1 + p2*p2;

                double ek = 0, eg = 0;
                for (int a = 0; a < NFIELDS; a++) {
                    ek += 0.5 * g->vel[a][idx] * g->vel[a][idx];
                    int kp = (k+1)%N, km = (k-1+N)%N;
                    double gx = (g->phi[a][idx+NN]-g->phi[a][idx-NN])/(2*dx);
                    double gy = (g->phi[a][idx+N] -g->phi[a][idx-N]) /(2*dx);
                    double gz = (g->phi[a][i*NN+j*N+kp]-g->phi[a][i*NN+j*N+km])/(2*dx);
                    eg += 0.5*(gx*gx + gy*gy + gz*gz);
                }
                double em = 0.5 * mass2 * rho;
                double P = p0*p1*p2;
                double ep = (mu/2.0)*P*P/(1.0+kappa*P*P);
                double epw = lpw*(p0*p1 + p1*p2 + p2*p0);
                E += (ek + eg + em + ep + epw) * dV;
            }
        }
    }
    return E;
}

/* Compute total phi^2 in shell */
static double phi2_in_shell(Grid *g, double r_min, double r_max) {
    int N = g->N, NN = N*N;
    double dx = g->dx, L = g->L, dV = dx*dx*dx;
    double rmin2 = r_min*r_min, rmax2 = r_max*r_max;

    double S = 0;
    #pragma omp parallel for reduction(+:S) schedule(static)
    for (int i = 1; i < N-1; i++) {
        double x = -L + i*dx;
        for (int j = 1; j < N-1; j++) {
            double y = -L + j*dx;
            double rp2 = x*x + y*y;
            if (rp2 < rmin2 || rp2 >= rmax2) continue;
            for (int k = 0; k < N; k++) {
                int idx = i*NN + j*N + k;
                double p0 = g->phi[0][idx], p1 = g->phi[1][idx], p2 = g->phi[2][idx];
                S += (p0*p0 + p1*p1 + p2*p2) * dV;
            }
        }
    }
    return S;
}

/* ================================================================ */

#define MAX_SNAP 200

int main(void) {
    printf("=== T3: Late-Time Radiation Rate ===\n\n");

    bimodal_init_params();
    printf("BIMODAL params:\n");
    for (int d = 0; d < NDIM; d++)
        printf("  %8s = %8.4f\n", PNAME[d], BIMODAL[d]);
    printf("\nmass = %.4f => mass^2 = %.4f\n", BIMODAL[14], BIMODAL[14]*BIMODAL[14]);

    int N = 192;
    double L = 60.0;
    printf("Grid: N=%d, L=%.1f, dx=%.4f\n", N, L, 2.0*L/(N-1));

    Grid *g = grid_alloc(N, L);
    printf("dt = %.6f\n", g->dt);
    printf("Damping: r_perp > %.1f (0.90*L)\n\n", 0.90*L);

    /* Initialize */
    init_braid(g, BIMODAL, -1);
    compute_forces(g, BIMODAL, -1);

    double T_max = 2000.0;
    double T_diag = 25.0;
    int steps_per_diag = (int)(T_diag / g->dt + 0.5);
    int total_steps = (int)(T_max / g->dt + 0.5);
    int n_diag = total_steps / steps_per_diag + 1;
    if (n_diag > MAX_SNAP) n_diag = MAX_SNAP;

    printf("Total steps: %d, diag every %d steps (%d snapshots)\n\n", total_steps, steps_per_diag, n_diag);

    /* Storage for time series */
    double t_arr[MAX_SNAP], E_core_arr[MAX_SNAP], E_shell_arr[MAX_SNAP];
    double E_total_arr[MAX_SNAP], fc_arr[MAX_SNAP], wind_arr[MAX_SNAP], peakP_arr[MAX_SNAP];
    int n_snap = 0;

    /* Print header */
    printf("%8s %12s %12s %12s %8s %8s %8s\n",
           "T", "E_core", "E_shell", "E_total", "fc", "|P|", "wind");
    printf("-------- ------------ ------------ ------------ -------- -------- --------\n");

    /* Initial snapshot */
    {
        double Ec = energy_in_shell(g, BIMODAL, 0.0, 8.0);
        double Es = energy_in_shell(g, BIMODAL, 8.0, 30.0);
        double Et = energy_in_shell(g, BIMODAL, 0.0, 1e6);
        double p2c = phi2_in_shell(g, 0.0, 8.0);
        double p2t = phi2_in_shell(g, 0.0, 1e6);
        double fc = p2c / (p2t + 1e-30);
        double w = compute_winding(g);
        Result res; compute_diagnostics(g, BIMODAL, &res);

        t_arr[0] = 0; E_core_arr[0] = Ec; E_shell_arr[0] = Es;
        E_total_arr[0] = Et; fc_arr[0] = fc; wind_arr[0] = w; peakP_arr[0] = res.peak_P;
        n_snap = 1;

        printf("%8.1f %12.2f %12.2f %12.2f %8.4f %8.4f %8.3f\n",
               0.0, Ec, Es, Et, fc, res.peak_P, w);
    }

    /* Main evolution */
    int step = 0;
    int blown = 0;
    double wall0 = omp_get_wtime();

    for (int diag = 1; diag <= n_diag && !blown; diag++) {
        int target = diag * steps_per_diag;
        if (target > total_steps) target = total_steps;

        for (; step < target; step++) {
            verlet_full_step(g, BIMODAL, -1);
            apply_damping_far(g);
        }

        if (check_blowup(g)) { blown = 1; printf("BLOWUP at step %d\n", step); break; }

        double T = step * g->dt;
        double Ec = energy_in_shell(g, BIMODAL, 0.0, 8.0);
        double Es = energy_in_shell(g, BIMODAL, 8.0, 30.0);
        double Et = energy_in_shell(g, BIMODAL, 0.0, 1e6);
        double p2c = phi2_in_shell(g, 0.0, 8.0);
        double p2t = phi2_in_shell(g, 0.0, 1e6);
        double fc = p2c / (p2t + 1e-30);
        double w = compute_winding(g);
        Result res; compute_diagnostics(g, BIMODAL, &res);

        if (n_snap < MAX_SNAP) {
            t_arr[n_snap] = T; E_core_arr[n_snap] = Ec; E_shell_arr[n_snap] = Es;
            E_total_arr[n_snap] = Et; fc_arr[n_snap] = fc; wind_arr[n_snap] = w;
            peakP_arr[n_snap] = res.peak_P;
            n_snap++;
        }

        printf("%8.1f %12.2f %12.2f %12.2f %8.4f %8.4f %8.3f\n",
               T, Ec, Es, Et, fc, res.peak_P, w);
        fflush(stdout);

        double wall = omp_get_wtime() - wall0;
        double frac = (double)step / total_steps;
        if (diag == 1 || diag % 10 == 0)
            printf("  [wall %.1fs, %.1f%% done, ETA %.0fs]\n", wall, 100*frac, wall*(1-frac)/(frac+1e-30));
    }

    double wall_total = omp_get_wtime() - wall0;
    printf("\nFinished %d steps in %.1f s (%.4f s/step)\n\n", step, wall_total, wall_total/step);

    /* Write TSV */
    FILE *fp = fopen("data/t3_results.tsv", "w");
    if (fp) {
        fprintf(fp, "T\tE_core\tE_shell\tE_total\tfc\tpeakP\twinding\n");
        for (int i = 0; i < n_snap; i++)
            fprintf(fp, "%.2f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.4f\n",
                    t_arr[i], E_core_arr[i], E_shell_arr[i], E_total_arr[i],
                    fc_arr[i], peakP_arr[i], wind_arr[i]);
        fclose(fp);
        printf("Wrote data/t3_results.tsv (%d rows)\n\n", n_snap);
    }

    /* ================================================================
       Analysis: dE_core/dt for T > 500, fit power law
       ================================================================ */
    printf("=== Radiation Analysis (T > 500) ===\n\n");

    /* Find snapshots with T > 500 */
    int i_start = -1;
    for (int i = 0; i < n_snap; i++) {
        if (t_arr[i] >= 500.0) { i_start = i; break; }
    }

    if (i_start < 0 || n_snap - i_start < 4) {
        printf("Not enough data after T=500 for analysis.\n");
        grid_free(g);
        return 0;
    }

    /* Compute dE_core/dt by centered finite difference */
    int n_rate = 0;
    double t_rate[MAX_SNAP], dEdt[MAX_SNAP];

    for (int i = i_start + 1; i < n_snap - 1; i++) {
        double dt_fd = t_arr[i+1] - t_arr[i-1];
        if (dt_fd < 1e-10) continue;
        t_rate[n_rate] = t_arr[i];
        dEdt[n_rate] = (E_core_arr[i+1] - E_core_arr[i-1]) / dt_fd;
        n_rate++;
    }

    if (n_rate < 3) {
        printf("Not enough rate data for fit.\n");
        grid_free(g);
        return 0;
    }

    printf("%8s %14s\n", "T", "dE_core/dt");
    printf("-------- --------------\n");
    for (int i = 0; i < n_rate; i++)
        printf("%8.1f %14.6f\n", t_rate[i], dEdt[i]);

    /* Characterize: mean and trend of |dE/dt| */
    double sum_abs = 0, sum_t = 0, sum_t2 = 0, sum_yt = 0, sum_y = 0;
    int n_neg = 0;
    for (int i = 0; i < n_rate; i++) {
        sum_abs += fabs(dEdt[i]);
        if (dEdt[i] < 0) n_neg++;
    }
    double mean_rate = sum_abs / n_rate;

    /* Linear fit of dE/dt vs t to detect trend */
    for (int i = 0; i < n_rate; i++) {
        sum_t += t_rate[i]; sum_t2 += t_rate[i]*t_rate[i];
        sum_y += dEdt[i]; sum_yt += dEdt[i]*t_rate[i];
    }
    double t_mean = sum_t / n_rate, y_mean = sum_y / n_rate;
    double slope_num = sum_yt - n_rate*t_mean*y_mean;
    double slope_den = sum_t2 - n_rate*t_mean*t_mean;
    double slope = slope_num / (slope_den + 1e-30);

    printf("\nMean |dE_core/dt| = %.6f\n", mean_rate);
    printf("Linear trend slope = %.2e (dE/dt per unit time)\n", slope);
    printf("Fraction negative = %d/%d\n", n_neg, n_rate);

    /* Try power-law fit: log(|dE/dt|) = a + b*log(t)  =>  |dE/dt| ~ t^b */
    double sx2 = 0, sxy = 0, sxl = 0, syl = 0;
    int n_fit = 0;
    for (int i = 0; i < n_rate; i++) {
        if (fabs(dEdt[i]) < 1e-20) continue;
        double lx = log(t_rate[i]);
        double ly = log(fabs(dEdt[i]));
        sxl += lx; syl += ly; sx2 += lx*lx; sxy += lx*ly;
        n_fit++;
    }
    if (n_fit >= 3) {
        double lx_m = sxl/n_fit, ly_m = syl/n_fit;
        double b_pw = (sxy - n_fit*lx_m*ly_m) / (sx2 - n_fit*lx_m*lx_m + 1e-30);
        double a_pw = ly_m - b_pw*lx_m;
        printf("\nPower-law fit: |dE/dt| ~ t^%.3f  (C=%.4e)\n", b_pw, exp(a_pw));

        if (b_pw < -2.0)
            printf("=> EFFECTIVELY STABLE (radiation decreases faster than 1/t^2)\n");
        else if (b_pw < -0.5)
            printf("=> MARGINAL (slow power-law decay of radiation rate)\n");
        else
            printf("=> STEADY RADIATION (braid is actively radiating)\n");
    }

    /* Energy retention */
    if (n_snap >= 2) {
        double E0 = E_core_arr[0];
        double Ef = E_core_arr[n_snap-1];
        printf("\nE_core retention: %.2f%% (initial=%.2f, final=%.2f)\n",
               100.0*Ef/(E0+1e-30), E0, Ef);
        printf("E_total initial=%.2f, final=%.2f (%.2f%% retained)\n",
               E_total_arr[0], E_total_arr[n_snap-1],
               100.0*E_total_arr[n_snap-1]/(E_total_arr[0]+1e-30));
    }

    /* Winding stability */
    printf("\nWinding: initial=%.3f, final=%.3f\n", wind_arr[0], wind_arr[n_snap-1]);

    printf("\n=== T3 Complete ===\n");

    grid_free(g);
    return 0;
}
