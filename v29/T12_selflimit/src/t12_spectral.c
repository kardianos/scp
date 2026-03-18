/*  t12_spectral.c — Frequency-resolved depletion analysis
 *
 *  The hypothesis: the braid absorbs at its resonant frequency omega_0
 *  and radiates at other frequencies. The depletion at omega_0 follows 1/r,
 *  but the TOTAL depletion (all frequencies) is partially cancelled by
 *  re-radiation, giving the observed ~1/r^1.7.
 *
 *  Method:
 *  Run two simulations (braid+bg vs bg-only control) for T=600.
 *  Record phi_a(t) at probe points at various radii.
 *  FFT the time series to get the spectral power at each radius.
 *  Compute the spectral depletion: delta_P(omega, r) = P_braid - P_control.
 *  If delta_P(omega_0, r) ~ 1/r: the hypothesis is confirmed.
 *
 *  Build: gcc -O3 -fopenmp -o t12_spectral src/t12_spectral.c -lm
 */

#include "../../src/braid_core.h"
#include <float.h>

/* FFT — simple DFT for moderate N (no library needed) */
static void dft_power(const double *signal, int N, double dt,
                      double *freq, double *power, int N_freq) {
    /* Compute power spectrum for frequencies 0..N_freq-1 */
    for (int f = 0; f < N_freq; f++) {
        freq[f] = f / (N * dt);  /* frequency in Hz (code units) */
        double re = 0, im = 0;
        for (int t = 0; t < N; t++) {
            double phase = -2.0 * PI * f * t / N;
            re += signal[t] * cos(phase);
            im += signal[t] * sin(phase);
        }
        power[f] = (re * re + im * im) / (N * N);
    }
}

/* M7 two-component force computation */
static double *B_phi[NFIELDS], *B_vel[NFIELDS], *B_acc[NFIELDS];
static double G_COUP = 0.01;

static void compute_forces_m7(Grid *g, const double *phys) {
    int N = g->N, NN = N*N, N3 = N*N*N;
    double idx2 = 1.0 / (g->dx * g->dx);
    double mu = phys[12], kappa = phys[13];
    double mass2 = phys[14] * phys[14];
    double lpw = phys[15];

    #pragma omp parallel for schedule(static)
    for (int idx = 0; idx < N3; idx++) {
        int i = idx / NN, j = (idx / N) % N, k = idx % N;
        if (i < 1 || i >= N-1 || j < 1 || j >= N-1) {
            for (int a = 0; a < NFIELDS; a++) {
                g->acc[a][idx] = 0;
                B_acc[a][idx] = 0;
            }
            continue;
        }
        int kp = (k+1)%N, km = (k-1+N)%N;
        int idx_kp = i*NN + j*N + kp, idx_km = i*NN + j*N + km;

        /* S fields (braid) */
        double s0 = g->phi[0][idx], s1 = g->phi[1][idx], s2 = g->phi[2][idx];
        double S2 = s0*s0 + s1*s1 + s2*s2;
        double P = s0 * s1 * s2;
        double denom = 1.0 + kappa * P * P;
        double mu_P_d2 = mu * P / (denom * denom);

        /* B fields (background) */
        double b0 = B_phi[0][idx], b1 = B_phi[1][idx], b2 = B_phi[2][idx];
        double B2 = b0*b0 + b1*b1 + b2*b2;

        for (int a = 0; a < NFIELDS; a++) {
            /* S Laplacian */
            double lap_s = (g->phi[a][idx+NN] + g->phi[a][idx-NN]
                          + g->phi[a][idx+N]  + g->phi[a][idx-N]
                          + g->phi[a][idx_kp] + g->phi[a][idx_km]
                          - 6.0 * g->phi[a][idx]) * idx2;
            double dPda = (a==0) ? s1*s2 : (a==1) ? s0*s2 : s0*s1;
            double f_triple = mu_P_d2 * dPda;
            double f_pw = lpw * (g->phi[(a+1)%3][idx] + g->phi[(a+2)%3][idx]);

            g->acc[a][idx] = lap_s - mass2 * g->phi[a][idx]
                           - f_triple - f_pw
                           - G_COUP * B2 * g->phi[a][idx];

            /* B Laplacian */
            double lap_b = (B_phi[a][idx+NN] + B_phi[a][idx-NN]
                          + B_phi[a][idx+N]  + B_phi[a][idx-N]
                          + B_phi[a][idx_kp] + B_phi[a][idx_km]
                          - 6.0 * B_phi[a][idx]) * idx2;

            B_acc[a][idx] = lap_b - mass2 * B_phi[a][idx]
                          - G_COUP * S2 * B_phi[a][idx];
        }
    }
}

static void verlet_step_m7(Grid *g, const double *phys) {
    int N3 = g->N * g->N * g->N;
    double hdt = 0.5 * g->dt, dt = g->dt;

    /* Kick S and B */
    for (int a = 0; a < NFIELDS; a++)
        for (int idx = 0; idx < N3; idx++) {
            g->vel[a][idx] += hdt * g->acc[a][idx];
            B_vel[a][idx]  += hdt * B_acc[a][idx];
        }
    /* Drift S and B */
    for (int a = 0; a < NFIELDS; a++)
        for (int idx = 0; idx < N3; idx++) {
            g->phi[a][idx] += dt * g->vel[a][idx];
            B_phi[a][idx]  += dt * B_vel[a][idx];
        }
    /* Force */
    compute_forces_m7(g, phys);
    /* Kick */
    for (int a = 0; a < NFIELDS; a++)
        for (int idx = 0; idx < N3; idx++) {
            g->vel[a][idx] += hdt * g->acc[a][idx];
            B_vel[a][idx]  += hdt * B_acc[a][idx];
        }
}

static void apply_edge_damping_m7(Grid *g) {
    int N = g->N, NN = N*N;
    double dx = g->dx, L = g->L;
    double rs = 0.85*L, re = 0.98*L, idr = 1.0/(re-rs+1e-30);
    for (int i = 0; i < N; i++) {
        double x = -L + i*dx;
        for (int j = 0; j < N; j++) {
            double y = -L + j*dx;
            double rp = sqrt(x*x+y*y);
            if (rp <= rs) continue;
            double f = (rp-rs)*idr; if (f>1) f=1;
            double d = 1.0 - 0.5*f*f;
            for (int kk = 0; kk < N; kk++) {
                int idx = i*NN + j*N + kk;
                for (int a = 0; a < NFIELDS; a++) {
                    g->phi[a][idx] *= d; g->vel[a][idx] *= d;
                    B_phi[a][idx]  *= d; B_vel[a][idx]  *= d;
                }
            }
        }
    }
}

/* ================================================================ */

#define N_PROBES 8
#define MAX_SAMPLES 1200

int main(int argc, char **argv) {
    omp_set_num_threads(16);
    bimodal_init_params();

    int N = 128;
    double L = 30.0;
    double T_total = 600.0;
    double A_bg = 0.1;

    printf("=== T12 Spectral Analysis: Frequency-Resolved Depletion ===\n");
    printf("N=%d, L=%.0f, T=%.0f, A_bg=%.2f\n", N, L, T_total, A_bg);

    double dx = 2.0*L/(N-1);
    double dt = 0.20*dx;
    double mass2 = BIMODAL[14]*BIMODAL[14];
    double k_bg = PI/L;
    double omega_bg = sqrt(k_bg*k_bg + mass2);

    printf("dx=%.4f dt=%.5f\n", dx, dt);
    printf("k_bg=%.4f omega_bg=%.4f (this is omega_0, the resonant frequency)\n",
           k_bg, omega_bg);

    /* Probe points at various radii along x-axis (y=0, z=z_mid) */
    double probe_r[N_PROBES] = {0, 2, 4, 6, 8, 12, 16, 20};
    int probe_i[N_PROBES];
    int ic = N/2, jc = N/2, kc = N/2;  /* center indices */
    int NN = N*N;

    for (int p = 0; p < N_PROBES; p++) {
        probe_i[p] = ic + (int)(probe_r[p] / dx + 0.5);
        if (probe_i[p] >= N-1) probe_i[p] = N-2;
        printf("Probe %d: r=%.0f -> grid i=%d\n", p, probe_r[p], probe_i[p]);
    }

    int n_steps = (int)(T_total / dt);
    /* Record every few steps to get enough temporal resolution for FFT */
    int record_every = 4;  /* record every 4 steps */
    int n_samples = n_steps / record_every;
    if (n_samples > MAX_SAMPLES) { record_every = n_steps / MAX_SAMPLES; n_samples = n_steps / record_every; }
    double dt_sample = dt * record_every;

    printf("Steps: %d, recording every %d -> %d samples, dt_sample=%.5f\n",
           n_steps, record_every, n_samples, dt_sample);
    printf("Nyquist freq: %.4f, freq resolution: %.6f\n",
           0.5/dt_sample, 1.0/(n_samples*dt_sample));

    /* Allocate probe time series — for BRAID and CONTROL runs */
    /* Store phi_0 at each probe for both S and B fields */
    double *ts_S_braid[N_PROBES], *ts_B_braid[N_PROBES];
    double *ts_S_ctrl[N_PROBES],  *ts_B_ctrl[N_PROBES];
    double *ts_rho_braid[N_PROBES], *ts_rho_ctrl[N_PROBES];
    for (int p = 0; p < N_PROBES; p++) {
        ts_S_braid[p] = calloc(n_samples, sizeof(double));
        ts_B_braid[p] = calloc(n_samples, sizeof(double));
        ts_S_ctrl[p]  = calloc(n_samples, sizeof(double));
        ts_B_ctrl[p]  = calloc(n_samples, sizeof(double));
        ts_rho_braid[p] = calloc(n_samples, sizeof(double));
        ts_rho_ctrl[p]  = calloc(n_samples, sizeof(double));
    }

    /* ============================================================
       RUN 1: BRAID + BACKGROUND (M7)
       ============================================================ */
    printf("\n=== Run 1: Braid + Background (M7) ===\n");

    Grid *g = grid_alloc(N, L);
    int N3 = N*N*N;
    for (int a = 0; a < NFIELDS; a++) {
        B_phi[a] = calloc(N3, sizeof(double));
        B_vel[a] = calloc(N3, sizeof(double));
        B_acc[a] = calloc(N3, sizeof(double));
    }

    /* Init S = braid */
    init_braid(g, BIMODAL, -1);

    /* Init B = background plane wave */
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            for (int kk = 0; kk < N; kk++) {
                int idx = i*NN + j*N + kk;
                double z = -L + kk*dx;
                for (int a = 0; a < NFIELDS; a++) {
                    double ph = k_bg*z + 2.0*PI*a/3.0;
                    B_phi[a][idx] = A_bg * cos(ph);
                    B_vel[a][idx] = omega_bg * A_bg * sin(ph);
                }
            }

    compute_forces_m7(g, BIMODAL);
    int sample = 0;

    for (int step = 0; step < n_steps; step++) {
        verlet_step_m7(g, BIMODAL);
        apply_edge_damping_m7(g);

        if (step % record_every == 0 && sample < n_samples) {
            for (int p = 0; p < N_PROBES; p++) {
                int idx = probe_i[p]*NN + jc*N + kc;
                /* Record field values */
                ts_S_braid[p][sample] = g->phi[0][idx];
                ts_B_braid[p][sample] = B_phi[0][idx];
                /* Record local energy density (S+B) */
                double rho = 0;
                for (int a = 0; a < NFIELDS; a++) {
                    rho += 0.5*g->vel[a][idx]*g->vel[a][idx]
                         + 0.5*mass2*g->phi[a][idx]*g->phi[a][idx];
                    rho += 0.5*B_vel[a][idx]*B_vel[a][idx]
                         + 0.5*mass2*B_phi[a][idx]*B_phi[a][idx];
                }
                ts_rho_braid[p][sample] = rho;
            }
            sample++;
        }

        if (step % (n_steps/10) == 0) {
            printf("  Braid: step %d/%d (%.0f%%)\n", step, n_steps, 100.0*step/n_steps);
            fflush(stdout);
        }
    }
    printf("  Braid run: %d samples collected\n", sample);
    int n_samples_actual = sample;

    /* Free braid grid, keep B arrays for reuse */
    grid_free(g);

    /* ============================================================
       RUN 2: CONTROL (background only, M7 with S=0)
       ============================================================ */
    printf("\n=== Run 2: Control (background only) ===\n");

    g = grid_alloc(N, L);
    /* S = zero (no braid) */
    grid_zero(g);
    /* B = same background */
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            for (int kk = 0; kk < N; kk++) {
                int idx = i*NN + j*N + kk;
                double z = -L + kk*dx;
                for (int a = 0; a < NFIELDS; a++) {
                    double ph = k_bg*z + 2.0*PI*a/3.0;
                    B_phi[a][idx] = A_bg * cos(ph);
                    B_vel[a][idx] = omega_bg * A_bg * sin(ph);
                }
            }
    memset(B_acc[0], 0, N3*sizeof(double));
    memset(B_acc[1], 0, N3*sizeof(double));
    memset(B_acc[2], 0, N3*sizeof(double));

    compute_forces_m7(g, BIMODAL);
    sample = 0;

    for (int step = 0; step < n_steps; step++) {
        verlet_step_m7(g, BIMODAL);
        apply_edge_damping_m7(g);

        if (step % record_every == 0 && sample < n_samples_actual) {
            for (int p = 0; p < N_PROBES; p++) {
                int idx = probe_i[p]*NN + jc*N + kc;
                ts_S_ctrl[p][sample] = g->phi[0][idx];
                ts_B_ctrl[p][sample] = B_phi[0][idx];
                double rho = 0;
                for (int a = 0; a < NFIELDS; a++) {
                    rho += 0.5*g->vel[a][idx]*g->vel[a][idx]
                         + 0.5*mass2*g->phi[a][idx]*g->phi[a][idx];
                    rho += 0.5*B_vel[a][idx]*B_vel[a][idx]
                         + 0.5*mass2*B_phi[a][idx]*B_phi[a][idx];
                }
                ts_rho_ctrl[p][sample] = rho;
            }
            sample++;
        }

        if (step % (n_steps/10) == 0) {
            printf("  Control: step %d/%d (%.0f%%)\n", step, n_steps, 100.0*step/n_steps);
            fflush(stdout);
        }
    }
    printf("  Control run: %d samples collected\n", sample);

    /* ============================================================
       SPECTRAL ANALYSIS
       ============================================================ */
    printf("\n=== Spectral Analysis ===\n");

    int N_freq = n_samples_actual / 2;  /* up to Nyquist */
    double *freq   = calloc(N_freq, sizeof(double));
    double *power  = calloc(N_freq, sizeof(double));
    double *power2 = calloc(N_freq, sizeof(double));

    /* Find the bin closest to omega_0 */
    double omega_0 = omega_bg;
    int bin_omega0 = (int)(omega_0 * n_samples_actual * dt_sample + 0.5);
    if (bin_omega0 >= N_freq) bin_omega0 = N_freq - 1;
    printf("omega_0 = %.4f, bin = %d\n", omega_0, bin_omega0);

    /* Output: spectral depletion at each probe radius */
    FILE *fp_spec = fopen("data/spectral_depletion.tsv", "w");
    fprintf(fp_spec, "r\tomega_0_power_braid\tomega_0_power_ctrl\tdelta_P_omega0\t"
                     "total_power_braid\ttotal_power_ctrl\tdelta_P_total\n");

    /* Also output full spectrum at a few key radii */
    FILE *fp_full[N_PROBES];
    for (int p = 0; p < N_PROBES; p++) {
        char fname[256];
        snprintf(fname, sizeof(fname), "data/spectrum_r%.0f.tsv", probe_r[p]);
        fp_full[p] = fopen(fname, "w");
        fprintf(fp_full[p], "freq\tP_B_braid\tP_B_ctrl\tdelta_P_B\tP_rho_braid\tP_rho_ctrl\tdelta_P_rho\n");
    }

    double delta_at_omega0[N_PROBES];
    double delta_total[N_PROBES];

    for (int p = 0; p < N_PROBES; p++) {
        /* FFT of B field (background component) for braid and control */
        double *P_B_braid = calloc(N_freq, sizeof(double));
        double *P_B_ctrl  = calloc(N_freq, sizeof(double));
        double *P_rho_braid = calloc(N_freq, sizeof(double));
        double *P_rho_ctrl  = calloc(N_freq, sizeof(double));

        dft_power(ts_B_braid[p], n_samples_actual, dt_sample, freq, P_B_braid, N_freq);
        dft_power(ts_B_ctrl[p],  n_samples_actual, dt_sample, freq, P_B_ctrl,  N_freq);
        dft_power(ts_rho_braid[p], n_samples_actual, dt_sample, freq, P_rho_braid, N_freq);
        dft_power(ts_rho_ctrl[p],  n_samples_actual, dt_sample, freq, P_rho_ctrl,  N_freq);

        /* Power at omega_0 (use a narrow band: bin_omega0 ± 2) */
        double P_w0_braid = 0, P_w0_ctrl = 0;
        double P_tot_braid = 0, P_tot_ctrl = 0;
        int bw = 3;  /* bandwidth */
        for (int f = 0; f < N_freq; f++) {
            P_tot_braid += P_B_braid[f];
            P_tot_ctrl  += P_B_ctrl[f];
            if (abs(f - bin_omega0) <= bw) {
                P_w0_braid += P_B_braid[f];
                P_w0_ctrl  += P_B_ctrl[f];
            }
        }

        delta_at_omega0[p] = P_w0_braid - P_w0_ctrl;
        delta_total[p] = P_tot_braid - P_tot_ctrl;

        fprintf(fp_spec, "%.1f\t%.8e\t%.8e\t%+.8e\t%.8e\t%.8e\t%+.8e\n",
                probe_r[p], P_w0_braid, P_w0_ctrl, delta_at_omega0[p],
                P_tot_braid, P_tot_ctrl, delta_total[p]);

        /* Full spectrum output */
        for (int f = 0; f < N_freq && f < 200; f++) {
            fprintf(fp_full[p], "%.6f\t%.8e\t%.8e\t%+.8e\t%.8e\t%.8e\t%+.8e\n",
                    freq[f], P_B_braid[f], P_B_ctrl[f],
                    P_B_braid[f] - P_B_ctrl[f],
                    P_rho_braid[f], P_rho_ctrl[f],
                    P_rho_braid[f] - P_rho_ctrl[f]);
        }

        free(P_B_braid); free(P_B_ctrl);
        free(P_rho_braid); free(P_rho_ctrl);
    }

    fclose(fp_spec);
    for (int p = 0; p < N_PROBES; p++) fclose(fp_full[p]);

    /* ============================================================
       FIT: delta_P(omega_0) vs r — power law
       ============================================================ */
    printf("\n=== Depletion at omega_0 vs radius ===\n");
    printf("  r      delta_P(w0)     delta_P(total)   ratio\n");

    for (int p = 0; p < N_PROBES; p++) {
        double ratio = (fabs(delta_total[p]) > 1e-30) ?
                       delta_at_omega0[p] / delta_total[p] : 0;
        printf("  %4.0f   %+12.6e   %+12.6e   %.3f\n",
               probe_r[p], delta_at_omega0[p], delta_total[p], ratio);
    }

    /* Power law fit on |delta_P(omega_0)| for r >= 4 */
    printf("\n=== Power-Law Fit: |delta_P(omega_0)| ~ A/r^alpha ===\n");
    {
        double sum_lnr = 0, sum_lnP = 0, sum_lr2 = 0, sum_lrlP = 0;
        int nfit = 0;
        for (int p = 0; p < N_PROBES; p++) {
            if (probe_r[p] < 4.0) continue;
            double absP = fabs(delta_at_omega0[p]);
            if (absP < 1e-30) continue;
            double lr = log(probe_r[p]), lP = log(absP);
            sum_lnr += lr; sum_lnP += lP;
            sum_lr2 += lr*lr; sum_lrlP += lr*lP;
            nfit++;
        }
        if (nfit >= 2) {
            double alpha = -(nfit*sum_lrlP - sum_lnr*sum_lnP) /
                            (nfit*sum_lr2 - sum_lnr*sum_lnr);
            double lnA = (sum_lnP + alpha*sum_lnr) / nfit;
            printf("  alpha_omega0 = %.4f (fit from %d points)\n", alpha, nfit);
            printf("  A = %.6e\n", exp(lnA));

            FILE *fp_fit = fopen("data/spectral_fit.tsv", "w");
            fprintf(fp_fit, "# Frequency-resolved depletion fit\n");
            fprintf(fp_fit, "quantity\talpha\tA\tnpoints\n");
            fprintf(fp_fit, "delta_P_omega0\t%.6f\t%.8e\t%d\n", alpha, exp(lnA), nfit);

            /* Also fit total */
            sum_lnr = sum_lnP = sum_lr2 = sum_lrlP = 0; nfit = 0;
            for (int p = 0; p < N_PROBES; p++) {
                if (probe_r[p] < 4.0) continue;
                double absP = fabs(delta_total[p]);
                if (absP < 1e-30) continue;
                double lr = log(probe_r[p]), lP = log(absP);
                sum_lnr += lr; sum_lnP += lP;
                sum_lr2 += lr*lr; sum_lrlP += lr*lP;
                nfit++;
            }
            if (nfit >= 2) {
                double alpha_tot = -(nfit*sum_lrlP - sum_lnr*sum_lnP) /
                                    (nfit*sum_lr2 - sum_lnr*sum_lnr);
                double lnA_tot = (sum_lnP + alpha_tot*sum_lnr) / nfit;
                printf("  alpha_total  = %.4f\n", alpha_tot);
                fprintf(fp_fit, "delta_P_total\t%.6f\t%.8e\t%d\n", alpha_tot, exp(lnA_tot), nfit);
            }
            fclose(fp_fit);
        }
    }

    printf("\n=== KEY QUESTION ===\n");
    printf("If alpha_omega0 < alpha_total: frequency-selective depletion\n");
    printf("   is shallower (closer to 1/r) than total depletion.\n");
    printf("   This confirms the harmonic absorption hypothesis.\n");
    printf("If alpha_omega0 >= alpha_total: no frequency selectivity.\n");
    printf("   The depletion is broadband, not resonant.\n");

    /* Cleanup */
    for (int p = 0; p < N_PROBES; p++) {
        free(ts_S_braid[p]); free(ts_B_braid[p]);
        free(ts_S_ctrl[p]);  free(ts_B_ctrl[p]);
        free(ts_rho_braid[p]); free(ts_rho_ctrl[p]);
    }
    for (int a = 0; a < NFIELDS; a++) {
        free(B_phi[a]); free(B_vel[a]); free(B_acc[a]);
    }
    free(freq); free(power); free(power2);
    grid_free(g);

    printf("\n=== T12 Spectral Analysis Complete ===\n");
    return 0;
}
