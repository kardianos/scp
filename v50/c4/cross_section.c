/* cross_section.c — Radial profile analysis of SFA frame
 *
 * Finds the particle centroid (|P|-weighted), then computes radial
 * shell profiles of key quantities: |phi|, |P|, |theta|, |curl(phi)|,
 * energy density, and the hardening term |theta|^2 * |curl(phi)|^2.
 *
 * Build: gcc -O3 -o cross_section cross_section.c -lzstd -lm
 * Usage: ./cross_section input.sfa [frame_idx]
 */

#define SFA_IMPLEMENTATION
#include "../../sfa/format/sfa.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define NBINS 40

int main(int argc, char **argv) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s input.sfa [frame_idx]\n", argv[0]);
        return 1;
    }

    SFA *sfa = sfa_open(argv[1]);
    if (!sfa) { fprintf(stderr, "Cannot open %s\n", argv[1]); return 1; }

    int frame = (argc > 2) ? atoi(argv[2]) : sfa->total_frames - 1;
    int N = sfa->Nx;
    double L = sfa->Lx;
    double dx = 2.0 * L / (N - 1);
    long N3 = (long)N * N * N;
    int NN = N * N;

    printf("File: %s  N=%d L=%.1f dx=%.4f  Frame %d/%d\n", argv[1], N, L, dx, frame, sfa->total_frames);

    /* Read frame */
    size_t frame_bytes = 0;
    for (uint32_t c = 0; c < sfa->n_columns; c++)
        frame_bytes += N3 * sfa_dtype_size[sfa->columns[c].dtype];
    void *buf = malloc(frame_bytes);
    if (sfa_read_frame(sfa, frame, buf) < 0) {
        fprintf(stderr, "Cannot read frame %d\n", frame);
        return 1;
    }

    double t = sfa_frame_time(sfa, frame);
    printf("Time: %.2f\n\n", t);

    /* Extract columns as float arrays */
    float *phi[3], *theta[3];
    size_t col_bytes = N3 * sizeof(float);
    for (int a = 0; a < 3; a++) {
        phi[a] = (float*)((char*)buf + a * col_bytes);
        theta[a] = (float*)((char*)buf + (3 + a) * col_bytes);
    }

    /* Pass 1: find centroid weighted by |P| */
    double cx = 0, cy = 0, cz = 0, wsum = 0;
    for (long idx = 0; idx < N3; idx++) {
        int i = idx / NN, j = (idx / N) % N, k = idx % N;
        double P = fabs((double)phi[0][idx] * phi[1][idx] * phi[2][idx]);
        double x = -L + i * dx, y = -L + j * dx, z = -L + k * dx;
        cx += P * x; cy += P * y; cz += P * z; wsum += P;
    }
    if (wsum > 0) { cx /= wsum; cy /= wsum; cz /= wsum; }
    printf("Centroid (|P|-weighted): (%.2f, %.2f, %.2f)\n\n", cx, cy, cz);

    /* Pass 2: radial shell profiles */
    double R_max = L * 0.9;  /* don't go to boundary */
    double dr = R_max / NBINS;
    double idx1 = 1.0 / (2.0 * dx);

    double bin_phi_rms[NBINS] = {0};
    double bin_P_abs[NBINS] = {0};
    double bin_theta_rms[NBINS] = {0};
    double bin_curl_rms[NBINS] = {0};
    double bin_harden[NBINS] = {0};
    double bin_mismatch[NBINS] = {0};
    int bin_count[NBINS] = {0};

    for (long idx = 0; idx < N3; idx++) {
        int i = idx / NN, j = (idx / N) % N, k = idx % N;
        double x = -L + i * dx - cx;
        double y = -L + j * dx - cy;
        double z = -L + k * dx - cz;
        double r = sqrt(x*x + y*y + z*z);
        int bin = (int)(r / dr);
        if (bin >= NBINS) continue;

        double p0 = phi[0][idx], p1 = phi[1][idx], p2 = phi[2][idx];
        double t0 = theta[0][idx], t1 = theta[1][idx], t2 = theta[2][idx];

        double phi_sq = p0*p0 + p1*p1 + p2*p2;
        double P = fabs(p0 * p1 * p2);
        double theta_sq = t0*t0 + t1*t1 + t2*t2;

        /* curl(phi) — need neighbors */
        int ip=(i+1)%N, im=(i-1+N)%N, jp=(j+1)%N, jm=(j-1+N)%N, kp=(k+1)%N, km=(k-1+N)%N;
        long n_jp=(long)i*NN+jp*N+k, n_jm=(long)i*NN+jm*N+k;
        long n_kp=(long)i*NN+j*N+kp, n_km=(long)i*NN+j*N+km;
        long n_ip=(long)ip*NN+j*N+k, n_im=(long)im*NN+j*N+k;

        double c0 = (phi[2][n_jp]-phi[2][n_jm]-phi[1][n_kp]+phi[1][n_km])*idx1;
        double c1 = (phi[0][n_kp]-phi[0][n_km]-phi[2][n_ip]+phi[2][n_im])*idx1;
        double c2 = (phi[1][n_ip]-phi[1][n_im]-phi[0][n_jp]+phi[0][n_jm])*idx1;
        double curl_sq = c0*c0 + c1*c1 + c2*c2;

        /* Cosserat mismatch |curl(phi)/2 - theta|² */
        double m0 = c0*0.5 - t0, m1 = c1*0.5 - t1, m2 = c2*0.5 - t2;
        double mismatch_sq = m0*m0 + m1*m1 + m2*m2;

        /* Hardening energy density: (β/2)|θ|²|∇×φ|² */
        double harden = theta_sq * curl_sq;

        bin_phi_rms[bin] += phi_sq;
        bin_P_abs[bin] += P;
        bin_theta_rms[bin] += theta_sq;
        bin_curl_rms[bin] += curl_sq;
        bin_harden[bin] += harden;
        bin_mismatch[bin] += mismatch_sq;
        bin_count[bin]++;
    }

    /* Print profiles */
    printf("%6s %8s %8s %10s %10s %10s %10s %8s\n",
           "r", "phi_rms", "|P|", "theta_rms", "curl_rms", "harden", "mismatch", "count");
    printf("------+--------+--------+----------+----------+----------+----------+--------\n");

    for (int b = 0; b < NBINS; b++) {
        if (bin_count[b] == 0) continue;
        double r = (b + 0.5) * dr;
        double n = bin_count[b];
        printf("%6.2f %8.4f %8.5f %10.5f %10.4f %10.6f %10.6f %8d\n",
               r,
               sqrt(bin_phi_rms[b] / n),
               bin_P_abs[b] / n,
               sqrt(bin_theta_rms[b] / n),
               sqrt(bin_curl_rms[b] / n),
               bin_harden[b] / n,
               bin_mismatch[b] / n,
               bin_count[b]);
    }

    /* Find shell: where is hardening maximum? */
    int max_bin = 0;
    double max_harden = 0;
    for (int b = 0; b < NBINS; b++) {
        if (bin_count[b] == 0) continue;
        double h = bin_harden[b] / bin_count[b];
        if (h > max_harden) { max_harden = h; max_bin = b; }
    }
    printf("\nPeak hardening at r=%.2f (bin %d), value=%.6f\n",
           (max_bin + 0.5) * dr, max_bin, max_harden);

    /* Find particle boundary: where |P| drops below 10% of peak */
    double P_peak = 0;
    for (int b = 0; b < NBINS; b++) {
        if (bin_count[b] == 0) continue;
        double p = bin_P_abs[b] / bin_count[b];
        if (p > P_peak) P_peak = p;
    }
    for (int b = 0; b < NBINS; b++) {
        if (bin_count[b] == 0) continue;
        double p = bin_P_abs[b] / bin_count[b];
        if (p < 0.1 * P_peak && b > 0) {
            printf("Particle boundary (|P| < 10%% peak) at r=%.2f\n", (b + 0.5) * dr);
            break;
        }
    }

    free(buf);
    sfa_close(sfa);
    return 0;
}
