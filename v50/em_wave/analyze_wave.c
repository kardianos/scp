/*  analyze_wave.c — Track θ_x wave packet propagation along z-axis
 *
 *  Reads the EM wave SFA and for each frame extracts:
 *  - θ_x centroid position along z (weighted by θ_x²)
 *  - θ_x peak position and amplitude
 *  - θ_x RMS amplitude
 *  - φ_x RMS (to see coupling/excitation of phi sector)
 *
 *  Uses the central y-z slice (i=N/2, j=N/2) for clean 1D profile.
 *
 *  Build: gcc -O3 -o analyze_wave analyze_wave.c -lzstd -lm
 *  Usage: ./analyze_wave em_wave.sfa
 */

#define SFA_IMPLEMENTATION
#include "../../sfa/format/sfa.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

static double f16_to_f64(uint16_t h) {
    uint16_t sign = h & 0x8000;
    int exp = (h >> 10) & 0x1F;
    uint16_t mant = h & 0x3FF;
    if (exp == 0) return sign ? -0.0 : 0.0;
    if (exp == 31) return sign ? -INFINITY : INFINITY;
    float f; uint32_t x = ((uint32_t)sign << 16) | ((uint32_t)(exp-15+127) << 23) | ((uint32_t)mant << 13);
    memcpy(&f, &x, 4);
    return (double)f;
}

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "Usage: %s input.sfa\n", argv[0]); return 1; }

    SFA *sfa = sfa_open(argv[1]);
    if (!sfa) { fprintf(stderr, "Cannot open: %s\n", argv[1]); return 1; }

    int N = sfa->Nx;
    double L = sfa->Lx;
    double dx = 2.0 * L / (N - 1);
    long N3 = (long)N * N * N;
    int NN = N * N;

    printf("# Grid: %d^3, L=%.1f, dx=%.4f, %u frames\n", N, L, dx, sfa->total_frames);
    printf("# t\tz_centroid\tz_peak\ttheta_x_peak\ttheta_x_rms\tphi_rms_perturbation\tspeed_est\n");

    void *buf = malloc(sfa->frame_bytes);
    double *theta_x = (double*)malloc(N3 * sizeof(double));
    double *phi_x = (double*)malloc(N3 * sizeof(double));

    /* Find column indices */
    int col_theta_x = -1, col_phi_x = -1;
    for (int c = 0; c < sfa->n_columns; c++) {
        if (sfa->columns[c].semantic == SFA_ANGLE && sfa->columns[c].component == 0)
            col_theta_x = c;
        if (sfa->columns[c].semantic == SFA_POSITION && sfa->columns[c].component == 0)
            col_phi_x = c;
    }
    if (col_theta_x < 0) { fprintf(stderr, "No theta_x column\n"); return 1; }

    double z_prev = 0;
    double t_prev = 0;

    for (uint32_t f = 0; f < sfa->total_frames; f++) {
        SFA_L2Entry loc;
        sfa_find_frame(sfa, f, &loc);
        sfa_read_frame(sfa, f, buf);

        /* Extract theta_x and phi_x */
        uint64_t off = 0;
        for (int c = 0; c < sfa->n_columns; c++) {
            int dtype = sfa->columns[c].dtype;
            int es = sfa_dtype_size[dtype];
            uint8_t *src = (uint8_t*)buf + off;

            double *target = NULL;
            if (c == col_theta_x) target = theta_x;
            else if (c == col_phi_x) target = phi_x;

            if (target) {
                if (dtype == SFA_F64) for(long i=0;i<N3;i++) target[i]=((double*)src)[i];
                else if (dtype == SFA_F32) for(long i=0;i<N3;i++) target[i]=(double)((float*)src)[i];
                else if (dtype == SFA_F16) for(long i=0;i<N3;i++) target[i]=f16_to_f64(((uint16_t*)src)[i]);
            }
            off += (uint64_t)N3 * es;
        }

        /* Compute z-centroid of θ_x² (using full 3D grid) */
        double sum_tx2 = 0, sum_tx2_z = 0;
        double peak_tx = 0, z_peak = 0;
        double sum_tx2_all = 0;

        for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
        for (int k = 0; k < N; k++) {
            long idx = (long)i*NN + j*N + k;
            double z = -L + k * dx;
            double tx = theta_x[idx];
            double tx2 = tx * tx;
            sum_tx2_all += tx2;

            /* Use central column (i=N/2, j=N/2) for 1D profile */
            if (i == N/2 && j == N/2) {
                sum_tx2 += tx2;
                sum_tx2_z += tx2 * z;
                if (fabs(tx) > fabs(peak_tx)) {
                    peak_tx = tx;
                    z_peak = z;
                }
            }
        }}}

        double z_centroid = (sum_tx2 > 0) ? sum_tx2_z / sum_tx2 : 0;
        double theta_x_rms = sqrt(sum_tx2_all / N3);

        /* Phi perturbation: subtract mean background amplitude */
        double sum_px2 = 0;
        double A_bg = 0.1;  /* expected background */
        for (long i = 0; i < N3; i++) {
            double dp = phi_x[i];
            sum_px2 += dp * dp;
        }
        double phi_rms = sqrt(sum_px2 / N3);
        /* Subtract expected bg contribution: A_bg²/2 */
        double phi_perturbation = sqrt(fabs(phi_rms*phi_rms - A_bg*A_bg/2.0));

        double speed = 0;
        if (f > 0 && loc.time > t_prev + 0.01) {
            speed = (z_centroid - z_prev) / (loc.time - t_prev);
        }

        printf("%.3f\t%.4f\t%.4f\t%.6f\t%.6e\t%.6e\t%.4f\n",
               loc.time, z_centroid, z_peak, peak_tx, theta_x_rms, phi_perturbation, speed);

        z_prev = z_centroid;
        t_prev = loc.time;
    }

    free(buf); free(theta_x); free(phi_x);
    sfa_close(sfa);
    return 0;
}
