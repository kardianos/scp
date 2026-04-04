/*  sfa_temporal_vec.c — Multi-frame temporal vectorization
 *
 *  Given an SFA file with multiple frames, fits 4D polynomial patches
 *  (x, y, z, t) that capture both spatial and temporal variation.
 *
 *  Approach:
 *  1. For each spatial block (8³ voxels), extract the time series
 *     across all frames → NxNxNxT data
 *  2. Fit a temporal polynomial (order 1-3) to each coefficient's
 *     time evolution
 *  3. Compare: per-frame vectorization (N_frames × N_patches × 64 coeffs)
 *     vs temporal (N_patches × 64 coeffs × (temporal_order+1))
 *
 *  For periodic signals (breathing), a Fourier fit (A*cos(ωt+φ))
 *  may be more efficient than polynomial.
 *
 *  Build: gcc -O3 -march=native -fopenmp -o sfa_temporal_vec \
 *         sfa_temporal_vec.c -lzstd -lm
 */

#define SFA_IMPLEMENTATION
#include "../format/sfa.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define PI 3.14159265358979323846

static float f16f(uint16_t h) {
    uint16_t s = h & 0x8000; int e = (h >> 10) & 0x1F; uint16_t m = h & 0x3FF;
    if (e == 0) return 0; if (e == 31) return s ? -1e30f : 1e30f;
    float f; uint32_t x = ((uint32_t)s << 16) | ((uint32_t)(e-15+127) << 23) | ((uint32_t)m << 13);
    memcpy(&f, &x, 4); return f;
}

static float col_f(void *buf, SFA *sfa, int c, long i) {
    long N3 = (long)sfa->Nx * sfa->Ny * sfa->Nz;
    uint64_t off = 0;
    for (int cc = 0; cc < c; cc++)
        off += (uint64_t)N3 * sfa_dtype_size[sfa->columns[cc].dtype];
    int dt = sfa->columns[c].dtype;
    uint8_t *src = (uint8_t*)buf + off;
    if (dt == SFA_F16) return f16f(((uint16_t*)src)[i]);
    if (dt == SFA_F32) return ((float*)src)[i];
    return (float)((double*)src)[i];
}

/* Fit a polynomial of given order to time series data.
 * Returns coefficients c[0..order] for c0 + c1*t + c2*t² + ... */
static void fit_temporal_poly(const double *times, const double *values, int n,
                               int order, double *coeffs, double *rms_err) {
    /* Least squares normal equations */
    int p = order + 1;
    double S[16] = {0}, Sy[8] = {0}; /* max order 3 → 4 params */

    for (int i = 0; i < n; i++) {
        double t = times[i];
        double y = values[i];
        double tk = 1;
        for (int k = 0; k < 2*p-1; k++) { S[k] += tk; tk *= t; }
        tk = 1;
        for (int k = 0; k < p; k++) { Sy[k] += tk * y; tk *= t; }
    }

    /* Solve via Gaussian elimination */
    double a[4][5];
    for (int r = 0; r < p; r++) {
        for (int c = 0; c < p; c++) a[r][c] = S[r+c];
        a[r][p] = Sy[r];
    }

    for (int col = 0; col < p; col++) {
        int piv = col;
        for (int row = col+1; row < p; row++)
            if (fabs(a[row][col]) > fabs(a[piv][col])) piv = row;
        if (piv != col) {
            for (int j = 0; j <= p; j++) {
                double tmp = a[col][j]; a[col][j] = a[piv][j]; a[piv][j] = tmp;
            }
        }
        if (fabs(a[col][col]) < 1e-15) continue;
        for (int row = col+1; row < p; row++) {
            double f = a[row][col] / a[col][col];
            for (int j = col; j <= p; j++) a[row][j] -= f * a[col][j];
        }
    }
    for (int row = p-1; row >= 0; row--) {
        double sum = a[row][p];
        for (int j = row+1; j < p; j++) sum -= a[row][j] * coeffs[j];
        coeffs[row] = (fabs(a[row][row]) > 1e-15) ? sum / a[row][row] : 0;
    }

    /* Compute RMS error */
    double se = 0;
    for (int i = 0; i < n; i++) {
        double t = times[i];
        double pred = 0, tk = 1;
        for (int k = 0; k < p; k++) { pred += coeffs[k] * tk; tk *= t; }
        double err = values[i] - pred;
        se += err * err;
    }
    *rms_err = sqrt(se / n);
}

/* Fit a sinusoidal model: A*cos(ω*t + φ) + offset
 * Uses DFT to find dominant frequency, then least-squares for amplitude/phase. */
static void fit_temporal_fourier(const double *times, const double *values, int n,
                                  double *amp, double *omega, double *phase, double *offset,
                                  double *rms_err) {
    /* Mean (offset) */
    double mean = 0;
    for (int i = 0; i < n; i++) mean += values[i];
    mean /= n;
    *offset = mean;

    /* DFT to find dominant frequency */
    double dt_avg = (n > 1) ? (times[n-1] - times[0]) / (n-1) : 1;
    double best_power = 0, best_omega = 0;
    int n_freq = n / 2;
    for (int f = 1; f <= n_freq; f++) {
        double w = 2 * PI * f / (n * dt_avg);
        double re = 0, im = 0;
        for (int i = 0; i < n; i++) {
            double v = values[i] - mean;
            re += v * cos(w * times[i]);
            im += v * sin(w * times[i]);
        }
        double power = re*re + im*im;
        if (power > best_power) {
            best_power = power;
            best_omega = w;
        }
    }
    *omega = best_omega;

    /* Least squares: values[i] ≈ A*cos(ω*t_i) + B*sin(ω*t_i) + C */
    double Scc = 0, Sss = 0, Scs = 0, Scy = 0, Ssy = 0;
    for (int i = 0; i < n; i++) {
        double c = cos(best_omega * times[i]);
        double s = sin(best_omega * times[i]);
        double y = values[i] - mean;
        Scc += c*c; Sss += s*s; Scs += c*s;
        Scy += c*y; Ssy += s*y;
    }

    double det = Scc*Sss - Scs*Scs;
    double A_cos = 0, B_sin = 0;
    if (fabs(det) > 1e-15) {
        A_cos = (Sss*Scy - Scs*Ssy) / det;
        B_sin = (Scc*Ssy - Scs*Scy) / det;
    }

    *amp = sqrt(A_cos*A_cos + B_sin*B_sin);
    *phase = atan2(-B_sin, A_cos);

    /* RMS error */
    double se = 0;
    for (int i = 0; i < n; i++) {
        double pred = mean + (*amp) * cos(best_omega * times[i] + (*phase));
        double err = values[i] - pred;
        se += err * err;
    }
    *rms_err = sqrt(se / n);
}

int main(int argc, char **argv) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s input.sfa [-col N] [-block_size B]\n", argv[0]);
        fprintf(stderr, "\nAnalyzes temporal coherence across SFA frames.\n");
        fprintf(stderr, "Reports compression ratios for polynomial and Fourier fits.\n");
        return 1;
    }

    int col = 0;      /* field column to analyze */
    int BS = 8;        /* block size */
    for (int i = 2; i < argc; i++) {
        if (!strcmp(argv[i], "-col") && i+1 < argc) col = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-block_size") && i+1 < argc) BS = atoi(argv[++i]);
    }

    SFA *sfa = sfa_open(argv[1]);
    if (!sfa) { fprintf(stderr, "Cannot open %s\n", argv[1]); return 1; }

    int N = sfa->Nx;
    long N3 = (long)N * N * N;
    int nf = sfa->total_frames;

    printf("sfa_temporal_vec: %s (N=%d, %d frames, col=%d)\n", argv[1], N, nf, col);

    if (nf < 3) {
        fprintf(stderr, "Need at least 3 frames for temporal analysis\n");
        sfa_close(sfa);
        return 1;
    }

    /* Load all frames for the target column */
    printf("Loading %d frames...\n", nf);
    float **frame_data = malloc(nf * sizeof(float*));
    double *frame_times = malloc(nf * sizeof(double));
    void *buf = malloc(sfa->frame_bytes);

    for (int fi = 0; fi < nf; fi++) {
        frame_data[fi] = malloc(N3 * sizeof(float));
        double t = sfa_frame_time(sfa, fi);
        frame_times[fi] = t;

        if (sfa_read_frame(sfa, fi, buf) != 0) {
            fprintf(stderr, "Cannot read frame %d\n", fi);
            return 1;
        }
        for (long i = 0; i < N3; i++)
            frame_data[fi][i] = col_f(buf, sfa, col, i);

        printf("  Frame %d: t=%.3f\n", fi, t);
    }
    free(buf);

    /* Normalize times to [0, 1] for numerical stability */
    double t0 = frame_times[0], tmax = frame_times[nf-1];
    double *tnorm = malloc(nf * sizeof(double));
    for (int fi = 0; fi < nf; fi++)
        tnorm[fi] = (nf > 1) ? (frame_times[fi] - t0) / (tmax - t0) : 0;

    /* Analyze temporal coherence at selected probe points */
    int n_blocks = (N / BS) * (N / BS) * (N / BS);
    int BN = N / BS;

    printf("\nAnalyzing %d blocks (%d³ × %d³)...\n", n_blocks, BN, BS);

    /* For each block center, extract the time series and fit */
    double poly1_total_rms = 0, poly2_total_rms = 0, poly3_total_rms = 0;
    double fourier_total_rms = 0;
    double raw_total_var = 0;
    int n_analyzed = 0;

    for (int bi = 0; bi < BN; bi++)
    for (int bj = 0; bj < BN; bj++)
    for (int bk = 0; bk < BN; bk++) {
        /* Block center voxel */
        int ci = bi * BS + BS/2;
        int cj = bj * BS + BS/2;
        int ck = bk * BS + BS/2;
        long idx = (long)ci * N * N + cj * N + ck;

        /* Extract time series at this voxel */
        double *ts = malloc(nf * sizeof(double));
        for (int fi = 0; fi < nf; fi++)
            ts[fi] = frame_data[fi][idx];

        /* Variance of time series */
        double mean = 0;
        for (int fi = 0; fi < nf; fi++) mean += ts[fi];
        mean /= nf;
        double var = 0;
        for (int fi = 0; fi < nf; fi++) { double d = ts[fi]-mean; var += d*d; }
        var /= nf;
        raw_total_var += var;

        /* Polynomial fits */
        double c1[2], c2[3], c3[4];
        double rms1, rms2, rms3;
        fit_temporal_poly(tnorm, ts, nf, 1, c1, &rms1);
        fit_temporal_poly(tnorm, ts, nf, 2, c2, &rms2);
        fit_temporal_poly(tnorm, ts, nf, 3, c3, &rms3);

        poly1_total_rms += rms1;
        poly2_total_rms += rms2;
        poly3_total_rms += rms3;

        /* Fourier fit */
        double amp, omega, phase, offset, rms_f;
        fit_temporal_fourier(frame_times, ts, nf, &amp, &omega, &phase, &offset, &rms_f);
        fourier_total_rms += rms_f;

        n_analyzed++;
        free(ts);
    }

    double avg_var = raw_total_var / n_analyzed;
    double avg_std = sqrt(avg_var);

    printf("\n=== Temporal Coherence Report ===\n");
    printf("  Frames: %d, time range: %.3f to %.3f\n", nf, t0, tmax);
    printf("  Blocks analyzed: %d (center voxel of each %d³ block)\n", n_analyzed, BS);
    printf("  Avg temporal std dev: %.6e\n", avg_std);
    printf("\n");

    printf("  %-20s  %12s  %12s  %8s\n", "Method", "Avg RMS", "Rel to std", "Params/block");
    printf("  %-20s  %12s  %12s  %8s\n", "------", "-------", "----------", "-----------");

    double r1 = poly1_total_rms / n_analyzed;
    double r2 = poly2_total_rms / n_analyzed;
    double r3 = poly3_total_rms / n_analyzed;
    double rf = fourier_total_rms / n_analyzed;

    printf("  %-20s  %12.6e  %12.2f%%  %8d\n", "Per-frame (raw)", avg_std, 100.0, nf);
    printf("  %-20s  %12.6e  %12.2f%%  %8d\n", "Linear (t)", r1, (avg_std>0?r1/avg_std*100:0), 2);
    printf("  %-20s  %12.6e  %12.2f%%  %8d\n", "Quadratic (t²)", r2, (avg_std>0?r2/avg_std*100:0), 3);
    printf("  %-20s  %12.6e  %12.2f%%  %8d\n", "Cubic (t³)", r3, (avg_std>0?r3/avg_std*100:0), 4);
    printf("  %-20s  %12.6e  %12.2f%%  %8d\n", "Fourier (cos(ωt+φ))", rf, (avg_std>0?rf/avg_std*100:0), 4);

    printf("\n=== Compression Analysis ===\n");
    int spatial_coeffs = 64; /* tricubic per block */
    printf("  Per-frame spatial: %d blocks × %d coeffs × %d frames = %d total\n",
           n_blocks, spatial_coeffs, nf, n_blocks * spatial_coeffs * nf);

    printf("  Temporal poly(1): %d blocks × %d spatial × 2 temporal = %d total (%.1f× vs per-frame)\n",
           n_blocks, spatial_coeffs, n_blocks * spatial_coeffs * 2,
           (double)(n_blocks * spatial_coeffs * nf) / (n_blocks * spatial_coeffs * 2));

    printf("  Temporal poly(3): %d blocks × %d spatial × 4 temporal = %d total (%.1f× vs per-frame)\n",
           n_blocks, spatial_coeffs, n_blocks * spatial_coeffs * 4,
           (double)(n_blocks * spatial_coeffs * nf) / (n_blocks * spatial_coeffs * 4));

    printf("  Temporal Fourier: %d blocks × %d spatial × 4 temporal = %d total (%.1f× vs per-frame)\n",
           n_blocks, spatial_coeffs, n_blocks * spatial_coeffs * 4,
           (double)(n_blocks * spatial_coeffs * nf) / (n_blocks * spatial_coeffs * 4));

    /* Cleanup */
    for (int fi = 0; fi < nf; fi++) free(frame_data[fi]);
    free(frame_data);
    free(frame_times);
    free(tnorm);
    sfa_close(sfa);
    return 0;
}
