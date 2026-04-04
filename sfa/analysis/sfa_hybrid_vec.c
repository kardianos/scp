/*  sfa_hybrid_vec.c — Hybrid temporal+spatial vectorization
 *
 *  Two-layer compression:
 *    Layer 1: Fourier temporal base — A*cos(ωt+φ) + offset per spatial block
 *             Captures the breathing mode (~85% of variation)
 *    Layer 2: Residual correction — short polynomial vectors on the
 *             per-frame residual (actual - Fourier prediction)
 *             The residual is small and smooth → compresses aggressively
 *
 *  Total params = Fourier(4 per spatial coeff) + Residual(adaptive, sparse)
 *
 *  Build: gcc -O3 -march=native -fopenmp -o sfa_hybrid_vec sfa_hybrid_vec.c -lzstd -lm
 */

#define SFA_IMPLEMENTATION
#include "../format/sfa.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define PI 3.14159265358979323846
#define MAX_FRAMES 256

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

/* ===== Fourier fit per voxel across time ===== */

typedef struct {
    float amp, omega, phase, offset;
} FourierCoeffs;

static void fit_fourier_voxel(const float **frames, int nf, const double *times,
                               long idx, FourierCoeffs *fc) {
    /* Mean */
    double mean = 0;
    for (int fi = 0; fi < nf; fi++) mean += frames[fi][idx];
    mean /= nf;
    fc->offset = (float)mean;

    /* DFT for dominant frequency */
    double dt_avg = (nf > 1) ? (times[nf-1] - times[0]) / (nf-1) : 1;
    double best_power = 0, best_omega = 0;
    for (int f = 1; f <= nf/2; f++) {
        double w = 2 * PI * f / (nf * dt_avg);
        double re = 0, im = 0;
        for (int fi = 0; fi < nf; fi++) {
            double v = frames[fi][idx] - mean;
            re += v * cos(w * times[fi]);
            im += v * sin(w * times[fi]);
        }
        double power = re*re + im*im;
        if (power > best_power) { best_power = power; best_omega = w; }
    }
    fc->omega = (float)best_omega;

    /* Least squares: A*cos(ωt) + B*sin(ωt) */
    double Scc=0, Sss=0, Scs=0, Scy=0, Ssy=0;
    for (int fi = 0; fi < nf; fi++) {
        double c = cos(best_omega * times[fi]);
        double s = sin(best_omega * times[fi]);
        double y = frames[fi][idx] - mean;
        Scc += c*c; Sss += s*s; Scs += c*s; Scy += c*y; Ssy += s*y;
    }
    double det = Scc*Sss - Scs*Scs;
    double Ac = 0, Bs = 0;
    if (fabs(det) > 1e-15) {
        Ac = (Sss*Scy - Scs*Ssy) / det;
        Bs = (Scc*Ssy - Scs*Scy) / det;
    }
    fc->amp = (float)sqrt(Ac*Ac + Bs*Bs);
    fc->phase = (float)atan2(-Bs, Ac);
}

static float fourier_predict(const FourierCoeffs *fc, double t) {
    return fc->offset + fc->amp * cosf(fc->omega * t + fc->phase);
}

/* ===== Spatial residual vectorization (simplified 1D along z) ===== */

typedef struct {
    float c[4]; /* polynomial coefficients up to cubic */
    int order;
    int start, span;
    float max_err;
} ResidVec;

static ResidVec fit_residual_cubic(const float *data, int i0, int i1) {
    ResidVec v = {0};
    v.order = 3;
    v.start = i0;
    v.span = i1 - i0 + 1;
    int n = v.span;

    /* Least squares cubic */
    double S[7]={0}, Sy[4]={0};
    for (int i = 0; i < n; i++) {
        double t = (n>1) ? (double)i/(n-1) : 0;
        double y = data[i0+i];
        double tk = 1;
        for (int k = 0; k < 7; k++) { S[k] += tk; tk *= t; }
        tk = 1;
        for (int k = 0; k < 4; k++) { Sy[k] += tk*y; tk *= t; }
    }
    double a[4][5];
    for (int r = 0; r < 4; r++) {
        for (int c = 0; c < 4; c++) a[r][c] = S[r+c];
        a[r][4] = Sy[r];
    }
    for (int col = 0; col < 4; col++) {
        int piv = col;
        for (int row = col+1; row < 4; row++)
            if (fabs(a[row][col]) > fabs(a[piv][col])) piv = row;
        if (piv != col)
            for (int j = 0; j < 5; j++) { double tmp=a[col][j]; a[col][j]=a[piv][j]; a[piv][j]=tmp; }
        if (fabs(a[col][col]) < 1e-15) continue;
        for (int row = col+1; row < 4; row++) {
            double f = a[row][col] / a[col][col];
            for (int j = col; j < 5; j++) a[row][j] -= f*a[col][j];
        }
    }
    double dc[4] = {0};
    for (int row = 3; row >= 0; row--) {
        double sum = a[row][4];
        for (int j = row+1; j < 4; j++) sum -= a[row][j]*dc[j];
        dc[row] = (fabs(a[row][row])>1e-15) ? sum/a[row][row] : 0;
    }
    for (int k = 0; k < 4; k++) v.c[k] = (float)dc[k];

    /* Error */
    v.max_err = 0;
    for (int i = 0; i < n; i++) {
        double t = (n>1) ? (double)i/(n-1) : 0;
        double pred = dc[0] + dc[1]*t + dc[2]*t*t + dc[3]*t*t*t;
        float err = fabsf(data[i0+i] - (float)pred);
        if (err > v.max_err) v.max_err = err;
    }
    return v;
}

/* Adaptive residual vectorization along one axis:
 * Start with large spans, split where error exceeds tolerance */
static int vectorize_residual_adaptive(const float *resid, int N, float tol,
                                        ResidVec *vecs, int max_vecs) {
    int nv = 0;
    int pos = 0;

    while (pos < N && nv < max_vecs) {
        /* Try largest possible span first, shrink if error too high */
        int best_span = 2;
        for (int span = N - pos; span >= 2; span /= 2) {
            if (pos + span > N) continue;
            ResidVec v = fit_residual_cubic(resid, pos, pos + span - 1);
            if (v.max_err <= tol || span <= 2) {
                vecs[nv++] = v;
                best_span = span;
                pos += span;
                goto next;
            }
        }
        /* Fallback: 2-point segment */
        vecs[nv++] = fit_residual_cubic(resid, pos, pos + 1);
        pos += 2;
        next:;
    }
    return nv;
}

int main(int argc, char **argv) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s input.sfa [-col N] [-tol 0.001]\n", argv[0]);
        return 1;
    }

    int col = 0;
    float tol = 0.001;
    for (int i = 2; i < argc; i++) {
        if (!strcmp(argv[i], "-col") && i+1<argc) col = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-tol") && i+1<argc) tol = atof(argv[++i]);
    }

    SFA *sfa = sfa_open(argv[1]);
    if (!sfa) { fprintf(stderr, "Cannot open %s\n", argv[1]); return 1; }

    int N = sfa->Nx;
    long N3 = (long)N*N*N;
    int nf = sfa->total_frames;
    if (nf < 3 || nf > MAX_FRAMES) {
        fprintf(stderr, "Need 3-%d frames (got %d)\n", MAX_FRAMES, nf);
        return 1;
    }

    printf("sfa_hybrid_vec: %s (N=%d, %d frames, col=%d, tol=%.4f)\n",
           argv[1], N, nf, col, tol);

    /* Load all frames */
    printf("Loading frames...\n");
    float **frames = malloc(nf * sizeof(float*));
    double times[MAX_FRAMES];
    void *buf = malloc(sfa->frame_bytes);

    for (int fi = 0; fi < nf; fi++) {
        frames[fi] = malloc(N3 * sizeof(float));
        times[fi] = sfa_frame_time(sfa, fi);
        if (sfa_read_frame(sfa, fi, buf) != 0) {
            fprintf(stderr, "Cannot read frame %d\n", fi);
            return 1;
        }
        for (long i = 0; i < N3; i++)
            frames[fi][i] = col_f(buf, sfa, col, i);
    }
    free(buf);

    /* ===== Layer 1: Fourier fit per voxel ===== */
    printf("Layer 1: Fourier temporal fit (%ld voxels × 4 params)...\n", N3);
    FourierCoeffs *fourier = malloc(N3 * sizeof(FourierCoeffs));

    double fourier_rms_total = 0;
    double raw_rms_total = 0;

    #pragma omp parallel for reduction(+:fourier_rms_total,raw_rms_total) schedule(static)
    for (long i = 0; i < N3; i++) {
        fit_fourier_voxel((const float**)frames, nf, times, i, &fourier[i]);

        /* Compute Fourier RMS for this voxel across time */
        double se_f = 0, se_r = 0, mean = fourier[i].offset;
        for (int fi = 0; fi < nf; fi++) {
            float pred = fourier_predict(&fourier[i], times[fi]);
            float actual = frames[fi][i];
            se_f += (actual - pred) * (actual - pred);
            se_r += (actual - mean) * (actual - mean);
        }
        fourier_rms_total += sqrt(se_f / nf);
        raw_rms_total += sqrt(se_r / nf);
    }

    double fourier_avg_rms = fourier_rms_total / N3;
    double raw_avg_rms = raw_rms_total / N3;
    printf("  Fourier residual: avg RMS = %.6e (%.1f%% of raw variation)\n",
           fourier_avg_rms, (raw_avg_rms > 0 ? fourier_avg_rms / raw_avg_rms * 100 : 0));

    /* ===== Layer 2: Vectorize the residual per frame ===== */
    printf("Layer 2: Spatial vectorization of residual (tol=%.4f)...\n", tol);

    int mid = N / 2; /* analyze z-axis through center */
    float *resid_slice = malloc(N * sizeof(float));
    ResidVec *rvecs = malloc(N * sizeof(ResidVec));

    int total_resid_vecs = 0;
    double resid_after_rms_total = 0;
    int total_raw_vecs_needed = 0;

    printf("\n  %5s  %8s  %8s  %8s  %8s  %10s\n",
           "Frame", "Resid_RMS", "Vecs", "Raw_Vecs", "Savings", "FinalRMS");

    for (int fi = 0; fi < nf; fi++) {
        /* Extract residual along z-axis at center */
        double resid_rms = 0;
        for (int k = 0; k < N; k++) {
            long idx = (long)mid * N*N + mid * N + k;
            float pred = fourier_predict(&fourier[idx], times[fi]);
            resid_slice[k] = frames[fi][idx] - pred;
            resid_rms += resid_slice[k] * resid_slice[k];
        }
        resid_rms = sqrt(resid_rms / N);

        /* Vectorize the residual adaptively */
        int nv = vectorize_residual_adaptive(resid_slice, N, tol, rvecs, N);
        total_resid_vecs += nv;

        /* Also vectorize the RAW data for comparison */
        float *raw_slice = malloc(N * sizeof(float));
        ResidVec *raw_vecs = malloc(N * sizeof(ResidVec));
        for (int k = 0; k < N; k++) {
            long idx = (long)mid * N*N + mid * N + k;
            raw_slice[k] = frames[fi][idx];
        }
        int nv_raw = vectorize_residual_adaptive(raw_slice, N, tol, raw_vecs, N);
        total_raw_vecs_needed += nv_raw;

        /* Reconstruct and measure final error */
        double final_rms = 0;
        for (int vi = 0; vi < nv; vi++) {
            ResidVec *v = &rvecs[vi];
            for (int i = 0; i < v->span; i++) {
                int k = v->start + i;
                if (k >= N) break;
                double t = (v->span > 1) ? (double)i / (v->span - 1) : 0;
                float resid_pred = v->c[0] + v->c[1]*t + v->c[2]*t*t + v->c[3]*t*t*t;

                long idx = (long)mid * N*N + mid * N + k;
                float fourier_pred = fourier_predict(&fourier[idx], times[fi]);
                float total_pred = fourier_pred + resid_pred;
                float actual = frames[fi][idx];
                float err = actual - total_pred;
                final_rms += err * err;
            }
        }
        final_rms = sqrt(final_rms / N);
        resid_after_rms_total += final_rms;

        double savings = (nv_raw > 0) ? (1.0 - (double)nv / nv_raw) * 100 : 0;
        printf("  %5d  %8.2e  %8d  %8d  %7.1f%%  %10.2e\n",
               fi, resid_rms, nv, nv_raw, savings, final_rms);

        free(raw_slice);
        free(raw_vecs);
    }

    printf("\n=== Hybrid Compression Summary ===\n");
    printf("  Layer 1 (Fourier): 4 params/voxel × %ld voxels = %ld params\n",
           N3, N3 * 4);
    printf("  Layer 2 (Residual): %d vectors × 4 coeffs = %d params (z-slice only)\n",
           total_resid_vecs, total_resid_vecs * 4);
    printf("  Raw vectorization: %d vectors × 4 coeffs = %d params (z-slice only)\n",
           total_raw_vecs_needed, total_raw_vecs_needed * 4);
    printf("  Residual vec savings: %.1f%% fewer vectors than raw\n",
           (total_raw_vecs_needed > 0 ?
            (1.0 - (double)total_resid_vecs / total_raw_vecs_needed) * 100 : 0));
    printf("  Final avg RMS: %.6e\n", resid_after_rms_total / nf);

    printf("\n=== Theoretical Full-Grid Compression ===\n");
    long raw_total = (long)N3 * nf; /* raw voxel count */
    long fourier_params = N3 * 4;   /* 4 per voxel for Fourier base */
    double resid_ratio = (total_raw_vecs_needed > 0) ?
        (double)total_resid_vecs / total_raw_vecs_needed : 1.0;
    /* For 3D: per-frame spatial = ~N3/8 patches × 64 coeffs */
    long spatial_per_frame = (N3 / 512) * 64;
    long raw_spatial_total = spatial_per_frame * nf;
    long hybrid_total = fourier_params + (long)(spatial_per_frame * nf * resid_ratio);

    printf("  Raw per-frame spatial: %ld params across %d frames\n", raw_spatial_total, nf);
    printf("  Hybrid (Fourier + residual): %ld params\n", hybrid_total);
    printf("  Compression vs raw: %.1f×\n", (double)raw_total / hybrid_total);
    printf("  Compression vs per-frame spatial: %.1f×\n", (double)raw_spatial_total / hybrid_total);

    /* Cleanup */
    for (int fi = 0; fi < nf; fi++) free(frames[fi]);
    free(frames); free(fourier); free(resid_slice); free(rvecs);
    sfa_close(sfa);
    return 0;
}
