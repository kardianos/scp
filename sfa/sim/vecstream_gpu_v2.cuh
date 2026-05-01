/*  vecstream_gpu_v2.cuh — GPU vecstream with multi-field patches + temporal streaming
 *
 *  Multi-field: 6 independent tricubic fits per 8³ block.
 *  Layout per patch: [phi_x(64) | phi_y(64) | phi_z(64) | theta_x(64) | theta_y(64) | theta_z(64)]
 *  Total: 6 × 64 = 384 coefficients per patch.
 *  Same layout as CPU kernel for reader compatibility.
 *
 *  Temporal: each patch coefficient maintains a running temporal model
 *  (mean + amplitude × cos(ωt + φ)). Per timestep, only the RESIDUAL
 *  (actual - predicted) is stored. Most residuals are near zero → sparse.
 *
 *  The temporal model is refit every N_refit frames for drift correction.
 */

#ifndef VECSTREAM_GPU_V2_CUH
#define VECSTREAM_GPU_V2_CUH

#include <cuda_runtime.h>
#include <stdint.h>
#include <math.h>

#define VS2_BS 8
#define VS2_NC 4
#define VS2_NCOEFFS 64  /* per scalar field per patch */
#define VS2_NFIELDS 6
#define VS2_MULTI_TOTAL (VS2_NFIELDS * VS2_NCOEFFS)  /* 6 × 64 = 384 */

/* Temporal model per coefficient: mean + amp*cos(omega*t + phase) */
typedef struct {
    float mean;
    float amp;
    float phase;
    /* omega is shared across all coefficients (breathing frequency) */
} TemporalModel;

/* ================================================================
   Kernel: Fit multi-field patches
   ================================================================
   One block per spatial patch (512 threads = 8³).
   Each thread reads all 6 field values at its voxel,
   computes the combined representation. */

/* Precomputed projection matrix: P = (V^T V)^{-1} V^T for tricubic fit */
__device__ __constant__ float d_vs2_P[VS2_BS][VS2_NC];      /* cubic: 8×4 */

static void vs2_gpu_init_projection(void) {
    /* Chebyshev basis: T0(s)=1, T1(s)=s, T2(s)=2s²-1, T3(s)=4s³-3s
     * where s = 2*t - 1 maps [0,1] -> [-1,1].
     * Condition number ~3.8 in 3D vs ~870k for monomials. */
    double V[VS2_BS][VS2_NC];
    for (int i = 0; i < VS2_BS; i++) {
        double t = (double)i / (VS2_BS - 1);
        double s = 2*t - 1;
        V[i][0]=1; V[i][1]=s; V[i][2]=2*s*s-1; V[i][3]=4*s*s*s-3*s;
    }
    double VtV[VS2_NC][VS2_NC] = {};
    for (int a = 0; a < VS2_NC; a++)
        for (int b = 0; b < VS2_NC; b++)
            for (int i = 0; i < VS2_BS; i++) VtV[a][b] += V[i][a]*V[i][b];

    double aug[VS2_NC][2*VS2_NC];
    for (int a = 0; a < VS2_NC; a++) {
        for (int b = 0; b < VS2_NC; b++) {
            aug[a][b] = VtV[a][b];
            aug[a][VS2_NC+b] = (a==b)?1.0:0.0;
        }
    }
    for (int col = 0; col < VS2_NC; col++) {
        double d = aug[col][col];
        for (int j = 0; j < 2*VS2_NC; j++) aug[col][j] /= d;
        for (int row = 0; row < VS2_NC; row++) {
            if (row==col) continue;
            double f = aug[row][col];
            for (int j = 0; j < 2*VS2_NC; j++) aug[row][j] -= f*aug[col][j];
        }
    }
    /* Projection matrix P = (V^T V)^{-1} V^T — computed in double, stored as float */
    float P[VS2_BS][VS2_NC];
    for (int i = 0; i < VS2_BS; i++)
        for (int a = 0; a < VS2_NC; a++) {
            double sum = 0;
            for (int m = 0; m < VS2_NC; m++) sum += aug[a][VS2_NC+m]*V[i][m];
            P[i][a] = (float)sum;
        }
    cudaMemcpyToSymbol(d_vs2_P, P, sizeof(P));
}

/* Multi-field patch fitting kernel: 6 independent tricubic fits per 8³ block.
 * Each block = one spatial patch (512 threads = 8³ voxels).
 * Thread tid = di*64 + dj*8 + dk.
 * Output layout: [phi_x(64) | phi_y(64) | phi_z(64) | theta_x(64) | theta_y(64) | theta_z(64)]
 * Same layout as CPU kernel for reader compatibility. */
__global__ void vs2_fit_multi_patches(
    const double *__restrict__ phi0,
    const double *__restrict__ phi1,
    const double *__restrict__ phi2,
    const double *__restrict__ theta0,
    const double *__restrict__ theta1,
    const double *__restrict__ theta2,
    int N, int BN,
    float *__restrict__ out_coeffs  /* n_patches × VS2_MULTI_TOTAL (384) */
) {
    int patch_idx = blockIdx.x;
    int bk = patch_idx % BN;
    int bj = (patch_idx / BN) % BN;
    int bi = patch_idx / (BN * BN);

    int di = threadIdx.x / 64;
    int dj = (threadIdx.x / 8) % 8;
    int dk = threadIdx.x % 8;

    int gi = bi*8+di, gj = bj*8+dj, gk = bk*8+dk;
    int valid = (gi < N && gj < N && gk < N);
    long idx = valid ? (long)gi*N*N + gj*N + gk : 0;

    /* Read all 6 fields */
    double fields[6];
    fields[0] = valid ? phi0[idx] : 0;
    fields[1] = valid ? phi1[idx] : 0;
    fields[2] = valid ? phi2[idx] : 0;
    fields[3] = valid ? theta0[idx] : 0;
    fields[4] = valid ? theta1[idx] : 0;
    fields[5] = valid ? theta2[idx] : 0;

    /* Shared memory: 6 fields × 64 coefficients — accumulate in float
     * (double atomicAdd on shared memory may silently fail on some GPUs) */
    __shared__ float s_coeffs[VS2_NFIELDS][VS2_NCOEFFS];

    /* Zero shared memory */
    for (int i = threadIdx.x; i < VS2_NFIELDS * VS2_NCOEFFS; i += blockDim.x)
        ((float*)s_coeffs)[i] = 0.0f;
    __syncthreads();

    /* Chebyshev basis evaluation: P already contains (V^T V)^{-1} V^T */
    for (int f = 0; f < VS2_NFIELDS; f++) {
        float val = (float)fields[f];
        for (int a = 0; a < VS2_NC; a++)
        for (int b = 0; b < VS2_NC; b++)
        for (int c = 0; c < VS2_NC; c++) {
            float w = d_vs2_P[di][a] * d_vs2_P[dj][b] * d_vs2_P[dk][c] * val;
            atomicAdd(&s_coeffs[f][a*16 + b*4 + c], w);
        }
    }
    __syncthreads();

    /* Write output */
    long base = (long)patch_idx * VS2_MULTI_TOTAL;
    for (int i = threadIdx.x; i < VS2_MULTI_TOTAL; i += blockDim.x)
        out_coeffs[base + i] = ((float*)s_coeffs)[i];
}

/* ================================================================
   Temporal streaming: predict + residual
   ================================================================ */

/* Predict coefficients from temporal model: predicted = mean + amp*cos(omega*t + phase) */
__global__ void vs2_temporal_predict(
    const float *__restrict__ temporal_mean,   /* n_total */
    const float *__restrict__ temporal_amp,    /* n_total */
    const float *__restrict__ temporal_phase,  /* n_total */
    float omega, float t,
    float *__restrict__ predicted,             /* n_total */
    int n_total
) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n_total) return;
    predicted[idx] = temporal_mean[idx] + temporal_amp[idx] * cosf(omega * t + temporal_phase[idx]);
}

/* Compute residual = actual - predicted, and count non-trivial residuals */
__global__ void vs2_temporal_residual(
    const float *__restrict__ actual,       /* n_total */
    const float *__restrict__ predicted,    /* n_total */
    float *__restrict__ residual,           /* n_total */
    float tol,
    uint32_t *__restrict__ nonzero_patches, /* patch indices with |residual| > tol */
    uint32_t *__restrict__ n_nonzero,       /* atomic counter */
    int n_patches, int coeffs_per_patch
) {
    int patch = blockIdx.x * blockDim.x + threadIdx.x;
    if (patch >= n_patches) return;

    long base = (long)patch * coeffs_per_patch;
    float max_r = 0;
    for (int c = 0; c < coeffs_per_patch; c++) {
        float r = actual[base + c] - predicted[base + c];
        residual[base + c] = r;
        float ar = fabsf(r);
        if (ar > max_r) max_r = ar;
    }

    if (max_r > tol) {
        uint32_t slot = atomicAdd(n_nonzero, 1);
        nonzero_patches[slot] = patch;
    }
}

/* Update temporal model with new data point: running Fourier accumulation */
__global__ void vs2_temporal_update(
    const float *__restrict__ actual,       /* n_total current coefficients */
    float *__restrict__ sum_cos,            /* n_total running cos sum */
    float *__restrict__ sum_sin,            /* n_total running sin sum */
    float *__restrict__ sum_mean,           /* n_total running mean sum */
    float cos_wt, float sin_wt,
    int n_total
) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n_total) return;
    float v = actual[idx];
    sum_mean[idx] += v;
    sum_cos[idx] += v * cos_wt;
    sum_sin[idx] += v * sin_wt;
}

/* Refit temporal model from accumulated sums */
__global__ void vs2_temporal_refit(
    const float *__restrict__ sum_cos,
    const float *__restrict__ sum_sin,
    const float *__restrict__ sum_mean,
    int count,
    float *__restrict__ out_mean,
    float *__restrict__ out_amp,
    float *__restrict__ out_phase,
    int n_total
) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n_total) return;

    float mean = sum_mean[idx] / count;
    float rc = sum_cos[idx] / count - mean; /* approximate */
    float rs = sum_sin[idx] / count;
    float amp = 2.0f * sqrtf(rc*rc + rs*rs);
    float phase = atan2f(-rs, rc);

    out_mean[idx] = mean;
    out_amp[idx] = amp;
    out_phase[idx] = phase;
}

/* Copy kernel */
__global__ void vs2_copy(const float *src, float *dst, int n) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < n) dst[i] = src[i];
}

#endif /* VECSTREAM_GPU_V2_CUH */
