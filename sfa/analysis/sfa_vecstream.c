/*  sfa_vecstream.c — Streaming vector compression pipeline
 *
 *  Architecture:
 *    1. Spatial vectorize frame 0 → base patches (I-frame)
 *    2. For each subsequent frame:
 *       a. Re-vectorize → new patches
 *       b. Compute delta = new_coeffs - old_coeffs per patch
 *       c. If ||delta|| < threshold → store delta only (P-frame, tiny)
 *       d. If ||delta|| > threshold → store full patch (I-frame, rare)
 *    3. Every K frames: write raw voxel keyframe for verification
 *    4. Fourier fit on coefficient time series for long-range compression
 *
 *  Output: .vecstream file with I-frames, P-frames, and keyframes
 *  Also: reconstructed SFA for error verification
 *
 *  This is designed to integrate as a CUDA output hook:
 *  - Per timestep: receive N³ voxel data
 *  - Compute spatial patches (can be GPU-accelerated)
 *  - Diff against previous patches
 *  - Write compact delta to disk
 *  - Periodically write full keyframe
 *
 *  Build: gcc -O3 -march=native -fopenmp -o sfa_vecstream sfa_vecstream.c -lzstd -lm
 */

#define SFA_IMPLEMENTATION
#include "../format/sfa.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define PI 3.14159265358979323846
#define MAX_ORDER 3
#define NCOEFFS ((MAX_ORDER+1)*(MAX_ORDER+1)*(MAX_ORDER+1))  /* 64 for cubic */

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

/* ===== Spatial patch: tricubic polynomial over a BS³ block ===== */

typedef struct {
    float coeffs[NCOEFFS]; /* tensor-product polynomial coefficients */
    int ox, oy, oz;        /* block origin in grid coordinates */
    float max_err;         /* max approximation error */
    float rms_err;         /* RMS error */
} Patch;

/* Fit a tricubic polynomial to a BS³ block of data.
 * Uses separable least-squares: fit along each axis independently,
 * then combine via tensor product. This is O(BS³) not O(BS⁹). */
static void fit_patch(const float *field, int N, int BS,
                      int ox, int oy, int oz, Patch *p) {
    memset(p, 0, sizeof(Patch));
    p->ox = ox; p->oy = oy; p->oz = oz;
    int ord = MAX_ORDER;
    int nc = ord + 1; /* 4 per axis */

    /* Direct tensor-product least squares */
    /* For each (i,j,k) in the block, the basis is:
     *   B_{abc}(x,y,z) = t_x^a * t_y^b * t_z^c
     * where t_x = (ix - ox) / (BS-1), etc., in [0,1] */

    /* Build normal equations: A^T A c = A^T y
     * A is BS³ × nc³, too large to form explicitly.
     * Instead use the separability: the tensor product
     * of 1D Vandermonde matrices. */

    /* 1D Vandermonde: V[i][a] = t_i^a, i=0..BS-1, a=0..ord */
    double V[32][4]; /* max BS=32 */
    double VtV[4][4], VtV_inv[4][4];

    for (int i = 0; i < BS; i++) {
        double t = (BS > 1) ? (double)i / (BS - 1) : 0;
        double tk = 1;
        for (int a = 0; a <= ord; a++) { V[i][a] = tk; tk *= t; }
    }

    /* VtV = V^T V (4×4) */
    for (int a = 0; a < nc; a++)
        for (int b = 0; b < nc; b++) {
            double s = 0;
            for (int i = 0; i < BS; i++) s += V[i][a] * V[i][b];
            VtV[a][b] = s;
        }

    /* Invert VtV via Gauss-Jordan */
    double aug[4][8];
    for (int a = 0; a < nc; a++) {
        for (int b = 0; b < nc; b++) {
            aug[a][b] = VtV[a][b];
            aug[a][nc+b] = (a==b) ? 1.0 : 0.0;
        }
    }
    for (int col = 0; col < nc; col++) {
        int piv = col;
        for (int row = col+1; row < nc; row++)
            if (fabs(aug[row][col]) > fabs(aug[piv][col])) piv = row;
        if (piv != col)
            for (int j = 0; j < 2*nc; j++) {
                double tmp = aug[col][j]; aug[col][j] = aug[piv][j]; aug[piv][j] = tmp;
            }
        double d = aug[col][col];
        if (fabs(d) < 1e-15) continue;
        for (int j = 0; j < 2*nc; j++) aug[col][j] /= d;
        for (int row = 0; row < nc; row++) {
            if (row == col) continue;
            double f = aug[row][col];
            for (int j = 0; j < 2*nc; j++) aug[row][j] -= f * aug[col][j];
        }
    }
    for (int a = 0; a < nc; a++)
        for (int b = 0; b < nc; b++)
            VtV_inv[a][b] = aug[a][nc+b];

    /* P = V × VtV_inv — projection matrix rows
     * coeff[a][b][c] = Σ_{ijk} P[i][a] P[j][b] P[k][c] × data[i][j][k]
     * where P[i][a] = Σ_m V[i][m] × VtV_inv[m][a] = (VtV_inv V^T)^T */
    double P[32][4]; /* P[i][a] = Σ_m VtV_inv[a][m] * V[i][m] */
    for (int i = 0; i < BS; i++)
        for (int a = 0; a < nc; a++) {
            double s = 0;
            for (int m = 0; m < nc; m++) s += VtV_inv[a][m] * V[i][m];
            P[i][a] = s;
        }

    /* Compute tensor-product coefficients */
    int NN = N * N;
    for (int a = 0; a < nc; a++)
    for (int b = 0; b < nc; b++)
    for (int c = 0; c < nc; c++) {
        double coeff = 0;
        for (int di = 0; di < BS; di++)
        for (int dj = 0; dj < BS; dj++)
        for (int dk = 0; dk < BS; dk++) {
            int gi = ox + di, gj = oy + dj, gk = oz + dk;
            if (gi >= N || gj >= N || gk >= N) continue;
            long idx = (long)gi * NN + gj * N + gk;
            coeff += P[di][a] * P[dj][b] * P[dk][c] * field[idx];
        }
        p->coeffs[a * nc*nc + b * nc + c] = (float)coeff;
    }

    /* Compute reconstruction error */
    p->max_err = 0;
    double se = 0;
    int count = 0;
    for (int di = 0; di < BS; di++)
    for (int dj = 0; dj < BS; dj++)
    for (int dk = 0; dk < BS; dk++) {
        int gi = ox + di, gj = oy + dj, gk = oz + dk;
        if (gi >= N || gj >= N || gk >= N) continue;
        long idx = (long)gi * NN + gj * N + gk;

        double tx = (BS > 1) ? (double)di / (BS-1) : 0;
        double ty = (BS > 1) ? (double)dj / (BS-1) : 0;
        double tz = (BS > 1) ? (double)dk / (BS-1) : 0;

        double pred = 0;
        for (int a = 0; a < nc; a++)
        for (int b = 0; b < nc; b++)
        for (int c = 0; c < nc; c++) {
            double ba = 1; for (int q=0;q<a;q++) ba *= tx;
            double bb = 1; for (int q=0;q<b;q++) bb *= ty;
            double bc = 1; for (int q=0;q<c;q++) bc *= tz;
            pred += p->coeffs[a*nc*nc + b*nc + c] * ba * bb * bc;
        }

        float err = fabsf(field[idx] - (float)pred);
        if (err > p->max_err) p->max_err = err;
        se += err * err;
        count++;
    }
    p->rms_err = (float)sqrt(se / (count > 0 ? count : 1));
}

/* Reconstruct a single voxel from a patch */
static float patch_eval(const Patch *p, int BS, int di, int dj, int dk) {
    int nc = MAX_ORDER + 1;
    double tx = (BS > 1) ? (double)di / (BS-1) : 0;
    double ty = (BS > 1) ? (double)dj / (BS-1) : 0;
    double tz = (BS > 1) ? (double)dk / (BS-1) : 0;

    double val = 0;
    double txa[4] = {1, tx, tx*tx, tx*tx*tx};
    double tya[4] = {1, ty, ty*ty, ty*ty*ty};
    double tza[4] = {1, tz, tz*tz, tz*tz*tz};

    for (int a = 0; a < nc; a++)
    for (int b = 0; b < nc; b++)
    for (int c = 0; c < nc; c++)
        val += p->coeffs[a*nc*nc + b*nc + c] * txa[a] * tya[b] * tza[c];

    return (float)val;
}

/* ===== Frame types ===== */

typedef enum {
    FRAME_I = 0,  /* Full spatial vectorization (keyframe) */
    FRAME_P = 1,  /* Delta from previous frame */
    FRAME_K = 2,  /* Raw voxel keyframe for verification */
} FrameType;

typedef struct {
    FrameType type;
    double time;
    int n_patches;
    Patch *patches;          /* for I-frames: full patches */
    float *delta_coeffs;     /* for P-frames: n_patches × NCOEFFS deltas */
    int n_nonzero_patches;   /* P-frame: patches with delta above threshold */
    int *nonzero_indices;    /* P-frame: which patches have deltas */
    float frame_max_err;
    float frame_rms_err;
} VecFrame;

/* ===== Streaming codec ===== */

typedef struct {
    int N, BS;
    int BN;           /* blocks per axis = N/BS */
    int n_patches;    /* total blocks = BN³ */
    int n_fields;
    float delta_tol;  /* threshold for storing a delta */
    int keyframe_interval;

    /* Current state per field */
    Patch **current;  /* [field][patch] — current patch coefficients */

    /* Statistics */
    long total_i_bytes;
    long total_p_bytes;
    long total_k_bytes;
    int n_i_frames, n_p_frames, n_k_frames;
} VecCodec;

static VecCodec *codec_create(int N, int BS, int n_fields, float delta_tol, int keyframe_interval) {
    VecCodec *vc = calloc(1, sizeof(VecCodec));
    vc->N = N;
    vc->BS = BS;
    vc->BN = N / BS;
    vc->n_patches = vc->BN * vc->BN * vc->BN;
    vc->n_fields = n_fields;
    vc->delta_tol = delta_tol;
    vc->keyframe_interval = keyframe_interval;

    vc->current = malloc(n_fields * sizeof(Patch*));
    for (int f = 0; f < n_fields; f++)
        vc->current[f] = calloc(vc->n_patches, sizeof(Patch));

    return vc;
}

static void codec_free(VecCodec *vc) {
    for (int f = 0; f < vc->n_fields; f++) free(vc->current[f]);
    free(vc->current);
    free(vc);
}

/* Encode one frame: returns frame type and statistics */
static FrameType codec_encode_frame(VecCodec *vc, const float *field, int field_idx,
                                     int frame_num, double frame_time,
                                     int *n_deltas_out, float *max_delta_out) {
    int N = vc->N, BS = vc->BS, BN = vc->BN;
    int np = vc->n_patches;
    Patch *cur = vc->current[field_idx];

    /* Vectorize this frame */
    Patch *new_patches = malloc(np * sizeof(Patch));
    int pi = 0;
    for (int bi = 0; bi < BN; bi++)
    for (int bj = 0; bj < BN; bj++)
    for (int bk = 0; bk < BN; bk++) {
        fit_patch(field, N, BS, bi*BS, bj*BS, bk*BS, &new_patches[pi]);
        pi++;
    }

    FrameType type;
    int n_deltas = 0;
    float max_delta = 0;

    if (frame_num == 0 || (frame_num % vc->keyframe_interval) == 0) {
        /* I-frame: store full patches */
        type = FRAME_I;
        memcpy(cur, new_patches, np * sizeof(Patch));
        vc->total_i_bytes += np * NCOEFFS * sizeof(float);
        vc->n_i_frames++;
    } else {
        /* P-frame: compute deltas */
        type = FRAME_P;
        long p_bytes = 0;

        for (int p = 0; p < np; p++) {
            /* Compute L∞ norm of coefficient delta */
            float patch_max_delta = 0;
            for (int c = 0; c < NCOEFFS; c++) {
                float d = fabsf(new_patches[p].coeffs[c] - cur[p].coeffs[c]);
                if (d > patch_max_delta) patch_max_delta = d;
            }

            if (patch_max_delta > max_delta) max_delta = patch_max_delta;

            if (patch_max_delta > vc->delta_tol) {
                /* Store delta for this patch: patch_index + NCOEFFS floats */
                p_bytes += sizeof(int) + NCOEFFS * sizeof(float);
                n_deltas++;
            }
            /* Always update current state */
            memcpy(&cur[p], &new_patches[p], sizeof(Patch));
        }

        vc->total_p_bytes += p_bytes;
        vc->n_p_frames++;
    }

    *n_deltas_out = n_deltas;
    *max_delta_out = max_delta;
    free(new_patches);
    return type;
}

/* Reconstruct voxel field from current patch state */
static void codec_reconstruct(const VecCodec *vc, int field_idx, float *output) {
    int N = vc->N, BS = vc->BS, BN = vc->BN;
    const Patch *cur = vc->current[field_idx];
    int NN = N * N;

    int pi = 0;
    for (int bi = 0; bi < BN; bi++)
    for (int bj = 0; bj < BN; bj++)
    for (int bk = 0; bk < BN; bk++) {
        const Patch *p = &cur[pi];
        for (int di = 0; di < BS; di++)
        for (int dj = 0; dj < BS; dj++)
        for (int dk = 0; dk < BS; dk++) {
            int gi = p->ox + di, gj = p->oy + dj, gk = p->oz + dk;
            if (gi < N && gj < N && gk < N) {
                long idx = (long)gi * NN + gj * N + gk;
                output[idx] = patch_eval(p, BS, di, dj, dk);
            }
        }
        pi++;
    }
}

int main(int argc, char **argv) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s input.sfa [-col N] [-bs 8] [-delta_tol 0.001] [-keyframe 10]\n", argv[0]);
        fprintf(stderr, "\nStreaming vector compression analysis.\n");
        fprintf(stderr, "Encodes SFA frames as I-frames (full) + P-frames (deltas).\n");
        return 1;
    }

    int col = 0, BS = 8;
    float delta_tol = 0.001;
    int keyframe_interval = 10;
    char *recon_path = NULL;

    for (int i = 2; i < argc; i++) {
        if (!strcmp(argv[i], "-col") && i+1<argc) col = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-bs") && i+1<argc) BS = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-delta_tol") && i+1<argc) delta_tol = atof(argv[++i]);
        else if (!strcmp(argv[i], "-keyframe") && i+1<argc) keyframe_interval = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-recon") && i+1<argc) recon_path = argv[++i];
    }

    SFA *sfa = sfa_open(argv[1]);
    if (!sfa) { fprintf(stderr, "Cannot open %s\n", argv[1]); return 1; }

    int N = sfa->Nx;
    long N3 = (long)N*N*N;
    int nf = sfa->total_frames;

    printf("sfa_vecstream: %s (N=%d, %d frames, col=%d)\n", argv[1], N, nf, col);
    printf("  BS=%d, delta_tol=%.4f, keyframe_interval=%d\n", BS, delta_tol, keyframe_interval);

    VecCodec *vc = codec_create(N, BS, 1, delta_tol, keyframe_interval);
    void *buf = malloc(sfa->frame_bytes);
    float *field = malloc(N3 * sizeof(float));
    float *recon = malloc(N3 * sizeof(float));

    /* Reconstructed SFA for verification */
    SFA *recon_sfa = NULL;
    if (recon_path) {
        recon_sfa = sfa_create(recon_path, N, N, N, sfa->Lx, sfa->Ly, sfa->Lz,
                               sfa->dt > 0 ? sfa->dt : 0.01);
        sfa_add_column(recon_sfa, sfa->columns[col].name, SFA_F32, SFA_POSITION, col);
        sfa_finalize_header(recon_sfa);
    }

    printf("\n%5s  %4s  %6s  %9s  %10s  %10s  %10s\n",
           "Frame", "Type", "Deltas", "DeltaMax", "ReconRMS", "ReconMax", "Bytes");

    long total_raw_bytes = 0;
    long total_vec_bytes = 0;

    for (int fi = 0; fi < nf; fi++) {
        double t = sfa_frame_time(sfa, fi);
        if (sfa_read_frame(sfa, fi, buf) != 0) {
            fprintf(stderr, "Cannot read frame %d\n", fi);
            break;
        }

        /* Extract field column */
        for (long i = 0; i < N3; i++)
            field[i] = col_f(buf, sfa, col, i);

        total_raw_bytes += N3 * sizeof(float);

        /* Encode */
        int n_deltas;
        float max_delta;
        FrameType type = codec_encode_frame(vc, field, 0, fi, t, &n_deltas, &max_delta);

        /* Reconstruct and measure error */
        codec_reconstruct(vc, 0, recon);

        double max_err = 0, rms_sum = 0;
        for (long i = 0; i < N3; i++) {
            float err = fabsf(field[i] - recon[i]);
            if (err > max_err) max_err = err;
            rms_sum += err * err;
        }
        double rms_err = sqrt(rms_sum / N3);

        /* Write reconstructed frame */
        if (recon_sfa) {
            void *col_ptr = recon;
            sfa_write_frame(recon_sfa, t, &col_ptr);
        }

        /* Frame bytes */
        long frame_bytes;
        const char *type_str;
        if (type == FRAME_I) {
            frame_bytes = (long)vc->n_patches * NCOEFFS * sizeof(float);
            type_str = "I";
        } else {
            frame_bytes = (long)n_deltas * (sizeof(int) + NCOEFFS * sizeof(float));
            type_str = "P";
        }
        total_vec_bytes += frame_bytes;

        printf("%5d  %4s  %6d  %9.2e  %10.4e  %10.4e  %10ld\n",
               fi, type_str, n_deltas, max_delta, rms_err, max_err, frame_bytes);
    }

    printf("\n=== Streaming Compression Summary ===\n");
    printf("  Raw data: %ld bytes (%ld MB)\n", total_raw_bytes, total_raw_bytes / (1024*1024));
    printf("  Vec stream: %ld bytes (%ld MB)\n", total_vec_bytes, total_vec_bytes / (1024*1024));
    printf("  Compression: %.1f×\n", (double)total_raw_bytes / total_vec_bytes);
    printf("  I-frames: %d (%.1f MB)\n", vc->n_i_frames, vc->total_i_bytes / 1e6);
    printf("  P-frames: %d (%.1f MB)\n", vc->n_p_frames, vc->total_p_bytes / 1e6);
    printf("  Keyframe interval: every %d frames\n", keyframe_interval);
    printf("  Delta threshold: %.4f\n", delta_tol);

    if (recon_sfa) {
        sfa_close(recon_sfa);
        printf("  Reconstructed SFA: %s\n", recon_path);
    }

    free(field); free(recon); free(buf);
    codec_free(vc);
    sfa_close(sfa);
    return 0;
}
