/*  sfa2vecstream.c — Convert SFA voxel files to vecstream format
 *
 *  Usage: sfa2vecstream input.sfa -o output.vecstream [-col 0,1,2,3,4,5]
 *         [-bs 8] [-delta_tol 0.01] [-keyframe 100]
 *
 *  Build: gcc -O3 -march=native -fopenmp -o sfa2vecstream sfa2vecstream.c -lzstd -lm
 */

#define SFA_IMPLEMENTATION
#include "../format/sfa.h"

#define VECSTREAM_IMPLEMENTATION
#include "../format/vecstream.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifdef _OPENMP
#include <omp.h>
#endif

/* Extract column c, voxel i as float from raw SFA frame buffer */
static float col_f(void *buf, SFA *sfa, int c, long i) {
    long N3 = (long)sfa->Nx * sfa->Ny * sfa->Nz;
    uint64_t off = 0;
    for (int cc = 0; cc < c; cc++)
        off += (uint64_t)N3 * sfa_dtype_size[sfa->columns[cc].dtype];
    int dt = sfa->columns[c].dtype;
    uint8_t *src = (uint8_t*)buf + off;
    if (dt == SFA_F16) return vs_f16_to_f32(((uint16_t*)src)[i]);
    if (dt == SFA_F32) return ((float*)src)[i];
    return (float)((double*)src)[i];
}

/* Parse comma-separated column list: "0,1,2" -> [0,1,2], return count */
static int parse_cols(const char *s, int *cols, int max) {
    int n = 0;
    const char *p = s;
    while (*p && n < max) {
        cols[n++] = atoi(p);
        while (*p && *p != ',') p++;
        if (*p == ',') p++;
    }
    return n;
}

int main(int argc, char **argv) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s input.sfa -o output.vecstream [-col 0,1,2,3,4,5]\n"
                        "       [-bs 8] [-delta_tol 0.01] [-keyframe 100]\n", argv[0]);
        return 1;
    }

    const char *input_path = argv[1];
    const char *output_path = NULL;
    int cols[64], n_cols = 0;
    int BS = 8;
    float delta_tol = 0.01f;
    int keyframe_interval = 100;

    for (int i = 2; i < argc; i++) {
        if (!strcmp(argv[i], "-o") && i+1 < argc) output_path = argv[++i];
        else if (!strcmp(argv[i], "-col") && i+1 < argc) n_cols = parse_cols(argv[++i], cols, 64);
        else if (!strcmp(argv[i], "-bs") && i+1 < argc) BS = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-delta_tol") && i+1 < argc) delta_tol = (float)atof(argv[++i]);
        else if (!strcmp(argv[i], "-keyframe") && i+1 < argc) keyframe_interval = atoi(argv[++i]);
    }

    /* Open SFA */
    SFA *sfa = sfa_open(input_path);
    if (!sfa) { fprintf(stderr, "Cannot open %s\n", input_path); return 1; }

    uint32_t Nx = sfa->Nx, Ny = sfa->Ny, Nz = sfa->Nz;
    long N3 = (long)Nx * Ny * Nz;
    int nf = sfa->total_frames;

    /* Default: all 6 position+angle fields, or all columns */
    if (n_cols == 0) {
        n_cols = (int)sfa->n_columns;
        if (n_cols > 6) n_cols = 6; /* first 6 by default */
        for (int i = 0; i < n_cols; i++) cols[i] = i;
    }

    /* Default output path */
    char default_out[1024];
    if (!output_path) {
        snprintf(default_out, sizeof(default_out), "%s.vecstream", input_path);
        output_path = default_out;
    }

    printf("sfa2vecstream: %s -> %s\n", input_path, output_path);
    printf("  Grid: %u x %u x %u (%ld voxels)\n", Nx, Ny, Nz, N3);
    printf("  Frames: %d, Fields: %d, BS: %d\n", nf, n_cols, BS);
    printf("  Delta tolerance: %.4f, Keyframe interval: %d\n", delta_tol, keyframe_interval);

    /* Create vecstream */
    VecStream *vs = vecstream_create(output_path, Nx, Ny, Nz,
                                     sfa->Lx, sfa->Ly, sfa->Lz,
                                     sfa->dt > 0 ? sfa->dt : 1.0,
                                     (uint16_t)n_cols, (uint16_t)BS);
    if (!vs) { fprintf(stderr, "Cannot create %s\n", output_path); sfa_close(sfa); return 1; }

    /* Block counts */
    uint32_t bpx = (Nx + BS - 1) / BS;
    uint32_t bpy = (Ny + BS - 1) / BS;
    uint32_t bpz = (Nz + BS - 1) / BS;
    uint32_t n_patches = bpx * bpy * bpz;

    /* Allocate working storage */
    void *sfa_buf = malloc(sfa->frame_bytes);
    float *field = (float *)malloc(N3 * sizeof(float));

    /* Per-field current patches for delta computation */
    VecPatch **current = (VecPatch **)malloc(n_cols * sizeof(VecPatch *));
    for (int f = 0; f < n_cols; f++)
        current[f] = (VecPatch *)calloc(n_patches, sizeof(VecPatch));

    /* Working arrays for P-frame construction */
    uint32_t *delta_indices = (uint32_t *)malloc(n_patches * sizeof(uint32_t));
    float *delta_coeffs = (float *)malloc((size_t)n_patches * VS_NCOEFFS * sizeof(float));

    /* Stats */
    long total_sfa_bytes = 0;
    long total_vs_bytes_estimate = 0;

    printf("\n%5s  %3s  %4s  %7s  %10s  %10s\n",
           "Frame", "Fld", "Type", "Deltas", "MaxDelta", "PayloadKB");

    for (int fi = 0; fi < nf; fi++) {
        double t = sfa_frame_time(sfa, fi);
        if (sfa_read_frame(sfa, fi, sfa_buf) != 0) {
            fprintf(stderr, "Cannot read frame %d\n", fi);
            break;
        }

        for (int fc = 0; fc < n_cols; fc++) {
            int col = cols[fc];

            /* Extract field into float array */
            for (long i = 0; i < N3; i++)
                field[i] = col_f(sfa_buf, sfa, col, i);

            total_sfa_bytes += N3 * sizeof(float);

            /* Decide frame type */
            int is_first = (fi == 0);
            int is_keyframe = (!is_first && keyframe_interval > 0 &&
                               (fi % keyframe_interval) == 0);

            if (is_first || is_keyframe) {
                /* Vectorize all patches */
                VecPatch *patches = (VecPatch *)malloc(n_patches * sizeof(VecPatch));
                uint32_t pi = 0;
                for (uint32_t bi = 0; bi < bpx; bi++)
                for (uint32_t bj = 0; bj < bpy; bj++)
                for (uint32_t bk = 0; bk < bpz; bk++) {
                    vecstream_fit_patch(field, Nx, Ny, Nz,
                                       bi * BS, bj * BS, bk * BS,
                                       BS, BS, BS, VS_MAX_ORDER, &patches[pi]);
                    pi++;
                }

                /* Write I-frame */
                vecstream_write_iframe(vs, t, (uint8_t)fc, patches, n_patches);

                /* Update current state */
                memcpy(current[fc], patches, n_patches * sizeof(VecPatch));

                long payload_kb = (long)(4 + n_patches * (16 + VS_NCOEFFS * 4)) / 1024;
                total_vs_bytes_estimate += payload_kb * 1024;

                printf("%5d  %3d  %4s  %7u  %10s  %10ld\n",
                       fi, fc, "I", n_patches, "-", payload_kb);

                /* Also write K-frame for verification on keyframes */
                if (is_keyframe) {
                    vecstream_write_kframe(vs, t, (uint8_t)fc, field, 1);
                    printf("%5d  %3d  %4s  %7s  %10s  %10ld\n",
                           fi, fc, "K", "-", "-", N3 * 4L / 1024);
                }

                free(patches);
            } else {
                /* Vectorize and compute deltas */
                uint32_t n_deltas = 0;
                float max_delta = 0;

                for (uint32_t pi = 0; pi < n_patches; pi++) {
                    /* Find block origin for this patch index */
                    uint32_t bk = pi % bpz;
                    uint32_t bj = (pi / bpz) % bpy;
                    uint32_t bi = pi / (bpz * bpy);

                    VecPatch new_patch;
                    vecstream_fit_patch(field, Nx, Ny, Nz,
                                       bi * BS, bj * BS, bk * BS,
                                       BS, BS, BS, VS_MAX_ORDER, &new_patch);

                    /* Compute L-inf delta */
                    float patch_max = 0;
                    for (int c = 0; c < VS_NCOEFFS; c++) {
                        float d = fabsf(new_patch.coeffs[c] - current[fc][pi].coeffs[c]);
                        if (d > patch_max) patch_max = d;
                    }
                    if (patch_max > max_delta) max_delta = patch_max;

                    if (patch_max > delta_tol) {
                        delta_indices[n_deltas] = pi;
                        for (int c = 0; c < VS_NCOEFFS; c++)
                            delta_coeffs[(size_t)n_deltas * VS_NCOEFFS + c] =
                                new_patch.coeffs[c] - current[fc][pi].coeffs[c];
                        n_deltas++;
                    }

                    /* Update current state */
                    memcpy(&current[fc][pi], &new_patch, sizeof(VecPatch));
                }

                /* Write P-frame */
                vecstream_write_pframe(vs, t, (uint8_t)fc, delta_indices, delta_coeffs,
                                       n_deltas, VS_NCOEFFS);

                long payload_kb = (long)(8 + n_deltas * (4 + VS_NCOEFFS * 4)) / 1024;
                total_vs_bytes_estimate += payload_kb * 1024;

                printf("%5d  %3d  %4s  %7u  %10.4e  %10ld\n",
                       fi, fc, "P", n_deltas, max_delta, payload_kb);
            }
        }
    }

    /* Close and report */
    vecstream_close(vs);
    sfa_close(sfa);

    /* Get actual file sizes */
    FILE *fp = fopen(output_path, "rb");
    long vs_size = 0;
    if (fp) {
        fseek(fp, 0, SEEK_END);
        vs_size = ftell(fp);
        fclose(fp);
    }

    FILE *fp2 = fopen(input_path, "rb");
    long sfa_size = 0;
    if (fp2) {
        fseek(fp2, 0, SEEK_END);
        sfa_size = ftell(fp2);
        fclose(fp2);
    }

    printf("\n=== Conversion Summary ===\n");
    printf("  SFA file: %ld bytes (%.1f MB)\n", sfa_size, sfa_size / 1e6);
    printf("  VecStream file: %ld bytes (%.1f MB)\n", vs_size, vs_size / 1e6);
    printf("  Compression: %.1fx vs SFA\n", sfa_size > 0 ? (double)sfa_size / vs_size : 0);
    printf("  Raw voxel data: %ld bytes (%.1f MB)\n", total_sfa_bytes, total_sfa_bytes / 1e6);
    printf("  Compression vs raw: %.1fx\n", total_sfa_bytes > 0 ? (double)total_sfa_bytes / vs_size : 0);

    /* Cleanup */
    free(sfa_buf); free(field);
    free(delta_indices); free(delta_coeffs);
    for (int f = 0; f < n_cols; f++) free(current[f]);
    free(current);

    return 0;
}
