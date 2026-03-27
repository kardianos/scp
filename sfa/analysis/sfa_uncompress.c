/*  sfa_uncompress.c — Convert compressed SFA to raw (mmap-friendly) format
 *
 *  Writes a new SFA file with uncompressed frame data laid out sequentially
 *  so frames can be accessed via mmap + pointer arithmetic with zero
 *  decompression overhead. Each frame is stored as raw column data in the
 *  native dtype (f16/f32/f64), no BSS, no zstd.
 *
 *  The output file has:
 *    - Same SFAH, CDEF, KVMD headers
 *    - codec=0 (raw) in flags
 *    - FRMD chunks with uncompressed data
 *    - Frames laid out contiguously for efficient mmap access
 *
 *  Optionally converts f16 columns to f32 for direct float* mapping.
 *
 *  Build: gcc -O3 -o sfa_uncompress sfa_uncompress.c -lzstd -lm
 *
 *  Usage:
 *    sfa_uncompress input.sfa -o output_raw.sfa
 *    sfa_uncompress input.sfa -o output_raw.sfa --f32   (convert f16 → f32)
 */

#define SFA_IMPLEMENTATION
#include "../format/sfa.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>

static inline float f16_to_f32(uint16_t h) {
    uint16_t sign = h & 0x8000;
    int exp = (h >> 10) & 0x1F;
    uint16_t mant = h & 0x3FF;
    if (exp == 0) return 0.0f;
    if (exp == 31) return sign ? -1e30f : 1e30f;
    float fv;
    uint32_t x = ((uint32_t)sign << 16) | ((uint32_t)(exp-15+127) << 23) | ((uint32_t)mant << 13);
    memcpy(&fv, &x, 4);
    return fv;
}

int main(int argc, char **argv) {
    if (argc < 4) {
        fprintf(stderr,
            "Usage: %s input.sfa -o output.sfa [--f32] [--colzstd]\n"
            "\n"
            "  --f32:      promote f16 columns to f32.\n"
            "  --colzstd:  per-column zstd compression (parallel decompression, ~30x ratio).\n"
            "              Without --colzstd, writes raw (uncompressed) for mmap.\n"
            "\n", argv[0]);
        return 1;
    }

    const char *in_path = argv[1];
    const char *out_path = NULL;
    int promote_f16 = 0;
    int use_colzstd = 0;

    for (int i = 2; i < argc; i++) {
        if (!strcmp(argv[i], "-o") && i+1 < argc) out_path = argv[++i];
        else if (!strcmp(argv[i], "--f32")) promote_f16 = 1;
        else if (!strcmp(argv[i], "--colzstd")) use_colzstd = 1;
    }
    if (!out_path) { fprintf(stderr, "Missing -o output.sfa\n"); return 1; }

    /* Open source */
    struct stat st;
    if (stat(in_path, &st) != 0) { fprintf(stderr, "Cannot stat %s\n", in_path); return 1; }
    uint64_t file_bytes = st.st_size;

    SFA *src = sfa_open(in_path);
    if (!src) { fprintf(stderr, "Cannot open %s\n", in_path); return 1; }

    int valid = sfa_count_valid_frames(src, file_bytes);
    if (valid == 0) { fprintf(stderr, "No valid frames in %s\n", in_path); sfa_close(src); return 1; }

    long N3 = (long)src->Nx * src->Ny * src->Nz;
    printf("Input: %s\n", in_path);
    printf("  Grid: %ux%ux%u (%ld voxels)\n", src->Nx, src->Ny, src->Nz, N3);
    printf("  Columns: %u, Frames: %d valid\n", src->n_columns, valid);

    /* Determine output column dtypes */
    uint8_t out_dtypes[64];
    int out_elem_sizes[64];
    uint64_t out_frame_bytes = 0;
    for (uint32_t c = 0; c < src->n_columns; c++) {
        uint8_t dt = src->columns[c].dtype;
        if (promote_f16 && dt == SFA_F16) {
            out_dtypes[c] = SFA_F32;
            out_elem_sizes[c] = 4;
        } else {
            out_dtypes[c] = dt;
            out_elem_sizes[c] = sfa_dtype_size[dt];
        }
        out_frame_bytes += (uint64_t)N3 * out_elem_sizes[c];
    }

    double est_size = (double)valid * out_frame_bytes / 1e9;
    printf("  Output frame size: %.1f MB (uncompressed)\n", out_frame_bytes / 1e6);
    printf("  Estimated output: %.1f GB (%d frames)\n", est_size, valid);
    if (promote_f16) printf("  f16 → f32 promotion enabled\n");

    /* Create output SFA */
    SFA *dst = sfa_create(out_path, src->Nx, src->Ny, src->Nz,
                          src->Lx, src->Ly, src->Lz, src->dt);
    if (!dst) { fprintf(stderr, "Cannot create %s\n", out_path); sfa_close(src); return 1; }

    /* Set codec */
    dst->flags = use_colzstd ? SFA_CODEC_COLZSTD : SFA_CODEC_RAW;

    /* Copy column definitions (with dtype promotion if needed) */
    dst->n_columns = src->n_columns;
    for (uint32_t c = 0; c < src->n_columns; c++) {
        dst->columns[c] = src->columns[c];
        dst->columns[c].dtype = out_dtypes[c];
    }

    /* Copy KVMD */
    SFA_KVMDSet kv[SFA_MAX_KVMD_SETS];
    int n_kv = sfa_read_kvmd(src, kv, SFA_MAX_KVMD_SETS);
    for (int i = 0; i < n_kv; i++) {
        const char *kptrs[128], *vptrs[128];
        for (int j = 0; j < kv[i].n_pairs; j++) {
            kptrs[j] = kv[i].keys[j];
            vptrs[j] = kv[i].values[j];
        }
        sfa_add_kvmd(dst, kv[i].set_id, 0xFFFFFFFF, 0, kptrs, vptrs, kv[i].n_pairs);
    }

    /* Add conversion metadata */
    {
        const char *keys[] = {"uncompressed_from", "promote_f16"};
        const char *vals[] = {in_path, promote_f16 ? "true" : "false"};
        sfa_add_kvmd(dst, 99, 0xFFFFFFFF, 0, keys, vals, 2);
    }

    /* Pre-allocate total_frames for JMPF sizing */
    dst->total_frames = valid;
    sfa_finalize_header(dst);

    /* Read source frame buffer */
    void *src_buf = malloc(src->frame_bytes);
    if (!src_buf) { fprintf(stderr, "Cannot allocate source buffer\n"); return 1; }

    /* Output column buffer (one column at a time) */
    void *col_buf = malloc((uint64_t)N3 * 4);  /* max 4 bytes per element for f32 */

    /* Per-column output pointers for sfa_write_frame */
    void **out_cols = (void**)malloc(src->n_columns * sizeof(void*));
    uint8_t *out_frame = (uint8_t*)malloc(out_frame_bytes);

    /* Process each frame */
    for (int f = 0; f < valid; f++) {
        double t = sfa_frame_time(src, f);

        if (sfa_read_frame(src, f, src_buf) != 0) {
            fprintf(stderr, "Failed to read frame %d\n", f);
            break;
        }

        /* Extract and optionally convert each column */
        uint64_t src_off = 0;
        uint64_t dst_off = 0;
        for (uint32_t c = 0; c < src->n_columns; c++) {
            uint8_t src_dt = src->columns[c].dtype;
            int src_es = sfa_dtype_size[src_dt];
            int dst_es = out_elem_sizes[c];
            uint8_t *src_col = (uint8_t*)src_buf + src_off;
            uint8_t *dst_col = out_frame + dst_off;

            if (promote_f16 && src_dt == SFA_F16) {
                /* f16 → f32 conversion */
                uint16_t *in16 = (uint16_t*)src_col;
                float *out32 = (float*)dst_col;
                for (long i = 0; i < N3; i++) {
                    out32[i] = f16_to_f32(in16[i]);
                }
            } else {
                /* Direct copy */
                memcpy(dst_col, src_col, (uint64_t)N3 * src_es);
            }

            out_cols[c] = dst_col;
            src_off += (uint64_t)N3 * src_es;
            dst_off += (uint64_t)N3 * dst_es;
        }

        /* Write uncompressed frame — use raw codec write */
        sfa_write_frame_ex(dst, t, out_cols, 0, out_frame_bytes);

        printf("\r  Frame %d/%d (t=%.2f)", f+1, valid, t);
        fflush(stdout);
    }
    printf("\n");

    uint32_t nf = dst->total_frames;
    sfa_close(dst);
    sfa_close(src);

    /* Report */
    struct stat out_st;
    stat(out_path, &out_st);
    printf("Output: %s\n", out_path);
    printf("  %u frames, %.1f GB\n", nf, out_st.st_size / 1e9);
    printf("  Codec: %s\n", use_colzstd ? "colzstd (per-column parallel)" : "raw (mmap-ready)");
    if (promote_f16) printf("  All columns are f32 (directly castable to float*)\n");

    free(src_buf);
    free(col_buf);
    free(out_cols);
    free(out_frame);
    return 0;
}
