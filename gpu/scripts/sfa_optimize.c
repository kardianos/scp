/*  sfa_optimize.c — Read an SFA file and rewrite with better compression
 *
 *  Reads frames from input SFA (f64 BSS+zstd), rewrites as f32+BSS+zstd.
 *  Runs on the remote GPU instance to avoid downloading 16 GB.
 *
 *  Build: gcc -O2 -o sfa_optimize sfa_optimize.c -lzstd -lm
 *  Usage: ./sfa_optimize input.sfa output.sfa [--f32] [--skip N]
 */

#define SFA_IMPLEMENTATION
#include "../../sfa/format/sfa.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

static void downcast_f64_to_f32(const double *src, float *dst, size_t n) {
    for (size_t i = 0; i < n; i++)
        dst[i] = (float)src[i];
}

int main(int argc, char **argv) {
    if (argc < 3) {
        printf("Usage: %s input.sfa output.sfa [--f32] [--skip N] [--maxframes N]\n", argv[0]);
        printf("  --f32       Downcast f64 to f32 (lossy, ~17x compression)\n");
        printf("  --skip N    Keep every Nth frame (default 1 = all)\n");
        printf("  --maxframes N  Stop after N output frames\n");
        return 1;
    }

    const char *inpath = argv[1];
    const char *outpath = argv[2];
    int use_f32 = 0;
    int skip = 1;
    int maxframes = 999999;

    for (int i = 3; i < argc; i++) {
        if (!strcmp(argv[i], "--f32")) use_f32 = 1;
        else if (!strcmp(argv[i], "--skip") && i+1 < argc) skip = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--maxframes") && i+1 < argc) maxframes = atoi(argv[++i]);
    }

    /* Open input */
    SFA *in = sfa_open(inpath);
    if (!in) { fprintf(stderr, "Cannot open %s\n", inpath); return 1; }

    printf("Input: %dx%dx%d, %d columns, %d frames\n",
           in->Nx, in->Ny, in->Nz, in->n_columns, in->total_frames);
    printf("Input frame size: %.1f MB (f64)\n", in->frame_bytes / 1e6);

    /* Create output */
    SFA *out = sfa_create(outpath, in->Nx, in->Ny, in->Nz,
                          in->Lx, in->Ly, in->Lz, in->dt);

    /* Copy column definitions, optionally changing dtype */
    for (uint32_t c = 0; c < in->n_columns; c++) {
        uint8_t dtype = in->columns[c].dtype;
        if (use_f32 && dtype == SFA_F64) dtype = SFA_F32;
        sfa_add_column(out, in->columns[c].name, dtype,
                       in->columns[c].semantic, in->columns[c].component);
    }
    sfa_finalize_header(out);

    /* Allocate buffers */
    void *in_buf = malloc(in->frame_bytes);
    uint64_t out_N = (uint64_t)in->Nx * in->Ny * in->Nz;

    /* For f32 mode: allocate f32 column buffers */
    float **f32_cols = NULL;
    if (use_f32) {
        f32_cols = (float**)malloc(in->n_columns * sizeof(float*));
        for (uint32_t c = 0; c < in->n_columns; c++)
            f32_cols[c] = (float*)malloc(out_N * sizeof(float));
    }

    int out_frames = 0;
    printf("Converting: skip=%d, f32=%d, maxframes=%d\n", skip, use_f32, maxframes);

    for (uint32_t f = 0; f < in->total_frames && out_frames < maxframes; f++) {
        if (f % skip != 0) continue;

        double t = sfa_frame_time(in, f);
        int rc = sfa_read_frame(in, f, in_buf);
        if (rc < 0) {
            printf("  Frame %d: read error (rc=%d), skipping\n", f, rc);
            continue;
        }

        if (use_f32) {
            /* Downcast each column from f64 to f32 */
            uint64_t off = 0;
            for (uint32_t c = 0; c < in->n_columns; c++) {
                int in_size = sfa_dtype_size[in->columns[c].dtype];
                if (in->columns[c].dtype == SFA_F64) {
                    downcast_f64_to_f32((double*)((uint8_t*)in_buf + off),
                                       f32_cols[c], out_N);
                }
                off += out_N * in_size;
            }
            void **cols = (void**)malloc(in->n_columns * sizeof(void*));
            for (uint32_t c = 0; c < in->n_columns; c++)
                cols[c] = f32_cols[c];
            sfa_write_frame(out, t, cols);
            free(cols);
        } else {
            /* Split the buffer into column pointers */
            void **cols = (void**)malloc(in->n_columns * sizeof(void*));
            uint64_t off = 0;
            for (uint32_t c = 0; c < in->n_columns; c++) {
                cols[c] = (uint8_t*)in_buf + off;
                off += out_N * sfa_dtype_size[in->columns[c].dtype];
            }
            sfa_write_frame(out, t, cols);
            free(cols);
        }

        out_frames++;
        if (out_frames % 10 == 0 || out_frames == 1)
            printf("  Frame %d/%d (input %d) t=%.1f\n",
                   out_frames, in->total_frames/skip, f, t);
    }

    sfa_close(out);
    sfa_close(in);
    free(in_buf);
    if (f32_cols) {
        for (uint32_t c = 0; c < in->n_columns; c++) free(f32_cols[c]);
        free(f32_cols);
    }

    /* Report sizes */
    FILE *fi = fopen(inpath, "rb"); fseek(fi, 0, SEEK_END);
    long in_size = ftell(fi); fclose(fi);
    FILE *fo = fopen(outpath, "rb"); fseek(fo, 0, SEEK_END);
    long out_size = ftell(fo); fclose(fo);

    printf("\nDone: %d frames written\n", out_frames);
    printf("Input:  %.1f MB (%d frames)\n", in_size/1e6, in->total_frames);
    printf("Output: %.1f MB (%d frames)\n", out_size/1e6, out_frames);
    printf("Ratio:  %.1fx smaller\n", (double)in_size / out_size);

    return 0;
}
