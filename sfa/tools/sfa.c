/*  sfa.c — SFA multi-tool
 *
 *  Unified command-line tool for SFA file operations.
 *  Subcommands dispatched by first argument.
 *
 *  Usage: sfa <command> [args...]
 *
 *  Commands:
 *    preview   Convert full SFA to uint8 preview for fast viewing
 *    info      (future) Print SFA file metadata
 *    extract   (future) Extract/truncate frames
 *
 *  Build: gcc -O3 -fopenmp -o sfa sfa.c -lzstd -lm
 */

#define SFA_IMPLEMENTATION
#include "../format/sfa.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>

/* ================================================================
   preview — Convert full SFA to uint8 preview
   ================================================================ */

static void usage_preview(void) {
    fprintf(stderr,
        "Usage: sfa preview <input.sfa> <output.sfa> [options]\n"
        "\n"
        "Create a small preview SFA with pre-computed uint8 channels.\n"
        "The output reuses the SFA container format and can be opened\n"
        "directly in volview.\n"
        "\n"
        "Options:\n"
        "  --N <size>        Preview grid size (default: 64)\n"
        "  --sample <n>      Sample every nth frame (default: 1 = all)\n"
        "  --channels <list> Channels: P,phi,theta,curl (default: P,phi,theta)\n"
        "\n"
        "Output columns (uint8, quantized 0-255):\n"
        "  P_abs     |phi_0 * phi_1 * phi_2|  (binding density)\n"
        "  phi_sq    phi_0^2 + phi_1^2 + phi_2^2  (field energy)\n"
        "  theta_sq  theta_0^2 + theta_1^2 + theta_2^2  (EM energy)\n"
        "\n"
        "All source KVMD parameters are preserved. Preview metadata\n"
        "(quantization ranges, source grid size) added as KVMD set 1.\n"
        "\n");
}

/* Downsample a float volume from src_N³ to dst_N³ via box averaging */
static void downsample_volume(const float *src, int src_N, float *dst, int dst_N) {
    double scale = (double)src_N / dst_N;
    int src_NN = src_N * src_N;
    int dst_NN = dst_N * dst_N;

    #pragma omp parallel for schedule(static) collapse(3)
    for (int di = 0; di < dst_N; di++) {
        for (int dj = 0; dj < dst_N; dj++) {
            for (int dk = 0; dk < dst_N; dk++) {
                int si0 = (int)(di * scale), si1 = (int)((di + 1) * scale);
                int sj0 = (int)(dj * scale), sj1 = (int)((dj + 1) * scale);
                int sk0 = (int)(dk * scale), sk1 = (int)((dk + 1) * scale);
                if (si1 > src_N) si1 = src_N;
                if (sj1 > src_N) sj1 = src_N;
                if (sk1 > src_N) sk1 = src_N;

                double sum = 0;
                int count = 0;
                for (int si = si0; si < si1; si++)
                    for (int sj = sj0; sj < sj1; sj++)
                        for (int sk = sk0; sk < sk1; sk++) {
                            sum += src[(long)si * src_NN + sj * src_N + sk];
                            count++;
                        }
                dst[(long)di * dst_NN + dj * dst_N + dk] = (count > 0) ? (float)(sum / count) : 0;
            }
        }
    }
}

/* Quantize float array to uint8 given range [0, max_val] → [0, 255] */
static void quantize_u8(const float *src, uint8_t *dst, long n, float max_val) {
    float scale = (max_val > 0) ? 255.0f / max_val : 0;
    for (long i = 0; i < n; i++) {
        float v = src[i] * scale;
        if (v < 0) v = 0;
        if (v > 255) v = 255;
        dst[i] = (uint8_t)(v + 0.5f);
    }
}

static int cmd_preview(int argc, char **argv) {
    if (argc < 2) { usage_preview(); return 1; }

    const char *input_path = argv[0];
    const char *output_path = argv[1];
    int preview_N = 64;
    int sample_every = 1;

    /* Parse options */
    for (int i = 2; i < argc; i++) {
        if (!strcmp(argv[i], "--N") && i + 1 < argc) preview_N = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--sample") && i + 1 < argc) sample_every = atoi(argv[++i]);
        else { fprintf(stderr, "Unknown option: %s\n", argv[i]); return 1; }
    }

    /* Open source */
    SFA *src = sfa_open(input_path);
    if (!src) { fprintf(stderr, "Cannot open %s\n", input_path); return 1; }

    int src_N = src->Nx;
    long src_N3 = (long)src_N * src_N * src_N;
    long dst_N3 = (long)preview_N * preview_N * preview_N;
    uint32_t n_frames = src->total_frames;
    double src_L = src->Lx;

    /* For streaming files with truncated frames, just skip bad reads
     * during processing — don't try to validate upfront. */

    printf("sfa preview: %s → %s\n", input_path, output_path);
    printf("  Source:  N=%d, %d frames, %.1f MB/frame\n",
           src_N, n_frames, src_N3 * 12.0 * 4 / 1e6);
    printf("  Preview: N=%d, sample every %d → %d frames\n",
           preview_N, sample_every, (n_frames + sample_every - 1) / sample_every);
    printf("  Channels: P_abs, phi_sq, theta_sq (uint8)\n\n");

    if (src->n_columns < 6) {
        fprintf(stderr, "Need at least 6 columns (3 phi + 3 theta)\n");
        sfa_close(src);
        return 1;
    }

    /* Pass 1: Scan a few frames to find quantization ranges */
    printf("Scanning quantization ranges...\n");
    float max_P = 0, max_phi = 0, max_theta = 0;

    /* Read frame data — sfa_read_frame returns decoded column-major data */
    size_t src_frame_bytes = src->frame_bytes;
    if (src_frame_bytes == 0) {
        /* Compute if not set */
        for (uint32_t c = 0; c < src->n_columns; c++)
            src_frame_bytes += src_N3 * sfa_dtype_size[src->columns[c].dtype];
    }
    printf("  Frame buffer: %.1f MB (%zu bytes)\n", src_frame_bytes / 1e6, src_frame_bytes);
    void *frame_buf = malloc(src_frame_bytes);
    if (!frame_buf) { fprintf(stderr, "Cannot allocate frame buffer\n"); return 1; }

    int scan_frames[] = {0, (int)(n_frames / 4), (int)(n_frames / 2),
                         (int)(3 * n_frames / 4), (int)(n_frames - 1)};
    int n_scan = 5;
    for (int s = 0; s < n_scan; s++) {
        int fi = scan_frames[s];
        if (fi < 0 || fi >= (int)n_frames) continue;
        if (sfa_read_frame(src, fi, frame_buf) < 0) continue;

        /* Extract columns — assume f32 for phi and theta */
        size_t col_bytes = src_N3 * sizeof(float);
        float *phi[3], *theta[3];
        for (int a = 0; a < 3; a++) {
            phi[a] = (float*)((char*)frame_buf + a * col_bytes);
            theta[a] = (float*)((char*)frame_buf + (3 + a) * col_bytes);
        }

        for (long i = 0; i < src_N3; i++) {
            float P = fabsf(phi[0][i] * phi[1][i] * phi[2][i]);
            float ps = phi[0][i]*phi[0][i] + phi[1][i]*phi[1][i] + phi[2][i]*phi[2][i];
            float ts = theta[0][i]*theta[0][i] + theta[1][i]*theta[1][i] + theta[2][i]*theta[2][i];
            if (P > max_P) max_P = P;
            if (ps > max_phi) max_phi = ps;
            if (ts > max_theta) max_theta = ts;
        }
    }

    /* Add 10% headroom to avoid clipping */
    max_P *= 1.1f; max_phi *= 1.1f; max_theta *= 1.1f;
    printf("  Quantization: P_max=%.5f  phi_sq_max=%.3f  theta_sq_max=%.5f\n\n",
           max_P, max_phi, max_theta);

    /* Create output SFA */
    SFA *out = sfa_create(output_path, preview_N, preview_N, preview_N,
                          src_L, src_L, src_L, src->dt);
    if (!out) {
        fprintf(stderr, "Cannot create output file %s (check directory exists)\n", output_path);
        sfa_close(src); free(frame_buf);
        return 1;
    }
    out->flags = SFA_CODEC_COLZSTD | SFA_FLAG_STREAMING;

    /* Copy source KVMD as set 0 */
    SFA_KVMDSet *kv = calloc(SFA_MAX_KVMD_SETS, sizeof(SFA_KVMDSet));
    int n_kv = sfa_read_kvmd(src, kv, SFA_MAX_KVMD_SETS);
    if (n_kv > 0) {
        SFA_KVMDSet *use = &kv[0];
        for (int i = 0; i < n_kv; i++)
            if (kv[i].first_frame == 0xFFFFFFFF) { use = &kv[i]; break; }
        /* Build pointer arrays from char[128][64] arrays */
        const char *kptrs[128], *vptrs[128];
        for (int i = 0; i < use->n_pairs && i < 128; i++) {
            kptrs[i] = use->keys[i];
            vptrs[i] = use->values[i];
        }
        sfa_add_kvmd(out, 0, 0xFFFFFFFF, 0xFFFFFFFF,
                     kptrs, vptrs, use->n_pairs);
    }

    /* Add preview KVMD as set 1 */
    {
        char vPN[32], vSN[32], vPM[32], vPHM[32], vTM[32], vSE[32];
        snprintf(vPN, 32, "%d", preview_N);
        snprintf(vSN, 32, "%d", src_N);
        snprintf(vPM, 32, "%.6f", max_P);
        snprintf(vPHM, 32, "%.6f", max_phi);
        snprintf(vTM, 32, "%.6f", max_theta);
        snprintf(vSE, 32, "%d", sample_every);
        const char *keys[] = {"preview", "source_N", "preview_N",
                              "quant_P_max", "quant_phi_sq_max", "quant_theta_sq_max",
                              "sample_every"};
        const char *vals[] = {"true", vSN, vPN, vPM, vPHM, vTM, vSE};
        sfa_add_kvmd(out, 1, 0xFFFFFFFF, 0xFFFFFFFF, keys, vals, 7);
    }

    /* Add 3 uint8 columns */
    sfa_add_column(out, "P_abs",    SFA_U8, SFA_POSITION, 0);
    sfa_add_column(out, "phi_sq",   SFA_U8, SFA_POSITION, 1);
    sfa_add_column(out, "theta_sq", SFA_U8, SFA_POSITION, 2);
    sfa_finalize_header(out);

    /* Allocate working buffers */
    float *derived[3];  /* P, phi_sq, theta_sq at source resolution */
    float *downsampled[3];  /* at preview resolution */
    uint8_t *quantized[3];
    for (int c = 0; c < 3; c++) {
        derived[c] = malloc(src_N3 * sizeof(float));
        downsampled[c] = malloc(dst_N3 * sizeof(float));
        quantized[c] = malloc(dst_N3);
    }

    /* Pass 2: Convert frames */
    int out_frames = 0;
    double wall0 = omp_get_wtime();

    for (uint32_t fi = 0; fi < n_frames; fi += sample_every) {
        if (sfa_read_frame(src, fi, frame_buf) < 0) {
            fprintf(stderr, "  Warning: cannot read frame %d, skipping\n", fi);
            continue;
        }
        double t = sfa_frame_time(src, fi);

        /* Extract source columns */
        size_t col_bytes = src_N3 * sizeof(float);
        float *phi[3], *theta[3];
        for (int a = 0; a < 3; a++) {
            phi[a] = (float*)((char*)frame_buf + a * col_bytes);
            theta[a] = (float*)((char*)frame_buf + (3 + a) * col_bytes);
        }

        /* Compute derived at full resolution */
        #pragma omp parallel for schedule(static)
        for (long i = 0; i < src_N3; i++) {
            derived[0][i] = fabsf(phi[0][i] * phi[1][i] * phi[2][i]);
            derived[1][i] = phi[0][i]*phi[0][i] + phi[1][i]*phi[1][i] + phi[2][i]*phi[2][i];
            derived[2][i] = theta[0][i]*theta[0][i] + theta[1][i]*theta[1][i] + theta[2][i]*theta[2][i];
        }

        /* Downsample */
        for (int c = 0; c < 3; c++)
            downsample_volume(derived[c], src_N, downsampled[c], preview_N);

        /* Quantize */
        float maxvals[3] = {max_P, max_phi, max_theta};
        for (int c = 0; c < 3; c++)
            quantize_u8(downsampled[c], quantized[c], dst_N3, maxvals[c]);

        /* Write frame */
        void *cols[3] = {quantized[0], quantized[1], quantized[2]};
        sfa_write_frame(out, t, cols);
        out_frames++;

        if (out_frames % 100 == 0 || fi + sample_every >= n_frames) {
            double elapsed = omp_get_wtime() - wall0;
            double pct = 100.0 * (fi + 1) / n_frames;
            printf("  Frame %d/%d (t=%.2f) [%.0f%% %.1fs]\r",
                   out_frames, (int)((n_frames + sample_every - 1) / sample_every),
                   t, pct, elapsed);
            fflush(stdout);
        }
    }
    printf("\n");

    sfa_close(out);
    sfa_close(src);

    /* Report */
    double wall = omp_get_wtime() - wall0;
    long out_size = 0;
    FILE *fcheck = fopen(output_path, "rb");
    if (fcheck) { fseek(fcheck, 0, SEEK_END); out_size = ftell(fcheck); fclose(fcheck); }

    printf("\nDone: %d frames, %.1f MB (%.1f KB/frame)\n",
           out_frames, out_size / 1e6, out_size / 1e3 / out_frames);
    printf("Compression: %.0f:1 (source frame %.1f MB → preview %.1f KB)\n",
           (src_N3 * 12.0 * 4) / (out_size / (double)out_frames),
           src_N3 * 12.0 * 4 / 1e6,
           out_size / 1e3 / out_frames);
    printf("Wall: %.1fs (%.1f frames/sec)\n", wall, out_frames / wall);

    for (int c = 0; c < 3; c++) { free(derived[c]); free(downsampled[c]); free(quantized[c]); }
    free(frame_buf);
    free(kv);
    return 0;
}

/* ================================================================
   Main dispatcher
   ================================================================ */

static void usage_main(void) {
    fprintf(stderr,
        "Usage: sfa <command> [args...]\n"
        "\n"
        "Commands:\n"
        "  preview   Create uint8 preview SFA for fast viewing\n"
        "\n"
        "Run 'sfa <command> --help' for command-specific usage.\n"
        "\n");
}

int main(int argc, char **argv) {
    if (argc < 2) { usage_main(); return 1; }

    const char *cmd = argv[1];

    if (!strcmp(cmd, "preview"))
        return cmd_preview(argc - 2, argv + 2);
    else {
        fprintf(stderr, "Unknown command: %s\n\n", cmd);
        usage_main();
        return 1;
    }
}
