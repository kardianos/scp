/*  sfa.c — SFA multi-tool
 *
 *  Unified command-line tool for SFA file operations.
 *  Subcommands dispatched by first argument.
 *
 *  Usage: sfa <command> [args...]
 *
 *  Commands:
 *    info      Print SFA file metadata and diagnostics
 *    extract   Extract a single frame to a new SFA seed file
 *    repair    Fix truncated SFA files (remove bad frames)
 *    convert   Convert between codecs (uncompress, f16→f32, recompress)
 *    preview   Create uint8 preview SFA for fast viewing
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
#include <sys/stat.h>

/* Shared helper: f16 → f64 conversion */
static inline double sfa_f16_to_f64(uint16_t h) {
    uint16_t sign = h & 0x8000; int exp = (h >> 10) & 0x1F; uint16_t mant = h & 0x3FF;
    if (exp == 0) return 0.0;
    if (exp == 31) return sign ? -1e30 : 1e30;
    float fv; uint32_t x = ((uint32_t)sign << 16) | ((uint32_t)(exp-15+127) << 23) | ((uint32_t)mant << 13);
    memcpy(&fv, &x, 4); return (double)fv;
}

/* Shared helper: f16 → f32 conversion */
static inline float sfa_f16_to_f32(uint16_t h) {
    return (float)sfa_f16_to_f64(h);
}

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

    /* Open source — fix streaming index first */
    sfa_fixup_index(input_path);
    SFA *src = sfa_open(input_path);
    if (!src) { fprintf(stderr, "Cannot open %s\n", input_path); return 1; }

    int src_N = src->Nx;
    long src_N3 = (long)src_N * src_N * src_N;
    long dst_N3 = (long)preview_N * preview_N * preview_N;
    uint32_t n_frames = src->total_frames;
    double src_L = src->Lx;

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

    /* Helper: extract column c from frame_buf as float value at index i */
    /* Column offsets (in bytes) within the frame buffer */
    uint64_t col_offsets[64];
    { uint64_t off = 0;
      for (uint32_t c = 0; c < src->n_columns; c++) {
          col_offsets[c] = off;
          off += (uint64_t)src_N3 * sfa_dtype_size[src->columns[c].dtype];
      }
    }

    /* Read a float from column c at voxel i, handling f16/f32/f64 */
    #define COL_F(buf, c, i) ( \
        (src->columns[c].dtype == SFA_F16) ? sfa_f16_to_f32(((uint16_t*)((uint8_t*)(buf) + col_offsets[c]))[i]) : \
        (src->columns[c].dtype == SFA_F32) ? ((float*)((uint8_t*)(buf) + col_offsets[c]))[i] : \
        (float)((double*)((uint8_t*)(buf) + col_offsets[c]))[i] )

    int scan_frames[] = {0, (int)(n_frames / 4), (int)(n_frames / 2),
                         (int)(3 * n_frames / 4), (int)(n_frames - 1)};
    int n_scan = 5;
    for (int s = 0; s < n_scan; s++) {
        int fi = scan_frames[s];
        if (fi < 0 || fi >= (int)n_frames) continue;
        if (sfa_read_frame(src, fi, frame_buf) < 0) continue;

        for (long i = 0; i < src_N3; i++) {
            float p0 = COL_F(frame_buf, 0, i);
            float p1 = COL_F(frame_buf, 1, i);
            float p2 = COL_F(frame_buf, 2, i);
            float t0 = COL_F(frame_buf, 3, i);
            float t1 = COL_F(frame_buf, 4, i);
            float t2 = COL_F(frame_buf, 5, i);
            float P = fabsf(p0 * p1 * p2);
            float ps = p0*p0 + p1*p1 + p2*p2;
            float ts = t0*t0 + t1*t1 + t2*t2;
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

        /* Compute derived at full resolution (handles f16/f32/f64) */
        #pragma omp parallel for schedule(static)
        for (long i = 0; i < src_N3; i++) {
            float p0 = COL_F(frame_buf, 0, i);
            float p1 = COL_F(frame_buf, 1, i);
            float p2 = COL_F(frame_buf, 2, i);
            float t0 = COL_F(frame_buf, 3, i);
            float t1 = COL_F(frame_buf, 4, i);
            float t2 = COL_F(frame_buf, 5, i);
            derived[0][i] = fabsf(p0 * p1 * p2);
            derived[1][i] = p0*p0 + p1*p1 + p2*p2;
            derived[2][i] = t0*t0 + t1*t1 + t2*t2;
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
   info — Print SFA file properties
   ================================================================ */

static const char *codec_name(uint32_t flags) {
    switch (flags & 0xF) {
    case 0: return "raw"; case 1: return "zstd"; case 2: return "BSS+zstd";
    case 3: return "f32+BSS+zstd"; case 4: return "f16+BSS+zstd";
    case 5: return "bq8+zstd"; case 6: return "idelta+BSS+zstd";
    case 7: return "colzstd"; default: return "unknown";
    }
}

static int cmd_info(int argc, char **argv) {
    if (argc < 1) { fprintf(stderr, "Usage: sfa info <file.sfa> [--json]\n"); return 1; }
    const char *path = argv[0];
    int json = (argc > 1 && !strcmp(argv[1], "--json"));

    struct stat st;
    if (stat(path, &st) != 0) { fprintf(stderr, "Cannot stat %s\n", path); return 1; }
    uint64_t file_bytes = st.st_size;

    SFA *s = sfa_open(path);
    if (!s) { fprintf(stderr, "Cannot open %s\n", path); return 1; }

    long N3 = (long)s->Nx * s->Ny * s->Nz;
    double dx = (s->Nx > 1) ? 2.0 * s->Lx / (s->Nx - 1) : 0;
    double dy = (s->Ny > 1) ? 2.0 * s->Ly / (s->Ny - 1) : 0;
    double dz = (s->Nz > 1) ? 2.0 * s->Lz / (s->Nz - 1) : 0;
    int is_cubic = (s->Nx == s->Ny && s->Ny == s->Nz);

    int has_vel=0;
    int col_dtype=-1, mixed_dtype=0, bytes_per_voxel=0;
    for (uint32_t c = 0; c < s->n_columns; c++) {
        if (col_dtype < 0) col_dtype = s->columns[c].dtype;
        else if (s->columns[c].dtype != col_dtype) mixed_dtype = 1;
        bytes_per_voxel += sfa_dtype_size[s->columns[c].dtype];
        if (s->columns[c].semantic == SFA_VELOCITY) has_vel = 1;
    }

    uint64_t frame_raw = (uint64_t)N3 * bytes_per_voxel;
    int n_valid = sfa_count_valid_frames(s, file_bytes);
    int truncated = (n_valid < (int)s->total_frames);

    double t_first=0, t_last=0;
    if (n_valid > 0) {
        t_first = sfa_frame_time(s, 0);
        t_last = sfa_frame_time(s, n_valid - 1);
    }

    SFA_KVMDSet kv[16];
    int n_kv = sfa_read_kvmd(s, kv, 16);

    if (json) {
        printf("{\"path\":\"%s\",\"size\":%lu,\"grid\":[%u,%u,%u],\"L\":[%.1f,%.1f,%.1f],",
               path, (unsigned long)file_bytes, s->Nx, s->Ny, s->Nz, s->Lx, s->Ly, s->Lz);
        printf("\"codec\":\"%s\",\"columns\":%u,\"dtype\":\"%s\",",
               codec_name(s->flags), s->n_columns, mixed_dtype?"mixed":sfa_dtype_name(col_dtype));
        printf("\"frames\":%u,\"valid\":%d,\"t_range\":[%.2f,%.2f]}\n",
               s->total_frames, n_valid, t_first, t_last);
    } else {
        printf("%-20s %s\n", "File:", path);
        printf("%-20s %.1f %s\n", "Size:",
               file_bytes>1e9?file_bytes/1e9:file_bytes>1e6?file_bytes/1e6:file_bytes/1e3,
               file_bytes>1e9?"GB":file_bytes>1e6?"MB":"KB");
        printf("%-20s v%u, %s%s\n", "Format:", s->version, codec_name(s->flags),
               (s->flags & SFA_FLAG_STREAMING) ? " [STREAMING]" : "");
        printf("%-20s %ux%ux%u = %ld voxels%s\n", "Grid:",
               s->Nx, s->Ny, s->Nz, N3, is_cubic?" (cubic)":"");
        printf("%-20s [%.1f, %.1f, %.1f]\n", "Domain (L):", s->Lx, s->Ly, s->Lz);
        printf("%-20s [%.4f, %.4f, %.4f]\n", "Resolution (dx):", dx, dy, dz);
        printf("%-20s %.10f\n", "Timestep (dt):", s->dt);
        printf("\n%-20s %u (%s%s)\n", "Columns:", s->n_columns,
               mixed_dtype?"mixed":sfa_dtype_name(col_dtype), has_vel?", restartable":"");
        for (uint32_t c = 0; c < s->n_columns; c++) {
            const char *sem_names[]={"pos","ang","vel","acc","nrg","bnd","tor","met","msk"};
            const char *sem=(s->columns[c].semantic<9)?sem_names[s->columns[c].semantic]:"cst";
            printf("  [%u] %-12s %s  %s/%d\n", c, s->columns[c].name,
                   sfa_dtype_name(s->columns[c].dtype), sem, s->columns[c].component);
        }
        printf("\n%-20s %u", "Frames:", s->total_frames);
        if (truncated) printf("  *** %d valid, %d bad ***", n_valid, (int)s->total_frames-n_valid);
        printf("\n");
        if (n_valid > 0)
            printf("%-20s %.2f → %.2f (duration %.2f)\n", "Time range:",
                   t_first, t_last, t_last - t_first);
        printf("%-20s %.2f:1\n", "Compression:", (double)frame_raw*n_valid/file_bytes);
        if (truncated)
            printf("\n!!! WARNING: FILE TRUNCATED — %d of %u frames valid !!!\n",
                   n_valid, s->total_frames);
        if (n_kv > 0) {
            printf("\n%-20s %d set(s)\n", "KVMD metadata:", n_kv);
            for (int i=0;i<n_kv;i++)
                for (int p=0;p<kv[i].n_pairs;p++)
                    printf("  %-16s = %s\n", kv[i].keys[p], kv[i].values[p]);
        }
    }
    sfa_close(s);
    return 0;
}

/* ================================================================
   extract — Extract a single frame to a new SFA seed
   ================================================================ */

static int cmd_extract(int argc, char **argv) {
    if (argc < 1) {
        fprintf(stderr, "Usage: sfa extract <input.sfa> -o <output.sfa> [--frame N]\n"
                "  Extracts frame N (default: last valid) to a single-frame SFA.\n");
        return 1;
    }
    const char *in_path = argv[0], *out_path = "extracted.sfa";
    int frame_idx = -1;
    for (int i=1;i<argc;i++) {
        if (!strcmp(argv[i],"-o") && i+1<argc) out_path=argv[++i];
        else if (!strcmp(argv[i],"--frame") && i+1<argc) frame_idx=atoi(argv[++i]);
    }

    struct stat st;
    if (stat(in_path, &st)!=0) { fprintf(stderr,"Cannot stat %s\n",in_path); return 1; }
    SFA *src = sfa_open(in_path);
    if (!src) { fprintf(stderr,"Cannot open %s\n",in_path); return 1; }

    int valid = sfa_count_valid_frames(src, st.st_size);
    if (valid==0) { fprintf(stderr,"No valid frames\n"); sfa_close(src); return 1; }
    if (frame_idx<0) frame_idx = valid-1;
    if (frame_idx>=valid) { frame_idx=valid-1; printf("Using last valid frame %d\n",frame_idx); }

    double t = sfa_frame_time(src, frame_idx);
    printf("Extracting frame %d (t=%.2f) from %s\n", frame_idx, t, in_path);

    void *buf = malloc(src->frame_bytes);
    if (!buf) { fprintf(stderr, "Cannot allocate %lu bytes\n", (unsigned long)src->frame_bytes); sfa_close(src); return 1; }
    if (sfa_read_frame(src, frame_idx, buf) != 0) {
        fprintf(stderr, "Failed to read frame %d\n", frame_idx); free(buf); sfa_close(src); return 1;
    }

    SFA *dst = sfa_create(out_path, src->Nx, src->Ny, src->Nz, src->Lx, src->Ly, src->Lz, src->dt);
    dst->n_columns = src->n_columns;
    for (uint32_t c=0;c<src->n_columns;c++) dst->columns[c] = src->columns[c];

    /* Copy KVMD */
    SFA_KVMDSet kv[SFA_MAX_KVMD_SETS];
    int n_kv = sfa_read_kvmd(src, kv, SFA_MAX_KVMD_SETS);
    for (int i=0;i<n_kv;i++) {
        const char *kp[128], *vp[128];
        for (int j=0;j<kv[i].n_pairs;j++) { kp[j]=kv[i].keys[j]; vp[j]=kv[i].values[j]; }
        sfa_add_kvmd(dst, kv[i].set_id, 0xFFFFFFFF, 0, kp, vp, kv[i].n_pairs);
    }
    { char fb[32],tb[32]; snprintf(fb,32,"%d",frame_idx); snprintf(tb,32,"%.6f",t);
      const char *k[]={"source_file","source_frame","source_time"};
      const char *v[]={in_path,fb,tb};
      sfa_add_kvmd(dst, 100, 0xFFFFFFFF, 0, k, v, 3);
    }
    sfa_finalize_header(dst);

    long N3 = (long)src->Nx*src->Ny*src->Nz;
    void **cols = malloc(src->n_columns*sizeof(void*));
    uint64_t off=0;
    for (uint32_t c=0;c<src->n_columns;c++) {
        cols[c]=(uint8_t*)buf+off;
        off+=(uint64_t)N3*sfa_dtype_size[src->columns[c].dtype];
    }
    sfa_write_frame(dst, t, cols);
    printf("Wrote %s: 1 frame at t=%.2f\n", out_path, t);

    free(cols); free(buf); sfa_close(src); sfa_close(dst);
    return 0;
}

/* ================================================================
   repair — Fix truncated SFA files
   ================================================================ */

#define SFAH_TOTAL_FRAMES_OFFSET 68

static int cmd_repair(int argc, char **argv) {
    if (argc < 1) { fprintf(stderr, "Usage: sfa repair <file.sfa>\n"); return 1; }
    const char *path = argv[0];

    struct stat st;
    if (stat(path,&st)!=0) { fprintf(stderr,"Cannot stat %s\n",path); return 1; }

    printf("Running index fixup scan...\n");
    int fixup = sfa_fixup_index(path);
    if (fixup>0) printf("  Recovered %d frame(s) from FRMD scan.\n", fixup);
    else printf("  No index fixup needed.\n");

    SFA *s = sfa_open(path);
    if (!s) { fprintf(stderr,"Cannot open %s\n",path); return 1; }

    int valid = sfa_count_valid_frames(s, st.st_size);
    int total = (int)s->total_frames;
    if (valid==total && total>0) {
        printf("File OK: all %d frames valid.\n", total);
        sfa_close(s); return 0;
    }
    if (valid==0) { fprintf(stderr,"No valid frames.\n"); sfa_close(s); return 1; }

    printf("Repairing: %d valid of %d total.\n", valid, total);
    double t_last = sfa_frame_time(s, valid-1);
    sfa_close(s);

    FILE *fp = fopen(path, "r+b");
    if (!fp) { fprintf(stderr, "Cannot open %s for writing\n", path); return 1; }
    uint32_t nv = (uint32_t)valid;
    fseek(fp, SFAH_TOTAL_FRAMES_OFFSET, SEEK_SET);
    if (fwrite(&nv, 4, 1, fp) != 1) {
        fprintf(stderr, "Failed to write total_frames\n"); fclose(fp); return 1;
    }

    /* Zero invalid JMPF entries */
    SFA *s2 = sfa_open(path);
    if (s2) {
        fseek(s2->fp, (long)s2->first_jtop_offset+12+4+4+8, SEEK_SET);
        uint64_t jmpf_off = 0;
        if (fread(&jmpf_off, 8, 1, s2->fp) != 1) {
            fprintf(stderr, "Warning: cannot read JMPF offset\n");
        } else {
            for (int f=valid;f<total;f++) {
                long eoff = (long)jmpf_off+12+8+(long)f*32;
                fseek(fp, eoff, SEEK_SET);
                char zeros[32]={0};
                if (fwrite(zeros,1,32,fp) != 32) {
                    fprintf(stderr, "Warning: failed to zero JMPF entry %d\n", f);
                    break;
                }
            }
        }
        sfa_close(s2);
    }
    fclose(fp);
    printf("Repaired: %d frames (t=0 to %.2f).\n", valid, t_last);
    return 0;
}

/* ================================================================
   convert — Convert between codecs (uncompress, f16→f32, recompress)
   ================================================================ */

static int cmd_convert(int argc, char **argv) {
    if (argc < 1) {
        fprintf(stderr, "Usage: sfa convert <input.sfa> -o <output.sfa> [--f32] [--colzstd]\n"
                "  --f32:      promote f16 columns to f32\n"
                "  --colzstd:  per-column zstd compression (default: raw/mmap)\n");
        return 1;
    }
    const char *in_path=argv[0], *out_path=NULL;
    int promote=0, colzstd=0;
    for (int i=1;i<argc;i++) {
        if (!strcmp(argv[i],"-o") && i+1<argc) out_path=argv[++i];
        else if (!strcmp(argv[i],"--f32")) promote=1;
        else if (!strcmp(argv[i],"--colzstd")) colzstd=1;
    }
    if (!out_path) { fprintf(stderr,"Missing -o output.sfa\n"); return 1; }

    struct stat st;
    if (stat(in_path,&st)!=0) { fprintf(stderr,"Cannot stat %s\n",in_path); return 1; }
    SFA *src = sfa_open(in_path);
    if (!src) { fprintf(stderr,"Cannot open %s\n",in_path); return 1; }

    int valid = sfa_count_valid_frames(src, st.st_size);
    if (valid==0) { fprintf(stderr,"No valid frames\n"); sfa_close(src); return 1; }
    long N3 = (long)src->Nx*src->Ny*src->Nz;

    printf("Converting: %s → %s\n", in_path, out_path);
    printf("  %ux%ux%u, %d frames, %s%s\n", src->Nx, src->Ny, src->Nz, valid,
           colzstd?"colzstd":"raw", promote?" +f16→f32":"");

    /* Determine output dtypes */
    uint8_t out_dt[64]; int out_es[64]; uint64_t out_fb=0;
    for (uint32_t c=0;c<src->n_columns;c++) {
        uint8_t d=src->columns[c].dtype;
        if (promote && d==SFA_F16) { out_dt[c]=SFA_F32; out_es[c]=4; }
        else { out_dt[c]=d; out_es[c]=sfa_dtype_size[d]; }
        out_fb+=(uint64_t)N3*out_es[c];
    }

    SFA *dst = sfa_create(out_path, src->Nx, src->Ny, src->Nz, src->Lx, src->Ly, src->Lz, src->dt);
    dst->flags = colzstd ? SFA_CODEC_COLZSTD : SFA_CODEC_RAW;
    dst->n_columns = src->n_columns;
    for (uint32_t c=0;c<src->n_columns;c++) { dst->columns[c]=src->columns[c]; dst->columns[c].dtype=out_dt[c]; }

    SFA_KVMDSet kv[SFA_MAX_KVMD_SETS];
    int n_kv = sfa_read_kvmd(src, kv, SFA_MAX_KVMD_SETS);
    for (int i=0;i<n_kv;i++) {
        const char *kp[128],*vp[128];
        for (int j=0;j<kv[i].n_pairs;j++) { kp[j]=kv[i].keys[j]; vp[j]=kv[i].values[j]; }
        sfa_add_kvmd(dst, kv[i].set_id, 0xFFFFFFFF, 0, kp, vp, kv[i].n_pairs);
    }
    dst->total_frames = valid;
    sfa_finalize_header(dst);

    void *sbuf = malloc(src->frame_bytes);
    uint8_t *obuf = malloc(out_fb);
    void **ocols = malloc(src->n_columns*sizeof(void*));

    for (int f=0;f<valid;f++) {
        double t = sfa_frame_time(src, f);
        if (sfa_read_frame(src,f,sbuf)!=0) { fprintf(stderr,"Read error frame %d\n",f); break; }

        uint64_t soff=0, doff=0;
        for (uint32_t c=0;c<src->n_columns;c++) {
            uint8_t sd=src->columns[c].dtype;
            int ses=sfa_dtype_size[sd];
            uint8_t *sc=(uint8_t*)sbuf+soff, *dc=obuf+doff;
            if (promote && sd==SFA_F16) {
                uint16_t *in16=(uint16_t*)sc; float *out32=(float*)dc;
                for (long i=0;i<N3;i++) out32[i]=sfa_f16_to_f32(in16[i]);
            } else memcpy(dc,sc,(uint64_t)N3*ses);
            ocols[c]=dc; soff+=(uint64_t)N3*ses; doff+=(uint64_t)N3*out_es[c];
        }
        sfa_write_frame_ex(dst, t, ocols, 0, out_fb);
        if ((f+1)%10==0 || f==valid-1) { printf("\r  %d/%d",f+1,valid); fflush(stdout); }
    }
    printf("\n");

    sfa_close(dst); sfa_close(src);
    struct stat ost; stat(out_path,&ost);
    printf("Done: %s (%.1f GB)\n", out_path, ost.st_size/1e9);
    free(sbuf); free(obuf); free(ocols);
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
        "  info      Print SFA file metadata and diagnostics\n"
        "  extract   Extract a single frame to a new SFA seed\n"
        "  repair    Fix truncated SFA files (remove bad frames)\n"
        "  convert   Convert between codecs (uncompress, f16→f32)\n"
        "  preview   Create uint8 preview SFA for fast viewing\n"
        "\n"
        "Run 'sfa <command> --help' for command-specific usage.\n"
        "\n");
}

int main(int argc, char **argv) {
    if (argc < 2) { usage_main(); return 1; }

    const char *cmd = argv[1];

    if (!strcmp(cmd, "info"))       return cmd_info(argc-2, argv+2);
    if (!strcmp(cmd, "extract"))    return cmd_extract(argc-2, argv+2);
    if (!strcmp(cmd, "repair"))     return cmd_repair(argc-2, argv+2);
    if (!strcmp(cmd, "convert"))    return cmd_convert(argc-2, argv+2);
    if (!strcmp(cmd, "preview"))    return cmd_preview(argc-2, argv+2);
    if (!strcmp(cmd, "help") || !strcmp(cmd, "--help")) { usage_main(); return 0; }

    fprintf(stderr, "Unknown command: %s\n\n", cmd);
    usage_main();
    return 1;
}
