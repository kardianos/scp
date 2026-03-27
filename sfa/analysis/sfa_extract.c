/*  sfa_extract.c — SFA frame extraction and truncation repair
 *
 *  Two modes:
 *    (a) repair:  Fix a truncated SFA file in-place — remove bad frames,
 *                 update header to reflect only valid frames.
 *    (b) extract: Copy frame N (or last valid frame) to a new single-frame
 *                 SFA file, preserving KVMD metadata. Useful for creating
 *                 restart seeds from simulation output.
 *
 *  Build: gcc -O3 -o sfa_extract sfa_extract.c -lzstd -lm
 *
 *  Usage:
 *    sfa_extract --repair input.sfa
 *    sfa_extract --extract input.sfa -o seed.sfa [--frame N]
 *    sfa_extract --extract input.sfa -o seed.sfa --frame -1   (last valid frame)
 */

#define SFA_IMPLEMENTATION
#include "../format/sfa.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

static inline double f16_to_f64(uint16_t h) {
    uint16_t sign = h & 0x8000;
    int exp = (h >> 10) & 0x1F;
    uint16_t mant = h & 0x3FF;
    if (exp == 0) return 0.0;
    if (exp == 31) return sign ? -1e30 : 1e30;
    float fv; uint32_t x = ((uint32_t)sign << 16) | ((uint32_t)(exp-15+127) << 23) | ((uint32_t)mant << 13);
    memcpy(&fv, &x, 4); return (double)fv;
}

/* Uses sfa_sfa_count_valid_frames() from sfa.h */

/* ---- REPAIR MODE ---- */
/* total_frames offset in SFAH: 12(chunk hdr) + 4(ver) + 4(flags) + 12(Nx,Ny,Nz)
   + 24(Lx,Ly,Lz) + 8(dt) + 4(n_columns) = 68 */
#define SFAH_TOTAL_FRAMES_OFFSET 68

static int do_repair(const char *path) {
    struct stat st;
    if (stat(path, &st) != 0) { fprintf(stderr, "Cannot stat %s\n", path); return 1; }
    uint64_t file_bytes = st.st_size;

    /* Step 1: Run sfa_fixup_index first to handle streaming files where
     * JMPF entries may be zeroed out (frames written but index not updated).
     * This also handles total_frames=0 from unfinalized streaming writes. */
    printf("Running index fixup scan...\n");
    int fixup_result = sfa_fixup_index(path);
    if (fixup_result > 0)
        printf("  Fixup recovered %d frame(s) from FRMD scan.\n", fixup_result);
    else if (fixup_result == 0)
        printf("  No index fixup needed.\n");
    else
        printf("  Fixup returned error (non-fatal, continuing).\n");

    /* Step 2: Open the (possibly fixed-up) file and count valid frames */
    SFA *s = sfa_open(path);
    if (!s) { fprintf(stderr, "Cannot open %s\n", path); return 1; }

    int valid = sfa_count_valid_frames(s, file_bytes);
    int total = (int)s->total_frames;

    if (valid == total && total > 0) {
        printf("File OK: all %d frames are valid. No repair needed.\n", total);
        sfa_close(s);
        return 0;
    }

    if (valid == 0) {
        fprintf(stderr, "No valid frames found — cannot repair.\n");
        sfa_close(s);
        return 1;
    }

    printf("Repairing: %d valid of %d total frames.\n", valid, total);
    double t_last = sfa_frame_time(s, valid - 1);
    printf("Last valid frame: %d (t=%.2f)\n", valid - 1, t_last);

    sfa_close(s);

    /* Step 3: Update total_frames in SFAH header */
    FILE *fp = fopen(path, "r+b");
    if (!fp) { fprintf(stderr, "Cannot open %s for writing\n", path); return 1; }

    uint32_t new_total = (uint32_t)valid;
    fseek(fp, SFAH_TOTAL_FRAMES_OFFSET, SEEK_SET);
    fwrite(&new_total, 4, 1, fp);
    printf("Updated total_frames: %d → %d\n", total, valid);

    /* Step 4: Zero out invalid JMPF entries */
    SFA *s2 = sfa_open(path);
    if (s2) {
        FILE *fp2 = s2->fp;
        fseek(fp2, (long)s2->first_jtop_offset + 12 + 4 + 4 + 8, SEEK_SET);
        uint64_t jmpf_off;
        fread(&jmpf_off, 8, 1, fp2);

        for (int f = valid; f < total; f++) {
            long entry_off = (long)jmpf_off + 12 + 8 + (long)f * 32;
            fseek(fp, entry_off, SEEK_SET);
            char zeros[32] = {0};
            fwrite(zeros, 1, 32, fp);
        }
        sfa_close(s2);
    }

    fclose(fp);
    printf("Repair complete. File now has %d valid frames (t=0 to %.2f).\n", valid, t_last);
    return 0;
}

/* ---- EXTRACT MODE ---- */
static int do_extract(const char *in_path, const char *out_path, int frame_idx) {
    struct stat st;
    if (stat(in_path, &st) != 0) { fprintf(stderr, "Cannot stat %s\n", in_path); return 1; }
    uint64_t file_bytes = st.st_size;

    SFA *src = sfa_open(in_path);
    if (!src) { fprintf(stderr, "Cannot open %s\n", in_path); return 1; }

    int valid = sfa_count_valid_frames(src, file_bytes);
    if (valid == 0) { fprintf(stderr, "No valid frames in %s\n", in_path); sfa_close(src); return 1; }

    /* Resolve frame index */
    if (frame_idx < 0) frame_idx = valid - 1;
    if (frame_idx >= valid) {
        fprintf(stderr, "Frame %d not valid (only %d valid frames). Using last valid (%d).\n",
                frame_idx, valid, valid - 1);
        frame_idx = valid - 1;
    }

    double t = sfa_frame_time(src, frame_idx);
    printf("Extracting frame %d (t=%.2f) from %s\n", frame_idx, t, in_path);
    printf("  Grid: %ux%ux%u  L=(%.1f,%.1f,%.1f)  %d columns\n",
           src->Nx, src->Ny, src->Nz, src->Lx, src->Ly, src->Lz, src->n_columns);

    /* Read the frame data */
    void *buf = malloc(src->frame_bytes);
    if (!buf) { fprintf(stderr, "Cannot allocate frame buffer (%lu bytes)\n",
                        (unsigned long)src->frame_bytes); sfa_close(src); return 1; }
    if (sfa_read_frame(src, frame_idx, buf) != 0) {
        fprintf(stderr, "Failed to read frame %d\n", frame_idx);
        free(buf); sfa_close(src); return 1;
    }

    /* Create output SFA with same grid/column structure */
    SFA *dst = sfa_create(out_path, src->Nx, src->Ny, src->Nz,
                          src->Lx, src->Ly, src->Lz, src->dt);
    if (!dst) {
        fprintf(stderr, "Cannot create %s\n", out_path);
        free(buf); sfa_close(src); return 1;
    }

    /* Copy column definitions */
    dst->n_columns = src->n_columns;
    for (uint32_t c = 0; c < src->n_columns; c++) {
        dst->columns[c] = src->columns[c];
    }

    /* Copy KVMD metadata from source using sfa_add_kvmd */
    SFA_KVMDSet kv[SFA_MAX_KVMD_SETS];
    int n_kv = sfa_read_kvmd(src, kv, SFA_MAX_KVMD_SETS);
    for (int i = 0; i < n_kv; i++) {
        const char *kptrs[128], *vptrs[128];
        for (int j = 0; j < kv[i].n_pairs; j++) {
            kptrs[j] = kv[i].keys[j];
            vptrs[j] = kv[i].values[j];
        }
        sfa_add_kvmd(dst, kv[i].set_id, 0xFFFFFFFF, 0,
                     kptrs, vptrs, kv[i].n_pairs);
    }

    /* Add extraction metadata */
    {
        char fbuf[32], tbuf[32];
        snprintf(fbuf, sizeof(fbuf), "%d", frame_idx);
        snprintf(tbuf, sizeof(tbuf), "%.6f", t);
        const char *keys[] = {"source_file", "source_frame", "source_time"};
        const char *vals[] = {in_path, fbuf, tbuf};
        sfa_add_kvmd(dst, 100, 0xFFFFFFFF, 0, keys, vals, 3);
    }

    /* Finalize header (writes SFAH, CDEF, KVMD, JTOP, JMPF) */
    sfa_finalize_header(dst);

    /* Decompose frame buffer into per-column arrays for sfa_write_frame */
    long N3 = (long)src->Nx * src->Ny * src->Nz;
    void **cols = (void**)malloc(src->n_columns * sizeof(void*));
    uint64_t col_off = 0;
    for (uint32_t c = 0; c < src->n_columns; c++) {
        int es = sfa_dtype_size[src->columns[c].dtype];
        cols[c] = (uint8_t*)buf + col_off;
        col_off += (uint64_t)N3 * es;
    }

    sfa_write_frame(dst, t, cols);

    printf("Wrote %s: 1 frame at t=%.2f\n", out_path, t);

    /* Print KVMD summary */
    if (n_kv > 0) {
        printf("  KVMD metadata preserved (%d sets, ", n_kv);
        int total_pairs = 0;
        for (int i = 0; i < n_kv; i++) total_pairs += kv[i].n_pairs;
        printf("%d pairs)\n", total_pairs);
    }

    free(cols);
    free(buf);
    sfa_close(src);
    sfa_close(dst);
    return 0;
}

/* ---- MAIN ---- */
int main(int argc, char **argv) {
    if (argc < 3) {
        fprintf(stderr,
            "Usage:\n"
            "  %s --repair input.sfa\n"
            "      Fix truncated SFA: remove bad frames, update header.\n"
            "\n"
            "  %s --extract input.sfa -o output.sfa [--frame N]\n"
            "      Extract frame N (default: last valid) to a new single-frame SFA.\n"
            "      Preserves KVMD metadata. Use --frame -1 for last valid frame.\n"
            "\n", argv[0], argv[0]);
        return 1;
    }

    if (!strcmp(argv[1], "--repair")) {
        return do_repair(argv[2]);
    }

    if (!strcmp(argv[1], "--extract")) {
        const char *in_path = argv[2];
        const char *out_path = "extracted.sfa";
        int frame_idx = -1;  /* default: last valid */

        for (int i = 3; i < argc; i++) {
            if (!strcmp(argv[i], "-o") && i+1 < argc) out_path = argv[++i];
            else if (!strcmp(argv[i], "--frame") && i+1 < argc) frame_idx = atoi(argv[++i]);
        }
        return do_extract(in_path, out_path, frame_idx);
    }

    fprintf(stderr, "Unknown mode: %s (use --repair or --extract)\n", argv[1]);
    return 1;
}
