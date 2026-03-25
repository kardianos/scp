/*  sfa_info.c — Print SFA file properties
 *
 *  Quick summary of an SFA file: grid, columns, frames, sizes, physics.
 *
 *  Build: cd sfa/analysis && make sfa_info
 *  Usage: ./sfa_info file.sfa [--json]
 */

#define SFA_IMPLEMENTATION
#include "../format/sfa.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>

static const char *codec_name(uint32_t flags) {
    switch (flags & 0xF) {
    case 0: return "raw";
    case 1: return "zstd";
    case 2: return "BSS+zstd";
    case 3: return "f32+BSS+zstd";
    case 4: return "f16+BSS+zstd";
    case 5: return "bq8+zstd";
    case 6: return "idelta+BSS+zstd";
    default: return "unknown";
    }
}

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "Usage: %s file.sfa [--json]\n", argv[0]); return 1; }
    const char *path = argv[1];
    int json = (argc > 2 && !strcmp(argv[2], "--json"));

    /* File size */
    struct stat st;
    if (stat(path, &st) != 0) { fprintf(stderr, "Cannot stat %s\n", path); return 1; }
    uint64_t file_bytes = st.st_size;

    SFA *s = sfa_open(path);
    if (!s) { fprintf(stderr, "Cannot open %s\n", path); return 1; }

    /* Derived quantities */
    long N3 = (long)s->Nx * s->Ny * s->Nz;
    double dx = (s->Nx > 1) ? 2.0 * s->Lx / (s->Nx - 1) : 0;
    double dy = (s->Ny > 1) ? 2.0 * s->Ly / (s->Ny - 1) : 0;
    double dz = (s->Nz > 1) ? 2.0 * s->Lz / (s->Nz - 1) : 0;
    int is_cubic = (s->Nx == s->Ny && s->Ny == s->Nz);
    int is_isotropic = (fabs(s->Lx - s->Ly) < 1e-6 && fabs(s->Ly - s->Lz) < 1e-6);

    /* Column dtype summary */
    int has_phi = 0, has_theta = 0, has_vel = 0, has_accel = 0;
    int col_dtype = -1;
    int mixed_dtype = 0;
    int bytes_per_voxel = 0;
    for (uint32_t c = 0; c < s->n_columns; c++) {
        if (col_dtype < 0) col_dtype = s->columns[c].dtype;
        else if (s->columns[c].dtype != col_dtype) mixed_dtype = 1;
        bytes_per_voxel += sfa_dtype_size[s->columns[c].dtype];
        if (s->columns[c].semantic == SFA_POSITION) has_phi = 1;
        if (s->columns[c].semantic == SFA_ANGLE) has_theta = 1;
        if (s->columns[c].semantic == SFA_VELOCITY) has_vel = 1;
        if (s->columns[c].semantic == SFA_ACCELERATION) has_accel = 1;
    }

    uint64_t frame_raw = (uint64_t)N3 * bytes_per_voxel;
    double compression = (s->total_frames > 0 && frame_raw > 0) ?
        (double)file_bytes / ((double)s->total_frames * frame_raw) : 0;

    /* Time range */
    double t_first = 0, t_last = 0, t_step = 0;
    if (s->total_frames > 0) {
        t_first = sfa_frame_time(s, 0);
        t_last = sfa_frame_time(s, s->total_frames - 1);
        if (s->total_frames > 1)
            t_step = (t_last - t_first) / (s->total_frames - 1);
    }

    /* Frame integrity check: verify each frame's data fits within the file */
    int n_valid_frames = 0, n_truncated_frames = 0, n_missing_frames = 0;
    int truncation_warning = 0;
    uint32_t first_bad_frame = s->total_frames;
    if (s->total_frames > 0 && s->fp) {
        /* Walk JTOP → JMPF to read per-frame offsets and sizes */
        FILE *fp = s->fp;
        long saved_pos = ftell(fp);

        /* Find JMPF via JTOP */
        fseek(fp, (long)s->first_jtop_offset + 12 + 4 + 4 + 8, SEEK_SET);
        uint64_t jmpf_off;
        fread(&jmpf_off, 8, 1, fp);

        /* Read JMPF entries */
        fseek(fp, (long)jmpf_off + 12 + 4 + 4, SEEK_SET);

        for (uint32_t f = 0; f < s->total_frames; f++) {
            double ftime;
            uint64_t foffset, fcomp_size;
            uint32_t fcrc, freserved;
            fread(&ftime, 8, 1, fp);
            fread(&foffset, 8, 1, fp);
            fread(&fcomp_size, 8, 1, fp);
            fread(&fcrc, 4, 1, fp);
            fread(&freserved, 4, 1, fp);

            if (foffset == 0 && fcomp_size == 0 && f > 0) {
                /* Empty JMPF slot — frame was never written */
                n_missing_frames++;
                if (f < first_bad_frame) first_bad_frame = f;
            } else {
                /* Check if the frame data fits within the file */
                /* Frame chunk = 12 bytes header + comp_size bytes data */
                uint64_t frame_end = foffset + 12 + fcomp_size;
                if (frame_end > file_bytes) {
                    n_truncated_frames++;
                    if (f < first_bad_frame) first_bad_frame = f;
                } else {
                    n_valid_frames++;
                }
            }
        }

        truncation_warning = (n_truncated_frames > 0 || n_missing_frames > 0);

        /* Recalculate compression using only valid frames */
        if (n_valid_frames > 0 && n_valid_frames != (int)s->total_frames) {
            compression = (n_valid_frames > 0 && frame_raw > 0) ?
                (double)file_bytes / ((double)n_valid_frames * frame_raw) : 0;
        }

        /* Update time range to reflect only valid frames */
        if (truncation_warning && first_bad_frame > 0) {
            /* Re-read the last valid frame's time */
            fseek(fp, (long)jmpf_off + 12 + 8 + (first_bad_frame - 1) * 32, SEEK_SET);
            fread(&t_last, 8, 1, fp);
            if (first_bad_frame > 1)
                t_step = (t_last - t_first) / (first_bad_frame - 1);
        }

        fseek(fp, saved_pos, SEEK_SET);
    }

    /* KVMD */
    SFA_KVMDSet kv[16];
    int n_kv = sfa_read_kvmd(s, kv, 16);
    int n_kv_pairs = 0;
    for (int i = 0; i < n_kv; i++) n_kv_pairs += kv[i].n_pairs;

    /* GDEF */
    SFA_GridDef grids[16];
    int n_grids = sfa_read_grids(s, grids, 16);

    /* Physics from KVMD */
    const char *kv_m = NULL, *kv_eta = NULL, *kv_mu = NULL, *kv_kappa = NULL;
    const char *kv_chirality = NULL, *kv_type = NULL;
    for (int i = 0; i < n_kv; i++) for (int p = 0; p < kv[i].n_pairs; p++) {
        if (!strcmp(kv[i].keys[p], "m")) kv_m = kv[i].values[p];
        if (!strcmp(kv[i].keys[p], "eta")) kv_eta = kv[i].values[p];
        if (!strcmp(kv[i].keys[p], "mu")) kv_mu = kv[i].values[p];
        if (!strcmp(kv[i].keys[p], "kappa")) kv_kappa = kv[i].values[p];
        if (!strcmp(kv[i].keys[p], "chirality")) kv_chirality = kv[i].values[p];
        if (!strcmp(kv[i].keys[p], "type")) kv_type = kv[i].values[p];
    }

    if (json) {
        printf("{\n");
        printf("  \"path\": \"%s\",\n", path);
        printf("  \"file_size_bytes\": %lu,\n", (unsigned long)file_bytes);
        printf("  \"file_size_human\": \"%.1f %s\",\n",
               file_bytes > 1e9 ? file_bytes/1e9 : file_bytes > 1e6 ? file_bytes/1e6 : file_bytes/1e3,
               file_bytes > 1e9 ? "GB" : file_bytes > 1e6 ? "MB" : "KB");
        printf("  \"version\": %u,\n", s->version);
        printf("  \"codec\": \"%s\",\n", codec_name(s->flags));
        printf("  \"streaming\": %s,\n", (s->flags & SFA_FLAG_STREAMING) ? "true" : "false");
        printf("  \"Nx\": %u, \"Ny\": %u, \"Nz\": %u,\n", s->Nx, s->Ny, s->Nz);
        printf("  \"N_total\": %ld,\n", N3);
        printf("  \"cubic\": %s,\n", is_cubic ? "true" : "false");
        printf("  \"Lx\": %.6f, \"Ly\": %.6f, \"Lz\": %.6f,\n", s->Lx, s->Ly, s->Lz);
        printf("  \"isotropic\": %s,\n", is_isotropic ? "true" : "false");
        printf("  \"dx\": %.6f, \"dy\": %.6f, \"dz\": %.6f,\n", dx, dy, dz);
        printf("  \"dt\": %.10f,\n", s->dt);
        printf("  \"n_columns\": %u,\n", s->n_columns);
        printf("  \"dtype\": \"%s\",\n", mixed_dtype ? "mixed" : sfa_dtype_name(col_dtype));
        printf("  \"bytes_per_voxel\": %d,\n", bytes_per_voxel);
        printf("  \"has_phi\": %s,\n", has_phi ? "true" : "false");
        printf("  \"has_theta\": %s,\n", has_theta ? "true" : "false");
        printf("  \"has_velocity\": %s,\n", has_vel ? "true" : "false");
        printf("  \"has_acceleration\": %s,\n", has_accel ? "true" : "false");
        printf("  \"restartable\": %s,\n", (has_phi && has_theta && has_vel) ? "true" : "false");
        printf("  \"total_frames\": %u,\n", s->total_frames);
        printf("  \"valid_frames\": %d,\n", n_valid_frames);
        printf("  \"truncated_frames\": %d,\n", n_truncated_frames);
        printf("  \"missing_frames\": %d,\n", n_missing_frames);
        printf("  \"integrity\": \"%s\",\n", truncation_warning ? "TRUNCATED" : "OK");
        printf("  \"t_first\": %.2f,\n", t_first);
        printf("  \"t_last\": %.2f,\n", t_last);
        printf("  \"t_step\": %.2f,\n", t_step);
        printf("  \"t_duration\": %.2f,\n", t_last - t_first);
        printf("  \"frame_raw_bytes\": %lu,\n", (unsigned long)frame_raw);
        printf("  \"frame_raw_human\": \"%.1f %s\",\n",
               frame_raw > 1e9 ? frame_raw/1e9 : frame_raw > 1e6 ? frame_raw/1e6 : frame_raw/1e3,
               frame_raw > 1e9 ? "GB" : frame_raw > 1e6 ? "MB" : "KB");
        printf("  \"compression_ratio\": %.2f,\n", compression);
        printf("  \"n_kvmd_sets\": %d,\n", n_kv);
        printf("  \"n_kvmd_pairs\": %d,\n", n_kv_pairs);
        printf("  \"n_grids\": %d,\n", n_grids);
        if (kv_m) printf("  \"physics_m\": \"%s\",\n", kv_m);
        if (kv_eta) printf("  \"physics_eta\": \"%s\",\n", kv_eta);
        if (kv_mu) printf("  \"physics_mu\": \"%s\",\n", kv_mu);
        if (kv_kappa) printf("  \"physics_kappa\": \"%s\",\n", kv_kappa);
        if (kv_chirality) printf("  \"chirality\": \"%s\",\n", kv_chirality);
        if (kv_type) printf("  \"type\": \"%s\",\n", kv_type);
        printf("  \"columns\": [\n");
        for (uint32_t c = 0; c < s->n_columns; c++) {
            const char *sem_names[] = {"position","angle","velocity","acceleration",
                                       "energy","binding","torsion","metric","mask"};
            const char *sem = (s->columns[c].semantic < 9) ? sem_names[s->columns[c].semantic] : "custom";
            printf("    {\"name\":\"%s\",\"dtype\":\"%s\",\"semantic\":\"%s\",\"component\":%d}%s\n",
                   s->columns[c].name, sfa_dtype_name(s->columns[c].dtype),
                   sem, s->columns[c].component, c < s->n_columns-1 ? "," : "");
        }
        printf("  ]\n}\n");
    } else {
        printf("%-20s %s\n", "File:", path);
        printf("%-20s %lu bytes (%.1f %s)\n", "Size:",
               (unsigned long)file_bytes,
               file_bytes > 1e9 ? file_bytes/1e9 : file_bytes > 1e6 ? file_bytes/1e6 : file_bytes/1e3,
               file_bytes > 1e9 ? "GB" : file_bytes > 1e6 ? "MB" : "KB");
        printf("%-20s v%u, %s%s\n", "Format:",
               s->version, codec_name(s->flags),
               (s->flags & SFA_FLAG_STREAMING) ? " [STREAMING]" : "");
        printf("\n");
        printf("%-20s %ux%ux%u = %ld voxels%s\n", "Grid:",
               s->Nx, s->Ny, s->Nz, N3, is_cubic ? " (cubic)" : "");
        printf("%-20s [%.1f, %.1f, %.1f]%s\n", "Domain (L):",
               s->Lx, s->Ly, s->Lz, is_isotropic ? " (isotropic)" : "");
        printf("%-20s [%.4f, %.4f, %.4f]\n", "Resolution (dx):", dx, dy, dz);
        printf("%-20s %.10f\n", "Timestep (dt):", s->dt);
        printf("\n");
        printf("%-20s %u (%s%s)\n", "Columns:",
               s->n_columns,
               mixed_dtype ? "mixed" : sfa_dtype_name(col_dtype),
               has_vel ? ", restartable" : "");
        for (uint32_t c = 0; c < s->n_columns; c++) {
            const char *sem_names[] = {"pos","ang","vel","acc","nrg","bnd","tor","met","msk"};
            const char *sem = (s->columns[c].semantic < 9) ? sem_names[s->columns[c].semantic] : "cst";
            printf("  [%u] %-12s %s  %s/%d\n", c, s->columns[c].name,
                   sfa_dtype_name(s->columns[c].dtype), sem, s->columns[c].component);
        }
        printf("\n");
        printf("%-20s %u", "Frames:", s->total_frames);
        if (truncation_warning)
            printf("  *** %d valid, %d truncated, %d missing ***",
                   n_valid_frames, n_truncated_frames, n_missing_frames);
        printf("\n");
        if (s->total_frames > 0) {
            printf("%-20s %.2f → %.2f (duration %.2f)\n", "Time range:",
                   t_first, t_last, t_last - t_first);
            if (n_valid_frames > 1)
                printf("%-20s %.2f (snap_dt)\n", "Frame interval:", t_step);
        }
        if (truncation_warning) {
            printf("\n");
            printf("!!! WARNING: FILE IS TRUNCATED !!!\n");
            printf("  Header claims %u frames but only %d are complete.\n",
                   s->total_frames, n_valid_frames);
            if (n_truncated_frames > 0)
                printf("  %d frame(s) have data extending past end of file (partial download?).\n",
                       n_truncated_frames);
            if (n_missing_frames > 0)
                printf("  %d frame(s) have empty index entries (never written?).\n",
                       n_missing_frames);
            printf("  First bad frame: %u (t=%.2f)\n", first_bad_frame,
                   first_bad_frame < s->total_frames ? sfa_frame_time(s, first_bad_frame) : -1.0);
            printf("  Safe to use frames 0-%u (t=0 to %.2f)\n",
                   first_bad_frame > 0 ? first_bad_frame - 1 : 0, t_last);
            printf("\n");
        }
        printf("%-20s %.1f %s (uncompressed)\n", "Frame size:",
               frame_raw > 1e9 ? frame_raw/1e9 : frame_raw > 1e6 ? frame_raw/1e6 : frame_raw/1e3,
               frame_raw > 1e9 ? "GB" : frame_raw > 1e6 ? "MB" : "KB");
        printf("%-20s %.2f:1\n", "Compression:", 1.0/compression);
        printf("%-20s %d bytes/voxel\n", "Voxel size:", bytes_per_voxel);
        printf("\n");
        printf("%-20s %s%s%s%s\n", "Content:",
               has_phi ? "φ " : "", has_theta ? "θ " : "",
               has_vel ? "v " : "", has_accel ? "a " : "");

        if (n_kv > 0) {
            printf("\n%-20s %d set(s), %d pair(s)\n", "KVMD metadata:", n_kv, n_kv_pairs);
            for (int i = 0; i < n_kv; i++)
                for (int p = 0; p < kv[i].n_pairs; p++)
                    printf("  %-16s = %s\n", kv[i].keys[p], kv[i].values[p]);
        }

        if (n_grids > 0) {
            printf("\n%-20s %d\n", "GDEF grids:", n_grids);
            for (int g = 0; g < n_grids; g++)
                printf("  grid %d: %ux%ux%u L=(%.1f,%.1f,%.1f) center=(%.1f,%.1f,%.1f) ghost=%d\n",
                       grids[g].grid_id, grids[g].Nx, grids[g].Ny, grids[g].Nz,
                       grids[g].Lx, grids[g].Ly, grids[g].Lz,
                       grids[g].cx, grids[g].cy, grids[g].cz, grids[g].ghost);
        }
    }

    sfa_close(s);
    return 0;
}
