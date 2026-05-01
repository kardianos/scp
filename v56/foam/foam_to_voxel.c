/*  foam_to_voxel.c — resample per-cell foam fields to a voxel grid SFA
 *
 *  Reads:
 *    foam_mesh.bin (for cell positions, L)
 *    foam_*.fsnp   (one or more per-cell field snapshots)
 *
 *  Writes:
 *    output.sfa    (NxNxN voxel SFA file with phi[3], theta[3] per voxel)
 *
 *  Algorithm: for each voxel, find the nearest cell using a uniform spatial
 *  grid index over cells. The voxel inherits that cell's field values
 *  (nearest-neighbor / Voronoi reconstruction). Each voxel sits in exactly
 *  one cell — that's the natural Voronoi rendering.
 *
 *  Build:
 *      gcc -O3 -fopenmp -o foam_to_voxel foam_to_voxel.c \
 *          -I../../sfa/format -lzstd -lm
 *
 *  Usage:
 *      ./foam_to_voxel foam_mesh.bin output.sfa N [snapshot1.fsnp ...]
 */

#define _GNU_SOURCE
#define SFA_IMPLEMENTATION
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <omp.h>
#include "../../sfa/format/sfa.h"

/* =====================================================================
 *  Mesh and grid index
 * ===================================================================== */

typedef struct {
    double L;
    uint32_t N_cells;
    double *cx, *cy, *cz;
    /* Grid index: each grid bin holds list of cell IDs */
    int Ng;
    double bin_size;
    int *bin_count;       /* Ng^3 ints */
    int *bin_offset;      /* Ng^3 ints */
    int *bin_cells;       /* sum(bin_count) ints, sorted by bin */
} FoamGrid;

static FoamGrid *foam_load(const char *path) {
    FILE *fp = fopen(path, "rb");
    if (!fp) { fprintf(stderr, "FATAL: cannot open mesh '%s'\n", path); exit(1); }
    char magic[4]; if (fread(magic, 1, 4, fp) != 4) exit(1);
    if (memcmp(magic, "FOAM", 4) != 0) { fprintf(stderr, "bad magic\n"); exit(1); }
    uint32_t version; if (fread(&version, 4, 1, fp) != 1) exit(1);
    FoamGrid *g = calloc(1, sizeof(FoamGrid));
    if (fread(&g->L, 8, 1, fp) != 1) exit(1);
    if (fread(&g->N_cells, 4, 1, fp) != 1) exit(1);
    uint32_t N_faces; if (fread(&N_faces, 4, 1, fp) != 1) exit(1);
    uint64_t reserved; if (fread(&reserved, 8, 1, fp) != 1) exit(1);

    g->cx = malloc(sizeof(double) * g->N_cells);
    g->cy = malloc(sizeof(double) * g->N_cells);
    g->cz = malloc(sizeof(double) * g->N_cells);
    for (uint32_t i = 0; i < g->N_cells; i++) {
        double rec[4];
        if (fread(rec, 8, 4, fp) != 4) exit(1);
        g->cx[i] = rec[0]; g->cy[i] = rec[1]; g->cz[i] = rec[2];
        /* skip volume */
    }
    fclose(fp);

    /* Build uniform-grid spatial index for nearest-neighbor lookup.
     * Bin size ≈ Voronoi cell scale; each bin then holds ~1-3 cells. */
    double box_vol = pow(2.0 * g->L, 3.0);
    double cell_vol = box_vol / g->N_cells;
    g->bin_size = pow(cell_vol, 1.0/3.0) * 1.5;
    g->Ng = (int)ceil((2.0 * g->L) / g->bin_size);
    g->bin_size = (2.0 * g->L) / g->Ng;
    long Ng3 = (long)g->Ng * g->Ng * g->Ng;

    printf("[foam] N_cells=%u L=%.2f Ng=%d bin_size=%.3f\n",
           g->N_cells, g->L, g->Ng, g->bin_size);

    g->bin_count  = calloc(Ng3, sizeof(int));
    g->bin_offset = calloc(Ng3, sizeof(int));

    /* First pass: count cells per bin */
    for (uint32_t i = 0; i < g->N_cells; i++) {
        int bx = (int)((g->cx[i] + g->L) / g->bin_size);
        int by = (int)((g->cy[i] + g->L) / g->bin_size);
        int bz = (int)((g->cz[i] + g->L) / g->bin_size);
        if (bx < 0) bx = 0; if (bx >= g->Ng) bx = g->Ng - 1;
        if (by < 0) by = 0; if (by >= g->Ng) by = g->Ng - 1;
        if (bz < 0) bz = 0; if (bz >= g->Ng) bz = g->Ng - 1;
        long b = (long)bx * g->Ng * g->Ng + by * g->Ng + bz;
        g->bin_count[b]++;
    }

    /* Prefix sum for offsets */
    int total = 0;
    for (long b = 0; b < Ng3; b++) {
        g->bin_offset[b] = total;
        total += g->bin_count[b];
    }
    g->bin_cells = malloc(sizeof(int) * total);

    /* Second pass: scatter cell IDs into bins */
    int *cursor = malloc(sizeof(int) * Ng3);
    memcpy(cursor, g->bin_offset, sizeof(int) * Ng3);
    for (uint32_t i = 0; i < g->N_cells; i++) {
        int bx = (int)((g->cx[i] + g->L) / g->bin_size);
        int by = (int)((g->cy[i] + g->L) / g->bin_size);
        int bz = (int)((g->cz[i] + g->L) / g->bin_size);
        if (bx < 0) bx = 0; if (bx >= g->Ng) bx = g->Ng - 1;
        if (by < 0) by = 0; if (by >= g->Ng) by = g->Ng - 1;
        if (bz < 0) bz = 0; if (bz >= g->Ng) bz = g->Ng - 1;
        long b = (long)bx * g->Ng * g->Ng + by * g->Ng + bz;
        g->bin_cells[cursor[b]++] = i;
    }
    free(cursor);

    return g;
}

static int find_nearest(const FoamGrid *g, double x, double y, double z) {
    /* Find nearest cell by searching the bin containing (x,y,z) and its
     * 26 neighbors (3x3x3 grid lookup). */
    int bx = (int)((x + g->L) / g->bin_size);
    int by = (int)((y + g->L) / g->bin_size);
    int bz = (int)((z + g->L) / g->bin_size);
    /* Periodic wrap */
    bx = ((bx % g->Ng) + g->Ng) % g->Ng;
    by = ((by % g->Ng) + g->Ng) % g->Ng;
    bz = ((bz % g->Ng) + g->Ng) % g->Ng;

    int best = -1;
    double best_d2 = 1e30;
    double S = 2.0 * g->L;

    for (int di = -1; di <= 1; di++) {
        int gi = ((bx + di) % g->Ng + g->Ng) % g->Ng;
        for (int dj = -1; dj <= 1; dj++) {
            int gj = ((by + dj) % g->Ng + g->Ng) % g->Ng;
            for (int dk = -1; dk <= 1; dk++) {
                int gk = ((bz + dk) % g->Ng + g->Ng) % g->Ng;
                long b = (long)gi * g->Ng * g->Ng + gj * g->Ng + gk;
                int n = g->bin_count[b];
                int o = g->bin_offset[b];
                for (int q = 0; q < n; q++) {
                    int c = g->bin_cells[o + q];
                    double dx = g->cx[c] - x;
                    double dy = g->cy[c] - y;
                    double dz = g->cz[c] - z;
                    /* Periodic distance */
                    if (dx >  0.5*S) dx -= S; if (dx < -0.5*S) dx += S;
                    if (dy >  0.5*S) dy -= S; if (dy < -0.5*S) dy += S;
                    if (dz >  0.5*S) dz -= S; if (dz < -0.5*S) dz += S;
                    double d2 = dx*dx + dy*dy + dz*dz;
                    if (d2 < best_d2) { best_d2 = d2; best = c; }
                }
            }
        }
    }
    return best;
}

/* =====================================================================
 *  Snapshot reader
 * ===================================================================== */

typedef struct {
    double t;
    uint32_t N_cells;
    double *phi[3];
    double *theta[3];
} Snapshot;

static Snapshot *snap_read(const char *path) {
    FILE *fp = fopen(path, "rb");
    if (!fp) { fprintf(stderr, "cannot open snap '%s'\n", path); return NULL; }
    char magic[4]; if (fread(magic, 1, 4, fp) != 4) { fclose(fp); return NULL; }
    if (memcmp(magic, "FSNP", 4) != 0) { fclose(fp); return NULL; }
    uint32_t version; if (fread(&version, 4, 1, fp) != 1) { fclose(fp); return NULL; }
    Snapshot *s = calloc(1, sizeof(Snapshot));
    if (fread(&s->N_cells, 4, 1, fp) != 1) { fclose(fp); free(s); return NULL; }
    if (fread(&s->t, 8, 1, fp) != 1) { fclose(fp); free(s); return NULL; }
    uint64_t flags; if (fread(&flags, 8, 1, fp) != 1) { fclose(fp); free(s); return NULL; }
    for (int a = 0; a < 3; a++) {
        s->phi[a]   = malloc(sizeof(double) * s->N_cells);
        s->theta[a] = malloc(sizeof(double) * s->N_cells);
    }
    for (uint32_t i = 0; i < s->N_cells; i++) {
        double rec[6];
        if (fread(rec, 8, 6, fp) != 6) { fclose(fp); /* leak ok */ return NULL; }
        for (int a = 0; a < 3; a++) {
            s->phi[a][i]   = rec[a];
            s->theta[a][i] = rec[3+a];
        }
    }
    fclose(fp);
    return s;
}

static void snap_free(Snapshot *s) {
    for (int a = 0; a < 3; a++) { free(s->phi[a]); free(s->theta[a]); }
    free(s);
}

/* =====================================================================
 *  Resample per-cell to NxNxN voxel grid via nearest-neighbor
 * ===================================================================== */

static void resample(const FoamGrid *g, const Snapshot *s, int N,
                     double *vox_phi[3], double *vox_theta[3]) {
    double L = g->L;
    double dx = 2.0 * L / (N - 1);
    long N3 = (long)N * N * N;
    int *nearest = malloc(sizeof(int) * N3);

    /* Find nearest cell per voxel */
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < N; i++) {
        double x = -L + i * dx;
        for (int j = 0; j < N; j++) {
            double y = -L + j * dx;
            for (int k = 0; k < N; k++) {
                double z = -L + k * dx;
                long idx = (long)i * N * N + j * N + k;
                nearest[idx] = find_nearest(g, x, y, z);
            }
        }
    }

    /* Assign field values from nearest cells */
    #pragma omp parallel for schedule(static)
    for (long idx = 0; idx < N3; idx++) {
        int c = nearest[idx];
        if (c < 0) {
            for (int a = 0; a < 3; a++) {
                vox_phi[a][idx] = 0;
                vox_theta[a][idx] = 0;
            }
            continue;
        }
        for (int a = 0; a < 3; a++) {
            vox_phi[a][idx]   = s->phi[a][c];
            vox_theta[a][idx] = s->theta[a][c];
        }
    }

    free(nearest);
}

/* =====================================================================
 *  Main
 * ===================================================================== */

int main(int argc, char **argv) {
    if (argc < 5) {
        fprintf(stderr, "Usage: %s foam_mesh.bin output.sfa N snap1.fsnp [snap2.fsnp ...]\n",
                argv[0]);
        return 1;
    }
    const char *mesh_path = argv[1];
    const char *out_path  = argv[2];
    int N = atoi(argv[3]);
    int n_snaps = argc - 4;

    int nthreads = 8;
    char *env = getenv("OMP_NUM_THREADS");
    if (env) nthreads = atoi(env);
    omp_set_num_threads(nthreads);

    FoamGrid *g = foam_load(mesh_path);
    long N3 = (long)N * N * N;

    double *vox_phi[3], *vox_theta[3];
    for (int a = 0; a < 3; a++) {
        vox_phi[a]   = malloc(sizeof(double) * N3);
        vox_theta[a] = malloc(sizeof(double) * N3);
    }

    /* Create SFA writer (F32 columns — passes float buffers, no f16 conversion) */
    SFA *sfa = sfa_create(out_path, N, N, N, g->L, g->L, g->L, 1.0);
    if (!sfa) { fprintf(stderr, "FATAL: cannot create '%s'\n", out_path); return 1; }
    sfa_add_column(sfa, "phi_x",   SFA_F32, SFA_POSITION, 0);
    sfa_add_column(sfa, "phi_y",   SFA_F32, SFA_POSITION, 1);
    sfa_add_column(sfa, "phi_z",   SFA_F32, SFA_POSITION, 2);
    sfa_add_column(sfa, "theta_x", SFA_F32, SFA_ANGLE, 0);
    sfa_add_column(sfa, "theta_y", SFA_F32, SFA_ANGLE, 1);
    sfa_add_column(sfa, "theta_z", SFA_F32, SFA_ANGLE, 2);
    if (sfa_finalize_header(sfa) != 0) {
        fprintf(stderr, "FATAL: finalize_header failed\n"); return 1;
    }

    /* Float buffers for writing (SFA_F32) */
    float *fbuf[6];
    for (int a = 0; a < 6; a++) fbuf[a] = malloc(sizeof(float) * N3);

    for (int s = 0; s < n_snaps; s++) {
        const char *snap_path = argv[4 + s];
        printf("[%d/%d] %s\n", s+1, n_snaps, snap_path);
        Snapshot *snp = snap_read(snap_path);
        if (!snp) { fprintf(stderr, "skip\n"); continue; }
        if (snp->N_cells != g->N_cells) {
            fprintf(stderr, "FATAL: snapshot N_cells=%u != mesh N_cells=%u\n",
                    snp->N_cells, g->N_cells);
            return 1;
        }
        resample(g, snp, N, vox_phi, vox_theta);

        /* Convert double → float */
        for (int a = 0; a < 3; a++) {
            for (long idx = 0; idx < N3; idx++) {
                fbuf[a][idx]   = (float)vox_phi[a][idx];
                fbuf[3+a][idx] = (float)vox_theta[a][idx];
            }
        }
        void *cols[6] = { fbuf[0], fbuf[1], fbuf[2], fbuf[3], fbuf[4], fbuf[5] };
        sfa_write_frame(sfa, snp->t, cols);

        printf("  t=%.2f -> frame written (N=%d^3)\n", snp->t, N);
        snap_free(snp);
    }

    for (int a = 0; a < 6; a++) free(fbuf[a]);
    sfa_close(sfa);

    for (int a = 0; a < 3; a++) { free(vox_phi[a]); free(vox_theta[a]); }
    return 0;
}
