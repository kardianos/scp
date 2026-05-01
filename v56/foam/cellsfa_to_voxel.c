/*  cellsfa_to_voxel.c — Convert a cell-native SFA (FMSH + FCEL frames)
 *  into a voxel SFA so existing tools (volview, sfa info, etc.) can read it.
 *
 *  Build:
 *      gcc -O3 -fopenmp -o cellsfa_to_voxel cellsfa_to_voxel.c \
 *          -I../../sfa/format -lzstd -lm
 *
 *  Usage:
 *      cellsfa_to_voxel <input_cell.sfa> <output_voxel.sfa> <voxel_N>
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

/* Spatial bin index over cell positions, identical pattern to foam_to_voxel. */
typedef struct {
    int Ng;
    double L;
    double bin_size;
    int *count;
    int *offset;
    int *cells;
    int  voxel_N;
    int *voxel_to_cell;
    double *cx, *cy, *cz;   /* cell centroids (for nearest-neighbour search) */
} ResIdx;

static ResIdx *res_build(const SFA_Mesh *m, int voxel_N) {
    ResIdx *r = calloc(1, sizeof(ResIdx));
    r->L = m->L;
    double box_vol = pow(2.0 * r->L, 3.0);
    double cell_vol = box_vol / m->N_cells;
    r->bin_size = pow(cell_vol, 1.0/3.0) * 1.5;
    r->Ng = (int)ceil((2.0 * r->L) / r->bin_size);
    r->bin_size = (2.0 * r->L) / r->Ng;
    long Ng3 = (long)r->Ng * r->Ng * r->Ng;
    r->count  = calloc(Ng3, sizeof(int));
    r->offset = calloc(Ng3, sizeof(int));

    /* Copy positions out for fast access; faster than going through struct */
    r->cx = malloc(sizeof(double) * m->N_cells);
    r->cy = malloc(sizeof(double) * m->N_cells);
    r->cz = malloc(sizeof(double) * m->N_cells);
    for (uint32_t i = 0; i < m->N_cells; i++) {
        r->cx[i] = m->cell_pos[3*i + 0];
        r->cy[i] = m->cell_pos[3*i + 1];
        r->cz[i] = m->cell_pos[3*i + 2];
    }

    for (uint32_t i = 0; i < m->N_cells; i++) {
        int bx = (int)((r->cx[i] + r->L) / r->bin_size);
        int by = (int)((r->cy[i] + r->L) / r->bin_size);
        int bz = (int)((r->cz[i] + r->L) / r->bin_size);
        if (bx < 0) bx = 0; if (bx >= r->Ng) bx = r->Ng - 1;
        if (by < 0) by = 0; if (by >= r->Ng) by = r->Ng - 1;
        if (bz < 0) bz = 0; if (bz >= r->Ng) bz = r->Ng - 1;
        r->count[(long)bx * r->Ng * r->Ng + by * r->Ng + bz]++;
    }
    int total = 0;
    for (long b = 0; b < Ng3; b++) { r->offset[b] = total; total += r->count[b]; }
    r->cells = malloc(sizeof(int) * total);
    int *cursor = malloc(sizeof(int) * Ng3);
    memcpy(cursor, r->offset, sizeof(int) * Ng3);
    for (uint32_t i = 0; i < m->N_cells; i++) {
        int bx = (int)((r->cx[i] + r->L) / r->bin_size);
        int by = (int)((r->cy[i] + r->L) / r->bin_size);
        int bz = (int)((r->cz[i] + r->L) / r->bin_size);
        if (bx < 0) bx = 0; if (bx >= r->Ng) bx = r->Ng - 1;
        if (by < 0) by = 0; if (by >= r->Ng) by = r->Ng - 1;
        if (bz < 0) bz = 0; if (bz >= r->Ng) bz = r->Ng - 1;
        r->cells[cursor[(long)bx * r->Ng * r->Ng + by * r->Ng + bz]++] = (int)i;
    }
    free(cursor);

    /* Pre-compute voxel→cell map */
    r->voxel_N = voxel_N;
    long N3 = (long)voxel_N * voxel_N * voxel_N;
    r->voxel_to_cell = malloc(sizeof(int) * N3);
    double L = r->L;
    double dx = 2.0 * L / (voxel_N - 1);
    double S = 2.0 * L;
    int Ng = r->Ng;

    printf("[resamp] voxel %d^3, mesh %u cells, bin_size=%.3f, Ng=%d\n",
           voxel_N, m->N_cells, r->bin_size, Ng);

    #pragma omp parallel for schedule(dynamic, 32)
    for (int i = 0; i < voxel_N; i++) {
        double x = -L + i * dx;
        int bx0 = ((int)((x + L) / r->bin_size) % Ng + Ng) % Ng;
        for (int j = 0; j < voxel_N; j++) {
            double y = -L + j * dx;
            int by0 = ((int)((y + L) / r->bin_size) % Ng + Ng) % Ng;
            for (int k = 0; k < voxel_N; k++) {
                double z = -L + k * dx;
                int bz0 = ((int)((z + L) / r->bin_size) % Ng + Ng) % Ng;
                int best = -1;
                double best_d2 = 1e30;
                for (int di = -1; di <= 1; di++) {
                    int gi = ((bx0 + di) % Ng + Ng) % Ng;
                    for (int dj = -1; dj <= 1; dj++) {
                        int gj = ((by0 + dj) % Ng + Ng) % Ng;
                        for (int dk = -1; dk <= 1; dk++) {
                            int gk = ((bz0 + dk) % Ng + Ng) % Ng;
                            long b = (long)gi * Ng * Ng + gj * Ng + gk;
                            int n = r->count[b], o = r->offset[b];
                            for (int q = 0; q < n; q++) {
                                int c = r->cells[o + q];
                                double cdx = r->cx[c] - x;
                                double cdy = r->cy[c] - y;
                                double cdz = r->cz[c] - z;
                                if (cdx >  0.5*S) cdx -= S; if (cdx < -0.5*S) cdx += S;
                                if (cdy >  0.5*S) cdy -= S; if (cdy < -0.5*S) cdy += S;
                                if (cdz >  0.5*S) cdz -= S; if (cdz < -0.5*S) cdz += S;
                                double d2 = cdx*cdx + cdy*cdy + cdz*cdz;
                                if (d2 < best_d2) { best_d2 = d2; best = c; }
                            }
                        }
                    }
                }
                r->voxel_to_cell[(long)i * voxel_N * voxel_N + j * voxel_N + k] = best;
            }
        }
    }
    return r;
}

static void res_free(ResIdx *r) {
    free(r->count); free(r->offset); free(r->cells); free(r->voxel_to_cell);
    free(r->cx); free(r->cy); free(r->cz);
    free(r);
}

int main(int argc, char **argv) {
    if (argc < 4) {
        fprintf(stderr, "Usage: %s <cell.sfa> <voxel.sfa> <voxel_N>\n", argv[0]);
        return 1;
    }
    const char *in_path  = argv[1];
    const char *out_path = argv[2];
    int N = atoi(argv[3]);
    int nthreads = 8;
    char *env = getenv("OMP_NUM_THREADS");
    if (env) nthreads = atoi(env);
    omp_set_num_threads(nthreads);

    SFA *in = sfa_open(in_path);
    if (!in) { fprintf(stderr, "cannot open %s\n", in_path); return 1; }
    printf("[in] %u frames, %u cols\n", in->total_frames, in->n_columns);

    SFA_Mesh mesh = {0};
    int have_mesh = 0;
    ResIdx *resamp = NULL;
    long N3 = (long)N * N * N;

    SFA *out = NULL;
    /* Output column count is decided when the first FCEL is encountered
     * (any value up to MAX_COLS — handles v55's 6-column and v56's
     * 8-column outputs without a recompile). */
    #define MAX_COLS 16
    float *vox_buf[MAX_COLS] = {0};
    int  out_ncols = 0;

    /* Per-cell state buffer (column-major: NCOL × N_cells f32 floats).
     * Allocated when first FCEL is seen. */
    float *cell_state = NULL;
    /* Temporal model (loaded from FCEL v2 frames). */
    float *t_mean = NULL, *t_amp = NULL, *t_phase = NULL;
    float t_omega = 0.0f;
    int   model_valid = 0;

    for (uint32_t f = 0; f < in->total_frames; f++) {
        SFA_L2Entry e;
        if (sfa_find_frame(in, f, &e) < 0) {
            fprintf(stderr, "find_frame %u failed\n", f); return 1;
        }
        if (e.frame_type == SFA_FRAME_MESH) {
            if (have_mesh) { sfa_mesh_free(&mesh); }
            if (resamp) { res_free(resamp); resamp = NULL; }
            if (sfa_read_mesh_frame(in, f, &mesh) != 0) {
                fprintf(stderr, "read_mesh_frame %u failed\n", f); return 1;
            }
            have_mesh = 1;
            resamp = res_build(&mesh, N);
            /* cell_state allocation is deferred until the first FCEL so
             * we know n_columns. */
            free(cell_state); cell_state = NULL;
            free(t_mean); free(t_amp); free(t_phase);
            t_mean = t_amp = t_phase = NULL; model_valid = 0;
            printf("[%u] FMSH t=%.3f  N_cells=%u  L=%.2f\n",
                   f, e.time, mesh.N_cells, mesh.L);
        } else if (e.frame_type == SFA_FRAME_CELL) {
            if (!have_mesh) {
                fprintf(stderr, "FCEL frame %u with no preceding mesh\n", f);
                return 1;
            }
            /* Always read the cell values (column-major into a temp buf). */
            uint32_t Nc, ncols;
            uint8_t  dtype;
            void *data = NULL;
            if (sfa_read_cell_frame(in, f, &Nc, &ncols, &dtype, &data) != 0) {
                fprintf(stderr, "read_cell_frame %u failed\n", f); return 1;
            }
            if (Nc != mesh.N_cells || dtype != SFA_F32 || ncols > MAX_COLS) {
                fprintf(stderr, "cell frame mismatch (Nc=%u, cols=%u, dtype=%d)\n",
                        Nc, ncols, dtype);
                free(data); return 1;
            }
            /* On first FCEL: allocate state and create output SFA */
            if (!out) {
                out_ncols = (int)ncols;
                cell_state = (float*)malloc(sizeof(float) * out_ncols * Nc);
                for (int k = 0; k < out_ncols; k++) {
                    vox_buf[k] = malloc(sizeof(float) * N3);
                }
                out = sfa_create(out_path, N, N, N, mesh.L, mesh.L, mesh.L, in->dt);
                for (int k = 0; k < out_ncols; k++) {
                    char name[16]; snprintf(name, 16, "C%d", k);
                    sfa_add_column(out, name, SFA_F32,
                                    (k < 4) ? SFA_CUSTOM : SFA_POSITION, k & 3);
                }
                sfa_finalize_header(out);
                printf("  output: %d columns × %d^3 voxels\n", out_ncols, N);
            }
            if ((int)ncols != out_ncols) {
                fprintf(stderr, "cell frame ncols=%u changed (out_ncols=%d)\n",
                        ncols, out_ncols);
                free(data); return 1;
            }
            memcpy(cell_state, data, (size_t)out_ncols * Nc * 4);
            free(data);
            /* If this FCEL embeds a temporal model (v2), load it. */
            uint32_t mNc, mNcols;
            float *nm=NULL, *nA=NULL, *nP=NULL;
            float nOmega = 0;
            if (sfa_read_cell_iframe_temporal(in, f, &nOmega, &mNc, &mNcols,
                                               &nm, &nA, &nP) == 0
                && mNc == Nc && (int)mNcols == out_ncols) {
                free(t_mean); free(t_amp); free(t_phase);
                t_mean = nm; t_amp = nA; t_phase = nP;
                t_omega = nOmega;
                model_valid = 1;
            }
            /* Output current state to voxel frame */
            for (int col = 0; col < out_ncols; col++) {
                float *src = cell_state + (long)col * Nc;
                #pragma omp parallel for schedule(static)
                for (long idx = 0; idx < N3; idx++) {
                    int c = resamp->voxel_to_cell[idx];
                    vox_buf[col][idx] = (c >= 0) ? src[c] : 0.0f;
                }
            }
            void *cols[MAX_COLS];
            for (int k = 0; k < out_ncols; k++) cols[k] = vox_buf[k];
            sfa_write_frame(out, e.time, cols);
            if (f % 10 == 0)
                printf("[%u] FCEL%s t=%.3f → voxel\n",
                       f, model_valid ? " v2 (+model)" : "", e.time);
        } else if (e.frame_type == SFA_FRAME_TEMPORAL_MODEL) {
            uint32_t mNc, mNcols;
            float *nm=NULL, *nA=NULL, *nP=NULL;
            float nOmega = 0;
            if (sfa_read_temporal_model_frame(in, f, &nOmega, &mNc, &mNcols,
                                                &nm, &nA, &nP) != 0) {
                fprintf(stderr, "read_temporal_model_frame %u failed\n", f); return 1;
            }
            free(t_mean); free(t_amp); free(t_phase);
            t_mean = nm; t_amp = nA; t_phase = nP;
            t_omega = nOmega;
            model_valid = 1;
            printf("[%u] FMTL t=%.3f  ω=%.3f (model frame)\n", f, e.time, nOmega);
        } else if (e.frame_type == SFA_FRAME_CELL_P) {
            if (!have_mesh || !model_valid) {
                fprintf(stderr, "FCEP frame %u without mesh+model loaded\n", f);
                return 1;
            }
            uint32_t Nc, ncols, n_changed;
            uint32_t *cell_ids = NULL;
            float    *deltas = NULL;
            if (sfa_read_cell_pframe(in, f, &Nc, &ncols, &n_changed,
                                       &cell_ids, &deltas) != 0) {
                fprintf(stderr, "read_cell_pframe %u failed\n", f); return 1;
            }
            if (Nc != mesh.N_cells || (int)ncols != out_ncols) {
                fprintf(stderr, "FCEP mismatch (Nc=%u out_ncols=%d ncols=%u)\n",
                        Nc, out_ncols, ncols);
                free(cell_ids); free(deltas); return 1;
            }
            /* Reconstruct: state[c, k] = predict(t)[c, k] (+ delta[c, k] if c in change list) */
            double t = e.time;
            #pragma omp parallel for schedule(static)
            for (uint32_t c = 0; c < Nc; c++) {
                long base = (long)c * out_ncols;
                for (int k = 0; k < out_ncols; k++) {
                    float pred = t_mean[base + k]
                               + t_amp[base + k]
                               * cosf(t_omega * (float)t + t_phase[base + k]);
                    cell_state[(long)k * Nc + c] = pred;
                }
            }
            for (uint32_t i = 0; i < n_changed; i++) {
                uint32_t c = cell_ids[i];
                for (int k = 0; k < out_ncols; k++) {
                    cell_state[(long)k * Nc + c] += deltas[(long)i * out_ncols + k];
                }
            }
            free(cell_ids); free(deltas);
            /* Output current state */
            for (int col = 0; col < out_ncols; col++) {
                float *src = cell_state + (long)col * Nc;
                #pragma omp parallel for schedule(static)
                for (long idx = 0; idx < N3; idx++) {
                    int c = resamp->voxel_to_cell[idx];
                    vox_buf[col][idx] = (c >= 0) ? src[c] : 0.0f;
                }
            }
            void *cols[MAX_COLS];
            for (int k = 0; k < out_ncols; k++) cols[k] = vox_buf[k];
            sfa_write_frame(out, e.time, cols);
            if (f % 10 == 0)
                printf("[%u] FCEP t=%.3f n_changed=%u → voxel\n",
                       f, e.time, n_changed);
        } else {
            printf("[%u] skip frame_type=%u\n", f, e.frame_type);
        }
    }
    free(cell_state);
    free(t_mean); free(t_amp); free(t_phase);

    if (out) sfa_close(out);
    if (resamp) res_free(resamp);
    if (have_mesh) sfa_mesh_free(&mesh);
    sfa_close(in);
    for (int i = 0; i < MAX_COLS; i++) if (vox_buf[i]) free(vox_buf[i]);
    return 0;
}
