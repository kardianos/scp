/*  sfa_vectorize3d.c --- 3D tensor-product polynomial vectorization of SFA field data
 *
 *  Pipeline:
 *    1. Read SFA voxel file
 *    2. Divide grid into 3D blocks (default 8^3)
 *    3. Fit tricubic polynomial (4x4x4 = 64 coefficients) per block
 *    4. Adaptive refinement: split high-error blocks, merge low-error neighbors
 *    5. Output: vector SFA file + reconstructed voxel SFA for comparison
 *
 *  Build: gcc -O3 -march=native -fopenmp -o sfa_vectorize3d sfa_vectorize3d.c -lzstd -lm
 *  Usage: ./sfa_vectorize3d input.sfa [-frame N] [-tol 0.01] [-block 8] [-max_iter 5]
 *         [-o_recon recon.sfa] [-o_vec vec.sfa] [-fields phi]
 */

#define SFA_IMPLEMENTATION
#include "../format/sfa.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>

/* ===================================================================
   Constants and configuration
   =================================================================== */

#define MAX_POLY_ORDER  4   /* 0..3 => constant, linear, quadratic, cubic */
#define MAX_COEFFS      64  /* 4^3 for tricubic */
#define MAX_BLOCKS       (256*256*256 / 8)  /* upper bound for 256^3 grid / min block */
#define MAX_REFINE_ITER  10
#define MIN_BLOCK_SIZE   4   /* don't split below 4^3 */

/* ===================================================================
   f16 conversion (from existing code)
   =================================================================== */

static float f16f(uint16_t h) {
    uint16_t s = h & 0x8000; int e = (h >> 10) & 0x1F; uint16_t m = h & 0x3FF;
    if (e == 0) return 0; if (e == 31) return s ? -1e30f : 1e30f;
    float f; uint32_t x = ((uint32_t)s << 16) | ((uint32_t)(e-15+127) << 23) | ((uint32_t)m << 13);
    memcpy(&f, &x, 4); return f;
}

static uint16_t f16h(float f) {
    uint32_t x; memcpy(&x, &f, 4);
    uint16_t s = (x >> 16) & 0x8000;
    int e = ((x >> 23) & 0xFF) - 127 + 15;
    uint16_t m = (x >> 13) & 0x3FF;
    if (e <= 0) return s;
    if (e >= 31) return s | 0x7C00;
    return s | (e << 10) | m;
}

/* Extract column c, voxel i as float */
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

/* ===================================================================
   3D Block Patch
   =================================================================== */

typedef struct {
    int bx, by, bz;         /* block origin in voxel coordinates */
    int sx, sy, sz;          /* block size in each dimension */
    int order;               /* polynomial order: 0=const, 1=linear, 2=quad, 3=cubic */
    int n_coeffs;            /* (order+1)^3 */
    double coeffs[MAX_COEFFS]; /* tensor-product polynomial coefficients */
    double max_error;        /* L-inf error in this block */
    double rms_error;        /* RMS error in this block */
    int active;              /* 1=active, 0=merged into another */
} Patch3D;

/* Dynamic array of patches */
typedef struct {
    Patch3D *patches;
    int count;
    int capacity;
} PatchList;

static void patchlist_init(PatchList *pl, int cap) {
    pl->count = 0;
    pl->capacity = cap;
    pl->patches = (Patch3D*)calloc(cap, sizeof(Patch3D));
}

static int patchlist_push(PatchList *pl, Patch3D *p) {
    if (pl->count >= pl->capacity) {
        pl->capacity *= 2;
        pl->patches = (Patch3D*)realloc(pl->patches, pl->capacity * sizeof(Patch3D));
    }
    pl->patches[pl->count] = *p;
    return pl->count++;
}

static void patchlist_free(PatchList *pl) {
    free(pl->patches);
    pl->patches = NULL;
    pl->count = pl->capacity = 0;
}

/* ===================================================================
   Least-squares tensor-product polynomial fit (3D)

   We fit f(u,v,w) = sum_{i,j,k=0}^{order} c[i*(O+1)^2 + j*(O+1) + k] * u^i * v^j * w^k

   where u,v,w in [0,1] parameterize the block.
   =================================================================== */

/* Solve Ax = b via Gaussian elimination with partial pivoting.
 * A is n x n, b is n x 1. Solution written into b. */
static int solve_linear(double *A, double *b, int n) {
    for (int col = 0; col < n; col++) {
        /* Partial pivot */
        int piv = col;
        for (int row = col + 1; row < n; row++)
            if (fabs(A[row * n + col]) > fabs(A[piv * n + col])) piv = row;
        if (piv != col) {
            for (int j = 0; j < n; j++) {
                double tmp = A[col * n + j]; A[col * n + j] = A[piv * n + j]; A[piv * n + j] = tmp;
            }
            double tmp = b[col]; b[col] = b[piv]; b[piv] = tmp;
        }
        double d = A[col * n + col];
        if (fabs(d) < 1e-15) return -1; /* singular */
        for (int row = col + 1; row < n; row++) {
            double f = A[row * n + col] / d;
            for (int j = col; j < n; j++) A[row * n + j] -= f * A[col * n + j];
            b[row] -= f * b[col];
        }
    }
    /* Back substitution */
    for (int row = n - 1; row >= 0; row--) {
        double sum = b[row];
        for (int j = row + 1; j < n; j++) sum -= A[row * n + j] * b[j];
        b[row] = sum / A[row * n + row];
    }
    return 0;
}

/* Precompute 1D Vandermonde basis powers: basis[order+1][n_pts]
 * basis[d][i] = t_i^d  where t_i = i/(n-1) */
static void vandermonde_1d(int n_pts, int order, double *basis) {
    int O1 = order + 1;
    for (int i = 0; i < n_pts; i++) {
        double t = (n_pts > 1) ? (double)i / (n_pts - 1) : 0.0;
        double tk = 1.0;
        for (int d = 0; d < O1; d++) {
            basis[d * n_pts + i] = tk;
            tk *= t;
        }
    }
}

/* Fit a 3D tensor-product polynomial to a block of data.
 * data: full 3D field array (Nx * Ny * Nz), row-major [x][y][z]
 * patch: defines block origin and size; coeffs are written here. */
static void fit_patch(const float *data, int Nx, int Ny, int Nz,
                      Patch3D *patch, int order) {
    int sx = patch->sx, sy = patch->sy, sz = patch->sz;
    int bx = patch->bx, by = patch->by, bz = patch->bz;
    int O1 = order + 1;
    int nc = O1 * O1 * O1;  /* number of coefficients */

    patch->order = order;
    patch->n_coeffs = nc;
    memset(patch->coeffs, 0, sizeof(patch->coeffs));

    /* Build 1D bases */
    double *bx_basis = (double*)malloc(O1 * sx * sizeof(double));
    double *by_basis = (double*)malloc(O1 * sy * sizeof(double));
    double *bz_basis = (double*)malloc(O1 * sz * sizeof(double));
    vandermonde_1d(sx, order, bx_basis);
    vandermonde_1d(sy, order, by_basis);
    vandermonde_1d(sz, order, bz_basis);

    /* Build normal equations: A * c = rhs
     * A[a][b] = sum_{ijk} phi_a(ijk) * phi_b(ijk)
     * rhs[a]  = sum_{ijk} phi_a(ijk) * data(ijk)
     * where phi_a = u^{ax} * v^{ay} * w^{az} is the a-th basis function */
    double *A = (double*)calloc(nc * nc, sizeof(double));
    double *rhs = (double*)calloc(nc, sizeof(double));

    for (int ix = 0; ix < sx; ix++) {
        int gx = bx + ix;
        if (gx >= Nx) continue;
        for (int iy = 0; iy < sy; iy++) {
            int gy = by + iy;
            if (gy >= Ny) continue;
            for (int iz = 0; iz < sz; iz++) {
                int gz = bz + iz;
                if (gz >= Nz) continue;

                float val = data[(long)gx * Ny * Nz + gy * Nz + gz];

                /* Compute tensor-product basis values for this voxel */
                double phi[MAX_COEFFS];
                int a = 0;
                for (int dx = 0; dx < O1; dx++)
                    for (int dy = 0; dy < O1; dy++)
                        for (int dz = 0; dz < O1; dz++)
                            phi[a++] = bx_basis[dx * sx + ix]
                                     * by_basis[dy * sy + iy]
                                     * bz_basis[dz * sz + iz];

                /* Accumulate normal equations */
                for (int ai = 0; ai < nc; ai++) {
                    rhs[ai] += phi[ai] * val;
                    for (int bi = ai; bi < nc; bi++)
                        A[ai * nc + bi] += phi[ai] * phi[bi];
                }
            }
        }
    }

    /* Symmetrize A */
    for (int ai = 0; ai < nc; ai++)
        for (int bi = 0; bi < ai; bi++)
            A[ai * nc + bi] = A[bi * nc + ai];

    /* Solve */
    if (solve_linear(A, rhs, nc) == 0) {
        memcpy(patch->coeffs, rhs, nc * sizeof(double));
    }

    free(A); free(rhs);
    free(bx_basis); free(by_basis); free(bz_basis);
}

/* Evaluate the patch polynomial at local coordinates (u,v,w) in [0,1]^3 */
static double eval_patch(const Patch3D *p, double u, double v, double w) {
    int O1 = p->order + 1;
    double val = 0;
    int a = 0;
    double upow[MAX_POLY_ORDER], vpow[MAX_POLY_ORDER], wpow[MAX_POLY_ORDER];
    upow[0] = 1; vpow[0] = 1; wpow[0] = 1;
    for (int d = 1; d < O1; d++) {
        upow[d] = upow[d-1] * u;
        vpow[d] = vpow[d-1] * v;
        wpow[d] = wpow[d-1] * w;
    }
    for (int dx = 0; dx < O1; dx++)
        for (int dy = 0; dy < O1; dy++)
            for (int dz = 0; dz < O1; dz++)
                val += p->coeffs[a++] * upow[dx] * vpow[dy] * wpow[dz];
    return val;
}

/* Compute error statistics for a patch */
static void score_patch(const float *data, int Nx, int Ny, int Nz, Patch3D *patch) {
    int sx = patch->sx, sy = patch->sy, sz = patch->sz;
    int bx = patch->bx, by = patch->by, bz = patch->bz;
    double max_err = 0, sum_sq = 0;
    int count = 0;

    for (int ix = 0; ix < sx; ix++) {
        int gx = bx + ix;
        if (gx >= Nx) continue;
        double u = (sx > 1) ? (double)ix / (sx - 1) : 0.0;
        for (int iy = 0; iy < sy; iy++) {
            int gy = by + iy;
            if (gy >= Ny) continue;
            double v = (sy > 1) ? (double)iy / (sy - 1) : 0.0;
            for (int iz = 0; iz < sz; iz++) {
                int gz = bz + iz;
                if (gz >= Nz) continue;
                double w = (sz > 1) ? (double)iz / (sz - 1) : 0.0;

                float orig = data[(long)gx * Ny * Nz + gy * Nz + gz];
                double approx = eval_patch(patch, u, v, w);
                double err = fabs(orig - approx);
                if (err > max_err) max_err = err;
                sum_sq += err * err;
                count++;
            }
        }
    }
    patch->max_error = max_err;
    patch->rms_error = (count > 0) ? sqrt(sum_sq / count) : 0;
}

/* ===================================================================
   Initial blocking: divide grid into uniform blocks
   =================================================================== */

static void create_initial_blocks(PatchList *pl, int Nx, int Ny, int Nz,
                                  int block_size, const float *data, int order) {
    int nbx = (Nx + block_size - 1) / block_size;
    int nby = (Ny + block_size - 1) / block_size;
    int nbz = (Nz + block_size - 1) / block_size;

    patchlist_init(pl, nbx * nby * nbz);

    for (int ix = 0; ix < nbx; ix++) {
        for (int iy = 0; iy < nby; iy++) {
            for (int iz = 0; iz < nbz; iz++) {
                Patch3D p = {0};
                p.bx = ix * block_size;
                p.by = iy * block_size;
                p.bz = iz * block_size;
                p.sx = (p.bx + block_size <= Nx) ? block_size : (Nx - p.bx);
                p.sy = (p.by + block_size <= Ny) ? block_size : (Ny - p.by);
                p.sz = (p.bz + block_size <= Nz) ? block_size : (Nz - p.bz);
                p.active = 1;
                patchlist_push(pl, &p);
            }
        }
    }

    /* Fit all patches in parallel */
    #pragma omp parallel for schedule(dynamic, 4)
    for (int i = 0; i < pl->count; i++) {
        fit_patch(data, Nx, Ny, Nz, &pl->patches[i], order);
        score_patch(data, Nx, Ny, Nz, &pl->patches[i]);
    }
}

/* ===================================================================
   Re-rasterization: reconstruct voxel data from patches
   =================================================================== */

static void rasterize_patches(const PatchList *pl, float *recon, int Nx, int Ny, int Nz) {
    long N3 = (long)Nx * Ny * Nz;
    memset(recon, 0, N3 * sizeof(float));

    #pragma omp parallel for schedule(dynamic, 4)
    for (int pi = 0; pi < pl->count; pi++) {
        const Patch3D *p = &pl->patches[pi];
        if (!p->active) continue;
        int sx = p->sx, sy = p->sy, sz = p->sz;
        int bx = p->bx, by = p->by, bz = p->bz;

        for (int ix = 0; ix < sx; ix++) {
            int gx = bx + ix;
            if (gx >= Nx) continue;
            double u = (sx > 1) ? (double)ix / (sx - 1) : 0.0;
            for (int iy = 0; iy < sy; iy++) {
                int gy = by + iy;
                if (gy >= Ny) continue;
                double v = (sy > 1) ? (double)iy / (sy - 1) : 0.0;
                for (int iz = 0; iz < sz; iz++) {
                    int gz = bz + iz;
                    if (gz >= Nz) continue;
                    double w = (sz > 1) ? (double)iz / (sz - 1) : 0.0;
                    recon[(long)gx * Ny * Nz + gy * Nz + gz] = (float)eval_patch(p, u, v, w);
                }
            }
        }
    }
}

/* ===================================================================
   Error scoring
   =================================================================== */

typedef struct {
    double linf;      /* max absolute error */
    double rms;       /* root-mean-square error */
    double mean_abs;  /* mean absolute error */
    double max_val;   /* max |data| for relative error */
    double rel_linf;  /* linf / max_val */
    int n_patches;
    int n_coeffs_total;
    long n_voxels;
    double compression; /* n_voxels / n_coeffs_total (as floats) */
} ErrorScore;

static ErrorScore compute_error(const float *orig, const float *recon,
                                int Nx, int Ny, int Nz,
                                const PatchList *pl) {
    ErrorScore s = {0};
    long N3 = (long)Nx * Ny * Nz;
    s.n_voxels = N3;

    double max_err = 0, sum_sq = 0, sum_abs = 0, max_val = 0;

    #pragma omp parallel for reduction(max:max_err,max_val) reduction(+:sum_sq,sum_abs)
    for (long i = 0; i < N3; i++) {
        double err = fabs(orig[i] - recon[i]);
        double av = fabs(orig[i]);
        if (err > max_err) max_err = err;
        if (av > max_val) max_val = av;
        sum_sq += err * err;
        sum_abs += err;
    }

    s.linf = max_err;
    s.rms = sqrt(sum_sq / N3);
    s.mean_abs = sum_abs / N3;
    s.max_val = max_val;
    s.rel_linf = (max_val > 1e-15) ? max_err / max_val : 0;

    /* Count active patches and total coefficients */
    int np = 0, nc = 0;
    for (int i = 0; i < pl->count; i++) {
        if (pl->patches[i].active) {
            np++;
            nc += pl->patches[i].n_coeffs;
        }
    }
    s.n_patches = np;
    s.n_coeffs_total = nc;
    s.compression = (nc > 0) ? (double)N3 / nc : 0;

    return s;
}

/* Per-region error: split into "core" (center 50%) and "edge" (outer shell) */
typedef struct {
    double core_rms, core_linf;
    double edge_rms, edge_linf;
} RegionError;

static RegionError compute_region_error(const float *orig, const float *recon,
                                        int Nx, int Ny, int Nz) {
    RegionError r = {0};
    int cx0 = Nx / 4, cx1 = 3 * Nx / 4;
    int cy0 = Ny / 4, cy1 = 3 * Ny / 4;
    int cz0 = Nz / 4, cz1 = 3 * Nz / 4;

    double core_sq = 0, edge_sq = 0;
    double core_max = 0, edge_max = 0;
    long core_n = 0, edge_n = 0;

    for (int ix = 0; ix < Nx; ix++)
    for (int iy = 0; iy < Ny; iy++)
    for (int iz = 0; iz < Nz; iz++) {
        long idx = (long)ix * Ny * Nz + iy * Nz + iz;
        double err = fabs(orig[idx] - recon[idx]);
        int in_core = (ix >= cx0 && ix < cx1 && iy >= cy0 && iy < cy1 && iz >= cz0 && iz < cz1);
        if (in_core) {
            core_sq += err * err;
            if (err > core_max) core_max = err;
            core_n++;
        } else {
            edge_sq += err * err;
            if (err > edge_max) edge_max = err;
            edge_n++;
        }
    }

    r.core_rms = (core_n > 0) ? sqrt(core_sq / core_n) : 0;
    r.core_linf = core_max;
    r.edge_rms = (edge_n > 0) ? sqrt(edge_sq / edge_n) : 0;
    r.edge_linf = edge_max;
    return r;
}

/* ===================================================================
   Adaptive refinement: split high-error patches, merge low-error neighbors
   =================================================================== */

/* Split a patch into up to 8 sub-patches (2x2x2).
 * Collects children into a temporary array to avoid realloc issues. */
static int split_patch(PatchList *pl, int patch_idx, const float *data,
                       int Nx, int Ny, int Nz, int order) {
    /* Copy parent data before any realloc */
    Patch3D parent = pl->patches[patch_idx];
    pl->patches[patch_idx].active = 0;

    int hx = parent.sx / 2, hy = parent.sy / 2, hz = parent.sz / 2;
    if (hx < MIN_BLOCK_SIZE) hx = parent.sx;
    if (hy < MIN_BLOCK_SIZE) hy = parent.sy;
    if (hz < MIN_BLOCK_SIZE) hz = parent.sz;

    int nx = (hx < parent.sx) ? 2 : 1;
    int ny = (hy < parent.sy) ? 2 : 1;
    int nz = (hz < parent.sz) ? 2 : 1;

    /* If we can't actually split, reactivate and return */
    if (nx * ny * nz <= 1) {
        pl->patches[patch_idx].active = 1;
        return 0;
    }

    Patch3D children[8];
    int nc = 0;

    for (int dx = 0; dx < nx; dx++)
    for (int dy = 0; dy < ny; dy++)
    for (int dz = 0; dz < nz; dz++) {
        Patch3D *child = &children[nc++];
        memset(child, 0, sizeof(*child));
        child->bx = parent.bx + dx * hx;
        child->by = parent.by + dy * hy;
        child->bz = parent.bz + dz * hz;
        child->sx = (dx == nx - 1) ? (parent.sx - dx * hx) : hx;
        child->sy = (dy == ny - 1) ? (parent.sy - dy * hy) : hy;
        child->sz = (dz == nz - 1) ? (parent.sz - dz * hz) : hz;
        child->active = 1;

        fit_patch(data, Nx, Ny, Nz, child, order);
        score_patch(data, Nx, Ny, Nz, child);
    }

    /* Now push all children (safe even if realloc happens) */
    for (int i = 0; i < nc; i++)
        patchlist_push(pl, &children[i]);

    return nc;
}

/* Try to merge two adjacent patches (same dimension along one axis) */
static int try_merge_patches(const float *data, int Nx, int Ny, int Nz,
                             Patch3D *a, Patch3D *b, int order, double tol,
                             Patch3D *merged) {
    if (!a->active || !b->active) return 0;

    /* Check adjacency: must share a face and have matching extents on the other two dims */
    int adj = 0;
    *merged = (Patch3D){0};
    merged->active = 1;

    /* Adjacent along X? */
    if (a->bx + a->sx == b->bx && a->by == b->by && a->bz == b->bz &&
        a->sy == b->sy && a->sz == b->sz) {
        merged->bx = a->bx; merged->by = a->by; merged->bz = a->bz;
        merged->sx = a->sx + b->sx; merged->sy = a->sy; merged->sz = a->sz;
        adj = 1;
    }
    /* Adjacent along Y? */
    else if (a->by + a->sy == b->by && a->bx == b->bx && a->bz == b->bz &&
             a->sx == b->sx && a->sz == b->sz) {
        merged->bx = a->bx; merged->by = a->by; merged->bz = a->bz;
        merged->sx = a->sx; merged->sy = a->sy + b->sy; merged->sz = a->sz;
        adj = 1;
    }
    /* Adjacent along Z? */
    else if (a->bz + a->sz == b->bz && a->bx == b->bx && a->by == b->by &&
             a->sx == b->sx && a->sy == b->sy) {
        merged->bx = a->bx; merged->by = a->by; merged->bz = a->bz;
        merged->sx = a->sx; merged->sy = a->sy; merged->sz = a->sz + b->sz;
        adj = 1;
    }

    if (!adj) return 0;

    /* Check merged size doesn't exceed 64 coefficients (order+1)^3 <= MAX_COEFFS */
    fit_patch(data, Nx, Ny, Nz, merged, order);
    score_patch(data, Nx, Ny, Nz, merged);

    return (merged->max_error <= tol);
}

/* One iteration of the adaptive refinement loop */
static int refine_iteration(PatchList *pl, const float *data,
                            int Nx, int Ny, int Nz, int order,
                            double tol, double merge_tol) {
    int changes = 0;

    /* Phase 1: Collect indices to split first, then split.
     * This avoids indexing issues when pl->patches reallocates. */
    int orig_count = pl->count;
    int *split_list = (int*)malloc(orig_count * sizeof(int));
    int n_split = 0;

    for (int i = 0; i < orig_count; i++) {
        if (!pl->patches[i].active) continue;
        if (pl->patches[i].max_error > tol) {
            if (pl->patches[i].sx > MIN_BLOCK_SIZE ||
                pl->patches[i].sy > MIN_BLOCK_SIZE ||
                pl->patches[i].sz > MIN_BLOCK_SIZE) {
                split_list[n_split++] = i;
            }
        }
    }

    for (int s = 0; s < n_split; s++) {
        if (split_patch(pl, split_list[s], data, Nx, Ny, Nz, order) > 0)
            changes++;
    }
    free(split_list);

    /* Phase 2: Try to merge adjacent patches with very low error.
     * Work on a snapshot of current count to avoid issues with new entries. */
    int n = pl->count;
    for (int i = 0; i < n; i++) {
        if (!pl->patches[i].active) continue;
        if (pl->patches[i].max_error > merge_tol) continue;

        for (int j = i + 1; j < n; j++) {
            if (!pl->patches[j].active) continue;
            if (pl->patches[j].max_error > merge_tol) continue;

            Patch3D merged;
            if (try_merge_patches(data, Nx, Ny, Nz,
                                  &pl->patches[i], &pl->patches[j],
                                  order, tol, &merged)) {
                pl->patches[i].active = 0;
                pl->patches[j].active = 0;
                patchlist_push(pl, &merged);
                changes++;
                break; /* Move on to next i */
            }
        }
    }

    return changes;
}

/* Full adaptive refinement loop */
static ErrorScore adaptive_refine(PatchList *pl, const float *data,
                                  float *recon, int Nx, int Ny, int Nz,
                                  int order, double target_rms,
                                  int max_iter) {
    ErrorScore best = {0};

    for (int iter = 0; iter < max_iter; iter++) {
        /* Rasterize and score */
        rasterize_patches(pl, recon, Nx, Ny, Nz);
        ErrorScore score = compute_error(data, recon, Nx, Ny, Nz, pl);

        printf("  iter %d: %d patches, %d coeffs, %.1fx compress, "
               "rms=%.6e linf=%.6e (%.4f%%)\n",
               iter, score.n_patches, score.n_coeffs_total,
               score.compression, score.rms, score.linf, score.rel_linf * 100);

        best = score;
        if (score.rms <= target_rms) {
            printf("  Converged: rms %.6e <= target %.6e\n", score.rms, target_rms);
            break;
        }

        double merge_tol = target_rms / 4.0;
        /* Use linf tolerance scaled from rms target --
         * allow max error ~ 10x the rms target */
        double split_tol = target_rms * 10.0;

        int changes = refine_iteration(pl, data, Nx, Ny, Nz, order,
                                       split_tol, merge_tol);
        if (changes == 0) {
            printf("  No further refinement possible.\n");
            break;
        }
    }

    return best;
}

/* ===================================================================
   SFA Vector Frame Format (Task 4)

   Vector frame layout within a chunk:
     uint32  magic = 0x56454333  ("VEC3")
     uint32  n_patches
     For each patch:
       uint8   fit_type   (0=polyline, 1=poly_coeffs, 2=sine, ...)
       uint8   field_term (which field: 0=phi_x, 1=phi_y, ...)
       uint16  n_coeffs
       int16   bx, by, bz   (block origin)
       int16   sx, sy, sz   (block size)
       uint8   order
       uint8   reserved
       float   max_error
       float   rms_error
       float   coeffs[n_coeffs]
   =================================================================== */

#define VEC3_MAGIC 0x56454333

/* Fit types */
#define FIT_POLYLINE     0
#define FIT_POLY_COEFFS  1
#define FIT_SINE         2
#define FIT_COSINE       3
#define FIT_CARRIER_WAVE 4
#define FIT_GAUSSIAN     5

typedef struct {
    uint32_t magic;
    uint32_t n_patches;
    uint8_t  field_term;
    uint8_t  reserved[3];
} VecFrameHeader;

typedef struct {
    uint8_t  fit_type;
    uint8_t  field_term;
    uint16_t n_coeffs;
    int16_t  bx, by, bz;
    int16_t  sx, sy, sz;
    uint8_t  order;
    uint8_t  reserved;
    float    max_error;
    float    rms_error;
    /* followed by n_coeffs floats */
} VecPatchHeader;

/* Write a vector frame to a file descriptor */
static int write_vector_frame(FILE *fp, const PatchList *pl, int field_term) {
    /* Count active patches */
    int n_active = 0;
    for (int i = 0; i < pl->count; i++)
        if (pl->patches[i].active) n_active++;

    VecFrameHeader hdr;
    hdr.magic = VEC3_MAGIC;
    hdr.n_patches = n_active;
    hdr.field_term = field_term;
    memset(hdr.reserved, 0, sizeof(hdr.reserved));
    fwrite(&hdr, sizeof(hdr), 1, fp);

    for (int i = 0; i < pl->count; i++) {
        const Patch3D *p = &pl->patches[i];
        if (!p->active) continue;

        VecPatchHeader ph;
        ph.fit_type = FIT_POLY_COEFFS;
        ph.field_term = field_term;
        ph.n_coeffs = p->n_coeffs;
        ph.bx = p->bx; ph.by = p->by; ph.bz = p->bz;
        ph.sx = p->sx; ph.sy = p->sy; ph.sz = p->sz;
        ph.order = p->order;
        ph.reserved = 0;
        ph.max_error = (float)p->max_error;
        ph.rms_error = (float)p->rms_error;
        fwrite(&ph, sizeof(ph), 1, fp);

        /* Write coefficients as float32 */
        float fc[MAX_COEFFS];
        for (int c = 0; c < p->n_coeffs; c++)
            fc[c] = (float)p->coeffs[c];
        fwrite(fc, sizeof(float), p->n_coeffs, fp);
    }

    return 0;
}

/* Read a vector frame from file */
static int read_vector_frame(FILE *fp, PatchList *pl, int *field_term_out) {
    VecFrameHeader hdr;
    if (fread(&hdr, sizeof(hdr), 1, fp) != 1) return -1;
    if (hdr.magic != VEC3_MAGIC) return -2;

    *field_term_out = hdr.field_term;
    patchlist_init(pl, hdr.n_patches);

    for (uint32_t i = 0; i < hdr.n_patches; i++) {
        VecPatchHeader ph;
        if (fread(&ph, sizeof(ph), 1, fp) != 1) return -3;

        Patch3D p = {0};
        p.bx = ph.bx; p.by = ph.by; p.bz = ph.bz;
        p.sx = ph.sx; p.sy = ph.sy; p.sz = ph.sz;
        p.order = ph.order;
        p.n_coeffs = ph.n_coeffs;
        p.max_error = ph.max_error;
        p.rms_error = ph.rms_error;
        p.active = 1;

        float fc[MAX_COEFFS];
        if ((int)fread(fc, sizeof(float), ph.n_coeffs, fp) != ph.n_coeffs) return -4;
        for (int c = 0; c < ph.n_coeffs; c++)
            p.coeffs[c] = fc[c];

        patchlist_push(pl, &p);
    }

    return 0;
}

/* ===================================================================
   Field extraction: get a flat float array for one field from SFA buffer
   =================================================================== */

static float *extract_field(void *buf, SFA *sfa, int col, int Nx, int Ny, int Nz) {
    long N3 = (long)Nx * Ny * Nz;
    float *data = (float*)malloc(N3 * sizeof(float));

    #pragma omp parallel for
    for (long i = 0; i < N3; i++)
        data[i] = col_f(buf, sfa, col, i);

    return data;
}

/* Compute derived quantities and return them as arrays.
 * Returns triple product P = phi_x * phi_y * phi_z */
static float *compute_triple_product(void *buf, SFA *sfa, int Nx, int Ny, int Nz) {
    long N3 = (long)Nx * Ny * Nz;
    float *P = (float*)malloc(N3 * sizeof(float));

    #pragma omp parallel for
    for (long i = 0; i < N3; i++)
        P[i] = col_f(buf, sfa, 0, i) * col_f(buf, sfa, 1, i) * col_f(buf, sfa, 2, i);

    return P;
}

/* Compute curl(phi) components */
static void compute_curl(void *buf, SFA *sfa, int Nx, int Ny, int Nz,
                         float *cx_out, float *cy_out, float *cz_out) {
    double dx = 2.0 * sfa->Lx / (Nx - 1);
    double idx2 = 1.0 / (2.0 * dx);
    long NN = (long)Ny * Nz;

    #pragma omp parallel for collapse(3) schedule(static)
    for (int i = 1; i < Nx-1; i++)
    for (int j = 1; j < Ny-1; j++)
    for (int k = 1; k < Nz-1; k++) {
        long idx = (long)i * NN + j * Nz + k;
        long ip = (long)(i+1)*NN + j*Nz + k;
        long im = (long)(i-1)*NN + j*Nz + k;
        long jp = (long)i*NN + (j+1)*Nz + k;
        long jm = (long)i*NN + (j-1)*Nz + k;
        long kp = (long)i*NN + j*Nz + (k+1);
        long km = (long)i*NN + j*Nz + (k-1);

        cx_out[idx] = (float)((col_f(buf,sfa,2,jp) - col_f(buf,sfa,2,jm)
                              - col_f(buf,sfa,1,kp) + col_f(buf,sfa,1,km)) * idx2);
        cy_out[idx] = (float)((col_f(buf,sfa,0,kp) - col_f(buf,sfa,0,km)
                              - col_f(buf,sfa,2,ip) + col_f(buf,sfa,2,im)) * idx2);
        cz_out[idx] = (float)((col_f(buf,sfa,1,ip) - col_f(buf,sfa,1,im)
                              - col_f(buf,sfa,0,jp) + col_f(buf,sfa,0,jm)) * idx2);
    }
}

/* ===================================================================
   Output: write reconstructed SFA file
   =================================================================== */

static int write_recon_sfa(const char *path, SFA *orig_sfa,
                           float **field_data, int n_fields,
                           const char **field_names) {
    int Nx = orig_sfa->Nx, Ny = orig_sfa->Ny, Nz = orig_sfa->Nz;
    long N3 = (long)Nx * Ny * Nz;

    SFA *out = sfa_create(path, Nx, Ny, Nz,
                          orig_sfa->Lx, orig_sfa->Ly, orig_sfa->Lz,
                          orig_sfa->dt);
    if (!out) return -1;

    /* Match the original columns for the fields we have */
    for (int c = 0; c < n_fields; c++) {
        uint8_t sem = SFA_CUSTOM, comp = c;
        /* Try to match original column metadata */
        if (c < (int)orig_sfa->n_columns) {
            sem = orig_sfa->columns[c].semantic;
            comp = orig_sfa->columns[c].component;
        }
        sfa_add_column(out, field_names[c], SFA_F32, sem, comp);
    }
    out->flags = SFA_CODEC_BSS | SFA_FLAG_STREAMING;
    sfa_finalize_header(out);

    /* Build column pointers */
    void *cols[16];
    for (int c = 0; c < n_fields; c++)
        cols[c] = field_data[c];

    sfa_write_frame(out, 0.0, cols);
    sfa_close(out);
    return 0;
}

/* Write vector frame to a standalone binary file */
static int write_vector_file(const char *path, const PatchList *pl_array,
                             int *field_terms, int n_fields,
                             int Nx, int Ny, int Nz) {
    FILE *fp = fopen(path, "wb");
    if (!fp) return -1;

    /* File header */
    char magic[8] = "VEC3D\0\0\0";
    fwrite(magic, 1, 8, fp);
    uint32_t ver = 1;
    fwrite(&ver, 4, 1, fp);
    uint32_t nf = n_fields;
    fwrite(&nf, 4, 1, fp);
    uint32_t dims[3] = {Nx, Ny, Nz};
    fwrite(dims, 4, 3, fp);

    /* Write each field's patches */
    for (int f = 0; f < n_fields; f++) {
        write_vector_frame(fp, &pl_array[f], field_terms[f]);
    }

    fclose(fp);
    return 0;
}

/* Read vector file and reconstruct */
static int read_vector_file(const char *path, PatchList *pl_array,
                            int *field_terms, int *n_fields_out,
                            int dims_out[3]) {
    FILE *fp = fopen(path, "rb");
    if (!fp) return -1;

    char magic[8];
    if (fread(magic, 1, 8, fp) != 8 || memcmp(magic, "VEC3D\0\0\0", 8) != 0) {
        fclose(fp); return -2;
    }
    uint32_t ver, nf;
    uint32_t dims[3];
    if (fread(&ver, 4, 1, fp) != 1) { fclose(fp); return -3; }
    if (fread(&nf, 4, 1, fp) != 1) { fclose(fp); return -3; }
    if (fread(dims, 4, 3, fp) != 3) { fclose(fp); return -3; }

    dims_out[0] = dims[0]; dims_out[1] = dims[1]; dims_out[2] = dims[2];
    *n_fields_out = nf;

    for (uint32_t f = 0; f < nf; f++) {
        int ft;
        int rc = read_vector_frame(fp, &pl_array[f], &ft);
        if (rc != 0) { fclose(fp); return rc; }
        field_terms[f] = ft;
    }

    fclose(fp);
    return 0;
}

/* ===================================================================
   Main entry point
   =================================================================== */

static const char *field_names[] = {
    "phi_x", "phi_y", "phi_z", "theta_x", "theta_y", "theta_z",
    "P", "curl_x", "curl_y", "curl_z"
};

int main(int argc, char **argv) {
    if (argc < 2) {
        fprintf(stderr,
            "Usage: %s input.sfa [options]\n"
            "\nOptions:\n"
            "  -frame N        Frame to process (default: 0)\n"
            "  -tol T          Target RMS error tolerance (default: 0.001)\n"
            "  -block B        Initial block size (default: 8)\n"
            "  -order O        Polynomial order 1-3 (default: 3)\n"
            "  -max_iter N     Max refinement iterations (default: 5)\n"
            "  -fields SPEC    Field selection: phi, theta, all, derived (default: all)\n"
            "  -o_recon PATH   Output reconstructed SFA (default: <input>_recon.sfa)\n"
            "  -o_vec PATH     Output vector file (default: <input>.vec3d)\n"
            "  -no_refine      Skip adaptive refinement (just initial fit)\n"
            "\nVectorizes 3D SFA field data into tensor-product polynomial patches.\n",
            argv[0]);
        return 1;
    }

    /* Parse arguments */
    char *sfa_path = argv[1];
    int frame = 0;
    double tol = 0.001;
    int block_size = 8;
    int order = 3;
    int max_iter = 5;
    char *field_spec = "all";
    char *recon_path = NULL;
    char *vec_path = NULL;
    int do_refine = 1;

    for (int i = 2; i < argc; i++) {
        if (!strcmp(argv[i], "-frame") && i+1 < argc) frame = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-tol") && i+1 < argc) tol = atof(argv[++i]);
        else if (!strcmp(argv[i], "-block") && i+1 < argc) block_size = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-order") && i+1 < argc) order = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-max_iter") && i+1 < argc) max_iter = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-fields") && i+1 < argc) field_spec = argv[++i];
        else if (!strcmp(argv[i], "-o_recon") && i+1 < argc) recon_path = argv[++i];
        else if (!strcmp(argv[i], "-o_vec") && i+1 < argc) vec_path = argv[++i];
        else if (!strcmp(argv[i], "-no_refine")) do_refine = 0;
    }

    if (order < 1) order = 1;
    if (order > 3) order = 3;
    if (block_size < MIN_BLOCK_SIZE) block_size = MIN_BLOCK_SIZE;

    /* Default output paths */
    char recon_buf[512], vec_buf[512];
    if (!recon_path) {
        snprintf(recon_buf, sizeof(recon_buf), "%.*s_recon.sfa",
                 (int)(strrchr(sfa_path, '.') ? (strrchr(sfa_path, '.') - sfa_path) : strlen(sfa_path)),
                 sfa_path);
        recon_path = recon_buf;
    }
    if (!vec_path) {
        snprintf(vec_buf, sizeof(vec_buf), "%.*s.vec3d",
                 (int)(strrchr(sfa_path, '.') ? (strrchr(sfa_path, '.') - sfa_path) : strlen(sfa_path)),
                 sfa_path);
        vec_path = vec_buf;
    }

    /* Open SFA file */
    SFA *sfa = sfa_open(sfa_path);
    if (!sfa) { fprintf(stderr, "Cannot open %s\n", sfa_path); return 1; }

    int Nx = sfa->Nx, Ny = sfa->Ny, Nz = sfa->Nz;
    long N3 = (long)Nx * Ny * Nz;

    printf("sfa_vectorize3d: %s\n", sfa_path);
    printf("  Grid: %d x %d x %d = %ld voxels\n", Nx, Ny, Nz, N3);
    printf("  Frame: %d, Block: %d^3, Order: %d (tricubic=%d coeffs)\n",
           frame, block_size, order, (order+1)*(order+1)*(order+1));
    printf("  Target RMS: %.6e, Max iterations: %d\n", tol, max_iter);
    printf("  Fields: %s\n", field_spec);

    /* Read frame data */
    void *buf = malloc(sfa->frame_bytes);
    if (!buf || sfa_read_frame(sfa, frame, buf) != 0) {
        fprintf(stderr, "Cannot read frame %d\n", frame);
        return 1;
    }

    /* Determine which fields to process */
    int process_cols[16];
    int n_process = 0;
    const char *proc_names[16];

    if (strcmp(field_spec, "phi") == 0 || strcmp(field_spec, "all") == 0) {
        for (int c = 0; c < 3 && c < (int)sfa->n_columns; c++) {
            process_cols[n_process] = c;
            proc_names[n_process] = field_names[c];
            n_process++;
        }
    }
    if (strcmp(field_spec, "theta") == 0 || strcmp(field_spec, "all") == 0) {
        for (int c = 3; c < 6 && c < (int)sfa->n_columns; c++) {
            process_cols[n_process] = c;
            proc_names[n_process] = field_names[c];
            n_process++;
        }
    }
    if (strcmp(field_spec, "derived") == 0 || strcmp(field_spec, "all") == 0) {
        /* Will be handled separately below */
    }
    /* If "phi" was specified and nothing matched (shouldn't happen), default to col 0-2 */
    if (n_process == 0 && strcmp(field_spec, "derived") != 0) {
        int ncols = (sfa->n_columns < 3) ? sfa->n_columns : 3;
        for (int c = 0; c < ncols; c++) {
            process_cols[n_process] = c;
            proc_names[n_process] = field_names[c];
            n_process++;
        }
    }

    printf("  Processing %d field(s)\n\n", n_process);

    /* Allocate reconstruction buffer */
    float *recon = (float*)malloc(N3 * sizeof(float));

    /* Arrays for patch lists (one per field) and reconstructed fields */
    PatchList *all_pl = (PatchList*)calloc(n_process, sizeof(PatchList));
    float **recon_fields = (float**)malloc(n_process * sizeof(float*));

    /* Process each field */
    for (int fi = 0; fi < n_process; fi++) {
        int col = process_cols[fi];
        printf("=== Field: %s (column %d) ===\n", proc_names[fi], col);

        /* Extract field data */
        float *data = extract_field(buf, sfa, col, Nx, Ny, Nz);

        /* Compute data range for context */
        float dmin = data[0], dmax = data[0];
        double dsum = 0;
        for (long i = 0; i < N3; i++) {
            if (data[i] < dmin) dmin = data[i];
            if (data[i] > dmax) dmax = data[i];
            dsum += data[i];
        }
        printf("  Data range: [%.6f, %.6f], mean: %.6f\n", dmin, dmax, dsum / N3);

        /* Create initial blocks and fit */
        printf("  Creating initial %d^3 blocks...\n", block_size);
        create_initial_blocks(&all_pl[fi], Nx, Ny, Nz, block_size, data, order);

        /* Rasterize and score initial fit */
        rasterize_patches(&all_pl[fi], recon, Nx, Ny, Nz);
        ErrorScore initial = compute_error(data, recon, Nx, Ny, Nz, &all_pl[fi]);
        printf("  Initial: %d patches, %d coeffs, %.1fx compress, "
               "rms=%.6e linf=%.6e (%.4f%%)\n",
               initial.n_patches, initial.n_coeffs_total,
               initial.compression, initial.rms, initial.linf,
               initial.rel_linf * 100);

        /* Adaptive refinement */
        ErrorScore final_score;
        if (do_refine && initial.rms > tol) {
            printf("  Running adaptive refinement...\n");
            final_score = adaptive_refine(&all_pl[fi], data, recon, Nx, Ny, Nz,
                                          order, tol, max_iter);
        } else {
            final_score = initial;
        }

        /* Final rasterization for output */
        rasterize_patches(&all_pl[fi], recon, Nx, Ny, Nz);

        /* Region-level error */
        RegionError reg = compute_region_error(data, recon, Nx, Ny, Nz);
        printf("  Region errors:\n");
        printf("    Core (center 50%%): rms=%.6e linf=%.6e\n", reg.core_rms, reg.core_linf);
        printf("    Edge (outer):      rms=%.6e linf=%.6e\n", reg.edge_rms, reg.edge_linf);

        /* Save reconstructed field */
        recon_fields[fi] = (float*)malloc(N3 * sizeof(float));
        memcpy(recon_fields[fi], recon, N3 * sizeof(float));

        printf("  Final: %d patches, %.1fx compression, rms=%.6e\n\n",
               final_score.n_patches, final_score.compression, final_score.rms);

        free(data);
    }

    /* Write reconstructed SFA file */
    printf("Writing reconstructed SFA: %s\n", recon_path);
    if (write_recon_sfa(recon_path, sfa, recon_fields, n_process, proc_names) != 0) {
        fprintf(stderr, "Failed to write reconstructed SFA\n");
    } else {
        printf("  Done.\n");
    }

    /* Write vector file */
    printf("Writing vector file: %s\n", vec_path);
    int *field_terms = (int*)malloc(n_process * sizeof(int));
    for (int fi = 0; fi < n_process; fi++)
        field_terms[fi] = process_cols[fi];
    if (write_vector_file(vec_path, all_pl, field_terms, n_process, Nx, Ny, Nz) != 0) {
        fprintf(stderr, "Failed to write vector file\n");
    } else {
        /* Report file sizes */
        FILE *vfp = fopen(vec_path, "rb");
        if (vfp) {
            fseek(vfp, 0, SEEK_END);
            long vec_size = ftell(vfp);
            fclose(vfp);
            printf("  Vector file: %ld bytes (%.1fx smaller than %ld raw voxels)\n",
                   vec_size, (double)(N3 * n_process * 4) / vec_size,
                   N3 * n_process);
        }
    }

    /* Verification: read back vector file and re-rasterize to confirm roundtrip */
    printf("\n=== Verification: roundtrip through vector file ===\n");
    {
        PatchList verify_pl[16];
        int verify_ft[16];
        int verify_nf = 0;
        int verify_dims[3];

        if (read_vector_file(vec_path, verify_pl, verify_ft, &verify_nf, verify_dims) == 0) {
            printf("  Read back %d fields, dims %d x %d x %d\n",
                   verify_nf, verify_dims[0], verify_dims[1], verify_dims[2]);

            for (int fi = 0; fi < verify_nf && fi < n_process; fi++) {
                float *verify_recon = (float*)malloc(N3 * sizeof(float));
                rasterize_patches(&verify_pl[fi], verify_recon, Nx, Ny, Nz);

                /* Compare with our original reconstruction */
                double max_diff = 0;
                for (long i = 0; i < N3; i++) {
                    double d = fabs(recon_fields[fi][i] - verify_recon[i]);
                    if (d > max_diff) max_diff = d;
                }
                printf("  Field %d (%s): roundtrip max diff = %.6e %s\n",
                       fi, proc_names[fi], max_diff,
                       (max_diff < 1e-5) ? "(OK)" : "(MISMATCH)");

                free(verify_recon);
                patchlist_free(&verify_pl[fi]);
            }
        } else {
            printf("  Failed to read back vector file!\n");
        }
    }

    /* Summary report */
    printf("\n=== Summary ===\n");
    printf("  Input: %s (%d x %d x %d, frame %d)\n", sfa_path, Nx, Ny, Nz, frame);
    printf("  Output (reconstructed): %s\n", recon_path);
    printf("  Output (vectors): %s\n", vec_path);
    printf("  Fields processed: %d\n", n_process);
    for (int fi = 0; fi < n_process; fi++) {
        int np = 0, nc = 0;
        for (int i = 0; i < all_pl[fi].count; i++) {
            if (all_pl[fi].patches[i].active) {
                np++;
                nc += all_pl[fi].patches[i].n_coeffs;
            }
        }
        printf("  %s: %d patches, %d coeffs, %.1fx compression\n",
               proc_names[fi], np, nc, (double)N3 / nc);
    }

    /* Cleanup */
    for (int fi = 0; fi < n_process; fi++) {
        patchlist_free(&all_pl[fi]);
        free(recon_fields[fi]);
    }
    free(all_pl);
    free(recon_fields);
    free(field_terms);
    free(recon);
    free(buf);
    sfa_close(sfa);

    return 0;
}
