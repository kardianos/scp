/*  sfa_vectorize.c — Progressive polynomial vectorization of SFA field data
 *
 *  Takes SFA voxel data and progressively compresses it into polynomial
 *  basis vectors by merging neighboring data points:
 *
 *  Level 0: Raw voxel values (N³ points per field)
 *  Level 1: 2-point vectors (linear interpolation between neighbors)
 *  Level 2: 3-point polynomials (quadratic fit of merged pairs)
 *  Level 3: 4-point polynomials (cubic fit of merged quadratics)
 *
 *  At each level, vectors that are similar (within tolerance) are merged,
 *  reducing the representation size while tracking precision.
 *
 *  Each physical term is vectorized separately:
 *    φ_x, φ_y, φ_z, θ_x, θ_y, θ_z
 *  Plus derived quantities computed from the simulation relationships:
 *    P = φ_x·φ_y·φ_z (triple product)
 *    curl(φ)_x, curl(φ)_y, curl(φ)_z
 *    |curl(φ)|²
 *    M_x, M_y, M_z (Cosserat mismatch = curl(φ)/2 - θ)
 *
 *  Build: gcc -O3 -march=native -fopenmp -o sfa_vectorize sfa_vectorize.c -lzstd -lm
 *  Usage: ./sfa_vectorize input.sfa [-frame N] [-tol 0.01] [-max_level 3]
 */

#define SFA_IMPLEMENTATION
#include "../format/sfa.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* ===== Data types ===== */

/* A polynomial vector: coefficients for polynomial fit along a 1D segment.
 * For level k, we have (k+1) coefficients: c[0] + c[1]*t + c[2]*t² + c[3]*t³
 * where t ∈ [0, 1] parameterizes the segment. */
#define MAX_POLY_ORDER 4

typedef struct {
    double c[MAX_POLY_ORDER]; /* polynomial coefficients */
    int order;                /* 1=linear, 2=quadratic, 3=cubic */
    int start_idx;            /* starting voxel index in the original data */
    int span;                 /* number of original voxels covered */
    double max_error;         /* max approximation error within this segment */
    double rms_error;         /* RMS error within this segment */
} PolyVec;

/* A collection of polynomial vectors for one field term along one axis */
typedef struct {
    PolyVec *vecs;
    int count;
    int capacity;
} VecList;

/* Statistics for a vectorization level */
typedef struct {
    int n_vectors;
    double total_rms;
    double max_error;
    double compression_ratio; /* original_points / n_vectors */
} LevelStats;

/* ===== f16 conversion ===== */

static float f16f(uint16_t h) {
    uint16_t s = h & 0x8000; int e = (h >> 10) & 0x1F; uint16_t m = h & 0x3FF;
    if (e == 0) return 0; if (e == 31) return s ? -1e30f : 1e30f;
    float f; uint32_t x = ((uint32_t)s << 16) | ((uint32_t)(e-15+127) << 23) | ((uint32_t)m << 13);
    memcpy(&f, &x, 4); return f;
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

/* ===== VecList operations ===== */

static void veclist_init(VecList *vl, int initial_cap) {
    vl->count = 0;
    vl->capacity = initial_cap;
    vl->vecs = malloc(initial_cap * sizeof(PolyVec));
}

static void veclist_push(VecList *vl, PolyVec *v) {
    if (vl->count >= vl->capacity) {
        vl->capacity *= 2;
        vl->vecs = realloc(vl->vecs, vl->capacity * sizeof(PolyVec));
    }
    vl->vecs[vl->count++] = *v;
}

static void veclist_free(VecList *vl) {
    free(vl->vecs);
    vl->vecs = NULL;
    vl->count = vl->capacity = 0;
}

/* ===== Polynomial fitting ===== */

/* Evaluate polynomial at parameter t ∈ [0,1] */
static double poly_eval(const PolyVec *v, double t) {
    double val = v->c[0];
    double tn = t;
    for (int i = 1; i <= v->order; i++) {
        val += v->c[i] * tn;
        tn *= t;
    }
    return val;
}

/* Fit a linear (2-point) polynomial to data[i0..i1] along one axis.
 * Returns the PolyVec and computes error statistics. */
static PolyVec fit_linear(const float *data, int i0, int i1) {
    PolyVec v = {0};
    v.order = 1;
    v.start_idx = i0;
    v.span = i1 - i0 + 1;
    v.c[0] = data[i0];
    v.c[1] = data[i1] - data[i0];

    /* Compute errors */
    v.max_error = 0;
    double sum_sq = 0;
    int n = v.span;
    for (int i = 0; i < n; i++) {
        double t = (n > 1) ? (double)i / (n - 1) : 0;
        double approx = v.c[0] + v.c[1] * t;
        double err = fabs(data[i0 + i] - approx);
        if (err > v.max_error) v.max_error = err;
        sum_sq += err * err;
    }
    v.rms_error = sqrt(sum_sq / n);
    return v;
}

/* Fit a quadratic (3-point) polynomial to data[i0..i2].
 * Uses least-squares if span > 3, exact fit if span == 3. */
static PolyVec fit_quadratic(const float *data, int i0, int i2) {
    PolyVec v = {0};
    v.order = 2;
    v.start_idx = i0;
    v.span = i2 - i0 + 1;
    int n = v.span;

    /* Least squares: minimize Σ(c0 + c1*t + c2*t² - y)² */
    /* Using normal equations with t_i = i/(n-1) */
    double S[5] = {0}, Sy[3] = {0}; /* S[k] = Σt^k, Sy[k] = Σt^k*y */
    for (int i = 0; i < n; i++) {
        double t = (n > 1) ? (double)i / (n - 1) : 0;
        double y = data[i0 + i];
        double tk = 1;
        for (int k = 0; k < 5; k++) { S[k] += tk; tk *= t; }
        tk = 1;
        for (int k = 0; k < 3; k++) { Sy[k] += tk * y; tk *= t; }
    }

    /* Solve 3x3 normal equations via Cramer's rule */
    double a[3][4] = {
        {S[0], S[1], S[2], Sy[0]},
        {S[1], S[2], S[3], Sy[1]},
        {S[2], S[3], S[4], Sy[2]}
    };

    /* Gaussian elimination */
    for (int col = 0; col < 3; col++) {
        /* Pivot */
        int piv = col;
        for (int row = col+1; row < 3; row++)
            if (fabs(a[row][col]) > fabs(a[piv][col])) piv = row;
        if (piv != col) {
            for (int j = 0; j < 4; j++) {
                double tmp = a[col][j]; a[col][j] = a[piv][j]; a[piv][j] = tmp;
            }
        }
        if (fabs(a[col][col]) < 1e-15) continue;
        for (int row = col+1; row < 3; row++) {
            double f = a[row][col] / a[col][col];
            for (int j = col; j < 4; j++) a[row][j] -= f * a[col][j];
        }
    }
    /* Back substitution */
    for (int row = 2; row >= 0; row--) {
        double sum = a[row][3];
        for (int j = row+1; j < 3; j++) sum -= a[row][j] * v.c[j];
        v.c[row] = (fabs(a[row][row]) > 1e-15) ? sum / a[row][row] : 0;
    }

    /* Compute errors */
    v.max_error = 0;
    double sum_sq = 0;
    for (int i = 0; i < n; i++) {
        double t = (n > 1) ? (double)i / (n - 1) : 0;
        double approx = v.c[0] + v.c[1]*t + v.c[2]*t*t;
        double err = fabs(data[i0 + i] - approx);
        if (err > v.max_error) v.max_error = err;
        sum_sq += err * err;
    }
    v.rms_error = sqrt(sum_sq / n);
    return v;
}

/* Fit a cubic (4-point) polynomial to data[i0..i3]. */
static PolyVec fit_cubic(const float *data, int i0, int i3) {
    PolyVec v = {0};
    v.order = 3;
    v.start_idx = i0;
    v.span = i3 - i0 + 1;
    int n = v.span;

    /* Least squares: 4x4 normal equations */
    double S[7] = {0}, Sy[4] = {0};
    for (int i = 0; i < n; i++) {
        double t = (n > 1) ? (double)i / (n - 1) : 0;
        double y = data[i0 + i];
        double tk = 1;
        for (int k = 0; k < 7; k++) { S[k] += tk; tk *= t; }
        tk = 1;
        for (int k = 0; k < 4; k++) { Sy[k] += tk * y; tk *= t; }
    }

    double a[4][5] = {
        {S[0], S[1], S[2], S[3], Sy[0]},
        {S[1], S[2], S[3], S[4], Sy[1]},
        {S[2], S[3], S[4], S[5], Sy[2]},
        {S[3], S[4], S[5], S[6], Sy[3]}
    };

    for (int col = 0; col < 4; col++) {
        int piv = col;
        for (int row = col+1; row < 4; row++)
            if (fabs(a[row][col]) > fabs(a[piv][col])) piv = row;
        if (piv != col) {
            for (int j = 0; j < 5; j++) {
                double tmp = a[col][j]; a[col][j] = a[piv][j]; a[piv][j] = tmp;
            }
        }
        if (fabs(a[col][col]) < 1e-15) continue;
        for (int row = col+1; row < 4; row++) {
            double f = a[row][col] / a[col][col];
            for (int j = col; j < 5; j++) a[row][j] -= f * a[col][j];
        }
    }
    for (int row = 3; row >= 0; row--) {
        double sum = a[row][4];
        for (int j = row+1; j < 4; j++) sum -= a[row][j] * v.c[j];
        v.c[row] = (fabs(a[row][row]) > 1e-15) ? sum / a[row][row] : 0;
    }

    v.max_error = 0;
    double sum_sq = 0;
    for (int i = 0; i < n; i++) {
        double t = (n > 1) ? (double)i / (n - 1) : 0;
        double approx = v.c[0] + v.c[1]*t + v.c[2]*t*t + v.c[3]*t*t*t;
        double err = fabs(data[i0 + i] - approx);
        if (err > v.max_error) v.max_error = err;
        sum_sq += err * err;
    }
    v.rms_error = sqrt(sum_sq / n);
    return v;
}

/* ===== Progressive merging ===== */

/* Try to merge two adjacent PolyVecs at the next polynomial level.
 * Returns 1 if merge succeeded (error within tolerance), 0 otherwise. */
static int try_merge(const float *data, const PolyVec *a, const PolyVec *b,
                     int target_order, double tol, PolyVec *merged) {
    int i0 = a->start_idx;
    int i1 = b->start_idx + b->span - 1;

    switch (target_order) {
        case 1: *merged = fit_linear(data, i0, i1); break;
        case 2: *merged = fit_quadratic(data, i0, i1); break;
        case 3: *merged = fit_cubic(data, i0, i1); break;
        default: return 0;
    }

    return merged->max_error <= tol;
}

/* Vectorize a 1D data array progressively.
 * Returns stats for each level. */
static void vectorize_1d(const float *data, int n, double tol, int max_level,
                         LevelStats *stats) {
    VecList current, next;
    veclist_init(&current, n / 2 + 1);
    veclist_init(&next, n / 4 + 1);

    /* Level 1: initial linear segments (pairs of adjacent points) */
    for (int i = 0; i < n - 1; i += 2) {
        int end = (i + 1 < n) ? i + 1 : i;
        PolyVec v = fit_linear(data, i, end);
        veclist_push(&current, &v);
    }
    if (n % 2 == 1) {
        /* Odd point at the end — single-point "vector" */
        PolyVec v = {0};
        v.order = 0;
        v.start_idx = n - 1;
        v.span = 1;
        v.c[0] = data[n - 1];
        veclist_push(&current, &v);
    }

    stats[0].n_vectors = current.count;
    stats[0].compression_ratio = (double)n / current.count;
    double max_err = 0, sum_rms = 0;
    for (int i = 0; i < current.count; i++) {
        if (current.vecs[i].max_error > max_err) max_err = current.vecs[i].max_error;
        sum_rms += current.vecs[i].rms_error;
    }
    stats[0].max_error = max_err;
    stats[0].total_rms = sum_rms / current.count;

    /* Levels 2,3: progressive merging */
    for (int level = 1; level < max_level; level++) {
        int target_order = level + 1; /* level 1→quadratic(2), level 2→cubic(3) */
        if (target_order > 3) target_order = 3;

        next.count = 0;
        int i = 0;
        while (i < current.count) {
            if (i + 1 < current.count) {
                PolyVec merged;
                if (try_merge(data, &current.vecs[i], &current.vecs[i+1],
                              target_order, tol, &merged)) {
                    veclist_push(&next, &merged);
                    i += 2; /* consumed both */
                    continue;
                }
            }
            /* Can't merge — keep as-is */
            veclist_push(&next, &current.vecs[i]);
            i++;
        }

        /* Swap current and next */
        VecList tmp = current;
        current = next;
        next = tmp;
        next.count = 0;

        stats[level].n_vectors = current.count;
        stats[level].compression_ratio = (double)n / current.count;
        max_err = 0; sum_rms = 0;
        for (int j = 0; j < current.count; j++) {
            if (current.vecs[j].max_error > max_err) max_err = current.vecs[j].max_error;
            sum_rms += current.vecs[j].rms_error;
        }
        stats[level].max_error = max_err;
        stats[level].total_rms = sum_rms / current.count;
    }

    veclist_free(&current);
    veclist_free(&next);
}

/* ===== Derived quantities ===== */

/* Compute curl(φ) and related quantities from raw field data */
static void compute_derived(void *buf, SFA *sfa, int N,
                            float *P_data,      /* N³ */
                            float *curl_x,      /* N³ */
                            float *curl_y,
                            float *curl_z,
                            float *curl_sq,     /* |curl|² */
                            float *M_x,         /* mismatch */
                            float *M_y,
                            float *M_z) {
    long NN = N * N;
    double dx = 2.0 * sfa->Lx / (N - 1);
    double idx2 = 1.0 / (2.0 * dx);

    #pragma omp parallel for collapse(3) schedule(static)
    for (int i = 1; i < N-1; i++)
    for (int j = 1; j < N-1; j++)
    for (int k = 1; k < N-1; k++) {
        long idx = (long)i * NN + j * N + k;

        /* Triple product P */
        float p0 = col_f(buf, sfa, 0, idx);
        float p1 = col_f(buf, sfa, 1, idx);
        float p2 = col_f(buf, sfa, 2, idx);
        P_data[idx] = p0 * p1 * p2;

        /* Curl via central differences */
        long ip = (long)(i+1)*NN + j*N + k;
        long im = (long)(i-1)*NN + j*N + k;
        long jp = (long)i*NN + (j+1)*N + k;
        long jm = (long)i*NN + (j-1)*N + k;
        long kp = (long)i*NN + j*N + (k+1);
        long km = (long)i*NN + j*N + (k-1);

        /* curl(φ)_x = dφ_z/dy - dφ_y/dz */
        float cx = (col_f(buf,sfa,2,jp) - col_f(buf,sfa,2,jm)
                  - col_f(buf,sfa,1,kp) + col_f(buf,sfa,1,km)) * idx2;
        /* curl(φ)_y = dφ_x/dz - dφ_z/dx */
        float cy = (col_f(buf,sfa,0,kp) - col_f(buf,sfa,0,km)
                  - col_f(buf,sfa,2,ip) + col_f(buf,sfa,2,im)) * idx2;
        /* curl(φ)_z = dφ_y/dx - dφ_x/dy */
        float cz = (col_f(buf,sfa,1,ip) - col_f(buf,sfa,1,im)
                  - col_f(buf,sfa,0,jp) + col_f(buf,sfa,0,jm)) * idx2;

        curl_x[idx] = cx;
        curl_y[idx] = cy;
        curl_z[idx] = cz;
        curl_sq[idx] = cx*cx + cy*cy + cz*cz;

        /* θ values */
        float tx = col_f(buf, sfa, 3, idx);
        float ty = col_f(buf, sfa, 4, idx);
        float tz = col_f(buf, sfa, 5, idx);

        /* Mismatch M = curl(φ)/2 - θ */
        M_x[idx] = cx * 0.5f - tx;
        M_y[idx] = cy * 0.5f - ty;
        M_z[idx] = cz * 0.5f - tz;
    }
}

/* ===== Reconstruction and scoring ===== */

/* Reconstruct a 1D array from polynomial vectors */
static void reconstruct_1d(const VecList *vl, float *recon, int n) {
    memset(recon, 0, n * sizeof(float));
    for (int vi = 0; vi < vl->count; vi++) {
        const PolyVec *v = &vl->vecs[vi];
        for (int i = 0; i < v->span; i++) {
            int idx = v->start_idx + i;
            if (idx >= n) break;
            double t = (v->span > 1) ? (double)i / (v->span - 1) : 0;
            recon[idx] = (float)poly_eval(v, t);
        }
    }
}

/* Score the reconstruction quality */
typedef struct {
    double max_err;   /* L∞ */
    double rms_err;   /* L² */
    double mean_err;  /* L¹ */
    double rel_max;   /* max error / max |data| */
    int n_vectors;
    double compression;
} ReconScore;

static ReconScore score_recon(const float *orig, const float *recon, int n, int n_vectors) {
    ReconScore s = {0};
    double max_val = 0;
    double sum_sq = 0, sum_abs = 0;

    for (int i = 0; i < n; i++) {
        double err = fabs(orig[i] - recon[i]);
        if (err > s.max_err) s.max_err = err;
        sum_sq += err * err;
        sum_abs += err;
        double av = fabs(orig[i]);
        if (av > max_val) max_val = av;
    }
    s.rms_err = sqrt(sum_sq / n);
    s.mean_err = sum_abs / n;
    s.rel_max = (max_val > 1e-15) ? s.max_err / max_val : 0;
    s.n_vectors = n_vectors;
    s.compression = (double)n / n_vectors;
    return s;
}

/* Full vectorize + reconstruct + score pipeline for 1D data.
 * Returns the vector list (caller must free) and the score. */
static VecList vectorize_and_score(const float *data, int n, double tol, int max_level,
                                    ReconScore *score_out) {
    VecList current, next;
    veclist_init(&current, n / 2 + 1);
    veclist_init(&next, n / 4 + 1);

    /* Level 1: linear pairs */
    for (int i = 0; i < n - 1; i += 2) {
        PolyVec v = fit_linear(data, i, i + 1);
        veclist_push(&current, &v);
    }
    if (n % 2 == 1) {
        PolyVec v = {0};
        v.order = 0; v.start_idx = n-1; v.span = 1; v.c[0] = data[n-1];
        veclist_push(&current, &v);
    }

    /* Progressive merging */
    for (int level = 1; level < max_level; level++) {
        int target_order = level + 1;
        if (target_order > 3) target_order = 3;

        next.count = 0;
        int i = 0;
        while (i < current.count) {
            if (i + 1 < current.count) {
                PolyVec merged;
                if (try_merge(data, &current.vecs[i], &current.vecs[i+1],
                              target_order, tol, &merged)) {
                    veclist_push(&next, &merged);
                    i += 2;
                    continue;
                }
            }
            veclist_push(&next, &current.vecs[i]);
            i++;
        }
        VecList tmp = current; current = next; next = tmp; next.count = 0;
    }
    veclist_free(&next);

    /* Reconstruct and score */
    float *recon = malloc(n * sizeof(float));
    reconstruct_1d(&current, recon, n);
    *score_out = score_recon(data, recon, n, current.count);
    free(recon);

    return current;
}

/* Iterative tolerance optimizer: find the best tolerance for a target compression */
static double optimize_tolerance(const float *data, int n, int max_level,
                                  double target_compression, double target_max_err,
                                  int max_iter) {
    double tol_lo = 1e-8, tol_hi = 1.0;
    double best_tol = 0.01;
    ReconScore best_score = {0};

    for (int iter = 0; iter < max_iter; iter++) {
        double tol = sqrt(tol_lo * tol_hi); /* geometric mean */
        ReconScore score;
        VecList vl = vectorize_and_score(data, n, tol, max_level, &score);
        veclist_free(&vl);

        if (score.compression >= target_compression && score.max_err <= target_max_err) {
            /* Good — try tighter */
            best_tol = tol;
            best_score = score;
            tol_hi = tol;
        } else if (score.compression < target_compression) {
            /* Not enough compression — relax tolerance */
            tol_lo = tol;
        } else {
            /* Too much error — tighten tolerance */
            tol_hi = tol;
        }

        if (tol_hi / tol_lo < 1.01) break; /* converged */
    }

    printf("  Optimized: tol=%.6e → %d vectors (%.1f×), max_err=%.6e, rms=%.6e\n",
           best_tol, best_score.n_vectors, best_score.compression,
           best_score.max_err, best_score.rms_err);
    return best_tol;
}

/* ===== Main ===== */

int main(int argc, char **argv) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s input.sfa [-frame N] [-tol 0.01] [-max_level 3] [-axis x|y|z]\n", argv[0]);
        fprintf(stderr, "\nProgressively vectorizes SFA field data into polynomial basis vectors.\n");
        fprintf(stderr, "Analyzes each field term along the chosen axis, reporting compression\n");
        fprintf(stderr, "ratios and approximation errors at each polynomial level.\n");
        return 1;
    }

    char *sfa_path = argv[1];
    int frame = 0;
    double tol = 0.01;
    int max_level = 3;
    int axis = 0; /* 0=x, 1=y, 2=z */

    for (int i = 2; i < argc; i++) {
        if (!strcmp(argv[i], "-frame") && i+1 < argc) frame = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-tol") && i+1 < argc) tol = atof(argv[++i]);
        else if (!strcmp(argv[i], "-max_level") && i+1 < argc) max_level = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-axis") && i+1 < argc) {
            i++;
            if (argv[i][0] == 'y' || argv[i][0] == 'Y') axis = 1;
            else if (argv[i][0] == 'z' || argv[i][0] == 'Z') axis = 2;
        }
    }
    if (max_level > 3) max_level = 3;

    SFA *sfa = sfa_open(sfa_path);
    if (!sfa) { fprintf(stderr, "Cannot open %s\n", sfa_path); return 1; }

    int N = sfa->Nx;
    long N3 = (long)N * N * N;

    printf("sfa_vectorize: %s (N=%d, frame=%d)\n", sfa_path, N, frame);
    printf("  tol=%.4f max_level=%d axis=%c\n", tol, max_level, "xyz"[axis]);

    void *buf = malloc(sfa->frame_bytes);
    if (!buf || sfa_read_frame(sfa, frame, buf) != 0) {
        fprintf(stderr, "Cannot read frame %d\n", frame);
        return 1;
    }

    /* Allocate derived quantity arrays */
    float *P_data  = calloc(N3, sizeof(float));
    float *curl_x  = calloc(N3, sizeof(float));
    float *curl_y  = calloc(N3, sizeof(float));
    float *curl_z  = calloc(N3, sizeof(float));
    float *curl_sq = calloc(N3, sizeof(float));
    float *M_x     = calloc(N3, sizeof(float));
    float *M_y     = calloc(N3, sizeof(float));
    float *M_z     = calloc(N3, sizeof(float));

    printf("Computing derived quantities...\n");
    compute_derived(buf, sfa, N, P_data, curl_x, curl_y, curl_z, curl_sq, M_x, M_y, M_z);

    /* Define the terms to vectorize */
    struct {
        const char *name;
        float *data;     /* if non-NULL, use this precomputed array */
        int sfa_col;     /* if data==NULL, extract this SFA column */
    } terms[] = {
        {"phi_x",     NULL, 0},
        {"phi_y",     NULL, 1},
        {"phi_z",     NULL, 2},
        {"theta_x",   NULL, 3},
        {"theta_y",   NULL, 4},
        {"theta_z",   NULL, 5},
        {"P",         P_data, -1},
        {"curl_x",    curl_x, -1},
        {"curl_y",    curl_y, -1},
        {"curl_z",    curl_z, -1},
        {"|curl|^2",  curl_sq, -1},
        {"M_x",       M_x, -1},
        {"M_y",       M_y, -1},
        {"M_z",       M_z, -1},
    };
    int n_terms = sizeof(terms) / sizeof(terms[0]);

    /* Extract 1D slices along the chosen axis through the center */
    int mid = N / 2;

    printf("\n%-12s  %6s  %10s  %10s  %10s\n",
           "Term", "Level", "Vectors", "MaxErr", "Compress");
    printf("%-12s  %6s  %10s  %10s  %10s\n",
           "----", "-----", "-------", "------", "--------");

    int do_optimize = 0;
    double target_compress = 16.0;
    double target_err = 0.01;
    for (int i = 2; i < argc; i++) {
        if (!strcmp(argv[i], "-optimize")) do_optimize = 1;
        else if (!strcmp(argv[i], "-target_compress") && i+1 < argc) target_compress = atof(argv[++i]);
        else if (!strcmp(argv[i], "-target_err") && i+1 < argc) target_err = atof(argv[++i]);
    }

    for (int ti = 0; ti < n_terms; ti++) {
        /* Extract 1D slice */
        float *slice = malloc(N * sizeof(float));

        for (int p = 0; p < N; p++) {
            long idx;
            switch (axis) {
                case 0: idx = (long)p * N * N + mid * N + mid; break;
                case 1: idx = (long)mid * N * N + p * N + mid; break;
                default: idx = (long)mid * N * N + mid * N + p; break;
            }

            if (terms[ti].data) {
                slice[p] = terms[ti].data[idx];
            } else {
                slice[p] = col_f(buf, sfa, terms[ti].sfa_col, idx);
            }
        }

        if (do_optimize) {
            printf("%-12s: optimizing (target %.1f× compress, %.4e max err)...\n",
                   terms[ti].name, target_compress, target_err);
            double opt_tol = optimize_tolerance(slice, N, max_level,
                                                target_compress, target_err, 30);

            /* Final vectorize + score with optimal tolerance */
            ReconScore score;
            VecList vl = vectorize_and_score(slice, N, opt_tol, max_level, &score);

            printf("  Final: %d vectors (%.1f×), max=%.6e rms=%.6e rel=%.4f%%\n\n",
                   score.n_vectors, score.compression,
                   score.max_err, score.rms_err, score.rel_max * 100);
            veclist_free(&vl);
        } else {
            /* Standard level-by-level report */
            LevelStats stats[3];
            vectorize_1d(slice, N, tol, max_level, stats);

            for (int lv = 0; lv < max_level; lv++) {
                printf("%-12s  L%-5d  %10d  %10.6f  %10.1f×\n",
                       terms[ti].name, lv + 1,
                       stats[lv].n_vectors,
                       stats[lv].max_error,
                       stats[lv].compression_ratio);
            }
            printf("\n");
        }
        free(slice);
    }

    /* Cleanup */
    free(buf);
    free(P_data); free(curl_x); free(curl_y); free(curl_z);
    free(curl_sq); free(M_x); free(M_y); free(M_z);
    sfa_close(sfa);

    return 0;
}
