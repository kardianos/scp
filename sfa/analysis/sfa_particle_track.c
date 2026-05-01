/*  sfa_particle_track.c -- Per-particle per-timestep tracking with chirality
 *
 *  Reads an SFA file frame-by-frame, detects particles via cluster analysis
 *  (flood-fill on |P| = |phi0*phi1*phi2|), and outputs per-particle-per-timestep
 *  diagnostics including chirality observables.
 *
 *  Output (TSV to stdout, one line per detected particle per timestep):
 *    t  pid  mass  cx  cy  cz  rms_r  phi_max  P_peak  E_pot
 *    H_cross  H_self  C_asym  theta_rms  curl_rms  nvox
 *
 *  Build:
 *    gcc -O3 -fopenmp -o sfa_particle_track sfa_particle_track.c -lzstd -lm
 *
 *  Usage:
 *    sfa_particle_track input.sfa [--stride N] [--t_start T] [--t_end T] [--voxel_only]
 */

#define SFA_IMPLEMENTATION
#include "../format/sfa.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>
#include <float.h>

/* ---- Physics defaults ---- */
static double MU    = -41.345;
static double KAPPA = 50.0;

/* ---- Cluster detection parameters ---- */
#define MAX_CLUSTERS       256
#define MIN_CLUSTER_VOXELS  50
#define MIN_CLUSTER_MASS    1.0
#define P_THRESHOLD_FRAC    0.01   /* cluster threshold = frac * P_max_global */

/* ---- Tracking ---- */
#define MAX_TRACKED    256

typedef struct {
    int    pid;
    double cx, cy, cz;    /* last known centroid */
    int    alive;          /* still tracked */
} TrackedParticle;

/* ---- Per-cluster diagnostics ---- */
typedef struct {
    int    id;             /* cluster label (1-based, local to this frame) */
    int    nvox;
    double mass;           /* sum |P| * dV */
    double cx, cy, cz;    /* P-weighted centroid */
    double rms_r;          /* RMS radius from centroid */
    double phi_max;        /* max |phi| in cluster */
    double P_peak;         /* max |P| in cluster */
    double E_pot;          /* sum V(P) * dV */
    double H_cross;        /* sum theta . curl(phi) * dV */
    double H_self;         /* sum phi . curl(phi) * dV */
    double C_asym;         /* sum (|curl(phi)|^2/4 - |theta|^2) * dV */
    double theta_rms;      /* RMS theta in cluster */
    double curl_rms;       /* RMS |curl(phi)| in cluster */
    int    pid;            /* assigned particle ID (after tracking) */
} ClusterDiag;

/* ---- Extract one column from raw SFA buffer as float* ---- */
static float *extract_column_f32(const void *buf, const SFA *sfa, int col_idx) {
    uint64_t N3 = sfa->N_total;
    uint64_t off = 0;
    for (int c = 0; c < col_idx; c++)
        off += N3 * sfa_dtype_size[sfa->columns[c].dtype];

    int dtype = sfa->columns[col_idx].dtype;
    const uint8_t *src = (const uint8_t*)buf + off;

    float *arr = (float*)malloc(N3 * sizeof(float));
    if (!arr) { fprintf(stderr, "malloc failed\n"); exit(1); }

    if (dtype == SFA_F32) {
        memcpy(arr, src, N3 * sizeof(float));
    } else if (dtype == SFA_F64) {
        const double *d = (const double*)src;
        for (uint64_t i = 0; i < N3; i++) arr[i] = (float)d[i];
    } else if (dtype == SFA_F16) {
        const uint16_t *h = (const uint16_t*)src;
        for (uint64_t i = 0; i < N3; i++)
            arr[i] = sfa_f16_to_f32(h[i]);
    }
    return arr;
}

/* ---- BFS flood-fill cluster detection (26-connectivity) ---- */
static int bfs_clusters_26(const float *P_arr, int N, double threshold,
                           int *labels, int *n_clusters_out)
{
    long N3 = (long)N * N * N;
    int NN = N * N;
    memset(labels, 0, N3 * sizeof(int));

    long *queue = (long*)malloc(N3 * sizeof(long));
    if (!queue) return -1;

    int n_clusters = 0;

    for (long idx = 0; idx < N3; idx++) {
        if (P_arr[idx] < threshold || labels[idx] != 0) continue;

        n_clusters++;
        if (n_clusters > MAX_CLUSTERS) { n_clusters--; break; }

        int front = 0, back = 0;
        queue[back++] = idx;
        labels[idx] = n_clusters;

        while (front < back) {
            long cur = queue[front++];
            int ci = (int)(cur / NN);
            int cj = (int)((cur / N) % N);
            int ck = (int)(cur % N);

            /* 26-connected neighbors (periodic BC) */
            for (int di = -1; di <= 1; di++)
            for (int dj = -1; dj <= 1; dj++)
            for (int dk = -1; dk <= 1; dk++) {
                if (di == 0 && dj == 0 && dk == 0) continue;
                int ni = (ci + di + N) % N;
                int nj = (cj + dj + N) % N;
                int nk = (ck + dk + N) % N;
                long nidx = (long)ni * NN + nj * N + nk;
                if (labels[nidx] == 0 && P_arr[nidx] >= threshold) {
                    labels[nidx] = n_clusters;
                    queue[back++] = nidx;
                }
            }
        }
    }

    *n_clusters_out = n_clusters;
    free(queue);
    return 0;
}

/* ---- Compute curl(phi) via central differences (periodic BC) ---- */
static void compute_curl(const float *fx, const float *fy, const float *fz,
                         float *curl_x, float *curl_y, float *curl_z,
                         int N, double inv2dx)
{
    int NN = N * N;
    long N3 = (long)N * N * N;

    #pragma omp parallel for schedule(static)
    for (long idx = 0; idx < N3; idx++) {
        int i = (int)(idx / NN);
        int j = (int)((idx / N) % N);
        int k = (int)(idx % N);

        int ip = (i + 1) % N, im = (i - 1 + N) % N;
        int jp = (j + 1) % N, jm = (j - 1 + N) % N;
        int kp = (k + 1) % N, km = (k - 1 + N) % N;

        /* Indices for neighbors along each axis */
        long ixp = (long)ip * NN + j * N + k;
        long ixm = (long)im * NN + j * N + k;
        long iyp = (long)i * NN + jp * N + k;
        long iym = (long)i * NN + jm * N + k;
        long izp = (long)i * NN + j * N + kp;
        long izm = (long)i * NN + j * N + km;

        /* curl = (dFz/dy - dFy/dz, dFx/dz - dFz/dx, dFy/dx - dFx/dy) */
        curl_x[idx] = (float)(((double)fz[iyp] - fz[iym]
                              - (double)fy[izp] + fy[izm]) * inv2dx);
        curl_y[idx] = (float)(((double)fx[izp] - fx[izm]
                              - (double)fz[ixp] + fz[ixm]) * inv2dx);
        curl_z[idx] = (float)(((double)fy[ixp] - fy[ixm]
                              - (double)fx[iyp] + fx[iym]) * inv2dx);
    }
}

/* ---- Scan JTOP chain to get real total frames (multi-JMPF safe) ---- */
static uint32_t scan_total_frames(SFA *s) {
    uint64_t jt = s->first_jtop_offset;
    uint32_t total = 0;
    FILE *fp = s->fp;

    while (jt) {
        fseek(fp, (long)jt + 12, SEEK_SET);
        uint32_t mx, cur;
        uint64_t nxt;
        if (fread(&mx, 4, 1, fp) != 1) break;
        if (fread(&cur, 4, 1, fp) != 1) break;
        if (fread(&nxt, 8, 1, fp) != 1) break;

        for (uint32_t i = 0; i < cur; i++) {
            uint64_t jo;
            uint32_t ff, fc;
            if (fread(&jo, 8, 1, fp) != 1) break;
            if (fread(&ff, 4, 1, fp) != 1) break;
            if (fread(&fc, 4, 1, fp) != 1) break;
            total += fc;
        }
        jt = nxt;
    }

    if (total > s->total_frames)
        s->total_frames = total;

    return s->total_frames;
}

/* ---- Main ---- */
int main(int argc, char **argv)
{
    if (argc < 2) {
        fprintf(stderr, "Usage: %s input.sfa [--stride N] [--t_start T] [--t_end T] [--voxel_only]\n", argv[0]);
        return 1;
    }

    const char *input_path = argv[1];
    int stride = 1;
    double t_start = -1e30;
    double t_end   =  1e30;
    int voxel_only = 0;

    for (int i = 2; i < argc; i++) {
        if (strcmp(argv[i], "--stride") == 0 && i + 1 < argc)
            stride = atoi(argv[++i]);
        else if (strcmp(argv[i], "--t_start") == 0 && i + 1 < argc)
            t_start = atof(argv[++i]);
        else if (strcmp(argv[i], "--t_end") == 0 && i + 1 < argc)
            t_end = atof(argv[++i]);
        else if (strcmp(argv[i], "--voxel_only") == 0)
            voxel_only = 1;
    }

    /* ---- Open SFA ---- */
    SFA *sfa = sfa_open(input_path);
    if (!sfa) {
        fprintf(stderr, "Failed to open %s\n", input_path);
        return 1;
    }

    /* Scan JTOP chain for real frame count (handles multi-JMPF) */
    scan_total_frames(sfa);

    int N = (int)sfa->Nx;
    double L = sfa->Lx;
    double dx = 2.0 * L / (N - 1);
    double dV = dx * dx * dx;
    double inv2dx = 1.0 / (2.0 * dx);
    long N3 = (long)N * N * N;
    int NN = N * N;

    if (sfa->n_columns < 6) {
        fprintf(stderr, "Need at least 6 columns (phi[3] + theta[3]), got %u\n", sfa->n_columns);
        sfa_close(sfa);
        return 1;
    }

    uint32_t total = sfa->total_frames;
    fprintf(stderr, "sfa_particle_track: %s  N=%d L=%.1f dx=%.4f frames=%u\n",
            input_path, N, L, dx, total);

    /* ---- Allocate frame buffer ---- */
    void *buf = malloc(sfa->frame_bytes);
    if (!buf) {
        fprintf(stderr, "Failed to allocate frame buffer (%lu bytes)\n",
                (unsigned long)sfa->frame_bytes);
        sfa_close(sfa);
        return 1;
    }

    /* ---- Allocate work arrays ---- */
    float *phi[3], *theta[3];
    float *curl_x, *curl_y, *curl_z;
    float *P_arr;
    int   *labels;

    for (int a = 0; a < 3; a++) {
        phi[a]   = (float*)malloc(N3 * sizeof(float));
        theta[a] = (float*)malloc(N3 * sizeof(float));
    }
    curl_x = (float*)malloc(N3 * sizeof(float));
    curl_y = (float*)malloc(N3 * sizeof(float));
    curl_z = (float*)malloc(N3 * sizeof(float));
    P_arr  = (float*)malloc(N3 * sizeof(float));
    labels = (int*)calloc(N3, sizeof(int));

    if (!phi[0] || !phi[1] || !phi[2] || !theta[0] || !theta[1] || !theta[2] ||
        !curl_x || !curl_y || !curl_z || !P_arr || !labels) {
        fprintf(stderr, "Failed to allocate work arrays\n");
        sfa_close(sfa);
        return 1;
    }

    /* ---- Tracking state ---- */
    TrackedParticle tracked[MAX_TRACKED];
    int n_tracked = 0;
    int next_pid = 1;
    memset(tracked, 0, sizeof(tracked));

    /* ---- Print header ---- */
    printf("t\tpid\tmass\tcx\tcy\tcz\trms_r\tphi_max\tP_peak\tE_pot\t"
           "H_cross\tH_self\tC_asym\ttheta_rms\tcurl_rms\tnvox\n");

    /* ---- Determine if file has vec frames (need sequential read) ---- */
    int has_vec_frames = 0;
    for (uint32_t fi = 0; fi < total && fi < 10; fi++) {
        uint32_t ftype = sfa_frame_type(sfa, fi);
        if (ftype == SFA_FRAME_VEC_I || ftype == SFA_FRAME_VEC_P) {
            has_vec_frames = 1;
            break;
        }
    }

    /* ---- Process frames ---- */
    /* IMPORTANT: Vec P-frames use delta encoding -- they depend on having
     * read all previous frames in order (to maintain prev_coeffs state).
     * So we must read ALL frames sequentially, but only analyze the ones
     * matching our stride/time/type filters. For pure-voxel files we can
     * skip directly. */
    int frames_processed = 0;

    for (uint32_t fi = 0; fi < total; fi++) {
        /* Check frame type */
        uint32_t ftype = sfa_frame_type(sfa, fi);

        /* Decide if this frame should be analyzed */
        int should_analyze = 1;

        /* Stride filter */
        if (fi % (uint32_t)stride != 0)
            should_analyze = 0;

        /* Voxel-only filter */
        if (voxel_only && ftype != SFA_FRAME_VOXEL && ftype != SFA_FRAME_VEC_K)
            should_analyze = 0;

        /* Time filter (get time for all frames -- cheap) */
        double t = sfa_frame_time(sfa, fi);
        if (t < t_start || t > t_end)
            should_analyze = 0;

        /* For vec files: must read every frame to maintain reconstruction state.
         * For pure-voxel files: skip non-analyzed frames entirely. */
        if (!should_analyze && !has_vec_frames)
            continue;

        /* Read frame (sfa_read_frame handles decompression and vec->voxel reconstruction) */
        if (sfa_read_frame(sfa, fi, buf) != 0) {
            if (should_analyze)
                fprintf(stderr, "  frame %u: read failed, skipping\n", fi);
            continue;
        }

        if (!should_analyze)
            continue;

        /* ---- Extract fields ---- */
        /* Instead of malloc/free per frame, copy into pre-allocated arrays */
        {
            uint64_t off = 0;
            for (int col = 0; col < 6; col++) {
                int dtype = sfa->columns[col].dtype;
                uint64_t col_bytes = N3 * sfa_dtype_size[dtype];
                const uint8_t *src = (const uint8_t*)buf + off;
                float *dst = (col < 3) ? phi[col] : theta[col - 3];

                if (dtype == SFA_F32) {
                    memcpy(dst, src, N3 * sizeof(float));
                } else if (dtype == SFA_F64) {
                    const double *d = (const double*)src;
                    for (uint64_t v = 0; v < (uint64_t)N3; v++)
                        dst[v] = (float)d[v];
                } else if (dtype == SFA_F16) {
                    const uint16_t *h = (const uint16_t*)src;
                    for (uint64_t v = 0; v < (uint64_t)N3; v++)
                        dst[v] = sfa_f16_to_f32(h[v]);
                }
                off += col_bytes;
            }
        }

        /* ---- Compute P = |phi0 * phi1 * phi2| ---- */
        double P_max_global = 0.0;
        #pragma omp parallel for schedule(static) reduction(max:P_max_global)
        for (long idx = 0; idx < N3; idx++) {
            double pv = fabs((double)phi[0][idx] * phi[1][idx] * phi[2][idx]);
            P_arr[idx] = (float)pv;
            if (pv > P_max_global) P_max_global = pv;
        }

        if (P_max_global < 1e-30) {
            /* Empty frame — nothing above noise */
            continue;
        }

        double threshold = P_THRESHOLD_FRAC * P_max_global;

        /* ---- Compute curl(phi) ---- */
        compute_curl(phi[0], phi[1], phi[2], curl_x, curl_y, curl_z, N, inv2dx);

        /* ---- BFS cluster detection (26-connectivity) ---- */
        int n_clusters_raw = 0;
        bfs_clusters_26(P_arr, N, (float)threshold, labels, &n_clusters_raw);

        /* ---- Count clusters and filter small ones ---- */
        int *sizes  = (int*)calloc(n_clusters_raw + 1, sizeof(int));
        double *masses = (double*)calloc(n_clusters_raw + 1, sizeof(double));

        for (long idx = 0; idx < N3; idx++) {
            int lab = labels[idx];
            if (lab > 0 && lab <= n_clusters_raw) {
                sizes[lab]++;
                masses[lab] += P_arr[idx] * dV;
            }
        }

        /* Remap: keep clusters with enough voxels AND mass */
        int *remap = (int*)calloc(n_clusters_raw + 1, sizeof(int));
        int n_valid = 0;
        for (int c = 1; c <= n_clusters_raw; c++) {
            if (sizes[c] >= MIN_CLUSTER_VOXELS && masses[c] >= MIN_CLUSTER_MASS)
                remap[c] = ++n_valid;
        }

        /* Apply remap to labels */
        for (long idx = 0; idx < N3; idx++) {
            int lab = labels[idx];
            if (lab > 0 && lab <= n_clusters_raw)
                labels[idx] = remap[lab];
            else
                labels[idx] = 0;
        }

        int n_clusters = n_valid;
        free(sizes);
        free(masses);
        free(remap);

        if (n_clusters == 0)
            continue;

        /* ---- Compute per-cluster diagnostics ---- */
        ClusterDiag *cd = (ClusterDiag*)calloc(n_clusters + 1, sizeof(ClusterDiag));
        double *wt = (double*)calloc(n_clusters + 1, sizeof(double));
        double *sx = (double*)calloc(n_clusters + 1, sizeof(double));
        double *sy = (double*)calloc(n_clusters + 1, sizeof(double));
        double *sz = (double*)calloc(n_clusters + 1, sizeof(double));

        /* --- Pass 1: accumulate centroid, mass, P_peak, phi_max, E_pot,
         *             H_cross, H_self, C_asym, theta_rms, curl_rms --- */
        for (long idx = 0; idx < N3; idx++) {
            int lab = labels[idx];
            if (lab <= 0 || lab > n_clusters) continue;

            int i = (int)(idx / NN);
            int j = (int)((idx / N) % N);
            int k = (int)(idx % N);
            double x = -L + i * dx;
            double y = -L + j * dx;
            double z = -L + k * dx;

            double p0 = (double)phi[0][idx];
            double p1 = (double)phi[1][idx];
            double p2 = (double)phi[2][idx];
            double P  = (double)P_arr[idx];

            double t0 = (double)theta[0][idx];
            double t1 = (double)theta[1][idx];
            double t2 = (double)theta[2][idx];

            double cx_v = (double)curl_x[idx];
            double cy_v = (double)curl_y[idx];
            double cz_v = (double)curl_z[idx];

            /* Cluster bookkeeping */
            cd[lab].nvox++;
            cd[lab].mass += P * dV;

            /* Centroid (P-weighted) */
            sx[lab] += x * P;
            sy[lab] += y * P;
            sz[lab] += z * P;
            wt[lab] += P;

            /* P_peak */
            if (P > cd[lab].P_peak) cd[lab].P_peak = P;

            /* phi_max (max |component| in cluster) */
            double pm = fabs(p0);
            if (fabs(p1) > pm) pm = fabs(p1);
            if (fabs(p2) > pm) pm = fabs(p2);
            if (pm > cd[lab].phi_max) cd[lab].phi_max = pm;

            /* V(P) = (mu/2) * P^2 / (1 + kappa * P^2) */
            double P2 = P * P;
            cd[lab].E_pot += (MU / 2.0) * P2 / (1.0 + KAPPA * P2) * dV;

            /* H_cross = theta . curl(phi) */
            cd[lab].H_cross += (t0 * cx_v + t1 * cy_v + t2 * cz_v) * dV;

            /* H_self = phi . curl(phi) */
            cd[lab].H_self += (p0 * cx_v + p1 * cy_v + p2 * cz_v) * dV;

            /* C_asym = |curl(phi)|^2/4 - |theta|^2 */
            double curl_mag2 = cx_v * cx_v + cy_v * cy_v + cz_v * cz_v;
            double theta_mag2 = t0 * t0 + t1 * t1 + t2 * t2;
            cd[lab].C_asym += (curl_mag2 / 4.0 - theta_mag2) * dV;

            /* Accumulators for RMS */
            cd[lab].theta_rms += theta_mag2;
            cd[lab].curl_rms  += curl_mag2;
        }

        /* Compute centroids */
        for (int c = 1; c <= n_clusters; c++) {
            cd[c].id = c;
            if (wt[c] > 1e-30) {
                cd[c].cx = sx[c] / wt[c];
                cd[c].cy = sy[c] / wt[c];
                cd[c].cz = sz[c] / wt[c];
            }
        }

        /* --- Pass 2: RMS radius from centroid --- */
        double *r2sum = (double*)calloc(n_clusters + 1, sizeof(double));

        for (long idx = 0; idx < N3; idx++) {
            int lab = labels[idx];
            if (lab <= 0 || lab > n_clusters) continue;

            int i = (int)(idx / NN);
            int j = (int)((idx / N) % N);
            int k = (int)(idx % N);
            double x = -L + i * dx;
            double y = -L + j * dx;
            double z = -L + k * dx;
            double P = (double)P_arr[idx];

            double rx = x - cd[lab].cx;
            double ry = y - cd[lab].cy;
            double rz = z - cd[lab].cz;
            r2sum[lab] += (rx * rx + ry * ry + rz * rz) * P;
        }

        /* Finalize RMS quantities */
        for (int c = 1; c <= n_clusters; c++) {
            if (wt[c] > 1e-30)
                cd[c].rms_r = sqrt(r2sum[c] / wt[c]);
            if (cd[c].nvox > 0) {
                cd[c].theta_rms = sqrt(cd[c].theta_rms / cd[c].nvox);
                cd[c].curl_rms  = sqrt(cd[c].curl_rms  / cd[c].nvox);
            }
        }

        free(r2sum);

        /* ---- Track particle identity across frames ---- */
        /*  Match clusters to tracked particles by minimum centroid distance.
         *  Greedy assignment: closest pair first, no double-assignment. */

        int *assigned_pid = (int*)calloc(n_clusters + 1, sizeof(int));
        int *track_used   = (int*)calloc(MAX_TRACKED, sizeof(int));

        /* Build distance matrix (n_clusters x n_tracked) and assign greedily */
        if (n_tracked > 0 && n_clusters > 0) {
            /* Collect (cluster, track, dist) pairs and sort by distance */
            int n_pairs = n_clusters * n_tracked;
            typedef struct { int cl; int tr; double dist; } MatchPair;
            MatchPair *pairs = (MatchPair*)malloc(n_pairs * sizeof(MatchPair));
            int np = 0;

            for (int c = 1; c <= n_clusters; c++) {
                for (int tr = 0; tr < n_tracked; tr++) {
                    if (!tracked[tr].alive) continue;
                    double ddx = cd[c].cx - tracked[tr].cx;
                    double ddy = cd[c].cy - tracked[tr].cy;
                    double ddz = cd[c].cz - tracked[tr].cz;
                    pairs[np].cl = c;
                    pairs[np].tr = tr;
                    pairs[np].dist = ddx*ddx + ddy*ddy + ddz*ddz;
                    np++;
                }
            }

            /* Simple insertion sort (small n) */
            for (int a = 1; a < np; a++) {
                MatchPair tmp = pairs[a];
                int b = a - 1;
                while (b >= 0 && pairs[b].dist > tmp.dist) {
                    pairs[b + 1] = pairs[b];
                    b--;
                }
                pairs[b + 1] = tmp;
            }

            /* Greedy assignment: closest unmatched pair */
            int *cl_taken = (int*)calloc(n_clusters + 1, sizeof(int));
            for (int p = 0; p < np; p++) {
                int c = pairs[p].cl;
                int tr = pairs[p].tr;
                if (cl_taken[c] || track_used[tr]) continue;
                assigned_pid[c] = tracked[tr].pid;
                cl_taken[c] = 1;
                track_used[tr] = 1;
            }
            free(cl_taken);
            free(pairs);
        }

        /* Assign new PIDs to unmatched clusters */
        for (int c = 1; c <= n_clusters; c++) {
            if (assigned_pid[c] == 0) {
                assigned_pid[c] = next_pid++;
            }
            cd[c].pid = assigned_pid[c];
        }

        /* Update tracking state: replace tracked list with this frame's clusters */
        n_tracked = 0;
        for (int c = 1; c <= n_clusters; c++) {
            if (n_tracked >= MAX_TRACKED) break;
            tracked[n_tracked].pid   = cd[c].pid;
            tracked[n_tracked].cx    = cd[c].cx;
            tracked[n_tracked].cy    = cd[c].cy;
            tracked[n_tracked].cz    = cd[c].cz;
            tracked[n_tracked].alive = 1;
            n_tracked++;
        }

        free(assigned_pid);
        free(track_used);

        /* ---- Output ---- */
        for (int c = 1; c <= n_clusters; c++) {
            printf("%.6f\t%d\t%.6e\t%.4f\t%.4f\t%.4f\t%.4f\t%.6f\t%.6e\t%.6e\t"
                   "%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%d\n",
                   t,
                   cd[c].pid,
                   cd[c].mass,
                   cd[c].cx, cd[c].cy, cd[c].cz,
                   cd[c].rms_r,
                   cd[c].phi_max,
                   cd[c].P_peak,
                   cd[c].E_pot,
                   cd[c].H_cross,
                   cd[c].H_self,
                   cd[c].C_asym,
                   cd[c].theta_rms,
                   cd[c].curl_rms,
                   cd[c].nvox);
        }

        free(cd);
        free(wt);
        free(sx);
        free(sy);
        free(sz);

        frames_processed++;
        if (frames_processed % 50 == 0)
            fprintf(stderr, "  processed %d frames (last t=%.2f)\n", frames_processed, t);
    }

    fprintf(stderr, "sfa_particle_track: %d frames processed, %d particles seen\n",
            frames_processed, next_pid - 1);

    /* ---- Cleanup ---- */
    free(buf);
    for (int a = 0; a < 3; a++) { free(phi[a]); free(theta[a]); }
    free(curl_x); free(curl_y); free(curl_z);
    free(P_arr);
    free(labels);
    sfa_close(sfa);

    return 0;
}
