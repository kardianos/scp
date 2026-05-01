/*  lissajous_analysis.c -- Per-cluster particle analysis for Lissajous chirality sweep
 *
 *  Reads the LAST frame of an SFA file, computes triple product P = |phi0*phi1*phi2|,
 *  finds connected clusters via BFS flood fill, and reports per-cluster metrics as TSV.
 *
 *  Output columns (one line per cluster, TSV to stdout):
 *    cluster_id  n_voxels  mass  cx  cy  cz  rms_r  bbox_dx  bbox_dy  bbox_dz  phi_max  P_peak  E_pot
 *
 *  Build:
 *    gcc -O3 -fopenmp -o lissajous_analysis lissajous_analysis.c -lzstd -lm
 *
 *  Usage:
 *    ./lissajous_analysis input.sfa [--threshold 0.001] [--header]
 */

#define SFA_IMPLEMENTATION
#include "../../sfa/format/sfa.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

/* Physics parameters */
#define MU      (-41.345)
#define KAPPA   50.0

/* Cluster detection */
#define MAX_CLUSTERS      128
#define MIN_CLUSTER_VOXELS 8
#define DEFAULT_THRESHOLD  0.001

/* ---- Extract one column from raw SFA buffer as float* ---- */
static float *extract_column_f32(void *buf, SFA *sfa, int col_idx) {
    uint64_t N3 = sfa->N_total;
    uint64_t off = 0;
    for (int c = 0; c < col_idx; c++)
        off += N3 * sfa_dtype_size[sfa->columns[c].dtype];

    int dtype = sfa->columns[col_idx].dtype;
    uint8_t *src = (uint8_t*)buf + off;

    float *arr = (float*)malloc(N3 * sizeof(float));
    if (!arr) { fprintf(stderr, "malloc failed\n"); exit(1); }

    if (dtype == SFA_F32) {
        memcpy(arr, src, N3 * sizeof(float));
    } else if (dtype == SFA_F64) {
        double *d = (double*)src;
        for (uint64_t i = 0; i < N3; i++) arr[i] = (float)d[i];
    } else if (dtype == SFA_F16) {
        uint16_t *h = (uint16_t*)src;
        for (uint64_t i = 0; i < N3; i++)
            arr[i] = sfa_f16_to_f32(h[i]);
    }
    return arr;
}

/* ---- BFS cluster detection ---- */
static int bfs_clusters(float *phi0, float *phi1, float *phi2,
                        int N, double threshold,
                        int *labels, int *n_clusters_out)
{
    long N3 = (long)N*N*N;
    int NN = N*N;
    memset(labels, 0, N3 * sizeof(int));

    long *queue = (long*)malloc(N3 * sizeof(long));
    if (!queue) return -1;

    int n_clusters = 0;

    for (long idx = 0; idx < N3; idx++) {
        double P = fabs((double)phi0[idx] * phi1[idx] * phi2[idx]);
        if (P < threshold || labels[idx] != 0) continue;

        n_clusters++;
        if (n_clusters > MAX_CLUSTERS) { n_clusters--; break; }

        int front = 0, back = 0;
        queue[back++] = idx;
        labels[idx] = n_clusters;

        while (front < back) {
            long cur = queue[front++];
            int i = (int)(cur / NN), j = (int)((cur / N) % N), k = (int)(cur % N);

            /* 6-connected neighbors (periodic BC) */
            int di[] = {1,-1,0,0,0,0};
            int dj[] = {0,0,1,-1,0,0};
            int dk[] = {0,0,0,0,1,-1};
            for (int d = 0; d < 6; d++) {
                int ni = (i+di[d]+N)%N, nj = (j+dj[d]+N)%N, nk = (k+dk[d]+N)%N;
                long nidx = (long)ni*NN + nj*N + nk;
                if (labels[nidx] == 0) {
                    double nP = fabs((double)phi0[nidx] * phi1[nidx] * phi2[nidx]);
                    if (nP >= threshold) {
                        labels[nidx] = n_clusters;
                        queue[back++] = nidx;
                    }
                }
            }
        }
    }

    *n_clusters_out = n_clusters;
    free(queue);
    return 0;
}

/* ---- Per-cluster statistics ---- */
typedef struct {
    int    id;
    int    n_voxels;
    double mass;          /* sum |P| * dV */
    double cx, cy, cz;   /* P-weighted centroid */
    double rms_r;         /* RMS radius from centroid */
    double bbox_dx, bbox_dy, bbox_dz; /* bounding box extent */
    double phi_max;       /* max |phi| in cluster (any component) */
    double P_peak;        /* max |P| in cluster */
    double E_pot;         /* sum V(P) * dV */
} ClusterStats;

int main(int argc, char **argv)
{
    if (argc < 2) {
        fprintf(stderr, "Usage: %s input.sfa [--threshold T] [--header]\n", argv[0]);
        return 1;
    }

    const char *input_path = argv[1];
    double threshold = DEFAULT_THRESHOLD;
    int print_header = 0;

    for (int i = 2; i < argc; i++) {
        if (strcmp(argv[i], "--threshold") == 0 && i+1 < argc)
            threshold = atof(argv[++i]);
        else if (strcmp(argv[i], "--header") == 0)
            print_header = 1;
    }

    /* Open SFA */
    SFA *sfa = sfa_open(input_path);
    if (!sfa) {
        fprintf(stderr, "Failed to open %s\n", input_path);
        return 1;
    }

    int N = (int)sfa->Nx;
    double L = sfa->Lx;
    double dx = 2.0 * L / (N - 1);
    double dV = dx * dx * dx;

    if (sfa->n_columns < 6) {
        fprintf(stderr, "Need at least 6 columns (phi + theta), got %u\n", sfa->n_columns);
        sfa_close(sfa);
        return 1;
    }

    /* Read last frame */
    uint32_t last_frame = sfa->total_frames - 1;
    void *buf = malloc(sfa->frame_bytes);
    if (!buf) {
        fprintf(stderr, "Failed to allocate frame buffer (%lu bytes)\n",
                (unsigned long)sfa->frame_bytes);
        sfa_close(sfa);
        return 1;
    }

    if (sfa_read_frame(sfa, last_frame, buf) != 0) {
        fprintf(stderr, "Failed to read frame %u\n", last_frame);
        free(buf);
        sfa_close(sfa);
        return 1;
    }

    /* Extract phi columns */
    float *phi[3];
    for (int a = 0; a < 3; a++)
        phi[a] = extract_column_f32(buf, sfa, a);

    free(buf);
    sfa_close(sfa);

    long N3 = (long)N*N*N;
    int NN = N*N;

    /* BFS cluster detection */
    int *labels = (int*)calloc(N3, sizeof(int));
    int n_clusters = 0;
    bfs_clusters(phi[0], phi[1], phi[2], N, threshold, labels, &n_clusters);

    /* Count cluster sizes and filter small ones */
    int *sizes = (int*)calloc(n_clusters + 1, sizeof(int));
    for (long idx = 0; idx < N3; idx++) {
        if (labels[idx] > 0 && labels[idx] <= n_clusters)
            sizes[labels[idx]]++;
    }

    /* Remap: only keep clusters with >= MIN_CLUSTER_VOXELS */
    int *remap = (int*)calloc(n_clusters + 1, sizeof(int));
    int n_valid = 0;
    for (int c = 1; c <= n_clusters; c++) {
        if (sizes[c] >= MIN_CLUSTER_VOXELS) {
            remap[c] = ++n_valid;
        }
    }
    /* Apply remap */
    for (long idx = 0; idx < N3; idx++) {
        int lab = labels[idx];
        if (lab > 0 && lab <= n_clusters)
            labels[idx] = remap[lab];
        else
            labels[idx] = 0;
    }
    n_clusters = n_valid;

    /* Compute per-cluster statistics */
    ClusterStats *cs = (ClusterStats*)calloc(n_clusters + 1, sizeof(ClusterStats));

    /* First pass: centroid (P-weighted), mass, P_peak, phi_max, E_pot */
    double *wt  = (double*)calloc(n_clusters + 1, sizeof(double));
    double *sx  = (double*)calloc(n_clusters + 1, sizeof(double));
    double *sy  = (double*)calloc(n_clusters + 1, sizeof(double));
    double *sz  = (double*)calloc(n_clusters + 1, sizeof(double));

    for (long idx = 0; idx < N3; idx++) {
        int lab = labels[idx];
        if (lab <= 0 || lab > n_clusters) continue;

        int i = (int)(idx / NN), j = (int)((idx / N) % N), k = (int)(idx % N);
        double x = -L + i*dx, y = -L + j*dx, z = -L + k*dx;

        double p0 = (double)phi[0][idx];
        double p1 = (double)phi[1][idx];
        double p2 = (double)phi[2][idx];
        double P = fabs(p0 * p1 * p2);

        /* Cluster mass */
        cs[lab].mass += P * dV;
        cs[lab].n_voxels++;

        /* Centroid accumulation */
        sx[lab] += x * P;
        sy[lab] += y * P;
        sz[lab] += z * P;
        wt[lab] += P;

        /* P_peak */
        if (P > cs[lab].P_peak) cs[lab].P_peak = P;

        /* phi_max (max absolute value of any component) */
        double pm = fabs(p0);
        if (fabs(p1) > pm) pm = fabs(p1);
        if (fabs(p2) > pm) pm = fabs(p2);
        if (pm > cs[lab].phi_max) cs[lab].phi_max = pm;

        /* V(P) = (mu/2) * P^2 / (1 + kappa * P^2) */
        double P2 = P * P;
        double Vp = (MU / 2.0) * P2 / (1.0 + KAPPA * P2);
        cs[lab].E_pot += Vp * dV;
    }

    /* Compute centroids */
    for (int c = 1; c <= n_clusters; c++) {
        cs[c].id = c;
        if (wt[c] > 1e-30) {
            cs[c].cx = sx[c] / wt[c];
            cs[c].cy = sy[c] / wt[c];
            cs[c].cz = sz[c] / wt[c];
        }
    }

    /* Second pass: RMS radius and bounding box */
    /* Track min/max coordinates per cluster */
    double *xmin = (double*)malloc((n_clusters+1) * sizeof(double));
    double *xmax = (double*)malloc((n_clusters+1) * sizeof(double));
    double *ymin = (double*)malloc((n_clusters+1) * sizeof(double));
    double *ymax = (double*)malloc((n_clusters+1) * sizeof(double));
    double *zmin = (double*)malloc((n_clusters+1) * sizeof(double));
    double *zmax = (double*)malloc((n_clusters+1) * sizeof(double));
    double *r2sum = (double*)calloc(n_clusters+1, sizeof(double));

    for (int c = 0; c <= n_clusters; c++) {
        xmin[c] = ymin[c] = zmin[c] = 1e30;
        xmax[c] = ymax[c] = zmax[c] = -1e30;
    }

    for (long idx = 0; idx < N3; idx++) {
        int lab = labels[idx];
        if (lab <= 0 || lab > n_clusters) continue;

        int i = (int)(idx / NN), j = (int)((idx / N) % N), k = (int)(idx % N);
        double x = -L + i*dx, y = -L + j*dx, z = -L + k*dx;

        double P = fabs((double)phi[0][idx] * phi[1][idx] * phi[2][idx]);

        double rx = x - cs[lab].cx;
        double ry = y - cs[lab].cy;
        double rz = z - cs[lab].cz;
        r2sum[lab] += (rx*rx + ry*ry + rz*rz) * P;

        if (x < xmin[lab]) xmin[lab] = x;
        if (x > xmax[lab]) xmax[lab] = x;
        if (y < ymin[lab]) ymin[lab] = y;
        if (y > ymax[lab]) ymax[lab] = y;
        if (z < zmin[lab]) zmin[lab] = z;
        if (z > zmax[lab]) zmax[lab] = z;
    }

    for (int c = 1; c <= n_clusters; c++) {
        if (wt[c] > 1e-30)
            cs[c].rms_r = sqrt(r2sum[c] / wt[c]);
        cs[c].bbox_dx = (xmax[c] > xmin[c]) ? (xmax[c] - xmin[c]) : 0;
        cs[c].bbox_dy = (ymax[c] > ymin[c]) ? (ymax[c] - ymin[c]) : 0;
        cs[c].bbox_dz = (zmax[c] > zmin[c]) ? (zmax[c] - zmin[c]) : 0;
    }

    /* Print header if requested */
    if (print_header) {
        printf("cluster_id\tn_voxels\tmass\tcx\tcy\tcz\trms_r\tbbox_dx\tbbox_dy\tbbox_dz\tphi_max\tP_peak\tE_pot\n");
    }

    /* Print per-cluster results */
    for (int c = 1; c <= n_clusters; c++) {
        printf("%d\t%d\t%.6e\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.6f\t%.6e\t%.6e\n",
               cs[c].id,
               cs[c].n_voxels,
               cs[c].mass,
               cs[c].cx, cs[c].cy, cs[c].cz,
               cs[c].rms_r,
               cs[c].bbox_dx, cs[c].bbox_dy, cs[c].bbox_dz,
               cs[c].phi_max,
               cs[c].P_peak,
               cs[c].E_pot);
    }

    /* If no clusters, print nothing (0 lines) */

    /* Cleanup */
    for (int a = 0; a < 3; a++) free(phi[a]);
    free(labels);
    free(sizes);
    free(remap);
    free(cs);
    free(wt); free(sx); free(sy); free(sz);
    free(xmin); free(xmax); free(ymin); free(ymax); free(zmin); free(zmax);
    free(r2sum);

    return 0;
}
