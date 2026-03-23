/*  sfa_frag.c — Braid fragmentation detector for SFA files
 *
 *  Reads 6-field Cosserat braid SFA files (phi_0..2, theta_0..2),
 *  computes the triple product P = phi_0 * phi_1 * phi_2, and runs
 *  connected-component labeling to detect when a braid structure
 *  fragments into separate clusters.
 *
 *  Build: gcc -O3 -march=native -fopenmp -o sfa_frag src/sfa_frag.c -lzstd -lm
 *  Usage: ./sfa_frag <file.sfa> [-threshold T] [-every N] [-o output.tsv]
 */

#define SFA_IMPLEMENTATION
#include "../../sfa/format/sfa.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>

/* ----------------------------------------------------------------
   Physics constants (must match simulation)
   ---------------------------------------------------------------- */
#define MU     (-41.345)
#define KAPPA  50.0

/* Minimum cluster size (grid points) to filter noise */
#define MIN_CLUSTER_SIZE 8

/* Maximum clusters we track per frame */
#define MAX_CLUSTERS 256

/* ----------------------------------------------------------------
   Cluster data
   ---------------------------------------------------------------- */
typedef struct {
    int    id;
    long   count;          /* number of grid points */
    double mass;           /* integral of rho * dV */
    double volume;         /* count * dV */
    double cx, cy, cz;    /* centroid (rho-weighted) */
    double rms_r;          /* RMS distance from centroid */
    double max_P;          /* peak |P| in cluster */
    double E_pot;          /* integral of V(P) over cluster */
} Cluster;

/* ----------------------------------------------------------------
   BFS queue (simple ring buffer)
   ---------------------------------------------------------------- */
typedef struct {
    long *data;
    long  head, tail, cap;
} Queue;

static void queue_init(Queue *q, long cap) {
    q->data = (long *)malloc(cap * sizeof(long));
    q->head = q->tail = 0;
    q->cap  = cap;
}

static void queue_free(Queue *q) {
    free(q->data);
}

static inline void queue_push(Queue *q, long val) {
    q->data[q->tail++] = val;
    if (q->tail >= q->cap) q->tail = 0;
}

static inline long queue_pop(Queue *q) {
    long val = q->data[q->head++];
    if (q->head >= q->cap) q->head = 0;
    return val;
}

static inline int queue_empty(Queue *q) {
    return q->head == q->tail;
}

/* ----------------------------------------------------------------
   Potential V(P) = (mu/2) * P^2 / (1 + kappa * P^2)
   ---------------------------------------------------------------- */
static inline double V_pot(double P) {
    double P2 = P * P;
    return (MU / 2.0) * P2 / (1.0 + KAPPA * P2);
}

/* ----------------------------------------------------------------
   Usage
   ---------------------------------------------------------------- */
static void usage(const char *prog) {
    fprintf(stderr,
        "Usage: %s <file.sfa> [-threshold T] [-every N] [-o output.tsv]\n"
        "\n"
        "  -threshold T  |P| threshold for active region (default: auto = 0.01 * max|P|)\n"
        "  -every N      analyze every Nth frame (default: 1)\n"
        "  -o FILE       output TSV file (default: stdout)\n",
        prog);
}

/* ----------------------------------------------------------------
   Main
   ---------------------------------------------------------------- */
int main(int argc, char **argv) {
    if (argc < 2) { usage(argv[0]); return 1; }

    /* Parse arguments */
    const char *sfa_path   = NULL;
    const char *out_path   = NULL;
    double      threshold  = -1.0;   /* negative = auto */
    int         every      = 1;

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-threshold") == 0 && i + 1 < argc) {
            threshold = atof(argv[++i]);
        } else if (strcmp(argv[i], "-every") == 0 && i + 1 < argc) {
            every = atoi(argv[++i]);
            if (every < 1) every = 1;
        } else if (strcmp(argv[i], "-o") == 0 && i + 1 < argc) {
            out_path = argv[++i];
        } else if (argv[i][0] != '-') {
            sfa_path = argv[i];
        } else {
            fprintf(stderr, "Unknown option: %s\n", argv[i]);
            usage(argv[0]);
            return 1;
        }
    }

    if (!sfa_path) {
        fprintf(stderr, "Error: no SFA file specified\n");
        usage(argv[0]);
        return 1;
    }

    /* Open SFA */
    SFA *sfa = sfa_open(sfa_path);
    if (!sfa) { fprintf(stderr, "Error: cannot open %s\n", sfa_path); return 1; }

    uint32_t Nx = sfa->Nx, Ny = sfa->Ny, Nz = sfa->Nz;
    uint64_t N3 = sfa->N_total;
    double Lx = sfa->Lx, Ly = sfa->Ly, Lz = sfa->Lz;
    double dx = 2.0 * Lx / (Nx - 1);
    double dy = 2.0 * Ly / (Ny - 1);
    double dz = 2.0 * Lz / (Nz - 1);
    double dV = dx * dy * dz;

    fprintf(stderr, "SFA: %s\n", sfa_path);
    fprintf(stderr, "  Grid: %u x %u x %u = %lu points\n",
            Nx, Ny, Nz, (unsigned long)N3);
    fprintf(stderr, "  Domain: [%.2f, %.2f] x [%.2f, %.2f] x [%.2f, %.2f]\n",
            -Lx, Lx, -Ly, Ly, -Lz, Lz);
    fprintf(stderr, "  dx=%.4f  dy=%.4f  dz=%.4f  dV=%.6e\n", dx, dy, dz, dV);
    fprintf(stderr, "  Columns: %u, Frames: %u\n", sfa->n_columns, sfa->total_frames);
    fprintf(stderr, "  Frame bytes (uncompressed): %lu\n",
            (unsigned long)sfa->frame_bytes);

    if (sfa->n_columns < 6) {
        fprintf(stderr, "Error: expected >= 6 columns (phi0..2, theta0..2), got %u\n",
                sfa->n_columns);
        sfa_close(sfa);
        return 1;
    }

    /* Print column names */
    fprintf(stderr, "  Column layout:\n");
    for (uint32_t c = 0; c < sfa->n_columns; c++) {
        fprintf(stderr, "    [%u] %-12s %s\n", c, sfa->columns[c].name,
                sfa_dtype_name(sfa->columns[c].dtype));
    }

    /* Allocate frame buffer */
    void *frame_buf = malloc(sfa->frame_bytes);
    if (!frame_buf) {
        fprintf(stderr, "Error: cannot allocate frame buffer (%lu bytes)\n",
                (unsigned long)sfa->frame_bytes);
        sfa_close(sfa);
        return 1;
    }

    /* Allocate work arrays */
    double *P_arr   = (double *)malloc(N3 * sizeof(double));
    double *rho_arr = (double *)malloc(N3 * sizeof(double));
    int    *labels  = (int    *)malloc(N3 * sizeof(int));
    if (!P_arr || !rho_arr || !labels) {
        fprintf(stderr, "Error: cannot allocate work arrays\n");
        sfa_close(sfa);
        return 1;
    }

    /* BFS queue -- worst case all points */
    Queue bfs;
    queue_init(&bfs, (long)N3 + 16);

    /* Open output */
    FILE *out_fp = stdout;
    if (out_path) {
        out_fp = fopen(out_path, "w");
        if (!out_fp) {
            fprintf(stderr, "Error: cannot open %s for writing\n", out_path);
            sfa_close(sfa);
            return 1;
        }
    }

    /* TSV header */
    fprintf(out_fp,
        "frame\ttime\tn_clusters\tcluster_id\tmass\tvolume\tcx\tcy\tcz\trms_r\tmax_P\tE_pot\n");

    /* Summary tracking */
    int max_clusters_seen = 0;
    int total_frames_analyzed = 0;

    /* Cluster count history (one entry per analyzed frame) */
    int    history_cap   = (int)((sfa->total_frames + every - 1) / every) + 1;
    int    *hist_frame   = (int    *)malloc(history_cap * sizeof(int));
    double *hist_time    = (double *)malloc(history_cap * sizeof(double));
    int    *hist_count   = (int    *)malloc(history_cap * sizeof(int));
    int     history_len  = 0;

    /* Precompute strides (constant across frames) */
    long stride_i = (long)Ny * Nz;
    long stride_j = (long)Nz;

    /* ----------------------------------------------------------------
       Frame loop
       ---------------------------------------------------------------- */
    for (uint32_t fi = 0; fi < sfa->total_frames; fi += (uint32_t)every) {
        int rc = sfa_read_frame(sfa, fi, frame_buf);
        if (rc < 0) {
            fprintf(stderr, "Warning: failed to read frame %u (rc=%d), skipping\n",
                    fi, rc);
            continue;
        }

        double frame_time = sfa_frame_time(sfa, fi);

        /* Extract phi_0, phi_1, phi_2 from frame buffer.
         * Column dtype determines element size (f32=4, f64=8). */
        int elem_sz = sfa_dtype_size[sfa->columns[0].dtype];
        int is_f32 = (sfa->columns[0].dtype == SFA_F32);

        /* Compute P = phi_0 * phi_1 * phi_2 and rho = sum(phi_a^2) */
        double max_absP = 0.0;

        #pragma omp parallel for reduction(max:max_absP) schedule(static)
        for (long n = 0; n < (long)N3; n++) {
            double p0, p1, p2;
            if (is_f32) {
                p0 = ((float *)((uint8_t *)frame_buf + 0*N3*elem_sz))[n];
                p1 = ((float *)((uint8_t *)frame_buf + 1*N3*elem_sz))[n];
                p2 = ((float *)((uint8_t *)frame_buf + 2*N3*elem_sz))[n];
            } else {
                p0 = ((double *)((uint8_t *)frame_buf + 0*N3*elem_sz))[n];
                p1 = ((double *)((uint8_t *)frame_buf + 1*N3*elem_sz))[n];
                p2 = ((double *)((uint8_t *)frame_buf + 2*N3*elem_sz))[n];
            }
            P_arr[n]   = p0 * p1 * p2;
            rho_arr[n] = p0 * p0 + p1 * p1 + p2 * p2;
            double ap  = fabs(P_arr[n]);
            if (ap > max_absP) max_absP = ap;
        }

        /* Set threshold on first frame if auto */
        double thr = threshold;
        if (thr < 0.0) {
            thr = 0.01 * max_absP;
            if (fi == 0) {
                fprintf(stderr, "  Auto threshold: %.6e (0.01 * max|P| = %.6e)\n",
                        thr, max_absP);
            }
        }

        /* Clear labels */
        #pragma omp parallel for schedule(static)
        for (long n = 0; n < (long)N3; n++)
            labels[n] = 0;

        /* ---- BFS connected-component labeling (serial) ----
         * 6-connectivity with periodic boundary wrapping. */
        int n_raw_clusters = 0;

        for (long n = 0; n < (long)N3; n++) {
            if (labels[n] != 0) continue;
            if (fabs(P_arr[n]) <= thr) continue;

            n_raw_clusters++;
            if (n_raw_clusters > MAX_CLUSTERS) {
                fprintf(stderr,
                    "Warning: frame %u exceeds %d clusters, truncating\n",
                    fi, MAX_CLUSTERS);
                n_raw_clusters = MAX_CLUSTERS;
                break;
            }

            int cid = n_raw_clusters;

            /* BFS flood fill */
            bfs.head = bfs.tail = 0;
            labels[n] = cid;
            queue_push(&bfs, n);

            while (!queue_empty(&bfs)) {
                long cur = queue_pop(&bfs);

                /* Decompose linear index: cur = i*Ny*Nz + j*Nz + k */
                long i   = cur / stride_i;
                long rem = cur - i * stride_i;
                long j   = rem / stride_j;
                long k   = rem - j * stride_j;

                /* 6-connected neighbors with periodic wrapping */
                long ni, nj, nk, nbr;

                #define TRY_NEIGHBOR(NI, NJ, NK) \
                    nbr = (NI)*stride_i + (NJ)*stride_j + (NK); \
                    if (labels[nbr] == 0 && fabs(P_arr[nbr]) > thr) { \
                        labels[nbr] = cid; queue_push(&bfs, nbr); \
                    }

                ni = (i + 1 < (long)Nx) ? i + 1 : 0;
                TRY_NEIGHBOR(ni, j, k);

                ni = (i > 0) ? i - 1 : (long)Nx - 1;
                TRY_NEIGHBOR(ni, j, k);

                nj = (j + 1 < (long)Ny) ? j + 1 : 0;
                TRY_NEIGHBOR(i, nj, k);

                nj = (j > 0) ? j - 1 : (long)Ny - 1;
                TRY_NEIGHBOR(i, nj, k);

                nk = (k + 1 < (long)Nz) ? k + 1 : 0;
                TRY_NEIGHBOR(i, j, nk);

                nk = (k > 0) ? k - 1 : (long)Nz - 1;
                TRY_NEIGHBOR(i, j, nk);

                #undef TRY_NEIGHBOR
            }
        }

        /* ---- Post-BFS: filter small clusters, compute per-cluster stats ----
         *
         * Strategy: count points per raw label, build a map from raw label
         * to compacted cluster index (skipping clusters with < MIN_CLUSTER_SIZE),
         * then do two passes over the grid:
         *   Pass 1: accumulate mass, weighted centroid, max_P, E_pot
         *   Pass 2: RMS radius from centroid (needs centroid from pass 1) */

        /* Count points per raw label */
        long label_count[MAX_CLUSTERS + 1];
        memset(label_count, 0, sizeof(label_count));
        for (long n = 0; n < (long)N3; n++) {
            int lab = labels[n];
            if (lab > 0 && lab <= MAX_CLUSTERS)
                label_count[lab]++;
        }

        /* Build map: raw label -> compacted index (0-based), or -1 if filtered */
        int label_to_new[MAX_CLUSTERS + 1];
        int n_clusters = 0;
        for (int lab = 1; lab <= n_raw_clusters; lab++) {
            if (label_count[lab] >= MIN_CLUSTER_SIZE) {
                label_to_new[lab] = n_clusters;
                n_clusters++;
            } else {
                label_to_new[lab] = -1;
            }
        }

        /* Initialize cluster structs */
        Cluster clusters[MAX_CLUSTERS];
        for (int c = 0; c < n_clusters; c++) {
            memset(&clusters[c], 0, sizeof(Cluster));
            clusters[c].id = c + 1;
        }

        /* Pass 1: mass, centroid (rho-weighted), max_P, E_pot */
        for (long n = 0; n < (long)N3; n++) {
            int lab = labels[n];
            if (lab <= 0 || lab > MAX_CLUSTERS) continue;
            int ci = label_to_new[lab];
            if (ci < 0) continue;

            Cluster *cl = &clusters[ci];

            long i   = n / stride_i;
            long rem = n - i * stride_i;
            long j   = rem / stride_j;
            long k   = rem - j * stride_j;

            double x = -Lx + i * dx;
            double y = -Ly + j * dy;
            double z = -Lz + k * dz;

            double rho_val = rho_arr[n];
            double P_val   = P_arr[n];

            cl->count++;
            cl->mass  += rho_val * dV;
            cl->E_pot += V_pot(P_val) * dV;
            cl->cx    += rho_val * x * dV;
            cl->cy    += rho_val * y * dV;
            cl->cz    += rho_val * z * dV;

            double ap = fabs(P_val);
            if (ap > cl->max_P) cl->max_P = ap;
        }

        /* Normalize centroids and compute volume */
        for (int c = 0; c < n_clusters; c++) {
            if (clusters[c].mass > 0.0) {
                clusters[c].cx /= clusters[c].mass;
                clusters[c].cy /= clusters[c].mass;
                clusters[c].cz /= clusters[c].mass;
            }
            clusters[c].volume = clusters[c].count * dV;
        }

        /* Pass 2: RMS radius from centroid (rho-weighted) */
        for (long n = 0; n < (long)N3; n++) {
            int lab = labels[n];
            if (lab <= 0 || lab > MAX_CLUSTERS) continue;
            int ci = label_to_new[lab];
            if (ci < 0) continue;

            Cluster *cl = &clusters[ci];

            long i   = n / stride_i;
            long rem = n - i * stride_i;
            long j   = rem / stride_j;
            long k   = rem - j * stride_j;

            double x = -Lx + i * dx;
            double y = -Ly + j * dy;
            double z = -Lz + k * dz;

            double ddx = x - cl->cx;
            double ddy = y - cl->cy;
            double ddz = z - cl->cz;
            double r2  = ddx * ddx + ddy * ddy + ddz * ddz;

            cl->rms_r += rho_arr[n] * r2 * dV;
        }

        for (int c = 0; c < n_clusters; c++) {
            if (clusters[c].mass > 0.0)
                clusters[c].rms_r = sqrt(clusters[c].rms_r / clusters[c].mass);
        }

        /* ---- Output rows ---- */
        if (n_clusters == 0) {
            fprintf(out_fp, "%u\t%.6f\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n",
                    fi, frame_time);
        } else {
            for (int c = 0; c < n_clusters; c++) {
                Cluster *cl = &clusters[c];
                fprintf(out_fp,
                    "%u\t%.6f\t%d\t%d\t%.6e\t%.6e\t%.4f\t%.4f\t%.4f\t%.4f\t%.6e\t%.6e\n",
                    fi, frame_time, n_clusters, cl->id,
                    cl->mass, cl->volume,
                    cl->cx, cl->cy, cl->cz,
                    cl->rms_r, cl->max_P, cl->E_pot);
            }
        }

        /* Track summary */
        if (n_clusters > max_clusters_seen)
            max_clusters_seen = n_clusters;

        hist_frame[history_len] = (int)fi;
        hist_time[history_len]  = frame_time;
        hist_count[history_len] = n_clusters;
        history_len++;
        total_frames_analyzed++;

        if (fi % 10 == 0 || n_clusters > 1) {
            fprintf(stderr, "  frame %u (t=%.2f): %d cluster%s, max|P|=%.4e\n",
                    fi, frame_time, n_clusters,
                    n_clusters != 1 ? "s" : "", max_absP);
        }
    }

    /* ----------------------------------------------------------------
       Summary to stderr
       ---------------------------------------------------------------- */
    fprintf(stderr, "\n=== Fragmentation Summary ===\n");
    fprintf(stderr, "  Total frames analyzed: %d\n", total_frames_analyzed);
    fprintf(stderr, "  Max clusters in any frame: %d\n", max_clusters_seen);
    fprintf(stderr, "  Fragmentation detected: %s\n",
            max_clusters_seen > 1 ? "YES" : "NO");

    /* Compact cluster count history (only print when count changes) */
    fprintf(stderr, "  Cluster history: ");
    int prev_count = -1;
    for (int h = 0; h < history_len; h++) {
        if (hist_count[h] != prev_count) {
            if (h > 0) fprintf(stderr, ", ");
            fprintf(stderr, "t=%.1f: %d", hist_time[h], hist_count[h]);
            prev_count = hist_count[h];
        }
    }
    fprintf(stderr, "\n");

    /* Cleanup */
    if (out_fp != stdout) fclose(out_fp);
    free(hist_frame);
    free(hist_time);
    free(hist_count);
    queue_free(&bfs);
    free(labels);
    free(rho_arr);
    free(P_arr);
    free(frame_buf);
    sfa_close(sfa);

    return 0;
}
