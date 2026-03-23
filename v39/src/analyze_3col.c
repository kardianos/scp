/*  analyze_3col.c — Radial profile + cluster analysis for 3-column SFA files
 *
 *  Reads phi_x, phi_y, phi_z from SFA, computes:
 *    1. Radial profiles of rho(r), P(r), V(r) per frame
 *    2. Connected-component BFS cluster detection on |P| > threshold
 *    3. Shell-integrated binding energy concentration
 *
 *  Build: gcc -O3 -march=native -fopenmp -o analyze_3col src/analyze_3col.c -lzstd -lm
 *  Usage: ./analyze_3col <file.sfa> [-threshold T] [-every N]
 */

#define SFA_IMPLEMENTATION
#include "../../sfa/format/sfa.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>

#define PI 3.14159265358979323846
#define MU    (-41.345)
#define KAPPA  50.0

#define N_SHELLS 40
#define MAX_CLUSTERS 256
#define MIN_CLUSTER_SIZE 8

/* ---- Potential ---- */
static inline double V_pot(double P) {
    double P2 = P * P;
    return (MU / 2.0) * P2 / (1.0 + KAPPA * P2);
}

/* ---- BFS Queue ---- */
typedef struct { long *data; long head, tail, cap; } Queue;
static void queue_init(Queue *q, long cap) { q->data=malloc(cap*sizeof(long)); q->head=q->tail=0; q->cap=cap; }
static void queue_free(Queue *q) { free(q->data); }
static inline void queue_push(Queue *q, long v) { q->data[q->tail++]=v; if(q->tail>=q->cap) q->tail=0; }
static inline long queue_pop(Queue *q) { long v=q->data[q->head++]; if(q->head>=q->cap) q->head=0; return v; }
static inline int queue_empty(Queue *q) { return q->head==q->tail; }

/* ---- Cluster ---- */
typedef struct {
    int id; long count;
    double mass, volume;
    double cx, cy, cz, rms_r;
    double max_P, E_pot;
} Cluster;

/* ---- Main ---- */
int main(int argc, char **argv) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <file.sfa> [-threshold T] [-every N]\n", argv[0]);
        return 1;
    }

    const char *sfa_path = argv[1];
    double threshold = -1;  /* auto */
    int every = 1;

    for (int i = 2; i < argc; i++) {
        if (!strcmp(argv[i], "-threshold") && i+1 < argc) threshold = atof(argv[++i]);
        else if (!strcmp(argv[i], "-every") && i+1 < argc) every = atoi(argv[++i]);
    }

    SFA *sfa = sfa_open(sfa_path);
    if (!sfa) { fprintf(stderr, "Failed to open %s\n", sfa_path); return 1; }

    int N = sfa->Nx;
    long N3 = (long)N * N * N;
    int NN = N * N;
    double L = sfa->Lx;
    double dx = 2.0 * L / (N - 1);
    double dV = dx * dx * dx;
    double R_max = L;
    double dr = R_max / N_SHELLS;

    printf("SFA: %s\n", sfa_path);
    printf("  Grid: %d^3 = %ld points\n", N, N3);
    printf("  Domain: [-%.1f, %.1f], dx=%.4f\n", L, L, dx);
    printf("  Columns: %d, Frames: %d\n", sfa->n_columns, sfa->total_frames);
    printf("  frame_bytes: %lu\n\n", (unsigned long)sfa->frame_bytes);

    if (sfa->n_columns != 3) {
        fprintf(stderr, "Expected 3 columns, got %d\n", sfa->n_columns);
        sfa_close(sfa);
        return 1;
    }

    /* Allocate frame buffer */
    void *buf = malloc(sfa->frame_bytes);
    int *labels = calloc(N3, sizeof(int));
    Queue bfs_q;
    queue_init(&bfs_q, N3/4 + 1024);

    /* Print headers */
    printf("=== RADIAL PROFILES ===\n");
    printf("frame\tt\tshell\tr_mid\trho2\tP_mean\tVpot\tvol\n");

    printf("\n=== CLUSTER ANALYSIS ===\n");
    printf("frame\tt\tn_clusters\tmax_|P|\ttotal_bind\t"
           "cl0_count\tcl0_cx\tcl0_cy\tcl0_cz\tcl0_rms\tcl0_maxP\tcl0_Epot\t"
           "cl1_count\tcl1_cx\tcl1_cy\tcl1_cz\n");

    printf("\n=== FRAME SUMMARIES ===\n");
    printf("frame\tt\tP_max\tP_int\tE_pot\t"
           "rho2_core\trho2_outer\t"
           "bind_core\tbind_outer\t"
           "R_half\tn_clusters\n");

    /* Process frames */
    for (uint32_t fi = 0; fi < sfa->total_frames; fi += every) {
        double t = sfa_frame_time(sfa, fi);
        if (sfa_read_frame(sfa, fi, buf) < 0) {
            fprintf(stderr, "  Frame %d: read error\n", fi);
            continue;
        }

        /* Extract 3 float columns */
        float *phi[3];
        for (int a = 0; a < 3; a++)
            phi[a] = (float*)((char*)buf + a * N3 * sizeof(float));

        /* ---- 1. Radial profiles ---- */
        double rho_shell[N_SHELLS], P_shell[N_SHELLS], Vpot_shell[N_SHELLS];
        double vol_shell[N_SHELLS];
        memset(rho_shell, 0, sizeof(rho_shell));
        memset(P_shell, 0, sizeof(P_shell));
        memset(Vpot_shell, 0, sizeof(Vpot_shell));
        memset(vol_shell, 0, sizeof(vol_shell));

        double P_max_global = 0, P_int_global = 0, Epot_global = 0;

        for (long idx = 0; idx < N3; idx++) {
            int i = (int)(idx / NN), j = (int)((idx / N) % N), k = (int)(idx % N);
            double x = -L + i * dx, y = -L + j * dx, z = -L + k * dx;
            double r = sqrt(x*x + y*y + z*z);

            double p0 = phi[0][idx], p1 = phi[1][idx], p2 = phi[2][idx];
            double rho2 = p0*p0 + p1*p1 + p2*p2;
            double P = p0 * p1 * p2;
            double absP = fabs(P);
            double vp = V_pot(P);

            if (absP > P_max_global) P_max_global = absP;
            P_int_global += absP * dV;
            Epot_global += vp * dV;

            int shell = (int)(r / dr);
            if (shell < N_SHELLS) {
                rho_shell[shell] += rho2 * dV;
                P_shell[shell] += P * dV;
                Vpot_shell[shell] += vp * dV;
                vol_shell[shell] += dV;
            }
        }

        /* Print radial profile */
        printf("\n");
        for (int s = 0; s < N_SHELLS; s++) {
            if (vol_shell[s] < 1e-20) continue;
            double r_mid = (s + 0.5) * dr;
            printf("RAD\t%d\t%.2f\t%d\t%.3f\t%.6e\t%.6e\t%.6e\t%.4f\n",
                   fi, t, s, r_mid,
                   rho_shell[s]/vol_shell[s],
                   P_shell[s]/vol_shell[s],
                   Vpot_shell[s]/vol_shell[s],
                   vol_shell[s]);
        }

        /* Core/outer split (at r = L/3) */
        double rho2_core = 0, rho2_outer = 0;
        double bind_core = 0, bind_outer = 0;
        double R_core = L / 3.0;
        int core_shells = (int)(R_core / dr);
        for (int s = 0; s < N_SHELLS; s++) {
            if (s < core_shells) {
                rho2_core += rho_shell[s];
                bind_core += Vpot_shell[s];
            } else {
                rho2_outer += rho_shell[s];
                bind_outer += Vpot_shell[s];
            }
        }

        /* Half-energy radius */
        double cum = 0, R_half = R_max;
        for (int s = 0; s < N_SHELLS; s++) {
            cum += fabs(Vpot_shell[s]);
            if (cum >= 0.5 * fabs(Epot_global)) { R_half = (s + 0.5) * dr; break; }
        }

        /* ---- 2. Cluster detection (BFS on |P| > threshold) ---- */
        double thr = threshold;
        if (thr < 0) thr = 0.01 * P_max_global;  /* auto */

        memset(labels, 0, N3 * sizeof(int));
        int n_clusters = 0;
        Cluster clusters[MAX_CLUSTERS];
        memset(clusters, 0, sizeof(clusters));

        for (long idx = 0; idx < N3; idx++) {
            double P = fabs(phi[0][idx] * phi[1][idx] * phi[2][idx]);
            if (P < thr || labels[idx] != 0) continue;

            n_clusters++;
            if (n_clusters > MAX_CLUSTERS) break;
            int cid = n_clusters;
            Cluster *cl = &clusters[cid - 1];
            cl->id = cid;

            /* BFS */
            bfs_q.head = bfs_q.tail = 0;
            labels[idx] = cid;
            queue_push(&bfs_q, idx);

            while (!queue_empty(&bfs_q)) {
                long cur = queue_pop(&bfs_q);
                int ci = (int)(cur / NN), cj = (int)((cur / N) % N), ck = (int)(cur % N);
                double cx = -L + ci * dx, cy = -L + cj * dx, cz = -L + ck * dx;

                double pp = phi[0][cur] * phi[1][cur] * phi[2][cur];
                double abspp = fabs(pp);
                double rho2 = phi[0][cur]*phi[0][cur] + phi[1][cur]*phi[1][cur] + phi[2][cur]*phi[2][cur];

                cl->count++;
                cl->mass += rho2 * dV;
                cl->cx += cx * rho2 * dV;
                cl->cy += cy * rho2 * dV;
                cl->cz += cz * rho2 * dV;
                if (abspp > cl->max_P) cl->max_P = abspp;
                cl->E_pot += V_pot(pp) * dV;

                /* 6-connected neighbors */
                int nbr[6][3] = {
                    {ci+1,cj,ck},{ci-1,cj,ck},
                    {ci,cj+1,ck},{ci,cj-1,ck},
                    {ci,cj,ck+1},{ci,cj,ck-1}
                };
                for (int n = 0; n < 6; n++) {
                    int ni = nbr[n][0], nj = nbr[n][1], nk = nbr[n][2];
                    if (ni < 0 || ni >= N || nj < 0 || nj >= N || nk < 0 || nk >= N) continue;
                    long nidx = (long)ni * NN + nj * N + nk;
                    if (labels[nidx] != 0) continue;
                    double nP = fabs(phi[0][nidx] * phi[1][nidx] * phi[2][nidx]);
                    if (nP < thr) continue;
                    labels[nidx] = cid;
                    queue_push(&bfs_q, nidx);
                }
            }

            /* Finalize centroid */
            if (cl->mass > 1e-30) {
                cl->cx /= cl->mass;
                cl->cy /= cl->mass;
                cl->cz /= cl->mass;
            }
        }

        /* Compute RMS radius for each cluster */
        for (long idx = 0; idx < N3; idx++) {
            if (labels[idx] <= 0) continue;
            Cluster *cl = &clusters[labels[idx] - 1];
            int ci = (int)(idx / NN), cj = (int)((idx / N) % N), ck = (int)(idx % N);
            double cx = -L + ci * dx - cl->cx;
            double cy = -L + cj * dx - cl->cy;
            double cz = -L + ck * dx - cl->cz;
            double rho2 = phi[0][idx]*phi[0][idx] + phi[1][idx]*phi[1][idx] + phi[2][idx]*phi[2][idx];
            cl->rms_r += (cx*cx + cy*cy + cz*cz) * rho2 * dV;
        }
        for (int c = 0; c < n_clusters; c++) {
            if (clusters[c].mass > 1e-30)
                clusters[c].rms_r = sqrt(clusters[c].rms_r / clusters[c].mass);
        }

        /* Sort clusters by count descending */
        for (int i = 0; i < n_clusters - 1; i++)
            for (int j = i + 1; j < n_clusters; j++)
                if (clusters[j].count > clusters[i].count) {
                    Cluster tmp = clusters[i]; clusters[i] = clusters[j]; clusters[j] = tmp;
                }

        /* Filter small clusters */
        int n_real = 0;
        double total_bind = 0;
        for (int c = 0; c < n_clusters; c++) {
            if (clusters[c].count >= MIN_CLUSTER_SIZE) {
                n_real++;
                total_bind += clusters[c].E_pot;
            }
        }

        /* Print cluster info */
        printf("CLU\t%d\t%.2f\t%d\t%.6e\t%.6e",
               fi, t, n_real, P_max_global, total_bind);
        for (int c = 0; c < 2 && c < n_real; c++) {
            printf("\t%ld\t%.3f\t%.3f\t%.3f\t%.3f\t%.4e\t%.4e",
                   clusters[c].count,
                   clusters[c].cx, clusters[c].cy, clusters[c].cz,
                   clusters[c].rms_r,
                   clusters[c].max_P, clusters[c].E_pot);
        }
        printf("\n");

        /* Frame summary */
        printf("SUM\t%d\t%.2f\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.3f\t%d\n",
               fi, t, P_max_global, P_int_global, Epot_global,
               rho2_core, rho2_outer, bind_core, bind_outer, R_half, n_real);

        fflush(stdout);
    }

    queue_free(&bfs_q);
    free(labels);
    free(buf);
    sfa_close(sfa);
    return 0;
}
