/*  phase_track.c — Track spatial distribution of intermediate phase group (φ≈0.7 rad)
 *
 *  For the V42 deuterium run, at T=500 a third phase group appears at
 *  atan2(phi_1, phi_0) ≈ 0.70 rad, intermediate between the two anti-phase
 *  groups (-0.79 and +2.40 rad). This program tracks where these
 *  intermediate-phase voxels are located relative to the baryon clusters.
 *
 *  Phase = atan2(phi_1, phi_0) at each active voxel (|P| > threshold).
 *  The three phase groups are:
 *    A: phase ≈ -0.79 rad
 *    B: phase ≈ +2.40 rad  (anti-phase to A, since -0.79+π ≈ 2.35)
 *    C: phase ≈ +0.70 rad  (intermediate — the "bond" candidate)
 *
 *  Build: gcc -O3 -fopenmp -o phase_track phase_track.c -lzstd -lm
 *  Usage: ./phase_track input.sfa
 */

#define SFA_IMPLEMENTATION
#include "../../sfa/format/sfa.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>

#define PI 3.14159265358979323846

/* Phase windows (atan2 in radians) */
#define PHASE_A_CENTER  (-0.79)
#define PHASE_B_CENTER  (2.40)
#define PHASE_C_CENTER  (0.70)   /* the intermediate phase */
#define PHASE_WINDOW    (0.5)    /* ±0.5 rad around center */

/* |P| threshold for "active" voxels */
#define P_THRESHOLD     0.002
/* Higher threshold for baryon cluster BFS */
#define CLUSTER_P_THRESH 0.01

/* Circular distance between two angles in [-π, π] */
static inline double angle_dist(double a, double b) {
    double d = a - b;
    while (d > PI) d -= 2*PI;
    while (d < -PI) d += 2*PI;
    return fabs(d);
}

/* BFS cluster info */
typedef struct {
    double cx, cy, cz;   /* centroid in physical coords */
    double mass;          /* total |P| */
    long   count;         /* voxel count */
    double phase;         /* atan2 phase at centroid */
} ClusterInfo;

/* Find clusters by BFS on |P|, return sorted by mass descending */
static int find_clusters_bfs(double *phi[3], int N, double L, double threshold,
                             ClusterInfo *clusters, int max_clusters) {
    long N3 = (long)N*N*N;
    int NN = N*N;
    double dx = 2.0*L/(N-1);
    char *visited = (char*)calloc(N3, 1);
    long *queue = (long*)malloc(N3 * sizeof(long));
    int n_clusters = 0;

    for (long idx = 0; idx < N3 && n_clusters < max_clusters; idx++) {
        double P = fabs(phi[0][idx] * phi[1][idx] * phi[2][idx]);
        if (P < threshold || visited[idx]) continue;

        double sx=0, sy=0, sz=0, sm=0;
        long cnt = 0;
        int front = 0, back = 0;
        queue[back++] = idx;
        visited[idx] = 1;

        while (front < back) {
            long cur = queue[front++];
            int i = (int)(cur/NN), j = (int)((cur/N)%N), k = (int)(cur%N);
            double Pv = fabs(phi[0][cur] * phi[1][cur] * phi[2][cur]);
            double x = -L + i*dx, y = -L + j*dx, z = -L + k*dx;
            sx += x*Pv; sy += y*Pv; sz += z*Pv; sm += Pv; cnt++;

            int di[] = {1,-1,0,0,0,0};
            int dj[] = {0,0,1,-1,0,0};
            int dk[] = {0,0,0,0,1,-1};
            for (int d = 0; d < 6; d++) {
                int ni = (i+di[d]+N)%N, nj = (j+dj[d]+N)%N, nk = (k+dk[d]+N)%N;
                long nidx = (long)ni*NN + nj*N + nk;
                if (!visited[nidx]) {
                    double nP = fabs(phi[0][nidx]*phi[1][nidx]*phi[2][nidx]);
                    if (nP >= threshold) {
                        visited[nidx] = 1;
                        queue[back++] = nidx;
                    }
                }
            }
        }

        if (sm > 0) {
            clusters[n_clusters].cx = sx/sm;
            clusters[n_clusters].cy = sy/sm;
            clusters[n_clusters].cz = sz/sm;
            clusters[n_clusters].mass = sm;
            clusters[n_clusters].count = cnt;
            /* Measure phase at centroid */
            int ic = (int)((clusters[n_clusters].cx + L)/dx + 0.5);
            int jc = (int)((clusters[n_clusters].cy + L)/dx + 0.5);
            int kc = (int)((clusters[n_clusters].cz + L)/dx + 0.5);
            if(ic<0)ic=0; if(ic>=N)ic=N-1;
            if(jc<0)jc=0; if(jc>=N)jc=N-1;
            if(kc<0)kc=0; if(kc>=N)kc=N-1;
            long cidx = (long)ic*NN + jc*N + kc;
            clusters[n_clusters].phase = atan2(phi[1][cidx], phi[0][cidx]);
            n_clusters++;
        }
    }

    free(visited);
    free(queue);

    /* Sort by mass descending */
    for (int i = 0; i < n_clusters-1; i++)
        for (int j = i+1; j < n_clusters; j++)
            if (clusters[j].mass > clusters[i].mass) {
                ClusterInfo tmp = clusters[i];
                clusters[i] = clusters[j];
                clusters[j] = tmp;
            }

    return n_clusters;
}

/* f16 decode helper */
static inline double f16_to_f64(uint16_t h) {
    uint16_t sign = h & 0x8000;
    int exp = (h >> 10) & 0x1F;
    uint16_t mant = h & 0x3FF;
    if (exp == 0) return 0.0;
    if (exp == 31) return sign ? -1e30 : 1e30;
    float fv;
    uint32_t x = ((uint32_t)sign << 16) | ((uint32_t)(exp-15+127) << 23) | ((uint32_t)mant << 13);
    memcpy(&fv, &x, 4);
    return (double)fv;
}

/* Extract column from frame buffer */
static double *extract_col(void *buf, SFA *sfa, int target_sem, int target_comp, long N3) {
    double *arr = (double*)calloc(N3, sizeof(double));
    uint64_t off = 0;
    for (uint32_t c = 0; c < sfa->n_columns; c++) {
        int dtype = sfa->columns[c].dtype, sem = sfa->columns[c].semantic, comp = sfa->columns[c].component;
        int es = sfa_dtype_size[dtype];
        if (sem == target_sem && comp == target_comp) {
            uint8_t *src = (uint8_t*)buf + off;
            if (dtype == SFA_F64) memcpy(arr, src, N3*8);
            else if (dtype == SFA_F32) for(long i=0;i<N3;i++) arr[i]=(double)((float*)src)[i];
            else if (dtype == SFA_F16) for(long i=0;i<N3;i++) arr[i]=f16_to_f64(((uint16_t*)src)[i]);
            return arr;
        }
        off += (uint64_t)N3 * es;
    }
    return arr;
}

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "Usage: %s input.sfa\n", argv[0]); return 1; }

    SFA *sfa = sfa_open(argv[1]);
    if (!sfa) { fprintf(stderr, "Cannot open %s\n", argv[1]); return 1; }

    int N = sfa->Nx;
    double L = sfa->Lx;
    long N3 = (long)N*N*N;
    int NN = N*N;
    double dx = 2.0*L/(N-1);
    double dV = dx*dx*dx;

    fprintf(stderr, "Opened: N=%d, L=%.1f, frames=%u, frame_bytes=%lu (%.1f GB)\n",
            N, L, sfa->total_frames, (unsigned long)sfa->frame_bytes,
            sfa->frame_bytes/1e9);
    fprintf(stderr, "Phase groups: A=%.2f, B=%.2f, C(intermed)=%.2f rad, window=±%.2f\n",
            PHASE_A_CENTER, PHASE_B_CENTER, PHASE_C_CENTER, PHASE_WINDOW);

    void *buf = malloc(sfa->frame_bytes);
    if (!buf) { fprintf(stderr, "Cannot allocate %.1f GB\n", sfa->frame_bytes/1e9); return 1; }

    /* TSV output */
    FILE *tsv = fopen("phase_tracking.tsv", "w");
    fprintf(tsv, "t\t"
                 "n_clusters\t"
                 "phaseA_count\tphaseA_vol\tphaseA_cx\tphaseA_cy\tphaseA_cz\t"
                 "phaseB_count\tphaseB_vol\tphaseB_cx\tphaseB_cy\tphaseB_cz\t"
                 "phaseC_count\tphaseC_vol\tphaseC_cx\tphaseC_cy\tphaseC_cz\t"
                 "phaseC_Pint\t"
                 "baryon1_cx\tbaryon1_cy\tbaryon1_cz\tbaryon1_phase\tbaryon1_mass\t"
                 "baryon2_cx\tbaryon2_cy\tbaryon2_cz\tbaryon2_phase\tbaryon2_mass\t"
                 "bond_mid_x\tbond_mid_y\tbond_mid_z\tbond_length\t"
                 "phaseC_bond_dist\tphaseC_frac_in_bond\n");

    /* Detailed per-frame output */
    FILE *detail = fopen("phase_detail.txt", "w");

    fprintf(stderr, "\nProcessing frames...\n");

    for (uint32_t f = 0; f < sfa->total_frames; f++) {
        double t = sfa_frame_time(sfa, f);
        fprintf(stderr, "  Frame %u (t=%.1f)...", f, t);
        double t0 = omp_get_wtime();

        sfa_read_frame(sfa, f, buf);

        double *phi[3];
        phi[0] = extract_col(buf, sfa, SFA_POSITION, 0, N3);
        phi[1] = extract_col(buf, sfa, SFA_POSITION, 1, N3);
        phi[2] = extract_col(buf, sfa, SFA_POSITION, 2, N3);

        /* 1) Find baryon clusters */
        #define MAX_CL 100
        ClusterInfo clusters[MAX_CL];
        int n_cl = find_clusters_bfs(phi, N, L, CLUSTER_P_THRESH, clusters, MAX_CL);

        /* Pick two largest as baryons */
        ClusterInfo b1 = {0}, b2 = {0};
        int n_baryons = 0;
        if (n_cl >= 1) { b1 = clusters[0]; n_baryons = 1; }
        if (n_cl >= 2) { b2 = clusters[1]; n_baryons = 2; }

        double bond_mx=0, bond_my=0, bond_mz=0, bond_len=0;
        if (n_baryons >= 2) {
            bond_mx = 0.5*(b1.cx + b2.cx);
            bond_my = 0.5*(b1.cy + b2.cy);
            bond_mz = 0.5*(b1.cz + b2.cz);
            double ddx = b2.cx-b1.cx, ddy = b2.cy-b1.cy, ddz = b2.cz-b1.cz;
            bond_len = sqrt(ddx*ddx + ddy*ddy + ddz*ddz);
        }

        /* 2) Compute atan2 phase at every active voxel, bin into three groups */
        /* Group A: phase near -0.79, Group B: near +2.40, Group C: near +0.70 */
        long countA=0, countB=0, countC=0, countOther=0;
        double sxA=0,syA=0,szA=0,swA=0;
        double sxB=0,syB=0,szB=0,swB=0;
        double sxC=0,syC=0,szC=0,swC=0;
        double PintC = 0;  /* integrated |P| in phase-C region */
        long countC_in_bond = 0;

        /* Phase histogram (36 bins over [-π, π]) */
        #define NHIST 36
        long phase_hist[NHIST] = {0};
        double phase_Pint_hist[NHIST] = {0};

        #pragma omp parallel for reduction(+:countA,countB,countC,countOther,sxA,syA,szA,swA,sxB,syB,szB,swB,sxC,syC,szC,swC,PintC,countC_in_bond)
        for (long idx = 0; idx < N3; idx++) {
            double p0 = phi[0][idx], p1 = phi[1][idx], p2 = phi[2][idx];
            double P = fabs(p0 * p1 * p2);
            if (P < P_THRESHOLD) continue;

            double phase = atan2(p1, p0);
            int hbin = (int)((phase + PI) / (2*PI) * NHIST);
            if (hbin < 0) hbin = 0; if (hbin >= NHIST) hbin = NHIST-1;

            /* Cannot do array reduction easily in C, skip histogram in parallel */

            double x = -L + (idx / NN) * dx;
            double y = -L + ((idx / N) % N) * dx;
            double z = -L + (idx % N) * dx;

            double dA = angle_dist(phase, PHASE_A_CENTER);
            double dB = angle_dist(phase, PHASE_B_CENTER);
            double dC = angle_dist(phase, PHASE_C_CENTER);

            if (dC < PHASE_WINDOW) {
                countC++;
                sxC += x*P; syC += y*P; szC += z*P; swC += P;
                PintC += P * dV;
                if (n_baryons >= 2 && bond_len > 0) {
                    double rdx = x-bond_mx, rdy = y-bond_my, rdz = z-bond_mz;
                    if (sqrt(rdx*rdx+rdy*rdy+rdz*rdz) < 0.5*bond_len)
                        countC_in_bond++;
                }
            } else if (dA < PHASE_WINDOW) {
                countA++;
                sxA += x*P; syA += y*P; szA += z*P; swA += P;
            } else if (dB < PHASE_WINDOW) {
                countB++;
                sxB += x*P; syB += y*P; szB += z*P; swB += P;
            } else {
                countOther++;
            }
        }

        /* Serial phase histogram pass */
        for (long idx = 0; idx < N3; idx++) {
            double P = fabs(phi[0][idx]*phi[1][idx]*phi[2][idx]);
            if (P < P_THRESHOLD) continue;
            double phase = atan2(phi[1][idx], phi[0][idx]);
            int hbin = (int)((phase + PI) / (2*PI) * NHIST);
            if (hbin < 0) hbin = 0; if (hbin >= NHIST) hbin = NHIST-1;
            phase_hist[hbin]++;
            phase_Pint_hist[hbin] += P * dV;
        }

        double cxA = swA>0 ? sxA/swA : 0, cyA = swA>0 ? syA/swA : 0, czA = swA>0 ? szA/swA : 0;
        double cxB = swB>0 ? sxB/swB : 0, cyB = swB>0 ? syB/swB : 0, czB = swB>0 ? szB/swB : 0;
        double cxC = swC>0 ? sxC/swC : 0, cyC = swC>0 ? syC/swC : 0, czC = swC>0 ? szC/swC : 0;
        double volA = countA*dV, volB = countB*dV, volC = countC*dV;

        double phaseC_bond_dist = 0;
        if (n_baryons >= 2 && swC > 0) {
            double ddx = cxC-bond_mx, ddy = cyC-bond_my, ddz = czC-bond_mz;
            phaseC_bond_dist = sqrt(ddx*ddx + ddy*ddy + ddz*ddz);
        }
        double fracC_in_bond = countC > 0 ? (double)countC_in_bond / countC : 0;

        /* TSV row */
        fprintf(tsv, "%.1f\t%d\t"
                     "%ld\t%.1f\t%.2f\t%.2f\t%.2f\t"
                     "%ld\t%.1f\t%.2f\t%.2f\t%.2f\t"
                     "%ld\t%.1f\t%.2f\t%.2f\t%.2f\t"
                     "%.4f\t"
                     "%.2f\t%.2f\t%.2f\t%.3f\t%.4f\t"
                     "%.2f\t%.2f\t%.2f\t%.3f\t%.4f\t"
                     "%.2f\t%.2f\t%.2f\t%.2f\t"
                     "%.2f\t%.4f\n",
                t, n_cl,
                countA, volA, cxA, cyA, czA,
                countB, volB, cxB, cyB, czB,
                countC, volC, cxC, cyC, czC,
                PintC,
                b1.cx, b1.cy, b1.cz, b1.phase, b1.mass,
                b2.cx, b2.cy, b2.cz, b2.phase, b2.mass,
                bond_mx, bond_my, bond_mz, bond_len,
                phaseC_bond_dist, fracC_in_bond);

        /* Detail output */
        fprintf(detail, "\n=== t=%.1f ===\n", t);
        fprintf(detail, "Clusters found: %d (top-2 masses: %.4f, %.4f)\n",
                n_cl, b1.mass, b2.mass);
        fprintf(detail, "Baryon 1: (%.2f, %.2f, %.2f) phase=%.3f mass=%.4f count=%ld\n",
                b1.cx, b1.cy, b1.cz, b1.phase, b1.mass, b1.count);
        if (n_baryons >= 2)
            fprintf(detail, "Baryon 2: (%.2f, %.2f, %.2f) phase=%.3f mass=%.4f count=%ld\n",
                    b2.cx, b2.cy, b2.cz, b2.phase, b2.mass, b2.count);
        fprintf(detail, "Bond: mid=(%.2f,%.2f,%.2f) length=%.2f\n",
                bond_mx, bond_my, bond_mz, bond_len);
        fprintf(detail, "\nPhase groups (|P|>%.3f):\n", P_THRESHOLD);
        fprintf(detail, "  A (≈-0.79): %ld voxels, vol=%.1f, centroid=(%.2f,%.2f,%.2f)\n",
                countA, volA, cxA, cyA, czA);
        fprintf(detail, "  B (≈+2.40): %ld voxels, vol=%.1f, centroid=(%.2f,%.2f,%.2f)\n",
                countB, volB, cxB, cyB, czB);
        fprintf(detail, "  C (≈+0.70): %ld voxels, vol=%.1f, centroid=(%.2f,%.2f,%.2f), P_int=%.4f\n",
                countC, volC, cxC, cyC, czC, PintC);
        fprintf(detail, "  Other:      %ld voxels\n", countOther);
        if (n_baryons >= 2) {
            fprintf(detail, "  C→bond_mid dist: %.2f (bond_len=%.2f)\n",
                    phaseC_bond_dist, bond_len);
            fprintf(detail, "  C fraction within bond region: %.4f\n", fracC_in_bond);
        }

        /* Print phase histogram */
        fprintf(detail, "\nPhase histogram (atan2, 36 bins over [-π,π]):\n");
        fprintf(detail, "  bin_center   count   P_int\n");
        for (int b = 0; b < NHIST; b++) {
            if (phase_hist[b] > 0) {
                double center = -PI + (b+0.5) * 2*PI / NHIST;
                fprintf(detail, "  %+6.2f     %7ld   %.4f\n",
                        center, phase_hist[b], phase_Pint_hist[b]);
            }
        }

        /* Print top clusters with phases */
        fprintf(detail, "\nTop clusters (by mass):\n");
        int top = n_cl < 10 ? n_cl : 10;
        for (int i = 0; i < top; i++)
            fprintf(detail, "  #%d: (%.1f,%.1f,%.1f) phase=%.3f mass=%.4f count=%ld\n",
                    i+1, clusters[i].cx, clusters[i].cy, clusters[i].cz,
                    clusters[i].phase, clusters[i].mass, clusters[i].count);

        double elapsed = omp_get_wtime() - t0;
        fprintf(stderr, " done (%.1fs) A=%ld B=%ld C=%ld other=%ld baryons=%d bond=%.1f\n",
                elapsed, countA, countB, countC, countOther, n_baryons, bond_len);

        for (int c = 0; c < 3; c++) free(phi[c]);
    }

    fclose(tsv);
    fclose(detail);
    sfa_close(sfa);
    free(buf);

    fprintf(stderr, "\nDone. Output: phase_tracking.tsv, phase_detail.txt\n");
    return 0;
}
