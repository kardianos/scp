/*  cluster_profile.c -- Per-cluster field characterization tool
 *
 *  For each selected frame of an SFA, detects clusters via BFS on |P|,
 *  then computes comprehensive per-cluster metrics in each cluster's
 *  local (centroid-centered) frame:
 *
 *    - Radial profiles:  rho(r), |P|(r), theta_rms(r), E_kin(r), E_pot(r),
 *                        |v|(r), |a|(r)   [30 bins from r=0 to R_max]
 *    - Energy decomposition per cluster
 *    - Phase coherence / winding
 *    - Velocity statistics (drift, dispersion, div(v))
 *    - Triple product structure (peak |P|, core volume)
 *    - Surface characterization (|P| at half-max, gradient at surface)
 *
 *  Output: one JSON file per input, containing an array of frame objects
 *          each with an array of cluster objects.
 *
 *  Build:
 *    gcc -O3 -fopenmp -o cluster_profile cluster_profile.c \
 *        -lzstd -lm
 *
 *  Usage:
 *    ./cluster_profile input.sfa [--frames 0,10,20,40] [--json out.json]
 *
 *  Frame indices are SFA frame numbers (0-based). Default: 0, 10, 20, 40
 *  which correspond to t=0, 50, 100, 200 for dt_save=5.
 */

#define SFA_IMPLEMENTATION
#include "../format/sfa.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>
#include <float.h>

/* ---- Physics defaults (overridden from KVMD) ---- */
static double MU = -41.345, KAPPA = 50.0, MASS2 = 2.25;
static double MTHETA2 = 0.0, ETA = 0.5;

/* ---- Tuning ---- */
#define MAX_CLUSTERS   64
#define RADIAL_BINS    30
#define P_THRESHOLD_FRAC 0.01   /* cluster threshold = frac * P_max_global */
#define MIN_CLUSTER_VOXELS 27   /* 3^3 minimum */

/* ---- Structures ---- */

typedef struct {
    double r_mid;           /* bin center radius */
    double rho;             /* mean field energy density sum(phi^2) */
    double P_abs;           /* mean |P| */
    double theta_rms;       /* rms theta in shell */
    double E_kin;           /* total kinetic energy in shell */
    double E_pot;           /* total V(P) energy in shell */
    double v_mag;           /* mean |v| in shell */
    double a_mag;           /* mean estimated |acceleration| in shell */
    int    count;           /* voxels in this bin */
} RadialBin;

typedef struct {
    int    id;              /* cluster index (1-based) */
    int    n_voxels;        /* number of voxels in cluster */

    /* Centroid (global coordinates) */
    double cx, cy, cz;

    /* Radial profiles (RADIAL_BINS entries) */
    RadialBin profile[RADIAL_BINS];
    double R_max;           /* max radius of cluster (from centroid) */

    /* Energy decomposition */
    double E_kin_phi;       /* phi kinetic energy */
    double E_kin_theta;     /* theta kinetic energy */
    double E_grad_phi;      /* phi gradient energy */
    double E_grad_theta;    /* theta gradient energy */
    double E_mass;          /* mass term energy */
    double E_pot;           /* V(P) potential energy */
    double E_coupling;      /* eta * phi . curl(theta) coupling */
    double E_total;

    /* Triple product structure */
    double P_peak;          /* max |P| in cluster */
    double P_int;           /* integral |P| dV in cluster */
    double core_volume;     /* volume where |P| > 0.5 * P_peak */

    /* Surface characterization */
    double R_half;          /* half-max radius (|P| = 0.5 * P_peak) */
    double P_at_Rhalf;      /* actual |P| at R_half (should be ~0.5*P_peak) */
    double dPdr_surface;    /* d|P|/dr at R_half (surface gradient) */

    /* Theta statistics */
    double theta_rms_total; /* overall theta rms in cluster */
    double theta_rms_core;  /* theta rms within r < R_half */

    /* Velocity field statistics */
    double v_drift[3];      /* mean velocity (bulk drift) */
    double v_disp;          /* velocity dispersion (random component) */
    double div_v;           /* mean divergence of velocity field */

    /* Phase coherence */
    double phase_corr_z;    /* phi_0 * phi_1 phase correlation along z */

    /* Moment of inertia eigenvalues (for aspect ratio) */
    double I_diag[3];       /* diagonal MOI (Ixx, Iyy, Izz) */
    double aspect;          /* I_max / I_min */
} ClusterData;

typedef struct {
    double time;
    int    frame_idx;
    int    n_clusters;
    double P_max_global;
    double threshold;
    ClusterData clusters[MAX_CLUSTERS];
} FrameResult;

/* ---- BFS cluster detection (returns label array, 0 = background) ---- */
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
        int count = 0;

        while (front < back) {
            long cur = queue[front++];
            count++;
            int i = (int)(cur / NN), j = (int)((cur / N) % N), k = (int)(cur % N);

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

        /* Remove tiny clusters */
        if (count < MIN_CLUSTER_VOXELS) {
            /* Erase this cluster label */
            for (long q = 0; q < back; q++) {
                /* Re-traverse queue to unlabel; but queue was overwritten...
                   Instead, just leave label and filter later. */
            }
            /* Actually, we already labeled them. Mark as label = -1 to remove. */
            /* Simpler: just keep them and filter by size later. */
        }
    }

    *n_clusters_out = n_clusters;
    free(queue);
    return 0;
}

/* ---- Extract one column from the raw SFA buffer as float* ---- */
static float *extract_column_f32(void *buf, SFA *sfa, int col_idx) {
    uint64_t N3 = sfa->N_total;
    uint64_t off = 0;
    for (int c = 0; c < col_idx; c++)
        off += N3 * sfa_dtype_size[sfa->columns[c].dtype];

    int dtype = sfa->columns[col_idx].dtype;
    uint8_t *src = (uint8_t*)buf + off;

    float *arr = (float*)malloc(N3 * sizeof(float));
    if (dtype == SFA_F32) {
        memcpy(arr, src, N3 * sizeof(float));
    } else if (dtype == SFA_F64) {
        double *d = (double*)src;
        for (uint64_t i = 0; i < N3; i++) arr[i] = (float)d[i];
    } else if (dtype == SFA_F16) {
        uint16_t *h = (uint16_t*)src;
        for (uint64_t i = 0; i < N3; i++) {
            int e = (h[i]>>10)&0x1F; uint16_t m = h[i]&0x3FF;
            uint16_t s = h[i]&0x8000;
            if (e==0) { arr[i]=0; } else {
                uint32_t x = ((uint32_t)s<<16)|((uint32_t)(e-15+127)<<23)|((uint32_t)m<<13);
                memcpy(&arr[i],&x,4);
            }
        }
    }
    return arr;
}

/* ---- Compute cluster centroids (weighted by |P|) ---- */
static void compute_centroids(float *phi[3], int *labels, int n_clusters,
                              int N, double L, double dx,
                              ClusterData *cd)
{
    long N3 = (long)N*N*N;
    int NN = N*N;

    /* Zero accumulators */
    double wt[MAX_CLUSTERS+1];
    double sx[MAX_CLUSTERS+1], sy[MAX_CLUSTERS+1], sz[MAX_CLUSTERS+1];
    int    cnt[MAX_CLUSTERS+1];
    memset(wt, 0, sizeof(wt));
    memset(sx, 0, sizeof(sx));
    memset(sy, 0, sizeof(sy));
    memset(sz, 0, sizeof(sz));
    memset(cnt, 0, sizeof(cnt));

    for (long idx = 0; idx < N3; idx++) {
        int lab = labels[idx];
        if (lab <= 0 || lab > n_clusters) continue;
        int i = (int)(idx / NN), j = (int)((idx / N) % N), k = (int)(idx % N);
        double x = -L + i*dx, y = -L + j*dx, z = -L + k*dx;
        double P = fabs((double)phi[0][idx] * phi[1][idx] * phi[2][idx]);
        sx[lab] += x*P; sy[lab] += y*P; sz[lab] += z*P;
        wt[lab] += P;
        cnt[lab]++;
    }

    for (int c = 1; c <= n_clusters; c++) {
        cd[c-1].id = c;
        cd[c-1].n_voxels = cnt[c];
        if (wt[c] > 1e-30) {
            cd[c-1].cx = sx[c] / wt[c];
            cd[c-1].cy = sy[c] / wt[c];
            cd[c-1].cz = sz[c] / wt[c];
        }
    }
}

/* ---- Process one frame ---- */
static void process_frame(SFA *sfa, void *buf, int N, double L,
                          FrameResult *res)
{
    long N3 = (long)N*N*N;
    int NN = N*N;
    double dx = 2.0*L / (N-1);
    double dV = dx*dx*dx;
    double idx1 = 1.0 / (2.0*dx);

    /* Extract all 12 columns as f32 */
    float *phi[3], *theta[3], *phi_vel[3], *theta_vel[3];
    for (int a = 0; a < 3; a++) {
        phi[a]       = extract_column_f32(buf, sfa, a);       /* cols 0-2 */
        theta[a]     = extract_column_f32(buf, sfa, a+3);     /* cols 3-5 */
        phi_vel[a]   = extract_column_f32(buf, sfa, a+6);     /* cols 6-8 */
        theta_vel[a] = extract_column_f32(buf, sfa, a+9);     /* cols 9-11 */
    }

    /* Global P_max */
    double P_max_global = 0;
    #pragma omp parallel for reduction(max:P_max_global) schedule(static)
    for (long idx = 0; idx < N3; idx++) {
        double P = fabs((double)phi[0][idx] * phi[1][idx] * phi[2][idx]);
        if (P > P_max_global) P_max_global = P;
    }
    res->P_max_global = P_max_global;
    double threshold = P_THRESHOLD_FRAC * P_max_global;
    if (threshold < 1e-20) threshold = 1e-20;
    res->threshold = threshold;

    /* BFS cluster detection */
    int *labels = (int*)calloc(N3, sizeof(int));
    int n_clusters = 0;
    bfs_clusters(phi[0], phi[1], phi[2], N, threshold, labels, &n_clusters);

    /* Compute cluster sizes and filter tiny ones */
    int *sizes = (int*)calloc(n_clusters + 1, sizeof(int));
    for (long idx = 0; idx < N3; idx++) {
        if (labels[idx] > 0 && labels[idx] <= n_clusters)
            sizes[labels[idx]]++;
    }
    /* Remap: keep only clusters with >= MIN_CLUSTER_VOXELS */
    int *remap = (int*)calloc(n_clusters + 1, sizeof(int));
    int n_valid = 0;
    for (int c = 1; c <= n_clusters; c++) {
        if (sizes[c] >= MIN_CLUSTER_VOXELS) {
            n_valid++;
            remap[c] = n_valid;
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
    if (n_clusters > MAX_CLUSTERS) n_clusters = MAX_CLUSTERS;
    res->n_clusters = n_clusters;
    free(sizes); free(remap);

    if (n_clusters == 0) {
        free(labels);
        for (int a = 0; a < 3; a++) {
            free(phi[a]); free(theta[a]); free(phi_vel[a]); free(theta_vel[a]);
        }
        return;
    }

    /* Compute centroids */
    compute_centroids(phi, labels, n_clusters, N, L, dx, res->clusters);

    /* ---- Per-cluster analysis ---- */
    for (int c = 0; c < n_clusters; c++) {
        ClusterData *cd = &res->clusters[c];

        /* Find R_max for this cluster */
        double Rmax = 0;
        for (long idx = 0; idx < N3; idx++) {
            if (labels[idx] != c+1) continue;
            int i = (int)(idx / NN), j = (int)((idx / N) % N), k = (int)(idx % N);
            double rx = -L + i*dx - cd->cx;
            double ry = -L + j*dx - cd->cy;
            double rz = -L + k*dx - cd->cz;
            double r = sqrt(rx*rx + ry*ry + rz*rz);
            if (r > Rmax) Rmax = r;
        }
        cd->R_max = Rmax;
        if (Rmax < dx) Rmax = dx;
        double dr = Rmax / RADIAL_BINS;

        /* Initialize profile bins */
        for (int b = 0; b < RADIAL_BINS; b++) {
            memset(&cd->profile[b], 0, sizeof(RadialBin));
            cd->profile[b].r_mid = (b + 0.5) * dr;
        }

        /* Accumulators */
        double ek_phi = 0, ek_theta = 0;
        double eg_phi = 0, eg_theta = 0;
        double e_mass = 0, e_pot = 0, e_coupling = 0;
        double P_peak = 0, P_int = 0, core_vol = 0;
        double vdrift[3] = {0,0,0}, v2sum = 0;
        double divv_sum = 0;
        double theta2_total = 0, theta2_core = 0;
        int theta_core_cnt = 0;
        double Ixx=0, Iyy=0, Izz=0;
        double phase_zsum = 0; int phase_zcnt = 0;

        /* First pass: find P_peak for core volume calculation */
        for (long idx = 0; idx < N3; idx++) {
            if (labels[idx] != c+1) continue;
            double P = fabs((double)phi[0][idx] * phi[1][idx] * phi[2][idx]);
            if (P > P_peak) P_peak = P;
        }
        cd->P_peak = P_peak;
        double P_half = 0.5 * P_peak;

        /* Second pass: compute everything */
        for (long idx = 0; idx < N3; idx++) {
            if (labels[idx] != c+1) continue;

            int i = (int)(idx / NN), j = (int)((idx / N) % N), k = (int)(idx % N);
            double x = -L + i*dx, y = -L + j*dx, z = -L + k*dx;
            double rx = x - cd->cx, ry = y - cd->cy, rz = z - cd->cz;
            double r = sqrt(rx*rx + ry*ry + rz*rz);

            /* Radial bin index */
            int bin = (int)(r / dr);
            if (bin >= RADIAL_BINS) bin = RADIAL_BINS - 1;

            /* Field values */
            double p0 = phi[0][idx], p1 = phi[1][idx], p2 = phi[2][idx];
            double P = p0 * p1 * p2;
            double Pa = fabs(P);
            double rho = p0*p0 + p1*p1 + p2*p2;

            /* Neighbor indices (periodic) */
            int ip = (i+1)%N, im = (i-1+N)%N;
            int jp = (j+1)%N, jm = (j-1+N)%N;
            int kp = (k+1)%N, km = (k-1+N)%N;
            long nip = (long)ip*NN+j*N+k, nim = (long)im*NN+j*N+k;
            long njp = (long)i*NN+jp*N+k, njm = (long)i*NN+jm*N+k;
            long nkp = (long)i*NN+j*N+kp, nkm = (long)i*NN+j*N+km;

            /* Gradient energy per component */
            double eg_this = 0;
            for (int a = 0; a < 3; a++) {
                double gx = (phi[a][nip] - phi[a][nim]) * idx1;
                double gy = (phi[a][njp] - phi[a][njm]) * idx1;
                double gz = (phi[a][nkp] - phi[a][nkm]) * idx1;
                eg_this += 0.5 * (gx*gx + gy*gy + gz*gz);
            }
            eg_phi += eg_this * dV;

            /* Mass energy */
            double em_this = 0.5 * MASS2 * rho;
            e_mass += em_this * dV;

            /* Potential V(P) */
            double P2 = P*P;
            double ep_this = (MU/2.0) * P2 / (1.0 + KAPPA*P2);
            e_pot += ep_this * dV;

            /* Kinetic energy */
            double ek_phi_this = 0, ek_theta_this = 0;
            for (int a = 0; a < 3; a++) {
                double pv = phi_vel[a][idx];
                double tv = theta_vel[a][idx];
                ek_phi_this += 0.5 * pv*pv;
                ek_theta_this += 0.5 * tv*tv;
            }
            ek_phi += ek_phi_this * dV;
            ek_theta += ek_theta_this * dV;

            /* Theta gradient energy */
            double etg_this = 0;
            for (int a = 0; a < 3; a++) {
                double gx = (theta[a][nip] - theta[a][nim]) * idx1;
                double gy = (theta[a][njp] - theta[a][njm]) * idx1;
                double gz = (theta[a][nkp] - theta[a][nkm]) * idx1;
                etg_this += 0.5 * (gx*gx + gy*gy + gz*gz);
            }
            eg_theta += etg_this * dV;

            /* Coupling energy: -eta * phi . curl(theta) */
            double ct0 = (theta[2][njp]-theta[2][njm] - theta[1][nkp]+theta[1][nkm])*idx1;
            double ct1 = (theta[0][nkp]-theta[0][nkm] - theta[2][nip]+theta[2][nim])*idx1;
            double ct2 = (theta[1][nip]-theta[1][nim] - theta[0][njp]+theta[0][njm])*idx1;
            double ec_this = -ETA * (p0*ct0 + p1*ct1 + p2*ct2);
            e_coupling += ec_this * dV;

            /* Triple product integral */
            P_int += Pa * dV;

            /* Core volume: where |P| > 0.5 * P_peak */
            if (Pa > P_half) core_vol += dV;

            /* Velocity statistics */
            double vx = phi_vel[0][idx], vy = phi_vel[1][idx], vz = phi_vel[2][idx];
            double vmag = sqrt(vx*vx + vy*vy + vz*vz);
            vdrift[0] += vx; vdrift[1] += vy; vdrift[2] += vz;
            v2sum += vmag*vmag;

            /* Divergence of velocity: dv_x/dx + dv_y/dy + dv_z/dz */
            double dvx = (phi_vel[0][nip] - phi_vel[0][nim]) * idx1;
            double dvy = (phi_vel[1][njp] - phi_vel[1][njm]) * idx1;
            double dvz = (phi_vel[2][nkp] - phi_vel[2][nkm]) * idx1;
            divv_sum += (dvx + dvy + dvz);

            /* Theta rms */
            double t2 = 0;
            for (int a = 0; a < 3; a++) {
                double tv = theta[a][idx];
                t2 += tv*tv;
            }
            theta2_total += t2;

            /* Moment of inertia (diagonal approx) */
            double w = rho * dV;
            Ixx += (ry*ry + rz*rz)*w;
            Iyy += (rx*rx + rz*rz)*w;
            Izz += (rx*rx + ry*ry)*w;

            /* Approximate acceleration: compute force from the equation of motion.
             * F_a = laplacian(phi_a) - m^2 * phi_a - dV/dphi_a + eta * curl(theta)_a
             * Instead of computing the full Laplacian (expensive with 7-point stencil
             * already available), use the simpler estimate: |a| ~ |F|/1 since mass
             * density = 1 in code units.  Actually just compute the Laplacian. */
            double a_mag2 = 0;
            for (int a = 0; a < 3; a++) {
                /* 7-point Laplacian (2nd order) */
                double lap = (phi[a][nip] + phi[a][nim] + phi[a][njp] + phi[a][njm]
                             + phi[a][nkp] + phi[a][nkm] - 6.0*phi[a][idx]) / (dx*dx);

                /* dV/dphi_a: V(P) = (mu/2)*P^2/(1+kappa*P^2)
                 * dV/dphi_a = mu * P * dP/dphi_a / (1+kappa*P^2)^2
                 * dP/dphi_a = product of other two fields */
                double dPda;
                if (a == 0) dPda = p1*p2;
                else if (a == 1) dPda = p0*p2;
                else dPda = p0*p1;
                double denom = (1.0 + KAPPA*P2);
                double dVda = MU * P * dPda / (denom*denom);

                /* curl(theta)_a (already computed for a=0,1,2 as ct0,ct1,ct2) */
                double curl_a;
                if (a == 0) curl_a = ct0;
                else if (a == 1) curl_a = ct1;
                else curl_a = ct2;

                double F_a = lap - MASS2 * phi[a][idx] - dVda + ETA * curl_a;
                a_mag2 += F_a * F_a;
            }

            /* Fill radial profile bins */
            RadialBin *rb = &cd->profile[bin];
            rb->rho += rho;
            rb->P_abs += Pa;
            rb->theta_rms += t2;
            rb->E_kin += (ek_phi_this + ek_theta_this) * dV;
            rb->E_pot += ep_this * dV;
            rb->v_mag += vmag;
            rb->a_mag += sqrt(a_mag2);
            rb->count++;

            /* Core theta rms */
            if (r < cd->R_max * 0.5) {
                theta2_core += t2;
                theta_core_cnt++;
            }

            /* Phase coherence along z (for voxels near z-axis in cluster frame) */
            if (fabs(rx) < 2*dx && fabs(ry) < 2*dx) {
                /* Sign of p0*p1 encodes relative phase */
                phase_zsum += p0 * p1;
                phase_zcnt++;
            }
        }

        /* Store energies */
        cd->E_kin_phi = ek_phi;
        cd->E_kin_theta = ek_theta;
        cd->E_grad_phi = eg_phi;
        cd->E_grad_theta = eg_theta;
        cd->E_mass = e_mass;
        cd->E_pot = e_pot;
        cd->E_coupling = e_coupling;
        cd->E_total = ek_phi + ek_theta + eg_phi + eg_theta + e_mass + e_pot + e_coupling;

        cd->P_int = P_int;
        cd->core_volume = core_vol;

        /* Theta rms */
        int nv = cd->n_voxels;
        if (nv > 0) {
            cd->theta_rms_total = sqrt(theta2_total / (3.0 * nv));
            cd->theta_rms_core = (theta_core_cnt > 0) ?
                                 sqrt(theta2_core / (3.0 * theta_core_cnt)) : 0;
        }

        /* Velocity stats */
        if (nv > 0) {
            cd->v_drift[0] = vdrift[0] / nv;
            cd->v_drift[1] = vdrift[1] / nv;
            cd->v_drift[2] = vdrift[2] / nv;
            double vd2 = cd->v_drift[0]*cd->v_drift[0] + cd->v_drift[1]*cd->v_drift[1]
                       + cd->v_drift[2]*cd->v_drift[2];
            cd->v_disp = sqrt(fmax(0, v2sum / nv - vd2));
            cd->div_v = divv_sum / nv;
        }

        /* Moment of inertia / aspect ratio */
        cd->I_diag[0] = Ixx;
        cd->I_diag[1] = Iyy;
        cd->I_diag[2] = Izz;
        double Imax = Ixx, Imin = Ixx;
        if (Iyy > Imax) Imax = Iyy; if (Iyy < Imin) Imin = Iyy;
        if (Izz > Imax) Imax = Izz; if (Izz < Imin) Imin = Izz;
        cd->aspect = (Imin > 1e-30) ? Imax / Imin : 999;

        /* Phase coherence */
        cd->phase_corr_z = (phase_zcnt > 0) ? phase_zsum / phase_zcnt : 0;

        /* Normalize radial profile bins to means */
        for (int b = 0; b < RADIAL_BINS; b++) {
            RadialBin *rb = &cd->profile[b];
            if (rb->count > 0) {
                rb->rho /= rb->count;
                rb->P_abs /= rb->count;
                rb->theta_rms = sqrt(rb->theta_rms / (3.0 * rb->count));
                /* E_kin and E_pot remain as totals in this shell */
                rb->v_mag /= rb->count;
                rb->a_mag /= rb->count;
            }
        }

        /* Surface characterization: find R_half */
        cd->R_half = cd->R_max;
        cd->P_at_Rhalf = 0;
        cd->dPdr_surface = 0;
        for (int b = 1; b < RADIAL_BINS; b++) {
            if (cd->profile[b].P_abs < P_half && cd->profile[b-1].P_abs >= P_half) {
                /* Interpolate R_half */
                double r0 = cd->profile[b-1].r_mid;
                double r1 = cd->profile[b].r_mid;
                double p0v = cd->profile[b-1].P_abs;
                double p1v = cd->profile[b].P_abs;
                if (fabs(p0v - p1v) > 1e-30) {
                    double frac = (p0v - P_half) / (p0v - p1v);
                    cd->R_half = r0 + frac * (r1 - r0);
                } else {
                    cd->R_half = 0.5*(r0+r1);
                }
                cd->P_at_Rhalf = P_half;
                cd->dPdr_surface = (p1v - p0v) / (r1 - r0);
                break;
            }
        }
    }

    free(labels);
    for (int a = 0; a < 3; a++) {
        free(phi[a]); free(theta[a]); free(phi_vel[a]); free(theta_vel[a]);
    }
}

/* ---- JSON output ---- */
static void write_json(FILE *fp, FrameResult *results, int n_frames, const char *sfa_path,
                       int N, double L)
{
    fprintf(fp, "{\n");
    fprintf(fp, "  \"sfa\": \"%s\",\n", sfa_path);
    fprintf(fp, "  \"N\": %d,\n", N);
    fprintf(fp, "  \"L\": %.2f,\n", L);
    fprintf(fp, "  \"physics\": {\"mu\": %.6f, \"kappa\": %.1f, \"m2\": %.4f, \"eta\": %.3f},\n",
            MU, KAPPA, MASS2, ETA);
    fprintf(fp, "  \"radial_bins\": %d,\n", RADIAL_BINS);
    fprintf(fp, "  \"frames\": [\n");

    for (int f = 0; f < n_frames; f++) {
        FrameResult *fr = &results[f];
        fprintf(fp, "    {\n");
        fprintf(fp, "      \"time\": %.2f,\n", fr->time);
        fprintf(fp, "      \"frame_idx\": %d,\n", fr->frame_idx);
        fprintf(fp, "      \"P_max_global\": %.8e,\n", fr->P_max_global);
        fprintf(fp, "      \"threshold\": %.8e,\n", fr->threshold);
        fprintf(fp, "      \"n_clusters\": %d,\n", fr->n_clusters);
        fprintf(fp, "      \"clusters\": [\n");

        for (int c = 0; c < fr->n_clusters; c++) {
            ClusterData *cd = &fr->clusters[c];
            fprintf(fp, "        {\n");
            fprintf(fp, "          \"id\": %d,\n", cd->id);
            fprintf(fp, "          \"n_voxels\": %d,\n", cd->n_voxels);
            fprintf(fp, "          \"centroid\": [%.4f, %.4f, %.4f],\n", cd->cx, cd->cy, cd->cz);
            fprintf(fp, "          \"R_max\": %.4f,\n", cd->R_max);

            /* Energy decomposition */
            fprintf(fp, "          \"energy\": {\n");
            fprintf(fp, "            \"E_kin_phi\": %.6e,\n", cd->E_kin_phi);
            fprintf(fp, "            \"E_kin_theta\": %.6e,\n", cd->E_kin_theta);
            fprintf(fp, "            \"E_grad_phi\": %.6e,\n", cd->E_grad_phi);
            fprintf(fp, "            \"E_grad_theta\": %.6e,\n", cd->E_grad_theta);
            fprintf(fp, "            \"E_mass\": %.6e,\n", cd->E_mass);
            fprintf(fp, "            \"E_pot\": %.6e,\n", cd->E_pot);
            fprintf(fp, "            \"E_coupling\": %.6e,\n", cd->E_coupling);
            fprintf(fp, "            \"E_total\": %.6e\n", cd->E_total);
            fprintf(fp, "          },\n");

            /* Triple product */
            fprintf(fp, "          \"triple_product\": {\n");
            fprintf(fp, "            \"P_peak\": %.8e,\n", cd->P_peak);
            fprintf(fp, "            \"P_int\": %.6e,\n", cd->P_int);
            fprintf(fp, "            \"core_volume\": %.4f\n", cd->core_volume);
            fprintf(fp, "          },\n");

            /* Surface */
            fprintf(fp, "          \"surface\": {\n");
            fprintf(fp, "            \"R_half\": %.4f,\n", cd->R_half);
            fprintf(fp, "            \"P_at_Rhalf\": %.8e,\n", cd->P_at_Rhalf);
            fprintf(fp, "            \"dPdr_surface\": %.8e\n", cd->dPdr_surface);
            fprintf(fp, "          },\n");

            /* Theta */
            fprintf(fp, "          \"theta\": {\n");
            fprintf(fp, "            \"rms_total\": %.8e,\n", cd->theta_rms_total);
            fprintf(fp, "            \"rms_core\": %.8e\n", cd->theta_rms_core);
            fprintf(fp, "          },\n");

            /* Velocity */
            fprintf(fp, "          \"velocity\": {\n");
            fprintf(fp, "            \"drift\": [%.6e, %.6e, %.6e],\n",
                    cd->v_drift[0], cd->v_drift[1], cd->v_drift[2]);
            fprintf(fp, "            \"dispersion\": %.6e,\n", cd->v_disp);
            fprintf(fp, "            \"div_v\": %.6e\n", cd->div_v);
            fprintf(fp, "          },\n");

            /* Phase coherence */
            fprintf(fp, "          \"phase_corr_z\": %.6e,\n", cd->phase_corr_z);

            /* Moment of inertia / aspect */
            fprintf(fp, "          \"moment_of_inertia\": [%.6e, %.6e, %.6e],\n",
                    cd->I_diag[0], cd->I_diag[1], cd->I_diag[2]);
            fprintf(fp, "          \"aspect\": %.4f,\n", cd->aspect);

            /* Radial profiles */
            fprintf(fp, "          \"radial_profile\": [\n");
            for (int b = 0; b < RADIAL_BINS; b++) {
                RadialBin *rb = &cd->profile[b];
                fprintf(fp, "            {\"r\": %.4f, \"rho\": %.6e, \"P_abs\": %.6e, "
                        "\"theta_rms\": %.6e, \"E_kin\": %.6e, \"E_pot\": %.6e, "
                        "\"v_mag\": %.6e, \"a_mag\": %.6e, \"count\": %d}%s\n",
                        rb->r_mid, rb->rho, rb->P_abs, rb->theta_rms,
                        rb->E_kin, rb->E_pot, rb->v_mag, rb->a_mag, rb->count,
                        (b < RADIAL_BINS-1) ? "," : "");
            }
            fprintf(fp, "          ]\n");

            fprintf(fp, "        }%s\n", (c < fr->n_clusters-1) ? "," : "");
        }

        fprintf(fp, "      ]\n");
        fprintf(fp, "    }%s\n", (f < n_frames-1) ? "," : "");
    }

    fprintf(fp, "  ]\n");
    fprintf(fp, "}\n");
}

/* ---- main ---- */
int main(int argc, char **argv)
{
    if (argc < 2) {
        fprintf(stderr, "Usage: %s input.sfa [--frames 0,10,20,40] [--json out.json]\n", argv[0]);
        return 1;
    }

    const char *sfa_path = argv[1];
    const char *json_path = NULL;
    int frame_list[256];
    int n_frames = 0;

    /* Parse args */
    for (int i = 2; i < argc; i++) {
        if (!strcmp(argv[i], "--json") && i+1 < argc) {
            json_path = argv[++i];
        } else if (!strcmp(argv[i], "--frames") && i+1 < argc) {
            i++;
            /* Parse comma-separated frame list */
            char *s = argv[i];
            while (*s && n_frames < 256) {
                frame_list[n_frames++] = atoi(s);
                while (*s && *s != ',') s++;
                if (*s == ',') s++;
            }
        }
    }

    /* Default frames: 0, 10, 20, 40 (t=0, 50, 100, 200) */
    if (n_frames == 0) {
        frame_list[0] = 0;
        frame_list[1] = 10;
        frame_list[2] = 20;
        frame_list[3] = 40;
        n_frames = 4;
    }

    /* Open SFA */
    SFA *sfa = sfa_open(sfa_path);
    if (!sfa) { fprintf(stderr, "Cannot open %s\n", sfa_path); return 1; }

    int N = sfa->Nx;
    double L = sfa->Lx;

    /* Read KVMD for physics parameters */
    SFA_KVMDSet kv[16];
    int n_kv = sfa_read_kvmd(sfa, kv, 16);
    for (int s = 0; s < n_kv; s++) {
        for (int p = 0; p < kv[s].n_pairs; p++) {
            if (!strcmp(kv[s].keys[p], "mu"))    MU = atof(kv[s].values[p]);
            if (!strcmp(kv[s].keys[p], "kappa")) KAPPA = atof(kv[s].values[p]);
            if (!strcmp(kv[s].keys[p], "m"))   { double m = atof(kv[s].values[p]); MASS2 = m*m; }
            if (!strcmp(kv[s].keys[p], "eta"))   ETA = atof(kv[s].values[p]);
            if (!strcmp(kv[s].keys[p], "m_theta")){ double mt = atof(kv[s].values[p]); MTHETA2 = mt*mt; }
        }
    }

    fprintf(stderr, "cluster_profile: %s (N=%d L=%.1f %u frames)\n", sfa_path, N, L, sfa->total_frames);
    fprintf(stderr, "Physics: mu=%.3f kappa=%.1f m^2=%.4f eta=%.3f\n", MU, KAPPA, MASS2, ETA);

    /* Clamp frame indices to available range */
    for (int i = 0; i < n_frames; i++) {
        if (frame_list[i] >= (int)sfa->total_frames)
            frame_list[i] = sfa->total_frames - 1;
        if (frame_list[i] < 0) frame_list[i] = 0;
    }

    /* Allocate frame buffer */
    void *buf = malloc(sfa->frame_bytes);
    if (!buf) { fprintf(stderr, "Cannot allocate %lu bytes for frame buffer\n",
                        (unsigned long)sfa->frame_bytes); return 1; }

    if (n_frames <= 0) { fprintf(stderr, "No frames to process\n"); return 1; }
    FrameResult *results = (FrameResult*)calloc((unsigned)n_frames, sizeof(FrameResult));

    for (int fi = 0; fi < n_frames; fi++) {
        int f = frame_list[fi];
        fprintf(stderr, "Processing frame %d ...\n", f);

        sfa_read_frame(sfa, f, buf);
        double t = sfa_frame_time(sfa, f);

        results[fi].time = t;
        results[fi].frame_idx = f;

        process_frame(sfa, buf, N, L, &results[fi]);

        fprintf(stderr, "  t=%.1f: %d clusters, P_max=%.4e\n",
                t, results[fi].n_clusters, results[fi].P_max_global);
        for (int c = 0; c < results[fi].n_clusters; c++) {
            ClusterData *cd = &results[fi].clusters[c];
            fprintf(stderr, "    Cluster %d: %d voxels, centroid=(%.1f,%.1f,%.1f), "
                    "P_peak=%.3e, E_pot=%.1f, R_half=%.2f, theta_rms=%.4e\n",
                    cd->id, cd->n_voxels, cd->cx, cd->cy, cd->cz,
                    cd->P_peak, cd->E_pot, cd->R_half, cd->theta_rms_total);
        }
    }

    /* Write JSON */
    FILE *jfp = stdout;
    if (json_path) {
        jfp = fopen(json_path, "w");
        if (!jfp) { fprintf(stderr, "Cannot open %s for writing\n", json_path); jfp = stdout; }
    }
    write_json(jfp, results, n_frames, sfa_path, N, L);
    if (json_path && jfp != stdout) {
        fclose(jfp);
        fprintf(stderr, "JSON written to %s\n", json_path);
    }

    free(buf);
    free(results);
    sfa_close(sfa);
    return 0;
}
