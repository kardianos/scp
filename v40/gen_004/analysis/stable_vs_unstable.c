/*  stable_vs_unstable.c -- Compare early vs late clusters to identify
 *  stability signatures.
 *
 *  Takes ONE SFA file and two frame indices (early, late).
 *  1. Detects clusters at both times
 *  2. Matches early->late clusters by centroid proximity
 *  3. Classifies matched pairs as GREW (stable) or SHRANK (unstable)
 *  4. Computes statistical comparison of per-cluster metrics
 *
 *  Build:
 *    gcc -O3 -fopenmp -o stable_vs_unstable stable_vs_unstable.c \
 *        -I/home/d/code/scp -lzstd -lm
 *
 *  Usage:
 *    ./stable_vs_unstable input.sfa --early 0 --late 40 [--json out.json]
 */

#define SFA_IMPLEMENTATION
#include "sfa/format/sfa.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>
#include <float.h>

/* ---- Physics defaults ---- */
static double MU = -41.345, KAPPA = 50.0, MASS2 = 2.25;
static double MTHETA2 = 0.0, ETA = 0.5;

#define MAX_CLUSTERS   64
#define RADIAL_BINS    30
#define P_THRESHOLD_FRAC 0.01
#define MIN_CLUSTER_VOXELS 27
#define MATCH_RADIUS   10.0   /* max centroid distance for a match */

/* ---- Compact cluster info ---- */
typedef struct {
    int    id;
    int    n_voxels;
    double cx, cy, cz;       /* centroid */
    double R_max;
    double R_half;            /* half-max radius */

    /* Energy */
    double E_kin_phi, E_kin_theta;
    double E_grad_phi, E_grad_theta;
    double E_mass, E_pot, E_coupling, E_total;

    /* Triple product */
    double P_peak, P_int, core_volume;

    /* Theta */
    double theta_rms, theta_rms_core;

    /* Velocity */
    double v_drift_mag, v_disp, div_v;

    /* Structure */
    double aspect;
    double rho_peak;          /* peak field energy density */
    double dPdr_surface;

    /* Radial profiles (means) */
    double prof_rho[RADIAL_BINS];
    double prof_P[RADIAL_BINS];
    double prof_theta[RADIAL_BINS];
    double prof_Ekin[RADIAL_BINS];
    double prof_Epot[RADIAL_BINS];
    double prof_vmag[RADIAL_BINS];
    double prof_amag[RADIAL_BINS];
    int    prof_count[RADIAL_BINS];
} ClusterInfo;

typedef struct {
    double time;
    int    frame_idx;
    int    n_clusters;
    double P_max_global;
    ClusterInfo clusters[MAX_CLUSTERS];
} FrameInfo;

/* ---- Matched pair ---- */
typedef struct {
    int early_id, late_id;    /* cluster indices (0-based in array) */
    double distance;          /* centroid distance */
    double dP_peak;           /* P_peak(late) - P_peak(early) */
    double dP_int;            /* P_int(late) - P_int(early) */
    double dn_voxels;         /* n_voxels change */
    int    grew;              /* 1=grew (stable), 0=shrank (unstable) */
} MatchedPair;

/* ---- Extract f32 column ---- */
static float *extract_col(void *buf, SFA *sfa, int col_idx) {
    uint64_t N3 = sfa->N_total;
    uint64_t off = 0;
    for (int c = 0; c < col_idx; c++)
        off += N3 * sfa_dtype_size[sfa->columns[c].dtype];
    int dtype = sfa->columns[col_idx].dtype;
    uint8_t *src = (uint8_t*)buf + off;
    float *arr = (float*)malloc(N3 * sizeof(float));
    if (dtype == SFA_F32) memcpy(arr, src, N3*sizeof(float));
    else if (dtype == SFA_F64) {
        double *d = (double*)src;
        for (uint64_t i = 0; i < N3; i++) arr[i] = (float)d[i];
    } else if (dtype == SFA_F16) {
        uint16_t *h = (uint16_t*)src;
        for (uint64_t i = 0; i < N3; i++) {
            int e = (h[i]>>10)&0x1F; uint16_t m = h[i]&0x3FF;
            uint16_t s = h[i]&0x8000;
            if (e==0) arr[i]=0; else {
                uint32_t x = ((uint32_t)s<<16)|((uint32_t)(e-15+127)<<23)|((uint32_t)m<<13);
                memcpy(&arr[i],&x,4);
            }
        }
    }
    return arr;
}

/* ---- BFS ---- */
static int bfs_clusters(float *phi0, float *phi1, float *phi2,
                        int N, double threshold, int *labels)
{
    long N3 = (long)N*N*N;
    int NN = N*N;
    memset(labels, 0, N3 * sizeof(int));
    long *queue = (long*)malloc(N3 * sizeof(long));
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
            int i = (int)(cur/NN), j = (int)((cur/N)%N), k = (int)(cur%N);
            int di[]={1,-1,0,0,0,0}, dj[]={0,0,1,-1,0,0}, dk[]={0,0,0,0,1,-1};
            for (int d = 0; d < 6; d++) {
                int ni = (i+di[d]+N)%N, nj = (j+dj[d]+N)%N, nk = (k+dk[d]+N)%N;
                long nidx = (long)ni*NN + nj*N + nk;
                if (labels[nidx] == 0) {
                    double nP = fabs((double)phi0[nidx]*phi1[nidx]*phi2[nidx]);
                    if (nP >= threshold) { labels[nidx] = n_clusters; queue[back++] = nidx; }
                }
            }
        }
    }
    free(queue);

    /* Filter small clusters */
    int *sizes = (int*)calloc(n_clusters+1, sizeof(int));
    for (long idx = 0; idx < N3; idx++)
        if (labels[idx] > 0) sizes[labels[idx]]++;
    int *remap = (int*)calloc(n_clusters+1, sizeof(int));
    int n_valid = 0;
    for (int c = 1; c <= n_clusters; c++)
        if (sizes[c] >= MIN_CLUSTER_VOXELS) remap[c] = ++n_valid;
    for (long idx = 0; idx < N3; idx++) {
        int lab = labels[idx];
        labels[idx] = (lab > 0 && lab <= n_clusters) ? remap[lab] : 0;
    }
    free(sizes); free(remap);
    return n_valid;
}

/* ---- Analyze one frame into FrameInfo ---- */
static void analyze_frame(SFA *sfa, void *buf, int N, double L, FrameInfo *fi)
{
    long N3 = (long)N*N*N;
    int NN = N*N;
    double dx = 2.0*L/(N-1), dV = dx*dx*dx, idx1 = 1.0/(2.0*dx);

    float *phi[3], *theta[3], *phi_vel[3], *theta_vel[3];
    for (int a = 0; a < 3; a++) {
        phi[a]       = extract_col(buf, sfa, a);
        theta[a]     = extract_col(buf, sfa, a+3);
        phi_vel[a]   = extract_col(buf, sfa, a+6);
        theta_vel[a] = extract_col(buf, sfa, a+9);
    }

    /* Global P_max */
    double P_max = 0;
    for (long idx = 0; idx < N3; idx++) {
        double P = fabs((double)phi[0][idx]*phi[1][idx]*phi[2][idx]);
        if (P > P_max) P_max = P;
    }
    fi->P_max_global = P_max;
    double threshold = P_THRESHOLD_FRAC * P_max;
    if (threshold < 1e-20) threshold = 1e-20;

    int *labels = (int*)calloc(N3, sizeof(int));
    int nc = bfs_clusters(phi[0], phi[1], phi[2], N, threshold, labels);
    fi->n_clusters = nc;

    /* Per-cluster analysis */
    for (int c = 0; c < nc; c++) {
        ClusterInfo *ci = &fi->clusters[c];
        memset(ci, 0, sizeof(ClusterInfo));
        ci->id = c + 1;

        /* Centroid */
        double sx=0, sy=0, sz=0, wt=0;
        int cnt = 0;
        for (long idx = 0; idx < N3; idx++) {
            if (labels[idx] != c+1) continue;
            cnt++;
            int i = (int)(idx/NN), j = (int)((idx/N)%N), k = (int)(idx%N);
            double x = -L+i*dx, y = -L+j*dx, z = -L+k*dx;
            double P = fabs((double)phi[0][idx]*phi[1][idx]*phi[2][idx]);
            sx += x*P; sy += y*P; sz += z*P; wt += P;
        }
        ci->n_voxels = cnt;
        if (wt > 1e-30) { ci->cx = sx/wt; ci->cy = sy/wt; ci->cz = sz/wt; }

        /* Full metrics */
        double Rmax = 0, P_peak = 0;
        double ek_phi=0, ek_theta=0, eg_phi=0, eg_theta=0;
        double e_mass=0, e_pot=0, e_coupling=0;
        double P_int=0, core_vol=0, rho_peak=0;
        double theta2=0, theta2_core=0; int t_core_cnt=0;
        double vdrift[3]={0}, v2sum=0, divv_sum=0;
        double Ixx=0, Iyy=0, Izz=0;

        /* First pass: find P_peak and Rmax */
        for (long idx = 0; idx < N3; idx++) {
            if (labels[idx] != c+1) continue;
            double P = fabs((double)phi[0][idx]*phi[1][idx]*phi[2][idx]);
            if (P > P_peak) P_peak = P;
            int i = (int)(idx/NN), j = (int)((idx/N)%N), k = (int)(idx%N);
            double rx = -L+i*dx-ci->cx, ry = -L+j*dx-ci->cy, rz = -L+k*dx-ci->cz;
            double r = sqrt(rx*rx+ry*ry+rz*rz);
            if (r > Rmax) Rmax = r;
        }
        ci->R_max = Rmax;
        ci->P_peak = P_peak;
        double P_half = 0.5 * P_peak;
        double dr = (Rmax > dx) ? Rmax / RADIAL_BINS : dx;

        /* Second pass */
        for (long idx = 0; idx < N3; idx++) {
            if (labels[idx] != c+1) continue;
            int i = (int)(idx/NN), j = (int)((idx/N)%N), k = (int)(idx%N);
            double x = -L+i*dx, y = -L+j*dx, z = -L+k*dx;
            double rx = x-ci->cx, ry = y-ci->cy, rz = z-ci->cz;
            double r = sqrt(rx*rx+ry*ry+rz*rz);

            double p0=phi[0][idx], p1=phi[1][idx], p2=phi[2][idx];
            double P=p0*p1*p2, Pa=fabs(P), P2=P*P;
            double rho = p0*p0+p1*p1+p2*p2;
            if (rho > rho_peak) rho_peak = rho;

            int ip=(i+1)%N, im=(i-1+N)%N, jp=(j+1)%N, jm=(j-1+N)%N, kp=(k+1)%N, km=(k-1+N)%N;
            long nip=(long)ip*NN+j*N+k, nim=(long)im*NN+j*N+k;
            long njp=(long)i*NN+jp*N+k, njm=(long)i*NN+jm*N+k;
            long nkp=(long)i*NN+j*N+kp, nkm=(long)i*NN+j*N+km;

            /* Energies */
            for (int a=0;a<3;a++) {
                ek_phi += 0.5*phi_vel[a][idx]*phi_vel[a][idx]*dV;
                ek_theta += 0.5*theta_vel[a][idx]*theta_vel[a][idx]*dV;
                double gx=(phi[a][nip]-phi[a][nim])*idx1;
                double gy=(phi[a][njp]-phi[a][njm])*idx1;
                double gz=(phi[a][nkp]-phi[a][nkm])*idx1;
                eg_phi += 0.5*(gx*gx+gy*gy+gz*gz)*dV;
                e_mass += 0.5*MASS2*phi[a][idx]*phi[a][idx]*dV;
                double tgx=(theta[a][nip]-theta[a][nim])*idx1;
                double tgy=(theta[a][njp]-theta[a][njm])*idx1;
                double tgz=(theta[a][nkp]-theta[a][nkm])*idx1;
                eg_theta += 0.5*(tgx*tgx+tgy*tgy+tgz*tgz)*dV;
                theta2 += theta[a][idx]*theta[a][idx];
            }
            e_pot += (MU/2.0)*P2/(1.0+KAPPA*P2)*dV;
            double ct0=(theta[2][njp]-theta[2][njm]-theta[1][nkp]+theta[1][nkm])*idx1;
            double ct1=(theta[0][nkp]-theta[0][nkm]-theta[2][nip]+theta[2][nim])*idx1;
            double ct2=(theta[1][nip]-theta[1][nim]-theta[0][njp]+theta[0][njm])*idx1;
            e_coupling -= ETA*(p0*ct0+p1*ct1+p2*ct2)*dV;

            P_int += Pa*dV;
            if (Pa > P_half) core_vol += dV;

            double vx=phi_vel[0][idx], vy=phi_vel[1][idx], vz=phi_vel[2][idx];
            vdrift[0]+=vx; vdrift[1]+=vy; vdrift[2]+=vz;
            v2sum += vx*vx+vy*vy+vz*vz;
            double dvx=(phi_vel[0][nip]-phi_vel[0][nim])*idx1;
            double dvy=(phi_vel[1][njp]-phi_vel[1][njm])*idx1;
            double dvz=(phi_vel[2][nkp]-phi_vel[2][nkm])*idx1;
            divv_sum += dvx+dvy+dvz;

            double w = rho*dV;
            Ixx+=(ry*ry+rz*rz)*w; Iyy+=(rx*rx+rz*rz)*w; Izz+=(rx*rx+ry*ry)*w;

            if (r < Rmax*0.5) { theta2_core += theta[0][idx]*theta[0][idx]+theta[1][idx]*theta[1][idx]+theta[2][idx]*theta[2][idx]; t_core_cnt++; }

            /* Radial profile */
            int bin = (int)(r/dr);
            if (bin >= RADIAL_BINS) bin = RADIAL_BINS-1;
            ci->prof_rho[bin] += rho;
            ci->prof_P[bin] += Pa;
            ci->prof_theta[bin] += theta[0][idx]*theta[0][idx]+theta[1][idx]*theta[1][idx]+theta[2][idx]*theta[2][idx];
            ci->prof_Ekin[bin] += (0.5*(vx*vx+vy*vy+vz*vz) + 0.5*(theta_vel[0][idx]*theta_vel[0][idx]+theta_vel[1][idx]*theta_vel[1][idx]+theta_vel[2][idx]*theta_vel[2][idx]))*dV;
            ci->prof_Epot[bin] += (MU/2.0)*P2/(1.0+KAPPA*P2)*dV;
            ci->prof_vmag[bin] += sqrt(vx*vx+vy*vy+vz*vz);

            /* Acceleration estimate: Laplacian + forces */
            double a2 = 0;
            for (int a=0;a<3;a++) {
                double lap = (phi[a][nip]+phi[a][nim]+phi[a][njp]+phi[a][njm]+phi[a][nkp]+phi[a][nkm]-6.0*phi[a][idx])/(dx*dx);
                double dPda = (a==0)?p1*p2:(a==1)?p0*p2:p0*p1;
                double denom = (1.0+KAPPA*P2);
                double dVda = MU*P*dPda/(denom*denom);
                double curl_a = (a==0)?ct0:(a==1)?ct1:ct2;
                double F_a = lap - MASS2*phi[a][idx] - dVda + ETA*curl_a;
                a2 += F_a*F_a;
            }
            ci->prof_amag[bin] += sqrt(a2);
            ci->prof_count[bin]++;
        }

        ci->E_kin_phi = ek_phi; ci->E_kin_theta = ek_theta;
        ci->E_grad_phi = eg_phi; ci->E_grad_theta = eg_theta;
        ci->E_mass = e_mass; ci->E_pot = e_pot; ci->E_coupling = e_coupling;
        ci->E_total = ek_phi+ek_theta+eg_phi+eg_theta+e_mass+e_pot+e_coupling;
        ci->P_int = P_int; ci->core_volume = core_vol; ci->rho_peak = rho_peak;
        ci->theta_rms = (cnt>0) ? sqrt(theta2/(3.0*cnt)) : 0;
        ci->theta_rms_core = (t_core_cnt>0) ? sqrt(theta2_core/(3.0*t_core_cnt)) : 0;

        if (cnt > 0) {
            vdrift[0]/=cnt; vdrift[1]/=cnt; vdrift[2]/=cnt;
            ci->v_drift_mag = sqrt(vdrift[0]*vdrift[0]+vdrift[1]*vdrift[1]+vdrift[2]*vdrift[2]);
            ci->v_disp = sqrt(fmax(0, v2sum/cnt - ci->v_drift_mag*ci->v_drift_mag));
            ci->div_v = divv_sum/cnt;
        }

        double Imax=Ixx,Imin=Ixx;
        if(Iyy>Imax)Imax=Iyy; if(Iyy<Imin)Imin=Iyy;
        if(Izz>Imax)Imax=Izz; if(Izz<Imin)Imin=Izz;
        ci->aspect = (Imin>1e-30)?Imax/Imin:999;

        /* Surface: find R_half from profile */
        ci->R_half = ci->R_max;
        for (int b = 0; b < RADIAL_BINS; b++) {
            if (ci->prof_count[b] > 0) {
                ci->prof_rho[b] /= ci->prof_count[b];
                ci->prof_P[b] /= ci->prof_count[b];
                ci->prof_theta[b] = sqrt(ci->prof_theta[b]/(3.0*ci->prof_count[b]));
                ci->prof_vmag[b] /= ci->prof_count[b];
                ci->prof_amag[b] /= ci->prof_count[b];
                /* E_kin and E_pot stay as totals */
            }
        }
        for (int b = 1; b < RADIAL_BINS; b++) {
            if (ci->prof_P[b] < P_half && ci->prof_P[b-1] >= P_half) {
                double r0 = (b-0.5)*dr, r1 = (b+0.5)*dr;
                ci->R_half = 0.5*(r0+r1);
                ci->dPdr_surface = (ci->prof_P[b]-ci->prof_P[b-1])/dr;
                break;
            }
        }
    }

    free(labels);
    for (int a=0;a<3;a++) { free(phi[a]); free(theta[a]); free(phi_vel[a]); free(theta_vel[a]); }
}

/* ---- Statistical summary ---- */
typedef struct {
    double mean, std, min, max;
    int    n;
} Stat;

static Stat compute_stat(double *vals, int n) {
    Stat s = {0,0,DBL_MAX,-DBL_MAX,n};
    if (n == 0) { s.min=0; s.max=0; return s; }
    for (int i=0;i<n;i++) {
        s.mean += vals[i];
        if (vals[i]<s.min) s.min=vals[i];
        if (vals[i]>s.max) s.max=vals[i];
    }
    s.mean /= n;
    for (int i=0;i<n;i++) s.std += (vals[i]-s.mean)*(vals[i]-s.mean);
    s.std = (n>1) ? sqrt(s.std/(n-1)) : 0;
    return s;
}

static void print_stat(FILE *fp, const char *name, Stat s, int last) {
    fprintf(fp, "      \"%s\": {\"mean\": %.6e, \"std\": %.6e, \"min\": %.6e, \"max\": %.6e, \"n\": %d}%s\n",
            name, s.mean, s.std, s.min, s.max, s.n, last?"":",");
}

/* ---- main ---- */
int main(int argc, char **argv)
{
    if (argc < 2) {
        fprintf(stderr, "Usage: %s input.sfa --early F1 --late F2 [--json out.json]\n", argv[0]);
        return 1;
    }

    const char *sfa_path = argv[1];
    const char *json_path = NULL;
    int early_frame = 0, late_frame = 40;

    for (int i = 2; i < argc; i++) {
        if (!strcmp(argv[i],"--early") && i+1<argc) early_frame = atoi(argv[++i]);
        if (!strcmp(argv[i],"--late")  && i+1<argc) late_frame  = atoi(argv[++i]);
        if (!strcmp(argv[i],"--json")  && i+1<argc) json_path   = argv[++i];
    }

    SFA *sfa = sfa_open(sfa_path);
    if (!sfa) { fprintf(stderr, "Cannot open %s\n", sfa_path); return 1; }

    int N = sfa->Nx; double L = sfa->Lx;

    /* Read physics from KVMD */
    SFA_KVMDSet kv[16];
    int n_kv = sfa_read_kvmd(sfa, kv, 16);
    for (int s=0;s<n_kv;s++) for(int p=0;p<kv[s].n_pairs;p++) {
        if(!strcmp(kv[s].keys[p],"mu"))    MU=atof(kv[s].values[p]);
        if(!strcmp(kv[s].keys[p],"kappa")) KAPPA=atof(kv[s].values[p]);
        if(!strcmp(kv[s].keys[p],"m"))   { double m=atof(kv[s].values[p]); MASS2=m*m; }
        if(!strcmp(kv[s].keys[p],"eta"))   ETA=atof(kv[s].values[p]);
    }

    if (early_frame >= (int)sfa->total_frames) early_frame = 0;
    if (late_frame  >= (int)sfa->total_frames) late_frame  = sfa->total_frames - 1;

    fprintf(stderr, "stable_vs_unstable: %s (N=%d L=%.1f)\n", sfa_path, N, L);
    fprintf(stderr, "Early frame: %d, Late frame: %d\n", early_frame, late_frame);

    void *buf = malloc(sfa->frame_bytes);

    /* Analyze early frame */
    FrameInfo early, late;
    memset(&early, 0, sizeof(early));
    memset(&late, 0, sizeof(late));

    sfa_read_frame(sfa, early_frame, buf);
    early.time = sfa_frame_time(sfa, early_frame);
    early.frame_idx = early_frame;
    fprintf(stderr, "Analyzing early frame %d (t=%.1f)...\n", early_frame, early.time);
    analyze_frame(sfa, buf, N, L, &early);
    fprintf(stderr, "  %d clusters found\n", early.n_clusters);

    sfa_read_frame(sfa, late_frame, buf);
    late.time = sfa_frame_time(sfa, late_frame);
    late.frame_idx = late_frame;
    fprintf(stderr, "Analyzing late frame %d (t=%.1f)...\n", late_frame, late.time);
    analyze_frame(sfa, buf, N, L, &late);
    fprintf(stderr, "  %d clusters found\n", late.n_clusters);

    free(buf);

    /* Match early -> late clusters by nearest centroid */
    int n_early = early.n_clusters, n_late = late.n_clusters;
    MatchedPair matches[MAX_CLUSTERS];
    int n_matches = 0;
    int *late_matched = (int*)calloc(n_late, sizeof(int));

    for (int e = 0; e < n_early && n_matches < MAX_CLUSTERS; e++) {
        ClusterInfo *ec = &early.clusters[e];
        double best_dist = MATCH_RADIUS;
        int best_l = -1;
        for (int l = 0; l < n_late; l++) {
            if (late_matched[l]) continue;
            ClusterInfo *lc = &late.clusters[l];
            double d = sqrt((ec->cx-lc->cx)*(ec->cx-lc->cx) +
                           (ec->cy-lc->cy)*(ec->cy-lc->cy) +
                           (ec->cz-lc->cz)*(ec->cz-lc->cz));
            if (d < best_dist) { best_dist = d; best_l = l; }
        }
        if (best_l >= 0) {
            late_matched[best_l] = 1;
            MatchedPair *mp = &matches[n_matches++];
            mp->early_id = e;
            mp->late_id = best_l;
            mp->distance = best_dist;
            mp->dP_peak = late.clusters[best_l].P_peak - ec->P_peak;
            mp->dP_int = late.clusters[best_l].P_int - ec->P_int;
            mp->dn_voxels = late.clusters[best_l].n_voxels - ec->n_voxels;
            /* Classify: grew = P_int increased OR voxel count increased */
            mp->grew = (mp->dP_int > 0 || mp->dn_voxels > 0) ? 1 : 0;
        }
    }
    free(late_matched);

    /* Unmatched early clusters = disappeared (definitely unstable) */
    /* Unmatched late clusters = newly formed */

    /* Count stable vs unstable */
    int n_stable = 0, n_unstable = 0;
    for (int m = 0; m < n_matches; m++) {
        if (matches[m].grew) n_stable++; else n_unstable++;
    }
    int n_disappeared = n_early - n_matches;
    int n_new = n_late - n_matches;  /* some late clusters unmatched */
    /* Count how many late clusters were matched */
    int late_match_count = 0;
    for (int m = 0; m < n_matches; m++) late_match_count++;
    n_new = n_late - late_match_count;

    fprintf(stderr, "\nMatching results:\n");
    fprintf(stderr, "  Matched pairs: %d\n", n_matches);
    fprintf(stderr, "  Stable (grew):  %d\n", n_stable);
    fprintf(stderr, "  Unstable (shrank): %d\n", n_unstable);
    fprintf(stderr, "  Disappeared: %d\n", n_disappeared);
    fprintf(stderr, "  Newly formed: %d\n", n_new);

    /* Collect metric arrays for stable vs unstable groups */
    #define NMETRICS 16
    const char *metric_names[NMETRICS] = {
        "n_voxels", "P_peak", "P_int", "core_volume",
        "E_pot", "E_total", "E_kin_phi", "E_grad_phi",
        "theta_rms", "theta_rms_core", "v_drift_mag", "v_disp",
        "div_v", "aspect", "rho_peak", "R_half"
    };

    double stable_vals[NMETRICS][MAX_CLUSTERS];
    double unstable_vals[NMETRICS][MAX_CLUSTERS];

    for (int m = 0; m < n_matches; m++) {
        ClusterInfo *ci = &early.clusters[matches[m].early_id];
        double *dst = matches[m].grew ? stable_vals[0] : unstable_vals[0];
        int idx = matches[m].grew ? 0 : 0;  /* We need separate counters */
        /* Use a simpler approach: collect into flat arrays */
        (void)dst; (void)idx;
    }

    /* Simpler: just collect directly */
    int si = 0, ui = 0;
    for (int m = 0; m < n_matches; m++) {
        ClusterInfo *ci = &early.clusters[matches[m].early_id];
        int *pi = matches[m].grew ? &si : &ui;
        double (*arr)[MAX_CLUSTERS] = matches[m].grew ? stable_vals : unstable_vals;
        int k = *pi;
        arr[0][k] = ci->n_voxels;
        arr[1][k] = ci->P_peak;
        arr[2][k] = ci->P_int;
        arr[3][k] = ci->core_volume;
        arr[4][k] = ci->E_pot;
        arr[5][k] = ci->E_total;
        arr[6][k] = ci->E_kin_phi;
        arr[7][k] = ci->E_grad_phi;
        arr[8][k] = ci->theta_rms;
        arr[9][k] = ci->theta_rms_core;
        arr[10][k] = ci->v_drift_mag;
        arr[11][k] = ci->v_disp;
        arr[12][k] = ci->div_v;
        arr[13][k] = ci->aspect;
        arr[14][k] = ci->rho_peak;
        arr[15][k] = ci->R_half;
        (*pi)++;
    }

    /* Write JSON */
    FILE *jfp = stdout;
    if (json_path) {
        jfp = fopen(json_path, "w");
        if (!jfp) { fprintf(stderr, "Cannot write %s\n", json_path); jfp = stdout; }
    }

    fprintf(jfp, "{\n");
    fprintf(jfp, "  \"sfa\": \"%s\",\n", sfa_path);
    fprintf(jfp, "  \"N\": %d, \"L\": %.2f,\n", N, L);
    fprintf(jfp, "  \"physics\": {\"mu\": %.6f, \"kappa\": %.1f, \"m2\": %.4f, \"eta\": %.3f},\n",
            MU, KAPPA, MASS2, ETA);
    fprintf(jfp, "  \"early\": {\"frame\": %d, \"time\": %.2f, \"n_clusters\": %d, \"P_max\": %.8e},\n",
            early_frame, early.time, n_early, early.P_max_global);
    fprintf(jfp, "  \"late\": {\"frame\": %d, \"time\": %.2f, \"n_clusters\": %d, \"P_max\": %.8e},\n",
            late_frame, late.time, n_late, late.P_max_global);
    fprintf(jfp, "  \"matching\": {\n");
    fprintf(jfp, "    \"n_matched\": %d,\n", n_matches);
    fprintf(jfp, "    \"n_stable\": %d,\n", n_stable);
    fprintf(jfp, "    \"n_unstable\": %d,\n", n_unstable);
    fprintf(jfp, "    \"n_disappeared\": %d,\n", n_disappeared);
    fprintf(jfp, "    \"n_new\": %d,\n", n_new);
    fprintf(jfp, "    \"match_radius\": %.1f\n", MATCH_RADIUS);
    fprintf(jfp, "  },\n");

    /* Matched pairs detail */
    fprintf(jfp, "  \"pairs\": [\n");
    for (int m = 0; m < n_matches; m++) {
        MatchedPair *mp = &matches[m];
        ClusterInfo *ec = &early.clusters[mp->early_id];
        ClusterInfo *lc = &late.clusters[mp->late_id];
        fprintf(jfp, "    {\n");
        fprintf(jfp, "      \"early_centroid\": [%.2f, %.2f, %.2f],\n", ec->cx, ec->cy, ec->cz);
        fprintf(jfp, "      \"late_centroid\": [%.2f, %.2f, %.2f],\n", lc->cx, lc->cy, lc->cz);
        fprintf(jfp, "      \"distance\": %.4f,\n", mp->distance);
        fprintf(jfp, "      \"grew\": %s,\n", mp->grew?"true":"false");
        fprintf(jfp, "      \"early\": {\"n_voxels\": %d, \"P_peak\": %.6e, \"P_int\": %.6e, "
                "\"E_pot\": %.6e, \"theta_rms\": %.6e, \"R_half\": %.4f, "
                "\"v_disp\": %.6e, \"div_v\": %.6e, \"aspect\": %.4f},\n",
                ec->n_voxels, ec->P_peak, ec->P_int, ec->E_pot, ec->theta_rms,
                ec->R_half, ec->v_disp, ec->div_v, ec->aspect);
        fprintf(jfp, "      \"late\": {\"n_voxels\": %d, \"P_peak\": %.6e, \"P_int\": %.6e, "
                "\"E_pot\": %.6e, \"theta_rms\": %.6e, \"R_half\": %.4f, "
                "\"v_disp\": %.6e, \"div_v\": %.6e, \"aspect\": %.4f},\n",
                lc->n_voxels, lc->P_peak, lc->P_int, lc->E_pot, lc->theta_rms,
                lc->R_half, lc->v_disp, lc->div_v, lc->aspect);
        fprintf(jfp, "      \"delta\": {\"dP_peak\": %.6e, \"dP_int\": %.6e, \"dn_voxels\": %.0f}\n",
                mp->dP_peak, mp->dP_int, mp->dn_voxels);
        fprintf(jfp, "    }%s\n", (m<n_matches-1)?",":"");
    }
    fprintf(jfp, "  ],\n");

    /* Statistical comparison: stable vs unstable (early-time metrics) */
    fprintf(jfp, "  \"stable_stats\": {\n");
    for (int mi = 0; mi < NMETRICS; mi++) {
        Stat s = compute_stat(stable_vals[mi], si);
        print_stat(jfp, metric_names[mi], s, mi==NMETRICS-1);
    }
    fprintf(jfp, "  },\n");
    fprintf(jfp, "  \"unstable_stats\": {\n");
    for (int mi = 0; mi < NMETRICS; mi++) {
        Stat s = compute_stat(unstable_vals[mi], ui);
        print_stat(jfp, metric_names[mi], s, mi==NMETRICS-1);
    }
    fprintf(jfp, "  },\n");

    /* Per-cluster early-time data for all clusters (for correlation analysis) */
    fprintf(jfp, "  \"early_clusters\": [\n");
    for (int e = 0; e < n_early; e++) {
        ClusterInfo *ci = &early.clusters[e];
        /* Find if this was matched and whether it grew */
        int matched = 0, grew = -1;
        for (int m = 0; m < n_matches; m++) {
            if (matches[m].early_id == e) { matched = 1; grew = matches[m].grew; break; }
        }
        fprintf(jfp, "    {\"id\": %d, \"n_voxels\": %d, \"centroid\": [%.2f,%.2f,%.2f], "
                "\"P_peak\": %.6e, \"P_int\": %.6e, \"E_pot\": %.6e, "
                "\"theta_rms\": %.6e, \"v_disp\": %.6e, \"div_v\": %.6e, "
                "\"aspect\": %.4f, \"R_half\": %.4f, \"rho_peak\": %.6e, "
                "\"matched\": %s, \"grew\": %s}%s\n",
                ci->id, ci->n_voxels, ci->cx, ci->cy, ci->cz,
                ci->P_peak, ci->P_int, ci->E_pot, ci->theta_rms,
                ci->v_disp, ci->div_v, ci->aspect, ci->R_half, ci->rho_peak,
                matched?"true":"false",
                (grew==1)?"true":(grew==0)?"false":"null",
                (e<n_early-1)?",":"");
    }
    fprintf(jfp, "  ],\n");

    /* Late clusters */
    fprintf(jfp, "  \"late_clusters\": [\n");
    for (int l = 0; l < n_late; l++) {
        ClusterInfo *ci = &late.clusters[l];
        int matched = 0;
        for (int m = 0; m < n_matches; m++)
            if (matches[m].late_id == l) { matched = 1; break; }
        fprintf(jfp, "    {\"id\": %d, \"n_voxels\": %d, \"centroid\": [%.2f,%.2f,%.2f], "
                "\"P_peak\": %.6e, \"P_int\": %.6e, \"E_pot\": %.6e, "
                "\"theta_rms\": %.6e, \"v_disp\": %.6e, \"div_v\": %.6e, "
                "\"aspect\": %.4f, \"R_half\": %.4f, \"rho_peak\": %.6e, "
                "\"matched\": %s}%s\n",
                ci->id, ci->n_voxels, ci->cx, ci->cy, ci->cz,
                ci->P_peak, ci->P_int, ci->E_pot, ci->theta_rms,
                ci->v_disp, ci->div_v, ci->aspect, ci->R_half, ci->rho_peak,
                matched?"true":"false",
                (l<n_late-1)?",":"");
    }
    fprintf(jfp, "  ]\n");

    fprintf(jfp, "}\n");

    if (json_path && jfp != stdout) {
        fclose(jfp);
        fprintf(stderr, "JSON written to %s\n", json_path);
    }

    sfa_close(sfa);
    return 0;
}
