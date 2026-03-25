/*  freq_phase.c — Frequency and phase analysis of SFA time series
 *
 *  For multi-frame SFA files, computes:
 *  1. Time series of per-cluster phase (carrier wave phase extraction)
 *  2. FFT-free frequency analysis via autocorrelation of E_pot, P_int, theta_rms
 *  3. Inter-cluster phase relationships (are clusters phase-locked?)
 *  4. Breathing mode period detection (peak-to-peak in E_pot)
 *
 *  Build: cd sfa/analysis && make freq_phase
 *  Usage: ./freq_phase input.sfa [--json output.json]
 */

#define SFA_IMPLEMENTATION
#include "../format/sfa.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>

#define PI 3.14159265358979323846
#define MAX_FRAMES 200
#define MAX_CLUSTERS 32

static const double MU = -41.345, KAPPA = 50.0, MASS2 = 2.25;

static inline double f16_to_f64(uint16_t h) {
    uint16_t sign = h & 0x8000; int exp = (h >> 10) & 0x1F; uint16_t mant = h & 0x3FF;
    if (exp == 0) return 0.0; if (exp == 31) return sign ? -1e30 : 1e30;
    float fv; uint32_t x = ((uint32_t)sign << 16) | ((uint32_t)(exp-15+127) << 23) | ((uint32_t)mant << 13);
    memcpy(&fv, &x, 4); return (double)fv;
}

/* Extract field array from SFA frame buffer, handling all dtypes */
static double *extract_column(void *buf, SFA *sfa, int target_sem, int target_comp, long N3) {
    double *arr = (double*)calloc(N3, sizeof(double));
    uint64_t off = 0;
    for (int c = 0; c < sfa->n_columns; c++) {
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
    return arr;  /* zeros if not found */
}

/* BFS cluster detection */
static int find_clusters(double *phi[3], int N, double threshold,
                         int *labels, double centroids[][3], int max_cl) {
    long N3 = (long)N*N*N; int NN = N*N;
    double dx = 1.0;  /* normalized */
    memset(labels, -1, N3 * sizeof(int));
    long *queue = (long*)malloc(N3 * sizeof(long));
    int n_clusters = 0;

    for (long idx = 0; idx < N3; idx++) {
        double P = fabs(phi[0][idx]*phi[1][idx]*phi[2][idx]);
        if (P < threshold || labels[idx] >= 0) continue;
        if (n_clusters >= max_cl) break;
        int cid = n_clusters++;
        double cx=0,cy=0,cz=0,wt=0;
        int front=0, back=0;
        queue[back++] = idx; labels[idx] = cid;
        while (front < back) {
            long cur = queue[front++];
            int i=(int)(cur/NN),j=(int)((cur/N)%N),k=(int)(cur%N);
            double w = fabs(phi[0][cur]*phi[1][cur]*phi[2][cur]);
            cx+=i*w; cy+=j*w; cz+=k*w; wt+=w;
            int di[]={1,-1,0,0,0,0}, dj[]={0,0,1,-1,0,0}, dk[]={0,0,0,0,1,-1};
            for(int d=0;d<6;d++){
                int ni=(i+di[d]+N)%N, nj=(j+dj[d]+N)%N, nk=(k+dk[d]+N)%N;
                long nidx=(long)ni*NN+nj*N+nk;
                if(labels[nidx]<0 && fabs(phi[0][nidx]*phi[1][nidx]*phi[2][nidx])>=threshold){
                    labels[nidx]=cid; queue[back++]=nidx;
                }
            }
        }
        if(wt>0){centroids[cid][0]=cx/wt;centroids[cid][1]=cy/wt;centroids[cid][2]=cz/wt;}
    }
    free(queue);
    return n_clusters;
}

/* Measure local phase of the carrier wave at a cluster centroid */
static double measure_phase(double *phi[3], int N, double cx, double cy, double cz) {
    int NN = N*N;
    int ic=(int)(cx+0.5), jc=(int)(cy+0.5), kc=(int)(cz+0.5);
    if(ic<0)ic=0;if(ic>=N)ic=N-1;
    if(jc<0)jc=0;if(jc>=N)jc=N-1;
    if(kc<0)kc=0;if(kc>=N)kc=N-1;
    long idx = (long)ic*NN + jc*N + kc;
    /* Phase from atan2(phi_1, phi_0) at centroid */
    return atan2(phi[1][idx], phi[0][idx]);
}

/* Autocorrelation for period detection */
static double find_period(double *series, int n, double dt) {
    if (n < 4) return 0;
    /* Remove mean */
    double mean = 0;
    for(int i=0;i<n;i++) mean+=series[i];
    mean /= n;
    double *s = (double*)malloc(n*sizeof(double));
    for(int i=0;i<n;i++) s[i]=series[i]-mean;

    /* Autocorrelation */
    int max_lag = n/2;
    double *ac = (double*)calloc(max_lag, sizeof(double));
    double ac0 = 0;
    for(int i=0;i<n;i++) ac0+=s[i]*s[i];
    if(ac0<1e-30){free(s);free(ac);return 0;}
    for(int lag=0;lag<max_lag;lag++){
        for(int i=0;i<n-lag;i++) ac[lag]+=s[i]*s[i+lag];
        ac[lag]/=ac0;
    }

    /* Find first peak after first zero crossing */
    int zero_cross = -1;
    for(int i=1;i<max_lag;i++) if(ac[i]<0){zero_cross=i;break;}
    double best_peak = 0; int best_lag = 0;
    if(zero_cross>0){
        for(int i=zero_cross;i<max_lag;i++){
            if(ac[i]>best_peak){best_peak=ac[i];best_lag=i;}
            if(ac[i]<best_peak*0.5 && best_lag>0) break;  /* past the peak */
        }
    }

    free(s); free(ac);
    return (best_lag > 0) ? best_lag * dt : 0;
}

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "Usage: %s input.sfa [--json out.json]\n", argv[0]); return 1; }
    const char *sfa_path = argv[1];
    const char *json_path = NULL;
    for(int i=2;i<argc;i++) if(!strcmp(argv[i],"--json")&&i+1<argc) json_path=argv[++i];

    SFA *sfa = sfa_open(sfa_path);
    if (!sfa) { fprintf(stderr, "Cannot open %s\n", sfa_path); return 1; }
    int N = sfa->Nx; double L = sfa->Lx;
    long N3 = (long)N*N*N;
    double dx = 2.0*L/(N-1);

    int nframes = sfa->total_frames;
    if (nframes > MAX_FRAMES) nframes = MAX_FRAMES;
    printf("Freq/phase analysis: %s (N=%d, %d frames)\n\n", sfa_path, N, nframes);

    /* Time series storage */
    double *ts_time = (double*)calloc(nframes, sizeof(double));
    double *ts_epot = (double*)calloc(nframes, sizeof(double));
    double *ts_pint = (double*)calloc(nframes, sizeof(double));
    double *ts_trms = (double*)calloc(nframes, sizeof(double));
    double *ts_phimax = (double*)calloc(nframes, sizeof(double));
    int *ts_nclusters = (int*)calloc(nframes, sizeof(int));

    /* Per-cluster phase tracking (up to 8 clusters) */
    double cluster_phases[MAX_FRAMES][MAX_CLUSTERS];
    double cluster_centroids[MAX_FRAMES][MAX_CLUSTERS][3];
    double cluster_Pint[MAX_FRAMES][MAX_CLUSTERS];
    memset(cluster_phases, 0, sizeof(cluster_phases));
    memset(cluster_Pint, 0, sizeof(cluster_Pint));

    void *buf = malloc(sfa->frame_bytes);
    int *labels = (int*)malloc(N3 * sizeof(int));

    for (int f = 0; f < nframes; f++) {
        sfa_read_frame(sfa, f, buf);
        ts_time[f] = sfa_frame_time(sfa, f);

        double *phi[3], *theta[3];
        for(int a=0;a<3;a++){
            phi[a] = extract_column(buf, sfa, SFA_POSITION, a, N3);
            theta[a] = extract_column(buf, sfa, SFA_ANGLE, a, N3);
        }

        /* Global metrics */
        double epot=0, pint=0, trms=0, phimax=0;
        double Pmax = 0;
        for(long i=0;i<N3;i++){
            double P = phi[0][i]*phi[1][i]*phi[2][i];
            double P2=P*P, Pabs=fabs(P);
            epot += (MU/2.0)*P2/(1.0+KAPPA*P2)*dx*dx*dx;
            pint += Pabs*dx*dx*dx;
            for(int a=0;a<3;a++){
                trms += theta[a][i]*theta[a][i];
                double ap=fabs(phi[a][i]); if(ap>phimax) phimax=ap;
            }
            if(Pabs>Pmax) Pmax=Pabs;
        }
        ts_epot[f] = epot;
        ts_pint[f] = pint;
        ts_trms[f] = sqrt(trms/(3.0*N3));
        ts_phimax[f] = phimax;

        /* Cluster detection */
        double cents[MAX_CLUSTERS][3];
        int nc = find_clusters(phi, N, 0.01*Pmax, labels, cents, MAX_CLUSTERS);
        ts_nclusters[f] = nc;

        /* Per-cluster phase and P_int */
        for(int c=0;c<nc && c<MAX_CLUSTERS;c++){
            cluster_centroids[f][c][0] = cents[c][0];
            cluster_centroids[f][c][1] = cents[c][1];
            cluster_centroids[f][c][2] = cents[c][2];
            cluster_phases[f][c] = measure_phase(phi, N, cents[c][0], cents[c][1], cents[c][2]);
            /* Per-cluster P_int */
            double cp = 0;
            for(long i=0;i<N3;i++){
                if(labels[i]==c) cp += fabs(phi[0][i]*phi[1][i]*phi[2][i])*dx*dx*dx;
            }
            cluster_Pint[f][c] = cp;
        }

        printf("  f=%2d t=%6.1f Ep=%8.1f Pint=%6.1f trms=%.4f phi=%.3f ncl=%d",
               f, ts_time[f], epot, pint, ts_trms[f], phimax, nc);
        if(nc>0) printf(" phases=[");
        for(int c=0;c<nc&&c<4;c++) printf("%.2f%s", cluster_phases[f][c], c<nc-1&&c<3?",":"");
        if(nc>0) printf("]");
        printf("\n");

        for(int a=0;a<3;a++){free(phi[a]);free(theta[a]);}
    }

    /* Period analysis */
    double dt_frames = (nframes > 1) ? ts_time[1] - ts_time[0] : 1;
    double period_epot = find_period(ts_epot, nframes, dt_frames);
    double period_pint = find_period(ts_pint, nframes, dt_frames);
    double period_trms = find_period(ts_trms, nframes, dt_frames);
    double period_phi  = find_period(ts_phimax, nframes, dt_frames);

    printf("\n=== FREQUENCY ANALYSIS ===\n");
    printf("Breathing period (E_pot):    %.1f time units (freq=%.4f)\n",
           period_epot, period_epot>0?1.0/period_epot:0);
    printf("Triple product period:       %.1f time units (freq=%.4f)\n",
           period_pint, period_pint>0?1.0/period_pint:0);
    printf("Theta oscillation period:    %.1f time units (freq=%.4f)\n",
           period_trms, period_trms>0?1.0/period_trms:0);
    printf("Phi amplitude period:        %.1f time units (freq=%.4f)\n",
           period_phi, period_phi>0?1.0/period_phi:0);

    /* Inter-cluster phase analysis */
    printf("\n=== PHASE RELATIONSHIPS ===\n");
    /* For each pair of frames, compute phase differences between clusters */
    if (nframes > 1 && ts_nclusters[0] >= 2) {
        printf("Phase differences between largest clusters (frame 0 vs last):\n");
        int nc0 = ts_nclusters[0], ncL = ts_nclusters[nframes-1];
        printf("  Frame 0 (%d clusters): ", nc0);
        for(int c=0;c<nc0&&c<4;c++) printf("φ%d=%.3f ", c, cluster_phases[0][c]);
        printf("\n");
        printf("  Frame %d (%d clusters): ", nframes-1, ncL);
        for(int c=0;c<ncL&&c<4;c++) printf("φ%d=%.3f ", c, cluster_phases[nframes-1][c]);
        printf("\n");

        /* Check if phase differences are preserved */
        if(nc0>=2 && ncL>=2){
            double dp0 = cluster_phases[0][1] - cluster_phases[0][0];
            double dpL = cluster_phases[nframes-1][1] - cluster_phases[nframes-1][0];
            printf("  Δφ(0→1) at t=0: %.3f, at t=%.0f: %.3f (drift=%.3f)\n",
                   dp0, ts_time[nframes-1], dpL, dpL-dp0);
        }
        if(nc0>=3 && ncL>=3){
            double dp0 = cluster_phases[0][2] - cluster_phases[0][0];
            double dpL = cluster_phases[nframes-1][2] - cluster_phases[nframes-1][0];
            printf("  Δφ(0→2) at t=0: %.3f, at t=%.0f: %.3f (drift=%.3f)\n",
                   dp0, ts_time[nframes-1], dpL, dpL-dp0);
        }
    }

    /* Cluster count stability */
    printf("\n=== MULTI-STRUCTURE PATTERNS ===\n");
    printf("Cluster count over time: ");
    for(int f=0;f<nframes;f++) printf("%d%s", ts_nclusters[f], f<nframes-1?" → ":"");
    printf("\n");

    /* Per-cluster P_int evolution */
    if(ts_nclusters[0] >= 2){
        printf("\nPer-cluster P_int evolution (top 3 clusters):\n");
        printf("  %-6s", "t");
        for(int c=0;c<3;c++) printf("  cl%d_Pint", c);
        printf("\n");
        for(int f=0;f<nframes;f++){
            printf("  %-6.0f", ts_time[f]);
            for(int c=0;c<3&&c<ts_nclusters[f];c++) printf("  %8.1f", cluster_Pint[f][c]);
            printf("\n");
        }
    }

    /* JSON output */
    if (json_path) {
        FILE *jfp = fopen(json_path, "w");
        fprintf(jfp, "{\n  \"sfa\": \"%s\",\n  \"N\": %d, \"n_frames\": %d,\n", sfa_path, N, nframes);
        fprintf(jfp, "  \"periods\": {\n");
        fprintf(jfp, "    \"E_pot\": %.2f,\n    \"P_int\": %.2f,\n    \"theta_rms\": %.2f,\n    \"phi_max\": %.2f\n  },\n",
                period_epot, period_pint, period_trms, period_phi);
        fprintf(jfp, "  \"frames\": [\n");
        for(int f=0;f<nframes;f++){
            fprintf(jfp, "    {\"t\":%.1f,\"Ep\":%.2f,\"Pint\":%.2f,\"trms\":%.5f,\"phi_max\":%.4f,\"n_clusters\":%d,\"phases\":[",
                    ts_time[f],ts_epot[f],ts_pint[f],ts_trms[f],ts_phimax[f],ts_nclusters[f]);
            for(int c=0;c<ts_nclusters[f]&&c<4;c++)
                fprintf(jfp,"%.4f%s",cluster_phases[f][c],c<ts_nclusters[f]-1&&c<3?",":"");
            fprintf(jfp,"],\"cluster_Pint\":[");
            for(int c=0;c<ts_nclusters[f]&&c<4;c++)
                fprintf(jfp,"%.2f%s",cluster_Pint[f][c],c<ts_nclusters[f]-1&&c<3?",":"");
            fprintf(jfp,"]}%s\n", f<nframes-1?",":"");
        }
        fprintf(jfp, "  ]\n}\n");
        fclose(jfp);
        printf("\nJSON written to %s\n", json_path);
    }

    free(buf); free(labels);
    free(ts_time); free(ts_epot); free(ts_pint); free(ts_trms); free(ts_phimax); free(ts_nclusters);
    sfa_close(sfa);
    return 0;
}
