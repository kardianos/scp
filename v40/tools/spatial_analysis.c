/*  spatial_analysis.c — Analyze spatial structure of SFA frames
 *
 *  For each frame: compute centroid, rms radius, octant energy distribution,
 *  theta stability (rms variation), and cluster count via connected components.
 *
 *  Build: gcc -O3 -fopenmp -o spatial_analysis spatial_analysis.c -lzstd -lm
 *  Usage: ./spatial_analysis input.sfa
 */

#define SFA_IMPLEMENTATION
#include "../../sfa/format/sfa.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>

static double MU = -41.345, KAPPA = 50.0;

/* Simple BFS cluster detection on |P| field */
static int count_clusters(double *phi[3], int N, double threshold) {
    long N3 = (long)N*N*N;
    int NN = N*N;
    char *visited = (char*)calloc(N3, 1);
    long *queue = (long*)malloc(N3 * sizeof(long));
    int n_clusters = 0;

    for (long idx = 0; idx < N3; idx++) {
        double P = fabs(phi[0][idx] * phi[1][idx] * phi[2][idx]);
        if (P < threshold || visited[idx]) continue;

        /* BFS from this seed */
        n_clusters++;
        int front = 0, back = 0;
        queue[back++] = idx;
        visited[idx] = 1;

        while (front < back) {
            long cur = queue[front++];
            int i = (int)(cur/NN), j = (int)((cur/N)%N), k = (int)(cur%N);

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
    }

    free(visited);
    free(queue);
    return n_clusters;
}

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "Usage: %s input.sfa\n", argv[0]); return 1; }

    SFA *sfa = sfa_open(argv[1]);
    if (!sfa) { fprintf(stderr, "Cannot open %s\n", argv[1]); return 1; }
    int N = sfa->Nx; double L = sfa->Lx;
    long N3 = (long)N*N*N; int NN = N*N;
    double dx = 2.0*L/(N-1), dV = dx*dx*dx;

    /* Read KVMD for physics params */
    SFA_KVMDSet kv[16];
    int n_kv = sfa_read_kvmd(sfa, kv, 16);
    for (int s=0;s<n_kv;s++) for(int p=0;p<kv[s].n_pairs;p++) {
        if(!strcmp(kv[s].keys[p],"mu")) MU=atof(kv[s].values[p]);
        if(!strcmp(kv[s].keys[p],"kappa")) KAPPA=atof(kv[s].values[p]);
    }

    printf("Spatial analysis: %s (N=%d L=%.1f %u frames)\n\n", argv[1], N, L, sfa->total_frames);
    printf("%-6s %8s %8s %8s %7s %7s %7s %7s %9s %8s %5s\n",
           "t", "cx", "cy", "cz", "R_rms", "θ_rms", "θ_var", "P_max", "E_pot", "asp", "clust");

    void *buf = malloc(sfa->frame_bytes);

    for (uint32_t f = 0; f < sfa->total_frames; f++) {
        sfa_read_frame(sfa, f, buf);
        double t = sfa_frame_time(sfa, f);

        /* Extract phi and theta arrays */
        double *phi[3]={NULL}, *theta[3]={NULL};
        uint64_t off = 0;
        for (int c=0; c<sfa->n_columns; c++) {
            int dtype=sfa->columns[c].dtype, sem=sfa->columns[c].semantic, comp=sfa->columns[c].component;
            int es=sfa_dtype_size[dtype]; uint8_t *src=(uint8_t*)buf+off;
            double *arr = (double*)malloc(N3*sizeof(double));
            if(dtype==SFA_F64) memcpy(arr,src,N3*8);
            else if(dtype==SFA_F32) for(long i=0;i<N3;i++) arr[i]=(double)((float*)src)[i];
            else { free(arr); arr=NULL; }
            if(arr) {
                if(sem==SFA_POSITION&&comp<3) phi[comp]=arr;
                else if(sem==SFA_ANGLE&&comp<3) theta[comp]=arr;
                else free(arr);
            }
            off += (uint64_t)N3*es;
        }
        if (!phi[0]) { printf("Frame %u: no phi data\n", f); continue; }

        /* Centroid and R_rms (weighted by |P|) */
        double cm[3]={0}, wt=0, Pmax=0;
        for (long idx=0; idx<N3; idx++) {
            double P = fabs(phi[0][idx]*phi[1][idx]*phi[2][idx]);
            if (P > Pmax) Pmax = P;
            int i=(int)(idx/NN),j=(int)((idx/N)%N),k=(int)(idx%N);
            double x=-L+i*dx, y=-L+j*dx, z=-L+k*dx;
            cm[0]+=x*P; cm[1]+=y*P; cm[2]+=z*P; wt+=P;
        }
        if(wt>0){cm[0]/=wt;cm[1]/=wt;cm[2]/=wt;}

        double R2_sum=0, R_wt=0;
        for (long idx=0; idx<N3; idx++) {
            double P = fabs(phi[0][idx]*phi[1][idx]*phi[2][idx]);
            int i=(int)(idx/NN),j=(int)((idx/N)%N),k=(int)(idx%N);
            double x=-L+i*dx-cm[0], y=-L+j*dx-cm[1], z=-L+k*dx-cm[2];
            R2_sum += (x*x+y*y+z*z)*P; R_wt += P;
        }
        double R_rms = (R_wt>0) ? sqrt(R2_sum/R_wt) : 0;

        /* Theta rms and variation */
        double trms=0, trms2=0;
        if (theta[0]) {
            for (long i=0;i<N3;i++) for(int a=0;a<3;a++) {
                double v = theta[a][i];
                trms += v*v;
                trms2 += v*v*v*v;  /* for variance */
            }
            trms = sqrt(trms/(3*N3));
            trms2 = sqrt(trms2/(3*N3));
        }
        double tvar = trms2 - trms*trms;  /* simplified variance proxy */

        /* E_pot */
        double ep = 0;
        for (long idx=0; idx<N3; idx++) {
            double P = phi[0][idx]*phi[1][idx]*phi[2][idx];
            ep += (MU/2.0)*P*P/(1.0+KAPPA*P*P)*dV;
        }

        /* Aspect ratio */
        double I[3]={0};
        for (long idx=0; idx<N3; idx++) {
            double P = fabs(phi[0][idx]*phi[1][idx]*phi[2][idx]);
            int i=(int)(idx/NN),j=(int)((idx/N)%N),k=(int)(idx%N);
            double x=-L+i*dx-cm[0], y=-L+j*dx-cm[1], z=-L+k*dx-cm[2];
            I[0]+=(y*y+z*z)*P; I[1]+=(x*x+z*z)*P; I[2]+=(x*x+y*y)*P;
        }
        double Imax=I[0],Imin=I[0];
        for(int a=1;a<3;a++){if(I[a]>Imax)Imax=I[a];if(I[a]<Imin)Imin=I[a];}
        double asp = (Imin>1e-30) ? Imax/Imin : 999;

        /* Cluster count */
        int clust = count_clusters(phi, N, 0.01 * Pmax);

        printf("%-6.1f %8.2f %8.2f %8.2f %7.2f %7.4f %7.4f %7.4f %9.1f %8.2f %5d\n",
               t, cm[0], cm[1], cm[2], R_rms, trms, tvar, Pmax, ep, asp, clust);

        for(int a=0;a<3;a++){free(phi[a]); if(theta[a])free(theta[a]);}
    }

    free(buf);
    sfa_close(sfa);
    return 0;
}
