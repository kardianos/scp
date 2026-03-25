/*  shell_analysis.c — Analyze radial shell structure of phi and theta fields
 *
 *  For each frame: compute radial profiles of phi energy, theta energy,
 *  and triple product in 50 thin shells from r=0 to r=L.
 *  Identifies breakaway structures and theta-dominated shells.
 *
 *  Build: gcc -O3 -fopenmp -o shell_analysis shell_analysis.c -I/home/d/code/scp -lzstd -lm
 */

#define SFA_IMPLEMENTATION
#include "sfa/format/sfa.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

static const double MU = -41.345, KAPPA = 50.0, MASS2 = 2.25;

static inline double f16_to_f64(uint16_t h) {
    uint16_t sign = h & 0x8000;
    int exp = (h >> 10) & 0x1F;
    uint16_t mant = h & 0x3FF;
    if (exp == 0) return 0.0;
    if (exp == 31) return sign ? -1e30 : 1e30;
    float fv; uint32_t x = ((uint32_t)sign << 16) | ((uint32_t)(exp-15+127) << 23) | ((uint32_t)mant << 13);
    memcpy(&fv, &x, 4); return (double)fv;
}

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "Usage: %s input.sfa [frame_idx]\n", argv[0]); return 1; }
    int target_frame = -1;
    if (argc > 2) target_frame = atoi(argv[2]);

    SFA *sfa = sfa_open(argv[1]);
    if (!sfa) { fprintf(stderr, "Cannot open %s\n", argv[1]); return 1; }
    int N = sfa->Nx; double L = sfa->Lx;
    long N3 = (long)N*N*N; int NN = N*N;
    double dx = 2.0*L/(N-1);

    printf("# Shell analysis: %s (N=%d L=%.1f %u frames)\n", argv[1], N, L, sfa->total_frames);

    if (target_frame < 0) target_frame = sfa->total_frames - 1;
    if (target_frame >= (int)sfa->total_frames) target_frame = sfa->total_frames - 1;

    void *buf = malloc(sfa->frame_bytes);
    sfa_read_frame(sfa, target_frame, buf);
    double t = sfa_frame_time(sfa, target_frame);

    /* Extract fields */
    double *phi[3] = {NULL}, *theta[3] = {NULL};
    for (int a = 0; a < 3; a++) {
        phi[a] = (double*)calloc(N3, sizeof(double));
        theta[a] = (double*)calloc(N3, sizeof(double));
    }

    uint64_t off = 0;
    for (int c = 0; c < sfa->n_columns; c++) {
        int dtype = sfa->columns[c].dtype;
        int sem = sfa->columns[c].semantic;
        int comp = sfa->columns[c].component;
        int es = sfa_dtype_size[dtype];
        uint8_t *src = (uint8_t*)buf + off;

        double *target = NULL;
        if (sem == SFA_POSITION && comp < 3) target = phi[comp];
        else if (sem == SFA_ANGLE && comp < 3) target = theta[comp];

        if (target) {
            if (dtype == SFA_F64) for(long i=0;i<N3;i++) target[i] = ((double*)src)[i];
            else if (dtype == SFA_F32) for(long i=0;i<N3;i++) target[i] = (double)((float*)src)[i];
            else if (dtype == SFA_F16) for(long i=0;i<N3;i++) target[i] = f16_to_f64(((uint16_t*)src)[i]);
        }
        off += (uint64_t)N3 * es;
    }
    free(buf);

    /* Find centroid weighted by |P| */
    double cm[3] = {0}, wt = 0;
    for (long idx = 0; idx < N3; idx++) {
        double P = fabs(phi[0][idx]*phi[1][idx]*phi[2][idx]);
        int i=(int)(idx/NN), j=(int)((idx/N)%N), k=(int)(idx%N);
        double x=-L+i*dx, y=-L+j*dx, z=-L+k*dx;
        cm[0]+=x*P; cm[1]+=y*P; cm[2]+=z*P; wt+=P;
    }
    if (wt > 0) { cm[0]/=wt; cm[1]/=wt; cm[2]/=wt; }
    printf("# Frame %d, t=%.1f, centroid=(%.2f, %.2f, %.2f)\n\n", target_frame, t, cm[0], cm[1], cm[2]);

    /* Radial shell analysis: 50 shells from r=0 to r=L */
    int NBINS = 50;
    double dr = L / NBINS;

    printf("# r_mid  phi_rms  theta_rms  |P|_avg  E_pot_dens  phi_E  theta_E  theta/phi  n_voxels  phi_max  theta_max\n");

    for (int b = 0; b < NBINS; b++) {
        double r_lo = b * dr, r_hi = (b + 1) * dr;
        double r_mid = (r_lo + r_hi) / 2.0;
        double phi_sum2 = 0, theta_sum2 = 0, P_sum = 0, Ep_sum = 0;
        double phi_E_sum = 0, theta_E_sum = 0;
        double phi_max = 0, theta_max = 0;
        long count = 0;

        for (long idx = 0; idx < N3; idx++) {
            int i=(int)(idx/NN), j=(int)((idx/N)%N), k=(int)(idx%N);
            double x=-L+i*dx-cm[0], y=-L+j*dx-cm[1], z=-L+k*dx-cm[2];
            double r = sqrt(x*x + y*y + z*z);
            if (r < r_lo || r >= r_hi) continue;
            count++;

            double p2 = 0, t2 = 0;
            for (int a = 0; a < 3; a++) {
                p2 += phi[a][idx] * phi[a][idx];
                t2 += theta[a][idx] * theta[a][idx];
                double ap = fabs(phi[a][idx]); if (ap > phi_max) phi_max = ap;
                double at = fabs(theta[a][idx]); if (at > theta_max) theta_max = at;
            }
            phi_sum2 += p2;
            theta_sum2 += t2;
            phi_E_sum += 0.5 * MASS2 * p2;  /* mass energy density proxy */
            theta_E_sum += 0.5 * t2;  /* theta field energy proxy (massless) */

            double P = phi[0][idx]*phi[1][idx]*phi[2][idx];
            double P2 = P*P;
            P_sum += fabs(P);
            Ep_sum += (MU/2.0)*P2/(1.0+KAPPA*P2);
        }

        if (count > 0) {
            double phi_rms = sqrt(phi_sum2 / count);
            double theta_rms = sqrt(theta_sum2 / count);
            double P_avg = P_sum / count;
            double Ep_dens = Ep_sum / count;
            double ratio = (phi_rms > 1e-10) ? theta_rms / phi_rms : 0;

            printf("%6.2f  %8.5f  %9.5f  %8.5f  %10.4f  %8.4f  %8.4f  %10.4f  %8ld  %8.5f  %9.5f\n",
                   r_mid, phi_rms, theta_rms, P_avg, Ep_dens,
                   phi_E_sum/count, theta_E_sum/count, ratio, count, phi_max, theta_max);
        }
    }

    /* Identify breakaway structures: look for local maxima in phi_rms or theta_rms
       at r > 2*R_core that are above background */
    printf("\n# Breakaway structure detection:\n");
    printf("# Looking for local maxima in phi/theta energy beyond the core...\n");

    /* Re-scan for disconnected concentrations */
    double bg_phi = 0, bg_theta = 0;
    long bg_count = 0;
    for (long idx = 0; idx < N3; idx++) {
        int i=(int)(idx/NN), j=(int)((idx/N)%N), k=(int)(idx%N);
        double x=-L+i*dx-cm[0], y=-L+j*dx-cm[1], z=-L+k*dx-cm[2];
        double r = sqrt(x*x + y*y + z*z);
        if (r > L * 0.7) {
            for (int a=0;a<3;a++) { bg_phi += phi[a][idx]*phi[a][idx]; bg_theta += theta[a][idx]*theta[a][idx]; }
            bg_count++;
        }
    }
    if (bg_count > 0) { bg_phi = sqrt(bg_phi/(3*bg_count)); bg_theta = sqrt(bg_theta/(3*bg_count)); }
    printf("# Background levels (r > 0.7*L): phi_rms=%.5f theta_rms=%.5f\n", bg_phi, bg_theta);

    /* Find concentrations above 3× background */
    printf("# Scanning for concentrations > 3x background...\n");
    double threshold_phi = 3.0 * bg_phi;
    double threshold_theta = 3.0 * bg_theta;

    /* Simple grid scan: divide space into 8x8x8 blocks, find those above threshold */
    int BLK = 8;
    int bsize = N / BLK;
    for (int bi=0; bi<BLK; bi++) for (int bj=0; bj<BLK; bj++) for (int bk=0; bk<BLK; bk++) {
        double block_phi2 = 0, block_theta2 = 0;
        long bcnt = 0;
        double bx = 0, by = 0, bz = 0;
        for (int i=bi*bsize; i<(bi+1)*bsize && i<N; i++)
        for (int j=bj*bsize; j<(bj+1)*bsize && j<N; j++)
        for (int k=bk*bsize; k<(bk+1)*bsize && k<N; k++) {
            long idx = (long)i*NN + j*N + k;
            for (int a=0;a<3;a++) { block_phi2+=phi[a][idx]*phi[a][idx]; block_theta2+=theta[a][idx]*theta[a][idx]; }
            bx += -L+i*dx; by += -L+j*dx; bz += -L+k*dx;
            bcnt++;
        }
        if (bcnt == 0) continue;
        double prms = sqrt(block_phi2/(3*bcnt));
        double trms = sqrt(block_theta2/(3*bcnt));
        bx/=bcnt; by/=bcnt; bz/=bcnt;
        double rdist = sqrt((bx-cm[0])*(bx-cm[0])+(by-cm[1])*(by-cm[1])+(bz-cm[2])*(bz-cm[2]));

        if (rdist > 5.0 && (prms > threshold_phi || trms > threshold_theta)) {
            char type[32] = "";
            if (prms > threshold_phi && trms > threshold_theta) strcpy(type, "phi+theta");
            else if (prms > threshold_phi) strcpy(type, "phi-dominated");
            else strcpy(type, "THETA-dominated");
            printf("#   Block at (%.1f,%.1f,%.1f) r=%.1f: phi_rms=%.4f theta_rms=%.4f [%s]\n",
                   bx, by, bz, rdist, prms, trms, type);
        }
    }

    for (int a=0;a<3;a++) { free(phi[a]); free(theta[a]); }
    sfa_close(sfa);
    return 0;
}
