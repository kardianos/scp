/*  gen_pair_seed.c — Place two Lissajous seeds on one grid, separated.
 *
 *  Takes two 64³ seeds and stamps them into an N³ grid at different positions.
 *  Adds analytical background.
 *
 *  Build: gcc -O3 -fopenmp -o gen_pair_seed gen_pair_seed.c -I../../sfa/format -lzstd -lm
 *  Usage: gen_pair_seed seedA.sfa seedB.sfa output.sfa N L separation
 */
#define SFA_IMPLEMENTATION
#include "../../sfa/format/sfa.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define PI 3.14159265358979323846

static void load_seed(const char *path, int *N, double *L, double phi[3][262144], double theta[3][262144], double vel[3][262144]) {
    SFA *s = sfa_open(path);
    if (!s) { fprintf(stderr, "Cannot open %s\n", path); exit(1); }
    *N = s->Nx; *L = s->Lx;
    long N3 = (long)(*N)*(*N)*(*N);
    void *buf = malloc(s->frame_bytes);
    sfa_read_frame(s, 0, buf);

    uint64_t off = 0;
    for (uint32_t col = 0; col < s->n_columns && col < 12; col++) {
        int dtype = s->columns[col].dtype;
        int sem = s->columns[col].semantic;
        int comp = s->columns[col].component;
        int es = sfa_dtype_size[dtype];
        uint8_t *src = (uint8_t*)buf + off;
        double *target = NULL;
        if (sem == 0 && comp < 3) target = phi[comp];
        else if (sem == 1 && comp < 3) target = theta[comp];
        else if (sem == 2 && comp < 3) target = vel[comp];
        if (target) {
            for (long i = 0; i < N3; i++) {
                if (dtype == 1) target[i] = (double)((float*)src)[i];
                else if (dtype == 0) {
                    uint16_t h; memcpy(&h, src+i*2, 2);
                    int exp = (h>>10)&0x1F; uint16_t mant = h&0x3FF; uint16_t sign = h&0x8000;
                    if (exp==0) target[i]=0;
                    else { float f; uint32_t x=((uint32_t)sign<<16)|((uint32_t)(exp-15+127)<<23)|((uint32_t)mant<<13);
                           memcpy(&f,&x,4); target[i]=(double)f; }
                }
            }
        }
        off += (uint64_t)N3 * es;
    }
    free(buf); sfa_close(s);
}

int main(int argc, char **argv) {
    if (argc < 6) {
        fprintf(stderr, "Usage: %s seedA.sfa seedB.sfa output.sfa N L [separation]\n", argv[0]);
        return 1;
    }
    int outN = atoi(argv[4]);
    double outL = atof(argv[5]);
    double sep = argc > 6 ? atof(argv[6]) : 15.0;
    double dx = 2.0 * outL / (outN - 1);
    long outN3 = (long)outN * outN * outN;
    int outNN = outN * outN;

    /* Load seeds */
    int sN; double sL;
    double phiA[3][262144], thetaA[3][262144], velA[3][262144];
    double phiB[3][262144], thetaB[3][262144], velB[3][262144];
    load_seed(argv[1], &sN, &sL, phiA, thetaA, velA);
    printf("Seed A: %s (%d³, L=%.1f)\n", argv[1], sN, sL);
    load_seed(argv[2], &sN, &sL, phiB, thetaB, velB);
    printf("Seed B: %s (%d³, L=%.1f)\n", argv[2], sN, sL);

    int sNN = sN * sN;
    double sDx = 2.0 * sL / (sN - 1);
    double m = 1.5, A_bg = 0.1;
    double k_bg = PI / outL;
    double omega_bg = sqrt(k_bg*k_bg + m*m);
    double delta[3] = {0, 3.0005, 4.4325};

    /* Allocate output */
    double *phi[3], *theta[3], *pvel[3], *tvel[3];
    for (int a = 0; a < 3; a++) {
        phi[a] = (double*)calloc(outN3, sizeof(double));
        theta[a] = (double*)calloc(outN3, sizeof(double));
        pvel[a] = (double*)calloc(outN3, sizeof(double));
        tvel[a] = (double*)calloc(outN3, sizeof(double));
    }

    /* Fill background */
    for (int i = 0; i < outN; i++)
    for (int j = 0; j < outN; j++)
    for (int k = 0; k < outN; k++) {
        long idx = (long)i*outNN + j*outN + k;
        double z = -outL + k*dx;
        for (int a = 0; a < 3; a++) {
            phi[a][idx] = A_bg * cos(k_bg*z + delta[a]);
            pvel[a][idx] = omega_bg * A_bg * sin(k_bg*z + delta[a]);
        }
    }

    /* Stamp seed A at center - sep/2 along x */
    int cxA = outN/2 - (int)(sep/(2*dx));
    int cyA = outN/2, czA = outN/2;
    int half = sN/2;
    printf("Stamping A at grid (%d,%d,%d), sep=%.1f\n", cxA, cyA, czA, sep);

    for (int ti = 0; ti < sN; ti++) {
        int gi = cxA + ti - half; if (gi < 0 || gi >= outN) continue;
        for (int tj = 0; tj < sN; tj++) {
            int gj = cyA + tj - half; if (gj < 0 || gj >= outN) continue;
            for (int tk = 0; tk < sN; tk++) {
                int gk = czA + tk - half; if (gk < 0 || gk >= outN) continue;
                long tidx = (long)ti*sNN + tj*sN + tk;
                long gidx = (long)gi*outNN + gj*outN + gk;
                double z = -sL + tk*sDx;
                for (int a = 0; a < 3; a++) {
                    double bg = A_bg * cos(k_bg*z + delta[a]);
                    double bgv = omega_bg * A_bg * sin(k_bg*z + delta[a]);
                    phi[a][gidx] += phiA[a][tidx] - bg;
                    pvel[a][gidx] += velA[a][tidx] - bgv;
                    theta[a][gidx] += thetaA[a][tidx];
                    tvel[a][gidx] = 0;  /* theta vel from seed isn't meaningful */
                }
            }
        }
    }

    /* Stamp seed B at center + sep/2 along x */
    int cxB = outN/2 + (int)(sep/(2*dx));
    printf("Stamping B at grid (%d,%d,%d)\n", cxB, cyA, czA);

    for (int ti = 0; ti < sN; ti++) {
        int gi = cxB + ti - half; if (gi < 0 || gi >= outN) continue;
        for (int tj = 0; tj < sN; tj++) {
            int gj = cyA + tj - half; if (gj < 0 || gj >= outN) continue;
            for (int tk = 0; tk < sN; tk++) {
                int gk = czA + tk - half; if (gk < 0 || gk >= outN) continue;
                long tidx = (long)ti*sNN + tj*sN + tk;
                long gidx = (long)gi*outNN + gj*outN + gk;
                double z = -sL + tk*sDx;
                for (int a = 0; a < 3; a++) {
                    double bg = A_bg * cos(k_bg*z + delta[a]);
                    double bgv = omega_bg * A_bg * sin(k_bg*z + delta[a]);
                    phi[a][gidx] += phiB[a][tidx] - bg;
                    pvel[a][gidx] += velB[a][tidx] - bgv;
                    theta[a][gidx] += thetaB[a][tidx];
                }
            }
        }
    }

    /* Write output */
    SFA *out = sfa_create(argv[3], outN, outN, outN, outL, outL, outL, 1.0);
    out->flags = 4; /* BSS */
    sfa_add_column(out, "phi_x", 1, 0, 0);
    sfa_add_column(out, "phi_y", 1, 0, 1);
    sfa_add_column(out, "phi_z", 1, 0, 2);
    sfa_add_column(out, "theta_x", 1, 1, 0);
    sfa_add_column(out, "theta_y", 1, 1, 1);
    sfa_add_column(out, "theta_z", 1, 1, 2);
    sfa_add_column(out, "phi_vx", 1, 2, 0);
    sfa_add_column(out, "phi_vy", 1, 2, 1);
    sfa_add_column(out, "phi_vz", 1, 2, 2);
    sfa_add_column(out, "theta_vx", 1, 2, 3);
    sfa_add_column(out, "theta_vy", 1, 2, 4);
    sfa_add_column(out, "theta_vz", 1, 2, 5);
    sfa_finalize_header(out);

    void *cols[12];
    float *cbufs[12];
    for (int c = 0; c < 12; c++) cbufs[c] = (float*)malloc(outN3 * sizeof(float));
    for (long i = 0; i < outN3; i++) {
        for (int a = 0; a < 3; a++) {
            cbufs[a][i] = (float)phi[a][i];
            cbufs[3+a][i] = (float)theta[a][i];
            cbufs[6+a][i] = (float)pvel[a][i];
            cbufs[9+a][i] = (float)tvel[a][i];
        }
    }
    for (int c = 0; c < 12; c++) cols[c] = cbufs[c];
    sfa_write_frame(out, 0.0, cols);
    sfa_close(out);

    printf("Wrote: %s (%d³, L=%.1f, sep=%.1f)\n", argv[3], outN, outL, sep);
    for (int c = 0; c < 12; c++) free(cbufs[c]);
    for (int a = 0; a < 3; a++) { free(phi[a]); free(theta[a]); free(pvel[a]); free(tvel[a]); }
    return 0;
}
