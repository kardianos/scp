/*  gen_lissajous_seed.c — Generate a 3-axis superposition seed
 *
 *  Takes a base blob seed (z-oriented carrier) and creates a composite by
 *  superimposing three rotated copies (carriers along x, y, z) with
 *  independent phase shifts.
 *
 *  The rotation maps the z-carrier to each axis by permuting spatial coords:
 *    z-blob: φ_a(x,y,z) as-is
 *    x-blob: φ_a(z,y,x) — swap x↔z so carrier runs along x
 *    y-blob: φ_a(x,z,y) — swap y↔z so carrier runs along y
 *
 *  Phase shifts are applied to the carrier: cos(k*axis + δ_a + Δ)
 *  by shifting the blob's position along its carrier axis.
 *
 *  Build: gcc -O3 -march=native -fopenmp -o gen_lissajous_seed gen_lissajous_seed.c -lzstd -lm
 *  Usage: gen_lissajous_seed base.sfa output.sfa Δ0 Δ1 Δ2 [--theta_gain G]
 */
#define SFA_IMPLEMENTATION
#include "../../sfa/format/sfa.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define PI 3.14159265358979323846

int main(int argc, char **argv) {
    if (argc < 6) {
        fprintf(stderr, "Usage: %s base.sfa output.sfa Δ0 Δ1 Δ2 [--theta_gain G]\n", argv[0]);
        fprintf(stderr, "  Δ0,Δ1,Δ2: phase shifts in radians for x,y,z blobs\n");
        fprintf(stderr, "  --theta_gain G: theta = -G * curl(phi) (default 2.5)\n");
        return 1;
    }

    const char *inpath = argv[1];
    const char *outpath = argv[2];
    double delta_phase[3] = { atof(argv[3]), atof(argv[4]), atof(argv[5]) };
    double theta_gain = -2.5;  /* negative = anti-aligned, our best chirality */

    for (int i = 6; i < argc; i++) {
        if (!strcmp(argv[i], "--theta_gain") && i+1 < argc) theta_gain = atof(argv[++i]);
    }

    /* Open base seed */
    SFA *in = sfa_open(inpath);
    if (!in) { fprintf(stderr, "Cannot open %s\n", inpath); return 1; }
    int N = in->Nx;
    long N3 = (long)N*N*N;
    int NN = N*N;
    double L = in->Lx;
    double dx = 2.0 * L / (N - 1);

    /* Read base seed fields */
    void *buf = malloc(in->frame_bytes);
    sfa_read_frame(in, 0, buf);

    double *base_phi[3], *base_vel[3];
    for (int a = 0; a < 3; a++) {
        base_phi[a] = (double*)calloc(N3, sizeof(double));
        base_vel[a] = (double*)calloc(N3, sizeof(double));
    }

    uint64_t off = 0;
    for (uint32_t col = 0; col < in->n_columns; col++) {
        int dtype = in->columns[col].dtype;
        int sem = in->columns[col].semantic;
        int comp = in->columns[col].component;
        int es = sfa_dtype_size[dtype];
        uint8_t *src = (uint8_t*)buf + off;
        double *target = NULL;
        if (sem == SFA_POSITION && comp < 3) target = base_phi[comp];
        else if (sem == SFA_VELOCITY && comp < 3) target = base_vel[comp];
        if (target) {
            for (long i = 0; i < N3; i++) {
                if (dtype == SFA_F32) target[i] = (double)((float*)src)[i];
                else if (dtype == SFA_F16) {
                    uint16_t h; memcpy(&h, src + i*2, 2);
                    uint16_t sign = h & 0x8000;
                    int exp = (h >> 10) & 0x1F;
                    uint16_t mant = h & 0x3FF;
                    if (exp == 0) target[i] = 0;
                    else { float f; uint32_t x = ((uint32_t)sign<<16)|((uint32_t)(exp-15+127)<<23)|((uint32_t)mant<<13);
                           memcpy(&f, &x, 4); target[i] = (double)f; }
                } else target[i] = ((double*)src)[i];
            }
        }
        off += (uint64_t)N3 * es;
    }
    free(buf); sfa_close(in);

    /* The base seed has a z-oriented carrier. To create x and y versions,
     * we permute the spatial indices:
     *   z-blob at (i,j,k): field value from base at (i,j,k)
     *   x-blob at (i,j,k): field value from base at (k,j,i)  [swap x↔z]
     *   y-blob at (i,j,k): field value from base at (i,k,j)  [swap y↔z]
     *
     * Phase shift Δ is applied by shifting the carrier axis coordinate:
     *   For z-blob: shift k by round(Δ/(k_bg * dx))
     *   For x-blob: shift i by the same amount
     *   For y-blob: shift j by the same amount
     *
     * k_bg = PI/L, so one full wavelength = 2*L grid points = N-1.
     * Phase shift Δ corresponds to shift of Δ/(2π) * (N-1) grid points.
     */
    double k_bg = PI / L;
    int shift[3];  /* grid-point shift for each axis */
    for (int ax = 0; ax < 3; ax++) {
        shift[ax] = (int)round(delta_phase[ax] / (2.0 * PI) * (N - 1));
        /* Wrap to [0, N) */
        shift[ax] = ((shift[ax] % N) + N) % N;
    }

    printf("Lissajous seed: Δ=(%.4f, %.4f, %.4f) rad, shifts=(%d, %d, %d) voxels\n",
           delta_phase[0], delta_phase[1], delta_phase[2], shift[0], shift[1], shift[2]);

    /* Superimpose three rotated+shifted copies */
    double *phi[3], *phi_vel[3];
    for (int a = 0; a < 3; a++) {
        phi[a] = (double*)calloc(N3, sizeof(double));
        phi_vel[a] = (double*)calloc(N3, sizeof(double));
    }

    for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
    for (int k = 0; k < N; k++) {
        long idx = (long)i*NN + j*N + k;

        /* z-blob: base at (i, j, (k+shift_z)%N) */
        {
            int ks = (k + shift[2]) % N;
            long src_idx = (long)i*NN + j*N + ks;
            for (int a = 0; a < 3; a++) {
                phi[a][idx]     += base_phi[a][src_idx];
                phi_vel[a][idx] += base_vel[a][src_idx];
            }
        }

        /* x-blob: base at ((k+shift_x)%N, j, i) — swap x↔z in source lookup */
        {
            int is = (k + shift[0]) % N;  /* source z = destination x, shifted */
            long src_idx = (long)is*NN + j*N + i;  /* source (z',j,x) where x=dest_k→src_i */
            /* Wait — need to be careful. The base carrier runs along z (index k).
             * For x-blob, we want the carrier to run along x (index i).
             * So at output voxel (i,j,k), the x-blob's value comes from base at
             * a position where the carrier coord = i (shifted).
             * Base is indexed as (i_src, j_src, k_src) with carrier along k_src.
             * So k_src maps to our i, and i_src maps to our k: */
            int k_src = (i + shift[0]) % N;  /* carrier axis = i, shifted */
            long si = (long)k*NN + j*N + k_src;  /* (k→i_src, j→j_src, k_src→k_src) */
            for (int a = 0; a < 3; a++) {
                phi[a][idx]     += base_phi[a][si];
                phi_vel[a][idx] += base_vel[a][si];
            }
        }

        /* y-blob: carrier along y (index j) */
        {
            int j_src = (j + shift[1]) % N;  /* carrier axis = j, shifted */
            long si = (long)i*NN + k*N + j_src;  /* (i→i_src, k→j_src, j_src→k_src) */
            for (int a = 0; a < 3; a++) {
                phi[a][idx]     += base_phi[a][si];
                phi_vel[a][idx] += base_vel[a][si];
            }
        }
    }

    /* Compute curl(phi) and set theta = theta_gain * curl(phi) */
    double *theta[3], *theta_vel[3];
    for (int a = 0; a < 3; a++) {
        theta[a] = (double*)calloc(N3, sizeof(double));
        theta_vel[a] = (double*)calloc(N3, sizeof(double));
    }

    double idx2 = 1.0 / (2.0 * dx);
    #pragma omp parallel for schedule(static)
    for (long vi = 0; vi < N3; vi++) {
        int i = (int)(vi / NN), j = (int)((vi / N) % N), k = (int)(vi % N);
        int ip = (i+1)%N, im = (i-1+N)%N;
        int jp = (j+1)%N, jm = (j-1+N)%N;
        int kp = (k+1)%N, km = (k-1+N)%N;
        double cx = (phi[2][(long)i*NN+jp*N+k] - phi[2][(long)i*NN+jm*N+k]) * idx2
                  - (phi[1][(long)i*NN+j*N+kp] - phi[1][(long)i*NN+j*N+km]) * idx2;
        double cy = (phi[0][(long)i*NN+j*N+kp] - phi[0][(long)i*NN+j*N+km]) * idx2
                  - (phi[2][(long)ip*NN+j*N+k] - phi[2][(long)im*NN+j*N+k]) * idx2;
        double cz = (phi[1][(long)ip*NN+j*N+k] - phi[1][(long)im*NN+j*N+k]) * idx2
                  - (phi[0][(long)i*NN+jp*N+k] - phi[0][(long)i*NN+jm*N+k]) * idx2;
        theta[0][vi] = theta_gain * cx;
        theta[1][vi] = theta_gain * cy;
        theta[2][vi] = theta_gain * cz;
    }

    /* Stats */
    double phi_max = 0, P_int = 0, theta_rms = 0;
    for (long i = 0; i < N3; i++) {
        double ps = phi[0][i]*phi[0][i] + phi[1][i]*phi[1][i] + phi[2][i]*phi[2][i];
        if (sqrt(ps) > phi_max) phi_max = sqrt(ps);
        P_int += fabs(phi[0][i] * phi[1][i] * phi[2][i]);
        theta_rms += theta[0][i]*theta[0][i] + theta[1][i]*theta[1][i] + theta[2][i]*theta[2][i];
    }
    double dx3 = dx*dx*dx;
    P_int *= dx3;
    theta_rms = sqrt(theta_rms / N3);
    printf("  phi_max=%.4f P_int=%.4f theta_rms=%.6f\n", phi_max, P_int, theta_rms);

    /* Write output */
    SFA *out = sfa_create(outpath, N, N, N, L, L, L, 1.0);
    out->flags = SFA_CODEC_BSS;
    sfa_add_column(out, "phi_x",    SFA_F32, SFA_POSITION, 0);
    sfa_add_column(out, "phi_y",    SFA_F32, SFA_POSITION, 1);
    sfa_add_column(out, "phi_z",    SFA_F32, SFA_POSITION, 2);
    sfa_add_column(out, "theta_x",  SFA_F32, SFA_ANGLE,    0);
    sfa_add_column(out, "theta_y",  SFA_F32, SFA_ANGLE,    1);
    sfa_add_column(out, "theta_z",  SFA_F32, SFA_ANGLE,    2);
    sfa_add_column(out, "phi_vx",   SFA_F32, SFA_VELOCITY, 0);
    sfa_add_column(out, "phi_vy",   SFA_F32, SFA_VELOCITY, 1);
    sfa_add_column(out, "phi_vz",   SFA_F32, SFA_VELOCITY, 2);
    sfa_add_column(out, "theta_vx", SFA_F32, SFA_VELOCITY, 3);
    sfa_add_column(out, "theta_vy", SFA_F32, SFA_VELOCITY, 4);
    sfa_add_column(out, "theta_vz", SFA_F32, SFA_VELOCITY, 5);
    sfa_finalize_header(out);

    void *cols[12];
    float *cbufs[12];
    for (int c = 0; c < 12; c++) cbufs[c] = (float*)malloc(N3 * sizeof(float));
    for (long i = 0; i < N3; i++) {
        for (int a = 0; a < 3; a++) {
            cbufs[a][i]   = (float)phi[a][i];
            cbufs[3+a][i] = (float)theta[a][i];
            cbufs[6+a][i] = (float)phi_vel[a][i];
            cbufs[9+a][i] = (float)theta_vel[a][i];
        }
    }
    for (int c = 0; c < 12; c++) cols[c] = cbufs[c];
    sfa_write_frame(out, 0.0, cols);
    sfa_close(out);

    printf("Wrote: %s\n", outpath);

    for (int a = 0; a < 3; a++) {
        free(phi[a]); free(phi_vel[a]); free(theta[a]); free(theta_vel[a]);
        free(base_phi[a]); free(base_vel[a]);
    }
    for (int c = 0; c < 12; c++) free(cbufs[c]);
    return 0;
}
