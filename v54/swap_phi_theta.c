/*  swap_phi_theta.c — Swap phi and theta channels in an SFA seed
 *
 *  Takes a phi-dominant seed and creates a theta-dominant one.
 *  phi_out = theta_in (weak, curl-derived)
 *  theta_out = phi_in (strong, structured)
 *
 *  Build: gcc -O3 -o swap_phi_theta swap_phi_theta.c -lzstd -lm
 *  Usage: swap_phi_theta input.sfa output.sfa [scale]
 *    scale: multiply the swapped fields (default 1.0)
 */
#define SFA_IMPLEMENTATION
#include "../sfa/format/sfa.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

int main(int argc, char **argv) {
    if (argc < 3) {
        fprintf(stderr, "Usage: %s input.sfa output.sfa [scale]\n", argv[0]);
        return 1;
    }
    const char *inpath = argv[1];
    const char *outpath = argv[2];
    double scale = argc > 3 ? atof(argv[3]) : 1.0;

    SFA *in = sfa_open(inpath);
    if (!in) { fprintf(stderr, "Cannot open %s\n", inpath); return 1; }
    int N = in->Nx;
    long N3 = (long)N*N*N;
    double L = in->Lx;

    void *buf = malloc(in->frame_bytes);
    sfa_read_frame(in, 0, buf);

    /* Parse columns: find phi, theta, phi_vel, theta_vel */
    float *phi[3]={0}, *theta[3]={0}, *phi_v[3]={0}, *theta_v[3]={0};
    uint64_t off = 0;
    for (uint32_t col = 0; col < in->n_columns; col++) {
        int dtype = in->columns[col].dtype;
        int sem = in->columns[col].semantic;
        int comp = in->columns[col].component;
        int es = sfa_dtype_size[dtype];
        float *ptr = NULL;
        if (dtype == SFA_F32) ptr = (float*)((uint8_t*)buf + off);

        if (ptr) {
            if (sem == SFA_POSITION && comp < 3) phi[comp] = ptr;
            else if (sem == SFA_ANGLE && comp < 3) theta[comp] = ptr;
            else if (sem == SFA_VELOCITY && comp < 3) phi_v[comp] = ptr;
            else if (sem == SFA_VELOCITY && comp >= 3 && comp < 6) theta_v[comp-3] = ptr;
        }
        off += (uint64_t)N3 * es;
    }

    if (!phi[0] || !theta[0]) {
        fprintf(stderr, "Could not find phi and theta columns\n");
        return 1;
    }

    /* Swap: new_phi = old_theta * scale, new_theta = old_phi * scale */
    float *new_phi[3], *new_theta[3], *new_pv[3], *new_tv[3];
    for (int a = 0; a < 3; a++) {
        new_phi[a] = malloc(N3 * sizeof(float));
        new_theta[a] = malloc(N3 * sizeof(float));
        new_pv[a] = calloc(N3, sizeof(float));
        new_tv[a] = calloc(N3, sizeof(float));
        for (long i = 0; i < N3; i++) {
            new_phi[a][i] = theta[a] ? scale * theta[a][i] : 0;
            new_theta[a][i] = scale * phi[a][i];
        }
        if (phi_v[a] && theta_v[a]) {
            for (long i = 0; i < N3; i++) {
                new_pv[a][i] = scale * theta_v[a][i];
                new_tv[a][i] = scale * phi_v[a][i];
            }
        }
    }

    /* Stats */
    double phi_max = 0, theta_max = 0;
    for (long i = 0; i < N3; i++) {
        for (int a = 0; a < 3; a++) {
            double v = fabs(new_phi[a][i]); if (v > phi_max) phi_max = v;
            v = fabs(new_theta[a][i]); if (v > theta_max) theta_max = v;
        }
    }
    printf("Swapped: phi_max=%.4f theta_max=%.4f (scale=%.2f)\n",
           phi_max, theta_max, scale);

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

    void *cols[12] = {
        new_phi[0], new_phi[1], new_phi[2],
        new_theta[0], new_theta[1], new_theta[2],
        new_pv[0], new_pv[1], new_pv[2],
        new_tv[0], new_tv[1], new_tv[2]
    };
    sfa_write_frame(out, 0.0, cols);
    sfa_close(out);
    printf("Wrote: %s\n", outpath);

    free(buf); sfa_close(in);
    for (int a = 0; a < 3; a++) {
        free(new_phi[a]); free(new_theta[a]);
        free(new_pv[a]); free(new_tv[a]);
    }
    return 0;
}
