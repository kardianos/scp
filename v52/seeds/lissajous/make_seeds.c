/*  make_seeds.c — Generate 5 chirality variants from a base blob seed.
 *
 *  Takes the base seed (phi blob with whatever theta it has) and creates:
 *    seed_theta_zero.sfa      — theta = 0
 *    seed_theta_pos_mod.sfa   — theta = +G * curl(phi), G moderate
 *    seed_theta_pos_ext.sfa   — theta = +G * curl(phi), G extreme
 *    seed_theta_neg_mod.sfa   — theta = -G * curl(phi), G moderate
 *    seed_theta_neg_ext.sfa   — theta = -G * curl(phi), G extreme
 *
 *  curl(phi) is the natural theta source in the Cosserat equation.
 *  Sign of theta relative to curl(phi) = chirality.
 *
 *  Build: gcc -O3 -march=native -fopenmp -o make_seeds make_seeds.c -lzstd -lm
 */
#define SFA_IMPLEMENTATION
#include "../../sfa/format/sfa.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

static void compute_curl(double **phi, int N, double dx, double **curl_out) {
    int NN = N * N;
    double idx = 1.0 / (2.0 * dx);
    #pragma omp parallel for schedule(static)
    for (long idx_v = 0; idx_v < (long)N*N*N; idx_v++) {
        int i = (int)(idx_v / NN), j = (int)((idx_v / N) % N), k = (int)(idx_v % N);
        int ip = (i+1)%N, im = (i-1+N)%N;
        int jp = (j+1)%N, jm = (j-1+N)%N;
        int kp = (k+1)%N, km = (k-1+N)%N;
        /* curl(phi)_x = dphi_z/dy - dphi_y/dz */
        curl_out[0][idx_v] = (phi[2][(long)i*NN+jp*N+k] - phi[2][(long)i*NN+jm*N+k]) * idx
                           - (phi[1][(long)i*NN+j*N+kp] - phi[1][(long)i*NN+j*N+km]) * idx;
        /* curl(phi)_y = dphi_x/dz - dphi_z/dx */
        curl_out[1][idx_v] = (phi[0][(long)i*NN+j*N+kp] - phi[0][(long)i*NN+j*N+km]) * idx
                           - (phi[2][(long)ip*NN+j*N+k] - phi[2][(long)im*NN+j*N+k]) * idx;
        /* curl(phi)_z = dphi_y/dx - dphi_x/dy */
        curl_out[2][idx_v] = (phi[1][(long)ip*NN+j*N+k] - phi[1][(long)im*NN+j*N+k]) * idx
                           - (phi[0][(long)i*NN+jp*N+k] - phi[0][(long)i*NN+jm*N+k]) * idx;
    }
}

static void write_seed(const char *path, int N, double L, double dt,
                       double **phi, double **theta, double **phi_vel, double **theta_vel,
                       int dtype) {
    SFA *out = sfa_create(path, N, N, N, L, L, L, dt);
    out->flags = SFA_CODEC_BSS;
    uint8_t sd = (dtype == 1) ? SFA_F32 : SFA_F16;
    sfa_add_column(out, "phi_x",    sd, SFA_POSITION, 0);
    sfa_add_column(out, "phi_y",    sd, SFA_POSITION, 1);
    sfa_add_column(out, "phi_z",    sd, SFA_POSITION, 2);
    sfa_add_column(out, "theta_x",  sd, SFA_ANGLE,    0);
    sfa_add_column(out, "theta_y",  sd, SFA_ANGLE,    1);
    sfa_add_column(out, "theta_z",  sd, SFA_ANGLE,    2);
    sfa_add_column(out, "phi_vx",   sd, SFA_VELOCITY, 0);
    sfa_add_column(out, "phi_vy",   sd, SFA_VELOCITY, 1);
    sfa_add_column(out, "phi_vz",   sd, SFA_VELOCITY, 2);
    sfa_add_column(out, "theta_vx", sd, SFA_VELOCITY, 3);
    sfa_add_column(out, "theta_vy", sd, SFA_VELOCITY, 4);
    sfa_add_column(out, "theta_vz", sd, SFA_VELOCITY, 5);
    sfa_finalize_header(out);

    long N3 = (long)N*N*N;
    void *cols[12];
    float *bufs[12];
    for (int c = 0; c < 12; c++) bufs[c] = (float*)malloc(N3 * sizeof(float));
    for (long i = 0; i < N3; i++) {
        for (int a = 0; a < 3; a++) {
            bufs[a][i]   = (float)phi[a][i];
            bufs[3+a][i] = (float)theta[a][i];
            bufs[6+a][i] = (float)phi_vel[a][i];
            bufs[9+a][i] = (float)theta_vel[a][i];
        }
    }
    for (int c = 0; c < 12; c++) cols[c] = bufs[c];
    sfa_write_frame(out, 0.0, cols);
    sfa_close(out);
    for (int c = 0; c < 12; c++) free(bufs[c]);
}

int main(int argc, char **argv) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s base_seed.sfa [outdir]\n", argv[0]);
        return 1;
    }
    const char *inpath = argv[1];
    const char *outdir = argc > 2 ? argv[2] : ".";

    /* Open base seed */
    SFA *in = sfa_open(inpath);
    if (!in) { fprintf(stderr, "Cannot open %s\n", inpath); return 1; }
    int N = in->Nx;
    long N3 = (long)N*N*N;
    double L = in->Lx;
    double dx = 2.0 * L / (N - 1);
    double dt = in->dt;
    printf("Base seed: %s (%d^3, L=%.2f)\n", inpath, N, L);

    /* Read frame */
    void *buf = malloc(in->frame_bytes);
    sfa_read_frame(in, 0, buf);

    /* Extract fields */
    double *phi[3], *theta_base[3], *phi_vel[3], *theta_vel[3];
    for (int a = 0; a < 3; a++) {
        phi[a]        = (double*)calloc(N3, sizeof(double));
        theta_base[a] = (double*)calloc(N3, sizeof(double));
        phi_vel[a]    = (double*)calloc(N3, sizeof(double));
        theta_vel[a]  = (double*)calloc(N3, sizeof(double));
    }

    uint64_t off = 0;
    for (uint32_t col = 0; col < in->n_columns; col++) {
        int dtype = in->columns[col].dtype;
        int sem = in->columns[col].semantic;
        int comp = in->columns[col].component;
        int es = sfa_dtype_size[dtype];
        uint8_t *src = (uint8_t*)buf + off;
        double *target = NULL;
        if (sem == SFA_POSITION && comp < 3) target = phi[comp];
        else if (sem == SFA_ANGLE && comp < 3) target = theta_base[comp];
        else if (sem == SFA_VELOCITY && comp < 3) target = phi_vel[comp];
        else if (sem == SFA_VELOCITY && comp >= 3 && comp < 6) target = theta_vel[comp-3];
        if (target) {
            for (long i = 0; i < N3; i++) {
                if (dtype == SFA_F32) target[i] = (double)((float*)src)[i];
                else if (dtype == SFA_F16) {
                    uint16_t h; memcpy(&h, src + i*2, 2);
                    uint16_t sign = h & 0x8000;
                    int exp = (h >> 10) & 0x1F;
                    uint16_t mant = h & 0x3FF;
                    if (exp == 0) { target[i] = sign ? -0.0 : 0.0; }
                    else {
                        float f; uint32_t x = ((uint32_t)sign<<16)|((uint32_t)(exp-15+127)<<23)|((uint32_t)mant<<13);
                        memcpy(&f, &x, 4); target[i] = (double)f;
                    }
                } else target[i] = ((double*)src)[i];
            }
        }
        off += (uint64_t)N3 * es;
    }
    free(buf); sfa_close(in);

    /* Compute curl(phi) */
    double *curl_phi[3];
    for (int a = 0; a < 3; a++) curl_phi[a] = (double*)calloc(N3, sizeof(double));
    compute_curl(phi, N, dx, curl_phi);

    /* Stats */
    double curl_rms = 0, phi_max = 0, theta_base_rms = 0;
    for (long i = 0; i < N3; i++) {
        double cs = curl_phi[0][i]*curl_phi[0][i] + curl_phi[1][i]*curl_phi[1][i] + curl_phi[2][i]*curl_phi[2][i];
        curl_rms += cs;
        double ps = phi[0][i]*phi[0][i] + phi[1][i]*phi[1][i] + phi[2][i]*phi[2][i];
        if (sqrt(ps) > phi_max) phi_max = sqrt(ps);
        double ts = theta_base[0][i]*theta_base[0][i] + theta_base[1][i]*theta_base[1][i] + theta_base[2][i]*theta_base[2][i];
        theta_base_rms += ts;
    }
    curl_rms = sqrt(curl_rms / N3);
    theta_base_rms = sqrt(theta_base_rms / N3);
    printf("  phi_max = %.4f\n", phi_max);
    printf("  curl(phi)_rms = %.6f\n", curl_rms);
    printf("  theta_base_rms = %.6f\n", theta_base_rms);

    /* Define gain factors */
    /* Moderate: G such that theta_rms ~ curl_rms * eta (natural coupling scale) */
    /* Extreme: 5x moderate */
    double G_mod = 0.5;   /* theta = 0.5 * curl(phi) — half the natural Cosserat coupling */
    double G_ext = 2.5;   /* theta = 2.5 * curl(phi) — strong chirality */
    printf("\nGain factors: moderate=%.2f, extreme=%.2f\n", G_mod, G_ext);

    /* Generate 5 seeds */
    double *theta[3], *tvel_zero[3];
    for (int a = 0; a < 3; a++) {
        theta[a] = (double*)calloc(N3, sizeof(double));
        tvel_zero[a] = (double*)calloc(N3, sizeof(double));
    }

    struct { const char *name; double gain; } seeds[] = {
        {"seed_theta_zero",    0.0},
        {"seed_theta_pos_mod", G_mod},
        {"seed_theta_pos_ext", G_ext},
        {"seed_theta_neg_mod", -G_mod},
        {"seed_theta_neg_ext", -G_ext},
    };

    for (int s = 0; s < 5; s++) {
        double G = seeds[s].gain;
        char path[512];
        snprintf(path, sizeof(path), "%s/%s.sfa", outdir, seeds[s].name);

        /* Set theta = G * curl(phi) */
        double trms = 0;
        for (long i = 0; i < N3; i++) {
            for (int a = 0; a < 3; a++)
                theta[a][i] = G * curl_phi[a][i];
            trms += theta[0][i]*theta[0][i] + theta[1][i]*theta[1][i] + theta[2][i]*theta[2][i];
        }
        trms = sqrt(trms / N3);

        write_seed(path, N, L, dt, phi, theta, phi_vel, tvel_zero, 1);
        printf("  %s: G=%.2f theta_rms=%.6f\n", seeds[s].name, G, trms);
    }

    /* Cleanup */
    for (int a = 0; a < 3; a++) {
        free(phi[a]); free(theta_base[a]); free(phi_vel[a]); free(theta_vel[a]);
        free(curl_phi[a]); free(theta[a]); free(tvel_zero[a]);
    }
    printf("\nDone. 5 seeds in %s/\n", outdir);
    return 0;
}
