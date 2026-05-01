/*  gen_temporal_seed_v2.c — Extract temporal mean perturbation + add background back
 *
 *  The temporal mean of polynomial coefficients captures the time-averaged field.
 *  But the background carrier wave cos(k*z+δ)*cos(ωt) has zero time average.
 *  So: temporal mean ≈ soliton perturbation only (E_total ≈ 4, useless).
 *
 *  Fix: evaluate temporal mean (= soliton perturbation), then ADD the analytical
 *  background back. Set velocities to the background velocity field.
 *
 *  Build: gcc -O3 -march=native -fopenmp -o gen_temporal_seed_v2 gen_temporal_seed_v2.c -lzstd -lm
 *  Usage: gen_temporal_seed_v2 input.sfa output.sfa
 */
#define SFA_IMPLEMENTATION
#include "../format/sfa.h"
#include "../../sfa/sim/scp_config.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

int main(int argc, char **argv) {
    if (argc < 3) {
        fprintf(stderr, "Usage: %s input.sfa output.sfa\n", argv[0]);
        return 1;
    }
    const char *inpath = argv[1], *outpath = argv[2];

    /* Open input */
    sfa_fixup_index(inpath);
    SFA *in = sfa_open(inpath);
    if (!in) { fprintf(stderr, "Cannot open %s\n", inpath); return 1; }

    int N = in->Nx;
    long N3 = (long)N * N * N;
    int NN = N * N;
    double L = in->Lx;
    double dx = 2.0 * L / (N - 1);
    printf("Input: %s (%d^3, L=%.2f, %d frames)\n", inpath, N, L, in->total_frames);

    /* Read KVMD to get physics parameters */
    double A_bg = 0.1, m = 1.5;
    double delta[3] = {0, 3.0005, 4.4325};
    SFA_KVMDSet kv[SFA_MAX_KVMD_SETS];
    int nkv = sfa_read_kvmd(in, kv, SFA_MAX_KVMD_SETS);
    for (int s = 0; s < nkv; s++) {
        for (int p = 0; p < kv[s].n_pairs; p++) {
            if (!strcmp(kv[s].keys[p], "A_bg")) A_bg = atof(kv[s].values[p]);
            if (!strcmp(kv[s].keys[p], "m")) m = atof(kv[s].values[p]);
            if (!strcmp(kv[s].keys[p], "delta")) sscanf(kv[s].values[p], "%lf,%lf,%lf", &delta[0], &delta[1], &delta[2]);
        }
    }
    double m2 = m * m;
    double k_bg = PI / L;
    double omega_bg = sqrt(k_bg * k_bg + m2);
    printf("Background: A_bg=%.3f, m=%.3f, k_bg=%.4f, omega_bg=%.4f\n", A_bg, m, k_bg, omega_bg);
    printf("Delta: %.4f, %.4f, %.4f\n", delta[0], delta[1], delta[2]);

    /* Find last I-frame with temporal model */
    int found_iframe = -1;
    for (int fi = in->total_frames - 1; fi >= 0; fi--) {
        uint32_t ft = sfa_frame_type(in, fi);
        if (ft == SFA_FRAME_VEC_I) {
            found_iframe = fi;
            break;
        }
    }
    if (found_iframe < 0) {
        fprintf(stderr, "No VEC_I frame found\n");
        return 1;
    }
    printf("Last I-frame: %d\n", found_iframe);

    /* Read the I-frame to populate temporal model */
    void *buf = malloc(in->frame_bytes);
    if (sfa_read_frame(in, found_iframe, buf) < 0) {
        fprintf(stderr, "Cannot read frame %d\n", found_iframe);
        return 1;
    }

    if (!in->vec_temporal_valid) {
        fprintf(stderr, "No temporal model in I-frame %d\n", found_iframe);
        return 1;
    }

    uint32_t np = in->vec_n_patches;
    uint16_t nc = in->vec_n_coeffs;
    long n_total = (long)np * nc;
    int n_fields = nc / 64;
    int pf_nc = 64;
    int bs = 8;
    int BN = N / bs;

    printf("Temporal model: %d patches, %d coeffs (%d fields x %d), omega=%.4f\n",
           np, nc, n_fields, pf_nc, in->vec_temporal_omega);

    /* Allocate 6 field arrays */
    double *phi[3], *theta[3], *phi_vel[3], *theta_vel[3];
    for (int a = 0; a < 3; a++) {
        phi[a] = (double*)calloc(N3, sizeof(double));
        theta[a] = (double*)calloc(N3, sizeof(double));
        phi_vel[a] = (double*)calloc(N3, sizeof(double));
        theta_vel[a] = (double*)calloc(N3, sizeof(double));
    }

    /* Evaluate temporal MEAN coefficients to voxel grid */
    float *mean = in->vec_temporal_mean;
    printf("Evaluating temporal mean to voxels...\n");

    #pragma omp parallel for schedule(dynamic)
    for (uint32_t pi = 0; pi < np; pi++) {
        int bk = pi % BN, bj = (pi / BN) % BN, bi = pi / (BN * BN);
        int ox = bi * bs, oy = bj * bs, oz = bk * bs;

        for (int di = 0; di < bs; di++) {
            int gi = ox + di; if (gi >= N) continue;
            double tx = (double)di / (bs - 1);
            double txa[4] = {1, tx, tx*tx, tx*tx*tx};
            for (int dj = 0; dj < bs; dj++) {
                int gj = oy + dj; if (gj >= N) continue;
                double ty = (double)dj / (bs - 1);
                double tya[4] = {1, ty, ty*ty, ty*ty*ty};
                for (int dk = 0; dk < bs; dk++) {
                    int gk = oz + dk; if (gk >= N) continue;
                    double tz = (double)dk / (bs - 1);
                    double tza[4] = {1, tz, tz*tz, tz*tz*tz};

                    long idx = (long)gi * NN + gj * N + gk;

                    for (int f = 0; f < n_fields && f < 6; f++) {
                        float *fc = mean + (long)pi * nc + f * pf_nc;
                        double val = 0;
                        for (int a = 0; a < 4; a++)
                        for (int b = 0; b < 4; b++) {
                            double ab = txa[a] * tya[b];
                            for (int c = 0; c < 4; c++)
                                val += fc[a*16 + b*4 + c] * ab * tza[c];
                        }
                        if (f < 3) phi[f][idx] = val;
                        else       theta[f-3][idx] = val;
                    }
                }
            }
        }
    }

    /* Report perturbation stats before adding background */
    for (int a = 0; a < 3; a++) {
        double mx = 0, mn = 0, sum = 0;
        for (long i = 0; i < N3; i++) {
            if (phi[a][i] > mx) mx = phi[a][i];
            if (phi[a][i] < mn) mn = phi[a][i];
            sum += phi[a][i];
        }
        printf("  perturbation phi_%d: min=%.6f max=%.6f mean=%.6e\n", a, mn, mx, sum/N3);
    }

    /* ADD analytical background */
    printf("Adding analytical background (A_bg=%.3f)...\n", A_bg);
    double E_kin_bg = 0, E_pot_bg = 0;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            for (int k = 0; k < N; k++) {
                long idx = (long)i * NN + j * N + k;
                double z = -L + k * dx;
                for (int a = 0; a < 3; a++) {
                    double ph_bg = k_bg * z + delta[a];
                    double bg_phi = A_bg * cos(ph_bg);
                    double bg_vel = omega_bg * A_bg * sin(ph_bg);
                    phi[a][idx] += bg_phi;
                    phi_vel[a][idx] = bg_vel;
                    E_kin_bg += 0.5 * bg_vel * bg_vel;
                    E_pot_bg += 0.5 * m2 * bg_phi * bg_phi;
                }
                /* theta velocities: zero (background has no theta) */
            }
        }
    }

    /* Compute seed diagnostics */
    double phi_max = 0, P_int = 0, E_total = 0;
    double E_kin = 0, E_mass = 0;
    for (long i = 0; i < N3; i++) {
        double p0 = phi[0][i], p1 = phi[1][i], p2 = phi[2][i];
        double P = fabs(p0 * p1 * p2);
        double ps = p0*p0 + p1*p1 + p2*p2;
        double ts = theta[0][i]*theta[0][i] + theta[1][i]*theta[1][i] + theta[2][i]*theta[2][i];
        double vs = phi_vel[0][i]*phi_vel[0][i] + phi_vel[1][i]*phi_vel[1][i] + phi_vel[2][i]*phi_vel[2][i];
        if (sqrt(ps) > phi_max) phi_max = sqrt(ps);
        P_int += P;
        E_kin += 0.5 * vs;
        E_mass += 0.5 * m2 * ps;
    }
    double dx3 = dx * dx * dx;
    E_total = (E_kin + E_mass) * dx3;  /* rough estimate, no gradient/potential terms */
    P_int *= dx3;

    printf("\nSeed diagnostics:\n");
    printf("  phi_max = %.4f\n", phi_max);
    printf("  P_int   = %.4e\n", P_int);
    printf("  E_kin   = %.4e (bg: %.4e)\n", E_kin * dx3, E_kin_bg * dx3);
    printf("  E_mass  = %.4e (bg: %.4e)\n", E_mass * dx3, E_pot_bg * dx3);
    printf("  E_total ~ %.4e (kin+mass, no gradient)\n", E_total);

    /* Write output seed */
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

    /* Cast to float and write */
    float *col = (float*)malloc(N3 * sizeof(float));
    void *cols[12];
    float *col_bufs[12];
    for (int c = 0; c < 12; c++) {
        col_bufs[c] = (float*)malloc(N3 * sizeof(float));
        cols[c] = col_bufs[c];
    }
    for (long i = 0; i < N3; i++) {
        col_bufs[0][i] = (float)phi[0][i];
        col_bufs[1][i] = (float)phi[1][i];
        col_bufs[2][i] = (float)phi[2][i];
        col_bufs[3][i] = (float)theta[0][i];
        col_bufs[4][i] = (float)theta[1][i];
        col_bufs[5][i] = (float)theta[2][i];
        col_bufs[6][i] = (float)phi_vel[0][i];
        col_bufs[7][i] = (float)phi_vel[1][i];
        col_bufs[8][i] = (float)phi_vel[2][i];
        col_bufs[9][i] = (float)theta_vel[0][i];
        col_bufs[10][i] = (float)theta_vel[1][i];
        col_bufs[11][i] = (float)theta_vel[2][i];
    }
    sfa_write_frame(out, 0.0, cols);
    sfa_close(out);

    printf("\nWrote seed: %s (%d^3, 12 cols f32)\n", outpath, N);

    /* Cleanup */
    for (int a = 0; a < 3; a++) {
        free(phi[a]); free(theta[a]); free(phi_vel[a]); free(theta_vel[a]);
    }
    for (int c = 0; c < 12; c++) free(col_bufs[c]);
    free(col); free(buf);
    sfa_close(in);
    return 0;
}
