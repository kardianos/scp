/*  gen_base_blob.c — Generate a single z-oriented braid seed at any N/L
 *
 *  This creates the base seed that gen_lissajous_seed.c needs as input.
 *  Physics: same as init_braid in scp_init.h
 *
 *  Build: gcc -O3 -march=native -fopenmp -o gen_base_blob gen_base_blob.c -lzstd -lm
 *  Usage: gen_base_blob N L output.sfa [m] [A] [R] [A_bg] [ellip] [d0] [d1] [d2]
 */
#define SFA_IMPLEMENTATION
#include "../sfa/format/sfa.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define PI 3.14159265358979323846

int main(int argc, char **argv) {
    if (argc < 4) {
        fprintf(stderr, "Usage: %s N L output.sfa [m] [A] [R] [A_bg] [ellip] [d0] [d1] [d2]\n", argv[0]);
        return 1;
    }
    int N = atoi(argv[1]);
    double L = atof(argv[2]);
    const char *outpath = argv[3];
    double m = argc > 4 ? atof(argv[4]) : 1.5;
    double A = argc > 5 ? atof(argv[5]) : 0.4;
    double R = argc > 6 ? atof(argv[6]) : 3.0;
    double A_bg = argc > 7 ? atof(argv[7]) : 0.1;
    double ellip = argc > 8 ? atof(argv[8]) : 0.33;
    double delta[3] = {0, 3.0005, 4.4325};
    if (argc > 9)  delta[0] = atof(argv[9]);
    if (argc > 10) delta[1] = atof(argv[10]);
    if (argc > 11) delta[2] = atof(argv[11]);

    long N3 = (long)N*N*N;
    int NN = N*N;
    double dx = 2.0 * L / (N - 1);
    double m2 = m*m;
    double kw = PI/L, omega = sqrt(kw*kw + m2);
    double sx = 1+ellip, sy = 1-ellip;
    double inv2R2 = 1.0/(2*R*R);
    double k_bg = PI/L, omega_bg = sqrt(k_bg*k_bg + m2);

    printf("Base blob: N=%d L=%.2f dx=%.4f m=%.2f A=%.2f R=%.1f A_bg=%.2f ellip=%.2f\n",
           N, L, dx, m, A, R, A_bg, ellip);
    printf("  delta=(%.4f, %.4f, %.4f)\n", delta[0], delta[1], delta[2]);

    float *phi[3], *phi_vel[3];
    for (int a = 0; a < 3; a++) {
        phi[a] = (float*)calloc(N3, sizeof(float));
        phi_vel[a] = (float*)calloc(N3, sizeof(float));
    }

    #pragma omp parallel for schedule(static) collapse(2)
    for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
    for (int k = 0; k < N; k++) {
        double x = -L + i*dx;
        double y = -L + j*dx;
        double z = -L + k*dx;
        long idx = (long)i*NN + j*N + k;
        double r2e = x*x/(sx*sx) + y*y/(sy*sy);
        double env = exp(-r2e * inv2R2);
        for (int a = 0; a < 3; a++) {
            double ph = kw*z + delta[a];
            double ph_bg = k_bg*z + 2*PI*a/3.0;
            phi[a][idx] = (float)(A*env*cos(ph) + A_bg*cos(ph_bg));
            phi_vel[a][idx] = (float)(omega*A*env*sin(ph) + omega_bg*A_bg*sin(ph_bg));
        }
    }

    SFA *out = sfa_create(outpath, N, N, N, L, L, L, 1.0);
    out->flags = SFA_CODEC_BSS;
    sfa_add_column(out, "phi_x",  SFA_F32, SFA_POSITION, 0);
    sfa_add_column(out, "phi_y",  SFA_F32, SFA_POSITION, 1);
    sfa_add_column(out, "phi_z",  SFA_F32, SFA_POSITION, 2);
    sfa_add_column(out, "phi_vx", SFA_F32, SFA_VELOCITY, 0);
    sfa_add_column(out, "phi_vy", SFA_F32, SFA_VELOCITY, 1);
    sfa_add_column(out, "phi_vz", SFA_F32, SFA_VELOCITY, 2);
    sfa_finalize_header(out);

    void *cols[6] = {phi[0], phi[1], phi[2], phi_vel[0], phi_vel[1], phi_vel[2]};
    sfa_write_frame(out, 0.0, cols);
    sfa_close(out);

    printf("Wrote: %s\n", outpath);
    for (int a = 0; a < 3; a++) { free(phi[a]); free(phi_vel[a]); }
    return 0;
}
