/*  gen_em_wave.c — Generate a transverse EM wave packet seed
 *
 *  Creates a θ_x Gaussian wave packet propagating in +z direction
 *  on top of the standard φ background. This tests free-space photon
 *  propagation in the Cosserat theory.
 *
 *  The wave packet:
 *    θ_x(z) = A_theta × exp(-(z-z0)²/(2σ²)) × cos(k0 z)
 *    θ_x_vel(z) = A_theta × ω × exp(...) × sin(k0 z)   [rightward]
 *
 *  With ω = k0 (massless dispersion, c=1), this is a right-moving
 *  transverse wave polarized in x, propagating in z.
 *
 *  Build: gcc -O3 -o gen_em_wave gen_em_wave.c -lzstd -lm
 *  Usage: ./gen_em_wave -o wave_seed.sfa [-N 128] [-L 30] [-A_theta 0.05]
 *         [-sigma 4] [-wavelength 6] [-z0 -10]
 */

#define SFA_IMPLEMENTATION
#include "../format/sfa.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define NFIELDS 3
#define PI 3.14159265358979323846

int main(int argc, char **argv) {
    int N = 128;
    double L = 30.0;
    double A_bg = 0.1;
    double A_theta = 0.05;       /* wave packet amplitude */
    double sigma = 4.0;          /* Gaussian envelope width */
    double wavelength = 6.0;     /* carrier wavelength */
    double z0 = -10.0;           /* wave packet center */
    double m2 = 2.25;
    double dt_factor = 0.025;
    char outpath[512] = "em_wave_seed.sfa";

    for (int i = 1; i < argc - 1; i += 2) {
        const char *k = argv[i], *v = argv[i+1];
        if      (!strcmp(k,"-N"))          N = atoi(v);
        else if (!strcmp(k,"-L"))          L = atof(v);
        else if (!strcmp(k,"-A_bg"))       A_bg = atof(v);
        else if (!strcmp(k,"-A_theta"))    A_theta = atof(v);
        else if (!strcmp(k,"-sigma"))      sigma = atof(v);
        else if (!strcmp(k,"-wavelength")) wavelength = atof(v);
        else if (!strcmp(k,"-z0"))         z0 = atof(v);
        else if (!strcmp(k,"-o"))          strncpy(outpath, v, 511);
    }

    double k0 = 2.0 * PI / wavelength;
    double omega = k0;  /* massless: ω = k (c = 1) */

    long N3 = (long)N * N * N;
    double dx = 2.0 * L / (N - 1);
    double dt = dt_factor * dx;
    int NN = N * N;

    fprintf(stderr, "gen_em_wave: N=%d L=%.1f A_theta=%.3f sigma=%.1f wavelength=%.1f z0=%.1f\n",
            N, L, A_theta, sigma, wavelength, z0);
    fprintf(stderr, "  k0=%.4f omega=%.4f  (massless dispersion)\n", k0, omega);
    fprintf(stderr, "  Polarization: theta_x,  Propagation: +z\n");

    double *phi[3], *theta[3], *phi_vel[3], *theta_vel[3];
    for (int a = 0; a < 3; a++) {
        phi[a]       = (double*)calloc(N3, sizeof(double));
        phi_vel[a]   = (double*)calloc(N3, sizeof(double));
        theta[a]     = (double*)calloc(N3, sizeof(double));
        theta_vel[a] = (double*)calloc(N3, sizeof(double));
    }

    /* Standard φ background + θ wave packet */
    double k_bg = PI / L;
    double omega_bg = sqrt(k_bg * k_bg + m2);

    for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
    for (int kk = 0; kk < N; kk++) {
        long idx = (long)i*NN + j*N + kk;
        double z = -L + kk * dx;

        /* φ background */
        for (int a = 0; a < NFIELDS; a++) {
            double ph = k_bg * z + 2.0 * PI * a / 3.0;
            phi[a][idx] = A_bg * cos(ph);
            phi_vel[a][idx] = omega_bg * A_bg * sin(ph);
        }

        /* θ_x wave packet: Gaussian × carrier, propagating in +z */
        double dz = z - z0;
        double env = A_theta * exp(-dz * dz / (2.0 * sigma * sigma));
        theta[0][idx]     = env * cos(k0 * z);
        theta_vel[0][idx] = env * omega * sin(k0 * z);
        /* theta_y, theta_z = 0 (linearly polarized in x) */
    }}}

    /* Write SFA */
    uint8_t sfa_dtype = SFA_F64;
    SFA *sfa = sfa_create(outpath, N, N, N, L, L, L, dt);

    sfa_add_column(sfa, "phi_x",    sfa_dtype, SFA_POSITION, 0);
    sfa_add_column(sfa, "phi_y",    sfa_dtype, SFA_POSITION, 1);
    sfa_add_column(sfa, "phi_z",    sfa_dtype, SFA_POSITION, 2);
    sfa_add_column(sfa, "theta_x",  sfa_dtype, SFA_ANGLE,    0);
    sfa_add_column(sfa, "theta_y",  sfa_dtype, SFA_ANGLE,    1);
    sfa_add_column(sfa, "theta_z",  sfa_dtype, SFA_ANGLE,    2);
    sfa_add_column(sfa, "phi_vx",   sfa_dtype, SFA_VELOCITY, 0);
    sfa_add_column(sfa, "phi_vy",   sfa_dtype, SFA_VELOCITY, 1);
    sfa_add_column(sfa, "phi_vz",   sfa_dtype, SFA_VELOCITY, 2);
    sfa_add_column(sfa, "theta_vx", sfa_dtype, SFA_VELOCITY, 3);
    sfa_add_column(sfa, "theta_vy", sfa_dtype, SFA_VELOCITY, 4);
    sfa_add_column(sfa, "theta_vz", sfa_dtype, SFA_VELOCITY, 5);
    sfa_finalize_header(sfa);

    void *cols[12] = {
        phi[0], phi[1], phi[2],
        theta[0], theta[1], theta[2],
        phi_vel[0], phi_vel[1], phi_vel[2],
        theta_vel[0], theta_vel[1], theta_vel[2]
    };
    sfa_write_frame(sfa, 0.0, cols);
    sfa_close(sfa);

    fprintf(stderr, "  Wrote: %s\n", outpath);

    for (int a = 0; a < 3; a++) {
        free(phi[a]); free(phi_vel[a]);
        free(theta[a]); free(theta_vel[a]);
    }
    return 0;
}
