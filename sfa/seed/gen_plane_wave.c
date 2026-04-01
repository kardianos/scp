/*  gen_plane_wave.c — Generate a polariton plane wave
 *
 *  Creates a true plane wave (uniform in x,y, propagating in z).
 *  This is what a photon looks like at the particle scale: flat sheets
 *  of oscillating θ_y + δφ_x propagating at v = sqrt(1 - η²/m²).
 *
 *  The wave fills the entire box — no Gaussian envelope.
 *  Use with periodic BC in z for steady-state, or absorbing for transient.
 *
 *  Build: gcc -O3 -o gen_plane_wave gen_plane_wave.c -lzstd -lm
 *  Usage: ./gen_plane_wave -o plane_wave.sfa [-N 128] [-L 30]
 *         [-A 0.005] [-wavelength 60] [-eta 0.5] [-m 1.5]
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
    double A_theta = 0.005;      /* small amplitude for linearity */
    double wavelength = 60.0;    /* full box = one half-wavelength */
    double m2 = 2.25;
    double eta = 0.5;
    double dt_factor = 0.025;
    char outpath[512] = "plane_wave.sfa";

    for (int i = 1; i < argc - 1; i += 2) {
        const char *k = argv[i], *v = argv[i+1];
        if      (!strcmp(k,"-N"))          N = atoi(v);
        else if (!strcmp(k,"-L"))          L = atof(v);
        else if (!strcmp(k,"-A"))          A_theta = atof(v);
        else if (!strcmp(k,"-wavelength")) wavelength = atof(v);
        else if (!strcmp(k,"-eta"))        eta = atof(v);
        else if (!strcmp(k,"-m"))        { double m = atof(v); m2 = m*m; }
        else if (!strcmp(k,"-o"))          strncpy(outpath, v, 511);
    }

    double k0 = 2.0 * PI / wavelength;
    double v_phase = sqrt(1.0 - eta*eta/m2);
    double omega = v_phase * k0;
    double mixing = -eta * k0 / m2;

    long N3 = (long)N * N * N;
    double dx = 2.0 * L / (N - 1);
    double dt = dt_factor * dx;
    int NN = N * N;

    fprintf(stderr, "gen_plane_wave: N=%d L=%.1f A=%.4f wavelength=%.1f\n",
            N, L, A_theta, wavelength);
    fprintf(stderr, "  k0=%.5f omega=%.5f v_phase=%.4f mixing=%.4f\n",
            k0, omega, v_phase, mixing);
    fprintf(stderr, "  This is a PLANE WAVE — uniform in x,y, propagating in +z\n");

    double *phi[3], *theta[3], *phi_vel[3], *theta_vel[3];
    for (int a = 0; a < 3; a++) {
        phi[a]       = (double*)calloc(N3, sizeof(double));
        phi_vel[a]   = (double*)calloc(N3, sizeof(double));
        theta[a]     = (double*)calloc(N3, sizeof(double));
        theta_vel[a] = (double*)calloc(N3, sizeof(double));
    }

    for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
    for (int kk = 0; kk < N; kk++) {
        long idx = (long)i*NN + j*N + kk;
        double z = -L + kk * dx;

        /* Plane wave polariton: uniform in x,y */
        theta[1][idx]     = A_theta * cos(k0 * z);
        theta_vel[1][idx] = A_theta * omega * sin(k0 * z);

        phi[0][idx]     = mixing * A_theta * cos(k0 * z);
        phi_vel[0][idx] = mixing * A_theta * omega * sin(k0 * z);
    }}}

    SFA *sfa = sfa_create(outpath, N, N, N, L, L, L, dt);
    uint8_t sfa_dtype = SFA_F64;
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
        phi[0], phi[1], phi[2], theta[0], theta[1], theta[2],
        phi_vel[0], phi_vel[1], phi_vel[2],
        theta_vel[0], theta_vel[1], theta_vel[2]
    };
    sfa_write_frame(sfa, 0.0, cols);
    sfa_close(sfa);

    fprintf(stderr, "  Wrote: %s\n", outpath);
    for (int a = 0; a < 3; a++) {
        free(phi[a]); free(phi_vel[a]); free(theta[a]); free(theta_vel[a]);
    }
    return 0;
}
