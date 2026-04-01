/*  gen_beam.c — Generate a directed polariton beam
 *
 *  A cylindrical beam: Gaussian in x,y (beam waist), extended plane-wave
 *  in z. This is what a laser looks like — directed, collimated,
 *  with visible wave fronts.
 *
 *  θ_y(x,y,z) = A × exp(-(x²+y²)/(2w²)) × cos(k₀z)
 *  δφ_x(x,y,z) = mixing × θ_y
 *
 *  The beam waist w controls how collimated it is.
 *  Large w → plane wave. Small w → diffracts quickly.
 *
 *  Build: gcc -O3 -o gen_beam gen_beam.c -lzstd -lm
 *  Usage: ./gen_beam -o beam.sfa [-N 128] [-L 30] [-A 0.01]
 *         [-waist 6] [-wavelength 20] [-eta 0.5] [-m 1.5]
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
    double A_theta = 0.01;
    double waist = 6.0;          /* beam waist (Gaussian width in x,y) */
    double wavelength = 20.0;    /* carrier wavelength along z */
    double m2 = 2.25;
    double eta = 0.5;
    double dt_factor = 0.025;
    char outpath[512] = "beam.sfa";

    for (int i = 1; i < argc - 1; i += 2) {
        const char *k = argv[i], *v = argv[i+1];
        if      (!strcmp(k,"-N"))          N = atoi(v);
        else if (!strcmp(k,"-L"))          L = atof(v);
        else if (!strcmp(k,"-A"))          A_theta = atof(v);
        else if (!strcmp(k,"-waist"))      waist = atof(v);
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

    /* Rayleigh range: z_R = π w² / λ — how far the beam stays collimated */
    double z_R = PI * waist * waist / wavelength;

    fprintf(stderr, "gen_beam: N=%d L=%.1f A=%.4f waist=%.1f wavelength=%.1f\n",
            N, L, A_theta, waist, wavelength);
    fprintf(stderr, "  k0=%.4f omega=%.4f v_phase=%.4f mixing=%.4f\n",
            k0, omega, v_phase, mixing);
    fprintf(stderr, "  Rayleigh range: z_R=%.1f (beam stays collimated for ±%.1f)\n",
            z_R, z_R);
    fprintf(stderr, "  Grid pts across beam: %.0f, grid pts per wavelength: %.0f\n",
            2*waist/dx, wavelength/dx);

    double *phi[3], *theta[3], *phi_vel[3], *theta_vel[3];
    for (int a = 0; a < 3; a++) {
        phi[a]       = (double*)calloc(N3, sizeof(double));
        phi_vel[a]   = (double*)calloc(N3, sizeof(double));
        theta[a]     = (double*)calloc(N3, sizeof(double));
        theta_vel[a] = (double*)calloc(N3, sizeof(double));
    }

    double inv2w2 = 1.0 / (2.0 * waist * waist);

    for (int i = 0; i < N; i++) { double x = -L + i * dx;
    for (int j = 0; j < N; j++) { double y = -L + j * dx;
    for (int kk = 0; kk < N; kk++) {
        long idx = (long)i*NN + j*N + kk;
        double z = -L + kk * dx;

        /* Slab beam: Gaussian in x only, uniform in y, plane wave in z.
         * This avoids cylindrical edge effects and diffraction rings.
         * The wavefronts are lines in y extending across the box. */
        double env = exp(-x*x * inv2w2);

        /* θ_y (main), δφ_x (coupled) */
        theta[1][idx]     = A_theta * env * cos(k0 * z);
        theta_vel[1][idx] = A_theta * env * omega * sin(k0 * z);

        phi[0][idx]     = mixing * A_theta * env * cos(k0 * z);
        phi_vel[0][idx] = mixing * A_theta * env * omega * sin(k0 * z);
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
