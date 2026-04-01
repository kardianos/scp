/*  gen_polariton.c — Generate a polariton eigenmode wave packet
 *
 *  Creates the CORRECT hybridized θ-φ mode, not a pure θ pulse.
 *  The photon-like eigenmode has:
 *    - θ_y component (main polarization)
 *    - δφ_x component with ratio δφ/δθ = -ηk/(m²) (from dispersion)
 *    - Phase velocity v = sqrt(1 - η²/m²)
 *
 *  For a z-propagating, x-polarized wave, the curl coupling pairs
 *  δφ_x with δθ_y (cross-polarization).
 *
 *  Build: gcc -O3 -o gen_polariton gen_polariton.c -lzstd -lm
 *  Usage: ./gen_polariton -o polariton.sfa [-N 128] [-L 30] [-A 0.02]
 *         [-sigma 5] [-wavelength 8] [-z0 -12] [-eta 0.5] [-m 1.5]
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
    double A_theta = 0.02;       /* wave packet amplitude (keep small for linearity) */
    double sigma = 5.0;          /* Gaussian envelope width */
    double wavelength = 8.0;     /* carrier wavelength */
    double z0 = -12.0;           /* wave packet center */
    double m2 = 2.25;            /* phi mass squared */
    double eta = 0.5;            /* curl coupling */
    double dt_factor = 0.025;
    char outpath[512] = "polariton_seed.sfa";

    for (int i = 1; i < argc - 1; i += 2) {
        const char *k = argv[i], *v = argv[i+1];
        if      (!strcmp(k,"-N"))          N = atoi(v);
        else if (!strcmp(k,"-L"))          L = atof(v);
        else if (!strcmp(k,"-A_bg"))       A_bg = atof(v);
        else if (!strcmp(k,"-A"))          A_theta = atof(v);
        else if (!strcmp(k,"-sigma"))      sigma = atof(v);
        else if (!strcmp(k,"-wavelength")) wavelength = atof(v);
        else if (!strcmp(k,"-z0"))         z0 = atof(v);
        else if (!strcmp(k,"-eta"))        eta = atof(v);
        else if (!strcmp(k,"-m"))        { double m = atof(v); m2 = m*m; }
        else if (!strcmp(k,"-o"))          strncpy(outpath, v, 511);
    }

    double k0 = 2.0 * PI / wavelength;
    double v_phase = sqrt(1.0 - eta*eta/m2);  /* polariton phase velocity */
    double omega = v_phase * k0;              /* ω = v·k (NOT k for polariton) */

    /* Mixing ratio: δφ_x / δθ_y = -η·k / m² at low k */
    double mixing = -eta * k0 / m2;

    long N3 = (long)N * N * N;
    double dx = 2.0 * L / (N - 1);
    double dt = dt_factor * dx;
    int NN = N * N;

    fprintf(stderr, "gen_polariton: N=%d L=%.1f A=%.4f sigma=%.1f wavelength=%.1f z0=%.1f\n",
            N, L, A_theta, sigma, wavelength, z0);
    fprintf(stderr, "  eta=%.3f m²=%.3f v_phase=%.4f (%.1f%% of c)\n", eta, m2, v_phase, v_phase*100);
    fprintf(stderr, "  k0=%.4f omega=%.4f mixing(δφ_x/δθ_y)=%.4f\n", k0, omega, mixing);
    fprintf(stderr, "  Polariton: θ_y + %.1f%% δφ_x, propagating in +z\n", fabs(mixing)*100);

    double *phi[3], *theta[3], *phi_vel[3], *theta_vel[3];
    for (int a = 0; a < 3; a++) {
        phi[a]       = (double*)calloc(N3, sizeof(double));
        phi_vel[a]   = (double*)calloc(N3, sizeof(double));
        theta[a]     = (double*)calloc(N3, sizeof(double));
        theta_vel[a] = (double*)calloc(N3, sizeof(double));
    }

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

        /* Polariton eigenmode wave packet:
         *   θ_y = A cos(k0 z) envelope     (main component)
         *   φ_x = mixing * A cos(k0 z) envelope  (coupled component)
         *
         * Velocities for +z propagation (ω = v_phase * k0):
         *   θ_y_vel = A ω sin(k0 z) envelope
         *   φ_x_vel = mixing * A ω sin(k0 z) envelope
         */
        double dz = z - z0;
        double env = A_theta * exp(-dz * dz / (2.0 * sigma * sigma));

        /* θ_y (main polarization — curl pairs x↔y for z-propagation) */
        theta[1][idx]     = env * cos(k0 * z);
        theta_vel[1][idx] = env * omega * sin(k0 * z);

        /* δφ_x (coupled component, with correct mixing ratio) */
        phi[0][idx]     += mixing * env * cos(k0 * z);
        phi_vel[0][idx] += mixing * env * omega * sin(k0 * z);
    }}}

    /* Write SFA */
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
