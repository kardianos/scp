/*  gen_diffuse.c — Smooth random chiral initial conditions for annealing
 *
 *  Generates a continuous 3D field with smooth θ domains of alternating
 *  chirality. No sharp edges — uses superposition of low-frequency
 *  sinusoidal modes to create a smooth random pattern.
 *
 *  φ: background carrier wave + smooth amplitude modulation
 *  θ: smooth random field with positive and negative regions (chirality domains)
 *
 *  The idea: start with a "warm" field where θ has structure at
 *  multiple scales, then let the equations sort it into particles
 *  via the V(P) binding and Cosserat dynamics.
 *
 *  Build: gcc -O3 -march=native -fopenmp -o gen_diffuse gen_diffuse.c -lzstd -lm
 */

#define SFA_IMPLEMENTATION
#include "../format/sfa.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define PI 3.14159265358979323846
#define NCOLS 12
#define MAX_MODES 64

typedef struct {
    double kx, ky, kz;  /* wavevector */
    double amp;          /* amplitude */
    double phase;        /* random phase */
} Mode;

/* Generate a smooth random 3D scalar field from superposition of sinusoidal modes.
 * Uses low-k modes (large wavelength) for smooth structure with no sharp edges. */
static void gen_random_modes(Mode *modes, int n_modes, double k_min, double k_max,
                              double amp, unsigned int *rng) {
    for (int i = 0; i < n_modes; i++) {
        /* Random wavevector with magnitude between k_min and k_max */
        double u1 = (double)(*rng = *rng * 1103515245 + 12345) / (double)0xFFFFFFFF;
        double u2 = (double)(*rng = *rng * 1103515245 + 12345) / (double)0xFFFFFFFF;
        double u3 = (double)(*rng = *rng * 1103515245 + 12345) / (double)0xFFFFFFFF;
        double u4 = (double)(*rng = *rng * 1103515245 + 12345) / (double)0xFFFFFFFF;

        /* Spherical random direction */
        double costh = 2*u1 - 1;
        double sinth = sqrt(1 - costh*costh);
        double phi = 2*PI*u2;
        double k_mag = k_min + (k_max - k_min) * u3;

        modes[i].kx = k_mag * sinth * cos(phi);
        modes[i].ky = k_mag * sinth * sin(phi);
        modes[i].kz = k_mag * costh;
        modes[i].amp = amp / sqrt(n_modes);  /* normalize by sqrt(N) for consistent RMS */
        modes[i].phase = 2*PI*u4;
    }
}

/* Evaluate smooth random field at point (x,y,z) */
static double eval_modes(const Mode *modes, int n_modes, double x, double y, double z) {
    double val = 0;
    for (int i = 0; i < n_modes; i++) {
        val += modes[i].amp * sin(modes[i].kx*x + modes[i].ky*y + modes[i].kz*z + modes[i].phase);
    }
    return val;
}

int main(int argc, char **argv) {
    int N = 384;
    double L = 50.0;
    double A_bg = 0.1;
    double theta_amp = 0.05;   /* RMS amplitude of θ perturbation */
    double phi_amp = 0.05;     /* RMS amplitude of φ modulation */
    double k_min = 0.3;        /* minimum mode wavenumber (large-scale structure) */
    double k_max = 1.5;        /* maximum mode wavenumber (carrier-scale) */
    int n_modes = 32;          /* number of random modes per field */
    double k_wave = 1.5;
    double delta[3] = {0.0, 3.0005, 4.4325};
    double Delta[3];
    unsigned int seed = 42;
    char outpath[512] = "diffuse_seed.sfa";
    int precision = 0; /* f16 default for GPU */

    Delta[0] = 0; Delta[1] = 2*PI/3; Delta[2] = 4*PI/3;

    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "-N") && i+1<argc) N = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-L") && i+1<argc) L = atof(argv[++i]);
        else if (!strcmp(argv[i], "-o") && i+1<argc) strncpy(outpath, argv[++i], 511);
        else if (!strcmp(argv[i], "-theta_amp") && i+1<argc) theta_amp = atof(argv[++i]);
        else if (!strcmp(argv[i], "-phi_amp") && i+1<argc) phi_amp = atof(argv[++i]);
        else if (!strcmp(argv[i], "-k_min") && i+1<argc) k_min = atof(argv[++i]);
        else if (!strcmp(argv[i], "-k_max") && i+1<argc) k_max = atof(argv[++i]);
        else if (!strcmp(argv[i], "-n_modes") && i+1<argc) n_modes = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-A_bg") && i+1<argc) A_bg = atof(argv[++i]);
        else if (!strcmp(argv[i], "-seed") && i+1<argc) seed = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-precision") && i+1<argc) {
            i++; precision = !strcmp(argv[i], "f16") ? 0 : 1;
        }
        else if (!strcmp(argv[i], "-h")) {
            fprintf(stderr, "Usage: %s [-N grid] [-L domain] [-theta_amp val] [-phi_amp val] "
                    "[-k_min val] [-k_max val] [-n_modes N] [-seed N] [-o output.sfa]\n", argv[0]);
            return 0;
        }
    }

    if (n_modes > MAX_MODES) n_modes = MAX_MODES;

    double dx = 2.0 * L / N;
    double dt = 0.025 * dx;
    long N3 = (long)N * N * N;
    int NN = N * N;

    printf("gen_diffuse: N=%d L=%.1f dx=%.4f\n", N, L, dx);
    printf("  theta_amp=%.3f phi_amp=%.3f k_range=[%.2f, %.2f] n_modes=%d\n",
           theta_amp, phi_amp, k_min, k_max, n_modes);
    printf("  A_bg=%.3f k_carrier=%.2f seed=%u\n", A_bg, k_wave, seed);

    /* Generate random mode sets:
     * - 3 sets for θ_x, θ_y, θ_z (independent random fields)
     * - 3 sets for φ amplitude modulation (per component)
     * - 3 sets for φ velocity modulation
     */
    unsigned int rng = seed;
    Mode theta_modes[3][MAX_MODES];
    Mode phi_mod_modes[3][MAX_MODES];
    Mode phi_vel_modes[3][MAX_MODES];

    for (int a = 0; a < 3; a++) {
        gen_random_modes(theta_modes[a], n_modes, k_min, k_max, theta_amp, &rng);
        gen_random_modes(phi_mod_modes[a], n_modes, k_min, k_max, phi_amp, &rng);
        gen_random_modes(phi_vel_modes[a], n_modes, k_min, k_max, phi_amp * k_wave, &rng);
    }

    /* Allocate fields */
    double *fields[NCOLS];
    for (int c = 0; c < NCOLS; c++) {
        fields[c] = calloc(N3, sizeof(double));
        if (!fields[c]) { fprintf(stderr, "OOM\n"); return 1; }
    }

    /* Fill grid */
    double P_max = 0, theta_max = 0, phi_max = 0;
    double omega = sqrt(k_wave * k_wave + 2.25); /* dispersion: ω² = k² + m² */

    #pragma omp parallel for collapse(3) reduction(max:P_max,theta_max,phi_max) schedule(static)
    for (int iz = 0; iz < N; iz++)
    for (int iy = 0; iy < N; iy++)
    for (int ix = 0; ix < N; ix++) {
        long idx = (long)iz * NN + (long)iy * N + ix;
        double x = -L + (ix + 0.5) * dx;
        double y = -L + (iy + 0.5) * dx;
        double z = -L + (iz + 0.5) * dx;

        /* φ_a = (A_bg + δA_a(x,y,z)) × cos(k×z + δ_a + Δ_a) */
        for (int a = 0; a < 3; a++) {
            double dA = eval_modes(phi_mod_modes[a], n_modes, x, y, z);
            double A_local = A_bg + dA;
            double phase = k_wave * z + delta[a] + Delta[a];
            fields[a][idx] = A_local * cos(phase);

            /* Velocity: ω × (A_bg + δA) × sin(phase) + perturbation */
            double dV = eval_modes(phi_vel_modes[a], n_modes, x, y, z);
            fields[6 + a][idx] = (A_bg * omega + dV) * sin(phase);

            double af = fabs(fields[a][idx]);
            if (af > phi_max) phi_max = af;
        }

        /* θ_a = smooth random field (independent per component) */
        for (int a = 0; a < 3; a++) {
            fields[3 + a][idx] = eval_modes(theta_modes[a], n_modes, x, y, z);
            /* θ velocity = 0 (let equations drive it) */

            double tf = fabs(fields[3 + a][idx]);
            if (tf > theta_max) theta_max = tf;
        }

        /* Track P_max */
        double P = fabs(fields[0][idx] * fields[1][idx] * fields[2][idx]);
        if (P > P_max) P_max = P;
    }

    printf("\nDiagnostics:\n");
    printf("  P_max = %.6e\n", P_max);
    printf("  phi_max = %.6e\n", phi_max);
    printf("  theta_max = %.6e\n", theta_max);

    /* Write SFA */
    uint8_t sfa_dtype = precision == 0 ? SFA_F16 : SFA_F32;
    SFA *sfa = sfa_create(outpath, N, N, N, L, L, L, dt);

    char buf[8][64];
    snprintf(buf[0], 64, "%.4f", theta_amp);
    snprintf(buf[1], 64, "%.4f", phi_amp);
    snprintf(buf[2], 64, "%.2f-%.2f", k_min, k_max);
    snprintf(buf[3], 64, "%d", n_modes);
    snprintf(buf[4], 64, "%u", seed);
    const char *keys[] = {"generator", "theta_amp", "phi_amp", "k_range", "n_modes", "random_seed"};
    const char *vals[] = {"gen_diffuse", buf[0], buf[1], buf[2], buf[3], buf[4]};
    sfa_add_kvmd(sfa, 0, 0xFFFFFFFF, 0xFFFFFFFF, keys, vals, 6);

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

    /* Convert and write */
    void *cols[NCOLS];
    int es = precision == 0 ? 2 : 4;
    for (int c = 0; c < NCOLS; c++) {
        cols[c] = malloc(N3 * es);
        if (precision == 1) {
            float *p = (float *)cols[c];
            for (long i = 0; i < N3; i++) p[i] = (float)fields[c][i];
        } else {
            uint16_t *p = (uint16_t *)cols[c];
            for (long i = 0; i < N3; i++) {
                float f = (float)fields[c][i];
                uint32_t x; memcpy(&x, &f, 4);
                uint16_t sign = (x >> 16) & 0x8000;
                int exp = ((x >> 23) & 0xFF) - 127 + 15;
                uint16_t mant = (x >> 13) & 0x3FF;
                p[i] = (exp <= 0) ? sign : (exp >= 31) ? (sign|0x7C00) : (sign|(exp<<10)|mant);
            }
        }
    }
    sfa_write_frame(sfa, 0.0, cols);
    for (int c = 0; c < NCOLS; c++) free(cols[c]);
    sfa_close(sfa);

    printf("\nWrote %s (%d^3, %d modes per field)\n", outpath, N, n_modes);
    for (int c = 0; c < NCOLS; c++) free(fields[c]);
    return 0;
}
