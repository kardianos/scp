/*  gen_spherical_seed.c — Generate proton/neutron seeds from spherical ansatz
 *
 *  The proton is modeled as a PHASE SOLITON: a localized region where the
 *  carrier phase shifts relative to the background. The ansatz is:
 *
 *    phi_a(r) = A(r) * cos(k*r + delta_a + chi_a * psi(r))
 *
 *  where:
 *    A(r) = A_bg + dA * exp(-r^2/(2*R_A^2))        -- amplitude envelope
 *    psi(r) = psi_0 * exp(-r^2/(2*R_psi^2))        -- phase shift
 *    chi_a = chirality per component (+1=Up, -1=Down)
 *    k = 1.5 (carrier wavenumber, from dispersion k^2 + m^2 = 4.5)
 *
 *  UUD (proton):  chi = (+1, +1, -1)
 *  UDD (neutron): chi = (+1, -1, -1)
 *
 *  Theta is initialized to algebraic equilibrium:
 *    theta_i = G * curl(phi)_i
 *    G = (alpha+eta) / (2*alpha + beta*|curl(phi)|^2)
 *
 *  Build: gcc -O3 -march=native -fopenmp -o gen_spherical_seed \
 *         gen_spherical_seed.c -lzstd -lm
 *
 *  Usage:
 *    ./gen_spherical_seed -N 128 -L 20 -o proton.sfa
 *    ./gen_spherical_seed -chirality UDD -o neutron.sfa
 *    ./gen_spherical_seed -dA 0.15 -R_A 3.5 -psi0 3.665 -R_psi 3.0 -o test.sfa
 */

#define SFA_IMPLEMENTATION
#include "../format/sfa.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>

#define PI 3.14159265358979323846
#define NCOLS 12   /* 3 phi + 3 theta + 3 phi_vel + 3 theta_vel */

/* ===== Default parameters ===== */
static int    N        = 128;
static double L        = 20.0;     /* half-domain */
static double A_bg     = 0.1;
static double dA       = 0.1;      /* amplitude boost at center */
static double R_A      = 3.5;      /* amplitude envelope radius */
static double psi0     = 3.665;    /* phase shift at center (210 deg, optimal UUD) */
static double R_psi    = 3.5;      /* phase envelope radius */
static double k_carrier= 1.5;     /* carrier wavenumber */
static double delta[3] = {0.0, 3.0005, 4.4325};  /* phase offsets */
static int    chi[3]   = {+1, +1, -1};            /* UUD chirality */
static double alpha_p  = 0.1;
static double eta_p    = 0.5;
static double beta_p   = 0.5;
static double cx = 0.0, cy = 0.0, cz = 0.0;  /* center of proton */
static char   outpath[512] = "spherical_seed.sfa";
static int    precision = 1;  /* 0=f16, 1=f32, 2=f64 */

static void usage(const char *prog) {
    fprintf(stderr,
        "Usage: %s [options] -o <output.sfa>\n"
        "\n"
        "Options:\n"
        "  -N <grid_size>    Grid size (default 128)\n"
        "  -L <half_domain>  Half-domain size (default 20)\n"
        "  -o <file.sfa>     Output file\n"
        "  -chirality UUD|UDD  Baryon type (default UUD=proton)\n"
        "  -dA <value>       Amplitude boost at center (default 0.1)\n"
        "  -R_A <value>      Amplitude envelope radius (default 3.5)\n"
        "  -psi0 <value>     Phase shift at center in radians (default 3.665)\n"
        "  -R_psi <value>    Phase envelope radius (default 3.5)\n"
        "  -k <value>        Carrier wavenumber (default 1.5)\n"
        "  -A_bg <value>     Background amplitude (default 0.1)\n"
        "  -center <x>,<y>,<z>  Proton center (default 0,0,0)\n"
        "  -precision f16|f32|f64  Output precision (default f32)\n"
        "  -delta <d0>,<d1>,<d2>  Phase offsets (default 0,3.0005,4.4325)\n",
        prog);
}

int main(int argc, char **argv) {
    /* Parse arguments */
    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "-N") && i+1 < argc) N = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-L") && i+1 < argc) L = atof(argv[++i]);
        else if (!strcmp(argv[i], "-o") && i+1 < argc) strncpy(outpath, argv[++i], 511);
        else if (!strcmp(argv[i], "-chirality") && i+1 < argc) {
            i++;
            if (!strcmp(argv[i], "UDD")) {
                chi[0] = +1; chi[1] = -1; chi[2] = -1;
                psi0 = 3.578;  /* optimal for UDD */
            } else if (!strcmp(argv[i], "UUD")) {
                chi[0] = +1; chi[1] = +1; chi[2] = -1;
                psi0 = 3.665;
            } else {
                fprintf(stderr, "Unknown chirality: %s (use UUD or UDD)\n", argv[i]);
                return 1;
            }
        }
        else if (!strcmp(argv[i], "-dA") && i+1 < argc) dA = atof(argv[++i]);
        else if (!strcmp(argv[i], "-R_A") && i+1 < argc) R_A = atof(argv[++i]);
        else if (!strcmp(argv[i], "-psi0") && i+1 < argc) psi0 = atof(argv[++i]);
        else if (!strcmp(argv[i], "-R_psi") && i+1 < argc) R_psi = atof(argv[++i]);
        else if (!strcmp(argv[i], "-k") && i+1 < argc) k_carrier = atof(argv[++i]);
        else if (!strcmp(argv[i], "-A_bg") && i+1 < argc) A_bg = atof(argv[++i]);
        else if (!strcmp(argv[i], "-center") && i+1 < argc) {
            sscanf(argv[++i], "%lf,%lf,%lf", &cx, &cy, &cz);
        }
        else if (!strcmp(argv[i], "-precision") && i+1 < argc) {
            i++;
            if (!strcmp(argv[i], "f16")) precision = 0;
            else if (!strcmp(argv[i], "f32")) precision = 1;
            else if (!strcmp(argv[i], "f64")) precision = 2;
        }
        else if (!strcmp(argv[i], "-delta") && i+1 < argc) {
            sscanf(argv[++i], "%lf,%lf,%lf", &delta[0], &delta[1], &delta[2]);
        }
        else if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "--help")) {
            usage(argv[0]); return 0;
        }
        else {
            fprintf(stderr, "Unknown option: %s\n", argv[i]);
            usage(argv[0]); return 1;
        }
    }

    double dx = 2.0 * L / N;
    double dt = 0.025 * dx;

    printf("gen_spherical_seed: N=%d L=%.1f dx=%.4f\n", N, L, dx);
    printf("  chirality: (%+d, %+d, %+d) = %s\n",
           chi[0], chi[1], chi[2],
           (chi[1] == 1) ? "UUD (proton)" : "UDD (neutron)");
    printf("  A_bg=%.3f dA=%.3f R_A=%.2f\n", A_bg, dA, R_A);
    printf("  psi0=%.4f (%.1f deg) R_psi=%.2f k=%.2f\n",
           psi0, psi0*180/PI, R_psi, k_carrier);
    printf("  delta=(%.4f, %.4f, %.4f)\n", delta[0], delta[1], delta[2]);
    printf("  center=(%.1f, %.1f, %.1f)\n", cx, cy, cz);
    printf("  output: %s (%s)\n", outpath,
           (const char*[]){"f16","f32","f64"}[precision]);

    /* Allocate per-column arrays (f64 internally) */
    long N3 = (long)N * N * N;
    double *fields[NCOLS];
    for (int c = 0; c < NCOLS; c++) {
        fields[c] = calloc(N3, sizeof(double));
        if (!fields[c]) { fprintf(stderr, "OOM\n"); return 1; }
    }

    /* Compute the Cosserat coupling prefactor for optimal psi */
    double CP_optimal = cos(delta[0]+chi[0]*psi0)
                      * cos(delta[1]+chi[1]*psi0)
                      * cos(delta[2]+chi[2]*psi0);
    printf("  C_P(psi0) = %.6f (triple product cosine factor)\n", CP_optimal);

    /* Track diagnostics */
    double P_max = 0, phi_sq_total = 0, theta_sq_total = 0;
    double curl_sq_max = 0;

    /* Fill grid */
    #pragma omp parallel for collapse(3) reduction(max:P_max,curl_sq_max) \
        reduction(+:phi_sq_total,theta_sq_total) schedule(static)
    for (int iz = 0; iz < N; iz++)
    for (int iy = 0; iy < N; iy++)
    for (int ix = 0; ix < N; ix++) {
        long idx = (long)iz * N * N + (long)iy * N + ix;

        /* World coordinates */
        double x = -L + (ix + 0.5) * dx - cx;
        double y = -L + (iy + 0.5) * dx - cy;
        double z = -L + (iz + 0.5) * dx - cz;
        double r2 = x*x + y*y + z*z;
        double r  = sqrt(r2 + 1e-30);  /* avoid division by zero */

        /* Radial profiles */
        double A_r   = A_bg + dA * exp(-r2 / (2.0 * R_A * R_A));
        double psi_r = psi0 * exp(-r2 / (2.0 * R_psi * R_psi));

        /* Derivatives */
        double dA_dr   = -dA * (r / (R_A * R_A)) * exp(-r2 / (2.0 * R_A * R_A));
        double dpsi_dr = -psi0 * (r / (R_psi * R_psi)) * exp(-r2 / (2.0 * R_psi * R_psi));

        /* phi_a values and radial derivatives */
        double phi[3], phi_prime[3];
        for (int a = 0; a < 3; a++) {
            double arg = k_carrier * r + delta[a] + chi[a] * psi_r;
            double ca = cos(arg);
            double sa = sin(arg);

            phi[a] = A_r * ca;
            /* d(phi_a)/dr = dA/dr * cos(arg) - A * (k + chi*dpsi/dr) * sin(arg) */
            phi_prime[a] = dA_dr * ca - A_r * (k_carrier + chi[a] * dpsi_dr) * sa;
        }

        /* Store phi (columns 0,1,2) */
        for (int a = 0; a < 3; a++)
            fields[a][idx] = phi[a];

        /* Triple product */
        double P = phi[0] * phi[1] * phi[2];
        double absP = fabs(P);
        if (absP > P_max) P_max = absP;

        /* Curl of phi: curl_i = sum_{j,k} eps_{ijk} (x_j/r) * phi_k'
         * d(phi_a)/dx_j = phi_a'(r) * x_j/r for spherically symmetric components
         */
        double nx = x/r, ny = y/r, nz = z/r;

        double curl_x = ny * phi_prime[2] - nz * phi_prime[1];
        double curl_y = nz * phi_prime[0] - nx * phi_prime[2];
        double curl_z = nx * phi_prime[1] - ny * phi_prime[0];

        double S = curl_x*curl_x + curl_y*curl_y + curl_z*curl_z;
        if (S > curl_sq_max) curl_sq_max = S;

        /* Theta: algebraic equilibrium
         * G = (alpha+eta) / (2*alpha + beta*S) */
        double G = (alpha_p + eta_p) / (2.0 * alpha_p + beta_p * S);

        double theta[3];
        theta[0] = G * curl_x;
        theta[1] = G * curl_y;
        theta[2] = G * curl_z;

        /* Store theta (columns 3,4,5) */
        for (int a = 0; a < 3; a++)
            fields[3 + a][idx] = theta[a];

        /* Velocities = 0 (columns 6-11 already zeroed by calloc) */

        /* Diagnostics */
        for (int a = 0; a < 3; a++) {
            phi_sq_total += phi[a] * phi[a];
            theta_sq_total += theta[a] * theta[a];
        }
    }

    printf("\nDiagnostics:\n");
    printf("  P_max = %.6e\n", P_max);
    printf("  |curl|^2_max = %.6f\n", curl_sq_max);
    printf("  phi_rms = %.6f\n", sqrt(phi_sq_total / (3.0 * N3)));
    printf("  theta_rms = %.6f\n", sqrt(theta_sq_total / (3.0 * N3)));

    /* Write SFA file */
    uint8_t sfa_dtype = (precision == 0) ? SFA_F16 : (precision == 1) ? SFA_F32 : SFA_F64;
    SFA *sfa = sfa_create(outpath, N, N, N, L, L, L, dt);

    /* Embed metadata */
    char buf[16][64];
    snprintf(buf[0], 64, "%d", N);
    snprintf(buf[1], 64, "%.6f", L);
    snprintf(buf[2], 64, "%.6f", A_bg);
    snprintf(buf[3], 64, "%.6f,%.6f,%.6f", delta[0], delta[1], delta[2]);
    snprintf(buf[4], 64, "%+d,%+d,%+d", chi[0], chi[1], chi[2]);
    snprintf(buf[5], 64, "%.4f", dA);
    snprintf(buf[6], 64, "%.4f", R_A);
    snprintf(buf[7], 64, "%.4f", psi0);
    snprintf(buf[8], 64, "%.4f", R_psi);
    snprintf(buf[9], 64, "%.4f", k_carrier);

    const char *keys[] = {"generator", "N", "L", "A_bg", "delta", "chirality",
                          "dA", "R_A", "psi0", "R_psi", "k"};
    const char *vals[] = {"gen_spherical_seed", buf[0], buf[1], buf[2], buf[3],
                          buf[4], buf[5], buf[6], buf[7], buf[8], buf[9]};
    sfa_add_kvmd(sfa, 0, 0xFFFFFFFF, 0xFFFFFFFF, keys, vals, 11);

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

    /* Convert to output precision and write */
    if (precision == 2) {
        void *cols[NCOLS];
        for (int c = 0; c < NCOLS; c++) cols[c] = fields[c];
        sfa_write_frame(sfa, 0.0, cols);
    } else {
        void *cols[NCOLS];
        int es = (precision == 0) ? 2 : 4;
        for (int c = 0; c < NCOLS; c++) {
            double *src = fields[c];
            cols[c] = malloc(N3 * es);
            if (precision == 1) {
                float *p = (float *)cols[c];
                for (long i = 0; i < N3; i++) p[i] = (float)src[i];
            } else {
                /* f16 */
                uint16_t *p = (uint16_t *)cols[c];
                for (long i = 0; i < N3; i++) {
                    float f = (float)src[i];
                    uint32_t x; memcpy(&x, &f, 4);
                    uint16_t sign = (x >> 16) & 0x8000;
                    int exp = ((x >> 23) & 0xFF) - 127 + 15;
                    uint16_t mant = (x >> 13) & 0x3FF;
                    p[i] = (exp <= 0) ? sign : (exp >= 31) ? (sign | 0x7C00) : (sign | (exp << 10) | mant);
                }
            }
        }
        sfa_write_frame(sfa, 0.0, cols);
        for (int c = 0; c < NCOLS; c++) free(cols[c]);
    }

    sfa_close(sfa);
    printf("\nWrote %s (%d^3 = %ld voxels, %d columns)\n",
           outpath, N, N3, NCOLS);

    for (int c = 0; c < NCOLS; c++) free(fields[c]);
    return 0;
}
