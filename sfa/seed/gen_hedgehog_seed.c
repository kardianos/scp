/*  gen_hedgehog_seed.c — Hedgehog-encoded proton/neutron seed
 *
 *  Each field component oscillates along its OWN axis:
 *    φ_a(x,y,z) = A(r) × cos(χ_a × k × coord_a + δ_a + Δ_a)
 *
 *  where coord_0 = x, coord_1 = y, coord_2 = z.
 *  This creates angular topology: P = φ₀φ₁φ₂ is confined because
 *  the three independent cosines average to zero over angles.
 *
 *  Unlike three separate tubes, this uses a SINGLE spherical envelope
 *  A(r) with self-consistent θ from curl(φ).
 *
 *  Build: gcc -O3 -march=native -fopenmp -o gen_hedgehog_seed \
 *         gen_hedgehog_seed.c -lzstd -lm
 */

#define SFA_IMPLEMENTATION
#include "../format/sfa.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>

#define PI 3.14159265358979323846
#define NCOLS 12

static int    N        = 128;
static double L        = 20.0;
static double A_bg     = 0.1;
static double dA       = 0.3;      /* amplitude boost — must be large for binding */
static double R_env    = 3.5;      /* envelope radius */
static double k_wave   = 1.5;
static double delta[3] = {0.0, 3.0005, 4.4325};
static double Delta[3];            /* carrier phases: 0, 2π/3, 4π/3 */
static int    chi[3]   = {+1, +1, -1};  /* UUD */
static double alpha_p  = 0.1;
static double eta_p    = 0.5;
static double beta_p   = 0.5;
static double cx = 0, cy = 0, cz = 0;
static char   outpath[512] = "hedgehog_seed.sfa";
static int    precision = 1;

static void usage(const char *prog) {
    fprintf(stderr,
        "Usage: %s [options] -o <output.sfa>\n"
        "  -N <grid>       Grid size (default 128)\n"
        "  -L <half_dom>   Half-domain (default 20)\n"
        "  -chirality UUD|UDD\n"
        "  -dA <value>     Amplitude boost (default 0.3)\n"
        "  -R <value>      Envelope radius (default 3.5)\n"
        "  -A_bg <value>   Background amplitude (default 0.1)\n"
        "  -k <value>      Carrier wavenumber (default 1.5)\n"
        "  -center <x>,<y>,<z>\n"
        "  -precision f16|f32\n"
        "  -delta <d0>,<d1>,<d2>\n",
        prog);
}

int main(int argc, char **argv) {
    Delta[0] = 0; Delta[1] = 2*PI/3; Delta[2] = 4*PI/3;

    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "-N") && i+1<argc) N = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-L") && i+1<argc) L = atof(argv[++i]);
        else if (!strcmp(argv[i], "-o") && i+1<argc) strncpy(outpath, argv[++i], 511);
        else if (!strcmp(argv[i], "-chirality") && i+1<argc) {
            i++;
            if (!strcmp(argv[i], "UDD")) { chi[0]=+1; chi[1]=-1; chi[2]=-1; }
            else { chi[0]=+1; chi[1]=+1; chi[2]=-1; } /* UUD default */
        }
        else if (!strcmp(argv[i], "-dA") && i+1<argc) dA = atof(argv[++i]);
        else if (!strcmp(argv[i], "-R") && i+1<argc) R_env = atof(argv[++i]);
        else if (!strcmp(argv[i], "-A_bg") && i+1<argc) A_bg = atof(argv[++i]);
        else if (!strcmp(argv[i], "-k") && i+1<argc) k_wave = atof(argv[++i]);
        else if (!strcmp(argv[i], "-center") && i+1<argc) {
            sscanf(argv[++i], "%lf,%lf,%lf", &cx, &cy, &cz);
        }
        else if (!strcmp(argv[i], "-precision") && i+1<argc) {
            i++; precision = !strcmp(argv[i], "f16") ? 0 : 1;
        }
        else if (!strcmp(argv[i], "-delta") && i+1<argc) {
            sscanf(argv[++i], "%lf,%lf,%lf", &delta[0], &delta[1], &delta[2]);
        }
        else if (!strcmp(argv[i], "-h")) { usage(argv[0]); return 0; }
    }

    double dx = 2.0 * L / N;
    double dt = 0.025 * dx;
    int NN = N * N;
    long N3 = (long)N * N * N;

    printf("gen_hedgehog_seed: N=%d L=%.1f dx=%.4f\n", N, L, dx);
    printf("  chirality: (%+d,%+d,%+d) = %s\n", chi[0], chi[1], chi[2],
           chi[1]==1 ? "UUD" : "UDD");
    printf("  A_bg=%.3f dA=%.3f R=%.2f k=%.2f\n", A_bg, dA, R_env, k_wave);

    double *fields[NCOLS];
    for (int c = 0; c < NCOLS; c++) {
        fields[c] = calloc(N3, sizeof(double));
        if (!fields[c]) { fprintf(stderr, "OOM\n"); return 1; }
    }

    double P_max = 0, Ep_total = 0, curl_sq_max = 0;
    double phi_sq = 0, theta_sq = 0;
    double mu = -41.345, kappa = 50.0;

    #pragma omp parallel for collapse(3) \
        reduction(max:P_max,curl_sq_max) reduction(+:phi_sq,theta_sq,Ep_total) \
        schedule(static)
    for (int iz = 0; iz < N; iz++)
    for (int iy = 0; iy < N; iy++)
    for (int ix = 0; ix < N; ix++) {
        long idx = (long)iz * NN + (long)iy * N + ix;

        double x = -L + (ix + 0.5) * dx - cx;
        double y = -L + (iy + 0.5) * dx - cy;
        double z = -L + (iz + 0.5) * dx - cz;
        double r2 = x*x + y*y + z*z;
        double r = sqrt(r2 + 1e-30);

        /* Spherical envelope */
        double env = exp(-r2 / (2.0 * R_env * R_env));
        double A_r = A_bg + dA * env;
        double dA_dr = -dA * (r / (R_env * R_env)) * env;

        /* Coordinates along each axis */
        double coord[3] = {x, y, z};

        /* φ_a = A(r) × cos(χ_a × k × coord_a + δ_a + Δ_a) */
        double phi[3];
        for (int a = 0; a < 3; a++) {
            double arg = chi[a] * k_wave * coord[a] + delta[a] + Delta[a];
            phi[a] = A_r * cos(arg);
        }

        /* Store phi */
        for (int a = 0; a < 3; a++)
            fields[a][idx] = phi[a];

        /* Triple product */
        double P = phi[0] * phi[1] * phi[2];
        double absP = fabs(P);
        if (absP > P_max) P_max = absP;

        /* V(P) */
        double P2 = P * P;
        double Vp = (mu / 2.0) * P2 / (1.0 + kappa * P2);
        Ep_total += Vp * dx * dx * dx;

        /* Curl of phi: need spatial derivatives
         * ∂φ_a/∂x_j = (∂A/∂x_j) cos(arg_a) - A × χ_a × k × δ_{aj} × sin(arg_a)
         * where ∂A/∂x_j = dA_dr × x_j/r
         */
        double dphi[3][3]; /* dphi[a][j] = ∂φ_a/∂x_j */
        double nx = x/r, ny = y/r, nz = z/r;
        double n_hat[3] = {nx, ny, nz};

        for (int a = 0; a < 3; a++) {
            double arg = chi[a] * k_wave * coord[a] + delta[a] + Delta[a];
            double ca = cos(arg);
            double sa = sin(arg);
            for (int j = 0; j < 3; j++) {
                /* Envelope gradient contribution */
                dphi[a][j] = dA_dr * n_hat[j] * ca;
                /* Carrier wave gradient: only when j == a */
                if (j == a) {
                    dphi[a][j] -= A_r * chi[a] * k_wave * sa;
                }
            }
        }

        /* curl(φ)_i = ε_{ijk} ∂φ_k/∂x_j */
        double curl[3];
        curl[0] = dphi[2][1] - dphi[1][2];  /* ∂φ_z/∂y - ∂φ_y/∂z */
        curl[1] = dphi[0][2] - dphi[2][0];  /* ∂φ_x/∂z - ∂φ_z/∂x */
        curl[2] = dphi[1][0] - dphi[0][1];  /* ∂φ_y/∂x - ∂φ_x/∂y */

        double S = curl[0]*curl[0] + curl[1]*curl[1] + curl[2]*curl[2];
        if (S > curl_sq_max) curl_sq_max = S;

        /* θ = G × curl(φ), G = (α+η)/(2α + β|curl|²) */
        double G = (alpha_p + eta_p) / (2.0 * alpha_p + beta_p * S);

        for (int a = 0; a < 3; a++)
            fields[3 + a][idx] = G * curl[a];

        /* Diagnostics */
        for (int a = 0; a < 3; a++) {
            phi_sq += phi[a] * phi[a];
            theta_sq += fields[3+a][idx] * fields[3+a][idx];
        }
    }

    printf("\nDiagnostics:\n");
    printf("  P_max = %.6e\n", P_max);
    printf("  E_pot = %.1f\n", Ep_total);
    printf("  |curl|^2_max = %.6f\n", curl_sq_max);
    printf("  phi_rms = %.6f\n", sqrt(phi_sq / (3.0 * N3)));
    printf("  theta_rms = %.6f\n", sqrt(theta_sq / (3.0 * N3)));

    /* Write SFA */
    uint8_t sfa_dtype = precision == 0 ? SFA_F16 : SFA_F32;
    SFA *sfa = sfa_create(outpath, N, N, N, L, L, L, dt);

    char buf[16][64];
    snprintf(buf[0], 64, "%+d,%+d,%+d", chi[0], chi[1], chi[2]);
    snprintf(buf[1], 64, "%.4f", dA);
    snprintf(buf[2], 64, "%.4f", R_env);
    const char *keys[] = {"generator", "chirality", "dA", "R_env"};
    const char *vals[] = {"gen_hedgehog_seed", buf[0], buf[1], buf[2]};
    sfa_add_kvmd(sfa, 0, 0xFFFFFFFF, 0xFFFFFFFF, keys, vals, 4);

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

    if (precision == 1) {
        void *cols[NCOLS];
        for (int c = 0; c < NCOLS; c++) {
            float *p = malloc(N3 * sizeof(float));
            for (long i = 0; i < N3; i++) p[i] = (float)fields[c][i];
            cols[c] = p;
        }
        sfa_write_frame(sfa, 0.0, cols);
        for (int c = 0; c < NCOLS; c++) free(cols[c]);
    } else {
        void *cols[NCOLS];
        for (int c = 0; c < NCOLS; c++) {
            uint16_t *p = malloc(N3 * 2);
            for (long i = 0; i < N3; i++) {
                float f = (float)fields[c][i];
                uint32_t x; memcpy(&x, &f, 4);
                uint16_t sign = (x >> 16) & 0x8000;
                int exp = ((x >> 23) & 0xFF) - 127 + 15;
                uint16_t mant = (x >> 13) & 0x3FF;
                p[i] = (exp <= 0) ? sign : (exp >= 31) ? (sign|0x7C00) : (sign|(exp<<10)|mant);
            }
            cols[c] = p;
        }
        sfa_write_frame(sfa, 0.0, cols);
        for (int c = 0; c < NCOLS; c++) free(cols[c]);
    }

    sfa_close(sfa);
    printf("\nWrote %s (%d^3, %d cols)\n", outpath, N, NCOLS);
    for (int c = 0; c < NCOLS; c++) free(fields[c]);
    return 0;
}
