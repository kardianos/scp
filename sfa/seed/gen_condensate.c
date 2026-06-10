/* gen_condensate.c — uniform charged condensate + long-wavelength noise seed
 * for the v66 complexified kernel (complex_phi=1, 12 fields + 12 velocities
 * = 24 columns). For Affleck-Dine fragmentation experiments (v67).
 *
 * Condensate (all voxels, all a = 0,1,2):
 *     u_a    = A (1 + xi(x))
 *     v_a    = 0
 *     udot_a = 0
 *     vdot_a = omega_c A (1 + xi(x))
 * where omega_c = sqrt(m2 + 2 Vt'(A^6) A^4) is the uniform-rotation frequency
 * of the homogeneous condensate, Vt(s) = (mu/2) s / (1 + kappa s),
 * Vt'(s) = (mu/2) / (1 + kappa s)^2, with standard params
 * m2 = 2.25, mu = -41.345, kappa = 50.
 *
 * xi(x) is a real random field of rms noise_amp built from ~200 random
 * plane waves with wavevectors snapped to the periodic-grid reciprocal
 * lattice k = (2 pi / (N dx)) n, n integer, restricted to
 * |k| <= 8 pi / (2 L) (long-wavelength seed noise). The SAME xi multiplies
 * u_a and vdot_a so the local phase (and local omega) stays coherent.
 * Theta sector (and its velocities) zero.
 *
 * Output: ONE SFA frame, f32, 24 columns named/typed exactly as the kernel
 * registers them (scp_sim.c sfa_add_column block / v66 SPEC §7.1).
 *
 * Build:
 *   gcc -O3 -march=native -fopenmp -o bin/gen_condensate \
 *       sfa/seed/gen_condensate.c -lzstd -lm
 *
 * Usage:
 *   gen_condensate N L A noise_amp rngseed out.sfa
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>

#define SFA_IMPLEMENTATION
#include "../format/sfa.h"

#define NCOLS   24
#define NWAVES  200

/* standard potential parameters (CLAUDE.md / v66 SPEC) */
#define M2     2.25
#define MU    (-41.345)
#define KAPPA  50.0

/* xorshift64* — deterministic, portable */
static uint64_t rng_state;
static uint64_t rng_next(void) {
    uint64_t x = rng_state;
    x ^= x >> 12; x ^= x << 25; x ^= x >> 27;
    rng_state = x;
    return x * 0x2545F4914F6CDD1DULL;
}
static double rng_uniform(void) {            /* [0,1) */
    return (double)(rng_next() >> 11) * (1.0 / 9007199254740992.0);
}

int main(int argc, char **argv) {
    if (argc != 7) {
        fprintf(stderr,
            "Usage: %s N L A noise_amp rngseed out.sfa\n"
            "  uniform condensate u_a=A, vdot_a=omega_c*A, multiplicative\n"
            "  long-wavelength noise of rms noise_amp on the modulus\n",
            argv[0]);
        return 1;
    }
    int      N         = atoi(argv[1]);
    double   L         = atof(argv[2]);
    double   A         = atof(argv[3]);
    double   noise_amp = atof(argv[4]);
    uint64_t rngseed   = (uint64_t)strtoull(argv[5], NULL, 10);
    const char *outpath = argv[6];

    if (N < 2 || L <= 0.0) {
        fprintf(stderr, "FATAL: bad N=%d or L=%g\n", N, L);
        return 1;
    }
    if (!(A > 0.0 && A < 0.7)) {
        fprintf(stderr, "FATAL: A=%g outside (0, 0.7)\n", A);
        return 1;
    }
    if (noise_amp < 0.0) {
        fprintf(stderr, "FATAL: noise_amp=%g < 0\n", noise_amp);
        return 1;
    }

    /* condensate rotation frequency: omega_c^2 = m2 + 2 Vt'(A^6) A^4 */
    double s0     = pow(A, 6.0);
    double Vtp    = (MU / 2.0) / ((1.0 + KAPPA * s0) * (1.0 + KAPPA * s0));
    double w2     = M2 + 2.0 * Vtp * s0 / (A * A);   /* 2 Vt'(A^6) A^4 */
    if (w2 <= 0.0) {
        fprintf(stderr, "FATAL: omega_c^2 = %g <= 0 (A=%g)\n", w2, A);
        return 1;
    }
    double omega_c = sqrt(w2);
    double rhoQ    = 3.0 * omega_c * A * A;          /* charge density */

    long   N3 = (long)N * N * N;
    double dx = 2.0 * L / (N - 1);   /* grid: x = -L + i*dx (scp_init.h) */
    double dV = dx * dx * dx;

    printf("gen_condensate: N=%d L=%.4f A=%.4f noise_amp=%.4g seed=%llu -> %s\n",
           N, L, A, noise_amp, (unsigned long long)rngseed, outpath);
    printf("  s0=A^6=%.6g  Vt'(s0)=%.6f  omega_c^2=%.6f  omega_c=%.6f\n",
           s0, Vtp, w2, omega_c);
    printf("  charge density 3*omega_c*A^2 = %.6f\n", rhoQ);

    /* ---- noise field xi: NWAVES box-commensurate plane waves ----------- */
    /* periodic grid period is N*dx (index wraps i=N -> i=0), so the
     * commensurate wavevectors are k = kfund * n, kfund = 2*pi/(N*dx).
     * Restrict |k| <= kmax = 8*pi/(2L) = 4*pi/L (long-wavelength seed). */
    rng_state = rngseed ? rngseed : 0x9E3779B97F4A7C15ULL;
    for (int w = 0; w < 8; w++) rng_next();          /* warm up */

    double kfund = 2.0 * M_PI / ((double)N * dx);
    double kmax  = 8.0 * M_PI / (2.0 * L);
    int    nmax  = (int)floor(kmax / kfund);
    double kx[NWAVES], ky[NWAVES], kz[NWAVES], amp[NWAVES], ph[NWAVES];
    int nw = 0;
    while (nw < NWAVES) {
        int nx = (int)floor(rng_uniform() * (2 * nmax + 1)) - nmax;
        int ny = (int)floor(rng_uniform() * (2 * nmax + 1)) - nmax;
        int nz = (int)floor(rng_uniform() * (2 * nmax + 1)) - nmax;
        double kk = kfund * sqrt((double)(nx*nx + ny*ny + nz*nz));
        if ((nx == 0 && ny == 0 && nz == 0) || kk > kmax) continue;
        kx[nw] = kfund * nx; ky[nw] = kfund * ny; kz[nw] = kfund * nz;
        /* Gaussian amplitude (Box-Muller) + uniform phase */
        double u1 = rng_uniform(), u2 = rng_uniform();
        if (u1 < 1e-300) u1 = 1e-300;
        amp[nw] = sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
        ph[nw]  = 2.0 * M_PI * rng_uniform();
        nw++;
    }
    printf("  noise: %d plane waves, kfund=%.6f, kmax=%.6f (|n|<=%d)\n",
           nw, kfund, kmax, nmax);

    double *xi = (double *)malloc((size_t)N3 * sizeof(double));
    if (!xi) { fprintf(stderr, "FATAL: alloc xi\n"); return 1; }

    double xsum = 0.0, x2sum = 0.0;
    #pragma omp parallel for reduction(+:xsum,x2sum) schedule(static)
    for (long idx = 0; idx < N3; idx++) {
        int i = (int)(idx / ((long)N * N));
        int j = (int)((idx / N) % N);
        int k = (int)(idx % N);
        double x = -L + i * dx, y = -L + j * dx, z = -L + k * dx;
        double v = 0.0;
        for (int w = 0; w < NWAVES; w++)
            v += amp[w] * cos(kx[w]*x + ky[w]*y + kz[w]*z + ph[w]);
        xi[idx] = v;
        xsum  += v;
        x2sum += v * v;
    }
    double mean = xsum / N3;
    double rms  = sqrt(x2sum / N3 - mean * mean);
    double scale = (rms > 1e-30 && noise_amp > 0.0) ? noise_amp / rms : 0.0;
    #pragma omp parallel for schedule(static)
    for (long idx = 0; idx < N3; idx++)
        xi[idx] = (xi[idx] - mean) * scale;
    printf("  noise: raw mean=%.4g rms=%.4g -> rescaled to rms=%.4g\n",
           mean, rms, noise_amp);

    /* ---- fill columns ---------------------------------------------------
     * column order MUST match kernel registration (scp_sim.c / SPEC §7.1):
     *  0-2  phi_x/y/z       (u_a)        12-14 phiim_x/y/z     (v_a)
     *  3-5  theta_x/y/z     (tu_a)       15-17 thetaim_x/y/z   (tv_a)
     *  6-8  phi_vx/vy/vz    (udot_a)     18-20 phiim_vx/vy/vz  (vdot_a)
     *  9-11 theta_vx/vy/vz  (tudot_a)    21-23 thetaim_vx/...  (tvdot_a)  */
    float *cols[NCOLS];
    for (int c = 0; c < NCOLS; c++) {
        cols[c] = (float *)calloc((size_t)N3, sizeof(float));
        if (!cols[c]) { fprintf(stderr, "FATAL: alloc column %d\n", c); return 1; }
    }

    double qsum = 0.0;   /* sum_a (u_a vdot_a - v_a udot_a) summed over grid */
    #pragma omp parallel for reduction(+:qsum) schedule(static)
    for (long idx = 0; idx < N3; idx++) {
        double f  = A * (1.0 + xi[idx]);
        double fd = omega_c * f;             /* vdot, same xi: phase coherent */
        for (int a = 0; a < 3; a++) {
            cols[0 + a][idx]  = (float)f;    /* u_a    */
            cols[18 + a][idx] = (float)fd;   /* vdot_a */
            /* v_a (12-14), udot_a (6-8), theta sector: stay zero */
        }
        qsum += 3.0 * f * fd;                /* u*vd - v*ud, v=ud=0 */
    }
    double Q_total = qsum * dV;
    printf("  Q_total = %.6f  (uniform estimate 3*omega_c*A^2*N^3*dV = %.6f)\n",
           Q_total, rhoQ * (double)N3 * dV);

    /* ---- write SFA ------------------------------------------------------ */
    double dt = 0.025 * dx;  /* standard dt factor (metadata only) */
    SFA *sfa = sfa_create(outpath, N, N, N, L, L, L, dt);
    if (!sfa) { fprintf(stderr, "FATAL: cannot create '%s'\n", outpath); return 1; }

    char bN[32], bL[32], bA[32], bn[32], bs[32], bw[32], bq[32];
    snprintf(bN, sizeof(bN), "%d", N);
    snprintf(bL, sizeof(bL), "%.6f", L);
    snprintf(bA, sizeof(bA), "%.6f", A);
    snprintf(bn, sizeof(bn), "%.6g", noise_amp);
    snprintf(bs, sizeof(bs), "%llu", (unsigned long long)rngseed);
    snprintf(bw, sizeof(bw), "%.6f", omega_c);
    snprintf(bq, sizeof(bq), "%.6f", Q_total);
    const char *keys[] = { "generator", "N", "L", "A", "noise_amp",
                           "rngseed", "omega_c", "Q_total" };
    const char *vals[] = { "gen_condensate", bN, bL, bA, bn, bs, bw, bq };
    sfa_add_kvmd(sfa, 0, 0xFFFFFFFF, 0xFFFFFFFF, keys, vals, 8);

    /* exact names/semantics/components as the kernel registers (24 cols) */
    sfa_add_column(sfa, "phi_x",      SFA_F32, SFA_POSITION, 0);
    sfa_add_column(sfa, "phi_y",      SFA_F32, SFA_POSITION, 1);
    sfa_add_column(sfa, "phi_z",      SFA_F32, SFA_POSITION, 2);
    sfa_add_column(sfa, "theta_x",    SFA_F32, SFA_ANGLE,    0);
    sfa_add_column(sfa, "theta_y",    SFA_F32, SFA_ANGLE,    1);
    sfa_add_column(sfa, "theta_z",    SFA_F32, SFA_ANGLE,    2);
    sfa_add_column(sfa, "phi_vx",     SFA_F32, SFA_VELOCITY, 0);
    sfa_add_column(sfa, "phi_vy",     SFA_F32, SFA_VELOCITY, 1);
    sfa_add_column(sfa, "phi_vz",     SFA_F32, SFA_VELOCITY, 2);
    sfa_add_column(sfa, "theta_vx",   SFA_F32, SFA_VELOCITY, 3);
    sfa_add_column(sfa, "theta_vy",   SFA_F32, SFA_VELOCITY, 4);
    sfa_add_column(sfa, "theta_vz",   SFA_F32, SFA_VELOCITY, 5);
    sfa_add_column(sfa, "phiim_x",    SFA_F32, SFA_POSITION, 3);
    sfa_add_column(sfa, "phiim_y",    SFA_F32, SFA_POSITION, 4);
    sfa_add_column(sfa, "phiim_z",    SFA_F32, SFA_POSITION, 5);
    sfa_add_column(sfa, "thetaim_x",  SFA_F32, SFA_ANGLE,    3);
    sfa_add_column(sfa, "thetaim_y",  SFA_F32, SFA_ANGLE,    4);
    sfa_add_column(sfa, "thetaim_z",  SFA_F32, SFA_ANGLE,    5);
    sfa_add_column(sfa, "phiim_vx",   SFA_F32, SFA_VELOCITY, 6);
    sfa_add_column(sfa, "phiim_vy",   SFA_F32, SFA_VELOCITY, 7);
    sfa_add_column(sfa, "phiim_vz",   SFA_F32, SFA_VELOCITY, 8);
    sfa_add_column(sfa, "thetaim_vx", SFA_F32, SFA_VELOCITY, 9);
    sfa_add_column(sfa, "thetaim_vy", SFA_F32, SFA_VELOCITY, 10);
    sfa_add_column(sfa, "thetaim_vz", SFA_F32, SFA_VELOCITY, 11);
    sfa_finalize_header(sfa);

    void *frame_cols[NCOLS];
    for (int c = 0; c < NCOLS; c++) frame_cols[c] = cols[c];
    sfa_write_frame(sfa, 0.0, frame_cols);
    sfa_close(sfa);

    printf("gen_condensate: wrote %s (%d^3, f32, 24 columns, 1 frame)\n",
           outpath, N);

    for (int c = 0; c < NCOLS; c++) free(cols[c]);
    free(xi);
    return 0;
}
