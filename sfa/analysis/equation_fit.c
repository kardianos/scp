/*  equation_fit.c — Fit equation-based model to proton template and measure error
 *
 *  Model: φ_a(x,y,z) = A(r) × cos(χ_a × k × z + δ_a + Δ_a + χ_a × ψ(r))
 *  where A(r) = A_bg + dA × exp(-r²/2R_A²)
 *        ψ(r) = ψ₀ × exp(-r²/2R_ψ²)
 *        r = distance from proton center
 *
 *  Fits the model parameters to minimize error against actual SFA data.
 *  Reports: parameter values, max error, RMS error, compression ratio.
 *
 *  Build: gcc -O3 -march=native -fopenmp -o equation_fit equation_fit.c -lzstd -lm
 */

#define SFA_IMPLEMENTATION
#include "../format/sfa.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define PI 3.14159265358979323846

static float f16f(uint16_t h) {
    uint16_t s = h & 0x8000; int e = (h >> 10) & 0x1F; uint16_t m = h & 0x3FF;
    if (e == 0) return 0; if (e == 31) return s ? -1e30f : 1e30f;
    float f; uint32_t x = ((uint32_t)s << 16) | ((uint32_t)(e-15+127) << 23) | ((uint32_t)m << 13);
    memcpy(&f, &x, 4); return f;
}

static float col_f(void *buf, SFA *sfa, int c, long i) {
    long N3 = (long)sfa->Nx * sfa->Ny * sfa->Nz;
    uint64_t off = 0;
    for (int cc = 0; cc < c; cc++)
        off += (uint64_t)N3 * sfa_dtype_size[sfa->columns[cc].dtype];
    int dt = sfa->columns[c].dtype;
    uint8_t *src = (uint8_t*)buf + off;
    if (dt == SFA_F16) return f16f(((uint16_t*)src)[i]);
    if (dt == SFA_F32) return ((float*)src)[i];
    return (float)((double*)src)[i];
}

/* Model parameters */
typedef struct {
    double cx, cy, cz;  /* proton center */
    double A_bg;         /* background amplitude */
    double k;            /* carrier wavenumber */
    double dA, R_A;      /* envelope: A_bg + dA*exp(-r²/2R²) */
    double psi0, R_psi;  /* phase shift */
    double delta[3];     /* phase offsets per component */
    double Delta[3];     /* carrier phases */
    int chi[3];          /* chirality */
} EqModel;

/* Evaluate model at a point */
static double model_eval(const EqModel *m, int a, double x, double y, double z) {
    double rx = x - m->cx, ry = y - m->cy, rz = z - m->cz;
    double r2 = rx*rx + ry*ry + rz*rz;
    double A_r = m->A_bg + m->dA * exp(-r2 / (2.0 * m->R_A * m->R_A));
    double psi_r = m->psi0 * exp(-r2 / (2.0 * m->R_psi * m->R_psi));
    double arg = m->chi[a] * m->k * z + m->delta[a] + m->Delta[a] + m->chi[a] * psi_r;
    return A_r * cos(arg);
}

/* Compute error of model against SFA data for one field component */
static void compute_error(const EqModel *m, int a, void *buf, SFA *sfa,
                          double *max_err, double *rms_err, double *max_val) {
    int N = sfa->Nx;
    long N3 = (long)N * N * N;
    double L = sfa->Lx;
    double dx = 2.0 * L / N;

    double me = 0, se = 0, mv = 0;

    #pragma omp parallel for collapse(3) reduction(max:me,mv) reduction(+:se)
    for (int iz = 0; iz < N; iz++)
    for (int iy = 0; iy < N; iy++)
    for (int ix = 0; ix < N; ix++) {
        long idx = (long)iz * N * N + (long)iy * N + ix;
        double x = -L + (ix + 0.5) * dx;
        double y = -L + (iy + 0.5) * dx;
        double z = -L + (iz + 0.5) * dx;

        double actual = col_f(buf, sfa, a, idx);
        double predicted = model_eval(m, a, x, y, z);
        double err = fabs(actual - predicted);
        double av = fabs(actual);

        if (err > me) me = err;
        if (av > mv) mv = av;
        se += err * err;
    }

    *max_err = me;
    *rms_err = sqrt(se / N3);
    *max_val = mv;
}

/* Simple grid search over key parameters to find best fit.
 * We know most parameters from theory — only dA, R_A, psi0, R_psi, and center need fitting. */
static void optimize_model(EqModel *m, void *buf, SFA *sfa) {
    int N = sfa->Nx;
    double L = sfa->Lx;
    double dx = 2.0 * L / N;

    /* Find proton center by maximum |P| */
    double best_P = 0;
    int bi = N/2, bj = N/2, bk = N/2;
    for (int iz = N/4; iz < 3*N/4; iz++)
    for (int iy = N/4; iy < 3*N/4; iy++)
    for (int ix = N/4; ix < 3*N/4; ix++) {
        long idx = (long)iz * N * N + (long)iy * N + ix;
        float p0 = col_f(buf, sfa, 0, idx);
        float p1 = col_f(buf, sfa, 1, idx);
        float p2 = col_f(buf, sfa, 2, idx);
        double P = fabs(p0 * p1 * p2);
        if (P > best_P) { best_P = P; bi = ix; bj = iy; bk = iz; }
    }
    m->cx = -L + (bi + 0.5) * dx;
    m->cy = -L + (bj + 0.5) * dx;
    m->cz = -L + (bk + 0.5) * dx;
    printf("  Center: (%.2f, %.2f, %.2f) [P_max=%.6e at voxel %d,%d,%d]\n",
           m->cx, m->cy, m->cz, best_P, bi, bj, bk);

    /* Grid search over dA, R_A, psi0, R_psi */
    double best_rms = 1e30;
    double best_dA = 0, best_RA = 3.5, best_psi = 0, best_Rpsi = 3.5;

    double dA_vals[] = {0.0, 0.02, 0.05, 0.1, 0.15, 0.2, 0.3, 0.5};
    double R_vals[] = {1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0};
    double psi_vals[] = {0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0};

    int n_dA = sizeof(dA_vals)/sizeof(double);
    int n_R = sizeof(R_vals)/sizeof(double);
    int n_psi = sizeof(psi_vals)/sizeof(double);

    printf("  Grid search: %d × %d × %d × %d = %d evaluations...\n",
           n_dA, n_R, n_psi, n_R, n_dA * n_R * n_psi * n_R);

    for (int idA = 0; idA < n_dA; idA++)
    for (int iRA = 0; iRA < n_R; iRA++)
    for (int ipsi = 0; ipsi < n_psi; ipsi++)
    for (int iRp = 0; iRp < n_R; iRp++) {
        m->dA = dA_vals[idA];
        m->R_A = R_vals[iRA];
        m->psi0 = psi_vals[ipsi];
        m->R_psi = R_vals[iRp];

        /* Quick error on phi_x only (component 0) */
        double me, re, mv;
        compute_error(m, 0, buf, sfa, &me, &re, &mv);

        if (re < best_rms) {
            best_rms = re;
            best_dA = m->dA;
            best_RA = m->R_A;
            best_psi = m->psi0;
            best_Rpsi = m->R_psi;
        }
    }

    m->dA = best_dA;
    m->R_A = best_RA;
    m->psi0 = best_psi;
    m->R_psi = best_Rpsi;

    printf("  Best fit: dA=%.3f R_A=%.2f psi0=%.3f R_psi=%.2f (rms=%.6e)\n",
           best_dA, best_RA, best_psi, best_Rpsi, best_rms);
}

int main(int argc, char **argv) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s input.sfa [-frame N]\n", argv[0]);
        return 1;
    }

    int frame = 0;
    for (int i = 2; i < argc; i++)
        if (!strcmp(argv[i], "-frame") && i+1 < argc) frame = atoi(argv[++i]);

    SFA *sfa = sfa_open(argv[1]);
    if (!sfa) { fprintf(stderr, "Cannot open %s\n", argv[1]); return 1; }

    int N = sfa->Nx;
    long N3 = (long)N * N * N;
    printf("equation_fit: %s (N=%d, frame=%d)\n", argv[1], N, frame);

    void *buf = malloc(sfa->frame_bytes);
    if (sfa_read_frame(sfa, frame, buf) != 0) {
        fprintf(stderr, "Cannot read frame %d\n", frame);
        return 1;
    }

    /* Initialize model with known parameters */
    EqModel m = {0};
    m.A_bg = 0.1;
    m.k = 1.5;
    m.delta[0] = 0; m.delta[1] = 3.0005; m.delta[2] = 4.4325;
    m.Delta[0] = 0; m.Delta[1] = 2*PI/3; m.Delta[2] = 4*PI/3;
    m.chi[0] = +1; m.chi[1] = +1; m.chi[2] = -1; /* UUD */

    printf("\nOptimizing model parameters...\n");
    optimize_model(&m, buf, sfa);

    printf("\n=== Equation model: 18 parameters ===\n");
    printf("  center = (%.2f, %.2f, %.2f)\n", m.cx, m.cy, m.cz);
    printf("  A_bg = %.3f, k = %.3f\n", m.A_bg, m.k);
    printf("  dA = %.3f, R_A = %.2f\n", m.dA, m.R_A);
    printf("  psi0 = %.3f (%.1f deg), R_psi = %.2f\n", m.psi0, m.psi0*180/PI, m.R_psi);
    printf("  delta = (%.4f, %.4f, %.4f)\n", m.delta[0], m.delta[1], m.delta[2]);
    printf("  chi = (%+d, %+d, %+d)\n", m.chi[0], m.chi[1], m.chi[2]);

    /* Evaluate error for all 3 phi components */
    printf("\n=== Error analysis (phi fields) ===\n");
    printf("%-8s  %12s  %12s  %12s  %10s\n", "Field", "MaxErr", "RMS", "MaxVal", "Rel%");

    double total_rms = 0;
    for (int a = 0; a < 3; a++) {
        double me, re, mv;
        compute_error(&m, a, buf, sfa, &me, &re, &mv);
        printf("phi_%c     %12.6e  %12.6e  %12.6e  %10.2f%%\n",
               'x'+a, me, re, mv, (mv > 0 ? me/mv*100 : 0));
        total_rms += re;
    }

    printf("\nMean RMS across phi fields: %.6e\n", total_rms / 3);
    printf("\nCompression: 18 params describe %ld values = %.0f× compression\n",
           N3 * 3, (double)(N3 * 3) / 18);

    /* Compare to generic tricubic */
    int n_blocks = (N/8) * (N/8) * (N/8);
    int tricubic_params = n_blocks * 64 * 3;
    printf("Generic tricubic: %d params (%.0f× compression)\n",
           tricubic_params, (double)(N3 * 3) / tricubic_params);
    printf("Equation advantage: %.0f× fewer params than tricubic\n",
           (double)tricubic_params / 18);

    free(buf);
    sfa_close(sfa);
    return 0;
}
