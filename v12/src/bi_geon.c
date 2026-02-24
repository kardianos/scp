/*
 * bi_geon.c — 3D Born-Infeld Hopfion in Yukawa Metric
 *
 * Combines V11 (3D BI evolution) with V9 (Yukawa gravity):
 *   - D-field BI formulation (analytical inversion, no Newton solver)
 *   - FFT Yukawa solver for 3D metric: (nabla^2 - mu^2) Phi = -kappa rho
 *   - ADM evolution: dB/dt = -curl(alpha E), dD/dt = curl(alpha H)
 *   - Hopf map initial data for B field
 *
 * PLAN.md Sections: 5 (Phase 3), 9 (Equations)
 *
 * Usage:
 *   bi_geon -evolve [-b 1.0] [-a 1.0] [-kappa 35] [-mu 6.47]
 *           [-N 128] [-L 6.0] [-T 10.0] [-metric_interval 10]
 *   bi_geon -flat    (no gravity, V11 comparison)
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include <fftw3.h>

/* ---- Parameters ---- */
static double b_param     = 1.0;    /* BI field strength */
static double a_param     = 1.0;    /* hopfion size */
static double kappa_param = 35.0;   /* Yukawa coupling */
static double mu_param    = 6.47;   /* Yukawa mass (1/range) */
static int    N           = 128;    /* grid points per dim */
static double L           = 6.0;    /* box half-width */
static double T_final     = 10.0;   /* total time */
static int    metric_interval = 10; /* metric update frequency (timesteps) */
static int    gravity_on  = 1;      /* 0 = flat space (V11 comparison) */

/* ---- Derived ---- */
static double h, dt;
static int    N3;

/* ---- State: D and B fields ---- */
static double *Dx, *Dy, *Dz;
static double *Bx, *By, *Bz;

/* ---- Working arrays ---- */
static double *Ex, *Ey, *Ez;     /* E = dH/dD (analytical) */
static double *Hfx, *Hfy, *Hfz; /* H = dH/dB (analytical) */
static double *aEx, *aEy, *aEz; /* alpha*E product (for curl) */
static double *aHx, *aHy, *aHz; /* alpha*H product (for curl) */

/* ---- Metric ---- */
static double *Phi;              /* gravitational potential */
static double *alpha_lapse;      /* lapse function = 1/sqrt(1 - 2*Phi) */
static double *rho_BI;           /* BI energy density */

/* ---- FFT arrays ---- */
static fftw_plan fft_forward, fft_backward;
static double *fft_in;
static fftw_complex *fft_out;

/* ---- RK4 scratch ---- */
static double *Dx0, *Dy0, *Dz0, *Bx0, *By0, *Bz0;  /* state backup */
static double *kDx, *kDy, *kDz, *kBx, *kBy, *kBz;  /* RHS accumulator */
static double *rdx, *rdy, *rdz, *rbx, *rby, *rbz;   /* temp RHS */

/* ---- Grid ---- */
static inline int idx(int i, int j, int k) {
    i = ((i % N) + N) % N;
    j = ((j % N) + N) % N;
    k = ((k % N) + N) % N;
    return (i * N + j) * N + k;
}

static inline double grid_x(int i) { return -L + (i + 0.5) * h; }
static inline double grid_y(int j) { return -L + (j + 0.5) * h; }
static inline double grid_z(int k) { return -L + (k + 0.5) * h; }

/* ---- Setup ---- */
static void setup_grid(void) {
    h  = 2.0 * L / N;
    N3 = N * N * N;
    dt = 0.4 * h / sqrt(3.0);
    printf("# Grid: N=%d, L=%.1f, h=%.4f, dt=%.6f, N3=%d\n", N, L, h, dt, N3);
    printf("# Memory: %.1f MB (fields) + %.1f MB (FFT)\n",
           42.0 * N3 * 8.0 / 1e6,
           (N3 + (N*N*(N/2+1))) * 8.0 * 2.0 / 1e6);
    printf("# Gravity: %s (kappa=%.2f, mu=%.4f)\n",
           gravity_on ? "ON" : "OFF", kappa_param, mu_param);
}

static void alloc_fields(void) {
    /* State */
    Dx = calloc(N3, sizeof(double));  Dy = calloc(N3, sizeof(double));
    Dz = calloc(N3, sizeof(double));
    Bx = calloc(N3, sizeof(double));  By = calloc(N3, sizeof(double));
    Bz = calloc(N3, sizeof(double));

    /* Working */
    Ex = calloc(N3, sizeof(double));  Ey = calloc(N3, sizeof(double));
    Ez = calloc(N3, sizeof(double));
    Hfx = calloc(N3, sizeof(double)); Hfy = calloc(N3, sizeof(double));
    Hfz = calloc(N3, sizeof(double));
    aEx = calloc(N3, sizeof(double)); aEy = calloc(N3, sizeof(double));
    aEz = calloc(N3, sizeof(double));
    aHx = calloc(N3, sizeof(double)); aHy = calloc(N3, sizeof(double));
    aHz = calloc(N3, sizeof(double));

    /* Metric */
    Phi = calloc(N3, sizeof(double));
    alpha_lapse = malloc(N3 * sizeof(double));
    rho_BI = calloc(N3, sizeof(double));
    for (int p = 0; p < N3; p++) alpha_lapse[p] = 1.0;

    /* FFT */
    int N_fft = N * N * (N / 2 + 1);
    fft_in  = fftw_malloc(N3 * sizeof(double));
    fft_out = fftw_malloc(N_fft * sizeof(fftw_complex));
    fft_forward  = fftw_plan_dft_r2c_3d(N, N, N, fft_in, fft_out, FFTW_MEASURE);
    fft_backward = fftw_plan_dft_c2r_3d(N, N, N, fft_out, fft_in, FFTW_MEASURE);

    /* RK4 scratch */
    Dx0 = calloc(N3, sizeof(double)); Dy0 = calloc(N3, sizeof(double));
    Dz0 = calloc(N3, sizeof(double));
    Bx0 = calloc(N3, sizeof(double)); By0 = calloc(N3, sizeof(double));
    Bz0 = calloc(N3, sizeof(double));
    kDx = calloc(N3, sizeof(double)); kDy = calloc(N3, sizeof(double));
    kDz = calloc(N3, sizeof(double));
    kBx = calloc(N3, sizeof(double)); kBy = calloc(N3, sizeof(double));
    kBz = calloc(N3, sizeof(double));
    rdx = calloc(N3, sizeof(double)); rdy = calloc(N3, sizeof(double));
    rdz = calloc(N3, sizeof(double));
    rbx = calloc(N3, sizeof(double)); rby = calloc(N3, sizeof(double));
    rbz = calloc(N3, sizeof(double));

    if (!Dx || !Dy || !Dz || !Bx || !By || !Bz ||
        !Ex || !Ey || !Ez || !Hfx || !Hfy || !Hfz ||
        !aEx || !aEy || !aEz || !aHx || !aHy || !aHz ||
        !Phi || !alpha_lapse || !rho_BI ||
        !fft_in || !fft_out ||
        !Dx0 || !Dy0 || !Dz0 || !Bx0 || !By0 || !Bz0 ||
        !kDx || !kDy || !kDz || !kBx || !kBy || !kBz ||
        !rdx || !rdy || !rdz || !rbx || !rby || !rbz) {
        fprintf(stderr, "Allocation failed\n");
        exit(1);
    }
}

static void free_fields(void) {
    free(Dx); free(Dy); free(Dz);
    free(Bx); free(By); free(Bz);
    free(Ex); free(Ey); free(Ez);
    free(Hfx); free(Hfy); free(Hfz);
    free(aEx); free(aEy); free(aEz);
    free(aHx); free(aHy); free(aHz);
    free(Phi); free(alpha_lapse); free(rho_BI);
    fftw_free(fft_in); fftw_free(fft_out);
    fftw_destroy_plan(fft_forward);
    fftw_destroy_plan(fft_backward);
    free(Dx0); free(Dy0); free(Dz0);
    free(Bx0); free(By0); free(Bz0);
    free(kDx); free(kDy); free(kDz);
    free(kBx); free(kBy); free(kBz);
    free(rdx); free(rdy); free(rdz);
    free(rbx); free(rby); free(rbz);
}

/* ==================================================================
 * Hopf map initial data (pure magnetic, D=0)
 * PLAN.md Section 5: Hopfion B field from V11
 * ================================================================== */
static void init_hopfion_B(void) {
    double a = a_param;
    double a2 = a * a;
    /* C0 chosen so peak |B| ~ b_param (from V11) */
    double C0 = b_param * a2 / 16.0;

    #pragma omp parallel for schedule(static) collapse(3)
    for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
    for (int k = 0; k < N; k++) {
        double x = grid_x(i);
        double y = grid_y(j);
        double z = grid_z(k);
        double r2 = x*x + y*y + z*z;
        double s = a2 + r2;
        double s2 = s * s;

        /* Stereographic embedding S3 -> R3 */
        double u1 = (a2 - r2) / s;
        double v1 = 2.0 * a * z / s;
        double u2 = 2.0 * a * x / s;
        double v2 = 2.0 * a * y / s;

        /* Partial derivatives */
        double du1_dx = -4.0 * a2 * x / s2;
        double du1_dy = -4.0 * a2 * y / s2;
        double du1_dz = -4.0 * a2 * z / s2;

        double dv1_dx = -4.0 * a * x * z / s2;
        double dv1_dy = -4.0 * a * y * z / s2;
        double dv1_dz =  2.0 * a * (s - 2.0 * z * z) / s2;

        double du2_dx =  2.0 * a * (s - 2.0 * x * x) / s2;
        double du2_dy = -4.0 * a * x * y / s2;
        double du2_dz = -4.0 * a * x * z / s2;

        double dv2_dx = -4.0 * a * x * y / s2;
        double dv2_dy =  2.0 * a * (s - 2.0 * y * y) / s2;
        double dv2_dz = -4.0 * a * y * z / s2;

        /* CP1 Hopf map: n = (n1, n2, n3) on S2 */
        double n1 = 2.0 * (u1*u2 + v1*v2);
        double n2 = 2.0 * (v1*u2 - u1*v2);
        double n3 = u1*u1 + v1*v1 - u2*u2 - v2*v2;

        /* dn/dx_j by chain rule */
        double dn1_dx = 2.0*(du1_dx*u2 + u1*du2_dx + dv1_dx*v2 + v1*dv2_dx);
        double dn1_dy = 2.0*(du1_dy*u2 + u1*du2_dy + dv1_dy*v2 + v1*dv2_dy);
        double dn1_dz = 2.0*(du1_dz*u2 + u1*du2_dz + dv1_dz*v2 + v1*dv2_dz);

        double dn2_dx = 2.0*(dv1_dx*u2 + v1*du2_dx - du1_dx*v2 - u1*dv2_dx);
        double dn2_dy = 2.0*(dv1_dy*u2 + v1*du2_dy - du1_dy*v2 - u1*dv2_dy);
        double dn2_dz = 2.0*(dv1_dz*u2 + v1*du2_dz - du1_dz*v2 - u1*dv2_dz);

        double dn3_dx = 2.0*(u1*du1_dx + v1*dv1_dx - u2*du2_dx - v2*dv2_dx);
        double dn3_dy = 2.0*(u1*du1_dy + v1*dv1_dy - u2*du2_dy - v2*dv2_dy);
        double dn3_dz = 2.0*(u1*du1_dz + v1*dv1_dz - u2*du2_dz - v2*dv2_dz);

        /* B = C0 * n . (dn/dy x dn/dz, dn/dz x dn/dx, dn/dx x dn/dy) */
        double cy_x = dn2_dy * dn3_dz - dn3_dy * dn2_dz;
        double cy_y = dn3_dy * dn1_dz - dn1_dy * dn3_dz;
        double cy_z = dn1_dy * dn2_dz - dn2_dy * dn1_dz;

        double cz_x = dn2_dz * dn3_dx - dn3_dz * dn2_dx;
        double cz_y = dn3_dz * dn1_dx - dn1_dz * dn3_dx;
        double cz_z = dn1_dz * dn2_dx - dn2_dz * dn1_dx;

        double cx_x = dn2_dx * dn3_dy - dn3_dx * dn2_dy;
        double cx_y = dn3_dx * dn1_dy - dn1_dx * dn3_dy;
        double cx_z = dn1_dx * dn2_dy - dn2_dx * dn1_dy;

        int p = idx(i, j, k);
        Bx[p] = C0 * (n1 * cy_x + n2 * cy_y + n3 * cy_z);
        By[p] = C0 * (n1 * cz_x + n2 * cz_y + n3 * cz_z);
        Bz[p] = C0 * (n1 * cx_x + n2 * cx_y + n3 * cx_z);

        /* D = 0 initially (pure magnetic) */
        Dx[p] = Dy[p] = Dz[p] = 0.0;
    }
}

/* ==================================================================
 * D-field BI constitutive relations (ANALYTICAL — no Newton solver)
 * PLAN.md Section 9: Analytical BI Constitutive Relations
 *
 * R = 1 + (D^2 + B^2)/b^2 + |D x B|^2/b^4
 * E = dH/dD = (D(1 + B^2/b^2) - B(B.D)/b^2) / sqrt(R)
 * H = dH/dB = (B(1 + D^2/b^2) - D(D.B)/b^2) / sqrt(R)
 * ================================================================== */
static inline void bi_constitutive_Dfield(
    double dx_f, double dy_f, double dz_f,
    double bx_f, double by_f, double bz_f,
    double b,
    double *ex_out, double *ey_out, double *ez_out,
    double *hx_out, double *hy_out, double *hz_out,
    double *rho_out)
{
    double b2 = b * b;
    double b4 = b2 * b2;

    double D2 = dx_f*dx_f + dy_f*dy_f + dz_f*dz_f;
    double B2 = bx_f*bx_f + by_f*by_f + bz_f*bz_f;
    double DdB = dx_f*bx_f + dy_f*by_f + dz_f*bz_f;

    /* Cross product D x B */
    double DxBx = dy_f*bz_f - dz_f*by_f;
    double DxBy = dz_f*bx_f - dx_f*bz_f;
    double DxBz = dx_f*by_f - dy_f*bx_f;
    double P2 = DxBx*DxBx + DxBy*DxBy + DxBz*DxBz;

    /* R = 1 + (D^2 + B^2)/b^2 + |D x B|^2/b^4 */
    double R = 1.0 + (D2 + B2) / b2 + P2 / b4;
    /* R >= 1 always (sum of 1 + non-negative terms) */
    double sqrtR = sqrt(R);
    double inv_sqrtR = 1.0 / sqrtR;

    /* E = (D(1 + B^2/b^2) - B(B.D)/b^2) / sqrt(R) */
    double fac_D = (1.0 + B2 / b2) * inv_sqrtR;
    double fac_B = (DdB / b2) * inv_sqrtR;
    *ex_out = dx_f * fac_D - bx_f * fac_B;
    *ey_out = dy_f * fac_D - by_f * fac_B;
    *ez_out = dz_f * fac_D - bz_f * fac_B;

    /* H = (B(1 + D^2/b^2) - D(D.B)/b^2) / sqrt(R) */
    double fac_B2 = (1.0 + D2 / b2) * inv_sqrtR;
    double fac_D2 = (DdB / b2) * inv_sqrtR;
    *hx_out = bx_f * fac_B2 - dx_f * fac_D2;
    *hy_out = by_f * fac_B2 - dy_f * fac_D2;
    *hz_out = bz_f * fac_B2 - dz_f * fac_D2;

    /* Energy density: rho_BI = b^2 (sqrt(R) - 1) */
    *rho_out = b2 * (sqrtR - 1.0);
}

/* ==================================================================
 * 4th-order spatial derivatives (same as V11)
 * ================================================================== */
static inline double deriv4(const double *f, int i, int j, int k, int dim) {
    int di = (dim == 0), dj = (dim == 1), dk = (dim == 2);
    double fm2 = f[idx(i - 2*di, j - 2*dj, k - 2*dk)];
    double fm1 = f[idx(i - di,   j - dj,   k - dk)];
    double fp1 = f[idx(i + di,   j + dj,   k + dk)];
    double fp2 = f[idx(i + 2*di, j + 2*dj, k + 2*dk)];
    return (-fp2 + 8.0*fp1 - 8.0*fm1 + fm2) / (12.0 * h);
}

/* Curl at a point: curl(F) */
static inline void curl_at(
    const double *fx, const double *fy, const double *fz,
    int i, int j, int k,
    double *cx, double *cy, double *cz)
{
    *cx = deriv4(fz, i, j, k, 1) - deriv4(fy, i, j, k, 2);
    *cy = deriv4(fx, i, j, k, 2) - deriv4(fz, i, j, k, 0);
    *cz = deriv4(fy, i, j, k, 0) - deriv4(fx, i, j, k, 1);
}

/* ==================================================================
 * FFT Yukawa solver
 * PLAN.md Section 9: 3D Yukawa via FFT
 *
 * (nabla^2 - mu^2) Phi = -kappa * rho_BI
 * => Phi(k) = -kappa * rho(k) / (k^2 + mu^2)
 * ================================================================== */
static void solve_yukawa_fft(void) {
    double dk = 2.0 * M_PI / (N * h);  /* wavenumber spacing */
    double mu2 = mu_param * mu_param;
    double norm = 1.0 / (double)N3;    /* FFTW normalization */
    int Nk = N / 2 + 1;                /* complex array size in last dim */

    /* Copy rho_BI into FFT input */
    memcpy(fft_in, rho_BI, N3 * sizeof(double));

    /* Forward FFT: rho(x) -> rho(k) */
    fftw_execute(fft_forward);

    /* Multiply by Green's function in k-space */
    #pragma omp parallel for schedule(static) collapse(3)
    for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
    for (int kk = 0; kk < Nk; kk++) {
        int ii = (i <= N/2) ? i : i - N;
        int jj = (j <= N/2) ? j : j - N;

        double kx = ii * dk;
        double ky = jj * dk;
        double kz = kk * dk;
        double k2 = kx*kx + ky*ky + kz*kz;

        int p = (i * N + j) * Nk + kk;

        if (k2 + mu2 < 1e-30) {
            /* Zero mode: set to zero (no constant offset) */
            fft_out[p][0] = 0.0;
            fft_out[p][1] = 0.0;
        } else {
            /* Phi(k) = kappa * rho(k) / (k^2 + mu^2) */
            /* Note: sign convention. Source is -kappa*rho on RHS, */
            /* so Phi(k) = kappa * rho(k) / (k^2 + mu^2) for attractive well */
            double G = kappa_param / (k2 + mu2) * norm;
            fft_out[p][0] *= -G;
            fft_out[p][1] *= -G;
        }
    }

    /* Backward FFT: Phi(k) -> Phi(x) */
    fftw_execute(fft_backward);

    /* Copy result to Phi array */
    memcpy(Phi, fft_in, N3 * sizeof(double));
}

/* ==================================================================
 * Metric update: compute lapse from potential
 * PLAN.md Section 5: Safeguard A — Pade lapse
 *
 * alpha = 1 / sqrt(1 - 2*Phi)   (Pade, always > 0)
 * ================================================================== */
static void update_metric(void) {
    if (!gravity_on) return;

    /* Compute energy density at all points */
    double b = b_param;
    #pragma omp parallel for schedule(static) collapse(3)
    for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
    for (int k = 0; k < N; k++) {
        int p = idx(i, j, k);
        double b2 = b * b;
        double b4 = b2 * b2;

        double D2 = Dx[p]*Dx[p] + Dy[p]*Dy[p] + Dz[p]*Dz[p];
        double B2 = Bx[p]*Bx[p] + By[p]*By[p] + Bz[p]*Bz[p];

        double DxBx_l = Dy[p]*Bz[p] - Dz[p]*By[p];
        double DxBy_l = Dz[p]*Bx[p] - Dx[p]*Bz[p];
        double DxBz_l = Dx[p]*By[p] - Dy[p]*Bx[p];
        double P2 = DxBx_l*DxBx_l + DxBy_l*DxBy_l + DxBz_l*DxBz_l;

        double R = 1.0 + (D2 + B2) / b2 + P2 / b4;
        rho_BI[p] = b2 * (sqrt(R) - 1.0);
    }

    /* Solve Yukawa equation in Fourier space */
    solve_yukawa_fft();

    /* Compute Pade lapse: alpha = 1/sqrt(1 - 2*Phi) */
    #pragma omp parallel for schedule(static)
    for (int p = 0; p < N3; p++) {
        double arg = 1.0 - 2.0 * Phi[p];
        if (arg < 0.01) arg = 0.01;  /* safety floor */
        alpha_lapse[p] = 1.0 / sqrt(arg);
    }
}

/* ==================================================================
 * Sponge absorbing layer (from V11)
 * Applied to D and B in outer 15% of domain
 * ================================================================== */
static void sponge_absorb(void) {
    double sponge_width = 0.15 * L;
    double sigma_max = 5.0;

    #pragma omp parallel for schedule(static) collapse(3)
    for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
    for (int k = 0; k < N; k++) {
        double x = grid_x(i);
        double y = grid_y(j);
        double z = grid_z(k);

        double dx_b = fmin(L - fabs(x), fmin(L - fabs(y), L - fabs(z)));
        if (dx_b >= sponge_width) continue;

        double s = (sponge_width - dx_b) / sponge_width;
        double damp = exp(-sigma_max * s * s * dt);

        int p = idx(i, j, k);
        Dx[p] *= damp;  Dy[p] *= damp;  Dz[p] *= damp;
        Bx[p] *= damp;  By[p] *= damp;  Bz[p] *= damp;
    }
}

/* ==================================================================
 * Compute RHS: dD/dt = curl(alpha*H), dB/dt = -curl(alpha*E)
 * PLAN.md Section 9: ADM Evolution Equations
 *
 * CRITICAL: alpha goes INSIDE the curl, not outside.
 * curl(alpha*E) = alpha*curl(E) + grad(alpha) x E
 * The grad(alpha) x E term is the gravitational lensing force.
 * ================================================================== */
static void compute_rhs(double b,
    double *d_x, double *d_y, double *d_z,
    double *b_x, double *b_y, double *b_z,
    double *rhs_dx, double *rhs_dy, double *rhs_dz,
    double *rhs_bx, double *rhs_by, double *rhs_bz)
{
    /* Step 1: Compute E,H from (D,B) analytically, then form alpha*E, alpha*H */
    #pragma omp parallel for schedule(static) collapse(3)
    for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
    for (int k = 0; k < N; k++) {
        int p = idx(i, j, k);
        double ex, ey, ez, hx_l, hy_l, hz_l, rho_l;

        bi_constitutive_Dfield(
            d_x[p], d_y[p], d_z[p],
            b_x[p], b_y[p], b_z[p], b,
            &ex, &ey, &ez,
            &hx_l, &hy_l, &hz_l,
            &rho_l);

        Ex[p] = ex;  Ey[p] = ey;  Ez[p] = ez;
        Hfx[p] = hx_l;  Hfy[p] = hy_l;  Hfz[p] = hz_l;

        /* Form alpha*E and alpha*H products for curl */
        double al = alpha_lapse[p];
        aEx[p] = al * ex;  aEy[p] = al * ey;  aEz[p] = al * ez;
        aHx[p] = al * hx_l; aHy[p] = al * hy_l; aHz[p] = al * hz_l;
    }

    /* Step 2: curl(alpha*H) -> rhs_D, -curl(alpha*E) -> rhs_B */
    #pragma omp parallel for schedule(static) collapse(3)
    for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
    for (int k = 0; k < N; k++) {
        int p = idx(i, j, k);

        double ch_x, ch_y, ch_z;
        curl_at(aHx, aHy, aHz, i, j, k, &ch_x, &ch_y, &ch_z);

        double ce_x, ce_y, ce_z;
        curl_at(aEx, aEy, aEz, i, j, k, &ce_x, &ce_y, &ce_z);

        rhs_dx[p] =  ch_x;   /* dD/dt = +curl(alpha*H) */
        rhs_dy[p] =  ch_y;
        rhs_dz[p] =  ch_z;
        rhs_bx[p] = -ce_x;   /* dB/dt = -curl(alpha*E) */
        rhs_by[p] = -ce_y;
        rhs_bz[p] = -ce_z;
    }
}

/* ==================================================================
 * RK4 step (adapted from V11)
 * ================================================================== */
static void rk4_step(double b) {
    /* Save state backup */
    memcpy(Dx0, Dx, N3 * sizeof(double));
    memcpy(Dy0, Dy, N3 * sizeof(double));
    memcpy(Dz0, Dz, N3 * sizeof(double));
    memcpy(Bx0, Bx, N3 * sizeof(double));
    memcpy(By0, By, N3 * sizeof(double));
    memcpy(Bz0, Bz, N3 * sizeof(double));

    /* k1 */
    compute_rhs(b, Dx, Dy, Dz, Bx, By, Bz,
                rdx, rdy, rdz, rbx, rby, rbz);
    #pragma omp parallel for schedule(static)
    for (int p = 0; p < N3; p++) {
        kDx[p] = rdx[p] / 6.0;  kDy[p] = rdy[p] / 6.0;  kDz[p] = rdz[p] / 6.0;
        kBx[p] = rbx[p] / 6.0;  kBy[p] = rby[p] / 6.0;  kBz[p] = rbz[p] / 6.0;
    }
    #pragma omp parallel for schedule(static)
    for (int p = 0; p < N3; p++) {
        Dx[p] = Dx0[p] + 0.5*dt*rdx[p];  Dy[p] = Dy0[p] + 0.5*dt*rdy[p];
        Dz[p] = Dz0[p] + 0.5*dt*rdz[p];
        Bx[p] = Bx0[p] + 0.5*dt*rbx[p];  By[p] = By0[p] + 0.5*dt*rby[p];
        Bz[p] = Bz0[p] + 0.5*dt*rbz[p];
    }

    /* k2 */
    compute_rhs(b, Dx, Dy, Dz, Bx, By, Bz,
                rdx, rdy, rdz, rbx, rby, rbz);
    #pragma omp parallel for schedule(static)
    for (int p = 0; p < N3; p++) {
        kDx[p] += rdx[p] / 3.0;  kDy[p] += rdy[p] / 3.0;  kDz[p] += rdz[p] / 3.0;
        kBx[p] += rbx[p] / 3.0;  kBy[p] += rby[p] / 3.0;  kBz[p] += rbz[p] / 3.0;
    }
    #pragma omp parallel for schedule(static)
    for (int p = 0; p < N3; p++) {
        Dx[p] = Dx0[p] + 0.5*dt*rdx[p];  Dy[p] = Dy0[p] + 0.5*dt*rdy[p];
        Dz[p] = Dz0[p] + 0.5*dt*rdz[p];
        Bx[p] = Bx0[p] + 0.5*dt*rbx[p];  By[p] = By0[p] + 0.5*dt*rby[p];
        Bz[p] = Bz0[p] + 0.5*dt*rbz[p];
    }

    /* k3 */
    compute_rhs(b, Dx, Dy, Dz, Bx, By, Bz,
                rdx, rdy, rdz, rbx, rby, rbz);
    #pragma omp parallel for schedule(static)
    for (int p = 0; p < N3; p++) {
        kDx[p] += rdx[p] / 3.0;  kDy[p] += rdy[p] / 3.0;  kDz[p] += rdz[p] / 3.0;
        kBx[p] += rbx[p] / 3.0;  kBy[p] += rby[p] / 3.0;  kBz[p] += rbz[p] / 3.0;
    }
    #pragma omp parallel for schedule(static)
    for (int p = 0; p < N3; p++) {
        Dx[p] = Dx0[p] + dt*rdx[p];  Dy[p] = Dy0[p] + dt*rdy[p];
        Dz[p] = Dz0[p] + dt*rdz[p];
        Bx[p] = Bx0[p] + dt*rbx[p];  By[p] = By0[p] + dt*rby[p];
        Bz[p] = Bz0[p] + dt*rbz[p];
    }

    /* k4 */
    compute_rhs(b, Dx, Dy, Dz, Bx, By, Bz,
                rdx, rdy, rdz, rbx, rby, rbz);
    #pragma omp parallel for schedule(static)
    for (int p = 0; p < N3; p++) {
        kDx[p] += rdx[p] / 6.0;  kDy[p] += rdy[p] / 6.0;  kDz[p] += rdz[p] / 6.0;
        kBx[p] += rbx[p] / 6.0;  kBy[p] += rby[p] / 6.0;  kBz[p] += rbz[p] / 6.0;
        Dx[p] = Dx0[p] + dt*kDx[p];  Dy[p] = Dy0[p] + dt*kDy[p];
        Dz[p] = Dz0[p] + dt*kDz[p];
        Bx[p] = Bx0[p] + dt*kBx[p];  By[p] = By0[p] + dt*kBy[p];
        Bz[p] = Bz0[p] + dt*kBz[p];
    }

    /* Sponge layer */
    sponge_absorb();
}

/* ==================================================================
 * Diagnostics
 * PLAN.md Section 5: Key observables
 * ================================================================== */
static void compute_diagnostics(double b,
    double *E_total, double *R_eff, double *B_peak, double *E_peak,
    double *divB_max, double *divD_max, double *alpha_min, double *Phi_min)
{
    double Etot = 0.0, R2w = 0.0, Ew = 0.0;
    double Bpk = 0.0, Epk = 0.0;
    double divB_mx = 0.0, divD_mx = 0.0;
    double al_min = 1e30, phi_min = 1e30;
    double b2 = b * b;
    double b4 = b2 * b2;
    double hv = h * h * h;

    #pragma omp parallel for schedule(static) collapse(3) \
        reduction(+:Etot,R2w,Ew) reduction(max:Bpk,Epk,divB_mx,divD_mx) \
        reduction(min:al_min,phi_min)
    for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
    for (int k = 0; k < N; k++) {
        int p = idx(i, j, k);

        double D2 = Dx[p]*Dx[p] + Dy[p]*Dy[p] + Dz[p]*Dz[p];
        double B2 = Bx[p]*Bx[p] + By[p]*By[p] + Bz[p]*Bz[p];
        double DxBx_l = Dy[p]*Bz[p] - Dz[p]*By[p];
        double DxBy_l = Dz[p]*Bx[p] - Dx[p]*Bz[p];
        double DxBz_l = Dx[p]*By[p] - Dy[p]*Bx[p];
        double P2 = DxBx_l*DxBx_l + DxBy_l*DxBy_l + DxBz_l*DxBz_l;

        double R = 1.0 + (D2 + B2) / b2 + P2 / b4;
        double u = b2 * (sqrt(R) - 1.0);

        Etot += u * hv;

        double x = grid_x(i);
        double y = grid_y(j);
        double z = grid_z(k);
        double r2 = x*x + y*y + z*z;
        R2w += r2 * u * hv;
        Ew  += u * hv;

        if (sqrt(B2) > Bpk) Bpk = sqrt(B2);
        if (sqrt(D2) > Epk) Epk = sqrt(D2);  /* actually |D|, not |E| */

        /* Divergence constraints */
        double dBx_dx = deriv4(Bx, i, j, k, 0);
        double dBy_dy = deriv4(By, i, j, k, 1);
        double dBz_dz = deriv4(Bz, i, j, k, 2);
        double divB = fabs(dBx_dx + dBy_dy + dBz_dz);
        if (divB > divB_mx) divB_mx = divB;

        double dDx_dx = deriv4(Dx, i, j, k, 0);
        double dDy_dy = deriv4(Dy, i, j, k, 1);
        double dDz_dz = deriv4(Dz, i, j, k, 2);
        double divD = fabs(dDx_dx + dDy_dy + dDz_dz);
        if (divD > divD_mx) divD_mx = divD;

        /* Metric extrema */
        if (alpha_lapse[p] < al_min) al_min = alpha_lapse[p];
        if (Phi[p] < phi_min) phi_min = Phi[p];
    }

    *E_total  = Etot;
    *R_eff    = (Ew > 1e-30) ? sqrt(R2w / Ew) : 0.0;
    *B_peak   = Bpk;
    *E_peak   = Epk;
    *divB_max = divB_mx;
    *divD_max = divD_mx;
    *alpha_min = al_min;
    *Phi_min   = phi_min;
}

/* ==================================================================
 * Main evolution loop
 * ================================================================== */
static void run_evolve(double b) {
    printf("\n# === V12 BI-GEON: b=%.4f, a=%.2f, kappa=%.2f, mu=%.4f ===\n",
           b, a_param, kappa_param, mu_param);
    printf("# === N=%d, L=%.1f, T=%.1f, gravity=%s, metric_update=%d steps ===\n",
           N, L, T_final, gravity_on ? "ON" : "OFF", metric_interval);

    /* Initialize Hopfion B field */
    init_hopfion_B();

    /* Initial metric solve */
    if (gravity_on) {
        update_metric();
    }

    /* Initial diagnostics */
    double Etot, Reff, Bpk, Dpk, divBmax, divDmax, al_min, phi_min;
    compute_diagnostics(b, &Etot, &Reff, &Bpk, &Dpk, &divBmax, &divDmax,
                        &al_min, &phi_min);
    double E0 = Etot;

    printf("#\n# t         E_total      dE/E0       R_eff    B_peak     D_peak"
           "     divB_max   divD_max   alpha_min  Phi_min\n");
    printf("  %-9.4f  %-12.6f %-+11.3e %-8.4f %-10.6f %-10.6f %-10.2e %-10.2e %-10.6f %-10.6f\n",
           0.0, Etot, 0.0, Reff, Bpk, Dpk, divBmax, divDmax, al_min, phi_min);

    /* Open data file */
    char fname[256];
    snprintf(fname, sizeof(fname), "data/evolve_b%.2f_k%.1f_mu%.2f_%s.dat",
             b, kappa_param, mu_param, gravity_on ? "grav" : "flat");
    FILE *fp = fopen(fname, "w");
    if (fp) {
        fprintf(fp, "# V12 BI-Geon evolution\n");
        fprintf(fp, "# b=%.4f a=%.4f kappa=%.4f mu=%.4f N=%d L=%.1f gravity=%d\n",
                b, a_param, kappa_param, mu_param, N, L, gravity_on);
        fprintf(fp, "# t  E_total  dE/E0  R_eff  B_peak  D_peak  divB_max  divD_max  alpha_min  Phi_min\n");
        fprintf(fp, "%.6f  %.10f  %.6e  %.6f  %.10f  %.10f  %.4e  %.4e  %.6f  %.6f\n",
                0.0, Etot, 0.0, Reff, Bpk, Dpk, divBmax, divDmax, al_min, phi_min);
    }

    int nsteps = (int)(T_final / dt + 0.5);
    int out_interval = nsteps / 50;
    if (out_interval < 1) out_interval = 1;

    for (int step = 1; step <= nsteps; step++) {
        /* Update metric periodically */
        if (gravity_on && step % metric_interval == 0) {
            update_metric();
        }

        rk4_step(b);
        double t = step * dt;

        if (step % out_interval == 0 || step == nsteps) {
            compute_diagnostics(b, &Etot, &Reff, &Bpk, &Dpk, &divBmax, &divDmax,
                                &al_min, &phi_min);
            double dE = (Etot - E0) / E0;
            printf("  %-9.4f  %-12.6f %-+11.3e %-8.4f %-10.6f %-10.6f %-10.2e %-10.2e %-10.6f %-10.6f\n",
                   t, Etot, dE, Reff, Bpk, Dpk, divBmax, divDmax, al_min, phi_min);
            if (fp) {
                fprintf(fp, "%.6f  %.10f  %.6e  %.6f  %.10f  %.10f  %.4e  %.4e  %.6f  %.6f\n",
                        t, Etot, dE, Reff, Bpk, Dpk, divBmax, divDmax, al_min, phi_min);
            }

            /* Early termination checks */
            if (al_min < 0.01) {
                printf("# WARNING: alpha_min < 0.01 — near horizon, stopping.\n");
                break;
            }
            if (fabs(dE) > 0.5) {
                printf("# WARNING: energy conservation violated (dE/E=%.3f), stopping.\n", dE);
                break;
            }
        }
    }

    if (fp) fclose(fp);
    printf("#\n# Output: %s\n", fname);

    /* Final summary */
    printf("#\n# === FINAL SUMMARY ===\n");
    printf("# E_total = %.6f (initial: %.6f, dE/E0 = %.4e)\n", Etot, E0, (Etot-E0)/E0);
    printf("# R_eff   = %.4f\n", Reff);
    printf("# B_peak  = %.6f\n", Bpk);
    printf("# D_peak  = %.6f\n", Dpk);
    printf("# divB    = %.4e\n", divBmax);
    printf("# divD    = %.4e\n", divDmax);
    if (gravity_on) {
        printf("# alpha_min = %.6f (Phi_min = %.6f)\n", al_min, phi_min);
        printf("# well_depth = %.4f\n", 1.0 - 1.0/(al_min*al_min));
    }
    printf("# Stable:  R_eff %s (initial R ~ %.2f)\n",
           Reff < 1.5 * Reff ? "YES" : "NO", Reff);
}

/* ==================================================================
 * Main
 * ================================================================== */
int main(int argc, char **argv) {
    int mode = 0;

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-evolve") == 0)      mode = 1;
        else if (strcmp(argv[i], "-flat") == 0)  { mode = 1; gravity_on = 0; }
        else if (strcmp(argv[i], "-b") == 0 && i+1 < argc) b_param = atof(argv[++i]);
        else if (strcmp(argv[i], "-a") == 0 && i+1 < argc) a_param = atof(argv[++i]);
        else if (strcmp(argv[i], "-kappa") == 0 && i+1 < argc) kappa_param = atof(argv[++i]);
        else if (strcmp(argv[i], "-mu") == 0 && i+1 < argc) mu_param = atof(argv[++i]);
        else if (strcmp(argv[i], "-N") == 0 && i+1 < argc) N = atoi(argv[++i]);
        else if (strcmp(argv[i], "-L") == 0 && i+1 < argc) L = atof(argv[++i]);
        else if (strcmp(argv[i], "-T") == 0 && i+1 < argc) T_final = atof(argv[++i]);
        else if (strcmp(argv[i], "-metric_interval") == 0 && i+1 < argc)
            metric_interval = atoi(argv[++i]);
        else {
            fprintf(stderr, "Unknown arg: %s\n", argv[i]);
            fprintf(stderr, "Usage: bi_geon [-evolve|-flat] [-b B] [-a A] [-kappa K] [-mu M]\n");
            fprintf(stderr, "               [-N N] [-L L] [-T T] [-metric_interval M]\n");
            return 1;
        }
    }

    if (mode == 0) {
        fprintf(stderr, "Must specify -evolve or -flat\n");
        return 1;
    }

    fftw_init_threads();
    fftw_plan_with_nthreads(omp_get_max_threads());

    setup_grid();
    alloc_fields();

    run_evolve(b_param);

    free_fields();
    fftw_cleanup_threads();
    return 0;
}
