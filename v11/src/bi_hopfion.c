/*
 * bi_hopfion.c — Born-Infeld magnetic hopfion evolution
 *
 * Tests whether BI nonlinearity stabilizes EM configurations with Hopf topology
 * in flat space. Pure magnetic initial data from Hopf map pullback.
 *
 * Time integration: RK4 method-of-lines on state (D, B).
 * E derived from (D, B) via Newton inversion at each RHS evaluation.
 *
 * Usage:
 *   bi_hopfion -evolve [-b 1.0] [-a 1.0] [-N 128] [-L 6.0] [-T 10.0]
 *   bi_hopfion -bscan  [-a 1.0] [-N 128] [-L 6.0] [-T 10.0]
 *   bi_hopfion -maxwell [-a 1.0] [-N 128] [-L 6.0] [-T 10.0]
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>

/* ---- Parameters ---- */
static double b_param = 1.0;    /* BI field strength */
static double a_param = 1.0;    /* hopfion size */
static int    N       = 128;    /* grid points per dim */
static double L       = 6.0;    /* box half-width */
static double T_final = 10.0;   /* total time */

/* ---- Derived ---- */
static double h, dt;
static int    N3;

/* ---- State: Dx,Dy,Dz, Bx,By,Bz ---- */
static double *Dx, *Dy, *Dz;
static double *Bx, *By, *Bz;

/* ---- Working arrays: Ex,Ey,Ez, Hx,Hy,Hz ---- */
static double *Ex, *Ey, *Ez;
static double *Hx, *Hy, *Hz;

/* ---- RK4 scratch: backup + accumulator + RHS temp ---- */
static double *Dx0, *Dy0, *Dz0, *Bx0, *By0, *Bz0;  /* state backup */
static double *kDx, *kDy, *kDz, *kBx, *kBy, *kBz;  /* RHS accumulator */
static double *rdx, *rdy, *rdz, *rbx, *rby, *rbz;   /* temp RHS */

/* ---- Grid ---- */
static inline int idx(int i, int j, int k) {
    /* periodic BC */
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
    dt = 0.4 * h / sqrt(3.0);  /* CFL for RK4 + 4th-order spatial */
    printf("Grid: N=%d, L=%.1f, h=%.4f, dt=%.6f, N3=%d\n", N, L, h, dt, N3);
    printf("Memory: %.1f MB\n", 30.0 * N3 * 8.0 / 1e6);
}

static void alloc_fields(void) {
    Ex = calloc(N3, sizeof(double));  Ey = calloc(N3, sizeof(double));
    Ez = calloc(N3, sizeof(double));
    Bx = calloc(N3, sizeof(double));  By = calloc(N3, sizeof(double));
    Bz = calloc(N3, sizeof(double));
    Dx = calloc(N3, sizeof(double));  Dy = calloc(N3, sizeof(double));
    Dz = calloc(N3, sizeof(double));
    Hx = calloc(N3, sizeof(double));  Hy = calloc(N3, sizeof(double));
    Hz = calloc(N3, sizeof(double));
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
    if (!Ex || !Ey || !Ez || !Bx || !By || !Bz ||
        !Dx || !Dy || !Dz || !Hx || !Hy || !Hz ||
        !Dx0 || !Dy0 || !Dz0 || !Bx0 || !By0 || !Bz0 ||
        !kDx || !kDy || !kDz || !kBx || !kBy || !kBz ||
        !rdx || !rdy || !rdz || !rbx || !rby || !rbz) {
        fprintf(stderr, "Allocation failed\n");
        exit(1);
    }
}

static void free_fields(void) {
    free(Ex); free(Ey); free(Ez);
    free(Bx); free(By); free(Bz);
    free(Dx); free(Dy); free(Dz);
    free(Hx); free(Hy); free(Hz);
    free(Dx0); free(Dy0); free(Dz0);
    free(Bx0); free(By0); free(Bz0);
    free(kDx); free(kDy); free(kDz);
    free(kBx); free(kBy); free(kBz);
    free(rdx); free(rdy); free(rdz);
    free(rbx); free(rby); free(rbz);
}

/* ---- Hopf map initial data (pure magnetic, E=D=0) ---- */
static void init_hopfion_B(void) {
    double a = a_param;
    double a2 = a * a;
    /* C0 chosen so peak |B| = b_param */
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

        /* Partial derivatives of (u1,v1,u2,v2) w.r.t. (x,y,z) */
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

        /* n = CP1 Hopf map: S3 -> S2 */
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
        /* Cross products for each B component */
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

        /* E = D = 0 initially */
        Ex[p] = Ey[p] = Ez[p] = 0.0;
        Dx[p] = Dy[p] = Dz[p] = 0.0;
    }
}

/* ---- Born-Infeld constitutive relations ---- */
/* Forward: (E,B) -> (D,H) */
static inline void bi_constitutive(
    double ex, double ey, double ez,
    double bx, double by, double bz,
    double b,
    double *dx, double *dy, double *dz,
    double *hx, double *hy, double *hz)
{
    double B2 = bx*bx + by*by + bz*bz;
    double E2 = ex*ex + ey*ey + ez*ez;
    double EdB = ex*bx + ey*by + ez*bz;
    double b2 = b * b;

    double arg = 1.0 + (B2 - E2) / b2 - (EdB * EdB) / (b2 * b2);
    if (arg < 1e-12) arg = 1e-12;  /* safety clamp */
    double Gamma = sqrt(arg);
    double invG = 1.0 / Gamma;

    /* D = (E + (E.B)B/b²) / Gamma */
    *dx = (ex + EdB * bx / b2) * invG;
    *dy = (ey + EdB * by / b2) * invG;
    *dz = (ez + EdB * bz / b2) * invG;

    /* H = (B - (E.B)E/b²) / Gamma */
    *hx = (bx - EdB * ex / b2) * invG;
    *hy = (by - EdB * ey / b2) * invG;
    *hz = (bz - EdB * ez / b2) * invG;
}

/* ---- Newton inversion: (D_target, B) -> E ---- */
static inline void newton_invert(
    double Dt_x, double Dt_y, double Dt_z,
    double bx, double by, double bz,
    double b,
    double *ex_out, double *ey_out, double *ez_out)
{
    /* Initial guess: E = D_target */
    double ex = Dt_x, ey = Dt_y, ez = Dt_z;
    double b2 = b * b;
    double b4 = b2 * b2;

    for (int iter = 0; iter < 8; iter++) {
        double B2 = bx*bx + by*by + bz*bz;
        double E2 = ex*ex + ey*ey + ez*ez;
        double EdB = ex*bx + ey*by + ez*bz;

        double arg = 1.0 + (B2 - E2) / b2 - (EdB * EdB) / b4;
        if (arg < 1e-12) arg = 1e-12;
        double Gamma = sqrt(arg);
        double invG = 1.0 / Gamma;

        double dx = (ex + EdB * bx / b2) * invG;
        double dy = (ey + EdB * by / b2) * invG;
        double dz = (ez + EdB * bz / b2) * invG;

        /* Residual G = D(E,B) - D_target */
        double gx = dx - Dt_x;
        double gy = dy - Dt_y;
        double gz = dz - Dt_z;

        double res2 = gx*gx + gy*gy + gz*gz;
        if (res2 < 1e-28) break;

        /* Jacobian J_{ij} = [delta_{ij} + (B_iB_j + D_iD_j)/b²] / Gamma */
        double J[3][3];
        double Bv[3] = {bx, by, bz};
        double Dv[3] = {dx, dy, dz};

        for (int m = 0; m < 3; m++)
        for (int n = 0; n < 3; n++) {
            J[m][n] = ((m == n ? 1.0 : 0.0) + (Bv[m]*Bv[n] + Dv[m]*Dv[n]) / b2) * invG;
        }

        /* Solve J * dE = -G by Gaussian elimination (3x3) */
        double rhs[3] = {-gx, -gy, -gz};

        for (int row = 0; row < 3; row++) {
            int piv = row;
            double pmax = fabs(J[row][row]);
            for (int r2 = row + 1; r2 < 3; r2++) {
                if (fabs(J[r2][row]) > pmax) {
                    pmax = fabs(J[r2][row]);
                    piv = r2;
                }
            }
            if (piv != row) {
                double tmp;
                for (int c = 0; c < 3; c++) {
                    tmp = J[row][c]; J[row][c] = J[piv][c]; J[piv][c] = tmp;
                }
                tmp = rhs[row]; rhs[row] = rhs[piv]; rhs[piv] = tmp;
            }

            double diag = J[row][row];
            if (fabs(diag) < 1e-30) break;

            for (int r2 = row + 1; r2 < 3; r2++) {
                double fac = J[r2][row] / diag;
                for (int c = row + 1; c < 3; c++)
                    J[r2][c] -= fac * J[row][c];
                rhs[r2] -= fac * rhs[row];
            }
        }

        double dE[3];
        for (int row = 2; row >= 0; row--) {
            double sum = rhs[row];
            for (int c = row + 1; c < 3; c++)
                sum -= J[row][c] * dE[c];
            dE[row] = sum / J[row][row];
        }

        ex += dE[0];
        ey += dE[1];
        ez += dE[2];
    }

    *ex_out = ex;
    *ey_out = ey;
    *ez_out = ez;
}

/* ---- 4th-order central derivative ---- */
static inline double deriv4(const double *f, int i, int j, int k, int dim) {
    int di = (dim == 0), dj = (dim == 1), dk = (dim == 2);
    double fm2 = f[idx(i - 2*di, j - 2*dj, k - 2*dk)];
    double fm1 = f[idx(i - di,   j - dj,   k - dk)];
    double fp1 = f[idx(i + di,   j + dj,   k + dk)];
    double fp2 = f[idx(i + 2*di, j + 2*dj, k + 2*dk)];
    return (-fp2 + 8.0*fp1 - 8.0*fm1 + fm2) / (12.0 * h);
}

/* ---- Curl at a point using 4th-order derivatives ---- */
static inline void curl_at(
    const double *fx, const double *fy, const double *fz,
    int i, int j, int k,
    double *cx, double *cy, double *cz)
{
    *cx = deriv4(fz, i, j, k, 1) - deriv4(fy, i, j, k, 2);
    *cy = deriv4(fx, i, j, k, 2) - deriv4(fz, i, j, k, 0);
    *cz = deriv4(fy, i, j, k, 0) - deriv4(fx, i, j, k, 1);
}

/* ---- Sponge absorbing layer (applied to state D, B) ---- */
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

/* ---- Compute RHS: given state (D, B), compute dD/dt and dB/dt ---- */
/*
 * BI equations: dD/dt = curl(H),  dB/dt = -curl(E)
 * Steps: 1. Invert (D,B) -> E at all points
 *        2. Compute H(E,B) at all points
 *        3. Compute curl(H) -> rhs_D, -curl(E) -> rhs_B
 *
 * rhs_D[x,y,z] and rhs_B[x,y,z] are written to the 6 output arrays.
 */
static void compute_rhs(double b,
    double *d_x, double *d_y, double *d_z,
    double *b_x, double *b_y, double *b_z,
    double *rhs_dx, double *rhs_dy, double *rhs_dz,
    double *rhs_bx, double *rhs_by, double *rhs_bz)
{
    /* Step 1: Newton-invert (D,B) -> E, then compute H(E,B) */
    #pragma omp parallel for schedule(static) collapse(3)
    for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
    for (int k = 0; k < N; k++) {
        int p = idx(i, j, k);
        double ex, ey, ez;
        newton_invert(d_x[p], d_y[p], d_z[p],
                      b_x[p], b_y[p], b_z[p], b,
                      &ex, &ey, &ez);
        Ex[p] = ex;
        Ey[p] = ey;
        Ez[p] = ez;

        double hx_l, hy_l, hz_l, dx_l, dy_l, dz_l;
        bi_constitutive(ex, ey, ez,
                        b_x[p], b_y[p], b_z[p], b,
                        &dx_l, &dy_l, &dz_l,
                        &hx_l, &hy_l, &hz_l);
        Hx[p] = hx_l;
        Hy[p] = hy_l;
        Hz[p] = hz_l;
    }

    /* Step 2: compute curl(H) -> rhs_D, -curl(E) -> rhs_B */
    #pragma omp parallel for schedule(static) collapse(3)
    for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
    for (int k = 0; k < N; k++) {
        int p = idx(i, j, k);

        double ch_x, ch_y, ch_z;
        curl_at(Hx, Hy, Hz, i, j, k, &ch_x, &ch_y, &ch_z);

        double ce_x, ce_y, ce_z;
        curl_at(Ex, Ey, Ez, i, j, k, &ce_x, &ce_y, &ce_z);

        rhs_dx[p] = ch_x;   /* dD/dt = curl(H) */
        rhs_dy[p] = ch_y;
        rhs_dz[p] = ch_z;
        rhs_bx[p] = -ce_x;  /* dB/dt = -curl(E) */
        rhs_by[p] = -ce_y;
        rhs_bz[p] = -ce_z;
    }
}

/* ---- RK4 step ---- */
/*
 * State: (Dx,Dy,Dz,Bx,By,Bz) at time n
 * Uses scratch: (Dx0,...,Bz0) for backup, (kDx,...,kBz) for RHS & accumulator
 *               (Ex,Ey,Ez,Hx,Hy,Hz) as working arrays
 */
static void rk4_step(double b) {
    /* Save state backup */
    memcpy(Dx0, Dx, N3 * sizeof(double));
    memcpy(Dy0, Dy, N3 * sizeof(double));
    memcpy(Dz0, Dz, N3 * sizeof(double));
    memcpy(Bx0, Bx, N3 * sizeof(double));
    memcpy(By0, By, N3 * sizeof(double));
    memcpy(Bz0, Bz, N3 * sizeof(double));

    /* ---- k1 = f(y_n) ---- */
    compute_rhs(b, Dx, Dy, Dz, Bx, By, Bz,
                rdx, rdy, rdz, rbx, rby, rbz);

    /* accumulator = k1/6 */
    #pragma omp parallel for schedule(static)
    for (int p = 0; p < N3; p++) {
        kDx[p] = rdx[p] / 6.0;  kDy[p] = rdy[p] / 6.0;  kDz[p] = rdz[p] / 6.0;
        kBx[p] = rbx[p] / 6.0;  kBy[p] = rby[p] / 6.0;  kBz[p] = rbz[p] / 6.0;
    }

    /* y_1 = y_n + dt/2 * k1 */
    #pragma omp parallel for schedule(static)
    for (int p = 0; p < N3; p++) {
        Dx[p] = Dx0[p] + 0.5 * dt * rdx[p];
        Dy[p] = Dy0[p] + 0.5 * dt * rdy[p];
        Dz[p] = Dz0[p] + 0.5 * dt * rdz[p];
        Bx[p] = Bx0[p] + 0.5 * dt * rbx[p];
        By[p] = By0[p] + 0.5 * dt * rby[p];
        Bz[p] = Bz0[p] + 0.5 * dt * rbz[p];
    }

    /* ---- k2 = f(y_1) ---- */
    compute_rhs(b, Dx, Dy, Dz, Bx, By, Bz,
                rdx, rdy, rdz, rbx, rby, rbz);

    /* accumulator += k2/3 */
    #pragma omp parallel for schedule(static)
    for (int p = 0; p < N3; p++) {
        kDx[p] += rdx[p] / 3.0;  kDy[p] += rdy[p] / 3.0;  kDz[p] += rdz[p] / 3.0;
        kBx[p] += rbx[p] / 3.0;  kBy[p] += rby[p] / 3.0;  kBz[p] += rbz[p] / 3.0;
    }

    /* y_2 = y_n + dt/2 * k2 */
    #pragma omp parallel for schedule(static)
    for (int p = 0; p < N3; p++) {
        Dx[p] = Dx0[p] + 0.5 * dt * rdx[p];
        Dy[p] = Dy0[p] + 0.5 * dt * rdy[p];
        Dz[p] = Dz0[p] + 0.5 * dt * rdz[p];
        Bx[p] = Bx0[p] + 0.5 * dt * rbx[p];
        By[p] = By0[p] + 0.5 * dt * rby[p];
        Bz[p] = Bz0[p] + 0.5 * dt * rbz[p];
    }

    /* ---- k3 = f(y_2) ---- */
    compute_rhs(b, Dx, Dy, Dz, Bx, By, Bz,
                rdx, rdy, rdz, rbx, rby, rbz);

    /* accumulator += k3/3 */
    #pragma omp parallel for schedule(static)
    for (int p = 0; p < N3; p++) {
        kDx[p] += rdx[p] / 3.0;  kDy[p] += rdy[p] / 3.0;  kDz[p] += rdz[p] / 3.0;
        kBx[p] += rbx[p] / 3.0;  kBy[p] += rby[p] / 3.0;  kBz[p] += rbz[p] / 3.0;
    }

    /* y_3 = y_n + dt * k3 */
    #pragma omp parallel for schedule(static)
    for (int p = 0; p < N3; p++) {
        Dx[p] = Dx0[p] + dt * rdx[p];
        Dy[p] = Dy0[p] + dt * rdy[p];
        Dz[p] = Dz0[p] + dt * rdz[p];
        Bx[p] = Bx0[p] + dt * rbx[p];
        By[p] = By0[p] + dt * rby[p];
        Bz[p] = Bz0[p] + dt * rbz[p];
    }

    /* ---- k4 = f(y_3) ---- */
    compute_rhs(b, Dx, Dy, Dz, Bx, By, Bz,
                rdx, rdy, rdz, rbx, rby, rbz);

    /* accumulator += k4/6;  y_{n+1} = y_n + dt * accumulator */
    #pragma omp parallel for schedule(static)
    for (int p = 0; p < N3; p++) {
        kDx[p] += rdx[p] / 6.0;  kDy[p] += rdy[p] / 6.0;  kDz[p] += rdz[p] / 6.0;
        kBx[p] += rbx[p] / 6.0;  kBy[p] += rby[p] / 6.0;  kBz[p] += rbz[p] / 6.0;

        Dx[p] = Dx0[p] + dt * kDx[p];
        Dy[p] = Dy0[p] + dt * kDy[p];
        Dz[p] = Dz0[p] + dt * kDz[p];
        Bx[p] = Bx0[p] + dt * kBx[p];
        By[p] = By0[p] + dt * kBy[p];
        Bz[p] = Bz0[p] + dt * kBz[p];
    }

    /* Recover E from final state for diagnostics */
    #pragma omp parallel for schedule(static) collapse(3)
    for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
    for (int k = 0; k < N; k++) {
        int p = idx(i, j, k);
        newton_invert(Dx[p], Dy[p], Dz[p],
                      Bx[p], By[p], Bz[p], b,
                      &Ex[p], &Ey[p], &Ez[p]);
    }

    /* Sponge layer */
    sponge_absorb();
}

/* ---- Diagnostics ---- */
static void compute_diagnostics(double b, double *E_total, double *R_eff,
                                 double *B_peak, double *E_peak,
                                 double *divB_max, double *divD_max)
{
    double Etot = 0.0, R2w = 0.0, Ew = 0.0;
    double Bpk = 0.0, Epk = 0.0;
    double divB_mx = 0.0, divD_mx = 0.0;
    double b2 = b * b;
    double hv = h * h * h;

    #pragma omp parallel for schedule(static) collapse(3) \
        reduction(+:Etot,R2w,Ew) reduction(max:Bpk,Epk,divB_mx,divD_mx)
    for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
    for (int k = 0; k < N; k++) {
        int p = idx(i, j, k);
        double ex = Ex[p], ey = Ey[p], ez = Ez[p];
        double bx_f = Bx[p], by_f = By[p], bz_f = Bz[p];

        double B2 = bx_f*bx_f + by_f*by_f + bz_f*bz_f;
        double E2 = ex*ex + ey*ey + ez*ez;
        double EdB = ex*bx_f + ey*by_f + ez*bz_f;

        double arg = 1.0 + (B2 - E2) / b2 - (EdB * EdB) / (b2 * b2);
        if (arg < 1e-12) arg = 1e-12;
        double Gamma = sqrt(arg);

        /* Cancellation-free energy: u = (B²(1+B²/b²) + E² + (E·B)²/b²) / (Gamma*(1+B²/b²+Gamma)) */
        double Y = 1.0 + B2 / b2;
        double u = (B2 * Y + E2 + EdB * EdB / b2) / (Gamma * (Y + Gamma));

        Etot += u * hv;

        double x = grid_x(i);
        double y = grid_y(j);
        double z = grid_z(k);
        double r2 = x*x + y*y + z*z;
        R2w += r2 * u * hv;
        Ew  += u * hv;

        double Bmag = sqrt(B2);
        double Emag = sqrt(E2);
        if (Bmag > Bpk) Bpk = Bmag;
        if (Emag > Epk) Epk = Emag;

        /* div B */
        double dBx_dx = deriv4(Bx, i, j, k, 0);
        double dBy_dy = deriv4(By, i, j, k, 1);
        double dBz_dz = deriv4(Bz, i, j, k, 2);
        double divB = fabs(dBx_dx + dBy_dy + dBz_dz);
        if (divB > divB_mx) divB_mx = divB;

        /* div D */
        double dDx_dx = deriv4(Dx, i, j, k, 0);
        double dDy_dy = deriv4(Dy, i, j, k, 1);
        double dDz_dz = deriv4(Dz, i, j, k, 2);
        double divD = fabs(dDx_dx + dDy_dy + dDz_dz);
        if (divD > divD_mx) divD_mx = divD;
    }

    *E_total = Etot;
    *R_eff = (Ew > 1e-30) ? sqrt(R2w / Ew) : 0.0;
    *B_peak = Bpk;
    *E_peak = Epk;
    *divB_max = divB_mx;
    *divD_max = divD_mx;
}

/* ---- Run modes ---- */

static void run_evolve(double b) {
    printf("\n=== EVOLVE: b=%.4f, a=%.2f, N=%d, L=%.1f, T=%.1f ===\n",
           b, a_param, N, L, T_final);

    init_hopfion_B();

    double Etot, Reff, Bpk, Epk, divBmax, divDmax;
    compute_diagnostics(b, &Etot, &Reff, &Bpk, &Epk, &divBmax, &divDmax);
    printf("t=0.000: E=%.6f  R_eff=%.4f  B_peak=%.6f  E_peak=%.6f  divB=%.2e  divD=%.2e\n",
           Etot, Reff, Bpk, Epk, divBmax, divDmax);

    double E0 = Etot;
    int nsteps = (int)(T_final / dt + 0.5);
    int out_interval = nsteps / 50;
    if (out_interval < 1) out_interval = 1;

    char fname[256];
    snprintf(fname, sizeof(fname), "data/evolve_b%.4f.dat", b);
    FILE *fp = fopen(fname, "w");
    if (fp) {
        fprintf(fp, "# t  E_total  R_eff  B_peak  E_peak  divB_max  divD_max  dE/E0\n");
        fprintf(fp, "%.6f  %.10f  %.6f  %.10f  %.10f  %.4e  %.4e  %.4e\n",
                0.0, Etot, Reff, Bpk, Epk, divBmax, divDmax, 0.0);
    }

    for (int step = 1; step <= nsteps; step++) {
        rk4_step(b);
        double t = step * dt;

        if (step % out_interval == 0 || step == nsteps) {
            compute_diagnostics(b, &Etot, &Reff, &Bpk, &Epk, &divBmax, &divDmax);
            double dE = (Etot - E0) / E0;
            printf("t=%.3f: E=%.6f(%.2e)  R=%.4f  B=%.6f  E=%.6f  dB=%.2e  dD=%.2e\n",
                   t, Etot, dE, Reff, Bpk, Epk, divBmax, divDmax);
            if (fp) {
                fprintf(fp, "%.6f  %.10f  %.6f  %.10f  %.10f  %.4e  %.4e  %.4e\n",
                        t, Etot, Reff, Bpk, Epk, divBmax, divDmax, dE);
            }
        }
    }

    if (fp) fclose(fp);
    printf("Done. Output: %s\n", fname);
}

static void run_bscan(void) {
    printf("\n=== B-SCAN: a=%.2f, N=%d, L=%.1f, T=%.1f ===\n",
           a_param, N, L, T_final);

    double b_values[] = {1.0, 2.0, 5.0, 10.0, 100.0, 1e10};
    int nb = sizeof(b_values) / sizeof(b_values[0]);

    FILE *fp = fopen("data/bscan.dat", "w");
    if (fp) fprintf(fp, "# b  E_init  E_final  R_init  R_final  B_peak_init  B_peak_final\n");

    for (int ib = 0; ib < nb; ib++) {
        double b = b_values[ib];
        printf("\n--- b = %.4e ---\n", b);

        init_hopfion_B();

        double Etot, Reff, Bpk, Epk, divBmax, divDmax;
        compute_diagnostics(b, &Etot, &Reff, &Bpk, &Epk, &divBmax, &divDmax);
        double E_init = Etot, R_init = Reff, B_init = Bpk;
        printf("  t=0: E=%.6f  R=%.4f  Bpk=%.6f\n", Etot, Reff, Bpk);

        int nsteps = (int)(T_final / dt + 0.5);
        int out_interval = nsteps / 10;
        if (out_interval < 1) out_interval = 1;

        char fname[256];
        snprintf(fname, sizeof(fname), "data/bscan_b%.4e.dat", b);
        FILE *fp2 = fopen(fname, "w");
        if (fp2) {
            fprintf(fp2, "# t  E_total  R_eff  B_peak  E_peak\n");
            fprintf(fp2, "%.6f  %.10f  %.6f  %.10f  %.10f\n",
                    0.0, Etot, Reff, Bpk, Epk);
        }

        for (int step = 1; step <= nsteps; step++) {
            rk4_step(b);

            if (step % out_interval == 0 || step == nsteps) {
                compute_diagnostics(b, &Etot, &Reff, &Bpk, &Epk, &divBmax, &divDmax);
                double t = step * dt;
                printf("  t=%.2f: E=%.6f  R=%.4f  Bpk=%.6f\n",
                       t, Etot, Reff, Bpk);
                if (fp2) {
                    fprintf(fp2, "%.6f  %.10f  %.6f  %.10f  %.10f\n",
                            t, Etot, Reff, Bpk, Epk);
                }
            }
        }

        if (fp2) fclose(fp2);

        if (fp) {
            fprintf(fp, "%.6e  %.10f  %.10f  %.6f  %.6f  %.10f  %.10f\n",
                    b, E_init, Etot, R_init, Reff, B_init, Bpk);
        }
    }

    if (fp) fclose(fp);
    printf("\nB-scan summary: data/bscan.dat\n");
}

static void run_maxwell_test(void) {
    double b_maxwell = 1e10;
    printf("\n=== MAXWELL TEST (b=%.0e, should match free EM dispersal) ===\n", b_maxwell);

    init_hopfion_B();

    double Etot, Reff, Bpk, Epk, divBmax, divDmax;
    compute_diagnostics(b_maxwell, &Etot, &Reff, &Bpk, &Epk, &divBmax, &divDmax);
    double E0 = Etot, R0 = Reff;
    printf("t=0.000: E=%.10f  R_eff=%.6f  B_peak=%.10f  divB=%.4e\n",
           Etot, Reff, Bpk, divBmax);

    FILE *fp = fopen("data/maxwell_test.dat", "w");
    if (fp) {
        fprintf(fp, "# t  E_total  R_eff  B_peak  E_peak  divB_max  dE/E0\n");
        fprintf(fp, "%.6f  %.10f  %.6f  %.10f  %.10f  %.4e  %.4e\n",
                0.0, Etot, Reff, Bpk, Epk, divBmax, 0.0);
    }

    int nsteps = (int)(T_final / dt + 0.5);
    int out_interval = nsteps / 20;
    if (out_interval < 1) out_interval = 1;

    for (int step = 1; step <= nsteps; step++) {
        rk4_step(b_maxwell);
        double t = step * dt;

        if (step % out_interval == 0 || step == nsteps) {
            compute_diagnostics(b_maxwell, &Etot, &Reff, &Bpk, &Epk, &divBmax, &divDmax);
            double dE = (Etot - E0) / E0;
            printf("t=%.3f: E=%.10f(%.2e)  R=%.6f  Bpk=%.6f  divB=%.4e\n",
                   t, Etot, dE, Reff, Bpk, divBmax);
            if (fp) {
                fprintf(fp, "%.6f  %.10f  %.6f  %.10f  %.10f  %.4e  %.4e\n",
                        t, Etot, Reff, Bpk, Epk, divBmax, dE);
            }
        }
    }

    if (fp) fclose(fp);

    printf("\n--- Validation ---\n");
    printf("Energy conservation: dE/E0 = %.4e (target: < 1e-4)\n",
           fabs(Etot - E0) / E0);
    printf("R_eff growth: R_final/R_init = %.4f (expect growth for free dispersal)\n",
           Reff / R0);
    printf("div B: %.4e\n", divBmax);
    printf("Maxwell test output: data/maxwell_test.dat\n");
}

/* ---- Main ---- */
int main(int argc, char **argv) {
    int mode = 0;

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-evolve") == 0)      mode = 1;
        else if (strcmp(argv[i], "-bscan") == 0)  mode = 2;
        else if (strcmp(argv[i], "-maxwell") == 0) mode = 3;
        else if (strcmp(argv[i], "-b") == 0 && i+1 < argc) b_param = atof(argv[++i]);
        else if (strcmp(argv[i], "-a") == 0 && i+1 < argc) a_param = atof(argv[++i]);
        else if (strcmp(argv[i], "-N") == 0 && i+1 < argc) N = atoi(argv[++i]);
        else if (strcmp(argv[i], "-L") == 0 && i+1 < argc) L = atof(argv[++i]);
        else if (strcmp(argv[i], "-T") == 0 && i+1 < argc) T_final = atof(argv[++i]);
        else {
            fprintf(stderr, "Unknown arg: %s\n", argv[i]);
            fprintf(stderr, "Usage: bi_hopfion [-evolve|-bscan|-maxwell] [-b B] [-a A] [-N N] [-L L] [-T T]\n");
            return 1;
        }
    }

    if (mode == 0) {
        fprintf(stderr, "Must specify -evolve, -bscan, or -maxwell\n");
        return 1;
    }

    setup_grid();
    alloc_fields();

    switch (mode) {
        case 1: run_evolve(b_param); break;
        case 2: run_bscan(); break;
        case 3: run_maxwell_test(); break;
    }

    free_fields();
    return 0;
}
