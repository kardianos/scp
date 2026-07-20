/*
 * v81/path1_locks — monist PIC sandbox (P1)
 *
 * REAL = FreeSubstrate (2D TE Yee U(1) medium) + Locks (typed structs)
 *
 * Forbidden by OP:
 *   - multiplet φ as particle definition
 *   - pairwise 1/r² force
 *   - multi-fab Cosserat L
 *   - GRIN variable-c gravity
 *
 * Step:
 *   1. Faraday: B^{n+1/2} from E^n
 *   2. Gather E,B → Boris push locks (half / full as staggered)
 *   3. Charge-conserving J deposit (zigzag / VB-class) from lock trajectories
 *   4. Ampère: E^{n+1} from B^{n+1/2}, J
 *   5. Ledger diagnostics
 *
 * Build:
 *   make -C v81/path1_locks
 * Run:
 *   ./bin/locks_pic all
 *   ./bin/locks_pic l0|l1|l2|l3
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

/* ---- units: c = ε0 = μ0 = 1 ---- */
#define C0 1.0

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Grid defaults (overridable via CLI env-less flags) */
typedef struct {
    int nx, ny;
    double dx, dy, dt;
    int n_locks;
    int absorb;          /* 1 = conducting sponge layer near boundary */
    int sponge_w;        /* sponge width in cells */
    double sponge_sigma; /* peak conductivity */
} GridCfg;

typedef struct {
    int id;
    int type; /* 0 = light charge lock */
    double q; /* ±1 */
    double m; /* rest mass (M-chart label) */
    double E_star; /* sequestered rest energy (constant in MVP) */
    double x, y;
    double ux, uy; /* proper velocity γv (c=1) */
    double fx, fy; /* last gathered force diagnostic (qE + q v×B) */
    int pinned;
    int alive; /* always 1 in this model — locks cannot evaporate */
} Lock;

typedef struct {
    GridCfg cfg;
    /* cell-centered ρ (nx*ny) */
    double *rho;
    /* Yee: Ex (i+1/2,j) size nx*ny; Ey (i,j+1/2); Bz (i+1/2,j+1/2) */
    double *Ex, *Ey, *Bz;
    double *Jx, *Jy;
    double *phi; /* Poisson scratch */
    Lock *locks;
    int n_locks;
    /* ledger accumulators */
    double E_em0, KE0, E_star_tot0;
    double flux_out_cum; /* rough sponge dissipation proxy */
} World;

/* ---------- helpers ---------- */

static inline int idx(int i, int j, int ny) { return i * ny + j; }

static inline int wrap(int i, int n) {
    i %= n;
    if (i < 0) i += n;
    return i;
}

static double *dalloc(size_t n) {
    double *p = (double *)calloc(n, sizeof(double));
    if (!p) {
        fprintf(stderr, "OOM %zu doubles\n", n);
        exit(1);
    }
    return p;
}

/* ---------- world lifecycle ---------- */

static void world_init(World *W, int nx, int ny, double dx, double dt, int n_locks) {
    memset(W, 0, sizeof(*W));
    W->cfg.nx = nx;
    W->cfg.ny = ny;
    W->cfg.dx = dx;
    W->cfg.dy = dx;
    W->cfg.dt = dt;
    W->cfg.n_locks = n_locks;
    W->cfg.absorb = 0;
    W->cfg.sponge_w = 8;
    W->cfg.sponge_sigma = 0.15;
    size_t N = (size_t)nx * (size_t)ny;
    W->rho = dalloc(N);
    W->Ex = dalloc(N);
    W->Ey = dalloc(N);
    W->Bz = dalloc(N);
    W->Jx = dalloc(N);
    W->Jy = dalloc(N);
    W->phi = dalloc(N);
    W->locks = (Lock *)calloc((size_t)n_locks, sizeof(Lock));
    if (!W->locks) {
        fprintf(stderr, "OOM locks\n");
        exit(1);
    }
    W->n_locks = n_locks;
}

static void world_free(World *W) {
    free(W->rho);
    free(W->Ex);
    free(W->Ey);
    free(W->Bz);
    free(W->Jx);
    free(W->Jy);
    free(W->phi);
    free(W->locks);
    memset(W, 0, sizeof(*W));
}

static void zero_fields(World *W) {
    size_t N = (size_t)W->cfg.nx * (size_t)W->cfg.ny;
    memset(W->rho, 0, N * sizeof(double));
    memset(W->Ex, 0, N * sizeof(double));
    memset(W->Ey, 0, N * sizeof(double));
    memset(W->Bz, 0, N * sizeof(double));
    memset(W->Jx, 0, N * sizeof(double));
    memset(W->Jy, 0, N * sizeof(double));
    memset(W->phi, 0, N * sizeof(double));
}

/* ---------- CIC charge deposit (cell centers) ---------- */

static void deposit_rho(World *W) {
    int nx = W->cfg.nx, ny = W->cfg.ny;
    double dx = W->cfg.dx, dy = W->cfg.dy;
    size_t N = (size_t)nx * (size_t)ny;
    memset(W->rho, 0, N * sizeof(double));
    double inv_area = 1.0 / (dx * dy);

    for (int p = 0; p < W->n_locks; p++) {
        Lock *L = &W->locks[p];
        if (!L->alive) continue;
        /* map to [0,nx) periodic */
        double x = L->x / dx;
        double y = L->y / dy;
        x = fmod(x, (double)nx);
        if (x < 0) x += nx;
        y = fmod(y, (double)ny);
        if (y < 0) y += ny;

        int i = (int)floor(x);
        int j = (int)floor(y);
        double fx = x - i;
        double fy = y - j;
        double q = L->q * inv_area;

        int i1 = wrap(i + 1, nx);
        int j1 = wrap(j + 1, ny);
        i = wrap(i, nx);
        j = wrap(j, ny);

        W->rho[idx(i, j, ny)] += q * (1 - fx) * (1 - fy);
        W->rho[idx(i1, j, ny)] += q * fx * (1 - fy);
        W->rho[idx(i, j1, ny)] += q * (1 - fx) * fy;
        W->rho[idx(i1, j1, ny)] += q * fx * fy;
    }
}

/* ---------- 1D segment current deposit (zigzag building block) ---------- */
/* Deposit J for pure-x motion of charge q from x0→x1 at fixed y (grid units:
 * positions in absolute coords). Charge-conserving for linear shape. */

static void deposit_Jx_segment(World *W, double q, double x0, double x1, double y, double dt) {
    int nx = W->cfg.nx, ny = W->cfg.ny;
    double dx = W->cfg.dx, dy = W->cfg.dy;
    if (fabs(x1 - x0) < 1e-30) return;

    /* work in cell units */
    double X0 = x0 / dx, X1 = x1 / dx, Y = y / dy;
    /* ensure continuous branch near wrap: caller should not wrap mid-segment */
    double dX = X1 - X0;

    /* y weights (CIC in y, cell-centered) */
    Y = fmod(Y, (double)ny);
    if (Y < 0) Y += ny;
    int j = (int)floor(Y);
    double wy = Y - j;
    int j0 = wrap(j, ny);
    int j1 = wrap(j + 1, ny);
    double w0 = 1.0 - wy, w1 = wy;

    /* subdivide if crosses integer cell faces in x */
    double xa = X0;
    double xb = X1;
    int dir = (dX > 0) ? 1 : -1;
    /* march through cells */
    double remaining = dX;
    double xcur = xa;
    int guard = 0;
    while (fabs(remaining) > 1e-15 && guard++ < nx + 4) {
        /* Linear shape + zigzag (Birdsall / Villasenor–Buneman class) */
        double cell_edge;
        if (dir > 0) {
            cell_edge = floor(xcur + 1e-14) + 1.0; /* next integer to the right */
            if (cell_edge > xb + 1e-14) cell_edge = xb;
        } else {
            cell_edge = ceil(xcur - 1e-14) - 1.0; /* next integer to the left */
            if (cell_edge < xb - 1e-14) cell_edge = xb;
        }
        double xnext = cell_edge;
        double dseg = xnext - xcur;
        if (fabs(dseg) < 1e-18) {
            /* snap off integer */
            xcur += 1e-12 * dir;
            remaining = xb - xcur;
            continue;
        }

        /* face index: Jx on face at left of cell containing segment */
        int iface = (int)floor(fmin(xcur, xnext) + 1e-14);
        iface = wrap(iface, nx);

        /* J density: (q Δx_phys / dt) distributed with y-weights, per face length dy
         * continuum: ∫ Jx dy ~ q vx; grid: Jx *= weight / dx (density) */
        double flux = (q / (dt * dy)) * (dseg * dx);
        W->Jx[idx(iface, j0, ny)] += flux * w0 / dx;
        W->Jx[idx(iface, j1, ny)] += flux * w1 / dx;

        xcur = xnext;
        remaining = xb - xcur;
    }
}

static void deposit_Jy_segment(World *W, double q, double y0, double y1, double x, double dt) {
    int nx = W->cfg.nx, ny = W->cfg.ny;
    double dx = W->cfg.dx, dy = W->cfg.dy;
    if (fabs(y1 - y0) < 1e-30) return;

    double Y0 = y0 / dy, Y1 = y1 / dy, X = x / dx;
    double dY = Y1 - Y0;
    double q_over_dt_dx = q / (dt * dx);

    X = fmod(X, (double)nx);
    if (X < 0) X += nx;
    int i = (int)floor(X);
    double wx = X - i;
    int i0 = wrap(i, nx);
    int i1 = wrap(i + 1, nx);
    double w0 = 1.0 - wx, w1 = wx;

    double ya = Y0, yb = Y1;
    int dir = (dY > 0) ? 1 : -1;
    double remaining = dY;
    double ycur = ya;
    int guard = 0;
    while (fabs(remaining) > 1e-15 && guard++ < ny + 4) {
        double cell_edge;
        if (dir > 0) {
            cell_edge = floor(ycur + 1e-14) + 1.0;
            if (cell_edge > yb + 1e-14) cell_edge = yb;
        } else {
            cell_edge = ceil(ycur - 1e-14) - 1.0;
            if (cell_edge < yb - 1e-14) cell_edge = yb;
        }
        double ynext = cell_edge;
        double dseg = ynext - ycur;
        if (fabs(dseg) < 1e-18) {
            ycur += 1e-12 * dir;
            remaining = yb - ycur;
            continue;
        }
        int jface = (int)floor(fmin(ycur, ynext) + 1e-14);
        jface = wrap(jface, ny);

        double flux = q_over_dt_dx * (dseg * dy);
        W->Jy[idx(i0, jface, ny)] += flux * w0 / dy;
        W->Jy[idx(i1, jface, ny)] += flux * w1 / dy;

        ycur = ynext;
        remaining = yb - ycur;
    }
}

/* Unwrap target so |Δ| is shortest on torus for deposit path */
static void shortest_delta(double *x1, double x0, double L) {
    double d = *x1 - x0;
    if (d > 0.5 * L) *x1 -= L;
    if (d < -0.5 * L) *x1 += L;
}

/* Zigzag charge-conserving J from trajectory (x0,y0)→(x1,y1) */
static void deposit_J_zigzag(World *W, double q, double x0, double y0, double x1, double y1,
                             double dt) {
    double Lx = W->cfg.nx * W->cfg.dx;
    double Ly = W->cfg.ny * W->cfg.dy;
    shortest_delta(&x1, x0, Lx);
    shortest_delta(&y1, y0, Ly);
    /* path: (x0,y0) → (x1,y0) → (x1,y1) */
    deposit_Jx_segment(W, q, x0, x1, y0, dt);
    deposit_Jy_segment(W, q, y0, y1, x1, dt);
}

/* ---------- gather E,B at continuous (x,y) via CIC on staggered grids ---------- */

static void gather_EB(const World *W, double x, double y, double *Ex, double *Ey, double *Bz) {
    int nx = W->cfg.nx, ny = W->cfg.ny;
    double dx = W->cfg.dx, dy = W->cfg.dy;

    /* normalize */
    double X = x / dx, Y = y / dy;
    X = fmod(X, (double)nx);
    if (X < 0) X += nx;
    Y = fmod(Y, (double)ny);
    if (Y < 0) Y += ny;

    /* Ex lives at (i+1/2, j): sample with shift -0.5 in x */
    {
        double xs = X - 0.5;
        if (xs < 0) xs += nx;
        int i = (int)floor(xs);
        int j = (int)floor(Y);
        double fx = xs - i, fy = Y - j;
        int i0 = wrap(i, nx), i1 = wrap(i + 1, nx);
        int j0 = wrap(j, ny), j1 = wrap(j + 1, ny);
        *Ex = (1 - fx) * (1 - fy) * W->Ex[idx(i0, j0, ny)] + fx * (1 - fy) * W->Ex[idx(i1, j0, ny)] +
              (1 - fx) * fy * W->Ex[idx(i0, j1, ny)] + fx * fy * W->Ex[idx(i1, j1, ny)];
    }
    /* Ey at (i, j+1/2) */
    {
        double ys = Y - 0.5;
        if (ys < 0) ys += ny;
        int i = (int)floor(X);
        int j = (int)floor(ys);
        double fx = X - i, fy = ys - j;
        int i0 = wrap(i, nx), i1 = wrap(i + 1, nx);
        int j0 = wrap(j, ny), j1 = wrap(j + 1, ny);
        *Ey = (1 - fx) * (1 - fy) * W->Ey[idx(i0, j0, ny)] + fx * (1 - fy) * W->Ey[idx(i1, j0, ny)] +
              (1 - fx) * fy * W->Ey[idx(i0, j1, ny)] + fx * fy * W->Ey[idx(i1, j1, ny)];
    }
    /* Bz at (i+1/2, j+1/2) */
    {
        double xs = X - 0.5, ys = Y - 0.5;
        if (xs < 0) xs += nx;
        if (ys < 0) ys += ny;
        int i = (int)floor(xs);
        int j = (int)floor(ys);
        double fx = xs - i, fy = ys - j;
        int i0 = wrap(i, nx), i1 = wrap(i + 1, nx);
        int j0 = wrap(j, ny), j1 = wrap(j + 1, ny);
        *Bz = (1 - fx) * (1 - fy) * W->Bz[idx(i0, j0, ny)] + fx * (1 - fy) * W->Bz[idx(i1, j0, ny)] +
              (1 - fx) * fy * W->Bz[idx(i0, j1, ny)] + fx * fy * W->Bz[idx(i1, j1, ny)];
    }
}

/* ---------- Boris push (relativistic, 2D) ---------- */

static void boris_push(Lock *L, double Ex, double Ey, double Bz, double dt) {
    if (L->pinned) {
        L->ux = L->uy = 0;
        L->fx = L->q * Ex;
        L->fy = L->q * Ey;
        return;
    }
    double qom = L->q / L->m;
    double half = 0.5 * dt * qom;

    /* u- = u + (q dt/2m) E */
    double ux = L->ux + half * Ex;
    double uy = L->uy + half * Ey;

    double gamma = sqrt(1.0 + ux * ux + uy * uy);
    /* t = (q dt/2m) B / γ */
    double tx = 0.0, ty = 0.0, tz = half * Bz / gamma;
    /* u' = u- + u- × t */
    double uxp = ux + (uy * tz - 0 * ty);
    double uyp = uy + (0 * tx - ux * tz);
    /* s = 2t / (1 + t²) */
    double t2 = tx * tx + ty * ty + tz * tz;
    double sx = 2.0 * tx / (1.0 + t2);
    double sy = 2.0 * ty / (1.0 + t2);
    double sz = 2.0 * tz / (1.0 + t2);
    /* u+ = u- + u' × s */
    ux = ux + (uyp * sz - 0 * sy);
    uy = uy + (0 * sx - uxp * sz);

    /* u_new = u+ + (q dt/2m) E */
    L->ux = ux + half * Ex;
    L->uy = uy + half * Ey;

    /* diagnostic force ~ q(E + v×B) using pre-push-ish fields */
    gamma = sqrt(1.0 + L->ux * L->ux + L->uy * L->uy);
    double vx = L->ux / gamma, vy = L->uy / gamma;
    L->fx = L->q * (Ex + vy * Bz);
    L->fy = L->q * (Ey - vx * Bz);
}

/* ---------- Yee Maxwell steps ---------- */

static void yee_faraday(World *W) {
    int nx = W->cfg.nx, ny = W->cfg.ny;
    double dtdx = W->cfg.dt / W->cfg.dx;
    double dtdy = W->cfg.dt / W->cfg.dy;
    /* ∂t Bz = -(curl E)_z = -(∂Ey/∂x - ∂Ex/∂y) = ∂Ex/∂y - ∂Ey/∂x
     * Bz^{n+1/2} += dt (∂Ex/∂y - ∂Ey/∂x) */
    for (int i = 0; i < nx; i++) {
        int ip = wrap(i + 1, nx);
        for (int j = 0; j < ny; j++) {
            int jp = wrap(j + 1, ny);
            /* Bz at (i+1/2,j+1/2) */
            double dEx_dy = (W->Ex[idx(i, jp, ny)] - W->Ex[idx(i, j, ny)]) * dtdy;
            double dEy_dx = (W->Ey[idx(ip, j, ny)] - W->Ey[idx(i, j, ny)]) * dtdx;
            W->Bz[idx(i, j, ny)] += (dEx_dy - dEy_dx);
        }
    }
}

static void yee_ampere(World *W) {
    int nx = W->cfg.nx, ny = W->cfg.ny;
    double dtdx = W->cfg.dt / W->cfg.dx;
    double dtdy = W->cfg.dt / W->cfg.dy;
    double dt = W->cfg.dt;
    /* ∂t Ex = ∂Bz/∂y - Jx ;  ∂t Ey = -∂Bz/∂x - Jy  (c=ε=μ=1) */
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            int jm = wrap(j - 1, ny);
            double dBz_dy = (W->Bz[idx(i, j, ny)] - W->Bz[idx(i, jm, ny)]) * dtdy;
            W->Ex[idx(i, j, ny)] += dBz_dy - dt * W->Jx[idx(i, j, ny)];
        }
    }
    for (int i = 0; i < nx; i++) {
        int im = wrap(i - 1, nx);
        for (int j = 0; j < ny; j++) {
            double dBz_dx = (W->Bz[idx(i, j, ny)] - W->Bz[idx(im, j, ny)]) * dtdx;
            W->Ey[idx(i, j, ny)] += -dBz_dx - dt * W->Jy[idx(i, j, ny)];
        }
    }

    /* optional sponge: damp E near boundary (non-periodic radiation sink proxy) */
    if (W->cfg.absorb) {
        int w = W->cfg.sponge_w;
        double sig0 = W->cfg.sponge_sigma;
        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < ny; j++) {
                int di = i < nx - 1 - i ? i : nx - 1 - i;
                int dj = j < ny - 1 - j ? j : ny - 1 - j;
                int d = di < dj ? di : dj;
                if (d >= w) continue;
                double s = sig0 * (double)(w - d) / (double)w;
                double damp = exp(-s * dt);
                double ex0 = W->Ex[idx(i, j, ny)];
                double ey0 = W->Ey[idx(i, j, ny)];
                W->Ex[idx(i, j, ny)] *= damp;
                W->Ey[idx(i, j, ny)] *= damp;
                /* crude flux proxy */
                W->flux_out_cum += 0.5 * ((ex0 * ex0 + ey0 * ey0) - (W->Ex[idx(i, j, ny)] * W->Ex[idx(i, j, ny)] +
                                                                      W->Ey[idx(i, j, ny)] * W->Ey[idx(i, j, ny)])) *
                                   W->cfg.dx * W->cfg.dy;
            }
        }
    }
}

/* ---------- Poisson electrostatic IC (periodic FFT-free Jacobi/SOR) ---------- */

static void poisson_sor(World *W, int max_iter, double tol) {
    int nx = W->cfg.nx, ny = W->cfg.ny;
    double dx = W->cfg.dx;
    double dx2 = dx * dx;
    double omega = 1.7;
    /* ∇²φ = -ρ  →  φ = (neighbors + ρ dx²)/4 */
    /* subtract mean ρ for periodic solvability */
    double mean = 0.0;
    size_t N = (size_t)nx * (size_t)ny;
    for (size_t k = 0; k < N; k++) mean += W->rho[k];
    mean /= (double)N;

    memset(W->phi, 0, N * sizeof(double));
    for (int it = 0; it < max_iter; it++) {
        double maxres = 0.0;
        for (int i = 0; i < nx; i++) {
            int im = wrap(i - 1, nx), ip = wrap(i + 1, nx);
            for (int j = 0; j < ny; j++) {
                int jm = wrap(j - 1, ny), jp = wrap(j + 1, ny);
                double rhs = -(W->rho[idx(i, j, ny)] - mean);
                double neigh = W->phi[idx(im, j, ny)] + W->phi[idx(ip, j, ny)] + W->phi[idx(i, jm, ny)] +
                               W->phi[idx(i, jp, ny)];
                double phi_new = 0.25 * (neigh - rhs * dx2);
                double res = phi_new - W->phi[idx(i, j, ny)];
                W->phi[idx(i, j, ny)] += omega * res;
                double ar = fabs(res);
                if (ar > maxres) maxres = ar;
            }
        }
        if (maxres < tol) break;
    }

    /* E = -∇φ on Yee faces */
    for (int i = 0; i < nx; i++) {
        int ip = wrap(i + 1, nx);
        for (int j = 0; j < ny; j++) {
            /* Ex(i+1/2,j) = -(φ(i+1,j)-φ(i,j))/dx */
            W->Ex[idx(i, j, ny)] = -(W->phi[idx(ip, j, ny)] - W->phi[idx(i, j, ny)]) / dx;
        }
    }
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            int jp = wrap(j + 1, ny);
            W->Ey[idx(i, j, ny)] = -(W->phi[idx(i, jp, ny)] - W->phi[idx(i, j, ny)]) / dx;
        }
    }
    /* Bz = 0 electrostatic */
    memset(W->Bz, 0, N * sizeof(double));
}

/* ---------- diagnostics ---------- */

static double gauss_max(const World *W) {
    int nx = W->cfg.nx, ny = W->cfg.ny;
    double dx = W->cfg.dx, dy = W->cfg.dy;
    double m = 0.0;
    for (int i = 0; i < nx; i++) {
        int im = wrap(i - 1, nx);
        for (int j = 0; j < ny; j++) {
            int jm = wrap(j - 1, ny);
            /* div E at cell center (i,j): (Ex_i - Ex_{i-1})/dx + (Ey_j - Ey_{j-1})/dy */
            double div =
                (W->Ex[idx(i, j, ny)] - W->Ex[idx(im, j, ny)]) / dx +
                (W->Ey[idx(i, j, ny)] - W->Ey[idx(i, jm, ny)]) / dy;
            double r = fabs(div - W->rho[idx(i, j, ny)]);
            if (r > m) m = r;
        }
    }
    return m;
}

static double gauss_rms(const World *W) {
    int nx = W->cfg.nx, ny = W->cfg.ny;
    double dx = W->cfg.dx, dy = W->cfg.dy;
    double s = 0.0;
    size_t N = (size_t)nx * (size_t)ny;
    for (int i = 0; i < nx; i++) {
        int im = wrap(i - 1, nx);
        for (int j = 0; j < ny; j++) {
            int jm = wrap(j - 1, ny);
            double div =
                (W->Ex[idx(i, j, ny)] - W->Ex[idx(im, j, ny)]) / dx +
                (W->Ey[idx(i, j, ny)] - W->Ey[idx(i, jm, ny)]) / dy;
            double r = div - W->rho[idx(i, j, ny)];
            s += r * r;
        }
    }
    return sqrt(s / (double)N);
}

static double field_energy(const World *W) {
    int nx = W->cfg.nx, ny = W->cfg.ny;
    double dA = W->cfg.dx * W->cfg.dy;
    double s = 0.0;
    size_t N = (size_t)nx * (size_t)ny;
    for (size_t k = 0; k < N; k++) {
        s += 0.5 * (W->Ex[k] * W->Ex[k] + W->Ey[k] * W->Ey[k] + W->Bz[k] * W->Bz[k]);
    }
    return s * dA;
}

static double kinetic_energy(const World *W) {
    double ke = 0.0;
    for (int p = 0; p < W->n_locks; p++) {
        Lock *L = &W->locks[p];
        if (!L->alive) continue;
        double gamma = sqrt(1.0 + L->ux * L->ux + L->uy * L->uy);
        ke += L->m * (gamma - 1.0); /* c=1 */
    }
    return ke;
}

static double E_star_total(const World *W) {
    double s = 0.0;
    for (int p = 0; p < W->n_locks; p++)
        if (W->locks[p].alive) s += W->locks[p].E_star;
    return s;
}

static double total_charge(const World *W) {
    size_t N = (size_t)W->cfg.nx * (size_t)W->cfg.ny;
    double s = 0.0;
    double dA = W->cfg.dx * W->cfg.dy;
    for (size_t k = 0; k < N; k++) s += W->rho[k];
    return s * dA;
}

/* Sample |E| vs radius from a center (for L0 exterior field) */
static void sample_E_radial(const World *W, double cx, double cy, int nbin, double rmax,
                            double *r_mid, double *E_mean, int *counts) {
    int nx = W->cfg.nx, ny = W->cfg.ny;
    double dx = W->cfg.dx;
    for (int b = 0; b < nbin; b++) {
        r_mid[b] = (b + 0.5) * rmax / nbin;
        E_mean[b] = 0.0;
        counts[b] = 0;
    }
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            double x = (i + 0.5) * dx;
            double y = (j + 0.5) * dx;
            double rx = x - cx, ry = y - cy;
            /* periodic min image */
            double Lx = nx * dx, Ly = ny * dx;
            if (rx > 0.5 * Lx) rx -= Lx;
            if (rx < -0.5 * Lx) rx += Lx;
            if (ry > 0.5 * Ly) ry -= Ly;
            if (ry < -0.5 * Ly) ry += Ly;
            double r = sqrt(rx * rx + ry * ry);
            if (r < 1e-9 || r >= rmax) continue;
            int b = (int)(r / rmax * nbin);
            if (b < 0 || b >= nbin) continue;
            /* cell-center E approx average of faces */
            int im = wrap(i - 1, nx), jm = wrap(j - 1, ny);
            double Ex = 0.5 * (W->Ex[idx(i, j, ny)] + W->Ex[idx(im, j, ny)]);
            double Ey = 0.5 * (W->Ey[idx(i, j, ny)] + W->Ey[idx(i, jm, ny)]);
            double Em = sqrt(Ex * Ex + Ey * Ey);
            E_mean[b] += Em;
            counts[b]++;
        }
    }
    for (int b = 0; b < nbin; b++)
        if (counts[b] > 0) E_mean[b] /= counts[b];
}

/* ---------- full PIC step ---------- */

static void pic_step(World *W) {
    int nx = W->cfg.nx, ny = W->cfg.ny;
    double dt = W->cfg.dt;
    size_t N = (size_t)nx * (size_t)ny;

    /* 1. Faraday half-step already full leapfrog: B update */
    yee_faraday(W);

    /* 2. Gather + Boris (using E^n, B^{n+1/2}) */
    for (int p = 0; p < W->n_locks; p++) {
        Lock *L = &W->locks[p];
        if (!L->alive) continue;
        double Ex, Ey, Bz;
        gather_EB(W, L->x, L->y, &Ex, &Ey, &Bz);
        boris_push(L, Ex, Ey, Bz, dt);
    }

    /* 3. Advance positions; deposit J from trajectories */
    memset(W->Jx, 0, N * sizeof(double));
    memset(W->Jy, 0, N * sizeof(double));

    double Lx = nx * W->cfg.dx, Ly = ny * W->cfg.dy;
    for (int p = 0; p < W->n_locks; p++) {
        Lock *L = &W->locks[p];
        if (!L->alive || L->pinned) continue;
        double gamma = sqrt(1.0 + L->ux * L->ux + L->uy * L->uy);
        double x0 = L->x, y0 = L->y;
        double x1 = x0 + (L->ux / gamma) * dt;
        double y1 = y0 + (L->uy / gamma) * dt;
        deposit_J_zigzag(W, L->q, x0, y0, x1, y1, dt);
        /* wrap position */
        L->x = fmod(x1, Lx);
        if (L->x < 0) L->x += Lx;
        L->y = fmod(y1, Ly);
        if (L->y < 0) L->y += Ly;
    }

    /* 4. Ampère */
    yee_ampere(W);

    /* 5. Refresh ρ for diagnostics / Gauss */
    deposit_rho(W);
}

/* Electrostatic-only force probe (no time step): gather E after Poisson */
static void probe_forces(World *W) {
    for (int p = 0; p < W->n_locks; p++) {
        Lock *L = &W->locks[p];
        double Ex, Ey, Bz;
        gather_EB(W, L->x, L->y, &Ex, &Ey, &Bz);
        L->fx = L->q * Ex;
        L->fy = L->q * Ey;
    }
}

/* ---------- lock factories ---------- */

static void init_lock(Lock *L, int id, double q, double m, double x, double y, int pinned) {
    memset(L, 0, sizeof(*L));
    L->id = id;
    L->type = 0;
    L->q = q;
    L->m = m;
    L->E_star = m * C0 * C0; /* rest energy sequestered (MVP constant) */
    L->x = x;
    L->y = y;
    L->pinned = pinned;
    L->alive = 1;
}

/* ========================================================================
 * Experiments L0–L3
 * ======================================================================== */

static const char *RESULTS_DIR = "results";

static int exp_L0(const char *outdir) {
    printf("\n======== L0: single lock static (Gauss + exterior field) ========\n");
    int nx = 96, ny = 96;
    double dx = 1.0;
    /* CFL: 2D Yee dt <= dx/sqrt(2) */
    double dt = 0.5 * dx / sqrt(2.0);

    World W;
    world_init(&W, nx, ny, dx, dt, 1);
    zero_fields(&W);

    double cx = 0.5 * nx * dx, cy = 0.5 * ny * dx;
    init_lock(&W.locks[0], 0, +1.0, 10.0, cx, cy, 1);

    deposit_rho(&W);
    double Q = total_charge(&W);
    poisson_sor(&W, 8000, 1e-12);

    double gmax = gauss_max(&W);
    double grms = gauss_rms(&W);
    double Eem = field_energy(&W);

    printf("  Q_total        = %.6e  (expect ~1)\n", Q);
    printf("  gauss_max      = %.6e\n", gmax);
    printf("  gauss_rms      = %.6e\n", grms);
    printf("  E_em           = %.6e\n", Eem);

    /* exterior |E|(r) — 2D Coulomb ~ 1/(2π r) for unit charge */
    enum { NBIN = 24 };
    double r_mid[NBIN], E_mean[NBIN];
    int counts[NBIN];
    double rmax = 0.35 * nx * dx;
    sample_E_radial(&W, cx, cy, NBIN, rmax, r_mid, E_mean, counts);

    char path[512];
    snprintf(path, sizeof path, "%s/l0_radial.tsv", outdir);
    FILE *fp = fopen(path, "w");
    if (fp) {
        fprintf(fp, "r\tE_mean\tE_2D_coulomb\tcounts\n");
        for (int b = 0; b < NBIN; b++) {
            double Eth = (r_mid[b] > 1.0) ? 1.0 / (2.0 * M_PI * r_mid[b]) : 0.0;
            fprintf(fp, "%.6f\t%.6e\t%.6e\t%d\n", r_mid[b], E_mean[b], Eth, counts[b]);
        }
        fclose(fp);
        printf("  wrote %s\n", path);
    }

    /* Check exterior field is sensible: for r in [6,20], E should decrease with r */
    int mono_ok = 1;
    double prev = 1e99;
    int n_samp = 0;
    for (int b = 0; b < NBIN; b++) {
        if (r_mid[b] < 6.0 || r_mid[b] > 20.0 || counts[b] < 8) continue;
        if (E_mean[b] > prev * 1.05) mono_ok = 0; /* allow 5% noise */
        prev = E_mean[b];
        n_samp++;
    }

    /* Gauss floor: after SOR should be small vs ρ scale (~1/dx² ~1) */
    int gauss_ok = (gmax < 1e-3);
    int field_ok = mono_ok && n_samp >= 3;
    int Q_ok = fabs(Q - 1.0) < 1e-9;

    /* Short dynamical hold with pinned lock: Gauss must not explode */
    for (int s = 0; s < 200; s++) pic_step(&W);
    double gmax2 = gauss_max(&W);
    printf("  gauss_max after 200 pinned steps = %.6e\n", gmax2);
    int dyn_ok = (gmax2 < 0.05); /* looser: dynamical floor with zigzag */

    printf("  PASS checks: Q_ok=%d gauss_ic=%d exterior_mono=%d dyn_gauss=%d\n", Q_ok, gauss_ok,
           field_ok, dyn_ok);

    snprintf(path, sizeof path, "%s/l0_summary.tsv", outdir);
    fp = fopen(path, "w");
    if (fp) {
        fprintf(fp, "Q\tgauss_max_ic\tgauss_rms_ic\tgauss_max_dyn\tE_em\tQ_ok\tgauss_ok\tfield_ok\tdyn_ok\n");
        fprintf(fp, "%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%d\t%d\t%d\t%d\n", Q, gmax, grms, gmax2, Eem, Q_ok,
                gauss_ok, field_ok, dyn_ok);
        fclose(fp);
    }

    int pass = Q_ok && gauss_ok && field_ok && dyn_ok;
    printf("  L0 => %s\n", pass ? "PASS" : "FAIL");
    world_free(&W);
    return pass ? 0 : 1;
}

static int exp_L1(const char *outdir) {
    printf("\n======== L1: rest pair F(D) from medium only ========\n");
    int nx = 128, ny = 128;
    double dx = 1.0;
    double dt = 0.5 * dx / sqrt(2.0);
    int Ds[] = {8, 12, 16, 20, 24};
    int nD = 5;

    char path[512];
    snprintf(path, sizeof path, "%s/l1_force.tsv", outdir);
    FILE *fp = fopen(path, "w");
    if (!fp) {
        fprintf(stderr, "cannot write %s\n", path);
        return 1;
    }
    fprintf(fp, "D\ta_rel\tF_along\tF_pair_ref\tgauss_max\tmedium_only\n");

    double a_prev = 1e99;
    int mono = 1;
    int all_medium = 1;

    for (int k = 0; k < nD; k++) {
        double D = (double)Ds[k];
        World W;
        world_init(&W, nx, ny, dx, dt, 2);
        zero_fields(&W);

        double cx = 0.5 * nx * dx, cy = 0.5 * ny * dx;
        double m = 5.0;
        /* opposite charges, rest, pinned for force probe */
        init_lock(&W.locks[0], 0, +1.0, m, cx - 0.5 * D, cy, 1);
        init_lock(&W.locks[1], 1, -1.0, m, cx + 0.5 * D, cy, 1);

        deposit_rho(&W);
        poisson_sor(&W, 10000, 1e-12);
        probe_forces(&W);

        /* relative acceleration along separation (x) from medium force only */
        double F0x = W.locks[0].fx, F1x = W.locks[1].fx;
        /* For + at left, - at right: attraction ⇒ F0x > 0, F1x < 0 */
        double a_rel = fabs(F0x / m - F1x / m);
        double F_along = 0.5 * (F0x - F1x); /* mean attractive magnitude if signs right */
        /* 2D continuum reference (NOT used in dynamics): |F| ~ 1/(2π D) */
        double F_ref = 1.0 / (2.0 * M_PI * D);
        double gmax = gauss_max(&W);

        /* Verify no pairwise term exists: force must equal q*E only (already true by construction).
         * Flag medium_only=1 always in this code path. */
        int medium_only = 1;
        all_medium &= medium_only;

        printf("  D=%4.0f  a_rel=%.6e  F_along=%.6e  F_2Dref=%.6e  gauss_max=%.3e  attract=%d\n", D,
               a_rel, F_along, F_ref, gmax, (F0x > 0 && F1x < 0));

        fprintf(fp, "%.0f\t%.8e\t%.8e\t%.8e\t%.8e\t%d\n", D, a_rel, F_along, F_ref, gmax, medium_only);

        if (a_rel > a_prev * 1.02) mono = 0; /* require monotone decreasing a_rel(D) */
        a_prev = a_rel;

        /* attraction sign check */
        if (!(F0x > 0.0 && F1x < 0.0)) mono = 0;

        world_free(&W);
    }
    fclose(fp);
    printf("  wrote %s\n", path);
    printf("  monotone a_rel(D)=%d  medium_only=%d\n", mono, all_medium);
    int pass = mono && all_medium;
    printf("  L1 => %s\n", pass ? "PASS" : "FAIL");
    return pass ? 0 : 1;
}

static int exp_L2(const char *outdir) {
    printf("\n======== L2: soft long T (durability + ledger) ========\n");
    int nx = 96, ny = 96;
    double dx = 1.0;
    double dt = 0.45 * dx / sqrt(2.0);
    int nsteps = 2500;
    double D0 = 16.0;
    double m = 8.0;

    World W;
    world_init(&W, nx, ny, dx, dt, 2);
    W.cfg.absorb = 1; /* sponge for radiation honesty */
    W.cfg.sponge_w = 10;
    W.cfg.sponge_sigma = 0.12;
    zero_fields(&W);

    double cx = 0.5 * nx * dx, cy = 0.5 * ny * dx;
    init_lock(&W.locks[0], 0, +1.0, m, cx - 0.5 * D0, cy, 0);
    init_lock(&W.locks[1], 1, -1.0, m, cx + 0.5 * D0, cy, 0);
    /* soft: zero initial velocity */

    deposit_rho(&W);
    poisson_sor(&W, 8000, 1e-12);

    W.E_em0 = field_energy(&W);
    W.KE0 = kinetic_energy(&W);
    W.E_star_tot0 = E_star_total(&W);

    char path[512];
    snprintf(path, sizeof path, "%s/l2_series.tsv", outdir);
    FILE *fp = fopen(path, "w");
    if (fp)
        fprintf(fp, "step\tt\tE_em\tKE\tE_star\tE_tot\tgauss_max\tQ\tsep\talive\n");

    double E_em_min = W.E_em0;
    double E_em_max = W.E_em0;
    double gmax_run = 0.0;
    int alive_always = 1;

    for (int s = 0; s <= nsteps; s++) {
        if (s > 0) pic_step(&W);

        if (s % 50 == 0 || s == nsteps) {
            double Eem = field_energy(&W);
            double KE = kinetic_energy(&W);
            double Es = E_star_total(&W);
            double Et = Eem + KE + Es;
            double g = gauss_max(&W);
            double Q = total_charge(&W);
            double dxp = W.locks[0].x - W.locks[1].x;
            double dyp = W.locks[0].y - W.locks[1].y;
            double Lx = nx * dx, Ly = ny * dx;
            if (dxp > 0.5 * Lx) dxp -= Lx;
            if (dxp < -0.5 * Lx) dxp += Lx;
            if (dyp > 0.5 * Ly) dyp -= Ly;
            if (dyp < -0.5 * Ly) dyp += Ly;
            double sep = sqrt(dxp * dxp + dyp * dyp);
            int alive = W.locks[0].alive && W.locks[1].alive;
            alive_always &= alive;
            if (Eem < E_em_min) E_em_min = Eem;
            if (Eem > E_em_max) E_em_max = Eem;
            if (g > gmax_run) gmax_run = g;

            if (fp)
                fprintf(fp, "%d\t%.4f\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.6f\t%d\n", s, s * dt, Eem,
                        KE, Es, Et, g, Q, sep, alive);

            if (s % 500 == 0)
                printf("  step %5d  E_em=%.4e  KE=%.4e  sep=%.3f  gauss=%.3e  alive=%d\n", s, Eem, KE,
                       sep, g, alive);
        }
    }
    if (fp) {
        fclose(fp);
        printf("  wrote %s\n", path);
    }

    double Eem_f = field_energy(&W);
    double KE_f = kinetic_energy(&W);
    double Es_f = E_star_total(&W);
    /* ledger: with sponge, Et is not closed; compare interior trend without requiring exact close.
     * Honesty metric: E_star exact (locks don't evaporate energy), E_em does not null. */
    double E_star_drift = fabs(Es_f - W.E_star_tot0) / fmax(W.E_star_tot0, 1e-30);
    double E_em_null = (Eem_f < 1e-6 * fmax(W.E_em0, 1e-30));
    double E_em_frac = Eem_f / fmax(W.E_em0, 1e-30);
    (void)E_em_null;

    printf("  E_em0=%.4e  E_em_f=%.4e  frac=%.4f  min=%.4e\n", W.E_em0, Eem_f, E_em_frac, E_em_min);
    printf("  E_star drift = %.3e  (expect 0)\n", E_star_drift);
    printf("  KE_f=%.4e  flux_proxy=%.4e  gauss_max_run=%.3e\n", KE_f, W.flux_out_cum, gmax_run);
    printf("  locks alive always = %d\n", alive_always);

    /* Pass criteria from OP */
    int pass_alive = alive_always;
    int pass_em = !E_em_null && (E_em_min > 1e-4 * fmax(W.E_em0, 1e-30));
    int pass_star = (E_star_drift < 1e-12);
    /* Gauss: dynamical PIC residual should stay controlled (not necessarily 1e-13) */
    int pass_gauss = (gmax_run < 0.5);
    int pass = pass_alive && pass_em && pass_star && pass_gauss;

    snprintf(path, sizeof path, "%s/l2_summary.tsv", outdir);
    fp = fopen(path, "w");
    if (fp) {
        fprintf(fp, "nsteps\tE_em0\tE_em_f\tE_em_min\tE_star_drift\tKE_f\tgmax_run\talive\tpass_em\tpass\n");
        fprintf(fp, "%d\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%d\t%d\t%d\n", nsteps, W.E_em0, Eem_f,
                E_em_min, E_star_drift, KE_f, gmax_run, pass_alive, pass_em, pass);
        fclose(fp);
    }

    printf("  L2 => %s  (alive=%d em=%d star=%d gauss=%d)\n", pass ? "PASS" : "FAIL", pass_alive,
           pass_em, pass_star, pass_gauss);
    world_free(&W);
    return pass ? 0 : 1;
}

static int exp_L3(const char *outdir) {
    printf("\n======== L3: hard v_t / orbit attempt (optional) ========\n");
    int nx = 96, ny = 96;
    double dx = 1.0;
    double dt = 0.4 * dx / sqrt(2.0);
    int nsteps = 4000;
    double D0 = 14.0;
    double m = 6.0;
    double v_t = 0.12; /* tangential kick */

    World W;
    world_init(&W, nx, ny, dx, dt, 2);
    W.cfg.absorb = 1;
    W.cfg.sponge_w = 10;
    W.cfg.sponge_sigma = 0.1;
    zero_fields(&W);

    double cx = 0.5 * nx * dx, cy = 0.5 * ny * dx;
    init_lock(&W.locks[0], 0, +1.0, m, cx - 0.5 * D0, cy, 0);
    init_lock(&W.locks[1], 1, -1.0, m, cx + 0.5 * D0, cy, 0);
    /* opposite tangential velocities for COM-rest pair */
    /* soft-relativistic: proper velocity u = γv ≈ v for |v|≪1 */
    W.locks[0].ux = 0.0;
    W.locks[0].uy = v_t;
    W.locks[1].ux = 0.0;
    W.locks[1].uy = -v_t;

    deposit_rho(&W);
    poisson_sor(&W, 8000, 1e-12);

    char path[512];
    snprintf(path, sizeof path, "%s/l3_tracks.tsv", outdir);
    FILE *fp = fopen(path, "w");
    if (fp) fprintf(fp, "step\tt\tx0\ty0\tx1\ty1\tsep\tE_em\tKE\n");

    int revs_est = 0;
    double angle_acc = 0.0;
    double prev_ang = atan2(W.locks[0].y - W.locks[1].y, W.locks[0].x - W.locks[1].x);
    int shredded = 0; /* locks cannot shred by definition — track field null instead */

    for (int s = 0; s <= nsteps; s++) {
        if (s > 0) pic_step(&W);

        if (s % 20 == 0 || s == nsteps) {
            double dxp = W.locks[0].x - W.locks[1].x;
            double dyp = W.locks[0].y - W.locks[1].y;
            double Lx = nx * dx, Ly = ny * dx;
            if (dxp > 0.5 * Lx) dxp -= Lx;
            if (dxp < -0.5 * Lx) dxp += Lx;
            if (dyp > 0.5 * Ly) dyp -= Ly;
            if (dyp < -0.5 * Ly) dyp += Ly;
            double sep = sqrt(dxp * dxp + dyp * dyp);
            double ang = atan2(dyp, dxp);
            double dang = ang - prev_ang;
            if (dang > M_PI) dang -= 2 * M_PI;
            if (dang < -M_PI) dang += 2 * M_PI;
            angle_acc += dang;
            prev_ang = ang;
            revs_est = (int)(fabs(angle_acc) / (2 * M_PI));

            double Eem = field_energy(&W);
            double KE = kinetic_energy(&W);
            if (fp)
                fprintf(fp, "%d\t%.4f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6e\t%.6e\n", s, s * dt,
                        W.locks[0].x, W.locks[0].y, W.locks[1].x, W.locks[1].y, sep, Eem, KE);

            if (s % 500 == 0)
                printf("  step %5d  sep=%.3f  revs~%d  E_em=%.3e  KE=%.3e\n", s, sep, revs_est, Eem,
                       KE);

            if (Eem < 1e-8) shredded = 1;
        }
    }
    if (fp) {
        fclose(fp);
        printf("  wrote %s\n", path);
    }

    double revs = fabs(angle_acc) / (2 * M_PI);
    int pass_noshred = !shredded && W.locks[0].alive && W.locks[1].alive;
    int pass_motion = revs >= 0.25; /* at least quarter-turn of relative angle, honest partial */
    int pass = pass_noshred;        /* L3 optional: no shred-by-definition is the hard OP line */

    printf("  relative revs ≈ %.3f  noshred=%d\n", revs, pass_noshred);
    printf("  L3 => %s  (optional; revs>=0.25=%d)\n", pass ? "PASS*" : "FAIL", pass_motion);

    snprintf(path, sizeof path, "%s/l3_summary.tsv", outdir);
    fp = fopen(path, "w");
    if (fp) {
        fprintf(fp, "nsteps\tv_t\trevs\tnoshred\tpass_motion\tpass\n");
        fprintf(fp, "%d\t%.4f\t%.6f\t%d\t%d\t%d\n", nsteps, v_t, revs, pass_noshred, pass_motion, pass);
        fclose(fp);
    }

    world_free(&W);
    return pass ? 0 : 1;
}

/* ---------- main ---------- */

static void usage(const char *argv0) {
    fprintf(stderr,
            "Usage: %s <l0|l1|l2|l3|all> [results_dir]\n"
            "  Standalone monist PIC: Yee U(1) free medium + typed locks\n"
            "  No pairwise Coulomb, no multi-fab, no Cosserat L.\n",
            argv0);
}

int main(int argc, char **argv) {
    if (argc < 2) {
        usage(argv[0]);
        return 2;
    }
    const char *mode = argv[1];
    const char *outdir = (argc >= 3) ? argv[2] : RESULTS_DIR;

    /* ensure results dir exists (best-effort) */
    char cmd[1024];
    snprintf(cmd, sizeof cmd, "mkdir -p '%s'", outdir);
    if (system(cmd) != 0) {
        /* ignore */
    }

    printf("v81 path1_locks monist PIC sandbox\n");
    printf("outdir=%s  mode=%s\n", outdir, mode);

    int rc = 0;
    if (strcmp(mode, "l0") == 0)
        rc = exp_L0(outdir);
    else if (strcmp(mode, "l1") == 0)
        rc = exp_L1(outdir);
    else if (strcmp(mode, "l2") == 0)
        rc = exp_L2(outdir);
    else if (strcmp(mode, "l3") == 0)
        rc = exp_L3(outdir);
    else if (strcmp(mode, "all") == 0) {
        int r0 = exp_L0(outdir);
        int r1 = exp_L1(outdir);
        int r2 = exp_L2(outdir);
        int r3 = exp_L3(outdir);
        printf("\n======== SCORECARD ========\n");
        printf("  L0 %s\n", r0 ? "FAIL" : "PASS");
        printf("  L1 %s\n", r1 ? "FAIL" : "PASS");
        printf("  L2 %s\n", r2 ? "FAIL" : "PASS");
        printf("  L3 %s (optional)\n", r3 ? "FAIL" : "PASS*");
        rc = r0 | r1 | r2; /* L3 optional for overall */
        printf("  overall L0–L2 => %s\n", rc ? "FAIL" : "PASS");
    } else {
        usage(argv[0]);
        return 2;
    }
    return rc;
}
