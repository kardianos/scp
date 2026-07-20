/* scp_locks.h — Stage-3 typed lock carriers on free U(1) gauge medium (v81 P1)
 *
 * REAL = FreeSubstrate (existing gauge A,E) + Locks (structs, not multiplet φ).
 * When n_locks==0: no allocation, no deposit, no push (nuclear path unchanged).
 *
 * Include once from scp_sim.c after Grid is defined. Provides:
 *   ScpLock, locks_alloc/free, locks_load_file, locks_deposit_rho,
 *   locks_add_J_to_Eacc, locks_push_and_move, locks_energy, locks_write_track
 *
 * Force:   F = −(g q) E + … (sign matched to projected lattice E; opposite attract)
 * Soft core (optional): short-range form-factor repulsion r < lock_soft_r
 * Gauss:   div E − g (ρ_matter + ρ_lock) = 0
 * Ampère:  E −= dt g j_phys after push (Cosserat J_lat = −j)
 */
#ifndef SCP_LOCKS_H
#define SCP_LOCKS_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* ScpLock is defined in scp_sim.c before Grid; this header is included after Grid. */

static inline int locks_wrap_i(int i, int N) {
    i %= N;
    if (i < 0) i += N;
    return i;
}

static inline long locks_idx(int i, int j, int k, int N) {
    return (long)i * N * N + (long)j * N + k;
}

/* Physical → continuous cell index in [0,N) for site-centered fields.
 * Site i lives at x_i = -L + i*dx, dx=2L/(N-1). Map x → (x+L)/dx. */
static inline void locks_phys_to_cell(const Grid *g, double x, double y, double z,
                                     double *ci, double *cj, double *ck) {
    double Nf = (double)g->N;
    double inv = 1.0 / g->dx;
    *ci = (x + g->L) * inv;
    *cj = (y + g->L) * inv;
    *ck = (z + g->L) * inv;
    /* fold into [0,N) using period ~ N (periodic indices) */
    *ci = fmod(*ci, Nf); if (*ci < 0) *ci += Nf;
    *cj = fmod(*cj, Nf); if (*cj < 0) *cj += Nf;
    *ck = fmod(*ck, Nf); if (*ck < 0) *ck += Nf;
}

static void locks_alloc(Grid *g, int n) {
    g->n_locks = n;
    g->locks = NULL;
    g->lock_rho = NULL;
    g->lock_J[0] = g->lock_J[1] = g->lock_J[2] = NULL;
    g->lock_bag = NULL;
    g->locks_track[0] = '\0';
    if (n <= 0) return;
    g->locks = (ScpLock *)calloc((size_t)n, sizeof(ScpLock));
    g->lock_rho = (double *)calloc((size_t)g->N3, sizeof(double));
    for (int d = 0; d < 3; d++)
        g->lock_J[d] = (double *)calloc((size_t)g->N3, sizeof(double));
    if (g->lock_bag_mode == 2)
        g->lock_bag = (double *)calloc((size_t)g->N3, sizeof(double));
    if (!g->locks || !g->lock_rho || !g->lock_J[0] || !g->lock_J[1] || !g->lock_J[2] ||
        (g->lock_bag_mode == 2 && !g->lock_bag)) {
        fprintf(stderr, "FATAL: locks alloc failed (n=%d)\n", n);
        exit(1);
    }
    printf("Locks: allocated %d carriers + rho/J%s\n", n,
           g->lock_bag_mode == 2 ? " + grid bag B(x)" : "");
}

static void locks_free(Grid *g) {
    if (!g) return;
    free(g->locks);
    free(g->lock_rho);
    free(g->lock_J[0]); free(g->lock_J[1]); free(g->lock_J[2]);
    free(g->lock_bag);
    g->locks = NULL;
    g->lock_rho = NULL;
    g->lock_J[0] = g->lock_J[1] = g->lock_J[2] = NULL;
    g->lock_bag = NULL;
    g->n_locks = 0;
}

/* TSV: id type q m E_star x y z ux uy uz pinned  (comments # allowed) */
static int locks_load_file(Grid *g, const char *path) {
    if (!path || !path[0] || g->n_locks <= 0) return 0;
    FILE *fp = fopen(path, "r");
    if (!fp) {
        fprintf(stderr, "FATAL: cannot open locks_file '%s'\n", path);
        exit(1);
    }
    char line[1024];
    int loaded = 0;
    while (fgets(line, sizeof(line), fp) && loaded < g->n_locks) {
        char *s = line;
        while (*s == ' ' || *s == '\t') s++;
        if (*s == '#' || *s == '\n' || *s == '\0') continue;
        ScpLock *L = &g->locks[loaded];
        int id = loaded, type = 0, pinned = 0;
        double q = 1, m = 1, Es = 1, x = 0, y = 0, z = 0, ux = 0, uy = 0, uz = 0;
        int n = sscanf(s, "%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %d",
                       &id, &type, &q, &m, &Es, &x, &y, &z, &ux, &uy, &uz, &pinned);
        if (n < 8) {
            /* compact: q m x y z [ux uy uz pinned [type]] */
            n = sscanf(s, "%lf %lf %lf %lf %lf %lf %lf %lf %d %d",
                       &q, &m, &x, &y, &z, &ux, &uy, &uz, &pinned, &type);
            if (n < 5) continue;
            id = loaded; Es = m;
        }
        if (type == 1) q = 0.0;
        L->id = id; L->type = type; L->q = q; L->m = m;
        L->E_star = (Es > 0) ? Es : m;
        L->x[0] = x; L->x[1] = y; L->x[2] = z;
        L->u[0] = ux; L->u[1] = uy; L->u[2] = uz;
        L->x_prev[0] = x; L->x_prev[1] = y; L->x_prev[2] = z;
        L->pinned = pinned; L->alive = 1;
        loaded++;
    }
    fclose(fp);
    if (loaded != g->n_locks) {
        fprintf(stderr, "FATAL: locks_file has %d rows, n_locks=%d\n", loaded, g->n_locks);
        exit(1);
    }
    printf("Locks: loaded %d from %s\n", loaded, path);
    for (int i = 0; i < loaded; i++) {
        ScpLock *L = &g->locks[i];
        printf("  lock[%d] type=%d q=%+.3g m=%.3g E*=%.3g x=(%.3f,%.3f,%.3f) pinned=%d\n",
               L->id, L->type, L->q, L->m, L->E_star, L->x[0], L->x[1], L->x[2], L->pinned);
    }
    return loaded;
}

/* Inline: "q,m,x,y,z[,ux,uy,uz,pinned[,type]]"
 * type: 0=charge lock (EM), 1=anti-lock / lock-Higgs bag (no EM ρ; bag pull on charges) */
static void locks_set_from_csv(Grid *g, int i, const char *csv) {
    if (i < 0 || i >= g->n_locks || !csv) return;
    ScpLock *L = &g->locks[i];
    double q = 1, m = 1, x = 0, y = 0, z = 0, ux = 0, uy = 0, uz = 0;
    int pinned = 0, type = 0;
    int n = sscanf(csv, "%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%d,%d",
                   &q, &m, &x, &y, &z, &ux, &uy, &uz, &pinned, &type);
    if (n < 5) return;
    if (type == 1) q = 0.0; /* anti-lock: neutral under U(1) */
    L->id = i; L->type = type; L->q = q; L->m = m; L->E_star = m;
    L->x[0] = x; L->x[1] = y; L->x[2] = z;
    L->u[0] = ux; L->u[1] = uy; L->u[2] = uz;
    L->x_prev[0] = x; L->x_prev[1] = y; L->x_prev[2] = z;
    L->pinned = pinned; L->alive = 1;
}

/* CIC deposit of charge density: ∫ ρ dV = Σ q  (charge locks only; anti-locks skip) */
static void locks_deposit_rho(Grid *g) {
    if (g->n_locks <= 0 || !g->lock_rho) return;
    const int N = g->N;
    const long N3 = g->N3;
    const double dV = g->dx * g->dx * g->dx;
    const double inv_dV = 1.0 / dV;
    memset(g->lock_rho, 0, (size_t)N3 * sizeof(double));
    for (int p = 0; p < g->n_locks; p++) {
        ScpLock *L = &g->locks[p];
        if (!L->alive || L->type == 1 || L->q == 0.0) continue;
        double ci, cj, ck;
        locks_phys_to_cell(g, L->x[0], L->x[1], L->x[2], &ci, &cj, &ck);
        int i0 = (int)floor(ci), j0 = (int)floor(cj), k0 = (int)floor(ck);
        double fx = ci - i0, fy = cj - j0, fz = ck - k0;
        double qv = L->q * inv_dV;
        for (int di = 0; di < 2; di++)
        for (int dj = 0; dj < 2; dj++)
        for (int dk = 0; dk < 2; dk++) {
            int ii = locks_wrap_i(i0 + di, N);
            int jj = locks_wrap_i(j0 + dj, N);
            int kk = locks_wrap_i(k0 + dk, N);
            double w = (di ? fx : 1 - fx) * (dj ? fy : 1 - fy) * (dk ? fz : 1 - fz);
            g->lock_rho[locks_idx(ii, jj, kk, N)] += qv * w;
        }
    }
}

/* --- Charge-conserving zigzag / VB current (1D face segments + transverse CIC) ---
 * Discrete continuity: (ρ^{n+1}-ρ^n)/dt + div_h J = 0 with same CIC shape as deposit_rho.
 * J_d lives on the same link as E_d (site → site+d). */

static void locks_shortest_delta_cell(double *x1, double x0, double Nf) {
    double d = *x1 - x0;
    if (d >  0.5 * Nf) *x1 -= Nf;
    if (d < -0.5 * Nf) *x1 += Nf;
}

/* Deposit Jx for motion x0→x1 (cell units) at fixed (y,z) cell coords.
 * Continuity with Cosserat div_h j = Σ(J_d[x]-J_d[x-d])/dx and ρ = q/dV * CIC:
 * charge q*dseg through face ⇒ J * (dx²) * dt = q*dseg ⇒ J = q*dseg/(dt dx²).
 * (Previous port had an extra dx in the numerator — Gauss blew up under motion.) */
static void locks_dep_Jx_seg(Grid *g, double q, double X0, double X1, double Y, double Z, double dt) {
    const int N = g->N;
    const double dx = g->dx;
    if (fabs(X1 - X0) < 1e-18) return;
    double dX = X1 - X0;
    int dir = (dX > 0) ? 1 : -1;
    Y = fmod(Y, (double)N); if (Y < 0) Y += N;
    Z = fmod(Z, (double)N); if (Z < 0) Z += N;
    int j0 = (int)floor(Y), k0 = (int)floor(Z);
    double wy = Y - j0, wz = Z - k0;
    double wj[2] = { 1 - wy, wy }, wk[2] = { 1 - wz, wz };
    double xcur = X0, remaining = dX;
    int guard = 0;
    const double inv_area_dt = 1.0 / (dx * dx * dt);
    while (fabs(remaining) > 1e-15 && guard++ < N + 8) {
        double cell_edge;
        if (dir > 0) {
            cell_edge = floor(xcur + 1e-14) + 1.0;
            if (cell_edge > X1 + 1e-14) cell_edge = X1;
        } else {
            cell_edge = ceil(xcur - 1e-14) - 1.0;
            if (cell_edge < X1 - 1e-14) cell_edge = X1;
        }
        double dseg = cell_edge - xcur;
        if (fabs(dseg) < 1e-18) { xcur += 1e-12 * dir; remaining = X1 - xcur; continue; }
        int iface = locks_wrap_i((int)floor(fmin(xcur, cell_edge) + 1e-14), N);
        double flux = q * dseg * inv_area_dt;
        for (int dj = 0; dj < 2; dj++)
        for (int dk = 0; dk < 2; dk++) {
            int jj = locks_wrap_i(j0 + dj, N);
            int kk = locks_wrap_i(k0 + dk, N);
            g->lock_J[0][locks_idx(iface, jj, kk, N)] += flux * wj[dj] * wk[dk];
        }
        xcur = cell_edge;
        remaining = X1 - xcur;
    }
}

static void locks_dep_Jy_seg(Grid *g, double q, double Y0, double Y1, double X, double Z, double dt) {
    const int N = g->N;
    const double dx = g->dx;
    if (fabs(Y1 - Y0) < 1e-18) return;
    double dY = Y1 - Y0;
    int dir = (dY > 0) ? 1 : -1;
    X = fmod(X, (double)N); if (X < 0) X += N;
    Z = fmod(Z, (double)N); if (Z < 0) Z += N;
    int i0 = (int)floor(X), k0 = (int)floor(Z);
    double wx = X - i0, wz = Z - k0;
    double wi[2] = { 1 - wx, wx }, wk[2] = { 1 - wz, wz };
    double ycur = Y0, remaining = dY;
    int guard = 0;
    const double inv_area_dt = 1.0 / (dx * dx * dt);
    while (fabs(remaining) > 1e-15 && guard++ < N + 8) {
        double cell_edge;
        if (dir > 0) {
            cell_edge = floor(ycur + 1e-14) + 1.0;
            if (cell_edge > Y1 + 1e-14) cell_edge = Y1;
        } else {
            cell_edge = ceil(ycur - 1e-14) - 1.0;
            if (cell_edge < Y1 - 1e-14) cell_edge = Y1;
        }
        double dseg = cell_edge - ycur;
        if (fabs(dseg) < 1e-18) { ycur += 1e-12 * dir; remaining = Y1 - ycur; continue; }
        int jface = locks_wrap_i((int)floor(fmin(ycur, cell_edge) + 1e-14), N);
        double flux = q * dseg * inv_area_dt;
        for (int di = 0; di < 2; di++)
        for (int dk = 0; dk < 2; dk++) {
            int ii = locks_wrap_i(i0 + di, N);
            int kk = locks_wrap_i(k0 + dk, N);
            g->lock_J[1][locks_idx(ii, jface, kk, N)] += flux * wi[di] * wk[dk];
        }
        ycur = cell_edge;
        remaining = Y1 - ycur;
    }
}

static void locks_dep_Jz_seg(Grid *g, double q, double Z0, double Z1, double X, double Y, double dt) {
    const int N = g->N;
    const double dx = g->dx;
    if (fabs(Z1 - Z0) < 1e-18) return;
    double dZ = Z1 - Z0;
    int dir = (dZ > 0) ? 1 : -1;
    X = fmod(X, (double)N); if (X < 0) X += N;
    Y = fmod(Y, (double)N); if (Y < 0) Y += N;
    int i0 = (int)floor(X), j0 = (int)floor(Y);
    double wx = X - i0, wy = Y - j0;
    double wi[2] = { 1 - wx, wx }, wj[2] = { 1 - wy, wy };
    double zcur = Z0, remaining = dZ;
    int guard = 0;
    const double inv_area_dt = 1.0 / (dx * dx * dt);
    while (fabs(remaining) > 1e-15 && guard++ < N + 8) {
        double cell_edge;
        if (dir > 0) {
            cell_edge = floor(zcur + 1e-14) + 1.0;
            if (cell_edge > Z1 + 1e-14) cell_edge = Z1;
        } else {
            cell_edge = ceil(zcur - 1e-14) - 1.0;
            if (cell_edge < Z1 - 1e-14) cell_edge = Z1;
        }
        double dseg = cell_edge - zcur;
        if (fabs(dseg) < 1e-18) { zcur += 1e-12 * dir; remaining = Z1 - zcur; continue; }
        int kface = locks_wrap_i((int)floor(fmin(zcur, cell_edge) + 1e-14), N);
        double flux = q * dseg * inv_area_dt;
        for (int di = 0; di < 2; di++)
        for (int dj = 0; dj < 2; dj++) {
            int ii = locks_wrap_i(i0 + di, N);
            int jj = locks_wrap_i(j0 + dj, N);
            g->lock_J[2][locks_idx(ii, jj, kface, N)] += flux * wi[di] * wj[dj];
        }
        zcur = cell_edge;
        remaining = Z1 - zcur;
    }
}

/* Zigzag path x_prev → x in physical coords; charge-conserving with CIC ρ. */
static void locks_deposit_J_esirkepov(Grid *g) {
    if (g->n_locks <= 0 || !g->lock_J[0]) return;
    const int N = g->N;
    const long N3 = g->N3;
    const double dt = g->dt;
    const double Nf = (double)N;
    for (int d = 0; d < 3; d++)
        memset(g->lock_J[d], 0, (size_t)N3 * sizeof(double));
    for (int p = 0; p < g->n_locks; p++) {
        ScpLock *L = &g->locks[p];
        if (!L->alive || L->pinned || L->type == 1 || L->q == 0.0) continue;
        double c0[3], c1[3];
        locks_phys_to_cell(g, L->x_prev[0], L->x_prev[1], L->x_prev[2], &c0[0], &c0[1], &c0[2]);
        locks_phys_to_cell(g, L->x[0], L->x[1], L->x[2], &c1[0], &c1[1], &c1[2]);
        for (int d = 0; d < 3; d++)
            locks_shortest_delta_cell(&c1[d], c0[d], Nf);
        /* path: x then y then z */
        locks_dep_Jx_seg(g, L->q, c0[0], c1[0], c0[1], c0[2], dt);
        locks_dep_Jy_seg(g, L->q, c0[1], c1[1], c1[0], c0[2], dt);
        locks_dep_Jz_seg(g, L->q, c0[2], c1[2], c1[0], c1[1], dt);
    }
}

/* Apply Ampère for just-completed lock motion.
 * Cosserat convention: ∂t E includes +g J_lat with ∂t ρ = div J_lat (SPEC Gauss-
 * conserving current). Physical PIC continuity is ∂t ρ + div j = 0, so J_lat = −j.
 * Deposit physical j (trajectory zigzag), then E += dt * g * (−j) = −dt g j. */
static void locks_apply_J_to_E(Grid *g, double dt) {
    if (g->n_locks <= 0 || !g->Efield[0] || !g->lock_J[0]) return;
    locks_deposit_J_esirkepov(g);
    const long N3 = g->N3;
    const double fac = -g->g_gauge * dt;   /* −g * j_phys */
    #pragma omp parallel for schedule(static)
    for (long i = 0; i < N3; i++) {
        g->Efield[0][i] += fac * g->lock_J[0][i];
        g->Efield[1][i] += fac * g->lock_J[1][i];
        g->Efield[2][i] += fac * g->lock_J[2][i];
    }
}

/* Gather E on links (half-cell stagger) + B from plaquette sines at cell centers. */
static void locks_gather_EB(const Grid *g, double x, double y, double z,
                            double E[3], double B[3]) {
    const int N = g->N;
    const double GG = g->g_gauge;
    const double inv_ga2 = (GG != 0.0) ? 1.0 / (GG * g->dx * g->dx) : 0.0;
    double ci, cj, ck;
    locks_phys_to_cell(g, x, y, z, &ci, &cj, &ck);

    /* E_d lives on link from site toward +d → sample at ci - 0.5 in direction d */
    for (int d = 0; d < 3; d++) {
        double s[3] = { ci, cj, ck };
        s[d] -= 0.5;
        if (s[d] < 0) s[d] += N;
        int i0 = (int)floor(s[0]), j0 = (int)floor(s[1]), k0 = (int)floor(s[2]);
        double fx = s[0] - i0, fy = s[1] - j0, fz = s[2] - k0;
        double acc = 0;
        for (int di = 0; di < 2; di++)
        for (int dj = 0; dj < 2; dj++)
        for (int dk = 0; dk < 2; dk++) {
            int ii = locks_wrap_i(i0 + di, N);
            int jj = locks_wrap_i(j0 + dj, N);
            int kk = locks_wrap_i(k0 + dk, N);
            double w = (di ? fx : 1 - fx) * (dj ? fy : 1 - fy) * (dk ? fz : 1 - fz);
            acc += w * g->Efield[d][locks_idx(ii, jj, kk, N)];
        }
        E[d] = acc;
    }

    /* B from plaquette: p0=(0,1)->Bz, p1=(1,2)->Bx, p2=(2,0)->By; B ≈ sinθ_P/(g a²) */
    {
        int i0 = (int)floor(ci), j0 = (int)floor(cj), k0 = (int)floor(ck);
        double fx = ci - i0, fy = cj - j0, fz = ck - k0;
        double bx = 0, by = 0, bz = 0;
        if (g->plaq_s[0] && inv_ga2 != 0.0) {
            for (int di = 0; di < 2; di++)
            for (int dj = 0; dj < 2; dj++)
            for (int dk = 0; dk < 2; dk++) {
                int ii = locks_wrap_i(i0 + di, N);
                int jj = locks_wrap_i(j0 + dj, N);
                int kk = locks_wrap_i(k0 + dk, N);
                long id = locks_idx(ii, jj, kk, N);
                double w = (di ? fx : 1 - fx) * (dj ? fy : 1 - fy) * (dk ? fz : 1 - fz);
                /* plaq_s may be stale mid-step; use th if plaq not ready — prefer plaq_s */
                bz += w * g->plaq_s[0][id] * inv_ga2;
                bx += w * g->plaq_s[1][id] * inv_ga2;
                by += w * g->plaq_s[2][id] * inv_ga2;
            }
        }
        B[0] = bx; B[1] = by; B[2] = bz;
    }
}

/* Refresh plaquette sines (needed if push runs without force pass). */
static void locks_refresh_plaquettes(Grid *g) {
    if (g->n_locks <= 0 || !g->th[0] || !g->plaq_s[0]) return;
    const int N = g->N, NN = N * N;
    const long N3 = g->N3;
    #pragma omp parallel for schedule(static)
    for (long idx = 0; idx < N3; idx++) {
        int i = (int)(idx / NN), j = (int)((idx / N) % N), k = (int)(idx % N);
        long np[3];
        np[0] = (long)((i + 1) % N) * NN + (long)j * N + k;
        np[1] = (long)i * NN + (long)((j + 1) % N) * N + k;
        np[2] = (long)i * NN + (long)j * N + (k + 1) % N;
        const int pa[3] = {0, 1, 2}, pb[3] = {1, 2, 0};
        for (int p = 0; p < 3; p++) {
            int a = pa[p], b = pb[p];
            double ang = g->th[a][idx] + g->th[b][np[a]]
                       - g->th[a][np[b]] - g->th[b][idx];
            g->plaq_s[p][idx] = sin(ang);
        }
    }
}

/* Soft form-factor core: charge–charge short-range repulsion (regularization). */
static void locks_soft_core_force(Grid *g, int p, double Fsc[3]) {
    Fsc[0] = Fsc[1] = Fsc[2] = 0;
    if (g->lock_soft_r <= 0 || g->lock_soft_k <= 0) return;
    ScpLock *A = &g->locks[p];
    if (A->type != 0) return; /* only charge locks */
    const double r0 = g->lock_soft_r;
    const double k = g->lock_soft_k;
    const double box = 2.0 * g->L;
    for (int q = 0; q < g->n_locks; q++) {
        if (q == p || !g->locks[q].alive || g->locks[q].type != 0) continue;
        ScpLock *B = &g->locks[q];
        double dx = A->x[0] - B->x[0];
        double dy = A->x[1] - B->x[1];
        double dz = A->x[2] - B->x[2];
        if (dx >  0.5 * box) dx -= box; if (dx < -0.5 * box) dx += box;
        if (dy >  0.5 * box) dy -= box; if (dy < -0.5 * box) dy += box;
        if (dz >  0.5 * box) dz -= box; if (dz < -0.5 * box) dz += box;
        double r2 = dx*dx + dy*dy + dz*dz;
        if (r2 < 1e-20 || r2 >= r0 * r0) continue;
        double r = sqrt(r2);
        double mag = k * (1.0 / r - 1.0 / r0); /* repel */
        Fsc[0] += mag * dx / r;
        Fsc[1] += mag * dy / r;
        Fsc[2] += mag * dz / r;
    }
}

/* Rebuild grid bag B(x): CIC deposit from type=1 anti-locks (weight ∝ E_star),
 * optional Jacobi smooth. Co-field — NOT Cosserat multiplet Higgs. */
static void locks_rebuild_grid_bag(Grid *g) {
    if (g->lock_bag_mode != 2 || !g->lock_bag || g->n_locks <= 0) return;
    const int N = g->N;
    const long N3 = g->N3;
    const double dV = g->dx * g->dx * g->dx;
    const double inv_dV = 1.0 / dV;
    const double sig = (g->lock_bag_r > 0) ? g->lock_bag_r : (2.0 * g->dx);
    const double inv2s2 = 1.0 / (2.0 * sig * sig);
    const int R = (int)ceil(3.0 * sig / g->dx);
    memset(g->lock_bag, 0, (size_t)N3 * sizeof(double));
    for (int p = 0; p < g->n_locks; p++) {
        ScpLock *L = &g->locks[p];
        if (!L->alive || L->type != 1) continue;
        double wtot = (L->E_star > 0) ? L->E_star : 1.0;
        double ci, cj, ck;
        locks_phys_to_cell(g, L->x[0], L->x[1], L->x[2], &ci, &cj, &ck);
        int i0 = (int)floor(ci + 1e-12);
        int j0 = (int)floor(cj + 1e-12);
        int k0 = (int)floor(ck + 1e-12);
        for (int di = -R; di <= R; di++)
        for (int dj = -R; dj <= R; dj++)
        for (int dk = -R; dk <= R; dk++) {
            int ii = locks_wrap_i(i0 + di, N);
            int jj = locks_wrap_i(j0 + dj, N);
            int kk = locks_wrap_i(k0 + dk, N);
            double dx = (di - (ci - i0)) * g->dx;
            double dy = (dj - (cj - j0)) * g->dx;
            double dz = (dk - (ck - k0)) * g->dx;
            double r2 = dx*dx + dy*dy + dz*dz;
            g->lock_bag[locks_idx(ii, jj, kk, N)] += wtot * inv_dV * exp(-r2 * inv2s2);
        }
    }
    int ns = g->lock_bag_smooth;
    if (ns < 0) ns = 0;
    if (ns > 0) {
        double *tmp = (double *)malloc((size_t)N3 * sizeof(double));
        if (!tmp) return;
        for (int s = 0; s < ns; s++) {
            #pragma omp parallel for schedule(static)
            for (long idx = 0; idx < N3; idx++) {
                int i = (int)(idx / (N * N)), j = (int)((idx / N) % N), k = (int)(idx % N);
                long ip = locks_idx(locks_wrap_i(i+1,N), j, k, N);
                long im = locks_idx(locks_wrap_i(i-1,N), j, k, N);
                long jp = locks_idx(i, locks_wrap_i(j+1,N), k, N);
                long jm = locks_idx(i, locks_wrap_i(j-1,N), k, N);
                long kp = locks_idx(i, j, locks_wrap_i(k+1,N), N);
                long km = locks_idx(i, j, locks_wrap_i(k-1,N), N);
                tmp[idx] = 0.5 * g->lock_bag[idx]
                    + (1.0/12.0) * (g->lock_bag[ip]+g->lock_bag[im]+g->lock_bag[jp]
                                   +g->lock_bag[jm]+g->lock_bag[kp]+g->lock_bag[km]);
            }
            memcpy(g->lock_bag, tmp, (size_t)N3 * sizeof(double));
        }
        free(tmp);
    }
}

/* Gather −κ ∇B at physical (x,y,z) via CIC of central differences on B. */
static void locks_gather_grad_bag(const Grid *g, double x, double y, double z, double gB[3]) {
    gB[0] = gB[1] = gB[2] = 0;
    if (!g->lock_bag) return;
    const int N = g->N;
    const double inv2dx = 0.5 / g->dx;
    double ci, cj, ck;
    locks_phys_to_cell(g, x, y, z, &ci, &cj, &ck);
    int i0 = (int)floor(ci), j0 = (int)floor(cj), k0 = (int)floor(ck);
    double fx = ci - i0, fy = cj - j0, fz = ck - k0;
    double gx = 0, gy = 0, gz = 0;
    for (int di = 0; di < 2; di++)
    for (int dj = 0; dj < 2; dj++)
    for (int dk = 0; dk < 2; dk++) {
        int ii = locks_wrap_i(i0 + di, N);
        int jj = locks_wrap_i(j0 + dj, N);
        int kk = locks_wrap_i(k0 + dk, N);
        double w = (di ? fx : 1 - fx) * (dj ? fy : 1 - fy) * (dk ? fz : 1 - fz);
        long ip = locks_idx(locks_wrap_i(ii+1,N), jj, kk, N);
        long im = locks_idx(locks_wrap_i(ii-1,N), jj, kk, N);
        long jp = locks_idx(ii, locks_wrap_i(jj+1,N), kk, N);
        long jm = locks_idx(ii, locks_wrap_i(jj-1,N), kk, N);
        long kp = locks_idx(ii, jj, locks_wrap_i(kk+1,N), N);
        long km = locks_idx(ii, jj, locks_wrap_i(kk-1,N), N);
        gx += w * (g->lock_bag[ip] - g->lock_bag[im]) * inv2dx;
        gy += w * (g->lock_bag[jp] - g->lock_bag[jm]) * inv2dx;
        gz += w * (g->lock_bag[kp] - g->lock_bag[km]) * inv2dx;
    }
    gB[0] = gx; gB[1] = gy; gB[2] = gz;
}

/* Bag force on charge lock p.
 * mode 1: pairwise anti-lock (legacy form-factor).
 * mode 2: F = −k ∇B from grid co-field (rebuild bag first outside). */
static void locks_bag_force(Grid *g, int p, double Fbag[3]) {
    Fbag[0] = Fbag[1] = Fbag[2] = 0;
    if (g->lock_bag_k <= 0) return;
    int mode = g->lock_bag_mode;
    /* back-compat: bag_k+r set without mode → pairwise */
    if (mode == 0 && g->lock_bag_r > 0) mode = 1;
    if (mode <= 0) return;
    ScpLock *A = &g->locks[p];
    if (A->type != 0) return;
    const double k = g->lock_bag_k;

    if (mode == 2) {
        double gB[3];
        locks_gather_grad_bag(g, A->x[0], A->x[1], A->x[2], gB);
        Fbag[0] = -k * gB[0];
        Fbag[1] = -k * gB[1];
        Fbag[2] = -k * gB[2];
        return;
    }

    /* mode 1: pairwise */
    if (g->lock_bag_r <= 0) return;
    const double r0 = g->lock_bag_r;
    const double box = 2.0 * g->L;
    for (int q = 0; q < g->n_locks; q++) {
        if (!g->locks[q].alive || g->locks[q].type != 1) continue;
        ScpLock *B = &g->locks[q];
        double dx = A->x[0] - B->x[0];
        double dy = A->x[1] - B->x[1];
        double dz = A->x[2] - B->x[2];
        if (dx >  0.5 * box) dx -= box; if (dx < -0.5 * box) dx += box;
        if (dy >  0.5 * box) dy -= box; if (dy < -0.5 * box) dy += box;
        if (dz >  0.5 * box) dz -= box; if (dz < -0.5 * box) dz += box;
        double r2 = dx*dx + dy*dy + dz*dz;
        if (r2 < 1e-20 || r2 >= r0 * r0) continue;
        double r = sqrt(r2);
        double mag = -k * (1.0 / r - 1.0 / r0);
        Fbag[0] += mag * dx / r;
        Fbag[1] += mag * dy / r;
        Fbag[2] += mag * dz / r;
    }
}

/* Relativistic Boris (charge) / force push (anti-lock); soft core + bag. */
static void locks_push_and_move(Grid *g, double dt) {
    if (g->n_locks <= 0) return;
    const double GG = g->g_gauge;
    const double box = 2.0 * g->L;

    locks_refresh_plaquettes(g);
    if (g->lock_bag_mode == 2)
        locks_rebuild_grid_bag(g);

    for (int p = 0; p < g->n_locks; p++) {
        ScpLock *L = &g->locks[p];
        if (!L->alive) continue;

        L->x_prev[0] = L->x[0];
        L->x_prev[1] = L->x[1];
        L->x_prev[2] = L->x[2];

        double E[3] = {0,0,0}, B[3] = {0,0,0}, Fsc[3], Fbag[3];
        locks_soft_core_force(g, p, Fsc);
        locks_bag_force(g, p, Fbag);

        if (L->pinned) {
            L->u[0] = L->u[1] = L->u[2] = 0;
            L->f[0] = Fsc[0] + Fbag[0];
            L->f[1] = Fsc[1] + Fbag[1];
            L->f[2] = Fsc[2] + Fbag[2];
            continue;
        }

        /* --- anti-lock: no EM; optional free motion under zero force (or pin) --- */
        if (L->type == 1 || L->q == 0.0) {
            /* bag does not self-force; anti-locks stay put unless given velocity */
            L->f[0] = L->f[1] = L->f[2] = 0;
            double u2 = L->u[0]*L->u[0] + L->u[1]*L->u[1] + L->u[2]*L->u[2];
            double gamma = sqrt(1.0 + u2);
            L->x[0] += dt * L->u[0] / gamma;
            L->x[1] += dt * L->u[1] / gamma;
            L->x[2] += dt * L->u[2] / gamma;
            for (int d = 0; d < 3; d++) {
                while (L->x[d] >  g->L) L->x[d] -= box;
                while (L->x[d] <= -g->L) L->x[d] += box;
            }
            continue;
        }

        locks_gather_EB(g, L->x[0], L->x[1], L->x[2], E, B);

        /* F_em = −(g q) E; plus soft core + bag pull from anti-locks */
        double qem = -GG * L->q;
        double Ftot[3] = {
            qem * E[0] + Fsc[0] + Fbag[0],
            qem * E[1] + Fsc[1] + Fbag[1],
            qem * E[2] + Fsc[2] + Fbag[2]
        };
        double Eeff[3] = {
            E[0] + (Fsc[0] + Fbag[0]) / qem,
            E[1] + (Fsc[1] + Fbag[1]) / qem,
            E[2] + (Fsc[2] + Fbag[2]) / qem
        };
        L->f[0] = Ftot[0]; L->f[1] = Ftot[1]; L->f[2] = Ftot[2];

        double qom = qem / L->m;
        double half = 0.5 * dt * qom;
        double ux = L->u[0] + half * Eeff[0];
        double uy = L->u[1] + half * Eeff[1];
        double uz = L->u[2] + half * Eeff[2];
        double gamma = sqrt(1.0 + ux*ux + uy*uy + uz*uz);
        double tx = half * B[0] / gamma;
        double ty = half * B[1] / gamma;
        double tz = half * B[2] / gamma;
        double uxp = ux + (uy * tz - uz * ty);
        double uyp = uy + (uz * tx - ux * tz);
        double uzp = uz + (ux * ty - uy * tx);
        double t2 = tx*tx + ty*ty + tz*tz;
        double sx = 2.0 * tx / (1.0 + t2);
        double sy = 2.0 * ty / (1.0 + t2);
        double sz = 2.0 * tz / (1.0 + t2);
        ux = ux + (uyp * sz - uzp * sy);
        uy = uy + (uzp * sx - uxp * sz);
        uz = uz + (uxp * sy - uyp * sx);
        ux += half * Eeff[0];
        uy += half * Eeff[1];
        uz += half * Eeff[2];
        L->u[0] = ux; L->u[1] = uy; L->u[2] = uz;

        gamma = sqrt(1.0 + ux*ux + uy*uy + uz*uz);
        double vx = ux / gamma, vy = uy / gamma, vz = uz / gamma;
        L->f[0] = qem * (E[0] + vy * B[2] - vz * B[1]) + Fsc[0] + Fbag[0];
        L->f[1] = qem * (E[1] + vz * B[0] - vx * B[2]) + Fsc[1] + Fbag[1];
        L->f[2] = qem * (E[2] + vx * B[1] - vy * B[0]) + Fsc[2] + Fbag[2];

        L->x[0] += dt * vx;
        L->x[1] += dt * vy;
        L->x[2] += dt * vz;
        for (int d = 0; d < 3; d++) {
            while (L->x[d] >  g->L) L->x[d] -= box;
            while (L->x[d] <= -g->L) L->x[d] += box;
        }
    }
    locks_deposit_rho(g);
}

static void locks_energy(const Grid *g, double *E_star, double *E_kin, int *alive) {
    double es = 0, ek = 0;
    int na = 0;
    for (int p = 0; p < g->n_locks; p++) {
        const ScpLock *L = &g->locks[p];
        if (!L->alive) continue;
        na++;
        es += L->E_star;
        double u2 = L->u[0]*L->u[0] + L->u[1]*L->u[1] + L->u[2]*L->u[2];
        double gamma = sqrt(1.0 + u2);
        ek += L->m * (gamma - 1.0);
    }
    if (E_star) *E_star = es;
    if (E_kin) *E_kin = ek;
    if (alive) *alive = na;
}

static void locks_write_track(Grid *g, double t, FILE *fp) {
    if (!fp || g->n_locks <= 0) return;
    for (int p = 0; p < g->n_locks; p++) {
        ScpLock *L = &g->locks[p];
        fprintf(fp, "%.6f\t%d\t%d\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t"
                    "%.12e\t%.12e\t%.12e\t%d\t%d\n",
                t, L->id, L->type, L->x[0], L->x[1], L->x[2],
                L->u[0], L->u[1], L->u[2],
                L->f[0], L->f[1], L->f[2],
                L->pinned, L->alive);
    }
    fflush(fp);
}

/* Total lock charge (sum q) for net-charge checks. */
static double locks_Q_total(const Grid *g) {
    double Q = 0;
    for (int p = 0; p < g->n_locks; p++)
        if (g->locks[p].alive) Q += g->locks[p].q;
    return Q;
}

#endif /* SCP_LOCKS_H */
