/*
 * hopfion_init.c — Initialization routines for hopfion/skyrmion configurations
 *
 * Key capability: compose multiple topological objects on the same grid.
 * The product ansatz q₁×q₂/ρ₀ preserves topology and creates linking.
 */

#include "hopfion_init.h"
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <omp.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ========== Profile I/O ========== */

RadialProfile *profile_load(const char *filename)
{
    FILE *fp = fopen(filename, "r");
    if (!fp) { fprintf(stderr, "Cannot open %s\n", filename); return NULL; }

    int has_rho = 0;
    int n = 0;
    char line[1024];
    while (fgets(line, sizeof(line), fp)) {
        if (line[0] == '#') {
            if (strstr(line, "rho(r)") || strstr(line, "Finite-lambda"))
                has_rho = 1;
            continue;
        }
        double rv, fv;
        if (sscanf(line, "%lf %lf", &rv, &fv) >= 2) n++;
    }
    rewind(fp);

    RadialProfile *p = malloc(sizeof(RadialProfile));
    p->r = malloc(n * sizeof(double));
    p->f = malloc(n * sizeof(double));
    p->rho = has_rho ? malloc(n * sizeof(double)) : NULL;
    p->n = n;

    int i = 0;
    while (fgets(line, sizeof(line), fp) && i < n) {
        if (line[0] == '#') continue;
        double rv, fv, fpv, rhov;
        int ncols = sscanf(line, "%lf %lf %lf %lf", &rv, &fv, &fpv, &rhov);
        if (ncols >= 2) {
            p->r[i] = rv;
            p->f[i] = fv;
            if (has_rho && ncols >= 4) p->rho[i] = rhov;
            i++;
        }
    }
    fclose(fp);

    p->dr = (n > 1) ? p->r[1] - p->r[0] : 0.001;
    p->r_max = p->r[n-1];

    printf("Loaded profile: %s (%d points, r=[0, %.1f]%s)\n",
           filename, n, p->r_max, has_rho ? ", finite-lambda" : ", sigma-model");
    return p;
}

void profile_free(RadialProfile *p)
{
    if (p) { free(p->r); free(p->f); if (p->rho) free(p->rho); free(p); }
}

double profile_interp_f(const RadialProfile *p, double r)
{
    if (r < 0) r = -r;
    if (r >= p->r_max) return 0.0;
    double fi = r / p->dr;
    int i = (int)fi;
    if (i >= p->n - 1) return 0.0;
    double t = fi - i;
    return (1-t)*p->f[i] + t*p->f[i+1];
}

double profile_interp_rho(const RadialProfile *p, double rho0, double r)
{
    if (!p->rho) return rho0;
    if (r < 0) r = -r;
    if (r >= p->r_max) return rho0;
    double fi = r / p->dr;
    int i = (int)fi;
    if (i >= p->n - 1) return rho0;
    double t = fi - i;
    return (1-t)*p->rho[i] + t*p->rho[i+1];
}

/* ========== Hedgehog Skyrmion ========== */

void init_skyrmion(SphericalGrid *g, const RadialProfile *prof, double rho0,
                   double cx, double cy, double cz)
{
    printf("Initializing B=1 hedgehog Skyrmion at (%.2f, %.2f, %.2f)\n", cx, cy, cz);

    #pragma omp parallel for schedule(static)
    for (int n = 0; n < g->n_active; n++) {
        int ix = g->active[n];
        int i, j, k;
        sg_unflatten(g->N, ix, &i, &j, &k);

        double x, y, z;
        sg_pos(g, i, j, k, &x, &y, &z);

        double dx = x - cx, dy = y - cy, dz = z - cz;
        double r = sqrt(dx*dx + dy*dy + dz*dz);
        double f = profile_interp_f(prof, r);
        double rho = profile_interp_rho(prof, rho0, r);

        double cf = cos(f), sf = sin(f);
        g->psi[ix].s = rho * cf;
        if (r > 1e-12) {
            double sr = rho * sf / r;
            g->psi[ix].f1 = sr * dx;
            g->psi[ix].f2 = sr * dy;
            g->psi[ix].f3 = sr * dz;
        } else {
            g->psi[ix].f1 = 0;
            g->psi[ix].f2 = 0;
            g->psi[ix].f3 = 0;
        }
        /* Weight sector = 0 */
        g->psi[ix].j1 = 0;
        g->psi[ix].j2 = 0;
        g->psi[ix].j3 = 0;
        g->psi[ix].p  = 0;
    }
}

/* ========== Hopfion (H=1 axially symmetric) ========== */

/* The Hopf map from R³ to S² with size parameter a:
 *   z₁ = (a² - r² + 2iaz) / (a² + r²)
 *   z₂ = 2a(x + iy) / (a² + r²)
 *   n̂ = (2Re(z₁z̄₂), 2Im(z₁z̄₂), |z₁|²-|z₂|²)
 *
 * The quaternion field from n̂ uses a "hedgehog-like" map:
 *   q = ρ₀ [cos(g/2) + sin(g/2)(n̂·σ)]
 * where g(r) is a profile function controlling the radial envelope.
 *
 * For the Faddeev-Skyrme model, g(r) satisfies a specific ODE.
 * For initialization, we use a Gaussian-like profile:
 *   g(r) = π × exp(-(r/a)⁴)    [compact support approximation]
 * This gives g(0) = π (anti-vacuum at center) and g(∞) = 0 (vacuum).
 */
void init_hopfion(SphericalGrid *g, double rho0, double a,
                  double cx, double cy, double cz)
{
    printf("Initializing H=1 hopfion at (%.2f, %.2f, %.2f), a=%.2f\n", cx, cy, cz, a);

    #pragma omp parallel for schedule(static)
    for (int n = 0; n < g->n_active; n++) {
        int ix = g->active[n];
        int ii, jj, kk;
        sg_unflatten(g->N, ix, &ii, &jj, &kk);

        double xp, yp, zp;
        sg_pos(g, ii, jj, kk, &xp, &yp, &zp);

        double x = xp - cx, y = yp - cy, z = zp - cz;
        double r2 = x*x + y*y + z*z;
        double denom = a*a + r2;

        /* Hopf map: compute z₁, z₂ ∈ C */
        double z1_re = (a*a - r2) / denom;
        double z1_im = 2*a*z / denom;
        double z2_re = 2*a*x / denom;
        double z2_im = 2*a*y / denom;

        /* n̂ = stereographic projection from (z₁,z₂) */
        double n1 = 2*(z1_re*z2_re + z1_im*z2_im);          /* 2 Re(z₁z̄₂) */
        double n2 = 2*(z1_im*z2_re - z1_re*z2_im);          /* 2 Im(z₁z̄₂) */
        double n3 = z1_re*z1_re + z1_im*z1_im - z2_re*z2_re - z2_im*z2_im;

        /* Normalize (should already be unit, but numerics) */
        double nn = sqrt(n1*n1 + n2*n2 + n3*n3);
        if (nn > 1e-15) { n1 /= nn; n2 /= nn; n3 /= nn; }
        else { n1 = 0; n2 = 0; n3 = 1; }

        /* Radial envelope: g(r) controls how much rotation is applied */
        double r = sqrt(r2);
        double ra = r / a;
        double gval = M_PI * exp(-ra*ra*ra*ra);  /* compact profile */

        /* Quaternion from n̂ and angle g:
         *   q = ρ₀ [cos(g/2) + sin(g/2)(n₁σ₁ + n₂σ₂ + n₃σ₃)] */
        double cg = cos(0.5 * gval);
        double sg_val = sin(0.5 * gval);

        g->psi[ix].s  = rho0 * cg;
        g->psi[ix].f1 = rho0 * sg_val * n1;
        g->psi[ix].f2 = rho0 * sg_val * n2;
        g->psi[ix].f3 = rho0 * sg_val * n3;
        g->psi[ix].j1 = 0;
        g->psi[ix].j2 = 0;
        g->psi[ix].j3 = 0;
        g->psi[ix].p  = 0;
    }
}

/* ========== (n,m) Hopfion ========== */

/* The (n,m) torus knot hopfion uses the parameterization:
 *   The preimage of the south pole traces out a (n,m)-torus knot.
 *   Hopf index = n × m.
 *
 * Construction: Use the modified Hopf map with winding numbers:
 *   z₁ = ((a²-r²+2iaz) / (a²+r²))^n
 *   z₂ = (2a(x+iy) / (a²+r²))^m      [note: power m, not 1]
 *   normalize: (z₁,z₂) → (z₁,z₂)/|(z₁,z₂)|
 *   then n̂ = standard S² projection
 */
void init_hopfion_nm(SphericalGrid *g, double rho0, double a,
                     int n_tor, int m_pol,
                     double cx, double cy, double cz)
{
    printf("Initializing (%d,%d)-hopfion (H=%d) at (%.2f, %.2f, %.2f), a=%.2f\n",
           n_tor, m_pol, n_tor*m_pol, cx, cy, cz, a);

    #pragma omp parallel for schedule(static)
    for (int n = 0; n < g->n_active; n++) {
        int ix = g->active[n];
        int ii, jj, kk;
        sg_unflatten(g->N, ix, &ii, &jj, &kk);

        double xp, yp, zp;
        sg_pos(g, ii, jj, kk, &xp, &yp, &zp);

        double x = xp - cx, y = yp - cy, z = zp - cz;
        double r2 = x*x + y*y + z*z;
        double denom = a*a + r2;

        /* Base Hopf map components */
        double w1_re = (a*a - r2) / denom;
        double w1_im = 2*a*z / denom;
        double w2_re = 2*a*x / denom;
        double w2_im = 2*a*y / denom;

        /* Raise to powers: z₁ = w₁^n, z₂ = w₂^m */
        /* Complex power: start with 1+0i and multiply n times */
        double z1_re = 1, z1_im = 0;
        for (int p = 0; p < n_tor; p++) {
            double re = z1_re * w1_re - z1_im * w1_im;
            double im = z1_re * w1_im + z1_im * w1_re;
            z1_re = re; z1_im = im;
        }

        double z2_re = 1, z2_im = 0;
        for (int p = 0; p < m_pol; p++) {
            double re = z2_re * w2_re - z2_im * w2_im;
            double im = z2_re * w2_im + z2_im * w2_re;
            z2_re = re; z2_im = im;
        }

        /* Normalize to unit sphere in C² */
        double norm = sqrt(z1_re*z1_re + z1_im*z1_im + z2_re*z2_re + z2_im*z2_im);
        if (norm > 1e-15) {
            z1_re /= norm; z1_im /= norm;
            z2_re /= norm; z2_im /= norm;
        }

        /* S² projection */
        double n1 = 2*(z1_re*z2_re + z1_im*z2_im);
        double n2 = 2*(z1_im*z2_re - z1_re*z2_im);
        double n3 = z1_re*z1_re + z1_im*z1_im - z2_re*z2_re - z2_im*z2_im;

        double nn = sqrt(n1*n1 + n2*n2 + n3*n3);
        if (nn > 1e-15) { n1 /= nn; n2 /= nn; n3 /= nn; }
        else { n1 = 0; n2 = 0; n3 = 1; }

        /* Radial envelope */
        double r = sqrt(r2);
        double ra = r / a;
        double gval = M_PI * exp(-ra*ra*ra*ra);

        double cg = cos(0.5 * gval);
        double sg_val = sin(0.5 * gval);

        g->psi[ix].s  = rho0 * cg;
        g->psi[ix].f1 = rho0 * sg_val * n1;
        g->psi[ix].f2 = rho0 * sg_val * n2;
        g->psi[ix].f3 = rho0 * sg_val * n3;
        g->psi[ix].j1 = 0;
        g->psi[ix].j2 = 0;
        g->psi[ix].j3 = 0;
        g->psi[ix].p  = 0;
    }
}

/* ========== Composed states ========== */

void init_two_skyrmions(SphericalGrid *g, const RadialProfile *prof, double rho0,
                        double z0, double vz)
{
    printf("Initializing two B=1 Skyrmions at z=±%.2f, v=±%.3fc\n", z0, vz);

    /* Use product ansatz: q = q₁(x-z₀ẑ) × q₂(x+z₀ẑ) / ρ₀ */
    #pragma omp parallel for schedule(static)
    for (int n = 0; n < g->n_active; n++) {
        int ix = g->active[n];
        int ii, jj, kk;
        sg_unflatten(g->N, ix, &ii, &jj, &kk);

        double x, y, z;
        sg_pos(g, ii, jj, kk, &x, &y, &z);

        /* Soliton 1 at (0,0,+z0) */
        double dx1 = x, dy1 = y, dz1 = z - z0;
        double r1 = sqrt(dx1*dx1 + dy1*dy1 + dz1*dz1);
        double f1 = profile_interp_f(prof, r1);
        double rho1 = profile_interp_rho(prof, rho0, r1);
        double cf1 = cos(f1), sf1 = sin(f1);

        double q1_s = rho1 * cf1;
        double q1_f1 = 0, q1_f2 = 0, q1_f3 = 0;
        if (r1 > 1e-12) {
            double sr1 = rho1 * sf1 / r1;
            q1_f1 = sr1 * dx1;
            q1_f2 = sr1 * dy1;
            q1_f3 = sr1 * dz1;
        }

        /* Soliton 2 at (0,0,-z0) */
        double dx2 = x, dy2 = y, dz2 = z + z0;
        double r2 = sqrt(dx2*dx2 + dy2*dy2 + dz2*dz2);
        double f2 = profile_interp_f(prof, r2);
        double rho2 = profile_interp_rho(prof, rho0, r2);
        double cf2 = cos(f2), sf2 = sin(f2);

        double q2_s = rho2 * cf2;
        double q2_f1 = 0, q2_f2 = 0, q2_f3 = 0;
        if (r2 > 1e-12) {
            double sr2 = rho2 * sf2 / r2;
            q2_f1 = sr2 * dx2;
            q2_f2 = sr2 * dy2;
            q2_f3 = sr2 * dz2;
        }

        /* Product: q = q₁ × q₂ / ρ₀ */
        double inv_rho0 = 1.0 / rho0;
        g->psi[ix].s  = inv_rho0 * (q1_s*q2_s  - q1_f1*q2_f1 - q1_f2*q2_f2 - q1_f3*q2_f3);
        g->psi[ix].f1 = inv_rho0 * (q1_s*q2_f1 + q1_f1*q2_s  - q1_f2*q2_f3 + q1_f3*q2_f2);
        g->psi[ix].f2 = inv_rho0 * (q1_s*q2_f2 + q1_f1*q2_f3 + q1_f2*q2_s  - q1_f3*q2_f1);
        g->psi[ix].f3 = inv_rho0 * (q1_s*q2_f3 - q1_f1*q2_f2 + q1_f2*q2_f1 + q1_f3*q2_s);
        g->psi[ix].j1 = 0;
        g->psi[ix].j2 = 0;
        g->psi[ix].j3 = 0;
        g->psi[ix].p  = 0;

        /* Velocity: v(x) = v_z × ∂q/∂z contribution from each soliton */
        /* dq₁/dz at this point (chain rule through product) */
        if (vz != 0) {
            /* Approximate: v = -v_z ∂(q₁q₂/ρ₀)/∂z for soliton 1 moving +z
             * + v_z ∂(q₁q₂/ρ₀)/∂z for soliton 2 moving -z */
            /* For well-separated solitons, the cross terms are small */
            /* Use finite difference */
            double eps = g->h;
            double z_plus = z + eps, z_minus = z - eps;

            /* Recompute at z±eps for finite-diff velocity */
            /* Soliton 1 contribution */
            double dz1p = z_plus - z0, dz1m = z_minus - z0;
            double r1p = sqrt(dx1*dx1 + dy1*dy1 + dz1p*dz1p);
            double r1m = sqrt(dx1*dx1 + dy1*dy1 + dz1m*dz1m);
            double f1p = profile_interp_f(prof, r1p);
            double f1m = profile_interp_f(prof, r1m);
            double rho1p = profile_interp_rho(prof, rho0, r1p);
            double rho1m = profile_interp_rho(prof, rho0, r1m);

            double dq1_s  = (rho1p*cos(f1p) - rho1m*cos(f1m)) / (2*eps);
            double dq1_f3 = 0;
            if (r1p > 1e-12 && r1m > 1e-12) {
                dq1_f3 = (rho1p*sin(f1p)*dz1p/r1p - rho1m*sin(f1m)*dz1m/r1m) / (2*eps);
            }

            /* Simplified: only z-component of velocity for head-on collision */
            g->vel[ix].s  = -vz * dq1_s * q2_s * inv_rho0
                            + vz * q1_s * dq1_s * inv_rho0;  /* approximate */
            /* Full velocity computation deferred — set to zero for now */
            g->vel[ix] = mv_zero();
        }
    }
}

void init_skyrmion_hopfion(SphericalGrid *g, const RadialProfile *prof,
                            double rho0, double hopfion_a, double hopfion_d)
{
    printf("Initializing Skyrmion + hopfion ring (a=%.2f, d=%.2f)\n",
           hopfion_a, hopfion_d);

    /* Step 1: Initialize Skyrmion at origin */
    init_skyrmion(g, prof, rho0, 0, 0, 0);

    /* Step 2: Compose with hopfion ring displaced by d along z */
    /* We need a temporary copy for the hopfion, then multiply */
    long total = (long)g->N * g->N * g->N;
    Multivector *skyrmion_field = malloc(total * sizeof(Multivector));
    memcpy(skyrmion_field, g->psi, total * sizeof(Multivector));

    /* Initialize hopfion on the grid */
    init_hopfion(g, rho0, hopfion_a, 0, 0, hopfion_d);

    /* Compose: q_total = q_skyrmion × q_hopfion / ρ₀ */
    double inv_rho0 = 1.0 / rho0;
    #pragma omp parallel for schedule(static)
    for (int n = 0; n < g->n_active; n++) {
        int ix = g->active[n];
        Multivector qs = skyrmion_field[ix];
        Multivector qh = g->psi[ix];

        /* Quaternion product (bulk only) */
        g->psi[ix].s  = inv_rho0 * (qs.s*qh.s  - qs.f1*qh.f1 - qs.f2*qh.f2 - qs.f3*qh.f3);
        g->psi[ix].f1 = inv_rho0 * (qs.s*qh.f1 + qs.f1*qh.s  - qs.f2*qh.f3 + qs.f3*qh.f2);
        g->psi[ix].f2 = inv_rho0 * (qs.s*qh.f2 + qs.f1*qh.f3 + qs.f2*qh.s  - qs.f3*qh.f1);
        g->psi[ix].f3 = inv_rho0 * (qs.s*qh.f3 - qs.f1*qh.f2 + qs.f2*qh.f1 + qs.f3*qh.s);
        g->psi[ix].j1 = 0;
        g->psi[ix].j2 = 0;
        g->psi[ix].j3 = 0;
        g->psi[ix].p  = 0;
    }

    free(skyrmion_field);
}

void init_boost(SphericalGrid *g, double rho0,
                double cx, double cy, double cz,
                double vx, double vy, double vz)
{
    /* v(x) ≈ -(v·∇)q for non-relativistic boost
     * Only applies near the soliton center (weight by proximity) */
    (void)rho0; (void)cx; (void)cy; (void)cz;

    #pragma omp parallel for schedule(static)
    for (int n = 0; n < g->n_active; n++) {
        int ix = g->active[n];
        int i, j, k;
        sg_unflatten(g->N, ix, &i, &j, &k);

        double r = sg_radius(g, i, j, k);
        if (r > g->R - 2.5*g->h) continue;

        /* Compute ∂q/∂x_i using 4th-order central differences */
        double c1 = 1.0 / (12.0 * g->h);
        Multivector dqdx = mv_zero(), dqdy = mv_zero(), dqdz = mv_zero();

        /* X */
        Multivector xm2 = sg_get(g, i-2, j, k);
        Multivector xm1 = sg_get(g, i-1, j, k);
        Multivector xp1 = sg_get(g, i+1, j, k);
        Multivector xp2 = sg_get(g, i+2, j, k);

        #define DERIV4(comp) c1*(xm2.comp - 8*xm1.comp + 8*xp1.comp - xp2.comp)
        dqdx.s = DERIV4(s); dqdx.f1 = DERIV4(f1); dqdx.f2 = DERIV4(f2); dqdx.f3 = DERIV4(f3);
        #undef DERIV4

        Multivector ym2 = sg_get(g, i, j-2, k);
        Multivector ym1 = sg_get(g, i, j-1, k);
        Multivector yp1 = sg_get(g, i, j+1, k);
        Multivector yp2 = sg_get(g, i, j+2, k);

        #define DERIV4(comp) c1*(ym2.comp - 8*ym1.comp + 8*yp1.comp - yp2.comp)
        dqdy.s = DERIV4(s); dqdy.f1 = DERIV4(f1); dqdy.f2 = DERIV4(f2); dqdy.f3 = DERIV4(f3);
        #undef DERIV4

        Multivector zm2 = sg_get(g, i, j, k-2);
        Multivector zm1 = sg_get(g, i, j, k-1);
        Multivector zp1 = sg_get(g, i, j, k+1);
        Multivector zp2 = sg_get(g, i, j, k+2);

        #define DERIV4(comp) c1*(zm2.comp - 8*zm1.comp + 8*zp1.comp - zp2.comp)
        dqdz.s = DERIV4(s); dqdz.f1 = DERIV4(f1); dqdz.f2 = DERIV4(f2); dqdz.f3 = DERIV4(f3);
        #undef DERIV4

        /* v(x) = -(vx ∂q/∂x + vy ∂q/∂y + vz ∂q/∂z) */
        g->vel[ix].s  = -(vx*dqdx.s  + vy*dqdy.s  + vz*dqdz.s);
        g->vel[ix].f1 = -(vx*dqdx.f1 + vy*dqdy.f1 + vz*dqdz.f1);
        g->vel[ix].f2 = -(vx*dqdx.f2 + vy*dqdy.f2 + vz*dqdz.f2);
        g->vel[ix].f3 = -(vx*dqdx.f3 + vy*dqdy.f3 + vz*dqdz.f3);
        g->vel[ix].j1 = 0;
        g->vel[ix].j2 = 0;
        g->vel[ix].j3 = 0;
        g->vel[ix].p  = 0;
    }
}
