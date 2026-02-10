/*
 * hopfion_field.c — Energy functional and force for composed hopfion dynamics
 *
 * Key new feature: L₆ = -λ₆ (B⁰)² term where B⁰ is the baryon density.
 * This is the "most topological" local term — it makes the winding density
 * directly dynamical.
 *
 * Derivatives: 4th-order central differences (matching existing hopfion_search code).
 * The 9-point stencil {1,-8,0,8,-1}/(12h) for first derivatives.
 * The consistent Laplacian uses {1,-16,64,16,-130,16,64,-16,1}/(144h²) per axis
 * (exact discrete gradient of discrete E₂ energy).
 *
 * Skyrme force: analytical 3-pass algorithm ported from field.c.
 * L₆ force: numerical finite differences of (B⁰)² energy.
 */

#include "hopfion_field.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ========== Derivative helpers ========== */

/* 4th-order central difference: df/dx_i at (i,j,k)
 * Uses safe access (returns vacuum outside sphere) */
static inline Multivector deriv_x(const SphericalGrid *g, int i, int j, int k)
{
    double c1 = 1.0 / (12.0 * g->h);
    Multivector m2 = sg_get(g, i-2, j, k);
    Multivector m1 = sg_get(g, i-1, j, k);
    Multivector p1 = sg_get(g, i+1, j, k);
    Multivector p2 = sg_get(g, i+2, j, k);

    Multivector r;
    r.s  = c1 * (m2.s  - 8*m1.s  + 8*p1.s  - p2.s);
    r.f1 = c1 * (m2.f1 - 8*m1.f1 + 8*p1.f1 - p2.f1);
    r.f2 = c1 * (m2.f2 - 8*m1.f2 + 8*p1.f2 - p2.f2);
    r.f3 = c1 * (m2.f3 - 8*m1.f3 + 8*p1.f3 - p2.f3);
    r.j1 = c1 * (m2.j1 - 8*m1.j1 + 8*p1.j1 - p2.j1);
    r.j2 = c1 * (m2.j2 - 8*m1.j2 + 8*p1.j2 - p2.j2);
    r.j3 = c1 * (m2.j3 - 8*m1.j3 + 8*p1.j3 - p2.j3);
    r.p  = c1 * (m2.p  - 8*m1.p  + 8*p1.p  - p2.p);
    return r;
}

static inline Multivector deriv_y(const SphericalGrid *g, int i, int j, int k)
{
    double c1 = 1.0 / (12.0 * g->h);
    Multivector m2 = sg_get(g, i, j-2, k);
    Multivector m1 = sg_get(g, i, j-1, k);
    Multivector p1 = sg_get(g, i, j+1, k);
    Multivector p2 = sg_get(g, i, j+2, k);

    Multivector r;
    r.s  = c1 * (m2.s  - 8*m1.s  + 8*p1.s  - p2.s);
    r.f1 = c1 * (m2.f1 - 8*m1.f1 + 8*p1.f1 - p2.f1);
    r.f2 = c1 * (m2.f2 - 8*m1.f2 + 8*p1.f2 - p2.f2);
    r.f3 = c1 * (m2.f3 - 8*m1.f3 + 8*p1.f3 - p2.f3);
    r.j1 = c1 * (m2.j1 - 8*m1.j1 + 8*p1.j1 - p2.j1);
    r.j2 = c1 * (m2.j2 - 8*m1.j2 + 8*p1.j2 - p2.j2);
    r.j3 = c1 * (m2.j3 - 8*m1.j3 + 8*p1.j3 - p2.j3);
    r.p  = c1 * (m2.p  - 8*m1.p  + 8*p1.p  - p2.p);
    return r;
}

static inline Multivector deriv_z(const SphericalGrid *g, int i, int j, int k)
{
    double c1 = 1.0 / (12.0 * g->h);
    Multivector m2 = sg_get(g, i, j, k-2);
    Multivector m1 = sg_get(g, i, j, k-1);
    Multivector p1 = sg_get(g, i, j, k+1);
    Multivector p2 = sg_get(g, i, j, k+2);

    Multivector r;
    r.s  = c1 * (m2.s  - 8*m1.s  + 8*p1.s  - p2.s);
    r.f1 = c1 * (m2.f1 - 8*m1.f1 + 8*p1.f1 - p2.f1);
    r.f2 = c1 * (m2.f2 - 8*m1.f2 + 8*p1.f2 - p2.f2);
    r.f3 = c1 * (m2.f3 - 8*m1.f3 + 8*p1.f3 - p2.f3);
    r.j1 = c1 * (m2.j1 - 8*m1.j1 + 8*p1.j1 - p2.j1);
    r.j2 = c1 * (m2.j2 - 8*m1.j2 + 8*p1.j2 - p2.j2);
    r.j3 = c1 * (m2.j3 - 8*m1.j3 + 8*p1.j3 - p2.j3);
    r.p  = c1 * (m2.p  - 8*m1.p  + 8*p1.p  - p2.p);
    return r;
}

/* Consistent Laplacian: exact discrete gradient of E₂ energy.
 *
 * The E₂ energy uses 4th-order central diffs D_d for first derivatives.
 * The exact discrete gradient of E₂ w.r.t. q(x) is Σ_d D_d(D_d q)(x),
 * which gives a 9-point stencil per direction:
 *
 *   weights at offsets -4..+4: {1, -16, 64, 16, -130, 16, 64, -16, 1} / (144 h²)
 *
 * This ensures F = -∇E exactly on the lattice (energy-force consistency).
 */
static inline Multivector consistent_laplacian(const SphericalGrid *g, int i, int j, int k)
{
    double c = 1.0 / (144.0 * g->h * g->h);
    static const double w[9] = {1, -16, 64, 16, -130, 16, 64, -16, 1};

    double rs = 0, rf1 = 0, rf2 = 0, rf3 = 0;
    double rj1 = 0, rj2 = 0, rj3 = 0, rp = 0;

    /* X direction */
    for (int m = -4; m <= 4; m++) {
        Multivector p = sg_get(g, i+m, j, k);
        double wm = w[m+4];
        rs += wm*p.s;  rf1 += wm*p.f1; rf2 += wm*p.f2; rf3 += wm*p.f3;
        rj1 += wm*p.j1; rj2 += wm*p.j2; rj3 += wm*p.j3; rp += wm*p.p;
    }
    /* Y direction */
    for (int m = -4; m <= 4; m++) {
        Multivector p = sg_get(g, i, j+m, k);
        double wm = w[m+4];
        rs += wm*p.s;  rf1 += wm*p.f1; rf2 += wm*p.f2; rf3 += wm*p.f3;
        rj1 += wm*p.j1; rj2 += wm*p.j2; rj3 += wm*p.j3; rp += wm*p.p;
    }
    /* Z direction */
    for (int m = -4; m <= 4; m++) {
        Multivector p = sg_get(g, i, j, k+m);
        double wm = w[m+4];
        rs += wm*p.s;  rf1 += wm*p.f1; rf2 += wm*p.f2; rf3 += wm*p.f3;
        rj1 += wm*p.j1; rj2 += wm*p.j2; rj3 += wm*p.j3; rp += wm*p.p;
    }

    Multivector r;
    r.s = c*rs; r.f1 = c*rf1; r.f2 = c*rf2; r.f3 = c*rf3;
    r.j1 = c*rj1; r.j2 = c*rj2; r.j3 = c*rj3; r.p = c*rp;
    return r;
}

/* ========== Left-invariant current R_i = q̃ ∂_i q ========== */

/* Quaternion-only operations (bulk sector, 4 components) */
typedef struct { double s, f1, f2, f3; } Quat4;

static inline Quat4 q4_from_mv(Multivector m) {
    return (Quat4){m.s, m.f1, m.f2, m.f3};
}

static inline Quat4 q4_rev(Quat4 a) {
    return (Quat4){a.s, -a.f1, -a.f2, -a.f3};
}

static inline Quat4 q4_mul(Quat4 a, Quat4 b) {
    return (Quat4){
        a.s*b.s  - a.f1*b.f1 - a.f2*b.f2 - a.f3*b.f3,
        a.s*b.f1 + a.f1*b.s  - a.f2*b.f3 + a.f3*b.f2,
        a.s*b.f2 + a.f1*b.f3 + a.f2*b.s  - a.f3*b.f1,
        a.s*b.f3 - a.f1*b.f2 + a.f2*b.f1 + a.f3*b.s
    };
}

static inline double q4_norm2(Quat4 a) {
    return a.s*a.s + a.f1*a.f1 + a.f2*a.f2 + a.f3*a.f3;
}

static inline Quat4 q4_add(Quat4 a, Quat4 b) {
    return (Quat4){a.s+b.s, a.f1+b.f1, a.f2+b.f2, a.f3+b.f3};
}

static inline Quat4 q4_sub(Quat4 a, Quat4 b) {
    return (Quat4){a.s-b.s, a.f1-b.f1, a.f2-b.f2, a.f3-b.f3};
}

static inline Quat4 q4_scale(double c, Quat4 a) {
    return (Quat4){c*a.s, c*a.f1, c*a.f2, c*a.f3};
}

static inline Quat4 q4_zero(void) {
    return (Quat4){0, 0, 0, 0};
}

/* Clifford scalar product for quaternions: <AB~>_0 = a.s*b.s - a.f1*b.f1 - ... */
static inline double q4_scalar_prod(Quat4 a, Quat4 b) {
    return a.s*b.s - a.f1*b.f1 - a.f2*b.f2 - a.f3*b.f3;
}

/* Left multiplication by basis quaternions:
 * eps[0] = 1, eps[1] = e23, eps[2] = e31, eps[3] = e12 */
static inline Quat4 q4_left_basis(int a, Quat4 q) {
    switch(a) {
    case 0: return q;
    case 1: return (Quat4){-q.f1,  q.s,   q.f3, -q.f2}; /* e23 * q */
    case 2: return (Quat4){-q.f2, -q.f3,  q.s,   q.f1}; /* e31 * q */
    case 3: return (Quat4){-q.f3,  q.f2, -q.f1,  q.s};  /* e12 * q */
    default: return q4_zero();
    }
}

/* ========== Baryon density ========== */

/* B⁰ = -(1/2π²) ε^{ijk} Tr(R_i R_j R_k) / (2|q|⁶)
 * where R_i = q̃ ∂_i q (left current)
 *
 * Expanded: B⁰ = -(1/2π²|q|⁶) ε^{ijk} ⟨R_i [R_j, R_k]⟩₀
 *
 * For unit quaternion (|q|=1), this gives ε^{ijk} ⟨(q̃∂_i q)[(q̃∂_j q),(q̃∂_k q)]⟩₀
 * which is 12× the baryon density.
 *
 * The standard formula for B⁰ on S³ (unit quaternion) is:
 *   B⁰ = -(1/2π²) ε_{abc} ε^{ijk} (q̃∂_i q)^a (q̃∂_j q)^b (q̃∂_k q)^c
 * where (q̃∂_i q)^a are the imaginary quaternion components.
 */
double hf_baryon_density(const SphericalGrid *g, int i, int j, int k)
{
    Multivector psi = sg_get(g, i, j, k);
    Quat4 q = q4_from_mv(psi);
    double q2 = q4_norm2(q);
    if (q2 < 1e-20) return 0.0;

    Quat4 qr = q4_rev(q);  /* q̃ */

    /* Compute R_i = q̃ ∂_i q (only imaginary parts matter for B⁰) */
    Multivector dx = deriv_x(g, i, j, k);
    Multivector dy = deriv_y(g, i, j, k);
    Multivector dz = deriv_z(g, i, j, k);

    Quat4 Rx = q4_mul(qr, q4_from_mv(dx));
    Quat4 Ry = q4_mul(qr, q4_from_mv(dy));
    Quat4 Rz = q4_mul(qr, q4_from_mv(dz));

    /* B⁰ = -(1/2π²) det[Rx_im, Ry_im, Rz_im] / |q|⁶
     * det = Rx1*(Ry2*Rz3 - Ry3*Rz2) - Rx2*(Ry1*Rz3 - Ry3*Rz1) + Rx3*(Ry1*Rz2 - Ry2*Rz1) */
    double det = Rx.f1 * (Ry.f2 * Rz.f3 - Ry.f3 * Rz.f2)
               - Rx.f2 * (Ry.f1 * Rz.f3 - Ry.f3 * Rz.f1)
               + Rx.f3 * (Ry.f1 * Rz.f2 - Ry.f2 * Rz.f1);

    double q6 = q2 * q2 * q2;
    return -det / (2.0 * M_PI * M_PI * q6);
}

double hf_baryon_charge(const SphericalGrid *g)
{
    double h3 = g->h * g->h * g->h;
    double Q = 0.0;

    #pragma omp parallel for reduction(+:Q) schedule(static)
    for (int n = 0; n < g->n_active; n++) {
        int ix = g->active[n];
        int i, j, k;
        sg_unflatten(g->N, ix, &i, &j, &k);

        /* Skip boundary cells (stencil needs ±2 neighbors) */
        double r = sg_radius(g, i, j, k);
        if (r > g->R - 4.5*g->h) continue;

        Q += hf_baryon_density(g, i, j, k) * h3;
    }
    return Q;
}

/* ========== Energy computation ========== */

/* E₄ (Skyrme) energy density at cell (i,j,k) */
static double skyrme_energy_density(const SphericalGrid *g, const HopfionParams *p,
                                     int i, int j, int k)
{
    Multivector psi = sg_get(g, i, j, k);
    Quat4 q = q4_from_mv(psi);
    Quat4 qr = q4_rev(q);

    Multivector dx = deriv_x(g, i, j, k);
    Multivector dy = deriv_y(g, i, j, k);
    Multivector dz = deriv_z(g, i, j, k);

    /* Left currents R_i = q̃ ∂_i q */
    Quat4 Rx = q4_mul(qr, q4_from_mv(dx));
    Quat4 Ry = q4_mul(qr, q4_from_mv(dy));
    Quat4 Rz = q4_mul(qr, q4_from_mv(dz));

    /* Commutators [R_i, R_j]: only 3 independent pairs (xy, xz, yz) */
    /* [A,B]_quaternion = AB - BA. For pure imaginary quaternions:
     * [A,B] has scalar part = 0 and vector part = 2(A×B)
     * But R_i has a scalar part too, so use full formula */
    Quat4 Cxy = {
        0,  /* scalar: A1B1+... cancels */
        2*(Rx.f2*Ry.f3 - Rx.f3*Ry.f2),
        2*(Rx.f3*Ry.f1 - Rx.f1*Ry.f3),
        2*(Rx.f1*Ry.f2 - Rx.f2*Ry.f1)
    };
    Quat4 Cxz = {
        0,
        2*(Rx.f2*Rz.f3 - Rx.f3*Rz.f2),
        2*(Rx.f3*Rz.f1 - Rx.f1*Rz.f3),
        2*(Rx.f1*Rz.f2 - Rx.f2*Rz.f1)
    };
    Quat4 Cyz = {
        0,
        2*(Ry.f2*Rz.f3 - Ry.f3*Rz.f2),
        2*(Ry.f3*Rz.f1 - Ry.f1*Rz.f3),
        2*(Ry.f1*Rz.f2 - Ry.f2*Rz.f1)
    };

    /* E4 = (1/4e²) Σ |[R_i,R_j]|² */
    double norm_xy = q4_norm2(Cxy);
    double norm_xz = q4_norm2(Cxz);
    double norm_yz = q4_norm2(Cyz);

    return (norm_xy + norm_xz + norm_yz) / (4.0 * p->e_skyrme * p->e_skyrme);
}

HopfionEnergy hf_energy(const SphericalGrid *g, const HopfionParams *p)
{
    double h3 = g->h * g->h * g->h;
    double E2 = 0, E4 = 0, E6 = 0, EV = 0, Epi = 0, ED = 0;

    #pragma omp parallel for reduction(+:E2,E4,E6,EV,Epi,ED) schedule(static)
    for (int n = 0; n < g->n_active; n++) {
        int ix = g->active[n];
        int i, j, k;
        sg_unflatten(g->N, ix, &i, &j, &k);

        double r = sg_radius(g, i, j, k);
        if (r > g->R - 4.5*g->h) continue;  /* skip near-boundary */

        Multivector psi = g->psi[ix];

        /* E₂: gradient energy = (1/2)|∂q|² */
        Multivector dx = deriv_x(g, i, j, k);
        Multivector dy = deriv_y(g, i, j, k);
        Multivector dz = deriv_z(g, i, j, k);

        double grad2 = mv_dot(dx, dx) + mv_dot(dy, dy) + mv_dot(dz, dz);
        E2 += 0.5 * grad2 * h3;

        /* E₄: Skyrme */
        E4 += skyrme_energy_density(g, p, i, j, k) * h3;

        /* E₆: Sextic (BPS) = λ₆ (B⁰)² */
        if (p->lambda6 > 0) {
            double B0 = hf_baryon_density(g, i, j, k);
            E6 += p->lambda6 * B0 * B0 * h3;
        }

        /* V: Mexican hat */
        double q2 = mv_bulk_norm2(psi);
        double dq = q2 - p->rho0 * p->rho0;
        EV += (p->lambda / 4.0) * dq * dq * h3;

        /* V_π: Pion mass */
        if (p->m_pi_sq > 0) {
            double qnorm = sqrt(q2);
            if (qnorm > 1e-15) {
                Epi += p->m_pi_sq * p->rho0 * p->rho0 * (1.0 - psi.s / qnorm) * h3;
            }
        }

        /* V_D: Degenerate mass */
        double w2 = mv_weight_norm2(psi);
        ED += 0.5 * p->mu * p->mu * w2 * h3;
    }

    HopfionEnergy en;
    en.E2 = E2;
    en.E4 = E4;
    en.E6 = E6;
    en.EV = EV;
    en.Epi = Epi;
    en.ED = ED;
    en.Etotal = E2 + E4 + E6 + EV + Epi + ED;
    return en;
}

double hf_kinetic_energy(const SphericalGrid *g, const HopfionParams *p)
{
    double h3 = g->h * g->h * g->h;
    double Ek = 0;
    double c2_inv = 1.0 / (p->c * p->c);

    #pragma omp parallel for reduction(+:Ek) schedule(static)
    for (int n = 0; n < g->n_active; n++) {
        int ix = g->active[n];
        Ek += 0.5 * c2_inv * mv_dot(g->vel[ix], g->vel[ix]) * h3;
    }
    return Ek;
}

/* ========== Skyrme force (analytical) ==========
 *
 * Ported from field.c. Three-pass algorithm:
 *
 * Pass 1: Precompute at each cell:
 *   A_d = q̃ ∂_d q  (right-currents)
 *   C_{d1,d2} = [A_d1, A_d2]  (commutators)
 *   G_d = Σ_{d'≠d} [A_{d'}, C_{d,d'}]  (Skyrme tensor)
 *
 * Pass 2: Precompute pi_{a,d}(x) = <q̃(x) ε_a, G_d(x)>₀
 *   where ε_a is the a-th basis quaternion.
 *
 * Pass 3: Compute force per component a:
 *   F4_a = (1/2e²) Σ_d { σ_a <ε_a·∂_d q, G_d>₀ - D_d(pi_{a,d}) }
 *   where σ_a = reversion sign (+1, -1, -1, -1).
 */

/* Precomputed Skyrme data per grid point */
typedef struct {
    Quat4 A[3];     /* right-currents */
    Quat4 G[3];     /* Skyrme G-tensor */
} SkyrmePre;

/* Quaternion-only derivative (bulk sector) by direction index */
static inline Quat4 q4_deriv(const SphericalGrid *g, int i, int j, int k, int dir)
{
    double c1 = 1.0 / (12.0 * g->h);
    Multivector m2, m1, p1, p2;

    switch (dir) {
    case 0:
        m2 = sg_get(g, i-2, j, k); m1 = sg_get(g, i-1, j, k);
        p1 = sg_get(g, i+1, j, k); p2 = sg_get(g, i+2, j, k);
        break;
    case 1:
        m2 = sg_get(g, i, j-2, k); m1 = sg_get(g, i, j-1, k);
        p1 = sg_get(g, i, j+1, k); p2 = sg_get(g, i, j+2, k);
        break;
    default:
        m2 = sg_get(g, i, j, k-2); m1 = sg_get(g, i, j, k-1);
        p1 = sg_get(g, i, j, k+1); p2 = sg_get(g, i, j, k+2);
        break;
    }

    Quat4 r;
    r.s  = c1 * (m2.s  - 8*m1.s  + 8*p1.s  - p2.s);
    r.f1 = c1 * (m2.f1 - 8*m1.f1 + 8*p1.f1 - p2.f1);
    r.f2 = c1 * (m2.f2 - 8*m1.f2 + 8*p1.f2 - p2.f2);
    r.f3 = c1 * (m2.f3 - 8*m1.f3 + 8*p1.f3 - p2.f3);
    return r;
}

/* Safe access to pi array: returns 0 for out-of-bounds or inactive cells
 * (vacuum has zero currents → pi = 0) */
static inline double pi_get(const double *pi, const SphericalGrid *g,
                             int i, int j, int k, int a, int d)
{
    if (i < 0 || i >= g->N || j < 0 || j >= g->N || k < 0 || k >= g->N)
        return 0.0;
    int ix = sg_idx(g->N, i, j, k);
    if (!g->mask[ix])
        return 0.0;
    return pi[12*ix + 3*a + d];
}

static void precompute_skyrme(const SphericalGrid *g, SkyrmePre *pre)
{
    int N = g->N;

    #pragma omp parallel for schedule(static)
    for (int n = 0; n < g->n_active; n++) {
        int ix = g->active[n];
        int i, j, k;
        sg_unflatten(N, ix, &i, &j, &k);

        double r = sg_radius(g, i, j, k);
        if (r > g->R - 4.5*g->h) {
            for (int d = 0; d < 3; d++) {
                pre[ix].A[d] = q4_zero();
                pre[ix].G[d] = q4_zero();
            }
            continue;
        }

        Quat4 q = q4_from_mv(g->psi[ix]);
        Quat4 qr = q4_rev(q);

        /* Right-currents A_d = q̃ ∂_d q */
        for (int d = 0; d < 3; d++) {
            Quat4 dq = q4_deriv(g, i, j, k, d);
            pre[ix].A[d] = q4_mul(qr, dq);
        }

        /* Commutators C_{d1,d2} = [A_d1, A_d2]
         * C[0]=C_{01}, C[1]=C_{02}, C[2]=C_{12} */
        Quat4 C[3];
        C[0] = q4_sub(q4_mul(pre[ix].A[0], pre[ix].A[1]),
                       q4_mul(pre[ix].A[1], pre[ix].A[0]));
        C[1] = q4_sub(q4_mul(pre[ix].A[0], pre[ix].A[2]),
                       q4_mul(pre[ix].A[2], pre[ix].A[0]));
        C[2] = q4_sub(q4_mul(pre[ix].A[1], pre[ix].A[2]),
                       q4_mul(pre[ix].A[2], pre[ix].A[1]));

        /* G_d = Σ_{d'≠d} [A_{d'}, C_{d,d'}] */
        for (int d = 0; d < 3; d++) {
            pre[ix].G[d] = q4_zero();
            for (int dp = 0; dp < 3; dp++) {
                if (dp == d) continue;
                /* Get C_{d,dp} (antisymmetric) */
                Quat4 Cdd;
                if (d < dp) {
                    int ci = (d==0 && dp==1) ? 0 : (d==0 && dp==2) ? 1 : 2;
                    Cdd = C[ci];
                } else {
                    int ci = (dp==0 && d==1) ? 0 : (dp==0 && d==2) ? 1 : 2;
                    Cdd = q4_scale(-1.0, C[ci]);
                }
                /* [A_{d'}, C_{d,d'}] */
                Quat4 term = q4_sub(q4_mul(pre[ix].A[dp], Cdd),
                                    q4_mul(Cdd, pre[ix].A[dp]));
                pre[ix].G[d] = q4_add(pre[ix].G[d], term);
            }
        }
    }
}

/* ========== L₆ Force (numerical finite differences) ==========
 *
 * F₆^a(x) = -λ₆ d/dq_a [ Σ_y (B⁰(y))² h³ ]
 *
 * Perturbing q_a(x) by ε only affects B⁰(y) at cells y within ±2 of x
 * (because B⁰ uses ±2 stencil derivatives). So the force at x requires
 * recomputing B⁰ at the local 5³ neighborhood.
 */
static void compute_l6_force(SphericalGrid *g, const HopfionParams *p)
{
    if (p->lambda6 <= 0) return;

    int N = g->N;
    double h3 = g->h * g->h * g->h;
    double eps = 1e-5 * g->h;  /* perturbation size, scaled by h */

    /* NOTE: This loop MUST be serial. Each iteration perturbs g->psi[ix]
     * in place and reads neighbors' fields via hf_baryon_density. Parallel
     * execution causes data races (thread A's perturbation corrupts
     * thread B's baryon density computation at overlapping stencil cells). */
    for (int n = 0; n < g->n_active; n++) {
        int ix = g->active[n];
        int ci, cj, ck;
        sg_unflatten(N, ix, &ci, &cj, &ck);

        double r = sg_radius(g, ci, cj, ck);
        if (r > g->R - 4.5*g->h) continue;

        /* For each bulk component a ∈ {s, f1, f2, f3} */
        double f6[4] = {0, 0, 0, 0};
        double *comp_ptr[4];
        Multivector *psi_x = &g->psi[ix];
        comp_ptr[0] = &psi_x->s;
        comp_ptr[1] = &psi_x->f1;
        comp_ptr[2] = &psi_x->f2;
        comp_ptr[3] = &psi_x->f3;

        for (int a = 0; a < 4; a++) {
            double orig = *comp_ptr[a];

            /* Compute Σ (B⁰)² at affected neighbors for +ε */
            *comp_ptr[a] = orig + eps;
            double sum_plus = 0;
            for (int di = -2; di <= 2; di++)
            for (int dj = -2; dj <= 2; dj++)
            for (int dk = -2; dk <= 2; dk++) {
                int ni = ci+di, nj = cj+dj, nk = ck+dk;
                if (!sg_is_active(g, ni, nj, nk)) continue;
                double rn = sg_radius(g, ni, nj, nk);
                if (rn > g->R - 4.5*g->h) continue;
                double B0 = hf_baryon_density(g, ni, nj, nk);
                sum_plus += B0 * B0;
            }

            /* Compute Σ (B⁰)² at affected neighbors for -ε */
            *comp_ptr[a] = orig - eps;
            double sum_minus = 0;
            for (int di = -2; di <= 2; di++)
            for (int dj = -2; dj <= 2; dj++)
            for (int dk = -2; dk <= 2; dk++) {
                int ni = ci+di, nj = cj+dj, nk = ck+dk;
                if (!sg_is_active(g, ni, nj, nk)) continue;
                double rn = sg_radius(g, ni, nj, nk);
                if (rn > g->R - 4.5*g->h) continue;
                double B0 = hf_baryon_density(g, ni, nj, nk);
                sum_minus += B0 * B0;
            }

            /* Restore and compute force */
            *comp_ptr[a] = orig;
            f6[a] = -p->lambda6 * (sum_plus - sum_minus) * h3 / (2.0 * eps);
        }

        /* Add L₆ force (only bulk components) */
        g->force[ix].s  += f6[0];
        g->force[ix].f1 += f6[1];
        g->force[ix].f2 += f6[2];
        g->force[ix].f3 += f6[3];
    }
}

/* ========== Force computation ========== */

void hf_force(SphericalGrid *g, const HopfionParams *p)
{
    int N = g->N;
    int N3 = N * N * N;
    double inv_2e2 = 1.0 / (2.0 * p->e_skyrme * p->e_skyrme);
    double rho0_sq = p->rho0 * p->rho0;

    /* --- Pass 1: Precompute Skyrme data (A_d, G_d) --- */
    SkyrmePre *pre = (SkyrmePre *)calloc((size_t)N3, sizeof(SkyrmePre));
    if (!pre) { fprintf(stderr, "hf_force: malloc SkyrmePre failed\n"); exit(1); }
    precompute_skyrme(g, pre);

    /* --- Pass 2: Precompute pi_{a,d}(x) = <q̃(x)·ε_a, G_d(x)>₀ --- */
    double *pi = (double *)calloc((size_t)N3 * 12, sizeof(double));
    if (!pi) { fprintf(stderr, "hf_force: malloc pi failed\n"); exit(1); }

    #pragma omp parallel for schedule(static)
    for (int n = 0; n < g->n_active; n++) {
        int ix = g->active[n];
        int i, j, k;
        sg_unflatten(N, ix, &i, &j, &k);

        double r = sg_radius(g, i, j, k);
        if (r > g->R - 4.5*g->h) continue;

        Quat4 qr = q4_rev(q4_from_mv(g->psi[ix]));
        for (int d = 0; d < 3; d++) {
            for (int a = 0; a < 4; a++) {
                /* q̃ · ε_a (right multiplication by basis element) */
                Quat4 qr_ea;
                switch (a) {
                case 0: qr_ea = qr; break;
                case 1: qr_ea = q4_mul(qr, (Quat4){0,1,0,0}); break;
                case 2: qr_ea = q4_mul(qr, (Quat4){0,0,1,0}); break;
                case 3: qr_ea = q4_mul(qr, (Quat4){0,0,0,1}); break;
                default: qr_ea = q4_zero();
                }
                /* Scalar part of (qr_ea · G_d) */
                Quat4 prod = q4_mul(qr_ea, pre[ix].G[d]);
                pi[12*ix + 3*a + d] = prod.s;
            }
        }
    }

    /* --- Pass 3: Compute all forces --- */
    static const double sigma[4] = {1.0, -1.0, -1.0, -1.0};

    #pragma omp parallel for schedule(static)
    for (int n = 0; n < g->n_active; n++) {
        int ix = g->active[n];
        int i, j, k;
        sg_unflatten(N, ix, &i, &j, &k);

        double r = sg_radius(g, i, j, k);
        if (r > g->R - 4.5*g->h) {
            g->force[ix] = mv_zero();
            continue;
        }

        Multivector psi = g->psi[ix];
        Multivector F = mv_zero();

        /* --- F_E2: consistent Laplacian --- */
        Multivector lap = consistent_laplacian(g, i, j, k);
        F = mv_add(F, lap);

        /* --- F_V: Mexican hat force --- */
        double q2 = mv_bulk_norm2(psi);
        double dq = q2 - rho0_sq;
        if (p->lambda > 0) {
            double coeff = -p->lambda * dq;
            Multivector fv = mv_zero();
            fv.s  = coeff * psi.s;
            fv.f1 = coeff * psi.f1;
            fv.f2 = coeff * psi.f2;
            fv.f3 = coeff * psi.f3;
            F = mv_add(F, fv);
        }

        /* --- F_π: pion mass force --- */
        if (p->m_pi_sq > 0 && q2 > 1e-20) {
            double inv_norm3 = 1.0 / (q2 * sqrt(q2));
            double coeff = p->m_pi_sq * rho0_sq;
            /* F_π = -∂V_π/∂q, same formula as field.c */
            F.s  += coeff * (q2 - psi.s * psi.s) * inv_norm3;
            F.f1 -= coeff * psi.s * psi.f1 * inv_norm3;
            F.f2 -= coeff * psi.s * psi.f2 * inv_norm3;
            F.f3 -= coeff * psi.s * psi.f3 * inv_norm3;
        }

        /* --- F_D: degenerate mass force --- */
        if (p->mu > 0) {
            double mu2 = p->mu * p->mu;
            F.j1 -= mu2 * psi.j1;
            F.j2 -= mu2 * psi.j2;
            F.j3 -= mu2 * psi.j3;
            F.p  -= mu2 * psi.p;
        }

        /* --- F_E4: Analytical Skyrme force (from field.c) --- */
        /* Term 1: σ_a <ε_a · ∂_d q, G_d>₀ summed over d */
        double f4[4] = {0, 0, 0, 0};
        for (int d = 0; d < 3; d++) {
            Quat4 dq = q4_deriv(g, i, j, k, d);
            for (int a = 0; a < 4; a++) {
                Quat4 ea_dq = q4_left_basis(a, dq);
                Quat4 prod = q4_mul(ea_dq, pre[ix].G[d]);
                f4[a] += sigma[a] * prod.s;
            }
        }

        /* Term 2: -D_d(pi_{a,d}) using 4th-order central difference */
        double inv12h = 1.0 / (12.0 * g->h);
        for (int a = 0; a < 4; a++) {
            for (int d = 0; d < 3; d++) {
                double pm2, pm1, pp1, pp2;
                switch (d) {
                case 0:
                    pm2 = pi_get(pi, g, i-2, j, k, a, d);
                    pm1 = pi_get(pi, g, i-1, j, k, a, d);
                    pp1 = pi_get(pi, g, i+1, j, k, a, d);
                    pp2 = pi_get(pi, g, i+2, j, k, a, d);
                    break;
                case 1:
                    pm2 = pi_get(pi, g, i, j-2, k, a, d);
                    pm1 = pi_get(pi, g, i, j-1, k, a, d);
                    pp1 = pi_get(pi, g, i, j+1, k, a, d);
                    pp2 = pi_get(pi, g, i, j+2, k, a, d);
                    break;
                default:
                    pm2 = pi_get(pi, g, i, j, k-2, a, d);
                    pm1 = pi_get(pi, g, i, j, k-1, a, d);
                    pp1 = pi_get(pi, g, i, j, k+1, a, d);
                    pp2 = pi_get(pi, g, i, j, k+2, a, d);
                    break;
                }
                double dpi = (-pp2 + 8*pp1 - 8*pm1 + pm2) * inv12h;
                f4[a] -= dpi;
            }
        }

        /* Apply Skyrme force with prefactor 1/(2e²) */
        F.s  += inv_2e2 * f4[0];
        F.f1 += inv_2e2 * f4[1];
        F.f2 += inv_2e2 * f4[2];
        F.f3 += inv_2e2 * f4[3];

        g->force[ix] = F;
    }

    free(pi);
    free(pre);

    /* --- L₆ force (numerical, separate pass) --- */
    compute_l6_force(g, p);
}

/* ========== Sigma-model projection ========== */

void hf_sigma_project(SphericalGrid *g, double rho0)
{
    #pragma omp parallel for schedule(static)
    for (int n = 0; n < g->n_active; n++) {
        int ix = g->active[n];
        Multivector *psi = &g->psi[ix];
        double q2 = mv_bulk_norm2(*psi);
        if (q2 > 1e-20) {
            double scale = rho0 / sqrt(q2);
            psi->s  *= scale;
            psi->f1 *= scale;
            psi->f2 *= scale;
            psi->f3 *= scale;
        }
    }
}

/* ========== Effective Metric ========== */

/*
 * BLV effective metric for Skyrme theory:
 *
 * The L₂ contribution to the effective metric is trivially η^{μν} (flat).
 * The L₄ (Skyrme) contribution depends on the background currents A_d = q̃∂_dq.
 *
 * For a general field configuration, the effective spatial metric at cell x is:
 *   g_4^{ij}(x) = (c₄/|q|⁴) [ Σ_{k≠i,k≠j} |A_k|² δ^{ij} - A_i·A_j (for i≠j) ]
 *
 * For the temporal component:
 *   g_4^{00}(x) = (c₄/|q|⁴) Σ_k |A_k|²
 *
 * where c₄ = 2ρ₀²/e² is the Skyrme coupling.
 *
 * We spherically average to get radial profiles g_00(r) and g_rr(r).
 * g_rr(r) is the radial-radial component of the effective metric.
 *
 * For a hedgehog, this gives:
 *   P(r) = total radial stiffness = 2r² + 4c₄ sin²f
 *   m(r) = inertia = r² + 2c₄ sin²f
 * and the effective line element ds² ~ -m(r)dt² + P(r)dr² + ...
 *
 * The ratio P/m = 2 for sigma-model (constant ρ₀) → no time dilation.
 */

EffectiveMetric *hf_effective_metric(const SphericalGrid *g, const HopfionParams *p,
                                      int n_bins)
{
    EffectiveMetric *em = (EffectiveMetric *)malloc(sizeof(EffectiveMetric));
    if (!em) return NULL;

    em->n_bins = n_bins;
    em->dr = g->R / n_bins;
    em->r   = (double *)calloc(n_bins, sizeof(double));
    em->g00 = (double *)calloc(n_bins, sizeof(double));
    em->grr = (double *)calloc(n_bins, sizeof(double));
    em->P   = (double *)calloc(n_bins, sizeof(double));
    em->m   = (double *)calloc(n_bins, sizeof(double));

    if (!em->r || !em->g00 || !em->grr || !em->P || !em->m) {
        hf_metric_free(em);
        return NULL;
    }

    for (int b = 0; b < n_bins; b++)
        em->r[b] = (b + 0.5) * em->dr;

    /* Accumulate contributions per radial bin */
    double *count = (double *)calloc(n_bins, sizeof(double));
    /* Per-bin accumulators for metric components */
    double *sum_g00_l2 = (double *)calloc(n_bins, sizeof(double));
    double *sum_g00_l4 = (double *)calloc(n_bins, sizeof(double));
    double *sum_grr_l2 = (double *)calloc(n_bins, sizeof(double));
    double *sum_grr_l4 = (double *)calloc(n_bins, sizeof(double));
    /* L₆ accumulator: β₆ = (2λ₆/π³) sin⁴f / r², adds equally to P and m */
    double *sum_beta6 = (double *)calloc(n_bins, sizeof(double));

    if (!count || !sum_g00_l2 || !sum_g00_l4 || !sum_grr_l2 || !sum_grr_l4
        || !sum_beta6) {
        free(count); free(sum_g00_l2); free(sum_g00_l4);
        free(sum_grr_l2); free(sum_grr_l4); free(sum_beta6);
        hf_metric_free(em);
        return NULL;
    }

    int N = g->N;
    double c4 = 2.0 * p->rho0 * p->rho0 / (p->e_skyrme * p->e_skyrme);

    for (int n = 0; n < g->n_active; n++) {
        int ix = g->active[n];
        int i, j, k;
        sg_unflatten(N, ix, &i, &j, &k);

        double x, y, z;
        sg_pos(g, i, j, k, &x, &y, &z);
        double r = sqrt(x*x + y*y + z*z);

        if (r > g->R - 4.5*g->h || r < 1e-10) continue;

        int bin = (int)(r / em->dr);
        if (bin < 0 || bin >= n_bins) continue;

        /* Compute right-currents A_d = q̃ ∂_d q */
        Quat4 q = q4_from_mv(g->psi[ix]);
        double q2 = q4_norm2(q);
        if (q2 < 1e-20) continue;
        double q4_inv = 1.0 / (q2 * q2);

        Quat4 qr = q4_rev(q);
        Quat4 A[3];
        for (int d = 0; d < 3; d++) {
            Quat4 dq = q4_deriv(g, i, j, k, d);
            A[d] = q4_mul(qr, dq);
        }

        /* |A_d|² for each direction */
        double A2[3];
        for (int d = 0; d < 3; d++)
            A2[d] = q4_norm2(A[d]);
        double A2_total = A2[0] + A2[1] + A2[2];

        /* L₂ contribution: g^{μν} = η^{μν} → g00=1, grr=1 (flat) */
        sum_g00_l2[bin] += 1.0;
        sum_grr_l2[bin] += 1.0;

        /* L₄ contribution to g^{00}: proportional to Σ_k |A_k|² */
        sum_g00_l4[bin] += c4 * q4_inv * A2_total;

        /* L₄ contribution to g^{rr}: need radial-radial component.
         * g_4^{rr} = (c₄/|q|⁴) [ Σ_{k≠r} |A_k|² ]
         *          = (c₄/|q|⁴) [ A2_total - |A_r|² ]
         * where A_r = (x/r)A_x + (y/r)A_y + (z/r)A_z
         * and |A_r|² = Σ_{ij} (x_i x_j / r²) <A_i, A_j>₀
         *
         * Actually, the full Skyrme effective metric in direction (i,j) is:
         * g_4^{ij} = (c₄/|q|⁴) [ (A2_total - A_i·A_i) δ^{ij} - A_i·A_j (i≠j) ]
         *          ... this is wrong. Let me be more careful.
         *
         * The variation of E₄ w.r.t. metric (Babichev-Langlois-Vernieri):
         * g_4^{ij} = (c₄/|q|⁴) Σ_{k} (|A_k|² δ^{ij} - 2 <A_i·A_k> <A_j·A_k> / ...  )
         *
         * For the Skyrme term E₄ = (1/4e²) Σ_{ij} |[A_i,A_j]|², the effective
         * metric tensor for small perturbations δφ propagating on the background is:
         *
         * From the Sturm-Liouville analysis (normal_modes.c), for a hedgehog:
         *   P(r) = 2r² + 4c₄ sin²f  (radial stiffness coefficient)
         *   m(r) = r² + 2c₄ sin²f   (inertia coefficient)
         *
         * For a general 3D field, we compute P and m by projecting:
         *   P(r) = spatial stiffness along radial direction
         *   m(r) = temporal inertia
         *
         * The L₂ contribution to both is r² (from the r² in spherical Laplacian).
         * The L₄ contribution is c₄ × (angular current content).
         *
         * For a general field: Σ_k |A_k|² is the total current squared.
         * The angular current content (perpendicular to r̂) contributes to P and m.
         *
         * Decompose: |A_r|² = radial current, |A_⊥|² = A2_total - |A_r|²
         */
        double r_hat[3] = {x/r, y/r, z/r};

        /* A_r (scalar-valued projection of each A component onto radial) */
        /* Actually A_d are quaternions; the "radial current" in the
         * Sturm-Liouville sense is the current along r̂:
         * A_r = Σ_d r̂_d A_d (quaternion sum) */
        Quat4 A_radial = q4_zero();
        for (int d = 0; d < 3; d++)
            A_radial = q4_add(A_radial, q4_scale(r_hat[d], A[d]));
        double A_r2 = q4_norm2(A_radial);
        double A_perp2 = A2_total - A_r2;

        /* For hedgehog: A_r has only imaginary part ~ f'(r) r̂·σ,
         * |A_r|² = f'² → after r² factor, gives r²f'² contribution to E₂.
         * A_⊥: gives sin²f/r² type terms.
         *
         * Sturm-Liouville coefficients (generalizing the 1D result):
         * P(r) = 2×(L₂ stiffness) + L₄ contribution
         * For L₂: stiffness in radial direction = 1 per cell (already in Laplacian)
         * For L₄: P gets 4c₄ × (|A_⊥|² / |q|⁴) from Skyrme perpendicular currents
         *         m gets 2c₄ × (|A_⊥|² / |q|⁴)
         * (These correspond to P(r)=2r²+4c₄sin²f and m(r)=r²+2c₄sin²f
         *  for the hedgehog where |A_⊥|²/|q|⁴ = sin²f/r²)
         */
        double perp_contrib = c4 * q4_inv * A_perp2;

        /* L₆ contribution: β₆ = (2λ₆/π³) sin⁴f / r²
         * sin²f = 4 sin²(f/2) cos²(f/2) = 4|f⃗|²s²/|q|⁴
         * where s = q.s (scalar part), |f⃗|² = q.f1²+q.f2²+q.f3² */
        if (p->lambda6 > 0) {
            double fvec2 = q.f1*q.f1 + q.f2*q.f2 + q.f3*q.f3;
            double sin2f = 4.0 * fvec2 * q.s * q.s / (q2 * q2);
            double sin4f = sin2f * sin2f;
            double beta6 = (2.0 * p->lambda6 / (M_PI * M_PI * M_PI)) * sin4f / (r * r);
            sum_beta6[bin] += beta6;
        }

        /* Accumulate raw metric components */
        sum_grr_l4[bin] += perp_contrib;
        count[bin] += 1.0;
    }

    /* Average and build metric profiles */
    for (int b = 0; b < n_bins; b++) {
        if (count[b] < 1.0) {
            em->g00[b] = 1.0;
            em->grr[b] = 1.0;
            em->P[b] = 0.0;
            em->m[b] = 0.0;
            continue;
        }

        double inv_n = 1.0 / count[b];
        double g00_l4_avg = sum_g00_l4[b] * inv_n;
        double grr_l4_avg = sum_grr_l4[b] * inv_n;

        double r = em->r[b];
        double r2 = r * r;

        /* Sturm-Liouville coefficients for hedgehog radial perturbations.
         * Both P and m use PERPENDICULAR (angular) currents only, because for
         * the hedgehog breathing mode, [A₀, A_r] = 0 (temporal perturbation
         * is aligned with radial current in isospin space). */

        /* L₆ contribution: β₆ adds equally to P and m (ratio 1:1, not 2:1)
         * because L₆ = λ₆ B^μ B_μ is Lorentz-invariant with B^r ∝ ḟ sin²f/r²
         * having the same structure as B⁰ ∝ f'sin²f/r². */
        double beta6_avg = sum_beta6[b] * inv_n;

        /* P(r) = 2r² + 4c₄×<|A_⊥|²/|q|⁴> + β₆  (radial stiffness) */
        em->P[b] = 2.0 * r2 + 4.0 * grr_l4_avg + beta6_avg;

        /* m(r) = r² + 2c₄×<|A_⊥|²/|q|⁴> + β₆  (temporal inertia) */
        em->m[b] = r2 + 2.0 * grr_l4_avg + beta6_avg;

        /* BLV metric: g^{00} uses total current (for general perturbations),
         * g^{rr} uses perpendicular current */
        em->g00[b] = r2 + 2.0 * g00_l4_avg;
        em->grr[b] = em->P[b];
    }

    free(count);
    free(sum_g00_l2); free(sum_g00_l4);
    free(sum_grr_l2); free(sum_grr_l4);
    free(sum_beta6);

    return em;
}

void hf_metric_free(EffectiveMetric *em)
{
    if (!em) return;
    free(em->r);
    free(em->g00);
    free(em->grr);
    free(em->P);
    free(em->m);
    free(em);
}
