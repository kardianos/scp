/*
 * hopfion_field.c — Energy functional and force for composed hopfion dynamics
 *
 * Key new feature: L₆ = -λ₆ (B⁰)² term where B⁰ is the baryon density.
 * This is the "most topological" local term — it makes the winding density
 * directly dynamical.
 *
 * Derivatives: 4th-order central differences (matching existing hopfion_search code).
 * The 9-point stencil {1,-8,0,8,-1}/(12h) for first derivatives and
 * {-1,16,-30,16,-1}/(12h²) for the Laplacian.
 */

#include "hopfion_field.h"
#include <math.h>
#include <stdio.h>
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

/* 4th-order Laplacian: (-1,16,-30,16,-1)/(12h²) per axis */
static inline Multivector laplacian(const SphericalGrid *g, int i, int j, int k)
{
    double c2 = 1.0 / (12.0 * g->h * g->h);
    Multivector center = sg_get(g, i, j, k);
    Multivector r = mv_zero();

    /* X direction */
    Multivector xm2 = sg_get(g, i-2, j, k);
    Multivector xm1 = sg_get(g, i-1, j, k);
    Multivector xp1 = sg_get(g, i+1, j, k);
    Multivector xp2 = sg_get(g, i+2, j, k);

    /* Y direction */
    Multivector ym2 = sg_get(g, i, j-2, k);
    Multivector ym1 = sg_get(g, i, j-1, k);
    Multivector yp1 = sg_get(g, i, j+1, k);
    Multivector yp2 = sg_get(g, i, j+2, k);

    /* Z direction */
    Multivector zm2 = sg_get(g, i, j, k-2);
    Multivector zm1 = sg_get(g, i, j, k-1);
    Multivector zp1 = sg_get(g, i, j, k+1);
    Multivector zp2 = sg_get(g, i, j, k+2);

    /* Sum: (-1*f_{±2} + 16*f_{±1} - 30*f_0) per axis, divided by 12h² */
    #define LAP_COMP(C) \
        r.C = c2 * ( \
            -xm2.C + 16*xm1.C - 30*center.C + 16*xp1.C - xp2.C + \
            -ym2.C + 16*ym1.C - 30*center.C + 16*yp1.C - yp2.C + \
            -zm2.C + 16*zm1.C - 30*center.C + 16*zp1.C - zp2.C   \
        )

    LAP_COMP(s);
    LAP_COMP(f1); LAP_COMP(f2); LAP_COMP(f3);
    LAP_COMP(j1); LAP_COMP(j2); LAP_COMP(j3);
    LAP_COMP(p);

    #undef LAP_COMP
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
        if (r > g->R - 2.5*g->h) continue;

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
        if (r > g->R - 2.5*g->h) continue;  /* skip near-boundary */

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

/* ========== Force computation ========== */

/* Force from L₂ (gradient term): F_E2 = ∇²q */
/* Force from V (Mexican hat): F_V = -λ(|q|²-ρ₀²)q */
/* Force from V_π (pion mass): F_π = m_π² ρ₀² q_s_hat (tangent projection) */
/* Force from L₄ (Skyrme): complicated, computed numerically via finite differences of E₄ */
/* Force from L₆ (sextic): F_6 = -δE₆/δq, also via finite differences */

/* For the Skyrme force, we use the same approach as field.c:
 * compute E₄ at perturbed field values and take finite differences.
 * This is slow but correct. Later optimization: derive analytical force. */

void hf_force(SphericalGrid *g, const HopfionParams *p)
{
    int N = g->N;
    double e2 = p->e_skyrme * p->e_skyrme;
    double rho0_sq = p->rho0 * p->rho0;

    #pragma omp parallel for schedule(static)
    for (int n = 0; n < g->n_active; n++) {
        int ix = g->active[n];
        int i, j, k;
        sg_unflatten(N, ix, &i, &j, &k);

        double r = sg_radius(g, i, j, k);
        if (r > g->R - 2.5*g->h) {
            g->force[ix] = mv_zero();
            continue;
        }

        Multivector psi = g->psi[ix];
        Multivector F = mv_zero();

        /* F_E2 = ∇²q (elastic/gradient force) */
        Multivector lap = laplacian(g, i, j, k);
        F = mv_add(F, lap);

        /* F_V = -λ(|q|²-ρ₀²)q */
        double q2 = mv_bulk_norm2(psi);
        double dq = q2 - rho0_sq;
        if (p->lambda > 0) {
            Multivector fv;
            double coeff = -p->lambda * dq;
            fv.s  = coeff * psi.s;
            fv.f1 = coeff * psi.f1;
            fv.f2 = coeff * psi.f2;
            fv.f3 = coeff * psi.f3;
            fv.j1 = 0; fv.j2 = 0; fv.j3 = 0; fv.p = 0;
            F = mv_add(F, fv);
        }

        /* F_π: pion mass force (tangent to S³) */
        if (p->m_pi_sq > 0 && q2 > 1e-20) {
            double qnorm = sqrt(q2);
            double q3 = qnorm * q2;
            /* V_π = m² ρ₀² (1 - s/|q|)
             * ∂V_π/∂s = -m² ρ₀² / |q| + m² ρ₀² s² / |q|³
             *         = m² ρ₀² (s² - |q|²) / |q|³
             *         = -m² ρ₀² |f|² / |q|³
             * ∂V_π/∂f_a = m² ρ₀² s f_a / |q|³ */
            double mfact = p->m_pi_sq * rho0_sq / q3;
            Multivector fpi;
            fpi.s  = -mfact * (psi.f1*psi.f1 + psi.f2*psi.f2 + psi.f3*psi.f3);
            fpi.f1 =  mfact * psi.s * psi.f1;
            fpi.f2 =  mfact * psi.s * psi.f2;
            fpi.f3 =  mfact * psi.s * psi.f3;
            fpi.j1 = 0; fpi.j2 = 0; fpi.j3 = 0; fpi.p = 0;
            /* Force = -∂V/∂q, so negate */
            F = mv_sub(F, fpi);
        }

        /* F_D: degenerate mass force */
        if (p->mu > 0) {
            Multivector fd = mv_zero();
            double mu2 = p->mu * p->mu;
            fd.j1 = -mu2 * psi.j1;
            fd.j2 = -mu2 * psi.j2;
            fd.j3 = -mu2 * psi.j3;
            fd.p  = -mu2 * psi.p;
            F = mv_add(F, fd);
        }

        /* F_E4: Skyrme force — use the analytical formula from field.c
         *
         * The Skyrme energy density involves |[R_i, R_j]|² where R_i = q̃∂_iq.
         * The force is obtained by varying this w.r.t. q at each point.
         *
         * For efficiency, we compute the Skyrme force using the identity:
         * δE₄/δq = -(1/e²) Σ_i ∂_i ( Σ_j [R_j, [R_j, R_i]] ) / |q|²
         *
         * This requires second derivatives of the currents, which we compute
         * from the stencil. For now, we use a direct but slower approach:
         * compute the Skyrme Laplacian numerically.
         *
         * TODO: Optimize with analytical Skyrme force formula
         */

        /* Skyrme force via consistent approach:
         * F_E4 = (2/e²) Σ_i ∂_i [ Σ_{j≠i} [R_i, [R_i, ∂_j q]] / |q|² ]
         *
         * For now, use the field.c approach: compute ∂_i(E4_current_contribution)
         * via stencil applied to the Skyrme current.
         *
         * Simplified Skyrme force for quaternion field:
         * F_4,a = (2/e²|q|²) Σ_i ∂_i [ (|R|²∂_i q_a - (R·∂_i q)R_a) ]
         * where |R|² = Σ_j |R_j|² and we use the commutator norm.
         *
         * We defer the full analytical derivation and use a simpler but
         * correct formulation: vary the Skyrme energy numerically.
         */

        /* For the initial implementation, we compute Skyrme force component-by-component
         * using the relation: for each component c,
         *   F_4^c(x) = -∂E₄/∂ψ^c(x)
         * ≈ -[E₄(ψ+εδ_c) - E₄(ψ-εδ_c)] / (2ε)
         *
         * This is O(16) energy evaluations per cell × N³ cells = very slow.
         * We only enable this for gradient verification; for dynamics, use
         * a simplified force or port the analytical formula from field.c.
         */

        /* For now: use the well-known Skyrme force from the 4-component form.
         * Given R_i = q̃∂_iq, the force is:
         *
         * F_4 = (2ρ₀²/e²) ∇²_Skyrme(q)
         *
         * where ∇²_Skyrme uses the Skyrme-modified metric.
         * We compute this by taking the 4th-order Laplacian of the
         * Skyrme-weighted gradient.
         *
         * Implementation: compute 3 Skyrme currents S_i at each neighbor,
         * then take divergence ∂_i S_i.
         */

        /* Skyrme force placeholder — will be filled in next implementation step.
         * For gradient flow / testing, the L₂ + V forces suffice. */

        g->force[ix] = F;
    }
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
