/*
 * hopfion_field.h — Energy functional and force computation for hopfion simulator
 *
 * Lagrangian terms:
 *   L₂ = (1/2)|∂q|²                          (quadratic gradient)
 *   L₄ = (1/4e²) Σ|[R_μ,R_ν]|²              (Skyrme quartic)
 *   L₆ = -λ₆ (B⁰)²                          (sextic/BPS baryon current squared)
 *   V  = (λ/4)(|q|²-ρ₀²)²                   (Mexican hat potential)
 *   V_π = m_π² ρ₀² (1 - q·s/|q|)            (pion mass)
 *
 * The force F = -δE/δq is computed for leapfrog time evolution.
 */

#ifndef HOPFION_FIELD_H
#define HOPFION_FIELD_H

#include "spherical_grid.h"

/* Physical parameters */
typedef struct {
    double rho0;        /* vacuum density */
    double e_skyrme;    /* Skyrme coupling (dimensionless) */
    double lambda;      /* bulk self-coupling (Mexican hat) */
    double lambda6;     /* sextic (BPS) coupling */
    double m_pi_sq;     /* pion mass squared */
    double mu;          /* degenerate sector mass */
    double c;           /* speed of light */
} HopfionParams;

/* Energy decomposition */
typedef struct {
    double E2;          /* gradient (quadratic) */
    double E4;          /* Skyrme (quartic) */
    double E6;          /* sextic (BPS) */
    double EV;          /* Mexican hat potential */
    double Epi;         /* pion mass potential */
    double ED;          /* degenerate sector mass */
    double Etotal;
} HopfionEnergy;

/* Compute total energy and decomposition */
HopfionEnergy hf_energy(const SphericalGrid *g, const HopfionParams *p);

/* Compute kinetic energy (1/2c²)|dq/dt|² */
double hf_kinetic_energy(const SphericalGrid *g, const HopfionParams *p);

/* Compute force F = -δE/δq at all active cells.
 * Result stored in g->force[]. */
void hf_force(SphericalGrid *g, const HopfionParams *p);

/* Compute baryon density B⁰(x) at cell (i,j,k).
 * Uses 4th-order central differences for derivatives. */
double hf_baryon_density(const SphericalGrid *g, int i, int j, int k);

/* Compute total baryon (Skyrmion) charge B = ∫ B⁰ d³x */
double hf_baryon_charge(const SphericalGrid *g);

/* Project field onto sigma-model constraint |q| = rho0
 * (only bulk quaternion; degenerate sector left unchanged) */
void hf_sigma_project(SphericalGrid *g, double rho0);

/* ========== Effective Metric ========== */

/* Radial profile of the BLV effective metric.
 *
 * For L₂ + L₄ Skyrme theory, the effective acoustic metric governing
 * small fluctuations around a background field φ₀ is:
 *
 *   g_eff^{μν} ∝ -∂²L/∂(∂_μφ^A)∂(∂_νφ^B)  evaluated on φ₀
 *
 * For a hedgehog soliton, this reduces to radial functions:
 *   P(r) = 2r² + 4c₄ sin²f(r)   (spatial radial)
 *   m(r) = r² + 2c₄ sin²f(r)    (temporal/inertia)
 *
 * Key result: For the sigma model (|q|=ρ₀), P(r)/m(r) = 2 identically,
 * meaning NO radial time dilation. Emergent gravity requires finite-λ
 * (varying ρ(r)) or L₆ to break this identity.
 */
typedef struct {
    int n_bins;         /* number of radial bins */
    double dr;          /* bin width */
    double *r;          /* bin centers: r[i] = (i+0.5)*dr */
    double *g00;        /* -g_00 effective metric (temporal) */
    double *grr;        /* g_rr effective metric (radial spatial) */
    double *P;          /* Sturm-Liouville P(r) = radial stiffness */
    double *m;          /* mass/inertia m(r) */
} EffectiveMetric;

/* Compute spherically-averaged effective metric from a 3D field configuration.
 * The metric is computed from L₂ + L₄ contributions.
 * Returns allocated EffectiveMetric (caller must free with hf_metric_free). */
EffectiveMetric *hf_effective_metric(const SphericalGrid *g, const HopfionParams *p,
                                      int n_bins);

/* Free effective metric */
void hf_metric_free(EffectiveMetric *em);

#endif /* HOPFION_FIELD_H */
