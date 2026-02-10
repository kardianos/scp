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

#endif /* HOPFION_FIELD_H */
