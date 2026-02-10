/*
 * topology.h — Topological invariant computation
 *
 * Computes:
 *   B: Skyrmion charge (baryon number) — maps S³→S³, π₃(S³)=Z
 *   H: Hopf invariant — maps S³→S², π₃(S²)=Z
 *
 * The Hopf invariant requires the CP¹ projection: q → n̂ = f/|f|
 * where q = s + f₁σ₁ + f₂σ₂ + f₃σ₃ and n̂ ∈ S².
 */

#ifndef TOPOLOGY_H
#define TOPOLOGY_H

#include "spherical_grid.h"

/* Skyrmion charge B = ∫ B⁰ d³x (already in hopfion_field.h, declared here too) */
double topo_skyrmion_charge(const SphericalGrid *g);

/* CP¹ unit vector at cell (i,j,k): n̂ = (f₁,f₂,f₃)/|f|
 * Returns 0 if |f| ≈ 0 (at vacuum s=ρ₀, f=0, n̂ undefined) */
typedef struct { double n1, n2, n3; } UnitVec;
UnitVec topo_cp1_vector(const SphericalGrid *g, int i, int j, int k);

/* Hopf invariant via the Whitehead integral:
 * H = (1/16π²) ∫ ε^{ijk} F_{ij} A_k d³x
 * where F_{ij} = n̂·(∂_i n̂ × ∂_j n̂) and ∇×A = F (solved via FFT).
 *
 * Returns the Hopf charge H (should be integer for topological configurations).
 * This is O(N³ log N) due to the FFT solve.
 *
 * NOTE: This requires FFTW3 for the Poisson solve. If FFTW is not available,
 * use topo_hopf_charge_direct() which is O(N⁶) but needs no external library.
 */
double topo_hopf_charge(const SphericalGrid *g);

/* Direct (slow) Hopf charge via the linking number formula.
 * Samples two preimage curves and computes their linking number.
 * Only accurate for simple configurations; use topo_hopf_charge() for general case. */
double topo_hopf_charge_direct(const SphericalGrid *g, UnitVec target1, UnitVec target2);

/* Compute the area-form pullback F_{ij} = n̂·(∂_i n̂ × ∂_j n̂)
 * at cell (i,j,k). Returns F_{12}, F_{13}, F_{23} (antisymmetric). */
void topo_area_form(const SphericalGrid *g, int i, int j, int k,
                    double *F12, double *F13, double *F23);

/* Integrated (absolute) area form = topological charge density of the CP¹ map */
double topo_cp1_charge_density(const SphericalGrid *g, int i, int j, int k);

#endif /* TOPOLOGY_H */
