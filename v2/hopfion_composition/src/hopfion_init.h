/*
 * hopfion_init.h — Hopfion and Skyrmion initialization
 *
 * Provides routines to initialize:
 *   1. Single Skyrmion (hedgehog, from 1D radial profile)
 *   2. Single Hopfion (axially symmetric, Hopf map)
 *   3. (n,m) Hopfion (torus knot, Hopf index = n×m)
 *   4. Composed state: multiple hopfions at different positions
 *   5. Hybrid: Skyrmion core + hopfion ring (FeGe-like)
 */

#ifndef HOPFION_INIT_H
#define HOPFION_INIT_H

#include "spherical_grid.h"

/* Radial profile for 1D Skyrmion initialization */
typedef struct {
    double *r;
    double *f;
    double *rho;    /* NULL for σ-model (|q|=ρ₀) */
    int n;
    double dr;
    double r_max;
} RadialProfile;

RadialProfile *profile_load(const char *filename);
void profile_free(RadialProfile *p);
double profile_interp_f(const RadialProfile *p, double r);
double profile_interp_rho(const RadialProfile *p, double rho0, double r);

/* ========== Single object initialization ========== */

/* Initialize a B=1 hedgehog Skyrmion centered at (cx,cy,cz)
 * using a 1D radial profile f(r). */
void init_skyrmion(SphericalGrid *g, const RadialProfile *prof, double rho0,
                   double cx, double cy, double cz);

/* Initialize an H=1 axially symmetric hopfion centered at (cx,cy,cz)
 * with size parameter 'a' (controls ring radius).
 *
 * The standard Hopf map: S³ → S²
 *   z₁ = (a² - r² + 2iaz) / (a² + r²)
 *   z₂ = 2a(x + iy) / (a² + r²)
 *   n̂ = (2 Re(z₁z̄₂), 2 Im(z₁z̄₂), |z₁|²-|z₂|²)
 *
 * The quaternion field is constructed from n̂ via the inverse Hopf map:
 *   q = ρ₀ (cos(g/2) + sin(g/2) n̂·σ)
 * where g(r) is a radial profile satisfying g(0)=π, g(∞)=0.
 */
void init_hopfion(SphericalGrid *g, double rho0, double a,
                  double cx, double cy, double cz);

/* Initialize an (n,m) hopfion with Hopf index H = n×m.
 * The n̂ field wraps n times toroidally and m times poloidally.
 * For (1,1): standard hopfion (unknotted ring)
 * For (2,1): H=2 (linked unknots)
 * For (2,3): H=6 (trefoil-like) */
void init_hopfion_nm(SphericalGrid *g, double rho0, double a,
                     int n_tor, int m_pol,
                     double cx, double cy, double cz);

/* ========== Composed state initialization ========== */

/* Compose two quaternion fields via the product ansatz:
 * q_total = q₁ × q₂ / ρ₀
 * This preserves topology: B_total = B₁ + B₂ (approximately, for separation >> size)
 * For hopfions: the preimage curves of q₁ and q₂ will link, contributing to H_total.
 */
void init_compose_product(SphericalGrid *g, double rho0);
/* Note: call individual init functions first to set up individual fields,
 * then compose. Use a temporary grid for the second field. */

/* Initialize two Skyrmions at ±z₀ on the z-axis with initial velocities ±v_z */
void init_two_skyrmions(SphericalGrid *g, const RadialProfile *prof, double rho0,
                        double z0, double vz);

/* Initialize a Skyrmion + hopfion ring composite:
 * Skyrmion (B=1) at center, hopfion ring at distance d from center.
 * This mimics the (Q,H) coupled states observed in FeGe. */
void init_skyrmion_hopfion(SphericalGrid *g, const RadialProfile *prof,
                            double rho0, double hopfion_a, double hopfion_d);

/* ========== Velocity initialization ========== */

/* Boost a localized field configuration centered at (cx,cy,cz)
 * to velocity (vx,vy,vz). Sets g->vel appropriately.
 * Uses the non-relativistic approximation v(x) ≈ -v·∇q for |v| << c. */
void init_boost(SphericalGrid *g, double rho0,
                double cx, double cy, double cz,
                double vx, double vy, double vz);

#endif /* HOPFION_INIT_H */
