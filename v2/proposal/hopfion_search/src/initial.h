/*
 * initial.h — Initial field configurations for soliton search
 */

#ifndef INITIAL_H
#define INITIAL_H

#include "field.h"

/* Set field to vacuum: Psi = rho0 everywhere */
void init_vacuum(Field *f, const Params *p);

/* Hedgehog ansatz for B=1 Skyrmion (topological charge Q=1).
 *
 * q(r) = rho0 * [cos(f(r)) + sin(f(r)) * (r_hat . bivectors)]
 * where f(r) = pi * R^2 / (r^2 + R^2) is the profile function.
 *
 * f(0) = pi (anti-vacuum at center), f(inf) = 0 (vacuum at infinity).
 * R = soliton size parameter.
 *
 * Degenerate sector: J = P = 0 initially.
 */
void init_hedgehog(Field *f, const Params *p, double R);

/* Add small random perturbation to all 8 components (for testing) */
void init_perturb(Field *f, double amplitude, unsigned int seed);

#endif /* INITIAL_H */
