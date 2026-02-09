/*
 * coupling.h — Degenerate sector coupling for CHPT solitons
 *
 * Three new energy terms that give the degenerate sector (j1,j2,j3,p)
 * non-trivial dynamics and coupling to the bulk soliton:
 *
 * 1. E_{2,D} = (1/2) sum_d |nabla p|^2
 *    Free kinetic energy for weight components.
 *
 * 2. E_{4,C} = (1/4e^2) sum_{i<j} |weight([R_i, R_j])|^2
 *    Skyrme cross-terms from full 8D right-currents. Repulsive near soliton.
 *
 * 3. E_int = (g^2/2) |q|^2 |nabla p|^2
 *    Bulk-degenerate gradient coupling. Attractive well at finite lambda.
 *
 * Combined: Lennard-Jones-like potential (short-range repulsion + long-range attraction).
 */

#ifndef COUPLING_H
#define COUPLING_H

#include "field.h"

/* Coupling energy decomposition */
typedef struct {
    double E2D;     /* degenerate gradient energy */
    double E4C;     /* Skyrme cross-coupling energy */
    double Eint;    /* bulk-degenerate interaction energy */
    double Etotal;  /* sum */
} CouplingEnergy;

/* Compute coupling energy and its decomposition.
 * g_coupling is the interaction strength for E_int.
 * All terms use e_skyrme from params for the Skyrme prefactor. */
CouplingEnergy coupling_energy(const Field *f, const Params *p, double g_coupling);

/* Compute the negative gradient of the coupling energy: force = -dE_coupling/dPsi.
 * The force is ADDED to the existing force array (not zeroed first).
 * This allows calling field_gradient() then coupling_gradient() and summing. */
void coupling_gradient(const Field *f, const Params *p, double g_coupling,
                       Multivector *force);

#endif /* COUPLING_H */
