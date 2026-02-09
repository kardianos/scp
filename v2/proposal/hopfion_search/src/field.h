/*
 * field.h — 3D field grid and energy functional for CHPT
 *
 * The field Psi(x) is stored on a regular 3D grid with periodic BCs.
 * Energy terms: E2 (gradient), E4 (Skyrme), EV (bulk potential), ED (degenerate potential)
 */

#ifndef FIELD_H
#define FIELD_H

#include "clifford.h"

/* Physical parameters */
typedef struct {
    double rho0;    /* vacuum density */
    double lambda;  /* bulk self-coupling */
    double e_skyrme;/* Skyrme coupling */
    double mu;      /* degenerate mass */
    double g_coupling; /* bulk-degenerate gradient coupling strength */
    double c;       /* speed of light (set to 1) */
} Params;

/* 3D grid */
typedef struct {
    int N;          /* grid points per dimension (cubic grid N x N x N) */
    double L;       /* half-width: domain is [-L, L]^3 */
    double h;       /* grid spacing = 2L/N */
    Multivector *psi; /* field values: psi[idx(i,j,k)] */
} Field;

/* Flat index from 3D coordinates (periodic BCs handled by caller) */
static inline int idx(int N, int i, int j, int k)
{
    return ((i % N + N) % N) * N * N
         + ((j % N + N) % N) * N
         + ((k % N + N) % N);
}

/* Energy decomposition */
typedef struct {
    double E2;      /* gradient energy */
    double E4;      /* Skyrme energy */
    double EV;      /* bulk potential energy */
    double ED;      /* degenerate potential energy */
    double Etotal;  /* sum */
} Energy;

/* Allocate/free field */
Field *field_alloc(int N, double L);
void field_free(Field *f);

/* Compute total energy and its decomposition */
Energy field_energy(const Field *f, const Params *p);

/* Compute the negative gradient of the energy: force = -dE/dPsi
 * This is used for gradient flow: dPsi/dtau = force
 * The force array must be pre-allocated with N^3 entries.
 */
void field_gradient(const Field *f, const Params *p, Multivector *force);

/* Compute topological charge (baryon number / degree of map) */
double field_topological_charge(const Field *f, const Params *p);

#endif /* FIELD_H */
