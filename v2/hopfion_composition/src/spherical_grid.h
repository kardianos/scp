/*
 * spherical_grid.h — Masked-Cartesian spherical grid
 *
 * A uniform Cartesian grid [-L,L]^3 where only cells with |x| <= R are active.
 * The active cells are enumerated in a flat list for efficient iteration.
 * Cells outside the sphere are clamped to vacuum.
 *
 * This saves ~48% memory vs a full cube and provides uniform boundary distance.
 */

#ifndef SPHERICAL_GRID_H
#define SPHERICAL_GRID_H

#include "clifford.h"
#include <stdlib.h>
#include <math.h>
#include <string.h>

/* Grid parameters */
typedef struct {
    int N;              /* grid points per dimension (underlying cube) */
    double L;           /* half-width of cube: domain is [-L,L]^3 */
    double h;           /* grid spacing = 2L/N */
    double R;           /* sphere radius (R <= L) */
    double R_sponge;    /* inner radius of sponge layer */

    /* Full cube storage (dense, for stencil access) */
    Multivector *psi;   /* field values: psi[flat_idx(i,j,k)] */
    Multivector *vel;   /* time derivative: dq/dt */
    Multivector *force; /* force: -dE/dpsi */

    /* Active cell enumeration */
    int *active;        /* active[n] = flat index of n-th active cell */
    int n_active;       /* number of active cells */

    /* Mask: 1 = active (inside sphere), 0 = boundary/outside */
    unsigned char *mask;

    /* Sponge damping coefficient: gamma[flat_idx] */
    double *sponge;     /* 0 inside R_sponge, ramps to 1 at R */

    /* Vacuum state */
    Multivector vacuum;
} SphericalGrid;

/* Flat index from (i,j,k) — NO periodic wrapping, caller must bounds-check */
static inline int sg_idx(int N, int i, int j, int k)
{
    return i * N * N + j * N + k;
}

/* Position of cell (i,j,k) */
static inline void sg_pos(const SphericalGrid *g, int i, int j, int k,
                           double *x, double *y, double *z)
{
    *x = -g->L + (i + 0.5) * g->h;
    *y = -g->L + (j + 0.5) * g->h;
    *z = -g->L + (k + 0.5) * g->h;
}

/* Radius of cell (i,j,k) */
static inline double sg_radius(const SphericalGrid *g, int i, int j, int k)
{
    double x, y, z;
    sg_pos(g, i, j, k, &x, &y, &z);
    return sqrt(x*x + y*y + z*z);
}

/* Check if (i,j,k) is in bounds AND active */
static inline int sg_is_active(const SphericalGrid *g, int i, int j, int k)
{
    if (i < 0 || i >= g->N || j < 0 || j >= g->N || k < 0 || k >= g->N)
        return 0;
    return g->mask[sg_idx(g->N, i, j, k)];
}

/* Safe field access: returns vacuum for out-of-bounds or inactive cells */
static inline Multivector sg_get(const SphericalGrid *g, int i, int j, int k)
{
    if (i < 0 || i >= g->N || j < 0 || j >= g->N || k < 0 || k >= g->N)
        return g->vacuum;
    int ix = sg_idx(g->N, i, j, k);
    if (!g->mask[ix])
        return g->vacuum;
    return g->psi[ix];
}

/* Allocate grid */
SphericalGrid *sg_alloc(int N, double L, double R, double sponge_width);

/* Free grid */
void sg_free(SphericalGrid *g);

/* Initialize all cells to vacuum */
void sg_set_vacuum(SphericalGrid *g, double rho0);

/* Save/load field snapshot */
int sg_save(const SphericalGrid *g, const char *filename);
int sg_load(SphericalGrid *g, const char *filename);

/* Extract (i,j,k) from flat index */
static inline void sg_unflatten(int N, int flat, int *i, int *j, int *k)
{
    *i = flat / (N * N);
    *j = (flat / N) % N;
    *k = flat % N;
}

#endif /* SPHERICAL_GRID_H */
