/*
 * topology.c — Topological invariant computation
 *
 * Skyrmion charge (B): standard formula from baryon density
 * Hopf charge (H): via area-form pullback F_{ij} = n̂·(∂_i n̂ × ∂_j n̂)
 *                   and Gauss linking integral
 */

#include "topology.h"
#include "hopfion_field.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ========== Skyrmion charge ========== */

double topo_skyrmion_charge(const SphericalGrid *g)
{
    return hf_baryon_charge(g);
}

/* ========== CP¹ projection ========== */

UnitVec topo_cp1_vector(const SphericalGrid *g, int i, int j, int k)
{
    Multivector psi = sg_get(g, i, j, k);
    double f2 = psi.f1*psi.f1 + psi.f2*psi.f2 + psi.f3*psi.f3;
    UnitVec n = {0, 0, 1};  /* default: north pole when f=0 */

    if (f2 > 1e-20) {
        double fnorm = sqrt(f2);
        n.n1 = psi.f1 / fnorm;
        n.n2 = psi.f2 / fnorm;
        n.n3 = psi.f3 / fnorm;
    }
    return n;
}

/* ========== Area form pullback ========== */

/* 4th-order central difference of n̂ in each direction */
static void cp1_derivatives(const SphericalGrid *g, int i, int j, int k,
                            double dn_dx[3], double dn_dy[3], double dn_dz[3])
{
    double c1 = 1.0 / (12.0 * g->h);

    /* X derivatives */
    UnitVec nxm2 = topo_cp1_vector(g, i-2, j, k);
    UnitVec nxm1 = topo_cp1_vector(g, i-1, j, k);
    UnitVec nxp1 = topo_cp1_vector(g, i+1, j, k);
    UnitVec nxp2 = topo_cp1_vector(g, i+2, j, k);
    dn_dx[0] = c1 * (nxm2.n1 - 8*nxm1.n1 + 8*nxp1.n1 - nxp2.n1);
    dn_dx[1] = c1 * (nxm2.n2 - 8*nxm1.n2 + 8*nxp1.n2 - nxp2.n2);
    dn_dx[2] = c1 * (nxm2.n3 - 8*nxm1.n3 + 8*nxp1.n3 - nxp2.n3);

    /* Y derivatives */
    UnitVec nym2 = topo_cp1_vector(g, i, j-2, k);
    UnitVec nym1 = topo_cp1_vector(g, i, j-1, k);
    UnitVec nyp1 = topo_cp1_vector(g, i, j+1, k);
    UnitVec nyp2 = topo_cp1_vector(g, i, j+2, k);
    dn_dy[0] = c1 * (nym2.n1 - 8*nym1.n1 + 8*nyp1.n1 - nyp2.n1);
    dn_dy[1] = c1 * (nym2.n2 - 8*nym1.n2 + 8*nyp1.n2 - nyp2.n2);
    dn_dy[2] = c1 * (nym2.n3 - 8*nym1.n3 + 8*nyp1.n3 - nyp2.n3);

    /* Z derivatives */
    UnitVec nzm2 = topo_cp1_vector(g, i, j, k-2);
    UnitVec nzm1 = topo_cp1_vector(g, i, j, k-1);
    UnitVec nzp1 = topo_cp1_vector(g, i, j, k+1);
    UnitVec nzp2 = topo_cp1_vector(g, i, j, k+2);
    dn_dz[0] = c1 * (nzm2.n1 - 8*nzm1.n1 + 8*nzp1.n1 - nzp2.n1);
    dn_dz[1] = c1 * (nzm2.n2 - 8*nzm1.n2 + 8*nzp1.n2 - nzp2.n2);
    dn_dz[2] = c1 * (nzm2.n3 - 8*nzm1.n3 + 8*nzp1.n3 - nzp2.n3);
}

void topo_area_form(const SphericalGrid *g, int i, int j, int k,
                    double *F12, double *F13, double *F23)
{
    UnitVec n = topo_cp1_vector(g, i, j, k);
    double dn_dx[3], dn_dy[3], dn_dz[3];
    cp1_derivatives(g, i, j, k, dn_dx, dn_dy, dn_dz);

    /* F_{ij} = n̂ · (∂_i n̂ × ∂_j n̂)
     * F_{12} = n̂ · (∂_x n̂ × ∂_y n̂) */
    double cross_xy[3] = {
        dn_dx[1]*dn_dy[2] - dn_dx[2]*dn_dy[1],
        dn_dx[2]*dn_dy[0] - dn_dx[0]*dn_dy[2],
        dn_dx[0]*dn_dy[1] - dn_dx[1]*dn_dy[0]
    };
    *F12 = n.n1*cross_xy[0] + n.n2*cross_xy[1] + n.n3*cross_xy[2];

    double cross_xz[3] = {
        dn_dx[1]*dn_dz[2] - dn_dx[2]*dn_dz[1],
        dn_dx[2]*dn_dz[0] - dn_dx[0]*dn_dz[2],
        dn_dx[0]*dn_dz[1] - dn_dx[1]*dn_dz[0]
    };
    *F13 = n.n1*cross_xz[0] + n.n2*cross_xz[1] + n.n3*cross_xz[2];

    double cross_yz[3] = {
        dn_dy[1]*dn_dz[2] - dn_dy[2]*dn_dz[1],
        dn_dy[2]*dn_dz[0] - dn_dy[0]*dn_dz[2],
        dn_dy[0]*dn_dz[1] - dn_dy[1]*dn_dz[0]
    };
    *F23 = n.n1*cross_yz[0] + n.n2*cross_yz[1] + n.n3*cross_yz[2];
}

double topo_cp1_charge_density(const SphericalGrid *g, int i, int j, int k)
{
    double F12, F13, F23;
    topo_area_form(g, i, j, k, &F12, &F13, &F23);
    /* This is not really a "charge density" in the Hopf sense,
     * but rather |F|² which measures how much S² area is being covered.
     * The actual Hopf charge requires the A∧F integral. */
    return sqrt(F12*F12 + F13*F13 + F23*F23);
}

/* ========== Hopf charge (FFT-free version) ========== */

/* Hopf charge via discretized A∧F integral.
 *
 * Strategy: We compute the "magnetic field" B_i = (1/2)ε_{ijk} F_{jk}
 * and solve ∇²A = -B (vector Poisson equation) on the grid using
 * Jacobi iteration. Then H = (1/4π²) ∫ A·B d³x.
 *
 * For an FFT-free implementation, we use the relaxation method.
 * This is O(N³ × N_iter) where N_iter ~ N for convergence.
 * Slow but correct and requires no external libraries.
 *
 * For production use, replace with FFTW-based solver.
 */
double topo_hopf_charge(const SphericalGrid *g)
{
    int N = g->N;
    long total = (long)N * N * N;
    double h = g->h;
    double h3 = h * h * h;

    /* Allocate B field (divergence-free "magnetic" field from area form) */
    double *Bx = calloc(total, sizeof(double));
    double *By = calloc(total, sizeof(double));
    double *Bz = calloc(total, sizeof(double));

    /* Allocate vector potential A */
    double *Ax = calloc(total, sizeof(double));
    double *Ay = calloc(total, sizeof(double));
    double *Az = calloc(total, sizeof(double));

    if (!Bx || !By || !Bz || !Ax || !Ay || !Az) {
        fprintf(stderr, "topo_hopf_charge: out of memory\n");
        free(Bx); free(By); free(Bz);
        free(Ax); free(Ay); free(Az);
        return 0.0;
    }

    /* Step 1: Compute B_i = (1/2)ε_{ijk} F_{jk} at each active cell */
    #pragma omp parallel for schedule(static)
    for (int n = 0; n < g->n_active; n++) {
        int ix = g->active[n];
        int i, j, k;
        sg_unflatten(N, ix, &i, &j, &k);

        double r = sg_radius(g, i, j, k);
        if (r > g->R - 2.5*h) continue;

        double F12, F13, F23;
        topo_area_form(g, i, j, k, &F12, &F13, &F23);

        /* B_x = F_{23}, B_y = -F_{13} = F_{31}, B_z = F_{12} */
        Bx[ix] = F23;
        By[ix] = -F13;
        Bz[ix] = F12;
    }

    /* Step 2: Solve ∇²A = -B via Jacobi relaxation.
     * Use 2nd-order Laplacian (6-point stencil) for speed. */
    int max_iter = 200;  /* limited iterations for speed */
    double omega = 1.5;  /* SOR over-relaxation factor */

    for (int iter = 0; iter < max_iter; iter++) {
        #pragma omp parallel for schedule(static)
        for (int n = 0; n < g->n_active; n++) {
            int ix = g->active[n];
            int i, j, k;
            sg_unflatten(N, ix, &i, &j, &k);

            double r = sg_radius(g, i, j, k);
            if (r > g->R - 1.5*h) continue;

            if (!sg_is_active(g, i-1,j,k) || !sg_is_active(g, i+1,j,k) ||
                !sg_is_active(g, i,j-1,k) || !sg_is_active(g, i,j+1,k) ||
                !sg_is_active(g, i,j,k-1) || !sg_is_active(g, i,j,k+1))
                continue;

            /* 6-point Laplacian: (sum of 6 neighbors - 6*center) / h² */
            int ixm = sg_idx(N, i-1, j, k);
            int ixp = sg_idx(N, i+1, j, k);
            int iym = sg_idx(N, i, j-1, k);
            int iyp = sg_idx(N, i, j+1, k);
            int izm = sg_idx(N, i, j, k-1);
            int izp = sg_idx(N, i, j, k+1);

            double new_Ax = (Ax[ixm]+Ax[ixp]+Ax[iym]+Ax[iyp]+Ax[izm]+Ax[izp] + h*h*Bx[ix]) / 6.0;
            double new_Ay = (Ay[ixm]+Ay[ixp]+Ay[iym]+Ay[iyp]+Ay[izm]+Ay[izp] + h*h*By[ix]) / 6.0;
            double new_Az = (Az[ixm]+Az[ixp]+Az[iym]+Az[iyp]+Az[izm]+Az[izp] + h*h*Bz[ix]) / 6.0;

            Ax[ix] = (1-omega)*Ax[ix] + omega*new_Ax;
            Ay[ix] = (1-omega)*Ay[ix] + omega*new_Ay;
            Az[ix] = (1-omega)*Az[ix] + omega*new_Az;
        }
    }

    /* Step 3: H = (1/4π²) ∫ A·B d³x */
    double H = 0;
    #pragma omp parallel for reduction(+:H) schedule(static)
    for (int n = 0; n < g->n_active; n++) {
        int ix = g->active[n];
        H += (Ax[ix]*Bx[ix] + Ay[ix]*By[ix] + Az[ix]*Bz[ix]) * h3;
    }
    H /= (4.0 * M_PI * M_PI);

    free(Bx); free(By); free(Bz);
    free(Ax); free(Ay); free(Az);
    return H;
}

/* Slow direct linking number — not implemented yet */
double topo_hopf_charge_direct(const SphericalGrid *g, UnitVec target1, UnitVec target2)
{
    (void)g; (void)target1; (void)target2;
    fprintf(stderr, "topo_hopf_charge_direct: not yet implemented\n");
    return 0.0;
}
