/*
 * initial.c — Initial field configurations
 */

#define _USE_MATH_DEFINES
#include <math.h>
#include <stdlib.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#include "initial.h"

void init_vacuum(Field *f, const Params *p)
{
    int N = f->N;
    for (int i = 0; i < N*N*N; i++) {
        f->psi[i] = mv_zero();
        f->psi[i].s = p->rho0;
    }
}

void init_hedgehog(Field *f, const Params *p, double R)
{
    int N = f->N;
    double rho0 = p->rho0;

    for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
    for (int k = 0; k < N; k++) {
        int ix = idx(N, i, j, k);

        /* Physical coordinates (centered) */
        double x = -f->L + (i + 0.5) * f->h;
        double y = -f->L + (j + 0.5) * f->h;
        double z = -f->L + (k + 0.5) * f->h;
        double r2 = x*x + y*y + z*z;
        double r = sqrt(r2);

        /* Profile function: f(r) = pi * R^2 / (r^2 + R^2) */
        double profile = M_PI * R * R / (r2 + R * R);

        double cos_f = cos(profile);
        double sin_f = sin(profile);

        /* Hedgehog: q = rho0 * [cos(f) + sin(f) * (r_hat . bivectors)]
         * r_hat = (x, y, z) / r
         * Bivector part: x/r * e23 + y/r * e31 + z/r * e12
         *
         * At r=0: profile=pi, so q = rho0*(-1) = -rho0 (anti-vacuum)
         * At r=inf: profile=0, so q = rho0*(+1) = +rho0 (vacuum)
         */
        f->psi[ix].s = rho0 * cos_f;

        if (r > 1e-12) {
            double sr = sin_f / r;
            f->psi[ix].f1 = rho0 * sr * x;  /* e23 component */
            f->psi[ix].f2 = rho0 * sr * y;  /* e31 component */
            f->psi[ix].f3 = rho0 * sr * z;  /* e12 component */
        } else {
            /* At origin: sin(pi)/0 is 0/0, but limit is well-defined
             * sin(f(r))/r -> pi*R^2 / R^2 * ... Actually at r=0:
             * f(0) = pi, sin(pi) = 0, so the bivector part is 0. */
            f->psi[ix].f1 = 0;
            f->psi[ix].f2 = 0;
            f->psi[ix].f3 = 0;
        }

        /* Degenerate sector: zero */
        f->psi[ix].j1 = 0;
        f->psi[ix].j2 = 0;
        f->psi[ix].j3 = 0;
        f->psi[ix].p  = 0;
    }
}

void init_perturb(Field *f, double amplitude, unsigned int seed)
{
    srand(seed);
    int N = f->N;
    for (int i = 0; i < N*N*N; i++) {
        double *comp = (double *)&f->psi[i];
        for (int a = 0; a < 8; a++) {
            comp[a] += amplitude * (2.0 * rand() / RAND_MAX - 1.0);
        }
    }
}
