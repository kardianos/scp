/*
 * verify.c — Verify gradient computation by comparing with finite differences.
 *
 * For each energy term, perturb a single field component and check that
 * the analytical gradient matches the numerical gradient.
 */

#define _USE_MATH_DEFINES
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "field.h"
#include "initial.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

int main(void)
{
    /* Small grid for fast verification */
    int N = 16;
    double L = 4.0;
    /* Test 1: without Skyrme (e very large) to verify E2+EV+ED gradient */
    printf("=== Test 1: No Skyrme (e=1e6) ===\n");
    Params params = {1.0, 10.0, 1e6, 1.0, 1.0};  /* rho0, lambda, e, mu, c */

    Field *f = field_alloc(N, L);
    init_hedgehog(f, &params, 1.5);

    /* Add some degenerate sector perturbation */
    init_perturb(f, 0.01, 42);

    int N3 = N * N * N;
    Multivector *force = (Multivector *)malloc((size_t)N3 * sizeof(Multivector));

    /* Compute analytical gradient */
    field_gradient(f, &params, force);

    /* Test at a few representative points */
    int test_points[] = {idx(N, N/4, N/4, N/4), idx(N, N/2, N/2, N/2),
                         idx(N, 3*N/4, N/4, N/2)};
    int ntest = 3;
    double eps = 1e-5;

    printf("=== Gradient Verification ===\n");
    printf("Comparing analytical force with finite-difference d(-E)/dPsi\n\n");

    int pass = 1;
    for (int t = 0; t < ntest; t++) {
        int ix = test_points[t];
        printf("Point %d (flat index %d):\n", t, ix);
        printf("  Comp | Analytical |   Numerical  | Rel Error\n");
        printf("  -----|------------|--------------|----------\n");

        double *comp = (double *)&f->psi[ix];
        double *fcomp = (double *)&force[ix];

        for (int a = 0; a < 8; a++) {
            double orig = comp[a];

            /* E(+eps) */
            comp[a] = orig + eps;
            Energy ep = field_energy(f, &params);

            /* E(-eps) */
            comp[a] = orig - eps;
            Energy em = field_energy(f, &params);

            /* Restore */
            comp[a] = orig;

            /* Numerical gradient: force = -(E+ - E-) / (2*eps) */
            double num = -(ep.Etotal - em.Etotal) / (2.0 * eps);

            /* Analytical: need to convert from "force density" to "force per point"
             * force[ix] stores -dE_density/dPsi = -(1/h³) dE/dPsi
             * But wait - let me check what the actual convention is.
             *
             * field_gradient computes: force = -d(sum_of_densities)/dPsi
             * E = h³ * sum_of_densities
             * dE/dPsi = h³ * d(sum)/dPsi
             * So: force = -(1/h³) * dE/dPsi
             * And: -dE/dPsi = h³ * force
             *
             * The numerical gradient -(E+ - E-)/(2eps) = -dE/dPsi = h³ * force_analytical
             */
            double h3 = f->h * f->h * f->h;
            double anal = h3 * fcomp[a];

            double rel_err = (fabs(num) > 1e-12) ?
                fabs(anal - num) / fabs(num) : fabs(anal - num);

            const char *names[] = {"s ", "f1", "f2", "f3", "j1", "j2", "j3", "P "};
            printf("  %s  | %+.5e | %+.5e | %.3e %s\n",
                   names[a], anal, num, rel_err, rel_err > 0.05 ? " ***FAIL***" : "");

            if (rel_err > 0.05 && fabs(num) > 1e-10) pass = 0;
        }
        printf("\n");
    }

    /* Also verify topological charge */
    double Q = field_topological_charge(f, &params);
    printf("Topological charge Q = %.6f (expect ~1.0 for hedgehog)\n", Q);

    /* Also check energy decomposition */
    Energy en = field_energy(f, &params);
    printf("\nEnergy: E=%.6e (E2=%.4e E4=%.4e EV=%.4e ED=%.4e)\n",
           en.Etotal, en.E2, en.E4, en.EV, en.ED);

    free(force);
    field_free(f);

    /* ===== Test 2: with Skyrme (e=2) to verify E4 gradient ===== */
    printf("\n=== Test 2: With Skyrme (e=2.0) ===\n");
    Params params2 = {1.0, 10.0, 2.0, 1.0, 1.0};

    Field *f2 = field_alloc(N, L);
    init_hedgehog(f2, &params2, 1.5);
    init_perturb(f2, 0.01, 42);

    Multivector *force2 = (Multivector *)malloc((size_t)N3 * sizeof(Multivector));
    field_gradient(f2, &params2, force2);

    int test_points2[] = {idx(N, N/4, N/4, N/4), idx(N, N/2, N/2, N/2),
                          idx(N, 3*N/4, N/4, N/2)};

    printf("Comparing analytical force with finite-difference d(-E)/dPsi\n\n");

    for (int t = 0; t < ntest; t++) {
        int ix = test_points2[t];
        printf("Point %d (flat index %d):\n", t, ix);
        printf("  Comp | Analytical |   Numerical  | Rel Error\n");
        printf("  -----|------------|--------------|----------\n");

        double *comp2 = (double *)&f2->psi[ix];
        double *fcomp2 = (double *)&force2[ix];

        for (int a = 0; a < 8; a++) {
            double orig = comp2[a];
            comp2[a] = orig + eps;
            Energy ep2 = field_energy(f2, &params2);
            comp2[a] = orig - eps;
            Energy em2 = field_energy(f2, &params2);
            comp2[a] = orig;

            double num = -(ep2.Etotal - em2.Etotal) / (2.0 * eps);
            double h3 = f2->h * f2->h * f2->h;
            double anal = h3 * fcomp2[a];

            double rel_err = (fabs(num) > 1e-12) ?
                fabs(anal - num) / fabs(num) : fabs(anal - num);

            const char *names[] = {"s ", "f1", "f2", "f3", "j1", "j2", "j3", "P "};
            printf("  %s  | %+.5e | %+.5e | %.3e %s\n",
                   names[a], anal, num, rel_err, rel_err > 0.05 ? " ***FAIL***" : "");

            if (rel_err > 0.05 && fabs(num) > 1e-10) pass = 0;
        }
        printf("\n");
    }

    Energy en2 = field_energy(f2, &params2);
    printf("Energy: E=%.6e (E2=%.4e E4=%.4e EV=%.4e ED=%.4e)\n",
           en2.Etotal, en2.E2, en2.E4, en2.EV, en2.ED);

    printf("\n%s\n", pass ? "=== ALL TESTS PASSED ===" : "=== SOME TESTS FAILED ===");

    free(force2);
    field_free(f2);
    return pass ? 0 : 1;
}
