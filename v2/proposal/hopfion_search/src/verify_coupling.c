/*
 * verify_coupling.c — Verify coupling gradient by comparing with finite differences.
 *
 * Tests all three coupling energy terms:
 *   E_{2,D}: degenerate gradient energy
 *   E_{4,C}: Skyrme cross-coupling energy
 *   E_int:   bulk-degenerate gradient coupling
 *
 * For each, perturbs a single field component and checks that
 * the analytical gradient matches the numerical gradient.
 */

#define _USE_MATH_DEFINES
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "field.h"
#include "coupling.h"
#include "initial.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Compute total energy including coupling */
static double total_coupling_energy(const Field *f, const Params *p, double g_coupling)
{
    CouplingEnergy ce = coupling_energy(f, p, g_coupling);
    return ce.Etotal;
}

int main(void)
{
    int N = 16;
    double L = 4.0;
    double eps = 1e-5;
    int pass = 1;

    /* === Test 1: E_{2,D} only (large e, g=0 suppresses E_{4,C} and E_int) === */
    printf("=== Test 1: E_{2,D} (degenerate gradient, e=1e6, g=0) ===\n");
    Params params1 = {1.0, 10.0, 1e6, 1.0, 0.0, 1.0};  /* rho0, lambda, e, mu, g, c */

    Field *f1 = field_alloc(N, L);
    init_hedgehog(f1, &params1, 1.5);
    init_perturb(f1, 0.05, 42);  /* Need non-trivial weight sector */

    int N3 = N * N * N;
    Multivector *force1 = (Multivector *)calloc((size_t)N3, sizeof(Multivector));

    /* Compute analytical coupling gradient (only adds to force, so zero first) */
    coupling_gradient(f1, &params1, 0.0, force1);

    int test_points[] = {idx(N, N/4, N/4, N/4), idx(N, N/2, N/2, N/2),
                         idx(N, 3*N/4, N/4, N/2)};
    int ntest = 3;

    printf("Comparing analytical coupling force with finite-difference d(-E_coupling)/dPsi\n\n");

    for (int t = 0; t < ntest; t++) {
        int ix = test_points[t];
        printf("Point %d (flat index %d):\n", t, ix);
        printf("  Comp | Analytical |   Numerical  | Rel Error\n");
        printf("  -----|------------|--------------|----------\n");

        double *comp = (double *)&f1->psi[ix];
        double *fcomp = (double *)&force1[ix];

        for (int a = 0; a < 8; a++) {
            double orig = comp[a];
            comp[a] = orig + eps;
            double ep = total_coupling_energy(f1, &params1, 0.0);
            comp[a] = orig - eps;
            double em = total_coupling_energy(f1, &params1, 0.0);
            comp[a] = orig;

            double num = -(ep - em) / (2.0 * eps);
            double h3 = f1->h * f1->h * f1->h;
            double anal = h3 * fcomp[a];

            double rel_err = (fabs(num) > 1e-12) ?
                fabs(anal - num) / fabs(num) : fabs(anal - num);

            const char *names[] = {"s ", "f1", "f2", "f3", "j1", "j2", "j3", "P "};
            printf("  %s  | %+.5e | %+.5e | %.3e %s\n",
                   names[a], anal, num, rel_err,
                   (rel_err > 0.05 && fabs(num) > 1e-8) ? " ***FAIL***" : "");

            if (rel_err > 0.05 && fabs(num) > 1e-8) pass = 0;
        }
        printf("\n");
    }

    CouplingEnergy ce1 = coupling_energy(f1, &params1, 0.0);
    printf("Coupling energy: E2D=%.4e, E4C=%.4e, Eint=%.4e\n\n",
           ce1.E2D, ce1.E4C, ce1.Eint);

    free(force1);
    field_free(f1);

    /* === Test 2: All coupling terms (e=2, g=1) === */
    printf("=== Test 2: Full coupling (e=2, g=1) ===\n");
    Params params2 = {1.0, 10.0, 2.0, 1.0, 1.0, 1.0};
    double g2 = 1.0;

    Field *f2 = field_alloc(N, L);
    init_hedgehog(f2, &params2, 1.5);
    init_perturb(f2, 0.05, 42);

    Multivector *force2 = (Multivector *)calloc((size_t)N3, sizeof(Multivector));
    coupling_gradient(f2, &params2, g2, force2);

    printf("Comparing analytical coupling force with finite-difference\n\n");

    for (int t = 0; t < ntest; t++) {
        int ix = test_points[t];
        printf("Point %d (flat index %d):\n", t, ix);
        printf("  Comp | Analytical |   Numerical  | Rel Error\n");
        printf("  -----|------------|--------------|----------\n");

        double *comp = (double *)&f2->psi[ix];
        double *fcomp = (double *)&force2[ix];

        for (int a = 0; a < 8; a++) {
            double orig = comp[a];
            comp[a] = orig + eps;
            double ep = total_coupling_energy(f2, &params2, g2);
            comp[a] = orig - eps;
            double em = total_coupling_energy(f2, &params2, g2);
            comp[a] = orig;

            double num = -(ep - em) / (2.0 * eps);
            double h3 = f2->h * f2->h * f2->h;
            double anal = h3 * fcomp[a];

            double rel_err = (fabs(num) > 1e-12) ?
                fabs(anal - num) / fabs(num) : fabs(anal - num);

            const char *names[] = {"s ", "f1", "f2", "f3", "j1", "j2", "j3", "P "};
            printf("  %s  | %+.5e | %+.5e | %.3e %s\n",
                   names[a], anal, num, rel_err,
                   (rel_err > 0.05 && fabs(num) > 1e-8) ? " ***FAIL***" : "");

            if (rel_err > 0.05 && fabs(num) > 1e-8) pass = 0;
        }
        printf("\n");
    }

    CouplingEnergy ce2 = coupling_energy(f2, &params2, g2);
    printf("Coupling energy: E2D=%.4e, E4C=%.4e, Eint=%.4e\n\n",
           ce2.E2D, ce2.E4C, ce2.Eint);

    /* === Test 3: Combined field.c + coupling.c gradient === */
    printf("=== Test 3: Combined gradient (field + coupling) ===\n");

    Multivector *force_combined = (Multivector *)calloc((size_t)N3, sizeof(Multivector));
    field_gradient(f2, &params2, force_combined);
    coupling_gradient(f2, &params2, g2, force_combined);

    printf("Comparing combined analytical force with finite-difference\n\n");

    for (int t = 0; t < ntest; t++) {
        int ix = test_points[t];
        printf("Point %d (flat index %d):\n", t, ix);
        printf("  Comp | Analytical |   Numerical  | Rel Error\n");
        printf("  -----|------------|--------------|----------\n");

        double *comp = (double *)&f2->psi[ix];
        double *fcomp = (double *)&force_combined[ix];

        for (int a = 0; a < 8; a++) {
            double orig = comp[a];

            comp[a] = orig + eps;
            Energy ep_f = field_energy(f2, &params2);
            double ep_c = total_coupling_energy(f2, &params2, g2);

            comp[a] = orig - eps;
            Energy em_f = field_energy(f2, &params2);
            double em_c = total_coupling_energy(f2, &params2, g2);

            comp[a] = orig;

            double num = -((ep_f.Etotal + ep_c) - (em_f.Etotal + em_c)) / (2.0 * eps);
            double h3 = f2->h * f2->h * f2->h;
            double anal = h3 * fcomp[a];

            double rel_err = (fabs(num) > 1e-12) ?
                fabs(anal - num) / fabs(num) : fabs(anal - num);

            const char *names[] = {"s ", "f1", "f2", "f3", "j1", "j2", "j3", "P "};
            printf("  %s  | %+.5e | %+.5e | %.3e %s\n",
                   names[a], anal, num, rel_err,
                   (rel_err > 0.05 && fabs(num) > 1e-8) ? " ***FAIL***" : "");

            if (rel_err > 0.05 && fabs(num) > 1e-8) pass = 0;
        }
        printf("\n");
    }

    printf("\n%s\n", pass ? "=== ALL TESTS PASSED ===" : "=== SOME TESTS FAILED ===");

    free(force_combined);
    free(force2);
    field_free(f2);
    return pass ? 0 : 1;
}
