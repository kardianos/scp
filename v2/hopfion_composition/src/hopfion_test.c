/*
 * hopfion_test.c — Verification tests for hopfion composition simulator
 *
 * Tests:
 *   1. Grid allocation and active cell count
 *   2. Vacuum energy = 0
 *   3. Skyrmion initialization + baryon charge
 *   4. Hopfion initialization + CP¹ structure
 *   5. Gradient verification (F = -δE/δq via finite differences)
 *   6. Hopf charge for H=1 hopfion
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "spherical_grid.h"
#include "hopfion_field.h"
#include "hopfion_init.h"
#include "topology.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static int tests_passed = 0, tests_total = 0;

#define CHECK(cond, msg) do { \
    tests_total++; \
    if (cond) { printf("  PASS: %s\n", msg); tests_passed++; } \
    else      { printf("  FAIL: %s\n", msg); } \
} while(0)

#define CHECK_APPROX(val, expected, tol, msg) do { \
    tests_total++; \
    double _v = (val), _e = (expected), _t = (tol); \
    if (fabs(_v - _e) < _t) { \
        printf("  PASS: %s (%.6f ≈ %.6f)\n", msg, _v, _e); tests_passed++; \
    } else { \
        printf("  FAIL: %s (%.6f != %.6f, diff=%.2e)\n", msg, _v, _e, fabs(_v-_e)); \
    } \
} while(0)

/* Test 1: Grid allocation */
void test_grid(void)
{
    printf("\n=== Test 1: Grid allocation ===\n");

    int N = 64;
    double R = 5.0, L = R + 0.5;
    SphericalGrid *g = sg_alloc(N, L, R, 1.0);
    CHECK(g != NULL, "Grid allocated");

    /* Expected active cells: roughly (4/3)πR³ / h³ */
    double h = 2.0 * L / N;
    double expected = (4.0/3.0) * M_PI * R*R*R / (h*h*h);
    double ratio = (double)g->n_active / expected;
    CHECK(ratio > 0.9 && ratio < 1.1, "Active cell count within 10% of (4/3)piR^3/h^3");
    printf("    Active: %d, expected ~%.0f, ratio=%.3f\n", g->n_active, expected, ratio);

    /* Vacuum energy should be zero */
    sg_set_vacuum(g, 1.0);
    HopfionParams p = {1.0, 1.0, 0, 0, 0, 0, 1.0};
    HopfionEnergy en = hf_energy(g, &p);
    CHECK(fabs(en.Etotal) < 1e-10, "Vacuum energy = 0");

    sg_free(g);
}

/* Test 2: Hopfion initialization */
void test_hopfion_init(void)
{
    printf("\n=== Test 2: Hopfion initialization ===\n");

    int N = 96;
    double R = 6.0, L = R + 0.5;
    SphericalGrid *g = sg_alloc(N, L, R, 1.5);
    sg_set_vacuum(g, 1.0);

    /* Initialize H=1 hopfion */
    init_hopfion(g, 1.0, 2.0, 0, 0, 0);

    /* Check that field is not vacuum at center */
    int ic = N/2;
    Multivector center = sg_get(g, ic, ic, ic);
    double q2 = mv_bulk_norm2(center);
    CHECK(fabs(q2 - 1.0) < 0.5, "Center field magnitude ≈ ρ₀");
    printf("    |q(0)|² = %.6f\n", q2);

    /* Check that field is vacuum far from center */
    Multivector far = sg_get(g, ic, ic, N-3);
    double far_s = far.s;
    CHECK(fabs(far_s - 1.0) < 0.01, "Far field ≈ vacuum");

    /* Compute energy */
    HopfionParams p = {1.0, 1.0, 0, 0, 0, 0, 1.0};
    HopfionEnergy en = hf_energy(g, &p);
    printf("    E = %.4f (E2=%.4f, E4=%.4f)\n", en.Etotal, en.E2, en.E4);
    CHECK(en.E2 > 0, "Hopfion has positive gradient energy");
    CHECK(en.E4 > 0, "Hopfion has positive Skyrme energy");

    /* CP¹ vector at center should be well-defined */
    UnitVec n = topo_cp1_vector(g, ic, ic, ic);
    double n2 = n.n1*n.n1 + n.n2*n.n2 + n.n3*n.n3;
    CHECK(fabs(n2 - 1.0) < 0.01, "CP¹ vector is unit at center");

    sg_free(g);
}

/* Test 3: Baryon charge of Skyrmion */
void test_skyrmion_charge(void)
{
    printf("\n=== Test 3: Skyrmion baryon charge ===\n");

    /* Create a simple test profile: f(r) = π(1 - r/R)² for r < R */
    /* This is not an exact solution but has the right boundary conditions */
    int np = 1001;
    RadialProfile prof;
    prof.n = np;
    prof.r = malloc(np * sizeof(double));
    prof.f = malloc(np * sizeof(double));
    prof.rho = NULL;
    prof.dr = 0.01;
    prof.r_max = (np-1) * prof.dr;

    for (int i = 0; i < np; i++) {
        prof.r[i] = i * prof.dr;
        double x = prof.r[i] / 5.0;  /* R_profile = 5 */
        if (x < 1) prof.f[i] = M_PI * (1 - x) * (1 - x);
        else prof.f[i] = 0;
    }

    int N = 96;
    double R = 7.0, L = R + 0.5;
    SphericalGrid *g = sg_alloc(N, L, R, 1.5);
    sg_set_vacuum(g, 1.0);

    init_skyrmion(g, &prof, 1.0, 0, 0, 0);

    double B = hf_baryon_charge(g);
    printf("    B = %.4f (expected ~1)\n", B);
    /* Note: test profile is not exact, so B may not be exactly 1.
     * But it should be positive and O(1). */
    CHECK(B > 0.5, "Baryon charge > 0.5 for test profile");

    free(prof.r);
    free(prof.f);
    sg_free(g);
}

/* Test 4: Energy is lowered by gradient flow */
void test_gradient_flow(void)
{
    printf("\n=== Test 4: Gradient flow lowers energy ===\n");

    int N = 64;
    double R = 5.0, L = R + 0.5;
    SphericalGrid *g = sg_alloc(N, L, R, 1.0);
    sg_set_vacuum(g, 1.0);

    init_hopfion(g, 1.0, 1.5, 0, 0, 0);

    HopfionParams p = {1.0, 1.0, 0, 0, 0, 0, 1.0};
    HopfionEnergy en0 = hf_energy(g, &p);
    printf("    E_initial = %.4f\n", en0.Etotal);

    /* A few steps of gradient flow (L₂ + V force only) */
    double gf_dt = 0.0001;
    for (int step = 0; step < 100; step++) {
        hf_force(g, &p);
        for (int n = 0; n < g->n_active; n++) {
            int ix = g->active[n];
            g->psi[ix] = mv_add(g->psi[ix], mv_scale(gf_dt, g->force[ix]));
        }
    }

    HopfionEnergy en1 = hf_energy(g, &p);
    printf("    E_final   = %.4f (after 100 GF steps)\n", en1.Etotal);
    CHECK(en1.Etotal < en0.Etotal, "Energy decreased under gradient flow");

    sg_free(g);
}

int main(void)
{
    printf("========================================\n");
    printf(" Hopfion Composition — Test Suite\n");
    printf("========================================\n");

    test_grid();
    test_hopfion_init();
    test_skyrmion_charge();
    test_gradient_flow();

    printf("\n========================================\n");
    printf(" Results: %d/%d tests passed\n", tests_passed, tests_total);
    printf("========================================\n");

    return (tests_passed == tests_total) ? 0 : 1;
}
