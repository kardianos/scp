/*
 * verify_gradient.c — Gradient consistency check for all force terms
 *
 * For each energy term (L₂, L₄, L₆, V, V_π, V_D), verifies that the
 * analytical force matches the numerical finite-difference gradient:
 *
 *   F_analytical × h³  ≈  -(E(q+εδ_a) - E(q-εδ_a)) / (2ε)
 *
 * Convention: hf_force stores -(1/h³) dE/dψ, so h³ × force = -dE/dψ.
 *
 * Tests:
 *   1. L₂ only (large e, no potential) — verifies consistent Laplacian
 *   2. L₂ + L₄ (Skyrme, e=2) — verifies analytical Skyrme force
 *   3. L₂ + L₄ + V (Mexican hat, λ=10) — verifies potential force
 *   4. L₂ + L₄ + V + V_π (pion mass) — verifies pion mass force
 *   5. L₂ + L₄ + L₆ (sextic) — verifies numerical L₆ force
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "spherical_grid.h"
#include "hopfion_field.h"
#include "hopfion_init.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Helper: perturb a hedgehog-like field so it's not on any special point */
static void perturb_field(SphericalGrid *g, double amp, int seed)
{
    srand(seed);
    for (int n = 0; n < g->n_active; n++) {
        int ix = g->active[n];
        Multivector *psi = &g->psi[ix];
        psi->s  += amp * (2.0*rand()/RAND_MAX - 1.0);
        psi->f1 += amp * (2.0*rand()/RAND_MAX - 1.0);
        psi->f2 += amp * (2.0*rand()/RAND_MAX - 1.0);
        psi->f3 += amp * (2.0*rand()/RAND_MAX - 1.0);
        psi->j1 += amp * (2.0*rand()/RAND_MAX - 1.0);
        psi->j2 += amp * (2.0*rand()/RAND_MAX - 1.0);
        psi->j3 += amp * (2.0*rand()/RAND_MAX - 1.0);
        psi->p  += amp * (2.0*rand()/RAND_MAX - 1.0);
    }
}

/* Run gradient verification on a given configuration.
 * Returns number of failures. */
static int verify_gradient(SphericalGrid *g, const HopfionParams *p,
                           const char *test_name, double tol)
{
    int N = g->N;
    double h3 = g->h * g->h * g->h;
    double eps = 1e-5;

    /* Compute analytical force */
    hf_force(g, p);

    /* Select test points: pick cells that are well inside the sphere.
     * The consistent Laplacian reads ±4 neighbors, and the energy sums
     * over cells within R-4.5h. For the Laplacian to be the exact gradient
     * of the truncated energy, ALL ±4 neighbors of the test cell must be
     * within the energy region. This requires r < R - 4.5h - 4h = R - 8.5h.
     * We use R - 9h for safety margin. */
    double r_safe = g->R - 9.0 * g->h;
    int ic = N / 2;
    int test_ijk[][3] = {
        {ic - N/8, ic,       ic},        /* offset from center */
        {ic,       ic + N/10, ic - N/10}, /* diagonal offset */
        {ic + N/8, ic - N/10, ic + N/12}, /* another offset */
    };
    int ntest = 3;

    printf("\n=== %s ===\n", test_name);
    printf("  eps = %.0e, h = %.4f, r_safe = %.3f\n\n", eps, g->h, r_safe);

    int failures = 0;
    int total_checks = 0;

    for (int t = 0; t < ntest; t++) {
        int i = test_ijk[t][0], j = test_ijk[t][1], k = test_ijk[t][2];
        int ix = sg_idx(N, i, j, k);

        /* Make sure this cell is active and well inside the consistency zone */
        if (!sg_is_active(g, i, j, k)) {
            printf("  Point %d (%d,%d,%d): SKIP (not active)\n\n", t, i, j, k);
            continue;
        }
        double r = sg_radius(g, i, j, k);
        if (r > r_safe) {
            printf("  Point %d (%d,%d,%d): SKIP (r=%.2f > r_safe=%.2f)\n\n",
                   t, i, j, k, r, r_safe);
            continue;
        }

        printf("  Point %d (%d,%d,%d), r=%.3f:\n", t, i, j, k, r);
        printf("    Comp | Analytical    |   Numerical   | Rel Error\n");
        printf("    -----|---------------|---------------|----------\n");

        double *comp = (double *)&g->psi[ix];
        double *fcomp = (double *)&g->force[ix];

        /* Test all 8 components (4 bulk + 4 degenerate) */
        int ncomp = 8;
        const char *names[] = {"s ", "f1", "f2", "f3", "j1", "j2", "j3", "P "};

        for (int a = 0; a < ncomp; a++) {
            double orig = comp[a];

            /* E(+eps) */
            comp[a] = orig + eps;
            HopfionEnergy ep = hf_energy(g, p);

            /* E(-eps) */
            comp[a] = orig - eps;
            HopfionEnergy em = hf_energy(g, p);

            /* Restore */
            comp[a] = orig;

            /* Numerical: -(E+ - E-)/(2ε) */
            double num = -(ep.Etotal - em.Etotal) / (2.0 * eps);

            /* Analytical: h³ × force */
            double anal = h3 * fcomp[a];

            /* Relative error */
            double denom = fabs(num) > 1e-12 ? fabs(num) : 1e-12;
            double rel_err = fabs(anal - num) / denom;

            /* For very small values, use absolute error */
            int is_small = (fabs(num) < 1e-10 && fabs(anal) < 1e-10);
            const char *status = "";
            if (!is_small && rel_err > tol) {
                status = " ***FAIL***";
                failures++;
            }
            total_checks++;

            printf("    %s  | %+.6e | %+.6e | %.3e%s\n",
                   names[a], anal, num, rel_err, status);
        }
        printf("\n");
    }

    printf("  Result: %d/%d checks passed",
           total_checks - failures, total_checks);
    if (failures > 0)
        printf(" (%d FAILURES)", failures);
    printf("\n");

    return failures;
}

/* Create a test Skyrmion profile: f(r) = π(1-r/R_prof)² for r < R_prof */
static RadialProfile make_test_profile(double R_prof)
{
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
        double x = prof.r[i] / R_prof;
        prof.f[i] = (x < 1.0) ? M_PI * (1 - x) * (1 - x) : 0.0;
    }
    return prof;
}

int main(void)
{
    printf("========================================\n");
    printf(" Gradient Consistency Verification\n");
    printf("========================================\n");

    /* Use a small grid for speed, but big enough for the stencils */
    int N = 48;
    double R = 5.0, L = R + 0.5;
    double sponge_w = 0.5;

    int total_failures = 0;

    /* --- Test 1: L₂ only (very large e kills L₄, no potential) --- */
    {
        SphericalGrid *g = sg_alloc(N, L, R, sponge_w);
        sg_set_vacuum(g, 1.0);

        RadialProfile prof = make_test_profile(3.5);
        init_skyrmion(g, &prof, 1.0, 0, 0, 0);
        perturb_field(g, 0.02, 42);

        HopfionParams p = {1.0, 1e6, 0, 0, 0, 0, 1.0};
        /* rho0=1, e=1e6 (kills L₄), lambda=0, lambda6=0, m_pi_sq=0, mu=0, c=1 */

        total_failures += verify_gradient(g, &p, "Test 1: L₂ only (e=1e6)", 0.01);

        free(prof.r); free(prof.f);
        sg_free(g);
    }

    /* --- Test 2: L₂ + L₄ (Skyrme force, e=2) --- */
    {
        SphericalGrid *g = sg_alloc(N, L, R, sponge_w);
        sg_set_vacuum(g, 1.0);

        RadialProfile prof = make_test_profile(3.5);
        init_skyrmion(g, &prof, 1.0, 0, 0, 0);
        perturb_field(g, 0.02, 42);

        HopfionParams p = {1.0, 2.0, 0, 0, 0, 0, 1.0};
        /* rho0=1, e=2, lambda=0, lambda6=0, m_pi_sq=0, mu=0, c=1 */

        total_failures += verify_gradient(g, &p, "Test 2: L₂ + L₄ (e=2.0)", 0.01);

        free(prof.r); free(prof.f);
        sg_free(g);
    }

    /* --- Test 3: L₂ + L₄ + V (Mexican hat, λ=10) --- */
    {
        SphericalGrid *g = sg_alloc(N, L, R, sponge_w);
        sg_set_vacuum(g, 1.0);

        RadialProfile prof = make_test_profile(3.5);
        init_skyrmion(g, &prof, 1.0, 0, 0, 0);
        perturb_field(g, 0.02, 42);

        HopfionParams p = {1.0, 2.0, 10.0, 0, 0, 0, 1.0};
        /* rho0=1, e=2, lambda=10, lambda6=0, m_pi_sq=0, mu=0, c=1 */

        total_failures += verify_gradient(g, &p, "Test 3: L₂ + L₄ + V (λ=10)", 0.01);

        free(prof.r); free(prof.f);
        sg_free(g);
    }

    /* --- Test 4: L₂ + L₄ + V + V_π (pion mass) --- */
    {
        SphericalGrid *g = sg_alloc(N, L, R, sponge_w);
        sg_set_vacuum(g, 1.0);

        RadialProfile prof = make_test_profile(3.5);
        init_skyrmion(g, &prof, 1.0, 0, 0, 0);
        perturb_field(g, 0.02, 42);

        HopfionParams p = {1.0, 2.0, 10.0, 0, 0.16, 0, 1.0};
        /* rho0=1, e=2, lambda=10, lambda6=0, m_pi_sq=0.16, mu=0, c=1 */

        total_failures += verify_gradient(g, &p, "Test 4: L₂ + L₄ + V + Vπ (m²π=0.16)", 0.01);

        free(prof.r); free(prof.f);
        sg_free(g);
    }

    /* --- Test 5: Full (L₂ + L₄ + V + V_π + V_D) --- */
    {
        SphericalGrid *g = sg_alloc(N, L, R, sponge_w);
        sg_set_vacuum(g, 1.0);

        RadialProfile prof = make_test_profile(3.5);
        init_skyrmion(g, &prof, 1.0, 0, 0, 0);
        perturb_field(g, 0.02, 42);

        HopfionParams p = {1.0, 2.0, 10.0, 0, 0.16, 2.0, 1.0};
        /* rho0=1, e=2, lambda=10, lambda6=0, m_pi_sq=0.16, mu=2.0, c=1 */

        total_failures += verify_gradient(g, &p, "Test 5: Full (L₂+L₄+V+Vπ+VD, μ=2)", 0.01);

        free(prof.r); free(prof.f);
        sg_free(g);
    }

    /* --- Test 6: L₂ + L₄ + L₆ (sextic, small grid) --- */
    /* L₆ force is expensive (numerical FD), so use small grid.
     * Need R large enough that R-9h is beyond soliton core. */
    {
        int N6 = 40;
        double R6 = 4.5, L6 = R6 + 0.5;
        SphericalGrid *g = sg_alloc(N6, L6, R6, sponge_w);
        sg_set_vacuum(g, 1.0);

        RadialProfile prof = make_test_profile(2.5);
        init_skyrmion(g, &prof, 1.0, 0, 0, 0);
        perturb_field(g, 0.02, 42);

        HopfionParams p = {1.0, 2.0, 0, 0.5, 0, 0, 1.0};
        /* rho0=1, e=2, lambda=0, lambda6=0.5, m_pi_sq=0, mu=0, c=1 */

        printf("\n  (L₆ test uses smaller grid N=%d — force is O(N³×125) per component)\n", N6);

                /* L₆ force is itself a numerical FD, so expect ~10⁻³ accuracy.
         * Use 2% tolerance (allowing ~1% from force FD + ~1% from verification FD). */
        total_failures += verify_gradient(g, &p, "Test 6: L₂ + L₄ + L₆ (λ₆=0.5)", 0.02);

        free(prof.r); free(prof.f);
        sg_free(g);
    }

    printf("\n========================================\n");
    if (total_failures == 0)
        printf(" ALL GRADIENT CHECKS PASSED\n");
    else
        printf(" %d GRADIENT CHECKS FAILED\n", total_failures);
    printf("========================================\n");

    return (total_failures == 0) ? 0 : 1;
}
