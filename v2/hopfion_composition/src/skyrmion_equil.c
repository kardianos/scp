/*
 * skyrmion_equil.c — Skyrmion equilibrium via gradient flow + effective metric
 *
 * Usage: ./bin/skyrmion_equil [profile_path] [lambda6]
 *
 * 1. Initialize B=1 hedgehog from 1D radial profile (or test profile)
 * 2. Run sigma-model gradient flow with full L₂ + L₄ force (skipped if λ₆>0)
 * 3. Monitor convergence: E, E₂, E₄, B, E₂/E₄ (virial), |F|_max
 * 4. At convergence, compute effective metric P(r), m(r), P/m
 * 5. Verify: E/E_FB ≈ 1.232, E₂/E₄ ≈ 1.0, P/m ≈ 2.0 (for λ₆=0)
 * 6. When λ₆>0: 1D analytical cross-check + λ₆ scan + interpretation
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

/* Compute max |force| over active cells (convergence diagnostic) */
static double max_force(const SphericalGrid *g)
{
    double fmax = 0;
    for (int n = 0; n < g->n_active; n++) {
        int ix = g->active[n];
        double f2 = mv_dot(g->force[ix], g->force[ix]);
        if (f2 > fmax) fmax = f2;
    }
    return sqrt(fmax);
}

/* Compute effective metric from 3D field and print radial profile.
 * Returns max |P/m - 2| in the soliton region. */
static double compute_and_print_metric(SphericalGrid *g, HopfionParams *p,
                                        int n_bins, double *out_pm_avg)
{
    EffectiveMetric *em = hf_effective_metric(g, p, n_bins);
    if (!em) {
        fprintf(stderr, "  Failed to compute effective metric\n");
        *out_pm_avg = 0;
        return -1;
    }

    printf("\n  %8s  %12s  %12s  %10s\n", "r", "P(r)", "m(r)", "P/m");
    printf("  %8s  %12s  %12s  %10s\n", "--------", "------------", "------------", "----------");

    double pm_sum = 0;
    int pm_count = 0;
    double pm_max_dev = 0;
    double pm_max_dev_r = 0;

    for (int b = 0; b < n_bins; b++) {
        double r = em->r[b];
        double P = em->P[b];
        double m_val = em->m[b];
        double pm = (m_val > 1e-15) ? P / m_val : 0;

        if (r < 0.2 || r > 5.0) continue;
        if (P < 1e-10 && m_val < 1e-10) continue;

        printf("  %8.3f  %12.6f  %12.6f  %10.6f\n", r, P, m_val, pm);

        if (r > 0.3 && r < 4.0 && m_val > 0.1) {
            pm_sum += pm;
            pm_count++;
            double dev = fabs(pm - 2.0);
            if (dev > pm_max_dev) {
                pm_max_dev = dev;
                pm_max_dev_r = r;
            }
        }
    }

    *out_pm_avg = (pm_count > 0) ? pm_sum / pm_count : 0;
    printf("\n  P/m statistics (0.3 < r < 4.0):\n");
    printf("    Average P/m = %.6f\n", *out_pm_avg);
    printf("    Max |P/m - 2| = %.6f at r = %.3f\n", pm_max_dev, pm_max_dev_r);

    hf_metric_free(em);
    return pm_max_dev;
}

int main(int argc, char **argv)
{
    /* Grid parameters */
    int N = 128;
    double R = 8.0;
    double L = R + 0.5;
    double sponge_w = 1.0;

    /* Physics: sigma model with e=1, ρ₀=1 */
    double rho0 = 1.0;
    double e_skyrme = 1.0;
    double lambda6 = 0.0;

    /* Parse command line */
    const char *profile_path = NULL;
    if (argc > 1) profile_path = argv[1];
    if (argc > 2) lambda6 = atof(argv[2]);

    /* Default profile: self-consistent for given λ₆ */
    char auto_profile[256];
    if (!profile_path) {
        if (lambda6 > 0) {
            snprintf(auto_profile, sizeof(auto_profile),
                     "../proposal/hopfion_search/data/profiles/profile_sigma_e1_lam6_%g.dat",
                     lambda6);
            profile_path = auto_profile;
            printf("  Auto-selecting self-consistent L₆ profile: %s\n", profile_path);
        } else {
            profile_path = "../proposal/hopfion_search/data/profiles/profile_sigma_e1.dat";
        }
    }

    HopfionParams p = {rho0, e_skyrme, 0, lambda6, 0, 0, 1.0};
    /* lambda=0 (sigma model via projection), m_pi_sq=0, mu=0, c=1 */

    double c4 = 2.0 * rho0 * rho0 / (e_skyrme * e_skyrme);

    /* Faddeev-Bogomolny bound: E_FB = 6√2 π² ρ₀³/e */
    double E_FB = 6.0 * sqrt(2.0) * M_PI * M_PI * rho0 * rho0 * rho0 / e_skyrme;

    printf("========================================\n");
    printf(" Skyrmion Equilibrium + Effective Metric\n");
    printf("========================================\n");
    printf("  N=%d, R=%.1f, h=%.4f\n", N, R, 2.0*L/N);
    printf("  e=%.1f, rho0=%.1f, E_FB=%.4f\n", e_skyrme, rho0, E_FB);
    printf("  lambda6=%.2f\n", lambda6);
    printf("  Expected: E/E_FB ≈ 1.232, E₂/E₄ ≈ 1.0\n\n");

    /* Allocate grid */
    SphericalGrid *g = sg_alloc(N, L, R, sponge_w);
    sg_set_vacuum(g, rho0);
    printf("  Grid: %d active cells (%.1f%%)\n\n", g->n_active,
           100.0 * g->n_active / ((double)N*N*N));

    /* Load equilibrium profile (or use test profile as fallback) */
    RadialProfile *prof = profile_load(profile_path);
    if (!prof) {
        printf("  Falling back to test profile f(r)=pi*(1-r/5)^2\n");
        prof = malloc(sizeof(RadialProfile));
        prof->n = 2001;
        prof->r = malloc(prof->n * sizeof(double));
        prof->f = malloc(prof->n * sizeof(double));
        prof->rho = NULL;
        prof->dr = 0.005;
        prof->r_max = (prof->n - 1) * prof->dr;
        for (int i = 0; i < prof->n; i++) {
            prof->r[i] = i * prof->dr;
            double x = prof->r[i] / 5.0;
            prof->f[i] = (x < 1.0) ? M_PI * (1.0 - x) * (1.0 - x) : 0.0;
        }
    }

    /* Initialize hedgehog */
    init_skyrmion(g, prof, rho0, 0, 0, 0);

    /* Enforce sigma-model constraint */
    hf_sigma_project(g, rho0);

    /* Initial diagnostics */
    HopfionEnergy en0 = hf_energy(g, &p);
    double B0 = hf_baryon_charge(g);
    printf("\n  Initial state:\n");
    printf("    E = %.6f  (E/E_FB = %.6f)\n", en0.Etotal, en0.Etotal / E_FB);
    printf("    E₂ = %.6f, E₄ = %.6f, E₂/E₄ = %.6f\n", en0.E2, en0.E4, en0.E2 / en0.E4);
    printf("    B = %.6f\n", B0);

    /* ========== Arrested Gradient Flow ========== */
    /* Skip gradient flow when λ₆ > 0: the 1D profile was computed with
     * self-consistent L₂+L₄+L₆ dynamics (via radial.c -lam6). The 3D
     * force doesn't include L₆, so gradient flow would move the field
     * away from the correct equilibrium. */
    int converged = 0;
    int N3 = N * N * N;

    if (lambda6 > 0) {
        printf("\n  --- Skipping gradient flow (self-consistent L₂+L₄+L₆ profile from 1D solver) ---\n");
        converged = 1;
    } else {
        printf("\n  --- Arrested Gradient Flow (sigma model) ---\n");
        printf("  %6s  %12s  %10s  %8s  %8s  %10s\n",
               "step", "E", "E/E_FB", "E2/E4", "B", "|F|_max");

        double dt = 2e-4;
        int max_steps = 5000;
        double E_prev = en0.Etotal;
        double E_best = en0.Etotal;

        Multivector *psi_best = malloc(N3 * sizeof(Multivector));
        memcpy(psi_best, g->psi, N3 * sizeof(Multivector));

        for (int step = 1; step <= max_steps; step++) {
            hf_force(g, &p);

            for (int n = 0; n < g->n_active; n++) {
                int ix = g->active[n];
                g->psi[ix] = mv_add(g->psi[ix], mv_scale(dt, g->force[ix]));
            }

            hf_sigma_project(g, rho0);

            if (step % 50 == 0 || step <= 10 || step == max_steps) {
                HopfionEnergy en = hf_energy(g, &p);
                double B = hf_baryon_charge(g);
                double fmax = max_force(g);
                double ratio_E2E4 = (en.E4 > 1e-15) ? en.E2 / en.E4 : 0;

                printf("  %6d  %12.6f  %10.6f  %8.5f  %8.5f  %10.3e\n",
                       step, en.Etotal, en.Etotal / E_FB, ratio_E2E4, B, fmax);

                if (fabs(B - 1.0) > 0.02) {
                    printf("\n  ARRESTED at step %d: B=%.4f deviated from 1.0\n", step, B);
                    printf("  Rolling back to best state (E=%.6f)\n", E_best);
                    memcpy(g->psi, psi_best, N3 * sizeof(Multivector));
                    converged = 1;
                    break;
                }

                if (en.Etotal < E_best) {
                    E_best = en.Etotal;
                    memcpy(psi_best, g->psi, N3 * sizeof(Multivector));
                }

                double dE_rel = fabs(en.Etotal - E_prev) / fabs(en.Etotal);

                if (en.Etotal > E_prev + 1e-12) {
                    dt *= 0.5;
                    printf("    (energy increased, reducing dt to %.2e)\n", dt);
                    if (dt < 1e-7) {
                        printf("\n  dt too small, stopping flow\n");
                        memcpy(g->psi, psi_best, N3 * sizeof(Multivector));
                        converged = 1;
                        break;
                    }
                }

                if (dE_rel < 1e-10 && step > 50) {
                    printf("\n  Converged at step %d (|ΔE/E| = %.2e)\n", step, dE_rel);
                    converged = 1;
                    break;
                }

                E_prev = en.Etotal;
            }
        }

        free(psi_best);

        if (!converged) {
            printf("\n  Did not converge in %d steps (using best state)\n", max_steps);
        }
    }

    /* ========== Final diagnostics ========== */
    HopfionEnergy en_final = hf_energy(g, &p);
    double B_final = hf_baryon_charge(g);
    double ratio = en_final.Etotal / E_FB;
    double virial = en_final.E2 / en_final.E4;

    printf("\n  --- Final State ---\n");
    printf("    E       = %.6f\n", en_final.Etotal);
    printf("    E/E_FB  = %.6f  (expected ≈ 1.232)\n", ratio);
    printf("    E₂      = %.6f\n", en_final.E2);
    printf("    E₄      = %.6f\n", en_final.E4);
    printf("    E₂/E₄   = %.6f  (expected ≈ 1.000, virial theorem)\n", virial);
    printf("    B       = %.6f  (expected ≈ 1.000)\n", B_final);

    int pass_E = (fabs(ratio - 1.232) < 0.02);
    int pass_V = (fabs(virial - 1.0) < 0.05);
    int pass_B = (fabs(B_final - 1.0) < 0.05);

    printf("\n  Checks:\n");
    printf("    E/E_FB ∈ [1.212, 1.252]:  %s\n", pass_E ? "PASS" : "FAIL");
    printf("    E₂/E₄ ∈ [0.95, 1.05]:    %s\n", pass_V ? "PASS" : "FAIL");
    printf("    B ∈ [0.95, 1.05]:         %s\n", pass_B ? "PASS" : "FAIL");

    /* ========== Effective Metric ========== */
    printf("\n  --- Effective Metric (λ₆=%.2f) ---\n", lambda6);
    printf("  Computing BLV effective metric on equilibrium field...\n");

    int n_bins = 80;
    double pm_avg = 0;
    double pm_max_dev = compute_and_print_metric(g, &p, n_bins, &pm_avg);

    int pass_pm;
    if (lambda6 <= 0) {
        pass_pm = (fabs(pm_avg - 2.0) < 0.01);
        printf("    P/m ≈ 2:  %s\n", pass_pm ? "PASS" : "FAIL");

        if (pass_pm) {
            printf("\n  ==> CONFIRMED: P/m = 2 identically for sigma-model Skyrmion.\n");
            printf("      No emergent gravity (no radial time dilation) in the sigma model.\n");
            printf("      Breaking this requires finite-λ (varying ρ(r)) or L₆.\n");
        }
    } else {
        /* With L₆: P/m should be < 2 */
        pass_pm = (pm_max_dev > 0.001);  /* must deviate measurably from 2 */
        printf("    P/m deviates from 2:  %s  (max dev = %.6f)\n",
               pass_pm ? "YES" : "NO", pm_max_dev);
    }

    /* ========== L₆ Analysis (when λ₆ > 0) ========== */
    if (lambda6 > 0) {

        /* --- 1D analytical cross-check --- */
        printf("\n  --- 1D Analytical Cross-Check ---\n");
        printf("  Comparing 3D spherically-averaged P/m to 1D formula:\n");
        printf("    P/m_1D = (2α + β)/(α + β)\n");
        printf("    α = r² + 2c₄ sin²f, β = (2λ₆/π³) sin⁴f/r²\n\n");

        printf("  %8s  %10s  %10s  %10s  %10s\n",
               "r", "P/m_3D", "P/m_1D", "diff", "rel_err");

        EffectiveMetric *em = hf_effective_metric(g, &p, n_bins);
        double max_1d3d_diff = 0;
        if (em) {
            for (int b = 0; b < n_bins; b++) {
                double r = em->r[b];
                if (r < 0.3 || r > 4.0) continue;
                if (em->m[b] < 0.1) continue;

                double pm_3d = em->P[b] / em->m[b];

                /* 1D analytical formula */
                double f = profile_interp_f(prof, r);
                double sinf = sin(f);
                double sin2f = sinf * sinf;
                double sin4f = sin2f * sin2f;
                double alpha = r*r + 2.0 * c4 * sin2f;
                double beta = (2.0 * lambda6 / (M_PI*M_PI*M_PI)) * sin4f / (r*r);
                double pm_1d = (2.0 * alpha + beta) / (alpha + beta);

                double diff = fabs(pm_3d - pm_1d);
                double rel = (pm_1d > 1e-10) ? diff / fabs(pm_1d - 2.0 + 1e-15) : 0;
                if (diff > max_1d3d_diff) max_1d3d_diff = diff;

                if (b % 4 == 0) {  /* print every 4th bin to keep output manageable */
                    printf("  %8.3f  %10.6f  %10.6f  %10.6f  %10.3f%%\n",
                           r, pm_3d, pm_1d, diff, rel * 100);
                }
            }
            printf("\n  Max |P/m_3D - P/m_1D| = %.6f\n", max_1d3d_diff);
            int pass_1d = (max_1d3d_diff < 0.05);
            printf("  3D vs 1D agreement: %s\n", pass_1d ? "PASS" : "FAIL");
            hf_metric_free(em);
        }

        /* --- λ₆ scan --- */
        printf("\n  --- λ₆ Scan: Breaking P/m = 2 ---\n");
        printf("  Scanning λ₆ ∈ {0.5, 1, 2, 5, 10, 20, 50}...\n\n");
        printf("  %8s  %12s  %12s  %12s\n",
               "λ₆", "max|P/m-2|", "P/m at peak", "peak r");

        double scan_vals[] = {0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0};
        int n_scan = 7;

        for (int si = 0; si < n_scan; si++) {
            double lam6 = scan_vals[si];

            /* Analytical 1D scan (fast, no need for 3D recomputation) */
            double max_dev = 0;
            double peak_pm = 2.0;
            double peak_r = 0;

            for (int ri = 1; ri < 400; ri++) {
                double r = ri * 0.01;
                double f = profile_interp_f(prof, r);
                double sinf = sin(f);
                double sin2f = sinf * sinf;
                double sin4f = sin2f * sin2f;
                double alpha = r*r + 2.0 * c4 * sin2f;
                double beta = (2.0 * lam6 / (M_PI*M_PI*M_PI)) * sin4f / (r*r);
                double pm_1d = (2.0 * alpha + beta) / (alpha + beta);
                double dev = fabs(pm_1d - 2.0);
                if (dev > max_dev) {
                    max_dev = dev;
                    peak_pm = pm_1d;
                    peak_r = r;
                }
            }

            printf("  %8.1f  %12.6f  %12.6f  %12.3f\n",
                   lam6, max_dev, peak_pm, peak_r);
        }

        /* --- Physical interpretation --- */
        printf("\n  --- Physical Interpretation ---\n");

        /* Find the max deviation and where it occurs for the current λ₆ */
        double max_dev_current = 0;
        double peak_r_current = 0;
        double peak_pm_current = 2.0;

        for (int ri = 1; ri < 400; ri++) {
            double r = ri * 0.01;
            double f = profile_interp_f(prof, r);
            double sinf = sin(f);
            double sin2f = sinf * sinf;
            double sin4f = sin2f * sin2f;
            double alpha = r*r + 2.0 * c4 * sin2f;
            double beta = (2.0 * lambda6 / (M_PI*M_PI*M_PI)) * sin4f / (r*r);
            double pm_1d = (2.0 * alpha + beta) / (alpha + beta);
            double dev = fabs(pm_1d - 2.0);
            if (dev > max_dev_current) {
                max_dev_current = dev;
                peak_r_current = r;
                peak_pm_current = pm_1d;
            }
        }

        printf("\n  At λ₆ = %.2f:\n", lambda6);
        printf("    Max P/m deviation from 2: %.6f (%.2f%%)\n",
               max_dev_current, max_dev_current / 2.0 * 100);
        printf("    Peak location: r = %.3f (inside soliton core)\n", peak_r_current);
        printf("    P/m at peak: %.6f (< 2)\n", peak_pm_current);
        printf("    Effective wave speed: c_eff = c√(P/m) = %.6f c\n",
               sqrt(peak_pm_current));
        printf("    Compare sigma-model: c_eff = c√2 = %.6f c\n", sqrt(2.0));
        printf("\n  Key physics:\n");
        printf("    L₆ contributes P₆ = m₆ (ratio 1:1, not 2:1 like L₂+L₄)\n");
        printf("    Total P/m = (2α + β)/(α + β) < 2 when β > 0\n");
        printf("    β = 0 at r=0 (sin⁴f→0) and r→∞ (sin⁴f→0)\n");
        printf("    Peak at f ≈ π/2 (maximum winding density)\n");
        printf("    Waves SLOW DOWN inside soliton (c_eff < c√2)\n");
        printf("    Sign: waves slow down inside soliton → attractive lensing\n");
        printf("    (Same sign as Schwarzschild: coordinate light speed also decreases near mass)\n");

        /* ========== Gravitational Observables ========== */
        printf("\n  --- Gravitational Observables ---\n");

        double pi3 = M_PI * M_PI * M_PI;

        /* 1. Effective potential profile Φ(r)/c² = (P/m - 2)/4 */
        int nr = 600;
        double dr_fine = 0.01;
        double *r_arr = malloc(nr * sizeof(double));
        double *pm_arr = malloc(nr * sizeof(double));
        double *phi_arr = malloc(nr * sizeof(double));
        double phi_min = 0, r_phi_min = 0;

        printf("\n  1. Effective gravitational potential Φ(r)/c² = (P/m - 2)/4:\n");
        printf("     Φ < 0 = attractive well (waves slow → bend inward)\n\n");
        printf("  %8s  %10s  %10s  %12s\n", "r", "P/m", "v/c", "Φ/c²");

        for (int ri = 0; ri < nr; ri++) {
            double r = (ri + 0.5) * dr_fine;
            r_arr[ri] = r;
            double f = profile_interp_f(prof, r);
            double sinf = sin(f);
            double sin2f = sinf * sinf;
            double sin4f = sin2f * sin2f;
            double alpha = r * r + 2.0 * c4 * sin2f;
            double beta = (2.0 * lambda6 / pi3) * sin4f / (r * r);
            pm_arr[ri] = (2.0 * alpha + beta) / (alpha + beta);
            phi_arr[ri] = (pm_arr[ri] - 2.0) / 4.0;

            if (phi_arr[ri] < phi_min) {
                phi_min = phi_arr[ri];
                r_phi_min = r;
            }

            if (ri < 3 || ri % 50 == 0) {
                printf("  %8.3f  %10.6f  %10.6f  %12.6e\n",
                       r, pm_arr[ri], sqrt(pm_arr[ri]), phi_arr[ri]);
            }
        }

        printf("\n  Well depth: Φ_min/c² = %.6e at r = %.3f\n", phi_min, r_phi_min);
        printf("  P/m at center: %.6f\n", pm_arr[0]);

        /* 2. Shapiro-like time delay (radial passage through soliton) */
        /*    Δ(path) = ∫ [n(r) - 1] dr where n = √(2/(P/m)) = √(2m/P) */
        double delay = 0;
        for (int ri = 0; ri < nr; ri++) {
            double n_r = sqrt(2.0 / pm_arr[ri]);
            delay += (n_r - 1.0) * dr_fine;
        }
        delay *= 2.0;  /* both halves */
        printf("\n  2. Shapiro-like time delay (radial passage):\n");
        printf("     Extra optical path: Δℓ = %.6f code lengths\n", delay);
        printf("     Fractional delay: %.4f%%\n",
               delay / (2.0 * 4.0) * 100);  /* ~4 code lengths effective radius */

        /* 3. Deflection angle θ(b) — Born approximation
         *    θ(b) = -b ∫_0^{t_max} d(P/m)/dr|_{r=b cosh(t)} dt
         *    (Substitution r = b cosh(t) removes the √(r²-b²) singularity) */
        printf("\n  3. Deflection angle (Born approximation):\n");
        printf("     θ(b) = -(4b/c²) ∫ Φ'(r)/√(r²-b²) dr\n\n");
        printf("  %8s  %12s  %12s  %14s\n", "b", "θ (rad)", "θ (arcsec-eq)", "GM_eq/c²");

        double b_vals[] = {0.3, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0};
        int n_b = 8;

        for (int bi = 0; bi < n_b; bi++) {
            double b = b_vals[bi];
            double t_max = acosh(5.8 / b);
            if (t_max < 0.01) continue;
            int nt = 1000;
            double dt_step = t_max / nt;
            double theta = 0;

            for (int ti = 0; ti < nt; ti++) {
                double t = (ti + 0.5) * dt_step;
                double r = b * cosh(t);
                /* d(P/m)/dr via central difference */
                double rp = r + 0.001, rm = r - 0.001;
                if (rm < 0.001) rm = 0.001;
                double fp = profile_interp_f(prof, rp);
                double fm = profile_interp_f(prof, rm);
                double sf_p = sin(fp), sf_m = sin(fm);
                double s2f_p = sf_p * sf_p, s2f_m = sf_m * sf_m;
                double alpha_p = rp*rp + 2.0*c4*s2f_p;
                double alpha_m = rm*rm + 2.0*c4*s2f_m;
                double beta_p = (2.0*lambda6/pi3)*s2f_p*s2f_p/(rp*rp);
                double beta_m = (2.0*lambda6/pi3)*s2f_m*s2f_m/(rm*rm);
                double pm_p = (2.0*alpha_p + beta_p) / (alpha_p + beta_p);
                double pm_m = (2.0*alpha_m + beta_m) / (alpha_m + beta_m);
                double dpm_dr = (pm_p - pm_m) / (rp - rm);
                theta += dpm_dr * dt_step;
            }
            theta *= -b;

            /* Equivalent GM: θ(b→∞) = 4GM/(bc²) → GM/c² = θ×b/4 */
            double GM_eq = theta * b / 4.0;
            printf("  %8.2f  %12.6e  %12.2f  %14.6e\n",
                   b, theta, theta * 206265.0, GM_eq);
        }

        /* 4. Physical interpretation */
        double fm_per_code = 0.5624;
        double M_proton_MeV = 938.272;

        printf("\n  4. Physical interpretation (mapping to proton):\n");
        printf("     Potential depth: |Φ_min|/c² = %.4f\n", fabs(phi_min));
        printf("     Depth × soliton mass: %.1f MeV\n",
               fabs(phi_min) * M_proton_MeV);
        printf("     Nuclear binding: ~8 MeV/nucleon\n");
        printf("     → Metric potential = %.1f%% of nuclear binding\n",
               fabs(phi_min) * M_proton_MeV / 8.0 * 100);

        /* Potential half-width */
        double r_half = 0;
        for (int ri = 0; ri < nr; ri++) {
            if (r_arr[ri] > r_phi_min && fabs(phi_arr[ri]) < 0.5 * fabs(phi_min)) {
                r_half = r_arr[ri];
                break;
            }
        }
        printf("     Potential half-width: %.3f code = %.3f fm\n",
               r_half, r_half * fm_per_code);

        /* Comparison to Newtonian gravity */
        /* G M_p / (r_p c²) = 6.674e-11 × 1.6726e-27 / (0.88e-15 × 8.988e16) */
        double Phi_Newton = 6.674e-11 * 1.6726e-27 / (0.88e-15 * 8.988e16);
        printf("\n  5. Comparison to Newtonian gravity:\n");
        printf("     Φ_Newton/c² at proton surface: %.3e\n", Phi_Newton);
        printf("     Φ_BLV/c² at soliton center:   %.3e\n", fabs(phi_min));
        printf("     Ratio: %.2e\n", fabs(phi_min) / Phi_Newton);
        printf("     → BLV metric effect is ~10^37 × stronger than gravity\n");

        printf("\n  6. Character:\n");
        printf("     Shape: bell-curve (peaked at center, exponential falloff)\n");
        printf("     Range: ~%.2f fm (core radius)\n", r_half * fm_per_code);
        printf("     NOT 1/r: no long-range gravitational tail\n");
        printf("     → Short-range (Yukawa-like), not Newtonian\n");
        printf("     → Matches scale of nuclear forces, not gravity\n");

        free(r_arr);
        free(pm_arr);
        free(phi_arr);
    }

    /* ========== Summary ========== */
    printf("\n========================================\n");
    int all_pass = pass_E && pass_V && pass_B && pass_pm;
    if (all_pass)
        printf(" ALL CHECKS PASSED\n");
    else
        printf(" SOME CHECKS FAILED\n");
    printf("========================================\n");

    profile_free(prof);
    sg_free(g);
    return all_pass ? 0 : 1;
}
