/*
 * torsion.c — Tests for density conservation and field torsion hypotheses
 *
 * Tests:
 *   A. σ-model torsion decay: frame deviation, connection, curvature
 *   B. Finite-λ density conservation: compute deficit ΔQ = ∫(1-ρ)4πr²dr
 *   C. Constrained soliton: add chemical potential μ, enforce ∫ρ = V·ρ₀
 *   D. Comparison with Path 3 (B⁰p → 1/r scalar coupling)
 *
 * Build: make torsion
 * Run:   ../bin/torsion [-lambda L] [-Rmax R]
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define NMAX 100001
#define PI 3.14159265358979323846

/* ========== Global grid ========== */
static int    NG;
static double hG, RmaxG;
static double rG[NMAX], fG[NMAX], fpG[NMAX], rhoG[NMAX];
static double S0G[NMAX];   /* S₀(r) = f'² + 2sin²f/r² */

/* ========== Hedgehog ODE (σ-model: ρ=1, c₄=2, e=1, ρ₀=1) ========== */

static double f_ode_rhs(double r, double fv, double fpv)
{
    if (r < 1e-10) return 0.0;
    double sf = sin(fv), s2f = sin(2.0*fv);
    double c4 = 2.0;
    double P = r*r + 2.0*c4*sf*sf;
    double num = -2.0*r*fpv - c4*fpv*fpv*s2f
                 + s2f*(1.0 + c4*sf*sf/(r*r));
    return num / P;
}

static void rk4_shoot(double a, int N, double h,
                      double *rv, double *fv, double *fpv)
{
    rv[0]  = 0.0;
    fv[0]  = PI;
    fpv[0] = -a;
    for (int i = 0; i < N; i++) {
        double r = rv[i], f = fv[i], fp = fpv[i];
        double k1f  = fp;
        double k1fp = f_ode_rhs(r, f, fp);

        double r2 = r + 0.5*h;
        double k2f  = fp + 0.5*h*k1fp;
        double k2fp = f_ode_rhs(r2, f + 0.5*h*k1f, fp + 0.5*h*k1fp);

        double k3f  = fp + 0.5*h*k2fp;
        double k3fp = f_ode_rhs(r2, f + 0.5*h*k2f, fp + 0.5*h*k2fp);

        double r4 = r + h;
        double k4f  = fp + h*k3fp;
        double k4fp = f_ode_rhs(r4, f + h*k3f, fp + h*k3fp);

        rv[i+1]  = r + h;
        fv[i+1]  = f  + (h/6.0)*(k1f  + 2*k2f  + 2*k3f  + k4f);
        fpv[i+1] = fp + (h/6.0)*(k1fp + 2*k2fp + 2*k3fp + k4fp);
    }
}

static double solve_profile(double Rmax, int N)
{
    double h = Rmax / N;
    double a_lo = 0.5, a_hi = 3.0;

    /* Bisect on slope a */
    double rv[NMAX], fv[NMAX], fpv[NMAX];
    for (int iter = 0; iter < 60; iter++) {
        double a = 0.5*(a_lo + a_hi);
        rk4_shoot(a, N, h, rv, fv, fpv);
        if (fv[N] > 0) a_lo = a; else a_hi = a;
    }
    double a = 0.5*(a_lo + a_hi);
    rk4_shoot(a, N, h, rv, fv, fpv);

    /* Copy to global */
    NG = N; hG = h; RmaxG = Rmax;
    for (int i = 0; i <= N; i++) {
        rG[i]  = rv[i];
        fG[i]  = fv[i];
        fpG[i] = fpv[i];
        rhoG[i] = 1.0;
    }

    /* Compute S₀ */
    for (int i = 0; i <= N; i++) {
        if (rG[i] < 1e-10) {
            S0G[i] = 3.0*a*a;
        } else {
            double sf = sin(fG[i]);
            S0G[i] = fpG[i]*fpG[i] + 2.0*sf*sf/(rG[i]*rG[i]);
        }
    }
    return a;
}

/* ========== ρ BVP solver (Thomas + Newton) ========== */

/* Solve: ρ'' + (2/r)ρ' = ρ·S₀ + λ·ρ·(ρ²-1) - μ
 * BCs:   ρ'(0) = 0,  ρ(R) = rho_inf
 */
static void solve_rho(double lambda, double mu, double rho_inf,
                      int max_newton, double tol)
{
    double al[NMAX], bl[NMAX], cl_arr[NMAX], res[NMAX], dx[NMAX];

    /* Initial guess */
    for (int i = 0; i <= NG; i++) rhoG[i] = 1.0;

    for (int newton = 0; newton < max_newton; newton++) {
        double max_res = 0;

        /* Interior points i = 1..NG-1 */
        for (int i = 1; i < NG; i++) {
            double r  = rG[i];
            double rh = rhoG[i];
            double F  = rh*S0G[i] + lambda*rh*(rh*rh - 1.0) - mu;
            double dF = S0G[i] + lambda*(3.0*rh*rh - 1.0);

            double rpp = (rhoG[i-1] - 2.0*rhoG[i] + rhoG[i+1])/(hG*hG);
            double rp  = (rhoG[i+1] - rhoG[i-1])/(2.0*hG);
            double residual = rpp + 2.0*rp/r - F;

            al[i]     = 1.0/(hG*hG) - 1.0/(hG*r);
            bl[i]     = -2.0/(hG*hG) - dF;
            cl_arr[i] = 1.0/(hG*hG) + 1.0/(hG*r);
            res[i]    = -residual;  /* Newton: RHS = -G(ρ) */

            if (fabs(residual) > max_res) max_res = fabs(residual);
        }

        /* r = 0: L'Hôpital gives 3ρ'' = F(ρ₀) */
        {
            double rh = rhoG[0];
            double F  = rh*S0G[0] + lambda*rh*(rh*rh - 1.0) - mu;
            double dF = S0G[0] + lambda*(3.0*rh*rh - 1.0);
            double rpp3 = 6.0*(rhoG[1] - rhoG[0])/(hG*hG); /* 3·(2(ρ₁-ρ₀)/h²) */
            double residual = rpp3 - F;

            al[0]     = 0.0;
            bl[0]     = -6.0/(hG*hG) - dF;
            cl_arr[0] = 6.0/(hG*hG);
            res[0]    = -residual;  /* Newton: RHS = -G(ρ) */

            if (fabs(residual) > max_res) max_res = fabs(residual);
        }

        /* r = R: ρ(R) = rho_inf */
        al[NG] = 0.0; bl[NG] = 1.0; cl_arr[NG] = 0.0;
        res[NG] = rho_inf - rhoG[NG];

        /* Thomas solve */
        for (int i = 1; i <= NG; i++) {
            double m = al[i] / bl[i-1];
            bl[i]  -= m * cl_arr[i-1];
            res[i] -= m * res[i-1];
        }
        dx[NG] = res[NG] / bl[NG];
        for (int i = NG-1; i >= 0; i--)
            dx[i] = (res[i] - cl_arr[i]*dx[i+1]) / bl[i];

        /* Damped update */
        for (int i = 0; i <= NG; i++) {
            double step = dx[i];
            if (step >  0.05) step =  0.05;
            if (step < -0.05) step = -0.05;
            rhoG[i] += step;
            if (rhoG[i] < 0.3) rhoG[i] = 0.3;
        }

        if (max_res < tol) {
            /* printf("  rho BVP converged in %d Newton iterations (res=%.2e)\n",
                    newton+1, max_res); */
            return;
        }
    }
    printf("  WARNING: rho BVP did not converge (tol=%.1e)\n", tol);
}

/* ========== Integrals ========== */

/* ∫ρ(r)·4πr²dr using Simpson's rule */
static double density_integral(void)
{
    double sum = 0;
    /* Trapezoidal (sufficient for smooth data) */
    for (int i = 0; i < NG; i++) {
        double r0 = rG[i],   rh0 = rhoG[i];
        double r1 = rG[i+1], rh1 = rhoG[i+1];
        sum += 0.5*(rh0*r0*r0 + rh1*r1*r1)*hG;
    }
    return 4.0*PI*sum;
}

/* Vacuum integral = (4/3)πR³ */
static double vacuum_integral(void)
{
    return (4.0/3.0)*PI*RmaxG*RmaxG*RmaxG;
}

/* ∫(1-ρ)·4πr²dr  (density deficit) */
static double deficit_integral(void)
{
    double sum = 0;
    for (int i = 0; i < NG; i++) {
        double r0 = rG[i],   d0 = 1.0 - rhoG[i];
        double r1 = rG[i+1], d1 = 1.0 - rhoG[i+1];
        sum += 0.5*(d0*r0*r0 + d1*r1*r1)*hG;
    }
    return 4.0*PI*sum;
}

/* ∫B⁰·4πr²dr = topological charge (should be 1.0) */
static double charge_integral(void)
{
    double sum = 0;
    for (int i = 0; i < NG; i++) {
        double r0 = rG[i], r1 = rG[i+1];
        double B0, B1;
        double a = fabs(fpG[0]);
        if (r0 < 1e-10) B0 = a*a*a/(2.0*PI*PI);
        else B0 = fabs(sin(fG[i])*sin(fG[i])*fpG[i])/(2.0*PI*PI*r0*r0);
        if (r1 < 1e-10) B1 = a*a*a/(2.0*PI*PI);
        else B1 = fabs(sin(fG[i+1])*sin(fG[i+1])*fpG[i+1])/(2.0*PI*PI*r1*r1);
        sum += 0.5*(B0*r0*r0 + B1*r1*r1)*hG;
    }
    return 4.0*PI*sum;
}

/* ========== Power-law fit: log(y) = A - n·log(r) for r in [r1,r2] ========== */

static void fit_power_law(int i_start, int i_end,
                          double *rv, double *yv,
                          double *C_out, double *n_out)
{
    /* Linear regression: ln(y) = ln(C) - n·ln(r) */
    double sx = 0, sy = 0, sxx = 0, sxy = 0;
    int count = 0;
    for (int i = i_start; i <= i_end; i++) {
        if (yv[i] <= 0 || rv[i] <= 0) continue;
        double x = log(rv[i]);
        double y = log(yv[i]);
        sx += x; sy += y; sxx += x*x; sxy += x*y;
        count++;
    }
    if (count < 2) { *C_out = 0; *n_out = 0; return; }
    double n = (count*sxy - sx*sy) / (count*sxx - sx*sx);
    double lnC = (sy - n*sx) / count;
    *n_out = -n;  /* y ~ C·r^{-n} means ln(y) = ln(C) - n·ln(r) */
    *C_out = exp(lnC);
}

/* ========== TEST A: σ-model torsion decay ========== */

static void test_torsion_decay(void)
{
    printf("\n");
    printf("================================================================\n");
    printf("  TEST A: σ-MODEL TORSION DECAY RATES\n");
    printf("================================================================\n\n");

    double a = fabs(fpG[0]);
    printf("  Profile: a = %.4f,  f(R=%.0f) = %.2e\n\n", a, RmaxG, fG[NG]);

    /* Compute quantities at each grid point */
    double frame_dev[NMAX];   /* D(r) = 2|sin(f/4)| */
    double connection[NMAX];  /* |ω(r)| = √S₀ */
    double baryon[NMAX];      /* B⁰(r) */
    double path3_pot[NMAX];   /* p(r) = g·Q_enc/(4πr) — Path 3 potential */

    double Q_enc = 0; /* enclosed topological charge */

    for (int i = 0; i <= NG; i++) {
        double r = rG[i], f = fG[i], fp = fpG[i];

        frame_dev[i]  = 2.0*fabs(sin(f/4.0));
        connection[i] = sqrt(S0G[i]);

        if (r < 1e-10) {
            baryon[i] = a*a*a/(2.0*PI*PI);
        } else {
            baryon[i] = fabs(sin(f)*sin(f)*fp)/(2.0*PI*PI*r*r);
        }

        /* Accumulate enclosed charge for Path 3 */
        if (i > 0) {
            double B0 = baryon[i-1], B1 = baryon[i];
            double r0 = rG[i-1], r1 = rG[i];
            Q_enc += 0.5*4.0*PI*(B0*r0*r0 + B1*r1*r1)*hG;
        }
        if (r < 1e-10) path3_pot[i] = 0;
        else path3_pot[i] = Q_enc/(4.0*PI*r);  /* g_top=1, κ²=1 */
    }

    /* Print sample points */
    printf("  %-8s %-12s %-12s %-12s %-12s %-12s\n",
           "r", "f(r)", "D(r)", "|ω(r)|", "B⁰(r)", "p_3(r)");
    printf("  %-8s %-12s %-12s %-12s %-12s %-12s\n",
           "", "(profile)", "(frame dev)", "(connection)", "(baryon)", "(Path 3)");
    printf("  -------------------------------------------------------------------\n");

    int sample_r[] = {1, 2, 3, 4, 5, 7, 10, 15, 20};
    int n_samples = 9;
    for (int s = 0; s < n_samples; s++) {
        int target_r = sample_r[s];
        int idx = (int)(target_r / hG + 0.5);
        if (idx > NG) break;
        printf("  %-8.1f %-12.4e %-12.4e %-12.4e %-12.4e %-12.4e\n",
               rG[idx], fG[idx], frame_dev[idx], connection[idx],
               baryon[idx], path3_pot[idx]);
    }

    /* Fit decay exponents at large r */
    int i_fit_start = (int)(5.0/hG + 0.5);
    int i_fit_end   = (int)(15.0/hG + 0.5);
    if (i_fit_end > NG) i_fit_end = NG - 1;

    double C, n;
    printf("\n  Decay fits (r = 5 to 15):\n");

    fit_power_law(i_fit_start, i_fit_end, rG, frame_dev, &C, &n);
    printf("  Frame deviation D(r):  D ~ %.2f / r^%.2f", C, n);
    printf("   (expected: ~2/r² from f~2/r²)\n");

    fit_power_law(i_fit_start, i_fit_end, rG, connection, &C, &n);
    printf("  Connection |ω(r)|:     |ω| ~ %.2f / r^%.2f", C, n);
    printf("   (expected: ~1/r³)\n");

    fit_power_law(i_fit_start, i_fit_end, rG, baryon, &C, &n);
    printf("  Baryon density B⁰(r):  B⁰ ~ %.3f / r^%.2f", C, n);
    printf("   (expected: ~1/r⁹)\n");

    fit_power_law(i_fit_start, i_fit_end, rG, path3_pot, &C, &n);
    printf("  Path 3 potential p(r): p ~ %.4f / r^%.2f", C, n);
    printf("   (expected: Q/(4πr) = 1/r)\n");

    /* Key comparison */
    printf("\n  CONCLUSION:\n");
    printf("  The frame torsion (connection) decays as 1/r³ — short range.\n");
    printf("  Only the Path 3 constraint potential gives 1/r (long range).\n");
    printf("  Frame deviation 1/r² is intermediate but NOT a force source.\n");
}

/* ========== TEST B: Finite-λ density conservation ========== */

static void test_density_conservation(double lambda)
{
    printf("\n");
    printf("================================================================\n");
    printf("  TEST B: DENSITY CONSERVATION AT FINITE λ\n");
    printf("================================================================\n\n");

    printf("  λ = %.0f,  R_max = %.1f,  N = %d,  h = %.4f\n\n", lambda, RmaxG, NG, hG);

    /* Save sigma-model ρ for later */
    double rho_sigma[NMAX];
    for (int i = 0; i <= NG; i++) rho_sigma[i] = 1.0;

    /* Solve ρ BVP (unconstrained: μ=0, rho_inf=1) */
    solve_rho(lambda, 0.0, 1.0, 200, 1e-10);

    printf("  ρ profile (unconstrained, μ = 0):\n");
    printf("  %-8s %-12s %-12s\n", "r", "ρ(r)", "1-ρ(r)");
    printf("  ----------------------------------\n");
    printf("  %-8.3f %-12.6f %-12.6f\n", rG[0], rhoG[0], 1.0-rhoG[0]);
    int r_vals[] = {1, 2, 3, 5, 10, 20};
    for (int s = 0; s < 6; s++) {
        int idx = (int)(r_vals[s]/hG + 0.5);
        if (idx > NG) break;
        printf("  %-8.3f %-12.6f %-12.2e\n", rG[idx], rhoG[idx], 1.0-rhoG[idx]);
    }

    double Q_dens = density_integral();
    double Q_vac  = vacuum_integral();
    double deficit = deficit_integral();

    printf("\n  Global integrals:\n");
    printf("  ∫ρ·4πr²dr  = %.6f  (soliton)\n", Q_dens);
    printf("  ∫1·4πr²dr  = %.6f  (vacuum)\n", Q_vac);
    printf("  Deficit ΔQ  = %.6f\n", deficit);
    printf("  ΔQ / Q_vac  = %.2e  (fractional)\n", deficit/Q_vac);

    /* Yukawa decay of (1-ρ) */
    double mass_rho = sqrt(2.0*lambda);
    printf("\n  ρ mode mass: m_ρ = √(2λ) = %.2f\n", mass_rho);
    printf("  Yukawa range: 1/m_ρ = %.4f code = %.4f fm\n",
           1.0/mass_rho, 0.5624/mass_rho);

    /* Fit exponential decay of (1-ρ) at intermediate r */
    printf("\n  Exponential decay fit (ln(1-ρ) vs r):\n");
    int i_start = (int)(0.5/hG), i_end = (int)(2.0/hG);
    if (lambda > 5000) { i_start = (int)(0.2/hG); i_end = (int)(0.8/hG); }
    if (i_end > NG - 1) i_end = NG - 1;
    double sx = 0, sy = 0, sxx = 0, sxy = 0; int count = 0;
    for (int i = i_start; i <= i_end; i++) {
        double d = 1.0 - rhoG[i];
        if (d < 1e-15) continue;
        double x = rG[i], y = log(d);
        sx += x; sy += y; sxx += x*x; sxy += x*y;
        count++;
    }
    if (count > 2) {
        double slope = (count*sxy - sx*sy)/(count*sxx - sx*sx);
        printf("  Fitted decay rate: κ = %.2f  (expected: √(2λ) = %.2f)\n",
               -slope, mass_rho);
    }

    printf("\n  CONCLUSION:\n");
    printf("  The soliton DEPLETES density at its core (ρ₀ = %.4f < 1).\n", rhoG[0]);
    printf("  Total density deficit: ΔQ = %.4f (finite, localized).\n", deficit);
    printf("  The deficit is exponentially localized (Yukawa, range %.4f).\n",
           1.0/mass_rho);
    printf("  NO long-range 1/r² depletion zone exists.\n");
    printf("  Narrative spec's \"density concentration\" is REVERSED: core has LESS density.\n");

    /* Restore ρ=1 for subsequent tests */
    for (int i = 0; i <= NG; i++) rhoG[i] = rho_sigma[i];
}

/* ========== TEST C: Constrained soliton ========== */

static void test_constrained(double lambda)
{
    printf("\n");
    printf("================================================================\n");
    printf("  TEST C: DENSITY-CONSTRAINED SOLITON (λ = %.0f)\n", lambda);
    printf("================================================================\n\n");

    double Q_vac = vacuum_integral();
    printf("  Vacuum density integral: Q_vac = ∫1·4πr²dr = %.4f\n", Q_vac);

    /* Step 1: unconstrained (μ=0) to get baseline deficit */
    solve_rho(lambda, 0.0, 1.0, 200, 1e-10);
    double Q_uncon = density_integral();
    double deficit = Q_vac - Q_uncon;
    printf("  Unconstrained: Q_dens = %.4f,  deficit = %.4f\n", Q_uncon, deficit);

    /* Step 2: estimate μ needed */
    /* At far field: ρ_∞ = 1 + μ/(2λ). Total integral shift ≈ (4/3)πR³·μ/(2λ) */
    double mu_est = deficit * 2.0 * lambda / ((4.0/3.0)*PI*RmaxG*RmaxG*RmaxG);
    printf("  Estimated μ for conservation: %.6e\n", mu_est);
    printf("  Expected far-field shift: ρ_∞ - 1 = μ/(2λ) = %.6e\n\n",
           mu_est/(2.0*lambda));

    /* Step 3: bisect on μ to enforce ∫ρ = Q_vac */
    double mu_lo = 0.0, mu_hi = mu_est * 5.0;

    for (int iter = 0; iter < 40; iter++) {
        double mu = 0.5*(mu_lo + mu_hi);
        double rho_inf = 1.0 + mu/(2.0*lambda);
        /* Re-initialize */
        for (int i = 0; i <= NG; i++) rhoG[i] = 1.0;
        solve_rho(lambda, mu, rho_inf, 200, 1e-10);
        double Q = density_integral();
        if (Q < Q_vac) mu_lo = mu; else mu_hi = mu;
    }

    double mu_star = 0.5*(mu_lo + mu_hi);
    double rho_inf = 1.0 + mu_star/(2.0*lambda);

    /* Re-solve at converged μ */
    for (int i = 0; i <= NG; i++) rhoG[i] = 1.0;
    solve_rho(lambda, mu_star, rho_inf, 200, 1e-10);
    double Q_con = density_integral();

    printf("  Converged: μ* = %.6e\n", mu_star);
    printf("  ρ_∞ = 1 + %.6e  (far-field shift)\n", rho_inf - 1.0);
    printf("  ∫ρ = %.6f  vs  Q_vac = %.6f  (error: %.2e)\n\n",
           Q_con, Q_vac, fabs(Q_con - Q_vac)/Q_vac);

    /* Step 4: print constrained ρ profile */
    printf("  Constrained ρ profile:\n");
    printf("  %-8s %-14s %-14s %-14s\n", "r", "ρ(r)", "ρ-ρ_∞", "ρ-1");
    printf("  -------------------------------------------------------\n");
    printf("  %-8.3f %-14.8f %-14.6e %-14.6e\n",
           rG[0], rhoG[0], rhoG[0]-rho_inf, rhoG[0]-1.0);
    double check_r[] = {0.5, 1.0, 2.0, 3.0, 5.0, 10.0, 15.0, 20.0, 25.0};
    for (int s = 0; s < 9; s++) {
        int idx = (int)(check_r[s]/hG + 0.5);
        if (idx > NG) break;
        printf("  %-8.3f %-14.8f %-14.6e %-14.6e\n",
               rG[idx], rhoG[idx], rhoG[idx]-rho_inf, rhoG[idx]-1.0);
    }

    /* Step 5: check far-field approach to ρ_∞ */
    printf("\n  Far-field analysis:\n");
    double mass_rho = sqrt(2.0*lambda*(3.0*rho_inf*rho_inf - 1.0));
    printf("  ρ mode mass around ρ_∞: m = %.2f\n", mass_rho);
    printf("  Yukawa range: 1/m = %.5f code = %.5f fm\n",
           1.0/mass_rho, 0.5624/mass_rho);

    /* Check: is far-field (ρ - ρ_∞) decaying exponentially or as 1/r²? */
    printf("\n  Approach to asymptote (ρ - ρ_∞) vs r:\n");
    printf("  %-8s %-14s %-14s\n", "r", "ρ-ρ_∞", "expected Yukawa");
    printf("  -----------------------------------------------\n");
    /* Find amplitude of Yukawa from intermediate r */
    int i_mid = (int)(3.0/hG + 0.5);
    if (lambda > 5000) i_mid = (int)(1.0/hG + 0.5);
    double A_yuk = 0;
    if (i_mid <= NG && rG[i_mid] > 0) {
        double dev = rhoG[i_mid] - rho_inf;
        if (fabs(dev) > 1e-15)
            A_yuk = dev * rG[i_mid] * exp(mass_rho * rG[i_mid]);
    }

    double check_r2[] = {1.0, 2.0, 3.0, 5.0, 8.0, 10.0, 15.0, 20.0};
    int nc2 = 8;
    if (lambda > 5000) {
        double temp[] = {0.2, 0.4, 0.6, 0.8, 1.0, 1.5, 2.0, 3.0};
        for (int s = 0; s < 8; s++) check_r2[s] = temp[s];
    }
    for (int s = 0; s < nc2; s++) {
        int idx = (int)(check_r2[s]/hG + 0.5);
        if (idx > NG) break;
        double dev = rhoG[idx] - rho_inf;
        double yuk = (rG[idx] > 0) ? A_yuk * exp(-mass_rho*rG[idx])/rG[idx] : 0;
        printf("  %-8.2f %-14.6e %-14.6e\n", rG[idx], dev, yuk);
    }

    /* Gravitational force comparison */
    double rho_shift = rho_inf - 1.0;
    printf("\n  CONCLUSION:\n");
    printf("  The constraint shifts the far-field density by Δρ = %.2e\n", rho_shift);
    printf("  (= %.2e of ρ₀ — essentially unmeasurable).\n", rho_shift);
    printf("  The approach to ρ_∞ is EXPONENTIAL (Yukawa, range %.4f).\n",
           1.0/mass_rho);
    printf("  There is NO 1/r² depletion zone.\n");
    printf("  The constraint redistributes ~%.2e density units to infinity.\n", deficit);
    printf("  Over volume (4/3)πR³ = %.0f, this is Δρ/ρ₀ = %.2e.\n",
           Q_vac, rho_shift);
    printf("  The density mode is MASSIVE (m=%.1f), so perturbations are Yukawa.\n",
           mass_rho);
    printf("  For 1/r depletion, would need MASSLESS density mode (flat potential).\n");

    /* Restore ρ=1 */
    for (int i = 0; i <= NG; i++) rhoG[i] = 1.0;
}

/* ========== TEST D: Comparison table ========== */

static void test_comparison(void)
{
    printf("\n");
    printf("================================================================\n");
    printf("  TEST D: MECHANISM COMPARISON\n");
    printf("================================================================\n\n");

    double a = fabs(fpG[0]);

    /* At representative distances */
    printf("  Quantity decay at various distances (σ-model):\n\n");
    printf("  %-6s %-12s %-12s %-12s %-12s %-12s\n",
           "r", "f(r)", "D=2sin(f/4)", "|ω|=√S₀", "B⁰", "p₃=Q/(4πr)");
    printf("  %-6s %-12s %-12s %-12s %-12s %-12s\n",
           "", "(profile)", "(frame dev)", "(torsion)", "(source)", "(potential)");
    printf("  -----------------------------------------------------------------------\n");

    double Q_total = charge_integral();

    for (double r_target = 1.0; r_target <= 20.5; r_target += 1.0) {
        int idx = (int)(r_target/hG + 0.5);
        if (idx > NG) break;
        double r = rG[idx], f = fG[idx], fp = fpG[idx];
        double D = 2.0*fabs(sin(f/4.0));
        double omega = sqrt(S0G[idx]);
        double B0;
        if (r < 1e-10) B0 = a*a*a/(2.0*PI*PI);
        else B0 = fabs(sin(f)*sin(f)*fp)/(2.0*PI*PI*r*r);
        double p3 = Q_total/(4.0*PI*r);

        printf("  %-6.0f %-12.4e %-12.4e %-12.4e %-12.4e %-12.4e\n",
               r, f, D, omega, B0, p3);
    }

    /* Summary */
    printf("\n  At r = 10 (≈ 5.6 fm from soliton center):\n");
    int i10 = (int)(10.0/hG + 0.5);
    if (i10 <= NG) {
        double r = rG[i10];
        double D = 2.0*fabs(sin(fG[i10]/4.0));
        double omega = sqrt(S0G[i10]);
        double p3 = Q_total/(4.0*PI*r);

        printf("  Frame deviation:  D = %.4e  (1/r² decay)\n", D);
        printf("  Frame torsion:    |ω| = %.4e  (1/r³ decay)\n", omega);
        printf("  Path 3 potential: p = %.4e  (1/r decay)\n", p3);
        printf("  Ratio p/D = %.1f  — Path 3 is %.0f× stronger at this distance\n",
               p3/D, p3/D);
        printf("  Ratio p/|ω| = %.1f  — Path 3 is %.0f× stronger\n",
               p3/omega, p3/omega);
    }

    /* Boosted soliton analysis */
    printf("\n  BOOSTED SOLITON (velocity v along ẑ):\n\n");
    printf("  At transverse distance d, soliton at origin:\n");
    printf("  Static frame deviation:  D(d) ~ 2/d²\n");
    printf("  Boost velocity field:    |v_boost| = v·|∂q/∂z| ~ v·sinf/(2d)\n");
    printf("  For large d: |v_boost| ~ v/d³  (from f~2/d², ∂r̂/∂z~1/d)\n\n");

    printf("  %-6s %-14s %-14s %-14s\n",
           "d", "D(d) static", "|v_boost|/v", "p₃(d) Path 3");
    printf("  -----------------------------------------------------------\n");
    for (double d = 2.0; d <= 20.5; d += 2.0) {
        int idx = (int)(d/hG + 0.5);
        if (idx > NG) break;
        double f = fG[idx];
        double D = 2.0*fabs(sin(f/4.0));
        /* v_boost at (d,0,0): ∂q̂/∂z = sin(f/2)·(ẑ/d)·σ component */
        double v_boost = fabs(sin(f/2.0)) / d;  /* relative to v */
        double p3 = Q_total/(4.0*PI*d);

        printf("  %-6.0f %-14.4e %-14.4e %-14.4e\n", d, D, v_boost, p3);
    }
    printf("\n  The boost creates a transverse velocity field that decays as ~1/d³.\n");
    printf("  This is NOT 1/d (unlike EM magnetic field from wire current).\n");
    printf("  The 1/d magnetic-like effect only arises through the Path 3\n");
    printf("  constraint sector: □p = g_top·B⁰ gives gravitomagnetic 1/r.\n");

    printf("\n  ================================================================\n");
    printf("  OVERALL SUMMARY\n");
    printf("  ================================================================\n\n");
    printf("  Hypothesis 1 (density conservation → 1/r² gravity):\n");
    printf("    NEGATIVE. Density is NOT conserved. Soliton has LESS density\n");
    printf("    at core (not more). Even with forced conservation, the density\n");
    printf("    mode is massive → Yukawa depletion, not 1/r².\n\n");

    printf("  Hypothesis 2 (frame torsion → non-local 1/r effect):\n");
    printf("    PARTIALLY CONFIRMED. The soliton IS a frame twist (not density\n");
    printf("    concentration). The twist extends non-locally. But the bulk\n");
    printf("    frame torsion decays as 1/r³ (too fast for gravity).\n\n");

    printf("  Hypothesis 3 (wire analogy — moving soliton → 1/r torsion vector):\n");
    printf("    NEGATIVE for bulk field (1/d³, not 1/d).\n");
    printf("    POSITIVE for constraint sector (Path 3): the degenerate field p\n");
    printf("    sourced by B⁰ gives BOTH 1/r scalar potential (static) AND\n");
    printf("    1/r vector potential (moving soliton). This IS the wire analogy.\n\n");

    printf("  KEY FINDING:\n");
    printf("    The Path 3 constraint channel (p field sourced by B⁰) IS the\n");
    printf("    \"non-local torsion\" the theory needs. It gives 1/r (static)\n");
    printf("    and gravitomagnetic effects (moving). The frame torsion\n");
    printf("    interpretation: p = scalar component of Cl⁺(3,0,1) frame torsion.\n");
    printf("    But the coupling g_top remains a free parameter.\n\n");

    printf("  WHAT WOULD FIX IT:\n");
    printf("    For density-based gravity: need massless density mode.\n");
    printf("    This requires a FLAT-BOTTOM potential V(ρ) with degenerate\n");
    printf("    vacuum (V''(ρ₀) = 0), instead of the standard Mexican hat\n");
    printf("    V = (λ/4)(ρ²-ρ₀²)² which gives massive ρ mode (m = √(2λ)ρ₀).\n");
    printf("    A BPS potential V = μ² (1-cosf) has this property for the\n");
    printf("    angular mode — but NOT for the radial (density) mode.\n");
}

/* ========== MAIN ========== */

int main(int argc, char **argv)
{
    double lambda = 100.0;
    double Rmax   = 30.0;
    int    N      = 30000;

    /* Parse args */
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-lambda") == 0 && i+1 < argc)
            lambda = atof(argv[++i]);
        else if (strcmp(argv[i], "-Rmax") == 0 && i+1 < argc)
            Rmax = atof(argv[++i]);
        else if (strcmp(argv[i], "-N") == 0 && i+1 < argc)
            N = atoi(argv[++i]);
    }

    if (N > NMAX-1) { fprintf(stderr, "N too large (max %d)\n", NMAX-1); return 1; }

    printf("Torsion Analysis: Testing Density Conservation & Field Torsion\n");
    printf("Parameters: λ = %.0f,  R_max = %.1f,  N = %d\n", lambda, Rmax, N);

    /* Solve σ-model profile */
    printf("\nSolving σ-model hedgehog...\n");
    double a = solve_profile(Rmax, N);
    printf("  Converged: a = %.6f,  Q = %.6f\n", a, charge_integral());

    /* Run all tests */
    test_torsion_decay();
    test_density_conservation(lambda);
    test_constrained(lambda);
    test_comparison();

    return 0;
}
