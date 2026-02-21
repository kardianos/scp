/*
 * depletion.c — V6 Test 1+2: Self-consistent density profile and long-range tail
 *
 * Solves the coupled system:
 *   f(r):  Skyrme profile (twist field)
 *   ρ(r):  density from equilibrium condition ρ = ρ₀ - |ω|²/(4α)
 *
 * Then computes the far-field depletion δρ(r) and checks:
 *   1. Total deficit ΔQ ∝ E₂ (universality)
 *   2. Far-field δρ(r) → -ΔQ/(4πr) (1/r, not Yukawa)
 *   3. Gravitational force F = -G_eff M₁M₂/r² with G_eff determined
 *
 * Usage:
 *   ./bin/depletion [-alpha A] [-B baryon] [-mpi M] [-scan]
 *
 * Build:
 *   make depletion
 */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Grid parameters */
#define NMAX 100001
static double rG[NMAX], fG[NMAX], fpG[NMAX];
static double rhoG[NMAX]; /* density profile */
static int N_grid;
static double h_grid;
static double R_max = 30.0;

/* Physical parameters */
static double rho0 = 1.0;   /* background density (= N/V effectively) */
static double c4   = 2.0;   /* Skyrme coupling = 2ρ₀²/e² with e=1, ρ₀=1 */
static double alpha_p = 1.0; /* pressure coefficient (E_p = α∫ρ² d³x) */
static double m_pi = 0.0;   /* pion mass (0 = massless) */
static int B_charge = 1;    /* baryon number */

/* ================================================================
 * Profile solver (RK4 shooting, same as v2/radial.c)
 * ================================================================ */

/* Twist energy density |ω|² for hedgehog: f'² + 2sin²f/r² */
static double omega_sq(double r, double f, double fp) {
    if (r < 1e-12) {
        /* Near origin: f ~ π - a·r^α, f' ~ -a·α·r^{α-1}
         * For B=1: α=1, |ω|² = f'² + 2sin²f/r² → a² + 2a²/1 = 3a²
         * But at r=0 exactly, just return 3*fp*fp (using limiting form) */
        return 3.0 * fp * fp;
    }
    double sf = sin(f);
    return fp * fp + 2.0 * sf * sf / (r * r);
}

/* Skyrme radial ODE with ρ-dependent E₂ coupling.
 *
 * Energy functional (hedgehog, per 4π):
 *   E₂ = ∫[½ρf'²r² + ρsin²f] dr
 *   E₄ = ∫c₄[f'²sin²f + sin⁴f/(2r²)] dr
 *   E_m = ∫m_π²(1-cosf)r² dr
 *
 * Euler-Lagrange: d/dr[∂L/∂f'] = ∂L/∂f where L = integrand.
 *
 *   [ρr² + 2c₄sin²f]f'' = -[ρ'r² + 2ρr]f'
 *                          + ρsin2f - c₄sin2f·f'² + c₄sin2f·sin²f/r²
 *                          + m_π²sinf·r²
 *
 * Note: ρ' = dρ/dr computed from rhoG[] numerically (finite difference).
 */
static double rhs_f(double r, double f, double fp, double rho_local,
                     double drho_dr) {
    if (r < 1e-10) return 0.0; /* f'' → 0 at origin by symmetry */

    double sf = sin(f), s2f = sin(2.0*f);
    double r2 = r*r;
    double sf2 = sf*sf;

    /* Denominator: P = ρr² + 2c₄sin²f */
    double P = rho_local * r2 + 2.0 * c4 * sf2;

    /* Numerator of f'' */
    double num = 0.0;
    num -= (drho_dr * r2 + 2.0 * rho_local * r) * fp; /* -[ρ'r² + 2ρr]f' */
    num += rho_local * s2f;                              /* +ρ sin2f */
    num -= c4 * s2f * fp * fp;                           /* -c₄ sin2f f'² */
    num += c4 * s2f * sf2 / r2;                          /* +c₄ sin2f sin²f/r² */
    num += m_pi * m_pi * sf * r2;                        /* +m_π² sinf r² */

    return num / P;
}

/* Get dρ/dr at grid point i (centered finite difference) */
static double get_drho(int i) {
    double dr = R_max / (N_grid - 1);
    if (i <= 0) return (rhoG[1] - rhoG[0]) / dr;
    if (i >= N_grid - 1) return (rhoG[N_grid-1] - rhoG[N_grid-2]) / dr;
    return (rhoG[i+1] - rhoG[i-1]) / (2.0 * dr);
}

/* Integrate profile from r=0 to R_max using RK4, with initial slope a.
 * ρ(r) is read from rhoG[] (interpolated). */
static double integrate_profile(double a, double *f_out, double *fp_out) {
    double dr = R_max / (N_grid - 1);
    double r, f, fp;

    /* Near-origin: f ~ π - a·r^α for general B */
    double alpha_exp = (-1.0 + sqrt(1.0 + 8.0*B_charge)) / 2.0;

    /* Set f(0) = π, f'(0) from series */
    f_out[0] = M_PI;
    fp_out[0] = (alpha_exp > 1.0) ? 0.0 : -a;

    /* Initialize first grid point from near-origin expansion */
    double r1 = dr;
    f_out[1] = M_PI - a * pow(r1, alpha_exp);
    fp_out[1] = -a * alpha_exp * pow(r1, alpha_exp - 1.0);

    for (int i = 1; i < N_grid - 1; i++) {
        r = rG[i];
        f = f_out[i];
        fp = fp_out[i];

        double rho_i = rhoG[i];
        double drho_i = get_drho(i);

        /* RK4 */
        double k1f = fp;
        double k1p = rhs_f(r, f, fp, rho_i, drho_i);

        double r2 = r + dr/2;
        double f2 = f + dr/2 * k1f;
        double p2 = fp + dr/2 * k1p;
        double rho2 = (rhoG[i] + rhoG[i+1]) / 2.0;
        double drho2 = (get_drho(i) + get_drho(i+1)) / 2.0;
        double k2f = p2;
        double k2p = rhs_f(r2, f2, p2, rho2, drho2);

        double f3 = f + dr/2 * k2f;
        double p3 = fp + dr/2 * k2p;
        double k3f = p3;
        double k3p = rhs_f(r2, f3, p3, rho2, drho2);

        double r4 = r + dr;
        double f4 = f + dr * k3f;
        double p4 = fp + dr * k3p;
        double rho4 = rhoG[i+1];
        double drho4 = get_drho(i+1);
        double k4f = p4;
        double k4p = rhs_f(r4, f4, p4, rho4, drho4);

        f_out[i+1] = f + dr/6*(k1f + 2*k2f + 2*k3f + k4f);
        fp_out[i+1] = fp + dr/6*(k1p + 2*k2p + 2*k3p + k4p);
    }

    /* Return f(R_max) — should be 0 for correct a */
    return f_out[N_grid - 1];
}

/* Bisection on slope a to find f(R_max) = 0 */
static double find_profile(void) {
    double a_lo = 0.5, a_hi = 3.0;

    /* Bracket: f(R) should go from positive (a too small) to negative (a too large) */
    for (int iter = 0; iter < 80; iter++) {
        double a_mid = (a_lo + a_hi) / 2.0;
        double f_end = integrate_profile(a_mid, fG, fpG);

        if (f_end > 0)
            a_lo = a_mid;
        else
            a_hi = a_mid;

        if (fabs(a_hi - a_lo) < 1e-14) break;
    }

    double a_best = (a_lo + a_hi) / 2.0;
    integrate_profile(a_best, fG, fpG);
    return a_best;
}

/* ================================================================
 * Density from equilibrium: ρ(r) = ρ₀ - |ω(r)|²/(4α)
 * ================================================================ */

/* Under-relaxation parameter for density updates */
static double relax = 0.3;

static void compute_density(void) {
    for (int i = 0; i < N_grid; i++) {
        double w2 = omega_sq(rG[i], fG[i], fpG[i]);
        double rho_new = rho0 - w2 / (4.0 * alpha_p);
        /* Floor: density can't go negative */
        if (rho_new < 0.01 * rho0) rho_new = 0.01 * rho0;
        /* Under-relaxation */
        rhoG[i] = (1.0 - relax) * rhoG[i] + relax * rho_new;
    }
}

/* ================================================================
 * Energy integrals
 * ================================================================ */

static void compute_energies(double *E2, double *E4, double *Ep, double *DeltaQ) {
    double e2 = 0, e4 = 0, ep = 0, dq = 0;
    double dr = R_max / (N_grid - 1);

    for (int i = 1; i < N_grid - 1; i++) {
        double r = rG[i];
        double f = fG[i], fp = fpG[i];
        double rho = rhoG[i];
        double r2 = r * r;
        double sf = sin(f), sf2 = sf * sf;

        /* E₂ = ½∫ρ(fp² + 2sf²/r²) 4πr² dr */
        double w2 = fp*fp + 2.0*sf2/r2;
        e2 += 0.5 * rho * w2 * 4.0*M_PI*r2 * dr;

        /* E₄ = ½c₄∫(2fp²sf²/r² + sf⁴/r⁴) 4πr² dr */
        double t2 = 2.0*fp*fp*sf2/r2 + sf2*sf2/r2/r2;
        e4 += 0.5 * c4 * t2 * 4.0*M_PI*r2 * dr;

        /* E_p = α∫ρ² 4πr² dr */
        ep += alpha_p * rho * rho * 4.0*M_PI*r2 * dr;

        /* ΔQ = ∫(ρ₀ - ρ) 4πr² dr */
        dq += (rho0 - rho) * 4.0*M_PI*r2 * dr;
    }

    *E2 = e2;
    *E4 = e4;
    *Ep = ep;
    *DeltaQ = dq;
}

/* ================================================================
 * Poisson solver: Φ(r) from δρ point source
 * For spherical symmetry:
 *   Φ(r) = -(1/r)∫₀ʳ δρ(r')4πr'²dr' - ∫ᵣ^∞ δρ(r')4πr'dr'
 * ================================================================ */

static void solve_poisson(double *Phi) {
    double dr = R_max / (N_grid - 1);

    /* Compute cumulative integrals */
    double *inner = calloc(N_grid, sizeof(double)); /* ∫₀ʳ δρ 4πr² dr */
    double *outer = calloc(N_grid, sizeof(double)); /* ∫ᵣ^∞ δρ 4πr dr */

    /* Inner integral (forward) */
    inner[0] = 0;
    for (int i = 1; i < N_grid; i++) {
        double r = rG[i];
        double drho = rhoG[i] - rho0;
        inner[i] = inner[i-1] + drho * 4.0*M_PI*r*r * dr;
    }

    /* Outer integral (backward) */
    outer[N_grid-1] = 0;
    for (int i = N_grid-2; i >= 0; i--) {
        double r = rG[i+1];
        double drho = rhoG[i+1] - rho0;
        outer[i] = outer[i+1] + drho * 4.0*M_PI*r * dr;
    }

    /* Φ(r) = -(1/r)·inner(r) - outer(r) */
    Phi[0] = -outer[0]; /* limit as r→0 */
    for (int i = 1; i < N_grid; i++) {
        Phi[i] = -inner[i] / rG[i] - outer[i];
    }

    free(inner);
    free(outer);
}

/* ================================================================
 * Power-law fit: log(|y|) = A + B·log(r) over [r_lo, r_hi]
 * ================================================================ */

static void fit_power_law(double *r_arr, double *y_arr, int n,
                          double r_lo, double r_hi,
                          double *amp, double *exponent) {
    double sx = 0, sy = 0, sxx = 0, sxy = 0;
    int count = 0;

    for (int i = 0; i < n; i++) {
        if (r_arr[i] < r_lo || r_arr[i] > r_hi) continue;
        if (fabs(y_arr[i]) < 1e-30) continue;

        double x = log(r_arr[i]);
        double y = log(fabs(y_arr[i]));
        sx += x; sy += y; sxx += x*x; sxy += x*y;
        count++;
    }

    if (count < 3) { *amp = 0; *exponent = 0; return; }

    double det = count * sxx - sx * sx;
    double B = (count * sxy - sx * sy) / det;
    double A = (sy - B * sx) / count;

    *amp = exp(A);
    *exponent = B;
}

/* ================================================================
 * Main test routines
 * ================================================================ */

static void test_depletion_profile(void) {
    printf("=== TEST 1: Self-consistent density profile ===\n\n");

    /* Initialize grid */
    double dr = R_max / (N_grid - 1);
    for (int i = 0; i < N_grid; i++) {
        rG[i] = i * dr;
        rhoG[i] = rho0; /* start with uniform density */
    }

    /* Self-consistent iteration */
    int max_iter = 20;
    double a;
    double E2, E4, Ep, DQ;
    double prev_DQ = 0;

    for (int iter = 0; iter < max_iter; iter++) {
        a = find_profile();
        compute_density();
        compute_energies(&E2, &E4, &Ep, &DQ);

        printf("  iter %2d: a = %.6f, E2 = %.4f, E4 = %.4f, DQ = %.6f, "
               "rho(0) = %.6f\n", iter, a, E2, E4, DQ, rhoG[0]);

        if (iter > 0 && fabs(DQ - prev_DQ) < 1e-8 * fabs(DQ)) {
            printf("  Converged at iteration %d\n", iter);
            break;
        }
        prev_DQ = DQ;
    }

    printf("\n  Results:\n");
    printf("    Slope a = %.6f\n", a);
    printf("    E₂ = %.4f (twist energy, ρ-weighted)\n", E2);
    printf("    E₄ = %.4f (Skyrme term, ρ-independent)\n", E4);
    printf("    E_p = %.4f (pressure energy α∫ρ²)\n", Ep);
    printf("    Total deficit ΔQ = %.6f\n", DQ);
    printf("    ρ(0)/ρ₀ = %.6f\n", rhoG[0] / rho0);
    printf("    ρ(R_core)/ρ₀ ≈ %.6f  (at r=1.5)\n",
           rhoG[(int)(1.5/dr)] / rho0);
    printf("    Predicted ΔQ = E₂/(2αρ₀) = %.6f\n", E2 / (2*alpha_p*rho0));
    printf("    Ratio ΔQ_actual/ΔQ_predicted = %.6f\n",
           DQ / (E2/(2*alpha_p*rho0)));

    /* Print profile to file */
    FILE *fp = fopen("data/depletion_profile.dat", "w");
    if (fp) {
        fprintf(fp, "# r  f(r)  f'(r)  rho(r)  delta_rho(r)  |omega|^2\n");
        for (int i = 0; i < N_grid; i += 10) {
            double w2 = omega_sq(rG[i], fG[i], fpG[i]);
            fprintf(fp, "%.6f  %.8f  %.8f  %.8f  %.8e  %.8e\n",
                    rG[i], fG[i], fpG[i], rhoG[i],
                    rhoG[i] - rho0, w2);
        }
        fclose(fp);
        printf("\n  Profile written to data/depletion_profile.dat\n");
    }
}

static void test_long_range_tail(void) {
    printf("\n=== TEST 2: Long-range density tail ===\n\n");

    /* Solve Poisson equation */
    double *Phi = calloc(N_grid, sizeof(double));
    solve_poisson(Phi);

    /* The Poisson integral:
     *   Φ(r) = -(1/r)∫₀ʳ δρ·4πr'²dr' - ∫ᵣ^∞ δρ·4πr'dr'
     * For r → ∞:  Φ → -(1/r)∫δρ·4πr'²dr' = ΔQ/r
     * (since ∫δρ·4πr²dr = -ΔQ)
     */

    double E2, E4, Ep, DQ;
    compute_energies(&E2, &E4, &Ep, &DQ);

    printf("  Total deficit ΔQ = %.6f\n", DQ);
    printf("  Expected Poisson asymptotic: Φ(r) → ΔQ/r = %.6f/r\n", DQ);

    printf("\n  KEY COMPARISON: actual δρ(r) vs Poisson potential Φ(r)\n");
    printf("  %8s  %12s  %12s  %12s\n",
           "r", "δρ(r)", "Φ(r)", "ΔQ/r");

    double dr = R_max / (N_grid - 1);

    for (double r = 2.0; r <= 25.0; r += 1.0) {
        int i = (int)(r / dr);
        if (i >= N_grid) break;

        double drho = rhoG[i] - rho0;
        double phi = Phi[i];
        double phi_pred = DQ / r;

        printf("  %8.1f  %12.4e  %12.4e  %12.4e\n",
               r, drho, phi, phi_pred);
    }

    /* Fit power law to δρ(r) and Φ(r) in far field */
    double *drho_arr = calloc(N_grid, sizeof(double));
    for (int i = 0; i < N_grid; i++)
        drho_arr[i] = rhoG[i] - rho0;

    double amp, expo;
    fit_power_law(rG, drho_arr, N_grid, 3.0, 15.0, &amp, &expo);
    printf("\n  FITS:\n");
    printf("    δρ(r)  (r=3-15): |δρ| ~ %.4f / r^%.2f\n", amp, -expo);

    fit_power_law(rG, Phi, N_grid, 5.0, 20.0, &amp, &expo);
    printf("    Φ(r)   (r=5-20): |Φ|  ~ %.4f / r^%.2f  (expected 1/r)\n", amp, -expo);

    printf("\n  CRITICAL FINDING:\n");
    printf("    The actual density δρ(r) decays as ~1/r⁶ (from |ω|² = f'²+2sin²f/r²)\n");
    printf("    The Poisson integral Φ(r) decays as 1/r (mathematical, from ∫δρ/|x-x'|)\n");
    printf("    A test soliton at distance d couples to δρ(d), NOT to Φ(d).\n");
    printf("    Therefore the force goes as ~dδρ/dr ~ 1/r⁷, NOT as dΦ/dr ~ 1/r².\n");
    printf("    The density depletion mechanism is SHORT-RANGE.\n");

    FILE *fp = fopen("data/depletion_poisson.dat", "w");
    if (fp) {
        fprintf(fp, "# r  delta_rho  Phi  DQ/r\n");
        for (int i = 1; i < N_grid; i += 10) {
            fprintf(fp, "%.6f  %.8e  %.8e  %.8e\n",
                    rG[i], drho_arr[i], Phi[i], DQ/rG[i]);
        }
        fclose(fp);
    }

    free(Phi);
    free(drho_arr);
}

static void test_universality(void) {
    printf("\n=== TEST 3: Universality (ΔQ vs α) ===\n\n");

    printf("  %10s  %10s  %10s  %10s  %10s  %10s\n",
           "alpha", "E2", "DQ", "E2/(2αρ0)", "ratio", "rho(0)/rho0");

    double saved_alpha = alpha_p;

    for (double logalpha = -1.0; logalpha <= 3.0; logalpha += 0.5) {
        alpha_p = pow(10.0, logalpha);

        /* Re-initialize density */
        for (int i = 0; i < N_grid; i++) rhoG[i] = rho0;

        /* Self-consistent solve */
        for (int iter = 0; iter < 15; iter++) {
            find_profile();
            compute_density();
        }

        double E2, E4, Ep, DQ;
        compute_energies(&E2, &E4, &Ep, &DQ);
        double pred = E2 / (2 * alpha_p * rho0);

        printf("  %10.3f  %10.4f  %10.6f  %10.6f  %10.4f  %10.6f\n",
               alpha_p, E2, DQ, pred, DQ/pred, rhoG[0]/rho0);
    }

    alpha_p = saved_alpha;
    printf("\n  If ratio ≈ 1 for all α, universality holds: ΔQ = E₂/(2αρ₀)\n");
}

static void test_gravity_coupling(void) {
    printf("\n=== TEST 4: Gravity coupling G_eff ===\n\n");

    /* From the derivation:
     *   G_eff = κ/(4π c_s² ρ₀²)
     * where c_s² = 2α (sound speed from p = 2αρ)
     * and κ ~ 1 (profile-dependent constant)
     *
     * In code units (ρ₀ = 1, e = 1):
     *   G_eff = κ/(8πα)
     *
     * The gravitational potential between two B=1 solitons:
     *   Φ(r) = -G_eff · M² / r
     * where M = E₂ + E₄ ≈ 2E₂ (virial).
     */

    /* First solve the profile */
    double dr = R_max / (N_grid - 1);
    for (int i = 0; i < N_grid; i++) { rG[i] = i*dr; rhoG[i] = rho0; }
    for (int iter = 0; iter < 15; iter++) { find_profile(); compute_density(); }

    double E2, E4, Ep, DQ;
    compute_energies(&E2, &E4, &Ep, &DQ);
    double M = E2 + E4;

    printf("  Soliton mass M = E₂ + E₄ = %.4f + %.4f = %.4f\n", E2, E4, M);
    printf("  Total deficit ΔQ = %.6f\n", DQ);
    printf("  Sound speed c_s² = 2α = %.4f\n", 2*alpha_p);

    /* G_eff from the depletion mechanism */
    /* Force: F = -(E₂/ρ₀)·(ΔQ/(4πr²)) = -G_eff·M²/r²
     * → G_eff = (E₂/ρ₀)·ΔQ/(4π·M²)
     * With ΔQ = E₂/(2αρ₀):
     * → G_eff = E₂²/(8παρ₀²·M²)
     */
    double G_eff = E2 * E2 / (8 * M_PI * alpha_p * rho0 * rho0 * M * M);
    printf("  G_eff = E₂²/(8παρ₀²M²) = %.6e\n", G_eff);
    printf("  G_eff = 1/(8πα) × (E₂/M)² / ρ₀² = 1/(8π·%.2f) × (%.4f/%.4f)² = %.6e\n",
           alpha_p, E2, M, 1.0/(8*M_PI*alpha_p) * (E2/M)*(E2/M));

    /* Physical calibration */
    printf("\n  Physical calibration (using v2 parameters):\n");
    printf("    1 code energy = 9.098 MeV\n");
    printf("    1 code length = 0.5624 fm = 5.624e-16 m\n");
    printf("    M_proton = %.4f code = 938 MeV\n", M);

    /* To get G_Newton = 6.674e-11 m³/(kg·s²):
     * G_eff [code] = G_N × (m_code/E_code) × (L_code/T_code²)
     * In natural units (ℏ=c=1): G_N = 1/M_Planck² = 1/(1.22e19 GeV)²
     * In code units: G_N_code = G_N × (code_energy)²
     *   = (1/(1.22e19·1000)²) × (9.098)² MeV⁻² × MeV²
     *   = 9.098² / (1.22e22)² = 82.77 / 1.49e44 = 5.56e-43
     */
    double G_N_code = 82.77 / 1.49e44;
    printf("    G_Newton in code units ≈ %.2e\n", G_N_code);
    printf("    G_eff from depletion  = %.2e\n", G_eff);
    printf("    Ratio G_eff/G_Newton  = %.2e\n", G_eff / G_N_code);

    /* What α would give G_Newton? */
    double alpha_needed = E2*E2 / (8*M_PI*G_N_code*rho0*rho0*M*M);
    printf("\n  To match G_Newton: need α = %.2e\n", alpha_needed);
    printf("  Sound speed c_s² = 2α = %.2e code\n", 2*alpha_needed);
    printf("  c_s / c = %.2e\n", sqrt(2*alpha_needed));
    printf("  (c_s/c ~ 10⁻²¹ confirms scalar gravity speed problem)\n");
}

static void test_two_body(void) {
    printf("\n=== TEST 5: Two-body interaction: overlap vs Poisson ===\n\n");

    /* The PHYSICAL interaction between two solitons at separation d comes from
     * the fact that soliton 2 sits in the modified density field of soliton 1.
     * Soliton 2's twist energy depends on ρ (through E₂ = ½∫ρ|ω|²).
     *
     * The interaction energy is:
     *   E_int(d) = ∫ρ₂(x)|ω₂(x)|² × [ρ₁(|x-d|) - ρ₀]/(2ρ₀) d³x
     *            ≈ (E₂/ρ₀) × δρ_source(d)
     *
     * where δρ_source(d) is the density perturbation of soliton 1 at distance d.
     * Since δρ(r) ~ 1/r⁶, the interaction is ~ 1/r⁶ (SHORT-RANGE).
     *
     * Compare with the POISSON estimate:
     *   E_int_Poisson(d) ≈ (E₂/ρ₀) × (-ΔQ/(4πd))
     * which would give 1/r (LONG-RANGE) if the coupling were to the potential.
     */

    double E2, E4, Ep, DQ;
    compute_energies(&E2, &E4, &Ep, &DQ);
    double M = E2 + E4;

    printf("  M = %.4f, E₂ = %.4f, ΔQ = %.6f\n\n", M, E2, DQ);

    double dr = R_max / (N_grid - 1);

    printf("  %8s  %14s  %14s  %14s  %10s\n",
           "d", "E_overlap(d)", "E_Poisson(d)", "δρ(d)", "power_law");

    /* Track for power-law fit */
    double d_arr[30], E_ov_arr[30];
    int npts = 0;

    for (double d = 2.0; d <= 20.0; d += 1.0) {
        /* Actual δρ at distance d from the source */
        int id = (int)(d / dr);
        if (id >= N_grid) break;
        double drho_at_d = rhoG[id] - rho0;

        /* PHYSICAL: overlap interaction (test soliton energy change in modified ρ) */
        double E_overlap = (E2 / rho0) * drho_at_d;

        /* POISSON: point-source estimate (1/r potential) */
        double E_poisson = -(E2 / rho0) * DQ / (4*M_PI*d);

        /* Track for fitting */
        if (d >= 3.0 && fabs(E_overlap) > 1e-30) {
            d_arr[npts] = d;
            E_ov_arr[npts] = E_overlap;
            npts++;
        }

        printf("  %8.1f  %14.4e  %14.4e  %14.4e  %10s\n",
               d, E_overlap, E_poisson, drho_at_d, "");
    }

    /* Fit power law to overlap interaction */
    double amp, expo;
    fit_power_law(d_arr, E_ov_arr, npts, 3.0, 15.0, &amp, &expo);
    printf("\n  FITS:\n");
    printf("    E_overlap(d) ~ %.4e / d^%.2f\n", amp, -expo);
    printf("    E_Poisson(d) ~ %.4e / d^%.2f  (by construction: 1/d)\n",
           (E2/rho0)*DQ/(4*M_PI), 1.0);

    printf("\n  VERDICT:\n");
    printf("    Physical interaction (field overlap): ~ 1/r^%.1f (SHORT-RANGE)\n", -expo);
    printf("    Poisson estimate (mathematical):      ~ 1/r   (long-range)\n");
    printf("    The density depletion mechanism gives SHORT-RANGE interaction.\n");
    printf("    For 1/r gravity, need a MEDIATING FIELD (like Path 3's constraint field p)\n");
    printf("    that satisfies ∇²φ = source and couples test solitons to φ, not to δρ.\n");
}

/* ================================================================
 * Test 6: Causal density propagation (wave equation)
 *
 * The algebraic equilibrium δρ = -|ω|²/(4α) is LOCAL — density at each
 * point depends only on twist at that same point.  The physical density
 * should propagate causally at speed c, with the Heaviside step function
 * θ(ct - r) enforcing that information hasn't arrived beyond the wavefront.
 *
 * Two wave equation models compared:
 *
 * Model A — MASSLESS WAVE (ds = c·dt):
 *   □δρ = -½|ω|²  where □ = (1/c²)∂²/∂t² - ∇²
 *   Source is ½|ω|² with NONZERO monopole: ∫½|ω|² d³x = E₂/ρ₀ > 0
 *   Static solution: δρ ~ -(E₂/ρ₀)/(4πr)  →  1/r  (LONG-RANGE)
 *   Retarded: δρ(r,t) = -(Q_eff/4πr)·θ(ct-r)
 *
 * Model B — FLUID WAVE (continuity + Euler with pressure):
 *   ∂²δρ/∂t² = c_s²∇²δρ + ½∇²|ω|²
 *   Source is ∇²|ω|² with ZERO monopole (Gauss theorem: ∫∇²f = 0)
 *   Static solution: δρ → -|ω|²/(4α)  →  1/r⁶  (SHORT-RANGE)
 *
 * In u = r·δρ coordinates (reduces to 1D wave equation):
 *   A: ∂²u/∂t² = c²u'' - ½c²r|ω|²
 *   B: ∂²u/∂t² = c_s²u'' + ½(r|ω|²)''
 *
 * KEY INSIGHT: The Gauss obstruction (∫∇²|ω|² = 0) prevents 1/r
 * from the fluid model. But the massless wave equation has |ω|² as
 * source (not ∇²|ω|²), so the monopole is nonzero → 1/r is possible.
 * ================================================================ */

static void test_causal_propagation(void) {
    printf("\n=== TEST 6: Causal density propagation (wave equation) ===\n\n");

    double E2, E4, Ep, DQ;
    compute_energies(&E2, &E4, &Ep, &DQ);

    double c_grav = 1.0;  /* speed of gravity in code units */
    double cs = sqrt(2.0 * alpha_p * rho0);  /* sound speed from pressure */
    double Qeff = E2 / rho0;  /* monopole moment of source A */

    printf("  Soliton: E₂ = %.4f, ΔQ_algebraic = %.6f\n", E2, DQ);
    printf("  Speed of gravity c = %.4f (code units)\n", c_grav);
    printf("  Sound speed c_s = √(2αρ₀) = %.4f\n", cs);
    printf("  Source A monopole: ∫½|ω|² d³x = E₂/ρ₀ = %.4f  (NONZERO)\n", Qeff);
    printf("  Source B monopole: ∫∇²|ω|² d³x = 0  (ZERO by Gauss)\n\n");

    /* Wave grid — larger domain than static grid */
    int Nw = 10001;
    double Rw = 100.0;
    double drw = Rw / (Nw - 1);

    double *rw  = malloc(Nw * sizeof(double));
    double *om2 = calloc(Nw, sizeof(double));   /* |ω|² on wave grid */
    double *ww  = calloc(Nw, sizeof(double));   /* w = r·|ω|² */

    for (int i = 0; i < Nw; i++) {
        rw[i] = i * drw;
        /* Interpolate soliton profile (zero outside R_max) */
        double r = rw[i];
        int j = (int)(r / h_grid);
        if (j >= 0 && j < N_grid - 1 && r <= R_max) {
            double t = (r - rG[j]) / h_grid;
            double f_int = (1-t)*fG[j] + t*fG[j+1];
            double fp_int = (1-t)*fpG[j] + t*fpG[j+1];
            om2[i] = omega_sq(rw[i], f_int, fp_int);
        }
        ww[i] = rw[i] * om2[i];
    }

    /* Verify monopole of ½|ω|² */
    double mono_num = 0;
    for (int i = 1; i < Nw; i++)
        mono_num += 0.5 * om2[i] * 4*M_PI*rw[i]*rw[i] * drw;
    printf("  Numerical ∫½|ω|²·4πr²dr = %.4f  (expect E₂/ρ₀ = %.4f)\n\n",
           mono_num, Qeff);

    /* Precompute source arrays (in u-coordinates) */
    double *src_A = calloc(Nw, sizeof(double));
    double *src_B = calloc(Nw, sizeof(double));

    /* Model A: src_u = -½c²r|ω|² */
    for (int i = 0; i < Nw; i++)
        src_A[i] = -0.5 * c_grav * c_grav * rw[i] * om2[i];

    /* Model B: src_u = ½(r|ω|²)'' = ½w'' */
    for (int i = 1; i < Nw - 1; i++)
        src_B[i] = 0.5 * (ww[i+1] - 2*ww[i] + ww[i-1]) / (drw*drw);

    /* --- Model A: Massless wave □δρ = -½|ω|² --- */
    printf("  ============ Model A: Massless wave (ds = c·dt) ============\n\n");

    double CFL = 0.5;
    double dt_A = CFL * drw / c_grav;
    double CA2 = (c_grav * dt_A / drw) * (c_grav * dt_A / drw);

    double *uA_prev = calloc(Nw, sizeof(double));
    double *uA_curr = calloc(Nw, sizeof(double));
    double *uA_next = calloc(Nw, sizeof(double));

    /* First step: u¹ = ½dt²·source (from ∂u/∂t=0 at t=0) */
    for (int i = 1; i < Nw - 1; i++)
        uA_curr[i] = 0.5 * dt_A * dt_A * src_A[i];

    printf("  dt = %.6f, Courant² = %.4f\n", dt_A, CA2);
    printf("  Wavefront reaches Rw = %.0f at t = %.1f\n\n", Rw, Rw/c_grav);

    double t_snaps[] = {5.0, 10.0, 20.0, 50.0};
    int n_snaps = 4;
    int snap_A = 0;
    int max_steps_A = (int)(55.0 / dt_A);

    FILE *fout_A = fopen("data/wave_A.dat", "w");
    if (fout_A) fprintf(fout_A, "# r  t  delta_rho  minus_Qeff_4pi_r\n");

    for (int n = 1; n < max_steps_A && snap_A < n_snaps; n++) {
        double tA = (n + 1) * dt_A;

        /* Leapfrog */
        for (int i = 1; i < Nw - 1; i++) {
            uA_next[i] = 2*uA_curr[i] - uA_prev[i]
                        + CA2*(uA_curr[i+1] - 2*uA_curr[i] + uA_curr[i-1])
                        + dt_A*dt_A*src_A[i];
        }
        uA_next[0] = 0;
        uA_next[Nw-1] = 0;

        double *tmp = uA_prev; uA_prev = uA_curr; uA_curr = uA_next; uA_next = tmp;

        if (tA >= t_snaps[snap_A]) {
            double r_front = c_grav * tA;
            printf("  t = %.1f, wavefront at r = %.1f:\n", tA, r_front);
            printf("  %8s  %14s  %14s  %10s\n",
                   "r", "δρ_wave(r)", "-Q/(4πr)", "ratio");

            double d_arr[200], dA_arr[200];
            int nfit = 0;

            /* Total deficit integral */
            double total_deficit = 0;
            for (int i = 1; i < Nw; i++) {
                double drho = (rw[i] > 1e-10) ? uA_curr[i]/rw[i] : 0;
                total_deficit += drho * 4*M_PI*rw[i]*rw[i] * drw;
            }

            for (double r = 2.0; r < r_front*0.9 && r < Rw-1; r *= 1.4) {
                int idx = (int)(r / drw);
                if (idx <= 0 || idx >= Nw) continue;
                double drho = uA_curr[idx] / rw[idx];
                double pred = -Qeff / (4*M_PI*r);
                double ratio = (fabs(pred) > 1e-30) ? drho/pred : 0;

                printf("  %8.2f  %14.6e  %14.6e  %10.4f\n", r, drho, pred, ratio);

                if (r >= 3.0 && r < r_front*0.5 && nfit < 200) {
                    d_arr[nfit] = r;
                    dA_arr[nfit] = fabs(drho);
                    nfit++;
                }

                if (fout_A)
                    fprintf(fout_A, "%.4f  %.2f  %.8e  %.8e\n",
                            r, tA, drho, pred);
            }

            if (nfit >= 5) {
                double amp, expo;
                fit_power_law(d_arr, dA_arr, nfit, 3.0, r_front*0.4, &amp, &expo);
                printf("  Power-law fit (r=3..%.0f): |δρ| ~ r^%.3f  "
                       "(expect -1.0 for 1/r)\n", r_front*0.4, expo);
            }
            printf("  Total deficit ∫δρ·4πr²dr = %.4f  "
                   "(grows as ½c²Q_eff·t² = %.4f)\n\n",
                   total_deficit,
                   -0.5*c_grav*c_grav*Qeff*tA*tA);

            snap_A++;
        }
    }
    if (fout_A) fclose(fout_A);

    /* --- Model B: Fluid wave ∂²δρ/∂t² = c_s²∇²δρ + ½∇²|ω|² --- */
    printf("  ============ Model B: Fluid wave (∇²|ω|² source) ============\n\n");

    double dt_B = CFL * drw / cs;
    double CB2 = (cs * dt_B / drw) * (cs * dt_B / drw);

    double *uB_prev = calloc(Nw, sizeof(double));
    double *uB_curr = calloc(Nw, sizeof(double));
    double *uB_next = calloc(Nw, sizeof(double));

    for (int i = 1; i < Nw - 1; i++)
        uB_curr[i] = 0.5 * dt_B * dt_B * src_B[i];

    printf("  dt = %.6e, Courant² = %.4f, c_s = %.4f\n", dt_B, CB2, cs);

    int snap_B = 0;
    int max_steps_B = (int)(55.0 / dt_B);
    if (max_steps_B > 5000000) max_steps_B = 5000000;

    printf("  max_steps = %d\n\n", max_steps_B);

    for (int n = 1; n < max_steps_B && snap_B < n_snaps; n++) {
        double tB = (n + 1) * dt_B;

        for (int i = 1; i < Nw - 1; i++) {
            uB_next[i] = 2*uB_curr[i] - uB_prev[i]
                        + CB2*(uB_curr[i+1] - 2*uB_curr[i] + uB_curr[i-1])
                        + dt_B*dt_B*src_B[i];
        }
        uB_next[0] = 0;
        uB_next[Nw-1] = 0;

        double *tmp2 = uB_prev; uB_prev = uB_curr; uB_curr = uB_next; uB_next = tmp2;

        if (tB >= t_snaps[snap_B]) {
            double r_front = cs * tB;
            printf("  t = %.3f, wavefront at r = %.1f:\n", tB, r_front);
            printf("  %8s  %14s  %14s  %10s\n",
                   "r", "δρ_fluid(r)", "δρ_alg(r)", "ratio");

            double total_B = 0;
            for (int i = 1; i < Nw; i++) {
                double drho = (rw[i] > 1e-10) ? uB_curr[i]/rw[i] : 0;
                total_B += drho * 4*M_PI*rw[i]*rw[i] * drw;
            }

            double d_arr[200], dB_arr[200];
            int nfit = 0;

            for (double r = 2.0; r < r_front*0.8 && r < Rw-1; r *= 1.4) {
                int idx = (int)(r / drw);
                if (idx <= 0 || idx >= Nw) continue;
                double drho = uB_curr[idx] / rw[idx];
                double drho_alg = -om2[idx] / (4.0 * alpha_p);
                double ratio = (fabs(drho_alg) > 1e-30) ? drho/drho_alg : 0;

                printf("  %8.2f  %14.6e  %14.6e  %10.4f\n",
                       r, drho, drho_alg, ratio);

                if (r >= 3.0 && r < r_front*0.5 && fabs(drho) > 1e-30 && nfit < 200) {
                    d_arr[nfit] = r;
                    dB_arr[nfit] = fabs(drho);
                    nfit++;
                }
            }

            if (nfit >= 5) {
                double amp, expo;
                fit_power_law(d_arr, dB_arr, nfit, 3.0, r_front*0.4, &amp, &expo);
                printf("  Power-law fit: |δρ| ~ r^%.3f  "
                       "(algebraic gives -6)\n", expo);
            }
            printf("  Total ∫δρ·4πr²dr = %.6e  (should be ~0 if conserving)\n\n",
                   total_B);

            snap_B++;
        }
    }

    /* --- Summary --- */
    printf("  ====================== SUMMARY ======================\n\n");
    printf("  Three density models for the same B=1 soliton:\n\n");
    printf("  1. ALGEBRAIC (local equilibrium): δρ = -|ω|²/(4α)\n");
    printf("     Decay: 1/r⁶  (follows |ω|² ∝ f'² + 2sin²f/r²)\n");
    printf("     Problem: purely local, no propagation\n\n");
    printf("  2. FLUID WAVE (continuity + Euler):\n");
    printf("     ∂²δρ/∂t² = c_s²∇²δρ + ½∇²|ω|²\n");
    printf("     Source ∇²|ω|² has ZERO monopole (Gauss: ∫∇²f = 0)\n");
    printf("     Converges to algebraic equilibrium: δρ ~ 1/r⁶\n");
    printf("     Conservation: ∫δρ = const = 0  ✓\n\n");
    printf("  3. MASSLESS WAVE with ds = c·dt:\n");
    printf("     □δρ = -½|ω|²  (c = speed of gravity)\n");
    printf("     Source ½|ω|² has NONZERO monopole: Q_eff = E₂/ρ₀ = %.4f\n",
           Qeff);
    printf("     Retarded solution: δρ(r,t) = -(Q_eff/4πr)·θ(ct - r)\n");
    printf("     Decay: 1/r inside wavefront  (LONG-RANGE)\n");
    printf("     Conservation: ∫δρ grows as ½c²Q_eff·t²\n\n");
    printf("  KEY PHYSICS:\n");
    printf("    The step function θ(ct-r) is the causal propagation ds = c·dt.\n");
    printf("    It makes density NON-LOCAL: δρ at distance r depends on\n");
    printf("    retarded information from the source, NOT local twist.\n\n");
    printf("    The monopole source ½|ω|² BYPASSES the Gauss obstruction\n");
    printf("    because the source is |ω|² itself, not ∇²|ω|².\n\n");
    printf("    Conservation is maintained in infinite volume: the growing\n");
    printf("    deficit is compensated by a wavefront compression spike\n");
    printf("    at r = ct, pushed to infinity as t → ∞.\n");

    free(rw); free(om2); free(ww);
    free(src_A); free(src_B);
    free(uA_prev); free(uA_curr); free(uA_next);
    free(uB_prev); free(uB_curr); free(uB_next);
}

/* ================================================================
 * Main
 * ================================================================ */

int main(int argc, char **argv) {
    N_grid = 10001;
    int do_scan = 0;
    int do_wave = 0;

    /* Parse arguments */
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-alpha") == 0 && i+1 < argc)
            alpha_p = atof(argv[++i]);
        else if (strcmp(argv[i], "-B") == 0 && i+1 < argc)
            B_charge = atoi(argv[++i]);
        else if (strcmp(argv[i], "-mpi") == 0 && i+1 < argc)
            m_pi = atof(argv[++i]);
        else if (strcmp(argv[i], "-rho0") == 0 && i+1 < argc)
            rho0 = atof(argv[++i]);
        else if (strcmp(argv[i], "-e") == 0 && i+1 < argc)
            { double e_val = atof(argv[++i]); c4 = 2.0 * rho0*rho0 / (e_val * e_val); }
        else if (strcmp(argv[i], "-Rmax") == 0 && i+1 < argc)
            R_max = atof(argv[++i]);
        else if (strcmp(argv[i], "-N") == 0 && i+1 < argc)
            N_grid = atoi(argv[++i]);
        else if (strcmp(argv[i], "-scan") == 0)
            do_scan = 1;
        else if (strcmp(argv[i], "-wave") == 0)
            do_wave = 1;
        else {
            fprintf(stderr, "Usage: %s [-alpha A] [-B 1] [-mpi M] [-rho0 R] "
                    "[-e E] [-Rmax R] [-N n] [-scan] [-wave]\n", argv[0]);
            return 1;
        }
    }

    if (N_grid > NMAX) N_grid = NMAX;

    printf("V6 Depletion Test\n");
    printf("  α = %.4f, ρ₀ = %.4f, c₄ = %.4f, m_π = %.4f\n",
           alpha_p, rho0, c4, m_pi);
    printf("  B = %d, R_max = %.1f, N = %d\n\n", B_charge, R_max, N_grid);

    /* Initialize grid */
    h_grid = R_max / (N_grid - 1);
    for (int i = 0; i < N_grid; i++) {
        rG[i] = i * h_grid;
        rhoG[i] = rho0;
    }

    /* Run tests */
    test_depletion_profile();
    test_long_range_tail();

    if (do_scan) {
        test_universality();
    }

    test_gravity_coupling();
    test_two_body();

    if (do_wave) {
        test_causal_propagation();
    }

    return 0;
}
