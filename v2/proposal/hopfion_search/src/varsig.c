/*
 * varsig.c — Self-Referential Mapping: Cl(3,0,1) ↔ SO(3,1)
 *
 * Correct self-referential mapping: constant ε₀² (NOT position-dependent).
 * At e₀² = -ε₀², the norm constraint is |q|² + ε₀²p² = ρ₀² (with j=0).
 * Through the constraint, the Skyrmion strain S₀(r) acts as an attractive
 * potential for the pseudoscalar p.
 *
 * Energy perturbation at small p:
 *   δE = (ε₀²/2) ∫[|∇p|² - S₀(r)p²] dV
 *
 * Key: E₄ is ρ-independent ([L_i,L_j] = [l̂_i,l̂_j], ∂ρ/ρ commutes with all).
 * So the effective potential comes ONLY from E₂.
 *
 * Eigenvalue problem:
 *   -∇²p - S₀(r)p = E₀p
 *
 * If E₀ < 0: Skyrmion is UNSTABLE to spontaneous pseudoscalar generation.
 * The pseudoscalar halo p(r) has range 1/κ where κ = √|E₀|.
 *
 * Using u(r) = rp(r) (s-wave, ℓ=0):
 *   -u'' - S₀(r)u = E₀u,  u(0) = 0
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ========== Profile I/O ========== */

typedef struct {
    int n;
    double *r, *f, *fp;
    double dr;
} Profile;

static int read_profile(const char *fname, Profile *prof) {
    FILE *fp = fopen(fname, "r");
    if (!fp) { fprintf(stderr, "Cannot open %s\n", fname); return -1; }

    int cap = 4096;
    prof->r  = malloc(cap * sizeof(double));
    prof->f  = malloc(cap * sizeof(double));
    prof->fp = malloc(cap * sizeof(double));
    int n = 0;
    char line[512];

    while (fgets(line, sizeof(line), fp)) {
        if (line[0] == '#') continue;
        double rv, fval, fpval = 0;
        double d1, d2, d3;
        int nc = sscanf(line, "%lf %lf %lf %lf %lf %lf",
                        &rv, &fval, &fpval, &d1, &d2, &d3);
        if (nc < 2) continue;
        if (n >= cap) {
            cap *= 2;
            prof->r  = realloc(prof->r,  cap * sizeof(double));
            prof->f  = realloc(prof->f,  cap * sizeof(double));
            prof->fp = realloc(prof->fp, cap * sizeof(double));
        }
        prof->r[n] = rv;
        prof->f[n] = fval;
        prof->fp[n] = fpval;
        n++;
    }
    fclose(fp);
    prof->n = n;
    prof->dr = (n > 1) ? prof->r[1] - prof->r[0] : 0.001;

    /* Check if f' was provided */
    int has_fp = 0;
    for (int i = 0; i < n; i++)
        if (fabs(prof->fp[i]) > 1e-30) { has_fp = 1; break; }

    if (!has_fp) {
        fprintf(stderr, "Computing f'(r) numerically...\n");
        for (int i = 0; i < n; i++) {
            if (i >= 2 && i < n-2)
                prof->fp[i] = (-prof->f[i+2] + 8*prof->f[i+1]
                               - 8*prof->f[i-1] + prof->f[i-2]) / (12*prof->dr);
            else if (i == 0)
                prof->fp[i] = (-3*prof->f[0] + 4*prof->f[1] - prof->f[2]) / (2*prof->dr);
            else if (i == n-1)
                prof->fp[i] = (prof->f[n-3] - 4*prof->f[n-2] + 3*prof->f[n-1]) / (2*prof->dr);
            else
                prof->fp[i] = (prof->f[i+1] - prof->f[i-1]) / (2*prof->dr);
        }
    }

    return 0;
}

static void free_profile(Profile *p) {
    free(p->r); free(p->f); free(p->fp);
}

/* ========== Physical quantities ========== */

/* Skyrmion strain: S₀(r) = f'² + 2sin²f/r² */
static double strain_S0(double r, double f, double fp, double a) {
    if (r > 1e-10) {
        double sf = sin(f);
        return fp * fp + 2.0 * sf * sf / (r * r);
    } else {
        return 3.0 * a * a;  /* L'Hôpital limit */
    }
}

/* Skyrme angular density: T₀(r) = 2f'²sin²f/r² + sin⁴f/r⁴ */
static double skyrme_T0(double r, double f, double fp, double a) {
    if (r > 1e-10) {
        double sf = sin(f);
        double sf2 = sf * sf;
        return 2.0 * fp * fp * sf2 / (r * r) + sf2 * sf2 / (r * r * r * r);
    } else {
        return 3.0 * a * a * a * a;  /* Both terms → a⁴ limit */
    }
}

/* Baryon density: B⁰(r) = -f'sin²f/(2π²r²) */
static double baryon_density(double r, double f, double fp, double a) {
    if (r > 1e-10) {
        double sf = sin(f);
        return -(1.0/(2*M_PI*M_PI)) * sf*sf * fp / (r*r);
    } else {
        return a*a*a / (2*M_PI*M_PI);
    }
}

/* ========== Eigenvalue solver ========== */

/*
 * Solve: -u'' - W(r)u = E₀u  where u(r) = rp(r), u(0) = 0
 *
 * W(r) is the attractive potential (S₀ or S₀ + correction).
 *
 * Method: Numerov integration from r=0 outward.
 * For u'' = -[W(r) + E₀]u, Numerov gives:
 *   u_{n+1} = [2u_n(1 - 5h²k_n²/12) - u_{n-1}(1 + h²k_{n-1}²/12)]
 *             / (1 + h²k_{n+1}²/12)
 * where k²(r) = W(r) + E₀.
 */
static int solve_eigenfunction(const double *W, const double *r_grid, int n,
                                double dr, double E0, double *u)
{
    /* Start: u(0) = 0, u(dr) = dr (normalization) */
    u[0] = 0.0;
    u[1] = dr;

    int nodes = 0;
    double h2 = dr * dr;

    for (int i = 1; i < n - 1; i++) {
        double k2_prev = W[i-1] + E0;
        double k2_curr = W[i]   + E0;
        double k2_next = W[i+1] + E0;

        double numer = 2.0 * u[i] * (1.0 - 5.0*h2*k2_curr/12.0)
                     - u[i-1] * (1.0 + h2*k2_prev/12.0);
        double denom = 1.0 + h2*k2_next/12.0;

        u[i+1] = numer / denom;

        /* Count sign changes (nodes) */
        if (u[i+1] * u[i] < 0) nodes++;
    }

    return nodes;
}

/* Find eigenvalue by bisection on node count */
static double find_eigenvalue(const double *W, const double *r_grid, int n,
                               double dr, int target_nodes,
                               double E_lo, double E_hi, double tol)
{
    double *u = calloc(n, sizeof(double));

    for (int iter = 0; iter < 200; iter++) {
        double E_mid = 0.5 * (E_lo + E_hi);
        int nodes = solve_eigenfunction(W, r_grid, n, dr, E_mid, u);

        if (nodes > target_nodes)
            E_hi = E_mid;  /* too many nodes → we're above eigenvalue */
        else
            E_lo = E_mid;  /* too few nodes → we're below eigenvalue */

        if (fabs(E_hi - E_lo) < tol) break;
    }

    free(u);
    return 0.5 * (E_lo + E_hi);
}

/* ========== Self-consistent solver ========== */

/*
 * At given ε₀², find the self-consistent p(r) profile.
 *
 * The equilibrium p(r) is the eigenfunction of the LOWEST eigenvalue E₀
 * of -∇²p - S₀p = E₀p, with amplitude set by the nonlinear constraint:
 *
 *   |q|² = ρ₀² - ε₀²p²  →  p² < ρ₀²/ε₀²
 *
 * For small ε₀², the amplitude is small and p₀(r) ≈ (eigenfunction) × small const.
 * The energy correction is δE = ε₀² E₀ × ∫p²dV (negative if E₀ < 0).
 */

/* ========== Main ========== */

int main(int argc, char *argv[])
{
    const char *profile_file = "data/profiles/profile_sigma_e1.dat";
    double e_param = 1.0;
    double rho0 = 1.0;
    int verbose = 0;

    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "-profile") && i+1 < argc) profile_file = argv[++i];
        else if (!strcmp(argv[i], "-e") && i+1 < argc) e_param = atof(argv[++i]);
        else if (!strcmp(argv[i], "-rho0") && i+1 < argc) rho0 = atof(argv[++i]);
        else if (!strcmp(argv[i], "-v")) verbose = 1;
        else {
            fprintf(stderr, "Usage: %s [-profile PATH] [-e E] [-rho0 R] [-v]\n", argv[0]);
            return 1;
        }
    }

    double c4 = 2.0 * rho0 * rho0 / (e_param * e_param);

    /* Read profile */
    Profile prof;
    if (read_profile(profile_file, &prof) < 0) return 1;
    int n = prof.n;
    double dr = prof.dr;

    printf("============================================================\n");
    printf(" Self-Referential Mapping: Cl(3,0,1) <-> SO(3,1)\n");
    printf(" Constraint Eigenvalue Problem\n");
    printf("============================================================\n\n");
    printf("Profile: %s (%d points, R_max=%.2f, dr=%.6f)\n",
           profile_file, n, prof.r[n-1], dr);
    printf("Parameters: e=%.4f, rho0=%.4f, c4=%.4f\n\n", e_param, rho0, c4);

    double a = -prof.fp[0];
    printf("Shooting parameter: a = -f'(0) = %.6f\n", a);
    printf("f(0) = %.6f (expect pi = %.6f)\n\n", prof.f[0], M_PI);

    /* ===== Section 1: Compute potential W(r) ===== */
    printf("--- Constraint-Mediated Potential ---\n\n");
    printf("At e0^2 = -eps0^2, norm constraint: |q|^2 + eps0^2*p^2 = rho0^2\n");
    printf("Energy perturbation: dE = (eps0^2/2) int[|grad p|^2 - W(r)p^2] dV\n");
    printf("E4 is rho-independent ([L_i,L_j]=[l_i,l_j], d(rho)/rho commutes)\n");
    printf("So W(r) = S0(r) from E2 only.\n\n");

    double *S0 = calloc(n, sizeof(double));
    double *T0 = calloc(n, sizeof(double));
    double *B0 = calloc(n, sizeof(double));

    for (int i = 0; i < n; i++) {
        S0[i] = strain_S0(prof.r[i], prof.f[i], prof.fp[i], a);
        T0[i] = skyrme_T0(prof.r[i], prof.f[i], prof.fp[i], a);
        B0[i] = baryon_density(prof.r[i], prof.f[i], prof.fp[i], a);
    }

    /* Compute E₂ and E₄ for verification */
    double E2 = 0, E4 = 0;
    for (int i = 1; i < n; i++) {
        double r0 = prof.r[i-1], r1 = prof.r[i];
        E2 += 0.5*dr * (2*M_PI*rho0*rho0*S0[i-1]*r0*r0
                       + 2*M_PI*rho0*rho0*S0[i]*r1*r1);
        E4 += 0.5*dr * (4*M_PI*(rho0*rho0*rho0*rho0/(e_param*e_param))*T0[i-1]*r0*r0
                       + 4*M_PI*(rho0*rho0*rho0*rho0/(e_param*e_param))*T0[i]*r1*r1);
    }
    double E_sol = E2 + E4;

    printf("E2 = %.4f, E4 = %.4f, E_sol = %.4f\n", E2, E4, E_sol);
    printf("Virial: E2/E4 = %.6f (expect 1.000)\n\n", E2/E4);

    printf("Potential W(r) = S0(r) = f'^2 + 2sin^2(f)/r^2:\n");
    printf("  W(0) = 3a^2 = %.6f\n", 3*a*a);
    printf("  Depth x range^2 estimate: W(0)*R_core^2 ~ %.1f >> pi^2/4 = %.2f\n",
           3*a*a * 4.0, M_PI*M_PI/4);
    printf("  => Bound states expected (deep potential well)\n\n");

    if (verbose) {
        printf("  %8s  %12s  %12s  %12s  %12s\n",
               "r", "f(r)", "S0(r)", "T0(r)", "B0(r)");
        for (int i = 0; i < n; i += (n > 200 ? n/40 : 1)) {
            printf("  %8.4f  %12.6f  %12.6f  %12.6f  %12.6e\n",
                   prof.r[i], prof.f[i], S0[i], T0[i], B0[i]);
        }
        printf("\n");
    }

    /* ===== Section 2: Eigenvalue scan ===== */
    printf("--- Eigenvalue Scan: -u'' - S0(r)u = E0*u ---\n\n");

    /* First, scan E₀ from very negative to zero, count nodes */
    printf("Node count vs E0 (u = r*p, Numerov method):\n");
    printf("  %12s  %8s\n", "E0", "nodes");

    double *u_tmp = calloc(n, sizeof(double));

    double E_scan_vals[] = {-10.0, -8.0, -6.0, -5.0, -4.0, -3.5, -3.0,
                            -2.5, -2.0, -1.5, -1.0, -0.5, -0.2, -0.1, -0.01, 0.0};
    int n_scan = sizeof(E_scan_vals)/sizeof(E_scan_vals[0]);

    for (int ie = 0; ie < n_scan; ie++) {
        int nodes = solve_eigenfunction(S0, prof.r, n, dr, E_scan_vals[ie], u_tmp);
        printf("  %12.4f  %8d\n", E_scan_vals[ie], nodes);
    }

    /* ===== Section 3: Find ground state eigenvalue ===== */
    printf("\n--- Ground State (0 nodes) ---\n\n");

    /* The ground state has 0 nodes. Find E₀ by bisection.
     * E₀ is at the transition: nodes goes from 0 (below E₀) to 1 (above E₀). */

    double E_lo = -10.0, E_hi = 0.0;

    /* Scan upward, find where node count transitions from 0 to ≥1 */
    {
        int prev_nodes = solve_eigenfunction(S0, prof.r, n, dr, -10.0, u_tmp);
        for (double E_try = -9.9; E_try <= 0.01; E_try += 0.1) {
            int nodes = solve_eigenfunction(S0, prof.r, n, dr, E_try, u_tmp);
            if (prev_nodes == 0 && nodes >= 1) {
                E_lo = E_try - 0.1;  /* 0 nodes: below eigenvalue */
                E_hi = E_try;        /* 1 node:  above eigenvalue */
                break;
            }
            prev_nodes = nodes;
        }
    }

    /* Bisect for ground state: 0 nodes */
    double E_ground = find_eigenvalue(S0, prof.r, n, dr, 0, E_lo, E_hi, 1e-10);

    printf("Ground state eigenvalue: E0 = %.10f\n", E_ground);
    printf("Decay length: 1/kappa = 1/sqrt(|E0|) = %.6f code lengths = %.4f fm\n",
           1.0/sqrt(fabs(E_ground)), 1.0/sqrt(fabs(E_ground)) * 0.5624);

    if (E_ground < 0) {
        printf("\n*** E0 < 0: INSTABILITY! ***\n");
        printf("The Skyrmion is unstable to spontaneous pseudoscalar generation.\n");
        printf("At ANY eps0^2 > 0, the soliton develops a p(r) halo.\n");
    } else {
        printf("\n*** E0 >= 0: STABLE ***\n");
        printf("The p = 0 solution is stable. No spontaneous pseudoscalar.\n");
    }

    /* Compute and normalize the ground state eigenfunction */
    solve_eigenfunction(S0, prof.r, n, dr, E_ground, u_tmp);

    /* Convert u → p = u/r and normalize: ∫p² 4πr² dr = 1 */
    double *p_eigen = calloc(n, sizeof(double));
    p_eigen[0] = u_tmp[1] / prof.r[1];  /* extrapolate from i=1 */
    for (int i = 1; i < n; i++)
        p_eigen[i] = u_tmp[i] / prof.r[i];

    double norm2 = 0;
    for (int i = 1; i < n; i++) {
        double r0 = prof.r[i-1], r1 = prof.r[i];
        norm2 += 0.5*dr * (4*M_PI*p_eigen[i-1]*p_eigen[i-1]*r0*r0
                          + 4*M_PI*p_eigen[i]*p_eigen[i]*r1*r1);
    }

    if (norm2 > 0) {
        double norm_fac = 1.0 / sqrt(norm2);
        for (int i = 0; i < n; i++)
            p_eigen[i] *= norm_fac;
    }

    printf("\nNormalized eigenfunction p0(r) [int p0^2 4pi r^2 dr = 1]:\n");
    printf("  %8s  %14s  %14s  %14s\n", "r", "p0(r)", "p0^2*4pi*r^2", "S0(r)");
    for (int i = 0; i < n; i += (n > 200 ? n/30 : 1)) {
        double r = prof.r[i];
        printf("  %8.4f  %14.6e  %14.6e  %14.6e\n",
               r, p_eigen[i],
               (r > 1e-10) ? p_eigen[i]*p_eigen[i]*4*M_PI*r*r : 0,
               S0[i]);
    }

    /* ===== Section 4: Find excited states ===== */
    printf("\n--- Excited States ---\n\n");

    /* Find first few excited states (1, 2, 3 nodes) */
    for (int target = 1; target <= 5; target++) {
        /* Scan for bracket: transition from target nodes to target+1 nodes */
        double e_lo = -10.0, e_hi = 0.0;
        int found = 0;

        int prev_nodes = solve_eigenfunction(S0, prof.r, n, dr, -10.0, u_tmp);
        for (double E_try = -9.99; E_try <= 0.01; E_try += 0.01) {
            int nodes = solve_eigenfunction(S0, prof.r, n, dr, E_try, u_tmp);
            if (prev_nodes <= target && nodes > target) {
                e_lo = E_try - 0.01;  /* target nodes: below eigenvalue */
                e_hi = E_try;         /* target+1 nodes: above eigenvalue */
                found = 1;
                break;
            }
            prev_nodes = nodes;
        }

        if (found) {
            double E_n = find_eigenvalue(S0, prof.r, n, dr, target, e_lo, e_hi, 1e-10);
            printf("  n=%d (%d nodes): E_%d = %.10f, kappa = %.6f, range = %.4f\n",
                   target, target, target, E_n,
                   (E_n < 0) ? sqrt(fabs(E_n)) : 0,
                   (E_n < 0) ? 1.0/sqrt(fabs(E_n)) : 1e10);
        } else {
            printf("  n=%d: no bound state found (E > 0)\n", target);
            break;
        }
    }

    /* ===== Section 5: Physical consequences ===== */
    printf("\n--- Physical Consequences ---\n\n");

    double kappa_ground = (E_ground < 0) ? sqrt(fabs(E_ground)) : 0;
    double range_ground = (kappa_ground > 0) ? 1.0/kappa_ground : 1e10;

    printf("Self-referential mapping: constant eps0^2 with constraint coupling\n\n");

    if (E_ground < 0) {
        printf("1. INSTABILITY: E0 = %.6f < 0\n", E_ground);
        printf("   The Skyrmion spontaneously generates a pseudoscalar halo p(r).\n");
        printf("   This is the CONSTRAINT-MEDIATED coupling through |q|^2 + eps0^2*p^2 = rho0^2.\n\n");

        printf("2. Halo profile: p(r) ~ p0(r) [eigenfunction], decays as e^{-%.4f r}/r\n",
               kappa_ground);
        printf("   Range: %.4f code = %.4f fm (%.1f%% of R_rms=1.50)\n",
               range_ground, range_ground * 0.5624,
               range_ground/1.50*100);
        printf("   This is SHORT-RANGE (exponential, not 1/r)\n\n");

        printf("3. Energy correction:\n");
        printf("   dE = eps0^2 * E0 * int[p0^2] dV = eps0^2 * (%.6f)\n", E_ground);
        printf("   Per unit eps0^2: dE/eps0^2 = %.6f (negative = attractive)\n\n", E_ground);

        printf("4. Effective ρ reduction at core:\n");
        printf("   |q(0)|^2 = rho0^2 - eps0^2 * p0(0)^2\n");
        printf("   p0(0) = %.6e (normalized), so delta_rho/rho ~ eps0^2 * p0(0)^2/(2*rho0^2)\n\n",
               p_eigen[0]);

        /* Compare with baryon density */
        printf("5. Source structure comparison:\n");
        printf("   %8s  %14s  %14s  %14s\n",
               "r", "p0^2*r^2", "B0*r^2", "S0*r^2");
        for (int i = 0; i < n; i += (n > 200 ? n/20 : 1)) {
            double r = prof.r[i];
            if (r < 0.01) continue;
            printf("   %8.4f  %14.6e  %14.6e  %14.6e\n",
                   r, p_eigen[i]*p_eigen[i]*r*r,
                   B0[i]*r*r, S0[i]*r*r);
        }

        printf("\n6. Two-soliton interaction from pseudoscalar exchange:\n");
        printf("   Since p(r) ~ e^{-kappa*r}/r (exponential tail):\n");
        printf("   U(D) ~ e^{-kappa*D}/D  (Yukawa, NOT 1/r)\n");
        printf("   kappa = %.6f => range = %.4f code = %.4f fm\n",
               kappa_ground, range_ground, range_ground * 0.5624);
        printf("   This is a SHORT-RANGE nuclear-type force, not gravity.\n");

    } else {
        printf("1. STABLE: E0 = %.6f >= 0\n", E_ground);
        printf("   p = 0 is an energy minimum. No spontaneous pseudoscalar.\n");
        printf("   The self-referential mapping produces NO new physics.\n");
    }

    /* ===== Section 6: Comparison with Path 3 ===== */
    printf("\n--- Comparison with Path 3 (B0*p coupling) ---\n\n");

    printf("Path 3 (degenerate.c): -kappa^2 nabla^2 p = g_top B0\n");
    printf("  Source: B0(r) = baryon density (external coupling g_top)\n");
    printf("  Solution: p ~ g_top/(4*pi*kappa^2*r) (1/r, Newtonian)\n");
    printf("  Range: infinite (massless) or Yukawa (massive)\n\n");

    printf("Self-referential (this): -nabla^2 p - S0(r)p = E0*p\n");
    printf("  Source: NONE (eigenvalue problem, not sourced equation)\n");
    printf("  Solution: p ~ e^{-kappa*r}/r (bound state, localized)\n");
    printf("  Range: 1/kappa = %.4f (finite)\n\n", range_ground);

    printf("KEY DIFFERENCE: Path 3 has a SOURCE (B0, giving 1/r far field).\n");
    printf("The self-referential mapping has no source term — the constraint\n");
    printf("coupling is an EIGENVALUE problem, not a Poisson equation.\n");
    printf("The resulting p-halo is localized (Yukawa-like), not 1/r.\n\n");

    printf("To get 1/r gravity, a LINEAR source term for p is required.\n");
    printf("This cannot come from:\n");
    printf("  - L2 scalar part: cross-terms <(dq)(d~d)>_0 = 0 (grade mismatch)\n");
    printf("  - L4 (Skyrme): first-order coupling vanishes (dF_{ij} is grade 4)\n");
    printf("  - Constraint: gives eigenvalue equation, not sourced equation\n");
    printf("  - WZW: couples B^mu*w_mu, vanishes on statics\n");
    printf("The coupling constant g_top remains a FREE PARAMETER.\n");

    /* ===== Section 7: Verify with ODE residual ===== */
    printf("\n--- Verification: ODE residual ---\n\n");

    /* Check that p_eigen satisfies -p'' - 2p'/r - S₀p = E₀p */
    double max_resid = 0, avg_resid = 0;
    int n_resid = 0;
    for (int i = 2; i < n-2; i++) {
        double r = prof.r[i];
        if (r < 2*dr) continue;  /* skip r~0 where numerical derivatives are poor */

        double pp = (-p_eigen[i+2] + 16*p_eigen[i+1] - 30*p_eigen[i]
                     + 16*p_eigen[i-1] - p_eigen[i-2]) / (12*dr*dr);
        double pv = (-p_eigen[i+2] + 8*p_eigen[i+1]
                     - 8*p_eigen[i-1] + p_eigen[i-2]) / (12*dr);

        double lhs = -pp - 2*pv/r - S0[i]*p_eigen[i];
        double rhs = E_ground * p_eigen[i];
        double resid = fabs(lhs - rhs);
        double scale = fabs(rhs) + 1e-30;

        if (resid/scale > max_resid) max_resid = resid/scale;
        avg_resid += resid/scale;
        n_resid++;
    }
    avg_resid /= (n_resid > 0 ? n_resid : 1);

    printf("ODE: -p'' - 2p'/r - S0*p = E0*p\n");
    printf("Max relative residual: %.6e\n", max_resid);
    printf("Avg relative residual: %.6e\n", avg_resid);

    /* ===== Final summary ===== */
    printf("\n============================================================\n");
    printf(" Summary: Self-Referential Mapping Cl(3,0,1) <-> SO(3,1)\n");
    printf("============================================================\n\n");

    printf("The correct self-referential mapping uses CONSTANT eps0^2\n");
    printf("(not position-dependent). At e0^2 = -eps0^2, the norm constraint\n");
    printf("|q|^2 + eps0^2*p^2 = rho0^2 couples the pseudoscalar p to the\n");
    printf("Skyrmion strain S0(r) through the E2 energy functional.\n\n");

    printf("Eigenvalue problem: -nabla^2 p - S0(r)p = E0*p\n");
    printf("  Ground state: E0 = %.10f\n", E_ground);

    if (E_ground < 0) {
        printf("  UNSTABLE: Skyrmion generates pseudoscalar halo\n");
        printf("  Range: %.4f code = %.4f fm (short-range, NOT 1/r)\n",
               range_ground, range_ground * 0.5624);
        printf("  Two-soliton interaction: Yukawa ~ e^{-%.4f D}/D\n\n", kappa_ground);
        printf("CONCLUSION: The constraint coupling produces a localized\n");
        printf("pseudoscalar halo, analogous to the pion cloud. This is a\n");
        printf("nuclear-scale effect, NOT long-range gravity. For 1/r gravity,\n");
        printf("a linear source term (g_top B0) is still required.\n");
    } else {
        printf("  STABLE: p = 0 is energetically favored\n");
        printf("  No new physics from self-referential mapping.\n");
    }

    /* Cleanup */
    free(S0);
    free(T0);
    free(B0);
    free(u_tmp);
    free(p_eigen);
    free_profile(&prof);

    return 0;
}
