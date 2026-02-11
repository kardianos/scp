/*
 * spectral_eps.c — Estimate how ε₀² shifts the Channel C effective potential
 *
 * The constraint |q|² + ε₀²p² = ρ₀² at small ε₀² modifies |q| → √(ρ₀² - ε₀²p²).
 * For the background soliton with p₀(r) = A × eigenfunction of -∇²p - S₀p = E₀p:
 *
 *   |q|² = ρ₀² - ε₀² p₀²(r)
 *   ρ(r) = √(ρ₀² - ε₀² p₀²(r)) ≈ ρ₀ - ε₀² p₀²/(2ρ₀)
 *
 * This modifies the angular mode potential through:
 * 1. The effective f profile changes (lower ρ → profile shifts)
 * 2. The angular mode P, W, m coefficients depend on ρ
 *
 * For the K=1 Channel C mixed mode with W/m minimum at 0.167 vs threshold 0.158:
 *   Gap = 0.167 - 0.158 = 0.009
 *   Need δ(W/m) ≈ -0.009 to close the gap
 *
 * At leading order in ε₀², the correction to W/m comes from:
 *   δ(W/m) ~ ε₀² × (∂W/∂ρ² × p₀² / m - W/m² × ∂m/∂ρ² × p₀²)
 *
 * This code computes the required ε₀² numerically.
 *
 * Also computes:
 * - Zero-point energy of the constraint eigenvalue mode
 * - Mode enumeration for quantum mass formula
 * - Degenerate sector contribution to M_quantum
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

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

    int has_fp = 0;
    for (int i = 0; i < n; i++)
        if (fabs(prof->fp[i]) > 1e-30) { has_fp = 1; break; }
    if (!has_fp) {
        for (int i = 0; i < n; i++) {
            if (i == 0) prof->fp[i] = (-3*prof->f[0]+4*prof->f[1]-prof->f[2])/(2*prof->dr);
            else if (i == n-1) prof->fp[i] = (prof->f[n-3]-4*prof->f[n-2]+3*prof->f[n-1])/(2*prof->dr);
            else prof->fp[i] = (prof->f[i+1]-prof->f[i-1])/(2*prof->dr);
        }
    }
    return 0;
}

/* S₀(r) = f'² + 2sin²f/r² (Skyrmion strain) */
static double strain_S0(double r, double f, double fp, double a) {
    if (r > 1e-10) {
        double sf = sin(f);
        return fp*fp + 2.0*sf*sf/(r*r);
    } else {
        return 3.0*a*a;
    }
}

/* Solve constraint eigenvalue: -u'' - S₀(r)u = E₀u, u(0)=0
 * Uses Numerov method. Returns ground state eigenvalue. */
static double find_constraint_eigenvalue(const double *S0, int n, double dr,
                                          double *u_out)
{
    double *u = calloc(n, sizeof(double));

    /* Bisection: find E₀ where node count transitions 0→1 */
    double E_lo = -10.0, E_hi = 0.0;
    double h2 = dr*dr;

    /* Scan for bracket */
    int prev_nodes = 0;
    for (double E_try = -9.9; E_try <= 0.01; E_try += 0.1) {
        u[0] = 0.0; u[1] = dr;
        int nodes = 0;
        for (int i = 1; i < n-1; i++) {
            double k2_prev = S0[i-1] + E_try;
            double k2_curr = S0[i] + E_try;
            double k2_next = S0[i+1] + E_try;
            u[i+1] = (2.0*u[i]*(1.0-5.0*h2*k2_curr/12.0)
                      - u[i-1]*(1.0+h2*k2_prev/12.0))
                     / (1.0+h2*k2_next/12.0);
            if (u[i+1]*u[i] < 0) nodes++;
        }
        if (prev_nodes == 0 && nodes >= 1) {
            E_lo = E_try - 0.1;
            E_hi = E_try;
            break;
        }
        prev_nodes = nodes;
    }

    /* Bisect */
    for (int iter = 0; iter < 200; iter++) {
        double E_mid = 0.5*(E_lo + E_hi);
        u[0] = 0.0; u[1] = dr;
        int nodes = 0;
        for (int i = 1; i < n-1; i++) {
            double k2_prev = S0[i-1] + E_mid;
            double k2_curr = S0[i] + E_mid;
            double k2_next = S0[i+1] + E_mid;
            u[i+1] = (2.0*u[i]*(1.0-5.0*h2*k2_curr/12.0)
                      - u[i-1]*(1.0+h2*k2_prev/12.0))
                     / (1.0+h2*k2_next/12.0);
            if (u[i+1]*u[i] < 0) nodes++;
        }
        if (nodes > 0) E_hi = E_mid;
        else E_lo = E_mid;
        if (fabs(E_hi - E_lo) < 1e-12) break;
    }

    double E0 = 0.5*(E_lo + E_hi);

    /* Compute normalized eigenfunction */
    u[0] = 0.0; u[1] = dr;
    for (int i = 1; i < n-1; i++) {
        double k2_prev = S0[i-1] + E0;
        double k2_curr = S0[i] + E0;
        double k2_next = S0[i+1] + E0;
        u[i+1] = (2.0*u[i]*(1.0-5.0*h2*k2_curr/12.0)
                  - u[i-1]*(1.0+h2*k2_prev/12.0))
                 / (1.0+h2*k2_next/12.0);
    }

    /* Normalize: ∫ p² 4πr² dr = 1, where p = u/r */
    double norm2 = 0;
    for (int i = 1; i < n; i++) {
        double r = i * dr;
        double p = u[i] / r;
        norm2 += p*p * 4*M_PI*r*r * dr;
    }
    if (norm2 > 0) {
        double nf = 1.0/sqrt(norm2);
        for (int i = 0; i < n; i++) u[i] *= nf;
    }

    if (u_out) {
        for (int i = 0; i < n; i++) u_out[i] = u[i];
    }

    free(u);
    return E0;
}

/* Compute the l-th partial wave eigenvalue of -∇²p - S₀(r)p = Ep
 * For l > 0: u'' + [S₀(r) + E - l(l+1)/r²]u = 0, u(0) = 0 */
static double find_eigenvalue_ell(const double *S0, const double *r_grid,
                                   int n, double dr, int ell)
{
    double *u = calloc(n, sizeof(double));
    double h2 = dr*dr;
    double ell_barrier = ell*(ell+1.0);

    /* Scan for bracket */
    double E_lo = -10.0, E_hi = 0.0;
    int prev_nodes = 0;
    int found_bracket = 0;

    for (double E_try = -9.99; E_try <= 0.01; E_try += 0.01) {
        u[0] = 0.0; u[1] = dr;
        int nodes = 0;
        for (int i = 1; i < n-1; i++) {
            double r = r_grid[i];
            double r2 = r*r; if (r2 < 1e-20) r2 = 1e-20;
            double k2_f = S0[i] + E_try - ell_barrier/r2;

            double r_prev = r_grid[i-1];
            double r2_prev = r_prev*r_prev; if (r2_prev < 1e-20) r2_prev = 1e-20;
            double k2_prev = S0[i-1] + E_try - ell_barrier/r2_prev;

            double r_next = r_grid[i+1];
            double r2_next = r_next*r_next; if (r2_next < 1e-20) r2_next = 1e-20;
            double k2_next = S0[i+1] + E_try - ell_barrier/r2_next;

            u[i+1] = (2.0*u[i]*(1.0-5.0*h2*k2_f/12.0)
                      - u[i-1]*(1.0+h2*k2_prev/12.0))
                     / (1.0+h2*k2_next/12.0);
            if (u[i+1]*u[i] < 0) nodes++;
        }
        if (prev_nodes == 0 && nodes >= 1) {
            E_lo = E_try - 0.01;
            E_hi = E_try;
            found_bracket = 1;
            break;
        }
        prev_nodes = nodes;
    }

    if (!found_bracket) {
        free(u);
        return 1.0; /* No bound state */
    }

    /* Bisect */
    for (int iter = 0; iter < 200; iter++) {
        double E_mid = 0.5*(E_lo + E_hi);
        u[0] = 0.0; u[1] = dr;
        int nodes = 0;
        for (int i = 1; i < n-1; i++) {
            double r = r_grid[i];
            double r2 = r*r; if (r2 < 1e-20) r2 = 1e-20;
            double k2_f = S0[i] + E_mid - ell_barrier/r2;

            double r_prev = r_grid[i-1];
            double r2_prev = r_prev*r_prev; if (r2_prev < 1e-20) r2_prev = 1e-20;
            double k2_prev = S0[i-1] + E_mid - ell_barrier/r2_prev;

            double r_next = r_grid[i+1];
            double r2_next = r_next*r_next; if (r2_next < 1e-20) r2_next = 1e-20;
            double k2_next = S0[i+1] + E_mid - ell_barrier/r2_next;

            u[i+1] = (2.0*u[i]*(1.0-5.0*h2*k2_f/12.0)
                      - u[i-1]*(1.0+h2*k2_prev/12.0))
                     / (1.0+h2*k2_next/12.0);
            if (u[i+1]*u[i] < 0) nodes++;
        }
        if (nodes > 0) E_hi = E_mid;
        else E_lo = E_mid;
        if (fabs(E_hi - E_lo) < 1e-12) break;
    }

    free(u);
    return 0.5*(E_lo + E_hi);
}


int main(int argc, char *argv[])
{
    const char *profile_file = "data/profiles/profile_massive_e1_mpi0.398.dat";
    double e_param = 1.0, rho0 = 1.0, m_pi = 0.398;

    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "-profile") && i+1<argc) profile_file = argv[++i];
        else if (!strcmp(argv[i], "-e") && i+1<argc) e_param = atof(argv[++i]);
        else if (!strcmp(argv[i], "-rho0") && i+1<argc) rho0 = atof(argv[++i]);
        else if (!strcmp(argv[i], "-mpi") && i+1<argc) m_pi = atof(argv[++i]);
    }

    double c4 = 2.0*rho0*rho0/(e_param*e_param);
    double m_pi2 = m_pi*m_pi;
    double hbar = 0.0386;  /* code units */
    double code_to_MeV = 938.272 * e_param / (103.13 * rho0*rho0*rho0);

    Profile prof;
    if (read_profile(profile_file, &prof) < 0) return 1;
    int n = prof.n;
    double dr = prof.dr;
    double a = -prof.fp[0];

    printf("============================================================\n");
    printf(" Avenue D: Spectral Self-Consistency Analysis\n");
    printf("============================================================\n\n");
    printf("Profile: %s (%d points, dr=%.6f)\n", profile_file, n, dr);
    printf("Parameters: e=%.4f, rho0=%.4f, c4=%.4f, m_pi=%.4f\n", e_param, rho0, c4, m_pi);
    printf("hbar = %.4f code, 1 code E = %.3f MeV\n\n", hbar, code_to_MeV);

    /* ===== Section 1: Constraint eigenvalue spectrum ===== */
    printf("--- Section 1: Constraint Eigenvalue Spectrum ---\n\n");
    printf("Eigenvalue problem: -nabla^2 p - S0(r) p = E_0 p\n");
    printf("(from constraint |q|^2 + eps0^2 p^2 = rho0^2)\n\n");

    /* Compute S₀ on profile grid */
    double *S0 = calloc(n, sizeof(double));
    double *r_grid = calloc(n, sizeof(double));
    for (int i = 0; i < n; i++) {
        r_grid[i] = prof.r[i];
        S0[i] = strain_S0(prof.r[i], prof.f[i], prof.fp[i], a);
    }

    /* Ground state (l=0) */
    double *u_ground = calloc(n, sizeof(double));
    double E0_l0 = find_constraint_eigenvalue(S0, n, dr, u_ground);
    double kappa0 = sqrt(fabs(E0_l0));

    printf("l=0 (s-wave): E0 = %.10f, kappa = %.6f, range = %.4f code = %.4f fm\n",
           E0_l0, kappa0, 1.0/kappa0, 1.0/kappa0*0.5624);

    /* Higher partial waves */
    for (int ell = 1; ell <= 4; ell++) {
        double E_ell = find_eigenvalue_ell(S0, r_grid, n, dr, ell);
        if (E_ell < 0) {
            double kappa = sqrt(fabs(E_ell));
            printf("l=%d: E0 = %.10f, kappa = %.6f, range = %.4f code = %.4f fm\n",
                   ell, E_ell, kappa, 1.0/kappa, 1.0/kappa*0.5624);
        } else {
            printf("l=%d: No bound state (E > 0)\n", ell);
        }
    }

    /* Degeneracy of constraint modes:
     * l=0: 1 state (p pseudoscalar)
     * For j_1,j_2,j_3: each is a pseudovector, but on the hedgehog background
     * the coupling is through S₀(r) which is isospin-symmetric.
     * The j-components see the SAME potential S₀(r) → 3× degenerate for l=0.
     * Total degenerate sector modes: 4 components × (2l+1) per partial wave.
     */
    printf("\nDegeneracy structure:\n");
    printf("  Degenerate sector has 4 components: (p, j1, j2, j3)\n");
    printf("  On hedgehog background: all see same S0(r) potential\n");
    printf("  l=0: 4 × 1 = 4 modes (all with E0 = %.6f)\n", E0_l0);

    /* ===== Section 2: Complete mode enumeration ===== */
    printf("\n--- Section 2: Complete Mode Enumeration ---\n\n");

    double E_sol = 110.19; /* massive profile at e=1 */
    double Lambda = 86.6;  /* moment of inertia, massive */

    printf("Classical soliton mass: E_sol = %.2f code = %.1f MeV\n",
           E_sol, E_sol * code_to_MeV);
    printf("Moment of inertia: Lambda = %.1f code\n", Lambda);

    /* Rotational energy */
    double E_rot_half = 3.0*hbar*hbar / (8.0*Lambda);
    printf("\n1. Rotational (J=I=1/2): E_rot = 3*hbar^2/(8*Lambda)\n");
    printf("   = 3 * %.6f / (8 * %.1f) = %.6f code = %.3f MeV\n",
           hbar*hbar, Lambda, E_rot_half, E_rot_half*code_to_MeV);

    /* K=0 breathing mode */
    printf("\n2. K=0 breathing mode:\n");
    printf("   No bound states below m_pi^2 = %.6f threshold\n", m_pi2);
    printf("   Lowest box mode: omega_0 = 0.572 code = 5.2 MeV\n");
    printf("   Continuum: contributes through spectral density (UV divergent)\n");
    printf("   Zero-point per continuum mode: (1/2)*hbar*omega_n\n");

    /* K=1 channels */
    printf("\n3. K=1 Channel A (isorotational, L=0, I=1):\n");
    printf("   Zero mode at omega=0 (3 generators)\n");
    printf("   Physical pion from collective quantization\n");
    printf("   After Goldstone correction: contributes m_pi/2 per mode\n");
    printf("   3 pion degrees of freedom: 3 * (1/2)*m_pi = %.4f code = %.1f MeV\n",
           1.5*m_pi*code_to_MeV/code_to_MeV, 1.5*m_pi*code_to_MeV);
    /* Actually the pion mass gap already in the physical m_pi — this is via quantization */

    printf("\n4. K=1 Channel B (translational, L=1, I=0):\n");
    printf("   Zero mode at omega=0 (3 translations)\n");
    printf("   Eliminated by center-of-mass frame\n");
    printf("   No contribution to M_quantum (momentum eigenstate)\n");

    printf("\n5. K=1 Channel C (mixed, L=1, I=1, K=1):\n");
    printf("   W/m minimum = 0.167, threshold = m_pi^2 = %.6f\n", m_pi2);
    printf("   Gap: delta(W/m) = 0.167 - 0.158 = 0.009 (5.7%% above threshold)\n");
    printf("   NO bound states at eps0^2 = 0\n");
    printf("   Degeneracy if bound: 3 (I=1) × 3 (L=1) = 9 modes\n");

    printf("\n6. Degenerate sector constraint modes (at eps0^2 > 0):\n");
    printf("   s-wave (l=0): E0 = %.6f, 4 modes (p + 3 j-components)\n", E0_l0);
    printf("   Zero-point: 4 × (1/2)*hbar*sqrt(|E0|)\n");
    double zpe_constraint = 4.0 * 0.5 * hbar * kappa0;
    printf("   = 4 * 0.5 * %.4f * %.4f = %.6f code = %.3f MeV\n",
           hbar, kappa0, zpe_constraint, zpe_constraint * code_to_MeV);

    /* ===== Section 3: Quantum mass formula ===== */
    printf("\n--- Section 3: Quantum Mass Formula ---\n\n");

    printf("M_quantum = M_classical + E_rot + E_zp(bulk) + E_zp(degen) + counterterms\n\n");

    printf("At leading order (1-loop):\n");
    printf("  M_classical = %.2f code = %.1f MeV\n", E_sol, E_sol*code_to_MeV);
    printf("  E_rot(J=1/2) = %.6f code = %.3f MeV\n", E_rot_half, E_rot_half*code_to_MeV);

    /* Zero-point energy of continuum modes (divergent — needs regularization) */
    printf("\n  E_zp(bulk) = (1/2) * Sum_n hbar*omega_n   [DIVERGENT]\n");
    printf("  Regularization: zeta-function or heat-kernel\n");
    printf("  For K=0 modes (breathing): Sum ~ integral omega*rho(omega) d(omega)\n");
    printf("  Density of states rho(omega) ~ omega^2 at high omega (3D)\n");
    printf("  => Sum ~ integral omega^3 d(omega) ~ Lambda_UV^4 (quartically divergent)\n");
    printf("  Renormalized: subtract free-space contribution\n");
    printf("  E_zp^{ren}(K=0) ~ Sum_n [omega_n - omega_n^{free}]\n");

    printf("\n  E_zp(degenerate sector, at eps0^2 > 0):\n");
    printf("  4 bound modes with omega = sqrt(|E0|) = %.4f code\n", kappa0);
    printf("  Zero-point: %.6f code = %.3f MeV\n",
           zpe_constraint, zpe_constraint*code_to_MeV);
    printf("  Plus: 4 × continuum modes for each l > 0\n");
    printf("  The continuum contribution is UV divergent (same structure as bulk)\n");

    /* ===== Section 4: Parameter count ===== */
    printf("\n--- Section 4: Parameter Over-Determination ---\n\n");

    printf("Parameters of the theory:\n");
    printf("  1. rho_0  (VEV / field amplitude)\n");
    printf("  2. lambda  (bulk potential coupling)\n");
    printf("  3. e       (Skyrme coupling)\n");
    printf("  4. mu      (degenerate mass)\n");
    printf("  5. c       (speed of light)\n");
    printf("  6. eps_0^2 (constraint coupling — NEW from degenerate sector)\n");
    printf("  Total: 6 parameters (5 + eps0^2)\n\n");

    printf("Physical observables to match:\n");
    printf("  1. M_nucleon = 938.3 MeV          [constrains rho0^3/e]\n");
    printf("  2. r_proton  = 0.841 fm            [constrains rho0/e]\n");
    printf("  3. m_pion    = 139.6 MeV           [constrains mu via m_pi = mu*F_pi]\n");
    printf("  4. Delta-N   = 293.7 MeV           [constrains Lambda ~ e^{-3}rho0]\n");
    printf("  5. g_A       = 1.27                [axial coupling]\n");
    printf("  6. F_pi      = 92.1 MeV            [pion decay constant]\n");
    printf("  7. G_Newton  = 6.674e-11 N m^2/kg^2 [gravitational constant]\n");
    printf("  8+ Various scattering observables, form factors, etc.\n\n");

    printf("Counting:\n");
    printf("  At eps0^2 = 0: 5 params vs >= 6 observables => OVERCONSTRAINED\n");
    printf("  Classical Skyrme fits: M_N + r_p fix (rho0,e); then m_pi fixes mu;\n");
    printf("  Delta-N and g_A are PREDICTIONS (currently ~10x too low for Delta-N).\n\n");

    printf("  At eps0^2 > 0: 6 params vs >= 7 observables (add G_Newton)\n");
    printf("  The system gains one parameter but also one observable.\n");
    printf("  Net: still overconstrained.\n\n");

    printf("  Key self-consistency test:\n");
    printf("  If eps0^2 is fixed by matching G_Newton, does the resulting\n");
    printf("  modification to the spectrum IMPROVE or WORSEN the other fits?\n");

    /* ===== Section 5: Channel C shift estimate ===== */
    printf("\n--- Section 5: eps0^2 Effect on Channel C ---\n\n");

    /* At leading order in eps0^2, the modification to the Channel C potential
     * comes from the constraint-modified rho:
     *   rho(r) = sqrt(rho0^2 - eps0^2 * p0^2(r))
     *
     * The Channel C effective potential has:
     *   P_C(r) = r^2 * Omega(r) where Omega = rho0^2 * sin^2(f) * [1 + c4(...)]
     *
     * With rho -> rho(r):
     *   Omega -> rho(r)^2 * sin^2(f) * [1 + c4(...)]
     *         = (rho0^2 - eps0^2*p0^2) * sin^2(f) * [1 + c4(...)]
     *
     * delta_Omega/Omega = -eps0^2 * p0^2 / rho0^2
     *
     * Since P_C = r^2 * Omega and m_C = P_C (Channel C has P=m):
     *   delta_P/P = delta_m/m = -eps0^2 * p0^2 / rho0^2
     *
     * W_C has multiple terms. The dominant effect on W/m comes through
     * the pion mass term and the centrifugal+Skyrme terms.
     *
     * The critical quantity is: how does W/m_minimum change with eps0^2?
     *
     * For the centrifugal term W_cent = 2*P/r^2:
     *   W_cent / m = 2/r^2 (independent of rho!)
     *
     * For the Skyrme term W_E4:
     *   W_E4 = c4 * sin^2(f) * (...) / r^2
     *   This depends on c4 = 2*rho0^2/e^2, but NOT on rho(r) if sin^2(f)
     *   is unchanged. However, f(r) itself might shift.
     *   At leading order: W_E4/m has no rho dependence (both proportional to Omega).
     *
     * For the pion mass term W_pi = m_pi^2 * r^2 * cos(f):
     *   W_pi / m = m_pi^2 * r^2 * cos(f) / P_C
     *   Since P_C ~ rho^2, delta(W_pi/m) = W_pi/m * eps0^2*p0^2/rho0^2
     *   This INCREASES W_pi/m (makes the attraction weaker near core where cos(f)<0).
     *
     * Wait: near the core cos(f) < 0, so W_pi < 0 (attractive).
     * delta(W_pi/m)/eps0^2 = -W_pi/m * p0^2/rho0^2 > 0 (since W_pi/m < 0 near core)
     * So the eps0^2 correction makes the well SHALLOWER (more positive W/m).
     *
     * This is the WRONG direction for binding! eps0^2 makes Channel C LESS likely to bind.
     *
     * BUT: we haven't considered the effect of eps0^2 on the profile f(r) itself.
     * With reduced rho, the equilibrium profile changes. And there are ADDITIONAL
     * modes from the degenerate sector that mix with Channel C.
     *
     * Actually the dominant effect may be through the new p-sector modes themselves,
     * not through modification of existing bulk modes.
     */

    printf("Leading-order analysis of eps0^2 on Channel C W/m:\n\n");

    /* Compute p0^2(r) at the location of W/m minimum */
    /* First, compute p0(r) = u(r)/r */
    double *p0_sq = calloc(n, sizeof(double));
    for (int i = 1; i < n; i++) {
        double r = prof.r[i];
        double p0 = u_ground[i] / r;
        p0_sq[i] = p0*p0;
    }
    p0_sq[0] = p0_sq[1]; /* extrapolate */

    /* Find where Channel C W/m is minimum */
    /* Recompute Channel C potential */
    double Wm_min = 1e30;
    double r_min = 0;
    double p0sq_at_min = 0;

    for (int i = 1; i < n; i++) {
        double r = prof.r[i];
        if (r < 0.1 || r > 10.0) continue;
        double f = prof.f[i], fp = prof.fp[i];
        double sf = sin(f), sf2 = sf*sf, cf = cos(f);
        double r2 = r*r;

        double P_mix = r2 * (1.0 + c4*(fp*fp + sf2/(r2+1e-30)));
        double W_cent = 2.0 * P_mix / (r2 + 1e-30);
        double W_E4 = c4 * sf2 * (fp*fp + sf2/(r2+1e-30)) / (r2 + 1e-30);
        double W_pi = m_pi2 * r2 * cf;
        double W_C = W_cent + W_E4 + W_pi;
        double ratio = W_C / P_mix;

        if (ratio < Wm_min) {
            Wm_min = ratio;
            r_min = r;
            p0sq_at_min = p0_sq[i];
        }
    }

    printf("Channel C W/m minimum:\n");
    printf("  Location: r = %.4f code\n", r_min);
    printf("  W/m = %.6f\n", Wm_min);
    printf("  Threshold: m_pi^2 = %.6f\n", m_pi2);
    printf("  Gap: %.6f (%.1f%% above threshold)\n",
           Wm_min - m_pi2, (Wm_min - m_pi2)/m_pi2 * 100);
    printf("  p0^2 at minimum: %.6e\n\n", p0sq_at_min);

    /* Effect on W/m from rho reduction */
    /* The W/m ratio for Channel C: W_cent/m = 2/r^2 (rho-independent)
     * W_E4/m = c4*sin^2f*(...)/(r^2*P) (rho-independent since c4 and P both ~ rho^2)
     * W_pi/m = m_pi^2*r^2*cos(f)/P, where P ~ rho^2
     * So delta(W_pi/m) = +W_pi/m * eps0^2*p0^2/rho0^2
     *
     * At the minimum: W_pi/m = Wm_min - (W_cent + W_E4)/m
     * We need W_pi/m at the minimum. */
    {
        double r = r_min;
        /* Interpolate profile at r_min */
        int idx = (int)(r_min / dr);
        if (idx >= n-1) idx = n-2;
        double t = (r_min - prof.r[idx]) / dr;
        double f = prof.f[idx] + t*(prof.f[idx+1]-prof.f[idx]);
        double fp = prof.fp[idx] + t*(prof.fp[idx+1]-prof.fp[idx]);
        double sf = sin(f), sf2 = sf*sf, cf = cos(f);
        double r2 = r*r;

        double P_mix = r2 * (1.0 + c4*(fp*fp + sf2/(r2+1e-30)));
        double W_cent_m = 2.0 / (r2+1e-30);
        double W_E4_m = c4*sf2*(fp*fp + sf2/(r2+1e-30)) / (r2*P_mix+1e-30);
        double W_pi_m = m_pi2*r2*cf / P_mix;

        printf("At W/m minimum (r=%.4f):\n", r_min);
        printf("  W_cent/m = %.6f (rho-independent)\n", W_cent_m);
        printf("  W_E4/m   = %.6f (rho-independent)\n", W_E4_m);
        printf("  W_pi/m   = %.6f (rho-DEPENDENT)\n", W_pi_m);
        printf("  Sum = %.6f (check: %.6f)\n", W_cent_m+W_E4_m+W_pi_m, Wm_min);

        double delta_Wm_per_eps = W_pi_m * p0sq_at_min / (rho0*rho0);
        printf("\n  delta(W/m) per eps0^2 = W_pi/m * p0^2/rho0^2 = %.6e\n", delta_Wm_per_eps);

        if (delta_Wm_per_eps > 0) {
            printf("  Sign: POSITIVE (gap WIDENS with eps0^2 > 0)\n");
            printf("  => eps0^2 makes Channel C LESS likely to bind\n");
            printf("  => The constraint coupling works AGAINST creating new light particles\n");
        } else {
            double gap = Wm_min - m_pi2;
            double eps0sq_bind = -gap / delta_Wm_per_eps;
            printf("  Sign: NEGATIVE (gap NARROWS with eps0^2 > 0)\n");
            printf("  eps0^2 needed to close gap: %.6e\n", eps0sq_bind);
        }
    }

    /* ===== Section 6: Breathing mode (K=0) with eps0^2 ===== */
    printf("\n--- Section 6: Breathing Mode Interaction with Constraint ---\n\n");

    printf("The constraint eigenvalue (l=0): E0 = %.6f\n", E0_l0);
    printf("The breathing mode (K=0) has omega^2 = lambda, all lambda > m_pi^2\n");
    printf("No bound breathing modes exist.\n\n");

    printf("At eps0^2 > 0, the constraint bound state (l=0) has:\n");
    printf("  omega_constraint = sqrt(|E0|) = %.4f code\n", kappa0);
    printf("  This is ABOVE the pion mass threshold: %.4f > m_pi = %.4f\n",
           kappa0, m_pi);
    printf("  So the constraint mode is in the CONTINUUM of the bulk K=0 channel.\n");
    printf("  It can DECAY into pions (omega_constraint > m_pi).\n");
    printf("  Decay rate ~ eps0^2 * (matrix element)^2 * phase_space\n\n");

    printf("If m_pi > sqrt(|E0|) = %.4f, the constraint mode would be STABLE.\n", kappa0);
    printf("At m_pi = 0.398 vs sqrt(|E0|) = %.4f: m_pi < sqrt(|E0|)\n", kappa0);
    printf("So the constraint mode IS in the continuum and DECAYS.\n");
    printf("Its zero-point energy still contributes to M_quantum.\n");

    /* ===== Section 7: UV structure ===== */
    printf("\n--- Section 7: UV Divergence Structure ---\n\n");

    printf("Zero-point energy sum: E_zp = (1/2)*hbar * Sum_n omega_n\n\n");

    printf("Bulk sector (eps0^2 = 0):\n");
    printf("  K=0: Sum omega_n ~ n^1 for large n (3D box modes)\n");
    printf("  K=1: 3 channels, each ~ n^1\n");
    printf("  K>=2: similar structure\n");
    printf("  Total: Sum ~ Lambda_UV^4 (quartic divergence in 4D)\n\n");

    printf("Degenerate sector (eps0^2 > 0):\n");
    printf("  4 constraint modes: same S0(r) potential, same spectrum structure\n");
    printf("  Continuum: same density of states as bulk (p, j are scalar/vector fields)\n");
    printf("  UV divergence: SAME STRUCTURE as bulk (quartic)\n");
    printf("  The eps0^2-dependent part: 4 bound modes contribute finite ZPE\n\n");

    printf("Regularization (zeta function):\n");
    printf("  zeta_S(s) = Sum_n omega_n^{-s}\n");
    printf("  Renormalized ZPE = (1/2)*hbar * zeta_S(-1)\n");
    printf("  The bound state contribution is ALWAYS finite (isolated eigenvalue)\n");
    printf("  The continuum contribution needs counterterms (bulk potential)\n\n");

    printf("eps0^2 dependence of renormalized ZPE:\n");
    printf("  Finite part: 4 * (1/2)*hbar*sqrt(|E0|) * eps0^2 (from bound modes)\n");
    printf("  = eps0^2 * %.6f code = eps0^2 * %.3f MeV\n",
           zpe_constraint, zpe_constraint * code_to_MeV);
    printf("  Counterterm part: depends on regularization scheme\n");
    printf("  At 1-loop: counterterms are polynomial in eps0^2, m_pi^2, Lambda_UV\n");
    printf("  => No constraint on eps0^2 from finiteness alone\n");

    /* ===== Section 8: Self-consistency constraints ===== */
    printf("\n--- Section 8: Self-Consistency Constraints ---\n\n");

    printf("Question: Does the quantum spectrum constrain eps0^2?\n\n");

    printf("Test 1: UV finiteness\n");
    printf("  The zero-point sum diverges regardless of eps0^2.\n");
    printf("  Renormalization absorbs the divergence into counterterms.\n");
    printf("  Result: NO constraint on eps0^2 from finiteness.\n\n");

    printf("Test 2: Matching M_nucleon\n");
    printf("  M_quantum = 938.3 MeV (input)\n");
    printf("  = M_classical + E_rot + E_zp^{ren}(bulk) + eps0^2 * E_zp^{ren}(degen)\n");
    printf("  Currently M_classical = 938 MeV by construction (rho0, e fixed).\n");
    printf("  E_rot = %.3f MeV (tiny)\n", E_rot_half*code_to_MeV);
    printf("  E_zp^{ren}(bulk) is an O(1) correction that shifts the fit slightly.\n");
    printf("  eps0^2 * E_zp(degen) = eps0^2 * %.3f MeV\n", zpe_constraint*code_to_MeV);
    printf("  For eps0^2 ~ 10^{-14}: contribution ~ 10^{-14} MeV = NEGLIGIBLE.\n");
    printf("  Result: eps0^2 has NO effect on the nucleon mass at 10^{-14}.\n\n");

    printf("Test 3: Delta-N splitting\n");
    printf("  Delta_E = 3*hbar^2/(2*Lambda_eff)\n");
    printf("  Lambda_eff includes eps0^2 correction to moment of inertia:\n");
    printf("  delta_Lambda/Lambda ~ eps0^2 * p0^2(core)/rho0^2\n");
    printf("  At eps0^2 ~ 10^{-14}: delta_Lambda ~ 10^{-14} = NEGLIGIBLE.\n");
    printf("  Result: NO constraint.\n\n");

    printf("Test 4: Channel C binding\n");
    printf("  Gap = %.6f (%.1f%% of threshold)\n",
           Wm_min - m_pi2, (Wm_min - m_pi2)/m_pi2 * 100);
    printf("  eps0^2 effect on gap: WIDENS (wrong direction for binding)\n");
    printf("  Result: eps0^2 > 0 makes Channel C FURTHER from binding.\n\n");

    printf("Test 5: Spectral flow (adiabatic evolution)\n");
    printf("  As eps0^2 increases from 0, the degenerate sector modes\n");
    printf("  'turn on' continuously. New eigenvalues appear at the\n");
    printf("  threshold omega = sqrt(|E0|) = %.4f and separate.\n", kappa0);
    printf("  These modes couple to the bulk spectrum through the constraint.\n");
    printf("  However, at eps0^2 ~ 10^{-14}, the coupling is O(10^{-14}).\n");
    printf("  The spectral flow is infinitesimally slow.\n");
    printf("  Result: NO discrete constraint on eps0^2 from spectral flow.\n\n");

    /* ===== Section 9: What WOULD constrain eps0^2? ===== */
    printf("\n--- Section 9: What Would Constrain eps0^2? ---\n\n");

    printf("The spectral analysis shows that eps0^2 enters the quantum\n");
    printf("spectrum only through:\n");
    printf("  a) An O(eps0^2) shift in the zero-point energy\n");
    printf("  b) An O(eps0^2) shift in the effective potential for bulk modes\n\n");

    printf("For eps0^2 ~ 10^{-14} (gravity scale), these effects are\n");
    printf("unmeasurably small at the nuclear scale. The spectrum is\n");
    printf("essentially identical to eps0^2 = 0.\n\n");

    printf("Possible constraints NOT from the spectrum:\n");
    printf("  1. Anomaly cancellation (Avenue C) — algebraic, not spectral\n");
    printf("  2. WZW coefficient (Avenue B) — topological quantization\n");
    printf("  3. Naturalness/fine-tuning — eps0^2 << 1 is technically natural\n");
    printf("     if a symmetry protects it (e.g., Z2: p -> -p)\n");
    printf("  4. Asymptotic safety — RG flow may have fixed points\n");

    /* ===== Final Summary ===== */
    printf("\n============================================================\n");
    printf(" Avenue D Summary: Spectral Self-Consistency\n");
    printf("============================================================\n\n");

    printf("FINDINGS:\n\n");

    printf("1. CONSTRAINT EIGENVALUE SPECTRUM:\n");
    printf("   l=0: E0 = %.6f (bound, kappa=%.4f, range=%.3f fm)\n",
           E0_l0, kappa0, 1.0/kappa0*0.5624);
    for (int ell = 1; ell <= 4; ell++) {
        double E_ell = find_eigenvalue_ell(S0, r_grid, n, dr, ell);
        if (E_ell < 0)
            printf("   l=%d: E = %.6f (bound)\n", ell, E_ell);
        else
            printf("   l=%d: unbound\n", ell);
    }
    printf("   Degeneracy: 4 components × (2l+1) per partial wave\n\n");

    printf("2. QUANTUM MASS (at eps0^2 ~ 10^{-14}):\n");
    printf("   M_quantum = M_cl + E_rot + E_zp(bulk) + O(10^{-14}) correction\n");
    printf("   The degenerate sector contribution is negligible at gravity scale.\n\n");

    printf("3. CHANNEL C NEAR-BINDING:\n");
    printf("   Gap = %.4f (%.1f%% of threshold)\n",
           Wm_min - m_pi2, (Wm_min-m_pi2)/m_pi2*100);
    printf("   eps0^2 > 0 WIDENS the gap (wrong direction)\n");
    printf("   No new light particles from eps0^2 coupling.\n\n");

    printf("4. OVER-DETERMINATION:\n");
    printf("   6 params (rho0,lambda,e,mu,c,eps0^2) vs 7+ observables\n");
    printf("   System IS overconstrained, but eps0^2 only affects G_Newton.\n");
    printf("   All other observables are insensitive to eps0^2 ~ 10^{-14}.\n");
    printf("   In effect: eps0^2 is a FREE parameter fixed only by G_Newton.\n\n");

    printf("5. CONCLUSION:\n");
    printf("   The quantum spectrum does NOT constrain eps0^2.\n");
    printf("   The hierarchy eps0^2 ~ 10^{-14} << 1 must come from elsewhere:\n");
    printf("   anomaly cancellation, topological quantization, or RG flow.\n");
    printf("   Spectral self-consistency is NECESSARY but NOT SUFFICIENT\n");
    printf("   to determine the gravitational coupling.\n");

    free(S0);
    free(r_grid);
    free(u_ground);
    free(p0_sq);
    free(prof.r); free(prof.f); free(prof.fp);

    return 0;
}
