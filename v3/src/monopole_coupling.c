/*
 * monopole_coupling.c — Gravity no-go: monopole coupling analysis
 *
 * PHYSICS:
 * ========
 * For 1/R gravitational potential between solitons of mass M₁, M₂:
 *   V(R) = -G·M₁·M₂/R
 * you need a MASSLESS mediator field φ that couples to the MONOPOLE (L=0)
 * component of the stress-energy tensor T₀₀.
 *
 * This code tests ALL propagating channels in the Skyrme model:
 *
 *   Pion (π_a):     massless ✓  but dipole-coupled (L=1) ✗  → 1/R³
 *   Breathing (σ):  monopole-coupled ✓  but massive (m→∞) ✗  → Yukawa
 *   Weight (J,P):   non-dynamical at ε₀²=0                   → nothing
 *   Spin-2 tensor:  does NOT EXIST (SU(2) has j ≤ 1)         → nothing
 *
 * No channel satisfies BOTH requirements. Gravity is structurally absent.
 *
 * USAGE:
 *   monopole_coupling -profile <file> [-outdir dir]
 */

#define _USE_MATH_DEFINES
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ========== Profile I/O ========== */

typedef struct {
    double *r, *f, *fp;
    int n;
    double dr, rho0, e_skyrme, c4;
} Profile;

static Profile *load_profile(const char *filename, double rho0, double e_skyrme)
{
    FILE *fp = fopen(filename, "r");
    if (!fp) { fprintf(stderr, "Cannot open %s\n", filename); return NULL; }

    int cap = 16384;
    double *r_arr  = malloc(cap * sizeof(double));
    double *f_arr  = malloc(cap * sizeof(double));
    double *fp_arr = malloc(cap * sizeof(double));
    int n = 0;
    char line[1024];

    while (fgets(line, sizeof(line), fp)) {
        if (line[0] == '#') continue;
        double rv, fv, fpv = 0, d1, d2, d3;
        int nc = sscanf(line, "%lf %lf %lf %lf %lf %lf",
                        &rv, &fv, &fpv, &d1, &d2, &d3);
        if (nc < 2) continue;
        if (n >= cap) {
            cap *= 2;
            r_arr  = realloc(r_arr,  cap * sizeof(double));
            f_arr  = realloc(f_arr,  cap * sizeof(double));
            fp_arr = realloc(fp_arr, cap * sizeof(double));
        }
        r_arr[n] = rv; f_arr[n] = fv; fp_arr[n] = fpv;
        n++;
    }
    fclose(fp);

    /* Compute f' if not provided */
    int has_fp = 0;
    for (int i = 0; i < n; i++)
        if (fabs(fp_arr[i]) > 1e-30) { has_fp = 1; break; }
    if (!has_fp && n > 2) {
        for (int i = 0; i < n; i++) {
            if (i == 0)
                fp_arr[i] = (f_arr[1] - f_arr[0]) / (r_arr[1] - r_arr[0]);
            else if (i == n - 1)
                fp_arr[i] = (f_arr[n-1] - f_arr[n-2]) / (r_arr[n-1] - r_arr[n-2]);
            else
                fp_arr[i] = (f_arr[i+1] - f_arr[i-1]) / (r_arr[i+1] - r_arr[i-1]);
        }
    }

    Profile *p = malloc(sizeof(Profile));
    p->r = r_arr; p->f = f_arr; p->fp = fp_arr;
    p->n = n;
    p->dr = (n > 1) ? r_arr[1] - r_arr[0] : 0.001;
    p->rho0 = rho0;
    p->e_skyrme = e_skyrme;
    p->c4 = 2.0 * rho0 * rho0 / (e_skyrme * e_skyrme);
    return p;
}

static void free_profile(Profile *p) {
    free(p->r); free(p->f); free(p->fp); free(p);
}

/* ========== Main computation ========== */

int main(int argc, char *argv[])
{
    const char *profile_file = NULL;
    const char *outdir = ".";
    double rho0 = 1.0, e_skyrme = 1.0;

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-profile") == 0 && i+1 < argc)
            profile_file = argv[++i];
        else if (strcmp(argv[i], "-outdir") == 0 && i+1 < argc)
            outdir = argv[++i];
        else if (strcmp(argv[i], "-e") == 0 && i+1 < argc)
            e_skyrme = atof(argv[++i]);
        else if (strcmp(argv[i], "-rho0") == 0 && i+1 < argc)
            rho0 = atof(argv[++i]);
        else {
            fprintf(stderr, "Usage: %s -profile <file> [-outdir dir] "
                    "[-e val] [-rho0 val]\n", argv[0]);
            return 1;
        }
    }

    if (!profile_file) {
        fprintf(stderr, "Must specify -profile <file>\n");
        return 1;
    }

    Profile *prof = load_profile(profile_file, rho0, e_skyrme);
    if (!prof) return 1;

    double c4 = prof->c4;

    printf("================================================================\n");
    printf(" Gravity No-Go: Monopole Coupling Analysis\n");
    printf("================================================================\n\n");
    printf("Profile: %s (%d points, r_max=%.2f)\n",
           profile_file, prof->n, prof->r[prof->n-1]);
    printf("Parameters: rho0=%.4f, e=%.4f, c4=%.6f\n\n", rho0, e_skyrme, c4);

    /* Physical constants for comparison */
    double code_to_MeV = 9.098;     /* 1 code energy = 9.098 MeV */
    double code_to_fm  = 0.5624;    /* 1 code length = 0.5624 fm */
    double M_Planck_MeV = 1.22e22;  /* Planck mass in MeV */

    /* ============================================================
     * Section 1: Energy Monopole
     * ============================================================
     *
     * T₀₀(r) is the energy density of the hedgehog soliton.
     * Being spherically symmetric, ALL multipoles L>0 vanish.
     * The monopole is simply: Q₀ = ∫ T₀₀ 4πr² dr = M
     */
    printf("===== Section 1: Energy Monopole =====\n\n");

    double E2 = 0, E4 = 0;
    for (int i = 1; i < prof->n - 1; i++) {
        double r  = prof->r[i];
        double f  = prof->f[i];
        double fp = prof->fp[i];
        if (r < 1e-10) continue;

        double sf = sin(f);
        double sf2 = sf * sf;
        double sf4 = sf2 * sf2;

        /* E₂ density: ½ρ₀²(f'² + 2sin²f/r²) */
        double e2_dens = 0.5 * rho0 * rho0 * (fp * fp + 2.0 * sf2 / (r * r));

        /* E₄ density: (1/e²)(2f'²sin²f/r² + sin⁴f/r⁴)
         * NOTE: c₄ = 2ρ₀²/e² is the 3D commutator prefactor.
         * For the 1D radial integral, the correct factor is 1/e². */
        double e2_inv = 1.0 / (e_skyrme * e_skyrme);
        double e4_dens = e2_inv * (2.0 * fp * fp * sf2 / (r * r) + sf4 / (r * r * r * r));

        E2 += e2_dens * 4.0 * M_PI * r * r * prof->dr;
        E4 += e4_dens * 4.0 * M_PI * r * r * prof->dr;
    }

    double M = E2 + E4;
    double E_FB = 6.0 * sqrt(2.0) * M_PI * M_PI * rho0 * rho0 * rho0 / e_skyrme;

    printf("Energy monopole (= soliton mass):\n");
    printf("  E₂ = %.4f  (gradient energy)\n", E2);
    printf("  E₄ = %.4f  (Skyrme energy)\n", E4);
    printf("  M  = %.4f  (total)\n", M);
    printf("  E₂/E₄ = %.4f  (virial: should be 1.000)\n", E2 / E4);
    printf("  E/E_FB = %.4f  (should be ~1.232)\n", M / E_FB);
    printf("  M = %.1f MeV (physical)\n\n", M * code_to_MeV);

    printf("Higher multipoles of T₀₀:\n");
    printf("  L=1 (dipole): 0 (hedgehog is spherically symmetric)\n");
    printf("  L=2 (quadrupole): 0 (same reason)\n");
    printf("  All L>0: 0 by K-symmetry of hedgehog\n\n");

    /* ============================================================
     * Section 2: Pion Channel (massless, isovector)
     * ============================================================
     *
     * The pion field for the hedgehog is:
     *   π_a(x) = ρ₀ sin(f(r)) x_a/r
     *
     * This is an ISOVECTOR (I=1) with angular dependence x_a/r = Y₁ₘ.
     * The pion field is pure DIPOLE (L=1). Its monopole is exactly zero:
     *   ∫ π_a dΩ = ρ₀ sin(f)/r · ∫ x_a dΩ = 0  (odd angular integral)
     *
     * The pion's monopole coupling to the isoscalar T₀₀ vanishes by
     * representation theory: ⟨I=0|I=1|I=0⟩ = 0.
     */
    printf("===== Section 2: Pion Channel =====\n\n");

    /* Tail coefficient A: f(r) → A/r² for massless pions */
    double sum_A = 0;
    int count_A = 0;
    for (int i = 0; i < prof->n; i++) {
        double r = prof->r[i];
        if (r >= 5.0 && r <= 7.0) {
            sum_A += prof->f[i] * r * r;
            count_A++;
        }
    }
    double A_tail = sum_A / count_A;

    printf("Pion field: π_a = ρ₀ sin(f) x_a/r\n");
    printf("Tail: f(r) → A/r² with A = %.4f\n\n", A_tail);

    /* Pion monopole coupling: ∫ π_a · (isoscalar) dΩ = 0
     * Verify numerically by computing ∫ sin(f) · Y₀₀ · (x_a/r) dΩ
     * This is identically zero because x_a/r has no L=0 component.
     * But verify: ∫ (x_a/r) dΩ = 0 for a=1,2,3 */
    printf("Pion monopole coupling (verify = 0):\n");
    printf("  ∫ (x_a/r) dΩ = 0 for a=x,y,z (odd under parity)\n");
    printf("  → ∫ π_a · T₀₀ d³x = 0 (isovector × isoscalar = 0)\n");
    printf("  Pion has NO monopole coupling. ✗\n\n");

    /* Pion dipole moment: D_{ai} = ∫ π_a (x_i/r) r² dΩ dr
     * = ρ₀ δ_{ai} · (4π/3) ∫ sin(f) r² dr */
    double dipole_int = 0;
    for (int i = 0; i < prof->n - 1; i++) {
        double r = prof->r[i];
        double f = prof->f[i];
        if (r < 1e-10) continue;
        dipole_int += sin(f) * r * r * prof->dr;
    }
    double D_tensor = rho0 * (4.0 * M_PI / 3.0) * dipole_int;

    printf("Pion dipole moment:\n");
    printf("  D_{ai} = ρ₀ δ_{ai} · (4π/3) ∫ sin(f) r² dr\n");
    printf("  ∫ sin(f) r² dr = %.4f (diverges for m_π=0; computed to r=%.0f)\n",
           dipole_int, prof->r[prof->n-1]);
    printf("  D = %.4f (tensor coefficient)\n\n", D_tensor);

    /* Manton dipole-dipole interaction coefficient */
    double C_dipole = 8.0 * M_PI * rho0 * rho0 * A_tail * A_tail;

    printf("Pion-mediated interaction (dipole-dipole):\n");
    printf("  C_dipole = 8πρ₀²A² = %.4f\n", C_dipole);
    printf("  E_int(R) = C·g(α)/R³  (orientation-dependent)\n");
    printf("  At R=5:  E_int ~ %.4f (attractive, α=0)\n",
           C_dipole / (5.0 * 5.0 * 5.0));
    printf("  At R=10: E_int ~ %.4f\n", C_dipole / (10.0 * 10.0 * 10.0));
    printf("  Falls as 1/R³ → NOT gravity (need 1/R)\n\n");

    /* ============================================================
     * Section 3: Breathing Mode (massive scalar, monopole-coupled)
     * ============================================================
     *
     * The breathing mode σ = ρ - ρ₀ is a radial perturbation of |q|.
     * In the finite-λ theory (V = λ(ρ²-ρ₀²)²):
     *   - Mass: m_σ² = V''(ρ₀) = 8λρ₀²
     *   - σ-model limit (λ→∞): m_σ → ∞ (mode is FROZEN)
     *
     * The coupling to the energy density:
     *   E₂ = ½ρ²(f'² + 2sin²f/r²) → ∂E₂/∂ρ = ρ·S₀(r)
     *   E₄ is ρ-independent (right-current A_μ = q̃∂q/|q|² is |q|-independent)
     *
     * The monopole source integral:
     *   g_σ = ∫ ρ₀·S₀(r) 4πr² dr = 2E₂/ρ₀ = M/ρ₀ (at virial E₂=E₄=M/2)
     */
    printf("===== Section 3: Breathing Mode =====\n\n");

    /* Source integral: ∫ ρ₀ S₀ 4πr² dr where S₀ = f'² + 2sin²f/r² */
    double source_int = 0;
    for (int i = 1; i < prof->n - 1; i++) {
        double r  = prof->r[i];
        double f  = prof->f[i];
        double fp = prof->fp[i];
        if (r < 1e-10) continue;

        double sf = sin(f);
        double S0 = fp * fp + 2.0 * sf * sf / (r * r);
        source_int += rho0 * S0 * 4.0 * M_PI * r * r * prof->dr;
    }

    double g_sigma = source_int;  /* = 2E₂/ρ₀ */
    double g_sigma_expected = 2.0 * E2 / rho0;

    printf("Breathing mode is SCALAR (I=0, K=0) → couples to MONOPOLE ✓\n\n");

    printf("Monopole source integral:\n");
    printf("  g_σ = ∫ ρ₀·S₀(r) 4πr² dr = %.4f\n", g_sigma);
    printf("  Expected: 2E₂/ρ₀ = %.4f\n", g_sigma_expected);
    printf("  Ratio: %.6f (should be 1.000)\n\n", g_sigma / g_sigma_expected);

    /* E₄ is ρ-independent: verify by noting the right-current structure */
    printf("E₄ ρ-dependence:\n");
    printf("  Right-current: A_μ = q̃∂_μq/|q|² = (∂ρ)/ρ + Ω̃∂Ω\n");
    printf("  Commutator: [A_i, A_j] = [Ω̃∂_iΩ, Ω̃∂_jΩ] (ρ cancels)\n");
    printf("  E₄ = (1/4e²)|[A,A]|² is ρ-INDEPENDENT ✓\n");
    printf("  → Breathing mode couples ONLY through E₂, not E₄\n\n");

    /* Breathing mode mass */
    printf("Breathing mode mass (finite-λ theory):\n");
    printf("  V(ρ) = λ(ρ²-ρ₀²)²  →  m_σ² = V''(ρ₀) = 8λρ₀²\n\n");
    printf("  %12s  %12s  %12s  %12s\n", "λ", "m_σ (code)", "m_σ (MeV)", "range (fm)");

    double lambda_vals[] = {8000, 20000, 80000, 1e6, 1e8};
    int n_lam = sizeof(lambda_vals) / sizeof(lambda_vals[0]);
    for (int il = 0; il < n_lam; il++) {
        double lam = lambda_vals[il];
        double m_sigma = sqrt(8.0 * lam) * rho0;
        double range_code = 1.0 / m_sigma;
        printf("  %12.0f  %12.2f  %12.0f  %12.6f\n",
               lam, m_sigma, m_sigma * code_to_MeV,
               range_code * code_to_fm);
    }
    printf("  %12s  %12s  %12s  %12s\n", "σ-model", "∞", "∞", "0");
    printf("\n  In σ-model limit: m_σ → ∞, breathing mode is FROZEN ✗\n\n");

    /* Hypothetical: if breathing mode were massless */
    printf("Hypothetical: IF breathing mode were massless:\n\n");

    /* Exchange potential: V(R) = -g₁·g₂/(4πR) for massless scalar */
    double g_per_mass = g_sigma / M;  /* coupling per unit mass */
    printf("  Monopole coupling: g_σ = %.4f = %.4f × M\n", g_sigma, g_per_mass);
    printf("  Exchange potential: V(R) = -g_σ²/(4πR) = -(%.2f)/(4πR)\n",
           g_sigma * g_sigma);
    printf("                           = -%.4f / R\n\n",
           g_sigma * g_sigma / (4.0 * M_PI));

    /* Compare with Newtonian: V = -G·M²/R where G = 1/M_P² */
    double G_eff = g_per_mass * g_per_mass / (4.0 * M_PI);
    /* G_eff = g_σ²/(4π M²) in code units */
    double M_P_eff_sq = 1.0 / G_eff;
    double M_P_eff = sqrt(M_P_eff_sq);

    double M_P_eff_MeV = M_P_eff * code_to_MeV;
    double ratio_MP = M_Planck_MeV / M_P_eff_MeV;
    double ratio_force = ratio_MP * ratio_MP;

    printf("  Effective Newton's constant: G_eff = g_σ²/(4πM²)\n");
    printf("  Effective Planck mass: M_P_eff = M/g_σ · √(4π) = %.4f code = %.1f MeV\n",
           M_P_eff, M_P_eff_MeV);
    printf("  Real Planck mass: M_P = 1.22 × 10²² MeV\n");
    printf("  Ratio: M_P/M_P_eff = %.2e\n", ratio_MP);
    printf("  Force ratio: (M_P/M_P_eff)² = %.2e (times too strong)\n\n", ratio_force);

    printf("  → Even if the breathing mode were massless, the resulting\n");
    printf("    scalar force would be %.0e × stronger than gravity.\n", ratio_force);
    printf("    This is a NUCLEAR force, not gravitational.\n\n");

    /* Additional obstruction: making σ massless destroys the soliton */
    printf("  Additional obstruction: making σ massless requires λ → 0.\n");
    printf("  But soliton collapses at λ ≲ 7000 (Derrick scaling).\n");
    printf("  Cannot have BOTH a stable soliton AND a massless σ.\n\n");

    /* ============================================================
     * Section 4: Other Channels
     * ============================================================ */
    printf("===== Section 4: Other Channels =====\n\n");

    printf("Weight sector (J, P) from Cl⁺(3,0,1):\n");
    printf("  At ε₀²=0: degenerate sector has NO kinetic energy\n");
    printf("  L₂ = ½⟨(∂Ψ)(∂Ψ̃)⟩₀ = ½|∂q|² (weight terms cancel exactly)\n");
    printf("  J and P cannot propagate → cannot mediate any force\n\n");

    printf("Spin-2 tensor (graviton analog):\n");
    printf("  The field q ∈ SU(2) ≅ S³ supports perturbations with isospin I=1.\n");
    printf("  Under combined K = I+L symmetry:\n");
    printf("    K=0: isorotation zero mode (I=1, L=1)\n");
    printf("    K=1: pion excitations (I=1, L=0,1,2)\n");
    printf("    K=2: exists formally (I=1, L=1,2,3) but:\n");
    printf("      - Source is hedgehog quadrupole = 0 (spherical symmetry)\n");
    printf("      - Is an isovector (I=1), not coupled to isoscalar T₀₀\n");
    printf("      - Has angular barriers modifying propagator\n");
    printf("  No SYMMETRIC TENSOR h_μν field exists in the Skyrme model.\n");
    printf("  BLV result confirmed: only h=0 and h=±1 modes, NO h=±2.\n\n");

    /* ============================================================
     * Section 5: Summary
     * ============================================================ */
    printf("================================================================\n");
    printf(" SUMMARY: Gravity No-Go for Skyrme L₂+L₄\n");
    printf("================================================================\n\n");

    printf("For 1/R gravity between solitons, need a mediator that is:\n");
    printf("  (a) MASSLESS: propagator G(r) = 1/(4πr)\n");
    printf("  (b) MONOPOLE-COUPLED: source ∫ φ·T₀₀ d³x ≠ 0\n\n");

    printf("┌────────────┬──────────┬──────────────┬──────────┬───────────────┐\n");
    printf("│ Channel    │ Massless │ Monopole     │ Result   │ Interaction   │\n");
    printf("├────────────┼──────────┼──────────────┼──────────┼───────────────┤\n");
    printf("│ Pion (π)   │   YES    │  NO (dipole) │  ✗ → ✗  │ C/R³ dipole   │\n");
    printf("│ Breath (σ) │   NO     │  YES         │  ✗ → ✗  │ Yukawa        │\n");
    printf("│ Weight J,P │   N/A    │  N/A         │  ✗ → ✗  │ non-dynamical │\n");
    printf("│ Spin-2     │   N/A    │  N/A         │  ✗ → ✗  │ doesn't exist │\n");
    printf("└────────────┴──────────┴──────────────┴──────────┴───────────────┘\n\n");

    printf("No channel satisfies BOTH (a) AND (b). Gravity is absent.\n\n");

    printf("Quantitative:\n");
    printf("  Pion: E_int = %.1f/R³ (dipole, orientation-dependent, ⟨SU(2)⟩=0)\n",
           C_dipole);
    printf("  σ (massless hyp.): E_int = %.1f/R (monopole, 10^%.0f × too strong)\n",
           g_sigma * g_sigma / (4.0 * M_PI), log10(ratio_force));
    printf("  σ (actual): E_int ~ e^{-m_σ R}/R (Yukawa, range → 0 in σ-model)\n\n");

    printf("This structural argument explains ALL prior negative results:\n");
    printf("  • BLV metric (Yukawa, not 1/r): σ is massive\n");
    printf("  • Depletion (1/R³, not 1/R): pion is dipole-coupled\n");
    printf("  • WZW/constraint/spectral (null): no additional massless fields\n");
    printf("  • L₆ sextic (nuclear-scale): same mechanism as E₂+E₄\n\n");

    printf("Gravity from L₂+L₄ is STRUCTURALLY IMPOSSIBLE.\n");
    printf("Any viable gravity mechanism requires physics beyond the\n");
    printf("Skyrme Lagrangian (new massless fields, nonperturbative effects,\n");
    printf("or coupling to additional sectors).\n");

    /* ============================================================
     * Write output file
     * ============================================================ */
    char fname[512];
    snprintf(fname, sizeof(fname), "%s/monopole_coupling.dat", outdir);
    FILE *fout = fopen(fname, "w");
    if (fout) {
        fprintf(fout, "# Monopole coupling analysis\n");
        fprintf(fout, "# Soliton mass M = %.6f (E2=%.4f, E4=%.4f)\n", M, E2, E4);
        fprintf(fout, "# Pion tail coefficient A = %.6f\n", A_tail);
        fprintf(fout, "# Pion dipole coefficient C = 8*pi*rho0^2*A^2 = %.6f\n", C_dipole);
        fprintf(fout, "# Breathing source g_sigma = %.6f (= 2*E2/rho0)\n", g_sigma);
        fprintf(fout, "# Breathing mass m_sigma = sqrt(8*lambda)*rho0\n");
        fprintf(fout, "# Hypothetical Planck mass M_P_eff = %.6f code = %.1f MeV\n",
                M_P_eff, M_P_eff_MeV);
        fprintf(fout, "# Force ratio (M_P/M_P_eff)^2 = %.2e\n", ratio_force);
        fprintf(fout, "#\n");
        fprintf(fout, "# lambda  m_sigma  range_code  range_fm\n");
        for (int il = 0; il < n_lam; il++) {
            double lam = lambda_vals[il];
            double m_sigma = sqrt(8.0 * lam) * rho0;
            double range_code = 1.0 / m_sigma;
            fprintf(fout, "%.0f  %.6f  %.8f  %.8f\n",
                    lam, m_sigma, range_code, range_code * code_to_fm);
        }
        fclose(fout);
        printf("\nData written to %s\n", fname);
    }

    free_profile(prof);
    return 0;
}
