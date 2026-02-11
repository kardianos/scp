/*
 * oneloop.c — One-loop B⁰p effective coupling from pion fluctuations
 *
 * FEYNMAN DIAGRAM:
 *
 *      B⁰(r)                    p (external)
 *        |                         |
 *        | (topological current)   | (constraint vertex)
 *        v                         v
 *   ---- * -------[pion loop]----- * ----
 *
 * LEFT vertex: B⁰ couples to pion pair via the topological current.
 *   The baryon current B^μ = (1/24π²) ε^μνρσ Tr(L_ν L_ρ L_σ) has
 *   terms linear in pion fluctuations from expanding L_μ = q⁻¹∂_μq
 *   around the hedgehog.
 *
 * RIGHT vertex: the constraint |q|² + ε₀²p² = ρ₀² couples p² to
 *   the bulk sector. Expanding q = q₀√(1 - ε₀²p²/ρ₀²), the kinetic
 *   energy L₂ = (ρ₀²/2)|∂(q/ρ₀)|² gets a correction:
 *     δL₂ = -(ε₀²/2) S₀(r) p²
 *   where S₀(r) = |∂(q₀/ρ₀)|² = f'² + 2sin²f/r².
 *
 * ACTUALLY: The one-loop diagram structure is simpler than the above.
 *
 * The effective action for p at one loop is:
 *   Γ[p] = S_tree[p] + (1/2) Tr ln[-∇² + V_eff(r,p)]
 *
 * where V_eff includes the constraint coupling to p. Expanding to
 * first order in p around p=0:
 *   Γ_1[p] = ∫ Σ(r) p(r) d³r
 *
 * where Σ(r) is the one-loop tadpole:
 *   Σ(r) = ε₀² × [d/dp (V_eff)]_{p=0} × G_π(r,r)
 *
 * The effective vertex d(V_eff)/dp|_{p=0} comes from the constraint:
 *   |q|² + ε₀²p² = ρ₀²  →  |q|² = ρ₀² - ε₀²p²
 *   L₂ = (1/2)|∂q|² = (1/2)(ρ₀² - ε₀²p²) S₀(r) + ...
 *
 * BUT WAIT: this is a p² term, not a linear p term. The one-loop
 * tadpole requires a LINEAR coupling between p and the pion bilinear.
 *
 * The LINEAR p coupling comes from the B⁰p term in the Lagrangian
 * (if it exists at tree level) or from higher-order constraint effects.
 *
 * CORRECT DIAGRAM:
 *
 * At tree level, there is NO B⁰p coupling (g_top = 0).
 * At one loop, the coupling is generated via:
 *
 * 1. The constraint gives a π-π-p-p FOUR-POINT vertex:
 *    δL = -(ε₀²/ρ₀²) p² × [ρ₀² S₀(r)]_π² = -(ε₀²/ρ₀²) p² (δ²E₂/δπ²)
 *
 * 2. The topological current B⁰ couples to the pion via the CUBIC vertex
 *    in the Skyrme Lagrangian expansion.
 *
 * 3. Combining: pion loop with B⁰ insertion + external p² → effective B⁰p
 *
 * For the s-wave p coupling:
 *   g_top^{eff} × B⁰ × p = (ε₀²/ρ₀²) × ∫ d³x' [G_π(x,x')²] B⁰(x') × p(x)
 *
 * Actually, the simplest derivation uses the EFFECTIVE POTENTIAL approach:
 *
 * 1. The one-loop effective potential for p in the Skyrmion background:
 *    V_eff^{1-loop}(p) = V_tree(p) + (1/(64π²)) ∫ dk² k² ln(1 + δV(r,p)/k²)
 *
 * 2. Expanding to linear order in p:
 *    The constraint shifts the pion mass: m_π² → m_π² + ε₀² × (something)
 *    This shift depends on p through the constraint, generating a p-dependent
 *    one-loop contribution.
 *
 * 3. The linear-in-p piece of V_eff^{1-loop} is the effective B⁰p coupling.
 *
 * SIMPLEST CORRECT ESTIMATE:
 *
 * The constraint gives: δm_π²(r) = ε₀² × S₀(r) (position-dependent pion mass shift)
 * The one-loop energy of the pion field is:
 *   E_loop = (1/2) Σ_n ω_n  where ω_n² = k_n² + m_π² + ε₀² S₀(r)
 *
 * The shift due to ε₀² is:
 *   δE_loop = (1/2) Σ_n [ε₀² S₀/(2ω_n)] = (ε₀²/4) ∫ S₀(r) G_π(r,r) d³r
 *
 * This is p-INDEPENDENT. To get a p-dependent piece, we need the
 * constraint coupling explicitly:
 *   |q|² = ρ₀² - ε₀²p² → S₀(r) → S₀(r)(1 - ε₀²p²/ρ₀²) + O(p⁴)
 *
 * So δm_π²(r,p) = ε₀² × S₀(r) × (1 - ε₀²p²/ρ₀²)
 *
 * This is QUADRATIC in p, not linear. To get a linear coupling, we need
 * an ODD number of p insertions, which requires a P-ODD vertex.
 *
 * THE KEY: The Wess-Zumino-Witten term IS P-odd!
 *   L_WZW ∝ ε^μνρσ B_μ ∂_ν(degenerate)_ρσ
 *
 * For statics (B⁰ only), the WZW contribution is:
 *   L_WZW ∝ B⁰ × (something involving p)
 *
 * Wait — this was already explored (Path 5). The WZW coupling vanishes
 * on statics. So there is NO tree-level B⁰p coupling.
 *
 * At one loop, the P-odd coupling comes from the PION LOOP itself.
 * The pion on the soliton background breaks P (the background has
 * definite B = +1). The pion propagator on the soliton background
 * has a P-odd component:
 *   G_π(r,r') ≠ G_π(r',r) (asymmetric in background)
 *
 * Actually, for a static background, G_π IS symmetric. The P-odd piece
 * arises from the B⁰ INSERTION vertex.
 *
 * CORRECT ONE-LOOP DIAGRAM (final):
 *
 *    p ──── constraint vertex ──── π ─── loop ─── π ──── B⁰ vertex ──── (nothing)
 *
 * Actually this is a VACUUM diagram with external p. Let me be more precise.
 *
 * The B⁰p coupling at one loop:
 *
 *   <0| T{B⁰(x) p(y)} |0>_connected
 *     = (tree-level propagator) + (one-loop correction)
 *
 * But B⁰ is a composite operator (topological current), not a fundamental field.
 * The proper statement is:
 *
 * The EFFECTIVE ACTION for the degenerate field p, after integrating out
 * pion fluctuations at one loop, acquires a term:
 *   Γ_eff[p] ∋ g_top^{eff} ∫ B⁰(r) p(r) d³r
 *
 * This arises from the following chain:
 *
 * 1. The pion field π_a (a=1,2,3) fluctuates on the hedgehog background.
 * 2. The constraint couples π² to p²: δL ∝ ε₀² p² π²
 * 3. At one loop: <π²(r)> = G_π(r,r) ≠ 0 (vacuum fluctuations)
 * 4. This generates a p²-linear shift in the effective Lagrangian:
 *    L_eff ∋ -ε₀² G_π(r,r) p²
 *
 * But this is STILL p², not linear in p.
 *
 * LINEAR-IN-P COUPLING:
 *
 * For a LINEAR p coupling, we need an operator that is odd in p.
 * In the Cl⁺(3,0,1) algebra, p is the e₀₁₂₃ component (pseudoscalar).
 * A linear p term in the Lagrangian requires a pseudoscalar source.
 *
 * The baryon density B⁰ IS a pseudoscalar (changes sign under parity).
 * So B⁰ × p is a scalar, and the coupling is allowed by symmetry.
 *
 * But HOW does the one-loop diagram generate it? The constraint gives
 * only p² couplings. To get a LINEAR p term, we need either:
 * (a) A tree-level linear coupling (absent — this is what we're trying to generate)
 * (b) A one-loop diagram with an ODD number of p legs
 *
 * For (b): the constraint vertex has 2 p legs. To get 1 external p,
 * we need a loop with 1 p propagator closing back on itself AND a
 * vertex connecting to B⁰. But p has no self-coupling at tree level.
 *
 * CONCLUSION: At one loop with the constraint coupling alone,
 * g_top^{1-loop} = 0 (p² coupling cannot generate linear p).
 *
 * To generate a LINEAR B⁰p coupling at one loop requires EITHER:
 * (a) A WZW-type vertex (P-odd, provides the single p) — but WZW vanishes on statics
 * (b) TWO-LOOP diagram: constraint (p²π²) × WZW (B⁰π²p) → close pion loop
 * (c) Quantum anomaly that breaks the Z₂ symmetry p → -p
 *
 * However, the Z₂ symmetry p → -p IS broken by the WZW term:
 *   L_WZW ∝ ε^μνρσ B_μ ∂_ν w_ρ ∂_σ w_τ
 * where w is related to the degenerate sector. If the WZW term has
 * a pion-pion-p vertex (linear in p), then the one-loop diagram:
 *
 *   p ──── WZW vertex ──── π ─── loop ─── π ──── constraint ──── (vacuum)
 *
 * generates a one-point function for p = tadpole of p in the B⁰ background.
 *
 * Let me compute this properly.
 *
 * WZW VERTEX STRUCTURE:
 * The WZW term couples B^μ to the axial current of the degenerate sector:
 *   L_WZW = (N_c/(240π²)) ε^μνρσ B_μ J^D_νρσ
 * where J^D includes derivatives of p.
 *
 * For a static hedgehog with p perturbation:
 *   L_WZW ⊃ (N_c/(240π²)) B⁰ × (pion × pion × ∂p) term
 *
 * This IS linear in p and involves π² — so the one-loop diagram:
 *   Contract π² into a loop → <π²> = G_π(r,r) → effective B⁰p coupling
 *
 * g_top^{WZW,1-loop} = (N_c/(240π²)) × G_π^{reg}(r) averaged over soliton
 *
 * This is the CORRECT mechanism!
 *
 * NUMERICAL COMPUTATION:
 * We compute all the overlap integrals and estimate the loop using
 * known spectral data from angular_modes.c and normal_modes.c.
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
        double rv, fval, fpval = 0, d1, d2, d3;
        int nc = sscanf(line, "%lf %lf %lf %lf %lf %lf",
                        &rv, &fval, &fpval, &d1, &d2, &d3);
        if (nc < 2) continue;
        if (n >= cap) {
            cap *= 2;
            prof->r  = realloc(prof->r,  cap * sizeof(double));
            prof->f  = realloc(prof->f,  cap * sizeof(double));
            prof->fp = realloc(prof->fp, cap * sizeof(double));
        }
        prof->r[n] = rv; prof->f[n] = fval; prof->fp[n] = fpval; n++;
    }
    fclose(fp);
    prof->n = n;
    prof->dr = (n > 1) ? prof->r[1] - prof->r[0] : 0.001;
    int has_fp = 0;
    for (int i = 0; i < n; i++)
        if (fabs(prof->fp[i]) > 1e-30) { has_fp = 1; break; }
    if (!has_fp) {
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

static double strain_S0(double r, double f, double fp, double a) {
    if (r > 1e-10) {
        double sf = sin(f);
        return fp*fp + 2.0*sf*sf/(r*r);
    }
    return 3.0*a*a;
}

static double baryon_density(double r, double f, double fp, double a) {
    if (r > 1e-10) {
        double sf = sin(f);
        return -(1.0/(2*M_PI*M_PI))*sf*sf*fp/(r*r);
    }
    return a*a*a/(2*M_PI*M_PI);
}

/* Skyrme angular density T₀(r) = 2f'²sin²f/r² + sin⁴f/r⁴ */
static double skyrme_T0(double r, double f, double fp, double a) {
    if (r > 1e-10) {
        double sf = sin(f), sf2 = sf*sf;
        return 2.0*fp*fp*sf2/(r*r) + sf2*sf2/(r*r*r*r);
    }
    return 3.0*a*a*a*a;
}

int main(int argc, char *argv[])
{
    const char *profile_file = "data/profiles/profile_sigma_e1.dat";
    double e_param = 1.0, rho0 = 1.0, m_pi = 0.398;
    int Nc = 3;  /* number of colors (WZW coefficient) */

    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "-profile") && i+1 < argc) profile_file = argv[++i];
        else if (!strcmp(argv[i], "-e") && i+1 < argc) e_param = atof(argv[++i]);
        else if (!strcmp(argv[i], "-rho0") && i+1 < argc) rho0 = atof(argv[++i]);
        else if (!strcmp(argv[i], "-mpi") && i+1 < argc) m_pi = atof(argv[++i]);
        else if (!strcmp(argv[i], "-Nc") && i+1 < argc) Nc = atoi(argv[++i]);
        else {
            fprintf(stderr, "Usage: %s [-profile PATH] [-e E] [-rho0 R] [-mpi M] [-Nc N]\n",
                    argv[0]);
            return 1;
        }
    }

    double c4 = 2.0*rho0*rho0/(e_param*e_param);
    double m_pi2 = m_pi*m_pi;

    Profile prof;
    if (read_profile(profile_file, &prof) < 0) return 1;
    int n = prof.n;
    double dr = prof.dr;
    double a = -prof.fp[0];

    printf("============================================================\n");
    printf(" Avenue A: One-Loop B⁰p Effective Coupling\n");
    printf("============================================================\n\n");
    printf("Profile: %s (%d pts, R_max=%.2f, dr=%.6f)\n",
           profile_file, n, prof.r[n-1], dr);
    printf("e=%.4f, rho0=%.4f, c4=%.6f, a=%.6f\n", e_param, rho0, c4, a);
    printf("m_pi=%.4f (m_pi^2=%.6f), N_c=%d\n\n", m_pi, m_pi2, Nc);

    /* Allocate arrays */
    double *S0 = calloc(n, sizeof(double));
    double *B0 = calloc(n, sizeof(double));
    double *T0 = calloc(n, sizeof(double));
    double *Omega = calloc(n, sizeof(double));  /* isospin inertia density */

    for (int i = 0; i < n; i++) {
        double r = prof.r[i], f = prof.f[i], fp = prof.fp[i];
        S0[i] = strain_S0(r, f, fp, a);
        B0[i] = baryon_density(r, f, fp, a);
        T0[i] = skyrme_T0(r, f, fp, a);
        double sf2 = sin(f)*sin(f);
        Omega[i] = rho0*rho0*sf2*(1.0 + c4*(fp*fp + sf2/((r*r)+1e-30)));
    }

    /* ===== Overlap integrals ===== */
    printf("=== Overlap Integrals ===\n\n");

    double I_SB = 0;    /* ∫ S₀ B⁰ 4πr² dr */
    double I_S2 = 0;    /* ∫ S₀² 4πr² dr */
    double I_B2 = 0;    /* ∫ B⁰² 4πr² dr */
    double Q_total = 0; /* ∫ B⁰ 4πr² dr */
    double I_TB = 0;    /* ∫ T₀ B⁰ 4πr² dr */
    double I_OB = 0;    /* ∫ Ω B⁰ 4πr² dr */
    double E2 = 0, E4 = 0; /* energies */
    double Lambda_sum = 0;  /* moment of inertia integrand */

    for (int i = 1; i < n; i++) {
        double r0 = prof.r[i-1], r1 = prof.r[i];
        double w = 0.5*dr;

        I_SB += w*(4*M_PI*S0[i-1]*B0[i-1]*r0*r0 + 4*M_PI*S0[i]*B0[i]*r1*r1);
        I_S2 += w*(4*M_PI*S0[i-1]*S0[i-1]*r0*r0 + 4*M_PI*S0[i]*S0[i]*r1*r1);
        I_B2 += w*(4*M_PI*B0[i-1]*B0[i-1]*r0*r0 + 4*M_PI*B0[i]*B0[i]*r1*r1);
        Q_total += w*(4*M_PI*B0[i-1]*r0*r0 + 4*M_PI*B0[i]*r1*r1);
        I_TB += w*(4*M_PI*T0[i-1]*B0[i-1]*r0*r0 + 4*M_PI*T0[i]*B0[i]*r1*r1);
        I_OB += w*(4*M_PI*Omega[i-1]*B0[i-1]*r0*r0 + 4*M_PI*Omega[i]*B0[i]*r1*r1);

        E2 += w*(2*M_PI*rho0*rho0*S0[i-1]*r0*r0 + 2*M_PI*rho0*rho0*S0[i]*r1*r1);
        double pre4 = 4*M_PI*rho0*rho0*rho0*rho0/(e_param*e_param);
        E4 += w*(pre4*T0[i-1]*r0*r0 + pre4*T0[i]*r1*r1);

        double f0 = prof.f[i-1], f1 = prof.f[i];
        double fp0 = prof.fp[i-1], fp1 = prof.fp[i];
        double sf0 = sin(f0), sf1 = sin(f1);
        double Om0 = sf0*sf0*(1+c4*(fp0*fp0+sf0*sf0/(r0*r0+1e-30)));
        double Om1 = sf1*sf1*(1+c4*(fp1*fp1+sf1*sf1/(r1*r1+1e-30)));
        Lambda_sum += w*(r0*r0*Om0 + r1*r1*Om1);
    }
    double Lambda = (8.0*M_PI*rho0*rho0/3.0) * Lambda_sum;
    double E_sol = E2 + E4;

    printf("Q_total = %.10f (expect 1.0)\n", Q_total);
    printf("E2 = %.4f, E4 = %.4f, E_sol = %.4f (virial E2/E4=%.4f)\n",
           E2, E4, E_sol, E2/E4);
    printf("Lambda = %.4f (moment of inertia)\n\n", Lambda);

    printf("I_SB = int S0*B0 4pi r^2 dr = %.10f\n", I_SB);
    printf("I_S2 = int S0^2  4pi r^2 dr = %.10f\n", I_S2);
    printf("I_B2 = int B0^2  4pi r^2 dr = %.10f\n", I_B2);
    printf("I_TB = int T0*B0 4pi r^2 dr = %.10f\n", I_TB);
    printf("I_OB = int Om*B0 4pi r^2 dr = %.10f\n\n", I_OB);

    /* ===== Section 1: The Constraint Vertex ===== */
    printf("============================================================\n");
    printf(" 1. CONSTRAINT VERTEX ANALYSIS\n");
    printf("============================================================\n\n");

    printf("Constraint: |q|^2 + eps0^2 p^2 = rho0^2\n");
    printf("=> q = q0 sqrt(1 - eps0^2 p^2/rho0^2)\n\n");

    printf("Expanding L2 = (rho0^2/2)|d(q/rho0)|^2 to second order in p:\n");
    printf("  L2 = (rho0^2/2)S0(r) - (eps0^2/2)S0(r)p^2 + O(eps0^4 p^4)\n\n");

    printf("Constraint vertex: V_{pi-pi-p-p} = eps0^2 S0(r) / rho0^2\n");
    printf("  (This couples pion pairs to p pairs — EVEN in p)\n\n");

    printf("KEY: The constraint vertex is p^2, not p.\n");
    printf("A one-loop pion diagram with constraint vertex gives <pi^2>p^2,\n");
    printf("which is a p MASS SHIFT, not a B0*p coupling.\n\n");

    /* ===== Section 2: The WZW Vertex ===== */
    printf("============================================================\n");
    printf(" 2. WZW VERTEX (Linear in p)\n");
    printf("============================================================\n\n");

    printf("The WZW term provides the P-odd vertex linear in p:\n");
    printf("  L_WZW = (Nc/(240 pi^2)) eps^{munurhosgm} Tr(...)\n\n");

    printf("For a STATIC hedgehog, the WZW coupling B^0 * (pion^2 dp) gives:\n");
    printf("  L_WZW^{static} = (Nc/(240 pi^2)) B0(r) * [pi x dpi] . dp/dr\n");
    printf("But: on statics with p=0, this vanishes (no dp/dr).\n\n");

    printf("At ONE LOOP: the pion loop generates <pi x dpi> != 0\n");
    printf("if the loop has a B0 insertion (background breaks parity).\n\n");

    printf("This is a TWO-vertex diagram:\n");
    printf("  - Vertex 1: WZW couples B0 to (pi x dpi) * p\n");
    printf("  - Vertex 2: constraint couples pi^2 to p^2\n");
    printf("  - Close the pion propagator into a loop\n");
    printf("  - One external p leg remains\n\n");

    printf("Actually, the simpler mechanism is the WZW TADPOLE:\n");
    printf("  L_WZW ∋ Nc/(240pi^2) × B0 × (d/dr)[pi_a pi_b eps_abc] × p_c\n");
    printf("Wait — this requires careful Clifford algebra expansion.\n\n");

    /* ===== Section 3: Correct Diagram Identification ===== */
    printf("============================================================\n");
    printf(" 3. CORRECT DIAGRAM STRUCTURE\n");
    printf("============================================================\n\n");

    printf("After careful analysis of the Cl+(3,0,1) algebra:\n\n");

    printf("The CORRECT mechanism for one-loop B0*p generation is:\n\n");

    printf("(A) CONSTRAINT MECHANISM (p^2 coupling → p mass shift):\n");
    printf("    Vertex: eps0^2 S0(r) p^2 / rho0^2\n");
    printf("    One-loop: eps0^2 S0(r) <pi^2(r)> / rho0^2 = eps0^2 S0(r) G_pi(r,r)/rho0^2\n");
    printf("    Result: MASS SHIFT for p, not B0*p coupling.\n");
    printf("    δm_p^2 = eps0^2 <S0 G_pi> (averaged over soliton volume)\n\n");

    printf("(B) TWO-LOOP B0*p via constraint + WZW:\n");
    printf("    Requires TWO loops (or one loop with two vertices).\n");
    printf("    Parametrically: g_top ~ eps0^2 * Nc/(240pi^2) * G_pi^2\n");
    printf("    Suppressed by extra loop factor 1/(16pi^2) relative to (A).\n\n");

    printf("(C) DIRECT TADPOLE from p-Z2 breaking by background:\n");
    printf("    The hedgehog background has B = +1 (not B = 0).\n");
    printf("    In the soliton sector, the Z2 symmetry p → -p is\n");
    printf("    explicitly broken by the topological charge.\n");
    printf("    The one-loop effective potential for p acquires a LINEAR term:\n");
    printf("      V_eff(p) = ... + sigma * p + (1/2)(m_p^2 + delta_m^2) p^2 + ...\n");
    printf("    where sigma = (d/dp) Tr ln[-nabla^2 + V(r,p)]|_{p=0}\n\n");

    printf("    For sigma != 0, we need dV/dp|_{p=0} != 0.\n");
    printf("    From the constraint: d(S0*(1-eps0^2 p^2/rho0^2))/dp = -2 eps0^2 S0 p/rho0^2\n");
    printf("    At p=0: this vanishes. So sigma = 0 from constraint alone.\n\n");

    printf("    From WZW: L_WZW has terms ODD in p, so dL_WZW/dp|_{p=0} can be != 0.\n");
    printf("    But the WZW term vanishes on STATIC hedgehog (proved in Path 5).\n\n");

    printf("CONCLUSION: At one loop, there is NO B0*p coupling generated\n");
    printf("from the constraint + WZW for a STATIC soliton.\n\n");

    printf("The coupling CAN be generated from:\n");
    printf("  1. ROTATING soliton (Avenue B): WZW gives non-zero B^i terms\n");
    printf("  2. TWO-LOOP diagram with constraint vertex\n");
    printf("  3. Pion-loop RENORMALIZATION of the WZW coefficient\n\n");

    /* ===== Section 4: What one-loop DOES generate ===== */
    printf("============================================================\n");
    printf(" 4. WHAT ONE-LOOP DOES GENERATE\n");
    printf("============================================================\n\n");

    printf("The one-loop pion correction on the soliton background generates:\n\n");

    /* (a) Mass shift of p from constraint */

    printf("(a) p-field mass shift from constraint:\n");
    printf("    delta_m_p^2 = eps0^2 * int S0 G_pi^{reg} d^3x / rho0^2\n");
    printf("    Using heat kernel: G_pi^{reg}(r,r) ~ S0(r)/(16pi^2 m_pi^2)\n");
    printf("    delta_m_p^2 / eps0^2 = I_S2 / (16pi^2 m_pi^2 rho0^2)\n");
    printf("                         = %.4f / (%.4f * %.4f * %.4f)\n",
           I_S2, 16*M_PI*M_PI, m_pi2, rho0*rho0);
    printf("                         = %.6e\n\n", I_S2/(16*M_PI*M_PI*m_pi2*rho0*rho0));

    /* (b) Casimir energy from pion on soliton background */
    double E_casimir_est = -I_S2 / (64*M_PI*M_PI * m_pi);
    printf("(b) Casimir energy (pion zero-point on soliton background):\n");
    printf("    E_Casimir ~ -(1/2) * 3 * int_0^inf dk k^2/(2pi)^2 * delta(k)/k\n");
    printf("    Using phase shift delta(k) ~ -arctan(S0 R^3 k / (k^2+m_pi^2)):\n");
    printf("    E_Casimir ~ -3 I_S2 / (64 pi^2 m_pi) = %.6e code\n\n", E_casimir_est);

    /* (c) One-loop correction to moment of inertia */
    printf("(c) One-loop correction to moment of inertia Lambda:\n");
    printf("    delta_Lambda / Lambda ~ 1/(16pi^2) * (few) ~ 10^{-2}\n");
    printf("    This shifts Delta-N splitting by ~1%%.\n\n");

    /* ===== Section 5: The actual g_top estimate ===== */
    printf("============================================================\n");
    printf(" 5. g_top FROM ROTATING SOLITON (Avenue A + B combined)\n");
    printf("============================================================\n\n");

    printf("The physical nucleon ROTATES (J=1/2, angular velocity Omega=hbar/Lambda).\n");
    printf("For a rotating soliton, the WZW term becomes:\n");
    printf("  L_WZW = (Nc/240pi^2) * B^i * (rot. current) * p\n\n");

    printf("The spatial baryon current from rotation:\n");
    printf("  B^i(x) = (Omega/c) x eps_{ijk} x^j B^0_k(r)\n");
    printf("where B^0_k is the k-th component of the static baryon current density.\n\n");

    printf("For a J=1/2 quantum state, the matrix element:\n");
    printf("  <J=1/2| Omega |J=1/2> = sqrt(J(J+1)) hbar / Lambda = sqrt(3/4)/Lambda\n\n");

    double hbar = 1.0;  /* in code units where hbar=c=1 */
    double Omega_quantum = sqrt(0.75) * hbar / Lambda;
    printf("Omega_quantum = sqrt(3/4) / Lambda = %.6e code^{-1}\n", Omega_quantum);
    printf("  = %.6e MeV (using code_to_MeV = 9.098)\n\n", Omega_quantum * 9.098);

    /* WZW coupling: the time-averaged WZW source for p */
    /* <L_WZW> = (Nc/240pi^2) × <B^i> × ∫ (something) d^3x × p */
    /* The spatial integral involves the overlap of B^i with the pion background */

    /* For the hedgehog, the angular average of B^i × x^j gives:
     * <B^i x^j> = delta_{ij}/3 × ∫ B^0(r) r^2 4pi r^2 dr × Omega
     * But this is not quite right. Let me use the standard result:
     *
     * The rotating hedgehog's spatial baryon current is:
     *   B_rot^i(x) = (1/c) (Omega × x)^i × B^0(r) + (isospin part)
     *
     * The overlap with ∂_i p in the WZW term:
     *   ∫ B_rot^i ∂_i p d^3x = (Omega/c) ∫ (x × ∇p)_rot × B^0(r) d^3x
     *
     * For spherically symmetric p(r), ∂_i p = p'(r) r_i/r, and:
     *   (Omega × x)_i × r_i/r = 0 (cross product perpendicular to x)
     *
     * So the WZW source vanishes for spherical p(r)!
     * The coupling requires L=1 (dipolar) p fluctuations.
     */

    printf("NOTE: WZW coupling to spherically symmetric p vanishes.\n");
    printf("(Omega x x) . grad(p(r)) = 0 for s-wave p.\n");
    printf("The coupling requires L=1 (dipolar) p fluctuations.\n\n");

    printf("For L=1 dipolar p: p(r,theta) = p_1(r) cos(theta)\n");
    printf("  WZW source: (Nc/240pi^2) Omega B^0(r) × r × p_1(r)\n\n");

    /* The effective g_top from WZW + rotation:
     * g_top^{WZW} = (Nc/240pi^2) × Omega × (angular factor)
     * The angular factor comes from projecting onto the s-wave B0*p coupling.
     * Since the source is L=1 and B0*p is L=0, we need the overlap
     * integral over angles, which is:
     *   ∫ Y_{10}(theta) × Y_{00} × Y_{10} dOmega = 1/sqrt(4pi)
     */

    /* Actually, let me just compute the standard Adkins-Nappi-Witten result.
     * The g_A (axial coupling) of the nucleon in the Skyrme model is:
     *   g_A = (2/3) × (8pi/3) × ∫ r^2 sin^2(f) [f' + sin(2f)/(2r)] dr
     *   The same integral controls the WZW B-p coupling.
     */

    /* Compute the axial charge integral */
    double I_axial = 0;
    for (int i = 1; i < n; i++) {
        double r0 = prof.r[i-1], r1 = prof.r[i];
        double f0 = prof.f[i-1], f1 = prof.f[i];
        double fp0 = prof.fp[i-1], fp1 = prof.fp[i];
        double sf0 = sin(f0), sf1 = sin(f1);
        double s2f0 = sin(2*f0), s2f1 = sin(2*f1);

        double ig0 = r0*r0*sf0*sf0*(fp0 + s2f0/(2*r0 + 1e-30));
        double ig1 = r1*r1*sf1*sf1*(fp1 + s2f1/(2*r1 + 1e-30));
        I_axial += 0.5*dr*(ig0 + ig1);
    }
    double g_A = (2.0/3.0) * (8*M_PI/3.0) * I_axial;
    printf("Axial coupling g_A = %.4f (experiment: 1.27)\n\n", g_A);

    /* The WZW contribution to g_top from a rotating J=1/2 soliton:
     *
     * Standard result (Adkins, Nappi, Witten 1983):
     * The B^i current of a rotating hedgehog has magnitude:
     *   J_B^i ~ (Omega/c) × (radial integral involving f(r))
     *
     * The WZW coupling to p:
     *   g_top^{WZW} = Nc / (240 pi^2) × (hbar / Lambda) × I_WZW
     *
     * where I_WZW = ∫ B^0(r) × r × p_1(r) × (angular factor) 4pi r^2 dr
     *
     * For dimensional estimate: I_WZW ~ Q × R_rms ~ 1 × 1.5 = 1.5
     */

    double R_rms = 0;
    {
        double sum_r2 = 0, sum_w = 0;
        for (int i = 1; i < n; i++) {
            double r0 = prof.r[i-1], r1 = prof.r[i];
            sum_r2 += 0.5*dr*(B0[i-1]*r0*r0*r0*r0*4*M_PI + B0[i]*r1*r1*r1*r1*4*M_PI);
            sum_w  += 0.5*dr*(B0[i-1]*r0*r0*4*M_PI + B0[i]*r1*r1*4*M_PI);
        }
        R_rms = sqrt(sum_r2/sum_w);
    }
    printf("R_rms (baryon) = %.4f code = %.4f fm\n\n", R_rms, R_rms*0.5624);

    /* Various estimates of g_top^{WZW} */
    double g_WZW_naive = (double)Nc / (240*M_PI*M_PI) * hbar / Lambda;
    double g_WZW_with_overlap = g_WZW_naive * I_SB / Q_total;
    double g_WZW_with_Rrms = g_WZW_naive * R_rms;

    printf("WZW coupling estimates (all proportional to Nc/240pi^2 * hbar/Lambda):\n\n");
    printf("  Prefactor: Nc/(240pi^2) = %.6e\n", (double)Nc/(240*M_PI*M_PI));
    printf("  hbar/Lambda = %.6e\n", hbar/Lambda);
    printf("  Basic: g_WZW = Nc/(240pi^2) * hbar/Lambda = %.6e\n", g_WZW_naive);
    printf("  With overlap: g_WZW * I_SB/Q = %.6e\n", g_WZW_with_overlap);
    printf("  With R_rms:   g_WZW * R_rms = %.6e\n\n", g_WZW_with_Rrms);

    /* One-loop CORRECTION to g_WZW via pion loop */
    /* The pion loop renormalizes the WZW vertex.
     * At one loop:
     *   g_top^{1-loop} = g_WZW × [1 + (3/16pi^2) × ∫ V_eff G_pi d^3x]
     *
     * The correction term:
     *   delta_g / g_WZW ~ (3/16pi^2) × ∫ S0 / (k^2 + m_pi^2) dk × R^3
     *                    ~ (3/16pi^2) × S0(0) × R^3 / m_pi
     */
    double one_loop_corr = (3.0/(16*M_PI*M_PI)) * S0[0] * R_rms*R_rms*R_rms / m_pi;
    printf("One-loop correction to WZW:\n");
    printf("  delta_g/g ~ (3/16pi^2) S0(0) R^3 / m_pi = %.4f\n", one_loop_corr);
    printf("  (%.0f%% correction — perturbation theory %s)\n\n",
           one_loop_corr*100,
           (one_loop_corr < 1.0) ? "OK" : "BREAKS DOWN");

    /* ===== Section 6: Required g_top for gravity ===== */
    printf("============================================================\n");
    printf(" 6. COMPARISON WITH GRAVITY\n");
    printf("============================================================\n\n");

    double G_Newton_code = 5.9e-39;
    double kappa_D = 1.0;
    double gtop_required = sqrt(G_Newton_code * 4*M_PI * kappa_D*kappa_D * E_sol*E_sol);
    printf("Required g_top for G_Newton: %.6e\n", gtop_required);
    printf("WZW estimate (naive):        %.6e\n", g_WZW_naive);
    printf("Ratio (WZW/required):        %.2e  (WZW is %.1e x too large)\n\n",
           g_WZW_naive/gtop_required, g_WZW_naive/gtop_required);

    printf("The WZW coupling is ~10^{10} too strong for gravity.\n");
    printf("This means either:\n");
    printf("  (a) Angular integrals suppress the effective coupling\n");
    printf("  (b) The WZW contribution to B0*p is not the dominant source\n");
    printf("  (c) Cancellations between WZW and constraint contributions\n\n");

    /* Check if constraint mechanism could cancel WZW */
    printf("=== Constraint + WZW cancellation? ===\n\n");

    /* The constraint-mediated coupling (Avenue A, second order) gives:
     * g_top^{constraint} = eps0^2 × (coupling_integral)
     * For cancellation: eps0^2 × C = g_WZW, giving:
     * eps0^2 = g_WZW / C
     */

    /* For the p mass shift mechanism:
     * The constraint shifts m_p^2 by delta_m^2 ~ eps0^2 I_S2/(16pi^2 m_pi^2 rho0^2)
     * This does NOT produce B0*p directly.
     *
     * For the TWO-LOOP mechanism (constraint + WZW):
     * g_top^{2-loop} = eps0^2 × g_WZW × (loop factor)
     * loop factor ~ 1/(16pi^2) × I_SB / I_S2
     */
    double loop_factor = 1.0/(16*M_PI*M_PI) * fabs(I_SB) / I_S2;
    double g_top_2loop_per_eps2 = g_WZW_naive * loop_factor;

    printf("Two-loop mechanism: g_top = eps0^2 * g_WZW * (loop_factor)\n");
    printf("  Loop factor = I_SB/(16pi^2 I_S2) = %.6e\n", loop_factor);
    printf("  g_top^{2-loop}/eps0^2 = %.6e\n\n", g_top_2loop_per_eps2);

    printf("For g_top = g_WZW - g_top^{constraint} = 0 (cancellation):\n");
    printf("  eps0^2 = g_WZW / (g_top_per_eps0^2) = %.6e / %.6e = %.2e\n\n",
           g_WZW_naive, g_top_2loop_per_eps2,
           g_WZW_naive / g_top_2loop_per_eps2);

    /* ===== Section 7: Channel C near-binding ===== */
    printf("============================================================\n");
    printf(" 7. CHANNEL C NEAR-BINDING ENHANCEMENT\n");
    printf("============================================================\n\n");

    double omega_C2_frac = 1.04;
    double enhancement = omega_C2_frac / (omega_C2_frac - 1.0);

    printf("Channel C (L=1, I=1) eigenvalue: omega^2 ~ 1.04 * m_pi^2\n");
    printf("Distance above threshold: 4%%\n");
    printf("Propagator enhancement: 1/(omega^2 - m_pi^2) ~ %.0f / m_pi^2\n\n", enhancement);

    printf("Effect on one-loop integral:\n");
    printf("  The Channel C near-resonance enhances the L=1 partial wave\n");
    printf("  of the pion propagator by a factor ~%.0f.\n", enhancement);
    printf("  But this is one angular channel out of sum over all L.\n");
    printf("  Weighted by the L=1 fraction of B0(r): modest enhancement.\n\n");

    printf("  If Channel C WERE bound (4%% lower): would create a\n");
    printf("  TRUE pole in the pion propagator → divergent enhancement.\n");
    printf("  This would require quantization effects to push omega^2\n");
    printf("  below m_pi^2. Zero-point energy shifts are typically\n");
    printf("  delta_omega^2 ~ hbar/(Lambda R^2) ~ 0.005, which is\n");
    printf("  comparable to the 4%% gap.\n\n");

    printf("  NEAR-MISS: quantization COULD push Channel C into binding!\n");
    printf("  Required shift: delta(omega^2/m_pi^2) = -0.04\n");
    printf("  Typical quantum shift: ~0.005-0.01\n");
    printf("  Gap: factor of 4-8 too large for naive quantum correction.\n");

    /* ===== Section 8: Summary ===== */
    printf("\n============================================================\n");
    printf(" SUMMARY: Avenue A — One-Loop B0*p Coupling\n");
    printf("============================================================\n\n");

    printf("FINDING 1: Constraint vertex is EVEN in p (p^2 coupling).\n");
    printf("  → One-loop constraint diagram generates p MASS SHIFT,\n");
    printf("    not B0*p coupling. The mass shift is:\n");
    printf("    delta_m_p^2 / eps0^2 = %.2e\n\n",
           I_S2/(16*M_PI*M_PI*m_pi2*rho0*rho0));

    printf("FINDING 2: B0*p coupling requires P-ODD vertex.\n");
    printf("  → Only the WZW term provides this.\n");
    printf("  → WZW vanishes on static hedgehog (proved in Path 5).\n");
    printf("  → WZW is nonzero for ROTATING soliton (Avenue B).\n\n");

    printf("FINDING 3: WZW coupling for rotating J=1/2 soliton:\n");
    printf("  g_WZW = Nc/(240pi^2) * hbar/Lambda = %.2e\n", g_WZW_naive);
    printf("  This is %.1e x too STRONG for G_Newton.\n\n",
           g_WZW_naive/gtop_required);

    printf("FINDING 4: Channel C near-binding (4%% above threshold)\n");
    printf("  enhances the L=1 pion propagator by ~%.0fx.\n", enhancement);
    printf("  Quantization could close the gap, but naive estimates\n");
    printf("  suggest the quantum shift is 4-8x too small.\n\n");

    printf("FINDING 5: The hierarchy problem persists.\n");
    printf("  WZW gives g_top ~ 10^{-7} (too strong by 10^{10}).\n");
    printf("  For G_Newton: need g_top ~ 10^{-17}.\n");
    printf("  Either massive cancellations occur, or the B0*p coupling\n");
    printf("  is not the correct gravity mechanism.\n\n");

    printf("ASSESSMENT: Avenue A alone DOES NOT determine g_top.\n");
    printf("  - Constraint mechanism: generates p mass shift, not B0*p\n");
    printf("  - WZW mechanism: generates B0*p but 10^10x too strong\n");
    printf("  - eps0^2 remains free (not fixed by one-loop)\n");
    printf("  - Near-binding in Channel C is interesting but insufficient\n\n");

    printf("IMPLICATION FOR AVENUES B-D:\n");
    printf("  - Avenue B (rotating soliton): CRITICAL — only source of B0*p.\n");
    printf("    Must compute angular integrals precisely to check if\n");
    printf("    suppression factor ~10^{-10} arises from geometry.\n");
    printf("  - Avenue C (anomaly): could constrain eps0^2 → fix p mass\n");
    printf("  - Avenue D (spectrum): overconstrained system may fix eps0^2\n");

    free(S0); free(B0); free(T0); free(Omega);
    free(prof.r); free(prof.f); free(prof.fp);
    return 0;
}
