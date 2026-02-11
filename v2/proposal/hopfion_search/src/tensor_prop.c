/*
 * tensor_prop.c — Scalar propagator on anisotropic BLV metric
 *
 * Tests whether Path 3 (scalar 1/r from B⁰p coupling) combined with
 * Path 4 (anisotropic BLV metric from non-hedgehog topology) can produce
 * effective spin-2 gravity at long range.
 *
 * METHOD: Born approximation for Green's function on perturbed metric.
 *
 * The scalar p satisfies: g^{ij}_eff ∂_i∂_j p = source
 * On hedgehog: g^{ij}_eff = 2δ^{ij} (flat, isotropic)
 * On hopfion:  g^{ij}_eff = 2δ^{ij} + δg^{ij}(r)  (anisotropic)
 *
 * The Green's function G = G₀ + δG where:
 *
 *   δG(r₁, r₂) = -(1/2) ∫ G₀(r₁,r') δg^{ij}(r') ∂'_i∂'_j G₀(r',r₂) d³r'
 *
 * For two solitons at separation D, the tensor correction to the interaction is:
 *
 *   δU/U₀ = -3Q(3cos²θ - 1)/(4πD²)
 *
 * where Q = 4π ∫₀^∞ Δ(r) r dr is the "metric quadrupole moment" and
 * Δ(r) is the traceless anisotropy amplitude.
 *
 * KEY RESULT: Tensor correction goes as 1/D² relative to scalar monopole.
 * In absolute terms: δU ~ 1/D³ (quadrupole), while U₀ ~ 1/D (monopole).
 * True spin-2 gravity has tensor structure IN the 1/D monopole (double light
 * deflection). Our model cannot reproduce this.
 *
 * PHYSICS DERIVATION:
 *
 * Wave equation: (2δ^{ij} + δg^{ij})∂_i∂_j p = J(r)
 *
 * Born approximation (first order in δg):
 *   p = p₀ + δp, where ∇²p₀ = J/2 and
 *   ∇²δp = -(1/2)δg^{ij} ∂_i∂_j p₀
 *
 * For axially symmetric (along ẑ) traceless anisotropy:
 *   δg^{ij}(r) = Δ(r)(3ê_zⁱê_zʲ - δ^{ij})
 *
 * The correction to G(0, Dñ̂) for D >> R_core:
 *   δG = -3(3cos²θ-1)Q / (32π²D³)
 *
 * where θ = angle(D̂, ẑ) and Q = 4π ∫₀^∞ Δ(r) r dr
 *
 * Two-soliton interaction:
 *   U₀(D) = g²/(8πD)           [scalar monopole, from Path 3]
 *   δU(D,θ) = -3g²Q(3cos²θ-1)/(32π²D³)  [tensor quadrupole]
 *
 * Ratio:
 *   δU/U₀ = -3Q(3cos²θ-1)/(4πD²)
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
    return 0;
}

static void free_profile(Profile *p) {
    free(p->r); free(p->f); free(p->fp);
}

/* ========== Physical quantities ========== */

/* Baryon density: B⁰(r) = -f'sin²f/(2π²r²) */
static double baryon_density(double r, double f, double fp, double a) {
    if (r > 1e-10) {
        double sf = sin(f);
        return -(1.0/(2*M_PI*M_PI)) * sf*sf * fp / (r*r);
    } else {
        return a*a*a / (2*M_PI*M_PI);
    }
}

/* Strain: S₀(r) = f'² + 2sin²f/r² */
static double strain_S0(double r, double f, double fp, double a) {
    if (r > 1e-10) {
        double sf = sin(f);
        return fp*fp + 2.0*sf*sf/(r*r);
    } else {
        return 3.0*a*a;
    }
}

/* Skyrme T₀(r) */
static double skyrme_T0(double r, double f, double fp, double a) {
    if (r > 1e-10) {
        double sf = sin(f);
        double sf2 = sf*sf;
        return 2.0*fp*fp*sf2/(r*r) + sf2*sf2/(r*r*r*r);
    } else {
        return 3.0*a*a*a*a;
    }
}

/* ========== BLV metric for hedgehog (isotropic control) ========== */

/*
 * For the hedgehog, P(r) = 2S₀ + 4c₄T₀, m(r) = S₀ + 2c₄T₀
 * P/m = 2 exactly (algebraic identity).
 *
 * For non-hedgehog (hopfion), the BLV metric becomes direction-dependent.
 * Path 4 measured:
 *   Hedgehog:           g^{zz}/g^{xx} = 1.000 (isotropic)
 *   Single hopfion:     g^{zz}/g^{xx} = 1.75  (peak)
 *   Two hopfions:       g^{zz}/g^{xx} = 0.47-1.87
 *   Skyrmion+hopfion:   g^{zz}/g^{xx} = 2.51  (strongest)
 *
 * We model the traceless anisotropy as:
 *   δg^{ij}(r) = Δ(r) × (3ê_zⁱê_zʲ - δ^{ij})
 * where Δ(r) is a profile calibrated to these measurements.
 */

/* ========== Model anisotropy profiles ========== */

/* Gaussian model: Δ(r) = Δ₀ × exp(-r²/σ²) */
static double anisotropy_gaussian(double r, double Delta0, double sigma) {
    return Delta0 * exp(-r*r / (sigma*sigma));
}

/* Soliton-matched model: Δ(r) proportional to Skyrme density T₀ */
static double anisotropy_skyrme(double r, double f, double fp, double a,
                                 double Delta0) {
    double T = skyrme_T0(r, f, fp, a);
    double T0_origin = 3.0*a*a*a*a;
    return Delta0 * T / T0_origin;
}

/* ========== Born integral computation ========== */

/*
 * Full Born correction to the interaction potential between two solitons
 * at separation D.
 *
 * The correction involves two integrals:
 *
 * I₁ = ∫ B⁰(r₁) G₀(r₁,r') d³r₁  [source of soliton 1 weighted by G₀]
 *
 * I₂ = ∫ δg^{ij}(r') ∂'_i∂'_j G₀(r',r₂) d³r'  [metric perturbation × propagator to soliton 2]
 *
 * For point-source approximation (valid when D >> R_core):
 *   I₁ → G₀(0,r') = -1/(4πr')
 *   I₂ → -(3D̂_iD̂_j - δ_{ij})/(4πD³) × 4π ∫ Δ(r) r dr
 *
 * The "metric quadrupole moment":
 *   Q = 4π ∫₀^∞ Δ(r) r dr
 *
 * Non-point-source (finite extent) correction:
 *   Q_full = ∫ B⁰(r₁) × [∫ Δ(r') G₀(r₁,r') F(r',D) d³r'] d³r₁
 * where F captures the ∂∂G₀ structure. This is harder to evaluate
 * but the point-source Q is valid for D > 2R_core.
 */

/* Compute the metric quadrupole moment Q for different models */
static double compute_Q_gaussian(double Delta0, double sigma) {
    /* Q = 4π ∫₀^∞ Δ₀ exp(-r²/σ²) r dr = 4π Δ₀ σ²/2 = 2π Δ₀ σ² */
    return 2.0 * M_PI * Delta0 * sigma * sigma;
}

static double compute_Q_numerical(const Profile *prof, double Delta0,
                                   double a, int use_skyrme_shape) {
    double Q = 0;
    double T0_origin = 3.0*a*a*a*a;
    for (int i = 1; i < prof->n; i++) {
        double r = prof->r[i];
        double Delta;
        if (use_skyrme_shape) {
            double T = skyrme_T0(r, prof->f[i], prof->fp[i], a);
            Delta = Delta0 * T / T0_origin;
        } else {
            Delta = anisotropy_gaussian(r, Delta0, 1.5);
        }
        Q += 4*M_PI * Delta * r * prof->dr;
    }
    return Q;
}

/* ========== Full Born integral (beyond point source) ========== */

/*
 * For finite-size solitons, the Born integral involves:
 *
 *   δU(D) = (g_top²/2) ∫∫ B⁰(r₁) G₀(r₁,r') δg^{ij}(r') ∂'_i∂'_j G₀(r',r₂) B⁰(r₂) d³r₁ d³r₂ d³r'
 *
 * For soliton 2 at Dẑ, point-source approx for soliton 2 but full structure for soliton 1:
 *
 *   δU(D) = (g_top²/2) × Q₂ × ∫ B⁰(r₁) × K(r₁, D) d³r₁
 *
 * where K(r₁, D) = ∫ G₀(r₁,r') δg^{ij}(r') ∂'_i∂'_j G₀(r',Dẑ) d³r'
 *
 * This is expensive (6D integral). Instead, we evaluate the radial form:
 *
 * For s-wave source B⁰(r) and axial anisotropy, the s-wave projection of K gives:
 *
 *   K_s(r₁, D) = (1/2) ∫₋₁¹ K(r₁ê_z, D) d(cosα) × (4π)  [angle average]
 *
 * Using the multipole expansion...
 *
 * Actually, for an s-wave source weighted by a quadrupole metric perturbation,
 * the leading correction is exactly the point-source result (by symmetry):
 *
 *   The s-wave component of the source ∫B⁰ d³r = Q_baryon = 1 (total charge)
 *   The quadrupole component of the source is zero (B⁰ is spherically symmetric)
 *
 * So the finite-size correction to the quadrupole moment is ZERO for a hedgehog
 * source. The point-source approximation is EXACT for the tensor correction.
 *
 * This is because: a spherically symmetric source cannot couple to a quadrupole
 * perturbation at leading order. The correction appears at order (R_core/D)⁴.
 */

/* ========== Main ========== */

int main(int argc, char *argv[])
{
    const char *profile_file = "data/profiles/profile_sigma_e1.dat";
    double e_param = 1.0;
    double rho0 = 1.0;

    /* Model parameters (calibrated to Path 4 results) */
    double Delta0 = 1.0;    /* traceless anisotropy amplitude (from Path 4) */
    double sigma = 1.5;     /* anisotropy range (code lengths) */
    int use_skyrme = 0;     /* 0: Gaussian, 1: Skyrme-shaped anisotropy */

    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "-profile") && i+1 < argc) profile_file = argv[++i];
        else if (!strcmp(argv[i], "-e") && i+1 < argc) e_param = atof(argv[++i]);
        else if (!strcmp(argv[i], "-Delta0") && i+1 < argc) Delta0 = atof(argv[++i]);
        else if (!strcmp(argv[i], "-sigma") && i+1 < argc) sigma = atof(argv[++i]);
        else if (!strcmp(argv[i], "-skyrme")) use_skyrme = 1;
        else {
            fprintf(stderr, "Usage: %s [-profile P] [-Delta0 D] [-sigma S] [-skyrme]\n",
                    argv[0]);
            return 1;
        }
    }

    double c4 = 2.0 * rho0*rho0 / (e_param*e_param);

    /* Read profile */
    Profile prof;
    if (read_profile(profile_file, &prof) < 0) return 1;
    int n = prof.n;
    double a = -prof.fp[0];

    printf("================================================================\n");
    printf("  Scalar Propagator on Anisotropic BLV Metric\n");
    printf("  Path 3 + Path 4 Combination Test\n");
    printf("================================================================\n\n");

    printf("Profile: %s (%d points)\n", profile_file, n);
    printf("Soliton: a=%.4f, c4=%.4f, e=%.4f\n\n", a, c4, e_param);

    /* ===== Section 1: Soliton properties ===== */
    printf("--- Soliton Properties ---\n\n");

    double E2 = 0, E4 = 0, Q_baryon = 0;
    for (int i = 1; i < n; i++) {
        double r = prof.r[i], r0 = prof.r[i-1];
        double f = prof.f[i], fp = prof.fp[i];
        double f0 = prof.f[i-1], fp0 = prof.fp[i-1];

        double S = strain_S0(r, f, fp, a);
        double S0 = strain_S0(r0, f0, fp0, a);
        double T = skyrme_T0(r, f, fp, a);
        double T0v = skyrme_T0(r0, f0, fp0, a);
        double B = baryon_density(r, f, fp, a);
        double B0 = baryon_density(r0, f0, fp0, a);

        E2 += 0.5*prof.dr * (2*M_PI*rho0*rho0*(S0*r0*r0 + S*r*r));
        E4 += 0.5*prof.dr * (4*M_PI*rho0*rho0*rho0*rho0/(e_param*e_param)
                             *(T0v*r0*r0 + T*r*r));
        Q_baryon += 0.5*prof.dr * (4*M_PI*(B0*r0*r0 + B*r*r));
    }

    double E_sol = E2 + E4;
    printf("E2 = %.4f, E4 = %.4f, E_sol = %.4f\n", E2, E4, E_sol);
    printf("Virial E2/E4 = %.6f, Q = %.6f\n\n", E2/E4, Q_baryon);

    /* ===== Section 2: BLV metric (hedgehog = isotropic control) ===== */
    printf("--- BLV Metric: Hedgehog (Isotropic Control) ---\n\n");

    printf("For hedgehog: g^{ij}_eff = P(r)/m(r) * delta^{ij}\n");
    printf("P/m = 2 exactly (algebraic identity for L2+L4)\n\n");

    printf("  %8s  %12s  %12s  %8s\n", "r", "P(r)", "m(r)", "P/m");
    for (int i = 0; i < n; i += n/15) {
        double r = prof.r[i];
        double f = prof.f[i], fp = prof.fp[i];
        double S = strain_S0(r, f, fp, a);
        double T = skyrme_T0(r, f, fp, a);

        double P = 2*S + 4*c4*T;
        double m = S + 2*c4*T;

        printf("  %8.4f  %12.6f  %12.6f  %8.4f\n", r, P, m, (m > 1e-30) ? P/m : 2.0);
    }

    /* ===== Section 3: Model anisotropy and quadrupole moment ===== */
    printf("\n--- Metric Anisotropy (from Path 4 data) ---\n\n");

    printf("Path 4 measurements:\n");
    printf("  Hedgehog:          g^{zz}/g^{xx} = 1.000 (isotropic)\n");
    printf("  Single hopfion:    g^{zz}/g^{xx} = 1.75  (75%% peak)\n");
    printf("  Skyrmion+hopfion:  g^{zz}/g^{xx} = 2.51  (151%% peak)\n\n");

    printf("Model: delta_g^{ij}(r) = Delta(r) * (3 z_i z_j - delta_ij)\n");
    if (use_skyrme) {
        printf("  Profile: Skyrme-shaped, Delta(r) = Delta0 * T0(r)/T0(0)\n");
    } else {
        printf("  Profile: Gaussian, Delta(r) = Delta0 * exp(-r^2/sigma^2)\n");
    }
    printf("  Delta0 = %.3f (traceless anisotropy amplitude)\n", Delta0);
    printf("  sigma  = %.3f code lengths (anisotropy range)\n\n", sigma);

    /* Calibration from Path 4:
     * Skyrmion+hopfion: g^{zz}/g^{xx} = 2.51
     * On hedgehog background, g^{ij} = 2δ^{ij}. So g^{zz} = 2 + 2Δ, g^{xx} = 2 - Δ
     * Ratio: (2+2Δ)/(2-Δ) = 2.51 → 2+2Δ = 2.51(2-Δ) = 5.02 - 2.51Δ
     *       → 4.51Δ = 3.02 → Δ = 0.67
     *
     * Single hopfion: (2+2Δ)/(2-Δ) = 1.75 → 2+2Δ = 3.5-1.75Δ → 3.75Δ = 1.5 → Δ = 0.40
     */
    double Delta_skyrmion_hopfion = 3.02 / 4.51;  /* 0.670 */
    double Delta_single_hopfion = 1.5 / 3.75;      /* 0.400 */

    printf("Calibrated anisotropy amplitudes:\n");
    printf("  Single hopfion:   Delta0 = %.3f\n", Delta_single_hopfion);
    printf("  Skyrmion+hopfion: Delta0 = %.3f\n", Delta_skyrmion_hopfion);
    printf("  Using: Delta0 = %.3f\n\n", Delta0);

    /* Compute quadrupole moment */
    double Q_gauss = compute_Q_gaussian(Delta0, sigma);
    double Q_num = compute_Q_numerical(&prof, Delta0, a, use_skyrme);

    printf("Metric quadrupole moment Q = 4pi int Delta(r) r dr:\n");
    if (!use_skyrme) {
        printf("  Analytical (Gaussian): Q = 2*pi*Delta0*sigma^2 = %.4f\n", Q_gauss);
    }
    printf("  Numerical:             Q = %.4f\n\n", Q_num);

    double Q = use_skyrme ? Q_num : Q_gauss;

    /* ===== Section 4: Born approximation results ===== */
    printf("--- Born Approximation: Tensor Correction ---\n\n");

    printf("Scalar monopole:      U0(D) = g^2 / (8*pi*D)\n");
    printf("Tensor quadrupole:    dU(D,theta) = -3*g^2*Q*(3cos^2(theta)-1) / (32*pi^2*D^3)\n");
    printf("Ratio:                dU/U0 = -3*Q*(3cos^2(theta)-1) / (4*pi*D^2)\n\n");

    printf("Tensor correction vs distance (theta=0, along anisotropy axis):\n");
    printf("  %10s  %12s  %12s  %12s  %12s\n",
           "D (code)", "D (fm)", "|dU/U0|", "dU/U0 max", "dU/U0 avg");

    double distances[] = {2.0, 3.0, 5.0, 7.0, 10.0, 15.0, 20.0, 50.0, 100.0,
                          1000.0, 1e6, 1e10, 1e20};
    int n_dist = sizeof(distances) / sizeof(distances[0]);

    for (int id = 0; id < n_dist; id++) {
        double D = distances[id];

        /* theta=0: 3cos²θ-1 = 2 (maximum correction, along z-axis) */
        double ratio_max = 3.0 * Q * 2.0 / (4*M_PI * D*D);

        /* theta=π/2: 3cos²θ-1 = -1 (perpendicular) */
        double ratio_perp = 3.0 * Q * 1.0 / (4*M_PI * D*D);

        /* angle-averaged |3cos²θ-1|: <|3x²-1|> = 4/3 for x=cosθ ∈ [-1,1] */
        double ratio_avg = 3.0 * Q * (4.0/3.0) / (4*M_PI * D*D);

        printf("  %10.1e  %12.4e  %12.4e  %12.4e  %12.4e\n",
               D, D*0.5624, ratio_max, ratio_max, ratio_avg);
    }

    /* ===== Section 5: Comparison with GR ===== */
    printf("\n--- Comparison with General Relativity ---\n\n");

    printf("In GR, gravity is mediated by a SPIN-2 graviton.\n");
    printf("This means the tensor structure is in the MONOPOLE (1/r) itself:\n");
    printf("  - Newtonian potential: Phi = -GM/r         [monopole, spin-0 component]\n");
    printf("  - Light deflection:    theta = 4GM/(c^2 b) [DOUBLE the Newtonian value]\n");
    printf("  - The factor of 2 comes from spin-2 (tensor) nature of the mediator\n");
    printf("  - This doubling occurs at ALL distances, not just near the source\n\n");

    printf("In our model (Path 3 scalar + Path 4 tensor):\n");
    printf("  - Monopole:   U0 ~ 1/D          [scalar, spin-0]\n");
    printf("  - Quadrupole: dU ~ 1/D^3         [tensor, from core anisotropy]\n");
    printf("  - Light deflection: theta_scalar = 2GM/(c^2 b)  [HALF of GR]\n");
    printf("  - Tensor correction to deflection: ~ (R_core/b)^2 << 1\n\n");

    printf("The fundamental problem:\n");
    printf("  GR spin-2:  tensor structure is in the PROPAGATOR (graviton carries spin-2)\n");
    printf("  Our model:  tensor structure is in the BACKGROUND (localized, core-scale)\n");
    printf("  A scalar mediator cannot acquire spin-2 from a localized background.\n");
    printf("  The tensor correction dies as (R_core/D)^2 = (%.1f fm / D)^2\n\n",
           sigma * 0.5624);

    /* Quantitative comparison */
    printf("At what distance does the tensor correction become negligible?\n");
    double D_1pct = sqrt(3.0 * Q * 2.0 / (4*M_PI * 0.01));
    double D_01pct = sqrt(3.0 * Q * 2.0 / (4*M_PI * 0.001));
    printf("  |dU/U0| < 1%%:   D > %.1f code = %.1f fm\n", D_1pct, D_1pct*0.5624);
    printf("  |dU/U0| < 0.1%%: D > %.1f code = %.1f fm\n", D_01pct, D_01pct*0.5624);

    /* ===== Section 6: Can we do better? ===== */
    printf("\n--- Can Any Modification Fix This? ---\n\n");

    printf("For the tensor correction to be 1/D (not 1/D^3), we would need:\n\n");

    printf("1. LONG-RANGE metric anisotropy: Delta(r) ~ 1/r^2 (not localized)\n");
    printf("   Then Q diverges, and the Born integral changes character.\n");
    printf("   But: BLV metric anisotropy requires |nabla q|^2 ~ 1/r^4 at large r,\n");
    printf("   while the Skyrmion has f ~ e^{-mu r}/r (exponential). IMPOSSIBLE.\n\n");

    printf("2. TENSOR mediator: replace scalar p with a spin-2 field h_{ij}\n");
    printf("   This requires a massless rank-2 field in the Lagrangian.\n");
    printf("   Cl+(3,0,1) has: scalars (s,p) and vectors (f_i, j_i),\n");
    printf("   but NO symmetric rank-2 tensor field. Would need Cl+(3,0,1) x Cl+(3,0,1)\n");
    printf("   or similar extension.\n\n");

    printf("3. EMERGENT spin-2 from scalar-scalar coupling:\n");
    printf("   Two scalar fields with specific coupling can mimic spin-2.\n");
    printf("   (Massive gravity bi-metric theories work this way.)\n");
    printf("   Our p and s fields could in principle do this, but both are\n");
    printf("   constrained: s = sqrt(rho0^2 - |f|^2), p is pseudoscalar.\n");
    printf("   The constraint removes the necessary degrees of freedom.\n\n");

    printf("4. NON-LINEAR backreaction (self-consistent loop):\n");
    printf("   p modifies q via constraint -> modified BLV metric -> modified p -> ...\n");
    printf("   But Path 6 showed: constraint eigenmode is LOCALIZED (0.46 fm).\n");
    printf("   The sourced 1/r component (Path 3) doesn't modify q.\n");
    printf("   No self-consistent long-range backreaction.\n\n");

    /* ===== Section 7: Skyrmion+hopfion specific calculation ===== */
    printf("--- Skyrmion+Hopfion Specific Numbers ---\n\n");

    double configs[][2] = {
        {Delta_single_hopfion, 1.5},    /* single hopfion */
        {Delta_skyrmion_hopfion, 1.5},  /* skyrmion+hopfion (strongest) */
    };
    const char *config_names[] = {
        "Single hopfion (Delta=0.40)",
        "Skyrmion+hopfion (Delta=0.67, strongest)",
    };

    for (int ic = 0; ic < 2; ic++) {
        double D0 = configs[ic][0];
        double sig = configs[ic][1];
        double Qc = 2.0 * M_PI * D0 * sig * sig;

        printf("%s:\n", config_names[ic]);
        printf("  Q = %.4f\n", Qc);
        printf("  Tensor correction (theta=0):\n");

        double test_D[] = {3.0, 5.0, 10.0, 20.0, 100.0};
        for (int j = 0; j < 5; j++) {
            double D = test_D[j];
            double ratio = 3.0 * Qc * 2.0 / (4*M_PI * D*D);
            printf("    D=%5.1f (%5.1f fm): dU/U0 = %.4e", D, D*0.5624, ratio);
            if (ratio > 0.01) printf(" (%.1f%%)", ratio*100);
            printf("\n");
        }
        printf("\n");
    }

    /* ===== Section 8: What about nuclear-scale physics? ===== */
    printf("--- Nuclear-Scale Relevance ---\n\n");

    printf("Although the tensor correction is negligible for gravity (macroscopic D),\n");
    printf("it IS significant at nuclear distances (D ~ 2-5 fm = 3.5-9 code):\n\n");

    printf("At D = 5 code (2.8 fm), skyrmion+hopfion:\n");
    double D_nuc = 5.0;
    double Q_sky_hop = 2.0 * M_PI * Delta_skyrmion_hopfion * 1.5*1.5;
    double ratio_nuc = 3.0 * Q_sky_hop * 2.0 / (4*M_PI * D_nuc*D_nuc);
    printf("  dU/U0 = %.1f%% (significant tensor nuclear force!)\n", ratio_nuc*100);
    printf("  This is a QUADRUPOLE (P2) correction to the nucleon-nucleon potential.\n");
    printf("  Sign: attractive along z, repulsive in equatorial plane.\n");
    printf("  Angular dependence: (3cos^2(theta) - 1)\n\n");

    printf("This tensor correction at nuclear distances is physically real.\n");
    printf("It represents the ORIENTATION-DEPENDENT force between non-spherical solitons.\n");
    printf("In nuclear physics, this is analogous to the tensor force from pion exchange,\n");
    printf("which is known to be crucial for deuteron binding and nuclear structure.\n");

    /* ===== Final summary ===== */
    printf("\n================================================================\n");
    printf("  CONCLUSION\n");
    printf("================================================================\n\n");

    printf("Combining Path 3 (scalar 1/r) with Path 4 (tensor metric) gives:\n\n");

    printf("  U(D,theta) = g^2/(8*pi*D) * [1 - %.2f*(3cos^2(theta)-1)/D^2 + O(1/D^4)]\n\n",
           3.0*Q_sky_hop*2.0/(4*M_PI));

    printf("1. MONOPOLE (1/D): scalar, spin-0, from B0*p coupling (Path 3)\n");
    printf("   -> Correct 1/r potential, correct sign, but HALF light deflection\n\n");

    printf("2. QUADRUPOLE (1/D^3): tensor, from BLV metric anisotropy (Path 4)\n");
    printf("   -> Correct tensor structure, but DIES as (R_core/D)^2\n");
    printf("   -> Nuclear-scale effect (%.1f%% at 2.8 fm), not gravity\n\n", ratio_nuc*100);

    printf("3. CANNOT reproduce GR: spin-2 gravity requires tensor structure\n");
    printf("   IN the 1/D monopole, not as a 1/D^3 correction.\n");
    printf("   A scalar mediator on a localized anisotropic background\n");
    printf("   cannot produce long-range spin-2 effects.\n\n");

    printf("The Cl+(3,0,1) Skyrmion framework provides:\n");
    printf("  [x] Topological solitons (nucleons)\n");
    printf("  [x] Pion cloud (Path 6 eigenvalue, 0.46 fm)\n");
    printf("  [x] 1/r monopole potential (Path 3, scalar)\n");
    printf("  [x] Nuclear tensor force (Path 4 + this, quadrupole)\n");
    printf("  [ ] Spin-2 gravity (requires tensor mediator, absent in algebra)\n");

    free_profile(&prof);
    return 0;
}
