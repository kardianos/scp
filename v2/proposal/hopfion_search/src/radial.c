/*
 * radial.c — Hedgehog Skyrmion Profile via Shooting Method
 *
 * Solves the Euler-Lagrange ODE for the hedgehog profile f(r):
 *
 *   f''(r² + 2c₄ sin²f + c₆ sin⁴f/r²) + 2f'r + c₄ f'² sin(2f)
 *     + c₆ f'² sin²f sin(2f)/r² - 2c₆ f' sin⁴f/r³
 *     - sin(2f)(1 + c₄ sin²f/r²) - m_π² r² sinf = 0
 *
 * where c₄ = 2ρ₀²/e², c₆ = λ₆/(2π⁴ρ₀²), m_π = pion mass in 1/length units.
 * Boundary conditions: f(0) = π, f(∞) = 0.
 *
 * Near r=0: f(r) = π - a·r + O(r³), so we shoot from r=ε with f'(0)=-a
 * and adjust a until f(R_max) ≈ 0.
 *
 * Energy in hedgehog ansatz (sigma model |q| = ρ₀):
 *   E₂ = 2π ρ₀² ∫₀^∞ [f'²r² + 2sin²f] dr
 *   E₄ = (4π ρ₀⁴/e²) ∫₀^∞ [2f'²sin²f + sin⁴f/r²] dr
 *   E₆ = (λ₆/π³) ∫₀^∞ f'² sin⁴f/r² dr
 *   E_V = 4π ρ₀² m_π² ∫₀^∞ r²(1-cosf) dr
 *
 * The L₆ = λ₆(B⁰)² term where B⁰ = -f'sin²f/(2π²r²) is the baryon density.
 *
 * Derrick scaling: E₂ ~ λ, E₄ ~ λ⁻¹, E₆ ~ λ⁻³, E_V ~ λ³
 * Virial: E₂ - E₄ - 3E₆ + 3E_V = 0
 * Asymptotic: f ~ exp(-m_π r)/r for m_π > 0 (vs 1/r² for m_π = 0)
 */

#define _USE_MATH_DEFINES
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ---------- ODE right-hand side ---------- */

/* The ODE: f'' = F(r, f, f')
 *
 * From EL equation of E = E2 + E4 + E6:
 *
 * The L₂+L₄ part gives (after dividing by 4πρ₀²):
 *   f''(r² + 2c4*sin²f) + 2f'r + c4*f'²sin(2f)
 *   - sin(2f) - c4*sin(2f)*sin²f/r² = 0
 *
 * The L₆ EL terms, also divided by 4πρ₀² (with c6 = λ₆/(2π⁴ρ₀²)):
 *   f'' coefficient (denom): += c6*sin⁴f/r²
 *   numerator:               -= c6*f'²*sin²f*sin(2f)/r²  +  2*c6*f'*sin⁴f/r³
 *
 * Full ODE:
 *   f''(r² + 2c4*sin²f + c6*sin⁴f/r²) + 2f'r + c4*f'²sin(2f)
 *   - c6*f'²*sin²f*sin(2f)/r² + 2*c6*f'*sin⁴f/r³
 *   - sin(2f)(1 + c4*sin²f/r²) - m_π²*r²*sinf = 0
 *
 * Near r=0: sinf ~ ar, so c6*sin⁴f/r² ~ c6*a⁴*r² → 0. Leading behavior unchanged.
 */
static double f_rhs(double r, double f, double fp, double c4, double c6, double m_pi2)
{
    double sf = sin(f);
    double s2f = sin(2.0 * f);
    double sf2 = sf * sf;
    double sf4 = sf2 * sf2;
    double r2 = r * r;
    double denom = r2 + 2.0 * c4 * sf2;

    if (r > 1e-14)
        denom += c6 * sf4 / r2;

    if (denom < 1e-30) {
        return 0.0;
    }

    double num = s2f * (1.0 + c4 * sf2 / r2)
               + m_pi2 * r2 * sf
               - 2.0 * fp * r
               - c4 * fp * fp * s2f;

    /* L₆ contributions: EL₆/(4πρ₀²) adds -c₆f'²sin²fsin(2f)/r² + 2c₆f'sin⁴f/r³ */
    if (c6 > 0 && r > 1e-14) {
        num -= c6 * fp * fp * sf2 * s2f / r2;
        num += 2.0 * c6 * fp * sf4 / (r2 * r);
    }

    return num / denom;
}

/* ---------- RK4 integration ---------- */

typedef struct {
    int N;
    double dr;
    double *r;
    double *f;
    double *fp;
} Profile;

static Profile *profile_alloc(int N, double R_max)
{
    Profile *p = (Profile *)malloc(sizeof(Profile));
    p->N = N;
    p->dr = R_max / N;
    p->r = (double *)calloc((size_t)(N + 1), sizeof(double));
    p->f = (double *)calloc((size_t)(N + 1), sizeof(double));
    p->fp = (double *)calloc((size_t)(N + 1), sizeof(double));
    for (int i = 0; i <= N; i++)
        p->r[i] = i * p->dr;
    return p;
}

static void profile_free(Profile *p) {
    if (p) { free(p->r); free(p->f); free(p->fp); free(p); }
}

/* Integrate f'' = F(r,f,f') from r=eps using RK4 */
static double shoot(Profile *p, double a, double c4, double c6, double m_pi2)
{
    double dr = p->dr;

    /* Initial conditions at r=eps (first grid point) */
    p->r[0] = 0.0;
    p->f[0] = M_PI;
    p->fp[0] = -a;

    /* At r=0, f''(0) = 0 by symmetry (m_π² term is O(r³), doesn't change this).
     * c₆ terms are also O(r²) at the origin. Start with standard linear step. */
    p->r[1] = dr;
    p->f[1] = M_PI - a * dr;
    p->fp[1] = -a;

    for (int i = 1; i < p->N; i++) {
        double r0 = p->r[i];
        double f0 = p->f[i];
        double v0 = p->fp[i];
        double h = dr;

        /* RK4 for the system f' = v, v' = F(r, f, v) */
        double k1f = v0;
        double k1v = f_rhs(r0, f0, v0, c4, c6, m_pi2);

        double k2f = v0 + 0.5 * h * k1v;
        double k2v = f_rhs(r0 + 0.5*h, f0 + 0.5*h*k1f, v0 + 0.5*h*k1v, c4, c6, m_pi2);

        double k3f = v0 + 0.5 * h * k2v;
        double k3v = f_rhs(r0 + 0.5*h, f0 + 0.5*h*k2f, v0 + 0.5*h*k2v, c4, c6, m_pi2);

        double k4f = v0 + h * k3v;
        double k4v = f_rhs(r0 + h, f0 + h*k3f, v0 + h*k3v, c4, c6, m_pi2);

        p->f[i+1] = f0 + (h / 6.0) * (k1f + 2*k2f + 2*k3f + k4f);
        p->fp[i+1] = v0 + (h / 6.0) * (k1v + 2*k2v + 2*k3v + k4v);
        p->r[i+1] = r0 + h;
    }

    return p->f[p->N];  /* f(R_max) — should be 0 */
}

/* ---------- Energy computation ---------- */

typedef struct { double E2, E4, E6, EV, Etotal, Q; } Energy;

static Energy compute_energy(const Profile *p, double e, double rho0,
                              double lambda6, double m_pi2)
{
    double inv_e2 = 1.0 / (e * e);
    double rho2 = rho0 * rho0;
    double rho4 = rho2 * rho2;
    double dr = p->dr;
    double pi3 = M_PI * M_PI * M_PI;

    double e2_sum = 0.0, e4_sum = 0.0, e6_sum = 0.0, ev_sum = 0.0, q_sum = 0.0;

    for (int i = 0; i <= p->N; i++) {
        double r = p->r[i];
        double fi = p->f[i];
        double fpi = p->fp[i];
        double sf = sin(fi);
        double sf2 = sf * sf;

        double w = dr;
        if (i == 0 || i == p->N) w *= 0.5;

        e2_sum += (fpi * fpi * r * r + 2.0 * sf2) * w;
        if (r > 1e-14) {
            e4_sum += (2.0 * fpi * fpi * sf2 + sf2 * sf2 / (r * r)) * w;
            e6_sum += fpi * fpi * sf2 * sf2 / (r * r) * w;
        }
        ev_sum += r * r * (1.0 - cos(fi)) * w;
        q_sum += (-fpi * sf2) * w;
    }

    Energy en;
    en.E2 = 2.0 * M_PI * rho2 * e2_sum;
    en.E4 = 4.0 * M_PI * inv_e2 * rho4 * e4_sum;
    en.E6 = (lambda6 / pi3) * e6_sum;
    en.EV = 4.0 * M_PI * rho2 * m_pi2 * ev_sum;
    en.Etotal = en.E2 + en.E4 + en.E6 + en.EV;
    en.Q = (2.0 / M_PI) * q_sum;
    return en;
}

/* ---------- Main ---------- */

int main(int argc, char **argv)
{
    int Nr = 10000;
    double R_max = 30.0;
    double e_skyrme = 4.0;
    double rho0 = 1.0;
    double m_pi = 0.0;      /* pion mass in 1/(code length) units */
    double lambda6 = 0.0;   /* sextic (BPS) coupling */
    const char *outfile = "profile.dat";

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-Nr") == 0 && i+1 < argc) Nr = atoi(argv[++i]);
        else if (strcmp(argv[i], "-Rmax") == 0 && i+1 < argc) R_max = atof(argv[++i]);
        else if (strcmp(argv[i], "-e") == 0 && i+1 < argc) e_skyrme = atof(argv[++i]);
        else if (strcmp(argv[i], "-rho0") == 0 && i+1 < argc) rho0 = atof(argv[++i]);
        else if (strcmp(argv[i], "-mpi") == 0 && i+1 < argc) m_pi = atof(argv[++i]);
        else if (strcmp(argv[i], "-lam6") == 0 && i+1 < argc) lambda6 = atof(argv[++i]);
        else if (strcmp(argv[i], "-o") == 0 && i+1 < argc) outfile = argv[++i];
        else {
            fprintf(stderr, "Usage: %s [-Nr N] [-Rmax R] [-e E] [-rho0 R] [-mpi M] [-lam6 L] [-o outfile]\n", argv[0]);
            return 1;
        }
    }

    /* c4 = ratio of E4 prefactor to E2 prefactor in hedgehog energy:
     * E2 prefactor = 2*pi*rho0^2, E4 prefactor = 4*pi*rho0^4/e^2
     * ratio = 2*rho0^2/e^2 */
    double c4 = 2.0 * rho0 * rho0 / (e_skyrme * e_skyrme);
    /* c6 for ODE: L₆ EL divided by 4πρ₀² normalization of L₂+L₄ terms */
    double pi4 = M_PI * M_PI * M_PI * M_PI;
    double c6 = lambda6 / (2.0 * pi4 * rho0 * rho0);
    double m_pi2 = m_pi * m_pi;

    printf("=== Hedgehog Skyrmion — Shooting Method ===\n");
    printf("Grid: Nr=%d, dr=%.6f, R_max=%.1f\n", Nr, R_max / Nr, R_max);
    printf("Parameters: e=%.4f (c4=%.6f), rho0=%.4f\n", e_skyrme, c4, rho0);
    if (lambda6 > 0)
        printf("Sextic: lambda6=%.6f (c6=%.6f)\n", lambda6, c6);
    if (m_pi > 0)
        printf("Pion mass: m_pi=%.6f (1/length), m_pi²=%.6f, decay length=%.3f\n",
               m_pi, m_pi2, 1.0/m_pi);
    printf("\n");

    Profile *p = profile_alloc(Nr, R_max);

    /* Bisection on the shooting parameter a = -f'(0)
     * a scales as 1/sqrt(c4) ~ e/sqrt(2), so adapt the search range */
    double a_lo = 0.1, a_hi = 20.0 / sqrt(c4);

    /* Find bracket: f(R_max) should change sign as a varies */
    double f_lo = shoot(p, a_lo, c4, c6, m_pi2);
    double f_hi = shoot(p, a_hi, c4, c6, m_pi2);
    printf("Bracket: a_lo=%.4f -> f(Rmax)=%.6f, a_hi=%.4f -> f(Rmax)=%.6f\n",
           a_lo, f_lo, a_hi, f_hi);

    if (f_lo * f_hi > 0) {
        /* Try to find bracket by scanning */
        printf("Scanning for bracket...\n");
        int found = 0;
        double step = a_hi / 500.0;
        for (double a = step; a <= a_hi; a += step) {
            double fend = shoot(p, a, c4, c6, m_pi2);
            if (fend * f_lo < 0) {
                a_hi = a;
                f_hi = fend;
                a_lo = a - step;
                f_lo = shoot(p, a_lo, c4, c6, m_pi2);
                found = 1;
                printf("Found bracket: a in [%.4f, %.4f]\n", a_lo, a_hi);
                break;
            }
        }
        if (!found) {
            fprintf(stderr, "Could not find shooting bracket\n");
            profile_free(p);
            return 1;
        }
    }

    /* Bisection to find a */
    printf("\nBisection search for a = -f'(0):\n");
    double a_best = 0;
    for (int iter = 0; iter < 100; iter++) {
        double a_mid = 0.5 * (a_lo + a_hi);
        double f_mid = shoot(p, a_mid, c4, c6, m_pi2);

        if (iter < 10 || iter % 10 == 0)
            printf("  iter %3d: a=%.12f, f(Rmax)=%+.6e\n", iter, a_mid, f_mid);

        if (fabs(f_mid) < 1e-12 || (a_hi - a_lo) < 1e-14) {
            a_best = a_mid;
            printf("  Converged at iter %d: a=%.14f, f(Rmax)=%+.6e\n", iter, a_mid, f_mid);
            break;
        }

        if (f_lo * f_mid < 0) {
            a_hi = a_mid;
            f_hi = f_mid;
        } else {
            a_lo = a_mid;
            f_lo = f_mid;
        }
        a_best = a_mid;
    }

    /* Final integration with best a */
    shoot(p, a_best, c4, c6, m_pi2);

    /* Compute energy */
    Energy en = compute_energy(p, e_skyrme, rho0, lambda6, m_pi2);

    printf("\n=== Result ===\n");
    printf("a = -f'(0) = %.14f\n", a_best);
    printf("E_total = %.12e\n", en.Etotal);
    printf("E2      = %.12e\n", en.E2);
    printf("E4      = %.12e\n", en.E4);
    if (lambda6 > 0)
        printf("E6      = %.12e\n", en.E6);
    printf("E_V     = %.12e\n", en.EV);
    printf("Q       = %.8f\n", en.Q);
    /* Derrick virial: E₂ - E₄ - 3E₆ + 3E_V = 0
     * (E₂~λ, E₄~λ⁻¹, E₆~λ⁻³, E_V~λ³) */
    double virial = en.E2 - en.E4 - 3.0 * en.E6 + 3.0 * en.EV;
    if (m_pi > 0) {
        printf("E2/E4   = %.6f\n", en.E2 / en.E4);
        printf("Virial: E2 - E4 - 3*E6 + 3*EV = %.6e (should be ~0)\n", virial);
    } else {
        printf("E2/(E4+3*E6) = %.6f (should be ~1 at virial)\n",
               en.E2 / (en.E4 + 3.0 * en.E6));
        printf("Virial: E2 - E4 - 3*E6 = %.6e (should be ~0)\n", virial);
    }

    /* The FB bound in our normalization: E_FB = 6*sqrt(2)*pi^2*rho0^3/e
     * Derived by rescaling r = R*sqrt(c4) to map to standard Skyrme ODE,
     * where E_std_FB = 12*pi^2, giving E_FB = sqrt(2)/(2e) * 12*pi^2 * rho0^3 */
    double E_FB = 6.0 * sqrt(2.0) * M_PI * M_PI * rho0 * rho0 * rho0 / e_skyrme;
    printf("\nFaddeev-Bogomolny bound (6*sqrt(2)*pi^2*rho0^3/e): %.6e\n", E_FB);
    printf("Ratio E/E_FB = %.6f (standard: ~1.232 for m_pi=0)\n", en.Etotal / E_FB);

    /* CHPT mass formula: E = E₂+E₄+E₆+E_V, virial: E₂ = E₄+3E₆-3E_V
     * → E = 2E₄ + 4E₆ - 2E_V */
    printf("\nCHPT mass formula: Mc^2 = 2*E4 + 4*E6 - 2*EV = %.12e\n",
           2.0*en.E4 + 4.0*en.E6 - 2.0*en.EV);
    if (m_pi > 0)
        printf("  (Virial: E2 = E4+3*E6-3*EV, so Mc^2 = E_total)\n");
    else
        printf("  (Virial: E2 = E4+3*E6, so Mc^2 = 2*E4+4*E6 = E_total)\n");

    /* Profile output */
    FILE *fp = fopen(outfile, "w");
    if (fp) {
        fprintf(fp, "# Hedgehog Skyrmion profile: e=%.4f rho0=%.4f m_pi=%.6f c4=%.6f\n",
                e_skyrme, rho0, m_pi, c4);
        fprintf(fp, "# r  f(r)  f'(r)  baryon_density  energy_density_E2  energy_density_E4\n");
        for (int i = 0; i <= p->N; i++) {
            double r = p->r[i];
            double fi = p->f[i];
            double fpi = p->fp[i];
            double sf = sin(fi);
            double sf2 = sf * sf;
            double bdens = (r > 1e-14) ? -fpi * sf2 / (M_PI * r * r) : 0.0;
            double ed2 = fpi * fpi * r * r + 2.0 * sf2;
            double ed4 = (r > 1e-14) ?
                (2.0 * fpi * fpi * sf2 + sf2 * sf2 / (r * r)) / (e_skyrme * e_skyrme) : 0.0;
            fprintf(fp, "%.8e  %.12e  %.12e  %.12e  %.12e  %.12e\n",
                    r, fi, fpi, bdens, ed2, ed4);
        }
        fclose(fp);
        printf("\nProfile written to %s\n", outfile);
    }

    profile_free(p);
    return 0;
}
