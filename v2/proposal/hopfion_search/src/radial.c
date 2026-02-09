/*
 * radial.c — Hedgehog Skyrmion Profile via Shooting Method
 *
 * Solves the Euler-Lagrange ODE for the hedgehog profile f(r):
 *
 *   f''(r² + c₄ sin²f) + 2f'r + c₄ f'² sin(2f) - sin(2f) - c₄ sin²(2f)/(2r²) = 0
 *
 * where c₄ = 2ρ₀²/e². Boundary conditions: f(0) = π, f(∞) = 0.
 *
 * Near r=0: f(r) = π - a·r + O(r³), so we shoot from r=ε with f'(0)=-a
 * and adjust a until f(R_max) ≈ 0.
 *
 * Energy in hedgehog ansatz (sigma model |q| = ρ₀):
 *   E₂ = 2π ρ₀² ∫₀^∞ [f'²r² + 2sin²f] dr
 *   E₄ = (4π ρ₀⁴/e²) ∫₀^∞ [2f'²sin²f + sin⁴f/r²] dr
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
 * From EL equation of E = E2 + E4:
 *
 * The Lagrangian density (integrand of E/(2π)):
 *   L(r,f,f') = f'²r² + 2sin²f + c4*(2f'²sin²f + sin⁴f/r²)
 *
 * EL: d/dr(∂L/∂f') = ∂L/∂f
 *   ∂L/∂f' = 2f'r² + c4*4f'sin²f = 2f'(r² + c4*2sin²f)
 *   d/dr(∂L/∂f') = 2f''(r² + c4*2sin²f) + 2f'(2r + c4*4f'sin(f)cos(f))
 *                 = 2f''(r² + c4*2sin²f) + 4f'r + c4*8f'²sin(f)cos(f)
 *   ∂L/∂f = 2sin(2f) + c4*(4f'²sin(f)cos(f) + 4sin³(f)cos(f)/r²)
 *          = 2sin(2f) + c4*4sin(f)cos(f)*(f'² + sin²(f)/r²)
 *
 * Setting d/dr(∂L/∂f') - ∂L/∂f = 0:
 *   2f''(r² + 2c4*sin²f) + 4f'r + 8c4*f'²sinf*cosf
 *   - 2sin(2f) - 4c4*sinf*cosf*(f'² + sin²f/r²) = 0
 *
 * Simplify (using sin(2f) = 2sinf*cosf):
 *   f'² terms: 8c4*sinf*cosf - 4c4*sinf*cosf = 4c4*sinf*cosf = 2c4*f'²sin(2f)
 *   2f''(r² + 2c4*sin²f) + 4f'r + 2c4*f'²sin(2f)
 *   - 2sin(2f) - 2c4*sin(2f)*sin²f/r² = 0
 *
 * Divide by 2:
 *   f''(r² + 2c4*sin²f) + 2f'r + c4*f'²sin(2f)
 *   - sin(2f) - c4*sin(2f)*sin²f/r² = 0
 *
 * So: f'' = [sin(2f)(1 + c4*sin²f/r²) - 2f'r - c4*f'²sin(2f)]
 *           / (r² + 2c4*sin²f)
 */
static double f_rhs(double r, double f, double fp, double c4)
{
    double sf = sin(f);
    double s2f = sin(2.0 * f);
    double sf2 = sf * sf;
    double denom = r * r + 2.0 * c4 * sf2;

    if (denom < 1e-30) {
        /* At r=0: use L'Hopital / series expansion.
         * f(r) ≈ π - ar, so sf ≈ ar, sf2 ≈ a²r², s2f ≈ -2a²r²
         * denom ≈ r² + 2c4*a²r² = r²(1 + 2c4*a²)
         * numerator ≈ s2f(1 + c4*a²) - 2fp*r ≈ -2a²r²(1+c4*a²) + 2ar
         * f'' ≈ [-2a²(1+c4*a²) + 2a/r] / (1+2c4*a²)
         * At r=0 this diverges, so we use the limit:
         * From the series: f = π - ar - br³/6 + ..., f'' = -br
         * The ODE at r=0 reduces to: 0 = 0 (trivially satisfied)
         * At small r, the second-order term gives the cubic correction.
         */
        return 0.0;
    }

    double num = s2f * (1.0 + c4 * sf2 / (r * r)) - 2.0 * fp * r - c4 * fp * fp * s2f;
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
static double shoot(Profile *p, double a, double c4)
{
    double dr = p->dr;

    /* Initial conditions at r=eps (first grid point) */
    /* f(r) = pi - a*r + higher order terms near r=0 */
    /* From the ODE, the next term is: f = pi - ar + br³/3 where
     * b can be determined from the ODE. For simplicity, start with linear. */
    p->r[0] = 0.0;
    p->f[0] = M_PI;
    p->fp[0] = -a;

    /* Use special start for first step from r=0 */
    /* At r=0, f''(0) needs care. From the ODE expansion:
     * f = π - ar + cr³ + ...,  sin(f) = ar - ...,  sin²(f) = a²r² - ...
     * The ODE gives: f''(0) = 0 by symmetry.
     * So start with f(dr) = π - a*dr, f'(dr) = -a (to first order). */
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
        double k1v = f_rhs(r0, f0, v0, c4);

        double k2f = v0 + 0.5 * h * k1v;
        double k2v = f_rhs(r0 + 0.5*h, f0 + 0.5*h*k1f, v0 + 0.5*h*k1v, c4);

        double k3f = v0 + 0.5 * h * k2v;
        double k3v = f_rhs(r0 + 0.5*h, f0 + 0.5*h*k2f, v0 + 0.5*h*k2v, c4);

        double k4f = v0 + h * k3v;
        double k4v = f_rhs(r0 + h, f0 + h*k3f, v0 + h*k3v, c4);

        p->f[i+1] = f0 + (h / 6.0) * (k1f + 2*k2f + 2*k3f + k4f);
        p->fp[i+1] = v0 + (h / 6.0) * (k1v + 2*k2v + 2*k3v + k4v);
        p->r[i+1] = r0 + h;
    }

    return p->f[p->N];  /* f(R_max) — should be 0 */
}

/* ---------- Energy computation ---------- */

typedef struct { double E2, E4, Etotal, Q; } Energy;

static Energy compute_energy(const Profile *p, double e, double rho0)
{
    double inv_e2 = 1.0 / (e * e);
    double rho2 = rho0 * rho0;
    double rho4 = rho2 * rho2;
    double dr = p->dr;

    double e2_sum = 0.0, e4_sum = 0.0, q_sum = 0.0;

    for (int i = 0; i <= p->N; i++) {
        double r = p->r[i];
        double fi = p->f[i];
        double fpi = p->fp[i];
        double sf = sin(fi);
        double sf2 = sf * sf;

        double w = dr;
        if (i == 0 || i == p->N) w *= 0.5;

        e2_sum += (fpi * fpi * r * r + 2.0 * sf2) * w;
        if (r > 1e-14)
            e4_sum += (2.0 * fpi * fpi * sf2 + sf2 * sf2 / (r * r)) * w;
        q_sum += (-fpi * sf2) * w;
    }

    Energy en;
    en.E2 = 2.0 * M_PI * rho2 * e2_sum;
    en.E4 = 4.0 * M_PI * inv_e2 * rho4 * e4_sum;
    en.Etotal = en.E2 + en.E4;
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
    const char *outfile = "profile.dat";

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-Nr") == 0 && i+1 < argc) Nr = atoi(argv[++i]);
        else if (strcmp(argv[i], "-Rmax") == 0 && i+1 < argc) R_max = atof(argv[++i]);
        else if (strcmp(argv[i], "-e") == 0 && i+1 < argc) e_skyrme = atof(argv[++i]);
        else if (strcmp(argv[i], "-rho0") == 0 && i+1 < argc) rho0 = atof(argv[++i]);
        else if (strcmp(argv[i], "-o") == 0 && i+1 < argc) outfile = argv[++i];
        else {
            fprintf(stderr, "Usage: %s [-Nr N] [-Rmax R] [-e E] [-rho0 R] [-o outfile]\n", argv[0]);
            return 1;
        }
    }

    /* c4 = ratio of E4 prefactor to E2 prefactor in hedgehog energy:
     * E2 prefactor = 2*pi*rho0^2, E4 prefactor = 4*pi*rho0^4/e^2
     * ratio = 2*rho0^2/e^2 */
    double c4 = 2.0 * rho0 * rho0 / (e_skyrme * e_skyrme);

    printf("=== Hedgehog Skyrmion — Shooting Method ===\n");
    printf("Grid: Nr=%d, dr=%.6f, R_max=%.1f\n", Nr, R_max / Nr, R_max);
    printf("Parameters: e=%.4f (c4=%.6f), rho0=%.4f\n", e_skyrme, c4, rho0);
    printf("\n");

    Profile *p = profile_alloc(Nr, R_max);

    /* Bisection on the shooting parameter a = -f'(0)
     * a scales as 1/sqrt(c4) ~ e/sqrt(2), so adapt the search range */
    double a_lo = 0.1, a_hi = 20.0 / sqrt(c4);

    /* Find bracket: f(R_max) should change sign as a varies */
    double f_lo = shoot(p, a_lo, c4);
    double f_hi = shoot(p, a_hi, c4);
    printf("Bracket: a_lo=%.4f -> f(Rmax)=%.6f, a_hi=%.4f -> f(Rmax)=%.6f\n",
           a_lo, f_lo, a_hi, f_hi);

    if (f_lo * f_hi > 0) {
        /* Try to find bracket by scanning */
        printf("Scanning for bracket...\n");
        int found = 0;
        double step = a_hi / 500.0;
        for (double a = step; a <= a_hi; a += step) {
            double fend = shoot(p, a, c4);
            if (fend * f_lo < 0) {
                a_hi = a;
                f_hi = fend;
                a_lo = a - step;
                f_lo = shoot(p, a_lo, c4);
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
        double f_mid = shoot(p, a_mid, c4);

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
    shoot(p, a_best, c4);

    /* Compute energy */
    Energy en = compute_energy(p, e_skyrme, rho0);

    printf("\n=== Result ===\n");
    printf("a = -f'(0) = %.14f\n", a_best);
    printf("E_total = %.12e\n", en.Etotal);
    printf("E2      = %.12e\n", en.E2);
    printf("E4      = %.12e\n", en.E4);
    printf("Q       = %.8f\n", en.Q);
    printf("E2/E4   = %.6f (should be ~1 at Bogomolny)\n", en.E2 / en.E4);
    printf("Virial: E2 - E4 = %.6e (should be ~0)\n", en.E2 - en.E4);

    /* The FB bound in our normalization: E_FB = 6*sqrt(2)*pi^2*rho0^3/e
     * Derived by rescaling r = R*sqrt(c4) to map to standard Skyrme ODE,
     * where E_std_FB = 12*pi^2, giving E_FB = sqrt(2)/(2e) * 12*pi^2 * rho0^3 */
    double E_FB = 6.0 * sqrt(2.0) * M_PI * M_PI * rho0 * rho0 * rho0 / e_skyrme;
    printf("\nFaddeev-Bogomolny bound (6*sqrt(2)*pi^2*rho0^3/e): %.6e\n", E_FB);
    printf("Ratio E/E_FB = %.6f (standard: ~1.232)\n", en.Etotal / E_FB);

    /* CHPT mass formula */
    printf("\nCHPT mass formula: Mc^2 = 2*E4 = %.12e\n", 2.0 * en.E4);
    printf("  (On sigma model without potential, EV=ED=0)\n");
    printf("  (With Derrick virial E2=E4, so Mc^2 = 2*E4 = E_total)\n");

    /* Profile output */
    FILE *fp = fopen(outfile, "w");
    if (fp) {
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
