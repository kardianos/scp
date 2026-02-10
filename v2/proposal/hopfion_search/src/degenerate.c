/*
 * degenerate.c — Degenerate sector dynamics: 1/r gravity from B⁰p coupling
 *
 * Solves the static sourced Poisson/Yukawa equation for the pseudoscalar
 * field p (e₀₁₂₃ component) coupled to the baryon density B⁰:
 *
 *   -κ_D² ∇²p + μ_D² p = g_top B⁰(r)
 *
 * where B⁰ = -f'sin²f/(2π²r²) is the topological charge density.
 *
 * Solutions:
 *   μ_D = 0 (massless): p(r→∞) → g_top B/(4πκ_D² r) — 1/r gravitational
 *   μ_D > 0 (massive):  p(r→∞) → g_top B e^{-μr/κ}/(4πκ_D² r) — Yukawa
 *
 * Two-soliton interaction at distance D:
 *   U(D) = g_top² B₁B₂ e^{-mD}/(4πκ_D² D)  where m = μ_D/κ_D
 *   → G_eff = g_top²/(4πκ_D² E_sol²) for gravity
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ========== Profile I/O (from bound.c pattern) ========== */

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

/* ========== Baryon density ========== */

/*
 * B⁰(r) = -f'sin²f/(2π²r²)
 * At r=0: L'Hôpital with f(r)~π-ar gives sin²f/r² → a², f'→ -a
 *   B⁰(0) = a³/(2π²)
 */
static double baryon_density(double r, double f, double fp, double a) {
    if (r > 1e-10) {
        double sf = sin(f);
        return -(1.0/(2*M_PI*M_PI)) * sf*sf * fp / (r*r);
    } else {
        return a*a*a / (2*M_PI*M_PI);
    }
}

/* ========== Massless solver (μ_D = 0): direct Coulomb integration ========== */

/*
 * Equation: -κ_D² ∇²p = g_top B⁰
 * Radial: -κ_D² [p'' + 2p'/r] = g_top B⁰(r)
 *
 * Solution by enclosed-charge method (same as maxwell.c):
 *   E_p(r) = -p'(r) = (g_top/(4πκ_D²)) Q(r)/r²
 *   where Q(r) = 4π ∫₀ʳ B⁰(r')r'² dr'
 *
 * Potential:
 *   p(r) = p(R_max) + ∫_r^R_max E_p(r') dr'
 *   p(R_max) = g_top B/(4πκ_D² R_max)
 */
static void solve_massless(const Profile *prof, double g_top, double kappa2,
                           double *B0, double *Q_enc, double *p_field, int n)
{
    double dr = prof->dr;
    double a = -prof->fp[0];

    /* Compute B⁰ and enclosed charge Q(r) */
    for (int i = 0; i < n; i++) {
        B0[i] = baryon_density(prof->r[i], prof->f[i], prof->fp[i], a);
    }

    Q_enc[0] = 0;
    for (int i = 1; i < n; i++) {
        Q_enc[i] = Q_enc[i-1] + 0.5*dr * (
            4*M_PI * B0[i-1] * prof->r[i-1]*prof->r[i-1] +
            4*M_PI * B0[i]   * prof->r[i]*prof->r[i]
        );
    }

    /* p'(r) = -(g_top/(4πκ²)) Q(r)/r² */
    double coeff = g_top / (4*M_PI * kappa2);

    /* Boundary: p(R_max) = g_top × Q(∞) / (4πκ² R_max) */
    double R_max = prof->r[n-1];
    p_field[n-1] = coeff * Q_enc[n-1] / R_max;

    /* Integrate inward: p(r) = p(r+dr) + ∫_r^{r+dr} (g_top/(4πκ²)) Q/r² dr */
    for (int i = n-2; i >= 0; i--) {
        double r1 = prof->r[i+1];
        double r0 = prof->r[i];
        double E1 = (r1 > 1e-10) ? coeff * Q_enc[i+1] / (r1*r1) : 0;
        double E0 = (r0 > 1e-10) ? coeff * Q_enc[i]   / (r0*r0) : 0;
        p_field[i] = p_field[i+1] + 0.5*dr * (E0 + E1);
    }
}

/* ========== Massive solver (μ_D > 0): Green's function method ========== */

/*
 * Equation: -κ_D² ∇²p + μ_D² p = g_top B⁰
 * Radial: -κ_D² [p'' + 2p'/r] + μ_D² p = g_top B⁰(r)
 *
 * Green's function for [-∇² + m²]G = δ³(x):
 *   G(r) = e^{-mr}/(4πr)  where m = μ_D/κ_D
 *
 * Solution: p(r) = (g_top/κ_D²) ∫ G(|x-x'|) B⁰(r') d³x'
 *
 * For spherical source (angle-averaged):
 *   p(r) = (g_top/(κ_D² m r)) [e^{-mr} I_sinh(r) + sinh(mr) I_exp(r)]
 *
 * where:
 *   I_sinh(r) = ∫₀ʳ B⁰(r')r' sinh(mr') dr'
 *   I_exp(r)  = ∫_r^∞ B⁰(r')r' e^{-mr'} dr'
 *
 * At r=0: p(0) = (g_top/κ_D²) ∫₀^∞ B⁰(r')r' e^{-mr'}/r' dr'×r'
 *              = (g_top/κ_D²) I_exp(0)  [finite]
 *
 * Overflow guard: use factored forms
 *   e^{-mr}sinh(mr') = (1/2)[e^{-m(r-r')} - e^{-m(r+r')}]  for r'<r
 *   sinh(mr)e^{-mr'} = (1/2)[e^{-m(r'-r)} - e^{-m(r'+r)}]  for r'>r
 * This avoids overflow when mr is large.
 */
static void solve_massive(const Profile *prof, double g_top, double kappa2,
                          double mu_D, double *B0, double *p_field, int n)
{
    double dr = prof->dr;
    double a = -prof->fp[0];
    double m = mu_D / sqrt(kappa2);  /* effective mass parameter */

    /* Compute B⁰ */
    for (int i = 0; i < n; i++) {
        B0[i] = baryon_density(prof->r[i], prof->f[i], prof->fp[i], a);
    }

    /* Forward sweep: I_sinh(r_i) = ∫₀^{r_i} B⁰(r')r' sinh(mr') dr' */
    double *I_sinh = calloc(n, sizeof(double));
    I_sinh[0] = 0;
    for (int i = 1; i < n; i++) {
        double r0 = prof->r[i-1], r1 = prof->r[i];
        double f0 = B0[i-1] * r0 * sinh(m*r0);
        double f1 = B0[i]   * r1 * sinh(m*r1);
        /* Guard overflow: if m*r > 500, sinh overflows but product may be ~0 */
        if (m*r0 > 500 || fabs(B0[i-1]) < 1e-300) f0 = 0;
        if (m*r1 > 500 || fabs(B0[i])   < 1e-300) f1 = 0;
        I_sinh[i] = I_sinh[i-1] + 0.5*dr*(f0 + f1);
    }

    /* Backward sweep: I_exp(r_i) = ∫_{r_i}^∞ B⁰(r')r' e^{-mr'} dr' */
    double *I_exp = calloc(n, sizeof(double));
    I_exp[n-1] = 0;
    for (int i = n-2; i >= 0; i--) {
        double r0 = prof->r[i], r1 = prof->r[i+1];
        double f0 = B0[i]   * r0 * exp(-m*r0);
        double f1 = B0[i+1] * r1 * exp(-m*r1);
        I_exp[i] = I_exp[i+1] + 0.5*dr*(f0 + f1);
    }

    /* Assemble p(r) = (g_top/(κ² m r)) [e^{-mr} I_sinh(r) + sinh(mr) I_exp(r)] */
    double prefac = g_top / (kappa2 * m);
    for (int i = 0; i < n; i++) {
        double r = prof->r[i];
        if (r < 1e-10) {
            /* At r→0: p(0) = (g_top/κ²) I_exp(0)
             * because lim_{r→0} [e^{-mr}I_sinh(r) + sinh(mr)I_exp(r)]/(mr)
             *   = lim I_sinh(r)/(mr) + I_exp(r) [sinh(mr)/(mr)→1]
             *   = 0 + I_exp(0) */
            p_field[i] = (g_top/kappa2) * I_exp[0];
        } else {
            double emr = exp(-m*r);
            double smr = sinh(m*r);
            /* Use safe form when mr is large */
            if (m*r > 500) {
                /* e^{-mr}×I_sinh ≈ 0, sinh(mr)×I_exp ≈ (e^{mr}/2)×I_exp */
                p_field[i] = prefac * 0.5 * exp(m*r) * I_exp[i] / r;
                /* But if I_exp[i] is already e^{-mr}×(...), this is fine */
            } else {
                p_field[i] = prefac * (emr * I_sinh[i] + smr * I_exp[i]) / r;
            }
        }
    }

    free(I_sinh);
    free(I_exp);
}

/* ========== Self-energy ========== */

/*
 * E_self = (g_top/2) ∫ B⁰(r) p(r) 4πr² dr
 *
 * For massless: E_self = g_top²/(8πκ²) ∫₀^∞ Q(r)²/r² dr × (1/r) ...
 * Actually just integrate directly.
 */
static double compute_self_energy(const Profile *prof, const double *B0,
                                  const double *p_field, double g_top, int n)
{
    double dr = prof->dr;
    double E = 0;
    for (int i = 1; i < n; i++) {
        double r0 = prof->r[i-1], r1 = prof->r[i];
        double f0 = B0[i-1] * p_field[i-1] * r0*r0;
        double f1 = B0[i]   * p_field[i]   * r1*r1;
        E += 0.5*dr*(f0 + f1);
    }
    return g_top * 2*M_PI * E;  /* 4π × (1/2) × g_top */
}

/* ========== ODE residual check ========== */

/*
 * Check: -κ² [p'' + 2p'/r] + μ² p = g_top B⁰
 * Use 4th-order finite differences for p'' and p'.
 */
static double ode_residual(const Profile *prof, const double *p_field,
                           const double *B0, double kappa2, double mu2,
                           double g_top, int n)
{
    double dr = prof->dr;
    double max_res = 0, max_rhs = 0;

    for (int i = 2; i < n-2; i++) {
        double r = prof->r[i];
        if (r < 1e-6) continue;

        /* 4th-order central diffs */
        double pp = (-p_field[i+2] + 16*p_field[i+1] - 30*p_field[i]
                     + 16*p_field[i-1] - p_field[i-2]) / (12*dr*dr);
        double p1 = (-p_field[i+2] + 8*p_field[i+1]
                     - 8*p_field[i-1] + p_field[i-2]) / (12*dr);

        double laplacian = pp + 2*p1/r;
        double lhs = -kappa2 * laplacian + mu2 * p_field[i];
        double rhs = g_top * B0[i];
        double res = fabs(lhs - rhs);
        if (res > max_res) max_res = res;
        if (fabs(rhs) > max_rhs) max_rhs = fabs(rhs);
    }

    return (max_rhs > 1e-30) ? max_res / max_rhs : max_res;
}

/* ========== Two-soliton interaction ========== */

/*
 * Point-source approximation (valid for D >> R_soliton):
 *   U_point(D) = g_top² B₁B₂ e^{-mD} / (4πκ² D)
 *
 * where m = μ_D/κ_D (0 for massless).
 *
 * Full overlap integral (angular average):
 *   U(D) = g_top ∫ B⁰(r) p_2(|x-x₂|) d³x
 *
 * where p_2 is the field from soliton 2 centered at distance D.
 * For spherical solitons, angular average:
 *   p_2_avg(r,D) = (1/2) ∫₋₁¹ p(√(r²+D²+2rDcosθ)) d(cosθ)
 *
 * Then: U(D) = g_top ∫₀^∞ B⁰(r) p_2_avg(r,D) 4πr² dr
 *
 * We evaluate p_2_avg by Gauss-Legendre quadrature over cosθ.
 */
static double interaction_point(double D, double g_top, double kappa2,
                                double m, int B1, int B2)
{
    if (D < 1e-10) return 0;
    return g_top * g_top * B1 * B2 * exp(-m*D) / (4*M_PI * kappa2 * D);
}

/* Interpolate p(r) from tabulated values */
static double interp_p(const Profile *prof, const double *p_field, double r, int n)
{
    if (r <= prof->r[0]) return p_field[0];
    if (r >= prof->r[n-1]) {
        /* Extrapolate: for massless, p ~ 1/r; for massive, p ~ e^{-mr}/r */
        /* Use last two points */
        return p_field[n-1] * prof->r[n-1] / r;  /* 1/r extrapolation */
    }
    double idx = r / prof->dr;
    int i0 = (int)idx;
    if (i0 >= n-1) i0 = n-2;
    double t = idx - i0;
    return (1-t)*p_field[i0] + t*p_field[i0+1];
}

static double interaction_full(const Profile *prof, const double *B0,
                               const double *p_field, double D,
                               double g_top, int n, int n_gl)
{
    double dr = prof->dr;

    /* Gauss-Legendre nodes and weights for cosθ ∈ [-1,1] */
    double *gl_x = malloc(n_gl * sizeof(double));
    double *gl_w = malloc(n_gl * sizeof(double));
    for (int i = 0; i < n_gl; i++) {
        double x = cos(M_PI * (i + 0.75) / (n_gl + 0.5));
        for (int iter = 0; iter < 20; iter++) {
            double p0 = 1, p1 = x;
            for (int j = 2; j <= n_gl; j++) {
                double p2 = ((2*j-1)*x*p1 - (j-1)*p0) / j;
                p0 = p1; p1 = p2;
            }
            double dp = n_gl * (x*p1 - p0) / (x*x - 1);
            x -= p1 / dp;
        }
        gl_x[i] = x;
        double p0 = 1, p1 = x;
        for (int j = 2; j <= n_gl; j++) {
            double p2 = ((2*j-1)*x*p1 - (j-1)*p0) / j;
            p0 = p1; p1 = p2;
        }
        gl_w[i] = 2.0 * (1 - x*x) / ((double)n_gl*n_gl * p0*p0);
    }

    /* U(D) = g_top ∫₀^∞ B⁰(r) <p(|x-D ẑ|)>_angle 4πr² dr */
    double U = 0;
    for (int i = 0; i < n; i++) {
        double r = prof->r[i];
        if (fabs(B0[i]) < 1e-30) continue;

        /* Angular average of p(√(r²+D²+2rDcosθ)) */
        double p_avg = 0;
        for (int j = 0; j < n_gl; j++) {
            double ct = gl_x[j];
            double dist2 = r*r + D*D + 2*r*D*ct;
            if (dist2 < 0) dist2 = 0;
            double dist = sqrt(dist2);
            p_avg += gl_w[j] * interp_p(prof, p_field, dist, n);
        }
        p_avg *= 0.5;  /* ∫₋₁¹ dμ / 2 normalization */

        double w = dr * ((i == 0 || i == n-1) ? 0.5 : 1.0);
        U += w * B0[i] * p_avg * 4*M_PI * r*r;
    }

    free(gl_x);
    free(gl_w);
    return g_top * U;
}

/* ========== Main ========== */

int main(int argc, char *argv[])
{
    const char *profile_file = "data/profiles/profile_sigma_e1.dat";
    double g_top = 1.0;
    double kappa_D = 1.0;
    double mu_D = 0.0;
    int scan = 0;
    double D_min = 2.0, D_max = 20.0;
    int n_D = 10;

    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "-profile") && i+1 < argc) profile_file = argv[++i];
        else if (!strcmp(argv[i], "-gtop") && i+1 < argc) g_top = atof(argv[++i]);
        else if (!strcmp(argv[i], "-kappa") && i+1 < argc) kappa_D = atof(argv[++i]);
        else if (!strcmp(argv[i], "-mu") && i+1 < argc) mu_D = atof(argv[++i]);
        else if (!strcmp(argv[i], "-scan")) scan = 1;
        else if (!strcmp(argv[i], "-Dmin") && i+1 < argc) D_min = atof(argv[++i]);
        else if (!strcmp(argv[i], "-Dmax") && i+1 < argc) D_max = atof(argv[++i]);
        else if (!strcmp(argv[i], "-nD") && i+1 < argc) n_D = atoi(argv[++i]);
        else {
            fprintf(stderr, "Usage: %s [-profile PATH] [-gtop G] [-kappa K] [-mu M]\n"
                    "       [-scan] [-Dmin D] [-Dmax D] [-nD N]\n", argv[0]);
            return 1;
        }
    }

    double kappa2 = kappa_D * kappa_D;
    double mu2 = mu_D * mu_D;
    double m = (kappa_D > 1e-30) ? mu_D / kappa_D : 0;  /* effective mass */

    /* Read soliton profile */
    Profile prof;
    if (read_profile(profile_file, &prof) < 0) return 1;
    int n = prof.n;

    printf("============================================================\n");
    printf(" Degenerate Sector: B⁰p Coupling → Gravitational Potential\n");
    printf("============================================================\n\n");
    printf("Profile: %s (%d points, R_max=%.2f, dr=%.6f)\n",
           profile_file, n, prof.r[n-1], prof.dr);
    printf("Parameters: g_top=%.6e, κ_D=%.4f, μ_D=%.6f\n", g_top, kappa_D, mu_D);
    printf("Effective mass: m = μ_D/κ_D = %.6f\n", m);
    printf("Regime: %s\n\n", mu_D < 1e-10 ? "MASSLESS (1/r potential)"
                                           : "MASSIVE (Yukawa potential)");

    /* Allocate arrays */
    double *B0      = calloc(n, sizeof(double));
    double *Q_enc   = calloc(n, sizeof(double));
    double *p_field = calloc(n, sizeof(double));

    /* Compute baryon density (needed for both solvers) */
    double a = -prof.fp[0];
    for (int i = 0; i < n; i++)
        B0[i] = baryon_density(prof.r[i], prof.f[i], prof.fp[i], a);

    /* Solve */
    if (mu_D < 1e-10) {
        solve_massless(&prof, g_top, kappa2, B0, Q_enc, p_field, n);
    } else {
        solve_massive(&prof, g_top, kappa2, mu_D, B0, p_field, n);
        /* Also compute Q_enc for diagnostics */
        Q_enc[0] = 0;
        for (int i = 1; i < n; i++)
            Q_enc[i] = Q_enc[i-1] + 0.5*prof.dr * (
                4*M_PI * B0[i-1]*prof.r[i-1]*prof.r[i-1] +
                4*M_PI * B0[i]*prof.r[i]*prof.r[i]
            );
    }

    /* ===== Verification ===== */
    printf("--- Verification ---\n");

    /* 1. Total baryon charge */
    double Q_total = Q_enc[n-1];
    printf("Q(∞) = %.6f  (expected: 1.000)\n", Q_total);

    /* 2. 1/r tail verification (massless) */
    if (mu_D < 1e-10) {
        printf("\n1/r tail check: p(r)×r should plateau to g_top/(4πκ²) = %.6e\n",
               g_top / (4*M_PI*kappa2));
        printf("  %8s  %14s  %14s  %10s\n", "r", "p(r)", "p(r)*r", "expected");
        double expected = g_top * Q_total / (4*M_PI * kappa2);
        for (int i = n/2; i < n; i += n/10) {
            printf("  %8.3f  %14.6e  %14.6e  %10.6e\n",
                   prof.r[i], p_field[i], p_field[i]*prof.r[i], expected);
        }
    } else {
        /* Yukawa tail: p(r) × r × e^{mr} should plateau */
        printf("\nYukawa tail: p(r)×r×e^{mr} should plateau to g_top/(4πκ²)\n");
        printf("  %8s  %14s  %14s  %10s\n", "r", "p(r)", "p*r*e^{mr}", "expected");
        double expected = g_top * Q_total / (4*M_PI * kappa2);
        for (int i = n/2; i < n; i += n/10) {
            double r = prof.r[i];
            if (m*r < 500)
                printf("  %8.3f  %14.6e  %14.6e  %10.6e\n",
                       r, p_field[i], p_field[i]*r*exp(m*r), expected);
        }
    }

    /* 3. ODE residual */
    double res = ode_residual(&prof, p_field, B0, kappa2, mu2, g_top, n);
    printf("\nODE residual (max |LHS-RHS| / max |RHS|): %.6e\n", res);

    /* ===== Self-energy ===== */
    printf("\n--- Self-energy ---\n");
    double E_self = compute_self_energy(&prof, B0, p_field, g_top, n);
    double E_sol = 103.13;  /* code units (e=1, ρ₀=1) */
    double E_code_to_MeV = 9.098;
    printf("E_self = %.6e (code units)\n", E_self);
    printf("E_self = %.6e MeV\n", E_self * E_code_to_MeV);
    printf("ΔE/E_sol = %.6e\n", E_self / E_sol);

    /* ===== p(r) profile ===== */
    printf("\n--- p(r) profile ---\n");
    printf("  %8s  %14s  %14s  %14s\n", "r", "B0(r)", "p(r)", "p(r)*r");
    for (int i = 0; i < n; i += (n > 200 ? n/50 : 1)) {
        printf("  %8.4f  %14.6e  %14.6e  %14.6e\n",
               prof.r[i], B0[i], p_field[i], p_field[i]*prof.r[i]);
    }

    /* ===== Two-soliton interaction ===== */
    printf("\n--- Two-soliton interaction (B₁=B₂=1) ---\n");
    printf("  %8s  %14s  %14s  %10s\n", "D", "U_full(D)", "U_point(D)", "ratio");
    int n_gl = 32;  /* Gauss-Legendre points for angular average */

    double D_test[] = {2.0, 3.0, 5.0, 7.0, 10.0, 15.0, 20.0};
    int n_test = sizeof(D_test)/sizeof(D_test[0]);
    for (int j = 0; j < n_test; j++) {
        double D = D_test[j];
        if (D > prof.r[n-1]) break;

        double U_full = interaction_full(&prof, B0, p_field, D, g_top, n, n_gl);
        double U_pt = interaction_point(D, g_top, kappa2, m, 1, 1);
        double ratio = (fabs(U_pt) > 1e-30) ? U_full / U_pt : 0;
        printf("  %8.2f  %14.6e  %14.6e  %10.6f\n", D, U_full, U_pt, ratio);
    }

    /* ===== G_eff and required g_top for G_Newton ===== */
    printf("\n--- Effective gravitational constant ---\n");
    printf("G_eff = g_top² / (4πκ²E_sol²)\n");
    double G_eff = g_top*g_top / (4*M_PI * kappa2 * E_sol*E_sol);
    printf("G_eff = %.6e (code units, with g_top=%.2e)\n", G_eff, g_top);

    /* Physical G_Newton in code units */
    /* G_N = 6.674e-11 m³/(kg·s²) */
    /* In code: G_N × E_sol²/(ℏc) ≈ 5.9e-39 */
    double G_Newton_code = 5.9e-39;
    printf("G_Newton (code) = %.2e\n", G_Newton_code);

    /* Required g_top for G_Newton: g_top² = G_N × 4πκ²E_sol² */
    double gtop_required = sqrt(G_Newton_code * 4*M_PI * kappa2 * E_sol*E_sol);
    printf("Required g_top for G_Newton: %.6e\n", gtop_required);
    printf("Hierarchy ratio: g_top_grav/g_top_nuclear ~ %.2e\n", gtop_required);

    /* ===== μ_D scan (if requested) ===== */
    if (scan) {
        printf("\n--- μ_D scan: transition from 1/r to Yukawa ---\n");
        printf("  %10s  %10s  %14s  %14s  %14s  %14s\n",
               "mu_D", "m=mu/kappa", "p(0)", "p(5)", "p(10)", "E_self");

        double mu_vals[] = {0, 0.001, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2,
                            0.398, 0.5, 1.0, 2.0};
        int n_mu = sizeof(mu_vals)/sizeof(mu_vals[0]);

        for (int im = 0; im < n_mu; im++) {
            double mu_v = mu_vals[im];
            double *p_v = calloc(n, sizeof(double));
            double *B0_v = calloc(n, sizeof(double));
            double *Q_v  = calloc(n, sizeof(double));

            for (int i = 0; i < n; i++)
                B0_v[i] = B0[i];

            if (mu_v < 1e-10)
                solve_massless(&prof, g_top, kappa2, B0_v, Q_v, p_v, n);
            else
                solve_massive(&prof, g_top, kappa2, mu_v, B0_v, p_v, n);

            double E_s = compute_self_energy(&prof, B0_v, p_v, g_top, n);

            /* Find p at r=5 and r=10 */
            int i5 = (int)(5.0/prof.dr);
            int i10 = (int)(10.0/prof.dr);
            if (i5 >= n) i5 = n-1;
            if (i10 >= n) i10 = n-1;

            printf("  %10.4f  %10.4f  %14.6e  %14.6e  %14.6e  %14.6e\n",
                   mu_v, mu_v/kappa_D, p_v[0], p_v[i5], p_v[i10], E_s);

            free(p_v); free(B0_v); free(Q_v);
        }
    }

    /* Custom D range interaction table */
    if (D_min < D_max && n_D > 0) {
        printf("\n--- Interaction table (D from %.1f to %.1f, %d points) ---\n",
               D_min, D_max, n_D);
        printf("  %8s  %14s  %14s  %10s  %14s\n",
               "D", "U_full", "U_point", "ratio", "F=-dU/dD");

        double dD = (n_D > 1) ? (D_max - D_min) / (n_D - 1) : 0;
        double *U_arr = malloc(n_D * sizeof(double));
        double *D_arr = malloc(n_D * sizeof(double));

        for (int j = 0; j < n_D; j++) {
            double D = D_min + j * dD;
            D_arr[j] = D;
            if (D > prof.r[n-1]) {
                U_arr[j] = interaction_point(D, g_top, kappa2, m, 1, 1);
            } else {
                U_arr[j] = interaction_full(&prof, B0, p_field, D, g_top, n, n_gl);
            }
        }

        for (int j = 0; j < n_D; j++) {
            double D = D_arr[j];
            double U_pt = interaction_point(D, g_top, kappa2, m, 1, 1);
            double ratio = (fabs(U_pt) > 1e-30) ? U_arr[j] / U_pt : 0;

            /* Force by central finite difference */
            double force = 0;
            if (j > 0 && j < n_D-1)
                force = -(U_arr[j+1] - U_arr[j-1]) / (2*dD);
            else if (j == 0 && n_D > 1)
                force = -(U_arr[1] - U_arr[0]) / dD;
            else if (j == n_D-1 && n_D > 1)
                force = -(U_arr[n_D-1] - U_arr[n_D-2]) / dD;

            printf("  %8.3f  %14.6e  %14.6e  %10.6f  %14.6e\n",
                   D, U_arr[j], U_pt, ratio, force);
        }

        free(U_arr); free(D_arr);
    }

    /* ===== Physical summary ===== */
    printf("\n============================================================\n");
    printf(" Physical Summary\n");
    printf("============================================================\n\n");

    printf("Lagrangian: L_D = (κ²/2c²)|ṗ|² - (κ²/2)|∇p|² - (μ²/2)p² + g_top B⁰ p\n\n");

    if (mu_D < 1e-10) {
        printf("MASSLESS case (μ_D=0):\n");
        printf("  p(r→∞) → g_top B / (4πκ² r) — Newtonian 1/r potential\n");
        printf("  U(D) = g_top² B₁B₂ / (4πκ² D) — gravitational attraction\n");
        printf("  G_eff = g_top² / (4πκ² E_sol²) = %.6e\n", G_eff);
        printf("  Required g_top for G_Newton: %.6e\n", gtop_required);
        printf("  → Reflects gravitational hierarchy: e₀² ~ 10⁻³⁸\n");
    } else {
        double range_fm = kappa_D / mu_D * 0.5624;  /* code_L × fm/code_L */
        printf("MASSIVE case (μ_D=%.4f):\n", mu_D);
        printf("  p(r→∞) → g_top B e^{-mr} / (4πκ² r) — Yukawa potential\n");
        printf("  Range: 1/m = κ/μ = %.4f code = %.4f fm\n", kappa_D/mu_D, range_fm);
        printf("  U(D) = g_top² B₁B₂ e^{-mD} / (4πκ² D)\n");
        if (fabs(mu_D - 0.398) < 0.01)
            printf("  → Physical pion mass: nuclear-range force (~1.4 fm)\n");
    }

    printf("\np(0) = %.6e\n", p_field[0]);
    printf("E_self = %.6e code = %.4f MeV\n", E_self, E_self * E_code_to_MeV);
    printf("ΔM/M_soliton = %.6e\n", E_self / E_sol);

    /* Cleanup */
    free(B0); free(Q_enc); free(p_field);
    free_profile(&prof);

    return 0;
}
