/*
 * soler.c — Soler model (nonlinear Dirac) soliton solver
 *
 * Model E from MODELS.md
 * Lagrangian: L = ψ̄(iγ∂ - m)ψ + λ(ψ̄ψ)²
 *
 * Ansatz: ψ = [f(r)χ_κ, ig(r)(σ·r̂)χ_κ] e^{-iωt}
 * Scalar density: ψ̄ψ = f² - g²
 *
 * EOM: iγ∂ψ = (m - 2λ(ψ̄ψ))ψ  →  m_eff = m - 2λS
 * Radial equations for ground state κ = -1 (λ_code = 2λ_Lagrangian):
 *   f' = (m + ω - λ(f²-g²)) g
 *   g' + 2g/r = (m - ω - λ(f²-g²)) f
 *
 * f = upper/large component: f(0) = f₀ (shooting parameter)
 * g = lower/small component: g(0) = 0
 *
 * Asymptotic: f,g ~ e^{-κr} with κ = √(m²-ω²) for ω < m.
 * g approaches 0 from below (standard relativistic behavior).
 *
 * Shooting: bisect on f₀ until f → 0 without zero crossing.
 *   f₀ too small → f diverges to +∞ (growing mode, no zero crossing)
 *   f₀ too large → f crosses zero (overshoots decay)
 *
 * Usage:
 *   soler -omega 0.9                 (single frequency)
 *   soler -scan                      (scan omega from 0.1 to 0.99)
 */

#define _DEFAULT_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define NMAX 200001

static double f_arr[NMAX], g_arr[NMAX];

/* Parameters */
static double par_m = 1.0;       /* fermion mass */
static double par_omega = 0.9;   /* frequency */
static double par_lambda = 1.0;  /* self-coupling */

/* RHS of Soler system for κ = -1 */
static void rhs(double r, double f, double g, double *df, double *dg)
{
    double S = f * f - g * g;  /* scalar density ψ̄ψ */
    *df = (par_m + par_omega - par_lambda * S) * g;
    if (r < 1e-12) {
        /* Taylor: g ~ g₁r, so g' + 2g/r → 3g₁ = (m-ω-λS)f at r=0 */
        *dg = (par_m - par_omega - par_lambda * S) * f / 3.0;
    } else {
        *dg = -2.0 * g / r + (par_m - par_omega - par_lambda * S) * f;
    }
}

/*
 * RK4 integration from r=0 to rmax.
 * Returns last valid index (N if complete, < N if diverges).
 */
static int integrate(double f0, int N, double rmax)
{
    double h = rmax / N;

    f_arr[0] = f0;
    g_arr[0] = 0;

    for (int i = 0; i < N; i++) {
        double r = i * h;
        double f = f_arr[i];
        double g = g_arr[i];

        double k1f, k1g, k2f, k2g, k3f, k3g, k4f, k4g;
        rhs(r, f, g, &k1f, &k1g);
        rhs(r + 0.5 * h, f + 0.5 * h * k1f, g + 0.5 * h * k1g, &k2f, &k2g);
        rhs(r + 0.5 * h, f + 0.5 * h * k2f, g + 0.5 * h * k2g, &k3f, &k3g);
        rhs(r + h, f + h * k3f, g + h * k3g, &k4f, &k4g);

        f_arr[i + 1] = f + h * (k1f + 2 * k2f + 2 * k3f + k4f) / 6.0;
        g_arr[i + 1] = g + h * (k1g + 2 * k2g + 2 * k3g + k4g) / 6.0;

        if (fabs(f_arr[i + 1]) > 1e10 || fabs(g_arr[i + 1]) > 1e10)
            return i + 1;
    }
    return N;
}

/*
 * Classify solution: -1 if f crosses zero in the physical region
 * (f₀ too large), +1 if f stays positive (f₀ too small).
 *
 * Only check up to the minimum of f²+g² — beyond that,
 * the growing exponential mode dominates and is unphysical.
 */
static int classify(int nv)
{
    /* Find minimum of f²+g² (end of physical region) */
    double amp_min = f_arr[0] * f_arr[0] + g_arr[0] * g_arr[0];
    int imin = 0;
    for (int i = 1; i <= nv; i++) {
        double amp = f_arr[i] * f_arr[i] + g_arr[i] * g_arr[i];
        if (amp < amp_min) {
            amp_min = amp;
            imin = i;
        }
    }

    /* Check for zero crossing only in physical region */
    for (int i = 1; i <= imin; i++) {
        if (f_arr[i] < 0) return -1;
    }
    return +1;
}

/*
 * Find f₀ by bisection for current par_omega.
 * Returns f₀ > 0 on success, -1 on failure.
 */
static double find_f0(int N, double rmax)
{
    /* Scan to find bracket: first +1 → -1 transition */
    double last_pos = -1;
    double first_neg = -1;
    int found_pos = 0;

    for (double f_scan = 0.001; f_scan <= 30.0; f_scan *= 1.01) {
        int nv = integrate(f_scan, N, rmax);
        int s = classify(nv);
        if (s == +1) {
            last_pos = f_scan;
            found_pos = 1;
        } else if (s == -1 && found_pos) {
            first_neg = f_scan;
            break;
        }
    }

    if (last_pos < 0 || first_neg < 0) return -1;

    double f_lo = last_pos;   /* +1: f stays positive */
    double f_hi = first_neg;  /* -1: f crosses zero */

    /* Bisection: 80 iterations for ~24 digits */
    for (int iter = 0; iter < 80; iter++) {
        double f_mid = 0.5 * (f_lo + f_hi);
        int nv = integrate(f_mid, N, rmax);
        int s = classify(nv);
        if (s == +1)
            f_lo = f_mid;
        else
            f_hi = f_mid;
    }

    double f0_best = 0.5 * (f_lo + f_hi);
    integrate(f0_best, N, rmax);

    /* Clip growing tail: find minimum of f²+g² in valid range, zero beyond */
    double amp_min = f_arr[0] * f_arr[0] + g_arr[0] * g_arr[0];
    int imin = 0;
    for (int i = 1; i <= N; i++) {
        double amp = f_arr[i] * f_arr[i] + g_arr[i] * g_arr[i];
        if (amp < amp_min) {
            amp_min = amp;
            imin = i;
        }
    }
    for (int i = imin + 1; i <= N; i++) {
        f_arr[i] = 0;
        g_arr[i] = 0;
    }

    return f0_best;
}

/*
 * Compute particle number Q = 4π∫(f²+g²)r²dr and R_rms.
 * Integrates only up to clip point (where f,g have been zeroed).
 */
static void compute_props(int N, double rmax,
                          double *Q_out, double *R_out)
{
    double h = rmax / N;
    double Q = 0, R2 = 0;

    for (int i = 0; i <= N; i++) {
        double r = i * h;
        double rho = f_arr[i] * f_arr[i] + g_arr[i] * g_arr[i];

        /* Simpson's rule weight */
        double w;
        if (i == 0 || i == N)
            w = 1.0 / 3.0;
        else if (i % 2 == 1)
            w = 4.0 / 3.0;
        else
            w = 2.0 / 3.0;

        Q += w * rho * r * r * h;
        R2 += w * rho * r * r * r * r * h;
    }

    *Q_out = 4.0 * M_PI * Q;
    *R_out = (Q > 0) ? sqrt(R2 / Q) : 0;
}

static void run_single(double omega, double rmax, int N, const char *outdir)
{
    par_omega = omega;
    double kappa = sqrt(par_m * par_m - par_omega * par_omega);

    /* Adjust domain: need at least 10 decay lengths */
    if (10.0 / kappa > rmax) rmax = 10.0 / kappa;
    if (rmax > 500) rmax = 500;

    printf("# Soler model: omega=%.6f, m=%.6f, lambda=%.6f\n",
           par_omega, par_m, par_lambda);
    printf("# rmax=%.1f, N=%d, h=%.6f, kappa=%.6f\n",
           rmax, N, rmax / N, kappa);

    double f0 = find_f0(N, rmax);
    if (f0 < 0) {
        printf("# FAILED: No bracket found for omega=%.6f\n", par_omega);
        return;
    }

    double Q, R;
    compute_props(N, rmax, &Q, &R);

    double E_bind = par_omega - par_m;
    double M = par_omega * Q;

    printf("#\n# RESULT: Soler soliton (nonlinear Dirac)\n");
    printf("#   f(0)      = %.8f  (upper component amplitude)\n", f0);
    printf("#   omega     = %.8f\n", par_omega);
    printf("#   E_bind    = %.8f  (omega - m, binding energy per particle)\n", E_bind);
    printf("#   Q         = %.6f  (particle number 4π∫(f²+g²)r²dr)\n", Q);
    printf("#   M=omega*Q = %.6f  (total mass-energy)\n", M);
    printf("#   R_rms     = %.6f\n", R);
    printf("#   kappa     = %.6f  (decay rate √(m²-ω²))\n", kappa);
    printf("#   S(0)=f₀²  = %.6f  (central scalar density ψ̄ψ)\n", f0 * f0);
    printf("#\n# POSITIVE: Nonlinear Dirac self-interaction provides confinement.\n");

    if (outdir) {
        char fname[256];
        snprintf(fname, sizeof(fname), "%s/soler_omega%.3f.dat", outdir, par_omega);
        FILE *fp = fopen(fname, "w");
        if (fp) {
            fprintf(fp, "# Soler soliton: omega=%.6f m=%.6f lambda=%.6f f0=%.8f\n",
                    par_omega, par_m, par_lambda, f0);
            fprintf(fp, "# r  f  g  rho=f²+g²  S=f²-g²\n");
            double h = rmax / N;
            int step = (N > 5000) ? N / 5000 : 1;
            for (int i = 0; i <= N; i += step) {
                double r = i * h;
                double f = f_arr[i], g = g_arr[i];
                fprintf(fp, "%.6f %.10e %.10e %.10e %.10e\n",
                        r, f, g, f * f + g * g, f * f - g * g);
            }
            fclose(fp);
            printf("#   Profile: %s\n", fname);
        }
    }
}

static void run_scan(double rmax, int N, const char *outdir)
{
    double omega_vals[] = {0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99};
    int nvals = sizeof(omega_vals) / sizeof(omega_vals[0]);

    printf("# Soler model scan: m=%.6f, lambda=%.6f\n", par_m, par_lambda);
    printf("# %10s  %12s  %12s  %12s  %12s  %8s\n",
           "omega", "f(0)", "E_bind", "Q", "M=omega*Q", "R_rms");

    for (int k = 0; k < nvals; k++) {
        par_omega = omega_vals[k];
        double rmax_k = rmax;
        double kappa = sqrt(par_m * par_m - par_omega * par_omega);
        if (kappa > 0 && 10.0 / kappa > rmax_k)
            rmax_k = 10.0 / kappa;
        if (rmax_k > 500)
            rmax_k = 500;

        double f0 = find_f0(N, rmax_k);
        if (f0 < 0) {
            printf("  %10.4f  FAILED\n", par_omega);
            continue;
        }

        double Q, R;
        compute_props(N, rmax_k, &Q, &R);
        double E_bind = par_omega - par_m;

        printf("  %10.4f  %12.6f  %12.8f  %12.6f  %12.6f  %8.4f  BOUND\n",
               par_omega, f0, E_bind, Q, par_omega * Q, R);

        /* Save profile */
        if (outdir) {
            char fname[256];
            snprintf(fname, sizeof(fname), "%s/soler_omega%.3f.dat", outdir, par_omega);
            FILE *fp = fopen(fname, "w");
            if (fp) {
                fprintf(fp, "# Soler: omega=%.6f m=%.6f lambda=%.6f f0=%.8f\n",
                        par_omega, par_m, par_lambda, f0);
                fprintf(fp, "# r  f  g  rho  S\n");
                double h = rmax_k / N;
                int step = (N > 5000) ? N / 5000 : 1;
                for (int i = 0; i <= N; i += step) {
                    double r = i * h;
                    double f = f_arr[i], g = g_arr[i];
                    fprintf(fp, "%.6f %.10e %.10e %.10e %.10e\n",
                            r, f, g, f * f + g * g, f * f - g * g);
                }
                fclose(fp);
            }
        }
    }
}

static void usage(const char *prog)
{
    fprintf(stderr,
            "Usage: %s [-omega ω] [-m mass] [-lambda λ] [-rmax R] [-N n] "
            "[-scan] [-outdir dir]\n",
            prog);
}

int main(int argc, char *argv[])
{
    double omega = 0.9;
    double rmax = 50.0;
    int N = 50000;
    int scan = 0;
    const char *outdir = NULL;

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-omega") == 0 && i + 1 < argc)
            omega = atof(argv[++i]);
        else if (strcmp(argv[i], "-m") == 0 && i + 1 < argc)
            par_m = atof(argv[++i]);
        else if (strcmp(argv[i], "-lambda") == 0 && i + 1 < argc)
            par_lambda = atof(argv[++i]);
        else if (strcmp(argv[i], "-rmax") == 0 && i + 1 < argc)
            rmax = atof(argv[++i]);
        else if (strcmp(argv[i], "-N") == 0 && i + 1 < argc)
            N = atoi(argv[++i]);
        else if (strcmp(argv[i], "-scan") == 0)
            scan = 1;
        else if (strcmp(argv[i], "-outdir") == 0 && i + 1 < argc)
            outdir = argv[++i];
        else {
            usage(argv[0]);
            return 1;
        }
    }

    if (N > NMAX - 1) N = NMAX - 1;

    if (scan) {
        run_scan(rmax, N, outdir);
    } else {
        run_single(omega, rmax, N, outdir);
    }

    return 0;
}
