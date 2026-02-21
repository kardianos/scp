/*
 * charged_soler.c — Charged Soler model (nonlinear Dirac + Coulomb)
 *
 * Extends soler.c with U(1) gauge coupling: Dirac + (ψ̄ψ)² + EM.
 *
 * Lagrangian: L = ψ̄(iγ^μ D_μ - m)ψ + λ(ψ̄ψ)² - (1/4)F_μν F^{μν}
 * where D_μ = ∂_μ - ieA_μ.
 *
 * System (flat space, spherically symmetric):
 *   f' = (m + (ω+eV) - λ(f²-g²)) g         (Dirac, upper)
 *   g' + 2g/r = (m - (ω+eV) - λ(f²-g²)) f  (Dirac, lower)
 *   V'' + 2V'/r = -4πe(f²+g²)               (Maxwell, Gauss's law)
 *
 * Method: SCF iteration
 *   1. Solve Soler with current V(r) → get f(r), g(r)
 *   2. Compute charge density ρ_ch = e(f²+g²)
 *   3. Solve Poisson for V via Green's function
 *   4. Under-relax V, repeat until convergence
 *
 * The EM coupling shifts ω → ω + eV(r), making the effective frequency
 * position-dependent. The Coulomb self-energy is repulsive (same-sign charge).
 *
 * Usage:
 *   charged_soler -omega 0.9 -e 0.1           (single case)
 *   charged_soler -omega 0.9 -scan_e           (scan charge coupling)
 */

#define _DEFAULT_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define NMAX 200001

static double f_arr[NMAX], g_arr[NMAX];
static double V_arr[NMAX], V_new[NMAX];

/* Parameters */
static double par_m = 1.0;
static double par_omega = 0.9;
static double par_lambda = 1.0;
static double par_e = 0.1;       /* electromagnetic coupling */

/* RHS of Soler system with position-dependent V(r) */
static void rhs(double r, double f, double g, int i,
                double *df, double *dg)
{
    double S = f * f - g * g;
    double omega_eff = par_omega + par_e * V_arr[i];

    *df = (par_m + omega_eff - par_lambda * S) * g;
    if (r < 1e-12) {
        *dg = (par_m - omega_eff - par_lambda * S) * f / 3.0;
    } else {
        *dg = -2.0 * g / r + (par_m - omega_eff - par_lambda * S) * f;
    }
}

/*
 * RK4 integration with position-dependent V(r).
 * V is interpolated at half-steps using V[i] and V[i+1].
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

        /* Save/restore V for midpoint: use linear interp */
        double V_save = V_arr[i];
        double V_mid = 0.5 * (V_arr[i] + V_arr[i + 1]);

        double k1f, k1g, k2f, k2g, k3f, k3g, k4f, k4g;

        /* k1: at (r, i) */
        rhs(r, f, g, i, &k1f, &k1g);

        /* k2, k3: at (r+h/2) — temporarily set V_arr[i] to midpoint */
        V_arr[i] = V_mid;
        rhs(r + 0.5 * h, f + 0.5 * h * k1f, g + 0.5 * h * k1g,
            i, &k2f, &k2g);
        rhs(r + 0.5 * h, f + 0.5 * h * k2f, g + 0.5 * h * k2g,
            i, &k3f, &k3g);
        V_arr[i] = V_save;  /* restore */

        /* k4: at (r+h, i+1) */
        rhs(r + h, f + h * k3f, g + h * k3g, i + 1, &k4f, &k4g);

        f_arr[i + 1] = f + h * (k1f + 2 * k2f + 2 * k3f + k4f) / 6.0;
        g_arr[i + 1] = g + h * (k1g + 2 * k2g + 2 * k3g + k4g) / 6.0;

        if (fabs(f_arr[i + 1]) > 1e10 || fabs(g_arr[i + 1]) > 1e10)
            return i + 1;
    }
    return N;
}

/* Classify solution: -1 if f crosses zero in physical region */
static int classify(int nv)
{
    double amp_min = f_arr[0] * f_arr[0] + g_arr[0] * g_arr[0];
    int imin = 0;
    for (int i = 1; i <= nv; i++) {
        double amp = f_arr[i] * f_arr[i] + g_arr[i] * g_arr[i];
        if (amp < amp_min) {
            amp_min = amp;
            imin = i;
        }
    }
    for (int i = 1; i <= imin; i++) {
        if (f_arr[i] < 0) return -1;
    }
    return +1;
}

/* Find f₀ by bisection */
static double find_f0(int N, double rmax)
{
    double last_pos = -1, first_neg = -1;
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

    double f_lo = last_pos;
    double f_hi = first_neg;

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

    /* Clip growing tail */
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
 * Solve Poisson for V(r) via Green's function:
 *   V'' + 2V'/r = -4πe(f²+g²)
 *
 * Solution: V(r) = (e/r)M_enc(r) + e·I_ext(r)
 * where M_enc(r) = 4π∫₀ʳ (f²+g²)r'²dr'  (enclosed charge)
 *       I_ext(r) = 4π∫ᵣ^∞ (f²+g²)r' dr'  (external integral)
 *
 * Note: V > 0 for positive charge (repulsive self-energy).
 */
static void solve_poisson_em(int N, double rmax)
{
    double h = rmax / N;
    static double M_enc[NMAX], I_ext[NMAX];

    /* Enclosed charge: M_enc(r) = 4π∫₀ʳ ρ_ch r'² dr' */
    M_enc[0] = 0;
    for (int i = 1; i <= N; i++) {
        double r = i * h, rp = (i - 1) * h;
        double rho_p = f_arr[i - 1] * f_arr[i - 1] + g_arr[i - 1] * g_arr[i - 1];
        double rho_c = f_arr[i] * f_arr[i] + g_arr[i] * g_arr[i];
        M_enc[i] = M_enc[i - 1] + 0.5 * 4.0 * M_PI * (rho_p * rp * rp + rho_c * r * r) * h;
    }

    /* External integral: I_ext(r) = 4π∫ᵣ^∞ ρ_ch r' dr' */
    I_ext[N] = 0;
    for (int i = N - 1; i >= 0; i--) {
        double r = i * h, rn = (i + 1) * h;
        double rho_c = f_arr[i] * f_arr[i] + g_arr[i] * g_arr[i];
        double rho_n = f_arr[i + 1] * f_arr[i + 1] + g_arr[i + 1] * g_arr[i + 1];
        I_ext[i] = I_ext[i + 1] + 0.5 * 4.0 * M_PI * (rho_c * r + rho_n * rn) * h;
    }

    /* V(r) = (e/r)M_enc + e·I_ext */
    V_new[0] = par_e * I_ext[0];  /* V(0) = e·∫₀^∞ ρ r dr */
    for (int i = 1; i <= N; i++) {
        double r = i * h;
        V_new[i] = par_e * (M_enc[i] / r + I_ext[i]);
    }
}

/* Compute particle number Q and R_rms */
static void compute_props(int N, double rmax,
                          double *Q_out, double *R_out, double *E_coul_out)
{
    double h = rmax / N;
    double Q = 0, R2 = 0, E_coul = 0;

    for (int i = 0; i <= N; i++) {
        double r = i * h;
        double rho = f_arr[i] * f_arr[i] + g_arr[i] * g_arr[i];

        double w;
        if (i == 0 || i == N)
            w = 1.0 / 3.0;
        else if (i % 2 == 1)
            w = 4.0 / 3.0;
        else
            w = 2.0 / 3.0;

        Q += w * rho * r * r * h;
        R2 += w * rho * r * r * r * r * h;
        /* Coulomb energy: E_C = (1/2)∫ ρ_ch V d³x */
        E_coul += w * rho * V_arr[i] * r * r * h;
    }

    *Q_out = 4.0 * M_PI * Q;
    *R_out = (Q > 0) ? sqrt(R2 / Q) : 0;
    *E_coul_out = 4.0 * M_PI * 0.5 * par_e * E_coul;
}

static void run_single(double omega, double e_charge, double rmax, int N,
                        const char *outdir)
{
    par_omega = omega;
    par_e = e_charge;
    double kappa = sqrt(par_m * par_m - par_omega * par_omega);

    if (10.0 / kappa > rmax) rmax = 10.0 / kappa;
    if (rmax > 500) rmax = 500;

    printf("# Charged Soler: omega=%.6f, m=%.6f, lambda=%.6f, e=%.6f\n",
           par_omega, par_m, par_lambda, par_e);
    printf("# rmax=%.1f, N=%d, h=%.6f, kappa=%.6f\n",
           rmax, N, rmax / N, kappa);

    /* Initialize V = 0 */
    for (int i = 0; i <= N; i++)
        V_arr[i] = 0;

    double f0_prev = 0;
    double relax = 0.2;
    int max_iter = 100;
    int converged = 0;

    printf("# %4s %14s %14s %14s %12s %12s\n",
           "iter", "f(0)", "V(0)", "E_Coulomb", "Q", "R_rms");

    for (int iter = 0; iter < max_iter; iter++) {
        double f0 = find_f0(N, rmax);
        if (f0 < 0) {
            printf("# FAILED at iter %d: no bracket found\n", iter);
            printf("# This may indicate the charge is too large for binding.\n");
            return;
        }

        /* Compute charge density and solve Poisson for V */
        solve_poisson_em(N, rmax);

        /* Under-relax */
        for (int i = 0; i <= N; i++)
            V_arr[i] = relax * V_new[i] + (1.0 - relax) * V_arr[i];

        double Q, R, E_coul;
        compute_props(N, rmax, &Q, &R, &E_coul);

        printf("  %4d %14.8f %14.8f %14.8f %12.6f %12.4f\n",
               iter, f0, V_arr[0], E_coul, Q, R);
        fflush(stdout);

        if (iter > 5 && fabs(f0 - f0_prev) < 1e-8 * fabs(f0)) {
            converged = 1;
            printf("# Converged at iter %d\n", iter);
            break;
        }
        f0_prev = f0;
    }

    if (!converged) {
        printf("# WARNING: did not converge in %d iterations\n", max_iter);
    }

    /* Final properties */
    double f0_final = find_f0(N, rmax);
    double Q, R, E_coul;
    compute_props(N, rmax, &Q, &R, &E_coul);

    double E_bind = par_omega - par_m;
    double M = par_omega * Q;

    printf("\n# RESULT: Charged Soler soliton\n");
    printf("#   f(0)       = %.8f\n", f0_final);
    printf("#   omega      = %.8f\n", par_omega);
    printf("#   e (charge) = %.8f\n", par_e);
    printf("#   E_bind     = %.8f  (omega - m)\n", E_bind);
    printf("#   Q          = %.6f  (particle number)\n", Q);
    printf("#   M=omega*Q  = %.6f\n", M);
    printf("#   R_rms      = %.6f\n", R);
    printf("#   V(0)       = %.8f  (central Coulomb potential)\n", V_arr[0]);
    printf("#   E_Coulomb  = %.8f  (Coulomb self-energy)\n", E_coul);
    printf("#   E_C/M      = %.6f  (Coulomb fraction)\n",
           M > 0 ? E_coul / M : 0);

    if (outdir) {
        char fname[256];
        snprintf(fname, sizeof(fname), "%s/charged_soler_omega%.3f_e%.4f.dat",
                 outdir, par_omega, par_e);
        FILE *fp = fopen(fname, "w");
        if (fp) {
            fprintf(fp, "# Charged Soler: omega=%.6f m=%.6f lambda=%.6f e=%.6f\n",
                    par_omega, par_m, par_lambda, par_e);
            fprintf(fp, "# f0=%.8f Q=%.6f R=%.6f V(0)=%.8f E_C=%.8f\n",
                    f0_final, Q, R, V_arr[0], E_coul);
            fprintf(fp, "# r  f  g  rho  S  V\n");
            double h = rmax / N;
            int step = (N > 5000) ? N / 5000 : 1;
            for (int i = 0; i <= N; i += step) {
                double r = i * h;
                double f = f_arr[i], g = g_arr[i];
                fprintf(fp, "%.6f %.10e %.10e %.10e %.10e %.10e\n",
                        r, f, g, f * f + g * g, f * f - g * g, V_arr[i]);
            }
            fclose(fp);
            printf("#   Profile: %s\n", fname);
        }
    }
}

static void run_scan_e(double omega, double rmax, int N, const char *outdir)
{
    double e_vals[] = {0.0, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0};
    int nvals = sizeof(e_vals) / sizeof(e_vals[0]);

    printf("# Charged Soler scan: omega=%.6f, m=%.6f, lambda=%.6f\n",
           omega, par_m, par_lambda);
    printf("# %8s  %12s  %12s  %12s  %12s  %8s  %8s\n",
           "e", "f(0)", "V(0)", "E_Coulomb", "Q", "R_rms", "E_C/M");

    for (int k = 0; k < nvals; k++) {
        par_e = e_vals[k];
        par_omega = omega;
        double rmax_k = rmax;
        double kappa = sqrt(par_m * par_m - par_omega * par_omega);
        if (kappa > 0 && 10.0 / kappa > rmax_k)
            rmax_k = 10.0 / kappa;
        if (rmax_k > 500) rmax_k = 500;

        /* Initialize V = 0 */
        for (int i = 0; i <= N; i++)
            V_arr[i] = 0;

        double f0_prev = 0;
        int converged = 0;
        double f0 = -1;

        for (int iter = 0; iter < 100; iter++) {
            f0 = find_f0(N, rmax_k);
            if (f0 < 0) break;

            solve_poisson_em(N, rmax_k);
            for (int i = 0; i <= N; i++)
                V_arr[i] = 0.2 * V_new[i] + 0.8 * V_arr[i];

            if (iter > 5 && fabs(f0 - f0_prev) < 1e-8 * fabs(f0)) {
                converged = 1;
                break;
            }
            f0_prev = f0;
        }

        if (f0 < 0 || !converged) {
            printf("  %8.4f  FAILED (charge too large or no convergence)\n", par_e);
            continue;
        }

        double Q, R, E_coul;
        compute_props(N, rmax_k, &Q, &R, &E_coul);
        double M = par_omega * Q;

        printf("  %8.4f  %12.6f  %12.8f  %12.8f  %12.6f  %8.4f  %8.5f\n",
               par_e, f0, V_arr[0], E_coul, Q, R,
               M > 0 ? E_coul / M : 0);

        if (outdir) {
            char fname[256];
            snprintf(fname, sizeof(fname), "%s/charged_soler_omega%.3f_e%.4f.dat",
                     outdir, omega, par_e);
            FILE *fp = fopen(fname, "w");
            if (fp) {
                fprintf(fp, "# r  f  g  rho  S  V\n");
                double h = rmax_k / N;
                int step = (N > 5000) ? N / 5000 : 1;
                for (int i = 0; i <= N; i += step) {
                    double r = i * h;
                    fprintf(fp, "%.6f %.10e %.10e %.10e %.10e %.10e\n",
                            r, f_arr[i], g_arr[i],
                            f_arr[i] * f_arr[i] + g_arr[i] * g_arr[i],
                            f_arr[i] * f_arr[i] - g_arr[i] * g_arr[i],
                            V_arr[i]);
                }
                fclose(fp);
            }
        }
    }
}

static void usage(const char *prog)
{
    fprintf(stderr,
            "Usage: %s [-omega ω] [-m mass] [-lambda λ] [-e charge] "
            "[-rmax R] [-N n] [-scan_e] [-outdir dir]\n",
            prog);
}

int main(int argc, char *argv[])
{
    double omega = 0.9;
    double e_charge = 0.1;
    double rmax = 50.0;
    int N = 50000;
    int scan_e = 0;
    const char *outdir = NULL;

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-omega") == 0 && i + 1 < argc)
            omega = atof(argv[++i]);
        else if (strcmp(argv[i], "-m") == 0 && i + 1 < argc)
            par_m = atof(argv[++i]);
        else if (strcmp(argv[i], "-lambda") == 0 && i + 1 < argc)
            par_lambda = atof(argv[++i]);
        else if (strcmp(argv[i], "-e") == 0 && i + 1 < argc)
            e_charge = atof(argv[++i]);
        else if (strcmp(argv[i], "-rmax") == 0 && i + 1 < argc)
            rmax = atof(argv[++i]);
        else if (strcmp(argv[i], "-N") == 0 && i + 1 < argc)
            N = atoi(argv[++i]);
        else if (strcmp(argv[i], "-scan_e") == 0)
            scan_e = 1;
        else if (strcmp(argv[i], "-outdir") == 0 && i + 1 < argc)
            outdir = argv[++i];
        else {
            usage(argv[0]);
            return 1;
        }
    }

    if (N > NMAX - 1) N = NMAX - 1;

    if (scan_e) {
        run_scan_e(omega, rmax, N, outdir);
    } else {
        run_single(omega, e_charge, rmax, N, outdir);
    }

    return 0;
}
