/*
 * soler.c — Soler model soliton solver for v5 parameter space exploration
 *
 * Lagrangian: L = ψ̄(iγ∂ - m)ψ + (λ/2)(ψ̄ψ)²
 *
 * Radial equations (κ = -1, ground state):
 *   f' = (m + ω - λS) g
 *   g' + 2g/r = (m - ω - λS) f
 * where S = f² - g² (scalar density).
 *
 * Total mass: M = ωQ + E_self, where
 *   Q = 4π ∫ (f²+g²) r² dr       (fermion number)
 *   E_self = 2πλ ∫ S² r² dr       (self-interaction energy)
 *
 * Physical nucleon constraint: M * R / (ℏc) = 4.00
 * (In code units with m=1: M_code * R_code = 4.00)
 *
 * Usage:
 *   soler -omega 0.9                     (single solution)
 *   soler -omega 0.9 -lambda 5.0         (single, different λ)
 *   soler -scan                           (ω scan at fixed λ)
 *   soler -scan2d                         (2D scan over λ × ω)
 *   soler -findMR 4.0                     (find ω where M*R = target, at fixed λ)
 */

#define _DEFAULT_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define NMAX 200001

static double f_arr[NMAX], g_arr[NMAX];

/* Parameters */
static double par_m = 1.0;
static double par_omega = 0.9;
static double par_lambda = 1.0;

/* RHS of Soler radial system */
static void rhs(double r, double f, double g, double *df, double *dg)
{
    double S = f * f - g * g;
    *df = (par_m + par_omega - par_lambda * S) * g;
    if (r < 1e-12)
        *dg = (par_m - par_omega - par_lambda * S) * f / 3.0;
    else
        *dg = -2.0 * g / r + (par_m - par_omega - par_lambda * S) * f;
}

/* RK4 integration. Returns last valid index. */
static int integrate(double f0, int N, double rmax)
{
    double h = rmax / N;
    f_arr[0] = f0;
    g_arr[0] = 0;

    for (int i = 0; i < N; i++) {
        double r = i * h;
        double f = f_arr[i], g = g_arr[i];
        double k1f, k1g, k2f, k2g, k3f, k3g, k4f, k4g;

        rhs(r, f, g, &k1f, &k1g);
        rhs(r + 0.5*h, f + 0.5*h*k1f, g + 0.5*h*k1g, &k2f, &k2g);
        rhs(r + 0.5*h, f + 0.5*h*k2f, g + 0.5*h*k2g, &k3f, &k3g);
        rhs(r + h, f + h*k3f, g + h*k3g, &k4f, &k4g);

        f_arr[i+1] = f + h*(k1f + 2*k2f + 2*k3f + k4f)/6.0;
        g_arr[i+1] = g + h*(k1g + 2*k2g + 2*k3g + k4g)/6.0;

        if (fabs(f_arr[i+1]) > 1e10 || fabs(g_arr[i+1]) > 1e10)
            return i + 1;
    }
    return N;
}

/* Classify: -1 if f crosses zero (overshoot), +1 if f stays positive */
static int classify(int nv)
{
    double amp_min = f_arr[0]*f_arr[0] + g_arr[0]*g_arr[0];
    int imin = 0;
    for (int i = 1; i <= nv; i++) {
        double amp = f_arr[i]*f_arr[i] + g_arr[i]*g_arr[i];
        if (amp < amp_min) { amp_min = amp; imin = i; }
    }
    for (int i = 1; i <= imin; i++)
        if (f_arr[i] < 0) return -1;
    return +1;
}

/* Bisection on f₀. Returns f₀ or -1 on failure. */
static double find_f0(int N, double rmax)
{
    double last_pos = -1, first_neg = -1;
    int found_pos = 0;

    for (double f_scan = 0.001; f_scan <= 30.0; f_scan *= 1.01) {
        int nv = integrate(f_scan, N, rmax);
        int s = classify(nv);
        if (s == +1) { last_pos = f_scan; found_pos = 1; }
        else if (s == -1 && found_pos) { first_neg = f_scan; break; }
    }
    if (last_pos < 0 || first_neg < 0) return -1;

    double f_lo = last_pos, f_hi = first_neg;
    for (int iter = 0; iter < 80; iter++) {
        double f_mid = 0.5*(f_lo + f_hi);
        int nv = integrate(f_mid, N, rmax);
        int s = classify(nv);
        if (s == +1) f_lo = f_mid; else f_hi = f_mid;
    }

    double f0 = 0.5*(f_lo + f_hi);
    integrate(f0, N, rmax);

    /* Clip growing tail */
    double amp_min = f_arr[0]*f_arr[0] + g_arr[0]*g_arr[0];
    int imin = 0;
    for (int i = 1; i <= N; i++) {
        double amp = f_arr[i]*f_arr[i] + g_arr[i]*g_arr[i];
        if (amp < amp_min) { amp_min = amp; imin = i; }
    }
    for (int i = imin+1; i <= N; i++) { f_arr[i] = 0; g_arr[i] = 0; }
    return f0;
}

/*
 * Compute observables: Q, R_rms, E_self, M_total = ωQ + E_self.
 */
static void compute_props(int N, double rmax,
    double *Q_out, double *R_out, double *Eself_out, double *M_out)
{
    double h = rmax / N;
    double Q = 0, R2 = 0, ES = 0;

    for (int i = 0; i <= N; i++) {
        double r = i * h;
        double fi = f_arr[i], gi = g_arr[i];
        double rho = fi*fi + gi*gi;          /* number density */
        double S = fi*fi - gi*gi;            /* scalar density */

        double w;
        if (i == 0 || i == N) w = 1.0/3.0;
        else if (i % 2 == 1) w = 4.0/3.0;
        else w = 2.0/3.0;

        Q  += w * rho * r*r * h;
        R2 += w * rho * r*r*r*r * h;
        ES += w * S*S * r*r * h;
    }

    *Q_out = 4.0 * M_PI * Q;
    *R_out = (Q > 0) ? sqrt(R2 / Q) : 0;
    *Eself_out = 2.0 * M_PI * par_lambda * ES;
    *M_out = par_omega * (*Q_out) + (*Eself_out);
}

/*
 * Solve one point and return success (1) or failure (0).
 * Fills in f0, Q, R, E_self, M.
 */
static int solve_one(double omega, double lambda, int N, double rmax,
    double *f0_out, double *Q_out, double *R_out, double *Eself_out, double *M_out)
{
    par_omega = omega;
    par_lambda = lambda;

    double kappa = sqrt(par_m*par_m - omega*omega);
    double rmax_use = rmax;
    if (10.0/kappa > rmax_use) rmax_use = 10.0/kappa;
    if (rmax_use > 500) rmax_use = 500;

    double f0 = find_f0(N, rmax_use);
    if (f0 < 0) return 0;

    compute_props(N, rmax_use, Q_out, R_out, Eself_out, M_out);
    *f0_out = f0;
    return 1;
}

/* ---- Run modes ---- */

static void run_single(double omega, double lambda, int N, double rmax,
                        const char *outdir)
{
    double f0, Q, R, Eself, M;
    if (!solve_one(omega, lambda, N, rmax, &f0, &Q, &R, &Eself, &M)) {
        printf("# FAILED: No solution at omega=%.6f lambda=%.6f\n", omega, lambda);
        return;
    }

    double kappa = sqrt(par_m*par_m - omega*omega);
    double MR = M * R;

    printf("# Soler soliton: m=%.6f lambda=%.6f omega=%.6f\n", par_m, lambda, omega);
    printf("#   f(0)      = %.8f\n", f0);
    printf("#   Q         = %.6f  (fermion number)\n", Q);
    printf("#   omega*Q   = %.6f  (kinetic mass)\n", omega*Q);
    printf("#   E_self    = %.6f  (self-interaction energy)\n", Eself);
    printf("#   M_total   = %.6f  (omega*Q + E_self)\n", M);
    printf("#   R_rms     = %.6f\n", R);
    printf("#   M*R       = %.6f  (target: 4.00)\n", MR);
    printf("#   kappa     = %.6f  (decay rate)\n", kappa);
    printf("#   E_self/M  = %.4f  (self-interaction fraction)\n", Eself/M);
    printf("#   lambda*m² = %.6f  (dimensionless coupling)\n", lambda*par_m*par_m);

    if (outdir) {
        char fname[256];
        snprintf(fname, sizeof(fname), "%s/soler_lam%.3f_om%.3f.dat",
                 outdir, lambda, omega);
        double rmax_use = rmax;
        if (10.0/kappa > rmax_use) rmax_use = 10.0/kappa;
        if (rmax_use > 500) rmax_use = 500;
        double h = rmax_use / N;
        FILE *fp = fopen(fname, "w");
        if (fp) {
            fprintf(fp, "# Soler: omega=%.6f m=%.6f lambda=%.6f f0=%.8f M=%.6f R=%.6f\n",
                    omega, par_m, lambda, f0, M, R);
            fprintf(fp, "# r  f  g  rho=f²+g²  S=f²-g²\n");
            int step = (N > 5000) ? N/5000 : 1;
            for (int i = 0; i <= N; i += step) {
                double r = i * h;
                fprintf(fp, "%.6f %.10e %.10e %.10e %.10e\n",
                        r, f_arr[i], g_arr[i],
                        f_arr[i]*f_arr[i] + g_arr[i]*g_arr[i],
                        f_arr[i]*f_arr[i] - g_arr[i]*g_arr[i]);
            }
            fclose(fp);
            printf("#   Profile: %s\n", fname);
        }
    }
}

/* Omega scan at fixed lambda */
static void run_scan(double lambda, int N, double rmax)
{
    printf("# Soler omega scan: m=%.6f lambda=%.6f (lambda*m²=%.6f)\n",
           par_m, lambda, lambda*par_m*par_m);
    printf("# %8s %10s %10s %10s %10s %10s %10s %8s\n",
           "omega", "f0", "Q", "omega*Q", "E_self", "M_total", "R_rms", "M*R");

    /* Dense scan from omega_min to near m */
    for (double omega = 0.10; omega <= 0.995; omega += 0.01) {
        if (omega >= par_m) break;
        double f0, Q, R, Eself, M;
        if (!solve_one(omega, lambda, N, rmax, &f0, &Q, &R, &Eself, &M))
            continue;
        printf("  %8.4f %10.5f %10.3f %10.3f %10.3f %10.3f %10.4f %8.3f\n",
               omega, f0, Q, omega*Q, Eself, M, R, M*R);
    }
}

/* 2D scan over lambda × omega */
static void run_scan2d(int N, double rmax)
{
    double lambdas[] = {0.01, 0.02, 0.05, 0.1, 0.2, 0.5,
                        1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0};
    int nlam = sizeof(lambdas)/sizeof(lambdas[0]);

    printf("# Soler 2D scan: m=%.6f\n", par_m);
    printf("# %8s %8s %10s %10s %10s %10s %10s %8s %8s\n",
           "lambda", "omega", "f0", "Q", "E_self", "M_total", "R_rms", "M*R", "Es/M");

    for (int il = 0; il < nlam; il++) {
        double lam = lambdas[il];
        int found_any = 0;
        double MR_min = 1e30, Q_at_MRmin = 0, om_at_MRmin = 0;
        double M_at_MRmin = 0, R_at_MRmin = 0;

        for (double omega = 0.10; omega <= 0.995; omega += 0.01) {
            if (omega >= par_m) break;
            double f0, Q, R, Eself, M;
            if (!solve_one(omega, lam, N, rmax, &f0, &Q, &R, &Eself, &M))
                continue;
            found_any = 1;

            double MR = M * R;
            printf("  %8.4f %8.4f %10.5f %10.3f %10.3f %10.3f %10.4f %8.3f %8.4f\n",
                   lam, omega, f0, Q, Eself, M, R, MR, Eself/M);

            if (MR < MR_min) {
                MR_min = MR;
                Q_at_MRmin = Q;
                om_at_MRmin = omega;
                M_at_MRmin = M;
                R_at_MRmin = R;
            }
        }

        if (found_any) {
            fprintf(stderr, "lambda=%8.4f: M*R_min=%8.3f at omega=%.3f "
                    "(Q=%.1f M=%.3f R=%.4f)\n",
                    lam, MR_min, om_at_MRmin, Q_at_MRmin,
                    M_at_MRmin, R_at_MRmin);
        }
    }
}

/*
 * Find omega where M*R = target, at fixed lambda.
 * Bisects on omega.
 */
static void run_findMR(double target, double lambda, int N, double rmax)
{
    printf("# Finding omega where M*R = %.3f at lambda=%.6f\n", target, lambda);

    /* First scan to find bracket: M*R decreases then increases with omega,
     * so we need two brackets (if the minimum M*R is below target). */
    double omegas[100], MRs[100];
    int npts = 0;

    for (double omega = 0.10; omega <= 0.995 && npts < 100; omega += 0.01) {
        if (omega >= par_m) break;
        double f0, Q, R, Eself, M;
        if (!solve_one(omega, lambda, N, rmax, &f0, &Q, &R, &Eself, &M))
            continue;
        omegas[npts] = omega;
        MRs[npts] = M * R;
        npts++;
    }

    if (npts < 2) {
        printf("# FAILED: Not enough solution points\n");
        return;
    }

    /* Find minimum M*R */
    double MR_min = MRs[0];
    int imin = 0;
    for (int i = 1; i < npts; i++) {
        if (MRs[i] < MR_min) { MR_min = MRs[i]; imin = i; }
    }

    printf("# M*R minimum = %.3f at omega=%.4f\n", MR_min, omegas[imin]);

    if (MR_min > target) {
        printf("# FAILED: M*R_min=%.3f > target=%.3f. No solution exists.\n",
               MR_min, target);
        return;
    }

    /* Find brackets: one on the lower-omega side, one on the upper */
    for (int side = 0; side < 2; side++) {
        int i_lo = -1, i_hi = -1;
        if (side == 0) {
            /* Lower omega side: scan from 0 to imin */
            for (int i = 0; i < imin; i++) {
                if (MRs[i] > target && MRs[i+1] <= target)
                    { i_lo = i; i_hi = i+1; break; }
            }
        } else {
            /* Upper omega side: scan from imin to end */
            for (int i = imin; i < npts-1; i++) {
                if (MRs[i] <= target && MRs[i+1] > target)
                    { i_lo = i; i_hi = i+1; break; }
            }
        }

        if (i_lo < 0) continue;

        /* Bisect */
        double om_lo = omegas[i_lo], om_hi = omegas[i_hi];
        double MR_lo = MRs[i_lo], MR_hi = MRs[i_hi];

        for (int iter = 0; iter < 50; iter++) {
            double om_mid = 0.5*(om_lo + om_hi);
            double f0, Q, R, Eself, M;
            if (!solve_one(om_mid, lambda, N, rmax, &f0, &Q, &R, &Eself, &M)) {
                /* If solve fails, shrink toward the working side */
                if (MR_lo <= target) om_hi = om_mid;
                else om_lo = om_mid;
                continue;
            }
            double MR = M * R;
            if ((MR_lo > target && MR > target) || (MR_lo <= target && MR <= target)) {
                om_lo = om_mid; MR_lo = MR;
            } else {
                om_hi = om_mid; MR_hi = MR;
            }
        }

        /* Final solution at midpoint */
        double om_final = 0.5*(om_lo + om_hi);
        double f0, Q, R, Eself, M_total;
        if (solve_one(om_final, lambda, N, rmax, &f0, &Q, &R, &Eself, &M_total)) {
            printf("# SOLUTION %s: omega=%.8f\n", side==0 ? "LOW" : "HIGH", om_final);
            printf("#   f(0)    = %.8f\n", f0);
            printf("#   Q       = %.6f\n", Q);
            printf("#   M_total = %.6f  (omega*Q=%.6f + E_self=%.6f)\n",
                   M_total, om_final*Q, Eself);
            printf("#   R_rms   = %.6f\n", R);
            printf("#   M*R     = %.6f  (target=%.3f)\n", M_total*R, target);
            printf("#   E_self/M = %.4f\n", Eself/M_total);
        }
    }
}

static void usage(const char *prog)
{
    fprintf(stderr,
        "Usage: %s [options]\n"
        "  -omega ω       frequency (default 0.9)\n"
        "  -m mass         fermion mass (default 1.0)\n"
        "  -lambda λ       self-coupling (default 1.0)\n"
        "  -rmax R         integration domain (default 50)\n"
        "  -N n            grid points (default 50000)\n"
        "  -scan           omega scan at fixed lambda\n"
        "  -scan2d         2D scan over lambda × omega\n"
        "  -findMR target  find omega where M*R = target\n"
        "  -outdir dir     output directory for profiles\n",
        prog);
}

int main(int argc, char *argv[])
{
    double omega = 0.9, lambda = 1.0, rmax = 50.0;
    int N = 50000;
    int mode_scan = 0, mode_scan2d = 0, mode_findMR = 0;
    double MR_target = 4.0;
    const char *outdir = NULL;

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-omega") == 0 && i+1 < argc)
            omega = atof(argv[++i]);
        else if (strcmp(argv[i], "-m") == 0 && i+1 < argc)
            par_m = atof(argv[++i]);
        else if (strcmp(argv[i], "-lambda") == 0 && i+1 < argc)
            lambda = atof(argv[++i]);
        else if (strcmp(argv[i], "-rmax") == 0 && i+1 < argc)
            rmax = atof(argv[++i]);
        else if (strcmp(argv[i], "-N") == 0 && i+1 < argc)
            N = atoi(argv[++i]);
        else if (strcmp(argv[i], "-scan") == 0)
            mode_scan = 1;
        else if (strcmp(argv[i], "-scan2d") == 0)
            mode_scan2d = 1;
        else if (strcmp(argv[i], "-findMR") == 0 && i+1 < argc) {
            mode_findMR = 1;
            MR_target = atof(argv[++i]);
        }
        else if (strcmp(argv[i], "-outdir") == 0 && i+1 < argc)
            outdir = argv[++i];
        else { usage(argv[0]); return 1; }
    }

    if (N > NMAX-1) N = NMAX-1;

    if (mode_scan2d)
        run_scan2d(N, rmax);
    else if (mode_scan)
        run_scan(lambda, N, rmax);
    else if (mode_findMR)
        run_findMR(MR_target, lambda, N, rmax);
    else
        run_single(omega, lambda, N, rmax, outdir);

    return 0;
}
