/*
 * strong_geon.c — Salam Strong Geon: Yukawa Hartree Solver
 *
 * Self-trapping of EM standing wave via massive spin-2 mediator (f_2).
 * Extends V8 geon.c by replacing the long-range Einstein mass integration
 * with a short-range Yukawa BVP:
 *
 *   V8 (GR):    m'(r) = kappa r^2 rho         (1/r potential)
 *   V9 (Yukawa): M'' - mu^2 M = -kappa r rho  (e^{-mu r}/r potential)
 *
 * The metric, EM mode, and Hartree iteration are identical to V8.
 *
 * Physical motivation: f_2(1270) spin-2 meson mediates a tensor force
 * with G_strong ~ 10^{38} G_N at range 1/m_f2 ~ 0.15 fm.
 *
 * Build: cd src && make
 * Usage: strong_geon -perturb -mu <M> [options]
 *        strong_geon -hartree -kappa <K> -mu <M> [options]
 */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define NMAX 16001

/* ---- Parameters ---- */
static double kappa   = 1.0;       /* coupling constant */
static double mu      = 1.0;       /* Yukawa mass (1/range) */
static double b_BI    = 0.0;       /* Born-Infeld field strength (0 = Maxwell) */
static int    ell     = 1;         /* angular momentum (>= 1) */
static double seed_amp = 0.0;      /* seed well amplitude (0 = flat start) */
static int    nl_metric = 0;       /* 0=linear, 1=exponential, 2=Padé */

/* ---- Multi-mode ---- */
#define MAX_MODES 4
static int n_multi_modes = 0;
static int mode_l_arr[MAX_MODES];
static double mode_occ_arr[MAX_MODES];

static int    Ngrid   = 4001;
static double Rmax    = 10.0;

/* ---- Grid ---- */
static double dr;
static double rg[NMAX];

/* ---- Metric ---- */
static double M_eff[NMAX];    /* effective mass function: f = 1 - 2M/r */
static double ff[NMAX];       /* f(r) = 1 - 2M/r */
static double ffp[NMAX];      /* f'(r) */
static double ffpp[NMAX];     /* f''(r) */

/* ---- EM mode ---- */
static double u_em[NMAX];     /* radial function u(r) */
static double v_em[NMAX];     /* symmetrized: v = sqrt(f) * u */

/* ==== Grid ==== */
static void setup_grid(void)
{
    if (Ngrid > NMAX) Ngrid = NMAX;
    dr = Rmax / (Ngrid - 1);
    for (int i = 0; i < Ngrid; i++)
        rg[i] = (i + 0.5) * dr;  /* half-cell offset avoids r=0 */
}

/* ==== Metric from mass function ==== */
/*
 * nl_metric selects the nonlinear completion:
 *   0: linear   f = 1 - 2M/r  (floor at 0.01, has horizon)
 *   1: exp      f = exp(-2M/r) (always positive, no horizon)
 *   2: Padé     f = 1/(1+2M/r) (always positive, no horizon)
 * All agree to O(M/r) for weak fields.
 */
static void compute_metric(void)
{
    int n = Ngrid;

    for (int i = 0; i < n; i++) {
        double x = 2.0 * M_eff[i] / rg[i];  /* 2M/r */
        switch (nl_metric) {
        case 1:  /* exponential */
            ff[i] = exp(-x);
            break;
        case 2:  /* Padé */
            ff[i] = 1.0 / (1.0 + x);
            break;
        default: /* linear with floor */
            ff[i] = 1.0 - x;
            if (ff[i] < 0.01) ff[i] = 0.01;
            break;
        }
    }

    /* f' by centered differences */
    ffp[0] = (ff[1] - ff[0]) / dr;
    for (int i = 1; i < n - 1; i++)
        ffp[i] = (ff[i+1] - ff[i-1]) / (2.0 * dr);
    ffp[n-1] = (ff[n-1] - ff[n-2]) / dr;

    /* f'' */
    ffpp[0] = (ff[1] - 2.0*ff[0] + ff[0]) / (dr*dr);
    for (int i = 1; i < n - 1; i++)
        ffpp[i] = (ff[i+1] - 2.0*ff[i] + ff[i-1]) / (dr*dr);
    ffpp[n-1] = (ff[n-1] - 2.0*ff[n-2] + ff[n-3]) / (dr*dr);
}

/* ==== Effective potential W(r) for symmetrized variable v ==== */
static double W_eff_pot(int i)
{
    double r = rg[i];
    double f = ff[i];
    double fp = ffp[i];
    double fpp = ffpp[i];
    return f * (double)(ell*(ell+1)) / (r*r) + f*fpp/2.0 - fp*fp/4.0;
}

/* ==== Eigenvalue solver: implicit imaginary time ==== */
static double solve_eigenvalue(void)
{
    int n = Ngrid;
    double At[NMAX], Bt[NMAX], Ct[NMAX], Dt[NMAX];
    double v_new[NMAX];
    double dr2 = dr * dr;

    /* Initial guess: r^l * exp(-r^2/sigma^2) */
    double norm2 = 0;
    for (int i = 0; i < n; i++) {
        v_em[i] = pow(rg[i], ell) * exp(-rg[i]*rg[i] / (Rmax*Rmax/4.0));
        norm2 += v_em[i] * v_em[i] * dr;
    }
    double inv = 1.0 / sqrt(norm2);
    for (int i = 0; i < n; i++) v_em[i] *= inv;

    double dtau = 0.3;
    int niter = 8000;

    for (int iter = 0; iter < niter; iter++) {
        /* i=0: ghost v[-1] = -v[0] (odd parity) */
        double f2_0 = ff[0] * ff[0];
        At[0] = 0.0;
        Bt[0] = 1.0 + dtau * (3.0 * f2_0 / dr2 + W_eff_pot(0));
        Ct[0] = -dtau * f2_0 / dr2;
        Dt[0] = v_em[0];

        for (int i = 1; i < n - 1; i++) {
            double f2 = ff[i] * ff[i];
            At[i] = -dtau * f2 / dr2;
            Bt[i] = 1.0 + dtau * (2.0 * f2 / dr2 + W_eff_pot(i));
            Ct[i] = -dtau * f2 / dr2;
            Dt[i] = v_em[i];
        }

        /* Dirichlet at Rmax: v(Rmax) = 0 */
        At[n-1] = 0.0; Bt[n-1] = 1.0; Ct[n-1] = 0.0; Dt[n-1] = 0.0;

        /* Thomas algorithm */
        double cp[NMAX], dp[NMAX];
        cp[0] = Ct[0] / Bt[0];
        dp[0] = Dt[0] / Bt[0];
        for (int i = 1; i < n; i++) {
            double den = Bt[i] - At[i] * cp[i-1];
            cp[i] = Ct[i] / den;
            dp[i] = (Dt[i] - At[i] * dp[i-1]) / den;
        }
        v_new[n-1] = dp[n-1];
        for (int i = n - 2; i >= 0; i--)
            v_new[i] = dp[i] - cp[i] * v_new[i+1];

        /* Normalize */
        norm2 = 0.0;
        for (int i = 0; i < n; i++) norm2 += v_new[i] * v_new[i] * dr;
        if (norm2 < 1e-30) break;
        inv = 1.0 / sqrt(norm2);
        for (int i = 0; i < n; i++) v_em[i] = v_new[i] * inv;
    }

    /* Compute eigenvalue: omega^2 = <v|H|v>/<v|v> */
    double hv = 0.0, vv = 0.0;
    {
        double f2 = ff[0] * ff[0];
        double vdd = (v_em[1] - 3.0*v_em[0]) / dr2;  /* ghost */
        double Hv = -f2 * vdd + W_eff_pot(0) * v_em[0];
        hv += v_em[0] * Hv * dr;
        vv += v_em[0] * v_em[0] * dr;
    }
    for (int i = 1; i < n - 1; i++) {
        double f2 = ff[i] * ff[i];
        double vdd = (v_em[i+1] - 2.0*v_em[i] + v_em[i-1]) / dr2;
        double Hv = -f2 * vdd + W_eff_pot(i) * v_em[i];
        hv += v_em[i] * Hv * dr;
        vv += v_em[i] * v_em[i] * dr;
    }

    return hv / vv;
}

/* ==== Reconstruct u from v ==== */
static void reconstruct_u(void)
{
    for (int i = 0; i < Ngrid; i++) {
        double sf = sqrt(ff[i]);
        u_em[i] = (sf > 1e-10) ? v_em[i] / sf : 0.0;
    }
}

/* ==== Compute EM energy density ==== */
static void compute_rho(double omega2, double *rho)
{
    int n = Ngrid;
    double up[NMAX];
    double N_l = (double)(ell*(ell+1)) / (double)(2*ell+1);

    /* u' */
    up[0] = (u_em[1] + u_em[0]) / (2.0 * dr);
    for (int i = 1; i < n - 1; i++)
        up[i] = (u_em[i+1] - u_em[i-1]) / (2.0 * dr);
    up[n-1] = (u_em[n-1] - u_em[n-2]) / dr;

    for (int i = 0; i < n; i++) {
        double r = rg[i];
        double f = ff[i];
        double u = u_em[i];

        double rho_M = N_l / (8.0 * M_PI * r * r) *
                       (omega2 * u * u / f + f * up[i] * up[i]);

        if (b_BI > 0.0 && rho_M > b_BI * b_BI) {
            rho[i] = b_BI * b_BI * (sqrt(1.0 + rho_M / (b_BI*b_BI)) - 1.0);
        } else {
            rho[i] = rho_M;
        }
    }
}

/* ==== Solve Yukawa BVP: M'' - mu^2 M = -kappa r rho ==== */
/*
 * Thomas algorithm for tridiagonal system on the half-cell grid.
 *
 * Interior: M_{i-1} - (2+mu^2*h^2) M_i + M_{i+1} = -h^2 kappa r_i rho_i
 * i=0:   -(3+mu^2*h^2) M_0 + M_1 = -h^2 kappa r_0 rho_0  (antisymmetry ghost)
 * i=N-1: M_{N-2} + (exp(-mu*h)-2-mu^2*h^2) M_{N-1} = D  (Robin/Yukawa BC)
 *
 * Result: M_eff > 0 everywhere (attractive well), f = 1 - 2M/r < 1.
 */
static void solve_yukawa_bvp(double kap, double mu_val, double *rho, double *M_out)
{
    int n = Ngrid;
    double h = dr;
    double h2 = h * h;
    double mu2h2 = mu_val * mu_val * h2;
    double emh = exp(-mu_val * h);

    double At[NMAX], Bt[NMAX], Ct[NMAX], Dt[NMAX];

    /* i=0: ghost M_{-1} = -M_0 (antisymmetry) */
    At[0] = 0.0;
    Bt[0] = -(3.0 + mu2h2);
    Ct[0] = 1.0;
    Dt[0] = -h2 * kap * rg[0] * rho[0];

    /* Interior points */
    for (int i = 1; i < n - 1; i++) {
        At[i] = 1.0;
        Bt[i] = -(2.0 + mu2h2);
        Ct[i] = 1.0;
        Dt[i] = -h2 * kap * rg[i] * rho[i];
    }

    /* i=N-1: Robin BC, ghost M_N = M_{N-1} * exp(-mu*h) */
    At[n-1] = 1.0;
    Bt[n-1] = emh - 2.0 - mu2h2;
    Ct[n-1] = 0.0;
    Dt[n-1] = -h2 * kap * rg[n-1] * rho[n-1];

    /* Thomas forward sweep */
    double cp[NMAX], dp[NMAX];
    cp[0] = Ct[0] / Bt[0];
    dp[0] = Dt[0] / Bt[0];
    for (int i = 1; i < n; i++) {
        double den = Bt[i] - At[i] * cp[i-1];
        cp[i] = Ct[i] / den;
        dp[i] = (Dt[i] - At[i] * dp[i-1]) / den;
    }

    /* Back substitution */
    M_out[n-1] = dp[n-1];
    for (int i = n - 2; i >= 0; i--)
        M_out[i] = dp[i] - cp[i] * M_out[i+1];
}

/* ==== Integrate Einstein equation (V8 Poisson, for mu=0 comparison) ==== */
static void integrate_einstein(double kap, double *rho, double *m_new)
{
    int n = Ngrid;
    m_new[0] = 0.0;
    for (int i = 1; i < n; i++) {
        double r_mid = 0.5 * (rg[i] + rg[i-1]);
        double rho_mid = 0.5 * (rho[i] + rho[i-1]);
        m_new[i] = m_new[i-1] + kap * r_mid * r_mid * rho_mid * dr;
    }
}

/* ==== Compute potential from M_eff ==== */
static void compute_potential(double *Phi_out)
{
    for (int i = 0; i < Ngrid; i++)
        Phi_out[i] = -M_eff[i] / rg[i];
}

/* ==== First-order perturbation theory ==== */
static void run_perturbative(void)
{
    int n = Ngrid;

    printf("# V9 Strong Geon — First-Order Perturbation Theory\n");
    printf("# l=%d  mu=%.4f  N=%d  Rmax=%.1f\n", ell, mu, Ngrid, Rmax);
    printf("#\n");

    /* Step 1: Flat-space eigenmode */
    for (int i = 0; i < n; i++) M_eff[i] = 0.0;
    compute_metric();

    double omega2_0 = solve_eigenvalue();
    reconstruct_u();

    printf("# Flat-space eigenvalue: omega^2_0 = %.10f\n", omega2_0);

    /* Step 2: Compute source at unit kappa */
    double rho[NMAX], M_new[NMAX];
    compute_rho(omega2_0, rho);

    if (mu > 0.0) {
        solve_yukawa_bvp(1.0, mu, rho, M_new);
    } else {
        integrate_einstein(1.0, rho, M_new);
    }

    /* Diagnostics */
    double M_total = M_new[n-1];
    double M_peak = 0;
    int i_peak = 0;
    for (int i = 0; i < n; i++) {
        if (M_new[i] > M_peak) { M_peak = M_new[i]; i_peak = i; }
    }

    printf("# Mass function: M_peak = %.6e at r=%.4f, M(Rmax) = %.6e\n",
           M_peak, rg[i_peak], M_total);
    printf("# Yukawa range: 1/mu = %.4f code units\n", 1.0/mu);

    /* Step 3: First-order eigenvalue shift */
    /* delta_W = -2 M_new(r)/r * l(l+1)/r^2 (leading order) */
    double delta_omega2 = 0.0, vv = 0.0;
    for (int i = 0; i < n; i++) {
        double r = rg[i];
        double delta_W = -2.0 * M_new[i] / r * (double)(ell*(ell+1)) / (r*r);
        delta_omega2 += v_em[i] * delta_W * v_em[i] * dr;
        vv += v_em[i] * v_em[i] * dr;
    }
    delta_omega2 /= vv;

    printf("#\n# First-order eigenvalue shift (at kappa=1):\n");
    printf("#   delta(omega^2)/kappa = %.10e\n", delta_omega2);
    printf("#   frac shift/kappa     = %.10e\n", delta_omega2 / omega2_0);

    if (delta_omega2 < 0) {
        printf("#\n# POSITIVE FEEDBACK: Yukawa well LOWERS the eigenvalue.\n");
    } else {
        printf("#\n# NEGATIVE FEEDBACK: Yukawa well RAISES the eigenvalue.\n");
    }

    double kappa_cr = -omega2_0 / delta_omega2;
    printf("#\n# Critical coupling (perturbative): kappa_cr = %.6e\n", kappa_cr);

    /* Well depth at unit kappa */
    double well_depth = 2.0 * M_peak / rg[i_peak];
    printf("# Well depth at kappa=1: %.6e\n", well_depth);
    printf("# Well depth at kappa_cr: %.6e\n", well_depth * kappa_cr);

    /* Print profile */
    printf("#\n# Profile at kappa=1:\n");
    printf("# %10s  %14s  %14s  %14s  %14s  %14s\n",
           "r", "u(r)", "rho(r)", "M_eff(r)", "Phi(r)", "delta_W(r)");

    int step = (n > 200) ? n / 200 : 1;
    for (int i = 0; i < n; i += step) {
        double r = rg[i];
        double Phi = -M_new[i] / r;
        double delta_W = -2.0 * M_new[i] / r * (double)(ell*(ell+1))/(r*r);
        printf("%12.6f  %14.6e  %14.6e  %14.6e  %14.6e  %14.6e\n",
               rg[i], u_em[i], rho[i], M_new[i], Phi, delta_W);
    }
}

/* ==== Hartree iteration ==== */
static void run_hartree(void)
{
    int n = Ngrid;
    double rho[NMAX], M_new[NMAX];

    printf("# V9 Strong Geon — Hartree Solver\n");
    printf("# kappa=%.6e  mu=%.4f  l=%d  b_BI=%.4f  N=%d  Rmax=%.1f\n",
           kappa, mu, ell, b_BI, Ngrid, Rmax);
    printf("#\n");

    /* Initialize mass function */
    if (seed_amp > 0.0) {
        /* Seeded start: compact Yukawa well at origin */
        /* M_eff(r) = A * r * exp(-mu*r) peaks at r = 1/mu */
        for (int i = 0; i < n; i++)
            M_eff[i] = seed_amp * rg[i] * exp(-mu * rg[i]);
        printf("# SEEDED start: M_eff = %.4f * r * exp(-%.4f*r)\n",
               seed_amp, mu);
        printf("# Seed well peak at r=%.4f, depth=%.6f\n",
               1.0/mu, 2.0*seed_amp*exp(-1.0)/1.0);
    } else {
        for (int i = 0; i < n; i++) M_eff[i] = 0.0;
    }
    compute_metric();

    double relax = 0.3;
    int max_iter = 300;

    printf("# %4s  %14s  %14s  %14s  %14s  %14s  %14s\n",
           "iter", "omega2", "delta_omega2", "M_peak",
           "well_depth", "f_min", "delta_m");

    double omega2_flat = 0;

    for (int iter = 0; iter < max_iter; iter++) {
        double omega2 = solve_eigenvalue();
        if (iter == 0) omega2_flat = omega2;

        reconstruct_u();
        compute_rho(omega2, rho);

        if (mu > 0.0) {
            solve_yukawa_bvp(kappa, mu, rho, M_new);
        } else {
            integrate_einstein(kappa, rho, M_new);
        }

        /* Mass peak */
        double M_peak = 0;
        for (int i = 0; i < n; i++)
            if (M_new[i] > M_peak) M_peak = M_new[i];

        /* Minimum f(r) using selected metric */
        double f_min = 1.0;
        for (int i = 0; i < n; i++) {
            double x = 2.0 * M_new[i] / rg[i];
            double f_test;
            switch (nl_metric) {
            case 1:  f_test = exp(-x); break;
            case 2:  f_test = 1.0 / (1.0 + x); break;
            default: f_test = 1.0 - x; break;
            }
            if (f_test < f_min) f_min = f_test;
        }
        double well_depth = 1.0 - f_min;

        /* Under-relax */
        double delta = 0;
        for (int i = 0; i < n; i++) {
            double m_upd = relax * M_new[i] + (1.0 - relax) * M_eff[i];
            delta += (m_upd - M_eff[i]) * (m_upd - M_eff[i]) * dr;
            M_eff[i] = m_upd;
        }
        delta = sqrt(delta);
        compute_metric();

        printf("%4d  %14.8f  %14.6e  %14.6e  %14.6e  %14.8f  %14.6e\n",
               iter, omega2, omega2 - omega2_flat, M_peak,
               well_depth, f_min, delta);
        fflush(stdout);

        if (delta < 1e-14 && iter > 5) {
            printf("# CONVERGED at iter %d\n", iter + 1);
            break;
        }
    }

    /* Final result */
    double omega2_final = solve_eigenvalue();
    reconstruct_u();

    printf("#\n# === FINAL RESULT ===\n");
    printf("# omega^2 (flat)  = %.10f\n", omega2_flat);
    printf("# omega^2 (geon)  = %.10f\n", omega2_final);
    printf("# delta(omega^2)  = %.6e\n", omega2_final - omega2_flat);
    printf("# frac shift      = %.6e\n",
           (omega2_final - omega2_flat) / omega2_flat);

    /* Well depth and effective radius */
    double f_min_final = 1.0;
    int i_min = 0;
    double M_peak_final = 0;
    int i_mpeak = 0;
    for (int i = 0; i < n; i++) {
        if (ff[i] < f_min_final) { f_min_final = ff[i]; i_min = i; }
        if (M_eff[i] > M_peak_final) { M_peak_final = M_eff[i]; i_mpeak = i; }
    }
    printf("# M_peak          = %.6e  (at r=%.4f)\n", M_peak_final, rg[i_mpeak]);
    printf("# f_min           = %.10f  (at r=%.4f)\n", f_min_final, rg[i_min]);
    printf("# well depth      = %.6e\n", 1.0 - f_min_final);
    printf("# Phi_min         = %.6e\n", -(1.0 - f_min_final) / 2.0);

    /* Mode effective radius */
    double r2_avg = 0, norm = 0;
    for (int i = 0; i < n; i++) {
        r2_avg += rg[i] * rg[i] * u_em[i] * u_em[i] * dr;
        norm += u_em[i] * u_em[i] * dr;
    }
    double R_eff = sqrt(r2_avg / norm);
    printf("# R_eff (mode)    = %.6f\n", R_eff);
    printf("# 1/mu (Yukawa)   = %.6f\n", 1.0/mu);
    printf("# R_eff * mu      = %.6f  (mode vs range)\n", R_eff * mu);

    /* Print profile */
    printf("#\n# Profile:\n");
    printf("# %10s  %14s  %14s  %14s  %14s  %14s\n",
           "r", "u(r)", "M_eff(r)", "f(r)", "Phi(r)", "rho(r)");

    double rho_final[NMAX];
    compute_rho(omega2_final, rho_final);

    int step = (n > 400) ? n / 400 : 1;
    for (int i = 0; i < n; i += step) {
        double Phi = -M_eff[i] / rg[i];
        printf("%12.6f  %14.6e  %14.6e  %14.8f  %14.6e  %14.6e\n",
               rg[i], u_em[i], M_eff[i], ff[i], Phi, rho_final[i]);
    }
}

/* ==== Coupling scan: kappa from 0 to kappa_max ==== */
static void run_scan(double kappa_max, int nk)
{
    printf("# V9 Strong Geon Coupling Scan\n");
    printf("# l=%d  mu=%.4f  b_BI=%.4f  N=%d  Rmax=%.1f\n",
           ell, mu, b_BI, Ngrid, Rmax);
    printf("# Scanning kappa from 0 to %.4e in %d steps\n#\n", kappa_max, nk);
    printf("# %14s  %14s  %14s  %14s  %14s  %14s  %14s\n",
           "kappa", "omega2", "delta_omega2", "M_peak",
           "well_depth", "frac_shift", "R_eff");

    double omega2_flat = 0;

    for (int ik = 0; ik <= nk; ik++) {
        double kap = kappa_max * (double)ik / (double)nk;

        /* Reset to flat */
        for (int i = 0; i < Ngrid; i++) M_eff[i] = 0.0;
        compute_metric();

        /* Hartree iteration (quiet) */
        int max_iter = 300;
        double relax = 0.3;
        double rho[NMAX], M_new[NMAX];

        for (int iter = 0; iter < max_iter; iter++) {
            double omega2 = solve_eigenvalue();
            if (ik == 0 && iter == 0) omega2_flat = omega2;

            reconstruct_u();
            compute_rho(omega2, rho);

            if (mu > 0.0) {
                solve_yukawa_bvp(kap, mu, rho, M_new);
            } else {
                integrate_einstein(kap, rho, M_new);
            }

            double delta = 0;
            for (int i = 0; i < Ngrid; i++) {
                double m_upd = relax * M_new[i] + (1.0 - relax) * M_eff[i];
                delta += (m_upd - M_eff[i]) * (m_upd - M_eff[i]) * dr;
                M_eff[i] = m_upd;
            }
            compute_metric();
            delta = sqrt(delta);
            if (delta < 1e-14 && iter > 5) break;
        }

        double omega2 = solve_eigenvalue();
        reconstruct_u();

        /* M_peak */
        double M_peak = 0;
        for (int i = 0; i < Ngrid; i++)
            if (M_eff[i] > M_peak) M_peak = M_eff[i];

        /* f_min */
        double f_min = 1.0;
        for (int i = 0; i < Ngrid; i++)
            if (ff[i] < f_min) f_min = ff[i];

        /* R_eff */
        double r2_avg = 0, norm = 0;
        for (int i = 0; i < Ngrid; i++) {
            r2_avg += rg[i] * rg[i] * u_em[i] * u_em[i] * dr;
            norm += u_em[i] * u_em[i] * dr;
        }
        double R_eff = sqrt(r2_avg / norm);

        printf("%14.6e  %14.8f  %14.6e  %14.6e  %14.6e  %14.6e  %14.6e\n",
               kap, omega2, omega2 - omega2_flat, M_peak,
               1.0 - f_min,
               (omega2_flat > 1e-15) ? (omega2 - omega2_flat)/omega2_flat : 0.0,
               R_eff);
        fflush(stdout);
    }
}

/* ==== Mu scan: vary mu at fixed kappa ==== */
static void run_muscan(double mu_min, double mu_max, int nm)
{
    printf("# V9 Strong Geon — Mu Scan (perturbative)\n");
    printf("# l=%d  kappa=1  N=%d  Rmax=%.1f\n", ell, Ngrid, Rmax);
    printf("# Scanning mu from %.4f to %.4f in %d steps\n#\n",
           mu_min, mu_max, nm);
    printf("# %14s  %14s  %14s  %14s  %14s  %14s\n",
           "mu", "omega2_0", "delta_omega2", "frac_shift",
           "kappa_cr", "M_peak");

    for (int im = 0; im <= nm; im++) {
        double mu_val = mu_min + (mu_max - mu_min) * (double)im / (double)nm;

        /* Flat-space mode */
        for (int i = 0; i < Ngrid; i++) M_eff[i] = 0.0;
        compute_metric();

        double omega2_0 = solve_eigenvalue();
        reconstruct_u();

        /* Source at unit kappa */
        double rho[NMAX], M_new[NMAX];
        compute_rho(omega2_0, rho);

        if (mu_val > 1e-10) {
            solve_yukawa_bvp(1.0, mu_val, rho, M_new);
        } else {
            integrate_einstein(1.0, rho, M_new);
        }

        /* First-order shift */
        double delta_omega2 = 0.0, vv = 0.0;
        for (int i = 0; i < Ngrid; i++) {
            double r = rg[i];
            double delta_W = -2.0 * M_new[i] / r * (double)(ell*(ell+1)) / (r*r);
            delta_omega2 += v_em[i] * delta_W * v_em[i] * dr;
            vv += v_em[i] * v_em[i] * dr;
        }
        delta_omega2 /= vv;

        double kappa_cr_val = (delta_omega2 < 0) ? -omega2_0/delta_omega2 : 1e30;

        /* M_peak */
        double M_peak = 0;
        for (int i = 0; i < Ngrid; i++)
            if (M_new[i] > M_peak) M_peak = M_new[i];

        printf("%14.6e  %14.8f  %14.6e  %14.6e  %14.6e  %14.6e\n",
               mu_val, omega2_0, delta_omega2, delta_omega2/omega2_0,
               kappa_cr_val, M_peak);
        fflush(stdout);
    }
}

/* ==== Parse mode specification string ==== */
/*
 * Format: "l1:n1,l2:n2,..."  where l=angular momentum, n=occupancy
 * Examples: "1:3" = 3 modes at l=1
 *           "1:2,2:1" = 2 modes at l=1 + 1 mode at l=2
 */
static void parse_modes(const char *spec)
{
    const char *p = spec;
    n_multi_modes = 0;
    while (*p && n_multi_modes < MAX_MODES) {
        int l_val = (int)strtol(p, (char **)&p, 10);
        double n_val = 1.0;
        if (*p == ':') {
            p++;
            n_val = strtod(p, (char **)&p);
        }
        if (l_val < 1) l_val = 1;
        mode_l_arr[n_multi_modes] = l_val;
        mode_occ_arr[n_multi_modes] = n_val;
        n_multi_modes++;
        if (*p == ',') p++;
    }
}

/* ==== Multi-mode Hartree iteration ==== */
/*
 * Multiple EM modes (different l values and/or occupancies) share
 * the same gravitational potential. Total EM density = sum of individual
 * mode densities weighted by occupancy.
 *
 * For N identical modes (same l): equivalent to kappa -> N*kappa.
 * For different l: genuinely new — different radial profiles overlap.
 */
static void run_multimode(void)
{
    int n = Ngrid;
    double rho_total[NMAX], rho_mode[NMAX], M_new[NMAX];
    double omega2_m[MAX_MODES], omega2_flat_m[MAX_MODES];
    int save_ell = ell;

    printf("# V9 Strong Geon — Multi-Mode Hartree\n");
    printf("# kappa=%.6e  mu=%.4f  N=%d  Rmax=%.1f  nlmetric=%d\n",
           kappa, mu, Ngrid, Rmax, nl_metric);
    printf("# Modes (%d): ", n_multi_modes);
    double total_occ = 0;
    for (int m = 0; m < n_multi_modes; m++) {
        printf("l=%d(x%.0f) ", mode_l_arr[m], mode_occ_arr[m]);
        total_occ += mode_occ_arr[m];
    }
    printf("\n# Total occupancy: %.0f\n#\n", total_occ);

    /* Initialize */
    for (int i = 0; i < n; i++) M_eff[i] = 0.0;
    compute_metric();

    double relax = 0.3;
    int max_iter = 500;
    for (int m = 0; m < n_multi_modes; m++)
        omega2_flat_m[m] = 0;

    /* Header */
    printf("# %4s", "iter");
    for (int m = 0; m < n_multi_modes; m++)
        printf("    omega2_m%d  ", m);
    printf("  %14s  %14s  %14s\n", "M_peak", "f_min", "delta_m");

    for (int iter = 0; iter < max_iter; iter++) {
        for (int i = 0; i < n; i++) rho_total[i] = 0.0;

        /* Solve each mode in current metric */
        for (int m = 0; m < n_multi_modes; m++) {
            ell = mode_l_arr[m];
            omega2_m[m] = solve_eigenvalue();
            if (iter == 0) omega2_flat_m[m] = omega2_m[m];
            reconstruct_u();
            compute_rho(omega2_m[m], rho_mode);
            for (int i = 0; i < n; i++)
                rho_total[i] += mode_occ_arr[m] * rho_mode[i];
        }

        /* Solve Yukawa BVP with total density */
        if (mu > 0.0)
            solve_yukawa_bvp(kappa, mu, rho_total, M_new);
        else
            integrate_einstein(kappa, rho_total, M_new);

        /* Under-relax */
        double delta = 0;
        for (int i = 0; i < n; i++) {
            double m_upd = relax * M_new[i] + (1.0 - relax) * M_eff[i];
            delta += (m_upd - M_eff[i]) * (m_upd - M_eff[i]) * dr;
            M_eff[i] = m_upd;
        }
        compute_metric();
        delta = sqrt(delta);

        /* Diagnostics */
        double M_peak = 0, f_min_it = 1.0;
        for (int i = 0; i < n; i++) {
            if (M_eff[i] > M_peak) M_peak = M_eff[i];
            if (ff[i] < f_min_it) f_min_it = ff[i];
        }

        printf("%4d", iter);
        for (int m = 0; m < n_multi_modes; m++)
            printf("  %14.6f", omega2_m[m]);
        printf("  %14.6e  %14.8f  %14.6e\n", M_peak, f_min_it, delta);
        fflush(stdout);

        if (delta < 1e-14 && iter > 5) {
            printf("# CONVERGED at iter %d\n", iter + 1);
            break;
        }
    }

    /* === Final results === */
    printf("#\n# === FINAL RESULTS ===\n");

    double f_min_final = 1.0;
    int i_fmin = 0;
    double M_pk = 0;
    for (int i = 0; i < n; i++) {
        if (ff[i] < f_min_final) { f_min_final = ff[i]; i_fmin = i; }
        if (M_eff[i] > M_pk) M_pk = M_eff[i];
    }
    printf("# f_min           = %.10f  (at r=%.4f)\n", f_min_final, rg[i_fmin]);
    printf("# well_depth      = %.6e\n", 1.0 - f_min_final);
    printf("# M_peak          = %.6e\n", M_pk);
    printf("#\n");

    double E_total = 0;
    for (int m = 0; m < n_multi_modes; m++) {
        ell = mode_l_arr[m];
        omega2_m[m] = solve_eigenvalue();
        reconstruct_u();

        double omega_m = (omega2_m[m] > 0) ? sqrt(omega2_m[m]) : 0;
        double E_MeV = omega_m * 197.3;

        double r2a = 0, nrm = 0;
        for (int i = 0; i < n; i++) {
            r2a += rg[i] * rg[i] * u_em[i] * u_em[i] * dr;
            nrm += u_em[i] * u_em[i] * dr;
        }
        double R_eff = sqrt(r2a / nrm);

        printf("# Mode %d: l=%d  occ=%.0f  omega^2=%.10f  omega=%.6f  "
               "E=%.1f MeV  R_eff=%.4f  shift=%.4e\n",
               m, mode_l_arr[m], mode_occ_arr[m],
               omega2_m[m], omega_m, E_MeV, R_eff,
               (omega2_flat_m[m] > 0) ?
               (omega2_m[m] - omega2_flat_m[m]) / omega2_flat_m[m] : 0.0);

        E_total += mode_occ_arr[m] * E_MeV;
    }

    printf("#\n# Total geon energy: %.1f MeV\n", E_total);
    printf("# Proton mass:      938.3 MeV\n");
    printf("# Ratio E/M_p:      %.3f\n", E_total / 938.3);

    /* Equivalent single-mode comparison */
    printf("#\n# === COMPARISON ===\n");
    printf("# For identical modes (all same l), equivalent to kappa_eff = %.1f * kappa = %.4f\n",
           total_occ, total_occ * kappa);

    ell = save_ell;
}

/* ==== WKB Lifetime computation ==== */
/*
 * After Hartree converges in a compact box, compute the quasi-bound
 * state lifetime by WKB tunneling through the centrifugal barrier.
 *
 * The Regge-Wheeler equation in tortoise coordinates:
 *   d^2 Psi/dr*^2 + [omega^2 - V_RW] Psi = 0
 * where V_RW = f(r) * l(l+1) / r^2 and dr* = dr / f(r).
 *
 * WKB action through the barrier:
 *   S = int_{r1}^{r2} sqrt(V_RW - omega^2) / f(r) dr
 *
 * Width: Gamma ~ omega * exp(-2S)
 * Lifetime: tau = 1/Gamma
 */
#define NEXT 20001

static void run_lifetime(void)
{
    int n = Ngrid;
    double rho[NMAX], M_new[NMAX];
    double ll1 = (double)(ell * (ell + 1));

    printf("# V9 Strong Geon — WKB Lifetime Analysis\n");
    printf("# kappa=%.6e  mu=%.4f  l=%d  N=%d  Rmax=%.4f\n",
           kappa, mu, ell, Ngrid, Rmax);
    printf("#\n");

    /* === Step 1: Hartree iteration === */
    if (seed_amp > 0.0) {
        for (int i = 0; i < n; i++)
            M_eff[i] = seed_amp * rg[i] * exp(-mu * rg[i]);
    } else {
        for (int i = 0; i < n; i++) M_eff[i] = 0.0;
    }
    compute_metric();

    double relax = 0.3;
    int max_iter = 300;
    double omega2_flat = 0;

    for (int iter = 0; iter < max_iter; iter++) {
        double omega2_i = solve_eigenvalue();
        if (iter == 0) omega2_flat = omega2_i;
        reconstruct_u();
        compute_rho(omega2_i, rho);

        if (mu > 0.0)
            solve_yukawa_bvp(kappa, mu, rho, M_new);
        else
            integrate_einstein(kappa, rho, M_new);

        double delta = 0;
        for (int i = 0; i < n; i++) {
            double m_upd = relax * M_new[i] + (1.0 - relax) * M_eff[i];
            delta += (m_upd - M_eff[i]) * (m_upd - M_eff[i]) * dr;
            M_eff[i] = m_upd;
        }
        compute_metric();
        delta = sqrt(delta);
        if (delta < 1e-14 && iter > 5) break;
    }

    double omega2 = solve_eigenvalue();
    reconstruct_u();

    if (omega2 <= 0) {
        printf("# ERROR: omega^2 = %.6e <= 0\n", omega2);
        return;
    }
    double omega = sqrt(omega2);

    /* Report Hartree result */
    double f_min = 1.0;
    double M_pk = 0;
    for (int i = 0; i < n; i++) {
        if (ff[i] < f_min) f_min = ff[i];
        if (M_eff[i] > M_pk) M_pk = M_eff[i];
    }
    double r2a = 0, nrm = 0;
    for (int i = 0; i < n; i++) {
        r2a += rg[i] * rg[i] * u_em[i] * u_em[i] * dr;
        nrm += u_em[i] * u_em[i] * dr;
    }
    double R_eff = sqrt(r2a / nrm);

    printf("# === Hartree Result ===\n");
    printf("# omega^2_flat = %.10f\n", omega2_flat);
    printf("# omega^2      = %.10f\n", omega2);
    printf("# omega        = %.10f\n", omega);
    printf("# f_min        = %.10f\n", f_min);
    printf("# well_depth   = %.6f\n", 1.0 - f_min);
    printf("# M_peak       = %.6e\n", M_pk);
    printf("# R_eff        = %.6f\n", R_eff);
    printf("# frac_shift   = %.6e\n", (omega2 - omega2_flat) / omega2_flat);

    /* === Step 2: Build extended grid beyond box === */
    double r_out_flat = sqrt(ll1 / omega2);
    double R_ext = r_out_flat + 2.0;
    if (R_ext < 3.0 * Rmax) R_ext = 3.0 * Rmax;

    int N_ext = NEXT;
    double dr_e = R_ext / (N_ext - 1);
    double M_at_Rmax = M_eff[n - 1];

    static double r_e[NEXT], M_e[NEXT], f_e[NEXT], V_e[NEXT];

    for (int i = 0; i < N_ext; i++) {
        r_e[i] = (i + 0.5) * dr_e;

        if (r_e[i] <= rg[n - 1]) {
            /* Interpolate from Hartree grid */
            double j_f = r_e[i] / dr - 0.5;
            int j = (int)j_f;
            if (j < 0) j = 0;
            if (j >= n - 1) j = n - 2;
            double t = j_f - j;
            if (t < 0.0) t = 0.0;
            if (t > 1.0) t = 1.0;
            M_e[i] = (1.0 - t) * M_eff[j] + t * M_eff[j + 1];
        } else {
            /* Yukawa extrapolation beyond box */
            M_e[i] = M_at_Rmax * exp(-mu * (r_e[i] - Rmax));
        }

        {
            double x = 2.0 * M_e[i] / r_e[i];
            switch (nl_metric) {
            case 1:  f_e[i] = exp(-x); break;
            case 2:  f_e[i] = 1.0 / (1.0 + x); break;
            default: f_e[i] = 1.0 - x;
                if (f_e[i] < 0.01) f_e[i] = 0.01;
                break;
            }
        }

        /* Regge-Wheeler potential */
        V_e[i] = f_e[i] * ll1 / (r_e[i] * r_e[i]);
    }

    printf("#\n# === Extended Domain ===\n");
    printf("# R_ext        = %.4f\n", R_ext);
    printf("# N_ext        = %d\n", N_ext);
    printf("# dr_ext       = %.6e\n", dr_e);
    printf("# r_out (flat) = %.6f\n", r_out_flat);
    printf("# M_eff(Rmax)  = %.6e\n", M_at_Rmax);

    /* === Step 3: Find barrier structure === */
    /*
     * Phases of V_RW(r) vs omega^2:
     *   1: V > omega^2 (centrifugal, small r)
     *   2: V < omega^2 (well — mode lives here)
     *   3: V > omega^2 (barrier — tunneling)
     *   4: V < omega^2 (free — mode escapes)
     */
    int phase = 1;
    int i_ws = -1, i_we = -1;    /* well start/end */
    int i_bs = -1, i_be = -1;    /* barrier start/end */
    double V_bmax = 0;

    for (int i = 0; i < N_ext; i++) {
        double V = V_e[i];
        if (phase == 1 && V < omega2) {
            phase = 2;
            i_ws = i;
        } else if (phase == 2 && V >= omega2) {
            phase = 3;
            i_we = i - 1;
            i_bs = i;
            V_bmax = V;
        } else if (phase == 3) {
            if (V > V_bmax) V_bmax = V;
            if (V < omega2) {
                phase = 4;
                i_be = i - 1;
            }
        }
    }
    if (phase == 3) i_be = N_ext - 1;  /* barrier extends to grid edge */
    if (phase == 2) i_we = N_ext - 1;  /* well extends to grid edge */

    int has_barrier = (i_bs >= 0 && i_be > i_bs);

    printf("#\n# === Barrier Structure ===\n");
    if (i_ws >= 0)
        printf("# Well start:    r = %.6f\n", r_e[i_ws]);
    else
        printf("# Well start:    NOT FOUND\n");
    if (has_barrier) {
        printf("# Well end:      r = %.6f\n", r_e[i_we]);
        printf("# Barrier start: r = %.6f\n", r_e[i_bs]);
        printf("# Barrier end:   r = %.6f\n", r_e[i_be]);
        printf("# Barrier width: %.6f\n", r_e[i_be] - r_e[i_bs]);
        printf("# V_barr_max:    %.6e\n", V_bmax);
        printf("# V_max/omega^2: %.4f\n", V_bmax / omega2);
    } else {
        printf("# NO BARRIER found after well\n");
    }

    /* === Step 4: WKB tunneling integral === */
    /*
     * S = int sqrt(V_RW - omega^2) / f dr   (tortoise coordinate)
     */
    double S_wkb = 0.0;
    if (has_barrier) {
        for (int i = i_bs; i <= i_be; i++) {
            if (V_e[i] > omega2 && f_e[i] > 1e-10) {
                S_wkb += sqrt(V_e[i] - omega2) / f_e[i] * dr_e;
            }
        }
    }

    double P_tunnel = exp(-2.0 * S_wkb);
    double Gamma_c = omega * P_tunnel;
    double tau_c = (Gamma_c > 1e-300) ? 1.0 / Gamma_c : 1e300;

    /* Physical units (L_0 = 1 fm) */
    double hbar_c = 197.3;          /* MeV fm */
    double t_code = 3.336e-24;      /* s (= L_0/c = 1fm / 3e23 fm/s) */
    double tau_phys = tau_c * t_code;
    double Gamma_MeV = Gamma_c * hbar_c;
    double omega_MeV = omega * hbar_c;
    double tau_hadron = 1e-23;      /* hadronic timescale (s) */

    printf("#\n# === WKB LIFETIME RESULTS ===\n");
    printf("# WKB action S:    %.6f\n", S_wkb);
    printf("# exp(-2S):        %.6e\n", P_tunnel);
    printf("#\n# Code units:\n");
    printf("#   omega:         %.10f\n", omega);
    printf("#   Gamma:         %.6e\n", Gamma_c);
    printf("#   tau:           %.6e\n", tau_c);
    printf("#\n# Physical (L_0 = 1 fm):\n");
    printf("#   omega:         %.4f MeV\n", omega_MeV);
    printf("#   Gamma:         %.4e MeV\n", Gamma_MeV);
    printf("#   tau:           %.4e s\n", tau_phys);
    printf("#   tau/tau_had:   %.4e\n", tau_phys / tau_hadron);
    printf("#   Q = omega/Gamma: %.4e\n",
           Gamma_c > 1e-300 ? omega / Gamma_c : 1e300);

    if (!has_barrier) {
        double T_osc = 2.0 * M_PI / omega;
        printf("#\n# *** NO BARRIER — mode leaks freely without box ***\n");
        printf("# Free-leak time ~ 1 period = %.4e code = %.4e s\n",
               T_osc, T_osc * t_code);
    } else {
        if (tau_phys > 100 * tau_hadron)
            printf("#\n# VERY LONG-LIVED: tau = %.0f x tau_hadron\n",
                   tau_phys / tau_hadron);
        else if (tau_phys > tau_hadron)
            printf("#\n# LONG-LIVED: tau = %.1f x tau_hadron\n",
                   tau_phys / tau_hadron);
        else
            printf("#\n# SHORT-LIVED: tau = %.3f x tau_hadron\n",
                   tau_phys / tau_hadron);
    }

    /* === Print extended potential profile === */
    printf("#\n# Extended potential profile:\n");
    printf("# %12s  %14s  %14s  %14s  %14s\n",
           "r", "V_RW", "f(r)", "V-omega2", "WKB_integrand");

    int step = (N_ext > 500) ? N_ext / 500 : 1;
    for (int i = 0; i < N_ext; i += step) {
        double diff = V_e[i] - omega2;
        double integ = (diff > 0 && f_e[i] > 1e-10) ?
                        sqrt(diff) / f_e[i] : 0.0;
        printf("%14.6f  %14.6e  %14.8f  %14.6e  %14.6e\n",
               r_e[i], V_e[i], f_e[i], diff, integ);
    }
}

/* ==== Main ==== */
int main(int argc, char *argv[])
{
    int mode_hartree = 0;
    int mode_scan = 0;
    int mode_perturb = 0;
    int mode_muscan = 0;
    int mode_lifetime = 0;
    int mode_multimode = 0;
    double kappa_max = 100.0;
    int nk = 20;
    double mu_min = 0.1, mu_max = 10.0;
    int nm = 20;

    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "-hartree"))   { mode_hartree = 1; continue; }
        if (!strcmp(argv[i], "-scan"))      { mode_scan = 1; continue; }
        if (!strcmp(argv[i], "-perturb"))   { mode_perturb = 1; continue; }
        if (!strcmp(argv[i], "-muscan"))    { mode_muscan = 1; continue; }
        if (!strcmp(argv[i], "-lifetime"))  { mode_lifetime = 1; continue; }
        if (!strcmp(argv[i], "-multimode")) { mode_multimode = 1; continue; }
        if (i + 1 < argc) {
            if (!strcmp(argv[i], "-kappa"))    { kappa = atof(argv[++i]); continue; }
            if (!strcmp(argv[i], "-kmax"))     { kappa_max = atof(argv[++i]); continue; }
            if (!strcmp(argv[i], "-nk"))       { nk = atoi(argv[++i]); continue; }
            if (!strcmp(argv[i], "-mu"))       { mu = atof(argv[++i]); continue; }
            if (!strcmp(argv[i], "-mumin"))    { mu_min = atof(argv[++i]); continue; }
            if (!strcmp(argv[i], "-mumax"))    { mu_max = atof(argv[++i]); continue; }
            if (!strcmp(argv[i], "-nm"))       { nm = atoi(argv[++i]); continue; }
            if (!strcmp(argv[i], "-seed"))     { seed_amp = atof(argv[++i]); continue; }
            if (!strcmp(argv[i], "-bBI"))      { b_BI = atof(argv[++i]); continue; }
            if (!strcmp(argv[i], "-nlmetric")) { nl_metric = atoi(argv[++i]); continue; }
            if (!strcmp(argv[i], "-modes"))    { parse_modes(argv[++i]); continue; }
            if (!strcmp(argv[i], "-l"))        { ell = atoi(argv[++i]); continue; }
            if (!strcmp(argv[i], "-N"))        { Ngrid = atoi(argv[++i]); continue; }
            if (!strcmp(argv[i], "-Rmax"))     { Rmax = atof(argv[++i]); continue; }
        }
    }

    if (ell < 1) {
        fprintf(stderr, "ERROR: EM modes require l >= 1 (got l=%d)\n", ell);
        return 1;
    }

    setup_grid();

    if (mode_perturb)
        run_perturbative();
    else if (mode_muscan)
        run_muscan(mu_min, mu_max, nm);
    else if (mode_scan)
        run_scan(kappa_max, nk);
    else if (mode_lifetime)
        run_lifetime();
    else if (mode_multimode) {
        if (n_multi_modes == 0) {
            /* Default: 3 x l=1 (proton-like) */
            n_multi_modes = 1;
            mode_l_arr[0] = ell;
            mode_occ_arr[0] = 3;
        }
        run_multimode();
    }
    else if (mode_hartree)
        run_hartree();
    else {
        printf("Usage: strong_geon -perturb|-hartree|-scan|-muscan|-lifetime|-multimode [options]\n");
        printf("  -perturb         First-order perturbation theory\n");
        printf("  -hartree         Full Hartree self-consistent solve\n");
        printf("  -scan            Scan kappa from 0 to kmax\n");
        printf("  -muscan          Scan mu at fixed kappa=1 (perturbative)\n");
        printf("  -lifetime        WKB lifetime of quasi-bound geon\n");
        printf("  -multimode       Multi-mode Hartree (use with -modes)\n");
        printf("\nParameters:\n");
        printf("  -kappa <K>       Coupling constant [%.2e]\n", kappa);
        printf("  -mu <M>          Yukawa mass [%.4f]\n", mu);
        printf("  -kmax <K>        Max kappa for scan [%.2e]\n", kappa_max);
        printf("  -nk <N>          Number of kappa scan steps [%d]\n", nk);
        printf("  -mumin <M>       Min mu for scan [%.4f]\n", mu_min);
        printf("  -mumax <M>       Max mu for scan [%.4f]\n", mu_max);
        printf("  -nm <N>          Number of mu scan steps [%d]\n", nm);
        printf("  -seed <A>        Seed well amplitude [%.2f]\n", seed_amp);
        printf("  -bBI <b>         Born-Infeld field strength [%.2f]\n", b_BI);
        printf("  -nlmetric <N>    Metric: 0=linear, 1=exp, 2=Padé [%d]\n", nl_metric);
        printf("  -modes <spec>    Mode specification: l1:n1,l2:n2,... (e.g. 1:3 or 1:2,2:1)\n");
        printf("  -l <L>           Angular momentum [%d]\n", ell);
        printf("  -N <N>           Grid points [%d]\n", Ngrid);
        printf("  -Rmax <R>        Domain size [%.1f]\n", Rmax);
        printf("\nPhysical reference (f_2 meson):\n");
        printf("  mu_f2 = 6.47 fm^{-1} (m_f2 = 1.275 GeV)\n");
        printf("  kappa_strong ~ 35 code units (L_0 = 1 fm)\n");
    }

    return 0;
}
