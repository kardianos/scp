/*
 * tov_geon.c — TOV Strong Geon: Two-Potential Yukawa Hartree Solver
 *
 * V10: Resolves the V9 "Redshift Paradox" (Q9) by decoupling
 * the temporal and spatial metric functions.
 *
 * V9:  ds² = -f(r) dt² + dr²/f(r) + r² dΩ²    (single f)
 * V10: ds² = -A(r) dt² + dr²/C(r) + r² dΩ²    (independent A, C)
 *
 * Two Yukawa BVPs determine the metric:
 *   M_t'' - μ² M_t = -κ r ρ           (temporal, from energy density)
 *   M_s'' - μ² M_s = -κ r ρ_s         (spatial, from mixed source)
 *
 * where ρ_s = α_s·ρ + (1-α_s)·p_r:
 *   α_s = 1: V9 recovery (A = C = f, single potential)
 *   α_s = 0: Fierz-Pauli massive spin-2 (spatial sourced by p_r)
 *
 * Physics: For EM standing waves, p_r = ρ_E - ρ_B changes sign.
 * Near the origin (magnetic dominated), p_r < 0, making the spatial
 * metric C > 1 (space compresses). This prevents the catastrophic
 * spatial stretching that caused the V9 redshift paradox.
 *
 * EM wave equation: -P v'' + W v = ω² v
 * where P = A·C, W = A·l(l+1)/r² + P(β'/4 + β²/16), β = (ln P)'.
 * Identical form to V9 but with P = A·C instead of P = f².
 *
 * Build: cd src && make
 * Usage: tov_geon -hartree -kappa <K> -mu <M> -alpha <A> [options]
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
#define NEXT 20001

/* ---- Parameters ---- */
static double kappa   = 1.0;
static double mu      = 1.0;
static double alpha_s = 1.0;   /* 0=Fierz-Pauli, 1=V9 recovery */
static int    ell     = 1;
static int    nl_metric = 0;   /* 0=linear, 1=exp, 2=Padé */
static double seed_amp = 0.0;

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

/* ---- Two-potential metric ---- */
static double M_t[NMAX];    /* temporal mass function */
static double M_s[NMAX];    /* spatial mass function */
static double AA[NMAX];     /* A(r) = -g_tt */
static double CC[NMAX];     /* C(r) = 1/g_rr */
static double PP[NMAX];     /* P = A·C */
static double bb[NMAX];     /* β = (ln P)' = A'/A + C'/C */
static double bbp[NMAX];    /* β' */

/* ---- EM mode ---- */
static double u_em[NMAX];
static double v_em[NMAX];

/* ==== Grid ==== */
static void setup_grid(void)
{
    if (Ngrid > NMAX) Ngrid = NMAX;
    dr = Rmax / (Ngrid - 1);
    for (int i = 0; i < Ngrid; i++)
        rg[i] = (i + 0.5) * dr;
}

/* ==== Nonlinear metric function ==== */
static double nl_func(double x)
{
    switch (nl_metric) {
    case 1:  return exp(-x);
    case 2:  return (1.0 + x > 0.01) ? 1.0 / (1.0 + x) : 100.0;
    default: { double f = 1.0 - x; return (f > 0.01) ? f : 0.01; }
    }
}

/* ==== Compute two-potential metric from M_t, M_s ==== */
static void compute_metric(void)
{
    int n = Ngrid;
    double Ap[NMAX], Cp[NMAX];

    for (int i = 0; i < n; i++) {
        double x_t = 2.0 * M_t[i] / rg[i];
        double x_s = 2.0 * M_s[i] / rg[i];
        AA[i] = nl_func(x_t);
        CC[i] = nl_func(x_s);
        PP[i] = AA[i] * CC[i];
    }

    /* A' and C' by centered differences */
    Ap[0] = (AA[1] - AA[0]) / dr;
    Cp[0] = (CC[1] - CC[0]) / dr;
    for (int i = 1; i < n - 1; i++) {
        Ap[i] = (AA[i+1] - AA[i-1]) / (2.0 * dr);
        Cp[i] = (CC[i+1] - CC[i-1]) / (2.0 * dr);
    }
    Ap[n-1] = (AA[n-1] - AA[n-2]) / dr;
    Cp[n-1] = (CC[n-1] - CC[n-2]) / dr;

    /* β = A'/A + C'/C */
    for (int i = 0; i < n; i++)
        bb[i] = Ap[i] / AA[i] + Cp[i] / CC[i];

    /* β' by centered differences */
    bbp[0] = (bb[1] - bb[0]) / dr;
    for (int i = 1; i < n - 1; i++)
        bbp[i] = (bb[i+1] - bb[i-1]) / (2.0 * dr);
    bbp[n-1] = (bb[n-1] - bb[n-2]) / dr;
}

/* ==== Effective potential W(r) ==== */
static double W_eff_pot(int i)
{
    double r = rg[i];
    double ll1 = (double)(ell * (ell + 1));
    return AA[i] * ll1 / (r * r)
         + PP[i] * (bbp[i] / 4.0 + bb[i] * bb[i] / 16.0);
}

/* ==== Eigenvalue solver: implicit imaginary time ==== */
static double solve_eigenvalue(void)
{
    int n = Ngrid;
    double At[NMAX], Bt[NMAX], Ct[NMAX], Dt[NMAX];
    double v_new[NMAX];
    double dr2 = dr * dr;

    /* Initial guess: r^l * gaussian */
    double norm2 = 0;
    for (int i = 0; i < n; i++) {
        v_em[i] = pow(rg[i], ell) * exp(-rg[i] * rg[i] / (Rmax * Rmax / 4.0));
        norm2 += v_em[i] * v_em[i] * dr;
    }
    double inv = 1.0 / sqrt(norm2);
    for (int i = 0; i < n; i++) v_em[i] *= inv;

    double dtau = 0.3;
    int niter = 8000;

    for (int iter = 0; iter < niter; iter++) {
        /* i=0: ghost v[-1] = -v[0] (odd parity) */
        double P0 = PP[0];
        At[0] = 0.0;
        Bt[0] = 1.0 + dtau * (3.0 * P0 / dr2 + W_eff_pot(0));
        Ct[0] = -dtau * P0 / dr2;
        Dt[0] = v_em[0];

        for (int i = 1; i < n - 1; i++) {
            double Pi = PP[i];
            At[i] = -dtau * Pi / dr2;
            Bt[i] = 1.0 + dtau * (2.0 * Pi / dr2 + W_eff_pot(i));
            Ct[i] = -dtau * Pi / dr2;
            Dt[i] = v_em[i];
        }

        /* Dirichlet at Rmax */
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

    /* Compute eigenvalue: ω² = <v|H|v>/<v|v> */
    double hv = 0.0, vv = 0.0;
    {
        double vdd = (v_em[1] - 3.0 * v_em[0]) / dr2;  /* ghost */
        double Hv = -PP[0] * vdd + W_eff_pot(0) * v_em[0];
        hv += v_em[0] * Hv * dr;
        vv += v_em[0] * v_em[0] * dr;
    }
    for (int i = 1; i < n - 1; i++) {
        double vdd = (v_em[i+1] - 2.0 * v_em[i] + v_em[i-1]) / dr2;
        double Hv = -PP[i] * vdd + W_eff_pot(i) * v_em[i];
        hv += v_em[i] * Hv * dr;
        vv += v_em[i] * v_em[i] * dr;
    }

    return hv / vv;
}

/* ==== Reconstruct u from v: u = v / P^{1/4} ==== */
static void reconstruct_u(void)
{
    for (int i = 0; i < Ngrid; i++) {
        double p14 = pow(PP[i], 0.25);
        u_em[i] = (p14 > 1e-10) ? v_em[i] / p14 : 0.0;
    }
}

/* ==== Compute EM energy density ρ and radial pressure p_r ==== */
/*
 * For ds² = -A dt² + dr²/C + r² dΩ²:
 *   ρ = N_l/(8πr²) [ω²u²/A + C(u')²]
 *   p_r = N_l/(8πr²) [ω²u²/A - C(u')²]
 *
 * ρ > 0 always. p_r changes sign (negative near origin where B > E).
 */
static void compute_rho_pr(double omega2, double *rho, double *pr)
{
    int n = Ngrid;
    double up[NMAX];
    double N_l = (double)(ell * (ell + 1)) / (double)(2 * ell + 1);

    /* u' */
    up[0] = (u_em[1] + u_em[0]) / (2.0 * dr);
    for (int i = 1; i < n - 1; i++)
        up[i] = (u_em[i+1] - u_em[i-1]) / (2.0 * dr);
    up[n-1] = (u_em[n-1] - u_em[n-2]) / dr;

    for (int i = 0; i < n; i++) {
        double r = rg[i];
        double rho_E = N_l / (8.0 * M_PI * r * r) * omega2 * u_em[i] * u_em[i] / AA[i];
        double rho_B = N_l / (8.0 * M_PI * r * r) * CC[i] * up[i] * up[i];
        rho[i] = rho_E + rho_B;
        pr[i]  = rho_E - rho_B;
    }
}

/* ==== Solve Yukawa BVP: M'' - μ²M = -κ r source ==== */
static void solve_yukawa_bvp(double kap, double mu_val, double *source, double *M_out)
{
    int n = Ngrid;
    double h = dr;
    double h2 = h * h;
    double mu2h2 = mu_val * mu_val * h2;
    double emh = exp(-mu_val * h);

    double At[NMAX], Bt[NMAX], Ct[NMAX], Dt[NMAX];

    /* i=0: ghost M_{-1} = -M_0 */
    At[0] = 0.0;
    Bt[0] = -(3.0 + mu2h2);
    Ct[0] = 1.0;
    Dt[0] = -h2 * kap * rg[0] * source[0];

    for (int i = 1; i < n - 1; i++) {
        At[i] = 1.0;
        Bt[i] = -(2.0 + mu2h2);
        Ct[i] = 1.0;
        Dt[i] = -h2 * kap * rg[i] * source[i];
    }

    /* Robin BC at Rmax */
    At[n-1] = 1.0;
    Bt[n-1] = emh - 2.0 - mu2h2;
    Ct[n-1] = 0.0;
    Dt[n-1] = -h2 * kap * rg[n-1] * source[n-1];

    /* Thomas forward sweep */
    double cp[NMAX], dp[NMAX];
    cp[0] = Ct[0] / Bt[0];
    dp[0] = Dt[0] / Bt[0];
    for (int i = 1; i < n; i++) {
        double den = Bt[i] - At[i] * cp[i-1];
        cp[i] = Ct[i] / den;
        dp[i] = (Dt[i] - At[i] * dp[i-1]) / den;
    }
    M_out[n-1] = dp[n-1];
    for (int i = n - 2; i >= 0; i--)
        M_out[i] = dp[i] - cp[i] * M_out[i+1];
}

/* ==== Hartree iteration ==== */
static void run_hartree(void)
{
    int n = Ngrid;
    double rho[NMAX], pr[NMAX], rho_s[NMAX];
    double Mt_new[NMAX], Ms_new[NMAX];

    printf("# V10 TOV Geon — Two-Potential Hartree Solver\n");
    printf("# kappa=%.6e  mu=%.4f  alpha_s=%.4f  l=%d  nlmetric=%d  N=%d  Rmax=%.4f\n",
           kappa, mu, alpha_s, ell, nl_metric, Ngrid, Rmax);
    printf("#\n");

    /* Initialize */
    if (seed_amp > 0.0) {
        for (int i = 0; i < n; i++) {
            M_t[i] = seed_amp * rg[i] * exp(-mu * rg[i]);
            M_s[i] = alpha_s * M_t[i];
        }
        printf("# SEEDED start: amp=%.4f\n", seed_amp);
    } else {
        for (int i = 0; i < n; i++) { M_t[i] = 0.0; M_s[i] = 0.0; }
    }
    compute_metric();

    double relax = 0.3;
    int max_iter = 300;
    double omega2_flat = 0;

    printf("# %4s  %14s  %14s  %14s  %14s  %14s  %14s  %14s\n",
           "iter", "omega2", "A_min", "C_min", "C_max", "P_min", "well_A", "delta_m");

    for (int iter = 0; iter < max_iter; iter++) {
        double omega2 = solve_eigenvalue();
        if (iter == 0) omega2_flat = omega2;

        reconstruct_u();
        compute_rho_pr(omega2, rho, pr);

        /* Spatial source: ρ_s = α_s·ρ + (1-α_s)·p_r */
        for (int i = 0; i < n; i++)
            rho_s[i] = alpha_s * rho[i] + (1.0 - alpha_s) * pr[i];

        /* Two Yukawa BVPs */
        solve_yukawa_bvp(kappa, mu, rho, Mt_new);
        solve_yukawa_bvp(kappa, mu, rho_s, Ms_new);

        /* Under-relax both */
        double delta = 0;
        for (int i = 0; i < n; i++) {
            double mt = relax * Mt_new[i] + (1.0 - relax) * M_t[i];
            double ms = relax * Ms_new[i] + (1.0 - relax) * M_s[i];
            delta += (mt - M_t[i]) * (mt - M_t[i]) * dr;
            delta += (ms - M_s[i]) * (ms - M_s[i]) * dr;
            M_t[i] = mt;
            M_s[i] = ms;
        }
        delta = sqrt(delta);
        compute_metric();

        /* Diagnostics */
        double A_min = 1, C_min = 1, C_max = 1, P_min = 1;
        for (int i = 0; i < n; i++) {
            if (AA[i] < A_min) A_min = AA[i];
            if (CC[i] < C_min) C_min = CC[i];
            if (CC[i] > C_max) C_max = CC[i];
            if (PP[i] < P_min) P_min = PP[i];
        }

        printf("%4d  %14.8f  %14.8f  %14.8f  %14.8f  %14.8f  %14.6e  %14.6e\n",
               iter, omega2, A_min, C_min, C_max, P_min, 1.0 - A_min, delta);
        fflush(stdout);

        if (delta < 1e-14 && iter > 5) {
            printf("# CONVERGED at iter %d\n", iter + 1);
            break;
        }
    }

    /* === Final results === */
    double omega2_final = solve_eigenvalue();
    reconstruct_u();

    printf("#\n# === FINAL RESULT ===\n");
    printf("# omega^2 (flat)  = %.10f\n", omega2_flat);
    printf("# omega^2 (geon)  = %.10f\n", omega2_final);
    printf("# delta(omega^2)  = %.6e\n", omega2_final - omega2_flat);
    printf("# frac shift      = %.6e\n",
           (omega2_flat > 1e-15) ? (omega2_final - omega2_flat) / omega2_flat : 0.0);

    double A_min = 1, C_min = 1, C_max = 1, P_min = 1;
    int i_Amin = 0, i_Cmax = 0;
    double Mt_peak = 0, Ms_peak = 0, Ms_min = 0;
    for (int i = 0; i < n; i++) {
        if (AA[i] < A_min) { A_min = AA[i]; i_Amin = i; }
        if (CC[i] < C_min) C_min = CC[i];
        if (CC[i] > C_max) { C_max = CC[i]; i_Cmax = i; }
        if (PP[i] < P_min) P_min = PP[i];
        if (M_t[i] > Mt_peak) Mt_peak = M_t[i];
        if (M_s[i] > Ms_peak) Ms_peak = M_s[i];
        if (M_s[i] < Ms_min) Ms_min = M_s[i];
    }

    printf("# A_min           = %.10f  (at r=%.4f)\n", A_min, rg[i_Amin]);
    printf("# C_min           = %.10f\n", C_min);
    printf("# C_max           = %.10f  (at r=%.4f)\n", C_max, rg[i_Cmax]);
    printf("# P_min           = %.10f\n", P_min);
    printf("# well_depth_A    = %.6e  (temporal)\n", 1.0 - A_min);
    printf("# M_t_peak        = %.6e\n", Mt_peak);
    printf("# M_s_peak        = %.6e\n", Ms_peak);
    printf("# M_s_min         = %.6e  (%s)\n", Ms_min,
           Ms_min < -1e-15 ? "REPULSIVE spatial" : "attractive");
    printf("# C at A_min      = %.6f  (spatial metric at deepest temporal well)\n",
           CC[i_Amin]);

    /* Mode radius */
    double r2a = 0, nrm = 0;
    for (int i = 0; i < n; i++) {
        r2a += rg[i] * rg[i] * u_em[i] * u_em[i] * dr;
        nrm += u_em[i] * u_em[i] * dr;
    }
    double R_eff = sqrt(r2a / nrm);

    printf("# R_eff           = %.6f\n", R_eff);
    printf("# 1/mu            = %.6f\n", 1.0 / mu);

    /* Energy */
    double omega = (omega2_final > 0) ? sqrt(omega2_final) : 0;
    double E_MeV = omega * 197.3;
    printf("# omega           = %.10f\n", omega);
    printf("# E_mode          = %.1f MeV\n", E_MeV);
    printf("# M_proton        = 938.3 MeV\n");
    printf("# E/M_p           = %.4f\n", E_MeV / 938.3);

    /* Print profile */
    printf("#\n# Profile:\n");
    printf("# %10s  %14s  %14s  %14s  %14s  %14s  %14s  %14s\n",
           "r", "u(r)", "A(r)", "C(r)", "M_t(r)", "M_s(r)", "rho(r)", "p_r(r)");

    double rho_f[NMAX], pr_f[NMAX];
    compute_rho_pr(omega2_final, rho_f, pr_f);

    int step = (n > 400) ? n / 400 : 1;
    for (int i = 0; i < n; i += step) {
        printf("%12.6f  %14.6e  %14.8f  %14.8f  %14.6e  %14.6e  %14.6e  %14.6e\n",
               rg[i], u_em[i], AA[i], CC[i], M_t[i], M_s[i], rho_f[i], pr_f[i]);
    }
}

/* ==== Coupling scan ==== */
static void run_scan(double kappa_max, int nk)
{
    printf("# V10 TOV Geon — Coupling Scan\n");
    printf("# alpha_s=%.4f  mu=%.4f  l=%d  N=%d  Rmax=%.4f  nlmetric=%d\n",
           alpha_s, mu, ell, Ngrid, Rmax, nl_metric);
    printf("# Scanning kappa from 0 to %.4e in %d steps\n#\n", kappa_max, nk);
    printf("# %14s  %14s  %14s  %14s  %14s  %14s  %14s  %14s\n",
           "kappa", "omega2", "A_min", "C_min", "C_max", "well_A", "E_MeV", "R_eff");

    double omega2_flat = 0;

    for (int ik = 0; ik <= nk; ik++) {
        double kap = kappa_max * (double)ik / (double)nk;

        for (int i = 0; i < Ngrid; i++) { M_t[i] = 0; M_s[i] = 0; }
        compute_metric();

        int max_iter = 300;
        double relax = 0.3;
        double rho[NMAX], pr[NMAX], rho_s[NMAX];
        double Mt_new[NMAX], Ms_new[NMAX];

        for (int iter = 0; iter < max_iter; iter++) {
            double omega2 = solve_eigenvalue();
            if (ik == 0 && iter == 0) omega2_flat = omega2;

            reconstruct_u();
            compute_rho_pr(omega2, rho, pr);

            for (int i = 0; i < Ngrid; i++)
                rho_s[i] = alpha_s * rho[i] + (1.0 - alpha_s) * pr[i];

            solve_yukawa_bvp(kap, mu, rho, Mt_new);
            solve_yukawa_bvp(kap, mu, rho_s, Ms_new);

            double delta = 0;
            for (int i = 0; i < Ngrid; i++) {
                double mt = relax * Mt_new[i] + (1.0 - relax) * M_t[i];
                double ms = relax * Ms_new[i] + (1.0 - relax) * M_s[i];
                delta += (mt - M_t[i]) * (mt - M_t[i]) * dr;
                delta += (ms - M_s[i]) * (ms - M_s[i]) * dr;
                M_t[i] = mt;
                M_s[i] = ms;
            }
            compute_metric();
            if (sqrt(delta) < 1e-14 && iter > 5) break;
        }

        double omega2 = solve_eigenvalue();
        reconstruct_u();

        double A_min = 1, C_min = 1, C_max = 1;
        for (int i = 0; i < Ngrid; i++) {
            if (AA[i] < A_min) A_min = AA[i];
            if (CC[i] < C_min) C_min = CC[i];
            if (CC[i] > C_max) C_max = CC[i];
        }

        double r2a = 0, nrm = 0;
        for (int i = 0; i < Ngrid; i++) {
            r2a += rg[i] * rg[i] * u_em[i] * u_em[i] * dr;
            nrm += u_em[i] * u_em[i] * dr;
        }
        double R_eff = sqrt(r2a / nrm);
        double omega = (omega2 > 0) ? sqrt(omega2) : 0;

        printf("%14.6e  %14.8f  %14.8f  %14.8f  %14.8f  %14.6e  %14.4f  %14.6f\n",
               kap, omega2, A_min, C_min, C_max, 1.0 - A_min, omega * 197.3, R_eff);
        fflush(stdout);
    }
}

/* ==== Alpha_s scan ==== */
static void run_alphascan(int na)
{
    printf("# V10 TOV Geon — alpha_s Scan\n");
    printf("# kappa=%.6e  mu=%.4f  l=%d  N=%d  Rmax=%.4f  nlmetric=%d\n",
           kappa, mu, ell, Ngrid, Rmax, nl_metric);
    printf("# Scanning alpha_s from 0 to 1 in %d steps\n#\n", na);
    printf("# %14s  %14s  %14s  %14s  %14s  %14s  %14s  %14s\n",
           "alpha_s", "omega2", "A_min", "C_min", "C_max", "well_A", "E_MeV", "R_eff");

    for (int ia = 0; ia <= na; ia++) {
        double as = (double)ia / (double)na;

        for (int i = 0; i < Ngrid; i++) { M_t[i] = 0; M_s[i] = 0; }
        compute_metric();

        int max_iter = 300;
        double relax = 0.3;
        double rho[NMAX], pr[NMAX], rho_s[NMAX];
        double Mt_new[NMAX], Ms_new[NMAX];

        for (int iter = 0; iter < max_iter; iter++) {
            double omega2 = solve_eigenvalue();
            reconstruct_u();
            compute_rho_pr(omega2, rho, pr);

            for (int i = 0; i < Ngrid; i++)
                rho_s[i] = as * rho[i] + (1.0 - as) * pr[i];

            solve_yukawa_bvp(kappa, mu, rho, Mt_new);
            solve_yukawa_bvp(kappa, mu, rho_s, Ms_new);

            double delta = 0;
            for (int i = 0; i < Ngrid; i++) {
                double mt = relax * Mt_new[i] + (1.0 - relax) * M_t[i];
                double ms = relax * Ms_new[i] + (1.0 - relax) * M_s[i];
                delta += (mt - M_t[i]) * (mt - M_t[i]) * dr;
                delta += (ms - M_s[i]) * (ms - M_s[i]) * dr;
                M_t[i] = mt;
                M_s[i] = ms;
            }
            compute_metric();
            if (sqrt(delta) < 1e-14 && iter > 5) break;
        }

        double omega2 = solve_eigenvalue();
        reconstruct_u();

        double A_min = 1, C_min = 1, C_max = 1;
        for (int i = 0; i < Ngrid; i++) {
            if (AA[i] < A_min) A_min = AA[i];
            if (CC[i] < C_min) C_min = CC[i];
            if (CC[i] > C_max) C_max = CC[i];
        }

        double r2a = 0, nrm = 0;
        for (int i = 0; i < Ngrid; i++) {
            r2a += rg[i] * rg[i] * u_em[i] * u_em[i] * dr;
            nrm += u_em[i] * u_em[i] * dr;
        }
        double R_eff = sqrt(r2a / nrm);
        double omega = (omega2 > 0) ? sqrt(omega2) : 0;

        printf("%14.6e  %14.8f  %14.8f  %14.8f  %14.8f  %14.6e  %14.4f  %14.6f\n",
               as, omega2, A_min, C_min, C_max, 1.0 - A_min, omega * 197.3, R_eff);
        fflush(stdout);
    }
}

/* ==== Rmax scan (for Q9 test) ==== */
static void run_rmaxscan(double rmin, double rmax, int nr)
{
    printf("# V10 TOV Geon — Rmax Scan (Q9 Test)\n");
    printf("# kappa=%.6e  mu=%.4f  alpha_s=%.4f  l=%d  nlmetric=%d  N=%d\n",
           kappa, mu, alpha_s, ell, nl_metric, Ngrid);
    printf("# Scanning Rmax from %.2f to %.2f in %d steps\n#\n", rmin, rmax, nr);
    printf("# %14s  %14s  %14s  %14s  %14s  %14s  %14s  %14s\n",
           "Rmax", "omega2", "A_min", "C_min", "C_max", "well_A", "E_MeV", "R_eff");

    for (int ir = 0; ir <= nr; ir++) {
        Rmax = rmin + (rmax - rmin) * (double)ir / (double)nr;
        setup_grid();

        for (int i = 0; i < Ngrid; i++) { M_t[i] = 0; M_s[i] = 0; }
        compute_metric();

        int max_iter = 300;
        double relax = 0.3;
        double rho[NMAX], pr[NMAX], rho_s[NMAX];
        double Mt_new[NMAX], Ms_new[NMAX];

        for (int iter = 0; iter < max_iter; iter++) {
            double omega2 = solve_eigenvalue();
            reconstruct_u();
            compute_rho_pr(omega2, rho, pr);

            for (int i = 0; i < Ngrid; i++)
                rho_s[i] = alpha_s * rho[i] + (1.0 - alpha_s) * pr[i];

            solve_yukawa_bvp(kappa, mu, rho, Mt_new);
            solve_yukawa_bvp(kappa, mu, rho_s, Ms_new);

            double delta = 0;
            for (int i = 0; i < Ngrid; i++) {
                double mt = relax * Mt_new[i] + (1.0 - relax) * M_t[i];
                double ms = relax * Ms_new[i] + (1.0 - relax) * M_s[i];
                delta += (mt - M_t[i]) * (mt - M_t[i]) * dr;
                delta += (ms - M_s[i]) * (ms - M_s[i]) * dr;
                M_t[i] = mt;
                M_s[i] = ms;
            }
            compute_metric();
            if (sqrt(delta) < 1e-14 && iter > 5) break;
        }

        double omega2 = solve_eigenvalue();
        reconstruct_u();

        double A_min = 1, C_min = 1, C_max = 1;
        for (int i = 0; i < Ngrid; i++) {
            if (AA[i] < A_min) A_min = AA[i];
            if (CC[i] < C_min) C_min = CC[i];
            if (CC[i] > C_max) C_max = CC[i];
        }

        double r2a = 0, nrm = 0;
        for (int i = 0; i < Ngrid; i++) {
            r2a += rg[i] * rg[i] * u_em[i] * u_em[i] * dr;
            nrm += u_em[i] * u_em[i] * dr;
        }
        double R_eff = sqrt(r2a / nrm);
        double omega = (omega2 > 0) ? sqrt(omega2) : 0;

        printf("%14.6f  %14.8f  %14.8f  %14.8f  %14.8f  %14.6e  %14.4f  %14.6f\n",
               Rmax, omega2, A_min, C_min, C_max, 1.0 - A_min, omega * 197.3, R_eff);
        fflush(stdout);
    }
}

/* ==== WKB Lifetime ==== */
static void run_lifetime(void)
{
    int n = Ngrid;
    double rho[NMAX], pr[NMAX], rho_s[NMAX];
    double Mt_new[NMAX], Ms_new[NMAX];
    double ll1 = (double)(ell * (ell + 1));

    printf("# V10 TOV Geon — WKB Lifetime Analysis\n");
    printf("# kappa=%.6e  mu=%.4f  alpha_s=%.4f  l=%d  N=%d  Rmax=%.4f  nlmetric=%d\n",
           kappa, mu, alpha_s, ell, Ngrid, Rmax, nl_metric);
    printf("#\n");

    /* Hartree iteration */
    if (seed_amp > 0.0) {
        for (int i = 0; i < n; i++) {
            M_t[i] = seed_amp * rg[i] * exp(-mu * rg[i]);
            M_s[i] = alpha_s * M_t[i];
        }
    } else {
        for (int i = 0; i < n; i++) { M_t[i] = 0; M_s[i] = 0; }
    }
    compute_metric();

    double relax = 0.3;
    int max_iter = 300;
    double omega2_flat = 0;

    for (int iter = 0; iter < max_iter; iter++) {
        double omega2_i = solve_eigenvalue();
        if (iter == 0) omega2_flat = omega2_i;
        reconstruct_u();
        compute_rho_pr(omega2_i, rho, pr);

        for (int i = 0; i < n; i++)
            rho_s[i] = alpha_s * rho[i] + (1.0 - alpha_s) * pr[i];

        solve_yukawa_bvp(kappa, mu, rho, Mt_new);
        solve_yukawa_bvp(kappa, mu, rho_s, Ms_new);

        double delta = 0;
        for (int i = 0; i < n; i++) {
            double mt = relax * Mt_new[i] + (1.0 - relax) * M_t[i];
            double ms = relax * Ms_new[i] + (1.0 - relax) * M_s[i];
            delta += (mt - M_t[i]) * (mt - M_t[i]) * dr;
            delta += (ms - M_s[i]) * (ms - M_s[i]) * dr;
            M_t[i] = mt;
            M_s[i] = ms;
        }
        compute_metric();
        if (sqrt(delta) < 1e-14 && iter > 5) break;
    }

    double omega2 = solve_eigenvalue();
    reconstruct_u();

    if (omega2 <= 0) {
        printf("# ERROR: omega^2 = %.6e <= 0\n", omega2);
        return;
    }
    double omega = sqrt(omega2);

    /* Hartree report */
    double A_min_h = 1, C_max_h = 1;
    for (int i = 0; i < n; i++) {
        if (AA[i] < A_min_h) A_min_h = AA[i];
        if (CC[i] > C_max_h) C_max_h = CC[i];
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
    printf("# A_min        = %.10f\n", A_min_h);
    printf("# C_max        = %.10f\n", C_max_h);
    printf("# well_depth_A = %.6f\n", 1.0 - A_min_h);
    printf("# R_eff        = %.6f\n", R_eff);

    /* Extended grid beyond box */
    double r_out_flat = sqrt(ll1 / omega2);
    double R_ext = r_out_flat + 2.0;
    if (R_ext < 3.0 * Rmax) R_ext = 3.0 * Rmax;

    int N_ext = NEXT;
    double dr_e = R_ext / (N_ext - 1);
    double Mt_Rmax = M_t[n - 1];

    static double r_e[NEXT], A_e[NEXT], V_e[NEXT];

    for (int i = 0; i < N_ext; i++) {
        r_e[i] = (i + 0.5) * dr_e;
        double Mt_i;
        if (r_e[i] <= rg[n - 1]) {
            double j_f = r_e[i] / dr - 0.5;
            int j = (int)j_f;
            if (j < 0) j = 0;
            if (j >= n - 1) j = n - 2;
            double t = j_f - j;
            if (t < 0) t = 0;
            if (t > 1) t = 1;
            Mt_i = (1.0 - t) * M_t[j] + t * M_t[j + 1];
        } else {
            Mt_i = Mt_Rmax * exp(-mu * (r_e[i] - Rmax));
        }
        double x_t = 2.0 * Mt_i / r_e[i];
        A_e[i] = nl_func(x_t);
        V_e[i] = A_e[i] * ll1 / (r_e[i] * r_e[i]);
    }

    printf("#\n# === Extended Domain ===\n");
    printf("# R_ext        = %.4f\n", R_ext);
    printf("# r_out (flat) = %.6f\n", r_out_flat);

    /* Barrier structure */
    int phase = 1;
    int i_ws = -1, i_we = -1, i_bs = -1, i_be = -1;
    double V_bmax = 0;

    for (int i = 0; i < N_ext; i++) {
        double V = V_e[i];
        if (phase == 1 && V < omega2) { phase = 2; i_ws = i; }
        else if (phase == 2 && V >= omega2) { phase = 3; i_we = i - 1; i_bs = i; V_bmax = V; }
        else if (phase == 3) {
            if (V > V_bmax) V_bmax = V;
            if (V < omega2) { phase = 4; i_be = i - 1; }
        }
    }
    if (phase == 3) i_be = N_ext - 1;
    if (phase == 2) i_we = N_ext - 1;
    int has_barrier = (i_bs >= 0 && i_be > i_bs);

    printf("#\n# === Barrier Structure ===\n");
    if (has_barrier) {
        printf("# Well:    r = %.4f to %.4f\n", r_e[i_ws], r_e[i_we]);
        printf("# Barrier: r = %.4f to %.4f  (width %.4f)\n",
               r_e[i_bs], r_e[i_be], r_e[i_be] - r_e[i_bs]);
        printf("# V_max/omega^2 = %.4f\n", V_bmax / omega2);
    } else {
        printf("# NO BARRIER found\n");
    }

    /* WKB tunneling */
    double S_wkb = 0.0;
    if (has_barrier) {
        for (int i = i_bs; i <= i_be; i++) {
            if (V_e[i] > omega2 && A_e[i] > 1e-10) {
                /* In two-potential metric, the tortoise factor uses both A and C.
                 * Beyond the box, C→1, so factor is just 1/√A ≈ 1/A^{1/2}.
                 * More precisely: dr* = dr/√(AC). Beyond box, C=1. */
                S_wkb += sqrt(V_e[i] - omega2) / A_e[i] * dr_e;
            }
        }
    }

    double P_tunnel = exp(-2.0 * S_wkb);
    double Gamma_c = omega * P_tunnel;
    double tau_c = (Gamma_c > 1e-300) ? 1.0 / Gamma_c : 1e300;

    double hbar_c = 197.3;
    double t_code = 3.336e-24;
    double tau_phys = tau_c * t_code;
    double tau_hadron = 1e-23;

    printf("#\n# === WKB LIFETIME RESULTS ===\n");
    printf("# WKB action S:    %.6f\n", S_wkb);
    printf("# exp(-2S):        %.6e\n", P_tunnel);
    printf("# omega:           %.4f MeV\n", omega * hbar_c);
    printf("# Gamma:           %.4e MeV\n", Gamma_c * hbar_c);
    printf("# tau:             %.4e s\n", tau_phys);
    printf("# tau/tau_had:     %.4e\n", tau_phys / tau_hadron);
    printf("# Q = omega/Gamma: %.4e\n",
           Gamma_c > 1e-300 ? omega / Gamma_c : 1e300);

    if (!has_barrier) {
        printf("#\n# *** NO BARRIER — mode leaks freely ***\n");
    } else if (tau_phys > 100 * tau_hadron) {
        printf("#\n# VERY LONG-LIVED: tau = %.0f x tau_hadron\n",
               tau_phys / tau_hadron);
    }
}

/* ==== Multi-mode Hartree ==== */
static void run_multimode(void)
{
    int n = Ngrid;
    double rho_total[NMAX], pr_total[NMAX], rho_s[NMAX];
    double rho_mode[NMAX], pr_mode[NMAX];
    double Mt_new[NMAX], Ms_new[NMAX];
    double omega2_m[MAX_MODES], omega2_flat_m[MAX_MODES];
    int save_ell = ell;

    printf("# V10 TOV Geon — Multi-Mode Hartree\n");
    printf("# kappa=%.6e  mu=%.4f  alpha_s=%.4f  nlmetric=%d  N=%d  Rmax=%.4f\n",
           kappa, mu, alpha_s, nl_metric, Ngrid, Rmax);
    printf("# Modes (%d): ", n_multi_modes);
    double total_occ = 0;
    for (int m = 0; m < n_multi_modes; m++) {
        printf("l=%d(x%.0f) ", mode_l_arr[m], mode_occ_arr[m]);
        total_occ += mode_occ_arr[m];
    }
    printf("\n# Total occupancy: %.0f\n#\n", total_occ);

    for (int i = 0; i < n; i++) { M_t[i] = 0; M_s[i] = 0; }
    compute_metric();

    double relax = 0.3;
    int max_iter = 500;
    for (int m = 0; m < n_multi_modes; m++) omega2_flat_m[m] = 0;

    printf("# %4s", "iter");
    for (int m = 0; m < n_multi_modes; m++) printf("    omega2_m%d  ", m);
    printf("  %14s  %14s  %14s  %14s\n", "A_min", "C_max", "P_min", "delta_m");

    for (int iter = 0; iter < max_iter; iter++) {
        for (int i = 0; i < n; i++) { rho_total[i] = 0; pr_total[i] = 0; }

        for (int m = 0; m < n_multi_modes; m++) {
            ell = mode_l_arr[m];
            omega2_m[m] = solve_eigenvalue();
            if (iter == 0) omega2_flat_m[m] = omega2_m[m];
            reconstruct_u();
            compute_rho_pr(omega2_m[m], rho_mode, pr_mode);
            for (int i = 0; i < n; i++) {
                rho_total[i] += mode_occ_arr[m] * rho_mode[i];
                pr_total[i]  += mode_occ_arr[m] * pr_mode[i];
            }
        }

        for (int i = 0; i < n; i++)
            rho_s[i] = alpha_s * rho_total[i] + (1.0 - alpha_s) * pr_total[i];

        solve_yukawa_bvp(kappa, mu, rho_total, Mt_new);
        solve_yukawa_bvp(kappa, mu, rho_s, Ms_new);

        double delta = 0;
        for (int i = 0; i < n; i++) {
            double mt = relax * Mt_new[i] + (1.0 - relax) * M_t[i];
            double ms = relax * Ms_new[i] + (1.0 - relax) * M_s[i];
            delta += (mt - M_t[i]) * (mt - M_t[i]) * dr;
            delta += (ms - M_s[i]) * (ms - M_s[i]) * dr;
            M_t[i] = mt;
            M_s[i] = ms;
        }
        compute_metric();
        delta = sqrt(delta);

        double A_min = 1, C_max_it = 1, P_min = 1;
        for (int i = 0; i < n; i++) {
            if (AA[i] < A_min) A_min = AA[i];
            if (CC[i] > C_max_it) C_max_it = CC[i];
            if (PP[i] < P_min) P_min = PP[i];
        }

        printf("%4d", iter);
        for (int m = 0; m < n_multi_modes; m++) printf("  %14.6f", omega2_m[m]);
        printf("  %14.8f  %14.8f  %14.8f  %14.6e\n", A_min, C_max_it, P_min, delta);
        fflush(stdout);

        if (delta < 1e-14 && iter > 5) {
            printf("# CONVERGED at iter %d\n", iter + 1);
            break;
        }
    }

    /* Final results */
    printf("#\n# === FINAL RESULTS ===\n");

    double A_min = 1, C_min_f = 1, C_max_f = 1;
    for (int i = 0; i < n; i++) {
        if (AA[i] < A_min) A_min = AA[i];
        if (CC[i] < C_min_f) C_min_f = CC[i];
        if (CC[i] > C_max_f) C_max_f = CC[i];
    }
    printf("# A_min      = %.10f\n", A_min);
    printf("# C_min      = %.10f\n", C_min_f);
    printf("# C_max      = %.10f\n", C_max_f);
    printf("# well_A     = %.6e\n", 1.0 - A_min);

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
    printf("# Ratio E/M_p:      %.4f\n", E_total / 938.3);

    ell = save_ell;
}

/* ==== Parse mode specification ==== */
static void parse_modes(const char *spec)
{
    const char *p = spec;
    n_multi_modes = 0;
    while (*p && n_multi_modes < MAX_MODES) {
        int l_val = (int)strtol(p, (char **)&p, 10);
        double n_val = 1.0;
        if (*p == ':') { p++; n_val = strtod(p, (char **)&p); }
        if (l_val < 1) l_val = 1;
        mode_l_arr[n_multi_modes] = l_val;
        mode_occ_arr[n_multi_modes] = n_val;
        n_multi_modes++;
        if (*p == ',') p++;
    }
}

/* ==== Main ==== */
int main(int argc, char *argv[])
{
    int mode_hartree = 0, mode_scan = 0, mode_alphascan = 0;
    int mode_lifetime = 0, mode_multimode = 0, mode_rmaxscan = 0;
    double kappa_max = 100.0;
    int nk = 20, na = 20, nr = 20;
    double rmin = 0.3, rmax_scan = 5.0;

    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "-hartree"))    { mode_hartree = 1; continue; }
        if (!strcmp(argv[i], "-scan"))       { mode_scan = 1; continue; }
        if (!strcmp(argv[i], "-alphascan"))  { mode_alphascan = 1; continue; }
        if (!strcmp(argv[i], "-rmaxscan"))   { mode_rmaxscan = 1; continue; }
        if (!strcmp(argv[i], "-lifetime"))   { mode_lifetime = 1; continue; }
        if (!strcmp(argv[i], "-multimode"))  { mode_multimode = 1; continue; }
        if (i + 1 < argc) {
            if (!strcmp(argv[i], "-kappa"))    { kappa = atof(argv[++i]); continue; }
            if (!strcmp(argv[i], "-kmax"))     { kappa_max = atof(argv[++i]); continue; }
            if (!strcmp(argv[i], "-nk"))       { nk = atoi(argv[++i]); continue; }
            if (!strcmp(argv[i], "-mu"))       { mu = atof(argv[++i]); continue; }
            if (!strcmp(argv[i], "-alpha"))    { alpha_s = atof(argv[++i]); continue; }
            if (!strcmp(argv[i], "-na"))       { na = atoi(argv[++i]); continue; }
            if (!strcmp(argv[i], "-nr"))       { nr = atoi(argv[++i]); continue; }
            if (!strcmp(argv[i], "-rmin"))     { rmin = atof(argv[++i]); continue; }
            if (!strcmp(argv[i], "-rmax"))     { rmax_scan = atof(argv[++i]); continue; }
            if (!strcmp(argv[i], "-seed"))     { seed_amp = atof(argv[++i]); continue; }
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

    if (mode_alphascan)
        run_alphascan(na);
    else if (mode_rmaxscan)
        run_rmaxscan(rmin, rmax_scan, nr);
    else if (mode_scan)
        run_scan(kappa_max, nk);
    else if (mode_lifetime)
        run_lifetime();
    else if (mode_multimode) {
        if (n_multi_modes == 0) {
            n_multi_modes = 1;
            mode_l_arr[0] = ell;
            mode_occ_arr[0] = 3;
        }
        run_multimode();
    }
    else if (mode_hartree)
        run_hartree();
    else {
        printf("Usage: tov_geon -hartree|-scan|-alphascan|-rmaxscan|-lifetime|-multimode [options]\n");
        printf("  -hartree         Two-potential Hartree solver\n");
        printf("  -scan            Scan kappa at fixed alpha_s\n");
        printf("  -alphascan       Scan alpha_s at fixed kappa\n");
        printf("  -rmaxscan        Scan Rmax (Q9 test)\n");
        printf("  -lifetime        WKB lifetime of quasi-bound geon\n");
        printf("  -multimode       Multi-mode Hartree (use with -modes)\n");
        printf("\nParameters:\n");
        printf("  -kappa <K>       Coupling constant [%.2e]\n", kappa);
        printf("  -mu <M>          Yukawa mass [%.4f]\n", mu);
        printf("  -alpha <A>       Spatial source mixing: 0=Fierz-Pauli, 1=V9 [%.2f]\n", alpha_s);
        printf("  -kmax <K>        Max kappa for scan [%.2e]\n", kappa_max);
        printf("  -nk <N>          Scan steps [%d]\n", nk);
        printf("  -na <N>          Alpha scan steps [%d]\n", na);
        printf("  -nr <N>          Rmax scan steps [%d]\n", nr);
        printf("  -rmin <R>        Min Rmax for scan [%.2f]\n", rmin);
        printf("  -rmax <R>        Max Rmax for scan [%.2f]\n", rmax_scan);
        printf("  -seed <A>        Seed well amplitude [%.2f]\n", seed_amp);
        printf("  -nlmetric <N>    Metric: 0=linear, 1=exp, 2=Pade [%d]\n", nl_metric);
        printf("  -modes <spec>    l1:n1,l2:n2,... (e.g. 1:3)\n");
        printf("  -l <L>           Angular momentum [%d]\n", ell);
        printf("  -N <N>           Grid points [%d]\n", Ngrid);
        printf("  -Rmax <R>        Domain size [%.1f]\n", Rmax);
        printf("\nPhysical reference (f_2 meson):\n");
        printf("  mu_f2 = 6.47 fm^{-1},  kappa_strong ~ 35\n");
        printf("\nalpha_s values:\n");
        printf("  1.0 = V9 recovery (A=C, single potential)\n");
        printf("  0.0 = Fierz-Pauli (spatial sourced by p_r)\n");
    }

    return 0;
}
