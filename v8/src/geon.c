/*
 * geon.c — Spherical Geon Hartree Solver
 *
 * Self-gravitating EM standing wave in spherically symmetric GR.
 *
 * Metric:
 *   ds^2 = -f(r) dt^2 + dr^2/f(r) + r^2 dOmega^2
 *   f(r) = 1 - 2*m(r)/r
 *   (delta = 0 approximation for weak gravity)
 *
 * EM field: TM_l axial mode, radial function u(r), frequency omega.
 * Wave equation (Regge-Wheeler, in r coordinates):
 *   f^2 u'' + f f' u' + [omega^2 - f l(l+1)/r^2] u = 0
 *
 * Symmetrized variable: v = sqrt(f) * u  removes first derivative:
 *   -f^2 v'' + W(r) v = omega^2 v
 *   W(r) = f l(l+1)/r^2 + f f''/2 - f'^2/4
 *
 * Einstein equation (time-averaged T_00):
 *   m'(r) = kappa * r^2 * rho_EM(r)
 *
 * where rho_EM = N_l/(8 pi r^2) * [omega^2 u^2 / f + f (u')^2]
 * with N_l = l(l+1)/(2l+1) (angular average factor).
 *
 * Hartree iteration:
 *   1. From m(r), compute f(r) and W(r)
 *   2. Solve eigenvalue problem for v(r), omega^2
 *   3. Reconstruct u = v/sqrt(f), compute rho_EM
 *   4. Integrate m'(r) = kappa r^2 rho_EM
 *   5. Under-relax m, repeat
 *
 * Key parameter: kappa = 4*pi*G (gravitational coupling).
 *   kappa = 0: flat space, delocalized box mode
 *   kappa > 0: gravitational well, mode shifts down
 *   kappa_cr: critical coupling for self-consistent geon
 *
 * Also supports Born-Infeld nonlinearity in the EM sector.
 *
 * Build: cd src && make
 * Usage: geon -hartree -kappa <value> [options]
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
static double kappa   = 0.0;      /* 4*pi*G (gravitational coupling) */
static double b_BI    = 0.0;      /* Born-Infeld field strength (0 = Maxwell) */
static int    ell     = 1;        /* angular momentum (>= 1 for EM) */

static int    Ngrid   = 8001;
static double Rmax    = 20.0;

/* Mode normalization: integral u^2 dr = 1 */
/* Total EM energy: E_EM = omega^2 (in these units) */

/* ---- Grid ---- */
static double dr;
static double rg[NMAX];

/* ---- Metric ---- */
static double mass[NMAX];   /* m(r) */
static double ff[NMAX];     /* f(r) = 1 - 2m/r */
static double ffp[NMAX];    /* f'(r) */
static double ffpp[NMAX];   /* f''(r) */

/* ---- EM mode ---- */
static double u_em[NMAX];   /* radial function u(r) */
static double v_em[NMAX];   /* symmetrized: v = sqrt(f) * u */

/* ==== Grid ==== */
static void setup_grid(void)
{
    if (Ngrid > NMAX) Ngrid = NMAX;
    dr = Rmax / (Ngrid - 1);
    for (int i = 0; i < Ngrid; i++)
        rg[i] = (i + 0.5) * dr;  /* half-cell offset avoids r=0 */
}

/* ==== Metric from mass function ==== */
static void compute_metric(void)
{
    int n = Ngrid;

    for (int i = 0; i < n; i++) {
        ff[i] = 1.0 - 2.0 * mass[i] / rg[i];
        if (ff[i] < 0.01) ff[i] = 0.01;  /* prevent horizon */
    }

    /* f' by centered differences */
    ffp[0] = (ff[1] - ff[0]) / dr;
    for (int i = 1; i < n - 1; i++)
        ffp[i] = (ff[i+1] - ff[i-1]) / (2.0 * dr);
    ffp[n-1] = (ff[n-1] - ff[n-2]) / dr;

    /* f'' */
    ffpp[0] = (ff[1] - 2.0*ff[0] + ff[0]) / (dr*dr);  /* ghost: f[-1]=f[0] */
    for (int i = 1; i < n - 1; i++)
        ffpp[i] = (ff[i+1] - 2.0*ff[i] + ff[i-1]) / (dr*dr);
    ffpp[n-1] = (ff[n-1] - 2.0*ff[n-2] + ff[n-3]) / (dr*dr);
}

/* ==== Effective potential W(r) for symmetrized variable v ==== */
/*
 * From the Regge-Wheeler equation:
 *   -f^2 v'' + W(r) v = omega^2 v
 * where W(r) = f l(l+1)/r^2 + f f''/2 - f'^2/4
 */
static double W_eff(int i)
{
    double r = rg[i];
    double f = ff[i];
    double fp = ffp[i];
    double fpp = ffpp[i];
    return f * (double)(ell*(ell+1)) / (r*r) + f*fpp/2.0 - fp*fp/4.0;
}

/* ==== Eigenvalue solver: implicit imaginary time ==== */
/*
 * Solves: -f^2 v'' + W(r) v = omega^2 v
 * using implicit scheme: (I + dtau H) v^{n+1} = v^n
 * where H v = -f^2 v'' + W v (in discretized form).
 *
 * Returns omega^2.
 */
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
        Bt[0] = 1.0 + dtau * (3.0 * f2_0 / dr2 + W_eff(0));
        Ct[0] = -dtau * f2_0 / dr2;
        Dt[0] = v_em[0];

        for (int i = 1; i < n - 1; i++) {
            double f2 = ff[i] * ff[i];
            At[i] = -dtau * f2 / dr2;
            Bt[i] = 1.0 + dtau * (2.0 * f2 / dr2 + W_eff(i));
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
        double Hv = -f2 * vdd + W_eff(0) * v_em[0];
        hv += v_em[0] * Hv * dr;
        vv += v_em[0] * v_em[0] * dr;
    }
    for (int i = 1; i < n - 1; i++) {
        double f2 = ff[i] * ff[i];
        double vdd = (v_em[i+1] - 2.0*v_em[i] + v_em[i-1]) / dr2;
        double Hv = -f2 * vdd + W_eff(i) * v_em[i];
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
/*
 * Time-averaged EM energy density (spherically averaged):
 *   rho_EM(r) = N_l / (8 pi r^2) * [omega^2 u^2 / f + f (u')^2]
 * where N_l = l(l+1)/(2l+1).
 *
 * With Born-Infeld (b_BI > 0): saturate rho_EM at b_BI^2.
 */
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

        /* Maxwell energy density */
        double rho_M = N_l / (8.0 * M_PI * r * r) *
                       (omega2 * u * u / f + f * up[i] * up[i]);

        if (b_BI > 0.0 && rho_M > b_BI * b_BI) {
            /* BI saturation: rho_BI = b^2 [sqrt(1 + rho_M/b^2) - 1] */
            /* For rho_M >> b^2: rho_BI ~ b * sqrt(rho_M) */
            rho[i] = b_BI * b_BI * (sqrt(1.0 + rho_M / (b_BI*b_BI)) - 1.0);
        } else {
            rho[i] = rho_M;
        }
    }
}

/* ==== Integrate Einstein equation: m'(r) = kappa * r^2 * rho ==== */
static void integrate_einstein(double *rho, double *m_new)
{
    int n = Ngrid;
    m_new[0] = 0.0;  /* regularity at origin */
    for (int i = 1; i < n; i++) {
        /* Simpson-like: m_new[i] = m_new[i-1] + kappa * r^2 * rho * dr */
        double r_mid = 0.5 * (rg[i] + rg[i-1]);
        double rho_mid = 0.5 * (rho[i] + rho[i-1]);
        m_new[i] = m_new[i-1] + kappa * r_mid * r_mid * rho_mid * dr;
    }
}

/* ==== Hartree iteration ==== */
static void run_hartree(void)
{
    int n = Ngrid;
    double rho[NMAX], m_new[NMAX];

    printf("# V8 Geon Hartree Solver\n");
    printf("# kappa=%.6e  l=%d  b_BI=%.4f  N=%d  Rmax=%.1f\n",
           kappa, ell, b_BI, Ngrid, Rmax);
    printf("#\n");

    /* Start with flat space */
    for (int i = 0; i < n; i++) mass[i] = 0.0;
    compute_metric();

    double relax = 0.3;
    int max_iter = 200;

    printf("# %4s  %14s  %14s  %14s  %14s  %14s  %14s\n",
           "iter", "omega2", "delta_omega2", "M_total",
           "Phi_min", "f_min", "delta_m");

    double omega2_flat = 0;

    for (int iter = 0; iter < max_iter; iter++) {
        /* Step 1: Solve eigenvalue problem */
        double omega2 = solve_eigenvalue();

        if (iter == 0) omega2_flat = omega2;

        /* Step 2: Reconstruct u, compute rho */
        reconstruct_u();
        compute_rho(omega2, rho);

        /* Step 3: Integrate Einstein equations */
        integrate_einstein(rho, m_new);

        /* Total mass */
        double M_total = m_new[n-1];

        /* Minimum f(r) */
        double f_min = 1.0;
        double Phi_min = 0.0;
        for (int i = 0; i < n; i++) {
            double f_test = 1.0 - 2.0 * m_new[i] / rg[i];
            if (f_test < f_min) {
                f_min = f_test;
                Phi_min = -mass[i] / rg[i];  /* Newtonian approx */
            }
        }

        /* Under-relax mass function */
        double delta = 0;
        for (int i = 0; i < n; i++) {
            double m_upd = relax * m_new[i] + (1.0 - relax) * mass[i];
            delta += (m_upd - mass[i]) * (m_upd - mass[i]) * dr;
            mass[i] = m_upd;
        }
        delta = sqrt(delta);

        compute_metric();

        printf("%4d  %14.8f  %14.6e  %14.6e  %14.6e  %14.8f  %14.6e\n",
               iter, omega2, omega2 - omega2_flat, M_total,
               Phi_min, f_min, delta);
        fflush(stdout);

        if (delta < 1e-14 && iter > 5) {
            printf("# CONVERGED at iter %d\n", iter + 1);
            break;
        }
    }

    /* Final eigenvalue */
    double omega2_final = solve_eigenvalue();
    reconstruct_u();

    printf("#\n# === FINAL RESULT ===\n");
    printf("# omega^2 (flat)  = %.10f\n", omega2_flat);
    printf("# omega^2 (geon)  = %.10f\n", omega2_final);
    printf("# delta(omega^2)  = %.6e\n", omega2_final - omega2_flat);
    printf("# frac shift      = %.6e\n",
           (omega2_final - omega2_flat) / omega2_flat);
    printf("# M_total         = %.6e\n", mass[Ngrid-1]);

    /* Compute well depth */
    double f_min_final = 1.0;
    int i_min = 0;
    for (int i = 0; i < n; i++) {
        if (ff[i] < f_min_final) {
            f_min_final = ff[i];
            i_min = i;
        }
    }
    printf("# f_min           = %.10f  (at r=%.4f)\n", f_min_final, rg[i_min]);
    printf("# well depth      = %.6e  (1-f_min)\n", 1.0 - f_min_final);
    printf("# Phi_min/c^2     = %.6e\n", -(1.0 - f_min_final) / 2.0);

    /* Print profile */
    printf("#\n# Profile:\n");
    printf("# %10s  %14s  %14s  %14s  %14s  %14s\n",
           "r", "u(r)", "m(r)", "f(r)", "W_eff(r)", "rho(r)");

    double rho_final[NMAX];
    compute_rho(omega2_final, rho_final);

    int step = (n > 400) ? n / 400 : 1;
    for (int i = 0; i < n; i += step) {
        printf("%12.6f  %14.6e  %14.6e  %14.8f  %14.6e  %14.6e\n",
               rg[i], u_em[i], mass[i], ff[i], W_eff(i), rho_final[i]);
    }
}

/* ==== Coupling scan: kappa from 0 to kappa_max ==== */
static void run_scan(double kappa_max, int nk)
{
    printf("# V8 Geon Coupling Scan\n");
    printf("# l=%d  b_BI=%.4f  N=%d  Rmax=%.1f\n",
           ell, b_BI, Ngrid, Rmax);
    printf("# Scanning kappa from 0 to %.4e in %d steps\n#\n", kappa_max, nk);
    printf("# %14s  %14s  %14s  %14s  %14s  %14s\n",
           "kappa", "omega2", "delta_omega2", "M_total",
           "well_depth", "frac_shift");

    double omega2_flat = 0;

    for (int ik = 0; ik <= nk; ik++) {
        kappa = kappa_max * (double)ik / (double)nk;

        /* Reset to flat */
        for (int i = 0; i < Ngrid; i++) mass[i] = 0.0;
        compute_metric();

        /* Hartree iteration (quiet) */
        int max_iter = 200;
        double relax = 0.3;
        double rho[NMAX], m_new[NMAX];

        for (int iter = 0; iter < max_iter; iter++) {
            double omega2 = solve_eigenvalue();
            if (ik == 0 && iter == 0) omega2_flat = omega2;

            reconstruct_u();
            compute_rho(omega2, rho);
            integrate_einstein(rho, m_new);

            double delta = 0;
            for (int i = 0; i < Ngrid; i++) {
                double m_upd = relax * m_new[i] + (1.0 - relax) * mass[i];
                delta += (m_upd - mass[i]) * (m_upd - mass[i]) * dr;
                mass[i] = m_upd;
            }
            compute_metric();
            delta = sqrt(delta);
            if (delta < 1e-14 && iter > 5) break;
        }

        double omega2 = solve_eigenvalue();
        double M_total = mass[Ngrid-1];

        double f_min = 1.0;
        for (int i = 0; i < Ngrid; i++)
            if (ff[i] < f_min) f_min = ff[i];

        printf("%14.6e  %14.8f  %14.6e  %14.6e  %14.6e  %14.6e\n",
               kappa, omega2, omega2 - omega2_flat, M_total,
               1.0 - f_min,
               (omega2_flat > 1e-15) ? (omega2 - omega2_flat)/omega2_flat : 0.0);
        fflush(stdout);
    }
}

/* ==== First-order perturbation theory ==== */
/*
 * Compute the gravitational eigenvalue shift analytically:
 *   delta(omega^2) = -kappa * <u | delta_V | u> / <u | u>
 * where delta_V comes from the change in the effective potential
 * due to the self-gravitating mass distribution.
 *
 * This is the DEFINITIVE test of positive vs negative feedback.
 */
static void run_perturbative(void)
{
    int n = Ngrid;

    printf("# V8 Geon — First-Order Perturbation Theory\n");
    printf("# l=%d  N=%d  Rmax=%.1f\n", ell, Ngrid, Rmax);
    printf("#\n");

    /* Step 1: Flat-space eigenmode */
    for (int i = 0; i < n; i++) mass[i] = 0.0;
    compute_metric();  /* f=1 everywhere */

    double omega2_0 = solve_eigenvalue();
    reconstruct_u();

    printf("# Flat-space eigenvalue: omega^2_0 = %.10f\n", omega2_0);
    printf("# Expected: (x_{l,1}/Rmax)^2 = %.10f  (x_{1,1}=4.4934)\n",
           (4.4934 / Rmax) * (4.4934 / Rmax));

    /* Step 2: Compute self-gravitating mass function at unit kappa */
    double rho[NMAX], m_self[NMAX];
    double kappa_unit = 1.0;
    kappa = kappa_unit;
    compute_rho(omega2_0, rho);
    integrate_einstein(rho, m_self);

    double M_unit = m_self[n-1];
    printf("# Total mass at kappa=1: M = %.10e\n", M_unit);

    /* Step 3: Compute the potential perturbation */
    /* delta_f(r) = -2*m(r)/r, so delta_W ~ delta_f * l(l+1)/r^2 */
    /* V_eff = f * l(l+1)/r^2; delta_V = delta_f * l(l+1)/r^2 = -2m/r * l(l+1)/r^2 */
    /* First-order shift: delta(omega^2) = <v|delta_W|v>/<v|v> */

    /* For the symmetrized variable: delta_W includes also f''/2 and f'^2 terms.
     * At first order in m/r (weak field), the dominant term is:
     *   delta_W ~ -2m(r)/r * l(l+1)/r^2
     * The f''/2 and f'^2/4 corrections are second-order in m/r. */

    double delta_omega2 = 0.0, vv = 0.0;
    for (int i = 0; i < n; i++) {
        double r = rg[i];
        double delta_W = -2.0 * m_self[i] / r * (double)(ell*(ell+1)) / (r*r);
        delta_omega2 += v_em[i] * delta_W * v_em[i] * dr;
        vv += v_em[i] * v_em[i] * dr;
    }
    delta_omega2 /= vv;

    printf("#\n# First-order eigenvalue shift (at kappa=1):\n");
    printf("#   delta(omega^2) = %.10e\n", delta_omega2);
    printf("#   frac shift     = %.10e\n", delta_omega2 / omega2_0);

    if (delta_omega2 < 0) {
        printf("#\n# POSITIVE FEEDBACK CONFIRMED: gravity LOWERS the eigenvalue.\n");
        printf("# The EM mode is attracted to the gravitational well.\n");
    } else {
        printf("#\n# NEGATIVE FEEDBACK: gravity RAISES the eigenvalue.\n");
        printf("# The EM mode is repelled from the gravitational well.\n");
    }

    /* Step 4: Scaling analysis */
    printf("#\n# === Scaling Analysis ===\n");
    printf("# At kappa=1 (Planck units):\n");
    printf("#   Well depth     = %.6e\n", 2.0 * M_unit / rg[n/4]);
    printf("#   frac shift     = %.6e per unit kappa\n",
           delta_omega2 / omega2_0);

    /* Physical values (nuclear scale) */
    double G_N = 6.674e-11;          /* m^3/(kg s^2) */
    double hbar = 1.055e-34;         /* J s */
    double c = 3.0e8;                /* m/s */
    double M_Pl_GeV = 1.221e19;      /* GeV */
    double r_proton_fm = 0.88;       /* fm */
    double m_proton_GeV = 0.938;     /* GeV */
    double GeV_to_invfm = 5.068;     /* 1/fm per GeV (hbar*c = 0.1973 GeV*fm) */

    /* kappa_phys in natural units (hbar=c=1): kappa = 4*pi/M_Pl^2 */
    double kappa_phys = 4.0 * M_PI / (M_Pl_GeV * M_Pl_GeV);

    /* The eigenvalue shift scales linearly with kappa:
     * delta(omega^2)/omega^2 = (delta/omega2_0) * kappa
     * For binding: need |delta|/omega2_0 * kappa ~ 1 */
    double kappa_cr = -omega2_0 / delta_omega2;  /* critical kappa */

    printf("#\n# Critical coupling for self-consistent geon:\n");
    printf("#   kappa_cr = %.6e  (need |shift| ~ omega^2)\n", kappa_cr);
    printf("#\n# Physical gravitational coupling:\n");
    printf("#   kappa_phys = 4*pi/M_Pl^2 = %.6e GeV^{-2}\n", kappa_phys);
    printf("#\n# COUPLING GAP:\n");
    printf("#   kappa_cr / kappa_phys = %.6e\n", kappa_cr / kappa_phys);
    printf("#   This is the factor by which G must be enhanced\n");
    printf("#   for gravity to self-trap light at this scale.\n");

    /* Physical interpretation */
    /* In our code units, Rmax sets the scale. If Rmax = 1 fm: */
    double scale_fm = Rmax;  /* treat code length as fm */
    double scale_GeV = GeV_to_invfm / scale_fm;  /* GeV per code length */
    double omega_GeV = sqrt(omega2_0) * scale_GeV;
    double E_geon_GeV = omega_GeV;  /* single photon mode energy */

    printf("#\n# If code length = 1 fm:\n");
    printf("#   omega = %.4f GeV = %.1f MeV\n", omega_GeV, omega_GeV*1000);
    printf("#   E_photon = %.4f GeV (single mode)\n", E_geon_GeV);
    printf("#   R_S(proton) = 2*G*m_p/c^2 = %.2e fm\n",
           2.0 * G_N * m_proton_GeV * 1.783e-27 / (c*c) * 1e15);
    (void)hbar; (void)r_proton_fm;

    /* Energy density diagnostics */
    printf("#\n# Energy density profile (flat space):\n");
    printf("# %10s  %14s  %14s  %14s  %14s\n",
           "r", "u(r)", "rho(r)", "m_self(r)", "delta_V(r)");
    int step = (n > 200) ? n / 200 : 1;
    for (int i = 0; i < n; i += step) {
        double r = rg[i];
        double delta_V = -2.0 * m_self[i] / r * (double)(ell*(ell+1))/(r*r);
        printf("%12.6f  %14.6e  %14.6e  %14.6e  %14.6e\n",
               rg[i], u_em[i], rho[i], m_self[i], delta_V);
    }
}

/* ==== Main ==== */
int main(int argc, char *argv[])
{
    int mode_hartree = 0;
    int mode_scan = 0;
    int mode_perturb = 0;
    double kappa_max = 1.0;
    int nk = 20;

    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "-hartree"))   { mode_hartree = 1; continue; }
        if (!strcmp(argv[i], "-scan"))      { mode_scan = 1; continue; }
        if (!strcmp(argv[i], "-perturb"))   { mode_perturb = 1; continue; }
        if (i + 1 < argc) {
            if (!strcmp(argv[i], "-kappa"))    { kappa = atof(argv[++i]); continue; }
            if (!strcmp(argv[i], "-kmax"))     { kappa_max = atof(argv[++i]); continue; }
            if (!strcmp(argv[i], "-nk"))       { nk = atoi(argv[++i]); continue; }
            if (!strcmp(argv[i], "-bBI"))      { b_BI = atof(argv[++i]); continue; }
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
    else if (mode_scan)
        run_scan(kappa_max, nk);
    else if (mode_hartree)
        run_hartree();
    else {
        printf("Usage: geon -perturb|-hartree|-scan [options]\n");
        printf("  -perturb         First-order perturbation theory\n");
        printf("  -hartree         Full Hartree self-consistent solve\n");
        printf("  -scan            Scan kappa from 0 to kmax\n");
        printf("  -kappa <K>       Gravitational coupling [%.2e]\n", kappa);
        printf("  -kmax <K>        Max kappa for scan [%.2e]\n", kappa_max);
        printf("  -nk <N>          Number of scan steps [%d]\n", nk);
        printf("  -bBI <b>         Born-Infeld field strength [%.2f]\n", b_BI);
        printf("  -l <L>           Angular momentum [%d]\n", ell);
        printf("  -N <N>           Grid points [%d]\n", Ngrid);
        printf("  -Rmax <R>        Domain size [%.1f]\n", Rmax);
    }

    return 0;
}
