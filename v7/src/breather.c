/*
 * breather.c — V7 Vector Breather in Disformal Metric
 *
 * Self-trapped massless EM wave in Bekenstein disformal metric.
 * Key differences from selftrap.c (scalar proxy):
 *   - NO mass term m^2 (massless vector field)
 *   - Source: alpha_s * F^2 with Born-Infeld saturation
 *   - Mass emerges from self-trapping energy, not Lagrangian mass
 *
 * F^2 = B^2 - E^2 for the TM_l radial mode:
 *   B^2_eff ~ (u'-u/r)^2/r^2   (spatial gradient, cos^2 phase)
 *   E^2_eff ~ udot^2/r^2        (time derivative, sin^2 phase)
 *   Time-averaged: <F^2> = [(u'-u/r)^2 - omega^2*u^2]/(2r^2)
 *
 * BI saturation: source = alpha_s*F^2/sqrt(1+|F^2|/(2*b^2))
 * prevents the parametric resonance that crashes the scalar code.
 *
 * Uses u = r*f for the TM_l radial mode function f(r).
 *
 * Modes:
 *   -hartree    Hartree self-consistent (with optional m^2 bootstrap)
 *   -quench     Mass quench: converge at m>0, track as m->0
 *   -coupled    Coupled BI time evolution
 *   -fsign +/-1 Sign of F^2 source (+1 = standard dilaton, -1 = flipped)
 *
 * Build: cd src && make breather
 */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define NMAX 8001
#define NV   4
#define A_MIN 0.01

/* ---- Parameters ---- */
static double alpha_c = 0.3;    /* conformal coupling */
static double beta_d  = 1.0;    /* disformal coupling */
static double alpha_s = 5.0;    /* source coupling */
static double b_BI    = 10.0;   /* Born-Infeld field strength */
static double m2_boot = 1.0;    /* bootstrap mass^2 for Hartree */

static int    Ngrid = 4001;
static double Rmax  = 40.0;
static double Tmax  = 200.0;
static double cfl   = 0.4;

static double amp_p = 5.0;
static double wid_p = 3.0;
static double k_p   = 3.0;
static double r0_p  = 15.0;
static double sig0_f = 5.0;     /* initial sigma well depth */
static double wid_f  = 4.0;     /* initial sigma well width */

static int ell = 1;
static int f_sign = -1;         /* F^2 source sign: -1 = attractive, +1 = standard */

static int mode_hartree = 0;
static int mode_quench  = 0;
static int mode_coupled = 0;
static int do_snap      = 0;

/* ---- Grid ---- */
static double dr;
static double rg[NMAX];

/* ---- State ---- */
static double Y[NV][NMAX];
static double K1[NV][NMAX], K2[NV][NMAX], K3[NV][NMAX], K4[NV][NMAX];
static double Ytmp[NV][NMAX];

/* ---- Metric ---- */
static double MET_v2[NMAX], MET_Ct[NMAX], MET_Cr[NMAX], MET_Crp[NMAX];

/* ---- Sponge ---- */
static double sponge[NMAX];

/* ==== Grid ==== */
static void setup_grid(void)
{
    if (Ngrid > NMAX) Ngrid = NMAX;
    dr = Rmax / Ngrid;
    for (int i = 0; i < Ngrid; i++)
        rg[i] = (i + 0.5) * dr;

    double r_sp = 0.85 * Rmax;
    for (int i = 0; i < Ngrid; i++) {
        if (rg[i] > r_sp) {
            double x = (rg[i] - r_sp) / (Rmax - r_sp);
            sponge[i] = 10.0 * x * x;
        } else {
            sponge[i] = 0.0;
        }
    }
}

/* ==== Metric from w = r*sigma ==== */
static void compute_metric_from_w(double *w)
{
    int n = Ngrid;
    double sig[NMAX], sigp[NMAX], wp[NMAX];

    wp[0] = (w[1] + w[0]) / (2.0 * dr);
    for (int i = 1; i < n - 1; i++)
        wp[i] = (w[i+1] - w[i-1]) / (2.0 * dr);
    wp[n-1] = (w[n-1] - w[n-2]) / dr;

    for (int i = 0; i < n; i++) {
        sig[i]  = w[i] / rg[i];
        sigp[i] = (wp[i] - sig[i]) / rg[i];
    }

    for (int i = 0; i < n; i++) {
        double A = exp(2.0 * alpha_c * sig[i]);
        if (A < A_MIN) A = A_MIN;
        double X = sigp[i] * sigp[i];
        double D = A + beta_d * X;
        MET_Ct[i] = sqrt(A * D);
        MET_Cr[i] = A * sqrt(A / D);
        MET_v2[i] = A / D;
    }

    MET_Crp[0] = (MET_Cr[1] - MET_Cr[0]) / dr;
    for (int i = 1; i < n - 1; i++)
        MET_Crp[i] = (MET_Cr[i+1] - MET_Cr[i-1]) / (2.0 * dr);
    MET_Crp[n-1] = (MET_Cr[n-1] - MET_Cr[n-2]) / dr;
}

/* ==== BI-saturated source ==== */
static double bi_source(double F2)
{
    double X = fabs(F2) / (2.0 * b_BI * b_BI);
    return (double)f_sign * alpha_s * F2 / sqrt(1.0 + X);
}

/* ==== Thomas algorithm ==== */
static void thomas_solve(int n, double *a, double *b, double *c, double *d, double *x)
{
    double cp[NMAX], dp[NMAX];
    cp[0] = c[0] / b[0];
    dp[0] = d[0] / b[0];
    for (int i = 1; i < n; i++) {
        double den = b[i] - a[i] * cp[i-1];
        cp[i] = c[i] / den;
        dp[i] = (d[i] - a[i] * dp[i-1]) / den;
    }
    x[n-1] = dp[n-1];
    for (int i = n-2; i >= 0; i--)
        x[i] = dp[i] - cp[i] * x[i+1];
}

/* ==== Sigma BVP: w'' = -r*source (massless, Neumann BC) ==== */
static void solve_sigma_bvp(double *source_r, double *w_out)
{
    int n = Ngrid;
    double A[NMAX], B[NMAX], C[NMAX], D[NMAX];
    double dr2 = dr * dr;

    for (int i = 0; i < n; i++) {
        A[i] = (i > 0) ? 1.0/dr2 : 0.0;
        C[i] = (i < n-1) ? 1.0/dr2 : 0.0;
        B[i] = -2.0/dr2;
        D[i] = -rg[i] * source_r[i];
    }
    B[0] = -3.0/dr2;   /* ghost: w[-1] = -w[0] */

    /* Massless sigma: Neumann w'(Rmax) = 0 */
    A[n-1] = -1.0; B[n-1] = 1.0; C[n-1] = 0.0; D[n-1] = 0.0;

    thomas_solve(n, A, B, C, D, w_out);
}

/* ==== Eigenstate finder with variable v^2(r) ==== */
/*
 * Solves: -v^2(r)*u'' + l(l+1)*u/r^2 + m^2*A(r)*u = omega^2 * u
 * using implicit imaginary time evolution.
 *
 * For vector (m^2=0): only centrifugal + metric modification.
 * For scalar (m^2>0): adds mass well (bootstrap).
 *
 * If u_init != NULL, uses it as initial guess (for quench tracking).
 */
static double find_eigenstate(double *v2, double m2, double *A_arr,
                              double *u_init, double *u_out)
{
    int n = Ngrid;
    double u[NMAX], u_new[NMAX];
    double At[NMAX], Bt[NMAX], Ct[NMAX], Dt[NMAX];
    double dr2 = dr * dr;

    /* Initial guess */
    if (u_init) {
        memcpy(u, u_init, n * sizeof(double));
    } else {
        for (int i = 0; i < n; i++)
            u[i] = pow(rg[i], ell) * exp(-rg[i]*rg[i] / 25.0);
    }
    u[n-1] = 0.0;

    double norm2 = 0;
    for (int i = 0; i < n; i++) norm2 += u[i]*u[i]*dr;
    if (norm2 > 1e-30) {
        double inv = 1.0/sqrt(norm2);
        for (int i = 0; i < n; i++) u[i] *= inv;
    }

    double dtau = 0.5;
    int niter = 5000;

    for (int iter = 0; iter < niter; iter++) {
        /* Operator: H = -v^2 d^2/dr^2 + l(l+1)/r^2 + m^2*A
         * Tridiagonal: (1 + dtau*H) u_new = u */

        /* i=0: ghost u[-1]=-u[0] → u''=(u1-3u0)/dr2 */
        At[0] = 0.0;
        Bt[0] = 1.0 + dtau*(3.0*v2[0]/dr2 + (double)(ell*(ell+1))/(rg[0]*rg[0])
                            + m2*A_arr[0]);
        Ct[0] = -dtau * v2[0] / dr2;
        Dt[0] = u[0];

        for (int i = 1; i < n-1; i++) {
            At[i] = -dtau * v2[i] / dr2;
            Bt[i] = 1.0 + dtau*(2.0*v2[i]/dr2 + (double)(ell*(ell+1))/(rg[i]*rg[i])
                                + m2*A_arr[i]);
            Ct[i] = -dtau * v2[i] / dr2;
            Dt[i] = u[i];
        }

        At[n-1] = 0.0; Bt[n-1] = 1.0; Ct[n-1] = 0.0; Dt[n-1] = 0.0;

        thomas_solve(n, At, Bt, Ct, Dt, u_new);

        int bad = 0;
        for (int i = 0; i < n; i++)
            if (isnan(u_new[i]) || isinf(u_new[i])) { bad = 1; break; }
        if (bad) { dtau *= 0.5; continue; }

        norm2 = 0.0;
        for (int i = 0; i < n; i++) norm2 += u_new[i]*u_new[i]*dr;
        if (norm2 < 1e-30) break;
        double inv = 1.0/sqrt(norm2);
        for (int i = 0; i < n; i++) u[i] = u_new[i]*inv;
    }

    /* Eigenvalue: omega^2 = <u|H|u>/<u|u> */
    double hu = 0.0, uu = 0.0;
    {
        double udd = (u[1] - 3.0*u[0])/dr2;
        double Hu = -v2[0]*udd + (double)(ell*(ell+1))*u[0]/(rg[0]*rg[0])
                    + m2*A_arr[0]*u[0];
        hu += u[0]*Hu*dr;
        uu += u[0]*u[0]*dr;
    }
    for (int i = 1; i < n-1; i++) {
        double udd = (u[i+1] - 2.0*u[i] + u[i-1])/dr2;
        double Hu = -v2[i]*udd + (double)(ell*(ell+1))*u[i]/(rg[i]*rg[i])
                    + m2*A_arr[i]*u[i];
        hu += u[i]*Hu*dr;
        uu += u[i]*u[i]*dr;
    }

    memcpy(u_out, u, n * sizeof(double));
    return hu / uu;
}

/* ==== Compute derivatives of u ==== */
static void compute_up(double *u, double *up, int n)
{
    up[0] = (u[1] + u[0]) / (2.0 * dr);
    for (int i = 1; i < n-1; i++)
        up[i] = (u[i+1] - u[i-1]) / (2.0 * dr);
    up[n-1] = (u[n-1] - u[n-2]) / dr;
}

/* ==== Compute time-averaged F^2 source from eigenstate ==== */
/*
 * <F^2> = <B^2-E^2> = [(u'-u/r)^2 - omega^2*u^2] / (2*r^2)
 * with BI saturation.
 *
 * Also computes the source for a "massive scalar" comparison:
 *   source_m2 = -alpha_s * m^2 * psi^2 = -alpha_s * m^2 * u^2/r^2
 */
static void compute_F2_source(double *u, double omega2,
                              double *source_out, double *F2_out)
{
    int n = Ngrid;
    double up[NMAX];
    compute_up(u, up, n);

    for (int i = 0; i < n; i++) {
        double ri = rg[i];
        double q = up[i] - u[i]/ri;
        /* Time-averaged F^2 for standing wave */
        double F2 = 0.5*(q*q - omega2*u[i]*u[i]) / (ri*ri);
        F2_out[i] = F2;
        source_out[i] = bi_source(F2);
    }
}

/* ==== Compute massive-scalar source (for comparison/bootstrap) ==== */
static void compute_m2_source(double *u, double m2, double *source_out)
{
    int n = Ngrid;
    for (int i = 0; i < n; i++) {
        double psi2 = u[i]*u[i] / (rg[i]*rg[i]);
        source_out[i] = -alpha_s * m2 * psi2;
    }
}

/* ==== Hartree solver ==== */
/*
 * source_mode: 0 = F^2 source (vector), 1 = m^2*psi^2 source (scalar)
 * m2: mass^2 for eigenvalue problem (0 for massless vector)
 * verbose: print per-iteration output
 * u_init: initial eigenstate guess (NULL for default)
 * w_out, u_out: converged solution output
 *
 * Returns converged omega^2, or -1 if failed.
 */
static double hartree_solve(int source_mode, double m2, int verbose,
                            double *u_init, double *w_out, double *u_out)
{
    int n = Ngrid;
    double w[NMAX], u_eig[NMAX], source[NMAX], F2_arr[NMAX];
    double v2_arr[NMAX], A_arr[NMAX], w_old[NMAX];

    /* Bootstrap: start with trial eigenstate, compute initial sigma */
    if (u_init) {
        memcpy(u_eig, u_init, n * sizeof(double));
    } else {
        double norm2 = 0;
        for (int i = 0; i < n; i++) {
            u_eig[i] = pow(rg[i], ell) * exp(-rg[i]*rg[i]/(wid_p*wid_p));
            norm2 += u_eig[i]*u_eig[i]*dr;
        }
        double inv = 1.0/sqrt(norm2);
        for (int i = 0; i < n; i++) u_eig[i] *= inv;
    }

    /* Compute initial sigma from trial eigenstate */
    if (source_mode == 0) {
        double omega2_trial = (m2 > 0.01) ? 0.9*m2 : 0.5;
        compute_F2_source(u_eig, omega2_trial, source, F2_arr);
    } else {
        compute_m2_source(u_eig, m2, source);
    }
    solve_sigma_bvp(source, w);

    double relax = 0.3;
    int max_iter = 500;
    double omega2_final = 0;
    (void)omega2_final;

    for (int iter = 0; iter < max_iter; iter++) {
        memcpy(w_old, w, n * sizeof(double));

        /* Compute A(r) from sigma */
        double sig_min = 0, A_min = 1;
        for (int i = 0; i < n; i++) {
            double sig = w[i] / rg[i];
            double A = exp(2.0 * alpha_c * sig);
            if (A < A_MIN) A = A_MIN;
            A_arr[i] = A;
            if (sig < sig_min) sig_min = sig;
            if (A < A_min) A_min = A;
        }

        /* Compute v^2(r) from metric */
        compute_metric_from_w(w);
        for (int i = 0; i < n; i++)
            v2_arr[i] = MET_v2[i];

        /* Find eigenstate */
        double omega2 = find_eigenstate(v2_arr, m2, A_arr,
                                        (iter == 0 && u_init) ? u_init : NULL,
                                        u_eig);

        /* Compute source */
        if (source_mode == 0) {
            compute_F2_source(u_eig, omega2, source, F2_arr);
        } else {
            compute_m2_source(u_eig, m2, source);
        }

        /* Solve sigma BVP */
        double w_new[NMAX];
        solve_sigma_bvp(source, w_new);

        /* Under-relax and check convergence */
        double delta = 0;
        for (int i = 0; i < n; i++) {
            double wn = relax*w_new[i] + (1.0-relax)*w_old[i];
            delta += (wn - w_old[i])*(wn - w_old[i])*dr;
            w[i] = wn;
        }
        delta = sqrt(delta);

        double D_charge = w[n-2];

        /* F^2 diagnostics */
        double F2_core = 0, F2_far = 0;
        if (source_mode == 0) {
            int i_core = (int)(2.0/dr);
            int i_far = (int)(10.0/dr);
            if (i_core < n) F2_core = F2_arr[i_core];
            if (i_far < n) F2_far = F2_arr[i_far];
        }

        if (verbose) {
            printf("%4d  %12.6f  %12.6f  %12.6f  %12.6f  %12.4e",
                   iter, omega2, sig_min, A_min, D_charge, delta);
            if (source_mode == 0)
                printf("  %12.4e  %12.4e", F2_core, F2_far);
            printf("\n");
            fflush(stdout);
        }

        omega2_final = omega2;

        if (delta < 1e-10 && iter > 5) {
            if (verbose) printf("# CONVERGED at iter %d\n", iter+1);
            break;
        }
    }

    /* Copy output */
    if (w_out) memcpy(w_out, w, n * sizeof(double));
    if (u_out) memcpy(u_out, u_eig, n * sizeof(double));

    return omega2_final;
}

/* ==== Hartree mode: single run ==== */
static void run_hartree(void)
{
    int n = Ngrid;

    printf("# V7 Vector Breather — Hartree Solver\n");
    printf("# alpha_c=%.4f  alpha_s=%.4f  b_BI=%.4f  l=%d  f_sign=%d\n",
           alpha_c, alpha_s, b_BI, ell, f_sign);
    printf("# m2_boot=%.4f  N=%d  Rmax=%.1f\n", m2_boot, Ngrid, Rmax);
    printf("#\n");

    /* Phase 1: Bootstrap with massive scalar source */
    printf("# Phase 1: Massive scalar bootstrap (m^2=%.4f)\n", m2_boot);
    printf("# %4s  %12s  %12s  %12s  %12s  %12s\n",
           "iter", "omega2", "sig_min", "A_min", "D_charge", "delta_w");

    double w_boot[NMAX], u_boot[NMAX];
    double omega2_boot = hartree_solve(1, m2_boot, 1, NULL, w_boot, u_boot);

    printf("#\n# Bootstrap: omega2=%.6f, E_bind=%.6f\n",
           omega2_boot, m2_boot - omega2_boot);

    /* Phase 2: Switch to F^2 source, keeping m^2 in eigenvalue eq */
    printf("#\n# Phase 2: Switch to F^2 source (keeping m^2=%.4f in eigenvalue)\n",
           m2_boot);
    printf("# %4s  %12s  %12s  %12s  %12s  %12s  %12s  %12s\n",
           "iter", "omega2", "sig_min", "A_min", "D_charge", "delta_w",
           "F2_core", "F2_far");

    double w_f2[NMAX], u_f2[NMAX];
    double omega2_f2 = hartree_solve(0, m2_boot, 1, u_boot, w_f2, u_f2);

    printf("#\n# F^2 source with m^2: omega2=%.6f\n", omega2_f2);

    /* Phase 3: Remove mass (m^2=0), pure vector */
    printf("#\n# Phase 3: Remove mass (m^2=0), pure vector mode\n");
    printf("# %4s  %12s  %12s  %12s  %12s  %12s  %12s  %12s\n",
           "iter", "omega2", "sig_min", "A_min", "D_charge", "delta_w",
           "F2_core", "F2_far");

    double w_vec[NMAX], u_vec[NMAX];
    double omega2_vec = hartree_solve(0, 0.0, 1, u_f2, w_vec, u_vec);

    printf("#\n# === FINAL RESULT (m^2=0, pure vector) ===\n");
    printf("# omega^2     = %.6f\n", omega2_vec);
    printf("# D_charge    = %.6f\n", w_vec[n-2]);

    /* Compute final F^2 profile */
    double source[NMAX], F2_arr[NMAX];
    compute_F2_source(u_vec, omega2_vec, source, F2_arr);

    double F2_int = 0;
    for (int i = 0; i < n; i++)
        F2_int += F2_arr[i] * 4.0*M_PI*rg[i]*rg[i] * dr;
    printf("# integral F^2 = %.6e\n", F2_int);

    /* Print final profile */
    printf("#\n# Final profile:\n");
    printf("# %10s  %12s  %12s  %12s  %12s  %12s\n",
           "r", "sigma", "A", "f_eig", "v2", "F2");
    for (int i = 0; i < n; i += (n > 200 ? n/200 : 1)) {
        double sig = w_vec[i] / rg[i];
        double A = exp(2.0 * alpha_c * sig);
        if (A < A_MIN) A = A_MIN;
        printf("%10.4f  %12.6f  %12.6f  %12.6f  %12.6f  %12.4e\n",
               rg[i], sig, A, u_vec[i]/rg[i], MET_v2[i], F2_arr[i]);
    }
}

/* ==== Mass quench: track solution as m^2 decreases ==== */
static void run_quench(void)
{
    int n = Ngrid;

    printf("# V7 Vector Breather — Mass Quench\n");
    printf("# alpha_c=%.4f  alpha_s=%.4f  b_BI=%.4f  l=%d  f_sign=%d\n",
           alpha_c, alpha_s, b_BI, ell, f_sign);
    printf("# N=%d  Rmax=%.1f\n#\n", Ngrid, Rmax);

    /* Step 1: Converge with massive scalar source at m^2=1 */
    printf("# Step 1: Massive bootstrap (m^2=1.0, scalar source)\n");
    double w_curr[NMAX], u_curr[NMAX];
    double omega2 = hartree_solve(1, 1.0, 0, NULL, w_curr, u_curr);
    printf("# Bootstrap: omega2=%.6f, D=%.6f\n\n", omega2, w_curr[n-2]);

    /* Step 2: Switch to F^2 source at m^2=1 */
    printf("# Step 2: Switch to F^2 source at m^2=1\n");
    omega2 = hartree_solve(0, 1.0, 0, u_curr, w_curr, u_curr);
    printf("# F^2 at m^2=1: omega2=%.6f, D=%.6f\n\n", omega2, w_curr[n-2]);

    /* Step 3: Quench m^2 toward 0 */
    printf("# Step 3: Mass quench\n");
    printf("# %10s  %12s  %12s  %12s  %12s  %12s\n",
           "m^2", "omega^2", "sig_min", "A_min", "D_charge", "F2_int");

    double m2_vals[] = {1.0, 0.5, 0.2, 0.1, 0.05, 0.02, 0.01, 0.005, 0.002, 0.001, 0.0};
    int n_m2 = sizeof(m2_vals) / sizeof(m2_vals[0]);

    for (int im = 0; im < n_m2; im++) {
        double m2 = m2_vals[im];
        double w_new[NMAX], u_new[NMAX];

        omega2 = hartree_solve(0, m2, 0, u_curr, w_new, u_new);

        /* Compute diagnostics */
        double sig_min = 0, A_min = 1;
        for (int i = 0; i < n; i++) {
            double sig = w_new[i] / rg[i];
            double A = exp(2.0*alpha_c*sig);
            if (A < A_MIN) A = A_MIN;
            if (sig < sig_min) sig_min = sig;
            if (A < A_min) A_min = A;
        }

        double source[NMAX], F2_arr[NMAX];
        compute_F2_source(u_new, omega2, source, F2_arr);
        double F2_int = 0;
        for (int i = 0; i < n; i++)
            F2_int += F2_arr[i] * 4.0*M_PI*rg[i]*rg[i] * dr;

        printf("%10.4f  %12.6f  %12.6f  %12.6f  %12.6f  %12.4e\n",
               m2, omega2, sig_min, A_min, w_new[n-2], F2_int);
        fflush(stdout);

        /* Use as starting point for next m^2 */
        memcpy(w_curr, w_new, n*sizeof(double));
        memcpy(u_curr, u_new, n*sizeof(double));
    }
}

/* ==== Coupled time evolution with BI source ==== */
static void compute_metric_Y(double st[NV][NMAX])
{
    compute_metric_from_w(st[2]);
}

static void compute_rhs(double st[NV][NMAX], double rhs[NV][NMAX])
{
    double *u = st[0], *ud = st[1], *w = st[2], *wd = st[3];
    int n = Ngrid;

    compute_metric_Y(st);

    double up[NMAX], upp[NMAX], wpp[NMAX];

    /* u derivatives with ghost u[-1]=-u[0] */
    up[0] = (u[1] + u[0]) / (2.0*dr);
    for (int i = 1; i < n-1; i++)
        up[i] = (u[i+1] - u[i-1]) / (2.0*dr);
    up[n-1] = (u[n-1] - u[n-2]) / dr;

    upp[0] = (u[1] - 3.0*u[0]) / (dr*dr);
    for (int i = 1; i < n-1; i++)
        upp[i] = (u[i+1] - 2.0*u[i] + u[i-1]) / (dr*dr);
    upp[n-1] = (u[n-1] - 2.0*u[n-2] + u[n-3]) / (dr*dr);

    /* w'' with ghost w[-1]=-w[0] */
    wpp[0] = (w[1] - 3.0*w[0]) / (dr*dr);
    for (int i = 1; i < n-1; i++)
        wpp[i] = (w[i+1] - 2.0*w[i] + w[i-1]) / (dr*dr);
    wpp[n-1] = (w[n-1] - 2.0*w[n-2] + w[n-3]) / (dr*dr);

    for (int i = 0; i < n; i++) {
        double ri = rg[i];

        /* Wave: NO mass term
         * u_tt = v^2*u'' + (C_r'/C_t)(u'-u/r) - l(l+1)u/r^2 */
        double qu = up[i] - u[i]/ri;
        double udd = MET_v2[i] * upp[i]
                   + (MET_Crp[i] / MET_Ct[i]) * qu
                   - (double)(ell*(ell+1)) * u[i] / (ri*ri);

        /* Instantaneous F^2 = (u'-u/r)^2/r^2 - udot^2/r^2 */
        double F2 = (qu*qu - ud[i]*ud[i]) / (ri*ri);
        double src = bi_source(F2);

        /* Sigma: w_tt = w'' + r*source */
        double wdd = wpp[i] + ri * src;

        /* Sponge */
        udd -= sponge[i] * ud[i];
        wdd -= sponge[i] * wd[i];

        rhs[0][i] = ud[i];
        rhs[1][i] = udd;
        rhs[2][i] = wd[i];
        rhs[3][i] = wdd;
    }

    /* Sommerfeld BC */
    int nl = n-1;
    rhs[1][nl] = -(ud[nl] - ud[nl>0 ? nl-1 : 0])/dr - sponge[nl]*ud[nl];
    rhs[3][nl] = -(wd[nl] - wd[nl>0 ? nl-1 : 0])/dr - sponge[nl]*wd[nl];
}

static void rk4_step(double dt)
{
    int n = Ngrid;

    compute_rhs(Y, K1);
    for (int v = 0; v < NV; v++)
        for (int i = 0; i < n; i++)
            Ytmp[v][i] = Y[v][i] + 0.5*dt*K1[v][i];

    compute_rhs(Ytmp, K2);
    for (int v = 0; v < NV; v++)
        for (int i = 0; i < n; i++)
            Ytmp[v][i] = Y[v][i] + 0.5*dt*K2[v][i];

    compute_rhs(Ytmp, K3);
    for (int v = 0; v < NV; v++)
        for (int i = 0; i < n; i++)
            Ytmp[v][i] = Y[v][i] + dt*K3[v][i];

    compute_rhs(Ytmp, K4);
    for (int v = 0; v < NV; v++)
        for (int i = 0; i < n; i++)
            Y[v][i] += (dt/6.0) * (K1[v][i] + 2*K2[v][i]
                                  + 2*K3[v][i] + K4[v][i]);
}

static double energy_wave(double rmin, double rmax_e)
{
    /* Massless: H = 1/2 * int [ud^2 + v^2*u'^2 + l(l+1)u^2/r^2] dr */
    double *u = Y[0], *ud = Y[1];
    double E = 0.0;
    for (int i = 0; i < Ngrid; i++) {
        if (rg[i] < rmin || rg[i] > rmax_e) continue;
        double up_i;
        if (i == 0) up_i = (u[1]+u[0])/(2.0*dr);
        else if (i == Ngrid-1) up_i = (u[i]-u[i-1])/dr;
        else up_i = (u[i+1]-u[i-1])/(2.0*dr);
        double ri = rg[i];
        double KE = ud[i]*ud[i];
        double GE = MET_v2[i]*up_i*up_i;
        double CE = (double)(ell*(ell+1))*u[i]*u[i]/(ri*ri);
        E += (KE + GE + CE) * dr;
    }
    return 0.5*E;
}

static void run_coupled(void)
{
    double dt = cfl * dr;
    int nsteps = (int)(Tmax/dt + 0.5);
    int nout = nsteps/100; if (nout<1) nout=1;

    /* Initialize: standing wave + pre-formed sigma well */
    for (int i = 0; i < Ngrid; i++) {
        double ri = rg[i];
        Y[0][i] = amp_p * pow(ri, ell) * exp(-0.5*ri*ri/(wid_p*wid_p));
        Y[1][i] = 0.0;
        double sig = -fabs(sig0_f) * exp(-0.5*ri*ri/(wid_f*wid_f));
        Y[2][i] = ri * sig;
        Y[3][i] = 0.0;
    }

    printf("# V7 Vector Breather — Coupled BI Evolution\n");
    printf("# alpha_c=%.4f  beta_d=%.4f  alpha_s=%.4f  b_BI=%.4f  l=%d  f_sign=%d\n",
           alpha_c, beta_d, alpha_s, b_BI, ell, f_sign);
    printf("# N=%d  Rmax=%.1f  Tmax=%.1f  dt=%.6f\n", Ngrid, Rmax, Tmax, dt);
    printf("#\n");

    compute_metric_Y(Y);
    double E0 = energy_wave(0.0, Rmax);
    printf("# Initial E = %.6f\n#\n", E0);
    printf("# %8s  %10s  %10s  %10s  %10s  %10s\n",
           "t", "E_total", "E_core", "sig_min", "f_max", "F2_max");

    for (int step = 0; step <= nsteps; step++) {
        double t = step * dt;

        if (step % nout == 0 || step == nsteps) {
            compute_metric_Y(Y);
            double Et = energy_wave(0.0, Rmax);
            double Ec = energy_wave(0.0, 10.0);

            double sig_min = 0, f_max = 0, F2_max = 0;
            double *u = Y[0], *ud = Y[1];
            double up[NMAX];
            compute_up(u, up, Ngrid);

            for (int i = 0; i < Ngrid; i++) {
                double si = Y[2][i]/rg[i];
                if (si < sig_min) sig_min = si;
                double fi = fabs(u[i]/rg[i]);
                if (fi > f_max) f_max = fi;
                double qu = up[i] - u[i]/rg[i];
                double F2 = (qu*qu - ud[i]*ud[i])/(rg[i]*rg[i]);
                if (fabs(F2) > F2_max) F2_max = fabs(F2);
            }

            printf("%9.3f  %10.4f  %10.4f  %10.6f  %10.6f  %10.4f\n",
                   t, Et, Ec, sig_min, f_max, F2_max);
            fflush(stdout);
        }

        if (step % (nout*10) == 0 && step > 0) {
            if (isnan(Y[0][Ngrid/2]) || isnan(Y[2][Ngrid/2])) {
                fprintf(stderr, "NaN at step %d, t=%.4f\n", step, step*dt);
                break;
            }
        }

        if (step < nsteps) rk4_step(dt);
    }

    double Ef = energy_wave(0.0, 10.0);
    double frac = (E0 > 1e-15) ? Ef/E0 : 0;
    printf("#\n# RESULT: core retention = %.1f%% %s\n",
           100*frac, frac > 0.1 ? "TRAPPED" : "DISPERSED");
}

/* ==== Main ==== */
int main(int argc, char *argv[])
{
    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "-hartree"))  { mode_hartree = 1; continue; }
        if (!strcmp(argv[i], "-quench"))   { mode_quench = 1; continue; }
        if (!strcmp(argv[i], "-coupled"))  { mode_coupled = 1; continue; }
        if (!strcmp(argv[i], "-snap"))     { do_snap = 1; continue; }
        if (i + 1 < argc) {
            if (!strcmp(argv[i], "-alpha"))  { alpha_c = atof(argv[++i]); continue; }
            if (!strcmp(argv[i], "-beta"))   { beta_d  = atof(argv[++i]); continue; }
            if (!strcmp(argv[i], "-alphas")) { alpha_s = atof(argv[++i]); continue; }
            if (!strcmp(argv[i], "-bBI"))    { b_BI    = atof(argv[++i]); continue; }
            if (!strcmp(argv[i], "-m2"))     { m2_boot = atof(argv[++i]); continue; }
            if (!strcmp(argv[i], "-amp"))    { amp_p   = atof(argv[++i]); continue; }
            if (!strcmp(argv[i], "-k"))      { k_p     = atof(argv[++i]); continue; }
            if (!strcmp(argv[i], "-r0"))     { r0_p    = atof(argv[++i]); continue; }
            if (!strcmp(argv[i], "-w"))      { wid_p   = atof(argv[++i]); continue; }
            if (!strcmp(argv[i], "-sig0"))   { sig0_f  = atof(argv[++i]); continue; }
            if (!strcmp(argv[i], "-wsig"))   { wid_f   = atof(argv[++i]); continue; }
            if (!strcmp(argv[i], "-N"))      { Ngrid   = atoi(argv[++i]); continue; }
            if (!strcmp(argv[i], "-Rmax"))   { Rmax    = atof(argv[++i]); continue; }
            if (!strcmp(argv[i], "-Tmax"))   { Tmax    = atof(argv[++i]); continue; }
            if (!strcmp(argv[i], "-l"))      { ell     = atoi(argv[++i]); continue; }
            if (!strcmp(argv[i], "-fsign"))  { f_sign  = atoi(argv[++i]); continue; }
        }
    }

    if (ell < 1) {
        fprintf(stderr, "ERROR: Vector modes require l >= 1 (got l=%d)\n", ell);
        return 1;
    }

    setup_grid();
    memset(Y, 0, sizeof(Y));

    if (mode_quench)
        run_quench();
    else if (mode_hartree)
        run_hartree();
    else if (mode_coupled)
        run_coupled();
    else {
        printf("Usage: breather -hartree|-quench|-coupled [options]\n");
        printf("  -alpha  A      conformal coupling [%.2f]\n", alpha_c);
        printf("  -alphas S      source coupling    [%.2f]\n", alpha_s);
        printf("  -bBI   B       Born-Infeld b      [%.2f]\n", b_BI);
        printf("  -m2    M       bootstrap mass^2   [%.2f]\n", m2_boot);
        printf("  -l     L       angular momentum   [%d]\n", ell);
        printf("  -fsign +/-1    F^2 source sign    [%d]\n", f_sign);
        printf("  -N     N       grid points        [%d]\n", Ngrid);
        printf("  -Rmax  R       domain size        [%.1f]\n", Rmax);
    }

    return 0;
}
