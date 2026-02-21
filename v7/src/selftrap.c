/*
 * selftrap.c — V7 Disformal Self-Trapping Solver
 *
 * Tests whether a scalar wave self-traps via disformal metric feedback.
 * Wave field psi propagates in effective metric:
 *   g~_uv = e^{2*alpha*sigma} eta_uv + beta * d_u(sigma) d_v(sigma)
 * Inertia field sigma sourced by wave energy density.
 *
 * Uses u = r*psi, w = r*sigma for regularity at r=0.
 *
 * Equations (quasi-static metric approximation):
 *   u_tt = v^2_eff * u_rr + (C_r'/C_t)*(u_r - u/r)
 *   w_tt = w_rr - alpha_s * r * epsilon_psi
 *
 * where v^2_eff = C_r/C_t = A/(A + beta*sigma'^2)
 *   A = exp(2*alpha_c*sigma), clamped >= A_MIN
 *   C_t = sqrt(A*(A+beta*X)), C_r = A^{3/2}/sqrt(A+beta*X), X = sigma'^2
 *   epsilon_psi = (1/2)(C_t*psi_t^2 + C_r*psi_r^2)
 *
 * The disformal term beta*d_u(sigma)*d_v(sigma) slows radial propagation:
 *   v_radial = c / sqrt(1 + beta*sigma'^2 / A) < c
 * while tangential speed remains c. This creates a refractive shell
 * that can trap the wave.
 *
 * Build: cd src && make
 * Usage: selftrap [options]
 *   -fixed          Fixed sigma profile, evolve psi only
 *   -coupled        Coupled evolution (default)
 *   -scan           Scan alpha_s for trapping threshold
 *   -alpha A        Conformal coupling          [0.2]
 *   -beta B         Disformal coupling           [1.0]
 *   -alphas S       Source coupling sigma<-psi    [0.5]
 *   -amp A          Pulse amplitude               [5.0]
 *   -k K            Pulse wave number             [2.0]
 *   -r0 R           Pulse center                  [30.0]
 *   -w W            Pulse width                   [3.0]
 *   -sig0 S         Fixed sigma amplitude         [5.0]
 *   -wsig W         Fixed sigma width             [4.0]
 *   -N N            Grid points                   [4001]
 *   -Rmax R         Domain size                   [80.0]
 *   -Tmax T         Evolution time                [200.0]
 *   -l L            Angular mode number            [1]
 *   -snap           Output snapshot files
 *
 * For l >= 1, centrifugal barrier l(l+1)/r^2 confines the wave near
 * the center. Combined with the disformal refractive shell at r ~ wsig,
 * this creates a potential well that can trap the wave.
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
#define NV   4       /* state variables: u, u_dot, w, w_dot */
#define A_MIN 0.01   /* floor on conformal factor */

/* ---- Parameters ---- */
static double alpha_c = 0.2;
static double beta_d  = 1.0;
static double alpha_s = 0.5;
static double mu_sig  = 0.0;   /* sigma mass (0 = massless) */
static double m_psi   = 0.0;   /* wave field mass */

static int    Ngrid = 4001;
static double Rmax  = 80.0;
static double Tmax  = 200.0;
static double cfl   = 0.4;

static double amp_p = 5.0;
static double r0_p  = 30.0;
static double wid_p = 3.0;
static double k_p   = 2.0;

static double sig0_f = 5.0;
static double wid_f  = 4.0;

static int ell = 1;   /* angular mode number */

static int mode_fixed = 0;
static int mode_scan  = 0;
static int mode_veff    = 0;   /* print V_eff(r) and exit */
static int mode_eigen   = 0;   /* imaginary time evolution to find eigenstate */
static int mode_hartree = 0;   /* Hartree self-consistent solver */
static int do_snap    = 0;
static int init_standing = 0; /* standing wave in well (u_dot=0) */
static int source_sign = -1;  /* -1 = attractive (default), +1 = repulsive (variational) */

/* ---- Grid ---- */
static double dr;
static double rg[NMAX];

/* ---- State: Y[0]=u, Y[1]=u_dot, Y[2]=w, Y[3]=w_dot ---- */
static double Y[NV][NMAX];
static double dY[NV][NMAX];
static double K1[NV][NMAX], K2[NV][NMAX], K3[NV][NMAX], K4[NV][NMAX];
static double Ytmp[NV][NMAX];

/* ---- Metric coefficients ---- */
static double MET_Ct[NMAX], MET_Cr[NMAX], MET_Crp[NMAX];
static double MET_v2[NMAX];   /* v^2_eff = C_r/C_t */

/* ---- Sponge layer ---- */
static double sponge[NMAX];

/* ==================================================================== */
/*  Grid setup                                                           */
/* ==================================================================== */
static void setup_grid(void)
{
    if (Ngrid > NMAX) Ngrid = NMAX;
    dr = Rmax / Ngrid;
    for (int i = 0; i < Ngrid; i++)
        rg[i] = (i + 0.5) * dr;

    /* Sponge in last 15% of domain */
    double r_sp = 0.85 * Rmax;
    double gmax = 10.0;
    for (int i = 0; i < Ngrid; i++) {
        if (rg[i] > r_sp) {
            double x = (rg[i] - r_sp) / (Rmax - r_sp);
            sponge[i] = gmax * x * x;
        } else {
            sponge[i] = 0.0;
        }
    }
}

/* ==================================================================== */
/*  Metric: C_t, C_r, C_r', v^2 from current w array                   */
/* ==================================================================== */
static void compute_metric(double st[NV][NMAX])
{
    double *w = st[2];
    int n = Ngrid;

    /* sigma = w/r, sigma' = (w' - sigma)/r = (w' - w/r)/r */
    double sig[NMAX], sigp[NMAX], wp[NMAX];

    /* w' with ghost point w[-1] = -w[0] (odd symmetry) */
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

        MET_Ct[i] = sqrt(A * D);          /* sqrt(A) * sqrt(A+beta*X) */
        MET_Cr[i] = A * sqrt(A / D);      /* A^{3/2} / sqrt(A+beta*X) */
        MET_v2[i] = A / D;                /* phase velocity squared */
    }

    /* C_r' by central differences */
    MET_Crp[0] = (MET_Cr[1] - MET_Cr[0]) / dr;
    for (int i = 1; i < n - 1; i++)
        MET_Crp[i] = (MET_Cr[i+1] - MET_Cr[i-1]) / (2.0 * dr);
    MET_Crp[n-1] = (MET_Cr[n-1] - MET_Cr[n-2]) / dr;
}

/* ==================================================================== */
/*  RHS computation                                                      */
/* ==================================================================== */
static void compute_rhs(double st[NV][NMAX], double rhs[NV][NMAX])
{
    double *u = st[0], *ud = st[1], *w = st[2], *wd = st[3];
    int n = Ngrid;

    compute_metric(st);

    /* u' with ghost u[-1] = -u[0] */
    double up[NMAX], upp[NMAX];
    up[0] = (u[1] + u[0]) / (2.0 * dr);
    for (int i = 1; i < n - 1; i++)
        up[i] = (u[i+1] - u[i-1]) / (2.0 * dr);
    up[n-1] = (u[n-1] - u[n-2]) / dr;

    /* u'' with ghost u[-1] = -u[0] */
    upp[0] = (u[1] - 3.0 * u[0]) / (dr * dr);
    for (int i = 1; i < n - 1; i++)
        upp[i] = (u[i+1] - 2.0*u[i] + u[i-1]) / (dr * dr);
    upp[n-1] = (u[n-1] - 2.0*u[n-2] + u[n-3]) / (dr * dr);

    /* w'' with ghost w[-1] = -w[0] */
    double wpp[NMAX];
    wpp[0] = (w[1] - 3.0 * w[0]) / (dr * dr);
    for (int i = 1; i < n - 1; i++)
        wpp[i] = (w[i+1] - 2.0*w[i] + w[i-1]) / (dr * dr);
    wpp[n-1] = (w[n-1] - 2.0*w[n-2] + w[n-3]) / (dr * dr);

    for (int i = 0; i < n; i++) {
        double ri = rg[i];

        /* Conformal factor at this point */
        double sig_i = w[i] / ri;
        double A_i = exp(2.0 * alpha_c * sig_i);
        if (A_i < A_MIN) A_i = A_MIN;

        /* ---- Wave equation: u_tt = v^2 u'' + (C_r'/C_t)(u'-u/r) - l(l+1)u/r^2 - m^2*A*u ---- */
        /* Mass potential m^2*A(r) creates well where A < 1 (sigma < 0) */
        double qu = up[i] - u[i] / ri;
        double udd = MET_v2[i] * upp[i]
                   + (MET_Crp[i] / MET_Ct[i]) * qu
                   - (double)(ell * (ell + 1)) * u[i] / (ri * ri)
                   - m_psi * m_psi * A_i * u[i];

        /* ---- Inertia: w_tt = w'' - mu^2*w - alpha_s * A * r * eps_psi ---- */
        double eps = 0.5 * (MET_Ct[i] * ud[i]*ud[i]
                          + MET_Cr[i] * qu*qu) / (ri * ri);
        /* Klein-Gordon restoring force: -mu^2 * w prevents runaway */
        /* source_sign: -1 = energy-density source (attractive, ad hoc)
         *              +1 = variational/conformal source (repulsive, correct) */
        double wdd = wpp[i] - mu_sig*mu_sig * w[i]
                   + source_sign * alpha_s * A_i * ri * eps;

        /* Sponge damping near boundary */
        udd -= sponge[i] * ud[i];
        wdd -= sponge[i] * wd[i];

        rhs[0][i] = ud[i];
        rhs[1][i] = udd;
        rhs[2][i] = mode_fixed ? 0.0 : wd[i];
        rhs[3][i] = mode_fixed ? 0.0 : wdd;
    }

    /* Boundary i=N-1: outgoing Sommerfeld */
    int nl = n - 1;
    rhs[1][nl] = -(ud[nl] - ud[nl > 0 ? nl-1 : 0]) / dr
                 - sponge[nl] * ud[nl];
    if (!mode_fixed) {
        rhs[3][nl] = -(wd[nl] - wd[nl > 0 ? nl-1 : 0]) / dr
                     - sponge[nl] * wd[nl];
    }
}

/* ==================================================================== */
/*  RK4 step                                                             */
/* ==================================================================== */
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

/* ==================================================================== */
/*  Initialization                                                       */
/* ==================================================================== */
static void init_pulse(void)
{
    /* u = r*psi = amp * exp(-(r-r0)^2/(2w^2)) * sin(k*r) */
    /* Ingoing: u_dot = -u' (speed c=1 toward center) */
    for (int i = 0; i < Ngrid; i++) {
        double ri = rg[i];
        double env = amp_p * exp(-0.5*(ri - r0_p)*(ri - r0_p) / (wid_p*wid_p));
        double envp = env * (-(ri - r0_p) / (wid_p * wid_p));

        Y[0][i] = env * sin(k_p * ri);
        Y[1][i] = -(envp * sin(k_p * ri) + env * k_p * cos(k_p * ri));
        Y[2][i] = 0.0;
        Y[3][i] = 0.0;
    }
}

static void init_fixed_sigma(void)
{
    /* sigma = -|sig0| * exp(-r^2/(2*wsig^2))  (negative well) */
    for (int i = 0; i < Ngrid; i++) {
        double ri = rg[i];
        double sig = -fabs(sig0_f) * exp(-0.5 * ri * ri / (wid_f * wid_f));
        Y[2][i] = ri * sig;
        Y[3][i] = 0.0;
    }
}

static void init_standing_wave(void)
{
    /* Standing wave u = amp * r^(l+1) * exp(-r^2/(2*w^2)), u_dot = 0
     * Centered at origin, matches angular momentum barrier behavior.
     * For use with fixed sigma well + massive wave to test bound states. */
    for (int i = 0; i < Ngrid; i++) {
        double ri = rg[i];
        /* r^l factor (u=r*psi, psi ~ r^l, so u ~ r^{l+1}) */
        double rl = pow(ri, ell);
        Y[0][i] = amp_p * rl * exp(-0.5 * ri * ri / (wid_p * wid_p));
        Y[1][i] = 0.0;  /* standing: u_dot = 0 */
    }
}

/* ==================================================================== */
/*  Print effective potential V_eff(r) = m^2*A(r) + l(l+1)/r^2 + ...   */
/* ==================================================================== */
static void print_veff(void)
{
    init_fixed_sigma();
    compute_metric(Y);

    printf("# V_eff profile for massive wave in disformal well\n");
    printf("# alpha_c=%.4f  beta_d=%.4f  sig0=%.4f  wsig=%.4f\n",
           alpha_c, beta_d, sig0_f, wid_f);
    printf("# m_psi=%.4f  l=%d\n", m_psi, ell);
    printf("# %10s  %12s  %12s  %12s  %12s  %12s  %12s\n",
           "r", "sigma", "A", "v2_eff", "V_mass", "V_cent", "V_total");

    for (int i = 0; i < Ngrid && rg[i] < 30.0; i++) {
        double ri = rg[i];
        double sig = Y[2][i] / ri;
        double A = exp(2.0 * alpha_c * sig);
        if (A < A_MIN) A = A_MIN;

        double V_mass = m_psi * m_psi * A;
        double V_cent = (double)(ell * (ell + 1)) / (ri * ri);
        /* Also include the 1/r² correction from metric (C_r'/C_t)(u'-u/r) */
        double V_total = V_mass + V_cent;

        printf("%10.4f  %12.6f  %12.6f  %12.6f  %12.6f  %12.6f  %12.6f\n",
               ri, sig, A, MET_v2[i], V_mass, V_cent, V_total);
    }

    /* Find well depth */
    double V_inf = m_psi * m_psi;  /* A -> 1 at infinity */
    double V_min = 1e30;
    double r_min = 0;
    for (int i = 2; i < Ngrid; i++) {
        double ri = rg[i];
        double sig = Y[2][i] / ri;
        double A = exp(2.0 * alpha_c * sig);
        if (A < A_MIN) A = A_MIN;
        double V = m_psi * m_psi * A + (double)(ell*(ell+1)) / (ri*ri);
        if (V < V_min) { V_min = V; r_min = ri; }
    }
    printf("#\n# V_inf = %.6f (= m^2)\n", V_inf);
    printf("# V_min = %.6f at r = %.4f\n", V_min, r_min);
    printf("# Well depth = %.6f\n", V_inf - V_min);
    printf("# Bound states possible if well depth > 0 and well is wide enough\n");
}

/* ==================================================================== */
/*  Energy diagnostics                                                   */
/* ==================================================================== */
static double energy_psi(double rmin, double rmax_e)
{
    /* Full energy including mass and centrifugal terms:
     * H = (1/2) integral [ud^2 + v^2*u'^2 + (l(l+1)/r^2 + m^2*A)*u^2] dr
     * This is the conserved Hamiltonian for the wave equation. */
    double *u = Y[0], *ud = Y[1];
    double E = 0.0;

    for (int i = 0; i < Ngrid; i++) {
        if (rg[i] < rmin || rg[i] > rmax_e) continue;

        double up_i;
        if (i == 0)
            up_i = (u[1] + u[0]) / (2.0 * dr);
        else if (i == Ngrid - 1)
            up_i = (u[i] - u[i-1]) / dr;
        else
            up_i = (u[i+1] - u[i-1]) / (2.0 * dr);

        double ri = rg[i];
        double sig_i = Y[2][i] / ri;
        double A_i = exp(2.0 * alpha_c * sig_i);
        if (A_i < A_MIN) A_i = A_MIN;

        double KE = ud[i] * ud[i];
        double GE = MET_v2[i] * up_i * up_i;  /* gradient energy (radial) */
        double CE = (double)(ell*(ell+1)) * u[i]*u[i] / (ri*ri); /* centrifugal */
        double ME = m_psi * m_psi * A_i * u[i] * u[i];           /* mass */

        E += (KE + GE + CE + ME) * dr;
    }
    return 0.5 * E;
}

static double energy_sigma(void)
{
    /* E_sigma in flat space = (1/2)(sigma_t^2 + sigma_r^2) per 4*pi*r^2 dr
     * = (1/2)(wd^2 + qw^2) dr   (same derivation as psi, with C_t=C_r=1) */
    double *w = Y[2], *wd = Y[3];
    double E = 0.0;

    for (int i = 0; i < Ngrid; i++) {
        double wp_i;
        if (i == 0)
            wp_i = (w[1] + w[0]) / (2.0 * dr);
        else if (i == Ngrid - 1)
            wp_i = (w[i] - w[i-1]) / dr;
        else
            wp_i = (w[i+1] - w[i-1]) / (2.0 * dr);

        double qw = wp_i - w[i] / rg[i];
        E += (wd[i]*wd[i] + qw*qw) * dr;
    }
    return 0.5 * E;
}

/* ==================================================================== */
/*  Snapshot output                                                      */
/* ==================================================================== */
static void output_snapshot(int idx, double t)
{
    if (!do_snap) return;
    char fn[256];
    snprintf(fn, sizeof(fn), "snap_%04d.dat", idx);
    FILE *fp = fopen(fn, "w");
    if (!fp) return;
    fprintf(fp, "# t = %.6f\n", t);
    fprintf(fp, "# r  psi  psi_dot  sigma  n_refr  v_eff\n");
    for (int i = 0; i < Ngrid; i += 2) {
        double psi_i  = Y[0][i] / rg[i];
        double psid_i = Y[1][i] / rg[i];
        double sig_i  = Y[2][i] / rg[i];
        double v2 = MET_v2[i] > 1e-12 ? MET_v2[i] : 1e-12;
        fprintf(fp, "%.5f  %.8e  %.8e  %.8e  %.6f  %.6f\n",
                rg[i], psi_i, psid_i, sig_i, 1.0/sqrt(v2), sqrt(v2));
    }
    fclose(fp);
}

/* ==================================================================== */
/*  Main evolution loop                                                  */
/* ==================================================================== */
static void run_evolution(void)
{
    double dt = cfl * dr;
    int nsteps = (int)(Tmax / dt + 0.5);
    int nout = nsteps / 100;
    if (nout < 1) nout = 1;
    int nsnap = nsteps / 20;
    if (nsnap < 1) nsnap = 1;

    printf("# V7 Disformal Self-Trapping\n");
    printf("# alpha_c=%.4f  beta_d=%.4f  alpha_s=%.4f\n",
           alpha_c, beta_d, alpha_s);
    printf("# amp=%.2f  k=%.2f  r0=%.2f  w=%.2f  l=%d\n",
           amp_p, k_p, r0_p, wid_p, ell);
    printf("# N=%d  Rmax=%.1f  Tmax=%.1f  dt=%.6f  dr=%.6f\n",
           Ngrid, Rmax, Tmax, dt, dr);
    printf("# mu_sig=%.4f  m_psi=%.4f  source_sign=%d\n",
           mu_sig, m_psi, source_sign);
    printf("# mode=%s  init=%s\n",
           mode_fixed ? "fixed_sigma" : "coupled",
           init_standing ? "standing" : "pulse");
    printf("#\n");

    compute_metric(Y);
    double E0 = energy_psi(0.0, Rmax);
    double r_core = 10.0;

    printf("# Initial E_psi = %.6f\n#\n", E0);
    printf("# %8s  %10s  %10s  %10s  %10s  %10s  %8s\n",
           "t", "E_psi", "E_core", "E_sigma", "sig_min", "psi_max", "n_max");

    int snap_idx = 0;
    double E_core_late = 0.0;

    for (int step = 0; step <= nsteps; step++) {
        double t = step * dt;

        if (step % nout == 0 || step == nsteps) {
            compute_metric(Y);
            double Ep = energy_psi(0.0, Rmax);
            double Ec = energy_psi(0.0, r_core);
            double Es = mode_fixed ? 0.0 : energy_sigma();

            double sig_min = 0, psi_max = 0, n_max = 1;
            for (int i = 0; i < Ngrid; i++) {
                double si = Y[2][i] / rg[i];
                double pi = fabs(Y[0][i] / rg[i]);
                double v2 = MET_v2[i] > 1e-12 ? MET_v2[i] : 1e-12;
                double ni = 1.0 / sqrt(v2);
                if (si < sig_min) sig_min = si;
                if (pi > psi_max) psi_max = pi;
                if (ni > n_max) n_max = ni;
            }

            printf("%9.3f  %10.4f  %10.4f  %10.4f  %10.6f  %10.6f  %8.3f\n",
                   t, Ep, Ec, Es, sig_min, psi_max, n_max);
            fflush(stdout);

            if (t > 0.8 * Tmax) E_core_late = Ec;
        }

        if (do_snap && step % nsnap == 0) {
            compute_metric(Y);
            output_snapshot(snap_idx++, t);
        }

        /* Check for NaN */
        if (step % (nout * 10) == 0 && step > 0) {
            if (isnan(Y[0][Ngrid/2]) || isnan(Y[2][Ngrid/2])) {
                fprintf(stderr, "NaN detected at step %d, t=%.4f\n",
                        step, t);
                break;
            }
        }

        if (step < nsteps)
            rk4_step(dt);
    }

    double frac = (E0 > 1e-15) ? E_core_late / E0 : 0.0;
    printf("#\n");
    printf("# SUMMARY: E_core(late)/E0 = %.6f\n", frac);
    if (frac > 0.10)
        printf("# RESULT: TRAPPED (%.1f%%)\n", 100.0*frac);
    else
        printf("# RESULT: DISPERSED (%.1f%%)\n", 100.0*frac);
}

/* ==================================================================== */
/*  Scan mode: sweep alpha_s to find trapping threshold                  */
/* ==================================================================== */
static void run_scan(void)
{
    printf("# V7 Scan: alpha_s vs trapping\n");
    printf("# alpha_c=%.4f  beta_d=%.4f  amp=%.2f  k=%.2f\n",
           alpha_c, beta_d, amp_p, k_p);
    printf("# %10s  %10s  %10s  %8s  %s\n",
           "alpha_s", "E_frac", "sig_min", "n_max", "result");

    double as_vals[] = {0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0};
    int n_as = sizeof(as_vals) / sizeof(as_vals[0]);

    for (int ia = 0; ia < n_as; ia++) {
        alpha_s = as_vals[ia];

        memset(Y, 0, sizeof(Y));
        init_pulse();
        if (mode_fixed) init_fixed_sigma();

        compute_metric(Y);
        double E0 = energy_psi(0.0, Rmax);

        double dt = cfl * dr;
        int nsteps = (int)(Tmax / dt + 0.5);

        int crashed = 0;
        for (int step = 0; step < nsteps; step++) {
            rk4_step(dt);
            /* Check stability every 1000 steps */
            if (step % 1000 == 999) {
                if (isnan(Y[0][Ngrid/2]) || isnan(Y[2][Ngrid/2])) {
                    crashed = 1;
                    break;
                }
            }
        }

        double frac = 0, sig_min = 0, n_max = 1;
        if (!crashed) {
            compute_metric(Y);
            double Ec = energy_psi(0.0, 10.0);
            frac = (E0 > 1e-15) ? Ec / E0 : 0.0;

            for (int i = 0; i < Ngrid; i++) {
                double si = Y[2][i] / rg[i];
                double v2 = MET_v2[i] > 1e-12 ? MET_v2[i] : 1e-12;
                double ni = 1.0 / sqrt(v2);
                if (si < sig_min) sig_min = si;
                if (ni > n_max) n_max = ni;
            }
        }

        printf("%10.4f  %10.6f  %10.6f  %8.3f  %s\n",
               alpha_s, frac, sig_min, n_max,
               crashed ? "CRASHED" : (frac > 0.1 ? "TRAPPED" : "DISPERSED"));
        fflush(stdout);
    }
}

/* ==================================================================== */
/*  Hartree self-consistent solver                                       */
/*  Iterates: sigma -> V_eff -> psi eigenstate -> source -> sigma        */
/*  This finds static self-trapped solutions without time-evolution.     */
/* ==================================================================== */

/* Thomas algorithm: solve tridiagonal system a_i x_{i-1} + b_i x_i + c_i x_{i+1} = d_i */
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

/* Solve sigma BVP: -sigma'' - 2sigma'/r + mu^2*sigma = S(r)
 * In w=r*sigma form: -w'' + mu^2*w = r*S(r) + (w-related corrections)
 * Actually: sigma'' + 2sigma'/r - mu^2*sigma = -S(r)
 * Which is: (1/r^2) d/dr(r^2 sigma') - mu^2*sigma = -S(r)
 * In u=r*sigma form: u'' - mu^2*u = -r*S(r), u(0)=0, u(N)=0
 */
static void solve_sigma_bvp(double *source_r, double *w_out)
{
    int n = Ngrid;
    double A[NMAX], B[NMAX], C[NMAX], D[NMAX];

    double dr2 = dr * dr;
    for (int i = 0; i < n; i++) {
        A[i] = (i > 0) ? 1.0 / dr2 : 0.0;
        C[i] = (i < n-1) ? 1.0 / dr2 : 0.0;
        B[i] = -2.0 / dr2 - mu_sig * mu_sig;
        /* w'' - mu^2*w = -r*S(r) → w''=-r*S + mu^2*w → Ax=d form below */
        D[i] = -rg[i] * source_r[i];
    }
    /* BC: w(0) ~ 0 (extrapolate: w[-1]=-w[0]) */
    B[0] = -3.0 / dr2 - mu_sig * mu_sig;  /* ghost: w[-1]=-w[0] → w'' ≈ (w1-3w0)/dr2 */

    if (mu_sig > 1e-10) {
        /* Massive: w → 0 at infinity, Dirichlet BC */
        D[n-1] = 0.0;
        C[n-1] = 0.0;
        B[n-1] = 1.0;
        A[n-1] = 0.0;
    } else {
        /* Massless: w → D (const) at infinity, Neumann BC: w'(Rmax)=0
         * w[n-1] - w[n-2] = 0 → w[n-1] = w[n-2] */
        A[n-1] = -1.0;
        B[n-1] = 1.0;
        C[n-1] = 0.0;
        D[n-1] = 0.0;
    }

    thomas_solve(n, A, B, C, D, w_out);
}

/* Implicit imaginary time evolution to find ground state eigenfunction.
 * Solves: (1 + dtau*H)*u^{n+1} = u^n  (backward Euler)
 * where H = -d²/dr² + V_eff(r), tridiagonal.
 * This is unconditionally stable, handles singular centrifugal barrier.
 * Returns eigenvalue omega^2 = <u|H|u>/<u|u>.
 */
static double find_eigenstate(double *V_eff, double *u_out)
{
    int n = Ngrid;
    double u[NMAX], u_new[NMAX];
    double At[NMAX], Bt[NMAX], Ct[NMAX], Dt[NMAX];

    /* Initial guess: u ~ r^(l+1) * exp(-r^2/w^2) */
    for (int i = 0; i < n; i++) {
        double ri = rg[i];
        u[i] = pow(ri, ell) * exp(-ri*ri / 25.0);
    }
    u[n-1] = 0.0;

    /* Normalize initial guess */
    double norm2 = 0;
    for (int i = 0; i < n; i++) norm2 += u[i]*u[i]*dr;
    if (norm2 > 1e-30) {
        double inv = 1.0 / sqrt(norm2);
        for (int i = 0; i < n; i++) u[i] *= inv;
    }

    double dtau = 0.5;  /* Large step OK with implicit scheme */
    int niter = 5000;
    double dr2 = dr * dr;

    for (int iter = 0; iter < niter; iter++) {
        /* Set up tridiagonal system: (1 + dtau*H)*u_new = u
         * H = -d²/dr² + V_eff
         * -d²/dr² → tridiagonal with -1/dr², 2/dr², -1/dr² */

        /* Ghost point: u[-1] = -u[0] → u'' ≈ (u1 - 3u0)/dr² */
        At[0] = 0.0;
        Bt[0] = 1.0 + dtau * (3.0/dr2 + V_eff[0]);
        Ct[0] = -dtau / dr2;
        Dt[0] = u[0];

        for (int i = 1; i < n-1; i++) {
            At[i] = -dtau / dr2;
            Bt[i] = 1.0 + dtau * (2.0/dr2 + V_eff[i]);
            Ct[i] = -dtau / dr2;
            Dt[i] = u[i];
        }

        /* Boundary: u(N-1) = 0 */
        At[n-1] = 0.0;
        Bt[n-1] = 1.0;
        Ct[n-1] = 0.0;
        Dt[n-1] = 0.0;

        thomas_solve(n, At, Bt, Ct, Dt, u_new);

        /* Normalize */
        norm2 = 0.0;
        for (int i = 0; i < n; i++)
            norm2 += u_new[i]*u_new[i] * dr;
        if (norm2 < 1e-30) break;
        double inv = 1.0 / sqrt(norm2);
        for (int i = 0; i < n; i++)
            u[i] = u_new[i] * inv;
    }

    /* Compute eigenvalue: omega^2 = <u|H|u>/<u|u> */
    double hu = 0.0, uu = 0.0;
    /* i=0: ghost u[-1]=-u[0] → u'' = (u1-3u0)/dr² */
    {
        double udd = (u[1] - 3.0*u[0]) / dr2;
        hu += u[0] * (-udd + V_eff[0]*u[0]) * dr;
        uu += u[0]*u[0] * dr;
    }
    for (int i = 1; i < n-1; i++) {
        double udd = (u[i+1] - 2.0*u[i] + u[i-1]) / dr2;
        hu += u[i] * (-udd + V_eff[i]*u[i]) * dr;
        uu += u[i]*u[i] * dr;
    }
    double omega2 = hu / uu;

    for (int i = 0; i < n; i++)
        u_out[i] = u[i];

    return omega2;
}

static void run_hartree(void)
{
    int n = Ngrid;
    double w[NMAX], psi_u[NMAX], source[NMAX], V_eff[NMAX];
    double w_old[NMAX];

    printf("# V7 Hartree Self-Consistent Solver\n");
    printf("# alpha_c=%.4f  beta_d=%.4f  alpha_s=%.4f\n",
           alpha_c, beta_d, alpha_s);
    printf("# m_psi=%.4f  mu_sig=%.4f  l=%d  N=%d  Rmax=%.1f\n",
           m_psi, mu_sig, ell, Ngrid, Rmax);
    printf("#\n");
    printf("# %4s  %12s  %12s  %12s  %12s  %12s  %12s  %s\n",
           "iter", "omega2", "sig_min", "A_min", "V_min", "psi_max", "delta_w", "status");

    /* Bootstrap: start with a TRIAL psi (Gaussian), compute sigma from it,
     * then iterate. This avoids the chicken-and-egg problem of starting
     * with sigma=0 (which has V >= m^2, no bound states). */

    /* Trial psi_u = r^(l+1) * exp(-r^2/w^2), normalized */
    double trial_w = 3.0;  /* width of trial wavefunction */
    double norm2 = 0;
    for (int i = 0; i < n; i++) {
        double ri = rg[i];
        psi_u[i] = pow(ri, ell) * exp(-ri*ri / (trial_w*trial_w));
        norm2 += psi_u[i]*psi_u[i] * dr;
    }
    double inv_norm = 1.0 / sqrt(norm2);
    for (int i = 0; i < n; i++)
        psi_u[i] *= inv_norm;

    /* Compute initial sigma from trial psi */
    for (int i = 0; i < n; i++) {
        double psi2 = psi_u[i]*psi_u[i] / (rg[i]*rg[i]);
        source[i] = -alpha_s * m_psi*m_psi * psi2;
    }
    solve_sigma_bvp(source, w);

    double relax = 0.3;  /* under-relaxation factor */
    int max_iter = 500;
    double omega2_final = 0, sig_min_final = 0, A_min_final = 1;
    double V_min_final = 1, psi_max_final = 0;

    for (int iter = 0; iter < max_iter; iter++) {
        /* Save old w for convergence check */
        memcpy(w_old, w, n * sizeof(double));

        /* Compute sigma = w/r and A(sigma) */
        double sig_min = 0, A_min = 1;
        for (int i = 0; i < n; i++) {
            double sig = w[i] / rg[i];
            double A = exp(2.0 * alpha_c * sig);
            if (A < A_MIN) A = A_MIN;
            if (sig < sig_min) sig_min = sig;
            if (A < A_min) A_min = A;

            /* Effective potential: V_eff(r) = m^2*A + l(l+1)/r^2 */
            V_eff[i] = m_psi*m_psi * A + (double)(ell*(ell+1)) / (rg[i]*rg[i]);
        }

        /* Find eigenstate in this potential */
        double omega2 = find_eigenstate(V_eff, psi_u);

        /* Check if bound (omega^2 < m^2) */
        int is_bound = (omega2 < m_psi*m_psi);

        /* Compute source S(r) from psi^2 */
        /* Source = -alpha_s * m^2 * psi^2 → sigma < 0 → A < 1 → well */
        double psi_max = 0;
        for (int i = 0; i < n; i++) {
            double psi2 = psi_u[i]*psi_u[i] / (rg[i]*rg[i]);
            source[i] = -alpha_s * m_psi*m_psi * psi2;
            if (fabs(psi_u[i]/rg[i]) > psi_max) psi_max = fabs(psi_u[i]/rg[i]);
        }

        /* Solve sigma BVP with this source */
        double w_new[NMAX];
        solve_sigma_bvp(source, w_new);

        /* Under-relaxation mixing */
        double delta = 0;
        for (int i = 0; i < n; i++) {
            double wn = relax * w_new[i] + (1.0 - relax) * w_old[i];
            delta += (wn - w_old[i]) * (wn - w_old[i]) * dr;
            w[i] = wn;
        }
        delta = sqrt(delta);

        double V_min = 1e30;
        for (int i = 0; i < n; i++)
            if (V_eff[i] < V_min) V_min = V_eff[i];

        printf("%4d  %12.6f  %12.6f  %12.6f  %12.6f  %12.6f  %12.4e  %s\n",
               iter, omega2, sig_min, A_min, V_min, psi_max, delta,
               is_bound ? "BOUND" : "unbound");
        fflush(stdout);

        if (!is_bound && iter > 20) {
            printf("# Lost bound state at iter %d — well too shallow\n", iter);
            break;
        }

        omega2_final = omega2;
        sig_min_final = sig_min;
        A_min_final = A_min;
        V_min_final = V_min;
        psi_max_final = psi_max;

        if (delta < 1e-10 && iter > 5) {
            printf("# CONVERGED after %d iterations\n", iter+1);
            break;
        }
    }

    /* Extract dilaton charge D = lim_{r->inf} r*sigma = w(large r) */
    double D_charge = w[n-2];  /* w = r*sigma → D at large r for massless */
    double Q_source = 0;       /* total source integral */
    for (int i = 0; i < n; i++) {
        double psi2 = psi_u[i]*psi_u[i] / (rg[i]*rg[i]);
        Q_source += alpha_s * m_psi*m_psi * psi2 * 4.0*M_PI*rg[i]*rg[i] * dr;
    }

    printf("#\n# === CONVERGED SELF-TRAPPED SOLUTION ===\n");
    printf("# omega^2     = %.6f  (binding energy = %.6f)\n",
           omega2_final, m_psi*m_psi - omega2_final);
    printf("# sigma_min   = %.6f  (at core)\n", sig_min_final);
    printf("# A_min       = %.6f  (conformal factor at core)\n", A_min_final);
    printf("# V_min       = %.6f  (effective potential minimum)\n", V_min_final);
    printf("# psi_max     = %.6f\n", psi_max_final);
    printf("# D_charge    = %.6f  (dilaton charge, sigma ~ D/r at large r)\n", D_charge);
    printf("# Q_source    = %.6f  (total source, D = -Q/(4pi) = %.6f)\n",
           Q_source, -Q_source/(4*M_PI));
    printf("# D/M ratio   = %.6f  (M ~ omega, D/omega)\n",
           D_charge / sqrt(omega2_final));

    /* Print final profile */
    printf("#\n# Final self-consistent profile:\n");
    printf("# %10s  %12s  %12s  %12s  %12s  %12s\n",
           "r", "sigma", "A", "psi", "V_eff", "r*sigma");
    for (int i = 0; i < n; i += (n > 200 ? n/200 : 1)) {
        double sig = w[i] / rg[i];
        double A = exp(2.0 * alpha_c * sig);
        if (A < A_MIN) A = A_MIN;
        double psi = psi_u[i] / rg[i];
        double Vef = m_psi*m_psi * A + (double)(ell*(ell+1))/(rg[i]*rg[i]);
        printf("%10.4f  %12.6f  %12.6f  %12.6f  %12.6f  %12.6f\n",
               rg[i], sig, A, psi, Vef, w[i]);
    }
}

/* ==================================================================== */
/*  Main                                                                 */
/* ==================================================================== */
int main(int argc, char *argv[])
{
    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "-fixed"))    { mode_fixed = 1; continue; }
        if (!strcmp(argv[i], "-coupled"))  { mode_fixed = 0; continue; }
        if (!strcmp(argv[i], "-scan"))     { mode_scan = 1; continue; }
        if (!strcmp(argv[i], "-snap"))     { do_snap = 1; continue; }
        if (!strcmp(argv[i], "-veff"))     { mode_veff = 1; continue; }
        if (!strcmp(argv[i], "-eigen"))    { mode_eigen = 1; continue; }
        if (!strcmp(argv[i], "-hartree"))  { mode_hartree = 1; continue; }
        if (!strcmp(argv[i], "-standing")) { init_standing = 1; continue; }
        if (!strcmp(argv[i], "-repulsive")){ source_sign = +1; continue; }
        if (i + 1 < argc) {
            if (!strcmp(argv[i], "-alpha"))  { alpha_c = atof(argv[++i]); continue; }
            if (!strcmp(argv[i], "-beta"))   { beta_d  = atof(argv[++i]); continue; }
            if (!strcmp(argv[i], "-alphas")) { alpha_s = atof(argv[++i]); continue; }
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
            if (!strcmp(argv[i], "-mu"))     { mu_sig  = atof(argv[++i]); continue; }
            if (!strcmp(argv[i], "-mpsi"))   { m_psi   = atof(argv[++i]); continue; }
        }
    }

    setup_grid();
    memset(Y, 0, sizeof(Y));

    if (mode_veff) {
        print_veff();
        return 0;
    }

    if (mode_hartree) {
        run_hartree();
        return 0;
    }

    if (init_standing) {
        init_standing_wave();
    } else {
        init_pulse();
    }
    if (mode_fixed) init_fixed_sigma();

    if (mode_scan)
        run_scan();
    else
        run_evolution();

    return 0;
}
