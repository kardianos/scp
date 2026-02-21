/*
 * ebid.c — Einstein-Born-Infeld-Dilaton soliton solver
 *
 * Solves the coupled (m, delta, phi) system for the EBId soliton.
 * Uses the exact Tamaki-Torii (PRD 62, 061501, 2000) field equations.
 *
 * Metric: ds^2 = -f*e^{-2*delta} dt^2 + dr^2/f + r^2 dOmega^2
 *   f = 1 - 2m/r
 * Note: e^{-2 delta} (not e^{+2 delta}).
 *
 * ODEs (Tamaki-Torii Eqs. 8-10, purely electric):
 *   m'     = -U + (r^2/2)*f*(phi')^2
 *   delta' = -r*(phi')^2
 *   phi''  = -(2/r)*phi' - (2/f)*[(m/r + U)*phi'/r - X]
 *
 * where (electric case, Eqs. 11-12):
 *   H = sqrt(1 + Q^2/(b*r^4))             [phi-independent!]
 *   U = e^{2*gamma*phi} * b * r^2 * (1-H)  [always <= 0]
 *   X = b * e^{2*gamma*phi} * (H-1)         [always >= 0]
 *
 * Note: -U = r^2 * X  (convenient identity).
 *
 * Near-origin (TT Eq. 19, electric soliton):
 *   phi(r) ~ -(1/(2*gamma)) * ln(4*gamma^2*Q*sqrt(b)*|ln r|)
 *   m(r)   ~ r / (4*|ln r|)   [for standard coupling]
 *   phi diverges to -infinity logarithmically at r=0.
 *
 * Key result (Clement & Gal'tsov, PRD 62, 124013, 2000):
 *   The soliton family is CONTINUOUS, parameterized by a constant c.
 *   Any c > c_cr gives an asymptotically flat solution with phi->0.
 *   No fine-tuning (bisection) is needed — just integrate outward.
 *
 * Integration uses logarithmic radial variable s = ln(r).
 *
 * Refs: Tamaki & Torii, PRD 62 (2000) 061501 (gr-qc/0004071)
 *       Clement & Gal'tsov, PRD 62 (2000) 124013 (hep-th/0007228)
 */

#define _DEFAULT_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define NMAX 500001

/* Solution arrays */
static double r_arr[NMAX];
static double m_arr[NMAX];
static double del_arr[NMAX];
static double phi_arr[NMAX];
static double psi_arr[NMAX];   /* dφ/dr */

/* Parameters */
static double par_b     = 1.0;    /* BI field strength */
static double par_Q     = 1.0;    /* electric charge */
static double par_gamma = 1.0;    /* dilaton coupling (1 = string theory) */

/* Grid */
static double r_start;

/* Bailout */
#define BAIL_NONE     0
#define BAIL_HORIZON  1
#define BAIL_PHI_BIG  2
#define BAIL_PSI_BIG  3
#define BAIL_M_BIG    4
#define BAIL_NAN      5
static int bail_reason;

/* --------------------------------------------------------------- */
/*  TT field equations in log-radial variable s = ln(r)            */
/* --------------------------------------------------------------- */

static double compute_H(double r)
{
    double r4 = r * r * r * r;
    return sqrt(1.0 + par_Q * par_Q / (par_b * r4));
}

/*
 * RHS of the 4-ODE system in s = ln(r).
 * y[0]=m, y[1]=delta, y[2]=phi, y[3]=psi=dφ/dr
 * Computes dy/ds = r * dy/dr.
 *
 * Uses the EXACT Tamaki-Torii equations (Eqs. 8-10).
 */
static void rhs_log(double s, const double *y, double *dy)
{
    double r = exp(s);
    double m   = y[0];
    double phi = y[2];
    double psi = y[3];   /* dφ/dr */

    double f = 1.0 - 2.0 * m / r;
    double H = compute_H(r);
    double e2gp = exp(2.0 * par_gamma * phi);

    /* U and X (TT Eqs. 11-12, electric case) */
    double U = e2gp * par_b * r * r * (1.0 - H);   /* U <= 0 */
    double X = par_b * e2gp * (H - 1.0);             /* X >= 0 */

    /* TT Eq. 8: m' = -U + (r^2/2)*f*psi^2 */
    double mp = -U + (r * r / 2.0) * f * psi * psi;

    /* TT Eq. 9: delta' = -r*psi^2  (NOTE: negative!) */
    double dp = -r * psi * psi;

    /* TT Eq. 10: phi'' = -(2/r)*psi - (2/f)*[(m/r+U)*psi/r - X] */
    double phipp;
    if (fabs(f) < 1e-15) {
        phipp = 0.0;
    } else {
        phipp = -(2.0 / r) * psi
                - (2.0 / f) * ((m / r + U) * psi / r - X);
    }

    /* dy/ds = r * dy/dr */
    dy[0] = r * mp;
    dy[1] = r * dp;
    dy[2] = r * psi;
    dy[3] = r * phipp;
}

/* --------------------------------------------------------------- */
/*  RK4 integrator on log grid                                     */
/* --------------------------------------------------------------- */

static int integrate(double phi0, double psi0, double m0,
                     int Ngrid, double rmax)
{
    double s_start = log(r_start);
    double s_max   = log(rmax);
    double ds = (s_max - s_start) / Ngrid;
    bail_reason = BAIL_NONE;

    r_arr[0]   = r_start;
    m_arr[0]   = m0;
    del_arr[0] = 0.0;    /* will shift so delta(inf)=0 at end */
    phi_arr[0] = phi0;
    psi_arr[0] = psi0;

    for (int i = 0; i < Ngrid; i++) {
        double s = s_start + i * ds;
        double y[4] = { m_arr[i], del_arr[i], phi_arr[i], psi_arr[i] };
        double k1[4], k2[4], k3[4], k4[4], yt[4];

        rhs_log(s, y, k1);
        for (int j = 0; j < 4; j++) yt[j] = y[j] + 0.5 * ds * k1[j];
        rhs_log(s + 0.5 * ds, yt, k2);
        for (int j = 0; j < 4; j++) yt[j] = y[j] + 0.5 * ds * k2[j];
        rhs_log(s + 0.5 * ds, yt, k3);
        for (int j = 0; j < 4; j++) yt[j] = y[j] + ds * k3[j];
        rhs_log(s + ds, yt, k4);

        m_arr[i+1]   = y[0] + (ds/6.0)*(k1[0]+2*k2[0]+2*k3[0]+k4[0]);
        del_arr[i+1] = y[1] + (ds/6.0)*(k1[1]+2*k2[1]+2*k3[1]+k4[1]);
        phi_arr[i+1] = y[2] + (ds/6.0)*(k1[2]+2*k2[2]+2*k3[2]+k4[2]);
        psi_arr[i+1] = y[3] + (ds/6.0)*(k1[3]+2*k2[3]+2*k3[3]+k4[3]);

        double rnext = exp(s + ds);
        r_arr[i+1] = rnext;
        double f = 1.0 - 2.0 * m_arr[i+1] / rnext;

        if (f < 0.0) {
            bail_reason = BAIL_HORIZON;
            return i + 1;
        }
        if (!isfinite(phi_arr[i+1]) || !isfinite(m_arr[i+1]) ||
            !isfinite(psi_arr[i+1]) || !isfinite(del_arr[i+1])) {
            bail_reason = BAIL_NAN;
            return i + 1;
        }
        if (fabs(phi_arr[i+1]) > 100.0) {
            bail_reason = BAIL_PHI_BIG;
            return i + 1;
        }
        if (fabs(psi_arr[i+1]) > 1e15) {
            bail_reason = BAIL_PSI_BIG;
            return i + 1;
        }
        if (fabs(m_arr[i+1]) > 1e10) {
            bail_reason = BAIL_M_BIG;
            return i + 1;
        }
    }
    return Ngrid;
}

/* --------------------------------------------------------------- */
/*  Near-origin initialization (TT Eq. 19 + CG sub-leading)       */
/* --------------------------------------------------------------- */

/*
 * Leading-order Tamaki-Torii asymptotic at r_start:
 *   phi ~ -(1/(2g))ln(C*u)   where C = 4*g^2*Q*sqrt(b), u = |ln r|
 *   psi ~ 1/(2*g*r*u)
 *   m   ~ r/(4u)             (for g = 1/(8*pi*gamma^2) ≈ TT normalization)
 *
 * The constant c in CG expansion shifts the sub-leading terms.
 */
static int init_quiet = 0;   /* suppress verbose output during bisection */

static void init_asymptotic(double *phi0, double *psi0, double *m0,
                            double c_param)
{
    double g = par_gamma;
    double u = -log(r_start);   /* u = |ln r_start| > 0 */
    double C = 4.0 * g * g * par_Q * sqrt(par_b);

    /* Leading-order */
    *phi0 = -(1.0 / (2.0 * g)) * log(C * u);
    *psi0 = 1.0 / (2.0 * g * r_start * u);

    /* Sub-leading correction from CG expansion */
    double g_grav = 1.0 / (4.0 * M_PI * g * g); /* G/gamma^2 in TT units */
    double tau = log(r_start);  /* tau < 0 */
    double L = (2.0 * g_grav - 1.0) * log(-tau) + c_param;

    *phi0 += -(1.0 / (2.0 * g)) * L / u;
    *m0 = g_grav * r_start / u;

    if (!init_quiet) {
        fprintf(stderr, "# Init: r_start=%.2e u=%.4f g_grav=%.6f\n",
                r_start, u, g_grav);
        fprintf(stderr, "# Init: phi0=%.6f psi0=%.6e m0=%.6e (c=%.4f L=%.4f)\n",
                *phi0, *psi0, *m0, c_param, L);
    }
}

/* --------------------------------------------------------------- */
/*  Soliton finder: bisect on c for phi(rmax) -> 0                 */
/* --------------------------------------------------------------- */

/*
 * Integrate with family parameter c and return phi_inf (extrapolated).
 * Sets *bailed = 1 if integration didn't reach rmax.
 * Returns phi_inf when bail_reason == BAIL_NONE.
 */
static double try_c(double c, int Ngrid, double rmax, int *bailed)
{
    double phi0, psi0, m0;
    init_asymptotic(&phi0, &psi0, &m0, c);

    integrate(phi0, psi0, m0, Ngrid, rmax);

    if (bail_reason != BAIL_NONE) {
        *bailed = 1;
        return 0.0;
    }
    *bailed = 0;

    /* Extrapolate phi_inf from phi = phi_inf + D/r at two far-field points. */
    int i1 = (int)(0.93 * Ngrid);
    int i2 = (int)(0.99 * Ngrid);
    double r1 = r_arr[i1], r2 = r_arr[i2];
    double p1 = phi_arr[i1], p2 = phi_arr[i2];
    double phi_inf = (p2 * r2 - p1 * r1) / (r2 - r1);

    return phi_inf;
}

/*
 * Find the soliton by bisecting on c (CG family parameter).
 * Returns the optimal c, or NAN on failure.
 * The soliton has phi(rmax) = 0.
 */
static double find_soliton(int Ngrid, double rmax)
{
    fprintf(stderr, "# Finding soliton: bisecting on c for phi_inf=0\n");
    init_quiet = 1;

    /* Scan to find bracket: c_lo with phi_inf>0, c_hi with phi_inf<0 */
    double c_lo = -5.0, c_hi = -5.0;
    int found_pos = 0, found_neg = 0;

    /*
     * Two-phase bracket search:
     * Phase 1: scan c upward, skip bail solutions. Among non-bailing,
     *          find first c_pos (phi_inf>0) and first c_neg (phi_inf<0).
     * Phase 2: if no c_pos found but c_neg exists, the soliton is at the
     *          bail/non-bail transition. Bracket: [last_bail_c, c_neg].
     */
    double last_bail_c = -1.0;
    int had_bail = 0;

    for (double c = 0.0; c <= 100.0; c += 0.5) {
        int bailed;
        double phi_inf = try_c(c, Ngrid, rmax, &bailed);

        if (bailed) {
            last_bail_c = c;
            had_bail = 1;
            continue;
        }

        if (phi_inf > 0.0) {
            c_lo = c;
            found_pos = 1;
        }
        if (phi_inf < 0.0) {
            if (found_pos) {
                /* Normal bracket: phi_inf went from + to - */
                c_hi = c;
                found_neg = 1;
            } else if (had_bail) {
                /* Soliton at bail transition: bail → negative phi_inf */
                c_lo = last_bail_c;
                c_hi = c;
                found_pos = 1;
                found_neg = 1;
            }
            break;
        }
    }

    if (!found_pos || !found_neg) {
        fprintf(stderr, "# ERROR: no bracket found (pos=%d neg=%d bail=%d)\n",
                found_pos, found_neg, had_bail);
        return NAN;
    }

    fprintf(stderr, "# Bracket: c in [%.6f, %.6f]\n", c_lo, c_hi);

    /* Bisection: 60 iterations.
     * If solution bails, treat as "c too low" (soliton is above).
     * If solution reaches rmax with phi_inf > 0, c too low.
     * If phi_inf < 0, c too high. */
    for (int iter = 0; iter < 60; iter++) {
        double c_mid = 0.5 * (c_lo + c_hi);
        int bailed;
        double phi_inf = try_c(c_mid, Ngrid, rmax, &bailed);

        if (bailed || phi_inf > 0.0)
            c_lo = c_mid;
        else
            c_hi = c_mid;
    }

    double c_best = 0.5 * (c_lo + c_hi);

    /* Final integration (verbose) */
    init_quiet = 0;
    double phi0, psi0, m0;
    init_asymptotic(&phi0, &psi0, &m0, c_best);
    integrate(phi0, psi0, m0, Ngrid, rmax);

    fprintf(stderr, "# Soliton found: c=%.10f\n", c_best);
    return c_best;
}

/* --------------------------------------------------------------- */
/*  Observables                                                    */
/* --------------------------------------------------------------- */

static void compute_observables(int npts,
                                double *M_out, double *D_out,
                                double *delta_shift)
{
    /* ADM mass: m(rmax) */
    *M_out = m_arr[npts];

    /* Dilaton charge: extrapolate from two-point fit phi = phi_inf + D/r.
     * Use far-field points (93% and 99% of log grid). */
    int i1 = (int)(0.93 * npts);
    int i2 = (int)(0.99 * npts);
    double r1 = r_arr[i1], r2 = r_arr[i2];
    double p1 = phi_arr[i1], p2 = phi_arr[i2];

    /* From phi = phi_inf + D/r:
     *   p1 = phi_inf + D/r1,  p2 = phi_inf + D/r2
     *   D = (p1-p2)*r1*r2/(r2-r1) */
    if (fabs(r2 - r1) > 1e-15) {
        double phi_inf = (p2 * r2 - p1 * r1) / (r2 - r1);
        *D_out = (p1 - p2) * r1 * r2 / (r2 - r1);
        /* Report phi_inf for diagnostics */
        fprintf(stderr, "#   phi_inf   = %.8e (should be ~0 for soliton)\n",
                phi_inf);
    } else {
        *D_out = r2 * p2;
    }

    /* Delta shift */
    *delta_shift = del_arr[npts];
}

/* --------------------------------------------------------------- */
/*  Output                                                         */
/* --------------------------------------------------------------- */

static void print_profile(int npts, double phi0, double psi0, FILE *fp)
{
    int step = (npts > 5000) ? npts / 5000 : 1;

    double M, D, delta_shift;
    compute_observables(npts, &M, &D, &delta_shift);

    fprintf(fp, "# EBId soliton (TT convention): Q=%.6f b=%.6f gamma=%.6f\n",
            par_Q, par_b, par_gamma);
    fprintf(fp, "# phi0=%.10f  psi0=%.10e  M=%.8f  D=%.8f\n",
            phi0, psi0, M, D);
    fprintf(fp, "# M/|D|=%.6f  (BPS: M^2+D^2=%.6f, Q^2=%.6f)\n",
            (fabs(D) > 1e-15) ? M / fabs(D) : 0.0,
            M * M + D * D, par_Q * par_Q);
    fprintf(fp, "# delta_shift=%.8f  (physical delta = computed - %.8f)\n",
            delta_shift, delta_shift);
    fprintf(fp, "# r  m(r)  delta(r)  phi(r)  psi(r)  f(r)\n");

    for (int i = 0; i <= npts; i += step) {
        double r = r_arr[i];
        double f = 1.0 - 2.0 * m_arr[i] / r;

        fprintf(fp, "%.10e %.10e %.10e %.10e %.10e %.10e\n",
                r, m_arr[i], del_arr[i] - delta_shift, phi_arr[i],
                psi_arr[i], f);
    }
}

/* --------------------------------------------------------------- */
/*  Main                                                           */
/* --------------------------------------------------------------- */

static void usage(const char *prog)
{
    fprintf(stderr, "Usage: %s [options]\n", prog);
    fprintf(stderr, "  -Qe <charge>    Electric charge (default 1.0)\n");
    fprintf(stderr, "  -b <beta>       BI field strength (default 1.0)\n");
    fprintf(stderr, "  -gamma <g>      Dilaton coupling (default 1.0)\n");
    fprintf(stderr, "  -c <c_param>    CG family parameter (skip auto-find)\n");
    fprintf(stderr, "  -find           Auto-find soliton (default if no -c)\n");
    fprintf(stderr, "  -N <grid>       Grid points (default 300000)\n");
    fprintf(stderr, "  -rmax <rmax>    Outer boundary (default 200)\n");
    fprintf(stderr, "  -rstart <rs>    Inner boundary (default 1e-8)\n");
    fprintf(stderr, "  -scan           Scan over Q_e\n");
    fprintf(stderr, "  -cscan          Scan over c (family parameter)\n");
    fprintf(stderr, "  -outdir <dir>   Write profile to file\n");
}

int main(int argc, char **argv)
{
    int Ngrid = 300000;
    double rmax = 200.0;
    double r_start_user = -1.0;
    double c_param = 0.0;
    int mode_scan = 0;
    int mode_cscan = 0;
    int mode_find = 0;
    int c_set = 0;
    char *outdir = NULL;

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-Qe") == 0 && i + 1 < argc)
            par_Q = atof(argv[++i]);
        else if (strcmp(argv[i], "-b") == 0 && i + 1 < argc)
            par_b = atof(argv[++i]);
        else if (strcmp(argv[i], "-gamma") == 0 && i + 1 < argc)
            par_gamma = atof(argv[++i]);
        else if (strcmp(argv[i], "-c") == 0 && i + 1 < argc) {
            c_param = atof(argv[++i]);
            c_set = 1;
        }
        else if (strcmp(argv[i], "-N") == 0 && i + 1 < argc)
            Ngrid = atoi(argv[++i]);
        else if (strcmp(argv[i], "-rmax") == 0 && i + 1 < argc)
            rmax = atof(argv[++i]);
        else if (strcmp(argv[i], "-rstart") == 0 && i + 1 < argc)
            r_start_user = atof(argv[++i]);
        else if (strcmp(argv[i], "-scan") == 0)
            mode_scan = 1;
        else if (strcmp(argv[i], "-cscan") == 0)
            mode_cscan = 1;
        else if (strcmp(argv[i], "-find") == 0)
            mode_find = 1;
        else if (strcmp(argv[i], "-outdir") == 0 && i + 1 < argc)
            outdir = argv[++i];
        else {
            usage(argv[0]);
            return 1;
        }
    }

    if (Ngrid > NMAX - 1) Ngrid = NMAX - 1;

    r_start = (r_start_user > 0) ? r_start_user : 1e-8;

    /* ----------------------------------------------------------- */
    /*  Scan over c (family parameter)                             */
    /* ----------------------------------------------------------- */
    if (mode_cscan) {
        printf("# EBId c-scan: Q=%.6f b=%.6f gamma=%.6f rstart=%.2e\n",
               par_Q, par_b, par_gamma, r_start);
        printf("# c         phi0         psi0         M            D            "
               "M/|D|      M^2+D^2    delta_inf\n");

        double c_values[] = {-2, -1, 0, 1, 2, 3, 5, 8, 10, 15, 20};
        int nc = sizeof(c_values) / sizeof(c_values[0]);

        for (int ic = 0; ic < nc; ic++) {
            double phi0, psi0, m0;
            init_asymptotic(&phi0, &psi0, &m0, c_values[ic]);

            int nv = integrate(phi0, psi0, m0, Ngrid, rmax);
            int npts = (nv < Ngrid) ? nv : Ngrid;

            double M, D, ds;
            compute_observables(npts, &M, &D, &ds);

            printf("%6.2f  %12.6f  %12.4e  %12.6f  %12.6f  %10.6f  %10.6f  %10.6f  bail=%d\n",
                   c_values[ic], phi0, psi0, M, D,
                   (fabs(D) > 1e-15) ? M / fabs(D) : 0.0,
                   M * M + D * D, ds, bail_reason);
            fflush(stdout);
        }
        return 0;
    }

    /* ----------------------------------------------------------- */
    /*  Scan over Q_e                                              */
    /* ----------------------------------------------------------- */
    if (mode_scan) {
        printf("# EBId Q-scan (auto-find): b=%.6f gamma=%.6f rstart=%.2e\n",
               par_b, par_gamma, r_start);
        printf("# Q_e        c_best       M            D            M/|D|      "
               "M^2+D^2    Q^2        delta_inf\n");

        double Q_values[] = {0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0};
        int nQ = sizeof(Q_values) / sizeof(Q_values[0]);

        for (int iq = 0; iq < nQ; iq++) {
            par_Q = Q_values[iq];
            double c_best = find_soliton(Ngrid, rmax);

            double M, D, ds;
            compute_observables(Ngrid, &M, &D, &ds);

            printf("%.4f  %12.6f  %12.6f  %12.6f  %10.6f  %10.6f  %10.6f  %10.6f\n",
                   par_Q, c_best, M, D,
                   (fabs(D) > 1e-15) ? M / fabs(D) : 0.0,
                   M * M + D * D, par_Q * par_Q, ds);
            fflush(stdout);
        }
        return 0;
    }

    /* ----------------------------------------------------------- */
    /*  Default: single soliton                                    */
    /* ----------------------------------------------------------- */
    fprintf(stderr, "# EBId solver (TT convention): Q=%.4f b=%.4f gamma=%.4f\n",
            par_Q, par_b, par_gamma);
    fprintf(stderr, "# Ngrid=%d rmax=%.1f r_start=%.4e (log grid)\n",
            Ngrid, rmax, r_start);

    double phi0, psi0, m0;
    int npts;

    if (mode_find || !c_set) {
        /* Auto-find soliton by bisecting on c */
        double c_best = find_soliton(Ngrid, rmax);
        if (isnan(c_best)) {
            fprintf(stderr, "# FAILED to find soliton\n");
            return 1;
        }
        c_param = c_best;
        npts = Ngrid;  /* find_soliton leaves arrays filled */
        phi0 = phi_arr[0];
        psi0 = psi_arr[0];
        m0   = m_arr[0];
    } else {
        init_asymptotic(&phi0, &psi0, &m0, c_param);
        int nv = integrate(phi0, psi0, m0, Ngrid, rmax);
        npts = (nv < Ngrid) ? nv : Ngrid;
    }

    if (bail_reason != BAIL_NONE) {
        fprintf(stderr, "# WARNING: integration bailed at reason=%d\n",
                bail_reason);
        fprintf(stderr, "#   r=%.6e phi=%.6f m=%.6e f=%.6e\n",
                r_arr[npts], phi_arr[npts], m_arr[npts],
                1.0 - 2.0 * m_arr[npts] / r_arr[npts]);
    }

    double M, D, delta_shift;
    compute_observables(npts, &M, &D, &delta_shift);

    fprintf(stderr, "# Solution at rmax=%.1f (c=%.10f):\n", r_arr[npts], c_param);
    fprintf(stderr, "#   phi0      = %.10f\n", phi0);
    fprintf(stderr, "#   psi0      = %.10e\n", psi0);
    fprintf(stderr, "#   M (ADM)   = %.8f\n", M);
    fprintf(stderr, "#   D (dil)   = %.8f\n", D);
    fprintf(stderr, "#   M/|D|     = %.6f\n",
            (fabs(D) > 1e-15) ? M / fabs(D) : 0.0);
    fprintf(stderr, "#   M^2+D^2   = %.6f  (Q^2=%.6f)\n",
            M * M + D * D, par_Q * par_Q);
    fprintf(stderr, "#   delta_inf = %.8f\n", delta_shift);
    fprintf(stderr, "#   phi(rmax) = %.8e\n", phi_arr[npts]);
    fprintf(stderr, "#   psi(rmax) = %.8e\n", psi_arr[npts]);

    print_profile(npts, phi0, psi0, stdout);

    if (outdir) {
        char fname[512];
        snprintf(fname, sizeof(fname),
                 "%s/ebid_Q%.3f_b%.3f_c%.1f.dat",
                 outdir, par_Q, par_b, c_param);
        FILE *fp = fopen(fname, "w");
        if (fp) {
            print_profile(npts, phi0, psi0, fp);
            fclose(fp);
            fprintf(stderr, "# Profile written to %s\n", fname);
        }
    }

    return 0;
}
