/*
 * Scalar Q-ball solver — Phase 1 of v4 research program
 *
 * Model: L = |dot{phi}|^2 - |grad phi|^2 - V(|phi|^2)
 *   V(s) = m^2 s - mu^2 s^2 + lam s^3     (s = |phi|^2)
 *   Default: m^2=1, mu^2=2, lam=1  ->  V(s) = s(1-s)^2
 *
 * Q-ball ansatz: phi(x,t) = f(r) e^{i omega t}
 * Profile ODE: f'' + (2/r)f' + [omega^2 - V'(f^2)]f = 0
 *   where V'(s) = 1 - 4s + 3s^2
 *
 * Shooting method: At large r, f ~ A e^{-kappa r}/r + B e^{+kappa r}/r
 *   where kappa = sqrt(m^2 - omega^2).
 *   Q-ball has B = 0 (pure decaying). We bisect on f(0) to make B = 0.
 *
 * Key subtlety: V'(0) = m^2 > omega^2, so the linear force at f=0 is
 *   repulsive — the profile never crosses zero. Instead, the sign of B
 *   is detected via the log-derivative: f'/f + kappa at large r.
 *   - B > 0: log-derivative > -kappa (growing mode pulls f upward)
 *   - B < 0: log-derivative < -kappa (growing mode pulls f downward)
 *   - B = 0: log-derivative = -kappa (pure exponential decay)
 *
 * Usage:
 *   qball -outdir data              (scan over omega)
 *   qball -omega 0.5 -outdir data   (single omega, output profile)
 *   qball -diagnose 0.5             (diagnostic: scan f0 for given omega)
 */

#define _USE_MATH_DEFINES
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define NR_MAX 100001

static double f_arr[NR_MAX];
static double fp_arr[NR_MAX];

/* Model: V(s) = M2*s - MU2*s^2 + LAM*s^3 */
static const double M2 = 1.0, MU2 = 2.0, LAM = 1.0;
static const double MASS = 1.0; /* m = sqrt(M2) */

/* DBI parameter: b = BI field limit.
 * b_dbi = 0 means standard (no DBI).
 * Finite b_dbi gives P(X) = b^2(sqrt(1+2X/b^2)-1) kinetic term.
 * Constraint: omega*f0 < b/sqrt(2).
 * Effective speed at center: v/c = sqrt(1 - 2*omega^2*f0^2/b^2). */
static double b_dbi = 0.0;

/* Angular momentum quantum number: 0 = spherical, 1 = dipolar (ell=1 Y_1^0) */
static int ell = 0;

static double V_of_s(double s)
{
    return M2 * s - MU2 * s * s + LAM * s * s * s;
}

static double dV_ds(double s)
{
    return M2 - 2.0 * MU2 * s + 3.0 * LAM * s * s;
}

/*
 * Angular-averaged potential for ell=1 (Y_1^0 harmonic).
 *
 * With phi = f(r) Y_1^0, |phi|^2 = f^2 * (3/(4pi)) cos^2(theta).
 * Projecting V(|phi|^2) onto the ell=1 angular mode:
 *   V_eff(s) = s - 9s^2/(10*pi) + 27s^3/(112*pi^2)
 *   V'_eff(s) = 1 - 9s/(5*pi) + 81s^2/(112*pi^2)
 * where s = f^2.
 */
static double Veff_of_s(double s)
{
    if (ell == 0) return V_of_s(s);
    return s - 9.0 * s * s / (10.0 * M_PI)
             + 27.0 * s * s * s / (112.0 * M_PI * M_PI);
}

static double dVeff_ds(double s)
{
    if (ell == 0) return dV_ds(s);
    return 1.0 - 9.0 * s / (5.0 * M_PI)
               + 81.0 * s * s / (112.0 * M_PI * M_PI);
}

/*
 * Compute g' = f'' for the profile ODE.
 *
 * Standard:  g' = -2g/r + (V'(f^2) - omega^2) f
 *
 * DBI:       P(X) = b^2(sqrt(1+2X/b^2) - 1),  X = g^2 - omega^2*f^2
 *   P' = 1/Gamma,  P'' = -1/(b^2 Gamma^3),  Gamma = sqrt(1+2X/b^2)
 *   g' = (b^2 Gamma^3 / (b^2 - 2*omega^2*f^2)) *
 *        [-2g/(r*Gamma) + (V' - omega^2/Gamma + 2*omega^2*g^2/(b^2*Gamma^3)) * f]
 *
 * Note the sign: P'' = -1/(b^2 Gamma^3), so the 2P''omega^2g^2 term is
 *   2*(-1/(b^2*Gamma^3))*omega^2*g^2 = -2*omega^2*g^2/(b^2*Gamma^3)
 *   but in the numerator it appears with + (since EOM has P''*stuff + V'f = 0
 *   rearranged as g' = ...).
 */
static double compute_gp(double f, double g, double r, double omega)
{
    double w2 = omega * omega;
    double Vp = dVeff_ds(f * f);

    if (b_dbi <= 0.0) {
        /* Standard: g' = -2g/r + ell(ell+1)f/r^2 + (V' - omega^2)*f */
        double gp = (Vp - w2) * f;
        if (r > 1e-15) {
            gp += -2.0 * g / r;
            if (ell > 0) gp += ell * (ell + 1) * f / (r * r);
        }
        return gp;
    }

    /* DBI mode */
    double b2 = b_dbi * b_dbi;
    double X = g * g - w2 * f * f;
    double Gam2 = 1.0 + 2.0 * X / b2;
    if (Gam2 < 1e-10) Gam2 = 1e-10; /* safety clamp */
    double Gam = sqrt(Gam2);
    double Gam3 = Gam * Gam2;

    double denom = b2 - 2.0 * w2 * f * f;
    if (fabs(denom) < 1e-30) denom = 1e-30; /* safety */

    /* Numerator terms */
    double n_fric = 0.0;
    if (r > 1e-15) n_fric = -2.0 * g / (r * Gam);

    double n_field = (Vp - w2 / Gam + 2.0 * w2 * g * g / (b2 * Gam3)) * f;

    return (b2 * Gam3 / denom) * (n_fric + n_field);
}

/*
 * Integrate profile ODE from r=0 to r=rmax using RK4.
 * Stores f_arr[], fp_arr[].
 * Returns: 0 on success, 1 if |f| > 1e10 (diverged).
 */
static int shoot(double f0, double omega, int N, double rmax)
{
    double h = rmax / N;

    if (ell > 0) {
        /* ell>0: f(0)=0, f'(0)=a (f0 = slope parameter).
         * Taylor: f = a*r + c*r^3, where c = (m^2-omega^2)*a/(2*(2*ell+3)).
         * For ell=1: c = (1-omega^2)*a/10. */
        double a = f0;
        double c = (M2 - omega * omega) * a / (2.0 * (2 * ell + 3));
        f_arr[0] = 0.0;
        fp_arr[0] = a;
        f_arr[1] = a * h + c * h * h * h;
        fp_arr[1] = a + 3.0 * c * h * h;
    } else {
        /* ell=0: Taylor start f(r) ~ f0 + f2*r^2 */
        double Vp0 = dV_ds(f0 * f0);
        double f2;
        if (b_dbi > 0.0) {
            double Gam0 = sqrt(1.0 - 2.0 * omega * omega * f0 * f0
                               / (b_dbi * b_dbi));
            f2 = (Gam0 * Vp0 - omega * omega) * f0 / 6.0;
        } else {
            f2 = (Vp0 - omega * omega) * f0 / 6.0;
        }
        f_arr[0] = f0;
        fp_arr[0] = 0.0;
        f_arr[1] = f0 + f2 * h * h;
        fp_arr[1] = 2.0 * f2 * h;
    }

    for (int i = 1; i < N; i++) {
        double r = i * h;
        double f = f_arr[i], g = fp_arr[i];

        double k1f = g;
        double k1g = compute_gp(f, g, r, omega);

        double rm = r + 0.5 * h;
        double fm = f + 0.5 * h * k1f;
        double gm = g + 0.5 * h * k1g;
        double k2f = gm;
        double k2g = compute_gp(fm, gm, rm, omega);

        fm = f + 0.5 * h * k2f;
        gm = g + 0.5 * h * k2g;
        double k3f = gm;
        double k3g = compute_gp(fm, gm, rm, omega);

        double re = r + h;
        fm = f + h * k3f;
        gm = g + h * k3g;
        double k4f = gm;
        double k4g = compute_gp(fm, gm, re, omega);

        f_arr[i + 1] = f + (h / 6.0) * (k1f + 2 * k2f + 2 * k3f + k4f);
        fp_arr[i + 1] = g + (h / 6.0) * (k1g + 2 * k2g + 2 * k3g + k4g);

        if (fabs(f_arr[i + 1]) > 1e10) return 1; /* diverged */
    }
    return 0;
}

/*
 * Detect zero-crossing: check if f(r) goes negative anywhere.
 * Returns +1 if f stays positive (undershoot / f0 too small).
 * Returns -1 if f goes negative (overshoot / f0 too large).
 * Also returns f_min via pointer.
 *
 * For V(s) = s(1-s)^2: the Q-ball profile starts at f0 near the
 * false vacuum, decreases through the core, and must approach f=0
 * at infinity. The mechanical potential has a valley around f=sqrt(s_lo)
 * where V'(f^2) = omega^2. With too little energy the profile gets
 * stuck in this valley (undershoot). With too much it crosses f=0
 * and goes negative (overshoot). The Q-ball is the boundary.
 */
static int zero_crossing(int N, double *fmin_out)
{
    double fmin = f_arr[0];
    for (int i = 1; i <= N; i++) {
        if (f_arr[i] < fmin) fmin = f_arr[i];
        if (f_arr[i] < 0) {
            if (fmin_out) *fmin_out = fmin;
            return -1; /* overshoot */
        }
    }
    if (fmin_out) *fmin_out = fmin;
    return +1; /* undershoot */
}

/*
 * Find Q-ball f(0) for given omega by bisecting on zero-crossing.
 *
 * For f0 near sqrt(s_lo): profile stays near sqrt(s_lo) forever (undershoot).
 * For f0 near f_fv: profile rolls hard, overshoots past f=0 (overshoot).
 * The Q-ball is the critical f0 between these behaviors.
 *
 * Bracket: a = undershoot (f stays positive), b = overshoot (f goes negative).
 */
static double find_f0(double omega, int N, double rmax)
{
    if (ell > 0) {
        /* ell>0: shooting parameter is slope a = f'(0).
         * Scan a from small to large, find undershoot/overshoot transition. */
        double a_min = 0.01, a_max = 10.0;
        int nscan = 500;

        /* Find bracket: first overshoot */
        double prev_a = a_min;
        int prev_s = 0;
        shoot(a_min, omega, N, rmax);
        double fm;
        prev_s = zero_crossing(N, &fm);

        double a_lo = a_min, a_hi = -1;
        int s_lo = prev_s, s_hi = 0;

        for (int j = 1; j <= nscan; j++) {
            double a = a_min + j * (a_max - a_min) / nscan;
            int div = shoot(a, omega, N, rmax);
            if (div) {
                /* Diverged = overshoot */
                a_hi = a;
                s_hi = -1;
                a_lo = prev_a;
                s_lo = prev_s;
                break;
            }
            int sj = zero_crossing(N, &fm);
            if (sj != prev_s) {
                a_lo = prev_a;
                s_lo = prev_s;
                a_hi = a;
                s_hi = sj;
                break;
            }
            prev_a = a;
            prev_s = sj;
        }

        if (a_hi < 0) {
            fprintf(stderr, "  omega=%.4f ell=%d: no overshoot found\n",
                    omega, ell);
            return -1;
        }

        /* Ensure a_lo is undershoot, a_hi is overshoot */
        if (s_lo < 0 && s_hi > 0) {
            double tmp = a_lo; a_lo = a_hi; a_hi = tmp;
        }

        /* Bisect */
        for (int iter = 0; iter < 200; iter++) {
            double mid = 0.5 * (a_lo + a_hi);
            if (a_hi - a_lo < 1e-14 * mid) break;
            int div = shoot(mid, omega, N, rmax);
            if (div) { a_hi = mid; continue; }
            int sm = zero_crossing(N, &fm);
            if (sm > 0) a_lo = mid;
            else a_hi = mid;
        }

        double a_best = 0.5 * (a_lo + a_hi);
        shoot(a_best, omega, N, rmax);
        return a_best;
    }

    /* ell=0: False vacuum roots: V'(s) = omega^2 => 3s^2 - 4s + (1-omega^2) = 0 */
    double disc = 16.0 - 12.0 * (1.0 - omega * omega);
    if (disc < 0) {
        fprintf(stderr, "No false vacuum for omega=%.4f\n", omega);
        return -1;
    }
    double s_lo = (4.0 - sqrt(disc)) / 6.0;
    double s_hi = (4.0 + sqrt(disc)) / 6.0;
    double f_lo = sqrt(fmax(s_lo, 0.0));
    double f_fv = sqrt(s_hi);

    /* Find bracket: a = undershoot, b = overshoot */
    double a = f_lo + 0.01 * (f_fv - f_lo);
    double b = f_fv - 0.001 * (f_fv - f_lo);

    /* Verify a is undershoot */
    shoot(a, omega, N, rmax);
    double fmin_a;
    int sa = zero_crossing(N, &fmin_a);

    /* Verify b is overshoot */
    shoot(b, omega, N, rmax);
    double fmin_b;
    int sb = zero_crossing(N, &fmin_b);

    if (sa == sb) {
        /* Both same sign — scan to find transition */
        int nscan = 500;
        double prev_f0 = a;
        int prev_s = sa;
        for (int j = 1; j <= nscan; j++) {
            double f0 = a + j * (b - a) / nscan;
            shoot(f0, omega, N, rmax);
            double fmin_j;
            int sj = zero_crossing(N, &fmin_j);
            if (sj != prev_s) {
                a = prev_f0;
                b = f0;
                sa = prev_s;
                sb = sj;
                goto bisect;
            }
            prev_f0 = f0;
            prev_s = sj;
        }
        /* No transition found — both endpoints same sign */
        if (sa > 0) {
            /* Both undershoot — try extending b closer to f_fv */
            b = f_fv - 1e-6;
            shoot(b, omega, N, rmax);
            sb = zero_crossing(N, &fmin_b);
            if (sb > 0) {
                fprintf(stderr, "  omega=%.4f: no overshoot found (f_fv=%.6f)\n", omega, f_fv);
                /* Return best undershoot (closest to crossing) */
                double best_f0 = a, best_fmin = 1e10;
                for (int j = 0; j <= 200; j++) {
                    double f0 = a + j * (b - a) / 200;
                    shoot(f0, omega, N, rmax);
                    double fm;
                    zero_crossing(N, &fm);
                    if (fm < best_fmin) { best_fmin = fm; best_f0 = f0; }
                }
                fprintf(stderr, "  Using best f0=%.8f (fmin=%.3e)\n", best_f0, best_fmin);
                shoot(best_f0, omega, N, rmax);
                return best_f0;
            }
        }
        if (sa > 0 && sb < 0) goto bisect;
        fprintf(stderr, "  omega=%.4f: bracketing failed\n", omega);
        return -1;
    }

    /* Ensure a is undershoot (+1), b is overshoot (-1) */
    if (sa < 0 && sb > 0) {
        double tmp = a; a = b; b = tmp;
        int tmps = sa; sa = sb; sb = tmps;
    }

bisect:
    /* Bisect: find f0 at the zero-crossing transition */
    for (int iter = 0; iter < 200; iter++) {
        double mid = 0.5 * (a + b);
        if (b - a < 1e-14 * mid) break;

        shoot(mid, omega, N, rmax);
        double fmin_m;
        int sm = zero_crossing(N, &fmin_m);

        if (sm > 0)
            a = mid;  /* undershoot */
        else
            b = mid;  /* overshoot */
    }

    double f0 = 0.5 * (a + b);
    shoot(f0, omega, N, rmax);
    return f0;
}

/*
 * Compute Q-ball properties from stored profile.
 */
static void compute_props(double omega, int N, double rmax,
                          double *E_out, double *Q_out, double *Rrms_out,
                          double *E_kin, double *E_grad, double *E_pot)
{
    double h = rmax / N;
    double e_sum = 0, q_sum = 0, i2_sum = 0, i4_sum = 0;
    double ek_sum = 0, eg_sum = 0, ep_sum = 0;
    double w2 = omega * omega;

    for (int i = 0; i <= N; i++) {
        double r = i * h;
        double f = f_arr[i], fp = fp_arr[i];
        double f2 = f * f, r2 = r * r;

        double e_kin_d, e_grad_d, e_pot_d, e_dens, q_dens;

        if (b_dbi > 0.0) {
            /* DBI: T_00 = 2*omega^2*f^2/Gamma + b^2*(Gamma-1) + V(f^2)
             *      j^0  = 2*omega*f^2/Gamma
             * where Gamma = sqrt(1 + 2*X/b^2), X = f'^2 - omega^2*f^2
             * Derived from L = -P(X) - V with P(X) = b^2(Gamma-1). */
            double b2 = b_dbi * b_dbi;
            double X = fp * fp - w2 * f2;
            double Gam = sqrt(1.0 + 2.0 * X / b2);
            e_kin_d = 2.0 * w2 * f2 / Gam;
            e_grad_d = b2 * (Gam - 1.0); /* DBI kinetic contribution */
            e_pot_d = V_of_s(f2);
            e_dens = e_kin_d + e_grad_d + e_pot_d;
            q_dens = 2.0 * omega * f2 / Gam;
        } else {
            /* Standard: T_00 = omega^2*f^2 + f'^2 + V_eff(f^2) */
            e_kin_d = w2 * f2;
            e_grad_d = fp * fp;
            if (ell > 0 && r > 1e-15)
                e_grad_d += ell * (ell + 1) * f2 / (r2);
            e_pot_d = Veff_of_s(f2);
            e_dens = e_kin_d + e_grad_d + e_pot_d;
            q_dens = 2.0 * omega * f2;
        }

        double w = (i == 0 || i == N) ? 0.5 : 1.0;
        e_sum += w * e_dens * r2 * h;
        q_sum += w * q_dens * r2 * h;
        i2_sum += w * f2 * r2 * h;
        i4_sum += w * f2 * r2 * r2 * h;
        ek_sum += w * e_kin_d * r2 * h;
        eg_sum += w * e_grad_d * r2 * h;
        ep_sum += w * e_pot_d * r2 * h;
    }

    /* ell=0: angular integral is 4*pi (no harmonic normalization).
     * ell>0: phi = f(r)*Y_l^0, angular integral of |Y_l^0|^2 = 1. */
    double ang_fac = (ell == 0) ? 4.0 * M_PI : 1.0;
    *E_out = ang_fac * e_sum;
    *Q_out = ang_fac * q_sum;
    *Rrms_out = (i2_sum > 0) ? sqrt(i4_sum / i2_sum) : 0;
    *E_kin = ang_fac * ek_sum;
    *E_grad = ang_fac * eg_sum;
    *E_pot = ang_fac * ep_sum;
}

/*
 * Diagnose mode: scan f0 for given omega, print delta at each f0.
 */
static void run_diagnose(double omega, int N, double rmax)
{
    double kappa = sqrt(MASS * MASS - omega * omega);
    double disc = 16.0 - 12.0 * (1.0 - omega * omega);
    double s_lo = (4.0 - sqrt(disc)) / 6.0;
    double s_hi = (4.0 + sqrt(disc)) / 6.0;
    double f_lo = sqrt(fmax(s_lo, 0.0));
    double f_fv = sqrt(s_hi);

    printf("# Diagnostic scan: omega=%.4f, kappa=%.4f, f_lo=%.4f, f_fv=%.4f\n",
           omega, kappa, f_lo, f_fv);
    printf("# %12s %14s %14s %6s\n",
           "f0", "f_min", "f(Rmax)", "cross");

    int nscan = 100;
    double a = f_lo * 0.5;
    double b = f_fv * 1.001;

    for (int j = 0; j <= nscan; j++) {
        double f0 = a + j * (b - a) / nscan;
        int div = shoot(f0, omega, N, rmax);

        double fmin;
        int cross = 0;
        if (!div) cross = zero_crossing(N, &fmin);
        else fmin = -999;

        printf("  %12.6f %14.6e %14.6e %6d %s\n",
               f0, fmin, f_arr[N], cross,
               div ? "DIVERGED" : "");
    }
}

int main(int argc, char **argv)
{
    int N = 20000;
    double rmax = 50.0;
    char outdir[256] = "data";
    int do_single = 0, do_diagnose = 0;
    double omega_single = 0.5;

    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "-omega") && i + 1 < argc) {
            omega_single = atof(argv[++i]);
            do_single = 1;
        } else if (!strcmp(argv[i], "-diagnose") && i + 1 < argc) {
            omega_single = atof(argv[++i]);
            do_diagnose = 1;
        } else if (!strcmp(argv[i], "-rmax") && i + 1 < argc) {
            rmax = atof(argv[++i]);
        } else if (!strcmp(argv[i], "-N") && i + 1 < argc) {
            N = atoi(argv[++i]);
            if (N >= NR_MAX) N = NR_MAX - 1;
        } else if (!strcmp(argv[i], "-b") && i + 1 < argc) {
            b_dbi = atof(argv[++i]);
        } else if (!strcmp(argv[i], "-l") && i + 1 < argc) {
            ell = atoi(argv[++i]);
        } else if (!strcmp(argv[i], "-outdir") && i + 1 < argc) {
            strncpy(outdir, argv[++i], sizeof(outdir) - 1);
        }
    }

    if (b_dbi > 0 && ell > 0) {
        fprintf(stderr, "DBI + ell>0 not supported (angular-averaged DBI not derived)\n");
        return 1;
    }

    if (do_diagnose) {
        run_diagnose(omega_single, N, rmax);
        return 0;
    }

    if (do_single) {
        double f0 = find_f0(omega_single, N, rmax);
        if (f0 < 0) { fprintf(stderr, "Failed\n"); return 1; }

        double E, Q, Rrms, Ek, Eg, Ep;
        compute_props(omega_single, N, rmax, &E, &Q, &Rrms, &Ek, &Eg, &Ep);

        printf("%s Q-ball (ell=%d): V(s) = s(1-s)^2, m = 1\n",
               b_dbi > 0 ? "DBI" : "Scalar", ell);
        if (b_dbi > 0) printf("  b_dbi = %.6f\n", b_dbi);
        printf("  omega = %.8f\n", omega_single);
        if (ell > 0)
            printf("  f'(0) = %.12f  (slope, shooting param)\n", f0);
        else
            printf("  f(0)  = %.12f\n", f0);
        printf("  E     = %.6f  (E_kin=%.4f  E_grad=%.4f  E_pot=%.4f)\n", E, Ek, Eg, Ep);
        printf("  Q     = %.6f\n", Q);
        printf("  E/Q   = %.8f  (stable if < %.1f)\n", Q > 0 ? E / Q : 0, MASS);
        printf("  R_rms = %.6f\n", Rrms);
        double kappa = sqrt(MASS * MASS - omega_single * omega_single);
        printf("  kappa = %.6f  (decay len = %.4f)\n", kappa, 1.0 / kappa);

        if (b_dbi > 0) {
            /* DBI effective speed at center */
            double Gam0_sq = 1.0 - 2.0 * omega_single * omega_single
                             * f0 * f0 / (b_dbi * b_dbi);
            double v_center = sqrt(fmax(Gam0_sq, 0.0));
            printf("  v(0)/c  = %.6f  (DBI speed at center)\n", v_center);
            printf("  Phi/c^2 = %.6f  (potential well depth)\n",
                   (1.0 - Gam0_sq) / 2.0);

            /* Compute speed profile v(r) */
            double h = rmax / N;
            printf("\n# DBI effective speed profile:\n");
            printf("# %10s %12s %12s %12s\n", "r", "f", "v/c", "Phi/c^2");
            int step = (N > 200) ? N / 200 : 1;
            for (int i = 0; i <= N; i += step) {
                double r = i * h;
                double fi = f_arr[i], gi = fp_arr[i];
                double Xi = gi * gi - omega_single * omega_single * fi * fi;
                double G2 = 1.0 + 2.0 * Xi / (b_dbi * b_dbi);
                if (G2 < 0) G2 = 0;
                double vi = sqrt(G2); /* Gamma = v/c for perturbations */
                double Phi = (1.0 - G2) / 2.0;
                printf("  %10.4f %12.6e %12.6f %12.6f\n", r, fi, vi, Phi);
            }
        }

        char fname[512];
        snprintf(fname, sizeof(fname), "%s/qball_profile_om%.3f%s%s.dat",
                 outdir, omega_single,
                 b_dbi > 0 ? "_dbi" : "",
                 ell > 0 ? "_l1" : "");
        FILE *fp = fopen(fname, "w");
        if (fp) {
            double h = rmax / N;
            fprintf(fp, "# r f fp\n");
            int step = (N > 2000) ? N / 2000 : 1;
            for (int i = 0; i <= N; i += step)
                fprintf(fp, "%.6f %.12e %.12e\n", i * h, f_arr[i], fp_arr[i]);
            fclose(fp);
            printf("  Profile: %s\n", fname);
        }

        /* Phase 4: Form factors F_ch(q^2) and F_M(q^2) */
        {
            double h_ff = rmax / N;
            double w2_ff = omega_single * omega_single;
            double b2 = (b_dbi > 0) ? b_dbi * b_dbi : 0;

            /* Compute normalization and radii */
            double ch_norm = 0, m_norm = 0;
            double ch_r2 = 0, m_r2 = 0, ch_r4 = 0, m_r4 = 0;
            for (int i = 0; i <= N; i++) {
                double r = i * h_ff;
                double fi = f_arr[i], gi = fp_arr[i];
                double f2i = fi * fi, r2 = r * r;
                double rho_ch, rho_m;

                if (b_dbi > 0.0) {
                    double Xi = gi * gi - w2_ff * f2i;
                    double Gam = sqrt(1.0 + 2.0 * Xi / b2);
                    rho_ch = 2.0 * omega_single * f2i / Gam;
                    rho_m = 2.0 * w2_ff * f2i / Gam + b2 * (Gam - 1.0) + V_of_s(f2i);
                } else {
                    rho_ch = 2.0 * omega_single * f2i;
                    rho_m = w2_ff * f2i + gi * gi + Veff_of_s(f2i);
                    if (ell > 0 && r > 1e-15)
                        rho_m += ell * (ell + 1) * f2i / r2;
                }

                double wt = (i == 0 || i == N) ? 0.5 : 1.0;
                ch_norm += wt * rho_ch * r2 * h_ff;
                m_norm  += wt * rho_m  * r2 * h_ff;
                ch_r2   += wt * rho_ch * r2 * r2 * h_ff;
                m_r2    += wt * rho_m  * r2 * r2 * h_ff;
                ch_r4   += wt * rho_ch * r2 * r2 * r2 * h_ff;
                m_r4    += wt * rho_m  * r2 * r2 * r2 * h_ff;
            }

            double R_ch = (ch_norm > 0) ? sqrt(ch_r2 / ch_norm) : 0;
            double R_M  = (m_norm > 0)  ? sqrt(m_r2 / m_norm) : 0;
            double R4_ch = (ch_norm > 0) ? ch_r4 / ch_norm : 0;
            double R4_M  = (m_norm > 0)  ? m_r4 / m_norm : 0;

            printf("\n# Phase 4: Form factors\n");
            printf("  R_ch  = %.6f  (charge radius)\n", R_ch);
            printf("  R_M   = %.6f  (mass radius)\n", R_M);
            printf("  <r^4>_ch = %.6f  <r^4>_M = %.6f\n", R4_ch, R4_M);

            /* Form factor output */
            char ff_fname[512];
            snprintf(ff_fname, sizeof(ff_fname),
                     "%s/qball_formfactor_om%.3f%s%s.dat",
                     outdir, omega_single,
                     b_dbi > 0 ? "_dbi" : "",
                     ell > 0 ? "_l1" : "");
            FILE *ff_fp = fopen(ff_fname, "w");
            if (ff_fp) fprintf(ff_fp, "# q F_ch F_M\n");

            printf("# %10s %14s %14s\n", "q", "F_ch(q)", "F_M(q)");
            for (int jq = 1; jq <= 100; jq++) {
                double q = jq * 0.1;
                double fch_sum = 0, fm_sum = 0;

                for (int i = 0; i <= N; i++) {
                    double r = i * h_ff;
                    double fi = f_arr[i], gi = fp_arr[i];
                    double f2i = fi * fi, r2 = r * r;
                    double rho_ch, rho_m;

                    if (b_dbi > 0.0) {
                        double Xi = gi * gi - w2_ff * f2i;
                        double Gam = sqrt(1.0 + 2.0 * Xi / b2);
                        rho_ch = 2.0 * omega_single * f2i / Gam;
                        rho_m = 2.0 * w2_ff * f2i / Gam
                                + b2 * (Gam - 1.0) + V_of_s(f2i);
                    } else {
                        rho_ch = 2.0 * omega_single * f2i;
                        rho_m = w2_ff * f2i + gi * gi + Veff_of_s(f2i);
                        if (ell > 0 && r > 1e-15)
                            rho_m += ell * (ell + 1) * f2i / r2;
                    }

                    double qr = q * r;
                    double j0 = (qr > 1e-10) ? sin(qr) / qr
                                : 1.0 - qr * qr / 6.0;
                    double wt = (i == 0 || i == N) ? 0.5 : 1.0;
                    fch_sum += wt * rho_ch * j0 * r2 * h_ff;
                    fm_sum  += wt * rho_m  * j0 * r2 * h_ff;
                }

                double Fch = (ch_norm > 0) ? fch_sum / ch_norm : 0;
                double Fm  = (m_norm > 0)  ? fm_sum / m_norm : 0;
                printf("  %10.4f %14.8f %14.8f\n", q, Fch, Fm);
                if (ff_fp) fprintf(ff_fp, "%.4f %.10f %.10f\n", q, Fch, Fm);
            }
            if (ff_fp) {
                fclose(ff_fp);
                printf("  Form factors: %s\n", ff_fname);
            }

            /* Phase 4: Bohr-Sommerfeld and stability */
            printf("\n# Phase 4: Stability and quantization\n");
            printf("  Binding: E - Q*m = %.6f  (%.1f%% of Q*m)\n",
                   E - Q * MASS, 100.0 * (E / Q - MASS));
            printf("  n = Q = %.1f  (semiclassical: large n -> excellent WKB)\n", Q);
            printf("  dE/dQ ~ omega = %.4f  (slope of mass-charge curve)\n",
                   omega_single);
        }

    } else {
        /* Dense scan over omega: thick-wall + thin-wall */
        printf("# %s Q-ball scan (ell=%d): V(s) = s(1-s)^2, m = 1\n",
               b_dbi > 0 ? "DBI" : "Scalar", ell);
        printf("# %8s %12s %14s %14s %12s %10s %10s %10s %10s\n",
               "omega", "f(0)", "E", "Q", "E/Q", "R_rms", "E_kin", "E_grad", "E_pot");

        char fname[512];
        snprintf(fname, sizeof(fname), "%s/qball_scan%s.dat",
                 outdir, ell > 0 ? "_l1" : "");
        FILE *fp = fopen(fname, "w");
        if (fp) fprintf(fp, "# omega f0 E Q E/Q R_rms E_kin E_grad E_pot\n");

        /* Dense scan from 0.60 to 0.99 */
        double omegas[] = {0.60, 0.65, 0.70, 0.72, 0.74, 0.76, 0.78,
                           0.80, 0.81, 0.82, 0.83, 0.84, 0.85, 0.86, 0.87, 0.88,
                           0.90, 0.92, 0.94, 0.96, 0.98, 0.99};
        int nom = sizeof(omegas) / sizeof(omegas[0]);

        for (int j = 0; j < nom; j++) {
            double om = omegas[j];
            /* Adaptive Rmax: keep kappa*Rmax ~ 30 for numerical stability.
             * Growing mode amplification e^{kappa*Rmax} must be < 10^{15}.
             * For high omega (small kappa), need larger Rmax for Q-ball to fit.
             * For low omega (large kappa), Rmax=50 is fine but Q-ball may be
             * larger than that (thin-wall) — mark these as unreliable. */
            double kappa = sqrt(1.0 - om * om);
            double rm = fmin(fmax(30.0 / kappa, 50.0), rmax);

            double f0 = find_f0(om, N, rm);
            if (f0 < 0) {
                printf("  %8.4f FAILED\n", om);
                continue;
            }

            double E, Q, Rrms, Ek, Eg, Ep;
            compute_props(om, N, rm, &E, &Q, &Rrms, &Ek, &Eg, &Ep);
            double EQ = (Q > 0) ? E / Q : 0;

            printf("  %8.4f %12.6f %14.4f %14.4f %12.6f %10.4f %10.4f %10.4f %10.4f\n",
                   om, f0, E, Q, EQ, Rrms, Ek, Eg, Ep);
            fflush(stdout);
            if (fp)
                fprintf(fp, "%.6f %.10f %.6f %.6f %.8f %.6f %.6f %.6f %.6f\n",
                        om, f0, E, Q, EQ, Rrms, Ek, Eg, Ep);
        }

        if (fp) { fclose(fp); printf("\n# Data: %s\n", fname); }
        printf("\n# Stability: E/Q < m = 1 means absolutely stable\n");
        printf("# E/Q -> omega for large Q (thin-wall) -> always stable for omega < 1\n");
    }

    return 0;
}
