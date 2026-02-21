/*
 * ed_soler.c — Einstein-Dirac-Soler solver (self-gravitating nonlinear Dirac)
 *
 * System: 4 coupled first-order ODEs for (α, β, A, T)
 * with algebraic scalar bilinear V = N·T·(α²-β²)/r²
 *
 * Metric:  ds² = -T⁻²dt² + A⁻¹dr² + r²dΩ²   (Finster convention)
 * Spinor:  ψ = e^{-iωt} (√T/r) [α(r) χ, iβ(r) χ']
 *
 * Dirac:  √A·α' = (N/(2r))·α - (ωT + m_eff)·β
 *         √A·β' = (ωT - m_eff)·α - (N/(2r))·β
 *
 * Einstein:
 *   rA' = 1 - A - 8πGNωT²(α²+β²) - 4πGr²λV²
 *   2rA(T'/T) = A-1 - 8πGNT√A(αβ'-βα') + 4πGr²λV²
 *
 * where m_eff = m - λV,  V = N·T·(α²-β²)/r².
 *
 * Convention: λ here is the Soler coupling in V = N·T·(α²-β²)/r²
 *   For single-fermion flat-space (soler.c), use λ_ed = λ_soler / N.
 *   Default: m=1, λ=0.5 (matches soler.c's λ=1 for N=2), N=2, G=0.
 *
 * Refs: Zhang & He, arXiv:2508.11714 (2025)
 *       Finster, Smoller, Yau, gr-qc/9801079 (1999)
 *
 * Usage:
 *   ed_soler -omega 0.9                  (single freq, flat space)
 *   ed_soler -omega 0.9 -G 1.0          (with gravity)
 *   ed_soler -scan                       (scan omega)
 *   ed_soler -scan -G 0.01              (scan with weak gravity)
 */

#define _DEFAULT_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define NMAX 200001

static double al_arr[NMAX];  /* α(r) */
static double be_arr[NMAX];  /* β(r) */
static double A_arr[NMAX];   /* A(r)  metric function */
static double T_arr[NMAX];   /* T(r)  metric function */

/* Parameters */
static double par_m     = 1.0;    /* fermion mass */
static double par_omega = 0.9;    /* frequency */
static double par_lambda = 0.5;   /* Soler coupling (= λ_soler/N) */
static double par_G     = 0.0;    /* Newton's constant (0 = flat space) */
static int    par_N     = 2;      /* number of fermions in filled shell */

/*
 * RHS of the coupled ODE system.
 * y[0]=α, y[1]=β, y[2]=A, y[3]=T
 * dy[0..3] = derivatives.
 */
static void rhs(double r, const double *y, double *dy)
{
    double alpha = y[0], beta = y[1], A = y[2], T = y[3];
    double sqrtA = sqrt(A);
    double N = par_N;

    /* Scalar bilinear and effective mass */
    double al2 = alpha * alpha, be2 = beta * beta;
    double V = N * T * (al2 - be2) / (r * r);
    double m_eff = par_m - par_lambda * V;

    /* Dirac equations */
    double dalpha = ((N * 0.5 / r) * alpha - (par_omega * T + m_eff) * beta) / sqrtA;
    double dbeta  = ((par_omega * T - m_eff) * alpha - (N * 0.5 / r) * beta) / sqrtA;

    dy[0] = dalpha;
    dy[1] = dbeta;

    /* Einstein equations */
    double G8piN = 8.0 * M_PI * par_G * N;
    double Vsq = V * V;

    /* A equation: rA' = 1 - A - 8πGNωT²(α²+β²) - 4πGr²λV² */
    dy[2] = (1.0 - A
             - G8piN * par_omega * T * T * (al2 + be2)
             - 4.0 * M_PI * par_G * r * r * par_lambda * Vsq) / r;

    /* T equation: 2rA(T'/T) = A-1 - 8πGNT√A(αβ'-βα') + 4πGr²λV²
     * We use α' and β' already computed above. */
    double cross = alpha * dbeta - beta * dalpha;  /* αβ'-βα' */
    dy[3] = T * (A - 1.0
                 - G8piN * T * sqrtA * cross
                 + 4.0 * M_PI * par_G * r * r * par_lambda * Vsq)
            / (2.0 * r * A);
}

/*
 * RK4 integration from r_start to r_max.
 * Initial conditions set from Taylor expansion at r_start.
 *
 * Taylor (N=2, T₀=1):
 *   α = α₁r, β = β₂r², A = 1 + A₂r², T = 1 + T₂r²
 *   where β₂ = (ω - m + 2λα₁²)α₁/3
 *         A₂ = -(16/3)πGα₁²(ω + λα₁²)     [with proper N=2 factors]
 *         T₂ = (4/3)πGα₁²(m - 2ω)
 */
static int integrate(double alpha1, int Ngrid, double rmax)
{
    double r_start = 1e-5;
    double h = (rmax - r_start) / Ngrid;
    double T0 = 1.0;  /* gauge choice */
    int N = par_N;

    /* Taylor expansion at r = r_start */
    double a1sq = alpha1 * alpha1;
    double V0 = N * T0 * a1sq;  /* V at r→0 (leading order) */
    double m_eff0 = par_m - par_lambda * V0;
    double beta2 = (par_omega * T0 - m_eff0) * alpha1 / (N + 1.0);

    /* A and T corrections (only needed for G > 0) */
    /* For N=2: A₂ and T₂ derived from the ODEs at r→0 */
    double A2 = 0, T2 = 0;
    if (par_G > 0) {
        /* From rA'=1-A-... at leading order:
         * 3A₂ = -8πGN·ω·T₀²·α₁² - 4πG·λ·V₀² */
        A2 = -(8.0 * M_PI * par_G * N * par_omega * T0 * T0 * a1sq
               + 4.0 * M_PI * par_G * par_lambda * V0 * V0) / 3.0;

        /* T₂ from the T equation at leading order */
        /* Using: cross₀ = αβ'-βα' at r→0:
         * α=α₁r, α'=α₁, β=β₂r², β'=2β₂r
         * cross = α₁r·2β₂r - β₂r²·α₁ = α₁β₂r²
         * So at leading order: 8πGNT√A·cross → 8πGN·α₁β₂r² */
        double term_cross = 8.0 * M_PI * par_G * N * T0 * alpha1 * beta2;
        double term_A = A2;  /* A-1 at leading order = A₂r² */
        double term_V = 4.0 * M_PI * par_G * par_lambda * V0 * V0;
        /* 4T₂ = A₂ - 8πGN·α₁β₂ + 4πGλV₀² */
        T2 = (term_A - term_cross + term_V) / 4.0;
    }

    double rs = r_start;
    al_arr[0] = alpha1 * rs;
    be_arr[0] = beta2 * rs * rs;
    A_arr[0]  = 1.0 + A2 * rs * rs;
    T_arr[0]  = T0 + T2 * rs * rs;

    for (int i = 0; i < Ngrid; i++) {
        double r = r_start + i * h;
        double y[4] = {al_arr[i], be_arr[i], A_arr[i], T_arr[i]};
        double k1[4], k2[4], k3[4], k4[4], yt[4];

        rhs(r, y, k1);

        for (int j = 0; j < 4; j++) yt[j] = y[j] + 0.5 * h * k1[j];
        rhs(r + 0.5 * h, yt, k2);

        for (int j = 0; j < 4; j++) yt[j] = y[j] + 0.5 * h * k2[j];
        rhs(r + 0.5 * h, yt, k3);

        for (int j = 0; j < 4; j++) yt[j] = y[j] + h * k3[j];
        rhs(r + h, yt, k4);

        for (int j = 0; j < 4; j++)
            yt[j] = y[j] + (h / 6.0) * (k1[j] + 2*k2[j] + 2*k3[j] + k4[j]);

        al_arr[i+1] = yt[0];
        be_arr[i+1] = yt[1];
        A_arr[i+1]  = yt[2];
        T_arr[i+1]  = yt[3];

        /* Bail on divergence */
        if (fabs(al_arr[i+1]) > 1e10 || fabs(be_arr[i+1]) > 1e10 ||
            A_arr[i+1] < 1e-6 || T_arr[i+1] < 1e-6 || T_arr[i+1] > 1e6)
            return i + 1;
    }
    return Ngrid;
}

/*
 * Classify solution: -1 if f=α/r crosses zero in the physical region
 * (α₁ too large), +1 if f stays positive (α₁ too small).
 *
 * Use f²+g² = (α/r)²+(β/r)² which starts at maximum and decays,
 * rather than α²+β² which starts at zero (would always give imin=0).
 */
static int classify(int nv, double rmax)
{
    double r_start = 1e-5;
    double h = (rmax - r_start) / nv;

    /* Find minimum of f²+g² after peak (defines end of physical region) */
    double amp_min = 1e30;
    int imin = nv;
    for (int i = 1; i <= nv; i++) {
        double r = r_start + i * h;
        double f = al_arr[i] / r, g = be_arr[i] / r;
        double amp = f * f + g * g;
        if (amp < amp_min) {
            amp_min = amp;
            imin = i;
        }
    }

    /* Check for f = α/r crossing zero in physical region */
    for (int i = 1; i <= imin; i++) {
        if (al_arr[i] < 0) return -1;
    }
    return +1;
}

/*
 * Find α₁ by bisection for current par_omega.
 */
static double find_alpha1(int Ngrid, double rmax)
{
    double last_pos = -1;
    double first_neg = -1;
    int found_pos = 0;

    for (double a_scan = 0.001; a_scan <= 30.0; a_scan *= 1.01) {
        int nv = integrate(a_scan, Ngrid, rmax);
        int s = classify(nv, rmax);
        if (s == +1) {
            last_pos = a_scan;
            found_pos = 1;
        } else if (s == -1 && found_pos) {
            first_neg = a_scan;
            break;
        }
    }

    if (last_pos < 0 || first_neg < 0) return -1;

    double a_lo = last_pos;
    double a_hi = first_neg;

    /* Bisection: 80 iterations for ~24 digits */
    for (int iter = 0; iter < 80; iter++) {
        double a_mid = 0.5 * (a_lo + a_hi);
        int nv = integrate(a_mid, Ngrid, rmax);
        int s = classify(nv, rmax);
        if (s == +1)
            a_lo = a_mid;
        else
            a_hi = a_mid;
    }

    double a1_best = 0.5 * (a_lo + a_hi);
    integrate(a1_best, Ngrid, rmax);

    /* Clip growing tail: use f²+g² = (α/r)²+(β/r)² */
    double r_start = 1e-5;
    double h = (rmax - r_start) / Ngrid;
    double amp_min = 1e30;
    int imin = Ngrid;
    for (int i = 1; i <= Ngrid; i++) {
        double r = r_start + i * h;
        double f = al_arr[i] / r, g = be_arr[i] / r;
        double amp = f * f + g * g;
        if (amp < amp_min) {
            amp_min = amp;
            imin = i;
        }
    }
    for (int i = imin + 1; i <= Ngrid; i++) {
        al_arr[i] = 0;
        be_arr[i] = 0;
    }

    return a1_best;
}

/*
 * Compute physical observables from the solution.
 *
 *  Q_norm  = 4π ∫ (α²+β²) T/√A dr   (normalization / particle number)
 *  R_rms   = √(∫r² ρ dr / ∫ρ dr)      where ρ = (α²+β²) T/√A
 *  M_ADM   = (1-A(r_max))·r_max/(2G)  (ADM mass)
 *  Phi_c   = 1 - 1/T(0)²              (central gravitational redshift)
 */
static void compute_props(int Ngrid, double rmax,
                          double *Q_out, double *R_out,
                          double *M_ADM_out, double *Phi_c_out,
                          double *E_self_out)
{
    double r_start = 1e-5;
    double h = (rmax - r_start) / Ngrid;
    double Q = 0, R2w = 0, Eself = 0;

    for (int i = 0; i <= Ngrid; i++) {
        double r = r_start + i * h;
        double al2 = al_arr[i] * al_arr[i];
        double be2 = be_arr[i] * be_arr[i];

        /* Skip clipped tail (α=β=0) */
        if (al2 + be2 < 1e-30) continue;

        double A = A_arr[i], T = T_arr[i];
        if (A < 1e-10 || !isfinite(A) || !isfinite(T)) continue;
        double sqrtA = sqrt(A);
        double rho = (al2 + be2) * T / sqrtA;

        /* Soler self-energy: 2πλ∫r²V²dr where V = N·T·(α²-β²)/r² */
        double V_sc = par_N * T * (al2 - be2) / (r * r);
        double eself = r * r * V_sc * V_sc;  /* r²V² */

        /* Simpson weights */
        double w;
        if (i == 0 || i == Ngrid)
            w = 1.0 / 3.0;
        else if (i % 2 == 1)
            w = 4.0 / 3.0;
        else
            w = 2.0 / 3.0;

        Q     += w * rho * h;
        R2w   += w * rho * r * r * h;
        Eself += w * eself * h;
    }

    *Q_out = 4.0 * M_PI * Q;
    *R_out = (Q > 0) ? sqrt(R2w / Q) : 0;
    *E_self_out = 2.0 * M_PI * par_lambda * Eself;

    /* ADM mass from A(r_max): A → 1 - 2GM/r */
    if (par_G > 0 && Ngrid > 0) {
        double r_end = rmax;
        double A_end = A_arr[Ngrid];
        *M_ADM_out = (1.0 - A_end) * r_end / (2.0 * par_G);
    } else {
        *M_ADM_out = 0;
    }

    /* Central potential: Φ/c² = 1 - T_∞²/T₀² (redshift between center and infinity) */
    double T_inf = T_arr[Ngrid];
    *Phi_c_out = 1.0 - (T_inf * T_inf) / (T_arr[0] * T_arr[0]);
}

static void run_single(double omega, double rmax, int Ngrid, const char *outdir)
{
    par_omega = omega;
    double kappa = sqrt(par_m * par_m - par_omega * par_omega);

    if (10.0 / kappa > rmax) rmax = 10.0 / kappa;
    if (rmax > 500) rmax = 500;

    printf("# Einstein-Dirac-Soler: omega=%.6f, m=%.6f, lambda=%.6f, G=%.6e, N=%d\n",
           par_omega, par_m, par_lambda, par_G, par_N);
    printf("# rmax=%.1f, Ngrid=%d, h=%.6f, kappa=%.6f\n",
           rmax, Ngrid, (rmax - 1e-5) / Ngrid, kappa);

    double alpha1 = find_alpha1(Ngrid, rmax);
    if (alpha1 < 0) {
        printf("# FAILED: No bracket found for omega=%.6f\n", par_omega);
        return;
    }

    double Q, R, M_ADM, Phi_c, E_self;
    compute_props(Ngrid, rmax, &Q, &R, &M_ADM, &Phi_c, &E_self);

    double M_total = par_N * par_omega * Q + E_self; /* N·ω·Q + self-energy */
    double T_inf = T_arr[Ngrid];
    double omega_phys = par_omega * T_inf;  /* physical frequency (rescaled to T_∞=1) */

    printf("#\n# RESULT: Einstein-Dirac-Soler soliton\n");
    printf("#   alpha1    = %.8f  (spinor amplitude at origin)\n", alpha1);
    printf("#   omega     = %.8f  (coordinate freq, gauge T₀=1)\n", par_omega);
    printf("#   omega_phys= %.8f  (physical freq, rescaled T_∞=1)\n", omega_phys);
    printf("#   E_bind    = %.8f  (ω_phys - m)\n", omega_phys - par_m);
    printf("#   Q_norm    = %.6f  (per-orbital charge)\n", Q);
    printf("#   M_total   = %.6f  (N·ω·Q + E_self)\n", M_total);
    printf("#   E_self    = %.6f  (Soler self-energy 2πλ∫r²V²dr)\n", E_self);
    printf("#   R_rms     = %.6f\n", R);
    printf("#   A(rmax)   = %.10f\n", A_arr[Ngrid]);
    printf("#   T(rmax)   = %.10f\n", T_inf);
    printf("#   T(0)      = %.10f\n", T_arr[0]);

    if (par_G > 0) {
        printf("#   M_ADM     = %.8f  ((1-A(rmax))·rmax/(2G))\n", M_ADM);
        printf("#   Phi_c/c2  = %.8e  (1 - T_∞²/T₀²)\n", Phi_c);
        printf("#   dM/M      = %.6e  ((M_ADM-M_total)/M_total)\n",
               M_total > 0 ? (M_ADM - M_total) / M_total : 0);
        printf("#   GM/R      = %.6e  (compactness)\n",
               par_G * M_ADM / R);
    }

    /* Save profile */
    if (outdir) {
        char fname[256];
        snprintf(fname, sizeof(fname), "%s/ed_soler_om%.3f_G%.1e.dat",
                 outdir, par_omega, par_G);
        FILE *fp = fopen(fname, "w");
        if (fp) {
            double r_start = 1e-5;
            double h = (rmax - r_start) / Ngrid;
            fprintf(fp, "# ED-Soler: omega=%.6f m=%.6f lambda=%.6f G=%.6e N=%d alpha1=%.8f\n",
                    par_omega, par_m, par_lambda, par_G, par_N, alpha1);
            fprintf(fp, "# r  alpha  beta  A  T  rho=(al2+be2)  V=N*T*(al2-be2)/r2\n");
            int step = (Ngrid > 5000) ? Ngrid / 5000 : 1;
            for (int i = 0; i <= Ngrid; i += step) {
                double r = r_start + i * h;
                double al = al_arr[i], be = be_arr[i];
                double rho = al * al + be * be;
                double V_sc = par_N * T_arr[i] * (al*al - be*be) / (r*r);
                fprintf(fp, "%.6f %.10e %.10e %.10e %.10e %.10e %.10e\n",
                        r, al, be, A_arr[i], T_arr[i], rho, V_sc);
            }
            fclose(fp);
            printf("#   Profile: %s\n", fname);
        }
    }
}

static void run_scan(double rmax, int Ngrid, const char *outdir)
{
    double omega_vals[] = {0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99};
    int nvals = sizeof(omega_vals) / sizeof(omega_vals[0]);

    printf("# Einstein-Dirac-Soler scan: m=%.6f, lambda=%.6f, G=%.6e, N=%d\n",
           par_m, par_lambda, par_G, par_N);
    printf("# %10s  %12s  %12s  %8s  %12s  %12s  %12s  %12s  %12s\n",
           "omega", "alpha1", "Q_norm", "R_rms", "M_total",
           "E_self", "M_ADM", "A(rmax)", "T(rmax)");

    for (int k = 0; k < nvals; k++) {
        par_omega = omega_vals[k];
        double rmax_k = rmax;
        double kappa = sqrt(par_m * par_m - par_omega * par_omega);
        if (kappa > 0 && 10.0 / kappa > rmax_k)
            rmax_k = 10.0 / kappa;
        if (rmax_k > 500)
            rmax_k = 500;

        double alpha1 = find_alpha1(Ngrid, rmax_k);
        if (alpha1 < 0) {
            printf("  %10.4f  FAILED\n", par_omega);
            continue;
        }

        double Q, R, M_ADM, Phi_c, E_self;
        compute_props(Ngrid, rmax_k, &Q, &R, &M_ADM, &Phi_c, &E_self);
        double M_total = par_N * par_omega * Q + E_self;

        printf("  %10.4f  %12.6f  %12.6f  %8.4f  %12.6f  %12.6f  %12.6e  %12.8f  %12.8f\n",
               par_omega, alpha1, Q, R, M_total, E_self,
               M_ADM, A_arr[Ngrid], T_arr[Ngrid]);

        /* Save profile */
        if (outdir) {
            double r_start = 1e-5;
            double h = (rmax_k - r_start) / Ngrid;
            char fname[256];
            snprintf(fname, sizeof(fname), "%s/ed_soler_om%.3f_G%.1e.dat",
                     outdir, par_omega, par_G);
            FILE *fp = fopen(fname, "w");
            if (fp) {
                fprintf(fp, "# ED-Soler: omega=%.6f m=%.6f lambda=%.6f G=%.6e alpha1=%.8f\n",
                        par_omega, par_m, par_lambda, par_G, alpha1);
                fprintf(fp, "# r  alpha  beta  A  T  rho  V\n");
                int step = (Ngrid > 5000) ? Ngrid / 5000 : 1;
                for (int i = 0; i <= Ngrid; i += step) {
                    double r = r_start + i * h;
                    double al = al_arr[i], be = be_arr[i];
                    fprintf(fp, "%.6f %.10e %.10e %.10e %.10e %.10e %.10e\n",
                            r, al, be, A_arr[i], T_arr[i],
                            al*al + be*be,
                            par_N * T_arr[i] * (al*al - be*be) / (r*r));
                }
                fclose(fp);
            }
        }
    }
}

static void usage(const char *prog)
{
    fprintf(stderr,
            "Usage: %s [-omega ω] [-m mass] [-lambda λ] [-G G] [-Nferm N]\n"
            "          [-rmax R] [-N n] [-scan] [-outdir dir]\n",
            prog);
}

int main(int argc, char *argv[])
{
    double omega = 0.9;
    double rmax = 50.0;
    int Ngrid = 100000;
    int scan = 0;
    const char *outdir = NULL;

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-omega") == 0 && i + 1 < argc)
            omega = atof(argv[++i]);
        else if (strcmp(argv[i], "-m") == 0 && i + 1 < argc)
            par_m = atof(argv[++i]);
        else if (strcmp(argv[i], "-lambda") == 0 && i + 1 < argc)
            par_lambda = atof(argv[++i]);
        else if (strcmp(argv[i], "-G") == 0 && i + 1 < argc)
            par_G = atof(argv[++i]);
        else if (strcmp(argv[i], "-Nferm") == 0 && i + 1 < argc)
            par_N = atoi(argv[++i]);
        else if (strcmp(argv[i], "-rmax") == 0 && i + 1 < argc)
            rmax = atof(argv[++i]);
        else if (strcmp(argv[i], "-N") == 0 && i + 1 < argc)
            Ngrid = atoi(argv[++i]);
        else if (strcmp(argv[i], "-scan") == 0)
            scan = 1;
        else if (strcmp(argv[i], "-outdir") == 0 && i + 1 < argc)
            outdir = argv[++i];
        else {
            usage(argv[0]);
            return 1;
        }
    }

    if (Ngrid > NMAX - 1) Ngrid = NMAX - 1;

    if (scan) {
        run_scan(rmax, Ngrid, outdir);
    } else {
        run_single(omega, rmax, Ngrid, outdir);
    }

    return 0;
}
