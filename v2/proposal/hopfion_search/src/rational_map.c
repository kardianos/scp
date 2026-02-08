/*
 * rational_map.c — Rational Map Ansatz Skyrmion Solver for B >= 1
 *
 * Generalizes the hedgehog profile solver to higher baryon numbers
 * using the rational map ansatz (Houghton, Manton, Sutcliffe 1998).
 *
 * The ansatz: q-hat = cos f(r) + sin f(r) n-hat(theta,phi) . sigma
 * where n-hat is determined by a degree-B rational map R: S^2 -> S^2.
 *
 * Energy:
 *   E2 = 2*pi*rho0^2 * integral(f'^2*r^2 + 2*B*sin^2(f)) dr
 *   E4 = (4*pi*rho0^4/e^2) * integral(2*B*f'^2*sin^2(f) + I*sin^4(f)/r^2) dr
 *
 * where B = degree of rational map, I = quartic angular integral.
 *
 * ODE:
 *   f'' = [sin(2f)*(B + I*c4*sin^2(f)/r^2) - 2*r*f' - B*c4*f'^2*sin(2f)]
 *         / (r^2 + 2*B*c4*sin^2(f))
 *
 * with c4 = 2*rho0^2/e^2, f(0) = pi, f(inf) = 0.
 *
 * Near-origin behavior:
 *   For the ODE to be regular at r=0, the profile must behave as
 *     f(r) = pi - a * r^alpha + ...
 *   where alpha = (-1 + sqrt(1 + 8*B)) / 2.
 *     B=1: alpha = 1  (linear, standard hedgehog)
 *     B=2: alpha = (-1+sqrt(17))/2 ~ 1.5616
 *     B=3: alpha = 2
 *     B=4: alpha = (-1+sqrt(33))/2 ~ 2.3723
 *
 *   The shooting parameter is 'a', the coefficient of r^alpha.
 *
 * The rational map integral I is computed numerically for known
 * optimal maps:
 *   B=1: R(z) = z                                      (spherical)
 *   B=2: R(z) = z^2                                    (axial)
 *   B=3: R(z) = (z^3 - sqrt(3)*i*z)/(sqrt(3)*i*z^2 - 1)  (tetrahedral)
 *   B=4: R(z) = (z^4 + 2*sqrt(3)*i*z^2 + 1)/(z^4 - 2*sqrt(3)*i*z^2 + 1) (cubic)
 */

#define _USE_MATH_DEFINES
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <complex.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ========== Rational Map Angular Integral I ========== */

/* Evaluate R(z) and R'(z) for the known optimal rational maps */
static void rational_map_B1(double complex z,
                            double complex *R, double complex *dR)
{
    *R = z;
    *dR = 1.0;
}

static void rational_map_B2(double complex z,
                            double complex *R, double complex *dR)
{
    *R = z * z;
    *dR = 2.0 * z;
}

static void rational_map_B3(double complex z,
                            double complex *R, double complex *dR)
{
    double complex s3i = I * sqrt(3.0);
    double complex z2 = z * z;
    double complex z3 = z2 * z;

    double complex num = z3 - s3i * z;
    double complex den = s3i * z2 - 1.0;

    *R = num / den;

    double complex dnum = 3.0 * z2 - s3i;
    double complex dden = 2.0 * s3i * z;

    *dR = (dnum * den - num * dden) / (den * den);
}

static void rational_map_B4(double complex z,
                            double complex *R, double complex *dR)
{
    double complex s3i2 = I * 2.0 * sqrt(3.0);
    double complex z2 = z * z;
    double complex z3 = z2 * z;
    double complex z4 = z2 * z2;

    double complex num = z4 + s3i2 * z2 + 1.0;
    double complex den = z4 - s3i2 * z2 + 1.0;

    *R = num / den;

    double complex dnum = 4.0 * z3 + 2.0 * s3i2 * z;
    double complex dden = 4.0 * z3 - 2.0 * s3i2 * z;

    *dR = (dnum * den - num * dden) / (den * den);
}

typedef void (*rational_map_fn)(double complex z,
                                double complex *R, double complex *dR);

/* Gauss-Legendre nodes and weights for n points on [-1,1] */
static void gauss_legendre(int n, double *nodes, double *weights)
{
    for (int i = 0; i < n; i++) {
        double x = cos(M_PI * (i + 0.75) / (n + 0.5));
        double dx;
        int iter = 0;
        do {
            double p0 = 1.0, p1 = x;
            for (int j = 2; j <= n; j++) {
                double p2 = ((2*j - 1) * x * p1 - (j - 1) * p0) / j;
                p0 = p1;
                p1 = p2;
            }
            double dp = n * (x * p1 - p0) / (x * x - 1.0);
            dx = -p1 / dp;
            x += dx;
            iter++;
        } while (fabs(dx) > 1e-15 && iter < 100);
        nodes[i] = x;
        double p0 = 1.0, p1_val = x;
        for (int j = 2; j <= n; j++) {
            double p2 = ((2*j - 1) * x * p1_val - (j - 1) * p0) / j;
            p0 = p1_val;
            p1_val = p2;
        }
        double dp = n * (x * p1_val - p0) / (x * x - 1.0);
        weights[i] = 2.0 / ((1.0 - x * x) * dp * dp);
    }
}

/*
 * Compute the angular integral I for a given rational map.
 *   I = (1/4pi) * integral_{S^2} b(z)^2 dOmega
 * where b(z) = (1+|z|^2)^2 |R'(z)|^2 / (1+|R(z)|^2)^2
 */
static double compute_angular_integral(rational_map_fn map, int B_degree)
{
    int N_theta = 400;
    int N_phi = 800;

    double *gl_nodes = (double *)malloc(N_theta * sizeof(double));
    double *gl_weights = (double *)malloc(N_theta * sizeof(double));
    gauss_legendre(N_theta, gl_nodes, gl_weights);

    double I_sum = 0.0;
    double B_check = 0.0;

    for (int i = 0; i < N_theta; i++) {
        double cos_theta = gl_nodes[i];
        double w_theta = gl_weights[i];

        double theta = acos(cos_theta);
        double half_theta = 0.5 * theta;
        double tan_half = tan(half_theta);

        double dphi = 2.0 * M_PI / N_phi;

        for (int j = 0; j < N_phi; j++) {
            double phi = j * dphi;

            double complex z = tan_half * cexp(I * phi);
            double abs_z2 = cabs(z) * cabs(z);

            double complex R_val, dR_val;
            map(z, &R_val, &dR_val);

            double abs_R2 = cabs(R_val) * cabs(R_val);
            double abs_dR2 = cabs(dR_val) * cabs(dR_val);

            double factor_z = (1.0 + abs_z2);
            double factor_R = (1.0 + abs_R2);
            double b = factor_z * factor_z * abs_dR2 / (factor_R * factor_R);

            I_sum += b * b * w_theta * dphi;
            B_check += b * w_theta * dphi;
        }
    }

    I_sum /= (4.0 * M_PI);
    B_check /= (4.0 * M_PI);

    printf("  B=%d: I = %.8f, B_check (should be %d) = %.8f\n",
           B_degree, I_sum, B_degree, B_check);

    free(gl_nodes);
    free(gl_weights);

    return I_sum;
}

/* ========== Radial ODE Solver ========== */

/*
 * ODE right-hand side:
 * f'' = [sin(2f)*(B + I*c4*sin^2(f)/r^2) - 2*r*f' - B*c4*f'^2*sin(2f)]
 *       / (r^2 + 2*B*c4*sin^2(f))
 */
static double f_rhs(double r, double f, double fp,
                    double c4, double B, double II)
{
    double sf = sin(f);
    double s2f = sin(2.0 * f);
    double sf2 = sf * sf;
    double r2 = r * r;
    double denom = r2 + 2.0 * B * c4 * sf2;

    if (denom < 1e-30)
        return 0.0;

    /* Guard against r=0 in sin^2(f)/r^2 term:
     * if f ~ pi - a*r^alpha, then sin^2(f)/r^2 ~ a^2*r^(2alpha-2),
     * and sin(2f)*sin^2(f)/r^2 ~ -2*a^3*r^(3alpha-2),
     * which is fine as long as 3*alpha-2 >= 0, i.e., alpha >= 2/3.
     * Since alpha >= 1 for B >= 1, this is always OK. But numerically
     * we still need to avoid 0/0. */
    double sf2_over_r2;
    if (r2 > 1e-30)
        sf2_over_r2 = sf2 / r2;
    else
        sf2_over_r2 = 0.0;

    double num = s2f * (B + II * c4 * sf2_over_r2)
               - 2.0 * fp * r
               - B * c4 * fp * fp * s2f;
    return num / denom;
}

typedef struct {
    int N;
    double dr;
    double *r;
    double *f;
    double *fp;
} Profile;

static Profile *profile_alloc(int N, double R_max)
{
    Profile *p = (Profile *)malloc(sizeof(Profile));
    p->N = N;
    p->dr = R_max / N;
    p->r = (double *)calloc((size_t)(N + 1), sizeof(double));
    p->f = (double *)calloc((size_t)(N + 1), sizeof(double));
    p->fp = (double *)calloc((size_t)(N + 1), sizeof(double));
    for (int i = 0; i <= N; i++)
        p->r[i] = i * p->dr;
    return p;
}

static void profile_free(Profile *p)
{
    if (p) { free(p->r); free(p->f); free(p->fp); free(p); }
}

/*
 * Near-origin behavior: f(r) = pi - a * r^alpha
 * where alpha = (-1 + sqrt(1 + 8*B)) / 2
 *
 * This comes from requiring regularity of the ODE at r=0.
 * Setting g(r) = pi - f(r) ~ a*r^alpha, the leading-order
 * balance in the ODE gives alpha*(alpha+1) = 2*B.
 */
static double origin_exponent(double B)
{
    return 0.5 * (-1.0 + sqrt(1.0 + 8.0 * B));
}

/*
 * Integrate the ODE via RK4 from r ~ r_start with initial conditions
 * determined by the near-origin series:
 *   f(r) = pi - a * r^alpha
 *   f'(r) = -a * alpha * r^(alpha-1)
 *
 * The shooting parameter is 'a' (coefficient of r^alpha).
 *
 * For B=1 (alpha=1): f(r) ~ pi - a*r, f'(r) ~ -a  (standard hedgehog)
 * For B>=2: alpha > 1, so f'(0) = 0.
 *
 * We start integration from r = N_start * dr (a few grid points in)
 * using the series initial conditions, then RK4 the rest.
 */
static double shoot(Profile *p, double a, double c4, double B, double II)
{
    double dr = p->dr;
    double alpha = origin_exponent(B);

    /* Number of initial points to set from series solution */
    int N_start = 4;

    /* Set initial points from series: f(r) = pi - a*r^alpha */
    for (int i = 0; i <= N_start; i++) {
        double r = i * dr;
        p->r[i] = r;
        if (r < 1e-30) {
            p->f[i] = M_PI;
            p->fp[i] = 0.0;  /* f'(0) = 0 for alpha > 1 (B >= 2) */
            if (B <= 1.0 + 1e-10) {
                /* B=1: alpha=1, f'(0) = -a */
                p->fp[i] = -a;
            }
        } else {
            double r_alpha = pow(r, alpha);
            double r_alpha_m1 = pow(r, alpha - 1.0);
            p->f[i] = M_PI - a * r_alpha;
            p->fp[i] = -a * alpha * r_alpha_m1;
        }
    }

    /* RK4 integration from N_start onward */
    for (int i = N_start; i < p->N; i++) {
        double r0 = p->r[i];
        double f0 = p->f[i];
        double v0 = p->fp[i];
        double h = dr;

        double k1f = v0;
        double k1v = f_rhs(r0, f0, v0, c4, B, II);

        double k2f = v0 + 0.5 * h * k1v;
        double k2v = f_rhs(r0 + 0.5*h, f0 + 0.5*h*k1f,
                           v0 + 0.5*h*k1v, c4, B, II);

        double k3f = v0 + 0.5 * h * k2v;
        double k3v = f_rhs(r0 + 0.5*h, f0 + 0.5*h*k2f,
                           v0 + 0.5*h*k2v, c4, B, II);

        double k4f = v0 + h * k3v;
        double k4v = f_rhs(r0 + h, f0 + h*k3f,
                           v0 + h*k3v, c4, B, II);

        p->f[i+1] = f0 + (h / 6.0) * (k1f + 2*k2f + 2*k3f + k4f);
        p->fp[i+1] = v0 + (h / 6.0) * (k1v + 2*k2v + 2*k3v + k4v);
        p->r[i+1] = r0 + h;
    }

    return p->f[p->N];
}

/* ========== Energy Computation ========== */

typedef struct {
    double E2, E4, Etotal, Q;
} Energy;

static Energy compute_energy(const Profile *p, double e_skyrme, double rho0,
                             double B, double II)
{
    double inv_e2 = 1.0 / (e_skyrme * e_skyrme);
    double rho2 = rho0 * rho0;
    double rho4 = rho2 * rho2;
    double dr = p->dr;

    double e2_sum = 0.0, e4_sum = 0.0, q_sum = 0.0;

    for (int i = 0; i <= p->N; i++) {
        double r = p->r[i];
        double fi = p->f[i];
        double fpi = p->fp[i];
        double sf = sin(fi);
        double sf2 = sf * sf;

        double w = dr;
        if (i == 0 || i == p->N) w *= 0.5;

        /* E2 integrand: f'^2 r^2 + 2*B*sin^2(f) */
        e2_sum += (fpi * fpi * r * r + 2.0 * B * sf2) * w;

        /* E4 integrand: 2*B*f'^2*sin^2(f) + I*sin^4(f)/r^2 */
        if (r > 1e-14)
            e4_sum += (2.0 * B * fpi * fpi * sf2
                     + II * sf2 * sf2 / (r * r)) * w;

        /* Topological charge integrand: (-f')*sin^2(f) */
        q_sum += (-fpi * sf2) * w;
    }

    Energy en;
    en.E2 = 2.0 * M_PI * rho2 * e2_sum;
    en.E4 = 4.0 * M_PI * inv_e2 * rho4 * e4_sum;
    en.Etotal = en.E2 + en.E4;
    en.Q = (2.0 * B / M_PI) * q_sum;
    return en;
}

/* ========== Shooting for one baryon number ========== */

typedef struct {
    int B;
    double I_angular;
    double alpha;
    double a_shooting;
    Energy en;
    Profile *prof;
} SkyrmionResult;

static SkyrmionResult solve_for_B(int B_val, double I_val,
                                  double e_skyrme, double rho0,
                                  int Nr, double R_max)
{
    SkyrmionResult res;
    memset(&res, 0, sizeof(res));
    res.B = B_val;
    res.I_angular = I_val;

    double c4 = 2.0 * rho0 * rho0 / (e_skyrme * e_skyrme);
    double alpha = origin_exponent((double)B_val);
    res.alpha = alpha;

    printf("\n===== B = %d, I = %.6f, c4 = %.6f, alpha = %.6f =====\n",
           B_val, I_val, c4, alpha);

    Profile *p = profile_alloc(Nr, R_max);

    /* Search range for shooting parameter a.
     * For B=1: a ~ 5.68 (slope at origin).
     * For higher B: a is the coefficient of r^alpha.
     * Since f must drop from pi to 0, a is of order pi/R_eff^alpha
     * where R_eff is the soliton size. */
    double a_lo = 0.01;
    double a_hi = 200.0;

    /* Scan for bracket: f(R_max) should change sign */
    printf("  Scanning for bracket (alpha=%.4f)...\n", alpha);

    /* Use adaptive scanning with increasing step size */
    double f_prev = shoot(p, a_lo, c4, (double)B_val, I_val);
    double a_prev = a_lo;
    int found = 0;
    double a_bracket_lo = a_lo, a_bracket_hi = a_lo;
    double f_bracket_lo = f_prev, f_bracket_hi = f_prev;

    /* Logarithmic scan to cover wide range */
    int n_scan = 2000;
    double log_lo = log(a_lo), log_hi = log(a_hi);
    double dlog = (log_hi - log_lo) / n_scan;

    for (int k = 1; k <= n_scan; k++) {
        double a = exp(log_lo + k * dlog);
        double fend = shoot(p, a, c4, (double)B_val, I_val);

        if (isnan(fend) || isinf(fend)) {
            a_prev = a;
            f_prev = fend;
            continue;
        }

        if (!isnan(f_prev) && !isinf(f_prev) && f_prev * fend < 0.0) {
            a_bracket_lo = a_prev;
            f_bracket_lo = f_prev;
            a_bracket_hi = a;
            f_bracket_hi = fend;
            found = 1;
            printf("  Found bracket: a in [%.6f, %.6f], "
                   "f_end in [%.6e, %.6e]\n",
                   a_bracket_lo, a_bracket_hi,
                   f_bracket_lo, f_bracket_hi);
            break;
        }

        a_prev = a;
        f_prev = fend;
    }

    if (!found) {
        fprintf(stderr, "  ERROR: Could not find bracket for B=%d\n", B_val);
        res.a_shooting = -1;
        res.prof = p;
        return res;
    }

    double a_lo2 = a_bracket_lo;
    double a_hi2 = a_bracket_hi;
    double f_lo2 = f_bracket_lo;

    /* Bisection */
    printf("  Bisecting...\n");
    double a_best = 0;
    for (int iter = 0; iter < 200; iter++) {
        double a_mid = 0.5 * (a_lo2 + a_hi2);
        double f_mid = shoot(p, a_mid, c4, (double)B_val, I_val);

        if (isnan(f_mid)) {
            a_hi2 = a_mid;
            continue;
        }

        if (iter < 10 || iter % 20 == 0)
            printf("    iter %3d: a=%.14f, f(Rmax)=%+.6e\n",
                   iter, a_mid, f_mid);

        if (fabs(f_mid) < 1e-12 || (a_hi2 - a_lo2) < 1e-15 * a_mid) {
            a_best = a_mid;
            printf("    Converged at iter %d: a=%.14f, f(Rmax)=%+.6e\n",
                   iter, a_mid, f_mid);
            break;
        }

        if (f_lo2 * f_mid < 0) {
            a_hi2 = a_mid;
        } else {
            a_lo2 = a_mid;
            f_lo2 = f_mid;
        }
        a_best = a_mid;
    }

    /* Final integration */
    shoot(p, a_best, c4, (double)B_val, I_val);
    res.a_shooting = a_best;
    res.en = compute_energy(p, e_skyrme, rho0, (double)B_val, I_val);
    res.prof = p;

    return res;
}

/* ========== Main ========== */

int main(int argc, char **argv)
{
    int Nr = 20000;
    double R_max = 30.0;
    double e_skyrme = 4.0;
    double rho0 = 1.0;
    int B_max = 4;

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-Nr") == 0 && i+1 < argc)
            Nr = atoi(argv[++i]);
        else if (strcmp(argv[i], "-Rmax") == 0 && i+1 < argc)
            R_max = atof(argv[++i]);
        else if (strcmp(argv[i], "-e") == 0 && i+1 < argc)
            e_skyrme = atof(argv[++i]);
        else if (strcmp(argv[i], "-rho0") == 0 && i+1 < argc)
            rho0 = atof(argv[++i]);
        else if (strcmp(argv[i], "-Bmax") == 0 && i+1 < argc)
            B_max = atoi(argv[++i]);
        else {
            fprintf(stderr,
                "Usage: %s [-Nr N] [-Rmax R] [-e E] [-rho0 R] [-Bmax B]\n",
                argv[0]);
            return 1;
        }
    }

    double c4 = 2.0 * rho0 * rho0 / (e_skyrme * e_skyrme);

    printf("============================================================\n");
    printf("  Rational Map Ansatz Skyrmion Solver (CHPT normalization)\n");
    printf("============================================================\n");
    printf("Grid: Nr=%d, dr=%.6f, R_max=%.1f\n", Nr, R_max / Nr, R_max);
    printf("Parameters: e=%.4f, rho0=%.4f, c4=%.6f\n", e_skyrme, rho0, c4);
    printf("Solving for B = 1..%d\n\n", B_max);

    /* Step 1: Compute angular integrals I for each B */
    printf("--- Computing angular integrals I ---\n");

    rational_map_fn maps[] = {
        NULL,
        rational_map_B1,
        rational_map_B2,
        rational_map_B3,
        rational_map_B4
    };

    double I_values[5] = {0.0, 0.0, 0.0, 0.0, 0.0};

    for (int B = 1; B <= B_max && B <= 4; B++) {
        I_values[B] = compute_angular_integral(maps[B], B);
    }

    printf("\nAngular integral summary:\n");
    for (int B = 1; B <= B_max && B <= 4; B++) {
        double alpha = origin_exponent((double)B);
        printf("  B=%d: I = %.6f, alpha = %.6f\n", B, I_values[B], alpha);
    }

    /* Step 2: Solve radial ODE for each B */
    SkyrmionResult results[5];
    memset(results, 0, sizeof(results));

    for (int B = 1; B <= B_max && B <= 4; B++) {
        results[B] = solve_for_B(B, I_values[B], e_skyrme, rho0, Nr, R_max);
    }

    /* Step 3: Summary table */
    printf("\n\n");
    printf("============================================================\n");
    printf("  RESULTS SUMMARY\n");
    printf("============================================================\n");
    printf("Parameters: e=%.4f, rho0=%.4f, c4=%.6f\n\n", e_skyrme, rho0, c4);

    printf("%-4s  %8s  %10s  %12s  %12s  %12s  %12s  %8s  %8s  %10s\n",
           "B", "alpha", "I", "a", "E_total", "E2", "E4",
           "Q", "E/E_FB", "E/(B*E1)");

    double E1_total = 0.0;

    for (int B = 1; B <= B_max && B <= 4; B++) {
        SkyrmionResult *res = &results[B];
        if (res->a_shooting < 0) {
            printf("%-4d  FAILED\n", B);
            continue;
        }

        double E_FB_B = 6.0 * sqrt(2.0) * M_PI * M_PI
                       * rho0 * rho0 * rho0 * (double)B / e_skyrme;
        double ratio_FB = res->en.Etotal / E_FB_B;

        if (B == 1) E1_total = res->en.Etotal;

        double binding = (E1_total > 0 && B > 0)
                       ? res->en.Etotal / ((double)B * E1_total) : 0.0;

        printf("%-4d  %8.4f  %10.4f  %12.8f  %12.6f  %12.6f  %12.6f"
               "  %8.4f  %8.4f  %10.6f\n",
               B, res->alpha, res->I_angular, res->a_shooting,
               res->en.Etotal, res->en.E2, res->en.E4,
               res->en.Q, ratio_FB, binding);
    }

    /* Step 4: Detailed output for each B */
    for (int B = 1; B <= B_max && B <= 4; B++) {
        SkyrmionResult *res = &results[B];
        if (res->a_shooting < 0) continue;

        double E_FB_B = 6.0 * sqrt(2.0) * M_PI * M_PI
                       * rho0 * rho0 * rho0 * (double)B / e_skyrme;

        printf("\n--- B=%d detailed ---\n", B);
        printf("  Origin exponent alpha = %.8f\n", res->alpha);
        printf("  Angular integral I = %.8f\n", res->I_angular);
        printf("  Shooting parameter a (coeff of r^alpha) = %.14f\n",
               res->a_shooting);
        printf("  E_total = %.12e\n", res->en.Etotal);
        printf("  E2      = %.12e\n", res->en.E2);
        printf("  E4      = %.12e\n", res->en.E4);
        printf("  Q       = %.8f (expected %d)\n", res->en.Q, B);
        printf("  E2/E4   = %.6f\n", res->en.E2 / res->en.E4);
        printf("  Virial: E2-E4 = %.6e\n", res->en.E2 - res->en.E4);
        printf("  E_FB(B=%d) = %.6e\n", B, E_FB_B);
        printf("  E/E_FB  = %.6f\n", res->en.Etotal / E_FB_B);
        if (B > 1 && E1_total > 0) {
            printf("  Binding energy per baryon: E(B)/(B*E(1)) = %.6f\n",
                   res->en.Etotal / ((double)B * E1_total));
            printf("  Binding energy: E(B) - B*E(1) = %.6f\n",
                   res->en.Etotal - (double)B * E1_total);
        }

        /* Write profile to file */
        char filename[64];
        snprintf(filename, sizeof(filename), "profile_B%d.dat", B);
        FILE *fp = fopen(filename, "w");
        if (fp) {
            fprintf(fp, "# Rational map ansatz profile for B=%d, I=%.8f\n",
                    B, res->I_angular);
            fprintf(fp, "# e=%.4f, rho0=%.4f, c4=%.6f, alpha=%.8f, "
                    "a=%.14f\n",
                    e_skyrme, rho0, c4, res->alpha, res->a_shooting);
            fprintf(fp, "# r  f(r)  f'(r)  baryon_density  "
                    "energy_density_E2  energy_density_E4\n");
            Profile *pp = res->prof;
            for (int i = 0; i <= pp->N; i++) {
                double r = pp->r[i];
                double fi = pp->f[i];
                double fpi = pp->fp[i];
                double sf = sin(fi);
                double sf2 = sf * sf;
                double bdens = (r > 1e-14) ?
                    -(double)B * fpi * sf2 / (M_PI * r * r) : 0.0;
                double ed2 = fpi * fpi * r * r + 2.0 * (double)B * sf2;
                double ed4 = (r > 1e-14) ?
                    (2.0 * (double)B * fpi * fpi * sf2
                     + res->I_angular * sf2 * sf2 / (r * r)) : 0.0;
                fprintf(fp, "%.8e  %.12e  %.12e  %.12e  %.12e  %.12e\n",
                        r, fi, fpi, bdens, ed2, ed4);
            }
            fclose(fp);
            printf("  Profile written to %s\n", filename);
        }
    }

    /* Cleanup */
    for (int B = 1; B <= B_max && B <= 4; B++) {
        if (results[B].prof) profile_free(results[B].prof);
    }

    printf("\nDone.\n");
    return 0;
}
