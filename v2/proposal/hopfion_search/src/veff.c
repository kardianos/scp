/*
 * veff.c — Effective potential for degenerate modes in modified Skyrme term
 *
 * Computes the angular-averaged coupling potential V_eff(r) that degenerate
 * perturbations δp experience near a B=1 hedgehog soliton, from the
 * modified Skyrme term |[R_μ,R_ν]|²₈ (full 8-component norm instead of
 * scalar extraction ⟨[R_μ,R_ν]²⟩₀).
 *
 * The coupling term at O(δp²) is:
 *   L₄,coupling = (1/4e²) Σ_{i<j} |[A_i, D_j] + [D_i, A_j]|²
 * where A_i = q̃₀∂_iq₀ (background right-current) and
 *       D_i = q̃₀∂_iδp + δp̃∂_iq₀ (linear perturbation current).
 *
 * We compute this for two perturbation types:
 *   1. Scalar (P mode): δp = g(r)·1    (ℓ=0)
 *   2. Vector (J mode): δp = g(r)r̂·σ  (ℓ=1, hedgehog-like)
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Quaternion: q = s + v1*i + v2*j + v3*k */
typedef struct { double s, v1, v2, v3; } Quat;

static Quat q_make(double s, double v1, double v2, double v3) {
    return (Quat){s, v1, v2, v3};
}

static Quat q_scale(Quat a, double c) {
    return q_make(c*a.s, c*a.v1, c*a.v2, c*a.v3);
}

static Quat q_add(Quat a, Quat b) {
    return q_make(a.s+b.s, a.v1+b.v1, a.v2+b.v2, a.v3+b.v3);
}

static Quat q_sub(Quat a, Quat b) {
    return q_make(a.s-b.s, a.v1-b.v1, a.v2-b.v2, a.v3-b.v3);
}

static Quat q_mul(Quat a, Quat b) {
    return q_make(
        a.s*b.s  - a.v1*b.v1 - a.v2*b.v2 - a.v3*b.v3,
        a.s*b.v1 + a.v1*b.s  + a.v2*b.v3 - a.v3*b.v2,
        a.s*b.v2 - a.v1*b.v3 + a.v2*b.s  + a.v3*b.v1,
        a.s*b.v3 + a.v1*b.v2 - a.v2*b.v1 + a.v3*b.s
    );
}

static Quat q_conj(Quat a) {
    return q_make(a.s, -a.v1, -a.v2, -a.v3);
}

static double q_norm2(Quat a) {
    return a.s*a.s + a.v1*a.v1 + a.v2*a.v2 + a.v3*a.v3;
}

/* Commutator [a,b] = ab - ba */
static Quat q_comm(Quat a, Quat b) {
    return q_sub(q_mul(a, b), q_mul(b, a));
}

/* Read profile file: columns r, f(r), [f'(r)] */
static int read_profile(const char *fname, int *npts, double **r_out,
                        double **f_out, double **fp_out)
{
    FILE *fp = fopen(fname, "r");
    if (!fp) { fprintf(stderr, "Cannot open %s\n", fname); return -1; }

    int cap = 4096;
    double *r = malloc(cap * sizeof(double));
    double *fv = malloc(cap * sizeof(double));
    double *fpv = malloc(cap * sizeof(double));
    int n = 0;
    char line[256];
    while (fgets(line, sizeof(line), fp)) {
        if (line[0] == '#') continue;
        double rv, fval, fpval = 0;
        int nc = sscanf(line, "%lf %lf %lf", &rv, &fval, &fpval);
        if (nc < 2) continue;
        if (n >= cap) {
            cap *= 2;
            r = realloc(r, cap * sizeof(double));
            fv = realloc(fv, cap * sizeof(double));
            fpv = realloc(fpv, cap * sizeof(double));
        }
        r[n] = rv; fv[n] = fval; fpv[n] = fpval;
        n++;
    }
    fclose(fp);
    *npts = n; *r_out = r; *f_out = fv; *fp_out = fpv;
    return 0;
}

/* Interpolate profile at given r */
static void interp_profile(double r, int npts, const double *rv,
                           const double *fv, const double *fpv,
                           double *f_out, double *fp_out)
{
    if (r <= rv[0]) { *f_out = fv[0]; *fp_out = fpv[0]; return; }
    if (r >= rv[npts-1]) { *f_out = fv[npts-1]; *fp_out = fpv[npts-1]; return; }

    /* Binary search */
    int lo = 0, hi = npts - 1;
    while (hi - lo > 1) {
        int mid = (lo + hi) / 2;
        if (rv[mid] <= r) lo = mid; else hi = mid;
    }

    double t = (r - rv[lo]) / (rv[hi] - rv[lo]);
    *f_out = fv[lo] + t * (fv[hi] - fv[lo]);
    *fp_out = fpv[lo] + t * (fpv[hi] - fpv[lo]);
}

/*
 * Compute the hedgehog right-current A_i at position x = r*xhat.
 * Returns A[0], A[1], A[2] (the three spatial right-currents as quaternions).
 *
 * A_i = ρ₀²[f'x̂_i(x̂·σ) + (sin(2f)/(2r))(σ_i - x̂_ix̂·σ) - (sin²f/r)(x̂×σ)_i]
 *
 * All three terms are pure quaternions (no scalar part).
 */
static void compute_A(double rho0, double f, double fp, double r,
                      const double xhat[3], Quat A[3])
{
    double sf = sin(f), cf = cos(f);
    double s2f = sin(2*f);  /* = 2*sf*cf */
    double sf2 = sf * sf;

    double inv_r = (r > 1e-10) ? 1.0/r : 0.0;

    for (int i = 0; i < 3; i++) {
        /* σ_i basis vector (pure imaginary quaternion) */
        double si[4] = {0, 0, 0, 0};
        si[i+1] = 1.0;

        /* x̂·σ as quaternion */
        Quat xdots = q_make(0, xhat[0], xhat[1], xhat[2]);

        /* Term 1: f' x̂_i (x̂·σ) */
        Quat t1 = q_scale(xdots, fp * xhat[i]);

        /* Term 2: (sin(2f)/(2r))(σ_i - x̂_i x̂·σ) */
        Quat sigma_i = q_make(si[0], si[1], si[2], si[3]);
        Quat t2 = q_scale(q_sub(sigma_i, q_scale(xdots, xhat[i])),
                          0.5 * s2f * inv_r);

        /* Term 3: -(sin²f/r)(x̂ × σ)_i
         * (x̂ × σ)_i = ε_{ijk} x̂_j σ_k (as quaternion) */
        double cross[3];
        /* Cross product x̂ × e_i (where e_i is unit vector in direction i):
         * Actually, (x̂ × σ)_i = ε_{ijk} x̂_j σ_k
         * This is a quaternion with vector part = x̂ × e_i */
        int j1 = (i+1)%3, j2 = (i+2)%3;
        cross[0] = cross[1] = cross[2] = 0;
        /* ε_{i,j1,j2} x̂_{j1} σ_{j2} - ε_{i,j2,j1} x̂_{j2} σ_{j1}
         * = x̂_{j1} σ_{j2} - x̂_{j2} σ_{j1} */
        /* Actually (x̂×σ)_i = Σ_{jk} ε_{ijk} x̂_j σ_k */
        double cv[3] = {0, 0, 0};
        for (int jj = 0; jj < 3; jj++)
            for (int kk = 0; kk < 3; kk++) {
                /* ε_{i,jj,kk} */
                int eps = 0;
                if (i == jj || i == kk || jj == kk) eps = 0;
                else if ((i==0&&jj==1&&kk==2)||(i==1&&jj==2&&kk==0)||(i==2&&jj==0&&kk==1)) eps = 1;
                else eps = -1;
                if (eps != 0) cv[kk] += eps * xhat[jj];
            }
        Quat xcross = q_make(0, cv[0], cv[1], cv[2]);
        Quat t3 = q_scale(xcross, -sf2 * inv_r);

        A[i] = q_scale(q_add(q_add(t1, t2), t3), rho0*rho0);
    }
}

/*
 * Compute D_i = q̃₀∂_i(δp) + δp̃∂_iq₀ for a scalar P-mode: δp = g(r)·1.
 *
 * ∂_i(δp) = g'(r)x̂_i · 1  (scalar quaternion)
 * δp̃ = g(r) · 1
 *
 * q̃₀ = ρ₀(cos f - sin f x̂·σ)
 * ∂_iq₀ = ρ₀[-sin f f' x̂_i + cos f f' x̂_i x̂·σ + (sin f/r)(σ_i - x̂_i x̂·σ)]
 */
static void compute_D_scalar(double rho0, double f, double fp, double r,
                             const double xhat[3],
                             double g, double gp, Quat D[3])
{
    double sf = sin(f), cf = cos(f);
    double inv_r = (r > 1e-10) ? 1.0/r : 0.0;
    Quat xdots = q_make(0, xhat[0], xhat[1], xhat[2]);

    /* q̃₀ = ρ₀(cos f - sin f x̂·σ) */
    Quat q0bar = q_add(q_make(rho0*cf, 0, 0, 0), q_scale(xdots, -rho0*sf));

    for (int i = 0; i < 3; i++) {
        /* ∂_i(δp) = g'x̂_i · 1 (scalar quaternion) */
        Quat dip = q_make(gp * xhat[i], 0, 0, 0);

        /* q̃₀ · ∂_i(δp) */
        Quat term1 = q_mul(q0bar, dip);

        /* δp̃ · ∂_iq₀ = g · ∂_iq₀ (just scaling) */
        /* ∂_iq₀ = ρ₀[-sf f' x̂_i + cf f' x̂_i x̂·σ + (sf/r)(σ_i - x̂_i x̂·σ)] */
        Quat sigma_i = q_make(0, (i==0)?1:0, (i==1)?1:0, (i==2)?1:0);
        Quat dq0_i = q_add(
            q_add(
                q_make(-rho0*sf*fp*xhat[i], 0, 0, 0),
                q_scale(xdots, rho0*cf*fp*xhat[i])
            ),
            q_scale(q_sub(sigma_i, q_scale(xdots, xhat[i])), rho0*sf*inv_r)
        );
        Quat term2 = q_scale(dq0_i, g);

        D[i] = q_add(term1, term2);
    }
}

/*
 * Compute D_i for a vector J-mode: δp = g(r)r̂·σ (hedgehog-like).
 *
 * ∂_i(δp) = g'x̂_i (x̂·σ) + (g/r)(σ_i - x̂_i x̂·σ)
 * δp̃ = -g(r)r̂·σ  (conjugate flips sign of vector part)
 */
static void compute_D_vector(double rho0, double f, double fp, double r,
                             const double xhat[3],
                             double g, double gp, Quat D[3])
{
    double sf = sin(f), cf = cos(f);
    double inv_r = (r > 1e-10) ? 1.0/r : 0.0;
    Quat xdots = q_make(0, xhat[0], xhat[1], xhat[2]);

    Quat q0bar = q_add(q_make(rho0*cf, 0, 0, 0), q_scale(xdots, -rho0*sf));

    for (int i = 0; i < 3; i++) {
        Quat sigma_i = q_make(0, (i==0)?1:0, (i==1)?1:0, (i==2)?1:0);

        /* ∂_i(δp) = g'x̂_i(x̂·σ) + (g/r)(σ_i - x̂_i x̂·σ) */
        Quat dip = q_add(
            q_scale(xdots, gp * xhat[i]),
            q_scale(q_sub(sigma_i, q_scale(xdots, xhat[i])), g * inv_r)
        );

        /* q̃₀ · ∂_i(δp) */
        Quat term1 = q_mul(q0bar, dip);

        /* δp̃ · ∂_iq₀ where δp̃ = -g(x̂·σ) (conjugate of g·x̂·σ) */
        Quat dpbar = q_scale(xdots, -g);

        /* ∂_iq₀ */
        Quat dq0_i = q_add(
            q_add(
                q_make(-rho0*sf*fp*xhat[i], 0, 0, 0),
                q_scale(xdots, rho0*cf*fp*xhat[i])
            ),
            q_scale(q_sub(sigma_i, q_scale(xdots, xhat[i])), rho0*sf*inv_r)
        );

        Quat term2 = q_mul(dpbar, dq0_i);

        D[i] = q_add(term1, term2);
    }
}

int main(int argc, char *argv[])
{
    const char *profile = "profile_B1.dat";
    double rho0 = 1.0, e = 4.0;
    int n_angular = 50;  /* Lebedev-like angular grid points per dimension */

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-profile") == 0 && i+1 < argc) profile = argv[++i];
        else if (strcmp(argv[i], "-rho0") == 0 && i+1 < argc) rho0 = atof(argv[++i]);
        else if (strcmp(argv[i], "-e") == 0 && i+1 < argc) e = atof(argv[++i]);
        else if (strcmp(argv[i], "-nang") == 0 && i+1 < argc) n_angular = atoi(argv[++i]);
    }

    /* Read profile */
    int npts;
    double *rv, *fv, *fpv;
    if (read_profile(profile, &npts, &rv, &fv, &fpv) < 0) return 1;

    /* If f'(r) not provided, compute numerically */
    int has_fp = 0;
    for (int i = 0; i < npts; i++)
        if (fabs(fpv[i]) > 1e-30) { has_fp = 1; break; }

    if (!has_fp) {
        fprintf(stderr, "Computing f'(r) numerically...\n");
        for (int i = 0; i < npts; i++) {
            if (i == 0) fpv[i] = (fv[1] - fv[0]) / (rv[1] - rv[0]);
            else if (i == npts-1) fpv[i] = (fv[npts-1] - fv[npts-2]) / (rv[npts-1] - rv[npts-2]);
            else fpv[i] = (fv[i+1] - fv[i-1]) / (rv[i+1] - rv[i-1]);
        }
    }

    double c4 = 2.0 * rho0 * rho0 / (e * e);
    double prefactor = 1.0 / (4.0 * e * e);  /* 1/(4e²) Skyrme coupling */

    printf("# Effective potential for degenerate modes near B=1 soliton\n");
    printf("# Profile: %s, rho0=%.3f, e=%.3f, c4=%.6f\n", profile, rho0, e, c4);
    printf("# Columns: r  VP_gamma  VP_alpha  VJ_gamma  VJ_alpha  E4_density\n");
    printf("# VP_gamma: P-mode mass potential (g=1,g'=0)\n");
    printf("# VP_alpha: P-mode kinetic coupling (g=0,g'=1)\n");
    printf("# VJ_gamma: J-mode mass potential (g=1,g'=0)\n");
    printf("# VJ_alpha: J-mode kinetic coupling (g=0,g'=1)\n");
    printf("# E4: standard Skyrme energy density (angular avg)\n");
    printf("# Radial eq: integral [alpha*g'^2 + gamma*g^2] r^2 dr\n");

    /* Angular quadrature: Gauss-Legendre in cos(θ), uniform in φ */
    int n_theta = n_angular, n_phi = 2 * n_angular;

    /* Precompute Gauss-Legendre nodes and weights for cos(θ) ∈ [-1,1] */
    double *cos_theta = malloc(n_theta * sizeof(double));
    double *w_theta = malloc(n_theta * sizeof(double));

    /* Simple Gauss-Legendre via Newton's method for Legendre roots */
    for (int i = 0; i < n_theta; i++) {
        /* Initial guess */
        double x = cos(M_PI * (i + 0.75) / (n_theta + 0.5));
        for (int iter = 0; iter < 20; iter++) {
            double p0 = 1, p1 = x;
            for (int j = 2; j <= n_theta; j++) {
                double p2 = ((2*j-1)*x*p1 - (j-1)*p0) / j;
                p0 = p1; p1 = p2;
            }
            double dp = n_theta * (x*p1 - p0) / (x*x - 1);
            x -= p1 / dp;
        }
        cos_theta[i] = x;
        double p0 = 1, p1 = x;
        for (int j = 2; j <= n_theta; j++) {
            double p2 = ((2*j-1)*x*p1 - (j-1)*p0) / j;
            p0 = p1; p1 = p2;
        }
        w_theta[i] = 2.0 / ((1 - x*x) * n_theta*n_theta * p0*p0);
    }

    double dphi = 2*M_PI / n_phi;

    /* Scan radial points */
    int nr_out = 500;
    double r_max = rv[npts-1];

    for (int ir = 0; ir < nr_out; ir++) {
        double r = (ir + 0.5) * r_max / nr_out;

        double f, fp;
        interp_profile(r, npts, rv, fv, fpv, &f, &fp);

        /* Compute three coupling coefficients separately:
         * γ(r): from (g=1, g'=0) — "mass-like" potential
         * α(r): from (g=0, g'=1) — kinetic modification
         * These enter the radial equation as: α(r)g'' + ... + γ(r)g */

        double sum_P_gamma = 0.0, sum_P_alpha = 0.0;
        double sum_J_gamma = 0.0, sum_J_alpha = 0.0;
        double sum_E4 = 0.0;

        for (int it = 0; it < n_theta; it++) {
            double ct = cos_theta[it];
            double st = sqrt(1.0 - ct*ct);

            for (int ip = 0; ip < n_phi; ip++) {
                double phi = ip * dphi;
                double xhat[3] = {st*cos(phi), st*sin(phi), ct};

                double wt = w_theta[it] * dphi / (4*M_PI);  /* normalize to 1 */

                /* Background right-currents */
                Quat A[3];
                compute_A(rho0, f, fp, r, xhat, A);

                /* Standard Skyrme: Σ_{i<j} |[A_i,A_j]|² */
                for (int ii = 0; ii < 3; ii++)
                    for (int jj = ii+1; jj < 3; jj++) {
                        Quat comm = q_comm(A[ii], A[jj]);
                        sum_E4 += wt * q_norm2(comm);
                    }

                /* Scalar P-mode: γ (g=1,g'=0) and α (g=0,g'=1) */
                Quat Dp_g[3], Dp_gp[3];
                compute_D_scalar(rho0, f, fp, r, xhat, 1.0, 0.0, Dp_g);
                compute_D_scalar(rho0, f, fp, r, xhat, 0.0, 1.0, Dp_gp);
                for (int ii = 0; ii < 3; ii++)
                    for (int jj = ii+1; jj < 3; jj++) {
                        Quat c1g = q_add(q_comm(A[ii], Dp_g[jj]),
                                         q_comm(Dp_g[ii], A[jj]));
                        sum_P_gamma += wt * q_norm2(c1g);
                        Quat c1gp = q_add(q_comm(A[ii], Dp_gp[jj]),
                                          q_comm(Dp_gp[ii], A[jj]));
                        sum_P_alpha += wt * q_norm2(c1gp);
                    }

                /* Vector J-mode: γ and α */
                Quat Dj_g[3], Dj_gp[3];
                compute_D_vector(rho0, f, fp, r, xhat, 1.0, 0.0, Dj_g);
                compute_D_vector(rho0, f, fp, r, xhat, 0.0, 1.0, Dj_gp);
                for (int ii = 0; ii < 3; ii++)
                    for (int jj = ii+1; jj < 3; jj++) {
                        Quat c1g = q_add(q_comm(A[ii], Dj_g[jj]),
                                         q_comm(Dj_g[ii], A[jj]));
                        sum_J_gamma += wt * q_norm2(c1g);
                        Quat c1gp = q_add(q_comm(A[ii], Dj_gp[jj]),
                                          q_comm(Dj_gp[ii], A[jj]));
                        sum_J_alpha += wt * q_norm2(c1gp);
                    }
            }
        }

        /* Multiply by prefactor (1/4e²) */
        double VP_gamma = prefactor * sum_P_gamma;  /* P-mode mass potential */
        double VP_alpha = prefactor * sum_P_alpha;  /* P-mode kinetic mod */
        double VJ_gamma = prefactor * sum_J_gamma;
        double VJ_alpha = prefactor * sum_J_alpha;
        double E4_dens = prefactor * sum_E4;

        printf("%10.6f  %14.8e  %14.8e  %14.8e  %14.8e  %14.8e\n",
               r, VP_gamma, VP_alpha, VJ_gamma, VJ_alpha, E4_dens);
    }

    free(rv); free(fv); free(fpv);
    free(cos_theta); free(w_theta);
    return 0;
}
