/*
 * bound.c — Bound state eigenvalue solver for degenerate perturbations
 *
 * Computes discrete bound state spectra for degenerate (p-sector)
 * perturbations around a B=1 hedgehog Skyrmion background.
 *
 * Two channels:
 *   P-mode (ℓ=0): δp = g(r)·1      (pseudoscalar)
 *   J-mode (ℓ=1): δp = g(r)·r̂·σ    (vector, hedgehog-like)
 *
 * Energy functional for radial mode g(r):
 *   E[g] = 4π ∫ r² [K(r) g'² + β(r) g g' + W(r) g²] dr
 *   T[g] = (ω²/2c²) · 4π ∫ r² g² dr
 *
 * Two coupling modes:
 *
 * MODE 1 (Skyrme, -skyrme): Three sources contribute to K, W, β:
 *   1. E_{2,D} = (1/2)|∇p|²           — free kinetic
 *   2. E_int  = (g_c²/2)ρ(r)²|∇p|²   — attractive (finite-λ)
 *   3. E_{4,C} = (1/4e²)Σ|F^w_{ij}|²  — Skyrme cross-coupling (always repulsive)
 *   4. V_D    = (μ²/2)|p|²            — degenerate mass
 *
 * MODE 2 (Geometric, -geom): Covariant derivative coupling:
 *   E_{cov} = (1/2)|D_i p|²  where D_i p = ∂_i p + g_cov [A_i, p]
 *   A_i = q̂₀⁻¹ ∂_i q̂₀ is the SU(2) connection from the hedgehog.
 *   J-mode (adjoint): full covariant coupling, cross-term can be NEGATIVE (attractive)
 *   P-mode (singlet): [A_i, scalar] = 0, no coupling
 *   Plus E_int and V_D as before.
 *
 * Eigenvalue equation (from δ(T-E)/δg = 0):
 *   (ω²/c²) g = -(1/r²) d/dr[r²·2K(r)·g'] + [2W(r) - β'(r) - 2β(r)/r]·g
 *
 * Bound state: ω² < c²μ², g → exp(-κr)/r at r→∞ where κ² = μ² - ω²/c².
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ========== Quaternion algebra (from veff.c) ========== */

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

static double q_norm2(Quat a) {
    return a.s*a.s + a.v1*a.v1 + a.v2*a.v2 + a.v3*a.v3;
}

static Quat q_comm(Quat a, Quat b) {
    return q_sub(q_mul(a, b), q_mul(b, a));
}

/* ========== Profile I/O ========== */

typedef struct {
    int n;
    double *r, *f, *fp, *rho;
} Profile;

/* Read profile: supports 3-col (r,f,f'), 4-col (r,f,f',ρ) and
 * 6-col sigma-model (r,f,f',baryon_dens,E2_dens,E4_dens) formats.
 * Only exactly 4 columns is treated as having ρ(r) data. */
static int read_profile(const char *fname, Profile *prof) {
    FILE *fp = fopen(fname, "r");
    if (!fp) { fprintf(stderr, "Cannot open %s\n", fname); return -1; }

    int cap = 4096;
    prof->r   = malloc(cap * sizeof(double));
    prof->f   = malloc(cap * sizeof(double));
    prof->fp  = malloc(cap * sizeof(double));
    prof->rho = malloc(cap * sizeof(double));
    int n = 0;
    char line[512];
    int has_rho = 0;
    int ncols_detected = 0;

    while (fgets(line, sizeof(line), fp)) {
        if (line[0] == '#') continue;
        double rv, fval, fpval = 0, rhoval = 1.0;
        double dummy1, dummy2;
        int nc = sscanf(line, "%lf %lf %lf %lf %lf %lf",
                        &rv, &fval, &fpval, &rhoval, &dummy1, &dummy2);
        if (nc < 2) continue;
        if (ncols_detected == 0) ncols_detected = nc;
        /* Only 4-column format has ρ as 4th col (finite-lambda profiles).
         * 6-column format has baryon_dens as 4th col — use ρ=1. */
        if (nc >= 5) rhoval = 1.0;
        if (nc == 4) has_rho = 1;
        if (n >= cap) {
            cap *= 2;
            prof->r   = realloc(prof->r,   cap * sizeof(double));
            prof->f   = realloc(prof->f,   cap * sizeof(double));
            prof->fp  = realloc(prof->fp,  cap * sizeof(double));
            prof->rho = realloc(prof->rho, cap * sizeof(double));
        }
        prof->r[n] = rv; prof->f[n] = fval;
        prof->fp[n] = fpval; prof->rho[n] = rhoval;
        n++;
    }
    fclose(fp);
    prof->n = n;

    /* Check if f' was provided (may be 0 in 2-col format) */
    int has_fp = 0;
    for (int i = 0; i < n; i++)
        if (fabs(prof->fp[i]) > 1e-30) { has_fp = 1; break; }

    if (!has_fp) {
        fprintf(stderr, "Computing f'(r) numerically...\n");
        for (int i = 0; i < n; i++) {
            if (i == 0) prof->fp[i] = (prof->f[1] - prof->f[0]) / (prof->r[1] - prof->r[0]);
            else if (i == n-1) prof->fp[i] = (prof->f[n-1] - prof->f[n-2]) / (prof->r[n-1] - prof->r[n-2]);
            else prof->fp[i] = (prof->f[i+1] - prof->f[i-1]) / (prof->r[i+1] - prof->r[i-1]);
        }
    }

    if (!has_rho) {
        if (ncols_detected >= 5)
            fprintf(stderr, "Detected %d-col sigma-model format — using rho=1\n", ncols_detected);
        else
            fprintf(stderr, "No rho(r) column — using rho=1 (sigma model)\n");
    }

    return 0;
}

static void free_profile(Profile *p) {
    free(p->r); free(p->f); free(p->fp); free(p->rho);
}

/* Interpolate profile at given r (linear) */
static void interp_profile(const Profile *prof, double r,
                           double *f_out, double *fp_out, double *rho_out)
{
    if (r <= prof->r[0]) {
        *f_out = prof->f[0]; *fp_out = prof->fp[0]; *rho_out = prof->rho[0];
        return;
    }
    if (r >= prof->r[prof->n-1]) {
        *f_out = prof->f[prof->n-1]; *fp_out = prof->fp[prof->n-1];
        *rho_out = prof->rho[prof->n-1];
        return;
    }
    int lo = 0, hi = prof->n - 1;
    while (hi - lo > 1) { int mid = (lo+hi)/2; if (prof->r[mid] <= r) lo = mid; else hi = mid; }
    double t = (r - prof->r[lo]) / (prof->r[hi] - prof->r[lo]);
    *f_out   = prof->f[lo]   + t * (prof->f[hi]   - prof->f[lo]);
    *fp_out  = prof->fp[lo]  + t * (prof->fp[hi]  - prof->fp[lo]);
    *rho_out = prof->rho[lo] + t * (prof->rho[hi] - prof->rho[lo]);
}

/* ========== Angular integration (from veff.c) ========== */

/* Compute hedgehog right-current A_i at position r*x̂ with local ρ */
static void compute_A(double rho, double f, double fp, double r,
                      const double xhat[3], Quat A[3])
{
    double sf = sin(f);
    double s2f = sin(2*f), sf2 = sf*sf;
    double rho2 = rho*rho;
    double inv_r = (r > 1e-10) ? 1.0/r : 0.0;

    for (int i = 0; i < 3; i++) {
        double si[4] = {0,0,0,0}; si[i+1] = 1.0;
        Quat xdots = q_make(0, xhat[0], xhat[1], xhat[2]);
        Quat t1 = q_scale(xdots, fp * xhat[i]);
        Quat sigma_i = q_make(si[0], si[1], si[2], si[3]);
        Quat t2 = q_scale(q_sub(sigma_i, q_scale(xdots, xhat[i])), 0.5*s2f*inv_r);

        int j1 = (i+1)%3, j2 = (i+2)%3;
        double cv[3] = {0,0,0};
        for (int jj = 0; jj < 3; jj++)
            for (int kk = 0; kk < 3; kk++) {
                int eps = 0;
                if (i==jj||i==kk||jj==kk) eps=0;
                else if ((i==0&&jj==1&&kk==2)||(i==1&&jj==2&&kk==0)||(i==2&&jj==0&&kk==1)) eps=1;
                else eps=-1;
                if (eps != 0) cv[kk] += eps * xhat[jj];
            }
        (void)j1; (void)j2;
        Quat t3 = q_scale(q_make(0, cv[0], cv[1], cv[2]), -sf2*inv_r);
        A[i] = q_scale(q_add(q_add(t1, t2), t3), rho2);
    }
}

/* D_i for scalar P-mode: δp = g·1 */
static void compute_D_scalar(double rho, double f, double fp, double r,
                             const double xhat[3], double g, double gp, Quat D[3])
{
    double sf = sin(f), cf = cos(f);
    double inv_r = (r > 1e-10) ? 1.0/r : 0.0;
    Quat xdots = q_make(0, xhat[0], xhat[1], xhat[2]);
    Quat q0bar = q_add(q_make(rho*cf, 0, 0, 0), q_scale(xdots, -rho*sf));

    for (int i = 0; i < 3; i++) {
        Quat dip = q_make(gp * xhat[i], 0, 0, 0);
        Quat term1 = q_mul(q0bar, dip);

        Quat sigma_i = q_make(0, (i==0)?1:0, (i==1)?1:0, (i==2)?1:0);
        Quat dq0_i = q_add(q_add(
            q_make(-rho*sf*fp*xhat[i], 0, 0, 0),
            q_scale(xdots, rho*cf*fp*xhat[i])),
            q_scale(q_sub(sigma_i, q_scale(xdots, xhat[i])), rho*sf*inv_r));
        Quat term2 = q_scale(dq0_i, g);
        D[i] = q_add(term1, term2);
    }
}

/* D_i for vector J-mode: δp = g·r̂·σ */
static void compute_D_vector(double rho, double f, double fp, double r,
                             const double xhat[3], double g, double gp, Quat D[3])
{
    double sf = sin(f), cf = cos(f);
    double inv_r = (r > 1e-10) ? 1.0/r : 0.0;
    Quat xdots = q_make(0, xhat[0], xhat[1], xhat[2]);
    Quat q0bar = q_add(q_make(rho*cf, 0, 0, 0), q_scale(xdots, -rho*sf));

    for (int i = 0; i < 3; i++) {
        Quat sigma_i = q_make(0, (i==0)?1:0, (i==1)?1:0, (i==2)?1:0);
        Quat dip = q_add(q_scale(xdots, gp*xhat[i]),
                         q_scale(q_sub(sigma_i, q_scale(xdots, xhat[i])), g*inv_r));
        Quat term1 = q_mul(q0bar, dip);

        Quat dpbar = q_scale(xdots, -g);
        Quat dq0_i = q_add(q_add(
            q_make(-rho*sf*fp*xhat[i], 0, 0, 0),
            q_scale(xdots, rho*cf*fp*xhat[i])),
            q_scale(q_sub(sigma_i, q_scale(xdots, xhat[i])), rho*sf*inv_r));
        Quat term2 = q_mul(dpbar, dq0_i);
        D[i] = q_add(term1, term2);
    }
}

/* Compute E_{4,C} contribution at one radial point via angular average.
 * Returns alpha (g'² coefficient), gamma (g² coefficient),
 * and beta (g·g' cross-coefficient) for the specified mode. */
static void compute_skyrme_coupling(double rho_val, double f, double fp, double r,
                                    double e_skyrme, int n_theta, int n_phi,
                                    const double *cos_theta, const double *w_theta,
                                    double dphi, int mode, /* 0=P, 1=J */
                                    double *alpha_out, double *gamma_out, double *beta_out)
{
    double prefactor = 1.0 / (4.0 * e_skyrme * e_skyrme);
    double sum_alpha = 0, sum_gamma = 0, sum_mixed = 0;

    for (int it = 0; it < n_theta; it++) {
        double ct = cos_theta[it];
        double st = sqrt(1.0 - ct*ct);
        for (int ip = 0; ip < n_phi; ip++) {
            double phi = ip * dphi;
            double xhat[3] = {st*cos(phi), st*sin(phi), ct};
            double wt = w_theta[it] * dphi / (4*M_PI);

            Quat A[3];
            compute_A(rho_val, f, fp, r, xhat, A);

            /* Compute |F^w|² for (g=0,g'=1), (g=1,g'=0), and (g=1,g'=1) */
            Quat D_gp[3], D_g[3], D_mix[3];
            if (mode == 0) {
                compute_D_scalar(rho_val, f, fp, r, xhat, 0.0, 1.0, D_gp);
                compute_D_scalar(rho_val, f, fp, r, xhat, 1.0, 0.0, D_g);
                compute_D_scalar(rho_val, f, fp, r, xhat, 1.0, 1.0, D_mix);
            } else {
                compute_D_vector(rho_val, f, fp, r, xhat, 0.0, 1.0, D_gp);
                compute_D_vector(rho_val, f, fp, r, xhat, 1.0, 0.0, D_g);
                compute_D_vector(rho_val, f, fp, r, xhat, 1.0, 1.0, D_mix);
            }

            for (int ii = 0; ii < 3; ii++)
                for (int jj = ii+1; jj < 3; jj++) {
                    Quat c_gp  = q_add(q_comm(A[ii], D_gp[jj]),  q_comm(D_gp[ii],  A[jj]));
                    Quat c_g   = q_add(q_comm(A[ii], D_g[jj]),   q_comm(D_g[ii],   A[jj]));
                    Quat c_mix = q_add(q_comm(A[ii], D_mix[jj]), q_comm(D_mix[ii], A[jj]));
                    sum_alpha += wt * q_norm2(c_gp);
                    sum_gamma += wt * q_norm2(c_g);
                    sum_mixed += wt * q_norm2(c_mix);
                }
        }
    }

    *alpha_out = prefactor * sum_alpha;
    *gamma_out = prefactor * sum_gamma;
    /* β from bilinearity: |F(g=1,g'=1)|² = α + γ + 2β(for g·g' cross-term) */
    *beta_out  = prefactor * (sum_mixed - sum_alpha - sum_gamma) / 2.0;
}

/* ========== Covariant derivative coupling (geometric mode) ========== */

/*
 * Compute covariant derivative kinetic term (1/2)|D_i p|² decomposed into
 * g'², g·g', g² coefficients, averaged over angles.
 *
 * D_i p = ∂_i p + g_cov [A_i, p]
 * A_i = q̂₀⁻¹ ∂_i q̂₀ = f'x̂_i(x̂·σ) + sin(2f)/(2r)(σ_i - x̂_i(x̂·σ)) - sin²f/r (x̂×ê_i)·σ
 *
 * P-mode (singlet): [A_i, scalar] = 0, so D_i = ∂_i (free)
 * J-mode (adjoint): [A_i, x̂·σ] = 2(A_i^vec × x̂)·σ (non-trivial coupling)
 *
 * Returns K, W, β INCLUDING the free kinetic term (no need to add 1/2 or 1/r² separately).
 */
static void compute_cov_coupling(double f, double fp, double r,
                                 double g_cov_param,
                                 int n_theta, int n_phi,
                                 const double *cos_theta, const double *w_theta,
                                 double dphi, int mode,
                                 double *K_out, double *W_out, double *beta_out)
{
    if (mode == 0) {
        /* P-mode: singlet under SU(2), no covariant coupling */
        *K_out = 0.5;   /* (1/2)|∂_i(g·1)|² = (1/2)g'² */
        *W_out = 0.0;
        *beta_out = 0.0;
        return;
    }

    /* J-mode: adjoint representation, full covariant coupling */
    double inv_r = (r > 1e-10) ? 1.0/r : 0.0;
    double sf = sin(f), cf = cos(f);
    double s2f = 2*sf*cf;  /* sin(2f) */
    double sf2 = sf*sf;    /* sin²(f) */

    double sum_K = 0, sum_W = 0, sum_beta = 0;

    for (int it = 0; it < n_theta; it++) {
        double ct = cos_theta[it];
        double st = sqrt(1.0 - ct*ct);
        for (int ip = 0; ip < n_phi; ip++) {
            double phi = ip * dphi;
            double xh[3] = {st*cos(phi), st*sin(phi), ct};
            double wt = w_theta[it] * dphi / (4*M_PI);

            double K_pt = 0, W_pt = 0, beta_pt = 0;

            for (int i = 0; i < 3; i++) {
                double ei[3] = {0,0,0};
                ei[i] = 1.0;

                /* V_i^{(1)} = x̂_i(x̂·σ)  — coefficient of g' in ∂_i p */
                Quat V1 = q_make(0, xh[0]*xh[i], xh[1]*xh[i], xh[2]*xh[i]);

                /* V_i^{(2)} = (1/r)(σ_i - x̂_i(x̂·σ))  — coefficient of g in ∂_i p */
                Quat V2 = q_make(0,
                    (ei[0] - xh[0]*xh[i]) * inv_r,
                    (ei[1] - xh[1]*xh[i]) * inv_r,
                    (ei[2] - xh[2]*xh[i]) * inv_r);

                /* Connection A_i (pure imaginary quaternion):
                 * A_i = f'x̂_i(x̂·σ) + s2f/(2r)(σ_i - x̂_i(x̂·σ)) - sf²/r(x̂×ê_i)·σ */
                double A_v[3];
                double cross_xe[3] = {
                    xh[1]*ei[2] - xh[2]*ei[1],
                    xh[2]*ei[0] - xh[0]*ei[2],
                    xh[0]*ei[1] - xh[1]*ei[0]
                };
                for (int k = 0; k < 3; k++) {
                    A_v[k] = fp * xh[i] * xh[k]                    /* radial */
                           + s2f*0.5*inv_r * (ei[k] - xh[k]*xh[i]) /* angular */
                           - sf2*inv_r * cross_xe[k];               /* cross */
                }

                /* C_i = [A_i, x̂·σ] = 2(A_vec × x̂)·σ */
                double C_v[3] = {
                    2*(A_v[1]*xh[2] - A_v[2]*xh[1]),
                    2*(A_v[2]*xh[0] - A_v[0]*xh[2]),
                    2*(A_v[0]*xh[1] - A_v[1]*xh[0])
                };

                /* W2_i = V_i^{(2)} + g_cov * C_i  (coefficient of g in D_i p) */
                double W2_v[3];
                for (int k = 0; k < 3; k++)
                    W2_v[k] = V2.v1*(k==0) + V2.v2*(k==1) + V2.v3*(k==2)
                              + g_cov_param * C_v[k];

                /* Accumulate |V1|², 2⟨V1,W2⟩, |W2|² */
                double V1_v[3] = {V1.v1, V1.v2, V1.v3};
                K_pt    += V1_v[0]*V1_v[0] + V1_v[1]*V1_v[1] + V1_v[2]*V1_v[2];
                beta_pt += 2*(V1_v[0]*W2_v[0] + V1_v[1]*W2_v[1] + V1_v[2]*W2_v[2]);
                W_pt    += W2_v[0]*W2_v[0] + W2_v[1]*W2_v[1] + W2_v[2]*W2_v[2];
            }

            sum_K    += wt * K_pt;
            sum_beta += wt * beta_pt;
            sum_W    += wt * W_pt;
        }
    }

    /* Factor of 1/2 from energy functional (1/2)|D_i p|² */
    *K_out    = sum_K / 2.0;
    *W_out    = sum_W / 2.0;
    *beta_out = sum_beta / 2.0;
}

/* ========== Gauss-Legendre nodes/weights ========== */

static void gauss_legendre(int n, double *nodes, double *weights) {
    for (int i = 0; i < n; i++) {
        double x = cos(M_PI * (i + 0.75) / (n + 0.5));
        for (int iter = 0; iter < 20; iter++) {
            double p0 = 1, p1 = x;
            for (int j = 2; j <= n; j++) {
                double p2 = ((2*j-1)*x*p1 - (j-1)*p0) / j;
                p0 = p1; p1 = p2;
            }
            double dp = n * (x*p1 - p0) / (x*x - 1);
            x -= p1 / dp;
        }
        nodes[i] = x;
        double p0 = 1, p1 = x;
        for (int j = 2; j <= n; j++) {
            double p2 = ((2*j-1)*x*p1 - (j-1)*p0) / j;
            p0 = p1; p1 = p2;
        }
        /* w_i = 2(1-x_i²) / (n² [P_{n-1}(x_i)]²) */
        weights[i] = 2.0 * (1 - x*x) / (n*n * p0*p0);
    }
}

/* ========== Potential construction ========== */

typedef struct {
    int nr;
    double dr, r_max;
    double *r, *K, *W, *beta;
    double *Weff;  /* W_eff = W - β'/(2) - β/r (for the ODE) */
} Potential;

/* coupling_type: 0 = Skyrme (E_{4,C}), 1 = Geometric (covariant derivative) */
static Potential *build_potential(const Profile *prof, double e_skyrme,
                                 double gc2, double mu2, int mode,
                                 int nr, double r_max, int n_ang,
                                 int coupling_type, double g_cov)
{
    Potential *pot = calloc(1, sizeof(Potential));
    pot->nr = nr;
    pot->dr = r_max / nr;
    pot->r_max = r_max;
    pot->r    = calloc(nr+1, sizeof(double));
    pot->K    = calloc(nr+1, sizeof(double));
    pot->W    = calloc(nr+1, sizeof(double));
    pot->beta = calloc(nr+1, sizeof(double));
    pot->Weff = calloc(nr+1, sizeof(double));

    /* Gauss-Legendre quadrature */
    int n_theta = n_ang, n_phi = 2*n_ang;
    double *cos_theta = malloc(n_theta * sizeof(double));
    double *w_theta   = malloc(n_theta * sizeof(double));
    gauss_legendre(n_theta, cos_theta, w_theta);
    double dphi = 2*M_PI / n_phi;

    for (int ir = 0; ir <= nr; ir++) {
        double r = (ir == 0) ? 1e-6 : ir * pot->dr;  /* avoid r=0 singularity */
        pot->r[ir] = r;

        double f, fp, rho;
        interp_profile(prof, r, &f, &fp, &rho);
        double rho2 = rho * rho;

        if (coupling_type == 1) {
            /* Geometric (covariant derivative) coupling */
            double K_cov, W_cov, beta_cov;
            compute_cov_coupling(f, fp, r, g_cov,
                                n_theta, n_phi, cos_theta, w_theta, dphi,
                                mode, &K_cov, &W_cov, &beta_cov);

            /* K_cov, W_cov already include free kinetic term */
            if (mode == 0) {
                pot->K[ir] = K_cov + gc2*rho2/2.0;
                pot->W[ir] = W_cov + mu2/2.0;
            } else {
                double inv_r2 = 1.0/(r*r);
                pot->K[ir] = K_cov + gc2*rho2/2.0;
                pot->W[ir] = W_cov + gc2*rho2*inv_r2 + mu2/2.0;
            }
            pot->beta[ir] = beta_cov;
        } else {
            /* Skyrme (E_{4,C}) coupling */
            double alpha, gamma_val, beta_val;
            compute_skyrme_coupling(rho, f, fp, r, e_skyrme,
                                   n_theta, n_phi, cos_theta, w_theta, dphi,
                                   mode, &alpha, &gamma_val, &beta_val);

            if (mode == 0) {
                pot->K[ir] = (1.0 + gc2*rho2)/2.0 + alpha;
                pot->W[ir] = gamma_val + mu2/2.0;
            } else {
                double inv_r2 = 1.0/(r*r);
                pot->K[ir] = (1.0 + gc2*rho2)/2.0 + alpha;
                pot->W[ir] = (1.0 + gc2*rho2)*inv_r2 + gamma_val + mu2/2.0;
            }
            pot->beta[ir] = beta_val;
        }
    }

    /* Compute V_eff for the momentum formulation:
     * V(r) = -β²/(2K) + 2W is the effective potential in the (g,p) system.
     * Binding requires V(r)/2 < μ²/2 (i.e., V/2 < threshold).
     * Store V/2 in Weff array for compatibility with diagnostics. */
    for (int ir = 0; ir <= nr; ir++) {
        double K = pot->K[ir];
        if (K < 1e-15) K = 1e-15;
        pot->Weff[ir] = pot->W[ir] - pot->beta[ir]*pot->beta[ir]/(4.0*K);
    }

    free(cos_theta); free(w_theta);
    return pot;
}

static void free_potential(Potential *p) {
    free(p->r); free(p->K); free(p->W); free(p->beta); free(p->Weff); free(p);
}

/* Interpolate potential arrays at arbitrary r (K, W, beta — no derivatives needed) */
static void interp_potential(const Potential *pot, double r,
                             double *K_out, double *W_out, double *beta_out)
{
    double idx = r / pot->dr;
    int i0 = (int)idx;
    if (i0 < 0) i0 = 0;
    if (i0 >= pot->nr) i0 = pot->nr - 1;
    int i1 = (i0 < pot->nr) ? i0+1 : i0;
    double t = idx - i0;
    if (t < 0) t = 0;
    if (t > 1) t = 1;

    *K_out    = (1-t)*pot->K[i0]    + t*pot->K[i1];
    *W_out    = (1-t)*pot->W[i0]    + t*pot->W[i1];
    *beta_out = (1-t)*pot->beta[i0] + t*pot->beta[i1];
}

/* ========== Shooting solver ========== */

/*
 * Momentum formulation to avoid noisy β' derivative.
 *
 * Energy: E[g] = 4π ∫ r² [K g'² + β g g' + W g²] dr
 *
 * Define conjugate momentum: p = r²(2K g' + β g)
 * Then the Euler-Lagrange system is:
 *   g' = p/(2K r²) - β g/(2K)
 *   p' = β p/(2K) + r²(-β²/(2K) + 2W - ω²/c²) g
 *
 * This uses only K, W, β (all smooth) — no derivatives needed.
 */

typedef struct {
    double omega2;       /* trial eigenvalue */
    double c2;           /* speed of light squared (=1) */
    const Potential *pot;
} ShootParams;

/* RHS for (g, p) system */
static void gp_rhs(double r, double g, double p,
                    const ShootParams *sp,
                    double *dg, double *dp)
{
    double K, W, beta;
    interp_potential(sp->pot, r, &K, &W, &beta);
    if (K < 1e-15) K = 1e-15;
    double r2 = r * r;
    if (r2 < 1e-20) r2 = 1e-20;

    *dg = p / (2.0*K*r2) - beta*g / (2.0*K);
    *dp = beta*p / (2.0*K) + r2 * (-beta*beta/(2.0*K) + 2.0*W - sp->omega2/sp->c2) * g;
}

/* Shoot from r=delta to r=R_max; returns g(R_max) and optionally fills wavefunction */
static double shoot(const ShootParams *sp, int mode,
                    double *g_out, double *gp_out, int nr_out)
{
    const Potential *pot = sp->pot;
    double dr = pot->dr;
    int nr = pot->nr;
    double delta = pot->r[0];

    /* Initial conditions */
    double g, p;
    if (mode == 0) {
        /* P-mode (ℓ=0): g(0)=1, g'(0)=0 */
        g = 1.0;
        /* p = r²(2K·g' + β·g) = δ²(0 + β(δ)·1) ≈ 0 */
        p = delta*delta * pot->beta[0];
    } else {
        /* J-mode (ℓ=1): g(δ)=δ, g'(δ)=1 */
        g = delta;
        /* p = δ²(2K(δ)·1 + β(δ)·δ) */
        p = delta*delta * (2.0*pot->K[0] + pot->beta[0]*delta);
    }

    if (g_out && nr_out > 0) {
        g_out[0] = g;
        /* Store g' for output: g' = p/(2Kr²) - βg/(2K) */
        double K0 = pot->K[0]; if (K0 < 1e-15) K0 = 1e-15;
        double d2 = delta*delta; if (d2 < 1e-20) d2 = 1e-20;
        gp_out[0] = p/(2.0*K0*d2) - pot->beta[0]*g/(2.0*K0);
    }

    /* RK4 integration of (g, p) system */
    for (int ir = 0; ir < nr; ir++) {
        double r = pot->r[ir];
        if (r < 1e-10) r = 1e-10;

        double k1g, k1p, k2g, k2p, k3g, k3p, k4g, k4p;
        double rmid = r + 0.5*dr, rend = r + dr;

        gp_rhs(r, g, p, sp, &k1g, &k1p);
        gp_rhs(rmid, g+0.5*dr*k1g, p+0.5*dr*k1p, sp, &k2g, &k2p);
        gp_rhs(rmid, g+0.5*dr*k2g, p+0.5*dr*k2p, sp, &k3g, &k3p);
        gp_rhs(rend, g+dr*k3g, p+dr*k3p, sp, &k4g, &k4p);

        g += (dr/6.0)*(k1g + 2*k2g + 2*k3g + k4g);
        p += (dr/6.0)*(k1p + 2*k2p + 2*k3p + k4p);

        if (g_out && ir+1 < nr_out) {
            g_out[ir+1] = g;
            /* Recover g' from momentum for output */
            double Ki, Wi, bi;
            interp_potential(pot, rend, &Ki, &Wi, &bi);
            if (Ki < 1e-15) Ki = 1e-15;
            double rr = rend*rend; if (rr < 1e-20) rr = 1e-20;
            gp_out[ir+1] = p/(2.0*Ki*rr) - bi*g/(2.0*Ki);
        }

        /* Overflow protection */
        if (fabs(g) > 1e30 || fabs(p) > 1e40) {
            double scale = 1e-20;
            g *= scale; p *= scale;
            if (g_out) {
                for (int j = 0; j <= ir+1 && j < nr_out; j++) {
                    g_out[j] *= scale;
                    gp_out[j] *= scale;
                }
            }
        }
    }

    return g;
}

/* Count nodes (sign changes) in wavefunction */
static int count_nodes(const double *g, int n) {
    int nodes = 0;
    for (int i = 1; i < n; i++)
        if (g[i-1]*g[i] < 0) nodes++;
    return nodes;
}

/* Find all eigenvalues in [omega2_lo, omega2_hi] by scanning for sign changes
 * of g(R_max), then refine each by bisection.
 * Returns number of eigenvalues found. eigenvalues[] must have room for max_ev entries. */
static int find_all_eigenvalues(const Potential *pot, int mode,
                                double omega2_lo, double omega2_hi,
                                double tol, double c2,
                                double *eigenvalues, int max_ev)
{
    ShootParams sp = { .omega2 = 0, .c2 = c2, .pot = pot };
    int nr = pot->nr;
    double *gw = malloc((nr+1) * sizeof(double));
    double *gpw = malloc((nr+1) * sizeof(double));

    /* Scan ω² range to find brackets (sign changes of g(R_max)) */
    int n_scan = 2000;
    double d_omega2 = (omega2_hi - omega2_lo) / n_scan;

    /* Store bracket endpoints */
    double brackets_lo[100], brackets_hi[100];
    int n_brackets = 0;

    sp.omega2 = omega2_lo;
    double g_prev = shoot(&sp, mode, NULL, NULL, 0);

    for (int is = 1; is <= n_scan; is++) {
        double om2 = omega2_lo + is * d_omega2;
        sp.omega2 = om2;
        double g_cur = shoot(&sp, mode, NULL, NULL, 0);

        if (g_prev * g_cur < 0 && n_brackets < 100) {
            brackets_lo[n_brackets] = om2 - d_omega2;
            brackets_hi[n_brackets] = om2;
            n_brackets++;
        }
        g_prev = g_cur;
    }

    /* Refine each bracket by bisection */
    int n_ev = 0;
    for (int ib = 0; ib < n_brackets && n_ev < max_ev; ib++) {
        double lo = brackets_lo[ib], hi = brackets_hi[ib];

        sp.omega2 = lo;
        double g_lo = shoot(&sp, mode, NULL, NULL, 0);

        for (int bi = 0; bi < 200; bi++) {
            double mid = 0.5*(lo + hi);
            sp.omega2 = mid;
            double g_mid = shoot(&sp, mode, NULL, NULL, 0);

            if (g_mid * g_lo > 0) { lo = mid; g_lo = g_mid; }
            else hi = mid;

            if (hi - lo < tol) break;
        }
        eigenvalues[n_ev++] = 0.5*(lo + hi);
    }

    free(gw); free(gpw);
    return n_ev;
}

/* ========== Variational energy check ========== */

/* Compute E[g]/T[g] = ω² from the wavefunction by direct quadrature */
static double variational_omega2(const Potential *pot, const double *g,
                                 const double *gp, int n, double c2)
{
    double E_num = 0, T_num = 0;
    double dr = pot->dr;

    for (int i = 0; i <= n && i <= pot->nr; i++) {
        double r = pot->r[i];
        double w = dr * ((i == 0 || i == pot->nr) ? 0.5 : 1.0);
        double r2 = r * r;

        T_num += w * r2 * g[i] * g[i];
        E_num += w * r2 * (pot->K[i]*gp[i]*gp[i] + pot->beta[i]*g[i]*gp[i]
                           + pot->W[i]*g[i]*g[i]);
    }

    /* E[g] = 4π ∫ r²[...]; T = (ω²/2c²) 4π ∫ r² g² */
    /* ω² = 2c² E_num / T_num */
    if (T_num < 1e-30) return 0;
    return 2.0 * c2 * E_num / T_num;
}

/* ========== Square-well verification test ========== */

static void test_square_well(void) {
    printf("\n===== Square Well Verification Test =====\n");

    /* Parameters: K=1/2 (constant), V(r) = -V0 for r<a, 0 for r>a, μ=0 */
    double V0 = 10.0, a = 2.0;
    double K_const = 0.5;
    double c2 = 1.0;

    /* The Schrödinger-like equation is:
     * (ω²/c²) g = -(1/r²)d/dr[r²·2K·g'] + 2W·g
     * With K=1/2, this is: (ω²/c²) g = -g'' - (2/r)g' + 2W·g
     * For ℓ=0: -g'' - (2/r)g' - 2V0·g = ω²·g  (inside)
     *          -g'' - (2/r)g' = ω²·g             (outside)
     *
     * Let u = rg. Then u'' + (2V0+ω²)u = 0 inside, u'' + ω²u = 0 outside.
     * Since ω² < 0 for bound states (below continuum at 0):
     *   Inside: u'' + (2V0 - |ω²|)u = 0 → k² = 2V0 - |ω²| > 0
     *   Outside: u'' - |ω²|u = 0 → κ² = |ω²|
     *
     * u(r<a) = sin(kr), u(r>a) = Ce^{-κr}
     * Matching: k·cos(ka)/sin(ka) = -κ → k·cot(ka) = -κ
     * With k² + κ² = 2V0.
     */

    /* Analytical eigenvalues: k cot(ka) = -κ, k² + κ² = 2V0 */
    printf("V0=%.1f, a=%.1f, K=%.1f\n", V0, a, K_const);
    printf("Analytical eigenvalues (k cot(ka) = -κ, k² + κ² = %.1f):\n", 2*V0);

    /* Scan for solutions, skipping cot(ka) singularities at ka = nπ */
    double analytic[10];
    int n_analytic = 0;
    double dk = 0.0001;
    for (int trial = 0; trial < 100000; trial++) {
        double k = 0.001 + trial * dk;
        if (k*k >= 2*V0) break;
        /* Skip near singularities of cot(ka): sin(ka) ≈ 0 */
        if (fabs(sin(k*a)) < 0.01) continue;
        double k2 = k + dk;
        if (k2*k2 >= 2*V0) continue;
        if (fabs(sin(k2*a)) < 0.01) continue;

        double kappa = sqrt(2*V0 - k*k);
        double lhs = k * cos(k*a) / sin(k*a) + kappa;
        double kappa2 = sqrt(2*V0 - k2*k2);
        double lhs2 = k2 * cos(k2*a) / sin(k2*a) + kappa2;
        if (lhs * lhs2 < 0) {
            double klo = k, khi = k2;
            for (int bi = 0; bi < 100; bi++) {
                double km = 0.5*(klo+khi);
                double kpm = sqrt(2*V0 - km*km);
                double val = km * cos(km*a) / sin(km*a) + kpm;
                if (val * lhs < 0) khi = km; else klo = km;
            }
            double kf = 0.5*(klo+khi);
            /* k² = ω² + 2V0, κ² = -ω²; k² + κ² = 2V0 */
            double energy = kf*kf - 2*V0;
            analytic[n_analytic++] = energy;
            printf("  n=%d: k=%.6f, κ=%.6f, ω²=%.6f (binding=%.6f)\n",
                   n_analytic-1, kf, sqrt(-energy), energy, -energy);
            if (n_analytic >= 10) break;
        }
    }

    /* Build potential for numerical solver */
    int nr_test = 2000;
    double r_max_test = 15.0;
    Potential *pot = calloc(1, sizeof(Potential));
    pot->nr = nr_test;
    pot->dr = r_max_test / nr_test;
    pot->r_max = r_max_test;
    pot->r    = calloc(nr_test+1, sizeof(double));
    pot->K    = calloc(nr_test+1, sizeof(double));
    pot->W    = calloc(nr_test+1, sizeof(double));
    pot->beta = calloc(nr_test+1, sizeof(double));
    pot->Weff = calloc(nr_test+1, sizeof(double));

    for (int i = 0; i <= nr_test; i++) {
        double r = (i == 0) ? 1e-6 : i * pot->dr;
        pot->r[i] = r;
        pot->K[i] = K_const;
        pot->W[i] = (r < a) ? -V0 : 0.0;
        pot->beta[i] = 0.0;
        pot->Weff[i] = pot->W[i];
    }

    printf("\nNumerical eigenvalues (scan + bisection):\n");
    ShootParams sp = { .omega2 = 0, .c2 = c2, .pot = pot };
    double *gw = malloc((nr_test+1)*sizeof(double));
    double *gpw = malloc((nr_test+1)*sizeof(double));

    double num_evals[20];
    int n_num = find_all_eigenvalues(pot, 0, -2*V0, 0, 1e-10, c2, num_evals, 20);
    printf("  Found %d numerical eigenvalues, %d analytical\n", n_num, n_analytic);

    int n_compare = (n_num < n_analytic) ? n_num : n_analytic;
    for (int n = 0; n < n_compare; n++) {
        sp.omega2 = num_evals[n];
        shoot(&sp, 0, gw, gpw, nr_test+1);
        double om2_var = variational_omega2(pot, gw, gpw, nr_test, c2);
        double err = fabs(num_evals[n] - analytic[n]) / fabs(analytic[n]) * 100;
        printf("  n=%d: ω²=%.8f (analytic=%.8f, err=%.4f%%), nodes=%d, var_ω²=%.6f\n",
               n, num_evals[n], analytic[n], err, count_nodes(gw, nr_test+1), om2_var);
    }
    /* Report any unmatched */
    for (int n = n_compare; n < n_num; n++)
        printf("  n=%d: ω²=%.8f (extra numerical)\n", n, num_evals[n]);
    for (int n = n_compare; n < n_analytic; n++)
        printf("  n=%d: ω²=%.8f (missing numerical)\n", n, analytic[n]);

    free(gw); free(gpw); free_potential(pot);
    printf("===== End Square Well Test =====\n\n");
}

/* ========== Main ========== */

int main(int argc, char *argv[]) {
    const char *profile_file = NULL;
    double e_skyrme = 2.0;
    double gc = 1.0;         /* coupling constant */
    double mu = 1.0;         /* degenerate mass */
    double c2 = 1.0;         /* c² = 1 in natural units */
    int n_ang = 30;          /* angular quadrature points */
    int nr = 1000;           /* radial grid points */
    double r_max = 15.0;     /* radial cutoff */
    int scan = 0;            /* parameter scan mode */
    int test = 0;            /* run square-well test */
    int verbose = 1;
    int coupling_type = 0;   /* 0=Skyrme, 1=geometric (covariant derivative) */
    double g_cov = 1.0;      /* covariant coupling strength */

    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "-profile") && i+1 < argc) profile_file = argv[++i];
        else if (!strcmp(argv[i], "-e") && i+1 < argc) e_skyrme = atof(argv[++i]);
        else if (!strcmp(argv[i], "-gc") && i+1 < argc) gc = atof(argv[++i]);
        else if (!strcmp(argv[i], "-mu") && i+1 < argc) mu = atof(argv[++i]);
        else if (!strcmp(argv[i], "-gcov") && i+1 < argc) g_cov = atof(argv[++i]);
        else if (!strcmp(argv[i], "-nang") && i+1 < argc) n_ang = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-nr") && i+1 < argc) nr = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-rmax") && i+1 < argc) r_max = atof(argv[++i]);
        else if (!strcmp(argv[i], "-geom")) coupling_type = 1;
        else if (!strcmp(argv[i], "-skyrme")) coupling_type = 0;
        else if (!strcmp(argv[i], "-scan")) scan = 1;
        else if (!strcmp(argv[i], "-test")) test = 1;
        else if (!strcmp(argv[i], "-quiet")) verbose = 0;
        else {
            fprintf(stderr, "Usage: %s [-profile FILE] [-e E] [-gc GC] [-mu MU] [-gcov G]\n"
                    "       [-geom|-skyrme] [-nang N] [-nr N] [-rmax R] [-scan] [-test]\n",
                    argv[0]);
            return 1;
        }
    }

    /* Square-well test mode */
    if (test) {
        test_square_well();
        if (!profile_file) return 0;
    }

    if (!profile_file) {
        fprintf(stderr, "Error: -profile FILE required (or use -test for verification)\n");
        return 1;
    }

    /* Read soliton background */
    Profile prof;
    if (read_profile(profile_file, &prof) < 0) return 1;

    printf("# Bound state solver for degenerate modes\n");
    printf("# Profile: %s (%d points, r_max=%.2f)\n", profile_file, prof.n, prof.r[prof.n-1]);
    printf("# Coupling: %s", coupling_type ? "GEOMETRIC (covariant derivative)" : "SKYRME (E_{4,C})");
    if (coupling_type)
        printf(", g_cov=%.4f", g_cov);
    else
        printf(", e=%.4f", e_skyrme);
    printf("\n");
    printf("# gc=%.4f, mu=%.4f, c²=%.4f\n", gc, mu, c2);
    printf("# Angular quadrature: n_theta=%d, n_phi=%d\n", n_ang, 2*n_ang);
    printf("# Radial grid: nr=%d, dr=%.6f, r_max=%.2f\n", nr, r_max/nr, r_max);

    double gc2 = gc * gc;
    double mu2 = mu * mu;

    if (!scan) {
        /* Single-parameter run: compute potential and find eigenvalues */
        for (int mode = 0; mode <= 1; mode++) {
            const char *mode_name = (mode == 0) ? "P-mode (l=0)" : "J-mode (l=1)";
            printf("\n# ======= %s =======\n", mode_name);

            /* Build effective potential */
            printf("# Building potential... (n_ang=%d per radial point)\n", n_ang);
            Potential *pot = build_potential(&prof, e_skyrme, gc2, mu2, mode, nr, r_max, n_ang,
                                            coupling_type, g_cov);

            /* Print potential profile */
            if (verbose) {
                printf("# Potential profile:\n");
                printf("# %10s %14s %14s %14s %14s\n", "r", "K(r)", "W(r)", "beta(r)", "W_eff(r)");
                for (int ir = 0; ir <= pot->nr; ir += pot->nr/50) {
                    printf("# %10.4f %14.8e %14.8e %14.8e %14.8e\n",
                           pot->r[ir], pot->K[ir], pot->W[ir], pot->beta[ir], pot->Weff[ir]);
                }
            }

            /* Continuum threshold: ω² = c²μ² */
            double omega2_max = c2 * mu2;
            /* Search from deep negative to just below threshold */
            /* For repulsive potential, eigenvalues may be ω² > 0 but below threshold */
            double omega2_min = -20.0;  /* generous lower bound */

            /* Scan for potential minimum to estimate search range */
            double weff_min = 1e30;
            for (int ir = 0; ir <= pot->nr; ir++)
                if (pot->Weff[ir] < weff_min) weff_min = pot->Weff[ir];

            /* The eigenvalue must satisfy ω² ≥ 2K_min × W_eff_min roughly */
            omega2_min = 2.0 * weff_min * c2;
            if (omega2_min > -20.0) omega2_min = -20.0;

            printf("# Continuum threshold: ω² = c²μ² = %.6f\n", omega2_max);
            printf("# W_eff minimum: %.6f\n", weff_min);
            printf("# Search range: ω² ∈ [%.4f, %.4f]\n", omega2_min, omega2_max);

            /* Find all bound states */
            printf("#\n# Bound states:\n");
            printf("# %4s %14s %14s %6s %14s\n",
                   "n", "omega2", "binding", "nodes", "var_omega2");

            double evals[20];
            int found = find_all_eigenvalues(pot, mode, omega2_min, omega2_max,
                                              1e-10, c2, evals, 20);

            double *gw = malloc((nr+1)*sizeof(double));
            double *gpw = malloc((nr+1)*sizeof(double));
            ShootParams sp = { .omega2 = 0, .c2 = c2, .pot = pot };

            for (int n = 0; n < found; n++) {
                sp.omega2 = evals[n];
                shoot(&sp, mode, gw, gpw, nr+1);
                int nodes = count_nodes(gw, nr+1);
                double om2_var = variational_omega2(pot, gw, gpw, nr, c2);
                double binding = omega2_max - evals[n];

                printf("  %4d %14.8f %14.8f %6d %14.8f\n",
                       n, evals[n], binding, nodes, om2_var);

                /* Save ground state wavefunction */
                if (n == 0 && verbose) {
                    double norm = 0;
                    double dr_val = pot->dr;
                    for (int ir = 0; ir <= pot->nr; ir++) {
                        double w = dr_val * ((ir == 0 || ir == pot->nr) ? 0.5 : 1.0);
                        norm += w * pot->r[ir]*pot->r[ir] * gw[ir]*gw[ir];
                    }
                    norm = sqrt(4*M_PI*norm);
                    if (norm > 1e-20) {
                        printf("#\n# Ground state wavefunction (normalized):\n");
                        printf("# %10s %14s %14s %14s\n", "r", "g(r)", "K(r)", "W_eff(r)");
                        for (int ir = 0; ir <= pot->nr; ir += pot->nr/100) {
                            printf("# %10.4f %14.8e %14.8e %14.8e\n",
                                   pot->r[ir], gw[ir]/norm, pot->K[ir], pot->Weff[ir]);
                        }
                    }
                }
            }

            if (found == 0) {
                printf("#  NO BOUND STATES FOUND\n");
                printf("# Potential analysis:\n");
                double W_max = -1e30, W_min = 1e30;
                double r_max_pot = 0, r_min_pot = 0;
                for (int ir = 0; ir <= pot->nr; ir++) {
                    if (pot->Weff[ir] > W_max) { W_max = pot->Weff[ir]; r_max_pot = pot->r[ir]; }
                    if (pot->Weff[ir] < W_min) { W_min = pot->Weff[ir]; r_min_pot = pot->r[ir]; }
                }
                printf("#   W_eff max = %.6f at r=%.4f (barrier)\n", W_max, r_max_pot);
                printf("#   W_eff min = %.6f at r=%.4f (well)\n", W_min, r_min_pot);
                printf("#   Threshold = %.6f (μ²/2)\n", mu2/2.0);
                double deficit = mu2/2.0 - W_min;  /* how much well exceeds threshold */
                if (W_min >= mu2/2.0)
                    printf("#   Well is ABOVE threshold by %.6f — need deeper well\n", -deficit);
                else
                    printf("#   Well is BELOW threshold by %.6f — might bind with better params\n", deficit);
            }

            free(gw); free(gpw);
            free_potential(pot);
        }
    } else {
        /* Parameter scan mode */
        if (coupling_type == 1) {
            /* Geometric mode: scan over (g_cov, mu) at fixed gc */
            printf("\n# Parameter scan (GEOMETRIC): bound state count vs (g_cov, mu)\n");
            printf("# gc=%.4f (fixed)\n", gc);
            printf("# %8s %8s %6s %14s %6s %14s\n",
                   "g_cov", "mu", "P_cnt", "P_ground_E", "J_cnt", "J_ground_E");
            printf("# -------  -------  -----  -------------  -----  -------------\n");

            double gcov_vals[] = {0.5, 1.0, 2.0, 5.0, 10.0};
            double mu_vals[] = {0.1, 0.5, 1.0, 2.0, 5.0};
            int n_gcov = sizeof(gcov_vals)/sizeof(gcov_vals[0]);
            int n_mu = sizeof(mu_vals)/sizeof(mu_vals[0]);

            for (int ig = 0; ig < n_gcov; ig++) {
                for (int im = 0; im < n_mu; im++) {
                    double gcov_v = gcov_vals[ig], mu_v = mu_vals[im];
                    double mu2_v = mu_v*mu_v;
                    double omega2_max = c2 * mu2_v;

                    int P_count = 0, J_count = 0;
                    double P_ground = 0, J_ground = 0;

                    for (int mode = 0; mode <= 1; mode++) {
                        Potential *pot = build_potential(&prof, e_skyrme, gc2, mu2_v,
                                                        mode, nr, r_max, n_ang,
                                                        coupling_type, gcov_v);

                        double omega2_min = -20.0;
                        double weff_min = 1e30;
                        for (int ir = 0; ir <= pot->nr; ir++)
                            if (pot->Weff[ir] < weff_min) weff_min = pot->Weff[ir];
                        double om_min = 2.0*weff_min*c2;
                        if (om_min > -20.0) om_min = -20.0;
                        if (om_min < omega2_min) omega2_min = om_min;

                        double ev_list[20];
                        int count = find_all_eigenvalues(pot, mode, omega2_min, omega2_max,
                                                         1e-8, c2, ev_list, 20);
                        double ground_e = (count > 0) ? ev_list[0] : 0;

                        if (mode == 0) { P_count = count; P_ground = ground_e; }
                        else           { J_count = count; J_ground = ground_e; }
                        free_potential(pot);
                    }

                    printf("  %8.2f %8.2f %6d %14.6f %6d %14.6f\n",
                           gcov_v, mu_v, P_count, P_ground, J_count, J_ground);
                }
            }
        } else {
        printf("\n# Parameter scan: bound state count vs (gc, mu)\n");
        printf("# %8s %8s %6s %14s %6s %14s\n",
               "gc", "mu", "P_cnt", "P_ground_E", "J_cnt", "J_ground_E");
        printf("# -------  -------  -----  -------------  -----  -------------\n");

        double gc_vals[] = {0.5, 1.0, 2.0, 5.0, 10.0};
        double mu_vals[] = {0.1, 0.5, 1.0, 2.0, 5.0};
        int n_gc = sizeof(gc_vals)/sizeof(gc_vals[0]);
        int n_mu = sizeof(mu_vals)/sizeof(mu_vals[0]);

        for (int ig = 0; ig < n_gc; ig++) {
            for (int im = 0; im < n_mu; im++) {
                double gc_v = gc_vals[ig], mu_v = mu_vals[im];
                double gc2_v = gc_v*gc_v, mu2_v = mu_v*mu_v;
                double omega2_max = c2 * mu2_v;

                int P_count = 0, J_count = 0;
                double P_ground = 0, J_ground = 0;

                for (int mode = 0; mode <= 1; mode++) {
                    Potential *pot = build_potential(&prof, e_skyrme, gc2_v, mu2_v,
                                                    mode, nr, r_max, n_ang,
                                                    coupling_type, g_cov);

                    double omega2_min = -20.0;
                    double weff_min = 1e30;
                    for (int ir = 0; ir <= pot->nr; ir++)
                        if (pot->Weff[ir] < weff_min) weff_min = pot->Weff[ir];
                    double om_min = 2.0*weff_min*c2;
                    if (om_min > -20.0) om_min = -20.0;
                    if (om_min < omega2_min) omega2_min = om_min;

                    /* Find eigenvalues */
                    double ev_list[20];
                    int count = find_all_eigenvalues(pot, mode, omega2_min, omega2_max,
                                                     1e-8, c2, ev_list, 20);
                    double ground_e = (count > 0) ? ev_list[0] : 0;

                    if (mode == 0) { P_count = count; P_ground = ground_e; }
                    else           { J_count = count; J_ground = ground_e; }
                    free_potential(pot);
                }

                printf("  %8.2f %8.2f %6d %14.6f %6d %14.6f\n",
                       gc_v, mu_v, P_count, P_ground, J_count, J_ground);
            }
        }
        } /* end Skyrme scan else */
    }

    free_profile(&prof);
    return 0;
}
