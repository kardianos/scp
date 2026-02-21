/*
 * nullrotor_metric.c — BLV effective metric for null-rotor propagation
 *                      on the B=1 hedgehog Skyrmion background
 *
 * PHYSICS:
 * ========
 * The Barceló-Liberati-Visser (BLV) effective metric describes how small
 * perturbations propagate on a nonlinear field background. For a Lagrangian
 * L(Ψ, ∂_μΨ), the effective metric is:
 *
 *   G^{μν}_{AB} = ∂²L / ∂(∂_μΨ_A) ∂(∂_νΨ_B)
 *
 * evaluated on the background field. Different perturbation polarizations A
 * see different effective metrics.
 *
 * For the Skyrme Lagrangian L = L₂ + L₄:
 *
 *   L₂ = (1/2) η^{μν} Σ_A (∂_μq_A)(∂_νq_A)
 *   L₄ = (c₄/4) Σ_{μ<ν} |[R_μ, R_ν]|²   where R_μ = (∂_μq) q̃ / |q|²
 *
 * L₂ gives the flat metric: G₂^{μν} = η^{μν} δ_{AB}
 * L₄ gives corrections that depend on the background and polarization.
 *
 * KEY RESULT (derived analytically):
 * On the hedgehog z-axis, q₀ = ρ₀(cosf + sinf σ₃), the L₄ stiffness
 * for spatial propagation involves commutators [e_A q̃₀, R_i^bg]:
 *
 *   P₄^{zz}_A = (c₄/2ρ₀⁴) × [|[e_A q̃₀, R_x]|² + |[e_A q̃₀, R_y]|²]
 *   P₄^{xx}_A = (c₄/2ρ₀⁴) × [|[e_A q̃₀, R_z]|² + |[e_A q̃₀, R_y]|²]
 *   m₄_A     = (c₄/2c²ρ₀⁴) × Σ_i |[e_A q̃₀, R_i]|²
 *
 * Phase velocity: v²(direction) = c² × (1 + P₄) / (1 + c²m₄)
 *
 * The v2 result P/m = 2 (from the Sturm-Liouville reduction) is an artifact
 * of the angular integration + volume element. The LOCAL BLV metric has
 * P ≠ m for off-diagonal propagation (e.g., σ₁ mode propagating radially).
 *
 * OUTPUTS:
 * ========
 * - Phase velocities v(r, polarization, direction) / c
 * - Effective gravitational potential Φ(r)
 * - Deflection angle θ(b) for null-rotor ray-tracing
 */

#define _USE_MATH_DEFINES
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ========== Profile I/O ========== */

typedef struct {
    double *r, *f, *fp;
    int n;
    double dr, rho0, e_skyrme, c4;
} Profile;

static Profile *load_profile(const char *filename, double rho0, double e_skyrme)
{
    FILE *fp = fopen(filename, "r");
    if (!fp) { fprintf(stderr, "Cannot open %s\n", filename); return NULL; }

    int cap = 16384;
    double *r_arr  = malloc(cap * sizeof(double));
    double *f_arr  = malloc(cap * sizeof(double));
    double *fp_arr = malloc(cap * sizeof(double));
    int n = 0;
    char line[1024];

    while (fgets(line, sizeof(line), fp)) {
        if (line[0] == '#') continue;
        double rv, fv, fpv = 0, d1, d2, d3;
        int nc = sscanf(line, "%lf %lf %lf %lf %lf %lf",
                        &rv, &fv, &fpv, &d1, &d2, &d3);
        if (nc < 2) continue;
        if (n >= cap) {
            cap *= 2;
            r_arr  = realloc(r_arr,  cap * sizeof(double));
            f_arr  = realloc(f_arr,  cap * sizeof(double));
            fp_arr = realloc(fp_arr, cap * sizeof(double));
        }
        r_arr[n] = rv; f_arr[n] = fv; fp_arr[n] = fpv;
        n++;
    }
    fclose(fp);

    /* Compute f' if not provided */
    int has_fp = 0;
    for (int i = 0; i < n; i++)
        if (fabs(fp_arr[i]) > 1e-30) { has_fp = 1; break; }
    if (!has_fp && n > 2) {
        for (int i = 0; i < n; i++) {
            if (i == 0)
                fp_arr[i] = (f_arr[1] - f_arr[0]) / (r_arr[1] - r_arr[0]);
            else if (i == n - 1)
                fp_arr[i] = (f_arr[n-1] - f_arr[n-2]) / (r_arr[n-1] - r_arr[n-2]);
            else
                fp_arr[i] = (f_arr[i+1] - f_arr[i-1]) / (r_arr[i+1] - r_arr[i-1]);
        }
    }

    Profile *p = malloc(sizeof(Profile));
    p->r = r_arr; p->f = f_arr; p->fp = fp_arr;
    p->n = n;
    p->dr = (n > 1) ? r_arr[1] - r_arr[0] : 0.001;
    p->rho0 = rho0;
    p->e_skyrme = e_skyrme;
    p->c4 = 2.0 * rho0 * rho0 / (e_skyrme * e_skyrme);
    return p;
}

static void interp_profile(const Profile *p, double r,
                            double *f_out, double *fp_out)
{
    if (r <= p->r[0]) { *f_out = p->f[0]; *fp_out = p->fp[0]; return; }
    if (r >= p->r[p->n-1]) {
        *f_out = p->f[p->n-1]; *fp_out = p->fp[p->n-1]; return;
    }
    int lo = 0, hi = p->n - 1;
    while (hi - lo > 1) {
        int mid = (lo + hi) / 2;
        if (p->r[mid] <= r) lo = mid; else hi = mid;
    }
    double t = (r - p->r[lo]) / (p->r[hi] - p->r[lo]);
    *f_out  = p->f[lo]  + t * (p->f[hi]  - p->f[lo]);
    *fp_out = p->fp[lo] + t * (p->fp[hi] - p->fp[lo]);
}

static void free_profile(Profile *p) {
    free(p->r); free(p->f); free(p->fp); free(p);
}

/* ========== Quaternion algebra (minimal) ========== */

/* Quaternion: (s, v1, v2, v3) = s + v1·σ₁ + v2·σ₂ + v3·σ₃
 * Product: (a,b)(c,d) = (ac - b·d, ad + cb + b×d)
 * Reverse: (s, v) → (s, -v)
 * Commutator: [a, b] = ab - ba  (scalar part always zero)
 */
typedef struct { double s, v1, v2, v3; } Quat;

static Quat q_mul(Quat a, Quat b) {
    return (Quat){
        a.s*b.s - a.v1*b.v1 - a.v2*b.v2 - a.v3*b.v3,
        a.s*b.v1 + a.v1*b.s + a.v2*b.v3 - a.v3*b.v2,
        a.s*b.v2 + a.v2*b.s + a.v3*b.v1 - a.v1*b.v3,
        a.s*b.v3 + a.v3*b.s + a.v1*b.v2 - a.v2*b.v1
    };
}

static Quat __attribute__((unused)) q_rev(Quat a) {
    return (Quat){ a.s, -a.v1, -a.v2, -a.v3 };
}

static Quat q_sub(Quat a, Quat b) {
    return (Quat){ a.s-b.s, a.v1-b.v1, a.v2-b.v2, a.v3-b.v3 };
}

static double q_norm2(Quat a) {
    return a.s*a.s + a.v1*a.v1 + a.v2*a.v2 + a.v3*a.v3;
}

static Quat q_comm(Quat a, Quat b) {
    /* [a, b] = ab - ba */
    Quat ab = q_mul(a, b);
    Quat ba = q_mul(b, a);
    return q_sub(ab, ba);
}

/* ========== BLV effective metric computation ========== */

/*
 * At the z-axis of the hedgehog (θ=0, r̂ = ẑ):
 *   q₀ = ρ₀(cosf + sinf σ₃)
 *   q̃₀ = ρ₀(cosf - sinf σ₃)
 *
 * Background right currents R_i = (∂_i q₀) q̃₀ / ρ₀²:
 *   R_z = f' σ₃
 *   R_x = (sinf/r)(cosf σ₁ + sinf σ₂)
 *   R_y = (sinf/r)(cosf σ₂ - sinf σ₁)
 *
 * Perturbation directions: e_A for A=0,1,2,3 → (1, σ₁, σ₂, σ₃)
 * Current perturbation: δR = (δ/ρ₀²) e_A q̃₀
 *
 * L₄ stiffness for z-propagation:
 *   P₄^{zz}_A = (c₄/2ρ₀⁴) [|[e_A q̃₀, R_x]|² + |[e_A q̃₀, R_y]|²]
 *
 * L₄ stiffness for x-propagation:
 *   P₄^{xx}_A = (c₄/2ρ₀⁴) [|[e_A q̃₀, R_z]|² + |[e_A q̃₀, R_y]|²]
 *
 * L₄ temporal inertia:
 *   m₄_A = (c₄/2ρ₀⁴) [|[e_A q̃₀, R_x]|² + |[e_A q̃₀, R_y]|² + |[e_A q̃₀, R_z]|²]
 *
 * Phase velocity: v² = c² × (1 + P₄) / (1 + P₄ + correction)
 * where the temporal inertia includes all three spatial directions.
 */

typedef struct {
    double P4_zz;   /* L₄ stiffness for z-propagation */
    double P4_xx;   /* L₄ stiffness for x-propagation */
    double m4;      /* L₄ temporal inertia (×c²) */
    double vz2;     /* v²_z / c² */
    double vx2;     /* v²_x / c² */
} BLVMetric;

static BLVMetric __attribute__((unused)) compute_blv_metric(double r, double f, double fp,
                                     double c4, double rho0)
{
    double sf = sin(f), cf = cos(f);
    double rho2 = rho0 * rho0;
    double rho4 = rho2 * rho2;
    double r_safe = (r < 1e-10) ? 1e-10 : r;

    /* Background right currents (pure quaternion part) */
    Quat Rz = {0, 0, 0, fp};
    Quat Rx = {0, sf*cf/r_safe, sf*sf/r_safe, 0};
    Quat Ry = {0, -sf*sf/r_safe, sf*cf/r_safe, 0};

    /* q̃₀ = ρ₀(cosf - sinf σ₃) */
    Quat qtilde = {rho0*cf, 0, 0, -rho0*sf};

    /* Perturbation directions */
    Quat eA[4] = {
        {1, 0, 0, 0},   /* e₀ = scalar */
        {0, 1, 0, 0},   /* e₁ = σ₁ */
        {0, 0, 1, 0},   /* e₂ = σ₂ */
        {0, 0, 0, 1}    /* e₃ = σ₃ */
    };

    BLVMetric met[4];

    for (int A = 0; A < 4; A++) {
        /* δR = e_A q̃₀ (quaternion product) */
        Quat dR = q_mul(eA[A], qtilde);

        /* Commutators [δR, R_i] */
        Quat comm_x = q_comm(dR, Rx);
        Quat comm_y = q_comm(dR, Ry);
        Quat comm_z = q_comm(dR, Rz);

        double Cx = q_norm2(comm_x);
        double Cy = q_norm2(comm_y);
        double Cz = q_norm2(comm_z);

        /* L₄ stiffness contributions (factor c₄/(2ρ₀⁴) from the formula) */
        double fac = c4 / (2.0 * rho4);

        met[A].P4_zz = fac * (Cx + Cy);
        met[A].P4_xx = fac * (Cz + Cy);
        met[A].m4    = fac * (Cx + Cy + Cz);  /* temporal inertia ×c² */

        /* Phase velocities:
         * v² = c² × (1 + P₄^{dir}) / (1 + m₄)
         * where m₄ already has the c² factor incorporated
         */
        met[A].vz2 = (1.0 + met[A].P4_zz) / (1.0 + met[A].m4);
        met[A].vx2 = (1.0 + met[A].P4_xx) / (1.0 + met[A].m4);
    }

    /* For output, return the physically most relevant mode.
     * We'll print all four in main(). Return the σ₁ (bivector transverse). */
    return met[1];
}

/* ========== Ray-tracing for deflection angle ========== */

/*
 * For a null-rotor passing the soliton at impact parameter b,
 * the deflection angle in the eikonal limit is:
 *
 *   θ(b) = -∫_{-∞}^{+∞} (∂n/∂r_⊥) ds
 *
 * where n(r) = c/v(r) is the refractive index and s is the path length.
 * For small deflections, the path is approximately straight (z-axis)
 * at distance b from the soliton.
 *
 * n(r) = c/v(r) = √((1 + m₄)/(1 + P₄))
 *
 * For a ray along z at distance b:
 *   r = √(b² + z²)
 *   ∂n/∂b = (dn/dr) × (b/r)
 *
 * θ(b) = -∫_{-∞}^{∞} (dn/dr)(b/r) dz
 *       = -b ∫₀^∞ (dn/dr)/r × 2dz   [symmetry]
 *       using z = b sinh(t): dz = b cosh(t) dt, r = b cosh(t)
 *       = -2 ∫₀^{t_max} (dn/dr)|_{r=b cosh(t)} dt
 */

static double compute_deflection(const Profile *prof, double b,
                                  int polarization, int prop_direction __attribute__((unused)))
{
    /* prop_direction: 0 = z (radial), 1 = x (transverse) */
    /* For a ray propagating along z at distance b from soliton,
     * the relevant speed is v_x (transverse to the soliton-ray line)
     * for the PERPENDICULAR component, but the dominant effect is from
     * the gradient of n along the radial direction from the soliton. */

    double c4 = prof->c4;
    double rho0 = prof->rho0;
    double rho4 = rho0*rho0*rho0*rho0;

    int nt = 2000;
    double t_max = 8.0;  /* cosh(8) ≈ 1490, so r_max ≈ 1490*b */
    double dt = t_max / nt;

    double theta = 0.0;

    for (int it = 0; it < nt; it++) {
        double t = (it + 0.5) * dt;
        double r = b * cosh(t);

        double f, fp;
        interp_profile(prof, r, &f, &fp);

        double sf = sin(f), cf = cos(f);
        double r_safe = (r < 1e-10) ? 1e-10 : r;

        /* Compute refractive index n(r) for the chosen polarization */
        Quat Rz = {0, 0, 0, fp};
        Quat Rx = {0, sf*cf/r_safe, sf*sf/r_safe, 0};
        Quat Ry = {0, -sf*sf/r_safe, sf*cf/r_safe, 0};
        Quat qtilde = {rho0*cf, 0, 0, -rho0*sf};

        Quat eA[4] = {{1,0,0,0}, {0,1,0,0}, {0,0,1,0}, {0,0,0,1}};
        Quat dR = q_mul(eA[polarization], qtilde);

        double Cx = q_norm2(q_comm(dR, Rx));
        double Cy = q_norm2(q_comm(dR, Ry));
        double Cz = q_norm2(q_comm(dR, Rz));

        double fac = c4 / (2.0 * rho4);
        double m4 = fac * (Cx + Cy + Cz);

        /* For a ray along z at distance b, the speed that matters for
         * bending is the TRANSVERSE speed v_⊥ at the closest approach.
         * The refractive index for the perpendicular direction is
         * n_⊥ = √((1 + m₄)/(1 + P₄^{perp})).
         *
         * For a ray along z, the perpendicular direction at any point
         * is in the x-y plane. By cylindrical symmetry around the z-axis,
         * we need the component along the radial (b̂) direction in the x-y plane.
         *
         * Actually, for a photon propagating along z at distance b from the
         * soliton, the relevant refractive index is for z-PROPAGATION at
         * that radial distance. The bending comes from the TRANSVERSE gradient
         * of this refractive index.
         *
         * n_z(r) = √((1 + m₄) / (1 + P₄^{zz}))  ... for z-propagation
         * But this is on the z-axis. For a general point at distance r from
         * the soliton center, we need to rotate the result.
         *
         * By spherical symmetry of the hedgehog, at any point at distance r:
         * n_parallel(r) = refractive index for propagation along r̂
         * n_perp(r) = refractive index for propagation perpendicular to r̂
         *
         * For a photon at distance r, propagating tangentially (perpendicular
         * to r̂ from soliton center), the bending is determined by dn_perp/dr.
         *
         * Using our z-axis results: n_parallel corresponds to v_z results,
         * and n_perp corresponds to v_x results, when the photon is on the z-axis.
         */

        /* For the propagation along z, at this distance r from center,
         * we use the fact that the photon trajectory is approximately along
         * a direction making angle α = atan(b/z) with the radial direction.
         * In the eikonal limit, the relevant refractive index is n(r) for
         * the propagation direction.
         *
         * For SIMPLICITY and to get the LEADING effect, use the average:
         * The isotropic part of n determines the monopole deflection.
         * n_iso = (n_z + 2*n_x) / 3  [angular average]
         *
         * Actually, the correct formula for ray bending is:
         * θ = -∫ ∇_⊥ n ds = -∫ (dn/dr)(b/r) dz
         *
         * This uses the refractive index for the propagation direction
         * (along the ray, approximately z) evaluated at distance r from center.
         */

        /* n for z-propagation: sqrt((1+m4)/(1+P4_zz)) */
        double P4_zz = fac * (Cx + Cy);

        /* Use the propagation direction n_z for the z-directed ray */
        double n_r = sqrt((1.0 + m4) / (1.0 + P4_zz));

        /* Now compute dn/dr by finite difference */
        double dr_fd = 0.001;
        double r_plus = r + dr_fd;
        double f_p, fp_p;
        interp_profile(prof, r_plus, &f_p, &fp_p);

        double sf_p = sin(f_p), cf_p = cos(f_p);
        Quat Rz_p = {0, 0, 0, fp_p};
        Quat Rx_p = {0, sf_p*cf_p/r_plus, sf_p*sf_p/r_plus, 0};
        Quat Ry_p = {0, -sf_p*sf_p/r_plus, sf_p*cf_p/r_plus, 0};
        Quat qt_p = {rho0*cf_p, 0, 0, -rho0*sf_p};
        Quat dR_p = q_mul(eA[polarization], qt_p);

        double Cx_p = q_norm2(q_comm(dR_p, Rx_p));
        double Cy_p = q_norm2(q_comm(dR_p, Ry_p));
        double Cz_p = q_norm2(q_comm(dR_p, Rz_p));
        double m4_p = fac * (Cx_p + Cy_p + Cz_p);
        double P4_zz_p = fac * (Cx_p + Cy_p);
        double n_rp = sqrt((1.0 + m4_p) / (1.0 + P4_zz_p));

        double dn_dr = (n_rp - n_r) / dr_fd;

        theta += dn_dr * dt;
    }

    /* Factor of -2 (symmetry of path, minus sign for bending toward higher n) */
    return -2.0 * theta;
}

/* ========== Main ========== */

int main(int argc, char *argv[])
{
    const char *profile_file = NULL;
    const char *outdir = ".";
    double rho0 = 1.0, e_skyrme = 1.0;

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-profile") == 0 && i+1 < argc)
            profile_file = argv[++i];
        else if (strcmp(argv[i], "-outdir") == 0 && i+1 < argc)
            outdir = argv[++i];
        else if (strcmp(argv[i], "-e") == 0 && i+1 < argc)
            e_skyrme = atof(argv[++i]);
        else if (strcmp(argv[i], "-rho0") == 0 && i+1 < argc)
            rho0 = atof(argv[++i]);
        else {
            fprintf(stderr, "Usage: %s -profile <file> [-outdir dir] "
                    "[-e val] [-rho0 val]\n", argv[0]);
            return 1;
        }
    }

    if (!profile_file) {
        fprintf(stderr, "Must specify -profile <file>\n");
        return 1;
    }

    Profile *prof = load_profile(profile_file, rho0, e_skyrme);
    if (!prof) return 1;

    double c4 = prof->c4;
    double rho4 = rho0*rho0*rho0*rho0;

    printf("===== Null-Rotor BLV Effective Metric =====\n");
    printf("Profile: %s (%d points, r_max=%.2f)\n",
           profile_file, prof->n, prof->r[prof->n-1]);
    printf("Parameters: rho0=%.4f, e=%.4f, c4=%.6f\n\n", rho0, e_skyrme, c4);

    /* ===== Step 1: Compute BLV metric at each radius ===== */
    printf("===== Phase velocities v/c at each radius =====\n");
    printf("# r  v_z(e0)  v_x(e0)  v_z(e1)  v_x(e1)  "
           "v_z(e2)  v_x(e2)  v_z(e3)  v_x(e3)\n");

    const char *pol_names[] = {"scalar(e0)", "sigma1(e1)",
                                "sigma2(e2)", "sigma3(e3)"};

    /* Output file for metric data */
    char fname[512];
    snprintf(fname, sizeof(fname), "%s/blv_metric.dat", outdir);
    FILE *fout = fopen(fname, "w");
    if (!fout) { fprintf(stderr, "Cannot open %s\n", fname); return 1; }
    fprintf(fout, "# BLV effective metric for null-rotor propagation\n");
    fprintf(fout, "# r  vz_e0/c  vx_e0/c  vz_e1/c  vx_e1/c  "
            "vz_e2/c  vx_e2/c  vz_e3/c  vx_e3/c  "
            "n_z(e1)  n_x(e1)  Phi_z(e1)/c2\n");

    /* Track extremes */
    double vmin_all = 1.0, vmax_all = 1.0;
    double phi_min = 0.0;  /* effective potential minimum */

    int nr_out = 500;
    double r_max = prof->r[prof->n - 1];
    double dr_out = r_max / nr_out;

    Quat eA[4] = {{1,0,0,0}, {0,1,0,0}, {0,0,1,0}, {0,0,0,1}};

    for (int ir = 0; ir <= nr_out; ir++) {
        double r = (ir == 0) ? 0.01 : ir * dr_out;
        double f, fp;
        interp_profile(prof, r, &f, &fp);

        double sf = sin(f), cf = cos(f);
        double r_safe = (r < 1e-10) ? 1e-10 : r;

        Quat Rz = {0, 0, 0, fp};
        Quat Rx = {0, sf*cf/r_safe, sf*sf/r_safe, 0};
        Quat Ry = {0, -sf*sf/r_safe, sf*cf/r_safe, 0};
        Quat qtilde = {rho0*cf, 0, 0, -rho0*sf};

        double vz[4], vx[4];

        for (int A = 0; A < 4; A++) {
            Quat dR = q_mul(eA[A], qtilde);
            double Cx = q_norm2(q_comm(dR, Rx));
            double Cy = q_norm2(q_comm(dR, Ry));
            double Cz = q_norm2(q_comm(dR, Rz));

            double fac = c4 / (2.0 * rho4);
            double P4_zz = fac * (Cx + Cy);
            double P4_xx = fac * (Cz + Cy);
            double m4    = fac * (Cx + Cy + Cz);

            vz[A] = sqrt((1.0 + P4_zz) / (1.0 + m4));
            vx[A] = sqrt((1.0 + P4_xx) / (1.0 + m4));

            if (vz[A] < vmin_all) vmin_all = vz[A];
            if (vx[A] < vmin_all) vmin_all = vx[A];
            if (vz[A] > vmax_all) vmax_all = vz[A];
            if (vx[A] > vmax_all) vmax_all = vx[A];
        }

        /* Refractive index and effective potential for σ₁ mode */
        double n_z1 = 1.0 / vz[1];
        double n_x1 = 1.0 / vx[1];
        double phi_z1 = (vz[1] * vz[1] - 1.0) / 2.0;  /* Φ/c² ≈ (v²-c²)/(2c²) */
        if (phi_z1 < phi_min) phi_min = phi_z1;

        fprintf(fout, "%.6f  %.10f  %.10f  %.10f  %.10f  "
                "%.10f  %.10f  %.10f  %.10f  "
                "%.10f  %.10f  %.10e\n",
                r, vz[0], vx[0], vz[1], vx[1],
                vz[2], vx[2], vz[3], vx[3],
                n_z1, n_x1, phi_z1);

        /* Print a few representative radii */
        if (ir % (nr_out/10) == 0 || ir == 0) {
            printf("r=%.3f: ", r);
            for (int A = 0; A < 4; A++)
                printf(" %s z=%.6f x=%.6f", pol_names[A], vz[A], vx[A]);
            printf("\n");
        }
    }
    fclose(fout);
    printf("\nMetric data written to %s\n", fname);

    /* ===== Summary ===== */
    printf("\n===== Summary =====\n");
    printf("Speed range: v/c ∈ [%.8f, %.8f]\n", vmin_all, vmax_all);
    printf("Minimum effective potential Φ_min/c² = %.6e\n", phi_min);

    if (fabs(vmin_all - 1.0) < 1e-10 && fabs(vmax_all - 1.0) < 1e-10) {
        printf("\n** NULL RESULT: v = c for all modes to machine precision.\n");
        printf("   No gravitational lensing from L₂ + L₄.\n");
    } else {
        printf("\n** NON-TRIVIAL RESULT: Phase velocity deviates from c!\n");
        printf("   Maximum deviation: %.6e\n",
               fmax(fabs(vmin_all - 1.0), fabs(vmax_all - 1.0)));

        /* ===== Step 2: Ray-tracing ===== */
        printf("\n===== Deflection angles =====\n");
        printf("# b  theta_e0  theta_e1  theta_e2  theta_e3  (radians)\n");

        snprintf(fname, sizeof(fname), "%s/blv_deflection.dat", outdir);
        FILE *fdefl = fopen(fname, "w");
        fprintf(fdefl, "# Deflection angle vs impact parameter\n");
        fprintf(fdefl, "# b  theta_e0  theta_e1  theta_e2  theta_e3\n");

        double b_values[] = {0.5, 1.0, 1.5, 2.0, 3.0, 5.0, 7.0, 10.0, 15.0, 20.0};
        int nb = sizeof(b_values) / sizeof(b_values[0]);

        for (int ib = 0; ib < nb; ib++) {
            double b = b_values[ib];
            double theta[4];
            for (int A = 0; A < 4; A++) {
                theta[A] = compute_deflection(prof, b, A, 0);
            }
            printf("b=%.1f: e0=%.6e  e1=%.6e  e2=%.6e  e3=%.6e rad\n",
                   b, theta[0], theta[1], theta[2], theta[3]);
            fprintf(fdefl, "%.6f  %.10e  %.10e  %.10e  %.10e\n",
                    b, theta[0], theta[1], theta[2], theta[3]);
        }
        fclose(fdefl);
        printf("Deflection data written to %s\n", fname);

        /* Convert to physical units */
        double code_to_fm = 0.5624;  /* 1 code length = 0.5624 fm */
        double code_to_MeV = 9.098;
        double E_sol = 103.13;
        double M_sol_MeV = E_sol * code_to_MeV;

        printf("\n===== Physical interpretation =====\n");
        printf("Soliton mass: %.1f MeV (= proton mass)\n", M_sol_MeV);
        printf("1 code length = %.4f fm\n", code_to_fm);

        /* GR prediction for comparison */
        double G_code = 5.9e-39;  /* Newton's G in code units */
        printf("\nGR prediction: θ_GR = 4GM/(c²b)\n");
        for (int ib = 0; ib < nb; ib++) {
            double b = b_values[ib];
            double theta_GR = 4.0 * G_code * E_sol / b;
            printf("b=%.1f: θ_GR = %.6e rad (our result / GR = %.3e)\n",
                   b, theta_GR,
                   compute_deflection(prof, b, 1, 0) / (theta_GR + 1e-100));
        }
    }

    /* ===== Cross-check: reproduce v2 SL P/m result ===== */
    printf("\n===== Cross-check: Sturm-Liouville P/m at each r =====\n");
    printf("(These are the INTEGRATED coefficients, not the local BLV metric)\n");
    printf("# r  P_SL  m_SL  P/m\n");

    for (int ir = 0; ir <= 10; ir++) {
        double r = (ir == 0) ? 0.01 : ir * r_max / 10.0;
        double f, fp;
        interp_profile(prof, r, &f, &fp);
        double sf2 = sin(f) * sin(f);
        double P_SL = 2*r*r + 4*c4*sf2;
        double m_SL = r*r + 2*c4*sf2;
        printf("r=%.3f: P=%.6f  m=%.6f  P/m=%.8f\n",
               r, P_SL, m_SL, P_SL / (m_SL + 1e-30));
    }

    free_profile(prof);
    return 0;
}
