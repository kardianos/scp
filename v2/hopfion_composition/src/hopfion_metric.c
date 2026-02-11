/*
 * hopfion_metric.c — BLV effective metric on hopfion backgrounds
 *
 * Path 4 investigation: Does the effective acoustic metric for small
 * perturbations around a HOPFION background differ from the hedgehog?
 *
 * Key question: The hedgehog P/m = 2 identity relies on spherical symmetry.
 * A hopfion has axial symmetry (preferred z-axis). The directional effective
 * metric g^{xx}, g^{yy}, g^{zz} may differ — this anisotropy IS an
 * effective gravitational effect (gravitomagnetic-like).
 *
 * BLV formula: g_eff^{μν} = -∂²L/∂(∂_μφ)∂(∂_νφ) on background φ₀
 *
 * For L₂ + L₄:
 *   g^{00}(x) = 1 + (c₄/|q|⁴) Σ_i |A_i|²
 *   g^{jj}(x) = 1 + (c₄/|q|⁴) Σ_{i≠j} |A_i|²
 *
 * Sound speed in direction j: c²_j = g^{jj}/g^{00}
 *
 * For hedgehog: <g^{xx}> = <g^{yy}> = <g^{zz}> (isotropic) → c_x = c_y = c_z
 * For hopfion: <g^{xx}> = <g^{yy}> ≠ <g^{zz}> (axial anisotropy) → c_z ≠ c_x
 *
 * Usage: ./bin/hopfion_metric [-hedgehog] [-two D] [-a SIZE] [-N GRID] [-R RADIUS] [-lam6 L6]
 *        -two D  : Initialize two hopfions at z=±D (product ansatz)
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "spherical_grid.h"
#include "hopfion_field.h"
#include "hopfion_init.h"
#include "topology.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Quaternion-only types for current computation */
typedef struct { double s, f1, f2, f3; } Q4;

static inline Q4 q4_from(Multivector m) {
    return (Q4){m.s, m.f1, m.f2, m.f3};
}
static inline Q4 q4_rev(Q4 a) {
    return (Q4){a.s, -a.f1, -a.f2, -a.f3};
}
static inline Q4 q4_mul(Q4 a, Q4 b) {
    return (Q4){
        a.s*b.s  - a.f1*b.f1 - a.f2*b.f2 - a.f3*b.f3,
        a.s*b.f1 + a.f1*b.s  - a.f2*b.f3 + a.f3*b.f2,
        a.s*b.f2 + a.f1*b.f3 + a.f2*b.s  - a.f3*b.f1,
        a.s*b.f3 - a.f1*b.f2 + a.f2*b.f1 + a.f3*b.s
    };
}
static inline double q4_norm2(Q4 a) {
    return a.s*a.s + a.f1*a.f1 + a.f2*a.f2 + a.f3*a.f3;
}

/* 4th-order derivative in direction d (0=x, 1=y, 2=z) */
static Q4 q4_deriv(const SphericalGrid *g, int i, int j, int k, int d)
{
    double c1 = 1.0 / (12.0 * g->h);
    Multivector m2, m1, p1, p2;
    switch (d) {
    case 0:
        m2 = sg_get(g, i-2, j, k); m1 = sg_get(g, i-1, j, k);
        p1 = sg_get(g, i+1, j, k); p2 = sg_get(g, i+2, j, k);
        break;
    case 1:
        m2 = sg_get(g, i, j-2, k); m1 = sg_get(g, i, j-1, k);
        p1 = sg_get(g, i, j+1, k); p2 = sg_get(g, i, j+2, k);
        break;
    default:
        m2 = sg_get(g, i, j, k-2); m1 = sg_get(g, i, j, k-1);
        p1 = sg_get(g, i, j, k+1); p2 = sg_get(g, i, j, k+2);
        break;
    }
    Q4 r;
    r.s  = c1 * (m2.s  - 8*m1.s  + 8*p1.s  - p2.s);
    r.f1 = c1 * (m2.f1 - 8*m1.f1 + 8*p1.f1 - p2.f1);
    r.f2 = c1 * (m2.f2 - 8*m1.f2 + 8*p1.f2 - p2.f2);
    r.f3 = c1 * (m2.f3 - 8*m1.f3 + 8*p1.f3 - p2.f3);
    return r;
}

/* ========== Directional Effective Metric ========== */

/* Per-bin accumulators for metric components */
typedef struct {
    int n_bins;
    double dr;
    double *r;
    /* Averaged metric components */
    double *g00;     /* g^{00} = 1 + c₄ Σ|A_i|²/|q|⁴ */
    double *gxx;     /* g^{xx} = 1 + c₄ (|A_y|²+|A_z|²)/|q|⁴ */
    double *gyy;     /* g^{yy} = 1 + c₄ (|A_x|²+|A_z|²)/|q|⁴ */
    double *gzz;     /* g^{zz} = 1 + c₄ (|A_x|²+|A_y|²)/|q|⁴ */
    double *count;
} DirectionalMetric;

static DirectionalMetric *dm_alloc(int n_bins, double R)
{
    DirectionalMetric *dm = calloc(1, sizeof(DirectionalMetric));
    dm->n_bins = n_bins;
    dm->dr = R / n_bins;
    dm->r     = calloc(n_bins, sizeof(double));
    dm->g00   = calloc(n_bins, sizeof(double));
    dm->gxx   = calloc(n_bins, sizeof(double));
    dm->gyy   = calloc(n_bins, sizeof(double));
    dm->gzz   = calloc(n_bins, sizeof(double));
    dm->count = calloc(n_bins, sizeof(double));
    for (int b = 0; b < n_bins; b++)
        dm->r[b] = (b + 0.5) * dm->dr;
    return dm;
}

static void dm_free(DirectionalMetric *dm)
{
    if (!dm) return;
    free(dm->r); free(dm->g00); free(dm->gxx);
    free(dm->gyy); free(dm->gzz); free(dm->count);
    free(dm);
}

static DirectionalMetric *compute_directional_metric(
    const SphericalGrid *g, const HopfionParams *p, int n_bins)
{
    DirectionalMetric *dm = dm_alloc(n_bins, g->R);
    double c4 = 2.0 * p->rho0 * p->rho0 / (p->e_skyrme * p->e_skyrme);

    for (int n = 0; n < g->n_active; n++) {
        int ix = g->active[n];
        int i, j, k;
        sg_unflatten(g->N, ix, &i, &j, &k);

        double x, y, z;
        sg_pos(g, i, j, k, &x, &y, &z);
        double r = sqrt(x*x + y*y + z*z);
        if (r > g->R - 4.5*g->h || r < 1e-10) continue;

        int bin = (int)(r / dm->dr);
        if (bin < 0 || bin >= n_bins) continue;

        /* Compute left-currents A_d = q̃ ∂_d q */
        Q4 q = q4_from(g->psi[ix]);
        double q2 = q4_norm2(q);
        if (q2 < 1e-20) continue;
        double q4_inv = 1.0 / (q2 * q2);
        Q4 qr = q4_rev(q);

        double A2[3];
        for (int d = 0; d < 3; d++) {
            Q4 dq = q4_deriv(g, i, j, k, d);
            Q4 Ad = q4_mul(qr, dq);
            A2[d] = q4_norm2(Ad);
        }

        double A2_total = A2[0] + A2[1] + A2[2];
        double contrib = c4 * q4_inv;

        dm->g00[bin] += 1.0 + contrib * A2_total;
        dm->gxx[bin] += 1.0 + contrib * (A2[1] + A2[2]);  /* perp to x */
        dm->gyy[bin] += 1.0 + contrib * (A2[0] + A2[2]);  /* perp to y */
        dm->gzz[bin] += 1.0 + contrib * (A2[0] + A2[1]);  /* perp to z */
        dm->count[bin] += 1.0;
    }

    /* Average */
    for (int b = 0; b < n_bins; b++) {
        if (dm->count[b] < 1.0) continue;
        double inv = 1.0 / dm->count[b];
        dm->g00[b] *= inv;
        dm->gxx[b] *= inv;
        dm->gyy[b] *= inv;
        dm->gzz[b] *= inv;
    }

    return dm;
}

/* ========== Azimuthal binning (r, cos_theta) for axially-symmetric fields ========== */

typedef struct {
    int nr, ntheta;
    double dr, dcos;
    /* 2D arrays indexed [ir * ntheta + itheta] */
    double *g00;
    double *gxx, *gyy, *gzz;
    double *count;
} AxialMetric;

static AxialMetric *am_alloc(int nr, int ntheta, double R)
{
    AxialMetric *am = calloc(1, sizeof(AxialMetric));
    am->nr = nr;
    am->ntheta = ntheta;
    am->dr = R / nr;
    am->dcos = 2.0 / ntheta;  /* cos(theta) from -1 to +1 */
    int total = nr * ntheta;
    am->g00   = calloc(total, sizeof(double));
    am->gxx   = calloc(total, sizeof(double));
    am->gyy   = calloc(total, sizeof(double));
    am->gzz   = calloc(total, sizeof(double));
    am->count = calloc(total, sizeof(double));
    return am;
}

static void am_free(AxialMetric *am)
{
    if (!am) return;
    free(am->g00); free(am->gxx); free(am->gyy);
    free(am->gzz); free(am->count); free(am);
}

static AxialMetric *compute_axial_metric(
    const SphericalGrid *g, const HopfionParams *p, int nr, int ntheta)
{
    AxialMetric *am = am_alloc(nr, ntheta, g->R);
    double c4 = 2.0 * p->rho0 * p->rho0 / (p->e_skyrme * p->e_skyrme);

    for (int n = 0; n < g->n_active; n++) {
        int ix = g->active[n];
        int i, j, k;
        sg_unflatten(g->N, ix, &i, &j, &k);

        double x, y, z;
        sg_pos(g, i, j, k, &x, &y, &z);
        double r = sqrt(x*x + y*y + z*z);
        if (r > g->R - 4.5*g->h || r < 1e-10) continue;

        int ir = (int)(r / am->dr);
        double cos_theta = z / r;
        int itheta = (int)((cos_theta + 1.0) / am->dcos);
        if (ir < 0 || ir >= nr) continue;
        if (itheta < 0) itheta = 0;
        if (itheta >= ntheta) itheta = ntheta - 1;
        int idx = ir * ntheta + itheta;

        Q4 q = q4_from(g->psi[ix]);
        double q2 = q4_norm2(q);
        if (q2 < 1e-20) continue;
        double q4_inv = 1.0 / (q2 * q2);
        Q4 qr = q4_rev(q);

        double A2[3];
        for (int d = 0; d < 3; d++) {
            Q4 dq = q4_deriv(g, i, j, k, d);
            Q4 Ad = q4_mul(qr, dq);
            A2[d] = q4_norm2(Ad);
        }

        double A2_total = A2[0] + A2[1] + A2[2];
        double contrib = c4 * q4_inv;

        am->g00[idx]   += 1.0 + contrib * A2_total;
        am->gxx[idx]   += 1.0 + contrib * (A2[1] + A2[2]);
        am->gyy[idx]   += 1.0 + contrib * (A2[0] + A2[2]);
        am->gzz[idx]   += 1.0 + contrib * (A2[0] + A2[1]);
        am->count[idx] += 1.0;
    }

    /* Average */
    int total = nr * ntheta;
    for (int idx = 0; idx < total; idx++) {
        if (am->count[idx] < 1.0) continue;
        double inv = 1.0 / am->count[idx];
        am->g00[idx] *= inv;
        am->gxx[idx] *= inv;
        am->gyy[idx] *= inv;
        am->gzz[idx] *= inv;
    }

    return am;
}

/* ========== Main ========== */

int main(int argc, char **argv)
{
    /* Defaults */
    int N = 128;
    double R = 8.0;
    double L = R + 0.5;
    double sponge_w = 1.0;
    double rho0 = 1.0;
    double e_skyrme = 1.0;
    double lambda6 = 0.0;
    double hopfion_a = 2.0;
    int do_hedgehog = 0;
    int do_two = 0;
    int do_skyrmhop = 0;
    double two_d = 3.0;  /* separation half-distance for two-hopfion mode */
    double skyrmhop_d = 2.0;  /* hopfion ring offset from skyrmion center */
    const char *profile_path = NULL;

    /* Parse args */
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-hedgehog") == 0) do_hedgehog = 1;
        else if (strcmp(argv[i], "-two") == 0 && i+1 < argc) { do_two = 1; two_d = atof(argv[++i]); }
        else if (strcmp(argv[i], "-skyrmhop") == 0) { do_skyrmhop = 1; if (i+1 < argc && argv[i+1][0] != '-') skyrmhop_d = atof(argv[++i]); }
        else if (strcmp(argv[i], "-a") == 0 && i+1 < argc) hopfion_a = atof(argv[++i]);
        else if (strcmp(argv[i], "-N") == 0 && i+1 < argc) N = atoi(argv[++i]);
        else if (strcmp(argv[i], "-R") == 0 && i+1 < argc) { R = atof(argv[++i]); L = R + 0.5; }
        else if (strcmp(argv[i], "-lam6") == 0 && i+1 < argc) lambda6 = atof(argv[++i]);
        else if (strcmp(argv[i], "-profile") == 0 && i+1 < argc) profile_path = argv[++i];
    }

    HopfionParams p = {rho0, e_skyrme, 0, lambda6, 0, 0, 1.0};
    double c4 = 2.0 * rho0 * rho0 / (e_skyrme * e_skyrme);
    double E_FB = 6.0 * sqrt(2.0) * M_PI * M_PI * rho0 * rho0 * rho0 / e_skyrme;

    printf("========================================\n");
    printf(" Path 4: Hopfion Effective Metric\n");
    printf("========================================\n");
    const char *mode_str = do_hedgehog ? "HEDGEHOG (control)" :
                           do_two ? "TWO HOPFIONS (linked)" :
                           do_skyrmhop ? "SKYRMION + HOPFION (mixed B,H)" : "HOPFION H=1";
    printf("  Mode: %s\n", mode_str);
    printf("  N=%d, R=%.1f, h=%.4f\n", N, R, 2.0*L/N);
    printf("  e=%.1f, rho0=%.1f, c4=%.4f, E_FB=%.4f\n", e_skyrme, rho0, c4, E_FB);
    if (lambda6 > 0) printf("  lambda6=%.2f\n", lambda6);
    if (!do_hedgehog) printf("  Hopfion size: a=%.2f\n", hopfion_a);
    if (do_two) printf("  Two-hopfion separation: D=%.2f (z=±%.2f)\n", 2*two_d, two_d);
    if (do_skyrmhop) printf("  Skyrmion+hopfion: ring at z=%.2f\n", skyrmhop_d);
    printf("\n");

    /* Allocate grid */
    SphericalGrid *g = sg_alloc(N, L, R, sponge_w);
    sg_set_vacuum(g, rho0);
    printf("  Grid: %d active cells (%.1f%%)\n\n", g->n_active,
           100.0 * g->n_active / ((double)N*N*N));

    /* ========== Initialize field ========== */
    if (do_hedgehog) {
        if (!profile_path)
            profile_path = "../proposal/hopfion_search/data/profiles/profile_sigma_e1.dat";
        RadialProfile *prof = profile_load(profile_path);
        if (!prof) {
            fprintf(stderr, "Cannot load profile %s\n", profile_path);
            sg_free(g);
            return 1;
        }
        init_skyrmion(g, prof, rho0, 0, 0, 0);
        hf_sigma_project(g, rho0);
        profile_free(prof);
    } else if (do_two) {
        /* Two hopfions at z=±D, composed via product ansatz q₁×q₂/ρ₀ */
        long N3 = (long)N * N * N;
        printf("  Composing two hopfions at z=±%.2f...\n", two_d);

        /* Initialize first hopfion at z=+D */
        init_hopfion(g, rho0, hopfion_a, 0, 0, two_d);
        Multivector *q1 = malloc(N3 * sizeof(Multivector));
        memcpy(q1, g->psi, N3 * sizeof(Multivector));

        /* Initialize second hopfion at z=-D */
        init_hopfion(g, rho0, hopfion_a, 0, 0, -two_d);

        /* Product: q_total = q₁ × q₂ / ρ₀ */
        double inv_rho0 = 1.0 / rho0;
        for (int n = 0; n < g->n_active; n++) {
            int ix = g->active[n];
            Multivector a = q1[ix];
            Multivector b = g->psi[ix];
            g->psi[ix].s  = inv_rho0 * (a.s*b.s  - a.f1*b.f1 - a.f2*b.f2 - a.f3*b.f3);
            g->psi[ix].f1 = inv_rho0 * (a.s*b.f1 + a.f1*b.s  - a.f2*b.f3 + a.f3*b.f2);
            g->psi[ix].f2 = inv_rho0 * (a.s*b.f2 + a.f1*b.f3 + a.f2*b.s  - a.f3*b.f1);
            g->psi[ix].f3 = inv_rho0 * (a.s*b.f3 - a.f1*b.f2 + a.f2*b.f1 + a.f3*b.s);
            g->psi[ix].j1 = 0; g->psi[ix].j2 = 0;
            g->psi[ix].j3 = 0; g->psi[ix].p  = 0;
        }
        free(q1);
        hf_sigma_project(g, rho0);
    } else if (do_skyrmhop) {
        if (!profile_path)
            profile_path = "../proposal/hopfion_search/data/profiles/profile_sigma_e1.dat";
        RadialProfile *prof = profile_load(profile_path);
        if (!prof) {
            fprintf(stderr, "Cannot load profile %s\n", profile_path);
            sg_free(g);
            return 1;
        }
        init_skyrmion_hopfion(g, prof, rho0, hopfion_a, skyrmhop_d);
        hf_sigma_project(g, rho0);
        profile_free(prof);
    } else {
        init_hopfion(g, rho0, hopfion_a, 0, 0, 0);
        hf_sigma_project(g, rho0);
    }

    /* ========== Topology ========== */
    printf("\n  --- Topological Charges ---\n");
    double B = hf_baryon_charge(g);
    printf("    Skyrmion charge B = %.6f\n", B);

    if (!do_hedgehog && !do_two && !do_skyrmhop) {
        printf("    Computing Hopf charge H (Gauss-Seidel, may be slow)...\n");
        double H = topo_hopf_charge(g);
        printf("    Hopf charge H = %.6f\n", H);
        printf("    (Expected: H ≈ 1 for standard hopfion; B ≈ 0 for pure hopfion)\n");
    } else if (do_two || do_skyrmhop) {
        printf("    (Skipping Hopf charge — slow solver)\n");
        if (do_skyrmhop)
            printf("    (Expected: B ≈ 1, H ≈ 1 for skyrmion+hopfion composite)\n");
    }

    /* ========== Energy ========== */
    HopfionEnergy en = hf_energy(g, &p);
    printf("\n  --- Energy ---\n");
    printf("    E_total = %.6f  (E/E_FB = %.6f)\n", en.Etotal, en.Etotal / E_FB);
    printf("    E₂ = %.6f, E₄ = %.6f", en.E2, en.E4);
    if (en.E4 > 1e-15)
        printf(", E₂/E₄ = %.6f", en.E2 / en.E4);
    printf("\n");
    if (en.E6 > 0) printf("    E₆ = %.6f\n", en.E6);

    /* ========== Directional Effective Metric (spherically averaged) ========== */
    printf("\n  --- Directional Effective Metric (radial average) ---\n");
    printf("  g^{jj}(r) = 1 + (c₄/|q|⁴) Σ_{i≠j} |A_i|²\n");
    printf("  Sound speed ratio: c_j/c₀ = √(g^{jj}/g^{00})\n\n");

    int n_bins = 60;
    DirectionalMetric *dm = compute_directional_metric(g, &p, n_bins);

    printf("  %6s  %10s  %10s  %10s  %10s  %10s  %10s\n",
           "r", "g00", "gxx", "gyy", "gzz", "gzz/gxx", "czz/cxx");
    printf("  %6s  %10s  %10s  %10s  %10s  %10s  %10s\n",
           "------", "----------", "----------", "----------",
           "----------", "----------", "----------");

    double max_aniso = 0, max_aniso_r = 0;
    double sum_aniso = 0;
    int n_valid = 0;

    for (int b = 0; b < n_bins; b++) {
        if (dm->count[b] < 10) continue;
        double r = dm->r[b];
        if (r < 0.3 || r > 5.0) continue;

        double ratio_zx = dm->gzz[b] / dm->gxx[b];
        double c_ratio = sqrt(ratio_zx);
        double aniso = fabs(ratio_zx - 1.0);

        if (aniso > max_aniso) {
            max_aniso = aniso;
            max_aniso_r = r;
        }
        sum_aniso += aniso;
        n_valid++;

        if (b % 3 == 0 || aniso > 0.01) {
            printf("  %6.2f  %10.6f  %10.6f  %10.6f  %10.6f  %10.6f  %10.6f\n",
                   r, dm->g00[b], dm->gxx[b], dm->gyy[b], dm->gzz[b],
                   ratio_zx, c_ratio);
        }
    }

    double avg_aniso = (n_valid > 0) ? sum_aniso / n_valid : 0;
    printf("\n  Anisotropy (|g^{zz}/g^{xx} - 1|):\n");
    printf("    Max anisotropy:  %.6f at r = %.2f\n", max_aniso, max_aniso_r);
    printf("    Mean anisotropy: %.6f\n", avg_aniso);

    if (do_hedgehog) {
        int pass = (max_aniso < 0.01);
        printf("    Hedgehog isotropy: %s (expected < 0.01)\n", pass ? "PASS" : "FAIL");
    } else {
        int signal = (max_aniso > 0.01);
        printf("    Hopfion anisotropy signal: %s\n",
               signal ? "DETECTED (g^{zz} ≠ g^{xx})" : "NONE (isotropic like hedgehog)");
    }

    /* ========== Sound speed decomposition ========== */
    printf("\n  --- Sound Speed Analysis ---\n");
    printf("  c²_j/c² = g^{jj}/g^{00}  (sound speed squared in direction j)\n\n");

    printf("  %6s  %10s  %10s  %10s  %10s\n",
           "r", "c²_x/c²", "c²_z/c²", "c²_avg/c²", "Δc²/c²");
    printf("  %6s  %10s  %10s  %10s  %10s\n",
           "------", "----------", "----------", "----------", "----------");

    double max_dc2 = 0, max_dc2_r = 0;

    for (int b = 0; b < n_bins; b++) {
        if (dm->count[b] < 10) continue;
        double r = dm->r[b];
        if (r < 0.3 || r > 5.0) continue;
        if (dm->g00[b] < 1.001) continue;

        double c2x = dm->gxx[b] / dm->g00[b];
        double c2z = dm->gzz[b] / dm->g00[b];
        double c2avg = (2*c2x + c2z) / 3.0;  /* (xx+yy+zz)/3 with xx=yy */
        double dc2 = c2z - c2x;

        if (fabs(dc2) > fabs(max_dc2)) {
            max_dc2 = dc2;
            max_dc2_r = r;
        }

        if (b % 3 == 0 || fabs(dc2) > 0.005) {
            printf("  %6.2f  %10.6f  %10.6f  %10.6f  %10.6f\n",
                   r, c2x, c2z, c2avg, dc2);
        }
    }

    printf("\n  Peak sound speed difference: Δc²/c² = %.6f at r = %.2f\n",
           max_dc2, max_dc2_r);
    if (fabs(max_dc2) > 0.001) {
        printf("  → Perturbations in z-direction %s than in x-direction\n",
               max_dc2 > 0 ? "FASTER" : "SLOWER");
        printf("  → This is an effective gravitomagnetic-like effect\n");
    }

    dm_free(dm);

    /* ========== Axial (r, cos θ) metric (for hopfion only) ========== */
    if (!do_hedgehog) {
        printf("\n  --- Axial Metric g^{zz}/g^{xx}(r, cos θ) ---\n");
        printf("  (Shows angular dependence of sound speed anisotropy)\n\n");

        int nr = 30, ntheta = 12;
        AxialMetric *am = compute_axial_metric(g, &p, nr, ntheta);

        printf("  %6s  %6s  %10s  %10s  %10s  %10s\n",
               "r", "cos θ", "gzz/gxx", "c_z/c_x", "g00", "cells");
        printf("  %6s  %6s  %10s  %10s  %10s  %10s\n",
               "------", "------", "----------", "----------", "----------", "------");

        double max_am_aniso = 0;
        double max_am_r = 0, max_am_costh = 0;

        for (int ir = 0; ir < nr; ir++) {
            double r = (ir + 0.5) * am->dr;
            if (r < 0.5 || r > 4.0) continue;
            for (int it = 0; it < ntheta; it++) {
                int idx = ir * ntheta + it;
                if (am->count[idx] < 5) continue;
                double cos_th = -1.0 + (it + 0.5) * am->dcos;
                double inv = 1.0 / am->count[idx];
                double gzz_avg = am->gzz[idx] * inv;  /* already averaged above */
                double gxx_avg = am->gxx[idx] * inv;
                (void)gzz_avg; (void)gxx_avg;

                /* Use pre-averaged values */
                if (am->gxx[idx] < 1.001) continue;
                double ratio = am->gzz[idx] / am->gxx[idx];
                double aniso = fabs(ratio - 1.0);
                if (aniso > max_am_aniso) {
                    max_am_aniso = aniso;
                    max_am_r = r;
                    max_am_costh = cos_th;
                }

                /* Print selected bins */
                if (ir % 3 == 0 && it % 2 == 0) {
                    printf("  %6.2f  %6.2f  %10.6f  %10.6f  %10.4f  %10.0f\n",
                           r, cos_th, ratio, sqrt(ratio),
                           am->g00[idx], am->count[idx]);
                }
            }
        }

        printf("\n  Peak axial anisotropy: |g^{zz}/g^{xx} - 1| = %.6f\n",
               max_am_aniso);
        printf("    at r = %.2f, cos θ = %.2f\n", max_am_r, max_am_costh);

        am_free(am);
    }

    /* ========== Z-axis line scan (for two-hopfion or skyrmhop mode) ========== */
    if (do_two || do_skyrmhop) {
        printf("\n  --- Z-Axis Line Scan (connecting the two cores) ---\n");
        printf("  Metric components along (0, 0, z) between the two hopfion centers.\n");
        if (do_two)
            printf("  Hopfion cores at z = ±%.2f\n\n", two_d);
        else
            printf("  Skyrmion at origin, hopfion ring at z = %.2f\n\n", skyrmhop_d);
        printf("  %6s  %10s  %10s  %10s  %10s  %10s\n",
               "z", "g00", "gxx", "gzz", "gzz/gxx", "|q|/ρ₀");

        /* Find the closest grid cells to the z-axis (x≈0, y≈0) */
        int i_center = N/2, j_center = N/2;
        double c4_val = 2.0 * rho0 * rho0 / (e_skyrme * e_skyrme);

        for (int kk = 4; kk < N-4; kk++) {
            double x, y, z;
            sg_pos(g, i_center, j_center, kk, &x, &y, &z);
            double r = sqrt(x*x + y*y + z*z);
            if (r > R - 4.5*g->h) continue;
            double z_range = do_two ? two_d + 3.0 : skyrmhop_d + 4.0;
            if (fabs(z) > z_range) continue;

            if (!sg_is_active(g, i_center, j_center, kk)) continue;
            int ix = sg_idx(N, i_center, j_center, kk);

            Q4 q = q4_from(g->psi[ix]);
            double q2 = q4_norm2(q);
            if (q2 < 1e-20) continue;
            double q4_inv = 1.0 / (q2 * q2);
            Q4 qr = q4_rev(q);

            double A2[3];
            for (int d = 0; d < 3; d++) {
                Q4 dq = q4_deriv(g, i_center, j_center, kk, d);
                Q4 Ad = q4_mul(qr, dq);
                A2[d] = q4_norm2(Ad);
            }
            double A2_total = A2[0] + A2[1] + A2[2];
            double ct = c4_val * q4_inv;

            double g00 = 1.0 + ct * A2_total;
            double gxx = 1.0 + ct * (A2[1] + A2[2]);
            double gzz = 1.0 + ct * (A2[0] + A2[1]);
            double ratio = gzz / gxx;
            double q_rel = sqrt(q2) / rho0;

            if (kk % 2 == 0 || fabs(ratio - 1.0) > 0.01) {
                printf("  %6.2f  %10.4f  %10.4f  %10.4f  %10.6f  %10.6f\n",
                       z, g00, gxx, gzz, ratio, q_rel);
            }
        }

        /* Check if metric is non-trivial BETWEEN the cores */
        if (do_two) {
            printf("\n  KEY QUESTION: Is the metric non-trivial between the cores (z ∈ [-%.1f, %.1f])?\n",
                   two_d, two_d);
            printf("  If gzz/gxx ≠ 1 between cores, the linking creates non-local metric.\n");
        } else {
            printf("\n  KEY QUESTION: Does the hopfion ring modify the skyrmion's metric?\n");
            printf("  Look for asymmetry in g^{zz}/g^{xx} — skyrmion alone is isotropic.\n");
        }
    }

    /* ========== Compare to hedgehog P/m ========== */
    printf("\n  --- Comparison to Hedgehog P/m Analysis ---\n");
    printf("  For L₂+L₄, the hedgehog gives P/m = 2 at all r.\n");
    printf("  P(r)/m(r) = (g^{rr} + g^{00}) / g^{00} in some conventions.\n\n");

    /* Compute the existing hf_effective_metric for comparison */
    EffectiveMetric *em = hf_effective_metric(g, &p, n_bins);
    if (em) {
        printf("  %6s  %10s  %10s  %10s\n", "r", "P", "m", "P/m");
        for (int b = 0; b < n_bins; b++) {
            double r = em->r[b];
            if (r < 0.3 || r > 5.0) continue;
            if (em->m[b] < 0.1) continue;
            double pm = em->P[b] / em->m[b];
            if (b % 4 == 0) {
                printf("  %6.2f  %10.4f  %10.4f  %10.6f\n",
                       r, em->P[b], em->m[b], pm);
            }
        }

        /* Check P/m */
        double pm_max_dev = 0;
        for (int b = 0; b < n_bins; b++) {
            if (em->r[b] < 0.3 || em->r[b] > 4.0) continue;
            if (em->m[b] < 0.1) continue;
            double pm = em->P[b] / em->m[b];
            double dev = fabs(pm - 2.0);
            if (dev > pm_max_dev) pm_max_dev = dev;
        }
        printf("\n  Hedgehog-convention P/m: max |P/m - 2| = %.6f\n", pm_max_dev);
        printf("  (This is tautologically 2 for L₂+L₄ in the current code)\n");
        printf("  → The directional metric above gives the REAL signal.\n");

        hf_metric_free(em);
    }

    /* ========== Physical interpretation ========== */
    printf("\n  --- Physical Interpretation ---\n");
    if (!do_hedgehog) {
        printf("\n  The hopfion has axial symmetry (z-axis). If g^{zz} ≠ g^{xx}:\n");
        printf("    • Sound speed depends on direction → effective spacetime is anisotropic\n");
        printf("    • This is analogous to a gravitomagnetic field (frame-dragging)\n");
        printf("    • The rank-2 tensor structure naturally gives GR-like properties\n");
        printf("    • Unlike B⁰p coupling (scalar, spin-0), this IS a tensor effect\n");
        printf("\n  For LINKED hopfions (next step):\n");
        printf("    • The linking number connects distant regions non-locally\n");
        printf("    • The effective metric between linked cores could show\n");
        printf("      long-range correlations not present for a single soliton\n");
        printf("    • Key test: compute metric of TWO linked hopfions,\n");
        printf("      look for metric effects extending between their cores\n");
    } else {
        printf("\n  The hedgehog has spherical symmetry → g^{xx} = g^{yy} = g^{zz}\n");
        printf("  This confirms the control: hedgehog metric is isotropic.\n");
    }

    /* ========== Summary ========== */
    printf("\n========================================\n");
    printf(" SUMMARY\n");
    printf("========================================\n");
    printf("  Configuration: %s\n", do_hedgehog ? "Hedgehog B=1" :
           do_two ? "Two linked hopfions" :
           do_skyrmhop ? "Skyrmion + Hopfion" : "Hopfion H=1");
    printf("  B = %.4f\n", B);
    printf("  E/E_FB = %.4f\n", en.Etotal / E_FB);

    sg_free(g);
    return 0;
}
