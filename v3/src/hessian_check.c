/*
 * hessian_check.c — Numerical verification of the L₄ Hessian
 *
 * Initializes a 3D B=1 hedgehog from a 1D radial profile, then:
 * 1. Computes the baseline force F(Ψ₀) (should be ~0 at equilibrium)
 * 2. For each test perturbation δΨ and several values of ε:
 *    - Computes F(Ψ₀ ± εδΨ)
 *    - Central-difference Hessian: H·δΨ = [F(+ε) - F(-ε)] / (2ε)
 *    - Checks convergence as ε → 0 (error should scale as ε²)
 * 3. Reports the converged Hessian structure and far-field behavior
 *
 * Test perturbations:
 *   (a) Breathing mode: δs = g(r), tests scalar Hessian
 *   (b) Bivector monopole: δf3 = g(r), tests bivector mass term
 *   (c) Bivector dipole: δf3 = g(r)·cos θ, tests angular coupling (EM-like)
 *
 * SUCCESS CRITERION #5: The numerical Hessian converges with O(ε²) error
 * for central differences (2nd order). The converged Hessian H·δΨ gives
 * the effective mass/potential for perturbations on the hedgehog background.
 *
 * PHYSICAL QUESTION: Does M²(r) → 0 as r → ∞? (Free Maxwell in far field)
 * Does it have structure in the core? What is the angular dependence?
 *
 * Usage: ./hessian_check [-N 64] [-L 6.0] [-e 1.0] [-profile PATH] [-outdir DIR]
 */

#define _USE_MATH_DEFINES
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>
#include "field.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ========== Profile I/O ========== */

typedef struct {
    int n;
    double *r, *f, *fp;
    double rmax, dr;
} Profile;

static Profile *profile_load(const char *filename)
{
    FILE *fp = fopen(filename, "r");
    if (!fp) { fprintf(stderr, "Cannot open %s\n", filename); return NULL; }

    int capacity = 25000;
    Profile *p = malloc(sizeof(Profile));
    p->r  = malloc(capacity * sizeof(double));
    p->f  = malloc(capacity * sizeof(double));
    p->fp = malloc(capacity * sizeof(double));
    p->n  = 0;

    char line[512];
    while (fgets(line, sizeof(line), fp)) {
        if (line[0] == '#' || line[0] == '\n') continue;
        double rv, fv, fpv = 0, d1, d2, d3;
        int nc = sscanf(line, "%lf %lf %lf %lf %lf %lf", &rv, &fv, &fpv, &d1, &d2, &d3);
        if (nc < 2) continue;
        if (p->n >= capacity) {
            capacity *= 2;
            p->r  = realloc(p->r,  capacity * sizeof(double));
            p->f  = realloc(p->f,  capacity * sizeof(double));
            p->fp = realloc(p->fp, capacity * sizeof(double));
        }
        p->r[p->n]  = rv;
        p->f[p->n]  = fv;
        p->fp[p->n] = fpv;
        p->n++;
    }
    fclose(fp);

    if (p->n < 2) { fprintf(stderr, "Profile too few points\n"); free(p); return NULL; }
    p->rmax = p->r[p->n - 1];
    p->dr = p->r[1] - p->r[0];
    return p;
}

static void profile_free(Profile *p) {
    if (p) { free(p->r); free(p->f); free(p->fp); free(p); }
}

static double profile_eval(const Profile *p, double r)
{
    if (r <= 0) return p->f[0];
    if (r >= p->rmax) return 0.0;
    double idx_f = r / p->dr;
    int i = (int)idx_f;
    if (i >= p->n - 1) i = p->n - 2;
    double t = idx_f - i;
    return (1.0 - t) * p->f[i] + t * p->f[i + 1];
}

/* ========== 3D hedgehog initialization ========== */

static void init_hedgehog(Field *fld, double rho0, const Profile *prof)
{
    int N = fld->N;

    #pragma omp parallel for collapse(3) schedule(static)
    for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
    for (int k = 0; k < N; k++) {
        int ix = idx(N, i, j, k);
        double x = -fld->L + (i + 0.5) * fld->h;
        double y = -fld->L + (j + 0.5) * fld->h;
        double z = -fld->L + (k + 0.5) * fld->h;
        double r = sqrt(x*x + y*y + z*z);

        double fval = profile_eval(prof, r);
        double cos_f = cos(fval);
        double sin_f = sin(fval);

        double nx, ny, nz;
        if (r > 1e-12) {
            double inv_r = 1.0 / r;
            nx = x * inv_r;
            ny = y * inv_r;
            nz = z * inv_r;
        } else {
            nx = 0; ny = 0; nz = 1;
        }

        fld->psi[ix].s  = rho0 * cos_f;
        fld->psi[ix].f1 = rho0 * sin_f * nx;
        fld->psi[ix].f2 = rho0 * sin_f * ny;
        fld->psi[ix].f3 = rho0 * sin_f * nz;
        fld->psi[ix].j1 = 0;
        fld->psi[ix].j2 = 0;
        fld->psi[ix].j3 = 0;
        fld->psi[ix].p  = 0;
    }
}

/* ========== Perturbation functions ========== */

/* Smooth radial envelope: g(r) = exp(-r²/σ²) */
static double radial_envelope(double r, double sigma)
{
    return exp(-r * r / (sigma * sigma));
}

/* Add perturbation to field: fld += epsilon * delta_psi
 * type 0: breathing (scalar): δs = g(r)
 * type 1: bivector monopole:  δf3 = g(r)
 * type 2: bivector dipole:    δf3 = g(r) * cos(θ) = g(r) * z/r
 * type 3: bivector quadrupole: δf3 = g(r) * (3cos²θ-1)/2
 */
static void add_perturbation(Field *fld, int type, double epsilon, double sigma)
{
    int N = fld->N;

    #pragma omp parallel for collapse(3) schedule(static)
    for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
    for (int k = 0; k < N; k++) {
        int ix = idx(N, i, j, k);
        double x = -fld->L + (i + 0.5) * fld->h;
        double y = -fld->L + (j + 0.5) * fld->h;
        double z = -fld->L + (k + 0.5) * fld->h;
        double r = sqrt(x*x + y*y + z*z);

        double g = radial_envelope(r, sigma);

        switch (type) {
        case 0: /* breathing: δs = ε g(r) */
            fld->psi[ix].s += epsilon * g;
            break;
        case 1: /* bivector monopole: δf3 = ε g(r) */
            fld->psi[ix].f3 += epsilon * g;
            break;
        case 2: /* bivector dipole: δf3 = ε g(r) cos θ */
            if (r > 1e-12)
                fld->psi[ix].f3 += epsilon * g * (z / r);
            break;
        case 3: /* bivector quadrupole: δf3 = ε g(r) P₂(cos θ) */
            if (r > 1e-12) {
                double ct = z / r;
                fld->psi[ix].f3 += epsilon * g * (1.5*ct*ct - 0.5);
            }
            break;
        }
    }
}

/* ========== Copy field ========== */

static void field_copy(const Field *src, Field *dst)
{
    int N3 = src->N * src->N * src->N;
    memcpy(dst->psi, src->psi, (size_t)N3 * sizeof(Multivector));
}

/* ========== Hessian analysis on a radial line ========== */

/* Sample the Hessian along the z-axis (θ=0) and x-axis (θ=π/2) */
static void sample_hessian_radial(const Multivector *hessian, const Field *fld,
                                   int component, /* 0=s, 1=f1, 2=f2, 3=f3 */
                                   double *r_out, double *hz_out, double *hx_out,
                                   int *n_out)
{
    int N = fld->N;
    int mid = N / 2;
    int count = 0;

    /* Along z-axis: (mid, mid, k) for k = mid..N-1 */
    for (int k = mid; k < N && count < N/2; k++) {
        int ix = idx(N, mid, mid, k);
        double z = -fld->L + (k + 0.5) * fld->h;
        double r = fabs(z);
        r_out[count] = r;

        switch (component) {
        case 0: hz_out[count] = hessian[ix].s;  break;
        case 1: hz_out[count] = hessian[ix].f1; break;
        case 2: hz_out[count] = hessian[ix].f2; break;
        case 3: hz_out[count] = hessian[ix].f3; break;
        }

        /* Along x-axis at same r: (mid+k-mid, mid, mid) */
        int ix2 = idx(N, k, mid, mid);
        switch (component) {
        case 0: hx_out[count] = hessian[ix2].s;  break;
        case 1: hx_out[count] = hessian[ix2].f1; break;
        case 2: hx_out[count] = hessian[ix2].f2; break;
        case 3: hx_out[count] = hessian[ix2].f3; break;
        }

        count++;
    }
    *n_out = count;
}

/* ========== Main ========== */

int main(int argc, char **argv)
{
    int N = 64;
    double L = 6.0;
    double rho0 = 1.0;
    double e_skyrme = 1.0;
    double sigma_pert = 1.5;  /* perturbation width */
    const char *profile_file = "data/profile_sigma_e1.dat";
    const char *output_dir = "data";

    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "-N") && i+1 < argc) N = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-L") && i+1 < argc) L = atof(argv[++i]);
        else if (!strcmp(argv[i], "-rho0") && i+1 < argc) rho0 = atof(argv[++i]);
        else if (!strcmp(argv[i], "-e") && i+1 < argc) e_skyrme = atof(argv[++i]);
        else if (!strcmp(argv[i], "-sigma") && i+1 < argc) sigma_pert = atof(argv[++i]);
        else if (!strcmp(argv[i], "-profile") && i+1 < argc) profile_file = argv[++i];
        else if (!strcmp(argv[i], "-outdir") && i+1 < argc) output_dir = argv[++i];
    }

    printf("============================================================\n");
    printf(" L₄ Hessian Check — Linearized Force on Hedgehog Background\n");
    printf("============================================================\n\n");

    printf("Grid: N=%d, L=%.2f, h=%.6f\n", N, L, 2.0*L/N);
    printf("Profile: %s\n", profile_file);
    printf("Parameters: ρ₀=%.4f, e=%.4f\n", rho0, e_skyrme);
    printf("Perturbation σ = %.4f\n\n", sigma_pert);

    /* Load profile */
    Profile *prof = profile_load(profile_file);
    if (!prof) return 1;
    printf("Profile loaded: %d points, R_max=%.2f\n\n", prof->n, prof->rmax);

    /* Set up parameters — sigma model (λ=0) for clean Hessian */
    Params par = {0};
    par.rho0 = rho0;
    par.e_skyrme = e_skyrme;
    par.lambda = 0.0;  /* sigma model: no potential, clean E₂+E₄ */
    par.mu = 0.0;
    par.m_pi_sq = 0.0;

    /* Allocate fields */
    Field *f0 = field_alloc(N, L);     /* background Ψ₀ */
    Field *fp = field_alloc(N, L);     /* Ψ₀ + ε δΨ */
    Field *fm = field_alloc(N, L);     /* Ψ₀ - ε δΨ */
    int N3 = N * N * N;
    Multivector *force0 = calloc(N3, sizeof(Multivector));
    Multivector *forcep = calloc(N3, sizeof(Multivector));
    Multivector *forcem = calloc(N3, sizeof(Multivector));
    Multivector *hessian = calloc(N3, sizeof(Multivector));
    Multivector *hessian_prev = calloc(N3, sizeof(Multivector));

    /* Initialize hedgehog */
    printf("--- Initializing 3D hedgehog ---\n");
    init_hedgehog(f0, rho0, prof);

    /* Compute baseline energy and force */
    Energy en = field_energy(f0, &par);
    printf("Energy: E₂=%.4f, E₄=%.4f, E_total=%.4f\n", en.E2, en.E4, en.Etotal);
    printf("E₂/E₄ = %.6f (virial: should be ~1.0)\n", en.E2/en.E4);

    double Q = field_topological_charge(f0, &par);
    printf("Topological charge: Q = %.6f\n", Q);

    field_gradient(f0, &par, force0);
    double max_f0 = 0, rms_f0 = 0;
    for (int i = 0; i < N3; i++) {
        double f2 = mv_dot(force0[i], force0[i]);
        if (f2 > max_f0) max_f0 = f2;
        rms_f0 += f2;
    }
    max_f0 = sqrt(max_f0);
    rms_f0 = sqrt(rms_f0 / N3);
    printf("Baseline force: max=%.6e, rms=%.6e (should be ~0 for equilibrium)\n\n",
           max_f0, rms_f0);

    /* ========== Hessian convergence test ========== */

    const char *pert_names[] = {
        "Breathing (δs)",
        "Bivector monopole (δf₃)",
        "Bivector dipole (δf₃·cosθ)",
        "Bivector quadrupole (δf₃·P₂)"
    };
    int n_pert_types = 4;

    /* Epsilon values for convergence test */
    double eps_vals[] = {0.1, 0.05, 0.025, 0.0125, 0.00625};
    int n_eps = 5;

    /* Allocate radial sampling arrays */
    int max_radial = N / 2;
    double *r_samp = malloc(max_radial * sizeof(double));
    double *hz_samp = malloc(max_radial * sizeof(double));
    double *hx_samp = malloc(max_radial * sizeof(double));

    for (int pt = 0; pt < n_pert_types; pt++) {
        printf("============================================================\n");
        printf(" Perturbation %d: %s\n", pt, pert_names[pt]);
        printf("============================================================\n\n");

        printf("  %12s  %14s  %14s  %14s  %10s\n",
               "epsilon", "||H||_max", "||H||_rms", "||ΔH||_rms", "order");

        double prev_delta_rms = 0;

        for (int ie = 0; ie < n_eps; ie++) {
            double eps = eps_vals[ie];

            /* Ψ₊ = Ψ₀ + ε δΨ */
            field_copy(f0, fp);
            add_perturbation(fp, pt, +eps, sigma_pert);

            /* Ψ₋ = Ψ₀ - ε δΨ */
            field_copy(f0, fm);
            add_perturbation(fm, pt, -eps, sigma_pert);

            /* Compute forces */
            field_gradient(fp, &par, forcep);
            field_gradient(fm, &par, forcem);

            /* Central-difference Hessian: H = [F(+ε) - F(-ε)] / (2ε) */
            double max_h = 0, rms_h = 0;
            double delta_rms = 0;
            double inv2eps = 1.0 / (2.0 * eps);

            for (int i = 0; i < N3; i++) {
                hessian[i].s  = (forcep[i].s  - forcem[i].s)  * inv2eps;
                hessian[i].f1 = (forcep[i].f1 - forcem[i].f1) * inv2eps;
                hessian[i].f2 = (forcep[i].f2 - forcem[i].f2) * inv2eps;
                hessian[i].f3 = (forcep[i].f3 - forcem[i].f3) * inv2eps;
                hessian[i].j1 = 0;
                hessian[i].j2 = 0;
                hessian[i].j3 = 0;
                hessian[i].p  = 0;

                double h2 = hessian[i].s*hessian[i].s + hessian[i].f1*hessian[i].f1
                           + hessian[i].f2*hessian[i].f2 + hessian[i].f3*hessian[i].f3;
                if (h2 > max_h) max_h = h2;
                rms_h += h2;

                if (ie > 0) {
                    double ds = hessian[i].s - hessian_prev[i].s;
                    double d1 = hessian[i].f1 - hessian_prev[i].f1;
                    double d2 = hessian[i].f2 - hessian_prev[i].f2;
                    double d3 = hessian[i].f3 - hessian_prev[i].f3;
                    delta_rms += ds*ds + d1*d1 + d2*d2 + d3*d3;
                }
            }
            max_h = sqrt(max_h);
            rms_h = sqrt(rms_h / N3);
            delta_rms = sqrt(delta_rms / N3);

            /* Convergence order: δH ~ ε^p → p = log(δH_new/δH_old) / log(ε_new/ε_old)
             * For central diff, expect p = 2 */
            double order = 0;
            if (ie > 1 && prev_delta_rms > 1e-30 && delta_rms > 1e-30)
                order = log(delta_rms / prev_delta_rms) / log(eps / eps_vals[ie-1]);

            printf("  %12.6f  %14.6e  %14.6e  %14.6e  %10.2f\n",
                   eps, max_h, rms_h, delta_rms, order);

            /* Save for next iteration */
            memcpy(hessian_prev, hessian, N3 * sizeof(Multivector));
            prev_delta_rms = delta_rms;
        }

        /* ===== Radial profile of converged Hessian (last ε) ===== */
        printf("\n  Converged Hessian radial profile (ε=%.6f):\n", eps_vals[n_eps-1]);

        /* For each component, sample along z-axis and x-axis */
        for (int comp = 0; comp < 4; comp++) {
            const char *comp_names[] = {"s", "f1", "f2", "f3"};
            int n_samp;
            sample_hessian_radial(hessian, f0, comp, r_samp, hz_samp, hx_samp, &n_samp);

            /* Check if this component has significant signal */
            double max_sig = 0;
            for (int i = 0; i < n_samp; i++) {
                if (fabs(hz_samp[i]) > max_sig) max_sig = fabs(hz_samp[i]);
                if (fabs(hx_samp[i]) > max_sig) max_sig = fabs(hx_samp[i]);
            }
            if (max_sig < 1e-10) continue;

            printf("\n  Component %s (max=%.4e):\n", comp_names[comp], max_sig);
            printf("  %8s  %14s  %14s  %14s\n", "r", "H_z(r)", "H_x(r)", "anisotropy");
            for (int i = 0; i < n_samp; i += (n_samp > 20 ? n_samp/20 : 1)) {
                double aniso = 0;
                if (fabs(hz_samp[i]) > 1e-15 && fabs(hx_samp[i]) > 1e-15)
                    aniso = hz_samp[i] / hx_samp[i];
                printf("  %8.4f  %14.6e  %14.6e  %14.4f\n",
                       r_samp[i], hz_samp[i], hx_samp[i], aniso);
            }
        }

        /* ===== Far-field check ===== */
        printf("\n  Far-field check (r > 3): does H → 0?\n");
        int mid = N/2;
        double max_far = 0;
        int count_far = 0;
        for (int k = 0; k < N; k++) {
            double z = -f0->L + (k + 0.5) * f0->h;
            double r = fabs(z);
            if (r < 3.0 || r > L - 1.0) continue;
            int ix = idx(N, mid, mid, k);
            double h2 = hessian[ix].s*hessian[ix].s + hessian[ix].f1*hessian[ix].f1
                       + hessian[ix].f2*hessian[ix].f2 + hessian[ix].f3*hessian[ix].f3;
            if (h2 > max_far) max_far = h2;
            count_far++;
        }
        max_far = sqrt(max_far);
        printf("  max |H| for r > 3: %.6e (over %d points)\n", max_far, count_far);
        printf("  max |H| at core:   %.6e\n",
               sqrt(hessian[idx(N,mid,mid,mid)].s*hessian[idx(N,mid,mid,mid)].s +
                    hessian[idx(N,mid,mid,mid)].f1*hessian[idx(N,mid,mid,mid)].f1 +
                    hessian[idx(N,mid,mid,mid)].f2*hessian[idx(N,mid,mid,mid)].f2 +
                    hessian[idx(N,mid,mid,mid)].f3*hessian[idx(N,mid,mid,mid)].f3));
        printf("  Ratio (far/core): %.6e\n\n", max_far /
               (1e-30 + sqrt(hessian[idx(N,mid,mid,mid)].s*hessian[idx(N,mid,mid,mid)].s +
                             hessian[idx(N,mid,mid,mid)].f1*hessian[idx(N,mid,mid,mid)].f1 +
                             hessian[idx(N,mid,mid,mid)].f2*hessian[idx(N,mid,mid,mid)].f2 +
                             hessian[idx(N,mid,mid,mid)].f3*hessian[idx(N,mid,mid,mid)].f3)));

        /* ===== Write radial profile to file ===== */
        char fname[512];
        snprintf(fname, sizeof(fname), "%s/hessian_pert%d.dat", output_dir, pt);
        FILE *out = fopen(fname, "w");
        if (out) {
            fprintf(out, "# Hessian radial profile for perturbation: %s\n", pert_names[pt]);
            fprintf(out, "# N=%d, L=%.2f, e=%.4f, sigma=%.4f, eps=%.6f\n",
                    N, L, e_skyrme, sigma_pert, eps_vals[n_eps-1]);
            fprintf(out, "# Columns: r  Hz_s  Hz_f1  Hz_f2  Hz_f3  Hx_s  Hx_f1  Hx_f2  Hx_f3\n");

            /* Sample along z-axis and x-axis */
            for (int k = mid; k < N; k++) {
                double z = -f0->L + (k + 0.5) * f0->h;
                double r = fabs(z);
                int iz = idx(N, mid, mid, k);
                int ix = idx(N, k, mid, mid);

                fprintf(out, "%.8e  %.10e  %.10e  %.10e  %.10e  %.10e  %.10e  %.10e  %.10e\n",
                        r,
                        hessian[iz].s, hessian[iz].f1, hessian[iz].f2, hessian[iz].f3,
                        hessian[ix].s, hessian[ix].f1, hessian[ix].f2, hessian[ix].f3);
            }
            fclose(out);
            printf("  Written: %s\n\n", fname);
        }
    }

    /* ========== Summary ========== */
    printf("============================================================\n");
    printf(" Summary\n");
    printf("============================================================\n\n");

    printf("Background:\n");
    printf("  E₂=%.4f, E₄=%.4f, Q=%.4f\n", en.E2, en.E4, Q);
    printf("  Baseline force: max=%.2e, rms=%.2e\n\n", max_f0, rms_f0);

    printf("Hessian convergence:\n");
    printf("  Central-difference Hessian should converge as O(ε²)\n");
    printf("  Convergence order ~2.0 confirms the force is smooth at equilibrium\n\n");

    printf("Physical interpretation:\n");
    printf("  H·δΨ = [linearized force at Ψ₀] × δΨ\n");
    printf("  For □δΨ + M²·δΨ = J_source:\n");
    printf("    M²(r) = -H(r)/δΨ(r) in the static limit\n");
    printf("    M²(r) → 0 as r → ∞ means free Maxwell in far field\n");
    printf("    M²(r) > 0 in core means perturbation is repelled (scattered)\n");
    printf("    M²(r) < 0 in core means perturbation is attracted (may bind)\n\n");

    printf("NULL RESULT CRITERION #2:\n");
    printf("  If H ≡ 0 for all bivector perturbations, the soliton does NOT\n");
    printf("  generate EM fields classically — sourced Maxwell requires quantum\n");
    printf("  collective coordinate quantization (as in standard ANW).\n");

    /* Cleanup */
    free(r_samp); free(hz_samp); free(hx_samp);
    free(force0); free(forcep); free(forcem);
    free(hessian); free(hessian_prev);
    field_free(f0); field_free(fp); field_free(fm);
    profile_free(prof);

    return 0;
}
