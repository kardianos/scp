/*
 * verify3d.c — 3D verification of radial Skyrmion profiles
 *
 * Reads a 1D radial profile f(r) from file, initializes a 3D hedgehog grid,
 * computes the 3D energy E₂, E₄, Q, and compares against the 1D values.
 * Also computes the gradient (force) to verify the configuration is
 * near a stationary point.
 *
 * Usage: ./verify3d [-N 64] [-L 6.0] [-e 4.0] [-profile profile_B1.dat]
 *                   [-B 1] [-relax 0] [-steps 100]
 *
 * For B>1, uses rational map ansatz angular functions with the radial profile.
 */

#define _USE_MATH_DEFINES
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>
#include "field.h"
#include "initial.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ---------- Profile table ---------- */

typedef struct {
    int n;          /* number of points */
    double *r;      /* radial positions */
    double *f;      /* profile f(r) */
    double *fp;     /* derivative f'(r) */
    double rmax;    /* maximum r in table */
    double dr;      /* uniform spacing (if uniform) */
} Profile;

static Profile *profile_load(const char *filename)
{
    FILE *fp = fopen(filename, "r");
    if (!fp) { fprintf(stderr, "Cannot open %s\n", filename); return NULL; }

    /* Count data lines (skip comments starting with #) */
    int capacity = 25000;
    Profile *p = malloc(sizeof(Profile));
    p->r  = malloc(capacity * sizeof(double));
    p->f  = malloc(capacity * sizeof(double));
    p->fp = malloc(capacity * sizeof(double));
    p->n  = 0;

    char line[512];
    while (fgets(line, sizeof(line), fp)) {
        if (line[0] == '#' || line[0] == '\n') continue;
        double rv, fv, fpv;
        if (sscanf(line, "%lf %lf %lf", &rv, &fv, &fpv) >= 3) {
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
    }
    fclose(fp);

    if (p->n < 2) { fprintf(stderr, "Profile has too few points\n"); free(p); return NULL; }
    p->rmax = p->r[p->n - 1];
    p->dr = p->r[1] - p->r[0];  /* assume uniform spacing */

    printf("Loaded profile: %d points, r=[0, %.3f], dr=%.6f\n", p->n, p->rmax, p->dr);
    return p;
}

static void profile_free(Profile *p)
{
    if (p) { free(p->r); free(p->f); free(p->fp); free(p); }
}

/* Linear interpolation on the profile table */
static double profile_eval(const Profile *p, double r)
{
    if (r <= 0) return p->f[0];
    if (r >= p->rmax) return 0.0;  /* f → 0 at infinity */

    /* Fast lookup for uniform grid */
    double idx_f = r / p->dr;
    int i = (int)idx_f;
    if (i >= p->n - 1) i = p->n - 2;
    double t = idx_f - i;
    return (1.0 - t) * p->f[i] + t * p->f[i + 1];
}

/* ---------- B>1 rational map angular functions ----------
 *
 * For B=1: hedgehog, r_hat = (sin θ cos φ, sin θ sin φ, cos θ)
 * For B>1: rational map R(z) maps S² → S², giving a non-radial angular
 * structure. The unit vector n̂ on the target S² is:
 *   n̂ = (2 Re(R), 2 Im(R), 1-|R|²) / (1+|R|²)
 * where z = tan(θ/2) e^{iφ} is the stereographic coordinate.
 *
 * Standard rational maps:
 *   B=1: R(z) = z                    (hedgehog)
 *   B=2: R(z) = z²                   (toroidal)
 *   B=3: R(z) = (z³ - √3 iz)/(√3 iz² - 1)   (tetrahedral)
 *   B=4: R(z) = (z⁴ + 2√3 iz² + 1)/(z⁴ - 2√3 iz² + 1)  (cubic)
 */

/* Compute the target unit vector n̂ from rational map for given (θ, φ).
 * Returns n = (nx, ny, nz) with |n| = 1. */
static void rational_map_n(int B, double theta, double phi,
                           double *nx, double *ny, double *nz)
{
    /* Stereographic coordinate z = tan(θ/2) e^{iφ} */
    double t2 = tan(0.5 * theta);
    double zr = t2 * cos(phi);
    double zi = t2 * sin(phi);

    double Rr, Ri;  /* R(z) = Rr + i*Ri */

    if (B == 1) {
        /* R(z) = z */
        Rr = zr; Ri = zi;
    } else if (B == 2) {
        /* R(z) = z² */
        Rr = zr*zr - zi*zi;
        Ri = 2*zr*zi;
    } else if (B == 3) {
        /* R(z) = (z³ - √3 iz) / (√3 iz² - 1) */
        double s3 = sqrt(3.0);
        /* z² */
        double z2r = zr*zr - zi*zi;
        double z2i = 2*zr*zi;
        /* z³ */
        double z3r = z2r*zr - z2i*zi;
        double z3i = z2r*zi + z2i*zr;
        /* √3 iz = √3 (-zi + i zr) */
        double s3izr = -s3*zi;
        double s3izi =  s3*zr;
        /* √3 iz² = √3 * i * z² = √3 (-z2i + i z2r) */
        double s3iz2r = -s3*z2i;
        double s3iz2i =  s3*z2r;
        /* numerator: z³ - √3 iz */
        double nr = z3r - s3izr;
        double ni = z3i - s3izi;
        /* denominator: √3 iz² - 1 */
        double dr = s3iz2r - 1.0;
        double di = s3iz2i;
        /* complex division */
        double d2 = dr*dr + di*di;
        if (d2 < 1e-30) { Rr = 1e15; Ri = 0; }
        else { Rr = (nr*dr + ni*di)/d2; Ri = (ni*dr - nr*di)/d2; }
    } else if (B == 4) {
        /* R(z) = (z⁴ + 2√3 iz² + 1) / (z⁴ - 2√3 iz² + 1) */
        double s3 = sqrt(3.0);
        double z2r = zr*zr - zi*zi;
        double z2i = 2*zr*zi;
        double z4r = z2r*z2r - z2i*z2i;
        double z4i = 2*z2r*z2i;
        /* 2√3 iz² */
        double tir = -2*s3*z2i;
        double tii =  2*s3*z2r;
        double nr = z4r + tir + 1.0;
        double ni = z4i + tii;
        double dr = z4r - tir + 1.0;
        double di = z4i - tii;
        double d2 = dr*dr + di*di;
        if (d2 < 1e-30) { Rr = 1e15; Ri = 0; }
        else { Rr = (nr*dr + ni*di)/d2; Ri = (ni*dr - nr*di)/d2; }
    } else {
        /* Default: hedgehog */
        Rr = zr; Ri = zi;
    }

    /* n̂ = (2 Re(R), 2 Im(R), 1-|R|²) / (1+|R|²) */
    double R2 = Rr*Rr + Ri*Ri;
    double inv = 1.0 / (1.0 + R2);
    *nx = 2.0 * Rr * inv;
    *ny = 2.0 * Ri * inv;
    *nz = (1.0 - R2) * inv;
}

/* ---------- Initialize 3D hedgehog from profile ---------- */

static void init_from_profile(Field *fld, const Params *par,
                              const Profile *prof, int B)
{
    int N = fld->N;
    double rho0 = par->rho0;

    #pragma omp parallel for collapse(3) schedule(static)
    for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
    for (int k = 0; k < N; k++) {
        int ix = idx(N, i, j, k);

        /* Physical coordinates (cell-centered) */
        double x = -fld->L + (i + 0.5) * fld->h;
        double y = -fld->L + (j + 0.5) * fld->h;
        double z = -fld->L + (k + 0.5) * fld->h;
        double r = sqrt(x*x + y*y + z*z);

        /* Profile function from table */
        double fval = profile_eval(prof, r);
        double cos_f = cos(fval);
        double sin_f = sin(fval);

        /* Unit direction vector */
        double nx, ny, nz;
        if (r > 1e-12) {
            if (B == 1) {
                /* Hedgehog: n̂ = r̂ */
                double inv_r = 1.0 / r;
                nx = x * inv_r;
                ny = y * inv_r;
                nz = z * inv_r;
            } else {
                /* Rational map ansatz */
                double theta = acos(z / r);
                double phi = atan2(y, x);
                rational_map_n(B, theta, phi, &nx, &ny, &nz);
            }
        } else {
            /* At origin: f(0)=π, sin(π)=0, so bivector part vanishes */
            nx = 0; ny = 0; nz = 1;
        }

        /* q = ρ₀[cos(f) + sin(f)(n̂ · σ)]
         * Bivector basis: e23 ↔ x, e31 ↔ y, e12 ↔ z */
        fld->psi[ix].s  = rho0 * cos_f;
        fld->psi[ix].f1 = rho0 * sin_f * nx;
        fld->psi[ix].f2 = rho0 * sin_f * ny;
        fld->psi[ix].f3 = rho0 * sin_f * nz;

        /* Degenerate sector: zero */
        fld->psi[ix].j1 = 0;
        fld->psi[ix].j2 = 0;
        fld->psi[ix].j3 = 0;
        fld->psi[ix].p  = 0;
    }
}

/* ---------- Gradient flow step ---------- */

static double max_force_norm(const Multivector *force, int N3)
{
    double maxf = 0;
    #pragma omp parallel for reduction(max:maxf) schedule(static)
    for (int i = 0; i < N3; i++) {
        double f2 = mv_dot(force[i], force[i]);
        if (f2 > maxf) maxf = f2;
    }
    return sqrt(maxf);
}

static double rms_force_norm(const Multivector *force, int N3)
{
    double sum = 0;
    #pragma omp parallel for reduction(+:sum) schedule(static)
    for (int i = 0; i < N3; i++) {
        sum += mv_dot(force[i], force[i]);
    }
    return sqrt(sum / N3);
}

static void gradient_flow_step(Field *f, const Multivector *force, double dt)
{
    int N3 = f->N * f->N * f->N;
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < N3; i++) {
        f->psi[i] = mv_add(f->psi[i], mv_scale(dt, force[i]));
    }
}

static void project_sigma_model(Field *f, double rho0)
{
    int N3 = f->N * f->N * f->N;
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < N3; i++) {
        double s = f->psi[i].s, b1 = f->psi[i].f1;
        double b2 = f->psi[i].f2, b3 = f->psi[i].f3;
        double r = sqrt(s*s + b1*b1 + b2*b2 + b3*b3);
        if (r > 1e-30) {
            double scale = rho0 / r;
            f->psi[i].s  *= scale;
            f->psi[i].f1 *= scale;
            f->psi[i].f2 *= scale;
            f->psi[i].f3 *= scale;
        }
    }
}

/* ---------- Main ---------- */

int main(int argc, char **argv)
{
    /* Default parameters */
    int N = 64;
    double L = 6.0;
    int B = 1;
    int do_relax = 0;
    int relax_steps = 200;
    double dt = 1e-4;
    const char *profile_file = "profile_B1.dat";

    /* Physical parameters */
    Params params;
    params.rho0 = 1.0;
    params.lambda = 100.0;  /* irrelevant for sigma model */
    params.e_skyrme = 4.0;
    params.mu = 1.0;
    params.g_coupling = 0.0;
    params.c = 1.0;

    /* Parse command line */
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-N") == 0 && i+1 < argc) N = atoi(argv[++i]);
        else if (strcmp(argv[i], "-L") == 0 && i+1 < argc) L = atof(argv[++i]);
        else if (strcmp(argv[i], "-e") == 0 && i+1 < argc) params.e_skyrme = atof(argv[++i]);
        else if (strcmp(argv[i], "-rho0") == 0 && i+1 < argc) params.rho0 = atof(argv[++i]);
        else if (strcmp(argv[i], "-B") == 0 && i+1 < argc) B = atoi(argv[++i]);
        else if (strcmp(argv[i], "-profile") == 0 && i+1 < argc) profile_file = argv[++i];
        else if (strcmp(argv[i], "-relax") == 0 && i+1 < argc) do_relax = atoi(argv[++i]);
        else if (strcmp(argv[i], "-steps") == 0 && i+1 < argc) relax_steps = atoi(argv[++i]);
        else if (strcmp(argv[i], "-dt") == 0 && i+1 < argc) dt = atof(argv[++i]);
        else {
            fprintf(stderr, "Usage: %s [-N 64] [-L 6.0] [-e 4.0] [-B 1] [-profile file] [-relax 0] [-steps 200] [-dt 1e-4]\n", argv[0]);
            return 1;
        }
    }

    printf("=== 3D Verification of Radial Skyrmion Profile ===\n");
    printf("Grid: N=%d, L=%.2f, h=%.6f\n", N, L, 2.0*L/N);
    printf("Parameters: rho0=%.3f, e=%.3f, B=%d\n", params.rho0, params.e_skyrme, B);
    printf("Profile: %s\n", profile_file);
    printf("Threads: %d\n", omp_get_max_threads());
    printf("\n");

    /* Load profile */
    Profile *prof = profile_load(profile_file);
    if (!prof) return 1;

    /* Allocate field */
    Field *fld = field_alloc(N, L);
    int N3 = N * N * N;

    if (do_relax) {
        printf("Memory: field=%.1f MB, total≈%.1f MB (with gradient arrays)\n",
               N3 * sizeof(Multivector) / 1e6,
               N3 * (sizeof(Multivector) * 2 + 24*8 + 12*8) / 1e6);
    } else {
        printf("Memory: field=%.1f MB (init+energy only)\n",
               N3 * sizeof(Multivector) / 1e6);
    }
    printf("\n");

    /* Initialize from profile */
    printf("Initializing %dD grid from profile...\n", 3);
    init_from_profile(fld, &params, prof, B);

    /* Compute energy */
    Energy en = field_energy(fld, &params);
    double Q = field_topological_charge(fld, &params);

    /* Reference values (1D) */
    double c4 = 2.0 * params.rho0 * params.rho0 / (params.e_skyrme * params.e_skyrme);
    double E_FB = 6.0 * sqrt(2.0) * M_PI * M_PI * params.rho0 * params.rho0 * params.rho0 / params.e_skyrme;

    printf("--- Initial 3D Energy ---\n");
    printf("E_total = %14.8f\n", en.Etotal);
    printf("E2      = %14.8f\n", en.E2);
    printf("E4      = %14.8f\n", en.E4);
    printf("EV      = %14.8f\n", en.EV);
    printf("ED      = %14.8f\n", en.ED);
    printf("Q       = %14.8f  (expected: %d)\n", Q, B);
    printf("E2/E4   = %14.8f  (expected: 1.0 for sigma model)\n", en.E2 / en.E4);
    printf("E/E_FB  = %14.8f\n", en.Etotal / E_FB);
    printf("c4      = %.6f\n", c4);
    printf("E_FB    = %.6f\n", E_FB);
    printf("\n");

    /* Reference values for comparison */
    double E_ref_1D;
    if (B == 1) E_ref_1D = 1.2322 * E_FB;
    else if (B == 2) E_ref_1D = 1.2083 * 2 * E_FB;
    else if (B == 3) E_ref_1D = 1.1843 * 3 * E_FB;
    else if (B == 4) E_ref_1D = 1.1370 * 4 * E_FB;
    else E_ref_1D = 0;

    if (E_ref_1D > 0) {
        printf("--- Comparison with 1D ---\n");
        printf("E_1D    = %14.8f  (reference)\n", E_ref_1D);
        printf("E_3D    = %14.8f  (this computation)\n", en.Etotal);
        printf("Error   = %14.8f  (%.4f%%)\n", en.Etotal - E_ref_1D,
               100.0 * (en.Etotal - E_ref_1D) / E_ref_1D);
        printf("Q error = %14.8f\n", Q - B);
        printf("\n");
    }

    /* Force computation and optional gradient flow relaxation */
    Multivector *force = NULL;

    double max_f = 0;

    if (do_relax) {
        force = (Multivector *)malloc((size_t)N3 * sizeof(Multivector));
        if (!force) { fprintf(stderr, "malloc failed\n"); return 1; }

        field_gradient(fld, &params, force);
        max_f = max_force_norm(force, N3);
        double rms_f = rms_force_norm(force, N3);

        printf("--- Force (stationarity check) ---\n");
        printf("|F|_max = %14.8e\n", max_f);
        printf("|F|_rms = %14.8e\n", rms_f);
        printf("(Should be small if configuration is near energy minimum)\n");
        printf("\n");
    }

    if (do_relax) {
        printf("--- Gradient Flow Relaxation (%d steps, sigma model) ---\n", relax_steps);
        printf("%6s  %14s  %14s  %14s  %14s  %10s  %12s\n",
               "step", "E_total", "E2", "E4", "Q", "E2/E4", "|F|_max");

        /* For sigma model relaxation, set lambda=0.
         * The sigma model projection enforces |q|=rho0 exactly,
         * making the EV force irrelevant. Keeping lambda>0 creates
         * huge but useless forces that slow convergence. */
        Params sigma_params = params;
        sigma_params.lambda = 0.0;

        double E_prev = en.Etotal;

        for (int step = 1; step <= relax_steps; step++) {
            field_gradient(fld, &sigma_params, force);
            max_f = max_force_norm(force, N3);

            /* Adaptive step with energy backtracking */
            double actual_dt = dt;
            Multivector *psi_save = NULL;

            /* Save field for backtracking */
            psi_save = (Multivector *)malloc((size_t)N3 * sizeof(Multivector));
            memcpy(psi_save, fld->psi, (size_t)N3 * sizeof(Multivector));

            int accepted = 0;
            for (int tries = 0; tries < 10; tries++) {
                if (tries > 0) {
                    memcpy(fld->psi, psi_save, (size_t)N3 * sizeof(Multivector));
                    actual_dt *= 0.5;
                }
                gradient_flow_step(fld, force, actual_dt);
                project_sigma_model(fld, params.rho0);

                Energy en_new = field_energy(fld, &params);
                if (en_new.Etotal <= E_prev + 1e-15) {
                    en = en_new;
                    E_prev = en.Etotal;
                    dt = fmin(actual_dt * 1.1, 0.1);  /* grow dt */
                    accepted = 1;
                    break;
                }
            }
            free(psi_save);
            if (!accepted) continue;

            if (step % 10 == 0 || step == 1 || step == relax_steps) {
                Q = field_topological_charge(fld, &params);
                printf("%6d  %14.8f  %14.8f  %14.8f  %14.8f  %10.6f  %12.6e\n",
                       step, en.Etotal, en.E2, en.E4, Q, en.E2/en.E4, max_f);
            }
        }
        printf("\n");

        /* Final comparison */
        printf("--- Final 3D Energy (after relaxation) ---\n");
        printf("E_total = %14.8f\n", en.Etotal);
        printf("E2      = %14.8f\n", en.E2);
        printf("E4      = %14.8f\n", en.E4);
        printf("Q       = %14.8f\n", Q);
        printf("E2/E4   = %14.8f\n", en.E2 / en.E4);
        printf("E/E_FB  = %14.8f\n", en.Etotal / E_FB);
        if (E_ref_1D > 0) {
            printf("Error   = %.4f%%\n", 100.0 * (en.Etotal - E_ref_1D) / E_ref_1D);
        }
        printf("|F|_max = %14.8e\n", max_f);
    }

    /* Cleanup */
    if (force) free(force);
    field_free(fld);
    profile_free(prof);

    return 0;
}
