/*
 * scatter.c — Soliton-soliton scattering simulation
 *
 * Time-dependent 3D evolution of two B=1 Skyrmions colliding head-on.
 *
 * PHYSICS:
 * The equation of motion for the quaternion field q(x,t) is:
 *   (1/c²) ∂²q/∂t² = F(q)
 * where F = -δE/δq is the force (negative energy gradient).
 *
 * F includes contributions from:
 *   F_E2 = ∇²q                          (gradient/elastic force)
 *   F_E4 = Skyrme force                  (topological stabilization)
 *   F_EV = -λ(|q|²-ρ₀²)q               (potential force)
 *
 * METHOD:
 * Leapfrog (Störmer-Verlet) time integration:
 *   v(t+dt/2) = v(t-dt/2) + dt × c² × F(q(t))
 *   q(t+dt)   = q(t)      + dt × v(t+dt/2)
 *
 * This is symplectic (energy-preserving) and 2nd-order accurate.
 *
 * INITIAL CONDITIONS:
 * Two B=1 solitons, centered at ±z₀ on the z-axis, moving towards
 * each other at speed ±v₀. We use the PRODUCT ANSATZ:
 *   q(x) = q₁(x - z₀ẑ) · q₂(x + z₀ẑ) / ρ₀
 * where q₁, q₂ are hedgehog solitons loaded from profile_B1.dat.
 *
 * The initial velocity field is:
 *   v(x) = ∂q/∂t|_{t=0} = (-v₀)(∂q₁/∂z)·q₂/ρ₀ + q₁·(+v₀)(∂q₂/∂z)/ρ₀
 *
 * For an ANTISOLITON (B=-1): q_anti(x) = q̃(x) (reverse orientation).
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

/* ========== Profile loading ========== */

typedef struct {
    double *r;
    double *f;
    double *rho;    /* radial amplitude ρ(r); NULL for σ-model profile */
    int n;
    double dr;
    double r_max;
} RadialProfile;

static RadialProfile *load_radial_profile(const char *filename)
{
    FILE *fp = fopen(filename, "r");
    if (!fp) { fprintf(stderr, "Cannot open %s\n", filename); return NULL; }

    /* Check header for finite-lambda indicator (has rho(r) column) */
    int has_rho = 0;
    int n = 0;
    char line[1024];
    while (fgets(line, sizeof(line), fp)) {
        if (line[0] == '#') {
            if (strstr(line, "rho(r)") || strstr(line, "Finite-lambda"))
                has_rho = 1;
            continue;
        }
        double rv, fv;
        if (sscanf(line, "%lf %lf", &rv, &fv) >= 2) n++;
    }
    rewind(fp);

    RadialProfile *p = malloc(sizeof(RadialProfile));
    p->r = malloc(n * sizeof(double));
    p->f = malloc(n * sizeof(double));
    p->rho = has_rho ? malloc(n * sizeof(double)) : NULL;
    p->n = n;

    int i = 0;
    while (fgets(line, sizeof(line), fp) && i < n) {
        if (line[0] == '#') continue;
        double rv, fv, fpv, rhov;
        int ncols = sscanf(line, "%lf %lf %lf %lf", &rv, &fv, &fpv, &rhov);
        if (ncols >= 2) {
            p->r[i] = rv;
            p->f[i] = fv;
            if (has_rho && ncols >= 4) p->rho[i] = rhov;
            i++;
        }
    }
    fclose(fp);

    p->dr = (n > 1) ? p->r[1] - p->r[0] : 0.001;
    p->r_max = p->r[n-1];
    return p;
}

static void free_radial_profile(RadialProfile *p)
{
    if (p) { free(p->r); free(p->f); if (p->rho) free(p->rho); free(p); }
}

/* Interpolate profile at arbitrary r (linear) */
static double interp_f(const RadialProfile *p, double r)
{
    if (r < 0) r = -r;
    if (r >= p->r_max) return 0.0;  /* f→0 at infinity */
    double fi = r / p->dr;
    int i = (int)fi;
    if (i >= p->n - 1) return 0.0;
    double t = fi - i;
    return (1-t)*p->f[i] + t*p->f[i+1];
}

/* Interpolate ρ(r); returns rho0 if no ρ profile available */
static double interp_rho(const RadialProfile *p, double rho0, double r)
{
    if (!p->rho) return rho0;
    if (r < 0) r = -r;
    if (r >= p->r_max) return rho0;  /* ρ→ρ₀ at infinity */
    double fi = r / p->dr;
    int i = (int)fi;
    if (i >= p->n - 1) return rho0;
    double t = fi - i;
    return (1-t)*p->rho[i] + t*p->rho[i+1];
}

/* Compute hedgehog quaternion q(x) = ρ₀[cos(f(r)) + sin(f(r)) r̂·σ]
 * centered at (cx, cy, cz) */
typedef struct { double s, f1, f2, f3; } Quat;

static Quat hedgehog_q(const RadialProfile *prof, double rho0,
                       double x, double y, double z,
                       double cx, double cy, double cz)
{
    double dx = x - cx, dy = y - cy, dz = z - cz;
    double r = sqrt(dx*dx + dy*dy + dz*dz);
    double f = interp_f(prof, r);
    double rho = interp_rho(prof, rho0, r);
    double cf = cos(f), sf = sin(f);

    Quat q;
    q.s = rho * cf;
    if (r > 1e-12) {
        double sr = rho * sf / r;
        q.f1 = sr * dx;
        q.f2 = sr * dy;
        q.f3 = sr * dz;
    } else {
        q.f1 = q.f2 = q.f3 = 0;
    }
    return q;
}

/* Quaternion product */
static inline Quat qmul(Quat a, Quat b) {
    return (Quat){
        a.s*b.s  - a.f1*b.f1 - a.f2*b.f2 - a.f3*b.f3,
        a.s*b.f1 + a.f1*b.s  - a.f2*b.f3 + a.f3*b.f2,
        a.s*b.f2 + a.f1*b.f3 + a.f2*b.s  - a.f3*b.f1,
        a.s*b.f3 - a.f1*b.f2 + a.f2*b.f1 + a.f3*b.s
    };
}

static inline Quat qrev(Quat a) {
    return (Quat){a.s, -a.f1, -a.f2, -a.f3};
}

static inline Quat qscale(double c, Quat a) {
    return (Quat){c*a.s, c*a.f1, c*a.f2, c*a.f3};
}

static inline Quat qadd(Quat a, Quat b) {
    return (Quat){a.s+b.s, a.f1+b.f1, a.f2+b.f2, a.f3+b.f3};
}

static inline Quat qsub(Quat a, Quat b) {
    return (Quat){a.s-b.s, a.f1-b.f1, a.f2-b.f2, a.f3-b.f3};
}

/* ========== Diagnostics ========== */

/* Compute center-of-charge along z-axis:
 * z_cm = ∫ z ρ_B(x) d³x / Q
 * Uses the topological charge density. */
static double charge_center_z(const Field *f, const Params *p, int sign)
{
    int N = f->N;
    double sum_zrho = 0, sum_rho = 0;

    #pragma omp parallel for collapse(3) reduction(+:sum_zrho,sum_rho) schedule(static)
    for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
    for (int k = 0; k < N; k++) {
        double z = -f->L + (k + 0.5) * f->h;
        int ix = idx(N, i, j, k);

        /* Approximate charge density from |q| variation */
        /* Use the actual scalar part deviation as a proxy for position */
        double q_s = f->psi[ix].s;
        /* The soliton center has q.s ≈ -ρ₀ (anti-vacuum) */
        double weight = (p->rho0 - q_s);  /* large at soliton center */
        if (weight < 0) weight = 0;
        weight = weight * weight;  /* sharpen */

        if (sign > 0 && z > 0) {
            sum_zrho += z * weight;
            sum_rho += weight;
        } else if (sign < 0 && z < 0) {
            sum_zrho += z * weight;
            sum_rho += weight;
        } else if (sign == 0) {
            sum_zrho += z * weight;
            sum_rho += weight;
        }
    }

    return (sum_rho > 1e-20) ? sum_zrho / sum_rho : 0;
}

/* Find the minimum of s-component along the z-axis (soliton core indicator) */
static double find_soliton_z(const Field *f, double rho0, int which)
{
    int N = f->N;
    /* Scan along z at x=y=0 (center of grid) */
    int ic = N/2;
    double best_z = 0;
    double best_s = rho0;

    if (which == 1) {
        /* Find minimum s with z > 0 (positive-z soliton) */
        for (int k = N/2; k < N; k++) {
            double s = f->psi[idx(N, ic, ic, k)].s;
            if (s < best_s) {
                best_s = s;
                best_z = -f->L + (k + 0.5) * f->h;
            }
        }
    } else {
        /* Find minimum s with z < 0 (negative-z soliton) */
        for (int k = 0; k < N/2; k++) {
            double s = f->psi[idx(N, ic, ic, k)].s;
            if (s < best_s) {
                best_s = s;
                best_z = -f->L + (k + 0.5) * f->h;
            }
        }
    }
    return best_z;
}

/* Find soliton core position in 3D via charge-density weighted centroid.
 * Scans full grid, splits into z>0 and z<0 halves.
 * Returns position via pointers. */
static void find_soliton_3d(const Field *f, double rho0,
                             double *x1, double *y1, double *z1_out,
                             double *x2, double *y2, double *z2_out)
{
    int N = f->N;
    double h = f->h, L = f->L;
    /* Weighted centroid using (rho0 - s)^2 as weight */
    double sx1=0, sy1=0, sz1=0, sw1=0;
    double sx2=0, sy2=0, sz2=0, sw2=0;

    for (int i = 2; i < N-2; i++) {
        double x = -L + (i + 0.5) * h;
        for (int j = 2; j < N-2; j++) {
            double y = -L + (j + 0.5) * h;
            for (int k = 2; k < N-2; k++) {
                double z = -L + (k + 0.5) * h;
                double s = f->psi[idx(N, i, j, k)].s;
                double w = rho0 - s;
                if (w < 0) w = 0;
                w = w * w;  /* sharpen */
                if (z > 0) {
                    sx1 += x * w; sy1 += y * w; sz1 += z * w; sw1 += w;
                } else {
                    sx2 += x * w; sy2 += y * w; sz2 += z * w; sw2 += w;
                }
            }
        }
    }
    *x1 = sw1 > 1e-20 ? sx1/sw1 : 0;
    *y1 = sw1 > 1e-20 ? sy1/sw1 : 0;
    *z1_out = sw1 > 1e-20 ? sz1/sw1 : 0;
    *x2 = sw2 > 1e-20 ? sx2/sw2 : 0;
    *y2 = sw2 > 1e-20 ? sy2/sw2 : 0;
    *z2_out = sw2 > 1e-20 ? sz2/sw2 : 0;
}

/* ========== Main simulation ========== */

int main(int argc, char *argv[])
{
    /* Default parameters */
    double v0 = 0.3;       /* initial speed (fraction of c) */
    double z0 = 4.0;       /* initial separation / 2 */
    double rho0 = 1.0;
    double e_skyrme = 4.0;
    double lambda = 10000.0;   /* finite λ for stability (σ-model unstable on lattice) */
    int sigma_model = 0;       /* full model by default */
    double c_light = 1.0;
    int N = 128;
    double L = 8.0;
    double T_max = 20.0;    /* total simulation time */
    double dt = 0.0;        /* auto-set from CFL */
    int output_interval = 50;
    int antisoliton = 0;    /* 0 = B+B, 1 = B+B̄ */
    int isorotate = 0;      /* 0 = attractive channel, 1 = repulsive (π-rotate soliton 2) */
    int single_soliton = 0; /* 1 = single soliton at origin (stability test) */
    const char *profile_file = NULL;  /* auto-detect based on lambda */
    double damp_gamma = 0;    /* velocity damping rate (0=off) */
    double settle_time = 0;   /* settling phase duration (damping on) */
    double rdamp = 0;         /* selective radial damping rate: damps v·q̂ component only */
    int bc_clamp = 3;         /* boundary layers clamped to vacuum (0=periodic, >0=Dirichlet) */

    /* Parse command line */
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-v") == 0 && i+1 < argc) v0 = atof(argv[++i]);
        else if (strcmp(argv[i], "-z") == 0 && i+1 < argc) z0 = atof(argv[++i]);
        else if (strcmp(argv[i], "-N") == 0 && i+1 < argc) N = atoi(argv[++i]);
        else if (strcmp(argv[i], "-L") == 0 && i+1 < argc) L = atof(argv[++i]);
        else if (strcmp(argv[i], "-T") == 0 && i+1 < argc) T_max = atof(argv[++i]);
        else if (strcmp(argv[i], "-dt") == 0 && i+1 < argc) dt = atof(argv[++i]);
        else if (strcmp(argv[i], "-e") == 0 && i+1 < argc) e_skyrme = atof(argv[++i]);
        else if (strcmp(argv[i], "-lambda") == 0 && i+1 < argc) {
            lambda = atof(argv[++i]);
            sigma_model = (lambda == 0) ? 1 : 0;
        }
        else if (strcmp(argv[i], "-sigma") == 0) { sigma_model = 1; lambda = 0; }
        else if (strcmp(argv[i], "-anti") == 0) antisoliton = 1;
        else if (strcmp(argv[i], "-out") == 0 && i+1 < argc) output_interval = atoi(argv[++i]);
        else if (strcmp(argv[i], "-profile") == 0 && i+1 < argc) profile_file = argv[++i];
        else if (strcmp(argv[i], "-damp") == 0 && i+1 < argc) damp_gamma = atof(argv[++i]);
        else if (strcmp(argv[i], "-settle") == 0 && i+1 < argc) settle_time = atof(argv[++i]);
        else if (strcmp(argv[i], "-isorot") == 0) isorotate = 1;
        else if (strcmp(argv[i], "-single") == 0) single_soliton = 1;
        else if (strcmp(argv[i], "-rdamp") == 0 && i+1 < argc) rdamp = atof(argv[++i]);
        else if (strcmp(argv[i], "-bc_clamp") == 0 && i+1 < argc) bc_clamp = atoi(argv[++i]);
        else {
            fprintf(stderr, "Usage: %s [-v speed] [-z sep/2] [-N grid] [-L box] "
                    "[-T time] [-dt step] [-lambda lam] [-sigma] [-anti] [-out interval] "
                    "[-profile file] [-damp gamma] [-settle time]\n", argv[0]);
            return 1;
        }
    }

    setvbuf(stdout, NULL, _IOLBF, 0);  /* line-buffered stdout */

    double h = 2.0 * L / N;
    if (dt <= 0) {
        double dt_cfl = 0.3 * h / c_light;  /* CFL for wave equation */
        if (lambda > 0) {
            /* Verlet stability for potential: dt < 2/omega
             * omega = sqrt(8λ)ρ₀c for radial oscillations around |q|=ρ₀
             * Use safety factor 0.5 */
            double omega = c_light * sqrt(8.0*lambda) * rho0;
            double dt_pot = 1.0 / omega;  /* 0.5 × (2/ω) */
            if (dt_pot < dt_cfl) dt_cfl = dt_pot;
        }
        dt = dt_cfl;
    }

    printf("==========================================================\n");
    printf(" Soliton-Soliton Scattering Simulation\n");
    printf("==========================================================\n\n");
    printf("Grid: N=%d, L=%.1f, h=%.4f\n", N, L, h);
    const char *channel_str = single_soliton ? "single B=1 (stability test)" :
                              antisoliton ? "B=+1 + B=-1 (soliton-antisoliton)" :
                              isorotate  ? "B=+1 + B=+1 (repulsive channel)" :
                                           "B=+1 + B=+1 (attractive channel)";
    printf("Solitons: %s, separation=%.1f, v=%.3fc\n", channel_str, 2*z0, v0);
    printf("Model: %s\n", sigma_model ? "σ-model (project |q|→ρ₀ each step)" : "full model");
    printf("Parameters: ρ₀=%.1f, e=%.1f, λ=%.0f, c=%.1f\n", rho0, e_skyrme, lambda, c_light);
    printf("Time: T_max=%.1f, dt=%.6f (CFL=%.3f)\n", T_max, dt, dt*c_light/h);
    if (settle_time > 0)
        printf("Settling: damp=%.1f for t<%.2f, then conservative\n", damp_gamma, settle_time);
    if (rdamp > 0)
        printf("Radial damping: γ=%.1f (kills breathing mode, preserves topology)\n", rdamp);
    printf("Boundary: %s\n", bc_clamp > 0 ? "Dirichlet (clamped)" : "periodic");
    printf("Threads: %d\n", omp_get_max_threads());

    size_t N3 = (size_t)N * N * N;
    double mem_GB = (double)(N3 * sizeof(Multivector) * 3 + N3 * sizeof(Multivector)) / 1e9;
    printf("Memory: ~%.1f GB (field + velocity + force)\n\n", mem_GB);

    /* Load radial profile */
    if (!profile_file) {
        /* Auto-select: use finite-λ profile if available, else σ-model */
        char auto_file[256];
        snprintf(auto_file, sizeof(auto_file), "profile_finlam_%.0f.dat", lambda);
        FILE *test = fopen(auto_file, "r");
        if (test) {
            fclose(test);
            profile_file = auto_file;
        } else {
            profile_file = "profile_B1.dat";
            printf("WARNING: No finite-λ profile for λ=%.0f, using σ-model profile.\n"
                   "         Generate with: ./finite_lambda_solver -save %.0f -o %s\n",
                   lambda, lambda, auto_file);
        }
    }
    RadialProfile *prof = load_radial_profile(profile_file);
    if (!prof) { fprintf(stderr, "Failed to load %s\n", profile_file); return 1; }
    printf("Loaded profile: %s (%d points, r=[0, %.1f]%s)\n",
           profile_file, prof->n, prof->r_max,
           prof->rho ? ", has ρ(r)" : ", σ-model");

    /* Allocate field and velocity */
    Params params = {rho0, lambda, e_skyrme, 0.0, c_light};
    Field *field = field_alloc(N, L);
    Multivector *vel = (Multivector *)calloc(N3, sizeof(Multivector));
    Multivector *force = (Multivector *)calloc(N3, sizeof(Multivector));
    if (!vel || !force) { fprintf(stderr, "Memory allocation failed\n"); return 1; }

    /* ========== Initialize: Product ansatz ========== */
    printf("\nInitializing two-soliton configuration...\n");

    #pragma omp parallel for collapse(3) schedule(static)
    for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
    for (int k = 0; k < N; k++) {
        double x = -L + (i + 0.5) * h;
        double y = -L + (j + 0.5) * h;
        double z = -L + (k + 0.5) * h;
        int ix = idx(N, i, j, k);

        Quat qp;
        Quat q1, q2;
        if (single_soliton) {
            /* Single hedgehog at origin */
            q1 = hedgehog_q(prof, rho0, x, y, z, 0, 0, 0);
            qp = q1;
        } else {
            /* Soliton 1: at +z₀, moving in -z direction */
            q1 = hedgehog_q(prof, rho0, x, y, z, 0, 0, +z0);

            /* Soliton 2: at -z₀, moving in +z direction */
            q2 = hedgehog_q(prof, rho0, x, y, z, 0, 0, -z0);
            if (antisoliton) {
                q2 = qrev(q2);
            } else if (isorotate) {
                q2.f2 = -q2.f2;
                q2.f3 = -q2.f3;
            }

            /* Product ansatz: q = q1 · q2 / ρ₀ */
            qp = qscale(1.0/rho0, qmul(q1, q2));
        }

        field->psi[ix].s  = qp.s;
        field->psi[ix].f1 = qp.f1;
        field->psi[ix].f2 = qp.f2;
        field->psi[ix].f3 = qp.f3;
        field->psi[ix].j1 = 0;
        field->psi[ix].j2 = 0;
        field->psi[ix].j3 = 0;
        field->psi[ix].p  = 0;

        /* Initial velocity: v = ∂q/∂t|_{t=0}
         * If settling phase is enabled, start at rest (boost after settling).
         * Otherwise, set velocity from product ansatz chain rule. */
        double v_init_speed = (settle_time > 0) ? 0.0 : v0;

        if (v_init_speed > 0) {
            /* For product ansatz q = q1·q2/ρ₀:
             * ∂q/∂t = (∂q1/∂t)·q2/ρ₀ + q1·(∂q2/∂t)/ρ₀
             * Soliton 1 at +z₀ moves -v₀ ẑ, Soliton 2 at -z₀ moves +v₀ ẑ */
            double dz_fd = 0.01;
            Quat q1p = hedgehog_q(prof, rho0, x, y, z+dz_fd, 0, 0, +z0);
            Quat q1m = hedgehog_q(prof, rho0, x, y, z-dz_fd, 0, 0, +z0);
            Quat dq1dz = qscale(0.5/dz_fd, qsub(q1p, q1m));

            Quat q2p_raw = hedgehog_q(prof, rho0, x, y, z+dz_fd, 0, 0, -z0);
            Quat q2m_raw = hedgehog_q(prof, rho0, x, y, z-dz_fd, 0, 0, -z0);
            if (antisoliton) {
                q2p_raw = qrev(q2p_raw);
                q2m_raw = qrev(q2m_raw);
            } else if (isorotate) {
                q2p_raw.f2 = -q2p_raw.f2; q2p_raw.f3 = -q2p_raw.f3;
                q2m_raw.f2 = -q2m_raw.f2; q2m_raw.f3 = -q2m_raw.f3;
            }
            Quat dq2dz = qscale(0.5/dz_fd, qsub(q2p_raw, q2m_raw));

            /* Chain rule: ∂f(x-vt)/∂t = -v·∇f.
             * Soliton 1 has v₁ = -v₀ẑ → ∂q₁/∂t = +v₀ ∂q₁/∂z
             * Soliton 2 has v₂ = +v₀ẑ → ∂q₂/∂t = -v₀ ∂q₂/∂z */
            Quat term1 = qscale(+v_init_speed/rho0, qmul(dq1dz, q2));
            Quat term2 = qscale(-v_init_speed/rho0, qmul(q1, dq2dz));
            Quat v_init = qadd(term1, term2);

            vel[ix].s  = v_init.s;
            vel[ix].f1 = v_init.f1;
            vel[ix].f2 = v_init.f2;
            vel[ix].f3 = v_init.f3;
        } else {
            vel[ix].s = vel[ix].f1 = vel[ix].f2 = vel[ix].f3 = 0;
        }
        vel[ix].j1 = 0;
        vel[ix].j2 = 0;
        vel[ix].j3 = 0;
        vel[ix].p  = 0;
    }

    /* Clamp boundary cells before first energy computation */
    if (bc_clamp > 0) {
        #pragma omp parallel for collapse(3)
        for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
        for (int k = 0; k < N; k++) {
            if (i < bc_clamp || i >= N-bc_clamp ||
                j < bc_clamp || j >= N-bc_clamp ||
                k < bc_clamp || k >= N-bc_clamp) {
                int ix = i*N*N + j*N + k;
                field->psi[ix].s  = rho0;
                field->psi[ix].f1 = 0;
                field->psi[ix].f2 = 0;
                field->psi[ix].f3 = 0;
                vel[ix].s = vel[ix].f1 = vel[ix].f2 = vel[ix].f3 = 0;
            }
        }
    }

    /* Initial diagnostics */
    Energy E0 = field_energy(field, &params);
    double Q0 = field_topological_charge(field, &params);

    /* Kinetic energy: T = (1/2c²) ∫ |v|² d³x */
    double T_kin = 0;
    #pragma omp parallel for reduction(+:T_kin)
    for (int ix = 0; ix < (int)N3; ix++) {
        T_kin += 0.5 * mv_dot(vel[ix], vel[ix]);
    }
    T_kin *= h*h*h / (c_light*c_light);

    double z1 = find_soliton_z(field, rho0, 1);
    double z2 = find_soliton_z(field, rho0, -1);

    printf("\n--- Initial state ---\n");
    printf("E_pot    = %.6f  (E2=%.4f, E4=%.4f, EV=%.4f)\n",
           E0.Etotal, E0.E2, E0.E4, E0.EV);
    printf("E_kin    = %.6f\n", T_kin);
    printf("E_total  = %.6f\n", E0.Etotal + T_kin);
    printf("Q        = %.6f  (expected: %d)\n", Q0, single_soliton ? 1 : (antisoliton ? 0 : 2));
    printf("z1       = %+.3f  (soliton 1)\n", z1);
    printf("z2       = %+.3f  (soliton 2)\n", z2);
    printf("sep      = %.3f\n", z1 - z2);

    /* Open output file */
    char outname[256];
    snprintf(outname, sizeof(outname), "scatter_%s_v%.2f.dat",
             antisoliton ? "anti" : "same", v0);
    FILE *fout = fopen(outname, "w");
    fprintf(fout, "# Soliton scattering: %s, v=%.3f, N=%d, L=%.1f, lambda=%.0f\n",
            antisoliton ? "B+B̄" : "B+B", v0, N, L, lambda);
    fprintf(fout, "# Columns: t  E_pot  E_kin  E_tot  Q  z1  z2  sep  E2  E4  EV\n");
    fprintf(fout, "%.6f  %.6f  %.6f  %.6f  %.6f  %+.4f  %+.4f  %.4f  %.4f  %.4f  %.4f\n",
            0.0, E0.Etotal, T_kin, E0.Etotal+T_kin, Q0, z1, z2, z1-z2,
            E0.E2, E0.E4, E0.EV);

    /* ========== Time evolution: Leapfrog ========== */
    printf("\n--- Time evolution ---\n");
    printf("step     t       E_pot     E_kin     E_tot     Q       sep     z1      z2\n");

    /* First half-step for velocity (to stagger v at t+dt/2) */
    field_gradient(field, &params, force);

    /* Project force tangentially for σ-model */
    if (sigma_model) {
        #pragma omp parallel for
        for (int ix = 0; ix < (int)N3; ix++) {
            double s = field->psi[ix].s, f1 = field->psi[ix].f1;
            double f2 = field->psi[ix].f2, f3 = field->psi[ix].f3;
            double n2 = s*s + f1*f1 + f2*f2 + f3*f3;
            if (n2 > 1e-20) {
                double fdq = (force[ix].s*s + force[ix].f1*f1
                            + force[ix].f2*f2 + force[ix].f3*f3) / n2;
                force[ix].s  -= fdq * s;
                force[ix].f1 -= fdq * f1;
                force[ix].f2 -= fdq * f2;
                force[ix].f3 -= fdq * f3;
            }
        }
    }

    /* Project initial velocity tangential to constraint */
    if (sigma_model) {
        #pragma omp parallel for
        for (int ix = 0; ix < (int)N3; ix++) {
            double s = field->psi[ix].s, f1 = field->psi[ix].f1;
            double f2 = field->psi[ix].f2, f3 = field->psi[ix].f3;
            double n2 = s*s + f1*f1 + f2*f2 + f3*f3;
            if (n2 > 1e-20) {
                double vdq = (vel[ix].s*s + vel[ix].f1*f1
                            + vel[ix].f2*f2 + vel[ix].f3*f3) / n2;
                vel[ix].s  -= vdq * s;
                vel[ix].f1 -= vdq * f1;
                vel[ix].f2 -= vdq * f2;
                vel[ix].f3 -= vdq * f3;
            }
        }
    }

    #pragma omp parallel for
    for (int ix = 0; ix < (int)N3; ix++) {
        /* v(dt/2) = v(0) + (dt/2) × c² × F(q(0)) */
        double fac = 0.5 * dt * c_light * c_light;
        vel[ix].s  += fac * force[ix].s;
        vel[ix].f1 += fac * force[ix].f1;
        vel[ix].f2 += fac * force[ix].f2;
        vel[ix].f3 += fac * force[ix].f3;
        /* Leave degenerate sector at zero */
    }

    int n_steps = (int)(T_max / dt);
    double t = 0;
    double E_total_0 = E0.Etotal + T_kin;
    int settled = (settle_time <= 0);  /* 1 once settling is done */
    int boost_applied = (v0 == 0 || settle_time <= 0);  /* 1 once boost given */

    for (int step = 1; step <= n_steps; step++) {
        t = step * dt;

        /* Check if settling just ended — apply velocity boost */
        if (!settled && t >= settle_time) {
            settled = 1;
            printf("--- Settling complete at t=%.3f ---\n", t);

            /* Recalibrate energy reference after settling */
            Energy E_settled = field_energy(field, &params);
            double T_settled = 0;
            #pragma omp parallel for reduction(+:T_settled)
            for (int ix = 0; ix < (int)N3; ix++)
                T_settled += 0.5 * mv_dot(vel[ix], vel[ix]);
            T_settled *= h*h*h / (c_light*c_light);
            E_total_0 = E_settled.Etotal + T_settled;
            double Q_settled = field_topological_charge(field, &params);
            printf("E_pot=%.4f, E_kin=%.4f, Q=%.4f\n",
                   E_settled.Etotal, T_settled, Q_settled);

            /* Apply velocity boost if v0 > 0 */
            if (v0 > 0 && !boost_applied) {
                boost_applied = 1;
                printf("Boosting solitons to v=%.3fc...\n", v0);
                /* Add ±v₀ẑ boost: v += ±v₀ ∂q/∂z */
                double dz_fd = h;
                #pragma omp parallel for collapse(3) schedule(static)
                for (int i = 0; i < N; i++)
                for (int j = 0; j < N; j++)
                for (int k = 0; k < N; k++) {
                    double z = -L + (k + 0.5) * h;
                    int ix = idx(N, i, j, k);
                    /* Finite diff ∂q/∂z from field */
                    int kp = (k+1 < N) ? k+1 : k;
                    int km = (k > 0) ? k-1 : k;
                    double dzeff = (kp - km) * h;
                    if (dzeff < 1e-15) continue;
                    Multivector dqdz;
                    dqdz.s  = (field->psi[idx(N,i,j,kp)].s  - field->psi[idx(N,i,j,km)].s)  / dzeff;
                    dqdz.f1 = (field->psi[idx(N,i,j,kp)].f1 - field->psi[idx(N,i,j,km)].f1) / dzeff;
                    dqdz.f2 = (field->psi[idx(N,i,j,kp)].f2 - field->psi[idx(N,i,j,km)].f2) / dzeff;
                    dqdz.f3 = (field->psi[idx(N,i,j,kp)].f3 - field->psi[idx(N,i,j,km)].f3) / dzeff;
                    /* Chain rule: ∂f(x-vt)/∂t = -v·∇f
                     * Soliton 1 at z>0: v₁=-v₀ẑ → ∂q₁/∂t = +v₀ ∂q/∂z
                     * Soliton 2 at z<0: v₂=+v₀ẑ → ∂q₂/∂t = -v₀ ∂q/∂z */
                    double sign = (z > 0) ? +v0 : -v0;
                    vel[ix].s  += sign * dqdz.s;
                    vel[ix].f1 += sign * dqdz.f1;
                    vel[ix].f2 += sign * dqdz.f2;
                    vel[ix].f3 += sign * dqdz.f3;
                }

                /* Update energy reference with boost */
                T_settled = 0;
                #pragma omp parallel for reduction(+:T_settled)
                for (int ix = 0; ix < (int)N3; ix++)
                    T_settled += 0.5 * mv_dot(vel[ix], vel[ix]);
                T_settled *= h*h*h / (c_light*c_light);
                E_total_0 = E_settled.Etotal + T_settled;
                printf("After boost: E_kin=%.4f, E_total=%.4f\n", T_settled, E_total_0);
            }
        }

        /* q(t+dt) = q(t) + dt × v(t+dt/2) */
        #pragma omp parallel for
        for (int ix = 0; ix < (int)N3; ix++) {
            field->psi[ix].s  += dt * vel[ix].s;
            field->psi[ix].f1 += dt * vel[ix].f1;
            field->psi[ix].f2 += dt * vel[ix].f2;
            field->psi[ix].f3 += dt * vel[ix].f3;
        }

        /* σ-model: project |q| → ρ₀ and velocity perpendicular to q */
        if (sigma_model) {
            #pragma omp parallel for
            for (int ix = 0; ix < (int)N3; ix++) {
                double s = field->psi[ix].s, f1 = field->psi[ix].f1;
                double f2 = field->psi[ix].f2, f3 = field->psi[ix].f3;
                double norm = sqrt(s*s + f1*f1 + f2*f2 + f3*f3);
                if (norm > 1e-15) {
                    double scale = rho0 / norm;
                    field->psi[ix].s  *= scale;
                    field->psi[ix].f1 *= scale;
                    field->psi[ix].f2 *= scale;
                    field->psi[ix].f3 *= scale;

                    /* Project velocity: v⊥ = v - (v·q̂)q̂ */
                    double qh_s = field->psi[ix].s/rho0, qh_f1 = field->psi[ix].f1/rho0;
                    double qh_f2 = field->psi[ix].f2/rho0, qh_f3 = field->psi[ix].f3/rho0;
                    double vdq = vel[ix].s*qh_s + vel[ix].f1*qh_f1
                               + vel[ix].f2*qh_f2 + vel[ix].f3*qh_f3;
                    vel[ix].s  -= vdq * qh_s;
                    vel[ix].f1 -= vdq * qh_f1;
                    vel[ix].f2 -= vdq * qh_f2;
                    vel[ix].f3 -= vdq * qh_f3;
                }
            }
        }

        /* Boundary clamping: set outer cells to vacuum q=(ρ₀,0,0,0), v=0
         * This gives effective Dirichlet BCs, eliminating the gradient
         * discontinuity from periodic wrapping of the hedgehog tail. */
        if (bc_clamp > 0) {
            #pragma omp parallel for collapse(3)
            for (int i = 0; i < N; i++)
            for (int j = 0; j < N; j++)
            for (int k = 0; k < N; k++) {
                if (i < bc_clamp || i >= N-bc_clamp ||
                    j < bc_clamp || j >= N-bc_clamp ||
                    k < bc_clamp || k >= N-bc_clamp) {
                    int ix = i*N*N + j*N + k;
                    field->psi[ix].s  = rho0;
                    field->psi[ix].f1 = 0;
                    field->psi[ix].f2 = 0;
                    field->psi[ix].f3 = 0;
                    vel[ix].s = vel[ix].f1 = vel[ix].f2 = vel[ix].f3 = 0;
                }
            }
        }

        /* Compute force F(q(t+dt)) */
        field_gradient(field, &params, force);

        /* σ-model: project force tangential to constraint |q|=ρ₀ */
        if (sigma_model) {
            #pragma omp parallel for
            for (int ix = 0; ix < (int)N3; ix++) {
                double s = field->psi[ix].s, f1 = field->psi[ix].f1;
                double f2 = field->psi[ix].f2, f3 = field->psi[ix].f3;
                double n2 = s*s + f1*f1 + f2*f2 + f3*f3;
                if (n2 > 1e-20) {
                    /* F_T = F - (F·q/|q|²)q */
                    double fdq = (force[ix].s*s + force[ix].f1*f1
                                + force[ix].f2*f2 + force[ix].f3*f3) / n2;
                    force[ix].s  -= fdq * s;
                    force[ix].f1 -= fdq * f1;
                    force[ix].f2 -= fdq * f2;
                    force[ix].f3 -= fdq * f3;
                }
            }
        }

        /* v(t+3dt/2) = v(t+dt/2) + dt × c² × F(q(t+dt)) */
        #pragma omp parallel for
        for (int ix = 0; ix < (int)N3; ix++) {
            double fac = dt * c_light * c_light;
            vel[ix].s  += fac * force[ix].s;
            vel[ix].f1 += fac * force[ix].f1;
            vel[ix].f2 += fac * force[ix].f2;
            vel[ix].f3 += fac * force[ix].f3;
        }

        /* Selective radial damping: damp v·q̂ (breathing mode) only.
         * Decomposes v = v_tangential + v_radial, damps v_radial by exp(-rdamp·dt).
         * This kills the breathing mode while preserving angular/translational dynamics. */
        if (rdamp > 0) {
            double rd = exp(-rdamp * dt);  /* radial damping factor per step */
            #pragma omp parallel for
            for (int ix = 0; ix < (int)N3; ix++) {
                double qs = field->psi[ix].s, qf1 = field->psi[ix].f1;
                double qf2 = field->psi[ix].f2, qf3 = field->psi[ix].f3;
                double n2 = qs*qs + qf1*qf1 + qf2*qf2 + qf3*qf3;
                if (n2 < 1e-20) continue;
                /* v_radial = (v·q/|q|²) q */
                double vdq = (vel[ix].s*qs + vel[ix].f1*qf1
                            + vel[ix].f2*qf2 + vel[ix].f3*qf3) / n2;
                /* v -= (1-rd) * v_radial */
                double fac = (1.0 - rd) * vdq;
                vel[ix].s  -= fac * qs;
                vel[ix].f1 -= fac * qf1;
                vel[ix].f2 -= fac * qf2;
                vel[ix].f3 -= fac * qf3;
            }
        }

        /* Velocity damping during settling phase */
        if (damp_gamma > 0 && t < settle_time) {
            double damp = 1.0 - damp_gamma * dt;
            if (damp < 0) damp = 0;
            #pragma omp parallel for
            for (int ix = 0; ix < (int)N3; ix++) {
                vel[ix].s  *= damp;
                vel[ix].f1 *= damp;
                vel[ix].f2 *= damp;
                vel[ix].f3 *= damp;
            }
        }

        /* Diagnostics */
        if (step % output_interval == 0 || step == n_steps) {
            Energy E = field_energy(field, &params);
            double Q = field_topological_charge(field, &params);

            /* Kinetic energy */
            double T = 0;
            #pragma omp parallel for reduction(+:T)
            for (int ix = 0; ix < (int)N3; ix++) {
                T += 0.5 * mv_dot(vel[ix], vel[ix]);
            }
            T *= h*h*h / (c_light*c_light);

            z1 = find_soliton_z(field, rho0, 1);
            z2 = find_soliton_z(field, rho0, -1);
            double sep = z1 - z2;

            /* 3D soliton tracking */
            double x1_3d, y1_3d, z1_3d, x2_3d, y2_3d, z2_3d;
            find_soliton_3d(field, rho0, &x1_3d, &y1_3d, &z1_3d,
                            &x2_3d, &y2_3d, &z2_3d);
            double dx = x1_3d - x2_3d, dy = y1_3d - y2_3d, dz = z1_3d - z2_3d;
            double sep3d = sqrt(dx*dx + dy*dy + dz*dz);

            double E_tot = E.Etotal + T;
            double dE = (E_tot - E_total_0) / E_total_0;

            printf("%5d  %7.3f  %9.4f  %9.4f  %9.4f  %7.4f  %7.3f  %+6.3f  %+6.3f",
                   step, t, E.Etotal, T, E_tot, Q, sep, z1, z2);
            printf("  3d:(%.2f,%.2f,%.2f)|(%.2f,%.2f,%.2f) r=%.3f",
                   x1_3d, y1_3d, z1_3d, x2_3d, y2_3d, z2_3d, sep3d);
            if (fabs(dE) > 0.01) printf("  ΔE=%.2e!", dE);
            printf("\n");

            fprintf(fout, "%.6f  %.6f  %.6f  %.6f  %.6f  %+.4f  %+.4f  %.4f  %.4f  %.4f  %.4f\n",
                    t, E.Etotal, T, E_tot, Q, z1, z2, sep,
                    E.E2, E.E4, E.EV);
            fflush(fout);

            /* Early termination if energy conservation is badly violated */
            if (fabs(dE) > 0.5) {
                printf("*** Energy conservation violated by >50%% — aborting ***\n");
                break;
            }

            /* Warn if topology degrades significantly from initial value */
            if (!antisoliton && Q < Q0 * 0.95) {
                printf("*** Topology loss detected (Q=%.4f, init=%.4f) ***\n", Q, Q0);
            }
        }
    }

    fclose(fout);
    printf("\nData written to %s\n", outname);

    /* Final state */
    Energy Ef = field_energy(field, &params);
    double Qf = field_topological_charge(field, &params);
    double Tf = 0;
    #pragma omp parallel for reduction(+:Tf)
    for (int ix = 0; ix < (int)N3; ix++) {
        Tf += 0.5 * mv_dot(vel[ix], vel[ix]);
    }
    Tf *= h*h*h / (c_light*c_light);

    printf("\n--- Final state (t=%.2f) ---\n", t);
    printf("E_pot    = %.6f\n", Ef.Etotal);
    printf("E_kin    = %.6f\n", Tf);
    printf("E_total  = %.6f  (initial: %.6f, ΔE/E=%.2e)\n",
           Ef.Etotal+Tf, E_total_0, (Ef.Etotal+Tf-E_total_0)/E_total_0);
    printf("Q        = %.6f  (initial: %.6f)\n", Qf, Q0);

    /* Cleanup */
    free(vel);
    free(force);
    field_free(field);
    free_radial_profile(prof);

    return 0;
}
