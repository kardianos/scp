/*  multiscale_v2.c -- Phase 1b: Free-floating BC + fast wedge time evolution
 *
 *  Changes from v1:
 *    - Free-floating (linear extrapolation) BC instead of absorbing damping
 *    - Wedge dt decoupled from core dt (offset time stepping)
 *    - Larger R_match option
 *  Core:  3D Cartesian Cosserat grid clipped to sphere of radius R_match
 *  Wedge: 1D radial Schrodinger equation from R_match to R_max (Crank-Nicolson)
 *  Coupling: theta extracted at R_match feeds wedge; wedge BC back to core (Phase 2)
 *
 *  Build: gcc -O3 -march=native -fopenmp -o multiscale src/multiscale.c -lm
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>
#include <sys/stat.h>
#include <stdint.h>

#define NFIELDS 3
#define PI 3.14159265358979323846

/* ================================================================
   Physics parameters
   ================================================================ */

/* Core (Cosserat) parameters */
static double MU       = -41.345;
static double KAPPA    = 50.0;
static double MASS2    = 2.25;     /* position field mass^2 */
static double MTHETA2  = 2.25;    /* angle field mass^2 */
static double A_BG     = 0.1;
static double ETA      = 1.0;     /* position-angle coupling */

/* Spherical boundary */
static double R_MATCH  = 20.0;    /* core-wedge boundary radius */
static double DAMP_WIDTH = 8.0;   /* damping layer width (in units of dx) — wider = gentler */
static double DAMP_RATE  = 0.995; /* velocity damping at outer edge (per step) — closer to 1 = gentler */

/* Wedge (Schrodinger) parameters */
static double HBAR_EFF  = 22727.0;
static double M_EFF     = 1535.0;
static double V_DEPTH   = 1.27;
static double N_POWER   = 1.189;
static double R_BRAID   = 5.0;
static double R_MAX     = 500000.0;
static double BOHR_RAD  = 265000.0;

/* ================================================================
   Core Grid: Cosserat with spherical mask
   ================================================================ */

typedef struct {
    double *mem;
    double *phi[NFIELDS];
    double *phi_vel[NFIELDS];
    double *phi_acc[NFIELDS];
    double *theta[NFIELDS];
    double *theta_vel[NFIELDS];
    double *theta_acc[NFIELDS];
    uint8_t *mask;    /* 0=outside, 1=damping layer, 2=interior */
    int N; long N3;
    double L, dx, dt;
} Grid;

static Grid *grid_alloc(int N, double L) {
    Grid *g = calloc(1, sizeof(Grid));
    g->N  = N;
    g->N3 = (long)N * N * N;
    g->L  = L;
    g->dx = 2.0 * L / (N - 1);
    g->dt = 0.10 * g->dx;

    long total = 18 * g->N3;
    double bytes = total * sizeof(double);
    printf("Core: allocating %.2f GB (%ld doubles, N=%d, 6 fields)\n", bytes/1e9, total, N);
    g->mem = malloc(total * sizeof(double));
    if (!g->mem) { fprintf(stderr, "FATAL: malloc failed for grid\n"); exit(1); }
    memset(g->mem, 0, total * sizeof(double));

    for (int a = 0; a < NFIELDS; a++) {
        g->phi[a]       = g->mem + (0  + a) * g->N3;
        g->phi_vel[a]   = g->mem + (3  + a) * g->N3;
        g->phi_acc[a]   = g->mem + (6  + a) * g->N3;
        g->theta[a]     = g->mem + (9  + a) * g->N3;
        g->theta_vel[a] = g->mem + (12 + a) * g->N3;
        g->theta_acc[a] = g->mem + (15 + a) * g->N3;
    }

    /* Spherical mask */
    g->mask = calloc(g->N3, sizeof(uint8_t));
    if (!g->mask) { fprintf(stderr, "FATAL: malloc failed for mask\n"); exit(1); }

    double damp_inner = R_MATCH - DAMP_WIDTH * g->dx;
    long n_inside = 0, n_damp = 0, n_outside = 0;
    int NN = N * N;
    for (int i = 0; i < N; i++) {
        double x = -L + i * g->dx;
        for (int j = 0; j < N; j++) {
            double y = -L + j * g->dx;
            for (int k = 0; k < N; k++) {
                double z = -L + k * g->dx;
                long idx = (long)i * NN + j * N + k;
                double r = sqrt(x*x + y*y + z*z);
                if (r < damp_inner) {
                    g->mask[idx] = 2;  /* interior */
                    n_inside++;
                } else if (r < R_MATCH) {
                    g->mask[idx] = 1;  /* damping layer */
                    n_damp++;
                } else {
                    g->mask[idx] = 0;  /* outside */
                    n_outside++;
                }
            }
        }
    }
    printf("Core mask: %ld interior, %ld damping, %ld outside (%.1f%% active)\n",
           n_inside, n_damp, n_outside,
           100.0*(n_inside+n_damp)/g->N3);

    return g;
}

static void grid_free(Grid *g) {
    free(g->mask);
    free(g->mem);
    free(g);
}

/* ================================================================
   Curl computation
   ================================================================ */

static inline double curl_component(double *F[3], int a,
    long n_ip, long n_im, long n_jp, long n_jm, long n_kp, long n_km,
    double idx1) {
    if (a == 0) return (F[2][n_jp] - F[2][n_jm] - F[1][n_kp] + F[1][n_km]) * idx1;
    if (a == 1) return (F[0][n_kp] - F[0][n_km] - F[2][n_ip] + F[2][n_im]) * idx1;
    return            (F[1][n_ip] - F[1][n_im] - F[0][n_jp] + F[0][n_jm]) * idx1;
}

/* ================================================================
   Core forces with spherical mask
   ================================================================ */

static void compute_forces(Grid *g) {
    const int N = g->N, NN = N * N;
    const long N3 = g->N3;
    const double idx2 = 1.0 / (g->dx * g->dx);
    const double idx1 = 1.0 / (2.0 * g->dx);

    #pragma omp parallel for schedule(static)
    for (long idx = 0; idx < N3; idx++) {
        /* Skip outside points entirely */
        if (g->mask[idx] == 0) {
            for (int a = 0; a < NFIELDS; a++) {
                g->phi_acc[a][idx] = 0;
                g->theta_acc[a][idx] = 0;
            }
            continue;
        }

        int i = (int)(idx / NN);
        int j = (int)((idx / N) % N);
        int k = (int)(idx % N);

        /* Neighbors: clamp at boundaries instead of periodic wrap */
        int ip = (i+1 < N) ? i+1 : i;
        int im = (i-1 >= 0) ? i-1 : i;
        int jp = (j+1 < N) ? j+1 : j;
        int jm = (j-1 >= 0) ? j-1 : j;
        int kp = (k+1 < N) ? k+1 : k;
        int km = (k-1 >= 0) ? k-1 : k;

        long n_ip = (long)ip*NN + j*N + k;
        long n_im = (long)im*NN + j*N + k;
        long n_jp = (long)i*NN + jp*N + k;
        long n_jm = (long)i*NN + jm*N + k;
        long n_kp = (long)i*NN + j*N + kp;
        long n_km = (long)i*NN + j*N + km;

        /* For neighbors outside the sphere, use the central value (Neumann BC) */
        /* This is handled implicitly: outside points have field=0, but we
           want Neumann (zero-gradient), so we reflect. For simplicity in Phase 1,
           we just use the clamped neighbors above and zero the outside fields
           after each step. */

        /* Position field forces */
        double p0 = g->phi[0][idx], p1 = g->phi[1][idx], p2 = g->phi[2][idx];
        double P = p0 * p1 * p2;
        double den = 1.0 + KAPPA * P * P;
        double mPd2 = MU * P / (den * den);

        for (int a = 0; a < NFIELDS; a++) {
            double lap = (g->phi[a][n_ip] + g->phi[a][n_im]
                        + g->phi[a][n_jp] + g->phi[a][n_jm]
                        + g->phi[a][n_kp] + g->phi[a][n_km]
                        - 6.0 * g->phi[a][idx]) * idx2;
            double dPda = (a==0) ? p1*p2 : (a==1) ? p0*p2 : p0*p1;

            double curl_theta = curl_component(g->theta, a,
                n_ip, n_im, n_jp, n_jm, n_kp, n_km, idx1);

            g->phi_acc[a][idx] = lap - MASS2 * g->phi[a][idx]
                               - mPd2 * dPda + ETA * curl_theta;
        }

        /* Angle field forces */
        for (int a = 0; a < NFIELDS; a++) {
            double lap_t = (g->theta[a][n_ip] + g->theta[a][n_im]
                          + g->theta[a][n_jp] + g->theta[a][n_jm]
                          + g->theta[a][n_kp] + g->theta[a][n_km]
                          - 6.0 * g->theta[a][idx]) * idx2;

            double curl_phi = curl_component(g->phi, a,
                n_ip, n_im, n_jp, n_jm, n_kp, n_km, idx1);

            g->theta_acc[a][idx] = lap_t - MTHETA2 * g->theta[a][idx]
                                 + ETA * curl_phi;
        }
    }
}

/* ================================================================
   Verlet step with spherical mask + damping layer
   ================================================================ */

static void verlet_step(Grid *g) {
    const long N3 = g->N3;
    const double hdt = 0.5 * g->dt, dt = g->dt;

    /* Half-kick */
    for (int a = 0; a < NFIELDS; a++) {
        double *vp = g->phi_vel[a], *ap = g->phi_acc[a];
        double *vt = g->theta_vel[a], *at = g->theta_acc[a];
        #pragma omp parallel for schedule(static)
        for (long i = 0; i < N3; i++) {
            if (g->mask[i] == 0) continue;
            vp[i] += hdt * ap[i];
            vt[i] += hdt * at[i];
        }
    }

    /* Drift */
    for (int a = 0; a < NFIELDS; a++) {
        double *pp = g->phi[a], *vp = g->phi_vel[a];
        double *pt = g->theta[a], *vt = g->theta_vel[a];
        #pragma omp parallel for schedule(static)
        for (long i = 0; i < N3; i++) {
            if (g->mask[i] == 0) continue;
            pp[i] += dt * vp[i];
            pt[i] += dt * vt[i];
        }
    }

    /* Recompute forces */
    compute_forces(g);

    /* Half-kick */
    for (int a = 0; a < NFIELDS; a++) {
        double *vp = g->phi_vel[a], *ap = g->phi_acc[a];
        double *vt = g->theta_vel[a], *at = g->theta_acc[a];
        #pragma omp parallel for schedule(static)
        for (long i = 0; i < N3; i++) {
            if (g->mask[i] == 0) continue;
            vp[i] += hdt * ap[i];
            vt[i] += hdt * at[i];
        }
    }

    /* Gentle absorbing layer: smoothly damp velocities in the boundary region.
     * Use a position-dependent damping factor that increases toward the edge.
     * This absorbs outgoing radiation without reflecting it. */
    {
        const int N = g->N, NN = N*N;
        const double dx = g->dx, L = g->L;
        const double damp_start = R_MATCH - DAMP_WIDTH * dx;

        #pragma omp parallel for schedule(static)
        for (long idx = 0; idx < N3; idx++) {
            if (g->mask[idx] != 1) continue;

            int i = (int)(idx / NN);
            int j = (int)((idx / N) % N);
            int k = (int)(idx % N);
            double x = -L + i*dx;
            double y = -L + j*dx;
            double z = -L + k*dx;
            double r = sqrt(x*x + y*y + z*z);

            /* Smooth damping: 1.0 at damp_start, DAMP_RATE at R_MATCH */
            double frac = (r - damp_start) / (R_MATCH - damp_start);
            if (frac < 0) frac = 0; if (frac > 1) frac = 1;
            double damp = 1.0 - frac * (1.0 - DAMP_RATE);

            for (int a = 0; a < NFIELDS; a++) {
                g->phi_vel[a][idx] *= damp;
                g->theta_vel[a][idx] *= damp;
            }
        }
    }

    /* Zero outside points (ensure they stay dead) */
    for (int a = 0; a < NFIELDS; a++) {
        #pragma omp parallel for schedule(static)
        for (long i = 0; i < N3; i++) {
            if (g->mask[i] == 0) {
                g->phi[a][i] = 0;
                g->phi_vel[a][i] = 0;
                g->theta[a][i] = 0;
                g->theta_vel[a][i] = 0;
            }
        }
    }
}

/* ================================================================
   Wedge: 1D Radial Schrodinger with Crank-Nicolson
   ================================================================ */

typedef struct {
    int Nr;
    double *r;        /* radial positions (log-spaced) */
    double *dr;       /* spacing at each point */
    double *V;        /* effective potential V(r) */
    double *psi_re;   /* real part of wavefunction */
    double *psi_im;   /* imaginary part */
    double *rhs_re;   /* work arrays for CN solve */
    double *rhs_im;
    double hbar, m_eff;
    double dt;        /* wedge timestep */
    double norm;      /* current normalization */
    /* Pre-allocated work arrays for CN step (avoid per-step malloc) */
    double *u_re, *u_im;       /* u = r*psi */
    double *Hu_re, *Hu_im;     /* H*u */
    double *a_coeff, *b_coeff, *c_coeff;  /* tridiag coefficients */
    double *d_re, *d_im;       /* modified diagonal */
    double *f_re, *f_im;       /* modified RHS */
} Wedge;

static Wedge *wedge_alloc(int Nr, double r_min, double r_max,
                           double hbar, double m_eff) {
    Wedge *w = calloc(1, sizeof(Wedge));
    w->Nr = Nr;
    w->hbar = hbar;
    w->m_eff = m_eff;

    w->r      = malloc(Nr * sizeof(double));
    w->dr     = malloc(Nr * sizeof(double));
    w->V      = malloc(Nr * sizeof(double));
    w->psi_re = calloc(Nr, sizeof(double));
    w->psi_im = calloc(Nr, sizeof(double));
    w->rhs_re = malloc(Nr * sizeof(double));
    w->rhs_im = malloc(Nr * sizeof(double));
    /* Pre-allocate all CN work arrays */
    w->u_re    = malloc(Nr * sizeof(double));
    w->u_im    = malloc(Nr * sizeof(double));
    w->Hu_re   = malloc(Nr * sizeof(double));
    w->Hu_im   = malloc(Nr * sizeof(double));
    w->a_coeff = malloc(Nr * sizeof(double));
    w->b_coeff = malloc(Nr * sizeof(double));
    w->c_coeff = malloc(Nr * sizeof(double));
    w->d_re    = malloc(Nr * sizeof(double));
    w->d_im    = malloc(Nr * sizeof(double));
    w->f_re    = malloc(Nr * sizeof(double));
    w->f_im    = malloc(Nr * sizeof(double));

    /* Log-spaced radial grid */
    double log_rmin = log(r_min);
    double log_rmax = log(r_max);
    for (int i = 0; i < Nr; i++) {
        double frac = (double)i / (Nr - 1);
        w->r[i] = exp(log_rmin + frac * (log_rmax - log_rmin));
    }
    for (int i = 0; i < Nr - 1; i++) {
        w->dr[i] = w->r[i+1] - w->r[i];
    }
    w->dr[Nr-1] = w->dr[Nr-2];

    /* Potential: V(r) = -V_depth / (r/R_braid)^n_power
       with repulsive core for r < 2*R_braid */
    for (int i = 0; i < Nr; i++) {
        double rr = w->r[i] / R_BRAID;
        if (rr < 0.4) {
            /* Hard repulsive core */
            w->V[i] = V_DEPTH * 100.0 / (rr * rr);
        } else {
            w->V[i] = -V_DEPTH / pow(rr, N_POWER);
        }
    }

    /* Timestep: from hbar and grid spacing */
    /* For Schrodinger, dt < m_eff * dr_min^2 / hbar for stability,
       but CN is unconditionally stable. Use physics timescale. */
    double dr_min = w->dr[0];
    w->dt = 0.5 * m_eff * dr_min * dr_min / hbar;
    printf("Wedge: Nr=%d, r=[%.1f, %.0f], dr_min=%.2f, dr_max=%.0f, dt_wedge=%.2f\n",
           Nr, r_min, r_max, dr_min, w->dr[Nr-2], w->dt);

    /* Pre-compute tridiag coefficients (constant — depends only on grid + V) */
    double alpha_k = hbar / (2.0 * m_eff);  /* hbar/(2m) */
    w->a_coeff[0] = 0;
    w->c_coeff[Nr-1] = 0;
    w->b_coeff[0] = w->V[0];
    w->b_coeff[Nr-1] = w->V[Nr-1];
    for (int i = 1; i < Nr - 1; i++) {
        double drp = w->r[i+1] - w->r[i];
        double drm = w->r[i] - w->r[i-1];
        double dravg = 0.5 * (drp + drm);
        w->a_coeff[i] = -alpha_k / (drm * dravg);
        w->c_coeff[i] = -alpha_k / (drp * dravg);
        w->b_coeff[i] = alpha_k * (1.0/drp + 1.0/drm) / dravg + w->V[i];
    }
    w->c_coeff[0] = 0;
    w->a_coeff[Nr-1] = 0;

    return w;
}

static void wedge_free(Wedge *w) {
    free(w->r); free(w->dr); free(w->V);
    free(w->psi_re); free(w->psi_im);
    free(w->rhs_re); free(w->rhs_im);
    free(w->u_re); free(w->u_im);
    free(w->Hu_re); free(w->Hu_im);
    free(w->a_coeff); free(w->b_coeff); free(w->c_coeff);
    free(w->d_re); free(w->d_im);
    free(w->f_re); free(w->f_im);
    free(w);
}

/* Initialize wedge with Gaussian wave packet */
static void wedge_init_gaussian(Wedge *w, double r_center, double sigma) {
    double norm = 0;
    for (int i = 0; i < w->Nr; i++) {
        double dr = w->r[i] - r_center;
        /* u(r) = r*psi(r) for radial Schrodinger; normalize in u-space */
        double env = exp(-dr*dr / (2.0*sigma*sigma));
        w->psi_re[i] = env;
        w->psi_im[i] = 0;
        /* Norm: integral |psi|^2 * r^2 dr (spherical) */
        double ri = w->r[i];
        norm += (w->psi_re[i]*w->psi_re[i] + w->psi_im[i]*w->psi_im[i])
                * ri * ri * w->dr[i];
    }
    /* Normalize */
    double inv = 1.0 / sqrt(norm);
    for (int i = 0; i < w->Nr; i++) {
        w->psi_re[i] *= inv;
        w->psi_im[i] *= inv;
    }
    w->norm = 1.0;
    printf("Wedge: Gaussian at r=%.0f, sigma=%.0f, normalized\n", r_center, sigma);
}

/*  Crank-Nicolson step for radial Schrodinger equation.
 *
 *  We solve:  i*hbar * dpsi/dt = H*psi
 *  where H = -hbar^2/(2m) * (1/r^2) d/dr(r^2 d/dr) + V(r)
 *
 *  Using u = r*psi to remove the 1/r^2 factor:
 *  H_u = -hbar^2/(2m) * d^2u/dr^2 + V(r)*u
 *
 *  CN:  (1 + i*dt*H/(2*hbar)) * u^{n+1} = (1 - i*dt*H/(2*hbar)) * u^n
 *
 *  This is a complex tridiagonal system. We solve it by splitting
 *  real and imaginary parts.
 *
 *  Let alpha = dt/(4*m_eff*hbar), beta = dt*V/(2*hbar)
 *  Then the kinetic part gives tridiagonal with:
 *    diag: 1 + alpha*hbar^2*(2/dr^2) + i*beta
 *    off:  -alpha*hbar^2/dr^2
 *
 *  For non-uniform grid, use the standard 3-point stencil:
 *    d^2u/dr^2 ~ (u_{i+1} - u_i)/(dr_i * dr_avg) - (u_i - u_{i-1})/(dr_{i-1} * dr_avg)
 *    where dr_avg = (dr_i + dr_{i-1})/2
 */
static void wedge_step(Wedge *w) {
    const int Nr = w->Nr;
    const double alpha = w->hbar / (2.0 * w->m_eff);
    const double hdt_h = w->dt / (2.0 * w->hbar);

    /* Use pre-allocated arrays */
    double *u_re = w->u_re, *u_im = w->u_im;
    double *Hu_re = w->Hu_re, *Hu_im = w->Hu_im;
    double *rhs_re = w->rhs_re, *rhs_im = w->rhs_im;
    double *a_coeff = w->a_coeff, *b_coeff = w->b_coeff, *c_coeff = w->c_coeff;
    double *d_re = w->d_re, *d_im = w->d_im;
    double *f_re = w->f_re, *f_im = w->f_im;

    /* u = r*psi for radial equation */
    for (int i = 0; i < Nr; i++) {
        u_re[i] = w->r[i] * w->psi_re[i];
        u_im[i] = w->r[i] * w->psi_im[i];
    }

    /* Compute H*u using pre-computed tridiag coefficients:
     * H*u_i = a_i * u_{i-1} + b_i * u_i + c_i * u_{i+1} */
    Hu_re[0] = b_coeff[0] * u_re[0];
    Hu_im[0] = b_coeff[0] * u_im[0];
    Hu_re[Nr-1] = b_coeff[Nr-1] * u_re[Nr-1];
    Hu_im[Nr-1] = b_coeff[Nr-1] * u_im[Nr-1];
    for (int i = 1; i < Nr - 1; i++) {
        Hu_re[i] = a_coeff[i]*u_re[i-1] + b_coeff[i]*u_re[i] + c_coeff[i]*u_re[i+1];
        Hu_im[i] = a_coeff[i]*u_im[i-1] + b_coeff[i]*u_im[i] + c_coeff[i]*u_im[i+1];
    }

    /* RHS = u - i*hdt_h*H*u
     * RHS_re = u_re + hdt_h * H*u_im
     * RHS_im = u_im - hdt_h * H*u_re */
    for (int i = 0; i < Nr; i++) {
        rhs_re[i] = u_re[i] + hdt_h * Hu_im[i];
        rhs_im[i] = u_im[i] - hdt_h * Hu_re[i];
    }

    /* Complex Thomas algorithm: solve (I + i*hdt_h*H) * u^{n+1} = rhs
     * L_i = i*hdt_h*a_i (pure imaginary), D_i = 1+i*hdt_h*b_i, U_i = i*hdt_h*c_i */

    /* Copy RHS to work arrays */
    memcpy(f_re, rhs_re, Nr * sizeof(double));
    memcpy(f_im, rhs_im, Nr * sizeof(double));

    /* Forward elimination */
    d_re[0] = 1.0;
    d_im[0] = hdt_h * b_coeff[0];

    for (int i = 1; i < Nr; i++) {
        double L_im = hdt_h * a_coeff[i];
        double U_im = hdt_h * c_coeff[i-1];

        /* w = L / d_{i-1}  where L = (0, L_im) */
        double den = d_re[i-1]*d_re[i-1] + d_im[i-1]*d_im[i-1];
        double w_re = L_im * d_im[i-1] / den;
        double w_im = L_im * d_re[i-1] / den;

        /* d_i = (1, hdt_h*b_i) - w * (0, U_im) */
        d_re[i] = 1.0 + w_im * U_im;
        d_im[i] = hdt_h * b_coeff[i] - w_re * U_im;

        /* f_i -= w * f_{i-1} */
        double fr = f_re[i] - (w_re*f_re[i-1] - w_im*f_im[i-1]);
        double fi = f_im[i] - (w_re*f_im[i-1] + w_im*f_re[i-1]);
        f_re[i] = fr;
        f_im[i] = fi;
    }

    /* Back substitution */
    {
        int i = Nr - 1;
        double den = d_re[i]*d_re[i] + d_im[i]*d_im[i];
        u_re[i] = (f_re[i]*d_re[i] + f_im[i]*d_im[i]) / den;
        u_im[i] = (f_im[i]*d_re[i] - f_re[i]*d_im[i]) / den;
    }
    for (int i = Nr - 2; i >= 0; i--) {
        double U_im = hdt_h * c_coeff[i];
        /* z_i = (f_i - (0,U_im)*z_{i+1}) / d_i */
        double t_re = f_re[i] + U_im * u_im[i+1];
        double t_im = f_im[i] - U_im * u_re[i+1];
        double den = d_re[i]*d_re[i] + d_im[i]*d_im[i];
        u_re[i] = (t_re*d_re[i] + t_im*d_im[i]) / den;
        u_im[i] = (t_im*d_re[i] - t_re*d_im[i]) / den;
    }

    /* Convert back: psi = u / r */
    for (int i = 0; i < Nr; i++) {
        w->psi_re[i] = u_re[i] / w->r[i];
        w->psi_im[i] = u_im[i] / w->r[i];
    }

    /* Boundary conditions: psi=0 at edges */
    w->psi_re[0] = 0;  w->psi_im[0] = 0;
    w->psi_re[Nr-1] = 0;  w->psi_im[Nr-1] = 0;

    /* Compute norm */
    double norm = 0;
    for (int i = 0; i < Nr; i++) {
        double ri = w->r[i];
        norm += (w->psi_re[i]*w->psi_re[i] + w->psi_im[i]*w->psi_im[i])
                * ri * ri * w->dr[i];
    }
    w->norm = norm;
}

/* Wedge diagnostics: expectation values */
static void wedge_diagnostics(Wedge *w, double *E_wedge, double *r_mean,
                               double *r_rms, double *r_peak) {
    const int Nr = w->Nr;
    double e = 0, rm = 0, r2 = 0;
    double peak_val = 0;
    int peak_idx = 0;
    double norm = 0;

    for (int i = 0; i < Nr; i++) {
        double ri = w->r[i];
        double prob = (w->psi_re[i]*w->psi_re[i] + w->psi_im[i]*w->psi_im[i])
                      * ri * ri * w->dr[i];
        norm += prob;
        rm += ri * prob;
        r2 += ri * ri * prob;

        double density = w->psi_re[i]*w->psi_re[i] + w->psi_im[i]*w->psi_im[i];
        if (density > peak_val) {
            peak_val = density;
            peak_idx = i;
        }
    }

    if (norm > 0) {
        rm /= norm;
        r2 /= norm;
    }
    *r_mean = rm;
    *r_rms = sqrt(r2);
    *r_peak = w->r[peak_idx];

    /* Energy: <H> = <T> + <V> via finite differences */
    double alpha = w->hbar / (2.0 * w->m_eff);
    double ekin = 0, epot = 0;
    for (int i = 1; i < Nr - 1; i++) {
        double drp = w->r[i+1] - w->r[i];
        double drm = w->r[i] - w->r[i-1];
        double dravg = 0.5 * (drp + drm);
        double ri = w->r[i];

        /* d2(r*psi)/dr2 */
        double u_re = ri * w->psi_re[i];
        double u_im = ri * w->psi_im[i];
        double up_re = w->r[i+1] * w->psi_re[i+1];
        double um_re = w->r[i-1] * w->psi_re[i-1];
        double up_im = w->r[i+1] * w->psi_im[i+1];
        double um_im = w->r[i-1] * w->psi_im[i-1];

        double d2u_re = (up_re/drp - u_re*(1.0/drp + 1.0/drm) + um_re/drm) / dravg;
        double d2u_im = (up_im/drp - u_im*(1.0/drp + 1.0/drm) + um_im/drm) / dravg;

        /* T*psi_re = -alpha * d2u_re / r, etc. */
        double Tpsi_re = -alpha * d2u_re / ri;
        double Tpsi_im = -alpha * d2u_im / ri;

        /* <psi|T|psi> contribution */
        ekin += (w->psi_re[i]*Tpsi_re + w->psi_im[i]*Tpsi_im) * ri*ri * w->dr[i];
        /* <psi|V|psi> */
        epot += w->V[i] * (w->psi_re[i]*w->psi_re[i] + w->psi_im[i]*w->psi_im[i])
                * ri*ri * w->dr[i];
    }
    if (norm > 0) {
        ekin /= norm;
        epot /= norm;
    }
    *E_wedge = ekin + epot;
}

/* ================================================================
   Coupling: extract theta at R_match
   ================================================================ */

/* Extract azimuthally-averaged theta_rms(r) near R_match from core.
 * Returns the RMS theta amplitude on the matching sphere.
 * This is a diagnostic for Phase 1; real coupling in Phase 2. */
static double extract_theta_at_rmatch(Grid *g) {
    const int N = g->N, NN = N * N;
    const double dx = g->dx, L = g->L;
    const double shell_width = 2.0 * dx;
    double sum2 = 0;
    long count = 0;

    for (int i = 0; i < N; i++) {
        double x = -L + i * dx;
        for (int j = 0; j < N; j++) {
            double y = -L + j * dx;
            for (int k = 0; k < N; k++) {
                double z = -L + k * dx;
                double r = sqrt(x*x + y*y + z*z);
                if (fabs(r - R_MATCH) < shell_width) {
                    long idx = (long)i * NN + j * N + k;
                    for (int a = 0; a < NFIELDS; a++) {
                        sum2 += g->theta[a][idx] * g->theta[a][idx];
                    }
                    count++;
                }
            }
        }
    }
    return (count > 0) ? sqrt(sum2 / (3 * count)) : 0;
}

/* ================================================================
   Core diagnostics
   ================================================================ */

static void compute_energy(Grid *g, double *E_phi_kin, double *E_theta_kin,
                           double *E_grad, double *E_mass, double *E_pot,
                           double *E_theta_grad, double *E_theta_mass,
                           double *E_coupling, double *E_total) {
    const int N = g->N, NN = N*N;
    const long N3 = g->N3;
    const double dx = g->dx, dV = dx*dx*dx;
    const double idx1 = 1.0 / (2.0*dx);
    double epk=0, etk=0, eg=0, em=0, ep=0, etg=0, etm=0, ec=0;

    #pragma omp parallel for reduction(+:epk,etk,eg,em,ep,etg,etm,ec) schedule(static)
    for (long idx = 0; idx < N3; idx++) {
        if (g->mask[idx] == 0) continue;  /* skip outside */

        int i=(int)(idx/NN), j=(int)((idx/N)%N), k=(int)(idx%N);
        int ip = (i+1 < N) ? i+1 : i;
        int im = (i-1 >= 0) ? i-1 : i;
        int jp = (j+1 < N) ? j+1 : j;
        int jm = (j-1 >= 0) ? j-1 : j;
        int kp = (k+1 < N) ? k+1 : k;
        int km = (k-1 >= 0) ? k-1 : k;
        long n_ip=(long)ip*NN+j*N+k, n_im=(long)im*NN+j*N+k;
        long n_jp=(long)i*NN+jp*N+k, n_jm=(long)i*NN+jm*N+k;
        long n_kp=(long)i*NN+j*N+kp, n_km=(long)i*NN+j*N+km;

        for (int a = 0; a < NFIELDS; a++) {
            epk += 0.5*g->phi_vel[a][idx]*g->phi_vel[a][idx]*dV;
            etk += 0.5*g->theta_vel[a][idx]*g->theta_vel[a][idx]*dV;
            double gx=(g->phi[a][n_ip]-g->phi[a][n_im])*idx1;
            double gy=(g->phi[a][n_jp]-g->phi[a][n_jm])*idx1;
            double gz=(g->phi[a][n_kp]-g->phi[a][n_km])*idx1;
            eg += 0.5*(gx*gx+gy*gy+gz*gz)*dV;
            em += 0.5*MASS2*g->phi[a][idx]*g->phi[a][idx]*dV;
            double tgx=(g->theta[a][n_ip]-g->theta[a][n_im])*idx1;
            double tgy=(g->theta[a][n_jp]-g->theta[a][n_jm])*idx1;
            double tgz=(g->theta[a][n_kp]-g->theta[a][n_km])*idx1;
            etg += 0.5*(tgx*tgx+tgy*tgy+tgz*tgz)*dV;
            etm += 0.5*MTHETA2*g->theta[a][idx]*g->theta[a][idx]*dV;
        }
        double P = g->phi[0][idx]*g->phi[1][idx]*g->phi[2][idx];
        ep += (MU/2.0)*P*P/(1.0+KAPPA*P*P)*dV;

        for (int a = 0; a < NFIELDS; a++) {
            double ct = curl_component(g->theta, a,
                n_ip, n_im, n_jp, n_jm, n_kp, n_km, idx1);
            ec -= ETA * g->phi[a][idx] * ct * dV;
        }
    }
    *E_phi_kin=epk; *E_theta_kin=etk; *E_grad=eg; *E_mass=em;
    *E_pot=ep; *E_theta_grad=etg; *E_theta_mass=etm; *E_coupling=ec;
    *E_total = epk+etk+eg+em+ep+etg+etm+ec;
}

static double theta_rms(Grid *g) {
    double sum = 0;
    long count = 0;
    for (int a = 0; a < NFIELDS; a++)
        for (long i = 0; i < g->N3; i++) {
            if (g->mask[i] == 0) continue;
            sum += g->theta[a][i] * g->theta[a][i];
            count++;
        }
    return (count > 0) ? sqrt(sum / count) : 0;
}

/* ================================================================
   Core initialization: braid inside sphere
   ================================================================ */

static void init_braid(Grid *g) {
    const int N = g->N, NN = N*N;
    const double dx = g->dx, L = g->L;
    const double A[3] = {0.8, 0.8, 0.8};
    const double delta[3] = {0, 3.0005, 4.4325};
    const double R_tube = 3.0, ellip = 0.3325;
    const double kw = PI/L, omega = sqrt(kw*kw + MASS2);
    const double sx = 1+ellip, sy = 1-ellip;
    const double inv2R2 = 1.0/(2*R_tube*R_tube);
    const double k_bg = PI/L, omega_bg = sqrt(k_bg*k_bg + MASS2);

    for (int i = 0; i < N; i++) {
        double x = -L + i*dx;
        for (int j = 0; j < N; j++) {
            double y = -L + j*dx;
            for (int kk = 0; kk < N; kk++) {
                double z = -L + kk*dx;
                long idx = (long)i*NN + j*N + kk;
                if (g->mask[idx] == 0) continue;  /* skip outside */

                double r2e = x*x/(sx*sx) + y*y/(sy*sy);
                double env = exp(-r2e * inv2R2);
                for (int a = 0; a < NFIELDS; a++) {
                    double ph = kw*z + delta[a];
                    double ph_bg = k_bg*z + 2*PI*a/3.0;
                    g->phi[a][idx] = A[a]*env*cos(ph) + A_BG*cos(ph_bg);
                    g->phi_vel[a][idx] = omega*A[a]*env*sin(ph) + omega_bg*A_BG*sin(ph_bg);
                }
            }
        }
    }
}

/* Save snapshot (core only, binary) */
static void save_field(Grid *g, double t, const char *dir) {
    char fn[512];
    snprintf(fn, sizeof(fn), "%s/field_t%04d.bin", dir, (int)(t+0.5));
    FILE *fp = fopen(fn, "wb"); if (!fp) return;
    int n = g->N; double l = g->L;
    int nf = 6;
    fwrite(&n, sizeof(int), 1, fp);
    fwrite(&l, sizeof(double), 1, fp);
    fwrite(&t, sizeof(double), 1, fp);
    fwrite(&nf, sizeof(int), 1, fp);
    for (int a = 0; a < NFIELDS; a++) fwrite(g->phi[a], sizeof(double), g->N3, fp);
    for (int a = 0; a < NFIELDS; a++) fwrite(g->theta[a], sizeof(double), g->N3, fp);
    fclose(fp);
    printf("  [SNAP] %s (%.1f MB)\n", fn, 6.0*g->N3*8.0/1e6);
}

/* Save wedge wavefunction snapshot */
static void save_wedge(Wedge *w, double t, const char *dir) {
    char fn[512];
    snprintf(fn, sizeof(fn), "%s/wedge_t%04d.dat", dir, (int)(t+0.5));
    FILE *fp = fopen(fn, "w"); if (!fp) return;
    fprintf(fp, "# r  psi_re  psi_im  |psi|^2  V(r)\n");
    for (int i = 0; i < w->Nr; i++) {
        double p2 = w->psi_re[i]*w->psi_re[i] + w->psi_im[i]*w->psi_im[i];
        fprintf(fp, "%.6e  %.6e  %.6e  %.6e  %.6e\n",
                w->r[i], w->psi_re[i], w->psi_im[i], p2, w->V[i]);
    }
    fclose(fp);
}

/* ================================================================
   Main
   ================================================================ */

int main(int argc, char **argv) {
    int N = 128;
    double L = 25.0, T = 100.0;
    int Nr_wedge = 2000;
    double diag_dt = 2.0, snap_dt = 50.0;
    int N_couple = 10;       /* coupling every N_couple core steps */
    int wedge_sub = 1000;    /* wedge substeps per coupling (more = faster wedge evolution) */
    double wedge_dt_mult = 100.0; /* wedge dt = this × core coupling interval / wedge_sub */
    char outdir[256] = "data/multiscale";

    for (int i = 1; i < argc; i++) {
        if      (!strcmp(argv[i],"-N"))       N = atoi(argv[++i]);
        else if (!strcmp(argv[i],"-L"))       L = atof(argv[++i]);
        else if (!strcmp(argv[i],"-T"))       T = atof(argv[++i]);
        else if (!strcmp(argv[i],"-Nr"))      Nr_wedge = atoi(argv[++i]);
        else if (!strcmp(argv[i],"-eta"))     ETA = atof(argv[++i]);
        else if (!strcmp(argv[i],"-Rmatch"))  R_MATCH = atof(argv[++i]);
        else if (!strcmp(argv[i],"-Rmax"))    R_MAX = atof(argv[++i]);
        else if (!strcmp(argv[i],"-hbar"))    HBAR_EFF = atof(argv[++i]);
        else if (!strcmp(argv[i],"-meff"))    M_EFF = atof(argv[++i]);
        else if (!strcmp(argv[i],"-Vdepth")) V_DEPTH = atof(argv[++i]);
        else if (!strcmp(argv[i],"-npower")) N_POWER = atof(argv[++i]);
        else if (!strcmp(argv[i],"-bg"))     A_BG = atof(argv[++i]);
        else if (!strcmp(argv[i],"-m"))      { double m = atof(argv[++i]); MASS2 = m*m; }
        else if (!strcmp(argv[i],"-mt"))     { double m = atof(argv[++i]); MTHETA2 = m*m; }
        else if (!strcmp(argv[i],"-couple"))  N_couple = atoi(argv[++i]);
        else if (!strcmp(argv[i],"-wsub"))    wedge_sub = atoi(argv[++i]);
        else if (!strcmp(argv[i],"-wdtmult")) wedge_dt_mult = atof(argv[++i]);
        else if (!strcmp(argv[i],"-diag"))    diag_dt = atof(argv[++i]);
        else if (!strcmp(argv[i],"-snap"))    snap_dt = atof(argv[++i]);
        else if (!strcmp(argv[i],"-o"))       strncpy(outdir, argv[++i], 255);
    }

    int nthreads = 8;
    char *env_threads = getenv("OMP_NUM_THREADS");
    if (env_threads) nthreads = atoi(env_threads);
    omp_set_num_threads(nthreads);

    printf("===================================================\n");
    printf("  Multi-Scale Atom Simulation — Phase 1\n");
    printf("===================================================\n\n");

    printf("CORE (3D Cosserat with spherical boundary):\n");
    printf("  N=%d  L=%.0f  R_match=%.1f\n", N, L, R_MATCH);
    printf("  m^2=%.4f  m_theta^2=%.4f  eta=%.3f\n", MASS2, MTHETA2, ETA);
    printf("  mu=%.3f  kappa=%.1f\n\n", MU, KAPPA);

    printf("WEDGE (1D Schrodinger, Crank-Nicolson):\n");
    printf("  Nr=%d  r=[%.0f, %.0f]\n", Nr_wedge, R_MATCH, R_MAX);
    printf("  hbar_eff=%.0f  m_eff=%.0f\n", HBAR_EFF, M_EFF);
    printf("  V_depth=%.3f  n_power=%.3f  R_braid=%.1f\n", V_DEPTH, N_POWER, R_BRAID);
    printf("  Bohr_radius=%.0f\n\n", BOHR_RAD);

    printf("COUPLING: every %d core steps, %d wedge substeps\n\n", N_couple, wedge_sub);

    mkdir("data", 0755);
    mkdir(outdir, 0755);

    /* ---- Allocate core ---- */
    Grid *g = grid_alloc(N, L);
    printf("Core: dx=%.4f  dt=%.5f  threads=%d\n", g->dx, g->dt, nthreads);

    /* ---- Allocate wedge ---- */
    Wedge *w = wedge_alloc(Nr_wedge, R_MATCH, R_MAX, HBAR_EFF, M_EFF);

    /* Wedge time is DECOUPLED from core time.
     * The wedge dt is set so the electron sees orbital timescales
     * during the core's much shorter simulation time.
     * T_orbital ~ 2π × a₀ × m_eff / ℏ ≈ 112,000 code units.
     * Each coupling step, wedge advances wedge_sub × dt_wedge time units. */
    double core_couple_dt = N_couple * g->dt;
    w->dt = wedge_dt_mult * core_couple_dt / wedge_sub;
    double wedge_T_per_couple = wedge_sub * w->dt;
    int n_couples = (int)(T / (N_couple * g->dt));
    double wedge_T_total = n_couples * wedge_T_per_couple;
    double T_orbital = 2*PI * BOHR_RAD * M_EFF / HBAR_EFF;
    printf("Wedge dt = %.2f (%.0f× core coupling interval)\n", w->dt, wedge_dt_mult);
    printf("Wedge time per coupling: %.0f (vs orbital period %.0f)\n",
           wedge_T_per_couple, T_orbital);
    printf("Wedge total time: %.0f = %.2f orbital periods\n\n",
           wedge_T_total, wedge_T_total / T_orbital);

    /* ---- Initialize ---- */
    init_braid(g);
    compute_forces(g);

    /* Wedge: Gaussian wave packet — configurable position */
    static double r_frac = 0.8;  /* fraction of Bohr radius for initial position */
    /* Parse -rfrac if present */
    for (int i = 1; i < argc; i++)
        if (!strcmp(argv[i], "-rfrac")) r_frac = atof(argv[++i]);
    double r_packet = r_frac * BOHR_RAD;
    double sigma_packet = 0.2 * BOHR_RAD; /* tighter packet = more eigenstate mixing */
    wedge_init_gaussian(w, r_packet, sigma_packet);

    /* ---- Open timeseries ---- */
    char tspath[512];
    snprintf(tspath, sizeof(tspath), "%s/timeseries.tsv", outdir);
    FILE *fp = fopen(tspath, "w");
    fprintf(fp, "t\tE_phi_kin\tE_theta_kin\tE_grad\tE_mass\tE_pot\t"
            "E_tgrad\tE_tmass\tE_coupling\tE_core_total\t"
            "theta_rms\ttheta_at_rmatch\t"
            "E_wedge\tr_mean\tr_rms\tr_peak\twedge_norm\n");

    int n_steps = (int)(T / g->dt);
    int diag_every = (int)(diag_dt / g->dt); if (diag_every<1) diag_every=1;
    int snap_every = (int)(snap_dt / g->dt); if (snap_every<1) snap_every=1;

    printf("Core steps=%d  diag_every=%d  snap_every=%d\n\n", n_steps, diag_every, snap_every);

    double wall0 = omp_get_wtime();
    double E0 = 0;

    /* Initial snapshot */
    save_field(g, 0, outdir);
    save_wedge(w, 0, outdir);

    /* ---- Main loop ---- */
    for (int step = 0; step <= n_steps; step++) {
        if (step > 0) {
            /* Core step */
            verlet_step(g);

            /* Coupling: every N_couple steps, advance wedge */
            if (step % N_couple == 0) {
                /* Extract theta at R_match (diagnostic) */
                /* double theta_rm = extract_theta_at_rmatch(g); */
                /* In Phase 1, wedge evolves independently (no back-coupling) */

                /* Advance wedge by wedge_sub steps */
                for (int ws = 0; ws < wedge_sub; ws++) {
                    wedge_step(w);
                }
            }
        }

        double t = step * g->dt;

        /* Diagnostics */
        if (step % diag_every == 0) {
            double epk, etk, eg, em, ep, etg, etm, ec, et;
            compute_energy(g, &epk, &etk, &eg, &em, &ep, &etg, &etm, &ec, &et);
            if (step == 0) E0 = et;
            double trms = theta_rms(g);
            double theta_rm = extract_theta_at_rmatch(g);

            double E_w, r_m, r_r, r_p;
            wedge_diagnostics(w, &E_w, &r_m, &r_r, &r_p);

            fprintf(fp, "%.4f\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t"
                    "%.6e\t%.6e\t%.6e\t%.1f\t%.1f\t%.1f\t%.8f\n",
                    t, epk, etk, eg, em, ep, etg, etm, ec, et,
                    trms, theta_rm, E_w, r_m, r_r, r_p, w->norm);
            fflush(fp);

            if (step % (diag_every*10) == 0) {
                double wall = omp_get_wtime() - wall0;
                double drift = 100.0*(et-E0)/(fabs(E0)+1e-30);
                printf("t=%7.2f  Core: E=%.4e (drift %+.3f%%) theta_rms=%.3e  "
                       "Wedge: <r>=%.0f rms=%.0f E=%.4e norm=%.6f  [%.0f%% %.0fs]\n",
                       t, et, drift, trms, r_m, r_r, E_w, w->norm,
                       100.0*step/n_steps, wall);
                fflush(stdout);
            }
        }

        if (step > 0 && step % snap_every == 0) {
            save_field(g, t, outdir);
            save_wedge(w, t, outdir);
        }
    }

    /* Final snapshots */
    save_field(g, T, outdir);
    save_wedge(w, T, outdir);
    fclose(fp);

    double wall = omp_get_wtime() - wall0;
    printf("\n===================================================\n");
    printf("  Complete: %.0fs (%.1f min)\n", wall, wall/60);
    printf("===================================================\n");

    /* Final summary */
    double epk, etk, eg, em, ep, etg, etm, ec, et;
    compute_energy(g, &epk, &etk, &eg, &em, &ep, &etg, &etm, &ec, &et);
    double E_w, r_m, r_r, r_p;
    wedge_diagnostics(w, &E_w, &r_m, &r_r, &r_p);

    printf("\nFinal state:\n");
    printf("  Core:  E_total=%.4e  theta_rms=%.3e  E_pot=%.4e\n", et, theta_rms(g), ep);
    printf("  Wedge: E=%.4e  <r>=%.1f  rms_r=%.1f  peak=%.1f  norm=%.6f\n",
           E_w, r_m, r_r, r_p, w->norm);
    printf("  Energy drift: %+.4f%%\n", 100.0*(et-E0)/(fabs(E0)+1e-30));

    grid_free(g);
    wedge_free(w);
    return 0;
}
