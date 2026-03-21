/*  v37_knot.c — Knot-geometry 6-field Cosserat simulation
 *
 *  Geometries:
 *    trefoil_R  — right-handed (2,3) torus knot
 *    trefoil_L  — left-handed (2,3) torus knot (mirror of trefoil_R)
 *    figure8    — figure-eight knot (4_1, amphichiral)
 *
 *  Initialization: one continuous tube following the knot curve, all three
 *  fields carry phase offsets along the arc length.
 *
 *  Equations (same as v37_crossed_v2):
 *    d^2 phi_a/dt^2 = Lap(phi_a) - m^2 phi_a - V'(P) + eta * curl(theta)_a
 *    d^2 theta_a/dt^2 = Lap(theta_a) - m_theta^2 theta_a + eta * curl(phi)_a
 *
 *  V(P) = (mu/2) P^2 / (1 + kappa P^2),  P = phi_0 phi_1 phi_2
 *
 *  Build: gcc -O3 -march=native -fopenmp -o v37_knot src/v37_knot.c -lzstd -lm
 *  Run:   ./v37_knot -geom trefoil_R -N 128 -L 25 -T 300 -snap 5 -o data/trefoil_R
 */

#define SFA_IMPLEMENTATION
#include "../../sfa/format/sfa.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>
#include <sys/stat.h>

#define NFIELDS 3
#define PI 3.14159265358979323846
#define N_CURVE 1000

/* ================================================================
   Geometry type
   ================================================================ */

enum { GEOM_TREFOIL_R = 0, GEOM_TREFOIL_L = 1, GEOM_FIGURE8 = 2 };

static int    GEOM_TYPE  = -1;  /* must be set via -geom */
static char   GEOM_NAME[32] = "";

/* Physics parameters */
static double MU      = -41.345;
static double KAPPA   = 50.0;
static double MASS2   = 2.25;     /* position field mass^2 */
static double MTHETA2 = 0.0;      /* angle field mass^2 (default massless) */
static double A_BG    = 0.1;
static double ETA     = 0.5;      /* coupling strength */

/* Knot geometry parameters */
static double R_MAJOR  = 3.0;     /* torus major radius (trefoil) */
static double R_MINOR  = 1.5;     /* torus minor radius (trefoil) */
static double F8_SCALE = 1.5;     /* figure-eight scale */
static double TUBE_R   = 2.0;     /* tube radius */
static double A_AMP    = 0.8;     /* field amplitude */
static int    N_OSC    = -1;      /* oscillations around loop (-1 = auto) */

/* Phase offsets (from v28 optimization) */
static double DELTA[3] = {0.0, 3.0005, 4.4325};

/* Absorbing boundary (gentler defaults than v2) */
static double DAMP_WIDTH = 4.0;
static double DAMP_RATE  = 0.005;

/* Death check: rolling average window */
#define DEATH_WINDOW 128
static double ep_history[DEATH_WINDOW];
static int ep_idx = 0, ep_count = 0;

/* ================================================================
   Grid
   ================================================================ */

typedef struct {
    double *mem;
    double *phi[NFIELDS];
    double *phi_vel[NFIELDS];
    double *phi_acc[NFIELDS];
    double *theta[NFIELDS];
    double *theta_vel[NFIELDS];
    double *theta_acc[NFIELDS];
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
    printf("Allocating %.2f GB (%ld doubles, N=%d, 6 fields)\n", bytes/1e9, total, N);
    g->mem = malloc(total * sizeof(double));
    if (!g->mem) { fprintf(stderr, "FATAL: malloc failed\n"); exit(1); }
    memset(g->mem, 0, total * sizeof(double));

    for (int a = 0; a < NFIELDS; a++) {
        g->phi[a]       = g->mem + (0  + a) * g->N3;
        g->phi_vel[a]   = g->mem + (3  + a) * g->N3;
        g->phi_acc[a]   = g->mem + (6  + a) * g->N3;
        g->theta[a]     = g->mem + (9  + a) * g->N3;
        g->theta_vel[a] = g->mem + (12 + a) * g->N3;
        g->theta_acc[a] = g->mem + (15 + a) * g->N3;
    }
    return g;
}

static void grid_free(Grid *g) { free(g->mem); free(g); }

/* ================================================================
   Curl helpers
   ================================================================ */

static inline double curl_component(double *F[3], int a,
    long n_ip, long n_im, long n_jp, long n_jm, long n_kp, long n_km,
    double idx1) {
    if (a == 0) return (F[2][n_jp] - F[2][n_jm] - F[1][n_kp] + F[1][n_km]) * idx1;
    if (a == 1) return (F[0][n_kp] - F[0][n_km] - F[2][n_ip] + F[2][n_im]) * idx1;
    return            (F[1][n_ip] - F[1][n_im] - F[0][n_jp] + F[0][n_jm]) * idx1;
}

/* ================================================================
   Forces
   ================================================================ */

static void compute_forces(Grid *g) {
    const int N = g->N, NN = N * N;
    const long N3 = g->N3;
    const double idx2 = 1.0 / (g->dx * g->dx);
    const double idx1 = 1.0 / (2.0 * g->dx);

    #pragma omp parallel for schedule(static)
    for (long idx = 0; idx < N3; idx++) {
        int i = (int)(idx / NN);
        int j = (int)((idx / N) % N);
        int k = (int)(idx % N);

        int ip = (i+1)%N, im = (i-1+N)%N;
        int jp = (j+1)%N, jm = (j-1+N)%N;
        int kp = (k+1)%N, km = (k-1+N)%N;

        long n_ip = (long)ip*NN + j*N + k;
        long n_im = (long)im*NN + j*N + k;
        long n_jp = (long)i*NN + jp*N + k;
        long n_jm = (long)i*NN + jm*N + k;
        long n_kp = (long)i*NN + j*N + kp;
        long n_km = (long)i*NN + j*N + km;

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
   Absorbing boundary damping
   ================================================================ */

static void apply_boundary_damping(Grid *g) {
    const int N = g->N, NN = N * N;
    const double L = g->L, dx = g->dx;

    #pragma omp parallel for schedule(static)
    for (long idx = 0; idx < g->N3; idx++) {
        int i = (int)(idx / (NN));
        int j = (int)((idx / N) % N);
        int k = (int)(idx % N);
        double x = -L + i*dx, y = -L + j*dx, z = -L + k*dx;

        /* Distance from nearest edge */
        double dx_edge = fmin(fabs(x + L), fabs(x - L));
        double dy_edge = fmin(fabs(y + L), fabs(y - L));
        double dz_edge = fmin(fabs(z + L), fabs(z - L));
        double d_edge = fmin(fmin(dx_edge, dy_edge), dz_edge);

        if (d_edge < DAMP_WIDTH) {
            double s = (DAMP_WIDTH - d_edge) / DAMP_WIDTH;  /* 0 at inner edge, 1 at boundary */
            double damp = 1.0 - DAMP_RATE * s * s;  /* quadratic profile */
            for (int a = 0; a < NFIELDS; a++) {
                g->phi_vel[a][idx] *= damp;
                g->theta_vel[a][idx] *= damp;
            }
        }
    }
}

/* ================================================================
   Verlet integrator (with absorbing boundary)
   ================================================================ */

static void verlet_step(Grid *g) {
    const long N3 = g->N3;
    const double hdt = 0.5 * g->dt, dt = g->dt;

    for (int a = 0; a < NFIELDS; a++) {
        double *vp = g->phi_vel[a], *ap = g->phi_acc[a];
        double *vt = g->theta_vel[a], *at = g->theta_acc[a];
        for (long i = 0; i < N3; i++) { vp[i] += hdt*ap[i]; vt[i] += hdt*at[i]; }
    }
    for (int a = 0; a < NFIELDS; a++) {
        double *pp = g->phi[a], *vp = g->phi_vel[a];
        double *pt = g->theta[a], *vt = g->theta_vel[a];
        for (long i = 0; i < N3; i++) { pp[i] += dt*vp[i]; pt[i] += dt*vt[i]; }
    }
    compute_forces(g);
    for (int a = 0; a < NFIELDS; a++) {
        double *vp = g->phi_vel[a], *ap = g->phi_acc[a];
        double *vt = g->theta_vel[a], *at = g->theta_acc[a];
        for (long i = 0; i < N3; i++) { vp[i] += hdt*ap[i]; vt[i] += hdt*at[i]; }
    }

    /* Apply absorbing boundary damping after each step */
    apply_boundary_damping(g);
}

/* ================================================================
   Knot curve generation
   ================================================================ */

typedef struct {
    double x[N_CURVE];
    double y[N_CURVE];
    double z[N_CURVE];
    double arclen[N_CURVE];  /* cumulative arc length at each sample */
    double L_curve;          /* total curve length (closed loop) */
    double bbox_min[3], bbox_max[3];
} KnotCurve;

static void knot_curve_trefoil(KnotCurve *kc, double R, double r, int left_handed) {
    for (int i = 0; i < N_CURVE; i++) {
        double t = 2.0 * PI * i / N_CURVE;
        kc->x[i] = (R + r * cos(3.0 * t)) * cos(2.0 * t);
        kc->y[i] = (R + r * cos(3.0 * t)) * sin(2.0 * t);
        kc->z[i] = (left_handed ? -1.0 : 1.0) * r * sin(3.0 * t);
    }
}

static void knot_curve_figure8(KnotCurve *kc, double scale) {
    for (int i = 0; i < N_CURVE; i++) {
        double t = 2.0 * PI * i / N_CURVE;
        kc->x[i] = scale * (2.0 + cos(2.0 * t)) * cos(3.0 * t);
        kc->y[i] = scale * (2.0 + cos(2.0 * t)) * sin(3.0 * t);
        kc->z[i] = scale * sin(4.0 * t);
    }
}

static void knot_compute_arclength(KnotCurve *kc) {
    kc->arclen[0] = 0.0;
    for (int i = 1; i < N_CURVE; i++) {
        double dx = kc->x[i] - kc->x[i-1];
        double dy = kc->y[i] - kc->y[i-1];
        double dz = kc->z[i] - kc->z[i-1];
        kc->arclen[i] = kc->arclen[i-1] + sqrt(dx*dx + dy*dy + dz*dz);
    }
    /* Close the loop */
    double dx = kc->x[0] - kc->x[N_CURVE-1];
    double dy = kc->y[0] - kc->y[N_CURVE-1];
    double dz = kc->z[0] - kc->z[N_CURVE-1];
    kc->L_curve = kc->arclen[N_CURVE-1] + sqrt(dx*dx + dy*dy + dz*dz);

    /* Bounding box */
    kc->bbox_min[0] = kc->bbox_max[0] = kc->x[0];
    kc->bbox_min[1] = kc->bbox_max[1] = kc->y[0];
    kc->bbox_min[2] = kc->bbox_max[2] = kc->z[0];
    for (int i = 1; i < N_CURVE; i++) {
        if (kc->x[i] < kc->bbox_min[0]) kc->bbox_min[0] = kc->x[i];
        if (kc->x[i] > kc->bbox_max[0]) kc->bbox_max[0] = kc->x[i];
        if (kc->y[i] < kc->bbox_min[1]) kc->bbox_min[1] = kc->y[i];
        if (kc->y[i] > kc->bbox_max[1]) kc->bbox_max[1] = kc->y[i];
        if (kc->z[i] < kc->bbox_min[2]) kc->bbox_min[2] = kc->z[i];
        if (kc->z[i] > kc->bbox_max[2]) kc->bbox_max[2] = kc->z[i];
    }
}

/* Count approximate self-crossings in the knot projection (xy plane) */
static int knot_count_crossings(KnotCurve *kc) {
    int crossings = 0;
    /* Check all pairs of segments for intersection in xy projection */
    for (int i = 0; i < N_CURVE; i++) {
        int i1 = (i + 1) % N_CURVE;
        double ax = kc->x[i], ay = kc->y[i];
        double bx = kc->x[i1], by = kc->y[i1];

        for (int j = i + 2; j < N_CURVE; j++) {
            if (i == 0 && j == N_CURVE - 1) continue;  /* skip adjacent */
            int j1 = (j + 1) % N_CURVE;
            double cx = kc->x[j], cy = kc->y[j];
            double dx = kc->x[j1], dy = kc->y[j1];

            /* 2D segment intersection test */
            double denom = (bx - ax) * (dy - cy) - (by - ay) * (dx - cx);
            if (fabs(denom) < 1e-12) continue;
            double t = ((cx - ax) * (dy - cy) - (cy - ay) * (dx - cx)) / denom;
            double u = ((cx - ax) * (by - ay) - (cy - ay) * (bx - ax)) / denom;
            if (t > 0.0 && t < 1.0 && u > 0.0 && u < 1.0)
                crossings++;
        }
    }
    return crossings;
}

/* ================================================================
   Initialization: Knot tube
   ================================================================ */

static void init_knot(Grid *g, KnotCurve *kc) {
    const int N = g->N, NN = N*N;
    const double dx = g->dx, L = g->L;
    const double tube_r = TUBE_R;
    const double inv2t2 = 1.0 / (2.0 * tube_r * tube_r);
    const double k_knot = 2.0 * PI * N_OSC / kc->L_curve;
    const double k_bg = PI / L;
    const double omega_bg = sqrt(k_bg * k_bg + MASS2);

    printf("  Init: knot tube (geom=%s)\n", GEOM_NAME);
    printf("  tube_r=%.1f, A=%.2f, n_osc=%d\n", tube_r, A_AMP, N_OSC);
    printf("  k_knot=%.4f (L_curve=%.3f)\n", k_knot, kc->L_curve);
    printf("  delta = {%.4f, %.4f, %.4f}\n", DELTA[0], DELTA[1], DELTA[2]);
    printf("  Background: per-axis directional, A_bg=%.2f\n", A_BG);
    printf("  Absorbing boundary: width=%.1f rate=%.4f\n", DAMP_WIDTH, DAMP_RATE);

    /* Pre-compute segment lengths for arc-length interpolation on closing segment */
    double seg_len[N_CURVE];
    for (int i = 0; i < N_CURVE - 1; i++) {
        double ddx = kc->x[i+1] - kc->x[i];
        double ddy = kc->y[i+1] - kc->y[i];
        double ddz = kc->z[i+1] - kc->z[i];
        seg_len[i] = sqrt(ddx*ddx + ddy*ddy + ddz*ddz);
    }
    {
        double ddx = kc->x[0] - kc->x[N_CURVE-1];
        double ddy = kc->y[0] - kc->y[N_CURVE-1];
        double ddz = kc->z[0] - kc->z[N_CURVE-1];
        seg_len[N_CURVE-1] = sqrt(ddx*ddx + ddy*ddy + ddz*ddz);
    }

    #pragma omp parallel for schedule(static) collapse(3)
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            for (int kk = 0; kk < N; kk++) {
                double px = -L + i * dx;
                double py = -L + j * dx;
                double pz = -L + kk * dx;
                long idx = (long)i * NN + j * N + kk;

                /* Find nearest point on the knot curve */
                double d_min = 1e30;
                double s_min = 0.0;

                for (int c = 0; c < N_CURVE; c++) {
                    int c1 = (c + 1) % N_CURVE;
                    double ax = kc->x[c],  ay = kc->y[c],  az = kc->z[c];
                    double bx = kc->x[c1], by = kc->y[c1], bz = kc->z[c1];

                    /* Vector AB and AP */
                    double abx = bx - ax, aby = by - ay, abz = bz - az;
                    double apx = px - ax, apy = py - ay, apz = pz - az;

                    double ab2 = abx*abx + aby*aby + abz*abz;
                    double t;
                    if (ab2 < 1e-30) {
                        t = 0.0;
                    } else {
                        t = (apx*abx + apy*aby + apz*abz) / ab2;
                        if (t < 0.0) t = 0.0;
                        if (t > 1.0) t = 1.0;
                    }

                    /* Closest point on segment */
                    double cx = ax + t * abx;
                    double cy = ay + t * aby;
                    double cz = az + t * abz;
                    double ddx = px - cx, ddy = py - cy, ddz = pz - cz;
                    double dist = sqrt(ddx*ddx + ddy*ddy + ddz*ddz);

                    if (dist < d_min) {
                        d_min = dist;
                        s_min = kc->arclen[c] + t * seg_len[c];
                    }
                }

                /* Tube envelope */
                double tube_env = exp(-d_min * d_min * inv2t2);

                for (int a = 0; a < NFIELDS; a++) {
                    double phase = k_knot * s_min + DELTA[a];
                    double val = A_AMP * tube_env * cos(phase);
                    double vel = 0.0;  /* stationary start */

                    /* Per-field directional background */
                    double bg_dir;
                    if (a == 0) bg_dir = pz;
                    else if (a == 1) bg_dir = px;
                    else bg_dir = py;
                    double ph_bg = k_bg * bg_dir + 2.0 * PI * a / 3.0;
                    val += A_BG * cos(ph_bg);
                    vel += omega_bg * A_BG * sin(ph_bg);

                    g->phi[a][idx] = val;
                    g->phi_vel[a][idx] = vel;
                }
            }
        }
    }
}

/* ================================================================
   Diagnostics (same as v37_crossed_v2)
   ================================================================ */

static void compute_energy(Grid *g, double *E_phi_kin, double *E_theta_kin,
                           double *E_grad, double *E_mass, double *E_pot,
                           double *E_theta_grad, double *E_theta_mass,
                           double *E_coupling, double *E_total) {
    const int N = g->N, NN = N*N;
    const long N3 = g->N3;
    const double dx = g->dx, dV = dx*dx*dx;
    const double idx1 = 1.0/(2.0*dx);
    double epk=0, etk=0, eg=0, em=0, ep=0, etg=0, etm=0, ec=0;

    #pragma omp parallel for reduction(+:epk,etk,eg,em,ep,etg,etm,ec) schedule(static)
    for (long idx = 0; idx < N3; idx++) {
        int i=(int)(idx/NN), j=(int)((idx/N)%N), k=(int)(idx%N);
        int ip=(i+1)%N, im=(i-1+N)%N;
        int jp=(j+1)%N, jm=(j-1+N)%N;
        int kp=(k+1)%N, km=(k-1+N)%N;
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
    for (int a = 0; a < NFIELDS; a++)
        for (long i = 0; i < g->N3; i++)
            sum += g->theta[a][i] * g->theta[a][i];
    return sqrt(sum / (3 * g->N3));
}

/* Triple product |P| integrated (measures knot survival) */
static double P_integrated(Grid *g) {
    const long N3 = g->N3;
    const double dV = g->dx * g->dx * g->dx;
    double total = 0;
    #pragma omp parallel for reduction(+:total)
    for (long i = 0; i < N3; i++) {
        double P = g->phi[0][i] * g->phi[1][i] * g->phi[2][i];
        total += fabs(P) * dV;
    }
    return total;
}

/* Moment of inertia tensor from field amplitude */
static void compute_inertia_tensor(Grid *g, double Ixx[3][3], double *E_tot_out) {
    const int N = g->N, NN = N*N;
    const long N3 = g->N3;
    const double dx = g->dx, L = g->L;
    const double dV = dx*dx*dx;

    double cm0=0, cm1=0, cm2=0, E_sum = 0;

    #pragma omp parallel for reduction(+:cm0,cm1,cm2,E_sum) schedule(static)
    for (long idx = 0; idx < N3; idx++) {
        int i=(int)(idx/NN), j=(int)((idx/N)%N), k=(int)(idx%N);
        double x = -L + i*dx, y = -L + j*dx, z = -L + k*dx;

        double w = 0;
        for (int a = 0; a < NFIELDS; a++)
            w += g->phi[a][idx] * g->phi[a][idx];
        w *= dV;

        cm0 += x * w; cm1 += y * w; cm2 += z * w;
        E_sum += w;
    }
    if (E_sum > 0) { cm0/=E_sum; cm1/=E_sum; cm2/=E_sum; }
    *E_tot_out = E_sum;

    double I00=0,I01=0,I02=0,I11=0,I12=0,I22=0;
    #pragma omp parallel for reduction(+:I00,I01,I02,I11,I12,I22) schedule(static)
    for (long idx = 0; idx < N3; idx++) {
        int i=(int)(idx/NN), j=(int)((idx/N)%N), k=(int)(idx%N);
        double x = -L + i*dx - cm0;
        double y = -L + j*dx - cm1;
        double z = -L + k*dx - cm2;

        double w = 0;
        for (int a = 0; a < NFIELDS; a++)
            w += g->phi[a][idx] * g->phi[a][idx];
        w *= dV;

        double r2 = x*x + y*y + z*z;
        I00 += (r2 - x*x) * w;
        I11 += (r2 - y*y) * w;
        I22 += (r2 - z*z) * w;
        I01 -= x*y * w;
        I02 -= x*z * w;
        I12 -= y*z * w;
    }
    Ixx[0][0]=I00; Ixx[0][1]=I01; Ixx[0][2]=I02;
    Ixx[1][0]=I01; Ixx[1][1]=I11; Ixx[1][2]=I12;
    Ixx[2][0]=I02; Ixx[2][1]=I12; Ixx[2][2]=I22;
}

/* Eigenvalues of 3x3 symmetric matrix (analytical) */
static void eigen3x3(double A[3][3], double evals[3]) {
    double p1 = A[0][1]*A[0][1] + A[0][2]*A[0][2] + A[1][2]*A[1][2];
    double q = (A[0][0] + A[1][1] + A[2][2]) / 3.0;
    double p2 = (A[0][0]-q)*(A[0][0]-q) + (A[1][1]-q)*(A[1][1]-q) + (A[2][2]-q)*(A[2][2]-q) + 2*p1;
    double p = sqrt(p2 / 6.0);

    if (p < 1e-30) { evals[0]=evals[1]=evals[2]=q; return; }

    double B[3][3];
    for (int i=0;i<3;i++) for (int j=0;j<3;j++)
        B[i][j] = (A[i][j] - (i==j?q:0)) / p;
    double detB = B[0][0]*(B[1][1]*B[2][2]-B[1][2]*B[2][1])
                - B[0][1]*(B[1][0]*B[2][2]-B[1][2]*B[2][0])
                + B[0][2]*(B[1][0]*B[2][1]-B[1][1]*B[2][0]);
    double r = detB / 2.0;
    if (r < -1) r = -1; if (r > 1) r = 1;
    double phi = acos(r) / 3.0;

    evals[0] = q + 2*p*cos(phi);
    evals[2] = q + 2*p*cos(phi + 2*PI/3);
    evals[1] = 3*q - evals[0] - evals[2];

    /* Sort descending */
    for (int i=0;i<2;i++) for (int j=i+1;j<3;j++)
        if (evals[j] > evals[i]) { double t=evals[i]; evals[i]=evals[j]; evals[j]=t; }
}

/* Spherical shell energy profile */
static void shell_energy(Grid *g, int n_shells, double *E_shells, double R_max) {
    const int N = g->N, NN = N*N;
    const long N3 = g->N3;
    const double dx = g->dx, L = g->L, dV = dx*dx*dx;
    double dr = R_max / n_shells;

    memset(E_shells, 0, n_shells * sizeof(double));

    #pragma omp parallel
    {
        double *local = calloc(n_shells, sizeof(double));
        #pragma omp for schedule(static)
        for (long idx = 0; idx < N3; idx++) {
            int i=(int)(idx/NN), j=(int)((idx/N)%N), k=(int)(idx%N);
            double x = -L + i*dx, y = -L + j*dx, z = -L + k*dx;
            double r = sqrt(x*x + y*y + z*z);
            int shell = (int)(r / dr);
            if (shell >= n_shells) continue;

            double e_dens = 0;
            for (int a = 0; a < NFIELDS; a++)
                e_dens += 0.5*g->phi[a][idx]*g->phi[a][idx];
            double P = g->phi[0][idx]*g->phi[1][idx]*g->phi[2][idx];
            e_dens += (MU/2.0)*P*P/(1.0+KAPPA*P*P);
            local[shell] += e_dens * dV;
        }
        #pragma omp critical
        for (int s=0;s<n_shells;s++) E_shells[s] += local[s];
        free(local);
    }
}

/* ================================================================
   Fragmentation detection
   ================================================================ */

static int detect_clusters(Grid *g, double threshold, int stride) {
    const int N = g->N, NN = N*N;
    int Nc = (N + stride - 1) / stride;
    long Nc3 = (long)Nc * Nc * Nc;

    int *mask = calloc(Nc3, sizeof(int));
    for (int ci = 0; ci < Nc; ci++) {
        int i = ci * stride; if (i >= N) i = N-1;
        for (int cj = 0; cj < Nc; cj++) {
            int j = cj * stride; if (j >= N) j = N-1;
            for (int ck = 0; ck < Nc; ck++) {
                int k = ck * stride; if (k >= N) k = N-1;
                long idx = (long)i*NN + j*N + k;
                double P = fabs(g->phi[0][idx] * g->phi[1][idx] * g->phi[2][idx]);
                long cidx = (long)ci * Nc * Nc + cj * Nc + ck;
                mask[cidx] = (P > threshold) ? 1 : 0;
            }
        }
    }

    int n_clusters = 0;
    long *stack = malloc(Nc3 * sizeof(long));

    for (long cidx = 0; cidx < Nc3; cidx++) {
        if (mask[cidx] != 1) continue;
        n_clusters++;
        int sp = 0;
        stack[sp++] = cidx;
        mask[cidx] = -1;

        while (sp > 0) {
            long cur = stack[--sp];
            int ci = (int)(cur / (Nc*Nc));
            int cj = (int)((cur / Nc) % Nc);
            int ck = (int)(cur % Nc);

            int di[6] = {-1,1,0,0,0,0};
            int dj[6] = {0,0,-1,1,0,0};
            int dk[6] = {0,0,0,0,-1,1};
            for (int n = 0; n < 6; n++) {
                int ni = ci+di[n], nj = cj+dj[n], nk = ck+dk[n];
                if (ni < 0 || ni >= Nc || nj < 0 || nj >= Nc || nk < 0 || nk >= Nc) continue;
                long nidx = (long)ni*Nc*Nc + nj*Nc + nk;
                if (mask[nidx] == 1) {
                    mask[nidx] = -1;
                    stack[sp++] = nidx;
                }
            }
        }
    }

    free(mask);
    free(stack);
    return n_clusters;
}

/* ================================================================
   Time-averaged death check
   ================================================================ */

static void ep_record(double ep) {
    ep_history[ep_idx] = ep;
    ep_idx = (ep_idx + 1) % DEATH_WINDOW;
    if (ep_count < DEATH_WINDOW) ep_count++;
}

static double ep_average(void) {
    if (ep_count == 0) return 0;
    double sum = 0;
    for (int i = 0; i < ep_count; i++) sum += fabs(ep_history[i]);
    return sum / ep_count;
}

/* ================================================================
   Main
   ================================================================ */

int main(int argc, char **argv) {
    int N = 128;
    double L = 25.0, T = 300.0;
    double diag_dt = 1.0, snap_dt = 5.0;
    char outdir[256] = "data/knot";

    for (int i = 1; i < argc; i++) {
        if      (!strcmp(argv[i],"-N"))      N = atoi(argv[++i]);
        else if (!strcmp(argv[i],"-L"))      L = atof(argv[++i]);
        else if (!strcmp(argv[i],"-T"))      T = atof(argv[++i]);
        else if (!strcmp(argv[i],"-bg"))     A_BG = atof(argv[++i]);
        else if (!strcmp(argv[i],"-m"))      { double m = atof(argv[++i]); MASS2 = m*m; }
        else if (!strcmp(argv[i],"-mt"))     { double m = atof(argv[++i]); MTHETA2 = m*m; }
        else if (!strcmp(argv[i],"-eta"))    ETA = atof(argv[++i]);
        else if (!strcmp(argv[i],"-mu"))     MU = atof(argv[++i]);
        else if (!strcmp(argv[i],"-kappa"))  KAPPA = atof(argv[++i]);
        else if (!strcmp(argv[i],"-diag"))   diag_dt = atof(argv[++i]);
        else if (!strcmp(argv[i],"-snap"))   snap_dt = atof(argv[++i]);
        else if (!strcmp(argv[i],"-o"))      strncpy(outdir, argv[++i], 255);
        else if (!strcmp(argv[i],"-R_major"))  R_MAJOR = atof(argv[++i]);
        else if (!strcmp(argv[i],"-r_minor"))  R_MINOR = atof(argv[++i]);
        else if (!strcmp(argv[i],"-f8_scale")) F8_SCALE = atof(argv[++i]);
        else if (!strcmp(argv[i],"-tube_r"))   TUBE_R = atof(argv[++i]);
        else if (!strcmp(argv[i],"-A"))         A_AMP = atof(argv[++i]);
        else if (!strcmp(argv[i],"-n_osc"))     N_OSC = atoi(argv[++i]);
        else if (!strcmp(argv[i],"-damp_width")) DAMP_WIDTH = atof(argv[++i]);
        else if (!strcmp(argv[i],"-damp_rate"))   DAMP_RATE = atof(argv[++i]);
        else if (!strcmp(argv[i],"-geom")) {
            i++;
            if      (!strcmp(argv[i], "trefoil_R")) { GEOM_TYPE = GEOM_TREFOIL_R; strncpy(GEOM_NAME, "trefoil_R", 31); }
            else if (!strcmp(argv[i], "trefoil_L")) { GEOM_TYPE = GEOM_TREFOIL_L; strncpy(GEOM_NAME, "trefoil_L", 31); }
            else if (!strcmp(argv[i], "figure8"))   { GEOM_TYPE = GEOM_FIGURE8;   strncpy(GEOM_NAME, "figure8", 31); }
            else {
                fprintf(stderr, "Unknown geometry: %s (use trefoil_R, trefoil_L, figure8)\n", argv[i]);
                return 1;
            }
        }
        else {
            fprintf(stderr, "Unknown option: %s\n", argv[i]);
            return 1;
        }
    }

    if (GEOM_TYPE < 0) {
        fprintf(stderr, "Error: -geom required (trefoil_R, trefoil_L, figure8)\n");
        return 1;
    }

    /* Auto-select N_OSC if not set */
    if (N_OSC < 0) {
        N_OSC = (GEOM_TYPE == GEOM_FIGURE8) ? 4 : 3;
    }

    int nthreads = 4;
    char *env_threads = getenv("OMP_NUM_THREADS");
    if (env_threads) nthreads = atoi(env_threads);
    omp_set_num_threads(nthreads);

    printf("=== V37 Knot: geom=%s ===\n", GEOM_NAME);
    printf("Equations: 6-field Cosserat (3 phi + 3 theta)\n");
    printf("m^2=%.4f, m_theta^2=%.4f, eta=%.3f, mu=%.3f, kappa=%.1f, A_bg=%.2f\n",
           MASS2, MTHETA2, ETA, MU, KAPPA, A_BG);
    printf("N=%d L=%.0f T=%.0f\n", N, L, T);

    /* Generate knot curve */
    KnotCurve kc;
    switch (GEOM_TYPE) {
        case GEOM_TREFOIL_R:
            printf("Trefoil (right-handed): R_major=%.1f r_minor=%.1f\n", R_MAJOR, R_MINOR);
            knot_curve_trefoil(&kc, R_MAJOR, R_MINOR, 0);
            break;
        case GEOM_TREFOIL_L:
            printf("Trefoil (left-handed): R_major=%.1f r_minor=%.1f\n", R_MAJOR, R_MINOR);
            knot_curve_trefoil(&kc, R_MAJOR, R_MINOR, 1);
            break;
        case GEOM_FIGURE8:
            printf("Figure-eight: scale=%.1f\n", F8_SCALE);
            knot_curve_figure8(&kc, F8_SCALE);
            break;
    }
    knot_compute_arclength(&kc);

    int crossings = knot_count_crossings(&kc);
    printf("Knot curve: L_curve=%.3f, crossings=%d\n", kc.L_curve, crossings);
    printf("  bbox: x=[%.2f, %.2f] y=[%.2f, %.2f] z=[%.2f, %.2f]\n",
           kc.bbox_min[0], kc.bbox_max[0],
           kc.bbox_min[1], kc.bbox_max[1],
           kc.bbox_min[2], kc.bbox_max[2]);

    mkdir("data", 0755); mkdir(outdir, 0755);
    Grid *g = grid_alloc(N, L);
    printf("dx=%.4f dt=%.5f threads=%d\n\n", g->dx, g->dt, nthreads);

    /* Initialize knot tube */
    init_knot(g, &kc);

    compute_forces(g);

    /* SFA archive with f32 columns */
    char sfapath[512];
    snprintf(sfapath, sizeof(sfapath), "%s.sfa", outdir);
    SFA *sfa = sfa_create(sfapath, N, N, N, L, L, L, g->dt);
    sfa_add_column(sfa, "phi_x",   SFA_F32, SFA_POSITION, 0);
    sfa_add_column(sfa, "phi_y",   SFA_F32, SFA_POSITION, 1);
    sfa_add_column(sfa, "phi_z",   SFA_F32, SFA_POSITION, 2);
    sfa_add_column(sfa, "theta_x", SFA_F32, SFA_ANGLE,    0);
    sfa_add_column(sfa, "theta_y", SFA_F32, SFA_ANGLE,    1);
    sfa_add_column(sfa, "theta_z", SFA_F32, SFA_ANGLE,    2);
    sfa_finalize_header(sfa);
    printf("SFA: %s (f32 columns)\n\n", sfapath);

    /* Allocate f32 conversion buffers (once) */
    float *f32_buf[6];
    for (int c = 0; c < 6; c++) {
        f32_buf[c] = malloc(g->N3 * sizeof(float));
        if (!f32_buf[c]) { fprintf(stderr, "FATAL: f32_buf malloc failed\n"); exit(1); }
    }

    /* Helper: pointers to source double arrays */
    double *src_double[6];
    src_double[0] = g->phi[0]; src_double[1] = g->phi[1]; src_double[2] = g->phi[2];
    src_double[3] = g->theta[0]; src_double[4] = g->theta[1]; src_double[5] = g->theta[2];

    /* Write t=0 frame */
    {
        for (int c = 0; c < 6; c++)
            for (long ii = 0; ii < g->N3; ii++)
                f32_buf[c][ii] = (float)src_double[c][ii];
        void *cols[] = {f32_buf[0], f32_buf[1], f32_buf[2], f32_buf[3], f32_buf[4], f32_buf[5]};
        sfa_write_frame(sfa, 0, cols);
    }

    /* Initial diagnostics */
    {
        double Ixx[3][3], evals[3], E_int;
        compute_inertia_tensor(g, Ixx, &E_int);
        eigen3x3(Ixx, evals);
        double aspect = (evals[2] > 1e-30) ? evals[0]/evals[2] : 999;
        printf("Initial inertia eigenvalues: %.1f  %.1f  %.1f  (aspect=%.2f)\n",
               evals[0], evals[1], evals[2], aspect);
        double P_int = P_integrated(g);
        printf("Initial |P| integrated: %.4e\n\n", P_int);
    }

    /* Timeseries file */
    char tspath[512];
    snprintf(tspath, sizeof(tspath), "%s/timeseries.tsv", outdir);
    FILE *fp = fopen(tspath, "w");
    fprintf(fp, "t\tE_phi_kin\tE_theta_kin\tE_grad\tE_mass\tE_pot\tE_tgrad\tE_tmass\tE_coupling\tE_total\ttheta_rms\tP_int\taspect\tn_clusters\n");

    int n_steps = (int)(T / g->dt);
    int diag_every = (int)(diag_dt / g->dt); if (diag_every<1) diag_every=1;
    int snap_every = (int)(snap_dt / g->dt); if (snap_every<1) snap_every=1;

    printf("Steps=%d diag_every=%d snap_every=%d (snap_dt=%.1f)\n\n", n_steps, diag_every, snap_every, snap_dt);

    double wall0 = omp_get_wtime();
    double E0 = 0, P0 = 0, Ep0 = 0;

    /* Time-averaged death check parameters */
    double death_time_window = 50.0;
    int death_min_samples = (int)(death_time_window / diag_dt);
    if (death_min_samples < 20) death_min_samples = 20;
    if (death_min_samples > DEATH_WINDOW) death_min_samples = DEATH_WINDOW;

    int fragmented = 0;

    for (int step = 0; step <= n_steps; step++) {
        if (step > 0) verlet_step(g);
        double t = step * g->dt;

        /* SFA snapshot */
        if (step > 0 && step % snap_every == 0) {
            for (int c = 0; c < 6; c++)
                for (long ii = 0; ii < g->N3; ii++)
                    f32_buf[c][ii] = (float)src_double[c][ii];
            void *cols[] = {f32_buf[0], f32_buf[1], f32_buf[2], f32_buf[3], f32_buf[4], f32_buf[5]};
            sfa_write_frame(sfa, t, cols);
        }

        if (step % diag_every == 0) {
            double epk, etk, eg, em, ep, etg, etm, ec, et;
            compute_energy(g, &epk, &etk, &eg, &em, &ep, &etg, &etm, &ec, &et);
            double trms = theta_rms(g);
            double P_int = P_integrated(g);

            /* Inertia tensor + fragmentation at major intervals */
            int is_major = (step % (diag_every*50) == 0);
            static double last_aspect = 1.0;
            static int last_clusters = 1;
            double aspect = last_aspect;
            int n_clusters = last_clusters;

            if (is_major) {
                double Ixx[3][3], evals[3], E_int;
                compute_inertia_tensor(g, Ixx, &E_int);
                eigen3x3(Ixx, evals);
                aspect = (evals[2] > 1e-30) ? evals[0]/evals[2] : 999;
                last_aspect = aspect;

                /* Cluster detection */
                double P_thresh = 0.001;
                n_clusters = detect_clusters(g, P_thresh, 4);
                last_clusters = n_clusters;

                if (n_clusters > 3 && !fragmented) {
                    fragmented = 1;
                    printf("  *** FRAGMENTATION DETECTED at t=%.1f (clusters=%d) ***\n",
                           t, n_clusters);
                }
            }

            if (step == 0) { E0 = et; P0 = P_int; Ep0 = ep; }

            /* Record for time-averaged death check */
            ep_record(ep);

            fprintf(fp, "%.2f\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%.4f\t%d\n",
                    t, epk, etk, eg, em, ep, etg, etm, ec, et, trms, P_int, aspect, n_clusters);
            fflush(fp);

            if (is_major) {
                double wall = omp_get_wtime() - wall0;
                double drift = 100.0*(et-E0)/(fabs(E0)+1e-30);
                double ep_avg = ep_average();
                printf("t=%7.1f E=%.3e (drift %+.3f%%) Ep=%.1f Ep_avg=%.1f P_int=%.3e trms=%.3e asp=%.2f clust=%d [%.0f%% %.0fs]\n",
                       t, et, drift, ep, ep_avg, P_int, trms, aspect, n_clusters, 100.0*step/n_steps, wall);
                fflush(stdout);
            }

            /* TIME-AVERAGED death check */
            if (ep_count >= death_min_samples && fabs(Ep0) > 1e-10) {
                double ep_avg = ep_average();
                if (ep_avg < 0.01 * fabs(Ep0)) {
                    printf("\n  *** Time-averaged |E_pot| (%.3e) < 1%% of initial (%.3e) -- knot DEAD. Stopping. ***\n",
                           ep_avg, fabs(Ep0));
                    for (int c = 0; c < 6; c++)
                        for (long ii = 0; ii < g->N3; ii++)
                            f32_buf[c][ii] = (float)src_double[c][ii];
                    void *cols[] = {f32_buf[0], f32_buf[1], f32_buf[2], f32_buf[3], f32_buf[4], f32_buf[5]};
                    sfa_write_frame(sfa, t, cols);
                    break;
                }
            }
        }
    }

    /* Final frame */
    {
        double t_final = n_steps * g->dt;
        for (int c = 0; c < 6; c++)
            for (long ii = 0; ii < g->N3; ii++)
                f32_buf[c][ii] = (float)src_double[c][ii];
        void *cols[] = {f32_buf[0], f32_buf[1], f32_buf[2], f32_buf[3], f32_buf[4], f32_buf[5]};
        sfa_write_frame(sfa, t_final, cols);
    }
    uint32_t sfa_nframes = sfa->total_frames;
    sfa_close(sfa);
    printf("\nSFA closed: %s (%u frames)\n", sfapath, sfa_nframes);

    /* Final shell analysis */
    printf("\n--- Spherical Shell Energy Profile (final state) ---\n");
    int n_shells = 20;
    double R_max = L;
    double E_shells[20];
    shell_energy(g, n_shells, E_shells, R_max);
    double dr = R_max / n_shells;
    for (int s = 0; s < n_shells; s++)
        printf("  r=[%.1f, %.1f]: E_shell=%.4e\n", s*dr, (s+1)*dr, E_shells[s]);

    /* Final summary */
    {
        double epk, etk, eg, em, ep, etg, etm, ec, et;
        compute_energy(g, &epk, &etk, &eg, &em, &ep, &etg, &etm, &ec, &et);
        double trms = theta_rms(g);
        double P_int = P_integrated(g);
        double Ixx[3][3], evals[3], E_int;
        compute_inertia_tensor(g, Ixx, &E_int);
        eigen3x3(Ixx, evals);
        double aspect = (evals[2] > 1e-30) ? evals[0]/evals[2] : 999;
        int n_clusters = detect_clusters(g, 0.001, 4);

        printf("\n=== FINAL SUMMARY: knot geom=%s ===\n", GEOM_NAME);
        printf("E_total=%.4e (drift %.3f%%)\n", et, 100*(et-E0)/(fabs(E0)+1e-30));
        printf("E_pot=%.4e (%.1f%% of initial)\n", ep, 100*ep/(Ep0+1e-30));
        printf("|P| integrated=%.4e (%.1f%% of initial)\n", P_int, 100*P_int/(P0+1e-30));
        printf("theta_rms=%.4e\n", trms);
        printf("Aspect ratio=%.3f (eigenvalues: %.1f %.1f %.1f)\n",
               aspect, evals[0], evals[1], evals[2]);
        printf("Clusters=%d  Fragmented=%s\n",
               n_clusters, fragmented ? "YES" : "NO");
        printf("Survival: %s\n", (fabs(ep) > 0.1*fabs(Ep0)) ? "YES (instantaneous)" : "NO (instantaneous)");
        double ep_avg = ep_average();
        printf("Survival (avg): %s (Ep_avg=%.3e vs Ep0=%.3e)\n",
               (ep_avg > 0.1*fabs(Ep0)) ? "YES" : "NO", ep_avg, Ep0);
    }

    fclose(fp);

    /* Free f32 buffers */
    for (int c = 0; c < 6; c++) free(f32_buf[c]);

    double wall = omp_get_wtime() - wall0;
    printf("\n=== Complete: %.0fs (%.1f min) ===\n", wall, wall/60);

    grid_free(g);
    return 0;
}
