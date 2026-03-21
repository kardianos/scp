/*  v37_compact.c — Compact braid geometry search (6-field Cosserat)
 *
 *  Tests four compact braid topologies to replace the infinite helical tube:
 *    1. truncated  — Gaussian z-envelope on standard helix
 *    2. borromean   — Three interlocking rings in orthogonal planes
 *    3. loop        — Double helix bent into a torus (closed loop)
 *    4. triloop     — Three helical torus loops, one per plane
 *
 *  Equations (same as v33):
 *    d²phi_a/dt² = Lap(phi_a) - m² phi_a - V'(P) + eta * curl(theta)_a
 *    d²theta_a/dt² = Lap(theta_a) - m_theta² theta_a + eta * curl(phi)_a
 *
 *  V(P) = (mu/2) P² / (1 + kappa P²),  P = phi_0 phi_1 phi_2
 *
 *  Build: gcc -O3 -march=native -fopenmp -o v37_sfa src/v37_sfa.c -lzstd -lm
 *  Run:   ./v37_sfa -geom truncated -o data/truncated
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

/* Physics parameters */
static double MU      = -41.345;
static double KAPPA   = 50.0;
static double MASS2   = 2.25;     /* position field mass² */
static double MTHETA2 = 0.0;     /* angle field mass² (default massless) */
static double A_BG    = 0.1;
static double ETA     = 0.5;     /* coupling strength */

/* Geometry types */
enum GeomType { GEOM_TRUNCATED=0, GEOM_BORROMEAN, GEOM_LOOP, GEOM_TRILOOP, GEOM_BRAID3 };
static const char *geom_names[] = {"truncated", "borromean", "loop", "triloop", "braid3"};

/* Braid3-specific parameters */
static double R_HELIX = 1.0;     /* braid3: helix radius */
static double R_TUBE  = 2.0;     /* braid3: tube radius */
static double SIGMA_Z = 3.0;     /* braid3: z-envelope width */

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
   Verlet integrator
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
}

/* ================================================================
   Initialization: Geometry 1 — Truncated Helix
   ================================================================ */

static void init_truncated(Grid *g) {
    const int N = g->N, NN = N*N;
    const double dx = g->dx, L = g->L;
    const double A[3] = {0.8, 0.8, 0.8};
    const double delta[3] = {0, 3.0005, 4.4325};
    const double R_tube = 3.0, ellip = 0.3325;
    const double sigma_z = R_tube;  /* z-envelope width = tube radius */
    const double kw = PI/L, omega = sqrt(kw*kw + MASS2);
    const double sx = 1+ellip, sy = 1-ellip;
    const double inv2R2 = 1.0/(2*R_tube*R_tube);
    const double k_bg = PI/L, omega_bg = sqrt(k_bg*k_bg + MASS2);

    printf("  Init: truncated helix (sigma_z=%.1f)\n", sigma_z);

    for (int i = 0; i < N; i++) {
        double x = -L + i*dx;
        for (int j = 0; j < N; j++) {
            double y = -L + j*dx;
            for (int kk = 0; kk < N; kk++) {
                double z = -L + kk*dx;
                long idx = (long)i*NN + j*N + kk;
                double r2e = x*x/(sx*sx) + y*y/(sy*sy);
                double z_env = exp(-z*z / (2*sigma_z*sigma_z));
                double env = exp(-r2e * inv2R2) * z_env;
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

/* ================================================================
   Initialization: Geometry 2 — Borromean Rings
   ================================================================ */

static void init_borromean(Grid *g) {
    const int N = g->N, NN = N*N;
    const double dx = g->dx, L = g->L;
    const double R_ring = 4.0;
    const double tube_r = 1.5;
    const double amp = 0.8;
    const double delta[3] = {0, 2*PI/3, 4*PI/3};
    const double inv2t2 = 1.0/(2*tube_r*tube_r);

    printf("  Init: Borromean rings (R=%.1f, tube=%.1f)\n", R_ring, tube_r);

    for (int i = 0; i < N; i++) {
        double x = -L + i*dx;
        for (int j = 0; j < N; j++) {
            double y = -L + j*dx;
            for (int kk = 0; kk < N; kk++) {
                double z = -L + kk*dx;
                long idx = (long)i*NN + j*N + kk;

                /* Ring 0: xy-plane (axis = z) */
                double r_xy = sqrt(x*x + y*y);
                double d0 = sqrt((r_xy - R_ring)*(r_xy - R_ring) + z*z);
                double angle0 = atan2(y, x);
                double env0 = amp * exp(-d0*d0 * inv2t2);

                /* Ring 1: xz-plane (axis = y) */
                double r_xz = sqrt(x*x + z*z);
                double d1 = sqrt((r_xz - R_ring)*(r_xz - R_ring) + y*y);
                double angle1 = atan2(z, x);
                double env1 = amp * exp(-d1*d1 * inv2t2);

                /* Ring 2: yz-plane (axis = x) */
                double r_yz = sqrt(y*y + z*z);
                double d2 = sqrt((r_yz - R_ring)*(r_yz - R_ring) + x*x);
                double angle2 = atan2(z, y);
                double env2 = amp * exp(-d2*d2 * inv2t2);

                /* Each field is primarily one ring, with background */
                double k_bg = PI/L;
                g->phi[0][idx] = env0 * cos(angle0 + delta[0])
                               + 0.3*env1 * cos(angle1)
                               + 0.3*env2 * cos(angle2)
                               + A_BG * cos(k_bg*x);
                g->phi[1][idx] = env1 * cos(angle1 + delta[1])
                               + 0.3*env0 * cos(angle0)
                               + 0.3*env2 * cos(angle2)
                               + A_BG * cos(k_bg*y);
                g->phi[2][idx] = env2 * cos(angle2 + delta[2])
                               + 0.3*env1 * cos(angle1)
                               + 0.3*env0 * cos(angle0)
                               + A_BG * cos(k_bg*z);

                /* Velocities: standing-wave style (omega * A * sin) */
                double omega = sqrt(MASS2);  /* k~0 for these localized modes */
                g->phi_vel[0][idx] = omega * g->phi[0][idx] * 0.3;
                g->phi_vel[1][idx] = omega * g->phi[1][idx] * 0.3;
                g->phi_vel[2][idx] = omega * g->phi[2][idx] * 0.3;
            }
        }
    }
}

/* ================================================================
   Initialization: Geometry 3 — Double Helix Closed Loop (Torus)
   ================================================================ */

static void init_loop(Grid *g) {
    const int N = g->N, NN = N*N;
    const double dx = g->dx, L = g->L;
    const double R_major = 5.0;
    const double R_minor = 1.5;
    const double k_twist = 3.0;   /* three full twists around the loop */
    const double amp = 0.8;
    const double delta[3] = {0, 2*PI/3, 4*PI/3};
    const double inv2R2 = 1.0/(2*R_minor*R_minor);

    printf("  Init: closed loop (R_major=%.1f, R_minor=%.1f, k_twist=%.0f)\n",
           R_major, R_minor, k_twist);

    for (int i = 0; i < N; i++) {
        double x = -L + i*dx;
        for (int j = 0; j < N; j++) {
            double y = -L + j*dx;
            for (int kk = 0; kk < N; kk++) {
                double z = -L + kk*dx;
                long idx = (long)i*NN + j*N + kk;

                double r_cyl = sqrt(x*x + y*y);
                double theta_tor = atan2(y, x);
                double r_tube = sqrt((r_cyl - R_major)*(r_cyl - R_major) + z*z);
                double env = amp * exp(-r_tube*r_tube * inv2R2);

                for (int a = 0; a < NFIELDS; a++) {
                    double ph = k_twist * theta_tor + delta[a];
                    double k_bg = PI/L;
                    double ph_bg = k_bg * (x*cos(2*PI*a/3) + y*sin(2*PI*a/3));
                    g->phi[a][idx] = env * cos(ph) + A_BG * cos(ph_bg);
                    /* Velocity: oscillating mode */
                    double omega = sqrt(MASS2);
                    g->phi_vel[a][idx] = omega * env * sin(ph) * 0.5;
                }
            }
        }
    }
}

/* ================================================================
   Initialization: Geometry 4 — Triple Closed-Loop Helix
   ================================================================ */

static void init_triloop(Grid *g) {
    const int N = g->N, NN = N*N;
    const double dx = g->dx, L = g->L;
    const double R_major = 4.0;
    const double R_minor = 1.2;
    const double k_twist = 2.0;   /* two twists per loop */
    const double amp = 0.8;
    const double delta[3] = {0, 2*PI/3, 4*PI/3};
    const double inv2R2 = 1.0/(2*R_minor*R_minor);

    printf("  Init: triple closed-loop (R_major=%.1f, R_minor=%.1f, k_twist=%.0f)\n",
           R_major, R_minor, k_twist);

    for (int i = 0; i < N; i++) {
        double x = -L + i*dx;
        for (int j = 0; j < N; j++) {
            double y = -L + j*dx;
            for (int kk = 0; kk < N; kk++) {
                double z = -L + kk*dx;
                long idx = (long)i*NN + j*N + kk;

                /* Loop 0: torus in xy-plane (axis = z) */
                double r0_cyl = sqrt(x*x + y*y);
                double theta0 = atan2(y, x);
                double r0_tube = sqrt((r0_cyl - R_major)*(r0_cyl - R_major) + z*z);
                double env0 = amp * exp(-r0_tube*r0_tube * inv2R2);

                /* Loop 1: torus in xz-plane (axis = y) */
                double r1_cyl = sqrt(x*x + z*z);
                double theta1 = atan2(z, x);
                double r1_tube = sqrt((r1_cyl - R_major)*(r1_cyl - R_major) + y*y);
                double env1 = amp * exp(-r1_tube*r1_tube * inv2R2);

                /* Loop 2: torus in yz-plane (axis = x) */
                double r2_cyl = sqrt(y*y + z*z);
                double theta2 = atan2(z, y);
                double r2_tube = sqrt((r2_cyl - R_major)*(r2_cyl - R_major) + x*x);
                double env2 = amp * exp(-r2_tube*r2_tube * inv2R2);

                for (int a = 0; a < NFIELDS; a++) {
                    /* Each field associated primarily with one loop */
                    double val = 0, vel = 0;
                    double omega = sqrt(MASS2);

                    /* Primary loop contribution */
                    double ph0 = k_twist * theta0 + delta[a];
                    double ph1 = k_twist * theta1 + delta[a];
                    double ph2 = k_twist * theta2 + delta[a];

                    if (a == 0) {
                        val = env0 * cos(ph0) + 0.3*env1*cos(ph1) + 0.3*env2*cos(ph2);
                        vel = omega*(env0*sin(ph0) + 0.3*env1*sin(ph1) + 0.3*env2*sin(ph2))*0.5;
                    } else if (a == 1) {
                        val = env1 * cos(ph1) + 0.3*env0*cos(ph0) + 0.3*env2*cos(ph2);
                        vel = omega*(env1*sin(ph1) + 0.3*env0*sin(ph0) + 0.3*env2*sin(ph2))*0.5;
                    } else {
                        val = env2 * cos(ph2) + 0.3*env0*cos(ph0) + 0.3*env1*cos(ph1);
                        vel = omega*(env2*sin(ph2) + 0.3*env0*sin(ph0) + 0.3*env1*sin(ph1))*0.5;
                    }

                    /* Background */
                    double k_bg = PI/L;
                    double dir = (a==0) ? x : (a==1) ? y : z;
                    val += A_BG * cos(k_bg * dir);

                    g->phi[a][idx] = val;
                    g->phi_vel[a][idx] = vel;
                }
            }
        }
    }
}

/* ================================================================
   Initialization: Geometry 5 — Triple-strand truncated braid (braid3)
   ================================================================ */

static void init_braid3(Grid *g) {
    const int N = g->N, NN = N*N;
    const double dx = g->dx, L = g->L;
    const double A_amp = 0.8;
    const double delta[3] = {0, 3.0005, 4.4325};
    const double Rh = R_HELIX;
    const double r_tube = R_TUBE;
    const double sigma_z = SIGMA_Z;
    const double inv2t2 = 1.0 / (2.0 * r_tube * r_tube);
    const double inv2sz2 = 1.0 / (2.0 * sigma_z * sigma_z);
    const double kz = PI / L;
    const double omega = sqrt(kz * kz + MASS2);
    const double k_bg = PI / L;
    const double omega_bg = sqrt(k_bg * k_bg + MASS2);

    printf("  Init: braid3 (R_helix=%.1f, r_tube=%.1f, sigma_z=%.1f)\n",
           Rh, r_tube, sigma_z);

    for (int i = 0; i < N; i++) {
        double x = -L + i * dx;
        for (int j = 0; j < N; j++) {
            double y = -L + j * dx;
            for (int kk = 0; kk < N; kk++) {
                double z = -L + kk * dx;
                long idx = (long)i * NN + j * N + kk;

                /* z-envelope for truncation */
                double z_env = exp(-z * z * inv2sz2);

                for (int a = 0; a < NFIELDS; a++) {
                    double val = 0, vel = 0;

                    /* Accumulate from all three strands */
                    for (int b = 0; b < 3; b++) {
                        double cx_b = Rh * cos(kz * z + 2.0 * PI * b / 3.0);
                        double cy_b = Rh * sin(kz * z + 2.0 * PI * b / 3.0);
                        double dx_b = x - cx_b;
                        double dy_b = y - cy_b;
                        double d2_b = dx_b * dx_b + dy_b * dy_b;
                        double env_b = exp(-d2_b * inv2t2);

                        double phase = kz * z + delta[a];
                        val += A_amp * env_b * cos(phase);
                        vel += omega * A_amp * env_b * sin(phase);
                    }

                    /* Apply z-envelope */
                    val *= z_env;
                    vel *= z_env;

                    /* Add background */
                    double ph_bg = k_bg * z + 2.0 * PI * a / 3.0;
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
   Diagnostics
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

/* Triple product |P| integrated (measures braid survival) */
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

/* Moment of inertia tensor from field amplitude (fast: no gradient computation) */
static void compute_inertia_tensor(Grid *g, double Ixx[3][3], double *E_tot_out) {
    const int N = g->N, NN = N*N;
    const long N3 = g->N3;
    const double dx = g->dx, L = g->L;
    const double dV = dx*dx*dx;

    /* Use phi^2 as weight (fast, no gradients needed) */
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

/* Spherical shell energy profile (for depletion check) */
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

/* Save binary snapshot */
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
    printf("  [SNAP] %s (%.1f GB, 6 fields)\n", fn, 6.0*g->N3*8.0/1e9);
}

/* ================================================================
   Main
   ================================================================ */

int main(int argc, char **argv) {
    int N = 128;
    double L = 15.0, T = 500.0;
    double diag_dt = 5.0, snap_dt = 100.0;
    char outdir[256] = "data/compact";
    enum GeomType geom = GEOM_TRUNCATED;

    for (int i = 1; i < argc; i++) {
        if      (!strcmp(argv[i],"-N"))      N = atoi(argv[++i]);
        else if (!strcmp(argv[i],"-L"))      L = atof(argv[++i]);
        else if (!strcmp(argv[i],"-T"))      T = atof(argv[++i]);
        else if (!strcmp(argv[i],"-bg"))     A_BG = atof(argv[++i]);
        else if (!strcmp(argv[i],"-m"))      { double m = atof(argv[++i]); MASS2 = m*m; }
        else if (!strcmp(argv[i],"-mt"))     { double m = atof(argv[++i]); MTHETA2 = m*m; }
        else if (!strcmp(argv[i],"-eta"))    ETA = atof(argv[++i]);
        else if (!strcmp(argv[i],"-diag"))   diag_dt = atof(argv[++i]);
        else if (!strcmp(argv[i],"-snap"))   snap_dt = atof(argv[++i]);
        else if (!strcmp(argv[i],"-o"))      strncpy(outdir, argv[++i], 255);
        else if (!strcmp(argv[i],"-Rh"))     R_HELIX = atof(argv[++i]);
        else if (!strcmp(argv[i],"-rtube"))  R_TUBE = atof(argv[++i]);
        else if (!strcmp(argv[i],"-sigmaz")) SIGMA_Z = atof(argv[++i]);
        else if (!strcmp(argv[i],"-geom")) {
            i++;
            if      (!strcmp(argv[i],"truncated")) geom = GEOM_TRUNCATED;
            else if (!strcmp(argv[i],"borromean"))  geom = GEOM_BORROMEAN;
            else if (!strcmp(argv[i],"loop"))       geom = GEOM_LOOP;
            else if (!strcmp(argv[i],"triloop"))    geom = GEOM_TRILOOP;
            else if (!strcmp(argv[i],"braid3"))     geom = GEOM_BRAID3;
            else { fprintf(stderr, "Unknown geometry: %s\n", argv[i]); return 1; }
        }
    }

    int nthreads = 8;
    char *env_threads = getenv("OMP_NUM_THREADS");
    if (env_threads) nthreads = atoi(env_threads);
    omp_set_num_threads(nthreads);

    printf("=== V37 Compact Braid: %s ===\n", geom_names[geom]);
    printf("Equations: 6-field Cosserat (3 phi + 3 theta)\n");
    printf("m²=%.4f, m_theta²=%.4f, eta=%.3f, mu=%.3f, kappa=%.1f, A_bg=%.2f\n",
           MASS2, MTHETA2, ETA, MU, KAPPA, A_BG);
    printf("N=%d L=%.0f T=%.0f\n", N, L, T);

    mkdir("data", 0755); mkdir(outdir, 0755);
    Grid *g = grid_alloc(N, L);
    printf("dx=%.4f dt=%.5f threads=%d\n\n", g->dx, g->dt, nthreads);

    /* Initialize based on geometry */
    switch (geom) {
        case GEOM_TRUNCATED: init_truncated(g); break;
        case GEOM_BORROMEAN: init_borromean(g); break;
        case GEOM_LOOP:      init_loop(g);      break;
        case GEOM_TRILOOP:   init_triloop(g);   break;
        case GEOM_BRAID3:    init_braid3(g);    break;
    }

    compute_forces(g);

    /* SFA archive */
    char sfapath[512];
    snprintf(sfapath, sizeof(sfapath), "%s.sfa", outdir);
    SFA *sfa = sfa_create(sfapath, N, N, N, L, L, L, g->dt);
    sfa_add_column(sfa, "phi_x",   SFA_F64, SFA_POSITION, 0);
    sfa_add_column(sfa, "phi_y",   SFA_F64, SFA_POSITION, 1);
    sfa_add_column(sfa, "phi_z",   SFA_F64, SFA_POSITION, 2);
    sfa_add_column(sfa, "theta_x", SFA_F64, SFA_ANGLE,    0);
    sfa_add_column(sfa, "theta_y", SFA_F64, SFA_ANGLE,    1);
    sfa_add_column(sfa, "theta_z", SFA_F64, SFA_ANGLE,    2);
    sfa_finalize_header(sfa);
    printf("SFA: %s\n\n", sfapath);

    /* Write t=0 frame */
    {
        void *cols[] = {g->phi[0],g->phi[1],g->phi[2],g->theta[0],g->theta[1],g->theta[2]};
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
    fprintf(fp, "t\tE_phi_kin\tE_theta_kin\tE_grad\tE_mass\tE_pot\tE_tgrad\tE_tmass\tE_coupling\tE_total\ttheta_rms\tP_int\taspect\n");

    int n_steps = (int)(T / g->dt);
    int diag_every = (int)(diag_dt / g->dt); if (diag_every<1) diag_every=1;
    int snap_every = (int)(snap_dt / g->dt); if (snap_every<1) snap_every=1;

    printf("Steps=%d diag_every=%d snap_every=%d\n\n", n_steps, diag_every, snap_every);

    double wall0 = omp_get_wtime();
    double E0 = 0, P0 = 0, Ep0 = 0;

    for (int step = 0; step <= n_steps; step++) {
        if (step > 0) verlet_step(g);
        double t = step * g->dt;

        if (step % diag_every == 0) {
            double epk, etk, eg, em, ep, etg, etm, ec, et;
            compute_energy(g, &epk, &etk, &eg, &em, &ep, &etg, &etm, &ec, &et);
            double trms = theta_rms(g);
            double P_int = P_integrated(g);

            /* Inertia tensor only at major print intervals (expensive!) */
            int is_major = (step % (diag_every*10) == 0);
            static double last_aspect = 1.0;
            double aspect = last_aspect;
            if (is_major) {
                double Ixx[3][3], evals[3], E_int;
                compute_inertia_tensor(g, Ixx, &E_int);
                eigen3x3(Ixx, evals);
                aspect = (evals[2] > 1e-30) ? evals[0]/evals[2] : 999;
                last_aspect = aspect;
            }

            if (step == 0) { E0 = et; P0 = P_int; Ep0 = ep; }

            fprintf(fp, "%.2f\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%.4f\n",
                    t, epk, etk, eg, em, ep, etg, etm, ec, et, trms, P_int, aspect);
            fflush(fp);

            if (is_major) {
                double wall = omp_get_wtime() - wall0;
                double drift = 100.0*(et-E0)/(fabs(E0)+1e-30);
                printf("t=%7.1f E=%.3e (drift %+.3f%%) Ep=%.1f P_int=%.3e theta_rms=%.3e aspect=%.2f [%.0f%% %.0fs]\n",
                       t, et, drift, ep, P_int, trms, aspect, 100.0*step/n_steps, wall);
                fflush(stdout);
            }

            /* Early termination: if E_pot has gone to ~zero, braid is dead */
            if (step > 0 && fabs(Ep0) > 1e-10 && fabs(ep) < 0.01*fabs(Ep0)) {
                printf("\n  *** E_pot collapsed to <1%% of initial — braid DEAD. Stopping. ***\n");
                void *cols[] = {g->phi[0],g->phi[1],g->phi[2],g->theta[0],g->theta[1],g->theta[2]};
                sfa_write_frame(sfa, t, cols);
                break;
            }
        }

        if (step > 0 && step % snap_every == 0) {
            void *cols[] = {g->phi[0],g->phi[1],g->phi[2],g->theta[0],g->theta[1],g->theta[2]};
            sfa_write_frame(sfa, t, cols);
        }
    }

    /* Final frame */
    {
        double t_final = n_steps * g->dt;
        void *cols[] = {g->phi[0],g->phi[1],g->phi[2],g->theta[0],g->theta[1],g->theta[2]};
        sfa_write_frame(sfa, t_final, cols);
    }
    sfa_close(sfa);

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

        printf("\n=== FINAL SUMMARY: %s ===\n", geom_names[geom]);
        printf("E_total=%.4e (drift %.3f%%)\n", et, 100*(et-E0)/(fabs(E0)+1e-30));
        printf("E_pot=%.4e (%.1f%% of initial)\n", ep, 100*ep/(Ep0+1e-30));
        printf("|P| integrated=%.4e (%.1f%% of initial)\n", P_int, 100*P_int/(P0+1e-30));
        printf("theta_rms=%.4e\n", trms);
        printf("Aspect ratio=%.3f (eigenvalues: %.1f %.1f %.1f)\n",
               aspect, evals[0], evals[1], evals[2]);
        printf("Survival: %s\n", (fabs(ep) > 0.1*fabs(Ep0)) ? "YES" : "NO");
    }

    fclose(fp);

    double wall = omp_get_wtime() - wall0;
    printf("\n=== Complete: %.0fs (%.1f min) ===\n", wall, wall/60);

    grid_free(g);
    return 0;
}
