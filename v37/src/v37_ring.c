/*  v37_ring.c — Intertwined compact braid geometries (6-field Cosserat)
 *
 *  Tests three co-located braid topologies:
 *    1. ring1   — Single torus ring with phase-offset oscillations
 *    2. sph     — Spherical harmonic l=1 modes with radial oscillation
 *    3. braid3  — Triple-strand truncated braid (helical sub-tubes)
 *
 *  All three keep fields CO-LOCATED on the same tube/region, unlike the
 *  failed ring geometries that separated fields onto different planes.
 *
 *  Equations (same as v37_compact / v33):
 *    d²phi_a/dt² = Lap(phi_a) - m² phi_a - V'(P) + eta * curl(theta)_a
 *    d²theta_a/dt² = Lap(theta_a) - m_theta² theta_a + eta * curl(phi)_a
 *
 *  V(P) = (mu/2) P² / (1 + kappa P²),  P = phi_0 phi_1 phi_2
 *
 *  FIX: Death check uses TIME-AVERAGED E_pot over last 50 time units,
 *       not instantaneous E_pot (avoids false kills on breathing zeros).
 *
 *  Build: gcc -O3 -march=native -fopenmp -o v37_ring src/v37_ring.c -lm
 *  Run:   ./v37_ring -geom ring1 -Rmaj 4 -nosc 2 -Rmin 2.0 -o data/ring1_R4n2
 */

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
enum GeomType { GEOM_RING1=0, GEOM_SPH, GEOM_BRAID3 };
static const char *geom_names[] = {"ring1", "sph", "braid3"};

/* Geometry-specific parameters */
static double R_MAJOR = 4.0;     /* ring1: major radius */
static double R_MINOR = 2.0;     /* ring1: tube radius */
static double N_OSC   = 2.0;     /* ring1: oscillations around ring */
static double SIGMA_R = 5.0;     /* sph: radial envelope width */
static double K_RAD   = 0.0;     /* sph: radial wavenumber (default PI/8) */
static double R_HELIX = 1.0;     /* braid3: helix radius */
static double R_TUBE  = 2.0;     /* braid3: tube radius */
static double SIGMA_Z = 3.0;     /* braid3: z-envelope width */

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
   Initialization: Geometry A — Ring1 (single torus, co-located)
   ================================================================ */

static void init_ring1(Grid *g) {
    const int N = g->N, NN = N*N;
    const double dx = g->dx, L = g->L;
    const double Rmaj = R_MAJOR, Rmin = R_MINOR;
    const double n_osc = N_OSC;
    const double A_amp = 0.8;
    const double delta[3] = {0, 3.0005, 4.4325};
    const double inv2R2 = 1.0 / (2.0 * Rmin * Rmin);
    const double k_ring = n_osc;  /* wavenumber: n_osc oscillations per 2pi */
    const double k_bg = PI / L;
    const double omega_bg = sqrt(k_bg * k_bg + MASS2);
    /* omega for ring mode: effective k = n_osc/R_major, plus mass */
    const double omega_ring = sqrt(k_ring * k_ring / (Rmaj * Rmaj) + MASS2);

    printf("  Init: ring1 (R_major=%.1f, R_minor=%.1f, n_osc=%.0f)\n",
           Rmaj, Rmin, n_osc);
    printf("  k_eff=%.3f omega_ring=%.3f delta={0, %.4f, %.4f}\n",
           k_ring/Rmaj, omega_ring, delta[1], delta[2]);

    for (int i = 0; i < N; i++) {
        double x = -L + i * dx;
        for (int j = 0; j < N; j++) {
            double y = -L + j * dx;
            for (int kk = 0; kk < N; kk++) {
                double z = -L + kk * dx;
                long idx = (long)i * NN + j * N + kk;

                /* Toroidal coordinates */
                double r_cyl = sqrt(x * x + y * y);
                double theta_tor = atan2(y, x);
                double d_tube = sqrt((r_cyl - Rmaj) * (r_cyl - Rmaj) + z * z);

                /* Gaussian envelope around the torus tube */
                double env = A_amp * exp(-d_tube * d_tube * inv2R2);

                for (int a = 0; a < NFIELDS; a++) {
                    /* Phase: n_osc oscillations around the ring + field-specific offset */
                    double phase = k_ring * theta_tor + delta[a];
                    double ph_bg = k_bg * z + 2.0 * PI * a / 3.0;

                    g->phi[a][idx] = env * cos(phase) + A_BG * cos(ph_bg);

                    /* Velocity: omega * env * sin(phase) for oscillating mode */
                    g->phi_vel[a][idx] = omega_ring * env * sin(phase)
                                       + omega_bg * A_BG * sin(ph_bg);
                }
            }
        }
    }
}

/* ================================================================
   Initialization: Geometry B — Spherical harmonic modes
   ================================================================ */

static void init_sph(Grid *g) {
    const int N = g->N, NN = N*N;
    const double dx = g->dx, L = g->L;
    const double A_amp = 0.8;
    const double sigma = SIGMA_R;
    const double kr = K_RAD;
    const double inv2s2 = 1.0 / (2.0 * sigma * sigma);
    const double omega_r = sqrt(kr * kr + MASS2);
    const double k_bg = PI / L;
    const double omega_bg = sqrt(k_bg * k_bg + MASS2);

    printf("  Init: sph (sigma_r=%.1f, k_r=%.4f)\n", sigma, kr);
    printf("  omega_r=%.3f\n", omega_r);

    for (int i = 0; i < N; i++) {
        double x = -L + i * dx;
        for (int j = 0; j < N; j++) {
            double y = -L + j * dx;
            for (int kk = 0; kk < N; kk++) {
                double z = -L + kk * dx;
                long idx = (long)i * NN + j * N + kk;

                double r = sqrt(x * x + y * y + z * z);
                double r_safe = r + 1e-30;

                /* Spherical angles */
                double cos_theta = z / r_safe;
                double sin_theta = sqrt(x * x + y * y) / r_safe;
                double cos_phi = (sin_theta > 1e-20) ? x / (r_safe * sin_theta) : 1.0;
                double sin_phi = (sin_theta > 1e-20) ? y / (r_safe * sin_theta) : 0.0;

                /* Radial envelope x oscillation */
                double env = A_amp * exp(-r * r * inv2s2);
                double f_r = env * cos(kr * r);
                /* Velocity radial part: d/dt of cos(kr*r - omega*t) at t=0 = omega*sin(kr*r) */
                double f_r_vel = omega_r * env * sin(kr * r);

                /* Three l=1 spherical harmonics */
                /* phi_0 = Y_10 ~ cos(theta) */
                double ph_bg0 = k_bg * z;
                g->phi[0][idx] = f_r * cos_theta + A_BG * cos(ph_bg0);
                g->phi_vel[0][idx] = f_r_vel * cos_theta + omega_bg * A_BG * sin(ph_bg0);

                /* phi_1 = Y_11 ~ sin(theta)*cos(phi) */
                double ph_bg1 = k_bg * z + 2.0 * PI / 3.0;
                g->phi[1][idx] = f_r * sin_theta * cos_phi + A_BG * cos(ph_bg1);
                g->phi_vel[1][idx] = f_r_vel * sin_theta * cos_phi + omega_bg * A_BG * sin(ph_bg1);

                /* phi_2 = Y_1-1 ~ sin(theta)*sin(phi) */
                double ph_bg2 = k_bg * z + 4.0 * PI / 3.0;
                g->phi[2][idx] = f_r * sin_theta * sin_phi + A_BG * cos(ph_bg2);
                g->phi_vel[2][idx] = f_r_vel * sin_theta * sin_phi + omega_bg * A_BG * sin(ph_bg2);
            }
        }
    }
}

/* ================================================================
   Initialization: Geometry C — Triple-strand truncated braid
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
                        /* Center of strand b at this z */
                        double cx_b = Rh * cos(kz * z + 2.0 * PI * b / 3.0);
                        double cy_b = Rh * sin(kz * z + 2.0 * PI * b / 3.0);
                        double dx_b = x - cx_b;
                        double dy_b = y - cy_b;
                        double d2_b = dx_b * dx_b + dy_b * dy_b;
                        double env_b = exp(-d2_b * inv2t2);

                        /* Field a's phase at this z */
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
    double L = 15.0, T = 200.0;
    double diag_dt = 1.0, snap_dt = 50.0;
    char outdir[256] = "data/ring1";
    enum GeomType geom = GEOM_RING1;
    K_RAD = PI / 8.0;  /* default sph k_r */

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
        else if (!strcmp(argv[i],"-Rmaj"))   R_MAJOR = atof(argv[++i]);
        else if (!strcmp(argv[i],"-Rmin"))   R_MINOR = atof(argv[++i]);
        else if (!strcmp(argv[i],"-nosc"))   N_OSC = atof(argv[++i]);
        else if (!strcmp(argv[i],"-sigr"))   SIGMA_R = atof(argv[++i]);
        else if (!strcmp(argv[i],"-kr"))     K_RAD = atof(argv[++i]);
        else if (!strcmp(argv[i],"-Rh"))     R_HELIX = atof(argv[++i]);
        else if (!strcmp(argv[i],"-rtube"))  R_TUBE = atof(argv[++i]);
        else if (!strcmp(argv[i],"-sigz"))   SIGMA_Z = atof(argv[++i]);
        else if (!strcmp(argv[i],"-geom")) {
            i++;
            if      (!strcmp(argv[i],"ring1"))  geom = GEOM_RING1;
            else if (!strcmp(argv[i],"sph"))     geom = GEOM_SPH;
            else if (!strcmp(argv[i],"braid3"))  geom = GEOM_BRAID3;
            else { fprintf(stderr, "Unknown geometry: %s\n", argv[i]); return 1; }
        }
    }

    int nthreads = 4;
    char *env_threads = getenv("OMP_NUM_THREADS");
    if (env_threads) nthreads = atoi(env_threads);
    omp_set_num_threads(nthreads);

    printf("=== V37 Intertwined Ring: %s ===\n", geom_names[geom]);
    printf("Equations: 6-field Cosserat (3 phi + 3 theta)\n");
    printf("m²=%.4f, m_theta²=%.4f, eta=%.3f, mu=%.3f, kappa=%.1f, A_bg=%.2f\n",
           MASS2, MTHETA2, ETA, MU, KAPPA, A_BG);
    printf("N=%d L=%.0f T=%.0f\n", N, L, T);

    mkdir("data", 0755); mkdir(outdir, 0755);
    Grid *g = grid_alloc(N, L);
    printf("dx=%.4f dt=%.5f threads=%d\n\n", g->dx, g->dt, nthreads);

    /* Initialize based on geometry */
    switch (geom) {
        case GEOM_RING1:  init_ring1(g);  break;
        case GEOM_SPH:    init_sph(g);    break;
        case GEOM_BRAID3: init_braid3(g); break;
    }

    compute_forces(g);

    /* Save t=0 snapshot */
    save_field(g, 0, outdir);

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

    /* Time-averaged death check: need time corresponding to ~50 time units */
    double death_time_window = 50.0;
    int death_min_samples = (int)(death_time_window / diag_dt);
    if (death_min_samples < 20) death_min_samples = 20;
    if (death_min_samples > DEATH_WINDOW) death_min_samples = DEATH_WINDOW;

    for (int step = 0; step <= n_steps; step++) {
        if (step > 0) verlet_step(g);
        double t = step * g->dt;

        if (step % diag_every == 0) {
            double epk, etk, eg, em, ep, etg, etm, ec, et;
            compute_energy(g, &epk, &etk, &eg, &em, &ep, &etg, &etm, &ec, &et);
            double trms = theta_rms(g);
            double P_int = P_integrated(g);

            /* Inertia tensor at major intervals */
            int is_major = (step % (diag_every*50) == 0);
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

            /* Record for time-averaged death check */
            ep_record(ep);

            fprintf(fp, "%.2f\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%.4f\n",
                    t, epk, etk, eg, em, ep, etg, etm, ec, et, trms, P_int, aspect);
            fflush(fp);

            if (is_major) {
                double wall = omp_get_wtime() - wall0;
                double drift = 100.0*(et-E0)/(fabs(E0)+1e-30);
                double ep_avg = ep_average();
                printf("t=%7.1f E=%.3e (drift %+.3f%%) Ep=%.1f Ep_avg=%.1f P_int=%.3e theta_rms=%.3e aspect=%.2f [%.0f%% %.0fs]\n",
                       t, et, drift, ep, ep_avg, P_int, trms, aspect, 100.0*step/n_steps, wall);
                fflush(stdout);
            }

            /* TIME-AVERAGED death check (fixes breathing-mode false kills) */
            if (ep_count >= death_min_samples && fabs(Ep0) > 1e-10) {
                double ep_avg = ep_average();
                if (ep_avg < 0.01 * fabs(Ep0)) {
                    printf("\n  *** Time-averaged |E_pot| (%.3e) < 1%% of initial (%.3e) — braid DEAD. Stopping. ***\n",
                           ep_avg, fabs(Ep0));
                    save_field(g, t, outdir);
                    break;
                }
            }
        }

        if (step > 0 && step % snap_every == 0) {
            double t_snap = step * g->dt;
            save_field(g, t_snap, outdir);
        }
    }

    /* Final snapshot */
    double t_final = n_steps * g->dt;
    save_field(g, t_final, outdir);

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
        printf("Survival: %s\n", (fabs(ep) > 0.1*fabs(Ep0)) ? "YES (instantaneous)" : "NO (instantaneous)");
        double ep_avg = ep_average();
        printf("Survival (avg): %s (Ep_avg=%.3e vs Ep0=%.3e)\n",
               (ep_avg > 0.1*fabs(Ep0)) ? "YES" : "NO", ep_avg, Ep0);
    }

    fclose(fp);

    double wall = omp_get_wtime() - wall0;
    printf("\n=== Complete: %.0fs (%.1f min) ===\n", wall, wall/60);

    grid_free(g);
    return 0;
}
