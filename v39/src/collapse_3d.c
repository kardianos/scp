/*  collapse_3d.c — 6-field Cosserat collapse with configurable mass coupling
 *
 *  Full 6-field equation: 3 position (phi) + 3 angle (theta)
 *    d^2 phi_a/dt^2 = Lap(phi_a) - m_eff^2 * phi_a - dV/dphi_a + eta*curl(theta)_a
 *    d^2 theta_a/dt^2 = Lap(theta_a) - m_theta^2 * theta_a + eta*curl(phi)_a
 *
 *  V(P) = (mu/2) P^2 / (1 + kappa_eff * P^2),  P = phi_0 * phi_1 * phi_2
 *
 *  Mass coupling modes:
 *    0: constant        m_eff^2 = m^2
 *    1: inverse         m_eff^2 = alpha / (1 + beta * Sigma)
 *    3: density kappa   m_eff^2 = m^2, kappa_eff = kappa0 / (1 + gamma * Sigma)
 *       Mode 3 includes the CORRECT Lagrangian force (Term I + Term II)
 *
 *  where Sigma = phi_0^2 + phi_1^2 + phi_2^2.
 *
 *  Build: gcc -O3 -march=native -fopenmp -o collapse_3d src/collapse_3d.c -lzstd -lm
 *  Run:   OMP_NUM_THREADS=8 ./collapse_3d -mode 3 -kappa_gamma 10 \
 *             -N 128 -L 10 -T 200 -dt_factor 0.025 -diag 2 -snap 5 \
 *             -damp_width 3 -damp_rate 0.01 -o data/kappa6f
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

/* ================================================================
   Physics parameters
   ================================================================ */

static double MU      = -41.345;
static double KAPPA   = 50.0;
static double MASS2   = 2.25;     /* m^2 for phi fields */
static double MTHETA2 = 0.0;     /* m_theta^2 for theta fields (0 = massless) */
static double ETA     = 0.5;     /* phi-theta coupling strength */

/* Mode 1: inverse coupling */
static double INV_ALPHA = 2.25;
static double INV_BETA  = 5.0;

/* Mode 3: density-dependent kappa */
static double KAPPA_GAMMA = 2.0;

/* Coupling mode */
static int MODE = 0;

/* Absorbing boundary */
static double DAMP_WIDTH = 3.0;
static double DAMP_RATE  = 0.01;

/* Initialization */
static double A_INIT  = 0.8;
static double SIGMA_W = 3.0;
static double DELTA[3] = {0.0, 3.0005, 4.4325};
static double A_BG    = 0.1;

/* Timestep factor: dt = dt_factor * dx */
static double DT_FACTOR = 0.025;

/* ================================================================
   Effective mass and kappa
   ================================================================ */

static inline double meff2_at(double sigma_phi2)
{
    switch (MODE) {
    case 0: return MASS2;
    case 1: return INV_ALPHA / (1.0 + INV_BETA * sigma_phi2);
    case 3: return MASS2;
    default: return MASS2;
    }
}

static inline double kappa_eff_at(double sigma_phi2)
{
    if (MODE == 3)
        return KAPPA / (1.0 + KAPPA_GAMMA * sigma_phi2);
    return KAPPA;
}

/* ================================================================
   Grid: 18 arrays (6 fields x {val, vel, acc})
   ================================================================ */

typedef struct {
    double *mem;
    /* Position fields */
    double *phi[NFIELDS];
    double *phi_vel[NFIELDS];
    double *phi_acc[NFIELDS];
    /* Angle fields */
    double *theta[NFIELDS];
    double *theta_vel[NFIELDS];
    double *theta_acc[NFIELDS];
    int N; long N3;
    double L, dx, dt;
} Grid;

static Grid *grid_alloc(int N, double L, double dt_factor) {
    Grid *g = calloc(1, sizeof(Grid));
    g->N  = N;
    g->N3 = (long)N * N * N;
    g->L  = L;
    g->dx = 2.0 * L / (N - 1);
    g->dt = dt_factor * g->dx;

    long total = 18 * g->N3;  /* 6 fields x (val, vel, acc) */
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
   Curl computation helper
   curl(F)_0 = dF2/dy - dF1/dz
   curl(F)_1 = dF0/dz - dF2/dx
   curl(F)_2 = dF1/dx - dF0/dy
   ================================================================ */

static inline double curl_component(double *F[3], int a,
    long n_ip, long n_im, long n_jp, long n_jm, long n_kp, long n_km,
    double idx1) {
    if (a == 0) return (F[2][n_jp] - F[2][n_jm] - F[1][n_kp] + F[1][n_km]) * idx1;
    if (a == 1) return (F[0][n_kp] - F[0][n_km] - F[2][n_ip] + F[2][n_im]) * idx1;
    return            (F[1][n_ip] - F[1][n_im] - F[0][n_jp] + F[0][n_jm]) * idx1;
}

/* ================================================================
   Initialize: various modes
   ================================================================ */

static void load_seed(Grid *g, const char *path) {
    FILE *fp = fopen(path, "rb");
    if (!fp) { fprintf(stderr, "FATAL: cannot open seed '%s'\n", path); exit(1); }
    int N_s; double L_s, t_s; int nf;
    fread(&N_s, sizeof(int), 1, fp);
    fread(&L_s, sizeof(double), 1, fp);
    fread(&t_s, sizeof(double), 1, fp);
    fread(&nf, sizeof(int), 1, fp);
    printf("Seed: %s (N=%d L=%.1f nf=%d)\n", path, N_s, L_s, nf);
    if (N_s != g->N) { fprintf(stderr, "FATAL: seed N=%d != grid N=%d\n", N_s, g->N); exit(1); }
    for (int a = 0; a < 3; a++) fread(g->phi[a], sizeof(double), g->N3, fp);
    if (nf >= 6) {
        for (int a = 0; a < 3; a++) fread(g->theta[a], sizeof(double), g->N3, fp);
    } else {
        fseek(fp, 0, SEEK_CUR); /* no theta in seed */
    }
    for (int a = 0; a < 3; a++) fread(g->phi_vel[a], sizeof(double), g->N3, fp);
    if (nf >= 6) {
        for (int a = 0; a < 3; a++) fread(g->theta_vel[a], sizeof(double), g->N3, fp);
    }
    fclose(fp);
    double mx = 0;
    for (long i = 0; i < g->N3; i++) if (fabs(g->phi[0][i]) > mx) mx = fabs(g->phi[0][i]);
    printf("  phi_0 max = %.4f\n", mx);
}

static void init_plain_braid(Grid *g) {
    const int N = g->N, NN = N * N;
    const double dx = g->dx, L = g->L;
    const double A[3] = {0.8, 0.8, 0.8};
    const double delta[3] = {0, 3.0005, 4.4325};
    const double R_tube = 3.0, ellip = 0.3325;
    const double kw = PI/L, omega = sqrt(kw*kw + MASS2);
    const double sx = 1+ellip, sy = 1-ellip;
    const double inv2R2 = 1.0/(2*R_tube*R_tube);
    const double k_bg = PI/L, omega_bg = sqrt(k_bg*k_bg + MASS2);

    printf("Init: plain braid (R_tube=%.1f, ellip=%.4f, A=%.1f, A_bg=%.2f)\n",
           R_tube, ellip, A[0], A_BG);

    for (int i = 0; i < N; i++) {
        double x = -L + i*dx;
        for (int j = 0; j < N; j++) {
            double y = -L + j*dx;
            for (int kk = 0; kk < N; kk++) {
                double z = -L + kk*dx;
                long idx = (long)i*NN + j*N + kk;
                double r2e = x*x/(sx*sx) + y*y/(sy*sy);
                double env = exp(-r2e * inv2R2);
                for (int a = 0; a < NFIELDS; a++) {
                    double ph = kw*z + delta[a];
                    double ph_bg = k_bg*z + 2*PI*a/3.0;
                    g->phi[a][idx] = A[a]*env*cos(ph) + A_BG*cos(ph_bg);
                    g->phi_vel[a][idx] = omega*A[a]*env*sin(ph) + omega_bg*A_BG*sin(ph_bg);
                }
                /* theta starts at zero — sourced by curl(phi) */
            }
        }
    }
}

static void init_oscillon(Grid *g) {
    const int N = g->N, NN = N * N;
    const double L = g->L, dx = g->dx;

    printf("Initializing oscillon sphere: A=%.3f sigma=%.3f\n", A_INIT, SIGMA_W);
    printf("  delta = {%.4f, %.4f, %.4f}\n", DELTA[0], DELTA[1], DELTA[2]);

    for (int i = 0; i < N; i++) {
        double x = -L + i * dx;
        for (int j = 0; j < N; j++) {
            double y = -L + j * dx;
            for (int k = 0; k < N; k++) {
                double z = -L + k * dx;
                long idx = (long)i * NN + j * N + k;
                double r2 = x*x + y*y + z*z;
                double env = A_INIT * exp(-r2 / (2.0 * SIGMA_W * SIGMA_W));
                for (int a = 0; a < NFIELDS; a++) {
                    g->phi[a][idx] = env * cos(DELTA[a]);
                    g->phi_vel[a][idx] = 0.0;
                }
                /* theta starts at zero */
            }
        }
    }

    for (int a = 0; a < NFIELDS; a++) {
        double mn = 1e30, mx = -1e30, rms = 0;
        for (long i = 0; i < g->N3; i++) {
            double v = g->phi[a][i];
            if (v < mn) mn = v;
            if (v > mx) mx = v;
            rms += v * v;
        }
        rms = sqrt(rms / g->N3);
        printf("  phi_%d: [%.4f, %.4f] rms=%.4f\n", a, mn, mx, rms);
    }
}

/* ================================================================
   Forces: 6-field Cosserat with configurable mass/kappa coupling
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

        /* --- Position field forces --- */
        double p0 = g->phi[0][idx], p1 = g->phi[1][idx], p2 = g->phi[2][idx];
        double sigma_phi2 = p0*p0 + p1*p1 + p2*p2;
        double me2 = meff2_at(sigma_phi2);
        double keff = kappa_eff_at(sigma_phi2);

        double P = p0 * p1 * p2;
        double P2 = P * P;
        double den = 1.0 + keff * P2;
        /* Term I: dV/dP = mu*P*(1+gamma*Sigma)^2 / D^2
           (encoded as mu*P/den^2, which is equivalent — see THEORY doc §2.3) */
        double dVdP = MU * P / (den * den);

        /* Term II (mode 3 only): density feedback force
           = -mu*gamma*kappa_0*P^4/D^2 * phi_a
           where D = 1 + gamma*Sigma + kappa_0*P^2 */
        double term2_coeff = 0.0;
        if (MODE == 3 && KAPPA_GAMMA > 0.0) {
            double D = 1.0 + KAPPA_GAMMA * sigma_phi2 + KAPPA * P2;
            term2_coeff = MU * KAPPA_GAMMA * KAPPA * P2 * P2 / (D * D);
        }

        for (int a = 0; a < NFIELDS; a++) {
            /* Laplacian */
            double lap = (g->phi[a][n_ip] + g->phi[a][n_im]
                        + g->phi[a][n_jp] + g->phi[a][n_jm]
                        + g->phi[a][n_kp] + g->phi[a][n_km]
                        - 6.0 * g->phi[a][idx]) * idx2;
            double dPda = (a==0) ? p1*p2 : (a==1) ? p0*p2 : p0*p1;

            /* Coupling: + eta * curl(theta)_a */
            double curl_theta = curl_component(g->theta, a,
                n_ip, n_im, n_jp, n_jm, n_kp, n_km, idx1);

            g->phi_acc[a][idx] = lap - me2 * g->phi[a][idx]
                               - dVdP * dPda - term2_coeff * g->phi[a][idx]
                               + ETA * curl_theta;
        }

        /* --- Angle field forces --- */
        for (int a = 0; a < NFIELDS; a++) {
            double lap_t = (g->theta[a][n_ip] + g->theta[a][n_im]
                          + g->theta[a][n_jp] + g->theta[a][n_jm]
                          + g->theta[a][n_kp] + g->theta[a][n_km]
                          - 6.0 * g->theta[a][idx]) * idx2;

            /* Coupling: + eta * curl(phi)_a */
            double curl_phi = curl_component(g->phi, a,
                n_ip, n_im, n_jp, n_jm, n_kp, n_km, idx1);

            g->theta_acc[a][idx] = lap_t - MTHETA2 * g->theta[a][idx]
                                 + ETA * curl_phi;
        }
    }
}

/* ================================================================
   Absorbing boundary damping (spherical) — all 6 fields
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

        double r = sqrt(x*x + y*y + z*z);
        double R_damp = L - DAMP_WIDTH;

        if (r > R_damp && DAMP_WIDTH > 0 && DAMP_RATE > 0) {
            double s = (r - R_damp) / DAMP_WIDTH;
            if (s > 1.0) s = 1.0;
            double damp = 1.0 - DAMP_RATE * s * s;
            for (int a = 0; a < NFIELDS; a++) {
                g->phi_vel[a][idx] *= damp;
                g->theta_vel[a][idx] *= damp;
            }
        }
    }
}

/* ================================================================
   Symplectic Velocity Verlet (all 6 fields)
   ================================================================ */

static void verlet_step(Grid *g) {
    const long N3 = g->N3;
    const double hdt = 0.5 * g->dt, dt = g->dt;

    /* Half-kick: all velocities */
    for (int a = 0; a < NFIELDS; a++) {
        double *vp = g->phi_vel[a], *ap = g->phi_acc[a];
        double *vt = g->theta_vel[a], *at = g->theta_acc[a];
        #pragma omp parallel for schedule(static)
        for (long i = 0; i < N3; i++) { vp[i] += hdt*ap[i]; vt[i] += hdt*at[i]; }
    }
    /* Drift: all fields */
    for (int a = 0; a < NFIELDS; a++) {
        double *pp = g->phi[a], *vp = g->phi_vel[a];
        double *pt = g->theta[a], *vt = g->theta_vel[a];
        #pragma omp parallel for schedule(static)
        for (long i = 0; i < N3; i++) { pp[i] += dt*vp[i]; pt[i] += dt*vt[i]; }
    }
    /* Recompute forces */
    compute_forces(g);
    /* Half-kick */
    for (int a = 0; a < NFIELDS; a++) {
        double *vp = g->phi_vel[a], *ap = g->phi_acc[a];
        double *vt = g->theta_vel[a], *at = g->theta_acc[a];
        #pragma omp parallel for schedule(static)
        for (long i = 0; i < N3; i++) { vp[i] += hdt*ap[i]; vt[i] += hdt*at[i]; }
    }

    apply_boundary_damping(g);
}

/* ================================================================
   Diagnostics
   ================================================================ */

static void compute_energy(Grid *g, double *E_phi_kin, double *E_theta_kin,
                           double *E_grad, double *E_mass, double *E_pot,
                           double *E_theta_grad, double *E_theta_mass,
                           double *E_coupling, double *E_total,
                           double *phi_max_out, double *P_max_out,
                           double *meff_center_out) {
    const int N = g->N, NN = N*N;
    const long N3 = g->N3;
    const double dx = g->dx, dV = dx*dx*dx;
    const double idx1 = 1.0/(2.0*dx);
    double epk=0, etk=0, eg=0, em=0, ep=0, etg=0, etm=0, ec=0;
    double phimax=0, pmax=0;

    #pragma omp parallel for reduction(+:epk,etk,eg,em,ep,etg,etm,ec) \
        reduction(max:phimax,pmax) schedule(static)
    for (long idx = 0; idx < N3; idx++) {
        int i=(int)(idx/NN), j=(int)((idx/N)%N), k=(int)(idx%N);
        int ip=(i+1)%N, im=(i-1+N)%N;
        int jp=(j+1)%N, jm=(j-1+N)%N;
        int kp=(k+1)%N, km=(k-1+N)%N;
        long n_ip=(long)ip*NN+j*N+k, n_im=(long)im*NN+j*N+k;
        long n_jp=(long)i*NN+jp*N+k, n_jm=(long)i*NN+jm*N+k;
        long n_kp=(long)i*NN+j*N+kp, n_km=(long)i*NN+j*N+km;

        double p0 = g->phi[0][idx], p1 = g->phi[1][idx], p2 = g->phi[2][idx];
        double sigma_phi2 = p0*p0 + p1*p1 + p2*p2;
        double me2 = meff2_at(sigma_phi2);
        double keff = kappa_eff_at(sigma_phi2);

        for (int a = 0; a < NFIELDS; a++) {
            /* Phi kinetic */
            epk += 0.5*g->phi_vel[a][idx]*g->phi_vel[a][idx]*dV;
            /* Theta kinetic */
            etk += 0.5*g->theta_vel[a][idx]*g->theta_vel[a][idx]*dV;
            /* Phi gradient */
            double gx=(g->phi[a][n_ip]-g->phi[a][n_im])*idx1;
            double gy=(g->phi[a][n_jp]-g->phi[a][n_jm])*idx1;
            double gz=(g->phi[a][n_kp]-g->phi[a][n_km])*idx1;
            eg += 0.5*(gx*gx+gy*gy+gz*gz)*dV;
            /* Phi mass */
            em += 0.5*me2*g->phi[a][idx]*g->phi[a][idx]*dV;
            /* Theta gradient */
            double tgx=(g->theta[a][n_ip]-g->theta[a][n_im])*idx1;
            double tgy=(g->theta[a][n_jp]-g->theta[a][n_jm])*idx1;
            double tgz=(g->theta[a][n_kp]-g->theta[a][n_km])*idx1;
            etg += 0.5*(tgx*tgx+tgy*tgy+tgz*tgz)*dV;
            /* Theta mass */
            etm += 0.5*MTHETA2*g->theta[a][idx]*g->theta[a][idx]*dV;

            double ap = fabs(g->phi[a][idx]);
            if (ap > phimax) phimax = ap;
        }
        /* Triple product potential (with effective kappa) */
        double P = p0*p1*p2;
        double P2 = P*P;
        if (MODE == 3) {
            /* V = (mu/2)*P^2*(1+gamma*Sigma)/(1+gamma*Sigma+kappa_0*P^2) */
            double D = 1.0 + KAPPA_GAMMA*sigma_phi2 + KAPPA*P2;
            ep += (MU/2.0)*P2*(1.0 + KAPPA_GAMMA*sigma_phi2)/D * dV;
        } else {
            ep += (MU/2.0)*P2/(1.0+keff*P2)*dV;
        }

        double Pabs = fabs(P);
        if (Pabs > pmax) pmax = Pabs;

        /* Coupling energy: -eta * Sum_a phi_a * curl(theta)_a */
        for (int a = 0; a < NFIELDS; a++) {
            double ct = curl_component(g->theta, a,
                n_ip, n_im, n_jp, n_jm, n_kp, n_km, idx1);
            ec -= ETA * g->phi[a][idx] * ct * dV;
        }
    }
    *E_phi_kin=epk; *E_theta_kin=etk; *E_grad=eg; *E_mass=em;
    *E_pot=ep; *E_theta_grad=etg; *E_theta_mass=etm; *E_coupling=ec;
    *E_total = epk+etk+eg+em+ep+etg+etm+ec;
    *phi_max_out = phimax;
    *P_max_out = pmax;

    /* m_eff at grid center */
    int ic = N/2;
    long center_idx = (long)ic*NN + ic*N + ic;
    double sig_c = g->phi[0][center_idx]*g->phi[0][center_idx]
                 + g->phi[1][center_idx]*g->phi[1][center_idx]
                 + g->phi[2][center_idx]*g->phi[2][center_idx];
    *meff_center_out = sqrt(fabs(meff2_at(sig_c)));
}

/* |P| integrated */
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

/* theta_rms — verification that angle fields are active */
static double theta_rms(Grid *g) {
    double sum = 0;
    #pragma omp parallel for reduction(+:sum)
    for (long i = 0; i < g->N3; i++) {
        for (int a = 0; a < NFIELDS; a++)
            sum += g->theta[a][i] * g->theta[a][i];
    }
    return sqrt(sum / (3.0 * g->N3));
}

/* Moment of inertia tensor */
static void compute_inertia_tensor(Grid *g, double Ixx[3][3], double *E_tot_out) {
    const int N = g->N, NN = N*N;
    const long N3 = g->N3;
    const double dx = g->dx, L = g->L;
    const double dV = dx*dx*dx;

    double cm0=0, cm1=0, cm2=0, E_sum=0;

    #pragma omp parallel for reduction(+:cm0,cm1,cm2,E_sum) schedule(static)
    for (long idx = 0; idx < N3; idx++) {
        int i=(int)(idx/NN), j=(int)((idx/N)%N), k=(int)(idx%N);
        double x = -L + i*dx, y = -L + j*dx, z = -L + k*dx;
        double w = 0;
        for (int a = 0; a < NFIELDS; a++)
            w += g->phi[a][idx] * g->phi[a][idx];
        w *= dV;
        cm0 += x*w; cm1 += y*w; cm2 += z*w;
        E_sum += w;
    }
    if (E_sum > 0) { cm0/=E_sum; cm1/=E_sum; cm2/=E_sum; }
    *E_tot_out = E_sum;

    double I00=0,I01=0,I02=0,I11=0,I12=0,I22=0;
    #pragma omp parallel for reduction(+:I00,I01,I02,I11,I12,I22) schedule(static)
    for (long idx = 0; idx < N3; idx++) {
        int i=(int)(idx/NN), j=(int)((idx/N)%N), k=(int)(idx%N);
        double x=-L+i*dx-cm0, y=-L+j*dx-cm1, z=-L+k*dx-cm2;
        double w=0;
        for (int a=0;a<NFIELDS;a++) w += g->phi[a][idx]*g->phi[a][idx];
        w *= dV;
        double r2 = x*x+y*y+z*z;
        I00+=(r2-x*x)*w; I11+=(r2-y*y)*w; I22+=(r2-z*z)*w;
        I01-=x*y*w; I02-=x*z*w; I12-=y*z*w;
    }
    Ixx[0][0]=I00; Ixx[0][1]=I01; Ixx[0][2]=I02;
    Ixx[1][0]=I01; Ixx[1][1]=I11; Ixx[1][2]=I12;
    Ixx[2][0]=I02; Ixx[2][1]=I12; Ixx[2][2]=I22;
}

/* Eigenvalues of 3x3 symmetric matrix */
static void eigen3x3(double A[3][3], double evals[3]) {
    double p1 = A[0][1]*A[0][1] + A[0][2]*A[0][2] + A[1][2]*A[1][2];
    double q = (A[0][0] + A[1][1] + A[2][2]) / 3.0;
    double p2 = (A[0][0]-q)*(A[0][0]-q) + (A[1][1]-q)*(A[1][1]-q)
              + (A[2][2]-q)*(A[2][2]-q) + 2*p1;
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
    for (int i=0;i<2;i++) for (int j=i+1;j<3;j++)
        if (evals[j] > evals[i]) { double t=evals[i]; evals[i]=evals[j]; evals[j]=t; }
}

/* ================================================================
   f32 downcast buffer for SFA output
   ================================================================ */

static float *f32_buf = NULL;
static long f32_buf_size = 0;

static float *to_f32(double *src, long n) {
    if (n > f32_buf_size) {
        f32_buf_size = n;
        f32_buf = realloc(f32_buf, n * sizeof(float));
    }
    for (long i = 0; i < n; i++)
        f32_buf[i] = (float)src[i];
    return f32_buf;
}

/* ================================================================
   SFA snapshot helper (6 columns: phi_x/y/z + theta_x/y/z)
   ================================================================ */

static void sfa_snap(SFA *sfa, Grid *g, double t) {
    const long n = g->N3;
    float *cols[6];
    for (int i = 0; i < 6; i++)
        cols[i] = malloc(n * sizeof(float));

    /* Convert phi fields */
    for (int a = 0; a < 3; a++) {
        to_f32(g->phi[a], n);
        memcpy(cols[a], f32_buf, n * sizeof(float));
    }
    /* Convert theta fields */
    for (int a = 0; a < 3; a++) {
        to_f32(g->theta[a], n);
        memcpy(cols[3+a], f32_buf, n * sizeof(float));
    }

    sfa_write_frame(sfa, t, (void**)cols);

    for (int i = 0; i < 6; i++)
        free(cols[i]);
}

/* ================================================================
   Main
   ================================================================ */

int main(int argc, char **argv) {
    int N = 128;
    double L = 10.0, T = 200.0;
    double diag_dt = 2.0, snap_dt = 5.0;
    char outdir[256] = "data/kappa6f";
    char seedpath[256] = "";
    char initmode[32] = "oscillon";

    for (int i = 1; i < argc; i++) {
        if      (!strcmp(argv[i],"-N"))           N = atoi(argv[++i]);
        else if (!strcmp(argv[i],"-L"))           L = atof(argv[++i]);
        else if (!strcmp(argv[i],"-T"))           T = atof(argv[++i]);
        else if (!strcmp(argv[i],"-mode"))        MODE = atoi(argv[++i]);
        else if (!strcmp(argv[i],"-m"))         { double m = atof(argv[++i]); MASS2 = m*m; }
        else if (!strcmp(argv[i],"-mt"))        { double m = atof(argv[++i]); MTHETA2 = m*m; }
        else if (!strcmp(argv[i],"-eta"))         ETA = atof(argv[++i]);
        else if (!strcmp(argv[i],"-mu"))          MU = atof(argv[++i]);
        else if (!strcmp(argv[i],"-kappa"))       KAPPA = atof(argv[++i]);
        else if (!strcmp(argv[i],"-inv_alpha"))   INV_ALPHA = atof(argv[++i]);
        else if (!strcmp(argv[i],"-inv_beta"))    INV_BETA = atof(argv[++i]);
        else if (!strcmp(argv[i],"-kappa_gamma")) KAPPA_GAMMA = atof(argv[++i]);
        else if (!strcmp(argv[i],"-A"))           A_INIT = atof(argv[++i]);
        else if (!strcmp(argv[i],"-sigma"))       SIGMA_W = atof(argv[++i]);
        else if (!strcmp(argv[i],"-bg"))          A_BG = atof(argv[++i]);
        else if (!strcmp(argv[i],"-dt_factor"))   DT_FACTOR = atof(argv[++i]);
        else if (!strcmp(argv[i],"-diag"))        diag_dt = atof(argv[++i]);
        else if (!strcmp(argv[i],"-snap"))        snap_dt = atof(argv[++i]);
        else if (!strcmp(argv[i],"-o"))           strncpy(outdir, argv[++i], 255);
        else if (!strcmp(argv[i],"-damp_width"))  DAMP_WIDTH = atof(argv[++i]);
        else if (!strcmp(argv[i],"-damp_rate"))   DAMP_RATE = atof(argv[++i]);
        else if (!strcmp(argv[i],"-seed"))        strncpy(seedpath, argv[++i], 255);
        else if (!strcmp(argv[i],"-init"))        strncpy(initmode, argv[++i], 31);
        else {
            fprintf(stderr, "Unknown option: %s\n", argv[i]);
            fprintf(stderr, "Usage: %s [-mode 0-3] [-N 128] [-L 10] [-T 200] "
                    "[-eta 0.5] [-mt 0] ...\n", argv[0]);
            return 1;
        }
    }

    int nthreads = 4;
    char *env_threads = getenv("OMP_NUM_THREADS");
    if (env_threads) nthreads = atoi(env_threads);
    omp_set_num_threads(nthreads);

    const char *mode_names[] = {"constant", "inverse", "unused", "density-kappa"};
    printf("=== collapse_3d: 6-field Cosserat + Configurable Mass ===\n");
    printf("d²phi_a/dt² = Lap(phi_a) - m²phi_a - dV/dphi_a + eta*curl(theta)_a\n");
    printf("d²theta_a/dt² = Lap(theta_a) - m_t²theta_a + eta*curl(phi)_a\n\n");
    printf("Mode: %d (%s)\n", MODE, mode_names[MODE]);
    printf("  mu=%.3f kappa=%.1f m²=%.4f m_theta²=%.4f eta=%.3f\n",
           MU, KAPPA, MASS2, MTHETA2, ETA);
    switch (MODE) {
    case 1: printf("  inv_alpha=%.4f inv_beta=%.4f\n", INV_ALPHA, INV_BETA); break;
    case 3: printf("  kappa_gamma=%.4f (Term II density feedback: ON)\n", KAPPA_GAMMA); break;
    }
    printf("Absorbing BC: damp_width=%.1f damp_rate=%.4f\n", DAMP_WIDTH, DAMP_RATE);
    printf("N=%d L=%.1f T=%.0f dt_factor=%.4f\n", N, L, T, DT_FACTOR);

    mkdir("data", 0755); mkdir(outdir, 0755);
    Grid *g = grid_alloc(N, L, DT_FACTOR);
    printf("dx=%.4f dt=%.6f threads=%d\n\n", g->dx, g->dt, nthreads);

    /* Initialize */
    if (seedpath[0]) {
        load_seed(g, seedpath);
    } else if (!strcmp(initmode, "braid")) {
        init_plain_braid(g);
    } else {
        init_oscillon(g);
    }

    /* Compute initial forces */
    compute_forces(g);

    /* SFA archive: 6 columns (compressed f32) */
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
    printf("SFA: %s (6 f32 columns, BSS+zstd compressed)\n\n", sfapath);

    /* Write t=0 frame */
    sfa_snap(sfa, g, 0.0);

    /* Timeseries file */
    char tspath[512];
    snprintf(tspath, sizeof(tspath), "%s/timeseries.tsv", outdir);
    FILE *fp = fopen(tspath, "w");
    fprintf(fp, "t\tE_phi_kin\tE_theta_kin\tE_grad\tE_mass\tE_pot\tE_tgrad\tE_tmass\t"
                "E_coupling\tE_total\tphi_max\tP_max\tm_eff_c\tP_int\ttheta_rms\taspect\n");

    int n_steps = (int)(T / g->dt);
    int diag_every = (int)(diag_dt / g->dt); if (diag_every<1) diag_every=1;
    int snap_every = (int)(snap_dt / g->dt); if (snap_every<1) snap_every=1;

    printf("Steps=%d diag_every=%d snap_every=%d\n", n_steps, diag_every, snap_every);

    /* Initial energy diagnostic */
    double epk0, etk0, eg0, em0, ep0, etg0, etm0, ec0, et0, pm0, Pm0, mc0;
    compute_energy(g, &epk0, &etk0, &eg0, &em0, &ep0, &etg0, &etm0, &ec0, &et0,
                   &pm0, &Pm0, &mc0);
    double Pint0 = P_integrated(g);
    double trms0 = theta_rms(g);
    printf("INITIAL: E_total=%.4e  E_pot=%.4f  phi_max=%.4f  P_max=%.4e  "
           "P_int=%.4e  theta_rms=%.3e\n\n", et0, ep0, pm0, Pm0, Pint0, trms0);

    fprintf(fp, "%.2f\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t"
                "%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.4f\n",
            0.0, epk0, etk0, eg0, em0, ep0, etg0, etm0, ec0, et0,
            pm0, Pm0, mc0, Pint0, trms0, 1.0);

    double wall0 = omp_get_wtime();
    double E0 = et0;

    for (int step = 1; step <= n_steps; step++) {
        verlet_step(g);
        double t = step * g->dt;

        /* SFA snapshot */
        if (step % snap_every == 0)
            sfa_snap(sfa, g, t);

        /* Diagnostics */
        if (step % diag_every == 0) {
            double epk, etk, eg, em, ep, etg, etm, ec, et, phimax, pmax, meff_c;
            compute_energy(g, &epk, &etk, &eg, &em, &ep, &etg, &etm, &ec, &et,
                           &phimax, &pmax, &meff_c);
            double P_int = P_integrated(g);
            double trms = theta_rms(g);

            /* Aspect ratio at major intervals */
            double aspect = 1.0;
            int is_major = (step % (diag_every * 25) == 0);
            if (is_major) {
                double Ixx[3][3], evals[3], E_int;
                compute_inertia_tensor(g, Ixx, &E_int);
                eigen3x3(Ixx, evals);
                aspect = (evals[2] > 1e-30) ? evals[0]/evals[2] : 999;
            }

            fprintf(fp, "%.2f\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t"
                        "%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.4f\n",
                    t, epk, etk, eg, em, ep, etg, etm, ec, et,
                    phimax, pmax, meff_c, P_int, trms, aspect);
            fflush(fp);

            if (is_major) {
                double wall = omp_get_wtime() - wall0;
                double drift = 100.0*(et-E0)/(fabs(E0)+1e-30);
                double ms_per_step = 1000.0*wall/step;
                printf("t=%7.1f E=%.4e (drift %+.3f%%) Ep=%.1f phi=%.4f "
                       "P_max=%.3e theta_rms=%.3e asp=%.2f [%.0f%% %.1fs %.2fms/step]\n",
                       t, et, drift, ep, phimax, pmax, trms, aspect,
                       100.0*step/n_steps, wall, ms_per_step);
                fflush(stdout);
            }
        }
    }

    /* Final frame */
    double t_final = n_steps * g->dt;
    sfa_snap(sfa, g, t_final);

    uint32_t sfa_nframes = sfa->total_frames;
    sfa_close(sfa);
    printf("\nSFA closed: %s (%u frames)\n", sfapath, sfa_nframes);

    /* Final summary */
    {
        double epk, etk, eg, em, ep, etg, etm, ec, et, phimax, pmax, meff_c;
        compute_energy(g, &epk, &etk, &eg, &em, &ep, &etg, &etm, &ec, &et,
                       &phimax, &pmax, &meff_c);
        double P_int = P_integrated(g);
        double trms = theta_rms(g);
        double Ixx[3][3], evals[3], E_int;
        compute_inertia_tensor(g, Ixx, &E_int);
        eigen3x3(Ixx, evals);
        double aspect = (evals[2] > 1e-30) ? evals[0]/evals[2] : 999;

        double wall = omp_get_wtime() - wall0;
        double ms_per_step = (n_steps > 0) ? 1000.0*wall/n_steps : 0;

        printf("\n=== FINAL SUMMARY (6-field Cosserat) ===\n");
        printf("Mode: %d (%s)\n", MODE, mode_names[MODE]);
        if (MODE == 3) printf("  kappa_gamma=%.4f (Term II: ON)\n", KAPPA_GAMMA);
        printf("  eta=%.3f  m²=%.4f  m_theta²=%.4f\n", ETA, MASS2, MTHETA2);
        printf("E_total=%.6e (drift %.4f%%)\n", et, 100*(et-E0)/(fabs(E0)+1e-30));
        printf("  E_phi_kin=%.4f  E_theta_kin=%.4f  E_coupling=%.4f\n", epk, etk, ec);
        printf("  E_grad=%.4f  E_mass=%.4f  E_pot=%.4f\n", eg, em, ep);
        printf("  E_theta_grad=%.4f  E_theta_mass=%.4f\n", etg, etm);
        printf("E_pot ratio (final/initial): %.6f\n", (fabs(ep0)>1e-30) ? ep/ep0 : 0.0);
        printf("phi_max=%.6f (initial=%.6f, ratio=%.4f)\n",
               phimax, pm0, (pm0>1e-30) ? phimax/pm0 : 0.0);
        printf("P_max=%.6e  P_int=%.6e (initial=%.6e, ratio=%.4f)\n",
               pmax, P_int, Pint0, (Pint0>1e-30) ? P_int/Pint0 : 0.0);
        printf("theta_rms=%.6e (initial=%.6e)\n", trms, trms0);
        printf("Aspect=%.3f (eigenvalues: %.1f %.1f %.1f)\n",
               aspect, evals[0], evals[1], evals[2]);
        printf("\nVerification checklist:\n");
        printf("  [%s] 18 arrays allocated (6 fields x 3)\n", "OK");
        printf("  [%s] theta_rms grew from zero: %.3e -> %.3e\n",
               (trms > trms0 + 1e-10) ? "OK" : "FAIL", trms0, trms);
        printf("  [%s] Curl coupling active (eta=%.3f)\n",
               (ETA > 0) ? "OK" : "WARN", ETA);
        printf("  [%s] Collapse: E_pot more negative? %s\n",
               (ep < ep0*1.1) ? "OK" : "--",
               (ep < ep0*1.1) ? "YES" : "NO");

        printf("\nWall time: %.1fs (%.1f min)\n", wall, wall/60);
        printf("Performance: %.2f ms/step (%d steps)\n", ms_per_step, n_steps);
    }

    fclose(fp);
    if (f32_buf) free(f32_buf);
    grid_free(g);

    return 0;
}
