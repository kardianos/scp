/*  collapse_3d.c — 3D three-field oscillon with configurable mass coupling
 *
 *  Three real scalar fields phi_0, phi_1, phi_2 on a 3D grid:
 *    d^2 phi_a/dt^2 = Lap(phi_a) - m_eff^2 * phi_a - dV/dphi_a
 *
 *  V(P) = (mu/2) P^2 / (1 + kappa_eff * P^2),  P = phi_0 * phi_1 * phi_2
 *
 *  Mass coupling modes (same as oscillon_1d.c):
 *    0: constant        m_eff^2 = m^2
 *    1: inverse         m_eff^2 = alpha / (1 + beta * Sigma)
 *    3: density kappa   m_eff^2 = m^2, kappa_eff = kappa0 / (1 + gamma * Sigma)
 *
 *  where Sigma = phi_0^2 + phi_1^2 + phi_2^2.
 *
 *  Initialization: Gaussian oscillon sphere (no seed file needed)
 *    phi_a(r) = A * exp(-r^2 / (2*sigma^2)) * cos(delta_a)
 *    vel_a(r) = 0
 *
 *  Build: gcc -O3 -march=native -fopenmp -o collapse_3d src/collapse_3d.c -lzstd -lm
 *  Run:   OMP_NUM_THREADS=8 ./collapse_3d -mode 1 -inv_alpha 1.0 -inv_beta 5 \
 *             -N 128 -L 13 -T 200 -dt_factor 0.025 -diag 2 -snap 5 \
 *             -damp_width 3 -damp_rate 0.01 -o data/inv_3d
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
static double MASS2   = 2.25;     /* m^2 for mode 0,3 */

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
   Grid
   ================================================================ */

typedef struct {
    double *mem;
    double *phi[NFIELDS];
    double *phi_vel[NFIELDS];
    double *phi_acc[NFIELDS];
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

    long total = 9 * g->N3;  /* 3 fields x (phi, vel, acc) */
    double bytes = total * sizeof(double);
    printf("Allocating %.2f GB (%ld doubles, N=%d, 3 fields)\n", bytes/1e9, total, N);
    g->mem = malloc(total * sizeof(double));
    if (!g->mem) { fprintf(stderr, "FATAL: malloc failed\n"); exit(1); }
    memset(g->mem, 0, total * sizeof(double));

    for (int a = 0; a < NFIELDS; a++) {
        g->phi[a]     = g->mem + (0 + a) * g->N3;
        g->phi_vel[a] = g->mem + (3 + a) * g->N3;
        g->phi_acc[a] = g->mem + (6 + a) * g->N3;
    }
    return g;
}

static void grid_free(Grid *g) { free(g->mem); free(g); }

/* ================================================================
   Initialize: Gaussian oscillon sphere
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
    /* skip theta (3 arrays) */
    fseek(fp, 3 * g->N3 * sizeof(double), SEEK_CUR);
    for (int a = 0; a < 3; a++) fread(g->phi_vel[a], sizeof(double), g->N3, fp);
    /* skip theta vel */
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
    const double A_BG = 0.1;

    printf("Init: plain braid (R_tube=%.1f, ellip=%.4f, A=%.1f)\n", R_tube, ellip, A[0]);
    printf("  kw=%.4f omega=%.4f delta={0, %.4f, %.4f}\n", kw, omega, delta[1], delta[2]);

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
                    g->phi_vel[a][idx] = 0.0;  /* start stationary */
                }
            }
        }
    }

    /* Print field statistics */
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
   Forces
   ================================================================ */

static void compute_forces(Grid *g) {
    const int N = g->N, NN = N * N;
    const long N3 = g->N3;
    const double idx2 = 1.0 / (g->dx * g->dx);

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
        double sigma_phi2 = p0*p0 + p1*p1 + p2*p2;
        double me2 = meff2_at(sigma_phi2);
        double keff = kappa_eff_at(sigma_phi2);

        double P = p0 * p1 * p2;
        double den = 1.0 + keff * P * P;
        double dVdP = MU * P / (den * den);

        for (int a = 0; a < NFIELDS; a++) {
            double lap = (g->phi[a][n_ip] + g->phi[a][n_im]
                        + g->phi[a][n_jp] + g->phi[a][n_jm]
                        + g->phi[a][n_kp] + g->phi[a][n_km]
                        - 6.0 * g->phi[a][idx]) * idx2;
            double dPda = (a==0) ? p1*p2 : (a==1) ? p0*p2 : p0*p1;

            g->phi_acc[a][idx] = lap - me2 * g->phi[a][idx] - dVdP * dPda;
        }
    }
}

/* ================================================================
   Absorbing boundary damping (spherical)
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

        /* Spherical absorbing boundary: damp based on distance from center */
        double r = sqrt(x*x + y*y + z*z);
        double R_damp = L - DAMP_WIDTH;  /* damping starts at this radius */
        double d_from_edge = L - r;      /* distance from outer sphere */

        if (r > R_damp && DAMP_WIDTH > 0 && DAMP_RATE > 0) {
            double s = (r - R_damp) / DAMP_WIDTH;  /* 0 at R_damp, 1 at L */
            if (s > 1.0) s = 1.0;
            double damp = 1.0 - DAMP_RATE * s * s;
            for (int a = 0; a < NFIELDS; a++) {
                g->phi_vel[a][idx] *= damp;
            }
        }
    }
}

/* ================================================================
   Verlet integrator
   ================================================================ */

static void verlet_step(Grid *g) {
    const long N3 = g->N3;
    const double hdt = 0.5 * g->dt, dt = g->dt;

    /* half-kick */
    for (int a = 0; a < NFIELDS; a++) {
        double *vp = g->phi_vel[a], *ap = g->phi_acc[a];
        #pragma omp parallel for schedule(static)
        for (long i = 0; i < N3; i++) vp[i] += hdt * ap[i];
    }
    /* drift */
    for (int a = 0; a < NFIELDS; a++) {
        double *pp = g->phi[a], *vp = g->phi_vel[a];
        #pragma omp parallel for schedule(static)
        for (long i = 0; i < N3; i++) pp[i] += dt * vp[i];
    }
    /* force */
    compute_forces(g);
    /* half-kick */
    for (int a = 0; a < NFIELDS; a++) {
        double *vp = g->phi_vel[a], *ap = g->phi_acc[a];
        #pragma omp parallel for schedule(static)
        for (long i = 0; i < N3; i++) vp[i] += hdt * ap[i];
    }

    apply_boundary_damping(g);
}

/* ================================================================
   Diagnostics
   ================================================================ */

static void compute_energy(Grid *g, double *E_kin, double *E_grad, double *E_mass,
                           double *E_pot, double *E_total,
                           double *phi_max_out, double *P_max_out,
                           double *meff_center_out) {
    const int N = g->N, NN = N*N;
    const long N3 = g->N3;
    const double dx = g->dx, dV = dx*dx*dx;
    const double idx1 = 1.0/(2.0*dx);
    double ek=0, eg=0, em=0, ep=0, phimax=0, pmax=0;

    #pragma omp parallel for reduction(+:ek,eg,em,ep) reduction(max:phimax,pmax) schedule(static)
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
            ek += 0.5*g->phi_vel[a][idx]*g->phi_vel[a][idx]*dV;
            double gx=(g->phi[a][n_ip]-g->phi[a][n_im])*idx1;
            double gy=(g->phi[a][n_jp]-g->phi[a][n_jm])*idx1;
            double gz=(g->phi[a][n_kp]-g->phi[a][n_km])*idx1;
            eg += 0.5*(gx*gx+gy*gy+gz*gz)*dV;
            em += 0.5*me2*g->phi[a][idx]*g->phi[a][idx]*dV;

            double ap = fabs(g->phi[a][idx]);
            if (ap > phimax) phimax = ap;
        }
        double P = p0 * p1 * p2;
        double P2 = P * P;
        ep += 0.5*MU*P2/(1.0+keff*P2)*dV;

        double Pabs = fabs(P);
        if (Pabs > pmax) pmax = Pabs;
    }
    *E_kin=ek; *E_grad=eg; *E_mass=em; *E_pot=ep;
    *E_total = ek+eg+em+ep;
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

/* |P| integrated (measures structure survival) */
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

/* Moment of inertia (to check if structure is 1 cluster) */
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
   SFA snapshot helper (write 3 phi columns)
   ================================================================ */

static void sfa_snap(SFA *sfa, Grid *g, double t) {
    float *p0 = to_f32(g->phi[0], g->N3);
    float *c0 = malloc(g->N3 * sizeof(float)); memcpy(c0, p0, g->N3*sizeof(float));
    float *p1 = to_f32(g->phi[1], g->N3);
    float *c1 = malloc(g->N3 * sizeof(float)); memcpy(c1, p1, g->N3*sizeof(float));
    float *p2 = to_f32(g->phi[2], g->N3);
    float *c2 = malloc(g->N3 * sizeof(float)); memcpy(c2, p2, g->N3*sizeof(float));
    void *cols[] = {c0, c1, c2};
    sfa_write_frame(sfa, t, cols);
    free(c0); free(c1); free(c2);
}

/* ================================================================
   Main
   ================================================================ */

int main(int argc, char **argv) {
    int N = 128;
    double L = 13.0, T = 200.0;
    double diag_dt = 2.0, snap_dt = 5.0;
    char outdir[256] = "data/inv_3d";
    char seedpath[256] = "";
    char initmode[32] = "oscillon";

    for (int i = 1; i < argc; i++) {
        if      (!strcmp(argv[i],"-N"))           N = atoi(argv[++i]);
        else if (!strcmp(argv[i],"-L"))           L = atof(argv[++i]);
        else if (!strcmp(argv[i],"-T"))           T = atof(argv[++i]);
        else if (!strcmp(argv[i],"-mode"))        MODE = atoi(argv[++i]);
        else if (!strcmp(argv[i],"-m"))         { double m = atof(argv[++i]); MASS2 = m*m; }
        else if (!strcmp(argv[i],"-mu"))          MU = atof(argv[++i]);
        else if (!strcmp(argv[i],"-kappa"))       KAPPA = atof(argv[++i]);
        else if (!strcmp(argv[i],"-inv_alpha"))   INV_ALPHA = atof(argv[++i]);
        else if (!strcmp(argv[i],"-inv_beta"))    INV_BETA = atof(argv[++i]);
        else if (!strcmp(argv[i],"-kappa_gamma")) KAPPA_GAMMA = atof(argv[++i]);
        else if (!strcmp(argv[i],"-A"))           A_INIT = atof(argv[++i]);
        else if (!strcmp(argv[i],"-sigma"))       SIGMA_W = atof(argv[++i]);
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
            fprintf(stderr, "Usage: %s [-mode 0-3] [-N 128] [-L 13] [-T 200] ...\n", argv[0]);
            return 1;
        }
    }

    int nthreads = 4;
    char *env_threads = getenv("OMP_NUM_THREADS");
    if (env_threads) nthreads = atoi(env_threads);
    omp_set_num_threads(nthreads);

    const char *mode_names[] = {"constant", "inverse", "unused", "density-kappa"};
    printf("=== collapse_3d: 3D Oscillon with Configurable Mass ===\n");
    printf("Mode: %d (%s)\n", MODE, mode_names[MODE]);
    printf("  mu=%.3f kappa=%.1f\n", MU, KAPPA);
    switch (MODE) {
    case 0: printf("  m^2=%.4f (m=%.4f)\n", MASS2, sqrt(MASS2)); break;
    case 1: printf("  inv_alpha=%.4f inv_beta=%.4f\n", INV_ALPHA, INV_BETA); break;
    case 3: printf("  m^2=%.4f kappa_gamma=%.4f\n", MASS2, KAPPA_GAMMA); break;
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

    /* SFA archive */
    char sfapath[512];
    snprintf(sfapath, sizeof(sfapath), "%s.sfa", outdir);
    SFA *sfa = sfa_create(sfapath, N, N, N, L, L, L, g->dt);
    sfa_add_column(sfa, "phi_x", SFA_F32, SFA_POSITION, 0);
    sfa_add_column(sfa, "phi_y", SFA_F32, SFA_POSITION, 1);
    sfa_add_column(sfa, "phi_z", SFA_F32, SFA_POSITION, 2);
    sfa_finalize_header(sfa);
    printf("SFA: %s (f32 columns)\n\n", sfapath);

    /* Write t=0 frame */
    sfa_snap(sfa, g, 0.0);

    /* Timeseries file */
    char tspath[512];
    snprintf(tspath, sizeof(tspath), "%s/timeseries.tsv", outdir);
    FILE *fp = fopen(tspath, "w");
    fprintf(fp, "t\tE_kin\tE_grad\tE_mass\tE_pot\tE_total\tphi_max\tP_max\tm_eff_center\tP_int\taspect\n");

    int n_steps = (int)(T / g->dt);
    int diag_every = (int)(diag_dt / g->dt); if (diag_every<1) diag_every=1;
    int snap_every = (int)(snap_dt / g->dt); if (snap_every<1) snap_every=1;

    printf("Steps=%d diag_every=%d snap_every=%d\n", n_steps, diag_every, snap_every);
    printf("Estimated time per step: dt=%.6f, total sim time: %.0f\n\n", g->dt, T);

    /* Initial energy diagnostic */
    double Ek0, Eg0, Em0, Ep0, Et0, pm0, Pm0, mc0;
    compute_energy(g, &Ek0, &Eg0, &Em0, &Ep0, &Et0, &pm0, &Pm0, &mc0);
    double Pint0 = P_integrated(g);
    printf("INITIAL: E_kin=%.4f E_grad=%.4f E_mass=%.4f E_pot=%.4f E_total=%.4f\n",
           Ek0, Eg0, Em0, Ep0, Et0);
    printf("  phi_max=%.4f P_max=%.4e P_int=%.4e m_eff_center=%.4f\n\n",
           pm0, Pm0, Pint0, mc0);

    fprintf(fp, "%.2f\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.4f\n",
            0.0, Ek0, Eg0, Em0, Ep0, Et0, pm0, Pm0, mc0, Pint0, 1.0);

    double wall0 = omp_get_wtime();
    double E0 = Et0;

    /* Key milestone times */
    double milestone_times[] = {50.0, 100.0, 150.0, 200.0};
    int n_milestones = 4;
    int milestone_idx = 0;

    for (int step = 1; step <= n_steps; step++) {
        verlet_step(g);
        double t = step * g->dt;

        /* SFA snapshot */
        if (step % snap_every == 0) {
            sfa_snap(sfa, g, t);
        }

        /* Diagnostics */
        if (step % diag_every == 0) {
            double ek, eg, em, ep, et, phimax, pmax, meff_c;
            compute_energy(g, &ek, &eg, &em, &ep, &et, &phimax, &pmax, &meff_c);
            double P_int = P_integrated(g);

            /* Aspect ratio at major intervals */
            double aspect = 1.0;
            int is_major = (step % (diag_every * 25) == 0);
            if (is_major) {
                double Ixx[3][3], evals[3], E_int;
                compute_inertia_tensor(g, Ixx, &E_int);
                eigen3x3(Ixx, evals);
                aspect = (evals[2] > 1e-30) ? evals[0]/evals[2] : 999;
            }

            fprintf(fp, "%.2f\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.4f\n",
                    t, ek, eg, em, ep, et, phimax, pmax, meff_c, P_int, aspect);
            fflush(fp);

            if (is_major) {
                double wall = omp_get_wtime() - wall0;
                double drift = 100.0*(et-E0)/(fabs(E0)+1e-30);
                double ms_per_step = (step > 0) ? 1000.0*wall/step : 0;
                printf("t=%7.1f E=%.4e (drift %+.3f%%) Ek=%.3f Eg=%.3f Em=%.3f Ep=%.3f "
                       "phi_max=%.4f P_max=%.3e m_c=%.4f asp=%.2f [%.0f%% %.1fs %.2fms/step]\n",
                       t, et, drift, ek, eg, em, ep, phimax, pmax, meff_c, aspect,
                       100.0*step/n_steps, wall, ms_per_step);
                fflush(stdout);
            }

            /* Check milestones */
            while (milestone_idx < n_milestones && t >= milestone_times[milestone_idx] - 0.5*g->dt) {
                printf("\n  >>> MILESTONE t=%.0f: E_pot=%.6f (initial=%.6f, ratio=%.4f) "
                       "phi_max=%.4f P_int=%.4e <<<\n\n",
                       milestone_times[milestone_idx], ep, Ep0,
                       (fabs(Ep0)>1e-30) ? ep/Ep0 : 0.0, phimax, P_int);
                milestone_idx++;
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
        double ek, eg, em, ep, et, phimax, pmax, meff_c;
        compute_energy(g, &ek, &eg, &em, &ep, &et, &phimax, &pmax, &meff_c);
        double P_int = P_integrated(g);
        double Ixx[3][3], evals[3], E_int;
        compute_inertia_tensor(g, Ixx, &E_int);
        eigen3x3(Ixx, evals);
        double aspect = (evals[2] > 1e-30) ? evals[0]/evals[2] : 999;

        double wall = omp_get_wtime() - wall0;
        double ms_per_step = (n_steps > 0) ? 1000.0*wall/n_steps : 0;

        printf("\n=== FINAL SUMMARY ===\n");
        printf("Mode: %d (%s)\n", MODE, mode_names[MODE]);
        if (MODE == 1) printf("  inv_alpha=%.4f inv_beta=%.4f\n", INV_ALPHA, INV_BETA);
        printf("E_total=%.6e (drift %.4f%%)\n", et, 100*(et-E0)/(fabs(E0)+1e-30));
        printf("E_kin=%.6f  E_grad=%.6f  E_mass=%.6f  E_pot=%.6f\n", ek, eg, em, ep);
        printf("E_pot ratio (final/initial): %.6f\n", (fabs(Ep0)>1e-30) ? ep/Ep0 : 0.0);
        printf("phi_max=%.6f (initial=%.6f, ratio=%.4f)\n",
               phimax, pm0, (pm0>1e-30) ? phimax/pm0 : 0.0);
        printf("P_max=%.6e  P_int=%.6e (initial=%.6e, ratio=%.4f)\n",
               pmax, P_int, Pint0, (Pint0>1e-30) ? P_int/Pint0 : 0.0);
        printf("m_eff_center=%.6f\n", meff_c);
        printf("Aspect ratio=%.3f (eigenvalues: %.1f %.1f %.1f)\n",
               aspect, evals[0], evals[1], evals[2]);
        printf("\nCollapse test: E_pot grew more negative? %s\n",
               (ep < Ep0*1.1) ? "YES" : "NO");
        printf("Amplitude grew? %s (phi_max ratio = %.4f)\n",
               (phimax > pm0*1.05) ? "YES" : "NO", (pm0>1e-30) ? phimax/pm0 : 0.0);
        printf("Stayed as 1 cluster? %s (aspect=%.3f)\n",
               (aspect < 3.0) ? "LIKELY" : "FRAGMENTED", aspect);

        printf("\nWall time: %.1fs (%.1f min)\n", wall, wall/60);
        printf("Performance: %.2f ms/step (%d steps)\n", ms_per_step, n_steps);
    }

    fclose(fp);
    if (f32_buf) free(f32_buf);
    grid_free(g);

    return 0;
}
