/*  v37_seedrun.c — Seed-initialized Cosserat simulation (6-field)
 *
 *  Reads a binary seed file (from seed_to_c_init.py) as initial condition,
 *  then runs the standard 6-field Cosserat equations with absorbing BC,
 *  SFA output (f32 columns), and inline fragmentation detection.
 *
 *  Physics identical to v37_crossed_v2.c:
 *    d^2 phi_a/dt^2 = Lap(phi_a) - m^2 phi_a - V'(P) + eta * curl(theta)_a
 *    d^2 theta_a/dt^2 = Lap(theta_a) - m_theta^2 theta_a + eta * curl(phi)_a
 *    V(P) = (mu/2) P^2 / (1 + kappa P^2),  P = phi_0 phi_1 phi_2
 *
 *  Build: gcc -O3 -march=native -fopenmp -o v37_seedrun src/v37_seedrun.c -lzstd -lm
 *  Run:   ./v37_seedrun -seed init.bin -N 128 -L 15 -T 500 -snap 5 \
 *             -damp_width 4 -damp_rate 0.005 -o data/seed_test
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
static double MASS2   = 2.25;     /* position field mass^2 */
static double MTHETA2 = 0.0;     /* angle field mass^2 (default massless) */
static double ETA     = 0.5;     /* coupling strength */

/* Absorbing boundary */
static double DAMP_WIDTH = 4.0;  /* code units from edge */
static double DAMP_RATE  = 0.005; /* per timestep */

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
   Seed file loader
   ================================================================ */

static void load_seed(Grid *g, const char *path) {
    FILE *fp = fopen(path, "rb");
    if (!fp) {
        fprintf(stderr, "FATAL: cannot open seed file '%s'\n", path);
        exit(1);
    }

    int N_seed; double L_seed, t_seed; int nf;
    fread(&N_seed, sizeof(int), 1, fp);
    fread(&L_seed, sizeof(double), 1, fp);
    fread(&t_seed, sizeof(double), 1, fp);
    fread(&nf, sizeof(int), 1, fp);

    printf("Seed file: %s\n", path);
    printf("  N_seed=%d, L_seed=%.2f, t_seed=%.4f, nfields=%d\n",
           N_seed, L_seed, t_seed, nf);

    if (N_seed != g->N) {
        fprintf(stderr, "FATAL: seed N=%d does not match grid N=%d\n", N_seed, g->N);
        fprintf(stderr, "  Upscale in Python first: python seed_to_c_init.py ... --N_out %d\n", g->N);
        fclose(fp);
        exit(1);
    }
    if (nf != 12) {
        fprintf(stderr, "FATAL: seed nfields=%d, expected 12\n", nf);
        fclose(fp);
        exit(1);
    }

    /* Read 12 arrays: phi0-2, theta0-2, phivel0-2, thetavel0-2 */
    for (int a = 0; a < 3; a++) fread(g->phi[a], sizeof(double), g->N3, fp);
    for (int a = 0; a < 3; a++) fread(g->theta[a], sizeof(double), g->N3, fp);
    for (int a = 0; a < 3; a++) fread(g->phi_vel[a], sizeof(double), g->N3, fp);
    for (int a = 0; a < 3; a++) fread(g->theta_vel[a], sizeof(double), g->N3, fp);

    fclose(fp);

    /* Print field statistics */
    for (int a = 0; a < 3; a++) {
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
    for (int a = 0; a < 3; a++) {
        double mn = 1e30, mx = -1e30, rms = 0;
        for (long i = 0; i < g->N3; i++) {
            double v = g->phi_vel[a][i];
            if (v < mn) mn = v;
            if (v > mx) mx = v;
            rms += v * v;
        }
        rms = sqrt(rms / g->N3);
        printf("  vel_%d: [%.4f, %.4f] rms=%.4f\n", a, mn, mx, rms);
    }
    for (int a = 0; a < 3; a++) {
        double rms = 0;
        for (long i = 0; i < g->N3; i++)
            rms += g->theta[a][i] * g->theta[a][i];
        rms = sqrt(rms / g->N3);
        printf("  theta_%d rms=%.6f\n", a, rms);
    }
}

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

/* ================================================================
   Fragmentation detection
   ================================================================ */

/*  Flood-fill cluster count on |P| > threshold.
 *  Uses a coarsened grid (stride) for speed at major diagnostic intervals.
 *  Returns number of connected clusters.
 */
static int detect_clusters(Grid *g, double threshold, int stride) {
    const int N = g->N, NN = N*N;
    int Nc = (N + stride - 1) / stride;  /* coarsened grid size */
    long Nc3 = (long)Nc * Nc * Nc;

    /* Build above-threshold mask on coarsened grid */
    int *mask = calloc(Nc3, sizeof(int));   /* 0=below, 1=above, -1=visited */
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

    /* Flood-fill using a stack */
    int n_clusters = 0;
    long *stack = malloc(Nc3 * sizeof(long));

    for (long cidx = 0; cidx < Nc3; cidx++) {
        if (mask[cidx] != 1) continue;
        /* New cluster found: flood fill */
        n_clusters++;
        int sp = 0;
        stack[sp++] = cidx;
        mask[cidx] = -1;  /* mark visited */

        while (sp > 0) {
            long cur = stack[--sp];
            int ci = (int)(cur / (Nc*Nc));
            int cj = (int)((cur / Nc) % Nc);
            int ck = (int)(cur % Nc);

            /* 6-connected neighbors */
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

/* Quick connection ratio: |P| at origin vs along axes */
static double connection_ratio(Grid *g) {
    const int N = g->N, NN = N*N;
    const double dx = g->dx, L = g->L;
    const double probe_dist = 3.0;  /* probe distance along each axis */

    /* Origin */
    int i0 = N/2, j0 = N/2, k0 = N/2;
    long idx0 = (long)i0*NN + j0*N + k0;
    double P0 = fabs(g->phi[0][idx0] * g->phi[1][idx0] * g->phi[2][idx0]);

    if (P0 < 1e-30) return 0.0;

    /* Probe along each axis at +/- probe_dist */
    double min_P = P0;
    int offsets[6][3] = {
        {1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1}
    };
    for (int p = 0; p < 6; p++) {
        int di = (int)(probe_dist / dx + 0.5) * offsets[p][0];
        int dj = (int)(probe_dist / dx + 0.5) * offsets[p][1];
        int dk = (int)(probe_dist / dx + 0.5) * offsets[p][2];
        int ii = i0 + di, jj = j0 + dj, kk = k0 + dk;
        if (ii < 0 || ii >= N || jj < 0 || jj >= N || kk < 0 || kk >= N) continue;
        long idx = (long)ii*NN + jj*N + kk;
        double P = fabs(g->phi[0][idx] * g->phi[1][idx] * g->phi[2][idx]);
        if (P < min_P) min_P = P;
    }

    return min_P / P0;
}

/* ================================================================
   Net z-winding number (from z-component phase gradient)
   ================================================================ */

static double compute_z_winding(Grid *g) {
    const int N = g->N, NN = N*N;

    /* Integrate winding of (phi_0, phi_1) around z-axis in xy-plane at z=0 */
    int k_mid = N / 2;
    double winding = 0;

    /* Sample on a circle of radius R in the z=0 plane */
    double R_sample = 2.0;
    int n_pts = 64;
    double L = g->L, dx = g->dx;

    for (int p = 0; p < n_pts; p++) {
        double theta1 = 2.0 * PI * p / n_pts;
        double theta2 = 2.0 * PI * (p + 1) / n_pts;

        double x1 = R_sample * cos(theta1), y1 = R_sample * sin(theta1);
        double x2 = R_sample * cos(theta2), y2 = R_sample * sin(theta2);

        /* Interpolate phi_0, phi_1 at both points (nearest grid) */
        int i1 = (int)((x1 + L) / dx + 0.5); if (i1<0) i1=0; if (i1>=N) i1=N-1;
        int j1 = (int)((y1 + L) / dx + 0.5); if (j1<0) j1=0; if (j1>=N) j1=N-1;
        int i2 = (int)((x2 + L) / dx + 0.5); if (i2<0) i2=0; if (i2>=N) i2=N-1;
        int j2 = (int)((y2 + L) / dx + 0.5); if (j2<0) j2=0; if (j2>=N) j2=N-1;

        long idx_1 = (long)i1 * NN + j1 * N + k_mid;
        long idx_2 = (long)i2 * NN + j2 * N + k_mid;

        double p0_1 = g->phi[0][idx_1], p1_1 = g->phi[1][idx_1];
        double p0_2 = g->phi[0][idx_2], p1_2 = g->phi[1][idx_2];

        /* Winding = sum of angle differences */
        double ang1 = atan2(p1_1, p0_1);
        double ang2 = atan2(p1_2, p0_2);
        double dang = ang2 - ang1;
        while (dang > PI) dang -= 2*PI;
        while (dang < -PI) dang += 2*PI;
        winding += dang;
    }

    return winding / (2.0 * PI);
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
   Main
   ================================================================ */

int main(int argc, char **argv) {
    int N = 128;
    double L = 15.0, T = 500.0;
    double diag_dt = 1.0, snap_dt = 5.0;
    char outdir[256] = "data/seed_test";
    char seedpath[512] = "";

    for (int i = 1; i < argc; i++) {
        if      (!strcmp(argv[i],"-seed"))   strncpy(seedpath, argv[++i], 511);
        else if (!strcmp(argv[i],"-N"))      N = atoi(argv[++i]);
        else if (!strcmp(argv[i],"-L"))      L = atof(argv[++i]);
        else if (!strcmp(argv[i],"-T"))      T = atof(argv[++i]);
        else if (!strcmp(argv[i],"-m"))      { double m = atof(argv[++i]); MASS2 = m*m; }
        else if (!strcmp(argv[i],"-mt"))     { double m = atof(argv[++i]); MTHETA2 = m*m; }
        else if (!strcmp(argv[i],"-eta"))    ETA = atof(argv[++i]);
        else if (!strcmp(argv[i],"-mu"))     MU = atof(argv[++i]);
        else if (!strcmp(argv[i],"-kappa"))  KAPPA = atof(argv[++i]);
        else if (!strcmp(argv[i],"-diag"))   diag_dt = atof(argv[++i]);
        else if (!strcmp(argv[i],"-snap"))   snap_dt = atof(argv[++i]);
        else if (!strcmp(argv[i],"-o"))      strncpy(outdir, argv[++i], 255);
        else if (!strcmp(argv[i],"-damp_width")) DAMP_WIDTH = atof(argv[++i]);
        else if (!strcmp(argv[i],"-damp_rate"))  DAMP_RATE = atof(argv[++i]);
        else {
            fprintf(stderr, "Unknown option: %s\n", argv[i]);
            fprintf(stderr, "Usage: %s -seed init.bin [-N 128] [-L 15] [-T 500] [-snap 5] "
                    "[-damp_width 4] [-damp_rate 0.005] [-o data/seed_test]\n", argv[0]);
            return 1;
        }
    }

    if (seedpath[0] == '\0') {
        fprintf(stderr, "ERROR: must specify -seed <path>\n");
        fprintf(stderr, "Usage: %s -seed init.bin [-N 128] [-L 15] [-T 500] [-snap 5] "
                "[-damp_width 4] [-damp_rate 0.005] [-o data/seed_test]\n", argv[0]);
        return 1;
    }

    int nthreads = 4;
    char *env_threads = getenv("OMP_NUM_THREADS");
    if (env_threads) nthreads = atoi(env_threads);
    omp_set_num_threads(nthreads);

    printf("=== V37 Seed Run ===\n");
    printf("Equations: 6-field Cosserat (3 phi + 3 theta)\n");
    printf("m^2=%.4f, m_theta^2=%.4f, eta=%.3f, mu=%.3f, kappa=%.1f\n",
           MASS2, MTHETA2, ETA, MU, KAPPA);
    printf("Absorbing BC: damp_width=%.1f damp_rate=%.4f\n", DAMP_WIDTH, DAMP_RATE);
    printf("N=%d L=%.1f T=%.0f\n", N, L, T);

    mkdir("data", 0755); mkdir(outdir, 0755);
    Grid *g = grid_alloc(N, L);
    printf("dx=%.4f dt=%.5f threads=%d\n\n", g->dx, g->dt, nthreads);

    /* Load seed file as initial condition */
    load_seed(g, seedpath);

    compute_forces(g);

    /* SFA archive with f32 columns for smaller files */
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

    /* Write t=0 frame (f32) */
    {
        float *p0 = to_f32(g->phi[0], g->N3);
        float *c0 = malloc(g->N3 * sizeof(float)); memcpy(c0, p0, g->N3*sizeof(float));
        float *p1 = to_f32(g->phi[1], g->N3);
        float *c1 = malloc(g->N3 * sizeof(float)); memcpy(c1, p1, g->N3*sizeof(float));
        float *p2 = to_f32(g->phi[2], g->N3);
        float *c2 = malloc(g->N3 * sizeof(float)); memcpy(c2, p2, g->N3*sizeof(float));
        float *t0 = to_f32(g->theta[0], g->N3);
        float *ct0 = malloc(g->N3 * sizeof(float)); memcpy(ct0, t0, g->N3*sizeof(float));
        float *t1 = to_f32(g->theta[1], g->N3);
        float *ct1 = malloc(g->N3 * sizeof(float)); memcpy(ct1, t1, g->N3*sizeof(float));
        float *t2 = to_f32(g->theta[2], g->N3);
        float *ct2 = malloc(g->N3 * sizeof(float)); memcpy(ct2, t2, g->N3*sizeof(float));
        /* Debug: verify f32 data before SFA write */
        float fmax0 = 0;
        for (long i = 0; i < g->N3; i++) if (fabsf(c0[i]) > fmax0) fmax0 = fabsf(c0[i]);
               fmax0, c0[0], g->phi[0][0]);
        void *cols[] = {c0, c1, c2, ct0, ct1, ct2};
        sfa_write_frame(sfa, 0, cols);
        free(c0); free(c1); free(c2); free(ct0); free(ct1); free(ct2);
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
        printf("Initial |P| integrated: %.4e\n", P_int);
        double winding = compute_z_winding(g);
        printf("Initial z-winding: %.4f\n", winding);
        double conn = connection_ratio(g);
        printf("Initial connection ratio: %.4f\n\n", conn);
    }

    /* Timeseries file */
    char tspath[512];
    snprintf(tspath, sizeof(tspath), "%s/timeseries.tsv", outdir);
    FILE *fp = fopen(tspath, "w");
    fprintf(fp, "t\tE_phi_kin\tE_theta_kin\tE_grad\tE_mass\tE_pot\tE_tgrad\tE_tmass\tE_coupling\tE_total\ttheta_rms\tP_int\taspect\tz_winding\tconn_ratio\tn_clusters\n");

    int n_steps = (int)(T / g->dt);
    int diag_every = (int)(diag_dt / g->dt); if (diag_every<1) diag_every=1;
    int snap_every = (int)(snap_dt / g->dt); if (snap_every<1) snap_every=1;

    printf("Steps=%d diag_every=%d snap_every=%d (snap_dt=%.1f)\n\n", n_steps, diag_every, snap_every, snap_dt);

    double wall0 = omp_get_wtime();
    double E0 = 0, P0 = 0, Ep0 = 0;

    /* Time-averaged death check: need time corresponding to ~50 time units */
    double death_time_window = 50.0;
    int death_min_samples = (int)(death_time_window / diag_dt);
    if (death_min_samples < 20) death_min_samples = 20;
    if (death_min_samples > DEATH_WINDOW) death_min_samples = DEATH_WINDOW;

    int fragmented = 0;

    for (int step = 0; step <= n_steps; step++) {
        if (step > 0) verlet_step(g);
        double t = step * g->dt;

        /* SFA snapshot (f32) */
        if (step > 0 && step % snap_every == 0) {
            float *p0 = to_f32(g->phi[0], g->N3);
            float *c0 = malloc(g->N3 * sizeof(float)); memcpy(c0, p0, g->N3*sizeof(float));
            float *p1 = to_f32(g->phi[1], g->N3);
            float *c1 = malloc(g->N3 * sizeof(float)); memcpy(c1, p1, g->N3*sizeof(float));
            float *p2 = to_f32(g->phi[2], g->N3);
            float *c2 = malloc(g->N3 * sizeof(float)); memcpy(c2, p2, g->N3*sizeof(float));
            float *t0 = to_f32(g->theta[0], g->N3);
            float *ct0 = malloc(g->N3 * sizeof(float)); memcpy(ct0, t0, g->N3*sizeof(float));
            float *t1 = to_f32(g->theta[1], g->N3);
            float *ct1 = malloc(g->N3 * sizeof(float)); memcpy(ct1, t1, g->N3*sizeof(float));
            float *t2 = to_f32(g->theta[2], g->N3);
            float *ct2 = malloc(g->N3 * sizeof(float)); memcpy(ct2, t2, g->N3*sizeof(float));
            void *cols[] = {c0, c1, c2, ct0, ct1, ct2};
            sfa_write_frame(sfa, t, cols);
            free(c0); free(c1); free(c2); free(ct0); free(ct1); free(ct2);
        }

        if (step % diag_every == 0) {
            double epk, etk, eg, em, ep, etg, etm, ec, et;
            compute_energy(g, &epk, &etk, &eg, &em, &ep, &etg, &etm, &ec, &et);
            double trms = theta_rms(g);
            double P_int = P_integrated(g);

            /* Inertia tensor + fragmentation at major intervals */
            int is_major = (step % (diag_every*50) == 0);
            static double last_aspect = 1.0;
            static double last_winding = 0.0;
            static double last_conn = 1.0;
            static int last_clusters = 1;
            double aspect = last_aspect;
            double winding = last_winding;
            double conn = last_conn;
            int n_clusters = last_clusters;

            if (is_major) {
                double Ixx[3][3], evals[3], E_int;
                compute_inertia_tensor(g, Ixx, &E_int);
                eigen3x3(Ixx, evals);
                aspect = (evals[2] > 1e-30) ? evals[0]/evals[2] : 999;
                last_aspect = aspect;
                winding = compute_z_winding(g);
                last_winding = winding;

                /* Fragmentation detection */
                conn = connection_ratio(g);
                last_conn = conn;

                /* Full cluster detection: use stride=4 for speed */
                double P_thresh = 0.001;  /* threshold for |P| */
                n_clusters = detect_clusters(g, P_thresh, 4);
                last_clusters = n_clusters;

                if (conn < 0.1 && !fragmented) {
                    fragmented = 1;
                    printf("  *** FRAGMENTATION DETECTED at t=%.1f (conn_ratio=%.4f, clusters=%d) ***\n",
                           t, conn, n_clusters);
                }
            }

            if (step == 0) { E0 = et; P0 = P_int; Ep0 = ep; }

            /* Record for time-averaged death check */
            ep_record(ep);

            fprintf(fp, "%.2f\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%.4f\t%.4f\t%.4f\t%d\n",
                    t, epk, etk, eg, em, ep, etg, etm, ec, et, trms, P_int, aspect, winding, conn, n_clusters);
            fflush(fp);

            if (is_major) {
                double wall = omp_get_wtime() - wall0;
                double drift = 100.0*(et-E0)/(fabs(E0)+1e-30);
                double ep_avg = ep_average();
                printf("t=%7.1f E=%.3e (drift %+.3f%%) Ep=%.1f Ep_avg=%.1f P_int=%.3e trms=%.3e asp=%.2f wind=%.3f conn=%.3f clust=%d [%.0f%% %.0fs]\n",
                       t, et, drift, ep, ep_avg, P_int, trms, aspect, winding, conn, n_clusters, 100.0*step/n_steps, wall);
                fflush(stdout);
            }

            /* TIME-AVERAGED death check (fixes breathing-mode false kills) */
            if (ep_count >= death_min_samples && fabs(Ep0) > 1e-10) {
                double ep_avg = ep_average();
                if (ep_avg < 0.01 * fabs(Ep0)) {
                    printf("\n  *** Time-averaged |E_pot| (%.3e) < 1%% of initial (%.3e) -- structure DEAD. Stopping. ***\n",
                           ep_avg, fabs(Ep0));
                    /* Write final frame */
                    float *p0 = to_f32(g->phi[0], g->N3);
                    float *c0 = malloc(g->N3 * sizeof(float)); memcpy(c0, p0, g->N3*sizeof(float));
                    float *p1 = to_f32(g->phi[1], g->N3);
                    float *c1 = malloc(g->N3 * sizeof(float)); memcpy(c1, p1, g->N3*sizeof(float));
                    float *p2 = to_f32(g->phi[2], g->N3);
                    float *c2 = malloc(g->N3 * sizeof(float)); memcpy(c2, p2, g->N3*sizeof(float));
                    float *t0f = to_f32(g->theta[0], g->N3);
                    float *ct0 = malloc(g->N3 * sizeof(float)); memcpy(ct0, t0f, g->N3*sizeof(float));
                    float *t1f = to_f32(g->theta[1], g->N3);
                    float *ct1 = malloc(g->N3 * sizeof(float)); memcpy(ct1, t1f, g->N3*sizeof(float));
                    float *t2f = to_f32(g->theta[2], g->N3);
                    float *ct2 = malloc(g->N3 * sizeof(float)); memcpy(ct2, t2f, g->N3*sizeof(float));
                    void *cols[] = {c0, c1, c2, ct0, ct1, ct2};
                    sfa_write_frame(sfa, t, cols);
                    free(c0); free(c1); free(c2); free(ct0); free(ct1); free(ct2);
                    break;
                }
            }
        }
    }

    /* Final frame */
    {
        double t_final = n_steps * g->dt;
        float *p0 = to_f32(g->phi[0], g->N3);
        float *c0 = malloc(g->N3 * sizeof(float)); memcpy(c0, p0, g->N3*sizeof(float));
        float *p1 = to_f32(g->phi[1], g->N3);
        float *c1 = malloc(g->N3 * sizeof(float)); memcpy(c1, p1, g->N3*sizeof(float));
        float *p2 = to_f32(g->phi[2], g->N3);
        float *c2 = malloc(g->N3 * sizeof(float)); memcpy(c2, p2, g->N3*sizeof(float));
        float *t0f = to_f32(g->theta[0], g->N3);
        float *ct0 = malloc(g->N3 * sizeof(float)); memcpy(ct0, t0f, g->N3*sizeof(float));
        float *t1f = to_f32(g->theta[1], g->N3);
        float *ct1 = malloc(g->N3 * sizeof(float)); memcpy(ct1, t1f, g->N3*sizeof(float));
        float *t2f = to_f32(g->theta[2], g->N3);
        float *ct2 = malloc(g->N3 * sizeof(float)); memcpy(ct2, t2f, g->N3*sizeof(float));
        void *cols[] = {c0, c1, c2, ct0, ct1, ct2};
        sfa_write_frame(sfa, t_final, cols);
        free(c0); free(c1); free(c2); free(ct0); free(ct1); free(ct2);
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
        double winding = compute_z_winding(g);
        double conn = connection_ratio(g);
        int n_clusters = detect_clusters(g, 0.001, 4);

        printf("\n=== FINAL SUMMARY: seed run ===\n");
        printf("Seed: %s\n", seedpath);
        printf("E_total=%.4e (drift %.3f%%)\n", et, 100*(et-E0)/(fabs(E0)+1e-30));
        printf("E_pot=%.4e (%.1f%% of initial)\n", ep, 100*ep/(Ep0+1e-30));
        printf("|P| integrated=%.4e (%.1f%% of initial)\n", P_int, 100*P_int/(P0+1e-30));
        printf("theta_rms=%.4e\n", trms);
        printf("z_winding=%.4f\n", winding);
        printf("Aspect ratio=%.3f (eigenvalues: %.1f %.1f %.1f)\n",
               aspect, evals[0], evals[1], evals[2]);
        printf("Connection ratio=%.4f  Clusters=%d  Fragmented=%s\n",
               conn, n_clusters, fragmented ? "YES" : "NO");
        printf("Survival: %s\n", (fabs(ep) > 0.1*fabs(Ep0)) ? "YES (instantaneous)" : "NO (instantaneous)");
        double ep_avg = ep_average();
        printf("Survival (avg): %s (Ep_avg=%.3e vs Ep0=%.3e)\n",
               (ep_avg > 0.1*fabs(Ep0)) ? "YES" : "NO", ep_avg, Ep0);
    }

    fclose(fp);

    if (f32_buf) free(f32_buf);

    double wall = omp_get_wtime() - wall0;
    printf("\n=== Complete: %.0fs (%.1f min) ===\n", wall, wall/60);

    grid_free(g);
    return 0;
}
