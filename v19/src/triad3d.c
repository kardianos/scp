/*
 * v19/src/triad3d.c
 *
 * 3D Skyrme sigma-model simulator for three-mode resonance testing.
 * Adapted from v18/src/skyrme3d.c with new B=0 initializations,
 * per-mode energy tracking, mode decomposition, time-averaged shell
 * accumulator, core energy fraction, and DFT spectrum.
 *
 * Physics:
 *   Lagrangian L = (1/2)|dq|^2 + (c4/4)[Tr(D)^2 - Tr(D^2)]
 *   on S^3 (unit quaternion field q, |q|=1).
 *   c4 = 2 rho0^2 / e^2.  Code units: e=1, rho0=1 => c4=2.
 *
 *   EOM: q_tt = Lap(q) + c4 * div(G) - (q.F + |v|^2) q
 *   where G_i = Tr(D) dq_i - D_{ij} dq_j
 *
 * Tests:
 *   1: B=0 hedgehog bump (Experiment 1A)
 *   2: Three-mode kicked triad (Experiment 1B)
 *   3: Three-mode displaced triad (Experiment 1C)
 *   4: Single-mode control (Experiment 1D)
 *   5: No-Skyrme bump (Experiment 1E)
 *   6: No-Skyrme triad (Experiment 1F)
 *
 * Compile: gcc -O3 -fopenmp -o triad3d v19/src/triad3d.c -lm
 * Run:     ./triad3d -test 2 -N 200 -L 16 -steps 6000
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>

/* --- Simulation state --- */
static int    N;            /* grid points per dimension */
static double L, dx, dt;   /* box size, spacing, timestep */
static double c4;           /* Skyrme coupling = 2*rho0^2/e^2 */
static int    ntot;          /* N^3 */

static double *R[4], *V[4]; /* quaternion field & velocity */
static double *G[3][4];     /* Skyrme stress tensor G_d^a */
static double *acc[4];      /* acceleration buffer */
static double *damp;        /* absorbing boundary mask */

/* --- Parameters for mode initialization --- */
static double param_A     = 0.8;
static double param_sigma = 1.5;
static double param_R0    = 3.0;
static double param_omega1 = 0.0;  /* 0 = use default pi*c/R0 */
static double param_omega2 = 0.0;
static double param_omega3 = 0.0;
static double param_A1    = 0.4;
static double param_A2    = 0.4;
static double param_A3    = 0.56;

/* --- Periodic index --- */
static inline int pbc(int i) { return ((i % N) + N) % N; }
static inline int flat(int i, int j, int k) {
    return pbc(i)*N*N + pbc(j)*N + pbc(k);
}
static inline double coord(int i) { return -L/2 + i*dx; }

/* --- sinc function: sin(x)/x, sinc(0) = 1 --- */
static inline double sinc(double x) {
    if (fabs(x) < 1e-12) return 1.0;
    return sin(x) / x;
}

/* --- Memory --- */
static void alloc_arrays(void) {
    ntot = N*N*N;
    for (int c = 0; c < 4; c++) {
        R[c]   = calloc(ntot, sizeof(double));
        V[c]   = calloc(ntot, sizeof(double));
        acc[c] = calloc(ntot, sizeof(double));
        for (int d = 0; d < 3; d++)
            G[d][c] = calloc(ntot, sizeof(double));
    }
    damp = malloc(ntot * sizeof(double));
}

static void free_arrays(void) {
    for (int c = 0; c < 4; c++) {
        free(R[c]); free(V[c]); free(acc[c]);
        for (int d = 0; d < 3; d++) free(G[d][c]);
    }
    free(damp);
}

/* --- Absorbing boundary layer ---
 * damp=1 in interior, damp -> 0.05 at edges.
 * Applied to velocity each step to absorb outgoing radiation. */
static void build_damping(double frac) {
    int b = (int)(N * frac);
    if (b < 2) b = 2;
    #pragma omp parallel for collapse(3)
    for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
    for (int k = 0; k < N; k++) {
        double m = 1.0;
        int ii[3] = {i, j, k};
        for (int d = 0; d < 3; d++) {
            if (ii[d] < b)
                m *= (double)ii[d] / b;
            else if (ii[d] >= N - b)
                m *= (double)(N - 1 - ii[d]) / b;
        }
        damp[i*N*N + j*N + k] = 1.0 - 0.95*(1.0 - m);
    }
}

/* --- Radial envelope g_n(r) = sinc(n*pi*r/R0) * exp(-r^2/(2*sigma^2)) --- */
static double g_envelope(int n, double r, double R0, double sigma) {
    double arg = n * M_PI * r / R0;
    return sinc(arg) * exp(-r*r / (2.0*sigma*sigma));
}

/* ====================== INITIALIZATIONS ====================== */

/* Test 1: B=0 Hedgehog Bump (Experiment 1A) */
static void init_bump(double A, double sigma) {
    #pragma omp parallel for collapse(3)
    for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
    for (int k = 0; k < N; k++) {
        double x = coord(i), y = coord(j), z = coord(k);
        double r = sqrt(x*x + y*y + z*z);
        int p = i*N*N + j*N + k;

        double env = A * exp(-r*r / (2.0*sigma*sigma));
        double q1, q2, q3;

        if (r < 1e-12) {
            q1 = q2 = q3 = env / sqrt(3.0);
        } else {
            q1 = env * x / r;
            q2 = env * y / r;
            q3 = env * z / r;
        }

        double sum2 = q1*q1 + q2*q2 + q3*q3;
        if (sum2 > 0.999) {
            double scale = sqrt(0.999 / sum2);
            q1 *= scale; q2 *= scale; q3 *= scale;
            sum2 = 0.999;
        }

        R[0][p] = sqrt(1.0 - sum2);
        R[1][p] = q1;
        R[2][p] = q2;
        R[3][p] = q3;
        V[0][p] = V[1][p] = V[2][p] = V[3][p] = 0.0;
    }
}

/* Test 2: Three-Mode Kicked Triad (Experiment 1B) */
static void init_kicked_triad(double A, double sigma, double R0,
                              double w1, double w2, double w3) {
    #pragma omp parallel for collapse(3)
    for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
    for (int k = 0; k < N; k++) {
        double x = coord(i), y = coord(j), z = coord(k);
        double r = sqrt(x*x + y*y + z*z);
        int p = i*N*N + j*N + k;

        /* Position: vacuum */
        R[0][p] = 1.0;
        R[1][p] = R[2][p] = R[3][p] = 0.0;

        /* Velocity: three orthogonal kicks */
        double g1 = g_envelope(1, r, R0, sigma);
        double g2 = g_envelope(2, r, R0, sigma);
        double g3 = g_envelope(3, r, R0, sigma);

        if (r < 1e-12) {
            V[1][p] = w1 * A * g1;
            V[2][p] = w2 * A * g2;
            V[3][p] = w3 * A * g3;
        } else {
            V[1][p] = w1 * A * g1 * x / r;
            V[2][p] = w2 * A * g2 * y / r;
            V[3][p] = w3 * A * g3 * z / r;
        }
        V[0][p] = 0.0;  /* tangent since q_a=0 */
    }
}

/* Test 3: Three-Mode Displaced Triad (Experiment 1C) */
static void init_displaced_triad(double A1, double A2, double A3,
                                 double sigma, double R0) {
    #pragma omp parallel for collapse(3)
    for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
    for (int k = 0; k < N; k++) {
        double x = coord(i), y = coord(j), z = coord(k);
        double r = sqrt(x*x + y*y + z*z);
        int p = i*N*N + j*N + k;

        double g1 = g_envelope(1, r, R0, sigma);
        double g2 = g_envelope(2, r, R0, sigma);
        double g3 = g_envelope(3, r, R0, sigma);

        double q1, q2, q3;
        if (r < 1e-12) {
            q1 = A1 * g1;
            q2 = A2 * g2;
            q3 = A3 * g3;
        } else {
            q1 = A1 * g1 * x / r;
            q2 = A2 * g2 * y / r;
            q3 = A3 * g3 * z / r;
        }

        double sum2 = q1*q1 + q2*q2 + q3*q3;
        if (sum2 > 0.999) {
            double scale = sqrt(0.999 / sum2);
            q1 *= scale; q2 *= scale; q3 *= scale;
            sum2 = 0.999;
        }

        R[0][p] = sqrt(1.0 - sum2);
        R[1][p] = q1;
        R[2][p] = q2;
        R[3][p] = q3;
        V[0][p] = V[1][p] = V[2][p] = V[3][p] = 0.0;
    }
}

/* Test 4: Single-Mode Control (Experiment 1D) */
static void init_single_mode(double A, double sigma, double R0) {
    #pragma omp parallel for collapse(3)
    for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
    for (int k = 0; k < N; k++) {
        double x = coord(i), y = coord(j), z = coord(k);
        double r = sqrt(x*x + y*y + z*z);
        int p = i*N*N + j*N + k;

        double g1 = g_envelope(1, r, R0, sigma);
        double q1;
        if (r < 1e-12) {
            q1 = A * g1;
        } else {
            q1 = A * g1 * x / r;
        }

        if (q1*q1 > 0.999) q1 = (q1 > 0 ? 1.0 : -1.0) * sqrt(0.999);

        R[0][p] = sqrt(1.0 - q1*q1);
        R[1][p] = q1;
        R[2][p] = 0.0;
        R[3][p] = 0.0;
        V[0][p] = V[1][p] = V[2][p] = V[3][p] = 0.0;
    }
}

/* ====================== CORE PHYSICS (from v18) ====================== */

/* --- Project R onto S^3, V onto tangent space --- */
static void project(void) {
    #pragma omp parallel for
    for (int p = 0; p < ntot; p++) {
        double norm = sqrt(R[0][p]*R[0][p] + R[1][p]*R[1][p] +
                          R[2][p]*R[2][p] + R[3][p]*R[3][p]);
        for (int c = 0; c < 4; c++) R[c][p] /= norm;
        double dot = 0;
        for (int c = 0; c < 4; c++) dot += R[c][p] * V[c][p];
        for (int c = 0; c < 4; c++) V[c][p] -= dot * R[c][p];
    }
}

/* --- Compute Skyrme stress G[d][c] at all grid points ---
 * D_{ij} = d_i q . d_j q     (strain tensor)
 * G_i^a = Tr(D) d_i q^a  -  sum_j D_{ij} d_j q^a */
static void compute_G(void) {
    double idx2 = 1.0 / (2.0 * dx);

    #pragma omp parallel for collapse(3)
    for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
    for (int k = 0; k < N; k++) {
        int p  = i*N*N + j*N + k;
        int nb[3][2];
        nb[0][0] = flat(i+1,j,k); nb[0][1] = flat(i-1,j,k);
        nb[1][0] = flat(i,j+1,k); nb[1][1] = flat(i,j-1,k);
        nb[2][0] = flat(i,j,k+1); nb[2][1] = flat(i,j,k-1);

        double dq[3][4];
        for (int d = 0; d < 3; d++)
            for (int c = 0; c < 4; c++)
                dq[d][c] = (R[c][nb[d][0]] - R[c][nb[d][1]]) * idx2;

        double D[3][3];
        for (int a = 0; a < 3; a++)
        for (int b = a; b < 3; b++) {
            double s = 0;
            for (int c = 0; c < 4; c++) s += dq[a][c] * dq[b][c];
            D[a][b] = s;
            if (b > a) D[b][a] = s;
        }
        double trD = D[0][0] + D[1][1] + D[2][2];

        for (int d = 0; d < 3; d++)
        for (int c = 0; c < 4; c++) {
            double Ddq = 0;
            for (int e = 0; e < 3; e++) Ddq += D[d][e] * dq[e][c];
            G[d][c][p] = trD * dq[d][c] - Ddq;
        }
    }
}

/* --- Compute acceleration ---
 * a = Lap(q) + c4*div(G) - (q.F + |v|^2) q */
static void compute_accel(void) {
    double idx2  = 1.0 / (2.0 * dx);
    double idxsq = 1.0 / (dx * dx);

    compute_G();

    #pragma omp parallel for collapse(3)
    for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
    for (int k = 0; k < N; k++) {
        int p = i*N*N + j*N + k;
        int nb[3][2];
        nb[0][0] = flat(i+1,j,k); nb[0][1] = flat(i-1,j,k);
        nb[1][0] = flat(i,j+1,k); nb[1][1] = flat(i,j-1,k);
        nb[2][0] = flat(i,j,k+1); nb[2][1] = flat(i,j,k-1);

        double F[4];
        for (int c = 0; c < 4; c++) {
            double lap = (R[c][nb[0][0]] + R[c][nb[0][1]] +
                         R[c][nb[1][0]] + R[c][nb[1][1]] +
                         R[c][nb[2][0]] + R[c][nb[2][1]] - 6*R[c][p]) * idxsq;

            double divG = 0;
            for (int d = 0; d < 3; d++)
                divG += (G[d][c][nb[d][0]] - G[d][c][nb[d][1]]) * idx2;

            F[c] = lap + c4 * divG;
        }

        double qdotF = 0, v2 = 0;
        for (int c = 0; c < 4; c++) {
            qdotF += R[c][p] * F[c];
            v2    += V[c][p] * V[c][p];
        }
        for (int c = 0; c < 4; c++)
            acc[c][p] = F[c] - (qdotF + v2) * R[c][p];
    }
}

/* --- Velocity-Verlet symplectic step --- */
static void step_vv(void) {
    /* Half-kick */
    compute_accel();
    #pragma omp parallel for
    for (int p = 0; p < ntot; p++)
        for (int c = 0; c < 4; c++)
            V[c][p] += 0.5 * dt * acc[c][p];

    /* Drift */
    #pragma omp parallel for
    for (int p = 0; p < ntot; p++)
        for (int c = 0; c < 4; c++)
            R[c][p] += dt * V[c][p];
    project();

    /* Half-kick at new position */
    compute_accel();
    #pragma omp parallel for
    for (int p = 0; p < ntot; p++)
        for (int c = 0; c < 4; c++)
            V[c][p] += 0.5 * dt * acc[c][p];

    /* Absorbing boundary */
    #pragma omp parallel for
    for (int p = 0; p < ntot; p++) {
        double d = damp[p];
        for (int c = 0; c < 4; c++) V[c][p] *= d;
    }
}

/* ====================== DIAGNOSTICS ====================== */

/* --- Energy diagnostics ---
 * E_kin = (1/2) int |V|^2
 * E_2   = (1/2) int Tr(D)
 * E_4   = (c4/4) int [Tr(D)^2 - Tr(D^2)] */
static void compute_energy(double *Ekin, double *E2, double *E4) {
    double ek = 0, e2 = 0, e4 = 0;
    double idx2 = 1.0 / (2.0 * dx);
    double h3 = dx*dx*dx;

    #pragma omp parallel for collapse(3) reduction(+:ek,e2,e4)
    for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
    for (int k = 0; k < N; k++) {
        int p = i*N*N + j*N + k;
        int nb[3][2];
        nb[0][0] = flat(i+1,j,k); nb[0][1] = flat(i-1,j,k);
        nb[1][0] = flat(i,j+1,k); nb[1][1] = flat(i,j-1,k);
        nb[2][0] = flat(i,j,k+1); nb[2][1] = flat(i,j,k-1);

        double v2 = 0;
        for (int c = 0; c < 4; c++) v2 += V[c][p] * V[c][p];
        ek += 0.5 * v2;

        double dq[3][4];
        for (int d = 0; d < 3; d++)
            for (int c = 0; c < 4; c++)
                dq[d][c] = (R[c][nb[d][0]] - R[c][nb[d][1]]) * idx2;

        double D[3][3];
        for (int a = 0; a < 3; a++)
        for (int b = a; b < 3; b++) {
            double s = 0;
            for (int c = 0; c < 4; c++) s += dq[a][c] * dq[b][c];
            D[a][b] = s;
            if (b > a) D[b][a] = s;
        }

        double trD  = D[0][0] + D[1][1] + D[2][2];
        double trD2 = 0;
        for (int a = 0; a < 3; a++)
            for (int b = 0; b < 3; b++)
                trD2 += D[a][b] * D[a][b];

        e2 += 0.5 * trD;
        e4 += 0.25 * c4 * (trD*trD - trD2);
    }

    *Ekin = ek * h3;
    *E2   = e2 * h3;
    *E4   = e4 * h3;
}

/* --- Baryon number (topological charge) --- */
static double compute_baryon(void) {
    double Q = 0;
    double idx2 = 1.0 / (2.0 * dx);
    double h3 = dx*dx*dx;

    #pragma omp parallel for collapse(3) reduction(+:Q)
    for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
    for (int k = 0; k < N; k++) {
        int p = i*N*N + j*N + k;
        int nb[3][2];
        nb[0][0] = flat(i+1,j,k); nb[0][1] = flat(i-1,j,k);
        nb[1][0] = flat(i,j+1,k); nb[1][1] = flat(i,j-1,k);
        nb[2][0] = flat(i,j,k+1); nb[2][1] = flat(i,j,k-1);

        double q[4], dq[3][4];
        for (int c = 0; c < 4; c++) {
            q[c] = R[c][p];
            for (int d = 0; d < 3; d++)
                dq[d][c] = (R[c][nb[d][0]] - R[c][nb[d][1]]) * idx2;
        }

        double l[3][3];
        for (int d = 0; d < 3; d++) {
            double *v = dq[d];
            l[d][0] = q[0]*v[1] - q[1]*v[0] - q[2]*v[3] + q[3]*v[2];
            l[d][1] = q[0]*v[2] + q[1]*v[3] - q[2]*v[0] - q[3]*v[1];
            l[d][2] = q[0]*v[3] - q[1]*v[2] + q[2]*v[1] - q[3]*v[0];
        }

        double det = l[0][0]*(l[1][1]*l[2][2] - l[1][2]*l[2][1])
                   - l[0][1]*(l[1][0]*l[2][2] - l[1][2]*l[2][0])
                   + l[0][2]*(l[1][0]*l[2][1] - l[1][1]*l[2][0]);
        Q += det;
    }

    return -Q * h3 / (2.0 * M_PI * M_PI);
}

/* --- Per-mode energy tracking (approximate, ignores cross-terms) ---
 * E_kin_a = (h^3/2) sum v_a^2
 * E_2_a   = (h^3/2) sum |grad q_a|^2
 * E_a     = E_kin_a + E_2_a */
static void compute_mode_energy(double E_mode[3]) {
    double ek[3] = {0,0,0}, e2[3] = {0,0,0};
    double idx2 = 1.0 / (2.0 * dx);
    double h3 = dx*dx*dx;

    #pragma omp parallel for collapse(3) \
        reduction(+:ek[:3]) reduction(+:e2[:3])
    for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
    for (int k = 0; k < N; k++) {
        int p = i*N*N + j*N + k;
        int nb[3][2];
        nb[0][0] = flat(i+1,j,k); nb[0][1] = flat(i-1,j,k);
        nb[1][0] = flat(i,j+1,k); nb[1][1] = flat(i,j-1,k);
        nb[2][0] = flat(i,j,k+1); nb[2][1] = flat(i,j,k-1);

        /* Per-component kinetic energy */
        for (int a = 0; a < 3; a++)
            ek[a] += 0.5 * V[a+1][p] * V[a+1][p];

        /* Per-component gradient energy */
        for (int a = 0; a < 3; a++) {
            double g2 = 0;
            for (int d = 0; d < 3; d++) {
                double gr = (R[a+1][nb[d][0]] - R[a+1][nb[d][1]]) * idx2;
                g2 += gr * gr;
            }
            e2[a] += 0.5 * g2;
        }
    }

    for (int a = 0; a < 3; a++)
        E_mode[a] = (ek[a] + e2[a]) * h3;
}

/* --- Core energy fraction: energy within r < R_core / E_total --- */
static double compute_core_fraction(double R_core) {
    double e_core = 0, e_total = 0;
    double idx2 = 1.0 / (2.0 * dx);
    double R2 = R_core * R_core;

    #pragma omp parallel for collapse(3) reduction(+:e_core,e_total)
    for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
    for (int k = 0; k < N; k++) {
        int p = i*N*N + j*N + k;
        int nb[3][2];
        nb[0][0] = flat(i+1,j,k); nb[0][1] = flat(i-1,j,k);
        nb[1][0] = flat(i,j+1,k); nb[1][1] = flat(i,j-1,k);
        nb[2][0] = flat(i,j,k+1); nb[2][1] = flat(i,j,k-1);

        double v2 = 0;
        for (int c = 0; c < 4; c++) v2 += V[c][p] * V[c][p];

        double dq[3][4];
        for (int d = 0; d < 3; d++)
            for (int c = 0; c < 4; c++)
                dq[d][c] = (R[c][nb[d][0]] - R[c][nb[d][1]]) * idx2;

        double D[3][3];
        for (int a = 0; a < 3; a++)
        for (int b = a; b < 3; b++) {
            double s = 0;
            for (int c = 0; c < 4; c++) s += dq[a][c] * dq[b][c];
            D[a][b] = s;
            if (b > a) D[b][a] = s;
        }
        double trD  = D[0][0] + D[1][1] + D[2][2];
        double trD2 = 0;
        for (int a = 0; a < 3; a++)
            for (int b = 0; b < 3; b++)
                trD2 += D[a][b] * D[a][b];

        double ed = 0.5*v2 + 0.5*trD + 0.25*c4*(trD*trD - trD2);
        e_total += ed;

        double x = coord(i), y = coord(j), z = coord(k);
        if (x*x + y*y + z*z < R2)
            e_core += ed;
    }

    /* No need to multiply by h3 since ratio cancels */
    return (e_total > 1e-30) ? e_core / e_total : 0.0;
}

/* --- Mode decomposition: project field onto initial mode shapes --- */
static void compute_mode_projection(double c_out[3], double sigma, double R0) {
    double c1 = 0, c2 = 0, c3 = 0;
    double h3 = dx*dx*dx;

    #pragma omp parallel for collapse(3) reduction(+:c1,c2,c3)
    for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
    for (int k = 0; k < N; k++) {
        double x = coord(i), y = coord(j), z = coord(k);
        double r = sqrt(x*x + y*y + z*z);
        int p = i*N*N + j*N + k;

        double g1 = g_envelope(1, r, R0, sigma);
        double g2 = g_envelope(2, r, R0, sigma);
        double g3 = g_envelope(3, r, R0, sigma);

        if (r < 1e-12) {
            /* At origin, x/r -> unit vector; use average projection */
            c1 += R[1][p] * g1;
            c2 += R[2][p] * g2;
            c3 += R[3][p] * g3;
        } else {
            c1 += R[1][p] * g1 * x / r;
            c2 += R[2][p] * g2 * y / r;
            c3 += R[3][p] * g3 * z / r;
        }
    }

    c_out[0] = c1 * h3;
    c_out[1] = c2 * h3;
    c_out[2] = c3 * h3;
}

/* --- Simple linear regression: y = slope*x + intercept --- */
static void linregress(const double *x, const double *y, int n,
                       double *slope, double *intercept, double *r2) {
    double sx=0, sy=0, sxx=0, sxy=0;
    for (int i=0; i<n; i++) { sx+=x[i]; sy+=y[i]; sxx+=x[i]*x[i]; sxy+=x[i]*y[i]; }
    double den = n*sxx - sx*sx;
    if (fabs(den) < 1e-30) { *slope = 0; *intercept = 0; *r2 = 0; return; }
    *slope = (n*sxy - sx*sy) / den;
    *intercept = (sy*sxx - sx*sxy) / den;
    double ymean = sy/n, sstot = 0, ssres = 0;
    for (int i=0; i<n; i++) {
        double yp = *slope * x[i] + *intercept;
        sstot += (y[i]-ymean)*(y[i]-ymean);
        ssres += (y[i]-yp)*(y[i]-yp);
    }
    *r2 = (sstot > 0) ? 1.0 - ssres/sstot : 0.0;
}

/* ====================== MAIN ====================== */
int main(int argc, char **argv) {
    /* Defaults */
    int    n = 200, steps = 6000, test = 1;
    double l = 16.0, timestep = 0.0, skyrme = 2.0;
    char   outdir[256] = "v19/data";

    for (int i = 1; i < argc; i++) {
        if      (!strcmp(argv[i],"-N")      && i+1<argc) n = atoi(argv[++i]);
        else if (!strcmp(argv[i],"-L")      && i+1<argc) l = atof(argv[++i]);
        else if (!strcmp(argv[i],"-dt")     && i+1<argc) timestep = atof(argv[++i]);
        else if (!strcmp(argv[i],"-c4")     && i+1<argc) skyrme = atof(argv[++i]);
        else if (!strcmp(argv[i],"-steps")  && i+1<argc) steps = atoi(argv[++i]);
        else if (!strcmp(argv[i],"-test")   && i+1<argc) test = atoi(argv[++i]);
        else if (!strcmp(argv[i],"-A")      && i+1<argc) param_A = atof(argv[++i]);
        else if (!strcmp(argv[i],"-sigma")  && i+1<argc) param_sigma = atof(argv[++i]);
        else if (!strcmp(argv[i],"-R0")     && i+1<argc) param_R0 = atof(argv[++i]);
        else if (!strcmp(argv[i],"-omega1") && i+1<argc) param_omega1 = atof(argv[++i]);
        else if (!strcmp(argv[i],"-omega2") && i+1<argc) param_omega2 = atof(argv[++i]);
        else if (!strcmp(argv[i],"-omega3") && i+1<argc) param_omega3 = atof(argv[++i]);
        else if (!strcmp(argv[i],"-A1")     && i+1<argc) param_A1 = atof(argv[++i]);
        else if (!strcmp(argv[i],"-A2")     && i+1<argc) param_A2 = atof(argv[++i]);
        else if (!strcmp(argv[i],"-A3")     && i+1<argc) param_A3 = atof(argv[++i]);
        else if (!strcmp(argv[i],"-o")      && i+1<argc) strncpy(outdir,argv[++i],255);
    }

    /* Set up grid */
    N  = n;
    L  = l;
    dx = L / (N - 1);
    c4 = skyrme;

    /* CFL: dt <= dx/sqrt(3) for wave equation; use 0.25*CFL for safety */
    double cfl_max = dx / sqrt(3.0);
    if (timestep <= 0) timestep = 0.25 * cfl_max;
    if (timestep > cfl_max) {
        fprintf(stderr, "WARNING: dt=%.6f > CFL=%.6f, clamping to 0.25*CFL\n",
                timestep, cfl_max);
        timestep = 0.25 * cfl_max;
    }
    dt = timestep;

    /* Default frequencies: omega_n = n*pi*c/R0 */
    double w1 = (param_omega1 > 0) ? param_omega1 : M_PI / param_R0;
    double w2 = (param_omega2 > 0) ? param_omega2 : 2.0 * M_PI / param_R0;
    double w3 = (param_omega3 > 0) ? param_omega3 : 3.0 * M_PI / param_R0;

    /* Tests 5,6 override c4 to 0 */
    if (test == 5 || test == 6) c4 = 0.0;

    double mem_gb = (double)(4+4+4+12+1) * N*N*(double)N * 8.0 / 1e9;
    printf("=== Triad 3D (v19) ===\n");
    printf("N=%d  L=%.1f  dx=%.5f  dt=%.6f  c4=%.2f\n", N, L, dx, dt, c4);
    printf("Test %d  steps=%d  threads=%d  mem=%.2f GB\n",
           test, steps, omp_get_max_threads(), mem_gb);

    setlinebuf(stdout);

    alloc_arrays();
    build_damping(0.12);

    /* Core radius for f_core diagnostic */
    double R_core = 2.0 * param_sigma;

    /* Initialize field based on test number */
    const char *test_names[] = {
        "", "B=0 Hedgehog Bump", "Three-Mode Kicked Triad",
        "Three-Mode Displaced Triad", "Single-Mode Control",
        "No-Skyrme Bump", "No-Skyrme Kicked Triad"
    };
    if (test < 1 || test > 6) {
        fprintf(stderr, "Unknown test %d (valid: 1-6)\n", test);
        return 1;
    }
    printf("\n--- Test %d: %s ---\n", test, test_names[test]);

    switch (test) {
        case 1:
            printf("A=%.2f sigma=%.2f\n", param_A, param_sigma);
            init_bump(param_A, param_sigma);
            break;
        case 2:
            printf("A=%.2f sigma=%.2f R0=%.2f w1=%.4f w2=%.4f w3=%.4f\n",
                   param_A, param_sigma, param_R0, w1, w2, w3);
            init_kicked_triad(param_A, param_sigma, param_R0, w1, w2, w3);
            break;
        case 3:
            printf("A1=%.2f A2=%.2f A3=%.2f sigma=%.2f R0=%.2f\n",
                   param_A1, param_A2, param_A3, param_sigma, param_R0);
            init_displaced_triad(param_A1, param_A2, param_A3,
                                 param_sigma, param_R0);
            break;
        case 4:
            printf("A=%.2f sigma=%.2f R0=%.2f\n",
                   param_A3, param_sigma, param_R0);
            init_single_mode(param_A3, param_sigma, param_R0);
            break;
        case 5:
            printf("A=%.2f sigma=%.2f (c4=0 override)\n", param_A, param_sigma);
            init_bump(param_A, param_sigma);
            break;
        case 6:
            printf("A=%.2f sigma=%.2f R0=%.2f (c4=0 override)\n",
                   param_A, param_sigma, param_R0);
            init_kicked_triad(param_A, param_sigma, param_R0, w1, w2, w3);
            break;
    }

    /* Initial diagnostics */
    double Ek, E2, E4;
    compute_energy(&Ek, &E2, &E4);
    double B = compute_baryon();
    double Etot = Ek + E2 + E4;
    double E_mode[3];
    compute_mode_energy(E_mode);
    double f_core = compute_core_fraction(R_core);

    printf("Initial: Ek=%.4f E2=%.4f E4=%.4f Etot=%.4f  B=%.6f\n",
           Ek, E2, E4, Etot, B);
    printf("Mode energies: E1=%.4f E2=%.4f E3=%.4f  f_core=%.4f\n",
           E_mode[0], E_mode[1], E_mode[2], f_core);

    /* --- Output files --- */
    char fname_ts[512], fname_modes[512];
    snprintf(fname_ts, sizeof(fname_ts), "%s/test%d_timeseries.tsv", outdir, test);
    snprintf(fname_modes, sizeof(fname_modes), "%s/test%d_modes.tsv", outdir, test);

    FILE *fp_ts = fopen(fname_ts, "w");
    if (!fp_ts) {
        fprintf(stderr, "Cannot open %s. Creating output directory.\n", fname_ts);
        char cmd[512];
        snprintf(cmd, sizeof(cmd), "mkdir -p %s", outdir);
        if (system(cmd) != 0) { fprintf(stderr, "mkdir failed\n"); return 1; }
        fp_ts = fopen(fname_ts, "w");
        if (!fp_ts) { perror("fopen timeseries"); return 1; }
    }
    FILE *fp_modes = fopen(fname_modes, "w");
    if (!fp_modes) { perror("fopen modes"); return 1; }

    fprintf(fp_ts, "step\ttime\tE_kin\tE2\tE4\tE_total\tB\tE1_mode\tE2_mode\tE3_mode\tf_core\n");
    fprintf(fp_ts, "%d\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n",
            0, 0.0, Ek, E2, E4, Etot, B, E_mode[0], E_mode[1], E_mode[2], f_core);

    fprintf(fp_modes, "step\ttime\tc1\tc2\tc3\tPhi\n");
    double c_proj[3];
    compute_mode_projection(c_proj, param_sigma, param_R0);
    fprintf(fp_modes, "%d\t%.6f\t%.8e\t%.8e\t%.8e\t%.4f\n",
            0, 0.0, c_proj[0], c_proj[1], c_proj[2], 0.0);

    /* --- Time-averaged shell accumulator --- */
    #define NBINS 100
    double tavg_delta[NBINS], tavg_rhoE[NBINS];
    int    tavg_cnt[NBINS];
    for (int b = 0; b < NBINS; b++) {
        tavg_delta[b] = tavg_rhoE[b] = 0.0;
        tavg_cnt[b] = 0;
    }
    int n_tavg = 0;
    int settle_step = (int)(0.2 * steps);
    double rmax_bin = L / 2.0;
    double dr_bin = rmax_bin / NBINS;

    /* --- q1(0,t) recording for FFT --- */
    int print_every = steps / 60;
    if (print_every < 1) print_every = 1;
    int max_records = steps / print_every + 2;
    double *q1_record = malloc(max_records * sizeof(double));
    double *t_record  = malloc(max_records * sizeof(double));
    int n_record = 0;

    /* Record initial q1(0) */
    int origin_idx = flat(N/2, N/2, N/2);
    q1_record[n_record] = R[1][origin_idx];
    t_record[n_record]  = 0.0;
    n_record++;

    /* Mode projection history for phase computation */
    double *c_hist[3];
    for (int a = 0; a < 3; a++)
        c_hist[a] = malloc(max_records * sizeof(double));
    double *t_hist = malloc(max_records * sizeof(double));
    int n_hist = 0;
    for (int a = 0; a < 3; a++) c_hist[a][0] = c_proj[a];
    t_hist[0] = 0.0;
    n_hist = 1;

    /* --- Time evolution --- */
    double t_start = omp_get_wtime();

    for (int s = 1; s <= steps; s++) {
        step_vv();

        if (s % print_every == 0 || s == steps) {
            double t = s * dt;

            compute_energy(&Ek, &E2, &E4);
            B = compute_baryon();
            Etot = Ek + E2 + E4;
            compute_mode_energy(E_mode);
            f_core = compute_core_fraction(R_core);
            compute_mode_projection(c_proj, param_sigma, param_R0);

            /* Record q1(origin) for FFT */
            if (n_record < max_records) {
                q1_record[n_record] = R[1][origin_idx];
                t_record[n_record]  = t;
                n_record++;
            }

            /* Record mode projection for phase lock */
            if (n_hist < max_records) {
                for (int a = 0; a < 3; a++) c_hist[a][n_hist] = c_proj[a];
                t_hist[n_hist] = t;
                n_hist++;
            }

            /* Compute phase lock diagnostic Phi = phi3 - phi1 - phi2 */
            double Phi = 0.0;
            if (n_hist >= 3) {
                int m = n_hist - 1;
                double dt_diag = t_hist[m] - t_hist[m-1];
                if (dt_diag > 1e-15) {
                    double phi[3];
                    for (int a = 0; a < 3; a++) {
                        double c_now  = c_hist[a][m];
                        double c_prev = c_hist[a][m-1];
                        double cdot = (c_now - c_prev) / dt_diag;
                        double wa = (a == 0) ? w1 : (a == 1) ? w2 : w3;
                        phi[a] = atan2(-cdot / wa, c_now);
                    }
                    Phi = phi[2] - phi[0] - phi[1];
                    /* Wrap to [-pi, pi] */
                    while (Phi >  M_PI) Phi -= 2*M_PI;
                    while (Phi < -M_PI) Phi += 2*M_PI;
                }
            }

            /* Write time series */
            fprintf(fp_ts, "%d\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n",
                    s, t, Ek, E2, E4, Etot, B,
                    E_mode[0], E_mode[1], E_mode[2], f_core);
            fflush(fp_ts);

            /* Write mode decomposition */
            fprintf(fp_modes, "%d\t%.6f\t%.8e\t%.8e\t%.8e\t%.4f\n",
                    s, t, c_proj[0], c_proj[1], c_proj[2], Phi);
            fflush(fp_modes);

            /* Time-averaged shell accumulation (after settling) */
            if (s > settle_step) {
                double idx2_loc = 1.0 / (2.0 * dx);

                /* Accumulate shell averages: delta = 1-q0, rhoE = energy density */
                /* Use per-bin counters for this diagnostic step */
                double step_delta[NBINS], step_rhoE[NBINS];
                int    step_cnt[NBINS];
                for (int b = 0; b < NBINS; b++) {
                    step_delta[b] = step_rhoE[b] = 0.0;
                    step_cnt[b] = 0;
                }

                for (int ii = 0; ii < N; ii++)
                for (int jj = 0; jj < N; jj++)
                for (int kk = 0; kk < N; kk++) {
                    double x = coord(ii), y = coord(jj), z = coord(kk);
                    double r = sqrt(x*x + y*y + z*z);
                    int b = (int)(r / dr_bin);
                    if (b < 0 || b >= NBINS) continue;

                    int pp = ii*N*N + jj*N + kk;
                    step_delta[b] += 1.0 - R[0][pp];

                    /* Energy density */
                    int nbl[3][2];
                    nbl[0][0] = flat(ii+1,jj,kk); nbl[0][1] = flat(ii-1,jj,kk);
                    nbl[1][0] = flat(ii,jj+1,kk); nbl[1][1] = flat(ii,jj-1,kk);
                    nbl[2][0] = flat(ii,jj,kk+1); nbl[2][1] = flat(ii,jj,kk-1);

                    double v2l = 0;
                    for (int c = 0; c < 4; c++) v2l += V[c][pp]*V[c][pp];

                    double dql[3][4];
                    for (int d = 0; d < 3; d++)
                        for (int c = 0; c < 4; c++)
                            dql[d][c] = (R[c][nbl[d][0]] - R[c][nbl[d][1]]) * idx2_loc;

                    double Dl[3][3];
                    for (int aa = 0; aa < 3; aa++)
                    for (int bb = aa; bb < 3; bb++) {
                        double ss = 0;
                        for (int c = 0; c < 4; c++) ss += dql[aa][c]*dql[bb][c];
                        Dl[aa][bb] = ss;
                        if (bb > aa) Dl[bb][aa] = ss;
                    }
                    double trDl = Dl[0][0]+Dl[1][1]+Dl[2][2];
                    double trD2l = 0;
                    for (int aa = 0; aa < 3; aa++)
                        for (int bb = 0; bb < 3; bb++)
                            trD2l += Dl[aa][bb]*Dl[aa][bb];

                    double ed = 0.5*v2l + 0.5*trDl + 0.25*c4*(trDl*trDl - trD2l);
                    step_rhoE[b] += ed;
                    step_cnt[b]++;
                }

                /* Add shell averages to accumulator */
                for (int b = 0; b < NBINS; b++) {
                    if (step_cnt[b] > 0) {
                        tavg_delta[b] += step_delta[b] / step_cnt[b];
                        tavg_rhoE[b]  += step_rhoE[b] / step_cnt[b];
                        tavg_cnt[b]++;
                    }
                }
                n_tavg++;
            }

            /* Console output */
            double elapsed = omp_get_wtime() - t_start;
            double rate = s / elapsed;
            printf("Step %5d | t=%7.3f | Et=%8.3f | B=%8.6f | E1=%6.3f E2=%6.3f E3=%6.3f | fc=%5.3f | Phi=%6.3f | %.0f/s\n",
                   s, t, Etot, B,
                   E_mode[0], E_mode[1], E_mode[2], f_core, Phi, rate);
        }
    }

    fclose(fp_ts);
    fclose(fp_modes);

    double elapsed = omp_get_wtime() - t_start;
    printf("\nDone: %d steps in %.1f s (%.1f steps/s)\n", steps, elapsed, steps/elapsed);
    printf("Timeseries: %s\n", fname_ts);
    printf("Modes: %s\n", fname_modes);

    /* ====================== POST-PROCESSING ====================== */

    /* --- Write time-averaged radial profile --- */
    char fname_tavg[512];
    snprintf(fname_tavg, sizeof(fname_tavg), "%s/test%d_tavg_profile.tsv", outdir, test);
    FILE *fp_tavg = fopen(fname_tavg, "w");
    if (fp_tavg) {
        fprintf(fp_tavg, "r\tavg_delta\tavg_rhoE\n");
        for (int b = 0; b < NBINS; b++) {
            if (tavg_cnt[b] > 0) {
                double r_mid = (b + 0.5) * dr_bin;
                fprintf(fp_tavg, "%.4f\t%.8e\t%.8e\n",
                        r_mid,
                        tavg_delta[b] / tavg_cnt[b],
                        tavg_rhoE[b] / tavg_cnt[b]);
            }
        }
        fclose(fp_tavg);
        printf("Time-averaged profile: %s\n", fname_tavg);
    }

    /* --- Far-field power-law fits on time-averaged profile --- */
    char fname_far[512];
    snprintf(fname_far, sizeof(fname_far), "%s/test%d_farfield.txt", outdir, test);
    FILE *fp_far = fopen(fname_far, "w");
    if (fp_far) {
        double fit_rmin = 2.0 * R_core;  /* Outside the core */
        double fit_rmax = L/2 - 0.12*L;  /* Inside absorbing layer */

        double logr[NBINS], logd[NBINS], logE[NBINS];
        int npts = 0;
        for (int b = 0; b < NBINS; b++) {
            if (tavg_cnt[b] < 1) continue;
            double r_mid = (b + 0.5) * dr_bin;
            if (r_mid < fit_rmin || r_mid > fit_rmax) continue;
            double avg_d = tavg_delta[b] / tavg_cnt[b];
            double avg_e = tavg_rhoE[b] / tavg_cnt[b];
            if (avg_d < 1e-20 || avg_e < 1e-20) continue;
            logr[npts] = log(r_mid);
            logd[npts] = log(avg_d);
            logE[npts] = log(avg_e);
            npts++;
        }

        fprintf(fp_far, "Far-field power-law fits (r in [%.1f, %.1f], %d bins, %d time samples)\n\n",
                fit_rmin, fit_rmax, npts, n_tavg);

        if (npts >= 5) {
            double sl, inter, r2;

            linregress(logr, logd, npts, &sl, &inter, &r2);
            fprintf(fp_far, "avg_delta ~ r^{%.3f}  R^2=%.4f  (free radiation: -2)\n", sl, r2);
            printf("\nFar-field: avg_delta ~ r^{%.3f} R^2=%.4f\n", sl, r2);

            linregress(logr, logE, npts, &sl, &inter, &r2);
            fprintf(fp_far, "avg_rhoE  ~ r^{%.3f}  R^2=%.4f  (free radiation: -2)\n", sl, r2);
            printf("Far-field: avg_rhoE  ~ r^{%.3f} R^2=%.4f\n", sl, r2);
        } else {
            fprintf(fp_far, "Not enough data for fit (%d points, need >= 5)\n", npts);
            printf("Not enough far-field data for power-law fit (%d points)\n", npts);
        }
        fclose(fp_far);
        printf("Far-field fits: %s\n", fname_far);
    }

    /* --- DFT of q1(origin, t) --- */
    char fname_spec[512];
    snprintf(fname_spec, sizeof(fname_spec), "%s/test%d_spectrum.tsv", outdir, test);
    FILE *fp_spec = fopen(fname_spec, "w");
    if (fp_spec && n_record >= 3) {
        double T_total = t_record[n_record - 1] - t_record[0];
        double dt_diag = T_total / (n_record - 1);
        double omega_max = M_PI / dt_diag;
        double d_omega = 2.0 * M_PI / T_total;
        int n_freq = (int)(omega_max / d_omega) + 1;
        if (n_freq > 2000) n_freq = 2000;

        fprintf(fp_spec, "omega\tpower\n");

        for (int f = 0; f < n_freq; f++) {
            double omega = f * d_omega;
            double re = 0, im = 0;
            for (int m = 0; m < n_record; m++) {
                double phase = omega * t_record[m];
                re += q1_record[m] * cos(phase);
                im -= q1_record[m] * sin(phase);
            }
            double power = (re*re + im*im) / (n_record * n_record);
            fprintf(fp_spec, "%.6f\t%.8e\n", omega, power);
        }
        fclose(fp_spec);
        printf("Spectrum: %s\n", fname_spec);

        /* Report peak frequency */
        double peak_omega = 0, peak_power = 0;
        for (int f = 1; f < n_freq; f++) {
            double omega = f * d_omega;
            double re = 0, im = 0;
            for (int m = 0; m < n_record; m++) {
                double phase = omega * t_record[m];
                re += q1_record[m] * cos(phase);
                im -= q1_record[m] * sin(phase);
            }
            double power = (re*re + im*im) / (n_record * n_record);
            if (power > peak_power) {
                peak_power = power;
                peak_omega = omega;
            }
        }
        printf("Spectrum peak: omega=%.4f (w1=%.4f, w2=%.4f, w3=%.4f)\n",
               peak_omega, w1, w2, w3);
    } else if (fp_spec) {
        fprintf(fp_spec, "omega\tpower\n");
        fclose(fp_spec);
    }

    /* --- Cleanup --- */
    free(q1_record);
    free(t_record);
    for (int a = 0; a < 3; a++) free(c_hist[a]);
    free(t_hist);
    free_arrays();

    return 0;
}
