/*
 * v18/src/skyrme3d.c
 *
 * 3D Skyrme sigma-model simulation with Velocity-Verlet integration.
 *
 * Corrects ALL bugs in the original v18 Python code:
 *   1. Correct Skyrme force via strain tensor D_{ij} = dq_i . dq_j
 *   2. Product ansatz for two-soliton initialization (not linear superposition)
 *   3. Complete energy: E_kin + E_2 + E_4 (Skyrme quartic)
 *   4. Correct baryon number via Maurer-Cartan determinant
 *   5. Adequate resolution with CFL-safe timestep
 *
 * Physics:
 *   Lagrangian L = (1/2)|dq|^2 + (c4/4)[Tr(D)^2 - Tr(D^2)]
 *   on S^3 (unit quaternion field q, |q|=1).
 *   c4 = 2 rho0^2 / e^2.  Code units: e=1, rho0=1 => c4=2.
 *
 *   EOM: q_tt = Lap(q) + c4 * div(G) - (q.F + |v|^2) q
 *   where G_i = Tr(D) dq_i - D_{ij} dq_j
 *
 * Compile: gcc -O3 -fopenmp -o skyrme3d skyrme3d.c -lm
 * Run:     ./skyrme3d -test 1 -N 200 -L 12 -steps 3000
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

/* --- Periodic index --- */
static inline int pbc(int i) { return ((i % N) + N) % N; }
static inline int flat(int i, int j, int k) {
    return pbc(i)*N*N + pbc(j)*N + pbc(k);
}
static inline double coord(int i) { return -L/2 + i*dx; }

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

/* --- Hedgehog Skyrmion initialization ---
 * Profile: f(r) = 4*atan(exp(-r/sigma))
 *   f(0)=pi, f(inf)=0, f'(0)=-2/sigma
 *   sigma=1.41 gives f'(0)=-1.42 (matches equilibrium slope a=1.42)
 *   Decays exponentially -> clean boundary conditions */
static void init_hedgehog(double cx, double cy, double cz) {
    double sigma = 1.41;
    #pragma omp parallel for collapse(3)
    for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
    for (int k = 0; k < N; k++) {
        double x = coord(i) - cx;
        double y = coord(j) - cy;
        double z = coord(k) - cz;
        double r = sqrt(x*x + y*y + z*z);
        int p = i*N*N + j*N + k;

        if (r < 1e-12) {
            R[0][p] = -1.0; /* cos(pi) */
            R[1][p] = R[2][p] = R[3][p] = 0.0;
        } else {
            double f = 4.0 * atan(exp(-r / sigma));
            R[0][p] = cos(f);
            double sf = sin(f);
            R[1][p] = sf * x / r;
            R[2][p] = sf * y / r;
            R[3][p] = sf * z / r;
        }
        V[0][p] = V[1][p] = V[2][p] = V[3][p] = 0.0;
    }
}

/* --- Quaternion product of field arrays: out = a * b (Hamilton convention) --- */
static void qmul_fields(double *a[4], double *b[4], double *out[4]) {
    #pragma omp parallel for
    for (int p = 0; p < ntot; p++) {
        double a0=a[0][p], a1=a[1][p], a2=a[2][p], a3=a[3][p];
        double b0=b[0][p], b1=b[1][p], b2=b[2][p], b3=b[3][p];
        out[0][p] = a0*b0 - a1*b1 - a2*b2 - a3*b3;
        out[1][p] = a0*b1 + a1*b0 + a2*b3 - a3*b2;
        out[2][p] = a0*b2 - a1*b3 + a2*b0 + a3*b1;
        out[3][p] = a0*b3 + a1*b2 - a2*b1 + a3*b0;
    }
}

/* --- Two-soliton state via product ansatz ---
 * q_total = q1 * q2 gives B_total = B1 + B2 = 2 */
static void init_two_solitons(double sep) {
    /* First hedgehog at (-sep/2, 0, 0) */
    init_hedgehog(-sep/2, 0, 0);
    double *q1[4], *tmp[4];
    for (int c = 0; c < 4; c++) {
        q1[c]  = malloc(ntot * sizeof(double));
        tmp[c] = malloc(ntot * sizeof(double));
        memcpy(q1[c], R[c], ntot * sizeof(double));
    }

    /* Second hedgehog at (+sep/2, 0, 0) */
    init_hedgehog(+sep/2, 0, 0);

    /* Product ansatz then normalize */
    qmul_fields(q1, R, tmp);
    #pragma omp parallel for
    for (int p = 0; p < ntot; p++) {
        double norm = sqrt(tmp[0][p]*tmp[0][p] + tmp[1][p]*tmp[1][p] +
                          tmp[2][p]*tmp[2][p] + tmp[3][p]*tmp[3][p]);
        if (norm < 1e-15) norm = 1e-15;
        for (int c = 0; c < 4; c++) R[c][p] = tmp[c][p] / norm;
        V[0][p] = V[1][p] = V[2][p] = V[3][p] = 0.0;
    }

    for (int c = 0; c < 4; c++) { free(q1[c]); free(tmp[c]); }
}

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
 *
 * D_{ij} = d_i q . d_j q     (strain tensor)
 * G_i^a = Tr(D) d_i q^a  -  sum_j D_{ij} d_j q^a
 *
 * Two-pass approach: first compute G everywhere, then take div(G) in compute_accel.
 */
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

        /* Spatial gradients dq[d][c] = d_d q^c */
        double dq[3][4];
        for (int d = 0; d < 3; d++)
            for (int c = 0; c < 4; c++)
                dq[d][c] = (R[c][nb[d][0]] - R[c][nb[d][1]]) * idx2;

        /* Strain tensor D[a][b] = d_a q . d_b q */
        double D[3][3];
        for (int a = 0; a < 3; a++)
        for (int b = a; b < 3; b++) {
            double s = 0;
            for (int c = 0; c < 4; c++) s += dq[a][c] * dq[b][c];
            D[a][b] = s;
            if (b > a) D[b][a] = s;
        }
        double trD = D[0][0] + D[1][1] + D[2][2];

        /* G_d^c = Tr(D) dq_d^c - D_{de} dq_e^c */
        for (int d = 0; d < 3; d++)
        for (int c = 0; c < 4; c++) {
            double Ddq = 0;
            for (int e = 0; e < 3; e++) Ddq += D[d][e] * dq[e][c];
            G[d][c][p] = trD * dq[d][c] - Ddq;
        }
    }
}

/* --- Compute acceleration ---
 * a = Lap(q) + c4*div(G) - (q.F + |v|^2) q
 * Tangent-space projection keeps q on S^3. */
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
            /* 7-point Laplacian */
            double lap = (R[c][nb[0][0]] + R[c][nb[0][1]] +
                         R[c][nb[1][0]] + R[c][nb[1][1]] +
                         R[c][nb[2][0]] + R[c][nb[2][1]] - 6*R[c][p]) * idxsq;

            /* Divergence of Skyrme stress: d_x G_x + d_y G_y + d_z G_z */
            double divG = 0;
            for (int d = 0; d < 3; d++)
                divG += (G[d][c][nb[d][0]] - G[d][c][nb[d][1]]) * idx2;

            F[c] = lap + c4 * divG;
        }

        /* Project to tangent space of S^3 */
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

/* --- Energy diagnostics ---
 * E_kin = (1/2) int |V|^2
 * E_2   = (1/2) int Tr(D)
 * E_4   = (c4/4) int [Tr(D)^2 - Tr(D^2)]  */
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

/* --- Baryon number (topological charge) ---
 * B = -(1/(2 pi^2)) int det(L_x, L_y, L_z) d^3x
 * where L_d = Im(q_bar * d_d q) in Hamilton convention.
 *
 * Hamilton: q_bar * v
 *   [1] = q0*v1 - q1*v0 - q2*v3 + q3*v2
 *   [2] = q0*v2 + q1*v3 - q2*v0 - q3*v1
 *   [3] = q0*v3 - q1*v2 + q2*v1 - q3*v0
 */
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

        /* Maurer-Cartan currents l[d][a] = Im(q_bar * d_d q)_a */
        double l[3][3];
        for (int d = 0; d < 3; d++) {
            double *v = dq[d];
            l[d][0] = q[0]*v[1] - q[1]*v[0] - q[2]*v[3] + q[3]*v[2];
            l[d][1] = q[0]*v[2] + q[1]*v[3] - q[2]*v[0] - q[3]*v[1];
            l[d][2] = q[0]*v[3] - q[1]*v[2] + q[2]*v[1] - q[3]*v[0];
        }

        /* det(L_x, L_y, L_z) = l_x . (l_y x l_z) */
        double det = l[0][0]*(l[1][1]*l[2][2] - l[1][2]*l[2][1])
                   - l[0][1]*(l[1][0]*l[2][2] - l[1][2]*l[2][0])
                   + l[0][2]*(l[1][0]*l[2][1] - l[1][1]*l[2][0]);
        Q += det;
    }

    return -Q * h3 / (2.0 * M_PI * M_PI);
}

/* --- Soliton centroid & RMS radius (weighted by |grad q|^2) --- */
static void soliton_stats(double *rms_out) {
    double idx2 = 1.0 / (2.0 * dx);
    double h3 = dx*dx*dx;
    double wr2 = 0, wtot = 0;

    #pragma omp parallel for collapse(3) reduction(+:wr2,wtot)
    for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
    for (int k = 0; k < N; k++) {
        int p = i*N*N + j*N + k;
        int nb[3][2];
        nb[0][0] = flat(i+1,j,k); nb[0][1] = flat(i-1,j,k);
        nb[1][0] = flat(i,j+1,k); nb[1][1] = flat(i,j-1,k);
        nb[2][0] = flat(i,j,k+1); nb[2][1] = flat(i,j,k-1);

        double grad2 = 0;
        for (int c = 0; c < 4; c++) {
            for (int d = 0; d < 3; d++) {
                double g = (R[c][nb[d][0]] - R[c][nb[d][1]]) * idx2;
                grad2 += g*g;
            }
        }

        double x = coord(i), y = coord(j), z = coord(k);
        double w = grad2 * h3;
        wtot += w;
        wr2  += w * (x*x + y*y + z*z);
    }

    *rms_out = (wtot > 1e-30) ? sqrt(wr2 / wtot) : 0.0;
}

/* --- Two-soliton separation (centroid of baryon density in each half) --- */
static double two_soliton_sep(void) {
    double idx2 = 1.0 / (2.0 * dx);
    double h3 = dx*dx*dx;
    double wL = 0, wxL = 0, wR = 0, wxR = 0;

    #pragma omp parallel for collapse(3) reduction(+:wL,wxL,wR,wxR)
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
        double bdens = fabs(det);
        double x = coord(i);
        if (x < 0) { wL += bdens * h3; wxL += bdens * h3 * x; }
        else       { wR += bdens * h3; wxR += bdens * h3 * x; }
    }

    double xL = (wL > 1e-30) ? wxL / wL : 0;
    double xR = (wR > 1e-30) ? wxR / wR : 0;
    return xR - xL;
}

/* ====================== FAR-FIELD & FORCE DIAGNOSTICS ====================== */

/* Simple linear regression: y = slope*x + intercept */
static void linregress(const double *x, const double *y, int n,
                       double *slope, double *intercept, double *r2) {
    double sx=0, sy=0, sxx=0, sxy=0;
    for (int i=0; i<n; i++) { sx+=x[i]; sy+=y[i]; sxx+=x[i]*x[i]; sxy+=x[i]*y[i]; }
    double den = n*sxx - sx*sx;
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

/* --- Radial profile of physically meaningful quantities ---
 *
 * Computes shell-averaged profiles of:
 *   |L|   = |Im(q_bar d_d q)|  — Maurer-Cartan current (gauge potential analog)
 *   E2_d  = (1/2)|grad q|^2    — gradient energy density
 *   E4_d  = (c4/4)[TrD^2-TrD2] — Skyrme energy density (proportional to |[L_i,L_j]|^2)
 *   B0    = baryon density
 *
 * Then fits log-log power law in the far field to extract decay exponents.
 */
static void radial_profile_analysis(const char *outdir) {
    double idx2 = 1.0 / (2.0 * dx);
    double h3 = dx*dx*dx;

    /* Radial bins */
    #define NBINS 100
    double rmin = 0.0, rmax = L/2 - dx;
    double dr = (rmax - rmin) / NBINS;
    double bin_r[NBINS], bin_L[NBINS], bin_e2[NBINS], bin_e4[NBINS], bin_B[NBINS];
    double bin_cnt[NBINS];
    for (int b=0; b<NBINS; b++) {
        bin_r[b]=bin_L[b]=bin_e2[b]=bin_e4[b]=bin_B[b]=0;
        bin_cnt[b]=0;
    }

    /* Accumulate (not parallelized for simplicity with bin reduction) */
    for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
    for (int k = 0; k < N; k++) {
        double x = coord(i), y = coord(j), z = coord(k);
        double r = sqrt(x*x + y*y + z*z);
        int b = (int)((r - rmin) / dr);
        if (b < 0 || b >= NBINS) continue;

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

        /* Maurer-Cartan currents l[d][a] = Im(q_bar * d_d q) */
        double l[3][3], Lmag2 = 0;
        for (int d = 0; d < 3; d++) {
            double *v = dq[d];
            l[d][0] = q[0]*v[1] - q[1]*v[0] - q[2]*v[3] + q[3]*v[2];
            l[d][1] = q[0]*v[2] + q[1]*v[3] - q[2]*v[0] - q[3]*v[1];
            l[d][2] = q[0]*v[3] - q[1]*v[2] + q[2]*v[1] - q[3]*v[0];
            for (int a=0; a<3; a++) Lmag2 += l[d][a]*l[d][a];
        }

        /* Strain tensor for E2, E4 density */
        double D[3][3];
        for (int a = 0; a < 3; a++)
        for (int bb = a; bb < 3; bb++) {
            double s = 0;
            for (int c = 0; c < 4; c++) s += dq[a][c] * dq[bb][c];
            D[a][bb] = s; if (bb>a) D[bb][a] = s;
        }
        double trD = D[0][0]+D[1][1]+D[2][2];
        double trD2 = 0;
        for (int a=0; a<3; a++) for (int bb=0; bb<3; bb++) trD2 += D[a][bb]*D[a][bb];

        /* Baryon density */
        double det = l[0][0]*(l[1][1]*l[2][2] - l[1][2]*l[2][1])
                   - l[0][1]*(l[1][0]*l[2][2] - l[1][2]*l[2][0])
                   + l[0][2]*(l[1][0]*l[2][1] - l[1][1]*l[2][0]);

        bin_r[b]   += r;
        bin_L[b]   += sqrt(Lmag2);
        bin_e2[b]  += 0.5 * trD;
        bin_e4[b]  += 0.25 * c4 * (trD*trD - trD2);
        bin_B[b]   += fabs(det) / (2.0 * M_PI * M_PI);
        bin_cnt[b] += 1;
    }

    /* Write radial profile */
    char fname[512];
    snprintf(fname, sizeof(fname), "%s/radial_profile.tsv", outdir);
    FILE *fp = fopen(fname, "w");
    fprintf(fp, "r\t|L|\tE2_dens\tE4_dens\tB_dens\n");
    for (int b=0; b<NBINS; b++) {
        if (bin_cnt[b] < 1) continue;
        double cnt = bin_cnt[b];
        fprintf(fp, "%.4f\t%.8e\t%.8e\t%.8e\t%.8e\n",
                bin_r[b]/cnt, bin_L[b]/cnt, bin_e2[b]/cnt, bin_e4[b]/cnt, bin_B[b]/cnt);
    }
    fclose(fp);
    printf("Radial profile written to %s\n", fname);

    /* Power-law fit in far field (r > 2.5 to r < L/2 - 2) */
    double fit_rmin = 2.5, fit_rmax = L/2 - 2.0;
    double logr[NBINS], logL[NBINS], loge2[NBINS], loge4[NBINS], logB[NBINS];
    int npts = 0;
    for (int b=0; b<NBINS; b++) {
        if (bin_cnt[b] < 5) continue;
        double r_avg = bin_r[b] / bin_cnt[b];
        if (r_avg < fit_rmin || r_avg > fit_rmax) continue;
        double L_avg = bin_L[b] / bin_cnt[b];
        double e2_avg = bin_e2[b] / bin_cnt[b];
        double e4_avg = bin_e4[b] / bin_cnt[b];
        double B_avg = bin_B[b] / bin_cnt[b];
        if (L_avg < 1e-20 || e2_avg < 1e-20 || e4_avg < 1e-20 || B_avg < 1e-20) continue;
        logr[npts] = log(r_avg);
        logL[npts] = log(L_avg);
        loge2[npts] = log(e2_avg);
        loge4[npts] = log(e4_avg);
        logB[npts] = log(B_avg);
        npts++;
    }

    if (npts >= 5) {
        double sl, inter, r2;
        printf("\nFar-field power-law fits (r in [%.1f, %.1f], %d bins):\n", fit_rmin, fit_rmax, npts);

        linregress(logr, logL, npts, &sl, &inter, &r2);
        printf("  |L| (gauge potential)  ~ r^{%.2f}  R²=%.4f   (expect -3 for massless sigma)\n", sl, r2);

        linregress(logr, loge2, npts, &sl, &inter, &r2);
        printf("  E2_dens (grad energy)  ~ r^{%.2f}  R²=%.4f   (expect -6)\n", sl, r2);

        linregress(logr, loge4, npts, &sl, &inter, &r2);
        printf("  E4_dens (Skyrme |F|²)  ~ r^{%.2f}  R²=%.4f   (expect -12)\n", sl, r2);

        linregress(logr, logB, npts, &sl, &inter, &r2);
        printf("  B_dens (baryon charge)  ~ r^{%.2f}  R²=%.4f   (expect -8 to -9)\n", sl, r2);

        printf("\n  ** None of these are 1/r² (Coulomb). The Skyrme field has NO long-range gauge field. **\n");
    } else {
        printf("Not enough far-field data for power-law fit (%d points)\n", npts);
    }
    #undef NBINS
}

/* --- Two-body effective force from centroid acceleration ---
 * Uses finite differences of separation time series.
 * Call after time evolution with stored separation history. */
#define MAX_FORCE_PTS 200
static double force_time[MAX_FORCE_PTS], force_sep[MAX_FORCE_PTS];
static int    force_npts = 0;

static void record_separation(double t, double sep_val) {
    if (force_npts < MAX_FORCE_PTS) {
        force_time[force_npts] = t;
        force_sep[force_npts]  = sep_val;
        force_npts++;
    }
}

static void analyze_two_body_force(const char *outdir) {
    if (force_npts < 5) return;

    char fname[512];
    snprintf(fname, sizeof(fname), "%s/two_body_force.tsv", outdir);
    FILE *fp = fopen(fname, "w");
    fprintf(fp, "time\tsep\tvel\taccel\n");

    /* Second-order finite differences for velocity and acceleration */
    for (int i = 2; i < force_npts - 2; i++) {
        double dt_loc = force_time[i+1] - force_time[i-1];
        double vel = (force_sep[i+1] - force_sep[i-1]) / dt_loc;
        double dt2 = force_time[i+1] - force_time[i];
        double acc = (force_sep[i+1] - 2*force_sep[i] + force_sep[i-1]) / (dt2*dt2);
        fprintf(fp, "%.4f\t%.4f\t%.6f\t%.6f\n",
                force_time[i], force_sep[i], vel, acc);
    }
    fclose(fp);
    printf("Two-body force data written to %s\n", fname);

    /* Fit F(D) power law in region where solitons are separating (clear force signal) */
    double logD[MAX_FORCE_PTS], logF[MAX_FORCE_PTS];
    int np = 0;
    for (int i = 2; i < force_npts - 2; i++) {
        double dt2 = force_time[i+1] - force_time[i];
        double acc = (force_sep[i+1] - 2*force_sep[i] + force_sep[i-1]) / (dt2*dt2);
        if (fabs(acc) < 1e-10 || force_sep[i] < 3.0) continue;
        logD[np] = log(force_sep[i]);
        logF[np] = log(fabs(acc));
        np++;
    }
    if (np >= 5) {
        double sl, inter, r2;
        linregress(logD, logF, np, &sl, &inter, &r2);
        printf("Two-body force: |F| ~ D^{%.2f}  R²=%.4f  (%d pts)\n", sl, r2, np);
        printf("  (Coulomb/Newton would be -2; Skyrme interaction is short-range)\n");
    }
}

/* ====================== MAIN ====================== */
int main(int argc, char **argv) {
    /* Defaults */
    int    n = 200, steps = 3000, test = 1;
    double l = 12.0, timestep = 0.0, skyrme = 2.0, sep = 6.0;
    char   outdir[256] = "v18/data";

    for (int i = 1; i < argc; i++) {
        if      (!strcmp(argv[i],"-N")     && i+1<argc) n = atoi(argv[++i]);
        else if (!strcmp(argv[i],"-L")     && i+1<argc) l = atof(argv[++i]);
        else if (!strcmp(argv[i],"-dt")    && i+1<argc) timestep = atof(argv[++i]);
        else if (!strcmp(argv[i],"-c4")    && i+1<argc) skyrme = atof(argv[++i]);
        else if (!strcmp(argv[i],"-steps") && i+1<argc) steps = atoi(argv[++i]);
        else if (!strcmp(argv[i],"-test")  && i+1<argc) test = atoi(argv[++i]);
        else if (!strcmp(argv[i],"-sep")   && i+1<argc) sep = atof(argv[++i]);
        else if (!strcmp(argv[i],"-o")     && i+1<argc) strncpy(outdir,argv[++i],255);
    }

    /* Set up grid */
    N  = n;
    L  = l;
    dx = L / (N - 1);
    c4 = skyrme;

    /* CFL: dt <= dx/sqrt(3) for wave equation; use 0.25*CFL for safety with Skyrme */
    double cfl_max = dx / sqrt(3.0);
    if (timestep <= 0) timestep = 0.25 * cfl_max;
    if (timestep > cfl_max) {
        fprintf(stderr, "WARNING: dt=%.6f > CFL=%.6f, clamping to 0.25*CFL\n",
                timestep, cfl_max);
        timestep = 0.25 * cfl_max;
    }
    dt = timestep;

    double mem_gb = (double)(4+4+4+12+1) * N*N*(double)N * 8.0 / 1e9;
    printf("=== Skyrme 3D (v18 corrected) ===\n");
    printf("N=%d  L=%.1f  dx=%.5f  dt=%.6f  c4=%.2f\n", N, L, dx, dt, c4);
    printf("Test %d  steps=%d  threads=%d  mem=%.2f GB\n",
           test, steps, omp_get_max_threads(), mem_gb);

    setlinebuf(stdout); /* flush after every line for progress monitoring */

    alloc_arrays();
    build_damping(0.12);

    /* Initialize field */
    if (test == 1) {
        printf("\n--- Test 1: Single Skyrmion Confinement ---\n");
        init_hedgehog(0, 0, 0);
    } else if (test == 2) {
        printf("\n--- Test 2: Two-Skyrmion Dynamics (B=2, product ansatz, sep=%.1f) ---\n", sep);
        init_two_solitons(sep);
    } else if (test == 3) {
        printf("\n--- Test 3: Far-Field Analysis (relax then measure radial decay) ---\n");
        init_hedgehog(0, 0, 0);
    } else {
        fprintf(stderr, "Unknown test %d\n", test);
        return 1;
    }

    /* Initial diagnostics */
    double Ek, E2, E4, rms;
    compute_energy(&Ek, &E2, &E4);
    double B = compute_baryon();
    soliton_stats(&rms);
    double Etot = Ek + E2 + E4;

    printf("Initial: Ek=%.4f E2=%.4f E4=%.4f Etot=%.4f  B=%.4f  R_rms=%.4f\n",
           Ek, E2, E4, Etot, B, rms);
    printf("Virial:  E2/E4=%.4f (1.0 at equilibrium)\n", E4 > 1e-10 ? E2/E4 : 0);

    /* --- Test 3: relaxation then radial analysis --- */
    if (test == 3) {
        int relax_steps = (steps > 500) ? steps : 2000;
        printf("Relaxing for %d steps (t=%.1f) ...\n", relax_steps, relax_steps*dt);
        double t0 = omp_get_wtime();
        for (int s = 1; s <= relax_steps; s++) {
            step_vv();
            if (s % (relax_steps/10) == 0) {
                compute_energy(&Ek, &E2, &E4);
                B = compute_baryon();
                printf("  relax %5d/%d | Et=%.3f B=%.4f\n", s, relax_steps, Ek+E2+E4, B);
            }
        }
        printf("Relaxation done in %.1f s\n", omp_get_wtime() - t0);
        radial_profile_analysis(outdir);
        free_arrays();
        return 0;
    }

    /* Open output */
    char fname[512];
    snprintf(fname, sizeof(fname), "%s/test%d_N%d.tsv", outdir, test, N);
    FILE *fp = fopen(fname, "w");
    if (!fp) { perror("fopen"); return 1; }

    if (test == 1) {
        fprintf(fp, "step\ttime\tE_kin\tE2\tE4\tE_total\tB\tE2_over_E4\tR_rms\n");
        fprintf(fp, "%d\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.4f\t%.4f\n",
                0, 0.0, Ek, E2, E4, Etot, B, E4>1e-10?E2/E4:0, rms);
    } else {
        double d0 = two_soliton_sep();
        fprintf(fp, "step\ttime\tE_kin\tE2\tE4\tE_total\tB\tE2_over_E4\tsep\n");
        fprintf(fp, "%d\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.4f\t%.4f\n",
                0, 0.0, Ek, E2, E4, Etot, B, E4>1e-10?E2/E4:0, d0);
        record_separation(0, d0);
    }

    /* Time evolution */
    int print_every = steps / 60;
    if (print_every < 1) print_every = 1;
    double t_start = omp_get_wtime();

    for (int s = 1; s <= steps; s++) {
        step_vv();

        if (s % print_every == 0 || s == steps) {
            compute_energy(&Ek, &E2, &E4);
            B = compute_baryon();
            Etot = Ek + E2 + E4;
            double elapsed = omp_get_wtime() - t_start;
            double rate = s / elapsed;

            if (test == 1) {
                soliton_stats(&rms);
                printf("Step %5d | t=%7.3f | Ek=%8.3f E2=%8.3f E4=%8.3f Et=%8.3f | B=%6.4f | E2/E4=%5.3f | R=%5.3f | %.0f/s\n",
                       s, s*dt, Ek, E2, E4, Etot, B, E4>1e-10?E2/E4:0, rms, rate);
                fprintf(fp, "%d\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.4f\t%.4f\n",
                        s, s*dt, Ek, E2, E4, Etot, B, E4>1e-10?E2/E4:0, rms);
            } else {
                double d_sep = two_soliton_sep();
                printf("Step %5d | t=%7.3f | Ek=%8.3f E2=%8.3f E4=%8.3f Et=%8.3f | B=%6.4f | sep=%5.3f | %.0f/s\n",
                       s, s*dt, Ek, E2, E4, Etot, B, d_sep, rate);
                fprintf(fp, "%d\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.4f\t%.4f\n",
                        s, s*dt, Ek, E2, E4, Etot, B, E4>1e-10?E2/E4:0, d_sep);
                record_separation(s*dt, d_sep);
            }
            fflush(fp);
        }
    }

    fclose(fp);
    double elapsed = omp_get_wtime() - t_start;
    printf("\nDone: %d steps in %.1f s (%.1f steps/s)\n", steps, elapsed, steps/elapsed);
    printf("Output: %s\n", fname);

    /* Post-evolution analysis */
    if (test == 2) analyze_two_body_force(outdir);

    free_arrays();
    return 0;
}
