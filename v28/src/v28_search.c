/*  v28_search.c — Automated braided soliton search
 *
 *  Targets: l=2>20%, |Phi_T|>0, winding conserved, fc>0.3, stable
 *  Method:  LHS screening → CMA-ES refinement → full validation
 *
 *  Build:  gcc -O3 -fopenmp -o v28_search v28_search.c -lm
 *  Run:    ./v28_search [-t1pop 1024] [-t2gen 5] [-threads 16] [-o data/]
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include <omp.h>
#include <time.h>

/* ================================================================
   Section 1: Constants and Configuration
   ================================================================ */

#define NDIM    16
#define NFIELDS 3
#define PI      3.14159265358979323846
#define MAX_N   128

/* Default search parameters (overridable via CLI) */
static int T1_POP   = 1024;   /* Tier 1 population */
static int T1_N     = 64;     /* Tier 1 grid */
static int T1_TIME  = 100;    /* Tier 1 sim time */

static int T2_NPOP  = 4;      /* CMA-ES populations */
static int T2_NGEN  = 5;      /* generations per pop */
static int T2_LAM   = 16;     /* candidates per gen */
static int T2_N     = 80;     /* Tier 2 grid */
static int T2_TIME  = 200;    /* Tier 2 sim time */

static int T3_TOP   = 10;     /* Tier 3 candidates */
static int T3_N     = 128;    /* Tier 3 grid */
static int T3_TIME  = 500;    /* Tier 3 sim time */

static double DOMAIN_L = 20.0;
static int NTHREADS = 16;
static char OUTDIR[256] = "data";

/* Parameter bounds: [lo, hi] in physical space */
/*  0: A1        1: A2        2: A3        3: delta2    4: delta3
    5: R_tube    6: ellip     7: ell_ang   8: k_fac     9: A_bg
   10: R_disp   11: ell_rot   12: mu       13: kappa    14: mass
   15: lam_pw
   R_disp: per-field strand displacement (0=all centered, >0=displaced 120deg apart)
   ell_rot: per-field ellipse rotation (0=same, 1=each rotated 2*pi*a/3) */
static const double LO[NDIM] = {
    0.1, 0.1, 0.1,   0.0, 0.0,
    1.5, 0.0, 0.0,   0.5, 0.0,
    0.0, 0.0,       -100.0, 5.0,  0.0,
    0.0
};
static const double HI[NDIM] = {
    2.0, 2.0, 2.0,   2*PI, 2*PI,
    6.0, 0.8, PI,    4.0,  0.5,
    3.0, 1.0,        -5.0, 100.0, 1.5,
    1.0
};
static const char *PNAME[NDIM] = {
    "A1","A2","A3","delta2","delta3",
    "R_tube","ellip","ell_ang","k_fac","A_bg",
    "R_disp","ell_rot","mu","kappa","mass",
    "lam_pw"
};

/* ================================================================
   Section 2: Types
   ================================================================ */

typedef struct {
    double params[NDIM];   /* normalized [0,1] */
    double phys[NDIM];     /* physical values */
    double fitness;
    double l2_frac;        /* quadrupole metric */
    double torsion_flux;   /* |Phi_T| */
    double fc;             /* field concentration */
    double peak_P;         /* max |triple product| */
    double winding_init;
    double winding_final;
    double energy;
    double transverse_l2;  /* xy-plane anisotropy */
    int stable;
} Result;

typedef struct {
    double *phi[NFIELDS];
    double *vel[NFIELDS];
    double *acc[NFIELDS];
    int N;
    double L, dx, dt;
} Grid;

typedef struct {
    double mean[NDIM];
    double sigma;
    double var[NDIM];      /* diagonal covariance */
    double ps[NDIM];       /* sigma evolution path */
    double pc[NDIM];       /* covariance path */
} CMA;

/* ================================================================
   Section 3: RNG (thread-safe splitmix64 + xoshiro256**)
   ================================================================ */

static __thread unsigned long long rng_s[4];

static unsigned long long splitmix64(unsigned long long *state) {
    unsigned long long z = (*state += 0x9e3779b97f4a7c15ULL);
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}

static void rng_seed(unsigned long long seed) {
    unsigned long long s = seed;
    for (int i = 0; i < 4; i++) rng_s[i] = splitmix64(&s);
}

static inline unsigned long long rotl(unsigned long long x, int k) {
    return (x << k) | (x >> (64 - k));
}

static unsigned long long rng_next(void) {
    unsigned long long *s = rng_s;
    unsigned long long result = rotl(s[1] * 5, 7) * 9;
    unsigned long long t = s[1] << 17;
    s[2] ^= s[0]; s[3] ^= s[1]; s[1] ^= s[2]; s[0] ^= s[3];
    s[2] ^= t; s[3] = rotl(s[3], 45);
    return result;
}

static double rng_uniform(void) {
    return (rng_next() >> 11) * 0x1.0p-53;
}

static double rng_normal(void) {
    /* Box-Muller */
    double u1 = rng_uniform(), u2 = rng_uniform();
    if (u1 < 1e-30) u1 = 1e-30;
    return sqrt(-2.0 * log(u1)) * cos(2.0 * PI * u2);
}

/* ================================================================
   Section 4: Parameter encoding / decoding
   ================================================================ */

static void decode_params(const double *norm, double *phys) {
    for (int d = 0; d < NDIM; d++)
        phys[d] = LO[d] + norm[d] * (HI[d] - LO[d]);
}

static void encode_params(const double *phys, double *norm) {
    for (int d = 0; d < NDIM; d++) {
        double range = HI[d] - LO[d];
        norm[d] = (fabs(range) > 1e-30) ? (phys[d] - LO[d]) / range : 0.5;
    }
}

static double clamp01(double x) { return x < 0 ? 0 : (x > 1 ? 1 : x); }

/* ================================================================
   Section 5: Grid memory
   ================================================================ */

static Grid *alloc_grid(double L) {
    Grid *g = calloc(1, sizeof(Grid));
    int N3 = MAX_N * MAX_N * MAX_N;
    for (int a = 0; a < NFIELDS; a++) {
        g->phi[a] = calloc(N3, sizeof(double));
        g->vel[a] = calloc(N3, sizeof(double));
        g->acc[a] = calloc(N3, sizeof(double));
    }
    g->L = L;
    return g;
}

static void free_grid(Grid *g) {
    for (int a = 0; a < NFIELDS; a++) {
        free(g->phi[a]); free(g->vel[a]); free(g->acc[a]);
    }
    free(g);
}

static void setup_grid(Grid *g, int N) {
    g->N = N;
    g->dx = 2.0 * g->L / (N - 1);
    g->dt = 0.20 * g->dx;
}

/* ================================================================
   Section 6: Initialization
   ================================================================ */

static void init_braid(Grid *g, const double *phys) {
    double A[3]     = {phys[0], phys[1], phys[2]};
    double delta[3] = {0.0, phys[3], phys[4]};
    double R_tube   = phys[5];
    double ellip    = phys[6];
    double ell_ang  = phys[7];
    double k_fac    = phys[8];
    double A_bg     = phys[9];
    double R_disp   = phys[10];
    double ell_rot  = phys[11];
    double mass     = phys[14];

    int N = g->N, NN = N * N;
    double dx = g->dx, L = g->L;
    double k = k_fac * PI / L;
    double omega = sqrt(k * k + mass * mass);
    double inv_2R2 = 1.0 / (2.0 * R_tube * R_tube);

    /* Per-field centers (displaced by R_disp at 120deg intervals) */
    double cx[3], cy[3];
    for (int a = 0; a < 3; a++) {
        double ang = 2.0 * PI * a / 3.0;
        cx[a] = R_disp * cos(ang);
        cy[a] = R_disp * sin(ang);
    }

    /* Per-field ellipse orientations */
    double ell_angles[3];
    for (int a = 0; a < 3; a++) {
        if (ell_rot > 0.5)
            ell_angles[a] = ell_ang + 2.0 * PI * a / 3.0;
        else
            ell_angles[a] = ell_ang;
    }

    double sx = 1.0 + ellip, sy = 1.0 - ellip;

    /* Zero arrays */
    int N3 = N * N * N;
    for (int a = 0; a < NFIELDS; a++) {
        memset(g->phi[a], 0, N3 * sizeof(double));
        memset(g->vel[a], 0, N3 * sizeof(double));
        memset(g->acc[a], 0, N3 * sizeof(double));
    }

    for (int i = 0; i < N; i++) {
        double x = -L + i * dx;
        for (int j = 0; j < N; j++) {
            double y = -L + j * dx;

            for (int kk = 0; kk < N; kk++) {
                double z = -L + kk * dx;
                int idx = i * NN + j * N + kk;

                for (int a = 0; a < NFIELDS; a++) {
                    /* Per-field displaced coordinates */
                    double xc = x - cx[a];
                    double yc = y - cy[a];

                    /* Per-field elliptical envelope */
                    double ca = cos(ell_angles[a]), sa = sin(ell_angles[a]);
                    double xr = xc * ca + yc * sa;
                    double yr = -xc * sa + yc * ca;
                    double r2e = xr * xr / (sx * sx) + yr * yr / (sy * sy);
                    double env = exp(-r2e * inv_2R2);

                    double ph = k * z + delta[a];
                    double amp = A[a] * env + A_bg;

                    g->phi[a][idx] = amp * cos(ph);
                    g->vel[a][idx] = omega * amp * sin(ph);
                }
            }
        }
    }
}

/* ================================================================
   Section 7: Force computation
   ================================================================ */

static void compute_forces(Grid *g, const double *phys, int parallel) {
    int N = g->N, NN = N * N, N3 = N * N * N;
    double idx2 = 1.0 / (g->dx * g->dx);
    double mu    = phys[12];
    double kappa = phys[13];
    double mass2 = phys[14] * phys[14];
    double lpw   = phys[15];

    #pragma omp parallel for if(parallel) schedule(static)
    for (int idx = 0; idx < N3; idx++) {
        int i = idx / NN;
        int j = (idx / N) % N;
        int k = idx % N;

        /* Skip boundary (margin=1 in x,y) */
        if (i < 1 || i >= N - 1 || j < 1 || j >= N - 1) {
            g->acc[0][idx] = g->acc[1][idx] = g->acc[2][idx] = 0;
            continue;
        }

        /* z-periodic neighbors */
        int kp = (k + 1) % N;
        int km = (k - 1 + N) % N;
        int idx_kp = i * NN + j * N + kp;
        int idx_km = i * NN + j * N + km;

        double p0 = g->phi[0][idx];
        double p1 = g->phi[1][idx];
        double p2 = g->phi[2][idx];
        double P  = p0 * p1 * p2;
        double denom = 1.0 + kappa * P * P;
        double mu_P_d2 = mu * P / (denom * denom);

        for (int a = 0; a < NFIELDS; a++) {
            /* Laplacian (7-point stencil) */
            double lap = (g->phi[a][idx + NN] + g->phi[a][idx - NN]
                        + g->phi[a][idx + N]  + g->phi[a][idx - N]
                        + g->phi[a][idx_kp]   + g->phi[a][idx_km]
                        - 6.0 * g->phi[a][idx]) * idx2;

            /* Triple product force: -dV/dphi_a */
            double dPda;
            if (a == 0)      dPda = p1 * p2;
            else if (a == 1) dPda = p0 * p2;
            else             dPda = p0 * p1;
            double f_triple = mu_P_d2 * dPda;

            /* Pairwise force */
            double f_pw = lpw * (g->phi[(a+1)%3][idx] + g->phi[(a+2)%3][idx]);

            g->acc[a][idx] = lap - mass2 * g->phi[a][idx] - f_triple - f_pw;
        }
    }
}

/* ================================================================
   Section 8: Time stepping
   ================================================================ */

static void verlet_step(Grid *g) {
    int N3 = g->N * g->N * g->N;
    double dt = g->dt;
    double hdt = 0.5 * dt;

    for (int a = 0; a < NFIELDS; a++) {
        double *p = g->phi[a], *v = g->vel[a], *ac = g->acc[a];
        for (int idx = 0; idx < N3; idx++) {
            v[idx] += hdt * ac[idx];
            p[idx] += dt  * v[idx];
        }
    }
}

static void verlet_finish(Grid *g) {
    int N3 = g->N * g->N * g->N;
    double hdt = 0.5 * g->dt;
    for (int a = 0; a < NFIELDS; a++) {
        double *v = g->vel[a], *ac = g->acc[a];
        for (int idx = 0; idx < N3; idx++)
            v[idx] += hdt * ac[idx];
    }
}

static void apply_damping(Grid *g) {
    int N = g->N, NN = N * N;
    double dx = g->dx, L = g->L;
    double r_start = 0.70 * L;
    double r_end   = 0.95 * L;
    double inv_dr  = 1.0 / (r_end - r_start + 1e-30);

    for (int i = 0; i < N; i++) {
        double x = -L + i * dx;
        for (int j = 0; j < N; j++) {
            double y = -L + j * dx;
            double rp = sqrt(x * x + y * y);
            if (rp <= r_start) continue;

            double f = (rp - r_start) * inv_dr;
            if (f > 1.0) f = 1.0;
            double damp = 1.0 - 0.98 * f * f;

            for (int kk = 0; kk < N; kk++) {
                int idx = i * NN + j * N + kk;
                for (int a = 0; a < NFIELDS; a++) {
                    g->phi[a][idx] *= damp;
                    g->vel[a][idx] *= damp;
                }
            }
        }
    }
}

/* ================================================================
   Section 9: Diagnostics
   ================================================================ */

static int check_blowup(Grid *g) {
    int N3 = g->N * g->N * g->N;
    for (int a = 0; a < NFIELDS; a++) {
        for (int idx = 0; idx < N3; idx += 37) { /* sparse check */
            if (fabs(g->phi[a][idx]) > 50.0) return 1;
        }
    }
    return 0;
}

static double compute_winding(Grid *g) {
    int N = g->N, NN = N * N;
    int ic = N / 2, jc = N / 2;
    double total = 0;

    for (int k = 0; k < N - 1; k++) {
        int idx0 = ic * NN + jc * N + k;
        int idx1 = ic * NN + jc * N + k + 1;
        double re0 = g->phi[0][idx0], im0 = g->phi[1][idx0];
        double re1 = g->phi[0][idx1], im1 = g->phi[1][idx1];
        double dth = atan2(im1 * re0 - re1 * im0, re1 * re0 + im1 * im0);
        total += dth;
    }
    /* Periodic wrap: last to first */
    {
        int idx0 = ic * NN + jc * N + (N - 1);
        int idx1 = ic * NN + jc * N + 0;
        double re0 = g->phi[0][idx0], im0 = g->phi[1][idx0];
        double re1 = g->phi[0][idx1], im1 = g->phi[1][idx1];
        double dth = atan2(im1 * re0 - re1 * im0, re1 * re0 + im1 * im0);
        total += dth;
    }
    return total / (2.0 * PI);
}

static void compute_diagnostics(Grid *g, const double *phys, Result *res) {
    int N = g->N, NN = N * N, N3 = N * N * N;
    double dx = g->dx, L = g->L, dV = dx * dx * dx;
    double mu = phys[12], kappa = phys[13], mass2 = phys[14] * phys[14];
    double lpw = phys[15];
    double R_core = 8.0;
    if (R_core > L / 3.0) R_core = L / 3.0;
    double Rc2 = R_core * R_core;

    double phi2_total = 0, phi2_core = 0;
    double E_total = 0;
    double peak_P = 0;

    /* Quadrupole tensor components */
    double Ixx = 0, Iyy = 0, Izz = 0, Ixy = 0, Ixz = 0, Iyz = 0;
    double M0 = 0, R2sum = 0;

    for (int i = 1; i < N - 1; i++) {
        double x = -L + i * dx;
        for (int j = 1; j < N - 1; j++) {
            double y = -L + j * dx;
            double rp2 = x * x + y * y;
            for (int k = 0; k < N; k++) {
                double z = -L + k * dx;
                int idx = i * NN + j * N + k;

                double p0 = g->phi[0][idx];
                double p1 = g->phi[1][idx];
                double p2 = g->phi[2][idx];
                double rho = p0 * p0 + p1 * p1 + p2 * p2;

                phi2_total += rho;
                if (rp2 < Rc2) phi2_core += rho;

                double P = p0 * p1 * p2;
                double absP = fabs(P);
                if (absP > peak_P) peak_P = absP;

                /* Energy (kinetic + gradient + mass + potential + pairwise) */
                double ek = 0, eg = 0;
                for (int a = 0; a < NFIELDS; a++) {
                    ek += 0.5 * g->vel[a][idx] * g->vel[a][idx];
                    /* Gradient via central diff (skip boundaries) */
                    int kp = (k + 1) % N, km = (k - 1 + N) % N;
                    double gx = (g->phi[a][idx + NN] - g->phi[a][idx - NN]) / (2 * dx);
                    double gy = (g->phi[a][idx + N]  - g->phi[a][idx - N])  / (2 * dx);
                    double gz = (g->phi[a][i*NN+j*N+kp] - g->phi[a][i*NN+j*N+km]) / (2*dx);
                    eg += 0.5 * (gx * gx + gy * gy + gz * gz);
                }
                double em = 0.5 * mass2 * rho;
                double ep = (mu / 2.0) * P * P / (1.0 + kappa * P * P);
                double epw = lpw * (p0*p1 + p1*p2 + p2*p0);
                E_total += (ek + eg + em + ep + epw) * dV;

                /* Quadrupole moments */
                double r2 = rp2 + z * z;
                M0    += rho;
                R2sum += rho * r2;
                Ixx   += rho * x * x;
                Iyy   += rho * y * y;
                Izz   += rho * z * z;
                Ixy   += rho * x * y;
                Ixz   += rho * x * z;
                Iyz   += rho * y * z;
            }
        }
    }

    /* Traceless quadrupole */
    double Qxx = 3*Ixx - R2sum, Qyy = 3*Iyy - R2sum, Qzz = 3*Izz - R2sum;
    double Qxy = 3*Ixy, Qxz = 3*Ixz, Qyz = 3*Iyz;
    double Q2 = Qxx*Qxx + Qyy*Qyy + Qzz*Qzz
              + 2*(Qxy*Qxy + Qxz*Qxz + Qyz*Qyz);
    double R2_mean = R2sum / (M0 + 1e-30);
    res->l2_frac = sqrt(Q2) / (M0 * R2_mean + 1e-30);

    /* Transverse (xy) anisotropy specifically */
    res->transverse_l2 = fabs(Ixx - Iyy) / (Ixx + Iyy + 1e-30);

    /* Field concentration */
    res->fc = phi2_core / (phi2_total + 1e-30);

    /* Torsion flux: integrate ALL three torsion components, take max
       omega_xy = 0.5*(d_x phi_1 - d_y phi_0)  [flux through z-slice]
       omega_xz = 0.5*(d_x phi_2 - d_z phi_0)  [flux through y-slice]
       omega_yz = 0.5*(d_y phi_2 - d_z phi_1)  [flux through x-slice]
       Also compute total torsion energy for localization check */
    double Phi_xy = 0, Phi_xz = 0, Phi_yz = 0;
    double torsion_energy = 0;
    {
        int k_mid = N / 2, j_mid = N / 2, i_mid = N / 2;
        for (int i = 2; i < N - 2; i++) {
            for (int j = 2; j < N - 2; j++) {
                /* xy torsion at z=z_mid */
                {
                    int idx = i * NN + j * N + k_mid;
                    double dx_p1 = (g->phi[1][idx+NN] - g->phi[1][idx-NN]) / (2*dx);
                    double dy_p0 = (g->phi[0][idx+N]  - g->phi[0][idx-N])  / (2*dx);
                    Phi_xy += 0.5 * (dx_p1 - dy_p0) * dx * dx;
                }
                /* xz torsion at y=y_mid */
                {
                    int idx = i * NN + j_mid * N + j; /* reuse j as k-index */
                    int kp = (j + 1) % N, km = (j - 1 + N) % N;
                    int idx_kp = i * NN + j_mid * N + kp;
                    int idx_km = i * NN + j_mid * N + km;
                    double dx_p2 = (g->phi[2][idx+NN] - g->phi[2][idx-NN]) / (2*dx);
                    double dz_p0 = (g->phi[0][idx_kp] - g->phi[0][idx_km]) / (2*dx);
                    Phi_xz += 0.5 * (dx_p2 - dz_p0) * dx * dx;
                }
                /* yz torsion at x=x_mid */
                {
                    int idx = i_mid * NN + i * N + j; /* reuse i as j, j as k */
                    int kp = (j + 1) % N, km = (j - 1 + N) % N;
                    int idx_kp = i_mid * NN + i * N + kp;
                    int idx_km = i_mid * NN + i * N + km;
                    double dy_p2 = (g->phi[2][idx+N]  - g->phi[2][idx-N])  / (2*dx);
                    double dz_p1 = (g->phi[1][idx_kp] - g->phi[1][idx_km]) / (2*dx);
                    Phi_yz += 0.5 * (dy_p2 - dz_p1) * dx * dx;
                }
            }
        }
        /* Also: total torsion energy in core */
        for (int i = 2; i < N - 2; i++) {
            double x = -L + i * dx;
            for (int j = 2; j < N - 2; j++) {
                double y = -L + j * dx;
                if (x*x + y*y > Rc2) continue;
                for (int k = 0; k < N; k++) {
                    int idx = i * NN + j * N + k;
                    int kp = (k+1)%N, km = (k-1+N)%N;
                    int ikp = i*NN+j*N+kp, ikm = i*NN+j*N+km;
                    /* All 3 torsion components */
                    double dxp1 = (g->phi[1][idx+NN]-g->phi[1][idx-NN])/(2*dx);
                    double dyp0 = (g->phi[0][idx+N] -g->phi[0][idx-N]) /(2*dx);
                    double dxp2 = (g->phi[2][idx+NN]-g->phi[2][idx-NN])/(2*dx);
                    double dzp0 = (g->phi[0][ikp]   -g->phi[0][ikm])   /(2*dx);
                    double dyp2 = (g->phi[2][idx+N] -g->phi[2][idx-N]) /(2*dx);
                    double dzp1 = (g->phi[1][ikp]   -g->phi[1][ikm])   /(2*dx);
                    double w01 = 0.5*(dxp1 - dyp0);
                    double w02 = 0.5*(dxp2 - dzp0);
                    double w12 = 0.5*(dyp2 - dzp1);
                    torsion_energy += (w01*w01 + w02*w02 + w12*w12) * dV;
                }
            }
        }
    }
    res->torsion_flux = fmax(fabs(Phi_xy), fmax(fabs(Phi_xz), fabs(Phi_yz)));

    res->peak_P = peak_P;
    res->energy = E_total;
    res->winding_final = compute_winding(g);
}

/* ================================================================
   Section 10: Fitness function
   ================================================================ */

static void compute_fitness(Result *r) {
    if (!r->stable) { r->fitness = 0; return; }
    if (r->fc < 0.02) { r->fitness = 0; return; }

    /* Transverse l2 is what matters for spin-2 gravity (xy-plane asymmetry).
       Total l2 is dominated by trivial cylindrical elongation. */
    double f = 0;
    f += 4.0 * fmin(r->transverse_l2 / 0.15, 1.0); /* trans_l2: 0→0.15 → 0→4 */
    f += 3.0 * fmin(r->torsion_flux / 0.3, 1.0);   /* Phi_T: 0→0.3 → 0→3 */
    f += 2.0 * fmin(r->fc / 0.5, 1.0);              /* fc: 0→0.5 → 0→2 */
    f += 1.0 * fmin(r->peak_P / 2.0, 1.0);          /* |P|: 0→2 → 0→1 */
    f += 0.5 * fmin(r->l2_frac / 1.0, 1.0);         /* total l2: mild bonus */

    /* Winding conservation bonus */
    double dw = fabs(r->winding_final - r->winding_init);
    if (fabs(round(r->winding_init) - r->winding_init) < 0.2) {
        f += 1.5 * exp(-5.0 * dw);  /* bonus for integer winding conserved */
    }

    r->fitness = f;
}

/* ================================================================
   Section 11: Evaluation wrapper
   ================================================================ */

static void evaluate(Grid *g, const double *norm_params, int N_grid,
                     double T_sim, Result *res, int parallel_grid) {
    /* Decode */
    double phys[NDIM];
    decode_params(norm_params, phys);
    memcpy(res->params, norm_params, NDIM * sizeof(double));
    memcpy(res->phys, phys, NDIM * sizeof(double));

    /* Setup and init */
    setup_grid(g, N_grid);
    init_braid(g, phys);

    res->winding_init = compute_winding(g);
    res->stable = 1;

    /* Time evolution */
    int Nsteps = (int)(T_sim / g->dt);
    int check_interval = Nsteps / 4;
    if (check_interval < 1) check_interval = 1;

    for (int step = 0; step < Nsteps; step++) {
        compute_forces(g, phys, parallel_grid);
        verlet_step(g);
        compute_forces(g, phys, parallel_grid);
        verlet_finish(g);
        apply_damping(g);

        /* Early termination */
        if (step > 0 && (step % check_interval == 0)) {
            if (check_blowup(g)) {
                res->stable = 0;
                res->fitness = 0;
                return;
            }
        }
    }

    compute_diagnostics(g, phys, res);
    compute_fitness(res);
}

/* ================================================================
   Section 12: LHS generation + anchors
   ================================================================ */

/* Fisher-Yates shuffle */
static void shuffle_int(int *arr, int n) {
    for (int i = n - 1; i > 0; i--) {
        int j = (int)(rng_uniform() * (i + 1));
        if (j > i) j = i;
        int tmp = arr[i]; arr[i] = arr[j]; arr[j] = tmp;
    }
}

static void latin_hypercube(double (*samples)[NDIM], int n) {
    int *perm = malloc(n * sizeof(int));
    for (int d = 0; d < NDIM; d++) {
        for (int i = 0; i < n; i++) perm[i] = i;
        shuffle_int(perm, n);
        for (int i = 0; i < n; i++) {
            double u = (perm[i] + rng_uniform()) / n;  /* uniform in [0,1] */
            samples[i][d] = u;
        }
    }
    free(perm);
}

static void set_anchors(double (*samples)[NDIM]) {
    /* Anchor configs in physical space, then encode to [0,1] */
    double anchors[][NDIM] = {
        /* 0: V27 optimal (m=0, 120deg, no displacement) */
        {0.8, 0.8, 0.8, 2*PI/3, 4*PI/3, 3.0, 0.0, 0.0,
         1.0, 0.0, 0.0, 0.0, -50, 50, 0.0, 0.0},
        /* 1: V27 regime 1 (m=1, 120deg) */
        {0.8, 0.8, 0.8, 2*PI/3, 4*PI/3, 3.0, 0.0, 0.0,
         1.0, 0.0, 0.0, 0.0, -20, 20, 1.0, 0.0},
        /* 2: 180deg + strand displacement (KEY: breaks torsion symmetry) */
        {0.8, 0.8, 0.8, 0.0, PI, 3.0, 0.0, 0.0,
         1.0, 0.0, 1.5, 0.0, -50, 50, 0.0, 0.0},
        /* 3: 120deg + displacement + per-field ellipse rotation */
        {0.8, 0.8, 0.8, 2*PI/3, 4*PI/3, 3.0, 0.4, 0.0,
         1.0, 0.0, 1.0, 1.0, -50, 50, 0.0, 0.0},
        /* 4: Asymmetric amplitudes + displacement */
        {1.2, 0.6, 0.3, 2*PI/3, 4*PI/3, 3.0, 0.0, 0.0,
         1.0, 0.0, 1.5, 0.0, -50, 50, 0.0, 0.0},
        /* 5: Elliptical + rotated per field (max torsion setup) */
        {0.8, 0.8, 0.8, 2*PI/3, 4*PI/3, 3.0, 0.6, 0.0,
         1.0, 0.0, 0.0, 1.0, -50, 50, 0.0, 0.0},
        /* 6: Large displacement + m=1 (true 3-strand braid) */
        {0.8, 0.8, 0.8, 2*PI/3, 4*PI/3, 2.0, 0.0, 0.0,
         1.0, 0.0, 2.5, 0.0, -50, 50, 1.0, 0.0},
        /* 7: 180deg + elliptical + rotated + mass */
        {1.0, 1.0, 0.8, 0.0, PI, 3.0, 0.5, 0.0,
         1.0, 0.0, 1.0, 1.0, -30, 30, 0.7, 0.3},
    };
    int n_anchors = 8;

    for (int i = 0; i < n_anchors; i++) {
        double norm[NDIM];
        encode_params(anchors[i], norm);
        for (int d = 0; d < NDIM; d++)
            samples[i][d] = clamp01(norm[d]);
    }
}

/* ================================================================
   Section 13: CMA-ES (sep-CMA-ES, diagonal covariance)
   ================================================================ */

static void cma_init(CMA *es, const double *seed_mean, double sigma0) {
    memcpy(es->mean, seed_mean, NDIM * sizeof(double));
    es->sigma = sigma0;
    for (int d = 0; d < NDIM; d++) {
        es->var[d] = 1.0;
        es->ps[d]  = 0.0;
        es->pc[d]  = 0.0;
    }
}

static void cma_sample(CMA *es, double (*pop)[NDIM], int lam) {
    for (int i = 0; i < lam; i++) {
        for (int d = 0; d < NDIM; d++) {
            double z = rng_normal();
            pop[i][d] = es->mean[d] + es->sigma * sqrt(es->var[d]) * z;
            pop[i][d] = clamp01(pop[i][d]);
        }
    }
}

static int cmp_result_desc(const void *a, const void *b) {
    double fa = ((const Result *)a)->fitness;
    double fb = ((const Result *)b)->fitness;
    return (fb > fa) - (fb < fa);
}

static void cma_update(CMA *es, Result *pop_sorted, int lam) {
    int mu = lam / 2;
    double n = NDIM;

    /* Weights */
    double w[mu], wsum = 0;
    for (int i = 0; i < mu; i++) {
        w[i] = log(mu + 0.5) - log(i + 1.0);
        wsum += w[i];
    }
    for (int i = 0; i < mu; i++) w[i] /= wsum;

    double mu_eff = 0;
    for (int i = 0; i < mu; i++) mu_eff += w[i] * w[i];
    mu_eff = 1.0 / mu_eff;

    /* Learning rates */
    double cs  = (mu_eff + 2) / (n + mu_eff + 5);
    double ds  = 1 + 2 * fmax(0, sqrt((mu_eff - 1) / (n + 1)) - 1) + cs;
    double cc  = (4 + mu_eff / n) / (n + 4 + 2 * mu_eff / n);
    double c1  = 2.0 / ((n + 1.3) * (n + 1.3) + mu_eff);
    double cmu_val = fmin(1 - c1,
                     2 * (mu_eff - 2 + 1.0 / mu_eff)
                     / ((n + 2) * (n + 2) + mu_eff));
    double chiN = sqrt(n) * (1 - 1.0 / (4 * n) + 1.0 / (21 * n * n));

    /* New mean */
    double m_old[NDIM], m_new[NDIM];
    memcpy(m_old, es->mean, sizeof(m_old));
    for (int d = 0; d < NDIM; d++) {
        m_new[d] = 0;
        for (int i = 0; i < mu; i++)
            m_new[d] += w[i] * pop_sorted[i].params[d];
    }
    memcpy(es->mean, m_new, sizeof(m_new));

    /* Step-size path */
    double ps_norm = 0;
    for (int d = 0; d < NDIM; d++) {
        double inv_sd = 1.0 / (sqrt(es->var[d]) + 1e-30);
        es->ps[d] = (1 - cs) * es->ps[d]
                   + sqrt(cs * (2 - cs) * mu_eff)
                   * (m_new[d] - m_old[d]) / (es->sigma * sqrt(es->var[d]) + 1e-30);
        ps_norm += es->ps[d] * es->ps[d];
    }
    ps_norm = sqrt(ps_norm);

    /* Sigma update */
    es->sigma *= exp(cs / ds * (ps_norm / chiN - 1));
    if (es->sigma < 1e-10) es->sigma = 1e-10;
    if (es->sigma > 2.0)   es->sigma = 2.0;

    /* Covariance path */
    for (int d = 0; d < NDIM; d++) {
        es->pc[d] = (1 - cc) * es->pc[d]
                   + sqrt(cc * (2 - cc) * mu_eff)
                   * (m_new[d] - m_old[d]) / es->sigma;
    }

    /* Diagonal covariance update */
    for (int d = 0; d < NDIM; d++) {
        double rank_mu_sum = 0;
        for (int i = 0; i < mu; i++) {
            double y = (pop_sorted[i].params[d] - m_old[d]) / es->sigma;
            rank_mu_sum += w[i] * y * y;
        }
        es->var[d] = (1 - c1 - cmu_val) * es->var[d]
                    + c1 * es->pc[d] * es->pc[d]
                    + cmu_val * rank_mu_sum;
        if (es->var[d] < 1e-20) es->var[d] = 1e-20;
        if (es->var[d] > 10.0)  es->var[d] = 10.0;
    }
}

/* ================================================================
   Section 14: I/O
   ================================================================ */

static void write_header(FILE *fp) {
    for (int d = 0; d < NDIM; d++) fprintf(fp, "%s\t", PNAME[d]);
    fprintf(fp, "fitness\tl2_frac\ttrans_l2\ttorsion\tfc\tpeak_P\t"
                "wind_init\twind_final\tenergy\tstable\n");
}

static void write_result(FILE *fp, const Result *r) {
    for (int d = 0; d < NDIM; d++) fprintf(fp, "%.4f\t", r->phys[d]);
    fprintf(fp, "%.4f\t%.6f\t%.6f\t%.6f\t%.4f\t%.4f\t%.3f\t%.3f\t%.1f\t%d\n",
            r->fitness, r->l2_frac, r->transverse_l2, r->torsion_flux,
            r->fc, r->peak_P, r->winding_init, r->winding_final,
            r->energy, r->stable);
}

static void print_top(Result *res, int n, int show) {
    qsort(res, n, sizeof(Result), cmp_result_desc);
    printf("\n  Rank  Fitness   l2      trans   torFlux fc      peakP   wind    mass    mu      delta2  delta3  ellip\n");
    printf("  ----  -------   ------  ------  ------- ------  ------  ------  ------  ------  ------  ------  ------\n");
    int m = show < n ? show : n;
    for (int i = 0; i < m; i++) {
        Result *r = &res[i];
        printf("  %3d   %7.3f   %.4f  %.4f  %.4f  %.4f  %.4f  %+.3f  %.2f  %7.1f  %.2f  %.2f  %.2f\n",
               i+1, r->fitness, r->l2_frac, r->transverse_l2, r->torsion_flux,
               r->fc, r->peak_P, r->winding_final,
               r->phys[14], r->phys[12], r->phys[3], r->phys[4], r->phys[6]);
    }
}

/* ================================================================
   Section 15: Correlation analysis
   ================================================================ */

static void analyze_correlations(Result *res, int n) {
    printf("\n=== Parameter-Metric Correlations (Pearson r) ===\n");
    printf("%-10s %8s %8s %8s %8s %8s %8s\n",
           "Param", "fitness", "l2", "trans_l2", "torsion", "fc", "peakP");
    printf("%-10s %8s %8s %8s %8s %8s %8s\n",
           "-----", "-------", "------", "--------", "-------", "------", "------");

    /* Only use stable results */
    int ns = 0;
    for (int i = 0; i < n; i++) if (res[i].stable) ns++;
    if (ns < 10) { printf("(too few stable results for correlation)\n"); return; }

    for (int d = 0; d < NDIM; d++) {
        /* Compute Pearson correlation with each metric */
        double metrics[6];  /* fitness, l2, trans, torsion, fc, peakP */
        double sx = 0, sx2 = 0;
        double sy[6] = {0}, sxy[6] = {0}, sy2[6] = {0};

        for (int i = 0; i < n; i++) {
            if (!res[i].stable) continue;
            double x = res[i].phys[d];
            sx += x; sx2 += x * x;
            double ys[6] = {res[i].fitness, res[i].l2_frac, res[i].transverse_l2,
                           res[i].torsion_flux, res[i].fc, res[i].peak_P};
            for (int m = 0; m < 6; m++) {
                sy[m]  += ys[m];
                sxy[m] += x * ys[m];
                sy2[m] += ys[m] * ys[m];
            }
        }

        printf("%-10s", PNAME[d]);
        for (int m = 0; m < 6; m++) {
            double num = ns * sxy[m] - sx * sy[m];
            double den = sqrt((ns * sx2 - sx * sx) * (ns * sy2[m] - sy[m] * sy[m]));
            double r = (den > 1e-30) ? num / den : 0;
            printf(" %+7.3f ", r);
        }
        printf("\n");
    }
}

/* ================================================================
   Section 16: Main
   ================================================================ */

static void parse_args(int argc, char **argv) {
    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "-t1pop")   && i+1 < argc) T1_POP   = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-t2gen")   && i+1 < argc) T2_NGEN  = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-t2lam")   && i+1 < argc) T2_LAM   = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-t3top")   && i+1 < argc) T3_TOP   = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-threads") && i+1 < argc) NTHREADS = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-o")       && i+1 < argc) strncpy(OUTDIR, argv[++i], 255);
        else if (!strcmp(argv[i], "-L")       && i+1 < argc) DOMAIN_L = atof(argv[++i]);
        else if (!strcmp(argv[i], "-t1n")     && i+1 < argc) T1_N     = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-t2n")     && i+1 < argc) T2_N     = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-t3n")     && i+1 < argc) T3_N     = atoi(argv[++i]);
        else {
            printf("Usage: %s [-t1pop N] [-t2gen N] [-t2lam N] [-t3top N]\n"
                   "          [-threads N] [-o dir] [-L domain]\n"
                   "          [-t1n N] [-t2n N] [-t3n N]\n", argv[0]);
            exit(0);
        }
    }
}

int main(int argc, char **argv) {
    parse_args(argc, argv);
    omp_set_num_threads(NTHREADS);

    printf("V28 Braided Soliton Search\n");
    printf("Threads: %d, Domain: L=%.0f\n", NTHREADS, DOMAIN_L);
    printf("Tier1: pop=%d N=%d T=%d\n", T1_POP, T1_N, T1_TIME);
    printf("Tier2: %d pops x %d gen x %d lam, N=%d T=%d\n",
           T2_NPOP, T2_NGEN, T2_LAM, T2_N, T2_TIME);
    printf("Tier3: top %d, N=%d T=%d\n\n", T3_TOP, T3_N, T3_TIME);

    /* Allocate per-thread grids */
    Grid **grids = malloc(NTHREADS * sizeof(Grid *));
    for (int t = 0; t < NTHREADS; t++)
        grids[t] = alloc_grid(DOMAIN_L);

    /* ============================================================
       TIER 1: LHS Screening
       ============================================================ */
    printf("=== TIER 1: Screening %d candidates at N=%d, T=%d ===\n", T1_POP, T1_N, T1_TIME);
    fflush(stdout);
    double t_start = omp_get_wtime();

    double (*t1_params)[NDIM] = malloc(T1_POP * sizeof(*t1_params));
    Result *t1_results = calloc(T1_POP, sizeof(Result));

    /* Generate LHS samples, then overwrite first 8 with anchors */
    rng_seed(42);
    latin_hypercube(t1_params, T1_POP);
    set_anchors(t1_params);

    int t1_done = 0;
    #pragma omp parallel for schedule(dynamic, 1) num_threads(NTHREADS)
    for (int c = 0; c < T1_POP; c++) {
        int tid = omp_get_thread_num();
        /* Each thread needs its own RNG state */
        if (tid == 0 && c == 0) rng_seed(42);
        else rng_seed(42 + c * 137);  /* deterministic per-candidate seed */

        evaluate(grids[tid], t1_params[c], T1_N, T1_TIME, &t1_results[c], 0);

        #pragma omp atomic
        t1_done++;

        if (t1_done % (T1_POP / 8) == 0) {
            #pragma omp critical
            {
                printf("  Tier1: %d/%d done (%.0fs)\n", t1_done, T1_POP,
                       omp_get_wtime() - t_start);
                fflush(stdout);
            }
        }
    }

    double t1_time = omp_get_wtime() - t_start;
    int t1_stable = 0;
    for (int i = 0; i < T1_POP; i++) if (t1_results[i].stable) t1_stable++;
    printf("  Tier1 complete: %.1fs, %d/%d stable\n", t1_time, t1_stable, T1_POP);

    /* Write Tier 1 results */
    {
        char path[512];
        snprintf(path, sizeof(path), "%s/tier1_results.tsv", OUTDIR);
        FILE *fp = fopen(path, "w");
        if (fp) {
            write_header(fp);
            for (int i = 0; i < T1_POP; i++) write_result(fp, &t1_results[i]);
            fclose(fp);
            printf("  Written: %s\n", path);
        }
    }

    print_top(t1_results, T1_POP, 20);

    /* ============================================================
       TIER 2: CMA-ES Refinement
       ============================================================ */
    printf("\n=== TIER 2: CMA-ES refinement (%d pops x %d gen x %d) at N=%d ===\n",
           T2_NPOP, T2_NGEN, T2_LAM, T2_N);
    fflush(stdout);
    double t2_start = omp_get_wtime();

    /* Sort Tier 1 by fitness, take top T2_NPOP as seeds */
    qsort(t1_results, T1_POP, sizeof(Result), cmp_result_desc);

    int t2_total = T2_NPOP * T2_NGEN * T2_LAM;
    Result *t2_results = calloc(t2_total, sizeof(Result));
    int t2_idx = 0;

    for (int pop = 0; pop < T2_NPOP; pop++) {
        CMA es;
        cma_init(&es, t1_results[pop].params, 0.3);
        printf("  Pop %d: seeded from Tier1 rank %d (fitness=%.3f)\n",
               pop, pop+1, t1_results[pop].fitness);

        for (int gen = 0; gen < T2_NGEN; gen++) {
            double (*pop_params)[NDIM] = malloc(T2_LAM * sizeof(*pop_params));
            Result *pop_results = calloc(T2_LAM, sizeof(Result));

            cma_sample(&es, pop_params, T2_LAM);

            #pragma omp parallel for schedule(dynamic, 1) num_threads(NTHREADS)
            for (int c = 0; c < T2_LAM; c++) {
                int tid = omp_get_thread_num();
                rng_seed(1000 + pop * 10000 + gen * 100 + c);
                evaluate(grids[tid], pop_params[c], T2_N, T2_TIME,
                         &pop_results[c], 0);
            }

            /* Sort and update CMA-ES */
            qsort(pop_results, T2_LAM, sizeof(Result), cmp_result_desc);
            cma_update(&es, pop_results, T2_LAM);

            printf("    Pop %d Gen %d: best=%.3f, sigma=%.4f\n",
                   pop, gen, pop_results[0].fitness, es.sigma);
            fflush(stdout);

            /* Store results */
            for (int c = 0; c < T2_LAM && t2_idx < t2_total; c++)
                t2_results[t2_idx++] = pop_results[c];

            free(pop_params);
            free(pop_results);
        }
    }

    double t2_time = omp_get_wtime() - t2_start;
    printf("  Tier2 complete: %.1fs, %d evaluations\n", t2_time, t2_idx);

    /* Write Tier 2 results */
    {
        char path[512];
        snprintf(path, sizeof(path), "%s/tier2_results.tsv", OUTDIR);
        FILE *fp = fopen(path, "w");
        if (fp) {
            write_header(fp);
            for (int i = 0; i < t2_idx; i++) write_result(fp, &t2_results[i]);
            fclose(fp);
            printf("  Written: %s\n", path);
        }
    }

    print_top(t2_results, t2_idx, 20);

    /* ============================================================
       TIER 3: Validation (full resolution)
       ============================================================ */
    printf("\n=== TIER 3: Validating top %d at N=%d, T=%d ===\n",
           T3_TOP, T3_N, T3_TIME);
    fflush(stdout);
    double t3_start = omp_get_wtime();

    /* Merge Tier 1 and Tier 2 results, take top T3_TOP */
    int n_all = T1_POP + t2_idx;
    Result *all_results = malloc(n_all * sizeof(Result));
    memcpy(all_results, t1_results, T1_POP * sizeof(Result));
    memcpy(all_results + T1_POP, t2_results, t2_idx * sizeof(Result));
    qsort(all_results, n_all, sizeof(Result), cmp_result_desc);

    /* Deduplicate: skip candidates with very similar params */
    Result *t3_candidates = calloc(T3_TOP, sizeof(Result));
    int t3_count = 0;
    for (int i = 0; i < n_all && t3_count < T3_TOP; i++) {
        int dup = 0;
        for (int j = 0; j < t3_count; j++) {
            double dist2 = 0;
            for (int d = 0; d < NDIM; d++) {
                double diff = all_results[i].params[d] - t3_candidates[j].params[d];
                dist2 += diff * diff;
            }
            if (dist2 < 0.01) { dup = 1; break; }
        }
        if (!dup) t3_candidates[t3_count++] = all_results[i];
    }

    Result *t3_results = calloc(t3_count, sizeof(Result));

    /* Sequential with grid-parallel OMP */
    for (int c = 0; c < t3_count; c++) {
        printf("  Tier3 %d/%d: evaluating...", c + 1, t3_count);
        fflush(stdout);
        evaluate(grids[0], t3_candidates[c].params, T3_N, T3_TIME,
                 &t3_results[c], 1);
        printf(" fitness=%.3f l2=%.4f tor=%.4f fc=%.4f\n",
               t3_results[c].fitness, t3_results[c].l2_frac,
               t3_results[c].torsion_flux, t3_results[c].fc);
        fflush(stdout);
    }

    double t3_time = omp_get_wtime() - t3_start;
    printf("  Tier3 complete: %.1fs\n", t3_time);

    /* Write Tier 3 results */
    {
        char path[512];
        snprintf(path, sizeof(path), "%s/tier3_results.tsv", OUTDIR);
        FILE *fp = fopen(path, "w");
        if (fp) {
            write_header(fp);
            for (int i = 0; i < t3_count; i++) write_result(fp, &t3_results[i]);
            fclose(fp);
            printf("  Written: %s\n", path);
        }
    }

    printf("\n=== FINAL RESULTS (Tier 3, N=%d, T=%d) ===\n", T3_N, T3_TIME);
    print_top(t3_results, t3_count, t3_count);

    /* ============================================================
       Analysis
       ============================================================ */
    printf("\n=== ANALYSIS ===\n");

    /* Correlations from all stable Tier 1+2 results */
    analyze_correlations(all_results, n_all);

    /* Per-metric best */
    printf("\n--- Best per metric (from all tiers) ---\n");
    Result *best_l2 = &all_results[0], *best_tor = &all_results[0];
    Result *best_fc = &all_results[0], *best_P = &all_results[0];
    for (int i = 0; i < n_all; i++) {
        if (!all_results[i].stable) continue;
        if (all_results[i].l2_frac > best_l2->l2_frac) best_l2 = &all_results[i];
        if (all_results[i].torsion_flux > best_tor->torsion_flux) best_tor = &all_results[i];
        if (all_results[i].fc > best_fc->fc) best_fc = &all_results[i];
        if (all_results[i].peak_P > best_P->peak_P) best_P = &all_results[i];
    }

    printf("Best l2=%.4f:   mu=%.1f kap=%.1f m=%.2f d2=%.2f d3=%.2f el=%.2f A=%.1f/%.1f/%.1f\n",
           best_l2->l2_frac, best_l2->phys[12], best_l2->phys[13], best_l2->phys[14],
           best_l2->phys[3], best_l2->phys[4], best_l2->phys[6],
           best_l2->phys[0], best_l2->phys[1], best_l2->phys[2]);
    printf("Best tor=%.4f:  mu=%.1f kap=%.1f m=%.2f d2=%.2f d3=%.2f el=%.2f A=%.1f/%.1f/%.1f\n",
           best_tor->torsion_flux, best_tor->phys[12], best_tor->phys[13], best_tor->phys[14],
           best_tor->phys[3], best_tor->phys[4], best_tor->phys[6],
           best_tor->phys[0], best_tor->phys[1], best_tor->phys[2]);
    printf("Best fc=%.4f:   mu=%.1f kap=%.1f m=%.2f d2=%.2f d3=%.2f el=%.2f A=%.1f/%.1f/%.1f\n",
           best_fc->fc, best_fc->phys[12], best_fc->phys[13], best_fc->phys[14],
           best_fc->phys[3], best_fc->phys[4], best_fc->phys[6],
           best_fc->phys[0], best_fc->phys[1], best_fc->phys[2]);
    printf("Best |P|=%.4f: mu=%.1f kap=%.1f m=%.2f d2=%.2f d3=%.2f el=%.2f A=%.1f/%.1f/%.1f\n",
           best_P->peak_P, best_P->phys[12], best_P->phys[13], best_P->phys[14],
           best_P->phys[3], best_P->phys[4], best_P->phys[6],
           best_P->phys[0], best_P->phys[1], best_P->phys[2]);

    /* Total time */
    double total = omp_get_wtime() - t_start;
    printf("\nTotal time: %.1f min (Tier1: %.1fs, Tier2: %.1fs, Tier3: %.1fs)\n",
           total / 60, t1_time, t2_time, t3_time);

    /* Cleanup */
    free(t1_params); free(t1_results);
    free(t2_results); free(all_results);
    free(t3_candidates); free(t3_results);
    for (int t = 0; t < NTHREADS; t++) free_grid(grids[t]);
    free(grids);

    return 0;
}
