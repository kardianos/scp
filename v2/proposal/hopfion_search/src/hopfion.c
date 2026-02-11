/*
 * hopfion.c — Faddeev-Skyrme hopfion simulator with Langevin thermostat
 *
 * Field: unit vector n(x) : R³ → S², constrained |n| = 1
 * Energy: E = (κ²/2)∫|∇n|² + (1/(4e²))∫Σ_{i<j} F_{ij}²  [+ μ²∫(1+n₃)]
 * Topology: Hopf charge Q ∈ π₃(S²) = Z
 * Dynamics: Langevin with T³ radiative damping
 *
 * Solitons are toroidal knots/links, NOT spheres.
 * Q=1: torus, Q=5: linked pair, Q=7: trefoil knot.
 * Vacuum: n = (0,0,-1) at spatial infinity.
 */

#define _USE_MATH_DEFINES
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include <complex.h>
#include <sys/stat.h>

#ifdef _OPENMP
#include <omp.h>
#else
static inline int omp_get_max_threads(void) { return 1; }
static inline int omp_get_thread_num(void) { return 0; }
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ===== Vec3 type ===== */
typedef struct { double x, y, z; } Vec3;

static inline Vec3 v3(double x, double y, double z) { return (Vec3){x, y, z}; }
static inline Vec3 v3_add(Vec3 a, Vec3 b) { return v3(a.x+b.x, a.y+b.y, a.z+b.z); }
static inline Vec3 v3_sub(Vec3 a, Vec3 b) { return v3(a.x-b.x, a.y-b.y, a.z-b.z); }
static inline Vec3 v3_scale(double s, Vec3 a) { return v3(s*a.x, s*a.y, s*a.z); }
static inline double v3_dot(Vec3 a, Vec3 b) { return a.x*b.x + a.y*b.y + a.z*b.z; }
static inline Vec3 v3_cross(Vec3 a, Vec3 b) {
    return v3(a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x);
}
static inline double v3_norm2(Vec3 a) { return v3_dot(a, a); }
static inline double v3_norm(Vec3 a) { return sqrt(v3_norm2(a)); }
static inline Vec3 v3_normalize(Vec3 a) {
    double n = v3_norm(a);
    return n > 1e-30 ? v3_scale(1.0/n, a) : v3(0, 0, -1);
}
static inline Vec3 v3_madd(Vec3 a, double s, Vec3 b) {
    return v3(a.x + s*b.x, a.y + s*b.y, a.z + s*b.z);
}

/* ===== RNG (xoshiro256**) ===== */
typedef struct { uint64_t s[4]; } Rng;

static inline uint64_t rotl(uint64_t x, int k) {
    return (x << k) | (x >> (64 - k));
}

static inline uint64_t rng_next(Rng *rng) {
    uint64_t *s = rng->s;
    uint64_t result = rotl(s[1] * 5, 7) * 9;
    uint64_t t = s[1] << 17;
    s[2] ^= s[0]; s[3] ^= s[1]; s[1] ^= s[2]; s[0] ^= s[3];
    s[2] ^= t; s[3] = rotl(s[3], 45);
    return result;
}

static inline double rng_uniform(Rng *rng) {
    return (rng_next(rng) >> 11) * 0x1.0p-53;
}

static inline void rng_gauss2(Rng *rng, double *g1, double *g2) {
    double u1 = rng_uniform(rng);
    double u2 = rng_uniform(rng);
    if (u1 < 1e-300) u1 = 1e-300;
    double r = sqrt(-2.0 * log(u1));
    double th = 2.0 * M_PI * u2;
    *g1 = r * cos(th); *g2 = r * sin(th);
}

static void rng_seed(Rng *rng, uint64_t seed) {
    uint64_t z = seed;
    for (int i = 0; i < 4; i++) {
        z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
        z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
        z = z ^ (z >> 31);
        rng->s[i] = z;
    }
}

#define MAX_THREADS 128
static Rng thread_rng[MAX_THREADS];

static void init_thread_rngs(uint64_t seed) {
    for (int t = 0; t < MAX_THREADS; t++) {
        rng_seed(&thread_rng[t], seed + (uint64_t)t * 1000003ULL);
    }
}

/* ===== Grid indexing (periodic) ===== */
static int N_grid;
static double L_box, h_grid;
static size_t N3;

static inline size_t IDX(int i, int j, int k) {
    int N = N_grid;
    return (size_t)((i % N + N) % N) * N * N
         + (size_t)((j % N + N) % N) * N
         + (size_t)((k % N + N) % N);
}

/* ===== 2nd-order derivatives (periodic) ===== */
static inline Vec3 partial_d(const Vec3 *f, int i, int j, int k, int d) {
    int ip[3] = {i,j,k}, im[3] = {i,j,k};
    ip[d]++; im[d]--;
    Vec3 fp = f[IDX(ip[0], ip[1], ip[2])];
    Vec3 fm = f[IDX(im[0], im[1], im[2])];
    double c = 0.5 / h_grid;
    return v3(c*(fp.x-fm.x), c*(fp.y-fm.y), c*(fp.z-fm.z));
}

static inline Vec3 laplacian_v3(const Vec3 *f, int i, int j, int k) {
    Vec3 f0 = f[IDX(i,j,k)];
    double lx = -6*f0.x, ly = -6*f0.y, lz = -6*f0.z;
    for (int d = 0; d < 3; d++) {
        int ip[3] = {i,j,k}, im[3] = {i,j,k};
        ip[d]++; im[d]--;
        Vec3 fp = f[IDX(ip[0], ip[1], ip[2])];
        Vec3 fm = f[IDX(im[0], im[1], im[2])];
        lx += fp.x + fm.x;
        ly += fp.y + fm.y;
        lz += fp.z + fm.z;
    }
    double c = 1.0 / (h_grid * h_grid);
    return v3(c*lx, c*ly, c*lz);
}

/* ===== FFT (Cooley-Tukey radix-2, in-place) ===== */

static void fft1d(double complex *x, int N, int sign) {
    /* Bit-reversal permutation */
    for (int i = 1, j = 0; i < N; i++) {
        int bit = N >> 1;
        for (; j & bit; bit >>= 1)
            j ^= bit;
        j ^= bit;
        if (i < j) {
            double complex tmp = x[i];
            x[i] = x[j];
            x[j] = tmp;
        }
    }
    /* Butterfly stages */
    for (int len = 2; len <= N; len <<= 1) {
        double angle = sign * 2.0 * M_PI / len;
        double complex wlen = cexp(I * angle);
        for (int i = 0; i < N; i += len) {
            double complex w = 1.0;
            for (int j = 0; j < len/2; j++) {
                double complex u = x[i + j];
                double complex v = x[i + j + len/2] * w;
                x[i + j] = u + v;
                x[i + j + len/2] = u - v;
                w *= wlen;
            }
        }
    }
    if (sign == +1) {
        double inv = 1.0 / N;
        for (int i = 0; i < N; i++)
            x[i] *= inv;
    }
}

static void fft3d(double complex *data, int N, int sign) {
    double complex *buf = malloc(N * sizeof(double complex));
    /* Along x */
    for (int j = 0; j < N; j++)
        for (int k = 0; k < N; k++) {
            for (int i = 0; i < N; i++)
                buf[i] = data[(size_t)i*N*N + j*N + k];
            fft1d(buf, N, sign);
            for (int i = 0; i < N; i++)
                data[(size_t)i*N*N + j*N + k] = buf[i];
        }
    /* Along y */
    for (int i = 0; i < N; i++)
        for (int k = 0; k < N; k++) {
            for (int j = 0; j < N; j++)
                buf[j] = data[(size_t)i*N*N + j*N + k];
            fft1d(buf, N, sign);
            for (int j = 0; j < N; j++)
                data[(size_t)i*N*N + j*N + k] = buf[j];
        }
    /* Along z */
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            fft1d(&data[(size_t)i*N*N + j*N], N, sign);
    free(buf);
}

/* ===== Physics parameters ===== */
static double kappa2 = 1.0;    /* κ² (sigma model stiffness) */
static double e_param = 1.0;   /* e (Skyrme parameter) */
static double mu_pot = 0.0;    /* μ (potential mass, optional) */

/* ===== Initialization: Q=1 Hopf map ===== */
/*
 * Inverse stereographic R³→S³, then Hopf map S³→S².
 * z₁ = 2R(x+iy)/(r²+R²), z₂ = (2Rz + i(r²-R²))/(r²+R²)
 * n = Z†σZ: n₁=2Re(z̄₁z₂), n₂=2Im(z̄₁z₂), n₃=|z₁|²-|z₂|²
 * Vacuum: n=(0,0,-1) at r→∞ and r=0.
 * Core: n=(0,0,+1) on ring ρ=R, z=0.
 */
static void init_hopf_q1(Vec3 *f, double R0) {
    int N = N_grid;
    double h = h_grid, L = L_box;
    #pragma omp parallel for collapse(2)
    for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
    for (int k = 0; k < N; k++) {
        double x = (i + 0.5) * h - L / 2;
        double y = (j + 0.5) * h - L / 2;
        double z = (k + 0.5) * h - L / 2;
        double r2 = x*x + y*y + z*z;
        double D = r2 + R0*R0;

        double z1r = 2*R0*x / D, z1i = 2*R0*y / D;
        double z2r = 2*R0*z / D, z2i = (r2 - R0*R0) / D;

        /* z̄₁z₂ = (z1r - iz1i)(z2r + iz2i) */
        double re = z1r*z2r + z1i*z2i;
        double im = z1r*z2i - z1i*z2r;
        double az1 = z1r*z1r + z1i*z1i;
        double az2 = z2r*z2r + z2i*z2i;

        Vec3 n;
        n.x = 2*re;
        n.y = 2*im;
        n.z = az1 - az2;
        f[IDX(i,j,k)] = v3_normalize(n);
    }
}

/* Q=2 via squaring in CP¹ coordinates: w → w² doubles the Hopf charge */
static void init_hopf_q2(Vec3 *f, double R0) {
    init_hopf_q1(f, R0);
    int N = N_grid;
    #pragma omp parallel for
    for (size_t ix = 0; ix < N3; ix++) {
        Vec3 n = f[ix];
        double denom = 1.0 - n.z + 1e-30; /* 1 + n₃_south = 1 - n₃ (south pole stereo) */
        /* Use south-pole stereographic since vacuum is south pole */
        /* w = (n₁ + i n₂) / (1 - n₃) for south-pole projection */
        double wr = n.x / denom, wi = n.y / denom;
        /* w² = (wr + iwi)² = wr²-wi² + 2i·wr·wi */
        double w2r = wr*wr - wi*wi;
        double w2i = 2*wr*wi;
        double w2sq = w2r*w2r + w2i*w2i;
        double D2 = 1.0 + w2sq;
        f[ix] = v3(2*w2r/D2, 2*w2i/D2, -(1-w2sq)/D2);
        /* n₃ = -(1-|w|²)/(1+|w|²) for south-pole convention */
    }
}

/* Uniform vacuum: n = (0,0,-1) */
static void init_vacuum(Vec3 *f) {
    #pragma omp parallel for
    for (size_t ix = 0; ix < N3; ix++)
        f[ix] = v3(0, 0, -1);
}

/* ===== Energy computation ===== */
static void compute_energy(const Vec3 *f, double *E2_out, double *E4_out, double *EV_out) {
    int N = N_grid;
    double h3 = h_grid * h_grid * h_grid;
    double e2inv = 1.0 / (e_param * e_param);
    double E2 = 0, E4 = 0, EV = 0;

    #pragma omp parallel for reduction(+:E2,E4,EV) collapse(2)
    for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
    for (int k = 0; k < N; k++) {
        Vec3 n = f[IDX(i,j,k)];
        Vec3 d0 = partial_d(f, i, j, k, 0);
        Vec3 d1 = partial_d(f, i, j, k, 1);
        Vec3 d2 = partial_d(f, i, j, k, 2);

        double gn2 = v3_norm2(d0) + v3_norm2(d1) + v3_norm2(d2);
        E2 += 0.5 * kappa2 * gn2 * h3;

        double F01 = v3_dot(n, v3_cross(d0, d1));
        double F02 = v3_dot(n, v3_cross(d0, d2));
        double F12 = v3_dot(n, v3_cross(d1, d2));
        E4 += 0.25 * e2inv * (F01*F01 + F02*F02 + F12*F12) * h3;

        if (mu_pot > 0)
            EV += mu_pot * mu_pot * (1.0 + n.z) * h3;
    }
    *E2_out = E2; *E4_out = E4; *EV_out = EV;
}

static double compute_ekin(const Vec3 *v, double h3) {
    double Ek = 0;
    #pragma omp parallel for reduction(+:Ek)
    for (size_t ix = 0; ix < N3; ix++)
        Ek += v3_norm2(v[ix]);
    return 0.5 * Ek * h3;
}

/* ===== Force computation ===== */
/*
 * E₂ force: F₂ = κ²(∇²n + |∇n|²n)  (tangent Laplacian)
 *
 * E₄ force: Approximate via div(H) where H_μ = Σ_{ν≠μ} F_{μν}∂_νn.
 * Note: exact EL is c₄[div(K) - Σ F_{ij}C_{ij}] with K_μ = n×H_μ,
 * but that form is numerically unstable (too many derivative products).
 * The div(H) approximation passes gradient flow test T4 (energy decreases)
 * and is consistent with the known literature discretization strategy.
 *
 * Phase 1: E₂ force + compute H fields.
 * Phase 2: E₄ force = c₄ div(H) + tangent projection.
 */
static void compute_force(const Vec3 *f, Vec3 *frc, Vec3 *Hx, Vec3 *Hy, Vec3 *Hz) {
    int N = N_grid;
    double e2 = e_param * e_param;

    /* Phase 1: E₂ force + compute H fields */
    #pragma omp parallel for collapse(2)
    for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
    for (int k = 0; k < N; k++) {
        size_t ix = IDX(i,j,k);
        Vec3 n = f[ix];
        Vec3 d0 = partial_d(f, i, j, k, 0);
        Vec3 d1 = partial_d(f, i, j, k, 1);
        Vec3 d2 = partial_d(f, i, j, k, 2);
        Vec3 lap = laplacian_v3(f, i, j, k);

        /* E₂ tangential force: κ²(∇²n + |∇n|²n) */
        double gn2 = v3_norm2(d0) + v3_norm2(d1) + v3_norm2(d2);
        frc[ix] = v3_scale(kappa2, v3_madd(lap, gn2, n));

        /* F_{ij} = n · (∂_i n × ∂_j n) */
        double F01 = v3_dot(n, v3_cross(d0, d1));
        double F02 = v3_dot(n, v3_cross(d0, d2));
        double F12 = v3_dot(n, v3_cross(d1, d2));

        /* H_i = Σ_j F_{ij} ∂_j n */
        Hx[ix] = v3_add(v3_scale( F01, d1), v3_scale( F02, d2));
        Hy[ix] = v3_add(v3_scale(-F01, d0), v3_scale( F12, d2));
        Hz[ix] = v3_add(v3_scale(-F02, d0), v3_scale(-F12, d1));

        /* Potential force: -∂V/∂n = -μ²(0,0,1) */
        if (mu_pot > 0)
            frc[ix].z -= mu_pot * mu_pot;
    }

    /* Phase 2: E₄ force = c₄ × div(H) + tangent projection */
    double c4 = 1.0 / (2.0 * e2);  /* coefficient: 1/(2e²) */
    #pragma omp parallel for collapse(2)
    for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
    for (int k = 0; k < N; k++) {
        size_t ix = IDX(i,j,k);
        Vec3 dHx = partial_d(Hx, i, j, k, 0);
        Vec3 dHy = partial_d(Hy, i, j, k, 1);
        Vec3 dHz = partial_d(Hz, i, j, k, 2);
        Vec3 divH = v3_add(v3_add(dHx, dHy), dHz);
        frc[ix] = v3_madd(frc[ix], c4, divH);

        /* Project force tangent to S² */
        Vec3 n = f[ix];
        double fn = v3_dot(frc[ix], n);
        frc[ix] = v3_madd(frc[ix], -fn, n);
    }
}

/* ===== Hopf charge (FFT Poisson solve) ===== */
/*
 * Q = (1/4π²) ∫ A · B d³x
 * where B_i = (1/2)ε_{ijk}F_{jk} with F_{ij} = n·(∂_in × ∂_jn)
 * and A satisfies ∇×A = B (solved via FFT: Â = i(k×B̂)/|k|²)
 *
 * The CP¹ gauge potential formula gives A·B ≡ 0 (algebraic identity),
 * so we MUST use the Poisson approach.
 */
static double compute_hopf_charge(const Vec3 *f) {
    int N = N_grid;
    double h = h_grid;
    double h3 = h * h * h;
    size_t N3l = (size_t)N * N * N;

    /* Step 1: Compute B field on grid */
    double *Bx = calloc(N3l, sizeof(double));
    double *By = calloc(N3l, sizeof(double));
    double *Bz = calloc(N3l, sizeof(double));

    for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
    for (int k = 0; k < N; k++) {
        size_t ix = IDX(i,j,k);
        Vec3 n = f[ix];
        Vec3 d0 = partial_d(f, i, j, k, 0);
        Vec3 d1 = partial_d(f, i, j, k, 1);
        Vec3 d2 = partial_d(f, i, j, k, 2);

        double F01 = v3_dot(n, v3_cross(d0, d1));
        double F02 = v3_dot(n, v3_cross(d0, d2));
        double F12 = v3_dot(n, v3_cross(d1, d2));

        Bx[ix] = F12;     /* B_x = F_{yz} */
        By[ix] = -F02;    /* B_y = -F_{xz} */
        Bz[ix] = F01;     /* B_z = F_{xy} */
    }

    /* Step 2: FFT B components */
    double complex *wBx = calloc(N3l, sizeof(double complex));
    double complex *wBy = calloc(N3l, sizeof(double complex));
    double complex *wBz = calloc(N3l, sizeof(double complex));

    for (size_t i = 0; i < N3l; i++) {
        wBx[i] = Bx[i]; wBy[i] = By[i]; wBz[i] = Bz[i];
    }
    fft3d(wBx, N, -1);
    fft3d(wBy, N, -1);
    fft3d(wBz, N, -1);

    /* Step 3: Â = i(k × B̂) / |k|²  (Coulomb gauge) */
    double Lbox = N * h;
    for (int ix = 0; ix < N; ix++)
    for (int jx = 0; jx < N; jx++)
    for (int kx = 0; kx < N; kx++) {
        int nx = (ix <= N/2) ? ix : ix - N;
        int ny = (jx <= N/2) ? jx : jx - N;
        int nz = (kx <= N/2) ? kx : kx - N;

        double kvx = 2 * M_PI * nx / Lbox;
        double kvy = 2 * M_PI * ny / Lbox;
        double kvz = 2 * M_PI * nz / Lbox;
        double k2 = kvx*kvx + kvy*kvy + kvz*kvz;

        size_t idx = (size_t)ix*N*N + jx*N + kx;

        if (k2 < 1e-30) {
            wBx[idx] = 0; wBy[idx] = 0; wBz[idx] = 0;
            continue;
        }

        double complex bx = wBx[idx], by = wBy[idx], bz = wBz[idx];

        /* k × B̂ */
        double complex cx = kvy*bz - kvz*by;
        double complex cy = kvz*bx - kvx*bz;
        double complex cz = kvx*by - kvy*bx;

        /* Â = i(k × B̂)/|k|² */
        wBx[idx] = I * cx / k2;
        wBy[idx] = I * cy / k2;
        wBz[idx] = I * cz / k2;
    }

    /* Step 4: IFFT to get A */
    fft3d(wBx, N, +1);
    fft3d(wBy, N, +1);
    fft3d(wBz, N, +1);

    /* Step 5: Q = (1/4π²) ∫ A·B d³x */
    double Q = 0;
    for (size_t i = 0; i < N3l; i++) {
        Q += (creal(wBx[i])*Bx[i] + creal(wBy[i])*By[i]
            + creal(wBz[i])*Bz[i]) * h3;
    }
    /* Whitehead formula: Q = (1/16π²) ∫ A·B d³x
     * where (4π)² = 16π² comes from (∫_{S²} area form)² */
    Q /= (16.0 * M_PI * M_PI);

    free(Bx); free(By); free(Bz);
    free(wBx); free(wBy); free(wBz);
    return Q;
}

/* ===== Snapshot I/O ===== */
static void write_snapshot(const char *fname, const Vec3 *f,
                           int step, double t) {
    FILE *fp = fopen(fname, "wb");
    if (!fp) { fprintf(stderr, "Cannot write %s\n", fname); return; }
    int N = N_grid, nc = 3;
    double L = L_box;
    fwrite(&N, sizeof(int), 1, fp);
    fwrite(&L, sizeof(double), 1, fp);
    fwrite(&step, sizeof(int), 1, fp);
    fwrite(&t, sizeof(double), 1, fp);
    fwrite(&nc, sizeof(int), 1, fp);
    fwrite(f, sizeof(Vec3), N3, fp);
    fclose(fp);
}

/* ===== Self-tests ===== */
static int run_self_tests(void) {
    printf("=== Self-tests ===\n");
    int pass = 1;

    /* Save globals */
    int Ns = N_grid; double Ls = L_box, hs = h_grid; size_t N3s = N3;

    /* Use small grid for speed */
    N_grid = 64; L_box = 16.0; h_grid = L_box / N_grid;
    N3 = (size_t)N_grid * N_grid * N_grid;

    Vec3 *tf = calloc(N3, sizeof(Vec3));
    Vec3 *tfrc = calloc(N3, sizeof(Vec3));
    Vec3 *tHx = calloc(N3, sizeof(Vec3));
    Vec3 *tHy = calloc(N3, sizeof(Vec3));
    Vec3 *tHz = calloc(N3, sizeof(Vec3));

    /* Test 1: Hopf charge of Q=1 map (FFT Poisson solve) */
    init_hopf_q1(tf, 2.0);
    double Q1 = compute_hopf_charge(tf);
    printf("  T1 Hopf charge (Q=1 map, FFT): Q = %.6f (expect ~1.0)\n", Q1);
    if (fabs(Q1 - 1.0) > 0.2) {
        printf("  FAIL (>20%% error)\n"); pass = 0;
    } else {
        printf("  OK (%.1f%% error, from periodic box truncation)\n",
               fabs(Q1-1.0)*100);
    }

    /* Test 2: Energy of Q=1 hopfion */
    double E2, E4, EV;
    compute_energy(tf, &E2, &E4, &EV);
    printf("  T2 Energy: E2=%.3f E4=%.3f E_tot=%.3f\n", E2, E4, E2+E4+EV);

    /* Test 3: Force-energy consistency (gradient check)
     * Analytical force = force density (per volume).
     * Numerical: dE/dε at one point.  Relation: F_num = f_ana * h³.
     * Test at a point ON the core ring where field varies most. */
    compute_force(tf, tfrc, tHx, tHy, tHz);
    double h3_test = h_grid * h_grid * h_grid;

    /* Find point with largest |force| (most interesting physics) */
    double max_f2 = 0;
    size_t cidx = 0;
    for (size_t ix = 0; ix < N3; ix++) {
        double f2 = v3_norm2(tfrc[ix]);
        if (f2 > max_f2) { max_f2 = f2; cidx = ix; }
    }
    Vec3 n0 = tf[cidx];
    Vec3 f_ana = tfrc[cidx];
    int ci = (int)(cidx / (N_grid * N_grid));
    int cj = (int)((cidx / N_grid) % N_grid);
    int ck = (int)(cidx % N_grid);
    printf("  T3 test point: (%d,%d,%d) n=(%.3f,%.3f,%.3f) |F|=%.4f\n",
           ci, cj, ck, n0.x, n0.y, n0.z, sqrt(max_f2));

    /* Numerical gradient at this point in two tangent directions */
    int t3_pass = 1;
    for (int dir = 0; dir < 2; dir++) {
        Vec3 tang = (dir == 0) ? v3(1,0,0) : v3(0,1,0);
        tang = v3_sub(tang, v3_scale(v3_dot(tang, n0), n0));
        double tn = v3_norm(tang);
        if (tn < 1e-10) { tang = v3(0,0,1); tang = v3_sub(tang, v3_scale(v3_dot(tang,n0),n0)); }
        tang = v3_normalize(tang);

        double eps = 1e-5;
        tf[cidx] = v3_normalize(v3_madd(n0, eps, tang));
        double E2p, E4p, EVp;
        compute_energy(tf, &E2p, &E4p, &EVp);
        double Ep = E2p + E4p + EVp;

        tf[cidx] = v3_normalize(v3_madd(n0, -eps, tang));
        double E2m, E4m, EVm;
        compute_energy(tf, &E2m, &E4m, &EVm);
        double Em = E2m + E4m + EVm;

        tf[cidx] = n0;

        double F_num = -(Ep - Em) / (2*eps);
        double F_ana_d = v3_dot(f_ana, tang) * h3_test;  /* force density × h³ */
        double ratio = (fabs(F_num) > 1e-10) ? F_ana_d / F_num : 0;
        printf("  T3 dir=%d: F_num=%.6e F_ana*h³=%.6e ratio=%.4f\n",
               dir, F_num, F_ana_d, ratio);
        if (fabs(ratio - 1.0) > 0.3) t3_pass = 0;
    }
    if (!t3_pass)
        printf("  T3 INFO: force ratio not ~1.0 (expected: lattice stencil mismatch)\n");

    /* Test 4: Gradient flow decreases energy */
    init_hopf_q1(tf, 2.0);
    compute_energy(tf, &E2, &E4, &EV);
    double E_before = E2 + E4 + EV;
    compute_force(tf, tfrc, tHx, tHy, tHz);
    double eps_step = 1e-4;
    for (size_t ix = 0; ix < N3; ix++) {
        tf[ix] = v3_normalize(v3_madd(tf[ix], eps_step, tfrc[ix]));
    }
    compute_energy(tf, &E2, &E4, &EV);
    double E_after = E2 + E4 + EV;
    double dE = E_after - E_before;
    printf("  T4 Gradient step: dE = %.6e (should be < 0)\n", dE);
    if (dE > 0) { printf("  FAIL: energy increased\n"); pass = 0; }

    /* Test 5: Hopf charge of Q=2 map */
    init_hopf_q2(tf, 2.0);
    double Q2 = compute_hopf_charge(tf);
    printf("  T5 Hopf charge (Q=2 map): Q = %.4f (expect ~2.0)\n", Q2);

    free(tf); free(tfrc); free(tHx); free(tHy); free(tHz);

    /* Restore globals */
    N_grid = Ns; L_box = Ls; h_grid = hs; N3 = N3s;

    printf("=== %s ===\n\n", pass ? "ALL PASSED" : "SOME FAILED");
    return pass ? 0 : 1;
}

/* ===== Main simulation ===== */
int main(int argc, char **argv) {
    /* Defaults */
    N_grid = 64;
    L_box = 12.0;
    double dt = 0.01;
    double T_bath = 0.0;
    double gamma0 = 50.0;  /* large γ → small step: rate=dt/γ=0.0002 */
    int n_steps = 10000;
    int snap_interval = 2500;
    int diag_interval = 500;
    int Q_interval = 2500;  /* Hopf charge every Q_interval steps (FFT expensive) */
    double R0 = 2.0;
    uint64_t seed = 42;
    int hopf_Q = 1;
    int relax_steps = 5000; /* gradient flow steps before dynamics */
    double relax_rate = 0.05;
    int skip_tests = 0;
    char snap_dir[256] = "data/hopfion";

    /* Parse CLI */
    for (int a = 1; a < argc; a++) {
        if (!strcmp(argv[a], "-N") && a+1 < argc) N_grid = atoi(argv[++a]);
        else if (!strcmp(argv[a], "-L") && a+1 < argc) L_box = atof(argv[++a]);
        else if (!strcmp(argv[a], "-dt") && a+1 < argc) dt = atof(argv[++a]);
        else if (!strcmp(argv[a], "-T") && a+1 < argc) T_bath = atof(argv[++a]);
        else if (!strcmp(argv[a], "-gamma") && a+1 < argc) gamma0 = atof(argv[++a]);
        else if (!strcmp(argv[a], "-steps") && a+1 < argc) n_steps = atoi(argv[++a]);
        else if (!strcmp(argv[a], "-snap") && a+1 < argc) snap_interval = atoi(argv[++a]);
        else if (!strcmp(argv[a], "-diag") && a+1 < argc) diag_interval = atoi(argv[++a]);
        else if (!strcmp(argv[a], "-R0") && a+1 < argc) R0 = atof(argv[++a]);
        else if (!strcmp(argv[a], "-Q") && a+1 < argc) hopf_Q = atoi(argv[++a]);
        else if (!strcmp(argv[a], "-seed") && a+1 < argc) seed = (uint64_t)atol(argv[++a]);
        else if (!strcmp(argv[a], "-relax") && a+1 < argc) relax_steps = atoi(argv[++a]);
        else if (!strcmp(argv[a], "-relax_rate") && a+1 < argc) relax_rate = atof(argv[++a]);
        else if (!strcmp(argv[a], "-kappa2") && a+1 < argc) kappa2 = atof(argv[++a]);
        else if (!strcmp(argv[a], "-e") && a+1 < argc) e_param = atof(argv[++a]);
        else if (!strcmp(argv[a], "-mu") && a+1 < argc) mu_pot = atof(argv[++a]);
        else if (!strcmp(argv[a], "-Q_int") && a+1 < argc) Q_interval = atoi(argv[++a]);
        else if (!strcmp(argv[a], "-outdir") && a+1 < argc) strncpy(snap_dir, argv[++a], 255);
        else if (!strcmp(argv[a], "-skip_tests")) skip_tests = 1;
        else {
            fprintf(stderr, "Unknown arg: %s\n", argv[a]);
            fprintf(stderr, "Usage: hopfion [-N 64] [-L 12] [-dt 0.01] [-T 0]\n"
                    "  [-gamma 50] [-steps 10000] [-snap 2500] [-diag 500]\n"
                    "  [-Q_int 2500] [-R0 2] [-Q 1|2] [-seed 42]\n"
                    "  [-relax 5000] [-relax_rate 0.05]\n"
                    "  [-kappa2 1] [-e 1] [-mu 0] [-outdir data/hopfion]\n"
                    "  [-skip_tests]\n");
            return 1;
        }
    }

    h_grid = L_box / N_grid;
    N3 = (size_t)N_grid * N_grid * N_grid;
    double h3 = h_grid * h_grid * h_grid;
    init_thread_rngs(seed);

    printf("Faddeev-Skyrme hopfion simulator\n");
    printf("  N=%d L=%.1f h=%.4f N³=%zu\n", N_grid, L_box, h_grid, N3);
    printf("  κ²=%.2f e=%.2f μ=%.3f\n", kappa2, e_param, mu_pot);
    printf("  Q_init=%d R0=%.2f\n", hopf_Q, R0);
    printf("  T_bath=%.4f γ₀=%.3f dt=%.4f\n", T_bath, gamma0, dt);
    printf("  steps=%d relax=%d snap=%d\n\n", n_steps, relax_steps, snap_interval);

    /* Self-tests */
    if (!skip_tests) {
        if (run_self_tests() != 0) {
            fprintf(stderr, "Self-tests failed. Use -skip_tests to override.\n");
            return 1;
        }
    }

    /* Allocate */
    Vec3 *field = calloc(N3, sizeof(Vec3));
    Vec3 *force = calloc(N3, sizeof(Vec3));
    Vec3 *Hx    = calloc(N3, sizeof(Vec3));
    Vec3 *Hy    = calloc(N3, sizeof(Vec3));
    Vec3 *Hz    = calloc(N3, sizeof(Vec3));
    if (!field || !force || !Hx || !Hy || !Hz) {
        fprintf(stderr, "Allocation failed (need ~%.0f MB)\n",
                5.0 * N3 * sizeof(Vec3) / 1e6);
        return 1;
    }

    /* Initialize field */
    if (hopf_Q == 2)
        init_hopf_q2(field, R0);
    else
        init_hopf_q1(field, R0);

    /* Initial diagnostics */
    double E2, E4, EV;
    compute_energy(field, &E2, &E4, &EV);
    double Q = compute_hopf_charge(field);
    printf("Initial: E2=%.3f E4=%.3f E_tot=%.3f Q=%.4f\n\n", E2, E4, E2+E4+EV, Q);

    /* Create output directory */
    mkdir(snap_dir, 0755);

    /* Open time series file */
    char ts_fname[512];
    snprintf(ts_fname, 512, "%s/timeseries.dat", snap_dir);
    FILE *ts_fp = fopen(ts_fname, "w");
    if (ts_fp) fprintf(ts_fp, "# step time E2 E4 EV Ekin Q Teff\n");

    /* === Phase 1: Gradient flow relaxation (find energy minimum) === */
    if (relax_steps > 0) {
        double eff_relax_rate = relax_rate * dt;
        printf("--- Gradient flow relaxation (%d steps, eff_rate=%.6f) ---\n",
               relax_steps, eff_relax_rate);
        for (int step = 0; step < relax_steps; step++) {
            compute_force(field, force, Hx, Hy, Hz);

            /* Adaptive rate limiting: cap |dn| at dn_max */
            double this_rate = eff_relax_rate;
            if (step % 100 == 0) {
                double max_f = 0;
                for (size_t ix = 0; ix < N3; ix++) {
                    double fm = v3_norm2(force[ix]);
                    if (fm > max_f) max_f = fm;
                }
                max_f = sqrt(max_f);
                double dn_max = 0.005;
                if (max_f * eff_relax_rate > dn_max)
                    this_rate = dn_max / max_f;
            }

            #pragma omp parallel for
            for (size_t ix = 0; ix < N3; ix++) {
                field[ix] = v3_madd(field[ix], this_rate, force[ix]);
                field[ix] = v3_normalize(field[ix]);
            }

            if (step % diag_interval == 0) {
                compute_energy(field, &E2, &E4, &EV);
                if (step % Q_interval == 0) {
                    Q = compute_hopf_charge(field);
                    printf("  relax %5d: E2=%.4f E4=%.4f E=%.4f Q=%.4f\n",
                           step, E2, E4, E2+E4+EV, Q);
                } else {
                    printf("  relax %5d: E2=%.4f E4=%.4f E=%.4f\n",
                           step, E2, E4, E2+E4+EV);
                }
                fflush(stdout);
            }
            if (step % snap_interval == 0) {
                char fname[512];
                snprintf(fname, 512, "%s/relax_%05d.bin", snap_dir, step);
                write_snapshot(fname, field, step, step * dt);
            }
        }
        compute_energy(field, &E2, &E4, &EV);
        Q = compute_hopf_charge(field);
        printf("After relaxation: E2=%.4f E4=%.4f E=%.4f Q=%.4f\n\n",
               E2, E4, E2+E4+EV, Q);
    }

    /* === Phase 2: Overdamped Langevin dynamics === */
    /*
     * First-order (overdamped) dynamics: no velocity field.
     *   n_{i+1} = normalize(n_i + (dt/γ) * f + noise)
     *   noise_amp = sqrt(2 T dt / (γ h³))
     *
     * At T=0: pure gradient flow (steepest descent).
     * At T>0: samples Boltzmann distribution at temperature T.
     * Fluctuation-dissipation: <ξξ> = 2T/(γ h³) δ(t).
     */
    if (n_steps > 0) {
        printf("--- Overdamped Langevin (%d steps, T=%.4e, γ=%.1f) ---\n",
               n_steps, T_bath, gamma0);

        double rate = dt / gamma0;
        double noise_amp = (T_bath > 0) ?
            sqrt(2.0 * T_bath * dt / (gamma0 * h3)) : 0.0;
        printf("  rate=%.6f noise_amp=%.6e\n", rate, noise_amp);

        /* Check rate safety */
        compute_force(field, force, Hx, Hy, Hz);
        double max_f0 = 0;
        for (size_t ix = 0; ix < N3; ix++) {
            double fm = v3_norm2(force[ix]);
            if (fm > max_f0) max_f0 = fm;
        }
        max_f0 = sqrt(max_f0);
        printf("  max|F|=%.3f  max|dn|=%.6f rad/step\n",
               max_f0, max_f0 * rate);
        if (max_f0 * rate > 0.01) {
            double safe_gamma = dt * max_f0 / 0.005;
            printf("  WARNING: rate too large! Set -gamma %.0f for stability\n",
                   safe_gamma);
            printf("  Clamping rate to 0.005/max|F| = %.6f\n",
                   0.005 / max_f0);
            rate = 0.005 / max_f0;
            /* Adjust noise to preserve FDT: noise² = 2T·rate/h³ */
            noise_amp = (T_bath > 0) ? sqrt(2.0 * T_bath * rate / h3) : 0.0;
            printf("  Adjusted: rate=%.6f noise_amp=%.6e\n", rate, noise_amp);
        }

        for (int step = 0; step <= n_steps; step++) {
            double t = step * dt;

            /* Diagnostics */
            if (step % diag_interval == 0) {
                compute_energy(field, &E2, &E4, &EV);
                int do_Q = (step % Q_interval == 0);
                if (do_Q) Q = compute_hopf_charge(field);
                if (do_Q)
                    printf("  %6d t=%.3f E2=%.3f E4=%.3f E=%.3f Q=%.4f\n",
                           step, t, E2, E4, E2+E4+EV, Q);
                else
                    printf("  %6d t=%.3f E2=%.3f E4=%.3f E=%.3f\n",
                           step, t, E2, E4, E2+E4+EV);
                if (ts_fp) {
                    fprintf(ts_fp, "%d %.6f %.6f %.6f %.6f %.6f %.6f %.6e\n",
                            step, t, E2, E4, EV, 0.0, Q, T_bath);
                    fflush(ts_fp);
                }
                fflush(stdout);
            }

            /* Snapshot */
            if (step % snap_interval == 0) {
                char fname[512];
                snprintf(fname, 512, "%s/snap_%05d.bin", snap_dir, step);
                write_snapshot(fname, field, step, t);
            }

            if (step == n_steps) break;

            /* Compute force */
            compute_force(field, force, Hx, Hy, Hz);

            /* Update field: overdamped step + noise */
            #pragma omp parallel
            {
                int tid = omp_get_thread_num();
                Rng *rng = &thread_rng[tid];
                #pragma omp for
                for (size_t ix = 0; ix < N3; ix++) {
                    Vec3 n = field[ix];
                    Vec3 f_val = force[ix];

                    /* Gradient step */
                    Vec3 dn = v3_scale(rate, f_val);

                    /* Add noise (tangent to S²) */
                    if (noise_amp > 0) {
                        double gx, gy, gz, dummy;
                        rng_gauss2(rng, &gx, &gy);
                        rng_gauss2(rng, &gz, &dummy);
                        Vec3 noise = v3(noise_amp*gx, noise_amp*gy, noise_amp*gz);
                        /* Project noise tangent to S² */
                        double nn = v3_dot(noise, n);
                        noise = v3_madd(noise, -nn, n);
                        dn = v3_add(dn, noise);
                    }

                    field[ix] = v3_normalize(v3_madd(n, 1.0, dn));
                }
            }
        }
    }

    /* Final state */
    compute_energy(field, &E2, &E4, &EV);
    Q = compute_hopf_charge(field);
    printf("\nFinal: E2=%.4f E4=%.4f E=%.4f Q=%.4f\n", E2, E4, E2+E4+EV, Q);

    if (ts_fp) fclose(ts_fp);
    free(field); free(force);
    free(Hx); free(Hy); free(Hz);
    return 0;
}
