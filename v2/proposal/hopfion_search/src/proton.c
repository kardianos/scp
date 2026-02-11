/*
 * proton.c — Langevin field dynamics simulator for Skyrmion thermal bath
 *
 * Tests the hypothesis that solitons require T>0 thermal fluctuations
 * to sustain themselves as self-reinforcing dynamic patterns.
 *
 * Modes:
 *   hedgehog: B=1 hedgehog from 1D profile (control)
 *   quarks:   Three individually-unstable field lumps at triangle vertices
 *   b3:       B=3 rational map (three visible quarks)
 *   extract:  Single quark extracted from B=3 (dissolves alone)
 *
 * Physics: Skyrme L2+L4 with Langevin thermostat (damping + noise).
 * Method: Leapfrog time-stepping with tangent-space projection on S3.
 *
 * Self-contained: no dependencies beyond -lm and -fopenmp.
 * Compile: gcc -O3 -fopenmp -std=c11 -o proton proton.c -lm
 */

#define _USE_MATH_DEFINES
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdint.h>
#include <time.h>

#ifdef _OPENMP
#include <omp.h>
#else
static inline int omp_get_max_threads(void) { return 1; }
static inline int omp_get_thread_num(void) { return 0; }
static inline int omp_get_num_threads(void) { return 1; }
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ========== Quaternion type and helpers ========== */

typedef struct { double s, f1, f2, f3; } Q4;

static inline Q4 q4_mul(Q4 a, Q4 b) {
    return (Q4){
        a.s*b.s  - a.f1*b.f1 - a.f2*b.f2 - a.f3*b.f3,
        a.s*b.f1 + a.f1*b.s  - a.f2*b.f3 + a.f3*b.f2,
        a.s*b.f2 + a.f1*b.f3 + a.f2*b.s  - a.f3*b.f1,
        a.s*b.f3 - a.f1*b.f2 + a.f2*b.f1 + a.f3*b.s
    };
}

static inline Q4 q4_rev(Q4 a) {
    return (Q4){a.s, -a.f1, -a.f2, -a.f3};
}

static inline Q4 q4_scale(double c, Q4 a) {
    return (Q4){c*a.s, c*a.f1, c*a.f2, c*a.f3};
}

static inline Q4 q4_add(Q4 a, Q4 b) {
    return (Q4){a.s+b.s, a.f1+b.f1, a.f2+b.f2, a.f3+b.f3};
}

static inline Q4 q4_sub(Q4 a, Q4 b) {
    return (Q4){a.s-b.s, a.f1-b.f1, a.f2-b.f2, a.f3-b.f3};
}

static inline double q4_norm2(Q4 a) {
    return a.s*a.s + a.f1*a.f1 + a.f2*a.f2 + a.f3*a.f3;
}

static inline double q4_dot(Q4 a, Q4 b) {
    return a.s*b.s + a.f1*b.f1 + a.f2*b.f2 + a.f3*b.f3;
}

static inline Q4 q4_zero(void) { return (Q4){0,0,0,0}; }

/* ========== Grid index with periodic BCs ========== */

static inline int idx(int N, int i, int j, int k) {
    return ((i % N + N) % N) * N * N
         + ((j % N + N) % N) * N
         + ((k % N + N) % N);
}

/* ========== RNG: xoshiro256** (thread-safe) ========== */

typedef struct { uint64_t s[4]; } Rng;

static inline uint64_t rotl(uint64_t x, int k) {
    return (x << k) | (x >> (64 - k));
}

static inline uint64_t rng_next(Rng *rng) {
    uint64_t *s = rng->s;
    uint64_t result = rotl(s[1] * 5, 7) * 9;
    uint64_t t = s[1] << 17;
    s[2] ^= s[0]; s[3] ^= s[1]; s[1] ^= s[2]; s[0] ^= s[3];
    s[2] ^= t;
    s[3] = rotl(s[3], 45);
    return result;
}

static inline double rng_uniform(Rng *rng) {
    return (rng_next(rng) >> 11) * 0x1.0p-53;
}

/* Box-Muller: generate two independent N(0,1) samples */
static inline void rng_gauss2(Rng *rng, double *g1, double *g2) {
    double u1 = rng_uniform(rng);
    double u2 = rng_uniform(rng);
    if (u1 < 1e-300) u1 = 1e-300;
    double r = sqrt(-2.0 * log(u1));
    double th = 2.0 * M_PI * u2;
    *g1 = r * cos(th);
    *g2 = r * sin(th);
}

static void rng_seed(Rng *rng, uint64_t seed) {
    /* SplitMix64 to initialize state */
    for (int i = 0; i < 4; i++) {
        seed += 0x9e3779b97f4a7c15ULL;
        uint64_t z = seed;
        z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
        z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
        rng->s[i] = z ^ (z >> 31);
    }
}

/* ========== Profile loading ========== */

typedef struct {
    double *r, *f;
    int n;
    double dr, r_max;
} RadialProfile;

static RadialProfile *load_profile(const char *filename) {
    FILE *fp = fopen(filename, "r");
    if (!fp) { fprintf(stderr, "Cannot open %s\n", filename); return NULL; }

    int n = 0;
    char line[1024];
    while (fgets(line, sizeof(line), fp)) {
        if (line[0] == '#') continue;
        double rv, fv;
        if (sscanf(line, "%lf %lf", &rv, &fv) >= 2) n++;
    }
    rewind(fp);

    RadialProfile *p = malloc(sizeof(RadialProfile));
    p->r = malloc(n * sizeof(double));
    p->f = malloc(n * sizeof(double));
    p->n = n;

    int i = 0;
    while (fgets(line, sizeof(line), fp) && i < n) {
        if (line[0] == '#') continue;
        double rv, fv;
        if (sscanf(line, "%lf %lf", &rv, &fv) >= 2) {
            p->r[i] = rv;
            p->f[i] = fv;
            i++;
        }
    }
    fclose(fp);

    p->dr = (n > 1) ? p->r[1] - p->r[0] : 0.001;
    p->r_max = p->r[n-1];
    return p;
}

static void free_profile(RadialProfile *p) {
    if (p) { free(p->r); free(p->f); free(p); }
}

static double interp_f(const RadialProfile *p, double r) {
    if (r < 0) r = -r;
    if (r >= p->r_max) return 0.0;
    double fi = r / p->dr;
    int i = (int)fi;
    if (i >= p->n - 1) return 0.0;
    double t = fi - i;
    return (1-t)*p->f[i] + t*p->f[i+1];
}

/* ========== Hedgehog quaternion ========== */

static Q4 hedgehog_q(const RadialProfile *prof, double rho0,
                     double x, double y, double z,
                     double cx, double cy, double cz)
{
    double dx = x-cx, dy = y-cy, dz = z-cz;
    double r = sqrt(dx*dx + dy*dy + dz*dz);
    double f = interp_f(prof, r);
    double cf = cos(f), sf = sin(f);
    Q4 q;
    q.s = rho0 * cf;
    if (r > 1e-12) {
        double sr = rho0 * sf / r;
        q.f1 = sr * dx;
        q.f2 = sr * dy;
        q.f3 = sr * dz;
    } else {
        q.f1 = q.f2 = q.f3 = 0;
    }
    return q;
}

/* ========== Rational map for B=3 ========== */

static void rational_map_n(int B, double theta, double phi,
                           double *nx, double *ny, double *nz)
{
    double t2 = tan(0.5 * theta);
    double zr = t2 * cos(phi);
    double zi = t2 * sin(phi);

    double Rr, Ri;

    if (B == 1) {
        Rr = zr; Ri = zi;
    } else if (B == 2) {
        Rr = zr*zr - zi*zi;
        Ri = 2*zr*zi;
    } else if (B == 3) {
        double s3 = sqrt(3.0);
        double z2r = zr*zr - zi*zi;
        double z2i = 2*zr*zi;
        double z3r = z2r*zr - z2i*zi;
        double z3i = z2r*zi + z2i*zr;
        double s3izr = -s3*zi;
        double s3izi =  s3*zr;
        double s3iz2r = -s3*z2i;
        double s3iz2i =  s3*z2r;
        double nr = z3r - s3izr;
        double ni = z3i - s3izi;
        double dr = s3iz2r - 1.0;
        double di = s3iz2i;
        double d2 = dr*dr + di*di;
        if (d2 < 1e-30) { Rr = 1e15; Ri = 0; }
        else { Rr = (nr*dr + ni*di)/d2; Ri = (ni*dr - nr*di)/d2; }
    } else if (B == 4) {
        double s3 = sqrt(3.0);
        double z2r = zr*zr - zi*zi;
        double z2i = 2*zr*zi;
        double z4r = z2r*z2r - z2i*z2i;
        double z4i = 2*z2r*z2i;
        double tir = -2*s3*z2i;
        double tii =  2*s3*z2r;
        double nr = z4r + tir + 1.0;
        double ni = z4i + tii;
        double dr = z4r - tir + 1.0;
        double di = z4i - tii;
        double d2 = dr*dr + di*di;
        if (d2 < 1e-30) { Rr = 1e15; Ri = 0; }
        else { Rr = (nr*dr + ni*di)/d2; Ri = (ni*dr - nr*di)/d2; }
    } else {
        Rr = zr; Ri = zi;
    }

    double R2 = Rr*Rr + Ri*Ri;
    double inv = 1.0 / (1.0 + R2);
    *nx = 2.0 * Rr * inv;
    *ny = 2.0 * Ri * inv;
    *nz = (1.0 - R2) * inv;
}

/* ========== 4th-order derivative helper ========== */

static inline Q4 q4_deriv(const Q4 *q, int N, int i, int j, int k, int dir, double inv12h) {
    int im2, im1, ip1, ip2;
    switch (dir) {
    case 0:
        im2 = idx(N,i-2,j,k); im1 = idx(N,i-1,j,k);
        ip1 = idx(N,i+1,j,k); ip2 = idx(N,i+2,j,k);
        break;
    case 1:
        im2 = idx(N,i,j-2,k); im1 = idx(N,i,j-1,k);
        ip1 = idx(N,i,j+1,k); ip2 = idx(N,i,j+2,k);
        break;
    default:
        im2 = idx(N,i,j,k-2); im1 = idx(N,i,j,k-1);
        ip1 = idx(N,i,j,k+1); ip2 = idx(N,i,j,k+2);
        break;
    }
    Q4 r;
    r.s  = (-q[ip2].s  + 8*q[ip1].s  - 8*q[im1].s  + q[im2].s)  * inv12h;
    r.f1 = (-q[ip2].f1 + 8*q[ip1].f1 - 8*q[im1].f1 + q[im2].f1) * inv12h;
    r.f2 = (-q[ip2].f2 + 8*q[ip1].f2 - 8*q[im1].f2 + q[im2].f2) * inv12h;
    r.f3 = (-q[ip2].f3 + 8*q[ip1].f3 - 8*q[im1].f3 + q[im2].f3) * inv12h;
    return r;
}

/* Consistent 9-point Laplacian per direction: {1,-16,64,16,-130,16,64,-16,1}/(144h^2) */
static inline Q4 q4_laplacian(const Q4 *q, int N, int i, int j, int k, double inv144h2) {
    static const double w[9] = {1, -16, 64, 16, -130, 16, 64, -16, 1};
    double rs = 0, rf1 = 0, rf2 = 0, rf3 = 0;

    for (int m = -4; m <= 4; m++) {
        double wm = w[m+4];
        int ix;
        ix = idx(N, i+m, j, k);
        rs += wm*q[ix].s; rf1 += wm*q[ix].f1; rf2 += wm*q[ix].f2; rf3 += wm*q[ix].f3;
        ix = idx(N, i, j+m, k);
        rs += wm*q[ix].s; rf1 += wm*q[ix].f1; rf2 += wm*q[ix].f2; rf3 += wm*q[ix].f3;
        ix = idx(N, i, j, k+m);
        rs += wm*q[ix].s; rf1 += wm*q[ix].f1; rf2 += wm*q[ix].f2; rf3 += wm*q[ix].f3;
    }
    return (Q4){inv144h2*rs, inv144h2*rf1, inv144h2*rf2, inv144h2*rf3};
}

/* ========== Left-multiply by basis quaternion ========== */

static inline Q4 q4_left_basis(int a, Q4 q) {
    switch(a) {
    case 0: return q;
    case 1: return (Q4){-q.f1,  q.s,   q.f3, -q.f2};
    case 2: return (Q4){-q.f2, -q.f3,  q.s,   q.f1};
    case 3: return (Q4){-q.f3,  q.f2, -q.f1,  q.s};
    default: return q4_zero();
    }
}

/* ========== Skyrme predata and force ========== */

typedef struct {
    Q4 A[3];
    Q4 G[3];
} SkyrmePre;

static void compute_skyrme_predata(const Q4 *q, int N, double h, SkyrmePre *pre) {
    double inv12h = 1.0 / (12.0 * h);
    #pragma omp parallel for collapse(3) schedule(static)
    for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
    for (int k = 0; k < N; k++) {
        int ix = idx(N,i,j,k);
        Q4 qr = q4_rev(q[ix]);
        Q4 dq[3];
        for (int d = 0; d < 3; d++)
            dq[d] = q4_deriv(q, N, i, j, k, d, inv12h);
        for (int d = 0; d < 3; d++)
            pre[ix].A[d] = q4_mul(qr, dq[d]);

        Q4 C[3];
        C[0] = q4_sub(q4_mul(pre[ix].A[0], pre[ix].A[1]),
                       q4_mul(pre[ix].A[1], pre[ix].A[0]));
        C[1] = q4_sub(q4_mul(pre[ix].A[0], pre[ix].A[2]),
                       q4_mul(pre[ix].A[2], pre[ix].A[0]));
        C[2] = q4_sub(q4_mul(pre[ix].A[1], pre[ix].A[2]),
                       q4_mul(pre[ix].A[2], pre[ix].A[1]));

        for (int d = 0; d < 3; d++) {
            pre[ix].G[d] = q4_zero();
            for (int dp = 0; dp < 3; dp++) {
                if (dp == d) continue;
                Q4 Cdd;
                if (d < dp) {
                    int ci = (d==0 && dp==1) ? 0 : (d==0 && dp==2) ? 1 : 2;
                    Cdd = C[ci];
                } else {
                    int ci = (dp==0 && d==1) ? 0 : (dp==0 && d==2) ? 1 : 2;
                    Cdd = q4_scale(-1.0, C[ci]);
                }
                Q4 term = q4_sub(q4_mul(pre[ix].A[dp], Cdd),
                                 q4_mul(Cdd, pre[ix].A[dp]));
                pre[ix].G[d] = q4_add(pre[ix].G[d], term);
            }
        }
    }
}

typedef struct { double E2, E4; } EnergyResult;

static EnergyResult compute_energy(const Q4 *q, int N, double h, double e_skyrme) {
    double h3 = h*h*h;
    double inv12h = 1.0 / (12.0 * h);
    double inv_4e2 = 1.0 / (4.0 * e_skyrme * e_skyrme);
    double e2_sum = 0, e4_sum = 0;

    #pragma omp parallel for collapse(3) reduction(+:e2_sum,e4_sum) schedule(static)
    for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
    for (int k = 0; k < N; k++) {
        int ix = idx(N,i,j,k);
        Q4 qr = q4_rev(q[ix]);
        Q4 dq[3], A[3];
        for (int d = 0; d < 3; d++) {
            dq[d] = q4_deriv(q, N, i, j, k, d, inv12h);
            e2_sum += 0.5 * q4_norm2(dq[d]);
            A[d] = q4_mul(qr, dq[d]);
        }
        for (int d1 = 0; d1 < 3; d1++)
        for (int d2 = d1+1; d2 < 3; d2++) {
            Q4 comm = q4_sub(q4_mul(A[d1], A[d2]), q4_mul(A[d2], A[d1]));
            Q4 comm2 = q4_mul(comm, comm);
            e4_sum -= inv_4e2 * comm2.s;
        }
    }

    EnergyResult en;
    en.E2 = e2_sum * h3;
    en.E4 = e4_sum * h3;
    return en;
}

static void compute_skyrme_force(const Q4 *q, int N, double h, double e_skyrme, Q4 *force) {
    size_t N3 = (size_t)N*N*N;
    double inv12h = 1.0 / (12.0 * h);
    double inv_2e2 = 1.0 / (2.0 * e_skyrme * e_skyrme);
    double inv144h2 = 1.0 / (144.0 * h * h);
    double sigma[4] = {1.0, -1.0, -1.0, -1.0};

    SkyrmePre *pre = malloc(N3 * sizeof(SkyrmePre));
    compute_skyrme_predata(q, N, h, pre);

    double *pi = calloc(N3 * 12, sizeof(double));

    #pragma omp parallel for collapse(3) schedule(static)
    for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
    for (int k = 0; k < N; k++) {
        int ix = idx(N,i,j,k);
        Q4 qr = q4_rev(q[ix]);
        for (int d = 0; d < 3; d++) {
            for (int a = 0; a < 4; a++) {
                Q4 qr_ea;
                switch (a) {
                case 0: qr_ea = qr; break;
                case 1: qr_ea = q4_mul(qr, (Q4){0,1,0,0}); break;
                case 2: qr_ea = q4_mul(qr, (Q4){0,0,1,0}); break;
                case 3: qr_ea = q4_mul(qr, (Q4){0,0,0,1}); break;
                default: qr_ea = q4_zero();
                }
                Q4 prod = q4_mul(qr_ea, pre[ix].G[d]);
                pi[12*ix + 3*a + d] = prod.s;
            }
        }
    }

    #pragma omp parallel for collapse(3) schedule(static)
    for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
    for (int k = 0; k < N; k++) {
        int ix = idx(N,i,j,k);

        Q4 lap = q4_laplacian(q, N, i, j, k, inv144h2);
        force[ix].s  = lap.s;
        force[ix].f1 = lap.f1;
        force[ix].f2 = lap.f2;
        force[ix].f3 = lap.f3;

        double f4[4] = {0, 0, 0, 0};
        for (int d = 0; d < 3; d++) {
            Q4 dq = q4_deriv(q, N, i, j, k, d, inv12h);
            for (int a = 0; a < 4; a++) {
                Q4 ea_dq = q4_left_basis(a, dq);
                Q4 prod = q4_mul(ea_dq, pre[ix].G[d]);
                f4[a] += sigma[a] * prod.s;
            }
        }
        for (int a = 0; a < 4; a++) {
            for (int d = 0; d < 3; d++) {
                int im2, im1, ip1, ip2;
                switch (d) {
                case 0:
                    im2 = idx(N,i-2,j,k); im1 = idx(N,i-1,j,k);
                    ip1 = idx(N,i+1,j,k); ip2 = idx(N,i+2,j,k);
                    break;
                case 1:
                    im2 = idx(N,i,j-2,k); im1 = idx(N,i,j-1,k);
                    ip1 = idx(N,i,j+1,k); ip2 = idx(N,i,j+2,k);
                    break;
                default:
                    im2 = idx(N,i,j,k-2); im1 = idx(N,i,j,k-1);
                    ip1 = idx(N,i,j,k+1); ip2 = idx(N,i,j,k+2);
                    break;
                }
                double dpi = (-pi[12*ip2 + 3*a + d] + 8*pi[12*ip1 + 3*a + d]
                              -8*pi[12*im1 + 3*a + d] + pi[12*im2 + 3*a + d]) * inv12h;
                f4[a] -= dpi;
            }
        }

        force[ix].s  += inv_2e2 * f4[0];
        force[ix].f1 += inv_2e2 * f4[1];
        force[ix].f2 += inv_2e2 * f4[2];
        force[ix].f3 += inv_2e2 * f4[3];
    }

    free(pi);
    free(pre);
}

/* ========== Topological charge ========== */

static double compute_Q(const Q4 *q, int N, double h) {
    double h3 = h*h*h;
    double inv12h = 1.0 / (12.0 * h);
    double sum = 0;

    #pragma omp parallel for collapse(3) reduction(+:sum) schedule(static)
    for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
    for (int k = 0; k < N; k++) {
        int ix = idx(N,i,j,k);
        Q4 qr = q4_rev(q[ix]);
        Q4 dq[3], A[3];
        for (int d = 0; d < 3; d++) {
            dq[d] = q4_deriv(q, N, i, j, k, d, inv12h);
            A[d] = q4_mul(qr, dq[d]);
        }
        Q4 A012 = q4_mul(q4_mul(A[0], A[1]), A[2]);
        double norm2 = q4_norm2(q[ix]);
        double norm6 = norm2 * norm2 * norm2;
        if (norm6 > 1e-30)
            sum += A012.s / norm6;
    }
    return -sum * h3 / (2.0 * M_PI * M_PI);
}

/* ========== Binary snapshot output ========== */

static void write_snapshot(const char *fname, const Q4 *q, int N, double L,
                           int step, double time_val) {
    FILE *fp = fopen(fname, "wb");
    if (!fp) { fprintf(stderr, "Cannot write %s\n", fname); return; }
    fwrite(&N, sizeof(int), 1, fp);
    fwrite(&L, sizeof(double), 1, fp);
    fwrite(&step, sizeof(int), 1, fp);
    fwrite(&time_val, sizeof(double), 1, fp);
    fwrite(q, sizeof(Q4), (size_t)N*N*N, fp);
    fclose(fp);
}

/* ========== Self-tests ========== */

static int selftest_force(int N, double L, double e_skyrme, double rho0,
                          const RadialProfile *prof) {
    printf("  Force consistency test...\n");
    double h = 2.0 * L / N;
    size_t N3 = (size_t)N*N*N;

    Q4 *q_test = calloc(N3, sizeof(Q4));
    Q4 *f_test = calloc(N3, sizeof(Q4));

    /* Initialize hedgehog */
    #pragma omp parallel for collapse(3) schedule(static)
    for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
    for (int k = 0; k < N; k++) {
        double x = -L + (i + 0.5) * h;
        double y = -L + (j + 0.5) * h;
        double z = -L + (k + 0.5) * h;
        int ix = idx(N,i,j,k);
        q_test[ix] = hedgehog_q(prof, rho0, x, y, z, 0, 0, 0);
        double n = sqrt(q4_norm2(q_test[ix]));
        if (n > 1e-15) q_test[ix] = q4_scale(rho0/n, q_test[ix]);
    }

    /* Compute energy */
    EnergyResult E0 = compute_energy(q_test, N, h, e_skyrme);

    /* Compute force */
    memset(f_test, 0, N3 * sizeof(Q4));
    compute_skyrme_force(q_test, N, h, e_skyrme, f_test);

    /* Verify: perturb field, check dE ~ -F.dq */
    double eps = 1e-5;
    int test_i = N/2, test_j = N/2, test_k = N/2;
    int test_ix = idx(N, test_i, test_j, test_k);

    /* Perturb f1 component */
    Q4 q_save = q_test[test_ix];
    q_test[test_ix].f1 += eps;
    double n = sqrt(q4_norm2(q_test[test_ix]));
    q_test[test_ix] = q4_scale(rho0/n, q_test[test_ix]);

    EnergyResult E1 = compute_energy(q_test, N, h, e_skyrme);
    double dE = (E1.E2 + E1.E4) - (E0.E2 + E0.E4);

    /* Restore */
    q_test[test_ix] = q_save;

    double h3 = h*h*h;

    /* For sigma model: F_tangent . dq_tangent should approximate -dE/h^3 */
    /* Finite diff gives dE/dq_a ~ -F_a * h^3  =>  dE ~ -F.dq * h^3 */
    double F_dq = f_test[test_ix].f1 * eps;  /* approximate */
    double predicted_dE = -F_dq * h3;

    double rel_err = (dE != 0) ? fabs(dE - predicted_dE) / fabs(dE) : fabs(predicted_dE);

    printf("    dE(actual)=%.6e, dE(predicted)=%.6e, rel_err=%.2e %s\n",
           dE, predicted_dE, rel_err, rel_err < 0.1 ? "OK" : "MARGINAL");

    free(q_test);
    free(f_test);
    return rel_err < 0.5;  /* generous for finite-diff */
}

static int selftest_noise(void) {
    printf("  Noise statistics test...\n");
    Rng rng;
    rng_seed(&rng, 42);

    int N_samples = 100000;
    double sum = 0, sum2 = 0;
    for (int i = 0; i < N_samples; i++) {
        double g1, g2;
        rng_gauss2(&rng, &g1, &g2);
        sum += g1 + g2;
        sum2 += g1*g1 + g2*g2;
    }
    double mean = sum / (2*N_samples);
    double var = sum2 / (2*N_samples) - mean*mean;

    printf("    mean=%.4e (expect 0), var=%.4f (expect 1) %s\n",
           mean, var,
           (fabs(mean) < 0.01 && fabs(var - 1.0) < 0.02) ? "OK" : "FAIL");
    return (fabs(mean) < 0.01 && fabs(var - 1.0) < 0.02);
}

static int selftest_conservation(int N, double L, double dt, double e_skyrme,
                                 double rho0, const RadialProfile *prof) {
    printf("  Energy conservation test (T=0, gamma=0, 100 steps)...\n");
    double h = 2.0 * L / N;
    size_t N3 = (size_t)N*N*N;

    Q4 *q_test = calloc(N3, sizeof(Q4));
    Q4 *v_test = calloc(N3, sizeof(Q4));
    Q4 *f_test = calloc(N3, sizeof(Q4));

    /* Initialize hedgehog */
    #pragma omp parallel for collapse(3) schedule(static)
    for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
    for (int k = 0; k < N; k++) {
        double x = -L + (i + 0.5) * h;
        double y = -L + (j + 0.5) * h;
        double z = -L + (k + 0.5) * h;
        int ix = idx(N,i,j,k);
        q_test[ix] = hedgehog_q(prof, rho0, x, y, z, 0, 0, 0);
        double n = sqrt(q4_norm2(q_test[ix]));
        if (n > 1e-15) q_test[ix] = q4_scale(rho0/n, q_test[ix]);
        v_test[ix] = q4_zero();
    }

    /* Give small initial velocity perturbation */
    int mid = N/2;
    v_test[idx(N,mid,mid,mid)].f1 = 0.1;
    /* Project tangent */
    {
        int ix = idx(N,mid,mid,mid);
        double vdq = q4_dot(v_test[ix], q_test[ix]) / q4_norm2(q_test[ix]);
        v_test[ix] = q4_sub(v_test[ix], q4_scale(vdq, q_test[ix]));
    }

    /* Initial energy */
    EnergyResult E0 = compute_energy(q_test, N, h, e_skyrme);
    double h3 = h*h*h;
    double Ekin0 = 0;
    for (size_t ix = 0; ix < N3; ix++)
        Ekin0 += 0.5 * q4_norm2(v_test[ix]);
    Ekin0 *= h3;
    double Etot0 = E0.E2 + E0.E4 + Ekin0;

    /* Leapfrog half-step */
    memset(f_test, 0, N3 * sizeof(Q4));
    compute_skyrme_force(q_test, N, h, e_skyrme, f_test);
    #pragma omp parallel for
    for (size_t ix = 0; ix < N3; ix++) {
        double n2 = q4_norm2(q_test[ix]);
        if (n2 > 1e-20) {
            double fdq = q4_dot(f_test[ix], q_test[ix]) / n2;
            f_test[ix] = q4_sub(f_test[ix], q4_scale(fdq, q_test[ix]));
        }
    }
    #pragma omp parallel for
    for (size_t ix = 0; ix < N3; ix++) {
        v_test[ix].s  += 0.5*dt * f_test[ix].s;
        v_test[ix].f1 += 0.5*dt * f_test[ix].f1;
        v_test[ix].f2 += 0.5*dt * f_test[ix].f2;
        v_test[ix].f3 += 0.5*dt * f_test[ix].f3;
    }

    /* 100 leapfrog steps */
    for (int step = 0; step < 100; step++) {
        #pragma omp parallel for
        for (size_t ix = 0; ix < N3; ix++) {
            q_test[ix].s  += dt * v_test[ix].s;
            q_test[ix].f1 += dt * v_test[ix].f1;
            q_test[ix].f2 += dt * v_test[ix].f2;
            q_test[ix].f3 += dt * v_test[ix].f3;
        }
        /* Sigma-project */
        #pragma omp parallel for
        for (size_t ix = 0; ix < N3; ix++) {
            double n = sqrt(q4_norm2(q_test[ix]));
            if (n > 1e-15) q_test[ix] = q4_scale(rho0/n, q_test[ix]);
            double vdq = q4_dot(v_test[ix], q_test[ix]) / (rho0*rho0);
            v_test[ix] = q4_sub(v_test[ix], q4_scale(vdq, q_test[ix]));
        }

        memset(f_test, 0, N3 * sizeof(Q4));
        compute_skyrme_force(q_test, N, h, e_skyrme, f_test);
        #pragma omp parallel for
        for (size_t ix = 0; ix < N3; ix++) {
            double n2 = q4_norm2(q_test[ix]);
            if (n2 > 1e-20) {
                double fdq = q4_dot(f_test[ix], q_test[ix]) / n2;
                f_test[ix] = q4_sub(f_test[ix], q4_scale(fdq, q_test[ix]));
            }
        }
        #pragma omp parallel for
        for (size_t ix = 0; ix < N3; ix++) {
            v_test[ix].s  += dt * f_test[ix].s;
            v_test[ix].f1 += dt * f_test[ix].f1;
            v_test[ix].f2 += dt * f_test[ix].f2;
            v_test[ix].f3 += dt * f_test[ix].f3;
        }
    }

    /* Final energy */
    EnergyResult E1 = compute_energy(q_test, N, h, e_skyrme);
    double Ekin1 = 0;
    for (size_t ix = 0; ix < N3; ix++)
        Ekin1 += 0.5 * q4_norm2(v_test[ix]);
    Ekin1 *= h3;
    double Etot1 = E1.E2 + E1.E4 + Ekin1;

    double dE_rel = fabs(Etot1 - Etot0) / (fabs(Etot0) + 1e-30);
    printf("    E0=%.6f, E100=%.6f, dE/E=%.2e %s\n",
           Etot0, Etot1, dE_rel, dE_rel < 0.01 ? "OK" : "FAIL");

    free(q_test); free(v_test); free(f_test);
    return dE_rel < 0.01;
}

/* ========== Initialization modes ========== */

enum Mode { MODE_HEDGEHOG, MODE_QUARKS, MODE_B3, MODE_EXTRACT };

static void init_hedgehog(Q4 *q, Q4 *vel, int N, double L, double h,
                          double rho0, const RadialProfile *prof) {
    #pragma omp parallel for collapse(3) schedule(static)
    for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
    for (int k = 0; k < N; k++) {
        double x = -L + (i + 0.5) * h;
        double y = -L + (j + 0.5) * h;
        double z = -L + (k + 0.5) * h;
        int ix = idx(N,i,j,k);
        q[ix] = hedgehog_q(prof, rho0, x, y, z, 0, 0, 0);
        double n = sqrt(q4_norm2(q[ix]));
        if (n > 1e-15) q[ix] = q4_scale(rho0/n, q[ix]);
        vel[ix] = q4_zero();
    }
}

static void init_three_quarks(Q4 *q, Q4 *vel, int N, double L, double h,
                              double rho0, double D, double A_amp, double sigma) {
    /* Three Gaussian perturbations at equilateral triangle vertices in xy-plane.
     * Each in a different SU(2) direction: f1, f2, f3.
     * g(r) = A * exp(-r^2/(2*sigma^2)), q = (cos g, sin g * e_dir) */
    double s3 = sqrt(3.0);
    double qx[3] = {+D,      -D/2.0,       -D/2.0};
    double qy[3] = { 0, +D*s3/2.0, -D*s3/2.0};
    double qz[3] = { 0,           0,           0};

    #pragma omp parallel for collapse(3) schedule(static)
    for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
    for (int k = 0; k < N; k++) {
        double x = -L + (i + 0.5) * h;
        double y = -L + (j + 0.5) * h;
        double z = -L + (k + 0.5) * h;
        int ix = idx(N,i,j,k);

        /* Start with vacuum */
        Q4 qtot = {rho0, 0, 0, 0};

        /* Compose three quark perturbations via quaternion product */
        for (int qk = 0; qk < 3; qk++) {
            double dx = x - qx[qk];
            double dy = y - qy[qk];
            double dz = z - qz[qk];
            double r2 = dx*dx + dy*dy + dz*dz;
            double g = A_amp * exp(-r2 / (2.0 * sigma * sigma));

            Q4 qpert;
            qpert.s = rho0 * cos(g);
            double sg = rho0 * sin(g);
            switch (qk) {
            case 0: qpert.f1 = sg; qpert.f2 = 0;  qpert.f3 = 0;  break;
            case 1: qpert.f1 = 0;  qpert.f2 = sg; qpert.f3 = 0;  break;
            case 2: qpert.f1 = 0;  qpert.f2 = 0;  qpert.f3 = sg; break;
            default: qpert = (Q4){rho0, 0, 0, 0}; break;
            }

            /* Product ansatz: qtot = qtot * qpert / rho0 */
            qtot = q4_scale(1.0/rho0, q4_mul(qtot, qpert));
        }

        /* Normalize to sigma model */
        double n = sqrt(q4_norm2(qtot));
        if (n > 1e-15) qtot = q4_scale(rho0/n, qtot);
        q[ix] = qtot;
        vel[ix] = q4_zero();
    }
}

static void init_b3_ratmap(Q4 *q, Q4 *vel, int N, double L, double h,
                           double rho0, const RadialProfile *prof) {
    #pragma omp parallel for collapse(3) schedule(static)
    for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
    for (int k = 0; k < N; k++) {
        double x = -L + (i + 0.5) * h;
        double y = -L + (j + 0.5) * h;
        double z = -L + (k + 0.5) * h;
        int ix = idx(N,i,j,k);
        double r = sqrt(x*x + y*y + z*z);

        double fval = interp_f(prof, r);
        double cos_f = cos(fval);
        double sin_f = sin(fval);

        double nx_d, ny_d, nz_d;
        if (r > 1e-12) {
            double theta = acos(z / r);
            double phi = atan2(y, x);
            rational_map_n(3, theta, phi, &nx_d, &ny_d, &nz_d);
        } else {
            nx_d = 0; ny_d = 0; nz_d = 1;
        }

        q[ix].s  = rho0 * cos_f;
        q[ix].f1 = rho0 * sin_f * nx_d;
        q[ix].f2 = rho0 * sin_f * ny_d;
        q[ix].f3 = rho0 * sin_f * nz_d;

        double n = sqrt(q4_norm2(q[ix]));
        if (n > 1e-15) q[ix] = q4_scale(rho0/n, q[ix]);
        vel[ix] = q4_zero();
    }
}

static void init_quark_extract(Q4 *q, Q4 *vel, int N, double L, double h,
                               double rho0, const RadialProfile *prof,
                               double R_cut_in) {
    /* First initialize B=3, then mask to vacuum outside a sphere
     * around one of the three baryon density peaks */
    init_b3_ratmap(q, vel, N, L, h, rho0, prof);

    /* The B=3 tetrahedral map has peaks roughly along the three
     * vertices of a tetrahedron. The first peak is near (0, 0, R_peak)
     * where R_peak ~ 1.5 (core size). Extract by masking:
     * keep field within sphere of radius R_cut around peak,
     * smoothly transition to vacuum outside. */
    double peak_x = 0, peak_y = 0, peak_z = 0;  /* approximate peak location */
    double R_cut = R_cut_in;  /* extraction radius */
    double R_smooth = 0.5;  /* smoothing width */

    /* Find actual peak: scan for maximum |pi| */
    double max_pi = 0;
    for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
    for (int k = 0; k < N; k++) {
        int ix = idx(N,i,j,k);
        double pi2 = q[ix].f1*q[ix].f1 + q[ix].f2*q[ix].f2 + q[ix].f3*q[ix].f3;
        if (pi2 > max_pi) {
            max_pi = pi2;
            peak_x = -L + (i + 0.5) * h;
            peak_y = -L + (j + 0.5) * h;
            peak_z = -L + (k + 0.5) * h;
        }
    }
    printf("  Extract peak at (%.2f, %.2f, %.2f), |pi|_max=%.4f\n",
           peak_x, peak_y, peak_z, sqrt(max_pi));

    /* Mask: smooth step from 1 (inside) to 0 (outside) */
    #pragma omp parallel for collapse(3) schedule(static)
    for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
    for (int k = 0; k < N; k++) {
        double x = -L + (i + 0.5) * h;
        double y = -L + (j + 0.5) * h;
        double z = -L + (k + 0.5) * h;
        double dx = x - peak_x, dy = y - peak_y, dz = z - peak_z;
        double r = sqrt(dx*dx + dy*dy + dz*dz);
        int ix = idx(N,i,j,k);

        /* Smooth step: 1 inside R_cut, 0 outside R_cut+R_smooth */
        double mask;
        if (r < R_cut) {
            mask = 1.0;
        } else if (r < R_cut + R_smooth) {
            double t = (r - R_cut) / R_smooth;
            mask = 0.5 * (1.0 + cos(M_PI * t));  /* cosine taper */
        } else {
            mask = 0.0;
        }

        /* Interpolate between current q and vacuum (rho0, 0, 0, 0) */
        q[ix].s  = mask * q[ix].s  + (1.0 - mask) * rho0;
        q[ix].f1 = mask * q[ix].f1;
        q[ix].f2 = mask * q[ix].f2;
        q[ix].f3 = mask * q[ix].f3;

        /* Re-normalize */
        double n = sqrt(q4_norm2(q[ix]));
        if (n > 1e-15) q[ix] = q4_scale(rho0/n, q[ix]);
    }
}

/* ========== Velocity rescaling thermostat ========== */

/* Compute kinetic energy: E_kin = (h^3/2) * sum|v|^2 */
static double compute_ekin(const Q4 *vel, size_t N3, double h3) {
    double Ekin = 0;
    #pragma omp parallel for reduction(+:Ekin)
    for (size_t ix = 0; ix < N3; ix++)
        Ekin += q4_norm2(vel[ix]);
    return 0.5 * Ekin * h3;
}

/* Effective temperature: T_eff = (2/(3*N^3)) * E_kin
 * At equilibrium, each of 3 tangent DOF per site has energy T/2,
 * so E_kin = (3/2)*N^3*T, giving T = (2/3)*E_kin/N^3. */
static double compute_teff(double Ekin, size_t N3) {
    return (N3 > 0) ? (2.0/3.0) * Ekin / (double)N3 : 0.0;
}

/* Rescale velocities to enforce T_eff = T_bath exactly.
 * v -> v * sqrt(T_bath / T_eff).
 * This is the Berendsen velocity-rescaling thermostat. */
static void velocity_rescale(Q4 *vel, const Q4 *field, size_t N3,
                             double h3, double T_bath) {
    if (T_bath <= 0) return;
    double Ekin = compute_ekin(vel, N3, h3);
    double T_eff = compute_teff(Ekin, N3);
    if (T_eff < 1e-30) return;
    double scale = sqrt(T_bath / T_eff);

    #pragma omp parallel for
    for (size_t ix = 0; ix < N3; ix++) {
        vel[ix] = q4_scale(scale, vel[ix]);
        /* Re-project tangent after rescaling */
        double n2 = q4_norm2(field[ix]);
        if (n2 > 1e-20) {
            double vdq = q4_dot(vel[ix], field[ix]) / n2;
            vel[ix] = q4_sub(vel[ix], q4_scale(vdq, field[ix]));
        }
    }
}

/* ========== MAIN ========== */

int main(int argc, char *argv[]) {
    /* Default parameters */
    int N = 128;
    double L = 8.0;
    double e_skyrme = 1.0;
    double rho0 = 1.0;
    double dt = 0.02;
    double T_bath = 0.0;        /* temperature */
    double gamma_damp = 0.5;    /* Langevin damping */
    int n_steps = 5000;
    int snap_interval = 100;
    int diag_interval = 10;
    const char *profile_file = NULL;
    const char *mode_str = "hedgehog";
    double D_quark = 2.0;      /* quark separation */
    double A_quark = 1.5;      /* quark amplitude */
    double sigma_quark = 1.0;  /* quark width */
    double R_cut_extract = 1.5; /* extraction radius for extract mode */
    uint64_t seed = 42;
    int run_selftests = 1;
    int n_therm = 500;          /* thermalization steps before soliton insertion */
    const char *snap_dir = "data/snapshots";
    const char *tseries_file = NULL;

    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "-mode") && i+1 < argc) mode_str = argv[++i];
        else if (!strcmp(argv[i], "-N") && i+1 < argc) N = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-L") && i+1 < argc) L = atof(argv[++i]);
        else if (!strcmp(argv[i], "-e") && i+1 < argc) e_skyrme = atof(argv[++i]);
        else if (!strcmp(argv[i], "-dt") && i+1 < argc) dt = atof(argv[++i]);
        else if (!strcmp(argv[i], "-T") && i+1 < argc) T_bath = atof(argv[++i]);
        else if (!strcmp(argv[i], "-gamma") && i+1 < argc) gamma_damp = atof(argv[++i]);
        else if (!strcmp(argv[i], "-steps") && i+1 < argc) n_steps = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-snap") && i+1 < argc) snap_interval = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-diag") && i+1 < argc) diag_interval = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-profile") && i+1 < argc) profile_file = argv[++i];
        else if (!strcmp(argv[i], "-D") && i+1 < argc) D_quark = atof(argv[++i]);
        else if (!strcmp(argv[i], "-A") && i+1 < argc) A_quark = atof(argv[++i]);
        else if (!strcmp(argv[i], "-sigma") && i+1 < argc) sigma_quark = atof(argv[++i]);
        else if (!strcmp(argv[i], "-seed") && i+1 < argc) seed = (uint64_t)atol(argv[++i]);
        else if (!strcmp(argv[i], "-Rcut") && i+1 < argc) R_cut_extract = atof(argv[++i]);
        else if (!strcmp(argv[i], "-notest")) run_selftests = 0;
        else if (!strcmp(argv[i], "-therm") && i+1 < argc) n_therm = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-snapdir") && i+1 < argc) snap_dir = argv[++i];
        else if (!strcmp(argv[i], "-tseries") && i+1 < argc) tseries_file = argv[++i];
        else {
            fprintf(stderr, "Usage: %s -mode hedgehog|quarks|b3|extract\n"
                    "  [-N 128] [-L 8] [-dt 0.02] [-T 1.0] [-gamma 0.5]\n"
                    "  [-steps 5000] [-snap 100] [-diag 10] [-profile PATH]\n"
                    "  [-D 2.0] [-A 1.5] [-sigma 1.0] [-seed 42]\n"
                    "  [-notest] [-therm 500] [-snapdir dir] [-tseries file]\n", argv[0]);
            return 1;
        }
    }

    setvbuf(stdout, NULL, _IOLBF, 0);

    /* Parse mode */
    enum Mode mode;
    if (!strcmp(mode_str, "hedgehog")) mode = MODE_HEDGEHOG;
    else if (!strcmp(mode_str, "quarks")) mode = MODE_QUARKS;
    else if (!strcmp(mode_str, "b3")) mode = MODE_B3;
    else if (!strcmp(mode_str, "extract")) mode = MODE_EXTRACT;
    else {
        fprintf(stderr, "Unknown mode: %s\n", mode_str);
        return 1;
    }

    double h = 2.0 * L / N;
    size_t N3 = (size_t)N * N * N;
    double h3 = h * h * h;

    printf("============================================================\n");
    printf(" Proton Field Dynamics — Langevin Simulation\n");
    printf("============================================================\n\n");
    printf("Mode: %s\n", mode_str);
    printf("Grid: N=%d, L=%.1f, h=%.6f\n", N, L, h);
    printf("Parameters: e=%.1f, rho0=%.1f, dt=%.4f\n", e_skyrme, rho0, dt);
    printf("Langevin: T=%.4f, gamma=%.4f\n", T_bath, gamma_damp);
    printf("Steps: %d, snap every %d, diag every %d\n", n_steps, snap_interval, diag_interval);
    printf("CFL: dt*sqrt(3)/h = %.4f (should be < 1)\n", dt*sqrt(3.0)/h);
    printf("Threads: %d\n", omp_get_max_threads());
    printf("Memory: ~%.1f MB\n",
           (double)(N3 * 3 * sizeof(Q4)) / 1e6);
    if (T_bath > 0)
        printf("Noise amplitude: sqrt(2*gamma*T*dt/h^3) = %.4e\n",
               sqrt(2.0 * gamma_damp * T_bath * dt / (h*h*h)));
    printf("RNG seed: %lu\n", (unsigned long)seed);
    printf("\n");

    /* Load profile if needed */
    RadialProfile *prof = NULL;
    if (mode == MODE_HEDGEHOG || mode == MODE_B3 || mode == MODE_EXTRACT) {
        if (!profile_file) profile_file = "data/profiles/profile_sigma_e1.dat";
        prof = load_profile(profile_file);
        if (!prof) { fprintf(stderr, "Failed to load %s\n", profile_file); return 1; }
        printf("Profile: %s (%d points, R_max=%.1f)\n\n", profile_file, prof->n, prof->r_max);
    }

    /* Self-tests */
    if (run_selftests) {
        printf("--- Self-tests ---\n");
        int ok1 = selftest_noise();
        int ok2 = 1, ok3 = 1;
        if (prof) {
            /* Use small grid for speed */
            int N_test = 64;
            double L_test = 6.0;
            ok2 = selftest_force(N_test, L_test, e_skyrme, rho0, prof);
            ok3 = selftest_conservation(N_test, L_test, dt, e_skyrme, rho0, prof);
        }
        if (!ok1 || !ok2 || !ok3)
            printf("WARNING: Some self-tests did not pass cleanly\n");
        printf("\n");
    }

    /* Allocate field, velocity, force */
    Q4 *field = calloc(N3, sizeof(Q4));
    Q4 *vel   = calloc(N3, sizeof(Q4));
    Q4 *force = calloc(N3, sizeof(Q4));

    if (!field || !vel || !force) {
        fprintf(stderr, "Memory allocation failed (need ~%.1f MB)\n",
                (double)(N3 * 3 * sizeof(Q4)) / 1e6);
        return 1;
    }

    /* Initialize per-thread RNGs */
    int max_threads = omp_get_max_threads();
    Rng *rngs = malloc(max_threads * sizeof(Rng));
    for (int t = 0; t < max_threads; t++)
        rng_seed(&rngs[t], seed + (uint64_t)t * 1000003ULL);

    /* Langevin noise amplitude: sqrt(2 * gamma * T * dt / h^3) */
    double noise_amp = (T_bath > 0 && gamma_damp > 0)
                     ? sqrt(2.0 * gamma_damp * T_bath * dt / h3) : 0.0;

    /* ========== PHASE 0: Thermalize vacuum ========== */
    /* Start with uniform vacuum field, run Langevin until T_eff = T_bath.
     * This creates the pre-existing thermal bath that the soliton is
     * immersed into. Periodic BCs (no Dirichlet clamping). */

    /* Initialize vacuum */
    #pragma omp parallel for
    for (size_t ix = 0; ix < N3; ix++) {
        field[ix] = (Q4){rho0, 0, 0, 0};
        vel[ix] = q4_zero();
    }

    if (T_bath > 0 && n_therm > 0) {
        printf("--- Phase 0: Thermalizing vacuum at T=%.4f (%d steps) ---\n",
               T_bath, n_therm);

        /* Seed velocity field with thermal noise at correct temperature.
         * Each tangent component gets variance T/h^3, so |v_comp| ~ sqrt(T/h^3).
         * Generate 4-component noise, then project tangent to S^3. */
        double v_thermal = sqrt(T_bath / h3);
        #pragma omp parallel for
        for (size_t ix = 0; ix < N3; ix++) {
            int tid = omp_get_thread_num();
            Rng *rng = &rngs[tid];
            double g1, g2, g3, g4;
            rng_gauss2(rng, &g1, &g2);
            rng_gauss2(rng, &g3, &g4);
            vel[ix].s  = v_thermal * g1;
            vel[ix].f1 = v_thermal * g2;
            vel[ix].f2 = v_thermal * g3;
            vel[ix].f3 = v_thermal * g4;
            /* Project tangent to field direction */
            double n2 = q4_norm2(field[ix]);
            if (n2 > 1e-20) {
                double vdq = q4_dot(vel[ix], field[ix]) / n2;
                vel[ix] = q4_sub(vel[ix], q4_scale(vdq, field[ix]));
            }
        }
        /* Rescale to exact T_bath */
        velocity_rescale(vel, field, N3, h3, T_bath);

        double Ekin_init = compute_ekin(vel, N3, h3);
        double T_init = compute_teff(Ekin_init, N3);
        printf("  Initial thermal seed: T_eff=%.6e (target %.4e)\n", T_init, T_bath);

        /* Run Langevin thermalization on vacuum */
        for (int step = 1; step <= n_therm; step++) {
            /* Position update */
            #pragma omp parallel for
            for (size_t ix = 0; ix < N3; ix++) {
                field[ix].s  += dt * vel[ix].s;
                field[ix].f1 += dt * vel[ix].f1;
                field[ix].f2 += dt * vel[ix].f2;
                field[ix].f3 += dt * vel[ix].f3;
            }
            /* Sigma-project field and velocity */
            #pragma omp parallel for
            for (size_t ix = 0; ix < N3; ix++) {
                double n = sqrt(q4_norm2(field[ix]));
                if (n > 1e-15) field[ix] = q4_scale(rho0/n, field[ix]);
                double vdq = q4_dot(vel[ix], field[ix]) / (rho0*rho0);
                vel[ix] = q4_sub(vel[ix], q4_scale(vdq, field[ix]));
            }

            /* Force */
            memset(force, 0, N3 * sizeof(Q4));
            compute_skyrme_force(field, N, h, e_skyrme, force);
            #pragma omp parallel for
            for (size_t ix = 0; ix < N3; ix++) {
                double n2 = q4_norm2(field[ix]);
                if (n2 > 1e-20) {
                    double fdq = q4_dot(force[ix], field[ix]) / n2;
                    force[ix] = q4_sub(force[ix], q4_scale(fdq, field[ix]));
                }
            }

            /* Velocity update with Langevin noise + damping */
            #pragma omp parallel for
            for (size_t ix = 0; ix < N3; ix++) {
                int tid = omp_get_thread_num();
                Rng *rng = &rngs[tid];

                vel[ix].s  += dt * force[ix].s;
                vel[ix].f1 += dt * force[ix].f1;
                vel[ix].f2 += dt * force[ix].f2;
                vel[ix].f3 += dt * force[ix].f3;

                if (gamma_damp > 0) {
                    double damp = 1.0 - gamma_damp * dt;
                    if (damp < 0) damp = 0;
                    vel[ix] = q4_scale(damp, vel[ix]);
                }

                if (noise_amp > 0) {
                    double g1, g2, g3, g4;
                    rng_gauss2(rng, &g1, &g2);
                    rng_gauss2(rng, &g3, &g4);
                    vel[ix].s  += noise_amp * g1;
                    vel[ix].f1 += noise_amp * g2;
                    vel[ix].f2 += noise_amp * g3;
                    vel[ix].f3 += noise_amp * g4;
                }

                double n2 = q4_norm2(field[ix]);
                if (n2 > 1e-20) {
                    double vdq = q4_dot(vel[ix], field[ix]) / n2;
                    vel[ix] = q4_sub(vel[ix], q4_scale(vdq, field[ix]));
                }
            }

            /* Velocity rescale: enforce T_eff = T_bath exactly */
            velocity_rescale(vel, field, N3, h3, T_bath);

            if (step % 100 == 0 || step == n_therm) {
                double Ekin = compute_ekin(vel, N3, h3);
                double Teff = compute_teff(Ekin, N3);
                EnergyResult Eth = compute_energy(field, N, h, e_skyrme);
                printf("  therm step %4d: T_eff=%.6e E2=%.2f E4=%.2f\n",
                       step, Teff, Eth.E2, Eth.E4);
            }
        }

        double Ekin_final = compute_ekin(vel, N3, h3);
        double T_final = compute_teff(Ekin_final, N3);
        printf("  Thermalization complete: T_eff=%.6e (target %.4e)\n\n", T_final, T_bath);
    } else if (T_bath <= 0) {
        printf("T=0: no thermalization needed.\n\n");
    }

    /* ========== PHASE 1: Insert soliton into thermal field ========== */
    /* Multiply the thermalized field by the soliton ansatz (product ansatz).
     * This preserves the thermal fluctuations while adding the soliton. */
    printf("--- Phase 1: Insert soliton ---\n");

    if (T_bath > 0) {
        /* Save thermal field, then compose with soliton via product ansatz */
        Q4 *q_soliton = calloc(N3, sizeof(Q4));

        switch (mode) {
        case MODE_HEDGEHOG:
            init_hedgehog(q_soliton, vel /* dummy */, N, L, h, rho0, prof);
            printf("Inserting hedgehog B=1 into thermal bath\n");
            break;
        case MODE_QUARKS:
            init_three_quarks(q_soliton, vel /* dummy */, N, L, h, rho0,
                              D_quark, A_quark, sigma_quark);
            printf("Inserting three quarks (D=%.1f, A=%.1f) into thermal bath\n",
                   D_quark, A_quark);
            break;
        case MODE_B3:
            init_b3_ratmap(q_soliton, vel /* dummy */, N, L, h, rho0, prof);
            printf("Inserting B=3 rational map into thermal bath\n");
            break;
        case MODE_EXTRACT:
            init_quark_extract(q_soliton, vel /* dummy */, N, L, h, rho0, prof, R_cut_extract);
            printf("Inserting extracted quark into thermal bath\n");
            break;
        }

        /* Product ansatz: field = q_soliton * q_thermal / rho0
         * This composes the soliton rotation with the thermal fluctuations. */
        #pragma omp parallel for
        for (size_t ix = 0; ix < N3; ix++) {
            Q4 q_th = field[ix];  /* thermal field */
            Q4 q_sol = q_soliton[ix];  /* soliton */
            field[ix] = q4_scale(1.0/rho0, q4_mul(q_sol, q_th));
            /* Re-normalize to sigma model */
            double n = sqrt(q4_norm2(field[ix]));
            if (n > 1e-15) field[ix] = q4_scale(rho0/n, field[ix]);
            /* Re-project velocity tangent to new field direction */
            double n2 = q4_norm2(field[ix]);
            if (n2 > 1e-20) {
                double vdq = q4_dot(vel[ix], field[ix]) / n2;
                vel[ix] = q4_sub(vel[ix], q4_scale(vdq, field[ix]));
            }
        }

        /* Rescale velocity back to T_bath after re-projection */
        velocity_rescale(vel, field, N3, h3, T_bath);
        free(q_soliton);
    } else {
        /* T=0: just initialize directly (no thermal field) */
        switch (mode) {
        case MODE_HEDGEHOG:
            init_hedgehog(field, vel, N, L, h, rho0, prof);
            printf("Hedgehog B=1 initialized (T=0)\n");
            break;
        case MODE_QUARKS:
            init_three_quarks(field, vel, N, L, h, rho0, D_quark, A_quark, sigma_quark);
            printf("Three quarks: D=%.1f, A=%.1f, sigma=%.1f (T=0)\n",
                   D_quark, A_quark, sigma_quark);
            break;
        case MODE_B3:
            init_b3_ratmap(field, vel, N, L, h, rho0, prof);
            printf("B=3 rational map initialized (T=0)\n");
            break;
        case MODE_EXTRACT:
            init_quark_extract(field, vel, N, L, h, rho0, prof, R_cut_extract);
            printf("Single quark extracted from B=3 (T=0)\n");
            break;
        }
    }

    /* Initial diagnostics */
    EnergyResult E0 = compute_energy(field, N, h, e_skyrme);
    double Q0 = compute_Q(field, N, h);
    double Ekin0 = compute_ekin(vel, N3, h3);
    double Teff0 = compute_teff(Ekin0, N3);
    printf("After insertion: E2=%.4f, E4=%.4f, E=%.4f, Q=%.4f, T_eff=%.4e\n\n",
           E0.E2, E0.E4, E0.E2+E0.E4, Q0, Teff0);

    /* Open time series file */
    FILE *ts_fp = NULL;
    if (tseries_file) {
        ts_fp = fopen(tseries_file, "w");
    } else {
        /* Default name based on mode and temperature */
        char ts_name[256];
        snprintf(ts_name, sizeof(ts_name), "data/tseries_%s_T%.2f.dat",
                 mode_str, T_bath);
        ts_fp = fopen(ts_name, "w");
        if (!ts_fp) {
            /* Try without data/ prefix */
            snprintf(ts_name, sizeof(ts_name), "tseries_%s_T%.2f.dat",
                     mode_str, T_bath);
            ts_fp = fopen(ts_name, "w");
        }
    }
    if (ts_fp) {
        fprintf(ts_fp, "# proton.c time series: mode=%s T=%.4f gamma=%.4f N=%d L=%.1f e=%.1f\n",
                mode_str, T_bath, gamma_damp, N, L, e_skyrme);
        fprintf(ts_fp, "# step  time  E2  E4  E_kin  Q  T_eff\n");
    }

    /* Write initial snapshot */
    {
        char snap_name[256];
        snprintf(snap_name, sizeof(snap_name), "%s/snapshot_%04d.bin", snap_dir, 0);
        /* Try to create directory (may fail if exists, that's ok) */
        char mkdir_cmd[300];
        snprintf(mkdir_cmd, sizeof(mkdir_cmd), "mkdir -p %s", snap_dir);
        if (system(mkdir_cmd) != 0) { /* ignore */ }
        write_snapshot(snap_name, field, N, L, 0, 0.0);
        printf("Wrote %s\n", snap_name);
    }

    /* ========== Leapfrog + Langevin time loop ========== */
    printf("\n--- Time evolution ---\n");
    printf("%6s %8s %10s %10s %10s %8s %10s\n",
           "step", "time", "E2", "E4", "E_kin", "Q", "T_eff");

    /* Step 0: compute force, do leapfrog half-step for velocity */
    memset(force, 0, N3 * sizeof(Q4));
    compute_skyrme_force(field, N, h, e_skyrme, force);

    /* Project force tangent */
    #pragma omp parallel for
    for (size_t ix = 0; ix < N3; ix++) {
        double n2 = q4_norm2(field[ix]);
        if (n2 > 1e-20) {
            double fdq = q4_dot(force[ix], field[ix]) / n2;
            force[ix] = q4_sub(force[ix], q4_scale(fdq, field[ix]));
        }
    }

    /* Half-step velocity: v(dt/2) = v(0) + (dt/2)*F - (gamma*dt/2)*v + noise_half */
    #pragma omp parallel for
    for (size_t ix = 0; ix < N3; ix++) {
        vel[ix].s  += 0.5*dt * force[ix].s;
        vel[ix].f1 += 0.5*dt * force[ix].f1;
        vel[ix].f2 += 0.5*dt * force[ix].f2;
        vel[ix].f3 += 0.5*dt * force[ix].f3;

        /* Langevin damping (half-step) */
        if (gamma_damp > 0) {
            double damp = 1.0 - 0.5*gamma_damp*dt;
            if (damp < 0) damp = 0;
            vel[ix] = q4_scale(damp, vel[ix]);
        }

        /* Project velocity tangent */
        double n2 = q4_norm2(field[ix]);
        if (n2 > 1e-20) {
            double vdq = q4_dot(vel[ix], field[ix]) / n2;
            vel[ix] = q4_sub(vel[ix], q4_scale(vdq, field[ix]));
        }
    }

    for (int step = 1; step <= n_steps; step++) {
        double t_now = step * dt;

        /* 1. Position update: q(t+dt) = q(t) + dt * v(t+dt/2) */
        #pragma omp parallel for
        for (size_t ix = 0; ix < N3; ix++) {
            field[ix].s  += dt * vel[ix].s;
            field[ix].f1 += dt * vel[ix].f1;
            field[ix].f2 += dt * vel[ix].f2;
            field[ix].f3 += dt * vel[ix].f3;
        }

        /* 2. Project |q| -> rho0, project v tangential */
        #pragma omp parallel for
        for (size_t ix = 0; ix < N3; ix++) {
            double n = sqrt(q4_norm2(field[ix]));
            if (n > 1e-15) field[ix] = q4_scale(rho0/n, field[ix]);
            double vdq = q4_dot(vel[ix], field[ix]) / (rho0*rho0);
            vel[ix] = q4_sub(vel[ix], q4_scale(vdq, field[ix]));
        }

        /* 3. Compute force */
        memset(force, 0, N3 * sizeof(Q4));
        compute_skyrme_force(field, N, h, e_skyrme, force);

        /* 5. Project force tangential */
        #pragma omp parallel for
        for (size_t ix = 0; ix < N3; ix++) {
            double n2 = q4_norm2(field[ix]);
            if (n2 > 1e-20) {
                double fdq = q4_dot(force[ix], field[ix]) / n2;
                force[ix] = q4_sub(force[ix], q4_scale(fdq, field[ix]));
            }
        }

        /* 6. Velocity update: v(t+3dt/2) = v(t+dt/2) + dt*F - gamma*dt*v + noise */
        #pragma omp parallel for
        for (size_t ix = 0; ix < N3; ix++) {
            int tid = omp_get_thread_num();
            Rng *rng = &rngs[tid];

            /* Deterministic: force + damping */
            vel[ix].s  += dt * force[ix].s;
            vel[ix].f1 += dt * force[ix].f1;
            vel[ix].f2 += dt * force[ix].f2;
            vel[ix].f3 += dt * force[ix].f3;

            /* Langevin damping */
            if (gamma_damp > 0) {
                double damp = 1.0 - gamma_damp * dt;
                if (damp < 0) damp = 0;
                vel[ix] = q4_scale(damp, vel[ix]);
            }

            /* Langevin noise */
            if (noise_amp > 0) {
                double g1, g2, g3, g4;
                rng_gauss2(rng, &g1, &g2);
                rng_gauss2(rng, &g3, &g4);
                vel[ix].s  += noise_amp * g1;
                vel[ix].f1 += noise_amp * g2;
                vel[ix].f2 += noise_amp * g3;
                vel[ix].f3 += noise_amp * g4;
            }

            /* 7. Project velocity tangential */
            double n2 = q4_norm2(field[ix]);
            if (n2 > 1e-20) {
                double vdq = q4_dot(vel[ix], field[ix]) / n2;
                vel[ix] = q4_sub(vel[ix], q4_scale(vdq, field[ix]));
            }
        }

        /* 8. Velocity rescale: enforce T_eff = T_bath exactly */
        velocity_rescale(vel, field, N3, h3, T_bath);

        /* Diagnostics */
        if (step % diag_interval == 0 || step == 1) {
            EnergyResult E = compute_energy(field, N, h, e_skyrme);
            double Q = compute_Q(field, N, h);

            double Ekin = compute_ekin(vel, N3, h3);
            double T_eff = compute_teff(Ekin, N3);

            printf("%6d %8.3f %10.4f %10.4f %10.4f %8.4f %10.4e\n",
                   step, t_now, E.E2, E.E4, Ekin, Q, T_eff);

            if (ts_fp) {
                fprintf(ts_fp, "%d %.6f %.6f %.6f %.6f %.6f %.6e\n",
                        step, t_now, E.E2, E.E4, Ekin, Q, T_eff);
                fflush(ts_fp);
            }

            if (isnan(E.E2) || isnan(E.E4)) {
                printf("NaN detected — stopping.\n");
                break;
            }
        }

        /* Snapshots */
        if (snap_interval > 0 && step % snap_interval == 0) {
            char snap_name[256];
            snprintf(snap_name, sizeof(snap_name), "%s/snapshot_%04d.bin",
                     snap_dir, step);
            write_snapshot(snap_name, field, N, L, step, t_now);
        }
    }

    printf("\n--- Done ---\n");
    {
        EnergyResult E_final = compute_energy(field, N, h, e_skyrme);
        double Q_final = compute_Q(field, N, h);
        printf("Final: E2=%.4f, E4=%.4f, E=%.4f, Q=%.4f\n",
               E_final.E2, E_final.E4, E_final.E2+E_final.E4, Q_final);
    }

    /* Clean up */
    if (ts_fp) fclose(ts_fp);
    free(field);
    free(vel);
    free(force);
    free(rngs);
    if (prof) free_profile(prof);

    return 0;
}
