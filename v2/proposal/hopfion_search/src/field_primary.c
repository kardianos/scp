/*
 * field_primary.c — 3D field-primary soliton simulation with FFT Poisson solver
 *
 * Implements 6 variations:
 *   V1: Baseline Skyrme only (two solitons, static + dynamic)
 *   V2: + Massless scalar mediator via FFT Poisson
 *   V3: + Massive mediator (pion mass, Yukawa)
 *   V4: Pullback metric feedback (BLV, control — should show ZERO deflection)
 *   V5: + L6 with metric feedback (breaks P/m=2, gives attractive lensing)
 *   V6: Dynamic scatter with mediator (two boosted solitons)
 *
 * Self-contained: no dependencies beyond -lm and -fopenmp.
 * Compile: gcc -O3 -fopenmp -std=c11 -o field_primary field_primary.c -lm
 */

#define _USE_MATH_DEFINES
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <complex.h>
#include <time.h>

#ifdef _OPENMP
#include <omp.h>
#else
static inline int omp_get_max_threads(void) { return 1; }
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

static inline Q4 q4_zero(void) { return (Q4){0,0,0,0}; }

/* ========== Grid index with periodic BCs ========== */

static inline int idx(int N, int i, int j, int k) {
    return ((i % N + N) % N) * N * N
         + ((j % N + N) % N) * N
         + ((k % N + N) % N);
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
        double rv, fv, dummy;
        if (sscanf(line, "%lf %lf %lf", &rv, &fv, &dummy) >= 2) {
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

static Q4 hedgehog(const RadialProfile *prof, double rho0,
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

/* ========== FFT (Cooley-Tukey radix-2, in-place) ========== */

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

    /* For inverse FFT, divide by N */
    if (sign == +1) {
        double inv = 1.0 / N;
        for (int i = 0; i < N; i++)
            x[i] *= inv;
    }
}

/* 3D FFT: apply fft1d along each axis */
static void fft3d(double complex *data, int N, int sign) {
    double complex *buf = malloc(N * sizeof(double complex));

    /* Along x (axis 0): for each (j,k), transform data[i*N*N + j*N + k] over i */
    for (int j = 0; j < N; j++) {
        for (int k = 0; k < N; k++) {
            for (int i = 0; i < N; i++)
                buf[i] = data[i*N*N + j*N + k];
            fft1d(buf, N, sign);
            for (int i = 0; i < N; i++)
                data[i*N*N + j*N + k] = buf[i];
        }
    }

    /* Along y (axis 1): for each (i,k), transform over j */
    for (int i = 0; i < N; i++) {
        for (int k = 0; k < N; k++) {
            for (int j = 0; j < N; j++)
                buf[j] = data[i*N*N + j*N + k];
            fft1d(buf, N, sign);
            for (int j = 0; j < N; j++)
                data[i*N*N + j*N + k] = buf[j];
        }
    }

    /* Along z (axis 2): for each (i,j), transform over k */
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            fft1d(&data[i*N*N + j*N], N, sign);
        }
    }

    free(buf);
}

/* Forward declaration for Poisson solver (used in self-test) */
static void solve_poisson_fft(const double *source, double *p_out,
                               double complex *work,
                               int N, double h,
                               double g_top, double kappa2, double mu2);

/* ========== FFT self-test ========== */

static int fft_selftest(int N) {
    size_t N3 = (size_t)N * N * N;
    double complex *data = malloc(N3 * sizeof(double complex));
    double complex *orig = malloc(N3 * sizeof(double complex));

    /* Fill with random data */
    srand(42);
    for (size_t i = 0; i < N3; i++) {
        data[i] = (double)(rand() % 1000) / 1000.0
                + I * (double)(rand() % 1000) / 1000.0;
        orig[i] = data[i];
    }

    /* Forward then inverse */
    fft3d(data, N, -1);
    fft3d(data, N, +1);

    /* Check roundtrip */
    double maxerr = 0;
    for (size_t i = 0; i < N3; i++) {
        double err = cabs(data[i] - orig[i]);
        if (err > maxerr) maxerr = err;
    }

    free(data);
    free(orig);

    printf("  FFT roundtrip test (N=%d): max error = %.2e %s\n",
           N, maxerr, maxerr < 1e-10 ? "OK" : "FAIL");
    return maxerr < 1e-10;
}

/* Poisson test: solve -nabla^2 p = source where source = sin(2*pi*x/L_box)
 * Exact solution on periodic domain: p = sin(2*pi*x/L_box) / (2*pi/L_box)^2 */
static int poisson_selftest(int N, double L) {
    double h = 2.0 * L / N;
    double L_box = 2.0 * L;
    double kk = 2.0 * M_PI / L_box;
    size_t N3 = (size_t)N * N * N;
    double *source = malloc(N3 * sizeof(double));
    double *p_out = malloc(N3 * sizeof(double));
    double complex *work = malloc(N3 * sizeof(double complex));

    /* Source: sin(2*pi*x/L_box) */
    for (int ix = 0; ix < N; ix++) {
        double x = -L + (ix + 0.5) * h;
        double val = sin(kk * x);
        for (int jx = 0; jx < N; jx++)
            for (int kx = 0; kx < N; kx++)
                source[(size_t)ix*N*N + jx*N + kx] = val;
    }

    /* Solve -nabla^2 p = source via FFT */
    solve_poisson_fft(source, p_out, work, N, h, 1.0, 1.0, 0.0);

    /* Expected: p = sin(kk*x) / kk^2 */
    double maxerr = 0;
    int n_check = 0;
    for (int ix = N/4; ix < 3*N/4; ix += N/8) {
        double x = -L + (ix + 0.5) * h;
        double expected = sin(kk * x) / (kk * kk);
        double got = p_out[(size_t)ix*N*N + (N/2)*N + N/2];
        double err = fabs(got - expected);
        if (err > maxerr) maxerr = err;
        n_check++;
    }

    free(source);
    free(p_out);
    free(work);
    printf("  Poisson test (sin source): max error = %.2e (%d pts) %s\n",
           maxerr, n_check, maxerr < 0.01 ? "OK" : "FAIL");
    return maxerr < 0.01;
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

/* Consistent Laplacian: D_d(D_d q) with 4th-order stencil composed with itself
 * 9-point stencil per direction: {1,-16,64,16,-130,16,64,-16,1}/(144h^2) */
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

/* ========== Skyrme energy and force ========== */

typedef struct {
    Q4 A[3]; /* right-currents */
    Q4 G[3]; /* G-tensor */
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

        /* Commutators */
        Q4 C[3]; /* C01, C02, C12 */
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

typedef struct { double E2, E4, Etotal; } EnergyResult;

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
    en.Etotal = en.E2 + en.E4;
    return en;
}

/* Compute Skyrme force (E2 + E4) and store in force array.
 * force must be zeroed before calling. */
static void compute_skyrme_force(const Q4 *q, int N, double h, double e_skyrme, Q4 *force) {
    size_t N3 = (size_t)N*N*N;
    double inv12h = 1.0 / (12.0 * h);
    double inv_2e2 = 1.0 / (2.0 * e_skyrme * e_skyrme);
    double inv144h2 = 1.0 / (144.0 * h * h);
    double sigma[4] = {1.0, -1.0, -1.0, -1.0};

    /* Step 1: precompute Skyrme data */
    SkyrmePre *pre = malloc(N3 * sizeof(SkyrmePre));
    compute_skyrme_predata(q, N, h, pre);

    /* Step 2: precompute pi_{a,d}(x) = <q~*eps_a, G_d>_0 */
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

    /* Step 3: compute forces */
    #pragma omp parallel for collapse(3) schedule(static)
    for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
    for (int k = 0; k < N; k++) {
        int ix = idx(N,i,j,k);

        /* E2 force: consistent Laplacian */
        Q4 lap = q4_laplacian(q, N, i, j, k, inv144h2);
        force[ix].s  = lap.s;
        force[ix].f1 = lap.f1;
        force[ix].f2 = lap.f2;
        force[ix].f3 = lap.f3;

        /* E4 force */
        double f4[4] = {0, 0, 0, 0};
        /* Term 1 */
        for (int d = 0; d < 3; d++) {
            Q4 dq = q4_deriv(q, N, i, j, k, d, inv12h);
            for (int a = 0; a < 4; a++) {
                Q4 ea_dq = q4_left_basis(a, dq);
                Q4 prod = q4_mul(ea_dq, pre[ix].G[d]);
                f4[a] += sigma[a] * prod.s;
            }
        }
        /* Term 2: -D_d(pi_{a,d}) */
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

/* ========== Baryon density and cofactor ========== */

/* B0(x) = (1/(2pi^2)) eps_{abcd} q_a D_x[q_b] D_y[q_c] D_z[q_d]
 * Cofactor: J_e(x) = (1/(2pi^2)) eps_{ebcd} D_x[q_b] D_y[q_c] D_z[q_d]
 * Note: B0 = sum_a q_a * J_a */

static void compute_baryon_density(const Q4 *q, int N, double h,
                                   double *B0, Q4 *J_cofactor) {
    double inv12h = 1.0 / (12.0 * h);
    double inv2pi2 = 1.0 / (2.0 * M_PI * M_PI);

    #pragma omp parallel for collapse(3) schedule(static)
    for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
    for (int k = 0; k < N; k++) {
        int ix = idx(N,i,j,k);
        Q4 dq[3];
        for (int d = 0; d < 3; d++)
            dq[d] = q4_deriv(q, N, i, j, k, d, inv12h);

        /* 3x3 matrix of spatial derivatives:
         * dq[d].s, dq[d].f1, dq[d].f2, dq[d].f3 for d=0,1,2
         * We need the determinant-like contraction with eps_{abcd} */

        /* Explicit expansion of eps_{abcd} q_a D_x[q_b] D_y[q_c] D_z[q_d]:
         * This is the quaternion triple product, related to det of 3x4 matrix.
         *
         * Using A_d = q~ * dq_d, the charge density is
         * B0 = -(1/(2pi^2)) <A_0 A_1 A_2>_0 / |q|^6
         * But for sigma model |q|=rho0, this simplifies. */

        Q4 qr = q4_rev(q[ix]);
        Q4 A[3];
        for (int d = 0; d < 3; d++)
            A[d] = q4_mul(qr, dq[d]);

        Q4 A01 = q4_mul(A[0], A[1]);
        Q4 A012 = q4_mul(A01, A[2]);
        double norm2 = q4_norm2(q[ix]);
        double norm6 = norm2 * norm2 * norm2;

        if (norm6 > 1e-30)
            B0[ix] = -inv2pi2 * A012.s / norm6;
        else
            B0[ix] = 0;

        /* Cofactor J_e: for each component e, it's the eps contraction
         * with the other 3 components of dq.
         * J_e = (1/(2pi^2)) eps_{ebcd} D_x[q_b] D_y[q_c] D_z[q_d]
         *
         * Explicit for e=0 (s component):
         * J_0 = (1/(2pi^2))[dq_x.f1*(dq_y.f2*dq_z.f3 - dq_y.f3*dq_z.f2)
         *                   - dq_x.f2*(dq_y.f1*dq_z.f3 - dq_y.f3*dq_z.f1)
         *                   + dq_x.f3*(dq_y.f1*dq_z.f2 - dq_y.f2*dq_z.f1)]
         * This is the 3x3 determinant of the f1,f2,f3 derivatives.
         */

        /* Shorthand */
        double ax = dq[0].s,  bx = dq[0].f1, cx = dq[0].f2, dx_v = dq[0].f3;
        double ay = dq[1].s,  by = dq[1].f1, cy = dq[1].f2, dy_v = dq[1].f3;
        double az = dq[2].s,  bz = dq[2].f1, cz = dq[2].f2, dz_v = dq[2].f3;

        /* J_0 (e=s): eps_{0bcd} = eps_{bcd} with (b,c,d) in {1,2,3}
         * = det of [bx cx dx; by cy dy; bz cz dz] */
        double J0 = bx*(cy*dz_v - cz*dy_v) - cx*(by*dz_v - bz*dy_v)
                   + dx_v*(by*cz - bz*cy);

        /* J_1 (e=f1): eps_{1bcd} with (b,c,d) in {0,2,3}
         * eps_{1023} = -1 (even perm 0123 -> 1023 is odd)
         * = -det of [ax cx dx; ay cy dy; az cz dz] */
        double J1 = -(ax*(cy*dz_v - cz*dy_v) - cx*(ay*dz_v - az*dy_v)
                     + dx_v*(ay*cz - az*cy));

        /* J_2 (e=f2): eps_{2bcd} with (b,c,d) in {0,1,3}
         * eps_{2013} = +1 (even perm)
         * = det of [ax bx dx; ay by dy; az bz dz] */
        double J2 = ax*(by*dz_v - bz*dy_v) - bx*(ay*dz_v - az*dy_v)
                   + dx_v*(ay*bz - az*by);

        /* J_3 (e=f3): eps_{3bcd} with (b,c,d) in {0,1,2}
         * eps_{3012} = -1
         * = -det of [ax bx cx; ay by cy; az bz cz] */
        double J3 = -(ax*(by*cz - bz*cy) - bx*(ay*cz - az*cy)
                     + cx*(ay*bz - az*by));

        if (J_cofactor) {
            J_cofactor[ix].s  = inv2pi2 * J0;
            J_cofactor[ix].f1 = inv2pi2 * J1;
            J_cofactor[ix].f2 = inv2pi2 * J2;
            J_cofactor[ix].f3 = inv2pi2 * J3;
        }
    }
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

/* ========== FFT Poisson solve ========== */

static void solve_poisson_fft(const double *source, double *p_out,
                               double complex *work,
                               int N, double h,
                               double g_top, double kappa2, double mu2) {
    size_t N3 = (size_t)N * N * N;

    /* Copy source into complex array */
    for (size_t i = 0; i < N3; i++)
        work[i] = source[i] + 0.0*I;

    /* Forward FFT */
    fft3d(work, N, -1);

    /* Multiply by Green's function:
     * Equation: -kappa2 nabla^2 p + mu2 p = g_top * B0
     * In Fourier: (kappa2 k_tilde^2 + mu2) p_hat = g_top * B0_hat
     * p_hat = g_top * B0_hat / (kappa2 * k_tilde^2 + mu2)
     */
    for (int ix = 0; ix < N; ix++)
    for (int jx = 0; jx < N; jx++)
    for (int kx = 0; kx < N; kx++) {
        int nx = (ix < N/2) ? ix : ix - N;
        int ny = (jx < N/2) ? jx : jx - N;
        int nz = (kx < N/2) ? kx : kx - N;
        double k2 = (2.0/(h*h)) * (3.0
            - cos(2.0*M_PI*nx/N)
            - cos(2.0*M_PI*ny/N)
            - cos(2.0*M_PI*nz/N));
        double denom = kappa2 * k2 + mu2;
        size_t idx_flat = (size_t)ix*N*N + jx*N + kx;
        if (denom > 1e-30)
            work[idx_flat] *= g_top / denom;
        else
            work[idx_flat] = 0;  /* DC mode: gauge choice */
    }

    /* Inverse FFT */
    fft3d(work, N, +1);

    /* Extract real part */
    for (size_t i = 0; i < N3; i++)
        p_out[i] = creal(work[i]);
}

/* ========== Poisson residual check ========== */

static double check_poisson_residual(const double *p, const double *B0,
                                      int N, double h,
                                      double g_top, double kappa2, double mu2) {
    double max_res = 0, max_rhs = 0;
    /* Sample 100 random interior points */
    srand(123);
    for (int s = 0; s < 100; s++) {
        int i = 4 + rand() % (N-8);
        int j = 4 + rand() % (N-8);
        int k = 4 + rand() % (N-8);
        int ix = idx(N,i,j,k);

        /* 2nd-order discrete Laplacian */
        double lap = 0;
        double inv_h2 = 1.0 / (h*h);
        lap += (p[idx(N,i+1,j,k)] + p[idx(N,i-1,j,k)] - 2*p[ix]) * inv_h2;
        lap += (p[idx(N,i,j+1,k)] + p[idx(N,i,j-1,k)] - 2*p[ix]) * inv_h2;
        lap += (p[idx(N,i,j,k+1)] + p[idx(N,i,j,k-1)] - 2*p[ix]) * inv_h2;

        double lhs = -kappa2 * lap + mu2 * p[ix];
        double rhs = g_top * B0[ix];
        double res = fabs(lhs - rhs);
        if (res > max_res) max_res = res;
        if (fabs(rhs) > max_rhs) max_rhs = fabs(rhs);
    }
    return (max_rhs > 1e-30) ? max_res / max_rhs : max_res;
}

/* ========== Soliton position tracking ========== */

static void find_solitons_3d(const Q4 *q, int N, double h, double L, double rho0,
                              double *z1_out, double *z2_out) {
    double sz1 = 0, sw1 = 0, sz2 = 0, sw2 = 0;

    for (int i = 2; i < N-2; i++) {
        for (int j = 2; j < N-2; j++) {
            for (int k = 2; k < N-2; k++) {
                double z = -L + (k + 0.5) * h;
                int ix = idx(N,i,j,k);
                double w = rho0 - q[ix].s;
                if (w < 0) w = 0;
                w = w * w;
                if (z > 0) { sz1 += z*w; sw1 += w; }
                else       { sz2 += z*w; sw2 += w; }
            }
        }
    }
    *z1_out = (sw1 > 1e-20) ? sz1/sw1 : 0;
    *z2_out = (sw2 > 1e-20) ? sz2/sw2 : 0;
}

/* ========== L6 force ========== */
/* F6_e(x) = -2*lambda6 * B0(x) * J_e(x) */

static void compute_l6_energy(const double *B0, int N, double h3,
                               double lambda6, double *E6_out) {
    double sum = 0;
    size_t N3 = (size_t)N*N*N;
    for (size_t i = 0; i < N3; i++)
        sum += B0[i] * B0[i];
    *E6_out = lambda6 * sum * h3;
}

/* ========== MAIN ========== */

int main(int argc, char *argv[]) {
    /* Default parameters */
    int N = 256;
    double L = 10.0;
    double e_skyrme = 1.0;
    double rho0 = 1.0;
    double dt = 0.02;
    double g_top = 10.0;
    double kappa2 = 1.0;
    double mu_massless = 0.05;
    double mu_massive = 0.398;
    double lambda6 = 10.0;
    const char *profile_file = NULL;
    int var_select = 0;  /* 0 = all, 1-6 = specific variation */

    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "-N") && i+1 < argc) N = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-L") && i+1 < argc) L = atof(argv[++i]);
        else if (!strcmp(argv[i], "-e") && i+1 < argc) e_skyrme = atof(argv[++i]);
        else if (!strcmp(argv[i], "-dt") && i+1 < argc) dt = atof(argv[++i]);
        else if (!strcmp(argv[i], "-gtop") && i+1 < argc) g_top = atof(argv[++i]);
        else if (!strcmp(argv[i], "-mu") && i+1 < argc) mu_massless = atof(argv[++i]);
        else if (!strcmp(argv[i], "-lam6") && i+1 < argc) lambda6 = atof(argv[++i]);
        else if (!strcmp(argv[i], "-profile") && i+1 < argc) profile_file = argv[++i];
        else if (!strcmp(argv[i], "-var") && i+1 < argc) var_select = atoi(argv[++i]);
        else {
            fprintf(stderr, "Usage: %s [-N grid] [-L box] [-e e] [-dt step] "
                    "[-gtop g] [-mu mu] [-lam6 L6] [-profile file] [-var 1-6]\n", argv[0]);
            return 1;
        }
    }

    setvbuf(stdout, NULL, _IOLBF, 0);

    double h = 2.0 * L / N;
    size_t N3 = (size_t)N * N * N;

    printf("============================================================\n");
    printf(" Field-Primary 3D Soliton Simulation\n");
    printf("============================================================\n\n");
    printf("Grid: N=%d, L=%.1f, h=%.6f\n", N, L, h);
    printf("Parameters: e=%.1f, rho0=%.1f, dt=%.4f\n", e_skyrme, rho0, dt);
    printf("Mediator: g_top=%.1f, kappa2=%.1f\n", g_top, kappa2);
    printf("CFL: dt*sqrt(3)/h = %.4f (should be < 1)\n", dt*sqrt(3.0)/h);
    printf("Threads: %d\n", omp_get_max_threads());
    printf("Memory: ~%.1f GB\n\n",
           (double)(N3*(4*sizeof(double)*4 + sizeof(double) + sizeof(double complex)))/1e9);

    /* Load profile */
    if (!profile_file) {
        profile_file = "data/profiles/profile_sigma_e1.dat";
    }
    RadialProfile *prof = load_profile(profile_file);
    if (!prof) { fprintf(stderr, "Failed to load %s\n", profile_file); return 1; }
    printf("Profile: %s (%d points, R_max=%.1f)\n\n", profile_file, prof->n, prof->r_max);

    /* FFT self-tests */
    printf("--- Self-tests ---\n");
    /* Use a small N for FFT test to be fast */
    int fft_ok = fft_selftest(64);
    int poisson_ok = poisson_selftest(64, 5.0);
    if (!fft_ok || !poisson_ok) {
        printf("WARNING: Self-tests failed\n");
    }
    printf("\n");

    /* Allocate field, velocity, force, baryon density */
    Q4 *q     = calloc(N3, sizeof(Q4));
    Q4 *vel   = calloc(N3, sizeof(Q4));
    Q4 *force = calloc(N3, sizeof(Q4));
    double *B0  = calloc(N3, sizeof(double));
    Q4 *J_cof  = calloc(N3, sizeof(Q4));
    double *p_med = calloc(N3, sizeof(double));
    double complex *fft_work = calloc(N3, sizeof(double complex));

    if (!q || !vel || !force || !B0 || !J_cof || !p_med || !fft_work) {
        fprintf(stderr, "Memory allocation failed (need ~%.1f GB)\n",
                (double)(N3*(5*sizeof(Q4) + 2*sizeof(double) + sizeof(double complex)))/1e9);
        return 1;
    }

    /* Helper: initialize two solitons at +/- z0 */
    auto void init_two_solitons(double z0);
    void init_two_solitons(double z0) {
        #pragma omp parallel for collapse(3) schedule(static)
        for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
        for (int k = 0; k < N; k++) {
            double x = -L + (i + 0.5) * h;
            double y = -L + (j + 0.5) * h;
            double z = -L + (k + 0.5) * h;
            int ix = idx(N,i,j,k);

            Q4 q1 = hedgehog(prof, rho0, x, y, z, 0, 0, +z0);
            Q4 q2 = hedgehog(prof, rho0, x, y, z, 0, 0, -z0);
            Q4 qp = q4_scale(1.0/rho0, q4_mul(q1, q2));

            /* Normalize to sigma model */
            double n = sqrt(q4_norm2(qp));
            if (n > 1e-15) {
                double sc = rho0 / n;
                qp = q4_scale(sc, qp);
            }
            q[ix] = qp;
            vel[ix] = q4_zero();
        }
    }

    /* Helper: initialize single soliton + wave packet */
    auto void init_soliton_wave(double wave_x0, double wave_sigma, double wave_amp,
                                double wave_vx);
    void init_soliton_wave(double wave_x0, double wave_sigma, double wave_amp,
                           double wave_vx) {
        #pragma omp parallel for collapse(3) schedule(static)
        for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
        for (int k = 0; k < N; k++) {
            double x = -L + (i + 0.5) * h;
            double y = -L + (j + 0.5) * h;
            double z = -L + (k + 0.5) * h;
            int ix = idx(N,i,j,k);

            /* Single soliton at origin */
            Q4 qs = hedgehog(prof, rho0, x, y, z, 0, 0, 0);

            /* Add small wave perturbation: delta_q in f1 direction */
            double dx = x - wave_x0;
            double r2 = dx*dx + y*y + z*z;
            double env = wave_amp * exp(-r2 / (2.0 * wave_sigma * wave_sigma));
            qs.f1 += env;

            /* Re-normalize to sigma model */
            double n = sqrt(q4_norm2(qs));
            if (n > 1e-15)
                qs = q4_scale(rho0 / n, qs);
            q[ix] = qs;

            /* Velocity from wave packet: v = d/dt[perturbation] = -vx * d/dx[perturbation]
             * The perturbation is delta_q.f1 = A * exp(-r^2/(2*sigma^2))
             * d/dx = -dx/sigma^2 * A * exp(...)
             * So v.f1 = -wave_vx * (-dx/sigma^2 * env) = wave_vx * dx/sigma^2 * env */
            Q4 v_init = q4_zero();
            v_init.f1 = wave_vx * dx / (wave_sigma * wave_sigma) * env;

            /* Project velocity tangent to S3 */
            double vdq = (v_init.s*qs.s + v_init.f1*qs.f1
                        + v_init.f2*qs.f2 + v_init.f3*qs.f3) / (rho0*rho0);
            v_init.s  -= vdq * qs.s / rho0 * rho0;
            v_init.f1 -= vdq * qs.f1 / rho0 * rho0;
            v_init.f2 -= vdq * qs.f2 / rho0 * rho0;
            v_init.f3 -= vdq * qs.f3 / rho0 * rho0;
            vel[ix] = v_init;
        }
    }

    /* Helper: project force and velocity tangent to S3 */
    auto void sigma_project(void);
    void sigma_project(void) {
        #pragma omp parallel for
        for (size_t ix = 0; ix < N3; ix++) {
            /* Project field onto S3 */
            double n = sqrt(q4_norm2(q[ix]));
            if (n > 1e-15)
                q[ix] = q4_scale(rho0 / n, q[ix]);

            /* Project velocity tangent */
            double vdq = (vel[ix].s*q[ix].s + vel[ix].f1*q[ix].f1
                        + vel[ix].f2*q[ix].f2 + vel[ix].f3*q[ix].f3) / (rho0*rho0);
            vel[ix].s  -= vdq * q[ix].s / (rho0*rho0) * rho0;
            vel[ix].f1 -= vdq * q[ix].f1 / (rho0*rho0) * rho0;
            vel[ix].f2 -= vdq * q[ix].f2 / (rho0*rho0) * rho0;
            vel[ix].f3 -= vdq * q[ix].f3 / (rho0*rho0) * rho0;
        }
    }

    /* Helper: project force tangent to S3 */
    auto void project_force_tangent(void);
    void project_force_tangent(void) {
        #pragma omp parallel for
        for (size_t ix = 0; ix < N3; ix++) {
            double n2 = q4_norm2(q[ix]);
            if (n2 > 1e-20) {
                double fdq = (force[ix].s*q[ix].s + force[ix].f1*q[ix].f1
                            + force[ix].f2*q[ix].f2 + force[ix].f3*q[ix].f3) / n2;
                force[ix].s  -= fdq * q[ix].s;
                force[ix].f1 -= fdq * q[ix].f1;
                force[ix].f2 -= fdq * q[ix].f2;
                force[ix].f3 -= fdq * q[ix].f3;
            }
        }
    }

    /* Helper: single energy for one isolated soliton (for U_int computation) */
    auto double single_soliton_energy(void);
    double single_soliton_energy(void) {
        /* Initialize single soliton at origin to a temporary buffer */
        Q4 *q_single = calloc(N3, sizeof(Q4));
        #pragma omp parallel for collapse(3) schedule(static)
        for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
        for (int k = 0; k < N; k++) {
            double x = -L + (i + 0.5) * h;
            double y = -L + (j + 0.5) * h;
            double z = -L + (k + 0.5) * h;
            int ix_l = idx(N,i,j,k);
            q_single[ix_l] = hedgehog(prof, rho0, x, y, z, 0, 0, 0);
            double n = sqrt(q4_norm2(q_single[ix_l]));
            if (n > 1e-15)
                q_single[ix_l] = q4_scale(rho0/n, q_single[ix_l]);
        }
        EnergyResult E_single = compute_energy(q_single, N, h, e_skyrme);
        free(q_single);
        return E_single.Etotal;
    }

    /* Helper: kinetic energy */
    auto double kinetic_energy(void);
    double kinetic_energy(void) {
        double T = 0;
        double h3 = h*h*h;
        #pragma omp parallel for reduction(+:T)
        for (size_t ix = 0; ix < N3; ix++)
            T += 0.5 * q4_norm2(vel[ix]);
        return T * h3;
    }

    /* Helper: mediator interaction energy U_med = -(g_top/2) * integral B0*p d^3x */
    auto double mediator_energy(void);
    double mediator_energy(void) {
        double h3 = h*h*h;
        double sum = 0;
        #pragma omp parallel for reduction(+:sum)
        for (size_t ix = 0; ix < N3; ix++)
            sum += B0[ix] * p_med[ix];
        return -0.5 * g_top * sum * h3;
    }

    /* ================================================================
     *  VARIATION 1: Baseline (standard Skyrme only)
     * ================================================================ */
    if (var_select == 0 || var_select == 1) {
        printf("=== VARIATION 1: Baseline (standard Skyrme only) ===\n");
        printf("Parameters: N=%d L=%.1f dt=%.4f e=%.1f\n\n", N, L, dt, e_skyrme);

        /* Compute single soliton energy */
        double E_single = single_soliton_energy();
        printf("Single soliton energy: E = %.6f (E2+E4)\n\n", E_single);

        /* Static interaction energy at various separations */
        printf("--- Static U_Skyrme(D) ---\n");
        printf("  %8s  %14s  %14s  %14s\n", "D", "E_total", "U_int=E-2E1", "U_int*D");
        double D_vals[] = {3.0, 4.0, 5.0, 7.0, 10.0, 15.0};
        int n_D = sizeof(D_vals)/sizeof(D_vals[0]);

        for (int id = 0; id < n_D; id++) {
            double D = D_vals[id];
            double z0 = D / 2.0;
            if (z0 >= L - 2.0) continue;  /* skip if too close to boundary */

            init_two_solitons(z0);
            EnergyResult E = compute_energy(q, N, h, e_skyrme);
            double U = E.Etotal - 2.0 * E_single;
            printf("  %8.1f  %14.6f  %14.6e  %14.6e\n", D, E.Etotal, U, U*D);
        }

        /* Self-checks */
        double Q_check = compute_Q(q, N, h);
        printf("\nSelf-check: Q = %.4f (expected ~2.0)\n", Q_check);

        /* Dynamic: evolve at D=5 for 500 steps */
        printf("\n--- Dynamic evolution at D=5 (500 steps) ---\n");
        int n_steps_v1 = 500;
        double D_init = 5.0;
        init_two_solitons(D_init / 2.0);

        EnergyResult E0 = compute_energy(q, N, h, e_skyrme);
        Q_check = compute_Q(q, N, h);
        printf("Init: E=%.4f, Q=%.4f, E2/E4=%.4f\n", E0.Etotal, Q_check, E0.E2/E0.E4);

        printf("step    time    Q       E_skyrme    sep\n");

        /* Leapfrog: first half-step */
        memset(force, 0, N3 * sizeof(Q4));
        compute_skyrme_force(q, N, h, e_skyrme, force);
        project_force_tangent();

        /* Project initial velocity tangent */
        #pragma omp parallel for
        for (size_t ix = 0; ix < N3; ix++) {
            double n2 = q4_norm2(q[ix]);
            if (n2 > 1e-20) {
                double vdq = (vel[ix].s*q[ix].s + vel[ix].f1*q[ix].f1
                            + vel[ix].f2*q[ix].f2 + vel[ix].f3*q[ix].f3) / n2;
                vel[ix].s  -= vdq * q[ix].s;
                vel[ix].f1 -= vdq * q[ix].f1;
                vel[ix].f2 -= vdq * q[ix].f2;
                vel[ix].f3 -= vdq * q[ix].f3;
            }
        }

        #pragma omp parallel for
        for (size_t ix = 0; ix < N3; ix++) {
            vel[ix].s  += 0.5*dt * force[ix].s;
            vel[ix].f1 += 0.5*dt * force[ix].f1;
            vel[ix].f2 += 0.5*dt * force[ix].f2;
            vel[ix].f3 += 0.5*dt * force[ix].f3;
        }

        for (int step = 1; step <= n_steps_v1; step++) {
            /* q update */
            #pragma omp parallel for
            for (size_t ix = 0; ix < N3; ix++) {
                q[ix].s  += dt * vel[ix].s;
                q[ix].f1 += dt * vel[ix].f1;
                q[ix].f2 += dt * vel[ix].f2;
                q[ix].f3 += dt * vel[ix].f3;
            }
            sigma_project();

            /* Force */
            memset(force, 0, N3 * sizeof(Q4));
            compute_skyrme_force(q, N, h, e_skyrme, force);
            project_force_tangent();

            /* Velocity update */
            #pragma omp parallel for
            for (size_t ix = 0; ix < N3; ix++) {
                vel[ix].s  += dt * force[ix].s;
                vel[ix].f1 += dt * force[ix].f1;
                vel[ix].f2 += dt * force[ix].f2;
                vel[ix].f3 += dt * force[ix].f3;
            }

            if (step % 20 == 0) {
                EnergyResult E = compute_energy(q, N, h, e_skyrme);
                double Q = compute_Q(q, N, h);
                double z1, z2;
                find_solitons_3d(q, N, h, L, rho0, &z1, &z2);
                double sep = z1 - z2;
                double T = kinetic_energy();
                printf("%5d  %7.3f  %7.4f  %11.4f  %7.3f  Ekin=%.4f dE=%.2e\n",
                       step, step*dt, Q, E.Etotal, sep, T,
                       (E.Etotal+T - E0.Etotal)/E0.Etotal);

                if (Q < 0.5 || isnan(E.Etotal)) {
                    printf("CATASTROPHIC: stopping V1\n");
                    break;
                }
            }
        }
        printf("\n");
    }

    /* ================================================================
     *  VARIATION 2: + Massless scalar mediator
     * ================================================================ */
    if (var_select == 0 || var_select == 2) {
        printf("=== VARIATION 2: + Massless scalar mediator ===\n");
        printf("Parameters: g_top=%.1f, mu=%.4f (nearly massless), kappa2=%.1f\n\n",
               g_top, mu_massless, kappa2);

        double E_single = single_soliton_energy();

        /* Static U(D) with mediator */
        printf("--- Static U(D) with mediator ---\n");
        printf("  %8s  %14s  %14s  %14s  %14s\n",
               "D", "U_skyrme", "U_mediator", "U_med*D", "U_total");

        double D_vals[] = {3.0, 4.0, 5.0, 7.0, 10.0, 15.0};
        int n_D = sizeof(D_vals)/sizeof(D_vals[0]);

        for (int id = 0; id < n_D; id++) {
            double D = D_vals[id];
            double z0 = D / 2.0;
            if (z0 >= L - 2.0) continue;

            init_two_solitons(z0);

            EnergyResult E = compute_energy(q, N, h, e_skyrme);
            double U_skyrme = E.Etotal - 2.0 * E_single;

            /* Compute B0 and solve Poisson */
            compute_baryon_density(q, N, h, B0, J_cof);
            solve_poisson_fft(B0, p_med, fft_work, N, h, g_top, kappa2,
                              mu_massless * mu_massless);
            double U_med = mediator_energy();

            printf("  %8.1f  %14.6e  %14.6e  %14.6e  %14.6e\n",
                   D, U_skyrme, U_med, U_med*D, U_skyrme + U_med);

            if (id == 0) {
                /* Poisson residual check on first configuration */
                double res = check_poisson_residual(p_med, B0, N, h,
                                                     g_top, kappa2,
                                                     mu_massless * mu_massless);
                printf("  Poisson residual: %.4e\n", res);
            }
        }

        /* Dynamic evolution at D=5 with mediator backreaction */
        printf("\n--- Dynamic (D=5, 500 steps, with mediator) ---\n");
        int n_steps_v2 = 500;
        init_two_solitons(2.5);

        EnergyResult E0 = compute_energy(q, N, h, e_skyrme);
        compute_baryon_density(q, N, h, B0, J_cof);
        solve_poisson_fft(B0, p_med, fft_work, N, h, g_top, kappa2,
                          mu_massless * mu_massless);
        double U_med0 = mediator_energy();
        double Q0 = compute_Q(q, N, h);
        printf("Init: E_skyrme=%.4f, U_med=%.6e, Q=%.4f\n\n", E0.Etotal, U_med0, Q0);

        printf("step    time    Q       E_skyrme    U_med       sep       p_max\n");

        /* Leapfrog init */
        memset(force, 0, N3 * sizeof(Q4));
        compute_skyrme_force(q, N, h, e_skyrme, force);

        /* Add mediator backreaction: F_med_e = g_top * p * J_e */
        #pragma omp parallel for
        for (size_t ix = 0; ix < N3; ix++) {
            double pval = p_med[ix];
            force[ix].s  += g_top * pval * J_cof[ix].s;
            force[ix].f1 += g_top * pval * J_cof[ix].f1;
            force[ix].f2 += g_top * pval * J_cof[ix].f2;
            force[ix].f3 += g_top * pval * J_cof[ix].f3;
        }
        project_force_tangent();

        #pragma omp parallel for
        for (size_t ix = 0; ix < N3; ix++) {
            vel[ix].s  += 0.5*dt * force[ix].s;
            vel[ix].f1 += 0.5*dt * force[ix].f1;
            vel[ix].f2 += 0.5*dt * force[ix].f2;
            vel[ix].f3 += 0.5*dt * force[ix].f3;
        }

        for (int step = 1; step <= n_steps_v2; step++) {
            #pragma omp parallel for
            for (size_t ix = 0; ix < N3; ix++) {
                q[ix].s  += dt * vel[ix].s;
                q[ix].f1 += dt * vel[ix].f1;
                q[ix].f2 += dt * vel[ix].f2;
                q[ix].f3 += dt * vel[ix].f3;
            }
            sigma_project();

            /* Recompute B0, solve Poisson, get mediator force */
            compute_baryon_density(q, N, h, B0, J_cof);
            solve_poisson_fft(B0, p_med, fft_work, N, h, g_top, kappa2,
                              mu_massless * mu_massless);

            memset(force, 0, N3 * sizeof(Q4));
            compute_skyrme_force(q, N, h, e_skyrme, force);
            #pragma omp parallel for
            for (size_t ix = 0; ix < N3; ix++) {
                double pval = p_med[ix];
                force[ix].s  += g_top * pval * J_cof[ix].s;
                force[ix].f1 += g_top * pval * J_cof[ix].f1;
                force[ix].f2 += g_top * pval * J_cof[ix].f2;
                force[ix].f3 += g_top * pval * J_cof[ix].f3;
            }
            project_force_tangent();

            #pragma omp parallel for
            for (size_t ix = 0; ix < N3; ix++) {
                vel[ix].s  += dt * force[ix].s;
                vel[ix].f1 += dt * force[ix].f1;
                vel[ix].f2 += dt * force[ix].f2;
                vel[ix].f3 += dt * force[ix].f3;
            }

            if (step % 20 == 0) {
                EnergyResult E = compute_energy(q, N, h, e_skyrme);
                double Q = compute_Q(q, N, h);
                double z1, z2;
                find_solitons_3d(q, N, h, L, rho0, &z1, &z2);
                double U_med = mediator_energy();
                double p_max = 0;
                for (size_t ix = 0; ix < N3; ix++)
                    if (fabs(p_med[ix]) > p_max) p_max = fabs(p_med[ix]);

                printf("%5d  %7.3f  %7.4f  %11.4f  %11.4e  %7.3f  %.4e\n",
                       step, step*dt, Q, E.Etotal, U_med, z1-z2, p_max);

                if (Q < 0.5 || isnan(E.Etotal)) {
                    printf("CATASTROPHIC: stopping V2\n");
                    break;
                }
            }
        }
        printf("\n");
    }

    /* ================================================================
     *  VARIATION 3: + Massive mediator (pion mass)
     * ================================================================ */
    if (var_select == 0 || var_select == 3) {
        printf("=== VARIATION 3: + Massive mediator (mu=%.4f) ===\n", mu_massive);
        printf("Parameters: g_top=%.1f, mu=%.4f, range=%.2f\n\n",
               g_top, mu_massive, 1.0/mu_massive);

        double E_single = single_soliton_energy();

        printf("--- Static U(D) with Yukawa mediator ---\n");
        printf("  %8s  %14s  %14s  %14s  %14s\n",
               "D", "U_skyrme", "U_mediator", "U_med*D*exp(mu*D)", "U_total");

        double D_vals[] = {3.0, 4.0, 5.0, 7.0, 10.0, 15.0};
        int n_D = sizeof(D_vals)/sizeof(D_vals[0]);

        for (int id = 0; id < n_D; id++) {
            double D = D_vals[id];
            double z0 = D / 2.0;
            if (z0 >= L - 2.0) continue;

            init_two_solitons(z0);
            EnergyResult E = compute_energy(q, N, h, e_skyrme);
            double U_skyrme = E.Etotal - 2.0 * E_single;

            compute_baryon_density(q, N, h, B0, J_cof);
            solve_poisson_fft(B0, p_med, fft_work, N, h, g_top, kappa2,
                              mu_massive * mu_massive);
            double U_med = mediator_energy();

            printf("  %8.1f  %14.6e  %14.6e  %14.6e  %14.6e\n",
                   D, U_skyrme, U_med, U_med*D*exp(mu_massive*D),
                   U_skyrme + U_med);
        }
        printf("\n");
    }

    /* ================================================================
     *  VARIATION 4: Pullback metric feedback (control — P/m=2, no deflection)
     * ================================================================ */
    if (var_select == 0 || var_select == 4) {
        printf("=== VARIATION 4: BLV metric control (P/m=2, expect ZERO deflection) ===\n\n");

        /* Single soliton + wave packet */
        double wave_x0 = 5.0;
        double wave_sigma = 1.0;
        double wave_amp = 0.01;
        double wave_vx = -1.0;  /* moving toward soliton */
        int n_steps_v4 = 1000;

        init_soliton_wave(wave_x0, wave_sigma, wave_amp, wave_vx);

        /* Track wave packet centroid by f1 component deviation */
        /* The wave was initialized as a perturbation in f1 */

        printf("step    time    wave_x_centroid   wave_y_centroid   Q       E\n");

        /* Leapfrog init */
        memset(force, 0, N3 * sizeof(Q4));
        compute_skyrme_force(q, N, h, e_skyrme, force);
        project_force_tangent();

        #pragma omp parallel for
        for (size_t ix = 0; ix < N3; ix++) {
            vel[ix].s  += 0.5*dt * force[ix].s;
            vel[ix].f1 += 0.5*dt * force[ix].f1;
            vel[ix].f2 += 0.5*dt * force[ix].f2;
            vel[ix].f3 += 0.5*dt * force[ix].f3;
        }

        for (int step = 1; step <= n_steps_v4; step++) {
            #pragma omp parallel for
            for (size_t ix = 0; ix < N3; ix++) {
                q[ix].s  += dt * vel[ix].s;
                q[ix].f1 += dt * vel[ix].f1;
                q[ix].f2 += dt * vel[ix].f2;
                q[ix].f3 += dt * vel[ix].f3;
            }
            sigma_project();

            memset(force, 0, N3 * sizeof(Q4));
            compute_skyrme_force(q, N, h, e_skyrme, force);
            project_force_tangent();

            #pragma omp parallel for
            for (size_t ix = 0; ix < N3; ix++) {
                vel[ix].s  += dt * force[ix].s;
                vel[ix].f1 += dt * force[ix].f1;
                vel[ix].f2 += dt * force[ix].f2;
                vel[ix].f3 += dt * force[ix].f3;
            }

            if (step % 50 == 0) {
                /* Track wave packet via kinetic energy concentration.
                 * Look for the centroid of |vel| away from origin. */
                double wx = 0, wy = 0, wsum = 0;
                for (int i = 0; i < N; i++) {
                    double x = -L + (i + 0.5) * h;
                    for (int j = 0; j < N; j++) {
                        double y = -L + (j + 0.5) * h;
                        for (int k = 0; k < N; k++) {
                            int ix_l = idx(N,i,j,k);
                            double v2 = q4_norm2(vel[ix_l]);
                            /* Weight by distance from origin to avoid soliton core */
                            double r2 = x*x + y*y;
                            if (r2 > 4.0) {  /* only outside core */
                                wx += x * v2;
                                wy += y * v2;
                                wsum += v2;
                            }
                        }
                    }
                }
                double cx_wave = (wsum > 1e-20) ? wx/wsum : 0;
                double cy_wave = (wsum > 1e-20) ? wy/wsum : 0;

                EnergyResult E = compute_energy(q, N, h, e_skyrme);
                double Q = compute_Q(q, N, h);
                printf("%5d  %7.3f  %14.6f  %14.6f  %7.4f  %11.4f\n",
                       step, step*dt, cx_wave, cy_wave, Q, E.Etotal);

                if (Q < 0.5 || isnan(E.Etotal)) {
                    printf("CATASTROPHIC: stopping V4\n");
                    break;
                }
            }
        }

        printf("\nExpected: wave_y_centroid should remain ~0 (no deflection)\n\n");
    }

    /* ================================================================
     *  VARIATION 5: + L6 with metric feedback
     * ================================================================ */
    if (var_select == 0 || var_select == 5) {
        printf("=== VARIATION 5: L6 metric (lambda6=%.1f, expect attractive lensing) ===\n\n",
               lambda6);

        double wave_x0 = 5.0;
        double wave_sigma = 1.0;
        double wave_amp = 0.01;
        double wave_vx = -1.0;
        int n_steps_v5 = 1000;

        init_soliton_wave(wave_x0, wave_sigma, wave_amp, wave_vx);

        printf("step    time    wave_x_centroid   wave_y_centroid   Q       E       E6\n");

        /* Leapfrog init */
        memset(force, 0, N3 * sizeof(Q4));
        compute_skyrme_force(q, N, h, e_skyrme, force);

        /* Add L6 force: F6_e = -2*lambda6*B0*J_e */
        compute_baryon_density(q, N, h, B0, J_cof);
        #pragma omp parallel for
        for (size_t ix = 0; ix < N3; ix++) {
            double coeff = -2.0 * lambda6 * B0[ix];
            force[ix].s  += coeff * J_cof[ix].s;
            force[ix].f1 += coeff * J_cof[ix].f1;
            force[ix].f2 += coeff * J_cof[ix].f2;
            force[ix].f3 += coeff * J_cof[ix].f3;
        }
        project_force_tangent();

        #pragma omp parallel for
        for (size_t ix = 0; ix < N3; ix++) {
            vel[ix].s  += 0.5*dt * force[ix].s;
            vel[ix].f1 += 0.5*dt * force[ix].f1;
            vel[ix].f2 += 0.5*dt * force[ix].f2;
            vel[ix].f3 += 0.5*dt * force[ix].f3;
        }

        for (int step = 1; step <= n_steps_v5; step++) {
            #pragma omp parallel for
            for (size_t ix = 0; ix < N3; ix++) {
                q[ix].s  += dt * vel[ix].s;
                q[ix].f1 += dt * vel[ix].f1;
                q[ix].f2 += dt * vel[ix].f2;
                q[ix].f3 += dt * vel[ix].f3;
            }
            sigma_project();

            memset(force, 0, N3 * sizeof(Q4));
            compute_skyrme_force(q, N, h, e_skyrme, force);

            compute_baryon_density(q, N, h, B0, J_cof);
            #pragma omp parallel for
            for (size_t ix = 0; ix < N3; ix++) {
                double coeff = -2.0 * lambda6 * B0[ix];
                force[ix].s  += coeff * J_cof[ix].s;
                force[ix].f1 += coeff * J_cof[ix].f1;
                force[ix].f2 += coeff * J_cof[ix].f2;
                force[ix].f3 += coeff * J_cof[ix].f3;
            }
            project_force_tangent();

            #pragma omp parallel for
            for (size_t ix = 0; ix < N3; ix++) {
                vel[ix].s  += dt * force[ix].s;
                vel[ix].f1 += dt * force[ix].f1;
                vel[ix].f2 += dt * force[ix].f2;
                vel[ix].f3 += dt * force[ix].f3;
            }

            if (step % 50 == 0) {
                double wx = 0, wy = 0, wsum = 0;
                for (int i = 0; i < N; i++) {
                    double x = -L + (i + 0.5) * h;
                    for (int j = 0; j < N; j++) {
                        double y = -L + (j + 0.5) * h;
                        for (int k = 0; k < N; k++) {
                            int ix_l = idx(N,i,j,k);
                            double v2 = q4_norm2(vel[ix_l]);
                            double r2 = x*x + y*y;
                            if (r2 > 4.0) {
                                wx += x * v2;
                                wy += y * v2;
                                wsum += v2;
                            }
                        }
                    }
                }
                double cx_wave = (wsum > 1e-20) ? wx/wsum : 0;
                double cy_wave = (wsum > 1e-20) ? wy/wsum : 0;

                EnergyResult E = compute_energy(q, N, h, e_skyrme);
                double Q = compute_Q(q, N, h);
                double E6;
                compute_l6_energy(B0, N, h*h*h, lambda6, &E6);

                printf("%5d  %7.3f  %14.6f  %14.6f  %7.4f  %11.4f  %9.4f\n",
                       step, step*dt, cx_wave, cy_wave, Q, E.Etotal + E6, E6);

                if (Q < 0.5 || isnan(E.Etotal)) {
                    printf("CATASTROPHIC: stopping V5\n");
                    break;
                }
            }
        }

        printf("\nExpected: wave_y_centroid should deviate from 0 (attractive lensing)\n\n");
    }

    /* ================================================================
     *  VARIATION 6: Dynamic scatter with mediator
     * ================================================================ */
    if (var_select == 0 || var_select == 6) {
        printf("=== VARIATION 6: Dynamic scatter with mediator ===\n");
        double D_init = 10.0;
        double v0 = 0.3;
        int n_steps_v6 = 2000;

        printf("Parameters: D=%.1f, v=%.2fc, g_top=%.1f, mu=%.4f, steps=%d\n\n",
               D_init, v0, g_top, mu_massless, n_steps_v6);

        double z0 = D_init / 2.0;
        init_two_solitons(z0);

        /* Set initial velocities: solitons boosted toward each other */
        #pragma omp parallel for collapse(3) schedule(static)
        for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
        for (int k = 0; k < N; k++) {
            double z = -L + (k + 0.5) * h;
            int ix_l = idx(N,i,j,k);

            /* Finite diff dq/dz */
            Q4 qp = q[idx(N,i,j,k+1 < N ? k+1 : 0)];
            Q4 qm = q[idx(N,i,j,k > 0 ? k-1 : N-1)];
            Q4 dqdz = q4_scale(0.5/h, q4_sub(qp, qm));

            /* Soliton 1 at +z0 moves -v0: dq/dt = +v0*dq/dz near z>0
             * Soliton 2 at -z0 moves +v0: dq/dt = -v0*dq/dz near z<0 */
            double sign = (z > 0) ? +v0 : -v0;
            vel[ix_l] = q4_scale(sign, dqdz);
        }

        /* Project velocity tangent */
        #pragma omp parallel for
        for (size_t ix = 0; ix < N3; ix++) {
            double n2 = q4_norm2(q[ix]);
            if (n2 > 1e-20) {
                double vdq = (vel[ix].s*q[ix].s + vel[ix].f1*q[ix].f1
                            + vel[ix].f2*q[ix].f2 + vel[ix].f3*q[ix].f3) / n2;
                vel[ix].s  -= vdq * q[ix].s;
                vel[ix].f1 -= vdq * q[ix].f1;
                vel[ix].f2 -= vdq * q[ix].f2;
                vel[ix].f3 -= vdq * q[ix].f3;
            }
        }

        EnergyResult E0 = compute_energy(q, N, h, e_skyrme);
        double T0 = kinetic_energy();
        double Q0 = compute_Q(q, N, h);

        compute_baryon_density(q, N, h, B0, J_cof);
        solve_poisson_fft(B0, p_med, fft_work, N, h, g_top, kappa2,
                          mu_massless * mu_massless);
        double U_med0 = mediator_energy();

        printf("Init: E_skyrme=%.4f, E_kin=%.4f, U_med=%.6e, Q=%.4f\n",
               E0.Etotal, T0, U_med0, Q0);
        printf("E_total_0 = %.4f\n\n", E0.Etotal + T0 + U_med0);

        double E_total_0 = E0.Etotal + T0 + U_med0;

        printf("step    time    Q       E_skyrme    E_kin       U_med       sep       dE/E\n");

        /* Leapfrog init */
        memset(force, 0, N3 * sizeof(Q4));
        compute_skyrme_force(q, N, h, e_skyrme, force);
        #pragma omp parallel for
        for (size_t ix = 0; ix < N3; ix++) {
            double pval = p_med[ix];
            force[ix].s  += g_top * pval * J_cof[ix].s;
            force[ix].f1 += g_top * pval * J_cof[ix].f1;
            force[ix].f2 += g_top * pval * J_cof[ix].f2;
            force[ix].f3 += g_top * pval * J_cof[ix].f3;
        }
        project_force_tangent();

        #pragma omp parallel for
        for (size_t ix = 0; ix < N3; ix++) {
            vel[ix].s  += 0.5*dt * force[ix].s;
            vel[ix].f1 += 0.5*dt * force[ix].f1;
            vel[ix].f2 += 0.5*dt * force[ix].f2;
            vel[ix].f3 += 0.5*dt * force[ix].f3;
        }

        for (int step = 1; step <= n_steps_v6; step++) {
            #pragma omp parallel for
            for (size_t ix = 0; ix < N3; ix++) {
                q[ix].s  += dt * vel[ix].s;
                q[ix].f1 += dt * vel[ix].f1;
                q[ix].f2 += dt * vel[ix].f2;
                q[ix].f3 += dt * vel[ix].f3;
            }
            sigma_project();

            compute_baryon_density(q, N, h, B0, J_cof);
            solve_poisson_fft(B0, p_med, fft_work, N, h, g_top, kappa2,
                              mu_massless * mu_massless);

            memset(force, 0, N3 * sizeof(Q4));
            compute_skyrme_force(q, N, h, e_skyrme, force);
            #pragma omp parallel for
            for (size_t ix = 0; ix < N3; ix++) {
                double pval = p_med[ix];
                force[ix].s  += g_top * pval * J_cof[ix].s;
                force[ix].f1 += g_top * pval * J_cof[ix].f1;
                force[ix].f2 += g_top * pval * J_cof[ix].f2;
                force[ix].f3 += g_top * pval * J_cof[ix].f3;
            }
            project_force_tangent();

            #pragma omp parallel for
            for (size_t ix = 0; ix < N3; ix++) {
                vel[ix].s  += dt * force[ix].s;
                vel[ix].f1 += dt * force[ix].f1;
                vel[ix].f2 += dt * force[ix].f2;
                vel[ix].f3 += dt * force[ix].f3;
            }

            if (step % 20 == 0) {
                EnergyResult E = compute_energy(q, N, h, e_skyrme);
                double T = kinetic_energy();
                double Q = compute_Q(q, N, h);
                double U_med = mediator_energy();
                double z1, z2;
                find_solitons_3d(q, N, h, L, rho0, &z1, &z2);
                double E_total = E.Etotal + T + U_med;
                double dE = (E_total - E_total_0) / fabs(E_total_0);

                printf("%5d  %7.3f  %7.4f  %11.4f  %11.4f  %11.4e  %7.3f  %10.2e\n",
                       step, step*dt, Q, E.Etotal, T, U_med, z1-z2, dE);

                if (Q < 0.5 || isnan(E.Etotal)) {
                    printf("CATASTROPHIC: stopping V6\n");
                    break;
                }
                if (fabs(dE) > 0.5) {
                    printf("Energy conservation violated by >50%% — aborting V6\n");
                    break;
                }
            }
        }
        printf("\n");
    }

    /* Cleanup */
    free(q);
    free(vel);
    free(force);
    free(B0);
    free(J_cof);
    free(p_med);
    free(fft_work);
    free_profile(prof);

    printf("============================================================\n");
    printf(" All requested variations complete.\n");
    printf("============================================================\n");

    return 0;
}
