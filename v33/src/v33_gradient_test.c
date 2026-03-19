/*  v33_gradient_test.c — Braid in an imposed density gradient
 *
 *  x-direction: pinned boundaries (high ρ left, low ρ right)
 *  y,z-directions: free-floating outflow (linear extrapolation)
 *
 *  Test: does the braid drift along ∇ρ?
 *  If yes: the standard equation captures gravity from density gradients.
 *
 *  Build: gcc -O3 -march=native -fopenmp -o v33_gtest src/v33_gradient_test.c -lm
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>
#include <sys/stat.h>

#define NFIELDS 3
#define PI 3.14159265358979323846

static double MU = -41.345, KAPPA = 50.0, MASS2 = 2.25;
static double A_HIGH = 0.15;   /* background amplitude at x=-L (high ρ) */
static double A_LOW  = 0.05;   /* background amplitude at x=+L (low ρ) */

/* ================================================================
   Grid: ONE allocation
   ================================================================ */

typedef struct {
    double *mem;
    double *phi[NFIELDS], *vel[NFIELDS], *acc[NFIELDS];
    /* Pinned boundary values (saved at init) */
    double *pin_phi[NFIELDS], *pin_vel[NFIELDS];
    int N; long N3;
    double L, dx, dt;
} Grid;

static Grid *grid_alloc(int N, double L) {
    Grid *g = calloc(1, sizeof(Grid));
    g->N = N; g->N3 = (long)N*N*N;
    g->L = L; g->dx = 2.0*L/(N-1);
    g->dt = 0.12 * g->dx;

    /* Main fields: ONE allocation */
    long total = 9 * g->N3;
    printf("Allocating %.2f GB (%ld doubles, N=%d)\n", total*8.0/1e9, total, N);
    g->mem = malloc(total * sizeof(double));
    if (!g->mem) { fprintf(stderr, "FATAL: malloc failed\n"); exit(1); }
    memset(g->mem, 0, total * sizeof(double));
    for (int a = 0; a < NFIELDS; a++) {
        g->phi[a] = g->mem + (0+a) * g->N3;
        g->vel[a] = g->mem + (3+a) * g->N3;
        g->acc[a] = g->mem + (6+a) * g->N3;
    }

    /* Pinned boundary storage (separate, small — only boundary slabs) */
    for (int a = 0; a < NFIELDS; a++) {
        g->pin_phi[a] = calloc(g->N3, sizeof(double));
        g->pin_vel[a] = calloc(g->N3, sizeof(double));
    }

    return g;
}

static void grid_free(Grid *g) {
    free(g->mem);
    for (int a = 0; a < NFIELDS; a++) {
        free(g->pin_phi[a]); free(g->pin_vel[a]);
    }
    free(g);
}

/* ================================================================
   Initialization: braid + gradient background
   ================================================================ */

static void init_gradient(Grid *g) {
    const int N = g->N, NN = N*N;
    const double dx = g->dx, L = g->L;
    const double kw = PI/L, omega = sqrt(kw*kw + MASS2);

    /* Braid params */
    const double A[3] = {0.8, 0.8, 0.8};
    const double delta[3] = {0, 3.0005, 4.4325};
    const double R_tube = 3.0, ellip = 0.3325;
    const double sx = 1+ellip, sy = 1-ellip;
    const double inv2R2 = 1.0/(2*R_tube*R_tube);

    for (int i = 0; i < N; i++) {
        double x = -L + i*dx;

        /* Background amplitude varies linearly in x */
        double frac = (x + L) / (2*L);  /* 0 at x=-L, 1 at x=+L */
        double A_bg = A_HIGH * (1.0 - frac) + A_LOW * frac;
        double k_bg = PI/L;
        double omega_bg = sqrt(k_bg*k_bg + MASS2);

        for (int j = 0; j < N; j++) {
            double y = -L + j*dx;
            for (int kk = 0; kk < N; kk++) {
                double z = -L + kk*dx;
                long idx = (long)i*NN + j*N + kk;

                /* Braid (centered at origin) */
                double xr = x, yr = y;
                double r2e = xr*xr/(sx*sx) + yr*yr/(sy*sy);
                double env = exp(-r2e * inv2R2);

                for (int a = 0; a < NFIELDS; a++) {
                    double ph = kw*z + delta[a];
                    double ph_bg = k_bg*z + 2*PI*a/3.0;
                    g->phi[a][idx] = A[a]*env*cos(ph) + A_bg*cos(ph_bg);
                    g->vel[a][idx] = omega*A[a]*env*sin(ph) + omega_bg*A_bg*sin(ph_bg);
                }
            }
        }
    }

    /* Save initial state as pinned boundary values */
    for (int a = 0; a < NFIELDS; a++) {
        memcpy(g->pin_phi[a], g->phi[a], g->N3 * sizeof(double));
        memcpy(g->pin_vel[a], g->vel[a], g->N3 * sizeof(double));
    }
}

/* ================================================================
   Boundary conditions
   ================================================================ */

static void apply_bc(Grid *g) {
    const int N = g->N, NN = N*N;
    const int margin = 3;

    /* x-direction: PIN to saved boundary values */
    #pragma omp parallel for schedule(static)
    for (long idx = 0; idx < g->N3; idx++) {
        int i = (int)(idx / NN);
        if (i < margin || i >= N - margin) {
            for (int a = 0; a < NFIELDS; a++) {
                g->phi[a][idx] = g->pin_phi[a][idx];
                g->vel[a][idx] = g->pin_vel[a][idx];
                g->acc[a][idx] = 0;
            }
        }
    }

    /* y-direction: FREE-FLOATING (linear extrapolation from interior) */
    #pragma omp parallel for schedule(static)
    for (int i = margin; i < N - margin; i++) {
        for (int kk = 0; kk < N; kk++) {
            for (int a = 0; a < NFIELDS; a++) {
                /* j=0,1: extrapolate from j=2,3 */
                for (int b = 0; b < margin; b++) {
                    int j_in1 = margin;
                    int j_in2 = margin + 1;
                    long idx_b  = (long)i*NN + b*N + kk;
                    long idx_i1 = (long)i*NN + j_in1*N + kk;
                    long idx_i2 = (long)i*NN + j_in2*N + kk;
                    g->phi[a][idx_b] = 2*g->phi[a][idx_i1] - g->phi[a][idx_i2];
                    g->vel[a][idx_b] = 2*g->vel[a][idx_i1] - g->vel[a][idx_i2];
                    g->acc[a][idx_b] = 0;
                }
                /* j=N-1,N-2: extrapolate from j=N-3, N-4 */
                for (int b = 0; b < margin; b++) {
                    int j_out = N - 1 - b;
                    int j_in1 = N - 1 - margin;
                    int j_in2 = N - 2 - margin;
                    long idx_b  = (long)i*NN + j_out*N + kk;
                    long idx_i1 = (long)i*NN + j_in1*N + kk;
                    long idx_i2 = (long)i*NN + j_in2*N + kk;
                    g->phi[a][idx_b] = 2*g->phi[a][idx_i1] - g->phi[a][idx_i2];
                    g->vel[a][idx_b] = 2*g->vel[a][idx_i1] - g->vel[a][idx_i2];
                    g->acc[a][idx_b] = 0;
                }
            }
        }
    }

    /* z-direction: FREE-FLOATING (linear extrapolation) */
    #pragma omp parallel for schedule(static)
    for (int i = margin; i < N - margin; i++) {
        for (int j = margin; j < N - margin; j++) {
            for (int a = 0; a < NFIELDS; a++) {
                for (int b = 0; b < margin; b++) {
                    /* k=0,1,2 */
                    long idx_b  = (long)i*NN + j*N + b;
                    long idx_i1 = (long)i*NN + j*N + margin;
                    long idx_i2 = (long)i*NN + j*N + margin + 1;
                    g->phi[a][idx_b] = 2*g->phi[a][idx_i1] - g->phi[a][idx_i2];
                    g->vel[a][idx_b] = 2*g->vel[a][idx_i1] - g->vel[a][idx_i2];
                    g->acc[a][idx_b] = 0;
                    /* k=N-1,N-2,N-3 */
                    int k_out = N - 1 - b;
                    long idx_ob = (long)i*NN + j*N + k_out;
                    long idx_o1 = (long)i*NN + j*N + (N-1-margin);
                    long idx_o2 = (long)i*NN + j*N + (N-2-margin);
                    g->phi[a][idx_ob] = 2*g->phi[a][idx_o1] - g->phi[a][idx_o2];
                    g->vel[a][idx_ob] = 2*g->vel[a][idx_o1] - g->vel[a][idx_o2];
                    g->acc[a][idx_ob] = 0;
                }
            }
        }
    }
}

/* ================================================================
   Forces: standard equation, NON-periodic (clamped neighbors at edges)
   ================================================================ */

static void compute_forces(Grid *g) {
    const int N = g->N, NN = N*N;
    const long N3 = g->N3;
    const double idx2 = 1.0/(g->dx*g->dx);
    const int margin = 3;

    #pragma omp parallel for schedule(static)
    for (long idx = 0; idx < N3; idx++) {
        int i = (int)(idx/NN), j = (int)((idx/N)%N), k = (int)(idx%N);

        /* Skip boundary cells (handled by apply_bc) */
        if (i < margin || i >= N-margin ||
            j < margin || j >= N-margin ||
            k < margin || k >= N-margin) {
            continue;  /* acc already set to 0 by apply_bc */
        }

        /* Interior: standard 6-point Laplacian (NOT periodic) */
        long n_ip = (long)(i+1)*NN + j*N + k;
        long n_im = (long)(i-1)*NN + j*N + k;
        long n_jp = (long)i*NN + (j+1)*N + k;
        long n_jm = (long)i*NN + (j-1)*N + k;
        long n_kp = (long)i*NN + j*N + (k+1);
        long n_km = (long)i*NN + j*N + (k-1);

        double p0 = g->phi[0][idx], p1 = g->phi[1][idx], p2 = g->phi[2][idx];
        double P = p0*p1*p2;
        double den = 1.0 + KAPPA*P*P;
        double mPd2 = MU*P/(den*den);

        for (int a = 0; a < NFIELDS; a++) {
            double lap = (g->phi[a][n_ip] + g->phi[a][n_im]
                        + g->phi[a][n_jp] + g->phi[a][n_jm]
                        + g->phi[a][n_kp] + g->phi[a][n_km]
                        - 6.0*g->phi[a][idx]) * idx2;
            double dPda = (a==0)?p1*p2:(a==1)?p0*p2:p0*p1;
            g->acc[a][idx] = lap - MASS2*g->phi[a][idx] - mPd2*dPda;
        }
    }
}

/* ================================================================
   Verlet + BC
   ================================================================ */

static void verlet_step(Grid *g) {
    const long N3 = g->N3;
    const double hdt = 0.5*g->dt, dt = g->dt;

    for (int a=0;a<NFIELDS;a++) {
        double *v=g->vel[a], *ac=g->acc[a];
        for (long idx=0;idx<N3;idx++) v[idx] += hdt*ac[idx];
    }
    for (int a=0;a<NFIELDS;a++) {
        double *p=g->phi[a], *v=g->vel[a];
        for (long idx=0;idx<N3;idx++) p[idx] += dt*v[idx];
    }
    apply_bc(g);
    compute_forces(g);
    for (int a=0;a<NFIELDS;a++) {
        double *v=g->vel[a], *ac=g->acc[a];
        for (long idx=0;idx<N3;idx++) v[idx] += hdt*ac[idx];
    }
}

/* ================================================================
   Diagnostics: track braid center of mass in x
   ================================================================ */

static double measure_braid_x(Grid *g) {
    const long N3 = g->N3;
    const int NN = g->N*g->N, N = g->N;
    const double dx = g->dx, L = g->L;

    /* Find high-energy centroid */
    double avg = 0;
    for (long idx = 0; idx < N3; idx++) {
        double p2 = 0;
        for (int a=0;a<NFIELDS;a++) p2 += g->phi[a][idx]*g->phi[a][idx];
        avg += p2;
    }
    avg /= N3;
    double thresh = 5.0 * avg;

    double wx = 0, wt = 0;
    for (long idx = 0; idx < N3; idx++) {
        double p2 = 0;
        for (int a=0;a<NFIELDS;a++) p2 += g->phi[a][idx]*g->phi[a][idx];
        if (p2 < thresh) continue;
        int i = (int)(idx / NN);
        double x = -L + i*dx;
        wx += x * p2;
        wt += p2;
    }
    return (wt > 0) ? wx/wt : 0;
}

static void save_field(Grid *g, double t, const char *dir) {
    char fn[512];
    snprintf(fn, sizeof(fn), "%s/field_t%04d.bin", dir, (int)(t+0.5));
    FILE *fp = fopen(fn, "wb"); if(!fp) return;
    int n=g->N; double l=g->L;
    fwrite(&n,sizeof(int),1,fp); fwrite(&l,sizeof(double),1,fp); fwrite(&t,sizeof(double),1,fp);
    for (int a=0;a<NFIELDS;a++) fwrite(g->phi[a],sizeof(double),g->N3,fp);
    fclose(fp);
}

/* ================================================================
   Main
   ================================================================ */

int main(int argc, char **argv) {
    int N = 512; double L = 20.0, T = 500.0;
    char outdir[256] = "data/gtest";

    for (int i=1;i<argc;i++) {
        if (!strcmp(argv[i],"-N")) N=atoi(argv[++i]);
        else if (!strcmp(argv[i],"-L")) L=atof(argv[++i]);
        else if (!strcmp(argv[i],"-T")) T=atof(argv[++i]);
        else if (!strcmp(argv[i],"-ah")) A_HIGH=atof(argv[++i]);
        else if (!strcmp(argv[i],"-al")) A_LOW=atof(argv[++i]);
        else if (!strcmp(argv[i],"-o")) strncpy(outdir,argv[++i],255);
    }

    omp_set_num_threads(16);
    printf("=== V33 Gradient Test: Braid in ρ Gradient ===\n");
    printf("x: pinned (A_high=%.3f left, A_low=%.3f right)\n", A_HIGH, A_LOW);
    printf("y,z: free-floating (linear extrapolation)\n");
    printf("N=%d L=%.0f T=%.0f\n", N, L, T);
    printf("ρ_high/ρ_low ≈ %.1f (amplitude ratio² = %.1f)\n",
           (A_HIGH*A_HIGH)/(A_LOW*A_LOW), (A_HIGH/A_LOW)*(A_HIGH/A_LOW));

    mkdir("data", 0755); mkdir(outdir, 0755);

    Grid *g = grid_alloc(N, L);
    printf("dx=%.4f dt=%.5f\n\n", g->dx, g->dt);

    init_gradient(g);
    apply_bc(g);
    compute_forces(g);

    char tspath[512];
    snprintf(tspath, sizeof(tspath), "%s/timeseries.tsv", outdir);
    FILE *fp = fopen(tspath, "w");
    fprintf(fp, "t\tx_braid\tE_total\n");

    int n_steps = (int)(T/g->dt);
    int diag_every = n_steps/100; if (diag_every<1) diag_every=1;
    int snap_every = n_steps/5; if (snap_every<1) snap_every=1;

    printf("Steps=%d\n\n", n_steps);
    double wall0 = omp_get_wtime();
    double x0 = measure_braid_x(g);
    printf("Initial braid x = %.4f\n\n", x0);
    save_field(g, 0, outdir);

    for (int step = 0; step <= n_steps; step++) {
        if (step > 0) verlet_step(g);
        double t = step * g->dt;

        if (step % diag_every == 0) {
            double xb = measure_braid_x(g);
            /* Quick energy estimate */
            double E = 0;
            for (long idx = 0; idx < g->N3; idx++)
                for (int a=0;a<NFIELDS;a++)
                    E += 0.5*g->vel[a][idx]*g->vel[a][idx]*g->dx*g->dx*g->dx;

            fprintf(fp, "%.2f\t%.6f\t%.4e\n", t, xb, E);
            fflush(fp);

            if (step % (diag_every*10) == 0) {
                double wall = omp_get_wtime()-wall0;
                printf("t=%7.1f x_braid=%+.4f (drift=%+.4f) E_kin=%.2e [%.0f%% %.0fs]\n",
                       t, xb, xb-x0, E, 100.0*step/n_steps, wall);
            }
        }
        if (step > 0 && step % snap_every == 0)
            save_field(g, t, outdir);
    }

    fclose(fp);
    double xf = measure_braid_x(g);
    save_field(g, T, outdir);

    printf("\n=== RESULT ===\n");
    printf("x_braid: %.4f → %.4f (drift = %+.4f)\n", x0, xf, xf-x0);
    printf("If drift < 0: braid moved toward HIGH ρ (left) = GRAVITY\n");
    printf("If drift > 0: braid moved toward LOW ρ (right) = ANTI-GRAVITY\n");
    printf("If drift ≈ 0: no response to gradient\n");
    printf("\n=== Complete: %.0fs ===\n", omp_get_wtime()-wall0);

    grid_free(g);
    return 0;
}
