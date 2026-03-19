/*  v33_drag_test.c — Measure effective drag on a kicked braid in different backgrounds
 *
 *  Hypothesis: the energy cost of braid motion is proportional to ρ (field density),
 *  NOT to distance. So a braid decelerates SLOWER in low-ρ backgrounds.
 *
 *  Method:
 *    1. Initialize single braid in uniform background of amplitude A_bg
 *    2. Apply a velocity kick: vel[a] -= v_kick × ∂φ_a/∂x (Galilean boost)
 *    3. Track braid center-of-mass position at high frequency
 *    4. Fit velocity decay: v(t) = v0 × exp(-t/τ)
 *    5. Compare drag coefficient 1/τ across different A_bg values
 *
 *  Run:  ./drag -bg 0.15 -vk 0.1 -o data/drag_high &
 *        ./drag -bg 0.05 -vk 0.1 -o data/drag_low &
 *        ./drag -bg 0.10 -vk 0.1 -o data/drag_mid &
 *
 *  Build: gcc -O3 -march=native -fopenmp -o drag src/v33_drag_test.c -lm
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>
#include <sys/stat.h>

#define NFIELDS 3
#define PI 3.14159265358979323846

static double MU    = -41.345;
static double KAPPA = 50.0;
static double MASS2 = 2.25;
static double A_BG  = 0.1;

/* ================================================================
   Grid: ONE allocation
   ================================================================ */

typedef struct {
    double *mem;
    double *phi[NFIELDS], *vel[NFIELDS], *acc[NFIELDS];
    int N; long N3;
    double L, dx, dt;
} Grid;

static Grid *grid_alloc(int N, double L) {
    Grid *g = calloc(1, sizeof(Grid));
    g->N = N; g->N3 = (long)N*N*N;
    g->L = L; g->dx = 2.0*L/(N-1);
    g->dt = 0.12 * g->dx;
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
    return g;
}

static void grid_free(Grid *g) { free(g->mem); free(g); }

/* ================================================================
   Forces: standard equation, periodic BC
   ================================================================ */

static void compute_forces(Grid *g) {
    const int N = g->N, NN = N*N;
    const long N3 = g->N3;
    const double idx2 = 1.0/(g->dx*g->dx);

    #pragma omp parallel for schedule(static)
    for (long idx = 0; idx < N3; idx++) {
        int i = (int)(idx/NN), j = (int)((idx/N)%N), k = (int)(idx%N);
        int ip=(i+1)%N, im=(i-1+N)%N;
        int jp=(j+1)%N, jm=(j-1+N)%N;
        int kp=(k+1)%N, km=(k-1+N)%N;

        long n_ip=(long)ip*NN+j*N+k, n_im=(long)im*NN+j*N+k;
        long n_jp=(long)i*NN+jp*N+k, n_jm=(long)i*NN+jm*N+k;
        long n_kp=(long)i*NN+j*N+kp, n_km=(long)i*NN+j*N+km;

        double p0=g->phi[0][idx], p1=g->phi[1][idx], p2=g->phi[2][idx];
        double P = p0*p1*p2;
        double den = 1.0 + KAPPA*P*P;
        double mPd2 = MU*P/(den*den);

        for (int a = 0; a < NFIELDS; a++) {
            double lap = (g->phi[a][n_ip]+g->phi[a][n_im]
                        +g->phi[a][n_jp]+g->phi[a][n_jm]
                        +g->phi[a][n_kp]+g->phi[a][n_km]
                        -6.0*g->phi[a][idx]) * idx2;
            double dPda = (a==0)?p1*p2:(a==1)?p0*p2:p0*p1;
            g->acc[a][idx] = lap - MASS2*g->phi[a][idx] - mPd2*dPda;
        }
    }
}

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
    compute_forces(g);
    for (int a=0;a<NFIELDS;a++) {
        double *v=g->vel[a], *ac=g->acc[a];
        for (long idx=0;idx<N3;idx++) v[idx] += hdt*ac[idx];
    }
}

/* ================================================================
   Initialize: braid + uniform background + velocity kick
   ================================================================ */

static void init_kicked_braid(Grid *g, double v_kick) {
    const int N = g->N, NN = N*N;
    const double dx = g->dx, L = g->L;
    const double A[3] = {0.8, 0.8, 0.8};
    const double delta[3] = {0, 3.0005, 4.4325};
    const double R_tube = 3.0, ellip = 0.3325;
    const double kw = PI/L;
    const double omega = sqrt(kw*kw + MASS2);
    const double sx = 1+ellip, sy = 1-ellip;
    const double inv2R2 = 1.0/(2*R_tube*R_tube);
    const double k_bg = PI/L;
    const double omega_bg = sqrt(k_bg*k_bg + MASS2);

    /* Standard init */
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
                    g->vel[a][idx] = omega*A[a]*env*sin(ph) + omega_bg*A_BG*sin(ph_bg);
                }
            }
        }
    }

    /* Apply velocity kick: vel[a] -= v_kick × ∂φ_a/∂x
     * Only within the braid envelope (smooth cutoff) */
    printf("Applying velocity kick: v_kick = %.4f\n", v_kick);
    double kick_energy = 0;
    for (int i = 0; i < N; i++) {
        double x = -L + i*dx;
        for (int j = 0; j < N; j++) {
            double y = -L + j*dx;
            for (int kk = 0; kk < N; kk++) {
                long idx = (long)i*NN + j*N + kk;

                /* Envelope weight: only kick the braid region */
                double r2e = x*x/((1+ellip)*(1+ellip)) + y*y/((1-ellip)*(1-ellip));
                double w = exp(-r2e * inv2R2);
                if (w < 0.01) continue;  /* skip far-field */

                for (int a = 0; a < NFIELDS; a++) {
                    /* ∂φ_a/∂x via central difference */
                    int ip = (i+1)%N, im = (i-1+N)%N;
                    double dpdx = (g->phi[a][(long)ip*NN+j*N+kk]
                                 - g->phi[a][(long)im*NN+j*N+kk]) / (2*dx);
                    double dv = -v_kick * w * dpdx;
                    g->vel[a][idx] += dv;
                    kick_energy += dv * dv * dx*dx*dx;
                }
            }
        }
    }
    printf("Kick energy injected: %.4e\n", 0.5*kick_energy);
}

/* ================================================================
   Diagnostics: braid position and velocity
   ================================================================ */

static void measure_braid(Grid *g, double *cx, double *cy, double *cz,
                          double *px, double *braid_E) {
    const long N3 = g->N3;
    const int NN = g->N*g->N, N = g->N;
    const double dx = g->dx, L = g->L, dV = dx*dx*dx;

    /* Average for threshold */
    double avg = 0;
    for (long idx = 0; idx < N3; idx++) {
        double p2 = 0;
        for (int a=0;a<NFIELDS;a++) p2 += g->phi[a][idx]*g->phi[a][idx];
        avg += p2;
    }
    avg /= N3;
    double thresh = 5.0 * avg;

    /* Centroid + x-momentum */
    double wx=0, wy=0, wz=0, wt=0, mom_x=0, eB=0;
    for (long idx = 0; idx < N3; idx++) {
        double p2 = 0;
        for (int a=0;a<NFIELDS;a++) p2 += g->phi[a][idx]*g->phi[a][idx];
        if (p2 < thresh) continue;

        int i = (int)(idx/NN), j = (int)((idx/N)%N), k = (int)(idx%N);
        double x = -L + i*dx, y = -L + j*dx, z = -L + k*dx;
        wx += x*p2; wy += y*p2; wz += z*p2; wt += p2;

        /* x-momentum: Σ_a vel[a] × ∂φ_a/∂x */
        for (int a = 0; a < NFIELDS; a++) {
            int ip=(i+1)%N, im=(i-1+N)%N;
            double dpdx = (g->phi[a][(long)ip*NN+j*N+k]
                         - g->phi[a][(long)im*NN+j*N+k]) / (2*dx);
            mom_x -= g->vel[a][idx] * dpdx * dV;  /* field momentum: -Σ vel·∂φ/∂x */
        }

        /* Braid energy (kinetic only, quick proxy) */
        for (int a=0;a<NFIELDS;a++)
            eB += 0.5*g->vel[a][idx]*g->vel[a][idx]*dV;
    }

    *cx = (wt > 0) ? wx/wt : 0;
    *cy = (wt > 0) ? wy/wt : 0;
    *cz = (wt > 0) ? wz/wt : 0;
    *px = mom_x;
    *braid_E = eB;
}

static void compute_energy(Grid *g, double *E_total) {
    const int N = g->N, NN = N*N;
    const long N3 = g->N3;
    const double dx = g->dx, dV = dx*dx*dx;
    double et = 0;

    #pragma omp parallel for reduction(+:et) schedule(static)
    for (long idx = 0; idx < N3; idx++) {
        int i=(int)(idx/NN), j=(int)((idx/N)%N), k=(int)(idx%N);
        int ip=(i+1)%N, im=(i-1+N)%N;
        int jp=(j+1)%N, jm=(j-1+N)%N;
        int kp=(k+1)%N, km=(k-1+N)%N;
        for (int a = 0; a < NFIELDS; a++) {
            et += 0.5*g->vel[a][idx]*g->vel[a][idx]*dV;
            double gx=(g->phi[a][(long)ip*NN+j*N+k]-g->phi[a][(long)im*NN+j*N+k])/(2*dx);
            double gy=(g->phi[a][(long)i*NN+jp*N+k]-g->phi[a][(long)i*NN+jm*N+k])/(2*dx);
            double gz=(g->phi[a][(long)i*NN+j*N+kp]-g->phi[a][(long)i*NN+j*N+km])/(2*dx);
            et += 0.5*(gx*gx+gy*gy+gz*gz)*dV;
            et += 0.5*MASS2*g->phi[a][idx]*g->phi[a][idx]*dV;
        }
        double P = g->phi[0][idx]*g->phi[1][idx]*g->phi[2][idx];
        et += (MU/2.0)*P*P/(1.0+KAPPA*P*P)*dV;
    }
    *E_total = et;
}

/* ================================================================
   Main
   ================================================================ */

int main(int argc, char **argv) {
    int N = 128; double L = 20.0, T = 200.0;
    double v_kick = 0.1;
    char outdir[256] = "data/drag";

    for (int i=1;i<argc;i++) {
        if      (!strcmp(argv[i],"-N"))  N = atoi(argv[++i]);
        else if (!strcmp(argv[i],"-L"))  L = atof(argv[++i]);
        else if (!strcmp(argv[i],"-T"))  T = atof(argv[++i]);
        else if (!strcmp(argv[i],"-bg")) A_BG = atof(argv[++i]);
        else if (!strcmp(argv[i],"-vk")) v_kick = atof(argv[++i]);
        else if (!strcmp(argv[i],"-o"))  strncpy(outdir, argv[++i], 255);
    }

    omp_set_num_threads(8);  /* 8 per run, two in parallel */

    printf("=== V33 Drag Test: Kicked Braid in Uniform Background ===\n");
    printf("A_bg=%.3f  v_kick=%.3f  N=%d  L=%.0f  T=%.0f\n", A_BG, v_kick, N, L, T);
    printf("ρ_bg ∝ A_bg² = %.4f\n\n", A_BG*A_BG);

    mkdir("data", 0755); mkdir(outdir, 0755);
    Grid *g = grid_alloc(N, L);
    printf("dx=%.4f dt=%.5f\n\n", g->dx, g->dt);

    init_kicked_braid(g, v_kick);
    compute_forces(g);

    /* Output */
    char tspath[512];
    snprintf(tspath, sizeof(tspath), "%s/timeseries.tsv", outdir);
    FILE *fp = fopen(tspath, "w");
    fprintf(fp, "t\tx_braid\ty_braid\tp_x\tv_x_est\tE_braid\tE_total\n");

    int n_steps = (int)(T/g->dt);
    int diag_every = (int)(0.5/g->dt);  /* every 0.5 time units — high frequency */
    if (diag_every < 1) diag_every = 1;

    printf("Steps=%d  diag_every=%d (dt_diag=%.2f)\n\n", n_steps, diag_every, diag_every*g->dt);

    double wall0 = omp_get_wtime();
    double E0 = 0;
    double prev_x = 0, prev_t = 0;

    for (int step = 0; step <= n_steps; step++) {
        if (step > 0) verlet_step(g);
        double t = step * g->dt;

        if (step % diag_every == 0) {
            double cx, cy, cz, px, eB, et;
            measure_braid(g, &cx, &cy, &cz, &px, &eB);
            compute_energy(g, &et);
            if (step == 0) { E0 = et; prev_x = cx; prev_t = t; }

            /* Velocity estimate from position difference */
            double vx_est = 0;
            if (step > 0 && (t - prev_t) > 0)
                vx_est = (cx - prev_x) / (t - prev_t);

            fprintf(fp, "%.3f\t%.6f\t%.6f\t%.6e\t%.6e\t%.4e\t%.4e\n",
                    t, cx, cy, px, vx_est, eB, et);
            fflush(fp);

            if (step % (diag_every * 20) == 0) {
                double wall = omp_get_wtime() - wall0;
                double drift_pct = 100.0*(et-E0)/(fabs(E0)+1e-30);
                printf("t=%7.1f  x=%+.4f  p_x=%+.3e  v_x=%+.4e  E_drift=%+.2f%%  [%.0f%% %.0fs]\n",
                       t, cx, px, vx_est, drift_pct, 100.0*step/n_steps, wall);
                fflush(stdout);
            }
            prev_x = cx; prev_t = t;
        }
    }

    fclose(fp);
    double cx, cy, cz, px, eB;
    measure_braid(g, &cx, &cy, &cz, &px, &eB);

    printf("\n=== RESULT ===\n");
    printf("A_bg=%.3f  v_kick=%.3f\n", A_BG, v_kick);
    printf("Final x_braid = %.4f  (total displacement from ~0)\n", cx);
    printf("Final p_x     = %.4e\n", px);
    printf("Run time: %.0fs\n", omp_get_wtime() - wall0);

    printf("\nCompare across A_bg values:\n");
    printf("  If displacement is LARGER at low A_bg → drag is LOWER in depleted field ✓\n");
    printf("  If displacement is SAME → drag is independent of ρ\n");
    printf("  If displacement is SMALLER at low A_bg → drag is HIGHER in depleted field\n");

    grid_free(g);
    return 0;
}
