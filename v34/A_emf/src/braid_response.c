/*  braid_response.c — Experiments 2c+2d: Braid response to torsion vs amplitude waves
 *
 *  Phase 1: Initialize braid, settle for T_settle time
 *  Phase 2: Apply pulse (torsion or amplitude), run and track braid position
 *
 *  Modes:
 *    -mode torsion    — field-space rotation pulse
 *    -mode amplitude  — amplitude modulation pulse
 *  Winding:
 *    -winding +1  (default) — standard braid delta = (0, +3.0005, +4.4325)
 *    -winding -1             — reversed delta = (0, -3.0005, -4.4325)
 *
 *  Build: gcc -O3 -march=native -fopenmp -o braid_response src/braid_response.c -lm
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

typedef struct {
    double *mem;
    double *phi[NFIELDS];
    double *vel[NFIELDS];
    double *acc[NFIELDS];
    int N;
    long N3;
    double L, dx, dt;
} Grid;

static Grid *grid_alloc(int N, double L) {
    Grid *g = calloc(1, sizeof(Grid));
    g->N  = N;
    g->N3 = (long)N * N * N;
    g->L  = L;
    g->dx = 2.0 * L / (N - 1);
    g->dt = 0.12 * g->dx;

    long total = 9 * g->N3;
    g->mem = malloc(total * sizeof(double));
    if (!g->mem) { fprintf(stderr, "FATAL: malloc failed\n"); exit(1); }
    memset(g->mem, 0, total * sizeof(double));

    for (int a = 0; a < NFIELDS; a++) {
        g->phi[a] = g->mem + (0 + a) * g->N3;
        g->vel[a] = g->mem + (3 + a) * g->N3;
        g->acc[a] = g->mem + (6 + a) * g->N3;
    }
    return g;
}

static void grid_free(Grid *g) { free(g->mem); free(g); }

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

        double p0 = g->phi[0][idx];
        double p1 = g->phi[1][idx];
        double p2 = g->phi[2][idx];
        double P  = p0 * p1 * p2;
        double den = 1.0 + KAPPA * P * P;
        double mPd2 = MU * P / (den * den);

        for (int a = 0; a < NFIELDS; a++) {
            double lap = (g->phi[a][n_ip] + g->phi[a][n_im]
                        + g->phi[a][n_jp] + g->phi[a][n_jm]
                        + g->phi[a][n_kp] + g->phi[a][n_km]
                        - 6.0 * g->phi[a][idx]) * idx2;
            double dPda = (a==0) ? p1*p2 : (a==1) ? p0*p2 : p0*p1;
            g->acc[a][idx] = lap - MASS2 * g->phi[a][idx] - mPd2 * dPda;
        }
    }
}

static void verlet_step(Grid *g) {
    const long N3 = g->N3;
    const double hdt = 0.5 * g->dt;
    const double dt  = g->dt;

    for (int a = 0; a < NFIELDS; a++) {
        double *v = g->vel[a], *ac = g->acc[a];
        for (long idx = 0; idx < N3; idx++) v[idx] += hdt * ac[idx];
    }
    for (int a = 0; a < NFIELDS; a++) {
        double *p = g->phi[a], *v = g->vel[a];
        for (long idx = 0; idx < N3; idx++) p[idx] += dt * v[idx];
    }
    compute_forces(g);
    for (int a = 0; a < NFIELDS; a++) {
        double *v = g->vel[a], *ac = g->acc[a];
        for (long idx = 0; idx < N3; idx++) v[idx] += hdt * ac[idx];
    }
}

/* ================================================================
   Initialization
   ================================================================ */

static void init_braid(Grid *g, double x_cen, double y_cen, int winding) {
    const int N = g->N, NN = N * N;
    const double dx = g->dx, L = g->L;
    const double A[3] = {0.8, 0.8, 0.8};
    double delta[3] = {0, 3.0005, 4.4325};
    if (winding < 0) { delta[1] = -delta[1]; delta[2] = -delta[2]; }

    const double R_tube = 3.0, ellip = 0.3325;
    const double kw = PI / L;
    const double omega = sqrt(kw*kw + MASS2);
    const double sx = 1+ellip, sy = 1-ellip;
    const double inv2R2 = 1.0 / (2*R_tube*R_tube);
    const double k_bg = PI / L;
    const double omega_bg = sqrt(k_bg*k_bg + MASS2);

    for (int i = 0; i < N; i++) {
        double x = -L + i*dx;
        for (int j = 0; j < N; j++) {
            double y = -L + j*dx;
            for (int kk = 0; kk < N; kk++) {
                double z = -L + kk*dx;
                long idx = (long)i*NN + j*N + kk;

                double xc = x - x_cen, yc = y - y_cen;
                double r2e = xc*xc/(sx*sx) + yc*yc/(sy*sy);
                double env = exp(-r2e * inv2R2);

                for (int a = 0; a < NFIELDS; a++) {
                    double ph = kw*z + delta[a];
                    double ph_bg = k_bg*z + 2*PI*a/3.0;
                    g->phi[a][idx] = A[a]*env*cos(ph) + A_BG*cos(ph_bg);
                    g->vel[a][idx] = omega*A[a]*env*sin(ph)
                                   + omega_bg*A_BG*sin(ph_bg);
                }
            }
        }
    }
    printf("Braid initialized: center=(%.1f,%.1f) winding=%+d\n",
           x_cen, y_cen, winding);
    printf("  delta = (%.4f, %.4f, %.4f)\n", delta[0], delta[1], delta[2]);
}

/* ================================================================
   Pulse application (same as torsion_wave.c)
   ================================================================ */

static void apply_torsion_pulse(Grid *g, double x0, double sigma, double eps,
                                 double k_wave, double v_group) {
    const int N = g->N, NN = N * N;
    const double dx = g->dx, L = g->L;
    const double s2 = sigma * sigma;

    for (int i = 0; i < N; i++) {
        double x = -L + i * dx;
        double env = exp(-(x - x0)*(x - x0) / (2*s2));
        double dtheta = eps * env * cos(k_wave * x);
        double denv_dx = -(x - x0) / s2 * env;
        double ddtheta_dx = eps * (denv_dx * cos(k_wave * x) -
                                    env * k_wave * sin(k_wave * x));

        if (fabs(dtheta) < 1e-15 && fabs(ddtheta_dx) < 1e-15) continue;
        double ct = cos(dtheta), st = sin(dtheta);

        for (int j = 0; j < N; j++) {
            for (int kk = 0; kk < N; kk++) {
                long idx = (long)i*NN + j*N + kk;

                double p0 = g->phi[0][idx];
                double p1 = g->phi[1][idx];
                g->phi[0][idx] = p0 * ct - p1 * st;
                g->phi[1][idx] = p0 * st + p1 * ct;

                double v0 = g->vel[0][idx];
                double v1 = g->vel[1][idx];
                g->vel[0][idx] = v0 * ct - v1 * st - v_group * ddtheta_dx * (-p0*st - p1*ct);
                g->vel[1][idx] = v0 * st + v1 * ct - v_group * ddtheta_dx * ( p0*ct - p1*st);
            }
        }
    }
    printf("Applied torsion pulse: x0=%.1f sigma=%.2f eps=%.3f\n", x0, sigma, eps);
}

static void apply_amplitude_pulse(Grid *g, double x0, double sigma, double eps,
                                   double k_wave, double v_group) {
    const int N = g->N, NN = N * N;
    const double dx = g->dx, L = g->L;
    const double s2 = sigma * sigma;

    for (int i = 0; i < N; i++) {
        double x = -L + i * dx;
        double env = exp(-(x - x0)*(x - x0) / (2*s2));
        double mod = eps * env * cos(k_wave * x);
        double denv_dx = -(x - x0) / s2 * env;
        double dmod_dx = eps * (denv_dx * cos(k_wave * x) -
                                 env * k_wave * sin(k_wave * x));

        if (fabs(mod) < 1e-15) continue;

        for (int j = 0; j < N; j++) {
            for (int kk = 0; kk < N; kk++) {
                long idx = (long)i*NN + j*N + kk;
                for (int a = 0; a < NFIELDS; a++) {
                    g->vel[a][idx] += g->vel[a][idx] * mod
                                    - v_group * g->phi[a][idx] * dmod_dx;
                    g->phi[a][idx] *= (1.0 + mod);
                }
            }
        }
    }
    printf("Applied amplitude pulse: x0=%.1f sigma=%.2f eps=%.3f\n", x0, sigma, eps);
}

/* ================================================================
   Diagnostics
   ================================================================ */

static void compute_energy(Grid *g, double *E_kin, double *E_grad,
                           double *E_mass, double *E_pot, double *E_total) {
    const int N = g->N, NN = N*N;
    const long N3 = g->N3;
    const double dx = g->dx, dV = dx*dx*dx;
    double ek=0, eg=0, em=0, ep=0;

    #pragma omp parallel for reduction(+:ek,eg,em,ep) schedule(static)
    for (long idx = 0; idx < N3; idx++) {
        int i=(int)(idx/NN), j=(int)((idx/N)%N), k=(int)(idx%N);
        int ip=(i+1)%N,im=(i-1+N)%N,jp=(j+1)%N,jm=(j-1+N)%N;
        int kp=(k+1)%N,km=(k-1+N)%N;

        for (int a = 0; a < NFIELDS; a++) {
            ek += 0.5 * g->vel[a][idx] * g->vel[a][idx] * dV;
            double gx = (g->phi[a][(long)ip*NN+j*N+k]-g->phi[a][(long)im*NN+j*N+k])/(2*dx);
            double gy = (g->phi[a][(long)i*NN+jp*N+k]-g->phi[a][(long)i*NN+jm*N+k])/(2*dx);
            double gz = (g->phi[a][(long)i*NN+j*N+kp]-g->phi[a][(long)i*NN+j*N+km])/(2*dx);
            eg += 0.5*(gx*gx+gy*gy+gz*gz) * dV;
            em += 0.5*MASS2*g->phi[a][idx]*g->phi[a][idx] * dV;
        }
        double P = g->phi[0][idx]*g->phi[1][idx]*g->phi[2][idx];
        ep += (MU/2.0)*P*P/(1.0+KAPPA*P*P) * dV;
    }
    *E_kin=ek; *E_grad=eg; *E_mass=em; *E_pot=ep; *E_total=ek+eg+em+ep;
}

/* Find braid centroid (energy-weighted above threshold) */
static void find_braid_center(Grid *g, double *bx, double *by, double *bz) {
    const int N = g->N, NN = N * N;
    const long N3 = g->N3;
    const double dx = g->dx, L = g->L;

    double avg = 0;
    for (long idx = 0; idx < N3; idx++) {
        double p2 = 0;
        for (int a = 0; a < NFIELDS; a++)
            p2 += g->phi[a][idx] * g->phi[a][idx];
        avg += p2;
    }
    avg /= N3;
    double thresh = 5.0 * avg;

    double wx = 0, wy = 0, wz = 0, wt = 0;
    for (long idx = 0; idx < N3; idx++) {
        double p2 = 0;
        for (int a = 0; a < NFIELDS; a++)
            p2 += g->phi[a][idx] * g->phi[a][idx];
        if (p2 < thresh) continue;

        int i = (int)(idx / NN);
        int j = (int)((idx / N) % N);
        int k = (int)(idx % N);
        double x = -L + i * dx;
        double y = -L + j * dx;
        double z = -L + k * dx;

        wx += x * p2; wy += y * p2; wz += z * p2; wt += p2;
    }
    *bx = (wt > 0) ? wx/wt : 0;
    *by = (wt > 0) ? wy/wt : 0;
    *bz = (wt > 0) ? wz/wt : 0;
}

/* Compute peak amplitude (max Sigma_phi2) */
static double peak_amplitude(Grid *g) {
    const long N3 = g->N3;
    double maxp2 = 0;
    for (long idx = 0; idx < N3; idx++) {
        double p2 = 0;
        for (int a = 0; a < NFIELDS; a++)
            p2 += g->phi[a][idx] * g->phi[a][idx];
        if (p2 > maxp2) maxp2 = p2;
    }
    return maxp2;
}

/* ================================================================
   Main
   ================================================================ */

int main(int argc, char **argv) {
    int N = 128;
    double L = 30.0;
    double T_settle = 50.0;
    double T_run = 200.0;
    double diag_dt = 1.0;
    char mode[32] = "torsion";
    char outdir[256] = "data/braid_torsion";
    int winding = 1;

    /* Pulse parameters */
    double pulse_x0 = -15.0;
    double pulse_sigma = 2.0 * sqrt(2.0);
    double pulse_eps = 0.05;
    double pulse_k = 2*PI / 4.0;
    double pulse_vg = 1.0;

    for (int i = 1; i < argc; i++) {
        if      (!strcmp(argv[i],"-N"))       N = atoi(argv[++i]);
        else if (!strcmp(argv[i],"-L"))       L = atof(argv[++i]);
        else if (!strcmp(argv[i],"-Tsettle")) T_settle = atof(argv[++i]);
        else if (!strcmp(argv[i],"-Trun"))    T_run = atof(argv[++i]);
        else if (!strcmp(argv[i],"-mode"))    strncpy(mode, argv[++i], 31);
        else if (!strcmp(argv[i],"-o"))       strncpy(outdir, argv[++i], 255);
        else if (!strcmp(argv[i],"-diag"))    diag_dt = atof(argv[++i]);
        else if (!strcmp(argv[i],"-winding")) winding = atoi(argv[++i]);
        else if (!strcmp(argv[i],"-eps"))     pulse_eps = atof(argv[++i]);
        else if (!strcmp(argv[i],"-x0"))      pulse_x0 = atof(argv[++i]);
    }

    {
        char *env_threads = getenv("OMP_NUM_THREADS");
        int nthreads = env_threads ? atoi(env_threads) : 4;
        omp_set_num_threads(nthreads);
    }

    printf("=== Braid Response: mode=%s winding=%+d ===\n", mode, winding);
    printf("N=%d L=%.0f T_settle=%.0f T_run=%.0f\n", N, L, T_settle, T_run);

    /* Create output directory */
    {
        char tmp[512];
        strncpy(tmp, outdir, sizeof(tmp));
        for (char *p = tmp + 1; *p; p++) {
            if (*p == '/') { *p = '\0'; mkdir(tmp, 0755); *p = '/'; }
        }
        mkdir(tmp, 0755);
    }

    Grid *g = grid_alloc(N, L);
    printf("dx=%.4f dt=%.5f\n\n", g->dx, g->dt);

    /* Phase 1: Initialize and settle braid */
    printf("=== Phase 1: Settling braid (T=%.0f) ===\n", T_settle);
    init_braid(g, 0, 0, winding);
    compute_forces(g);

    int settle_steps = (int)(T_settle / g->dt);
    int diag_every = (int)(diag_dt / g->dt); if (diag_every < 1) diag_every = 1;

    double wall0 = omp_get_wtime();

    for (int step = 0; step < settle_steps; step++) {
        verlet_step(g);
        if (step % (settle_steps/5) == 0) {
            double t = step * g->dt;
            double bx, by, bz;
            find_braid_center(g, &bx, &by, &bz);
            double pk = peak_amplitude(g);
            printf("  settle t=%.1f center=(%.3f,%.3f,%.3f) peak_phi2=%.4f\n",
                   t, bx, by, bz, pk);
        }
    }

    double bx0, by0, bz0;
    find_braid_center(g, &bx0, &by0, &bz0);
    printf("Settled: center=(%.3f, %.3f, %.3f)\n\n", bx0, by0, bz0);

    /* Phase 2: Apply pulse and run */
    printf("=== Phase 2: Apply %s pulse, run T=%.0f ===\n", mode, T_run);

    if (strcmp(mode, "torsion") == 0) {
        apply_torsion_pulse(g, pulse_x0, pulse_sigma, pulse_eps, pulse_k, pulse_vg);
    } else if (strcmp(mode, "amplitude") == 0) {
        apply_amplitude_pulse(g, pulse_x0, pulse_sigma, pulse_eps, pulse_k, pulse_vg);
    } else {
        fprintf(stderr, "Unknown mode: %s\n", mode);
        return 1;
    }
    compute_forces(g);

    /* Output */
    char tspath[512];
    snprintf(tspath, sizeof(tspath), "%s/tracking.tsv", outdir);
    FILE *fp = fopen(tspath, "w");
    fprintf(fp, "t\tbraid_x\tbraid_y\tbraid_z\tpeak_phi2\tE_total\tE_kin\tE_pot\n");

    int run_steps = (int)(T_run / g->dt);
    printf("Run: %d steps, diag every %d\n\n", run_steps, diag_every);

    for (int step = 0; step <= run_steps; step++) {
        if (step > 0) verlet_step(g);
        double t = step * g->dt;

        if (step % diag_every == 0) {
            double bx, by, bz;
            find_braid_center(g, &bx, &by, &bz);
            double pk = peak_amplitude(g);
            double ek, eg, em, ep, et;
            compute_energy(g, &ek, &eg, &em, &ep, &et);

            fprintf(fp, "%.3f\t%.6f\t%.6f\t%.6f\t%.6e\t%.6e\t%.6e\t%.6e\n",
                    t, bx, by, bz, pk, et, ek, ep);
            fflush(fp);

            if (step % (diag_every * 50) == 0) {
                double wall = omp_get_wtime() - wall0;
                printf("t=%6.1f braid=(%.3f,%.3f,%.3f) dx=%.4f peak=%.4f [%.0fs]\n",
                       t, bx, by, bz, bx-bx0, pk, wall);
            }
        }
    }

    fclose(fp);
    double wall = omp_get_wtime() - wall0;
    printf("\n=== Complete: %.0fs (%.1f min) ===\n", wall, wall/60);
    printf("Output: %s\n", tspath);

    grid_free(g);
    return 0;
}
