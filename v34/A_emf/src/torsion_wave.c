/*  torsion_wave.c — Experiment 2b: Torsion wave packet propagation
 *
 *  Initialize a torsion pulse (field-space rotation) in the uniform background
 *  and track its propagation. Compare with an amplitude pulse (control).
 *
 *  Mode: -mode torsion  (default) — rotation in (phi_0, phi_1) plane
 *        -mode amplitude            — amplitude modulation phi_a *= (1+eps*env)
 *
 *  Build: gcc -O3 -march=native -fopenmp -o torsion_wave src/torsion_wave.c -lm
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>
#include <sys/stat.h>

#define NFIELDS 3
#define PI 3.14159265358979323846

/* Physics parameters */
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
    double bytes = total * sizeof(double);
    printf("Allocating %.2f GB (%ld doubles, N=%d)\n", bytes / 1e9, total, N);
    g->mem = malloc(total * sizeof(double));
    if (!g->mem) {
        fprintf(stderr, "FATAL: malloc failed for %.2f GB\n", bytes / 1e9);
        exit(1);
    }
    memset(g->mem, 0, total * sizeof(double));

    for (int a = 0; a < NFIELDS; a++) {
        g->phi[a] = g->mem + (0 + a) * g->N3;
        g->vel[a] = g->mem + (3 + a) * g->N3;
        g->acc[a] = g->mem + (6 + a) * g->N3;
    }
    return g;
}

static void grid_free(Grid *g) {
    free(g->mem);
    free(g);
}

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
        for (long idx = 0; idx < N3; idx++)
            v[idx] += hdt * ac[idx];
    }
    for (int a = 0; a < NFIELDS; a++) {
        double *p = g->phi[a], *v = g->vel[a];
        for (long idx = 0; idx < N3; idx++)
            p[idx] += dt * v[idx];
    }
    compute_forces(g);
    for (int a = 0; a < NFIELDS; a++) {
        double *v = g->vel[a], *ac = g->acc[a];
        for (long idx = 0; idx < N3; idx++)
            v[idx] += hdt * ac[idx];
    }
}

/* ================================================================
   Initialization: background + pulse
   ================================================================ */

static void init_background(Grid *g) {
    const int N = g->N, NN = N * N;
    const double dx = g->dx, L = g->L;
    const double k_bg = PI / L;
    const double omega_bg = sqrt(k_bg*k_bg + MASS2);

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            for (int kk = 0; kk < N; kk++) {
                long idx = (long)i*NN + j*N + kk;
                double z = -L + kk * dx;
                for (int a = 0; a < NFIELDS; a++) {
                    double ph = k_bg * z + 2*PI*a/3.0;
                    g->phi[a][idx] = A_BG * cos(ph);
                    g->vel[a][idx] = omega_bg * A_BG * sin(ph);
                }
            }
        }
    }
}

/* Apply torsion pulse: rotate (phi_0, phi_1) by delta_theta(x) */
static void apply_torsion_pulse(Grid *g, double x0, double sigma, double eps,
                                 double k_wave, double v_group) {
    const int N = g->N, NN = N * N;
    const double dx = g->dx, L = g->L;
    const double s2 = sigma * sigma;
    int count = 0;

    for (int i = 0; i < N; i++) {
        double x = -L + i * dx;
        double env = exp(-(x - x0)*(x - x0) / (2*s2));
        double dtheta = eps * env * cos(k_wave * x);

        /* d(dtheta)/dx for velocity correction */
        double denv_dx = -(x - x0) / s2 * env;
        double ddtheta_dx = eps * (denv_dx * cos(k_wave * x) -
                                    env * k_wave * sin(k_wave * x));

        if (fabs(dtheta) < 1e-15 && fabs(ddtheta_dx) < 1e-15) continue;

        double ct = cos(dtheta), st = sin(dtheta);
        count++;

        for (int j = 0; j < N; j++) {
            for (int kk = 0; kk < N; kk++) {
                long idx = (long)i*NN + j*N + kk;

                /* Rotate phi_0, phi_1 */
                double p0 = g->phi[0][idx];
                double p1 = g->phi[1][idx];
                g->phi[0][idx] = p0 * ct - p1 * st;
                g->phi[1][idx] = p0 * st + p1 * ct;

                /* Rotate velocities + add traveling-wave correction */
                double v0 = g->vel[0][idx];
                double v1 = g->vel[1][idx];
                g->vel[0][idx] = v0 * ct - v1 * st - v_group * ddtheta_dx * (-p0*st - p1*ct);
                g->vel[1][idx] = v0 * st + v1 * ct - v_group * ddtheta_dx * ( p0*ct - p1*st);
            }
        }
    }
    printf("Torsion pulse: x0=%.1f sigma=%.2f eps=%.3f k=%.3f v_g=%.2f (%d x-slices affected)\n",
           x0, sigma, eps, k_wave, v_group, count);
}

/* Apply amplitude pulse: phi_a *= (1 + eps * env * cos(k_wave * x)) */
static void apply_amplitude_pulse(Grid *g, double x0, double sigma, double eps,
                                   double k_wave, double v_group) {
    const int N = g->N, NN = N * N;
    const double dx = g->dx, L = g->L;
    const double s2 = sigma * sigma;
    int count = 0;

    for (int i = 0; i < N; i++) {
        double x = -L + i * dx;
        double env = exp(-(x - x0)*(x - x0) / (2*s2));
        double mod = eps * env * cos(k_wave * x);

        double denv_dx = -(x - x0) / s2 * env;
        double dmod_dx = eps * (denv_dx * cos(k_wave * x) -
                                 env * k_wave * sin(k_wave * x));

        if (fabs(mod) < 1e-15) continue;
        count++;

        for (int j = 0; j < N; j++) {
            for (int kk = 0; kk < N; kk++) {
                long idx = (long)i*NN + j*N + kk;
                for (int a = 0; a < NFIELDS; a++) {
                    /* Velocity correction for traveling wave */
                    g->vel[a][idx] += g->vel[a][idx] * mod
                                    - v_group * g->phi[a][idx] * dmod_dx;
                    g->phi[a][idx] *= (1.0 + mod);
                }
            }
        }
    }
    printf("Amplitude pulse: x0=%.1f sigma=%.2f eps=%.3f k=%.3f v_g=%.2f (%d x-slices affected)\n",
           x0, sigma, eps, k_wave, v_group, count);
}

/* ================================================================
   Diagnostics: 1D profiles along x (averaged over y,z)
   ================================================================ */

static void compute_1d_profiles(Grid *g, double *phi2_x, double *detJ_x,
                                 double *delta_P_x, int Nx) {
    const int N = g->N, NN = N * N;
    const double dx = g->dx;

    /* Compute background P for reference */
    /* P_bg = A_BG^3 * cos(kz) * cos(kz+2pi/3) * cos(kz+4pi/3) averaged over z
       The average of cos(a)cos(b)cos(c) with specific phases is complicated,
       so just compute it numerically */

    for (int i = 0; i < Nx; i++) {
        phi2_x[i] = 0;
        detJ_x[i] = 0;
        delta_P_x[i] = 0;
    }

    for (int i = 0; i < N; i++) {
        double sum_p2 = 0, sum_det = 0, sum_P = 0;
        int cnt = 0;

        for (int j = 0; j < N; j++) {
            for (int kk = 0; kk < N; kk++) {
                long idx = (long)i*NN + j*N + kk;

                /* Sigma_phi2 */
                double p2 = 0;
                for (int a = 0; a < NFIELDS; a++)
                    p2 += g->phi[a][idx] * g->phi[a][idx];
                sum_p2 += p2;

                /* P = phi_0 * phi_1 * phi_2 */
                double P = g->phi[0][idx] * g->phi[1][idx] * g->phi[2][idx];
                sum_P += P;

                /* det(J) — needs neighbors */
                int ip = (i+1)%N, im = (i-1+N)%N;
                int jp = (j+1)%N, jm = (j-1+N)%N;
                int kp = (kk+1)%N, km = (kk-1+N)%N;

                double J[3][3];
                for (int a = 0; a < NFIELDS; a++) {
                    J[a][0] = (g->phi[a][(long)ip*NN + j*N + kk] -
                               g->phi[a][(long)im*NN + j*N + kk]) / (2*dx);
                    J[a][1] = (g->phi[a][(long)i*NN + jp*N + kk] -
                               g->phi[a][(long)i*NN + jm*N + kk]) / (2*dx);
                    J[a][2] = (g->phi[a][(long)i*NN + j*N + kp] -
                               g->phi[a][(long)i*NN + j*N + km]) / (2*dx);
                }
                double det = J[0][0]*(J[1][1]*J[2][2] - J[1][2]*J[2][1])
                           - J[0][1]*(J[1][0]*J[2][2] - J[1][2]*J[2][0])
                           + J[0][2]*(J[1][0]*J[2][1] - J[1][1]*J[2][0]);
                sum_det += det;
                cnt++;
            }
        }

        phi2_x[i] = sum_p2 / cnt;
        detJ_x[i] = sum_det / cnt;
        delta_P_x[i] = sum_P / cnt;
    }
}

/* ================================================================
   Energy diagnostics
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

/* ================================================================
   Main
   ================================================================ */

int main(int argc, char **argv) {
    int N = 128;
    double L = 30.0, T = 60.0;
    double diag_dt = 1.0, prof_dt = 5.0;
    char mode[32] = "torsion";
    char outdir[256] = "data/torsion_wave";

    /* Pulse parameters */
    double pulse_x0 = -15.0;
    double pulse_sigma = 2.0 * sqrt(2.0);  /* sigma = 2*sqrt(2) */
    double pulse_eps = 0.05;
    double pulse_k = 2*PI / 4.0;   /* wavelength = 4 */
    double pulse_vg = 1.0;          /* group velocity guess */

    for (int i = 1; i < argc; i++) {
        if      (!strcmp(argv[i],"-N"))      N = atoi(argv[++i]);
        else if (!strcmp(argv[i],"-L"))      L = atof(argv[++i]);
        else if (!strcmp(argv[i],"-T"))      T = atof(argv[++i]);
        else if (!strcmp(argv[i],"-mode"))   strncpy(mode, argv[++i], 31);
        else if (!strcmp(argv[i],"-o"))      strncpy(outdir, argv[++i], 255);
        else if (!strcmp(argv[i],"-diag"))   diag_dt = atof(argv[++i]);
        else if (!strcmp(argv[i],"-prof"))   prof_dt = atof(argv[++i]);
        else if (!strcmp(argv[i],"-eps"))    pulse_eps = atof(argv[++i]);
        else if (!strcmp(argv[i],"-x0"))     pulse_x0 = atof(argv[++i]);
        else if (!strcmp(argv[i],"-vg"))     pulse_vg = atof(argv[++i]);
        else if (!strcmp(argv[i],"-kw"))     pulse_k = 2*PI / atof(argv[++i]);
    }

    {
        char *env_threads = getenv("OMP_NUM_THREADS");
        int nthreads = env_threads ? atoi(env_threads) : 4;
        omp_set_num_threads(nthreads);
    }

    printf("=== Exp 2b: Torsion Wave Propagation ===\n");
    printf("Mode: %s, N=%d L=%.0f T=%.0f\n", mode, N, L, T);

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

    /* Initialize background */
    init_background(g);

    /* Apply pulse */
    if (strcmp(mode, "torsion") == 0) {
        apply_torsion_pulse(g, pulse_x0, pulse_sigma, pulse_eps, pulse_k, pulse_vg);
    } else if (strcmp(mode, "amplitude") == 0) {
        apply_amplitude_pulse(g, pulse_x0, pulse_sigma, pulse_eps, pulse_k, pulse_vg);
    } else {
        fprintf(stderr, "Unknown mode: %s\n", mode);
        return 1;
    }

    compute_forces(g);

    /* Output files */
    char tspath[512], profpath[512];
    snprintf(tspath, sizeof(tspath), "%s/timeseries.tsv", outdir);
    snprintf(profpath, sizeof(profpath), "%s/profiles.tsv", outdir);
    FILE *fp_ts = fopen(tspath, "w");
    FILE *fp_prof = fopen(profpath, "w");
    fprintf(fp_ts, "t\tE_kin\tE_grad\tE_mass\tE_pot\tE_total\n");
    fprintf(fp_prof, "# 1D profiles along x (averaged over y,z)\n");
    fprintf(fp_prof, "t\tx\tphi2\tdetJ\tP\n");

    int n_steps  = (int)(T / g->dt);
    int diag_every = (int)(diag_dt / g->dt); if (diag_every < 1) diag_every = 1;
    int prof_every = (int)(prof_dt / g->dt); if (prof_every < 1) prof_every = 1;

    printf("Steps=%d diag_every=%d prof_every=%d\n\n", n_steps, diag_every, prof_every);

    double *phi2_x   = calloc(N, sizeof(double));
    double *detJ_x   = calloc(N, sizeof(double));
    double *deltaP_x = calloc(N, sizeof(double));

    double wall0 = omp_get_wtime();
    double E0 = 0;

    for (int step = 0; step <= n_steps; step++) {
        if (step > 0) verlet_step(g);
        double t = step * g->dt;

        if (step % diag_every == 0) {
            double ek, eg, em, ep, et;
            compute_energy(g, &ek, &eg, &em, &ep, &et);
            if (step == 0) E0 = et;
            fprintf(fp_ts, "%.3f\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\n",
                    t, ek, eg, em, ep, et);
            fflush(fp_ts);

            if (step % (diag_every * 20) == 0) {
                double wall = omp_get_wtime() - wall0;
                double pct = 100.0 * step / n_steps;
                double drift = 100.0 * (et - E0) / (fabs(E0) + 1e-30);
                printf("t=%6.1f E=%.4e (drift %+.4f%%) [%.0f%% %.0fs]\n",
                       t, et, drift, pct, wall);
            }
        }

        if (step % prof_every == 0) {
            compute_1d_profiles(g, phi2_x, detJ_x, deltaP_x, N);
            for (int i = 0; i < N; i++) {
                double x = -g->L + i * g->dx;
                fprintf(fp_prof, "%.3f\t%.4f\t%.8e\t%.8e\t%.8e\n",
                        t, x, phi2_x[i], detJ_x[i], deltaP_x[i]);
            }
            fflush(fp_prof);
        }
    }

    fclose(fp_ts);
    fclose(fp_prof);

    double wall = omp_get_wtime() - wall0;
    printf("\n=== Complete: %.0fs (%.1f min) ===\n", wall, wall/60);

    free(phi2_x); free(detJ_x); free(deltaP_x);
    grid_free(g);
    return 0;
}
