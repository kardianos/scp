/*  v33_cosserat_sfa.c — 6-field Cosserat with SFA archive output
 *
 *  Same physics as v33_cosserat.c, but writes a single .sfa file
 *  with every frame compressed (zstd), seekable, and self-describing.
 *
 *  Build: gcc -O3 -march=native -fopenmp -o cosserat_sfa \
 *         src/v33_cosserat_sfa.c -lzstd -lm
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>
#include <sys/stat.h>

#define SFA_IMPLEMENTATION
#include "../viewer/sfa.h"

#define NFIELDS 3
#define PI 3.14159265358979323846

/* Physics parameters */
static double MU     = -41.345;
static double KAPPA  = 50.0;
static double MASS2  = 2.25;    /* position field mass */
static double MTHETA2 = 2.25;  /* angle field mass (can differ) */
static double A_BG   = 0.1;
static double ETA    = 1.0;    /* position-angle coupling strength */

/* ================================================================
   Grid: ONE allocation for everything (18 arrays)
   ================================================================ */

typedef struct {
    double *mem;
    /* Position fields */
    double *phi[NFIELDS];      /* displacement */
    double *phi_vel[NFIELDS];
    double *phi_acc[NFIELDS];
    /* Angle fields */
    double *theta[NFIELDS];    /* rotation */
    double *theta_vel[NFIELDS];
    double *theta_acc[NFIELDS];
    int N; long N3;
    double L, dx, dt;
} Grid;

static Grid *grid_alloc(int N, double L) {
    Grid *g = calloc(1, sizeof(Grid));
    g->N  = N;
    g->N3 = (long)N * N * N;
    g->L  = L;
    g->dx = 2.0 * L / (N - 1);
    g->dt = 0.10 * g->dx;  /* slightly smaller dt for stability with coupling */

    /* ONE allocation: 18 arrays (9 position + 9 angle) */
    long total = 18 * g->N3;
    double bytes = total * sizeof(double);
    printf("Allocating %.2f GB (%ld doubles, N=%d, 6 fields)\n", bytes/1e9, total, N);
    g->mem = malloc(total * sizeof(double));
    if (!g->mem) { fprintf(stderr, "FATAL: malloc failed\n"); exit(1); }
    memset(g->mem, 0, total * sizeof(double));

    for (int a = 0; a < NFIELDS; a++) {
        g->phi[a]       = g->mem + (0  + a) * g->N3;
        g->phi_vel[a]   = g->mem + (3  + a) * g->N3;
        g->phi_acc[a]   = g->mem + (6  + a) * g->N3;
        g->theta[a]     = g->mem + (9  + a) * g->N3;
        g->theta_vel[a] = g->mem + (12 + a) * g->N3;
        g->theta_acc[a] = g->mem + (15 + a) * g->N3;
    }
    return g;
}

static void grid_free(Grid *g) { free(g->mem); free(g); }

/* ================================================================
   Curl computation helpers
   ================================================================ */

/* curl(F)_a at grid point idx, for field array F[3]
   curl(F)_0 = ∂F₂/∂y - ∂F₁/∂z
   curl(F)_1 = ∂F₀/∂z - ∂F₂/∂x
   curl(F)_2 = ∂F₁/∂x - ∂F₀/∂y  */
static inline double curl_component(double *F[3], int a,
    long n_ip, long n_im, long n_jp, long n_jm, long n_kp, long n_km,
    double idx1) {
    if (a == 0) return (F[2][n_jp] - F[2][n_jm] - F[1][n_kp] + F[1][n_km]) * idx1;
    if (a == 1) return (F[0][n_kp] - F[0][n_km] - F[2][n_ip] + F[2][n_im]) * idx1;
    return            (F[1][n_ip] - F[1][n_im] - F[0][n_jp] + F[0][n_jm]) * idx1;
}

/* ================================================================
   Forces: position + angle, with coupling
   ================================================================ */

static void compute_forces(Grid *g) {
    const int N = g->N, NN = N * N;
    const long N3 = g->N3;
    const double idx2 = 1.0 / (g->dx * g->dx);
    const double idx1 = 1.0 / (2.0 * g->dx);

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

        /* --- Position field forces --- */
        double p0 = g->phi[0][idx], p1 = g->phi[1][idx], p2 = g->phi[2][idx];
        double P = p0 * p1 * p2;
        double den = 1.0 + KAPPA * P * P;
        double mPd2 = MU * P / (den * den);

        for (int a = 0; a < NFIELDS; a++) {
            /* Standard: Laplacian + mass + V'(P) */
            double lap = (g->phi[a][n_ip] + g->phi[a][n_im]
                        + g->phi[a][n_jp] + g->phi[a][n_jm]
                        + g->phi[a][n_kp] + g->phi[a][n_km]
                        - 6.0 * g->phi[a][idx]) * idx2;
            double dPda = (a==0) ? p1*p2 : (a==1) ? p0*p2 : p0*p1;

            /* Coupling: + η × curl(θ)_a */
            double curl_theta = curl_component(g->theta, a,
                n_ip, n_im, n_jp, n_jm, n_kp, n_km, idx1);

            g->phi_acc[a][idx] = lap - MASS2 * g->phi[a][idx]
                               - mPd2 * dPda + ETA * curl_theta;
        }

        /* --- Angle field forces --- */
        for (int a = 0; a < NFIELDS; a++) {
            /* Laplacian of θ_a */
            double lap_t = (g->theta[a][n_ip] + g->theta[a][n_im]
                          + g->theta[a][n_jp] + g->theta[a][n_jm]
                          + g->theta[a][n_kp] + g->theta[a][n_km]
                          - 6.0 * g->theta[a][idx]) * idx2;

            /* Coupling: + η × curl(φ)_a */
            double curl_phi = curl_component(g->phi, a,
                n_ip, n_im, n_jp, n_jm, n_kp, n_km, idx1);

            g->theta_acc[a][idx] = lap_t - MTHETA2 * g->theta[a][idx]
                                 + ETA * curl_phi;
        }
    }
}

/* ================================================================
   Symplectic Velocity Verlet (all 6 fields)
   ================================================================ */

static void verlet_step(Grid *g) {
    const long N3 = g->N3;
    const double hdt = 0.5 * g->dt, dt = g->dt;

    /* Half-kick: all velocities */
    for (int a = 0; a < NFIELDS; a++) {
        double *vp = g->phi_vel[a], *ap = g->phi_acc[a];
        double *vt = g->theta_vel[a], *at = g->theta_acc[a];
        for (long i = 0; i < N3; i++) { vp[i] += hdt*ap[i]; vt[i] += hdt*at[i]; }
    }
    /* Drift: all fields */
    for (int a = 0; a < NFIELDS; a++) {
        double *pp = g->phi[a], *vp = g->phi_vel[a];
        double *pt = g->theta[a], *vt = g->theta_vel[a];
        for (long i = 0; i < N3; i++) { pp[i] += dt*vp[i]; pt[i] += dt*vt[i]; }
    }
    /* Recompute forces */
    compute_forces(g);
    /* Half-kick */
    for (int a = 0; a < NFIELDS; a++) {
        double *vp = g->phi_vel[a], *ap = g->phi_acc[a];
        double *vt = g->theta_vel[a], *at = g->theta_acc[a];
        for (long i = 0; i < N3; i++) { vp[i] += hdt*ap[i]; vt[i] += hdt*at[i]; }
    }
}

/* ================================================================
   Initialization
   ================================================================ */

static void init_braid(Grid *g, double x_cen, double y_cen) {
    const int N = g->N, NN = N*N;
    const double dx = g->dx, L = g->L;
    const double A[3] = {0.8, 0.8, 0.8};
    const double delta[3] = {0, 3.0005, 4.4325};
    const double R_tube = 3.0, ellip = 0.3325;
    const double kw = PI/L, omega = sqrt(kw*kw + MASS2);
    const double sx = 1+ellip, sy = 1-ellip;
    const double inv2R2 = 1.0/(2*R_tube*R_tube);
    const double k_bg = PI/L, omega_bg = sqrt(k_bg*k_bg + MASS2);

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
                    g->phi[a][idx] += A[a]*env*cos(ph) + A_BG*cos(ph_bg);
                    g->phi_vel[a][idx] += omega*A[a]*env*sin(ph) + omega_bg*A_BG*sin(ph_bg);
                }
                /* θ fields start at zero — will be sourced by curl(φ) */
            }
        }
    }
}

/* ================================================================
   Diagnostics
   ================================================================ */

static void compute_energy(Grid *g, double *E_phi_kin, double *E_theta_kin,
                           double *E_grad, double *E_mass, double *E_pot,
                           double *E_theta_grad, double *E_theta_mass,
                           double *E_coupling, double *E_total) {
    const int N = g->N, NN = N*N;
    const long N3 = g->N3;
    const double dx = g->dx, dV = dx*dx*dx;
    const double idx1 = 1.0/(2.0*dx);
    double epk=0, etk=0, eg=0, em=0, ep=0, etg=0, etm=0, ec=0;

    #pragma omp parallel for reduction(+:epk,etk,eg,em,ep,etg,etm,ec) schedule(static)
    for (long idx = 0; idx < N3; idx++) {
        int i=(int)(idx/NN), j=(int)((idx/N)%N), k=(int)(idx%N);
        int ip=(i+1)%N, im=(i-1+N)%N;
        int jp=(j+1)%N, jm=(j-1+N)%N;
        int kp=(k+1)%N, km=(k-1+N)%N;
        long n_ip=(long)ip*NN+j*N+k, n_im=(long)im*NN+j*N+k;
        long n_jp=(long)i*NN+jp*N+k, n_jm=(long)i*NN+jm*N+k;
        long n_kp=(long)i*NN+j*N+kp, n_km=(long)i*NN+j*N+km;

        for (int a = 0; a < NFIELDS; a++) {
            /* Position kinetic */
            epk += 0.5*g->phi_vel[a][idx]*g->phi_vel[a][idx]*dV;
            /* Angle kinetic */
            etk += 0.5*g->theta_vel[a][idx]*g->theta_vel[a][idx]*dV;
            /* Position gradient */
            double gx=(g->phi[a][n_ip]-g->phi[a][n_im])*idx1;
            double gy=(g->phi[a][n_jp]-g->phi[a][n_jm])*idx1;
            double gz=(g->phi[a][n_kp]-g->phi[a][n_km])*idx1;
            eg += 0.5*(gx*gx+gy*gy+gz*gz)*dV;
            /* Position mass */
            em += 0.5*MASS2*g->phi[a][idx]*g->phi[a][idx]*dV;
            /* Angle gradient */
            double tgx=(g->theta[a][n_ip]-g->theta[a][n_im])*idx1;
            double tgy=(g->theta[a][n_jp]-g->theta[a][n_jm])*idx1;
            double tgz=(g->theta[a][n_kp]-g->theta[a][n_km])*idx1;
            etg += 0.5*(tgx*tgx+tgy*tgy+tgz*tgz)*dV;
            /* Angle mass */
            etm += 0.5*MTHETA2*g->theta[a][idx]*g->theta[a][idx]*dV;
        }
        /* Triple product potential */
        double P = g->phi[0][idx]*g->phi[1][idx]*g->phi[2][idx];
        ep += (MU/2.0)*P*P/(1.0+KAPPA*P*P)*dV;

        /* Coupling energy: -η × φ · curl(θ) or equivalently -η × θ · curl(φ)
           E_coupling = -η ∫ Σ_a φ_a × curl(θ)_a dV  */
        for (int a = 0; a < NFIELDS; a++) {
            double ct = curl_component(g->theta, a,
                n_ip, n_im, n_jp, n_jm, n_kp, n_km, idx1);
            ec -= ETA * g->phi[a][idx] * ct * dV;
        }
    }
    *E_phi_kin=epk; *E_theta_kin=etk; *E_grad=eg; *E_mass=em;
    *E_pot=ep; *E_theta_grad=etg; *E_theta_mass=etm; *E_coupling=ec;
    *E_total = epk+etk+eg+em+ep+etg+etm+ec;
}

/* Braid centroid (position fields only) */
static double measure_braid_x(Grid *g) {
    const long N3 = g->N3;
    const int NN = g->N*g->N, N = g->N;
    const double dx = g->dx, L = g->L;
    double avg = 0;
    for (long idx = 0; idx < N3; idx++) {
        double p2 = 0;
        for (int a=0;a<NFIELDS;a++) p2 += g->phi[a][idx]*g->phi[a][idx];
        avg += p2;
    }
    avg /= N3;
    double wx = 0, wt = 0;
    for (long idx = 0; idx < N3; idx++) {
        double p2 = 0;
        for (int a=0;a<NFIELDS;a++) p2 += g->phi[a][idx]*g->phi[a][idx];
        if (p2 < 5.0*avg) continue;
        int i = (int)(idx / NN);
        double x = -L + i*dx;
        wx += x*p2; wt += p2;
    }
    return (wt > 0) ? wx/wt : 0;
}

/* θ field RMS (measure of angle field excitation) */
static double theta_rms(Grid *g) {
    double sum = 0;
    for (int a = 0; a < NFIELDS; a++)
        for (long i = 0; i < g->N3; i++)
            sum += g->theta[a][i] * g->theta[a][i];
    return sqrt(sum / (3 * g->N3));
}

/* SFA archive writer (replaces individual .bin snapshots) */
static SFA *sfa_archive = NULL;

static void sfa_write_snapshot(Grid *g, double t) {
    void *cols[] = {
        g->phi[0], g->phi[1], g->phi[2],
        g->theta[0], g->theta[1], g->theta[2]
    };
    sfa_write_frame(sfa_archive, t, cols);
}

/* ================================================================
   Main
   ================================================================ */

int main(int argc, char **argv) {
    int N = 128;
    double L = 20.0, T = 300.0;
    int n_braids = 1;
    double D = 20.0;
    double diag_dt = 5.0, snap_dt = 100.0;
    double dt_factor = 1.0;
    char outdir[256] = "data/cosserat";

    for (int i = 1; i < argc; i++) {
        if      (!strcmp(argv[i],"-N"))      N = atoi(argv[++i]);
        else if (!strcmp(argv[i],"-L"))      L = atof(argv[++i]);
        else if (!strcmp(argv[i],"-T"))      T = atof(argv[++i]);
        else if (!strcmp(argv[i],"-braids")) n_braids = atoi(argv[++i]);
        else if (!strcmp(argv[i],"-D"))      D = atof(argv[++i]);
        else if (!strcmp(argv[i],"-bg"))     A_BG = atof(argv[++i]);
        else if (!strcmp(argv[i],"-m"))      MASS2 = atof(argv[++i]) * atof(argv[i]);
        else if (!strcmp(argv[i],"-mt"))     MTHETA2 = atof(argv[++i]) * atof(argv[i]);
        else if (!strcmp(argv[i],"-eta"))    ETA = atof(argv[++i]);
        else if (!strcmp(argv[i],"-diag"))   diag_dt = atof(argv[++i]);
        else if (!strcmp(argv[i],"-snap"))   snap_dt = atof(argv[++i]);
        else if (!strcmp(argv[i],"-o"))      strncpy(outdir, argv[++i], 255);
        else if (!strcmp(argv[i],"-dt_factor")) dt_factor = atof(argv[++i]);
    }

    int nthreads = 8;
    char *env_threads = getenv("OMP_NUM_THREADS");
    if (env_threads) nthreads = atoi(env_threads);
    omp_set_num_threads(nthreads);

    printf("=== V33 Cosserat: 3 Position + 3 Angle Fields ===\n");
    printf("∂²φ_a/∂t² = ∇²φ_a - m²φ_a - V'(φ) + η×curl(θ)_a\n");
    printf("∂²θ_a/∂t² = ∇²θ_a - m_θ²θ_a       + η×curl(φ)_a\n\n");
    printf("m² = %.4f, m_θ² = %.4f, η = %.3f\n", MASS2, MTHETA2, ETA);
    printf("N=%d L=%.0f T=%.0f braids=%d D=%.0f bg=%.2f\n", N, L, T, n_braids, D, A_BG);

    mkdir("data", 0755);
    Grid *g = grid_alloc(N, L);
    g->dt *= dt_factor;
    printf("dx=%.4f dt=%.5f (factor=%.2f) threads=%d\n\n", g->dx, g->dt, dt_factor, nthreads);

    if (n_braids == 1) {
        init_braid(g, 0, 0);
    } else {
        for (int b = 0; b < n_braids; b++) {
            double bx = (n_braids==2) ? (b==0?-D/2:D/2) : D/2*cos(2*PI*b/n_braids);
            double by = (n_braids==2) ? 0 : D/2*sin(2*PI*b/n_braids);
            init_braid(g, bx, by);
        }
    }
    compute_forces(g);

    /* Set up SFA archive */
    char sfapath[512];
    snprintf(sfapath, sizeof(sfapath), "%s.sfa", outdir);
    sfa_archive = sfa_create(sfapath, N, N, N, L, L, L, g->dt);
    sfa_add_column(sfa_archive, "phi_x",   SFA_F64, SFA_POSITION, 0);
    sfa_add_column(sfa_archive, "phi_y",   SFA_F64, SFA_POSITION, 1);
    sfa_add_column(sfa_archive, "phi_z",   SFA_F64, SFA_POSITION, 2);
    sfa_add_column(sfa_archive, "theta_x", SFA_F64, SFA_ANGLE,    0);
    sfa_add_column(sfa_archive, "theta_y", SFA_F64, SFA_ANGLE,    1);
    sfa_add_column(sfa_archive, "theta_z", SFA_F64, SFA_ANGLE,    2);
    sfa_finalize_header(sfa_archive);
    printf("SFA archive: %s (6 columns, f64)\n\n", sfapath);

    /* Also keep text timeseries */
    char tspath[512];
    snprintf(tspath, sizeof(tspath), "%s_timeseries.tsv", outdir);
    FILE *fp = fopen(tspath, "w");
    fprintf(fp, "t\tE_phi_kin\tE_theta_kin\tE_grad\tE_mass\tE_pot\tE_tgrad\tE_tmass\tE_coupling\tE_total\ttheta_rms\n");

    int n_steps = (int)(T / g->dt);
    int diag_every = (int)(diag_dt / g->dt); if (diag_every<1) diag_every=1;

    printf("Steps=%d diag=%d (writing every diag frame to SFA)\n\n", n_steps, diag_every);

    double wall0 = omp_get_wtime();
    double E0 = 0;
    int frame_count = 0;

    /* Write initial frame */
    sfa_write_snapshot(g, 0);
    frame_count++;

    for (int step = 0; step <= n_steps; step++) {
        if (step > 0) verlet_step(g);
        double t = step * g->dt;

        if (step % diag_every == 0) {
            double epk, etk, eg, em, ep, etg, etm, ec, et;
            compute_energy(g, &epk, &etk, &eg, &em, &ep, &etg, &etm, &ec, &et);
            if (step == 0) E0 = et;
            double trms = theta_rms(g);

            fprintf(fp, "%.2f\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\n",
                    t, epk, etk, eg, em, ep, etg, etm, ec, et, trms);
            fflush(fp);

            /* Write every diag frame to SFA */
            if (step > 0) {
                sfa_write_snapshot(g, t);
                frame_count++;
            }

            if (step % (diag_every*10) == 0) {
                double wall = omp_get_wtime() - wall0;
                double drift = 100.0*(et-E0)/(fabs(E0)+1e-30);
                printf("t=%7.1f E=%.4e (drift %+.3f%%) θ_rms=%.3e Ep=%.1f frames=%d [%.0f%% %.0fs]\n",
                       t, et, drift, trms, ep, frame_count, 100.0*step/n_steps, wall);
                fflush(stdout);
            }
        }
    }

    fclose(fp);
    sfa_close(sfa_archive);

    double wall = omp_get_wtime() - wall0;
    printf("\n=== Complete: %.0fs (%.1f min) ===\n", wall, wall/60);
    printf("SFA: %s (%d frames)\n", sfapath, frame_count);

    /* Show file size */
    FILE *sz = fopen(sfapath, "rb");
    fseek(sz, 0, SEEK_END);
    long fsize = ftell(sz);
    fclose(sz);
    long raw = (long)frame_count * g->N3 * 6 * 8;
    printf("Size: %.1f MB (raw would be %.1f MB, ratio %.1f×)\n",
           fsize/1e6, raw/1e6, (double)raw/fsize);

    grid_free(g);
    return 0;
}
