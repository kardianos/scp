/*  v30_expand.c — Rotating FRW expansion: do braids form spontaneously?
 *
 *  Dense rotating field → cosmological expansion → vortex formation → braids?
 *
 *  EOM: ∂²φ/∂t² + 3H·∂φ/∂t = (1/a²)∇²φ - m²φ - ∂V/∂φ
 *
 *  Build: gcc -O3 -fopenmp -o v30_expand src/v30_expand.c -lm
 *  Run:   ./v30_expand -nw 2               (n_wind=2)
 *         ./v30_expand -nw 0               (control, no rotation)
 *         ./v30_expand -nw 1 -H0 0.03      (slower inflation)
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>
#include <sys/stat.h>

#define NFIELDS 3
#define PI 3.14159265358979323846

/* ================================================================
   Parameters
   ================================================================ */

static int    SIM_N     = 256;
static double SIM_L     = 100.0;
static double A0_INIT   = 3.0;    /* initial field amplitude (dense) */
static int    N_WIND    = 2;      /* azimuthal phase windings */
static double H0_INF    = 0.05;   /* Hubble rate during inflation */
static double T_INF     = 50.0;   /* inflation duration */
static double T_FORM    = 300.0;  /* formation phase duration */
static double MU        = -41.345;
static double KAPPA     = 50.0;
static double MASS      = 1.50;
static double DIAG_DT   = 5.0;    /* diagnostic interval */
static double SLICE_DT  = 10.0;   /* 2D slice interval */
static char   OUTDIR[256] = "data/nw2";

/* ================================================================
   Grid
   ================================================================ */

typedef struct {
    double *phi[NFIELDS];
    double *vel[NFIELDS];
    double *acc[NFIELDS];
    int N;
    double L, dx, dt;
} Grid;

static Grid *grid_alloc(int N, double L) {
    Grid *g = calloc(1, sizeof(Grid));
    int N3 = N * N * N;
    for (int a = 0; a < NFIELDS; a++) {
        g->phi[a] = calloc(N3, sizeof(double));
        g->vel[a] = calloc(N3, sizeof(double));
        g->acc[a] = calloc(N3, sizeof(double));
    }
    g->N = N; g->L = L;
    g->dx = 2.0 * L / (N - 1);
    g->dt = 0.15;  /* fixed, safe for dx~0.78 */
    if (g->dt > 0.2 * g->dx) g->dt = 0.2 * g->dx;
    return g;
}

static void grid_free(Grid *g) {
    for (int a = 0; a < NFIELDS; a++) {
        free(g->phi[a]); free(g->vel[a]); free(g->acc[a]);
    }
    free(g);
}

/* ================================================================
   Force computation with scale factor
   ================================================================ */

static void compute_forces(Grid *g, double inv_a2) {
    int N = g->N, NN = N*N, N3 = N*N*N;
    double idx2 = inv_a2 / (g->dx * g->dx);  /* 1/a² in Laplacian */
    double mass2 = MASS * MASS;

    #pragma omp parallel for schedule(static)
    for (int idx = 0; idx < N3; idx++) {
        int i = idx / NN, j = (idx / N) % N, k = idx % N;

        /* Fully periodic BC in all directions */
        int ip = (i+1) % N, im = (i-1+N) % N;
        int jp = (j+1) % N, jm = (j-1+N) % N;
        int kp = (k+1) % N, km = (k-1+N) % N;

        int idx_ip = ip*NN + j*N + k;
        int idx_im = im*NN + j*N + k;
        int idx_jp = i*NN + jp*N + k;
        int idx_jm = i*NN + jm*N + k;
        int idx_kp = i*NN + j*N + kp;
        int idx_km = i*NN + j*N + km;

        double p0 = g->phi[0][idx], p1 = g->phi[1][idx], p2 = g->phi[2][idx];
        double P = p0 * p1 * p2;
        double denom = 1.0 + KAPPA * P * P;
        double mu_P_d2 = MU * P / (denom * denom);

        for (int a = 0; a < NFIELDS; a++) {
            double lap = (g->phi[a][idx_ip] + g->phi[a][idx_im]
                        + g->phi[a][idx_jp] + g->phi[a][idx_jm]
                        + g->phi[a][idx_kp] + g->phi[a][idx_km]
                        - 6.0 * g->phi[a][idx]) * idx2;

            double dPda = (a==0) ? p1*p2 : (a==1) ? p0*p2 : p0*p1;
            double f_triple = mu_P_d2 * dPda;

            g->acc[a][idx] = lap - mass2 * g->phi[a][idx] - f_triple;
        }
    }
}

/* ================================================================
   Verlet step with Hubble friction
   ================================================================ */

static void verlet_step_frw(Grid *g, double a, double H) {
    int N3 = g->N * g->N * g->N;
    double dt = g->dt;
    double hdt = 0.5 * dt;
    double friction = 1.0 - 3.0 * H * hdt;  /* 1st order friction per half-step */
    if (friction < 0.5) friction = 0.5;       /* safety clamp */
    double inv_a2 = 1.0 / (a * a);

    /* Half kick + friction */
    for (int a_ = 0; a_ < NFIELDS; a_++) {
        double *v = g->vel[a_], *ac = g->acc[a_];
        for (int idx = 0; idx < N3; idx++) {
            v[idx] = v[idx] * friction + hdt * ac[idx];
        }
    }

    /* Drift */
    for (int a_ = 0; a_ < NFIELDS; a_++) {
        double *p = g->phi[a_], *v = g->vel[a_];
        for (int idx = 0; idx < N3; idx++)
            p[idx] += dt * v[idx];
    }

    /* Recompute forces */
    compute_forces(g, inv_a2);

    /* Half kick + friction */
    for (int a_ = 0; a_ < NFIELDS; a_++) {
        double *v = g->vel[a_], *ac = g->acc[a_];
        for (int idx = 0; idx < N3; idx++) {
            v[idx] = v[idx] * friction + hdt * ac[idx];
        }
    }
}

/* ================================================================
   Initialization: dense rotating field
   ================================================================ */

static void init_rotating_field(Grid *g) {
    int N = g->N, NN = N*N;
    double dx = g->dx, L = g->L;
    double k0 = PI / L;
    double omega0 = sqrt(k0*k0 + MASS*MASS);
    double Omega_rot = N_WIND * 0.1;  /* angular velocity of rigid-body rotation */

    double R_blob = 0.5 * L;  /* blob radius: field decays to ~0 at 0.5L */

    printf("  Init: A0=%.1f, n_wind=%d, k0=%.4f, omega0=%.4f\n",
           A0_INIT, N_WIND, k0, omega0);
    printf("  Method: Gaussian blob with atan2 phase winding (n_wind vortex lines)\n");
    printf("  R_blob=%.1f (field decays to ~0 by r=%.1f, well inside L=%.0f)\n",
           R_blob, 1.5*R_blob, L);

    /* Initialize: Gaussian blob × helical wave × azimuthal phase winding.
       The blob ensures field → 0 at large r, avoiding boundary singularity.
       The atan2 creates genuine vortex lines in the interior.
       n_wind=0: no vortices (smooth, control).
       n_wind≥1: n_wind vortex lines passing through the blob center. */

    for (int i = 0; i < N; i++) {
        double x = -L + i * dx;
        for (int j = 0; j < N; j++) {
            double y = -L + j * dx;
            double r_perp = sqrt(x*x + y*y);
            double envelope = exp(-r_perp*r_perp / (2.0*R_blob*R_blob));

            /* Azimuthal angle for phase winding */
            double theta = atan2(y, x);

            for (int k = 0; k < N; k++) {
                double z = -L + k * dx;
                int idx = i * NN + j * N + k;

                for (int a = 0; a < NFIELDS; a++) {
                    double phase = k0 * z + 2.0*PI*a/3.0;
                    if (N_WIND > 0 && r_perp > 0.5*dx) {
                        phase += N_WIND * theta;
                    }
                    g->phi[a][idx] = A0_INIT * envelope * cos(phase);
                    g->vel[a][idx] = omega0 * A0_INIT * envelope * sin(phase);
                }
            }
        }
    }
}

/* ================================================================
   Diagnostics
   ================================================================ */

static void compute_global(Grid *g, double inv_a2, double *E_total,
                           double *E_kin, double *E_grad, double *E_mass,
                           double *E_pot, double *max_rho, double *avg_rho) {
    int N = g->N, NN = N*N, N3 = N*N*N;
    double dx = g->dx, dV = dx*dx*dx;
    double mass2 = MASS*MASS;
    double ek=0, eg=0, em=0, ep=0, mr=0, sr=0;

    #pragma omp parallel for reduction(+:ek,eg,em,ep,sr) reduction(max:mr) schedule(static)
    for (int idx = 0; idx < N3; idx++) {
        int i = idx/NN, j = (idx/N)%N, k = idx%N;
        int kp=(k+1)%N, km=(k-1+N)%N;
        int ip=(i+1)%N, im=(i-1+N)%N;
        int jp=(j+1)%N, jm=(j-1+N)%N;

        double local_ek = 0, local_eg = 0, local_em = 0;
        for (int a = 0; a < NFIELDS; a++) {
            local_ek += 0.5 * g->vel[a][idx] * g->vel[a][idx];
            double gx = (g->phi[a][ip*NN+j*N+k] - g->phi[a][im*NN+j*N+k]) / (2*dx);
            double gy = (g->phi[a][i*NN+jp*N+k] - g->phi[a][i*NN+jm*N+k]) / (2*dx);
            double gz = (g->phi[a][i*NN+j*N+kp] - g->phi[a][i*NN+j*N+km]) / (2*dx);
            local_eg += 0.5 * inv_a2 * (gx*gx + gy*gy + gz*gz);
            local_em += 0.5 * mass2 * g->phi[a][idx] * g->phi[a][idx];
        }
        double P = g->phi[0][idx]*g->phi[1][idx]*g->phi[2][idx];
        double local_ep = (MU/2.0)*P*P/(1.0+KAPPA*P*P);
        double rho = local_ek + local_eg + local_em + local_ep;

        ek += local_ek * dV;
        eg += local_eg * dV;
        em += local_em * dV;
        ep += local_ep * dV;
        sr += rho;
        if (rho > mr) mr = rho;
    }

    *E_kin = ek; *E_grad = eg; *E_mass = em; *E_pot = ep;
    *E_total = ek + eg + em + ep;
    *max_rho = mr;
    *avg_rho = sr / N3;
}

/* Particle detection: find localized energy peaks */
static int find_particles(Grid *g, double inv_a2, FILE *fp, double t) {
    int N = g->N, NN = N*N, N3 = N*N*N;
    double dx = g->dx, L = g->L;
    double mass2 = MASS*MASS;

    /* Compute energy density */
    double *rho = calloc(N3, sizeof(double));

    #pragma omp parallel for schedule(static)
    for (int idx = 0; idx < N3; idx++) {
        int i=idx/NN, j=(idx/N)%N, k=idx%N;
        int kp=(k+1)%N, km=(k-1+N)%N;
        int ip=(i+1)%N, im=(i-1+N)%N;
        int jp=(j+1)%N, jm=(j-1+N)%N;
        double e = 0;
        for (int a = 0; a < NFIELDS; a++) {
            e += 0.5 * g->vel[a][idx] * g->vel[a][idx];
            double gx = (g->phi[a][ip*NN+j*N+k]-g->phi[a][im*NN+j*N+k])/(2*dx);
            double gy = (g->phi[a][i*NN+jp*N+k]-g->phi[a][i*NN+jm*N+k])/(2*dx);
            double gz = (g->phi[a][i*NN+j*N+kp]-g->phi[a][i*NN+j*N+km])/(2*dx);
            e += 0.5*inv_a2*(gx*gx+gy*gy+gz*gz);
            e += 0.5*mass2*g->phi[a][idx]*g->phi[a][idx];
        }
        double P = g->phi[0][idx]*g->phi[1][idx]*g->phi[2][idx];
        e += (MU/2.0)*P*P/(1.0+KAPPA*P*P);
        rho[idx] = e;
    }

    /* Find threshold: 5× average */
    double sum = 0;
    for (int idx = 0; idx < N3; idx++) sum += rho[idx];
    double avg = sum / N3;
    double thresh = 5.0 * avg;

    /* Find local maxima above threshold */
    int n_peaks = 0;
    int max_peaks = 1000;
    double *peak_x = calloc(max_peaks, sizeof(double));
    double *peak_y = calloc(max_peaks, sizeof(double));
    double *peak_z = calloc(max_peaks, sizeof(double));
    double *peak_e = calloc(max_peaks, sizeof(double));

    for (int i = 2; i < N-2; i++) {
        double x = -L + i*dx;
        for (int j = 2; j < N-2; j++) {
            double y = -L + j*dx;
            for (int k = 0; k < N; k++) {
                int idx = i*NN + j*N + k;
                if (rho[idx] < thresh) continue;

                /* Check if local max (6 neighbors) */
                int kp=(k+1)%N, km=(k-1+N)%N;
                if (rho[idx] > rho[(i+1)*NN+j*N+k] &&
                    rho[idx] > rho[(i-1)*NN+j*N+k] &&
                    rho[idx] > rho[i*NN+(j+1)*N+k] &&
                    rho[idx] > rho[i*NN+(j-1)*N+k] &&
                    rho[idx] > rho[i*NN+j*N+kp] &&
                    rho[idx] > rho[i*NN+j*N+km]) {

                    if (n_peaks < max_peaks) {
                        double z = -L + k*dx;
                        peak_x[n_peaks] = x;
                        peak_y[n_peaks] = y;
                        peak_z[n_peaks] = z;
                        peak_e[n_peaks] = rho[idx];
                        n_peaks++;
                    }
                }
            }
        }
    }

    /* Cluster nearby peaks (within 5 dx) */
    double cluster_r = 5.0 * dx;
    int *cluster_id = calloc(n_peaks, sizeof(int));
    for (int p = 0; p < n_peaks; p++) cluster_id[p] = p;
    for (int p = 0; p < n_peaks; p++) {
        for (int q = p+1; q < n_peaks; q++) {
            double dr = sqrt((peak_x[p]-peak_x[q])*(peak_x[p]-peak_x[q])
                           + (peak_y[p]-peak_y[q])*(peak_y[p]-peak_y[q])
                           + (peak_z[p]-peak_z[q])*(peak_z[p]-peak_z[q]));
            if (dr < cluster_r) {
                int old = cluster_id[q], new_ = cluster_id[p];
                for (int r = 0; r < n_peaks; r++)
                    if (cluster_id[r] == old) cluster_id[r] = new_;
            }
        }
    }

    /* Count unique clusters */
    int n_particles = 0;
    int *seen = calloc(n_peaks > 0 ? n_peaks : 1, sizeof(int));
    for (int p = 0; p < n_peaks; p++) {
        int cid = cluster_id[p];
        int found = 0;
        for (int s = 0; s < n_particles; s++)
            if (seen[s] == cid) { found = 1; break; }
        if (!found && n_particles < max_peaks)
            seen[n_particles++] = cid;
    }

    /* Write particle info */
    if (fp) {
        for (int c = 0; c < n_particles; c++) {
            /* Find centroid and total energy of cluster */
            double cx=0,cy=0,cz=0,ce=0;
            int cn=0;
            for (int p = 0; p < n_peaks; p++) {
                if (cluster_id[p] == seen[c]) {
                    cx += peak_x[p]*peak_e[p];
                    cy += peak_y[p]*peak_e[p];
                    cz += peak_z[p]*peak_e[p];
                    ce += peak_e[p];
                    cn++;
                }
            }
            if (ce > 0) { cx/=ce; cy/=ce; cz/=ce; }
            fprintf(fp, "%.2f\t%d\t%.2f\t%.2f\t%.2f\t%.4e\t%d\n",
                    t, c, cx, cy, cz, ce, cn);
        }
    }

    free(rho); free(peak_x); free(peak_y); free(peak_z); free(peak_e);
    free(cluster_id); free(seen);
    return n_particles;
}

/* Write 2D slice of energy density (subsampled 2:1) */
static void write_slice_xy(Grid *g, double inv_a2, const char *fname) {
    int N = g->N, NN = N*N;
    double dx = g->dx, L = g->L;
    double mass2 = MASS*MASS;
    int k_mid = N/2;
    int sub = 2;  /* subsample factor */

    FILE *fp = fopen(fname, "w");
    fprintf(fp, "x\ty\trho\tP\n");

    for (int i = 0; i < N; i += sub) {
        double x = -L + i*dx;
        for (int j = 0; j < N; j += sub) {
            double y = -L + j*dx;
            int idx = i*NN + j*N + k_mid;
            int kp=(k_mid+1)%N, km=(k_mid-1+N)%N;
            int ip=(i+1)%N, im=(i-1+N)%N;
            int jp=(j+1)%N, jm=(j-1+N)%N;

            double e = 0;
            for (int a = 0; a < NFIELDS; a++) {
                e += 0.5*g->vel[a][idx]*g->vel[a][idx];
                double gx = (g->phi[a][ip*NN+j*N+k_mid]-g->phi[a][im*NN+j*N+k_mid])/(2*dx);
                double gy = (g->phi[a][i*NN+jp*N+k_mid]-g->phi[a][i*NN+jm*N+k_mid])/(2*dx);
                double gz = (g->phi[a][i*NN+j*N+kp]-g->phi[a][i*NN+j*N+km])/(2*dx);
                e += 0.5*inv_a2*(gx*gx+gy*gy+gz*gz);
                e += 0.5*mass2*g->phi[a][idx]*g->phi[a][idx];
            }
            double P = g->phi[0][idx]*g->phi[1][idx]*g->phi[2][idx];
            e += (MU/2.0)*P*P/(1.0+KAPPA*P*P);

            fprintf(fp, "%.2f\t%.2f\t%.6e\t%.6e\n", x, y, e, P);
        }
    }
    fclose(fp);
}

/* Write full 3D energy density as binary */
static void write_full3d(Grid *g, double inv_a2, const char *fname) {
    int N = g->N, NN = N*N, N3 = N*N*N;
    double dx = g->dx, L = g->L;
    double mass2 = MASS*MASS;

    double *rho = malloc(N3 * sizeof(double));

    #pragma omp parallel for schedule(static)
    for (int idx = 0; idx < N3; idx++) {
        int i=idx/NN, j=(idx/N)%N, k=idx%N;
        int kp=(k+1)%N, km=(k-1+N)%N;
        int ip=(i+1)%N, im=(i-1+N)%N;
        int jp=(j+1)%N, jm=(j-1+N)%N;
        double e = 0;
        for (int a = 0; a < NFIELDS; a++) {
            e += 0.5*g->vel[a][idx]*g->vel[a][idx];
            double gx=(g->phi[a][ip*NN+j*N+k]-g->phi[a][im*NN+j*N+k])/(2*dx);
            double gy=(g->phi[a][i*NN+jp*N+k]-g->phi[a][i*NN+jm*N+k])/(2*dx);
            double gz=(g->phi[a][i*NN+j*N+kp]-g->phi[a][i*NN+j*N+km])/(2*dx);
            e += 0.5*inv_a2*(gx*gx+gy*gy+gz*gz);
            e += 0.5*mass2*g->phi[a][idx]*g->phi[a][idx];
        }
        double P = g->phi[0][idx]*g->phi[1][idx]*g->phi[2][idx];
        e += (MU/2.0)*P*P/(1.0+KAPPA*P*P);
        rho[idx] = e;
    }

    FILE *fp = fopen(fname, "wb");
    int n = N;
    double l = L, t_dummy = 0;
    fwrite(&n, sizeof(int), 1, fp);
    fwrite(&l, sizeof(double), 1, fp);
    fwrite(&t_dummy, sizeof(double), 1, fp);
    fwrite(rho, sizeof(double), N3, fp);
    fclose(fp);
    free(rho);
}

/* ================================================================
   Main
   ================================================================ */

int main(int argc, char **argv) {
    /* Parse args */
    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "-nw") && i+1<argc) N_WIND = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-N") && i+1<argc) SIM_N = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-L") && i+1<argc) SIM_L = atof(argv[++i]);
        else if (!strcmp(argv[i], "-A0") && i+1<argc) A0_INIT = atof(argv[++i]);
        else if (!strcmp(argv[i], "-H0") && i+1<argc) H0_INF = atof(argv[++i]);
        else if (!strcmp(argv[i], "-tinf") && i+1<argc) T_INF = atof(argv[++i]);
        else if (!strcmp(argv[i], "-tform") && i+1<argc) T_FORM = atof(argv[++i]);
        else if (!strcmp(argv[i], "-o") && i+1<argc) strncpy(OUTDIR, argv[++i], 255);
    }

    /* Auto output dir */
    if (OUTDIR[0] == 'd' && OUTDIR[1] == 'a') { /* default */
        snprintf(OUTDIR, sizeof(OUTDIR), "data/nw%d", N_WIND);
    }

    omp_set_num_threads(16);

    printf("=== V30: Rotating FRW Expansion ===\n");
    printf("N=%d, L=%.0f, A0=%.1f, n_wind=%d\n", SIM_N, SIM_L, A0_INIT, N_WIND);
    printf("H0=%.4f, t_inf=%.0f, t_form=%.0f\n", H0_INF, T_INF, T_FORM);
    printf("a_final = exp(%.2f) = %.2f\n", H0_INF*T_INF, exp(H0_INF*T_INF));
    printf("Output: %s/\n\n", OUTDIR);

    /* Create output directory */
    mkdir("data", 0755);
    mkdir(OUTDIR, 0755);

    Grid *g = grid_alloc(SIM_N, SIM_L);
    printf("dx=%.4f, dt=%.4f\n", g->dx, g->dt);

    /* Initialize */
    init_rotating_field(g);
    double a_scale = 1.0;
    compute_forces(g, 1.0 / (a_scale * a_scale));

    /* Open timeseries file */
    char tspath[512];
    snprintf(tspath, sizeof(tspath), "%s/timeseries.tsv", OUTDIR);
    FILE *fp_ts = fopen(tspath, "w");
    fprintf(fp_ts, "t\ta\tH\tE_total\tE_kin\tE_grad\tE_mass\tE_pot\tmax_rho\tavg_rho\tn_particles\n");

    /* Open particles file */
    char partpath[512];
    snprintf(partpath, sizeof(partpath), "%s/particles.tsv", OUTDIR);
    FILE *fp_part = fopen(partpath, "w");
    fprintf(fp_part, "t\tpid\tx\ty\tz\tpeak_rho\tn_subpeaks\n");

    double T_total = T_INF + T_FORM;
    int n_steps = (int)(T_total / g->dt);
    int steps_inf = (int)(T_INF / g->dt);
    int diag_steps = (int)(DIAG_DT / g->dt);
    int slice_steps = (int)(SLICE_DT / g->dt);
    if (diag_steps < 1) diag_steps = 1;
    if (slice_steps < 1) slice_steps = 1;

    /* Full 3D dump times */
    double dump_times[] = {0, T_INF, T_INF+50, T_INF+150, T_total};
    int n_dumps = 5;
    int next_dump = 0;

    printf("Total steps: %d (inflation: %d, formation: %d)\n\n",
           n_steps, steps_inf, n_steps - steps_inf);

    double wall_start = omp_get_wtime();

    for (int step = 0; step <= n_steps; step++) {
        double t = step * g->dt;

        /* Update scale factor and Hubble parameter */
        double H;
        if (step < steps_inf) {
            H = H0_INF;
            a_scale = exp(H0_INF * t);
        } else {
            H = 0;
            a_scale = exp(H0_INF * T_INF);  /* frozen */
        }
        double inv_a2 = 1.0 / (a_scale * a_scale);

        /* Diagnostics */
        if (step % diag_steps == 0) {
            double E_tot, E_kin, E_grad, E_mass, E_pot, max_rho, avg_rho;
            compute_global(g, inv_a2, &E_tot, &E_kin, &E_grad, &E_mass,
                          &E_pot, &max_rho, &avg_rho);

            /* Particle detection (only after inflation) */
            int n_part = 0;
            if (t >= T_INF) {
                n_part = find_particles(g, inv_a2, fp_part, t);
            }

            fprintf(fp_ts, "%.2f\t%.4f\t%.5f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.6e\t%.6e\t%d\n",
                    t, a_scale, H, E_tot, E_kin, E_grad, E_mass, E_pot,
                    max_rho, avg_rho, n_part);
            fflush(fp_ts);

            double wall = omp_get_wtime() - wall_start;
            double pct = 100.0 * step / n_steps;
            printf("t=%6.1f a=%6.2f E=%.1e max_rho=%.2e n_part=%d [%.0f%% %.0fs]\n",
                   t, a_scale, E_tot, max_rho, n_part, pct, wall);
            fflush(stdout);
        }

        /* 2D slices */
        if (step % slice_steps == 0) {
            char fname[512];
            snprintf(fname, sizeof(fname), "%s/slice_xy_t%03d.tsv",
                     OUTDIR, (int)(t+0.5));
            write_slice_xy(g, inv_a2, fname);
        }

        /* Full 3D dumps */
        if (next_dump < n_dumps && t >= dump_times[next_dump] - g->dt*0.5) {
            char fname[512];
            snprintf(fname, sizeof(fname), "%s/full3d_t%03d.bin",
                     OUTDIR, (int)(t+0.5));
            printf("  Writing full 3D dump: %s (%.0f MB)\n", fname,
                   SIM_N*SIM_N*SIM_N*8.0/1e6);
            write_full3d(g, inv_a2, fname);
            next_dump++;
        }

        /* Time step */
        if (step < n_steps) {
            verlet_step_frw(g, a_scale, H);
        }
    }

    fclose(fp_ts);
    fclose(fp_part);

    double wall_total = omp_get_wtime() - wall_start;
    printf("\n=== Complete: %.0f seconds ===\n", wall_total);

    grid_free(g);
    return 0;
}
