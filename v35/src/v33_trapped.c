/*  v33_trapped.c — 6-field Cosserat with trapped theta wave packet
 *
 *  Tests whether a theta wave packet can be trapped in the phi-braid's
 *  potential well (electron analog). NO theta self-interaction.
 *  The phi-braid provides the "box" through its curl coupling.
 *
 *  Key differences from v33_theta_self.c:
 *  - MU_THETA=0 always (no theta self-interaction)
 *  - Theta initialized as a Gaussian wave packet, not a braid
 *  - Theta packet centroid tracking diagnostic
 *
 *  Build: gcc -O3 -march=native -fopenmp -o v33_trapped \
 *         src/v33_trapped.c -lm
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>
#include <sys/stat.h>
#include <errno.h>

#define NFIELDS 3
#define PI 3.14159265358979323846

/* Physics parameters */
static double MU       = -41.345;
static double KAPPA    = 50.0;
static double MASS2    = 2.25;    /* position field mass */
static double MTHETA2  = 0.0;    /* angle field mass (0 = massless theta) */
static double A_BG     = 0.1;
static double ETA      = 0.5;    /* position-angle coupling strength */

/* Theta wave packet parameters */
static double THETA_X     = 10.0;  /* packet center x */
static double THETA_Y     = 0.0;   /* packet center y */
static double A_THETA     = 0.05;  /* packet amplitude */
static double SIGMA_THETA = 2.0;   /* packet Gaussian width */

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
    g->dt = 0.10 * g->dx;

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

static inline double curl_component(double *F[3], int a,
    long n_ip, long n_im, long n_jp, long n_jm, long n_kp, long n_km,
    double idx1) {
    if (a == 0) return (F[2][n_jp] - F[2][n_jm] - F[1][n_kp] + F[1][n_km]) * idx1;
    if (a == 1) return (F[0][n_kp] - F[0][n_km] - F[2][n_ip] + F[2][n_im]) * idx1;
    return            (F[1][n_ip] - F[1][n_im] - F[0][n_jp] + F[0][n_jm]) * idx1;
}

/* ================================================================
   Forces: position + angle, with curl coupling only (no theta self-int)
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
            double lap = (g->phi[a][n_ip] + g->phi[a][n_im]
                        + g->phi[a][n_jp] + g->phi[a][n_jm]
                        + g->phi[a][n_kp] + g->phi[a][n_km]
                        - 6.0 * g->phi[a][idx]) * idx2;
            double dPda = (a==0) ? p1*p2 : (a==1) ? p0*p2 : p0*p1;

            double curl_theta = curl_component(g->theta, a,
                n_ip, n_im, n_jp, n_jm, n_kp, n_km, idx1);

            g->phi_acc[a][idx] = lap - MASS2 * g->phi[a][idx]
                               - mPd2 * dPda + ETA * curl_theta;
        }

        /* --- Angle field forces (no self-interaction) --- */
        for (int a = 0; a < NFIELDS; a++) {
            double lap_t = (g->theta[a][n_ip] + g->theta[a][n_im]
                          + g->theta[a][n_jp] + g->theta[a][n_jm]
                          + g->theta[a][n_kp] + g->theta[a][n_km]
                          - 6.0 * g->theta[a][idx]) * idx2;

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
            }
        }
    }
}

/* Initialize a theta WAVE PACKET (Gaussian envelope, not braid) */
static void init_theta_packet(Grid *g, double x_cen, double y_cen,
                               double amp, double sigma) {
    if (amp <= 0.0) return;
    const int N = g->N, NN = N*N;
    const double dx = g->dx, L = g->L;
    const double kw = PI/L;
    /* Massless dispersion: omega = k*c = k (c=1) */
    /* If theta has mass, omega = sqrt(k^2 + m_theta^2) */
    const double omega_t = sqrt(kw*kw + MTHETA2);
    const double inv2sig2 = 1.0 / (2.0 * sigma * sigma);

    printf("Initializing theta wave packet: A=%.4f sigma=%.2f at (%.1f, %.1f)\n",
           amp, sigma, x_cen, y_cen);
    printf("  k_z = pi/L = %.4f, omega = %.4f\n", kw, omega_t);

    for (int i = 0; i < N; i++) {
        double x = -L + i*dx;
        for (int j = 0; j < N; j++) {
            double y = -L + j*dx;
            double r2 = (x - x_cen)*(x - x_cen) + (y - y_cen)*(y - y_cen);
            double env = amp * exp(-r2 * inv2sig2);
            for (int kk = 0; kk < N; kk++) {
                double z = -L + kk*dx;
                long idx = (long)i*NN + j*N + kk;
                for (int a = 0; a < NFIELDS; a++) {
                    double ph = kw*z + 2*PI*a/3.0;
                    g->theta[a][idx] += env * cos(ph);
                    g->theta_vel[a][idx] += omega_t * env * sin(ph);
                }
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
                           double *E_coupling,
                           double *E_total) {
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
            epk += 0.5*g->phi_vel[a][idx]*g->phi_vel[a][idx]*dV;
            etk += 0.5*g->theta_vel[a][idx]*g->theta_vel[a][idx]*dV;
            double gx=(g->phi[a][n_ip]-g->phi[a][n_im])*idx1;
            double gy=(g->phi[a][n_jp]-g->phi[a][n_jm])*idx1;
            double gz=(g->phi[a][n_kp]-g->phi[a][n_km])*idx1;
            eg += 0.5*(gx*gx+gy*gy+gz*gz)*dV;
            em += 0.5*MASS2*g->phi[a][idx]*g->phi[a][idx]*dV;
            double tgx=(g->theta[a][n_ip]-g->theta[a][n_im])*idx1;
            double tgy=(g->theta[a][n_jp]-g->theta[a][n_jm])*idx1;
            double tgz=(g->theta[a][n_kp]-g->theta[a][n_km])*idx1;
            etg += 0.5*(tgx*tgx+tgy*tgy+tgz*tgz)*dV;
            etm += 0.5*MTHETA2*g->theta[a][idx]*g->theta[a][idx]*dV;
        }
        /* phi triple product potential */
        double P = g->phi[0][idx]*g->phi[1][idx]*g->phi[2][idx];
        ep += (MU/2.0)*P*P/(1.0+KAPPA*P*P)*dV;

        /* Coupling energy */
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

/* theta field RMS */
static double theta_rms(Grid *g) {
    double sum = 0;
    for (int a = 0; a < NFIELDS; a++)
        for (long i = 0; i < g->N3; i++)
            sum += g->theta[a][i] * g->theta[a][i];
    return sqrt(sum / (3 * g->N3));
}

/* ================================================================
   Theta packet centroid tracking
   ================================================================
   Computes theta^2 weighted centroid in the xy-plane and total
   theta energy near that centroid. Uses theta_kin + theta_grad
   density as weight for robust tracking even as packet disperses.
   ================================================================ */

static void measure_theta_packet(Grid *g,
    double *cx, double *cy, double *cz,
    double *theta_E, double *theta_peak, double *theta_spread) {

    const int N = g->N, NN = N*N;
    const double dx = g->dx, L = g->L, dV = dx*dx*dx;
    const double idx1 = 1.0/(2.0*dx);

    /* Pass 1: compute theta^2 at each xy column (sum over z),
       use as weight for xy centroid */
    double wx_sum = 0, wy_sum = 0, wz_sum = 0, w_total = 0;
    double peak = 0;

    /* Use theta^2 as the weight (simpler, avoids gradient cost in pass 1) */
    for (int i = 0; i < N; i++) {
        double x = -L + i*dx;
        for (int j = 0; j < N; j++) {
            double y = -L + j*dx;
            for (int kk = 0; kk < N; kk++) {
                double z = -L + kk*dx;
                long idx = (long)i*NN + j*N + kk;
                double t2 = 0;
                for (int a = 0; a < NFIELDS; a++)
                    t2 += g->theta[a][idx] * g->theta[a][idx];
                if (t2 > peak) peak = t2;
                double w = t2 * dV;
                wx_sum += w * x;
                wy_sum += w * y;
                wz_sum += w * z;
                w_total += w;
            }
        }
    }

    if (w_total > 0) {
        *cx = wx_sum / w_total;
        *cy = wy_sum / w_total;
        *cz = wz_sum / w_total;
    } else {
        *cx = *cy = *cz = 0;
    }
    *theta_peak = sqrt(peak);

    /* Pass 2: compute theta energy near the centroid (within R_cut of cx,cy)
       and spread (RMS distance from centroid) */
    double R_cut = 3.0 * SIGMA_THETA;  /* capture region: 3 sigma */
    double e_local = 0, spread_sum = 0, spread_w = 0;

    for (int i = 0; i < N; i++) {
        double x = -L + i*dx;
        for (int j = 0; j < N; j++) {
            double y = -L + j*dx;
            double r2 = (x - *cx)*(x - *cx) + (y - *cy)*(y - *cy);
            /* Always compute spread, cut energy at R_cut */
            for (int kk = 0; kk < N; kk++) {
                long idx = (long)i*NN + j*N + kk;
                double t2 = 0;
                for (int a = 0; a < NFIELDS; a++)
                    t2 += g->theta[a][idx] * g->theta[a][idx];
                double w = t2 * dV;
                spread_sum += w * r2;
                spread_w += w;

                if (r2 < R_cut*R_cut) {
                    /* kinetic energy */
                    for (int a = 0; a < NFIELDS; a++)
                        e_local += 0.5*g->theta_vel[a][idx]*g->theta_vel[a][idx]*dV;
                    /* gradient energy */
                    int ip=(i+1)%N, im=(i-1+N)%N;
                    int jp=(j+1)%N, jm=(j-1+N)%N;
                    int kp=(kk+1)%N, km=(kk-1+N)%N;
                    for (int a = 0; a < NFIELDS; a++) {
                        double tgx=(g->theta[a][(long)ip*NN+j*N+kk]-g->theta[a][(long)im*NN+j*N+kk])*idx1;
                        double tgy=(g->theta[a][(long)i*NN+jp*N+kk]-g->theta[a][(long)i*NN+jm*N+kk])*idx1;
                        double tgz=(g->theta[a][(long)i*NN+j*N+kp]-g->theta[a][(long)i*NN+j*N+km])*idx1;
                        e_local += 0.5*(tgx*tgx+tgy*tgy+tgz*tgz)*dV;
                    }
                    /* mass energy */
                    if (MTHETA2 > 0) {
                        for (int a = 0; a < NFIELDS; a++)
                            e_local += 0.5*MTHETA2*g->theta[a][idx]*g->theta[a][idx]*dV;
                    }
                }
            }
        }
    }
    *theta_E = e_local;
    *theta_spread = (spread_w > 0) ? sqrt(spread_sum / spread_w) : 0;
}

/* ================================================================
   Recursive mkdir (like mkdir -p)
   ================================================================ */

static void mkdirs(const char *path) {
    char tmp[512];
    strncpy(tmp, path, sizeof(tmp)-1);
    tmp[sizeof(tmp)-1] = '\0';
    for (char *p = tmp + 1; *p; p++) {
        if (*p == '/') {
            *p = '\0';
            mkdir(tmp, 0755);
            *p = '/';
        }
    }
    mkdir(tmp, 0755);
}

/* ================================================================
   Main
   ================================================================ */

int main(int argc, char **argv) {
    int N = 80;
    double L = 25.0, T = 200.0;
    int n_braids = 1;
    double D = 20.0;
    double diag_dt = 1.0;
    double dt_factor = 1.0;
    char outdir[256] = "data/trapped/default";

    for (int i = 1; i < argc; i++) {
        if      (!strcmp(argv[i],"-N"))      N = atoi(argv[++i]);
        else if (!strcmp(argv[i],"-L"))      L = atof(argv[++i]);
        else if (!strcmp(argv[i],"-T"))      T = atof(argv[++i]);
        else if (!strcmp(argv[i],"-braids")) n_braids = atoi(argv[++i]);
        else if (!strcmp(argv[i],"-D"))      D = atof(argv[++i]);
        else if (!strcmp(argv[i],"-bg"))     A_BG = atof(argv[++i]);
        else if (!strcmp(argv[i],"-m"))      { double v = atof(argv[++i]); MASS2 = v*v; }
        else if (!strcmp(argv[i],"-mt"))     { double v = atof(argv[++i]); MTHETA2 = v*v; }
        else if (!strcmp(argv[i],"-eta"))    ETA = atof(argv[++i]);
        else if (!strcmp(argv[i],"-diag"))   diag_dt = atof(argv[++i]);
        else if (!strcmp(argv[i],"-o"))      strncpy(outdir, argv[++i], 255);
        else if (!strcmp(argv[i],"-dt_factor")) dt_factor = atof(argv[++i]);
        else if (!strcmp(argv[i],"-A_theta"))    A_THETA = atof(argv[++i]);
        else if (!strcmp(argv[i],"-theta_x"))    THETA_X = atof(argv[++i]);
        else if (!strcmp(argv[i],"-theta_y"))    THETA_Y = atof(argv[++i]);
        else if (!strcmp(argv[i],"-sigma_theta")) SIGMA_THETA = atof(argv[++i]);
    }

    int nthreads = 4;
    char *env_threads = getenv("OMP_NUM_THREADS");
    if (env_threads) nthreads = atoi(env_threads);
    omp_set_num_threads(nthreads);

    printf("=== V35 Trapped Theta Wave Packet ===\n");
    printf("Equations:\n");
    printf("  d^2 phi_a/dt^2 = lap(phi_a) - m^2 phi_a - V'(phi) + eta*curl(theta)_a\n");
    printf("  d^2 theta_a/dt^2 = lap(theta_a) - m_t^2 theta_a + eta*curl(phi)_a\n\n");
    printf("  V(P) = (mu/2) P^2/(1+kap P^2),  P = phi_0 phi_1 phi_2\n");
    printf("  NO theta self-interaction (trapping by phi-braid only)\n\n");
    printf("phi:   mu=%.3f  kap=%.1f  m^2=%.4f\n", MU, KAPPA, MASS2);
    printf("theta: m_t^2=%.4f\n", MTHETA2);
    printf("eta=%.3f  bg=%.2f\n", ETA, A_BG);
    printf("theta packet: A=%.4f sigma=%.2f at (%.1f, %.1f)\n",
           A_THETA, SIGMA_THETA, THETA_X, THETA_Y);
    printf("N=%d L=%.0f T=%.0f braids=%d D=%.0f\n\n", N, L, T, n_braids, D);

    mkdirs(outdir);
    Grid *g = grid_alloc(N, L);
    g->dt *= dt_factor;
    printf("dx=%.4f dt=%.5f (factor=%.2f) threads=%d\n\n", g->dx, g->dt, dt_factor, nthreads);

    /* Initialize phi-braid(s) */
    if (n_braids == 1) {
        init_braid(g, 0, 0);
    } else {
        for (int b = 0; b < n_braids; b++) {
            double bx = (n_braids==2) ? (b==0?-D/2:D/2) : D/2*cos(2*PI*b/n_braids);
            double by = (n_braids==2) ? 0 : D/2*sin(2*PI*b/n_braids);
            init_braid(g, bx, by);
        }
    }

    /* Initialize theta wave packet */
    init_theta_packet(g, THETA_X, THETA_Y, A_THETA, SIGMA_THETA);

    compute_forces(g);

    /* Timeseries file: energy + theta packet tracking */
    char tspath[512];
    snprintf(tspath, sizeof(tspath), "%s/timeseries.tsv", outdir);
    FILE *fp = fopen(tspath, "w");
    if (!fp) { fprintf(stderr, "FATAL: cannot open %s\n", tspath); exit(1); }
    fprintf(fp, "t\tE_phi_kin\tE_theta_kin\tE_grad\tE_mass\tE_pot\tE_tgrad\tE_tmass\tE_coupling\tE_total\ttheta_rms\tpkt_x\tpkt_y\tpkt_r\tpkt_E\tpkt_peak\tpkt_spread\n");

    /* Separate theta tracking file (higher time resolution) */
    char trackpath[512];
    snprintf(trackpath, sizeof(trackpath), "%s/theta_track.tsv", outdir);
    FILE *ftrack = fopen(trackpath, "w");
    if (!ftrack) { fprintf(stderr, "FATAL: cannot open %s\n", trackpath); exit(1); }
    fprintf(ftrack, "t\tpkt_x\tpkt_y\tpkt_r\tpkt_E\tpkt_peak\tpkt_spread\n");

    int n_steps = (int)(T / g->dt);
    int diag_every = (int)(diag_dt / g->dt); if (diag_every<1) diag_every=1;

    printf("Steps=%d diag_every=%d\n\n", n_steps, diag_every);

    double wall0 = omp_get_wtime();
    double E0 = 0;

    for (int step = 0; step <= n_steps; step++) {
        if (step > 0) verlet_step(g);
        double t = step * g->dt;

        if (step % diag_every == 0) {
            double epk, etk, eg, em, ep, etg, etm, ec, et;
            compute_energy(g, &epk, &etk, &eg, &em, &ep, &etg, &etm, &ec, &et);
            if (step == 0) E0 = et;
            double trms = theta_rms(g);

            /* Theta packet tracking */
            double px, py, pz, pE, ppeak, pspread;
            measure_theta_packet(g, &px, &py, &pz, &pE, &ppeak, &pspread);
            double pr = sqrt(px*px + py*py);

            fprintf(fp, "%.2f\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%.4f\t%.4f\t%.4f\t%.4e\t%.4e\t%.4f\n",
                    t, epk, etk, eg, em, ep, etg, etm, ec, et, trms,
                    px, py, pr, pE, ppeak, pspread);
            fflush(fp);

            fprintf(ftrack, "%.2f\t%.4f\t%.4f\t%.4f\t%.4e\t%.4e\t%.4f\n",
                    t, px, py, pr, pE, ppeak, pspread);
            fflush(ftrack);

            if (step % (diag_every*20) == 0) {
                double wall = omp_get_wtime() - wall0;
                double drift = 100.0*(et-E0)/(fabs(E0)+1e-30);
                printf("t=%7.1f E=%.4e (drift %+.3f%%) theta_rms=%.3e pkt=(%.1f,%.1f) r=%.2f E_pkt=%.2e spread=%.2f [%.0f%% %.0fs]\n",
                       t, et, drift, trms, px, py, pr, pE, pspread,
                       100.0*step/n_steps, wall);
                fflush(stdout);
            }
        }
    }

    fclose(fp);
    fclose(ftrack);

    /* Write final fields as raw binary */
    {
        char path[512];
        snprintf(path, sizeof(path), "%s/theta_final.bin", outdir);
        FILE *ft = fopen(path, "wb");
        if (ft) {
            int32_t nn = N;
            fwrite(&nn, sizeof(int32_t), 1, ft);
            for (int a = 0; a < NFIELDS; a++)
                fwrite(g->theta[a], sizeof(double), g->N3, ft);
            fclose(ft);
        }
        snprintf(path, sizeof(path), "%s/phi_final.bin", outdir);
        FILE *fp2 = fopen(path, "wb");
        if (fp2) {
            int32_t nn = N;
            fwrite(&nn, sizeof(int32_t), 1, fp2);
            for (int a = 0; a < NFIELDS; a++)
                fwrite(g->phi[a], sizeof(double), g->N3, fp2);
            fclose(fp2);
        }
    }

    double wall = omp_get_wtime() - wall0;
    printf("\n=== Complete: %.0fs (%.1f min) ===\n", wall, wall/60);
    printf("Output: %s/timeseries.tsv\n", outdir);
    printf("        %s/theta_track.tsv\n", outdir);

    grid_free(g);
    return 0;
}
