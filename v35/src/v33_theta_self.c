/*  v33_theta_self.c — 6-field Cosserat with θ self-interaction
 *
 *  Adds triple-product self-coupling to the θ field:
 *    V_θ(Q) = (μ_θ/2) Q² / (1 + κ_θ Q²),  Q = θ₀ θ₁ θ₂
 *
 *  This gives θ its own binding potential, enabling θ-solitons.
 *
 *  Build: gcc -O3 -march=native -fopenmp -o v33_theta_self \
 *         src/v33_theta_self.c -lm
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
static double MTHETA2  = 0.0;    /* angle field mass (0 for scan) */
static double A_BG     = 0.1;
static double ETA      = 0.5;    /* position-angle coupling strength */

/* θ self-interaction parameters */
static double MU_THETA    = 0.0;   /* θ triple-product coupling (0 = off) */
static double KAPPA_THETA = 50.0;  /* θ saturation parameter */

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
   Forces: position + angle, with coupling + θ self-interaction
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

        /* --- Angle field forces --- */
        /* Pre-compute θ self-interaction terms */
        double t0 = g->theta[0][idx], t1 = g->theta[1][idx], t2 = g->theta[2][idx];
        double Q = t0 * t1 * t2;
        double den_t = 1.0 + KAPPA_THETA * Q * Q;
        double mQd2 = MU_THETA * Q / (den_t * den_t);

        for (int a = 0; a < NFIELDS; a++) {
            double lap_t = (g->theta[a][n_ip] + g->theta[a][n_im]
                          + g->theta[a][n_jp] + g->theta[a][n_jm]
                          + g->theta[a][n_kp] + g->theta[a][n_km]
                          - 6.0 * g->theta[a][idx]) * idx2;

            double curl_phi = curl_component(g->phi, a,
                n_ip, n_im, n_jp, n_jm, n_kp, n_km, idx1);

            g->theta_acc[a][idx] = lap_t - MTHETA2 * g->theta[a][idx]
                                 + ETA * curl_phi;

            /* θ self-interaction: -∂V_θ/∂θ_a */
            if (MU_THETA != 0.0) {
                double dQda = (a==0) ? t1*t2 : (a==1) ? t0*t2 : t0*t1;
                g->theta_acc[a][idx] += -mQd2 * dQda;
            }
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

/* Initialize a θ-braid: same helical pattern as φ-braid but in θ fields */
static double A_THETA_INIT = 0.0;  /* θ-braid amplitude (0 = off) */
static double THETA_BRAID_XCEN = 10.0;  /* θ-braid x center */

static void init_theta_braid(Grid *g, double x_cen, double y_cen) {
    if (A_THETA_INIT <= 0.0) return;
    const int N = g->N, NN = N*N;
    const double dx = g->dx, L = g->L;
    const double A_t[3] = {A_THETA_INIT, A_THETA_INIT, A_THETA_INIT};
    const double delta[3] = {0, 3.0005, 4.4325};
    const double R_tube = 3.0, ellip = 0.3325;
    /* Use theta mass for dispersion if nonzero, else use phi mass */
    double m2eff = (MTHETA2 > 0) ? MTHETA2 : MASS2;
    const double kw = PI/L, omega = sqrt(kw*kw + m2eff);
    const double sx = 1+ellip, sy = 1-ellip;
    const double inv2R2 = 1.0/(2*R_tube*R_tube);

    printf("Initializing theta-braid: A_theta=%.4f at (%.1f, %.1f)\n",
           A_THETA_INIT, x_cen, y_cen);

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
                    g->theta[a][idx] += A_t[a]*env*cos(ph);
                    g->theta_vel[a][idx] += omega*A_t[a]*env*sin(ph);
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
                           double *E_coupling, double *E_pot_theta,
                           double *E_total) {
    const int N = g->N, NN = N*N;
    const long N3 = g->N3;
    const double dx = g->dx, dV = dx*dx*dx;
    const double idx1 = 1.0/(2.0*dx);
    double epk=0, etk=0, eg=0, em=0, ep=0, etg=0, etm=0, ec=0, ept=0;

    #pragma omp parallel for reduction(+:epk,etk,eg,em,ep,etg,etm,ec,ept) schedule(static)
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
        /* φ triple product potential */
        double P = g->phi[0][idx]*g->phi[1][idx]*g->phi[2][idx];
        ep += (MU/2.0)*P*P/(1.0+KAPPA*P*P)*dV;

        /* θ triple product potential */
        if (MU_THETA != 0.0) {
            double Q = g->theta[0][idx]*g->theta[1][idx]*g->theta[2][idx];
            ept += (MU_THETA/2.0)*Q*Q/(1.0+KAPPA_THETA*Q*Q)*dV;
        }

        /* Coupling energy */
        for (int a = 0; a < NFIELDS; a++) {
            double ct = curl_component(g->theta, a,
                n_ip, n_im, n_jp, n_jm, n_kp, n_km, idx1);
            ec -= ETA * g->phi[a][idx] * ct * dV;
        }
    }
    *E_phi_kin=epk; *E_theta_kin=etk; *E_grad=eg; *E_mass=em;
    *E_pot=ep; *E_theta_grad=etg; *E_theta_mass=etm; *E_coupling=ec;
    *E_pot_theta=ept;
    *E_total = epk+etk+eg+em+ep+etg+etm+ec+ept;
}

/* θ field RMS */
static double theta_rms(Grid *g) {
    double sum = 0;
    for (int a = 0; a < NFIELDS; a++)
        for (long i = 0; i < g->N3; i++)
            sum += g->theta[a][i] * g->theta[a][i];
    return sqrt(sum / (3 * g->N3));
}

/* max|Q| = max|θ₀ θ₁ θ₂| */
static double theta_Q_max(Grid *g) {
    double qmax = 0;
    for (long i = 0; i < g->N3; i++) {
        double q = fabs(g->theta[0][i] * g->theta[1][i] * g->theta[2][i]);
        if (q > qmax) qmax = q;
    }
    return qmax;
}

/* θ energy within cylindrical radius R_cut of a given center (x_c, y_c) */
static double theta_energy_local(Grid *g, double x_c, double y_c, double R_cut) {
    const int N = g->N, NN = N*N;
    const double dx = g->dx, L = g->L, dV = dx*dx*dx;
    double e_kin = 0, e_pot = 0;
    for (int i = 0; i < N; i++) {
        double x = -L + i*dx;
        for (int j = 0; j < N; j++) {
            double y = -L + j*dx;
            double r2 = (x-x_c)*(x-x_c) + (y-y_c)*(y-y_c);
            if (r2 > R_cut*R_cut) continue;
            for (int kk = 0; kk < N; kk++) {
                long idx = (long)i*NN + j*N + kk;
                for (int a = 0; a < NFIELDS; a++)
                    e_kin += 0.5*g->theta_vel[a][idx]*g->theta_vel[a][idx]*dV;
                if (MU_THETA != 0.0) {
                    double Q = g->theta[0][idx]*g->theta[1][idx]*g->theta[2][idx];
                    e_pot += (MU_THETA/2.0)*Q*Q/(1.0+KAPPA_THETA*Q*Q)*dV;
                }
            }
        }
    }
    return e_kin + e_pot;
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
    int N = 64;
    double L = 20.0, T = 100.0;
    int n_braids = 1;
    double D = 20.0;
    double diag_dt = 2.0;
    double dt_factor = 1.0;
    char outdir[256] = "data/theta_self/default";

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
        else if (!strcmp(argv[i],"-mu_t"))   MU_THETA = atof(argv[++i]);
        else if (!strcmp(argv[i],"-kap_t"))  KAPPA_THETA = atof(argv[++i]);
        else if (!strcmp(argv[i],"-diag"))   diag_dt = atof(argv[++i]);
        else if (!strcmp(argv[i],"-o"))      strncpy(outdir, argv[++i], 255);
        else if (!strcmp(argv[i],"-dt_factor")) dt_factor = atof(argv[++i]);
        else if (!strcmp(argv[i],"-A_theta"))  A_THETA_INIT = atof(argv[++i]);
        else if (!strcmp(argv[i],"-theta_x"))  THETA_BRAID_XCEN = atof(argv[++i]);
    }

    int nthreads = 4;
    char *env_threads = getenv("OMP_NUM_THREADS");
    if (env_threads) nthreads = atoi(env_threads);
    omp_set_num_threads(nthreads);

    printf("=== V35 Cosserat + θ Self-Interaction ===\n");
    printf("Equations:\n");
    printf("  d²phi_a/dt² = lap(phi_a) - m²phi_a - V'(phi) + eta*curl(theta)_a\n");
    printf("  d²theta_a/dt² = lap(theta_a) - m_t²theta_a + eta*curl(phi)_a - V_t'(theta)\n\n");
    printf("  V(P) = (mu/2) P²/(1+kap P²),  P = phi_0 phi_1 phi_2\n");
    printf("  V_t(Q) = (mu_t/2) Q²/(1+kap_t Q²),  Q = theta_0 theta_1 theta_2\n\n");
    printf("phi:   mu=%.3f  kap=%.1f  m²=%.4f\n", MU, KAPPA, MASS2);
    printf("theta: mu_t=%.1f  kap_t=%.1f  m_t²=%.4f\n", MU_THETA, KAPPA_THETA, MTHETA2);
    printf("eta=%.3f  bg=%.2f\n", ETA, A_BG);
    printf("N=%d L=%.0f T=%.0f braids=%d D=%.0f\n\n", N, L, T, n_braids, D);

    mkdirs(outdir);
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
    /* Optional theta-braid initialization */
    init_theta_braid(g, THETA_BRAID_XCEN, 0);
    compute_forces(g);

    /* Timeseries file */
    char tspath[512];
    snprintf(tspath, sizeof(tspath), "%s/timeseries.tsv", outdir);
    FILE *fp = fopen(tspath, "w");
    if (!fp) { fprintf(stderr, "FATAL: cannot open %s\n", tspath); exit(1); }
    fprintf(fp, "t\tE_phi_kin\tE_theta_kin\tE_grad\tE_mass\tE_pot\tE_tgrad\tE_tmass\tE_coupling\tE_pot_theta\tE_total\ttheta_rms\tQ_max\n");

    int n_steps = (int)(T / g->dt);
    int diag_every = (int)(diag_dt / g->dt); if (diag_every<1) diag_every=1;

    printf("Steps=%d diag_every=%d\n\n", n_steps, diag_every);

    double wall0 = omp_get_wtime();
    double E0 = 0;

    for (int step = 0; step <= n_steps; step++) {
        if (step > 0) verlet_step(g);
        double t = step * g->dt;

        if (step % diag_every == 0) {
            double epk, etk, eg, em, ep, etg, etm, ec, ept, et;
            compute_energy(g, &epk, &etk, &eg, &em, &ep, &etg, &etm, &ec, &ept, &et);
            if (step == 0) E0 = et;
            double trms = theta_rms(g);
            double qmax = theta_Q_max(g);

            fprintf(fp, "%.2f\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\n",
                    t, epk, etk, eg, em, ep, etg, etm, ec, ept, et, trms, qmax);
            fflush(fp);

            if (step % (diag_every*10) == 0) {
                double wall = omp_get_wtime() - wall0;
                double drift = 100.0*(et-E0)/(fabs(E0)+1e-30);
                printf("t=%7.1f E=%.4e (drift %+.3f%%) theta_rms=%.3e Q_max=%.3e Ept=%.2e [%.0f%% %.0fs]\n",
                       t, et, drift, trms, qmax, ept, 100.0*step/n_steps, wall);
                fflush(stdout);
            }
        }
    }

    fclose(fp);

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

    grid_free(g);
    return 0;
}
