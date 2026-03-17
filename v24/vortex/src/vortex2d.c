/*
 * vortex2d.c — 2D three-field vortex pair: massive scalars with triple-product coupling
 *
 * Lagrangian (2D):
 *   L = sum_a [ (1/2)(dt phi_a)^2 - (1/2)|grad phi_a|^2 - (m^2/2)phi_a^2 ]
 *     - (mu/2) P^2 / (1 + kappa P^2)
 *   P = phi_1 phi_2 phi_3
 *
 * Fields phi_1, phi_2 initialized as vortex-antivortex pair (Lamb-Oseen profile).
 * phi_3 initialized with small random noise to seed triple-product effects.
 *
 * Key question: does the vortex pair survive and propagate in the massive three-field system?
 *
 * Compile: gcc -O3 -Wall -o vortex2d src/vortex2d.c -lm
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Parameters */
static double mu     = -20.0;
static double kappa  = 20.0;
static double mass   = 1.0;
static double Gamma  = 1.0;   /* circulation strength */
static double a_core = 2.0;   /* vortex core radius */
static double d_sep  = 10.0;  /* pair separation */
static double noise  = 0.01;  /* phi3 noise amplitude */
static int    Nx     = 512;
static int    Ny     = 512;
static double Lx     = 40.0;  /* half-extent: domain is [-Lx, Lx] */
static double Ly     = 40.0;
static double tfinal = 2000.0;
static char   outdir[512] = "data";
static int    seed   = 42;

static void parse_args(int argc, char **argv)
{
    for (int i = 1; i < argc - 1; i += 2) {
        if      (!strcmp(argv[i], "-mu"))     mu     = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-kappa"))  kappa  = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-mass"))   mass   = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-Gamma"))  Gamma  = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-acore"))  a_core = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-dsep"))   d_sep  = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-noise"))  noise  = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-Nx"))     Nx     = atoi(argv[i+1]);
        else if (!strcmp(argv[i], "-Ny"))     Ny     = atoi(argv[i+1]);
        else if (!strcmp(argv[i], "-Lx"))     Lx     = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-Ly"))     Ly     = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-tfinal")) tfinal = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-seed"))   seed   = atoi(argv[i+1]);
        else if (!strcmp(argv[i], "-o"))      strncpy(outdir, argv[i+1], sizeof(outdir)-1);
        else { fprintf(stderr, "Unknown flag: %s\n", argv[i]); exit(1); }
    }
}

/* Index into flat 2D array */
#define IDX(i,j) ((i)*Ny + (j))

/* -dV/dphi_a where V = (mu/2)P^2/(1+kappa P^2), P = phi1*phi2*phi3 */
static inline double force_pot(double p1, double p2, double p3, int a)
{
    double P  = p1 * p2 * p3;
    double P2 = P * P;
    double denom2 = (1.0 + kappa * P2) * (1.0 + kappa * P2);
    double dP;
    switch (a) {
        case 0: dP = p2 * p3; break;
        case 1: dP = p1 * p3; break;
        case 2: dP = p1 * p2; break;
        default: dP = 0.0;
    }
    return -mu * P * dP / denom2;
}

/* Add single Lamb-Oseen vortex contribution at (x0, y0) with sign s */
static void add_vortex(double *phi1, double *phi2, double dx, double dy,
                       double x0, double y0, double sign)
{
    for (int i = 0; i < Nx; i++) {
        double x = -Lx + i * dx;
        for (int j = 0; j < Ny; j++) {
            double y = -Ly + j * dy;
            double rx = x - x0;
            double ry = y - y0;
            double r2 = rx*rx + ry*ry;
            if (r2 < 1e-30) r2 = 1e-30;
            double fac = sign * Gamma * (1.0 - exp(-r2 / (a_core * a_core))) / r2;
            phi1[IDX(i,j)] += -fac * ry;
            phi2[IDX(i,j)] += +fac * rx;
        }
    }
}

/* Simple LCG random in [-1,1] */
static unsigned int rng_state;
static double randf(void)
{
    rng_state = rng_state * 1103515245 + 12345;
    return 2.0 * ((rng_state >> 16) & 0x7FFF) / 32767.0 - 1.0;
}

/* Find vortex centroid in upper half-plane (y > 0) using |curl phi|^2 weighting */
static void find_vortex_centroid(const double *phi1, const double *phi2,
                                  double dx, double dy,
                                  double *cx, double *cy)
{
    double wx = 0, wy = 0, wtot = 0;
    for (int i = 1; i < Nx-1; i++) {
        double x = -Lx + i * dx;
        for (int j = Ny/2; j < Ny-1; j++) {  /* upper half only */
            double y = -Ly + j * dy;
            /* vorticity omega_z = dphi2/dx - dphi1/dy */
            double dp2dx = (phi2[IDX(i+1,j)] - phi2[IDX(i-1,j)]) / (2.0*dx);
            double dp1dy = (phi1[IDX(i,j+1)] - phi1[IDX(i,j-1)]) / (2.0*dy);
            double omega = dp2dx - dp1dy;
            double w = omega * omega;
            wx += w * x;
            wy += w * y;
            wtot += w;
        }
    }
    if (wtot > 1e-30) {
        *cx = wx / wtot;
        *cy = wy / wtot;
    } else {
        *cx = 0;
        *cy = d_sep / 2.0;
    }
}

/* Compute divergence integral = sum |div phi|^2 */
static double compute_div2(const double *phi1, const double *phi2,
                           double dx, double dy)
{
    double s = 0;
    for (int i = 1; i < Nx-1; i++)
        for (int j = 1; j < Ny-1; j++) {
            double dp1dx = (phi1[IDX(i+1,j)] - phi1[IDX(i-1,j)]) / (2.0*dx);
            double dp2dy = (phi2[IDX(i,j+1)] - phi2[IDX(i,j-1)]) / (2.0*dy);
            double div = dp1dx + dp2dy;
            s += div * div * dx * dy;
        }
    return s;
}

/* Compute peak vorticity in upper half-plane */
static double peak_vorticity(const double *phi1, const double *phi2,
                              double dx, double dy)
{
    double peak = 0;
    for (int i = 1; i < Nx-1; i++)
        for (int j = Ny/2; j < Ny-1; j++) {
            double dp2dx = (phi2[IDX(i+1,j)] - phi2[IDX(i-1,j)]) / (2.0*dx);
            double dp1dy = (phi1[IDX(i,j+1)] - phi1[IDX(i,j-1)]) / (2.0*dy);
            double omega = fabs(dp2dx - dp1dy);
            if (omega > peak) peak = omega;
        }
    return peak;
}

int main(int argc, char **argv)
{
    parse_args(argc, argv);

    double dx = 2.0 * Lx / (Nx - 1);
    double dy = 2.0 * Ly / (Ny - 1);
    double dx2 = dx * dx;
    double dy2 = dy * dy;
    double m2 = mass * mass;

    /* CFL for 2D: dt < 1/sqrt(1/dx^2 + 1/dy^2 + m^2) */
    double cfl = 1.0 / sqrt(1.0/dx2 + 1.0/dy2 + m2);
    double dt = 0.4 * cfl;
    int Nt = (int)(tfinal / dt) + 1;

    int N = Nx * Ny;

    printf("vortex2d: 2D vortex-antivortex pair\n");
    printf("  mu=%.3f kappa=%.4f mass=%.4f\n", mu, kappa, mass);
    printf("  Gamma=%.3f a_core=%.3f d_sep=%.3f noise=%.4f\n",
           Gamma, a_core, d_sep, noise);
    printf("  Nx=%d Ny=%d Lx=%.1f Ly=%.1f dx=%.5f dy=%.5f\n",
           Nx, Ny, Lx, Ly, dx, dy);
    printf("  dt=%.6f Nt=%d tfinal=%.0f\n", dt, Nt, tfinal);
    printf("  Grid points: %d  Memory: %.1f MB\n",
           N, 9.0 * N * sizeof(double) / (1024.0*1024.0));

    /* Allocate: 3 fields x (phi, vel, acc) */
    double *phi[3], *vel[3], *acc[3];
    for (int a = 0; a < 3; a++) {
        phi[a] = calloc(N, sizeof(double));
        vel[a] = calloc(N, sizeof(double));
        acc[a] = calloc(N, sizeof(double));
        if (!phi[a] || !vel[a] || !acc[a]) {
            fprintf(stderr, "Allocation failed\n");
            return 1;
        }
    }

    /* Absorbing boundary: damp in outer 25% */
    double *damp = malloc(N * sizeof(double));
    double x_abs = Lx * 0.75;
    double y_abs = Ly * 0.75;
    for (int i = 0; i < Nx; i++) {
        double x = -Lx + i * dx;
        double ax = fabs(x);
        double fx = (ax > x_abs) ? (ax - x_abs) / (Lx - x_abs) : 0.0;
        for (int j = 0; j < Ny; j++) {
            double y = -Ly + j * dy;
            double ay = fabs(y);
            double fy = (ay > y_abs) ? (ay - y_abs) / (Ly - y_abs) : 0.0;
            double f = (fx > fy) ? fx : fy;
            damp[IDX(i,j)] = 1.0 - 0.98 * f * f;
        }
    }

    /* Initialize vortex-antivortex pair:
     * Vortex at (0, +d/2) with +Gamma
     * Antivortex at (0, -d/2) with -Gamma
     * This pair propagates in +x direction */
    add_vortex(phi[0], phi[1], dx, dy, 0.0, +d_sep/2.0, +1.0);
    add_vortex(phi[0], phi[1], dx, dy, 0.0, -d_sep/2.0, -1.0);

    /* phi3 = small random noise */
    rng_state = (unsigned int)seed;
    for (int k = 0; k < N; k++)
        phi[2][k] = noise * randf();

    /* Measure initial fields */
    double peak_phi1 = 0, peak_phi2 = 0;
    for (int k = 0; k < N; k++) {
        if (fabs(phi[0][k]) > peak_phi1) peak_phi1 = fabs(phi[0][k]);
        if (fabs(phi[1][k]) > peak_phi2) peak_phi2 = fabs(phi[1][k]);
    }
    printf("  Initial peak |phi1|=%.6f |phi2|=%.6f\n", peak_phi1, peak_phi2);

    /* Compute acceleration macro */
    #define COMPUTE_ACC() do { \
        for (int a = 0; a < 3; a++) { \
            for (int i = 0; i < Nx; i++) { \
                acc[a][IDX(i,0)] = 0; \
                acc[a][IDX(i,Ny-1)] = 0; \
            } \
            for (int j = 0; j < Ny; j++) { \
                acc[a][IDX(0,j)] = 0; \
                acc[a][IDX(Nx-1,j)] = 0; \
            } \
            for (int i = 1; i < Nx-1; i++) \
                for (int j = 1; j < Ny-1; j++) { \
                    double lapl = (phi[a][IDX(i+1,j)] + phi[a][IDX(i-1,j)] - 2.0*phi[a][IDX(i,j)]) / dx2 \
                                + (phi[a][IDX(i,j+1)] + phi[a][IDX(i,j-1)] - 2.0*phi[a][IDX(i,j)]) / dy2; \
                    double fp = force_pot(phi[0][IDX(i,j)], phi[1][IDX(i,j)], phi[2][IDX(i,j)], a); \
                    acc[a][IDX(i,j)] = lapl - m2*phi[a][IDX(i,j)] + fp; \
                } \
        } \
    } while(0)

    COMPUTE_ACC();

    /* Open time series file */
    char tspath[600];
    snprintf(tspath, sizeof(tspath), "%s/vortex_ts.tsv", outdir);
    FILE *fts = fopen(tspath, "w");
    if (!fts) {
        fprintf(stderr, "Cannot open %s\n", tspath);
        return 1;
    }
    fprintf(fts, "time\tpeak_phi1\tpeak_phi2\tpeak_phi3\t"
                 "E_kin\tE_grad\tE_mass\tE_pot\tE_total\t"
                 "vortex_x\tvortex_y\tpeak_vort\tdiv2\n");

    int rec_every  = Nt / 10000;
    if (rec_every < 1) rec_every = 1;
    int print_every = Nt / 80;
    if (print_every < 1) print_every = 1;

    /* Snapshot times */
    double snap_times[] = {0, 10, 50, 100, 200, 500, 1000, 2000};
    int n_snaps = sizeof(snap_times) / sizeof(snap_times[0]);
    int next_snap = 0;

    for (int n = 0; n <= Nt; n++) {
        double t = n * dt;

        /* Snapshots */
        if (next_snap < n_snaps && t >= snap_times[next_snap] - 0.5*dt) {
            char snappath[600];
            snprintf(snappath, sizeof(snappath),
                     "%s/vortex_profile_t%04d.tsv", outdir, (int)snap_times[next_snap]);
            FILE *fsnap = fopen(snappath, "w");
            if (fsnap) {
                fprintf(fsnap, "x\ty\tphi1\tphi2\tphi3\n");
                /* Subsample: every 4th point */
                int step = 4;
                for (int i = 0; i < Nx; i += step)
                    for (int j = 0; j < Ny; j += step) {
                        double x = -Lx + i * dx;
                        double y = -Ly + j * dy;
                        fprintf(fsnap, "%.4f\t%.4f\t%.6e\t%.6e\t%.6e\n",
                                x, y, phi[0][IDX(i,j)], phi[1][IDX(i,j)], phi[2][IDX(i,j)]);
                    }
                fclose(fsnap);
                printf("  Snapshot t=%g written: %s\n", snap_times[next_snap], snappath);
            }
            next_snap++;
        }

        int do_rec   = (n % rec_every == 0);
        int do_print = (n % print_every == 0);

        if (do_rec || do_print) {
            double Ek = 0, Eg = 0, Em = 0, Ep = 0;
            double peak[3] = {0, 0, 0};

            for (int i = 1; i < Nx-1; i++)
                for (int j = 1; j < Ny-1; j++) {
                    int k = IDX(i,j);
                    for (int a = 0; a < 3; a++) {
                        Ek += 0.5 * vel[a][k] * vel[a][k] * dx * dy;
                        double dpx = (phi[a][IDX(i+1,j)] - phi[a][IDX(i-1,j)]) / (2.0*dx);
                        double dpy = (phi[a][IDX(i,j+1)] - phi[a][IDX(i,j-1)]) / (2.0*dy);
                        Eg += 0.5 * (dpx*dpx + dpy*dpy) * dx * dy;
                        Em += 0.5 * m2 * phi[a][k] * phi[a][k] * dx * dy;
                        if (fabs(phi[a][k]) > peak[a]) peak[a] = fabs(phi[a][k]);
                    }
                    double P = phi[0][k] * phi[1][k] * phi[2][k];
                    double P2 = P * P;
                    double V = 0.5 * mu * P2 / (1.0 + kappa * P2);
                    Ep += V * dx * dy;
                }

            double Et = Ek + Eg + Em + Ep;

            double vx, vy;
            find_vortex_centroid(phi[0], phi[1], dx, dy, &vx, &vy);
            double pvort = peak_vorticity(phi[0], phi[1], dx, dy);
            double d2 = compute_div2(phi[0], phi[1], dx, dy);

            if (do_rec)
                fprintf(fts, "%.4f\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t"
                             "%.4f\t%.4f\t%.6e\t%.6e\n",
                        t, peak[0], peak[1], peak[2],
                        Ek, Eg, Em, Ep, Et,
                        vx, vy, pvort, d2);
            if (do_print)
                printf("  t=%7.1f  pk=(%.4f,%.4f,%.4f)  E=%.3f(k=%.3f g=%.3f m=%.3f p=%.3f)"
                       "  vort@(%.2f,%.2f) |w|=%.4f\n",
                       t, peak[0], peak[1], peak[2],
                       Et, Ek, Eg, Em, Ep,
                       vx, vy, pvort);
        }

        if (n == Nt) break;

        /* Velocity Verlet */
        for (int a = 0; a < 3; a++)
            for (int k = 0; k < N; k++)
                vel[a][k] += 0.5 * dt * acc[a][k];
        for (int a = 0; a < 3; a++)
            for (int k = 0; k < N; k++)
                phi[a][k] += dt * vel[a][k];
        COMPUTE_ACC();
        for (int a = 0; a < 3; a++)
            for (int k = 0; k < N; k++)
                vel[a][k] += 0.5 * dt * acc[a][k];

        /* Absorbing boundary */
        for (int a = 0; a < 3; a++)
            for (int k = 0; k < N; k++) {
                vel[a][k] *= damp[k];
                phi[a][k] *= damp[k];
            }
    }

    fclose(fts);
    printf("\nOutput: %s\n", tspath);

    /* Cleanup */
    for (int a = 0; a < 3; a++) { free(phi[a]); free(vel[a]); free(acc[a]); }
    free(damp);
    return 0;
}
