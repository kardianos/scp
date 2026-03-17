/*
 * maxwell_d.c — Elastic dislocation analysis of three-field oscillon
 *
 * Phase 1: 1D oscillon strain analysis
 *   - Evolves triad oscillon (same physics as v21/triad1d.c)
 *   - Computes strain eps_{ax} = d_x phi_a
 *   - Dislocation density trivially zero in 1D (needs >=2 spatial dims)
 *
 * Phase 2: 2D oscillon strain topology
 *   - 2D triad oscillon on Nx x Ny grid
 *   - Strain tensor eps_{ij} = (1/2)(d_i phi_j + d_j phi_i)  [3 fields = "displacement"]
 *   - In 2D with 3 fields: phi_a(x,y), a=1,2,3
 *     eps_{ax,bx} involves d_x phi_a and d_x phi_b etc.
 *   - Incompatibility (dislocation density):
 *     For 2D: alpha_a = eps_{kl} d_k eps_{la}  (Nye tensor, 2D version)
 *     In 2D: alpha_a = d_x eps_{ya} - d_y eps_{xa}
 *            = d_x [(d_y phi_a + d_a phi_y?)]  ... but phi_a is a SCALAR field
 *
 *   KEY INTERPRETATION: phi_a(x,y) is the a-th component of displacement u_a.
 *   We have 3 displacement components and 2 spatial dimensions.
 *   Strain: eps_{ia} = d_i phi_a  (unsymmetrized distortion tensor)
 *   Incompatibility (2D Nye tensor):
 *     alpha_a = (curl eps)_a = d_x eps_{ya} - d_y eps_{xa}
 *             = d_x d_y phi_a - d_y d_x phi_a
 *   For smooth fields: this is ZERO (equality of mixed partials).
 *
 *   But the SYMMETRIZED strain incompatibility (Saint-Venant):
 *     eta = d_xx eps_{yy} + d_yy eps_{xx} - 2 d_xy eps_{xy}
 *   where eps_{ij} = (1/2)(d_i u_j + d_j u_i) sums over field components:
 *     eps_{xx} = sum_a (d_x phi_a)^2 ... no, that's different.
 *
 *   Actually with 3 fields as displacement:
 *     eps_{ij} = (1/2) sum_a (d_i phi_a)(d_j phi_a)   [metric strain, nonlinear]
 *   Or linear:
 *     eps_{ia} = d_i phi_a  (distortion, 2x3 tensor)
 *   Nye dislocation density:
 *     alpha_{3a} = eps_{3kl} d_k eps_{la} = d_1 eps_{2a} - d_2 eps_{1a}
 *                = d_x(d_y phi_a) - d_y(d_x phi_a)
 *   This is zero for smooth C^2 fields.
 *
 *   So we compute BOTH:
 *   (a) The distortion-based Nye tensor (should be zero for smooth oscillon)
 *   (b) The nonlinear metric strain Saint-Venant incompatibility
 *
 * Compile: gcc -O3 -Wall -o maxwell_d src/maxwell_d.c -lm
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ========== PARAMETERS ========== */

/* 1D parameters (v21 defaults) */
static double mu_1d    = -10.0;
static double kappa_1d = 10.0;

/* 2D parameters (stronger coupling for stability, from v21 3D findings) */
static double mu_2d    = -20.0;
static double kappa_2d = 20.0;

/* Shared */
static double mass   = 1.0;
static double A_init = 0.8;
static double sigma  = 3.0;

/* 1D grid */
static int    Nx_1d  = 2000;
static double xmax_1d = 60.0;
static double t1d    = 500.0;   /* enough for equilibration + measurement */

/* 2D grid */
static int    Nx_2d  = 256;
static int    Ny_2d  = 256;
static double L_2d   = 20.0;
static double t2d    = 200.0;

static char outdir[512] = "v24/maxwell_d/data";
static int  phase = 0;  /* 0 = both, 1 = 1D only, 2 = 2D only */

static void parse_args(int argc, char **argv)
{
    for (int i = 1; i < argc - 1; i += 2) {
        if      (!strcmp(argv[i], "-mu1d"))    mu_1d    = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-kappa1d")) kappa_1d = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-mu2d"))    mu_2d    = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-kappa2d")) kappa_2d = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-mass"))    mass     = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-A"))       A_init   = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-sigma"))   sigma    = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-Nx1d"))    Nx_1d    = atoi(argv[i+1]);
        else if (!strcmp(argv[i], "-xmax"))    xmax_1d  = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-t1d"))     t1d      = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-Nx2d"))    Nx_2d    = atoi(argv[i+1]);
        else if (!strcmp(argv[i], "-Ny2d"))    Ny_2d    = atoi(argv[i+1]);
        else if (!strcmp(argv[i], "-L2d"))     L_2d     = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-t2d"))     t2d      = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-phase"))   phase    = atoi(argv[i+1]);
        else if (!strcmp(argv[i], "-o"))       strncpy(outdir, argv[i+1], sizeof(outdir)-1);
        else { fprintf(stderr, "Unknown flag: %s\n", argv[i]); exit(1); }
    }
}

/* -dV/dphi_a where V = (mu/2)P^2/(1+kappa P^2), P = phi1*phi2*phi3 */
static inline double force_pot(double p1, double p2, double p3, int a,
                               double mu, double kap)
{
    double P  = p1 * p2 * p3;
    double P2 = P * P;
    double denom2 = (1.0 + kap * P2) * (1.0 + kap * P2);
    double dP;
    switch (a) {
        case 0: dP = p2 * p3; break;
        case 1: dP = p1 * p3; break;
        case 2: dP = p1 * p2; break;
        default: dP = 0.0;
    }
    return -mu * P * dP / denom2;
}


/* ================================================================
 *  PHASE 1: 1D strain analysis
 * ================================================================ */
static void run_phase1(void)
{
    printf("=== PHASE 1: 1D Oscillon Strain Analysis ===\n");
    printf("  mu=%.1f kappa=%.1f mass=%.3f A=%.3f sigma=%.3f\n",
           mu_1d, kappa_1d, mass, A_init, sigma);

    int Nx = Nx_1d;
    double dx = 2.0 * xmax_1d / (Nx - 1);
    double dx2 = dx * dx;
    double m2 = mass * mass;

    /* CFL */
    double kmax = M_PI / dx;
    double dt = 0.8 * 2.0 / sqrt(kmax*kmax + m2);
    int Nt = (int)(t1d / dt) + 1;

    printf("  Nx=%d dx=%.5f dt=%.6f Nt=%d tfinal=%.0f\n", Nx, dx, dt, Nt, t1d);

    double *phi[3], *vel[3], *acc[3];
    for (int a = 0; a < 3; a++) {
        phi[a] = calloc(Nx, sizeof(double));
        vel[a] = calloc(Nx, sizeof(double));
        acc[a] = calloc(Nx, sizeof(double));
    }

    /* Absorbing boundary */
    double *damp = malloc(Nx * sizeof(double));
    double x_abs = xmax_1d * 0.75;
    for (int i = 0; i < Nx; i++) {
        double x = -xmax_1d + i * dx;
        double ax = fabs(x);
        if (ax > x_abs) {
            double f = (ax - x_abs) / (xmax_1d - x_abs);
            damp[i] = 1.0 - 0.98 * f * f;
        } else {
            damp[i] = 1.0;
        }
    }

    /* Initialize: symmetric Gaussians */
    for (int a = 0; a < 3; a++)
        for (int i = 0; i < Nx; i++) {
            double x = -xmax_1d + i * dx;
            phi[a][i] = A_init * exp(-x * x / (2.0 * sigma * sigma));
        }

    /* Compute acceleration */
    #define ACC1D() do { \
        for (int a = 0; a < 3; a++) { \
            acc[a][0] = acc[a][1] = acc[a][Nx-2] = acc[a][Nx-1] = 0; \
            for (int i = 1; i < Nx - 1; i++) { \
                double lapl = (phi[a][i+1] - 2.0*phi[a][i] + phi[a][i-1]) / dx2; \
                double fp = force_pot(phi[0][i], phi[1][i], phi[2][i], a, mu_1d, kappa_1d); \
                acc[a][i] = lapl - m2*phi[a][i] + fp; \
            } \
        } \
    } while(0)

    ACC1D();

    /* Time series */
    char tspath[600];
    snprintf(tspath, sizeof(tspath), "%s/phase1_ts.tsv", outdir);
    FILE *fts = fopen(tspath, "w");
    if (!fts) { fprintf(stderr, "Cannot open %s\n", tspath); return; }
    fprintf(fts, "time\tphi1_0\tpeak1\tE_total\tE_pot\n");

    int ic = Nx / 2;
    int print_every = Nt / 20;
    if (print_every < 1) print_every = 1;
    int rec_every = Nt / 5000;
    if (rec_every < 1) rec_every = 1;

    /* Snapshot times for strain profile */
    double snap_times[] = {0, 100, 250, 400, 499};
    int n_snaps = 5;
    int next_snap = 0;

    for (int n = 0; n <= Nt; n++) {
        double t = n * dt;

        /* Strain profile snapshots */
        if (next_snap < n_snaps && t >= snap_times[next_snap] - 0.5*dt) {
            char spath[600];
            snprintf(spath, sizeof(spath), "%s/strain_1d_t%04d.tsv",
                     outdir, (int)snap_times[next_snap]);
            FILE *fs = fopen(spath, "w");
            if (fs) {
                fprintf(fs, "x\tphi1\tphi2\tphi3\teps_1x\teps_2x\teps_3x\t"
                            "eps_mag\tP\n");
                for (int i = 2; i < Nx - 2; i++) {
                    double x = -xmax_1d + i * dx;
                    /* 4th-order central difference for strain */
                    double eps[3];
                    for (int a = 0; a < 3; a++)
                        eps[a] = (-phi[a][i+2] + 8.0*phi[a][i+1]
                                  - 8.0*phi[a][i-1] + phi[a][i-2]) / (12.0*dx);
                    double emag = sqrt(eps[0]*eps[0] + eps[1]*eps[1] + eps[2]*eps[2]);
                    double P = phi[0][i] * phi[1][i] * phi[2][i];
                    fprintf(fs, "%.6f\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\n",
                            x, phi[0][i], phi[1][i], phi[2][i],
                            eps[0], eps[1], eps[2], emag, P);
                }
                fclose(fs);
                printf("  Strain snapshot t=%.0f: %s\n", snap_times[next_snap], spath);
            }
            next_snap++;
        }

        /* Time series */
        int do_rec = (n % rec_every == 0);
        int do_print = (n % print_every == 0);

        if (do_rec || do_print) {
            double Ek = 0, Eg = 0, Em = 0, Ep = 0;
            double peak = 0;
            for (int i = 1; i < Nx - 1; i++) {
                for (int a = 0; a < 3; a++) {
                    Ek += 0.5 * vel[a][i] * vel[a][i] * dx;
                    double dp = (phi[a][i+1] - phi[a][i-1]) / (2.0*dx);
                    Eg += 0.5 * dp * dp * dx;
                    Em += 0.5 * m2 * phi[a][i] * phi[a][i] * dx;
                    if (fabs(phi[a][i]) > peak) peak = fabs(phi[a][i]);
                }
                double P = phi[0][i] * phi[1][i] * phi[2][i];
                double P2 = P * P;
                Ep += 0.5 * mu_1d * P2 / (1.0 + kappa_1d * P2) * dx;
            }
            double Et = Ek + Eg + Em + Ep;

            if (do_rec)
                fprintf(fts, "%.4f\t%.6e\t%.6e\t%.6e\t%.6e\n",
                        t, phi[0][ic], peak, Et, Ep);
            if (do_print)
                printf("  t=%7.1f  phi0=%.4f  peak=%.4f  E=%.4f  Ep=%.4f\n",
                       t, phi[0][ic], peak, Et, Ep);
        }

        if (n == Nt) break;

        /* Velocity Verlet */
        for (int a = 0; a < 3; a++)
            for (int i = 1; i < Nx - 1; i++)
                vel[a][i] += 0.5 * dt * acc[a][i];
        for (int a = 0; a < 3; a++)
            for (int i = 1; i < Nx - 1; i++)
                phi[a][i] += dt * vel[a][i];
        ACC1D();
        for (int a = 0; a < 3; a++)
            for (int i = 1; i < Nx - 1; i++)
                vel[a][i] += 0.5 * dt * acc[a][i];

        /* Absorbing boundary */
        for (int a = 0; a < 3; a++)
            for (int i = 0; i < Nx; i++) {
                vel[a][i] *= damp[i];
                phi[a][i] *= damp[i];
            }
    }

    fclose(fts);

    /* Phase 1 summary: compute strain statistics at final time */
    printf("\n--- Phase 1 Strain Summary (t=%.0f) ---\n", t1d);
    double max_eps[3] = {0}, max_phi[3] = {0};
    for (int i = 2; i < Nx - 2; i++) {
        for (int a = 0; a < 3; a++) {
            double eps = (-phi[a][i+2] + 8.0*phi[a][i+1]
                          - 8.0*phi[a][i-1] + phi[a][i-2]) / (12.0*dx);
            if (fabs(eps) > max_eps[a]) max_eps[a] = fabs(eps);
            if (fabs(phi[a][i]) > max_phi[a]) max_phi[a] = fabs(phi[a][i]);
        }
    }
    printf("  Peak |phi|: (%.4f, %.4f, %.4f)\n", max_phi[0], max_phi[1], max_phi[2]);
    printf("  Peak |eps_{ax}|: (%.4f, %.4f, %.4f)\n", max_eps[0], max_eps[1], max_eps[2]);
    printf("  1D dislocation density: ZERO (trivially, needs >=2D)\n");
    printf("  Output: %s\n\n", tspath);

    for (int a = 0; a < 3; a++) { free(phi[a]); free(vel[a]); free(acc[a]); }
    free(damp);
}


/* ================================================================
 *  PHASE 2: 2D oscillon + strain incompatibility
 * ================================================================ */

#define IDX2(i,j) ((i)*Ny + (j))

static void run_phase2(void)
{
    int Nx = Nx_2d, Ny = Ny_2d;
    double Lx = L_2d, Ly = L_2d;

    printf("=== PHASE 2: 2D Oscillon Strain Topology ===\n");
    printf("  mu=%.1f kappa=%.1f mass=%.3f A=%.3f sigma=%.3f\n",
           mu_2d, kappa_2d, mass, A_init, sigma);

    double dx = 2.0 * Lx / (Nx - 1);
    double dy = 2.0 * Ly / (Ny - 1);
    double dx2 = dx * dx;
    double dy2 = dy * dy;
    double m2 = mass * mass;

    /* CFL */
    double cfl = 1.0 / sqrt(1.0/dx2 + 1.0/dy2 + m2);
    double dt = 0.4 * cfl;
    int Nt = (int)(t2d / dt) + 1;

    int N = Nx * Ny;

    printf("  Nx=%d Ny=%d Lx=%.1f Ly=%.1f dx=%.5f dy=%.5f\n",
           Nx, Ny, Lx, Ly, dx, dy);
    printf("  dt=%.6f Nt=%d tfinal=%.0f\n", dt, Nt, t2d);
    printf("  Memory: %.1f MB\n", 9.0 * N * sizeof(double) / (1024.0*1024.0));

    double *phi[3], *vel[3], *acc[3];
    for (int a = 0; a < 3; a++) {
        phi[a] = calloc(N, sizeof(double));
        vel[a] = calloc(N, sizeof(double));
        acc[a] = calloc(N, sizeof(double));
        if (!phi[a] || !vel[a] || !acc[a]) {
            fprintf(stderr, "Allocation failed\n"); return;
        }
    }

    /* Absorbing boundary */
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
            damp[IDX2(i,j)] = 1.0 - 0.98 * f * f;
        }
    }

    /* Initialize: all three fields as identical 2D Gaussians centered at origin
     * This is the radially symmetric oscillon (no vortex topology) */
    for (int a = 0; a < 3; a++)
        for (int i = 0; i < Nx; i++) {
            double x = -Lx + i * dx;
            for (int j = 0; j < Ny; j++) {
                double y = -Ly + j * dy;
                double r2 = x*x + y*y;
                phi[a][IDX2(i,j)] = A_init * exp(-r2 / (2.0 * sigma * sigma));
            }
        }

    /* Compute acceleration */
    #define ACC2D() do { \
        for (int a = 0; a < 3; a++) { \
            for (int i = 0; i < Nx; i++) { \
                acc[a][IDX2(i,0)] = 0; \
                acc[a][IDX2(i,Ny-1)] = 0; \
            } \
            for (int j = 0; j < Ny; j++) { \
                acc[a][IDX2(0,j)] = 0; \
                acc[a][IDX2(Nx-1,j)] = 0; \
            } \
            for (int i = 1; i < Nx-1; i++) \
                for (int j = 1; j < Ny-1; j++) { \
                    double lapl = (phi[a][IDX2(i+1,j)] + phi[a][IDX2(i-1,j)] - 2.0*phi[a][IDX2(i,j)]) / dx2 \
                                + (phi[a][IDX2(i,j+1)] + phi[a][IDX2(i,j-1)] - 2.0*phi[a][IDX2(i,j)]) / dy2; \
                    double fp = force_pot(phi[0][IDX2(i,j)], phi[1][IDX2(i,j)], phi[2][IDX2(i,j)], a, mu_2d, kappa_2d); \
                    acc[a][IDX2(i,j)] = lapl - m2*phi[a][IDX2(i,j)] + fp; \
                } \
        } \
    } while(0)

    ACC2D();

    /* Time series file */
    char tspath[600];
    snprintf(tspath, sizeof(tspath), "%s/phase2_ts.tsv", outdir);
    FILE *fts = fopen(tspath, "w");
    if (!fts) { fprintf(stderr, "Cannot open %s\n", tspath); return; }
    fprintf(fts, "time\tphi1_center\tpeak1\tE_total\tE_pot\tf_core\t"
                 "max_nye_1\tmax_nye_2\tmax_nye_3\tmax_nye_mag\t"
                 "max_SV_incompat\tburgers_1\tburgers_2\tburgers_3\n");

    int print_every = Nt / 20;
    if (print_every < 1) print_every = 1;
    int rec_every = Nt / 2000;
    if (rec_every < 1) rec_every = 1;

    int ic = Nx / 2, jc = Ny / 2;
    double core_r = 3.0 * sigma;

    /* Snapshot times */
    double snap_times[] = {0, 50, 100, 150, 199};
    int n_snaps = 5;
    int next_snap = 0;

    /* Temporary arrays for strain diagnostics */
    double *nye[3];       /* Nye dislocation density alpha_a */
    double *sv_incompat;  /* Saint-Venant incompatibility */
    for (int a = 0; a < 3; a++)
        nye[a] = calloc(N, sizeof(double));
    sv_incompat = calloc(N, sizeof(double));

    for (int n = 0; n <= Nt; n++) {
        double t = n * dt;

        int do_rec = (n % rec_every == 0);
        int do_print = (n % print_every == 0);
        int do_snap = (next_snap < n_snaps && t >= snap_times[next_snap] - 0.5*dt);

        if (do_rec || do_print || do_snap) {
            /* Compute strain-based diagnostics on interior */

            /* Distortion: beta_{ia} = d_i phi_a
             * Nye tensor (2D, out-of-plane component):
             *   alpha_a = d_x(d_y phi_a) - d_y(d_x phi_a)
             * For smooth fields this should be ~0 (up to discretization error).
             *
             * Saint-Venant incompatibility:
             *   eta = d_yy eps_xx + d_xx eps_yy - 2 d_xy eps_xy
             * where for each field component a:
             *   eps_xx^{(a)} = d_x phi_a,  eps_yy^{(a)} = d_y phi_a
             *   eps_xy^{(a)} = (1/2)(d_y phi_a ... )
             *
             * Actually the correct definition with multiple field components:
             * The distortion tensor is beta_{ia} = d_i phi_a  (i=spatial, a=field)
             * The (symmetric) strain is eps_{ij} = (1/2)(beta_{ia}delta_{ja} + beta_{ja}delta_{ia})
             *   ... only makes sense when field index = spatial index.
             *
             * With 3 fields and 2 spatial dims, the "elastic" interpretation:
             *   u_1(x,y) = phi_1, u_2(x,y) = phi_2 are in-plane displacements
             *   u_3(x,y) = phi_3 is out-of-plane displacement
             *
             * In-plane strain (2x2):
             *   eps_xx = d_x phi_1
             *   eps_yy = d_y phi_2
             *   eps_xy = (1/2)(d_y phi_1 + d_x phi_2)
             *
             * Out-of-plane:
             *   eps_xz = (1/2) d_x phi_3,  eps_yz = (1/2) d_y phi_3
             *
             * Saint-Venant (in-plane, 2D):
             *   eta = d_yy(eps_xx) + d_xx(eps_yy) - 2 d_xy(eps_xy)
             *       = d_yy(d_x phi_1) + d_xx(d_y phi_2) - d_xy(d_y phi_1 + d_x phi_2)
             *       = d_xyy phi_1 + d_xxy phi_2 - d_xyy phi_1 - d_xxy phi_2
             *       = 0  (for smooth fields!)
             *
             * The Nye tensor for the full 3D distortion:
             *   alpha_{3a} = d_1 beta_{2a} - d_2 beta_{1a}
             *              = d_x(d_y phi_a) - d_y(d_x phi_a) = 0 for smooth fields
             *
             * CONCLUSION: For smooth C^2 fields, ALL incompatibility measures vanish.
             * This is a fundamental result: oscillons (smooth field configurations)
             * have ZERO dislocation content. Only singular configurations (vortices,
             * domain walls with discontinuous derivatives) carry dislocations.
             *
             * We compute the Nye tensor numerically to verify this is zero to
             * discretization accuracy, and also compute the nonlinear Burgers integral.
             */

            /* Compute Nye tensor: alpha_a = d_x(d_y phi_a) - d_y(d_x phi_a) */
            for (int a = 0; a < 3; a++) {
                for (int i = 0; i < Nx; i++)
                    for (int j = 0; j < Ny; j++)
                        nye[a][IDX2(i,j)] = 0;

                for (int i = 2; i < Nx-2; i++)
                    for (int j = 2; j < Ny-2; j++) {
                        /* d_x(d_y phi) using centered differences */
                        /* d_y phi at (i+1,j) and (i-1,j) */
                        double dyphi_ip = (phi[a][IDX2(i+1,j+1)] - phi[a][IDX2(i+1,j-1)]) / (2.0*dy);
                        double dyphi_im = (phi[a][IDX2(i-1,j+1)] - phi[a][IDX2(i-1,j-1)]) / (2.0*dy);
                        double dx_dyphi = (dyphi_ip - dyphi_im) / (2.0*dx);

                        /* d_y(d_x phi) */
                        double dxphi_jp = (phi[a][IDX2(i+1,j+1)] - phi[a][IDX2(i-1,j+1)]) / (2.0*dx);
                        double dxphi_jm = (phi[a][IDX2(i+1,j-1)] - phi[a][IDX2(i-1,j-1)]) / (2.0*dx);
                        double dy_dxphi = (dxphi_jp - dxphi_jm) / (2.0*dy);

                        nye[a][IDX2(i,j)] = dx_dyphi - dy_dxphi;
                    }
            }

            /* Saint-Venant incompatibility:
             * eta = d_yy(d_x phi_1) + d_xx(d_y phi_2) - d_xy(d_y phi_1 + d_x phi_2)
             * Compute numerically with second-order finite differences */
            for (int i = 0; i < Nx; i++)
                for (int j = 0; j < Ny; j++)
                    sv_incompat[IDX2(i,j)] = 0;

            for (int i = 2; i < Nx-2; i++)
                for (int j = 2; j < Ny-2; j++) {
                    /* eps_xx = d_x phi_1 at various grid points */
                    double exx_c = (phi[0][IDX2(i+1,j)] - phi[0][IDX2(i-1,j)]) / (2.0*dx);
                    double exx_jp = (phi[0][IDX2(i+1,j+1)] - phi[0][IDX2(i-1,j+1)]) / (2.0*dx);
                    double exx_jm = (phi[0][IDX2(i+1,j-1)] - phi[0][IDX2(i-1,j-1)]) / (2.0*dx);
                    double d_yy_exx = (exx_jp - 2.0*exx_c + exx_jm) / dy2;

                    /* eps_yy = d_y phi_2 at various grid points */
                    double eyy_c = (phi[1][IDX2(i,j+1)] - phi[1][IDX2(i,j-1)]) / (2.0*dy);
                    double eyy_ip = (phi[1][IDX2(i+1,j+1)] - phi[1][IDX2(i+1,j-1)]) / (2.0*dy);
                    double eyy_im = (phi[1][IDX2(i-1,j+1)] - phi[1][IDX2(i-1,j-1)]) / (2.0*dy);
                    double d_xx_eyy = (eyy_ip - 2.0*eyy_c + eyy_im) / dx2;

                    /* eps_xy = (1/2)(d_y phi_1 + d_x phi_2) at various grid points */
                    /* Need d_xy(eps_xy) = mixed partial */
                    double exy_pp = 0.5*((phi[0][IDX2(i+1,j+2)] - phi[0][IDX2(i+1,j)]) / (2.0*dy)
                                       + (phi[1][IDX2(i+2,j+1)] - phi[1][IDX2(i,j+1)]) / (2.0*dx));
                    double exy_pm = 0.5*((phi[0][IDX2(i+1,j)] - phi[0][IDX2(i+1,j-2)]) / (2.0*dy)
                                       + (phi[1][IDX2(i+2,j-1)] - phi[1][IDX2(i,j-1)]) / (2.0*dx));
                    double exy_mp = 0.5*((phi[0][IDX2(i-1,j+2)] - phi[0][IDX2(i-1,j)]) / (2.0*dy)
                                       + (phi[1][IDX2(i,j+1)] - phi[1][IDX2(i-2,j+1)]) / (2.0*dx));
                    double exy_mm = 0.5*((phi[0][IDX2(i-1,j)] - phi[0][IDX2(i-1,j-2)]) / (2.0*dy)
                                       + (phi[1][IDX2(i,j-1)] - phi[1][IDX2(i-2,j-1)]) / (2.0*dx));
                    double d_xy_exy = (exy_pp - exy_pm - exy_mp + exy_mm) / (4.0*dx*dy);

                    sv_incompat[IDX2(i,j)] = d_yy_exx + d_xx_eyy - 2.0*d_xy_exy;
                }

            /* Burgers vector: line integral of distortion around a large contour
             * b_a = oint beta_{ia} dl_i = oint (d_i phi_a) dl_i
             * For a closed contour: b_a = 0 for smooth fields (by Stokes' theorem)
             * Compute on a square contour around the oscillon core */
            double burgers[3] = {0, 0, 0};
            int r_cont = (int)(2.0 * sigma / dx);  /* contour at ~2 sigma */
            if (r_cont > 2 && ic - r_cont > 1 && ic + r_cont < Nx - 1 &&
                jc - r_cont > 1 && jc + r_cont < Ny - 1) {
                for (int a = 0; a < 3; a++) {
                    double bval = 0;
                    /* Bottom edge: y = jc - r_cont, x goes right */
                    int jbot = jc - r_cont;
                    for (int i = ic - r_cont; i < ic + r_cont; i++) {
                        double dxphi = (phi[a][IDX2(i+1,jbot)] - phi[a][IDX2(i-1,jbot)]) / (2.0*dx);
                        bval += dxphi * dx;
                    }
                    /* Right edge: x = ic + r_cont, y goes up */
                    int iright = ic + r_cont;
                    for (int j = jbot; j < jc + r_cont; j++) {
                        double dyphi = (phi[a][IDX2(iright,j+1)] - phi[a][IDX2(iright,j-1)]) / (2.0*dy);
                        bval += dyphi * dy;
                    }
                    /* Top edge: y = jc + r_cont, x goes left */
                    int jtop = jc + r_cont;
                    for (int i = ic + r_cont; i > ic - r_cont; i--) {
                        double dxphi = (phi[a][IDX2(i+1,jtop)] - phi[a][IDX2(i-1,jtop)]) / (2.0*dx);
                        bval -= dxphi * dx;
                    }
                    /* Left edge: x = ic - r_cont, y goes down */
                    int ileft = ic - r_cont;
                    for (int j = jc + r_cont; j > jbot; j--) {
                        double dyphi = (phi[a][IDX2(ileft,j+1)] - phi[a][IDX2(ileft,j-1)]) / (2.0*dy);
                        bval -= dyphi * dy;
                    }
                    burgers[a] = bval;
                }
            }

            /* Compute statistics */
            double max_nye[3] = {0, 0, 0};
            double max_nye_mag = 0;
            double max_sv = 0;
            for (int i = 3; i < Nx-3; i++)
                for (int j = 3; j < Ny-3; j++) {
                    double nm = 0;
                    for (int a = 0; a < 3; a++) {
                        double v = fabs(nye[a][IDX2(i,j)]);
                        if (v > max_nye[a]) max_nye[a] = v;
                        nm += nye[a][IDX2(i,j)] * nye[a][IDX2(i,j)];
                    }
                    nm = sqrt(nm);
                    if (nm > max_nye_mag) max_nye_mag = nm;
                    double sv = fabs(sv_incompat[IDX2(i,j)]);
                    if (sv > max_sv) max_sv = sv;
                }

            /* Energy */
            double Ek = 0, Eg = 0, Em = 0, Ep = 0, Ecore = 0, Eall = 0;
            double peak = 0;
            for (int i = 1; i < Nx-1; i++)
                for (int j = 1; j < Ny-1; j++) {
                    int k = IDX2(i,j);
                    double x = -Lx + i * dx;
                    double y = -Ly + j * dy;
                    double e = 0;
                    for (int a = 0; a < 3; a++) {
                        double ek = 0.5 * vel[a][k] * vel[a][k];
                        double dpx = (phi[a][IDX2(i+1,j)] - phi[a][IDX2(i-1,j)]) / (2.0*dx);
                        double dpy = (phi[a][IDX2(i,j+1)] - phi[a][IDX2(i,j-1)]) / (2.0*dy);
                        double eg = 0.5 * (dpx*dpx + dpy*dpy);
                        double em = 0.5 * m2 * phi[a][k] * phi[a][k];
                        Ek += ek * dx * dy;
                        Eg += eg * dx * dy;
                        Em += em * dx * dy;
                        e += ek + eg + em;
                        if (fabs(phi[a][k]) > peak) peak = fabs(phi[a][k]);
                    }
                    double P = phi[0][k] * phi[1][k] * phi[2][k];
                    double P2 = P * P;
                    double V = 0.5 * mu_2d * P2 / (1.0 + kappa_2d * P2);
                    Ep += V * dx * dy;
                    e += V;
                    double dV = dx * dy;
                    Eall += e * dV;
                    double r = sqrt(x*x + y*y);
                    if (r < core_r) Ecore += e * dV;
                }
            double Et = Ek + Eg + Em + Ep;
            double fc = (Eall > 1e-20) ? Ecore / Eall : 0.0;

            if (do_rec)
                fprintf(fts, "%.4f\t%.6e\t%.6e\t%.6e\t%.6e\t%.4f\t"
                             "%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\n",
                        t, phi[0][IDX2(ic,jc)], peak, Et, Ep, fc,
                        max_nye[0], max_nye[1], max_nye[2], max_nye_mag,
                        max_sv, burgers[0], burgers[1], burgers[2]);

            if (do_print)
                printf("  t=%6.1f  phi0=%.4f  pk=%.4f  E=%.3f  fc=%.3f  "
                       "|nye|=%.2e  |SV|=%.2e  b=(%.2e,%.2e,%.2e)\n",
                       t, phi[0][IDX2(ic,jc)], peak, Et, fc,
                       max_nye_mag, max_sv,
                       burgers[0], burgers[1], burgers[2]);

            /* Write 2D snapshot */
            if (do_snap) {
                char spath[600];
                snprintf(spath, sizeof(spath), "%s/strain_2d_t%04d.tsv",
                         outdir, (int)snap_times[next_snap]);
                FILE *fs = fopen(spath, "w");
                if (fs) {
                    fprintf(fs, "x\ty\tphi1\tphi2\tphi3\t"
                                "eps_xx\teps_yy\teps_xy\t"
                                "nye1\tnye2\tnye3\tSV_incompat\n");
                    /* Subsample for output */
                    int step = 2;
                    for (int i = 2; i < Nx-2; i += step)
                        for (int j = 2; j < Ny-2; j += step) {
                            double x = -Lx + i * dx;
                            double y = -Ly + j * dy;
                            /* In-plane strain components */
                            double eps_xx = (phi[0][IDX2(i+1,j)] - phi[0][IDX2(i-1,j)]) / (2.0*dx);
                            double eps_yy = (phi[1][IDX2(i,j+1)] - phi[1][IDX2(i,j-1)]) / (2.0*dy);
                            double eps_xy = 0.5 * (
                                (phi[0][IDX2(i,j+1)] - phi[0][IDX2(i,j-1)]) / (2.0*dy) +
                                (phi[1][IDX2(i+1,j)] - phi[1][IDX2(i-1,j)]) / (2.0*dx));
                            fprintf(fs, "%.4f\t%.4f\t%.6e\t%.6e\t%.6e\t"
                                        "%.6e\t%.6e\t%.6e\t"
                                        "%.6e\t%.6e\t%.6e\t%.6e\n",
                                    x, y,
                                    phi[0][IDX2(i,j)], phi[1][IDX2(i,j)], phi[2][IDX2(i,j)],
                                    eps_xx, eps_yy, eps_xy,
                                    nye[0][IDX2(i,j)], nye[1][IDX2(i,j)], nye[2][IDX2(i,j)],
                                    sv_incompat[IDX2(i,j)]);
                        }
                    fclose(fs);
                    printf("  Snapshot t=%.0f: %s\n", snap_times[next_snap], spath);
                }
                next_snap++;
            }
        }

        if (n == Nt) break;

        /* Velocity Verlet */
        for (int a = 0; a < 3; a++)
            for (int k = 0; k < N; k++)
                vel[a][k] += 0.5 * dt * acc[a][k];
        for (int a = 0; a < 3; a++)
            for (int k = 0; k < N; k++)
                phi[a][k] += dt * vel[a][k];
        ACC2D();
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

    printf("\n--- Phase 2 Summary ---\n");
    printf("  Output: %s\n\n", tspath);

    for (int a = 0; a < 3; a++) { free(phi[a]); free(vel[a]); free(acc[a]); free(nye[a]); }
    free(damp); free(sv_incompat);
}


/* ================================================================ */
int main(int argc, char **argv)
{
    parse_args(argc, argv);

    printf("maxwell_d: Elastic Dislocation Analysis of Three-Field Oscillon\n\n");

    if (phase == 0 || phase == 1)
        run_phase1();
    if (phase == 0 || phase == 2)
        run_phase2();

    return 0;
}
