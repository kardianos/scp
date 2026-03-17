/*
 * proca_spin.c — V24-S6: Spin classification of the Proca mediator
 *
 * Three massive scalar fields with triple-product + pairwise + cross-gradient couplings.
 *
 * Lagrangian (2D):
 *   L = sum_a [ (1/2)(dt phi_a)^2 - (1/2)|grad phi_a|^2 - (m^2/2)phi_a^2 ]
 *     - (mu/2) P^2 / (1 + kappa P^2)                     [triple product]
 *     - (lambda/2) sum_{a<b} (phi_a - phi_b)^2            [pairwise]
 *     + (eta/2) sum_{a,b} (d_a phi_b)(d_b phi_a)          [cross-gradient, a=spatial, b=field]
 *
 * In 2D with 2 fields (phi_1, phi_2): field index = spatial index.
 *   Cross-gradient: eta * [ (dx phi_2)(dy phi_1) + (dy phi_1)(dx phi_2) ]
 *                 = eta * (dx phi_2)(dy phi_1)   [symmetric part, doubled]
 *
 * Mode decomposition at each point:
 *   S  = (phi_1 + phi_2) / sqrt(2)   [symmetric, compression]
 *   A  = (phi_1 - phi_2) / sqrt(2)   [antisymmetric, Proca candidate]
 *
 * Tests:
 *   1: 1D baseline — verify forward/backward symmetry of antisymmetric radiation
 *   2: 2D with pairwise only (lambda=0.99, eta=0) — measure angular pattern
 *   3: 2D with pairwise + cross-gradient (lambda=0.99, eta=0.5) — spin-2?
 *
 * Compile: gcc -O3 -Wall -o proca_spin src/proca_spin.c -lm
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Parameters */
static double mu      = -20.0;
static double kappa   = 20.0;
static double mass    = 1.0;
static double lambda  = 0.99;   /* pairwise coupling */
static double eta     = 0.0;    /* cross-gradient coupling */
static double A_init  = 0.8;
static double sigma   = 3.0;    /* Gaussian width */
static double eps     = 0.05;   /* antisymmetric perturbation amplitude */
static int    Nx      = 256;
static int    Ny      = 256;
static double L       = 20.0;   /* half-extent: domain [-L, L] */
static double tfinal  = 1000.0;
static int    test    = 2;
static double R_meas  = 15.0;   /* measurement radius */
static int    N_ang   = 72;     /* angular bins (every 5 degrees) */
static char   outdir[512] = "v24/proca_spin/data";

static void parse_args(int argc, char **argv)
{
    for (int i = 1; i < argc - 1; i += 2) {
        if      (!strcmp(argv[i], "-mu"))      mu      = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-kappa"))   kappa   = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-mass"))    mass    = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-lambda"))  lambda  = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-eta"))     eta     = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-A"))       A_init  = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-sigma"))   sigma   = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-eps"))     eps     = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-Nx"))      Nx      = atoi(argv[i+1]);
        else if (!strcmp(argv[i], "-Ny"))      Ny      = atoi(argv[i+1]);
        else if (!strcmp(argv[i], "-L"))       L       = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-tfinal"))  tfinal  = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-test"))    test    = atoi(argv[i+1]);
        else if (!strcmp(argv[i], "-Rmeas"))   R_meas  = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-Nang"))    N_ang   = atoi(argv[i+1]);
        else if (!strcmp(argv[i], "-o"))       strncpy(outdir, argv[i+1], sizeof(outdir)-1);
        else { fprintf(stderr, "Unknown flag: %s\n", argv[i]); exit(1); }
    }
}

#define IDX(i,j) ((i)*Ny + (j))

/* -dV/dphi_a for triple product: V = (mu/2)P^2/(1+kappa P^2), P = phi1*phi2*phi3 */
static inline double force_triple(double p1, double p2, double p3, int a)
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

/* ====================================================================
 * TEST 1: 1D baseline — two fields, pairwise coupling
 * ==================================================================== */
static void run_test1(void)
{
    int N1 = 2000;
    double xmax = 60.0;
    double dx = 2.0 * xmax / (N1 - 1);
    double dx2 = dx * dx;
    double m2 = mass * mass;

    double kmax = M_PI / dx;
    double dt = 0.8 * 2.0 / sqrt(kmax*kmax + m2);
    int Nt = (int)(tfinal / dt) + 1;

    printf("Test 1: 1D baseline — forward/backward symmetry\n");
    printf("  N=%d xmax=%.1f dx=%.5f dt=%.6f Nt=%d\n", N1, xmax, dx, dt, Nt);

    /* 3 fields in 1D */
    double *phi[3], *vel[3], *acc[3];
    for (int a = 0; a < 3; a++) {
        phi[a] = calloc(N1, sizeof(double));
        vel[a] = calloc(N1, sizeof(double));
        acc[a] = calloc(N1, sizeof(double));
    }

    /* Absorbing boundary: outer 25% */
    double *damp = malloc(N1 * sizeof(double));
    double x_abs = xmax * 0.75;
    for (int i = 0; i < N1; i++) {
        double x = -xmax + i * dx;
        double ax = fabs(x);
        if (ax > x_abs) {
            double f = (ax - x_abs) / (xmax - x_abs);
            damp[i] = 1.0 - 0.98 * f * f;
        } else {
            damp[i] = 1.0;
        }
    }

    /* Initialize: symmetric Gaussian oscillon + small antisymmetric perturbation */
    int ic = N1 / 2;
    for (int i = 0; i < N1; i++) {
        double x = -xmax + i * dx;
        double g = A_init * exp(-x*x / (2.0*sigma*sigma));
        phi[0][i] = g + eps * g;   /* phi_1 = (1+eps) * g */
        phi[1][i] = g - eps * g;   /* phi_2 = (1-eps) * g */
        phi[2][i] = g;             /* phi_3 = g (needed for triple product) */
    }

    /* Macro: compute 1D acceleration with triple + pairwise + mass */
    #define ACC1D() do { \
        for (int a = 0; a < 3; a++) { \
            acc[a][0] = acc[a][N1-1] = 0; \
            for (int i = 1; i < N1-1; i++) { \
                double lapl = (phi[a][i+1] - 2.0*phi[a][i] + phi[a][i-1]) / dx2; \
                double fp = force_triple(phi[0][i], phi[1][i], phi[2][i], a); \
                /* Pairwise: -dV_pw/dphi_a where V_pw = (lambda/2) sum_{a<b} (phi_a-phi_b)^2 */ \
                double fpw = 0; \
                for (int b = 0; b < 3; b++) { \
                    if (b != a) fpw -= lambda * (phi[a][i] - phi[b][i]); \
                } \
                acc[a][i] = lapl - m2*phi[a][i] + fp + fpw; \
            } \
        } \
    } while(0)

    ACC1D();

    int print_every = Nt / 40;
    if (print_every < 1) print_every = 1;

    /* Track antisymmetric mode at measurement points */
    /* Measure at +/- 40 code lengths (well into far field) */
    double x_meas = 40.0;
    int i_plus  = (int)((x_meas + xmax) / dx);
    int i_minus = (int)((-x_meas + xmax) / dx);
    if (i_plus >= N1-1) i_plus = N1-2;
    if (i_minus < 1) i_minus = 1;

    /* Accumulate |A|^2 at measurement points over second half */
    double A2_plus = 0, A2_minus = 0;
    int n_meas = 0;
    int meas_start = Nt / 2;

    for (int n = 0; n <= Nt; n++) {
        double t = n * dt;

        if (n >= meas_start) {
            double A_p = (phi[0][i_plus] - phi[1][i_plus]) / M_SQRT2;
            double A_m = (phi[0][i_minus] - phi[1][i_minus]) / M_SQRT2;
            A2_plus  += A_p * A_p;
            A2_minus += A_m * A_m;
            n_meas++;
        }

        if (n % print_every == 0) {
            double S0 = (phi[0][ic] + phi[1][ic]) / M_SQRT2;
            double A0 = (phi[0][ic] - phi[1][ic]) / M_SQRT2;
            printf("  t=%7.1f  S(0)=%+.4f  A(0)=%+.4f  phi3(0)=%+.4f\n",
                   t, S0, A0, phi[2][ic]);
        }

        if (n == Nt) break;

        /* Velocity Verlet */
        for (int a = 0; a < 3; a++)
            for (int i = 1; i < N1-1; i++)
                vel[a][i] += 0.5 * dt * acc[a][i];
        for (int a = 0; a < 3; a++)
            for (int i = 1; i < N1-1; i++)
                phi[a][i] += dt * vel[a][i];
        ACC1D();
        for (int a = 0; a < 3; a++)
            for (int i = 1; i < N1-1; i++)
                vel[a][i] += 0.5 * dt * acc[a][i];

        /* absorbing */
        for (int a = 0; a < 3; a++)
            for (int i = 0; i < N1; i++) {
                vel[a][i] *= damp[i];
                phi[a][i] *= damp[i];
            }
    }

    A2_plus  /= n_meas;
    A2_minus /= n_meas;

    printf("\n=== Test 1 Results ===\n");
    printf("  <|A|^2> at x=+%.0f: %.6e\n", x_meas, A2_plus);
    printf("  <|A|^2> at x=-%.0f: %.6e\n", x_meas, A2_minus);
    double ratio = (A2_plus > 1e-30 && A2_minus > 1e-30) ?
                   A2_plus / A2_minus : 0.0;
    printf("  Forward/backward ratio: %.4f\n", ratio);
    printf("  Symmetric (spin-0)? %s (ratio should be ~1.0)\n",
           (ratio > 0.8 && ratio < 1.2) ? "YES" : "NO");

    /* Write results */
    char path[600];
    snprintf(path, sizeof(path), "%s/test1_symmetry.tsv", outdir);
    FILE *f = fopen(path, "w");
    if (f) {
        fprintf(f, "x_meas\tA2_plus\tA2_minus\tratio\n");
        fprintf(f, "%.1f\t%.6e\t%.6e\t%.4f\n", x_meas, A2_plus, A2_minus, ratio);
        fclose(f);
        printf("  Output: %s\n", path);
    }

    for (int a = 0; a < 3; a++) { free(phi[a]); free(vel[a]); free(acc[a]); }
    free(damp);
}

/* ====================================================================
 * TEST 2 & 3: 2D radiation pattern
 * ==================================================================== */

/* Simple LCG */
static unsigned int rng_state = 42;
static double randf(void)
{
    rng_state = rng_state * 1103515245 + 12345;
    return 2.0 * ((rng_state >> 16) & 0x7FFF) / 32767.0 - 1.0;
}

static void run_test2d(int test_num)
{
    double dx = 2.0 * L / (Nx - 1);
    double dy = 2.0 * L / (Ny - 1);
    double dx2 = dx * dx;
    double dy2 = dy * dy;
    double m2 = mass * mass;

    /* CFL for 2D */
    double cfl = 1.0 / sqrt(1.0/dx2 + 1.0/dy2 + m2);
    double dt = 0.4 * cfl;
    int Nt = (int)(tfinal / dt) + 1;

    int N = Nx * Ny;

    const char *test_desc = (test_num == 2) ?
        "2D pairwise only (lambda=0.99, eta=0)" :
        "2D pairwise + cross-gradient (lambda=0.99, eta=0.5)";

    printf("Test %d: %s\n", test_num, test_desc);
    printf("  mu=%.3f kappa=%.4f mass=%.4f lambda=%.4f eta=%.4f\n",
           mu, kappa, mass, lambda, eta);
    printf("  Nx=%d Ny=%d L=%.1f dx=%.5f dy=%.5f\n", Nx, Ny, L, dx, dy);
    printf("  dt=%.6f Nt=%d tfinal=%.0f\n", dt, Nt, tfinal);
    printf("  R_meas=%.1f N_ang=%d eps=%.4f\n", R_meas, N_ang, eps);
    printf("  Grid: %d points, Memory: %.1f MB\n",
           N, 9.0 * N * sizeof(double) / (1024.0*1024.0));

    /* Allocate: 3 fields x (phi, vel, acc) */
    double *phi[3], *vel[3], *acc[3];
    for (int a = 0; a < 3; a++) {
        phi[a] = calloc(N, sizeof(double));
        vel[a] = calloc(N, sizeof(double));
        acc[a] = calloc(N, sizeof(double));
        if (!phi[a] || !vel[a] || !acc[a]) {
            fprintf(stderr, "Allocation failed\n");
            exit(1);
        }
    }

    /* Circular absorbing boundary: damp outside r > 0.75*L */
    double *damp = malloc(N * sizeof(double));
    double r_abs = L * 0.75;
    for (int i = 0; i < Nx; i++) {
        double x = -L + i * dx;
        for (int j = 0; j < Ny; j++) {
            double y = -L + j * dy;
            double r = sqrt(x*x + y*y);
            if (r > r_abs) {
                double f = (r - r_abs) / (L - r_abs);
                if (f > 1.0) f = 1.0;
                damp[IDX(i,j)] = 1.0 - 0.98 * f * f;
            } else {
                damp[IDX(i,j)] = 1.0;
            }
        }
    }

    /* Initialize: Gaussian oscillon for all three fields + antisymmetric perturbation */
    int ic = Nx / 2, jc = Ny / 2;
    for (int i = 0; i < Nx; i++) {
        double x = -L + i * dx;
        for (int j = 0; j < Ny; j++) {
            double y = -L + j * dy;
            double r2 = x*x + y*y;
            double g = A_init * exp(-r2 / (2.0*sigma*sigma));

            /* Symmetric base + antisymmetric perturbation:
             * phi_1 = g + eps*g, phi_2 = g - eps*g
             * This gives A = (phi_1-phi_2)/sqrt(2) = sqrt(2)*eps*g */
            phi[0][IDX(i,j)] = g * (1.0 + eps);
            phi[1][IDX(i,j)] = g * (1.0 - eps);
            phi[2][IDX(i,j)] = g;
        }
    }

    /* Precompute cross-gradient scratch arrays if eta != 0 */
    double *cg_force[3];
    int use_cg = (fabs(eta) > 1e-15);
    if (use_cg) {
        for (int a = 0; a < 3; a++)
            cg_force[a] = calloc(N, sizeof(double));
    }

    /* Compute cross-gradient force contribution.
     *
     * Cross-gradient energy density (2D, 2 fields, field index a, spatial index i):
     *   E_cg = -(eta/2) * [ (dx phi_2)(dy phi_1) + (dy phi_1)(dx phi_2) ]
     *        = -eta * (dx phi_2)(dy phi_1)
     *
     * The cross-gradient couples field index a to spatial index a.
     * General form: E_cg = (eta/2) sum_{a!=b} (d_a phi_b)(d_b phi_a)
     *   = eta * (dx phi_2)(dy phi_1)    [2D, only a=x,b=y and a=y,b=x, same term]
     *
     * Force on phi_1 from cross-gradient: need d/dt(dL/d(dt phi_1)) and dL/d(phi_1).
     * Since E_cg depends on spatial derivatives only:
     *   F_1 = -delta E_cg / delta phi_1
     * E_cg = eta * (dx phi_2)(dy phi_1)
     * delta/delta phi_1: integrate by parts on dy:
     *   F_1 = -eta * dx(d/dy phi_2)... wait, need to be careful.
     *
     * E_cg = eta * integral (dx phi_2)(dy phi_1) dx dy
     * Variation w.r.t. phi_1:
     *   delta E_cg = eta * integral (dx phi_2)(dy delta phi_1) dx dy
     *              = -eta * integral dy(dx phi_2) * delta phi_1 dx dy   [integrate by parts in y]
     *   => F_1 = +eta * d_y(d_x phi_2) = eta * d_{xy} phi_2
     *
     * Variation w.r.t. phi_2:
     *   delta E_cg = eta * integral (dx delta phi_2)(dy phi_1) dx dy
     *              = -eta * integral dx(dy phi_1) * delta phi_2 dx dy   [integrate by parts in x]
     *   => F_2 = +eta * d_x(d_y phi_1) = eta * d_{xy} phi_1
     *
     * So: F_a = eta * d_{xy} phi_{other}, where other = 3-a for a=1,2.
     * phi_3 does not participate in cross-gradient (only 2 spatial dims, 2 field indices).
     */

    /* Compute acceleration macro for 2D */
    /* This is called with phi[], acc[] in scope */
    #define COMPUTE_ACC_2D() do { \
        /* Cross-gradient force if enabled */ \
        if (use_cg) { \
            for (int a = 0; a < 3; a++) \
                memset(cg_force[a], 0, N * sizeof(double)); \
            /* F_0 = eta * d_{xy} phi_1 */ \
            /* F_1 = eta * d_{xy} phi_0 */ \
            for (int i = 1; i < Nx-1; i++) \
                for (int j = 1; j < Ny-1; j++) { \
                    /* d_{xy} phi_b = d/dx(d/dy phi_b) approx by cross-stencil */ \
                    /* d_{xy} = [phi(i+1,j+1) - phi(i+1,j-1) - phi(i-1,j+1) + phi(i-1,j-1)] / (4 dx dy) */ \
                    double dxy_1 = (phi[1][IDX(i+1,j+1)] - phi[1][IDX(i+1,j-1)] \
                                  - phi[1][IDX(i-1,j+1)] + phi[1][IDX(i-1,j-1)]) / (4.0*dx*dy); \
                    double dxy_0 = (phi[0][IDX(i+1,j+1)] - phi[0][IDX(i+1,j-1)] \
                                  - phi[0][IDX(i-1,j+1)] + phi[0][IDX(i-1,j-1)]) / (4.0*dx*dy); \
                    cg_force[0][IDX(i,j)] = eta * dxy_1; \
                    cg_force[1][IDX(i,j)] = eta * dxy_0; \
                    /* phi_3 gets no cross-gradient force */ \
                } \
        } \
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
                    double fp = force_triple(phi[0][IDX(i,j)], phi[1][IDX(i,j)], phi[2][IDX(i,j)], a); \
                    /* Pairwise force: -lambda * sum_{b!=a} (phi_a - phi_b) */ \
                    double fpw = 0; \
                    for (int b = 0; b < 3; b++) { \
                        if (b != a) fpw -= lambda * (phi[a][IDX(i,j)] - phi[b][IDX(i,j)]); \
                    } \
                    double fcg = use_cg ? cg_force[a][IDX(i,j)] : 0.0; \
                    acc[a][IDX(i,j)] = lapl - m2*phi[a][IDX(i,j)] + fp + fpw + fcg; \
                } \
        } \
    } while(0)

    COMPUTE_ACC_2D();

    int print_every = Nt / 40;
    if (print_every < 1) print_every = 1;

    /* Angular measurement bins */
    double *ang_A2    = calloc(N_ang, sizeof(double));  /* accumulated |A|^2 */
    double *ang_S2    = calloc(N_ang, sizeof(double));  /* accumulated |S|^2 */
    double *ang_count = calloc(N_ang, sizeof(double));
    int n_meas = 0;
    int meas_start = Nt / 2;  /* measure second half */

    /* Precompute angular bin assignments for grid points near r=R_meas */
    /* Use all points within [R_meas - dr, R_meas + dr] */
    double dr_band = 2.0 * dx;  /* band width */
    int n_ring = 0;
    int *ring_idx = NULL;
    int *ring_ang = NULL;

    /* Count ring points first */
    for (int i = 1; i < Nx-1; i++) {
        double x = -L + i * dx;
        for (int j = 1; j < Ny-1; j++) {
            double y = -L + j * dy;
            double r = sqrt(x*x + y*y);
            if (fabs(r - R_meas) < dr_band) n_ring++;
        }
    }
    ring_idx = malloc(n_ring * sizeof(int));
    ring_ang = malloc(n_ring * sizeof(int));

    int ri = 0;
    for (int i = 1; i < Nx-1; i++) {
        double x = -L + i * dx;
        for (int j = 1; j < Ny-1; j++) {
            double y = -L + j * dy;
            double r = sqrt(x*x + y*y);
            if (fabs(r - R_meas) < dr_band) {
                double theta = atan2(y, x);  /* [-pi, pi] */
                int bin = (int)((theta + M_PI) / (2.0*M_PI) * N_ang);
                if (bin >= N_ang) bin = N_ang - 1;
                if (bin < 0) bin = 0;
                ring_idx[ri] = IDX(i,j);
                ring_ang[ri] = bin;
                ri++;
            }
        }
    }
    printf("  Ring points at r=%.1f: %d\n", R_meas, n_ring);

    /* Time evolution */
    for (int n = 0; n <= Nt; n++) {
        double t = n * dt;

        /* Accumulate angular pattern in second half */
        if (n >= meas_start) {
            for (int ri2 = 0; ri2 < n_ring; ri2++) {
                int k = ring_idx[ri2];
                int bin = ring_ang[ri2];
                double A = (phi[0][k] - phi[1][k]) / M_SQRT2;
                double S = (phi[0][k] + phi[1][k]) / M_SQRT2;
                ang_A2[bin] += A * A;
                ang_S2[bin] += S * S;
                ang_count[bin] += 1.0;
            }
            n_meas++;
        }

        if (n % print_every == 0) {
            int k0 = IDX(ic, jc);
            double S0 = (phi[0][k0] + phi[1][k0]) / M_SQRT2;
            double A0 = (phi[0][k0] - phi[1][k0]) / M_SQRT2;

            /* Compute energy */
            double Ek = 0, Eg = 0, Em = 0, Ep = 0, Epw = 0;
            /* Sample a fraction for speed */
            int step = 2;
            for (int i = 1; i < Nx-1; i += step)
                for (int j = 1; j < Ny-1; j += step) {
                    int k = IDX(i,j);
                    for (int a = 0; a < 3; a++) {
                        Ek += 0.5 * vel[a][k] * vel[a][k];
                        double dpx = (phi[a][IDX(i+1,j)] - phi[a][IDX(i-1,j)]) / (2.0*dx);
                        double dpy = (phi[a][IDX(i,j+1)] - phi[a][IDX(i,j-1)]) / (2.0*dy);
                        Eg += 0.5 * (dpx*dpx + dpy*dpy);
                        Em += 0.5 * m2 * phi[a][k] * phi[a][k];
                    }
                    double P = phi[0][k] * phi[1][k] * phi[2][k];
                    double P2 = P * P;
                    Ep += 0.5 * mu * P2 / (1.0 + kappa * P2);
                    for (int a = 0; a < 3; a++)
                        for (int b = a+1; b < 3; b++) {
                            double d = phi[a][k] - phi[b][k];
                            Epw += 0.5 * lambda * d * d;
                        }
                }
            double area = dx * dy * step * step;
            Ek *= area; Eg *= area; Em *= area; Ep *= area; Epw *= area;
            double Et = Ek + Eg + Em + Ep + Epw;

            printf("  t=%7.1f  S(0)=%+.4f  A(0)=%+.6f  phi3(0)=%+.4f  "
                   "E=%.2f (k=%.2f g=%.2f m=%.2f p=%.2f pw=%.2f)\n",
                   t, S0, A0, phi[2][k0], Et, Ek, Eg, Em, Ep, Epw);
        }

        if (n == Nt) break;

        /* Velocity Verlet */
        for (int a = 0; a < 3; a++)
            for (int k = 0; k < N; k++)
                vel[a][k] += 0.5 * dt * acc[a][k];
        for (int a = 0; a < 3; a++)
            for (int k = 0; k < N; k++)
                phi[a][k] += dt * vel[a][k];
        COMPUTE_ACC_2D();
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

    /* Normalize angular bins */
    for (int b = 0; b < N_ang; b++) {
        if (ang_count[b] > 0) {
            ang_A2[b] /= ang_count[b];
            ang_S2[b] /= ang_count[b];
        }
    }

    /* Write angular pattern */
    char path[600];
    snprintf(path, sizeof(path), "%s/test%d_angular.tsv", outdir, test_num);
    FILE *f = fopen(path, "w");
    if (f) {
        fprintf(f, "theta_deg\tA2\tS2\n");
        for (int b = 0; b < N_ang; b++) {
            double theta = -180.0 + (b + 0.5) * 360.0 / N_ang;
            fprintf(f, "%.1f\t%.6e\t%.6e\n", theta, ang_A2[b], ang_S2[b]);
        }
        fclose(f);
        printf("\n  Angular pattern written: %s\n", path);
    }

    /* Analyze: compute multipole moments of the angular pattern */
    /* C_l = (1/N) sum_k A2(theta_k) * e^{i*l*theta_k} */
    printf("\n=== Test %d Results: Angular Multipole Analysis ===\n", test_num);

    double C0_re = 0, C0_im = 0;  /* l=0 monopole (spin-0) */
    double C1_re = 0, C1_im = 0;  /* l=1 dipole (spin-1) */
    double C2_re = 0, C2_im = 0;  /* l=2 quadrupole (spin-2) */
    double C3_re = 0, C3_im = 0;  /* l=3 */
    double C4_re = 0, C4_im = 0;  /* l=4 */

    int n_valid = 0;
    for (int b = 0; b < N_ang; b++) {
        if (ang_count[b] > 0) {
            double theta = -M_PI + (b + 0.5) * 2.0*M_PI / N_ang;
            double A2 = ang_A2[b];
            C0_re += A2; C0_im += 0;
            C1_re += A2 * cos(theta);   C1_im += A2 * sin(theta);
            C2_re += A2 * cos(2*theta); C2_im += A2 * sin(2*theta);
            C3_re += A2 * cos(3*theta); C3_im += A2 * sin(3*theta);
            C4_re += A2 * cos(4*theta); C4_im += A2 * sin(4*theta);
            n_valid++;
        }
    }

    if (n_valid > 0) {
        C0_re /= n_valid;
        C1_re /= n_valid; C1_im /= n_valid;
        C2_re /= n_valid; C2_im /= n_valid;
        C3_re /= n_valid; C3_im /= n_valid;
        C4_re /= n_valid; C4_im /= n_valid;
    }

    double P0 = C0_re;  /* monopole power */
    double P1 = sqrt(C1_re*C1_re + C1_im*C1_im);
    double P2 = sqrt(C2_re*C2_re + C2_im*C2_im);
    double P3 = sqrt(C3_re*C3_re + C3_im*C3_im);
    double P4 = sqrt(C4_re*C4_re + C4_im*C4_im);

    printf("  Multipole amplitudes (|C_l|/C_0):\n");
    printf("    l=0 (monopole/spin-0): %.6e (reference)\n", P0);
    if (P0 > 1e-30) {
        printf("    l=1 (dipole/spin-1):   %.6e  ratio=%.4f\n", P1, P1/P0);
        printf("    l=2 (quadrupole/spin-2): %.6e  ratio=%.4f\n", P2, P2/P0);
        printf("    l=3 (octupole):          %.6e  ratio=%.4f\n", P3, P3/P0);
        printf("    l=4:                     %.6e  ratio=%.4f\n", P4, P4/P0);
    }

    /* Classification */
    double r1 = P0 > 1e-30 ? P1/P0 : 0;
    double r2 = P0 > 1e-30 ? P2/P0 : 0;
    const char *classification;
    if (r1 < 0.1 && r2 < 0.1)
        classification = "ISOTROPIC (spin-0 scalar)";
    else if (r1 > 0.2 && r1 > r2)
        classification = "DIPOLAR (spin-1 vector)";
    else if (r2 > 0.2 && r2 > r1)
        classification = "QUADRUPOLAR (spin-2 tensor)";
    else
        classification = "MIXED / AMBIGUOUS";

    printf("\n  CLASSIFICATION: %s\n", classification);

    /* Also measure forward/backward asymmetry as in 1D */
    /* +x vs -x */
    int b_plus  = (int)(0.5 * N_ang);  /* theta=0 => +x */
    int b_minus = 0;                    /* theta=-180 => -x */
    if (b_plus >= N_ang) b_plus = N_ang - 1;
    printf("  Forward (theta=0): %.6e  Backward (theta=180): %.6e  ratio=%.4f\n",
           ang_A2[b_plus], ang_A2[b_minus],
           (ang_A2[b_minus] > 1e-30) ? ang_A2[b_plus]/ang_A2[b_minus] : 0.0);

    /* Min/max of angular pattern */
    double a_min = 1e30, a_max = 0;
    int b_min_idx = 0, b_max_idx = 0;
    for (int b = 0; b < N_ang; b++) {
        if (ang_A2[b] < a_min && ang_count[b] > 0) { a_min = ang_A2[b]; b_min_idx = b; }
        if (ang_A2[b] > a_max && ang_count[b] > 0) { a_max = ang_A2[b]; b_max_idx = b; }
    }
    double aniso = (a_min > 1e-30) ? a_max / a_min : 0;
    printf("  Max/min ratio (anisotropy): %.4f (max at %.0f deg, min at %.0f deg)\n",
           aniso,
           -180.0 + (b_max_idx + 0.5) * 360.0 / N_ang,
           -180.0 + (b_min_idx + 0.5) * 360.0 / N_ang);

    /* Write multipole results */
    snprintf(path, sizeof(path), "%s/test%d_multipoles.tsv", outdir, test_num);
    f = fopen(path, "w");
    if (f) {
        fprintf(f, "l\tC_re\tC_im\t|C_l|\tratio\n");
        fprintf(f, "0\t%.6e\t0\t%.6e\t1.0000\n", P0, P0);
        fprintf(f, "1\t%.6e\t%.6e\t%.6e\t%.4f\n", C1_re, C1_im, P1, P0>1e-30?P1/P0:0);
        fprintf(f, "2\t%.6e\t%.6e\t%.6e\t%.4f\n", C2_re, C2_im, P2, P0>1e-30?P2/P0:0);
        fprintf(f, "3\t%.6e\t%.6e\t%.6e\t%.4f\n", C3_re, C3_im, P3, P0>1e-30?P3/P0:0);
        fprintf(f, "4\t%.6e\t%.6e\t%.6e\t%.4f\n", C4_re, C4_im, P4, P0>1e-30?P4/P0:0);
        fclose(f);
        printf("  Multipoles written: %s\n", path);
    }

    /* Cleanup */
    for (int a = 0; a < 3; a++) { free(phi[a]); free(vel[a]); free(acc[a]); }
    if (use_cg) for (int a = 0; a < 3; a++) free(cg_force[a]);
    free(damp); free(ang_A2); free(ang_S2); free(ang_count);
    free(ring_idx); free(ring_ang);
}

int main(int argc, char **argv)
{
    parse_args(argc, argv);

    printf("proca_spin: V24-S6 — Spin Classification of Proca Mediator\n\n");

    switch (test) {
        case 1:
            run_test1();
            break;
        case 2:
            eta = 0.0;  /* pairwise only */
            run_test2d(2);
            break;
        case 3:
            if (fabs(eta) < 1e-15) eta = 0.5;  /* default cross-gradient */
            run_test2d(3);
            break;
        default:
            fprintf(stderr, "Unknown test %d (use 1, 2, or 3)\n", test);
            return 1;
    }

    return 0;
}
