/*
 * chain3d.c — Two-oscillon 3D stability test (V23-D Phase 2)
 *
 * Tests whether two oscillons at separation D remain bound in 3D,
 * where radiation pressure falls as 1/r² (geometric dilution).
 *
 * Phase 1: Equilibrate a single oscillon (N_eq grid, L_eq box, t_equil time)
 *          Save 3D profile at a breathing peak.
 * Phase 2: Place two copies at z = ±D/2 on a larger grid (N, L)
 *          Evolve for t_run, tracking z-centroids via energy weighting.
 *
 * Physics: triple-product coupling (v21 model)
 *   L = Σ_a [ ½(∂φ_a)² - ½m²φ_a² ] - V(P),  P = φ₁φ₂φ₃
 *   V = (μ/2) P² / (1 + κP²)
 *
 * Compile: gcc -O3 -fopenmp -Wall -o chain3d src/chain3d.c -lm
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

/* Parameters */
static double mu       = -20.0;
static double kappa    = 20.0;
static double mass     = 1.0;
static double A_init   = 0.8;
static double sigma    = 3.0;

/* Equilibration grid */
static int    N_eq     = 96;
static double L_eq     = 15.0;
static double t_equil  = 500.0;

/* Two-body grid */
static int    N        = 96;
static double L        = 25.0;
static double t_run    = 2000.0;
static double D_sep    = 16.0;

static double cfl_frac = 0.25;
static char   outdir[512] = "v23/phonon/data";

static void parse_args(int argc, char **argv)
{
    for (int i = 1; i < argc - 1; i += 2) {
        if      (!strcmp(argv[i], "-mu"))     mu       = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-kappa"))  kappa    = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-mass"))   mass     = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-A"))      A_init   = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-sigma"))  sigma    = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-Neq"))    N_eq     = atoi(argv[i+1]);
        else if (!strcmp(argv[i], "-Leq"))    L_eq     = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-teq"))    t_equil  = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-N"))      N        = atoi(argv[i+1]);
        else if (!strcmp(argv[i], "-L"))      L        = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-trun"))   t_run    = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-D"))      D_sep    = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-cfl"))    cfl_frac = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-o"))      strncpy(outdir, argv[i+1], sizeof(outdir)-1);
        else { fprintf(stderr, "Unknown flag: %s\n", argv[i]); exit(1); }
    }
}

/* ─── Index helpers ─── */
static inline long idx3(int i, int j, int k, int NN)
{
    return (long)i * NN * NN + (long)j * NN + (long)k;
}

/* ─── Allocate/free field arrays ─── */
typedef struct {
    double *phi[3], *vel[3], *acc[3];
    double *damp;
    int NN;
    double LL, dx, dx2, dt;
    long Ngrid;
} Grid;

static void grid_alloc(Grid *g, int NN, double LL)
{
    g->NN = NN;
    g->LL = LL;
    g->dx = 2.0 * LL / (NN - 1);
    g->dx2 = g->dx * g->dx;
    g->dt = cfl_frac * g->dx;
    g->Ngrid = (long)NN * NN * NN;
    for (int a = 0; a < 3; a++) {
        g->phi[a] = calloc(g->Ngrid, sizeof(double));
        g->vel[a] = calloc(g->Ngrid, sizeof(double));
        g->acc[a] = calloc(g->Ngrid, sizeof(double));
        if (!g->phi[a] || !g->vel[a] || !g->acc[a]) {
            fprintf(stderr, "Allocation failed\n"); exit(1);
        }
    }
    g->damp = malloc(g->Ngrid * sizeof(double));
    if (!g->damp) { fprintf(stderr, "Allocation failed\n"); exit(1); }

    double R_inner = LL * 0.70;
    double R_outer = LL * 0.95;
    double dx = g->dx;

    #pragma omp parallel for schedule(static)
    for (long idx = 0; idx < g->Ngrid; idx++) {
        int ii = idx / (NN * NN);
        int jj = (idx / NN) % NN;
        int kk = idx % NN;
        double x = -LL + ii * dx;
        double y = -LL + jj * dx;
        double z = -LL + kk * dx;
        double r = sqrt(x*x + y*y + z*z);
        if (r > R_inner) {
            double f = (r - R_inner) / (R_outer - R_inner);
            if (f > 1.0) f = 1.0;
            g->damp[idx] = 1.0 - 0.98 * f * f;
        } else {
            g->damp[idx] = 1.0;
        }
    }
}

static void grid_free(Grid *g)
{
    for (int a = 0; a < 3; a++) {
        free(g->phi[a]); free(g->vel[a]); free(g->acc[a]);
    }
    free(g->damp);
}

/* ─── Compute acceleration ─── */
static void compute_acc(Grid *g)
{
    int NN = g->NN;
    double dx2 = g->dx2;
    double m2 = mass * mass;

    for (int a = 0; a < 3; a++) {
        #pragma omp parallel for schedule(static)
        for (int i = 1; i < NN-1; i++) {
            for (int j = 1; j < NN-1; j++) {
                for (int k = 1; k < NN-1; k++) {
                    long idx = idx3(i,j,k,NN);

                    double lapl = (g->phi[a][idx3(i+1,j,k,NN)]
                                 + g->phi[a][idx3(i-1,j,k,NN)]
                                 + g->phi[a][idx3(i,j+1,k,NN)]
                                 + g->phi[a][idx3(i,j-1,k,NN)]
                                 + g->phi[a][idx3(i,j,k+1,NN)]
                                 + g->phi[a][idx3(i,j,k-1,NN)]
                                 - 6.0 * g->phi[a][idx]) / dx2;

                    double p0 = g->phi[0][idx];
                    double p1 = g->phi[1][idx];
                    double p2 = g->phi[2][idx];
                    double P = p0 * p1 * p2;
                    double P2 = P * P;
                    double denom2 = (1.0 + kappa * P2);
                    denom2 *= denom2;

                    double dP;
                    switch (a) {
                        case 0: dP = p1 * p2; break;
                        case 1: dP = p0 * p2; break;
                        default: dP = p0 * p1; break;
                    }
                    double dVdphi = mu * P * dP / denom2;

                    g->acc[a][idx] = lapl - m2 * g->phi[a][idx] - dVdphi;
                }
            }
        }

        /* Boundary: zero acceleration */
        #pragma omp parallel for schedule(static)
        for (long idx = 0; idx < g->Ngrid; idx++) {
            int ii = idx / (NN * NN);
            int jj = (idx / NN) % NN;
            int kk = idx % NN;
            if (ii == 0 || ii == NN-1 || jj == 0 || jj == NN-1 || kk == 0 || kk == NN-1)
                g->acc[a][idx] = 0.0;
        }
    }
}

/* ─── One Verlet step ─── */
static void verlet_step(Grid *g)
{
    long Ngrid = g->Ngrid;
    double dt = g->dt;

    /* Half-kick */
    for (int a = 0; a < 3; a++) {
        #pragma omp parallel for schedule(static)
        for (long idx = 0; idx < Ngrid; idx++)
            g->vel[a][idx] += 0.5 * dt * g->acc[a][idx];
    }

    /* Drift */
    for (int a = 0; a < 3; a++) {
        #pragma omp parallel for schedule(static)
        for (long idx = 0; idx < Ngrid; idx++)
            g->phi[a][idx] += dt * g->vel[a][idx];
    }

    /* Recompute acceleration */
    compute_acc(g);

    /* Half-kick */
    for (int a = 0; a < 3; a++) {
        #pragma omp parallel for schedule(static)
        for (long idx = 0; idx < Ngrid; idx++)
            g->vel[a][idx] += 0.5 * dt * g->acc[a][idx];
    }

    /* Absorbing boundary */
    for (int a = 0; a < 3; a++) {
        #pragma omp parallel for schedule(static)
        for (long idx = 0; idx < Ngrid; idx++) {
            g->vel[a][idx] *= g->damp[idx];
            g->phi[a][idx] *= g->damp[idx];
        }
    }
}

/* ─── Compute energy density at a grid point ─── */
static inline double energy_density(Grid *g, int i, int j, int k)
{
    int NN = g->NN;
    double dx = g->dx;
    double m2 = mass * mass;
    long idx = idx3(i,j,k,NN);
    double e = 0;
    for (int a = 0; a < 3; a++) {
        e += 0.5 * g->vel[a][idx] * g->vel[a][idx];
        double gx = (g->phi[a][idx3(i+1,j,k,NN)] - g->phi[a][idx3(i-1,j,k,NN)]) / (2*dx);
        double gy = (g->phi[a][idx3(i,j+1,k,NN)] - g->phi[a][idx3(i,j-1,k,NN)]) / (2*dx);
        double gz = (g->phi[a][idx3(i,j,k+1,NN)] - g->phi[a][idx3(i,j,k-1,NN)]) / (2*dx);
        e += 0.5 * (gx*gx + gy*gy + gz*gz);
        e += 0.5 * m2 * g->phi[a][idx] * g->phi[a][idx];
    }
    double P = g->phi[0][idx] * g->phi[1][idx] * g->phi[2][idx];
    double P2 = P * P;
    e += 0.5 * mu * P2 / (1.0 + kappa * P2);
    return e;
}

/* ─── Compute total energy ─── */
static double total_energy(Grid *g)
{
    int NN = g->NN;
    double dV = g->dx * g->dx * g->dx;
    double Etot = 0;

    #pragma omp parallel for reduction(+:Etot) schedule(static)
    for (int i = 1; i < NN-1; i++) {
        for (int j = 1; j < NN-1; j++) {
            for (int k = 1; k < NN-1; k++) {
                Etot += energy_density(g, i, j, k) * dV;
            }
        }
    }
    return Etot;
}

/* ─── Find peak |phi[0]| (proxy for oscillon amplitude) ─── */
static double peak_amplitude(Grid *g)
{
    double pk = 0;
    #pragma omp parallel for reduction(max:pk) schedule(static)
    for (long idx = 0; idx < g->Ngrid; idx++) {
        double a = fabs(g->phi[0][idx]);
        if (a > pk) pk = a;
    }
    return pk;
}

int main(int argc, char **argv)
{
    parse_args(argc, argv);

    double m2 = mass * mass;
    (void)m2;

    printf("=== chain3d: Two-oscillon 3D stability test ===\n");
    printf("  mu=%.1f kappa=%.1f mass=%.3f A=%.3f sigma=%.3f\n",
           mu, kappa, mass, A_init, sigma);
    printf("  D_sep=%.1f  t_equil=%.0f  t_run=%.0f\n", D_sep, t_equil, t_run);
    printf("  Threads: %d\n", omp_get_max_threads());
    fflush(stdout);

    /* ═══════════════════════════════════════════════
     * PHASE 1: Equilibrate single oscillon
     * ═══════════════════════════════════════════════ */
    printf("\n--- Phase 1: Equilibrate single oscillon ---\n");
    printf("  N_eq=%d  L_eq=%.1f  dx=%.4f\n", N_eq, L_eq, 2.0*L_eq/(N_eq-1));
    fflush(stdout);

    Grid geq;
    grid_alloc(&geq, N_eq, L_eq);

    printf("  Memory: %.1f MB per field, %.1f MB total\n",
           geq.Ngrid*8.0/1e6, geq.Ngrid*8.0*10/1e6);
    fflush(stdout);

    /* Initialize: spherical Gaussian */
    double dx_eq = geq.dx;
    #pragma omp parallel for schedule(static)
    for (long idx = 0; idx < geq.Ngrid; idx++) {
        int i = idx / (N_eq * N_eq);
        int j = (idx / N_eq) % N_eq;
        int k = idx % N_eq;
        double x = -L_eq + i * dx_eq;
        double y = -L_eq + j * dx_eq;
        double z = -L_eq + k * dx_eq;
        double r2 = x*x + y*y + z*z;
        double g = A_init * exp(-r2 / (2.0 * sigma * sigma));
        for (int a = 0; a < 3; a++)
            geq.phi[a][idx] = g;
    }

    compute_acc(&geq);

    int Nt_eq = (int)(t_equil / geq.dt) + 1;
    int print_eq = Nt_eq / 50;
    if (print_eq < 1) print_eq = 1;

    double wall_start = omp_get_wtime();

    /* Track center value to find breathing peak */
    long ic_eq = idx3(N_eq/2, N_eq/2, N_eq/2, N_eq);
    double prev_phi0 = geq.phi[0][ic_eq];
    double prev_prev_phi0 = prev_phi0;
    int last_peak_step = 0;
    double last_peak_val = 0;

    /* We'll save the profile at the LAST breathing peak before t_equil */
    /* Use a buffer for the snapshot */
    double *snap_phi[3];
    double *snap_vel[3];
    for (int a = 0; a < 3; a++) {
        snap_phi[a] = malloc(geq.Ngrid * sizeof(double));
        snap_vel[a] = malloc(geq.Ngrid * sizeof(double));
    }

    for (int n = 0; n < Nt_eq; n++) {
        verlet_step(&geq);
        double t = (n+1) * geq.dt;

        double cur = geq.phi[0][ic_eq];

        /* Detect breathing peak: prev > prev_prev AND prev > cur */
        if (prev_phi0 > prev_prev_phi0 && prev_phi0 > cur && prev_phi0 > 0.01) {
            last_peak_step = n;
            last_peak_val = prev_phi0;

            /* Save snapshot */
            for (int a = 0; a < 3; a++) {
                memcpy(snap_phi[a], geq.phi[a], geq.Ngrid * sizeof(double));
                memcpy(snap_vel[a], geq.vel[a], geq.Ngrid * sizeof(double));
            }
        }

        prev_prev_phi0 = prev_phi0;
        prev_phi0 = cur;

        if ((n+1) % print_eq == 0) {
            double E = total_energy(&geq);
            double pk = peak_amplitude(&geq);
            double wall = omp_get_wtime() - wall_start;
            double frac = (double)(n+1) / Nt_eq;
            double eta = (frac > 0.01) ? wall * (1.0-frac)/frac : 0;
            printf("  t=%7.1f  phi0=%.4f  pk=%.4f  E=%.2f  [%.0fs, ETA %.0fs]\n",
                   t, cur, pk, E, wall, eta);
            fflush(stdout);
        }
    }

    double wall_eq = omp_get_wtime() - wall_start;
    printf("  Phase 1 done: %.1f s. Last peak at step %d (val=%.4f)\n",
           wall_eq, last_peak_step, last_peak_val);
    fflush(stdout);

    if (last_peak_step == 0) {
        fprintf(stderr, "ERROR: No breathing peak detected during equilibration!\n");
        return 1;
    }

    /* ═══════════════════════════════════════════════
     * PHASE 2: Place two oscillons, evolve, track
     * ═══════════════════════════════════════════════ */
    printf("\n--- Phase 2: Two oscillons at D=%.1f ---\n", D_sep);
    printf("  N=%d  L=%.1f  dx=%.4f\n", N, L, 2.0*L/(N-1));
    fflush(stdout);

    Grid g2;
    grid_alloc(&g2, N, L);

    printf("  Memory: %.1f MB per field, %.1f MB total\n",
           g2.Ngrid*8.0/1e6, g2.Ngrid*8.0*10/1e6);
    fflush(stdout);

    /* Interpolate equilibrated profile onto two-body grid at z = ±D/2.
     * The equilibrated profile is on grid [-L_eq, L_eq] with spacing dx_eq.
     * We use trilinear interpolation. */
    double z_upper_init = +D_sep / 2.0;
    double z_lower_init = -D_sep / 2.0;
    double dx2_grid = g2.dx;

    printf("  Placing oscillons at z = %.1f and z = %.1f\n",
           z_upper_init, z_lower_init);
    printf("  Interpolating from equilibration grid (N=%d,L=%.1f) to production grid (N=%d,L=%.1f)\n",
           N_eq, L_eq, N, L);
    fflush(stdout);

    #pragma omp parallel for schedule(static)
    for (long idx = 0; idx < g2.Ngrid; idx++) {
        int i = idx / (N * N);
        int j = (idx / N) % N;
        int k = idx % N;
        double x = -L + i * dx2_grid;
        double y = -L + j * dx2_grid;
        double z = -L + k * dx2_grid;

        double val_phi[3] = {0, 0, 0};
        double val_vel[3] = {0, 0, 0};

        /* Contribution from upper oscillon (centered at z_upper_init) */
        double xu = x, yu = y, zu = z - z_upper_init;
        /* Map to equilibration grid indices */
        double fi_u = (xu + L_eq) / dx_eq;
        double fj_u = (yu + L_eq) / dx_eq;
        double fk_u = (zu + L_eq) / dx_eq;
        int ii_u = (int)floor(fi_u); int jj_u = (int)floor(fj_u); int kk_u = (int)floor(fk_u);

        if (ii_u >= 0 && ii_u < N_eq-1 && jj_u >= 0 && jj_u < N_eq-1 &&
            kk_u >= 0 && kk_u < N_eq-1) {
            double fx = fi_u - ii_u, fy = fj_u - jj_u, fz = fk_u - kk_u;
            for (int a = 0; a < 3; a++) {
                /* Trilinear interpolation */
                double c000 = snap_phi[a][idx3(ii_u,  jj_u,  kk_u,  N_eq)];
                double c100 = snap_phi[a][idx3(ii_u+1,jj_u,  kk_u,  N_eq)];
                double c010 = snap_phi[a][idx3(ii_u,  jj_u+1,kk_u,  N_eq)];
                double c110 = snap_phi[a][idx3(ii_u+1,jj_u+1,kk_u,  N_eq)];
                double c001 = snap_phi[a][idx3(ii_u,  jj_u,  kk_u+1,N_eq)];
                double c101 = snap_phi[a][idx3(ii_u+1,jj_u,  kk_u+1,N_eq)];
                double c011 = snap_phi[a][idx3(ii_u,  jj_u+1,kk_u+1,N_eq)];
                double c111 = snap_phi[a][idx3(ii_u+1,jj_u+1,kk_u+1,N_eq)];
                val_phi[a] += c000*(1-fx)*(1-fy)*(1-fz)
                            + c100*fx*(1-fy)*(1-fz)
                            + c010*(1-fx)*fy*(1-fz)
                            + c110*fx*fy*(1-fz)
                            + c001*(1-fx)*(1-fy)*fz
                            + c101*fx*(1-fy)*fz
                            + c011*(1-fx)*fy*fz
                            + c111*fx*fy*fz;

                double v000 = snap_vel[a][idx3(ii_u,  jj_u,  kk_u,  N_eq)];
                double v100 = snap_vel[a][idx3(ii_u+1,jj_u,  kk_u,  N_eq)];
                double v010 = snap_vel[a][idx3(ii_u,  jj_u+1,kk_u,  N_eq)];
                double v110 = snap_vel[a][idx3(ii_u+1,jj_u+1,kk_u,  N_eq)];
                double v001 = snap_vel[a][idx3(ii_u,  jj_u,  kk_u+1,N_eq)];
                double v101 = snap_vel[a][idx3(ii_u+1,jj_u,  kk_u+1,N_eq)];
                double v011 = snap_vel[a][idx3(ii_u,  jj_u+1,kk_u+1,N_eq)];
                double v111 = snap_vel[a][idx3(ii_u+1,jj_u+1,kk_u+1,N_eq)];
                val_vel[a] += v000*(1-fx)*(1-fy)*(1-fz)
                            + v100*fx*(1-fy)*(1-fz)
                            + v010*(1-fx)*fy*(1-fz)
                            + v110*fx*fy*(1-fz)
                            + v001*(1-fx)*(1-fy)*fz
                            + v101*fx*(1-fy)*fz
                            + v011*(1-fx)*fy*fz
                            + v111*fx*fy*fz;
            }
        }

        /* Contribution from lower oscillon (centered at z_lower_init) */
        double xl = x, yl = y, zl = z - z_lower_init;
        double fi_l = (xl + L_eq) / dx_eq;
        double fj_l = (yl + L_eq) / dx_eq;
        double fk_l = (zl + L_eq) / dx_eq;
        int ii_l = (int)floor(fi_l); int jj_l = (int)floor(fj_l); int kk_l = (int)floor(fk_l);

        if (ii_l >= 0 && ii_l < N_eq-1 && jj_l >= 0 && jj_l < N_eq-1 &&
            kk_l >= 0 && kk_l < N_eq-1) {
            double fx = fi_l - ii_l, fy = fj_l - jj_l, fz = fk_l - kk_l;
            for (int a = 0; a < 3; a++) {
                double c000 = snap_phi[a][idx3(ii_l,  jj_l,  kk_l,  N_eq)];
                double c100 = snap_phi[a][idx3(ii_l+1,jj_l,  kk_l,  N_eq)];
                double c010 = snap_phi[a][idx3(ii_l,  jj_l+1,kk_l,  N_eq)];
                double c110 = snap_phi[a][idx3(ii_l+1,jj_l+1,kk_l,  N_eq)];
                double c001 = snap_phi[a][idx3(ii_l,  jj_l,  kk_l+1,N_eq)];
                double c101 = snap_phi[a][idx3(ii_l+1,jj_l,  kk_l+1,N_eq)];
                double c011 = snap_phi[a][idx3(ii_l,  jj_l+1,kk_l+1,N_eq)];
                double c111 = snap_phi[a][idx3(ii_l+1,jj_l+1,kk_l+1,N_eq)];
                val_phi[a] += c000*(1-fx)*(1-fy)*(1-fz)
                            + c100*fx*(1-fy)*(1-fz)
                            + c010*(1-fx)*fy*(1-fz)
                            + c110*fx*fy*(1-fz)
                            + c001*(1-fx)*(1-fy)*fz
                            + c101*fx*(1-fy)*fz
                            + c011*(1-fx)*fy*fz
                            + c111*fx*fy*fz;

                double v000 = snap_vel[a][idx3(ii_l,  jj_l,  kk_l,  N_eq)];
                double v100 = snap_vel[a][idx3(ii_l+1,jj_l,  kk_l,  N_eq)];
                double v010 = snap_vel[a][idx3(ii_l,  jj_l+1,kk_l,  N_eq)];
                double v110 = snap_vel[a][idx3(ii_l+1,jj_l+1,kk_l,  N_eq)];
                double v001 = snap_vel[a][idx3(ii_l,  jj_l,  kk_l+1,N_eq)];
                double v101 = snap_vel[a][idx3(ii_l+1,jj_l,  kk_l+1,N_eq)];
                double v011 = snap_vel[a][idx3(ii_l,  jj_l+1,kk_l+1,N_eq)];
                double v111 = snap_vel[a][idx3(ii_l+1,jj_l+1,kk_l+1,N_eq)];
                val_vel[a] += v000*(1-fx)*(1-fy)*(1-fz)
                            + v100*fx*(1-fy)*(1-fz)
                            + v010*(1-fx)*fy*(1-fz)
                            + v110*fx*fy*(1-fz)
                            + v001*(1-fx)*(1-fy)*fz
                            + v101*fx*(1-fy)*fz
                            + v011*(1-fx)*fy*fz
                            + v111*fx*fy*fz;
            }
        }

        for (int a = 0; a < 3; a++) {
            g2.phi[a][idx] = val_phi[a];
            g2.vel[a][idx] = val_vel[a];
        }
    }

    /* Free equilibration grid and snapshots */
    grid_free(&geq);
    for (int a = 0; a < 3; a++) {
        free(snap_phi[a]);
        free(snap_vel[a]);
    }

    printf("  Initialization done. Computing initial energy...\n");
    fflush(stdout);

    compute_acc(&g2);

    double E0 = total_energy(&g2);
    printf("  Initial E_total = %.4f\n", E0);
    fflush(stdout);

    /* Output file */
    char tspath[600];
    snprintf(tspath, sizeof(tspath), "%s/chain3d_D%02d_ts.tsv",
             outdir, (int)D_sep);
    FILE *fts = fopen(tspath, "w");
    if (!fts) { fprintf(stderr, "Cannot open %s\n", tspath); return 1; }
    fprintf(fts, "time\tseparation\tE_total\tz_upper\tz_lower\t"
                 "pk_upper\tpk_lower\tE_upper\tE_lower\n");

    int Nt_run = (int)(t_run / g2.dt) + 1;
    int diag_every = Nt_run / 5000;
    if (diag_every < 1) diag_every = 1;
    int print_every = Nt_run / 100;
    if (print_every < 1) print_every = 1;

    double wall_p2 = omp_get_wtime();

    for (int n = 0; n <= Nt_run; n++) {
        double t = n * g2.dt;

        int do_diag  = (n % diag_every == 0);
        int do_print = (n % print_every == 0);

        if (do_diag || do_print) {
            /* Compute energy-weighted centroids, split at z=0 */
            double wz_up = 0, w_up = 0;
            double wz_lo = 0, w_lo = 0;
            double E_tot = 0, E_up = 0, E_lo = 0;
            double pk_up = 0, pk_lo = 0;
            double dV = g2.dx * g2.dx * g2.dx;

            #pragma omp parallel
            {
                double lwzu=0, lwu=0, lwzl=0, lwl=0;
                double lEt=0, lEu=0, lEl=0;
                double lpku=0, lpkl=0;

                #pragma omp for schedule(static) nowait
                for (int i = 1; i < N-1; i++) {
                    double x = -L + i * dx2_grid;
                    for (int j = 1; j < N-1; j++) {
                        double y = -L + j * dx2_grid;
                        for (int k = 1; k < N-1; k++) {
                            double z = -L + k * dx2_grid;
                            double e = energy_density(&g2, i, j, k);
                            double edv = e * dV;
                            lEt += edv;

                            if (z > 0) {
                                lEu += edv;
                                lwzu += z * edv;
                                lwu += edv;
                                double ap = fabs(g2.phi[0][idx3(i,j,k,N)]);
                                if (ap > lpku) lpku = ap;
                            } else {
                                lEl += edv;
                                lwzl += z * edv;
                                lwl += edv;
                                double ap = fabs(g2.phi[0][idx3(i,j,k,N)]);
                                if (ap > lpkl) lpkl = ap;
                            }
                        }
                    }
                }

                #pragma omp critical
                {
                    wz_up += lwzu; w_up += lwu;
                    wz_lo += lwzl; w_lo += lwl;
                    E_tot += lEt; E_up += lEu; E_lo += lEl;
                    if (lpku > pk_up) pk_up = lpku;
                    if (lpkl > pk_lo) pk_lo = lpkl;
                }
            }

            double z_up = (w_up > 1e-20) ? wz_up / w_up : z_upper_init;
            double z_lo = (w_lo > 1e-20) ? wz_lo / w_lo : z_lower_init;
            double sep = z_up - z_lo;

            if (do_diag)
                fprintf(fts, "%.4f\t%.6f\t%.4f\t%.6f\t%.6f\t%.4f\t%.4f\t%.4f\t%.4f\n",
                        t, sep, E_tot, z_up, z_lo, pk_up, pk_lo, E_up, E_lo);

            if (do_print) {
                double wall = omp_get_wtime() - wall_p2;
                double frac = (double)n / Nt_run;
                double eta = (frac > 0.01) ? wall * (1.0-frac)/frac : 0;
                printf("  t=%7.1f  sep=%.4f  z+=%+.3f z-=%+.3f  "
                       "E=%.2f  pk=%.3f/%.3f  [%.0fs, ETA %.0fs]\n",
                       t, sep, z_up, z_lo, E_tot, pk_up, pk_lo, wall, eta);
                fflush(stdout);
            }
        }

        if (n == Nt_run) break;

        verlet_step(&g2);
    }

    fclose(fts);

    double wall_total = omp_get_wtime() - wall_start;
    printf("\n=== SUMMARY ===\n");
    printf("  D_sep=%.1f  t_run=%.0f\n", D_sep, t_run);
    printf("  Wall time: %.1f s (%.1f min)\n", wall_total, wall_total/60);
    printf("  Output: %s\n", tspath);
    printf("===\n");

    grid_free(&g2);
    return 0;
}
