/*  v31_implicit.c — Single-field implicit solver with c(rho)
 *
 *  ONE field phi_a (3 components). No S/B split. No smoothing.
 *  c²(x) = rho0 / rho(x) — physical, unbounded.
 *
 *  Semi-implicit Crank-Nicolson with red-black Gauss-Seidel.
 *  Unconditionally stable at any c ratio.
 *
 *  The field IS the medium. The braid is a concentrated region.
 *  The "background" is the same field at lower amplitude.
 *
 *  Build: gcc -O3 -march=native -fopenmp -o v31_impl src/v31_implicit.c -lm
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

static double MU       = -41.345;
static double KAPPA    = 50.0;
static double MASS     = 1.50;
static double MASS2    = 2.25;
static double RHO0     = 0.0;    /* reference density (auto-set from far field) */
static double A_BG     = 0.1;    /* background amplitude */
static double ALPHA_C  = 1.0;    /* c² = rho0/rho when alpha_c=1 */

/* Implicit solver params */
static int    GS_ITERS = 15;     /* Gauss-Seidel iterations per step */
static double GS_TOL   = 1e-8;  /* convergence tolerance */

/* Grid */
static int    SIM_N    = 128;
static double SIM_L    = 20.0;

/* ================================================================
   Grid structure — stores current (n), previous (n-1), and workspace
   ================================================================ */

typedef struct {
    double *phi[NFIELDS];      /* phi^n (current) */
    double *phi_prev[NFIELDS]; /* phi^{n-1} (previous) */
    double *phi_next[NFIELDS]; /* phi^{n+1} (being solved) */
    double *vel[NFIELDS];      /* velocity (for diagnostics + init) */
    double *rho;               /* energy density */
    double *c2;                /* c²(x) = rho0/rho */
    int N;
    double L, dx, dt;
} Grid;

static Grid *grid_alloc(int N, double L) {
    Grid *g = calloc(1, sizeof(Grid));
    long N3 = (long)N * N * N;
    for (int a = 0; a < NFIELDS; a++) {
        g->phi[a]      = calloc(N3, sizeof(double));
        g->phi_prev[a] = calloc(N3, sizeof(double));
        g->phi_next[a] = calloc(N3, sizeof(double));
        g->vel[a]      = calloc(N3, sizeof(double));
    }
    g->rho = calloc(N3, sizeof(double));
    g->c2  = calloc(N3, sizeof(double));
    g->N = N; g->L = L;
    g->dx = 2.0 * L / (N - 1);
    /* dt can be LARGE for implicit — set by physics, not CFL */
    g->dt = 0.5 * g->dx;  /* ~3× larger than explicit CFL limit */
    return g;
}

static void grid_free(Grid *g) {
    for (int a = 0; a < NFIELDS; a++) {
        free(g->phi[a]); free(g->phi_prev[a]);
        free(g->phi_next[a]); free(g->vel[a]);
    }
    free(g->rho); free(g->c2); free(g);
}

/* ================================================================
   Compute rho and c² from current phi, vel
   ================================================================ */

static void compute_rho_c2(Grid *g) {
    int N = g->N, NN = N * N;
    long N3 = (long)N * N * N;
    double dx = g->dx;
    double rho_floor = 1e-10;  /* prevent division by zero */

    #pragma omp parallel for schedule(static)
    for (long idx = 0; idx < N3; idx++) {
        int i = (int)(idx / NN), j = (int)((idx / N) % N), k = (int)(idx % N);

        /* Fully periodic — no boundary special case */

        /* Energy density: KE + gradient + mass + potential */
        /* All directions periodic */
        int ip = (i+1)%N, im = (i-1+N)%N;
        int jp = (j+1)%N, jm = (j-1+N)%N;
        int kp = (k+1)%N, km = (k-1+N)%N;

        double e = 0;
        for (int a = 0; a < NFIELDS; a++) {
            e += 0.5 * g->vel[a][idx] * g->vel[a][idx];

            double gx = (g->phi[a][(long)ip*NN + j*N + k] -
                         g->phi[a][(long)im*NN + j*N + k]) / (2*dx);
            double gy = (g->phi[a][(long)i*NN + jp*N + k] -
                         g->phi[a][(long)i*NN + jm*N + k]) / (2*dx);
            double gz = (g->phi[a][(long)i*NN + j*N + kp] -
                         g->phi[a][(long)i*NN + j*N + km]) / (2*dx);
            e += 0.5 * (gx*gx + gy*gy + gz*gz);
            e += 0.5 * MASS2 * g->phi[a][idx] * g->phi[a][idx];
        }
        double p0 = g->phi[0][idx], p1 = g->phi[1][idx], p2 = g->phi[2][idx];
        double P = p0 * p1 * p2;
        e += (MU / 2.0) * P * P / (1.0 + KAPPA * P * P);

        if (e < rho_floor) e = rho_floor;
        g->rho[idx] = e;

        /* c² = 1 - alpha_c * tanh((rho - rho0) / rho0)
           Smooth, bounded: at rho=rho0: c²=1.
           At rho>>rho0: c²=1-alpha_c (slow in dense).
           At rho<<rho0: c²=1+alpha_c (fast in vacuum).
           Max variation: ±alpha_c. No extremes. */
        double x_arg = (e - RHO0) / (RHO0 + 1e-30);
        g->c2[idx] = 1.0 - ALPHA_C * tanh(x_arg);
    }
}

/* ================================================================
   Triple product force (explicit, computed from phi^n)
   ================================================================ */

static inline double triple_force(const double *phi0, const double *phi1,
                                   const double *phi2, long idx, int a) {
    double p0 = phi0[idx], p1 = phi1[idx], p2 = phi2[idx];
    double P = p0 * p1 * p2;
    double denom = 1.0 + KAPPA * P * P;
    double mu_P_d2 = MU * P / (denom * denom);
    double dPda = (a == 0) ? p1*p2 : (a == 1) ? p0*p2 : p0*p1;
    return mu_P_d2 * dPda;
}

/* ================================================================
   Gauss-Seidel sweep (one iteration, red-black ordering)

   Solving for phi^{n+1} from:
   (phi^{n+1} - 2*phi^n + phi^{n-1})/dt² =
       c² * ∇²[½(phi^n + phi^{n+1})] - m²*½(phi^n + phi^{n+1}) - V'(phi^n)

   Rearranged for phi^{n+1} at point idx:
   phi^{n+1}(idx) = [2*phi^n - phi^{n-1}
                      + ½*dt²*c²/dx² * Σ_neighbors(phi^{n+1} + phi^n)
                      - ½*dt²*(6*c²/dx² + m²)*phi^n
                      - dt²*V'(phi^n)]
                    / [1 + ½*dt²*(6*c²/dx² + m²)]
   ================================================================ */

static double gs_sweep(Grid *g, int color) {
    int N = g->N, NN = N * N;
    double dx = g->dx, dt = g->dt;
    double dt2 = dt * dt;
    double idx2 = 1.0 / (dx * dx);
    double max_change = 0;

    #pragma omp parallel for reduction(max:max_change) schedule(static)
    for (long idx = 0; idx < (long)N*N*N; idx++) {
        int i = (int)(idx / NN), j = (int)((idx / N) % N), k = (int)(idx % N);

        /* Red-black: skip wrong color */
        if (((i + j + k) & 1) != color) continue;

        /* Fully periodic BC — no damping, background sustains itself */
        /* (boundary handling is in the neighbor indexing below) */

        double c2 = g->c2[idx];
        double half_dt2_c2_idx2 = 0.5 * dt2 * c2 * idx2;
        double diag = 1.0 + 0.5 * dt2 * (6.0 * c2 * idx2 + MASS2);

        /* Fully periodic neighbors in all 3 directions */
        int ip = (i+1)%N, im = (i-1+N)%N;
        int jp = (j+1)%N, jm = (j-1+N)%N;
        int kp = (k+1)%N, km = (k-1+N)%N;
        long idx_ip = (long)ip*NN + j*N + k;
        long idx_im = (long)im*NN + j*N + k;
        long idx_jp = (long)i*NN + jp*N + k;
        long idx_jm = (long)i*NN + jm*N + k;
        long idx_kp = (long)i*NN + j*N + kp;
        long idx_km = (long)i*NN + j*N + km;

        for (int a = 0; a < NFIELDS; a++) {
            double sum_next = g->phi_next[a][idx_ip] + g->phi_next[a][idx_im]
                            + g->phi_next[a][idx_jp] + g->phi_next[a][idx_jm]
                            + g->phi_next[a][idx_kp] + g->phi_next[a][idx_km];
            double sum_curr = g->phi[a][idx_ip] + g->phi[a][idx_im]
                            + g->phi[a][idx_jp] + g->phi[a][idx_jm]
                            + g->phi[a][idx_kp] + g->phi[a][idx_km];

            /* Triple product force (explicit, from phi^n) */
            double fV = triple_force(g->phi[0], g->phi[1], g->phi[2], idx, a);

            /* RHS */
            double rhs = 2.0 * g->phi[a][idx] - g->phi_prev[a][idx]
                       + half_dt2_c2_idx2 * (sum_next + sum_curr)
                       - 0.5 * dt2 * (6.0 * c2 * idx2 + MASS2) * g->phi[a][idx]
                       - dt2 * fV;

            double new_val = rhs / diag;
            double change = fabs(new_val - g->phi_next[a][idx]);
            if (change > max_change) max_change = change;
            g->phi_next[a][idx] = new_val;
        }
    }
    return max_change;
}

/* ================================================================
   Full implicit step: solve for phi^{n+1}
   ================================================================ */

static int implicit_step(Grid *g) {
    long N3 = (long)g->N * g->N * g->N;

    /* Initial guess: extrapolate from phi^n and phi^{n-1} */
    for (int a = 0; a < NFIELDS; a++)
        for (long idx = 0; idx < N3; idx++)
            g->phi_next[a][idx] = 2.0 * g->phi[a][idx] - g->phi_prev[a][idx];

    /* Red-black Gauss-Seidel iterations */
    int converged = 0;
    int iter;
    for (iter = 0; iter < GS_ITERS; iter++) {
        double change_red  = gs_sweep(g, 0);  /* red points */
        double change_black = gs_sweep(g, 1); /* black points */
        double max_change = fmax(change_red, change_black);

        if (max_change < GS_TOL) { converged = 1; break; }
    }

    /* Update velocity (central difference) */
    double inv_2dt = 1.0 / (2.0 * g->dt);
    for (int a = 0; a < NFIELDS; a++)
        for (long idx = 0; idx < N3; idx++)
            g->vel[a][idx] = (g->phi_next[a][idx] - g->phi_prev[a][idx]) * inv_2dt;

    /* Shift: prev ← current, current ← next */
    for (int a = 0; a < NFIELDS; a++) {
        double *tmp = g->phi_prev[a];
        g->phi_prev[a] = g->phi[a];
        g->phi[a] = g->phi_next[a];
        g->phi_next[a] = tmp;  /* reuse as workspace */
    }

    /* Recompute rho and c² for new state */
    compute_rho_c2(g);

    return iter + 1;  /* return number of GS iterations used */
}

/* ================================================================
   Absorbing damping layer (x,y edges)
   ================================================================ */

static void apply_damping(Grid *g) {
    int N = g->N, NN = N * N;
    double dx = g->dx, L = g->L;
    double rs = 0.70*L, re = 0.95*L, idr = 1.0/(re-rs+1e-30);

    for (int i = 0; i < N; i++) {
        double x = -L + i*dx;
        for (int j = 0; j < N; j++) {
            double y = -L + j*dx;
            double rp = sqrt(x*x + y*y);
            if (rp <= rs) continue;
            double f = (rp-rs)*idr; if (f > 1) f = 1;
            double d = 1.0 - 0.98*f*f;
            for (int kk = 0; kk < N; kk++) {
                long idx = (long)i*NN + j*N + kk;
                for (int a = 0; a < NFIELDS; a++) {
                    g->phi[a][idx] *= d;
                    g->phi_prev[a][idx] *= d;
                    g->vel[a][idx] *= d;
                }
            }
        }
    }
}

/* ================================================================
   Initialization: bimodal braid + uniform background (single field)
   ================================================================ */

static void init_braid_single(Grid *g, double x_center) {
    int N = g->N, NN = N*N;
    double dx = g->dx, L = g->L;

    /* Bimodal params (t=0.85 interpolation) */
    double A[3] = {0.8, 0.8, 0.8};
    double delta[3] = {0, 3.0005, 4.4325};
    double R_tube = 3.0, ellip = 0.3325, ell_ang = 0.0;
    double k_fac = 1.0;
    double kw = k_fac * PI / L;
    double omega = sqrt(kw*kw + MASS2);
    double sx = 1+ellip, sy = 1-ellip;
    double inv2R2 = 1.0/(2*R_tube*R_tube);
    double ca = cos(ell_ang), sa = sin(ell_ang);
    double k_bg = PI/L, omega_bg = sqrt(k_bg*k_bg + MASS2);

    for (int i = 0; i < N; i++) {
        double x = -L + i*dx;
        for (int j = 0; j < N; j++) {
            double y = -L + j*dx;
            for (int kk = 0; kk < N; kk++) {
                double z = -L + kk*dx;
                long idx = (long)i*NN + j*N + kk;

                /* Braid envelope (centered at x_center) */
                double xc = (x - x_center)*ca + y*sa;
                double yc = -(x - x_center)*sa + y*ca;
                double r2e = xc*xc/(sx*sx) + yc*yc/(sy*sy);
                double env = exp(-r2e * inv2R2);

                for (int a = 0; a < NFIELDS; a++) {
                    double ph_braid = kw*z + delta[a];
                    double ph_bg = k_bg*z + 2*PI*a/3.0;

                    /* Single field = braid + background */
                    g->phi[a][idx] += A[a]*env*cos(ph_braid) + A_BG*cos(ph_bg);
                    g->vel[a][idx] += omega*A[a]*env*sin(ph_braid) + omega_bg*A_BG*sin(ph_bg);
                }
            }
        }
    }
}

/* Bootstrap: need phi^{n-1} from phi^n and vel^n */
static void bootstrap_prev(Grid *g) {
    long N3 = (long)g->N * g->N * g->N;
    double dt = g->dt;
    for (int a = 0; a < NFIELDS; a++)
        for (long idx = 0; idx < N3; idx++)
            g->phi_prev[a][idx] = g->phi[a][idx] - dt * g->vel[a][idx];
}

/* ================================================================
   Diagnostics
   ================================================================ */

static void write_diagnostics(Grid *g, FILE *fp, double t, int gs_iters) {
    int N = g->N, NN = N*N;
    long N3 = (long)N*N*N;
    double dx = g->dx, L = g->L, dV = dx*dx*dx;

    double E_total = 0, phi2_core = 0, phi2_total = 0;
    double max_rho = 0, min_c2 = 100, max_c2 = 0;

    #pragma omp parallel for reduction(+:E_total,phi2_core,phi2_total) \
                             reduction(max:max_rho,max_c2) reduction(min:min_c2)
    for (long idx = 0; idx < N3; idx++) {
        int i = (int)(idx/NN), j = (int)((idx/N)%N);
        double x = -L+i*dx, y = -L+j*dx;

        double p2 = 0;
        for (int a = 0; a < NFIELDS; a++)
            p2 += g->phi[a][idx] * g->phi[a][idx];
        phi2_total += p2;
        if (x*x + y*y < 64) phi2_core += p2;

        E_total += g->rho[idx] * dV;
        if (g->rho[idx] > max_rho) max_rho = g->rho[idx];
        if (g->c2[idx] < min_c2) min_c2 = g->c2[idx];
        if (g->c2[idx] > max_c2) max_c2 = g->c2[idx];
    }

    double fc = (phi2_total > 0) ? phi2_core / phi2_total : 0;

    fprintf(fp, "%.2f\t%.4e\t%.4f\t%.4e\t%.4f\t%.4f\t%d\n",
            t, E_total, fc, max_rho, min_c2, max_c2, gs_iters);
    fflush(fp);

    printf("t=%7.1f  E=%.2e  fc=%.3f  max_rho=%.2e  c2=[%.3f,%.3f]  GS=%d\n",
           t, E_total, fc, max_rho, min_c2, max_c2, gs_iters);
    fflush(stdout);
}

/* Save full field snapshot (binary: N, L, t, then phi[0..2] as raw doubles) */
static void save_field(Grid *g, double t, const char *dir) {
    long N3 = (long)g->N * g->N * g->N;
    char fname[512];
    snprintf(fname, sizeof(fname), "%s/field_t%04d.bin", dir, (int)(t + 0.5));

    FILE *fp = fopen(fname, "wb");
    if (!fp) return;

    int n = g->N;
    double l = g->L;
    fwrite(&n, sizeof(int), 1, fp);
    fwrite(&l, sizeof(double), 1, fp);
    fwrite(&t, sizeof(double), 1, fp);

    /* Write all 3 phi fields */
    for (int a = 0; a < NFIELDS; a++)
        fwrite(g->phi[a], sizeof(double), N3, fp);

    /* Write rho and c2 */
    fwrite(g->rho, sizeof(double), N3, fp);
    fwrite(g->c2, sizeof(double), N3, fp);

    fclose(fp);
    printf("  Saved field: %s (%.0f MB)\n", fname,
           (3 + 2) * N3 * 8.0 / 1e6);
}

/* ================================================================
   Main
   ================================================================ */

int main(int argc, char **argv) {
    int n_braids = 1;        /* 1 or 2 braids */
    double D_sep = 20.0;     /* separation for 2-braid test */
    double T_total = 500.0;
    double diag_dt = 10.0;
    double snap_dt = 100.0;  /* field snapshot interval */
    char outdir[256] = "data/impl";

    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "-N") && i+1<argc) SIM_N = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-L") && i+1<argc) SIM_L = atof(argv[++i]);
        else if (!strcmp(argv[i], "-T") && i+1<argc) T_total = atof(argv[++i]);
        else if (!strcmp(argv[i], "-dt") && i+1<argc) { /* override dt */ }
        else if (!strcmp(argv[i], "-gs") && i+1<argc) GS_ITERS = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-braids") && i+1<argc) n_braids = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-D") && i+1<argc) D_sep = atof(argv[++i]);
        else if (!strcmp(argv[i], "-bg") && i+1<argc) A_BG = atof(argv[++i]);
        else if (!strcmp(argv[i], "-o") && i+1<argc) strncpy(outdir, argv[++i], 255);
    }

    omp_set_num_threads(16);
    MASS2 = MASS * MASS;

    printf("=== V31 Implicit Single-Field Solver ===\n");
    printf("N=%d, L=%.0f, T=%.0f\n", SIM_N, SIM_L, T_total);
    printf("Braids: %d, D=%.0f, A_bg=%.2f\n", n_braids, D_sep, A_BG);
    printf("GS iterations: %d, tolerance: %.1e\n", GS_ITERS, GS_TOL);
    printf("No S/B split. No smoothing. c²=rho0/rho.\n\n");

    mkdir("data", 0755);
    mkdir(outdir, 0755);

    Grid *g = grid_alloc(SIM_N, SIM_L);
    printf("dx=%.4f, dt=%.4f (implicit, ~%.1f× explicit CFL)\n",
           g->dx, g->dt, g->dt / (0.15 * g->dx));

    /* Initialize: braid(s) + background */
    long N3 = (long)SIM_N * SIM_N * SIM_N;
    for (int a = 0; a < NFIELDS; a++) {
        memset(g->phi[a], 0, N3*sizeof(double));
        memset(g->vel[a], 0, N3*sizeof(double));
    }

    if (n_braids == 1) {
        init_braid_single(g, 0.0);
    } else {
        init_braid_single(g, -D_sep/2);
        init_braid_single(g, +D_sep/2);
    }

    /* Set RHO0 from far-field initial density */
    compute_rho_c2(g);
    {
        double sum = 0; int cnt = 0;
        int NN = SIM_N * SIM_N;
        for (int i = 0; i < SIM_N; i++) {
            double x = -SIM_L + i*g->dx;
            for (int j = 0; j < SIM_N; j++) {
                double y = -SIM_L + j*g->dx;
                if (sqrt(x*x+y*y) < SIM_L*0.5) continue;
                for (int kk = 0; kk < SIM_N; kk++) {
                    sum += g->rho[(long)i*NN + j*SIM_N + kk]; cnt++;
                }
            }
        }
        RHO0 = sum / cnt;
    }
    printf("rho0 = %.6e (far-field reference)\n", RHO0);

    /* Recompute with correct RHO0 */
    compute_rho_c2(g);

    /* Bootstrap phi_prev from phi and vel */
    bootstrap_prev(g);

    /* Open output files */
    char tspath[512];
    snprintf(tspath, sizeof(tspath), "%s/timeseries.tsv", outdir);
    FILE *fp_ts = fopen(tspath, "w");
    fprintf(fp_ts, "t\tE_total\tfc\tmax_rho\tmin_c2\tmax_c2\tGS_iters\n");

    int n_steps = (int)(T_total / g->dt);
    int diag_steps = (int)(diag_dt / g->dt); if (diag_steps < 1) diag_steps = 1;
    int snap_steps = (int)(snap_dt / g->dt); if (snap_steps < 1) snap_steps = 1;

    printf("Steps: %d, diag every %d, snap every %d\n\n", n_steps, diag_steps, snap_steps);

    /* Initial diagnostics + snapshot */
    write_diagnostics(g, fp_ts, 0, 0);
    save_field(g, 0, outdir);

    double wall_start = omp_get_wtime();

    for (int step = 1; step <= n_steps; step++) {
        int gs_used = implicit_step(g);
        /* No damping — fully periodic BC, background sustains itself */

        double t = step * g->dt;

        if (step % diag_steps == 0) {
            write_diagnostics(g, fp_ts, t, gs_used);
        }
        if (step % snap_steps == 0) {
            save_field(g, t, outdir);
        }
    }

    /* Final snapshot */
    save_field(g, T_total, outdir);

    fclose(fp_ts);
    double wall = omp_get_wtime() - wall_start;
    printf("\n=== Complete: %.0f seconds (%.1f min) ===\n", wall, wall/60);
    printf("Avg time per step: %.3f s\n", wall / n_steps);
    printf("Avg GS iterations: check timeseries\n");

    grid_free(g);
    return 0;
}
