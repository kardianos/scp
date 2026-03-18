/*  v30_cfield.c — Field-dependent speed of light
 *
 *  The field IS the medium. c depends on local field density:
 *     c_eff²(x) = c₀² × (ρ(x)/ρ_ref)^alpha
 *
 *  Dense field → faster propagation. Depleted field → slower propagation.
 *  Waves curve TOWARD depleted regions (gravitational lensing analog).
 *
 *  Setup:
 *  - LARGE domain (L=5000), N=512 (dx≈20)
 *  - Concentrated field blob at center (R_init=200)
 *  - Pinned boundaries (phi=0 at edges = empty space)
 *  - NO FRW expansion — the field expands naturally through its own dynamics
 *  - Asymmetric expansion from nonlinear c(ρ) creates structure
 *
 *  EOM: ∂²φ_a/∂t² = c_eff²(x) × ∇²φ_a - m²φ_a - ∂V/∂φ_a
 *
 *  Build: gcc -O3 -fopenmp -o v30_cfield src/v30_cfield.c -lm
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

static int    SIM_N      = 512;
static double SIM_L      = 5000.0;
static double A0_INIT    = 1.0;     /* initial blob amplitude (near quartic vacuum v=1) */
static double R_INIT     = 200.0;   /* initial blob radius */
static double C_ALPHA    = 0.5;     /* c_eff = (rho/rho_ref)^alpha */
static double RHO_REF    = 1.0;     /* reference density for c_eff=1 */
static double MU         = -41.345;
static double KAPPA      = 50.0;
static double MASS       = 0.0;   /* massless at large scale — mass emerges from c(rho) */
static double T_TOTAL    = 2000.0;
static double DIAG_DT    = 50.0;
static int    SMOOTH_R   = 2;       /* smoothing radius for rho (grid cells) */
static char   OUTDIR[256] = "data/cfield";

/* ================================================================
   Grid
   ================================================================ */

typedef struct {
    double *phi[NFIELDS];
    double *vel[NFIELDS];
    double *acc[NFIELDS];
    double *rho;       /* local energy density */
    double *c2_eff;    /* c_eff² at each point */
    int N;
    double L, dx, dt;
} Grid;

static Grid *grid_alloc(int N, double L) {
    Grid *g = calloc(1, sizeof(Grid));
    long N3 = (long)N * N * N;
    for (int a = 0; a < NFIELDS; a++) {
        g->phi[a] = calloc(N3, sizeof(double));
        g->vel[a] = calloc(N3, sizeof(double));
        g->acc[a] = calloc(N3, sizeof(double));
    }
    g->rho    = calloc(N3, sizeof(double));
    g->c2_eff = calloc(N3, sizeof(double));
    g->N = N; g->L = L;
    g->dx = 2.0 * L / (N - 1);
    /* CFL: need dt << min(dx/c_max, 1/omega_max).
       c_max=1 (c_eff≤1 in dense regions), omega=sqrt(k²+m²)≈m=1.5 for long waves.
       Period = 2π/1.5 ≈ 4.2. Need dt << 1 for stability. */
    g->dt = fmin(0.15 * g->dx, 0.5);  /* cap at 0.5 regardless of dx */
    return g;
}

static void grid_free(Grid *g) {
    for (int a = 0; a < NFIELDS; a++) {
        free(g->phi[a]); free(g->vel[a]); free(g->acc[a]);
    }
    free(g->rho); free(g->c2_eff);
    free(g);
}

/* ================================================================
   Compute local energy density and c_eff²
   ================================================================ */

static void compute_rho_and_ceff(Grid *g) {
    int N = g->N, NN = N*N;
    long N3 = (long)N*N*N;
    double dx = g->dx;
    double mass2 = MASS * MASS;

    /* Step 1: raw energy density */
    #pragma omp parallel for schedule(static)
    for (long idx = 0; idx < N3; idx++) {
        int i = idx / NN, j = (idx / N) % N, k = idx % N;

        /* Boundary: rho=0 at edges */
        if (i < 1 || i >= N-1 || j < 1 || j >= N-1 || k < 1 || k >= N-1) {
            g->rho[idx] = 0;
            continue;
        }

        double e = 0;
        for (int a = 0; a < NFIELDS; a++) {
            e += 0.5 * g->vel[a][idx] * g->vel[a][idx];
            double gx = (g->phi[a][idx+NN] - g->phi[a][idx-NN]) / (2*dx);
            double gy = (g->phi[a][idx+N]  - g->phi[a][idx-N])  / (2*dx);
            double gz = (g->phi[a][idx+1]  - g->phi[a][idx-1])  / (2*dx);
            e += 0.5 * (gx*gx + gy*gy + gz*gz);
            e += 0.5 * mass2 * g->phi[a][idx] * g->phi[a][idx];
        }
        g->rho[idx] = e;
    }

    /* Step 2: smooth rho (box filter, SMOOTH_R cells radius) */
    if (SMOOTH_R > 0) {
        double *rho_smooth = calloc(N3, sizeof(double));
        int sr = SMOOTH_R;

        #pragma omp parallel for schedule(static)
        for (long idx = 0; idx < N3; idx++) {
            int i = idx / NN, j = (idx / N) % N, k = idx % N;
            if (i < sr || i >= N-sr || j < sr || j >= N-sr || k < sr || k >= N-sr) {
                rho_smooth[idx] = g->rho[idx];
                continue;
            }
            double sum = 0;
            int count = 0;
            for (int di = -sr; di <= sr; di++)
                for (int dj = -sr; dj <= sr; dj++)
                    for (int dk = -sr; dk <= sr; dk++) {
                        long nidx = (long)(i+di)*NN + (j+dj)*N + (k+dk);
                        sum += g->rho[nidx];
                        count++;
                    }
            rho_smooth[idx] = sum / count;
        }

        memcpy(g->rho, rho_smooth, N3 * sizeof(double));
        free(rho_smooth);
    }

    /* Step 3: c_eff²
       Dense field → SLOWER c (like Schwarzschild: clocks tick slower near mass)
       Depleted field → c approaches c₀ = 1
       c_eff² = 1 / (1 + (rho/rho_ref)^alpha)
       At rho=0: c_eff² = 1 (vacuum speed)
       At rho=rho_ref: c_eff² = 0.5
       At rho>>rho_ref: c_eff² → 0 (frozen inside dense region)
       This creates gravitational lensing: waves slow near mass, curve toward it. */
    #pragma omp parallel for schedule(static)
    for (long idx = 0; idx < N3; idx++) {
        double ratio = g->rho[idx] / (RHO_REF + 1e-30);
        g->c2_eff[idx] = 1.0 / (1.0 + pow(ratio, C_ALPHA));
        if (g->c2_eff[idx] < 0.01) g->c2_eff[idx] = 0.01;  /* floor */
    }
}

/* ================================================================
   Force computation with field-dependent c
   ================================================================ */

static void compute_forces(Grid *g) {
    int N = g->N, NN = N*N;
    long N3 = (long)N*N*N;
    double idx2 = 1.0 / (g->dx * g->dx);
    double mass2 = MASS * MASS;
    int margin = 2;

    #pragma omp parallel for schedule(static)
    for (long idx = 0; idx < N3; idx++) {
        int i = idx / NN, j = (idx / N) % N, k = idx % N;

        /* Pinned boundary: acc=0, phi=0 at edges */
        if (i < margin || i >= N-margin || j < margin || j >= N-margin
            || k < margin || k >= N-margin) {
            for (int a = 0; a < NFIELDS; a++) {
                g->acc[a][idx] = 0;
                g->phi[a][idx] = 0;
                g->vel[a][idx] = 0;
            }
            continue;
        }

        double c2 = g->c2_eff[idx];

        /* Stable quartic self-interaction: V = (lambda/4)(phi²-v²)²
           This replaces the triple-product potential which is unstable at dx>>1.
           The quartic keeps field amplitudes bounded while allowing structure. */
        double phi2_total = 0;
        for (int a = 0; a < NFIELDS; a++)
            phi2_total += g->phi[a][idx] * g->phi[a][idx];
        double lambda_q = 0.1;  /* quartic coupling */
        double v2_vac = 1.0;    /* preferred field magnitude */

        for (int a = 0; a < NFIELDS; a++) {
            double lap = (g->phi[a][idx+NN] + g->phi[a][idx-NN]
                        + g->phi[a][idx+N]  + g->phi[a][idx-N]
                        + g->phi[a][idx+1]  + g->phi[a][idx-1]
                        - 6.0 * g->phi[a][idx]) * idx2;

            /* Quartic force: -dV/dphi_a = -lambda*(phi²-v²)*phi_a */
            double f_quartic = lambda_q * (phi2_total - v2_vac) * g->phi[a][idx];

            /* c_eff² multiplies the Laplacian */
            g->acc[a][idx] = c2 * lap - mass2 * g->phi[a][idx] - f_quartic;
        }
    }
}

/* ================================================================
   Verlet step
   ================================================================ */

static void verlet_step(Grid *g) {
    long N3 = (long)g->N * g->N * g->N;
    double hdt = 0.5 * g->dt, dt = g->dt;

    for (int a = 0; a < NFIELDS; a++) {
        double *p = g->phi[a], *v = g->vel[a], *ac = g->acc[a];
        for (long idx = 0; idx < N3; idx++)
            v[idx] += hdt * ac[idx];
    }
    for (int a = 0; a < NFIELDS; a++) {
        double *p = g->phi[a], *v = g->vel[a];
        for (long idx = 0; idx < N3; idx++)
            p[idx] += dt * v[idx];
    }

    compute_rho_and_ceff(g);
    compute_forces(g);

    for (int a = 0; a < NFIELDS; a++) {
        double *v = g->vel[a], *ac = g->acc[a];
        for (long idx = 0; idx < N3; idx++)
            v[idx] += hdt * ac[idx];
    }
}

/* ================================================================
   Initialization: concentrated blob with small asymmetry
   ================================================================ */

static void init_blob(Grid *g) {
    int N = g->N, NN = N*N;
    double dx = g->dx, L = g->L;
    double k0 = PI / R_INIT;  /* wavenumber matched to blob size */
    double omega0 = sqrt(k0*k0 + MASS*MASS);

    printf("  Init: A0=%.1f, R_init=%.0f, k0=%.4f, omega0=%.4f\n",
           A0_INIT, R_INIT, k0, omega0);

    /* Concentrated blob with helical structure + small random asymmetry */
    unsigned int seed = 42;
    for (int i = 0; i < N; i++) {
        double x = -L + i * dx;
        for (int j = 0; j < N; j++) {
            double y = -L + j * dx;
            for (int k = 0; k < N; k++) {
                double z = -L + k * dx;
                long idx = (long)i * NN + j * N + k;
                double r = sqrt(x*x + y*y + z*z);

                /* Gaussian envelope */
                double env = exp(-r*r / (2.0*R_INIT*R_INIT));

                /* Small asymmetry (1% noise) to seed instabilities */
                seed = seed * 1103515245 + 12345;
                double noise = ((seed >> 16) & 0x7fff) / 32768.0 - 0.5;

                for (int a = 0; a < NFIELDS; a++) {
                    double phase = k0 * z + 2.0*PI*a/3.0;
                    g->phi[a][idx] = A0_INIT * env * cos(phase) * (1.0 + 0.01*noise);
                    g->vel[a][idx] = omega0 * A0_INIT * env * sin(phase) * (1.0 + 0.01*noise);
                    seed = seed * 1103515245 + 12345;
                    noise = ((seed >> 16) & 0x7fff) / 32768.0 - 0.5;
                }
            }
        }
    }
}

/* ================================================================
   Diagnostics
   ================================================================ */

static void write_diagnostics(Grid *g, FILE *fp_ts, double t) {
    int N = g->N, NN = N*N;
    long N3 = (long)N*N*N;
    double dx = g->dx, L = g->L;

    double E_total = 0, max_rho = 0, max_c2 = 0;
    double cx = 0, cy = 0, cz = 0, total_rho = 0;
    double R_rms = 0;
    double phi2_inner = 0, phi2_total = 0;
    double R_inner = 0.1 * L;

    #pragma omp parallel for reduction(+:E_total,cx,cy,cz,total_rho,R_rms,phi2_inner,phi2_total) \
                             reduction(max:max_rho,max_c2) schedule(static)
    for (long idx = 0; idx < N3; idx++) {
        int i = idx / NN, j = (idx / N) % N, k = idx % N;
        double x = -L + i*dx, y = -L + j*dx, z = -L + k*dx;
        double r = sqrt(x*x + y*y + z*z);

        double rho = g->rho[idx];
        double c2  = g->c2_eff[idx];
        double dV  = dx*dx*dx;

        E_total += rho * dV;
        if (rho > max_rho) max_rho = rho;
        if (c2 > max_c2) max_c2 = c2;

        cx += x * rho; cy += y * rho; cz += z * rho;
        total_rho += rho;
        R_rms += r*r * rho;

        double p2 = 0;
        for (int a = 0; a < NFIELDS; a++)
            p2 += g->phi[a][idx] * g->phi[a][idx];
        phi2_total += p2;
        if (r < R_inner) phi2_inner += p2;
    }

    if (total_rho > 0) { cx /= total_rho; cy /= total_rho; cz /= total_rho; }
    double Rrms = (total_rho > 0) ? sqrt(R_rms / total_rho) : 0;
    double fc = (phi2_total > 0) ? phi2_inner / phi2_total : 0;

    /* Count particles (energy peaks above 5× average) */
    double avg_rho = total_rho / N3;
    double thresh = 5.0 * avg_rho;
    int n_peaks = 0;
    for (long idx = 0; idx < N3; idx += 37) /* sparse check */
        if (g->rho[idx] > thresh) n_peaks++;

    fprintf(fp_ts, "%.1f\t%.4e\t%.4e\t%.4e\t%.1f\t%.1f\t%.1f\t%.1f\t%.4f\t%d\t%.4e\n",
            t, E_total, max_rho, max_c2, cx, cy, cz, Rrms, fc, n_peaks, avg_rho);
    fflush(fp_ts);

    printf("t=%7.0f  E=%.2e  max_rho=%.2e  max_c2=%.1f  R_rms=%.0f  fc=%.3f  n_peak=%d\n",
           t, E_total, max_rho, max_c2, Rrms, fc, n_peaks);
    fflush(stdout);
}

/* Write 2D slice through z=0 */
static void write_slice(Grid *g, double t) {
    int N = g->N, NN = N*N;
    double dx = g->dx, L = g->L;
    int k_mid = N/2;
    int sub = 4;  /* subsample 4:1 for manageable file size */

    char fname[512];
    snprintf(fname, sizeof(fname), "%s/slice_t%04d.tsv", OUTDIR, (int)(t+0.5));
    FILE *fp = fopen(fname, "w");
    if (!fp) return;
    fprintf(fp, "x\ty\trho\tc2\n");

    for (int i = 0; i < N; i += sub) {
        double x = -L + i*dx;
        for (int j = 0; j < N; j += sub) {
            double y = -L + j*dx;
            long idx = (long)i*NN + j*N + k_mid;
            fprintf(fp, "%.0f\t%.0f\t%.6e\t%.4f\n",
                    x, y, g->rho[idx], g->c2_eff[idx]);
        }
    }
    fclose(fp);
}

/* ================================================================
   Main
   ================================================================ */

int main(int argc, char **argv) {
    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "-N") && i+1<argc) SIM_N = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-L") && i+1<argc) SIM_L = atof(argv[++i]);
        else if (!strcmp(argv[i], "-A0") && i+1<argc) A0_INIT = atof(argv[++i]);
        else if (!strcmp(argv[i], "-R") && i+1<argc) R_INIT = atof(argv[++i]);
        else if (!strcmp(argv[i], "-alpha") && i+1<argc) C_ALPHA = atof(argv[++i]);
        else if (!strcmp(argv[i], "-rhoref") && i+1<argc) RHO_REF = atof(argv[++i]);
        else if (!strcmp(argv[i], "-T") && i+1<argc) T_TOTAL = atof(argv[++i]);
        else if (!strcmp(argv[i], "-o") && i+1<argc) strncpy(OUTDIR, argv[++i], 255);
        else if (!strcmp(argv[i], "-mu") && i+1<argc) MU = atof(argv[++i]);
        else if (!strcmp(argv[i], "-m") && i+1<argc) MASS = atof(argv[++i]);
    }

    omp_set_num_threads(16);

    printf("=== V30 C-Field: Density-Dependent Speed of Light ===\n");
    printf("N=%d, L=%.0f, dx=%.2f\n", SIM_N, SIM_L, 2*SIM_L/(SIM_N-1));
    printf("A0=%.1f, R_init=%.0f\n", A0_INIT, R_INIT);
    printf("c_eff² = (rho/rho_ref)^alpha, alpha=%.2f, rho_ref=%.2f\n", C_ALPHA, RHO_REF);
    printf("mu=%.2f, kappa=%.1f, mass=%.2f\n", MU, KAPPA, MASS);
    printf("T_total=%.0f, diag_dt=%.0f\n", T_TOTAL, DIAG_DT);
    printf("Output: %s/\n", OUTDIR);

    long N3 = (long)SIM_N * SIM_N * SIM_N;
    printf("Grid: %ld points (%.1f GB)\n\n", N3, N3 * 11.0 * 8 / 1e9);

    mkdir("data", 0755);
    mkdir(OUTDIR, 0755);

    Grid *g = grid_alloc(SIM_N, SIM_L);
    printf("dx=%.2f, dt=%.4f\n", g->dx, g->dt);

    /* Initialize */
    init_blob(g);
    compute_rho_and_ceff(g);

    /* Set rho_ref from initial peak density if not specified */
    if (RHO_REF == 1.0) {
        double max_rho = 0;
        for (long idx = 0; idx < N3; idx++)
            if (g->rho[idx] > max_rho) max_rho = g->rho[idx];
        RHO_REF = max_rho;
        printf("Auto rho_ref = %.4e (initial peak)\n", RHO_REF);
        /* Recompute c_eff with correct rho_ref */
        compute_rho_and_ceff(g);
    }

    compute_forces(g);

    /* Open timeseries */
    char tspath[512];
    snprintf(tspath, sizeof(tspath), "%s/timeseries.tsv", OUTDIR);
    FILE *fp_ts = fopen(tspath, "w");
    fprintf(fp_ts, "t\tE_total\tmax_rho\tmax_c2\tcx\tcy\tcz\tR_rms\tfc\tn_peaks\tavg_rho\n");

    int n_steps = (int)(T_TOTAL / g->dt);
    int diag_steps = (int)(DIAG_DT / g->dt);
    int slice_steps = diag_steps * 2;  /* slices every 2× diag interval */
    if (diag_steps < 1) diag_steps = 1;

    printf("Steps: %d, diag every %d, slice every %d\n\n", n_steps, diag_steps, slice_steps);

    double wall0 = omp_get_wtime();

    /* Initial snapshot */
    write_diagnostics(g, fp_ts, 0);
    write_slice(g, 0);

    for (int step = 1; step <= n_steps; step++) {
        verlet_step(g);

        double t = step * g->dt;

        if (step % diag_steps == 0) {
            write_diagnostics(g, fp_ts, t);
        }
        if (step % slice_steps == 0) {
            write_slice(g, t);
        }
    }

    fclose(fp_ts);

    double wall = omp_get_wtime() - wall0;
    printf("\n=== Complete: %.0f seconds (%.1f min) ===\n", wall, wall/60);

    grid_free(g);
    return 0;
}
