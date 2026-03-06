/*
 * v20/src/triple3d.c
 *
 * 3D simulator for three coupled real scalar fields with saturating
 * triple-product potential. Uses OpenMP.
 *
 * Physics:
 *   L = sum_a [ (1/2)(dt phi_a)^2 - (1/2)|grad phi_a|^2 ] - V
 *   V = (mu/2) P^2 / (1 + kappa P^2),   P = phi_1 * phi_2 * phi_3
 *   EOM: phi_a_tt = laplacian(phi_a) - dV/dphi_a
 *
 * Compile: gcc -O3 -fopenmp -o triple3d v20/src/triple3d.c -lm
 * Run:     ./triple3d -test 1 -N 128 -L 20 -steps 6000 -mu -10 -kappa 0.1
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>

/* --- Grid state --- */
static int    N;
static double L, dx, dt;
static int    ntot;

/* Three scalar fields: position, velocity, acceleration */
static double *phi[3], *vel[3], *acc_buf[3];
static double *damp;  /* absorbing boundary mask */

/* Physics parameters */
static double param_mu    = -10.0;
static double param_kappa = 0.1;

/* Initialization parameters */
static double param_A     = 1.0;
static double param_sigma = 1.5;
static double param_sep   = 3.0;

/* --- Index helpers (NO periodic BC -- Dirichlet zero at boundary) --- */
static inline int idx(int i, int j, int k) {
    return i * N * N + j * N + k;
}
static inline double coord(int i) { return -L / 2 + i * dx; }

/* --- Memory --- */
static void alloc_arrays(void) {
    ntot = N * N * N;
    for (int a = 0; a < 3; a++) {
        phi[a]     = calloc(ntot, sizeof(double));
        vel[a]     = calloc(ntot, sizeof(double));
        acc_buf[a] = calloc(ntot, sizeof(double));
    }
    damp = malloc(ntot * sizeof(double));
}

static void free_arrays(void) {
    for (int a = 0; a < 3; a++) {
        free(phi[a]); free(vel[a]); free(acc_buf[a]);
    }
    free(damp);
}

/* --- Absorbing boundary layer (outer 12%) ---
 * damp=1 in interior, smoothly -> 0.05 at edges. */
static void build_damping(double frac) {
    int b = (int)(N * frac);
    if (b < 2) b = 2;
    #pragma omp parallel for collapse(3)
    for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
    for (int k = 0; k < N; k++) {
        double m = 1.0;
        int ii[3] = {i, j, k};
        for (int d = 0; d < 3; d++) {
            if (ii[d] < b)
                m *= (double)ii[d] / b;
            else if (ii[d] >= N - b)
                m *= (double)(N - 1 - ii[d]) / b;
        }
        damp[i * N * N + j * N + k] = 1.0 - 0.95 * (1.0 - m);
    }
}

/* ====================== INITIALIZATIONS ====================== */

/* Test 1: Three overlapping Gaussians (all centered at origin) */
static void init_overlapping(double A, double sigma) {
    #pragma omp parallel for collapse(3)
    for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
    for (int k = 0; k < N; k++) {
        double x = coord(i), y = coord(j), z = coord(k);
        double r2 = x * x + y * y + z * z;
        int p = idx(i, j, k);
        double val = A * exp(-r2 / (2.0 * sigma * sigma));
        phi[0][p] = val;
        phi[1][p] = val;
        phi[2][p] = val;
        vel[0][p] = vel[1][p] = vel[2][p] = 0.0;
    }
}

/* Test 2: Three separated Gaussians along x, y, z axes */
static void init_separated(double A, double sigma, double sep) {
    #pragma omp parallel for collapse(3)
    for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
    for (int k = 0; k < N; k++) {
        double x = coord(i), y = coord(j), z = coord(k);
        int p = idx(i, j, k);
        double s2 = 2.0 * sigma * sigma;
        double dx1 = x - sep, dy1 = y, dz1 = z;
        phi[0][p] = A * exp(-(dx1 * dx1 + dy1 * dy1 + dz1 * dz1) / s2);
        double dx2 = x, dy2 = y - sep, dz2 = z;
        phi[1][p] = A * exp(-(dx2 * dx2 + dy2 * dy2 + dz2 * dz2) / s2);
        double dx3 = x, dy3 = y, dz3 = z - sep;
        phi[2][p] = A * exp(-(dx3 * dx3 + dy3 * dy3 + dz3 * dz3) / s2);
        vel[0][p] = vel[1][p] = vel[2][p] = 0.0;
    }
}

/* Test 3: Single-field control (phi_2 = phi_3 = 0, V = 0) */
static void init_single(double A, double sigma) {
    #pragma omp parallel for collapse(3)
    for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
    for (int k = 0; k < N; k++) {
        double x = coord(i), y = coord(j), z = coord(k);
        double r2 = x * x + y * y + z * z;
        int p = idx(i, j, k);
        phi[0][p] = A * exp(-r2 / (2.0 * sigma * sigma));
        phi[1][p] = 0.0;
        phi[2][p] = 0.0;
        vel[0][p] = vel[1][p] = vel[2][p] = 0.0;
    }
}

/* Test 4: Two-field control (phi_3 = 0, V = 0) */
static void init_twofield(double A, double sigma) {
    #pragma omp parallel for collapse(3)
    for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
    for (int k = 0; k < N; k++) {
        double x = coord(i), y = coord(j), z = coord(k);
        double r2 = x * x + y * y + z * z;
        int p = idx(i, j, k);
        double val = A * exp(-r2 / (2.0 * sigma * sigma));
        phi[0][p] = val;
        phi[1][p] = val;
        phi[2][p] = 0.0;
        vel[0][p] = vel[1][p] = vel[2][p] = 0.0;
    }
}

/* ====================== CORE PHYSICS ====================== */

/* Compute acceleration: acc[a] = laplacian(phi[a]) - dV/dphi_a */
static void compute_accel(void) {
    double idxsq = 1.0 / (dx * dx);

    #pragma omp parallel for collapse(3)
    for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
    for (int k = 0; k < N; k++) {
        int p = idx(i, j, k);

        /* 7-point Laplacian with zero Dirichlet BC at boundaries */
        for (int a = 0; a < 3; a++) {
            double center = phi[a][p];
            double lap = -6.0 * center;
            lap += (i > 0)     ? phi[a][idx(i - 1, j, k)] : 0.0;
            lap += (i < N - 1) ? phi[a][idx(i + 1, j, k)] : 0.0;
            lap += (j > 0)     ? phi[a][idx(i, j - 1, k)] : 0.0;
            lap += (j < N - 1) ? phi[a][idx(i, j + 1, k)] : 0.0;
            lap += (k > 0)     ? phi[a][idx(i, j, k - 1)] : 0.0;
            lap += (k < N - 1) ? phi[a][idx(i, j, k + 1)] : 0.0;
            lap *= idxsq;

            /* Potential force: dV/dphi_a = mu * P * (product of other two) / (1+kappa*P^2)^2 */
            double p1 = phi[0][p], p2 = phi[1][p], p3 = phi[2][p];
            double P = p1 * p2 * p3;
            double kP2 = param_kappa * P * P;
            double denom = (1.0 + kP2);
            denom *= denom;  /* (1 + kappa*P^2)^2 */
            double other;
            if (a == 0)      other = p2 * p3;
            else if (a == 1) other = p1 * p3;
            else             other = p1 * p2;
            double dVda = param_mu * P * other / denom;

            acc_buf[a][p] = lap - dVda;
        }
    }
}

/* Velocity-Verlet: half-kick, drift, recompute accel, half-kick, damp */
static void step_vv(void) {
    /* Half-kick */
    compute_accel();
    #pragma omp parallel for
    for (int p = 0; p < ntot; p++)
        for (int a = 0; a < 3; a++)
            vel[a][p] += 0.5 * dt * acc_buf[a][p];

    /* Drift */
    #pragma omp parallel for
    for (int p = 0; p < ntot; p++)
        for (int a = 0; a < 3; a++)
            phi[a][p] += dt * vel[a][p];

    /* Recompute acceleration at new position */
    compute_accel();
    #pragma omp parallel for
    for (int p = 0; p < ntot; p++)
        for (int a = 0; a < 3; a++)
            vel[a][p] += 0.5 * dt * acc_buf[a][p];

    /* Absorbing boundary damping */
    #pragma omp parallel for
    for (int p = 0; p < ntot; p++) {
        double d = damp[p];
        for (int a = 0; a < 3; a++)
            vel[a][p] *= d;
    }
}

/* ====================== DIAGNOSTICS ====================== */

/*
 * Compute all energies and per-field peaks in one pass.
 * E_kin, E_grad, E_pot: total kinetic, gradient, potential
 * E_field[3]: per-field kinetic+gradient (no potential)
 * phi_peak[3]: max|phi_a|
 * f_core: energy fraction within r < R_core
 */
static void compute_diagnostics(double *Ekin, double *Egrad, double *Epot,
                                double E_field[3], double phi_peak[3],
                                double *f_core, double R_core) {
    double ek = 0, eg = 0, ep = 0;
    double ef[3] = {0, 0, 0};
    double pk[3] = {0, 0, 0};
    double e_core = 0, e_total = 0;
    double idx2 = 1.0 / (2.0 * dx);
    double h3 = dx * dx * dx;
    double R2 = R_core * R_core;

    #pragma omp parallel for collapse(3) \
        reduction(+:ek, eg, ep, ef[:3], e_core, e_total) \
        reduction(max:pk[:3])
    for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
    for (int k = 0; k < N; k++) {
        int p = idx(i, j, k);
        double x = coord(i), y = coord(j), z = coord(k);

        /* Per-field kinetic and gradient */
        double ek_loc = 0, eg_loc = 0;
        for (int a = 0; a < 3; a++) {
            double v2 = vel[a][p] * vel[a][p];
            ek_loc += 0.5 * v2;
            ef[a] += 0.5 * v2;

            /* Gradient via central differences (zero BC at edges) */
            double g2 = 0;
            double c = phi[a][p];
            double gx = ((i > 0 && i < N - 1)
                ? (phi[a][idx(i + 1, j, k)] - phi[a][idx(i - 1, j, k)]) * idx2
                : (i == 0 ? (phi[a][idx(1, j, k)] - 0.0) / dx
                          : (0.0 - phi[a][idx(N - 2, j, k)]) / dx));
            double gy = ((j > 0 && j < N - 1)
                ? (phi[a][idx(i, j + 1, k)] - phi[a][idx(i, j - 1, k)]) * idx2
                : (j == 0 ? (phi[a][idx(i, 1, k)] - 0.0) / dx
                          : (0.0 - phi[a][idx(i, N - 2, k)]) / dx));
            double gz = ((k > 0 && k < N - 1)
                ? (phi[a][idx(i, j, k + 1)] - phi[a][idx(i, j, k - 1)]) * idx2
                : (k == 0 ? (phi[a][idx(i, j, 1)] - 0.0) / dx
                          : (0.0 - phi[a][idx(i, j, N - 2)]) / dx));
            (void)c;
            g2 = gx * gx + gy * gy + gz * gz;
            eg_loc += 0.5 * g2;
            ef[a] += 0.5 * g2;

            /* Peak tracking */
            double absv = fabs(phi[a][p]);
            if (absv > pk[a]) pk[a] = absv;
        }
        ek += ek_loc;
        eg += eg_loc;

        /* Potential */
        double P = phi[0][p] * phi[1][p] * phi[2][p];
        double kP2 = param_kappa * P * P;
        double vp = 0.5 * param_mu * P * P / (1.0 + kP2);
        ep += vp;

        /* Core fraction */
        double ed = ek_loc + eg_loc + vp;
        e_total += ed;
        if (x * x + y * y + z * z < R2)
            e_core += ed;
    }

    *Ekin  = ek * h3;
    *Egrad = eg * h3;
    *Epot  = ep * h3;
    for (int a = 0; a < 3; a++) {
        E_field[a]  = ef[a] * h3;
        phi_peak[a] = pk[a];
    }
    *f_core = (fabs(e_total) > 1e-30) ? e_core / e_total : 0.0;
}

/* --- Linear regression: y = slope*x + intercept --- */
static void linregress(const double *x, const double *y, int n,
                       double *slope, double *intercept, double *r2) {
    double sx = 0, sy = 0, sxx = 0, sxy = 0;
    for (int i = 0; i < n; i++) {
        sx += x[i]; sy += y[i]; sxx += x[i] * x[i]; sxy += x[i] * y[i];
    }
    double den = n * sxx - sx * sx;
    if (fabs(den) < 1e-30) { *slope = 0; *intercept = 0; *r2 = 0; return; }
    *slope = (n * sxy - sx * sy) / den;
    *intercept = (sy * sxx - sx * sxy) / den;
    double ymean = sy / n, sstot = 0, ssres = 0;
    for (int i = 0; i < n; i++) {
        double yp = *slope * x[i] + *intercept;
        sstot += (y[i] - ymean) * (y[i] - ymean);
        ssres += (y[i] - yp) * (y[i] - yp);
    }
    *r2 = (sstot > 0) ? 1.0 - ssres / sstot : 0.0;
}

/* ====================== MAIN ====================== */
int main(int argc, char **argv) {
    /* Defaults */
    int    n = 128, steps = 6000, test = 1;
    double l = 20.0;
    char   outdir[256] = "v20/data";

    for (int i = 1; i < argc; i++) {
        if      (!strcmp(argv[i], "-N")      && i + 1 < argc) n = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-L")      && i + 1 < argc) l = atof(argv[++i]);
        else if (!strcmp(argv[i], "-steps")  && i + 1 < argc) steps = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-test")   && i + 1 < argc) test = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-mu")     && i + 1 < argc) param_mu = atof(argv[++i]);
        else if (!strcmp(argv[i], "-kappa")  && i + 1 < argc) param_kappa = atof(argv[++i]);
        else if (!strcmp(argv[i], "-A")      && i + 1 < argc) param_A = atof(argv[++i]);
        else if (!strcmp(argv[i], "-sigma")  && i + 1 < argc) param_sigma = atof(argv[++i]);
        else if (!strcmp(argv[i], "-sep")    && i + 1 < argc) param_sep = atof(argv[++i]);
        else if (!strcmp(argv[i], "-o")      && i + 1 < argc) strncpy(outdir, argv[++i], 255);
    }

    /* Set up grid */
    N  = n;
    L  = l;
    dx = L / (N - 1);
    ntot = N * N * N;

    /* CFL: dt = 0.25 * dx/sqrt(3) */
    dt = 0.25 * dx / sqrt(3.0);

    double mem_gb = 3.0 * 3.0 * (double)ntot * 8.0 / 1e9 + (double)ntot * 8.0 / 1e9;
    printf("=== Triple Product 3D (v20) ===\n");
    printf("N=%d  L=%.1f  dx=%.5f  dt=%.6f\n", N, L, dx, dt);
    printf("mu=%.4f  kappa=%.4f  A=%.2f  sigma=%.2f  sep=%.2f\n",
           param_mu, param_kappa, param_A, param_sigma, param_sep);
    printf("Test %d  steps=%d  threads=%d  mem=%.2f GB\n",
           test, steps, omp_get_max_threads(), mem_gb);

    setlinebuf(stdout);

    alloc_arrays();
    build_damping(0.12);

    double R_core = 2.0 * param_sigma;

    /* Initialize */
    const char *test_names[] = {
        "", "Overlapping Gaussians", "Separated Gaussians",
        "Single-Field Control", "Two-Field Control"
    };
    if (test < 1 || test > 4) {
        fprintf(stderr, "Unknown test %d (valid: 1-4)\n", test);
        return 1;
    }
    printf("\n--- Test %d: %s ---\n", test, test_names[test]);

    switch (test) {
        case 1: init_overlapping(param_A, param_sigma); break;
        case 2: init_separated(param_A, param_sigma, param_sep); break;
        case 3: init_single(param_A, param_sigma); break;
        case 4: init_twofield(param_A, param_sigma); break;
    }

    /* Initial diagnostics */
    double Ek, Eg, Ep, Ef[3], pk[3], fc;
    compute_diagnostics(&Ek, &Eg, &Ep, Ef, pk, &fc, R_core);
    double Et = Ek + Eg + Ep;

    printf("Initial: Ek=%.4f Eg=%.4f Ep=%.4f Et=%.4f\n", Ek, Eg, Ep, Et);
    printf("Per-field: E1=%.4f E2=%.4f E3=%.4f\n", Ef[0], Ef[1], Ef[2]);
    printf("Peaks: %.4f %.4f %.4f  f_core=%.4f\n", pk[0], pk[1], pk[2], fc);

    /* --- Output files --- */
    char fname_ts[512];
    snprintf(fname_ts, sizeof(fname_ts), "%s/test%d_timeseries.tsv", outdir, test);
    FILE *fp_ts = fopen(fname_ts, "w");
    if (!fp_ts) {
        char cmd[512];
        snprintf(cmd, sizeof(cmd), "mkdir -p %s", outdir);
        if (system(cmd) != 0) { fprintf(stderr, "mkdir failed\n"); return 1; }
        fp_ts = fopen(fname_ts, "w");
        if (!fp_ts) { perror("fopen timeseries"); return 1; }
    }

    fprintf(fp_ts, "step\ttime\tE_kin\tE_grad\tE_pot\tE_total\tE1\tE2\tE3\tf_core\tphi1_peak\tphi2_peak\tphi3_peak\n");
    fprintf(fp_ts, "%d\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n",
            0, 0.0, Ek, Eg, Ep, Et, Ef[0], Ef[1], Ef[2], fc, pk[0], pk[1], pk[2]);

    /* --- Time-averaged shell accumulator --- */
    #define NBINS 100
    double tavg_phi2[NBINS], tavg_rhoE[NBINS];
    int    tavg_cnt[NBINS];
    for (int b = 0; b < NBINS; b++) {
        tavg_phi2[b] = tavg_rhoE[b] = 0.0;
        tavg_cnt[b] = 0;
    }
    int n_tavg = 0;
    int settle_step = (int)(0.2 * steps);
    double rmax_bin = L / 2.0;
    double dr_bin = rmax_bin / NBINS;

    int print_every = steps / 60;
    if (print_every < 1) print_every = 1;

    /* --- Time evolution --- */
    double t_start = omp_get_wtime();

    for (int s = 1; s <= steps; s++) {
        step_vv();

        if (s % print_every == 0 || s == steps) {
            double t = s * dt;

            compute_diagnostics(&Ek, &Eg, &Ep, Ef, pk, &fc, R_core);
            Et = Ek + Eg + Ep;

            fprintf(fp_ts, "%d\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n",
                    s, t, Ek, Eg, Ep, Et, Ef[0], Ef[1], Ef[2], fc, pk[0], pk[1], pk[2]);
            fflush(fp_ts);

            /* Time-averaged shell accumulation (after settling) */
            if (s > settle_step) {
                double idx2_loc = 1.0 / (2.0 * dx);

                double step_phi2[NBINS], step_rhoE[NBINS];
                int    step_cnt[NBINS];
                for (int b = 0; b < NBINS; b++) {
                    step_phi2[b] = step_rhoE[b] = 0.0;
                    step_cnt[b] = 0;
                }

                for (int ii = 0; ii < N; ii++)
                for (int jj = 0; jj < N; jj++)
                for (int kk = 0; kk < N; kk++) {
                    double x = coord(ii), y = coord(jj), z = coord(kk);
                    double r = sqrt(x * x + y * y + z * z);
                    int b = (int)(r / dr_bin);
                    if (b < 0 || b >= NBINS) continue;
                    int pp = idx(ii, jj, kk);

                    /* Field amplitude squared */
                    double f2 = phi[0][pp] * phi[0][pp]
                              + phi[1][pp] * phi[1][pp]
                              + phi[2][pp] * phi[2][pp];
                    step_phi2[b] += f2;

                    /* Energy density */
                    double ek_loc = 0, eg_loc = 0;
                    for (int a = 0; a < 3; a++) {
                        ek_loc += 0.5 * vel[a][pp] * vel[a][pp];
                        double gx = ((ii > 0 && ii < N - 1)
                            ? (phi[a][idx(ii + 1, jj, kk)] - phi[a][idx(ii - 1, jj, kk)]) * idx2_loc
                            : 0.0);
                        double gy = ((jj > 0 && jj < N - 1)
                            ? (phi[a][idx(ii, jj + 1, kk)] - phi[a][idx(ii, jj - 1, kk)]) * idx2_loc
                            : 0.0);
                        double gz = ((kk > 0 && kk < N - 1)
                            ? (phi[a][idx(ii, jj, kk + 1)] - phi[a][idx(ii, jj, kk - 1)]) * idx2_loc
                            : 0.0);
                        eg_loc += 0.5 * (gx * gx + gy * gy + gz * gz);
                    }
                    double P = phi[0][pp] * phi[1][pp] * phi[2][pp];
                    double kP2 = param_kappa * P * P;
                    double vp = 0.5 * param_mu * P * P / (1.0 + kP2);
                    double ed = ek_loc + eg_loc + vp;

                    step_rhoE[b] += ed;
                    step_cnt[b]++;
                }

                for (int b = 0; b < NBINS; b++) {
                    if (step_cnt[b] > 0) {
                        tavg_phi2[b] += step_phi2[b] / step_cnt[b];
                        tavg_rhoE[b] += step_rhoE[b] / step_cnt[b];
                        tavg_cnt[b]++;
                    }
                }
                n_tavg++;
            }

            /* Console output */
            printf("Step %5d | t=%7.3f | Et=%9.3f | Epot=%9.3f | phi_peaks=%.4f %.4f %.4f | fc=%.3f\n",
                   s, t, Et, Ep, pk[0], pk[1], pk[2], fc);
        }
    }

    fclose(fp_ts);

    double elapsed = omp_get_wtime() - t_start;
    printf("\nDone: %d steps in %.1f s (%.1f steps/s)\n", steps, elapsed, steps / elapsed);
    printf("Timeseries: %s\n", fname_ts);

    /* ====================== POST-PROCESSING ====================== */

    /* --- Write time-averaged radial profile --- */
    char fname_tavg[512];
    snprintf(fname_tavg, sizeof(fname_tavg), "%s/test%d_tavg_profile.tsv", outdir, test);
    FILE *fp_tavg = fopen(fname_tavg, "w");
    if (fp_tavg) {
        fprintf(fp_tavg, "r\tavg_phi2\tavg_rhoE\n");
        for (int b = 0; b < NBINS; b++) {
            if (tavg_cnt[b] > 0) {
                double r_mid = (b + 0.5) * dr_bin;
                fprintf(fp_tavg, "%.4f\t%.8e\t%.8e\n",
                        r_mid,
                        tavg_phi2[b] / tavg_cnt[b],
                        tavg_rhoE[b] / tavg_cnt[b]);
            }
        }
        fclose(fp_tavg);
        printf("Time-averaged profile: %s\n", fname_tavg);
    }

    /* --- Far-field power-law fits --- */
    char fname_far[512];
    snprintf(fname_far, sizeof(fname_far), "%s/test%d_farfield.txt", outdir, test);
    FILE *fp_far = fopen(fname_far, "w");
    if (fp_far) {
        double fit_rmin = 2.0 * R_core;
        double fit_rmax = L / 2 - 0.12 * L;

        double logr[NBINS], logp[NBINS], logE[NBINS];
        int npts = 0;
        for (int b = 0; b < NBINS; b++) {
            if (tavg_cnt[b] < 1) continue;
            double r_mid = (b + 0.5) * dr_bin;
            if (r_mid < fit_rmin || r_mid > fit_rmax) continue;
            double avg_p = tavg_phi2[b] / tavg_cnt[b];
            double avg_e = tavg_rhoE[b] / tavg_cnt[b];
            if (avg_p < 1e-20 || fabs(avg_e) < 1e-20) continue;
            logr[npts] = log(r_mid);
            logp[npts] = log(avg_p);
            logE[npts] = log(fabs(avg_e));
            npts++;
        }

        fprintf(fp_far, "Far-field power-law fits (r in [%.1f, %.1f], %d bins, %d time samples)\n\n",
                fit_rmin, fit_rmax, npts, n_tavg);

        if (npts >= 5) {
            double sl, inter, r2;

            linregress(logr, logp, npts, &sl, &inter, &r2);
            fprintf(fp_far, "avg_phi2 ~ r^{%.3f}  R^2=%.4f  (free radiation: -2)\n", sl, r2);
            printf("\nFar-field: avg_phi2 ~ r^{%.3f} R^2=%.4f\n", sl, r2);

            linregress(logr, logE, npts, &sl, &inter, &r2);
            fprintf(fp_far, "avg_rhoE ~ r^{%.3f}  R^2=%.4f  (free radiation: -2)\n", sl, r2);
            printf("Far-field: avg_rhoE ~ r^{%.3f} R^2=%.4f\n", sl, r2);
        } else {
            fprintf(fp_far, "Not enough data for fit (%d points, need >= 5)\n", npts);
            printf("Not enough far-field data for power-law fit (%d points)\n", npts);
        }
        fclose(fp_far);
        printf("Far-field fits: %s\n", fname_far);
    }

    free_arrays();
    return 0;
}
