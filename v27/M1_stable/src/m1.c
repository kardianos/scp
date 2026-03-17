/*
 * m1.c — V27-M1: Stable Propagating Braid
 *
 * Based on V26-DynA (propagating helical wave braid).
 * Additions:
 *   - Pairwise coupling: lambda_pw * (phi1*phi2 + phi2*phi3 + phi3*phi1)
 *   - Configurable k_wave (number of twists)
 *   - Test modes: M1a (lambda_pw scan), M1b (L=10), M1c (k scan)
 *
 * Lagrangian:
 *   L = Sum_a [1/2(dt phi_a)^2 - 1/2|grad phi_a|^2 - 1/2 m^2 phi_a^2]
 *     - (mu/2) P^2 / (1 + kappa P^2)
 *     - lambda_pw * (phi_0*phi_1 + phi_1*phi_2 + phi_2*phi_0)
 *
 * Compile: gcc -O3 -fopenmp -Wall -o m1 src/m1.c -lm
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

/* --- Parameters --- */
static double mu_pot    = -20.0;
static double kappa     = 20.0;
static double mass      = 1.0;
static double A0        = 0.8;
static double R_tube    = 3.0;
static double lambda_pw = 0.0;   /* pairwise coupling */

static int    N         = 128;
static double L         = 20.0;
static double tfinal    = 500.0;
static double cfl_frac  = 0.20;
static char   outdir[512] = "data";

/* --- Index helpers --- */
#define IDX(i,j,k) ((long)(i)*N*N + (long)(j)*N + (long)(k))

/* --- Globals --- */
static double *phi[3], *vel[3], *acc[3];
static double *damp;
static double dx, dx2, m2, dt;
static long Ngrid;

/* --- Periodic z wrap --- */
static inline int wrap_z(int k)
{
    if (k < 0)   return k + N;
    if (k >= N)  return k - N;
    return k;
}

/* --- Initialization: Propagating Helical Braid --- */
static void init_braid(double k_w, double omega)
{
    printf("  k = %.6f, omega = %.6f, v_g = k/omega = %.4f\n",
           k_w, omega, k_w / omega);

    #pragma omp parallel for schedule(static)
    for (long idx = 0; idx < Ngrid; idx++) {
        int ii = idx / (N * N);
        int jj = (idx / N) % N;
        int kk = idx % N;
        double x = -L + ii * dx;
        double y = -L + jj * dx;
        double z = -L + kk * dx;

        double r_perp = sqrt(x * x + y * y);
        double envelope = A0 * exp(-r_perp * r_perp / (2.0 * R_tube * R_tube));

        for (int a = 0; a < 3; a++) {
            double phase = k_w * z + 2.0 * M_PI * a / 3.0;
            phi[a][idx] = envelope * cos(phase);
            vel[a][idx] = omega * envelope * sin(phase);
        }
    }
}

/* --- Compute acceleration (triple product + pairwise, periodic z) --- */
static void compute_acc(void)
{
    for (int a = 0; a < 3; a++) {
        int b = (a + 1) % 3;
        int c = (a + 2) % 3;

        #pragma omp parallel for schedule(static)
        for (int i = 2; i < N-2; i++) {
            for (int j = 2; j < N-2; j++) {
                for (int k = 0; k < N; k++) {
                    long idx = IDX(i, j, k);

                    int km1 = wrap_z(k - 1);
                    int kp1 = wrap_z(k + 1);

                    double lapl = (phi[a][IDX(i+1,j,k)] + phi[a][IDX(i-1,j,k)]
                                 + phi[a][IDX(i,j+1,k)] + phi[a][IDX(i,j-1,k)]
                                 + phi[a][IDX(i,j,kp1)] + phi[a][IDX(i,j,km1)]
                                 - 6.0 * phi[a][idx]) / dx2;

                    /* Triple-product potential force */
                    double p0 = phi[0][idx], p1 = phi[1][idx], p2 = phi[2][idx];
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
                    double dVdphi_triple = mu_pot * P * dP / denom2;

                    /* Pairwise coupling force:
                     * V_pw = lambda_pw * (phi0*phi1 + phi1*phi2 + phi2*phi0)
                     * dV_pw/dphi_a = lambda_pw * (phi_b + phi_c) */
                    double dVdphi_pw = lambda_pw * (phi[b][idx] + phi[c][idx]);

                    acc[a][idx] = lapl - m2 * phi[a][idx] - dVdphi_triple - dVdphi_pw;
                }
            }
        }

        /* Boundary: zero acceleration for x,y edges */
        #pragma omp parallel for schedule(static)
        for (long idx = 0; idx < Ngrid; idx++) {
            int ii = idx / (N * N);
            int jj = (idx / N) % N;
            if (ii < 2 || ii >= N-2 || jj < 2 || jj >= N-2) {
                acc[a][idx] = 0.0;
            }
        }
    }
}

/* --- Velocity Verlet step --- */
static void verlet_step(void)
{
    for (int a = 0; a < 3; a++) {
        #pragma omp parallel for schedule(static)
        for (long idx = 0; idx < Ngrid; idx++)
            vel[a][idx] += 0.5 * dt * acc[a][idx];
    }

    for (int a = 0; a < 3; a++) {
        #pragma omp parallel for schedule(static)
        for (long idx = 0; idx < Ngrid; idx++)
            phi[a][idx] += dt * vel[a][idx];
    }

    compute_acc();

    for (int a = 0; a < 3; a++) {
        #pragma omp parallel for schedule(static)
        for (long idx = 0; idx < Ngrid; idx++)
            vel[a][idx] += 0.5 * dt * acc[a][idx];
    }

    /* Absorbing boundary damping (x,y only) */
    for (int a = 0; a < 3; a++) {
        #pragma omp parallel for schedule(static)
        for (long idx = 0; idx < Ngrid; idx++) {
            vel[a][idx] *= damp[idx];
            phi[a][idx] *= damp[idx];
        }
    }
}

/* --- Diagnostics structure --- */
typedef struct {
    double Ek, Eg, Em, Ep, Epw, Et;
    double peak[3];
    double fc;
    double peak_P;
    double rho_center;
    double Pz;
    double l2_frac;
} Diag;

/* --- Legendre polynomials --- */
static double legendre_P_fn(int l, double x)
{
    switch (l) {
        case 0: return 1.0;
        case 1: return x;
        case 2: return 0.5*(3.0*x*x - 1.0);
        default: return 0.0;
    }
}

/* --- Compute diagnostics --- */
static Diag compute_diag(double core_radius)
{
    Diag d;
    memset(&d, 0, sizeof(d));

    double Ek = 0, Eg = 0, Em = 0, Ep = 0, Epw = 0;
    double Ecore = 0, Eall = 0;
    double peak[3] = {0, 0, 0};
    double peak_P = 0;
    double Pz_total = 0;

    /* l=2 multipole */
    double c_l[3] = {0, 0, 0};

    #pragma omp parallel
    {
        double lEk = 0, lEg = 0, lEm = 0, lEp = 0, lEpw = 0, lEc = 0, lEa = 0;
        double lpk[3] = {0, 0, 0};
        double lpkP = 0;
        double lPz = 0;
        double lc_l[3] = {0, 0, 0};

        #pragma omp for schedule(static) nowait
        for (int i = 2; i < N-2; i++) {
            double x = -L + i * dx;
            for (int j = 2; j < N-2; j++) {
                double y = -L + j * dx;
                for (int k = 0; k < N; k++) {
                    double z = -L + k * dx;
                    long idx = IDX(i, j, k);
                    double dV = dx * dx * dx;
                    double e_loc = 0;

                    int km1 = wrap_z(k - 1);
                    int kp1 = wrap_z(k + 1);

                    for (int a = 0; a < 3; a++) {
                        double v2 = vel[a][idx] * vel[a][idx];
                        lEk += 0.5 * v2 * dV;
                        e_loc += 0.5 * v2;

                        double gx = (phi[a][IDX(i+1,j,k)] - phi[a][IDX(i-1,j,k)]) / (2*dx);
                        double gy = (phi[a][IDX(i,j+1,k)] - phi[a][IDX(i,j-1,k)]) / (2*dx);
                        double gz = (phi[a][IDX(i,j,kp1)] - phi[a][IDX(i,j,km1)]) / (2*dx);
                        double grad2 = gx*gx + gy*gy + gz*gz;
                        lEg += 0.5 * grad2 * dV;
                        e_loc += 0.5 * grad2;

                        double mass_e = 0.5 * m2 * phi[a][idx] * phi[a][idx];
                        lEm += mass_e * dV;
                        e_loc += mass_e;

                        double ap = fabs(phi[a][idx]);
                        if (ap > lpk[a]) lpk[a] = ap;

                        /* z-momentum: -vel_a * dz(phi_a) */
                        lPz += -vel[a][idx] * gz * dV;
                    }

                    double p0 = phi[0][idx], p1 = phi[1][idx], p2 = phi[2][idx];
                    double P_val = p0 * p1 * p2;
                    double P2 = P_val * P_val;
                    double Vloc = 0.5 * mu_pot * P2 / (1.0 + kappa * P2);
                    lEp += Vloc * dV;
                    e_loc += Vloc;

                    /* Pairwise energy */
                    double Vpw = lambda_pw * (p0*p1 + p1*p2 + p2*p0);
                    lEpw += Vpw * dV;
                    e_loc += Vpw;

                    double absP = fabs(P_val);
                    if (absP > lpkP) lpkP = absP;

                    lEa += e_loc * dV;
                    double r = sqrt(x*x + y*y + z*z);
                    if (r < core_radius) lEc += e_loc * dV;

                    /* Multipole: energy density on shells */
                    if (r > 0.5) {
                        double cos_th = z / r;
                        for (int l = 0; l <= 2; l++)
                            lc_l[l] += e_loc * legendre_P_fn(l, cos_th) * dV;
                    }
                }
            }
        }

        #pragma omp critical
        {
            Ek += lEk; Eg += lEg; Em += lEm; Ep += lEp; Epw += lEpw;
            Ecore += lEc; Eall += lEa;
            Pz_total += lPz;
            for (int a = 0; a < 3; a++)
                if (lpk[a] > peak[a]) peak[a] = lpk[a];
            if (lpkP > peak_P) peak_P = lpkP;
            for (int l = 0; l <= 2; l++) c_l[l] += lc_l[l];
        }
    }

    d.Ek = Ek; d.Eg = Eg; d.Em = Em; d.Ep = Ep; d.Epw = Epw;
    d.Et = Ek + Eg + Em + Ep + Epw;
    for (int a = 0; a < 3; a++) d.peak[a] = peak[a];
    d.fc = (Eall > 1e-20) ? Ecore / Eall : 0.0;
    d.peak_P = peak_P;
    d.Pz = Pz_total;

    double sum_abs = fabs(c_l[0]) + fabs(c_l[1]) + fabs(c_l[2]);
    d.l2_frac = (sum_abs > 1e-30) ? fabs(c_l[2]) / sum_abs : 0.0;

    /* Center energy density */
    {
        long ic = IDX(N/2, N/2, N/2);
        int ci = N/2, cj = N/2, ck = N/2;
        int ckm1 = wrap_z(ck - 1);
        int ckp1 = wrap_z(ck + 1);
        double rho = 0;
        for (int a = 0; a < 3; a++) {
            rho += 0.5 * vel[a][ic] * vel[a][ic];
            double gx = (phi[a][IDX(ci+1,cj,ck)] - phi[a][IDX(ci-1,cj,ck)]) / (2*dx);
            double gy = (phi[a][IDX(ci,cj+1,ck)] - phi[a][IDX(ci,cj-1,ck)]) / (2*dx);
            double gz = (phi[a][IDX(ci,cj,ckp1)] - phi[a][IDX(ci,cj,ckm1)]) / (2*dx);
            rho += 0.5 * (gx*gx + gy*gy + gz*gz);
            rho += 0.5 * m2 * phi[a][ic] * phi[a][ic];
        }
        double p0 = phi[0][ic], p1 = phi[1][ic], p2 = phi[2][ic];
        double P_val = p0 * p1 * p2;
        double P2 = P_val * P_val;
        rho += 0.5 * mu_pot * P2 / (1.0 + kappa * P2);
        rho += lambda_pw * (p0*p1 + p1*p2 + p2*p0);
        d.rho_center = rho;
    }

    return d;
}

/* --- Run a single test configuration --- */
static void run_test(const char *label, double k_w, double Lval, double lpw)
{
    /* Update globals for this run */
    L = Lval;
    lambda_pw = lpw;
    dx  = 2.0 * L / (N - 1);
    dx2 = dx * dx;
    m2  = mass * mass;
    dt  = cfl_frac * dx;

    double omega = sqrt(k_w * k_w + mass * mass);
    double v_g = k_w / omega;

    printf("\n");
    printf("================================================================\n");
    printf("  M1 Test: %s\n", label);
    printf("  L=%.1f  lambda_pw=%.3f  k=%.4f  omega=%.4f  v_g=%.4f  dx=%.4f  dt=%.5f\n",
           L, lambda_pw, k_w, omega, v_g, dx, dt);
    printf("================================================================\n\n");

    /* Zero fields */
    for (int a = 0; a < 3; a++) {
        memset(phi[a], 0, Ngrid * sizeof(double));
        memset(vel[a], 0, Ngrid * sizeof(double));
        memset(acc[a], 0, Ngrid * sizeof(double));
    }

    /* Initialize */
    init_braid(k_w, omega);

    /* Set damping: absorbing in x,y only */
    double R_abs_inner = L * 0.70;
    double R_abs_outer = L * 0.95;
    #pragma omp parallel for schedule(static)
    for (long idx = 0; idx < Ngrid; idx++) {
        int ii = idx / (N * N);
        int jj = (idx / N) % N;
        double x = -L + ii * dx;
        double y = -L + jj * dx;
        double rperp = sqrt(x*x + y*y);
        if (rperp > R_abs_inner) {
            double f = (rperp - R_abs_inner) / (R_abs_outer - R_abs_inner);
            if (f > 1.0) f = 1.0;
            damp[idx] = 1.0 - 0.98 * f * f;
        } else {
            damp[idx] = 1.0;
        }
    }

    double core_radius = (L < 15.0) ? L * 0.4 : 8.0;

    int Nt = (int)(tfinal / dt) + 1;
    int print_every = Nt / 20;
    if (print_every < 1) print_every = 1;
    int diag_every = Nt / 2000;
    if (diag_every < 1) diag_every = 1;

    /* Output file */
    char path[600];
    snprintf(path, sizeof(path), "%s/m1_%s.tsv", outdir, label);
    FILE *fout = fopen(path, "w");
    if (!fout) { fprintf(stderr, "Cannot open %s\n", path); return; }
    fprintf(fout, "time\tE_total\tE_kin\tE_grad\tE_mass\tE_pot\tE_pw\t"
                  "fc\tpeak_P\tPz\tl2_frac\n");

    compute_acc();

    double wall_start = omp_get_wtime();

    /* Initial diagnostics */
    Diag d0 = compute_diag(core_radius);
    printf("  t=%7.1f  E=%.2f  fc=%.4f  |P|=%.6f  Pz=%.4f  l2=%.3f\n",
           0.0, d0.Et, d0.fc, d0.peak_P, d0.Pz, d0.l2_frac);
    fprintf(fout, "%.4f\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6f\t%.6e\t%.6e\t%.6f\n",
            0.0, d0.Et, d0.Ek, d0.Eg, d0.Em, d0.Ep, d0.Epw,
            d0.fc, d0.peak_P, d0.Pz, d0.l2_frac);

    for (int n = 1; n <= Nt; n++) {
        verlet_step();
        double t = n * dt;

        int do_diag  = (n % diag_every == 0);
        int do_print = (n % print_every == 0) || (n == Nt);

        if (do_diag || do_print) {
            Diag d = compute_diag(core_radius);

            if (do_diag)
                fprintf(fout, "%.4f\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t"
                              "%.6f\t%.6e\t%.6e\t%.6f\n",
                        t, d.Et, d.Ek, d.Eg, d.Em, d.Ep, d.Epw,
                        d.fc, d.peak_P, d.Pz, d.l2_frac);

            if (do_print) {
                double elapsed = omp_get_wtime() - wall_start;
                double frac = (double)n / Nt;
                double eta_t = (frac > 0.001) ? elapsed * (1.0-frac)/frac : 0;
                printf("  t=%7.1f  E=%.2f  fc=%.4f  |P|=%.6f  Pz=%.4f  l2=%.3f  [%.0fs, ETA %.0fs]\n",
                       t, d.Et, d.fc, d.peak_P, d.Pz, d.l2_frac, elapsed, eta_t);
                fflush(stdout);
            }
        }
    }

    fclose(fout);

    /* Final state summary */
    Diag dfinal = compute_diag(core_radius);
    double elapsed1 = omp_get_wtime() - wall_start;

    printf("\n  === %s FINAL (%.1f sec) ===\n", label, elapsed1);
    printf("  E=%.2f  fc=%.4f  |P|=%.6f  Pz=%.4f  l2=%.3f\n",
           dfinal.Et, dfinal.fc, dfinal.peak_P, dfinal.Pz, dfinal.l2_frac);
    printf("  PASS criteria: fc > 0.5 AND |P| > 0.1\n");
    printf("  fc > 0.5? %s    |P| > 0.1? %s\n",
           dfinal.fc > 0.5 ? "YES" : "NO",
           dfinal.peak_P > 0.1 ? "YES" : "NO");
    printf("  OVERALL: %s\n",
           (dfinal.fc > 0.5 && dfinal.peak_P > 0.1) ? "*** PASS ***" : "FAIL");
    fflush(stdout);
}

/* --- Main --- */
int main(int argc, char **argv)
{
    /* Parse args */
    int test_mode = 0;  /* 0=all, 1=M1a, 2=M1b, 3=M1c */
    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "-test") && i+1 < argc) {
            if (!strcmp(argv[i+1], "M1a") || !strcmp(argv[i+1], "m1a")) test_mode = 1;
            else if (!strcmp(argv[i+1], "M1b") || !strcmp(argv[i+1], "m1b")) test_mode = 2;
            else if (!strcmp(argv[i+1], "M1c") || !strcmp(argv[i+1], "m1c")) test_mode = 3;
            i++;
        }
        else if (!strcmp(argv[i], "-mu") && i+1 < argc)     { mu_pot = atof(argv[++i]); }
        else if (!strcmp(argv[i], "-kappa") && i+1 < argc)   { kappa = atof(argv[++i]); }
        else if (!strcmp(argv[i], "-mass") && i+1 < argc)    { mass = atof(argv[++i]); }
        else if (!strcmp(argv[i], "-A") && i+1 < argc)       { A0 = atof(argv[++i]); }
        else if (!strcmp(argv[i], "-N") && i+1 < argc)       { N = atoi(argv[++i]); }
        else if (!strcmp(argv[i], "-L") && i+1 < argc)       { L = atof(argv[++i]); }
        else if (!strcmp(argv[i], "-tfinal") && i+1 < argc)  { tfinal = atof(argv[++i]); }
        else if (!strcmp(argv[i], "-cfl") && i+1 < argc)     { cfl_frac = atof(argv[++i]); }
        else if (!strcmp(argv[i], "-Rtube") && i+1 < argc)   { R_tube = atof(argv[++i]); }
        else if (!strcmp(argv[i], "-o") && i+1 < argc)       { strncpy(outdir, argv[++i], sizeof(outdir)-1); }
        else { fprintf(stderr, "Unknown flag: %s\n", argv[i]); exit(1); }
    }

    /* Set grid for max domain (L=20) for allocation */
    double L_max = 20.0;
    dx  = 2.0 * L_max / (N - 1);
    dx2 = dx * dx;
    m2  = mass * mass;
    dt  = cfl_frac * dx;
    Ngrid = (long)N * N * N;

    printf("=== V27-M1: Stable Propagating Braid ===\n");
    printf("Parameters: mu=%.1f kappa=%.1f mass=%.3f A0=%.3f R_tube=%.1f\n",
           mu_pot, kappa, mass, A0, R_tube);
    printf("Grid: N=%d  Ngrid=%ld (%.1f M)  Memory: %.1f MB\n",
           N, Ngrid, Ngrid/1e6, Ngrid*8.0*10/1e6);
    printf("Threads: %d\n", omp_get_max_threads());
    printf("tfinal=%.0f  cfl=%.2f\n", tfinal, cfl_frac);
    if (test_mode == 0) printf("Running: ALL tests (M1a + M1b + M1c)\n");
    else if (test_mode == 1) printf("Running: M1a only (pairwise scan)\n");
    else if (test_mode == 2) printf("Running: M1b only (small domain)\n");
    else if (test_mode == 3) printf("Running: M1c only (k scan)\n");
    fflush(stdout);

    /* Allocate */
    for (int a = 0; a < 3; a++) {
        phi[a] = calloc(Ngrid, sizeof(double));
        vel[a] = calloc(Ngrid, sizeof(double));
        acc[a] = calloc(Ngrid, sizeof(double));
        if (!phi[a] || !vel[a] || !acc[a]) {
            fprintf(stderr, "Allocation failed\n");
            return 1;
        }
    }
    damp = malloc(Ngrid * sizeof(double));
    if (!damp) { fprintf(stderr, "Allocation failed for damp\n"); return 1; }

    /* DynA baseline: k = pi/L (0.5 twists per domain = 1 twist per period 2L) */
    double k_base = M_PI / 20.0;  /* pi/L for L=20 */

    /* ===== M1a: Pairwise coupling scan (L=20) ===== */
    if (test_mode == 0 || test_mode == 1) {
        printf("\n\n########################################\n");
        printf("# M1a: Pairwise coupling scan (L=20)\n");
        printf("########################################\n");

        double lpw_vals[] = {0.0, 0.3, 0.5, 0.7};
        int n_lpw = 4;
        for (int t = 0; t < n_lpw; t++) {
            char label[64];
            snprintf(label, sizeof(label), "M1a_lpw%.1f", lpw_vals[t]);
            run_test(label, k_base, 20.0, lpw_vals[t]);
        }
    }

    /* ===== M1b: Small domain L=10 ===== */
    if (test_mode == 0 || test_mode == 2) {
        printf("\n\n########################################\n");
        printf("# M1b: Small domain L=10\n");
        printf("########################################\n");

        /* k = pi/L for L=10: same physical wavelength = 1 twist per 2L=20 */
        double k_b = M_PI / 10.0;
        run_test("M1b_L10", k_b, 10.0, 0.0);
    }

    /* ===== M1c: Propagation speed scan ===== */
    if (test_mode == 0 || test_mode == 3) {
        printf("\n\n########################################\n");
        printf("# M1c: k scan (propagation speed)\n");
        printf("########################################\n");

        /* Best lambda_pw from M1a — use 0.5 as default */
        double lpw_best = 0.5;

        /* k values: pi/20, 2pi/20, 4pi/20, 6pi/20 */
        double k_vals[] = {M_PI/20.0, 2.0*M_PI/20.0, 4.0*M_PI/20.0, 6.0*M_PI/20.0};
        const char *k_labels[] = {"M1c_k0.5", "M1c_k1.0", "M1c_k2.0", "M1c_k3.0"};
        int n_k = 4;
        for (int t = 0; t < n_k; t++) {
            run_test(k_labels[t], k_vals[t], 20.0, lpw_best);
        }
    }

    /* Summary */
    printf("\n\n========================================\n");
    printf("=== V27-M1 Complete ===\n");
    printf("========================================\n");
    printf("Results in %s/m1_*.tsv\n", outdir);

    /* Cleanup */
    for (int a = 0; a < 3; a++) {
        free(phi[a]); free(vel[a]); free(acc[a]);
    }
    free(damp);

    return 0;
}
