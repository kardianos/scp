/*
 * triad3d.c — 3D three-body oscillon: three massive scalars + triple-product coupling
 *
 * Lagrangian:
 *   L = Σ_a [ ½(∂_t φ_a)² - ½|∇φ_a|² - ½m²φ_a² ] - V(P)
 *   V = (μ/2) P² / (1 + κ P²),   P = φ₁ φ₂ φ₃
 *
 * EOM: ∂²_t φ_a = ∇²φ_a - m²φ_a - ∂V/∂φ_a
 *
 * Integration: velocity Verlet (symplectic)
 * Boundary: absorbing spherical shell in outer region
 * Parallelism: OpenMP over grid loops
 *
 * Tests:
 *   1: Symmetric triad (A,A,A) — should form oscillon
 *   2: Asymmetric triad (A, 0.8A, 0.6A) — test equalization
 *   3: Two-field control (A,A,0) — P=0, should disperse
 *   4: Single-field control (A,0,0) — free massive, disperses
 *
 * Compile: gcc -O3 -fopenmp -Wall -o triad3d v21/src/triad3d.c -lm
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

/* ─── Parameters ─── */
static double mu       = -10.0;
static double kappa    = 10.0;
static double mass     = 1.0;
static double A_init   = 0.8;
static double sigma    = 3.0;
static int    N        = 200;     /* grid points per side */
static double L        = 30.0;    /* half-width: box is [-L, L]^3 */
static double tfinal   = 5000.0;
static int    test     = 1;
static double cfl_frac = 0.25;    /* dt = cfl_frac * dx */
static char   outdir[512] = "v21/data";

static void parse_args(int argc, char **argv)
{
    for (int i = 1; i < argc - 1; i += 2) {
        if      (!strcmp(argv[i], "-mu"))     mu       = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-kappa"))  kappa    = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-mass"))   mass     = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-A"))      A_init   = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-sigma"))  sigma    = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-N"))      N        = atoi(argv[i+1]);
        else if (!strcmp(argv[i], "-L"))      L        = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-tfinal")) tfinal   = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-test"))   test     = atoi(argv[i+1]);
        else if (!strcmp(argv[i], "-cfl"))    cfl_frac = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-o"))      strncpy(outdir, argv[i+1], sizeof(outdir)-1);
        else { fprintf(stderr, "Unknown flag: %s\n", argv[i]); exit(1); }
    }
}

/* ─── Index helpers ─── */
#define IDX(i,j,k) ((long)(i)*N*N + (long)(j)*N + (long)(k))

/* ─── Main ─── */
int main(int argc, char **argv)
{
    parse_args(argc, argv);

    double dx  = 2.0 * L / (N - 1);
    double dx2 = dx * dx;
    double m2  = mass * mass;
    double dt  = cfl_frac * dx;
    long   Ngrid = (long)N * N * N;
    int    Nt  = (int)(tfinal / dt) + 1;

    /* Initialize amplitudes per test */
    double Ainit[3] = {0, 0, 0};
    const char *test_desc;
    switch (test) {
        case 1: Ainit[0]=Ainit[1]=Ainit[2]=A_init;
                test_desc = "Symmetric triad (A,A,A)"; break;
        case 2: Ainit[0]=A_init; Ainit[1]=0.8*A_init; Ainit[2]=0.6*A_init;
                test_desc = "Asymmetric triad (A,0.8A,0.6A)"; break;
        case 3: Ainit[0]=A_init; Ainit[1]=A_init; Ainit[2]=0;
                test_desc = "Two-field control (A,A,0)"; break;
        case 4: Ainit[0]=A_init; Ainit[1]=0; Ainit[2]=0;
                test_desc = "Single-field control (A,0,0)"; break;
        default: fprintf(stderr, "Unknown test %d\n", test); return 1;
    }

    printf("triad3d test %d: %s\n", test, test_desc);
    printf("  mu=%.3f kappa=%.4f mass=%.4f A=%.3f sigma=%.3f\n",
           mu, kappa, mass, A_init, sigma);
    printf("  N=%d L=%.1f dx=%.5f dt=%.6f cfl=%.3f\n", N, L, dx, dt, cfl_frac);
    printf("  Ngrid=%ld (%.1f M)  Nt=%d  tfinal=%.0f\n",
           Ngrid, Ngrid/1e6, Nt, tfinal);
    printf("  Memory: %.1f MB per field array, %.1f MB total\n",
           Ngrid*8.0/1e6, Ngrid*8.0*10/1e6);  /* 3*(phi+vel+acc) + damp = 10 */
    printf("  Threads: %d\n", omp_get_max_threads());
    fflush(stdout);

    /* Allocate */
    double *phi[3], *vel[3], *acc[3];
    for (int a = 0; a < 3; a++) {
        phi[a] = calloc(Ngrid, sizeof(double));
        vel[a] = calloc(Ngrid, sizeof(double));
        acc[a] = calloc(Ngrid, sizeof(double));
        if (!phi[a] || !vel[a] || !acc[a]) {
            fprintf(stderr, "Allocation failed for field %d\n", a);
            return 1;
        }
    }

    /* Absorbing boundary: spherical damping layer */
    double *damp = malloc(Ngrid * sizeof(double));
    if (!damp) { fprintf(stderr, "Allocation failed for damp\n"); return 1; }

    double R_abs_inner = L * 0.70;  /* absorbing starts at 70% of half-width */
    double R_abs_outer = L * 0.95;  /* full damping at 95% */

    #pragma omp parallel for schedule(static)
    for (long idx = 0; idx < Ngrid; idx++) {
        int i = idx / (N * N);
        int j = (idx / N) % N;
        int k = idx % N;
        double x = -L + i * dx;
        double y = -L + j * dx;
        double z = -L + k * dx;
        double r = sqrt(x*x + y*y + z*z);
        if (r > R_abs_inner) {
            double f = (r - R_abs_inner) / (R_abs_outer - R_abs_inner);
            if (f > 1.0) f = 1.0;
            damp[idx] = 1.0 - 0.98 * f * f;  /* quadratic ramp to 0.02 */
        } else {
            damp[idx] = 1.0;
        }
    }

    /* Initialize: spherical Gaussians */
    #pragma omp parallel for schedule(static)
    for (long idx = 0; idx < Ngrid; idx++) {
        int i = idx / (N * N);
        int j = (idx / N) % N;
        int k = idx % N;
        double x = -L + i * dx;
        double y = -L + j * dx;
        double z = -L + k * dx;
        double r2 = x*x + y*y + z*z;
        double g = exp(-r2 / (2.0 * sigma * sigma));
        for (int a = 0; a < 3; a++)
            phi[a][idx] = Ainit[a] * g;
    }

    /* ─── Compute acceleration ─── */
    /* acc_a = ∇²φ_a - m²φ_a - ∂V/∂φ_a */
    /* ∂V/∂φ_a = μ P (∂P/∂φ_a) / (1+κP²)²  where P = φ₁φ₂φ₃ */
    void compute_acc(void) {
        for (int a = 0; a < 3; a++) {
            #pragma omp parallel for schedule(static)
            for (int i = 1; i < N-1; i++) {
                for (int j = 1; j < N-1; j++) {
                    for (int k = 1; k < N-1; k++) {
                        long idx = IDX(i,j,k);

                        /* 7-point Laplacian */
                        double lapl = (phi[a][IDX(i+1,j,k)] + phi[a][IDX(i-1,j,k)]
                                     + phi[a][IDX(i,j+1,k)] + phi[a][IDX(i,j-1,k)]
                                     + phi[a][IDX(i,j,k+1)] + phi[a][IDX(i,j,k-1)]
                                     - 6.0 * phi[a][idx]) / dx2;

                        /* Triple-product force */
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
                        double dVdphi = mu * P * dP / denom2;

                        acc[a][idx] = lapl - m2 * phi[a][idx] - dVdphi;
                    }
                }
            }

            /* Boundary: zero acceleration */
            #pragma omp parallel for schedule(static)
            for (long idx = 0; idx < Ngrid; idx++) {
                int i = idx / (N * N);
                int j = (idx / N) % N;
                int k = idx % N;
                if (i == 0 || i == N-1 || j == 0 || j == N-1 || k == 0 || k == N-1)
                    acc[a][idx] = 0.0;
            }
        }
    }

    compute_acc();

    /* ─── Output setup ─── */
    char tspath[600];
    snprintf(tspath, sizeof(tspath), "%s/triad3d_test%d_ts.tsv", outdir, test);
    FILE *fts = fopen(tspath, "w");
    if (!fts) { fprintf(stderr, "Cannot open %s\n", tspath); return 1; }
    fprintf(fts, "time\tphi1_0\tphi2_0\tphi3_0\tpeak1\tpeak2\tpeak3\t"
                 "E_kin\tE_grad\tE_mass\tE_pot\tE_total\tf_core\n");

    /* DFT storage for center value */
    int max_dft = 100000;
    double *phi0_hist = malloc(max_dft * sizeof(double));
    double *t_hist = malloc(max_dft * sizeof(double));
    int n_dft = 0;

    /* Diagnostic intervals */
    int diag_every = Nt / 5000;
    if (diag_every < 1) diag_every = 1;
    int print_every = Nt / 100;
    if (print_every < 1) print_every = 1;
    int dft_every = Nt / max_dft;
    if (dft_every < 1) dft_every = 1;

    /* Radial profile output times */
    double profile_times[] = {0, 100, 500, 1000, 2000, 3000, 4000, 5000};
    int n_profile_times = sizeof(profile_times)/sizeof(profile_times[0]);
    int next_profile = 0;

    double core_radius = 3.0 * sigma;
    long ic = IDX(N/2, N/2, N/2);  /* center index */

    /* ─── Time integration ─── */
    double wall_start = omp_get_wtime();

    for (int n = 0; n <= Nt; n++) {
        double t = n * dt;

        /* Record center for DFT */
        if (n % dft_every == 0 && n_dft < max_dft) {
            phi0_hist[n_dft] = phi[0][ic];
            t_hist[n_dft] = t;
            n_dft++;
        }

        int do_diag  = (n % diag_every == 0);
        int do_print = (n % print_every == 0);

        /* Radial profile snapshot */
        if (next_profile < n_profile_times && t >= profile_times[next_profile] - dt/2) {
            char profpath[600];
            snprintf(profpath, sizeof(profpath), "%s/triad3d_test%d_profile_t%04d.tsv",
                     outdir, test, (int)profile_times[next_profile]);
            FILE *fp = fopen(profpath, "w");
            if (fp) {
                fprintf(fp, "r\tphi1\tphi2\tphi3\trho_E\n");
                int nr = (int)(L / dx);
                double *sum_phi = calloc(3 * nr, sizeof(double));
                double *sum_rho = calloc(nr, sizeof(double));
                int *count = calloc(nr, sizeof(int));

                #pragma omp parallel
                {
                    double *lsum_phi = calloc(3 * nr, sizeof(double));
                    double *lsum_rho = calloc(nr, sizeof(double));
                    int *lcount = calloc(nr, sizeof(int));

                    #pragma omp for schedule(static) nowait
                    for (int i = 1; i < N-1; i++) {
                        double x = -L + i * dx;
                        for (int j = 1; j < N-1; j++) {
                            double y = -L + j * dx;
                            for (int k = 1; k < N-1; k++) {
                                double z = -L + k * dx;
                                double r = sqrt(x*x + y*y + z*z);
                                int ir = (int)(r / dx);
                                if (ir >= nr) continue;
                                long idx = IDX(i,j,k);

                                for (int a = 0; a < 3; a++)
                                    lsum_phi[a*nr + ir] += phi[a][idx];

                                /* energy density */
                                double e = 0;
                                for (int a = 0; a < 3; a++) {
                                    e += 0.5 * vel[a][idx] * vel[a][idx];
                                    double gx = (phi[a][IDX(i+1,j,k)] - phi[a][IDX(i-1,j,k)]) / (2*dx);
                                    double gy = (phi[a][IDX(i,j+1,k)] - phi[a][IDX(i,j-1,k)]) / (2*dx);
                                    double gz = (phi[a][IDX(i,j,k+1)] - phi[a][IDX(i,j,k-1)]) / (2*dx);
                                    e += 0.5 * (gx*gx + gy*gy + gz*gz);
                                    e += 0.5 * m2 * phi[a][idx] * phi[a][idx];
                                }
                                double P = phi[0][idx] * phi[1][idx] * phi[2][idx];
                                double P2 = P * P;
                                e += 0.5 * mu * P2 / (1.0 + kappa * P2);
                                lsum_rho[ir] += e;
                                lcount[ir]++;
                            }
                        }
                    }

                    #pragma omp critical
                    {
                        for (int ir = 0; ir < nr; ir++) {
                            for (int a = 0; a < 3; a++)
                                sum_phi[a*nr + ir] += lsum_phi[a*nr + ir];
                            sum_rho[ir] += lsum_rho[ir];
                            count[ir] += lcount[ir];
                        }
                    }
                    free(lsum_phi); free(lsum_rho); free(lcount);
                }

                for (int ir = 0; ir < nr; ir++) {
                    if (count[ir] == 0) continue;
                    double r = (ir + 0.5) * dx;
                    fprintf(fp, "%.6f\t%.6e\t%.6e\t%.6e\t%.6e\n",
                            r,
                            sum_phi[0*nr+ir] / count[ir],
                            sum_phi[1*nr+ir] / count[ir],
                            sum_phi[2*nr+ir] / count[ir],
                            sum_rho[ir] / count[ir]);
                }
                fclose(fp);
                free(sum_phi); free(sum_rho); free(count);
                printf("  [profile written: t=%.0f -> %s]\n", t, profpath);
            }
            next_profile++;
        }

        /* Full diagnostics */
        if (do_diag || do_print) {
            double Ek = 0, Eg = 0, Em = 0, Ep = 0;
            double Ecore = 0, Eall = 0;
            double peak[3] = {0, 0, 0};

            #pragma omp parallel
            {
                double lEk=0, lEg=0, lEm=0, lEp=0, lEc=0, lEa=0;
                double lpk[3] = {0,0,0};

                #pragma omp for schedule(static) nowait
                for (int i = 1; i < N-1; i++) {
                    double x = -L + i * dx;
                    for (int j = 1; j < N-1; j++) {
                        double y = -L + j * dx;
                        for (int k = 1; k < N-1; k++) {
                            double z = -L + k * dx;
                            long idx = IDX(i,j,k);
                            double dV = dx * dx * dx;

                            double e_loc = 0;
                            for (int a = 0; a < 3; a++) {
                                double v2 = vel[a][idx] * vel[a][idx];
                                lEk += 0.5 * v2 * dV;
                                e_loc += 0.5 * v2;

                                double gx = (phi[a][IDX(i+1,j,k)] - phi[a][IDX(i-1,j,k)]) / (2*dx);
                                double gy = (phi[a][IDX(i,j+1,k)] - phi[a][IDX(i,j-1,k)]) / (2*dx);
                                double gz = (phi[a][IDX(i,j,k+1)] - phi[a][IDX(i,j,k-1)]) / (2*dx);
                                double grad2 = gx*gx + gy*gy + gz*gz;
                                lEg += 0.5 * grad2 * dV;
                                e_loc += 0.5 * grad2;

                                double mass_e = 0.5 * m2 * phi[a][idx] * phi[a][idx];
                                lEm += mass_e * dV;
                                e_loc += mass_e;

                                double ap = fabs(phi[a][idx]);
                                if (ap > lpk[a]) lpk[a] = ap;
                            }

                            double P = phi[0][idx] * phi[1][idx] * phi[2][idx];
                            double P2 = P * P;
                            double Vloc = 0.5 * mu * P2 / (1.0 + kappa * P2);
                            lEp += Vloc * dV;
                            e_loc += Vloc;

                            lEa += e_loc * dV;
                            double r = sqrt(x*x + y*y + z*z);
                            if (r < core_radius)
                                lEc += e_loc * dV;
                        }
                    }
                }

                #pragma omp critical
                {
                    Ek += lEk; Eg += lEg; Em += lEm; Ep += lEp;
                    Ecore += lEc; Eall += lEa;
                    for (int a = 0; a < 3; a++)
                        if (lpk[a] > peak[a]) peak[a] = lpk[a];
                }
            }

            double Et = Ek + Eg + Em + Ep;
            double fc = (Eall > 1e-20) ? Ecore / Eall : 0.0;

            if (do_diag)
                fprintf(fts, "%.6f\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t"
                             "%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6f\n",
                        t, phi[0][ic], phi[1][ic], phi[2][ic],
                        peak[0], peak[1], peak[2], Ek, Eg, Em, Ep, Et, fc);

            if (do_print) {
                double wall_now = omp_get_wtime();
                double elapsed = wall_now - wall_start;
                double frac = (double)n / Nt;
                double eta = (frac > 0.001) ? elapsed * (1.0 - frac) / frac : 0;
                printf("  t=%7.1f  phi0=(%+.3f,%+.3f,%+.3f)  pk=(%.3f,%.3f,%.3f)  "
                       "E=%+.2f  Ep=%+.2f  fc=%.4f  [%.0fs elapsed, ETA %.0fs]\n",
                       t, phi[0][ic], phi[1][ic], phi[2][ic],
                       peak[0], peak[1], peak[2], Et, Ep, fc,
                       elapsed, eta);
                fflush(stdout);
            }
        }

        if (n == Nt) break;

        /* ─── Velocity Verlet ─── */

        /* 1. Half-kick */
        for (int a = 0; a < 3; a++) {
            #pragma omp parallel for schedule(static)
            for (long idx = 0; idx < Ngrid; idx++)
                vel[a][idx] += 0.5 * dt * acc[a][idx];
        }

        /* 2. Drift */
        for (int a = 0; a < 3; a++) {
            #pragma omp parallel for schedule(static)
            for (long idx = 0; idx < Ngrid; idx++)
                phi[a][idx] += dt * vel[a][idx];
        }

        /* 3. Recompute acceleration */
        compute_acc();

        /* 4. Half-kick */
        for (int a = 0; a < 3; a++) {
            #pragma omp parallel for schedule(static)
            for (long idx = 0; idx < Ngrid; idx++)
                vel[a][idx] += 0.5 * dt * acc[a][idx];
        }

        /* 5. Absorbing boundary damping */
        for (int a = 0; a < 3; a++) {
            #pragma omp parallel for schedule(static)
            for (long idx = 0; idx < Ngrid; idx++) {
                vel[a][idx] *= damp[idx];
                phi[a][idx] *= damp[idx];
            }
        }
    }

    double wall_total = omp_get_wtime() - wall_start;
    fclose(fts);

    /* ─── DFT of φ₁(origin, t) ─── */
    int dft_start = n_dft / 2;
    double peak_om = 0;
    if (n_dft - dft_start > 100) {
        char dftpath[600];
        snprintf(dftpath, sizeof(dftpath), "%s/triad3d_test%d_spectrum.tsv", outdir, test);
        FILE *fdft = fopen(dftpath, "w");
        fprintf(fdft, "omega\tpower\n");
        double T = t_hist[n_dft-1] - t_hist[dft_start];
        int nf = 500;
        double peak_pow = 0;
        for (int kk = 0; kk < nf; kk++) {
            double omega = 3.0 * mass * kk / nf;
            double re = 0, im = 0;
            for (int j = dft_start; j < n_dft; j++) {
                double dtj = (j > dft_start) ?
                    (t_hist[j]-t_hist[j-1]) : (t_hist[dft_start+1]-t_hist[dft_start]);
                re += phi0_hist[j] * cos(omega * t_hist[j]) * dtj;
                im += phi0_hist[j] * sin(omega * t_hist[j]) * dtj;
            }
            double pw = (re*re + im*im) / (T*T);
            fprintf(fdft, "%.6f\t%.6e\n", omega, pw);
            if (pw > peak_pow) { peak_pow = pw; peak_om = omega; }
        }
        fclose(fdft);
        printf("\nSpectrum: peak omega = %.4f (mass gap = %.4f)\n", peak_om, mass);
        printf("Oscillon (omega < m)? %s\n",
               (peak_om > 0.01 && peak_om < mass) ? "YES" : "NO");
    }

    printf("\n=== SUMMARY ===\n");
    printf("  Test %d: %s\n", test, test_desc);
    printf("  Wall time: %.1f seconds (%.1f min)\n", wall_total, wall_total/60);
    printf("  Output: %s\n", tspath);
    printf("===\n");

    /* Cleanup */
    for (int a = 0; a < 3; a++) {
        free(phi[a]); free(vel[a]); free(acc[a]);
    }
    free(damp); free(phi0_hist); free(t_hist);
    return 0;
}
