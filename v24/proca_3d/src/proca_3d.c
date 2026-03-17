/*
 * proca_3d.c — 3D oscillon with triple-product + pairwise coupling
 *
 * Based on v21/src/triad3d.c. Adds pairwise mass coupling:
 *   V_pw = lambda * (phi1*phi2 + phi2*phi3 + phi3*phi1)
 *   => acc[a] += -lambda*(phi[b]+phi[c])
 *
 * Phase 1: Single oscillon stability scan over lambda
 * Phase 2: Two-oscillon interaction at best lambda
 *
 * Compile: gcc -O3 -fopenmp -Wall -o proca_3d src/proca_3d.c -lm
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

/* ─── Parameters ─── */
static double mu       = -20.0;
static double kappa    = 20.0;
static double mass     = 1.0;
static double A_init   = 0.8;
static double sigma    = 3.0;
static double lambda   = 0.0;   /* pairwise coupling strength */
static int    N        = 96;
static double L        = 15.0;
static double tfinal   = 500.0;
static int    phase    = 1;     /* 1=single, 2=two-oscillon */
static double D_sep    = 20.0;  /* separation for phase 2 */
static double cfl_frac = 0.25;
static char   outdir[512] = "v24/proca_3d/data";

static void parse_args(int argc, char **argv)
{
    for (int i = 1; i < argc - 1; i += 2) {
        if      (!strcmp(argv[i], "-mu"))     mu       = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-kappa"))  kappa    = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-mass"))   mass     = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-A"))      A_init   = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-sigma"))  sigma    = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-lambda")) lambda   = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-N"))      N        = atoi(argv[i+1]);
        else if (!strcmp(argv[i], "-L"))      L        = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-tfinal")) tfinal   = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-phase"))  phase    = atoi(argv[i+1]);
        else if (!strcmp(argv[i], "-D"))      D_sep    = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-cfl"))    cfl_frac = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-o"))      strncpy(outdir, argv[i+1], sizeof(outdir)-1);
        else { fprintf(stderr, "Unknown flag: %s\n", argv[i]); exit(1); }
    }
}

/* ─── Index helpers ─── */
#define IDX(i,j,k) ((long)(i)*N*N + (long)(j)*N + (long)(k))

/* ─── Globals for field arrays ─── */
static double *phi_arr[3], *vel_arr[3], *acc_arr[3];
static double *damp_arr;
static double dx, dx2, m2, dt_g;
static long Ngrid;

static void compute_acc(void)
{
    double m2_loc = m2;
    double mu_loc = mu;
    double kappa_loc = kappa;
    double lam = lambda;
    double dx2_loc = dx2;
    int NN = N;

    for (int a = 0; a < 3; a++) {
        int b = (a + 1) % 3;
        int c = (a + 2) % 3;

        #pragma omp parallel for schedule(static)
        for (int i = 1; i < NN-1; i++) {
            for (int j = 1; j < NN-1; j++) {
                for (int k = 1; k < NN-1; k++) {
                    long idx = IDX(i,j,k);

                    /* 7-point Laplacian */
                    double lapl = (phi_arr[a][IDX(i+1,j,k)] + phi_arr[a][IDX(i-1,j,k)]
                                 + phi_arr[a][IDX(i,j+1,k)] + phi_arr[a][IDX(i,j-1,k)]
                                 + phi_arr[a][IDX(i,j,k+1)] + phi_arr[a][IDX(i,j,k-1)]
                                 - 6.0 * phi_arr[a][idx]) / dx2_loc;

                    /* Triple-product force */
                    double p0 = phi_arr[0][idx], p1 = phi_arr[1][idx], p2 = phi_arr[2][idx];
                    double P = p0 * p1 * p2;
                    double P2 = P * P;
                    double denom2 = (1.0 + kappa_loc * P2);
                    denom2 *= denom2;

                    double dP;
                    switch (a) {
                        case 0: dP = p1 * p2; break;
                        case 1: dP = p0 * p2; break;
                        default: dP = p0 * p1; break;
                    }
                    double dVdphi = mu_loc * P * dP / denom2;

                    /* Pairwise coupling: -lambda*(phi_b + phi_c) */
                    double pw = lam * (phi_arr[b][idx] + phi_arr[c][idx]);

                    acc_arr[a][idx] = lapl - m2_loc * phi_arr[a][idx] - dVdphi - pw;
                }
            }
        }

        /* Boundary: zero acceleration */
        #pragma omp parallel for schedule(static)
        for (long idx = 0; idx < Ngrid; idx++) {
            int i = idx / (N * N);
            int j = (idx / N) % N;
            int k = idx % N;
            if (i == 0 || i == NN-1 || j == 0 || j == NN-1 || k == 0 || k == NN-1)
                acc_arr[a][idx] = 0.0;
        }
    }
}

int main(int argc, char **argv)
{
    parse_args(argc, argv);

    dx  = 2.0 * L / (N - 1);
    dx2 = dx * dx;
    m2  = mass * mass;
    dt_g = cfl_frac * dx;
    double dt = dt_g;
    Ngrid = (long)N * N * N;
    int Nt = (int)(tfinal / dt) + 1;

    printf("proca_3d: Phase %d, lambda=%.4f\n", phase, lambda);
    printf("  mu=%.3f kappa=%.4f mass=%.4f A=%.3f sigma=%.3f\n",
           mu, kappa, mass, A_init, sigma);
    printf("  N=%d L=%.1f dx=%.5f dt=%.6f cfl=%.3f\n", N, L, dx, dt, cfl_frac);
    printf("  Ngrid=%ld (%.1f M)  Nt=%d  tfinal=%.0f\n",
           Ngrid, Ngrid/1e6, Nt, tfinal);
    printf("  Memory: %.1f MB per field array, %.1f MB total\n",
           Ngrid*8.0/1e6, Ngrid*8.0*10/1e6);
    printf("  Threads: %d\n", omp_get_max_threads());
    if (phase == 2)
        printf("  D_sep=%.1f (two-oscillon along z)\n", D_sep);
    fflush(stdout);

    /* Allocate */
    for (int a = 0; a < 3; a++) {
        phi_arr[a] = calloc(Ngrid, sizeof(double));
        vel_arr[a] = calloc(Ngrid, sizeof(double));
        acc_arr[a] = calloc(Ngrid, sizeof(double));
        if (!phi_arr[a] || !vel_arr[a] || !acc_arr[a]) {
            fprintf(stderr, "Allocation failed for field %d\n", a);
            return 1;
        }
    }

    /* Absorbing boundary */
    damp_arr = malloc(Ngrid * sizeof(double));
    if (!damp_arr) { fprintf(stderr, "Allocation failed for damp\n"); return 1; }

    double R_abs_inner = L * 0.70;
    double R_abs_outer = L * 0.95;

    if (phase == 1) {
        /* Spherical absorbing layer centered at origin */
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
                damp_arr[idx] = 1.0 - 0.98 * f * f;
            } else {
                damp_arr[idx] = 1.0;
            }
        }
    } else {
        /* Phase 2: box-edge absorbing (need room for two oscillons along z) */
        #pragma omp parallel for schedule(static)
        for (long idx = 0; idx < Ngrid; idx++) {
            int i = idx / (N * N);
            int j = (idx / N) % N;
            int k = idx % N;
            double x = -L + i * dx;
            double y = -L + j * dx;
            double z = -L + k * dx;
            /* Use distance to nearest box face */
            double d_edge = L - fabs(x);
            double dy = L - fabs(y);
            double dz = L - fabs(z);
            if (dy < d_edge) d_edge = dy;
            if (dz < d_edge) d_edge = dz;
            double w = L * 0.20; /* absorbing width = 20% of L */
            if (d_edge < w) {
                double f = 1.0 - d_edge / w;
                damp_arr[idx] = 1.0 - 0.98 * f * f;
            } else {
                damp_arr[idx] = 1.0;
            }
        }
    }

    /* Initialize: spherical Gaussians */
    if (phase == 1) {
        /* Single oscillon at origin */
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
                phi_arr[a][idx] = A_init * g;
        }
    } else {
        /* Two oscillons along z at ±D/2 */
        #pragma omp parallel for schedule(static)
        for (long idx = 0; idx < Ngrid; idx++) {
            int i = idx / (N * N);
            int j = (idx / N) % N;
            int k = idx % N;
            double x = -L + i * dx;
            double y = -L + j * dx;
            double z = -L + k * dx;
            double z1 = z - D_sep/2.0;
            double z2 = z + D_sep/2.0;
            double r2_1 = x*x + y*y + z1*z1;
            double r2_2 = x*x + y*y + z2*z2;
            double g1 = exp(-r2_1 / (2.0 * sigma * sigma));
            double g2 = exp(-r2_2 / (2.0 * sigma * sigma));
            for (int a = 0; a < 3; a++)
                phi_arr[a][idx] = A_init * (g1 + g2);
        }
    }

    compute_acc();

    /* ─── Output setup ─── */
    char tspath[600];
    if (phase == 1)
        snprintf(tspath, sizeof(tspath), "%s/proca3d_lam%.2f_ts.tsv", outdir, lambda);
    else
        snprintf(tspath, sizeof(tspath), "%s/proca3d_lam%.2f_D%.0f_ts.tsv",
                 outdir, lambda, D_sep);

    FILE *fts = fopen(tspath, "w");
    if (!fts) { fprintf(stderr, "Cannot open %s\n", tspath); return 1; }

    if (phase == 1) {
        fprintf(fts, "time\tphi1_0\tphi2_0\tphi3_0\tpeak1\tpeak2\tpeak3\t"
                     "E_kin\tE_grad\tE_mass\tE_pot\tE_pw\tE_total\tf_core\n");
    } else {
        fprintf(fts, "time\tpeak1\tpeak2\tpeak3\tE_total\tf_core\t"
                     "z_com1\tz_com2\tsep\n");
    }

    /* DFT storage for center value */
    int max_dft = 50000;
    double *phi0_hist = malloc(max_dft * sizeof(double));
    double *t_hist = malloc(max_dft * sizeof(double));
    int n_dft = 0;

    /* Diagnostic intervals */
    int diag_every = Nt / 5000;
    if (diag_every < 1) diag_every = 1;
    int print_every = Nt / 50;
    if (print_every < 1) print_every = 1;
    int dft_every = Nt / max_dft;
    if (dft_every < 1) dft_every = 1;

    double core_radius = 3.0 * sigma;
    long ic = IDX(N/2, N/2, N/2);

    /* ─── Time integration ─── */
    double wall_start = omp_get_wtime();

    for (int n = 0; n <= Nt; n++) {
        double t = n * dt;

        /* Record center for DFT */
        if (n % dft_every == 0 && n_dft < max_dft) {
            phi0_hist[n_dft] = phi_arr[0][ic];
            t_hist[n_dft] = t;
            n_dft++;
        }

        int do_diag  = (n % diag_every == 0);
        int do_print = (n % print_every == 0);

        if (do_diag || do_print) {
            double Ek = 0, Eg = 0, Em = 0, Ep = 0, Epw = 0;
            double Ecore = 0, Eall = 0;
            double peak[3] = {0, 0, 0};

            /* For phase 2: COM tracking */
            double z_w1 = 0, z_w2 = 0, w1 = 0, w2 = 0;

            #pragma omp parallel
            {
                double lEk=0, lEg=0, lEm=0, lEp=0, lEpw=0, lEc=0, lEa=0;
                double lpk[3] = {0,0,0};
                double lzw1=0, lzw2=0, lw1=0, lw2=0;

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
                                double v2 = vel_arr[a][idx] * vel_arr[a][idx];
                                lEk += 0.5 * v2 * dV;
                                e_loc += 0.5 * v2;

                                double gx = (phi_arr[a][IDX(i+1,j,k)] - phi_arr[a][IDX(i-1,j,k)]) / (2*dx);
                                double gy = (phi_arr[a][IDX(i,j+1,k)] - phi_arr[a][IDX(i,j-1,k)]) / (2*dx);
                                double gz = (phi_arr[a][IDX(i,j,k+1)] - phi_arr[a][IDX(i,j,k-1)]) / (2*dx);
                                double grad2 = gx*gx + gy*gy + gz*gz;
                                lEg += 0.5 * grad2 * dV;
                                e_loc += 0.5 * grad2;

                                double mass_e = 0.5 * m2 * phi_arr[a][idx] * phi_arr[a][idx];
                                lEm += mass_e * dV;
                                e_loc += mass_e;

                                double ap = fabs(phi_arr[a][idx]);
                                if (ap > lpk[a]) lpk[a] = ap;
                            }

                            double PP = phi_arr[0][idx] * phi_arr[1][idx] * phi_arr[2][idx];
                            double P2 = PP * PP;
                            double Vloc = 0.5 * mu * P2 / (1.0 + kappa * P2);
                            lEp += Vloc * dV;
                            e_loc += Vloc;

                            /* Pairwise energy: lambda*(phi0*phi1 + phi1*phi2 + phi2*phi0) */
                            double pw_e = lambda * (phi_arr[0][idx]*phi_arr[1][idx]
                                        + phi_arr[1][idx]*phi_arr[2][idx]
                                        + phi_arr[2][idx]*phi_arr[0][idx]);
                            lEpw += pw_e * dV;
                            e_loc += pw_e;

                            lEa += e_loc * dV;
                            double r = sqrt(x*x + y*y + z*z);
                            if (r < core_radius)
                                lEc += e_loc * dV;

                            /* Phase 2: COM of each soliton (split by z>0, z<0) */
                            if (phase == 2) {
                                double rho_e = 0;
                                for (int a = 0; a < 3; a++)
                                    rho_e += phi_arr[a][idx] * phi_arr[a][idx];
                                if (z > 0) {
                                    lzw1 += z * rho_e * dV;
                                    lw1 += rho_e * dV;
                                } else {
                                    lzw2 += z * rho_e * dV;
                                    lw2 += rho_e * dV;
                                }
                            }
                        }
                    }
                }

                #pragma omp critical
                {
                    Ek += lEk; Eg += lEg; Em += lEm; Ep += lEp; Epw += lEpw;
                    Ecore += lEc; Eall += lEa;
                    for (int a = 0; a < 3; a++)
                        if (lpk[a] > peak[a]) peak[a] = lpk[a];
                    z_w1 += lzw1; z_w2 += lzw2; w1 += lw1; w2 += lw2;
                }
            }

            double Et = Ek + Eg + Em + Ep + Epw;
            double fc = (Eall > 1e-20) ? Ecore / Eall : 0.0;

            if (do_diag) {
                if (phase == 1) {
                    fprintf(fts, "%.6f\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t"
                                 "%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6f\n",
                            t, phi_arr[0][ic], phi_arr[1][ic], phi_arr[2][ic],
                            peak[0], peak[1], peak[2], Ek, Eg, Em, Ep, Epw, Et, fc);
                } else {
                    double zc1 = (w1 > 1e-20) ? z_w1 / w1 : D_sep/2;
                    double zc2 = (w2 > 1e-20) ? z_w2 / w2 : -D_sep/2;
                    double sep = zc1 - zc2;
                    fprintf(fts, "%.6f\t%.6e\t%.6e\t%.6e\t%.6e\t%.6f\t"
                                 "%.6f\t%.6f\t%.6f\n",
                            t, peak[0], peak[1], peak[2], Et, fc,
                            zc1, zc2, sep);
                }
            }

            if (do_print) {
                double wall_now = omp_get_wtime();
                double elapsed = wall_now - wall_start;
                double frac = (double)n / Nt;
                double eta = (frac > 0.001) ? elapsed * (1.0 - frac) / frac : 0;
                if (phase == 1) {
                    printf("  t=%7.1f  phi0=(%+.3f,%+.3f,%+.3f)  pk=(%.3f,%.3f,%.3f)  "
                           "E=%+.2f  Epw=%+.2f  fc=%.4f  [%.0fs, ETA %.0fs]\n",
                           t, phi_arr[0][ic], phi_arr[1][ic], phi_arr[2][ic],
                           peak[0], peak[1], peak[2], Et, Epw, fc,
                           elapsed, eta);
                } else {
                    double zc1 = (w1 > 1e-20) ? z_w1 / w1 : D_sep/2;
                    double zc2 = (w2 > 1e-20) ? z_w2 / w2 : -D_sep/2;
                    printf("  t=%7.1f  pk=(%.3f,%.3f,%.3f)  E=%+.2f  fc=%.4f  "
                           "sep=%.3f  [%.0fs, ETA %.0fs]\n",
                           t, peak[0], peak[1], peak[2], Et, fc,
                           zc1-zc2, elapsed, eta);
                }
                fflush(stdout);
            }
        }

        if (n == Nt) break;

        /* ─── Velocity Verlet ─── */

        /* 1. Half-kick */
        for (int a = 0; a < 3; a++) {
            #pragma omp parallel for schedule(static)
            for (long idx = 0; idx < Ngrid; idx++)
                vel_arr[a][idx] += 0.5 * dt * acc_arr[a][idx];
        }

        /* 2. Drift */
        for (int a = 0; a < 3; a++) {
            #pragma omp parallel for schedule(static)
            for (long idx = 0; idx < Ngrid; idx++)
                phi_arr[a][idx] += dt * vel_arr[a][idx];
        }

        /* 3. Recompute acceleration */
        compute_acc();

        /* 4. Half-kick */
        for (int a = 0; a < 3; a++) {
            #pragma omp parallel for schedule(static)
            for (long idx = 0; idx < Ngrid; idx++)
                vel_arr[a][idx] += 0.5 * dt * acc_arr[a][idx];
        }

        /* 5. Absorbing boundary damping */
        for (int a = 0; a < 3; a++) {
            #pragma omp parallel for schedule(static)
            for (long idx = 0; idx < Ngrid; idx++) {
                vel_arr[a][idx] *= damp_arr[idx];
                phi_arr[a][idx] *= damp_arr[idx];
            }
        }
    }

    double wall_total = omp_get_wtime() - wall_start;
    fclose(fts);

    /* ─── DFT of phi_1(origin, t) ─── */
    int dft_start = n_dft / 2;
    double peak_om = 0;
    if (phase == 1 && n_dft - dft_start > 100) {
        char dftpath[600];
        snprintf(dftpath, sizeof(dftpath), "%s/proca3d_lam%.2f_spectrum.tsv", outdir, lambda);
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
    printf("  Phase %d, lambda=%.4f\n", phase, lambda);
    printf("  Wall time: %.1f seconds (%.1f min)\n", wall_total, wall_total/60);
    printf("  Output: %s\n", tspath);
    if (phase == 1)
        printf("  Peak omega: %.4f, gap margin: %.1f%%\n",
               peak_om, (1.0 - peak_om/mass)*100);
    printf("===\n");

    /* Cleanup */
    for (int a = 0; a < 3; a++) {
        free(phi_arr[a]); free(vel_arr[a]); free(acc_arr[a]);
    }
    free(damp_arr); free(phi0_hist); free(t_hist);
    return 0;
}
