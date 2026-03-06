/*
 * triad3d_grav.c — 3D three-body oscillon + massless scalar gravity (Φ)
 *
 * Physics:
 *   □Φ = -α ρ_source              (Φ sourced by oscillon energy density)
 *   m_eff² = m₀² - β Φ            (Φ modulates effective mass)
 *   ∂²φ_a/∂t² = ∇²φ_a - m_eff²·φ_a - ∂V/∂φ_a
 *
 * Sign conventions: Φ < 0 near source (attractive), -βΦ > 0 → deeper mass gap
 *
 * Φ is initialized from the Poisson solution (∇²Φ = αρ) to avoid transients.
 * ρ_source = FULL energy density (kinetic + gradient + mass + potential).
 * Since E_total is conserved, ∫ρ d³x = const → no monopole radiation.
 *
 * Compile: gcc -O3 -fopenmp -Wall -o triad3d_grav v21/src/triad3d_grav.c -lm
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

/* Parameters */
static double mu       = -20.0;
static double kappa    = 20.0;
static double mass0    = 1.0;     /* bare mass m₀ */
static double A_init   = 0.8;
static double sigma    = 3.0;
static int    N        = 100;
static double L        = 20.0;
static double tfinal   = 500.0;
static double cfl_frac = 0.25;
static char   outdir[512] = "v21/data";

/* Gravity coupling */
static double alpha_g  = 0.0;    /* Φ source coupling (0 = no gravity) */
static double beta_g   = 0.0;    /* mass modulation coupling (0 = no backreaction) */

static void parse_args(int argc, char **argv)
{
    for (int i = 1; i < argc - 1; i += 2) {
        if      (!strcmp(argv[i], "-mu"))     mu       = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-kappa"))  kappa    = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-mass"))   mass0    = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-A"))      A_init   = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-sigma"))  sigma    = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-N"))      N        = atoi(argv[i+1]);
        else if (!strcmp(argv[i], "-L"))      L        = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-tfinal")) tfinal   = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-cfl"))    cfl_frac = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-o"))      strncpy(outdir, argv[i+1], sizeof(outdir)-1);
        else if (!strcmp(argv[i], "-alpha"))  alpha_g  = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-beta"))   beta_g   = atof(argv[i+1]);
        else { fprintf(stderr, "Unknown flag: %s\n", argv[i]); exit(1); }
    }
}

#define IDX(i,j,k) ((long)(i)*N*N + (long)(j)*N + (long)(k))

int main(int argc, char **argv)
{
    parse_args(argc, argv);

    double dx  = 2.0 * L / (N - 1);
    double dx2 = dx * dx;
    double m02 = mass0 * mass0;
    double dt  = cfl_frac * dx;
    long   Ngrid = (long)N * N * N;
    int    Nt  = (int)(tfinal / dt) + 1;
    double dV  = dx * dx * dx;

    printf("triad3d_grav: oscillon + scalar gravity\n");
    printf("  mu=%.3f kappa=%.4f mass0=%.4f A=%.3f sigma=%.3f\n",
           mu, kappa, mass0, A_init, sigma);
    printf("  alpha=%.6e beta=%.6e\n", alpha_g, beta_g);
    printf("  N=%d L=%.1f dx=%.5f dt=%.6f\n", N, L, dx, dt);
    printf("  Ngrid=%ld (%.1f M)  Nt=%d  tfinal=%.0f\n",
           Ngrid, Ngrid/1e6, Nt, tfinal);
    printf("  Threads: %d\n", omp_get_max_threads());
    fflush(stdout);

    /* Allocate: 3 scalar fields + gravity field */
    double *phi[3], *vel[3], *acc[3];
    double *phi_g, *vel_g, *acc_g;  /* gravity scalar Φ */
    double *rho_src;                 /* source energy density */

    for (int a = 0; a < 3; a++) {
        phi[a] = calloc(Ngrid, sizeof(double));
        vel[a] = calloc(Ngrid, sizeof(double));
        acc[a] = calloc(Ngrid, sizeof(double));
    }
    phi_g   = calloc(Ngrid, sizeof(double));
    vel_g   = calloc(Ngrid, sizeof(double));
    acc_g   = calloc(Ngrid, sizeof(double));
    rho_src = calloc(Ngrid, sizeof(double));

    /* Absorbing boundary */
    double *damp = malloc(Ngrid * sizeof(double));
    double R_abs_inner = L * 0.70;
    double R_abs_outer = L * 0.95;

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
            damp[idx] = 1.0 - 0.98 * f * f;
        } else {
            damp[idx] = 1.0;
        }
    }

    /* Initialize: spherical Gaussians for all three fields */
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
            phi[a][idx] = A_init * g;
    }

    /* Initialize Φ from Poisson solution of initial ρ.
     * Compute radial ρ, solve ODE, interpolate back to 3D.
     * Uses Green's function: Φ(r) = -α[I1(r)/r + I2(r)]
     * where I1(r) = ∫₀ʳ ρ(r')r'²dr', I2(r) = ∫ᵣ^∞ ρ(r')r'dr'
     */
    if (alpha_g > 0) {
        /* Build radial ρ profile from initial Gaussians */
        int Nr = (int)(L / dx) + 1;
        double *rho_rad = calloc(Nr, sizeof(double));
        int *count = calloc(Nr, sizeof(int));

        #pragma omp parallel
        {
            double *lrho = calloc(Nr, sizeof(double));
            int *lcnt = calloc(Nr, sizeof(int));

            #pragma omp for schedule(static) nowait
            for (int i = 1; i < N-1; i++) {
                double x = -L + i * dx;
                for (int j = 1; j < N-1; j++) {
                    double y = -L + j * dx;
                    for (int k = 1; k < N-1; k++) {
                        double z = -L + k * dx;
                        long idx2 = IDX(i,j,k);
                        double r = sqrt(x*x + y*y + z*z);
                        int ir = (int)(r / dx);
                        if (ir >= Nr) continue;

                        /* ρ_source = kinetic + gradient + mass + potential
                         * (at t=0, vel=0 so kinetic contribution is zero) */
                        double e = 0;
                        for (int a = 0; a < 3; a++) {
                            double gx = (phi[a][IDX(i+1,j,k)] - phi[a][IDX(i-1,j,k)]) / (2*dx);
                            double gy = (phi[a][IDX(i,j+1,k)] - phi[a][IDX(i,j-1,k)]) / (2*dx);
                            double gz = (phi[a][IDX(i,j,k+1)] - phi[a][IDX(i,j,k-1)]) / (2*dx);
                            e += 0.5 * (gx*gx + gy*gy + gz*gz);
                            e += 0.5 * m02 * phi[a][idx2] * phi[a][idx2];
                        }
                        double P = phi[0][idx2] * phi[1][idx2] * phi[2][idx2];
                        double P2 = P * P;
                        e += 0.5 * mu * P2 / (1.0 + kappa * P2);

                        lrho[ir] += e;
                        lcnt[ir]++;
                    }
                }
            }

            #pragma omp critical
            {
                for (int ir = 0; ir < Nr; ir++) {
                    rho_rad[ir] += lrho[ir];
                    count[ir] += lcnt[ir];
                }
            }
            free(lrho); free(lcnt);
        }

        for (int ir = 0; ir < Nr; ir++)
            if (count[ir] > 0) rho_rad[ir] /= count[ir];

        /* Solve Poisson via Green's function */
        double *I1 = calloc(Nr, sizeof(double));
        double *I2 = calloc(Nr, sizeof(double));

        I1[0] = 0;
        for (int ir = 1; ir < Nr; ir++) {
            double r = ir * dx;
            I1[ir] = I1[ir-1] + rho_rad[ir] * r * r * dx;
        }
        I2[Nr-1] = 0;
        for (int ir = Nr-2; ir >= 0; ir--) {
            double r = (ir+1) * dx;
            I2[ir] = I2[ir+1] + rho_rad[ir+1] * r * dx;
        }

        double Phi0 = -alpha_g * I2[0];
        double Q_total = 0;
        for (int ir = 1; ir < Nr; ir++) {
            double r = ir * dx;
            Q_total += rho_rad[ir] * 4.0 * M_PI * r * r * dx;
        }

        printf("  Poisson init: Phi(0) = %.6e, Q = %.4f\n", Phi0, Q_total);
        printf("  Expected 1/r coefficient: %.6e\n", -alpha_g * Q_total / (4.0 * M_PI));

        /* Interpolate Φ(r) to 3D grid */
        #pragma omp parallel for schedule(static)
        for (long idx2 = 0; idx2 < Ngrid; idx2++) {
            int i = idx2 / (N * N);
            int j = (idx2 / N) % N;
            int k = idx2 % N;
            double x = -L + i * dx;
            double y = -L + j * dx;
            double z = -L + k * dx;
            double r = sqrt(x*x + y*y + z*z);
            int ir = (int)(r / dx);

            if (ir <= 0) {
                phi_g[idx2] = Phi0;
            } else if (ir >= Nr - 1) {
                double rr = r;
                if (rr < dx) rr = dx;
                phi_g[idx2] = -alpha_g * I1[Nr-1] / rr;
            } else {
                double rr = r;
                if (rr < 0.01 * dx) rr = 0.01 * dx;
                phi_g[idx2] = -alpha_g * (I1[ir] / rr + I2[ir]);
            }
        }

        free(rho_rad); free(count); free(I1); free(I2);
    }

    /* Compute source ρ and accelerations */
    void compute_acc(void) {
        /* First compute ρ_source at every point (FULL energy density) */
        #pragma omp parallel for schedule(static)
        for (int i = 1; i < N-1; i++) {
            for (int j = 1; j < N-1; j++) {
                for (int k = 1; k < N-1; k++) {
                    long idx = IDX(i,j,k);
                    double e = 0;
                    for (int a = 0; a < 3; a++) {
                        /* Kinetic energy */
                        e += 0.5 * vel[a][idx] * vel[a][idx];
                        /* Gradient energy */
                        double gx = (phi[a][IDX(i+1,j,k)] - phi[a][IDX(i-1,j,k)]) / (2*dx);
                        double gy = (phi[a][IDX(i,j+1,k)] - phi[a][IDX(i,j-1,k)]) / (2*dx);
                        double gz = (phi[a][IDX(i,j,k+1)] - phi[a][IDX(i,j,k-1)]) / (2*dx);
                        e += 0.5 * (gx*gx + gy*gy + gz*gz);
                        /* Mass energy */
                        e += 0.5 * m02 * phi[a][idx] * phi[a][idx];
                    }
                    double P = phi[0][idx] * phi[1][idx] * phi[2][idx];
                    double P2 = P * P;
                    e += 0.5 * mu * P2 / (1.0 + kappa * P2);
                    rho_src[idx] = e;
                }
            }
        }

        /* Φ field acceleration: ∇²Φ - α·ρ */
        if (alpha_g > 0) {
            #pragma omp parallel for schedule(static)
            for (int i = 1; i < N-1; i++) {
                for (int j = 1; j < N-1; j++) {
                    for (int k = 1; k < N-1; k++) {
                        long idx = IDX(i,j,k);
                        double lapl = (phi_g[IDX(i+1,j,k)] + phi_g[IDX(i-1,j,k)]
                                     + phi_g[IDX(i,j+1,k)] + phi_g[IDX(i,j-1,k)]
                                     + phi_g[IDX(i,j,k+1)] + phi_g[IDX(i,j,k-1)]
                                     - 6.0 * phi_g[idx]) / dx2;
                        acc_g[idx] = lapl - alpha_g * rho_src[idx];
                    }
                }
            }
            /* Boundary */
            #pragma omp parallel for schedule(static)
            for (long idx = 0; idx < Ngrid; idx++) {
                int i = idx / (N * N);
                int j = (idx / N) % N;
                int k = idx % N;
                if (i == 0 || i == N-1 || j == 0 || j == N-1 || k == 0 || k == N-1)
                    acc_g[idx] = 0.0;
            }
        }

        /* φ_a field accelerations */
        for (int a = 0; a < 3; a++) {
            #pragma omp parallel for schedule(static)
            for (int i = 1; i < N-1; i++) {
                for (int j = 1; j < N-1; j++) {
                    for (int k = 1; k < N-1; k++) {
                        long idx = IDX(i,j,k);

                        double lapl = (phi[a][IDX(i+1,j,k)] + phi[a][IDX(i-1,j,k)]
                                     + phi[a][IDX(i,j+1,k)] + phi[a][IDX(i,j-1,k)]
                                     + phi[a][IDX(i,j,k+1)] + phi[a][IDX(i,j,k-1)]
                                     - 6.0 * phi[a][idx]) / dx2;

                        /* Effective mass: m_eff² = m₀² - β·Φ */
                        double m_eff2 = m02 - beta_g * phi_g[idx];

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

                        acc[a][idx] = lapl - m_eff2 * phi[a][idx] - dVdphi;
                    }
                }
            }

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

    /* Output setup */
    char tspath[600];
    snprintf(tspath, sizeof(tspath), "%s/grav_ts.tsv", outdir);
    FILE *fts = fopen(tspath, "w");
    fprintf(fts, "time\tphi1_0\tphi_g_0\tpeak1\tpeak_g\t"
                 "E_kin\tE_grad\tE_mass\tE_pot\tE_total\tf_core\t"
                 "Phi_center\tPhi_r10\tQ_monopole\n");

    /* DFT storage */
    int max_dft = 100000;
    double *phi0_hist = malloc(max_dft * sizeof(double));
    double *t_hist = malloc(max_dft * sizeof(double));
    int n_dft = 0;

    int diag_every = Nt / 5000;
    if (diag_every < 1) diag_every = 1;
    int print_every = Nt / 100;
    if (print_every < 1) print_every = 1;
    int dft_every = Nt / max_dft;
    if (dft_every < 1) dft_every = 1;

    double core_radius = 3.0 * sigma;
    long ic = IDX(N/2, N/2, N/2);
    double wall_start = omp_get_wtime();

    /* Profile output times */
    double profile_times[] = {0, 50, 100, 200, 500};
    int n_profile_times = sizeof(profile_times)/sizeof(profile_times[0]);
    int next_profile = 0;

    for (int n = 0; n <= Nt; n++) {
        double t = n * dt;

        if (n % dft_every == 0 && n_dft < max_dft) {
            phi0_hist[n_dft] = phi[0][ic];
            t_hist[n_dft] = t;
            n_dft++;
        }

        int do_diag  = (n % diag_every == 0);
        int do_print = (n % print_every == 0);

        /* Radial profile snapshot (including Φ) */
        if (next_profile < n_profile_times && t >= profile_times[next_profile] - dt/2) {
            char profpath[600];
            snprintf(profpath, sizeof(profpath), "%s/grav_profile_t%04d.tsv",
                     outdir, (int)profile_times[next_profile]);
            FILE *fp = fopen(profpath, "w");
            if (fp) {
                fprintf(fp, "r\tphi1\tPhi_g\trho_E\n");
                int nr = (int)(L / dx);
                double *sum_phi1 = calloc(nr, sizeof(double));
                double *sum_phig = calloc(nr, sizeof(double));
                double *sum_rho  = calloc(nr, sizeof(double));
                int *cnt = calloc(nr, sizeof(int));

                #pragma omp parallel
                {
                    double *ls1 = calloc(nr, sizeof(double));
                    double *lsg = calloc(nr, sizeof(double));
                    double *lsr = calloc(nr, sizeof(double));
                    int *lc = calloc(nr, sizeof(int));

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
                                long idx2 = IDX(i,j,k);

                                ls1[ir] += phi[0][idx2];
                                lsg[ir] += phi_g[idx2];

                                /* Full energy density */
                                double e = 0;
                                for (int a = 0; a < 3; a++) {
                                    e += 0.5 * vel[a][idx2] * vel[a][idx2];
                                    double gx = (phi[a][IDX(i+1,j,k)] - phi[a][IDX(i-1,j,k)]) / (2*dx);
                                    double gy = (phi[a][IDX(i,j+1,k)] - phi[a][IDX(i,j-1,k)]) / (2*dx);
                                    double gz = (phi[a][IDX(i,j,k+1)] - phi[a][IDX(i,j,k-1)]) / (2*dx);
                                    e += 0.5 * (gx*gx + gy*gy + gz*gz);
                                    double m2_loc = m02 - beta_g * phi_g[idx2];
                                    e += 0.5 * m2_loc * phi[a][idx2] * phi[a][idx2];
                                }
                                double P = phi[0][idx2] * phi[1][idx2] * phi[2][idx2];
                                double P2 = P * P;
                                e += 0.5 * mu * P2 / (1.0 + kappa * P2);
                                lsr[ir] += e;
                                lc[ir]++;
                            }
                        }
                    }

                    #pragma omp critical
                    {
                        for (int ir = 0; ir < nr; ir++) {
                            sum_phi1[ir] += ls1[ir];
                            sum_phig[ir] += lsg[ir];
                            sum_rho[ir]  += lsr[ir];
                            cnt[ir] += lc[ir];
                        }
                    }
                    free(ls1); free(lsg); free(lsr); free(lc);
                }

                for (int ir = 0; ir < nr; ir++) {
                    if (cnt[ir] == 0) continue;
                    double r = (ir + 0.5) * dx;
                    fprintf(fp, "%.6f\t%.6e\t%.6e\t%.6e\n",
                            r, sum_phi1[ir]/cnt[ir], sum_phig[ir]/cnt[ir],
                            sum_rho[ir]/cnt[ir]);
                }
                fclose(fp);
                free(sum_phi1); free(sum_phig); free(sum_rho); free(cnt);
                printf("  [profile: t=%.0f -> %s]\n", t, profpath);
            }
            next_profile++;
        }

        /* Diagnostics */
        if (do_diag || do_print) {
            double Ek=0, Eg=0, Em=0, Ep=0, Ecore=0, Eall=0;
            double peak1=0, peak_g=0;
            double Q_mono = 0;  /* monopole moment ∫ρ d³x */

            /* Φ at center and at r≈10 */
            double Phi_c = phi_g[ic];
            /* Find Φ at r≈10 (sample a point on x-axis) */
            int ir10 = (int)(10.0 / dx) + N/2;
            double Phi_r10 = 0;
            if (ir10 >= 0 && ir10 < N)
                Phi_r10 = phi_g[IDX(ir10, N/2, N/2)];

            #pragma omp parallel
            {
                double lEk=0, lEg=0, lEm=0, lEp=0, lEc=0, lEa=0, lQ=0;
                double lpk1=0, lpkg=0;

                #pragma omp for schedule(static) nowait
                for (int i = 1; i < N-1; i++) {
                    double x = -L + i * dx;
                    for (int j = 1; j < N-1; j++) {
                        double y = -L + j * dx;
                        for (int k = 1; k < N-1; k++) {
                            double z = -L + k * dx;
                            long idx2 = IDX(i,j,k);
                            double r = sqrt(x*x + y*y + z*z);
                            double dvol = dV;

                            double e_loc = 0;
                            for (int a = 0; a < 3; a++) {
                                double v2 = vel[a][idx2] * vel[a][idx2];
                                lEk += 0.5 * v2 * dvol;
                                e_loc += 0.5 * v2;

                                double gx = (phi[a][IDX(i+1,j,k)] - phi[a][IDX(i-1,j,k)]) / (2*dx);
                                double gy = (phi[a][IDX(i,j+1,k)] - phi[a][IDX(i,j-1,k)]) / (2*dx);
                                double gz = (phi[a][IDX(i,j,k+1)] - phi[a][IDX(i,j,k-1)]) / (2*dx);
                                lEg += 0.5 * (gx*gx+gy*gy+gz*gz) * dvol;
                                e_loc += 0.5 * (gx*gx+gy*gy+gz*gz);

                                double m2_loc = m02 - beta_g * phi_g[idx2];
                                lEm += 0.5 * m2_loc * phi[a][idx2] * phi[a][idx2] * dvol;
                                e_loc += 0.5 * m2_loc * phi[a][idx2] * phi[a][idx2];

                                double ap = fabs(phi[a][idx2]);
                                if (a == 0 && ap > lpk1) lpk1 = ap;
                            }

                            double P = phi[0][idx2]*phi[1][idx2]*phi[2][idx2];
                            double P2 = P*P;
                            double Vloc = 0.5 * mu * P2 / (1.0 + kappa * P2);
                            lEp += Vloc * dvol;
                            e_loc += Vloc;

                            lEa += e_loc * dvol;
                            lQ += rho_src[idx2] * dvol;
                            if (r < core_radius) lEc += e_loc * dvol;

                            double apg = fabs(phi_g[idx2]);
                            if (apg > lpkg) lpkg = apg;
                        }
                    }
                }

                #pragma omp critical
                {
                    Ek += lEk; Eg += lEg; Em += lEm; Ep += lEp;
                    Ecore += lEc; Eall += lEa; Q_mono += lQ;
                    if (lpk1 > peak1) peak1 = lpk1;
                    if (lpkg > peak_g) peak_g = lpkg;
                }
            }

            double Et = Ek + Eg + Em + Ep;
            double fc = (Eall > 1e-20) ? Ecore / Eall : 0.0;

            if (do_diag)
                fprintf(fts, "%.6f\t%.6e\t%.6e\t%.6e\t%.6e\t"
                             "%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6f\t"
                             "%.6e\t%.6e\t%.6e\n",
                        t, phi[0][ic], phi_g[ic], peak1, peak_g,
                        Ek, Eg, Em, Ep, Et, fc,
                        Phi_c, Phi_r10, Q_mono);

            if (do_print) {
                double wall = omp_get_wtime() - wall_start;
                double frac = (double)n / Nt;
                double eta = (frac > 0.001) ? wall * (1.0 - frac) / frac : 0;
                printf("  t=%7.1f  pk1=%.3f  E=%+.2f  fc=%.4f  "
                       "Phi_c=%.4e  Phi_r10=%.4e  Q=%.2f  "
                       "[%.0fs, ETA %.0fs]\n",
                       t, peak1, Et, fc, Phi_c, Phi_r10, Q_mono,
                       wall, eta);
                fflush(stdout);
            }
        }

        if (n == Nt) break;

        /* Velocity Verlet */
        /* 1. Half-kick (all fields including Φ) */
        for (int a = 0; a < 3; a++) {
            #pragma omp parallel for schedule(static)
            for (long idx2 = 0; idx2 < Ngrid; idx2++)
                vel[a][idx2] += 0.5 * dt * acc[a][idx2];
        }
        if (alpha_g > 0) {
            #pragma omp parallel for schedule(static)
            for (long idx2 = 0; idx2 < Ngrid; idx2++)
                vel_g[idx2] += 0.5 * dt * acc_g[idx2];
        }

        /* 2. Drift */
        for (int a = 0; a < 3; a++) {
            #pragma omp parallel for schedule(static)
            for (long idx2 = 0; idx2 < Ngrid; idx2++)
                phi[a][idx2] += dt * vel[a][idx2];
        }
        if (alpha_g > 0) {
            #pragma omp parallel for schedule(static)
            for (long idx2 = 0; idx2 < Ngrid; idx2++)
                phi_g[idx2] += dt * vel_g[idx2];
        }

        /* 3. Recompute */
        compute_acc();

        /* 4. Half-kick */
        for (int a = 0; a < 3; a++) {
            #pragma omp parallel for schedule(static)
            for (long idx2 = 0; idx2 < Ngrid; idx2++)
                vel[a][idx2] += 0.5 * dt * acc[a][idx2];
        }
        if (alpha_g > 0) {
            #pragma omp parallel for schedule(static)
            for (long idx2 = 0; idx2 < Ngrid; idx2++)
                vel_g[idx2] += 0.5 * dt * acc_g[idx2];
        }

        /* 5. Absorbing boundary */
        for (int a = 0; a < 3; a++) {
            #pragma omp parallel for schedule(static)
            for (long idx2 = 0; idx2 < Ngrid; idx2++) {
                vel[a][idx2] *= damp[idx2];
                phi[a][idx2] *= damp[idx2];
            }
        }
        if (alpha_g > 0) {
            #pragma omp parallel for schedule(static)
            for (long idx2 = 0; idx2 < Ngrid; idx2++) {
                vel_g[idx2] *= damp[idx2];
                phi_g[idx2] *= damp[idx2];
            }
        }
    }

    fclose(fts);

    /* DFT */
    int dft_start = n_dft / 2;
    if (n_dft - dft_start > 100) {
        char dftpath[600];
        snprintf(dftpath, sizeof(dftpath), "%s/grav_spectrum.tsv", outdir);
        FILE *fdft = fopen(dftpath, "w");
        fprintf(fdft, "omega\tpower\n");
        double T = t_hist[n_dft-1] - t_hist[dft_start];
        int nf = 500;
        double peak_pow = 0, peak_om = 0;
        for (int kk = 0; kk < nf; kk++) {
            double omega = 3.0 * mass0 * kk / nf;
            double re = 0, im = 0;
            for (int jj = dft_start; jj < n_dft; jj++) {
                double dtj = (jj > dft_start) ?
                    (t_hist[jj]-t_hist[jj-1]) : (t_hist[dft_start+1]-t_hist[dft_start]);
                re += phi0_hist[jj] * cos(omega * t_hist[jj]) * dtj;
                im += phi0_hist[jj] * sin(omega * t_hist[jj]) * dtj;
            }
            double pw = (re*re + im*im) / (T*T);
            fprintf(fdft, "%.6f\t%.6e\n", omega, pw);
            if (pw > peak_pow) { peak_pow = pw; peak_om = omega; }
        }
        fclose(fdft);
        printf("\nSpectrum: peak omega = %.4f (mass gap = %.4f)\n", peak_om, mass0);
    }

    double wall_total = omp_get_wtime() - wall_start;
    printf("\n=== SUMMARY ===\n");
    printf("  Wall time: %.1f s (%.1f min)\n", wall_total, wall_total/60);
    printf("  Output: %s\n", tspath);
    printf("===\n");

    for (int a = 0; a < 3; a++) { free(phi[a]); free(vel[a]); free(acc[a]); }
    free(phi_g); free(vel_g); free(acc_g); free(rho_src); free(damp);
    free(phi0_hist); free(t_hist);
    return 0;
}
