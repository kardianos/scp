/*
 * two_oscillon.c — Two-oscillon gravitational interaction in 3D
 *
 * Places two oscillons at (0,0,+D/2) and (0,0,-D/2) and evolves with
 * scalar gravity coupling. Measures center-of-mass separation vs time
 * to detect gravitational attraction.
 *
 * Physics:
 *   Three scalar fields phi_a with triple-product coupling (v21 model)
 *   + massless scalar Phi sourced by full energy density (Option A)
 *   + backreaction m_eff^2 = m0^2 - beta*Phi
 *
 * Key design choices:
 *   - Box: [-L,L]^3 with L=60, N=400 (dx=0.30, same as v21 production)
 *   - Two Gaussians on z-axis separated by D (default 20)
 *   - Gravity initialized from Poisson solution of combined density
 *   - Absorbing boundary at 70-95% of L
 *   - Track centroid of each oscillon (z>0 and z<0 halves)
 *
 * Compile: gcc -O3 -fopenmp -Wall -o two_oscillon v22/src/two_oscillon.c -lm
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

/* Physics parameters */
static double mu       = -20.0;
static double kappa    = 20.0;
static double mass0    = 1.0;
static double A_init   = 0.8;
static double sigma    = 3.0;

/* Gravity */
static double alpha_g  = 0.1;
static double beta_g   = 0.1;

/* Grid */
static int    N        = 400;
static double L        = 60.0;
static double tfinal   = 500.0;
static double cfl_frac = 0.25;

/* Two-body setup */
static double D_sep    = 20.0;   /* initial separation (center-to-center) */
static double phase2   = 0.0;    /* breathing phase offset of lower oscillon (radians) */
static double omega_b  = 0.95;   /* breathing frequency for phase init */

/* Output */
static char   outdir[512] = "v22/data";
static char   tag[128]    = "";

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
        else if (!strcmp(argv[i], "-D"))      D_sep    = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-alpha"))  alpha_g  = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-beta"))   beta_g   = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-phase"))  phase2   = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-omegab")) omega_b = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-tag"))   strncpy(tag, argv[i+1], sizeof(tag)-1);
        else if (!strcmp(argv[i], "-o"))      strncpy(outdir, argv[i+1], sizeof(outdir)-1);
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

    printf("two_oscillon: gravitational interaction test\n");
    printf("  mu=%.3f kappa=%.4f mass0=%.4f A=%.3f sigma=%.3f\n",
           mu, kappa, mass0, A_init, sigma);
    printf("  alpha=%.6e beta=%.6e D=%.1f phase=%.4f\n", alpha_g, beta_g, D_sep, phase2);
    printf("  N=%d L=%.1f dx=%.5f dt=%.6f\n", N, L, dx, dt);
    printf("  Ngrid=%ld (%.1f M, %.1f GB)  Nt=%d  tfinal=%.0f\n",
           Ngrid, Ngrid/1e6, Ngrid*14.0*8/1e9, Nt, tfinal);
    printf("  Threads: %d\n", omp_get_max_threads());
    fflush(stdout);

    /* Allocate arrays */
    double *phi[3], *vel[3], *acc[3];
    double *phi_g, *vel_g, *acc_g;
    double *rho_src, *damp;

    for (int a = 0; a < 3; a++) {
        phi[a] = calloc(Ngrid, sizeof(double));
        vel[a] = calloc(Ngrid, sizeof(double));
        acc[a] = calloc(Ngrid, sizeof(double));
        if (!phi[a] || !vel[a] || !acc[a]) {
            fprintf(stderr, "Out of memory allocating phi/vel/acc[%d]\n", a);
            return 1;
        }
    }
    phi_g   = calloc(Ngrid, sizeof(double));
    vel_g   = calloc(Ngrid, sizeof(double));
    acc_g   = calloc(Ngrid, sizeof(double));
    rho_src = calloc(Ngrid, sizeof(double));
    damp    = malloc(Ngrid * sizeof(double));
    if (!phi_g || !vel_g || !acc_g || !rho_src || !damp) {
        fprintf(stderr, "Out of memory\n");
        return 1;
    }
    printf("  Memory allocated: %.1f GB\n", Ngrid*14.0*8/1e9);
    fflush(stdout);

    /* Absorbing boundary */
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

    /* Initialize: two oscillons on z-axis at z = +D/2 and z = -D/2
     * Upper oscillon at phase 0: phi = A*f(r), v = 0
     * Lower oscillon at phase delta: phi = A*cos(delta)*f(r), v = -A*omega*sin(delta)*f(r) */
    double z1 = +D_sep / 2.0;
    double z2 = -D_sep / 2.0;
    double cos_ph = cos(phase2);
    double sin_ph = sin(phase2);

    #pragma omp parallel for schedule(static)
    for (long idx = 0; idx < Ngrid; idx++) {
        int i = idx / (N * N);
        int j = (idx / N) % N;
        int k = idx % N;
        double x = -L + i * dx;
        double y = -L + j * dx;
        double z = -L + k * dx;

        double r1_2 = x*x + y*y + (z - z1)*(z - z1);
        double r2_2 = x*x + y*y + (z - z2)*(z - z2);
        double g1 = A_init * exp(-r1_2 / (2.0 * sigma * sigma));
        double g2 = A_init * exp(-r2_2 / (2.0 * sigma * sigma));

        for (int a = 0; a < 3; a++) {
            phi[a][idx] = g1 + cos_ph * g2;
            vel[a][idx] = -omega_b * sin_ph * g2;
        }
    }

    /* Initialize Phi from Poisson solution of combined density */
    if (alpha_g > 0) {
        printf("  Computing Poisson initialization...\n");
        fflush(stdout);

        /* Build radial rho from EACH oscillon center, solve individually,
         * then superpose (valid in linear regime) */
        int Nr = (int)(L / dx) + 1;

        /* For the combined system, compute 1D radial profile around origin
         * then use full 3D Poisson via radial average */
        double *rho_rad = calloc(Nr, sizeof(double));
        long *count = calloc(Nr, sizeof(long));

        /* Compute energy density at each point, bin by distance from origin */
        #pragma omp parallel
        {
            double *lrho = calloc(Nr, sizeof(double));
            long *lcnt = calloc(Nr, sizeof(long));

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

        /* Solve Poisson via Green's function (radial) */
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
        fflush(stdout);

        /* Interpolate Phi(r) to 3D grid */
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

    /* --- Compute acceleration --- */
    void compute_acc(void) {
        /* Source: full energy density */
        #pragma omp parallel for schedule(static)
        for (int i = 1; i < N-1; i++) {
            for (int j = 1; j < N-1; j++) {
                for (int k = 1; k < N-1; k++) {
                    long idx = IDX(i,j,k);
                    double e = 0;
                    for (int a = 0; a < 3; a++) {
                        e += 0.5 * vel[a][idx] * vel[a][idx];
                        double gx = (phi[a][IDX(i+1,j,k)] - phi[a][IDX(i-1,j,k)]) / (2*dx);
                        double gy = (phi[a][IDX(i,j+1,k)] - phi[a][IDX(i,j-1,k)]) / (2*dx);
                        double gz = (phi[a][IDX(i,j,k+1)] - phi[a][IDX(i,j,k-1)]) / (2*dx);
                        e += 0.5 * (gx*gx + gy*gy + gz*gz);
                        e += 0.5 * m02 * phi[a][idx] * phi[a][idx];
                    }
                    double P = phi[0][idx] * phi[1][idx] * phi[2][idx];
                    double P2 = P * P;
                    e += 0.5 * mu * P2 / (1.0 + kappa * P2);
                    rho_src[idx] = e;
                }
            }
        }

        /* Phi acceleration: laplacian(Phi) - alpha * rho */
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
            #pragma omp parallel for schedule(static)
            for (long idx = 0; idx < Ngrid; idx++) {
                int ii = idx / (N * N);
                int jj = (idx / N) % N;
                int kk = idx % N;
                if (ii == 0 || ii == N-1 || jj == 0 || jj == N-1 || kk == 0 || kk == N-1)
                    acc_g[idx] = 0.0;
            }
        }

        /* phi_a accelerations */
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

                        double m_eff2 = m02 - beta_g * phi_g[idx];

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
                int ii = idx / (N * N);
                int jj = (idx / N) % N;
                int kk = idx % N;
                if (ii == 0 || ii == N-1 || jj == 0 || jj == N-1 || kk == 0 || kk == N-1)
                    acc[a][idx] = 0.0;
            }
        }
    }

    compute_acc();

    /* --- Output setup --- */
    char tspath[600];
    if (tag[0])
        snprintf(tspath, sizeof(tspath), "%s/two_oscillon_%s_ts.tsv", outdir, tag);
    else
        snprintf(tspath, sizeof(tspath), "%s/two_oscillon_ts.tsv", outdir);
    FILE *fts = fopen(tspath, "w");
    fprintf(fts, "time\tE_total\tE_upper\tE_lower\tfc_upper\tfc_lower\t"
                 "z_upper\tz_lower\tseparation\t"
                 "pk_upper\tpk_lower\tPhi_mid\tQ_total\n");

    /* DFT storage */
    int max_dft = 100000;
    double *phi0_hist = malloc(max_dft * sizeof(double));
    double *t_hist = malloc(max_dft * sizeof(double));
    int n_dft = 0;

    int diag_every = Nt / 10000;
    if (diag_every < 1) diag_every = 1;
    int print_every = Nt / 200;
    if (print_every < 1) print_every = 1;
    int dft_every = Nt / max_dft;
    if (dft_every < 1) dft_every = 1;

    double core_radius = 3.0 * sigma;
    long ic_upper = IDX(N/2, N/2, N/2 + (int)(z1/dx));
    long ic_lower = IDX(N/2, N/2, N/2 + (int)(z2/dx));
    double wall_start = omp_get_wtime();

    /* Profile output at fixed intervals */
    double next_profile_t = 0;
    double profile_interval = 50.0;

    for (int n = 0; n <= Nt; n++) {
        double t = n * dt;

        /* DFT sampling (upper oscillon) */
        if (n % dft_every == 0 && n_dft < max_dft) {
            phi0_hist[n_dft] = phi[0][ic_upper];
            t_hist[n_dft] = t;
            n_dft++;
        }

        int do_diag  = (n % diag_every == 0);
        int do_print = (n % print_every == 0);

        /* Z-axis profile snapshot */
        if (t >= next_profile_t - dt/2) {
            char profpath[600];
            if (tag[0])
                snprintf(profpath, sizeof(profpath), "%s/two_osc_%s_profile_t%04d.tsv",
                         outdir, tag, (int)next_profile_t);
            else
                snprintf(profpath, sizeof(profpath), "%s/two_osc_profile_t%04d.tsv",
                         outdir, (int)next_profile_t);
            FILE *fp = fopen(profpath, "w");
            if (fp) {
                fprintf(fp, "z\tphi1\tPhi_g\trho_E\n");
                int imid = N/2;
                for (int k = 0; k < N; k++) {
                    double z = -L + k * dx;
                    long idx2 = IDX(imid, imid, k);
                    fprintf(fp, "%.6f\t%.6e\t%.6e\t%.6e\n",
                            z, phi[0][idx2], phi_g[idx2], rho_src[idx2]);
                }
                fclose(fp);
                printf("  [profile: t=%.0f -> %s]\n", t, profpath);
            }
            next_profile_t += profile_interval;
        }

        /* Diagnostics: track each oscillon separately */
        if (do_diag || do_print) {
            double E_total = 0, E_upper = 0, E_lower = 0;
            double Ecore_up = 0, Eall_up = 0, Ecore_lo = 0, Eall_lo = 0;
            double pk_upper = 0, pk_lower = 0;
            double Q_total = 0;

            /* Energy-weighted centroid for z > 0 and z < 0 */
            double wz_up = 0, w_up = 0;
            double wz_lo = 0, w_lo = 0;

            #pragma omp parallel
            {
                double lEt=0, lEu=0, lEl=0;
                double lEcu=0, lEau=0, lEcl=0, lEal=0;
                double lpku=0, lpkl=0, lQ=0;
                double lwzu=0, lwu=0, lwzl=0, lwl=0;

                #pragma omp for schedule(static) nowait
                for (int i = 1; i < N-1; i++) {
                    double x = -L + i * dx;
                    for (int j = 1; j < N-1; j++) {
                        double y = -L + j * dx;
                        for (int k = 1; k < N-1; k++) {
                            double z = -L + k * dx;
                            long idx2 = IDX(i,j,k);
                            double dvol = dV;

                            /* Energy density */
                            double e_loc = 0;
                            for (int a = 0; a < 3; a++) {
                                e_loc += 0.5 * vel[a][idx2] * vel[a][idx2];
                                double gx = (phi[a][IDX(i+1,j,k)] - phi[a][IDX(i-1,j,k)]) / (2*dx);
                                double gy = (phi[a][IDX(i,j+1,k)] - phi[a][IDX(i,j-1,k)]) / (2*dx);
                                double gz = (phi[a][IDX(i,j,k+1)] - phi[a][IDX(i,j,k-1)]) / (2*dx);
                                e_loc += 0.5 * (gx*gx + gy*gy + gz*gz);
                                double m2l = m02 - beta_g * phi_g[idx2];
                                e_loc += 0.5 * m2l * phi[a][idx2] * phi[a][idx2];
                            }
                            double P = phi[0][idx2]*phi[1][idx2]*phi[2][idx2];
                            double P2 = P*P;
                            e_loc += 0.5 * mu * P2 / (1.0 + kappa * P2);

                            double edv = e_loc * dvol;
                            lEt += edv;
                            lQ += rho_src[idx2] * dvol;

                            /* Assign to upper or lower oscillon */
                            if (z > 0) {
                                lEu += edv;
                                lEau += edv;
                                double r_up = sqrt(x*x + y*y + (z-z1)*(z-z1));
                                if (r_up < core_radius) lEcu += edv;

                                /* Centroid */
                                lwzu += z * edv;
                                lwu += edv;

                                double ap = fabs(phi[0][idx2]);
                                if (ap > lpku) lpku = ap;
                            } else {
                                lEl += edv;
                                lEal += edv;
                                double r_lo = sqrt(x*x + y*y + (z-z2)*(z-z2));
                                if (r_lo < core_radius) lEcl += edv;

                                lwzl += z * edv;
                                lwl += edv;

                                double ap = fabs(phi[0][idx2]);
                                if (ap > lpkl) lpkl = ap;
                            }
                        }
                    }
                }

                #pragma omp critical
                {
                    E_total += lEt; E_upper += lEu; E_lower += lEl;
                    Ecore_up += lEcu; Eall_up += lEau;
                    Ecore_lo += lEcl; Eall_lo += lEal;
                    Q_total += lQ;
                    wz_up += lwzu; w_up += lwu;
                    wz_lo += lwzl; w_lo += lwl;
                    if (lpku > pk_upper) pk_upper = lpku;
                    if (lpkl > pk_lower) pk_lower = lpkl;
                }
            }

            double fc_up = (Eall_up > 1e-20) ? Ecore_up / Eall_up : 0.0;
            double fc_lo = (Eall_lo > 1e-20) ? Ecore_lo / Eall_lo : 0.0;
            double z_up = (w_up > 1e-20) ? wz_up / w_up : z1;
            double z_lo = (w_lo > 1e-20) ? wz_lo / w_lo : z2;
            double sep = z_up - z_lo;

            /* Phi at midpoint (between the two oscillons) */
            double Phi_mid = phi_g[IDX(N/2, N/2, N/2)];

            if (do_diag)
                fprintf(fts, "%.6f\t%.4f\t%.4f\t%.4f\t%.6f\t%.6f\t"
                             "%.6f\t%.6f\t%.6f\t"
                             "%.4f\t%.4f\t%.6e\t%.4f\n",
                        t, E_total, E_upper, E_lower, fc_up, fc_lo,
                        z_up, z_lo, sep,
                        pk_upper, pk_lower, Phi_mid, Q_total);

            if (do_print) {
                double wall = omp_get_wtime() - wall_start;
                double frac = (double)n / Nt;
                double eta = (frac > 0.001) ? wall * (1.0 - frac) / frac : 0;
                printf("  t=%7.1f  sep=%.4f  z+=%+.3f z-=%+.3f  "
                       "E=%.1f (%.1f+%.1f)  fc=%.3f/%.3f  "
                       "Phi_mid=%.4e  [%.0fs, ETA %.0fs]\n",
                       t, sep, z_up, z_lo,
                       E_total, E_upper, E_lower, fc_up, fc_lo,
                       Phi_mid, wall, eta);
                fflush(stdout);
            }
        }

        if (n == Nt) break;

        /* --- Velocity Verlet --- */
        /* 1. Half-kick */
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

        /* 3. Recompute accelerations */
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

    /* DFT of upper oscillon */
    int dft_start = n_dft / 2;
    if (n_dft - dft_start > 100) {
        char dftpath[600];
        snprintf(dftpath, sizeof(dftpath), "%s/two_osc_spectrum.tsv", outdir);
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
