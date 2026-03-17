/*
 * dynB.c — V26-DynB: Rotating Braid (Spin from Angular Rotation)
 *
 * Initialization: phi_a = A(r_perp) * cos(theta + kz + 2*pi*a/3)
 *                 vel_a = +Omega * A(r_perp) * sin(theta + kz + 2*pi*a/3)
 *
 * where theta = atan2(y,x) is the azimuthal angle.
 *
 * Scans Omega = {0.05, 0.1, 0.2, 0.5}.
 * Lagrangian: L = sum_a [1/2(dt phi_a)^2 - 1/2|grad phi_a|^2 - 1/2 m^2 phi_a^2]
 *             - (mu/2) P^2/(1 + kappa P^2)
 *
 * Compile: gcc -O3 -fopenmp -Wall -o dynB src/dynB.c -lm
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

/* ─── Parameters ─── */
static double mu_pot    = -20.0;
static double kappa     = 20.0;
static double mass      = 1.0;
static double A0        = 0.8;
static double Omega     = 0.1;   /* rotation frequency */

/* Braid parameters */
static double R_tube    = 3.0;   /* tube radius */

static int    N         = 128;
static double L         = 20.0;
static double tfinal    = 500.0;
static double cfl_frac  = 0.20;
static char   outdir[512] = "data";

/* Derived */
static double dx, dx2, m2, dt;
static long   Ngrid;

/* ─── Index helpers ─── */
#define IDX(i,j,k) ((long)(i)*N*N + (long)(j)*N + (long)(k))

/* ─── Globals ─── */
static double *phi[3], *vel[3], *acc[3];
static double *damp;

/* ─── Periodic BC helper: wrap index in z direction ─── */
static inline int wrap_z(int k)
{
    if (k < 0)   return k + N;
    if (k >= N)  return k - N;
    return k;
}

/* ─── Initialization: Rotating Braid ─── */
static void init_rotating_braid(double Om)
{
    double k_twist = 2.0 * M_PI / (2.0 * L);  /* one full twist across domain */

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
        double theta = atan2(y, x);

        for (int a = 0; a < 3; a++) {
            double phase = theta + k_twist * z + 2.0 * M_PI * a / 3.0;
            phi[a][idx] = envelope * cos(phase);
            vel[a][idx] = Om * envelope * sin(phase);  /* rotation velocity */
        }
    }
}

/* ─── Compute acceleration (triple-product potential only) ─── */
static void compute_acc(void)
{
    for (int a = 0; a < 3; a++) {
        #pragma omp parallel for schedule(static)
        for (int i = 2; i < N-2; i++) {
            for (int j = 2; j < N-2; j++) {
                for (int k = 0; k < N; k++) {
                    long idx = IDX(i, j, k);

                    /* 7-point Laplacian with periodic z */
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
                    double dVdphi = mu_pot * P * dP / denom2;

                    acc[a][idx] = lapl - m2 * phi[a][idx] - dVdphi;
                }
            }
        }

        /* Boundary: zero acceleration in x,y (periodic in z) */
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

/* ─── Velocity Verlet step ─── */
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

    /* Absorbing boundary damping (x,y only; z is periodic) */
    for (int a = 0; a < 3; a++) {
        #pragma omp parallel for schedule(static)
        for (long idx = 0; idx < Ngrid; idx++) {
            vel[a][idx] *= damp[idx];
            phi[a][idx] *= damp[idx];
        }
    }
}

/* ─── Diagnostics ─── */
typedef struct {
    double Ek, Eg, Em, Ep, Et;
    double peak[3];
    double fc;
    double peak_P;
    double Lz;      /* z-component of angular momentum */
} Diag;

static Diag compute_diag(double core_radius)
{
    Diag d;
    memset(&d, 0, sizeof(d));

    double Ek=0, Eg=0, Em=0, Ep=0;
    double Ecore=0, Eall=0;
    double peak[3] = {0,0,0};
    double peak_P = 0;
    double Lz_total = 0;

    #pragma omp parallel
    {
        double lEk=0, lEg=0, lEm=0, lEp=0, lEc=0, lEa=0;
        double lpk[3] = {0,0,0};
        double lpkP = 0;
        double lLz = 0;

        #pragma omp for schedule(static) nowait
        for (int i = 2; i < N-2; i++) {
            double x = -L + i * dx;
            for (int j = 2; j < N-2; j++) {
                double y = -L + j * dx;
                for (int k = 0; k < N; k++) {
                    double z = -L + k * dx;
                    (void)z;
                    long idx = IDX(i,j,k);
                    double dV = dx * dx * dx;
                    double e_loc = 0;

                    for (int a = 0; a < 3; a++) {
                        double v2 = vel[a][idx] * vel[a][idx];
                        lEk += 0.5 * v2 * dV;
                        e_loc += 0.5 * v2;

                        int km1 = wrap_z(k-1);
                        int kp1 = wrap_z(k+1);
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

                        /* Angular momentum: L_z = sum_a integral (x * dphi_a/dy - y * dphi_a/dx) * vel_a */
                        /* Actually: L_z = sum_a integral (x * v_a * dphi_a/dy - y * v_a * dphi_a/dx) ... */
                        /* Canonical: L_z = -sum_a integral pi_a (x d_y phi_a - y d_x phi_a) */
                        /* pi_a = d_t phi_a = vel_a */
                        lLz += vel[a][idx] * (x * gy - y * gx) * dV;
                    }

                    double p0 = phi[0][idx], p1 = phi[1][idx], p2 = phi[2][idx];
                    double P = p0 * p1 * p2;
                    double P2 = P * P;
                    double Vloc = 0.5 * mu_pot * P2 / (1.0 + kappa * P2);
                    lEp += Vloc * dV;
                    e_loc += Vloc;

                    double absP = fabs(P);
                    if (absP > lpkP) lpkP = absP;

                    lEa += e_loc * dV;
                    double r = sqrt(x*x + y*y);
                    if (r < core_radius) lEc += e_loc * dV;
                }
            }
        }

        #pragma omp critical
        {
            Ek += lEk; Eg += lEg; Em += lEm; Ep += lEp;
            Ecore += lEc; Eall += lEa;
            Lz_total += lLz;
            for (int a = 0; a < 3; a++)
                if (lpk[a] > peak[a]) peak[a] = lpk[a];
            if (lpkP > peak_P) peak_P = lpkP;
        }
    }

    d.Ek = Ek; d.Eg = Eg; d.Em = Em; d.Ep = Ep;
    d.Et = Ek + Eg + Em + Ep;
    for (int a = 0; a < 3; a++) d.peak[a] = peak[a];
    d.fc = (Eall > 1e-20) ? Ecore / Eall : 0.0;
    d.peak_P = peak_P;
    d.Lz = Lz_total;

    return d;
}

/* ─── Legendre polynomials ─── */
static double legendre_P(int l, double x)
{
    switch (l) {
        case 0: return 1.0;
        case 1: return x;
        case 2: return 0.5*(3.0*x*x - 1.0);
        default: return 0.0;
    }
}

/* ─── Compute strain tensor at a grid point ─── */
static void compute_strain(int i, int j, int k, double eps[6])
{
    double dphi[3][3];
    int km1 = wrap_z(k-1);
    int kp1 = wrap_z(k+1);
    for (int a = 0; a < 3; a++) {
        dphi[a][0] = (phi[a][IDX(i+1,j,k)] - phi[a][IDX(i-1,j,k)]) / (2.0*dx);
        dphi[a][1] = (phi[a][IDX(i,j+1,k)] - phi[a][IDX(i,j-1,k)]) / (2.0*dx);
        dphi[a][2] = (phi[a][IDX(i,j,kp1)] - phi[a][IDX(i,j,km1)]) / (2.0*dx);
    }

    eps[0] = dphi[0][0];                              /* eps_xx */
    eps[1] = dphi[1][1];                              /* eps_yy */
    eps[2] = dphi[2][2];                              /* eps_zz */
    eps[3] = 0.5*(dphi[1][0] + dphi[0][1]);           /* eps_xy */
    eps[4] = 0.5*(dphi[2][0] + dphi[0][2]);           /* eps_xz */
    eps[5] = 0.5*(dphi[2][1] + dphi[1][2]);           /* eps_yz */
}

/* ─── Sample strain on spherical shell, compute multipoles ─── */
static void compute_strain_multipoles(double R_shell, double c_l_out[3])
{
    int N_sample = 200;
    double golden_ratio = (1.0 + sqrt(5.0)) / 2.0;

    double c_l[3] = {0, 0, 0};
    double sum_abs_sigma2 = 0;
    int n_valid = 0;

    for (int ns = 0; ns < N_sample; ns++) {
        double cos_th = 1.0 - 2.0 * (ns + 0.5) / N_sample;
        double sin_th = sqrt(1.0 - cos_th * cos_th);
        double phi_ang = 2.0 * M_PI * ns / golden_ratio;

        double x = R_shell * sin_th * cos(phi_ang);
        double y = R_shell * sin_th * sin(phi_ang);
        double z = R_shell * cos_th;

        int i = (int)((x + L) / dx + 0.5);
        int j_idx = (int)((y + L) / dx + 0.5);
        int kk = (int)((z + L) / dx + 0.5);
        /* wrap z for periodic BC */
        if (kk < 0) kk += N;
        if (kk >= N) kk -= N;

        if (i < 2 || i >= N-2 || j_idx < 2 || j_idx >= N-2)
            continue;

        double eps[6];
        compute_strain(i, j_idx, kk, eps);

        double theta_tr = eps[0] + eps[1] + eps[2];

        double sig[6];
        sig[0] = eps[0] - theta_tr / 3.0;
        sig[1] = eps[1] - theta_tr / 3.0;
        sig[2] = eps[2] - theta_tr / 3.0;
        sig[3] = eps[3];
        sig[4] = eps[4];
        sig[5] = eps[5];

        double sigma2 = sig[0]*sig[0] + sig[1]*sig[1] + sig[2]*sig[2]
                       + 2.0*(sig[3]*sig[3] + sig[4]*sig[4] + sig[5]*sig[5]);

        for (int l = 0; l <= 2; l++)
            c_l[l] += sigma2 * legendre_P(l, cos_th);
        sum_abs_sigma2 += sigma2;
        n_valid++;
    }

    for (int l = 0; l <= 2; l++)
        c_l_out[l] = c_l[l];

    double sum_abs = fabs(c_l[0]) + fabs(c_l[1]) + fabs(c_l[2]);
    printf("    Strain multipoles (R=%.1f, %d pts):\n", R_shell, n_valid);
    for (int l = 0; l <= 2; l++) {
        double frac = (sum_abs > 1e-30) ? fabs(c_l[l]) / sum_abs : 0.0;
        printf("      l=%d: coeff=%.6e  frac=%.4f\n", l, c_l[l], frac);
    }
    printf("    Mean |sigma|^2 = %.6e\n",
           (n_valid > 0) ? sum_abs_sigma2 / n_valid : 0.0);
}

/* ─── DFT to find peak frequency ─── */
static double find_peak_omega(double *hist, double *t_hist, int n_pts, int start,
                              double omega_max)
{
    double T = t_hist[n_pts-1] - t_hist[start];
    if (T < 10.0 || n_pts - start < 50) return 0.0;

    double mean = 0;
    for (int j = start; j < n_pts; j++) mean += hist[j];
    mean /= (n_pts - start);

    int nf = 500;
    double peak_pow = 0.0, peak_om = 0.0;
    for (int kk = 1; kk < nf; kk++) {
        double omega = omega_max * kk / nf;
        double re = 0, im = 0;
        for (int j = start; j < n_pts; j++) {
            double dtj = (j > start) ?
                (t_hist[j]-t_hist[j-1]) : (t_hist[start+1]-t_hist[start]);
            double val = hist[j] - mean;
            re += val * cos(omega * t_hist[j]) * dtj;
            im += val * sin(omega * t_hist[j]) * dtj;
        }
        double pw = re*re + im*im;
        if (pw > peak_pow) { peak_pow = pw; peak_om = omega; }
    }
    return peak_om;
}

/* ─── Run one Omega value ─── */
typedef struct {
    double Omega;
    double fc_init, fc_final;
    double peak_P_init, peak_P_final;
    double Lz_init, Lz_final;
    double Et_init, Et_final;
    double Ek_init, Ek_final;
    double cl_init[3], cl_final[3];
    double omega_peak;  /* DFT peak of center rho */
    int    survived;
} RunResult;

static RunResult run_omega(double Om)
{
    RunResult res;
    memset(&res, 0, sizeof(res));
    res.Omega = Om;

    printf("\n================================================================\n");
    printf("  Omega = %.3f\n", Om);
    printf("================================================================\n\n");

    /* Zero fields */
    for (int a = 0; a < 3; a++) {
        memset(phi[a], 0, Ngrid * sizeof(double));
        memset(vel[a], 0, Ngrid * sizeof(double));
        memset(acc[a], 0, Ngrid * sizeof(double));
    }

    /* Initialize rotating braid */
    init_rotating_braid(Om);

    /* Set damping (x,y only, z periodic) */
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

    double core_radius = 8.0;
    int Nt = (int)(tfinal / dt) + 1;

    /* DFT storage */
    int max_dft = 20000;
    double *rho_hist = malloc(max_dft * sizeof(double));
    double *Lz_hist  = malloc(max_dft * sizeof(double));
    double *t_hist   = malloc(max_dft * sizeof(double));
    int n_dft = 0;

    int dft_every = Nt / max_dft;
    if (dft_every < 1) dft_every = 1;
    int print_every = Nt / 10;
    if (print_every < 1) print_every = 1;
    int diag_every = Nt / 200;
    if (diag_every < 1) diag_every = 1;

    /* Output file */
    char path[600];
    snprintf(path, sizeof(path), "%s/dynB_Om%.3f.tsv", outdir, Om);
    FILE *fout = fopen(path, "w");
    if (!fout) { fprintf(stderr, "Cannot open %s\n", path); res.survived = 0; return res; }
    fprintf(fout, "time\tE_total\tE_kin\tE_grad\tE_mass\tE_pot\t"
                  "fc\tpeak0\tpeak1\tpeak2\tpeak_P\tLz\n");

    compute_acc();

    /* Initial diagnostics */
    Diag d0 = compute_diag(core_radius);
    res.fc_init = d0.fc;
    res.peak_P_init = d0.peak_P;
    res.Lz_init = d0.Lz;
    res.Et_init = d0.Et;
    res.Ek_init = d0.Ek;

    printf("  Initial: E=%.2f  Ek=%.2f  fc=%.4f  |P|=%.6f  Lz=%.4f\n",
           d0.Et, d0.Ek, d0.fc, d0.peak_P, d0.Lz);

    /* Initial strain multipoles */
    printf("  Initial multipoles:\n");
    compute_strain_multipoles(8.0, res.cl_init);

    double wall_start = omp_get_wtime();

    for (int n = 0; n <= Nt; n++) {
        double t = n * dt;

        /* DFT history */
        if (n % dft_every == 0 && n_dft < max_dft) {
            long ic = IDX(N/2, N/2, N/2);
            /* Quick center rho */
            double rho = 0;
            int ci = N/2, cj = N/2, ck = N/2;
            for (int a = 0; a < 3; a++) {
                rho += 0.5 * vel[a][ic] * vel[a][ic];
                double gx = (phi[a][IDX(ci+1,cj,ck)] - phi[a][IDX(ci-1,cj,ck)]) / (2*dx);
                double gy = (phi[a][IDX(ci,cj+1,ck)] - phi[a][IDX(ci,cj-1,ck)]) / (2*dx);
                int km1 = wrap_z(ck-1), kp1 = wrap_z(ck+1);
                double gz = (phi[a][IDX(ci,cj,kp1)] - phi[a][IDX(ci,cj,km1)]) / (2*dx);
                rho += 0.5 * (gx*gx + gy*gy + gz*gz);
                rho += 0.5 * m2 * phi[a][ic] * phi[a][ic];
            }
            double p0 = phi[0][ic], p1 = phi[1][ic], p2 = phi[2][ic];
            double P = p0 * p1 * p2;
            double P2 = P * P;
            rho += 0.5 * mu_pot * P2 / (1.0 + kappa * P2);

            rho_hist[n_dft] = rho;
            t_hist[n_dft] = t;
            n_dft++;
        }

        int do_diag  = (n % diag_every == 0);
        int do_print = (n % print_every == 0);

        if (do_diag || do_print) {
            Diag d = compute_diag(core_radius);

            if (do_diag)
                fprintf(fout, "%.4f\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t"
                              "%.6f\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\n",
                        t, d.Et, d.Ek, d.Eg, d.Em, d.Ep,
                        d.fc, d.peak[0], d.peak[1], d.peak[2], d.peak_P, d.Lz);

            if (do_print) {
                double elapsed = omp_get_wtime() - wall_start;
                double frac = (double)n / Nt;
                double eta_t = (frac > 0.001) ? elapsed * (1.0-frac)/frac : 0;
                printf("  t=%7.1f  E=%.2f  fc=%.4f  |P|=%.6f  Lz=%.4f  [%.0fs, ETA %.0fs]\n",
                       t, d.Et, d.fc, d.peak_P, d.Lz, elapsed, eta_t);
                fflush(stdout);
            }
        }

        if (n == Nt) break;
        verlet_step();
    }

    fclose(fout);

    /* Final diagnostics */
    Diag dfinal = compute_diag(core_radius);
    res.fc_final = dfinal.fc;
    res.peak_P_final = dfinal.peak_P;
    res.Lz_final = dfinal.Lz;
    res.Et_final = dfinal.Et;
    res.Ek_final = dfinal.Ek;
    res.survived = (dfinal.fc > 0.01 && dfinal.peak_P > 1e-6) ? 1 : 0;

    printf("  Final: E=%.2f  fc=%.4f  |P|=%.6f  Lz=%.4f  survived=%s\n",
           dfinal.Et, dfinal.fc, dfinal.peak_P, dfinal.Lz,
           res.survived ? "YES" : "NO");

    /* Final strain multipoles */
    printf("  Final multipoles:\n");
    compute_strain_multipoles(8.0, res.cl_final);

    /* DFT peak */
    if (n_dft > 100) {
        int start = n_dft / 2;
        res.omega_peak = find_peak_omega(rho_hist, t_hist, n_dft, start, 5.0);
        printf("  DFT peak omega (center rho, 2nd half): %.4f\n", res.omega_peak);
    }

    double elapsed = omp_get_wtime() - wall_start;
    printf("  Wall time: %.1f sec\n", elapsed);

    free(rho_hist);
    free(Lz_hist);
    free(t_hist);

    return res;
}

/* ─── Main ─── */
int main(int argc, char **argv)
{
    /* Parse optional overrides */
    for (int i = 1; i < argc - 1; i += 2) {
        if      (!strcmp(argv[i], "-mu"))     mu_pot   = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-kappa"))  kappa    = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-mass"))   mass     = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-A"))      A0       = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-N"))      N        = atoi(argv[i+1]);
        else if (!strcmp(argv[i], "-L"))      L        = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-tfinal")) tfinal   = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-cfl"))    cfl_frac = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-Rtube"))  R_tube   = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-o"))      strncpy(outdir, argv[i+1], sizeof(outdir)-1);
        else { fprintf(stderr, "Unknown flag: %s\n", argv[i]); exit(1); }
    }

    dx  = 2.0 * L / (N - 1);
    dx2 = dx * dx;
    m2  = mass * mass;
    dt  = cfl_frac * dx;
    Ngrid = (long)N * N * N;

    printf("=== V26-DynB: Rotating Braid ===\n");
    printf("Parameters:\n");
    printf("  mu=%.1f  kappa=%.1f  mass=%.3f  A0=%.3f  R_tube=%.1f\n",
           mu_pot, kappa, mass, A0, R_tube);
    printf("  N=%d  L=%.1f  dx=%.4f  dt=%.5f\n", N, L, dx, dt);
    printf("  Ngrid=%ld (%.1f M)  Memory: %.1f MB\n",
           Ngrid, Ngrid/1e6, Ngrid*8.0*10/1e6);
    printf("  Threads: %d\n", omp_get_max_threads());
    printf("  tfinal=%.0f\n", tfinal);
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

    /* Omega scan */
    double Omega_vals[] = {0.05, 0.1, 0.2, 0.5};
    int n_omega = 4;
    RunResult results[4];

    for (int io = 0; io < n_omega; io++) {
        results[io] = run_omega(Omega_vals[io]);
    }

    /* ─── Summary ─── */
    printf("\n\n");
    printf("================================================================\n");
    printf("  SUMMARY: V26-DynB Rotating Braid\n");
    printf("================================================================\n\n");
    printf("%-8s  %-8s  %-8s  %-10s  %-10s  %-10s  %-10s  %-8s  %-6s\n",
           "Omega", "fc_i", "fc_f", "|P|_i", "|P|_f",
           "Lz_i", "Lz_f", "omega_pk", "surv");
    printf("%-8s  %-8s  %-8s  %-10s  %-10s  %-10s  %-10s  %-8s  %-6s\n",
           "-----", "-----", "-----", "--------", "--------",
           "--------", "--------", "--------", "----");
    for (int io = 0; io < n_omega; io++) {
        RunResult *r = &results[io];
        printf("%-8.3f  %-8.4f  %-8.4f  %-10.6f  %-10.6f  %-10.4f  %-10.4f  %-8.4f  %-6s\n",
               r->Omega, r->fc_init, r->fc_final, r->peak_P_init, r->peak_P_final,
               r->Lz_init, r->Lz_final, r->omega_peak,
               r->survived ? "YES" : "NO");
    }

    printf("\nMultipole fractions (l=2) at final time:\n");
    printf("%-8s  %-12s  %-12s  %-12s  %-12s\n",
           "Omega", "c_l0", "c_l1", "c_l2", "l2_frac");
    for (int io = 0; io < n_omega; io++) {
        RunResult *r = &results[io];
        double sum = fabs(r->cl_final[0]) + fabs(r->cl_final[1]) + fabs(r->cl_final[2]);
        double l2f = (sum > 1e-30) ? fabs(r->cl_final[2]) / sum : 0.0;
        printf("%-8.3f  %-12.4e  %-12.4e  %-12.4e  %-12.4f\n",
               r->Omega, r->cl_final[0], r->cl_final[1], r->cl_final[2], l2f);
    }

    /* Write summary TSV */
    {
        char spath[600];
        snprintf(spath, sizeof(spath), "%s/dynB_summary.tsv", outdir);
        FILE *fs = fopen(spath, "w");
        if (fs) {
            fprintf(fs, "Omega\tfc_init\tfc_final\tpeakP_init\tpeakP_final\t"
                        "Lz_init\tLz_final\tEt_init\tEt_final\tEk_init\tEk_final\t"
                        "cl0_f\tcl1_f\tcl2_f\tl2_frac\tomega_peak\tsurvived\n");
            for (int io = 0; io < n_omega; io++) {
                RunResult *r = &results[io];
                double sum = fabs(r->cl_final[0]) + fabs(r->cl_final[1]) + fabs(r->cl_final[2]);
                double l2f = (sum > 1e-30) ? fabs(r->cl_final[2]) / sum : 0.0;
                fprintf(fs, "%.3f\t%.6f\t%.6f\t%.6e\t%.6e\t%.6f\t%.6f\t"
                            "%.4f\t%.4f\t%.4f\t%.4f\t"
                            "%.6e\t%.6e\t%.6e\t%.4f\t%.4f\t%d\n",
                        r->Omega, r->fc_init, r->fc_final,
                        r->peak_P_init, r->peak_P_final,
                        r->Lz_init, r->Lz_final,
                        r->Et_init, r->Et_final, r->Ek_init, r->Ek_final,
                        r->cl_final[0], r->cl_final[1], r->cl_final[2],
                        l2f, r->omega_peak, r->survived);
            }
            fclose(fs);
            printf("\nSummary written to %s/dynB_summary.tsv\n", outdir);
        }
    }

    /* Cleanup */
    for (int a = 0; a < 3; a++) {
        free(phi[a]); free(vel[a]); free(acc[a]);
    }
    free(damp);

    printf("\n=== V26-DynB Complete ===\n");
    return 0;
}
