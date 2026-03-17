/*
 * dynA.c — V26-DynA: Propagating Helical Wave Braid
 *
 * Based on v26.c mode 2 (twisted tube, triple product only).
 * Key change: initial velocity to make the braid PROPAGATE along z.
 *
 * Two runs:
 *   run 0: propagating braid (vel = omega * envelope * sin(phase))
 *   run 1: static braid control (vel = 0)
 *
 * Lagrangian:
 *   L = Sum_a [1/2(dt phi_a)^2 - 1/2|grad phi_a|^2 - 1/2 m^2 phi_a^2]
 *     - (mu/2) P^2 / (1 + kappa P^2)
 *
 * Compile: gcc -O3 -fopenmp -Wall -o dynA src/dynA.c -lm
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
static double R_tube    = 3.0;

static int    N         = 128;
static double L         = 20.0;
static double tfinal    = 500.0;
static double cfl_frac  = 0.20;
static char   outdir[512] = "data";

/* ─── Index helpers ─── */
#define IDX(i,j,k) ((long)(i)*N*N + (long)(j)*N + (long)(k))

/* ─── Globals ─── */
static double *phi[3], *vel[3], *acc[3];
static double *damp;
static double dx, dx2, m2, dt;
static long Ngrid;

/* ─── Periodic z wrap ─── */
static inline int wrap_z(int k)
{
    if (k < 0)   return k + N;
    if (k >= N)  return k - N;
    return k;
}

/* ─── Initialization: Propagating Helical Braid ─── */
static void init_braid(int propagating)
{
    double Lz = 2.0 * L;  /* full domain length in z */
    double k_wave = 2.0 * M_PI / Lz;
    double omega = sqrt(k_wave * k_wave + mass * mass);

    printf("  k = %.6f, omega = %.6f, v_g = k/omega = %.4f\n",
           k_wave, omega, k_wave / omega);

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
            double phase = k_wave * z + 2.0 * M_PI * a / 3.0;
            phi[a][idx] = envelope * cos(phase);
            if (propagating) {
                vel[a][idx] = omega * envelope * sin(phase);
            } else {
                vel[a][idx] = 0.0;
            }
        }
    }
}

/* ─── Compute acceleration (triple product only, periodic z) ─── */
static void compute_acc(void)
{
    for (int a = 0; a < 3; a++) {
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
                    double dVdphi = mu_pot * P * dP / denom2;

                    acc[a][idx] = lapl - m2 * phi[a][idx] - dVdphi;
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

    /* Absorbing boundary damping (x,y only) */
    for (int a = 0; a < 3; a++) {
        #pragma omp parallel for schedule(static)
        for (long idx = 0; idx < Ngrid; idx++) {
            vel[a][idx] *= damp[idx];
            phi[a][idx] *= damp[idx];
        }
    }
}

/* ─── Diagnostics structure ─── */
typedef struct {
    double Ek, Eg, Em, Ep, Et;
    double peak[3];
    double fc;
    double peak_P;
    double rho_center;
    double Pz;          /* z-momentum */
} Diag;

/* ─── Compute diagnostics ─── */
static Diag compute_diag(double core_radius)
{
    Diag d;
    memset(&d, 0, sizeof(d));

    double Ek = 0, Eg = 0, Em = 0, Ep = 0;
    double Ecore = 0, Eall = 0;
    double peak[3] = {0, 0, 0};
    double peak_P = 0;
    double Pz_total = 0;

    #pragma omp parallel
    {
        double lEk = 0, lEg = 0, lEm = 0, lEp = 0, lEc = 0, lEa = 0;
        double lpk[3] = {0, 0, 0};
        double lpkP = 0;
        double lPz = 0;

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

                    double absP = fabs(P_val);
                    if (absP > lpkP) lpkP = absP;

                    lEa += e_loc * dV;
                    double r = sqrt(x*x + y*y + z*z);
                    if (r < core_radius) lEc += e_loc * dV;
                }
            }
        }

        #pragma omp critical
        {
            Ek += lEk; Eg += lEg; Em += lEm; Ep += lEp;
            Ecore += lEc; Eall += lEa;
            Pz_total += lPz;
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
    d.Pz = Pz_total;

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
        d.rho_center = rho;
    }

    return d;
}

/* ─── Fast center density ─── */
static double center_rho(void)
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
    return rho;
}

/* ─── DFT to find peak frequency ─── */
static double find_peak_omega(double *hist, double *t_hist, int n_pts, int start,
                              double omega_max, double *peak_power_out)
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
    if (peak_power_out) *peak_power_out = peak_pow;
    return peak_om;
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
    int km1 = wrap_z(k - 1);
    int kp1 = wrap_z(k + 1);

    double dphi[3][3];
    for (int a = 0; a < 3; a++) {
        dphi[a][0] = (phi[a][IDX(i+1,j,k)] - phi[a][IDX(i-1,j,k)]) / (2.0*dx);
        dphi[a][1] = (phi[a][IDX(i,j+1,k)] - phi[a][IDX(i,j-1,k)]) / (2.0*dx);
        dphi[a][2] = (phi[a][IDX(i,j,kp1)] - phi[a][IDX(i,j,km1)]) / (2.0*dx);
    }

    eps[0] = dphi[0][0];
    eps[1] = dphi[1][1];
    eps[2] = dphi[2][2];
    eps[3] = 0.5*(dphi[1][0] + dphi[0][1]);
    eps[4] = 0.5*(dphi[2][0] + dphi[0][2]);
    eps[5] = 0.5*(dphi[2][1] + dphi[1][2]);
}

/* ─── Sample strain on spherical shell, compute multipoles ─── */
static void compute_strain_multipoles(double R_shell, const char *label)
{
    int N_sample = 200;
    double golden_ratio = (1.0 + sqrt(5.0)) / 2.0;

    double c_l[3] = {0, 0, 0};
    double norm = 0;
    int n_valid = 0;

    char path[600];
    snprintf(path, sizeof(path), "%s/dynA_%s_strain.tsv", outdir, label);
    FILE *fstrain = fopen(path, "w");
    snprintf(path, sizeof(path), "%s/dynA_%s_multipoles.tsv", outdir, label);
    FILE *fmulti = fopen(path, "w");

    if (!fstrain || !fmulti) {
        fprintf(stderr, "Cannot open strain output files\n");
        if (fstrain) fclose(fstrain);
        if (fmulti) fclose(fmulti);
        return;
    }

    fprintf(fstrain, "# Strain on shell r=%.1f (%s)\n", R_shell, label);
    fprintf(fstrain, "theta\tphi_ang\teps_xx\teps_yy\teps_zz\teps_xy\teps_xz\teps_yz\t"
                     "theta_tr\tsigma2\n");

    for (int ns = 0; ns < N_sample; ns++) {
        double cos_th = 1.0 - 2.0 * (ns + 0.5) / N_sample;
        double sin_th = sqrt(1.0 - cos_th * cos_th);
        double phi_ang = 2.0 * M_PI * ns / golden_ratio;

        double x = R_shell * sin_th * cos(phi_ang);
        double y = R_shell * sin_th * sin(phi_ang);
        double z = R_shell * cos_th;

        int i = (int)((x + L) / dx + 0.5);
        int j_idx = (int)((y + L) / dx + 0.5);
        int k = (int)((z + L) / dx + 0.5);

        if (i < 2 || i >= N-2 || j_idx < 2 || j_idx >= N-2 || k < 0 || k >= N)
            continue;

        double eps[6];
        compute_strain(i, j_idx, k, eps);

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

        fprintf(fstrain, "%.6f\t%.6f\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\n",
                acos(cos_th), phi_ang, eps[0], eps[1], eps[2], eps[3], eps[4], eps[5],
                theta_tr, sigma2);

        for (int l = 0; l <= 2; l++)
            c_l[l] += sigma2 * legendre_P(l, cos_th);
        norm += sigma2;
        n_valid++;
    }

    fprintf(fmulti, "# Multipole decomposition of |sigma|^2 on shell r=%.1f (%s)\n",
            R_shell, label);
    fprintf(fmulti, "l\tcoeff\tfraction\n");
    double sum_abs = fabs(c_l[0]) + fabs(c_l[1]) + fabs(c_l[2]);
    for (int l = 0; l <= 2; l++) {
        double frac = (sum_abs > 1e-30) ? fabs(c_l[l]) / sum_abs : 0.0;
        fprintf(fmulti, "%d\t%.6e\t%.6f\n", l, c_l[l], frac);
        printf("    l=%d: coefficient = %.6e, fraction = %.4f\n", l, c_l[l], frac);
    }
    printf("    Total |sigma|^2 on shell: %.6e (from %d valid points)\n",
           (n_valid > 0) ? norm / n_valid : 0.0, n_valid);

    fclose(fstrain);
    fclose(fmulti);
}

/* ─── Run one braid (propagating or static) ─── */
static void run_braid(int propagating)
{
    const char *label = propagating ? "propagating" : "static";
    printf("\n");
    printf("================================================================\n");
    printf("  DynA: %s braid\n", label);
    printf("================================================================\n\n");

    m2 = mass * mass;

    /* Zero fields */
    for (int a = 0; a < 3; a++) {
        memset(phi[a], 0, Ngrid * sizeof(double));
        memset(vel[a], 0, Ngrid * sizeof(double));
        memset(acc[a], 0, Ngrid * sizeof(double));
    }

    /* Initialize */
    init_braid(propagating);

    /* Set damping: absorbing in x,y only (periodic in z) */
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
    int max_dft = 50000;
    double *rho_hist = malloc(max_dft * sizeof(double));
    double *phi0_hist = malloc(max_dft * sizeof(double));
    double *t_hist = malloc(max_dft * sizeof(double));
    int n_dft = 0;

    int diag_every = Nt / 5000;
    if (diag_every < 1) diag_every = 1;
    int print_every = Nt / 30;
    if (print_every < 1) print_every = 1;
    int dft_every = Nt / max_dft;
    if (dft_every < 1) dft_every = 1;

    /* Output file */
    char path[600];
    snprintf(path, sizeof(path), "%s/dynA_%s_phase1.tsv", outdir, label);
    FILE *f1 = fopen(path, "w");
    if (!f1) { fprintf(stderr, "Cannot open %s\n", path); return; }
    fprintf(f1, "time\tE_total\tE_kin\tE_grad\tE_mass\tE_pot\t"
                "fc\tpeak0\tpeak1\tpeak2\tpeak_P\trho_center\tPz\n");

    compute_acc();

    double wall_start = omp_get_wtime();

    /* Initial diagnostics */
    {
        Diag d0 = compute_diag(core_radius);
        printf("  t=%7.1f  E=%.2f  fc=%.4f  pk=(%.4f,%.4f,%.4f)  |P|=%.6f  Pz=%.4f\n",
               0.0, d0.Et, d0.fc, d0.peak[0], d0.peak[1], d0.peak[2], d0.peak_P, d0.Pz);
    }

    for (int n = 0; n <= Nt; n++) {
        double t = n * dt;

        /* DFT history */
        if (n % dft_every == 0 && n_dft < max_dft) {
            long ic = IDX(N/2, N/2, N/2);
            rho_hist[n_dft] = center_rho();
            phi0_hist[n_dft] = phi[0][ic];
            t_hist[n_dft] = t;
            n_dft++;
        }

        int do_diag  = (n % diag_every == 0);
        int do_print = (n % print_every == 0);

        if (do_diag || do_print) {
            Diag d = compute_diag(core_radius);

            if (do_diag)
                fprintf(f1, "%.4f\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t"
                            "%.6f\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\n",
                        t, d.Et, d.Ek, d.Eg, d.Em, d.Ep,
                        d.fc, d.peak[0], d.peak[1], d.peak[2], d.peak_P, d.rho_center, d.Pz);

            if (do_print) {
                double elapsed = omp_get_wtime() - wall_start;
                double frac = (double)n / Nt;
                double eta_t = (frac > 0.001) ? elapsed * (1.0-frac)/frac : 0;
                printf("  t=%7.1f  E=%.2f  fc=%.4f  |P|=%.6f  Pz=%.4f  [%.0fs, ETA %.0fs]\n",
                       t, d.Et, d.fc, d.peak_P, d.Pz, elapsed, eta_t);
                fflush(stdout);
            }
        }

        if (n == Nt) break;
        verlet_step();
    }

    fclose(f1);

    /* Final state */
    Diag dfinal = compute_diag(core_radius);
    double elapsed1 = omp_get_wtime() - wall_start;

    printf("\nEvolution complete (%.1f sec)\n", elapsed1);
    printf("  Final: E=%.2f  fc=%.4f  pk=(%.4f,%.4f,%.4f)  |P|=%.6f  Pz=%.4f\n",
           dfinal.Et, dfinal.fc, dfinal.peak[0], dfinal.peak[1], dfinal.peak[2],
           dfinal.peak_P, dfinal.Pz);

    int survived = (dfinal.fc > 0.01 && dfinal.peak_P > 1e-6);
    printf("  Survived? %s (fc=%.4f, |P|max=%.6f)\n",
           survived ? "YES" : "NO", dfinal.fc, dfinal.peak_P);

    /* Phase 2: Non-Breathing Verification */
    printf("\n--- Phase 2: Non-Breathing Verification ---\n");

    if (n_dft < 100) {
        printf("  Not enough DFT samples (%d), skipping.\n", n_dft);
    } else {
        int start = n_dft / 2;

        double peak_power_rho = 0;
        double omega_rho = find_peak_omega(rho_hist, t_hist, n_dft, start, 5.0, &peak_power_rho);

        double peak_power_phi = 0;
        double omega_phi = find_peak_omega(phi0_hist, t_hist, n_dft, start, 5.0, &peak_power_phi);

        double mean_rho = 0;
        for (int j = start; j < n_dft; j++) mean_rho += rho_hist[j];
        mean_rho /= (n_dft - start);

        double var_rho = 0;
        for (int j = start; j < n_dft; j++) {
            double dr = rho_hist[j] - mean_rho;
            var_rho += dr * dr;
        }
        var_rho /= (n_dft - start);

        printf("  rho(center): mean=%.4e, variance=%.4e, relative var=%.4e\n",
               mean_rho, var_rho,
               (fabs(mean_rho) > 1e-30) ? var_rho / (mean_rho * mean_rho) : 0.0);
        printf("  Peak omega (rho_center): %.4f (power=%.4e)\n", omega_rho, peak_power_rho);
        printf("  Peak omega (phi_0_center): %.4f (power=%.4e)\n", omega_phi, peak_power_phi);

        int is_breathing = (omega_rho > 0.1 && var_rho / (mean_rho * mean_rho + 1e-30) > 0.01);
        printf("  Breathing? %s\n", is_breathing ? "YES" : "NO (static or dispersed)");

        /* Write DFT data */
        snprintf(path, sizeof(path), "%s/dynA_%s_dft.tsv", outdir, label);
        FILE *fdft = fopen(path, "w");
        if (fdft) {
            fprintf(fdft, "# DFT of rho(center,t) second half\n");
            fprintf(fdft, "omega\tpower_rho\tpower_phi\n");
            double mean_phi0 = 0;
            for (int j = start; j < n_dft; j++) mean_phi0 += phi0_hist[j];
            mean_phi0 /= (n_dft - start);

            int nf = 500;
            for (int kk = 1; kk < nf; kk++) {
                double omega = 5.0 * kk / nf;
                double re_r = 0, im_r = 0, re_p = 0, im_p = 0;
                for (int j = start; j < n_dft; j++) {
                    double dtj = (j > start) ?
                        (t_hist[j]-t_hist[j-1]) : (t_hist[start+1]-t_hist[start]);
                    double c = cos(omega * t_hist[j]);
                    double s = sin(omega * t_hist[j]);
                    re_r += (rho_hist[j] - mean_rho) * c * dtj;
                    im_r += (rho_hist[j] - mean_rho) * s * dtj;
                    re_p += (phi0_hist[j] - mean_phi0) * c * dtj;
                    im_p += (phi0_hist[j] - mean_phi0) * s * dtj;
                }
                fprintf(fdft, "%.6f\t%.6e\t%.6e\n", omega, re_r*re_r+im_r*im_r, re_p*re_p+im_p*im_p);
            }
            fclose(fdft);
        }

        /* Write rho time history */
        snprintf(path, sizeof(path), "%s/dynA_%s_rho_history.tsv", outdir, label);
        FILE *frho = fopen(path, "w");
        if (frho) {
            fprintf(frho, "time\trho_center\tphi0_center\n");
            for (int j = 0; j < n_dft; j++)
                fprintf(frho, "%.4f\t%.6e\t%.6e\n", t_hist[j], rho_hist[j], phi0_hist[j]);
            fclose(frho);
        }
    }

    /* Phase 3: Strain Multipoles */
    printf("\n--- Phase 3: Strain Multipoles ---\n");
    compute_strain_multipoles(8.0, label);
    printf("Phase 3 complete.\n");
    fflush(stdout);

    free(rho_hist);
    free(phi0_hist);
    free(t_hist);
}

/* ─── Main ─── */
int main(int argc, char **argv)
{
    /* Parse args */
    for (int i = 1; i < argc - 1; i += 2) {
        if      (!strcmp(argv[i], "-mu"))       mu_pot    = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-kappa"))    kappa     = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-mass"))     mass      = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-A"))        A0        = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-N"))        N         = atoi(argv[i+1]);
        else if (!strcmp(argv[i], "-L"))        L         = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-tfinal"))   tfinal    = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-cfl"))      cfl_frac  = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-Rtube"))    R_tube    = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-o"))        strncpy(outdir, argv[i+1], sizeof(outdir)-1);
        else { fprintf(stderr, "Unknown flag: %s\n", argv[i]); exit(1); }
    }

    dx  = 2.0 * L / (N - 1);
    dx2 = dx * dx;
    m2  = mass * mass;
    dt  = cfl_frac * dx;
    Ngrid = (long)N * N * N;

    double Lz = 2.0 * L;
    double k_wave = 2.0 * M_PI / Lz;
    double omega = sqrt(k_wave * k_wave + mass * mass);

    printf("=== V26-DynA: Propagating Helical Wave Braid ===\n");
    printf("Parameters:\n");
    printf("  mu=%.1f  kappa=%.1f  mass=%.3f  A0=%.3f  R_tube=%.1f\n",
           mu_pot, kappa, mass, A0, R_tube);
    printf("  k=%.6f  omega=%.6f  v_g=%.4f\n", k_wave, omega, k_wave/omega);
    printf("  N=%d  L=%.1f  dx=%.4f  dt=%.5f\n", N, L, dx, dt);
    printf("  Ngrid=%ld (%.1f M)  Memory: %.1f MB\n",
           Ngrid, Ngrid/1e6, Ngrid*8.0*10/1e6);
    printf("  Threads: %d\n", omp_get_max_threads());
    printf("  tfinal=%.0f  cfl=%.2f\n", tfinal, cfl_frac);
    printf("  BC: periodic z, absorbing x,y\n");
    printf("  Couplings: mu=%.1f, kappa=%.1f (triple product ONLY)\n", mu_pot, kappa);
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

    /* Run 1: Propagating braid */
    run_braid(1);

    /* Run 2: Static braid control */
    run_braid(0);

    /* Cleanup */
    for (int a = 0; a < 3; a++) {
        free(phi[a]); free(vel[a]); free(acc[a]);
    }
    free(damp);

    printf("\n=== V26-DynA Complete ===\n");
    return 0;
}
