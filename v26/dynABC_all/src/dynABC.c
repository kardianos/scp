/*
 * dynABC.c — V26-DynABC: Full Dynamic Braid (Propagating + Rotating + Massless)
 *
 * The ULTIMATE test: m=0, propagating at c along z, rotating at Omega around z.
 * Three braided displacement fields with triple-product potential only.
 * "Particle as process" — does a dynamical, massless, rotating, propagating
 * braided configuration maintain itself through continuous motion?
 *
 * Lagrangian:
 *   L = Sum_a [1/2(dt phi_a)^2 - 1/2|grad phi_a|^2] - (mu/2)P^2/(1+kappa*P^2)
 *
 * Initialization:
 *   phase_a = theta + k*z + 2*pi*a/3
 *   phi_a(x,0) = A(r_perp) * cos(phase_a)
 *   vel_a(x,0) = (k + Omega) * A(r_perp) * sin(phase_a)
 *
 * Two runs: DYNAMIC (with velocities) and STATIC control (vel=0).
 *
 * Measurements:
 *   1. Survival (fc, |P|)
 *   2. Non-breathing (DFT of center density)
 *   3. l=2 fraction (strain multipoles on shell R=8)
 *   4. Propagation speed (track phase along z)
 *   5. Angular momentum conservation (L_z)
 *
 * Compile: gcc -O3 -fopenmp -Wall -o dynABC src/dynABC.c -lm
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

/* --- Parameters --- */
static double mu_pot    = -20.0;
static double kappa_pot = 20.0;
static double mass      = 0.0;     /* MASSLESS */
static double A0        = 0.8;
static double R_tube    = 3.0;
static double Omega     = 0.1;     /* angular rotation frequency */

static int    N         = 128;
static double L         = 20.0;
static double tfinal    = 500.0;
static double cfl_frac  = 0.20;
static char   outdir[512] = "data";

static void parse_args(int argc, char **argv)
{
    for (int i = 1; i < argc - 1; i += 2) {
        if      (!strcmp(argv[i], "-mu"))       mu_pot    = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-kappa"))    kappa_pot = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-mass"))     mass      = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-A"))        A0        = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-N"))        N         = atoi(argv[i+1]);
        else if (!strcmp(argv[i], "-L"))        L         = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-tfinal"))   tfinal    = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-cfl"))      cfl_frac  = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-Rtube"))    R_tube    = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-Omega"))    Omega     = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-o"))        strncpy(outdir, argv[i+1], sizeof(outdir)-1);
        else { fprintf(stderr, "Unknown flag: %s\n", argv[i]); exit(1); }
    }
}

/* --- Index helpers --- */
#define IDX(i,j,k) ((long)(i)*N*N + (long)(j)*N + (long)(k))

/* --- Globals --- */
static double *phi[3], *vel[3], *acc[3];
static double *damp;
static double dx, dx2, m2, dt;
static long Ngrid;
static double Lz;  /* domain length in z = 2*L */
static double k_axial;  /* axial wave number = 2*pi/Lz */

/* Periodic z BC */
static inline int wrap_z(int k)
{
    if (k < 0)  return k + N;
    if (k >= N) return k - N;
    return k;
}

/* --- Initialization: Dynamic braid (propagating + rotating) --- */
static void init_dynamic_braid(int set_velocity)
{
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
            double phase = theta + k_axial * z + 2.0 * M_PI * a / 3.0;
            phi[a][idx] = envelope * cos(phase);
            if (set_velocity) {
                /* vel = d(phi)/dt for wave propagating at c along z and rotating at Omega */
                /* phi ~ envelope * cos(theta + k*z - (k+Omega)*t + offset) */
                /* d/dt phi = (k+Omega) * envelope * sin(phase) at t=0 */
                vel[a][idx] = (k_axial + Omega) * envelope * sin(phase);
            } else {
                vel[a][idx] = 0.0;
            }
        }
    }
}

/* --- Set absorbing damping (cylindrical: damp in x,y only, periodic z) --- */
static void set_damping(void)
{
    double R_abs_inner = L * 0.70;
    double R_abs_outer = L * 0.95;

    #pragma omp parallel for schedule(static)
    for (long idx = 0; idx < Ngrid; idx++) {
        int ii = idx / (N * N);
        int jj = (idx / N) % N;
        double x = -L + ii * dx;
        double y = -L + jj * dx;
        double rperp = sqrt(x * x + y * y);
        if (rperp > R_abs_inner) {
            double f = (rperp - R_abs_inner) / (R_abs_outer - R_abs_inner);
            if (f > 1.0) f = 1.0;
            damp[idx] = 1.0 - 0.98 * f * f;
        } else {
            damp[idx] = 1.0;
        }
    }
}

/* --- Compute acceleration (massless, triple-product only, periodic z) --- */
static void compute_acc(void)
{
    for (int a = 0; a < 3; a++) {
        #pragma omp parallel for schedule(static)
        for (int i = 2; i < N-2; i++) {
            for (int j = 2; j < N-2; j++) {
                for (int k = 0; k < N; k++) {
                    long idx = IDX(i, j, k);

                    int km1 = wrap_z(k-1);
                    int kp1 = wrap_z(k+1);

                    double lapl = (phi[a][IDX(i+1,j,k)] + phi[a][IDX(i-1,j,k)]
                                 + phi[a][IDX(i,j+1,k)] + phi[a][IDX(i,j-1,k)]
                                 + phi[a][IDX(i,j,kp1)] + phi[a][IDX(i,j,km1)]
                                 - 6.0 * phi[a][idx]) / dx2;

                    /* Triple-product potential force */
                    double p0 = phi[0][idx], p1 = phi[1][idx], p2 = phi[2][idx];
                    double P = p0 * p1 * p2;
                    double P2 = P * P;
                    double denom2 = (1.0 + kappa_pot * P2);
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

        /* Boundary: zero acceleration at x,y edges */
        #pragma omp parallel for schedule(static)
        for (long idx = 0; idx < Ngrid; idx++) {
            int ii = idx / (N * N);
            int jj = (idx / N) % N;
            if (ii < 2 || ii >= N-2 || jj < 2 || jj >= N-2)
                acc[a][idx] = 0.0;
        }
    }
}

/* --- Velocity Verlet step --- */
static void verlet_step(void)
{
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

/* --- Diagnostics structure --- */
typedef struct {
    double Ek, Eg, Em, Ep, Et;
    double peak[3];
    double fc;
    double peak_P;
    double rho_center;
    double Lz_angular;     /* z-component of angular momentum */
    double Pz_linear;      /* z-component of linear momentum */
} Diag;

/* --- Compute diagnostics --- */
static Diag compute_diag(double core_radius)
{
    Diag d;
    memset(&d, 0, sizeof(d));

    double Ek=0, Eg=0, Em=0, Ep=0;
    double Ecore=0, Eall=0;
    double peak[3] = {0,0,0};
    double peak_P = 0;
    double Lz_ang = 0;
    double Pz_lin = 0;

    /* Center density */
    {
        long ic = IDX(N/2, N/2, N/2);
        int ci = N/2, cj = N/2, ck = N/2;
        int km1 = wrap_z(ck-1), kp1 = wrap_z(ck+1);
        double rho = 0;
        for (int a = 0; a < 3; a++) {
            rho += 0.5 * vel[a][ic] * vel[a][ic];
            double gx = (phi[a][IDX(ci+1,cj,ck)] - phi[a][IDX(ci-1,cj,ck)]) / (2*dx);
            double gy = (phi[a][IDX(ci,cj+1,ck)] - phi[a][IDX(ci,cj-1,ck)]) / (2*dx);
            double gz = (phi[a][IDX(ci,cj,kp1)] - phi[a][IDX(ci,cj,km1)]) / (2*dx);
            rho += 0.5 * (gx*gx + gy*gy + gz*gz);
            rho += 0.5 * m2 * phi[a][ic] * phi[a][ic];
        }
        double p0 = phi[0][ic], p1 = phi[1][ic], p2 = phi[2][ic];
        double P = p0 * p1 * p2;
        double P2 = P * P;
        rho += 0.5 * mu_pot * P2 / (1.0 + kappa_pot * P2);
        d.rho_center = rho;
    }

    #pragma omp parallel
    {
        double lEk=0, lEg=0, lEm=0, lEp=0, lEc=0, lEa=0;
        double lpk[3] = {0,0,0};
        double lpkP = 0;
        double lLz = 0, lPz = 0;

        #pragma omp for schedule(static) nowait
        for (int i = 2; i < N-2; i++) {
            double x = -L + i * dx;
            for (int j = 2; j < N-2; j++) {
                double y = -L + j * dx;
                for (int k = 2; k < N-2; k++) {
                    double z = -L + k * dx;
                    long idx = IDX(i,j,k);
                    double dV = dx * dx * dx;
                    double e_loc = 0;

                    int km1 = wrap_z(k-1), kp1 = wrap_z(k+1);

                    for (int a = 0; a < 3; a++) {
                        double v_a = vel[a][idx];
                        double v2 = v_a * v_a;
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

                        /* Angular momentum: L_z = sum_a integral (x * v_a * gy - y * v_a * gx) dV
                         * Actually L_z = sum_a integral v_a * (x * d_y phi_a - y * d_x phi_a) dV
                         * (field angular momentum = integral T^{0i} x_j epsilon_{ijk} dV)
                         * More precisely: L_z = -sum_a integral (dt phi_a)(x dy phi_a - y dx phi_a) dV */
                        lLz += -v_a * (x * gy - y * gx) * dV;

                        /* Linear momentum along z: P_z = -sum_a integral (dt phi_a)(dz phi_a) dV */
                        lPz += -v_a * gz * dV;
                    }

                    double p0 = phi[0][idx], p1 = phi[1][idx], p2 = phi[2][idx];
                    double P = p0 * p1 * p2;
                    double P2 = P * P;
                    double Vloc = 0.5 * mu_pot * P2 / (1.0 + kappa_pot * P2);
                    lEp += Vloc * dV;
                    e_loc += Vloc;

                    double absP = fabs(P);
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
            for (int a = 0; a < 3; a++)
                if (lpk[a] > peak[a]) peak[a] = lpk[a];
            if (lpkP > peak_P) peak_P = lpkP;
            Lz_ang += lLz;
            Pz_lin += lPz;
        }
    }

    d.Ek = Ek; d.Eg = Eg; d.Em = Em; d.Ep = Ep;
    d.Et = Ek + Eg + Em + Ep;
    for (int a = 0; a < 3; a++) d.peak[a] = peak[a];
    d.fc = (Eall > 1e-20) ? Ecore / Eall : 0.0;
    d.peak_P = peak_P;
    d.Lz_angular = Lz_ang;
    d.Pz_linear = Pz_lin;

    return d;
}

/* --- Fast center density (no full grid scan) --- */
static double center_rho(void)
{
    long ic = IDX(N/2, N/2, N/2);
    int ci = N/2, cj = N/2, ck = N/2;
    int km1 = wrap_z(ck-1), kp1 = wrap_z(ck+1);
    double rho = 0;
    for (int a = 0; a < 3; a++) {
        rho += 0.5 * vel[a][ic] * vel[a][ic];
        double gx = (phi[a][IDX(ci+1,cj,ck)] - phi[a][IDX(ci-1,cj,ck)]) / (2*dx);
        double gy = (phi[a][IDX(ci,cj+1,ck)] - phi[a][IDX(ci,cj-1,ck)]) / (2*dx);
        double gz = (phi[a][IDX(ci,cj,kp1)] - phi[a][IDX(ci,cj,km1)]) / (2*dx);
        rho += 0.5 * (gx*gx + gy*gy + gz*gz);
        rho += 0.5 * m2 * phi[a][ic] * phi[a][ic];
    }
    double p0 = phi[0][ic], p1 = phi[1][ic], p2 = phi[2][ic];
    double P = p0 * p1 * p2;
    double P2 = P * P;
    rho += 0.5 * mu_pot * P2 / (1.0 + kappa_pot * P2);
    return rho;
}

/* --- DFT to find peak frequency --- */
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

/* --- Legendre polynomials --- */
static double legendre_P(int l, double x)
{
    switch (l) {
        case 0: return 1.0;
        case 1: return x;
        case 2: return 0.5*(3.0*x*x - 1.0);
        default: return 0.0;
    }
}

/* --- Compute strain tensor at a grid point --- */
static void compute_strain(int i, int j, int k, double eps[6])
{
    int km1 = wrap_z(k-1), kp1 = wrap_z(k+1);
    double dphi[3][3];
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

/* --- Strain multipoles on spherical shell --- */
static void compute_strain_multipoles(double R_shell, const char *label)
{
    int N_sample = 200;
    double golden_ratio = (1.0 + sqrt(5.0)) / 2.0;

    double c_l[3] = {0, 0, 0};
    double norm = 0;
    int n_valid = 0;

    char path[600];
    snprintf(path, sizeof(path), "%s/dynABC_%s_strain.tsv", outdir, label);
    FILE *fstrain = fopen(path, "w");
    snprintf(path, sizeof(path), "%s/dynABC_%s_multipoles.tsv", outdir, label);
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

        int ig = (int)((x + L) / dx + 0.5);
        int jg = (int)((y + L) / dx + 0.5);
        int kg = (int)((z + L) / dx + 0.5);

        if (ig < 2 || ig >= N-2 || jg < 2 || jg >= N-2 || kg < 2 || kg >= N-2)
            continue;

        double eps[6];
        compute_strain(ig, jg, kg, eps);

        double theta_tr = eps[0] + eps[1] + eps[2];

        double sig[6];
        sig[0] = eps[0] - theta_tr / 3.0;
        sig[1] = eps[1] - theta_tr / 3.0;
        sig[2] = eps[2] - theta_tr / 3.0;
        sig[3] = eps[3]; sig[4] = eps[4]; sig[5] = eps[5];

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

/* --- Track propagation: find z-position of phi_0 maximum at r_perp=0 --- */
static double track_phase_z(void)
{
    /* Scan along z-axis at (x,y) = center (i=N/2, j=N/2) */
    int ic = N/2, jc = N/2;
    double max_phi = -1e30;
    int k_max = 0;
    for (int k = 0; k < N; k++) {
        double val = phi[0][IDX(ic, jc, k)];
        if (val > max_phi) {
            max_phi = val;
            k_max = k;
        }
    }
    /* Parabolic interpolation for sub-grid precision */
    int km = wrap_z(k_max - 1);
    int kp = wrap_z(k_max + 1);
    double f_m = phi[0][IDX(ic, jc, km)];
    double f_0 = phi[0][IDX(ic, jc, k_max)];
    double f_p = phi[0][IDX(ic, jc, kp)];
    double denom = 2.0 * (2.0 * f_0 - f_m - f_p);
    double shift = 0.0;
    if (fabs(denom) > 1e-30)
        shift = (f_m - f_p) / denom;
    return -L + (k_max + shift) * dx;
}

/* --- Run one configuration (dynamic or static) --- */
static void run_config(const char *label, int set_velocity)
{
    printf("\n================================================================\n");
    printf("  %s: %s braid  (m=%.1f, mu=%.1f, kappa=%.1f, Omega=%.3f)\n",
           label, set_velocity ? "DYNAMIC" : "STATIC", mass, mu_pot, kappa_pot, Omega);
    printf("  k_axial=%.4f, v_init=(k+Omega)=%.4f\n", k_axial, k_axial + Omega);
    printf("================================================================\n\n");

    /* Zero fields */
    for (int a = 0; a < 3; a++) {
        memset(phi[a], 0, Ngrid * sizeof(double));
        memset(vel[a], 0, Ngrid * sizeof(double));
        memset(acc[a], 0, Ngrid * sizeof(double));
    }

    init_dynamic_braid(set_velocity);
    set_damping();

    double core_radius = 8.0;
    int Nt = (int)(tfinal / dt) + 1;

    /* DFT storage */
    int max_dft = 50000;
    double *rho_hist  = malloc(max_dft * sizeof(double));
    double *phi0_hist = malloc(max_dft * sizeof(double));
    double *t_hist    = malloc(max_dft * sizeof(double));
    int n_dft = 0;
    int dft_every = Nt / max_dft;
    if (dft_every < 1) dft_every = 1;

    int diag_every = Nt / 200;
    if (diag_every < 1) diag_every = 1;
    int print_every = Nt / 20;
    if (print_every < 1) print_every = 1;

    /* Phase tracking storage */
    int max_phase = 5000;
    double *z_phase = malloc(max_phase * sizeof(double));
    double *t_phase = malloc(max_phase * sizeof(double));
    int n_phase = 0;
    int phase_every = Nt / max_phase;
    if (phase_every < 1) phase_every = 1;

    /* Output file */
    char path[600];
    snprintf(path, sizeof(path), "%s/dynABC_%s_evolution.tsv", outdir, label);
    FILE *fout = fopen(path, "w");
    if (!fout) {
        fprintf(stderr, "Cannot open %s\n", path);
        free(rho_hist); free(phi0_hist); free(t_hist);
        free(z_phase); free(t_phase);
        return;
    }
    fprintf(fout, "time\tE_total\tE_kin\tE_grad\tE_mass\tE_pot\t"
                  "fc\tpeak0\tpeak1\tpeak2\tpeak_P\trho_center\t"
                  "Lz_angular\tPz_linear\tz_phase\n");

    compute_acc();

    double wall_start = omp_get_wtime();

    /* Initial diagnostics */
    Diag d0 = compute_diag(core_radius);
    double z0 = track_phase_z();
    printf("  t=%7.1f  E=%.2f  fc=%.4f  pk=(%.3f,%.3f,%.3f)  |P|=%.6f  Lz=%.4f  Pz=%.4f  z0=%.3f\n",
           0.0, d0.Et, d0.fc, d0.peak[0], d0.peak[1], d0.peak[2],
           d0.peak_P, d0.Lz_angular, d0.Pz_linear, z0);

    double Lz_initial = d0.Lz_angular;
    double Pz_initial = d0.Pz_linear;
    double E_initial  = d0.Et;

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

        /* Phase tracking */
        if (n % phase_every == 0 && n_phase < max_phase) {
            z_phase[n_phase] = track_phase_z();
            t_phase[n_phase] = t;
            n_phase++;
        }

        int do_diag  = (n % diag_every == 0);
        int do_print = (n % print_every == 0);

        if (do_diag || do_print) {
            Diag d = compute_diag(core_radius);
            double zp = track_phase_z();

            if (do_diag) {
                fprintf(fout, "%.4f\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t"
                              "%.6f\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t"
                              "%.6e\t%.6e\t%.6e\n",
                        t, d.Et, d.Ek, d.Eg, d.Em, d.Ep,
                        d.fc, d.peak[0], d.peak[1], d.peak[2],
                        d.peak_P, d.rho_center,
                        d.Lz_angular, d.Pz_linear, zp);
            }

            if (do_print) {
                double elapsed = omp_get_wtime() - wall_start;
                double frac = (double)n / Nt;
                double eta_t = (frac > 0.001) ? elapsed * (1.0-frac)/frac : 0;
                printf("  t=%7.1f  E=%.2f  fc=%.4f  pk=(%.3f,%.3f,%.3f)  |P|=%.6f  Lz=%.4f  Pz=%.4f  [%.0fs, ETA %.0fs]\n",
                       t, d.Et, d.fc, d.peak[0], d.peak[1], d.peak[2],
                       d.peak_P, d.Lz_angular, d.Pz_linear, elapsed, eta_t);
                fflush(stdout);
            }
        }

        if (n == Nt) break;
        verlet_step();
    }

    fclose(fout);

    /* Final diagnostics */
    Diag dfinal = compute_diag(core_radius);
    double elapsed1 = omp_get_wtime() - wall_start;

    printf("\n--- %s: Evolution complete (%.1f sec) ---\n", label, elapsed1);
    printf("  Final: E=%.2f  fc=%.4f  pk=(%.3f,%.3f,%.3f)  |P|max=%.6f\n",
           dfinal.Et, dfinal.fc, dfinal.peak[0], dfinal.peak[1], dfinal.peak[2], dfinal.peak_P);
    printf("  Lz: initial=%.4f  final=%.4f  ratio=%.4f\n",
           Lz_initial, dfinal.Lz_angular,
           (fabs(Lz_initial) > 1e-20) ? dfinal.Lz_angular / Lz_initial : 0.0);
    printf("  Pz: initial=%.4f  final=%.4f  ratio=%.4f\n",
           Pz_initial, dfinal.Pz_linear,
           (fabs(Pz_initial) > 1e-20) ? dfinal.Pz_linear / Pz_initial : 0.0);
    printf("  Energy: initial=%.4f  final=%.4f  ratio=%.4f\n",
           E_initial, dfinal.Et,
           (fabs(E_initial) > 1e-20) ? dfinal.Et / E_initial : 0.0);

    /* SURVIVAL test */
    int survived = (dfinal.fc > 0.3 && dfinal.peak_P > 0.01);
    printf("  SURVIVED? %s (fc=%.4f > 0.3?, |P|=%.6f > 0.01?)\n",
           survived ? "YES" : "NO", dfinal.fc, dfinal.peak_P);

    /* ===== Non-Breathing Analysis (DFT) ===== */
    printf("\n--- %s: Non-Breathing Analysis ---\n", label);

    double omega_rho = 0, peak_power_rho = 0;
    double omega_phi = 0, peak_power_phi = 0;
    double mean_rho = 0, var_rho = 0;

    if (n_dft >= 100) {
        int start = n_dft / 2;

        omega_rho = find_peak_omega(rho_hist, t_hist, n_dft, start, 5.0, &peak_power_rho);
        omega_phi = find_peak_omega(phi0_hist, t_hist, n_dft, start, 5.0, &peak_power_phi);

        for (int j = start; j < n_dft; j++) mean_rho += rho_hist[j];
        mean_rho /= (n_dft - start);

        for (int j = start; j < n_dft; j++) {
            double dr = rho_hist[j] - mean_rho;
            var_rho += dr * dr;
        }
        var_rho /= (n_dft - start);

        double rel_var = (fabs(mean_rho) > 1e-30) ? var_rho / (mean_rho * mean_rho) : 0.0;
        printf("  rho(center): mean=%.4e, variance=%.4e, relative var=%.4e\n",
               mean_rho, var_rho, rel_var);
        printf("  Peak omega (rho_center): %.4f (power=%.4e)\n", omega_rho, peak_power_rho);
        printf("  Peak omega (phi_0_center): %.4f (power=%.4e)\n", omega_phi, peak_power_phi);

        int is_breathing = (omega_rho > 0.1 && var_rho / (mean_rho * mean_rho + 1e-30) > 0.01);
        printf("  Breathing? %s\n", is_breathing ? "YES (oscillon-like)" : "NO (non-breathing)");

        /* Write DFT data */
        snprintf(path, sizeof(path), "%s/dynABC_%s_dft.tsv", outdir, label);
        FILE *fdft = fopen(path, "w");
        if (fdft) {
            fprintf(fdft, "# DFT of rho(center,t) and phi0(center,t), second half\n");
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
        snprintf(path, sizeof(path), "%s/dynABC_%s_rho_history.tsv", outdir, label);
        FILE *frho = fopen(path, "w");
        if (frho) {
            fprintf(frho, "time\trho_center\tphi0_center\n");
            for (int j = 0; j < n_dft; j++)
                fprintf(frho, "%.4f\t%.6e\t%.6e\n", t_hist[j], rho_hist[j], phi0_hist[j]);
            fclose(frho);
        }
    } else {
        printf("  Not enough DFT samples (%d), skipping.\n", n_dft);
    }

    /* ===== Strain Multipoles ===== */
    printf("\n--- %s: Strain Multipoles (R=8) ---\n", label);
    compute_strain_multipoles(8.0, label);

    /* ===== Propagation Speed ===== */
    printf("\n--- %s: Propagation Speed ---\n", label);
    if (n_phase >= 10) {
        /* Fit propagation speed from phase tracking via linear regression */
        /* Use second half to avoid transient */
        int start = n_phase / 2;
        int count = n_phase - start;

        /* Handle z wrapping: unwrap phase position */
        double *z_unwrap = malloc(n_phase * sizeof(double));
        z_unwrap[0] = z_phase[0];
        for (int j = 1; j < n_phase; j++) {
            double dz = z_phase[j] - z_phase[j-1];
            /* Unwrap: if jump > L, subtract 2L; if < -L, add 2L */
            if (dz > L) dz -= 2.0 * L;
            if (dz < -L) dz += 2.0 * L;
            z_unwrap[j] = z_unwrap[j-1] + dz;
        }

        /* Linear regression: z = v_prop * t + b */
        double sum_t = 0, sum_z = 0, sum_tz = 0, sum_t2 = 0;
        for (int j = start; j < n_phase; j++) {
            sum_t  += t_phase[j];
            sum_z  += z_unwrap[j];
            sum_tz += t_phase[j] * z_unwrap[j];
            sum_t2 += t_phase[j] * t_phase[j];
        }
        double v_prop = (count * sum_tz - sum_t * sum_z) / (count * sum_t2 - sum_t * sum_t);
        double z_mean = sum_z / count;
        double t_mean = sum_t / count;

        /* R^2 goodness of fit */
        double ss_res = 0, ss_tot = 0;
        double b = z_mean - v_prop * t_mean;
        for (int j = start; j < n_phase; j++) {
            double pred = v_prop * t_phase[j] + b;
            ss_res += (z_unwrap[j] - pred) * (z_unwrap[j] - pred);
            ss_tot += (z_unwrap[j] - z_mean) * (z_unwrap[j] - z_mean);
        }
        double R2 = (ss_tot > 1e-30) ? 1.0 - ss_res / ss_tot : 0.0;

        printf("  v_propagation = %.6f (expected %.6f = k*c/|k|)\n",
               v_prop, k_axial > 0 ? 1.0 : -1.0);
        printf("  R^2 = %.6f\n", R2);
        printf("  z(0) = %.3f, z(end) = %.3f, Delta_z = %.3f\n",
               z_unwrap[0], z_unwrap[n_phase-1], z_unwrap[n_phase-1] - z_unwrap[0]);

        /* Write phase tracking data */
        snprintf(path, sizeof(path), "%s/dynABC_%s_phase_track.tsv", outdir, label);
        FILE *fphase = fopen(path, "w");
        if (fphase) {
            fprintf(fphase, "time\tz_phase\tz_unwrap\n");
            for (int j = 0; j < n_phase; j++)
                fprintf(fphase, "%.4f\t%.6f\t%.6f\n", t_phase[j], z_phase[j], z_unwrap[j]);
            fclose(fphase);
        }

        free(z_unwrap);
    } else {
        printf("  Not enough phase tracking samples (%d)\n", n_phase);
    }

    /* ===== Angular Momentum Conservation ===== */
    printf("\n--- %s: Angular Momentum Conservation ---\n", label);
    printf("  L_z(0) = %.6f\n", Lz_initial);
    printf("  L_z(end) = %.6f\n", dfinal.Lz_angular);
    if (fabs(Lz_initial) > 1e-10) {
        printf("  Delta L_z / L_z(0) = %.6f\n",
               (dfinal.Lz_angular - Lz_initial) / Lz_initial);
    }
    printf("  P_z(0) = %.6f\n", Pz_initial);
    printf("  P_z(end) = %.6f\n", dfinal.Pz_linear);

    free(rho_hist);
    free(phi0_hist);
    free(t_hist);
    free(z_phase);
    free(t_phase);
}

/* --- Main --- */
int main(int argc, char **argv)
{
    parse_args(argc, argv);

    dx  = 2.0 * L / (N - 1);
    dx2 = dx * dx;
    m2  = mass * mass;
    dt  = cfl_frac * dx;
    Ngrid = (long)N * N * N;
    Lz  = 2.0 * L;
    k_axial = 2.0 * M_PI / Lz;

    printf("=== V26-DynABC: Full Dynamic Braid (Propagating + Rotating + Massless) ===\n");
    printf("Parameters:\n");
    printf("  mu=%.1f  kappa=%.1f  mass=%.3f  A0=%.3f  R_tube=%.1f\n",
           mu_pot, kappa_pot, mass, A0, R_tube);
    printf("  Omega=%.4f  k_axial=%.4f  v_init=(k+Omega)=%.4f\n",
           Omega, k_axial, k_axial + Omega);
    printf("  N=%d  L=%.1f  dx=%.4f  dt=%.5f\n", N, L, dx, dt);
    printf("  Ngrid=%ld (%.1f M)  Memory: %.1f MB\n",
           Ngrid, Ngrid/1e6, Ngrid*8.0*10/1e6);
    printf("  Threads: %d\n", omp_get_max_threads());
    printf("  tfinal=%.0f  Nt=%d\n", tfinal, (int)(tfinal / dt) + 1);
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

    /* ============================================
     * Run 1: DYNAMIC braid (propagating + rotating)
     * ============================================ */
    run_config("dynamic", 1);

    /* ============================================
     * Run 2: STATIC control (same phase, vel=0)
     * ============================================ */
    run_config("static", 0);

    /* ============================================
     * Summary comparison
     * ============================================ */
    printf("\n================================================================\n");
    printf("  SUMMARY: DYNAMIC vs STATIC comparison\n");
    printf("  See data files for detailed comparison.\n");
    printf("================================================================\n");

    /* Cleanup */
    for (int a = 0; a < 3; a++) {
        free(phi[a]); free(vel[a]); free(acc[a]);
    }
    free(damp);

    printf("\n=== V26-DynABC Complete ===\n");
    return 0;
}
