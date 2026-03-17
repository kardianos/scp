/*
 * v26_phase4.c — Torsion Constraint and Teleparallel Gravity
 *
 * Based on v26.c mode 2 (twisted tube, triple-product only, NO mass term).
 * Adds symmetric strain (kappa_S) and antisymmetric torsion (kappa_T)
 * coupling terms that decompose the cross-gradient into:
 *   epsilon_{ij} = 1/2(d_i phi_j + d_j phi_i)  [strain, spin-2]
 *   omega_{ij}   = 1/2(d_i phi_j - d_j phi_i)  [torsion, spin-1]
 *
 * Phase 4a: Torsion scan (kappa_T = 0..2, kappa_S = 0)
 * Phase 4b: Joint (kappa_S, kappa_T) scan
 * Phase 4c: Torsion flux quantization
 * Phase 4d: Self-consistent teleparallel metric
 *
 * Compile: gcc -O3 -fopenmp -Wall -o v26p4 src/v26_phase4.c -lm
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

/* --- Parameters --- */
static double mu_pot    = -20.0;
static double kappa_pot = 20.0;
static double mass      = 1.0;    /* mass term (v26 mode 2 keeps m=1) */
static double A0        = 0.8;
static double R_tube    = 3.0;

static double kappa_S   = 0.0;   /* strain coupling */
static double kappa_T   = 0.0;   /* torsion coupling */

static int    N         = 128;
static double L         = 20.0;
static double tfinal    = 500.0;
static double cfl_frac  = 0.20;
static int    phase     = 0;     /* 0=4a, 1=4b, 2=4c, 3=4d */
static char   outdir[512] = "data";

/* Phase 4d parameters */
static double alpha_g_max = 1.0;
static double alpha_g_ramp_time = 100.0;

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
        else if (!strcmp(argv[i], "-kS"))       kappa_S   = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-kT"))       kappa_T   = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-phase"))    phase     = atoi(argv[i+1]);
        else if (!strcmp(argv[i], "-o"))        strncpy(outdir, argv[i+1], sizeof(outdir)-1);
        else if (!strcmp(argv[i], "-ag"))       alpha_g_max = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-agramp"))   alpha_g_ramp_time = atof(argv[i+1]);
        else { fprintf(stderr, "Unknown flag: %s\n", argv[i]); exit(1); }
    }
}

/* --- Index helpers --- */
#define IDX(i,j,k) ((long)(i)*N*N + (long)(j)*N + (long)(k))

/* --- Globals --- */
static double *phi[3], *vel[3], *acc_arr[3];
static double *damp;
static double dx, dx2, m2, dt;
static long Ngrid;

/* Periodic z BC */
static inline int wrap_z(int k)
{
    if (k < 0)  return k + N;
    if (k >= N) return k - N;
    return k;
}

/* --- Initialization: Twisted Tube (mode 2 from v26.c) --- */
static void init_twisted_tube(void)
{
    double Lz = 2.0 * L;  /* domain length in z */

    #pragma omp parallel for schedule(static)
    for (long idx = 0; idx < Ngrid; idx++) {
        int ii = idx / (N * N);
        int jj = (idx / N) % N;
        int kk = idx % N;
        double x = -L + ii * dx;
        double y = -L + jj * dx;
        double z = -L + kk * dx;

        double r_perp = sqrt(x * x + y * y);
        double A_r = A0 * exp(-r_perp * r_perp / (2.0 * R_tube * R_tube));

        double twist = 2.0 * M_PI * z / Lz;

        phi[0][idx] = A_r * cos(twist);
        phi[1][idx] = A_r * cos(twist + 2.0 * M_PI / 3.0);
        phi[2][idx] = A_r * cos(twist + 4.0 * M_PI / 3.0);
    }
}

/* --- Set absorbing damping (cylindrical: damp in x,y only) --- */
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

/* --- Compute d_a(div phi) using mixed second derivatives --- */
/* div phi = d_x phi_0 + d_y phi_1 + d_z phi_2 */
/* d_a(div phi) = sum_b d_a(d_b phi_b) */
static double compute_d_a_divphi(int a, int i, int j, int k)
{
    double d_a_div = 0.0;

    /* b=0: d_a(d_x phi_0) */
    if (a == 0) {
        d_a_div += (phi[0][IDX(i+1,j,k)] - 2.0*phi[0][IDX(i,j,k)]
                  + phi[0][IDX(i-1,j,k)]) / dx2;
    } else if (a == 1) {
        d_a_div += (phi[0][IDX(i+1,j+1,k)] - phi[0][IDX(i+1,j-1,k)]
                  - phi[0][IDX(i-1,j+1,k)] + phi[0][IDX(i-1,j-1,k)])
                  / (4.0 * dx2);
    } else {
        int km1 = wrap_z(k-1), kp1 = wrap_z(k+1);
        d_a_div += (phi[0][IDX(i+1,j,kp1)] - phi[0][IDX(i+1,j,km1)]
                  - phi[0][IDX(i-1,j,kp1)] + phi[0][IDX(i-1,j,km1)])
                  / (4.0 * dx2);
    }

    /* b=1: d_a(d_y phi_1) */
    if (a == 1) {
        d_a_div += (phi[1][IDX(i,j+1,k)] - 2.0*phi[1][IDX(i,j,k)]
                  + phi[1][IDX(i,j-1,k)]) / dx2;
    } else if (a == 0) {
        d_a_div += (phi[1][IDX(i+1,j+1,k)] - phi[1][IDX(i+1,j-1,k)]
                  - phi[1][IDX(i-1,j+1,k)] + phi[1][IDX(i-1,j-1,k)])
                  / (4.0 * dx2);
    } else {
        int km1 = wrap_z(k-1), kp1 = wrap_z(k+1);
        d_a_div += (phi[1][IDX(i,j+1,kp1)] - phi[1][IDX(i,j+1,km1)]
                  - phi[1][IDX(i,j-1,kp1)] + phi[1][IDX(i,j-1,km1)])
                  / (4.0 * dx2);
    }

    /* b=2: d_a(d_z phi_2) */
    if (a == 2) {
        int km1 = wrap_z(k-1), kp1 = wrap_z(k+1);
        d_a_div += (phi[2][IDX(i,j,kp1)] - 2.0*phi[2][IDX(i,j,k)]
                  + phi[2][IDX(i,j,km1)]) / dx2;
    } else if (a == 0) {
        int km1 = wrap_z(k-1), kp1 = wrap_z(k+1);
        d_a_div += (phi[2][IDX(i+1,j,kp1)] - phi[2][IDX(i+1,j,km1)]
                  - phi[2][IDX(i-1,j,kp1)] + phi[2][IDX(i-1,j,km1)])
                  / (4.0 * dx2);
    } else {
        int km1 = wrap_z(k-1), kp1 = wrap_z(k+1);
        d_a_div += (phi[2][IDX(i,j+1,kp1)] - phi[2][IDX(i,j+1,km1)]
                  - phi[2][IDX(i,j-1,kp1)] + phi[2][IDX(i,j-1,km1)])
                  / (4.0 * dx2);
    }

    return d_a_div;
}

/* --- Compute acceleration (standard mode) --- */
static void compute_acc(void)
{
    /* Combined coefficients:
     * acc_a = (1 + kS/2 + kT/2) * lapl_a + (kS/2 - kT/2) * d_a(div phi) - dV/dphi_a
     * Torsion: +kT/2 * (lapl - d_a_div)
     * Strain:  +kS/2 * (lapl + d_a_div) */
    double lapl_coeff = 1.0 + 0.5 * kappa_S + 0.5 * kappa_T;
    double div_coeff  = 0.5 * (kappa_S - kappa_T);

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

                    /* Torsion + Strain coupling */
                    double d_a_div = 0.0;
                    if (fabs(div_coeff) > 1e-15) {
                        d_a_div = compute_d_a_divphi(a, i, j, k);
                    }

                    acc_arr[a][idx] = lapl_coeff * lapl - m2 * phi[a][idx]
                                   + div_coeff * d_a_div - dVdphi;
                }
            }
        }

        /* Boundary: zero acceleration for x,y edges */
        #pragma omp parallel for schedule(static)
        for (long idx = 0; idx < Ngrid; idx++) {
            int ii = idx / (N * N);
            int jj = (idx / N) % N;
            if (ii < 2 || ii >= N-2 || jj < 2 || jj >= N-2)
                acc_arr[a][idx] = 0.0;
        }
    }
}

/* --- Compute acceleration with self-consistent metric correction (Phase 4d) --- */
static void compute_acc_metric(double alpha_g)
{
    /* First compute standard acceleration */
    compute_acc();

    if (alpha_g < 1e-12) return;

    /* Add metric correction: -2 alpha_g * eps_{ij} * d_i d_j phi_a */
    for (int a = 0; a < 3; a++) {
        #pragma omp parallel for schedule(static)
        for (int i = 2; i < N-2; i++) {
            for (int j = 2; j < N-2; j++) {
                for (int k = 0; k < N; k++) {
                    long idx = IDX(i, j, k);
                    int km1 = wrap_z(k-1), kp1 = wrap_z(k+1);

                    /* Compute strain eps_{ij} = 0.5(d_i phi_j + d_j phi_i) at this point */
                    double dphi[3][3];
                    for (int b = 0; b < 3; b++) {
                        dphi[b][0] = (phi[b][IDX(i+1,j,k)] - phi[b][IDX(i-1,j,k)]) / (2.0*dx);
                        dphi[b][1] = (phi[b][IDX(i,j+1,k)] - phi[b][IDX(i,j-1,k)]) / (2.0*dx);
                        dphi[b][2] = (phi[b][IDX(i,j,kp1)] - phi[b][IDX(i,j,km1)]) / (2.0*dx);
                    }

                    double eps[3][3];
                    for (int p = 0; p < 3; p++)
                        for (int q = 0; q < 3; q++)
                            eps[p][q] = 0.5 * (dphi[q][p] + dphi[p][q]);

                    /* Second derivatives d_i d_j phi_a */
                    double d2phi[3][3];
                    /* d_xx */
                    d2phi[0][0] = (phi[a][IDX(i+1,j,k)] - 2.0*phi[a][idx]
                                 + phi[a][IDX(i-1,j,k)]) / dx2;
                    /* d_yy */
                    d2phi[1][1] = (phi[a][IDX(i,j+1,k)] - 2.0*phi[a][idx]
                                 + phi[a][IDX(i,j-1,k)]) / dx2;
                    /* d_zz */
                    d2phi[2][2] = (phi[a][IDX(i,j,kp1)] - 2.0*phi[a][idx]
                                 + phi[a][IDX(i,j,km1)]) / dx2;
                    /* d_xy */
                    d2phi[0][1] = d2phi[1][0] = (phi[a][IDX(i+1,j+1,k)] - phi[a][IDX(i+1,j-1,k)]
                                                - phi[a][IDX(i-1,j+1,k)] + phi[a][IDX(i-1,j-1,k)])
                                                / (4.0 * dx2);
                    /* d_xz */
                    d2phi[0][2] = d2phi[2][0] = (phi[a][IDX(i+1,j,kp1)] - phi[a][IDX(i+1,j,km1)]
                                                - phi[a][IDX(i-1,j,kp1)] + phi[a][IDX(i-1,j,km1)])
                                                / (4.0 * dx2);
                    /* d_yz */
                    d2phi[1][2] = d2phi[2][1] = (phi[a][IDX(i,j+1,kp1)] - phi[a][IDX(i,j+1,km1)]
                                                - phi[a][IDX(i,j-1,kp1)] + phi[a][IDX(i,j-1,km1)])
                                                / (4.0 * dx2);

                    /* Metric correction: -2 alpha_g sum_{ij} eps_{ij} d_i d_j phi_a */
                    double correction = 0.0;
                    for (int p = 0; p < 3; p++)
                        for (int q = 0; q < 3; q++)
                            correction += eps[p][q] * d2phi[p][q];

                    acc_arr[a][idx] -= 2.0 * alpha_g * correction;
                }
            }
        }
    }
}

/* --- Velocity Verlet step --- */
static void verlet_step_standard(void)
{
    for (int a = 0; a < 3; a++) {
        #pragma omp parallel for schedule(static)
        for (long idx = 0; idx < Ngrid; idx++)
            vel[a][idx] += 0.5 * dt * acc_arr[a][idx];
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
            vel[a][idx] += 0.5 * dt * acc_arr[a][idx];
    }

    for (int a = 0; a < 3; a++) {
        #pragma omp parallel for schedule(static)
        for (long idx = 0; idx < Ngrid; idx++) {
            vel[a][idx] *= damp[idx];
            phi[a][idx] *= damp[idx];
        }
    }
}

/* --- Velocity Verlet step with metric (Phase 4d) --- */
static void verlet_step_metric(double alpha_g)
{
    for (int a = 0; a < 3; a++) {
        #pragma omp parallel for schedule(static)
        for (long idx = 0; idx < Ngrid; idx++)
            vel[a][idx] += 0.5 * dt * acc_arr[a][idx];
    }

    for (int a = 0; a < 3; a++) {
        #pragma omp parallel for schedule(static)
        for (long idx = 0; idx < Ngrid; idx++)
            phi[a][idx] += dt * vel[a][idx];
    }

    compute_acc_metric(alpha_g);

    for (int a = 0; a < 3; a++) {
        #pragma omp parallel for schedule(static)
        for (long idx = 0; idx < Ngrid; idx++)
            vel[a][idx] += 0.5 * dt * acc_arr[a][idx];
    }

    for (int a = 0; a < 3; a++) {
        #pragma omp parallel for schedule(static)
        for (long idx = 0; idx < Ngrid; idx++) {
            vel[a][idx] *= damp[idx];
            phi[a][idx] *= damp[idx];
        }
    }
}

/* --- Diagnostics --- */
typedef struct {
    double Ek, Eg, Ep, Et;
    double peak[3];
    double fc;
    double peak_P;
    double rho_center;
    /* Torsion/strain at center */
    double omega_xy, omega_xz, omega_yz;
    double eps_xx, eps_yy, eps_zz, eps_xy, eps_xz, eps_yz;
    double torsion_flux_z;
} Diag;

/* --- Compute torsion and strain at a grid point --- */
static void compute_torsion_strain(int i, int j, int k,
                                   double omega_out[3], double eps_out[6])
{
    int km1 = wrap_z(k-1), kp1 = wrap_z(k+1);
    double dphi[3][3];
    for (int a = 0; a < 3; a++) {
        dphi[a][0] = (phi[a][IDX(i+1,j,k)] - phi[a][IDX(i-1,j,k)]) / (2.0*dx);
        dphi[a][1] = (phi[a][IDX(i,j+1,k)] - phi[a][IDX(i,j-1,k)]) / (2.0*dx);
        dphi[a][2] = (phi[a][IDX(i,j,kp1)] - phi[a][IDX(i,j,km1)]) / (2.0*dx);
    }

    /* omega_{ij} = 0.5(d_i phi_j - d_j phi_i), vorticity Omega_k = eps_{ijk} omega_{ij} */
    double w_01 = 0.5 * (dphi[1][0] - dphi[0][1]);  /* omega_xy */
    double w_02 = 0.5 * (dphi[2][0] - dphi[0][2]);  /* omega_xz */
    double w_12 = 0.5 * (dphi[2][1] - dphi[1][2]);  /* omega_yz */

    /* Omega_x = eps_{x,yz} omega_{yz} - eps_{x,zy} omega_{zy} = 2*omega_{yz} */
    /* Actually Omega_k = eps_{ijk} omega_{ij} = omega_{12} - omega_{21} + ... */
    /* Omega_z = omega_{xy} - omega_{yx} = 2*omega_{xy} since omega is antisymmetric */
    /* But standard: Omega_k = eps_{ijk} (d_i phi_j - d_j phi_i)/2  */
    /* Omega_z = (d_x phi_y - d_y phi_x) = 2*w_01 */
    omega_out[0] = 2.0 * w_12;  /* Omega_x */
    omega_out[1] = -2.0 * w_02; /* Omega_y = eps_{yxz}*w_xz = -(d_x phi_z - d_z phi_x) */
    omega_out[2] = 2.0 * w_01;  /* Omega_z */

    /* Strain */
    eps_out[0] = dphi[0][0];                              /* eps_xx */
    eps_out[1] = dphi[1][1];                              /* eps_yy */
    eps_out[2] = dphi[2][2];                              /* eps_zz */
    eps_out[3] = 0.5 * (dphi[1][0] + dphi[0][1]);        /* eps_xy */
    eps_out[4] = 0.5 * (dphi[2][0] + dphi[0][2]);        /* eps_xz */
    eps_out[5] = 0.5 * (dphi[2][1] + dphi[1][2]);        /* eps_yz */
}

/* --- Compute torsion flux Phi_T = integral of Omega_z over xy-plane at z=0 --- */
static double compute_torsion_flux(void)
{
    int k_mid = N / 2;
    double flux = 0.0;

    #pragma omp parallel for reduction(+:flux) schedule(static)
    for (int i = 2; i < N-2; i++) {
        for (int j = 2; j < N-2; j++) {
            double omega[3], eps[6];
            compute_torsion_strain(i, j, k_mid, omega, eps);
            flux += omega[2] * dx * dx;  /* Omega_z * dA */
        }
    }
    return flux;
}

/* --- Fast center density --- */
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

/* --- Compute full diagnostics --- */
static Diag compute_diag(double core_radius)
{
    Diag d;
    memset(&d, 0, sizeof(d));

    /* Center values */
    int ci = N/2, cj = N/2, ck = N/2;
    double omega_c[3], eps_c[6];
    compute_torsion_strain(ci, cj, ck, omega_c, eps_c);
    d.omega_xy = omega_c[2] / 2.0;  /* back to omega_{xy} from Omega_z */
    d.omega_xz = -omega_c[1] / 2.0;
    d.omega_yz = omega_c[0] / 2.0;
    d.eps_xx = eps_c[0]; d.eps_yy = eps_c[1]; d.eps_zz = eps_c[2];
    d.eps_xy = eps_c[3]; d.eps_xz = eps_c[4]; d.eps_yz = eps_c[5];

    d.rho_center = center_rho();
    d.torsion_flux_z = compute_torsion_flux();

    double Ek = 0, Eg = 0, Em = 0, Ep = 0;
    double Ecore = 0, Eall = 0;
    double peak[3] = {0,0,0};
    double peak_P = 0;

    #pragma omp parallel
    {
        double lEk = 0, lEg = 0, lEm = 0, lEp = 0, lEc = 0, lEa = 0;
        double lpk[3] = {0,0,0};
        double lpkP = 0;

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
        }
    }

    d.Ek = Ek; d.Eg = Eg; d.Ep = Ep;
    d.Et = Ek + Eg + Em + Ep;
    for (int a = 0; a < 3; a++) d.peak[a] = peak[a];
    d.fc = (Eall > 1e-20) ? Ecore / Eall : 0.0;
    d.peak_P = peak_P;
    return d;
}

/* --- DFT peak frequency finder --- */
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

/* --- Strain multipoles on spherical shell --- */
static void compute_strain_multipoles(double R_shell, FILE *fstrain, FILE *fmulti)
{
    int N_sample = 200;
    double golden_ratio = (1.0 + sqrt(5.0)) / 2.0;

    double c_l[3] = {0, 0, 0};
    double norm = 0;
    int n_valid = 0;

    fprintf(fstrain, "# Strain on shell r=%.1f\n", R_shell);
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

        double omega_pt[3], eps[6];
        compute_torsion_strain(ig, jg, kg, omega_pt, eps);

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

    fprintf(fmulti, "# Multipole decomposition of |sigma|^2 on shell r=%.1f\n", R_shell);
    fprintf(fmulti, "l\tcoeff\tfraction\n");
    double sum_abs = fabs(c_l[0]) + fabs(c_l[1]) + fabs(c_l[2]);
    for (int l = 0; l <= 2; l++) {
        double frac = (sum_abs > 1e-30) ? fabs(c_l[l]) / sum_abs : 0.0;
        fprintf(fmulti, "%d\t%.6e\t%.6f\n", l, c_l[l], frac);
        printf("    l=%d: coefficient = %.6e, fraction = %.4f\n", l, c_l[l], frac);
    }
    printf("    Total |sigma|^2 on shell: %.6e (from %d valid points)\n",
           (n_valid > 0) ? norm / n_valid : 0.0, n_valid);
}

/* --- Run single (kappa_S, kappa_T) configuration --- */
static void run_single(const char *label, double kS, double kT,
                       int use_metric, FILE *fsummary)
{
    kappa_S = kS;
    kappa_T = kT;

    printf("\n================================================================\n");
    printf("  Run: %s  (kappa_S=%.2f, kappa_T=%.2f, metric=%d)\n",
           label, kS, kT, use_metric);
    printf("================================================================\n\n");

    /* Zero fields */
    for (int a = 0; a < 3; a++) {
        memset(phi[a], 0, Ngrid * sizeof(double));
        memset(vel[a], 0, Ngrid * sizeof(double));
        memset(acc_arr[a], 0, Ngrid * sizeof(double));
    }

    init_twisted_tube();
    set_damping();

    /* Adjust dt for effective wave speed: c_eff = sqrt(1 + kS/2 + kT/2) */
    double c_eff = sqrt(1.0 + 0.5*kS + 0.5*kT);
    dt = cfl_frac * dx / c_eff;
    printf("  c_eff=%.3f  dt=%.5f\n", c_eff, dt);

    double core_radius = 8.0;
    int Nt = (int)(tfinal / dt) + 1;

    /* DFT storage */
    int max_dft = 50000;
    double *rho_hist = malloc(max_dft * sizeof(double));
    double *t_hist   = malloc(max_dft * sizeof(double));
    int n_dft = 0;
    int dft_every = Nt / max_dft;
    if (dft_every < 1) dft_every = 1;

    int print_every = Nt / 20;
    if (print_every < 1) print_every = 1;
    int diag_every = Nt / 2000;
    if (diag_every < 1) diag_every = 1;

    /* Output file */
    char path[600];
    snprintf(path, sizeof(path), "%s/v26p4_%s.tsv", outdir, label);
    FILE *fout = fopen(path, "w");
    if (!fout) { fprintf(stderr, "Cannot open %s\n", path); free(rho_hist); free(t_hist); return; }
    fprintf(fout, "time\tE_total\tE_kin\tE_grad\tE_pot\tfc\tpeak0\tpeak1\tpeak2\t"
                  "peak_P\trho_ctr\t"
                  "omega_xy\tomega_xz\tomega_yz\t"
                  "eps_xx\teps_yy\teps_zz\teps_xy\teps_xz\teps_yz\t"
                  "torsion_flux_z\n");

    if (use_metric)
        compute_acc_metric(0.0);
    else
        compute_acc();

    double wall_start = omp_get_wtime();

    /* Initial */
    Diag d0 = compute_diag(core_radius);
    printf("  t=%7.1f  E=%.2f  fc=%.4f  pk=(%.3f,%.3f,%.3f)  |P|=%.6f  Phi_T=%.4f\n",
           0.0, d0.Et, d0.fc, d0.peak[0], d0.peak[1], d0.peak[2], d0.peak_P, d0.torsion_flux_z);

    double fc_final = 0, peak_P_final = 0, E_final = 0;
    double omega_xy_final = 0, eps_xx_final = 0;
    double torsion_flux_final = 0;

    for (int n = 0; n <= Nt; n++) {
        double t = n * dt;

        /* DFT history */
        if (n % dft_every == 0 && n_dft < max_dft) {
            rho_hist[n_dft] = center_rho();
            t_hist[n_dft] = t;
            n_dft++;
        }

        int do_diag  = (n % diag_every == 0);
        int do_print = (n % print_every == 0);

        if (do_diag || do_print) {
            Diag d = compute_diag(core_radius);

            if (do_diag) {
                fprintf(fout, "%.4f\t%.6e\t%.6e\t%.6e\t%.6e\t%.6f\t%.6e\t%.6e\t%.6e\t"
                              "%.6e\t%.6e\t"
                              "%.6e\t%.6e\t%.6e\t"
                              "%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t"
                              "%.6e\n",
                        t, d.Et, d.Ek, d.Eg, d.Ep, d.fc,
                        d.peak[0], d.peak[1], d.peak[2], d.peak_P, d.rho_center,
                        d.omega_xy, d.omega_xz, d.omega_yz,
                        d.eps_xx, d.eps_yy, d.eps_zz, d.eps_xy, d.eps_xz, d.eps_yz,
                        d.torsion_flux_z);
            }

            if (do_print) {
                double elapsed = omp_get_wtime() - wall_start;
                double frac = (double)n / Nt;
                double eta_t = (frac > 0.001) ? elapsed * (1.0-frac)/frac : 0;
                printf("  t=%7.1f  E=%.2f  fc=%.4f  pk=(%.3f,%.3f,%.3f)  |P|=%.6f  Phi_T=%.4f  [%.0fs, ETA %.0fs]\n",
                       t, d.Et, d.fc, d.peak[0], d.peak[1], d.peak[2],
                       d.peak_P, d.torsion_flux_z, elapsed, eta_t);
                fflush(stdout);
            }

            /* Store final values */
            fc_final = d.fc;
            peak_P_final = d.peak_P;
            E_final = d.Et;
            omega_xy_final = d.omega_xy;
            eps_xx_final = d.eps_xx;
            torsion_flux_final = d.torsion_flux_z;
        }

        if (n == Nt) break;

        if (use_metric) {
            double alpha_g = alpha_g_max;
            if (t < alpha_g_ramp_time)
                alpha_g = alpha_g_max * t / alpha_g_ramp_time;
            verlet_step_metric(alpha_g);
        } else {
            verlet_step_standard();
        }
    }

    fclose(fout);
    double elapsed = omp_get_wtime() - wall_start;

    /* Phase 2: breathing analysis */
    printf("\n  Breathing analysis:\n");
    double peak_power_rho = 0;
    double omega_rho = 0;
    double mean_rho = 0, var_rho = 0;
    if (n_dft >= 100) {
        int start = n_dft / 2;
        omega_rho = find_peak_omega(rho_hist, t_hist, n_dft, start, 5.0, &peak_power_rho);

        for (int j = start; j < n_dft; j++) mean_rho += rho_hist[j];
        mean_rho /= (n_dft - start);

        for (int j = start; j < n_dft; j++) {
            double dr = rho_hist[j] - mean_rho;
            var_rho += dr * dr;
        }
        var_rho /= (n_dft - start);

        double rel_var = (fabs(mean_rho) > 1e-30) ? var_rho / (mean_rho * mean_rho) : 0.0;
        printf("    rho(center): mean=%.4e, rel_var=%.4e, peak_omega=%.4f\n",
               mean_rho, rel_var, omega_rho);
    }

    /* Phase 3: strain multipoles */
    printf("\n  Strain multipoles (R=8):\n");
    snprintf(path, sizeof(path), "%s/v26p4_%s_strain.tsv", outdir, label);
    FILE *fstrain = fopen(path, "w");
    snprintf(path, sizeof(path), "%s/v26p4_%s_multipoles.tsv", outdir, label);
    FILE *fmulti = fopen(path, "w");
    if (fstrain && fmulti) {
        compute_strain_multipoles(8.0, fstrain, fmulti);
        fclose(fstrain);
        fclose(fmulti);
    }

    /* Torsion flux at multiple z-planes */
    printf("\n  Torsion flux Phi_T(z=0) = %.6f  (%.4f * 2*pi)\n",
           torsion_flux_final, torsion_flux_final / (2.0 * M_PI));

    /* Summary line */
    int survived = (fc_final > 0.01 && peak_P_final > 1e-6);
    printf("\n  RESULT: survived=%s  fc=%.4f  |P|max=%.6f  E=%.2f  Phi_T=%.4f  (%.1fs)\n",
           survived ? "YES" : "NO", fc_final, peak_P_final, E_final,
           torsion_flux_final, elapsed);

    if (fsummary) {
        fprintf(fsummary, "%.2f\t%.2f\t%.2f\t%d\t%.4f\t%.6f\t%.4f\t%.6e\t%.6e\t%.6e\t%.4f\t%.4f\n",
                mass, kS, kT, survived, fc_final, peak_P_final, E_final,
                omega_xy_final, eps_xx_final, torsion_flux_final,
                omega_rho, (fabs(mean_rho) > 1e-30) ? var_rho / (mean_rho * mean_rho) : 0.0);
    }

    free(rho_hist);
    free(t_hist);
}

/* --- Phase 4c: Detailed torsion flux analysis --- */
static void run_phase4c(double kS_best, double kT_best, FILE *fsummary)
{
    printf("\n================================================================\n");
    printf("  Phase 4c: Torsion Flux Quantization\n");
    printf("  (kappa_S=%.2f, kappa_T=%.2f)\n", kS_best, kT_best);
    printf("================================================================\n\n");

    kappa_S = kS_best;
    kappa_T = kT_best;

    /* Zero and init */
    for (int a = 0; a < 3; a++) {
        memset(phi[a], 0, Ngrid * sizeof(double));
        memset(vel[a], 0, Ngrid * sizeof(double));
        memset(acc_arr[a], 0, Ngrid * sizeof(double));
    }
    init_twisted_tube();
    set_damping();

    /* Adjust dt for effective wave speed */
    double c_eff = sqrt(1.0 + 0.5*kS_best + 0.5*kT_best);
    dt = cfl_frac * dx / c_eff;

    /* Evolve to equilibrium */
    compute_acc();
    int Nt = (int)(tfinal / dt) + 1;
    int print_every = Nt / 10;
    if (print_every < 1) print_every = 1;

    printf("  Evolving to t=%.0f...\n", tfinal);
    for (int n = 0; n < Nt; n++) {
        verlet_step_standard();
        if (n % print_every == 0) {
            double t = n * dt;
            printf("    t=%.0f\n", t);
            fflush(stdout);
        }
    }

    /* Compute torsion flux at multiple z-planes */
    printf("\n  Torsion flux at different z-planes:\n");
    char path[600];
    snprintf(path, sizeof(path), "%s/v26p4_phase4c_flux.tsv", outdir);
    FILE *fflux = fopen(path, "w");
    if (fflux) {
        fprintf(fflux, "z\tOmega_z_flux\tOmega_x_flux\tOmega_y_flux\n");

        for (int kslice = 4; kslice < N-4; kslice += 4) {
            double z = -L + kslice * dx;
            double flux_z = 0, flux_x = 0, flux_y = 0;

            for (int i = 2; i < N-2; i++) {
                for (int j = 2; j < N-2; j++) {
                    double omega[3], eps[6];
                    compute_torsion_strain(i, j, kslice, omega, eps);
                    flux_z += omega[2] * dx * dx;
                    flux_x += omega[0] * dx * dx;
                    flux_y += omega[1] * dx * dx;
                }
            }

            fprintf(fflux, "%.4f\t%.6e\t%.6e\t%.6e\n", z, flux_z, flux_x, flux_y);
            if (kslice % 16 == 0) {
                printf("    z=%6.2f: Phi_z=%.4f (%.3f*2pi)  Phi_x=%.4f  Phi_y=%.4f\n",
                       z, flux_z, flux_z/(2*M_PI), flux_x, flux_y);
            }
        }
        fclose(fflux);
    }

    /* Also compute Omega field in a cross-section */
    snprintf(path, sizeof(path), "%s/v26p4_phase4c_vorticity_xy.tsv", outdir);
    FILE *fvor = fopen(path, "w");
    if (fvor) {
        int k_mid = N/2;
        fprintf(fvor, "# Vorticity cross-section at z=0 (k=%d)\n", k_mid);
        fprintf(fvor, "x\ty\tOmega_x\tOmega_y\tOmega_z\t|Omega|\n");
        for (int i = 2; i < N-2; i += 2) {
            double x = -L + i * dx;
            for (int j = 2; j < N-2; j += 2) {
                double y = -L + j * dx;
                double omega[3], eps[6];
                compute_torsion_strain(i, j, k_mid, omega, eps);
                double omag = sqrt(omega[0]*omega[0] + omega[1]*omega[1] + omega[2]*omega[2]);
                fprintf(fvor, "%.4f\t%.4f\t%.6e\t%.6e\t%.6e\t%.6e\n",
                        x, y, omega[0], omega[1], omega[2], omag);
            }
        }
        fclose(fvor);
    }

    if (fsummary) {
        fprintf(fsummary, "\nPhase 4c: Torsion flux quantization at (kS=%.2f, kT=%.2f)\n",
                kS_best, kT_best);
        fprintf(fsummary, "See v26p4_phase4c_flux.tsv for z-dependent flux\n");
    }
}

/* --- Main --- */
int main(int argc, char **argv)
{
    parse_args(argc, argv);

    dx  = 2.0 * L / (N - 1);
    dx2 = dx * dx;
    m2  = mass * mass;
    Ngrid = (long)N * N * N;
    /* dt is set per-run to account for kappa_S, kappa_T effective speed */

    printf("=== V26 Phase 4: Torsion Constraint and Teleparallel Gravity ===\n");
    printf("Parameters:\n");
    printf("  mu=%.1f  kappa=%.1f  mass=%.3f  A0=%.3f  R_tube=%.1f\n",
           mu_pot, kappa_pot, mass, A0, R_tube);
    printf("  N=%d  L=%.1f  dx=%.4f  dt=%.5f\n", N, L, dx, dt);
    printf("  Ngrid=%ld (%.1f M)  Memory: %.1f MB\n",
           Ngrid, Ngrid/1e6, Ngrid*8.0*10/1e6);
    printf("  Threads: %d\n", omp_get_max_threads());
    printf("  phase=%d  tfinal=%.0f\n", phase, tfinal);
    fflush(stdout);

    /* Allocate */
    for (int a = 0; a < 3; a++) {
        phi[a] = calloc(Ngrid, sizeof(double));
        vel[a] = calloc(Ngrid, sizeof(double));
        acc_arr[a] = calloc(Ngrid, sizeof(double));
        if (!phi[a] || !vel[a] || !acc_arr[a]) {
            fprintf(stderr, "Allocation failed\n");
            return 1;
        }
    }
    damp = malloc(Ngrid * sizeof(double));
    if (!damp) { fprintf(stderr, "Allocation failed for damp\n"); return 1; }

    /* Summary file */
    char path[600];
    snprintf(path, sizeof(path), "%s/v26p4_summary.tsv", outdir);
    FILE *fsummary = fopen(path, "w");
    if (fsummary)
        fprintf(fsummary, "mass\tkappa_S\tkappa_T\tsurvived\tfc\tpeak_P\tE_total\t"
                          "omega_xy_ctr\teps_xx_ctr\ttorsion_flux_z\t"
                          "breathing_omega\tbreathing_relvar\n");

    if (phase == 0 || phase == 99) {
        /* ============================================
         * Phase 4a: Torsion scan (kappa_T scan, kappa_S=0)
         * Keep mass=1 (v26 mode 2 baseline).
         * KEY Q: does kappa_T TIGHTEN the braid (increase fc)?
         * Use tfinal=200 for scan.
         * ============================================ */
        printf("\n\n########################################\n");
        printf("# PHASE 4a: Torsion Scan (m=1)          #\n");
        printf("########################################\n");

        double save_tf = tfinal;
        tfinal = 200.0;
        mass = 1.0; m2 = 1.0;

        double kT_vals[] = {0.0, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0};
        int n_kT = 7;
        for (int ik = 0; ik < n_kT; ik++) {
            char label[64];
            snprintf(label, sizeof(label), "4a_kT%.2f", kT_vals[ik]);
            run_single(label, 0.0, kT_vals[ik], 0, fsummary);
        }

        tfinal = save_tf;
    }

    if (phase == 1 || phase == 99) {
        /* ============================================
         * Phase 4b: Joint (kappa_S, kappa_T) scan, m=1
         * Does kappa_S != kappa_T improve l=2?
         * Use tfinal=200 for scan.
         * ============================================ */
        printf("\n\n########################################\n");
        printf("# PHASE 4b: Strain + Torsion Scan (m=1) #\n");
        printf("########################################\n");

        double save_tf = tfinal;
        tfinal = 200.0;
        mass = 1.0; m2 = 1.0;
        double pairs[][2] = {
            {0.5, 0.0}, {0.0, 0.5}, {0.5, 0.5},
            {1.0, 0.5}, {0.5, 1.0}, {0.5, 2.0}
        };
        int n_pairs = 6;
        for (int ip = 0; ip < n_pairs; ip++) {
            char label[64];
            snprintf(label, sizeof(label), "4b_kS%.2f_kT%.2f", pairs[ip][0], pairs[ip][1]);
            run_single(label, pairs[ip][0], pairs[ip][1], 0, fsummary);
        }
        tfinal = save_tf;
    }

    if (phase == 2 || phase == 99) {
        /* ============================================
         * Phase 4c: Torsion flux quantization
         * Use best kT from 4a (we'll use kT=2.0 default)
         * Full tfinal=500 for accurate flux measurement.
         * ============================================ */
        printf("\n\n########################################\n");
        printf("# PHASE 4c: Torsion Flux Quantization   #\n");
        printf("########################################\n");

        mass = 1.0; m2 = 1.0;
        run_phase4c(0.0, 2.0, fsummary);
    }

    if (phase == 3 || phase == 99) {
        /* ============================================
         * Phase 4d: Self-consistent metric
         * Full tfinal=500.
         * ============================================ */
        printf("\n\n########################################\n");
        printf("# PHASE 4d: Self-Consistent Metric      #\n");
        printf("########################################\n");

        mass = 1.0; m2 = 1.0;
        char label[64];
        snprintf(label, sizeof(label), "4d_metric_kS0.00_kT2.00");
        run_single(label, 0.0, 2.0, 1, fsummary);
    }

    if (fsummary) fclose(fsummary);

    /* Cleanup */
    for (int a = 0; a < 3; a++) {
        free(phi[a]); free(vel[a]); free(acc_arr[a]);
    }
    free(damp);

    printf("\n=== V26 Phase 4 Complete ===\n");
    return 0;
}
