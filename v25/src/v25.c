/*
 * v25.c — Three-Field Elastic Gravity in 3D
 *
 * Phases 1-4:
 *   1. 3D elastic oscillon with cross-gradient, Lame, pairwise couplings
 *   2. Strain field diagnostics: compression/shear decomposition, multipoles
 *   3. Self-consistent metric: h_{ij} = d_i phi_j + d_j phi_i modifies EOM
 *   4. Gravitational wave detection: TT strain at 4 angles after velocity kick
 *
 * Lagrangian:
 *   L = Sum_a [1/2(dt phi_a)^2 - 1/2|grad phi_a|^2 - 1/2 m^2 phi_a^2]
 *     - V(P)                               P = phi_1 phi_2 phi_3
 *     - lambda_pw(phi_1 phi_2 + phi_2 phi_3 + phi_3 phi_1)
 *     - 1/2 eta (d_i phi_j)(d_j phi_i)     cross-gradient
 *     - 1/2 lambda_L (div phi)^2            Lame compression
 *
 * Self-consistent metric (Phase 3):
 *   Laplacian: d^2 phi_a -> (delta_ij + alpha_g * h_ij) d_i d_j phi_a
 *   Mass:      m^2 -> m^2 (1 + alpha_g * theta)
 *   h_ij = d_i phi_j + d_j phi_i,  theta = h_kk / 2
 *
 * Compile: gcc -O3 -fopenmp -Wall -o v25 src/v25.c -lm
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
static double A_init    = 0.8;
static double sig_init  = 3.0;
static double lambda_pw = 0.5;
static double eta       = 0.1;
static double lambda_L  = 0.1;
static double alpha_g   = 0.001;  /* metric coupling */
static double v_kick    = 0.01;   /* Phase 4 velocity kick */

static int    N         = 96;
static double L         = 15.0;
static double tfinal    = 500.0;
static double cfl_frac  = 0.20;
static char   outdir[512] = "data";

static void parse_args(int argc, char **argv)
{
    for (int i = 1; i < argc - 1; i += 2) {
        if      (!strcmp(argv[i], "-mu"))       mu_pot    = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-kappa"))    kappa     = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-mass"))     mass      = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-A"))        A_init    = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-sigma"))    sig_init  = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-N"))        N         = atoi(argv[i+1]);
        else if (!strcmp(argv[i], "-L"))        L         = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-tfinal"))   tfinal    = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-cfl"))      cfl_frac  = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-lpw"))      lambda_pw = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-eta"))      eta       = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-lL"))       lambda_L  = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-ag"))       alpha_g   = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-vkick"))    v_kick    = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-o"))        strncpy(outdir, argv[i+1], sizeof(outdir)-1);
        else { fprintf(stderr, "Unknown flag: %s\n", argv[i]); exit(1); }
    }
}

/* ─── Index helpers ─── */
#define IDX(i,j,k) ((long)(i)*N*N + (long)(j)*N + (long)(k))

/* ─── Globals ─── */
static double *phi[3], *vel[3], *acc[3];
static double *damp;
static double dx, dx2, m2, dt;
static long Ngrid;

/* ─── Hermite smooth step for ramp-up ─── */
static double hermite_step(double t, double t0, double t1)
{
    if (t <= t0) return 0.0;
    if (t >= t1) return 1.0;
    double s = (t - t0) / (t1 - t0);
    return s * s * (3.0 - 2.0 * s);
}

/* ─── Compute acceleration for Phase 1 (elastic, no metric) ─── */
static void compute_acc_elastic(void)
{
    double eta_lL = eta + lambda_L;  /* combined cross-gradient + Lame coefficient */

    for (int a = 0; a < 3; a++) {
        #pragma omp parallel for schedule(static)
        for (int i = 2; i < N-2; i++) {
            for (int j = 2; j < N-2; j++) {
                for (int k = 2; k < N-2; k++) {
                    long idx = IDX(i,j,k);

                    /* 7-point Laplacian */
                    double lapl = (phi[a][IDX(i+1,j,k)] + phi[a][IDX(i-1,j,k)]
                                 + phi[a][IDX(i,j+1,k)] + phi[a][IDX(i,j-1,k)]
                                 + phi[a][IDX(i,j,k+1)] + phi[a][IDX(i,j,k-1)]
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

                    /* Pairwise coupling: -lambda_pw * (phi_b + phi_c) */
                    double pw_force = 0.0;
                    for (int b = 0; b < 3; b++) {
                        if (b != a) pw_force += phi[b][idx];
                    }
                    pw_force *= -lambda_pw;

                    /* Cross-gradient + Lame: (eta + lambda_L) * d_a(div_phi)
                     * div_phi = d_x phi_0 + d_y phi_1 + d_z phi_2
                     * d_a(div_phi) = d_a d_x phi_0 + d_a d_y phi_1 + d_a d_z phi_2
                     *
                     * For diagonal term (a=b): d_a^2 phi_a (standard second deriv)
                     * For off-diagonal (a!=b): mixed partial d_a d_b phi_b
                     *   using 4-point formula: [f(a+1,b+1)-f(a+1,b-1)-f(a-1,b+1)+f(a-1,b-1)]/(4 dx^2)
                     */
                    double d_a_div = 0.0;

                    /* b=0: d_a(d_x phi_0) */
                    if (a == 0) {
                        /* d_x^2 phi_0 */
                        d_a_div += (phi[0][IDX(i+1,j,k)] - 2.0*phi[0][idx]
                                  + phi[0][IDX(i-1,j,k)]) / dx2;
                    } else if (a == 1) {
                        /* d_y d_x phi_0: mixed */
                        d_a_div += (phi[0][IDX(i+1,j+1,k)] - phi[0][IDX(i+1,j-1,k)]
                                  - phi[0][IDX(i-1,j+1,k)] + phi[0][IDX(i-1,j-1,k)])
                                  / (4.0 * dx2);
                    } else {
                        /* d_z d_x phi_0: mixed */
                        d_a_div += (phi[0][IDX(i+1,j,k+1)] - phi[0][IDX(i+1,j,k-1)]
                                  - phi[0][IDX(i-1,j,k+1)] + phi[0][IDX(i-1,j,k-1)])
                                  / (4.0 * dx2);
                    }

                    /* b=1: d_a(d_y phi_1) */
                    if (a == 1) {
                        /* d_y^2 phi_1 */
                        d_a_div += (phi[1][IDX(i,j+1,k)] - 2.0*phi[1][idx]
                                  + phi[1][IDX(i,j-1,k)]) / dx2;
                    } else if (a == 0) {
                        /* d_x d_y phi_1: mixed */
                        d_a_div += (phi[1][IDX(i+1,j+1,k)] - phi[1][IDX(i+1,j-1,k)]
                                  - phi[1][IDX(i-1,j+1,k)] + phi[1][IDX(i-1,j-1,k)])
                                  / (4.0 * dx2);
                    } else {
                        /* d_z d_y phi_1: mixed */
                        d_a_div += (phi[1][IDX(i,j+1,k+1)] - phi[1][IDX(i,j+1,k-1)]
                                  - phi[1][IDX(i,j-1,k+1)] + phi[1][IDX(i,j-1,k-1)])
                                  / (4.0 * dx2);
                    }

                    /* b=2: d_a(d_z phi_2) */
                    if (a == 2) {
                        /* d_z^2 phi_2 */
                        d_a_div += (phi[2][IDX(i,j,k+1)] - 2.0*phi[2][idx]
                                  + phi[2][IDX(i,j,k-1)]) / dx2;
                    } else if (a == 0) {
                        /* d_x d_z phi_2: mixed */
                        d_a_div += (phi[2][IDX(i+1,j,k+1)] - phi[2][IDX(i+1,j,k-1)]
                                  - phi[2][IDX(i-1,j,k+1)] + phi[2][IDX(i-1,j,k-1)])
                                  / (4.0 * dx2);
                    } else {
                        /* d_y d_z phi_2: mixed */
                        d_a_div += (phi[2][IDX(i,j+1,k+1)] - phi[2][IDX(i,j+1,k-1)]
                                  - phi[2][IDX(i,j-1,k+1)] + phi[2][IDX(i,j-1,k-1)])
                                  / (4.0 * dx2);
                    }

                    acc[a][idx] = lapl - m2 * phi[a][idx] - dVdphi
                                + pw_force + eta_lL * d_a_div;
                }
            }
        }

        /* Boundary: zero acceleration (layers 0,1 and N-2,N-1) */
        #pragma omp parallel for schedule(static)
        for (long idx = 0; idx < Ngrid; idx++) {
            int ii = idx / (N * N);
            int jj = (idx / N) % N;
            int kk = idx % N;
            if (ii < 2 || ii >= N-2 || jj < 2 || jj >= N-2 || kk < 2 || kk >= N-2)
                acc[a][idx] = 0.0;
        }
    }
}

/* ─── Compute acceleration for Phase 3 (with self-consistent metric) ─── */
static void compute_acc_metric(double ag_eff)
{
    double eta_lL = eta + lambda_L;

    for (int a = 0; a < 3; a++) {
        #pragma omp parallel for schedule(static)
        for (int i = 2; i < N-2; i++) {
            for (int j = 2; j < N-2; j++) {
                for (int k = 2; k < N-2; k++) {
                    long idx = IDX(i,j,k);

                    /* Standard 7-point Laplacian */
                    double lapl = (phi[a][IDX(i+1,j,k)] + phi[a][IDX(i-1,j,k)]
                                 + phi[a][IDX(i,j+1,k)] + phi[a][IDX(i,j-1,k)]
                                 + phi[a][IDX(i,j,k+1)] + phi[a][IDX(i,j,k-1)]
                                 - 6.0 * phi[a][idx]) / dx2;

                    /* Metric correction: h_{ij} d_i d_j phi_a
                     * h_{ij} = d_i phi_j + d_j phi_i
                     * Diagonal only approximation (as specified):
                     *   h_{11} = 2 * d_x phi_0
                     *   h_{22} = 2 * d_y phi_1
                     *   h_{33} = 2 * d_z phi_2
                     * theta = (h_{11}+h_{22}+h_{33})/2 = d_x phi_0 + d_y phi_1 + d_z phi_2
                     */
                    double dx_phi0 = (phi[0][IDX(i+1,j,k)] - phi[0][IDX(i-1,j,k)]) / (2.0*dx);
                    double dy_phi1 = (phi[1][IDX(i,j+1,k)] - phi[1][IDX(i,j-1,k)]) / (2.0*dx);
                    double dz_phi2 = (phi[2][IDX(i,j,k+1)] - phi[2][IDX(i,j,k-1)]) / (2.0*dx);

                    double h11 = 2.0 * dx_phi0;
                    double h22 = 2.0 * dy_phi1;
                    double h33 = 2.0 * dz_phi2;
                    double theta = dx_phi0 + dy_phi1 + dz_phi2;  /* = tr(h)/2 */

                    /* Diagonal second derivatives of phi_a */
                    double d2x = (phi[a][IDX(i+1,j,k)] - 2.0*phi[a][idx]
                                + phi[a][IDX(i-1,j,k)]) / dx2;
                    double d2y = (phi[a][IDX(i,j+1,k)] - 2.0*phi[a][idx]
                                + phi[a][IDX(i,j-1,k)]) / dx2;
                    double d2z = (phi[a][IDX(i,j,k+1)] - 2.0*phi[a][idx]
                                + phi[a][IDX(i,j,k-1)]) / dx2;

                    /* Metric-corrected Laplacian: (delta_ij + ag * h_ij) d_i d_j phi_a */
                    double metric_corr = ag_eff * (h11 * d2x + h22 * d2y + h33 * d2z);

                    /* Mass correction: m^2 -> m^2(1 + ag * theta) */
                    double m2_eff = m2 * (1.0 + ag_eff * theta);

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

                    /* Pairwise coupling */
                    double pw_force = 0.0;
                    for (int b = 0; b < 3; b++) {
                        if (b != a) pw_force += phi[b][idx];
                    }
                    pw_force *= -lambda_pw;

                    /* Cross-gradient + Lame (same as elastic) */
                    double d_a_div = 0.0;

                    if (a == 0) {
                        d_a_div += (phi[0][IDX(i+1,j,k)] - 2.0*phi[0][idx]
                                  + phi[0][IDX(i-1,j,k)]) / dx2;
                    } else if (a == 1) {
                        d_a_div += (phi[0][IDX(i+1,j+1,k)] - phi[0][IDX(i+1,j-1,k)]
                                  - phi[0][IDX(i-1,j+1,k)] + phi[0][IDX(i-1,j-1,k)])
                                  / (4.0 * dx2);
                    } else {
                        d_a_div += (phi[0][IDX(i+1,j,k+1)] - phi[0][IDX(i+1,j,k-1)]
                                  - phi[0][IDX(i-1,j,k+1)] + phi[0][IDX(i-1,j,k-1)])
                                  / (4.0 * dx2);
                    }

                    if (a == 1) {
                        d_a_div += (phi[1][IDX(i,j+1,k)] - 2.0*phi[1][idx]
                                  + phi[1][IDX(i,j-1,k)]) / dx2;
                    } else if (a == 0) {
                        d_a_div += (phi[1][IDX(i+1,j+1,k)] - phi[1][IDX(i+1,j-1,k)]
                                  - phi[1][IDX(i-1,j+1,k)] + phi[1][IDX(i-1,j-1,k)])
                                  / (4.0 * dx2);
                    } else {
                        d_a_div += (phi[1][IDX(i,j+1,k+1)] - phi[1][IDX(i,j+1,k-1)]
                                  - phi[1][IDX(i,j-1,k+1)] + phi[1][IDX(i,j-1,k-1)])
                                  / (4.0 * dx2);
                    }

                    if (a == 2) {
                        d_a_div += (phi[2][IDX(i,j,k+1)] - 2.0*phi[2][idx]
                                  + phi[2][IDX(i,j,k-1)]) / dx2;
                    } else if (a == 0) {
                        d_a_div += (phi[2][IDX(i+1,j,k+1)] - phi[2][IDX(i+1,j,k-1)]
                                  - phi[2][IDX(i-1,j,k+1)] + phi[2][IDX(i-1,j,k-1)])
                                  / (4.0 * dx2);
                    } else {
                        d_a_div += (phi[2][IDX(i,j+1,k+1)] - phi[2][IDX(i,j+1,k-1)]
                                  - phi[2][IDX(i,j-1,k+1)] + phi[2][IDX(i,j-1,k-1)])
                                  / (4.0 * dx2);
                    }

                    acc[a][idx] = lapl + metric_corr - m2_eff * phi[a][idx]
                                - dVdphi + pw_force + eta_lL * d_a_div;
                }
            }
        }

        /* Boundary: zero acceleration */
        #pragma omp parallel for schedule(static)
        for (long idx = 0; idx < Ngrid; idx++) {
            int ii = idx / (N * N);
            int jj = (idx / N) % N;
            int kk = idx % N;
            if (ii < 2 || ii >= N-2 || jj < 2 || jj >= N-2 || kk < 2 || kk >= N-2)
                acc[a][idx] = 0.0;
        }
    }
}

/* ─── Compute energy diagnostics ─── */
typedef struct {
    double Ek, Eg, Em, Ep, Epw, Et;
    double peak[3];
    double fc;
    double phi_center[3];
} Diag;

static Diag compute_diag(double core_radius)
{
    Diag d;
    memset(&d, 0, sizeof(d));
    long ic = IDX(N/2, N/2, N/2);
    for (int a = 0; a < 3; a++)
        d.phi_center[a] = phi[a][ic];

    double Ek=0, Eg=0, Em=0, Ep=0, Epw=0;
    double Ecore=0, Eall=0;
    double peak[3] = {0,0,0};

    #pragma omp parallel
    {
        double lEk=0, lEg=0, lEm=0, lEp=0, lEpw=0, lEc=0, lEa=0;
        double lpk[3] = {0,0,0};

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

                    double p0 = phi[0][idx], p1 = phi[1][idx], p2 = phi[2][idx];
                    double P = p0 * p1 * p2;
                    double P2 = P * P;
                    double Vloc = 0.5 * mu_pot * P2 / (1.0 + kappa * P2);
                    lEp += Vloc * dV;
                    e_loc += Vloc;

                    double pw_e = lambda_pw * (p0*p1 + p1*p2 + p2*p0);
                    lEpw += pw_e * dV;
                    e_loc += pw_e;

                    lEa += e_loc * dV;
                    double r = sqrt(x*x + y*y + z*z);
                    if (r < core_radius) lEc += e_loc * dV;
                }
            }
        }

        #pragma omp critical
        {
            Ek += lEk; Eg += lEg; Em += lEm; Ep += lEp; Epw += lEpw;
            Ecore += lEc; Eall += lEa;
            for (int a = 0; a < 3; a++)
                if (lpk[a] > peak[a]) peak[a] = lpk[a];
        }
    }

    d.Ek = Ek; d.Eg = Eg; d.Em = Em; d.Ep = Ep; d.Epw = Epw;
    d.Et = Ek + Eg + Em + Ep + Epw;
    for (int a = 0; a < 3; a++) d.peak[a] = peak[a];
    d.fc = (Eall > 1e-20) ? Ecore / Eall : 0.0;

    return d;
}

/* ─── DFT to find peak frequency ─── */
static double find_peak_omega(double *hist, double *t_hist, int n_pts, int start,
                              double omega_max)
{
    double T = t_hist[n_pts-1] - t_hist[start];
    if (T < 10.0 || n_pts - start < 50) return 0.0;

    int nf = 500;
    double peak_pow = 0.0, peak_om = 0.0;
    for (int kk = 1; kk < nf; kk++) {
        double omega = omega_max * kk / nf;
        double re = 0, im = 0;
        for (int j = start; j < n_pts; j++) {
            double dtj = (j > start) ?
                (t_hist[j]-t_hist[j-1]) : (t_hist[start+1]-t_hist[start]);
            re += hist[j] * cos(omega * t_hist[j]) * dtj;
            im += hist[j] * sin(omega * t_hist[j]) * dtj;
        }
        double pw = re*re + im*im;
        if (pw > peak_pow) { peak_pow = pw; peak_om = omega; }
    }
    return peak_om;
}

/* ─── Compute strain tensor at a grid point ─── */
/* eps_{ij} = 0.5*(d_i phi_j + d_j phi_i) */
/* Returns 6 independent components: eps_xx, eps_yy, eps_zz, eps_xy, eps_xz, eps_yz */
static void compute_strain(int i, int j, int k, double eps[6])
{
    /* First derivatives d_i phi_a using centered differences */
    double dphi[3][3]; /* dphi[a][dir] = d_{dir} phi_a */
    for (int a = 0; a < 3; a++) {
        dphi[a][0] = (phi[a][IDX(i+1,j,k)] - phi[a][IDX(i-1,j,k)]) / (2.0*dx);
        dphi[a][1] = (phi[a][IDX(i,j+1,k)] - phi[a][IDX(i,j-1,k)]) / (2.0*dx);
        dphi[a][2] = (phi[a][IDX(i,j,k+1)] - phi[a][IDX(i,j,k-1)]) / (2.0*dx);
    }

    /* eps_{ij} = 0.5*(d_i phi_j + d_j phi_i)
     * field index = spatial index in elastic interpretation */
    eps[0] = dphi[0][0];                              /* eps_xx = d_x phi_0 */
    eps[1] = dphi[1][1];                              /* eps_yy = d_y phi_1 */
    eps[2] = dphi[2][2];                              /* eps_zz = d_z phi_2 */
    eps[3] = 0.5*(dphi[1][0] + dphi[0][1]);           /* eps_xy */
    eps[4] = 0.5*(dphi[2][0] + dphi[0][2]);           /* eps_xz */
    eps[5] = 0.5*(dphi[2][1] + dphi[1][2]);           /* eps_yz */
}

/* ─── Velocity Verlet step ─── */
static void verlet_step(int use_metric, double ag_eff)
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
    if (use_metric)
        compute_acc_metric(ag_eff);
    else
        compute_acc_elastic();

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

/* ─── Spherical harmonic helpers for multipole analysis ─── */
/* Y_{l,m} at (theta, phi) for real spherical harmonics, we use Legendre polynomials */
/* P_0(x) = 1, P_1(x) = x, P_2(x) = (3x^2-1)/2 */
static double legendre_P(int l, double x)
{
    switch (l) {
        case 0: return 1.0;
        case 1: return x;
        case 2: return 0.5*(3.0*x*x - 1.0);
        default: return 0.0;
    }
}

/* ─── Sample strain on spherical shell, compute multipoles ─── */
static void compute_strain_multipoles(double R_shell, FILE *fstrain, FILE *fmulti)
{
    /* Sample points on the sphere using Fibonacci spiral */
    int N_sample = 200;
    double golden_ratio = (1.0 + sqrt(5.0)) / 2.0;

    /* Accumulators for multipole decomposition of |sigma|^2 */
    /* We project onto P_l(cos theta) for l=0,1,2 */
    double c_l[3] = {0, 0, 0};   /* multipole coefficients */
    double norm = 0;

    fprintf(fstrain, "# Strain on shell r=%.1f\n", R_shell);
    fprintf(fstrain, "theta\tphi_ang\teps_xx\teps_yy\teps_zz\teps_xy\teps_xz\teps_yz\t"
                     "theta_tr\tsigma2\n");

    for (int ns = 0; ns < N_sample; ns++) {
        /* Fibonacci spiral point */
        double cos_th = 1.0 - 2.0 * (ns + 0.5) / N_sample;
        double sin_th = sqrt(1.0 - cos_th * cos_th);
        double phi_ang = 2.0 * M_PI * ns / golden_ratio;

        double x = R_shell * sin_th * cos(phi_ang);
        double y = R_shell * sin_th * sin(phi_ang);
        double z = R_shell * cos_th;

        /* Grid indices (nearest neighbor interpolation) */
        int i = (int)((x + L) / dx + 0.5);
        int j_idx = (int)((y + L) / dx + 0.5);
        int k = (int)((z + L) / dx + 0.5);

        if (i < 2 || i >= N-2 || j_idx < 2 || j_idx >= N-2 || k < 2 || k >= N-2)
            continue;

        double eps[6];
        compute_strain(i, j_idx, k, eps);

        /* Trace (compression = spin-0) */
        double theta_tr = eps[0] + eps[1] + eps[2];

        /* Traceless part sigma_{ij} = eps_{ij} - theta/3 * delta_{ij} */
        double sig[6];
        sig[0] = eps[0] - theta_tr / 3.0;
        sig[1] = eps[1] - theta_tr / 3.0;
        sig[2] = eps[2] - theta_tr / 3.0;
        sig[3] = eps[3];
        sig[4] = eps[4];
        sig[5] = eps[5];

        /* |sigma|^2 = sig_xx^2 + sig_yy^2 + sig_zz^2 + 2*(sig_xy^2+sig_xz^2+sig_yz^2) */
        double sigma2 = sig[0]*sig[0] + sig[1]*sig[1] + sig[2]*sig[2]
                       + 2.0*(sig[3]*sig[3] + sig[4]*sig[4] + sig[5]*sig[5]);

        fprintf(fstrain, "%.6f\t%.6f\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\n",
                acos(cos_th), phi_ang, eps[0], eps[1], eps[2], eps[3], eps[4], eps[5],
                theta_tr, sigma2);

        /* Multipole projection: integrate sigma2 * P_l(cos theta) */
        for (int l = 0; l <= 2; l++)
            c_l[l] += sigma2 * legendre_P(l, cos_th);
        norm += sigma2;
    }

    /* Normalize: fraction of total in each multipole */
    fprintf(fmulti, "# Multipole decomposition of |sigma|^2 on shell r=%.1f\n", R_shell);
    fprintf(fmulti, "l\tcoeff\tfraction\n");
    double sum_abs = fabs(c_l[0]) + fabs(c_l[1]) + fabs(c_l[2]);
    for (int l = 0; l <= 2; l++) {
        double frac = (sum_abs > 1e-30) ? fabs(c_l[l]) / sum_abs : 0.0;
        fprintf(fmulti, "%d\t%.6e\t%.6f\n", l, c_l[l], frac);
        printf("    l=%d: coefficient = %.6e, fraction = %.4f\n", l, c_l[l], frac);
    }
    printf("    Total |sigma|^2 on shell: %.6e\n", norm / N_sample);
}

/* ─── Compute TT strain at a point ─── */
/* For a point at (theta, phi_ang) on shell radius R from origin,
 * the TT projection removes the trace and the longitudinal components.
 * For simplicity, compute the full strain and extract h_+ and h_x
 * in the wave frame defined by the radial direction. */
static void compute_tt_strain(double theta, double phi_ang, double R,
                              double *h_plus, double *h_cross)
{
    double sin_th = sin(theta);
    double cos_th = cos(theta);
    double sin_ph = sin(phi_ang);
    double cos_ph = cos(phi_ang);

    /* Position on shell */
    double x = R * sin_th * cos_ph;
    double y = R * sin_th * sin_ph;
    double z = R * cos_th;

    /* Grid indices */
    int i = (int)((x + L) / dx + 0.5);
    int j = (int)((y + L) / dx + 0.5);
    int k = (int)((z + L) / dx + 0.5);

    *h_plus = 0.0;
    *h_cross = 0.0;

    if (i < 2 || i >= N-2 || j < 2 || j >= N-2 || k < 2 || k >= N-2)
        return;

    double eps[6];
    compute_strain(i, j, k, eps);

    /* Full 3x3 strain tensor */
    double e[3][3];
    e[0][0] = eps[0]; e[1][1] = eps[1]; e[2][2] = eps[2];
    e[0][1] = e[1][0] = eps[3];
    e[0][2] = e[2][0] = eps[4];
    e[1][2] = e[2][1] = eps[5];

    /* Radial unit vector */
    /* Build transverse basis:
     * e_theta = (cos_th*cos_ph, cos_th*sin_ph, -sin_th)
     * e_phi   = (-sin_ph, cos_ph, 0)
     */
    double et[3] = {cos_th*cos_ph, cos_th*sin_ph, -sin_th};
    double ep[3] = {-sin_ph, cos_ph, 0.0};

    /* TT projection: h_+ = e_{theta}^i e_{theta}^j eps_ij - e_{phi}^i e_{phi}^j eps_ij
     *                h_x = 2 e_{theta}^i e_{phi}^j eps_ij
     * These are the standard GW polarizations in the wave frame. */
    double htt = 0, hpp = 0, htp = 0;
    for (int a = 0; a < 3; a++)
        for (int b = 0; b < 3; b++) {
            htt += et[a] * et[b] * e[a][b];
            hpp += ep[a] * ep[b] * e[a][b];
            htp += et[a] * ep[b] * e[a][b];
        }

    *h_plus  = htt - hpp;
    *h_cross = 2.0 * htp;
}

/* ─── Main ─── */
int main(int argc, char **argv)
{
    parse_args(argc, argv);

    dx  = 2.0 * L / (N - 1);
    dx2 = dx * dx;
    m2  = mass * mass;
    dt  = cfl_frac * dx;
    Ngrid = (long)N * N * N;

    printf("=== V25: Three-Field Elastic Gravity in 3D ===\n");
    printf("Parameters:\n");
    printf("  mu=%.1f  kappa=%.1f  mass=%.3f  A=%.3f  sigma=%.3f\n",
           mu_pot, kappa, mass, A_init, sig_init);
    printf("  lambda_pw=%.3f  eta=%.3f  lambda_L=%.3f  alpha_g=%.4f\n",
           lambda_pw, eta, lambda_L, alpha_g);
    printf("  N=%d  L=%.1f  dx=%.4f  dt=%.5f\n", N, L, dx, dt);
    printf("  Ngrid=%ld (%.1f M)  Memory: %.1f MB\n",
           Ngrid, Ngrid/1e6, Ngrid*8.0*10/1e6);
    printf("  Threads: %d\n", omp_get_max_threads());
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

    /* Initialize: spherical Gaussians, all 3 fields equal */
    #pragma omp parallel for schedule(static)
    for (long idx = 0; idx < Ngrid; idx++) {
        int i = idx / (N * N);
        int j = (idx / N) % N;
        int k = idx % N;
        double x = -L + i * dx;
        double y = -L + j * dx;
        double z = -L + k * dx;
        double r2 = x*x + y*y + z*z;
        double g = exp(-r2 / (2.0 * sig_init * sig_init));
        for (int a = 0; a < 3; a++)
            phi[a][idx] = A_init * g;
    }

    double core_radius = 3.0 * sig_init;
    long ic = IDX(N/2, N/2, N/2);

    /* ================================================================
     * PHASE 1: 3D Elastic Oscillon
     * ================================================================ */
    printf("\n=== PHASE 1: 3D Elastic Oscillon (t=0..%.0f) ===\n", tfinal);
    fflush(stdout);

    int Nt1 = (int)(tfinal / dt) + 1;

    /* DFT storage */
    int max_dft = 50000;
    double *phi0_hist = malloc(max_dft * sizeof(double));
    double *t_hist = malloc(max_dft * sizeof(double));
    int n_dft = 0;

    int diag_every = Nt1 / 5000;
    if (diag_every < 1) diag_every = 1;
    int print_every = Nt1 / 50;
    if (print_every < 1) print_every = 1;
    int dft_every = Nt1 / max_dft;
    if (dft_every < 1) dft_every = 1;

    /* Output file */
    char path[600];
    snprintf(path, sizeof(path), "%s/v25_phase1_ts.tsv", outdir);
    FILE *f1 = fopen(path, "w");
    if (!f1) { fprintf(stderr, "Cannot open %s\n", path); return 1; }
    fprintf(f1, "time\tE_total\tE_kin\tE_grad\tE_mass\tE_pot\tE_pw\t"
                "fc\tpeak0\tpeak1\tpeak2\tphi0_ctr\n");

    compute_acc_elastic();

    double wall_start = omp_get_wtime();

    for (int n = 0; n <= Nt1; n++) {
        double t = n * dt;

        /* DFT history */
        if (n % dft_every == 0 && n_dft < max_dft) {
            phi0_hist[n_dft] = phi[0][ic];
            t_hist[n_dft] = t;
            n_dft++;
        }

        int do_diag  = (n % diag_every == 0);
        int do_print = (n % print_every == 0);

        if (do_diag || do_print) {
            Diag d = compute_diag(core_radius);

            if (do_diag)
                fprintf(f1, "%.4f\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t"
                            "%.6f\t%.6e\t%.6e\t%.6e\t%.6e\n",
                        t, d.Et, d.Ek, d.Eg, d.Em, d.Ep, d.Epw,
                        d.fc, d.peak[0], d.peak[1], d.peak[2], d.phi_center[0]);

            if (do_print) {
                double elapsed = omp_get_wtime() - wall_start;
                double frac = (double)n / Nt1;
                double eta_t = (frac > 0.001) ? elapsed * (1.0-frac)/frac : 0;
                printf("  t=%7.1f  pk=(%.4f,%.4f,%.4f)  E=%.2f  fc=%.4f  [%.0fs, ETA %.0fs]\n",
                       t, d.peak[0], d.peak[1], d.peak[2], d.Et, d.fc,
                       elapsed, eta_t);
                fflush(stdout);
            }
        }

        if (n == Nt1) break;
        verlet_step(0, 0.0);
    }

    fclose(f1);

    /* Compute breathing frequency from DFT */
    double peak_omega = find_peak_omega(phi0_hist, t_hist, n_dft, n_dft/2, 3.0*mass);
    printf("\nPhase 1 complete (%.1f sec)\n", omp_get_wtime() - wall_start);
    printf("  Peak omega = %.4f  (mass = %.4f, omega/m = %.4f)\n",
           peak_omega, mass, peak_omega / mass);
    printf("  Oscillon? %s\n",
           (peak_omega > 0.01 && peak_omega < mass) ? "YES" : "NO");

    {
        Diag d = compute_diag(core_radius);
        printf("  Final: E=%.2f  fc=%.4f  pk=(%.4f,%.4f,%.4f)\n",
               d.Et, d.fc, d.peak[0], d.peak[1], d.peak[2]);
    }
    fflush(stdout);

    /* ================================================================
     * PHASE 2: Strain Field Diagnostics
     * ================================================================ */
    printf("\n=== PHASE 2: Strain Field Diagnostics ===\n");
    fflush(stdout);

    double R_shell = 10.0;

    snprintf(path, sizeof(path), "%s/v25_phase2_strain.tsv", outdir);
    FILE *fstrain = fopen(path, "w");
    snprintf(path, sizeof(path), "%s/v25_phase2_multipoles.tsv", outdir);
    FILE *fmulti = fopen(path, "w");

    if (!fstrain || !fmulti) {
        fprintf(stderr, "Cannot open Phase 2 output files\n");
        return 1;
    }

    compute_strain_multipoles(R_shell, fstrain, fmulti);

    fclose(fstrain);
    fclose(fmulti);

    printf("Phase 2 complete.\n");
    fflush(stdout);

    /* ================================================================
     * PHASE 3: Self-Consistent Metric
     * ================================================================ */
    printf("\n=== PHASE 3: Self-Consistent Metric (t=%.0f..%.0f) ===\n",
           tfinal, 2.0*tfinal);
    fflush(stdout);

    /* Run Phase 3 for another tfinal time units with metric ramp-up */
    double t3_start = tfinal;
    double t3_ramp_end = tfinal + 100.0;  /* ramp alpha_g over 100 time units */
    double t3_end = 2.0 * tfinal;
    int Nt3 = (int)((t3_end - t3_start) / dt) + 1;

    snprintf(path, sizeof(path), "%s/v25_phase3_ts.tsv", outdir);
    FILE *f3 = fopen(path, "w");
    if (!f3) { fprintf(stderr, "Cannot open Phase 3 output\n"); return 1; }
    fprintf(f3, "time\tE_total\tfc\tpeak0\talpha_g_eff\tphi0_ctr\n");

    int diag3_every = Nt3 / 5000;
    if (diag3_every < 1) diag3_every = 1;
    int print3_every = Nt3 / 50;
    if (print3_every < 1) print3_every = 1;

    /* Reset DFT storage for Phase 3 */
    n_dft = 0;
    int dft3_every = Nt3 / max_dft;
    if (dft3_every < 1) dft3_every = 1;

    double wall3 = omp_get_wtime();

    /* Initial acceleration with zero metric coupling */
    compute_acc_metric(0.0);

    for (int n = 0; n <= Nt3; n++) {
        double t = t3_start + n * dt;
        double ag_eff = alpha_g * hermite_step(t, t3_start, t3_ramp_end);

        /* DFT history */
        if (n % dft3_every == 0 && n_dft < max_dft) {
            phi0_hist[n_dft] = phi[0][ic];
            t_hist[n_dft] = t;
            n_dft++;
        }

        int do_diag  = (n % diag3_every == 0);
        int do_print = (n % print3_every == 0);

        if (do_diag || do_print) {
            Diag d = compute_diag(core_radius);

            if (do_diag)
                fprintf(f3, "%.4f\t%.6e\t%.6f\t%.6e\t%.6e\t%.6e\n",
                        t, d.Et, d.fc, d.peak[0], ag_eff, d.phi_center[0]);

            if (do_print) {
                double elapsed = omp_get_wtime() - wall3;
                double frac = (double)n / Nt3;
                double eta_t = (frac > 0.001) ? elapsed*(1.0-frac)/frac : 0;
                printf("  t=%7.1f  ag=%.4f  pk=%.4f  E=%.2f  fc=%.4f  [%.0fs, ETA %.0fs]\n",
                       t, ag_eff, d.peak[0], d.Et, d.fc, elapsed, eta_t);
                fflush(stdout);
            }
        }

        if (n == Nt3) break;

        /* Verlet with metric */
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
        /* 3. Recompute with potentially updated ag_eff */
        double ag_next = alpha_g * hermite_step(t + dt, t3_start, t3_ramp_end);
        compute_acc_metric(ag_next);
        /* 4. Half-kick */
        for (int a = 0; a < 3; a++) {
            #pragma omp parallel for schedule(static)
            for (long idx = 0; idx < Ngrid; idx++)
                vel[a][idx] += 0.5 * dt * acc[a][idx];
        }
        /* 5. Absorbing boundary */
        for (int a = 0; a < 3; a++) {
            #pragma omp parallel for schedule(static)
            for (long idx = 0; idx < Ngrid; idx++) {
                vel[a][idx] *= damp[idx];
                phi[a][idx] *= damp[idx];
            }
        }
    }

    fclose(f3);

    double peak_omega3 = find_peak_omega(phi0_hist, t_hist, n_dft, n_dft/2, 3.0*mass);
    printf("\nPhase 3 complete (%.1f sec)\n", omp_get_wtime() - wall3);
    printf("  Peak omega = %.4f  (omega/m = %.4f)\n", peak_omega3, peak_omega3/mass);
    {
        Diag d = compute_diag(core_radius);
        printf("  Final: E=%.2f  fc=%.4f  pk=(%.4f,%.4f,%.4f)\n",
               d.Et, d.fc, d.peak[0], d.peak[1], d.peak[2]);
    }
    fflush(stdout);

    /* ================================================================
     * PHASE 4: Gravitational Wave Detection
     * ================================================================ */
    printf("\n=== PHASE 4: Gravitational Wave Detection ===\n");
    fflush(stdout);

    /* Expand grid for far-field detection if needed */
    /* We keep the same grid but use L=20 for Phase 4 if the user wants.
     * For simplicity, continue on same grid (L=15 gives R_shell=12 at ~38 grid pts from center).
     * That's within the undamped region (R_abs_inner = 0.70*15 = 10.5).
     * We need R_shell < R_abs_inner, so use R_shell = 8 or rebuild grid.
     * Actually: let's rebuild with L=20 for Phase 4. */

    double L_old = L;
    L = 20.0;
    dx = 2.0 * L / (N - 1);
    dx2 = dx * dx;
    dt = cfl_frac * dx;

    /* Save current fields on old grid, re-initialize on new grid */
    double *phi_save[3], *vel_save[3];
    double dx_old = 2.0 * L_old / (N - 1);
    for (int a = 0; a < 3; a++) {
        phi_save[a] = phi[a];
        vel_save[a] = vel[a];
        phi[a] = calloc(Ngrid, sizeof(double));
        vel[a] = calloc(Ngrid, sizeof(double));
        if (!phi[a] || !vel[a]) {
            fprintf(stderr, "Allocation failed for Phase 4\n");
            return 1;
        }
    }

    /* Interpolate from old grid to new grid (nearest neighbor) */
    #pragma omp parallel for schedule(static)
    for (long idx = 0; idx < Ngrid; idx++) {
        int ii = idx / (N * N);
        int jj = (idx / N) % N;
        int kk = idx % N;
        double x = -L + ii * dx;
        double y = -L + jj * dx;
        double z = -L + kk * dx;

        /* Find nearest old grid point */
        int io = (int)((x + L_old) / dx_old + 0.5);
        int jo = (int)((y + L_old) / dx_old + 0.5);
        int ko = (int)((z + L_old) / dx_old + 0.5);

        if (io >= 0 && io < N && jo >= 0 && jo < N && ko >= 0 && ko < N) {
            long oidx = IDX(io, jo, ko);
            for (int a = 0; a < 3; a++) {
                phi[a][idx] = phi_save[a][oidx];
                vel[a][idx] = vel_save[a][oidx];
            }
        }
    }

    for (int a = 0; a < 3; a++) {
        free(phi_save[a]);
        free(vel_save[a]);
    }

    /* Rebuild damping for new L */
    R_abs_inner = L * 0.70;
    R_abs_outer = L * 0.95;

    #pragma omp parallel for schedule(static)
    for (long idx = 0; idx < Ngrid; idx++) {
        int ii = idx / (N * N);
        int jj = (idx / N) % N;
        int kk = idx % N;
        double x = -L + ii * dx;
        double y = -L + jj * dx;
        double z = -L + kk * dx;
        double r = sqrt(x*x + y*y + z*z);
        if (r > R_abs_inner) {
            double f = (r - R_abs_inner) / (R_abs_outer - R_abs_inner);
            if (f > 1.0) f = 1.0;
            damp[idx] = 1.0 - 0.98 * f * f;
        } else {
            damp[idx] = 1.0;
        }
    }

    /* Apply velocity kick in z-direction to all 3 fields */
    printf("  Applying v_kick = %.4f in z-direction\n", v_kick);
    #pragma omp parallel for schedule(static)
    for (long idx = 0; idx < Ngrid; idx++) {
        int ii = idx / (N * N);
        int jj = (idx / N) % N;
        int kk = idx % N;
        /* Boost: vel_a += -v_kick * d_z(phi_a) (Lorentz boost approximation) */
        for (int a = 0; a < 3; a++) {
            if (kk >= 1 && kk < N-1) {
                double dz_phi = (phi[a][IDX(ii,jj,kk+1)] - phi[a][IDX(ii,jj,kk-1)]) / (2.0*dx);
                vel[a][idx] += -v_kick * dz_phi;
            }
        }
    }

    /* Recompute acceleration with full metric coupling */
    for (int a = 0; a < 3; a++) {
        memset(acc[a], 0, Ngrid * sizeof(double));
    }
    compute_acc_metric(alpha_g);

    /* Track TT strain at 4 angles, R=12 */
    double R_gw = 12.0;
    /* Angles: (theta, phi_ang) */
    double gw_theta[4] = {0.001,           M_PI/2.0,    M_PI/2.0,      M_PI/4.0};
    double gw_phi[4]   = {0.0,             0.0,         M_PI/4.0,      0.0};
    const char *gw_label[4] = {"z-axis(th=0)", "x-axis(th=pi/2,ph=0)",
                                "xy-45(th=pi/2,ph=pi/4)", "45-deg(th=pi/4,ph=0)"};

    double t4_end = 200.0;  /* Phase 4 run time */
    int Nt4 = (int)(t4_end / dt) + 1;

    snprintf(path, sizeof(path), "%s/v25_phase4_gw.tsv", outdir);
    FILE *f4 = fopen(path, "w");
    if (!f4) { fprintf(stderr, "Cannot open Phase 4 output\n"); return 1; }
    fprintf(f4, "time\th+_z\thx_z\th+_x\thx_x\th+_xy45\thx_xy45\th+_45\thx_45\n");

    int diag4_every = Nt4 / 2000;
    if (diag4_every < 1) diag4_every = 1;
    int print4_every = Nt4 / 20;
    if (print4_every < 1) print4_every = 1;

    double wall4 = omp_get_wtime();

    /* Track peak strain for summary */
    double max_hp[4] = {0,0,0,0}, max_hx[4] = {0,0,0,0};

    for (int n = 0; n <= Nt4; n++) {
        double t = n * dt;

        if (n % diag4_every == 0) {
            fprintf(f4, "%.4f", t);
            for (int g = 0; g < 4; g++) {
                double hp, hx;
                compute_tt_strain(gw_theta[g], gw_phi[g], R_gw, &hp, &hx);
                fprintf(f4, "\t%.6e\t%.6e", hp, hx);
                if (fabs(hp) > max_hp[g]) max_hp[g] = fabs(hp);
                if (fabs(hx) > max_hx[g]) max_hx[g] = fabs(hx);
            }
            fprintf(f4, "\n");
        }

        if (n % print4_every == 0) {
            Diag d = compute_diag(core_radius);
            double elapsed = omp_get_wtime() - wall4;
            double frac = (double)n / Nt4;
            double eta_t = (frac > 0.001) ? elapsed*(1.0-frac)/frac : 0;
            printf("  t=%6.1f  pk=%.4f  E=%.2f  fc=%.4f  [%.0fs, ETA %.0fs]\n",
                   t, d.peak[0], d.Et, d.fc, elapsed, eta_t);
            fflush(stdout);
        }

        if (n == Nt4) break;

        /* Verlet with full metric */
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
        compute_acc_metric(alpha_g);
        for (int a = 0; a < 3; a++) {
            #pragma omp parallel for schedule(static)
            for (long idx = 0; idx < Ngrid; idx++)
                vel[a][idx] += 0.5 * dt * acc[a][idx];
        }
        for (int a = 0; a < 3; a++) {
            #pragma omp parallel for schedule(static)
            for (long idx = 0; idx < Ngrid; idx++) {
                vel[a][idx] *= damp[idx];
                phi[a][idx] *= damp[idx];
            }
        }
    }

    fclose(f4);

    printf("\nPhase 4 complete (%.1f sec)\n", omp_get_wtime() - wall4);
    printf("  Peak TT strain at R=%.0f:\n", R_gw);
    for (int g = 0; g < 4; g++) {
        printf("    %s:  |h+|_max = %.4e   |hx|_max = %.4e\n",
               gw_label[g], max_hp[g], max_hx[g]);
    }

    /* Angular pattern analysis */
    printf("\n  Angular pattern:\n");
    if (max_hp[0] < 0.1 * max_hp[1] && max_hp[1] > 1e-10) {
        printf("    h+ suppressed along z-axis (forward) -> consistent with quadrupolar\n");
    } else if (max_hp[0] > 1e-10) {
        printf("    h+ NOT suppressed along z-axis -> not purely quadrupolar\n");
    } else {
        printf("    h+ too small to determine pattern\n");
    }
    if (max_hp[1] > max_hp[0] && max_hp[3] > 0.3*max_hp[1] && max_hp[3] < 0.9*max_hp[1]) {
        printf("    Intermediate h+ at 45 deg -> spin-2 pattern likely\n");
    }

    /* ================================================================
     * FINAL SUMMARY
     * ================================================================ */
    double wall_total = omp_get_wtime() - wall_start;

    printf("\n=== FINAL SUMMARY ===\n");
    printf("  Total wall time: %.1f sec (%.1f min)\n", wall_total, wall_total/60);
    printf("  Phase 1 (elastic oscillon):   omega=%.4f  omega/m=%.4f\n",
           peak_omega, peak_omega/mass);
    printf("  Phase 3 (self-consistent):    omega=%.4f  omega/m=%.4f\n",
           peak_omega3, peak_omega3/mass);
    printf("  Phase 4 (GW):  max |h+|_equatorial = %.4e\n", max_hp[1]);
    printf("===\n");

    /* Cleanup */
    for (int a = 0; a < 3; a++) {
        free(phi[a]); free(vel[a]); free(acc[a]);
    }
    free(damp); free(phi0_hist); free(t_hist);

    return 0;
}
