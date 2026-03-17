/*
 * binary.c — Binary Oscillon Gravitational Wave Emission (V25 Phase 5)
 *
 * Two oscillons orbiting each other in the xy-plane.
 * Measures TT strain at multiple angles to detect spin-2 (quadrupolar)
 * gravitational wave emission at 2*Omega (twice the orbital frequency).
 *
 * Based on v25.c (Phases 1-4). All couplings preserved:
 *   elastic eta, Lame lambda_L, pairwise lambda_pw,
 *   triple product mu/kappa, self-consistent metric alpha_g.
 *
 * Strategy:
 *   1. Equilibrate a single oscillon at the origin (t=0..t_equil)
 *   2. Copy the equilibrated profile to two locations at x=+/-D/2
 *   3. Apply tangential velocities v_y = +/-v_orb
 *   4. Evolve with self-consistent metric, tracking positions and GW strain
 *
 * Compile: gcc -O3 -fopenmp -Wall -o binary src/binary.c -lm
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
static double alpha_g   = 0.001;

static double D_sep     = 12.0;   /* initial separation */
static double v_orb     = 0.1;    /* orbital velocity */
static double t_equil   = 200.0;  /* equilibration time for single oscillon */
static double t_binary  = 500.0;  /* binary evolution time */

static int    N         = 128;
static double L         = 30.0;
static double cfl_frac  = 0.20;
static char   outdir[512] = "data";

/* For equilibration: smaller grid */
static int    N_eq      = 96;
static double L_eq      = 15.0;

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
        else if (!strcmp(argv[i], "-D"))        D_sep     = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-vorb"))     v_orb     = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-teq"))      t_equil   = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-tbin"))     t_binary  = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-lpw"))      lambda_pw = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-eta"))      eta       = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-lL"))       lambda_L  = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-ag"))       alpha_g   = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-cfl"))      cfl_frac  = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-o"))        strncpy(outdir, argv[i+1], sizeof(outdir)-1);
        else { fprintf(stderr, "Unknown flag: %s\n", argv[i]); exit(1); }
    }
}

/* ─── Grid globals ─── */
static double *phi[3], *vel[3], *acc[3];
static double *damp;
static double dx, dx2, m2, dt;
static long Ngrid;
static int Ncur;    /* current grid size (N_eq during equil, N during binary) */
static double Lcur; /* current box half-size */

#define IDX(i,j,k) ((long)(i)*Ncur*Ncur + (long)(j)*Ncur + (long)(k))

/* ─── Hermite smooth step ─── */
static double hermite_step(double t, double t0, double t1)
{
    if (t <= t0) return 0.0;
    if (t >= t1) return 1.0;
    double s = (t - t0) / (t1 - t0);
    return s * s * (3.0 - 2.0 * s);
}

/* ─── Compute acceleration with self-consistent metric ─── */
static void compute_acc_metric(double ag_eff)
{
    double eta_lL = eta + lambda_L;

    for (int a = 0; a < 3; a++) {
        #pragma omp parallel for schedule(static)
        for (int i = 2; i < Ncur-2; i++) {
            for (int j = 2; j < Ncur-2; j++) {
                for (int k = 2; k < Ncur-2; k++) {
                    long idx = IDX(i,j,k);

                    /* 7-point Laplacian */
                    double lapl = (phi[a][IDX(i+1,j,k)] + phi[a][IDX(i-1,j,k)]
                                 + phi[a][IDX(i,j+1,k)] + phi[a][IDX(i,j-1,k)]
                                 + phi[a][IDX(i,j,k+1)] + phi[a][IDX(i,j,k-1)]
                                 - 6.0 * phi[a][idx]) / dx2;

                    /* Metric correction (diagonal approximation) */
                    double metric_corr = 0.0;
                    double theta = 0.0;
                    if (ag_eff > 1e-15) {
                        double dx_phi0 = (phi[0][IDX(i+1,j,k)] - phi[0][IDX(i-1,j,k)]) / (2.0*dx);
                        double dy_phi1 = (phi[1][IDX(i,j+1,k)] - phi[1][IDX(i,j-1,k)]) / (2.0*dx);
                        double dz_phi2 = (phi[2][IDX(i,j,k+1)] - phi[2][IDX(i,j,k-1)]) / (2.0*dx);

                        double h11 = 2.0 * dx_phi0;
                        double h22 = 2.0 * dy_phi1;
                        double h33 = 2.0 * dz_phi2;
                        theta = dx_phi0 + dy_phi1 + dz_phi2;

                        double d2x = (phi[a][IDX(i+1,j,k)] - 2.0*phi[a][idx]
                                    + phi[a][IDX(i-1,j,k)]) / dx2;
                        double d2y = (phi[a][IDX(i,j+1,k)] - 2.0*phi[a][idx]
                                    + phi[a][IDX(i,j-1,k)]) / dx2;
                        double d2z = (phi[a][IDX(i,j,k+1)] - 2.0*phi[a][idx]
                                    + phi[a][IDX(i,j,k-1)]) / dx2;

                        metric_corr = ag_eff * (h11 * d2x + h22 * d2y + h33 * d2z);
                    }

                    double m2_eff = m2 * (1.0 + ag_eff * theta);

                    /* Triple-product potential */
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

                    /* Cross-gradient + Lame: (eta+lambda_L) * d_a(div phi) */
                    double d_a_div = 0.0;

                    /* b=0: d_a(d_x phi_0) */
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

                    /* b=1: d_a(d_y phi_1) */
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

                    /* b=2: d_a(d_z phi_2) */
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
            int ii = idx / (Ncur * Ncur);
            int jj = (idx / Ncur) % Ncur;
            int kk = idx % Ncur;
            if (ii < 2 || ii >= Ncur-2 || jj < 2 || jj >= Ncur-2 || kk < 2 || kk >= Ncur-2)
                acc[a][idx] = 0.0;
        }
    }
}

/* ─── Compute strain tensor at a grid point ─── */
static void compute_strain(int i, int j, int k, double eps[6])
{
    double dphi[3][3];
    for (int a = 0; a < 3; a++) {
        dphi[a][0] = (phi[a][IDX(i+1,j,k)] - phi[a][IDX(i-1,j,k)]) / (2.0*dx);
        dphi[a][1] = (phi[a][IDX(i,j+1,k)] - phi[a][IDX(i,j-1,k)]) / (2.0*dx);
        dphi[a][2] = (phi[a][IDX(i,j,k+1)] - phi[a][IDX(i,j,k-1)]) / (2.0*dx);
    }
    eps[0] = dphi[0][0];
    eps[1] = dphi[1][1];
    eps[2] = dphi[2][2];
    eps[3] = 0.5*(dphi[1][0] + dphi[0][1]);
    eps[4] = 0.5*(dphi[2][0] + dphi[0][2]);
    eps[5] = 0.5*(dphi[2][1] + dphi[1][2]);
}

/* ─── Compute TT strain h+ and hx at a point on sphere ─── */
static void compute_tt_strain(double theta, double phi_ang, double R,
                              double *h_plus, double *h_cross)
{
    double sin_th = sin(theta);
    double cos_th = cos(theta);
    double sin_ph = sin(phi_ang);
    double cos_ph = cos(phi_ang);

    double x = R * sin_th * cos_ph;
    double y = R * sin_th * sin_ph;
    double z = R * cos_th;

    int i = (int)((x + Lcur) / dx + 0.5);
    int j = (int)((y + Lcur) / dx + 0.5);
    int k = (int)((z + Lcur) / dx + 0.5);

    *h_plus = 0.0;
    *h_cross = 0.0;

    if (i < 2 || i >= Ncur-2 || j < 2 || j >= Ncur-2 || k < 2 || k >= Ncur-2)
        return;

    double eps[6];
    compute_strain(i, j, k, eps);

    double e[3][3];
    e[0][0] = eps[0]; e[1][1] = eps[1]; e[2][2] = eps[2];
    e[0][1] = e[1][0] = eps[3];
    e[0][2] = e[2][0] = eps[4];
    e[1][2] = e[2][1] = eps[5];

    /* Transverse basis vectors */
    double et[3] = {cos_th*cos_ph, cos_th*sin_ph, -sin_th};
    double ep[3] = {-sin_ph, cos_ph, 0.0};

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

/* ─── Find peak frequency via DFT ─── */
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

/* ─── Energy-weighted centroid for one half of domain ─── */
/* side=+1: x>0 (oscillon A), side=-1: x<0 (oscillon B) */
static void find_centroid(int side, double *cx, double *cy, double *cz, double *E_half)
{
    double sx = 0, sy = 0, sz = 0, se = 0;

    #pragma omp parallel
    {
        double lsx=0, lsy=0, lsz=0, lse=0;
        #pragma omp for schedule(static) nowait
        for (int i = 2; i < Ncur-2; i++) {
            double x = -Lcur + i * dx;
            if (side > 0 && x < 0) continue;
            if (side < 0 && x >= 0) continue;
            for (int j = 2; j < Ncur-2; j++) {
                double y = -Lcur + j * dx;
                for (int k = 2; k < Ncur-2; k++) {
                    double z = -Lcur + k * dx;
                    long idx = IDX(i,j,k);

                    double e_loc = 0;
                    for (int a = 0; a < 3; a++) {
                        double v2 = vel[a][idx] * vel[a][idx];
                        e_loc += 0.5 * v2;

                        double gx = (phi[a][IDX(i+1,j,k)] - phi[a][IDX(i-1,j,k)]) / (2*dx);
                        double gy = (phi[a][IDX(i,j+1,k)] - phi[a][IDX(i,j-1,k)]) / (2*dx);
                        double gz = (phi[a][IDX(i,j,k+1)] - phi[a][IDX(i,j,k-1)]) / (2*dx);
                        e_loc += 0.5 * (gx*gx + gy*gy + gz*gz);
                        e_loc += 0.5 * m2 * phi[a][idx] * phi[a][idx];
                    }
                    double p0 = phi[0][idx], p1 = phi[1][idx], p2 = phi[2][idx];
                    double P = p0*p1*p2, P2 = P*P;
                    e_loc += 0.5 * mu_pot * P2 / (1.0 + kappa*P2);
                    e_loc += lambda_pw * (p0*p1 + p1*p2 + p2*p0);

                    if (e_loc > 0) {
                        double dV = dx*dx*dx;
                        lsx += x * e_loc * dV;
                        lsy += y * e_loc * dV;
                        lsz += z * e_loc * dV;
                        lse += e_loc * dV;
                    }
                }
            }
        }
        #pragma omp critical
        { sx += lsx; sy += lsy; sz += lsz; se += lse; }
    }

    *E_half = se;
    if (se > 1e-10) {
        *cx = sx / se; *cy = sy / se; *cz = sz / se;
    } else {
        *cx = (side > 0) ? D_sep/2.0 : -D_sep/2.0;
        *cy = 0; *cz = 0;
    }
}

/* ─── Total energy diagnostic ─── */
static double compute_total_energy(void)
{
    double Etot = 0;
    #pragma omp parallel for reduction(+:Etot) schedule(static)
    for (int i = 2; i < Ncur-2; i++) {
        for (int j = 2; j < Ncur-2; j++) {
            for (int k = 2; k < Ncur-2; k++) {
                long idx = IDX(i,j,k);
                double dV = dx*dx*dx;
                double e_loc = 0;
                for (int a = 0; a < 3; a++) {
                    e_loc += 0.5 * vel[a][idx] * vel[a][idx];
                    double gx = (phi[a][IDX(i+1,j,k)] - phi[a][IDX(i-1,j,k)]) / (2*dx);
                    double gy = (phi[a][IDX(i,j+1,k)] - phi[a][IDX(i,j-1,k)]) / (2*dx);
                    double gz = (phi[a][IDX(i,j,k+1)] - phi[a][IDX(i,j,k-1)]) / (2*dx);
                    e_loc += 0.5 * (gx*gx + gy*gy + gz*gz);
                    e_loc += 0.5 * m2 * phi[a][idx] * phi[a][idx];
                }
                double p0 = phi[0][idx], p1 = phi[1][idx], p2 = phi[2][idx];
                double P = p0*p1*p2, P2 = P*P;
                e_loc += 0.5 * mu_pot * P2 / (1.0 + kappa*P2);
                e_loc += lambda_pw * (p0*p1 + p1*p2 + p2*p0);
                Etot += e_loc * dV;
            }
        }
    }
    return Etot;
}

/* ─── Multipole decomposition of |sigma^TT|^2 on shell ─── */
static void compute_binary_multipoles(double R_shell, double c_l[3], double *total_sigma2)
{
    int N_sample = 400;
    double golden_ratio = (1.0 + sqrt(5.0)) / 2.0;

    for (int l = 0; l < 3; l++) c_l[l] = 0;
    *total_sigma2 = 0;

    for (int ns = 0; ns < N_sample; ns++) {
        double cos_th = 1.0 - 2.0 * (ns + 0.5) / N_sample;
        double sin_th = sqrt(1.0 - cos_th * cos_th);
        double phi_ang = 2.0 * M_PI * ns / golden_ratio;

        double x = R_shell * sin_th * cos(phi_ang);
        double y = R_shell * sin_th * sin(phi_ang);
        double z = R_shell * cos_th;

        int i = (int)((x + Lcur) / dx + 0.5);
        int j_idx = (int)((y + Lcur) / dx + 0.5);
        int k = (int)((z + Lcur) / dx + 0.5);

        if (i < 2 || i >= Ncur-2 || j_idx < 2 || j_idx >= Ncur-2 || k < 2 || k >= Ncur-2)
            continue;

        double eps[6];
        compute_strain(i, j_idx, k, eps);

        double theta_tr = eps[0] + eps[1] + eps[2];
        double sig[6];
        sig[0] = eps[0] - theta_tr / 3.0;
        sig[1] = eps[1] - theta_tr / 3.0;
        sig[2] = eps[2] - theta_tr / 3.0;
        sig[3] = eps[3]; sig[4] = eps[4]; sig[5] = eps[5];

        double sigma2 = sig[0]*sig[0] + sig[1]*sig[1] + sig[2]*sig[2]
                       + 2.0*(sig[3]*sig[3] + sig[4]*sig[4] + sig[5]*sig[5]);

        for (int l = 0; l <= 2; l++)
            c_l[l] += sigma2 * legendre_P(l, cos_th);
        *total_sigma2 += sigma2;
    }

    /* Normalize by number of samples */
    *total_sigma2 /= N_sample;
}

/* ─── Allocate fields for current grid ─── */
static void alloc_fields(void)
{
    Ngrid = (long)Ncur * Ncur * Ncur;
    for (int a = 0; a < 3; a++) {
        phi[a] = calloc(Ngrid, sizeof(double));
        vel[a] = calloc(Ngrid, sizeof(double));
        acc[a] = calloc(Ngrid, sizeof(double));
        if (!phi[a] || !vel[a] || !acc[a]) {
            fprintf(stderr, "Allocation failed\n"); exit(1);
        }
    }
    damp = malloc(Ngrid * sizeof(double));
    if (!damp) { fprintf(stderr, "Alloc failed for damp\n"); exit(1); }
}

/* ─── Free fields ─── */
static void free_fields(void)
{
    for (int a = 0; a < 3; a++) {
        free(phi[a]); free(vel[a]); free(acc[a]);
        phi[a] = vel[a] = acc[a] = NULL;
    }
    free(damp); damp = NULL;
}

/* ─── Setup damping layer ─── */
static void setup_damping(void)
{
    double R_abs_inner = Lcur * 0.70;
    double R_abs_outer = Lcur * 0.95;

    #pragma omp parallel for schedule(static)
    for (long idx = 0; idx < Ngrid; idx++) {
        int i = idx / (Ncur * Ncur);
        int j = (idx / Ncur) % Ncur;
        int k = idx % Ncur;
        double x = -Lcur + i * dx;
        double y = -Lcur + j * dx;
        double z = -Lcur + k * dx;
        double r = sqrt(x*x + y*y + z*z);
        if (r > R_abs_inner) {
            double f = (r - R_abs_inner) / (R_abs_outer - R_abs_inner);
            if (f > 1.0) f = 1.0;
            damp[idx] = 1.0 - 0.98 * f * f;
        } else {
            damp[idx] = 1.0;
        }
    }
}

/* ─── Verlet step ─── */
static void verlet_step(double ag_eff)
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
    compute_acc_metric(ag_eff);
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

/* ─── Main ─── */
int main(int argc, char **argv)
{
    parse_args(argc, argv);

    m2 = mass * mass;

    printf("=== V25 Phase 5: Binary Oscillon GW Emission ===\n");
    printf("Parameters:\n");
    printf("  mu=%.1f  kappa=%.1f  mass=%.3f  A=%.3f  sigma=%.3f\n",
           mu_pot, kappa, mass, A_init, sig_init);
    printf("  lambda_pw=%.3f  eta=%.3f  lambda_L=%.3f  alpha_g=%.4f\n",
           lambda_pw, eta, lambda_L, alpha_g);
    printf("  D=%.1f  v_orb=%.4f\n", D_sep, v_orb);
    printf("  Omega=%.6f  T_orb=%.1f\n", v_orb/(D_sep/2.0),
           2.0*M_PI/(v_orb/(D_sep/2.0)));
    printf("  t_equil=%.0f  t_binary=%.0f\n", t_equil, t_binary);
    printf("  Threads: %d\n", omp_get_max_threads());
    fflush(stdout);

    double wall_start = omp_get_wtime();

    /* ================================================================
     * STEP 1: Equilibrate single oscillon on small grid
     * ================================================================ */
    printf("\n=== STEP 1: Equilibrate single oscillon (N=%d, L=%.0f, t=%.0f) ===\n",
           N_eq, L_eq, t_equil);
    fflush(stdout);

    Ncur = N_eq;
    Lcur = L_eq;
    dx = 2.0 * Lcur / (Ncur - 1);
    dx2 = dx * dx;
    dt = cfl_frac * dx;

    alloc_fields();
    setup_damping();

    printf("  dx=%.4f  dt=%.5f  Ngrid=%ld (%.1f M)\n", dx, dt, Ngrid, Ngrid/1e6);
    fflush(stdout);

    /* Initialize: spherical Gaussian */
    #pragma omp parallel for schedule(static)
    for (long idx = 0; idx < Ngrid; idx++) {
        int i = idx / (Ncur * Ncur);
        int j = (idx / Ncur) % Ncur;
        int k = idx % Ncur;
        double x = -Lcur + i * dx;
        double y = -Lcur + j * dx;
        double z = -Lcur + k * dx;
        double r2 = x*x + y*y + z*z;
        double g = exp(-r2 / (2.0 * sig_init * sig_init));
        for (int a = 0; a < 3; a++)
            phi[a][idx] = A_init * g;
    }

    compute_acc_metric(0.0);  /* no metric during equilibration */

    int Nt_eq = (int)(t_equil / dt) + 1;
    int print_eq = Nt_eq / 20;
    if (print_eq < 1) print_eq = 1;

    double wall_eq = omp_get_wtime();
    for (int n = 0; n <= Nt_eq; n++) {
        if (n % print_eq == 0) {
            double t = n * dt;
            long ic = IDX(Ncur/2, Ncur/2, Ncur/2);
            double phi0c = phi[0][ic];
            double elapsed = omp_get_wtime() - wall_eq;
            double frac = (double)n / Nt_eq;
            double eta_t = (frac > 0.001) ? elapsed*(1.0-frac)/frac : 0;
            printf("  t=%7.1f  phi0_ctr=%.6f  [%.0fs, ETA %.0fs]\n",
                   t, phi0c, elapsed, eta_t);
            fflush(stdout);
        }
        if (n == Nt_eq) break;
        verlet_step(0.0);
    }

    printf("  Equilibration complete (%.1f sec)\n", omp_get_wtime() - wall_eq);

    /* Save equilibrated profile on small grid */
    double *eq_phi[3], *eq_vel[3];
    for (int a = 0; a < 3; a++) {
        eq_phi[a] = malloc(Ngrid * sizeof(double));
        eq_vel[a] = malloc(Ngrid * sizeof(double));
        memcpy(eq_phi[a], phi[a], Ngrid * sizeof(double));
        memcpy(eq_vel[a], vel[a], Ngrid * sizeof(double));
    }
    int Ncur_eq = Ncur;
    double dx_eq = dx;

    free_fields();

    /* ================================================================
     * STEP 2: Place two oscillons on large grid
     * ================================================================ */
    printf("\n=== STEP 2: Place binary on large grid (N=%d, L=%.0f) ===\n", N, L);
    fflush(stdout);

    Ncur = N;
    Lcur = L;
    dx = 2.0 * Lcur / (Ncur - 1);
    dx2 = dx * dx;
    dt = cfl_frac * dx;

    alloc_fields();
    setup_damping();

    printf("  dx=%.4f  dt=%.5f  Ngrid=%ld (%.1f M)\n", dx, dt, Ngrid, Ngrid/1e6);
    printf("  R_abs_inner=%.1f  (measurement sphere R=20 is OK)\n", Lcur*0.70);
    fflush(stdout);

    /* Place two copies of equilibrated oscillon at x = +/- D/2 */
    #pragma omp parallel for schedule(static)
    for (long idx = 0; idx < Ngrid; idx++) {
        int i = idx / (Ncur * Ncur);
        int j = (idx / Ncur) % Ncur;
        int k = idx % Ncur;
        double x = -Lcur + i * dx;
        double y = -Lcur + j * dx;
        double z = -Lcur + k * dx;

        for (int a = 0; a < 3; a++) {
            phi[a][idx] = 0.0;
            vel[a][idx] = 0.0;
        }

        /* Oscillon A at (+D/2, 0, 0) */
        double xA = x - D_sep / 2.0;
        int ioA = (int)((xA + L_eq) / dx_eq + 0.5);
        int joA = (int)((y + L_eq) / dx_eq + 0.5);
        int koA = (int)((z + L_eq) / dx_eq + 0.5);

        if (ioA >= 0 && ioA < Ncur_eq && joA >= 0 && joA < Ncur_eq && koA >= 0 && koA < Ncur_eq) {
            long oidx = (long)ioA * Ncur_eq * Ncur_eq + (long)joA * Ncur_eq + koA;
            for (int a = 0; a < 3; a++) {
                phi[a][idx] += eq_phi[a][oidx];
                vel[a][idx] += eq_vel[a][oidx];
            }
        }

        /* Oscillon B at (-D/2, 0, 0) */
        double xB = x + D_sep / 2.0;
        int ioB = (int)((xB + L_eq) / dx_eq + 0.5);
        int joB = (int)((y + L_eq) / dx_eq + 0.5);
        int koB = (int)((z + L_eq) / dx_eq + 0.5);

        if (ioB >= 0 && ioB < Ncur_eq && joB >= 0 && joB < Ncur_eq && koB >= 0 && koB < Ncur_eq) {
            long oidx = (long)ioB * Ncur_eq * Ncur_eq + (long)joB * Ncur_eq + koB;
            for (int a = 0; a < 3; a++) {
                phi[a][idx] += eq_phi[a][oidx];
                vel[a][idx] += eq_vel[a][oidx];
            }
        }
    }

    /* Free equilibrated profiles */
    for (int a = 0; a < 3; a++) {
        free(eq_phi[a]); free(eq_vel[a]);
    }

    /* Apply tangential velocity: A gets +v_orb in y, B gets -v_orb in y */
    /* Boost: vel_a += -v * d_y(phi_a)  for oscillon at that location */
    printf("  Applying orbital velocities: v_orb = %.4f\n", v_orb);
    #pragma omp parallel for schedule(static)
    for (long idx = 0; idx < Ngrid; idx++) {
        int i = idx / (Ncur * Ncur);
        int j = (idx / Ncur) % Ncur;
        int k = idx % Ncur;
        double x = -Lcur + i * dx;

        if (j >= 1 && j < Ncur-1) {
            for (int a = 0; a < 3; a++) {
                double dy_phi = (phi[a][IDX(i,j+1,k)] - phi[a][IDX(i,j-1,k)]) / (2.0*dx);
                /* Near oscillon A (x>0): boost in +y direction -> vel += -v_orb * dy_phi */
                /* Near oscillon B (x<0): boost in -y direction -> vel += +v_orb * dy_phi */
                if (x > 0) {
                    vel[a][idx] += -v_orb * dy_phi;
                } else if (x < 0) {
                    vel[a][idx] += +v_orb * dy_phi;
                }
            }
        }
    }

    /* ================================================================
     * STEP 3: Evolve binary with metric, tracking positions and GW
     * ================================================================ */
    printf("\n=== STEP 3: Binary evolution (t=0..%.0f) ===\n", t_binary);
    printf("  Expected orbital freq Omega = %.6f, period = %.1f\n",
           v_orb/(D_sep/2.0), 2.0*M_PI/(v_orb/(D_sep/2.0)));
    printf("  GW freq 2*Omega = %.6f, GW period = %.1f\n",
           2.0*v_orb/(D_sep/2.0), M_PI/(v_orb/(D_sep/2.0)));
    fflush(stdout);

    /* Metric ramp-up over first 50 time units */
    double t_ramp = 50.0;

    /* TT strain measurement angles */
    /* Orbit is in xy-plane -> z-axis is orbital axis (poles) */
    #define N_ANGLES 6
    double gw_theta[N_ANGLES] = {
        0.001,          /* pole (z-axis), theta~0 */
        M_PI/2.0,       /* equator, phi=0 (x-axis) */
        M_PI/2.0,       /* equator, phi=pi/2 (y-axis) */
        M_PI/4.0,       /* 45 deg, phi=0 */
        M_PI - 0.001,   /* south pole */
        M_PI/2.0        /* equator, phi=pi/4 */
    };
    double gw_phi_ang[N_ANGLES] = {
        0.0, 0.0, M_PI/2.0, 0.0, 0.0, M_PI/4.0
    };
    const char *gw_label[N_ANGLES] = {
        "pole(th=0)", "eq_x(th=pi/2,ph=0)", "eq_y(th=pi/2,ph=pi/2)",
        "45deg(th=pi/4,ph=0)", "spole(th=pi)", "eq_45(th=pi/2,ph=pi/4)"
    };
    double R_gw = 20.0;

    int Nt_bin = (int)(t_binary / dt) + 1;

    /* Output files */
    char path[600];
    snprintf(path, sizeof(path), "%s/binary_v%.3f_ts.tsv", outdir, v_orb);
    FILE *fts = fopen(path, "w");
    if (!fts) { fprintf(stderr, "Cannot open %s\n", path); return 1; }
    fprintf(fts, "time\txA\tyA\tzA\txB\tyB\tzB\tsep\tE_total\tEA\tEB\n");

    snprintf(path, sizeof(path), "%s/binary_v%.3f_gw.tsv", outdir, v_orb);
    FILE *fgw = fopen(path, "w");
    if (!fgw) { fprintf(stderr, "Cannot open %s\n", path); return 1; }
    fprintf(fgw, "time");
    for (int g = 0; g < N_ANGLES; g++)
        fprintf(fgw, "\th+_%d\thx_%d", g, g);
    fprintf(fgw, "\n");

    snprintf(path, sizeof(path), "%s/binary_v%.3f_multipoles.tsv", outdir, v_orb);
    FILE *fmp = fopen(path, "w");
    if (!fmp) { fprintf(stderr, "Cannot open %s\n", path); return 1; }
    fprintf(fmp, "time\tc_l0\tc_l1\tc_l2\tfrac_l0\tfrac_l1\tfrac_l2\ttotal_sigma2\n");

    int diag_every = Nt_bin / 2000;
    if (diag_every < 1) diag_every = 1;
    int print_every = Nt_bin / 50;
    if (print_every < 1) print_every = 1;
    int multi_every = Nt_bin / 200;  /* multipoles less frequently */
    if (multi_every < 1) multi_every = 1;

    /* DFT storage for h+ at equator */
    int max_dft = 10000;
    double *hp_eq_hist = malloc(max_dft * sizeof(double));
    double *t_dft_hist = malloc(max_dft * sizeof(double));
    int n_dft = 0;
    int dft_every = Nt_bin / max_dft;
    if (dft_every < 1) dft_every = 1;

    /* Peak strain trackers */
    double max_hp[N_ANGLES], max_hx[N_ANGLES];
    memset(max_hp, 0, sizeof(max_hp));
    memset(max_hx, 0, sizeof(max_hx));

    compute_acc_metric(0.0);

    double wall_bin = omp_get_wtime();

    for (int n = 0; n <= Nt_bin; n++) {
        double t = n * dt;
        double ag_eff = alpha_g * hermite_step(t, 0.0, t_ramp);

        int do_diag  = (n % diag_every == 0);
        int do_print = (n % print_every == 0);
        int do_multi = (n % multi_every == 0);

        if (do_diag || do_print) {
            /* Track positions */
            double xA, yA, zA, EA, xB, yB, zB, EB;
            find_centroid(+1, &xA, &yA, &zA, &EA);
            find_centroid(-1, &xB, &yB, &zB, &EB);
            double sep = sqrt((xA-xB)*(xA-xB) + (yA-yB)*(yA-yB) + (zA-zB)*(zA-zB));
            double Etot = EA + EB;

            if (do_diag) {
                fprintf(fts, "%.4f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.4f\t%.4f\t%.4f\t%.4f\n",
                        t, xA, yA, zA, xB, yB, zB, sep, Etot, EA, EB);

                /* GW strain at all angles */
                fprintf(fgw, "%.4f", t);
                for (int g = 0; g < N_ANGLES; g++) {
                    double hp, hx;
                    compute_tt_strain(gw_theta[g], gw_phi_ang[g], R_gw, &hp, &hx);
                    fprintf(fgw, "\t%.6e\t%.6e", hp, hx);
                    if (fabs(hp) > max_hp[g]) max_hp[g] = fabs(hp);
                    if (fabs(hx) > max_hx[g]) max_hx[g] = fabs(hx);
                }
                fprintf(fgw, "\n");

                /* DFT of equatorial h+ */
                if (n % dft_every == 0 && n_dft < max_dft) {
                    double hp, hx;
                    compute_tt_strain(M_PI/2.0, 0.0, R_gw, &hp, &hx);
                    hp_eq_hist[n_dft] = hp;
                    t_dft_hist[n_dft] = t;
                    n_dft++;
                }
            }

            if (do_multi) {
                double c_l[3], total_s2;
                compute_binary_multipoles(R_gw, c_l, &total_s2);
                double sum_abs = fabs(c_l[0]) + fabs(c_l[1]) + fabs(c_l[2]);
                fprintf(fmp, "%.4f\t%.6e\t%.6e\t%.6e\t%.6f\t%.6f\t%.6f\t%.6e\n",
                        t, c_l[0], c_l[1], c_l[2],
                        (sum_abs > 0) ? fabs(c_l[0])/sum_abs : 0,
                        (sum_abs > 0) ? fabs(c_l[1])/sum_abs : 0,
                        (sum_abs > 0) ? fabs(c_l[2])/sum_abs : 0,
                        total_s2);
            }

            if (do_print) {
                double elapsed = omp_get_wtime() - wall_bin;
                double frac = (double)n / Nt_bin;
                double eta_t = (frac > 0.001) ? elapsed*(1.0-frac)/frac : 0;
                printf("  t=%7.1f  sep=%.2f  posA=(%.2f,%.2f)  E=%.1f  ag=%.4f  "
                       "[%.0fs, ETA %.0fs]\n",
                       t, sep, xA, yA, Etot, ag_eff, elapsed, eta_t);
                fflush(stdout);
            }
        }

        if (n == Nt_bin) break;

        /* Verlet with metric (ramp up) */
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
        double ag_next = alpha_g * hermite_step(t + dt, 0.0, t_ramp);
        compute_acc_metric(ag_next);
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

    fclose(fts);
    fclose(fgw);
    fclose(fmp);

    /* ================================================================
     * ANALYSIS
     * ================================================================ */
    double Omega = v_orb / (D_sep / 2.0);
    double T_orb = 2.0 * M_PI / Omega;

    printf("\n=== ANALYSIS ===\n");
    printf("  Orbital: Omega=%.6f  T_orb=%.1f  2*Omega=%.6f\n", Omega, T_orb, 2.0*Omega);
    printf("  Evolution time: %.0f  (~%.2f orbits)\n", t_binary, t_binary/T_orb);

    printf("\n  Peak TT strain at R=%.0f:\n", R_gw);
    for (int g = 0; g < N_ANGLES; g++) {
        printf("    %s:  |h+|=%.4e  |hx|=%.4e\n",
               gw_label[g], max_hp[g], max_hx[g]);
    }

    /* Angular pattern test */
    printf("\n  SPIN-2 TEST:\n");
    double hp_pole = max_hp[0];
    double hp_eq   = (max_hp[1] > max_hp[2]) ? max_hp[1] : max_hp[2];
    double hp_45   = max_hp[3];

    if (hp_eq > 1e-12) {
        double ratio_pole = hp_pole / hp_eq;
        double ratio_45   = hp_45 / hp_eq;
        printf("    |h+|_pole / |h+|_equator = %.4f", ratio_pole);
        if (ratio_pole < 0.2)
            printf("  -> SUPPRESSED at poles (spin-2 signature)\n");
        else
            printf("  -> NOT suppressed (isotropic or monopolar)\n");

        printf("    |h+|_45deg / |h+|_equator = %.4f", ratio_45);
        if (ratio_45 > 0.2 && ratio_45 < 0.9)
            printf("  -> intermediate (consistent with quadrupole)\n");
        else
            printf("\n");
    } else {
        printf("    h+ too small to test angular pattern\n");
    }

    /* DFT of equatorial h+ to find GW frequency */
    if (n_dft > 100) {
        double peak_om = find_peak_omega(hp_eq_hist, t_dft_hist, n_dft, n_dft/4,
                                         0.2);  /* scan up to omega=0.2 */
        printf("\n  GW frequency from DFT of h+ at equator:\n");
        printf("    omega_GW = %.6f\n", peak_om);
        printf("    2*Omega  = %.6f\n", 2.0*Omega);
        if (peak_om > 0.001) {
            printf("    Ratio omega_GW / (2*Omega) = %.3f", peak_om / (2.0*Omega));
            if (fabs(peak_om/(2.0*Omega) - 1.0) < 0.3)
                printf("  -> MATCHES quadrupole emission at 2*Omega!\n");
            else
                printf("  -> does NOT match 2*Omega\n");
        }
    }

    /* Final multipole snapshot */
    {
        double c_l[3], total_s2;
        compute_binary_multipoles(R_gw, c_l, &total_s2);
        double sum_abs = fabs(c_l[0]) + fabs(c_l[1]) + fabs(c_l[2]);
        printf("\n  Final multipole decomposition (R=%.0f):\n", R_gw);
        for (int l = 0; l <= 2; l++) {
            double frac = (sum_abs > 0) ? fabs(c_l[l])/sum_abs : 0;
            printf("    l=%d: %.6e  (%.1f%%)\n", l, c_l[l], 100*frac);
        }
    }

    double wall_total = omp_get_wtime() - wall_start;
    printf("\n=== DONE (%.1f sec = %.1f min) ===\n", wall_total, wall_total/60);

    /* Cleanup */
    free_fields();
    free(hp_eq_hist);
    free(t_dft_hist);

    return 0;
}
