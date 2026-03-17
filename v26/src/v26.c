/*
 * v26.c — Braided Solitons: Topological Binding from Knot Theory
 *
 * Phases 1-3:
 *   1. Braided configuration initialization and survival
 *   2. Non-breathing verification (DFT of center density)
 *   3. Strain multipoles on spherical shell
 *
 * Four run configurations:
 *   mode=0: Twisted tube (trefoil braid), full couplings
 *   mode=1: Borromean rings, full couplings
 *   mode=2: Twisted tube, triple product only (m=1, mu, kappa, rest zero)
 *   mode=3: Control oscillon (Gaussian init, full couplings)
 *
 * Lagrangian (identical to V25):
 *   L = Sum_a [1/2(dt phi_a)^2 - 1/2|grad phi_a|^2 - 1/2 m^2 phi_a^2]
 *     - V(P)                               P = phi_1 phi_2 phi_3
 *     - lambda_pw(phi_1 phi_2 + phi_2 phi_3 + phi_3 phi_1)
 *     - 1/2 eta (d_i phi_j)(d_j phi_i)     cross-gradient
 *     - 1/2 lambda_L (div phi)^2            Lame compression
 *
 * Compile: gcc -O3 -fopenmp -Wall -o v26 src/v26.c -lm
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
static double lambda_pw = 0.5;
static double eta       = 0.1;
static double lambda_L  = 0.1;

/* Braid parameters */
static double R_tube    = 3.0;    /* tube radius for twisted tube */
static double R0_ring   = 6.0;    /* ring radius for Borromean */
static double r_tube    = 2.0;    /* tube thickness for Borromean */
static double sig_init  = 3.0;    /* Gaussian sigma for oscillon control */

static int    N         = 128;
static double L         = 20.0;
static double tfinal    = 500.0;
static double cfl_frac  = 0.20;
static int    mode      = 0;      /* 0=twisted, 1=borromean, 2=twisted-tripleonly, 3=oscillon */
static char   outdir[512] = "data";

static void parse_args(int argc, char **argv)
{
    for (int i = 1; i < argc - 1; i += 2) {
        if      (!strcmp(argv[i], "-mu"))       mu_pot    = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-kappa"))    kappa     = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-mass"))     mass      = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-A"))        A0        = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-N"))        N         = atoi(argv[i+1]);
        else if (!strcmp(argv[i], "-L"))        L         = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-tfinal"))   tfinal    = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-cfl"))      cfl_frac  = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-lpw"))      lambda_pw = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-eta"))      eta       = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-lL"))       lambda_L  = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-Rtube"))    R_tube    = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-R0"))       R0_ring   = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-rtube"))    r_tube    = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-mode"))     mode      = atoi(argv[i+1]);
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

/* BC type: 0 = absorbing (spherical damping), 1 = periodic in z */
static int bc_periodic_z = 0;

/* ─── Mode names ─── */
static const char *mode_names[] = {
    "twisted_tube_full",
    "borromean_full",
    "twisted_tube_tripleonly",
    "oscillon_control"
};

/* ─── Initialization: Twisted Tube (Trefoil Braid) ─── */
static void init_twisted_tube(void)
{
    double L_twist = 2.0 * L;  /* one full twist across domain (periodic z) */

    #pragma omp parallel for schedule(static)
    for (long idx = 0; idx < Ngrid; idx++) {
        int ii = idx / (N * N);
        int jj = (idx / N) % N;
        int kk = idx % N;
        double x = -L + ii * dx;
        double y = -L + jj * dx;
        double z = -L + kk * dx;

        double r = sqrt(x * x + y * y);
        double A_r = A0 * exp(-r * r / (2.0 * R_tube * R_tube));

        double twist = 2.0 * M_PI * z / L_twist;

        phi[0][idx] = A_r * cos(twist + 0.0);
        phi[1][idx] = A_r * cos(twist + 2.0 * M_PI / 3.0);
        phi[2][idx] = A_r * cos(twist + 4.0 * M_PI / 3.0);
    }
}

/* ─── Initialization: Borromean Rings ─── */
static void init_borromean(void)
{
    #pragma omp parallel for schedule(static)
    for (long idx = 0; idx < Ngrid; idx++) {
        int ii = idx / (N * N);
        int jj = (idx / N) % N;
        int kk = idx % N;
        double x = -L + ii * dx;
        double y = -L + jj * dx;
        double z = -L + kk * dx;

        /* Ring 1 (phi_0): circle in xy-plane, radius R0 */
        double rho1 = sqrt(x * x + y * y);
        double d1 = sqrt((rho1 - R0_ring) * (rho1 - R0_ring) + z * z);
        phi[0][idx] = A0 * exp(-d1 * d1 / (2.0 * r_tube * r_tube));

        /* Ring 2 (phi_1): circle in xz-plane, radius R0 */
        double rho2 = sqrt(x * x + z * z);
        double d2 = sqrt((rho2 - R0_ring) * (rho2 - R0_ring) + y * y);
        phi[1][idx] = A0 * exp(-d2 * d2 / (2.0 * r_tube * r_tube));

        /* Ring 3 (phi_2): circle in yz-plane, radius R0 */
        double rho3 = sqrt(y * y + z * z);
        double d3 = sqrt((rho3 - R0_ring) * (rho3 - R0_ring) + x * x);
        phi[2][idx] = A0 * exp(-d3 * d3 / (2.0 * r_tube * r_tube));
    }
}

/* ─── Initialization: Oscillon control (spherical Gaussian) ─── */
static void init_oscillon(void)
{
    #pragma omp parallel for schedule(static)
    for (long idx = 0; idx < Ngrid; idx++) {
        int ii = idx / (N * N);
        int jj = (idx / N) % N;
        int kk = idx % N;
        double x = -L + ii * dx;
        double y = -L + jj * dx;
        double z = -L + kk * dx;
        double r2 = x * x + y * y + z * z;
        double g = exp(-r2 / (2.0 * sig_init * sig_init));
        for (int a = 0; a < 3; a++)
            phi[a][idx] = A0 * g;
    }
}

/* ─── Periodic BC helper: wrap index in z direction ─── */
static inline int wrap_z(int k)
{
    if (k < 0)   return k + N;
    if (k >= N)  return k - N;
    return k;
}

/* ─── Neighbor index with optional periodic z BC ─── */
static inline long neighbor(int i, int j, int k)
{
    if (bc_periodic_z) {
        k = wrap_z(k);
    }
    return IDX(i, j, k);
}

/* ─── Compute acceleration ─── */
static void compute_acc(void)
{
    double eta_lL = eta + lambda_L;

    for (int a = 0; a < 3; a++) {
        #pragma omp parallel for schedule(static)
        for (int i = 2; i < N-2; i++) {
            for (int j = 2; j < N-2; j++) {
                /* k range depends on BC */
                int k_lo = bc_periodic_z ? 0 : 2;
                int k_hi = bc_periodic_z ? N : N-2;

                for (int k = k_lo; k < k_hi; k++) {
                    long idx = IDX(i, j, k);

                    /* 7-point Laplacian with periodic z if needed */
                    int km1 = bc_periodic_z ? wrap_z(k-1) : k-1;
                    int kp1 = bc_periodic_z ? wrap_z(k+1) : k+1;

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

                    /* Pairwise coupling: -lambda_pw * (phi_b + phi_c) */
                    double pw_force = 0.0;
                    for (int b = 0; b < 3; b++) {
                        if (b != a) pw_force += phi[b][idx];
                    }
                    pw_force *= -lambda_pw;

                    /* Cross-gradient + Lame: (eta + lambda_L) * d_a(div_phi) */
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
                        int km1b = bc_periodic_z ? wrap_z(k-1) : k-1;
                        int kp1b = bc_periodic_z ? wrap_z(k+1) : k+1;
                        d_a_div += (phi[0][IDX(i+1,j,kp1b)] - phi[0][IDX(i+1,j,km1b)]
                                  - phi[0][IDX(i-1,j,kp1b)] + phi[0][IDX(i-1,j,km1b)])
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
                        int km1b = bc_periodic_z ? wrap_z(k-1) : k-1;
                        int kp1b = bc_periodic_z ? wrap_z(k+1) : k+1;
                        d_a_div += (phi[1][IDX(i,j+1,kp1b)] - phi[1][IDX(i,j+1,km1b)]
                                  - phi[1][IDX(i,j-1,kp1b)] + phi[1][IDX(i,j-1,km1b)])
                                  / (4.0 * dx2);
                    }

                    /* b=2: d_a(d_z phi_2) */
                    if (a == 2) {
                        int km1b = bc_periodic_z ? wrap_z(k-1) : k-1;
                        int kp1b = bc_periodic_z ? wrap_z(k+1) : k+1;
                        d_a_div += (phi[2][IDX(i,j,kp1b)] - 2.0*phi[2][idx]
                                  + phi[2][IDX(i,j,km1b)]) / dx2;
                    } else if (a == 0) {
                        int km1b = bc_periodic_z ? wrap_z(k-1) : k-1;
                        int kp1b = bc_periodic_z ? wrap_z(k+1) : k+1;
                        d_a_div += (phi[2][IDX(i+1,j,kp1b)] - phi[2][IDX(i+1,j,km1b)]
                                  - phi[2][IDX(i-1,j,kp1b)] + phi[2][IDX(i-1,j,km1b)])
                                  / (4.0 * dx2);
                    } else {
                        int km1b = bc_periodic_z ? wrap_z(k-1) : k-1;
                        int kp1b = bc_periodic_z ? wrap_z(k+1) : k+1;
                        d_a_div += (phi[2][IDX(i,j+1,kp1b)] - phi[2][IDX(i,j+1,km1b)]
                                  - phi[2][IDX(i,j-1,kp1b)] + phi[2][IDX(i,j-1,km1b)])
                                  / (4.0 * dx2);
                    }

                    acc[a][idx] = lapl - m2 * phi[a][idx] - dVdphi
                                + pw_force + eta_lL * d_a_div;
                }
            }
        }

        /* Boundary: zero acceleration for non-periodic directions */
        #pragma omp parallel for schedule(static)
        for (long idx = 0; idx < Ngrid; idx++) {
            int ii = idx / (N * N);
            int jj = (idx / N) % N;
            int kk = idx % N;
            if (ii < 2 || ii >= N-2 || jj < 2 || jj >= N-2) {
                acc[a][idx] = 0.0;
            }
            if (!bc_periodic_z && (kk < 2 || kk >= N-2)) {
                acc[a][idx] = 0.0;
            }
        }
    }
}

/* ─── Velocity Verlet step ─── */
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

/* ─── Diagnostics structure ─── */
typedef struct {
    double Ek, Eg, Em, Ep, Epw, Et;
    double peak[3];
    double fc;
    double phi_center[3];
    double peak_P;           /* max |triple product| */
    double rho_center;       /* energy density at center */
} Diag;

/* ─── Compute diagnostics ─── */
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
    double peak_P = 0;
    double rho_ctr = 0;

    /* Compute energy density at center */
    {
        int ci = N/2, cj = N/2, ck = N/2;
        for (int a = 0; a < 3; a++) {
            rho_ctr += 0.5 * vel[a][ic] * vel[a][ic];
            double gx = (phi[a][IDX(ci+1,cj,ck)] - phi[a][IDX(ci-1,cj,ck)]) / (2*dx);
            double gy = (phi[a][IDX(ci,cj+1,ck)] - phi[a][IDX(ci,cj-1,ck)]) / (2*dx);
            double gz = (phi[a][IDX(ci,cj,ck+1)] - phi[a][IDX(ci,cj,ck-1)]) / (2*dx);
            rho_ctr += 0.5 * (gx*gx + gy*gy + gz*gz);
            rho_ctr += 0.5 * m2 * phi[a][ic] * phi[a][ic];
        }
        double p0c = phi[0][ic], p1c = phi[1][ic], p2c = phi[2][ic];
        double Pc = p0c * p1c * p2c;
        double Pc2 = Pc * Pc;
        rho_ctr += 0.5 * mu_pot * Pc2 / (1.0 + kappa * Pc2);
        rho_ctr += lambda_pw * (p0c*p1c + p1c*p2c + p2c*p0c);
    }

    #pragma omp parallel
    {
        double lEk=0, lEg=0, lEm=0, lEp=0, lEpw=0, lEc=0, lEa=0;
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
            Ek += lEk; Eg += lEg; Em += lEm; Ep += lEp; Epw += lEpw;
            Ecore += lEc; Eall += lEa;
            for (int a = 0; a < 3; a++)
                if (lpk[a] > peak[a]) peak[a] = lpk[a];
            if (lpkP > peak_P) peak_P = lpkP;
        }
    }

    d.Ek = Ek; d.Eg = Eg; d.Em = Em; d.Ep = Ep; d.Epw = Epw;
    d.Et = Ek + Eg + Em + Ep + Epw;
    for (int a = 0; a < 3; a++) d.peak[a] = peak[a];
    d.fc = (Eall > 1e-20) ? Ecore / Eall : 0.0;
    d.peak_P = peak_P;
    d.rho_center = rho_ctr;

    return d;
}

/* ─── Fast center density (no full grid scan) ─── */
static double center_rho(void)
{
    long ic = IDX(N/2, N/2, N/2);
    int ci = N/2, cj = N/2, ck = N/2;
    double rho = 0;
    for (int a = 0; a < 3; a++) {
        rho += 0.5 * vel[a][ic] * vel[a][ic];
        double gx = (phi[a][IDX(ci+1,cj,ck)] - phi[a][IDX(ci-1,cj,ck)]) / (2*dx);
        double gy = (phi[a][IDX(ci,cj+1,ck)] - phi[a][IDX(ci,cj-1,ck)]) / (2*dx);
        double gz = (phi[a][IDX(ci,cj,ck+1)] - phi[a][IDX(ci,cj,ck-1)]) / (2*dx);
        rho += 0.5 * (gx*gx + gy*gy + gz*gz);
        rho += 0.5 * m2 * phi[a][ic] * phi[a][ic];
    }
    double p0 = phi[0][ic], p1 = phi[1][ic], p2 = phi[2][ic];
    double P = p0 * p1 * p2;
    double P2 = P * P;
    rho += 0.5 * mu_pot * P2 / (1.0 + kappa * P2);
    rho += lambda_pw * (p0*p1 + p1*p2 + p2*p0);
    return rho;
}

/* ─── DFT to find peak frequency ─── */
static double find_peak_omega(double *hist, double *t_hist, int n_pts, int start,
                              double omega_max, double *peak_power_out)
{
    double T = t_hist[n_pts-1] - t_hist[start];
    if (T < 10.0 || n_pts - start < 50) return 0.0;

    /* Compute mean and subtract */
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

/* ─── Compute strain tensor at a grid point ─── */
static void compute_strain(int i, int j, int k, double eps[6])
{
    double dphi[3][3];
    for (int a = 0; a < 3; a++) {
        dphi[a][0] = (phi[a][IDX(i+1,j,k)] - phi[a][IDX(i-1,j,k)]) / (2.0*dx);
        dphi[a][1] = (phi[a][IDX(i,j+1,k)] - phi[a][IDX(i,j-1,k)]) / (2.0*dx);
        dphi[a][2] = (phi[a][IDX(i,j,k+1)] - phi[a][IDX(i,j,k-1)]) / (2.0*dx);
    }

    eps[0] = dphi[0][0];                              /* eps_xx */
    eps[1] = dphi[1][1];                              /* eps_yy */
    eps[2] = dphi[2][2];                              /* eps_zz */
    eps[3] = 0.5*(dphi[1][0] + dphi[0][1]);           /* eps_xy */
    eps[4] = 0.5*(dphi[2][0] + dphi[0][2]);           /* eps_xz */
    eps[5] = 0.5*(dphi[2][1] + dphi[1][2]);           /* eps_yz */
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

/* ─── Sample strain on spherical shell, compute multipoles ─── */
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

        int i = (int)((x + L) / dx + 0.5);
        int j_idx = (int)((y + L) / dx + 0.5);
        int k = (int)((z + L) / dx + 0.5);

        if (i < 2 || i >= N-2 || j_idx < 2 || j_idx >= N-2 || k < 2 || k >= N-2)
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

/* ─── Run one configuration through all phases ─── */
static void run_config(int run_mode)
{
    const char *name = mode_names[run_mode];
    printf("\n");
    printf("================================================================\n");
    printf("  Configuration: %s (mode=%d)\n", name, run_mode);
    printf("================================================================\n\n");

    /* Set coupling parameters based on mode */
    double save_mass = mass, save_lpw = lambda_pw, save_eta = eta, save_lL = lambda_L;

    if (run_mode == 2) {
        /* Triple product only: keep m, mu, kappa; zero everything else */
        lambda_pw = 0.0;
        eta = 0.0;
        lambda_L = 0.0;
    }
    m2 = mass * mass;

    /* Set BC type */
    if (run_mode == 0 || run_mode == 2) {
        bc_periodic_z = 1;  /* twisted tube needs periodic z */
    } else {
        bc_periodic_z = 0;  /* absorbing for Borromean and oscillon */
    }

    /* Zero fields and velocities */
    for (int a = 0; a < 3; a++) {
        memset(phi[a], 0, Ngrid * sizeof(double));
        memset(vel[a], 0, Ngrid * sizeof(double));
        memset(acc[a], 0, Ngrid * sizeof(double));
    }

    /* Initialize fields */
    switch (run_mode) {
        case 0: case 2: init_twisted_tube(); break;
        case 1:         init_borromean(); break;
        case 3:         init_oscillon(); break;
    }

    /* Set damping profile */
    if (bc_periodic_z) {
        /* For periodic z: only damp in x,y at large radius from z-axis */
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
    } else {
        /* Spherical absorbing shell */
        double R_abs_inner = L * 0.70;
        double R_abs_outer = L * 0.95;
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
    }

    double core_radius = 8.0;  /* for core fraction */

    /* ================================================================
     * PHASE 1: Braided Configuration Survival
     * ================================================================ */
    printf("--- Phase 1: Evolution (t=0..%.0f) ---\n", tfinal);
    fflush(stdout);

    int Nt = (int)(tfinal / dt) + 1;

    /* DFT storage for Phase 2 */
    int max_dft = 50000;
    double *rho_hist = malloc(max_dft * sizeof(double));  /* rho(center, t) */
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
    snprintf(path, sizeof(path), "%s/v26_%s_phase1.tsv", outdir, name);
    FILE *f1 = fopen(path, "w");
    if (!f1) { fprintf(stderr, "Cannot open %s\n", path); return; }
    fprintf(f1, "time\tE_total\tE_kin\tE_grad\tE_mass\tE_pot\tE_pw\t"
                "fc\tpeak0\tpeak1\tpeak2\tpeak_P\trho_center\n");

    compute_acc();

    double wall_start = omp_get_wtime();

    /* Initial diagnostics */
    {
        Diag d0 = compute_diag(core_radius);
        printf("  t=%7.1f  E=%.2f  fc=%.4f  pk=(%.4f,%.4f,%.4f)  |P|max=%.6f  rho_ctr=%.4e\n",
               0.0, d0.Et, d0.fc, d0.peak[0], d0.peak[1], d0.peak[2], d0.peak_P, d0.rho_center);
    }

    for (int n = 0; n <= Nt; n++) {
        double t = n * dt;

        /* DFT history (fast, center only) */
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
                fprintf(f1, "%.4f\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t"
                            "%.6f\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\n",
                        t, d.Et, d.Ek, d.Eg, d.Em, d.Ep, d.Epw,
                        d.fc, d.peak[0], d.peak[1], d.peak[2], d.peak_P, d.rho_center);

            if (do_print) {
                double elapsed = omp_get_wtime() - wall_start;
                double frac = (double)n / Nt;
                double eta_t = (frac > 0.001) ? elapsed * (1.0-frac)/frac : 0;
                printf("  t=%7.1f  E=%.2f  fc=%.4f  pk=(%.4f,%.4f,%.4f)  |P|max=%.6f  [%.0fs, ETA %.0fs]\n",
                       t, d.Et, d.fc, d.peak[0], d.peak[1], d.peak[2], d.peak_P,
                       elapsed, eta_t);
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

    printf("\nPhase 1 complete (%.1f sec)\n", elapsed1);
    printf("  Final: E=%.2f  fc=%.4f  pk=(%.4f,%.4f,%.4f)  |P|max=%.6f\n",
           dfinal.Et, dfinal.fc, dfinal.peak[0], dfinal.peak[1], dfinal.peak[2], dfinal.peak_P);

    /* Determine survival */
    int survived = (dfinal.fc > 0.01 && dfinal.peak_P > 1e-6);
    printf("  Survived? %s (fc=%.4f, |P|max=%.6f)\n",
           survived ? "YES" : "NO", dfinal.fc, dfinal.peak_P);

    /* ================================================================
     * PHASE 2: Non-Breathing Verification
     * ================================================================ */
    printf("\n--- Phase 2: Non-Breathing Verification ---\n");

    if (n_dft < 100) {
        printf("  Not enough DFT samples (%d), skipping.\n", n_dft);
    } else {
        int start = n_dft / 2;  /* second half only */

        /* DFT of rho(center, t) */
        double peak_power_rho = 0;
        double omega_rho = find_peak_omega(rho_hist, t_hist, n_dft, start, 5.0, &peak_power_rho);

        /* DFT of phi_0(center, t) */
        double peak_power_phi = 0;
        double omega_phi = find_peak_omega(phi0_hist, t_hist, n_dft, start, 5.0, &peak_power_phi);

        /* Compute variance for comparison */
        double mean_rho = 0;
        for (int j = start; j < n_dft; j++) mean_rho += rho_hist[j];
        mean_rho /= (n_dft - start);

        /* Variance of rho */
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
        snprintf(path, sizeof(path), "%s/v26_%s_phase2_dft.tsv", outdir, name);
        FILE *fdft = fopen(path, "w");
        if (fdft) {
            fprintf(fdft, "# DFT of rho(center,t) second half\n");
            fprintf(fdft, "omega\tpower_rho\tpower_phi\n");
            /* Precompute mean_phi */
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
        snprintf(path, sizeof(path), "%s/v26_%s_rho_history.tsv", outdir, name);
        FILE *frho = fopen(path, "w");
        if (frho) {
            fprintf(frho, "time\trho_center\tphi0_center\n");
            for (int j = 0; j < n_dft; j++)
                fprintf(frho, "%.4f\t%.6e\t%.6e\n", t_hist[j], rho_hist[j], phi0_hist[j]);
            fclose(frho);
        }
    }

    /* ================================================================
     * PHASE 3: Strain Multipoles
     * ================================================================ */
    printf("\n--- Phase 3: Strain Multipoles ---\n");

    double R_shell = 8.0;

    snprintf(path, sizeof(path), "%s/v26_%s_phase3_strain.tsv", outdir, name);
    FILE *fstrain = fopen(path, "w");
    snprintf(path, sizeof(path), "%s/v26_%s_phase3_multipoles.tsv", outdir, name);
    FILE *fmulti = fopen(path, "w");

    if (!fstrain || !fmulti) {
        fprintf(stderr, "Cannot open Phase 3 output files\n");
    } else {
        compute_strain_multipoles(R_shell, fstrain, fmulti);
        fclose(fstrain);
        fclose(fmulti);
    }

    printf("Phase 3 complete.\n");
    fflush(stdout);

    /* Restore parameters */
    mass = save_mass;
    lambda_pw = save_lpw;
    eta = save_eta;
    lambda_L = save_lL;
    m2 = mass * mass;

    free(rho_hist);
    free(phi0_hist);
    free(t_hist);
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

    printf("=== V26: Braided Solitons — Topological Binding from Knot Theory ===\n");
    printf("Parameters:\n");
    printf("  mu=%.1f  kappa=%.1f  mass=%.3f  A0=%.3f\n",
           mu_pot, kappa, mass, A0);
    printf("  lambda_pw=%.3f  eta=%.3f  lambda_L=%.3f\n",
           lambda_pw, eta, lambda_L);
    printf("  R_tube=%.1f  R0_ring=%.1f  r_tube=%.1f\n", R_tube, R0_ring, r_tube);
    printf("  N=%d  L=%.1f  dx=%.4f  dt=%.5f\n", N, L, dx, dt);
    printf("  Ngrid=%ld (%.1f M)  Memory: %.1f MB\n",
           Ngrid, Ngrid/1e6, Ngrid*8.0*10/1e6);
    printf("  Threads: %d\n", omp_get_max_threads());
    printf("  mode=%d  tfinal=%.0f\n", mode, tfinal);
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

    if (mode >= 0 && mode <= 3) {
        /* Run single configuration */
        run_config(mode);
    } else if (mode == 99) {
        /* Run all four configurations sequentially */
        for (int m = 0; m <= 3; m++) {
            run_config(m);
        }
    } else {
        fprintf(stderr, "Unknown mode %d. Use 0-3 or 99 (all).\n", mode);
        return 1;
    }

    /* Cleanup */
    for (int a = 0; a < 3; a++) {
        free(phi[a]); free(vel[a]); free(acc[a]);
    }
    free(damp);

    printf("\n=== V26 Complete ===\n");
    return 0;
}
