/*
 * deformed.c — V25 Phase 6: Deformed Oscillon — Tidal Spin-2 Radiation
 *
 * Based on v25.c (Phases 1-4). Equilibrates ONCE, saves checkpoint, then
 * runs multiple tidal experiments from that checkpoint.
 *
 * Tidal potential:
 *   V_tidal = -1/2 eps_T * cos(Omega_T * t) * (x^2 - y^2) * Sum_a phi_a^2
 * Force on phi_a:
 *   acc_a += eps_T * cos(Omega_T*t) * (x^2 - y^2) * phi_a
 *
 * Measures: Q_22, Q_20, TT strain multipoles, h+ at equator/poles, Love number k2
 *
 * Compile: gcc -O3 -fopenmp -Wall -o deformed src/deformed.c -lm
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

/* --- Parameters --- */
static double mu_pot    = -20.0;
static double kappa     = 20.0;
static double mass      = 1.0;
static double A_init    = 0.8;
static double sig_init  = 3.0;
static double lambda_pw = 0.5;
static double eta       = 0.1;
static double lambda_L  = 0.1;
static double alpha_g   = 0.001;

/* Tidal parameters */
static double eps_T     = 0.01;
static double Omega_T   = 0.1;

static int    N         = 96;
static double L         = 15.0;
static double t_equil   = 500.0;
static double t_tidal   = 500.0;
static double cfl_frac  = 0.20;
static char   outdir[512] = "data";

static int    scan_mode = 0;  /* 0=single, 1=eps scan, 2=omega scan */

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
        else if (!strcmp(argv[i], "-t_equil"))  t_equil   = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-t_tidal"))  t_tidal   = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-cfl"))      cfl_frac  = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-lpw"))      lambda_pw = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-eta"))      eta       = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-lL"))       lambda_L  = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-ag"))       alpha_g   = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-eps"))      eps_T     = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-omega"))    Omega_T   = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-scan"))     scan_mode = atoi(argv[i+1]);
        else if (!strcmp(argv[i], "-o"))        strncpy(outdir, argv[i+1], sizeof(outdir)-1);
        else { fprintf(stderr, "Unknown flag: %s\n", argv[i]); exit(1); }
    }
}

/* --- Index helpers --- */
#define IDX(i,j,k) ((long)(i)*N*N + (long)(j)*N + (long)(k))

/* --- Globals --- */
static double *phi[3], *vel[3], *acc[3];
static double *phi_ckpt[3], *vel_ckpt[3];  /* checkpoint */
static double *damp;
static double dx, dx2, m2, dt;
static long Ngrid;

/* --- Hermite smooth step --- */
static double hermite_step(double t, double t0, double t1)
{
    if (t <= t0) return 0.0;
    if (t >= t1) return 1.0;
    double s = (t - t0) / (t1 - t0);
    return s * s * (3.0 - 2.0 * s);
}

/* --- Compute acceleration for elastic phase (no tidal) --- */
static void compute_acc_elastic(void)
{
    double eta_lL = eta + lambda_L;

    for (int a = 0; a < 3; a++) {
        #pragma omp parallel for schedule(static)
        for (int i = 2; i < N-2; i++) {
            for (int j = 2; j < N-2; j++) {
                for (int k = 2; k < N-2; k++) {
                    long idx = IDX(i,j,k);

                    double lapl = (phi[a][IDX(i+1,j,k)] + phi[a][IDX(i-1,j,k)]
                                 + phi[a][IDX(i,j+1,k)] + phi[a][IDX(i,j-1,k)]
                                 + phi[a][IDX(i,j,k+1)] + phi[a][IDX(i,j,k-1)]
                                 - 6.0 * phi[a][idx]) / dx2;

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

                    double pw_force = 0.0;
                    for (int b = 0; b < 3; b++) {
                        if (b != a) pw_force += phi[b][idx];
                    }
                    pw_force *= -lambda_pw;

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

                    acc[a][idx] = lapl - m2 * phi[a][idx] - dVdphi
                                + pw_force + eta_lL * d_a_div;
                }
            }
        }

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

/* --- Compute acceleration with tidal potential --- */
static void compute_acc_tidal(double tidal_amp)
{
    double eta_lL = eta + lambda_L;

    for (int a = 0; a < 3; a++) {
        #pragma omp parallel for schedule(static)
        for (int i = 2; i < N-2; i++) {
            double x = -L + i * dx;
            for (int j = 2; j < N-2; j++) {
                double y = -L + j * dx;
                double x2my2 = x*x - y*y;
                for (int k = 2; k < N-2; k++) {
                    long idx = IDX(i,j,k);

                    double lapl = (phi[a][IDX(i+1,j,k)] + phi[a][IDX(i-1,j,k)]
                                 + phi[a][IDX(i,j+1,k)] + phi[a][IDX(i,j-1,k)]
                                 + phi[a][IDX(i,j,k+1)] + phi[a][IDX(i,j,k-1)]
                                 - 6.0 * phi[a][idx]) / dx2;

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

                    double pw_force = 0.0;
                    for (int b = 0; b < 3; b++) {
                        if (b != a) pw_force += phi[b][idx];
                    }
                    pw_force *= -lambda_pw;

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

                    double tidal_force = tidal_amp * x2my2 * phi[a][idx];

                    acc[a][idx] = lapl - m2 * phi[a][idx] - dVdphi
                                + pw_force + eta_lL * d_a_div + tidal_force;
                }
            }
        }

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

/* --- Velocity Verlet steps --- */
static void verlet_step_elastic(void)
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
    compute_acc_elastic();
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

static void verlet_step_tidal(double tidal_amp_next)
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
    compute_acc_tidal(tidal_amp_next);
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

/* --- Energy diagnostics --- */
typedef struct {
    double Ek, Eg, Em, Ep, Epw, Et;
    double peak[3];
    double fc;
} Diag;

static Diag compute_diag(double core_radius)
{
    Diag d;
    memset(&d, 0, sizeof(d));

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

/* --- Quadrupole moments --- */
static double compute_Q22(void)
{
    double Q22 = 0.0;
    #pragma omp parallel for schedule(static) reduction(+:Q22)
    for (int i = 2; i < N-2; i++) {
        double x = -L + i * dx;
        for (int j = 2; j < N-2; j++) {
            double y = -L + j * dx;
            double x2my2 = x*x - y*y;
            for (int k = 2; k < N-2; k++) {
                long idx = IDX(i,j,k);
                double dV = dx * dx * dx;
                double rho = 0.0;
                for (int a = 0; a < 3; a++) {
                    rho += 0.5 * vel[a][idx] * vel[a][idx];
                    double gx = (phi[a][IDX(i+1,j,k)] - phi[a][IDX(i-1,j,k)]) / (2*dx);
                    double gy = (phi[a][IDX(i,j+1,k)] - phi[a][IDX(i,j-1,k)]) / (2*dx);
                    double gz = (phi[a][IDX(i,j,k+1)] - phi[a][IDX(i,j,k-1)]) / (2*dx);
                    rho += 0.5 * (gx*gx + gy*gy + gz*gz);
                    rho += 0.5 * m2 * phi[a][idx] * phi[a][idx];
                }
                double p0 = phi[0][idx], p1 = phi[1][idx], p2 = phi[2][idx];
                double P = p0 * p1 * p2;
                double P2 = P * P;
                rho += 0.5 * mu_pot * P2 / (1.0 + kappa * P2);
                rho += lambda_pw * (p0*p1 + p1*p2 + p2*p0);
                Q22 += x2my2 * rho * dV;
            }
        }
    }
    return Q22;
}

static double compute_Q20(void)
{
    double Q20 = 0.0;
    #pragma omp parallel for schedule(static) reduction(+:Q20)
    for (int i = 2; i < N-2; i++) {
        double x = -L + i * dx;
        for (int j = 2; j < N-2; j++) {
            double y = -L + j * dx;
            for (int k = 2; k < N-2; k++) {
                double z = -L + k * dx;
                long idx = IDX(i,j,k);
                double dV = dx * dx * dx;
                double w = 2.0*z*z - x*x - y*y;
                double rho = 0.0;
                for (int a = 0; a < 3; a++) {
                    rho += 0.5 * vel[a][idx] * vel[a][idx];
                    double gx = (phi[a][IDX(i+1,j,k)] - phi[a][IDX(i-1,j,k)]) / (2*dx);
                    double gy = (phi[a][IDX(i,j+1,k)] - phi[a][IDX(i,j-1,k)]) / (2*dx);
                    double gz = (phi[a][IDX(i,j,k+1)] - phi[a][IDX(i,j,k-1)]) / (2*dx);
                    rho += 0.5 * (gx*gx + gy*gy + gz*gz);
                    rho += 0.5 * m2 * phi[a][idx] * phi[a][idx];
                }
                double p0 = phi[0][idx], p1 = phi[1][idx], p2 = phi[2][idx];
                double P = p0 * p1 * p2;
                double P2 = P * P;
                rho += 0.5 * mu_pot * P2 / (1.0 + kappa * P2);
                rho += lambda_pw * (p0*p1 + p1*p2 + p2*p0);
                Q20 += w * rho * dV;
            }
        }
    }
    return Q20;
}

/* --- Oscillon RMS radius --- */
static double compute_R_rms(void)
{
    double num = 0.0, den = 0.0;
    #pragma omp parallel for schedule(static) reduction(+:num,den)
    for (int i = 2; i < N-2; i++) {
        double x = -L + i * dx;
        for (int j = 2; j < N-2; j++) {
            double y = -L + j * dx;
            for (int k = 2; k < N-2; k++) {
                double z = -L + k * dx;
                long idx = IDX(i,j,k);
                double dV = dx * dx * dx;
                double r2 = x*x + y*y + z*z;
                double rho = 0.0;
                for (int a = 0; a < 3; a++)
                    rho += phi[a][idx] * phi[a][idx];
                num += r2 * rho * dV;
                den += rho * dV;
            }
        }
    }
    return (den > 1e-20) ? sqrt(num / den) : 0.0;
}

/* --- Strain tensor at grid point --- */
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

/* --- TT strain at (theta, phi_ang) on shell R --- */
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

    int i = (int)((x + L) / dx + 0.5);
    int j = (int)((y + L) / dx + 0.5);
    int k = (int)((z + L) / dx + 0.5);

    *h_plus = 0.0;
    *h_cross = 0.0;

    if (i < 2 || i >= N-2 || j < 2 || j >= N-2 || k < 2 || k >= N-2)
        return;

    double eps[6];
    compute_strain(i, j, k, eps);

    double e[3][3];
    e[0][0] = eps[0]; e[1][1] = eps[1]; e[2][2] = eps[2];
    e[0][1] = e[1][0] = eps[3];
    e[0][2] = e[2][0] = eps[4];
    e[1][2] = e[2][1] = eps[5];

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

/* --- Compute strain multipoles on spherical shell --- */
static void compute_strain_multipoles(double R_shell, double c_l[3],
                                       double *total_sigma2)
{
    int N_sample = 200;
    double golden_ratio = (1.0 + sqrt(5.0)) / 2.0;

    c_l[0] = c_l[1] = c_l[2] = 0;
    *total_sigma2 = 0;

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

        for (int l = 0; l <= 2; l++)
            c_l[l] += sigma2 * legendre_P(l, cos_th);
        *total_sigma2 += sigma2;
    }
    *total_sigma2 /= N_sample;
}

/* --- Save checkpoint --- */
static void save_checkpoint(void)
{
    for (int a = 0; a < 3; a++) {
        memcpy(phi_ckpt[a], phi[a], Ngrid * sizeof(double));
        memcpy(vel_ckpt[a], vel[a], Ngrid * sizeof(double));
    }
}

/* --- Restore checkpoint --- */
static void restore_checkpoint(void)
{
    for (int a = 0; a < 3; a++) {
        memcpy(phi[a], phi_ckpt[a], Ngrid * sizeof(double));
        memcpy(vel[a], vel_ckpt[a], Ngrid * sizeof(double));
    }
}

/* ================================================================
 * Run tidal phase from checkpoint (fields already restored)
 * ================================================================ */
typedef struct {
    double eps_T, Omega_T;
    double Q22_max, Q22_rms;
    double Q20_max;
    double hp_equator_max, hp_pole_max;
    double hx_equator_max;
    double hp_45_max;
    double l0_frac, l1_frac, l2_frac;
    double R_osc;
    double k2;
    double E_final;
    double peak_final;
} TidalResult;

static TidalResult run_tidal_phase(double eps, double omega, double R_osc)
{
    TidalResult res;
    memset(&res, 0, sizeof(res));
    res.eps_T = eps;
    res.Omega_T = omega;
    res.R_osc = R_osc;

    double core_radius = 3.0 * sig_init;
    int Nt = (int)(t_tidal / dt) + 1;
    int diag_every = Nt / 2000;
    if (diag_every < 1) diag_every = 1;
    int print_every = Nt / 10;
    if (print_every < 1) print_every = 1;

    /* Output file */
    char path[600];
    if (omega > 1e-10) {
        snprintf(path, sizeof(path), "%s/tidal_driven_eps%.4f_omega%.3f.tsv",
                 outdir, eps, omega);
    } else {
        snprintf(path, sizeof(path), "%s/tidal_static_eps%.4f.tsv", outdir, eps);
    }
    FILE *fts = fopen(path, "w");
    if (!fts) { fprintf(stderr, "Cannot open %s\n", path); return res; }
    fprintf(fts, "time\tQ22\tQ20\thp_equator\thx_equator\thp_pole\thx_pole\t"
                 "hp_45\thx_45\tE_total\tpeak0\n");

    double t_ramp = 50.0;
    double R_gw = 12.0;

    /* Initial acc */
    double amp0 = (omega > 1e-10) ? eps * cos(omega * 0.0) * hermite_step(0.0, 0.0, t_ramp)
                                   : eps * hermite_step(0.0, 0.0, t_ramp);
    compute_acc_tidal(amp0);

    double Q22_max = 0, Q22_sum2 = 0;
    double Q20_max = 0;
    double hp_eq_max = 0, hp_pole_max = 0, hx_eq_max = 0, hp_45_max = 0;
    int n_samples = 0;

    double wall = omp_get_wtime();

    printf("  Tidal phase (eps=%.4f, Omega=%.3f, t=0..%.0f)...\n",
           eps, omega, t_tidal);

    for (int n = 0; n <= Nt; n++) {
        double t = n * dt;

        if (n % diag_every == 0) {
            double Q22 = compute_Q22();
            double Q20 = compute_Q20();

            double hp_eq, hx_eq;
            compute_tt_strain(M_PI/2.0, 0.0, R_gw, &hp_eq, &hx_eq);
            double hp_pole, hx_pole;
            compute_tt_strain(0.001, 0.0, R_gw, &hp_pole, &hx_pole);
            double hp_45, hx_45;
            compute_tt_strain(M_PI/4.0, 0.0, R_gw, &hp_45, &hx_45);

            Diag d = compute_diag(core_radius);

            fprintf(fts, "%.4f\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\n",
                    t, Q22, Q20, hp_eq, hx_eq, hp_pole, hx_pole,
                    hp_45, hx_45, d.Et, d.peak[0]);

            if (fabs(Q22) > Q22_max) Q22_max = fabs(Q22);
            if (fabs(Q20) > Q20_max) Q20_max = fabs(Q20);
            if (fabs(hp_eq) > hp_eq_max) hp_eq_max = fabs(hp_eq);
            if (fabs(hx_eq) > hx_eq_max) hx_eq_max = fabs(hx_eq);
            if (fabs(hp_pole) > hp_pole_max) hp_pole_max = fabs(hp_pole);
            if (fabs(hp_45) > hp_45_max) hp_45_max = fabs(hp_45);

            if (t > t_ramp) {
                Q22_sum2 += Q22 * Q22;
                n_samples++;
            }
        }

        if (n % print_every == 0) {
            double Q22 = compute_Q22();
            Diag d = compute_diag(core_radius);
            double elapsed = omp_get_wtime() - wall;
            double frac = (double)n / Nt;
            double eta_t = (frac > 0.001) ? elapsed*(1.0-frac)/frac : 0;
            printf("    t=%7.1f  Q22=%.4e  pk=%.4f  E=%.2f  [%.0fs, ETA %.0fs]\n",
                   t, Q22, d.peak[0], d.Et, elapsed, eta_t);
            fflush(stdout);
        }

        if (n == Nt) break;

        double t_next = (n + 1) * dt;
        double ramp_next = hermite_step(t_next, 0.0, t_ramp);
        double amp_next;
        if (omega > 1e-10) {
            amp_next = eps * cos(omega * t_next) * ramp_next;
        } else {
            amp_next = eps * ramp_next;
        }
        verlet_step_tidal(amp_next);
    }

    fclose(fts);

    /* Multipole analysis at end */
    double c_l[3], total_sig2;
    compute_strain_multipoles(R_gw, c_l, &total_sig2);
    double sum_abs = fabs(c_l[0]) + fabs(c_l[1]) + fabs(c_l[2]);
    res.l0_frac = (sum_abs > 1e-30) ? fabs(c_l[0]) / sum_abs : 0;
    res.l1_frac = (sum_abs > 1e-30) ? fabs(c_l[1]) / sum_abs : 0;
    res.l2_frac = (sum_abs > 1e-30) ? fabs(c_l[2]) / sum_abs : 0;

    res.Q22_max = Q22_max;
    res.Q22_rms = (n_samples > 0) ? sqrt(Q22_sum2 / n_samples) : 0;
    res.Q20_max = Q20_max;
    res.hp_equator_max = hp_eq_max;
    res.hp_pole_max = hp_pole_max;
    res.hx_equator_max = hx_eq_max;
    res.hp_45_max = hp_45_max;

    double R5 = R_osc * R_osc * R_osc * R_osc * R_osc;
    res.k2 = (eps > 1e-20 && R5 > 1e-20) ? Q22_max / (eps * R5) : 0;

    Diag d = compute_diag(core_radius);
    res.E_final = d.Et;
    res.peak_final = d.peak[0];

    printf("  Done (%.1f s). Q22_max=%.4e  h+_eq=%.4e  h+_pole=%.4e  "
           "l2=%.1f%%  k2=%.4e\n",
           omp_get_wtime() - wall, res.Q22_max, res.hp_equator_max,
           res.hp_pole_max, 100*res.l2_frac, res.k2);
    fflush(stdout);

    return res;
}

/* --- Print result summary --- */
static void print_result(TidalResult *r)
{
    printf("    eps=%.4f Omega=%.3f: Q22_max=%.4e Q22_rms=%.4e Q20_max=%.4e "
           "h+_eq=%.4e h+_pole=%.4e h+_45=%.4e "
           "l0=%.1f%% l1=%.1f%% l2=%.1f%% k2=%.4e E=%.2f\n",
           r->eps_T, r->Omega_T, r->Q22_max, r->Q22_rms, r->Q20_max,
           r->hp_equator_max, r->hp_pole_max, r->hp_45_max,
           100*r->l0_frac, 100*r->l1_frac, 100*r->l2_frac,
           r->k2, r->E_final);
}

/* ================================================================
 * Main
 * ================================================================ */
int main(int argc, char **argv)
{
    parse_args(argc, argv);

    dx  = 2.0 * L / (N - 1);
    dx2 = dx * dx;
    m2  = mass * mass;
    dt  = cfl_frac * dx;
    Ngrid = (long)N * N * N;

    printf("=== V25 Phase 6: Deformed Oscillon — Tidal Spin-2 Radiation ===\n");
    printf("Parameters:\n");
    printf("  mu=%.1f  kappa=%.1f  mass=%.3f  A=%.3f  sigma=%.3f\n",
           mu_pot, kappa, mass, A_init, sig_init);
    printf("  lambda_pw=%.3f  eta=%.3f  lambda_L=%.3f  alpha_g=%.4f\n",
           lambda_pw, eta, lambda_L, alpha_g);
    printf("  eps_T=%.4f  Omega_T=%.3f  scan_mode=%d\n", eps_T, Omega_T, scan_mode);
    printf("  N=%d  L=%.1f  dx=%.4f  dt=%.5f\n", N, L, dx, dt);
    printf("  Ngrid=%ld (%.1f M)  Memory: %.1f MB\n",
           Ngrid, Ngrid/1e6, Ngrid*8.0*10/1e6);
    printf("  t_equil=%.0f  t_tidal=%.0f\n", t_equil, t_tidal);
    printf("  Threads: %d\n", omp_get_max_threads());
    fflush(stdout);

    /* Allocate */
    for (int a = 0; a < 3; a++) {
        phi[a] = calloc(Ngrid, sizeof(double));
        vel[a] = calloc(Ngrid, sizeof(double));
        acc[a] = calloc(Ngrid, sizeof(double));
        phi_ckpt[a] = calloc(Ngrid, sizeof(double));
        vel_ckpt[a] = calloc(Ngrid, sizeof(double));
        if (!phi[a] || !vel[a] || !acc[a] || !phi_ckpt[a] || !vel_ckpt[a]) {
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

    double wall_total = omp_get_wtime();

    /* ================================================================
     * EQUILIBRATION (done once)
     * ================================================================ */
    printf("\n=== EQUILIBRATION (t=0..%.0f) ===\n", t_equil);
    fflush(stdout);

    /* Initialize: spherical Gaussians */
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

    int Nt_eq = (int)(t_equil / dt) + 1;
    int print_every = Nt_eq / 20;
    if (print_every < 1) print_every = 1;

    compute_acc_elastic();

    for (int n = 0; n <= Nt_eq; n++) {
        if (n % print_every == 0) {
            Diag d = compute_diag(core_radius);
            double t = n * dt;
            printf("  t=%7.1f  pk=%.4f  E=%.2f  fc=%.4f\n",
                   t, d.peak[0], d.Et, d.fc);
            fflush(stdout);
        }
        if (n == Nt_eq) break;
        verlet_step_elastic();
    }

    Diag d_eq = compute_diag(core_radius);
    double R_osc = compute_R_rms();
    double Q22_pre = compute_Q22();
    double Q20_pre = compute_Q20();

    printf("Equilibrated: E=%.2f  fc=%.4f  pk=%.4f  R_osc=%.3f\n",
           d_eq.Et, d_eq.fc, d_eq.peak[0], R_osc);
    printf("Pre-tidal: Q22=%.6e  Q20=%.6e  (should be ~0)\n", Q22_pre, Q20_pre);
    printf("Equilibration took %.1f sec\n", omp_get_wtime() - wall_total);
    fflush(stdout);

    /* Save checkpoint */
    save_checkpoint();

    /* ================================================================
     * TIDAL EXPERIMENTS
     * ================================================================ */

    /* Summary table file */
    char spath[600];
    snprintf(spath, sizeof(spath), "%s/tidal_summary.tsv", outdir);
    FILE *fsum = fopen(spath, "w");
    if (!fsum) { fprintf(stderr, "Cannot open %s\n", spath); return 1; }
    fprintf(fsum, "eps_T\tOmega_T\tQ22_max\tQ22_rms\tQ20_max\t"
                  "hp_equator\thp_pole\thx_equator\thp_45\t"
                  "l0_frac\tl1_frac\tl2_frac\tR_osc\tk2\tE_final\n");

    int n_runs = 0;
    TidalResult results[20];

    if (scan_mode == 0) {
        /* Single run */
        restore_checkpoint();
        results[0] = run_tidal_phase(eps_T, Omega_T, R_osc);
        n_runs = 1;

    } else if (scan_mode == 1) {
        /* Phase 6a: Static tidal + Phase 6b: Driven at 3 eps values */
        double eps_vals[] = {0.001, 0.01, 0.1};
        int n_eps = 3;

        printf("\n========================================\n");
        printf("Phase 6a: Static tidal deformation\n");
        printf("========================================\n");
        for (int ie = 0; ie < n_eps; ie++) {
            printf("\n--- Static: eps_T = %.4f ---\n", eps_vals[ie]);
            restore_checkpoint();
            results[n_runs] = run_tidal_phase(eps_vals[ie], 0.0, R_osc);
            print_result(&results[n_runs]);
            n_runs++;
        }

        printf("\n========================================\n");
        printf("Phase 6b: Oscillating tidal (Omega_T=%.3f)\n", Omega_T);
        printf("========================================\n");
        for (int ie = 0; ie < n_eps; ie++) {
            printf("\n--- Driven: eps_T = %.4f, Omega_T = %.3f ---\n",
                   eps_vals[ie], Omega_T);
            restore_checkpoint();
            results[n_runs] = run_tidal_phase(eps_vals[ie], Omega_T, R_osc);
            print_result(&results[n_runs]);
            n_runs++;
        }

    } else if (scan_mode == 2) {
        /* Phase 6c: Omega scan at fixed eps_T */
        double omega_vals[] = {0.05, 0.1, 0.2, 0.5};
        int n_omega = 4;

        printf("\n========================================\n");
        printf("Phase 6c: Resonance scan (eps_T=%.4f)\n", eps_T);
        printf("========================================\n");
        for (int io = 0; io < n_omega; io++) {
            printf("\n--- Driven: eps_T = %.4f, Omega_T = %.3f ---\n",
                   eps_T, omega_vals[io]);
            restore_checkpoint();
            results[n_runs] = run_tidal_phase(eps_T, omega_vals[io], R_osc);
            print_result(&results[n_runs]);
            n_runs++;
        }
    }

    /* Write summary */
    for (int i = 0; i < n_runs; i++) {
        TidalResult *r = &results[i];
        fprintf(fsum, "%.4f\t%.3f\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t"
                      "%.4f\t%.4f\t%.4f\t%.3f\t%.6e\t%.2f\n",
                r->eps_T, r->Omega_T, r->Q22_max, r->Q22_rms, r->Q20_max,
                r->hp_equator_max, r->hp_pole_max, r->hx_equator_max, r->hp_45_max,
                r->l0_frac, r->l1_frac, r->l2_frac, r->R_osc, r->k2, r->E_final);
    }
    fclose(fsum);

    /* ================================================================
     * FINAL SUMMARY
     * ================================================================ */
    double total_time = omp_get_wtime() - wall_total;

    printf("\n=== FINAL SUMMARY ===\n");
    printf("Total wall time: %.1f sec (%.1f min)\n", total_time, total_time/60);
    printf("R_osc = %.3f code units\n", R_osc);
    printf("\n%-8s %-8s %-12s %-12s %-12s %-12s %-8s %-8s %-8s %-12s\n",
           "eps_T", "Omega_T", "Q22_max", "hp_equator", "hp_pole", "hp_45",
           "l0%", "l1%", "l2%", "k2");
    printf("%-8s %-8s %-12s %-12s %-12s %-12s %-8s %-8s %-8s %-12s\n",
           "------", "------", "----------", "----------", "----------", "----------",
           "------", "------", "------", "----------");
    for (int i = 0; i < n_runs; i++) {
        TidalResult *r = &results[i];
        printf("%-8.4f %-8.3f %-12.4e %-12.4e %-12.4e %-12.4e %-8.1f %-8.1f %-8.1f %-12.4e\n",
               r->eps_T, r->Omega_T, r->Q22_max, r->hp_equator_max,
               r->hp_pole_max, r->hp_45_max,
               100*r->l0_frac, 100*r->l1_frac, 100*r->l2_frac, r->k2);
    }
    printf("===\n");

    /* Cleanup */
    for (int a = 0; a < 3; a++) {
        free(phi[a]); free(vel[a]); free(acc[a]);
        free(phi_ckpt[a]); free(vel_ckpt[a]);
    }
    free(damp);

    return 0;
}
