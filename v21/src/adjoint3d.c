/*
 * adjoint3d.c — Forward-mode tangent optimization for 3D three-body oscillon
 *
 * Optimizes (A, sigma, mu, kappa) to maximize f_core in 3D.
 * Uses forward-mode tangent differentiation through velocity Verlet.
 * All 4 tangent vectors propagated simultaneously (share Jacobian computation).
 *
 * Compile: gcc -O3 -fopenmp -Wall -o adjoint3d v21/src/adjoint3d.c -lm
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

/* Optimization parameters (updated each iteration) */
static double opt_A     = 0.8;
static double opt_sig   = 3.0;
static double opt_mu    = -20.0;
static double opt_kappa = 20.0;

/* Fixed parameters */
static double mass     = 1.0;
static int    N        = 100;
static double L        = 20.0;
static double tfinal   = 200.0;
static double cfl_frac = 0.25;

/* Optimizer settings */
static double lr       = 0.001;
static int    n_iter   = 20;
static int    opt_mu_kappa = 1;  /* 0: only optimize A,sigma; 1: all four */

static void parse_args(int argc, char **argv)
{
    for (int i = 1; i < argc - 1; i += 2) {
        if      (!strcmp(argv[i], "-A"))      opt_A     = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-sigma"))  opt_sig   = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-mu"))     opt_mu    = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-kappa"))  opt_kappa = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-mass"))   mass      = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-N"))      N         = atoi(argv[i+1]);
        else if (!strcmp(argv[i], "-L"))      L         = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-tfinal")) tfinal    = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-cfl"))    cfl_frac  = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-lr"))     lr        = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-niter"))  n_iter    = atoi(argv[i+1]);
        else if (!strcmp(argv[i], "-opt_mk")) opt_mu_kappa = atoi(argv[i+1]);
    }
}

#define IDX(i,j,k) ((long)(i)*N*N + (long)(j)*N + (long)(k))

/* Number of tangent parameters */
#define NPAR 4  /* A, sigma, mu, kappa */

typedef struct {
    double fc;
    double grad[NPAR];  /* dfc/dA, dfc/dsig, dfc/dmu, dfc/dkappa */
} GradResult;

static GradResult compute_gradient(double A, double sig, double mu, double kappa)
{
    double dx  = 2.0 * L / (N - 1);
    double dx2 = dx * dx;
    double m2  = mass * mass;
    double dt  = cfl_frac * dx;
    long   Ngrid = (long)N * N * N;
    int    Nt  = (int)(tfinal / dt) + 1;
    double core_radius = 3.0 * sig;
    double dV = dx * dx * dx;

    /* Allocate field arrays */
    double *phi[3], *vel[3], *acc[3];
    for (int a = 0; a < 3; a++) {
        phi[a] = calloc(Ngrid, sizeof(double));
        vel[a] = calloc(Ngrid, sizeof(double));
        acc[a] = calloc(Ngrid, sizeof(double));
    }

    /* Tangent arrays: dphi[p][a], dvel[p][a], dacc[p][a] for p in {A, sig, mu, kappa} */
    double *dphi[NPAR][3], *dvel[NPAR][3], *dacc[NPAR][3];
    for (int p = 0; p < NPAR; p++)
        for (int a = 0; a < 3; a++) {
            dphi[p][a] = calloc(Ngrid, sizeof(double));
            dvel[p][a] = calloc(Ngrid, sizeof(double));
            dacc[p][a] = calloc(Ngrid, sizeof(double));
        }

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

    /* Initialize fields and tangents */
    #pragma omp parallel for schedule(static)
    for (long idx = 0; idx < Ngrid; idx++) {
        int i = idx / (N * N);
        int j = (idx / N) % N;
        int k = idx % N;
        double x = -L + i * dx;
        double y = -L + j * dx;
        double z = -L + k * dx;
        double r2 = x*x + y*y + z*z;
        double g = exp(-r2 / (2.0 * sig * sig));

        for (int a = 0; a < 3; a++)
            phi[a][idx] = A * g;

        /* dphi/dA = g */
        for (int a = 0; a < 3; a++)
            dphi[0][a][idx] = g;

        /* dphi/dsig = A * r2/sig^3 * g */
        for (int a = 0; a < 3; a++)
            dphi[1][a][idx] = A * r2 / (sig * sig * sig) * g;

        /* dphi/dmu = dphi/dkappa = 0 (initial condition doesn't depend on coupling) */
    }

    /* --- Compute acceleration and tangent acceleration --- */
    /* Force: F_a = -mu * P * dP_a / (1+kappa*P^2)^2 */
    /* where P = phi0*phi1*phi2, dP_a = product of phi_b for b!=a */

    void compute_acc_and_tangent(void) {
        for (int a = 0; a < 3; a++) {
            #pragma omp parallel for schedule(static)
            for (int i = 1; i < N-1; i++) {
                for (int j = 1; j < N-1; j++) {
                    for (int k = 1; k < N-1; k++) {
                        long idx = IDX(i,j,k);

                        /* 7-point Laplacian */
                        double lapl = (phi[a][IDX(i+1,j,k)] + phi[a][IDX(i-1,j,k)]
                                     + phi[a][IDX(i,j+1,k)] + phi[a][IDX(i,j-1,k)]
                                     + phi[a][IDX(i,j,k+1)] + phi[a][IDX(i,j,k-1)]
                                     - 6.0 * phi[a][idx]) / dx2;

                        /* Force computation */
                        double p0 = phi[0][idx], p1 = phi[1][idx], p2 = phi[2][idx];
                        double P = p0 * p1 * p2;
                        double P2 = P * P;
                        double D = 1.0 + kappa * P2;
                        double D2 = D * D;
                        double D3 = D2 * D;

                        double dP_a;
                        switch (a) {
                            case 0: dP_a = p1 * p2; break;
                            case 1: dP_a = p0 * p2; break;
                            default: dP_a = p0 * p1; break;
                        }
                        double F_a = -mu * P * dP_a / D2;

                        acc[a][idx] = lapl - m2 * phi[a][idx] + F_a;

                        /* --- Tangent of acceleration --- */
                        /* dacc_a/dp = lapl(dphi_a/dp) - m2*(dphi_a/dp)
                         *           + sum_b dF_a/dphi_b * dphi_b/dp
                         *           + dF_a/dp_direct                */

                        /* Compute Jacobian dF_a/dphi_b once, reuse for all params */
                        double dP_b[3];
                        dP_b[0] = p1 * p2;
                        dP_b[1] = p0 * p2;
                        dP_b[2] = p0 * p1;

                        double dFdphi[3]; /* dF_a/dphi_b */
                        for (int b = 0; b < 3; b++) {
                            double d2P;
                            if (a == b) d2P = 0.0;
                            else d2P = phi[3-a-b][idx]; /* the third field */

                            double term1 = (dP_b[b] * dP_a + P * d2P) / D2;
                            double term2 = P * dP_a * 4.0 * kappa * P * dP_b[b] / D3;
                            dFdphi[b] = -mu * (term1 - term2);
                        }

                        /* Direct derivatives of F_a w.r.t. mu, kappa */
                        double dF_dmu = -P * dP_a / D2;  /* dF_a/dmu */
                        double dF_dkappa = 2.0 * mu * P * dP_a * P2 / D3; /* dF_a/dkappa */

                        for (int p = 0; p < NPAR; p++) {
                            /* Laplacian of tangent */
                            double dlapl = (dphi[p][a][IDX(i+1,j,k)] + dphi[p][a][IDX(i-1,j,k)]
                                          + dphi[p][a][IDX(i,j+1,k)] + dphi[p][a][IDX(i,j-1,k)]
                                          + dphi[p][a][IDX(i,j,k+1)] + dphi[p][a][IDX(i,j,k-1)]
                                          - 6.0 * dphi[p][a][idx]) / dx2;

                            /* Chain rule through force */
                            double dfp = 0;
                            for (int b = 0; b < 3; b++)
                                dfp += dFdphi[b] * dphi[p][b][idx];

                            /* Direct parameter dependence */
                            if (p == 2) dfp += dF_dmu;     /* p=2 is mu */
                            if (p == 3) dfp += dF_dkappa;  /* p=3 is kappa */

                            dacc[p][a][idx] = dlapl - m2 * dphi[p][a][idx] + dfp;
                        }
                    }
                }
            }

            /* Boundary: zero */
            #pragma omp parallel for schedule(static)
            for (long idx = 0; idx < Ngrid; idx++) {
                int ii = idx / (N * N);
                int jj = (idx / N) % N;
                int kk = idx % N;
                if (ii == 0 || ii == N-1 || jj == 0 || jj == N-1 || kk == 0 || kk == N-1) {
                    acc[a][idx] = 0.0;
                    for (int p = 0; p < NPAR; p++)
                        dacc[p][a][idx] = 0.0;
                }
            }
        }
    }

    compute_acc_and_tangent();

    /* --- Time integration --- */
    for (int n = 0; n < Nt; n++) {
        /* 1. Half-kick */
        for (int a = 0; a < 3; a++) {
            #pragma omp parallel for schedule(static)
            for (long idx = 0; idx < Ngrid; idx++) {
                vel[a][idx] += 0.5 * dt * acc[a][idx];
                for (int p = 0; p < NPAR; p++)
                    dvel[p][a][idx] += 0.5 * dt * dacc[p][a][idx];
            }
        }

        /* 2. Drift */
        for (int a = 0; a < 3; a++) {
            #pragma omp parallel for schedule(static)
            for (long idx = 0; idx < Ngrid; idx++) {
                phi[a][idx] += dt * vel[a][idx];
                for (int p = 0; p < NPAR; p++)
                    dphi[p][a][idx] += dt * dvel[p][a][idx];
            }
        }

        /* 3. Recompute */
        compute_acc_and_tangent();

        /* 4. Half-kick */
        for (int a = 0; a < 3; a++) {
            #pragma omp parallel for schedule(static)
            for (long idx = 0; idx < Ngrid; idx++) {
                vel[a][idx] += 0.5 * dt * acc[a][idx];
                for (int p = 0; p < NPAR; p++)
                    dvel[p][a][idx] += 0.5 * dt * dacc[p][a][idx];
            }
        }

        /* 5. Absorbing boundary */
        for (int a = 0; a < 3; a++) {
            #pragma omp parallel for schedule(static)
            for (long idx = 0; idx < Ngrid; idx++) {
                vel[a][idx] *= damp[idx];
                phi[a][idx] *= damp[idx];
                for (int p = 0; p < NPAR; p++) {
                    dvel[p][a][idx] *= damp[idx];
                    dphi[p][a][idx] *= damp[idx];
                }
            }
        }
    }

    /* --- Compute f_core and its gradient --- */
    double Ecore = 0, Eall = 0;
    double dEcore[NPAR] = {0}, dEall[NPAR] = {0};

    #pragma omp parallel
    {
        double lEc = 0, lEa = 0;
        double ldEc[NPAR] = {0}, ldEa[NPAR] = {0};

        #pragma omp for schedule(static) nowait
        for (int i = 1; i < N-1; i++) {
            double x = -L + i * dx;
            for (int j = 1; j < N-1; j++) {
                double y = -L + j * dx;
                for (int kk = 1; kk < N-1; kk++) {
                    double z = -L + kk * dx;
                    long idx = IDX(i,j,kk);
                    double r = sqrt(x*x + y*y + z*z);

                    /* Energy density */
                    double e = 0;
                    double de[NPAR] = {0};

                    for (int a = 0; a < 3; a++) {
                        /* Kinetic */
                        e += 0.5 * vel[a][idx] * vel[a][idx];
                        for (int p = 0; p < NPAR; p++)
                            de[p] += vel[a][idx] * dvel[p][a][idx];

                        /* Gradient */
                        double gx = (phi[a][IDX(i+1,j,kk)] - phi[a][IDX(i-1,j,kk)]) / (2*dx);
                        double gy = (phi[a][IDX(i,j+1,kk)] - phi[a][IDX(i,j-1,kk)]) / (2*dx);
                        double gz = (phi[a][IDX(i,j,kk+1)] - phi[a][IDX(i,j,kk-1)]) / (2*dx);
                        e += 0.5 * (gx*gx + gy*gy + gz*gz);

                        for (int p = 0; p < NPAR; p++) {
                            double dgx = (dphi[p][a][IDX(i+1,j,kk)] - dphi[p][a][IDX(i-1,j,kk)]) / (2*dx);
                            double dgy = (dphi[p][a][IDX(i,j+1,kk)] - dphi[p][a][IDX(i,j-1,kk)]) / (2*dx);
                            double dgz = (dphi[p][a][IDX(i,j,kk+1)] - dphi[p][a][IDX(i,j,kk-1)]) / (2*dx);
                            de[p] += gx*dgx + gy*dgy + gz*dgz;
                        }

                        /* Mass */
                        e += 0.5 * m2 * phi[a][idx] * phi[a][idx];
                        for (int p = 0; p < NPAR; p++)
                            de[p] += m2 * phi[a][idx] * dphi[p][a][idx];
                    }

                    /* Potential V = mu/2 * P^2 / (1+kappa*P^2) */
                    double p0 = phi[0][idx], p1 = phi[1][idx], p2 = phi[2][idx];
                    double P = p0 * p1 * p2;
                    double P2 = P * P;
                    double D = 1.0 + kappa * P2;
                    double V = 0.5 * mu * P2 / D;
                    e += V;

                    /* dV/dp for each parameter */
                    for (int p = 0; p < NPAR; p++) {
                        /* dP/dp = sum_a prod_{b!=a} phi_b * dphi_a/dp */
                        double dP = 0;
                        double dP_a[3] = {p1*p2, p0*p2, p0*p1};
                        for (int a = 0; a < 3; a++)
                            dP += dP_a[a] * dphi[p][a][idx];

                        /* dV/dP = mu * P / D - mu * kappa * P^3 / D^2
                         *       = mu * P / D^2 */
                        double dV_via_P = mu * P / (D * D) * dP;
                        de[p] += dV_via_P;

                        /* Direct: dV/dmu = P^2/(2D), dV/dkappa = -mu*P^4/(2*D^2) */
                        if (p == 2) de[p] += 0.5 * P2 / D;
                        if (p == 3) de[p] += -0.5 * mu * P2 * P2 / (D * D);
                    }

                    lEa += e * dV;
                    for (int p = 0; p < NPAR; p++)
                        ldEa[p] += de[p] * dV;

                    if (r < core_radius) {
                        lEc += e * dV;
                        for (int p = 0; p < NPAR; p++)
                            ldEc[p] += de[p] * dV;
                    }
                }
            }
        }

        #pragma omp critical
        {
            Ecore += lEc;
            Eall += lEa;
            for (int p = 0; p < NPAR; p++) {
                dEcore[p] += ldEc[p];
                dEall[p] += ldEa[p];
            }
        }
    }

    double fc = (Eall > 1e-20) ? Ecore / Eall : 0.0;
    GradResult res;
    res.fc = fc;
    for (int p = 0; p < NPAR; p++) {
        if (Eall > 1e-20)
            res.grad[p] = (dEcore[p] * Eall - Ecore * dEall[p]) / (Eall * Eall);
        else
            res.grad[p] = 0.0;
    }

    /* Cleanup */
    for (int a = 0; a < 3; a++) {
        free(phi[a]); free(vel[a]); free(acc[a]);
        for (int p = 0; p < NPAR; p++) {
            free(dphi[p][a]); free(dvel[p][a]); free(dacc[p][a]);
        }
    }
    free(damp);

    return res;
}

int main(int argc, char **argv)
{
    parse_args(argc, argv);

    printf("# 3D Adjoint Optimization (forward-mode tangent)\n");
    printf("# mass=%.2f N=%d L=%.1f tfinal=%.0f cfl=%.3f lr=%.4f\n",
           mass, N, L, tfinal, cfl_frac, lr);
    printf("# Starting: A=%.4f sigma=%.4f mu=%.4f kappa=%.4f\n",
           opt_A, opt_sig, opt_mu, opt_kappa);
    printf("# Optimizing: A, sigma%s\n",
           opt_mu_kappa ? ", mu, kappa" : "");
    printf("# iter\tA\tsigma\tmu\tkappa\tfc\tdfc/dA\tdfc/dsig\tdfc/dmu\tdfc/dkap\twall_s\n");
    fflush(stdout);

    for (int iter = 0; iter < n_iter; iter++) {
        /* True vacuum check */
        double m2_min = fabs(opt_mu) * pow(2.0/opt_kappa, 2.0/3.0) / 9.0;
        if (mass * mass < m2_min * 1.1) {
            printf("# WARNING: near false vacuum boundary (m2=%.3f, m2_min=%.3f)\n",
                   mass*mass, m2_min);
        }

        double t0 = omp_get_wtime();
        GradResult g = compute_gradient(opt_A, opt_sig, opt_mu, opt_kappa);
        double wall = omp_get_wtime() - t0;

        printf("%d\t%.6f\t%.6f\t%.4f\t%.4f\t%.6f\t%+.4e\t%+.4e\t%+.4e\t%+.4e\t%.1f\n",
               iter, opt_A, opt_sig, opt_mu, opt_kappa,
               g.fc, g.grad[0], g.grad[1], g.grad[2], g.grad[3], wall);
        fflush(stdout);

        /* Gradient ascent on fc (we want to maximize fc) */
        opt_A   += lr * g.grad[0];
        opt_sig += lr * g.grad[1];
        if (opt_mu_kappa) {
            /* Separate learning rates for physics params (they have different scales) */
            opt_mu    += lr * 0.1 * g.grad[2];
            opt_kappa += lr * 0.1 * g.grad[3];
        }

        /* Clamp */
        if (opt_A < 0.1) opt_A = 0.1;
        if (opt_A > 3.0) opt_A = 3.0;
        if (opt_sig < 1.0) opt_sig = 1.0;
        if (opt_sig > 6.0) opt_sig = 6.0;
        if (opt_mu > -2.0) opt_mu = -2.0;
        if (opt_mu < -100.0) opt_mu = -100.0;
        if (opt_kappa < 1.0) opt_kappa = 1.0;
        if (opt_kappa > 200.0) opt_kappa = 200.0;
    }

    printf("\n# Final: A=%.6f sigma=%.6f mu=%.4f kappa=%.4f\n",
           opt_A, opt_sig, opt_mu, opt_kappa);
    return 0;
}
