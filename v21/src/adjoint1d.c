/*
 * adjoint1d.c — Adjoint optimization for 1D three-body oscillon parameters
 *
 * Forward: velocity Verlet integration of three massive scalars + triple-product
 * Backward: adjoint (reverse-mode autodiff) through the time-stepper
 * Optimizes: A, sigma (initial condition shape) for given (mu, kappa, mass)
 *
 * Loss = -f_core(T) + lambda_omega * max(0, omega/m - 0.9)^2
 *        (maximize core fraction, penalize if omega too close to gap)
 *
 * Gradients: computed via adjoint of velocity Verlet with checkpointing
 *
 * Compile: gcc -O3 -Wall -o adjoint1d v21/src/adjoint1d.c -lm
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Physics parameters (fixed during optimization) */
static double mu     = -10.0;
static double kappa  = 10.0;
static double mass   = 1.0;

/* Optimization variables */
static double A_opt  = 0.8;
static double sig_opt = 3.0;

/* Grid */
static int    Nx     = 1000;
static double xmax   = 40.0;
static double tfinal = 300.0;

/* Optimizer */
static double lr     = 0.01;   /* learning rate */
static int    n_iter = 50;

static void parse_args(int argc, char **argv)
{
    for (int i = 1; i < argc - 1; i += 2) {
        if      (!strcmp(argv[i], "-mu"))     mu      = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-kappa"))  kappa   = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-mass"))   mass    = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-A"))      A_opt   = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-sigma"))  sig_opt = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-Nx"))     Nx      = atoi(argv[i+1]);
        else if (!strcmp(argv[i], "-xmax"))   xmax    = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-tfinal")) tfinal  = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-lr"))     lr      = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-niter"))  n_iter  = atoi(argv[i+1]);
    }
}

/* Force from potential: -dV/dphi_a */
static double force_pot(double p1, double p2, double p3, int a)
{
    double P  = p1 * p2 * p3;
    double P2 = P * P;
    double denom2 = (1.0 + kappa * P2);
    denom2 *= denom2;
    double dP;
    switch (a) {
        case 0: dP = p2 * p3; break;
        case 1: dP = p1 * p3; break;
        case 2: dP = p1 * p2; break;
        default: dP = 0.0;
    }
    return -mu * P * dP / denom2;
}

/* d(force_pot)/d(phi_b) — needed for adjoint */
static void dforce_dphi(double p1, double p2, double p3, int a, double out[3])
{
    double phi[3] = {p1, p2, p3};
    double P = p1 * p2 * p3;
    double P2 = P * P;
    double D = 1.0 + kappa * P2;
    double D2 = D * D;
    double D3 = D2 * D;

    /* dP/dphi_b */
    double dP[3];
    dP[0] = p2 * p3;
    dP[1] = p1 * p3;
    dP[2] = p1 * p2;

    /* force_a = -mu * P * dP_a / D^2 */
    /* d(force_a)/d(phi_b) = -mu * [dP_b * dP_a / D^2  +  P * d(dP_a)/d(phi_b) / D^2
     *                              - 2 * kappa * P * dP_b * P * dP_a / D^3 * 2P] */

    /* Simpler: let F_a = -mu * P * dP_a / D^2 */
    /* dF_a/dphi_b = -mu * [ dP_b * dP_a + P * d2P_ab ] / D^2
     *              + mu * P * dP_a * 2 * kappa * 2 * P * dP_b / D^3 */

    /* d2P_{ab}: d(dP_a)/d(phi_b) */
    /* dP_0 = p2*p3 -> d(dP_0)/dp0 = 0, d(dP_0)/dp1 = p3, d(dP_0)/dp2 = p2 */
    /* In general: d(dP_a)/d(phi_b) = product of phi_c for c != a and c != b */
    /*   = 0 if a == b, else phi_{3-a-b} if all different */

    for (int b = 0; b < 3; b++) {
        double d2P;
        if (a == b) {
            d2P = 0.0;
        } else {
            /* the third index */
            int c = 3 - a - b;
            d2P = phi[c];
        }

        double term1 = (dP[b] * dP[a] + P * d2P) / D2;
        double term2 = P * dP[a] * 4.0 * kappa * P * dP[b] / D3;

        out[b] = -mu * (term1 - term2);
    }
}

typedef struct {
    double loss;
    double fc;
    double omega_over_m;
    double dL_dA;
    double dL_dsig;
} GradResult;

static GradResult compute_gradient(double A, double sig)
{
    double dx = 2.0 * xmax / (Nx - 1);
    double dx2 = dx * dx;
    double m2 = mass * mass;

    double kmax = M_PI / dx;
    double dt = 0.8 * 2.0 / sqrt(kmax*kmax + m2);
    int Nt = (int)(tfinal / dt) + 1;

    double *phi[3], *vel[3], *acc[3];
    for (int a = 0; a < 3; a++) {
        phi[a] = calloc(Nx, sizeof(double));
        vel[a] = calloc(Nx, sizeof(double));
        acc[a] = calloc(Nx, sizeof(double));
    }

    /* Absorbing boundary */
    double *damp = malloc(Nx * sizeof(double));
    double x_abs = xmax * 0.75;
    for (int i = 0; i < Nx; i++) {
        double x = -xmax + i * dx;
        double ax = fabs(x);
        if (ax > x_abs) {
            double f = (ax - x_abs) / (xmax - x_abs);
            damp[i] = 1.0 - 0.98 * f * f;
        } else {
            damp[i] = 1.0;
        }
    }

    /* Initialize: Gaussians */
    /* phi_a(x) = A * exp(-x^2/(2*sig^2)) */
    /* dphi/dA = exp(-x^2/(2*sig^2)) */
    /* dphi/dsig = A * x^2/sig^3 * exp(-x^2/(2*sig^2)) */
    double *dphi_dA[3], *dphi_dsig[3];
    for (int a = 0; a < 3; a++) {
        dphi_dA[a] = calloc(Nx, sizeof(double));
        dphi_dsig[a] = calloc(Nx, sizeof(double));
    }

    for (int a = 0; a < 3; a++)
        for (int i = 0; i < Nx; i++) {
            double x = -xmax + i * dx;
            double g = exp(-x * x / (2.0 * sig * sig));
            phi[a][i] = A * g;
            dphi_dA[a][i] = g;
            dphi_dsig[a][i] = A * x * x / (sig * sig * sig) * g;
        }

    /* Forward-mode tangent propagation (cheaper than adjoint for 2 params) */
    /* Tangent variables: dphi_dA, dvel_dA, dphi_dsig, dvel_dsig */
    double *dvel_dA[3], *dvel_dsig[3];
    for (int a = 0; a < 3; a++) {
        dvel_dA[a] = calloc(Nx, sizeof(double));
        dvel_dsig[a] = calloc(Nx, sizeof(double));
    }

    /* Compute acceleration and its tangent */
    #define COMPUTE_ACC_FWD() do { \
        for (int a = 0; a < 3; a++) { \
            acc[a][0] = acc[a][1] = acc[a][Nx-2] = acc[a][Nx-1] = 0; \
            for (int i = 1; i < Nx - 1; i++) { \
                double lapl = (phi[a][i+1] - 2.0*phi[a][i] + phi[a][i-1]) / dx2; \
                double fp = force_pot(phi[0][i], phi[1][i], phi[2][i], a); \
                acc[a][i] = lapl - m2*phi[a][i] + fp; \
            } \
        } \
    } while(0)

    /* Tangent of acceleration w.r.t. parameter p (A or sigma) */
    /* dacc_a/dp = lapl(dphi_a/dp) - m2*(dphi_a/dp) + sum_b dF_a/dphi_b * dphi_b/dp */
    double *dacc_dA[3], *dacc_dsig[3];
    for (int a = 0; a < 3; a++) {
        dacc_dA[a] = calloc(Nx, sizeof(double));
        dacc_dsig[a] = calloc(Nx, sizeof(double));
    }

    #define COMPUTE_DACC(suffix) do { \
        for (int a = 0; a < 3; a++) { \
            dacc_##suffix[a][0] = dacc_##suffix[a][1] = 0; \
            dacc_##suffix[a][Nx-2] = dacc_##suffix[a][Nx-1] = 0; \
            for (int i = 1; i < Nx - 1; i++) { \
                double dlapl = (dphi_##suffix[a][i+1] - 2.0*dphi_##suffix[a][i] \
                               + dphi_##suffix[a][i-1]) / dx2; \
                double dFdphi[3]; \
                dforce_dphi(phi[0][i], phi[1][i], phi[2][i], a, dFdphi); \
                double dfp = 0; \
                for (int b = 0; b < 3; b++) \
                    dfp += dFdphi[b] * dphi_##suffix[b][i]; \
                dacc_##suffix[a][i] = dlapl - m2 * dphi_##suffix[a][i] + dfp; \
            } \
        } \
    } while(0)

    COMPUTE_ACC_FWD();
    COMPUTE_DACC(dA);
    COMPUTE_DACC(dsig);

    double core_r = 3.0 * sig;

    for (int n = 0; n < Nt; n++) {
        /* Velocity Verlet + tangent */
        /* Step 1: v += 0.5*dt*a */
        for (int a = 0; a < 3; a++)
            for (int i = 1; i < Nx - 1; i++) {
                vel[a][i] += 0.5 * dt * acc[a][i];
                dvel_dA[a][i] += 0.5 * dt * dacc_dA[a][i];
                dvel_dsig[a][i] += 0.5 * dt * dacc_dsig[a][i];
            }
        /* Step 2: x += dt*v */
        for (int a = 0; a < 3; a++)
            for (int i = 1; i < Nx - 1; i++) {
                phi[a][i] += dt * vel[a][i];
                dphi_dA[a][i] += dt * dvel_dA[a][i];
                dphi_dsig[a][i] += dt * dvel_dsig[a][i];
            }
        /* Step 3: recompute a */
        COMPUTE_ACC_FWD();
        COMPUTE_DACC(dA);
        COMPUTE_DACC(dsig);
        /* Step 4: v += 0.5*dt*a */
        for (int a = 0; a < 3; a++)
            for (int i = 1; i < Nx - 1; i++) {
                vel[a][i] += 0.5 * dt * acc[a][i];
                dvel_dA[a][i] += 0.5 * dt * dacc_dA[a][i];
                dvel_dsig[a][i] += 0.5 * dt * dacc_dsig[a][i];
            }
        /* Absorbing boundary (including tangent) */
        for (int a = 0; a < 3; a++)
            for (int i = 0; i < Nx; i++) {
                vel[a][i] *= damp[i];
                phi[a][i] *= damp[i];
                dvel_dA[a][i] *= damp[i];
                dphi_dA[a][i] *= damp[i];
                dvel_dsig[a][i] *= damp[i];
                dphi_dsig[a][i] *= damp[i];
            }
    }

    /* Compute loss and its gradient */
    /* Loss = -f_core = -Ecore/Eall */
    double Ecore = 0, Eall = 0;
    double dEcore_dA = 0, dEall_dA = 0;
    double dEcore_dsig = 0, dEall_dsig = 0;

    for (int i = 1; i < Nx - 1; i++) {
        double x = -xmax + i * dx;
        double e = 0, de_dA = 0, de_dsig = 0;

        for (int a = 0; a < 3; a++) {
            double dp = (phi[a][i+1] - phi[a][i-1]) / (2.0*dx);
            double ddp_dA = (dphi_dA[a][i+1] - dphi_dA[a][i-1]) / (2.0*dx);
            double ddp_dsig = (dphi_dsig[a][i+1] - dphi_dsig[a][i-1]) / (2.0*dx);

            e += 0.5*vel[a][i]*vel[a][i] + 0.5*dp*dp + 0.5*m2*phi[a][i]*phi[a][i];
            de_dA += vel[a][i]*dvel_dA[a][i] + dp*ddp_dA + m2*phi[a][i]*dphi_dA[a][i];
            de_dsig += vel[a][i]*dvel_dsig[a][i] + dp*ddp_dsig + m2*phi[a][i]*dphi_dsig[a][i];
        }

        double P = phi[0][i]*phi[1][i]*phi[2][i];
        double P2 = P*P;
        double V = 0.5*mu*P2/(1.0+kappa*P2);
        e += V;

        /* dV/dp = mu*P*dP/dp / (1+kP^2)^2 */
        double D = 1.0 + kappa*P2;
        for (int pp = 0; pp < 2; pp++) {  /* pp=0: dA, pp=1: dsig */
            double dP = 0;
            for (int a = 0; a < 3; a++) {
                double prod = 1.0;
                for (int b = 0; b < 3; b++) {
                    if (b == a)
                        prod *= (pp == 0) ? dphi_dA[b][i] : dphi_dsig[b][i];
                    else
                        prod *= phi[b][i];
                }
                dP += prod;
            }
            double dV = mu * P * dP / (D*D);
            if (pp == 0) de_dA += dV;
            else de_dsig += dV;
        }

        Eall += e * dx;
        dEall_dA += de_dA * dx;
        dEall_dsig += de_dsig * dx;
        if (fabs(x) < core_r) {
            Ecore += e * dx;
            dEcore_dA += de_dA * dx;
            dEcore_dsig += de_dsig * dx;
        }
    }

    double fc = (Eall > 1e-20) ? Ecore / Eall : 0.0;
    /* d(fc)/dp = (dEcore*Eall - Ecore*dEall) / Eall^2 */
    double dfc_dA = 0, dfc_dsig = 0;
    if (Eall > 1e-20) {
        dfc_dA = (dEcore_dA * Eall - Ecore * dEall_dA) / (Eall * Eall);
        dfc_dsig = (dEcore_dsig * Eall - Ecore * dEall_dsig) / (Eall * Eall);
    }

    /* Loss = -fc */
    double loss = -fc;
    double dL_dA = -dfc_dA;
    double dL_dsig = -dfc_dsig;

    GradResult res = {loss, fc, 0.0, dL_dA, dL_dsig};

    for (int a = 0; a < 3; a++) {
        free(phi[a]); free(vel[a]); free(acc[a]);
        free(dphi_dA[a]); free(dphi_dsig[a]);
        free(dvel_dA[a]); free(dvel_dsig[a]);
        free(dacc_dA[a]); free(dacc_dsig[a]);
    }
    free(damp);

    return res;
}

int main(int argc, char **argv)
{
    parse_args(argc, argv);

    printf("# Adjoint 1D Optimization\n");
    printf("# mu=%.2f kappa=%.2f mass=%.2f\n", mu, kappa, mass);
    printf("# Starting: A=%.4f sigma=%.4f\n", A_opt, sig_opt);
    printf("# iter\tA\tsigma\tloss\tfc\tdL/dA\tdL/dsig\n");

    for (int iter = 0; iter < n_iter; iter++) {
        GradResult g = compute_gradient(A_opt, sig_opt);

        printf("%d\t%.6f\t%.6f\t%.6f\t%.6f\t%+.6e\t%+.6e\n",
               iter, A_opt, sig_opt, g.loss, g.fc, g.dL_dA, g.dL_dsig);
        fflush(stdout);

        /* Gradient descent with clamping */
        A_opt -= lr * g.dL_dA;
        sig_opt -= lr * g.dL_dsig;

        /* Clamp to reasonable range */
        if (A_opt < 0.1) A_opt = 0.1;
        if (A_opt > 3.0) A_opt = 3.0;
        if (sig_opt < 0.5) sig_opt = 0.5;
        if (sig_opt > 6.0) sig_opt = 6.0;
    }

    printf("\n# Final: A=%.6f sigma=%.6f\n", A_opt, sig_opt);
    return 0;
}
