/*
 * scan1d.c — 1D parameter scan for three-body oscillon
 *
 * Runs triad1d physics for moderate time (t=500) across parameter grid.
 * Reports: omega/m ratio, f_core, energy retention, peak amplitude.
 * Goal: find parameters that minimize omega/m (deepest subgap).
 *
 * Compile: gcc -O3 -Wall -o scan1d v21/src/scan1d.c -lm
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    double mu, kappa, mass, A, sigma;
    double omega, omega_over_m, fc, E_retain, peak;
    int valid;
} Result;

static int Nx = 2000;
static double xmax = 60.0;

static double force_pot(double p1, double p2, double p3, int a,
                        double mu, double kappa)
{
    double P = p1 * p2 * p3;
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

static Result run_sim(double mu, double kappa, double mass, double A, double sigma)
{
    Result res = {mu, kappa, mass, A, sigma, 0,0,0,0,0, 0};

    /* True vacuum check */
    double m2_min = fabs(mu) * pow(2.0/kappa, 2.0/3.0) / 9.0;
    if (mass * mass < m2_min) {
        return res; /* false vacuum — skip */
    }

    double dx = 2.0 * xmax / (Nx - 1);
    double dx2 = dx * dx;
    double m2 = mass * mass;
    double tfinal = 500.0;

    double kmax = M_PI / dx;
    double dt = 0.8 * 2.0 / sqrt(kmax*kmax + m2);
    int Nt = (int)(tfinal / dt) + 1;

    double *phi[3], *vel[3], *acc[3];
    for (int a = 0; a < 3; a++) {
        phi[a] = calloc(Nx, sizeof(double));
        vel[a] = calloc(Nx, sizeof(double));
        acc[a] = calloc(Nx, sizeof(double));
    }

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

    for (int a = 0; a < 3; a++)
        for (int i = 0; i < Nx; i++) {
            double x = -xmax + i * dx;
            phi[a][i] = A * exp(-x * x / (2.0 * sigma * sigma));
        }

    #define COMPUTE_ACC_SCAN() do { \
        for (int a = 0; a < 3; a++) { \
            acc[a][0] = acc[a][1] = acc[a][Nx-2] = acc[a][Nx-1] = 0; \
            for (int i = 1; i < Nx - 1; i++) { \
                double lapl = (phi[a][i+1] - 2.0*phi[a][i] + phi[a][i-1]) / dx2; \
                double fp = force_pot(phi[0][i], phi[1][i], phi[2][i], a, mu, kappa); \
                acc[a][i] = lapl - m2*phi[a][i] + fp; \
            } \
        } \
    } while(0)

    COMPUTE_ACC_SCAN();

    /* DFT storage — sample center value in second half */
    int max_dft = 20000;
    double *phi0_hist = malloc(max_dft * sizeof(double));
    double *t_hist = malloc(max_dft * sizeof(double));
    int n_dft = 0;
    int dft_every = Nt / max_dft;
    if (dft_every < 1) dft_every = 1;

    int ic = Nx / 2;
    double core_r = 3.0 * sigma;
    double E0 = 0;

    /* Compute initial energy */
    for (int i = 1; i < Nx - 1; i++) {
        double x = -xmax + i * dx;
        for (int a = 0; a < 3; a++) {
            double dp = (phi[a][i+1] - phi[a][i-1]) / (2.0*dx);
            E0 += (0.5*vel[a][i]*vel[a][i] + 0.5*dp*dp + 0.5*m2*phi[a][i]*phi[a][i]) * dx;
        }
        double P = phi[0][i]*phi[1][i]*phi[2][i];
        double P2 = P*P;
        E0 += 0.5*mu*P2/(1.0+kappa*P2) * dx;
    }

    for (int n = 0; n <= Nt; n++) {
        double t = n * dt;

        if (n % dft_every == 0 && n_dft < max_dft) {
            phi0_hist[n_dft] = phi[0][ic];
            t_hist[n_dft] = t;
            n_dft++;
        }

        if (n == Nt) break;

        for (int a = 0; a < 3; a++)
            for (int i = 1; i < Nx - 1; i++)
                vel[a][i] += 0.5 * dt * acc[a][i];
        for (int a = 0; a < 3; a++)
            for (int i = 1; i < Nx - 1; i++)
                phi[a][i] += dt * vel[a][i];
        COMPUTE_ACC_SCAN();
        for (int a = 0; a < 3; a++)
            for (int i = 1; i < Nx - 1; i++)
                vel[a][i] += 0.5 * dt * acc[a][i];
        for (int a = 0; a < 3; a++)
            for (int i = 0; i < Nx; i++) {
                vel[a][i] *= damp[i];
                phi[a][i] *= damp[i];
            }
    }

    /* Final diagnostics */
    double Ef = 0, Ecore = 0, Eall = 0;
    double peak = 0;
    for (int i = 1; i < Nx - 1; i++) {
        double x = -xmax + i * dx;
        double e = 0;
        for (int a = 0; a < 3; a++) {
            double dp = (phi[a][i+1] - phi[a][i-1]) / (2.0*dx);
            e += 0.5*vel[a][i]*vel[a][i] + 0.5*dp*dp + 0.5*m2*phi[a][i]*phi[a][i];
            if (fabs(phi[a][i]) > peak) peak = fabs(phi[a][i]);
        }
        double P = phi[0][i]*phi[1][i]*phi[2][i];
        double P2 = P*P;
        e += 0.5*mu*P2/(1.0+kappa*P2);
        Ef += e * dx;
        Eall += e * dx;
        if (fabs(x) < core_r) Ecore += e * dx;
    }

    double fc = (Eall > 1e-20) ? Ecore / Eall : 0.0;

    /* DFT — second half only */
    int dft_start = n_dft / 2;
    double peak_omega = 0, peak_pow = 0;
    if (n_dft - dft_start > 50) {
        double T = t_hist[n_dft-1] - t_hist[dft_start];
        int nf = 300;
        for (int k = 1; k < nf; k++) {
            double omega = 3.0 * mass * k / nf;
            double re = 0, im = 0;
            for (int j = dft_start; j < n_dft; j++) {
                double dtj = (j > dft_start) ?
                    (t_hist[j]-t_hist[j-1]) : (t_hist[dft_start+1]-t_hist[dft_start]);
                re += phi0_hist[j] * cos(omega * t_hist[j]) * dtj;
                im += phi0_hist[j] * sin(omega * t_hist[j]) * dtj;
            }
            double pw = (re*re + im*im) / (T*T);
            if (pw > peak_pow) { peak_pow = pw; peak_omega = omega; }
        }
    }

    res.omega = peak_omega;
    res.omega_over_m = (mass > 0) ? peak_omega / mass : 99.0;
    res.fc = fc;
    res.E_retain = (E0 > 1e-20) ? Ef / E0 : 0.0;
    res.peak = peak;
    res.valid = (peak_omega > 0.01 && fc > 0.5);

    for (int a = 0; a < 3; a++) { free(phi[a]); free(vel[a]); free(acc[a]); }
    free(damp); free(phi0_hist); free(t_hist);

    return res;
}

int main(void)
{
    printf("# 1D Parameter Scan: Three-Body Oscillon\n");
    printf("# mu\tkappa\tmass\tA\tsigma\tomega\tw/m\tfc\tE_ret\tpeak\tvalid\n");

    /* Scan grid */
    double mus[]    = {-5, -10, -20, -30};
    double kappas[] = {5, 10, 20, 50};
    double masses[] = {0.8, 1.0, 1.5, 2.0};
    double As[]     = {0.5, 0.8, 1.0, 1.5, 2.0};
    double sigmas[] = {1.5, 2.0, 3.0, 4.0};

    int nmu = sizeof(mus)/sizeof(mus[0]);
    int nk = sizeof(kappas)/sizeof(kappas[0]);
    int nm = sizeof(masses)/sizeof(masses[0]);
    int nA = sizeof(As)/sizeof(As[0]);
    int ns = sizeof(sigmas)/sizeof(sigmas[0]);

    int total = nmu * nk * nm * nA * ns;
    int done = 0;

    for (int im = 0; im < nmu; im++)
    for (int ik = 0; ik < nk; ik++)
    for (int imm = 0; imm < nm; imm++)
    for (int iA = 0; iA < nA; iA++)
    for (int is = 0; is < ns; is++) {
        double mu = mus[im];
        double kappa = kappas[ik];
        double mass = masses[imm];
        double A = As[iA];
        double sig = sigmas[is];

        Result r = run_sim(mu, kappa, mass, A, sig);
        done++;

        if (r.omega > 0.01) {
            printf("%.1f\t%.1f\t%.2f\t%.2f\t%.2f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%d",
                   r.mu, r.kappa, r.mass, r.A, r.sigma,
                   r.omega, r.omega_over_m, r.fc, r.E_retain, r.peak, r.valid);
            if (r.valid && r.omega_over_m < 0.9)
                printf("\t<-- GOOD");
            printf("\n");
            fflush(stdout);
        }

        if (done % 50 == 0)
            fprintf(stderr, "Progress: %d/%d (%.0f%%)\n", done, total, 100.0*done/total);
    }

    fprintf(stderr, "Done: %d configurations\n", total);
    return 0;
}
