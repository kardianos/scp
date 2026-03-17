/*
 * eos.c — Equation of state diagnostic for the three-body oscillon
 *
 * Computes full stress-energy tensor T^{mu nu} at every grid point:
 *   T^00 = energy density rho = sum_a [v_a^2/2 + (dx phi_a)^2/2 + m^2 phi_a^2/2] + V
 *   T^11 = pressure P = sum_a [v_a^2/2 + (dx phi_a)^2/2 - m^2 phi_a^2/2] - V
 *   T^01 = momentum flux = -sum_a v_a * dx phi_a
 *
 * Time-averages over breathing cycles after equilibration.
 * Outputs P vs rho at each grid point, trace T^mu_mu, P/rho ratios.
 *
 * Two modes: (1) triple-product only, (2) with pairwise coupling lambda
 *
 * Compile: gcc -O3 -Wall -o eos src/eos.c -lm
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Parameters */
static double mu       = -20.0;
static double kappa    = 20.0;
static double mass     = 1.0;
static double A_init   = 0.8;
static double sigma_init = 3.0;
static int    Nx       = 4000;
static double xmax     = 100.0;
static double tfinal   = 10000.0;
static double lambda   = 0.0;    /* pairwise coupling */
static char   outdir[512] = "v24/fundamental/testA_eos/data";

static void parse_args(int argc, char **argv)
{
    for (int i = 1; i < argc - 1; i += 2) {
        if      (!strcmp(argv[i], "-mu"))     mu       = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-kappa"))  kappa    = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-mass"))   mass     = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-A"))      A_init   = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-sigma"))  sigma_init = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-Nx"))     Nx       = atoi(argv[i+1]);
        else if (!strcmp(argv[i], "-xmax"))   xmax     = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-tfinal")) tfinal   = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-lambda")) lambda   = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-o"))      strncpy(outdir, argv[i+1], sizeof(outdir)-1);
        else { fprintf(stderr, "Unknown flag: %s\n", argv[i]); exit(1); }
    }
}

/* Triple-product potential: V = (mu/2) P^2 / (1 + kappa P^2) */
static inline double pot_triple(double p1, double p2, double p3)
{
    double P  = p1 * p2 * p3;
    double P2 = P * P;
    return 0.5 * mu * P2 / (1.0 + kappa * P2);
}

/* Pairwise potential: V_pw = lambda * (p1*p2 + p2*p3 + p3*p1) */
static inline double pot_pairwise(double p1, double p2, double p3)
{
    return lambda * (p1*p2 + p2*p3 + p3*p1);
}

/* -dV/dphi_a for triple product */
static double force_triple(double p1, double p2, double p3, int a)
{
    double P  = p1 * p2 * p3;
    double P2 = P * P;
    double denom2 = (1.0 + kappa * P2) * (1.0 + kappa * P2);
    double dP;
    switch (a) {
        case 0: dP = p2 * p3; break;
        case 1: dP = p1 * p3; break;
        case 2: dP = p1 * p2; break;
        default: dP = 0.0;
    }
    return -mu * P * dP / denom2;
}

/* -dV_pw/dphi_a for pairwise coupling */
static double force_pairwise(double p1, double p2, double p3, int a)
{
    switch (a) {
        case 0: return -lambda * (p2 + p3);
        case 1: return -lambda * (p1 + p3);
        case 2: return -lambda * (p1 + p2);
        default: return 0.0;
    }
}

int main(int argc, char **argv)
{
    parse_args(argc, argv);

    double dx  = 2.0 * xmax / (Nx - 1);
    double dx2 = dx * dx;
    double m2  = mass * mass;

    /* CFL condition */
    double kmax = M_PI / dx;
    double dt = 0.8 * 2.0 / sqrt(kmax * kmax + m2);
    int Nt = (int)(tfinal / dt) + 1;

    printf("eos: Stress-energy tensor diagnostic\n");
    printf("  mu=%.1f kappa=%.1f mass=%.3f lambda=%.4f\n", mu, kappa, mass, lambda);
    printf("  Nx=%d xmax=%.0f dx=%.5f dt=%.6f Nt=%d tfinal=%.0f\n",
           Nx, xmax, dx, dt, Nt, tfinal);

    /* Allocate fields */
    double *phi[3], *vel[3], *acc[3];
    for (int a = 0; a < 3; a++) {
        phi[a] = calloc(Nx, sizeof(double));
        vel[a] = calloc(Nx, sizeof(double));
        acc[a] = calloc(Nx, sizeof(double));
    }

    /* Absorbing boundary: outer 25% */
    double *damp = malloc(Nx * sizeof(double));
    double x_abs = xmax * 0.75;
    for (int i = 0; i < Nx; i++) {
        double x = -xmax + i * dx;
        double ax = fabs(x);
        if (ax > x_abs)
            damp[i] = 1.0 - 0.98 * ((ax - x_abs) / (xmax - x_abs)) *
                                     ((ax - x_abs) / (xmax - x_abs));
        else
            damp[i] = 1.0;
    }

    /* Initialize: symmetric Gaussians */
    for (int a = 0; a < 3; a++)
        for (int i = 0; i < Nx; i++) {
            double x = -xmax + i * dx;
            phi[a][i] = A_init * exp(-x * x / (2.0 * sigma_init * sigma_init));
        }

    /* Compute acceleration */
    #define COMPUTE_ACC() do { \
        for (int a = 0; a < 3; a++) { \
            acc[a][0] = acc[a][1] = acc[a][Nx-2] = acc[a][Nx-1] = 0; \
            for (int i = 1; i < Nx - 1; i++) { \
                double lapl = (phi[a][i+1] - 2.0*phi[a][i] + phi[a][i-1]) / dx2; \
                double ft = force_triple(phi[0][i], phi[1][i], phi[2][i], a); \
                double fp = (lambda != 0.0) ? force_pairwise(phi[0][i], phi[1][i], phi[2][i], a) : 0.0; \
                acc[a][i] = lapl - m2 * phi[a][i] + ft + fp; \
            } \
        } \
    } while(0)

    COMPUTE_ACC();

    /* Time-averaging accumulators (only in the inner non-absorbing region) */
    /* We average from t_avg_start to tfinal */
    double t_avg_start = 5000.0;
    int n_avg_start = (int)(t_avg_start / dt);
    int n_avg = 0;

    double *avg_T00 = calloc(Nx, sizeof(double));
    double *avg_T11 = calloc(Nx, sizeof(double));
    double *avg_T01 = calloc(Nx, sizeof(double));
    double *avg_T00sq = calloc(Nx, sizeof(double));  /* for variance */
    double *avg_T11sq = calloc(Nx, sizeof(double));

    /* Also track peak phi at center for period measurement */
    int max_dft = 50000;
    double *phi0_hist = malloc(max_dft * sizeof(double));
    double *t_hist = malloc(max_dft * sizeof(double));
    int n_dft = 0;
    int dft_every = Nt / max_dft;
    if (dft_every < 1) dft_every = 1;

    int ic = Nx / 2;
    int print_every = Nt / 40;
    if (print_every < 1) print_every = 1;

    /* How often to accumulate averages — every 10 steps is fine */
    int avg_every = 10;

    for (int n = 0; n <= Nt; n++) {
        double t = n * dt;

        /* DFT recording */
        if (n % dft_every == 0 && n_dft < max_dft) {
            phi0_hist[n_dft] = phi[0][ic];
            t_hist[n_dft] = t;
            n_dft++;
        }

        /* Accumulate stress-energy averages */
        if (n >= n_avg_start && n % avg_every == 0) {
            for (int i = 1; i < Nx - 1; i++) {
                double kin = 0, grad = 0, mass_e = 0;
                double mom = 0;
                for (int a = 0; a < 3; a++) {
                    double v = vel[a][i];
                    double dp = (phi[a][i+1] - phi[a][i-1]) / (2.0 * dx);
                    kin += 0.5 * v * v;
                    grad += 0.5 * dp * dp;
                    mass_e += 0.5 * m2 * phi[a][i] * phi[a][i];
                    mom += v * dp;
                }
                double Vt = pot_triple(phi[0][i], phi[1][i], phi[2][i]);
                double Vp = (lambda != 0.0) ? pot_pairwise(phi[0][i], phi[1][i], phi[2][i]) : 0.0;
                double V = Vt + Vp;

                double T00 = kin + grad + mass_e + V;
                double T11 = kin + grad - mass_e - V;
                double T01 = -mom;

                avg_T00[i] += T00;
                avg_T11[i] += T11;
                avg_T01[i] += T01;
                avg_T00sq[i] += T00 * T00;
                avg_T11sq[i] += T11 * T11;
            }
            n_avg++;
        }

        /* Progress */
        if (n % print_every == 0) {
            double Et = 0;
            for (int i = 1; i < Nx - 1; i++) {
                for (int a = 0; a < 3; a++) {
                    Et += (0.5 * vel[a][i] * vel[a][i] +
                           0.5 * m2 * phi[a][i] * phi[a][i]) * dx;
                    double dp = (phi[a][i+1] - phi[a][i-1]) / (2.0 * dx);
                    Et += 0.5 * dp * dp * dx;
                }
                Et += pot_triple(phi[0][i], phi[1][i], phi[2][i]) * dx;
                if (lambda != 0.0)
                    Et += pot_pairwise(phi[0][i], phi[1][i], phi[2][i]) * dx;
            }
            printf("  t=%7.1f  phi0=(%+.4f,%+.4f,%+.4f)  E=%.4f  n_avg=%d\n",
                   t, phi[0][ic], phi[1][ic], phi[2][ic], Et, n_avg);
        }

        if (n == Nt) break;

        /* Velocity Verlet */
        for (int a = 0; a < 3; a++)
            for (int i = 1; i < Nx - 1; i++)
                vel[a][i] += 0.5 * dt * acc[a][i];
        for (int a = 0; a < 3; a++)
            for (int i = 1; i < Nx - 1; i++)
                phi[a][i] += dt * vel[a][i];
        COMPUTE_ACC();
        for (int a = 0; a < 3; a++)
            for (int i = 1; i < Nx - 1; i++)
                vel[a][i] += 0.5 * dt * acc[a][i];

        /* Absorbing boundary */
        for (int a = 0; a < 3; a++)
            for (int i = 0; i < Nx; i++) {
                vel[a][i] *= damp[i];
                phi[a][i] *= damp[i];
            }
    }

    printf("\nAveraging: %d samples from t=%.0f to t=%.0f\n",
           n_avg, t_avg_start, tfinal);

    /* Normalize averages */
    if (n_avg > 0) {
        for (int i = 0; i < Nx; i++) {
            avg_T00[i] /= n_avg;
            avg_T11[i] /= n_avg;
            avg_T01[i] /= n_avg;
            avg_T00sq[i] /= n_avg;
            avg_T11sq[i] /= n_avg;
        }
    }

    /* Write spatial profile of time-averaged T^{mu nu} */
    char path[600];
    const char *suffix = (lambda > 0) ? "_pw" : "";
    snprintf(path, sizeof(path), "%s/eos_profile%s.tsv", outdir, suffix);
    FILE *fp = fopen(path, "w");
    if (!fp) { fprintf(stderr, "Cannot open %s\n", path); return 1; }
    fprintf(fp, "x\trho\tP\tT01\tP_over_rho\ttrace\trho_std\tP_std\n");
    for (int i = 1; i < Nx - 1; i++) {
        double x = -xmax + i * dx;
        double rho = avg_T00[i];
        double P   = avg_T11[i];
        double T01 = avg_T01[i];
        double P_rho = (fabs(rho) > 1e-20) ? P / rho : 0.0;
        /* Trace: T^mu_mu with (+,-) signature = T^0_0 + T^1_1 = rho - P */
        /* But proposal says "T^mu_mu = -rho + P" (using (-,+) convention) */
        /* We output both: trace_plus = rho - P, trace_minus = -rho + P */
        double trace = rho - P;  /* = 2*(mass_energy + V) */
        double rho_std = sqrt(fabs(avg_T00sq[i] - rho * rho));
        double P_std   = sqrt(fabs(avg_T11sq[i] - P * P));
        fprintf(fp, "%.6f\t%.8e\t%.8e\t%.8e\t%.6f\t%.8e\t%.8e\t%.8e\n",
                x, rho, P, T01, P_rho, trace, rho_std, P_std);
    }
    fclose(fp);
    printf("Wrote: %s\n", path);

    /* Write P-vs-rho scatter data (one point per grid point in inner region) */
    snprintf(path, sizeof(path), "%s/eos_scatter%s.tsv", outdir, suffix);
    fp = fopen(path, "w");
    fprintf(fp, "rho\tP\tP_over_rho\ttrace\tx\n");
    for (int i = 1; i < Nx - 1; i++) {
        double x = -xmax + i * dx;
        if (fabs(x) > x_abs) continue;  /* skip absorbing layer */
        double rho = avg_T00[i];
        double P   = avg_T11[i];
        if (fabs(rho) < 1e-30 && fabs(P) < 1e-30) continue;
        double P_rho = (fabs(rho) > 1e-20) ? P / rho : 0.0;
        double trace = rho - P;
        fprintf(fp, "%.8e\t%.8e\t%.6f\t%.8e\t%.4f\n", rho, P, P_rho, trace, x);
    }
    fclose(fp);
    printf("Wrote: %s\n", path);

    /* Summary statistics */
    /* Core: |x| < 5 */
    double core_rho = 0, core_P = 0, core_T01 = 0, core_trace = 0;
    int core_n = 0;
    /* Vacuum: 30 < |x| < 60 (well outside core, inside non-absorbing) */
    double vac_rho = 0, vac_P = 0, vac_T01 = 0, vac_trace = 0;
    int vac_n = 0;
    /* Mid: 10 < |x| < 20 */
    double mid_rho = 0, mid_P = 0;
    int mid_n = 0;

    for (int i = 1; i < Nx - 1; i++) {
        double x = -xmax + i * dx;
        double ax = fabs(x);
        double rho = avg_T00[i];
        double P   = avg_T11[i];
        double T01 = avg_T01[i];
        double tr  = rho - P;

        if (ax < 5.0) {
            core_rho += rho; core_P += P; core_T01 += T01; core_trace += tr;
            core_n++;
        }
        if (ax > 30.0 && ax < 60.0) {
            vac_rho += rho; vac_P += P; vac_T01 += T01; vac_trace += tr;
            vac_n++;
        }
        if (ax > 10.0 && ax < 20.0) {
            mid_rho += rho; mid_P += P;
            mid_n++;
        }
    }

    printf("\n=== EQUATION OF STATE RESULTS (lambda=%.4f) ===\n", lambda);

    if (core_n > 0) {
        core_rho /= core_n; core_P /= core_n; core_T01 /= core_n; core_trace /= core_n;
        printf("\nCORE (|x|<5): <rho>=%.6e  <P>=%.6e  P/rho=%.6f\n",
               core_rho, core_P, (fabs(core_rho) > 1e-20) ? core_P/core_rho : 0.0);
        printf("  <T01>=%.6e  trace(rho-P)=%.6e\n", core_T01, core_trace);
        printf("  P = -rho? P/rho = %.6f (need -1.0 for vacuum EOS)\n",
               (fabs(core_rho) > 1e-20) ? core_P/core_rho : 0.0);
    }

    if (vac_n > 0) {
        vac_rho /= vac_n; vac_P /= vac_n; vac_T01 /= vac_n; vac_trace /= vac_n;
        printf("\nVACUUM (30<|x|<60): <rho>=%.6e  <P>=%.6e  P/rho=%.6f\n",
               vac_rho, vac_P, (fabs(vac_rho) > 1e-20) ? vac_P/vac_rho : 0.0);
        printf("  <T01>=%.6e  trace(rho-P)=%.6e\n", vac_T01, vac_trace);
        printf("  P = -rho? P/rho = %.6f (need -1.0 for vacuum EOS)\n",
               (fabs(vac_rho) > 1e-20) ? vac_P/vac_rho : 0.0);
    }

    if (mid_n > 0) {
        mid_rho /= mid_n; mid_P /= mid_n;
        printf("\nMID-RANGE (10<|x|<20): <rho>=%.6e  <P>=%.6e  P/rho=%.6f\n",
               mid_rho, mid_P, (fabs(mid_rho) > 1e-20) ? mid_P/mid_rho : 0.0);
    }

    /* Check center point specifically */
    printf("\nCENTER (x=0): rho=%.6e  P=%.6e  P/rho=%.6f  T01=%.6e  trace=%.6e\n",
           avg_T00[ic], avg_T11[ic],
           (fabs(avg_T00[ic]) > 1e-20) ? avg_T11[ic]/avg_T00[ic] : 0.0,
           avg_T01[ic], avg_T00[ic] - avg_T11[ic]);

    /* Check: total integrated T01 should be ~0 for static oscillon */
    double total_T01 = 0;
    for (int i = 1; i < Nx - 1; i++)
        total_T01 += avg_T01[i] * dx;
    printf("\nIntegrated <T01> = %.6e (should be ~0 for static)\n", total_T01);

    /* Total energy and pressure integral */
    double total_rho = 0, total_P = 0;
    for (int i = 1; i < Nx - 1; i++) {
        total_rho += avg_T00[i] * dx;
        total_P += avg_T11[i] * dx;
    }
    printf("Integrated <rho> = %.6f\n", total_rho);
    printf("Integrated <P>   = %.6f\n", total_P);
    printf("Virial check: int(P)/int(rho) = %.6f\n",
           (fabs(total_rho) > 1e-20) ? total_P / total_rho : 0.0);

    /* DFT of phi_1(0,t) — second half */
    int dft_start = n_dft / 2;
    double peak_om = 0;
    if (n_dft - dft_start > 100) {
        snprintf(path, sizeof(path), "%s/eos_spectrum%s.tsv", outdir, suffix);
        fp = fopen(path, "w");
        fprintf(fp, "omega\tpower\n");
        double T = t_hist[n_dft-1] - t_hist[dft_start];
        int nf = 500;
        double peak_pow = 0;
        for (int k = 0; k < nf; k++) {
            double omega = 3.0 * mass * k / nf;
            double re = 0, im = 0;
            for (int j = dft_start; j < n_dft; j++) {
                double dtj = (j > dft_start) ?
                    (t_hist[j] - t_hist[j-1]) : (t_hist[dft_start+1] - t_hist[dft_start]);
                re += phi0_hist[j] * cos(omega * t_hist[j]) * dtj;
                im += phi0_hist[j] * sin(omega * t_hist[j]) * dtj;
            }
            double pw = (re * re + im * im) / (T * T);
            fprintf(fp, "%.6f\t%.6e\n", omega, pw);
            if (pw > peak_pow) { peak_pow = pw; peak_om = omega; }
        }
        fclose(fp);
        printf("\nSpectrum: peak omega = %.4f (mass gap = %.4f)\n", peak_om, mass);
    }

    printf("\n=== END ===\n");

    /* Cleanup */
    for (int a = 0; a < 3; a++) { free(phi[a]); free(vel[a]); free(acc[a]); }
    free(damp); free(avg_T00); free(avg_T11); free(avg_T01);
    free(avg_T00sq); free(avg_T11sq);
    free(phi0_hist); free(t_hist);
    return 0;
}
