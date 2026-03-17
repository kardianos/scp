/*
 * confine.c — Test C: Confining potential for triad oscillon
 *
 * Three potential modes:
 *   mode 0: Saturating   V = (mu/2)P^2/(1+kappa*P^2)       [baseline]
 *   mode 1: Confining    V = -sigma*sqrt(P^2 + eps^2)       [linear confinement]
 *   mode 2: Conf+quartic V = -sigma*sqrt(P^2+eps^2) + (kc/2)*P^4  [stabilized]
 *
 * P = phi_1 * phi_2 * phi_3
 *
 * Based on v21/src/triad1d.c
 * Compile: gcc -O3 -Wall -o confine src/confine.c -lm
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Parameters */
static double mass   = 1.0;
static double A_init = 0.8;
static double sig_init = 3.0;   /* Gaussian width */
static int    Nx     = 4000;
static double xmax   = 100.0;
static double tfinal = 20000.0;
static char   outdir[512] = "v24/fundamental/testC_confine/data";

/* Potential parameters */
static int    pot_mode = 1;     /* 0=saturating, 1=confining, 2=conf+quartic */
static double mu       = -10.0;
static double kappa    = 0.1;
static double sigma_c  = 10.0;  /* string tension */
static double eps      = 1e-6;  /* regularization */
static double kc       = 0.0;   /* quartic stabilizer */

static void parse_args(int argc, char **argv)
{
    for (int i = 1; i < argc - 1; i += 2) {
        if      (!strcmp(argv[i], "-mode"))    pot_mode = atoi(argv[i+1]);
        else if (!strcmp(argv[i], "-mu"))      mu       = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-kappa"))   kappa    = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-sigma_c")) sigma_c  = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-eps"))     eps      = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-kc"))      kc       = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-mass"))    mass     = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-A"))       A_init   = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-sig"))     sig_init = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-Nx"))      Nx       = atoi(argv[i+1]);
        else if (!strcmp(argv[i], "-xmax"))    xmax     = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-tfinal"))  tfinal   = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-o"))       strncpy(outdir, argv[i+1], sizeof(outdir)-1);
        else { fprintf(stderr, "Unknown flag: %s\n", argv[i]); exit(1); }
    }
}

/*
 * Potential energy density V(P) at a single point
 */
static double V_pot(double P)
{
    double P2 = P * P;
    switch (pot_mode) {
        case 0: /* Saturating: V = (mu/2)*P^2/(1+kappa*P^2) */
            return 0.5 * mu * P2 / (1.0 + kappa * P2);
        case 1: /* Confining: V = -sigma*sqrt(P^2+eps^2) */
            return -sigma_c * sqrt(P2 + eps * eps);
        case 2: /* Confining + quartic: V = -sigma*sqrt(P^2+eps^2) + (kc/2)*P^4 */
            return -sigma_c * sqrt(P2 + eps * eps) + 0.5 * kc * P2 * P2;
        default:
            return 0.0;
    }
}

/*
 * Force: -dV/dphi_a
 * dP/dphi_a: for a=0, dP=phi2*phi3, etc.
 */
static double force_pot(double p1, double p2, double p3, int a)
{
    double P  = p1 * p2 * p3;
    double dP;
    switch (a) {
        case 0: dP = p2 * p3; break;
        case 1: dP = p1 * p3; break;
        case 2: dP = p1 * p2; break;
        default: dP = 0.0;
    }

    double P2 = P * P;
    switch (pot_mode) {
        case 0: { /* -dV/dphi_a = -mu*P*dP/(1+kappa*P^2)^2 */
            double denom2 = (1.0 + kappa * P2) * (1.0 + kappa * P2);
            return -mu * P * dP / denom2;
        }
        case 1: { /* -dV/dphi_a = sigma*P*dP/sqrt(P^2+eps^2) */
            double sq = sqrt(P2 + eps * eps);
            return sigma_c * P * dP / sq;
        }
        case 2: { /* confining + quartic */
            double sq = sqrt(P2 + eps * eps);
            return sigma_c * P * dP / sq - 2.0 * kc * P2 * P * dP;
        }
        default:
            return 0.0;
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
    double dt   = 0.8 * 2.0 / sqrt(kmax * kmax + m2);
    int    Nt   = (int)(tfinal / dt) + 1;

    const char *mode_desc;
    switch (pot_mode) {
        case 0: mode_desc = "Saturating V=(mu/2)P^2/(1+kP^2)"; break;
        case 1: mode_desc = "Confining V=-sigma*sqrt(P^2+eps^2)"; break;
        case 2: mode_desc = "Conf+quartic V=-sigma*sqrt(P^2+eps^2)+(kc/2)P^4"; break;
        default: mode_desc = "Unknown"; break;
    }

    printf("confine: mode %d — %s\n", pot_mode, mode_desc);
    printf("  mass=%.4f A=%.3f sig_init=%.3f\n", mass, A_init, sig_init);
    if (pot_mode == 0)
        printf("  mu=%.3f kappa=%.4f\n", mu, kappa);
    else if (pot_mode == 1)
        printf("  sigma_c=%.3f eps=%.2e\n", sigma_c, eps);
    else
        printf("  sigma_c=%.3f eps=%.2e kc=%.4f\n", sigma_c, eps, kc);
    printf("  Nx=%d xmax=%.1f dx=%.5f dt=%.6f Nt=%d tfinal=%.0f\n",
           Nx, xmax, dx, dt, Nt, tfinal);

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
        if (ax > x_abs) {
            double f = (ax - x_abs) / (xmax - x_abs);
            damp[i] = 1.0 - 0.98 * f * f;
        } else {
            damp[i] = 1.0;
        }
    }

    /* Initialize: Gaussians (symmetric triad) */
    for (int a = 0; a < 3; a++)
        for (int i = 0; i < Nx; i++) {
            double x = -xmax + i * dx;
            phi[a][i] = A_init * exp(-x * x / (2.0 * sig_init * sig_init));
        }

    /* Compute acceleration */
    #define COMPUTE_ACC() do { \
        for (int a = 0; a < 3; a++) { \
            acc[a][0] = acc[a][1] = acc[a][Nx-2] = acc[a][Nx-1] = 0; \
            for (int i = 1; i < Nx - 1; i++) { \
                double lapl = (phi[a][i+1] - 2.0*phi[a][i] + phi[a][i-1]) / dx2; \
                double fp = force_pot(phi[0][i], phi[1][i], phi[2][i], a); \
                acc[a][i] = lapl - m2*phi[a][i] + fp; \
            } \
        } \
    } while(0)

    COMPUTE_ACC();

    /* Output file */
    char tspath[600];
    snprintf(tspath, sizeof(tspath), "%s/confine_mode%d_sc%.1f_ts.tsv",
             outdir, pot_mode, (pot_mode == 0) ? mu : sigma_c);
    FILE *fts = fopen(tspath, "w");
    if (!fts) { fprintf(stderr, "Cannot open %s\n", tspath); return 1; }
    fprintf(fts, "time\tphi1_0\tpeak1\tE_kin\tE_grad\tE_mass\tE_pot\tE_total\t"
                 "fc\tdEdt\n");

    /* DFT storage */
    int max_dft = 50000;
    double *phi0_hist = malloc(max_dft * sizeof(double));
    double *t_hist    = malloc(max_dft * sizeof(double));
    int n_dft = 0;

    int rec_every   = Nt / 20000;
    if (rec_every < 1) rec_every = 1;
    int print_every = Nt / 40;
    if (print_every < 1) print_every = 1;
    int dft_every   = Nt / max_dft;
    if (dft_every < 1) dft_every = 1;

    double core_r = 3.0 * sig_init;
    int    ic     = Nx / 2;

    double E_prev = 0.0;
    double t_prev = 0.0;
    int    first_rec = 1;

    for (int n = 0; n <= Nt; n++) {
        double t = n * dt;

        if (n % dft_every == 0 && n_dft < max_dft) {
            phi0_hist[n_dft] = phi[0][ic];
            t_hist[n_dft] = t;
            n_dft++;
        }

        int do_rec   = (n % rec_every == 0);
        int do_print = (n % print_every == 0);

        if (do_rec || do_print) {
            double Ek = 0, Eg = 0, Em = 0, Ep = 0;
            double Ecore = 0, Eall = 0;
            double peak = 0;

            for (int i = 1; i < Nx - 1; i++) {
                double x = -xmax + i * dx;
                for (int a = 0; a < 3; a++) {
                    Ek += 0.5 * vel[a][i] * vel[a][i] * dx;
                    double dp = (phi[a][i+1] - phi[a][i-1]) / (2.0 * dx);
                    Eg += 0.5 * dp * dp * dx;
                    Em += 0.5 * m2 * phi[a][i] * phi[a][i] * dx;
                    if (fabs(phi[a][i]) > peak) peak = fabs(phi[a][i]);
                }

                double P = phi[0][i] * phi[1][i] * phi[2][i];
                double V = V_pot(P);
                Ep += V * dx;

                double e = V;
                for (int a = 0; a < 3; a++) {
                    e += 0.5 * vel[a][i] * vel[a][i];
                    double dp = (phi[a][i+1] - phi[a][i-1]) / (2.0 * dx);
                    e += 0.5 * dp * dp + 0.5 * m2 * phi[a][i] * phi[a][i];
                }
                Eall += e * dx;
                if (fabs(x) < core_r) Ecore += e * dx;
            }

            double Et = Ek + Eg + Em + Ep;
            double fc = (Eall > 1e-20) ? Ecore / Eall : 0.0;

            double dEdt = 0.0;
            if (!first_rec && t > t_prev) {
                dEdt = (Et - E_prev) / (t - t_prev);
            }
            E_prev = Et;
            t_prev = t;
            first_rec = 0;

            if (do_rec)
                fprintf(fts, "%.6f\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t"
                             "%.6f\t%.6e\n",
                        t, phi[0][ic], peak, Ek, Eg, Em, Ep, Et, fc, dEdt);
            if (do_print)
                printf("  t=%7.1f  phi0=(%+.4f)  pk=%.4f  E=%+.4f  Ep=%+.4f  "
                       "fc=%.3f  dE/dt=%.2e\n",
                       t, phi[0][ic], peak, Et, Ep, fc, dEdt);
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

        /* absorbing boundary */
        for (int a = 0; a < 3; a++)
            for (int i = 0; i < Nx; i++) {
                vel[a][i] *= damp[i];
                phi[a][i] *= damp[i];
            }
    }

    fclose(fts);

    /* DFT of phi_1(0,t) — second half */
    int dft_start = n_dft / 2;
    if (n_dft - dft_start > 100) {
        char dftpath[600];
        snprintf(dftpath, sizeof(dftpath), "%s/confine_mode%d_sc%.1f_spectrum.tsv",
                 outdir, pot_mode, (pot_mode == 0) ? mu : sigma_c);
        FILE *fdft = fopen(dftpath, "w");
        fprintf(fdft, "omega\tpower\n");
        double T = t_hist[n_dft-1] - t_hist[dft_start];
        int nf = 500;
        double peak_pow = 0, peak_om = 0;
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
            fprintf(fdft, "%.6f\t%.6e\n", omega, pw);
            if (pw > peak_pow) { peak_pow = pw; peak_om = omega; }
        }
        fclose(fdft);
        printf("\nSpectrum: peak omega = %.4f (mass gap = %.4f)\n", peak_om, mass);
        printf("Oscillon (omega < m)? %s\n",
               (peak_om > 0.01 && peak_om < mass) ? "YES" : "NO");
    }

    printf("Output: %s\n", tspath);

    for (int a = 0; a < 3; a++) { free(phi[a]); free(vel[a]); free(acc[a]); }
    free(damp); free(phi0_hist); free(t_hist);
    return 0;
}
