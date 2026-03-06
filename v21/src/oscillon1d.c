/*
 * oscillon1d.c — 1D oscillon search: massive scalar + biharmonic + saturating φ^6
 *
 * Symmetric ansatz: phi_1 = phi_2 = phi_3 = phi, so P = phi^3
 * Effective single-field Lagrangian:
 *   L = (3/2)(dt phi)^2 - (3/2)(dx phi)^2 - (3 beta/2)(d^2x phi)^2
 *     - (3/2) m^2 phi^2 - (mu/2) phi^6 / (1 + kappa phi^6)
 *
 * Simplify: divide by 3, redefine mu_eff = mu/3, keep m^2 same:
 *   L_eff = (1/2)(dt phi)^2 - (1/2)(dx phi)^2 - (beta/2)(d^2x phi)^2
 *         - (m^2/2) phi^2 - (mu_eff/2) phi^6 / (1 + kappa phi^6)
 *
 * EOM: d^2t phi = d^2x phi - beta d^4x phi - m^2 phi - dV_eff/dphi
 *   dV_eff/dphi = 3 mu_eff phi^5 / (1 + kappa phi^6)^2
 *
 * Dispersion: omega^2 = k^2 + beta k^4 + m^2
 *   Gap: omega_min = m (at k=0)
 *   Oscillon: oscillates at omega < m -> cannot radiate
 *
 * Compile: gcc -O3 -Wall -o oscillon1d v21/src/oscillon1d.c -lm
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static double mu_eff = -20.0;   /* effective coupling (mu/3) */
static double kappa  = 1.0;
static double beta   = 0.1;
static double mass   = 1.0;
static double A_init = 1.5;
static double sigma  = 2.0;
static int    Nx     = 1000;
static double xmax   = 30.0;
static double tfinal = 500.0;
static char   outdir[512] = "v21/data";

static void parse_args(int argc, char **argv)
{
    for (int i = 1; i < argc - 1; i += 2) {
        if      (!strcmp(argv[i], "-mu"))     mu_eff = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-kappa"))  kappa  = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-beta"))   beta   = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-mass"))   mass   = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-A"))      A_init = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-sigma"))  sigma  = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-Nx"))     Nx     = atoi(argv[i+1]);
        else if (!strcmp(argv[i], "-xmax"))   xmax   = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-tfinal")) tfinal = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-o"))      strncpy(outdir, argv[i+1], sizeof(outdir)-1);
        else { fprintf(stderr, "Unknown flag: %s\n", argv[i]); exit(1); }
    }
}

int main(int argc, char **argv)
{
    parse_args(argc, argv);

    double dx = 2.0 * xmax / (Nx - 1);
    double dx2 = dx * dx;
    double dx4 = dx2 * dx2;
    double m2 = mass * mass;

    /* CFL: omega_max^2 = (pi/dx)^2 + beta*(pi/dx)^4 + m^2
     * dt < 2/omega_max */
    double kmax = M_PI / dx;
    double omega_max = sqrt(kmax * kmax + beta * kmax * kmax * kmax * kmax + m2);
    double dt = 0.4 / omega_max;  /* safety factor 0.4 */
    int Nt = (int)(tfinal / dt) + 1;

    printf("oscillon1d: mu_eff=%.3f kappa=%.4f beta=%.4f mass=%.4f A=%.3f sigma=%.3f\n",
           mu_eff, kappa, beta, mass, A_init, sigma);
    printf("  Nx=%d xmax=%.1f dx=%.5f dt=%.6f omega_max=%.2f Nt=%d\n",
           Nx, xmax, dx, dt, omega_max, Nt);
    printf("  Mass gap: omega_min = m = %.4f\n", mass);
    printf("  Natural freq estimate (small amp): omega ~ sqrt(m^2 + mu_eff * ...) \n");

    double *phi = calloc(Nx, sizeof(double));
    double *vel = calloc(Nx, sizeof(double));
    double *acc = calloc(Nx, sizeof(double));
    double *damp = malloc(Nx * sizeof(double));

    /* Absorbing layer: outer 20% */
    double absorb_frac = 0.20;
    double x_abs = xmax * (1.0 - absorb_frac);
    for (int i = 0; i < Nx; i++) {
        double x = -xmax + i * dx;
        double ax = fabs(x);
        if (ax > x_abs) {
            double f = (ax - x_abs) / (xmax * absorb_frac);
            damp[i] = 1.0 - 0.98 * f * f;  /* quadratic ramp */
        } else {
            damp[i] = 1.0;
        }
    }

    /* Initialize: Gaussian bump, zero velocity */
    for (int i = 0; i < Nx; i++) {
        double x = -xmax + i * dx;
        phi[i] = A_init * exp(-x * x / (2.0 * sigma * sigma));
    }

    /* Compute acceleration */
    #define COMPUTE_ACC() do { \
        for (int i = 2; i < Nx - 2; i++) { \
            double lapl = (phi[i+1] - 2.0*phi[i] + phi[i-1]) / dx2; \
            double bih  = (phi[i-2] - 4.0*phi[i-1] + 6.0*phi[i] \
                          - 4.0*phi[i+1] + phi[i+2]) / dx4; \
            double p = phi[i]; \
            double p5 = p*p*p*p*p; \
            double p6 = p5*p; \
            double denom2 = (1.0 + kappa*p6); denom2 *= denom2; \
            double dVdp = 3.0 * mu_eff * p5 / denom2; \
            acc[i] = lapl - beta*bih - m2*p - dVdp; \
        } \
        acc[0] = acc[1] = acc[Nx-2] = acc[Nx-1] = 0.0; \
    } while(0)

    COMPUTE_ACC();

    /* Output */
    char tspath[600];
    snprintf(tspath, sizeof(tspath), "%s/oscillon1d_timeseries.tsv", outdir);
    FILE *fts = fopen(tspath, "w");
    if (!fts) { fprintf(stderr, "Cannot open %s\n", tspath); return 1; }
    fprintf(fts, "time\tphi0\tpeak\tE_kin\tE_grad\tE_bih\tE_mass\tE_pot\tE_total\tf_core\n");

    /* Record center value for DFT */
    int max_dft = 100000;
    double *phi0_hist = malloc(max_dft * sizeof(double));
    double *t_hist = malloc(max_dft * sizeof(double));
    int n_dft = 0;

    int rec_every = Nt / 10000;
    if (rec_every < 1) rec_every = 1;
    int print_every = Nt / 30;
    if (print_every < 1) print_every = 1;
    int dft_every = Nt / max_dft;
    if (dft_every < 1) dft_every = 1;

    double core_radius = 3.0 * sigma;
    int i_center = Nx / 2;

    for (int n = 0; n <= Nt; n++) {
        double t = n * dt;

        /* record phi(0) for DFT */
        if (n % dft_every == 0 && n_dft < max_dft) {
            phi0_hist[n_dft] = phi[i_center];
            t_hist[n_dft] = t;
            n_dft++;
        }

        int do_rec = (n % rec_every == 0);
        int do_print = (n % print_every == 0);

        if (do_rec || do_print) {
            double Ek = 0, Eg = 0, Eb = 0, Em = 0, Ep = 0;
            double Ecore = 0, Eall = 0;
            double peak = 0;

            for (int i = 2; i < Nx - 2; i++) {
                double x = -xmax + i * dx;
                Ek += 0.5 * vel[i] * vel[i] * dx;
                double dphi = (phi[i+1] - phi[i-1]) / (2.0*dx);
                Eg += 0.5 * dphi * dphi * dx;
                double d2phi = (phi[i+1] - 2.0*phi[i] + phi[i-1]) / dx2;
                Eb += 0.5 * beta * d2phi * d2phi * dx;
                Em += 0.5 * m2 * phi[i] * phi[i] * dx;
                double p = phi[i], p6 = p*p*p*p*p*p;
                double V = 0.5 * mu_eff * p6 / (1.0 + kappa*p6);
                Ep += V * dx;

                double e = 0.5*vel[i]*vel[i] + 0.5*dphi*dphi
                         + 0.5*beta*d2phi*d2phi + 0.5*m2*phi[i]*phi[i] + V;
                Eall += e * dx;
                if (fabs(x) < core_radius)
                    Ecore += e * dx;
                if (fabs(phi[i]) > peak) peak = fabs(phi[i]);
            }

            double Et = Ek + Eg + Eb + Em + Ep;
            double fc = (Eall > 1e-20) ? Ecore / Eall : 0.0;

            if (do_rec)
                fprintf(fts, "%.6f\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6f\n",
                        t, phi[i_center], peak, Ek, Eg, Eb, Em, Ep, Et, fc);
            if (do_print)
                printf("  t=%7.1f  phi0=%+.4f  pk=%.4f  E=%+.3f  Ek=%.2f  Eg=%.2f  Eb=%.2f  "
                       "Em=%.2f  Ep=%+.2f  fc=%.3f\n",
                       t, phi[i_center], peak, Et, Ek, Eg, Eb, Em, Ep, fc);
        }

        if (n == Nt) break;

        /* Velocity Verlet */
        for (int i = 2; i < Nx - 2; i++)
            vel[i] += 0.5 * dt * acc[i];
        for (int i = 2; i < Nx - 2; i++)
            phi[i] += dt * vel[i];
        COMPUTE_ACC();
        for (int i = 2; i < Nx - 2; i++)
            vel[i] += 0.5 * dt * acc[i];

        /* absorbing boundary */
        for (int i = 0; i < Nx; i++) {
            vel[i] *= damp[i];
            phi[i] *= damp[i];  /* also damp phi to push toward vacuum at boundary */
        }
    }

    fclose(fts);

    /* DFT of phi(0,t) */
    if (n_dft > 100) {
        char dftpath[600];
        snprintf(dftpath, sizeof(dftpath), "%s/oscillon1d_spectrum.tsv", outdir);
        FILE *fdft = fopen(dftpath, "w");
        if (fdft) {
            fprintf(fdft, "omega\tpower\n");
            double T = t_hist[n_dft-1] - t_hist[0];
            int n_omega = 500;
            double omega_max_dft = 5.0 * mass;
            for (int k = 0; k < n_omega; k++) {
                double omega = omega_max_dft * k / n_omega;
                double re = 0, im = 0;
                for (int j = 0; j < n_dft; j++) {
                    double dt_j = (j > 0) ? (t_hist[j] - t_hist[j-1]) : (t_hist[1] - t_hist[0]);
                    re += phi0_hist[j] * cos(omega * t_hist[j]) * dt_j;
                    im += phi0_hist[j] * sin(omega * t_hist[j]) * dt_j;
                }
                double power = (re*re + im*im) / (T*T);
                fprintf(fdft, "%.6f\t%.6e\n", omega, power);
            }
            fclose(fdft);
            printf("Spectrum written to %s\n", dftpath);
        }
    }

    printf("Timeseries written to %s\n", tspath);
    printf("\n=== OSCILLON TEST ===\n");
    printf("  Mass gap: omega = %.4f\n", mass);
    printf("  If phi(0) oscillates coherently with omega < m, this is an oscillon.\n");
    printf("  Check spectrum for peak below m = %.4f\n", mass);
    printf("===\n");

    free(phi); free(vel); free(acc); free(damp);
    free(phi0_hist); free(t_hist);
    return 0;
}
