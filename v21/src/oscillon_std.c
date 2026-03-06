/*
 * oscillon_std.c — Standard oscillon model for validation
 *
 * V(phi) = (m^2/2)phi^2 - (g/3)|phi|^3 + (lambda/4)phi^4
 * (The phi^3 is even: -g/3 |phi|^3 = -(g/3)phi^2 |phi| to preserve Z_2)
 * Actually, simpler: use the standard asymmetric model:
 *   V(phi) = (m^2/2)phi^2 (1 - phi/phi_0)^2 ... no, let's use:
 *
 * Standard sine-Gordon-like:
 *   V(phi) = m^2(1 - cos(phi))  [periodic, mass gap = m at phi=0]
 *   Force: -dV/dphi = -m^2 sin(phi)
 *   EOM: phi_tt = phi_xx - beta phi_xxxx - m^2 sin(phi)
 *   Oscillon: breather at omega < m
 *
 * Alternative: phi^4 with negative cubic (Segur-Kruskal model):
 *   V = (m^2/2)phi^2 - (g/4)phi^4
 *   Force: -m^2 phi + g phi^3
 *   Well-known to support oscillons for 1D.
 *   (Needs phi^6 or boundary for global stability.)
 *
 * We use: V = (m^2/2)phi^2 - (g/4)phi^4 + (h/6)phi^6
 *   With g>0, h>0: attractive phi^4, stabilizing phi^6.
 *   Oscillon: omega ~ m sqrt(1 - c*A^2) < m
 *
 * Compile: gcc -O3 -Wall -o oscillon_std v21/src/oscillon_std.c -lm
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static double mass   = 1.0;
static double g4     = 2.0;    /* phi^4 attractive coupling */
static double h6     = 1.0;    /* phi^6 stabilizing */
static double beta   = 0.0;    /* biharmonic (optional) */
static double A_init = 0.5;
static double sigma  = 3.0;
static int    Nx     = 2000;
static double xmax   = 60.0;
static double tfinal = 1000.0;
static char   outdir[512] = "v21/data";

static void parse_args(int argc, char **argv)
{
    for (int i = 1; i < argc - 1; i += 2) {
        if      (!strcmp(argv[i], "-mass"))   mass   = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-g4"))     g4     = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-h6"))     h6     = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-beta"))   beta   = atof(argv[i+1]);
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

    /* CFL */
    double kmax = M_PI / dx;
    double om_max2 = kmax*kmax + beta*kmax*kmax*kmax*kmax + m2;
    double dt = 0.8 * 2.0 / sqrt(om_max2);
    if (beta > 0) {
        double dt_bih = 0.25 * dx2 / sqrt(beta);
        if (dt_bih < dt) dt = dt_bih;
    }
    int Nt = (int)(tfinal / dt) + 1;

    printf("oscillon_std: m=%.3f g4=%.3f h6=%.3f beta=%.4f A=%.3f sigma=%.3f\n",
           mass, g4, h6, beta, A_init, sigma);
    printf("  Nx=%d xmax=%.1f dx=%.5f dt=%.6f Nt=%d\n", Nx, xmax, dx, dt, Nt);

    double *phi = calloc(Nx, sizeof(double));
    double *vel = calloc(Nx, sizeof(double));
    double *acc = calloc(Nx, sizeof(double));

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

    /* Initialize: Gaussian */
    for (int i = 0; i < Nx; i++) {
        double x = -xmax + i * dx;
        phi[i] = A_init * exp(-x * x / (2.0 * sigma * sigma));
    }

    /* Force = -dV/dphi = -m^2 phi + g4 phi^3 - h6 phi^5 */
    #define FORCE(p) (-(m2)*(p) + (g4)*(p)*(p)*(p) - (h6)*(p)*(p)*(p)*(p)*(p))

    #define COMPUTE_ACC() do { \
        acc[0] = acc[1] = acc[Nx-2] = acc[Nx-1] = 0; \
        for (int i = 2; i < Nx - 2; i++) { \
            double lapl = (phi[i+1] - 2.0*phi[i] + phi[i-1]) / dx2; \
            double bih = 0; \
            if (beta > 0) bih = (phi[i-2] - 4.0*phi[i-1] + 6.0*phi[i] \
                                - 4.0*phi[i+1] + phi[i+2]) / dx4; \
            acc[i] = lapl - beta*bih + FORCE(phi[i]); \
        } \
    } while(0)

    COMPUTE_ACC();

    /* Output */
    char tspath[600];
    snprintf(tspath, sizeof(tspath), "%s/oscillon_std_ts.tsv", outdir);
    FILE *fts = fopen(tspath, "w");
    fprintf(fts, "time\tphi0\tpeak\tE_total\tf_core\n");

    int max_dft = 50000;
    double *phi0_hist = malloc(max_dft * sizeof(double));
    double *t_hist = malloc(max_dft * sizeof(double));
    int n_dft = 0;

    int rec_every = Nt / 20000;
    if (rec_every < 1) rec_every = 1;
    int print_every = Nt / 30;
    if (print_every < 1) print_every = 1;
    int dft_every = Nt / max_dft;
    if (dft_every < 1) dft_every = 1;

    double core_r = 3.0 * sigma;
    int ic = Nx / 2;

    for (int n = 0; n <= Nt; n++) {
        double t = n * dt;

        if (n % dft_every == 0 && n_dft < max_dft) {
            phi0_hist[n_dft] = phi[ic];
            t_hist[n_dft] = t;
            n_dft++;
        }

        int do_rec = (n % rec_every == 0);
        int do_print = (n % print_every == 0);

        if (do_rec || do_print) {
            double Ek = 0, Eg = 0, Eb = 0, Ep = 0;
            double Ecore = 0, Eall = 0;
            double peak = 0;

            for (int i = 2; i < Nx - 2; i++) {
                double x = -xmax + i * dx;
                Ek += 0.5 * vel[i] * vel[i] * dx;
                double dp = (phi[i+1] - phi[i-1]) / (2.0*dx);
                Eg += 0.5 * dp * dp * dx;
                double p = phi[i];
                double V = m2*p*p/2 - g4*p*p*p*p/4 + h6*p*p*p*p*p*p/6;
                Ep += V * dx;
                if (beta > 0) {
                    double d2p = (phi[i+1] - 2.0*phi[i] + phi[i-1]) / dx2;
                    Eb += 0.5 * beta * d2p * d2p * dx;
                }

                double e = 0.5*vel[i]*vel[i] + 0.5*dp*dp + V;
                Eall += e * dx;
                if (fabs(x) < core_r) Ecore += e * dx;
                if (fabs(p) > peak) peak = fabs(p);
            }

            double Et = Ek + Eg + Eb + Ep;
            double fc = (Eall > 1e-20) ? Ecore / Eall : 0.0;

            if (do_rec)
                fprintf(fts, "%.6f\t%.6e\t%.6e\t%.6e\t%.6f\n",
                        t, phi[ic], peak, Et, fc);
            if (do_print)
                printf("  t=%7.1f  phi0=%+.4f  pk=%.4f  E=%+.4f  fc=%.3f  "
                       "Ek=%.3f  Eg=%.3f  Ep=%.3f\n",
                       t, phi[ic], peak, Et, fc, Ek, Eg, Ep);
        }

        if (n == Nt) break;

        /* Velocity Verlet */
        for (int i = 2; i < Nx - 2; i++) vel[i] += 0.5 * dt * acc[i];
        for (int i = 2; i < Nx - 2; i++) phi[i] += dt * vel[i];
        COMPUTE_ACC();
        for (int i = 2; i < Nx - 2; i++) vel[i] += 0.5 * dt * acc[i];

        for (int i = 0; i < Nx; i++) {
            vel[i] *= damp[i];
            phi[i] *= damp[i];
        }
    }

    fclose(fts);

    /* DFT — only use second half of time series */
    int dft_start = n_dft / 2;
    if (n_dft - dft_start > 100) {
        char dftpath[600];
        snprintf(dftpath, sizeof(dftpath), "%s/oscillon_std_spectrum.tsv", outdir);
        FILE *fdft = fopen(dftpath, "w");
        fprintf(fdft, "omega\tpower\n");
        double T = t_hist[n_dft-1] - t_hist[dft_start];
        int nf = 500;
        for (int k = 0; k < nf; k++) {
            double omega = 3.0 * mass * k / nf;
            double re = 0, im = 0;
            for (int j = dft_start; j < n_dft; j++) {
                double dt_j = (j > dft_start) ?
                    (t_hist[j] - t_hist[j-1]) : (t_hist[dft_start+1] - t_hist[dft_start]);
                double val = phi0_hist[j];
                re += val * cos(omega * t_hist[j]) * dt_j;
                im += val * sin(omega * t_hist[j]) * dt_j;
            }
            fprintf(fdft, "%.6f\t%.6e\n", omega, (re*re + im*im) / (T*T));
        }
        fclose(fdft);
        printf("Spectrum (2nd half) written to %s\n", dftpath);

        /* Find peak */
        double max_pow = 0, peak_omega = 0;
        FILE *fread = fopen(dftpath, "r");
        char line[256];
        fgets(line, sizeof(line), fread); /* header */
        while (fgets(line, sizeof(line), fread)) {
            double om, pw;
            sscanf(line, "%lf\t%lf", &om, &pw);
            if (pw > max_pow) { max_pow = pw; peak_omega = om; }
        }
        fclose(fread);
        printf("  Peak frequency: omega = %.4f (mass gap = %.4f)\n", peak_omega, mass);
        printf("  Oscillon? %s\n", (peak_omega < mass && peak_omega > 0.1) ? "YES (omega < m)" : "NO");
    }

    free(phi); free(vel); free(acc); free(damp);
    free(phi0_hist); free(t_hist);
    return 0;
}
