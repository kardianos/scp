/*
 * triad1d.c — 1D three-body oscillon: three massive scalars with triple-product coupling
 *
 * Lagrangian:
 *   L = sum_a [ (1/2)(dt phi_a)^2 - (1/2)(dx phi_a)^2 - (m^2/2)phi_a^2 ]
 *     - (mu/2) P^2 / (1 + kappa P^2)
 *   P = phi_1 phi_2 phi_3
 *
 * Each field individually: free massive Klein-Gordon (disperses).
 * Triple product coupling: attractive when all three overlap.
 * Mass gap m prevents radiation when oscillon frequency omega < m.
 *
 * Test modes:
 *   1: All three identical (symmetric, P = phi^3)
 *   2: Three different amplitudes (A, 0.8A, 0.6A)
 *   3: Two fields only (control — should disperse since P=0)
 *   4: Single field (control — free massive, disperses)
 *
 * Compile: gcc -O3 -Wall -o triad1d v21/src/triad1d.c -lm
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static double mu     = -10.0;
static double kappa  = 0.1;
static double mass   = 1.0;
static double A_init = 0.8;
static double sigma  = 3.0;
static int    Nx     = 2000;
static double xmax   = 60.0;
static double tfinal = 2000.0;
static int    test   = 1;
static char   outdir[512] = "v21/data";

static void parse_args(int argc, char **argv)
{
    for (int i = 1; i < argc - 1; i += 2) {
        if      (!strcmp(argv[i], "-mu"))     mu     = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-kappa"))  kappa  = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-mass"))   mass   = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-A"))      A_init = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-sigma"))  sigma  = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-Nx"))     Nx     = atoi(argv[i+1]);
        else if (!strcmp(argv[i], "-xmax"))   xmax   = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-tfinal")) tfinal = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-test"))   test   = atoi(argv[i+1]);
        else if (!strcmp(argv[i], "-o"))      strncpy(outdir, argv[i+1], sizeof(outdir)-1);
        else { fprintf(stderr, "Unknown flag: %s\n", argv[i]); exit(1); }
    }
}

/* -dV/dphi_a where V = (mu/2)P^2/(1+kappa P^2), P = phi1 phi2 phi3 */
static double force_pot(double p1, double p2, double p3, int a)
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

int main(int argc, char **argv)
{
    parse_args(argc, argv);

    double dx = 2.0 * xmax / (Nx - 1);
    double dx2 = dx * dx;
    double m2 = mass * mass;

    /* CFL */
    double kmax = M_PI / dx;
    double dt = 0.8 * 2.0 / sqrt(kmax*kmax + m2);
    int Nt = (int)(tfinal / dt) + 1;

    /* Initialize amplitudes per test */
    double Ainit[3] = {0, 0, 0};
    const char *test_desc;
    switch (test) {
        case 1: Ainit[0]=Ainit[1]=Ainit[2]=A_init;
                test_desc = "Symmetric triad (A,A,A)"; break;
        case 2: Ainit[0]=A_init; Ainit[1]=0.8*A_init; Ainit[2]=0.6*A_init;
                test_desc = "Asymmetric triad (A,0.8A,0.6A)"; break;
        case 3: Ainit[0]=A_init; Ainit[1]=A_init; Ainit[2]=0;
                test_desc = "Two-field control (A,A,0)"; break;
        case 4: Ainit[0]=A_init; Ainit[1]=0; Ainit[2]=0;
                test_desc = "Single-field control (A,0,0)"; break;
        default: fprintf(stderr, "Unknown test %d\n", test); return 1;
    }

    printf("triad1d test %d: %s\n", test, test_desc);
    printf("  mu=%.3f kappa=%.4f mass=%.4f A=%.3f sigma=%.3f\n",
           mu, kappa, mass, A_init, sigma);
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

    /* Initialize: Gaussians */
    for (int a = 0; a < 3; a++)
        for (int i = 0; i < Nx; i++) {
            double x = -xmax + i * dx;
            phi[a][i] = Ainit[a] * exp(-x * x / (2.0 * sigma * sigma));
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

    /* Output */
    char tspath[600];
    snprintf(tspath, sizeof(tspath), "%s/triad1d_test%d_ts.tsv", outdir, test);
    FILE *fts = fopen(tspath, "w");
    fprintf(fts, "time\tphi1_0\tphi2_0\tphi3_0\tpeak1\tpeak2\tpeak3\t"
                 "E_kin\tE_grad\tE_mass\tE_pot\tE_total\tf_core\n");

    /* DFT storage */
    int max_dft = 50000;
    double *phi0_hist = malloc(max_dft * sizeof(double));
    double *t_hist = malloc(max_dft * sizeof(double));
    int n_dft = 0;

    int rec_every = Nt / 20000;
    if (rec_every < 1) rec_every = 1;
    int print_every = Nt / 40;
    if (print_every < 1) print_every = 1;
    int dft_every = Nt / max_dft;
    if (dft_every < 1) dft_every = 1;

    double core_r = 3.0 * sigma;
    int ic = Nx / 2;

    for (int n = 0; n <= Nt; n++) {
        double t = n * dt;

        if (n % dft_every == 0 && n_dft < max_dft) {
            phi0_hist[n_dft] = phi[0][ic];
            t_hist[n_dft] = t;
            n_dft++;
        }

        int do_rec = (n % rec_every == 0);
        int do_print = (n % print_every == 0);

        if (do_rec || do_print) {
            double Ek = 0, Eg = 0, Em = 0, Ep = 0;
            double Ecore = 0, Eall = 0;
            double peak[3] = {0};

            for (int i = 1; i < Nx - 1; i++) {
                double x = -xmax + i * dx;
                for (int a = 0; a < 3; a++) {
                    Ek += 0.5 * vel[a][i] * vel[a][i] * dx;
                    double dp = (phi[a][i+1] - phi[a][i-1]) / (2.0*dx);
                    Eg += 0.5 * dp * dp * dx;
                    Em += 0.5 * m2 * phi[a][i] * phi[a][i] * dx;
                    if (fabs(phi[a][i]) > peak[a]) peak[a] = fabs(phi[a][i]);
                }

                double P = phi[0][i] * phi[1][i] * phi[2][i];
                double P2 = P * P;
                double V = 0.5 * mu * P2 / (1.0 + kappa * P2);
                Ep += V * dx;

                /* total density */
                double e = V;
                for (int a = 0; a < 3; a++) {
                    e += 0.5*vel[a][i]*vel[a][i];
                    double dp = (phi[a][i+1] - phi[a][i-1]) / (2.0*dx);
                    e += 0.5*dp*dp + 0.5*m2*phi[a][i]*phi[a][i];
                }
                Eall += e * dx;
                if (fabs(x) < core_r) Ecore += e * dx;
            }

            double Et = Ek + Eg + Em + Ep;
            double fc = (Eall > 1e-20) ? Ecore / Eall : 0.0;

            if (do_rec)
                fprintf(fts, "%.6f\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t"
                             "%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6f\n",
                        t, phi[0][ic], phi[1][ic], phi[2][ic],
                        peak[0], peak[1], peak[2], Ek, Eg, Em, Ep, Et, fc);
            if (do_print)
                printf("  t=%7.1f  p0=(%+.3f,%+.3f,%+.3f)  pk=(%.3f,%.3f,%.3f)  "
                       "E=%+.4f  Ep=%+.4f  fc=%.3f\n",
                       t, phi[0][ic], phi[1][ic], phi[2][ic],
                       peak[0], peak[1], peak[2], Et, Ep, fc);
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
        snprintf(dftpath, sizeof(dftpath), "%s/triad1d_test%d_spectrum.tsv", outdir, test);
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
                    (t_hist[j]-t_hist[j-1]) : (t_hist[dft_start+1]-t_hist[dft_start]);
                re += phi0_hist[j] * cos(omega * t_hist[j]) * dtj;
                im += phi0_hist[j] * sin(omega * t_hist[j]) * dtj;
            }
            double pw = (re*re + im*im) / (T*T);
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
