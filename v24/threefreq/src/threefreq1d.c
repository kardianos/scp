/*
 * threefreq1d.c — Three-frequency non-degenerate oscillon
 *
 * Three massive scalars with saturating triple-product coupling.
 * Each field initialized with a DIFFERENT oscillation frequency via
 * initial velocity: v_a(x,0) = -omega_a * A * g(x).
 *
 * Key question: do non-degenerate frequencies reduce radiation (dE/dt)?
 * Do frequencies converge back to degenerate or remain split?
 *
 * Based on v21/src/triad1d.c.
 *
 * Compile: gcc -O3 -Wall -o threefreq1d src/threefreq1d.c -lm
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Parameters */
static double mu     = -20.0;
static double kappa  = 20.0;
static double mass   = 1.0;
static double A_init = 0.8;
static double sigma  = 3.0;
static int    Nx     = 4000;
static double xmax   = 100.0;
static double tfinal = 10000.0;
static char   outdir[512] = "data";
static char   label[128]  = "control";

/* Three frequencies */
static double omega[3] = {0.87, 0.87, 0.87};

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
        else if (!strcmp(argv[i], "-o"))      strncpy(outdir, argv[i+1], sizeof(outdir)-1);
        else if (!strcmp(argv[i], "-label"))  strncpy(label, argv[i+1], sizeof(label)-1);
        else if (!strcmp(argv[i], "-w1"))     omega[0] = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-w2"))     omega[1] = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-w3"))     omega[2] = atof(argv[i+1]);
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

    printf("threefreq1d [%s]\n", label);
    printf("  omega = (%.4f, %.4f, %.4f)\n", omega[0], omega[1], omega[2]);
    printf("  mu=%.3f kappa=%.4f mass=%.4f A=%.3f sigma=%.3f\n",
           mu, kappa, mass, A_init, sigma);
    printf("  Nx=%d xmax=%.1f dx=%.5f dt=%.6f Nt=%d tfinal=%.0f\n",
           Nx, xmax, dx, dt, Nt, tfinal);

    /* Print frequency analysis */
    double O1 = omega[0]+omega[1]+omega[2];
    double O2 = omega[0]+omega[1]-omega[2];
    double O3 = omega[0]-omega[1]+omega[2];
    double O4 = -omega[0]+omega[1]+omega[2];
    printf("  Triple-product frequencies:\n");
    printf("    Omega_1 = %.4f (w1+w2+w3) %s\n", O1, fabs(O1)>mass?"ABOVE GAP":"below gap");
    printf("    Omega_2 = %.4f (w1+w2-w3) %s\n", O2, fabs(O2)>mass?"ABOVE GAP":"below gap");
    printf("    Omega_3 = %.4f (w1-w2+w3) %s\n", O3, fabs(O3)>mass?"ABOVE GAP":"below gap");
    printf("    Omega_4 = %.4f (-w1+w2+w3) %s\n", O4, fabs(O4)>mass?"ABOVE GAP":"below gap");

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

    /* Initialize: Gaussians with frequency-dependent initial velocity */
    for (int a = 0; a < 3; a++)
        for (int i = 0; i < Nx; i++) {
            double x = -xmax + i * dx;
            double g = exp(-x * x / (2.0 * sigma * sigma));
            phi[a][i] = A_init * g;
            vel[a][i] = -omega[a] * A_init * g;
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

    /* Output time series */
    char tspath[600];
    snprintf(tspath, sizeof(tspath), "%s/threefreq_%s_ts.tsv", outdir, label);
    FILE *fts = fopen(tspath, "w");
    if (!fts) { fprintf(stderr, "Cannot open %s\n", tspath); return 1; }
    fprintf(fts, "time\tphi1_0\tphi2_0\tphi3_0\tpeak1\tpeak2\tpeak3\t"
                 "E_kin\tE_grad\tE_mass\tE_pot\tE_total\tf_core\n");

    /* DFT storage for all three fields */
    int max_dft = 50000;
    double *phi_hist[3];
    for (int a = 0; a < 3; a++)
        phi_hist[a] = malloc(max_dft * sizeof(double));
    double *t_hist = malloc(max_dft * sizeof(double));
    /* Also store P(0,t) for triple-product DFT */
    double *P_hist = malloc(max_dft * sizeof(double));
    int n_dft = 0;

    int rec_every = Nt / 20000;
    if (rec_every < 1) rec_every = 1;
    int print_every = Nt / 40;
    if (print_every < 1) print_every = 1;
    int dft_every = Nt / max_dft;
    if (dft_every < 1) dft_every = 1;

    double core_r = 3.0 * sigma;
    int ic = Nx / 2;

    /* Windowed dE/dt: store E at regular intervals */
    int n_eslope = 0;
    int max_eslope = 200;
    double *E_slope_t = malloc(max_eslope * sizeof(double));
    double *E_slope_E = malloc(max_eslope * sizeof(double));
    int eslope_every = Nt / max_eslope;
    if (eslope_every < 1) eslope_every = 1;

    for (int n = 0; n <= Nt; n++) {
        double t = n * dt;

        if (n % dft_every == 0 && n_dft < max_dft) {
            for (int a = 0; a < 3; a++)
                phi_hist[a][n_dft] = phi[a][ic];
            P_hist[n_dft] = phi[0][ic] * phi[1][ic] * phi[2][ic];
            t_hist[n_dft] = t;
            n_dft++;
        }

        int do_rec = (n % rec_every == 0);
        int do_print = (n % print_every == 0);
        int do_eslope = (n % eslope_every == 0);

        if (do_rec || do_print || do_eslope) {
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

            if (do_eslope && n_eslope < max_eslope) {
                E_slope_t[n_eslope] = t;
                E_slope_E[n_eslope] = Et;
                n_eslope++;
            }

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
    printf("\nTime series: %s\n", tspath);

    /* === DFT analysis: per-field spectra (second half) === */
    int dft_start = n_dft / 2;
    int n_used = n_dft - dft_start;
    if (n_used > 100) {
        double T = t_hist[n_dft-1] - t_hist[dft_start];
        int nf = 600;
        double omega_max = 3.5 * mass;

        /* Per-field DFT */
        for (int a = 0; a < 3; a++) {
            char dftpath[600];
            snprintf(dftpath, sizeof(dftpath), "%s/threefreq_%s_spectrum%d.tsv",
                     outdir, label, a+1);
            FILE *fdft = fopen(dftpath, "w");
            fprintf(fdft, "omega\tpower\n");
            double peak_pow = 0, peak_om = 0;
            for (int k = 0; k < nf; k++) {
                double om = omega_max * k / nf;
                double re = 0, im = 0;
                for (int j = dft_start; j < n_dft; j++) {
                    double dtj = (j > dft_start) ?
                        (t_hist[j]-t_hist[j-1]) : (t_hist[dft_start+1]-t_hist[dft_start]);
                    re += phi_hist[a][j] * cos(om * t_hist[j]) * dtj;
                    im += phi_hist[a][j] * sin(om * t_hist[j]) * dtj;
                }
                double pw = (re*re + im*im) / (T*T);
                fprintf(fdft, "%.6f\t%.6e\n", om, pw);
                if (pw > peak_pow) { peak_pow = pw; peak_om = om; }
            }
            fclose(fdft);
            printf("  Field %d spectrum: peak omega = %.4f %s\n",
                   a+1, peak_om, peak_om < mass ? "(below gap)" : "(ABOVE gap)");
        }

        /* Triple-product DFT */
        {
            char dftpath[600];
            snprintf(dftpath, sizeof(dftpath), "%s/threefreq_%s_Pspectrum.tsv",
                     outdir, label);
            FILE *fdft = fopen(dftpath, "w");
            fprintf(fdft, "omega\tpower\n");
            double peak_pow = 0, peak_om = 0;
            for (int k = 0; k < nf; k++) {
                double om = omega_max * k / nf;
                double re = 0, im = 0;
                for (int j = dft_start; j < n_dft; j++) {
                    double dtj = (j > dft_start) ?
                        (t_hist[j]-t_hist[j-1]) : (t_hist[dft_start+1]-t_hist[dft_start]);
                    re += P_hist[j] * cos(om * t_hist[j]) * dtj;
                    im += P_hist[j] * sin(om * t_hist[j]) * dtj;
                }
                double pw = (re*re + im*im) / (T*T);
                fprintf(fdft, "%.6f\t%.6e\n", om, pw);
                if (pw > peak_pow) { peak_pow = pw; peak_om = om; }
            }
            fclose(fdft);
            printf("  P(t) spectrum: peak omega = %.4f %s\n",
                   peak_om, peak_om < mass ? "(below gap)" : "(ABOVE gap)");
        }
    }

    /* === Energy slope analysis === */
    {
        char slopepath[600];
        snprintf(slopepath, sizeof(slopepath), "%s/threefreq_%s_eslope.tsv",
                 outdir, label);
        FILE *fsl = fopen(slopepath, "w");
        fprintf(fsl, "t_mid\tdEdt\tE\n");

        /* Use windows of 10 points for local slope */
        int win = 10;
        double dEdt_early = 0, dEdt_late = 0;
        int n_early = 0, n_late = 0;
        double t_mid_boundary = tfinal * 0.5;

        for (int i = win; i < n_eslope; i++) {
            /* linear regression over window */
            double st = 0, sE = 0, st2 = 0, stE = 0;
            int nw = 0;
            for (int j = i - win; j <= i; j++) {
                st += E_slope_t[j];
                sE += E_slope_E[j];
                st2 += E_slope_t[j] * E_slope_t[j];
                stE += E_slope_t[j] * E_slope_E[j];
                nw++;
            }
            double slope = (nw*stE - st*sE) / (nw*st2 - st*st);
            double tmid = st / nw;
            double Emid = sE / nw;
            fprintf(fsl, "%.2f\t%.6e\t%.6e\n", tmid, slope, Emid);

            if (tmid < t_mid_boundary) { dEdt_early += slope; n_early++; }
            else { dEdt_late += slope; n_late++; }
        }
        fclose(fsl);

        if (n_early > 0 && n_late > 0) {
            printf("\n  Energy loss rate:\n");
            printf("    Early (t < %.0f): dE/dt = %.4e\n",
                   t_mid_boundary, dEdt_early / n_early);
            printf("    Late  (t > %.0f): dE/dt = %.4e\n",
                   t_mid_boundary, dEdt_late / n_late);
        }
    }

    /* === Summary line for aggregation === */
    {
        /* Final state: use last recorded energy */
        double Ef = 0;
        if (n_eslope > 0) Ef = E_slope_E[n_eslope - 1];
        double E0 = (n_eslope > 0) ? E_slope_E[0] : 0;

        /* Late dE/dt from last 20% */
        double dEdt_final = 0;
        int n_final = 0;
        int win = 10;
        for (int i = (int)(n_eslope * 0.8); i < n_eslope - win; i++) {
            double st = 0, sE = 0, st2 = 0, stE = 0;
            int nw = 0;
            for (int j = i; j < i + win && j < n_eslope; j++) {
                st += E_slope_t[j]; sE += E_slope_E[j];
                st2 += E_slope_t[j]*E_slope_t[j]; stE += E_slope_t[j]*E_slope_E[j];
                nw++;
            }
            if (nw > 1) {
                double slope = (nw*stE - st*sE) / (nw*st2 - st*st);
                dEdt_final += slope; n_final++;
            }
        }
        if (n_final > 0) dEdt_final /= n_final;

        char sumpath[600];
        snprintf(sumpath, sizeof(sumpath), "%s/threefreq_%s_summary.tsv", outdir, label);
        FILE *fsum = fopen(sumpath, "w");
        fprintf(fsum, "label\tw1\tw2\tw3\tmu\tkappa\tE0\tEf\tE_retained\tdEdt_late\n");
        fprintf(fsum, "%s\t%.4f\t%.4f\t%.4f\t%.1f\t%.1f\t%.4f\t%.4f\t%.4f\t%.4e\n",
                label, omega[0], omega[1], omega[2], mu, kappa,
                E0, Ef, (E0 > 1e-20) ? Ef/E0 : 0.0, dEdt_final);
        fclose(fsum);
        printf("\n  Summary: E0=%.4f  Ef=%.4f  retained=%.1f%%  dE/dt(late)=%.4e\n",
               E0, Ef, (E0>1e-20)?100*Ef/E0:0, dEdt_final);
    }

    printf("Output: %s\n", tspath);

    for (int a = 0; a < 3; a++) { free(phi[a]); free(vel[a]); free(acc[a]); free(phi_hist[a]); }
    free(damp); free(t_hist); free(P_hist);
    free(E_slope_t); free(E_slope_E);
    return 0;
}
