/*
 * mediator1d.c — 1D three-field oscillon with per-field masses
 *
 * V23-E Path A: Break S3 symmetry by giving field 3 a different mass.
 * Fields 1,2 have mass m1=m2=1.0 (massive), field 3 has mass m3 (scanned).
 *
 * Lagrangian:
 *   L = sum_a [ (1/2)(dt phi_a)^2 - (1/2)(dx phi_a)^2 - (m_a^2/2)phi_a^2 ]
 *     - (mu/2) P^2 / (1 + kappa P^2)
 *   P = phi_1 phi_2 phi_3
 *
 * When m3=0, field 3 is massless: it gets tachyonic effective mass inside the
 * oscillon core from the coupling, potentially creating a localized "charge"
 * that sources a long-range 1/|x| tail.
 *
 * Compile: gcc -O3 -Wall -o mediator1d src/mediator1d.c -lm
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static double mu     = -20.0;
static double kappa  = 20.0;
static double m1     = 1.0;    /* mass of fields 1,2 */
static double m3     = 1.0;    /* mass of field 3 (scanned) */
static double A_init = 0.8;
static double sigma  = 3.0;
static int    Nx     = 4000;
static double xmax   = 100.0;
static double tfinal = 10000.0;
static int    do_scan = 0;     /* 1 = run Phase 2 scan */
static char   outdir[512] = "v23/mediator/data";

static void parse_args(int argc, char **argv)
{
    for (int i = 1; i < argc - 1; i += 2) {
        if      (!strcmp(argv[i], "-mu"))     mu     = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-kappa"))  kappa  = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-m1"))     m1     = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-m3"))     m3     = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-A"))      A_init = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-sigma"))  sigma  = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-Nx"))     Nx     = atoi(argv[i+1]);
        else if (!strcmp(argv[i], "-xmax"))   xmax   = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-tfinal")) tfinal = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-scan"))   do_scan = atoi(argv[i+1]);
        else if (!strcmp(argv[i], "-o"))      strncpy(outdir, argv[i+1], sizeof(outdir)-1);
        else { fprintf(stderr, "Unknown flag: %s\n", argv[i]); exit(1); }
    }
}

/* -dV/dphi_a where V = (mu/2)P^2/(1+kappa P^2), P = phi1 phi2 phi3 */
static inline double force_pot(double p1, double p2, double p3, int a)
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

/* Potential energy density at a point */
static inline double V_pot(double p1, double p2, double p3)
{
    double P = p1 * p2 * p3;
    double P2 = P * P;
    return 0.5 * mu * P2 / (1.0 + kappa * P2);
}

typedef struct {
    int survived;
    double omega;
    double A1_peak;
    double A3_peak;
    double E_total;
    double fc;
    double Q;          /* "charge" of field 3 if measurable */
    double phi3_far;   /* field 3 value far from core */
} RunResult;

static RunResult run_single(double m3_val, int save_ts, int save_profile)
{
    RunResult res = {0};

    double mass_sq[3];
    mass_sq[0] = m1 * m1;
    mass_sq[1] = m1 * m1;
    mass_sq[2] = m3_val * m3_val;

    double dx = 2.0 * xmax / (Nx - 1);
    double dx2 = dx * dx;

    /* CFL: use largest mass for stability */
    double m2max = mass_sq[0];
    if (mass_sq[1] > m2max) m2max = mass_sq[1];
    if (mass_sq[2] > m2max) m2max = mass_sq[2];
    double kmax = M_PI / dx;
    double dt = 0.8 * 2.0 / sqrt(kmax * kmax + m2max + 1.0);  /* +1 safety */
    int Nt = (int)(tfinal / dt) + 1;

    printf("mediator1d: m3=%.4f\n", m3_val);
    printf("  mu=%.1f kappa=%.1f m1=%.4f m3=%.4f A=%.3f sigma=%.3f\n",
           mu, kappa, m1, m3_val, A_init, sigma);
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

    /* Initialize: Gaussians for all three fields */
    for (int a = 0; a < 3; a++)
        for (int i = 0; i < Nx; i++) {
            double x = -xmax + i * dx;
            phi[a][i] = A_init * exp(-x * x / (2.0 * sigma * sigma));
        }

    int ic = Nx / 2;  /* center index */

    /* Compute acceleration macro (per-field mass) */
    #define COMPUTE_ACC() do { \
        for (int a = 0; a < 3; a++) { \
            acc[a][0] = acc[a][1] = acc[a][Nx-2] = acc[a][Nx-1] = 0; \
            for (int i = 1; i < Nx - 1; i++) { \
                double lapl = (phi[a][i+1] - 2.0*phi[a][i] + phi[a][i-1]) / dx2; \
                double fp = force_pot(phi[0][i], phi[1][i], phi[2][i], a); \
                acc[a][i] = lapl - mass_sq[a]*phi[a][i] + fp; \
            } \
        } \
    } while(0)

    COMPUTE_ACC();

    /* Time series output */
    FILE *fts = NULL;
    if (save_ts) {
        char tspath[600];
        snprintf(tspath, sizeof(tspath), "%s/mediator_m3_%.4f_ts.tsv", outdir, m3_val);
        fts = fopen(tspath, "w");
        if (!fts) { fprintf(stderr, "Cannot open %s\n", tspath); exit(1); }
        fprintf(fts, "time\tphi1_0\tphi2_0\tphi3_0\tpeak1\tpeak2\tpeak3\t"
                     "E_kin\tE_grad\tE_mass\tE_pot\tE_total\tf_core\n");
    }

    /* DFT storage for frequency measurement */
    int max_dft = 50000;
    double *phi0_hist = malloc(max_dft * sizeof(double));
    double *phi3_hist = malloc(max_dft * sizeof(double));
    double *t_hist = malloc(max_dft * sizeof(double));
    int n_dft = 0;

    int rec_every = Nt / 20000;
    if (rec_every < 1) rec_every = 1;
    int print_every = Nt / 40;
    if (print_every < 1) print_every = 1;
    int dft_every = Nt / max_dft;
    if (dft_every < 1) dft_every = 1;

    double core_r = 3.0 * sigma;

    double last_E = 0, last_fc = 0;
    double last_peak1 = 0, last_peak3 = 0;

    for (int n = 0; n <= Nt; n++) {
        double t = n * dt;

        if (n % dft_every == 0 && n_dft < max_dft) {
            phi0_hist[n_dft] = phi[0][ic];
            phi3_hist[n_dft] = phi[2][ic];
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
                    Em += 0.5 * mass_sq[a] * phi[a][i] * phi[a][i] * dx;
                    if (fabs(phi[a][i]) > peak[a]) peak[a] = fabs(phi[a][i]);
                }

                double V = V_pot(phi[0][i], phi[1][i], phi[2][i]);
                Ep += V * dx;

                double e = V;
                for (int a = 0; a < 3; a++) {
                    e += 0.5*vel[a][i]*vel[a][i];
                    double dp = (phi[a][i+1] - phi[a][i-1]) / (2.0*dx);
                    e += 0.5*dp*dp + 0.5*mass_sq[a]*phi[a][i]*phi[a][i];
                }
                Eall += e * dx;
                if (fabs(x) < core_r) Ecore += e * dx;
            }

            double Et = Ek + Eg + Em + Ep;
            double fc = (Eall > 1e-20) ? Ecore / Eall : 0.0;

            last_E = Et;
            last_fc = fc;
            last_peak1 = peak[0];
            last_peak3 = peak[2];

            if (fts && do_rec)
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

    if (fts) fclose(fts);

    /* Survival criterion: f_core > 0.1 and peak1 > 0.01 */
    res.survived = (last_fc > 0.1 && last_peak1 > 0.01) ? 1 : 0;
    res.A1_peak = last_peak1;
    res.A3_peak = last_peak3;
    res.E_total = last_E;
    res.fc = last_fc;

    /* DFT for frequency (field 1, second half) */
    int dft_start = n_dft / 2;
    double peak_pow = 0, peak_om = 0;
    if (n_dft - dft_start > 100) {
        double T = t_hist[n_dft-1] - t_hist[dft_start];
        int nf = 500;
        for (int k = 1; k < nf; k++) {
            double omega = 3.0 * m1 * k / nf;
            double re = 0, im = 0;
            for (int j = dft_start; j < n_dft; j++) {
                double dtj = (j > dft_start) ?
                    (t_hist[j]-t_hist[j-1]) : (t_hist[dft_start+1]-t_hist[dft_start]);
                re += phi0_hist[j] * cos(omega * t_hist[j]) * dtj;
                im += phi0_hist[j] * sin(omega * t_hist[j]) * dtj;
            }
            double pw = (re*re + im*im) / (T*T);
            if (pw > peak_pow) { peak_pow = pw; peak_om = omega; }
        }
        res.omega = peak_om;

        /* Also save spectrum */
        if (save_ts) {
            char dftpath[600];
            snprintf(dftpath, sizeof(dftpath), "%s/mediator_m3_%.4f_spectrum.tsv", outdir, m3_val);
            FILE *fdft = fopen(dftpath, "w");
            if (fdft) {
                fprintf(fdft, "omega\tpower\n");
                for (int k = 0; k < nf; k++) {
                    double omega = 3.0 * m1 * k / nf;
                    double re = 0, im = 0;
                    for (int j = dft_start; j < n_dft; j++) {
                        double dtj = (j > dft_start) ?
                            (t_hist[j]-t_hist[j-1]) : (t_hist[dft_start+1]-t_hist[dft_start]);
                        re += phi0_hist[j] * cos(omega * t_hist[j]) * dtj;
                        im += phi0_hist[j] * sin(omega * t_hist[j]) * dtj;
                    }
                    double pw = (re*re + im*im) / (T*T);
                    fprintf(fdft, "%.6f\t%.6e\n", omega, pw);
                }
                fclose(fdft);
            }
        }
    }

    /* Field 3 spatial profile at final time (for tail analysis) */
    if (save_profile) {
        char profpath[600];
        snprintf(profpath, sizeof(profpath), "%s/mediator_m3_%.4f_profile.tsv", outdir, m3_val);
        FILE *fprof = fopen(profpath, "w");
        if (fprof) {
            fprintf(fprof, "x\tphi1\tphi2\tphi3\n");
            for (int i = 0; i < Nx; i++) {
                double x = -xmax + i * dx;
                fprintf(fprof, "%.6f\t%.6e\t%.6e\t%.6e\n",
                        x, phi[0][i], phi[1][i], phi[2][i]);
            }
            fclose(fprof);
        }

        /* Measure "charge" Q from field 3 tail:
         * In 1D, massless Green's function: phi3 ~ Q_eff * |x| / 2 + const
         * So dphi3/dx = Q_eff/2 * sign(x) at large |x|
         * Measure slope at |x| = 30-50 (well outside core) */
        if (m3_val < 0.01 && res.survived) {
            int i30 = (int)((30.0 + xmax) / dx);
            int i50 = (int)((50.0 + xmax) / dx);
            if (i30 > 0 && i50 < Nx-1) {
                double slope = (phi[2][i50] - phi[2][i30]) / ((i50 - i30) * dx);
                res.Q = 2.0 * slope;  /* Q = 2 * dphi3/dx */
                res.phi3_far = phi[2][i50];
            }
        }
    }

    printf("  -> survived=%d omega=%.4f A1=%.4f A3=%.4f E=%.4f fc=%.3f Q=%.4e\n\n",
           res.survived, res.omega, res.A1_peak, res.A3_peak, res.E_total, res.fc, res.Q);

    for (int a = 0; a < 3; a++) { free(phi[a]); free(vel[a]); free(acc[a]); }
    free(damp); free(phi0_hist); free(phi3_hist); free(t_hist);

    return res;
}

int main(int argc, char **argv)
{
    parse_args(argc, argv);

    if (!do_scan) {
        /* Single run with specified m3 */
        run_single(m3, 1, 1);
    } else {
        /* Phase 2 scan */
        double m3_vals[] = {1.0, 0.8, 0.6, 0.4, 0.2, 0.1, 0.05, 0.0};
        int n_vals = 8;

        printf("============================================================\n");
        printf("Phase 2: m3 scan\n");
        printf("  mu=%.1f kappa=%.1f m1=%.4f\n", mu, kappa, m1);
        printf("  Nx=%d xmax=%.1f tfinal=%.0f\n", Nx, xmax, tfinal);
        printf("============================================================\n\n");

        RunResult results[8];
        for (int k = 0; k < n_vals; k++) {
            results[k] = run_single(m3_vals[k], 1, 1);
        }

        /* Summary table */
        printf("\n============================================================\n");
        printf("SUMMARY TABLE\n");
        printf("============================================================\n");
        printf("m3      survived  omega   A1_peak  A3_peak  E_total    fc     Q\n");
        printf("------  --------  ------  -------  -------  ---------  -----  ----------\n");
        for (int k = 0; k < n_vals; k++) {
            printf("%.4f  %s       %.4f  %.5f  %.5f  %+.4e  %.3f  %+.4e\n",
                   m3_vals[k],
                   results[k].survived ? "YES" : "NO ",
                   results[k].omega,
                   results[k].A1_peak,
                   results[k].A3_peak,
                   results[k].E_total,
                   results[k].fc,
                   results[k].Q);
        }
        printf("============================================================\n");

        /* Save summary to file */
        char sumpath[600];
        snprintf(sumpath, sizeof(sumpath), "%s/scan_summary.tsv", outdir);
        FILE *fsum = fopen(sumpath, "w");
        if (fsum) {
            fprintf(fsum, "m3\tsurvived\tomega\tA1_peak\tA3_peak\tE_total\tfc\tQ\n");
            for (int k = 0; k < n_vals; k++) {
                fprintf(fsum, "%.4f\t%d\t%.6f\t%.6f\t%.6f\t%.6e\t%.6f\t%.6e\n",
                        m3_vals[k], results[k].survived, results[k].omega,
                        results[k].A1_peak, results[k].A3_peak,
                        results[k].E_total, results[k].fc, results[k].Q);
            }
            fclose(fsum);
        }
    }

    return 0;
}
