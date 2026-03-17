/*
 * uud180.c — UUD/UDD oscillon in the 180-degree anti-phase state
 *
 * Three massive scalars with triple-product coupling, per-field masses.
 * UUD (proton): phi1=+A, phi2=+A, phi3=-A (two up, one down)
 * UDD (neutron): phi1=+A, phi2=-A, phi3=-A (one up, two down)
 *
 * Lagrangian:
 *   L = sum_a [ (1/2)(dt phi_a)^2 - (1/2)(dx phi_a)^2 - (m_a^2/2)phi_a^2 ]
 *     - (mu/2) P^2 / (1 + kappa P^2)
 *   P = phi_1 phi_2 phi_3
 *
 * Scans m_D from 1.0 to 0.70 for both UUD and UDD configurations.
 * Outputs per-run time series and a summary mass_ordering.tsv.
 *
 * Compile: gcc -O3 -Wall -o uud180 src/uud180.c -lm
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static double mu     = -20.0;
static double kappa  = 20.0;
static double m_U    = 1.0;
static double sigma  = 3.0;
static double A_init = 0.8;
static int    Nx     = 4000;
static double xmax   = 100.0;
static double tfinal = 10000.0;
static char   outdir[512] = "data";

/* Per-field masses (squared) and signs */
static double m2[3];
static double signs[3];

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

typedef struct {
    double E_final;
    double E_pot_final;
    double fc_final;
    double peak_final[3];
    double omega_peak;
    int    stable;
} RunResult;

static RunResult run_one(const char *label, double m_D,
                         double s1, double s2, double s3)
{
    RunResult res;
    memset(&res, 0, sizeof(res));

    double dx = 2.0 * xmax / (Nx - 1);
    double dx2 = dx * dx;

    /* Per-field masses */
    m2[0] = (s1 > 0 ? m_U : m_D) * (s1 > 0 ? m_U : m_D);
    m2[1] = (s2 > 0 ? m_U : m_D) * (s2 > 0 ? m_U : m_D);
    m2[2] = (s3 > 0 ? m_U : m_D) * (s3 > 0 ? m_U : m_D);
    signs[0] = s1; signs[1] = s2; signs[2] = s3;

    /* CFL: use minimum mass for stability */
    double m_min = m_D < m_U ? m_D : m_U;
    double kmax = M_PI / dx;
    double dt = 0.8 * 2.0 / sqrt(kmax*kmax + m_min*m_min);
    int Nt = (int)(tfinal / dt) + 1;

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

    /* Initialize: signed Gaussians */
    for (int a = 0; a < 3; a++)
        for (int i = 0; i < Nx; i++) {
            double x = -xmax + i * dx;
            phi[a][i] = signs[a] * A_init * exp(-x * x / (2.0 * sigma * sigma));
        }

    /* Compute acceleration macro */
    #define COMPUTE_ACC() do { \
        for (int a = 0; a < 3; a++) { \
            acc[a][0] = acc[a][1] = acc[a][Nx-2] = acc[a][Nx-1] = 0; \
            for (int i = 1; i < Nx - 1; i++) { \
                double lapl = (phi[a][i+1] - 2.0*phi[a][i] + phi[a][i-1]) / dx2; \
                double fp = force_pot(phi[0][i], phi[1][i], phi[2][i], a); \
                acc[a][i] = lapl - m2[a]*phi[a][i] + fp; \
            } \
        } \
    } while(0)

    COMPUTE_ACC();

    /* Output file */
    char tspath[600];
    snprintf(tspath, sizeof(tspath), "%s/%s_mD%.2f_ts.tsv", outdir, label, m_D);
    FILE *fts = fopen(tspath, "w");
    if (!fts) { fprintf(stderr, "Cannot open %s\n", tspath); exit(1); }
    fprintf(fts, "time\tphi1_0\tphi2_0\tphi3_0\tpeak1\tpeak2\tpeak3\t"
                 "E_kin\tE_grad\tE_mass\tE_pot\tE_total\tf_core\n");

    /* DFT storage for all 3 fields */
    int max_dft = 50000;
    double *phi_hist[3];
    for (int a = 0; a < 3; a++)
        phi_hist[a] = malloc(max_dft * sizeof(double));
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

    double last_E = 0, last_fc = 0;
    double last_peak[3] = {0};
    double last_Ep = 0;

    for (int n = 0; n <= Nt; n++) {
        double t = n * dt;

        if (n % dft_every == 0 && n_dft < max_dft) {
            for (int a = 0; a < 3; a++)
                phi_hist[a][n_dft] = phi[a][ic];
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
                    Em += 0.5 * m2[a] * phi[a][i] * phi[a][i] * dx;
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
                    e += 0.5*dp*dp + 0.5*m2[a]*phi[a][i]*phi[a][i];
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
                printf("  [%s mD=%.2f] t=%7.1f  p0=(%+.3f,%+.3f,%+.3f)  "
                       "pk=(%.3f,%.3f,%.3f)  E=%+.4f  Ep=%+.4f  fc=%.3f\n",
                       label, m_D, t,
                       phi[0][ic], phi[1][ic], phi[2][ic],
                       peak[0], peak[1], peak[2], Et, Ep, fc);

            last_E = Et; last_Ep = Ep; last_fc = fc;
            for (int a = 0; a < 3; a++) last_peak[a] = peak[a];
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

    /* DFT of each field — second half */
    int dft_start = n_dft / 2;
    double omega_peak_all = 0;
    if (n_dft - dft_start > 100) {
        char dftpath[600];
        snprintf(dftpath, sizeof(dftpath), "%s/%s_mD%.2f_spectrum.tsv",
                 outdir, label, m_D);
        FILE *fdft = fopen(dftpath, "w");
        fprintf(fdft, "omega\tpower1\tpower2\tpower3\n");
        double T = t_hist[n_dft-1] - t_hist[dft_start];
        int nf = 500;
        double peak_pow = 0;
        for (int k = 0; k < nf; k++) {
            double omega = 3.0 * m_U / nf * k;  /* scan up to 3*m_U */
            double pw_tot = 0;
            fprintf(fdft, "%.6f", omega);
            for (int a = 0; a < 3; a++) {
                double re = 0, im = 0;
                for (int j = dft_start; j < n_dft; j++) {
                    double dtj = (j > dft_start) ?
                        (t_hist[j]-t_hist[j-1]) : (t_hist[dft_start+1]-t_hist[dft_start]);
                    re += phi_hist[a][j] * cos(omega * t_hist[j]) * dtj;
                    im += phi_hist[a][j] * sin(omega * t_hist[j]) * dtj;
                }
                double pw = (re*re + im*im) / (T*T);
                fprintf(fdft, "\t%.6e", pw);
                pw_tot += pw;
            }
            fprintf(fdft, "\n");
            if (pw_tot > peak_pow) { peak_pow = pw_tot; omega_peak_all = omega; }
        }
        fclose(fdft);
    }

    res.E_final = last_E;
    res.E_pot_final = last_Ep;
    res.fc_final = last_fc;
    for (int a = 0; a < 3; a++) res.peak_final[a] = last_peak[a];
    res.omega_peak = omega_peak_all;
    res.stable = (last_fc > 0.3) ? 1 : 0;

    printf("  [%s mD=%.2f] FINAL: E=%.4f  Ep=%.4f  fc=%.3f  omega=%.4f  %s\n\n",
           label, m_D, res.E_final, res.E_pot_final, res.fc_final,
           res.omega_peak, res.stable ? "STABLE" : "DISPERSED");

    for (int a = 0; a < 3; a++) { free(phi[a]); free(vel[a]); free(acc[a]); free(phi_hist[a]); }
    free(damp); free(t_hist);

    return res;
    #undef COMPUTE_ACC
}

int main(int argc, char **argv)
{
    /* Parse optional overrides */
    for (int i = 1; i < argc - 1; i += 2) {
        if      (!strcmp(argv[i], "-mu"))     mu     = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-kappa"))  kappa  = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-A"))      A_init = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-sigma"))  sigma  = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-Nx"))     Nx     = atoi(argv[i+1]);
        else if (!strcmp(argv[i], "-xmax"))   xmax   = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-tfinal")) tfinal = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-o"))      strncpy(outdir, argv[i+1], sizeof(outdir)-1);
        else { fprintf(stderr, "Unknown flag: %s\n", argv[i]); exit(1); }
    }

    printf("=== V24-180B: UUD/UDD in the 180-degree Anti-Phase State ===\n");
    printf("  mu=%.1f  kappa=%.1f  m_U=%.2f  A=%.2f  sigma=%.1f\n",
           mu, kappa, m_U, A_init, sigma);
    printf("  Nx=%d  xmax=%.0f  tfinal=%.0f\n\n", Nx, xmax, tfinal);

    double mD_vals[] = {1.00, 0.95, 0.90, 0.85, 0.80, 0.70};
    int n_mD = sizeof(mD_vals) / sizeof(mD_vals[0]);

    /* Summary file */
    char sumpath[600];
    snprintf(sumpath, sizeof(sumpath), "%s/mass_ordering.tsv", outdir);
    FILE *fsum = fopen(sumpath, "w");
    if (!fsum) { fprintf(stderr, "Cannot open %s\n", sumpath); exit(1); }
    fprintf(fsum, "m_D\tE_UUD\tE_UDD\tE_UDD-E_UUD\tfc_UUD\tfc_UDD\t"
                  "omega_UUD\tomega_UDD\tstable_UUD\tstable_UDD\n");

    for (int im = 0; im < n_mD; im++) {
        double mD = mD_vals[im];
        printf("========== m_D = %.2f ==========\n\n", mD);

        /* UUD: (+,+,-) — two up, one down */
        RunResult uud = run_one("uud180", mD, +1.0, +1.0, -1.0);

        /* UDD: (+,-,-) — one up, two down */
        RunResult udd = run_one("udd180", mD, +1.0, -1.0, -1.0);

        double dE = udd.E_final - uud.E_final;
        printf("  >>> m_D=%.2f: E_UUD=%.4f  E_UDD=%.4f  delta=%.4f  %s\n\n",
               mD, uud.E_final, udd.E_final, dE,
               dE > 0 ? "UDD > UUD (neutron heavier)" : "UUD > UDD (proton heavier)");

        fprintf(fsum, "%.2f\t%.6f\t%.6f\t%.6f\t%.4f\t%.4f\t%.4f\t%.4f\t%d\t%d\n",
                mD, uud.E_final, udd.E_final, dE,
                uud.fc_final, udd.fc_final,
                uud.omega_peak, udd.omega_peak,
                uud.stable, udd.stable);
    }

    fclose(fsum);
    printf("Summary written to %s/mass_ordering.tsv\n", outdir);

    return 0;
}
