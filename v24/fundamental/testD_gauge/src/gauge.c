/*
 * gauge.c — Test D: Local Breathing Gauge
 *
 * Three massive scalars φ_a with triple-product coupling (v21 oscillon)
 * plus a massless scalar field Ω(x,t) that modulates the effective mass:
 *
 *   m²_eff(x) = m² · (1 + g_Ω · Ω(x,t))
 *
 * Ω is sourced by the oscillon energy density:
 *   □Ω = ∂²Ω/∂t² - ∂²Ω/∂x² = g_source · ρ(x,t)
 *
 * Tests:
 *   1: Single oscillon + Ω coupling (measure Ω profile, ω shift)
 *   2: Two oscillons separated by D — does Ω mediate a force?
 *   3: Control: single oscillon, g_Ω=0 (Ω decoupled, reference ω)
 *
 * Compile: gcc -O3 -Wall -o gauge src/gauge.c -lm
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* --- Parameters --- */
static double mu      = -20.0;
static double kappa   = 20.0;
static double mass    = 1.0;
static double A_init  = 0.8;
static double sigma   = 3.0;
static double g_Omega = 0.1;
static double g_source = 0.01;
static int    Nx      = 4000;
static double xmax    = 100.0;
static double tfinal  = 10000.0;
static int    test    = 1;
static double sep     = 30.0;  /* separation for two-oscillon test */
static char   outdir[512] = "v24/fundamental/testD_gauge/data";

static void parse_args(int argc, char **argv)
{
    for (int i = 1; i < argc - 1; i += 2) {
        if      (!strcmp(argv[i], "-mu"))      mu      = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-kappa"))   kappa   = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-mass"))    mass    = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-A"))       A_init  = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-sigma"))   sigma   = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-gO"))      g_Omega = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-gs"))      g_source= atof(argv[i+1]);
        else if (!strcmp(argv[i], "-Nx"))      Nx      = atoi(argv[i+1]);
        else if (!strcmp(argv[i], "-xmax"))    xmax    = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-tfinal"))  tfinal  = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-test"))    test    = atoi(argv[i+1]);
        else if (!strcmp(argv[i], "-sep"))     sep     = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-o"))       strncpy(outdir, argv[i+1], sizeof(outdir)-1);
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

/* Energy density at grid point i (total: kinetic + gradient + mass + potential) */
static inline double energy_density(double *phi[3], double *vel[3],
                                     double *Om, double *vOm,
                                     double dx, double m2, int i)
{
    double e = 0.0;
    double m2_eff = m2 * (1.0 + g_Omega * Om[i]);

    for (int a = 0; a < 3; a++) {
        e += 0.5 * vel[a][i] * vel[a][i];
        double dp = (phi[a][i+1] - phi[a][i-1]) / (2.0 * dx);
        e += 0.5 * dp * dp;
        e += 0.5 * m2_eff * phi[a][i] * phi[a][i];
    }

    double P = phi[0][i] * phi[1][i] * phi[2][i];
    double P2 = P * P;
    e += 0.5 * mu * P2 / (1.0 + kappa * P2);

    /* Ω kinetic + gradient energy */
    e += 0.5 * vOm[i] * vOm[i];
    double dOm = (Om[i+1] - Om[i-1]) / (2.0 * dx);
    e += 0.5 * dOm * dOm;

    return e;
}

/* Source energy density for Ω (just the matter sector, no Ω self-energy) */
static inline double matter_density(double *phi[3], double *vel[3],
                                     double dx, double m2, int i)
{
    double e = 0.0;
    for (int a = 0; a < 3; a++) {
        e += 0.5 * vel[a][i] * vel[a][i];
        double dp = (phi[a][i+1] - phi[a][i-1]) / (2.0 * dx);
        e += 0.5 * dp * dp;
        e += 0.5 * m2 * phi[a][i] * phi[a][i];
    }
    double P = phi[0][i] * phi[1][i] * phi[2][i];
    double P2 = P * P;
    e += 0.5 * mu * P2 / (1.0 + kappa * P2);
    return e;
}

int main(int argc, char **argv)
{
    parse_args(argc, argv);

    double dx  = 2.0 * xmax / (Nx - 1);
    double dx2 = dx * dx;
    double m2  = mass * mass;

    /* CFL: c=1 for massless Ω, but also need stability for massive φ */
    double dt = 0.4 * dx;  /* CFL < 1 for wave equation */
    int Nt = (int)(tfinal / dt) + 1;

    const char *test_desc;
    switch (test) {
        case 1: test_desc = "Single oscillon + Omega"; break;
        case 2: test_desc = "Two oscillons (force test)"; break;
        case 3: test_desc = "Control: g_Omega=0"; break;
        default: fprintf(stderr, "Unknown test %d\n", test); return 1;
    }

    double g_Om_eff = (test == 3) ? 0.0 : g_Omega;

    printf("gauge test %d: %s\n", test, test_desc);
    printf("  mu=%.1f kappa=%.1f mass=%.3f A=%.3f sigma=%.3f\n",
           mu, kappa, mass, A_init, sigma);
    printf("  g_Omega=%.4f g_source=%.4f\n", g_Om_eff, g_source);
    printf("  Nx=%d xmax=%.1f dx=%.5f dt=%.6f Nt=%d tfinal=%.0f\n",
           Nx, xmax, dx, dt, Nt, tfinal);
    if (test == 2) printf("  separation=%.1f\n", sep);

    /* Allocate fields */
    double *phi[3], *vel[3], *acc[3];
    for (int a = 0; a < 3; a++) {
        phi[a] = calloc(Nx, sizeof(double));
        vel[a] = calloc(Nx, sizeof(double));
        acc[a] = calloc(Nx, sizeof(double));
    }
    double *Om   = calloc(Nx, sizeof(double));  /* Ω field */
    double *vOm  = calloc(Nx, sizeof(double));  /* Ω velocity */
    double *aOm  = calloc(Nx, sizeof(double));  /* Ω acceleration */

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
    if (test == 2) {
        /* Two oscillons at ±sep/2 */
        for (int a = 0; a < 3; a++)
            for (int i = 0; i < Nx; i++) {
                double x = -xmax + i * dx;
                double x1 = x - sep / 2.0;
                double x2 = x + sep / 2.0;
                phi[a][i] = A_init * exp(-x1*x1/(2.0*sigma*sigma))
                          + A_init * exp(-x2*x2/(2.0*sigma*sigma));
            }
    } else {
        /* Single oscillon at center */
        for (int a = 0; a < 3; a++)
            for (int i = 0; i < Nx; i++) {
                double x = -xmax + i * dx;
                phi[a][i] = A_init * exp(-x*x/(2.0*sigma*sigma));
            }
    }

    /* Ω starts at zero */

    /* Macro: compute accelerations for φ_a and Ω */
    #define COMPUTE_ACC() do { \
        double m2_loc; \
        for (int a = 0; a < 3; a++) { \
            acc[a][0] = acc[a][1] = acc[a][Nx-2] = acc[a][Nx-1] = 0; \
            for (int ii = 1; ii < Nx - 1; ii++) { \
                double lapl = (phi[a][ii+1] - 2.0*phi[a][ii] + phi[a][ii-1]) / dx2; \
                m2_loc = m2 * (1.0 + g_Om_eff * Om[ii]); \
                double fp = force_pot(phi[0][ii], phi[1][ii], phi[2][ii], a); \
                acc[a][ii] = lapl - m2_loc * phi[a][ii] + fp; \
            } \
        } \
        /* Ω: □Ω = g_source * ρ_matter */ \
        aOm[0] = aOm[1] = aOm[Nx-2] = aOm[Nx-1] = 0; \
        for (int ii = 1; ii < Nx - 1; ii++) { \
            double laplO = (Om[ii+1] - 2.0*Om[ii] + Om[ii-1]) / dx2; \
            double rho = matter_density(phi, vel, dx, m2, ii); \
            aOm[ii] = laplO + g_source * rho; \
        } \
    } while(0)

    COMPUTE_ACC();

    /* Output files */
    char tspath[600], profpath[600];
    snprintf(tspath, sizeof(tspath), "%s/gauge_test%d_ts.tsv", outdir, test);
    FILE *fts = fopen(tspath, "w");
    if (!fts) { fprintf(stderr, "Cannot open %s\n", tspath); return 1; }
    if (test == 2) {
        fprintf(fts, "time\tphi1_L\tphi1_R\tOm_L\tOm_R\tOm_mid\t"
                     "E_total\tE_core_L\tE_core_R\txcm_L\txcm_R\tsep\n");
    } else {
        fprintf(fts, "time\tphi1_0\tOm_0\tOm_10\tOm_20\tOm_40\t"
                     "E_total\tE_core\tf_core\tpeak_phi\tpeak_Om\n");
    }

    /* DFT storage */
    int max_dft = 50000;
    double *phi0_hist = malloc(max_dft * sizeof(double));
    double *t_hist    = malloc(max_dft * sizeof(double));
    int n_dft = 0;

    int rec_every   = Nt / 20000; if (rec_every < 1) rec_every = 1;
    int print_every = Nt / 40;    if (print_every < 1) print_every = 1;
    int dft_every   = Nt / max_dft; if (dft_every < 1) dft_every = 1;
    int prof_every  = Nt / 10;    if (prof_every < 1) prof_every = 1;

    double core_r = 3.0 * sigma;
    int ic = Nx / 2;

    /* Grid indices for Ω samples at specific distances */
    int idx_10 = ic + (int)(10.0 / dx);
    int idx_20 = ic + (int)(20.0 / dx);
    int idx_40 = ic + (int)(40.0 / dx);
    if (idx_10 >= Nx) idx_10 = Nx - 2;
    if (idx_20 >= Nx) idx_20 = Nx - 2;
    if (idx_40 >= Nx) idx_40 = Nx - 2;

    /* For two-oscillon test */
    int ic_L = ic - (int)(sep / (2.0 * dx));
    int ic_R = ic + (int)(sep / (2.0 * dx));
    if (ic_L < 1) ic_L = 1;
    if (ic_R >= Nx - 1) ic_R = Nx - 2;

    int prof_count = 0;

    printf("\nStarting evolution...\n");

    for (int n = 0; n <= Nt; n++) {
        double t = n * dt;

        /* DFT sampling */
        if (n % dft_every == 0 && n_dft < max_dft) {
            phi0_hist[n_dft] = phi[0][ic];
            t_hist[n_dft] = t;
            n_dft++;
        }

        int do_rec   = (n % rec_every == 0);
        int do_print = (n % print_every == 0);
        int do_prof  = (n % prof_every == 0);

        if (do_rec || do_print || do_prof) {
            if (test == 2) {
                /* Two-oscillon diagnostics */
                double Et = 0, EcL = 0, EcR = 0;
                double xcm_L = 0, wt_L = 0, xcm_R = 0, wt_R = 0;
                double xL = -sep / 2.0, xR = sep / 2.0;

                for (int i = 2; i < Nx - 2; i++) {
                    double x = -xmax + i * dx;
                    double e = energy_density(phi, vel, Om, vOm, dx, m2, i);
                    Et += e * dx;
                    if (fabs(x - xL) < core_r) EcL += e * dx;
                    if (fabs(x - xR) < core_r) EcR += e * dx;
                    /* Centroid: left half x<0, right half x>0 */
                    if (x < 0) { xcm_L += x * e * dx; wt_L += e * dx; }
                    else       { xcm_R += x * e * dx; wt_R += e * dx; }
                }
                if (wt_L > 1e-20) xcm_L /= wt_L;
                if (wt_R > 1e-20) xcm_R /= wt_R;
                double sep_meas = xcm_R - xcm_L;

                if (do_rec)
                    fprintf(fts, "%.6f\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6f\t%.6f\t%.6f\n",
                            t, phi[0][ic_L], phi[0][ic_R],
                            Om[ic_L], Om[ic_R], Om[ic],
                            Et, EcL, EcR, xcm_L, xcm_R, sep_meas);
                if (do_print)
                    printf("  t=%7.1f  phi_L=%+.3f phi_R=%+.3f  Om_L=%.3e Om_R=%.3e Om_mid=%.3e  "
                           "E=%.3f  sep=%.3f\n",
                           t, phi[0][ic_L], phi[0][ic_R],
                           Om[ic_L], Om[ic_R], Om[ic],
                           Et, sep_meas);
            } else {
                /* Single oscillon diagnostics */
                double Et = 0, Ecore = 0;
                double peak_phi = 0, peak_Om = 0;
                for (int i = 2; i < Nx - 2; i++) {
                    double x = -xmax + i * dx;
                    double e = energy_density(phi, vel, Om, vOm, dx, m2, i);
                    Et += e * dx;
                    if (fabs(x) < core_r) Ecore += e * dx;
                    for (int a = 0; a < 3; a++)
                        if (fabs(phi[a][i]) > peak_phi) peak_phi = fabs(phi[a][i]);
                    if (fabs(Om[i]) > peak_Om) peak_Om = fabs(Om[i]);
                }
                double fc = (Et > 1e-20) ? Ecore / Et : 0.0;

                if (do_rec)
                    fprintf(fts, "%.6f\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6f\t%.6e\t%.6e\n",
                            t, phi[0][ic], Om[ic],
                            (idx_10 < Nx) ? Om[idx_10] : 0.0,
                            (idx_20 < Nx) ? Om[idx_20] : 0.0,
                            (idx_40 < Nx) ? Om[idx_40] : 0.0,
                            Et, Ecore, fc, peak_phi, peak_Om);
                if (do_print)
                    printf("  t=%7.1f  phi0=%+.4f  Om(0)=%+.3e Om(10)=%+.3e Om(20)=%+.3e Om(40)=%+.3e  "
                           "E=%.4f fc=%.3f pk_phi=%.3f\n",
                           t, phi[0][ic], Om[ic], Om[idx_10], Om[idx_20], Om[idx_40],
                           Et, fc, peak_phi);
            }

            /* Snapshot Ω profile at selected times */
            if (do_prof) {
                snprintf(profpath, sizeof(profpath), "%s/gauge_test%d_profile_t%d.tsv",
                         outdir, test, (int)t);
                FILE *fp = fopen(profpath, "w");
                if (fp) {
                    fprintf(fp, "x\tphi1\tphi2\tphi3\tOmega\trho\n");
                    int stride = 4;
                    for (int i = 2; i < Nx - 2; i += stride) {
                        double x = -xmax + i * dx;
                        double rho = matter_density(phi, vel, dx, m2, i);
                        fprintf(fp, "%.6f\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\n",
                                x, phi[0][i], phi[1][i], phi[2][i], Om[i], rho);
                    }
                    fclose(fp);
                    prof_count++;
                }
            }
        }

        if (n == Nt) break;

        /* Velocity Verlet: half-kick */
        for (int a = 0; a < 3; a++)
            for (int i = 1; i < Nx - 1; i++)
                vel[a][i] += 0.5 * dt * acc[a][i];
        for (int i = 1; i < Nx - 1; i++)
            vOm[i] += 0.5 * dt * aOm[i];

        /* Drift */
        for (int a = 0; a < 3; a++)
            for (int i = 1; i < Nx - 1; i++)
                phi[a][i] += dt * vel[a][i];
        for (int i = 1; i < Nx - 1; i++)
            Om[i] += dt * vOm[i];

        /* Recompute accelerations */
        COMPUTE_ACC();

        /* Half-kick */
        for (int a = 0; a < 3; a++)
            for (int i = 1; i < Nx - 1; i++)
                vel[a][i] += 0.5 * dt * acc[a][i];
        for (int i = 1; i < Nx - 1; i++)
            vOm[i] += 0.5 * dt * aOm[i];

        /* Absorbing boundary */
        for (int a = 0; a < 3; a++)
            for (int i = 0; i < Nx; i++) {
                vel[a][i] *= damp[i];
                phi[a][i] *= damp[i];
            }
        for (int i = 0; i < Nx; i++) {
            vOm[i] *= damp[i];
            Om[i]  *= damp[i];
        }
    }

    fclose(fts);

    /* DFT of phi_1(0,t) — second half */
    int dft_start = n_dft / 2;
    double peak_om = 0, peak_pow = 0;
    if (n_dft - dft_start > 100) {
        char dftpath[600];
        snprintf(dftpath, sizeof(dftpath), "%s/gauge_test%d_spectrum.tsv", outdir, test);
        FILE *fdft = fopen(dftpath, "w");
        fprintf(fdft, "omega\tpower\n");
        double T = t_hist[n_dft-1] - t_hist[dft_start];
        int nf = 500;
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

    /* DFT of Ω(0,t) — to check if Ω is static or oscillating */
    if (test != 2) {
        /* Collect Ω time history at center from the profile snapshots is hard;
           instead do DFT of Ω at center using same sampling */
        /* Re-run a lightweight DFT on the time series data */
        /* We'll use the recorded Om_0 values from fts — but easier to just
           note Ω behavior from the time series output */
        printf("\nOmega analysis:\n");
        printf("  Ω(0) at end: check time series file\n");
    }

    printf("\nProfile snapshots written: %d\n", prof_count);
    printf("Output: %s\n", tspath);

    /* Cleanup */
    for (int a = 0; a < 3; a++) { free(phi[a]); free(vel[a]); free(acc[a]); }
    free(Om); free(vOm); free(aOm);
    free(damp); free(phi0_hist); free(t_hist);
    return 0;
}
