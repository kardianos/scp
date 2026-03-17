/*
 * rad180.c — Radiation comparison: 0-degree vs 180-degree oscillon
 *
 * Runs BOTH configurations simultaneously and compares:
 *   0-deg:  phi1=phi2=phi3 = +A*gaussian  (P = +f^3 cos^3 wt)
 *   180-deg: phi1=phi2=+A, phi3=-A         (P = -f^3 cos^3 wt)
 *
 * Also tests ASYMMETRIC 180-deg: phi1,phi2 at A+, phi3 at -A-
 *
 * Lagrangian:
 *   L = sum_a [ (1/2)(dt phi_a)^2 - (1/2)(dx phi_a)^2 - (m^2/2)phi_a^2 ]
 *     - (mu/2) P^2 / (1 + kappa P^2)
 *   P = phi_1 phi_2 phi_3
 *
 * Compile: gcc -O3 -Wall -o rad180 src/rad180.c -lm
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ===================== Parameters ===================== */
static double mu     = -20.0;
static double kappa  = 20.0;
static double mass   = 1.0;
static double A_init = 0.8;
static double sigma  = 3.0;
static int    Nx     = 4000;
static double xmax   = 100.0;
static double tfinal = 20000.0;
static char   outdir[512] = "data";

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

/* ===================== Oscillon State ===================== */
typedef struct {
    double *phi[3];
    double *vel[3];
    double *acc[3];
} State;

static void state_alloc(State *s, int n)
{
    for (int a = 0; a < 3; a++) {
        s->phi[a] = calloc(n, sizeof(double));
        s->vel[a] = calloc(n, sizeof(double));
        s->acc[a] = calloc(n, sizeof(double));
    }
}

static void state_free(State *s)
{
    for (int a = 0; a < 3; a++) {
        free(s->phi[a]); free(s->vel[a]); free(s->acc[a]);
    }
}

static void compute_acc(State *s, double dx, double m2)
{
    double dx2 = dx * dx;
    for (int a = 0; a < 3; a++) {
        s->acc[a][0] = s->acc[a][1] = s->acc[a][Nx-2] = s->acc[a][Nx-1] = 0;
        for (int i = 1; i < Nx - 1; i++) {
            double lapl = (s->phi[a][i+1] - 2.0*s->phi[a][i] + s->phi[a][i-1]) / dx2;
            double fp = force_pot(s->phi[0][i], s->phi[1][i], s->phi[2][i], a);
            s->acc[a][i] = lapl - m2 * s->phi[a][i] + fp;
        }
    }
}

static void step_vv(State *s, double dt, double dx, double m2, double *damp)
{
    for (int a = 0; a < 3; a++)
        for (int i = 1; i < Nx - 1; i++)
            s->vel[a][i] += 0.5 * dt * s->acc[a][i];
    for (int a = 0; a < 3; a++)
        for (int i = 1; i < Nx - 1; i++)
            s->phi[a][i] += dt * s->vel[a][i];
    compute_acc(s, dx, m2);
    for (int a = 0; a < 3; a++)
        for (int i = 1; i < Nx - 1; i++)
            s->vel[a][i] += 0.5 * dt * s->acc[a][i];
    /* absorbing boundary */
    for (int a = 0; a < 3; a++)
        for (int i = 0; i < Nx; i++) {
            s->vel[a][i] *= damp[i];
            s->phi[a][i] *= damp[i];
        }
}

/* Compute energies. Returns total energy, fills breakdown. */
static double compute_energy(State *s, double dx, double m2, double core_r,
                             double *Ek_out, double *Eg_out, double *Em_out,
                             double *Ep_out, double *Ecore_out)
{
    double Ek = 0, Eg = 0, Em = 0, Ep = 0, Ecore = 0, Eall = 0;
    for (int i = 1; i < Nx - 1; i++) {
        double x = -xmax + i * dx;
        double eloc = 0;
        for (int a = 0; a < 3; a++) {
            double v2 = 0.5 * s->vel[a][i] * s->vel[a][i];
            double dp = (s->phi[a][i+1] - s->phi[a][i-1]) / (2.0*dx);
            double g2 = 0.5 * dp * dp;
            double m2e = 0.5 * m2 * s->phi[a][i] * s->phi[a][i];
            Ek += v2 * dx;
            Eg += g2 * dx;
            Em += m2e * dx;
            eloc += v2 + g2 + m2e;
        }
        double P = s->phi[0][i] * s->phi[1][i] * s->phi[2][i];
        double P2 = P * P;
        double V = 0.5 * mu * P2 / (1.0 + kappa * P2);
        Ep += V * dx;
        eloc += V;

        Eall += eloc * dx;
        if (fabs(x) < core_r) Ecore += eloc * dx;
    }
    if (Ek_out) *Ek_out = Ek;
    if (Eg_out) *Eg_out = Eg;
    if (Em_out) *Em_out = Em;
    if (Ep_out) *Ep_out = Ep;
    if (Ecore_out) *Ecore_out = Ecore;
    return Ek + Eg + Em + Ep;
}

/* ===================== DFT of time series ===================== */
static void do_dft(const char *fname, double *data, double *times, int n_pts,
                   int start, double m, int n_freq)
{
    FILE *f = fopen(fname, "w");
    if (!f) { fprintf(stderr, "Cannot open %s\n", fname); return; }
    fprintf(f, "omega\tpower\n");
    double T = times[n_pts-1] - times[start];
    double peak_pow = 0, peak_om = 0;
    for (int k = 0; k < n_freq; k++) {
        double omega = 3.0 * m * k / n_freq;
        double re = 0, im = 0;
        for (int j = start; j < n_pts; j++) {
            double dtj = (j > start) ? (times[j]-times[j-1]) :
                         (times[start+1]-times[start]);
            re += data[j] * cos(omega * times[j]) * dtj;
            im += data[j] * sin(omega * times[j]) * dtj;
        }
        double pw = (re*re + im*im) / (T*T);
        fprintf(f, "%.6f\t%.6e\n", omega, pw);
        if (pw > peak_pow) { peak_pow = pw; peak_om = omega; }
    }
    fclose(f);
    printf("  DFT %s: peak omega = %.4f\n", fname, peak_om);
}

/* ===================== Run one configuration ===================== */
typedef struct {
    double Ainit[3];       /* initial amplitudes (sign included) */
    const char *label;
    const char *ts_file;
    const char *dft_file;
} Config;

typedef struct {
    /* Late-time dE/dt measurement */
    double E_at_t1;        /* E at t=10000 */
    double E_at_t2;        /* E at t=tfinal */
    double dEdt_late;      /* (E2-E1)/(t2-t1) */
    double Ecore_at_t2;    /* core fraction at end */
    double peak_omega;     /* fundamental frequency */
    double peak_amp[3];    /* field amplitudes at end */
    double E_total_final;
} RunResult;

static void run_config(Config *cfg, double dx, double dt, int Nt, double *damp,
                       double m2, RunResult *res)
{
    State s;
    state_alloc(&s, Nx);

    int ic = Nx / 2;
    double core_r = 3.0 * sigma;

    /* Initialize Gaussians */
    for (int a = 0; a < 3; a++)
        for (int i = 0; i < Nx; i++) {
            double x = -xmax + i * dx;
            s.phi[a][i] = cfg->Ainit[a] * exp(-x*x / (2.0*sigma*sigma));
        }

    compute_acc(&s, dx, m2);

    /* Time series file */
    char tspath[600];
    snprintf(tspath, sizeof(tspath), "%s/%s", outdir, cfg->ts_file);
    FILE *fts = fopen(tspath, "w");
    fprintf(fts, "time\tphi1_0\tphi2_0\tphi3_0\tpeak1\tpeak2\tpeak3\t"
                 "E_kin\tE_grad\tE_mass\tE_pot\tE_total\tE_core\tf_core\n");

    /* DFT history */
    int max_dft = 50000;
    double *phi1_hist = malloc(max_dft * sizeof(double));
    double *phi3_hist = malloc(max_dft * sizeof(double));
    double *t_hist = malloc(max_dft * sizeof(double));
    int n_dft = 0;

    int rec_every = Nt / 20000;
    if (rec_every < 1) rec_every = 1;
    int print_every = Nt / 20;
    if (print_every < 1) print_every = 1;
    int dft_every = Nt / max_dft;
    if (dft_every < 1) dft_every = 1;

    double E_t1 = 0;
    int got_t1 = 0;

    printf("\n=== Running: %s ===\n", cfg->label);
    printf("  Init amplitudes: (%.3f, %.3f, %.3f)\n",
           cfg->Ainit[0], cfg->Ainit[1], cfg->Ainit[2]);

    for (int n = 0; n <= Nt; n++) {
        double t = n * dt;

        if (n % dft_every == 0 && n_dft < max_dft) {
            phi1_hist[n_dft] = s.phi[0][ic];
            phi3_hist[n_dft] = s.phi[2][ic];
            t_hist[n_dft] = t;
            n_dft++;
        }

        int do_rec = (n % rec_every == 0);
        int do_print = (n % print_every == 0);

        if (do_rec || do_print) {
            double Ek, Eg, Em, Ep, Ecore;
            double Et = compute_energy(&s, dx, m2, core_r, &Ek, &Eg, &Em, &Ep, &Ecore);
            double fc = (Et > 1e-20) ? Ecore / Et : 0.0;

            double peak[3] = {0};
            for (int i = 1; i < Nx - 1; i++)
                for (int a = 0; a < 3; a++)
                    if (fabs(s.phi[a][i]) > peak[a]) peak[a] = fabs(s.phi[a][i]);

            /* Record E at t ~ 10000 */
            if (!got_t1 && t >= 10000.0) {
                E_t1 = Et;
                got_t1 = 1;
            }

            if (do_rec)
                fprintf(fts, "%.6f\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t"
                             "%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6f\n",
                        t, s.phi[0][ic], s.phi[1][ic], s.phi[2][ic],
                        peak[0], peak[1], peak[2], Ek, Eg, Em, Ep, Et, Ecore, fc);
            if (do_print)
                printf("  t=%7.0f  p0=(%+.4f,%+.4f,%+.4f)  pk=(%.4f,%.4f,%.4f)  "
                       "E=%.4f  fc=%.3f\n",
                       t, s.phi[0][ic], s.phi[1][ic], s.phi[2][ic],
                       peak[0], peak[1], peak[2], Et, fc);

            if (n == Nt) {
                res->E_at_t2 = Et;
                res->Ecore_at_t2 = fc;
                res->E_total_final = Et;
                for (int a = 0; a < 3; a++) res->peak_amp[a] = peak[a];
            }
        }

        if (n == Nt) break;
        step_vv(&s, dt, dx, m2, damp);
    }

    fclose(fts);

    res->E_at_t1 = E_t1;
    res->dEdt_late = (res->E_at_t2 - E_t1) / (tfinal - 10000.0);

    /* DFT of phi1 and phi3 — second half only */
    int dft_start = n_dft / 2;
    if (n_dft - dft_start > 100) {
        char dftpath[600];
        snprintf(dftpath, sizeof(dftpath), "%s/%s", outdir, cfg->dft_file);
        do_dft(dftpath, phi1_hist, t_hist, n_dft, dft_start, mass, 500);

        /* Also DFT of phi3 for asymmetric cases */
        char dft3path[600];
        snprintf(dft3path, sizeof(dft3path), "%s/dft_phi3_%s", outdir, cfg->dft_file);
        do_dft(dft3path, phi3_hist, t_hist, n_dft, dft_start, mass, 500);
    }

    /* Estimate peak frequency from phi1_hist (simple peak search) */
    {
        double T = t_hist[n_dft-1] - t_hist[dft_start];
        double best_pw = 0, best_om = 0;
        for (int k = 1; k < 500; k++) {
            double omega = 3.0 * mass * k / 500.0;
            double re = 0, im = 0;
            for (int j = dft_start; j < n_dft; j++) {
                double dtj = (j > dft_start) ? (t_hist[j]-t_hist[j-1]) :
                             (t_hist[dft_start+1]-t_hist[dft_start]);
                re += phi1_hist[j] * cos(omega * t_hist[j]) * dtj;
                im += phi1_hist[j] * sin(omega * t_hist[j]) * dtj;
            }
            double pw = (re*re + im*im) / (T*T);
            if (pw > best_pw) { best_pw = pw; best_om = omega; }
        }
        res->peak_omega = best_om;
    }

    printf("  dE/dt (t>10000) = %.6e\n", res->dEdt_late);
    printf("  Peak omega = %.4f (mass gap = %.4f)\n", res->peak_omega, mass);

    free(phi1_hist); free(phi3_hist); free(t_hist);
    state_free(&s);
}

/* ===================== main ===================== */
int main(int argc, char **argv)
{
    parse_args(argc, argv);

    double dx = 2.0 * xmax / (Nx - 1);
    double m2 = mass * mass;
    double kmax = M_PI / dx;
    double dt = 0.8 * 2.0 / sqrt(kmax*kmax + m2);
    int Nt = (int)(tfinal / dt) + 1;

    printf("rad180: Radiation comparison 0-deg vs 180-deg\n");
    printf("  mu=%.3f kappa=%.4f mass=%.4f A=%.3f sigma=%.3f\n",
           mu, kappa, mass, A_init, sigma);
    printf("  Nx=%d xmax=%.1f dx=%.5f dt=%.6f Nt=%d tfinal=%.0f\n",
           Nx, xmax, dx, dt, Nt, tfinal);

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

    /* =========== Configuration 1: 0-degree =========== */
    Config cfg0 = {
        .Ainit = {A_init, A_init, A_init},
        .label = "0-degree (+++)",
        .ts_file = "rad_0deg_ts.tsv",
        .dft_file = "rad_spectrum_0deg.tsv"
    };
    RunResult res0;
    run_config(&cfg0, dx, dt, Nt, damp, m2, &res0);

    /* =========== Configuration 2: 180-degree symmetric =========== */
    Config cfg180 = {
        .Ainit = {A_init, A_init, -A_init},
        .label = "180-degree symmetric (++-)",
        .ts_file = "rad_180deg_ts.tsv",
        .dft_file = "rad_spectrum_180deg.tsv"
    };
    RunResult res180;
    run_config(&cfg180, dx, dt, Nt, damp, m2, &res180);

    /* =========== Configuration 3: Asymmetric 180 — mild =========== */
    Config cfg_asym1 = {
        .Ainit = {1.2*A_init, 1.2*A_init, -0.8*A_init},
        .label = "180-degree asymmetric (1.2A, 1.2A, -0.8A)",
        .ts_file = "rad_asym_mild_ts.tsv",
        .dft_file = "rad_spectrum_asym_mild.tsv"
    };
    RunResult res_asym1;
    run_config(&cfg_asym1, dx, dt, Nt, damp, m2, &res_asym1);

    /* =========== Configuration 4: Asymmetric 180 — strong =========== */
    Config cfg_asym2 = {
        .Ainit = {1.5*A_init, 1.5*A_init, -0.5*A_init},
        .label = "180-degree asymmetric (1.5A, 1.5A, -0.5A)",
        .ts_file = "rad_asym_strong_ts.tsv",
        .dft_file = "rad_spectrum_asym_strong.tsv"
    };
    RunResult res_asym2;
    run_config(&cfg_asym2, dx, dt, Nt, damp, m2, &res_asym2);

    /* ===================== Summary ===================== */
    printf("\n");
    printf("=========================================\n");
    printf("         RADIATION COMPARISON SUMMARY\n");
    printf("=========================================\n");
    printf("\n");
    printf("%-35s %12s %12s %12s %12s\n",
           "Configuration", "E(t=10k)", "E(t=20k)", "dE/dt(late)", "omega");
    printf("%-35s %12s %12s %12s %12s\n",
           "-----------------------------------", "--------", "--------", "-----------", "-----");

    RunResult *results[] = {&res0, &res180, &res_asym1, &res_asym2};
    const char *labels[] = {
        "0-deg (+++)",
        "180-deg (++-)",
        "Asym 180 (1.2,1.2,-0.8)",
        "Asym 180 (1.5,1.5,-0.5)"
    };
    for (int i = 0; i < 4; i++) {
        printf("%-35s %12.4f %12.4f %12.6e %12.4f\n",
               labels[i], results[i]->E_at_t1, results[i]->E_at_t2,
               results[i]->dEdt_late, results[i]->peak_omega);
    }

    printf("\nFinal peak amplitudes:\n");
    printf("%-35s %10s %10s %10s\n", "Configuration", "|phi1|", "|phi2|", "|phi3|");
    for (int i = 0; i < 4; i++) {
        printf("%-35s %10.5f %10.5f %10.5f\n",
               labels[i], results[i]->peak_amp[0], results[i]->peak_amp[1],
               results[i]->peak_amp[2]);
    }

    /* Key comparison */
    printf("\n--- KEY COMPARISON ---\n");
    double ratio = (res180.dEdt_late != 0) ?
        res0.dEdt_late / res180.dEdt_late : 0.0;
    printf("dE/dt ratio (0deg / 180deg) = %.6f\n", ratio);
    printf("Radiation rates %s\n",
           (fabs(ratio - 1.0) < 0.01) ? "IDENTICAL (within 1%%)" :
           (fabs(ratio - 1.0) < 0.05) ? "SIMILAR (within 5%%)" :
           "DIFFERENT (>5%% apart)");

    /* Asymmetry persistence */
    printf("\n--- ASYMMETRY PERSISTENCE ---\n");
    for (int i = 2; i < 4; i++) {
        double asym = (results[i]->peak_amp[2] > 1e-10) ?
            results[i]->peak_amp[0] / results[i]->peak_amp[2] : 999.0;
        printf("%-35s  |phi1|/|phi3| = %.4f  (init: %.4f)\n",
               labels[i], asym,
               fabs(results[i] == &res_asym1 ? 1.2/0.8 : 1.5/0.5));
    }

    /* Write summary to file */
    char sumpath[600];
    snprintf(sumpath, sizeof(sumpath), "%s/summary.txt", outdir);
    FILE *fsum = fopen(sumpath, "w");
    fprintf(fsum, "# rad180 summary\n");
    fprintf(fsum, "# mu=%.3f kappa=%.4f mass=%.4f A=%.3f sigma=%.3f\n",
            mu, kappa, mass, A_init, sigma);
    fprintf(fsum, "# Nx=%d xmax=%.1f tfinal=%.0f\n\n", Nx, xmax, tfinal);
    fprintf(fsum, "%-35s %12s %12s %12s %12s %10s %10s %10s\n",
            "config", "E_t10k", "E_t20k", "dEdt_late", "omega",
            "pk1", "pk2", "pk3");
    for (int i = 0; i < 4; i++) {
        fprintf(fsum, "%-35s %12.6f %12.6f %12.6e %12.6f %10.5f %10.5f %10.5f\n",
                labels[i], results[i]->E_at_t1, results[i]->E_at_t2,
                results[i]->dEdt_late, results[i]->peak_omega,
                results[i]->peak_amp[0], results[i]->peak_amp[1],
                results[i]->peak_amp[2]);
    }
    fprintf(fsum, "\ndEdt_ratio(0/180) = %.6f\n", ratio);
    fclose(fsum);

    printf("\nOutput files in %s/\n", outdir);
    free(damp);
    return 0;
}
