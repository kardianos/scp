/*
 * char180.c — Characterize the 180-degree anti-phase oscillon
 *
 * Runs side-by-side: 0-degree (all +) vs 180-degree (phi3 flipped)
 * Then stability tests on the 180-degree state.
 *
 * Lagrangian:
 *   L = sum_a [ (1/2)(dt phi_a)^2 - (1/2)(dx phi_a)^2 - (m^2/2)phi_a^2 ]
 *     - (mu/2) P^2 / (1 + kappa P^2)
 *   P = phi_1 phi_2 phi_3
 *
 * Compile: gcc -O3 -Wall -o char180 src/char180.c -lm
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Parameters from proposal */
static double mu     = -20.0;
static double kappa  = 20.0;
static double mass   = 1.0;
static double A_init = 0.8;
static double sigma  = 3.0;
static int    Nx     = 4000;
static double xmax   = 100.0;
static double tfinal = 20000.0;
static char   outdir[512] = "data";

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
    double P  = p1 * p2 * p3;
    double P2 = P * P;
    return 0.5 * mu * P2 / (1.0 + kappa * P2);
}

/* ===================== Single Oscillon Simulation ===================== */

typedef struct {
    double *phi[3];
    double *vel[3];
    double *acc[3];
    double *damp;
    int Nx;
    double dx, dx2, dt, m2;
    double xmax;
} Sim;

static Sim *sim_create(int nx, double xm)
{
    Sim *s = calloc(1, sizeof(Sim));
    s->Nx = nx;
    s->xmax = xm;
    s->dx = 2.0 * xm / (nx - 1);
    s->dx2 = s->dx * s->dx;
    s->m2 = mass * mass;

    double kmax = M_PI / s->dx;
    s->dt = 0.8 * 2.0 / sqrt(kmax * kmax + s->m2);

    for (int a = 0; a < 3; a++) {
        s->phi[a] = calloc(nx, sizeof(double));
        s->vel[a] = calloc(nx, sizeof(double));
        s->acc[a] = calloc(nx, sizeof(double));
    }

    /* Absorbing boundary: outer 25% */
    s->damp = malloc(nx * sizeof(double));
    double x_abs = xm * 0.75;
    for (int i = 0; i < nx; i++) {
        double x = -xm + i * s->dx;
        double ax = fabs(x);
        if (ax > x_abs) {
            double f = (ax - x_abs) / (xm - x_abs);
            s->damp[i] = 1.0 - 0.98 * f * f;
        } else {
            s->damp[i] = 1.0;
        }
    }
    return s;
}

static void sim_free(Sim *s)
{
    for (int a = 0; a < 3; a++) {
        free(s->phi[a]); free(s->vel[a]); free(s->acc[a]);
    }
    free(s->damp);
    free(s);
}

static void sim_compute_acc(Sim *s)
{
    int nx = s->Nx;
    for (int a = 0; a < 3; a++) {
        s->acc[a][0] = s->acc[a][1] = s->acc[a][nx-2] = s->acc[a][nx-1] = 0;
        for (int i = 1; i < nx - 1; i++) {
            double lapl = (s->phi[a][i+1] - 2.0*s->phi[a][i] + s->phi[a][i-1]) / s->dx2;
            double fp = force_pot(s->phi[0][i], s->phi[1][i], s->phi[2][i], a);
            s->acc[a][i] = lapl - s->m2 * s->phi[a][i] + fp;
        }
    }
}

static void sim_step(Sim *s)
{
    int nx = s->Nx;
    double dt = s->dt;

    /* Velocity Verlet */
    for (int a = 0; a < 3; a++)
        for (int i = 1; i < nx - 1; i++)
            s->vel[a][i] += 0.5 * dt * s->acc[a][i];
    for (int a = 0; a < 3; a++)
        for (int i = 1; i < nx - 1; i++)
            s->phi[a][i] += dt * s->vel[a][i];
    sim_compute_acc(s);
    for (int a = 0; a < 3; a++)
        for (int i = 1; i < nx - 1; i++)
            s->vel[a][i] += 0.5 * dt * s->acc[a][i];

    /* absorbing boundary */
    for (int a = 0; a < 3; a++)
        for (int i = 0; i < nx; i++) {
            s->vel[a][i] *= s->damp[i];
            s->phi[a][i] *= s->damp[i];
        }
}

typedef struct {
    double Ek, Eg, Em, Ep, Et;
    double fc;
    double peak[3];
    double phi0[3];  /* center values */
    double P_center;
} Diag;

static Diag sim_diag(Sim *s)
{
    Diag d = {0};
    int ic = s->Nx / 2;
    double core_r = 3.0 * sigma;
    double Ecore = 0, Eall = 0;

    for (int a = 0; a < 3; a++)
        d.phi0[a] = s->phi[a][ic];
    d.P_center = s->phi[0][ic] * s->phi[1][ic] * s->phi[2][ic];

    for (int i = 1; i < s->Nx - 1; i++) {
        double x = -s->xmax + i * s->dx;
        for (int a = 0; a < 3; a++) {
            d.Ek += 0.5 * s->vel[a][i] * s->vel[a][i] * s->dx;
            double dp = (s->phi[a][i+1] - s->phi[a][i-1]) / (2.0 * s->dx);
            d.Eg += 0.5 * dp * dp * s->dx;
            d.Em += 0.5 * s->m2 * s->phi[a][i] * s->phi[a][i] * s->dx;
            if (fabs(s->phi[a][i]) > d.peak[a]) d.peak[a] = fabs(s->phi[a][i]);
        }
        double V = V_pot(s->phi[0][i], s->phi[1][i], s->phi[2][i]);
        d.Ep += V * s->dx;

        double e = V;
        for (int a = 0; a < 3; a++) {
            e += 0.5 * s->vel[a][i] * s->vel[a][i];
            double dp = (s->phi[a][i+1] - s->phi[a][i-1]) / (2.0 * s->dx);
            e += 0.5 * dp * dp + 0.5 * s->m2 * s->phi[a][i] * s->phi[a][i];
        }
        Eall += e * s->dx;
        if (fabs(x) < core_r) Ecore += e * s->dx;
    }
    d.Et = d.Ek + d.Eg + d.Em + d.Ep;
    d.fc = (Eall > 1e-20) ? Ecore / Eall : 0.0;
    return d;
}

/* ===================== DFT ===================== */

typedef struct {
    double *data;
    double *time;
    int n, cap;
} TimeSeries;

static TimeSeries ts_create(int cap)
{
    TimeSeries ts;
    ts.cap = cap;
    ts.n = 0;
    ts.data = malloc(cap * sizeof(double));
    ts.time = malloc(cap * sizeof(double));
    return ts;
}

static void ts_add(TimeSeries *ts, double t, double val)
{
    if (ts->n < ts->cap) {
        ts->time[ts->n] = t;
        ts->data[ts->n] = val;
        ts->n++;
    }
}

static void ts_free(TimeSeries *ts)
{
    free(ts->data); free(ts->time);
}

/* DFT of second half, write spectrum */
static double dft_peak(TimeSeries *ts, const char *path, double omega_max, int nfreq)
{
    int start = ts->n / 2;
    if (ts->n - start < 100) return -1.0;

    FILE *f = fopen(path, "w");
    fprintf(f, "omega\tpower\n");

    double T = ts->time[ts->n - 1] - ts->time[start];
    double peak_pow = 0, peak_om = 0;

    for (int k = 0; k < nfreq; k++) {
        double omega = omega_max * k / nfreq;
        double re = 0, im = 0;
        for (int j = start; j < ts->n; j++) {
            double dtj = (j > start) ?
                (ts->time[j] - ts->time[j-1]) :
                (ts->time[start+1] - ts->time[start]);
            re += ts->data[j] * cos(omega * ts->time[j]) * dtj;
            im += ts->data[j] * sin(omega * ts->time[j]) * dtj;
        }
        double pw = (re * re + im * im) / (T * T);
        fprintf(f, "%.6f\t%.6e\n", omega, pw);
        if (pw > peak_pow) { peak_pow = pw; peak_om = omega; }
    }
    fclose(f);
    return peak_om;
}

/* ===================== Run One Simulation ===================== */

typedef struct {
    double omega_phi;   /* peak frequency of phi_1(0,t) */
    double omega_P;     /* peak frequency of P(0,t) */
    double E_init, E_final;
    double dEdt;        /* linear fit */
    double fc_avg;
    double peak_avg;
} RunResult;

static RunResult run_sim(Sim *s, const char *label, const char *ts_path)
{
    RunResult res = {0};
    int Nt = (int)(tfinal / s->dt) + 1;

    sim_compute_acc(s);

    FILE *fts = fopen(ts_path, "w");
    fprintf(fts, "time\tphi1_0\tphi2_0\tphi3_0\tpeak1\tpeak2\tpeak3\t"
                 "E_kin\tE_grad\tE_mass\tE_pot\tE_total\tf_core\tP_center\n");

    int max_dft = 50000;
    TimeSeries ts_phi = ts_create(max_dft);
    TimeSeries ts_P   = ts_create(max_dft);

    int rec_every = Nt / 20000;
    if (rec_every < 1) rec_every = 1;
    int print_every = Nt / 40;
    if (print_every < 1) print_every = 1;
    int dft_every = Nt / max_dft;
    if (dft_every < 1) dft_every = 1;

    /* For dE/dt linear fit: collect (t, E) pairs */
    int nE = 0, maxE = 2000;
    double *tE = malloc(maxE * sizeof(double));
    double *EE = malloc(maxE * sizeof(double));
    int E_every = Nt / maxE;
    if (E_every < 1) E_every = 1;

    double fc_sum = 0;
    int fc_count = 0;
    double peak_sum = 0;
    int peak_count = 0;

    for (int n = 0; n <= Nt; n++) {
        double t = n * s->dt;

        if (n % dft_every == 0) {
            int ic = s->Nx / 2;
            ts_add(&ts_phi, t, s->phi[0][ic]);
            ts_add(&ts_P, t, s->phi[0][ic] * s->phi[1][ic] * s->phi[2][ic]);
        }

        int do_rec = (n % rec_every == 0);
        int do_print = (n % print_every == 0);

        if (do_rec || do_print || (n % E_every == 0)) {
            Diag d = sim_diag(s);

            if (n % E_every == 0 && nE < maxE) {
                tE[nE] = t;
                EE[nE] = d.Et;
                nE++;
            }

            if (n == 0) res.E_init = d.Et;
            if (n == Nt) res.E_final = d.Et;

            /* average fc and peak in second half */
            if (t > tfinal * 0.5) {
                fc_sum += d.fc; fc_count++;
                peak_sum += d.peak[0]; peak_count++;
            }

            if (do_rec)
                fprintf(fts, "%.6f\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t"
                             "%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6f\t%.6e\n",
                        t, d.phi0[0], d.phi0[1], d.phi0[2],
                        d.peak[0], d.peak[1], d.peak[2],
                        d.Ek, d.Eg, d.Em, d.Ep, d.Et, d.fc, d.P_center);
            if (do_print)
                printf("  [%s] t=%7.1f  p0=(%+.3f,%+.3f,%+.3f)  pk=(%.3f,%.3f,%.3f)  "
                       "E=%+.4f  Ep=%+.4f  fc=%.3f\n",
                       label, t, d.phi0[0], d.phi0[1], d.phi0[2],
                       d.peak[0], d.peak[1], d.peak[2], d.Et, d.Ep, d.fc);
        }

        if (n == Nt) break;
        sim_step(s);
    }
    fclose(fts);

    /* dE/dt via linear regression on second half */
    {
        int start = nE / 2;
        double st = 0, sE = 0, stE = 0, st2 = 0;
        int nn = nE - start;
        for (int i = start; i < nE; i++) {
            st += tE[i]; sE += EE[i]; stE += tE[i]*EE[i]; st2 += tE[i]*tE[i];
        }
        double denom = nn * st2 - st * st;
        if (fabs(denom) > 1e-30)
            res.dEdt = (nn * stE - st * sE) / denom;
    }

    res.fc_avg = fc_count > 0 ? fc_sum / fc_count : 0;
    res.peak_avg = peak_count > 0 ? peak_sum / peak_count : 0;

    /* DFT */
    char specpath[600];
    snprintf(specpath, sizeof(specpath), "%s/%s_phi_spectrum.tsv", outdir, label);
    res.omega_phi = dft_peak(&ts_phi, specpath, 3.0 * mass, 1000);

    snprintf(specpath, sizeof(specpath), "%s/%s_P_spectrum.tsv", outdir, label);
    res.omega_P = dft_peak(&ts_P, specpath, 5.0 * mass, 1500);

    ts_free(&ts_phi);
    ts_free(&ts_P);
    free(tE); free(EE);

    printf("\n  [%s] SUMMARY: omega_phi=%.4f  omega_P=%.4f  E_init=%.4f  E_final=%.4f  "
           "dE/dt=%.2e  fc=%.4f  peak=%.4f\n\n",
           label, res.omega_phi, res.omega_P, res.E_init, res.E_final,
           res.dEdt, res.fc_avg, res.peak_avg);

    return res;
}

/* ===================== Main ===================== */

int main(int argc, char **argv)
{
    /* Parse optional outdir */
    for (int i = 1; i < argc - 1; i += 2) {
        if (!strcmp(argv[i], "-o")) strncpy(outdir, argv[i+1], sizeof(outdir)-1);
        else if (!strcmp(argv[i], "-mu"))     mu     = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-kappa"))  kappa  = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-mass"))   mass   = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-A"))      A_init = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-sigma"))  sigma  = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-Nx"))     Nx     = atoi(argv[i+1]);
        else if (!strcmp(argv[i], "-xmax"))   xmax   = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-tfinal")) tfinal = atof(argv[i+1]);
    }

    printf("=== V24-180A: Characterize the 180-degree Anti-Phase Oscillon ===\n");
    printf("  mu=%.1f kappa=%.1f mass=%.1f A=%.3f sigma=%.1f\n", mu, kappa, mass, A_init, sigma);
    printf("  Nx=%d xmax=%.0f tfinal=%.0f\n\n", Nx, xmax, tfinal);

    /* ============ Phase 1+2: Side-by-side 0 vs 180 ============ */
    printf("====== PHASE 1+2: Side-by-side comparison ======\n\n");

    /* 0-degree */
    {
        printf("--- 0-degree oscillon ---\n");
        Sim *s = sim_create(Nx, xmax);
        for (int a = 0; a < 3; a++)
            for (int i = 0; i < s->Nx; i++) {
                double x = -s->xmax + i * s->dx;
                s->phi[a][i] = A_init * exp(-x * x / (2.0 * sigma * sigma));
            }
        char path[600];
        snprintf(path, sizeof(path), "%s/char180_0deg_ts.tsv", outdir);
        RunResult r0 = run_sim(s, "0deg", path);
        sim_free(s);

        printf("0-degree result:\n");
        printf("  omega_phi = %.4f  (omega/m = %.4f)\n", r0.omega_phi, r0.omega_phi/mass);
        printf("  omega_P   = %.4f\n", r0.omega_P);
        printf("  E_init    = %.4f\n", r0.E_init);
        printf("  E_final   = %.4f\n", r0.E_final);
        printf("  dE/dt     = %.4e\n", r0.dEdt);
        printf("  fc_avg    = %.4f\n", r0.fc_avg);
        printf("  peak_avg  = %.4f\n\n", r0.peak_avg);
    }

    /* 180-degree */
    RunResult r180;
    {
        printf("--- 180-degree oscillon ---\n");
        Sim *s = sim_create(Nx, xmax);
        for (int a = 0; a < 3; a++) {
            double sign = (a == 2) ? -1.0 : +1.0;
            for (int i = 0; i < s->Nx; i++) {
                double x = -s->xmax + i * s->dx;
                s->phi[a][i] = sign * A_init * exp(-x * x / (2.0 * sigma * sigma));
            }
        }
        char path[600];
        snprintf(path, sizeof(path), "%s/char180_180deg_ts.tsv", outdir);
        r180 = run_sim(s, "180deg", path);
        sim_free(s);

        printf("180-degree result:\n");
        printf("  omega_phi = %.4f  (omega/m = %.4f)\n", r180.omega_phi, r180.omega_phi/mass);
        printf("  omega_P   = %.4f\n", r180.omega_P);
        printf("  E_init    = %.4f\n", r180.E_init);
        printf("  E_final   = %.4f\n", r180.E_final);
        printf("  dE/dt     = %.4e\n", r180.dEdt);
        printf("  fc_avg    = %.4f\n", r180.fc_avg);
        printf("  peak_avg  = %.4f\n\n", r180.peak_avg);
    }

    /* ============ Phase 3: Stability Tests ============ */
    printf("====== PHASE 3: Stability Tests ======\n\n");

    /* 3a: Flip perturbation — add +epsilon to phi3 (push toward 0 state) */
    {
        printf("--- Stability 3a: flip perturbation (+0.1 on phi3) ---\n");
        Sim *s = sim_create(Nx, xmax);
        double eps = 0.1;
        for (int a = 0; a < 3; a++) {
            double sign = (a == 2) ? -1.0 : +1.0;
            for (int i = 0; i < s->Nx; i++) {
                double x = -s->xmax + i * s->dx;
                double g = exp(-x * x / (2.0 * sigma * sigma));
                s->phi[a][i] = sign * A_init * g;
                if (a == 2) s->phi[a][i] += eps * A_init * g; /* push toward + */
            }
        }
        char path[600];
        snprintf(path, sizeof(path), "%s/char180_stab_flip_ts.tsv", outdir);
        RunResult rf = run_sim(s, "flip", path);
        sim_free(s);

        printf("Flip perturbation result:\n");
        printf("  omega_phi = %.4f  E_final = %.4f  dE/dt = %.4e  fc = %.4f\n\n",
               rf.omega_phi, rf.E_final, rf.dEdt, rf.fc_avg);
    }

    /* 3b: Amplitude perturbation — scale all by 1.2 */
    {
        printf("--- Stability 3b: amplitude perturbation (1.2x) ---\n");
        Sim *s = sim_create(Nx, xmax);
        double scale = 1.2;
        for (int a = 0; a < 3; a++) {
            double sign = (a == 2) ? -1.0 : +1.0;
            for (int i = 0; i < s->Nx; i++) {
                double x = -s->xmax + i * s->dx;
                s->phi[a][i] = sign * scale * A_init * exp(-x * x / (2.0 * sigma * sigma));
            }
        }
        char path[600];
        snprintf(path, sizeof(path), "%s/char180_stab_amp_ts.tsv", outdir);
        RunResult ra = run_sim(s, "amp", path);
        sim_free(s);

        printf("Amplitude perturbation result:\n");
        printf("  omega_phi = %.4f  E_final = %.4f  dE/dt = %.4e  fc = %.4f\n\n",
               ra.omega_phi, ra.E_final, ra.dEdt, ra.fc_avg);
    }

    /* 3c: Asymmetric perturbation — phi1 at 0.9A, phi2 at 1.1A */
    {
        printf("--- Stability 3c: asymmetric perturbation (0.9A, 1.1A, -A) ---\n");
        Sim *s = sim_create(Nx, xmax);
        double scales[3] = {0.9, 1.1, -1.0};
        for (int a = 0; a < 3; a++) {
            for (int i = 0; i < s->Nx; i++) {
                double x = -s->xmax + i * s->dx;
                s->phi[a][i] = scales[a] * A_init * exp(-x * x / (2.0 * sigma * sigma));
            }
        }
        char path[600];
        snprintf(path, sizeof(path), "%s/char180_stab_asym_ts.tsv", outdir);
        RunResult ras = run_sim(s, "asym", path);
        sim_free(s);

        printf("Asymmetric perturbation result:\n");
        printf("  omega_phi = %.4f  E_final = %.4f  dE/dt = %.4e  fc = %.4f\n\n",
               ras.omega_phi, ras.E_final, ras.dEdt, ras.fc_avg);
    }

    printf("=== DONE ===\n");
    return 0;
}
