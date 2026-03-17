/*
 * phase120_1d.c — 120-degree phase-separated oscillon vs symmetric (0-degree) control
 *
 * Three massive scalars with triple-product coupling.
 * Two modes:
 *   A) Symmetric (0-phase): phi_a = f(x)*cos(wt), all in phase -> P = f^3 cos^3(wt)
 *   B) 120-phase: phi_a = f(x)*cos(wt + 2*pi*a/3) -> P = (f^3/4)*cos(3wt)
 *
 * The 120-degree state eliminates fundamental and 2nd harmonic from P,
 * potentially extending lifetime. If 3*omega < m, NO radiation at all.
 *
 * Protocol:
 *   1. Equilibrate symmetric oscillon for t_equil
 *   2. Fork: (a) continue symmetric as control, (b) rotate phases to 120-degree
 *   3. Evolve both for t_run, compare energy loss and phase stability
 *
 * Also: strong-binding scan (mu=-60, -100) to check if omega < m/3 achievable.
 *
 * Compile: gcc -O3 -Wall -o phase120_1d src/phase120_1d.c -lm
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
static double t_equil = 5000.0;
static double t_run   = 10000.0;
static char   outdir[512] = "v24/phase120/data";

static void parse_args(int argc, char **argv)
{
    for (int i = 1; i < argc - 1; i += 2) {
        if      (!strcmp(argv[i], "-mu"))      mu      = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-kappa"))   kappa   = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-mass"))    mass    = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-A"))       A_init  = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-sigma"))   sigma   = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-Nx"))      Nx      = atoi(argv[i+1]);
        else if (!strcmp(argv[i], "-xmax"))    xmax    = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-t_equil")) t_equil = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-t_run"))   t_run   = atof(argv[i+1]);
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

/* State: fields, velocities, accelerations for 3 components */
typedef struct {
    double *phi[3], *vel[3], *acc[3];
} State;

static State alloc_state(int n) {
    State s;
    for (int a = 0; a < 3; a++) {
        s.phi[a] = calloc(n, sizeof(double));
        s.vel[a] = calloc(n, sizeof(double));
        s.acc[a] = calloc(n, sizeof(double));
    }
    return s;
}

static void free_state(State *s) {
    for (int a = 0; a < 3; a++) {
        free(s->phi[a]); free(s->vel[a]); free(s->acc[a]);
    }
}

static void copy_state(State *dst, const State *src, int n) {
    for (int a = 0; a < 3; a++) {
        memcpy(dst->phi[a], src->phi[a], n * sizeof(double));
        memcpy(dst->vel[a], src->vel[a], n * sizeof(double));
        memcpy(dst->acc[a], src->acc[a], n * sizeof(double));
    }
}

static double dx_g, dx2_g, m2_g;

static void compute_acc(State *s, int n) {
    for (int a = 0; a < 3; a++) {
        s->acc[a][0] = s->acc[a][1] = s->acc[a][n-2] = s->acc[a][n-1] = 0;
        for (int i = 1; i < n - 1; i++) {
            double lapl = (s->phi[a][i+1] - 2.0*s->phi[a][i] + s->phi[a][i-1]) / dx2_g;
            double fp = force_pot(s->phi[0][i], s->phi[1][i], s->phi[2][i], a);
            s->acc[a][i] = lapl - m2_g * s->phi[a][i] + fp;
        }
    }
}

static void step_verlet(State *s, double dt, double *damp, int n) {
    for (int a = 0; a < 3; a++)
        for (int i = 1; i < n - 1; i++)
            s->vel[a][i] += 0.5 * dt * s->acc[a][i];
    for (int a = 0; a < 3; a++)
        for (int i = 1; i < n - 1; i++)
            s->phi[a][i] += dt * s->vel[a][i];
    compute_acc(s, n);
    for (int a = 0; a < 3; a++)
        for (int i = 1; i < n - 1; i++)
            s->vel[a][i] += 0.5 * dt * s->acc[a][i];
    /* absorbing boundary */
    for (int a = 0; a < 3; a++)
        for (int i = 0; i < n; i++) {
            s->vel[a][i] *= damp[i];
            s->phi[a][i] *= damp[i];
        }
}

typedef struct {
    double Ek, Eg, Em, Ep, Et;
    double peak[3];
    double fc;
    double phi0[3];
} Diag;

static Diag diagnose(State *s, int n, double core_r) {
    Diag d = {0};
    int ic = n / 2;
    double Ecore = 0, Eall = 0;

    for (int a = 0; a < 3; a++)
        d.phi0[a] = s->phi[a][ic];

    for (int i = 1; i < n - 1; i++) {
        double x = -xmax + i * dx_g;
        for (int a = 0; a < 3; a++) {
            d.Ek += 0.5 * s->vel[a][i] * s->vel[a][i] * dx_g;
            double dp = (s->phi[a][i+1] - s->phi[a][i-1]) / (2.0*dx_g);
            d.Eg += 0.5 * dp * dp * dx_g;
            d.Em += 0.5 * m2_g * s->phi[a][i] * s->phi[a][i] * dx_g;
            if (fabs(s->phi[a][i]) > d.peak[a]) d.peak[a] = fabs(s->phi[a][i]);
        }
        double P = s->phi[0][i] * s->phi[1][i] * s->phi[2][i];
        double P2 = P * P;
        double V = 0.5 * mu * P2 / (1.0 + kappa * P2);
        d.Ep += V * dx_g;

        double e = V;
        for (int a = 0; a < 3; a++) {
            e += 0.5*s->vel[a][i]*s->vel[a][i];
            double dp = (s->phi[a][i+1] - s->phi[a][i-1]) / (2.0*dx_g);
            e += 0.5*dp*dp + 0.5*m2_g*s->phi[a][i]*s->phi[a][i];
        }
        Eall += e * dx_g;
        if (fabs(x) < core_r) Ecore += e * dx_g;
    }
    d.Et = d.Ek + d.Eg + d.Em + d.Ep;
    d.fc = (Eall > 1e-20) ? Ecore / Eall : 0.0;
    return d;
}

/* Measure phase difference between two oscillating signals using their
 * instantaneous values and derivatives (Hilbert-like).
 * phase_a = atan2(-vel_a, omega*phi_a) approximately */
static double measure_phase_diff(double phi_a, double vel_a, double phi_b, double vel_b, double omega)
{
    if (omega < 0.01) return 0.0;
    double pa = atan2(-vel_a/omega, phi_a);
    double pb = atan2(-vel_b/omega, phi_b);
    double diff = pa - pb;
    /* wrap to [-pi, pi] */
    while (diff > M_PI) diff -= 2.0*M_PI;
    while (diff < -M_PI) diff += 2.0*M_PI;
    return diff;
}

/* DFT of a time series, return peak frequency */
static double do_dft(double *sig, double *times, int n, int start,
                     double omega_max, int nfreq, FILE *fout)
{
    double T = times[n-1] - times[start];
    if (T < 1.0 || n - start < 50) return 0.0;

    double peak_pow = 0, peak_om = 0;
    for (int k = 0; k < nfreq; k++) {
        double omega = omega_max * k / nfreq;
        double re = 0, im = 0;
        for (int j = start; j < n; j++) {
            double dtj = (j > start) ?
                (times[j]-times[j-1]) : (times[start+1]-times[start]);
            re += sig[j] * cos(omega * times[j]) * dtj;
            im += sig[j] * sin(omega * times[j]) * dtj;
        }
        double pw = (re*re + im*im) / (T*T);
        if (fout) fprintf(fout, "%.6f\t%.6e\n", omega, pw);
        if (pw > peak_pow) { peak_pow = pw; peak_om = omega; }
    }
    return peak_om;
}

/* DFT of P(t) — triple product at center */
static double do_dft_P(double *sig, double *times, int n, int start,
                       double omega_max, int nfreq, FILE *fout)
{
    double T = times[n-1] - times[start];
    if (T < 1.0 || n - start < 50) return 0.0;

    double peak_pow = 0, peak_om = 0;
    for (int k = 0; k < nfreq; k++) {
        double omega = omega_max * k / nfreq;
        double re = 0, im = 0;
        for (int j = start; j < n; j++) {
            double dtj = (j > start) ?
                (times[j]-times[j-1]) : (times[start+1]-times[start]);
            re += sig[j] * cos(omega * times[j]) * dtj;
            im += sig[j] * sin(omega * times[j]) * dtj;
        }
        double pw = (re*re + im*im) / (T*T);
        if (fout) fprintf(fout, "%.6f\t%.6e\n", omega, pw);
        if (pw > peak_pow && omega > 0.01) { peak_pow = pw; peak_om = omega; }
    }
    return peak_om;
}

int main(int argc, char **argv)
{
    parse_args(argc, argv);

    dx_g = 2.0 * xmax / (Nx - 1);
    dx2_g = dx_g * dx_g;
    m2_g = mass * mass;

    /* CFL */
    double kmax = M_PI / dx_g;
    double dt = 0.8 * 2.0 / sqrt(kmax*kmax + m2_g);
    int Nt_equil = (int)(t_equil / dt) + 1;
    int Nt_run   = (int)(t_run / dt) + 1;

    printf("phase120_1d: mu=%.1f kappa=%.1f mass=%.3f\n", mu, kappa, mass);
    printf("  Nx=%d xmax=%.1f dx=%.5f dt=%.6f\n", Nx, xmax, dx_g, dt);
    printf("  t_equil=%.0f (Nt=%d), t_run=%.0f (Nt=%d)\n",
           t_equil, Nt_equil, t_run, Nt_run);

    /* Absorbing boundary: outer 25% */
    double *damp = malloc(Nx * sizeof(double));
    double x_abs = xmax * 0.75;
    for (int i = 0; i < Nx; i++) {
        double x = -xmax + i * dx_g;
        double ax = fabs(x);
        if (ax > x_abs) {
            double f = (ax - x_abs) / (xmax - x_abs);
            damp[i] = 1.0 - 0.98 * f * f;
        } else {
            damp[i] = 1.0;
        }
    }

    /* Initialize symmetric state */
    State sym = alloc_state(Nx);
    for (int a = 0; a < 3; a++)
        for (int i = 0; i < Nx; i++) {
            double x = -xmax + i * dx_g;
            sym.phi[a][i] = A_init * exp(-x * x / (2.0 * sigma * sigma));
        }
    compute_acc(&sym, Nx);

    double core_r = 3.0 * sigma;

    /* ===== Phase 1: Equilibrate symmetric oscillon ===== */
    printf("\n=== Equilibrating symmetric oscillon for t=%.0f ===\n", t_equil);
    int print_eq = Nt_equil / 20;
    if (print_eq < 1) print_eq = 1;

    /* Track phi1(0,t) during equilibration to measure omega */
    int max_hist_eq = 50000;
    double *phi_eq_hist = malloc(max_hist_eq * sizeof(double));
    double *t_eq_hist   = malloc(max_hist_eq * sizeof(double));
    int n_eq_hist = 0;
    int eq_rec = Nt_equil / max_hist_eq;
    if (eq_rec < 1) eq_rec = 1;

    for (int n = 0; n < Nt_equil; n++) {
        double t = n * dt;

        if (n % eq_rec == 0 && n_eq_hist < max_hist_eq) {
            phi_eq_hist[n_eq_hist] = sym.phi[0][Nx/2];
            t_eq_hist[n_eq_hist] = t;
            n_eq_hist++;
        }

        if (n % print_eq == 0) {
            Diag d = diagnose(&sym, Nx, core_r);
            printf("  equil t=%7.1f  E=%.4f  Ep=%.4f  fc=%.3f  pk=(%.3f,%.3f,%.3f)\n",
                   t, d.Et, d.Ep, d.fc, d.peak[0], d.peak[1], d.peak[2]);
        }
        step_verlet(&sym, dt, damp, Nx);
    }

    /* Measure omega from equilibrated oscillon */
    int dft_start_eq = n_eq_hist * 3 / 4;  /* last quarter */
    double omega_eq = do_dft(phi_eq_hist, t_eq_hist, n_eq_hist, dft_start_eq,
                             3.0*mass, 500, NULL);
    printf("\nEquilibrated omega = %.4f (mass gap = %.4f, m/3 = %.4f)\n",
           omega_eq, mass, mass/3.0);
    printf("  3*omega = %.4f %s m\n", 3.0*omega_eq,
           3.0*omega_eq < mass ? "<" : ">=");

    free(phi_eq_hist); free(t_eq_hist);

    /* ===== Phase 2: Fork into control (0-phase) and 120-phase states ===== */
    printf("\n=== Forking: control (0-phase) vs 120-phase ===\n");

    /* Control = copy of equilibrated symmetric state */
    State ctrl = alloc_state(Nx);
    copy_state(&ctrl, &sym, Nx);

    /* 120-phase state: rotate phases
     * Current state: phi_a = f(x)*cos(omega*t_now + theta)
     *                vel_a = -omega*f(x)*sin(omega*t_now + theta)
     * For symmetric: all same phase theta.
     *
     * We want:
     *   phi_1 = f*cos(theta), vel_1 = -omega*f*sin(theta)  [unchanged]
     *   phi_2 = f*cos(theta + 2pi/3), vel_2 = -omega*f*sin(theta + 2pi/3)
     *   phi_3 = f*cos(theta + 4pi/3), vel_3 = -omega*f*sin(theta + 4pi/3)
     *
     * From current phi_1 = f*cos(theta), vel_1 = -omega*f*sin(theta):
     *   f = sqrt(phi_1^2 + vel_1^2/omega^2)
     *   theta = atan2(-vel_1/omega, phi_1)
     *
     * Then rotate field 2 by +2pi/3 and field 3 by +4pi/3.
     */
    State ph120 = alloc_state(Nx);
    copy_state(&ph120, &sym, Nx);

    double omega_use = (omega_eq > 0.01) ? omega_eq : 0.8;
    printf("Using omega=%.4f for phase rotation\n", omega_use);

    for (int i = 0; i < Nx; i++) {
        double p1 = sym.phi[0][i];
        double v1 = sym.vel[0][i];
        double f_amp = sqrt(p1*p1 + v1*v1/(omega_use*omega_use));
        double theta = atan2(-v1/omega_use, p1);

        /* Field 1: unchanged */
        ph120.phi[0][i] = p1;
        ph120.vel[0][i] = v1;

        /* Field 2: phase + 2pi/3 */
        double th2 = theta + 2.0*M_PI/3.0;
        ph120.phi[1][i] = f_amp * cos(th2);
        ph120.vel[1][i] = -omega_use * f_amp * sin(th2);

        /* Field 3: phase + 4pi/3 */
        double th3 = theta + 4.0*M_PI/3.0;
        ph120.phi[2][i] = f_amp * cos(th3);
        ph120.vel[2][i] = -omega_use * f_amp * sin(th3);
    }
    compute_acc(&ph120, Nx);
    compute_acc(&ctrl, Nx);

    /* Verify initial energies */
    Diag d_ctrl0 = diagnose(&ctrl, Nx, core_r);
    Diag d_ph120_0 = diagnose(&ph120, Nx, core_r);
    printf("Control  E(0) = %.4f, Ep = %.4f, fc = %.3f\n",
           d_ctrl0.Et, d_ctrl0.Ep, d_ctrl0.fc);
    printf("120-deg  E(0) = %.4f, Ep = %.4f, fc = %.3f\n",
           d_ph120_0.Et, d_ph120_0.Ep, d_ph120_0.fc);

    /* ===== Phase 3: Evolve both and record ===== */
    printf("\n=== Evolving for t=%.0f ===\n", t_run);

    char path_ctrl[600], path_ph120[600];
    snprintf(path_ctrl, sizeof(path_ctrl), "%s/phase120_control_ts.tsv", outdir);
    snprintf(path_ph120, sizeof(path_ph120), "%s/phase120_ts.tsv", outdir);
    FILE *f_ctrl = fopen(path_ctrl, "w");
    FILE *f_ph120 = fopen(path_ph120, "w");
    if (!f_ctrl || !f_ph120) {
        fprintf(stderr, "Cannot open output files in %s\n", outdir);
        return 1;
    }

    fprintf(f_ctrl,  "time\tphi1_0\tphi2_0\tphi3_0\tpeak1\tpeak2\tpeak3\t"
                     "E_kin\tE_grad\tE_mass\tE_pot\tE_total\tf_core\t"
                     "dtheta12\tdtheta23\tP_center\n");
    fprintf(f_ph120, "time\tphi1_0\tphi2_0\tphi3_0\tpeak1\tpeak2\tpeak3\t"
                     "E_kin\tE_grad\tE_mass\tE_pot\tE_total\tf_core\t"
                     "dtheta12\tdtheta23\tP_center\n");

    int rec_every = Nt_run / 20000;
    if (rec_every < 1) rec_every = 1;
    int print_every = Nt_run / 40;
    if (print_every < 1) print_every = 1;

    /* DFT history */
    int max_dft = 50000;
    double *phi1_ctrl_hist  = malloc(max_dft * sizeof(double));
    double *phi1_ph120_hist = malloc(max_dft * sizeof(double));
    double *P_ctrl_hist     = malloc(max_dft * sizeof(double));
    double *P_ph120_hist    = malloc(max_dft * sizeof(double));
    double *t_hist          = malloc(max_dft * sizeof(double));
    int n_dft = 0;
    int dft_every = Nt_run / max_dft;
    if (dft_every < 1) dft_every = 1;

    int ic = Nx / 2;

    for (int n = 0; n <= Nt_run; n++) {
        double t = n * dt;

        /* DFT sampling */
        if (n % dft_every == 0 && n_dft < max_dft) {
            phi1_ctrl_hist[n_dft]  = ctrl.phi[0][ic];
            phi1_ph120_hist[n_dft] = ph120.phi[0][ic];
            P_ctrl_hist[n_dft]     = ctrl.phi[0][ic] * ctrl.phi[1][ic] * ctrl.phi[2][ic];
            P_ph120_hist[n_dft]    = ph120.phi[0][ic] * ph120.phi[1][ic] * ph120.phi[2][ic];
            t_hist[n_dft] = t;
            n_dft++;
        }

        int do_rec   = (n % rec_every == 0);
        int do_print = (n % print_every == 0);

        if (do_rec || do_print) {
            Diag dc = diagnose(&ctrl, Nx, core_r);
            Diag dp = diagnose(&ph120, Nx, core_r);

            /* Phase differences */
            double dth12_c = measure_phase_diff(ctrl.phi[0][ic], ctrl.vel[0][ic],
                                                ctrl.phi[1][ic], ctrl.vel[1][ic], omega_use);
            double dth23_c = measure_phase_diff(ctrl.phi[1][ic], ctrl.vel[1][ic],
                                                ctrl.phi[2][ic], ctrl.vel[2][ic], omega_use);
            double dth12_p = measure_phase_diff(ph120.phi[0][ic], ph120.vel[0][ic],
                                                ph120.phi[1][ic], ph120.vel[1][ic], omega_use);
            double dth23_p = measure_phase_diff(ph120.phi[1][ic], ph120.vel[1][ic],
                                                ph120.phi[2][ic], ph120.vel[2][ic], omega_use);

            double P_c = ctrl.phi[0][ic] * ctrl.phi[1][ic] * ctrl.phi[2][ic];
            double P_p = ph120.phi[0][ic] * ph120.phi[1][ic] * ph120.phi[2][ic];

            if (do_rec) {
                fprintf(f_ctrl, "%.6f\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t"
                        "%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6f\t%.6f\t%.6f\t%.6e\n",
                        t, dc.phi0[0], dc.phi0[1], dc.phi0[2],
                        dc.peak[0], dc.peak[1], dc.peak[2],
                        dc.Ek, dc.Eg, dc.Em, dc.Ep, dc.Et, dc.fc,
                        dth12_c, dth23_c, P_c);
                fprintf(f_ph120, "%.6f\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t"
                        "%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6f\t%.6f\t%.6f\t%.6e\n",
                        t, dp.phi0[0], dp.phi0[1], dp.phi0[2],
                        dp.peak[0], dp.peak[1], dp.peak[2],
                        dp.Ek, dp.Eg, dp.Em, dp.Ep, dp.Et, dp.fc,
                        dth12_p, dth23_p, P_p);
            }

            if (do_print) {
                printf("  t=%7.1f  CTRL: E=%.4f Ep=%.4f fc=%.3f pk=%.3f  "
                       "120: E=%.4f Ep=%.4f fc=%.3f pk=%.3f  "
                       "dth12=%.2f dth23=%.2f\n",
                       t, dc.Et, dc.Ep, dc.fc, dc.peak[0],
                       dp.Et, dp.Ep, dp.fc, dp.peak[0],
                       dth12_p*180.0/M_PI, dth23_p*180.0/M_PI);
            }
        }

        if (n == Nt_run) break;

        step_verlet(&ctrl, dt, damp, Nx);
        step_verlet(&ph120, dt, damp, Nx);
    }

    fclose(f_ctrl);
    fclose(f_ph120);

    /* ===== DFT analysis ===== */
    printf("\n=== DFT Analysis ===\n");
    int dft_start = n_dft / 2;

    /* phi_1 spectrum for control */
    char dft_path[600];
    snprintf(dft_path, sizeof(dft_path), "%s/phase120_ctrl_phi_spectrum.tsv", outdir);
    FILE *fdft = fopen(dft_path, "w");
    fprintf(fdft, "omega\tpower\n");
    double omega_ctrl = do_dft(phi1_ctrl_hist, t_hist, n_dft, dft_start, 3.0*mass, 500, fdft);
    fclose(fdft);
    printf("Control phi_1 spectrum: peak omega = %.4f\n", omega_ctrl);

    /* phi_1 spectrum for 120-phase */
    snprintf(dft_path, sizeof(dft_path), "%s/phase120_ph120_phi_spectrum.tsv", outdir);
    fdft = fopen(dft_path, "w");
    fprintf(fdft, "omega\tpower\n");
    double omega_ph120 = do_dft(phi1_ph120_hist, t_hist, n_dft, dft_start, 3.0*mass, 500, fdft);
    fclose(fdft);
    printf("120-deg phi_1 spectrum: peak omega = %.4f\n", omega_ph120);

    /* P(t) spectrum for control */
    snprintf(dft_path, sizeof(dft_path), "%s/phase120_ctrl_P_spectrum.tsv", outdir);
    fdft = fopen(dft_path, "w");
    fprintf(fdft, "omega\tpower\n");
    double omP_ctrl = do_dft_P(P_ctrl_hist, t_hist, n_dft, dft_start, 4.0*mass, 600, fdft);
    fclose(fdft);
    printf("Control P(t) spectrum: peak omega = %.4f\n", omP_ctrl);

    /* P(t) spectrum for 120-phase */
    snprintf(dft_path, sizeof(dft_path), "%s/phase120_ph120_P_spectrum.tsv", outdir);
    fdft = fopen(dft_path, "w");
    fprintf(fdft, "omega\tpower\n");
    double omP_ph120 = do_dft_P(P_ph120_hist, t_hist, n_dft, dft_start, 4.0*mass, 600, fdft);
    fclose(fdft);
    printf("120-deg P(t) spectrum: peak omega = %.4f\n", omP_ph120);

    /* Final energies */
    Diag dc_final = diagnose(&ctrl, Nx, core_r);
    Diag dp_final = diagnose(&ph120, Nx, core_r);

    printf("\n=== Summary ===\n");
    printf("  omega_equil = %.4f (m/3 = %.4f, 3*omega = %.4f)\n",
           omega_eq, mass/3.0, 3.0*omega_eq);
    printf("  Control:  E(0)=%.4f -> E(end)=%.4f  (retained %.1f%%), fc=%.3f\n",
           d_ctrl0.Et, dc_final.Et,
           100.0*dc_final.Et/d_ctrl0.Et, dc_final.fc);
    printf("  120-deg:  E(0)=%.4f -> E(end)=%.4f  (retained %.1f%%), fc=%.3f\n",
           d_ph120_0.Et, dp_final.Et,
           100.0*dp_final.Et/d_ph120_0.Et, dp_final.fc);
    printf("  3*omega < m? %s\n", 3.0*omega_eq < mass ? "YES (sub-triple-harmonic)" : "NO");

    /* Summary file */
    char sumpath[600];
    snprintf(sumpath, sizeof(sumpath), "%s/phase120_summary.tsv", outdir);
    FILE *fsum = fopen(sumpath, "w");
    fprintf(fsum, "parameter\tvalue\n");
    fprintf(fsum, "mu\t%.1f\n", mu);
    fprintf(fsum, "kappa\t%.1f\n", kappa);
    fprintf(fsum, "mass\t%.3f\n", mass);
    fprintf(fsum, "omega_equil\t%.6f\n", omega_eq);
    fprintf(fsum, "3*omega\t%.6f\n", 3.0*omega_eq);
    fprintf(fsum, "m/3\t%.6f\n", mass/3.0);
    fprintf(fsum, "sub_triple\t%d\n", 3.0*omega_eq < mass ? 1 : 0);
    fprintf(fsum, "ctrl_E_init\t%.6f\n", d_ctrl0.Et);
    fprintf(fsum, "ctrl_E_final\t%.6f\n", dc_final.Et);
    fprintf(fsum, "ctrl_E_retained\t%.6f\n", dc_final.Et/d_ctrl0.Et);
    fprintf(fsum, "ctrl_fc_final\t%.6f\n", dc_final.fc);
    fprintf(fsum, "ctrl_omega\t%.6f\n", omega_ctrl);
    fprintf(fsum, "ph120_E_init\t%.6f\n", d_ph120_0.Et);
    fprintf(fsum, "ph120_E_final\t%.6f\n", dp_final.Et);
    fprintf(fsum, "ph120_E_retained\t%.6f\n", dp_final.Et/d_ph120_0.Et);
    fprintf(fsum, "ph120_fc_final\t%.6f\n", dp_final.fc);
    fprintf(fsum, "ph120_omega\t%.6f\n", omega_ph120);
    fprintf(fsum, "ctrl_P_peak_omega\t%.6f\n", omP_ctrl);
    fprintf(fsum, "ph120_P_peak_omega\t%.6f\n", omP_ph120);
    fclose(fsum);

    printf("\nOutput: %s, %s, %s\n", path_ctrl, path_ph120, sumpath);

    free_state(&sym); free_state(&ctrl); free_state(&ph120);
    free(damp);
    free(phi1_ctrl_hist); free(phi1_ph120_hist);
    free(P_ctrl_hist); free(P_ph120_hist); free(t_hist);
    return 0;
}
