/*
 * pt120.c — Pairwise + Triple Product 120-degree oscillon scanner
 *
 * V24-PT: Three massive scalars with BOTH:
 *   V_pw     = lambda * (phi1*phi2 + phi2*phi3 + phi3*phi1)  [pairwise]
 *   V_triple = (mu/2) P^2 / (1 + kappa P^2),  P = phi1*phi2*phi3  [triple]
 *
 * EOM for field a:
 *   d^2 phi_a/dt^2 = d^2 phi_a/dx^2 - m^2 phi_a
 *                    - lambda * (phi_b + phi_c)
 *                    - mu * P * dP/dphi_a / (1 + kappa P^2)^2
 *
 * Normal modes (linearized about phi=0):
 *   Symmetric:     phi_a = (1/sqrt(3))*(phi_S),  m^2_S = m^2 + 2*lambda
 *   Antisymmetric: two modes with  m^2_A = m^2 - lambda
 *
 * 120-degree phasing: phi_a(x,t) = f(x)*cos(omega*t + 2*pi*a/3)
 *   - Projects purely onto antisymmetric modes
 *   - Frequency omega < m_A = sqrt(m^2 - lambda)
 *   - Triple product P = (f^3/4)*cos(3*omega*t)
 *   - At spatial infinity: free fields with mass m
 *   - Radiation threshold: 3*omega < m => NO radiation
 *
 * Protocol for each lambda:
 *   1. Initialize directly as 120-degree Gaussian with omega ~ 0.8*m_A
 *   2. Let system settle for t_equil (absorbing boundary eats transients)
 *   3. Evolve for t_run, measuring omega, dE/dt, stability
 *   4. Record whether 3*omega < m (sub-threshold)
 *
 * Compile: gcc -O3 -Wall -o pt120 src/pt120.c -lm
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Default parameters */
static double mu      = -20.0;
static double kappa   = 20.0;
static double mass    = 1.0;
static double A_init  = 0.8;
static double sig     = 3.0;
static int    Nx      = 4000;
static double xmax    = 100.0;
static double t_equil = 3000.0;   /* settle transients */
static double t_run   = 20000.0;  /* main measurement run */
static char   outdir[512] = "v24/120triple/data";

/* Scan parameters */
static double lam_list[] = {0.80, 0.82, 0.84, 0.86, 0.87, 0.88, 0.89, 0.90, 0.92, 0.95};
static int    n_lam = 10;
static double lam_single = -1.0;  /* if >= 0, run single lambda only */

static void parse_args(int argc, char **argv)
{
    for (int i = 1; i < argc - 1; i += 2) {
        if      (!strcmp(argv[i], "-mu"))      mu      = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-kappa"))   kappa   = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-mass"))    mass    = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-A"))       A_init  = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-sigma"))   sig     = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-Nx"))      Nx      = atoi(argv[i+1]);
        else if (!strcmp(argv[i], "-xmax"))    xmax    = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-t_equil")) t_equil = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-t_run"))   t_run   = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-lambda"))  lam_single = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-o"))       strncpy(outdir, argv[i+1], sizeof(outdir)-1);
        else { fprintf(stderr, "Unknown flag: %s\n", argv[i]); exit(1); }
    }
}

/* Global grid parameters */
static double dx_g, dx2_g, m2_g;
static double lambda_g;

static inline double force_triple(double p1, double p2, double p3, int a)
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

static inline double force_pairwise(double p1, double p2, double p3, int a)
{
    switch (a) {
        case 0: return -lambda_g * (p2 + p3);
        case 1: return -lambda_g * (p1 + p3);
        case 2: return -lambda_g * (p1 + p2);
        default: return 0.0;
    }
}

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

static void compute_acc(State *s, int n) {
    for (int a = 0; a < 3; a++) {
        s->acc[a][0] = s->acc[a][1] = s->acc[a][n-2] = s->acc[a][n-1] = 0;
        for (int i = 1; i < n - 1; i++) {
            double lapl = (s->phi[a][i+1] - 2.0*s->phi[a][i] + s->phi[a][i-1]) / dx2_g;
            double ft = force_triple(s->phi[0][i], s->phi[1][i], s->phi[2][i], a);
            double fp = force_pairwise(s->phi[0][i], s->phi[1][i], s->phi[2][i], a);
            s->acc[a][i] = lapl - m2_g * s->phi[a][i] + ft + fp;
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
    for (int a = 0; a < 3; a++)
        for (int i = 0; i < n; i++) {
            s->vel[a][i] *= damp[i];
            s->phi[a][i] *= damp[i];
        }
}

typedef struct {
    double Ek, Eg, Em, Ep_triple, Ep_pw, Et;
    double peak[3];
    double fc;
    double phi0[3];
} Diag;

static Diag diagnose(State *s, int n, double core_r) {
    Diag d = {0};
    int ic = n / 2;
    double Ecore_abs = 0, Eall_abs = 0;

    for (int a = 0; a < 3; a++)
        d.phi0[a] = s->phi[a][ic];

    for (int i = 1; i < n - 1; i++) {
        double x = -xmax + i * dx_g;
        double ek_i = 0, eg_i = 0, em_i = 0;
        for (int a = 0; a < 3; a++) {
            ek_i += 0.5 * s->vel[a][i] * s->vel[a][i];
            double dp = (s->phi[a][i+1] - s->phi[a][i-1]) / (2.0*dx_g);
            eg_i += 0.5 * dp * dp;
            em_i += 0.5 * m2_g * s->phi[a][i] * s->phi[a][i];
            if (fabs(s->phi[a][i]) > d.peak[a]) d.peak[a] = fabs(s->phi[a][i]);
        }
        d.Ek += ek_i * dx_g;
        d.Eg += eg_i * dx_g;
        d.Em += em_i * dx_g;

        double p1 = s->phi[0][i], p2 = s->phi[1][i], p3 = s->phi[2][i];
        double P = p1 * p2 * p3;
        double P2 = P * P;
        double vt = 0.5 * mu * P2 / (1.0 + kappa * P2);
        double vpw = lambda_g * (p1*p2 + p2*p3 + p3*p1);
        d.Ep_triple += vt * dx_g;
        d.Ep_pw += vpw * dx_g;

        double e_abs = ek_i + eg_i + em_i + fabs(vt) + fabs(vpw);
        Eall_abs += e_abs * dx_g;
        if (fabs(x) < core_r) Ecore_abs += e_abs * dx_g;
    }
    d.Et = d.Ek + d.Eg + d.Em + d.Ep_triple + d.Ep_pw;
    d.fc = (Eall_abs > 1e-20) ? Ecore_abs / Eall_abs : 0.0;
    return d;
}

static double do_dft(double *sig_data, double *times, int n, int start,
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
            re += sig_data[j] * cos(omega * times[j]) * dtj;
            im += sig_data[j] * sin(omega * times[j]) * dtj;
        }
        double pw = (re*re + im*im) / (T*T);
        if (fout) fprintf(fout, "%.6f\t%.6e\n", omega, pw);
        if (pw > peak_pow && omega > 0.01) { peak_pow = pw; peak_om = omega; }
    }
    return peak_om;
}

static double measure_phase_diff(double phi_a, double vel_a, double phi_b, double vel_b, double omega)
{
    if (omega < 0.01) return 0.0;
    double pa = atan2(-vel_a/omega, phi_a);
    double pb = atan2(-vel_b/omega, phi_b);
    double diff = pa - pb;
    while (diff > M_PI) diff -= 2.0*M_PI;
    while (diff < -M_PI) diff += 2.0*M_PI;
    return diff;
}

typedef struct {
    double lambda;
    double m2_anti, m_anti;
    double omega;        /* measured late-time frequency */
    double three_omega;
    int    sub_thresh;   /* 3*omega < m */
    double E_init, E_settled, E_final;
    double E_retained;   /* E_final / E_settled */
    double dEdt_early, dEdt_late;
    double fc_final, peak_final;
    double phase12_final, phase23_final;
} Result;

static Result run_one(double lam, double *damp)
{
    Result res = {0};
    res.lambda = lam;
    res.m2_anti = m2_g - lam;
    res.m_anti = (res.m2_anti > 0) ? sqrt(res.m2_anti) : 0.0;

    lambda_g = lam;

    printf("\n======================================================\n");
    printf("  lambda=%.4f  m^2_A=%.4f  m_A=%.4f  m_S=%.4f\n",
           lam, res.m2_anti, res.m_anti, sqrt(m2_g + 2*lam));
    if (res.m2_anti <= 0.0) {
        printf("  TACHYONIC! Skipping.\n");
        return res;
    }

    /* Predicted threshold */
    double om_pred = 0.8 * res.m_anti;
    printf("  Predicted omega~%.3f, 3omega~%.3f %s m=%.3f\n",
           om_pred, 3*om_pred, 3*om_pred < mass ? "<" : ">=", mass);
    printf("======================================================\n");

    double kmax = M_PI / dx_g;
    double dt = 0.8 * 2.0 / sqrt(kmax*kmax + m2_g + 2.0*fabs(lam));
    int Nt_equil = (int)(t_equil / dt) + 1;
    int Nt_run   = (int)(t_run / dt) + 1;
    int ic = Nx / 2;
    double core_r = 3.0 * sig;

    printf("  dt=%.6f  Nt_equil=%d  Nt_run=%d\n", dt, Nt_equil, Nt_run);

    /* ===== Initialize directly as 120-degree state ===== */
    /* phi_a(x,0) = A * exp(-x^2/(2*sigma^2)) * cos(2*pi*a/3)
     * vel_a(x,0) = -omega_est * A * exp(-x^2/(2*sigma^2)) * sin(2*pi*a/3)
     *
     * Width sigma adapted to m_A: sigma ~ 1/sqrt(m_A^2 - omega^2)
     */
    State st = alloc_state(Nx);

    double omega_est = om_pred;
    if (omega_est < 0.05) omega_est = 0.1;
    double kappa_loc = sqrt(res.m2_anti - omega_est*omega_est);
    double sig_use = (kappa_loc > 0.05) ? 1.0/kappa_loc : sig;
    if (sig_use > 2.0*sig) sig_use = 2.0*sig;
    if (sig_use < 1.0) sig_use = 1.0;

    double phases[3] = {0.0, 2.0*M_PI/3.0, 4.0*M_PI/3.0};
    printf("  Init: A=%.3f sigma=%.3f omega_est=%.4f\n", A_init, sig_use, omega_est);

    for (int a = 0; a < 3; a++)
        for (int i = 0; i < Nx; i++) {
            double x = -xmax + i * dx_g;
            double f = A_init * exp(-x*x / (2.0*sig_use*sig_use));
            st.phi[a][i] = f * cos(phases[a]);
            st.vel[a][i] = -omega_est * f * sin(phases[a]);
        }
    compute_acc(&st, Nx);

    Diag d0 = diagnose(&st, Nx, core_r);
    res.E_init = d0.Et;
    printf("  E_init = %.4f  Epw = %.4f  Etr = %.4f\n",
           d0.Et, d0.Ep_pw, d0.Ep_triple);

    /* DFT/energy storage */
    int max_hist = 50000;
    double *phi_hist = malloc(max_hist * sizeof(double));
    double *t_hist   = malloc(max_hist * sizeof(double));
    int n_hist = 0;

    /* ===== Phase 1: Equilibrate (absorb transients) ===== */
    printf("  Equilibrating for t=%.0f ...\n", t_equil);
    int print_eq = Nt_equil / 10;
    if (print_eq < 1) print_eq = 1;
    int hist_every = Nt_equil / max_hist;
    if (hist_every < 1) hist_every = 1;

    for (int n = 0; n < Nt_equil; n++) {
        if (n % hist_every == 0 && n_hist < max_hist) {
            phi_hist[n_hist] = st.phi[0][ic];
            t_hist[n_hist] = n * dt;
            n_hist++;
        }
        if (n % print_eq == 0) {
            Diag d = diagnose(&st, Nx, core_r);
            printf("    equil t=%7.0f  E=%+.4f  pk=%.4f  fc=%.3f\n",
                   n*dt, d.Et, d.peak[0], d.fc);
        }
        step_verlet(&st, dt, damp, Nx);
    }

    /* Measure omega from equilibration */
    int dft_start_eq = n_hist * 3 / 4;
    double omega_eq = do_dft(phi_hist, t_hist, n_hist, dft_start_eq, 3.0*mass, 800, NULL);
    if (omega_eq > 0.05) omega_est = omega_eq;

    Diag d_eq = diagnose(&st, Nx, core_r);
    res.E_settled = d_eq.Et;
    printf("  After equil: omega=%.4f E=%.4f pk=%.4f fc=%.3f\n",
           omega_eq, d_eq.Et, d_eq.peak[0], d_eq.fc);

    /* ===== Phase 2: Main measurement run ===== */
    printf("  Main run for t=%.0f ...\n", t_run);

    n_hist = 0;
    int hist_every2 = Nt_run / max_hist;
    if (hist_every2 < 1) hist_every2 = 1;

    int max_ehist = 4000;
    double *E_hist  = malloc(max_ehist * sizeof(double));
    double *tE_hist = malloc(max_ehist * sizeof(double));
    int n_ehist = 0;
    int ehist_every = Nt_run / max_ehist;
    if (ehist_every < 1) ehist_every = 1;

    /* Time series file */
    char tspath[600];
    snprintf(tspath, sizeof(tspath), "%s/pt120_lam%.4f_ts.tsv", outdir, lam);
    FILE *fts = fopen(tspath, "w");
    if (!fts) { fprintf(stderr, "Cannot open %s\n", tspath); free_state(&st); return res; }
    fprintf(fts, "time\tphi1_0\tphi2_0\tphi3_0\tpeak1\tE_total\tE_pw\tE_triple\tf_core\t"
                 "dtheta12\tdtheta23\n");

    int rec_every = Nt_run / 20000;
    if (rec_every < 1) rec_every = 1;
    int print_every = Nt_run / 20;
    if (print_every < 1) print_every = 1;

    for (int n = 0; n <= Nt_run; n++) {
        double t = n * dt;

        if (n % hist_every2 == 0 && n_hist < max_hist) {
            phi_hist[n_hist] = st.phi[0][ic];
            t_hist[n_hist] = t;
            n_hist++;
        }

        if (n % ehist_every == 0 && n_ehist < max_ehist) {
            Diag d = diagnose(&st, Nx, core_r);
            E_hist[n_ehist] = d.Et;
            tE_hist[n_ehist] = t;
            n_ehist++;
        }

        int do_rec = (n % rec_every == 0);
        int do_print = (n % print_every == 0);

        if (do_rec || do_print) {
            Diag d = diagnose(&st, Nx, core_r);
            double dth12 = measure_phase_diff(st.phi[0][ic], st.vel[0][ic],
                                               st.phi[1][ic], st.vel[1][ic], omega_est);
            double dth23 = measure_phase_diff(st.phi[1][ic], st.vel[1][ic],
                                               st.phi[2][ic], st.vel[2][ic], omega_est);

            if (do_rec) {
                fprintf(fts, "%.4f\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6f\t"
                        "%.4f\t%.4f\n",
                        t, d.phi0[0], d.phi0[1], d.phi0[2],
                        d.peak[0], d.Et, d.Ep_pw, d.Ep_triple, d.fc,
                        dth12*180.0/M_PI, dth23*180.0/M_PI);
            }

            if (do_print) {
                printf("    t=%7.0f  pk=%.4f  E=%+.4f  Epw=%+.3f  fc=%.3f  "
                       "dth=(%.0f,%.0f)\n",
                       t, d.peak[0], d.Et, d.Ep_pw, d.fc,
                       dth12*180.0/M_PI, dth23*180.0/M_PI);
            }
        }

        if (n == Nt_run) break;
        step_verlet(&st, dt, damp, Nx);
    }
    fclose(fts);

    /* DFT: last half and last quarter */
    dft_start_eq = n_hist / 2;
    char dftpath[600];
    snprintf(dftpath, sizeof(dftpath), "%s/pt120_lam%.4f_spectrum.tsv", outdir, lam);
    FILE *fdft = fopen(dftpath, "w");
    fprintf(fdft, "omega\tpower\n");
    double omega_meas = do_dft(phi_hist, t_hist, n_hist, dft_start_eq, 3.0*mass, 800, fdft);
    fclose(fdft);

    int dft_late = n_hist * 3 / 4;
    double omega_late = do_dft(phi_hist, t_hist, n_hist, dft_late, 3.0*mass, 800, NULL);

    res.omega = (omega_late > 0.01) ? omega_late : omega_meas;
    res.three_omega = 3.0 * res.omega;
    res.sub_thresh = (res.three_omega < mass) ? 1 : 0;

    /* dE/dt */
    if (n_ehist > 20) {
        int n20 = n_ehist / 5;
        res.dEdt_early = (E_hist[n20] - E_hist[0]) / (tE_hist[n20] - tE_hist[0]);
        int i0 = n_ehist - n20;
        res.dEdt_late = (E_hist[n_ehist-1] - E_hist[i0]) / (tE_hist[n_ehist-1] - tE_hist[i0]);
    }

    /* Final */
    Diag df = diagnose(&st, Nx, core_r);
    res.E_final = df.Et;
    res.E_retained = (fabs(res.E_settled) > 1e-10) ? df.Et / res.E_settled : 0.0;
    res.fc_final = df.fc;
    res.peak_final = df.peak[0];

    res.phase12_final = measure_phase_diff(st.phi[0][ic], st.vel[0][ic],
                                            st.phi[1][ic], st.vel[1][ic],
                                            res.omega) * 180.0/M_PI;
    res.phase23_final = measure_phase_diff(st.phi[1][ic], st.vel[1][ic],
                                            st.phi[2][ic], st.vel[2][ic],
                                            res.omega) * 180.0/M_PI;

    printf("\n  === RESULT lam=%.4f ===\n", lam);
    printf("  omega=%.4f  3omega=%.4f  %s m=%.3f\n",
           res.omega, res.three_omega,
           res.sub_thresh ? "<" : ">=", mass);
    printf("  E: settled=%.4f final=%.4f ratio=%.4f\n",
           res.E_settled, res.E_final, res.E_retained);
    printf("  dE/dt: early=%.2e late=%.2e\n", res.dEdt_early, res.dEdt_late);
    printf("  fc=%.3f pk=%.4f phase=(%.0f,%.0f)\n",
           res.fc_final, res.peak_final, res.phase12_final, res.phase23_final);

    free_state(&st);
    free(phi_hist); free(t_hist);
    free(E_hist); free(tE_hist);

    return res;
}

int main(int argc, char **argv)
{
    parse_args(argc, argv);

    dx_g  = 2.0 * xmax / (Nx - 1);
    dx2_g = dx_g * dx_g;
    m2_g  = mass * mass;

    printf("pt120: Pairwise + Triple Product 120-degree oscillon scanner\n");
    printf("  mu=%.1f kappa=%.1f mass=%.3f A=%.3f sigma=%.3f\n",
           mu, kappa, mass, A_init, sig);
    printf("  Nx=%d xmax=%.1f dx=%.5f\n", Nx, xmax, dx_g);
    printf("  t_equil=%.0f  t_run=%.0f\n", t_equil, t_run);
    printf("  Threshold: 3*omega < m for NO radiation\n\n");

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

    double *lambdas;
    int nlam;
    if (lam_single >= 0.0) {
        lambdas = &lam_single;
        nlam = 1;
    } else {
        lambdas = lam_list;
        nlam = n_lam;
    }

    Result *results = malloc(nlam * sizeof(Result));

    for (int il = 0; il < nlam; il++) {
        results[il] = run_one(lambdas[il], damp);
    }

    /* Summary table */
    printf("\n\n================================================================\n");
    printf("                    SUMMARY TABLE\n");
    printf("================================================================\n");
    printf("lambda  m_A     omega   3omega  3w<m  E_settl   E_final   E_rat   dE/dt_late  fc    pk    ph12  ph23\n");
    printf("------  ------  ------  ------  ----  --------  --------  ------  ----------  ----  ----  ----  ----\n");

    for (int il = 0; il < nlam; il++) {
        Result *r = &results[il];
        printf("%.4f  %.4f  %.4f  %.4f  %s  %+8.3f  %+8.3f  %6.4f  %+.2e  %.3f %.4f  %4.0f  %4.0f\n",
               r->lambda, r->m_anti, r->omega, r->three_omega,
               r->sub_thresh ? "YES " : " NO ",
               r->E_settled, r->E_final, r->E_retained,
               r->dEdt_late,
               r->fc_final, r->peak_final,
               r->phase12_final, r->phase23_final);
    }

    /* Summary TSV */
    char sumpath[600];
    snprintf(sumpath, sizeof(sumpath), "%s/pt120_summary.tsv", outdir);
    FILE *fsum = fopen(sumpath, "w");
    fprintf(fsum, "lambda\tm2_anti\tm_anti\tomega\tthree_omega\tsub_thresh\t"
                  "E_init\tE_settled\tE_final\tE_retained\t"
                  "dEdt_early\tdEdt_late\tfc_final\tpeak_final\tphase12\tphase23\n");
    for (int il = 0; il < nlam; il++) {
        Result *r = &results[il];
        fprintf(fsum, "%.4f\t%.6f\t%.6f\t%.6f\t%.6f\t%d\t"
                "%.6f\t%.6f\t%.6f\t%.6f\t%.6e\t%.6e\t%.6f\t%.6f\t%.2f\t%.2f\n",
                r->lambda, r->m2_anti, r->m_anti, r->omega, r->three_omega,
                r->sub_thresh,
                r->E_init, r->E_settled, r->E_final, r->E_retained,
                r->dEdt_early, r->dEdt_late, r->fc_final, r->peak_final,
                r->phase12_final, r->phase23_final);
    }
    fclose(fsum);

    printf("\nSummary: %s\n", sumpath);

    free(damp);
    free(results);
    return 0;
}
