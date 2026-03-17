/*
 * uud120.c — UUD proton-like 120° oscillon with per-field masses
 *
 * Three massive scalars with:
 *   - Per-field masses: m_1²=m_2²=m_U², m_3²=m_D²
 *   - Pairwise coupling lambda (off-diagonal mass matrix)
 *   - Triple product: V = (mu/2)*P²/(1+kappa*P²), P = phi_1*phi_2*phi_3
 *
 * EOM: ddot(phi_a) = lapl(phi_a) - m_a² phi_a - lambda*(phi_b + phi_c) - dV/dphi_a
 *
 * Protocol:
 *   1. Equilibrate symmetric oscillon (triple product only, no pairwise, equal masses)
 *   2. Rotate to 120° phases, set m_3=m_D, enable pairwise lambda
 *   3. Evolve and measure stability
 *
 * Compile: gcc -O3 -Wall -o uud120 v24/120uud/src/uud120.c -lm
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static double mu      = -20.0;
static double kappa   = 20.0;
static double m_U     = 1.0;
static double lambda  = 0.05;
static double A_init  = 0.8;
static double sigma_g = 3.0;
static int    Nx      = 4000;
static double xmax    = 100.0;
static double t_equil = 3000.0;
static double t_run   = 5000.0;
static char   outdir[512] = "v24/120uud/data";

static void parse_args(int argc, char **argv)
{
    for (int i = 1; i < argc - 1; i += 2) {
        if      (!strcmp(argv[i], "-mu"))      mu      = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-kappa"))   kappa   = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-mU"))      m_U     = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-lambda"))  lambda  = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-A"))       A_init  = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-sigma"))   sigma_g = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-Nx"))      Nx      = atoi(argv[i+1]);
        else if (!strcmp(argv[i], "-xmax"))    xmax    = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-t_equil")) t_equil = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-t_run"))   t_run   = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-o"))       strncpy(outdir, argv[i+1], sizeof(outdir)-1);
        else { fprintf(stderr, "Unknown flag: %s\n", argv[i]); exit(1); }
    }
}

static double m2_field[3];
static double lambda_active;

static inline double force_triple(double p1, double p2, double p3, int a)
{
    double P = p1 * p2 * p3;
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
    double *phi[3], *vel[3], *acc[3];
} State;

static State alloc_state(int n)
{
    State s;
    for (int a = 0; a < 3; a++) {
        s.phi[a] = calloc(n, sizeof(double));
        s.vel[a] = calloc(n, sizeof(double));
        s.acc[a] = calloc(n, sizeof(double));
    }
    return s;
}

static void free_state(State *s)
{
    for (int a = 0; a < 3; a++) {
        free(s->phi[a]); free(s->vel[a]); free(s->acc[a]);
    }
}

static double dx_g, dx2_g;

static void compute_acc(State *s, int n)
{
    for (int a = 0; a < 3; a++) {
        s->acc[a][0] = s->acc[a][1] = s->acc[a][n-2] = s->acc[a][n-1] = 0;
        for (int i = 1; i < n - 1; i++) {
            double lapl = (s->phi[a][i+1] - 2.0*s->phi[a][i] + s->phi[a][i-1]) / dx2_g;
            double ft = force_triple(s->phi[0][i], s->phi[1][i], s->phi[2][i], a);
            /* Pairwise: -lambda*(sum of other fields) */
            double pw = 0;
            if (lambda_active > 1e-15) {
                for (int b = 0; b < 3; b++)
                    if (b != a) pw += s->phi[b][i];
            }
            s->acc[a][i] = lapl - m2_field[a]*s->phi[a][i] - lambda_active*pw + ft;
        }
    }
}

static void step_verlet(State *s, double dt, double *damp, int n)
{
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
    double Ek, Eg, Em, Ep, Epw, Et;
    double peak[3];
    double fc;
    double phi0[3], vel0[3];
} Diag;

static Diag diagnose(State *s, int n, double core_r)
{
    Diag d = {0};
    int ic = n / 2;
    double Ecore = 0, Eall = 0;

    for (int a = 0; a < 3; a++) {
        d.phi0[a] = s->phi[a][ic];
        d.vel0[a] = s->vel[a][ic];
    }

    for (int i = 1; i < n - 1; i++) {
        double x = -xmax + i * dx_g;
        for (int a = 0; a < 3; a++) {
            d.Ek += 0.5 * s->vel[a][i] * s->vel[a][i] * dx_g;
            double dp = (s->phi[a][i+1] - s->phi[a][i-1]) / (2.0*dx_g);
            d.Eg += 0.5 * dp * dp * dx_g;
            d.Em += 0.5 * m2_field[a] * s->phi[a][i] * s->phi[a][i] * dx_g;
            if (fabs(s->phi[a][i]) > d.peak[a]) d.peak[a] = fabs(s->phi[a][i]);
        }
        double P = s->phi[0][i] * s->phi[1][i] * s->phi[2][i];
        double P2 = P * P;
        d.Ep += 0.5 * mu * P2 / (1.0 + kappa * P2) * dx_g;

        double vpw = 0;
        for (int a = 0; a < 3; a++)
            for (int b = a+1; b < 3; b++)
                vpw += s->phi[a][i] * s->phi[b][i];
        d.Epw += lambda_active * vpw * dx_g;

        double e = 0.5 * mu * P2 / (1.0 + kappa * P2) + lambda_active * vpw;
        for (int a = 0; a < 3; a++) {
            e += 0.5*s->vel[a][i]*s->vel[a][i];
            double dp = (s->phi[a][i+1] - s->phi[a][i-1]) / (2.0*dx_g);
            e += 0.5*dp*dp + 0.5*m2_field[a]*s->phi[a][i]*s->phi[a][i];
        }
        Eall += e * dx_g;
        if (fabs(x) < core_r) Ecore += e * dx_g;
    }
    d.Et = d.Ek + d.Eg + d.Em + d.Ep + d.Epw;
    d.fc = (fabs(Eall) > 1e-20) ? Ecore / Eall : 0.0;
    return d;
}

static double measure_phase_diff(double phi_a, double vel_a,
                                  double phi_b, double vel_b, double omega)
{
    if (omega < 0.01) return 0.0;
    double pa = atan2(-vel_a/omega, phi_a);
    double pb = atan2(-vel_b/omega, phi_b);
    double diff = pa - pb;
    while (diff > M_PI) diff -= 2.0*M_PI;
    while (diff < -M_PI) diff += 2.0*M_PI;
    return diff;
}

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
        if (pw > peak_pow && omega > 0.01) { peak_pow = pw; peak_om = omega; }
    }
    return peak_om;
}

int main(int argc, char **argv)
{
    parse_args(argc, argv);

    dx_g  = 2.0 * xmax / (Nx - 1);
    dx2_g = dx_g * dx_g;

    double kmax = M_PI / dx_g;
    double m2eff_max = m_U*m_U + 2.0*lambda;
    double dt = 0.4 * 2.0 / sqrt(kmax*kmax + m2eff_max);
    int Nt_equil = (int)(t_equil / dt) + 1;
    int Nt_run   = (int)(t_run / dt) + 1;

    printf("uud120: UUD proton-like 120-deg oscillon\n");
    printf("  m_U=%.3f, lambda=%.3f, mu=%.1f, kappa=%.1f\n", m_U, lambda, mu, kappa);
    printf("  A=%.3f, sigma=%.3f, Nx=%d, xmax=%.1f, dx=%.5f, dt=%.6f\n",
           A_init, sigma_g, Nx, xmax, dx_g, dt);
    printf("  t_equil=%.0f (Nt=%d), t_run=%.0f (Nt=%d)\n",
           t_equil, Nt_equil, t_run, Nt_run);

    /* Absorbing boundary */
    double *damp = malloc(Nx * sizeof(double));
    double x_abs = xmax * 0.75;
    for (int i = 0; i < Nx; i++) {
        double x = -xmax + i * dx_g;
        double ax = fabs(x);
        damp[i] = (ax > x_abs) ? 1.0 - 0.98*((ax-x_abs)/(xmax-x_abs))*((ax-x_abs)/(xmax-x_abs)) : 1.0;
    }

    /* ===== Phase 1: Equilibrate symmetric oscillon (triple product only) ===== */
    printf("\n=== Phase 1: Equilibrating symmetric oscillon ===\n");
    m2_field[0] = m2_field[1] = m2_field[2] = m_U * m_U;
    lambda_active = 0.0;

    State sym = alloc_state(Nx);
    for (int a = 0; a < 3; a++)
        for (int i = 0; i < Nx; i++) {
            double x = -xmax + i * dx_g;
            sym.phi[a][i] = A_init * exp(-x*x / (2.0*sigma_g*sigma_g));
        }
    compute_acc(&sym, Nx);

    double core_r = 3.0 * sigma_g;
    int ic = Nx / 2;

    int max_hist = 50000;
    double *phi_eq_hist = malloc(max_hist * sizeof(double));
    double *t_eq_hist   = malloc(max_hist * sizeof(double));
    int n_eq = 0;
    int eq_rec = Nt_equil / max_hist; if (eq_rec < 1) eq_rec = 1;
    int print_eq = Nt_equil / 20; if (print_eq < 1) print_eq = 1;

    for (int n = 0; n < Nt_equil; n++) {
        double t = n * dt;
        if (n % eq_rec == 0 && n_eq < max_hist) {
            phi_eq_hist[n_eq] = sym.phi[0][ic];
            t_eq_hist[n_eq] = t;
            n_eq++;
        }
        if (n % print_eq == 0) {
            Diag d = diagnose(&sym, Nx, core_r);
            printf("  equil t=%7.1f  E=%.4f  Ep=%.4f  fc=%.3f  pk=%.3f\n",
                   t, d.Et, d.Ep, d.fc, d.peak[0]);
        }
        step_verlet(&sym, dt, damp, Nx);
    }

    int dft_start_eq = n_eq * 3 / 4;
    double omega_eq = do_dft(phi_eq_hist, t_eq_hist, n_eq, dft_start_eq,
                             3.0*m_U, 500, NULL);
    Diag d_eq = diagnose(&sym, Nx, core_r);
    double m_eff_120 = (m_U*m_U > lambda) ? sqrt(m_U*m_U - lambda) : 0.0;

    printf("\nEquilibrated: omega=%.4f, m=%.4f, m_eff(120)=%.4f\n",
           omega_eq, m_U, m_eff_120);
    printf("  3*omega=%.4f, omega/m_eff(120)=%.4f\n",
           3.0*omega_eq, (m_eff_120 > 0.01) ? omega_eq/m_eff_120 : 999.0);
    printf("  E=%.4f, fc=%.3f, pk=%.3f\n", d_eq.Et, d_eq.fc, d_eq.peak[0]);
    printf("  omega < m_eff(120)? %s\n",
           omega_eq < m_eff_120 ? "YES (stable 120-deg)" : "NO (120-deg will radiate)");

    free(phi_eq_hist); free(t_eq_hist);

    double omega_use = (omega_eq > 0.01) ? omega_eq : 0.5;

    /* ===== Phase 2+3: m_D scan ===== */
    double mD_vals[] = {1.00, 0.98, 0.96, 0.94, 0.92, 0.90, 0.85};
    int n_scan = sizeof(mD_vals) / sizeof(mD_vals[0]);

    char sumpath[600];
    snprintf(sumpath, sizeof(sumpath), "%s/uud_summary.tsv", outdir);
    FILE *fsum = fopen(sumpath, "w");
    fprintf(fsum, "m_D\tm_D_over_mU\tE_init\tE_final\tE_retained\tfc_final\t"
                  "omega1\tomega3\tdth12_deg\tdth23_deg\tsurvived\n");

    for (int si = 0; si < n_scan; si++) {
        double mD = mD_vals[si];
        double m_eff_D = (mD*mD > lambda) ? sqrt(mD*mD - lambda) : 0.0;

        printf("\n========================================\n");
        printf("UUD: m_U=%.3f, m_D=%.3f, lambda=%.3f\n", m_U, mD, lambda);
        printf("  m_eff_U(120)=%.4f, m_eff_D(120)=%.4f\n", m_eff_120, m_eff_D);
        printf("========================================\n");

        /* Set masses */
        m2_field[0] = m_U * m_U;
        m2_field[1] = m_U * m_U;
        m2_field[2] = mD * mD;
        lambda_active = lambda;

        /* Rotate equilibrated symmetric state to 120° */
        State uud = alloc_state(Nx);
        for (int i = 0; i < Nx; i++) {
            double p1 = sym.phi[0][i];
            double v1 = sym.vel[0][i];
            double f_amp = sqrt(p1*p1 + v1*v1/(omega_use*omega_use));
            double theta = atan2(-v1/omega_use, p1);

            uud.phi[0][i] = f_amp * cos(theta);
            uud.vel[0][i] = -omega_use * f_amp * sin(theta);

            double th2 = theta + 2.0*M_PI/3.0;
            uud.phi[1][i] = f_amp * cos(th2);
            uud.vel[1][i] = -omega_use * f_amp * sin(th2);

            double th3 = theta + 4.0*M_PI/3.0;
            uud.phi[2][i] = f_amp * cos(th3);
            uud.vel[2][i] = -omega_use * f_amp * sin(th3);
        }
        compute_acc(&uud, Nx);

        Diag d_init = diagnose(&uud, Nx, core_r);
        printf("  E_init=%.4f, Ep=%.4f, Epw=%.4f, fc=%.3f\n",
               d_init.Et, d_init.Ep, d_init.Epw, d_init.fc);

        /* Time series */
        char tspath[600];
        snprintf(tspath, sizeof(tspath), "%s/uud_mD%.2f_ts.tsv", outdir, mD);
        FILE *fts = fopen(tspath, "w");
        fprintf(fts, "time\tphi1_0\tphi2_0\tphi3_0\tpeak1\tpeak2\tpeak3\t"
                     "E_kin\tE_grad\tE_mass\tE_pot\tE_pw\tE_total\tf_core\t"
                     "dth12\tdth23\n");

        double *phi1h = malloc(max_hist * sizeof(double));
        double *phi3h = malloc(max_hist * sizeof(double));
        double *th    = malloc(max_hist * sizeof(double));
        int ndft = 0;
        int rec_ev = Nt_run / 20000; if (rec_ev < 1) rec_ev = 1;
        int prt_ev = Nt_run / 40;    if (prt_ev < 1) prt_ev = 1;
        int dft_ev = Nt_run / max_hist; if (dft_ev < 1) dft_ev = 1;

        for (int n = 0; n <= Nt_run; n++) {
            double t = n * dt;

            if (n % dft_ev == 0 && ndft < max_hist) {
                phi1h[ndft] = uud.phi[0][ic];
                phi3h[ndft] = uud.phi[2][ic];
                th[ndft] = t;
                ndft++;
            }

            int do_rec = (n % rec_ev == 0);
            int do_prt = (n % prt_ev == 0);

            if (do_rec || do_prt) {
                Diag d = diagnose(&uud, Nx, core_r);
                double dth12 = measure_phase_diff(uud.phi[0][ic], uud.vel[0][ic],
                                                   uud.phi[1][ic], uud.vel[1][ic], omega_use);
                double dth23 = measure_phase_diff(uud.phi[1][ic], uud.vel[1][ic],
                                                   uud.phi[2][ic], uud.vel[2][ic], omega_use);
                if (do_rec)
                    fprintf(fts, "%.6f\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t"
                            "%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6f\t%.6f\t%.6f\n",
                            t, d.phi0[0], d.phi0[1], d.phi0[2],
                            d.peak[0], d.peak[1], d.peak[2],
                            d.Ek, d.Eg, d.Em, d.Ep, d.Epw, d.Et, d.fc,
                            dth12, dth23);
                if (do_prt)
                    printf("  t=%7.1f  pk=(%.3f,%.3f,%.3f)  E=%+.4f  fc=%.3f  dth=(%.0f,%.0f)\n",
                           t, d.peak[0], d.peak[1], d.peak[2],
                           d.Et, d.fc, dth12*180/M_PI, dth23*180/M_PI);
            }

            if (n == Nt_run) break;
            step_verlet(&uud, dt, damp, Nx);
        }
        fclose(fts);

        Diag d_final = diagnose(&uud, Nx, core_r);
        double E_ret = (fabs(d_init.Et) > 1e-10) ? d_final.Et / d_init.Et : 0.0;
        double dth12_f = measure_phase_diff(uud.phi[0][ic], uud.vel[0][ic],
                                             uud.phi[1][ic], uud.vel[1][ic], omega_use);
        double dth23_f = measure_phase_diff(uud.phi[1][ic], uud.vel[1][ic],
                                             uud.phi[2][ic], uud.vel[2][ic], omega_use);

        int dft_start = ndft / 2;
        char dp[600];
        snprintf(dp, sizeof(dp), "%s/uud_mD%.2f_phi1_spectrum.tsv", outdir, mD);
        FILE *fd = fopen(dp, "w"); fprintf(fd, "omega\tpower\n");
        double om1 = do_dft(phi1h, th, ndft, dft_start, 3.0*m_U, 500, fd);
        fclose(fd);

        snprintf(dp, sizeof(dp), "%s/uud_mD%.2f_phi3_spectrum.tsv", outdir, mD);
        fd = fopen(dp, "w"); fprintf(fd, "omega\tpower\n");
        double om3 = do_dft(phi3h, th, ndft, dft_start, 3.0*m_U, 500, fd);
        fclose(fd);

        /* Survival: >50% energy retained, >30% in core, phase ~120 deg */
        int survived = (E_ret > 0.50 && d_final.fc > 0.30) ? 1 : 0;

        printf("\n  RESULT mD=%.3f: E=%+.4f -> %+.4f (%.1f%%), fc=%.3f\n",
               mD, d_init.Et, d_final.Et, 100.0*E_ret, d_final.fc);
        printf("  omega1=%.4f, omega3=%.4f, dth=(%.1f,%.1f)  %s\n",
               om1, om3, dth12_f*180/M_PI, dth23_f*180/M_PI,
               survived ? "SURVIVED" : "DIED");

        fprintf(fsum, "%.3f\t%.3f\t%.6f\t%.6f\t%.6f\t%.6f\t%.4f\t%.4f\t%.4f\t%.4f\t%d\n",
                mD, mD/m_U, d_init.Et, d_final.Et, E_ret,
                d_final.fc, om1, om3,
                dth12_f*180/M_PI, dth23_f*180/M_PI, survived);
        fflush(fsum);

        free_state(&uud);
        free(phi1h); free(phi3h); free(th);
    }

    fclose(fsum);
    printf("\nSummary: %s\n", sumpath);
    free_state(&sym); free(damp);
    return 0;
}
