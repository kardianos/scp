/*
 * udd120.c — UDD (neutron-like) 120-degree oscillon
 *
 * Three massive scalars with triple-product coupling + pairwise coupling.
 * Field 1 = "up" (mass m_U), Fields 2,3 = "down" (mass m_D).
 * 120-degree phase initialization: phases 0, 2pi/3, 4pi/3.
 *
 * EOM:
 *   d²phi_a/dt² = d²phi_a/dx² - m_a² phi_a - lambda*(sum_{b!=a} phi_b) - dV_triple/dphi_a
 *
 * V_triple = (mu/2) P² / (1 + kappa P²),  P = phi_1 phi_2 phi_3
 *
 * Mode: UDD (field 1 = up, fields 2,3 = down)
 *   or  UUD (fields 1,2 = up, field 3 = down) via -mode flag
 *
 * Compile: gcc -O3 -Wall -o udd120 src/udd120.c -lm
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Parameters */
static double mu      = -20.0;
static double kappa   = 20.0;
static double m_U     = 1.0;
static double m_D     = 1.0;
static double lambda  = 0.85;
static double A_init  = 0.8;
static double sigma_g = 3.0;
static int    Nx      = 4000;
static double xmax    = 100.0;
static double tfinal  = 10000.0;
static int    mode    = 0;  /* 0 = UDD (neutron), 1 = UUD (proton) */
static char   outdir[512] = "v24/120udd/data";
static char   prefix[64]  = "udd";

static void parse_args(int argc, char **argv)
{
    for (int i = 1; i < argc - 1; i += 2) {
        if      (!strcmp(argv[i], "-mu"))      mu      = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-kappa"))   kappa   = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-mU"))      m_U     = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-mD"))      m_D     = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-lambda"))  lambda  = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-A"))       A_init  = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-sigma"))   sigma_g = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-Nx"))      Nx      = atoi(argv[i+1]);
        else if (!strcmp(argv[i], "-xmax"))    xmax    = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-tfinal"))  tfinal  = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-mode"))    mode    = atoi(argv[i+1]);
        else if (!strcmp(argv[i], "-o"))       strncpy(outdir, argv[i+1], sizeof(outdir)-1);
        else if (!strcmp(argv[i], "-prefix"))  strncpy(prefix, argv[i+1], sizeof(prefix)-1);
        else { fprintf(stderr, "Unknown flag: %s\n", argv[i]); exit(1); }
    }
}

/* Per-field mass squared */
static double m2[3];

/* -dV/dphi_a where V = (mu/2)P^2/(1+kappa P^2), P = phi1 phi2 phi3 */
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

/* Pairwise coupling force: -d/dphi_a [ (lambda/2) sum_{b<c} (phi_b - phi_c)^2 ]
 * = ... actually the simpler form: -lambda * sum_{b!=a} phi_b
 * This comes from V_pw = lambda * sum_{a<b} phi_a phi_b
 * dV_pw/dphi_a = lambda * sum_{b!=a} phi_b
 * force = -dV_pw/dphi_a = -lambda * sum_{b!=a} phi_b
 */
static inline double force_pairwise(double p1, double p2, double p3, int a)
{
    switch (a) {
        case 0: return -lambda * (p2 + p3);
        case 1: return -lambda * (p1 + p3);
        case 2: return -lambda * (p1 + p2);
        default: return 0.0;
    }
}

static double dx_g, dx2_g;
static double *phi[3], *vel[3], *acc_f[3];

static void compute_acc(void)
{
    for (int a = 0; a < 3; a++) {
        acc_f[a][0] = acc_f[a][1] = acc_f[a][Nx-2] = acc_f[a][Nx-1] = 0;
        for (int i = 1; i < Nx - 1; i++) {
            double lapl = (phi[a][i+1] - 2.0*phi[a][i] + phi[a][i-1]) / dx2_g;
            double ft = force_triple(phi[0][i], phi[1][i], phi[2][i], a);
            double fp = force_pairwise(phi[0][i], phi[1][i], phi[2][i], a);
            acc_f[a][i] = lapl - m2[a]*phi[a][i] + ft + fp;
        }
    }
}

typedef struct {
    double Ek, Eg, Em, Ep, Epw, Et;
    double peak[3];
    double fc;
    double phi0[3];
} Diag;

static Diag diagnose(double core_r)
{
    Diag d = {0};
    int ic = Nx / 2;
    double Ecore = 0, Eall = 0;

    for (int a = 0; a < 3; a++)
        d.phi0[a] = phi[a][ic];

    for (int i = 1; i < Nx - 1; i++) {
        double x = -xmax + i * dx_g;
        for (int a = 0; a < 3; a++) {
            d.Ek += 0.5 * vel[a][i] * vel[a][i] * dx_g;
            double dp = (phi[a][i+1] - phi[a][i-1]) / (2.0*dx_g);
            d.Eg += 0.5 * dp * dp * dx_g;
            d.Em += 0.5 * m2[a] * phi[a][i] * phi[a][i] * dx_g;
            if (fabs(phi[a][i]) > d.peak[a]) d.peak[a] = fabs(phi[a][i]);
        }
        /* Triple-product potential */
        double P = phi[0][i] * phi[1][i] * phi[2][i];
        double P2 = P * P;
        double V = 0.5 * mu * P2 / (1.0 + kappa * P2);
        d.Ep += V * dx_g;

        /* Pairwise potential: lambda * sum_{a<b} phi_a phi_b */
        double Vpw = lambda * (phi[0][i]*phi[1][i] + phi[0][i]*phi[2][i] + phi[1][i]*phi[2][i]);
        d.Epw += Vpw * dx_g;

        /* total density */
        double e = V + Vpw;
        for (int a = 0; a < 3; a++) {
            e += 0.5*vel[a][i]*vel[a][i];
            double dp = (phi[a][i+1] - phi[a][i-1]) / (2.0*dx_g);
            e += 0.5*dp*dp + 0.5*m2[a]*phi[a][i]*phi[a][i];
        }
        Eall += e * dx_g;
        if (fabs(x) < core_r) Ecore += e * dx_g;
    }
    d.Et = d.Ek + d.Eg + d.Em + d.Ep + d.Epw;
    d.fc = (Eall > 1e-20) ? Ecore / Eall : 0.0;
    return d;
}

/* Measure phase difference */
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

/* DFT */
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

    /* Set per-field masses based on mode */
    if (mode == 0) {
        /* UDD: field 1 = up, fields 2,3 = down */
        m2[0] = m_U * m_U;
        m2[1] = m_D * m_D;
        m2[2] = m_D * m_D;
        snprintf(prefix, sizeof(prefix), "udd");
    } else {
        /* UUD: fields 1,2 = up, field 3 = down */
        m2[0] = m_U * m_U;
        m2[1] = m_U * m_U;
        m2[2] = m_D * m_D;
        snprintf(prefix, sizeof(prefix), "uud");
    }

    dx_g = 2.0 * xmax / (Nx - 1);
    dx2_g = dx_g * dx_g;

    /* CFL: use minimum mass for stability */
    double m_min = (m_U < m_D) ? m_U : m_D;
    double m2_min = m_min * m_min;
    double kmax = M_PI / dx_g;
    double dt = 0.8 * 2.0 / sqrt(kmax*kmax + m2_min + 4.0*lambda);
    /* add pairwise coupling to stability bound */
    int Nt = (int)(tfinal / dt) + 1;

    const char *mode_str = (mode == 0) ? "UDD (neutron)" : "UUD (proton)";
    printf("%s 120-degree oscillon\n", mode_str);
    printf("  mu=%.1f kappa=%.1f m_U=%.3f m_D=%.3f lambda=%.3f\n",
           mu, kappa, m_U, m_D, lambda);
    printf("  Nx=%d xmax=%.1f dx=%.5f dt=%.6f Nt=%d tfinal=%.0f\n",
           Nx, xmax, dx_g, dt, Nt, tfinal);

    for (int a = 0; a < 3; a++) {
        phi[a]   = calloc(Nx, sizeof(double));
        vel[a]   = calloc(Nx, sizeof(double));
        acc_f[a] = calloc(Nx, sizeof(double));
    }

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

    /* Initialize 120-degree state
     * Mode UDD: phi_1 = A_U * g(x) * cos(0)       [up, phase 0]
     *           phi_2 = A_D * g(x) * cos(2pi/3)    [down, phase 2pi/3]
     *           phi_3 = A_D * g(x) * cos(4pi/3)    [down, phase 4pi/3]
     * Mode UUD: phi_1 = A_U * g(x) * cos(0)        [up, phase 0]
     *           phi_2 = A_U * g(x) * cos(2pi/3)    [up, phase 2pi/3]
     *           phi_3 = A_D * g(x) * cos(4pi/3)    [down, phase 4pi/3]
     *
     * Velocities from d/dt[A*g(x)*cos(omega*t + theta)] at t=0:
     *   vel = -omega * A * g(x) * sin(theta)
     * Use omega_a = sqrt(m_a^2 - lambda) if real, else m_a
     */
    double omega_U = (m_U*m_U > lambda) ? sqrt(m_U*m_U - lambda) : m_U;
    double omega_D = (m_D*m_D > lambda) ? sqrt(m_D*m_D - lambda) : m_D;

    double phases[3], amps[3], omegas[3];
    phases[0] = 0.0;
    phases[1] = 2.0*M_PI/3.0;
    phases[2] = 4.0*M_PI/3.0;

    if (mode == 0) {
        /* UDD */
        amps[0] = A_init; omegas[0] = omega_U;
        amps[1] = A_init; omegas[1] = omega_D;
        amps[2] = A_init; omegas[2] = omega_D;
    } else {
        /* UUD */
        amps[0] = A_init; omegas[0] = omega_U;
        amps[1] = A_init; omegas[1] = omega_U;
        amps[2] = A_init; omegas[2] = omega_D;
    }

    printf("  omega_U=%.4f omega_D=%.4f\n", omega_U, omega_D);
    printf("  phases: (%.1f, %.1f, %.1f) deg\n",
           phases[0]*180/M_PI, phases[1]*180/M_PI, phases[2]*180/M_PI);

    for (int a = 0; a < 3; a++)
        for (int i = 0; i < Nx; i++) {
            double x = -xmax + i * dx_g;
            double g = exp(-x * x / (2.0 * sigma_g * sigma_g));
            phi[a][i] = amps[a] * g * cos(phases[a]);
            vel[a][i] = -omegas[a] * amps[a] * g * sin(phases[a]);
        }

    compute_acc();

    double core_r = 3.0 * sigma_g;

    /* Output: time series */
    char tspath[600];
    snprintf(tspath, sizeof(tspath), "%s/%s_mD%.2f_ts.tsv", outdir, prefix, m_D);
    FILE *fts = fopen(tspath, "w");
    if (!fts) { fprintf(stderr, "Cannot open %s\n", tspath); return 1; }
    fprintf(fts, "time\tphi1_0\tphi2_0\tphi3_0\tpeak1\tpeak2\tpeak3\t"
                 "E_kin\tE_grad\tE_mass\tE_pot\tE_pairwise\tE_total\tf_core\t"
                 "dtheta12\tdtheta13\tdtheta23\n");

    /* DFT storage */
    int max_dft = 50000;
    double *phi1_hist = malloc(max_dft * sizeof(double));
    double *phi2_hist = malloc(max_dft * sizeof(double));
    double *phi3_hist = malloc(max_dft * sizeof(double));
    double *P_hist    = malloc(max_dft * sizeof(double));
    double *t_hist    = malloc(max_dft * sizeof(double));
    int n_dft = 0;

    int rec_every = Nt / 20000;
    if (rec_every < 1) rec_every = 1;
    int print_every = Nt / 40;
    if (print_every < 1) print_every = 1;
    int dft_every = Nt / max_dft;
    if (dft_every < 1) dft_every = 1;

    int ic = Nx / 2;

    /* Track early and late energy for summary */
    double E_init = 0, E_t1000 = 0, E_t5000 = 0, E_final = 0;
    double fc_init = 0, fc_final = 0;
    double peak_init[3] = {0}, peak_final[3] = {0};
    int got_t1000 = 0, got_t5000 = 0;

    for (int n = 0; n <= Nt; n++) {
        double t = n * dt;

        /* DFT sampling */
        if (n % dft_every == 0 && n_dft < max_dft) {
            phi1_hist[n_dft] = phi[0][ic];
            phi2_hist[n_dft] = phi[1][ic];
            phi3_hist[n_dft] = phi[2][ic];
            P_hist[n_dft]    = phi[0][ic] * phi[1][ic] * phi[2][ic];
            t_hist[n_dft] = t;
            n_dft++;
        }

        int do_rec   = (n % rec_every == 0);
        int do_print = (n % print_every == 0);

        if (do_rec || do_print || n == 0) {
            Diag d = diagnose(core_r);

            /* Use average omega for phase measurement */
            double omega_avg = (omegas[0] + omegas[1] + omegas[2]) / 3.0;
            double dth12 = measure_phase_diff(phi[0][ic], vel[0][ic],
                                               phi[1][ic], vel[1][ic], omega_avg);
            double dth13 = measure_phase_diff(phi[0][ic], vel[0][ic],
                                               phi[2][ic], vel[2][ic], omega_avg);
            double dth23 = measure_phase_diff(phi[1][ic], vel[1][ic],
                                               phi[2][ic], vel[2][ic], omega_avg);

            if (do_rec)
                fprintf(fts, "%.6f\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t"
                        "%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6f\t"
                        "%.6f\t%.6f\t%.6f\n",
                        t, d.phi0[0], d.phi0[1], d.phi0[2],
                        d.peak[0], d.peak[1], d.peak[2],
                        d.Ek, d.Eg, d.Em, d.Ep, d.Epw, d.Et, d.fc,
                        dth12, dth13, dth23);

            if (do_print)
                printf("  t=%7.1f  p0=(%+.3f,%+.3f,%+.3f)  pk=(%.3f,%.3f,%.3f)  "
                       "E=%+.4f  Ep=%+.4f  Epw=%+.4f  fc=%.3f  "
                       "dth=(%.0f,%.0f,%.0f)\n",
                       t, d.phi0[0], d.phi0[1], d.phi0[2],
                       d.peak[0], d.peak[1], d.peak[2],
                       d.Et, d.Ep, d.Epw, d.fc,
                       dth12*180/M_PI, dth13*180/M_PI, dth23*180/M_PI);

            if (n == 0) {
                E_init = d.Et; fc_init = d.fc;
                for (int a = 0; a < 3; a++) peak_init[a] = d.peak[a];
            }
            if (!got_t1000 && t >= 1000.0) { E_t1000 = d.Et; got_t1000 = 1; }
            if (!got_t5000 && t >= 5000.0) { E_t5000 = d.Et; got_t5000 = 1; }
            E_final = d.Et; fc_final = d.fc;
            for (int a = 0; a < 3; a++) peak_final[a] = d.peak[a];
        }

        if (n == Nt) break;

        /* Velocity Verlet */
        for (int a = 0; a < 3; a++)
            for (int i = 1; i < Nx - 1; i++)
                vel[a][i] += 0.5 * dt * acc_f[a][i];
        for (int a = 0; a < 3; a++)
            for (int i = 1; i < Nx - 1; i++)
                phi[a][i] += dt * vel[a][i];
        compute_acc();
        for (int a = 0; a < 3; a++)
            for (int i = 1; i < Nx - 1; i++)
                vel[a][i] += 0.5 * dt * acc_f[a][i];

        /* absorbing boundary */
        for (int a = 0; a < 3; a++)
            for (int i = 0; i < Nx; i++) {
                vel[a][i] *= damp[i];
                phi[a][i] *= damp[i];
            }
    }

    fclose(fts);

    /* DFT analysis */
    printf("\n=== DFT Analysis ===\n");
    int dft_start = n_dft / 2;

    char dft_path[600];
    snprintf(dft_path, sizeof(dft_path), "%s/%s_mD%.2f_phi1_spectrum.tsv", outdir, prefix, m_D);
    FILE *fdft = fopen(dft_path, "w");
    fprintf(fdft, "omega\tpower\n");
    double omega_1 = do_dft(phi1_hist, t_hist, n_dft, dft_start, 3.0*m_U, 500, fdft);
    fclose(fdft);
    printf("  phi_1 (up) peak omega = %.4f\n", omega_1);

    snprintf(dft_path, sizeof(dft_path), "%s/%s_mD%.2f_phi2_spectrum.tsv", outdir, prefix, m_D);
    fdft = fopen(dft_path, "w");
    fprintf(fdft, "omega\tpower\n");
    double omega_2 = do_dft(phi2_hist, t_hist, n_dft, dft_start, 3.0*m_U, 500, fdft);
    fclose(fdft);
    printf("  phi_2 peak omega = %.4f\n", omega_2);

    snprintf(dft_path, sizeof(dft_path), "%s/%s_mD%.2f_phi3_spectrum.tsv", outdir, prefix, m_D);
    fdft = fopen(dft_path, "w");
    fprintf(fdft, "omega\tpower\n");
    double omega_3 = do_dft(phi3_hist, t_hist, n_dft, dft_start, 3.0*m_U, 500, fdft);
    fclose(fdft);
    printf("  phi_3 peak omega = %.4f\n", omega_3);

    snprintf(dft_path, sizeof(dft_path), "%s/%s_mD%.2f_P_spectrum.tsv", outdir, prefix, m_D);
    fdft = fopen(dft_path, "w");
    fprintf(fdft, "omega\tpower\n");
    double omega_P = do_dft(P_hist, t_hist, n_dft, dft_start, 4.0*m_U, 600, fdft);
    fclose(fdft);
    printf("  P(t) peak omega = %.4f\n", omega_P);

    /* Mass matrix eigenvalues (analytical) */
    printf("\n=== Mass Matrix Eigenvalues ===\n");
    double mU2 = m_U*m_U, mD2 = m_D*m_D;
    if (mode == 0) {
        /* UDD: M² = [[mU², lam, lam], [lam, mD², lam], [lam, lam, mD²]] */
        /* DD antisymmetric: eigenvalue = mD² - lambda (exact, doesn't involve U) */
        double eig_DD = mD2 - lambda;
        /* The other two are from 2x2 block: trace = mU²+mD²+2lam, det = ... */
        double tr = mU2 + mD2 + 2*lambda;
        double det = (mU2 + lambda) * (mD2 + lambda) - 2*lambda*lambda;
        double disc = tr*tr - 4*det;
        if (disc < 0) disc = 0;
        double eig_plus  = (tr + sqrt(disc)) / 2.0;
        double eig_minus = (tr - sqrt(disc)) / 2.0;
        printf("  UDD mass matrix eigenvalues:\n");
        printf("    DD antisymmetric: %.4f (omega=%.4f)\n", eig_DD, eig_DD > 0 ? sqrt(eig_DD) : 0.0);
        printf("    Mixed mode -:     %.4f (omega=%.4f)\n", eig_minus, eig_minus > 0 ? sqrt(eig_minus) : 0.0);
        printf("    Mixed mode +:     %.4f (omega=%.4f)\n", eig_plus, eig_plus > 0 ? sqrt(eig_plus) : 0.0);
    } else {
        /* UUD: M² = [[mU², lam, lam], [lam, mU², lam], [lam, lam, mD²]] */
        double eig_UU = mU2 - lambda;
        double tr = mU2 + mD2 + 2*lambda;
        double det = (mU2 + lambda) * (mD2 + lambda) - 2*lambda*lambda;
        double disc = tr*tr - 4*det;
        if (disc < 0) disc = 0;
        double eig_plus  = (tr + sqrt(disc)) / 2.0;
        double eig_minus = (tr - sqrt(disc)) / 2.0;
        printf("  UUD mass matrix eigenvalues:\n");
        printf("    UU antisymmetric: %.4f (omega=%.4f)\n", eig_UU, eig_UU > 0 ? sqrt(eig_UU) : 0.0);
        printf("    Mixed mode -:     %.4f (omega=%.4f)\n", eig_minus, eig_minus > 0 ? sqrt(eig_minus) : 0.0);
        printf("    Mixed mode +:     %.4f (omega=%.4f)\n", eig_plus, eig_plus > 0 ? sqrt(eig_plus) : 0.0);
    }

    /* Summary */
    printf("\n=== Summary ===\n");
    printf("  Mode: %s, m_U=%.3f, m_D=%.3f, lambda=%.3f\n",
           mode_str, m_U, m_D, lambda);
    printf("  E_init=%.4f  E(t=1000)=%.4f  E(t=5000)=%.4f  E_final=%.4f\n",
           E_init, E_t1000, E_t5000, E_final);
    printf("  E retained: %.1f%% (t=1000), %.1f%% (t=5000), %.1f%% (t=end)\n",
           100.0*E_t1000/E_init, 100.0*E_t5000/E_init, 100.0*E_final/E_init);
    printf("  fc: %.3f -> %.3f\n", fc_init, fc_final);
    printf("  peaks: (%.3f,%.3f,%.3f) -> (%.3f,%.3f,%.3f)\n",
           peak_init[0], peak_init[1], peak_init[2],
           peak_final[0], peak_final[1], peak_final[2]);
    printf("  omega_phi: (%.4f, %.4f, %.4f)\n", omega_1, omega_2, omega_3);

    /* Write summary line to stdout in TSV format for easy collection */
    printf("\nSUMMARY_TSV\t%s\t%.3f\t%.3f\t%.3f\t%.4f\t%.4f\t%.4f\t%.4f\t%.1f\t%.1f\t%.3f\t%.4f\t%.4f\t%.4f\n",
           (mode==0)?"UDD":"UUD", m_U, m_D, lambda,
           E_init, E_t1000, E_t5000, E_final,
           100.0*E_t5000/E_init, 100.0*E_final/E_init,
           fc_final, omega_1, omega_2, omega_3);

    printf("\nOutput: %s\n", tspath);

    for (int a = 0; a < 3; a++) { free(phi[a]); free(vel[a]); free(acc_f[a]); }
    free(damp);
    free(phi1_hist); free(phi2_hist); free(phi3_hist);
    free(P_hist); free(t_hist);
    return 0;
}
