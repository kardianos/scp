/*
 * oscillon_1d.c — 1D three-field oscillon with configurable mass coupling
 *
 * Three real scalar fields phi_0, phi_1, phi_2:
 *   d^2 phi_a/dt^2 = d^2 phi_a/dx^2 - m_eff^2 * phi_a - dV/dphi_a
 *
 * V(P) = (mu/2) P^2 / (1 + kappa_eff * P^2),  P = phi_0 * phi_1 * phi_2
 *
 * Mass coupling modes:
 *   0: constant        m_eff^2 = m^2
 *   1: inverse         m_eff^2 = alpha / (1 + beta * Sigma)
 *   2: hybrid          m_eff^2 = m_const^2 + alpha / (1 + beta * Sigma)
 *   3: density kappa   m_eff^2 = m^2, kappa_eff = kappa0 / (1 + gamma * Sigma)
 *
 * where Sigma = phi_0^2 + phi_1^2 + phi_2^2.
 *
 * Build: gcc -O3 -fopenmp -o oscillon_1d src/oscillon_1d.c -lm
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <omp.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ================================================================
   Parameters
   ================================================================ */

/* Physics */
static double mu      = -41.345;
static double kappa   = 50.0;
static double mass2   = 2.25;      /* m^2 for modes 0,2,3 */

/* Mode 1/2: inverse coupling */
static double inv_alpha = 2.25;
static double inv_beta  = 5.0;

/* Mode 2: hybrid constant part */
static double m_const2  = 1.0;

/* Mode 3: density-dependent kappa */
static double kappa_gamma = 2.0;

/* Coupling mode */
static int mode = 0;

/* Grid */
static int    N  = 2048;
static double L  = 50.0;

/* Time */
static double T_final  = 1000.0;
static double diag_dt  = 10.0;
static double snap_dt  = 50.0;

/* Initial condition */
static double A_init  = 0.8;
static double sigma_w = 3.0;
static double k_wave  = 0.0;
static double blob_W  = 5.0;
static double delta[3] = {0.0, 3.0005, 4.4325};

enum InitType { INIT_OSCILLON=0, INIT_BLOB, INIT_BRAID };
static int init_type = INIT_OSCILLON;

/* Output */
static char outdir[512] = "output";

/* ================================================================
   Command-line parsing
   ================================================================ */

static void usage(const char *prog)
{
    fprintf(stderr,
        "Usage: %s [options]\n"
        "  -mode N         Mass coupling mode (0=const, 1=inverse, 2=hybrid, 3=density-kappa)\n"
        "  -m F            Mass (sqrt of m^2) for modes 0,2,3 [%.3f]\n"
        "  -mu F           Potential mu [%.3f]\n"
        "  -kappa F        Potential kappa [%.1f]\n"
        "  -inv_alpha F    Inverse coupling alpha [%.3f]\n"
        "  -inv_beta F     Inverse coupling beta [%.3f]\n"
        "  -m_const F      Hybrid constant mass [%.3f]\n"
        "  -kappa_gamma F  Density kappa gamma [%.3f]\n"
        "  -N I            Grid points [%d]\n"
        "  -L F            Half-domain size [%.1f]\n"
        "  -T F            Final time [%.1f]\n"
        "  -A F            Initial amplitude [%.3f]\n"
        "  -sigma F        Gaussian width [%.3f]\n"
        "  -k F            Wave number [%.3f]\n"
        "  -blob_W F       Blob half-width [%.3f]\n"
        "  -diag_dt F      Diagnostic interval [%.1f]\n"
        "  -snap_dt F      Snapshot interval [%.1f]\n"
        "  -init S         Initial condition (oscillon, blob, braid)\n"
        "  -o DIR          Output directory [%s]\n",
        prog, sqrt(mass2), mu, kappa, inv_alpha, inv_beta,
        sqrt(m_const2), kappa_gamma, N, L, T_final,
        A_init, sigma_w, k_wave, blob_W, diag_dt, snap_dt, outdir);
    exit(1);
}

static void parse_args(int argc, char **argv)
{
    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "--help")) usage(argv[0]);
        if (i + 1 >= argc) { fprintf(stderr, "Missing argument for %s\n", argv[i]); exit(1); }

        if      (!strcmp(argv[i], "-mode"))        mode = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-m"))          { double m = atof(argv[++i]); mass2 = m*m; }
        else if (!strcmp(argv[i], "-mu"))           mu = atof(argv[++i]);
        else if (!strcmp(argv[i], "-kappa"))        kappa = atof(argv[++i]);
        else if (!strcmp(argv[i], "-inv_alpha"))    inv_alpha = atof(argv[++i]);
        else if (!strcmp(argv[i], "-inv_beta"))     inv_beta = atof(argv[++i]);
        else if (!strcmp(argv[i], "-m_const"))    { double m = atof(argv[++i]); m_const2 = m*m; }
        else if (!strcmp(argv[i], "-kappa_gamma"))  kappa_gamma = atof(argv[++i]);
        else if (!strcmp(argv[i], "-N"))            N = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-L"))            L = atof(argv[++i]);
        else if (!strcmp(argv[i], "-T"))            T_final = atof(argv[++i]);
        else if (!strcmp(argv[i], "-A"))            A_init = atof(argv[++i]);
        else if (!strcmp(argv[i], "-sigma"))        sigma_w = atof(argv[++i]);
        else if (!strcmp(argv[i], "-k"))            k_wave = atof(argv[++i]);
        else if (!strcmp(argv[i], "-blob_W"))       blob_W = atof(argv[++i]);
        else if (!strcmp(argv[i], "-diag_dt"))      diag_dt = atof(argv[++i]);
        else if (!strcmp(argv[i], "-snap_dt"))      snap_dt = atof(argv[++i]);
        else if (!strcmp(argv[i], "-o"))            strncpy(outdir, argv[++i], sizeof(outdir)-1);
        else if (!strcmp(argv[i], "-init")) {
            i++;
            if      (!strcmp(argv[i], "oscillon")) init_type = INIT_OSCILLON;
            else if (!strcmp(argv[i], "blob"))      init_type = INIT_BLOB;
            else if (!strcmp(argv[i], "braid"))     init_type = INIT_BRAID;
            else { fprintf(stderr, "Unknown init: %s\n", argv[i]); exit(1); }
        }
        else { fprintf(stderr, "Unknown flag: %s\n", argv[i]); exit(1); }
    }

    if (mode < 0 || mode > 3) {
        fprintf(stderr, "Invalid mode %d (must be 0-3)\n", mode);
        exit(1);
    }
}

/* ================================================================
   Effective mass and kappa
   ================================================================ */

static inline double meff2_at(double sigma_phi2)
{
    switch (mode) {
    case 0: return mass2;
    case 1: return inv_alpha / (1.0 + inv_beta * sigma_phi2);
    case 2: return m_const2 + inv_alpha / (1.0 + inv_beta * sigma_phi2);
    case 3: return mass2;
    default: return mass2;
    }
}

static inline double kappa_eff_at(double sigma_phi2)
{
    if (mode == 3)
        return kappa / (1.0 + kappa_gamma * sigma_phi2);
    return kappa;
}

/* ================================================================
   Potential force: -dV/dphi_a
   V(P) = (mu/2) P^2 / (1 + kappa_eff * P^2)
   dV/dP = mu * P / (1 + kappa_eff * P^2)^2
   dV/dphi_a = dV/dP * dP/dphi_a
   ================================================================ */

static inline double force_pot(double p0, double p1, double p2, int a, double keff)
{
    double P = p0 * p1 * p2;
    double P2 = P * P;
    double denom = 1.0 + keff * P2;
    double dVdP = mu * P / (denom * denom);
    double dPda;
    switch (a) {
    case 0: dPda = p1 * p2; break;
    case 1: dPda = p0 * p2; break;
    case 2: dPda = p0 * p1; break;
    default: dPda = 0.0;
    }
    return -dVdP * dPda;
}

/* ================================================================
   Main
   ================================================================ */

int main(int argc, char **argv)
{
    parse_args(argc, argv);

    /* Grid setup */
    double dx  = 2.0 * L / N;   /* periodic: N cells over [-L, L) */
    double dx2 = dx * dx;

    /* CFL: dt < dx / c, with some safety for the mass term */
    double dt_cfl = 0.4 * dx;
    /* Also limit by mass oscillation: dt < 2/m */
    double m_max = sqrt(mass2);
    if (mode == 1) m_max = sqrt(inv_alpha);
    if (mode == 2) m_max = sqrt(m_const2 + inv_alpha);
    if (m_max > 0) {
        double dt_mass = 0.4 * 2.0 / m_max;
        if (dt_mass < dt_cfl) dt_cfl = dt_mass;
    }
    double dt = dt_cfl;
    long Nt = (long)(T_final / dt) + 1;

    /* Diagnostic intervals in steps */
    long diag_steps = (long)(diag_dt / dt);
    if (diag_steps < 1) diag_steps = 1;
    long snap_steps = (long)(snap_dt / dt);
    if (snap_steps < 1) snap_steps = 1;

    const char *mode_names[] = {"constant", "inverse", "hybrid", "density-kappa"};
    const char *init_names[] = {"oscillon", "blob", "braid"};

    printf("oscillon_1d: mode=%d (%s), init=%s\n", mode, mode_names[mode], init_names[init_type]);
    printf("  mu=%.3f kappa=%.1f\n", mu, kappa);
    switch (mode) {
    case 0: printf("  m^2=%.4f (m=%.4f)\n", mass2, sqrt(mass2)); break;
    case 1: printf("  inv_alpha=%.4f inv_beta=%.4f\n", inv_alpha, inv_beta); break;
    case 2: printf("  m_const^2=%.4f inv_alpha=%.4f inv_beta=%.4f\n",
                   m_const2, inv_alpha, inv_beta); break;
    case 3: printf("  m^2=%.4f kappa_gamma=%.4f\n", mass2, kappa_gamma); break;
    }
    printf("  N=%d L=%.1f dx=%.6f dt=%.6f\n", N, L, dx, dt);
    printf("  T=%.1f (%ld steps), diag every %ld steps, snap every %ld steps\n",
           T_final, Nt, diag_steps, snap_steps);
    printf("  A=%.3f sigma=%.3f k=%.3f\n", A_init, sigma_w, k_wave);
    printf("  delta = {%.4f, %.4f, %.4f}\n", delta[0], delta[1], delta[2]);

    /* Create output directory */
    mkdir(outdir, 0755);

    /* Allocate fields: phi, vel, acc for 3 components */
    double *phi[3], *vel[3], *acc[3];
    for (int a = 0; a < 3; a++) {
        phi[a] = calloc(N, sizeof(double));
        vel[a] = calloc(N, sizeof(double));
        acc[a] = calloc(N, sizeof(double));
        if (!phi[a] || !vel[a] || !acc[a]) {
            fprintf(stderr, "FATAL: allocation failed\n");
            return 1;
        }
    }

    /* --------------------------------------------------------
       Initialize fields
       -------------------------------------------------------- */
    switch (init_type) {
    case INIT_OSCILLON:
        for (int a = 0; a < 3; a++)
            for (int i = 0; i < N; i++) {
                double x = -L + (i + 0.5) * dx;
                phi[a][i] = A_init * exp(-x*x / (2.0 * sigma_w * sigma_w))
                            * cos(k_wave * x + delta[a]);
            }
        /* Zero initial velocity */
        break;

    case INIT_BLOB:
        for (int a = 0; a < 3; a++)
            for (int i = 0; i < N; i++) {
                double x = -L + (i + 0.5) * dx;
                phi[a][i] = (fabs(x) < blob_W) ? A_init : 0.0;
            }
        break;

    case INIT_BRAID:
        /* Traveling wave: phi_a with initial velocity from omega * A * sin(...) */
        {
            /* Estimate omega from dispersion: omega^2 = k^2 + m_eff^2 */
            double meff2_0 = meff2_at(3.0 * A_init * A_init);
            double omega = sqrt(k_wave * k_wave + meff2_0);
            if (omega < 0.01) omega = 1.0;

            printf("  braid: omega=%.4f (k=%.3f, meff2=%.4f)\n", omega, k_wave, meff2_0);

            for (int a = 0; a < 3; a++)
                for (int i = 0; i < N; i++) {
                    double x = -L + (i + 0.5) * dx;
                    double env = exp(-x*x / (2.0 * sigma_w * sigma_w));
                    double phase = k_wave * x + delta[a];
                    phi[a][i] = A_init * env * cos(phase);
                    vel[a][i] = -A_init * env * omega * sin(phase);
                }
        }
        break;
    }

    /* --------------------------------------------------------
       Compute initial acceleration
       -------------------------------------------------------- */
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < N; i++) {
        int ip = (i + 1) % N;
        int im = (i - 1 + N) % N;
        double p0 = phi[0][i], p1 = phi[1][i], p2 = phi[2][i];
        double sig = p0*p0 + p1*p1 + p2*p2;
        double me2 = meff2_at(sig);
        double keff = kappa_eff_at(sig);
        for (int a = 0; a < 3; a++) {
            double lapl = (phi[a][ip] - 2.0*phi[a][i] + phi[a][im]) / dx2;
            double fp = force_pot(p0, p1, p2, a, keff);
            acc[a][i] = lapl - me2 * phi[a][i] + fp;
        }
    }

    /* --------------------------------------------------------
       Open timeseries output
       -------------------------------------------------------- */
    char tspath[600];
    snprintf(tspath, sizeof(tspath), "%s/timeseries.tsv", outdir);
    FILE *fts = fopen(tspath, "w");
    if (!fts) { fprintf(stderr, "Cannot open %s\n", tspath); return 1; }
    fprintf(fts, "t\tE_kin\tE_grad\tE_mass\tE_pot\tE_total\tphi_max\tP_max\tm_eff_center\n");

    int snap_count = 0;

    /* --------------------------------------------------------
       Time evolution: Velocity Verlet
       -------------------------------------------------------- */
    for (long n = 0; n <= Nt; n++) {
        double t = n * dt;

        /* ---- Diagnostics ---- */
        int do_diag = (n % diag_steps == 0);
        int do_snap = (n % snap_steps == 0);

        if (do_diag || do_snap) {
            double E_kin = 0, E_grad = 0, E_mass = 0, E_pot = 0;
            double phi_max = 0, P_max = 0;

            #pragma omp parallel for reduction(+:E_kin,E_grad,E_mass,E_pot) \
                                     reduction(max:phi_max,P_max) schedule(static)
            for (int i = 0; i < N; i++) {
                int ip = (i + 1) % N;
                double p0 = phi[0][i], p1 = phi[1][i], p2 = phi[2][i];
                double sig = p0*p0 + p1*p1 + p2*p2;
                double me2 = meff2_at(sig);
                double keff = kappa_eff_at(sig);

                double ek = 0, eg = 0, em = 0;
                for (int a = 0; a < 3; a++) {
                    double v = vel[a][i];
                    ek += 0.5 * v * v;
                    double dp = (phi[a][ip] - phi[a][i]) / dx;
                    eg += 0.5 * dp * dp;
                    em += 0.5 * me2 * phi[a][i] * phi[a][i];

                    double ap = fabs(phi[a][i]);
                    if (ap > phi_max) phi_max = ap;
                }

                double P = p0 * p1 * p2;
                double P2 = P * P;
                double ep = 0.5 * mu * P2 / (1.0 + keff * P2);
                E_pot += ep;

                double Pabs = fabs(P);
                if (Pabs > P_max) P_max = Pabs;

                E_kin  += ek;
                E_grad += eg;
                E_mass += em;
            }
            E_kin  *= dx;
            E_grad *= dx;
            E_mass *= dx;
            E_pot  *= dx;
            double E_total = E_kin + E_grad + E_mass + E_pot;

            /* m_eff at center */
            int ic = N / 2;
            double sig_c = phi[0][ic]*phi[0][ic] + phi[1][ic]*phi[1][ic] + phi[2][ic]*phi[2][ic];
            double me_center = sqrt(fabs(meff2_at(sig_c)));

            if (do_diag) {
                fprintf(fts, "%.4f\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\n",
                        t, E_kin, E_grad, E_mass, E_pot, E_total,
                        phi_max, P_max, me_center);
            }

            /* Print summary */
            if (do_diag) {
                printf("  t=%8.1f  E=%.4f (K=%.3f G=%.3f M=%.3f V=%.3f)  "
                       "phi_max=%.4f  P_max=%.4e  m_c=%.4f\n",
                       t, E_total, E_kin, E_grad, E_mass, E_pot,
                       phi_max, P_max, me_center);
            }

            /* Field snapshot */
            if (do_snap) {
                char snappath[600];
                snprintf(snappath, sizeof(snappath), "%s/field_t%04d.dat",
                         outdir, snap_count);
                FILE *fsnap = fopen(snappath, "w");
                if (fsnap) {
                    fprintf(fsnap, "# t=%.4f  mode=%d\n", t, mode);
                    fprintf(fsnap, "# x\tphi_0\tphi_1\tphi_2\n");
                    /* Write every point (or subsample for very large N) */
                    int stride = 1;
                    if (N > 8192) stride = N / 4096;
                    for (int i = 0; i < N; i += stride) {
                        double x = -L + (i + 0.5) * dx;
                        fprintf(fsnap, "%.6f\t%.8e\t%.8e\t%.8e\n",
                                x, phi[0][i], phi[1][i], phi[2][i]);
                    }
                    fclose(fsnap);
                }
                snap_count++;
            }
        }

        if (n == Nt) break;

        /* ---- Verlet: half-kick ---- */
        #pragma omp parallel for schedule(static)
        for (int i = 0; i < N; i++)
            for (int a = 0; a < 3; a++)
                vel[a][i] += 0.5 * dt * acc[a][i];

        /* ---- Verlet: drift ---- */
        #pragma omp parallel for schedule(static)
        for (int i = 0; i < N; i++)
            for (int a = 0; a < 3; a++)
                phi[a][i] += dt * vel[a][i];

        /* ---- Verlet: recompute acceleration ---- */
        #pragma omp parallel for schedule(static)
        for (int i = 0; i < N; i++) {
            int ip = (i + 1) % N;
            int im = (i - 1 + N) % N;
            double p0 = phi[0][i], p1 = phi[1][i], p2 = phi[2][i];
            double sig = p0*p0 + p1*p1 + p2*p2;
            double me2 = meff2_at(sig);
            double keff = kappa_eff_at(sig);
            for (int a = 0; a < 3; a++) {
                double lapl = (phi[a][ip] - 2.0*phi[a][i] + phi[a][im]) / dx2;
                double fp = force_pot(p0, p1, p2, a, keff);
                acc[a][i] = lapl - me2 * phi[a][i] + fp;
            }
        }

        /* ---- Verlet: half-kick ---- */
        #pragma omp parallel for schedule(static)
        for (int i = 0; i < N; i++)
            for (int a = 0; a < 3; a++)
                vel[a][i] += 0.5 * dt * acc[a][i];
    }

    fclose(fts);

    /* --------------------------------------------------------
       Final summary
       -------------------------------------------------------- */
    /* Recompute final energies */
    double E_kin = 0, E_grad = 0, E_mass = 0, E_pot = 0;
    double phi_max = 0;

    for (int i = 0; i < N; i++) {
        int ip = (i + 1) % N;
        double p0 = phi[0][i], p1 = phi[1][i], p2 = phi[2][i];
        double sig = p0*p0 + p1*p1 + p2*p2;
        double me2 = meff2_at(sig);
        double keff = kappa_eff_at(sig);
        for (int a = 0; a < 3; a++) {
            E_kin  += 0.5 * vel[a][i] * vel[a][i];
            double dp = (phi[a][ip] - phi[a][i]) / dx;
            E_grad += 0.5 * dp * dp;
            E_mass += 0.5 * me2 * phi[a][i] * phi[a][i];
            if (fabs(phi[a][i]) > phi_max) phi_max = fabs(phi[a][i]);
        }
        double P = p0 * p1 * p2;
        double P2 = P * P;
        E_pot += 0.5 * mu * P2 / (1.0 + keff * P2);
    }
    E_kin *= dx; E_grad *= dx; E_mass *= dx; E_pot *= dx;

    printf("\n=== FINAL (t=%.1f) ===\n", T_final);
    printf("  E_total = %.6f\n", E_kin + E_grad + E_mass + E_pot);
    printf("  E_kin=%.4f  E_grad=%.4f  E_mass=%.4f  E_pot=%.4f\n",
           E_kin, E_grad, E_mass, E_pot);
    printf("  phi_max = %.6f\n", phi_max);
    printf("  Snapshots: %d files in %s/\n", snap_count, outdir);
    printf("  Timeseries: %s\n", tspath);

    /* Cleanup */
    for (int a = 0; a < 3; a++) {
        free(phi[a]); free(vel[a]); free(acc[a]);
    }

    return 0;
}
