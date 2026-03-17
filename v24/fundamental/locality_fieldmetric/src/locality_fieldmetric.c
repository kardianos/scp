/*
 * locality_fieldmetric.c — Field-derived effective metric from energy density
 *
 * Option 5: g_eff(x,t) = 1 + beta * <rho(x,t)> / rho_0
 * where <rho> is EMA-smoothed energy density, rho_0 = average energy density.
 *
 * Backreaction: m^2_eff = m^2 * g_eff, c^2_eff = g_eff (wave speed modification)
 *
 * Three massive scalars with triple-product coupling (v21 oscillon model):
 *   L = sum_a [ (1/2)(dt phi_a)^2 - (c^2_eff/2)(dx phi_a)^2 - (m^2_eff/2)phi_a^2 ]
 *     - (mu/2) P^2 / (1 + kappa P^2),  P = phi_1 phi_2 phi_3
 *
 * Phase 1: Single oscillon, scan beta
 * Phase 2: 8-oscillon lattice with g_eff at best beta
 * Phase 3: Causality test — boost one oscillon, measure g_eff propagation delay
 *
 * Compile: gcc -O3 -Wall -o locality_fieldmetric src/locality_fieldmetric.c -lm
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Physics parameters */
static double mu      = -20.0;
static double kap     = 20.0;
static double mass    = 1.0;
static double A_init  = 0.8;
static double sig     = 3.0;

/* Grid */
static int    Nx      = 4000;
static double xmax    = 100.0;

/* Time */
static double t_equil = 5000.0;
static double t_test  = 5000.0;

/* Metric parameters */
static double beta_g  = -0.01;
static double tau_ema = 7.0;

/* Lattice */
static int    N_osc   = 1;
static double lam     = 0.5;     /* lattice: amplitude scaling */
static double d_sep   = 16.0;    /* oscillon separation */

/* Boost (Phase 3) */
static double v_boost = 0.0;     /* boost velocity for rightmost oscillon */
static int    boost_idx = -1;    /* which oscillon to boost (-1 = rightmost) */

/* Output */
static char   outdir[512] = "data";

static void parse_args(int argc, char **argv)
{
    for (int i = 1; i < argc - 1; i += 2) {
        if      (!strcmp(argv[i], "-mu"))      mu      = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-kappa"))   kap     = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-mass"))    mass    = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-A"))       A_init  = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-sigma"))   sig     = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-Nx"))      Nx      = atoi(argv[i+1]);
        else if (!strcmp(argv[i], "-xmax"))    xmax    = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-t_equil")) t_equil = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-t_test"))  t_test  = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-beta"))    beta_g  = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-tau"))     tau_ema = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-N_osc"))   N_osc   = atoi(argv[i+1]);
        else if (!strcmp(argv[i], "-lambda"))  lam     = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-d"))       d_sep   = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-vboost"))  v_boost = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-bidx"))    boost_idx = atoi(argv[i+1]);
        else if (!strcmp(argv[i], "-o"))       strncpy(outdir, argv[i+1], sizeof(outdir)-1);
        else { fprintf(stderr, "Unknown flag: %s\n", argv[i]); exit(1); }
    }
}

/* -dV/dphi_a where V = (mu/2)P^2/(1+kappa P^2), P = phi1 phi2 phi3 */
static inline double force_pot(double p1, double p2, double p3, int a)
{
    double P  = p1 * p2 * p3;
    double P2 = P * P;
    double denom2 = (1.0 + kap * P2) * (1.0 + kap * P2);
    double dP;
    switch (a) {
        case 0: dP = p2 * p3; break;
        case 1: dP = p1 * p3; break;
        case 2: dP = p1 * p2; break;
        default: dP = 0.0;
    }
    return -mu * P * dP / denom2;
}

/* Energy density at grid point i (needs phi, vel arrays and dx) */
static inline double energy_density(double **phi, double **vel, int i, double dx, double m2)
{
    double e = 0.0;
    for (int a = 0; a < 3; a++) {
        e += 0.5 * vel[a][i] * vel[a][i];
        double dp = (phi[a][i+1] - phi[a][i-1]) / (2.0 * dx);
        e += 0.5 * dp * dp;
        e += 0.5 * m2 * phi[a][i] * phi[a][i];
    }
    double P  = phi[0][i] * phi[1][i] * phi[2][i];
    double P2 = P * P;
    e += 0.5 * mu * P2 / (1.0 + kap * P2);
    return e;
}


int main(int argc, char **argv)
{
    parse_args(argc, argv);

    double dx  = 2.0 * xmax / (Nx - 1);
    double dx2 = dx * dx;
    double m2  = mass * mass;

    /* CFL */
    double kmax = M_PI / dx;
    double dt = 0.8 * 2.0 / sqrt(kmax * kmax + m2);
    /* With metric correction, effective c can increase; reduce dt for safety */
    double max_geff = 1.0 + fabs(beta_g) * 20.0; /* conservative bound */
    if (max_geff > 1.0)
        dt /= sqrt(max_geff);

    int Nt_equil = (int)(t_equil / dt) + 1;
    int Nt_test  = (int)(t_test / dt) + 1;
    int Nt_total = Nt_equil + Nt_test;

    printf("locality_fieldmetric: mu=%.1f kappa=%.1f mass=%.3f\n", mu, kap, mass);
    printf("  beta=%.6f tau=%.1f N_osc=%d lambda=%.3f d=%.1f\n",
           beta_g, tau_ema, N_osc, lam, d_sep);
    printf("  Nx=%d xmax=%.1f dx=%.5f dt=%.6f\n", Nx, xmax, dx, dt);
    printf("  t_equil=%.0f (%d steps), t_test=%.0f (%d steps)\n",
           t_equil, Nt_equil, t_test, Nt_test);
    if (v_boost != 0.0)
        printf("  v_boost=%.4f on oscillon %d\n", v_boost,
               boost_idx >= 0 ? boost_idx : N_osc - 1);

    /* Allocate fields */
    double *phi[3], *vel[3], *acc[3];
    for (int a = 0; a < 3; a++) {
        phi[a] = calloc(Nx, sizeof(double));
        vel[a] = calloc(Nx, sizeof(double));
        acc[a] = calloc(Nx, sizeof(double));
    }

    /* EMA energy density and g_eff */
    double *rho_ema = calloc(Nx, sizeof(double));
    double *g_eff   = malloc(Nx * sizeof(double));
    for (int i = 0; i < Nx; i++) g_eff[i] = 1.0;

    /* Absorbing boundary: outer 20% */
    double *damp = malloc(Nx * sizeof(double));
    double x_abs = xmax * 0.80;
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

    /* Initialize oscillons */
    double *osc_pos = malloc(N_osc * sizeof(double));
    if (N_osc == 1) {
        osc_pos[0] = 0.0;
    } else {
        /* Center the chain */
        double chain_len = (N_osc - 1) * d_sep;
        for (int k = 0; k < N_osc; k++)
            osc_pos[k] = -chain_len / 2.0 + k * d_sep;
    }

    for (int k = 0; k < N_osc; k++) {
        double amp = (N_osc > 1) ? lam : A_init;
        for (int a = 0; a < 3; a++)
            for (int i = 0; i < Nx; i++) {
                double x = -xmax + i * dx - osc_pos[k];
                phi[a][i] += amp * exp(-x * x / (2.0 * sig * sig));
            }
    }

    /* Apply boost if requested (Phase 3) */
    int boosted = (boost_idx >= 0) ? boost_idx : N_osc - 1;
    /* Boost is applied AFTER equilibration, but we need to track which */

    /* Compute acceleration (with g_eff) */
    #define COMPUTE_ACC_GEFF() do { \
        for (int a = 0; a < 3; a++) { \
            acc[a][0] = acc[a][1] = acc[a][Nx-2] = acc[a][Nx-1] = 0; \
            for (int i = 1; i < Nx - 1; i++) { \
                double ge = g_eff[i]; \
                double c2e = ge; \
                double m2e = m2 * ge; \
                double lapl = (phi[a][i+1] - 2.0*phi[a][i] + phi[a][i-1]) / dx2; \
                double fp = force_pot(phi[0][i], phi[1][i], phi[2][i], a); \
                acc[a][i] = c2e * lapl - m2e * phi[a][i] + fp; \
            } \
        } \
    } while(0)

    /* No metric during equilibration */
    #define COMPUTE_ACC_FREE() do { \
        for (int a = 0; a < 3; a++) { \
            acc[a][0] = acc[a][1] = acc[a][Nx-2] = acc[a][Nx-1] = 0; \
            for (int i = 1; i < Nx - 1; i++) { \
                double lapl = (phi[a][i+1] - 2.0*phi[a][i] + phi[a][i-1]) / dx2; \
                double fp = force_pot(phi[0][i], phi[1][i], phi[2][i], a); \
                acc[a][i] = lapl - m2 * phi[a][i] + fp; \
            } \
        } \
    } while(0)

    COMPUTE_ACC_FREE();

    /* Diagnostics storage */
    int max_dft = 50000;
    double *phi0_hist = malloc(max_dft * sizeof(double));
    double *t_hist    = malloc(max_dft * sizeof(double));
    int n_dft = 0;

    int rec_every   = Nt_total / 20000; if (rec_every < 1) rec_every = 1;
    int print_every = Nt_total / 60;    if (print_every < 1) print_every = 1;
    int dft_every   = Nt_total / max_dft; if (dft_every < 1) dft_every = 1;

    double core_r = 3.0 * sig;

    /* Output files */
    char tspath[600], gpath[600], causalpath[600];
    snprintf(tspath, sizeof(tspath), "%s/fieldmetric_ts.tsv", outdir);
    FILE *fts = fopen(tspath, "w");
    if (!fts) { fprintf(stderr, "Cannot open %s\n", tspath); return 1; }
    fprintf(fts, "time\tphase\tphi1_0\tpeak1\tE_total\tE_pot\tf_core\tg_eff_0\trho_0\n");

    /* For lattice: track g_eff at each oscillon position */
    FILE *fgeff = NULL;
    if (N_osc > 1) {
        snprintf(gpath, sizeof(gpath), "%s/fieldmetric_geff.tsv", outdir);
        fgeff = fopen(gpath, "w");
        fprintf(fgeff, "time\tphase");
        for (int k = 0; k < N_osc; k++) fprintf(fgeff, "\tg_%d\trho_%d", k, k);
        fprintf(fgeff, "\n");
    }

    /* For causality: track g_eff at all oscillon positions with fine time resolution */
    FILE *fcausal = NULL;
    if (v_boost != 0.0) {
        snprintf(causalpath, sizeof(causalpath), "%s/fieldmetric_causal.tsv", outdir);
        fcausal = fopen(causalpath, "w");
        fprintf(fcausal, "time");
        for (int k = 0; k < N_osc; k++) fprintf(fcausal, "\tg_%d\tx_%d", k, k);
        fprintf(fcausal, "\n");
    }

    /* Index of center (for single oscillon) or first oscillon pos */
    int ic = Nx / 2;

    /* rho_0 will be computed from equilibrated state */
    double rho_0 = 1.0;  /* placeholder, updated at end of equilibration */

    /* EMA alpha */
    double ema_alpha = dt / tau_ema;  /* exponential moving average factor */

    /* Track equilibration energy for rho_0 computation */

    printf("\n=== Phase: Equilibration (t=0 to %.0f, no metric) ===\n", t_equil);

    int phase = 0; /* 0=equil, 1=metric-on, 2=boost */

    for (int n = 0; n <= Nt_total; n++) {
        double t = n * dt;

        /* Phase transitions */
        if (n == Nt_equil && phase == 0) {
            phase = 1;
            printf("\n=== Phase: Metric ON (t=%.0f, beta=%.6f) ===\n", t, beta_g);

            /* Compute rho_0 = total energy / domain length */
            double E_tot = 0.0;
            for (int i = 1; i < Nx - 1; i++) {
                E_tot += energy_density(phi, vel, i, dx, m2) * dx;
            }
            /* Use physical domain (exclude absorbing region) */
            rho_0 = E_tot / (2.0 * x_abs);
            printf("  rho_0 = E_total / L_phys = %.6f / %.1f = %.6e\n",
                   E_tot, 2.0*x_abs, rho_0);
            printf("  At center: rho_center ~ %.4f, g_eff-1 ~ %.6f\n",
                   energy_density(phi, vel, ic, dx, m2),
                   beta_g * energy_density(phi, vel, ic, dx, m2) / rho_0);

            /* Initialize EMA to current energy density */
            for (int i = 1; i < Nx - 1; i++)
                rho_ema[i] = energy_density(phi, vel, i, dx, m2);

            /* Apply boost if in Phase 3 mode */
            if (v_boost != 0.0) {
                printf("  Applying boost v=%.4f to oscillon %d at x=%.1f\n",
                       v_boost, boosted, osc_pos[boosted]);
                /* Lorentz boost: phi(x) -> phi(gamma*(x - x0)), vel += -v*dphi/dx */
                int i_osc = (int)((osc_pos[boosted] + xmax) / dx);
                /* Simple: add velocity kick within oscillon region */
                for (int a = 0; a < 3; a++) {
                    for (int i = 2; i < Nx - 2; i++) {
                        double x = -xmax + i * dx;
                        double r = fabs(x - osc_pos[boosted]);
                        if (r < 4.0 * sig) {
                            double dphi = (phi[a][i+1] - phi[a][i-1]) / (2.0 * dx);
                            vel[a][i] += -v_boost * dphi;
                        }
                    }
                }
                (void)i_osc;
            }
        }

        /* Update EMA and g_eff if metric is on */
        if (phase >= 1 && n > Nt_equil) {
            for (int i = 1; i < Nx - 1; i++) {
                double rho_now = energy_density(phi, vel, i, dx, m2);
                rho_ema[i] += ema_alpha * (rho_now - rho_ema[i]);
                g_eff[i] = 1.0 + beta_g * rho_ema[i] / rho_0;
            }
        }

        /* DFT recording */
        if (n % dft_every == 0 && n_dft < max_dft) {
            phi0_hist[n_dft] = phi[0][ic];
            t_hist[n_dft] = t;
            n_dft++;
        }

        /* Diagnostics */
        int do_rec   = (n % rec_every == 0);
        int do_print = (n % print_every == 0);

        if (do_rec || do_print) {
            double Ek = 0, Eg = 0, Em = 0, Ep = 0;
            double Ecore = 0, Eall = 0;
            double peak1 = 0;

            for (int i = 1; i < Nx - 1; i++) {
                double x = -xmax + i * dx;
                for (int a = 0; a < 3; a++) {
                    Ek += 0.5 * vel[a][i] * vel[a][i] * dx;
                    double dp = (phi[a][i+1] - phi[a][i-1]) / (2.0*dx);
                    Eg += 0.5 * dp * dp * dx;
                    Em += 0.5 * m2 * phi[a][i] * phi[a][i] * dx;
                    if (a == 0 && fabs(phi[a][i]) > peak1) peak1 = fabs(phi[a][i]);
                }
                double P = phi[0][i] * phi[1][i] * phi[2][i];
                double P2 = P * P;
                Ep += 0.5 * mu * P2 / (1.0 + kap * P2) * dx;

                double e_loc = energy_density(phi, vel, i, dx, m2);
                Eall += e_loc * dx;
                if (N_osc == 1 && fabs(x) < core_r) Ecore += e_loc * dx;
                if (N_osc > 1 && fabs(x - osc_pos[0]) < core_r) Ecore += e_loc * dx;
            }
            double Et = Ek + Eg + Em + Ep;
            double fc = (Eall > 1e-20) ? Ecore / Eall : 0.0;
            double geff0 = g_eff[ic];
            double rho0  = rho_ema[ic];

            if (do_rec) {
                fprintf(fts, "%.4f\t%d\t%.6e\t%.6e\t%.6e\t%.6e\t%.6f\t%.8f\t%.6e\n",
                        t, phase, phi[0][ic], peak1, Et, Ep, fc, geff0, rho0);
            }
            if (do_print) {
                printf("  t=%7.1f [%d]  phi0=%+.4f  pk=%.4f  E=%.4f  Ep=%+.4f  "
                       "fc=%.3f  g=%.6f\n",
                       t, phase, phi[0][ic], peak1, Et, Ep, fc, geff0);
            }
        }

        /* Causality tracking (fine time resolution) */
        if (fcausal && phase >= 1 && n % (rec_every > 5 ? rec_every/5 : 1) == 0) {
            fprintf(fcausal, "%.4f", t);
            for (int k = 0; k < N_osc; k++) {
                int idx = (int)((osc_pos[k] + xmax) / dx);
                if (idx < 1) idx = 1;
                if (idx >= Nx-1) idx = Nx-2;
                fprintf(fcausal, "\t%.8f\t%.2f", g_eff[idx], osc_pos[k]);
            }
            fprintf(fcausal, "\n");
        }

        /* Lattice g_eff tracking */
        if (fgeff && n % rec_every == 0) {
            fprintf(fgeff, "%.4f\t%d", t, phase);
            for (int k = 0; k < N_osc; k++) {
                int idx = (int)((osc_pos[k] + xmax) / dx);
                if (idx < 1) idx = 1;
                if (idx >= Nx-1) idx = Nx-2;
                fprintf(fgeff, "\t%.8f\t%.6e", g_eff[idx], rho_ema[idx]);
            }
            fprintf(fgeff, "\n");
        }

        if (n == Nt_total) break;

        /* Velocity Verlet integration */
        if (phase == 0) {
            /* Free evolution during equilibration */
            for (int a = 0; a < 3; a++)
                for (int i = 1; i < Nx - 1; i++)
                    vel[a][i] += 0.5 * dt * acc[a][i];
            for (int a = 0; a < 3; a++)
                for (int i = 1; i < Nx - 1; i++)
                    phi[a][i] += dt * vel[a][i];
            COMPUTE_ACC_FREE();
            for (int a = 0; a < 3; a++)
                for (int i = 1; i < Nx - 1; i++)
                    vel[a][i] += 0.5 * dt * acc[a][i];
        } else {
            /* Evolution with g_eff metric */
            for (int a = 0; a < 3; a++)
                for (int i = 1; i < Nx - 1; i++)
                    vel[a][i] += 0.5 * dt * acc[a][i];
            for (int a = 0; a < 3; a++)
                for (int i = 1; i < Nx - 1; i++)
                    phi[a][i] += dt * vel[a][i];
            COMPUTE_ACC_GEFF();
            for (int a = 0; a < 3; a++)
                for (int i = 1; i < Nx - 1; i++)
                    vel[a][i] += 0.5 * dt * acc[a][i];
        }

        /* Absorbing boundary */
        for (int a = 0; a < 3; a++)
            for (int i = 0; i < Nx; i++) {
                vel[a][i] *= damp[i];
                phi[a][i] *= damp[i];
            }
    }

    fclose(fts);
    if (fgeff) fclose(fgeff);
    if (fcausal) fclose(fcausal);

    /* DFT: separate equilibrium and test phases */
    printf("\n=== Frequency Analysis ===\n");
    {
        /* Find equilibrium and test DFT ranges */
        int eq_start = -1, eq_end = -1, ts_start = -1, ts_end = -1;
        for (int j = 0; j < n_dft; j++) {
            if (t_hist[j] >= t_equil * 0.5 && eq_start < 0) eq_start = j;
            if (t_hist[j] < t_equil) eq_end = j;
            if (t_hist[j] >= t_equil + t_test * 0.5 && ts_start < 0) ts_start = j;
        }
        ts_end = n_dft - 1;

        double omega_eq = -1, omega_ts = -1;
        if (eq_start >= 0 && eq_end > eq_start + 100) {
            /* Build sub-arrays for equilibrium phase */
            int neq = eq_end - eq_start + 1;
            double *h = phi0_hist + eq_start;
            double *th = t_hist + eq_start;
            double T = th[neq-1] - th[0];
            double peak_pow = 0;
            for (int k = 1; k < 500; k++) {
                double omega = 3.0 * mass * k / 500.0;
                double re = 0, im = 0;
                for (int j = 0; j < neq; j++) {
                    double dtj = (j > 0) ? (th[j]-th[j-1]) : (th[1]-th[0]);
                    re += h[j] * cos(omega * th[j]) * dtj;
                    im += h[j] * sin(omega * th[j]) * dtj;
                }
                double pw = (re*re + im*im) / (T*T);
                if (pw > peak_pow) { peak_pow = pw; omega_eq = omega; }
            }
            printf("  Equilibrium omega = %.4f (mass gap = %.4f)\n", omega_eq, mass);
        }
        if (ts_start >= 0 && ts_end > ts_start + 100) {
            int nts = ts_end - ts_start + 1;
            double *h = phi0_hist + ts_start;
            double *th = t_hist + ts_start;
            double T = th[nts-1] - th[0];
            double peak_pow = 0;
            for (int k = 1; k < 500; k++) {
                double omega = 3.0 * mass * k / 500.0;
                double re = 0, im = 0;
                for (int j = 0; j < nts; j++) {
                    double dtj = (j > 0) ? (th[j]-th[j-1]) : (th[1]-th[0]);
                    re += h[j] * cos(omega * th[j]) * dtj;
                    im += h[j] * sin(omega * th[j]) * dtj;
                }
                double pw = (re*re + im*im) / (T*T);
                if (pw > peak_pow) { peak_pow = pw; omega_ts = omega; }
            }
            printf("  Test-phase omega  = %.4f (mass gap = %.4f)\n", omega_ts, mass);
        }

        if (omega_eq > 0 && omega_ts > 0) {
            double shift = omega_ts - omega_eq;
            printf("  Frequency shift   = %+.6f (%.2f%%)\n",
                   shift, 100.0 * shift / omega_eq);
        }

        /* Write spectrum file */
        char specpath[600];
        snprintf(specpath, sizeof(specpath), "%s/fieldmetric_spectrum.tsv", outdir);
        FILE *fspec = fopen(specpath, "w");
        if (fspec) {
            fprintf(fspec, "omega\tpower_equil\tpower_test\n");
            for (int k = 0; k < 500; k++) {
                double omega = 3.0 * mass * k / 500.0;
                double pw_eq = 0, pw_ts = 0;

                if (eq_start >= 0 && eq_end > eq_start + 100) {
                    int neq = eq_end - eq_start + 1;
                    double *h = phi0_hist + eq_start;
                    double *th = t_hist + eq_start;
                    double T = th[neq-1] - th[0];
                    double re = 0, im = 0;
                    for (int j = 0; j < neq; j++) {
                        double dtj = (j > 0) ? (th[j]-th[j-1]) : (th[1]-th[0]);
                        re += h[j] * cos(omega * th[j]) * dtj;
                        im += h[j] * sin(omega * th[j]) * dtj;
                    }
                    pw_eq = (re*re + im*im) / (T*T);
                }
                if (ts_start >= 0 && ts_end > ts_start + 100) {
                    int nts = ts_end - ts_start + 1;
                    double *h = phi0_hist + ts_start;
                    double *th = t_hist + ts_start;
                    double T = th[nts-1] - th[0];
                    double re = 0, im = 0;
                    for (int j = 0; j < nts; j++) {
                        double dtj = (j > 0) ? (th[j]-th[j-1]) : (th[1]-th[0]);
                        re += h[j] * cos(omega * th[j]) * dtj;
                        im += h[j] * sin(omega * th[j]) * dtj;
                    }
                    pw_ts = (re*re + im*im) / (T*T);
                }
                fprintf(fspec, "%.6f\t%.6e\t%.6e\n", omega, pw_eq, pw_ts);
            }
            fclose(fspec);
            printf("  Spectrum: %s\n", specpath);
        }
    }

    /* Final summary */
    printf("\n=== Summary ===\n");
    printf("  Output: %s\n", tspath);
    if (fgeff) printf("  g_eff tracking: %s\n", gpath);
    if (fcausal) printf("  Causality data: %s\n", causalpath);

    /* Cleanup */
    for (int a = 0; a < 3; a++) { free(phi[a]); free(vel[a]); free(acc[a]); }
    free(rho_ema); free(g_eff); free(damp); free(osc_pos);
    free(phi0_hist); free(t_hist);
    return 0;
}
