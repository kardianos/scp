/*
 * dissipative1d.c — Dissipative soliton: gain-loss balanced oscillon
 *
 * Modified EOM (three massive scalars with triple-product coupling):
 *   d²φ_a/dt² = ∇²φ_a - m²φ_a - dV/dφ_a - γ·dφ_a/dt + G·S_a
 *
 * Gain modes:
 *   Mode A (proposal): S_a = P²·φ_a / (1 + κ_g·P²)  [alignment gain]
 *   Mode B (cubic):    S_a = Σ²·φ_a / (1 + κ_g·Σ²)  where Σ²=Σφ_a²  [amplitude gain]
 *   Mode C (saturated anti-damping): S_a = dφ_a/dt / (1 + α·Σ²)  [core anti-damping]
 *
 * Phase 1: Find G_balance for each γ (binary search).
 * Phase 2: Perturb ±20%, check limit-cycle recovery.
 * Phase 3: Start from noise.
 *
 * Compile: gcc -O3 -Wall -o dissipative1d src/dissipative1d.c -lm
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Physical parameters */
static double mu     = -20.0;
static double kappa  = 20.0;
static double mass   = 1.0;
static double A_init = 0.8;
static double sigma  = 3.0;

/* Gain */
static double kappa_g = 100.0;
static int gain_mode  = 0;  /* 0=A(P²), 1=B(Σ²), 2=C(anti-damp) */

/* Grid */
static int    Nx     = 4000;
static double xmax   = 100.0;

/* Time */
static double t_equil = 3000.0;
static double t_diss  = 10000.0;

/* Dissipation/gain */
static double gamma_damp = 0.01;
static double G_gain     = 0.0;

static int phase = 0;
static char outdir[512] = "data";

static void parse_args(int argc, char **argv)
{
    for (int i = 1; i < argc - 1; i += 2) {
        if      (!strcmp(argv[i], "-mu"))      mu         = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-kappa"))   kappa      = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-mass"))    mass       = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-A"))       A_init     = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-sigma"))   sigma      = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-kappa_g")) kappa_g    = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-gmode"))   gain_mode  = atoi(argv[i+1]);
        else if (!strcmp(argv[i], "-Nx"))      Nx         = atoi(argv[i+1]);
        else if (!strcmp(argv[i], "-xmax"))    xmax       = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-t_equil")) t_equil    = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-t_diss"))  t_diss     = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-gamma"))   gamma_damp = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-G"))       G_gain     = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-phase"))   phase      = atoi(argv[i+1]);
        else if (!strcmp(argv[i], "-o"))       strncpy(outdir, argv[i+1], sizeof(outdir)-1);
        else { fprintf(stderr, "Unknown flag: %s\n", argv[i]); exit(1); }
    }
}

/* -dV/dφ_a */
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

/* Gain source for different modes */
static inline double gain_source(double p1, double p2, double p3,
                                 double v1, double v2, double v3,
                                 int a, int mode)
{
    double phi_a, vel_a;
    switch (a) {
        case 0: phi_a = p1; vel_a = v1; break;
        case 1: phi_a = p2; vel_a = v2; break;
        case 2: phi_a = p3; vel_a = v3; break;
        default: phi_a = 0; vel_a = 0;
    }

    if (mode == 0) {
        /* Mode A: P² alignment gain (saturated) */
        double P = p1 * p2 * p3;
        double P2 = P * P;
        return P2 * phi_a / (1.0 + kappa_g * P2);
    } else if (mode == 1) {
        /* Mode B: Σ² amplitude gain (saturated) */
        double S2 = p1*p1 + p2*p2 + p3*p3;
        return S2 * phi_a / (1.0 + kappa_g * S2);
    } else {
        /* Mode C: saturated anti-damping */
        double S2 = p1*p1 + p2*p2 + p3*p3;
        return vel_a / (1.0 + kappa_g * S2);
    }
}

/* Diagnostics */
typedef struct {
    double E_kin, E_grad, E_mass, E_pot, E_total;
    double E_core, E_all, f_core;
    double peak[3];
    double P_diss, P_gain;
    int is_nan;
} Diag;

static Diag compute_diag(double **phi, double **vel, double dx,
                         double m2, double core_r, double G, double gamma_d)
{
    Diag d;
    memset(&d, 0, sizeof(d));
    for (int i = 1; i < Nx - 1; i++) {
        double x = -xmax + i * dx;
        for (int a = 0; a < 3; a++) {
            if (!isfinite(phi[a][i]) || !isfinite(vel[a][i])) {
                d.is_nan = 1; return d;
            }
            d.E_kin += 0.5 * vel[a][i] * vel[a][i] * dx;
            double dp = (phi[a][i+1] - phi[a][i-1]) / (2.0*dx);
            d.E_grad += 0.5 * dp * dp * dx;
            d.E_mass += 0.5 * m2 * phi[a][i] * phi[a][i] * dx;
            if (fabs(phi[a][i]) > d.peak[a]) d.peak[a] = fabs(phi[a][i]);
        }

        double P = phi[0][i] * phi[1][i] * phi[2][i];
        double P2 = P * P;
        double V = 0.5 * mu * P2 / (1.0 + kappa * P2);
        d.E_pot += V * dx;

        for (int a = 0; a < 3; a++) {
            d.P_diss += gamma_d * vel[a][i] * vel[a][i] * dx;
            double gs = gain_source(phi[0][i], phi[1][i], phi[2][i],
                                    vel[0][i], vel[1][i], vel[2][i],
                                    a, gain_mode);
            d.P_gain += G * gs * vel[a][i] * dx;
        }

        double e = V;
        for (int a = 0; a < 3; a++) {
            e += 0.5*vel[a][i]*vel[a][i];
            double dp = (phi[a][i+1] - phi[a][i-1]) / (2.0*dx);
            e += 0.5*dp*dp + 0.5*m2*phi[a][i]*phi[a][i];
        }
        d.E_all += e * dx;
        if (fabs(x) < core_r) d.E_core += e * dx;
    }
    d.E_total = d.E_kin + d.E_grad + d.E_mass + d.E_pot;
    d.f_core = (d.E_all > 1e-20) ? d.E_core / d.E_all : 0.0;
    return d;
}

/* RNG */
static unsigned long rng_state = 12345;
static double rand_uniform(void)
{
    rng_state = rng_state * 6364136223846793005ULL + 1442695040888963407ULL;
    return (double)(rng_state >> 33) / (double)(1ULL << 31);
}
static double rand_gauss(void)
{
    double u1 = rand_uniform() + 1e-30;
    double u2 = rand_uniform();
    return sqrt(-2.0*log(u1)) * cos(2.0*M_PI*u2);
}

/*
 * Evolve fields. Returns 0 on success, 1 on blowup.
 */
static int evolve(double **phi, double **vel, int init_mode,
                  double gamma_d, double G, double t_total,
                  FILE *tsfile, double *E_out)
{
    double dx = 2.0 * xmax / (Nx - 1);
    double dx2 = dx * dx;
    double m2 = mass * mass;

    double kmax = M_PI / dx;
    double cfl_fac = (G > 0.0) ? 0.4 : 0.8;
    double dt = cfl_fac * 2.0 / sqrt(kmax*kmax + m2);
    int Nt = (int)(t_total / dt) + 1;

    double *acc[3];
    for (int a = 0; a < 3; a++)
        acc[a] = calloc(Nx, sizeof(double));

    /* Absorbing boundary */
    double *damp_bd = malloc(Nx * sizeof(double));
    double x_abs = xmax * 0.75;
    for (int i = 0; i < Nx; i++) {
        double x = -xmax + i * dx;
        double ax = fabs(x);
        if (ax > x_abs) {
            double f = (ax - x_abs) / (xmax - x_abs);
            damp_bd[i] = 1.0 - 0.98 * f * f;
        } else {
            damp_bd[i] = 1.0;
        }
    }

    if (init_mode == 0) {
        for (int a = 0; a < 3; a++)
            for (int i = 0; i < Nx; i++) {
                double x = -xmax + i * dx;
                phi[a][i] = A_init * exp(-x*x / (2.0*sigma*sigma));
                vel[a][i] = 0.0;
            }
    } else if (init_mode == 2) {
        double noise_amp = 0.05;
        for (int a = 0; a < 3; a++)
            for (int i = 0; i < Nx; i++) {
                phi[a][i] = noise_amp * rand_gauss();
                vel[a][i] = noise_amp * rand_gauss() * 0.5;
            }
    }

    /* Macro for acceleration computation */
    #define COMPUTE_ACC() do { \
        for (int a = 0; a < 3; a++) { \
            acc[a][0] = acc[a][Nx-1] = 0; \
            for (int i = 1; i < Nx - 1; i++) { \
                double lapl = (phi[a][i+1] - 2.0*phi[a][i] + phi[a][i-1]) / dx2; \
                double fp = force_pot(phi[0][i], phi[1][i], phi[2][i], a); \
                double gs = G * gain_source(phi[0][i], phi[1][i], phi[2][i], \
                                            vel[0][i], vel[1][i], vel[2][i], \
                                            a, gain_mode); \
                acc[a][i] = lapl - m2*phi[a][i] + fp - gamma_d*vel[a][i] + gs; \
            } \
        } \
    } while(0)

    COMPUTE_ACC();

    int rec_every = Nt / 10000;
    if (rec_every < 1) rec_every = 1;
    int ic = Nx / 2;
    double core_r = 3.0 * sigma;
    double E_final = 0;
    int blowup = 0;
    int check_every = Nt / 1000;
    if (check_every < 1) check_every = 1;

    for (int n = 0; n <= Nt; n++) {
        double t = n * dt;

        if (n % rec_every == 0) {
            Diag d = compute_diag(phi, vel, dx, m2, core_r, G, gamma_d);
            if (d.is_nan) { blowup = 1; break; }
            E_final = d.E_total;
            if (tsfile)
                fprintf(tsfile, "%.4f\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t"
                                "%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.4f\t%.6e\t%.6e\n",
                        t, phi[0][ic], phi[1][ic], phi[2][ic],
                        d.peak[0], d.peak[1], d.peak[2],
                        d.E_kin, d.E_grad, d.E_mass, d.E_pot, d.E_total,
                        d.f_core, d.P_diss, d.P_gain);
        }

        if (n % check_every == 0 && n > 0) {
            if (!isfinite(phi[0][ic]) || !isfinite(vel[0][ic]) || fabs(phi[0][ic]) > 1e6) {
                blowup = 1; break;
            }
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

        for (int a = 0; a < 3; a++)
            for (int i = 0; i < Nx; i++) {
                vel[a][i] *= damp_bd[i];
                phi[a][i] *= damp_bd[i];
            }
    }

    #undef COMPUTE_ACC

    for (int a = 0; a < 3; a++) free(acc[a]);
    free(damp_bd);

    if (E_out) *E_out = E_final;
    return blowup;
}

/* Copy state arrays */
static void copy_state(double **dst_phi, double **dst_vel,
                       double **src_phi, double **src_vel)
{
    for (int a = 0; a < 3; a++) {
        memcpy(dst_phi[a], src_phi[a], Nx * sizeof(double));
        memcpy(dst_vel[a], src_vel[a], Nx * sizeof(double));
    }
}

/* Allocate 3-field arrays */
static void alloc_fields(double **phi, double **vel)
{
    for (int a = 0; a < 3; a++) {
        phi[a] = calloc(Nx, sizeof(double));
        vel[a] = calloc(Nx, sizeof(double));
    }
}
static void free_fields(double **phi, double **vel)
{
    for (int a = 0; a < 3; a++) { free(phi[a]); free(vel[a]); }
}

/* ============================================================ */
/* Phase 1: Balance scan                                        */
/* ============================================================ */
static void phase1_balance_scan(void)
{
    const char *mode_names[] = {"A(P^2)", "B(Sigma^2)", "C(anti-damp)"};
    printf("\n=== PHASE 1: Gain-Loss Balance Scan [mode %s, kappa_g=%.1f] ===\n",
           mode_names[gain_mode], kappa_g);

    double gammas[] = {0.001, 0.01, 0.1};
    int n_gamma = 3;

    char sumpath[600];
    snprintf(sumpath, sizeof(sumpath), "%s/dissipative_summary.tsv", outdir);
    FILE *fsum = fopen(sumpath, "w");
    if (!fsum) { fprintf(stderr, "Cannot open %s\n", sumpath); return; }
    fprintf(fsum, "gamma\tG_balance\tE_initial\tE_final\tf_core\tpeak\tnote\tgain_mode\n");

    for (int gi = 0; gi < n_gamma; gi++) {
        double gam = gammas[gi];
        printf("\n--- gamma = %.4f ---\n", gam);

        double *phi_eq[3], *vel_eq[3];
        alloc_fields(phi_eq, vel_eq);

        printf("  Equilibrating (t=%.0f)...\n", t_equil);
        evolve(phi_eq, vel_eq, 0, 0.0, 0.0, t_equil, NULL, NULL);

        double dx = 2.0 * xmax / (Nx - 1);
        double m2 = mass * mass;
        double core_r = 3.0 * sigma;
        Diag d0 = compute_diag(phi_eq, vel_eq, dx, m2, core_r, 0, 0);
        double E0 = d0.E_total;
        printf("  Equilibrated: E=%.4f, peak=%.4f, fc=%.4f\n",
               E0, d0.peak[0], d0.f_core);

        /* Find upper bound for G */
        double G_lo = 0.0, G_hi = 1.0;
        double test_time = 500.0;

        int found_upper = 0;
        for (int attempt = 0; attempt < 20 && !found_upper; attempt++) {
            double *phi_t[3], *vel_t[3];
            alloc_fields(phi_t, vel_t);
            copy_state(phi_t, vel_t, phi_eq, vel_eq);
            double Ef;
            int blow = evolve(phi_t, vel_t, 1, gam, G_hi, test_time, NULL, &Ef);
            free_fields(phi_t, vel_t);
            if (blow || Ef > E0 * 1.1 || !isfinite(Ef)) {
                found_upper = 1;
                printf("  Upper bound: G=%.6f (blow=%d, E=%.4f)\n", G_hi, blow, Ef);
            } else {
                G_hi *= 3.0;
            }
        }

        if (!found_upper) {
            printf("  Cannot find upper bound. Skipping.\n");
            fprintf(fsum, "%.4f\t-1\t%.4f\t-1\t-1\t-1\tno_upper_bound\t%s\n",
                    gam, E0, mode_names[gain_mode]);
            free_fields(phi_eq, vel_eq);
            continue;
        }

        /* Binary search */
        double G_best = 0.0, E_best = 0.0, fc_best = 0.0, pk_best = 0.0;
        int best_balanced = 0;

        for (int iter = 0; iter < 30; iter++) {
            double G_mid = 0.5 * (G_lo + G_hi);

            double *phi_t[3], *vel_t[3];
            alloc_fields(phi_t, vel_t);
            copy_state(phi_t, vel_t, phi_eq, vel_eq);

            double Ef;
            int blow = evolve(phi_t, vel_t, 1, gam, G_mid, test_time, NULL, &Ef);
            Diag df = {0};
            if (!blow) df = compute_diag(phi_t, vel_t, dx, m2, core_r, G_mid, gam);
            free_fields(phi_t, vel_t);

            if (blow || df.is_nan || !isfinite(Ef)) {
                G_hi = G_mid;
                if (iter % 5 == 0)
                    printf("  iter=%2d: G=%.8f → BLOWUP\n", iter, G_mid);
            } else {
                double ratio = df.E_total / E0;

                if (iter % 5 == 0)
                    printf("  iter=%2d: G=%.8f, E=%.4f, ratio=%.4f, pk=%.4f, fc=%.4f\n",
                           iter, G_mid, df.E_total, ratio, df.peak[0], df.f_core);

                G_best = G_mid;
                E_best = df.E_total;
                fc_best = df.f_core;
                pk_best = df.peak[0];

                if (ratio > 0.85 && ratio < 1.15) {
                    best_balanced = 1;
                    /* Narrow search */
                    double w = (G_hi - G_lo) * 0.3;
                    G_lo = G_mid - w;
                    G_hi = G_mid + w;
                    if (G_lo < 0) G_lo = 0;
                } else if (ratio > 1.15) {
                    G_hi = G_mid;
                } else {
                    G_lo = G_mid;
                }
            }
        }

        printf("  G_balance = %.8f\n", G_best);

        /* Longer validation: 2000 t.u. */
        printf("  Validation run (t=2000)...\n");
        {
            double *phi_t[3], *vel_t[3];
            alloc_fields(phi_t, vel_t);
            copy_state(phi_t, vel_t, phi_eq, vel_eq);

            char tspath[600];
            snprintf(tspath, sizeof(tspath), "%s/dissipative_g%.4f_G%.6f_mode%d_ts.tsv",
                     outdir, gam, G_best, gain_mode);
            FILE *fts = fopen(tspath, "w");
            if (fts) {
                fprintf(fts, "time\tphi1_0\tphi2_0\tphi3_0\tpeak1\tpeak2\tpeak3\t"
                             "E_kin\tE_grad\tE_mass\tE_pot\tE_total\tf_core\tP_diss\tP_gain\n");
                double Ef;
                int blow = evolve(phi_t, vel_t, 1, gam, G_best, 2000.0, fts, &Ef);
                fclose(fts);
                Diag df = compute_diag(phi_t, vel_t, dx, m2, core_r, G_best, gam);
                if (blow || df.is_nan)
                    printf("  Validation: BLOWUP\n");
                else {
                    printf("  Validation: E=%.4f (E0=%.4f), ratio=%.4f, pk=%.4f, fc=%.4f\n",
                           df.E_total, E0, df.E_total/E0, df.peak[0], df.f_core);
                    E_best = df.E_total;
                    fc_best = df.f_core;
                    pk_best = df.peak[0];
                    if (df.E_total/E0 > 0.7 && df.E_total/E0 < 1.3 && df.f_core > 0.5)
                        best_balanced = 1;
                }
                printf("  Output: %s\n", tspath);
            }
            free_fields(phi_t, vel_t);
        }

        const char *note = best_balanced ? "balanced" :
                           (E_best < 0.01*E0) ? "decayed" :
                           (E_best > 5*E0) ? "blowup" : "marginal";

        fprintf(fsum, "%.4f\t%.8f\t%.4f\t%.4f\t%.4f\t%.4f\t%s\t%s\n",
                gam, G_best, E0, E_best, fc_best, pk_best, note,
                mode_names[gain_mode]);
        printf("  => gamma=%.4f, G=%.8f, status=%s\n\n", gam, G_best, note);

        free_fields(phi_eq, vel_eq);
    }

    fclose(fsum);
    printf("Phase 1 summary: %s\n", sumpath);
}

/* ============================================================ */
/* Phase 2                                                      */
/* ============================================================ */
static void phase2_stability(double gam, double G)
{
    printf("\n=== PHASE 2: Stability (gamma=%.4f, G=%.8f, mode=%d) ===\n",
           gam, G, gain_mode);

    double dx = 2.0 * xmax / (Nx - 1);
    double m2 = mass * mass;
    double core_r = 3.0 * sigma;

    double *phi_eq[3], *vel_eq[3];
    alloc_fields(phi_eq, vel_eq);

    printf("  Equilibrating...\n");
    evolve(phi_eq, vel_eq, 0, 0.0, 0.0, t_equil, NULL, NULL);

    printf("  Reaching dissipative steady state (t=2000)...\n");
    int blow = evolve(phi_eq, vel_eq, 1, gam, G, 2000.0, NULL, NULL);
    if (blow) {
        printf("  BLOWUP. Aborting.\n");
        free_fields(phi_eq, vel_eq);
        return;
    }

    Diag d0 = compute_diag(phi_eq, vel_eq, dx, m2, core_r, G, gam);
    printf("  Steady state: E=%.4f, peak=%.4f, fc=%.4f\n",
           d0.E_total, d0.peak[0], d0.f_core);

    double perturbs[] = {1.2, 0.8};
    const char *pnames[] = {"plus20", "minus20"};

    for (int pi = 0; pi < 2; pi++) {
        double fac = perturbs[pi];
        printf("\n  Perturbation: x%.1f\n", fac);

        double *phi_t[3], *vel_t[3];
        alloc_fields(phi_t, vel_t);
        for (int a = 0; a < 3; a++)
            for (int i = 0; i < Nx; i++) {
                phi_t[a][i] = phi_eq[a][i] * fac;
                vel_t[a][i] = vel_eq[a][i] * fac;
            }

        char tspath[600];
        snprintf(tspath, sizeof(tspath), "%s/dissipative_stability_%s_mode%d_ts.tsv",
                 outdir, pnames[pi], gain_mode);
        FILE *fts = fopen(tspath, "w");
        fprintf(fts, "time\tphi1_0\tphi2_0\tphi3_0\tpeak1\tpeak2\tpeak3\t"
                     "E_kin\tE_grad\tE_mass\tE_pot\tE_total\tf_core\tP_diss\tP_gain\n");

        double Ef;
        int bl = evolve(phi_t, vel_t, 1, gam, G, 5000.0, fts, &Ef);
        fclose(fts);

        if (bl) {
            printf("  BLOWUP\n");
        } else {
            Diag df = compute_diag(phi_t, vel_t, dx, m2, core_r, G, gam);
            printf("  Final: E=%.4f (target=%.4f), ratio=%.4f\n",
                   df.E_total, d0.E_total,
                   (d0.E_total > 1e-10) ? df.E_total/d0.E_total : 0.0);
            printf("  Recovery: %s\n",
                   (d0.E_total > 1e-10 && fabs(df.E_total/d0.E_total - 1.0) < 0.15)
                   ? "YES (limit cycle)" : "NO");
        }
        printf("  Output: %s\n", tspath);

        free_fields(phi_t, vel_t);
    }

    free_fields(phi_eq, vel_eq);
}

/* ============================================================ */
/* Phase 3                                                      */
/* ============================================================ */
static void phase3_from_noise(double gam, double G)
{
    printf("\n=== PHASE 3: From Noise (gamma=%.4f, G=%.8f, mode=%d) ===\n",
           gam, G, gain_mode);

    double G_tests[] = {G, G*2.0, G*5.0, G*10.0};
    int n_tests = 4;

    double dx = 2.0 * xmax / (Nx - 1);
    double m2 = mass * mass;
    double core_r = 3.0 * sigma;

    for (int ti = 0; ti < n_tests; ti++) {
        double Gt = G_tests[ti];
        printf("\n  --- G=%.6f (%.1fx balance) ---\n", Gt, Gt/G);

        double *phi_t[3], *vel_t[3];
        alloc_fields(phi_t, vel_t);

        char tspath[600];
        snprintf(tspath, sizeof(tspath), "%s/dissipative_fromNoise_G%.4f_mode%d_ts.tsv",
                 outdir, Gt, gain_mode);
        FILE *fts = fopen(tspath, "w");
        fprintf(fts, "time\tphi1_0\tphi2_0\tphi3_0\tpeak1\tpeak2\tpeak3\t"
                     "E_kin\tE_grad\tE_mass\tE_pot\tE_total\tf_core\tP_diss\tP_gain\n");

        rng_state = 42 + ti;
        double Ef;
        int blow = evolve(phi_t, vel_t, 2, gam, Gt, t_diss, fts, &Ef);
        fclose(fts);

        if (blow) {
            printf("  BLOWUP\n");
        } else {
            Diag df = compute_diag(phi_t, vel_t, dx, m2, core_r, Gt, gam);
            printf("  Final: E=%.4f, peak=%.4f, fc=%.4f\n",
                   df.E_total, df.peak[0], df.f_core);
            printf("  Formed oscillon? %s\n",
                   (df.f_core > 0.3 && df.peak[0] > 0.05) ? "YES" : "NO");
        }
        printf("  Output: %s\n", tspath);

        free_fields(phi_t, vel_t);
    }
}

/* ============================================================ */
/* Main                                                         */
/* ============================================================ */
int main(int argc, char **argv)
{
    parse_args(argc, argv);

    const char *mode_names[] = {"A(P^2)", "B(Sigma^2)", "C(anti-damp)"};
    printf("dissipative1d: mu=%.1f kappa=%.1f mass=%.1f kappa_g=%.1f mode=%s\n",
           mu, kappa, mass, kappa_g, mode_names[gain_mode]);
    printf("  Nx=%d xmax=%.1f t_equil=%.0f t_diss=%.0f\n", Nx, xmax, t_equil, t_diss);

    if (phase == 0 || phase == 1) {
        phase1_balance_scan();
    }

    /* Read G_balance */
    double gam_test = 0.01;
    double G_test = G_gain;

    if ((phase == 0 || phase == 2 || phase == 3) && G_test <= 0.0) {
        char sumpath[600];
        snprintf(sumpath, sizeof(sumpath), "%s/dissipative_summary.tsv", outdir);
        FILE *fs = fopen(sumpath, "r");
        if (fs) {
            char line[1024];
            if (fgets(line, sizeof(line), fs)) {
                while (fgets(line, sizeof(line), fs)) {
                    double g, Gb;
                    if (sscanf(line, "%lf\t%lf", &g, &Gb) == 2) {
                        if (fabs(g - gam_test) < 1e-6) {
                            G_test = Gb;
                            break;
                        }
                    }
                }
            }
            fclose(fs);
        }
        if (G_test <= 0.0) G_test = 100.0;
    }

    if (phase == 0 || phase == 2)
        phase2_stability(gam_test, G_test);

    if (phase == 0 || phase == 3)
        phase3_from_noise(gam_test, G_test);

    printf("\n=== DONE ===\n");
    return 0;
}
