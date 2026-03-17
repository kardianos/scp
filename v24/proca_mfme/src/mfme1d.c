/*
 * mfme1d.c — V24-P5: Combined MF+ME (Goldstone backreaction + Pairwise coupling)
 *
 * Five field arrays: phi_1, phi_2, phi_3, theta (Goldstone)
 *
 * Lagrangian:
 *   L = sum_a [ (1/2)(dt phi_a)^2 - (1/2)(dx phi_a)^2 - (m^2/2)phi_a^2 ]
 *     - (mu/2) P^2 / (1 + kappa P^2)           [triple-product]
 *     - lambda(phi_1*phi_2 + phi_2*phi_3 + phi_3*phi_1)  [pairwise]
 *     + (1/2)(dt theta)^2 - (1/2)(dx theta)^2   [Goldstone KE]
 *     + g * (dt theta) * P^2 / (1 + kappa*P^2)^2 [theta coupling]
 *
 * EOM:
 *   d_tt phi_a = d_xx phi_a - m^2 phi_a - lambda(phi_b + phi_c)
 *              - dV_triple/dphi_a
 *              + 2*g*(d_t theta)*P*(dP/dphi_a)/(1+kappa*P^2)^2   [backreaction]
 *
 *   d_tt theta = d_xx theta - g * d_t(P^2/(1+kappa*P^2)^2)        [source from P]
 *     simplified for small kappa*P^2: d_tt theta ~ d_xx theta - g * d_t(P^2)
 *     We use the FULL regulated form: source = d_t[ P^2/(1+kappa*P^2)^2 ]
 *
 * Phases:
 *   1: Baselines: (g=1,lam=0) then (g=0,lam=0.99) — each to tfinal
 *   2: Combined: (g=1,lam=0.99)
 *   3: Push lambda past 1.0 with g=1: scan {0.99, 0.995, 1.0, 1.01, 1.05}
 *
 * Compile: gcc -O3 -Wall -o mfme1d src/mfme1d.c -lm
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Physics parameters */
static double mu     = -20.0;
static double kappa  = 20.0;
static double mass   = 1.0;
static double A_init = 0.8;
static double sigma  = 3.0;
static double g_coup = 1.0;     /* theta coupling */
static double lambda = 0.0;     /* pairwise coupling */

/* Grid */
static int    Nx     = 8000;
static double xmax   = 200.0;
static double tfinal = 10000.0;

static int    run_phase = 0;  /* 0=single run, 1=baselines, 2=combined, 3=push scan */
static char   outdir[512] = "v24/proca_mfme/data";

static void parse_args(int argc, char **argv)
{
    for (int i = 1; i < argc - 1; i += 2) {
        if      (!strcmp(argv[i], "-mu"))     mu     = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-kappa"))  kappa  = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-mass"))   mass   = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-A"))      A_init = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-sigma"))  sigma  = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-g"))      g_coup = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-lambda")) lambda = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-Nx"))     Nx     = atoi(argv[i+1]);
        else if (!strcmp(argv[i], "-xmax"))   xmax   = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-tfinal")) tfinal = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-phase"))  run_phase = atoi(argv[i+1]);
        else if (!strcmp(argv[i], "-o"))      strncpy(outdir, argv[i+1], sizeof(outdir)-1);
        else { fprintf(stderr, "Unknown flag: %s\n", argv[i]); exit(1); }
    }
}

/* -dV_triple/dphi_a where V = (mu/2)P^2/(1+kappa*P^2), P = phi1*phi2*phi3 */
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

/* Regulated P^2 factor: R(P) = P^2 / (1 + kappa*P^2)^2 */
static inline double reg_P2(double P)
{
    double P2 = P * P;
    double d = 1.0 + kappa * P2;
    return P2 / (d * d);
}

/* ======================================================================
 * Single evolution run: returns results via struct
 * ====================================================================== */
typedef struct {
    double omega;       /* peak DFT frequency of phi_1(0,t) */
    double E_final;     /* total energy at end */
    double E_phi;       /* phi-sector energy at end */
    double E_theta;     /* theta-sector energy at end */
    double fc_final;    /* core fraction at end */
    double peak_phi;    /* peak |phi| at end */
    double max_theta;   /* max |theta| at end */
    int    survived;    /* 1 if fc > 0.5 at end */
} RunResult;

static RunResult run_evolution(double g, double lam, const char *label)
{
    RunResult res = {0};

    double dx  = 2.0 * xmax / (Nx - 1);
    double dx2 = dx * dx;
    double m2  = mass * mass;

    /* CFL: account for massless theta (c=1) and massive phi */
    double kmax = M_PI / dx;
    double dt_mass = 0.8 * 2.0 / sqrt(kmax * kmax + m2 + 2.0 * fabs(lam));
    double dt_wave = 0.8 * dx;
    double dt = (dt_mass < dt_wave) ? dt_mass : dt_wave;
    int Nt = (int)(tfinal / dt) + 1;

    double m2_anti = m2 - lam;
    double m_A = (m2_anti > 0) ? sqrt(m2_anti) : 0.0;

    printf("\n=== RUN: %s ===\n", label);
    printf("  g=%.4f lambda=%.4f m2_anti=%.4f m_A=%.4f\n", g, lam, m2_anti, m_A);
    printf("  mu=%.1f kappa=%.1f mass=%.4f A=%.3f sigma=%.3f\n",
           mu, kappa, mass, A_init, sigma);
    printf("  Nx=%d xmax=%.1f dx=%.5f dt=%.6f Nt=%d tfinal=%.0f\n",
           Nx, xmax, dx, dt, Nt, tfinal);
    fflush(stdout);

    /* Allocate: 3 phi + 1 theta */
    double *phi[3], *vel[3], *acc[3];
    double *theta, *vth, *ath;
    double *RP_cur, *RP_prev;  /* R(P) = P^2/(1+kappa*P^2)^2 at current/previous step */

    for (int a = 0; a < 3; a++) {
        phi[a] = calloc(Nx, sizeof(double));
        vel[a] = calloc(Nx, sizeof(double));
        acc[a] = calloc(Nx, sizeof(double));
    }
    theta   = calloc(Nx, sizeof(double));
    vth     = calloc(Nx, sizeof(double));
    ath     = calloc(Nx, sizeof(double));
    RP_cur  = calloc(Nx, sizeof(double));
    RP_prev = calloc(Nx, sizeof(double));

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

    /* Initialize: Gaussian oscillon for all three phi, theta=0 */
    for (int a = 0; a < 3; a++)
        for (int i = 0; i < Nx; i++) {
            double x = -xmax + i * dx;
            phi[a][i] = A_init * exp(-x * x / (2.0 * sigma * sigma));
        }

    /* Initialize RP */
    for (int i = 0; i < Nx; i++) {
        double P = phi[0][i] * phi[1][i] * phi[2][i];
        RP_cur[i] = reg_P2(P);
        RP_prev[i] = RP_cur[i];
    }

    /* Acceleration computation macro:
     * phi_a: lapl - m^2*phi_a - lambda*(phi_b + phi_c) - dV_triple/dphi_a
     *        + 2*g*(d_t theta)*P*(dP/dphi_a)/(1+kappa*P^2)^2  [backreaction]
     * theta: lapl_theta - g * d_t[ P^2/(1+kappa*P^2)^2 ]
     *        i.e. d_tt theta = d_xx theta - g*(RP_cur - RP_prev)/dt
     */

    int ic = Nx / 2;
    double core_r = 3.0 * sigma;

    /* Time series output */
    char tspath[600];
    snprintf(tspath, sizeof(tspath), "%s/mfme_%s_ts.tsv", outdir, label);
    FILE *fts = fopen(tspath, "w");
    if (!fts) { fprintf(stderr, "Cannot open %s\n", tspath); return res; }
    fprintf(fts, "time\tphi1_0\tpeak_phi\tE_phi\tE_theta\tE_total\tfc\t"
                 "theta_0\tmax_theta\n");

    /* DFT storage */
    int max_dft = 50000;
    double *phi0_hist = malloc(max_dft * sizeof(double));
    double *t_hist    = malloc(max_dft * sizeof(double));
    int n_dft = 0;

    int rec_every = Nt / 20000;
    if (rec_every < 1) rec_every = 1;
    int print_every = Nt / 40;
    if (print_every < 1) print_every = 1;
    int dft_every = Nt / max_dft;
    if (dft_every < 1) dft_every = 1;

    /* Main loop */
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

        if (do_rec || do_print) {
            double Ek = 0, Eg = 0, Em = 0, Ep = 0, Epw = 0;
            double Eth_kin = 0, Eth_grad = 0;
            double Ecore = 0, Eall = 0;
            double pk_phi = 0, mx_theta = 0;

            for (int i = 1; i < Nx - 1; i++) {
                double x = -xmax + i * dx;

                /* phi sector */
                for (int a = 0; a < 3; a++) {
                    Ek += 0.5 * vel[a][i] * vel[a][i] * dx;
                    double dp = (phi[a][i+1] - phi[a][i-1]) / (2.0*dx);
                    Eg += 0.5 * dp * dp * dx;
                    Em += 0.5 * m2 * phi[a][i] * phi[a][i] * dx;
                    if (fabs(phi[a][i]) > pk_phi) pk_phi = fabs(phi[a][i]);
                }

                /* Pairwise potential energy */
                Epw += lam * (phi[0][i]*phi[1][i] + phi[1][i]*phi[2][i]
                             + phi[2][i]*phi[0][i]) * dx;

                /* Triple potential energy */
                double P = phi[0][i] * phi[1][i] * phi[2][i];
                double P2 = P * P;
                double V = 0.5 * mu * P2 / (1.0 + kappa * P2);
                Ep += V * dx;

                /* theta sector */
                Eth_kin  += 0.5 * vth[i] * vth[i] * dx;
                double dth = (theta[i+1] - theta[i-1]) / (2.0*dx);
                Eth_grad += 0.5 * dth * dth * dx;
                if (fabs(theta[i]) > mx_theta) mx_theta = fabs(theta[i]);

                /* Total density for core fraction (phi sector only) */
                double e = V;
                for (int a = 0; a < 3; a++) {
                    e += 0.5*vel[a][i]*vel[a][i];
                    double dp2 = (phi[a][i+1] - phi[a][i-1]) / (2.0*dx);
                    e += 0.5*dp2*dp2 + 0.5*m2*phi[a][i]*phi[a][i];
                }
                e += lam * (phi[0][i]*phi[1][i] + phi[1][i]*phi[2][i]
                           + phi[2][i]*phi[0][i]);
                Eall += e * dx;
                if (fabs(x) < core_r) Ecore += e * dx;
            }

            double E_p = Ek + Eg + Em + Ep + Epw;
            double E_t = Eth_kin + Eth_grad;
            double Et  = E_p + E_t;
            double fc  = (Eall > 1e-20) ? Ecore / Eall : 0.0;

            if (do_rec)
                fprintf(fts, "%.6f\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6f\t%.6e\t%.6e\n",
                        t, phi[0][ic], pk_phi, E_p, E_t, Et, fc,
                        theta[ic], mx_theta);

            if (do_print)
                printf("  t=%7.0f  p0=%+.4f  pk=%.3f  Ephi=%.3f  Eth=%.3e  fc=%.3f  th0=%.2e\n",
                       t, phi[0][ic], pk_phi, E_p, E_t, fc, theta[ic]);

            /* Store final values */
            res.E_phi    = E_p;
            res.E_theta  = E_t;
            res.E_final  = Et;
            res.fc_final = fc;
            res.peak_phi = pk_phi;
            res.max_theta = mx_theta;
        }

        if (n == Nt) break;

        /* --- Velocity Verlet --- */

        /* 1. Half-kick velocities */
        /* Compute accelerations inline for half-kick */
        /* First full acceleration computation */
        for (int a = 0; a < 3; a++) {
            acc[a][0] = acc[a][1] = acc[a][Nx-2] = acc[a][Nx-1] = 0;
            int b = (a+1)%3, c = (a+2)%3;
            for (int i = 1; i < Nx - 1; i++) {
                double lapl = (phi[a][i+1] - 2.0*phi[a][i] + phi[a][i-1]) / dx2;
                double fp = force_triple(phi[0][i], phi[1][i], phi[2][i], a);

                /* Backreaction: +2*g*(d_t theta)*P*(dP/dphi_a)/(1+kappa*P^2)^2 */
                double P_loc = phi[0][i] * phi[1][i] * phi[2][i];
                double P2_loc = P_loc * P_loc;
                double denom2 = (1.0 + kappa * P2_loc);
                denom2 = denom2 * denom2;
                double dPda;
                switch (a) {
                    case 0: dPda = phi[1][i] * phi[2][i]; break;
                    case 1: dPda = phi[0][i] * phi[2][i]; break;
                    case 2: dPda = phi[0][i] * phi[1][i]; break;
                    default: dPda = 0.0;
                }
                double f_back = 2.0 * g * vth[i] * P_loc * dPda / denom2;

                acc[a][i] = lapl - m2*phi[a][i] - lam*(phi[b][i]+phi[c][i]) + fp + f_back;
            }
        }
        /* theta acceleration: d_tt theta = d_xx theta - g * d_t(RP) */
        ath[0] = ath[1] = ath[Nx-2] = ath[Nx-1] = 0;
        for (int i = 1; i < Nx - 1; i++) {
            double lapl_th = (theta[i+1] - 2.0*theta[i] + theta[i-1]) / dx2;
            double dRPdt = (RP_cur[i] - RP_prev[i]) / dt;
            ath[i] = lapl_th - g * dRPdt;
        }

        /* Half-kick */
        for (int a = 0; a < 3; a++)
            for (int i = 1; i < Nx - 1; i++)
                vel[a][i] += 0.5 * dt * acc[a][i];
        for (int i = 1; i < Nx - 1; i++)
            vth[i] += 0.5 * dt * ath[i];

        /* 2. Drift positions */
        for (int a = 0; a < 3; a++)
            for (int i = 1; i < Nx - 1; i++)
                phi[a][i] += dt * vel[a][i];
        for (int i = 1; i < Nx - 1; i++)
            theta[i] += dt * vth[i];

        /* 3. Update RP history */
        for (int i = 0; i < Nx; i++)
            RP_prev[i] = RP_cur[i];
        for (int i = 0; i < Nx; i++) {
            double P = phi[0][i] * phi[1][i] * phi[2][i];
            RP_cur[i] = reg_P2(P);
        }

        /* 4. Recompute accelerations at new positions */
        for (int a = 0; a < 3; a++) {
            acc[a][0] = acc[a][1] = acc[a][Nx-2] = acc[a][Nx-1] = 0;
            int b = (a+1)%3, c = (a+2)%3;
            for (int i = 1; i < Nx - 1; i++) {
                double lapl = (phi[a][i+1] - 2.0*phi[a][i] + phi[a][i-1]) / dx2;
                double fp = force_triple(phi[0][i], phi[1][i], phi[2][i], a);

                double P_loc = phi[0][i] * phi[1][i] * phi[2][i];
                double P2_loc = P_loc * P_loc;
                double denom2 = (1.0 + kappa * P2_loc);
                denom2 = denom2 * denom2;
                double dPda;
                switch (a) {
                    case 0: dPda = phi[1][i] * phi[2][i]; break;
                    case 1: dPda = phi[0][i] * phi[2][i]; break;
                    case 2: dPda = phi[0][i] * phi[1][i]; break;
                    default: dPda = 0.0;
                }
                double f_back = 2.0 * g * vth[i] * P_loc * dPda / denom2;

                acc[a][i] = lapl - m2*phi[a][i] - lam*(phi[b][i]+phi[c][i]) + fp + f_back;
            }
        }
        ath[0] = ath[1] = ath[Nx-2] = ath[Nx-1] = 0;
        for (int i = 1; i < Nx - 1; i++) {
            double lapl_th = (theta[i+1] - 2.0*theta[i] + theta[i-1]) / dx2;
            double dRPdt = (RP_cur[i] - RP_prev[i]) / dt;
            ath[i] = lapl_th - g * dRPdt;
        }

        /* 5. Second half-kick */
        for (int a = 0; a < 3; a++)
            for (int i = 1; i < Nx - 1; i++)
                vel[a][i] += 0.5 * dt * acc[a][i];
        for (int i = 1; i < Nx - 1; i++)
            vth[i] += 0.5 * dt * ath[i];

        /* 6. Absorbing boundary */
        for (int a = 0; a < 3; a++)
            for (int i = 0; i < Nx; i++) {
                vel[a][i] *= damp[i];
                phi[a][i] *= damp[i];
            }
        for (int i = 0; i < Nx; i++) {
            vth[i]   *= damp[i];
            theta[i] *= damp[i];
        }
    }

    fclose(fts);
    res.survived = (res.fc_final > 0.5) ? 1 : 0;

    /* DFT of phi_1(0,t) — second half */
    int dft_start = n_dft / 2;
    if (n_dft - dft_start > 100) {
        char dftpath[600];
        snprintf(dftpath, sizeof(dftpath), "%s/mfme_%s_spectrum.tsv", outdir, label);
        FILE *fdft = fopen(dftpath, "w");
        fprintf(fdft, "omega\tpower\n");
        double T = t_hist[n_dft-1] - t_hist[dft_start];
        int nf = 500;
        double peak_pow = 0, peak_om = 0;
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
        res.omega = peak_om;
        printf("\n  Spectrum: peak omega = %.4f (mass gap = %.4f)\n", peak_om, mass);
        printf("  Oscillon (omega < m)? %s\n",
               (peak_om > 0.01 && peak_om < mass) ? "YES" : "NO");
    }

    printf("  fc_final=%.4f  survived=%s  E_phi=%.3f  E_theta=%.3e\n",
           res.fc_final, res.survived ? "YES" : "NO", res.E_phi, res.E_theta);
    printf("  Output: %s\n", tspath);

    for (int a = 0; a < 3; a++) { free(phi[a]); free(vel[a]); free(acc[a]); }
    free(theta); free(vth); free(ath);
    free(RP_cur); free(RP_prev);
    free(damp); free(phi0_hist); free(t_hist);
    return res;
}

/* ======================================================================
 * Phase 1: Baselines
 * ====================================================================== */
static void phase1_baselines(void)
{
    printf("############################################################\n");
    printf("# PHASE 1: BASELINES                                       #\n");
    printf("############################################################\n");

    RunResult r_mf = run_evolution(1.0, 0.0, "baseline_MF_g1_lam0");
    RunResult r_me = run_evolution(0.0, 0.99, "baseline_ME_g0_lam099");

    printf("\n============================================================\n");
    printf("PHASE 1 SUMMARY\n");
    printf("============================================================\n");
    printf("  MF baseline (g=1, lam=0):\n");
    printf("    omega=%.4f  fc=%.4f  survived=%s\n", r_mf.omega, r_mf.fc_final, r_mf.survived ? "YES" : "NO");
    printf("    E_phi=%.3f  E_theta=%.3e  peak_phi=%.4f\n", r_mf.E_phi, r_mf.E_theta, r_mf.peak_phi);
    printf("  ME baseline (g=0, lam=0.99):\n");
    printf("    omega=%.4f  fc=%.4f  survived=%s\n", r_me.omega, r_me.fc_final, r_me.survived ? "YES" : "NO");
    printf("    E_phi=%.3f  E_theta=%.3e  peak_phi=%.4f\n", r_me.E_phi, r_me.E_theta, r_me.peak_phi);
    printf("    m_A=%.6f (predicted range=%.4f)\n",
           (1.0 - 0.99 > 0) ? sqrt(1.0 - 0.99) : 0.0,
           (1.0 - 0.99 > 0) ? 1.0/sqrt(1.0 - 0.99) : 1e30);
    printf("============================================================\n");
}

/* ======================================================================
 * Phase 2: Combined
 * ====================================================================== */
static void phase2_combined(void)
{
    printf("############################################################\n");
    printf("# PHASE 2: COMBINED MF+ME (g=1, lambda=0.99)               #\n");
    printf("############################################################\n");

    RunResult r = run_evolution(1.0, 0.99, "combined_g1_lam099");

    printf("\n============================================================\n");
    printf("PHASE 2 SUMMARY\n");
    printf("============================================================\n");
    printf("  Combined (g=1, lam=0.99):\n");
    printf("    omega=%.4f  fc=%.4f  survived=%s\n", r.omega, r.fc_final, r.survived ? "YES" : "NO");
    printf("    E_phi=%.3f  E_theta=%.3e  peak_phi=%.4f  max_theta=%.4e\n",
           r.E_phi, r.E_theta, r.peak_phi, r.max_theta);
    printf("    m_A=%.6f (predicted range=%.4f)\n", sqrt(0.01), 1.0/sqrt(0.01));
    printf("============================================================\n");
}

/* ======================================================================
 * Phase 3: Push lambda past 1.0 with g=1
 * ====================================================================== */
static void phase3_push(void)
{
    printf("############################################################\n");
    printf("# PHASE 3: PUSH lambda PAST 1.0 WITH g=1                   #\n");
    printf("############################################################\n");

    double lam_vals[] = {0.99, 0.995, 1.0, 1.01, 1.05};
    int n_lam = sizeof(lam_vals) / sizeof(lam_vals[0]);
    RunResult results[5];

    char sumpath[600];
    snprintf(sumpath, sizeof(sumpath), "%s/phase3_summary.tsv", outdir);
    FILE *fsum = fopen(sumpath, "w");
    fprintf(fsum, "lambda\tm2_anti\tm_A\tomega\tfc\tE_phi\tE_theta\tpeak_phi\tsurvived\n");

    for (int k = 0; k < n_lam; k++) {
        double lam = lam_vals[k];
        double m2a = mass * mass - lam;

        char label[128];
        snprintf(label, sizeof(label), "push_g1_lam%.3f", lam);

        /* For tachyonic case (m2_anti < 0), still run — the backreaction might stabilize */
        printf("\n--- lambda=%.4f  m2_anti=%.4f ---\n", lam, m2a);
        results[k] = run_evolution(1.0, lam, label);

        double m_A = (m2a > 0) ? sqrt(m2a) : 0.0;
        fprintf(fsum, "%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4e\t%.4f\t%d\n",
                lam, m2a, m_A, results[k].omega, results[k].fc_final,
                results[k].E_phi, results[k].E_theta, results[k].peak_phi,
                results[k].survived);
        fflush(fsum);
    }
    fclose(fsum);

    /* Also run ME-only (g=0) at the same lambda values for comparison */
    printf("\n############################################################\n");
    printf("# PHASE 3 CONTROL: Same lambdas with g=0 (no backreaction) #\n");
    printf("############################################################\n");

    char sumpath2[600];
    snprintf(sumpath2, sizeof(sumpath2), "%s/phase3_control.tsv", outdir);
    FILE *fsum2 = fopen(sumpath2, "w");
    fprintf(fsum2, "lambda\tm2_anti\tomega\tfc\tE_phi\tsurvived\n");

    for (int k = 0; k < n_lam; k++) {
        double lam = lam_vals[k];
        double m2a = mass * mass - lam;

        char label[128];
        snprintf(label, sizeof(label), "control_g0_lam%.3f", lam);

        printf("\n--- CONTROL: lambda=%.4f  m2_anti=%.4f ---\n", lam, m2a);
        RunResult rc = run_evolution(0.0, lam, label);

        fprintf(fsum2, "%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%d\n",
                lam, m2a, rc.omega, rc.fc_final, rc.E_phi, rc.survived);
        fflush(fsum2);
    }
    fclose(fsum2);

    printf("\n============================================================\n");
    printf("PHASE 3 SUMMARY\n");
    printf("============================================================\n");
    printf("  lambda    m2_anti   omega_g1  fc_g1   surv_g1  omega_g0  fc_g0   surv_g0\n");
    printf("  -------   -------   --------  ------  -------  --------  ------  -------\n");
    /* Note: control results would need separate storage; we print g=1 results */
    for (int k = 0; k < n_lam; k++) {
        double lam = lam_vals[k];
        double m2a = mass * mass - lam;
        printf("  %.4f   %+.4f   %.4f    %.4f  %s\n",
               lam, m2a, results[k].omega, results[k].fc_final,
               results[k].survived ? "YES" : "NO");
    }
    printf("============================================================\n");
    printf("Output: %s\n", sumpath);
}

/* ====================================================================== */
int main(int argc, char **argv)
{
    parse_args(argc, argv);

    switch (run_phase) {
        case 0:
            /* Single run with current g, lambda */
            run_evolution(g_coup, lambda, "single");
            break;
        case 1:
            phase1_baselines();
            break;
        case 2:
            phase2_combined();
            break;
        case 3:
            phase3_push();
            break;
        case 99:
            /* Full run: all phases sequentially */
            phase1_baselines();
            phase2_combined();
            phase3_push();
            break;
        default:
            fprintf(stderr, "Unknown phase %d (use 0-3 or 99 for all)\n", run_phase);
            return 1;
    }
    return 0;
}
