/*
 * proca_field.c — V24-S2: Mediator Field Characterization
 *
 * Decomposes three-field oscillon system into symmetric/antisymmetric modes:
 *   S(x)  = (phi_1 + phi_2 + phi_3) / sqrt(3)      (symmetric)
 *   A1(x) = (phi_1 - phi_2) / sqrt(2)               (antisymmetric mode 1)
 *   A2(x) = (phi_1 + phi_2 - 2*phi_3) / sqrt(6)    (antisymmetric mode 2)
 *
 * Lagrangian:
 *   L = sum_a [ (1/2)(dt phi_a)^2 - (1/2)(dx phi_a)^2 - (m^2/2)phi_a^2 ]
 *     - (mu/2) P^2 / (1 + kappa P^2)           [triple product]
 *     - lambda * (phi_1*phi_2 + phi_2*phi_3 + phi_3*phi_1)  [pairwise]
 *
 * Phase 1: Equilibrate single oscillon for t_equil
 * Phase 2: Place two equilibrated profiles at +-D/2, evolve t_run
 *          Output mode decomposition profiles and midpoint DFT
 * Phase 3: Repeat Phase 1+2 with lambda=0 as control
 *
 * Compile: gcc -O3 -Wall -o proca_field src/proca_field.c -lm
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
static double sig    = 3.0;
static double lambda = 0.99;

/* Grid */
static int    Nx     = 8000;
static double xmax   = 200.0;

/* Timing */
static double t_equil = 5000.0;
static double t_run   = 2000.0;
static double dt_diag = 50.0;

/* Two-oscillon */
static double D_sep  = 30.0;

static char outdir[512] = "v24/proca_field/data";

static void parse_args(int argc, char **argv)
{
    for (int i = 1; i < argc - 1; i += 2) {
        if      (!strcmp(argv[i], "-mu"))     mu     = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-kappa"))  kappa  = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-mass"))   mass   = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-A"))      A_init = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-sigma"))  sig    = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-lambda")) lambda = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-Nx"))     Nx     = atoi(argv[i+1]);
        else if (!strcmp(argv[i], "-xmax"))   xmax   = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-tequil")) t_equil= atof(argv[i+1]);
        else if (!strcmp(argv[i], "-trun"))   t_run  = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-D"))      D_sep  = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-dtdiag")) dt_diag= atof(argv[i+1]);
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

/* ======================================================================
 * Compute timestep from CFL condition
 * ====================================================================== */
static double compute_dt(double dx, double lam)
{
    double m2 = mass * mass;
    double kmax = M_PI / dx;
    double dt_mass = 0.8 * 2.0 / sqrt(kmax * kmax + m2 + 2.0 * fabs(lam));
    double dt_wave = 0.8 * dx;
    return (dt_mass < dt_wave) ? dt_mass : dt_wave;
}

/* ======================================================================
 * Phase 1: Equilibrate a single oscillon
 *   Returns saved phi[3][Nx] and vel[3][Nx] at a breathing peak
 * ====================================================================== */
typedef struct {
    double *phi[3];
    double *vel[3];
    double omega;
    double E_final;
    double peak_amp;
    int    Nx_save;
    double dx_save;
    double xmax_save;
} EquilProfile;

static EquilProfile equilibrate(double lam, const char *label)
{
    EquilProfile prof = {0};

    double dx  = 2.0 * xmax / (Nx - 1);
    double dx2 = dx * dx;
    double m2  = mass * mass;
    double dt  = compute_dt(dx, lam);
    int    Nt  = (int)(t_equil / dt) + 1;

    double m2_anti = m2 - lam;
    double m_A = (m2_anti > 0) ? sqrt(m2_anti) : 0.0;

    printf("\n=== EQUILIBRATE: %s ===\n", label);
    printf("  lambda=%.4f m2_anti=%.4f m_A=%.4f range=%.2f\n",
           lam, m2_anti, m_A, m_A > 0 ? 1.0/m_A : 1e30);
    printf("  Nx=%d xmax=%.1f dx=%.5f dt=%.6f Nt=%d t_equil=%.0f\n",
           Nx, xmax, dx, dt, Nt, t_equil);
    fflush(stdout);

    double *phi[3], *vel[3], *acc[3];
    for (int a = 0; a < 3; a++) {
        phi[a] = calloc(Nx, sizeof(double));
        vel[a] = calloc(Nx, sizeof(double));
        acc[a] = calloc(Nx, sizeof(double));
    }

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

    /* Initialize: Gaussians centered at x=0 */
    int ic = Nx / 2;
    for (int a = 0; a < 3; a++)
        for (int i = 0; i < Nx; i++) {
            double x = -xmax + i * dx;
            phi[a][i] = A_init * exp(-x * x / (2.0 * sig * sig));
        }

    /* DFT storage */
    int max_dft = 50000;
    double *phi0_hist = malloc(max_dft * sizeof(double));
    double *t_hist    = malloc(max_dft * sizeof(double));
    int n_dft = 0;
    int dft_every = Nt / max_dft;
    if (dft_every < 1) dft_every = 1;

    int print_every = Nt / 20;
    if (print_every < 1) print_every = 1;

    /* Track peak for saving profile at breathing maximum */
    double prev_phi0 = 0, prev2_phi0 = 0;
    int save_ready = 0;  /* only save after t > t_equil/2 */
    int saved = 0;

    /* Allocate save buffers */
    for (int a = 0; a < 3; a++) {
        prof.phi[a] = calloc(Nx, sizeof(double));
        prof.vel[a] = calloc(Nx, sizeof(double));
    }

    for (int n = 0; n <= Nt; n++) {
        double t = n * dt;

        if (n % dft_every == 0 && n_dft < max_dft) {
            phi0_hist[n_dft] = phi[0][ic];
            t_hist[n_dft] = t;
            n_dft++;
        }

        if (t > t_equil * 0.5) save_ready = 1;

        /* Detect breathing peak: phi0(ic) was at a local max */
        if (save_ready && n > 2) {
            double cur = phi[0][ic];
            if (prev_phi0 > prev2_phi0 && prev_phi0 > cur && prev_phi0 > 0.1) {
                /* Save at previous step (peak) — actually save current state,
                   close enough for a nearly sinusoidal oscillon */
                for (int a = 0; a < 3; a++) {
                    memcpy(prof.phi[a], phi[a], Nx * sizeof(double));
                    memcpy(prof.vel[a], vel[a], Nx * sizeof(double));
                }
                saved = 1;
            }
        }
        prev2_phi0 = prev_phi0;
        prev_phi0 = phi[0][ic];

        if (n % print_every == 0) {
            double pk = 0;
            for (int i = 0; i < Nx; i++)
                if (fabs(phi[0][i]) > pk) pk = fabs(phi[0][i]);
            printf("  t=%7.0f  phi0(0)=%+.4f  peak=%.4f\n", t, phi[0][ic], pk);
            fflush(stdout);
        }

        if (n == Nt) break;

        /* Velocity Verlet: compute acceleration */
        for (int a = 0; a < 3; a++) {
            acc[a][0] = acc[a][1] = acc[a][Nx-2] = acc[a][Nx-1] = 0;
            int b = (a+1)%3, c = (a+2)%3;
            for (int i = 1; i < Nx - 1; i++) {
                double lapl = (phi[a][i+1] - 2.0*phi[a][i] + phi[a][i-1]) / dx2;
                double fp = force_triple(phi[0][i], phi[1][i], phi[2][i], a);
                acc[a][i] = lapl - m2*phi[a][i] - lam*(phi[b][i]+phi[c][i]) + fp;
            }
        }

        /* Half-kick */
        for (int a = 0; a < 3; a++)
            for (int i = 1; i < Nx - 1; i++)
                vel[a][i] += 0.5 * dt * acc[a][i];

        /* Drift */
        for (int a = 0; a < 3; a++)
            for (int i = 1; i < Nx - 1; i++)
                phi[a][i] += dt * vel[a][i];

        /* Recompute acceleration */
        for (int a = 0; a < 3; a++) {
            acc[a][0] = acc[a][1] = acc[a][Nx-2] = acc[a][Nx-1] = 0;
            int b = (a+1)%3, c = (a+2)%3;
            for (int i = 1; i < Nx - 1; i++) {
                double lapl = (phi[a][i+1] - 2.0*phi[a][i] + phi[a][i-1]) / dx2;
                double fp = force_triple(phi[0][i], phi[1][i], phi[2][i], a);
                acc[a][i] = lapl - m2*phi[a][i] - lam*(phi[b][i]+phi[c][i]) + fp;
            }
        }

        /* Second half-kick */
        for (int a = 0; a < 3; a++)
            for (int i = 1; i < Nx - 1; i++)
                vel[a][i] += 0.5 * dt * acc[a][i];

        /* Absorbing boundary */
        for (int a = 0; a < 3; a++)
            for (int i = 0; i < Nx; i++) {
                vel[a][i] *= damp[i];
                phi[a][i] *= damp[i];
            }
    }

    /* If no peak was found, save final state */
    if (!saved) {
        printf("  WARNING: no breathing peak detected, saving final state\n");
        for (int a = 0; a < 3; a++) {
            memcpy(prof.phi[a], phi[a], Nx * sizeof(double));
            memcpy(prof.vel[a], vel[a], Nx * sizeof(double));
        }
    }

    /* Compute peak amplitude of saved profile */
    prof.peak_amp = 0;
    for (int i = 0; i < Nx; i++)
        if (fabs(prof.phi[0][i]) > prof.peak_amp)
            prof.peak_amp = fabs(prof.phi[0][i]);

    /* DFT for frequency */
    int dft_start = n_dft / 2;
    prof.omega = 0;
    if (n_dft - dft_start > 100) {
        double T = t_hist[n_dft-1] - t_hist[dft_start];
        int nf = 500;
        double peak_pow = 0;
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
            if (pw > peak_pow) { peak_pow = pw; prof.omega = omega; }
        }
        printf("  Equilibrated: omega=%.4f peak_amp=%.4f\n", prof.omega, prof.peak_amp);
    }

    prof.Nx_save   = Nx;
    prof.dx_save   = dx;
    prof.xmax_save = xmax;

    /* Compute final energy */
    prof.E_final = 0;
    for (int i = 1; i < Nx - 1; i++) {
        for (int a = 0; a < 3; a++) {
            prof.E_final += 0.5 * vel[a][i] * vel[a][i] * dx;
            double dp = (phi[a][i+1] - phi[a][i-1]) / (2.0*dx);
            prof.E_final += 0.5 * dp * dp * dx;
            prof.E_final += 0.5 * m2 * phi[a][i] * phi[a][i] * dx;
        }
        double P = phi[0][i] * phi[1][i] * phi[2][i];
        double P2 = P * P;
        prof.E_final += 0.5 * mu * P2 / (1.0 + kappa * P2) * dx;
        prof.E_final += lam * (phi[0][i]*phi[1][i] + phi[1][i]*phi[2][i]
                              + phi[2][i]*phi[0][i]) * dx;
    }
    printf("  E_final=%.4f\n", prof.E_final);

    for (int a = 0; a < 3; a++) { free(phi[a]); free(vel[a]); free(acc[a]); }
    free(damp); free(phi0_hist); free(t_hist);

    return prof;
}

/* ======================================================================
 * Phase 2: Two-oscillon evolution with mode decomposition
 * ====================================================================== */
typedef struct {
    double S_rms_mid;     /* RMS of S in midpoint region */
    double A1_rms_mid;    /* RMS of A1 in midpoint region */
    double A2_rms_mid;    /* RMS of A2 in midpoint region */
    double omega_S;       /* peak frequency of S at midpoint */
    double omega_A1;      /* peak frequency of A1 at midpoint */
    double omega_A2;      /* peak frequency of A2 at midpoint */
} TwoBodyResult;

static TwoBodyResult two_body(EquilProfile *prof, double lam, const char *label)
{
    TwoBodyResult res = {0};

    double dx  = prof->dx_save;
    double dx2 = dx * dx;
    double m2  = mass * mass;
    double dt  = compute_dt(dx, lam);
    int    Nt  = (int)(t_run / dt) + 1;

    double m2_anti = m2 - lam;
    double m_A = (m2_anti > 0) ? sqrt(m2_anti) : 0.0;

    printf("\n=== TWO-BODY: %s ===\n", label);
    printf("  lambda=%.4f D=%.1f m_A=%.4f range=%.2f\n",
           lam, D_sep, m_A, m_A > 0 ? 1.0/m_A : 1e30);
    printf("  Nx=%d xmax=%.1f dx=%.5f dt=%.6f Nt=%d t_run=%.0f\n",
           Nx, xmax, dx, dt, Nt, t_run);
    fflush(stdout);

    double *phi[3], *vel[3], *acc[3];
    for (int a = 0; a < 3; a++) {
        phi[a] = calloc(Nx, sizeof(double));
        vel[a] = calloc(Nx, sizeof(double));
        acc[a] = calloc(Nx, sizeof(double));
    }

    /* Absorbing boundary */
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

    /* Initialize: place two copies of the equilibrated profile at +-D/2 */
    int ic = Nx / 2;
    int shift = (int)(D_sep / (2.0 * dx) + 0.5);  /* grid points for D/2 */
    printf("  Center index=%d, shift=%d pts (D/2=%.2f)\n", ic, shift, shift*dx);

    for (int a = 0; a < 3; a++) {
        /* Profile is centered at ic in the saved arrays */
        for (int i = 0; i < Nx; i++) {
            double val = 0, vval = 0;

            /* Left oscillon: shift profile to i - shift relative to center */
            int jL = i - ic + shift + (prof->Nx_save / 2);
            if (jL >= 0 && jL < prof->Nx_save) {
                val  += prof->phi[a][jL];
                vval += prof->vel[a][jL];
            }

            /* Right oscillon: shift profile to i + shift relative to center */
            int jR = i - ic - shift + (prof->Nx_save / 2);
            if (jR >= 0 && jR < prof->Nx_save) {
                val  += prof->phi[a][jR];
                vval += prof->vel[a][jR];
            }

            phi[a][i] = val;
            vel[a][i] = vval;
        }
    }

    /* Add antisymmetric perturbation to seed A1 mode:
     * A1 = (phi_1 - phi_2)/sqrt(2), so add +eps to phi_1, -eps to phi_2
     * Gaussian centered at midpoint (x=0) with width sigma_pert */
    double eps_pert = 0.01;   /* small perturbation amplitude */
    double sig_pert = 3.0;    /* width of perturbation */
    printf("  Antisymmetric seed: eps=%.4f sigma=%.2f\n", eps_pert, sig_pert);
    for (int i = 0; i < Nx; i++) {
        double x = -xmax + i * dx;
        double g = eps_pert * exp(-x * x / (2.0 * sig_pert * sig_pert));
        phi[0][i] += g;    /* +eps in phi_1 */
        phi[1][i] -= g;    /* -eps in phi_2 */
        /* phi_3 unchanged: this seeds A1 = sqrt(2)*eps*Gaussian, A2 = 0 */
    }

    /* Midpoint DFT storage */
    int max_dft = (int)(t_run / dt) + 1;
    if (max_dft > 100000) max_dft = 100000;
    double *S_mid_hist  = malloc(max_dft * sizeof(double));
    double *A1_mid_hist = malloc(max_dft * sizeof(double));
    double *A2_mid_hist = malloc(max_dft * sizeof(double));
    double *t_hist      = malloc(max_dft * sizeof(double));
    int n_dft = 0;
    int dft_every = Nt / max_dft;
    if (dft_every < 1) dft_every = 1;

    /* Diagnostic interval */
    int diag_every = (int)(dt_diag / dt + 0.5);
    if (diag_every < 1) diag_every = 1;

    int print_every = Nt / 40;
    if (print_every < 1) print_every = 1;

    /* Profile output times: t=0, 100, 500, 1000, 2000 */
    double prof_times[] = {0.0, 100.0, 500.0, 1000.0, 2000.0};
    int n_prof_times = sizeof(prof_times) / sizeof(prof_times[0]);
    int prof_written[5] = {0};

    /* Time series file */
    char tspath[600];
    snprintf(tspath, sizeof(tspath), "%s/twobody_%s_ts.tsv", outdir, label);
    FILE *fts = fopen(tspath, "w");
    if (!fts) { fprintf(stderr, "Cannot open %s\n", tspath); return res; }
    fprintf(fts, "time\tS_mid\tA1_mid\tA2_mid\t"
                 "S_rms_between\tA1_rms_between\tA2_rms_between\t"
                 "E_total\tpeak_phi\n");

    /* Midpoint region: D/4 < |x| < 3D/4 from center
       Actually, midpoint between oscillons at +-D/2 is x=0.
       "Between" region: |x| < D/2 - core_width, say |x| < D/4 */
    int i_mid = ic;   /* midpoint = center of grid */
    double x_between_lo = -D_sep * 0.375;   /* inner 75% of gap */
    double x_between_hi =  D_sep * 0.375;

    /* RMS accumulators for final result */
    double S_rms_acc = 0, A1_rms_acc = 0, A2_rms_acc = 0;
    int n_rms = 0;

    for (int n = 0; n <= Nt; n++) {
        double t = n * dt;

        /* Mode decomposition at midpoint */
        double S_m  = (phi[0][i_mid] + phi[1][i_mid] + phi[2][i_mid]) / sqrt(3.0);
        double A1_m = (phi[0][i_mid] - phi[1][i_mid]) / sqrt(2.0);
        double A2_m = (phi[0][i_mid] + phi[1][i_mid] - 2.0*phi[2][i_mid]) / sqrt(6.0);

        /* DFT sampling */
        if (n % dft_every == 0 && n_dft < max_dft) {
            S_mid_hist[n_dft]  = S_m;
            A1_mid_hist[n_dft] = A1_m;
            A2_mid_hist[n_dft] = A2_m;
            t_hist[n_dft] = t;
            n_dft++;
        }

        int do_diag  = (n % diag_every == 0);
        int do_print = (n % print_every == 0);

        if (do_diag || do_print) {
            /* Compute RMS of S, A1, A2 in the between-oscillons region */
            double S_sq = 0, A1_sq = 0, A2_sq = 0;
            int n_between = 0;
            double pk_phi = 0;
            double Etot = 0;

            for (int i = 1; i < Nx - 1; i++) {
                double x = -xmax + i * dx;

                /* Peak amplitude */
                for (int a = 0; a < 3; a++)
                    if (fabs(phi[a][i]) > pk_phi) pk_phi = fabs(phi[a][i]);

                /* Energy density */
                double e = 0;
                for (int a = 0; a < 3; a++) {
                    e += 0.5 * vel[a][i] * vel[a][i];
                    double dp = (phi[a][i+1] - phi[a][i-1]) / (2.0*dx);
                    e += 0.5 * dp * dp;
                    e += 0.5 * m2 * phi[a][i] * phi[a][i];
                }
                double P = phi[0][i] * phi[1][i] * phi[2][i];
                double P2 = P * P;
                e += 0.5 * mu * P2 / (1.0 + kappa * P2);
                e += lam * (phi[0][i]*phi[1][i] + phi[1][i]*phi[2][i]
                           + phi[2][i]*phi[0][i]);
                Etot += e * dx;

                /* Between-region RMS */
                if (x >= x_between_lo && x <= x_between_hi) {
                    double S  = (phi[0][i] + phi[1][i] + phi[2][i]) / sqrt(3.0);
                    double A1 = (phi[0][i] - phi[1][i]) / sqrt(2.0);
                    double A2 = (phi[0][i] + phi[1][i] - 2.0*phi[2][i]) / sqrt(6.0);
                    S_sq  += S * S;
                    A1_sq += A1 * A1;
                    A2_sq += A2 * A2;
                    n_between++;
                }
            }

            double S_rms  = (n_between > 0) ? sqrt(S_sq / n_between) : 0;
            double A1_rms = (n_between > 0) ? sqrt(A1_sq / n_between) : 0;
            double A2_rms = (n_between > 0) ? sqrt(A2_sq / n_between) : 0;

            /* Accumulate for time-averaged RMS (second half) */
            if (t > t_run * 0.5) {
                S_rms_acc  += S_rms;
                A1_rms_acc += A1_rms;
                A2_rms_acc += A2_rms;
                n_rms++;
            }

            if (do_diag)
                fprintf(fts, "%.6f\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\n",
                        t, S_m, A1_m, A2_m, S_rms, A1_rms, A2_rms, Etot, pk_phi);

            if (do_print)
                printf("  t=%7.0f  S_m=%+.4e A1_m=%+.4e A2_m=%+.4e  "
                       "S_rms=%.3e A1_rms=%.3e  E=%.3f\n",
                       t, S_m, A1_m, A2_m, S_rms, A1_rms, Etot);
        }

        /* Output spatial profiles at selected times */
        for (int pt = 0; pt < n_prof_times; pt++) {
            if (!prof_written[pt] && t >= prof_times[pt] - 0.5*dt) {
                char ppath[600];
                snprintf(ppath, sizeof(ppath), "%s/field_decomp_%s_t%04d.tsv",
                         outdir, label, (int)prof_times[pt]);
                FILE *fp = fopen(ppath, "w");
                if (fp) {
                    fprintf(fp, "x\tphi1\tphi2\tphi3\tS\tA1\tA2\n");
                    /* Output every 4th point to keep file size manageable */
                    for (int i = 0; i < Nx; i += 4) {
                        double x = -xmax + i * dx;
                        double S  = (phi[0][i] + phi[1][i] + phi[2][i]) / sqrt(3.0);
                        double A1 = (phi[0][i] - phi[1][i]) / sqrt(2.0);
                        double A2 = (phi[0][i] + phi[1][i] - 2.0*phi[2][i]) / sqrt(6.0);
                        fprintf(fp, "%.4f\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\n",
                                x, phi[0][i], phi[1][i], phi[2][i], S, A1, A2);
                    }
                    fclose(fp);
                    printf("  Wrote profile: %s\n", ppath);
                }
                prof_written[pt] = 1;
            }
        }

        if (n == Nt) break;

        /* Velocity Verlet */
        /* Compute acceleration */
        for (int a = 0; a < 3; a++) {
            acc[a][0] = acc[a][1] = acc[a][Nx-2] = acc[a][Nx-1] = 0;
            int b = (a+1)%3, c = (a+2)%3;
            for (int i = 1; i < Nx - 1; i++) {
                double lapl = (phi[a][i+1] - 2.0*phi[a][i] + phi[a][i-1]) / dx2;
                double fp = force_triple(phi[0][i], phi[1][i], phi[2][i], a);
                acc[a][i] = lapl - m2*phi[a][i] - lam*(phi[b][i]+phi[c][i]) + fp;
            }
        }

        /* Half-kick */
        for (int a = 0; a < 3; a++)
            for (int i = 1; i < Nx - 1; i++)
                vel[a][i] += 0.5 * dt * acc[a][i];

        /* Drift */
        for (int a = 0; a < 3; a++)
            for (int i = 1; i < Nx - 1; i++)
                phi[a][i] += dt * vel[a][i];

        /* Recompute acceleration */
        for (int a = 0; a < 3; a++) {
            acc[a][0] = acc[a][1] = acc[a][Nx-2] = acc[a][Nx-1] = 0;
            int b = (a+1)%3, c = (a+2)%3;
            for (int i = 1; i < Nx - 1; i++) {
                double lapl = (phi[a][i+1] - 2.0*phi[a][i] + phi[a][i-1]) / dx2;
                double fp = force_triple(phi[0][i], phi[1][i], phi[2][i], a);
                acc[a][i] = lapl - m2*phi[a][i] - lam*(phi[b][i]+phi[c][i]) + fp;
            }
        }

        /* Second half-kick */
        for (int a = 0; a < 3; a++)
            for (int i = 1; i < Nx - 1; i++)
                vel[a][i] += 0.5 * dt * acc[a][i];

        /* Absorbing boundary */
        for (int a = 0; a < 3; a++)
            for (int i = 0; i < Nx; i++) {
                vel[a][i] *= damp[i];
                phi[a][i] *= damp[i];
            }
    }

    fclose(fts);

    /* Time-averaged RMS */
    if (n_rms > 0) {
        res.S_rms_mid  = S_rms_acc / n_rms;
        res.A1_rms_mid = A1_rms_acc / n_rms;
        res.A2_rms_mid = A2_rms_acc / n_rms;
    }

    /* DFT of midpoint time series — second half only */
    int dft_start = n_dft / 2;
    if (n_dft - dft_start > 100) {
        char dftpath[600];
        snprintf(dftpath, sizeof(dftpath), "%s/midpoint_spectrum_%s.tsv", outdir, label);
        FILE *fdft = fopen(dftpath, "w");
        fprintf(fdft, "omega\tpower_S\tpower_A1\tpower_A2\n");

        double T = t_hist[n_dft-1] - t_hist[dft_start];
        int nf = 500;
        double pk_S = 0, pk_A1 = 0, pk_A2 = 0;
        double om_S = 0, om_A1 = 0, om_A2 = 0;

        for (int k = 0; k < nf; k++) {
            double omega = 3.0 * mass * k / nf;
            double re_S = 0, im_S = 0;
            double re_A1 = 0, im_A1 = 0;
            double re_A2 = 0, im_A2 = 0;

            for (int j = dft_start; j < n_dft; j++) {
                double dtj = (j > dft_start) ?
                    (t_hist[j]-t_hist[j-1]) : (t_hist[dft_start+1]-t_hist[dft_start]);
                double cs = cos(omega * t_hist[j]) * dtj;
                double sn = sin(omega * t_hist[j]) * dtj;
                re_S  += S_mid_hist[j]  * cs;  im_S  += S_mid_hist[j]  * sn;
                re_A1 += A1_mid_hist[j] * cs;  im_A1 += A1_mid_hist[j] * sn;
                re_A2 += A2_mid_hist[j] * cs;  im_A2 += A2_mid_hist[j] * sn;
            }
            double pw_S  = (re_S*re_S   + im_S*im_S)   / (T*T);
            double pw_A1 = (re_A1*re_A1 + im_A1*im_A1) / (T*T);
            double pw_A2 = (re_A2*re_A2 + im_A2*im_A2) / (T*T);

            fprintf(fdft, "%.6f\t%.6e\t%.6e\t%.6e\n", omega, pw_S, pw_A1, pw_A2);

            if (pw_S  > pk_S)  { pk_S  = pw_S;  om_S  = omega; }
            if (pw_A1 > pk_A1) { pk_A1 = pw_A1; om_A1 = omega; }
            if (pw_A2 > pk_A2) { pk_A2 = pw_A2; om_A2 = omega; }
        }
        fclose(fdft);

        res.omega_S  = om_S;
        res.omega_A1 = om_A1;
        res.omega_A2 = om_A2;

        printf("\n  Midpoint DFT:\n");
        printf("    S:  peak omega=%.4f  power=%.4e\n", om_S, pk_S);
        printf("    A1: peak omega=%.4f  power=%.4e\n", om_A1, pk_A1);
        printf("    A2: peak omega=%.4f  power=%.4e\n", om_A2, pk_A2);
        printf("    mass gap m=%.4f, m_A=%.4f\n", mass, m_A);
    }

    printf("\n  Time-averaged RMS between oscillons (t > %.0f):\n", t_run * 0.5);
    printf("    S_rms =%.6e\n", res.S_rms_mid);
    printf("    A1_rms=%.6e\n", res.A1_rms_mid);
    printf("    A2_rms=%.6e\n", res.A2_rms_mid);

    for (int a = 0; a < 3; a++) { free(phi[a]); free(vel[a]); free(acc[a]); }
    free(damp);
    free(S_mid_hist); free(A1_mid_hist); free(A2_mid_hist); free(t_hist);

    return res;
}

/* ====================================================================== */
int main(int argc, char **argv)
{
    parse_args(argc, argv);

    printf("############################################################\n");
    printf("# V24-S2: Mediator Field Characterization                  #\n");
    printf("############################################################\n");
    printf("  mu=%.1f kappa=%.1f mass=%.4f A=%.3f sigma=%.3f\n",
           mu, kappa, mass, A_init, sig);
    printf("  lambda=%.4f D=%.1f\n", lambda, D_sep);
    printf("  Nx=%d xmax=%.1f t_equil=%.0f t_run=%.0f dt_diag=%.1f\n",
           Nx, xmax, t_equil, t_run, dt_diag);

    /* ============================================================
     * Step 1: Equilibrate at lambda=0.99
     * ============================================================ */
    printf("\n############################################################\n");
    printf("# STEP 1: Equilibrate at lambda=%.4f                      #\n", lambda);
    printf("############################################################\n");
    EquilProfile prof_lam = equilibrate(lambda, "lam099");

    /* ============================================================
     * Step 2: Two-oscillon with lambda=0.99
     * ============================================================ */
    printf("\n############################################################\n");
    printf("# STEP 2: Two-body at lambda=%.4f, D=%.1f                 #\n", lambda, D_sep);
    printf("############################################################\n");
    TwoBodyResult res_lam = two_body(&prof_lam, lambda, "lam099");

    /* ============================================================
     * Step 3: Control — equilibrate at lambda=0
     * ============================================================ */
    printf("\n############################################################\n");
    printf("# STEP 3: Control equilibration at lambda=0                #\n");
    printf("############################################################\n");
    EquilProfile prof_0 = equilibrate(0.0, "lam000");

    /* ============================================================
     * Step 4: Two-oscillon control with lambda=0
     * ============================================================ */
    printf("\n############################################################\n");
    printf("# STEP 4: Two-body control at lambda=0, D=%.1f            #\n", D_sep);
    printf("############################################################\n");
    TwoBodyResult res_0 = two_body(&prof_0, 0.0, "lam000");

    /* ============================================================
     * Summary
     * ============================================================ */
    printf("\n############################################################\n");
    printf("# SUMMARY                                                  #\n");
    printf("############################################################\n");
    printf("\n  Equilibrated profiles:\n");
    printf("    lambda=%.4f: omega=%.4f  peak=%.4f  E=%.4f\n",
           lambda, prof_lam.omega, prof_lam.peak_amp, prof_lam.E_final);
    printf("    lambda=0.00:  omega=%.4f  peak=%.4f  E=%.4f\n",
           prof_0.omega, prof_0.peak_amp, prof_0.E_final);

    printf("\n  Two-body RMS between oscillons (time-averaged, second half):\n");
    printf("    lambda=%.4f:  S_rms=%.6e  A1_rms=%.6e  A2_rms=%.6e\n",
           lambda, res_lam.S_rms_mid, res_lam.A1_rms_mid, res_lam.A2_rms_mid);
    printf("    lambda=0.00:   S_rms=%.6e  A1_rms=%.6e  A2_rms=%.6e\n",
           res_0.S_rms_mid, res_0.A1_rms_mid, res_0.A2_rms_mid);

    double ratio_S  = (res_0.S_rms_mid  > 1e-30) ? res_lam.S_rms_mid  / res_0.S_rms_mid  : 0;
    double ratio_A1 = (res_0.A1_rms_mid > 1e-30) ? res_lam.A1_rms_mid / res_0.A1_rms_mid : 0;
    double ratio_A2 = (res_0.A2_rms_mid > 1e-30) ? res_lam.A2_rms_mid / res_0.A2_rms_mid : 0;

    printf("\n  Ratios (lambda=%.4f / lambda=0):\n", lambda);
    printf("    S:  %.4f\n", ratio_S);
    printf("    A1: %.4f\n", ratio_A1);
    printf("    A2: %.4f\n", ratio_A2);

    printf("\n  Midpoint DFT peaks:\n");
    printf("    lambda=%.4f:  omega_S=%.4f  omega_A1=%.4f  omega_A2=%.4f\n",
           lambda, res_lam.omega_S, res_lam.omega_A1, res_lam.omega_A2);
    printf("    lambda=0.00:   omega_S=%.4f  omega_A1=%.4f  omega_A2=%.4f\n",
           res_0.omega_S, res_0.omega_A1, res_0.omega_A2);

    double m_A_lam = sqrt(fabs(mass*mass - lambda));
    printf("\n  Mass scales:\n");
    printf("    Symmetric mass m = %.4f (range = %.2f)\n", mass, 1.0/mass);
    printf("    Antisymmetric mass m_A = %.4f (range = %.2f)\n", m_A_lam, 1.0/m_A_lam);

    int A_nonzero = (res_lam.A1_rms_mid > 1e-10 || res_lam.A2_rms_mid > 1e-10);
    int A_enhanced = (ratio_A1 > 2.0 || ratio_A2 > 2.0);
    printf("\n  KEY QUESTIONS:\n");
    printf("    Is antisymmetric mode nonzero between oscillons at lambda=%.4f? %s\n",
           lambda, A_nonzero ? "YES" : "NO");
    printf("    Is it larger than at lambda=0? %s (ratio A1=%.2f, A2=%.2f)\n",
           A_enhanced ? "YES" : "NO", ratio_A1, ratio_A2);

    /* Write summary table */
    char sumpath[600];
    snprintf(sumpath, sizeof(sumpath), "%s/summary.tsv", outdir);
    FILE *fsum = fopen(sumpath, "w");
    if (fsum) {
        fprintf(fsum, "lambda\tS_rms_mid\tA1_rms_mid\tA2_rms_mid\t"
                      "omega_S\tomega_A1\tomega_A2\tomega_equil\tpeak_equil\tE_equil\n");
        fprintf(fsum, "%.4f\t%.6e\t%.6e\t%.6e\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n",
                lambda, res_lam.S_rms_mid, res_lam.A1_rms_mid, res_lam.A2_rms_mid,
                res_lam.omega_S, res_lam.omega_A1, res_lam.omega_A2,
                prof_lam.omega, prof_lam.peak_amp, prof_lam.E_final);
        fprintf(fsum, "0.0000\t%.6e\t%.6e\t%.6e\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n",
                res_0.S_rms_mid, res_0.A1_rms_mid, res_0.A2_rms_mid,
                res_0.omega_S, res_0.omega_A1, res_0.omega_A2,
                prof_0.omega, prof_0.peak_amp, prof_0.E_final);
        fclose(fsum);
        printf("\n  Output: %s\n", sumpath);
    }

    printf("\n############################################################\n");
    printf("# DONE                                                     #\n");
    printf("############################################################\n");

    /* Cleanup */
    for (int a = 0; a < 3; a++) {
        free(prof_lam.phi[a]); free(prof_lam.vel[a]);
        free(prof_0.phi[a]);   free(prof_0.vel[a]);
    }

    return 0;
}
