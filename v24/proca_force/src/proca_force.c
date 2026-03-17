/*
 * proca_force.c -- V24-S1: Clean Force Measurement with Pairwise Coupling
 *
 * Three massive scalars with triple-product AND pairwise coupling:
 *   L = sum_a [ (1/2)(dt phi_a)^2 - (1/2)(dx phi_a)^2 - (m^2/2)phi_a^2 ]
 *     - (mu/2) P^2 / (1 + kappa P^2)
 *     - lambda(phi_1 phi_2 + phi_2 phi_3 + phi_3 phi_1)
 *
 * Mass spectrum:
 *   antisymmetric mode: m^2_A = m^2 - lambda  (lighter, longer range)
 *   symmetric mode:     m^2_S = m^2 + 2*lambda (heavier)
 *
 * Protocol:
 *   Step 1: Equilibrate single oscillon WITH pairwise coupling for t=10000
 *           Save state at breathing peak. Repeat for each lambda.
 *   Step 2: Place two equilibrated profiles at +/- D/2. Evolve t=2000.
 *           Track sep(t). Measure initial force from quadratic fit t in [0,500].
 *   Step 3: Fit F(D) to Yukawa. Compare range with 1/m_A.
 *
 * Compile: gcc -O3 -Wall -o proca_force src/proca_force.c -lm
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* --- Parameters --- */
static double mu       = -20.0;
static double kappa    = 20.0;
static double mass     = 1.0;
static double A_init   = 0.8;
static double sig      = 3.0;
static char   outdir[512] = "data";

/* Equilibration grid */
static int    Nx_eq    = 4000;
static double xmax_eq  = 80.0;
static double t_equil  = 10000.0;

/* Interaction grid */
static int    Nx_int   = 8000;
static double xmax_int = 300.0;
static double t_run    = 2000.0;

/* Scans */
static double lam_list[] = {0.0, 0.5, 0.9, 0.99, 0.999};
static int    n_lam      = 5;
static double D_list[]   = {15.0, 20.0, 30.0, 40.0, 60.0, 80.0, 100.0};
static int    n_D        = 7;

/* -dV_triple/dphi_a where V = (mu/2)P^2/(1+kappa P^2), P = phi1 phi2 phi3 */
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

/* ===================================================================
 *  Saved profile data
 * =================================================================== */
typedef struct {
    int    N;
    double xmax;
    double dx;
    double omega;
    double A_peak;
    double E_total;
    double *phi[3];
    double *vel[3];
} profile_t;

static void profile_free(profile_t *p)
{
    for (int a = 0; a < 3; a++) {
        if (p->phi[a]) { free(p->phi[a]); p->phi[a] = NULL; }
        if (p->vel[a]) { free(p->vel[a]); p->vel[a] = NULL; }
    }
}

static int profile_save(const char *path, const profile_t *p)
{
    FILE *f = fopen(path, "w");
    if (!f) { fprintf(stderr, "Cannot open %s for writing\n", path); return -1; }
    fprintf(f, "# N=%d xmax=%.6f dx=%.6e omega=%.6f A_peak=%.6f E_total=%.6f\n",
            p->N, p->xmax, p->dx, p->omega, p->A_peak, p->E_total);
    fprintf(f, "# x  phi1  phi2  phi3  vel1  vel2  vel3\n");
    for (int i = 0; i < p->N; i++) {
        double x = -p->xmax + i * p->dx;
        fprintf(f, "%.8e  %.12e  %.12e  %.12e  %.12e  %.12e  %.12e\n",
                x, p->phi[0][i], p->phi[1][i], p->phi[2][i],
                p->vel[0][i], p->vel[1][i], p->vel[2][i]);
    }
    fclose(f);
    return 0;
}

/* ===================================================================
 *  Step 1: Equilibrate a single oscillon WITH pairwise coupling
 * =================================================================== */
static profile_t run_equilibrate(double lam)
{
    profile_t prof = {0};
    prof.N    = Nx_eq;
    prof.xmax = xmax_eq;
    prof.dx   = 2.0 * xmax_eq / (Nx_eq - 1);

    double dx  = prof.dx;
    double dx2 = dx * dx;
    double m2  = mass * mass;
    int    Nx  = Nx_eq;
    int    ic  = Nx / 2;

    double m2_anti = m2 - lam;
    double m2_sym  = m2 + 2.0 * lam;
    if (m2_anti <= 0) {
        printf("  WARNING: tachyonic at lambda=%.4f, m2_anti=%.4f\n", lam, m2_anti);
    }

    /* CFL based on heaviest mode */
    double kmax = M_PI / dx;
    double dt   = 0.8 * 2.0 / sqrt(kmax * kmax + (m2_sym > m2 ? m2_sym : m2));
    int    Nt   = (int)(t_equil / dt) + 1;

    printf("  Equilibrate: lam=%.4f Nx=%d xmax=%.0f dx=%.5f dt=%.6f Nt=%d\n",
           lam, Nx, xmax_eq, dx, dt, Nt);
    printf("    m2_anti=%.4f m_A=%.4f  m2_sym=%.4f m_S=%.4f\n",
           m2_anti, m2_anti > 0 ? sqrt(m2_anti) : 0.0, m2_sym, sqrt(m2_sym));

    /* Allocate fields */
    double *phi[3], *vel[3], *acc[3];
    for (int a = 0; a < 3; a++) {
        phi[a] = calloc(Nx, sizeof(double));
        vel[a] = calloc(Nx, sizeof(double));
        acc[a] = calloc(Nx, sizeof(double));
    }

    /* Absorbing boundary: outer 25% */
    double *damp = malloc(Nx * sizeof(double));
    double x_abs = xmax_eq * 0.75;
    for (int i = 0; i < Nx; i++) {
        double x = -xmax_eq + i * dx;
        double ax = fabs(x);
        if (ax > x_abs) {
            double f = (ax - x_abs) / (xmax_eq - x_abs);
            damp[i] = 1.0 - 0.98 * f * f;
        } else {
            damp[i] = 1.0;
        }
    }

    /* Initialize: symmetric Gaussians */
    for (int a = 0; a < 3; a++)
        for (int i = 0; i < Nx; i++) {
            double x = -xmax_eq + i * dx;
            phi[a][i] = A_init * exp(-x * x / (2.0 * sig * sig));
        }

    /* Compute acceleration with pairwise coupling */
    #define COMPUTE_ACC_EQ() do { \
        for (int aa = 0; aa < 3; aa++) { \
            acc[aa][0] = acc[aa][1] = acc[aa][Nx-2] = acc[aa][Nx-1] = 0; \
            int bb = (aa+1)%3, cc = (aa+2)%3; \
            for (int ii = 1; ii < Nx - 1; ii++) { \
                double lapl = (phi[aa][ii+1] - 2.0*phi[aa][ii] + phi[aa][ii-1]) / dx2; \
                double fp = force_triple(phi[0][ii], phi[1][ii], phi[2][ii], aa); \
                acc[aa][ii] = lapl - m2*phi[aa][ii] - lam*(phi[bb][ii]+phi[cc][ii]) + fp; \
            } \
        } \
    } while(0)

    COMPUTE_ACC_EQ();

    /* DFT storage for frequency measurement */
    int max_dft = 100000;
    double *phi0_hist = malloc(max_dft * sizeof(double));
    double *t_hist    = malloc(max_dft * sizeof(double));
    int n_dft = 0;
    int dft_every = Nt / max_dft;
    if (dft_every < 1) dft_every = 1;

    int print_every = Nt / 10;
    if (print_every < 1) print_every = 1;

    /* Snapshot storage for peak detection */
    double *snap_phi[3], *snap_vel[3];
    for (int a = 0; a < 3; a++) {
        snap_phi[a] = calloc(Nx, sizeof(double));
        snap_vel[a] = calloc(Nx, sizeof(double));
    }
    int have_snap = 0;
    double prev_phi0 = 0, prev_prev_phi0 = 0;
    double snap_phi0_val = 0;

    for (int n = 0; n <= Nt; n++) {
        double t = n * dt;

        /* Record for DFT */
        if (n % dft_every == 0 && n_dft < max_dft) {
            phi0_hist[n_dft] = phi[0][ic];
            t_hist[n_dft] = t;
            n_dft++;
        }

        /* Detect positive peaks in last 5% and save full state */
        if (t > t_equil * 0.95) {
            double cur = phi[0][ic];
            if (prev_phi0 > prev_prev_phi0 && prev_phi0 > cur && prev_phi0 > 0.01) {
                for (int a = 0; a < 3; a++) {
                    memcpy(snap_phi[a], phi[a], Nx * sizeof(double));
                    memcpy(snap_vel[a], vel[a], Nx * sizeof(double));
                }
                have_snap = 1;
                snap_phi0_val = prev_phi0;
            }
            prev_prev_phi0 = prev_phi0;
            prev_phi0 = cur;
        }

        /* Print progress */
        if (n % print_every == 0) {
            double Eall = 0, Ecore = 0, peak = 0;
            double core_r = 3.0 * sig;
            for (int i = 1; i < Nx - 1; i++) {
                double x = -xmax_eq + i * dx;
                double e = 0;
                for (int a = 0; a < 3; a++) {
                    e += 0.5 * vel[a][i] * vel[a][i];
                    double dp = (phi[a][i+1] - phi[a][i-1]) / (2.0*dx);
                    e += 0.5 * dp * dp + 0.5 * m2 * phi[a][i] * phi[a][i];
                    if (fabs(phi[a][i]) > peak) peak = fabs(phi[a][i]);
                }
                e += lam * (phi[0][i]*phi[1][i] + phi[1][i]*phi[2][i]
                           + phi[2][i]*phi[0][i]);
                double P = phi[0][i] * phi[1][i] * phi[2][i];
                e += 0.5 * mu * P * P / (1.0 + kappa * P * P);
                Eall += e * dx;
                if (fabs(x) < core_r) Ecore += e * dx;
            }
            double fc = (Eall > 1e-20) ? Ecore / Eall : 0.0;
            printf("    t=%7.0f  phi0=%+.4f  pk=%.4f  E=%+.4f  fc=%.4f\n",
                   t, phi[0][ic], peak, Eall, fc);
        }

        if (n == Nt) break;

        /* Velocity Verlet */
        for (int a = 0; a < 3; a++)
            for (int i = 1; i < Nx - 1; i++)
                vel[a][i] += 0.5 * dt * acc[a][i];
        for (int a = 0; a < 3; a++)
            for (int i = 1; i < Nx - 1; i++)
                phi[a][i] += dt * vel[a][i];
        COMPUTE_ACC_EQ();
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
    #undef COMPUTE_ACC_EQ

    /* Use saved snapshot */
    if (have_snap) {
        printf("  Using saved peak snapshot: phi1(0)=%.6f\n", snap_phi0_val);
        for (int a = 0; a < 3; a++) {
            memcpy(phi[a], snap_phi[a], Nx * sizeof(double));
            memcpy(vel[a], snap_vel[a], Nx * sizeof(double));
        }
    } else {
        printf("  WARNING: no positive peak found in last 5%%, using final state\n");
    }

    /* Compute energy at saved state */
    double E_save = 0, A_save = 0;
    for (int i = 1; i < Nx - 1; i++) {
        for (int a = 0; a < 3; a++) {
            E_save += 0.5 * vel[a][i] * vel[a][i] * dx;
            double dp = (phi[a][i+1] - phi[a][i-1]) / (2.0*dx);
            E_save += 0.5 * dp * dp * dx;
            E_save += 0.5 * m2 * phi[a][i] * phi[a][i] * dx;
            if (fabs(phi[a][i]) > A_save) A_save = fabs(phi[a][i]);
        }
        E_save += lam * (phi[0][i]*phi[1][i] + phi[1][i]*phi[2][i]
                        + phi[2][i]*phi[0][i]) * dx;
        double P = phi[0][i] * phi[1][i] * phi[2][i];
        E_save += 0.5 * mu * P * P / (1.0 + kappa * P * P) * dx;
    }

    /* Measure frequency from DFT (second half) */
    double peak_omega = 0, peak_pow_val = 0;
    int dft_start = n_dft / 2;
    if (n_dft - dft_start > 100) {
        double T = t_hist[n_dft-1] - t_hist[dft_start];
        int nf = 2000;
        for (int k = 1; k < nf; k++) {
            double omega = 2.0 * mass * k / nf;
            double re = 0, im = 0;
            for (int j = dft_start; j < n_dft; j++) {
                double dtj = (j > dft_start) ?
                    (t_hist[j]-t_hist[j-1]) : (t_hist[dft_start+1]-t_hist[dft_start]);
                re += phi0_hist[j] * cos(omega * t_hist[j]) * dtj;
                im += phi0_hist[j] * sin(omega * t_hist[j]) * dtj;
            }
            double pw = (re*re + im*im) / (T*T);
            if (pw > peak_pow_val) { peak_pow_val = pw; peak_omega = omega; }
        }
    }

    /* Fill profile struct */
    for (int a = 0; a < 3; a++) {
        prof.phi[a] = malloc(Nx * sizeof(double));
        prof.vel[a] = malloc(Nx * sizeof(double));
        memcpy(prof.phi[a], phi[a], Nx * sizeof(double));
        memcpy(prof.vel[a], vel[a], Nx * sizeof(double));
    }
    prof.omega   = peak_omega;
    prof.A_peak  = A_save;
    prof.E_total = E_save;

    printf("  omega=%.4f  A_peak=%.4f  E_total=%.4f\n", peak_omega, A_save, E_save);

    /* Cleanup */
    for (int a = 0; a < 3; a++) { free(phi[a]); free(vel[a]); free(acc[a]); }
    for (int a = 0; a < 3; a++) { free(snap_phi[a]); free(snap_vel[a]); }
    free(damp); free(phi0_hist); free(t_hist);

    return prof;
}

/* ===================================================================
 *  Interpolate profile onto interaction grid
 * =================================================================== */
static double interp_profile(const double *data, int N_src, double dx_src,
                              double xmax_src, double x_query)
{
    double xi = (x_query + xmax_src) / dx_src;
    int i0 = (int)floor(xi);
    if (i0 < 0 || i0 >= N_src - 1) return 0.0;
    double frac = xi - i0;
    return data[i0] * (1.0 - frac) + data[i0 + 1] * frac;
}

/* ===================================================================
 *  Energy centroid (split at x=0)
 * =================================================================== */
static void compute_energy_centroids(double *phi[], double *vel[], int Nx,
                                      double dx, double xm, double m2,
                                      double lam,
                                      double *x_left, double *x_right,
                                      double *E_left, double *E_right)
{
    double num_L = 0, den_L = 0, num_R = 0, den_R = 0;
    int ic = Nx / 2;

    for (int i = 1; i < Nx - 1; i++) {
        double x = -xm + i * dx;
        double e = 0;
        for (int a = 0; a < 3; a++) {
            e += 0.5 * vel[a][i] * vel[a][i];
            double dp = (phi[a][i+1] - phi[a][i-1]) / (2.0*dx);
            e += 0.5 * dp * dp;
            e += 0.5 * m2 * phi[a][i] * phi[a][i];
        }
        e += lam * (phi[0][i]*phi[1][i] + phi[1][i]*phi[2][i]
                    + phi[2][i]*phi[0][i]);
        double P = phi[0][i] * phi[1][i] * phi[2][i];
        e += 0.5 * mu * P * P / (1.0 + kappa * P * P);

        double w = fabs(e) * dx;
        if (i < ic) {
            num_L += x * w;
            den_L += w;
        } else {
            num_R += x * w;
            den_R += w;
        }
    }

    *x_left  = (den_L > 1e-30) ? num_L / den_L : 0;
    *x_right = (den_R > 1e-30) ? num_R / den_R : 0;
    *E_left  = den_L;
    *E_right = den_R;
}

/* ===================================================================
 *  Step 2: Two-oscillon interaction
 * =================================================================== */
typedef struct {
    double F_fit;      /* quadratic coefficient c in sep(t) = a + bt + ct^2 */
    double D_init;     /* initial separation from centroid */
    double D_final;    /* final separation */
    int    valid;
} interact_result_t;

static interact_result_t run_interact(double lam, double D, const profile_t *prof)
{
    interact_result_t res = {0};

    int    Nx  = Nx_int;
    double xm  = xmax_int;
    double dx  = 2.0 * xm / (Nx - 1);
    double dx2 = dx * dx;
    double m2  = mass * mass;

    double m2_sym = m2 + 2.0 * lam;
    double kmax = M_PI / dx;
    double dt   = 0.8 * 2.0 / sqrt(kmax * kmax + (m2_sym > m2 ? m2_sym : m2));
    int    Nt   = (int)(t_run / dt) + 1;

    printf("    Interact: lam=%.4f D=%.0f Nx=%d xmax=%.0f dx=%.5f dt=%.6f Nt=%d\n",
           lam, D, Nx, xm, dx, dt, Nt);

    double *phi[3], *vel[3], *acc[3];
    for (int a = 0; a < 3; a++) {
        phi[a] = calloc(Nx, sizeof(double));
        vel[a] = calloc(Nx, sizeof(double));
        acc[a] = calloc(Nx, sizeof(double));
    }

    /* Absorbing boundary */
    double *damp = malloc(Nx * sizeof(double));
    double x_abs = xm * 0.75;
    for (int i = 0; i < Nx; i++) {
        double x = -xm + i * dx;
        double ax = fabs(x);
        if (ax > x_abs) {
            double f = (ax - x_abs) / (xm - x_abs);
            damp[i] = 1.0 - 0.98 * f * f;
        } else {
            damp[i] = 1.0;
        }
    }

    /* Initialize from saved profile: two copies at +/- D/2 */
    for (int a = 0; a < 3; a++) {
        for (int i = 0; i < Nx; i++) {
            double x = -xm + i * dx;
            double f1 = interp_profile(prof->phi[a], prof->N, prof->dx, prof->xmax, x - D/2);
            double f2 = interp_profile(prof->phi[a], prof->N, prof->dx, prof->xmax, x + D/2);
            phi[a][i] = f1 + f2;

            double v1 = interp_profile(prof->vel[a], prof->N, prof->dx, prof->xmax, x - D/2);
            double v2 = interp_profile(prof->vel[a], prof->N, prof->dx, prof->xmax, x + D/2);
            vel[a][i] = v1 + v2;
        }
    }

    /* Compute acceleration with pairwise coupling */
    #define COMPUTE_ACC_INT() do { \
        for (int aa = 0; aa < 3; aa++) { \
            acc[aa][0] = acc[aa][1] = acc[aa][Nx-2] = acc[aa][Nx-1] = 0; \
            int bb = (aa+1)%3, cc = (aa+2)%3; \
            for (int ii = 1; ii < Nx - 1; ii++) { \
                double lapl = (phi[aa][ii+1] - 2.0*phi[aa][ii] + phi[aa][ii-1]) / dx2; \
                double fp = force_triple(phi[0][ii], phi[1][ii], phi[2][ii], aa); \
                acc[aa][ii] = lapl - m2*phi[aa][ii] - lam*(phi[bb][ii]+phi[cc][ii]) + fp; \
            } \
        } \
    } while(0)

    COMPUTE_ACC_INT();

    /* Track separation vs time */
    int max_track = 50000;
    double *sep_hist = malloc(max_track * sizeof(double));
    double *t_track  = malloc(max_track * sizeof(double));
    int n_track = 0;
    int track_every = Nt / max_track;
    if (track_every < 1) track_every = 1;

    /* Output TSV */
    char path[600];
    snprintf(path, sizeof(path), "%s/force_lam%.3f_D%.0f_ts.tsv", outdir, lam, D);
    FILE *fts = fopen(path, "w");
    if (fts) fprintf(fts, "time\tseparation\tx_right\tx_left\tE_total\n");

    int print_every = Nt / 10;
    if (print_every < 1) print_every = 1;

    for (int n = 0; n <= Nt; n++) {
        double t = n * dt;

        if (n % track_every == 0) {
            double xL, xR, EL, ER;
            compute_energy_centroids(phi, vel, Nx, dx, xm, m2, lam, &xL, &xR, &EL, &ER);
            double s = xR - xL;

            if (n_track < max_track) {
                sep_hist[n_track] = s;
                t_track[n_track]  = t;
                n_track++;
            }

            if (fts)
                fprintf(fts, "%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n", t, s, xR, xL, EL+ER);

            if (n % print_every == 0)
                printf("      t=%7.0f  sep=%.4f  xR=%.4f  xL=%.4f  E=%.4f\n",
                       t, s, xR, xL, EL+ER);
        }

        if (n == Nt) break;

        /* Velocity Verlet */
        for (int a = 0; a < 3; a++)
            for (int i = 1; i < Nx - 1; i++)
                vel[a][i] += 0.5 * dt * acc[a][i];
        for (int a = 0; a < 3; a++)
            for (int i = 1; i < Nx - 1; i++)
                phi[a][i] += dt * vel[a][i];
        COMPUTE_ACC_INT();
        for (int a = 0; a < 3; a++)
            for (int i = 1; i < Nx - 1; i++)
                vel[a][i] += 0.5 * dt * acc[a][i];
        for (int a = 0; a < 3; a++)
            for (int i = 0; i < Nx; i++) {
                vel[a][i] *= damp[i];
                phi[a][i] *= damp[i];
            }
    }
    #undef COMPUTE_ACC_INT

    if (fts) fclose(fts);

    /* --- Fit sep(t) using cycle-averaged separation --- */
    if (n_track > 100) {
        double T_avg = 10.0;  /* averaging window > breathing period */

        /* Compute cycle-averaged separation */
        int n_avg = 0;
        double *t_avg_arr = malloc(n_track * sizeof(double));
        double *s_avg_arr = malloc(n_track * sizeof(double));

        for (int j = 0; j < n_track; j++) {
            double tc = t_track[j];
            if (tc > t_run - T_avg) break;

            double sum_s = 0;
            int cnt = 0;
            for (int k = j; k < n_track; k++) {
                if (t_track[k] > tc + T_avg/2) break;
                if (t_track[k] >= tc - T_avg/2) {
                    sum_s += sep_hist[k];
                    cnt++;
                }
            }
            for (int k = j - 1; k >= 0; k--) {
                if (t_track[k] < tc - T_avg/2) break;
                sum_s += sep_hist[k];
                cnt++;
            }

            if (cnt > 0 && tc >= T_avg/2) {
                t_avg_arr[n_avg] = tc;
                s_avg_arr[n_avg] = sum_s / cnt;
                n_avg++;
            }
        }

        /* Fit averaged sep over t in [20, 500] to s = a + b*t + c*t^2 */
        double S[5] = {0};
        double Sy[3] = {0};
        int nfit = 0;

        for (int j = 0; j < n_avg; j++) {
            if (t_avg_arr[j] < 20.0) continue;
            if (t_avg_arr[j] > 500.0) break;
            double tj = t_avg_arr[j];
            double sj = s_avg_arr[j];
            double tk = 1.0;
            for (int k = 0; k < 5; k++) { S[k] += tk; tk *= tj; }
            Sy[0] += sj;
            Sy[1] += sj * tj;
            Sy[2] += sj * tj * tj;
            nfit++;
        }

        if (nfit >= 10) {
            /* 3x3 Gaussian elimination */
            double M[3][4] = {
                {S[0], S[1], S[2], Sy[0]},
                {S[1], S[2], S[3], Sy[1]},
                {S[2], S[3], S[4], Sy[2]}
            };

            for (int col = 0; col < 3; col++) {
                int piv = col;
                for (int r = col+1; r < 3; r++)
                    if (fabs(M[r][col]) > fabs(M[piv][col])) piv = r;
                if (piv != col)
                    for (int c = 0; c < 4; c++) {
                        double tmp = M[col][c]; M[col][c] = M[piv][c]; M[piv][c] = tmp;
                    }
                if (fabs(M[col][col]) < 1e-30) continue;
                for (int r = col+1; r < 3; r++) {
                    double fac = M[r][col] / M[col][col];
                    for (int c = col; c < 4; c++)
                        M[r][c] -= fac * M[col][c];
                }
            }
            double coeff[3] = {0};
            for (int r = 2; r >= 0; r--) {
                if (fabs(M[r][r]) < 1e-30) continue;
                coeff[r] = M[r][3];
                for (int c = r+1; c < 3; c++)
                    coeff[r] -= M[r][c] * coeff[c];
                coeff[r] /= M[r][r];
            }

            /* sep(t) = coeff[0] + coeff[1]*t + coeff[2]*t^2
             * d^2(sep)/dt^2 = 2*coeff[2]
             * Acceleration = 2*coeff[2] */
            res.F_fit = 2.0 * coeff[2];   /* acceleration of separation */
            res.D_init = coeff[0];
            res.valid = 1;

            printf("    Fit(avg): a=%.4f b=%.4e c=%.4e  => accel=%.4e\n",
                   coeff[0], coeff[1], coeff[2], 2.0 * coeff[2]);
        }

        free(t_avg_arr);
        free(s_avg_arr);
    }

    res.D_final = (n_track > 0) ? sep_hist[n_track - 1] : D;

    for (int a = 0; a < 3; a++) { free(phi[a]); free(vel[a]); free(acc[a]); }
    free(damp); free(sep_hist); free(t_track);

    return res;
}

/* ===================================================================
 *  Yukawa fit: ln|F| = ln|F0| - D/lambda  =>  linear fit in D
 * =================================================================== */
static void fit_yukawa(double *D_vals, double *F_vals, int nf,
                       double *F0_out, double *lambda_out)
{
    double sx = 0, sy = 0, sxx = 0, sxy = 0;
    int nfit = 0;
    for (int j = 0; j < nf; j++) {
        if (fabs(F_vals[j]) < 1e-30) continue;
        double lnF = log(fabs(F_vals[j]));
        sx  += D_vals[j];
        sy  += lnF;
        sxx += D_vals[j] * D_vals[j];
        sxy += D_vals[j] * lnF;
        nfit++;
    }
    if (nfit < 2) { *F0_out = 0; *lambda_out = 0; return; }
    double det = nfit * sxx - sx * sx;
    if (fabs(det) < 1e-30) { *F0_out = 0; *lambda_out = 0; return; }
    double slope     = (nfit * sxy - sx * sy) / det;
    double intercept = (sy * sxx - sx * sxy) / det;
    *lambda_out = -1.0 / slope;
    *F0_out     = exp(intercept);
}

/* =================================================================== */
int main(int argc, char **argv)
{
    /* Parse optional overrides */
    for (int i = 1; i < argc - 1; i += 2) {
        if      (!strcmp(argv[i], "-mu"))      mu      = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-kappa"))   kappa   = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-mass"))    mass    = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-A"))       A_init  = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-sigma"))   sig     = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-o"))       strncpy(outdir, argv[i+1], sizeof(outdir)-1);
        else if (!strcmp(argv[i], "-Nx_eq"))   Nx_eq   = atoi(argv[i+1]);
        else if (!strcmp(argv[i], "-Nx_int"))  Nx_int  = atoi(argv[i+1]);
        else if (!strcmp(argv[i], "-xmax_eq"))  xmax_eq  = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-xmax_int")) xmax_int = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-t_equil"))  t_equil  = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-t_run"))    t_run    = atof(argv[i+1]);
    }

    printf("=== V24-S1: Proca Force Measurement ===\n\n");
    printf("mu=%.1f kappa=%.1f mass=%.3f A=%.3f sigma=%.3f\n", mu, kappa, mass, A_init, sig);
    printf("Equil: Nx=%d xmax=%.0f t=%.0f\n", Nx_eq, xmax_eq, t_equil);
    printf("Interact: Nx=%d xmax=%.0f t=%.0f\n\n", Nx_int, xmax_int, t_run);

    /* ========== STEP 1: EQUILIBRATION ========== */
    printf("========== STEP 1: EQUILIBRATION ==========\n\n");

    profile_t profiles[10];
    double omegas[10], E_totals[10], A_peaks[10];
    int alive[10];
    memset(alive, 0, sizeof(alive));

    for (int li = 0; li < n_lam; li++) {
        double lam = lam_list[li];
        printf("\n--- Equilibrating lambda=%.4f ---\n", lam);
        profiles[li] = run_equilibrate(lam);

        omegas[li]   = profiles[li].omega;
        E_totals[li] = profiles[li].E_total;
        A_peaks[li]  = profiles[li].A_peak;

        /* Check if alive: significant amplitude and subgap frequency */
        alive[li] = (A_peaks[li] > 0.05 && omegas[li] > 0.01);
        printf("  => omega=%.4f  E=%.4f  A=%.4f  alive=%s\n",
               omegas[li], E_totals[li], A_peaks[li], alive[li] ? "YES" : "NO");

        /* Save profile */
        char path[600];
        snprintf(path, sizeof(path), "%s/profile_lam%.3f.dat", outdir, lam);
        profile_save(path, &profiles[li]);
        printf("  Saved: %s\n", path);
    }

    /* Equilibration summary */
    printf("\n\n========== EQUILIBRATION SUMMARY ==========\n\n");
    printf("%-8s  %-8s  %-8s  %-10s  %-10s  %-10s  %-10s  %-5s\n",
           "lambda", "m_A", "1/m_A", "omega", "E_total", "A_peak", "xi_pred", "alive");
    printf("--------  --------  --------  ----------  ----------  ----------  ----------  -----\n");
    for (int li = 0; li < n_lam; li++) {
        double lam = lam_list[li];
        double m2a = mass*mass - lam;
        double m_A = (m2a > 0) ? sqrt(m2a) : 0.0;
        double range_A = (m_A > 1e-10) ? 1.0 / m_A : 9999.0;
        double k2 = mass*mass - omegas[li]*omegas[li];
        double xi_pred = (k2 > 0) ? 1.0/sqrt(k2) : 9999.0;
        printf("%-8.4f  %-8.4f  %-8.4f  %-10.4f  %-10.4f  %-10.4f  %-10.4f  %-5s\n",
               lam, m_A, range_A, omegas[li], E_totals[li], A_peaks[li],
               xi_pred, alive[li] ? "YES" : "NO");
    }

    /* ========== STEP 2: TWO-OSCILLON INTERACTIONS ========== */
    printf("\n\n========== STEP 2: TWO-OSCILLON INTERACTIONS ==========\n\n");

    double F_table[10][10];
    double D_init_table[10][10];
    int    valid_table[10][10];
    memset(valid_table, 0, sizeof(valid_table));

    for (int li = 0; li < n_lam; li++) {
        if (!alive[li]) {
            printf("\n--- Skipping lambda=%.4f (not alive) ---\n", lam_list[li]);
            continue;
        }

        printf("\n--- lambda=%.4f: interactions ---\n", lam_list[li]);
        for (int di = 0; di < n_D; di++) {
            double D = D_list[di];
            printf("\n  D=%.0f:\n", D);
            interact_result_t ir = run_interact(lam_list[li], D, &profiles[li]);
            F_table[li][di] = ir.F_fit;
            D_init_table[li][di] = ir.D_init;
            valid_table[li][di] = ir.valid;
            printf("    F=%.4e  D_init=%.4f  D_final=%.4f\n", ir.F_fit, ir.D_init, ir.D_final);
        }
    }

    /* ========== STEP 3: YUKAWA FITS ========== */
    printf("\n\n========== STEP 3: YUKAWA FITS ==========\n\n");

    /* Force summary table */
    char rpath[600];
    snprintf(rpath, sizeof(rpath), "%s/force_summary.tsv", outdir);
    FILE *fsum = fopen(rpath, "w");
    if (fsum) {
        fprintf(fsum, "lambda\tD");
        for (int di = 0; di < n_D; di++)
            fprintf(fsum, "\tF_D%.0f", D_list[di]);
        fprintf(fsum, "\n");
    }

    /* Print detailed F(D) for each lambda */
    for (int li = 0; li < n_lam; li++) {
        if (!alive[li]) continue;

        printf("lambda=%.4f  F(D):\n", lam_list[li]);
        if (fsum) fprintf(fsum, "%.4f", lam_list[li]);
        for (int di = 0; di < n_D; di++) {
            printf("  D=%3.0f: F=%+.4e  %s\n", D_list[di], F_table[li][di],
                   valid_table[li][di] ? "" : "(invalid)");
            if (fsum) fprintf(fsum, "\t%.6e", F_table[li][di]);
        }
        if (fsum) fprintf(fsum, "\n");
        printf("\n");
    }
    if (fsum) fclose(fsum);

    /* Yukawa fits */
    printf("\n%-8s  %-8s  %-10s  %-10s  %-10s  %-10s  %-8s\n",
           "lambda", "m_A", "1/m_A", "lam_fit", "F0_fit", "ratio", "sign");
    printf("--------  --------  ----------  ----------  ----------  ----------  --------\n");

    char ypath[600];
    snprintf(ypath, sizeof(ypath), "%s/yukawa_fits.tsv", outdir);
    FILE *fy = fopen(ypath, "w");
    if (fy) fprintf(fy, "lambda\tm_A\trange_pred\tlam_fit\tF0_fit\tE_total\n");

    for (int li = 0; li < n_lam; li++) {
        if (!alive[li]) continue;

        /* Collect valid F(D) data */
        double Dv[10], Fv[10];
        int nf = 0;
        int sign_pos = 0, sign_neg = 0;
        for (int di = 0; di < n_D; di++) {
            if (valid_table[li][di]) {
                Dv[nf] = D_list[di];
                Fv[nf] = F_table[li][di];
                if (Fv[nf] > 0) sign_pos++;
                else sign_neg++;
                nf++;
            }
        }

        double F0, lam_fit;
        fit_yukawa(Dv, Fv, nf, &F0, &lam_fit);

        double m2a = mass*mass - lam_list[li];
        double m_A = (m2a > 0) ? sqrt(m2a) : 0.0;
        double range_pred = (m_A > 1e-10) ? 1.0 / m_A : 9999.0;

        const char *sign_str = (sign_pos > sign_neg) ? "repuls" :
                               (sign_neg > sign_pos) ? "attract" : "mixed";

        double ratio = (range_pred > 0.01 && range_pred < 9000) ? lam_fit / range_pred : 0.0;

        printf("%-8.4f  %-8.4f  %-10.4f  %-10.4f  %-10.4e  %-10.4f  %-8s\n",
               lam_list[li], m_A, range_pred, lam_fit, F0, ratio, sign_str);

        if (fy) fprintf(fy, "%.4f\t%.6f\t%.6f\t%.6f\t%.6e\t%.6f\n",
                        lam_list[li], m_A, range_pred, lam_fit, F0, E_totals[li]);
    }
    if (fy) fclose(fy);

    /* ========== STEP 4: KEY COMPARISON ========== */
    printf("\n\n========== STEP 4: KEY COMPARISON ==========\n\n");

    /* Find F at D=30 for lambda=0 and lambda=0.99 */
    int li_0 = -1, li_99 = -1, li_999 = -1;
    for (int li = 0; li < n_lam; li++) {
        if (fabs(lam_list[li] - 0.0) < 0.001) li_0 = li;
        if (fabs(lam_list[li] - 0.99) < 0.001) li_99 = li;
        if (fabs(lam_list[li] - 0.999) < 0.001) li_999 = li;
    }

    int di_30 = -1, di_60 = -1;
    for (int di = 0; di < n_D; di++) {
        if (fabs(D_list[di] - 30.0) < 0.1) di_30 = di;
        if (fabs(D_list[di] - 60.0) < 0.1) di_60 = di;
    }

    if (li_0 >= 0 && li_99 >= 0 && di_30 >= 0) {
        printf("CRITICAL TEST: F at D=30\n");
        printf("  lambda=0.00: F(D=30) = %+.4e  (range 1/m_A = %.2f)\n",
               F_table[li_0][di_30], 1.0/mass);
        printf("  lambda=0.99: F(D=30) = %+.4e  (range 1/m_A = %.2f)\n",
               F_table[li_99][di_30], 1.0/sqrt(mass*mass - 0.99));

        double ratio = 0;
        if (fabs(F_table[li_0][di_30]) > 1e-30)
            ratio = fabs(F_table[li_99][di_30]) / fabs(F_table[li_0][di_30]);
        printf("  |F(lam=0.99)| / |F(lam=0)| at D=30 = %.2f\n", ratio);
        printf("  Enhancement expected: exp((D/xi_0 - D/xi_99)) ~ exp(30*(1-0.1)) ~ exp(27) ~ 5e11\n");
        printf("  => Force at D=30 LARGER at lam=0.99? %s\n",
               (fabs(F_table[li_99][di_30]) > fabs(F_table[li_0][di_30])) ? "YES" : "NO");
    }

    if (li_999 >= 0 && di_60 >= 0) {
        printf("\n  lambda=0.999 at D=60: F = %+.4e  (range 1/m_A = %.2f)\n",
               F_table[li_999][di_60], 1.0/sqrt(mass*mass - 0.999));
        printf("  Detectable (|F| > 1e-10)? %s\n",
               (fabs(F_table[li_999][di_60]) > 1e-10) ? "YES" : "NO");
    }

    /* ========== SUMMARY TABLE ========== */
    printf("\n\n========== FULL RESULTS TABLE ==========\n\n");
    printf("%-8s", "lambda");
    for (int di = 0; di < n_D; di++)
        printf("  D=%-6.0f", D_list[di]);
    printf("  lam_fit  1/m_A\n");

    for (int li = 0; li < n_lam; li++) {
        if (!alive[li]) continue;
        printf("%-8.4f", lam_list[li]);
        for (int di = 0; di < n_D; di++) {
            if (valid_table[li][di])
                printf("  %+.2e", F_table[li][di]);
            else
                printf("  %8s", "---");
        }
        /* Yukawa fit */
        double Dv[10], Fv[10]; int nf = 0;
        for (int di = 0; di < n_D; di++)
            if (valid_table[li][di]) {
                Dv[nf] = D_list[di]; Fv[nf] = F_table[li][di]; nf++;
            }
        double F0, lf;
        fit_yukawa(Dv, Fv, nf, &F0, &lf);
        double m2a = mass*mass - lam_list[li];
        double range_pred = (m2a > 0) ? 1.0/sqrt(m2a) : 9999.0;
        printf("  %-8.2f  %-8.2f\n", lf, range_pred);
    }

    /* Cleanup */
    for (int li = 0; li < n_lam; li++)
        profile_free(&profiles[li]);

    printf("\n=== Done ===\n");
    return 0;
}
