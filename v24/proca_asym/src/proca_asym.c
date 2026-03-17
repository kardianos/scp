/*
 * proca_asym.c -- V24-S1b: Asymmetric Oscillon Proca Force Measurement
 *
 * KEY INSIGHT: Symmetric oscillons (phi1=phi2=phi3) don't source the Proca
 * (antisymmetric) channel. We need ASYMMETRIC oscillons to activate Proca exchange.
 *
 * Lagrangian:
 *   L = sum_a [ (1/2)(dt phi_a)^2 - (1/2)(dx phi_a)^2 - (m_a^2/2)phi_a^2 ]
 *     - (mu/2) P^2 / (1 + kappa P^2)
 *     - lambda(phi_1 phi_2 + phi_2 phi_3 + phi_3 phi_1)
 *   P = phi_1 phi_2 phi_3
 *
 * Mass spectrum with pairwise coupling:
 *   antisymmetric mode: m^2_A = m^2 - lambda  (lighter, longer range)
 *   symmetric mode:     m^2_S = m^2 + 2*lambda (heavier)
 *   Proca range = 1/m_A = 1/sqrt(m^2 - lambda) = 10 at lambda=0.99
 *
 * Options:
 *   SYM:  phi1=phi2=phi3 = A*g           (control, no Proca coupling)
 *   180:  phi1=phi2=+A*g, phi3=-A*g      (large antisymmetric content)
 *   UUD:  m1=m2=1.0, m3=0.95, all +A*g  (small asymmetry from mass diff)
 *   SEED: phi1=phi2=A*g, phi3=0.9*A*g    (10% amplitude asymmetry)
 *
 * Protocol:
 *   Phase 1: Equilibrate each option at lambda=0.99 for t=10000.
 *   Phase 2: For each at D={20,30,40,60,80}, measure force (quadratic fit t in [0,500]).
 *   Phase 3: Compare F(asym) vs F(sym) — is Proca channel active?
 *
 * Compile: gcc -O3 -Wall -o proca_asym src/proca_asym.c -lm
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
static double lam      = 0.99;      /* pairwise coupling */
static char   outdir[512] = "data";

/* Equilibration grid */
static int    Nx_eq    = 4000;
static double xmax_eq  = 80.0;
static double t_equil  = 10000.0;

/* Interaction grid */
static int    Nx_int   = 8000;
static double xmax_int = 300.0;
static double t_run    = 2000.0;

/* D scan */
static double D_list[] = {20.0, 30.0, 40.0, 60.0, 80.0};
static int    n_D      = 5;

/* Oscillon option IDs */
#define OPT_SYM   0
#define OPT_180   1
#define OPT_UUD   2
#define OPT_SEED  3
#define N_OPT     4

static const char *opt_names[] = {"SYM", "180", "UUD", "SEED"};

/* UUD mass for field 3 */
static double m3_uud = 0.95;

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
    double S_content;   /* symmetric mode content at peak */
    double A2_content;  /* antisymmetric mode content at peak */
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

/* Compute S, A1, A2 content (integrated) */
static void measure_mode_content(const double *phi0, const double *phi1,
                                  const double *phi2, int Nx, double dx,
                                  double *S_out, double *A1_out, double *A2_out)
{
    double S2 = 0, A12 = 0, A22 = 0;
    for (int i = 0; i < Nx; i++) {
        double p0 = phi0[i], p1 = phi1[i], p2 = phi2[i];
        double S  = (p0 + p1 + p2) / sqrt(3.0);
        double A1 = (p0 - p1) / sqrt(2.0);
        double A2 = (p0 + p1 - 2.0*p2) / sqrt(6.0);
        S2  += S  * S  * dx;
        A12 += A1 * A1 * dx;
        A22 += A2 * A2 * dx;
    }
    *S_out  = sqrt(S2);
    *A1_out = sqrt(A12);
    *A2_out = sqrt(A22);
}

/* ===================================================================
 *  Equilibrate a single oscillon
 * =================================================================== */
static profile_t run_equilibrate(int opt)
{
    profile_t prof = {0};
    prof.N    = Nx_eq;
    prof.xmax = xmax_eq;
    prof.dx   = 2.0 * xmax_eq / (Nx_eq - 1);

    double dx  = prof.dx;
    double dx2 = dx * dx;
    int    Nx  = Nx_eq;
    int    ic  = Nx / 2;

    /* Per-field masses */
    double m2[3];
    for (int a = 0; a < 3; a++) m2[a] = mass * mass;
    if (opt == OPT_UUD) m2[2] = m3_uud * m3_uud;

    /* CFL based on heaviest mode */
    double m2_sym  = mass * mass + 2.0 * lam;
    double kmax = M_PI / dx;
    double dt   = 0.8 * 2.0 / sqrt(kmax * kmax + m2_sym);
    int    Nt   = (int)(t_equil / dt) + 1;

    double m2_anti = mass * mass - lam;
    double m_A = (m2_anti > 0) ? sqrt(m2_anti) : 0.0;

    printf("  Equilibrate [%s]: Nx=%d xmax=%.0f dx=%.5f dt=%.6f Nt=%d\n",
           opt_names[opt], Nx, xmax_eq, dx, dt, Nt);
    printf("    m_A=%.4f  1/m_A=%.2f  m_S=%.4f\n",
           m_A, m_A > 1e-10 ? 1.0/m_A : 9999.0, sqrt(m2_sym));

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

    /* Initialize based on option */
    for (int i = 0; i < Nx; i++) {
        double x = -xmax_eq + i * dx;
        double g = exp(-x * x / (2.0 * sig * sig));
        switch (opt) {
            case OPT_SYM:
                phi[0][i] = phi[1][i] = phi[2][i] = A_init * g;
                break;
            case OPT_180:
                phi[0][i] = phi[1][i] = A_init * g;
                phi[2][i] = -A_init * g;
                break;
            case OPT_UUD:
                phi[0][i] = phi[1][i] = phi[2][i] = A_init * g;
                break;
            case OPT_SEED:
                phi[0][i] = phi[1][i] = A_init * g;
                phi[2][i] = 0.9 * A_init * g;
                break;
        }
    }

    /* Measure initial mode content */
    {
        double S0, A10, A20;
        measure_mode_content(phi[0], phi[1], phi[2], Nx, dx, &S0, &A10, &A20);
        printf("    Initial mode content: S=%.4f  A1=%.4f  A2=%.4f\n", S0, A10, A20);
    }

    /* Macro to compute acceleration with pairwise coupling and per-field masses */
    #define COMPUTE_ACC_EQ() do { \
        for (int aa = 0; aa < 3; aa++) { \
            acc[aa][0] = acc[aa][1] = acc[aa][Nx-2] = acc[aa][Nx-1] = 0; \
            int bb = (aa+1)%3, cc = (aa+2)%3; \
            for (int ii = 1; ii < Nx - 1; ii++) { \
                double lapl = (phi[aa][ii+1] - 2.0*phi[aa][ii] + phi[aa][ii-1]) / dx2; \
                double fp = force_triple(phi[0][ii], phi[1][ii], phi[2][ii], aa); \
                acc[aa][ii] = lapl - m2[aa]*phi[aa][ii] - lam*(phi[bb][ii]+phi[cc][ii]) + fp; \
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
            double cur = fabs(phi[0][ic]);
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
            double Eall = 0, peak = 0;
            for (int i = 1; i < Nx - 1; i++) {
                double e = 0;
                for (int a = 0; a < 3; a++) {
                    e += 0.5 * vel[a][i] * vel[a][i];
                    double dp = (phi[a][i+1] - phi[a][i-1]) / (2.0*dx);
                    e += 0.5 * dp * dp + 0.5 * m2[a] * phi[a][i] * phi[a][i];
                    if (fabs(phi[a][i]) > peak) peak = fabs(phi[a][i]);
                }
                e += lam * (phi[0][i]*phi[1][i] + phi[1][i]*phi[2][i]
                           + phi[2][i]*phi[0][i]);
                double P = phi[0][i] * phi[1][i] * phi[2][i];
                e += 0.5 * mu * P * P / (1.0 + kappa * P * P);
                Eall += e * dx;
            }
            printf("    t=%7.0f  phi0=(%+.4f,%+.4f,%+.4f)  pk=%.4f  E=%+.4f\n",
                   t, phi[0][ic], phi[1][ic], phi[2][ic], peak, Eall);
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
        printf("    Using saved peak snapshot: |phi1(0)|=%.6f\n", snap_phi0_val);
        for (int a = 0; a < 3; a++) {
            memcpy(phi[a], snap_phi[a], Nx * sizeof(double));
            memcpy(vel[a], snap_vel[a], Nx * sizeof(double));
        }
    } else {
        printf("    WARNING: no peak found in last 5%%, using final state\n");
    }

    /* Measure mode content at saved state */
    {
        double S0, A10, A20;
        measure_mode_content(phi[0], phi[1], phi[2], Nx, dx, &S0, &A10, &A20);
        prof.S_content  = S0;
        prof.A2_content = A20;
        printf("    Equilibrated mode content: S=%.4f  A1=%.4f  A2=%.4f\n", S0, A10, A20);
        printf("    A2/S ratio = %.4f\n", S0 > 1e-30 ? A20/S0 : 0.0);
    }

    /* Compute energy at saved state */
    double E_save = 0, A_save = 0;
    for (int i = 1; i < Nx - 1; i++) {
        for (int a = 0; a < 3; a++) {
            E_save += 0.5 * vel[a][i] * vel[a][i] * dx;
            double dp = (phi[a][i+1] - phi[a][i-1]) / (2.0*dx);
            E_save += 0.5 * dp * dp * dx;
            E_save += 0.5 * m2[a] * phi[a][i] * phi[a][i] * dx;
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

    printf("    omega=%.4f  A_peak=%.4f  E_total=%.4f\n", peak_omega, A_save, E_save);

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
                                      double dx, double xm, double m2_arr[3],
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
            e += 0.5 * m2_arr[a] * phi[a][i] * phi[a][i];
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
 *  Two-oscillon interaction
 * =================================================================== */
typedef struct {
    double F_fit;      /* acceleration = 2*c from sep(t) = a + bt + ct^2 */
    double D_init;     /* initial separation from fit */
    double D_final;    /* final separation */
    int    valid;
} interact_result_t;

static interact_result_t run_interact(int opt, double D, const profile_t *prof)
{
    interact_result_t res = {0};

    int    Nx  = Nx_int;
    double xm  = xmax_int;
    double dx  = 2.0 * xm / (Nx - 1);
    double dx2 = dx * dx;

    /* Per-field masses */
    double m2[3];
    for (int a = 0; a < 3; a++) m2[a] = mass * mass;
    if (opt == OPT_UUD) m2[2] = m3_uud * m3_uud;

    double m2_sym = mass * mass + 2.0 * lam;
    double kmax = M_PI / dx;
    double dt   = 0.8 * 2.0 / sqrt(kmax * kmax + m2_sym);
    int    Nt   = (int)(t_run / dt) + 1;

    printf("    Interact [%s] D=%.0f: Nx=%d dx=%.5f dt=%.6f Nt=%d\n",
           opt_names[opt], D, Nx, dx, dt, Nt);

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

    /* Initialize: two copies at +/- D/2 */
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
                acc[aa][ii] = lapl - m2[aa]*phi[aa][ii] - lam*(phi[bb][ii]+phi[cc][ii]) + fp; \
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

    int print_every = Nt / 5;
    if (print_every < 1) print_every = 1;

    for (int n = 0; n <= Nt; n++) {
        double t = n * dt;

        if (n % track_every == 0) {
            double xL, xR, EL, ER;
            compute_energy_centroids(phi, vel, Nx, dx, xm, m2, &xL, &xR, &EL, &ER);
            double s = xR - xL;

            if (n_track < max_track) {
                sep_hist[n_track] = s;
                t_track[n_track]  = t;
                n_track++;
            }

            if (n % print_every == 0)
                printf("      t=%7.0f  sep=%.4f  E=%.4f\n", t, s, EL+ER);
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

    /* Fit sep(t) using cycle-averaged separation */
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

        /* Quadratic fit over t in [20, 500] */
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

            res.F_fit = 2.0 * coeff[2];
            res.D_init = coeff[0];
            res.valid = 1;

            printf("    Fit: D0=%.4f  v=%.4e  accel=%.4e\n",
                   coeff[0], coeff[1], 2.0 * coeff[2]);
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
 *  Yukawa fit: ln|F| = ln|F0| - D/lambda_range
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
        if      (!strcmp(argv[i], "-mu"))       mu       = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-kappa"))    kappa    = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-mass"))     mass     = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-A"))        A_init   = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-sigma"))    sig      = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-lam"))      lam      = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-m3"))       m3_uud   = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-o"))        strncpy(outdir, argv[i+1], sizeof(outdir)-1);
        else if (!strcmp(argv[i], "-Nx_eq"))    Nx_eq    = atoi(argv[i+1]);
        else if (!strcmp(argv[i], "-Nx_int"))   Nx_int   = atoi(argv[i+1]);
        else if (!strcmp(argv[i], "-xmax_eq"))  xmax_eq  = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-xmax_int")) xmax_int = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-t_equil"))  t_equil  = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-t_run"))    t_run    = atof(argv[i+1]);
    }

    printf("=== V24-S1b: Asymmetric Oscillon Proca Force Measurement ===\n\n");
    printf("mu=%.1f kappa=%.1f mass=%.3f A=%.3f sigma=%.3f lambda=%.4f\n",
           mu, kappa, mass, A_init, sig, lam);
    printf("Equil: Nx=%d xmax=%.0f t=%.0f\n", Nx_eq, xmax_eq, t_equil);
    printf("Interact: Nx=%d xmax=%.0f t=%.0f\n", Nx_int, xmax_int, t_run);

    double m2_anti = mass*mass - lam;
    double m_A = (m2_anti > 0) ? sqrt(m2_anti) : 0.0;
    double range_A = (m_A > 1e-10) ? 1.0/m_A : 9999.0;
    printf("Proca mass m_A=%.4f  range 1/m_A=%.2f\n\n", m_A, range_A);

    /* ========== PHASE 1: EQUILIBRATION ========== */
    printf("========== PHASE 1: EQUILIBRATION ==========\n\n");

    profile_t profiles[N_OPT];
    int alive[N_OPT];

    for (int opt = 0; opt < N_OPT; opt++) {
        printf("\n--- Equilibrating option %s ---\n", opt_names[opt]);
        profiles[opt] = run_equilibrate(opt);

        /* Check if alive */
        alive[opt] = (profiles[opt].A_peak > 0.05);
        printf("  => omega=%.4f  E=%.4f  A_peak=%.4f  S=%.4f  A2=%.4f  alive=%s\n\n",
               profiles[opt].omega, profiles[opt].E_total, profiles[opt].A_peak,
               profiles[opt].S_content, profiles[opt].A2_content,
               alive[opt] ? "YES" : "NO");
    }

    /* Equilibration summary */
    printf("\n========== EQUILIBRATION SUMMARY ==========\n\n");
    printf("%-6s  %-8s  %-8s  %-8s  %-8s  %-8s  %-5s\n",
           "Opt", "omega", "E_total", "A_peak", "S", "A2", "alive");
    printf("------  --------  --------  --------  --------  --------  -----\n");
    for (int opt = 0; opt < N_OPT; opt++) {
        printf("%-6s  %-8.4f  %-8.4f  %-8.4f  %-8.4f  %-8.4f  %-5s\n",
               opt_names[opt], profiles[opt].omega, profiles[opt].E_total,
               profiles[opt].A_peak, profiles[opt].S_content, profiles[opt].A2_content,
               alive[opt] ? "YES" : "NO");
    }

    /* ========== PHASE 2: TWO-OSCILLON INTERACTIONS ========== */
    printf("\n\n========== PHASE 2: TWO-OSCILLON INTERACTIONS ==========\n\n");

    double F_table[N_OPT][10];
    double D_init_table[N_OPT][10];
    int    valid_table[N_OPT][10];
    memset(valid_table, 0, sizeof(valid_table));

    for (int opt = 0; opt < N_OPT; opt++) {
        if (!alive[opt]) {
            printf("\n--- Skipping %s (not alive) ---\n", opt_names[opt]);
            continue;
        }

        printf("\n--- %s: interactions ---\n", opt_names[opt]);
        for (int di = 0; di < n_D; di++) {
            double D = D_list[di];
            printf("\n  D=%.0f:\n", D);
            interact_result_t ir = run_interact(opt, D, &profiles[opt]);
            F_table[opt][di] = ir.F_fit;
            D_init_table[opt][di] = ir.D_init;
            valid_table[opt][di] = ir.valid;
            printf("    F=%.4e  D_init=%.4f  D_final=%.4f\n",
                   ir.F_fit, ir.D_init, ir.D_final);
        }
    }

    /* ========== PHASE 3: FORCE COMPARISON ========== */
    printf("\n\n========== PHASE 3: FORCE COMPARISON ==========\n\n");

    /* Force table */
    printf("Force F(D) for each oscillon type:\n\n");
    printf("%-6s", "D");
    for (int opt = 0; opt < N_OPT; opt++)
        if (alive[opt]) printf("  %-12s", opt_names[opt]);
    printf("\n------");
    for (int opt = 0; opt < N_OPT; opt++)
        if (alive[opt]) printf("  ------------");
    printf("\n");

    for (int di = 0; di < n_D; di++) {
        printf("%-6.0f", D_list[di]);
        for (int opt = 0; opt < N_OPT; opt++) {
            if (!alive[opt]) continue;
            if (valid_table[opt][di])
                printf("  %+12.4e", F_table[opt][di]);
            else
                printf("  %12s", "---");
        }
        printf("\n");
    }

    /* Yukawa fits per option */
    printf("\n\nYukawa fits: ln|F| = ln|F0| - D/xi\n\n");
    printf("%-6s  %-10s  %-10s  %-10s  %-8s\n",
           "Opt", "F0", "xi_fit", "1/m_A", "sign");
    printf("------  ----------  ----------  ----------  --------\n");

    for (int opt = 0; opt < N_OPT; opt++) {
        if (!alive[opt]) continue;

        double Dv[10], Fv[10]; int nf = 0;
        int sign_pos = 0, sign_neg = 0;
        for (int di = 0; di < n_D; di++) {
            if (valid_table[opt][di]) {
                Dv[nf] = D_list[di];
                Fv[nf] = F_table[opt][di];
                if (Fv[nf] > 0) sign_pos++; else sign_neg++;
                nf++;
            }
        }

        double F0, xi_fit;
        fit_yukawa(Dv, Fv, nf, &F0, &xi_fit);

        const char *sign_str = (sign_pos > sign_neg) ? "repuls" :
                               (sign_neg > sign_pos) ? "attract" : "mixed";

        printf("%-6s  %-10.4e  %-10.4f  %-10.4f  %-8s\n",
               opt_names[opt], F0, xi_fit, range_A, sign_str);
    }

    /* ========== KEY TEST: D=60 comparison ========== */
    printf("\n\n========== KEY TEST: PROCA CHANNEL DETECTION ==========\n\n");

    int di_60 = -1;
    for (int di = 0; di < n_D; di++)
        if (fabs(D_list[di] - 60.0) < 0.1) di_60 = di;

    if (di_60 >= 0 && alive[OPT_SYM]) {
        double F_sym = fabs(F_table[OPT_SYM][di_60]);
        printf("At D=60 (Proca range 1/m_A = %.1f):\n\n", range_A);
        printf("  %-6s  F(D=60) = %+.4e  (CONTROL: no Proca coupling)\n",
               opt_names[OPT_SYM],
               valid_table[OPT_SYM][di_60] ? F_table[OPT_SYM][di_60] : 0.0);

        for (int opt = 1; opt < N_OPT; opt++) {
            if (!alive[opt] || !valid_table[opt][di_60]) continue;
            double F_asym = fabs(F_table[opt][di_60]);
            double ratio = (F_sym > 1e-30) ? F_asym / F_sym : 0.0;
            printf("  %-6s  F(D=60) = %+.4e  |F_asym/F_sym| = %.2f  %s\n",
                   opt_names[opt], F_table[opt][di_60], ratio,
                   ratio > 3.0 ? "** PROCA ACTIVE **" :
                   ratio > 1.5 ? "* marginal *" : "(no enhancement)");
        }
    }

    /* Also check all D values for the most interesting comparison */
    printf("\n\nRatio |F(asym)/F(sym)| at each D:\n\n");
    printf("%-6s", "D");
    for (int opt = 1; opt < N_OPT; opt++)
        if (alive[opt]) printf("  %-12s", opt_names[opt]);
    printf("\n------");
    for (int opt = 1; opt < N_OPT; opt++)
        if (alive[opt]) printf("  ------------");
    printf("\n");

    for (int di = 0; di < n_D; di++) {
        printf("%-6.0f", D_list[di]);
        double F_sym = fabs(F_table[OPT_SYM][di]);
        for (int opt = 1; opt < N_OPT; opt++) {
            if (!alive[opt]) continue;
            if (valid_table[opt][di] && valid_table[OPT_SYM][di] && F_sym > 1e-30)
                printf("  %12.4f", fabs(F_table[opt][di]) / F_sym);
            else
                printf("  %12s", "---");
        }
        printf("\n");
    }

    /* ========== WRITE SUMMARY TSV ========== */
    {
        char path[600];
        snprintf(path, sizeof(path), "%s/force_summary.tsv", outdir);
        FILE *f = fopen(path, "w");
        if (f) {
            fprintf(f, "option\tD\tF_fit\tvalid\n");
            for (int opt = 0; opt < N_OPT; opt++) {
                if (!alive[opt]) continue;
                for (int di = 0; di < n_D; di++) {
                    fprintf(f, "%s\t%.0f\t%.6e\t%d\n",
                            opt_names[opt], D_list[di], F_table[opt][di],
                            valid_table[opt][di]);
                }
            }
            fclose(f);
            printf("\nForce summary: %s\n", path);
        }
    }

    /* ========== FINAL VERDICT ========== */
    printf("\n\n========== FINAL VERDICT ==========\n\n");

    int proca_detected = 0;
    if (di_60 >= 0 && alive[OPT_SYM] && valid_table[OPT_SYM][di_60]) {
        double F_sym = fabs(F_table[OPT_SYM][di_60]);
        for (int opt = 1; opt < N_OPT; opt++) {
            if (!alive[opt] || !valid_table[opt][di_60]) continue;
            double ratio = (F_sym > 1e-30) ?
                fabs(F_table[opt][di_60]) / F_sym : 0.0;
            if (ratio > 3.0) proca_detected = 1;
        }
    }

    if (proca_detected)
        printf("RESULT: Proca channel IS active. Asymmetric oscillons show enhanced\n"
               "long-range force compared to symmetric control.\n");
    else
        printf("RESULT: No clear Proca enhancement detected. Either:\n"
               "  (a) Asymmetric content equilibrates to zero (180 collapses to SYM)\n"
               "  (b) Proca coupling too weak at these separations\n"
               "  (c) Both channels show similar tail-overlap force\n");

    /* Cleanup */
    for (int opt = 0; opt < N_OPT; opt++)
        profile_free(&profiles[opt]);

    printf("\n=== Done ===\n");
    return 0;
}
