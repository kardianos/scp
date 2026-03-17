/*
 * boost1d.c — V24 Tests 1+2: Boosted oscillon stability and harmonic Doppler shift
 *
 * Phase 1: Equilibrate a stationary oscillon (μ=-20, κ=20, m=1.0, A=0.8, σ=3.0)
 *          Save the full rest-frame solution over several oscillation periods.
 * Phase 2: Lorentz-boost to velocities v=0.0..0.7, evolve, measure dE/dt
 *
 * KEY INSIGHT: The Lorentz boost requires the rest-frame field at (x', t') where
 *   x' = γ(x - x₀), t' = -γv(x - x₀)
 * Different lab-frame x points at t=0 correspond to DIFFERENT rest-frame times.
 * We must interpolate in both space AND time from the rest-frame solution.
 *
 * Compile: gcc -O3 -Wall -o boost1d v24/src/boost1d.c -lm
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ---- Parameters ---- */
static double mu     = -20.0;
static double kappa  = 20.0;
static double mass   = 1.0;
static double A_init = 0.8;
static double sigma  = 3.0;

/* Phase 1: equilibration */
static int    Nx_eq     = 4000;
static double xmax_eq   = 100.0;
static double t_equil   = 5000.0;

/* Phase 2: boosted evolution */
static double t_run     = 2000.0;
static double R_window  = 30.0;

static char outdir[512] = "v24/data";

/* ---- Potential force ---- */
static double force_pot(double p1, double p2, double p3, int a)
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

/* ---- Energy density ---- */
static double energy_density_1(double *phi[3], double *vel[3], int i, int Nx, double dx, double m2)
{
    double e = 0.0;
    for (int a = 0; a < 3; a++) {
        e += 0.5 * vel[a][i] * vel[a][i];
        double dp = (i > 0 && i < Nx-1) ?
            (phi[a][i+1] - phi[a][i-1]) / (2.0*dx) : 0.0;
        e += 0.5 * dp * dp;
        e += 0.5 * m2 * phi[a][i] * phi[a][i];
    }
    double P = phi[0][i] * phi[1][i] * phi[2][i];
    double P2 = P * P;
    e += 0.5 * mu * P2 / (1.0 + kappa * P2);
    return e;
}

/* ============================================================
 *  PHASE 1: Equilibrate and save multiple periods of rest-frame solution
 * ============================================================ */

typedef struct {
    /* Saved rest-frame solution over N_snap snapshots covering several periods */
    double **phi_snaps;    /* phi_snaps[snap_idx][grid_idx] */
    double **vel_snaps;    /* vel_snaps[snap_idx][grid_idx] */
    double *t_snaps;       /* time of each snapshot */
    int N_snap;            /* number of snapshots */
    int Nx;                /* grid size */
    double dx;             /* grid spacing */
    double xmax;           /* half-domain */
    double omega;          /* breathing frequency */
    double peak_amp;       /* peak amplitude */
    double E_equil;        /* energy at equilibration end */
} EquilResult;

static double interp_1d(const double *f, int N, double dx, double xmax, double x)
{
    double xi = (x + xmax) / dx;
    int i0 = (int)floor(xi);
    if (i0 < 0) return 0.0;
    if (i0 >= N - 1) return 0.0;
    double frac = xi - i0;
    return f[i0] * (1.0 - frac) + f[i0 + 1] * frac;
}

/* Interpolate rest-frame solution at (x', t') using bilinear interpolation */
static void interp_rest_frame(const EquilResult *eq, double xp, double tp,
                               double *phi_out, double *vel_out)
{
    /* Find bracketing time snapshots.
     * We use periodic interpolation within the saved window. */
    double t0 = eq->t_snaps[0];
    double t1 = eq->t_snaps[eq->N_snap - 1];
    double T_span = t1 - t0;

    /* Wrap tp into the saved time window */
    double tp_mod = tp - t0;
    /* If tp is outside the window, wrap modulo T_span */
    if (tp_mod < 0) tp_mod += T_span * (1 + (int)(-tp_mod / T_span));
    if (tp_mod >= T_span) tp_mod -= T_span * (int)(tp_mod / T_span);
    double tp_eff = t0 + tp_mod;

    /* Find bracketing snapshots */
    int j0 = 0;
    for (int j = 0; j < eq->N_snap - 1; j++) {
        if (eq->t_snaps[j+1] > tp_eff) { j0 = j; break; }
    }
    int j1 = j0 + 1;
    if (j1 >= eq->N_snap) j1 = eq->N_snap - 1;

    double dt_snap = eq->t_snaps[j1] - eq->t_snaps[j0];
    double frac_t = (dt_snap > 1e-30) ? (tp_eff - eq->t_snaps[j0]) / dt_snap : 0.0;

    /* Interpolate phi and vel in space at both time snapshots, then interpolate in time */
    double phi0 = interp_1d(eq->phi_snaps[j0], eq->Nx, eq->dx, eq->xmax, xp);
    double phi1 = interp_1d(eq->phi_snaps[j1], eq->Nx, eq->dx, eq->xmax, xp);
    *phi_out = phi0 * (1.0 - frac_t) + phi1 * frac_t;

    double vel0 = interp_1d(eq->vel_snaps[j0], eq->Nx, eq->dx, eq->xmax, xp);
    double vel1 = interp_1d(eq->vel_snaps[j1], eq->Nx, eq->dx, eq->xmax, xp);
    *vel_out = vel0 * (1.0 - frac_t) + vel1 * frac_t;
}

static EquilResult equilibrate(void)
{
    printf("=== PHASE 1: Equilibrating stationary oscillon ===\n");
    printf("  mu=%.1f kappa=%.1f mass=%.1f A=%.3f sigma=%.3f\n",
           mu, kappa, mass, A_init, sigma);

    int Nx = Nx_eq;
    double xmax = xmax_eq;
    double dx = 2.0 * xmax / (Nx - 1);
    double dx2 = dx * dx;
    double m2 = mass * mass;

    /* CFL */
    double kmax = M_PI / dx;
    double dt = 0.8 * 2.0 / sqrt(kmax * kmax + m2);
    int Nt = (int)(t_equil / dt) + 1;

    printf("  Nx=%d xmax=%.1f dx=%.5f dt=%.6f Nt=%d\n", Nx, xmax, dx, dt, Nt);

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

    /* Initialize: Gaussians (symmetric triad) */
    for (int a = 0; a < 3; a++)
        for (int i = 0; i < Nx; i++) {
            double x = -xmax + i * dx;
            phi[a][i] = A_init * exp(-x * x / (2.0 * sigma * sigma));
        }

    #define COMPUTE_ACC_EQ() do { \
        for (int a = 0; a < 3; a++) { \
            acc[a][0] = acc[a][1] = acc[a][Nx-2] = acc[a][Nx-1] = 0; \
            for (int i = 1; i < Nx - 1; i++) { \
                double lapl = (phi[a][i+1] - 2.0*phi[a][i] + phi[a][i-1]) / dx2; \
                double fp = force_pot(phi[0][i], phi[1][i], phi[2][i], a); \
                acc[a][i] = lapl - m2*phi[a][i] + fp; \
            } \
        } \
    } while(0)

    COMPUTE_ACC_EQ();

    /* DFT storage for frequency measurement */
    int max_dft = 50000;
    double *phi0_hist = malloc(max_dft * sizeof(double));
    double *t_hist_dft = malloc(max_dft * sizeof(double));
    int n_dft = 0;
    int dft_every = Nt / max_dft;
    if (dft_every < 1) dft_every = 1;

    int ic = Nx / 2;
    int print_every = Nt / 20;
    if (print_every < 1) print_every = 1;

    /* Run equilibration to t_equil */
    for (int n = 0; n <= Nt; n++) {
        double t = n * dt;

        if (n % dft_every == 0 && n_dft < max_dft) {
            phi0_hist[n_dft] = phi[0][ic];
            t_hist_dft[n_dft] = t;
            n_dft++;
        }

        if (n % print_every == 0) {
            double Et = 0;
            for (int i = 1; i < Nx - 1; i++)
                Et += energy_density_1(phi, vel, i, Nx, dx, m2) * dx;
            printf("  t=%7.0f  phi0=%.4f  E=%.4f\n", t, phi[0][ic], Et);
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

        for (int a = 0; a < 3; a++)
            for (int i = 0; i < Nx; i++) {
                vel[a][i] *= damp[i];
                phi[a][i] *= damp[i];
            }
    }

    /* DFT to find omega (second half) */
    int dft_start = n_dft / 2;
    double T_dft = t_hist_dft[n_dft-1] - t_hist_dft[dft_start];
    double peak_pow = 0, peak_om = 0;
    int nf = 1000;
    for (int k = 0; k < nf; k++) {
        double omega = 3.0 * mass * k / nf;
        double re = 0, im = 0;
        for (int j = dft_start; j < n_dft; j++) {
            double dtj = (j > dft_start) ?
                (t_hist_dft[j] - t_hist_dft[j-1]) : (t_hist_dft[dft_start+1] - t_hist_dft[dft_start]);
            re += phi0_hist[j] * cos(omega * t_hist_dft[j]) * dtj;
            im += phi0_hist[j] * sin(omega * t_hist_dft[j]) * dtj;
        }
        double pw = (re * re + im * im) / (T_dft * T_dft);
        if (pw > peak_pow) { peak_pow = pw; peak_om = omega; }
    }

    double E_final = 0;
    for (int i = 1; i < Nx - 1; i++)
        E_final += energy_density_1(phi, vel, i, Nx, dx, m2) * dx;

    printf("\n  Equilibration complete.\n");
    printf("  E_final = %.4f, omega = %.4f (mass gap = %.4f)\n", E_final, peak_om, mass);

    /* Now continue evolving and save snapshots covering several periods.
     * Period T = 2*pi/omega. Save ~100 snapshots per period over 5 periods.
     * For v=0.7, gamma=1.4, the time offset at the oscillon edge (|x|~15) is
     * gamma*v*15 = 1.4*0.7*15 = 14.7, which is ~2 periods. Save 5 periods. */
    double T_period = 2.0 * M_PI / peak_om;
    int n_periods = 5;
    double t_save_span = n_periods * T_period;
    int snaps_per_period = 100;
    int N_snap = n_periods * snaps_per_period;
    double dt_snap = t_save_span / N_snap;

    printf("  Saving %d snapshots over %.1f time units (%.1f periods)\n",
           N_snap, t_save_span, (double)n_periods);

    EquilResult res;
    res.phi_snaps = malloc(N_snap * sizeof(double*));
    res.vel_snaps = malloc(N_snap * sizeof(double*));
    res.t_snaps   = malloc(N_snap * sizeof(double));
    res.N_snap = N_snap;
    res.Nx = Nx;
    res.dx = dx;
    res.xmax = xmax;
    res.omega = peak_om;
    res.peak_amp = 0;
    res.E_equil = E_final;

    for (int s = 0; s < N_snap; s++) {
        res.phi_snaps[s] = malloc(Nx * sizeof(double));
        res.vel_snaps[s] = malloc(Nx * sizeof(double));
    }

    /* Continue evolving from current state, saving snapshots */
    double t_base = t_equil;
    int snap_idx = 0;
    double next_snap_t = t_base;

    /* Reset acceleration */
    COMPUTE_ACC_EQ();

    int Nt2 = (int)(t_save_span / dt) + 1;
    for (int n = 0; n <= Nt2; n++) {
        double t = t_base + n * dt;

        /* Save snapshot if due */
        if (snap_idx < N_snap && t >= next_snap_t) {
            memcpy(res.phi_snaps[snap_idx], phi[0], Nx * sizeof(double));
            memcpy(res.vel_snaps[snap_idx], vel[0], Nx * sizeof(double));
            res.t_snaps[snap_idx] = t - t_base;  /* time relative to start of save window */

            double pk = fabs(phi[0][ic]);
            if (pk > res.peak_amp) res.peak_amp = pk;

            snap_idx++;
            next_snap_t = t_base + snap_idx * dt_snap;
        }

        if (n == Nt2) break;

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

        for (int a = 0; a < 3; a++)
            for (int i = 0; i < Nx; i++) {
                vel[a][i] *= damp[i];
                phi[a][i] *= damp[i];
            }
    }

    printf("  Saved %d snapshots, peak amp = %.4f\n\n", snap_idx, res.peak_amp);

    /* Save first snapshot as profile for reference */
    {
        char profpath[600];
        snprintf(profpath, sizeof(profpath), "%s/equil_profile.tsv", outdir);
        FILE *fp = fopen(profpath, "w");
        fprintf(fp, "x\tphi\tvel\n");
        for (int i = 0; i < Nx; i++) {
            double x = -xmax + i * dx;
            fprintf(fp, "%.6f\t%.8e\t%.8e\n", x, res.phi_snaps[0][i], res.vel_snaps[0][i]);
        }
        fclose(fp);
    }

    for (int a = 0; a < 3; a++) { free(phi[a]); free(vel[a]); free(acc[a]); }
    free(damp); free(phi0_hist); free(t_hist_dft);

    return res;
}

/* ============================================================
 *  PHASE 2: Boost and evolve at velocity v
 *
 *  The Lorentz boost:
 *    Lab frame (x, t=0) ↔ Rest frame (x', t') where
 *      x' = γ(x - x₀)
 *      t' = -γv(x - x₀)
 *
 *    φ_lab(x, 0) = φ_rest(x', t')
 *    ∂_t φ_lab(x, 0) = γ[-v·∂_{x'}φ_rest + ∂_{t'}φ_rest] at (x',t')
 * ============================================================ */

typedef struct {
    double v, gamma, E_initial, E_final, dE_dt, omega, A_peak, harm_2w_amp;
} BoostResult;

static BoostResult run_boost(double v, const EquilResult *eq)
{
    double gamma = 1.0 / sqrt(1.0 - v * v);
    printf("=== Boosted evolution: v=%.2f, gamma=%.4f ===\n", v, gamma);

    /* Choose grid size based on velocity:
     * Need enough room for the oscillon to travel v*t_run + buffer.
     * Keep dx matched to equilibration grid for resolution. */
    double dx = eq->dx;  /* SAME resolution as equilibration */
    double dx2 = dx * dx;
    double m2 = mass * mass;

    /* Domain must contain: oscillon start + travel + tail + absorbing boundary.
     * Place oscillon at x=0 at t=0. It will be at x=v*t_run at the end.
     * Need xmax > v*t_run + R_window + absorbing zone.
     * Absorbing zone = 0.2 * xmax. So xmax * 0.8 > v*t_run + R_window + 50. */
    double travel = v * t_run;
    double needed = travel + R_window + 50.0;
    double xmax = needed / 0.7;
    if (xmax < 100.0) xmax = 100.0;  /* minimum domain */

    int Nx = (int)(2.0 * xmax / dx) + 1;

    /* CFL */
    double kmax = M_PI / dx;
    double dt = 0.8 * 2.0 / sqrt(kmax * kmax + m2);
    int Nt = (int)(t_run / dt) + 1;

    printf("  Nx=%d xmax=%.1f dx=%.5f dt=%.6f Nt=%d\n", Nx, xmax, dx, dt, Nt);
    printf("  Travel distance: %.1f, total domain: %.1f\n", travel, 2*xmax);

    double *phi[3], *vel[3], *acc[3];
    for (int a = 0; a < 3; a++) {
        phi[a] = calloc(Nx, sizeof(double));
        vel[a] = calloc(Nx, sizeof(double));
        acc[a] = calloc(Nx, sizeof(double));
    }

    /* Absorbing boundary: outer 20% */
    double *dampv = malloc(Nx * sizeof(double));
    double x_abs = xmax * 0.80;
    for (int i = 0; i < Nx; i++) {
        double x = -xmax + i * dx;
        double ax = fabs(x);
        if (ax > x_abs) {
            double f = (ax - x_abs) / (xmax - x_abs);
            dampv[i] = 1.0 - 0.98 * f * f;
        } else {
            dampv[i] = 1.0;
        }
    }

    /* Initialize boosted profile using spacetime interpolation of rest-frame solution.
     * x' = γ(x - 0), t' = -γv(x - 0)  [x₀ = 0]
     * Place center at x=0 in the lab frame. */
    double x0 = 0.0;

    /* Time offset range: t' ∈ [-γv*xmax_effective, +γv*xmax_effective]
     * For the oscillon, only |x'| < ~20 matters, so |x| < 20/γ.
     * t' range: |t'| < γv * 20/γ = 20v. For v=0.7: |t'| < 14.
     * Our saved snapshots span T_period * n_periods ≈ 7.2*5 = 36. OK. */

    double t_snap_span = eq->t_snaps[eq->N_snap - 1] - eq->t_snaps[0];
    printf("  Rest-frame snapshot span: %.2f\n", t_snap_span);
    printf("  Max time offset: |t'| < %.2f (at |x|=20/γ)\n", gamma * v * 20.0 / gamma);

    for (int i = 0; i < Nx; i++) {
        double x = -xmax + i * dx;
        double xp = gamma * (x - x0);
        double tp = -gamma * v * (x - x0);

        /* Get rest-frame phi and vel at (x', t') */
        double phi_rest, vel_rest;
        interp_rest_frame(eq, xp, tp, &phi_rest, &vel_rest);

        /* Also need spatial derivative at (x', t') for the Lorentz transform of velocity.
         * ∂_t φ_lab = γ(-v·∂_{x'} φ_rest + ∂_{t'} φ_rest)
         * We approximate ∂_{x'} φ_rest by finite difference in x'. */
        double phi_rest_plus, vel_rest_plus, phi_rest_minus, vel_rest_minus;
        double dxp = eq->dx;  /* same as rest-frame grid spacing */
        interp_rest_frame(eq, xp + dxp, tp, &phi_rest_plus, &vel_rest_plus);
        interp_rest_frame(eq, xp - dxp, tp, &phi_rest_minus, &vel_rest_minus);
        double dphi_dxp = (phi_rest_plus - phi_rest_minus) / (2.0 * dxp);

        double boosted_phi = phi_rest;
        double boosted_vel = gamma * (-v * dphi_dxp + vel_rest);

        for (int a = 0; a < 3; a++) {
            phi[a][i] = boosted_phi;
            vel[a][i] = boosted_vel;
        }
    }

    /* Compute acceleration */
    #define COMPUTE_ACC_RUN() do { \
        for (int a = 0; a < 3; a++) { \
            acc[a][0] = acc[a][1] = acc[a][Nx-2] = acc[a][Nx-1] = 0; \
            for (int i = 1; i < Nx - 1; i++) { \
                double lapl = (phi[a][i+1] - 2.0*phi[a][i] + phi[a][i-1]) / dx2; \
                double fp = force_pot(phi[0][i], phi[1][i], phi[2][i], a); \
                acc[a][i] = lapl - m2*phi[a][i] + fp; \
            } \
        } \
    } while(0)

    COMPUTE_ACC_RUN();

    /* Compute initial energy in window around x0 */
    double E_init_window = 0;
    for (int i = 1; i < Nx - 1; i++) {
        double x = -xmax + i * dx;
        if (fabs(x - x0) < R_window)
            E_init_window += energy_density_1(phi, vel, i, Nx, dx, m2) * dx;
    }
    printf("  E_initial (window) = %.6f\n", E_init_window);

    /* Time series output */
    char tspath[600];
    snprintf(tspath, sizeof(tspath), "%s/boost_v%.1f_ts.tsv", outdir, v);
    FILE *fts = fopen(tspath, "w");
    fprintf(fts, "time\tE_window\tpeak_amp\tx_center\n");

    /* DFT storage */
    int max_dft = 50000;
    double *phi_center_hist = malloc(max_dft * sizeof(double));
    double *t_dft_hist = malloc(max_dft * sizeof(double));
    int n_dft = 0;
    int dft_every = Nt / max_dft;
    if (dft_every < 1) dft_every = 1;

    /* Energy time series for dE/dt fit */
    int max_ets = 20000;
    double *E_ts = malloc(max_ets * sizeof(double));
    double *t_ets = malloc(max_ets * sizeof(double));
    int n_ets = 0;
    int ets_every = Nt / max_ets;
    if (ets_every < 1) ets_every = 1;

    int rec_every = Nt / 10000;
    if (rec_every < 1) rec_every = 1;
    int print_every = Nt / 20;
    if (print_every < 1) print_every = 1;

    double last_peak_amp = 0;

    for (int n = 0; n <= Nt; n++) {
        double t = n * dt;
        double x_expected = x0 + v * t;

        int do_rec = (n % rec_every == 0);
        int do_print = (n % print_every == 0);
        int do_dft = (n % dft_every == 0 && n_dft < max_dft);
        int do_ets = (n % ets_every == 0 && n_ets < max_ets);

        if (do_rec || do_print || do_dft || do_ets) {
            double wsum = 0, wxsum = 0, E_win = 0, pk = 0;
            int i_exp = (int)((x_expected + xmax) / dx);
            int i_lo = i_exp - (int)(R_window / dx);
            int i_hi = i_exp + (int)(R_window / dx);
            if (i_lo < 1) i_lo = 1;
            if (i_hi > Nx - 2) i_hi = Nx - 2;

            for (int i = i_lo; i <= i_hi; i++) {
                double x = -xmax + i * dx;
                double ed = energy_density_1(phi, vel, i, Nx, dx, m2);
                E_win += ed * dx;
                wsum += ed * dx;
                wxsum += ed * x * dx;
                double ap = fabs(phi[0][i]);
                if (ap > pk) pk = ap;
            }
            double x_center = (wsum > 1e-20) ? wxsum / wsum : x_expected;
            last_peak_amp = pk;

            if (do_rec)
                fprintf(fts, "%.6f\t%.8e\t%.8e\t%.6f\n", t, E_win, pk, x_center);

            if (do_print)
                printf("  t=%7.0f  x_c=%8.2f (exp=%8.2f)  E_win=%.4f  pk=%.4f\n",
                       t, x_center, x_expected, E_win, pk);

            if (do_dft) {
                double xi = (x_center + xmax) / dx;
                int ii = (int)xi;
                if (ii >= 1 && ii < Nx - 1) {
                    double frac = xi - ii;
                    phi_center_hist[n_dft] = phi[0][ii] * (1.0 - frac) + phi[0][ii+1] * frac;
                } else {
                    phi_center_hist[n_dft] = 0;
                }
                t_dft_hist[n_dft] = t;
                n_dft++;
            }

            if (do_ets) {
                E_ts[n_ets] = E_win;
                t_ets[n_ets] = t;
                n_ets++;
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
        COMPUTE_ACC_RUN();
        for (int a = 0; a < 3; a++)
            for (int i = 1; i < Nx - 1; i++)
                vel[a][i] += 0.5 * dt * acc[a][i];

        for (int a = 0; a < 3; a++)
            for (int i = 0; i < Nx; i++) {
                vel[a][i] *= dampv[i];
                phi[a][i] *= dampv[i];
            }
    }

    fclose(fts);

    /* ---- DFT analysis ---- */
    int dft_start2 = n_dft / 2;
    double T_dft2 = t_dft_hist[n_dft-1] - t_dft_hist[dft_start2];
    double pk_pow = 0, pk_om = 0;
    int nf2 = 2000;
    double om_max = 3.0 * mass;
    double *spectrum = calloc(nf2, sizeof(double));
    double *om_arr = calloc(nf2, sizeof(double));

    for (int k = 0; k < nf2; k++) {
        double omega = om_max * k / nf2;
        om_arr[k] = omega;
        double re = 0, im = 0;
        for (int j = dft_start2; j < n_dft; j++) {
            double dtj = (j > dft_start2) ?
                (t_dft_hist[j] - t_dft_hist[j-1]) : (t_dft_hist[dft_start2+1] - t_dft_hist[dft_start2]);
            re += phi_center_hist[j] * cos(omega * t_dft_hist[j]) * dtj;
            im += phi_center_hist[j] * sin(omega * t_dft_hist[j]) * dtj;
        }
        double pw = (re * re + im * im) / (T_dft2 * T_dft2);
        spectrum[k] = pw;
        if (pw > pk_pow) { pk_pow = pw; pk_om = omega; }
    }

    /* Find 2nd harmonic amplitude */
    double target_2w = 2.0 * pk_om;
    double harm_2w_pow = 0;
    for (int k = 0; k < nf2; k++) {
        if (fabs(om_arr[k] - target_2w) < om_max / nf2) {
            if (spectrum[k] > harm_2w_pow) harm_2w_pow = spectrum[k];
        }
    }
    double harm_2w_amp = (pk_pow > 1e-30) ? sqrt(harm_2w_pow / pk_pow) : 0;

    /* Save spectrum */
    {
        char specpath[600];
        snprintf(specpath, sizeof(specpath), "%s/boost_v%.1f_spectrum.tsv", outdir, v);
        FILE *fspec = fopen(specpath, "w");
        fprintf(fspec, "omega\tpower\n");
        for (int k = 0; k < nf2; k++)
            fprintf(fspec, "%.6f\t%.8e\n", om_arr[k], spectrum[k]);
        fclose(fspec);
    }
    free(spectrum); free(om_arr);

    /* ---- dE/dt linear fit (second half) ---- */
    int ets_start = n_ets / 2;
    double sum_t = 0, sum_E = 0, sum_tE = 0, sum_t2 = 0;
    int nfit = n_ets - ets_start;
    for (int j = ets_start; j < n_ets; j++) {
        sum_t  += t_ets[j];
        sum_E  += E_ts[j];
        sum_tE += t_ets[j] * E_ts[j];
        sum_t2 += t_ets[j] * t_ets[j];
    }
    double dE_dt = 0;
    if (nfit > 2) {
        double det = nfit * sum_t2 - sum_t * sum_t;
        if (fabs(det) > 1e-30)
            dE_dt = (nfit * sum_tE - sum_t * sum_E) / det;
    }
    double E_final_window = (n_ets > 0) ? E_ts[n_ets - 1] : 0;

    printf("  omega = %.4f, 2w = %.4f (gap = %.4f)\n", pk_om, 2*pk_om, mass);
    printf("  E_initial = %.4f, E_final = %.4f, dE/dt = %.6e\n",
           E_init_window, E_final_window, dE_dt);
    printf("  2w harmonic relative amplitude = %.6e\n", harm_2w_amp);

    double doppler_back = (v > 0.01) ? sqrt((1.0 - v) / (1.0 + v)) : 1.0;
    double freq_2w_back = 2.0 * pk_om * doppler_back;
    printf("  Backward Doppler of 2w: %.4f (gap=%.4f, below? %s)\n\n",
           freq_2w_back, mass, freq_2w_back < mass ? "YES" : "NO");

    BoostResult res2;
    res2.v = v;
    res2.gamma = gamma;
    res2.E_initial = E_init_window;
    res2.E_final = E_final_window;
    res2.dE_dt = dE_dt;
    res2.omega = pk_om;
    res2.A_peak = last_peak_amp;
    res2.harm_2w_amp = harm_2w_amp;

    for (int a = 0; a < 3; a++) { free(phi[a]); free(vel[a]); free(acc[a]); }
    free(dampv); free(phi_center_hist); free(t_dft_hist);
    free(E_ts); free(t_ets);

    return res2;
}

/* ============================================================ */
int main(int argc, char **argv)
{
    for (int i = 1; i < argc - 1; i += 2) {
        if      (!strcmp(argv[i], "-mu"))      mu      = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-kappa"))   kappa   = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-mass"))    mass    = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-A"))       A_init  = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-sigma"))   sigma   = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-o"))       strncpy(outdir, argv[i+1], sizeof(outdir)-1);
        else if (!strcmp(argv[i], "-t_equil")) t_equil = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-t_run"))   t_run   = atof(argv[i+1]);
        else { fprintf(stderr, "Unknown flag: %s\n", argv[i]); exit(1); }
    }

    printf("V24 Boost Test: mu=%.1f kappa=%.1f mass=%.1f A=%.3f sigma=%.3f\n\n",
           mu, kappa, mass, A_init, sigma);

    /* Phase 1 */
    EquilResult eq = equilibrate();

    /* Phase 2 */
    double velocities[] = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7};
    int nv = sizeof(velocities) / sizeof(velocities[0]);
    BoostResult *results = malloc(nv * sizeof(BoostResult));

    for (int iv = 0; iv < nv; iv++) {
        results[iv] = run_boost(velocities[iv], &eq);
    }

    /* ---- Summary ---- */
    printf("\n========== SUMMARY ==========\n");
    printf("v\tgamma\tE_init\tE_final\tdE/dt\t\tomega\tA_peak\tharm_2w\t2w_back\tbelow_gap?\n");
    for (int iv = 0; iv < nv; iv++) {
        BoostResult *r = &results[iv];
        double db = (r->v > 0.01) ? sqrt((1.0 - r->v) / (1.0 + r->v)) : 1.0;
        double f2b = 2.0 * r->omega * db;
        printf("%.1f\t%.4f\t%.4f\t%.4f\t%.6e\t%.4f\t%.4f\t%.4e\t%.4f\t%s\n",
               r->v, r->gamma, r->E_initial, r->E_final, r->dE_dt,
               r->omega, r->A_peak, r->harm_2w_amp, f2b,
               f2b < mass ? "YES" : "NO");
    }

    /* Save summary TSV */
    char sumpath[600];
    snprintf(sumpath, sizeof(sumpath), "%s/boost_summary.tsv", outdir);
    FILE *fsum = fopen(sumpath, "w");
    fprintf(fsum, "v\tgamma\tE_initial\tE_final\tdE_dt\tomega\tA_peak\tharm_2w_amp\tfreq_2w_back\tbelow_gap\n");
    for (int iv = 0; iv < nv; iv++) {
        BoostResult *r = &results[iv];
        double db = (r->v > 0.01) ? sqrt((1.0 - r->v) / (1.0 + r->v)) : 1.0;
        double f2b = 2.0 * r->omega * db;
        fprintf(fsum, "%.2f\t%.6f\t%.8e\t%.8e\t%.8e\t%.6f\t%.8e\t%.8e\t%.6f\t%d\n",
                r->v, r->gamma, r->E_initial, r->E_final, r->dE_dt,
                r->omega, r->A_peak, r->harm_2w_amp, f2b,
                f2b < mass ? 1 : 0);
    }
    fclose(fsum);
    printf("\nSummary written to %s\n", sumpath);

    /* Cleanup */
    for (int s = 0; s < eq.N_snap; s++) {
        free(eq.phi_snaps[s]);
        free(eq.vel_snaps[s]);
    }
    free(eq.phi_snaps); free(eq.vel_snaps); free(eq.t_snaps);
    free(results);

    return 0;
}
