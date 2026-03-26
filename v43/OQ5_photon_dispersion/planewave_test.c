/*
 * planewave_test.c — Plane-wave pulse propagation in Cosserat field
 *
 * Tests OQ5: photon dispersion relation in coupled phi-theta system.
 *
 * Physics: Full 6-field Cosserat
 *   d²phi_a/dt² = nabla²phi_a - m²phi_a - dV/dphi_a + eta*curl(theta)_a
 *   d²theta_a/dt² = nabla²theta_a + eta*curl(phi)_a
 *   V(P) = (mu/2)*P²/(1+kappa*P²),  P = phi_0*phi_1*phi_2
 *
 * Initial condition: plane-wave pulse in theta_z only
 *   theta_z = A * exp(-(x-x0)²/(2*sigma_x²)) * sin(k0*x)
 *   phi_a = A_bg * cos(k_bg*z + 2*pi*a/3)   (standing background)
 *
 * Two runs: eta=0 (control, v_g=1.0) and eta=0.5 (coupled, v_g~0.964)
 *
 * Coupled dispersion: (w²-k²)(w²-k²-m²) = eta²*k²
 *   Photon branch: w² = k² + m²/2 - sqrt(m⁴/4 + eta²*k²)
 *   At k0=2.0, eta=0.5: v_g = 0.964
 *
 * Build: gcc -O3 -march=native -fopenmp -o planewave_test planewave_test.c -lm
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>

/* Grid parameters */
#define N    128
#define L    25.0
#define DX   (2.0 * L / N)
#define NNN  (N * N * N)

/* Physics parameters */
#define M2     2.25
#define MU     (-41.345)
#define KAPPA  50.0

/* Time stepping */
#define DT       (0.1 * DX)
#define T_END    50.0
#define MEAS_DT  1.0

/* Initial condition */
#define AMP      0.01
#define K0       2.0
#define SIGMA_X  5.0      /* compact pulse: FWHM ~ 12 code units, fits in domain */
#define X0       (-15.0)
#define A_BG     0.0      /* zero background: coupling comes from induced phi */

/* Absorbing boundary: damp near ALL edges to prevent wraparound */
#define DAMP_START  (0.8 * L)   /* = 20.0 — gives 5 code units of damping */
#define DAMP_MAX    4.0         /* strong damping in boundary zone */

/* Phase offsets for phi background */
static const double delta[3] = {0.0, 3.0005, 4.4325};

/* 3D index */
#define IDX(i,j,k) ((i)*N*N + (j)*N + (k))

/* Periodic wrap */
static inline int pw(int i) {
    return ((i % N) + N) % N;
}

/* Grid coordinate: [-L, L) */
static inline double coord(int i) {
    return -L + (i + 0.5) * DX;
}

/* Field arrays */
static double *phi[3], *phi_dot[3], *phi_acc[3], *phi_acc_old[3];
static double *theta[3], *theta_dot[3], *theta_acc[3], *theta_acc_old[3];

void alloc_fields(void) {
    for (int f = 0; f < 3; f++) {
        phi[f]         = calloc(NNN, sizeof(double));
        phi_dot[f]     = calloc(NNN, sizeof(double));
        phi_acc[f]     = calloc(NNN, sizeof(double));
        phi_acc_old[f] = calloc(NNN, sizeof(double));
        theta[f]       = calloc(NNN, sizeof(double));
        theta_dot[f]   = calloc(NNN, sizeof(double));
        theta_acc[f]   = calloc(NNN, sizeof(double));
        theta_acc_old[f] = calloc(NNN, sizeof(double));
    }
}

void free_fields(void) {
    for (int f = 0; f < 3; f++) {
        free(phi[f]); free(phi_dot[f]); free(phi_acc[f]); free(phi_acc_old[f]);
        free(theta[f]); free(theta_dot[f]); free(theta_acc[f]); free(theta_acc_old[f]);
    }
}

void zero_fields(void) {
    for (int f = 0; f < 3; f++) {
        memset(phi[f],         0, NNN * sizeof(double));
        memset(phi_dot[f],     0, NNN * sizeof(double));
        memset(phi_acc[f],     0, NNN * sizeof(double));
        memset(phi_acc_old[f], 0, NNN * sizeof(double));
        memset(theta[f],       0, NNN * sizeof(double));
        memset(theta_dot[f],   0, NNN * sizeof(double));
        memset(theta_acc[f],   0, NNN * sizeof(double));
        memset(theta_acc_old[f], 0, NNN * sizeof(double));
    }
}

/* 7-point Laplacian (periodic) */
static inline double laplacian(double *u, int i, int j, int k) {
    double inv_dx2 = 1.0 / (DX * DX);
    return inv_dx2 * (
        u[IDX(pw(i+1), j, k)] + u[IDX(pw(i-1), j, k)]
      + u[IDX(i, pw(j+1), k)] + u[IDX(i, pw(j-1), k)]
      + u[IDX(i, j, pw(k+1))] + u[IDX(i, j, pw(k-1))]
      - 6.0 * u[IDX(i, j, k)]
    );
}

/* Central difference partial derivatives */
static inline double ddx(double *u, int i, int j, int k) {
    return (u[IDX(pw(i+1), j, k)] - u[IDX(pw(i-1), j, k)]) / (2.0 * DX);
}
static inline double ddy(double *u, int i, int j, int k) {
    return (u[IDX(i, pw(j+1), k)] - u[IDX(i, pw(j-1), k)]) / (2.0 * DX);
}
static inline double ddz(double *u, int i, int j, int k) {
    return (u[IDX(i, j, pw(k+1))] - u[IDX(i, j, pw(k-1))]) / (2.0 * DX);
}

/* Damping coefficient at position (absorbing boundary at all edges) */
static inline double damp_coeff(double x, double y, double z) {
    double d = 0.0;
    double ax = fabs(x), ay = fabs(y), az = fabs(z);
    if (ax > DAMP_START) {
        double s = (ax - DAMP_START) / (L - DAMP_START);
        if (s > 1.0) s = 1.0;
        if (s > d) d = s;
    }
    if (ay > DAMP_START) {
        double s = (ay - DAMP_START) / (L - DAMP_START);
        if (s > 1.0) s = 1.0;
        if (s > d) d = s;
    }
    if (az > DAMP_START) {
        double s = (az - DAMP_START) / (L - DAMP_START);
        if (s > 1.0) s = 1.0;
        if (s > d) d = s;
    }
    return DAMP_MAX * d * d;  /* quadratic ramp */
}

/* Current eta (set per run) */
static double eta_val = 0.0;

/* Background phi energy (computed at t=0 before pulse, for subtraction) */
static double E_phi_bg = 0.0;

/* Compute accelerations with absorbing damping */
void compute_accelerations(void) {
    #pragma omp parallel for collapse(3) schedule(static)
    for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
    for (int k = 0; k < N; k++) {
        int idx = IDX(i, j, k);
        double x = coord(i), y = coord(j), z = coord(k);
        double gamma = damp_coeff(x, y, z);

        /* Product P = phi_0 * phi_1 * phi_2 */
        double p0 = phi[0][idx], p1 = phi[1][idx], p2 = phi[2][idx];
        double P = p0 * p1 * p2;
        double P2 = P * P;
        double denom = 1.0 + KAPPA * P2;

        /* V'(P) = mu * P / (1 + kappa*P²)² */
        double dVdP = MU * P / (denom * denom);

        /* dP/dphi_a */
        double dPdphi[3];
        dPdphi[0] = p1 * p2;
        dPdphi[1] = p0 * p2;
        dPdphi[2] = p0 * p1;

        /* curl(theta) */
        double curl_theta[3];
        curl_theta[0] = ddy(theta[2], i, j, k) - ddz(theta[1], i, j, k);
        curl_theta[1] = ddz(theta[0], i, j, k) - ddx(theta[2], i, j, k);
        curl_theta[2] = ddx(theta[1], i, j, k) - ddy(theta[0], i, j, k);

        /* curl(phi) */
        double curl_phi[3];
        curl_phi[0] = ddy(phi[2], i, j, k) - ddz(phi[1], i, j, k);
        curl_phi[1] = ddz(phi[0], i, j, k) - ddx(phi[2], i, j, k);
        curl_phi[2] = ddx(phi[1], i, j, k) - ddy(phi[0], i, j, k);

        /* phi accelerations */
        for (int a = 0; a < 3; a++) {
            phi_acc[a][idx] = laplacian(phi[a], i, j, k)
                            - M2 * phi[a][idx]
                            - dVdP * dPdphi[a]
                            + eta_val * curl_theta[a]
                            - gamma * phi_dot[a][idx];
        }

        /* theta accelerations */
        for (int a = 0; a < 3; a++) {
            theta_acc[a][idx] = laplacian(theta[a], i, j, k)
                              + eta_val * curl_phi[a]
                              - gamma * theta_dot[a][idx];
        }
    }
}

/* Initialize plane-wave pulse in theta_z + phi background */
void initialize(double v_g) {
    double k_bg = M_PI / L;  /* one half-wavelength across box */

    #pragma omp parallel for collapse(3)
    for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
    for (int k = 0; k < N; k++) {
        int idx = IDX(i, j, k);
        double x = coord(i);
        double z = coord(k);

        /* Phi background: standing wave in z, velocity=0 at t=0 */
        for (int a = 0; a < 3; a++) {
            phi[a][idx] = A_BG * cos(k_bg * z + 2.0 * M_PI * a / 3.0 + delta[a]);
            phi_dot[a][idx] = 0.0;
        }

        /* Theta_z plane-wave pulse */
        double dx_x = x - X0;
        double env = exp(-dx_x * dx_x / (2.0 * SIGMA_X * SIGMA_X));
        theta[2][idx] = AMP * env * sin(K0 * x);
        theta_dot[2][idx] = -v_g * AMP * env *
            (-dx_x / (SIGMA_X * SIGMA_X) * sin(K0 * x) + K0 * cos(K0 * x));
    }
}

/* Compute phi energy */
double compute_E_phi(void) {
    double E_phi = 0;
    double dV = DX * DX * DX;

    #pragma omp parallel for collapse(3) reduction(+:E_phi)
    for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
    for (int k = 0; k < N; k++) {
        int idx = IDX(i, j, k);
        for (int a = 0; a < 3; a++) {
            double pdot = phi_dot[a][idx];
            double gx = ddx(phi[a], i, j, k);
            double gy = ddy(phi[a], i, j, k);
            double gz = ddz(phi[a], i, j, k);
            E_phi += 0.5 * (pdot*pdot + gx*gx + gy*gy + gz*gz + M2 * phi[a][idx] * phi[a][idx]) * dV;
        }
        double P = phi[0][idx] * phi[1][idx] * phi[2][idx];
        double P2 = P * P;
        E_phi += 0.5 * MU * P2 / (1.0 + KAPPA * P2) * dV;
    }
    return E_phi;
}

/* Measure: project theta_z onto x-axis, find envelope peak position and width.
 *
 * We project theta_z^2 onto x to get P(x) = sum_{j,k} theta_z(i,j,k)^2.
 * For a modulated pulse theta_z = A*env*sin(k0*x), P(x) ~ env^2 * sin^2(k0*x).
 * To extract the ENVELOPE position and width, we smooth P(x) by averaging
 * over a carrier wavelength (Hilbert-like envelope extraction).
 */
void measure(double t, FILE *fp) {
    /* Step 1: Project theta_z² onto x-axis */
    double Px[N];
    memset(Px, 0, sizeof(Px));

    for (int i = 0; i < N; i++) {
        double s = 0.0;
        #pragma omp parallel for collapse(2) reduction(+:s)
        for (int j = 0; j < N; j++)
        for (int k = 0; k < N; k++) {
            s += theta[2][IDX(i,j,k)] * theta[2][IDX(i,j,k)];
        }
        Px[i] = s;
    }

    /* Step 2: Smooth P(x) to extract envelope squared.
     * Average over ~1 carrier wavelength (lambda = 2*pi/k0 = 3.14, ~8 grid pts).
     * Use a window of +/- 4 grid points. */
    int hw = (int)(M_PI / (K0 * DX) + 0.5);  /* half-wavelength in grid points */
    if (hw < 2) hw = 2;
    double Px_smooth[N];
    for (int i = 0; i < N; i++) {
        double s = 0.0;
        int cnt = 0;
        for (int di = -hw; di <= hw; di++) {
            int ii = pw(i + di);
            s += Px[ii];
            cnt++;
        }
        Px_smooth[i] = s / cnt;
    }

    /* Step 3: Find peak of smoothed P(x) */
    int i_peak = 0;
    double p_max = Px_smooth[0];
    for (int i = 1; i < N; i++) {
        if (Px_smooth[i] > p_max) {
            p_max = Px_smooth[i];
            i_peak = i;
        }
    }

    /* Step 4: Weighted moments of smoothed P(x) for Gaussian fit */
    double sum_w = 0, sum_wx = 0, sum_wx2 = 0;
    double x_peak_ref = coord(i_peak);
    for (int i = 0; i < N; i++) {
        double x = coord(i);
        double dx = x - x_peak_ref;
        /* Periodic unwrap */
        while (dx >  L) dx -= 2.0 * L;
        while (dx < -L) dx += 2.0 * L;
        sum_w   += Px_smooth[i];
        sum_wx  += dx * Px_smooth[i];
        sum_wx2 += dx * dx * Px_smooth[i];
    }

    double x_peak = x_peak_ref;
    double sigma_meas = SIGMA_X;
    /* Peak amplitude: sqrt of smoothed peak / (N*N) gives RMS-like amplitude */
    double peak_amp = sqrt(p_max / (N * N));

    if (sum_w > 0) {
        x_peak = x_peak_ref + sum_wx / sum_w;
        double var = sum_wx2 / sum_w - (sum_wx / sum_w) * (sum_wx / sum_w);
        /* For smoothed env^2, the variance = sigma_x^2 (envelope variance) */
        if (var > 0) sigma_meas = sqrt(var);
    }

    /* Step 5: Compute energies */
    double E_theta = 0;
    double dV = DX * DX * DX;

    #pragma omp parallel for collapse(3) reduction(+:E_theta)
    for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
    for (int k = 0; k < N; k++) {
        int idx = IDX(i, j, k);
        for (int a = 0; a < 3; a++) {
            double tdot = theta_dot[a][idx];
            double gx = ddx(theta[a], i, j, k);
            double gy = ddy(theta[a], i, j, k);
            double gz = ddz(theta[a], i, j, k);
            E_theta += 0.5 * (tdot*tdot + gx*gx + gy*gy + gz*gz) * dV;
        }
    }

    double E_phi = compute_E_phi();
    /* Subtract background phi energy to see pulse-induced phi excitation */
    double dE_phi = E_phi - E_phi_bg;

    fprintf(fp, "%.4f\t%.6f\t%.6f\t%.8e\t%.8e\t%.8e\t%.8e\n",
            t, x_peak, sigma_meas, peak_amp, E_theta, E_phi, dE_phi);

    printf("  t=%6.2f  x_pk=%8.4f  sig=%7.4f  amp=%.4e  E_th=%.4e  dE_ph=%.4e\n",
           t, x_peak, sigma_meas, peak_amp, E_theta, dE_phi);
}

/* Run one simulation */
void run_simulation(double eta, double v_g_pred, const char *outfile,
                    double *out_vg, double *out_sigma0, double *out_sigma_final,
                    double *out_ephi_frac) {
    eta_val = eta;
    zero_fields();
    initialize(v_g_pred);

    /* Compute background phi energy (before any evolution) */
    E_phi_bg = compute_E_phi();

    FILE *fp = fopen(outfile, "w");
    fprintf(fp, "t\tx_peak\tsigma\tpeak_amp\tE_theta\tE_phi\tdE_phi\n");

    printf("\n=== Run: eta=%.2f, predicted v_g=%.4f ===\n", eta, v_g_pred);
    printf("  Output: %s\n", outfile);
    printf("  Background E_phi = %.4e\n", E_phi_bg);

    compute_accelerations();

    double t = 0.0;
    double next_meas = 0.0;
    int n_steps = (int)(T_END / DT) + 1;

    measure(t, fp);
    next_meas += MEAS_DT;

    double t_start = omp_get_wtime();

    for (int step = 0; step < n_steps; step++) {
        t = (step + 1) * DT;

        /* Velocity Verlet step 1: positions */
        #pragma omp parallel for
        for (int idx = 0; idx < NNN; idx++) {
            for (int a = 0; a < 3; a++) {
                phi[a][idx]   += phi_dot[a][idx] * DT + 0.5 * phi_acc[a][idx] * DT * DT;
                theta[a][idx] += theta_dot[a][idx] * DT + 0.5 * theta_acc[a][idx] * DT * DT;
            }
        }

        /* Save old accelerations */
        for (int a = 0; a < 3; a++) {
            memcpy(phi_acc_old[a], phi_acc[a], NNN * sizeof(double));
            memcpy(theta_acc_old[a], theta_acc[a], NNN * sizeof(double));
        }

        compute_accelerations();

        /* Velocity Verlet step 2: velocities */
        #pragma omp parallel for
        for (int idx = 0; idx < NNN; idx++) {
            for (int a = 0; a < 3; a++) {
                phi_dot[a][idx]   += 0.5 * (phi_acc_old[a][idx] + phi_acc[a][idx]) * DT;
                theta_dot[a][idx] += 0.5 * (theta_acc_old[a][idx] + theta_acc[a][idx]) * DT;
            }
        }

        if (t >= next_meas - DT * 0.5) {
            measure(t, fp);
            next_meas += MEAS_DT;
        }
    }

    double elapsed = omp_get_wtime() - t_start;
    printf("  Completed in %.1f seconds\n", elapsed);
    fclose(fp);

    /* Read back results for analysis */
    fp = fopen(outfile, "r");
    char line[1024];
    (void)fgets(line, sizeof(line), fp); /* skip header */

    double t_vals[200], x_vals[200], s_vals[200], amp_vals[200];
    double e_theta[200], de_phi[200];
    int n = 0;
    while (fgets(line, sizeof(line), fp) && n < 200) {
        double e_phi_raw;
        sscanf(line, "%lf %lf %lf %lf %lf %lf %lf",
               &t_vals[n], &x_vals[n], &s_vals[n], &amp_vals[n],
               &e_theta[n], &e_phi_raw, &de_phi[n]);
        n++;
    }
    fclose(fp);

    if (n < 5) {
        printf("  WARNING: too few data points (%d)\n", n);
        *out_vg = 0; *out_sigma0 = 0; *out_sigma_final = 0; *out_ephi_frac = 0;
        return;
    }

    /* Unwrap x positions */
    double x_unwrap[200];
    x_unwrap[0] = x_vals[0];
    for (int i = 1; i < n; i++) {
        double dx = x_vals[i] - x_vals[i-1];
        while (dx >  L) dx -= 2.0 * L;
        while (dx < -L) dx += 2.0 * L;
        x_unwrap[i] = x_unwrap[i-1] + dx;
    }

    /* Linear fit for velocity.
     * Only use data where the pulse is away from boundaries.
     * Pulse starts at x0=-15, damping at |x|>17.5.
     * At v~1, pulse reaches x=17.5 at t~32.5.
     * Use t=2 to t=30 (safely before damping zone). */
    double t_fit_min = 2.0, t_fit_max = 28.0;
    double st = 0, sx = 0, stt = 0, stx = 0, sw = 0;
    for (int i = 0; i < n; i++) {
        if (t_vals[i] >= t_fit_min && t_vals[i] <= t_fit_max) {
            st += t_vals[i]; sx += x_unwrap[i];
            stt += t_vals[i] * t_vals[i]; stx += t_vals[i] * x_unwrap[i];
            sw += 1.0;
        }
    }
    double v_meas = 0;
    if (sw > 2 && (sw * stt - st * st) > 0)
        v_meas = (sw * stx - st * sx) / (sw * stt - st * st);

    /* Sigma at early/mid times (before damping corrupts) */
    *out_sigma0 = s_vals[0];
    /* Find sigma at t~25 */
    int i_mid = 0;
    for (int i = 0; i < n; i++) {
        if (t_vals[i] <= 25.0) i_mid = i;
    }
    *out_sigma_final = s_vals[i_mid];

    /* E_phi / E_theta: average over t=15-25 (steady state, before boundary) */
    double avg_frac = 0;
    int n_avg = 0;
    for (int i = 0; i < n; i++) {
        if (t_vals[i] >= 15.0 && t_vals[i] <= 25.0 && e_theta[i] > 0) {
            avg_frac += de_phi[i] / e_theta[i];
            n_avg++;
        }
    }
    if (n_avg > 0) avg_frac /= n_avg;

    *out_vg = v_meas;
    *out_ephi_frac = avg_frac;
}

int main(void) {
    printf("===============================================\n");
    printf("OQ5 Plane-Wave Pulse Propagation Test\n");
    printf("===============================================\n");
    printf("Grid: N=%d, L=%.1f, dx=%.4f, domain=[%.1f, %.1f]\n",
           N, L, DX, -L, L);
    printf("Carrier: k0=%.1f, wavelength=%.4f, pts/wavelength=%.1f\n",
           K0, 2.0*M_PI/K0, 2.0*M_PI/K0/DX);
    printf("Pulse: A=%.3f, sigma_x=%.1f, x0=%.1f\n", AMP, SIGMA_X, X0);
    printf("Background: A_bg=%.2f, phi with delta offsets\n", A_BG);
    printf("Physics: m^2=%.2f, mu=%.3f, kappa=%.1f\n", M2, MU, KAPPA);
    printf("Time: T=%.1f, dt=%.6f, meas_dt=%.1f\n", T_END, DT, MEAS_DT);
    printf("Absorbing boundary: starts at |x|>%.1f\n", DAMP_START);
    printf("Memory: %.1f MB total\n", NNN * 8.0 * 24 / 1e6);

    /* Lattice-corrected dispersion predictions.
     * 7-point Laplacian: k_eff = (2/dx)*sin(k*dx/2) replaces k in dispersion.
     * Group velocity: d(omega)/dk, where omega depends on k_eff(k). */
    double k_eff = 2.0 * sin(K0 * DX / 2.0) / DX;
    double k_eff_deriv = cos(K0 * DX / 2.0);  /* d(k_eff)/dk */

    /* Continuum predictions */
    double m4_4 = M2 * M2 / 4.0;
    double eta2_k2_cont = 0.25 * K0 * K0;
    double sq_cont = sqrt(m4_4 + eta2_k2_cont);
    double omega2_cont = K0 * K0 + M2 / 2.0 - sq_cont;
    double omega_cont = sqrt(omega2_cont);
    /* d(omega^2)/dk = 2k - eta^2*k/sqrt(...) */
    double domega2_dk_cont = 2.0*K0 - 0.25*K0/sq_cont;
    double vg_cont = domega2_dk_cont / (2.0 * omega_cont);

    /* Lattice predictions: use k_eff everywhere k appears */
    double eta2_keff2 = 0.25 * k_eff * k_eff;
    double sq_lat = sqrt(m4_4 + eta2_keff2);
    double omega2_lat = k_eff * k_eff + M2 / 2.0 - sq_lat;
    double omega_lat = sqrt(omega2_lat);
    /* d(omega^2)/dk = d(omega^2)/d(k_eff) * d(k_eff)/dk */
    double domega2_dkeff = 2.0*k_eff - 0.25*k_eff/sq_lat;
    double vg_lat = domega2_dkeff * k_eff_deriv / (2.0 * omega_lat);

    /* eta=0 lattice: omega = k_eff, v_g = d(k_eff)/dk = cos(k*dx/2) */
    double vg0_lat = k_eff_deriv;

    printf("\nLattice dispersion (k0=%.1f, dx=%.4f):\n", K0, DX);
    printf("  k_eff = %.6f (vs k0=%.6f), d(k_eff)/dk = %.6f\n",
           k_eff, K0, k_eff_deriv);
    printf("  eta=0:   v_g_lattice = %.6f (continuum: 1.000)\n", vg0_lat);
    printf("  eta=0.5: omega_lat = %.6f, v_g_lattice = %.6f (continuum: %.6f)\n",
           omega_lat, vg_lat, vg_cont);
    double v_g_pred = vg_lat;

    alloc_fields();

    double vg0, sigma0_0, sigmaf_0, efrac_0;
    double vg1, sigma0_1, sigmaf_1, efrac_1;

    /* Run 1: eta=0 (control) — use lattice v_g for clean rightward-only pulse */
    run_simulation(0.0, vg0_lat, "planewave_eta0.tsv",
                   &vg0, &sigma0_0, &sigmaf_0, &efrac_0);

    /* Run 2: eta=0.5 (coupled) — use lattice-predicted v_g */
    run_simulation(0.5, vg_lat, "planewave_eta05.tsv",
                   &vg1, &sigma0_1, &sigmaf_1, &efrac_1);

    /* Summary */
    double pred_efrac = 0.25 * k_eff * k_eff / ((k_eff * k_eff + M2) * (k_eff * k_eff + M2));

    printf("\n");
    printf("=============================================\n");
    printf("       OQ5 PLANE-WAVE PULSE SUMMARY\n");
    printf("=============================================\n");
    printf("\n--- eta=0.0 (control, uncoupled) ---\n");
    printf("  Measured v_group:    %.6f  (fit t=2..25)\n", vg0);
    printf("  Continuum predicted: 1.000000\n");
    printf("  Lattice predicted:   %.6f\n", vg0_lat);
    printf("  Error vs lattice:    %.4f (%.2f%%)\n",
           fabs(vg0 - vg0_lat), fabs(vg0 - vg0_lat) / vg0_lat * 100.0);
    printf("  Sigma: %.4f -> %.4f (%.1f%% change, t=0..25)\n",
           sigma0_0, sigmaf_0,
           sigma0_0 > 0 ? 100.0 * (sigmaf_0 - sigma0_0) / sigma0_0 : 0.0);
    printf("  dE_phi/E_theta:      %.2e (should be ~0)\n", efrac_0);

    printf("\n--- eta=0.5 (coupled) ---\n");
    printf("  Measured v_group:    %.6f  (fit t=2..25)\n", vg1);
    printf("  Continuum predicted: %.6f\n", vg_cont);
    printf("  Lattice predicted:   %.6f\n", vg_lat);
    printf("  Error vs lattice:    %.4f (%.2f%%)\n",
           fabs(vg1 - vg_lat), fabs(vg1 - vg_lat) / vg_lat * 100.0);
    printf("  Sigma: %.4f -> %.4f (%.1f%% change, t=0..25)\n",
           sigma0_1, sigmaf_1,
           sigma0_1 > 0 ? 100.0 * (sigmaf_1 - sigma0_1) / sigma0_1 : 0.0);
    printf("  dE_phi/E_theta:      %.4f (predicted ~%.4f)\n", efrac_1, pred_efrac);

    printf("\n--- Coupling slowdown ---\n");
    double slowdown_meas = vg0 > 0 ? (vg0 - vg1) / vg0 : 0;
    double slowdown_lat = vg0_lat > 0 ? (vg0_lat - vg_lat) / vg0_lat : 0;
    printf("  Measured:  delta_v/v = %.4f (%.2f%%)\n", slowdown_meas, slowdown_meas * 100.0);
    printf("  Predicted: delta_v/v = %.4f (%.2f%%)\n", slowdown_lat, slowdown_lat * 100.0);

    printf("\n--- PASS/FAIL ---\n");
    int pass0 = fabs(vg0 - vg0_lat) / vg0_lat < 0.02;
    int pass1 = fabs(vg1 - vg_lat) / vg_lat < 0.02;
    printf("  eta=0 velocity:   %s (within 2%% of lattice prediction %.4f)\n",
           pass0 ? "PASS" : "FAIL", vg0_lat);
    printf("  eta=0.5 velocity: %s (within 2%% of lattice prediction %.4f)\n",
           pass1 ? "PASS" : "FAIL", vg_lat);
    printf("=============================================\n");

    free_fields();
    return 0;
}
