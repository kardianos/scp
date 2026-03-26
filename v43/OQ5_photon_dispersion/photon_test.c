/*
 * photon_test.c — Test pure θ wave packet propagation in empty Cosserat field
 *
 * Physics: 6-field Cosserat with φ=0 background.
 *   ∂²φ_a/∂t² = ∇²φ_a - m²φ_a - ∂V/∂φ_a + η curl(θ)_a
 *   ∂²θ_a/∂t² = ∇²θ_a                      + η curl(φ)_a
 *
 * With φ=0 everywhere and V'(0)=0, the θ equation reduces to □θ=0 (massless wave).
 * The φ equation gets a source η curl(θ), but with φ=0 initially and small θ,
 * any φ excitation should be negligible.
 *
 * Test: Gaussian θ_z pulse propagating in +x direction should move at c=1,
 * not disperse, and not excite φ.
 *
 * Build: gcc -O3 -march=native -fopenmp -o photon_test photon_test.c -lm
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>

/* Grid parameters */
#define N 128
#define L 40.0
#define DX (L / N)
#define NNN (N * N * N)

/* Physics parameters */
#define M2   2.25
#define MU   (-41.345)
#define KAPPA 50.0
#define ETA  0.5

/* Time stepping */
#define DT   (0.12 * DX)
#define T_END 60.0
#define MEAS_INTERVAL 1.0

/* Initial condition */
#define AMP   0.01
#define SIGMA 3.0
#define X0    (-15.0)

/* 3D index */
#define IDX(i,j,k) ((i)*N*N + (j)*N + (k))

/* Periodic wrap */
static inline int pw(int i) {
    return ((i % N) + N) % N;
}

/* Grid coordinate (centered at 0) */
static inline double coord(int i) {
    return -L/2.0 + (i + 0.5) * DX;
}

/* Periodic distance (for centroid tracking) */
static inline double periodic_dist(double x, double ref) {
    double d = x - ref;
    while (d > L/2.0)  d -= L;
    while (d < -L/2.0) d += L;
    return d;
}

/* Field arrays: phi[3][NNN], theta[3][NNN], velocities, accelerations */
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

/* 7-point Laplacian with periodic BCs */
static inline double laplacian(double *u, int i, int j, int k) {
    double inv_dx2 = 1.0 / (DX * DX);
    return inv_dx2 * (
        u[IDX(pw(i+1), j, k)] + u[IDX(pw(i-1), j, k)]
      + u[IDX(i, pw(j+1), k)] + u[IDX(i, pw(j-1), k)]
      + u[IDX(i, j, pw(k+1))] + u[IDX(i, j, pw(k-1))]
      - 6.0 * u[IDX(i, j, k)]
    );
}

/* Partial derivatives (central difference, periodic) */
static inline double ddx(double *u, int i, int j, int k) {
    return (u[IDX(pw(i+1), j, k)] - u[IDX(pw(i-1), j, k)]) / (2.0 * DX);
}
static inline double ddy(double *u, int i, int j, int k) {
    return (u[IDX(i, pw(j+1), k)] - u[IDX(i, pw(j-1), k)]) / (2.0 * DX);
}
static inline double ddz(double *u, int i, int j, int k) {
    return (u[IDX(i, j, pw(k+1))] - u[IDX(i, j, pw(k-1))]) / (2.0 * DX);
}

/* Compute accelerations for all fields */
void compute_accelerations(void) {
    double inv_dx2 = 1.0 / (DX * DX);

    #pragma omp parallel for collapse(3) schedule(static)
    for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
    for (int k = 0; k < N; k++) {
        int idx = IDX(i, j, k);

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

        /* curl(theta)_a = epsilon_{abc} d_b theta_c */
        double curl_theta[3];
        curl_theta[0] = ddy(theta[2], i, j, k) - ddz(theta[1], i, j, k);
        curl_theta[1] = ddz(theta[0], i, j, k) - ddx(theta[2], i, j, k);
        curl_theta[2] = ddx(theta[1], i, j, k) - ddy(theta[0], i, j, k);

        /* curl(phi)_a */
        double curl_phi[3];
        curl_phi[0] = ddy(phi[2], i, j, k) - ddz(phi[1], i, j, k);
        curl_phi[1] = ddz(phi[0], i, j, k) - ddx(phi[2], i, j, k);
        curl_phi[2] = ddx(phi[1], i, j, k) - ddy(phi[0], i, j, k);

        /* phi accelerations */
        for (int a = 0; a < 3; a++) {
            phi_acc[a][idx] = laplacian(phi[a], i, j, k)
                            - M2 * phi[a][idx]
                            - dVdP * dPdphi[a]
                            + ETA * curl_theta[a];
        }

        /* theta accelerations */
        for (int a = 0; a < 3; a++) {
            theta_acc[a][idx] = laplacian(theta[a], i, j, k)
                              + ETA * curl_phi[a];
        }
    }
}

/* Initialize: Gaussian θ_z pulse moving in +x */
void initialize(void) {
    #pragma omp parallel for collapse(3)
    for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
    for (int k = 0; k < N; k++) {
        int idx = IDX(i, j, k);
        double x = coord(i);
        double y = coord(j);
        double z = coord(k);

        double r2 = (x - X0)*(x - X0) + y*y + z*z;
        double gauss = AMP * exp(-r2 / (2.0 * SIGMA * SIGMA));

        /* theta_z = Gaussian */
        theta[2][idx] = gauss;

        /* d(theta_z)/dt for rightward-moving pulse at c=1:
         * For f(x - ct), df/dt = -c * df/dx = -c * f * (-(x-x0)/sigma^2)
         * = (x-x0)/sigma^2 * f   (with c=1)
         */
        theta_dot[2][idx] = gauss * (x - X0) / (SIGMA * SIGMA);
    }
}

/* Measure diagnostics */
void measure(double t, FILE *fp) {
    double sum_tz2 = 0, sum_x_tz2 = 0, sum_x2_tz2 = 0;
    double peak_amp = 0;
    double E_theta = 0, E_phi = 0;
    double dV = DX * DX * DX;

    /* First pass: find rough centroid for periodic unwrapping */
    /* Use sin/cos method for periodic centroid */
    double sin_sum = 0, cos_sum = 0, w_sum = 0;

    #pragma omp parallel for collapse(3) reduction(+:sin_sum,cos_sum,w_sum)
    for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
    for (int k = 0; k < N; k++) {
        int idx = IDX(i, j, k);
        double w = theta[2][idx] * theta[2][idx];
        double angle = 2.0 * M_PI * coord(i) / L;
        sin_sum += w * sin(angle);
        cos_sum += w * cos(angle);
        w_sum += w;
    }
    double x_cent = L / (2.0 * M_PI) * atan2(sin_sum, cos_sum);

    /* Second pass: compute width, energy, peak */
    #pragma omp parallel for collapse(3) reduction(+:sum_tz2,sum_x_tz2,sum_x2_tz2,E_theta,E_phi) reduction(max:peak_amp)
    for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
    for (int k = 0; k < N; k++) {
        int idx = IDX(i, j, k);
        double x = coord(i);

        /* θ_z² weighted moments */
        double tz2 = theta[2][idx] * theta[2][idx];
        double dx_p = periodic_dist(x, x_cent);
        sum_tz2 += tz2;
        sum_x_tz2 += dx_p * tz2;
        sum_x2_tz2 += dx_p * dx_p * tz2;

        double aval = fabs(theta[2][idx]);
        if (aval > peak_amp) peak_amp = aval;

        /* θ energy: ½[(∂θ/∂t)² + |∇θ|²] */
        for (int a = 0; a < 3; a++) {
            double tdot = theta_dot[a][idx];
            double gx = ddx(theta[a], i, j, k);
            double gy = ddy(theta[a], i, j, k);
            double gz = ddz(theta[a], i, j, k);
            E_theta += 0.5 * (tdot*tdot + gx*gx + gy*gy + gz*gz) * dV;
        }

        /* φ energy: ½[(∂φ/∂t)² + |∇φ|² + m²φ²] + V */
        for (int a = 0; a < 3; a++) {
            double pdot = phi_dot[a][idx];
            double gx = ddx(phi[a], i, j, k);
            double gy = ddy(phi[a], i, j, k);
            double gz = ddz(phi[a], i, j, k);
            E_phi += 0.5 * (pdot*pdot + gx*gx + gy*gy + gz*gz + M2 * phi[a][idx]*phi[a][idx]) * dV;
        }
        /* V(P) contribution */
        double P = phi[0][idx] * phi[1][idx] * phi[2][idx];
        double P2 = P * P;
        E_phi += 0.5 * MU * P2 / (1.0 + KAPPA * P2) * dV;
    }

    /* Refine centroid */
    if (sum_tz2 > 0) {
        x_cent += sum_x_tz2 / sum_tz2;
    }

    double sigma_meas = 0;
    if (sum_tz2 > 0) {
        /* Recompute width around refined centroid */
        double s2 = 0, sw = 0;
        #pragma omp parallel for collapse(3) reduction(+:s2,sw)
        for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
        for (int k = 0; k < N; k++) {
            int idx = IDX(i, j, k);
            double x = coord(i);
            double tz2 = theta[2][idx] * theta[2][idx];
            double dx_p = periodic_dist(x, x_cent);
            s2 += dx_p * dx_p * tz2;
            sw += tz2;
        }
        sigma_meas = sqrt(s2 / sw);
    }

    fprintf(fp, "%.4f\t%.6f\t%.6f\t%.8e\t%.8e\t%.8e\n",
            t, x_cent, sigma_meas, peak_amp, E_theta, E_phi);

    printf("t=%6.2f  x_c=%8.4f  σ=%7.4f  peak=%.4e  E_θ=%.4e  E_φ=%.4e\n",
           t, x_cent, sigma_meas, peak_amp, E_theta, E_phi);
}

int main(void) {
    printf("Cosserat θ-wave photon test\n");
    printf("N=%d, L=%.1f, dx=%.4f, dt=%.6f\n", N, L, DX, DT);
    printf("m²=%.2f, μ=%.3f, κ=%.1f, η=%.2f\n", M2, MU, KAPPA, ETA);
    printf("Pulse: A=%.3f, σ=%.1f, x0=%.1f\n", AMP, SIGMA, X0);
    printf("T_end=%.1f, steps=%d\n\n", T_END, (int)(T_END / DT));

    alloc_fields();
    initialize();

    FILE *fp = fopen("photon_test.tsv", "w");
    fprintf(fp, "t\tx_centroid\tsigma\tpeak_amplitude\tE_theta\tE_phi\n");

    /* Initial acceleration */
    compute_accelerations();

    /* Velocity Verlet integration */
    double t = 0.0;
    double next_meas = 0.0;
    int step = 0;
    int n_steps = (int)(T_END / DT) + 1;

    measure(t, fp);
    next_meas += MEAS_INTERVAL;

    double t_start = omp_get_wtime();

    for (step = 0; step < n_steps; step++) {
        t = (step + 1) * DT;

        /* Velocity Verlet step 1: update positions */
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

        /* Compute new accelerations */
        compute_accelerations();

        /* Velocity Verlet step 2: update velocities */
        #pragma omp parallel for
        for (int idx = 0; idx < NNN; idx++) {
            for (int a = 0; a < 3; a++) {
                phi_dot[a][idx]   += 0.5 * (phi_acc_old[a][idx] + phi_acc[a][idx]) * DT;
                theta_dot[a][idx] += 0.5 * (theta_acc_old[a][idx] + theta_acc[a][idx]) * DT;
            }
        }

        /* Measurement */
        if (t >= next_meas - DT * 0.5) {
            measure(t, fp);
            next_meas += MEAS_INTERVAL;
        }
    }

    double elapsed = omp_get_wtime() - t_start;
    printf("\nSimulation completed in %.1f seconds\n", elapsed);

    fclose(fp);

    /* Summary analysis */
    fp = fopen("photon_test.tsv", "r");
    char line[1024];
    fgets(line, sizeof(line), fp); /* skip header */

    double t_vals[200], x_vals[200], s_vals[200], e_theta[200], e_phi[200];
    int n = 0;
    while (fgets(line, sizeof(line), fp) && n < 200) {
        double pk;
        sscanf(line, "%lf %lf %lf %lf %lf %lf",
               &t_vals[n], &x_vals[n], &s_vals[n], &pk, &e_theta[n], &e_phi[n]);
        n++;
    }
    fclose(fp);

    if (n >= 10) {
        /* Linear fit for velocity: use middle portion to avoid edge effects */
        int i0 = 2, i1 = n - 3;
        /* Unwrap x positions for velocity fit */
        double x_unwrap[200];
        x_unwrap[0] = x_vals[0];
        for (int i = 1; i < n; i++) {
            double dx = x_vals[i] - x_vals[i-1];
            while (dx > L/2) dx -= L;
            while (dx < -L/2) dx += L;
            x_unwrap[i] = x_unwrap[i-1] + dx;
        }

        /* Least squares for velocity */
        double st = 0, sx = 0, stt = 0, stx = 0, sw = 0;
        for (int i = i0; i <= i1; i++) {
            st += t_vals[i]; sx += x_unwrap[i];
            stt += t_vals[i]*t_vals[i]; stx += t_vals[i]*x_unwrap[i];
            sw += 1.0;
        }
        double v_group = (sw*stx - st*sx) / (sw*stt - st*st);

        /* Dispersion: sigma change rate */
        double ds_dt = (s_vals[i1] - s_vals[i0]) / (t_vals[i1] - t_vals[i0]);

        /* Energy conservation */
        double E0 = e_theta[0];
        double E_final = e_theta[n-1];
        double dE_rel = (E_final - E0) / E0;

        /* Max phi energy */
        double max_ephi = 0;
        for (int i = 0; i < n; i++)
            if (e_phi[i] > max_ephi) max_ephi = e_phi[i];

        printf("\n========== SUMMARY ==========\n");
        printf("Measured group velocity: %.6f  (expected 1.0, error=%.2e)\n",
               v_group, fabs(v_group - 1.0));
        printf("Dispersion rate dσ/dt:  %.6f  (expected 0.0)\n", ds_dt);
        printf("Initial σ:              %.4f\n", s_vals[0]);
        printf("Final σ:                %.4f\n", s_vals[n-1]);
        printf("σ change:               %.2f%%\n", 100.0*(s_vals[n-1]-s_vals[0])/s_vals[0]);
        printf("E_θ conservation:       %.2e relative\n", dE_rel);
        printf("Max E_φ:                %.2e (coupling leakage)\n", max_ephi);
        printf("E_φ/E_θ ratio:          %.2e\n", max_ephi / E0);
        printf("=============================\n");
    }

    free_fields();
    return 0;
}
