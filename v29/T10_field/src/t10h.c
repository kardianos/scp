/*  t10h.c — Two Charged Complex Braids: Force Measurement
 *
 *  From T10G: the complex braid carries U(1) charge Q~3-14.
 *  Question: do two charged braids interact via a 1/r Coulomb force?
 *
 *  Setup: Two complex braids at separation D, with charges of same or opposite sign.
 *  Measure the force (change in separation) over time.
 *
 *  For comparison also measure the REAL-field case (no charge, just Yukawa repulsion).
 *
 *  Charge sign control: multiply the imaginary components of one braid by -1
 *  to flip its charge.
 *
 *  Build: gcc -O3 -fopenmp -o t10h src/t10h.c -lm
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>

#define NREAL 6
#define PI 3.14159265358979323846

typedef struct {
    double *phi[NREAL];
    double *vel[NREAL];
    double *acc[NREAL];
    int N;
    double L, dx, dt;
} CGrid;

static CGrid *cgrid_alloc(int N, double L) {
    CGrid *g = calloc(1, sizeof(CGrid));
    int N3 = N*N*N;
    for (int a = 0; a < NREAL; a++) {
        g->phi[a] = calloc(N3, sizeof(double));
        g->vel[a] = calloc(N3, sizeof(double));
        g->acc[a] = calloc(N3, sizeof(double));
    }
    g->N = N; g->L = L;
    g->dx = 2.0*L/(N-1);
    g->dt = 0.20 * g->dx;
    return g;
}

static void cgrid_free(CGrid *g) {
    for (int a = 0; a < NREAL; a++) {
        free(g->phi[a]); free(g->vel[a]); free(g->acc[a]);
    }
    free(g);
}

static void ccompute_forces(CGrid *g, double mu, double kappa, double mass2) {
    int N = g->N, NN = N*N, N3 = N*N*N;
    double idx2 = 1.0 / (g->dx * g->dx);

    #pragma omp parallel for schedule(static)
    for (int idx = 0; idx < N3; idx++) {
        int i = idx / NN, j = (idx / N) % N, k = idx % N;
        if (i < 1 || i >= N-1 || j < 1 || j >= N-1) {
            for (int a = 0; a < NREAL; a++) g->acc[a][idx] = 0;
            continue;
        }
        int kp = (k+1)%N, km = (k-1+N)%N;
        int idx_kp = i*NN + j*N + kp, idx_km = i*NN + j*N + km;

        double r0 = g->phi[0][idx], i0 = g->phi[1][idx];
        double r1 = g->phi[2][idx], i1 = g->phi[3][idx];
        double r2 = g->phi[4][idx], i2 = g->phi[5][idx];

        /* P = ψ₀ψ₁ψ₂ (complex) */
        double t_re = r0*r1 - i0*i1, t_im = r0*i1 + i0*r1;
        double reP = t_re*r2 - t_im*i2, imP = t_re*i2 + t_im*r2;
        double P2 = reP*reP + imP*imP;
        double D = 1.0 + kappa*P2;
        double D2 = D*D;
        double dVdPr = mu * reP / D2;
        double dVdPi = mu * imP / D2;

        /* dP/dψ_a */
        double dP0_re = r1*r2-i1*i2, dP0_im = r1*i2+i1*r2;
        double dP1_re = r0*r2-i0*i2, dP1_im = r0*i2+i0*r2;
        double dP2_re = r0*r1-i0*i1, dP2_im = r0*i1+i0*r1;

        double fV[NREAL];
        fV[0] = -(dVdPr*dP0_re + dVdPi*dP0_im);
        fV[1] = -(dVdPi*dP0_re - dVdPr*dP0_im);
        fV[2] = -(dVdPr*dP1_re + dVdPi*dP1_im);
        fV[3] = -(dVdPi*dP1_re - dVdPr*dP1_im);
        fV[4] = -(dVdPr*dP2_re + dVdPi*dP2_im);
        fV[5] = -(dVdPi*dP2_re - dVdPr*dP2_im);

        for (int a = 0; a < NREAL; a++) {
            double lap = (g->phi[a][idx+NN] + g->phi[a][idx-NN]
                        + g->phi[a][idx+N]  + g->phi[a][idx-N]
                        + g->phi[a][idx_kp] + g->phi[a][idx_km]
                        - 6.0*g->phi[a][idx]) * idx2;
            g->acc[a][idx] = lap - mass2*g->phi[a][idx] + fV[a];
        }
    }
}

static void cverlet_step(CGrid *g, double mu, double kappa, double mass2) {
    int N3 = g->N*g->N*g->N;
    double hdt = 0.5*g->dt;
    for (int a = 0; a < NREAL; a++)
        for (int idx = 0; idx < N3; idx++)
            g->vel[a][idx] += hdt * g->acc[a][idx];
    for (int a = 0; a < NREAL; a++)
        for (int idx = 0; idx < N3; idx++)
            g->phi[a][idx] += g->dt * g->vel[a][idx];
    ccompute_forces(g, mu, kappa, mass2);
    for (int a = 0; a < NREAL; a++)
        for (int idx = 0; idx < N3; idx++)
            g->vel[a][idx] += hdt * g->acc[a][idx];
}

static void capply_damping_xy(CGrid *g) {
    int N = g->N, NN = N*N;
    double dx = g->dx, L = g->L;
    double r_start = 0.70*L, r_end = 0.95*L;
    double inv_dr = 1.0/(r_end-r_start+1e-30);
    for (int i = 0; i < N; i++) {
        double x = -L + i*dx;
        for (int j = 0; j < N; j++) {
            double y = -L + j*dx;
            double rp = sqrt(x*x+y*y);
            if (rp <= r_start) continue;
            double f = (rp-r_start)*inv_dr;
            if (f > 1.0) f = 1.0;
            double damp = 1.0-0.98*f*f;
            for (int kk = 0; kk < N; kk++) {
                int idx = i*NN+j*N+kk;
                for (int a = 0; a < NREAL; a++) {
                    g->phi[a][idx] *= damp;
                    g->vel[a][idx] *= damp;
                }
            }
        }
    }
}

/* Initialize a single complex braid at position (x_off, 0, 0) */
static void init_complex_braid(CGrid *g, double x_off, double phase_twist,
                                double charge_sign, double mass2) {
    int N = g->N, NN = N*N;
    double dx = g->dx, L = g->L;
    double A = 0.8, R_tube = 3.0, ellip = 0.33;
    double k = PI / L;
    double omega = sqrt(k*k + mass2);
    double sx = 1.0+ellip, sy = 1.0-ellip;
    double inv_2R2 = 1.0/(2.0*R_tube*R_tube);
    double delta[3] = {0.0, 3.53*0.85, 4.92*0.85};

    for (int i = 0; i < N; i++) {
        double x = -L + i*dx - x_off;
        for (int j = 0; j < N; j++) {
            double y = -L + j*dx;
            for (int kk = 0; kk < N; kk++) {
                double z = -L + kk*dx;
                int idx = i*NN + j*N + kk;
                for (int c = 0; c < 3; c++) {
                    double xr = x/sx, yr = y/sy;
                    double r2e = xr*xr + yr*yr;
                    double env = A * exp(-r2e * inv_2R2);
                    double ph = k*z + delta[c];
                    /* Add to existing field (superposition) */
                    g->phi[2*c][idx] += env * cos(ph);
                    g->vel[2*c][idx] += omega * env * sin(ph);
                    g->phi[2*c+1][idx] += charge_sign * phase_twist * env * sin(ph);
                    g->vel[2*c+1][idx] += charge_sign * phase_twist * omega * env * cos(ph);
                }
            }
        }
    }
}

/* Find center of mass of energy density in the left/right half */
static double find_x_centroid(CGrid *g, double xmin, double xmax) {
    int N = g->N, NN = N*N;
    double dx = g->dx, L = g->L;
    double wx = 0, wtot = 0;
    for (int i = 1; i < N-1; i++) {
        double x = -L + i*dx;
        if (x < xmin || x > xmax) continue;
        for (int j = 1; j < N-1; j++) {
            double y = -L + j*dx;
            for (int k = 0; k < N; k++) {
                int idx = i*NN + j*N + k;
                double rho = 0;
                for (int a = 0; a < NREAL; a++)
                    rho += g->phi[a][idx]*g->phi[a][idx];
                wx += rho * x;
                wtot += rho;
            }
        }
    }
    return wx / (wtot + 1e-30);
}

int main(int argc, char **argv) {
    int N = 96;
    double L = 30.0; /* large box for two braids */
    double D = 15.0; /* separation */
    double T_equil = 100.0;
    double T_track = 200.0;

    double mu = -41.3, kappa = 50.0, mass2 = 2.25;
    double phase_twist = 0.3;

    printf("=== T10H: Two Charged Complex Braids ===\n\n");
    printf("Grid: N=%d, L=%.1f, D=%.1f\n", N, L, D);
    printf("mu=%.1f, kappa=%.1f, m²=%.2f\n\n", mu, kappa, mass2);

    FILE *fout = fopen("data/t10h_separation.tsv", "w");
    fprintf(fout, "config\tt\tx_left\tx_right\tseparation\n");

    /* Three test cases:
     * 1. Same charge (both +)
     * 2. Opposite charge (+ and -)
     * 3. Real only (no imaginary parts = no charge)
     */
    const char *configs[] = {"same_charge", "opposite_charge", "no_charge"};
    double charge_signs_R[] = {+1.0, -1.0, 0.0};  /* right braid charge sign */
    int nconfigs = 3;

    for (int ic = 0; ic < nconfigs; ic++) {
        printf("=== %s ===\n", configs[ic]);
        fflush(stdout);

        CGrid *g = cgrid_alloc(N, L);
        int N3 = N*N*N;
        for (int a = 0; a < NREAL; a++) {
            memset(g->phi[a], 0, N3*sizeof(double));
            memset(g->vel[a], 0, N3*sizeof(double));
        }

        /* Left braid at -D/2, charge +1 */
        init_complex_braid(g, -D/2, phase_twist, +1.0, mass2);
        /* Right braid at +D/2, charge given by config */
        double cs = charge_signs_R[ic];
        double pt = (ic == 2) ? 0.0 : phase_twist; /* no imaginary for real-only */
        init_complex_braid(g, +D/2, pt, cs, mass2);

        ccompute_forces(g, mu, kappa, mass2);

        /* Brief equilibration with damping */
        int equil_steps = (int)(T_equil / g->dt);
        for (int s = 0; s < equil_steps; s++) {
            cverlet_step(g, mu, kappa, mass2);
            capply_damping_xy(g);
        }

        /* Track separation */
        double x_left_0 = find_x_centroid(g, -L, 0);
        double x_right_0 = find_x_centroid(g, 0, L);
        double D0 = x_right_0 - x_left_0;
        printf("  Initial: left=%.2f, right=%.2f, D=%.2f\n", x_left_0, x_right_0, D0);

        int track_steps = (int)(T_track / g->dt);
        int meas_every = track_steps / 100;
        if (meas_every < 1) meas_every = 1;

        double last_D = D0;
        for (int s = 0; s < track_steps; s++) {
            cverlet_step(g, mu, kappa, mass2);
            capply_damping_xy(g);

            if (s % meas_every == 0) {
                double xl = find_x_centroid(g, -L, 0);
                double xr = find_x_centroid(g, 0, L);
                double sep = xr - xl;
                fprintf(fout, "%s\t%.3f\t%.4f\t%.4f\t%.4f\n",
                        configs[ic], (s+equil_steps)*g->dt, xl, xr, sep);
                last_D = sep;
            }
        }

        double dD = last_D - D0;
        printf("  Final: D=%.2f, ΔD=%.4f", last_D, dD);
        if (dD > 0.1) printf(" → REPULSION\n");
        else if (dD < -0.1) printf(" → ATTRACTION\n");
        else printf(" → no significant force\n");

        printf("\n");
        cgrid_free(g);
    }

    fclose(fout);

    printf("=== T10H Complete ===\n");
    printf("Data in data/t10h_separation.tsv\n");
    return 0;
}
