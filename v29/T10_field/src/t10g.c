/*  t10g.c — Complex Field Upgrade: First Test
 *
 *  Upgrade from 3 real scalars to 3 complex scalars:
 *    ψ_a = (φ_a + i χ_a) / √2
 *
 *  Lagrangian:
 *    L = Σ_a |∂ψ_a|² - m²|ψ_a|² - (μ/2)|P|²/(1+κ|P|²)
 *    P = ψ₀ ψ₁ ψ₂ (complex triple product)
 *    |P|² = (Re P)² + (Im P)²
 *
 *  This model has U(1)² global symmetry:
 *    ψ₀ → e^{iα}ψ₀, ψ₁ → e^{iβ}ψ₁, ψ₂ → e^{-i(α+β)}ψ₂
 *  preserving P. Two conserved charges → two "EM-like" currents.
 *
 *  Questions:
 *  1. Does the braid survive the upgrade?
 *  2. Are the conserved charges localized on the braid?
 *  3. Does phase rotation create a long-range current?
 *
 *  Implementation: 6 real fields (φ₀,χ₀,φ₁,χ₁,φ₂,χ₂)
 *  P = (φ₀+iχ₀)(φ₁+iχ₁)(φ₂+iχ₂)
 *
 *  Build: gcc -O3 -fopenmp -o t10g src/t10g.c -lm
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>

#define NREAL 6  /* 3 complex = 6 real */
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

/* Complex triple product: P = ψ₀ψ₁ψ₂ where ψ_a = (φ_{2a} + i φ_{2a+1})
 * Returns (Re P, Im P) */
static void complex_product(double r0, double i0, double r1, double i1,
                            double r2, double i2, double *reP, double *imP) {
    /* (r0+i*i0)*(r1+i*i1) = (r0*r1-i0*i1) + i*(r0*i1+i0*r1) */
    double t_re = r0*r1 - i0*i1;
    double t_im = r0*i1 + i0*r1;
    /* * (r2+i*i2) */
    *reP = t_re*r2 - t_im*i2;
    *imP = t_re*i2 + t_im*r2;
}

/* Compute forces for complex field model */
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

        /* Complex fields at this point */
        double r0 = g->phi[0][idx], i0 = g->phi[1][idx]; /* ψ₀ */
        double r1 = g->phi[2][idx], i1 = g->phi[3][idx]; /* ψ₁ */
        double r2 = g->phi[4][idx], i2 = g->phi[5][idx]; /* ψ₂ */

        /* Complex triple product P = ψ₀ψ₁ψ₂ */
        double reP, imP;
        complex_product(r0, i0, r1, i1, r2, i2, &reP, &imP);
        double P2 = reP*reP + imP*imP;

        /* V = (μ/2)|P|²/(1+κ|P|²) */
        double D = 1.0 + kappa * P2;
        double D2 = D * D;

        /* dV/dP_R = μ P_R / D² (where P_R = Re P) */
        /* dV/dP_I = μ P_I / D² */
        double dVdPr = mu * reP / D2;
        double dVdPi = mu * imP / D2;

        /* dP/dψ_0 = ψ₁ψ₂, etc.
         * d(Re P)/d(r0) = Re(ψ₁ψ₂), d(Re P)/d(i0) = -Im(ψ₁ψ₂)
         * d(Im P)/d(r0) = Im(ψ₁ψ₂), d(Im P)/d(i0) = Re(ψ₁ψ₂)
         */

        /* ψ₁ψ₂ */
        double dP0_re = r1*r2 - i1*i2;
        double dP0_im = r1*i2 + i1*r2;
        /* ψ₀ψ₂ */
        double dP1_re = r0*r2 - i0*i2;
        double dP1_im = r0*i2 + i0*r2;
        /* ψ₀ψ₁ */
        double dP2_re = r0*r1 - i0*i1;
        double dP2_im = r0*i1 + i0*r1;

        /* Force from potential: -dV/dφ_a */
        /* For ψ₀ = (r0, i0):
         * -dV/dr0 = -(dVdPr * d(ReP)/dr0 + dVdPi * d(ImP)/dr0)
         *         = -(dVdPr * dP0_re + dVdPi * dP0_im)
         * -dV/di0 = -(dVdPr * (-dP0_im) + dVdPi * dP0_re)
         *         = -(-dVdPr * dP0_im + dVdPi * dP0_re)
         */
        double fV[NREAL];
        fV[0] = -(dVdPr * dP0_re + dVdPi * dP0_im); /* -dV/dr0 */
        fV[1] = -(dVdPi * dP0_re - dVdPr * dP0_im); /* -dV/di0 */
        fV[2] = -(dVdPr * dP1_re + dVdPi * dP1_im);
        fV[3] = -(dVdPi * dP1_re - dVdPr * dP1_im);
        fV[4] = -(dVdPr * dP2_re + dVdPi * dP2_im);
        fV[5] = -(dVdPi * dP2_re - dVdPr * dP2_im);

        for (int a = 0; a < NREAL; a++) {
            double lap = (g->phi[a][idx+NN] + g->phi[a][idx-NN]
                        + g->phi[a][idx+N]  + g->phi[a][idx-N]
                        + g->phi[a][idx_kp] + g->phi[a][idx_km]
                        - 6.0 * g->phi[a][idx]) * idx2;
            g->acc[a][idx] = lap - mass2 * g->phi[a][idx] + fV[a];
        }
    }
}

static void cverlet_step(CGrid *g, double mu, double kappa, double mass2) {
    int N3 = g->N * g->N * g->N;
    double hdt = 0.5 * g->dt;
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
    double inv_dr = 1.0/(r_end - r_start + 1e-30);
    for (int i = 0; i < N; i++) {
        double x = -L + i*dx;
        for (int j = 0; j < N; j++) {
            double y = -L + j*dx;
            double rp = sqrt(x*x + y*y);
            if (rp <= r_start) continue;
            double f = (rp - r_start) * inv_dr;
            if (f > 1.0) f = 1.0;
            double damp = 1.0 - 0.98*f*f;
            for (int kk = 0; kk < N; kk++) {
                int idx = i*NN + j*N + kk;
                for (int a = 0; a < NREAL; a++) {
                    g->phi[a][idx] *= damp;
                    g->vel[a][idx] *= damp;
                }
            }
        }
    }
}

/* Measure conserved U(1) charge: Q_α = Σ Im(ψ*_a ∂_t ψ_a)
 * For charge associated with ψ₀ → e^{iα}ψ₀:
 * Q₀ = Im(ψ₀* ∂_t ψ₀) = r0*v_i0 - i0*v_r0 (at each point) */
static double measure_charge(CGrid *g, int field_idx) {
    int N = g->N, NN = N*N, N3 = N*N*N;
    double dx = g->dx, dV = dx*dx*dx;
    double Q = 0;
    int r_idx = 2*field_idx, i_idx = 2*field_idx+1;
    for (int idx = 0; idx < N3; idx++) {
        Q += (g->phi[r_idx][idx]*g->vel[i_idx][idx]
            - g->phi[i_idx][idx]*g->vel[r_idx][idx]) * dV;
    }
    return Q;
}

/* Measure charge density in core (r < R_core) */
static double measure_charge_core(CGrid *g, int field_idx, double R_core) {
    int N = g->N, NN = N*N;
    double dx = g->dx, L = g->L, dV = dx*dx*dx;
    double Q = 0;
    double Rc2 = R_core * R_core;
    int r_idx = 2*field_idx, i_idx = 2*field_idx+1;
    for (int i = 1; i < N-1; i++) {
        double x = -L + i*dx;
        for (int j = 1; j < N-1; j++) {
            double y = -L + j*dx;
            if (x*x + y*y > Rc2) continue;
            for (int k = 0; k < N; k++) {
                int idx = i*NN + j*N + k;
                Q += (g->phi[r_idx][idx]*g->vel[i_idx][idx]
                    - g->phi[i_idx][idx]*g->vel[r_idx][idx]) * dV;
            }
        }
    }
    return Q;
}

/* fc for complex field */
static double measure_fc(CGrid *g) {
    int N = g->N, NN = N*N;
    double dx = g->dx, L = g->L;
    double R_core = 8.0;
    if (R_core > L/3.0) R_core = L/3.0;
    double Rc2 = R_core * R_core;
    double core = 0, total = 0;
    for (int i = 1; i < N-1; i++) {
        double x = -L + i*dx;
        for (int j = 1; j < N-1; j++) {
            double y = -L + j*dx;
            double rp2 = x*x + y*y;
            for (int k = 0; k < N; k++) {
                int idx = i*NN + j*N + k;
                double rho = 0;
                for (int a = 0; a < NREAL; a++)
                    rho += g->phi[a][idx] * g->phi[a][idx];
                total += rho;
                if (rp2 < Rc2) core += rho;
            }
        }
    }
    return core / (total + 1e-30);
}

int main(int argc, char **argv) {
    int N = 96;
    double L = 20.0;
    double T_equil = 300.0;

    double mu = -41.3, kappa = 50.0, mass2 = 2.25;

    printf("=== T10G: Complex Field Upgrade ===\n\n");
    printf("Grid: N=%d, L=%.1f\n", N, L);
    printf("6 real fields (3 complex scalars)\n");
    printf("V = (mu/2)|P|^2/(1+kappa|P|^2), P = psi0*psi1*psi2\n");
    printf("mu=%.1f, kappa=%.1f, m^2=%.2f\n\n", mu, kappa, mass2);

    CGrid *g = cgrid_alloc(N, L);
    double dx = g->dx, dt = g->dt;
    int NN = N*N;

    /* Initialize: helical braid in real parts, small perturbation in imaginary parts
     * This is like the real braid with a small twist in phase space */
    printf("Initializing complex braid...\n");
    {
        double A = 0.8, R_tube = 3.0, ellip = 0.33;
        double k = PI / L;
        double omega = sqrt(k*k + mass2);
        double sx = 1.0 + ellip, sy = 1.0 - ellip;
        double inv_2R2 = 1.0 / (2.0*R_tube*R_tube);

        double delta[3] = {0.0, 3.53*0.85, 4.92*0.85}; /* bimodal phases */
        double phase_twist = 0.3; /* small imaginary component */

        for (int i = 0; i < N; i++) {
            double x = -L + i*dx;
            for (int j = 0; j < N; j++) {
                double y = -L + j*dx;
                for (int kk = 0; kk < N; kk++) {
                    double z = -L + kk*dx;
                    int idx = i*NN + j*N + kk;
                    for (int c = 0; c < 3; c++) {
                        double xr = x/(sx), yr = y/(sy);
                        double r2e = xr*xr + yr*yr;
                        double env = A * exp(-r2e * inv_2R2);
                        double ph = k * z + delta[c];
                        /* Real part */
                        g->phi[2*c][idx] = env * cos(ph);
                        g->vel[2*c][idx] = omega * env * sin(ph);
                        /* Imaginary part: small twist */
                        g->phi[2*c+1][idx] = phase_twist * env * sin(ph);
                        g->vel[2*c+1][idx] = phase_twist * omega * env * cos(ph);
                    }
                }
            }
        }
    }

    ccompute_forces(g, mu, kappa, mass2);

    /* Equilibrate */
    printf("Equilibrating (T=%.0f)...\n", T_equil);
    fflush(stdout);
    int equil_steps = (int)(T_equil / dt);

    /* Monitor charges during equilibration */
    for (int s = 0; s < equil_steps; s++) {
        cverlet_step(g, mu, kappa, mass2);
        capply_damping_xy(g);

        if (s % (equil_steps/5) == 0) {
            double Q0 = measure_charge(g, 0);
            double Q1 = measure_charge(g, 1);
            double Q2 = measure_charge(g, 2);
            double fc = measure_fc(g);
            printf("  t=%5.0f: fc=%.4f, Q0=%.2f, Q1=%.2f, Q2=%.2f\n",
                   s*dt, fc, Q0, Q1, Q2);
        }
    }

    printf("\nAfter equilibration:\n");
    double fc = measure_fc(g);
    double Q0 = measure_charge(g, 0);
    double Q1 = measure_charge(g, 1);
    double Q2 = measure_charge(g, 2);
    double Q0c = measure_charge_core(g, 0, 8.0);
    double Q1c = measure_charge_core(g, 1, 8.0);
    double Q2c = measure_charge_core(g, 2, 8.0);

    printf("  fc = %.4f\n", fc);
    printf("  Total charges: Q0=%.4f, Q1=%.4f, Q2=%.4f\n", Q0, Q1, Q2);
    printf("  Core charges:  Q0=%.4f, Q1=%.4f, Q2=%.4f\n", Q0c, Q1c, Q2c);
    printf("  Conservation: Q0+Q1+Q2 = %.6f (should be ~0)\n", Q0+Q1+Q2);

    /* Track charge conservation during free evolution */
    printf("\nFree evolution (T=100)...\n");
    FILE *fout = fopen("data/t10g_charges.tsv", "w");
    fprintf(fout, "t\tfc\tQ0\tQ1\tQ2\tQ_total\tQ0_core\n");

    double T_free = 100.0;
    int free_steps = (int)(T_free / dt);
    int meas_every = free_steps / 50;
    if (meas_every < 1) meas_every = 1;

    for (int s = 0; s < free_steps; s++) {
        cverlet_step(g, mu, kappa, mass2);
        capply_damping_xy(g);

        if (s % meas_every == 0) {
            double t = s * dt;
            double fc2 = measure_fc(g);
            double q0 = measure_charge(g, 0);
            double q1 = measure_charge(g, 1);
            double q2 = measure_charge(g, 2);
            double q0c = measure_charge_core(g, 0, 8.0);
            fprintf(fout, "%.3f\t%.4f\t%.4f\t%.4f\t%.4f\t%.6f\t%.4f\n",
                    t, fc2, q0, q1, q2, q0+q1+q2, q0c);
        }
    }

    fclose(fout);

    /* Final state */
    fc = measure_fc(g);
    Q0 = measure_charge(g, 0);
    Q1 = measure_charge(g, 1);
    Q2 = measure_charge(g, 2);

    printf("\nFinal state:\n");
    printf("  fc = %.4f\n", fc);
    printf("  Charges: Q0=%.4f, Q1=%.4f, Q2=%.4f\n", Q0, Q1, Q2);
    printf("  Sum = %.6f\n", Q0+Q1+Q2);

    if (fc > 0.5)
        printf("\n  BRAID SURVIVES complex upgrade!\n");
    else
        printf("\n  Braid does NOT survive complex upgrade.\n");

    if (fabs(Q0) > 0.1 || fabs(Q1) > 0.1)
        printf("  Non-zero charges detected → CHARGE-CARRYING SOLITON!\n");
    else
        printf("  Charges near zero → no net charge on soliton.\n");

    cgrid_free(g);

    printf("\n=== T10G Complete ===\n");
    printf("Data in data/t10g_charges.tsv\n");
    return 0;
}
