/*  t10c.c — Dispersion Relation Around a Condensate
 *
 *  Can mass EMERGE from the field configuration?
 *
 *  Linearize the EOM around three backgrounds:
 *  (a) Vacuum: φ_a = 0 → ω² = k² + m² (trivial)
 *  (b) Uniform condensate: φ_a = φ₀ for all a → modified dispersion
 *  (c) 1D modulated state: φ_a = φ₀ cos(Kz) → band structure
 *
 *  For (b): the Hessian around φ_a = φ₀ gives M²_ab = m²δ_ab + ∂²V/∂φ_a∂φ_b
 *  where ∂²V/∂φ_a∂φ_b evaluated at uniform condensate.
 *
 *  If any eigenvalue of M² is NEGATIVE → tachyonic instability → the uniform
 *  condensate is UNSTABLE → spontaneous structure formation! This would mean
 *  braids can form from a homogeneous field.
 *
 *  For (c): solve the Hill equation (Floquet theory) to get band structure.
 *  Gaps in the band structure = mass generation.
 *
 *  Build: gcc -O3 -fopenmp -o t10c src/t10c.c -lm
 */

#include "../../src/braid_core.h"
#include <float.h>

/* Compute Hessian at uniform condensate φ_a = v for all a */
static void hessian_uniform(double v, double mu, double kappa, double mass2,
                            double H[3][3]) {
    double P = v * v * v;
    double P2 = P * P;
    double D = 1.0 + kappa * P2;
    double D2 = D * D;
    double D3 = D2 * D;

    /* dP/dφ_a = v² (all same) */
    double dP = v * v;

    /* d²P/dφ_a dφ_b = v for a≠b, 0 for a=b */

    double f1 = mu * (1.0 - 3.0*kappa*P2) / D3;
    double f2 = mu * P / D2;

    for (int a = 0; a < 3; a++) {
        for (int b = 0; b < 3; b++) {
            H[a][b] = f1 * dP * dP; /* dP_a * dP_b (all same for uniform) */
            if (a != b) H[a][b] += f2 * v;
            if (a == b) H[a][b] += mass2;
        }
    }
}

/* Eigenvalues of 3x3 symmetric matrix */
static void eigen3(double A[3][3], double eig[3]) {
    double p1 = A[0][1]*A[0][1] + A[0][2]*A[0][2] + A[1][2]*A[1][2];
    double q = (A[0][0] + A[1][1] + A[2][2]) / 3.0;

    if (p1 < 1e-30) {
        eig[0] = A[0][0]; eig[1] = A[1][1]; eig[2] = A[2][2];
        if (eig[0] > eig[1]) { double t=eig[0]; eig[0]=eig[1]; eig[1]=t; }
        if (eig[1] > eig[2]) { double t=eig[1]; eig[1]=eig[2]; eig[2]=t; }
        if (eig[0] > eig[1]) { double t=eig[0]; eig[0]=eig[1]; eig[1]=t; }
        return;
    }

    double p2 = (A[0][0]-q)*(A[0][0]-q) + (A[1][1]-q)*(A[1][1]-q)
              + (A[2][2]-q)*(A[2][2]-q) + 2*p1;
    double p = sqrt(p2 / 6.0);
    double invp = 1.0 / p;

    double B[3][3];
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            B[i][j] = (A[i][j] - (i==j ? q : 0)) * invp;

    double detB = B[0][0]*(B[1][1]*B[2][2]-B[1][2]*B[2][1])
                - B[0][1]*(B[1][0]*B[2][2]-B[1][2]*B[2][0])
                + B[0][2]*(B[1][0]*B[2][1]-B[1][1]*B[2][0]);

    double r = detB / 2.0;
    if (r <= -1) r = -1; if (r >= 1) r = 1;
    double phi = acos(r) / 3.0;

    eig[0] = q + 2*p*cos(phi + 2*PI/3);
    eig[1] = q + 2*p*cos(phi);
    eig[2] = 3*q - eig[0] - eig[1];

    if (eig[0] > eig[1]) { double t=eig[0]; eig[0]=eig[1]; eig[1]=t; }
    if (eig[1] > eig[2]) { double t=eig[1]; eig[1]=eig[2]; eig[2]=t; }
    if (eig[0] > eig[1]) { double t=eig[0]; eig[0]=eig[1]; eig[1]=t; }
}

/* ================================================================
   Numerical dispersion: initialize plane wave on uniform condensate,
   measure ω from oscillation period
   ================================================================ */

/* 1D grid for fast dispersion measurement */
typedef struct {
    double phi[NFIELDS][4096]; /* 1D fields */
    double vel[NFIELDS][4096];
    double acc[NFIELDS][4096];
    int N;
    double L, dx, dt;
} Grid1D;

static void compute_forces_1d(Grid1D *g1, double mu, double kappa, double mass2,
                               double bg[3]) {
    int N = g1->N;
    double idx2 = 1.0 / (g1->dx * g1->dx);

    for (int i = 0; i < N; i++) {
        int ip = (i+1) % N, im = (i-1+N) % N;
        /* Total field = bg + perturbation */
        double f[3], df[3];
        for (int a = 0; a < 3; a++) {
            f[a] = bg[a] + g1->phi[a][i];
            df[a] = g1->phi[a][i];
        }

        double P = f[0]*f[1]*f[2];
        double denom = 1.0 + kappa*P*P;
        double mu_P_d2 = mu * P / (denom*denom);

        for (int a = 0; a < 3; a++) {
            double lap = (g1->phi[a][ip] + g1->phi[a][im] - 2.0*g1->phi[a][i]) * idx2;
            double dPda = (a==0) ? f[1]*f[2] : (a==1) ? f[0]*f[2] : f[0]*f[1];
            /* Subtract the background force (which should be zero for uniform bg) */
            double P_bg = bg[0]*bg[1]*bg[2];
            double denom_bg = 1.0 + kappa*P_bg*P_bg;
            double dPda_bg = (a==0)?bg[1]*bg[2]:(a==1)?bg[0]*bg[2]:bg[0]*bg[1];
            double f_bg = mu*P_bg/(denom_bg*denom_bg)*dPda_bg;

            g1->acc[a][i] = lap - mass2*g1->phi[a][i]
                          - (mu_P_d2*dPda - f_bg - mass2*bg[a]);
        }
    }
}

int main(int argc, char **argv) {
    double mu, kappa, mass2;

    bimodal_init_params();
    mu = BIMODAL[12];
    kappa = BIMODAL[13];
    mass2 = BIMODAL[14] * BIMODAL[14]; /* = 2.25 */

    printf("=== T10C: Dispersion Relation ===\n\n");
    printf("Parameters: mu=%.1f, kappa=%.1f, m²=%.2f\n\n", mu, kappa, mass2);

    FILE *fout = fopen("data/t10c_dispersion.tsv", "w");
    fprintf(fout, "condensate\tmode\tk\tomega_sq\teff_mass_sq\tstable\n");

    /* ============================================================
       Part 1: Analytical Hessian for uniform condensate
       ============================================================ */
    printf("=== Part 1: Hessian eigenvalues vs condensate amplitude ===\n\n");
    printf("  v\tM²_1\tM²_2\tM²_3\tstable?\n");

    double v_values[] = {0, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0, 5.0};
    int nv = 11;

    double v_tachyonic = -1; /* first v where any eigenvalue goes negative */

    for (int iv = 0; iv < nv; iv++) {
        double v = v_values[iv];
        double H[3][3];
        hessian_uniform(v, mu, kappa, mass2, H);
        double eig[3];
        eigen3(H, eig);

        int stable = (eig[0] >= 0) ? 1 : 0;
        printf("  %.1f\t%.4f\t%.4f\t%.4f\t%s\n", v,
               eig[0], eig[1], eig[2], stable ? "yes" : "TACHYONIC");

        for (int e = 0; e < 3; e++)
            fprintf(fout, "%.4f\tanalytic_%d\t0\t%.6f\t%.6f\t%d\n",
                    v, e, eig[e], eig[e]-mass2, stable);

        if (!stable && v_tachyonic < 0)
            v_tachyonic = v;
    }

    if (v_tachyonic > 0) {
        printf("\n  TACHYONIC INSTABILITY at v = %.2f!\n", v_tachyonic);
        printf("  → Uniform condensate is UNSTABLE above this amplitude\n");
        printf("  → Spontaneous structure formation is POSSIBLE\n");
    } else {
        printf("\n  No tachyonic instability found. Condensate is stable at all tested v.\n");
    }

    /* ============================================================
       Part 2: Find exact tachyonic threshold by bisection
       ============================================================ */
    printf("\n=== Part 2: Tachyonic threshold ===\n");
    {
        double vlo = 0, vhi = 10;
        for (int iter = 0; iter < 60; iter++) {
            double vmid = 0.5*(vlo+vhi);
            double H[3][3];
            hessian_uniform(vmid, mu, kappa, mass2, H);
            double eig[3];
            eigen3(H, eig);
            if (eig[0] < 0) vhi = vmid; else vlo = vmid;
        }
        printf("  Tachyonic threshold: v_c = %.6f\n", 0.5*(vlo+vhi));
        printf("  At v_c: ");
        double H[3][3];
        hessian_uniform(0.5*(vlo+vhi), mu, kappa, mass2, H);
        double eig[3];
        eigen3(H, eig);
        printf("M² = (%.6f, %.6f, %.6f)\n", eig[0], eig[1], eig[2]);

        /* Also find where the MAXIMUM eigenvalue equals zero (complete instability) */
        vlo = 0; vhi = 100;
        for (int iter = 0; iter < 60; iter++) {
            double vmid = 0.5*(vlo+vhi);
            double HH[3][3];
            hessian_uniform(vmid, mu, kappa, mass2, HH);
            double ee[3];
            eigen3(HH, ee);
            if (ee[2] < 0) vhi = vmid; else vlo = vmid;
        }
        printf("  All-tachyonic threshold: v_all = %.6f\n", 0.5*(vlo+vhi));
    }

    /* ============================================================
       Part 3: Which MODE goes tachyonic first?
       At uniform condensate, the symmetry tells us: the 3x3 Hessian has
       one eigenvalue = m² + 2f1·dP² + 2f2·v (the symmetric mode)
       and two degenerate eigenvalues (the antisymmetric modes).
       ============================================================ */
    printf("\n=== Part 3: Mode identification ===\n");
    {
        /* At v = 0.5 (near threshold) */
        double v = 0.5;
        double H[3][3];
        hessian_uniform(v, mu, kappa, mass2, H);
        printf("  Hessian at v=%.1f:\n", v);
        for (int a = 0; a < 3; a++)
            printf("    [%.6f, %.6f, %.6f]\n", H[a][0], H[a][1], H[a][2]);

        /* The symmetric mode: (1,1,1)/sqrt(3) */
        /* The two antisymmetric modes: (1,-1,0)/sqrt(2) and (1,1,-2)/sqrt(6) */
        double sym_M2 = 0;
        for (int a = 0; a < 3; a++)
            for (int b = 0; b < 3; b++)
                sym_M2 += H[a][b] / 3.0;

        double asym_M2 = H[0][0] - H[0][1]; /* (1,-1,0) mode */

        printf("  Symmetric mode M² (all fields in phase): %.6f\n", sym_M2);
        printf("  Antisymmetric mode M² (fields out of phase): %.6f\n", asym_M2);
        if (sym_M2 < asym_M2)
            printf("  → SYMMETRIC mode goes tachyonic first (braid-like instability)\n");
        else
            printf("  → ANTISYMMETRIC mode goes tachyonic first (texture-like instability)\n");
    }

    /* ============================================================
       Part 4: Numerical verification (1D simulation)
       Initialize perturbation on condensate, check if it grows
       ============================================================ */
    printf("\n=== Part 4: Numerical growth rate test ===\n");
    {
        int N1d = 512;
        double L1d = 20.0;
        double dx1d = 2.0*L1d / N1d;
        double dt1d = 0.1 * dx1d;
        double T_test = 50.0;
        int steps = (int)(T_test / dt1d);

        /* Test at v=0.5 (near threshold) and v=1.0 (above if tachyonic) */
        double test_v[] = {0.3, 0.5, 1.0, 2.0};
        int ntv = 4;
        double test_k = 0.5; /* low k to see tachyonic growth */
        double pert_amp = 1e-6;

        printf("  v\tk\tgrowth_rate\tresult\n");

        for (int iv = 0; iv < ntv; iv++) {
            double v = test_v[iv];
            double bg[3] = {v, v, v};

            Grid1D g1;
            g1.N = N1d; g1.L = L1d; g1.dx = dx1d; g1.dt = dt1d;

            /* Initialize small perturbation: symmetric mode */
            for (int i = 0; i < N1d; i++) {
                double z = -L1d + i * dx1d;
                for (int a = 0; a < 3; a++) {
                    g1.phi[a][i] = pert_amp * cos(test_k * z);
                    g1.vel[a][i] = 0;
                }
            }

            /* Compute initial forces */
            compute_forces_1d(&g1, mu, kappa, mass2, bg);

            double amp0 = pert_amp;
            double amp_final = 0;

            for (int s = 0; s < steps; s++) {
                /* Verlet */
                for (int a = 0; a < 3; a++)
                    for (int i = 0; i < N1d; i++)
                        g1.vel[a][i] += 0.5 * dt1d * g1.acc[a][i];
                for (int a = 0; a < 3; a++)
                    for (int i = 0; i < N1d; i++)
                        g1.phi[a][i] += dt1d * g1.vel[a][i];
                compute_forces_1d(&g1, mu, kappa, mass2, bg);
                for (int a = 0; a < 3; a++)
                    for (int i = 0; i < N1d; i++)
                        g1.vel[a][i] += 0.5 * dt1d * g1.acc[a][i];

                /* Track amplitude */
                double max_a = 0;
                for (int i = 0; i < N1d; i++) {
                    double a2 = 0;
                    for (int a = 0; a < 3; a++)
                        a2 += g1.phi[a][i] * g1.phi[a][i];
                    if (a2 > max_a) max_a = a2;
                }
                amp_final = sqrt(max_a);

                if (amp_final > 1.0) { /* clearly growing */
                    double gamma = log(amp_final / amp0) / (s * dt1d);
                    printf("  %.1f\t%.1f\tGROWTH γ=%.4f at t=%.1f\n",
                           v, test_k, gamma, s*dt1d);
                    fprintf(fout, "%.4f\tgrowth\t%.4f\t%.6f\t%.6f\t0\n",
                            v, test_k, -gamma*gamma, -gamma*gamma);
                    break;
                }
            }

            if (amp_final <= 1.0) {
                double growth = log(amp_final / (amp0 * sqrt(3.0))) / T_test;
                printf("  %.1f\t%.1f\t%.6f\t%s\n", v, test_k, growth,
                       growth > 0.01 ? "growing" : "stable/oscillating");
            }
        }
    }

    /* ============================================================
       Part 5: Effective mass vs condensate amplitude (full scan)
       ============================================================ */
    printf("\n=== Part 5: Effective mass² vs condensate v (fine scan) ===\n");
    {
        FILE *feff = fopen("data/t10c_effective_mass.tsv", "w");
        fprintf(feff, "v\tM2_min\tM2_mid\tM2_max\tP\n");

        for (int iv = 0; iv <= 100; iv++) {
            double v = iv * 0.05;
            double H[3][3];
            hessian_uniform(v, mu, kappa, mass2, H);
            double eig[3];
            eigen3(H, eig);
            double P = v*v*v;
            fprintf(feff, "%.4f\t%.6f\t%.6f\t%.6f\t%.6f\n",
                    v, eig[0], eig[1], eig[2], P);
        }
        fclose(feff);
        printf("  Written to data/t10c_effective_mass.tsv\n");
    }

    fclose(fout);

    printf("\n=== SUMMARY ===\n");
    printf("Key questions:\n");
    printf("  1. Is the uniform condensate UNSTABLE? → does a tachyonic mode exist?\n");
    printf("  2. If yes, at what amplitude? → that determines the condensate density\n");
    printf("  3. Which mode goes tachyonic? → determines what structures form\n");
    printf("  4. Does this explain the braid mass? → only if mass² shifts match\n");

    return 0;
}
