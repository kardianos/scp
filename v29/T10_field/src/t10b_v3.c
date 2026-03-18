/*  t10b_v3.c — Effective Metric from Linearized EOM (direct Hessian approach)
 *
 *  Instead of tracking dispersive wave packets, directly compute the
 *  linearized equation of motion at each point on the frozen background:
 *
 *    d²δφ_a/dt² = ∇²δφ_a - M²_ab(x) δφ_b
 *
 *  where M²_ab(x) = m²δ_ab + ∂²V/∂φ_a∂φ_b |_{bg}
 *
 *  The eigenvalues of M²_ab(x) give the local effective mass squared
 *  for each mode. The group velocity for a mode with mass² = M² and
 *  wave number k is:
 *
 *    v_g = k / sqrt(k² + M²)
 *
 *  Where M² is LARGER, waves are SLOWER → gravitational lensing.
 *  Where M² is SMALLER (or negative), waves are FASTER → anti-lensing.
 *
 *  This gives us the effective refractive index:
 *    n(x) = v_g(vacuum) / v_g(x) = sqrt(k²+M²(x)) / sqrt(k²+m²)
 *
 *  For GRAVITY we need n > 1 near the braid (slower = attractive).
 *
 *  ALSO: check if the Hessian is anisotropic (different eigenvalues
 *  for different modes). If so, the effective metric is NOT conformally
 *  flat — it has genuine tensor structure.
 *
 *  Build: gcc -O3 -fopenmp -o t10b_v3 src/t10b_v3.c -lm
 */

#include "../../src/braid_core.h"
#include <float.h>

/* ================================================================
   Compute the 3x3 Hessian matrix ∂²V/∂φ_a∂φ_b at a point
   V = (μ/2)P²/(1+κP²), P = φ₀φ₁φ₂
   ================================================================ */

static void compute_hessian(double b0, double b1, double b2,
                            double mu, double kappa,
                            double H[3][3]) {
    double P = b0 * b1 * b2;
    double P2 = P * P;
    double D = 1.0 + kappa * P2;
    double D2 = D * D;
    double D3 = D2 * D;

    /* dP/dφ_a */
    double dP[3] = {b1*b2, b0*b2, b0*b1};

    /* d²P/dφ_a dφ_b:
     * (0,0)=0, (0,1)=b2, (0,2)=b1
     * (1,0)=b2, (1,1)=0, (1,2)=b0
     * (2,0)=b1, (2,1)=b0, (2,2)=0
     */
    double d2P[3][3] = {{0,b2,b1},{b2,0,b0},{b1,b0,0}};

    /* ∂²V/∂φ_a∂φ_b = μ * [ dP_a dP_b (1-3κP²)/D³ + P d²P_ab/D² ] */
    double f1 = mu * (1.0 - 3.0*kappa*P2) / D3;
    double f2 = mu * P / D2;

    for (int a = 0; a < 3; a++)
        for (int b = 0; b < 3; b++)
            H[a][b] = f1 * dP[a] * dP[b] + f2 * d2P[a][b];
}

/* Eigenvalues of 3x3 symmetric matrix using analytical formula */
static void eigen3(double A[3][3], double eig[3]) {
    double p1 = A[0][1]*A[0][1] + A[0][2]*A[0][2] + A[1][2]*A[1][2];
    double q = (A[0][0] + A[1][1] + A[2][2]) / 3.0;

    if (p1 < 1e-30) {
        /* diagonal */
        eig[0] = A[0][0]; eig[1] = A[1][1]; eig[2] = A[2][2];
        /* sort */
        if (eig[0] > eig[1]) { double t=eig[0]; eig[0]=eig[1]; eig[1]=t; }
        if (eig[1] > eig[2]) { double t=eig[1]; eig[1]=eig[2]; eig[2]=t; }
        if (eig[0] > eig[1]) { double t=eig[0]; eig[0]=eig[1]; eig[1]=t; }
        return;
    }

    double p2 = (A[0][0]-q)*(A[0][0]-q) + (A[1][1]-q)*(A[1][1]-q)
              + (A[2][2]-q)*(A[2][2]-q) + 2*p1;
    double p = sqrt(p2 / 6.0);
    double invp = 1.0 / p;

    /* B = (1/p)(A - qI) */
    double B[3][3];
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            B[i][j] = (A[i][j] - (i==j ? q : 0)) * invp;

    /* det(B) */
    double detB = B[0][0]*(B[1][1]*B[2][2]-B[1][2]*B[2][1])
                - B[0][1]*(B[1][0]*B[2][2]-B[1][2]*B[2][0])
                + B[0][2]*(B[1][0]*B[2][1]-B[1][1]*B[2][0]);

    double r = detB / 2.0;
    if (r <= -1) r = -1; if (r >= 1) r = 1;
    double phi = acos(r) / 3.0;

    eig[0] = q + 2*p*cos(phi + 2*PI/3);
    eig[1] = q + 2*p*cos(phi);
    eig[2] = 3*q - eig[0] - eig[1];

    /* sort ascending */
    if (eig[0] > eig[1]) { double t=eig[0]; eig[0]=eig[1]; eig[1]=t; }
    if (eig[1] > eig[2]) { double t=eig[1]; eig[1]=eig[2]; eig[2]=t; }
    if (eig[0] > eig[1]) { double t=eig[0]; eig[0]=eig[1]; eig[1]=t; }
}

/* ================================================================
   Main
   ================================================================ */

int main(int argc, char **argv) {
    int N = 96;
    double L = 20.0;
    double T_equil = 200.0;
    double mass2 = 2.25;
    double k_test = 4.0; /* reference wave number for n_eff calculation */

    printf("=== T10B v3: Effective Metric from Hessian ===\n\n");
    printf("Grid: N=%d, L=%.1f\n", N, L);
    printf("Reference k=%.1f for refractive index\n", k_test);
    printf("Vacuum mass²=%.2f, omega_vac=%.4f, vg_vac=%.4f\n\n",
           mass2, sqrt(k_test*k_test+mass2), k_test/sqrt(k_test*k_test+mass2));

    /* Equilibrate braid */
    printf("Equilibrating braid...\n");
    fflush(stdout);
    bimodal_init_params();
    Grid *g = grid_alloc(N, L);
    init_braid(g, BIMODAL, -1);
    compute_forces(g, BIMODAL, 2.25);

    int equil_steps = (int)(T_equil / g->dt);
    for (int s = 0; s < equil_steps; s++) {
        verlet_full_step(g, BIMODAL, 2.25);
        apply_damping_xy(g);
    }
    Result res;
    snprintf(res.label, sizeof(res.label), "equil");
    compute_diagnostics(g, BIMODAL, &res);
    printf("  fc=%.4f, energy=%.1f\n\n", res.fc, res.energy);

    double mu = BIMODAL[12], kappa = BIMODAL[13];
    double dx = g->dx;
    int NN = N * N;

    /* Open output files */
    FILE *frad = fopen("data/t10b_v3_radial.tsv", "w");
    FILE *f2d  = fopen("data/t10b_v3_xy_slice.tsv", "w");
    FILE *fzax = fopen("data/t10b_v3_zaxis.tsv", "w");

    fprintf(frad, "r\tM2_min\tM2_mid\tM2_max\tn_eff_min\tn_eff_max\tanisotropy\tbg_rho\n");
    fprintf(f2d, "x\ty\tM2_min\tM2_max\tn_eff_avg\tbg_P\n");
    fprintf(fzax, "z\tM2_min\tM2_mid\tM2_max\tbg_rho\tbg_P\n");

    /* ============================================================
       1. Radial profile: average Hessian eigenvalues vs r (xy-plane)
       ============================================================ */
    printf("=== Radial profile of effective mass² ===\n");
    printf("  r\tM²_min\tM²_mid\tM²_max\tn_eff_min\tn_eff_max\tanisotropy\n");

    int nr = 40;
    double dr = L / nr;

    for (int ir = 0; ir < nr; ir++) {
        double r = (ir + 0.5) * dr;
        /* Average over azimuth at z=0 (mid-plane) */
        int k_mid = N / 2;
        double sum_eig[3] = {0,0,0};
        double sum_rho = 0, sum_P = 0;
        int count = 0;

        for (int i = 1; i < N-1; i++) {
            double x = -L + i * dx;
            for (int j = 1; j < N-1; j++) {
                double y = -L + j * dx;
                double rp = sqrt(x*x + y*y);
                if (fabs(rp - r) > dr * 0.6) continue;
                int idx = i * NN + j * N + k_mid;
                double b0 = g->phi[0][idx], b1 = g->phi[1][idx], b2 = g->phi[2][idx];

                double H[3][3];
                compute_hessian(b0, b1, b2, mu, kappa, H);
                /* Add the bare mass */
                for (int a = 0; a < 3; a++)
                    H[a][a] += mass2;

                double eig[3];
                eigen3(H, eig);

                sum_eig[0] += eig[0];
                sum_eig[1] += eig[1];
                sum_eig[2] += eig[2];
                sum_rho += b0*b0 + b1*b1 + b2*b2;
                sum_P += fabs(b0*b1*b2);
                count++;
            }
        }

        if (count > 0) {
            for (int e = 0; e < 3; e++) sum_eig[e] /= count;
            sum_rho /= count;
            sum_P /= count;

            /* Refractive index: n = sqrt(k²+M²) / sqrt(k²+m²_vac)
             * where m²_vac = 2.25 */
            double omega_vac = sqrt(k_test*k_test + mass2);
            double n_min = sqrt(k_test*k_test + sum_eig[0]) / omega_vac;
            double n_max = sqrt(k_test*k_test + sum_eig[2]) / omega_vac;
            double aniso = (sum_eig[2] - sum_eig[0]) / (fabs(sum_eig[1]) + 1e-30);

            printf("  %.2f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n",
                   r, sum_eig[0], sum_eig[1], sum_eig[2], n_min, n_max, aniso);

            fprintf(frad, "%.4f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n",
                    r, sum_eig[0], sum_eig[1], sum_eig[2], n_min, n_max, aniso, sum_rho);
        }
    }

    /* ============================================================
       2. 2D slice in xy-plane at z=0
       ============================================================ */
    printf("\n=== 2D xy-slice at z=0 ===\n");
    {
        int k_mid = N / 2;
        int step = 2; /* every other point for manageable output */
        for (int i = 1; i < N-1; i += step) {
            double x = -L + i * dx;
            for (int j = 1; j < N-1; j += step) {
                double y = -L + j * dx;
                int idx = i * NN + j * N + k_mid;
                double b0 = g->phi[0][idx], b1 = g->phi[1][idx], b2 = g->phi[2][idx];

                double H[3][3];
                compute_hessian(b0, b1, b2, mu, kappa, H);
                for (int a = 0; a < 3; a++) H[a][a] += mass2;

                double eig[3];
                eigen3(H, eig);

                double omega_vac = sqrt(k_test*k_test + mass2);
                double n_avg = (sqrt(k_test*k_test+eig[0]) + sqrt(k_test*k_test+eig[1])
                              + sqrt(k_test*k_test+eig[2])) / (3.0 * omega_vac);
                double P = b0 * b1 * b2;

                fprintf(f2d, "%.3f\t%.3f\t%.6f\t%.6f\t%.6f\t%.6e\n",
                        x, y, eig[0], eig[2], n_avg, P);
            }
        }
    }

    /* ============================================================
       3. Z-axis profile (through braid core at x=y=0)
       ============================================================ */
    printf("\n=== Z-axis profile (x=y=0) ===\n");
    printf("  z\tM²_min\tM²_mid\tM²_max\trho\tP\n");
    {
        int ic = N/2, jc = N/2;
        for (int kk = 0; kk < N; kk++) {
            double z = -L + kk * dx;
            int idx = ic * NN + jc * N + kk;
            double b0 = g->phi[0][idx], b1 = g->phi[1][idx], b2 = g->phi[2][idx];

            double H[3][3];
            compute_hessian(b0, b1, b2, mu, kappa, H);
            for (int a = 0; a < 3; a++) H[a][a] += mass2;

            double eig[3];
            eigen3(H, eig);

            double rho = b0*b0 + b1*b1 + b2*b2;
            double P = b0*b1*b2;

            if (kk % 4 == 0)
                printf("  %.2f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4e\n",
                       z, eig[0], eig[1], eig[2], rho, P);

            fprintf(fzax, "%.4f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6e\n",
                    z, eig[0], eig[1], eig[2], rho, P);
        }
    }

    fclose(frad);
    fclose(f2d);
    fclose(fzax);

    /* ============================================================
       Analysis: What does this mean for gravity?
       ============================================================ */
    printf("\n=== ANALYSIS ===\n\n");

    /* Re-scan center for key numbers */
    {
        int ic = N/2, jc = N/2, kc = N/2;
        int idx = ic*NN + jc*N + kc;
        double b0 = g->phi[0][idx], b1 = g->phi[1][idx], b2 = g->phi[2][idx];
        double H[3][3];
        compute_hessian(b0, b1, b2, mu, kappa, H);
        for (int a = 0; a < 3; a++) H[a][a] += mass2;
        double eig[3];
        eigen3(H, eig);

        double omega_vac = sqrt(k_test*k_test + mass2);
        printf("At braid center (x=y=z=0):\n");
        printf("  bg fields: (%.4f, %.4f, %.4f), rho=%.4f, P=%.4e\n",
               b0, b1, b2, b0*b0+b1*b1+b2*b2, b0*b1*b2);
        printf("  Effective mass² eigenvalues: (%.4f, %.4f, %.4f)\n", eig[0], eig[1], eig[2]);
        printf("  Vacuum mass² = %.4f\n", mass2);
        printf("  ΔM² = (%.4f, %.4f, %.4f)\n", eig[0]-mass2, eig[1]-mass2, eig[2]-mass2);
        printf("  n_eff at k=%.1f: (%.4f, %.4f, %.4f)\n", k_test,
               sqrt(k_test*k_test+eig[0])/omega_vac,
               sqrt(k_test*k_test+eig[1])/omega_vac,
               sqrt(k_test*k_test+eig[2])/omega_vac);

        if (eig[2] > mass2 + 0.01)
            printf("  → At least one mode is SLOWER (n>1) → attractive lensing\n");
        if (eig[0] < mass2 - 0.01)
            printf("  → At least one mode is FASTER (n<1) → repulsive lensing\n");
        if (eig[0] < 0)
            printf("  → TACHYONIC mode (M²<0) → instability / spontaneous structure\n");

        double aniso = (eig[2] - eig[0]) / (fabs(eig[1]) + 1e-30);
        printf("  Anisotropy: %.4f", aniso);
        if (aniso > 0.1)
            printf(" → TENSOR structure (not just scalar conformal factor)\n");
        else
            printf(" → approximately scalar (conformal)\n");
    }

    printf("\n  Far field (r >> R_core):\n");
    {
        int i_far = N-5, jc = N/2, kc = N/2;
        int idx = i_far*NN + jc*N + kc;
        double b0 = g->phi[0][idx], b1 = g->phi[1][idx], b2 = g->phi[2][idx];
        double H[3][3];
        compute_hessian(b0, b1, b2, mu, kappa, H);
        for (int a = 0; a < 3; a++) H[a][a] += mass2;
        double eig[3];
        eigen3(H, eig);
        printf("  bg fields: (%.4e, %.4e, %.4e)\n", b0, b1, b2);
        printf("  M² eigenvalues: (%.6f, %.6f, %.6f)\n", eig[0], eig[1], eig[2]);
        printf("  → should be ~%.4f (vacuum)\n", mass2);
    }

    /* ============================================================
       Time-averaged version: the braid oscillates, so the Hessian
       varies in time. Average over one oscillation period.
       ============================================================ */
    printf("\n=== Time-averaged Hessian (T=20, ~5 oscillations) ===\n");
    {
        double T_avg = 20.0;
        int avg_steps = (int)(T_avg / g->dt);
        int n_samples = 0;

        /* Accumulate radial Hessian averages */
        #define NRBIN 30
        double avg_eig[NRBIN][3];
        int avg_count[NRBIN];
        memset(avg_eig, 0, sizeof(avg_eig));
        memset(avg_count, 0, sizeof(avg_count));

        double dr2 = L / NRBIN;
        int sample_every = avg_steps / 20;
        if (sample_every < 1) sample_every = 1;

        for (int s = 0; s < avg_steps; s++) {
            verlet_full_step(g, BIMODAL, 2.25);
            apply_damping_xy(g);

            if (s % sample_every == 0) {
                n_samples++;
                int k_mid = N/2;
                for (int i = 1; i < N-1; i++) {
                    double x = -L + i * dx;
                    for (int j = 1; j < N-1; j++) {
                        double y = -L + j * dx;
                        double rp = sqrt(x*x + y*y);
                        int ir = (int)(rp / dr2);
                        if (ir >= NRBIN) continue;

                        int idx = i*NN + j*N + k_mid;
                        double b0 = g->phi[0][idx], b1 = g->phi[1][idx], b2 = g->phi[2][idx];
                        double H[3][3];
                        compute_hessian(b0, b1, b2, mu, kappa, H);
                        for (int a = 0; a < 3; a++) H[a][a] += mass2;
                        double eig[3];
                        eigen3(H, eig);

                        avg_eig[ir][0] += eig[0];
                        avg_eig[ir][1] += eig[1];
                        avg_eig[ir][2] += eig[2];
                        avg_count[ir]++;
                    }
                }
            }
        }

        printf("  Samples: %d\n", n_samples);
        printf("  r\t<M²_min>\t<M²_mid>\t<M²_max>\t<n_min>\t<n_max>\n");

        FILE *favg = fopen("data/t10b_v3_time_avg.tsv", "w");
        fprintf(favg, "r\tM2_min\tM2_mid\tM2_max\tn_min\tn_max\n");

        double omega_vac = sqrt(k_test*k_test + mass2);
        for (int ir = 0; ir < NRBIN; ir++) {
            if (avg_count[ir] == 0) continue;
            double r = (ir + 0.5) * dr2;
            for (int e = 0; e < 3; e++)
                avg_eig[ir][e] /= avg_count[ir];

            double n_min = sqrt(k_test*k_test + avg_eig[ir][0]) / omega_vac;
            double n_max = sqrt(k_test*k_test + avg_eig[ir][2]) / omega_vac;

            printf("  %.2f\t%.4f\t\t%.4f\t\t%.4f\t\t%.4f\t%.4f\n",
                   r, avg_eig[ir][0], avg_eig[ir][1], avg_eig[ir][2], n_min, n_max);
            fprintf(favg, "%.4f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n",
                    r, avg_eig[ir][0], avg_eig[ir][1], avg_eig[ir][2], n_min, n_max);
        }
        fclose(favg);
    }

    grid_free(g);

    printf("\n=== T10B v3 Complete ===\n");
    printf("Data in data/t10b_v3_*.tsv\n");
    return 0;
}
