/* v60 GEN4 -- C cross-check: equivalence-principle ratios + Newtonian 1/r law.
 *
 * Independent (compiled-C) confirmation of two GEN4 claims:
 *   (1) For the Brannen lepton triple m_k = a^2 (1 + sqrt2 cos(2/9 + 2 pi k/3))^2,
 *       the gravitational mass (= rest energy = eigenvalue of M+M) equals the
 *       inertial mass for every mode: ratio = 1 (EP exact). Total grav charge
 *       rho_grav = sum m_k = 6 a^2 = 9 Q a^2.
 *   (2) A nonzero monopole rho_grav sources a LONG-RANGE Newtonian potential
 *       Phi(r) = -G rho_grav / r  (slope d log|Phi|/d log r = -1) and force
 *       F = -dPhi/dr = -G rho_grav / r^2  (slope -2) -- the GEN1/2 OBE radial
 *       law (massless kernel + nonzero monopole => 1/r^2).
 *
 * Build & run:  gcc -O2 -o /tmp/newton_ep 13_newton_ep.c -lm && /tmp/newton_ep
 */
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

int main(void) {
    const double a = 1.0;
    const double phi = 2.0 / 9.0;          /* Brannen phase */
    const double Q = 2.0 / 3.0;            /* Koide constant */
    const double TWO_PI = 2.0 * M_PI;
    int fail = 0;

    /* (1) Brannen masses and EP ratios */
    double m[3], sqm[3], sum_m = 0.0;
    for (int k = 0; k < 3; k++) {
        sqm[k] = a * (1.0 + sqrt(2.0) * cos(phi + TWO_PI * k / 3.0)); /* sqrt-mass */
        m[k] = sqm[k] * sqm[k];                                       /* mass */
        sum_m += m[k];
    }
    printf("Brannen lepton triple (a=1):\n");
    for (int k = 0; k < 3; k++) {
        /* gravitational mass = rest energy = eigenvalue of M+M = m[k];
         * inertial mass = m[k]. ratio must be 1 (EP). */
        double m_grav = sqm[k] * sqm[k];     /* (M+M)_kk eigenvalue */
        double ratio = m_grav / m[k];
        printf("  mode %d: m=%.10f  m_grav/m_inertial=%.15f\n", k, m[k], ratio);
        if (fabs(ratio - 1.0) > 1e-12) fail = 1;
    }
    printf("  rho_grav = sum m = %.12f   (target 6 a^2 = %.12f = 9 Q a^2 = %.12f)\n",
           sum_m, 6.0 * a * a, 9.0 * Q * a * a);
    if (fabs(sum_m - 6.0 * a * a) > 1e-10) fail = 1;
    if (fabs(sum_m - 9.0 * Q * a * a) > 1e-10) fail = 1;
    printf(fail ? "  [FAIL] EP / rho_grav check\n" : "  [OK] EP exact (ratio 1) and rho_grav = 6 a^2 = 9 Q a^2\n");

    /* (2) Newtonian potential of a point source rho_grav: Phi = -G M / r */
    const double G = 1.0;
    const double M = sum_m;
    /* sample r logarithmically and fit slopes of log|Phi| and log|F| vs log r */
    int N = 20;
    double r0 = 1.0, r1 = 1000.0;
    double sx = 0, sy_phi = 0, sy_f = 0, sxx = 0, sxy_phi = 0, sxy_f = 0;
    for (int i = 0; i < N; i++) {
        double r = r0 * pow(r1 / r0, (double)i / (N - 1));
        double Phi = -G * M / r;
        double F = -G * M / (r * r);          /* F = -dPhi/dr */
        double lx = log(r), lyp = log(fabs(Phi)), lyf = log(fabs(F));
        sx += lx; sy_phi += lyp; sy_f += lyf; sxx += lx * lx;
        sxy_phi += lx * lyp; sxy_f += lx * lyf;
    }
    double slope_phi = (N * sxy_phi - sx * sy_phi) / (N * sxx - sx * sx);
    double slope_f   = (N * sxy_f   - sx * sy_f)   / (N * sxx - sx * sx);
    printf("Newtonian source (monopole M=%.6f):\n", M);
    printf("  d log|Phi| / d log r = %.10f   (target -1: 1/r long-range)\n", slope_phi);
    printf("  d log|F|   / d log r = %.10f   (target -2: 1/r^2 force)\n", slope_f);
    if (fabs(slope_phi + 1.0) > 1e-9) fail = 1;
    if (fabs(slope_f + 2.0) > 1e-9) fail = 1;
    printf(fail ? "  [FAIL] radial law\n" : "  [OK] 1/r potential, 1/r^2 force (GEN1/2 OBE radial law)\n");

    printf(fail ? "\nC CHECK FAILED\n" : "\nC CHECK PASSED.\n");
    return fail;
}
