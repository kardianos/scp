/*
 * poisson_phi.c — Solve ∇²Φ = α·ρ for the oscillon's static gravitational potential
 *
 * Reads a radial profile (r, phi1, phi2, phi3, rho_E) from TSV.
 * Solves the spherical Poisson equation:
 *   Φ'' + (2/r)Φ' = α·ρ(r)
 * with BC: Φ'(0)=0, Φ(R_max) = -(α·Q)/(4π·R_max)  where Q = ∫ρ·4πr²dr
 *
 * Outputs: r, Φ(r), Φ_1overr(r) = -(αQ)/(4πr) for comparison
 *
 * Compile: gcc -O3 -Wall -o poisson_phi v21/src/poisson_phi.c -lm
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAXPTS 500

int main(int argc, char **argv)
{
    const char *infile = "v21/data/triad3d_test1_profile_t1000.tsv";
    double alpha = 1.0;  /* coupling (will rescale afterward) */
    int Nr = 2000;       /* integration grid (finer than profile) */

    for (int i = 1; i < argc - 1; i += 2) {
        if      (!strcmp(argv[i], "-i"))     infile = argv[i+1];
        else if (!strcmp(argv[i], "-alpha")) alpha  = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-Nr"))    Nr     = atoi(argv[i+1]);
    }

    /* Read profile */
    FILE *fin = fopen(infile, "r");
    if (!fin) { fprintf(stderr, "Cannot open %s\n", infile); return 1; }

    double r_prof[MAXPTS], rho_prof[MAXPTS];
    int np = 0;
    char line[1024];
    fgets(line, sizeof(line), fin); /* skip header */
    while (fgets(line, sizeof(line), fin) && np < MAXPTS) {
        double r, p1, p2, p3, rho;
        if (sscanf(line, "%lf %lf %lf %lf %lf", &r, &p1, &p2, &p3, &rho) == 5) {
            r_prof[np] = r;
            rho_prof[np] = rho;
            np++;
        }
    }
    fclose(fin);

    if (np < 5) { fprintf(stderr, "Too few profile points: %d\n", np); return 1; }

    double R_max = r_prof[np-1];
    double dr = R_max / (Nr - 1);

    /* Interpolate ρ onto fine grid */
    double *rho = calloc(Nr, sizeof(double));
    for (int i = 0; i < Nr; i++) {
        double r = i * dr;
        /* Linear interpolation in profile data */
        if (r <= r_prof[0]) {
            rho[i] = rho_prof[0];
        } else if (r >= r_prof[np-1]) {
            rho[i] = 0.0;
        } else {
            int j = 0;
            while (j < np-1 && r_prof[j+1] < r) j++;
            double frac = (r - r_prof[j]) / (r_prof[j+1] - r_prof[j]);
            rho[i] = rho_prof[j] + frac * (rho_prof[j+1] - rho_prof[j]);
        }
    }

    /* Compute total source: Q = ∫ρ·4πr²dr */
    double Q = 0;
    for (int i = 0; i < Nr; i++) {
        double r = i * dr;
        if (i == 0) continue;
        Q += rho[i] * 4.0 * M_PI * r * r * dr;
    }
    printf("# Total source Q = integral(rho * 4pi r^2 dr) = %.6f\n", Q);
    printf("# alpha = %.6e\n", alpha);
    printf("# Expected far-field: Phi(r) -> -(alpha*Q)/(4*pi*r) = %.6e / r\n",
           -alpha * Q / (4.0 * M_PI));

    /* Solve Poisson equation using the standard method:
     * Φ(r) = -(α/r) ∫₀ʳ ρ(r')r'² dr' - α ∫ᵣ^∞ ρ(r')r' dr'
     * (Green's function for spherical Poisson equation)
     */
    double *Phi = calloc(Nr, sizeof(double));

    /* Build cumulative integrals */
    /* I1(r) = ∫₀ʳ ρ(r') r'² dr' */
    double *I1 = calloc(Nr, sizeof(double));
    I1[0] = 0;
    for (int i = 1; i < Nr; i++) {
        double r = i * dr;
        I1[i] = I1[i-1] + rho[i] * r * r * dr;
    }

    /* I2(r) = ∫ᵣ^∞ ρ(r') r' dr' */
    double *I2 = calloc(Nr, sizeof(double));
    I2[Nr-1] = 0;
    for (int i = Nr-2; i >= 0; i--) {
        double r = (i+1) * dr;
        I2[i] = I2[i+1] + rho[i+1] * r * dr;
    }

    /* Φ(r) = -α [I1(r)/r + I2(r)] */
    Phi[0] = -alpha * I2[0];  /* at r=0: I1/r → 0 (L'Hôpital), only I2 term */
    for (int i = 1; i < Nr; i++) {
        double r = i * dr;
        Phi[i] = -alpha * (I1[i] / r + I2[i]);
    }

    /* Output */
    printf("# r\tPhi\tPhi_1overr\trho\tr*Phi\n");
    for (int i = 0; i < Nr; i++) {
        double r = i * dr;
        double Phi_1r = (r > 0.01) ? -alpha * Q / (4.0 * M_PI * r) : Phi[0];
        printf("%.6f\t%.6e\t%.6e\t%.6e\t%.6e\n",
               r, Phi[i], Phi_1r, rho[i],
               (r > 0.01) ? r * Phi[i] : 0.0);
    }

    /* Summary */
    fprintf(stderr, "Source: Q = %.4f (total energy in profile)\n", Q);
    fprintf(stderr, "Phi(0) = %.6e\n", Phi[0]);
    fprintf(stderr, "Phi(R/2) = %.6e, -(aQ)/(4pi*R/2) = %.6e\n",
            Phi[Nr/2], -alpha*Q/(4.0*M_PI*R_max/2));
    fprintf(stderr, "Phi(R) = %.6e, -(aQ)/(4pi*R) = %.6e\n",
            Phi[Nr-1], -alpha*Q/(4.0*M_PI*R_max));
    fprintf(stderr, "r*Phi at R/2 = %.6e, at R = %.6e (should be const = %.6e)\n",
            (R_max/2)*Phi[Nr/2], R_max*Phi[Nr-1], -alpha*Q/(4.0*M_PI));

    free(rho); free(Phi); free(I1); free(I2);
    return 0;
}
