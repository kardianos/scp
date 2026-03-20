/*  torsion_profile.c — Experiment 2a: Torsion profile of a settled braid
 *
 *  Reads a settled braid snapshot, computes:
 *    - Jacobian J_ai = dphi_a/dx_i (3x3 matrix) at each point
 *    - det(J) = epsilon_ijk dphi_0/dx_i dphi_1/dx_j dphi_2/dx_k
 *    - Sigma_phi2 = phi_0^2 + phi_1^2 + phi_2^2
 *  Averages in spherical shells around the braid center.
 *
 *  Build: gcc -O3 -march=native -fopenmp -o torsion_profile src/torsion_profile.c -lm
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>

#define NFIELDS 3

int main(int argc, char **argv) {
    const char *infile = "/home/d/code/scp/v34/phonon_test/data/phonon/field_t0200.bin";
    const char *outfile = "data/torsion_profile.tsv";

    if (argc > 1) infile = argv[1];
    if (argc > 2) outfile = argv[2];

    /* Read snapshot */
    FILE *fp = fopen(infile, "rb");
    if (!fp) { fprintf(stderr, "Cannot open %s\n", infile); return 1; }

    int N;
    double L, t;
    fread(&N, sizeof(int), 1, fp);
    fread(&L, sizeof(double), 1, fp);
    fread(&t, sizeof(double), 1, fp);
    printf("Snapshot: N=%d L=%.1f t=%.1f\n", N, L, t);

    long N3 = (long)N * N * N;
    double *phi[NFIELDS];
    for (int a = 0; a < NFIELDS; a++) {
        phi[a] = malloc(N3 * sizeof(double));
        if (!phi[a]) { fprintf(stderr, "Malloc failed for field %d\n", a); return 1; }
        size_t nr = fread(phi[a], sizeof(double), N3, fp);
        if ((long)nr != N3) {
            fprintf(stderr, "Short read field %d: got %zu expected %ld\n", a, nr, N3);
            return 1;
        }
    }
    fclose(fp);
    printf("Loaded %ld doubles per field (%.1f MB total)\n", N3, 3.0*N3*8.0/1e6);

    double dx = 2.0 * L / (N - 1);
    int NN = N * N;

    /* Step 1: Find braid center (energy-weighted centroid above 5x average Sigma_phi2) */
    double avg_phi2 = 0;
    for (long idx = 0; idx < N3; idx++) {
        double p2 = 0;
        for (int a = 0; a < NFIELDS; a++)
            p2 += phi[a][idx] * phi[a][idx];
        avg_phi2 += p2;
    }
    avg_phi2 /= N3;
    double thresh = 5.0 * avg_phi2;
    printf("avg(Sigma_phi2) = %.6e, threshold = %.6e\n", avg_phi2, thresh);

    double cx = 0, cy = 0, cz = 0, wtot = 0;
    for (long idx = 0; idx < N3; idx++) {
        double p2 = 0;
        for (int a = 0; a < NFIELDS; a++)
            p2 += phi[a][idx] * phi[a][idx];
        if (p2 < thresh) continue;

        int i = (int)(idx / NN);
        int j = (int)((idx / N) % N);
        int k = (int)(idx % N);
        double x = -L + i * dx;
        double y = -L + j * dx;
        double z = -L + k * dx;

        cx += x * p2;
        cy += y * p2;
        cz += z * p2;
        wtot += p2;
    }
    if (wtot > 0) { cx /= wtot; cy /= wtot; cz /= wtot; }
    printf("Braid center: (%.3f, %.3f, %.3f)\n", cx, cy, cz);

    /* Step 2: Compute det(J), det(J)^2, Sigma_phi2 in spherical shells */
    double dr = 0.5;
    double R_max = L;  /* max shell radius */
    int nbins = (int)(R_max / dr) + 1;

    double *sum_detJ   = calloc(nbins, sizeof(double));
    double *sum_detJ2  = calloc(nbins, sizeof(double));
    double *sum_phi2   = calloc(nbins, sizeof(double));
    double *sum_grad2  = calloc(nbins, sizeof(double));  /* |grad phi|^2 for comparison */
    long   *counts     = calloc(nbins, sizeof(long));

    printf("Computing torsion profile (%d bins, dr=%.2f)...\n", nbins, dr);

    #pragma omp parallel
    {
        double *loc_detJ  = calloc(nbins, sizeof(double));
        double *loc_detJ2 = calloc(nbins, sizeof(double));
        double *loc_phi2  = calloc(nbins, sizeof(double));
        double *loc_grad2 = calloc(nbins, sizeof(double));
        long   *loc_cnt   = calloc(nbins, sizeof(long));

        #pragma omp for schedule(static)
        for (long idx = 0; idx < N3; idx++) {
            int i = (int)(idx / NN);
            int j = (int)((idx / N) % N);
            int k = (int)(idx % N);

            double x = -L + i * dx;
            double y = -L + j * dx;
            double z = -L + k * dx;

            double r = sqrt((x-cx)*(x-cx) + (y-cy)*(y-cy) + (z-cz)*(z-cz));
            int bin = (int)(r / dr);
            if (bin >= nbins) continue;

            /* Periodic neighbors */
            int ip = (i+1) % N, im = (i-1+N) % N;
            int jp = (j+1) % N, jm = (j-1+N) % N;
            int kp = (k+1) % N, km = (k-1+N) % N;

            /* Compute Jacobian J[a][dir] = dphi_a/dx_dir */
            double J[3][3];
            double g2 = 0;  /* |grad phi|^2 */
            for (int a = 0; a < NFIELDS; a++) {
                J[a][0] = (phi[a][(long)ip*NN + j*N + k] -
                           phi[a][(long)im*NN + j*N + k]) / (2*dx);
                J[a][1] = (phi[a][(long)i*NN + jp*N + k] -
                           phi[a][(long)i*NN + jm*N + k]) / (2*dx);
                J[a][2] = (phi[a][(long)i*NN + j*N + kp] -
                           phi[a][(long)i*NN + j*N + km]) / (2*dx);
                g2 += J[a][0]*J[a][0] + J[a][1]*J[a][1] + J[a][2]*J[a][2];
            }

            /* det(J) = epsilon_ijk J[0][i] J[1][j] J[2][k] */
            double det = J[0][0]*(J[1][1]*J[2][2] - J[1][2]*J[2][1])
                       - J[0][1]*(J[1][0]*J[2][2] - J[1][2]*J[2][0])
                       + J[0][2]*(J[1][0]*J[2][1] - J[1][1]*J[2][0]);

            /* Sigma_phi2 */
            double p2 = 0;
            for (int a = 0; a < NFIELDS; a++)
                p2 += phi[a][idx] * phi[a][idx];

            loc_detJ[bin]  += det;
            loc_detJ2[bin] += det * det;
            loc_phi2[bin]  += p2;
            loc_grad2[bin] += g2;
            loc_cnt[bin]++;
        }

        #pragma omp critical
        {
            for (int b = 0; b < nbins; b++) {
                sum_detJ[b]  += loc_detJ[b];
                sum_detJ2[b] += loc_detJ2[b];
                sum_phi2[b]  += loc_phi2[b];
                sum_grad2[b] += loc_grad2[b];
                counts[b]    += loc_cnt[b];
            }
        }

        free(loc_detJ); free(loc_detJ2); free(loc_phi2); free(loc_grad2); free(loc_cnt);
    }

    /* Compute background Sigma_phi2 from outermost non-empty shell */
    double phi2_bg = 0;
    int bg_count = 0;
    for (int b = nbins-1; b >= nbins-5; b--) {
        if (counts[b] > 0) {
            phi2_bg += sum_phi2[b] / counts[b];
            bg_count++;
        }
    }
    if (bg_count > 0) phi2_bg /= bg_count;
    printf("Background Sigma_phi2 = %.6e\n", phi2_bg);

    /* Write output */
    FILE *out = fopen(outfile, "w");
    if (!out) { fprintf(stderr, "Cannot open %s for writing\n", outfile); return 1; }
    fprintf(out, "# Torsion profile of settled braid (snapshot t=%.1f)\n", t);
    fprintf(out, "# Braid center: (%.3f, %.3f, %.3f)\n", cx, cy, cz);
    fprintf(out, "# Background Sigma_phi2 = %.6e\n", phi2_bg);
    fprintf(out, "r\tdetJ_mean\tdetJ2_mean\tphi2_mean\tdelta_phi2\tgrad2_mean\tcounts\n");

    for (int b = 0; b < nbins; b++) {
        if (counts[b] == 0) continue;
        double r = (b + 0.5) * dr;
        double dJ_mean  = sum_detJ[b]  / counts[b];
        double dJ2_mean = sum_detJ2[b] / counts[b];
        double p2_mean  = sum_phi2[b]  / counts[b];
        double g2_mean  = sum_grad2[b] / counts[b];
        double dp2      = p2_mean - phi2_bg;

        fprintf(out, "%.3f\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%ld\n",
                r, dJ_mean, dJ2_mean, p2_mean, dp2, g2_mean, counts[b]);
    }
    fclose(out);
    printf("Output written to %s\n", outfile);

    /* Print summary */
    printf("\n=== Torsion Profile Summary ===\n");
    printf("%-8s %-14s %-14s %-14s %-14s\n",
           "r", "<det(J)>", "<det(J)^2>", "delta(phi2)", "<grad2>");
    for (int b = 0; b < nbins && b < 40; b++) {
        if (counts[b] == 0) continue;
        double r = (b + 0.5) * dr;
        double dJ_mean  = sum_detJ[b]  / counts[b];
        double dJ2_mean = sum_detJ2[b] / counts[b];
        double dp2      = sum_phi2[b] / counts[b] - phi2_bg;
        double g2_mean  = sum_grad2[b] / counts[b];
        printf("%-8.2f %-14.6e %-14.6e %-14.6e %-14.6e\n",
               r, dJ_mean, dJ2_mean, dp2, g2_mean);
    }

    /* Cleanup */
    for (int a = 0; a < NFIELDS; a++) free(phi[a]);
    free(sum_detJ); free(sum_detJ2); free(sum_phi2); free(sum_grad2); free(counts);

    return 0;
}
