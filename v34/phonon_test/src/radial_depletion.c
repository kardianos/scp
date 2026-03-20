/*  radial_depletion.c — Measure spherical radial depletion profile around a braid
 *
 *  Reads a field snapshot from v33_G, finds the braid center via
 *  energy-weighted centroid, then computes rho(r) = <sum phi_a^2>
 *  averaged in spherical shells.
 *
 *  Build: gcc -O3 -march=native -o radial_depletion src/radial_depletion.c -lm
 *  Usage: ./radial_depletion <snapshot.bin> [r_max] [dr]
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define NFIELDS 3

int main(int argc, char **argv) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <snapshot.bin> [r_max=50] [dr=0.5]\n", argv[0]);
        return 1;
    }

    const char *fname = argv[1];
    double r_max = (argc > 2) ? atof(argv[2]) : 50.0;
    double dr    = (argc > 3) ? atof(argv[3]) : 0.5;

    /* Read snapshot header */
    FILE *fp = fopen(fname, "rb");
    if (!fp) { perror("fopen"); return 1; }

    int N;
    double L, t;
    fread(&N, sizeof(int), 1, fp);
    fread(&L, sizeof(double), 1, fp);
    fread(&t, sizeof(double), 1, fp);

    long N3 = (long)N * N * N;
    double dx = 2.0 * L / (N - 1);
    int NN = N * N;

    printf("Snapshot: N=%d L=%.1f t=%.1f dx=%.4f\n", N, L, t, dx);
    printf("Grid: %d^3 = %ld points, box [-%.0f, %.0f]\n", N, N3, L, L);

    /* Allocate and read fields */
    double **phi = malloc(NFIELDS * sizeof(double*));
    for (int a = 0; a < NFIELDS; a++) {
        phi[a] = malloc(N3 * sizeof(double));
        size_t nr = fread(phi[a], sizeof(double), N3, fp);
        if ((long)nr != N3) {
            fprintf(stderr, "ERROR: field %d: read %zu of %ld\n", a, nr, N3);
            return 1;
        }
    }
    fclose(fp);
    printf("Fields loaded: %.2f GB\n", 3.0 * N3 * 8.0 / 1e9);

    /* Step 1: Compute phi^2 at each point and find average */
    double *phi2 = malloc(N3 * sizeof(double));
    double avg_phi2 = 0;
    for (long idx = 0; idx < N3; idx++) {
        double p2 = 0;
        for (int a = 0; a < NFIELDS; a++)
            p2 += phi[a][idx] * phi[a][idx];
        phi2[idx] = p2;
        avg_phi2 += p2;
    }
    avg_phi2 /= N3;
    printf("Average phi^2 = %.6e\n", avg_phi2);

    /* Step 2: Find braid center (energy-weighted centroid, threshold 5x average) */
    double thresh = 5.0 * avg_phi2;
    double wx = 0, wy = 0, wz = 0, wtot = 0;
    long n_above = 0;

    for (int i = 0; i < N; i++) {
        double x = -L + i * dx;
        for (int j = 0; j < N; j++) {
            double y = -L + j * dx;
            for (int k = 0; k < N; k++) {
                long idx = (long)i * NN + j * N + k;
                if (phi2[idx] > thresh) {
                    double z = -L + k * dx;
                    double w = phi2[idx];
                    wx += w * x;
                    wy += w * y;
                    wz += w * z;
                    wtot += w;
                    n_above++;
                }
            }
        }
    }

    double cx = 0, cy = 0, cz = 0;
    if (wtot > 0) {
        cx = wx / wtot;
        cy = wy / wtot;
        cz = wz / wtot;
    }
    printf("Braid center: (%.3f, %.3f, %.3f)  [%ld pts above threshold]\n",
           cx, cy, cz, n_above);

    /* Step 3: Compute spherical shell averages */
    int nbins = (int)(r_max / dr) + 1;
    double *rho_sum = calloc(nbins, sizeof(double));
    long   *counts  = calloc(nbins, sizeof(long));

    for (int i = 0; i < N; i++) {
        double x = -L + i * dx - cx;
        for (int j = 0; j < N; j++) {
            double y = -L + j * dx - cy;
            for (int k = 0; k < N; k++) {
                double z = -L + k * dx - cz;
                double r = sqrt(x*x + y*y + z*z);
                int b = (int)(r / dr);
                if (b >= nbins) continue;
                long idx = (long)i * NN + j * N + k;
                rho_sum[b] += phi2[idx];
                counts[b]++;
            }
        }
    }

    /* Step 4: Compute background from far field (r > 40) */
    double bg_sum = 0;
    long bg_count = 0;
    int bg_start = (int)(40.0 / dr);
    for (int b = bg_start; b < nbins; b++) {
        if (counts[b] > 0) {
            bg_sum += rho_sum[b];  /* sum of phi2 values */
            bg_count += counts[b];
        }
    }
    double rho_bg = (bg_count > 0) ? bg_sum / bg_count : 0;
    printf("Background rho (r>40): %.6e  (%ld points)\n", rho_bg, bg_count);

    /* Step 5: Output radial profile */
    char outfn[512];
    snprintf(outfn, sizeof(outfn), "data/depletion_t%04d.tsv", (int)(t + 0.5));
    FILE *fout = fopen(outfn, "w");
    if (!fout) { perror("fopen output"); return 1; }

    fprintf(fout, "# Radial depletion profile at t=%.1f\n", t);
    fprintf(fout, "# Center: (%.3f, %.3f, %.3f)\n", cx, cy, cz);
    fprintf(fout, "# Background rho (r>40): %.6e\n", rho_bg);
    fprintf(fout, "r\trho\tdelta_rho\tcounts\n");

    printf("\n%8s %14s %14s %10s\n", "r", "rho(r)", "delta_rho", "counts");
    printf("-----------------------------------------------------\n");

    for (int b = 0; b < nbins; b++) {
        double r = (b + 0.5) * dr;
        double rho = (counts[b] > 0) ? rho_sum[b] / counts[b] : 0;
        double delta = rho - rho_bg;

        fprintf(fout, "%.2f\t%.8e\t%+.8e\t%ld\n", r, rho, delta, counts[b]);

        if (b % 2 == 0 || b < 20) {
            printf("%8.2f %14.6e %14.6e %10ld\n", r, rho, delta, counts[b]);
        }
    }

    fclose(fout);
    printf("\nOutput written to %s\n", outfn);

    /* Cleanup */
    for (int a = 0; a < NFIELDS; a++) free(phi[a]);
    free(phi);
    free(phi2);
    free(rho_sum);
    free(counts);

    return 0;
}
