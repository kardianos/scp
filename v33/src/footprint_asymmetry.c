/*  footprint_asymmetry.c — Measure braid's directional footprint in a ρ gradient
 *
 *  Reads field snapshot, finds braid center, extracts the x-directional
 *  energy density profile through the braid, and measures:
 *    1. Half-width on the high-ρ side (−x) vs low-ρ side (+x)
 *    2. Field-weighted center vs geometric center of footprint
 *    3. Asymmetry ratio R_low / R_high
 *
 *  If R_low > R_high: braid reaches further into depleted side (supports
 *  the "fixed intake amount, variable spatial reach" hypothesis).
 *
 *  Build: gcc -O3 -march=native -o footprint src/footprint_asymmetry.c -lm
 *  Usage: ./footprint data/gtest_quick/field_t0000.bin [A_high A_low]
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define NFIELDS 3

static double compute_rho(double *phi[NFIELDS], long idx) {
    double rho = 0;
    for (int a = 0; a < NFIELDS; a++)
        rho += phi[a][idx] * phi[a][idx];
    return rho;
}

int main(int argc, char **argv) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <snapshot.bin> [A_high A_low]\n", argv[0]);
        return 1;
    }

    double A_high = 0.15, A_low = 0.05;
    if (argc >= 4) { A_high = atof(argv[2]); A_low = atof(argv[3]); }

    /* --- Read snapshot --- */
    FILE *fp = fopen(argv[1], "rb");
    if (!fp) { perror("fopen"); return 1; }

    int N; double L, t;
    fread(&N, sizeof(int), 1, fp);
    fread(&L, sizeof(double), 1, fp);
    fread(&t, sizeof(double), 1, fp);
    long N3 = (long)N * N * N;
    double dx = 2.0 * L / (N - 1);

    printf("Snapshot: N=%d L=%.1f t=%.1f dx=%.4f\n", N, L, t, dx);
    printf("Expected gradient: A_high=%.3f (x=-L) → A_low=%.3f (x=+L)\n", A_high, A_low);

    double *phi[NFIELDS];
    for (int a = 0; a < NFIELDS; a++) {
        phi[a] = malloc(N3 * sizeof(double));
        fread(phi[a], sizeof(double), N3, fp);
    }
    fclose(fp);

    /* --- Find braid center (energy-weighted centroid above 5× average) --- */
    int NN = N * N;
    double avg_rho = 0;
    for (long idx = 0; idx < N3; idx++)
        avg_rho += compute_rho(phi, idx);
    avg_rho /= N3;
    double thresh = 5.0 * avg_rho;

    double cx = 0, cy = 0, cz = 0, wt = 0;
    for (long idx = 0; idx < N3; idx++) {
        double rho = compute_rho(phi, idx);
        if (rho < thresh) continue;
        int i = (int)(idx / NN), j = (int)((idx / N) % N), k = (int)(idx % N);
        double x = -L + i * dx, y = -L + j * dx, z = -L + k * dx;
        cx += x * rho; cy += y * rho; cz += z * rho; wt += rho;
    }
    cx /= wt; cy /= wt; cz /= wt;
    printf("Braid center: (%.3f, %.3f, %.3f)\n\n", cx, cy, cz);

    /* --- Extract x-directional profile through braid center ---
     * Average ρ in a cylinder of radius R_cyl around (cy, cz) */
    double R_cyl = 4.0;  /* should capture braid cross-section */
    int n_xbins = N;
    double *rho_x = calloc(n_xbins, sizeof(double));
    double *rho_bg_x = calloc(n_xbins, sizeof(double));  /* expected background */
    int *counts = calloc(n_xbins, sizeof(int));

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            double y = -L + j * dx;
            double dy = y - cy;
            for (int k = 0; k < N; k++) {
                double z = -L + k * dx;
                double dz = z - cz;
                double r_perp = sqrt(dy * dy + dz * dz);
                if (r_perp > R_cyl) continue;
                long idx = (long)i * NN + j * N + k;
                rho_x[i] += compute_rho(phi, idx);
                counts[i]++;
            }
        }
    }

    /* Compute expected background ρ at each x (from linear gradient) */
    for (int i = 0; i < N; i++) {
        double x = -L + i * dx;
        double frac = (x + L) / (2.0 * L);
        double A_bg = A_high * (1.0 - frac) + A_low * frac;
        /* background ρ ≈ 3 × A_bg² (three fields, each ~A_bg²/2 on average for cos²) */
        rho_bg_x[i] = 1.5 * A_bg * A_bg;
        if (counts[i] > 0)
            rho_x[i] /= counts[i];
    }

    /* --- Compute perturbation: δρ(x) = ρ(x) - ρ_bg(x) --- */
    double *delta_rho = calloc(n_xbins, sizeof(double));
    double max_delta = 0;
    int i_peak = 0;
    for (int i = 0; i < N; i++) {
        delta_rho[i] = rho_x[i] - rho_bg_x[i];
        if (delta_rho[i] > max_delta) {
            max_delta = delta_rho[i];
            i_peak = i;
        }
    }
    double x_peak = -L + i_peak * dx;
    printf("Peak perturbation at x=%.3f (δρ/ρ_bg = %.1f)\n", x_peak, max_delta / rho_bg_x[i_peak]);

    /* --- Measure half-width on each side --- */
    double half_max = 0.5 * max_delta;
    double tenth_max = 0.1 * max_delta;

    /* Search left (toward high ρ, x < x_peak) */
    double x_half_left = x_peak, x_tenth_left = x_peak;
    for (int i = i_peak - 1; i >= 0; i--) {
        if (delta_rho[i] < half_max && x_half_left == x_peak)
            x_half_left = -L + i * dx;
        if (delta_rho[i] < tenth_max && x_tenth_left == x_peak) {
            x_tenth_left = -L + i * dx;
            break;
        }
    }

    /* Search right (toward low ρ, x > x_peak) */
    double x_half_right = x_peak, x_tenth_right = x_peak;
    for (int i = i_peak + 1; i < N; i++) {
        if (delta_rho[i] < half_max && x_half_right == x_peak)
            x_half_right = -L + i * dx;
        if (delta_rho[i] < tenth_max && x_tenth_right == x_peak) {
            x_tenth_right = -L + i * dx;
            break;
        }
    }

    double R_high = fabs(x_peak - x_half_left);   /* half-width toward high ρ */
    double R_low  = fabs(x_half_right - x_peak);   /* half-width toward low ρ */
    double R10_high = fabs(x_peak - x_tenth_left);
    double R10_low  = fabs(x_tenth_right - x_peak);

    printf("\n=== FOOTPRINT ASYMMETRY (relative to perturbation peak) ===\n");
    printf("                   Toward HIGH ρ    Toward LOW ρ     Ratio (low/high)\n");
    printf("Half-width (50%%):  %8.3f          %8.3f          %.3f\n",
           R_high, R_low, R_low / (R_high + 1e-10));
    printf("Tenth-width (10%%): %8.3f          %8.3f          %.3f\n",
           R10_high, R10_low, R10_low / (R10_high + 1e-10));

    /* --- Compute field-weighted center vs geometric center of footprint --- */
    /* Field-weighted center = ∫ x·δρ(x) dx / ∫ δρ(x) dx (over positive δρ) */
    double wx_field = 0, w_field = 0;
    double x_min_foot = L, x_max_foot = -L;
    for (int i = 0; i < N; i++) {
        if (delta_rho[i] > tenth_max) {
            double x = -L + i * dx;
            wx_field += x * delta_rho[i];
            w_field += delta_rho[i];
            if (x < x_min_foot) x_min_foot = x;
            if (x > x_max_foot) x_max_foot = x;
        }
    }
    double x_field = wx_field / (w_field + 1e-30);
    double x_geom = 0.5 * (x_min_foot + x_max_foot);

    printf("\n=== CENTERING ANALYSIS ===\n");
    printf("Footprint extent:  [%.3f, %.3f]  (width=%.3f)\n",
           x_min_foot, x_max_foot, x_max_foot - x_min_foot);
    printf("Geometric center:  %.4f  (midpoint of footprint)\n", x_geom);
    printf("Field-wt center:   %.4f  (∫x·δρ / ∫δρ)\n", x_field);
    printf("Shift (geom→field): %+.4f  (%s)\n",
           x_field - x_geom,
           (x_field - x_geom > 0) ? "field center toward LOW ρ" :
           (x_field - x_geom < 0) ? "field center toward HIGH ρ" : "symmetric");

    printf("\nIf R_low/R_high > 1: braid reaches further into depleted side ✓\n");
    printf("If field center shifted toward HIGH ρ relative to geom center:\n");
    printf("  → braid centers its FIELD (not distance) = energy-minimization picture\n");
    printf("If field center shifted toward LOW ρ relative to geom center:\n");
    printf("  → braid centers its DISTANCE and overshoots = spatial-centering picture\n");

    /* --- Dump x-profile for plotting --- */
    char pfn[512];
    snprintf(pfn, sizeof(pfn), "%s.xprofile.tsv", argv[1]);
    FILE *pf = fopen(pfn, "w");
    fprintf(pf, "x\trho\trho_bg\tdelta_rho\n");
    for (int i = 0; i < N; i++) {
        double x = -L + i * dx;
        fprintf(pf, "%.4f\t%.6e\t%.6e\t%.6e\n", x, rho_x[i], rho_bg_x[i], delta_rho[i]);
    }
    fclose(pf);
    printf("\nProfile saved to %s\n", pfn);

    for (int a = 0; a < NFIELDS; a++) free(phi[a]);
    free(rho_x); free(rho_bg_x); free(delta_rho); free(counts);
    return 0;
}
