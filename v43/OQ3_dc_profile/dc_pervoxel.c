/*
 * dc_pervoxel.c — Per-voxel DC extraction for OQ3 corrected analysis
 *
 * Accumulates per-VOXEL time averages of theta_a across all frames,
 * then computes DC magnitude and AC RMS at each voxel, then bins
 * into radial shells to get DC_rms(r), AC_rms(r), DC_frac(r).
 *
 * This is the CORRECT methodology: DC is computed per-voxel FIRST,
 * then statistics are taken over radial shells. The original dc_profile.c
 * shell-averaged the signed theta before computing DC, which cancels
 * dipolar structure.
 *
 * Build: gcc -O3 -o dc_pervoxel dc_pervoxel.c -lzstd -lm
 */

#define SFA_IMPLEMENTATION
#include "../../sfa/format/sfa.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_RBINS 50
#define R_MAX 12.0
#define DR 0.5
#define TRANSIENT_FRAMES 25  /* skip first 25 frames (t < ~5) for no-transient analysis */

int main(int argc, char **argv)
{
    const char *sfa_path = "/home/d/code/scp/v34/torsion_coupling/data/braid_hires.sfa";
    if (argc > 1) sfa_path = argv[1];

    SFA *s = sfa_open(sfa_path);
    if (!s) { fprintf(stderr, "Failed to open %s\n", sfa_path); return 1; }

    uint32_t N = s->Nx;
    double L = s->Lx;
    /* Vertex-centered grid: dx = 2*L/(N-1), matching the simulation kernel */
    double dx = 2.0 * L / (N - 1);
    uint64_t N3 = (uint64_t)N * N * N;
    int nframes = s->total_frames;

    printf("SFA: N=%u, L=%.1f, dx=%.4f (vertex-centered), frames=%d, cols=%u\n",
           N, L, dx, nframes, s->n_columns);
    for (uint32_t c = 0; c < s->n_columns; c++)
        printf("  col %u: %s\n", c, s->columns[c].name);

    int nbins = (int)(R_MAX / DR);
    if (nbins > MAX_RBINS) nbins = MAX_RBINS;
    printf("Radial binning: dr=%.1f, r_max=%.1f, nbins=%d\n", DR, R_MAX, nbins);

    /* Allocate per-voxel accumulators: sum and sum-of-squares for theta_{x,y,z} */
    /* Full analysis (all frames) */
    double *sum_tx  = (double *)calloc(N3, sizeof(double));
    double *sum_ty  = (double *)calloc(N3, sizeof(double));
    double *sum_tz  = (double *)calloc(N3, sizeof(double));
    double *sum2_tx = (double *)calloc(N3, sizeof(double));
    double *sum2_ty = (double *)calloc(N3, sizeof(double));
    double *sum2_tz = (double *)calloc(N3, sizeof(double));
    /* No-transient analysis (skip first TRANSIENT_FRAMES frames) */
    double *nt_sum_tx  = (double *)calloc(N3, sizeof(double));
    double *nt_sum_ty  = (double *)calloc(N3, sizeof(double));
    double *nt_sum_tz  = (double *)calloc(N3, sizeof(double));
    double *nt_sum2_tx = (double *)calloc(N3, sizeof(double));
    double *nt_sum2_ty = (double *)calloc(N3, sizeof(double));
    double *nt_sum2_tz = (double *)calloc(N3, sizeof(double));

    if (!sum_tx || !sum_ty || !sum_tz || !sum2_tx || !sum2_ty || !sum2_tz ||
        !nt_sum_tx || !nt_sum_ty || !nt_sum_tz || !nt_sum2_tx || !nt_sum2_ty || !nt_sum2_tz) {
        fprintf(stderr, "Failed to allocate per-voxel accumulators\n");
        return 1;
    }

    /* Allocate frame buffer */
    double *buf = (double *)malloc(s->frame_bytes);
    if (!buf) { fprintf(stderr, "malloc frame buffer failed\n"); return 1; }

    /* --- Pass 1: Compute centroid from frame 0 phi fields --- */
    if (sfa_read_frame(s, 0, buf) != 0) {
        fprintf(stderr, "Failed to read frame 0\n"); return 1;
    }
    double *phi_x = buf + 0 * N3;
    double *phi_y = buf + 1 * N3;
    double *phi_z = buf + 2 * N3;

    double cx = 0, cy = 0, cz = 0, wtot = 0;
    for (uint32_t iz = 0; iz < N; iz++)
    for (uint32_t iy = 0; iy < N; iy++)
    for (uint32_t ix = 0; ix < N; ix++) {
        uint64_t idx = (uint64_t)iz * N * N + iy * N + ix;
        double P = fabs(phi_x[idx] * phi_y[idx] * phi_z[idx]);
        double x = -L + ix * dx;
        double y = -L + iy * dx;
        double z = -L + iz * dx;
        cx += P * x;
        cy += P * y;
        cz += P * z;
        wtot += P;
    }
    if (wtot > 0) { cx /= wtot; cy /= wtot; cz /= wtot; }
    printf("Centroid from frame 0: (%.3f, %.3f, %.3f), |P|_tot = %.4e\n",
           cx, cy, cz, wtot);

    /* Check centroid at mid and last frame */
    int check_frames[] = {nframes/2, nframes-1};
    for (int cf = 0; cf < 2; cf++) {
        int fr = check_frames[cf];
        if (sfa_read_frame(s, fr, buf) != 0) continue;
        phi_x = buf + 0 * N3;
        phi_y = buf + 1 * N3;
        phi_z = buf + 2 * N3;
        double ccx = 0, ccy = 0, ccz = 0, ww = 0;
        for (uint32_t iz = 0; iz < N; iz++)
        for (uint32_t iy = 0; iy < N; iy++)
        for (uint32_t ix = 0; ix < N; ix++) {
            uint64_t idx = (uint64_t)iz * N * N + iy * N + ix;
            double P = fabs(phi_x[idx] * phi_y[idx] * phi_z[idx]);
            double x = -L + ix * dx;
            double y = -L + iy * dx;
            double z = -L + iz * dx;
            ccx += P * x; ccy += P * y; ccz += P * z;
            ww += P;
        }
        if (ww > 0) { ccx /= ww; ccy /= ww; ccz /= ww; }
        double drift = sqrt((ccx-cx)*(ccx-cx) + (ccy-cy)*(ccy-cy) + (ccz-cz)*(ccz-cz));
        printf("Centroid at frame %d: (%.3f, %.3f, %.3f), drift=%.4f (%.2f grid cells)\n",
               fr, ccx, ccy, ccz, drift, drift / dx);
    }

    /* --- Pass 2: Accumulate per-voxel theta sums --- */
    printf("\nAccumulating per-voxel theta across %d frames...\n", nframes);
    int nt_count = 0;  /* number of non-transient frames */
    for (int frame = 0; frame < nframes; frame++) {
        if (sfa_read_frame(s, frame, buf) != 0) {
            fprintf(stderr, "Failed to read frame %d\n", frame);
            continue;
        }

        double *theta_x = buf + 3 * N3;
        double *theta_y = buf + 4 * N3;
        double *theta_z = buf + 5 * N3;

        for (uint64_t idx = 0; idx < N3; idx++) {
            double tx = theta_x[idx];
            double ty = theta_y[idx];
            double tz = theta_z[idx];
            sum_tx[idx] += tx;
            sum_ty[idx] += ty;
            sum_tz[idx] += tz;
            sum2_tx[idx] += tx * tx;
            sum2_ty[idx] += ty * ty;
            sum2_tz[idx] += tz * tz;
        }

        if (frame >= TRANSIENT_FRAMES) {
            nt_count++;
            for (uint64_t idx = 0; idx < N3; idx++) {
                double tx = theta_x[idx];
                double ty = theta_y[idx];
                double tz = theta_z[idx];
                nt_sum_tx[idx] += tx;
                nt_sum_ty[idx] += ty;
                nt_sum_tz[idx] += tz;
                nt_sum2_tx[idx] += tx * tx;
                nt_sum2_ty[idx] += ty * ty;
                nt_sum2_tz[idx] += tz * tz;
            }
        }

        if (frame % 50 == 0 || frame == nframes - 1)
            printf("  frame %d/%d\n", frame, nframes);
    }

    printf("Total frames: %d, non-transient frames: %d (skipped first %d)\n",
           nframes, nt_count, TRANSIENT_FRAMES);

    /* --- Compute per-voxel DC and AC, then bin into radial shells --- */
    /* Bins: DC_mag^2 sum, AC_rms^2 sum, count */
    double dc2_bin[MAX_RBINS], ac2_bin[MAX_RBINS], bin_count[MAX_RBINS];
    double nt_dc2_bin[MAX_RBINS], nt_ac2_bin[MAX_RBINS];
    /* Also accumulate total energy for global DC fraction */
    double total_dc_energy = 0, total_theta_energy = 0;
    double nt_total_dc_energy = 0, nt_total_theta_energy = 0;

    memset(dc2_bin, 0, sizeof(dc2_bin));
    memset(ac2_bin, 0, sizeof(ac2_bin));
    memset(bin_count, 0, sizeof(bin_count));
    memset(nt_dc2_bin, 0, sizeof(nt_dc2_bin));
    memset(nt_ac2_bin, 0, sizeof(nt_ac2_bin));

    for (uint32_t iz = 0; iz < N; iz++)
    for (uint32_t iy = 0; iy < N; iy++)
    for (uint32_t ix = 0; ix < N; ix++) {
        uint64_t idx = (uint64_t)iz * N * N + iy * N + ix;
        double x = -L + ix * dx - cx;
        double y = -L + iy * dx - cy;
        double z = -L + iz * dx - cz;
        double r = sqrt(x*x + y*y + z*z);

        int bin = (int)(r / DR);
        if (bin >= nbins) continue;  /* r >= R_MAX */

        /* Full analysis */
        double dcx = sum_tx[idx] / nframes;
        double dcy = sum_ty[idx] / nframes;
        double dcz = sum_tz[idx] / nframes;
        double dc_mag2 = dcx*dcx + dcy*dcy + dcz*dcz;

        double ac_var_x = sum2_tx[idx] / nframes - dcx*dcx;
        double ac_var_y = sum2_ty[idx] / nframes - dcy*dcy;
        double ac_var_z = sum2_tz[idx] / nframes - dcz*dcz;
        /* Clamp negative variances from floating point */
        if (ac_var_x < 0) ac_var_x = 0;
        if (ac_var_y < 0) ac_var_y = 0;
        if (ac_var_z < 0) ac_var_z = 0;
        double ac_rms2 = ac_var_x + ac_var_y + ac_var_z;

        dc2_bin[bin] += dc_mag2;
        ac2_bin[bin] += ac_rms2;
        bin_count[bin] += 1.0;

        total_dc_energy += dc_mag2;
        total_theta_energy += dc_mag2 + ac_rms2;

        /* No-transient analysis */
        double nt_dcx = nt_sum_tx[idx] / nt_count;
        double nt_dcy = nt_sum_ty[idx] / nt_count;
        double nt_dcz = nt_sum_tz[idx] / nt_count;
        double nt_dc_mag2 = nt_dcx*nt_dcx + nt_dcy*nt_dcy + nt_dcz*nt_dcz;

        double nt_ac_var_x = nt_sum2_tx[idx] / nt_count - nt_dcx*nt_dcx;
        double nt_ac_var_y = nt_sum2_ty[idx] / nt_count - nt_dcy*nt_dcy;
        double nt_ac_var_z = nt_sum2_tz[idx] / nt_count - nt_dcz*nt_dcz;
        if (nt_ac_var_x < 0) nt_ac_var_x = 0;
        if (nt_ac_var_y < 0) nt_ac_var_y = 0;
        if (nt_ac_var_z < 0) nt_ac_var_z = 0;
        double nt_ac_rms2 = nt_ac_var_x + nt_ac_var_y + nt_ac_var_z;

        nt_dc2_bin[bin] += nt_dc_mag2;
        nt_ac2_bin[bin] += nt_ac_rms2;

        nt_total_dc_energy += nt_dc_mag2;
        nt_total_theta_energy += nt_dc_mag2 + nt_ac_rms2;
    }

    /* --- Write output TSV --- */
    const char *outpath = "dc_pervoxel.tsv";
    FILE *fp = fopen(outpath, "w");
    if (!fp) { fprintf(stderr, "Cannot open %s\n", outpath); return 1; }

    fprintf(fp, "# Per-voxel DC extraction (OQ3 corrected)\n");
    fprintf(fp, "# Source: %s (%d frames, N=%u, L=%.1f, dx=%.4f vertex-centered)\n",
            sfa_path, nframes, N, L, dx);
    fprintf(fp, "# Full: all %d frames; No-transient: frames %d-%d (%d frames)\n",
            nframes, TRANSIENT_FRAMES, nframes-1, nt_count);
    fprintf(fp, "# DC_rms = sqrt(<DC_mag^2>_shell), AC_rms = sqrt(<AC_rms^2>_shell)\n");
    fprintf(fp, "# DC_frac = DC_rms / sqrt(DC_rms^2 + AC_rms^2)\n");
    fprintf(fp, "r\tDC_rms\tAC_rms\tDC_frac\tDC_rms_notrans\tDC_frac_notrans\tnvoxels\n");

    printf("\n--- Per-Voxel DC Radial Profile ---\n");
    printf("%6s  %12s  %12s  %10s  %12s  %10s  %8s\n",
           "r", "DC_rms", "AC_rms", "DC_frac", "DC_rms_nt", "DC_frac_nt", "nvoxels");

    for (int b = 0; b < nbins; b++) {
        double r_mid = (b + 0.5) * DR;
        if (bin_count[b] < 1.0) {
            fprintf(fp, "%.4f\t0\t0\t0\t0\t0\t0\n", r_mid);
            continue;
        }

        double dc_rms = sqrt(dc2_bin[b] / bin_count[b]);
        double ac_rms = sqrt(ac2_bin[b] / bin_count[b]);
        double total_rms = sqrt(dc_rms*dc_rms + ac_rms*ac_rms);
        double dc_frac = (total_rms > 0) ? dc_rms / total_rms : 0;

        double nt_dc_rms = sqrt(nt_dc2_bin[b] / bin_count[b]);
        double nt_ac_rms = sqrt(nt_ac2_bin[b] / bin_count[b]);
        double nt_total_rms = sqrt(nt_dc_rms*nt_dc_rms + nt_ac_rms*nt_ac_rms);
        double nt_dc_frac = (nt_total_rms > 0) ? nt_dc_rms / nt_total_rms : 0;

        fprintf(fp, "%.4f\t%.6e\t%.6e\t%.6f\t%.6e\t%.6f\t%.0f\n",
                r_mid, dc_rms, ac_rms, dc_frac, nt_dc_rms, nt_dc_frac, bin_count[b]);

        printf("%6.2f  %12.6e  %12.6e  %10.6f  %12.6e  %10.6f  %8.0f\n",
               r_mid, dc_rms, ac_rms, dc_frac, nt_dc_rms, nt_dc_frac, bin_count[b]);
    }

    fclose(fp);

    /* --- Summary at specific radii --- */
    printf("\n--- Summary: DC_frac at selected radii ---\n");
    double target_r[] = {2.0, 5.0, 8.0, 10.0};
    int ntargets = 4;
    for (int t = 0; t < ntargets; t++) {
        int b = (int)(target_r[t] / DR);
        if (b >= nbins || bin_count[b] < 1) {
            printf("  r=%5.1f: no data\n", target_r[t]);
            continue;
        }
        double dc_rms = sqrt(dc2_bin[b] / bin_count[b]);
        double ac_rms = sqrt(ac2_bin[b] / bin_count[b]);
        double total_rms = sqrt(dc_rms*dc_rms + ac_rms*ac_rms);
        double dc_frac = (total_rms > 0) ? dc_rms / total_rms : 0;

        double nt_dc_rms = sqrt(nt_dc2_bin[b] / bin_count[b]);
        double nt_ac_rms = sqrt(nt_ac2_bin[b] / bin_count[b]);
        double nt_total_rms = sqrt(nt_dc_rms*nt_dc_rms + nt_ac_rms*nt_ac_rms);
        double nt_dc_frac = (nt_total_rms > 0) ? nt_dc_rms / nt_total_rms : 0;

        printf("  r=%5.1f:  DC_frac=%.4f (%.2f%%)  no-transient=%.4f (%.2f%%)\n",
               target_r[t], dc_frac, dc_frac*100.0, nt_dc_frac, nt_dc_frac*100.0);
    }

    /* Global DC fraction within r < R_MAX */
    double global_dc_frac = (total_theta_energy > 0)
        ? sqrt(total_dc_energy / total_theta_energy) : 0;
    double nt_global_dc_frac = (nt_total_theta_energy > 0)
        ? sqrt(nt_total_dc_energy / nt_total_theta_energy) : 0;

    printf("\n--- Global DC fraction (r < %.0f) ---\n", R_MAX);
    printf("  All frames:      sqrt(DC_energy/total_energy) = %.4f (%.2f%%)\n",
           global_dc_frac, global_dc_frac*100.0);
    printf("  No transient:    sqrt(DC_energy/total_energy) = %.4f (%.2f%%)\n",
           nt_global_dc_frac, nt_global_dc_frac*100.0);

    printf("\n--- Comparison with V34 claim ---\n");
    printf("V34 claimed DC ~ 0.2%% of total theta amplitude.\n");
    printf("Original dc_profile.c (shell-averaged DC) found 5-12%% in core.\n");
    printf("This per-voxel analysis finds the values above.\n");
    printf("If DC_frac is constant with r, it may be a systematic artifact.\n");
    printf("If DC_frac decays with r, it has physical radial structure.\n");

    printf("\nWrote: %s\n", outpath);

    /* Cleanup */
    free(buf);
    free(sum_tx); free(sum_ty); free(sum_tz);
    free(sum2_tx); free(sum2_ty); free(sum2_tz);
    free(nt_sum_tx); free(nt_sum_ty); free(nt_sum_tz);
    free(nt_sum2_tx); free(nt_sum2_ty); free(nt_sum2_tz);
    sfa_close(s);
    return 0;
}
