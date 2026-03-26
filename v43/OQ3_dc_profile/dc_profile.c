/*
 * dc_profile.c — Compute time-averaged signed theta radial profile
 * from V34 braid_hires.sfa (6-field Cosserat, N=80, L=25, 264 frames)
 *
 * Build: gcc -O3 -o dc_profile dc_profile.c -lzstd -lm
 */

#define SFA_IMPLEMENTATION
#include "../../sfa/format/sfa.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define NBINS 50

int main(int argc, char **argv)
{
    const char *sfa_path = "/home/d/code/scp/v34/torsion_coupling/data/braid_hires.sfa";
    if (argc > 1) sfa_path = argv[1];

    SFA *s = sfa_open(sfa_path);
    if (!s) { fprintf(stderr, "Failed to open %s\n", sfa_path); return 1; }

    uint32_t N = s->Nx;
    double L = s->Lx;
    double h = 2.0 * L / N;
    uint64_t N3 = (uint64_t)N * N * N;
    int nframes = s->total_frames;

    printf("SFA: N=%u, L=%.1f, h=%.4f, frames=%d, cols=%u\n",
           N, L, h, nframes, s->n_columns);
    for (uint32_t c = 0; c < s->n_columns; c++)
        printf("  col %u: %s\n", c, s->columns[c].name);

    /* Allocate frame buffer */
    double *buf = (double *)malloc(s->frame_bytes);
    if (!buf) { fprintf(stderr, "malloc failed\n"); return 1; }

    /* Radial bins: accumulate over all frames */
    double dc_tx[NBINS], dc_ty[NBINS], dc_tz[NBINS];  /* signed theta sums */
    double rms_sum[NBINS];                               /* theta^2 sums */
    double count[NBINS];                                 /* voxel counts */
    memset(dc_tx, 0, sizeof(dc_tx));
    memset(dc_ty, 0, sizeof(dc_ty));
    memset(dc_tz, 0, sizeof(dc_tz));
    memset(rms_sum, 0, sizeof(rms_sum));
    memset(count, 0, sizeof(count));

    double dr = L / NBINS;  /* bin width: L/50 = 0.5 for L=25 */

    for (int frame = 0; frame < nframes; frame++) {
        if (sfa_read_frame(s, frame, buf) != 0) {
            fprintf(stderr, "Failed to read frame %d\n", frame);
            continue;
        }

        double *phi_x  = buf + 0 * N3;
        double *phi_y  = buf + 1 * N3;
        double *phi_z  = buf + 2 * N3;
        double *theta_x = buf + 3 * N3;
        double *theta_y = buf + 4 * N3;
        double *theta_z = buf + 5 * N3;

        /* Find centroid weighted by |P| = |phi_x * phi_y * phi_z| */
        double cx = 0, cy = 0, cz = 0, wtot = 0;
        for (uint32_t iz = 0; iz < N; iz++)
        for (uint32_t iy = 0; iy < N; iy++)
        for (uint32_t ix = 0; ix < N; ix++) {
            uint64_t idx = (uint64_t)iz * N * N + iy * N + ix;
            double P = fabs(phi_x[idx] * phi_y[idx] * phi_z[idx]);
            double x = -L + (ix + 0.5) * h;
            double y = -L + (iy + 0.5) * h;
            double z = -L + (iz + 0.5) * h;
            cx += P * x;
            cy += P * y;
            cz += P * z;
            wtot += P;
        }
        if (wtot > 0) { cx /= wtot; cy /= wtot; cz /= wtot; }

        if (frame % 50 == 0)
            printf("  frame %3d: centroid = (%.2f, %.2f, %.2f), |P|_tot = %.4e\n",
                   frame, cx, cy, cz, wtot);

        /* Bin into radial shells */
        for (uint32_t iz = 0; iz < N; iz++)
        for (uint32_t iy = 0; iy < N; iy++)
        for (uint32_t ix = 0; ix < N; ix++) {
            uint64_t idx = (uint64_t)iz * N * N + iy * N + ix;
            double x = -L + (ix + 0.5) * h - cx;
            double y = -L + (iy + 0.5) * h - cy;
            double z = -L + (iz + 0.5) * h - cz;
            double r = sqrt(x*x + y*y + z*z);

            int bin = (int)(r / dr);
            if (bin >= NBINS) continue;

            double tx = theta_x[idx];
            double ty = theta_y[idx];
            double tz = theta_z[idx];

            dc_tx[bin] += tx;
            dc_ty[bin] += ty;
            dc_tz[bin] += tz;
            rms_sum[bin] += tx*tx + ty*ty + tz*tz;
            count[bin] += 1.0;
        }
    }

    /* Compute time-averaged profiles and write TSV */
    const char *outpath = "dc_profile.tsv";
    FILE *fp = fopen(outpath, "w");
    if (!fp) { fprintf(stderr, "Cannot open %s\n", outpath); return 1; }

    fprintf(fp, "# Time-averaged signed theta radial profile\n");
    fprintf(fp, "# Source: %s (%d frames, N=%u, L=%.1f)\n", sfa_path, nframes, N, L);
    fprintf(fp, "# Columns:\n");
    fprintf(fp, "#   r_mid    - radial shell midpoint\n");
    fprintf(fp, "#   DC       - |mean(theta_vec)| = magnitude of time-averaged signed theta vector\n");
    fprintf(fp, "#   RMS      - sqrt(mean(theta^2)) = RMS of all theta components\n");
    fprintf(fp, "#   DC_RMS   - DC/RMS ratio (1.0 = perfectly coherent, 0.0 = pure noise)\n");
    fprintf(fp, "#   mean_tx  - time-averaged theta_x (signed)\n");
    fprintf(fp, "#   mean_ty  - time-averaged theta_y (signed)\n");
    fprintf(fp, "#   mean_tz  - time-averaged theta_z (signed)\n");
    fprintf(fp, "#   nvoxels  - total voxel-frame samples in this bin\n");
    fprintf(fp, "r_mid\tDC\tRMS\tDC_RMS\tmean_tx\tmean_ty\tmean_tz\tnvoxels\n");

    printf("\n--- Results ---\n");
    printf("%6s  %12s  %12s  %10s\n", "r_mid", "DC", "RMS", "DC/RMS");

    for (int b = 0; b < NBINS; b++) {
        double r_mid = (b + 0.5) * dr;
        if (count[b] < 1.0) {
            fprintf(fp, "%.4f\t0\t0\t0\t0\t0\t0\t0\n", r_mid);
            continue;
        }

        double mtx = dc_tx[b] / count[b];
        double mty = dc_ty[b] / count[b];
        double mtz = dc_tz[b] / count[b];
        double DC = sqrt(mtx*mtx + mty*mty + mtz*mtz);
        double RMS = sqrt(rms_sum[b] / count[b]);
        double ratio = (RMS > 0) ? DC / RMS : 0;

        fprintf(fp, "%.4f\t%.6e\t%.6e\t%.6f\t%.6e\t%.6e\t%.6e\t%.0f\n",
                r_mid, DC, RMS, ratio, mtx, mty, mtz, count[b]);

        printf("%6.2f  %12.6e  %12.6e  %10.6f\n", r_mid, DC, RMS, ratio);
    }

    fclose(fp);

    /* Print summary at specific radii */
    printf("\n--- Summary: DC/RMS at selected radii ---\n");
    double target_r[] = {2.0, 5.0, 10.0, 15.0, 20.0};
    int ntargets = 5;
    for (int t = 0; t < ntargets; t++) {
        int b = (int)(target_r[t] / dr);
        if (b >= NBINS || count[b] < 1) {
            printf("  r=%.0f: no data\n", target_r[t]);
            continue;
        }
        double mtx = dc_tx[b] / count[b];
        double mty = dc_ty[b] / count[b];
        double mtz = dc_tz[b] / count[b];
        double DC = sqrt(mtx*mtx + mty*mty + mtz*mtz);
        double RMS = sqrt(rms_sum[b] / count[b]);
        double ratio = (RMS > 0) ? DC / RMS : 0;
        printf("  r=%5.1f:  DC=%.4e  RMS=%.4e  DC/RMS=%.4f (%.2f%%)\n",
               target_r[t], DC, RMS, ratio, ratio * 100.0);
    }

    printf("\nInterpretation: DC/RMS << 1 means theta is AC-dominated (random sign).\n");
    printf("DC/RMS ~ 1 means theta has a coherent DC offset.\n");
    printf("The V34 claim was DC ~ 0.2%% of RMS.\n");

    printf("\nWrote: %s\n", outpath);

    free(buf);
    sfa_close(s);
    return 0;
}
