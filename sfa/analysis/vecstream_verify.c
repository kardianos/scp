/*  vecstream_verify.c — Verify vecstream reconstruction accuracy
 *
 *  Usage: vecstream_verify input.vecstream [-field 0]
 *
 *  Reconstructs field at each K-frame time from I+P chain,
 *  compares against embedded K-frame voxels, reports error metrics.
 *
 *  Build: gcc -O3 -march=native -fopenmp -o vecstream_verify vecstream_verify.c -lzstd -lm
 */

#define VECSTREAM_IMPLEMENTATION
#include "../format/vecstream.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

int main(int argc, char **argv) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s input.vecstream [-field 0]\n", argv[0]);
        return 1;
    }

    const char *input_path = argv[1];
    int target_field = 0;

    for (int i = 2; i < argc; i++) {
        if (!strcmp(argv[i], "-field") && i+1 < argc) target_field = atoi(argv[++i]);
    }

    VecStream *vs = vecstream_open(input_path);
    if (!vs) { fprintf(stderr, "Cannot open %s\n", input_path); return 1; }

    uint32_t Nx = vs->Nx, Ny = vs->Ny, Nz = vs->Nz;
    uint64_t N3 = (uint64_t)Nx * Ny * Nz;
    uint32_t total_frames = vecstream_frame_count(vs);

    printf("vecstream_verify: %s\n", input_path);
    printf("  Grid: %u x %u x %u (%lu voxels)\n", Nx, Ny, Nz, (unsigned long)N3);
    printf("  Total index entries: %u\n", total_frames);
    printf("  Target field: %d\n", target_field);

    /* Count frame types */
    int n_i = 0, n_p = 0, n_k = 0;
    for (uint32_t i = 0; i < total_frames; i++) {
        if (vecstream_frame_field(vs, i) != (uint8_t)target_field) continue;
        switch (vecstream_frame_type(vs, i)) {
            case VS_FRAME_I: n_i++; break;
            case VS_FRAME_P: n_p++; break;
            case VS_FRAME_K: n_k++; break;
        }
    }
    printf("  Field %d: %d I-frames, %d P-frames, %d K-frames\n",
           target_field, n_i, n_p, n_k);

    if (n_k == 0) {
        printf("\n  No K-frames found for field %d. Cannot verify.\n", target_field);
        printf("  (Reconstruction accuracy can only be checked against K-frames.)\n");

        /* Still report I-frame self-consistency */
        printf("\n  Checking I-frame self-consistency...\n");
        float *recon = (float *)malloc(N3 * sizeof(float));
        for (uint32_t i = 0; i < total_frames; i++) {
            if (vecstream_frame_field(vs, i) != (uint8_t)target_field) continue;
            if (vecstream_frame_type(vs, i) != VS_FRAME_I) continue;

            double t = vecstream_frame_time(vs, i);
            if (vecstream_reconstruct(vs, i, (uint8_t)target_field, recon) != 0) {
                printf("  Frame %u (t=%.4f): RECONSTRUCT FAILED\n", i, t);
                continue;
            }

            /* Compute field statistics */
            double sum = 0, sum2 = 0, mn = 1e30, mx = -1e30;
            for (uint64_t v = 0; v < N3; v++) {
                double val = recon[v];
                sum += val;
                sum2 += val * val;
                if (val < mn) mn = val;
                if (val > mx) mx = val;
            }
            double mean = sum / N3;
            double rms = sqrt(sum2 / N3);
            printf("  I-frame %u (t=%.4f): mean=%.4e, rms=%.4e, range=[%.4e, %.4e]\n",
                   i, t, mean, rms, mn, mx);
        }
        free(recon);
        vecstream_close(vs);
        return 0;
    }

    /* For each K-frame: reconstruct from I+P chain, compare with K-frame voxels */
    float *k_voxels = (float *)malloc(N3 * sizeof(float));
    float *recon_voxels = (float *)malloc(N3 * sizeof(float));

    printf("\n%8s  %10s  %12s  %12s  %12s  %12s  %s\n",
           "KFrame", "Time", "MaxErr", "RMSErr", "MeanErr", "FieldRMS", "Status");

    int n_verified = 0, n_passed = 0;
    double worst_max_err = 0;
    double worst_rms_err = 0;

    /* Find the I-frame or preceding reconstructable frame for each K-frame */
    for (uint32_t i = 0; i < total_frames; i++) {
        if (vecstream_frame_field(vs, i) != (uint8_t)target_field) continue;
        if (vecstream_frame_type(vs, i) != VS_FRAME_K) continue;

        double t = vecstream_frame_time(vs, i);

        /* Read K-frame raw voxels */
        if (vecstream_read_kframe(vs, i, k_voxels) != 0) {
            printf("%8u  %10.4f  %12s  %12s  %12s  %12s  KFRAME_READ_FAIL\n",
                   i, t, "-", "-", "-", "-");
            continue;
        }

        /* Find the corresponding I-frame or P-frame at the same time.
         * The K-frame is paired with an I-frame at the same time. We need
         * to reconstruct using the I+P chain that precedes the K-frame. */
        /* Find the I-frame at same time for this field */
        int recon_idx = -1;
        for (int j = (int)i - 1; j >= 0; j--) {
            if (vecstream_frame_field(vs, j) != (uint8_t)target_field) continue;
            if (fabs(vecstream_frame_time(vs, j) - t) < 1e-12) {
                if (vecstream_frame_type(vs, j) == VS_FRAME_I ||
                    vecstream_frame_type(vs, j) == VS_FRAME_P) {
                    recon_idx = j;
                    break;
                }
            }
        }

        if (recon_idx < 0) {
            printf("%8u  %10.4f  %12s  %12s  %12s  %12s  NO_MATCHING_IFRAME\n",
                   i, t, "-", "-", "-", "-");
            continue;
        }

        /* Reconstruct from I+P chain */
        if (vecstream_reconstruct(vs, (uint32_t)recon_idx,
                                  (uint8_t)target_field, recon_voxels) != 0) {
            printf("%8u  %10.4f  %12s  %12s  %12s  %12s  RECONSTRUCT_FAIL\n",
                   i, t, "-", "-", "-", "-");
            continue;
        }

        /* Compute error metrics */
        double max_err = 0, sum_err = 0, sum_err2 = 0;
        double sum_val2 = 0;
        for (uint64_t v = 0; v < N3; v++) {
            double err = fabs((double)recon_voxels[v] - (double)k_voxels[v]);
            if (err > max_err) max_err = err;
            sum_err += err;
            sum_err2 += err * err;
            sum_val2 += (double)k_voxels[v] * k_voxels[v];
        }
        double mean_err = sum_err / N3;
        double rms_err = sqrt(sum_err2 / N3);
        double field_rms = sqrt(sum_val2 / N3);

        /* Relative error: rms_err / field_rms */
        double rel_err = (field_rms > 0) ? rms_err / field_rms : rms_err;

        if (max_err > worst_max_err) worst_max_err = max_err;
        if (rms_err > worst_rms_err) worst_rms_err = rms_err;

        const char *status = (rel_err < 0.01) ? "PASS" : "WARN";
        if (rel_err >= 0.05) status = "FAIL";

        printf("%8u  %10.4f  %12.4e  %12.4e  %12.4e  %12.4e  %s (%.2f%%)\n",
               i, t, max_err, rms_err, mean_err, field_rms, status, rel_err * 100);

        n_verified++;
        if (rel_err < 0.01) n_passed++;
    }

    printf("\n=== Verification Summary ===\n");
    printf("  K-frames verified: %d\n", n_verified);
    printf("  Passed (<1%% relative error): %d\n", n_passed);
    printf("  Worst max error: %.4e\n", worst_max_err);
    printf("  Worst RMS error: %.4e\n", worst_rms_err);

    if (n_verified > 0 && n_passed == n_verified)
        printf("  Result: ALL PASSED\n");
    else if (n_verified == 0)
        printf("  Result: NO K-FRAMES TO VERIFY\n");
    else
        printf("  Result: %d/%d PASSED\n", n_passed, n_verified);

    free(k_voxels);
    free(recon_voxels);
    vecstream_close(vs);
    return 0;
}
