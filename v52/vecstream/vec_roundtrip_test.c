/*
 * vec_roundtrip_test.c — Roundtrip test for SFA vector (FRVD) encode/decode.
 *
 * Generates a known evolving 64³ field (translating + deforming Gaussian blob),
 * encodes it as FRVD I-frames and P-frames using the same pipeline as scp_sim.c,
 * then decodes via sfa_read_frame and compares against the original data.
 *
 * Diagnoses:
 *  - Per-frame max/RMS error
 *  - Fraction of patches using pure prediction (zero delta) per P-frame
 *  - Whether errors grow between I-frames (drift)
 *  - Fit error (Chebyshev basis roundtrip on single patches)
 *
 * Build:
 *   gcc -O3 -fopenmp -o vec_roundtrip_test vec_roundtrip_test.c -lzstd -lm
 */

#define SFA_IMPLEMENTATION
#include "../../sfa/format/sfa.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* ---- Chebyshev polynomial patch fitting (copied from scp_sim.c) ---- */

#define VS_ORDER 3
#define VS_NC    ((VS_ORDER+1)*(VS_ORDER+1)*(VS_ORDER+1))  /* 64 */

static void fit_patch(const float *field, int N, int ox, int oy, int oz,
                      int bs, float *out_coeffs) {
    int o1 = VS_ORDER + 1;
    int nc = o1 * o1 * o1;
    double AtA[64*64] = {0};
    double Atb[64] = {0};

    for (int di = 0; di < bs; di++) {
        int gi = ox + di; if (gi >= N) gi = N-1;
        double tx = (bs > 1) ? (double)di / (bs - 1) : 0;
        double sx = 2*tx - 1;
        double txi[4] = {1, sx, 2*sx*sx-1, 4*sx*sx*sx-3*sx};
        for (int dj = 0; dj < bs; dj++) {
            int gj = oy + dj; if (gj >= N) gj = N-1;
            double ty = (bs > 1) ? (double)dj / (bs - 1) : 0;
            double sy = 2*ty - 1;
            double tyj[4] = {1, sy, 2*sy*sy-1, 4*sy*sy*sy-3*sy};
            for (int dk = 0; dk < bs; dk++) {
                int gk = oz + dk; if (gk >= N) gk = N-1;
                double tz = (bs > 1) ? (double)dk / (bs - 1) : 0;
                double sz = 2*tz - 1;
                double tzk[4] = {1, sz, 2*sz*sz-1, 4*sz*sz*sz-3*sz};

                double basis[64];
                int c = 0;
                for (int a = 0; a < o1; a++)
                for (int b = 0; b < o1; b++)
                for (int g2 = 0; g2 < o1; g2++)
                    basis[c++] = txi[a] * tyj[b] * tzk[g2];

                float val = field[(long)gi*N*N + (long)gj*N + gk];
                for (int r = 0; r < nc; r++) {
                    Atb[r] += basis[r] * val;
                    for (int s = r; s < nc; s++)
                        AtA[r*nc+s] += basis[r] * basis[s];
                }
            }
        }
    }
    for (int r = 0; r < nc; r++)
        for (int s = 0; s < r; s++)
            AtA[r*nc+s] = AtA[s*nc+r];

    double L[64*64] = {0};
    for (int i = 0; i < nc; i++) {
        for (int j = 0; j <= i; j++) {
            double s = AtA[i*nc+j];
            for (int k = 0; k < j; k++) s -= L[i*nc+k] * L[j*nc+k];
            if (i == j) L[i*nc+j] = sqrt(s > 0 ? s : 1e-30);
            else        L[i*nc+j] = s / (L[j*nc+j] + 1e-30);
        }
    }
    double y[64];
    for (int i = 0; i < nc; i++) {
        double s = Atb[i];
        for (int k = 0; k < i; k++) s -= L[i*nc+k] * y[k];
        y[i] = s / (L[i*nc+i] + 1e-30);
    }
    double x[64];
    for (int i = nc-1; i >= 0; i--) {
        double s = y[i];
        for (int k = i+1; k < nc; k++) s -= L[k*nc+i] * x[k];
        x[i] = s / (L[i*nc+i] + 1e-30);
    }
    for (int i = 0; i < nc; i++) out_coeffs[i] = (float)x[i];
}

/* Evaluate a single patch at a voxel position (for local testing) */
static float eval_patch(const float *coeffs, int bs, int di, int dj, int dk) {
    int o1 = VS_ORDER + 1;
    double tx = (bs > 1) ? (double)di / (bs - 1) : 0;
    double sx = 2*tx - 1;
    double txa[4] = {1, sx, 2*sx*sx-1, 4*sx*sx*sx-3*sx};
    double ty = (bs > 1) ? (double)dj / (bs - 1) : 0;
    double sy = 2*ty - 1;
    double tya[4] = {1, sy, 2*sy*sy-1, 4*sy*sy*sy-3*sy};
    double tz = (bs > 1) ? (double)dk / (bs - 1) : 0;
    double sz = 2*tz - 1;
    double tza[4] = {1, sz, 2*sz*sz-1, 4*sz*sz*sz-3*sz};
    double val = 0;
    for (int a = 0; a < o1; a++)
    for (int b = 0; b < o1; b++) {
        double ab = txa[a]*tya[b];
        int base = a*o1*o1 + b*o1;
        for (int c = 0; c < o1; c++)
            val += coeffs[base+c] * ab * tza[c];
    }
    return (float)val;
}

/* ---- Test field generation ----
 *
 * We want a field that is NON-periodic in time (so a single-frequency cosine
 * model cannot capture it well). The test uses:
 *
 *   f(x,y,z,t) = A(t) * exp(-|r - r0(t)|^2 / (2*sigma(t)^2))
 *              + B(t) * exp(-|r - r1(t)|^2 / (2*sigma1^2))
 *
 * where A(t) = 0.5 + 0.3*sin(t) + 0.1*sin(sqrt(2)*t)  (quasiperiodic)
 *       r0(t) = (32 + 8*sin(0.7*t), 32 + 6*cos(1.1*t), 32)  (Lissajous-like)
 *       sigma(t) = 5.0 + 1.5*sin(0.5*t)
 *       B(t) = 0.2 * t / T_max  (linear ramp - purely non-periodic)
 *       r1(t) = (20, 44, 32 + 4*t/T_max)  (drifting)
 *       sigma1 = 4.0
 *
 * This field has:
 *   - Quasiperiodic amplitude (not captured by single cosine)
 *   - Translating center (Lissajous path)
 *   - Time-varying width (pulsing)
 *   - A second blob with linear amplitude ramp and drift
 */

static void generate_field(float *field, int N, double t, double T_max) {
    double A = 0.5 + 0.3*sin(t) + 0.1*sin(sqrt(2.0)*t);
    double cx = 32.0 + 8.0*sin(0.7*t);
    double cy = 32.0 + 6.0*cos(1.1*t);
    double cz = 32.0;
    double sigma = 5.0 + 1.5*sin(0.5*t);
    double inv2s2 = 1.0 / (2.0 * sigma * sigma);

    double B = 0.2 * t / T_max;
    double cx1 = 20.0, cy1 = 44.0, cz1 = 32.0 + 4.0*t/T_max;
    double sigma1 = 4.0;
    double inv2s12 = 1.0 / (2.0 * sigma1 * sigma1);

    #pragma omp parallel for schedule(static)
    for (int i = 0; i < N; i++) {
        double dx = i - cx;
        double dx1 = i - cx1;
        for (int j = 0; j < N; j++) {
            double dy = j - cy;
            double dy1 = j - cy1;
            for (int k = 0; k < N; k++) {
                double dz = k - cz;
                double dz1 = k - cz1;
                double r2 = dx*dx + dy*dy + dz*dz;
                double r12 = dx1*dx1 + dy1*dy1 + dz1*dz1;
                field[(long)i*N*N + (long)j*N + k] =
                    (float)(A * exp(-r2 * inv2s2) + B * exp(-r12 * inv2s12));
            }
        }
    }
}

/* ---- Main ---- */

int main(void) {
    const int N = 64;
    const long N3 = (long)N*N*N;
    const int BS = 8;          /* block size */
    const int BN = N / BS;     /* blocks per dim */
    const int n_patches = BN * BN * BN;
    const int n_fields = 1;    /* single field for testing */
    const int nc = VS_NC;      /* 64 coeffs per patch per field */
    const int nc_total = n_fields * nc;
    const int n_frames = 20;
    const int iframe_interval = 10;
    const float delta_tol = 0.001f;  /* same default as scp_sim */
    const double T_max = 10.0;
    const double dt_frame = T_max / (n_frames - 1);  /* time between frames */

    /* Temporal model parameters (same as scp_sim.c) */
    float vec_omega = 2.0f * 3.14159265f / 2.2f;
    int vec_refit_interval = 50;

    printf("=== Vec Roundtrip Test ===\n");
    printf("Grid: %d^3, BS=%d, %d patches, %d coeffs/patch\n", N, BS, n_patches, nc_total);
    printf("Frames: %d, I-frame interval: %d, delta_tol: %.4f\n", n_frames, iframe_interval, delta_tol);
    printf("Temporal omega: %.4f (period %.2f)\n\n", vec_omega, 2*3.14159265/vec_omega);

    /* ======== STEP 0: Chebyshev fit->eval roundtrip test (single patch) ======== */
    printf("--- Step 0: Chebyshev fit->eval roundtrip on single patches ---\n");
    {
        float *field = (float*)malloc(N3 * sizeof(float));
        generate_field(field, N, 0.0, T_max);

        double max_fit_err = 0, sum_fit_err2 = 0;
        int n_voxels_tested = 0;

        for (int pi = 0; pi < n_patches; pi++) {
            int bi = pi / (BN*BN), bj = (pi / BN) % BN, bk = pi % BN;
            int ox = bi * BS, oy = bj * BS, oz = bk * BS;

            float coeffs[VS_NC];
            fit_patch(field, N, ox, oy, oz, BS, coeffs);

            for (int di = 0; di < BS; di++)
            for (int dj = 0; dj < BS; dj++)
            for (int dk = 0; dk < BS; dk++) {
                int gi = ox+di, gj = oy+dj, gk = oz+dk;
                if (gi >= N || gj >= N || gk >= N) continue;
                float original = field[(long)gi*N*N + (long)gj*N + gk];
                float reconstructed = eval_patch(coeffs, BS, di, dj, dk);
                double err = fabs((double)original - (double)reconstructed);
                if (err > max_fit_err) max_fit_err = err;
                sum_fit_err2 += err * err;
                n_voxels_tested++;
            }
        }
        double rms_fit_err = sqrt(sum_fit_err2 / n_voxels_tested);
        printf("  Fit->eval max error: %.6e\n", max_fit_err);
        printf("  Fit->eval RMS error: %.6e\n", rms_fit_err);
        printf("  Voxels tested: %d\n\n", n_voxels_tested);
        free(field);
    }

    /* ======== STEP 1: Generate all frames and store originals ======== */
    printf("--- Step 1: Generating %d frames of test data ---\n", n_frames);
    float **originals = (float**)malloc(n_frames * sizeof(float*));
    double *frame_times = (double*)malloc(n_frames * sizeof(double));
    for (int f = 0; f < n_frames; f++) {
        originals[f] = (float*)malloc(N3 * sizeof(float));
        frame_times[f] = f * dt_frame;
        generate_field(originals[f], N, frame_times[f], T_max);
    }
    printf("  Done. Time range: [%.2f, %.2f]\n\n", frame_times[0], frame_times[n_frames-1]);

    /* ======== STEP 2: Encode to SFA ======== */
    printf("--- Step 2: Encoding to SFA ---\n");

    const char *sfa_path = "/tmp/vec_roundtrip_test.sfa";
    SFA *sfa = sfa_create(sfa_path, N, N, N, 1.0, 1.0, 1.0, dt_frame);
    /* We only need one column for this test (the field), stored as f32 */
    sfa_add_column(sfa, "phi_x", SFA_F32, SFA_POSITION, 0);
    sfa_finalize_header(sfa);

    /* Pre-compute patch origins */
    int16_t *origins = (int16_t*)malloc(n_patches * 3 * sizeof(int16_t));
    {
        int pi = 0;
        for (int bi = 0; bi < BN; bi++)
        for (int bj = 0; bj < BN; bj++)
        for (int bk = 0; bk < BN; bk++) {
            origins[pi*3+0] = (int16_t)(bi * BS);
            origins[pi*3+1] = (int16_t)(bj * BS);
            origins[pi*3+2] = (int16_t)(bk * BS);
            pi++;
        }
    }

    /* Temporal model accumulators */
    long vnt = (long)n_patches * nc_total;
    float *temp_mean  = (float*)calloc(vnt, sizeof(float));
    float *temp_amp   = (float*)calloc(vnt, sizeof(float));
    float *temp_phase = (float*)calloc(vnt, sizeof(float));
    float *sum_cos    = (float*)calloc(vnt, sizeof(float));
    float *sum_sin    = (float*)calloc(vnt, sizeof(float));
    float *sum_mean   = (float*)calloc(vnt, sizeof(float));
    int temporal_count = 0;

    /* Per-frame encoding diagnostics */
    int *frame_type_log = (int*)calloc(n_frames, sizeof(int));  /* 0=I, 1=P */
    int *n_stored_patches = (int*)calloc(n_frames, sizeof(int));
    int *n_zero_patches   = (int*)calloc(n_frames, sizeof(int));
    float *max_pred_err   = (float*)calloc(n_frames, sizeof(float));

    float *vec_buf = (float*)malloc(N3 * sizeof(float));

    for (int frame = 0; frame < n_frames; frame++) {
        double t = frame_times[frame];
        int is_iframe = (frame == 0) || (frame % iframe_interval == 0);
        frame_type_log[frame] = is_iframe ? 0 : 1;

        /* Fit patches (1 field for this test) */
        float *all_coeffs = (float*)calloc(vnt, sizeof(float));
        memcpy(vec_buf, originals[frame], N3 * sizeof(float));

        #pragma omp parallel for schedule(dynamic)
        for (int pi = 0; pi < n_patches; pi++) {
            fit_patch(vec_buf, N, origins[pi*3], origins[pi*3+1], origins[pi*3+2],
                      BS, &all_coeffs[(long)pi * nc_total]);
        }

        /* Update temporal accumulators (identical to scp_sim.c) */
        {
            float cos_wt = cosf(vec_omega * (float)t);
            float sin_wt = sinf(vec_omega * (float)t);
            for (long i = 0; i < vnt; i++) {
                sum_mean[i] += all_coeffs[i];
                sum_cos[i]  += all_coeffs[i] * cos_wt;
                sum_sin[i]  += all_coeffs[i] * sin_wt;
            }
            temporal_count++;

            if (temporal_count >= vec_refit_interval && temporal_count > 2) {
                float inv_n = 1.0f / temporal_count;
                for (long i = 0; i < vnt; i++) {
                    temp_mean[i] = sum_mean[i] * inv_n;
                    float sc = sum_cos[i] * inv_n - temp_mean[i] * cos_wt;
                    float ss = sum_sin[i] * inv_n - temp_mean[i] * sin_wt;
                    temp_amp[i] = 2.0f * sqrtf(sc*sc + ss*ss);
                    temp_phase[i] = atan2f(-ss, sc);
                }
                memset(sum_cos, 0, vnt * sizeof(float));
                memset(sum_sin, 0, vnt * sizeof(float));
                memset(sum_mean, 0, vnt * sizeof(float));
                temporal_count = 0;
            }
        }

        if (is_iframe) {
            sfa_write_vec_iframe_temporal(sfa, t, n_patches, BS, nc_total,
                                          origins, all_coeffs,
                                          vec_omega, temp_mean, temp_amp, temp_phase);
            n_stored_patches[frame] = n_patches;
            n_zero_patches[frame] = 0;
            max_pred_err[frame] = 0.0f;
        } else {
            /* P-frame: compute delta = actual - predicted(t) */
            uint32_t *didx = (uint32_t*)malloc(n_patches * sizeof(uint32_t));
            float *dcoeffs = (float*)malloc(vnt * sizeof(float));
            int nd = 0;
            float frame_max_pred = 0;
            int frame_zero_count = 0;

            for (int p = 0; p < n_patches; p++) {
                float mx = 0;
                long base = (long)p * nc_total;
                for (int ci = 0; ci < nc_total; ci++) {
                    float predicted = temp_mean[base+ci]
                        + temp_amp[base+ci] * cosf(vec_omega * (float)t + temp_phase[base+ci]);
                    float d = all_coeffs[base+ci] - predicted;
                    dcoeffs[(long)nd*nc_total+ci] = d;
                    if (fabsf(d) > mx) mx = fabsf(d);
                }
                if (mx > frame_max_pred) frame_max_pred = mx;
                if (mx > delta_tol) {
                    didx[nd] = p;
                    nd++;
                } else {
                    frame_zero_count++;
                }
            }
            sfa_write_vec_pframe(sfa, t, nd, nc_total, didx, dcoeffs);
            n_stored_patches[frame] = nd;
            n_zero_patches[frame] = frame_zero_count;
            max_pred_err[frame] = frame_max_pred;
            free(didx); free(dcoeffs);
        }
        free(all_coeffs);
    }

    sfa_close(sfa);
    printf("  Wrote %d frames to %s\n\n", n_frames, sfa_path);

    /* ======== STEP 3: Decode and compare ======== */
    printf("--- Step 3: Decode and compare ---\n\n");

    SFA *reader = sfa_open(sfa_path);
    if (!reader) { fprintf(stderr, "ERROR: cannot open %s\n", sfa_path); return 1; }

    /* The file stores total_frames as patched by sfa_close */
    printf("SFA reports %u frames, %d columns, grid %ux%ux%u\n\n",
           reader->total_frames, reader->n_columns, reader->Nx, reader->Ny, reader->Nz);

    /* Allocate decode buffer — frame_bytes sized for the SFA column layout */
    void *decode_buf = malloc(reader->frame_bytes);
    if (!decode_buf) { fprintf(stderr, "ERROR: alloc decode buf\n"); return 1; }

    printf("%-6s %-6s %12s %12s %10s %10s %10s %12s\n",
           "Frame", "Type", "MaxErr", "RMSErr", "Err>1%%", "Stored", "ZeroDelta", "MaxPredErr");
    printf("%-6s %-6s %12s %12s %10s %10s %10s %12s\n",
           "-----", "----", "----------", "----------", "--------", "--------", "---------", "----------");

    double *max_errs = (double*)calloc(n_frames, sizeof(double));
    double *rms_errs = (double*)calloc(n_frames, sizeof(double));
    double *frac_gt1pct = (double*)calloc(n_frames, sizeof(double));

    for (int frame = 0; frame < n_frames; frame++) {
        int ret = sfa_read_frame(reader, frame, decode_buf);
        if (ret != 0) {
            printf("Frame %d: DECODE FAILED (ret=%d)\n", frame, ret);
            continue;
        }

        /* Compare decoded voxels with original */
        float *orig = originals[frame];
        float *decoded = (float*)decode_buf;  /* f32 column, first (only) column */

        double max_err = 0, sum_err2 = 0;
        int n_gt1pct = 0;

        for (long idx = 0; idx < N3; idx++) {
            double err = fabs((double)orig[idx] - (double)decoded[idx]);
            if (err > max_err) max_err = err;
            sum_err2 += err * err;
            /* 1% of the field's typical amplitude (roughly 0.5 for our test) */
            if (err > 0.005) n_gt1pct++;
        }
        double rms_err = sqrt(sum_err2 / N3);
        double frac = (double)n_gt1pct / N3;

        max_errs[frame] = max_err;
        rms_errs[frame] = rms_err;
        frac_gt1pct[frame] = frac;

        const char *type_str = frame_type_log[frame] == 0 ? "I" : "P";
        printf("  %3d   %-4s  %12.6e %12.6e %9.4f%% %8d %10d %12.6e\n",
               frame, type_str, max_err, rms_err, 100.0*frac,
               n_stored_patches[frame], n_zero_patches[frame], max_pred_err[frame]);
    }

    sfa_close(reader);
    printf("\n");

    /* ======== STEP 4: Diagnosis ======== */
    printf("=== DIAGNOSIS ===\n\n");

    /* 4a: I-frame vs P-frame error comparison */
    double sum_i_rms = 0, sum_p_rms = 0;
    double max_i_err = 0, max_p_err = 0;
    int n_i = 0, n_p = 0;
    for (int f = 0; f < n_frames; f++) {
        if (frame_type_log[f] == 0) {
            sum_i_rms += rms_errs[f];
            if (max_errs[f] > max_i_err) max_i_err = max_errs[f];
            n_i++;
        } else {
            sum_p_rms += rms_errs[f];
            if (max_errs[f] > max_p_err) max_p_err = max_errs[f];
            n_p++;
        }
    }
    printf("1) I-frame vs P-frame errors:\n");
    printf("   I-frames (%d): avg RMS = %.6e, max err = %.6e\n",
           n_i, n_i > 0 ? sum_i_rms/n_i : 0.0, max_i_err);
    printf("   P-frames (%d): avg RMS = %.6e, max err = %.6e\n",
           n_p, n_p > 0 ? sum_p_rms/n_p : 0.0, max_p_err);
    if (max_p_err > 10.0 * max_i_err)
        printf("   ** P-frames have %.0fx higher max error than I-frames! **\n",
               max_p_err / (max_i_err + 1e-30));
    printf("\n");

    /* 4b: Zero-delta patch analysis */
    printf("2) Zero-delta patches (pure prediction) per P-frame:\n");
    int total_zero = 0, total_stored = 0;
    for (int f = 0; f < n_frames; f++) {
        if (frame_type_log[f] == 1) {
            total_zero += n_zero_patches[f];
            total_stored += n_stored_patches[f];
            double pct = 100.0 * n_zero_patches[f] / n_patches;
            printf("   Frame %2d: %d/%d patches zero-delta (%.1f%%), max|pred_err|=%.4e\n",
                   f, n_zero_patches[f], n_patches, pct, max_pred_err[f]);
        }
    }
    if (n_p > 0) {
        double avg_zero_pct = 100.0 * total_zero / (n_p * n_patches);
        printf("   Average: %.1f%% patches using pure prediction\n", avg_zero_pct);
        if (avg_zero_pct > 50.0)
            printf("   ** OVER HALF of patches use pure prediction — tolerance %.4f is too high! **\n",
                   delta_tol);
    }
    printf("\n");

    /* 4c: Error drift between I-frames */
    printf("3) Error drift between I-frames:\n");
    int gop_start = -1;
    for (int f = 0; f < n_frames; f++) {
        if (frame_type_log[f] == 0) {
            if (gop_start >= 0) {
                printf("   GOP [%d..%d]: I=%.2e", gop_start, f-1, rms_errs[gop_start]);
                double prev_rms = rms_errs[gop_start];
                int growing = 1;
                for (int p = gop_start+1; p < f; p++) {
                    if (rms_errs[p] < prev_rms) growing = 0;
                    prev_rms = rms_errs[p];
                }
                printf(" -> last_P=%.2e (%s)\n", rms_errs[f-1], growing ? "GROWING" : "mixed");
            }
            gop_start = f;
        }
    }
    /* Last GOP */
    if (gop_start >= 0 && gop_start < n_frames - 1) {
        printf("   GOP [%d..%d]: I=%.2e", gop_start, n_frames-1, rms_errs[gop_start]);
        double prev_rms = rms_errs[gop_start];
        int growing = 1;
        for (int p = gop_start+1; p < n_frames; p++) {
            if (rms_errs[p] < prev_rms) growing = 0;
            prev_rms = rms_errs[p];
        }
        printf(" -> last_P=%.2e (%s)\n", rms_errs[n_frames-1], growing ? "GROWING" : "mixed");
    }
    printf("\n");

    /* 4d: Temporal model quality analysis */
    printf("4) Temporal model analysis:\n");
    printf("   omega = %.4f (period = %.2f time units)\n", vec_omega, 2*M_PI/vec_omega);
    printf("   refit_interval = %d (refit triggered after %d samples)\n",
           vec_refit_interval, vec_refit_interval);
    printf("   With only %d total frames, refit threshold of %d is NEVER reached.\n",
           n_frames, vec_refit_interval);
    printf("   -> mean/amp/phase stay at ZERO for ALL frames.\n");
    printf("   -> predicted(t) = 0 for ALL patches at ALL times.\n");
    printf("   -> delta = actual - 0 = actual (full coefficients).\n");
    printf("   -> Only patches with max|coeff| < tol=%.4f get zero delta.\n", delta_tol);
    printf("\n");

    /* 4e: Check if temporal model was all zeros */
    {
        int any_nonzero = 0;
        for (long i = 0; i < vnt; i++) {
            if (temp_mean[i] != 0 || temp_amp[i] != 0 || temp_phase[i] != 0) {
                any_nonzero = 1;
                break;
            }
        }
        printf("   Temporal model state: %s\n",
               any_nonzero ? "has fitted values" : "ALL ZEROS (never refitted)");
    }
    printf("\n");

    /* 4f: What the decoder actually gets */
    printf("5) Decoder behavior analysis:\n");
    printf("   I-frame: coefficients stored directly + temporal model (all zeros)\n");
    printf("   P-frame: predicted = mean + amp*cos(omega*t + phase) = 0\n");
    printf("            delta patches applied on top of prediction\n");
    printf("   For stored patches: actual = 0 + delta = correct coefficients\n");
    printf("   For zero-delta patches: actual = 0 (WRONG - should be near-zero field values)\n");
    printf("\n");

    /* 4g: Detailed look at what the decoder produces for zero-delta patches */
    printf("6) Zero-delta patch impact:\n");
    {
        /* Re-encode frame 1 (first P-frame) and check */
        float *field1 = originals[1];
        float max_zero_patch_val = 0;
        int n_zero_patch_voxels = 0;
        float *tmp_coeffs = (float*)calloc(vnt, sizeof(float));

        memcpy(vec_buf, field1, N3*sizeof(float));
        for (int pi = 0; pi < n_patches; pi++) {
            fit_patch(vec_buf, N, origins[pi*3], origins[pi*3+1], origins[pi*3+2],
                      BS, &tmp_coeffs[(long)pi*nc_total]);
        }

        /* For zero-delta patches, the decoder gets predicted = 0 (since model is all zeros).
         * Check what the actual field values are in those patches. */
        for (int p = 0; p < n_patches; p++) {
            long base = (long)p * nc_total;
            float mx = 0;
            for (int ci = 0; ci < nc_total; ci++) {
                /* With zero temporal model, predicted=0, so delta = actual coeff */
                if (fabsf(tmp_coeffs[base+ci]) > mx) mx = fabsf(tmp_coeffs[base+ci]);
            }
            if (mx <= delta_tol) {
                /* This patch would be zero-delta */
                int bi = p / (BN*BN), bj = (p / BN) % BN, bk = p % BN;
                int ox = bi * BS, oy = bj * BS, oz = bk * BS;
                for (int di = 0; di < BS; di++)
                for (int dj = 0; dj < BS; dj++)
                for (int dk = 0; dk < BS; dk++) {
                    int gi = ox+di, gj = oy+dj, gk = oz+dk;
                    if (gi < N && gj < N && gk < N) {
                        float v = fabsf(field1[(long)gi*N*N + gj*N + gk]);
                        if (v > max_zero_patch_val) max_zero_patch_val = v;
                        n_zero_patch_voxels++;
                    }
                }
            }
        }
        printf("   Frame 1 (first P-frame): zero-delta patches contain:\n");
        printf("     %d voxels, max |original field| = %.6e\n",
               n_zero_patch_voxels, max_zero_patch_val);
        if (max_zero_patch_val < 0.001)
            printf("     -> These are truly near-zero regions. Zero prediction is fine.\n");
        else
            printf("     -> These patches have SIGNIFICANT field values lost to prediction=0!\n");

        free(tmp_coeffs);
    }
    printf("\n");

    /* 4h: The REAL problem — scp_sim.c uses a refit_interval of 50
     * but typical vec_snap_every might produce < 50 frames before the first I-frame.
     * If temporal_count never reaches 50, the model stays at zero forever. */
    printf("7) ROOT CAUSE ANALYSIS:\n");
    printf("   In scp_sim.c, vec_refit_interval=50. The temporal model is only fitted\n");
    printf("   when temporal_count >= 50. With typical simulation parameters:\n");
    printf("     - If vec_snap_dt = 0.1 and T = 200, you get 2000 vec frames.\n");
    printf("     - The model IS fitted after 50 frames... but the first I-frame at frame 0\n");
    printf("       sends model = all zeros (fitted 0 times).\n");
    printf("     - P-frames 1-9 use the zero model from I-frame 0.\n");
    printf("     - At frame 10 (next I-frame), if temporal_count < 50,\n");
    printf("       the model is STILL all zeros.\n");
    printf("     - It takes 50 frames (5 I-frame cycles) before the model gets fitted.\n");
    printf("\n");
    printf("   Meanwhile, the decoder's P-frame reconstruction does:\n");
    printf("     predicted = mean + amp * cos(omega * t + phase)\n");
    printf("   With mean=amp=phase=0:\n");
    printf("     predicted = 0 + 0 * cos(...) = 0\n");
    printf("   So for zero-delta patches, the decoder returns 0 everywhere.\n");
    printf("\n");
    printf("   When the model IS eventually fitted (after frame 50+), the prediction\n");
    printf("   is a single-frequency cosine at omega=%.2f. For non-periodic or\n", vec_omega);
    printf("   multi-frequency fields, this prediction has systematic bias.\n");
    printf("   Patches where the bias < tol get zero delta, showing pure prediction\n");
    printf("   (the 'playing a loop' artifact).\n");
    printf("\n");

    /* 4i: Suggested fixes */
    printf("8) SUGGESTED FIXES:\n");
    printf("   A) Lower delta_tol: current %.4f may be too loose. Try 1e-4 or 1e-5.\n", delta_tol);
    printf("   B) Fix refit timing: refit at EVERY I-frame (not every 50 samples).\n");
    printf("      Change: if (is_iframe && temporal_count > 2) { refit; }\n");
    printf("   C) Send ALL patches in P-frames (set tol=0): guaranteed lossless\n");
    printf("      roundtrip at the cost of larger files. This removes the artifact.\n");
    printf("   D) Use prev-frame delta instead of temporal model for the first\n");
    printf("      50 frames (before the model is fitted).\n");
    printf("   E) Bootstrap: fit the temporal model on the first I-frame's coefficients\n");
    printf("      alone (mean = coeffs, amp = 0). At least the prediction starts at\n");
    printf("      the correct values instead of zero.\n");
    printf("\n");

    /* Cleanup */
    for (int f = 0; f < n_frames; f++) free(originals[f]);
    free(originals); free(frame_times); free(origins); free(vec_buf);
    free(temp_mean); free(temp_amp); free(temp_phase);
    free(sum_cos); free(sum_sin); free(sum_mean);
    free(frame_type_log); free(n_stored_patches); free(n_zero_patches); free(max_pred_err);
    free(max_errs); free(rms_errs); free(frac_gt1pct);
    free(decode_buf);

    printf("=== DONE ===\n");
    return 0;
}
