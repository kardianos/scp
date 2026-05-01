/*  gen_temporal_seed.c — Extract temporal-mean seed from FRVD simulation output
 *
 *  Reads an SFA file containing FRVD I-frames with temporal models (produced by
 *  scp_sim with vec_snap_dt > 0). Finds the last I-frame, extracts the temporal
 *  mean polynomial coefficients, evaluates them to a voxel grid, and writes the
 *  result as a standard SFA seed file (6 fields, zero velocities).
 *
 *  The temporal model in each I-frame encodes: coeff(t) = mean + amp*cos(omega*t + phase)
 *  The "mean" represents the time-averaged equilibrium field. Extracting it produces
 *  a seed that should have minimal breathing oscillation when re-simulated.
 *
 *  Build: gcc -O3 -fopenmp -o gen_temporal_seed gen_temporal_seed.c -lzstd -lm
 *
 *  Usage: gen_temporal_seed input.sfa output_seed.sfa [-amp_report]
 *
 *  The -amp_report flag prints per-field amplitude statistics to help assess
 *  how much the field was oscillating (lower = more stable seed).
 */

#define SFA_IMPLEMENTATION
#include "../format/sfa.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifdef _OPENMP
#include <omp.h>
#endif

/* ---- Evaluate tricubic polynomial coefficients to a voxel grid ----
 *
 * The coefficient layout per patch is: 6 fields * 64 coefficients = 384 floats.
 * Each field's 64 coefficients are a 4x4x4 tricubic (order 3).
 * Polynomial basis: t = di/(bs-1), product of [1,t,t^2,t^3] in each dimension.
 */
static void eval_coeffs_to_voxels(
    const float *coeffs,        /* n_patches * nc_total floats */
    const int16_t *origins,     /* n_patches * 3 (ox, oy, oz) */
    uint32_t n_patches,
    int bs,                     /* block size (8) */
    uint16_t nc_total,          /* total coeffs per patch (384 = 6*64) */
    int N,                      /* grid dimension */
    int n_fields,               /* number of fields (6) */
    float **field_voxels)       /* output: n_fields arrays of N^3 floats */
{
    int pf_nc = nc_total / n_fields;  /* coeffs per field per patch (64) */
    int order = 3;
    if (pf_nc <= 27) order = 2;
    if (pf_nc <= 8)  order = 1;
    int o1 = order + 1;

    /* Zero output */
    long N3 = (long)N * N * N;
    for (int f = 0; f < n_fields; f++)
        memset(field_voxels[f], 0, N3 * sizeof(float));

    #ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic)
    #endif
    for (uint32_t pi = 0; pi < n_patches; pi++) {
        int ox = origins[pi*3], oy = origins[pi*3+1], oz = origins[pi*3+2];
        const float *pc = coeffs + (long)pi * nc_total;

        for (int di = 0; di < bs; di++) {
            int gi = ox + di;
            if (gi < 0 || gi >= N) continue;
            double tx = (bs > 1) ? (double)di / (bs - 1) : 0;
            double txa[4] = {1, tx, tx*tx, tx*tx*tx};
            for (int dj = 0; dj < bs; dj++) {
                int gj = oy + dj;
                if (gj < 0 || gj >= N) continue;
                double ty = (bs > 1) ? (double)dj / (bs - 1) : 0;
                double tya[4] = {1, ty, ty*ty, ty*ty*ty};
                for (int dk = 0; dk < bs; dk++) {
                    int gk = oz + dk;
                    if (gk < 0 || gk >= N) continue;
                    double tz = (bs > 1) ? (double)dk / (bs - 1) : 0;
                    double tza[4] = {1, tz, tz*tz, tz*tz*tz};
                    long idx = (long)gi * N * N + gj * N + gk;

                    for (int f = 0; f < n_fields; f++) {
                        const float *fc = pc + f * pf_nc;
                        double val = 0;
                        for (int a = 0; a < o1; a++)
                        for (int b = 0; b < o1; b++) {
                            double ab = txa[a] * tya[b];
                            int base = a * o1 * o1 + b * o1;
                            for (int c = 0; c < o1; c++)
                                val += fc[base + c] * ab * tza[c];
                        }
                        field_voxels[f][idx] = (float)val;
                    }
                }
            }
        }
    }
}

static void usage(const char *prog) {
    fprintf(stderr,
        "Usage: %s input.sfa output_seed.sfa [-amp_report]\n"
        "\n"
        "Extracts the temporal-mean field from the last FRVD I-frame in input.sfa\n"
        "and writes it as a standard SFA seed (6 fields, zero velocities).\n"
        "\n"
        "Options:\n"
        "  -amp_report   Print per-field amplitude statistics\n",
        prog);
}

int main(int argc, char **argv) {
    if (argc < 3) { usage(argv[0]); return 1; }

    const char *input_path = argv[1];
    const char *output_path = argv[2];
    int amp_report = 0;
    for (int i = 3; i < argc; i++) {
        if (!strcmp(argv[i], "-amp_report")) amp_report = 1;
    }

    /* Open input SFA */
    SFA *sfa = sfa_open(input_path);
    if (!sfa) {
        fprintf(stderr, "ERROR: cannot open '%s'\n", input_path);
        return 1;
    }
    printf("Input: %s\n", input_path);
    printf("  Grid: %ux%ux%u, L=(%.2f,%.2f,%.2f), %u columns, %u frames\n",
           sfa->Nx, sfa->Ny, sfa->Nz, sfa->Lx, sfa->Ly, sfa->Lz,
           sfa->n_columns, sfa->total_frames);

    /* Scan backwards to find the last I-frame with temporal model */
    int last_iframe = -1;
    for (int f = (int)sfa->total_frames - 1; f >= 0; f--) {
        uint32_t ft = sfa_frame_type(sfa, f);
        if (ft == SFA_FRAME_VEC_I) {
            last_iframe = f;
            break;
        }
    }
    if (last_iframe < 0) {
        fprintf(stderr, "ERROR: no FRVD I-frame found in '%s'\n", input_path);
        fprintf(stderr, "  (Run simulation with vec_snap_dt > 0 to produce FRVD frames)\n");
        sfa_close(sfa);
        return 1;
    }
    printf("  Last I-frame: index %d (time=%.2f)\n",
           last_iframe, sfa_frame_time(sfa, last_iframe));

    /* Read the I-frame to trigger temporal model parsing in sfa_read_frame.
     * We need the raw FRVD payload to get the temporal model directly, but
     * sfa_read_frame already parses it and stores it in the SFA struct. */
    void *voxbuf = malloc(sfa->frame_bytes);
    if (!voxbuf) {
        fprintf(stderr, "ERROR: malloc failed for frame buffer (%lu bytes)\n",
                (unsigned long)sfa->frame_bytes);
        sfa_close(sfa);
        return 1;
    }
    if (sfa_read_frame(sfa, last_iframe, voxbuf) != 0) {
        fprintf(stderr, "ERROR: failed to read I-frame %d\n", last_iframe);
        free(voxbuf);
        sfa_close(sfa);
        return 1;
    }
    free(voxbuf);  /* We don't need the voxel evaluation; we want the temporal model */

    /* Check that temporal model was loaded */
    if (!sfa->vec_temporal_valid) {
        fprintf(stderr, "ERROR: I-frame %d has no temporal model (flags & 0x01 not set)\n",
                last_iframe);
        fprintf(stderr, "  The simulation may not have run long enough to build a temporal model.\n");
        sfa_close(sfa);
        return 1;
    }

    uint32_t n_patches = sfa->vec_n_patches;
    uint16_t nc_total = sfa->vec_n_coeffs;
    float omega = sfa->vec_temporal_omega;
    long n_total = (long)n_patches * nc_total;

    printf("  Temporal model: %u patches, %u coeffs/patch, omega=%.4f\n",
           n_patches, nc_total, omega);

    /* Determine fields and block size */
    int n_fields = 6;
    int pf_nc = nc_total / n_fields;
    if (nc_total != (uint16_t)(n_fields * pf_nc)) {
        /* Try other field counts */
        if (nc_total % 64 == 0) { n_fields = nc_total / 64; pf_nc = 64; }
        else if (nc_total % 27 == 0) { n_fields = nc_total / 27; pf_nc = 27; }
        else if (nc_total % 8 == 0)  { n_fields = nc_total / 8;  pf_nc = 8; }
        else {
            fprintf(stderr, "ERROR: cannot determine field count from nc_total=%u\n", nc_total);
            sfa_close(sfa);
            return 1;
        }
    }
    printf("  Fields: %d, coeffs/field: %d\n", n_fields, pf_nc);

    /* Compute block size from n_patches and grid size */
    int N = sfa->Nx;
    int BN_cand;
    int bs = 0;
    for (int try_bs = 4; try_bs <= 16; try_bs++) {
        BN_cand = N / try_bs;
        if (BN_cand * try_bs == N && (uint32_t)(BN_cand * BN_cand * BN_cand) == n_patches) {
            bs = try_bs;
            break;
        }
    }
    if (bs == 0) {
        fprintf(stderr, "ERROR: cannot determine block size (N=%d, n_patches=%u)\n",
                N, n_patches);
        sfa_close(sfa);
        return 1;
    }
    printf("  Block size: %d (grid %d / %d blocks per dim)\n", bs, N, N/bs);

    /* Amplitude statistics */
    if (amp_report) {
        printf("\n  Amplitude statistics (temporal oscillation per field):\n");
        for (int f = 0; f < n_fields && f < 6; f++) {
            const char *names[] = {"phi_x", "phi_y", "phi_z", "theta_x", "theta_y", "theta_z"};
            double amp_sum = 0, amp_max = 0, mean_sum = 0;
            long count = 0;
            for (uint32_t pi = 0; pi < n_patches; pi++) {
                for (int ci = 0; ci < pf_nc; ci++) {
                    long idx = (long)pi * nc_total + f * pf_nc + ci;
                    double a = fabs(sfa->vec_temporal_amp[idx]);
                    double m = fabs(sfa->vec_temporal_mean[idx]);
                    amp_sum += a;
                    mean_sum += m;
                    if (a > amp_max) amp_max = a;
                    count++;
                }
            }
            printf("    %-10s  mean_avg=%.6f  amp_avg=%.6f  amp_max=%.6f  amp/mean=%.4f\n",
                   names[f], mean_sum / count, amp_sum / count, amp_max,
                   mean_sum > 0 ? amp_sum / mean_sum : 0);
        }
        printf("\n");
    }

    /* Reconstruct patch origins from grid geometry */
    int BN = N / bs;
    int16_t *origins = (int16_t*)malloc(n_patches * 3 * sizeof(int16_t));
    {
        int pi = 0;
        for (int bi = 0; bi < BN; bi++)
        for (int bj = 0; bj < BN; bj++)
        for (int bk = 0; bk < BN; bk++) {
            origins[pi*3+0] = (int16_t)(bi * bs);
            origins[pi*3+1] = (int16_t)(bj * bs);
            origins[pi*3+2] = (int16_t)(bk * bs);
            pi++;
        }
    }

    /* Evaluate temporal mean to voxels */
    long N3 = (long)N * N * N;
    float *field_voxels[6];
    for (int f = 0; f < n_fields; f++)
        field_voxels[f] = (float*)calloc(N3, sizeof(float));

    printf("Evaluating temporal mean to %d^3 voxel grid...\n", N);
    eval_coeffs_to_voxels(sfa->vec_temporal_mean, origins, n_patches,
                          bs, nc_total, N, n_fields, field_voxels);

    /* Report field statistics */
    for (int f = 0; f < n_fields && f < 6; f++) {
        const char *names[] = {"phi_x", "phi_y", "phi_z", "theta_x", "theta_y", "theta_z"};
        float fmin = 1e30f, fmax = -1e30f;
        double fsum = 0;
        for (long i = 0; i < N3; i++) {
            float v = field_voxels[f][i];
            if (v < fmin) fmin = v;
            if (v > fmax) fmax = v;
            fsum += v;
        }
        printf("  %-10s  min=%.6f  max=%.6f  mean=%.6e\n",
               names[f], fmin, fmax, fsum / N3);
    }

    /* Write output SFA seed: 6 fields + 6 zero velocities = 12 columns, f32 */
    printf("\nWriting seed: %s\n", output_path);
    SFA *out = sfa_create(output_path, N, N, N, sfa->Lx, sfa->Ly, sfa->Lz, sfa->dt);
    if (!out) {
        fprintf(stderr, "ERROR: cannot create output '%s'\n", output_path);
        for (int f = 0; f < n_fields; f++) free(field_voxels[f]);
        free(origins);
        sfa_close(sfa);
        return 1;
    }

    sfa_add_column(out, "phi_x",    SFA_F32, SFA_POSITION, 0);
    sfa_add_column(out, "phi_y",    SFA_F32, SFA_POSITION, 1);
    sfa_add_column(out, "phi_z",    SFA_F32, SFA_POSITION, 2);
    sfa_add_column(out, "theta_x",  SFA_F32, SFA_ANGLE,    0);
    sfa_add_column(out, "theta_y",  SFA_F32, SFA_ANGLE,    1);
    sfa_add_column(out, "theta_z",  SFA_F32, SFA_ANGLE,    2);
    sfa_add_column(out, "phi_vx",   SFA_F32, SFA_VELOCITY, 0);
    sfa_add_column(out, "phi_vy",   SFA_F32, SFA_VELOCITY, 1);
    sfa_add_column(out, "phi_vz",   SFA_F32, SFA_VELOCITY, 2);
    sfa_add_column(out, "theta_vx", SFA_F32, SFA_VELOCITY, 3);
    sfa_add_column(out, "theta_vy", SFA_F32, SFA_VELOCITY, 4);
    sfa_add_column(out, "theta_vz", SFA_F32, SFA_VELOCITY, 5);
    sfa_finalize_header(out);

    /* Build column data: 6 fields from temporal mean, 6 zero velocity arrays */
    float *zero_vel = (float*)calloc(N3, sizeof(float));
    void *cols[12];
    for (int f = 0; f < n_fields && f < 6; f++)
        cols[f] = field_voxels[f];
    /* If fewer than 6 fields, zero-fill remaining */
    for (int f = n_fields; f < 6; f++)
        cols[f] = zero_vel;
    /* Velocities: all zero (cold start from equilibrium) */
    for (int f = 0; f < 6; f++)
        cols[6 + f] = zero_vel;

    sfa_write_frame(out, 0.0, cols);
    sfa_close(out);

    printf("  Wrote 1 frame: %d^3 grid, L=%.2f, 12 columns (f32)\n", N, sfa->Lx);
    printf("  Velocities: zero (cold start from temporal mean)\n");

    /* Cleanup */
    free(zero_vel);
    for (int f = 0; f < n_fields; f++) free(field_voxels[f]);
    free(origins);
    sfa_close(sfa);

    printf("\nDone. Use as seed with: init=sfa, init_sfa=%s\n", output_path);
    return 0;
}
