/*  extract_phases.c — Extract proton templates at 4 breathing phases + time-averaged
 *
 *  Reads specific frames from proton_formation.sfa corresponding to breathing
 *  cycle extrema and midpoints, extracts 64^3 sub-cube templates centered on
 *  the |P|-weighted centroid, and computes a time-averaged template across
 *  one complete breathing cycle (frames 264–321).
 *
 *  Outputs:
 *    proton_peak.sfa      — frame 321, most concentrated (P_int=164.6)
 *    proton_trough.sfa    — frame 264, most expanded (P_int=49.7)
 *    proton_rising.sfa    — frame 313, contracting through average (P_int=108.5)
 *    proton_falling.sfa   — frame 292, expanding through average (P_int=99.1)
 *    proton_averaged.sfa  — time-averaged over frames 264–321
 *    proton_phases.tsv    — radial profile comparison table
 *
 *  Build: gcc -O3 -o extract_phases extract_phases.c -lzstd -lm
 */

#define SFA_IMPLEMENTATION
#include "../../sfa/format/sfa.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* ---- Constants ---- */
#define TEMPLATE_N  64
#define NR_BINS     48
#define DR          0.25

/* Phase frames */
#define FRAME_PEAK     321
#define FRAME_TROUGH   264
#define FRAME_RISING   313
#define FRAME_FALLING  292
#define FRAME_AVG_START 264
#define FRAME_AVG_END   321

#define N_PHASES 4

static const int phase_frames[N_PHASES] = {FRAME_PEAK, FRAME_TROUGH, FRAME_RISING, FRAME_FALLING};
static const char *phase_names[N_PHASES] = {"peak", "trough", "rising", "falling"};
static const char *phase_files[N_PHASES] = {
    "proton_peak.sfa", "proton_trough.sfa", "proton_rising.sfa", "proton_falling.sfa"
};

/* ---- f16 to float conversion ---- */
static float f16_to_f32(uint16_t h) {
    int e = (h >> 10) & 0x1F;
    uint16_t m = h & 0x3FF;
    uint16_t s = h & 0x8000;
    if (e == 0) return 0.0f;
    uint32_t x = ((uint32_t)s << 16) | ((uint32_t)(e - 15 + 127) << 23) | ((uint32_t)m << 13);
    float fv;
    memcpy(&fv, &x, 4);
    return fv;
}

/* ---- Extract double arrays from SFA frame buffer ---- */
static void extract_fields(SFA *sfa, void *buf, long N3,
                           double **phi, double **theta,
                           double **phi_vel, double **theta_vel) {
    uint64_t off = 0;
    for (uint32_t c = 0; c < sfa->n_columns; c++) {
        int dtype = sfa->columns[c].dtype;
        int sem   = sfa->columns[c].semantic;
        int comp  = sfa->columns[c].component;
        int es    = sfa_dtype_size[dtype];
        uint8_t *src = (uint8_t *)buf + off;

        double *arr = (double *)malloc(N3 * sizeof(double));
        if (dtype == SFA_F64)
            memcpy(arr, src, N3 * 8);
        else if (dtype == SFA_F32)
            for (long i = 0; i < N3; i++) arr[i] = (double)((float *)src)[i];
        else if (dtype == SFA_F16)
            for (long i = 0; i < N3; i++) arr[i] = (double)f16_to_f32(((uint16_t *)src)[i]);
        else
            for (long i = 0; i < N3; i++) arr[i] = 0.0;

        if (sem == SFA_POSITION && comp < 3) phi[comp] = arr;
        else if (sem == SFA_ANGLE && comp < 3) theta[comp] = arr;
        else if (sem == SFA_VELOCITY && comp < 3) phi_vel[comp] = arr;
        else if (sem == SFA_VELOCITY && comp >= 3 && comp < 6) theta_vel[comp - 3] = arr;
        else free(arr);

        off += (uint64_t)N3 * es;
    }
}

static void free_fields(double **phi, double **theta,
                        double **phi_vel, double **theta_vel) {
    for (int a = 0; a < 3; a++) {
        if (phi[a])       { free(phi[a]);       phi[a] = NULL; }
        if (theta[a])     { free(theta[a]);     theta[a] = NULL; }
        if (phi_vel[a])   { free(phi_vel[a]);   phi_vel[a] = NULL; }
        if (theta_vel[a]) { free(theta_vel[a]); theta_vel[a] = NULL; }
    }
}

/* ---- Find |P|-weighted centroid ---- */
static void find_centroid(double **phi, int N, double L,
                          double *cx, double *cy, double *cz) {
    double h = 2.0 * L / N;
    double wsum = 0, wx = 0, wy = 0, wz = 0;
    for (int iz = 0; iz < N; iz++)
    for (int iy = 0; iy < N; iy++)
    for (int ix = 0; ix < N; ix++) {
        long idx = (long)iz * N * N + (long)iy * N + ix;
        double P = fabs(phi[0][idx] * phi[1][idx] * phi[2][idx]);
        if (P < 1e-12) continue;
        double x = -L + (ix + 0.5) * h;
        double y = -L + (iy + 0.5) * h;
        double z = -L + (iz + 0.5) * h;
        wx += P * x;  wy += P * y;  wz += P * z;
        wsum += P;
    }
    *cx = wx / wsum;
    *cy = wy / wsum;
    *cz = wz / wsum;
}

/* ---- Compute sub-cube start indices from centroid ---- */
static void subcube_start(double cx, double cy, double cz,
                          int N, double L, int TN,
                          int *x0, int *y0, int *z0) {
    double h = 2.0 * L / N;
    int ci = (int)((cx + L) / h);
    int cj = (int)((cy + L) / h);
    int ck = (int)((cz + L) / h);
    *x0 = ci - TN / 2; if (*x0 < 0) *x0 = 0;
    *y0 = cj - TN / 2; if (*y0 < 0) *y0 = 0;
    *z0 = ck - TN / 2; if (*z0 < 0) *z0 = 0;
    if (*x0 + TN > N) *x0 = N - TN;
    if (*y0 + TN > N) *y0 = N - TN;
    if (*z0 + TN > N) *z0 = N - TN;
}

/* ---- Radial profile data ---- */
typedef struct {
    double phi_rms_sum, theta_rms_sum, P_abs_sum;
    int count;
} ShellBin;

/* ---- Compute radial profile on a 64^3 sub-cube (in sub-cube coords) ---- */
static void compute_radial_profile_subcube(float **cols, int TN, double L_sub,
                                           ShellBin *bins) {
    double h = 2.0 * L_sub / TN;
    memset(bins, 0, NR_BINS * sizeof(ShellBin));
    int half = TN / 2;

    for (int iz = 0; iz < TN; iz++)
    for (int iy = 0; iy < TN; iy++)
    for (int ix = 0; ix < TN; ix++) {
        long idx = (long)iz * TN * TN + (long)iy * TN + ix;
        double x = (ix - half + 0.5) * h;
        double y = (iy - half + 0.5) * h;
        double z = (iz - half + 0.5) * h;
        double r = sqrt(x*x + y*y + z*z);

        int bin = (int)(r / DR);
        if (bin >= NR_BINS) continue;

        double p0 = cols[0][idx], p1 = cols[1][idx], p2 = cols[2][idx];
        double phi_rms = sqrt(p0*p0 + p1*p1 + p2*p2);
        double P_abs = fabs(p0 * p1 * p2);

        double t0 = cols[3][idx], t1 = cols[4][idx], t2 = cols[5][idx];
        double theta_rms = sqrt(t0*t0 + t1*t1 + t2*t2);

        bins[bin].phi_rms_sum += phi_rms;
        bins[bin].theta_rms_sum += theta_rms;
        bins[bin].P_abs_sum += P_abs;
        bins[bin].count++;
    }
}

/* ---- Extract sub-cube from full-grid fields into float arrays ---- */
static void extract_subcube(double **phi, double **theta,
                            double **phi_vel, double **theta_vel,
                            int N, int x0, int y0, int z0, int TN,
                            float **cols) {
    for (int tz = 0; tz < TN; tz++)
    for (int ty = 0; ty < TN; ty++)
    for (int tx = 0; tx < TN; tx++) {
        long tidx = (long)tz * TN * TN + (long)ty * TN + tx;
        int sx = x0 + tx, sy = y0 + ty, sz = z0 + tz;
        long sidx = (long)sz * N * N + (long)sy * N + sx;

        cols[0][tidx]  = (float)phi[0][sidx];
        cols[1][tidx]  = (float)phi[1][sidx];
        cols[2][tidx]  = (float)phi[2][sidx];
        cols[3][tidx]  = (float)theta[0][sidx];
        cols[4][tidx]  = (float)theta[1][sidx];
        cols[5][tidx]  = (float)theta[2][sidx];
        cols[6][tidx]  = phi_vel[0]   ? (float)phi_vel[0][sidx]   : 0.0f;
        cols[7][tidx]  = phi_vel[1]   ? (float)phi_vel[1][sidx]   : 0.0f;
        cols[8][tidx]  = phi_vel[2]   ? (float)phi_vel[2][sidx]   : 0.0f;
        cols[9][tidx]  = theta_vel[0] ? (float)theta_vel[0][sidx] : 0.0f;
        cols[10][tidx] = theta_vel[1] ? (float)theta_vel[1][sidx] : 0.0f;
        cols[11][tidx] = theta_vel[2] ? (float)theta_vel[2][sidx] : 0.0f;
    }
}

/* ---- Write a 64^3 sub-cube to SFA ---- */
static void write_template_sfa(const char *path, float **cols, int TN, double L_sub) {
    SFA *out = sfa_create(path, TN, TN, TN, L_sub, L_sub, L_sub, 1.0);
    sfa_add_column(out, "phi_x",    SFA_F32, SFA_POSITION, 0);
    sfa_add_column(out, "phi_y",    SFA_F32, SFA_POSITION, 1);
    sfa_add_column(out, "phi_z",    SFA_F32, SFA_POSITION, 2);
    sfa_add_column(out, "theta_x",  SFA_F32, SFA_ANGLE, 0);
    sfa_add_column(out, "theta_y",  SFA_F32, SFA_ANGLE, 1);
    sfa_add_column(out, "theta_z",  SFA_F32, SFA_ANGLE, 2);
    sfa_add_column(out, "phi_vx",   SFA_F32, SFA_VELOCITY, 0);
    sfa_add_column(out, "phi_vy",   SFA_F32, SFA_VELOCITY, 1);
    sfa_add_column(out, "phi_vz",   SFA_F32, SFA_VELOCITY, 2);
    sfa_add_column(out, "theta_vx", SFA_F32, SFA_VELOCITY, 3);
    sfa_add_column(out, "theta_vy", SFA_F32, SFA_VELOCITY, 4);
    sfa_add_column(out, "theta_vz", SFA_F32, SFA_VELOCITY, 5);
    sfa_finalize_header(out);

    void *ptrs[12];
    for (int c = 0; c < 12; c++) ptrs[c] = cols[c];
    sfa_write_frame(out, 0.0, ptrs);
    sfa_close(out);
}

/* ---- Main ---- */
int main(int argc, char **argv) {
    const char *sfa_path = "proton_formation.sfa";
    if (argc > 1) sfa_path = argv[1];

    /* Open SFA */
    SFA *sfa = sfa_open(sfa_path);
    if (!sfa) { fprintf(stderr, "Cannot open %s\n", sfa_path); return 1; }

    int N = sfa->Nx;
    double L = sfa->Lx;
    long N3 = (long)N * N * N;
    int TN = TEMPLATE_N;
    double h = 2.0 * L / N;
    double L_sub = TN * h / 2.0;

    printf("SFA: %s\n", sfa_path);
    printf("  Grid: %d^3, L=%.1f, h=%.4f, frames=%d, columns=%d\n",
           N, L, h, sfa->total_frames, sfa->n_columns);
    printf("  Template: %d^3, L_sub=%.3f\n", TN, L_sub);
    printf("  Columns:");
    for (uint32_t c = 0; c < sfa->n_columns; c++)
        printf(" %s(%s)", sfa->columns[c].name, sfa_dtype_name(sfa->columns[c].dtype));
    printf("\n\n");

    long TN3 = (long)TN * TN * TN;
    void *buf = malloc(sfa->frame_bytes);

    /* Allocate sub-cube float arrays (reused for each phase) */
    float *cols[12];
    for (int c = 0; c < 12; c++)
        cols[c] = (float *)malloc(TN3 * sizeof(float));

    /* Radial profiles for each phase + averaged */
    ShellBin phase_bins[N_PHASES + 1][NR_BINS];  /* +1 for averaged */

    /* ================================================================
       PART 1: Extract 4 phase templates
       ================================================================ */
    for (int p = 0; p < N_PHASES; p++) {
        int frame = phase_frames[p];
        printf("--- Phase: %s (frame %d) ---\n", phase_names[p], frame);

        if (sfa_read_frame(sfa, frame, buf) < 0) {
            fprintf(stderr, "Failed to read frame %d\n", frame);
            continue;
        }
        double t = sfa_frame_time(sfa, frame);
        printf("  t = %.4f\n", t);

        /* Extract fields */
        double *phi[3] = {NULL}, *theta[3] = {NULL};
        double *phi_vel[3] = {NULL}, *theta_vel[3] = {NULL};
        extract_fields(sfa, buf, N3, phi, theta, phi_vel, theta_vel);

        /* Centroid */
        double cx, cy, cz;
        find_centroid(phi, N, L, &cx, &cy, &cz);
        printf("  Centroid: (%.3f, %.3f, %.3f)\n", cx, cy, cz);

        /* Sub-cube extraction */
        int x0, y0, z0;
        subcube_start(cx, cy, cz, N, L, TN, &x0, &y0, &z0);
        printf("  Sub-cube start: (%d, %d, %d)\n", x0, y0, z0);

        extract_subcube(phi, theta, phi_vel, theta_vel, N, x0, y0, z0, TN, cols);

        /* Write SFA */
        write_template_sfa(phase_files[p], cols, TN, L_sub);
        printf("  Wrote %s\n", phase_files[p]);

        /* Radial profile */
        compute_radial_profile_subcube(cols, TN, L_sub, phase_bins[p]);

        /* Print summary */
        double P_total = 0;
        for (int b = 0; b < NR_BINS; b++)
            P_total += phase_bins[p][b].P_abs_sum;
        printf("  P_total (sum |P|) = %.4f\n\n", P_total * h * h * h);

        free_fields(phi, theta, phi_vel, theta_vel);
    }

    /* ================================================================
       PART 2: Time-averaged template (frames 264–321)
       ================================================================ */
    int n_avg_frames = FRAME_AVG_END - FRAME_AVG_START + 1;
    printf("=== Time-averaged template (frames %d-%d, %d frames) ===\n",
           FRAME_AVG_START, FRAME_AVG_END, n_avg_frames);

    /* Accumulator arrays (double precision) */
    double *acc[12];
    for (int c = 0; c < 12; c++) {
        acc[c] = (double *)calloc(TN3, sizeof(double));
    }

    int frames_read = 0;
    for (int frame = FRAME_AVG_START; frame <= FRAME_AVG_END; frame++) {
        if (sfa_read_frame(sfa, frame, buf) < 0) {
            fprintf(stderr, "Failed to read frame %d, skipping\n", frame);
            continue;
        }

        /* Extract fields */
        double *phi[3] = {NULL}, *theta[3] = {NULL};
        double *phi_vel[3] = {NULL}, *theta_vel[3] = {NULL};
        extract_fields(sfa, buf, N3, phi, theta, phi_vel, theta_vel);

        /* Find centroid for THIS frame (re-center each frame) */
        double cx, cy, cz;
        find_centroid(phi, N, L, &cx, &cy, &cz);

        /* Extract sub-cube centered on this frame's centroid */
        int x0, y0, z0;
        subcube_start(cx, cy, cz, N, L, TN, &x0, &y0, &z0);
        extract_subcube(phi, theta, phi_vel, theta_vel, N, x0, y0, z0, TN, cols);

        /* Accumulate */
        for (int c = 0; c < 12; c++)
            for (long i = 0; i < TN3; i++)
                acc[c][i] += (double)cols[c][i];

        frames_read++;
        free_fields(phi, theta, phi_vel, theta_vel);

        if (frame % 10 == 0)
            printf("  Frame %d: centroid=(%.3f, %.3f, %.3f)\n", frame, cx, cy, cz);
    }

    printf("  Averaged %d frames\n", frames_read);

    /* Divide by frame count to get average, store in float cols */
    double inv_n = 1.0 / frames_read;
    for (int c = 0; c < 12; c++)
        for (long i = 0; i < TN3; i++)
            cols[c][i] = (float)(acc[c][i] * inv_n);

    /* Write averaged template */
    write_template_sfa("proton_averaged.sfa", cols, TN, L_sub);
    printf("  Wrote proton_averaged.sfa\n");

    /* Radial profile of averaged template */
    compute_radial_profile_subcube(cols, TN, L_sub, phase_bins[N_PHASES]);

    /* Clean up accumulators */
    for (int c = 0; c < 12; c++) free(acc[c]);

    /* ================================================================
       PART 3: Comparison table — proton_phases.tsv
       ================================================================ */
    FILE *tsv = fopen("proton_phases.tsv", "w");
    fprintf(tsv, "r\t"
            "phi_rms_peak\tphi_rms_trough\tphi_rms_rising\tphi_rms_falling\tphi_rms_averaged\t"
            "theta_rms_peak\ttheta_rms_trough\ttheta_rms_rising\ttheta_rms_falling\ttheta_rms_averaged\t"
            "P_abs_peak\tP_abs_trough\tP_abs_rising\tP_abs_falling\tP_abs_averaged\n");

    printf("\n=== Radial Profile Comparison ===\n");
    printf("%6s  %10s %10s %10s %10s %10s  |  %10s %10s %10s %10s %10s  |  %10s %10s %10s %10s %10s\n",
           "r",
           "phi_peak", "phi_tro", "phi_ris", "phi_fal", "phi_avg",
           "tht_peak", "tht_tro", "tht_ris", "tht_fal", "tht_avg",
           "P_peak", "P_tro", "P_ris", "P_fal", "P_avg");
    printf("------  ---------- ---------- ---------- ---------- ----------  |  ---------- ---------- ---------- ---------- ----------  |  ---------- ---------- ---------- ---------- ----------\n");

    for (int b = 0; b < NR_BINS; b++) {
        double r = (b + 0.5) * DR;

        double phi_rms[5], theta_rms[5], P_abs[5];
        for (int p = 0; p < N_PHASES + 1; p++) {
            double n_inv = (phase_bins[p][b].count > 0) ? 1.0 / phase_bins[p][b].count : 0;
            phi_rms[p]   = phase_bins[p][b].phi_rms_sum * n_inv;
            theta_rms[p] = phase_bins[p][b].theta_rms_sum * n_inv;
            P_abs[p]     = phase_bins[p][b].P_abs_sum * n_inv;
        }

        fprintf(tsv, "%.4f", r);
        for (int p = 0; p < 5; p++) fprintf(tsv, "\t%.6e", phi_rms[p]);
        for (int p = 0; p < 5; p++) fprintf(tsv, "\t%.6e", theta_rms[p]);
        for (int p = 0; p < 5; p++) fprintf(tsv, "\t%.6e", P_abs[p]);
        fprintf(tsv, "\n");

        if (b < 28) {
            printf("%6.2f  %10.4e %10.4e %10.4e %10.4e %10.4e  |  %10.4e %10.4e %10.4e %10.4e %10.4e  |  %10.4e %10.4e %10.4e %10.4e %10.4e\n",
                   r,
                   phi_rms[0], phi_rms[1], phi_rms[2], phi_rms[3], phi_rms[4],
                   theta_rms[0], theta_rms[1], theta_rms[2], theta_rms[3], theta_rms[4],
                   P_abs[0], P_abs[1], P_abs[2], P_abs[3], P_abs[4]);
        }
    }

    fclose(tsv);
    printf("\nWrote proton_phases.tsv (%d bins, dr=%.2f)\n", NR_BINS, DR);

    /* Summary: peak-to-trough ratios at core */
    printf("\n=== Breathing Amplitude Summary ===\n");
    for (int b = 0; b < 8; b++) {
        double r = (b + 0.5) * DR;
        double n_peak = (phase_bins[0][b].count > 0) ? 1.0 / phase_bins[0][b].count : 0;
        double n_trough = (phase_bins[1][b].count > 0) ? 1.0 / phase_bins[1][b].count : 0;
        double n_avg = (phase_bins[4][b].count > 0) ? 1.0 / phase_bins[4][b].count : 0;

        double P_peak = phase_bins[0][b].P_abs_sum * n_peak;
        double P_trough = phase_bins[1][b].P_abs_sum * n_trough;
        double P_avg = phase_bins[4][b].P_abs_sum * n_avg;

        double phi_peak = phase_bins[0][b].phi_rms_sum * n_peak;
        double phi_trough = phase_bins[1][b].phi_rms_sum * n_trough;

        double ratio_P = (P_trough > 1e-15) ? P_peak / P_trough : 0;
        double ratio_phi = (phi_trough > 1e-15) ? phi_peak / phi_trough : 0;

        printf("  r=%.2f: phi_peak/trough=%.3f, P_peak/trough=%.3f, P_avg=%.4e\n",
               r, ratio_phi, ratio_P, P_avg);
    }

    /* Clean up */
    for (int c = 0; c < 12; c++) free(cols[c]);
    free(buf);
    sfa_close(sfa);

    printf("\nDone. Generated 5 templates + comparison table.\n");
    return 0;
}
