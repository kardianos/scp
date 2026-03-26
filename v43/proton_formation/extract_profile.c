/*  extract_profile.c — Extract radial profile + 3D template from converged proton
 *
 *  Reads the last frame of proton_formation.sfa, finds the proton centroid,
 *  computes radial shell averages, and extracts a 64^3 sub-cube template.
 *
 *  Also reads frame 0 to compare initial vs converged profiles, and reads
 *  the diagnostics TSV to characterize formation dynamics.
 *
 *  Build: gcc -O3 -o extract_profile extract_profile.c -lzstd -lm
 */

#define SFA_IMPLEMENTATION
#include "../../sfa/format/sfa.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

/* ---- f16 → float conversion ---- */
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

/* ---- Extract double arrays from a raw SFA frame buffer ---- */
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

static void free_fields(double **phi, double **theta, double **phi_vel, double **theta_vel) {
    for (int a = 0; a < 3; a++) {
        if (phi[a]) free(phi[a]);
        if (theta[a]) free(theta[a]);
        if (phi_vel[a]) free(phi_vel[a]);
        if (theta_vel[a]) free(theta_vel[a]);
    }
}

/* ---- Radial profile data ---- */
#define NR_BINS 48
#define DR 0.25
#define TEMPLATE_N 64

typedef struct {
    double phi_rms_sum, theta_rms_sum, P_abs_sum;
    double phi_v_rms_sum, theta_v_rms_sum, v_radial_sum;
    double phi_max, theta_max;
    int count;
} ShellBin;

/* ---- Find centroid weighted by |P| ---- */
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
        wx += P * x;
        wy += P * y;
        wz += P * z;
        wsum += P;
    }
    *cx = wx / wsum;
    *cy = wy / wsum;
    *cz = wz / wsum;
}

/* ---- Compute radial profile ---- */
static void compute_radial_profile(double **phi, double **theta,
                                   double **phi_vel, double **theta_vel,
                                   int N, double L,
                                   double cx, double cy, double cz,
                                   ShellBin *bins) {
    double h = 2.0 * L / N;
    memset(bins, 0, NR_BINS * sizeof(ShellBin));

    for (int iz = 0; iz < N; iz++)
    for (int iy = 0; iy < N; iy++)
    for (int ix = 0; ix < N; ix++) {
        long idx = (long)iz * N * N + (long)iy * N + ix;
        double x = -L + (ix + 0.5) * h - cx;
        double y = -L + (iy + 0.5) * h - cy;
        double z = -L + (iz + 0.5) * h - cz;
        double r = sqrt(x * x + y * y + z * z);

        int bin = (int)(r / DR);
        if (bin >= NR_BINS) continue;

        double p0 = phi[0][idx], p1 = phi[1][idx], p2 = phi[2][idx];
        double phi_rms = sqrt(p0 * p0 + p1 * p1 + p2 * p2);
        double P_abs = fabs(p0 * p1 * p2);

        double t0 = theta[0][idx], t1 = theta[1][idx], t2 = theta[2][idx];
        double theta_rms = sqrt(t0 * t0 + t1 * t1 + t2 * t2);

        double pv_rms = 0, tv_rms = 0, v_rad = 0;
        if (phi_vel[0]) {
            double vx = phi_vel[0][idx], vy = phi_vel[1][idx], vz = phi_vel[2][idx];
            pv_rms = sqrt(vx * vx + vy * vy + vz * vz);
            if (r > 1e-8)
                v_rad = (vx * x + vy * y + vz * z) / r;
        }
        if (theta_vel[0]) {
            double tvx = theta_vel[0][idx], tvy = theta_vel[1][idx], tvz = theta_vel[2][idx];
            tv_rms = sqrt(tvx * tvx + tvy * tvy + tvz * tvz);
        }

        bins[bin].phi_rms_sum += phi_rms;
        bins[bin].theta_rms_sum += theta_rms;
        bins[bin].P_abs_sum += P_abs;
        bins[bin].phi_v_rms_sum += pv_rms;
        bins[bin].theta_v_rms_sum += tv_rms;
        bins[bin].v_radial_sum += v_rad;
        if (phi_rms > bins[bin].phi_max) bins[bin].phi_max = phi_rms;
        if (theta_rms > bins[bin].theta_max) bins[bin].theta_max = theta_rms;
        bins[bin].count++;
    }
}

/* ---- Extract 3D sub-cube template ---- */
static void extract_template(SFA *sfa, void *buf,
                             int N, double L,
                             double cx, double cy, double cz) {
    long N3 = (long)N * N * N;
    double h = 2.0 * L / N;
    int TN = TEMPLATE_N;

    /* Center voxel indices */
    int ci = (int)((cx + L) / h);
    int cj = (int)((cy + L) / h);
    int ck = (int)((cz + L) / h);

    /* Sub-cube start (clamped) */
    int x0 = ci - TN / 2; if (x0 < 0) x0 = 0;
    int y0 = cj - TN / 2; if (y0 < 0) y0 = 0;
    int z0 = ck - TN / 2; if (z0 < 0) z0 = 0;
    if (x0 + TN > N) x0 = N - TN;
    if (y0 + TN > N) y0 = N - TN;
    if (z0 + TN > N) z0 = N - TN;

    printf("Template sub-cube: start=(%d,%d,%d), size=%d^3\n", x0, y0, z0, TN);

    /* Extract fields to double arrays */
    double *phi[3] = {NULL}, *theta[3] = {NULL};
    double *phi_vel[3] = {NULL}, *theta_vel[3] = {NULL};
    extract_fields(sfa, buf, N3, phi, theta, phi_vel, theta_vel);

    /* Template domain: L_template = TN/2 * h */
    double L_tmpl = TN * h / 2.0;
    double dt_tmpl = 1.0; /* placeholder */

    /* Create template SFA with f32 columns */
    SFA *out = sfa_create("proton_template.sfa", TN, TN, TN,
                          L_tmpl, L_tmpl, L_tmpl, dt_tmpl);
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

    long TN3 = (long)TN * TN * TN;
    float *col_data[12];
    for (int c = 0; c < 12; c++)
        col_data[c] = (float *)malloc(TN3 * sizeof(float));

    /* Copy sub-cube data */
    for (int tz = 0; tz < TN; tz++)
    for (int ty = 0; ty < TN; ty++)
    for (int tx = 0; tx < TN; tx++) {
        long tidx = (long)tz * TN * TN + (long)ty * TN + tx;
        int sx = x0 + tx, sy = y0 + ty, sz = z0 + tz;
        long sidx = (long)sz * N * N + (long)sy * N + sx;

        col_data[0][tidx] = (float)phi[0][sidx];
        col_data[1][tidx] = (float)phi[1][sidx];
        col_data[2][tidx] = (float)phi[2][sidx];
        col_data[3][tidx] = (float)theta[0][sidx];
        col_data[4][tidx] = (float)theta[1][sidx];
        col_data[5][tidx] = (float)theta[2][sidx];
        col_data[6][tidx]  = phi_vel[0]  ? (float)phi_vel[0][sidx]  : 0.0f;
        col_data[7][tidx]  = phi_vel[1]  ? (float)phi_vel[1][sidx]  : 0.0f;
        col_data[8][tidx]  = phi_vel[2]  ? (float)phi_vel[2][sidx]  : 0.0f;
        col_data[9][tidx]  = theta_vel[0] ? (float)theta_vel[0][sidx] : 0.0f;
        col_data[10][tidx] = theta_vel[1] ? (float)theta_vel[1][sidx] : 0.0f;
        col_data[11][tidx] = theta_vel[2] ? (float)theta_vel[2][sidx] : 0.0f;
    }

    void *ptrs[12];
    for (int c = 0; c < 12; c++) ptrs[c] = col_data[c];
    sfa_write_frame(out, 0.0, ptrs);
    sfa_close(out);

    for (int c = 0; c < 12; c++) free(col_data[c]);
    free_fields(phi, theta, phi_vel, theta_vel);

    printf("Wrote proton_template.sfa: %d^3 grid, L=%.3f, 12 columns (f32)\n", TN, L_tmpl);
}

/* ---- Diagnostics analysis ---- */
typedef struct {
    double time;
    double E_total, E_pot, P_int, theta_rms;
    double E_phi_kin, E_theta_kin, E_grad, E_mass;
} DiagRow;

static int read_diag(const char *path, DiagRow **rows_out) {
    FILE *fp = fopen(path, "r");
    if (!fp) { fprintf(stderr, "Cannot open %s\n", path); return 0; }

    int cap = 1024, n = 0;
    DiagRow *rows = (DiagRow *)malloc(cap * sizeof(DiagRow));
    char line[4096];
    fgets(line, sizeof(line), fp); /* skip header */

    while (fgets(line, sizeof(line), fp)) {
        DiagRow *r = &rows[n];
        double E_tgrad, E_tmass, E_coupling;
        if (sscanf(line, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%*f\t%*f\t%lf\t%lf",
                   &r->time, &r->E_phi_kin, &r->E_theta_kin, &r->E_grad, &r->E_mass,
                   &r->E_pot, &E_tgrad, &E_tmass, &E_coupling, &r->E_total,
                   &r->P_int, &r->theta_rms) >= 12) {
            n++;
            if (n >= cap) { cap *= 2; rows = realloc(rows, cap * sizeof(DiagRow)); }
        }
    }
    fclose(fp);
    *rows_out = rows;
    return n;
}

static void analyze_formation(DiagRow *rows, int n, FILE *md) {
    if (n < 10) {
        fprintf(md, "Too few diagnostic rows (%d) for analysis.\n", n);
        return;
    }

    /* Find convergence time: when |dP_int/dt| < threshold for 10 consecutive steps */
    double threshold = 0.5; /* P_int change per time unit */
    int converge_idx = -1;
    for (int i = 10; i < n - 10; i++) {
        int stable = 1;
        for (int j = 0; j < 10; j++) {
            double dp = fabs(rows[i + j + 1].P_int - rows[i + j].P_int);
            double dt = rows[i + j + 1].time - rows[i + j].time;
            if (dt > 0 && dp / dt > threshold) { stable = 0; break; }
        }
        if (stable) { converge_idx = i; break; }
    }

    /* Find breathing period from E_pot oscillations in converged region */
    double breathing_period = 0;
    int start = (converge_idx > 0) ? converge_idx : n / 2;
    /* Find zero-crossings of (E_pot - E_pot_mean) */
    double E_pot_mean = 0;
    for (int i = start; i < n; i++) E_pot_mean += rows[i].E_pot;
    E_pot_mean /= (n - start);

    int n_crossings = 0;
    double first_cross = 0, last_cross = 0;
    for (int i = start + 1; i < n; i++) {
        double a = rows[i - 1].E_pot - E_pot_mean;
        double b = rows[i].E_pot - E_pot_mean;
        if (a * b < 0) {
            double t_cross = rows[i - 1].time + (rows[i].time - rows[i - 1].time) *
                             fabs(a) / (fabs(a) + fabs(b));
            if (n_crossings == 0) first_cross = t_cross;
            last_cross = t_cross;
            n_crossings++;
        }
    }
    if (n_crossings >= 2)
        breathing_period = 2.0 * (last_cross - first_cross) / (n_crossings - 1);

    /* Equilibrium values (average over last 20% of frames) */
    int eq_start = (int)(0.8 * n);
    double eq_P_int = 0, eq_E_total = 0, eq_E_pot = 0;
    double eq_E_kin = 0, eq_theta_rms = 0;
    int eq_n = n - eq_start;
    for (int i = eq_start; i < n; i++) {
        eq_P_int += rows[i].P_int;
        eq_E_total += rows[i].E_total;
        eq_E_pot += rows[i].E_pot;
        eq_E_kin += rows[i].E_phi_kin;
        eq_theta_rms += rows[i].theta_rms;
    }
    eq_P_int /= eq_n;
    eq_E_total /= eq_n;
    eq_E_pot /= eq_n;
    eq_E_kin /= eq_n;
    eq_theta_rms /= eq_n;

    /* P_int range in equilibrium */
    double P_min = 1e30, P_max = -1e30;
    double E_pot_min = 1e30, E_pot_max = -1e30;
    for (int i = eq_start; i < n; i++) {
        if (rows[i].P_int < P_min) P_min = rows[i].P_int;
        if (rows[i].P_int > P_max) P_max = rows[i].P_int;
        if (rows[i].E_pot < E_pot_min) E_pot_min = rows[i].E_pot;
        if (rows[i].E_pot > E_pot_max) E_pot_max = rows[i].E_pot;
    }

    fprintf(md, "# Proton Formation Analysis\n\n");
    fprintf(md, "## Formation Dynamics\n\n");
    fprintf(md, "| Quantity | Value |\n");
    fprintf(md, "|---|---|\n");
    fprintf(md, "| Total simulation time | %.2f |\n", rows[n - 1].time);
    fprintf(md, "| Number of diagnostics rows | %d |\n", n);
    if (converge_idx >= 0)
        fprintf(md, "| Convergence time (P_int stable) | %.2f |\n", rows[converge_idx].time);
    else
        fprintf(md, "| Convergence time | not reached |\n");
    fprintf(md, "| Breathing period (from E_pot) | %.3f |\n", breathing_period);
    fprintf(md, "| Equilibrium P_int | %.2f (range: %.2f -- %.2f) |\n", eq_P_int, P_min, P_max);
    fprintf(md, "| Equilibrium E_total | %.2f |\n", eq_E_total);
    fprintf(md, "| Equilibrium E_pot | %.2f (range: %.2f -- %.2f) |\n", eq_E_pot, E_pot_min, E_pot_max);
    fprintf(md, "| Equilibrium E_phi_kin | %.2f |\n", eq_E_kin);
    fprintf(md, "| Equilibrium theta_rms | %.6f |\n", eq_theta_rms);
    fprintf(md, "\n");

    /* Early vs late comparison */
    fprintf(md, "## Time Evolution Summary\n\n");
    fprintf(md, "| t | E_total | E_pot | P_int | theta_rms |\n");
    fprintf(md, "|---|---|---|---|---|\n");
    int steps[] = {0, n/10, n/4, n/2, 3*n/4, n-1};
    for (int s = 0; s < 6; s++) {
        int i = steps[s];
        fprintf(md, "| %.2f | %.2f | %.2f | %.2f | %.6f |\n",
                rows[i].time, rows[i].E_total, rows[i].E_pot,
                rows[i].P_int, rows[i].theta_rms);
    }
    fprintf(md, "\n");
}

/* ---- Main ---- */
int main(int argc, char **argv) {
    const char *sfa_path = "proton_formation.sfa";
    const char *diag_path = "proton_formation_diag.tsv";

    if (argc > 1) sfa_path = argv[1];
    if (argc > 2) diag_path = argv[2];

    /* ---- Open SFA ---- */
    SFA *sfa = sfa_open(sfa_path);
    if (!sfa) { fprintf(stderr, "Cannot open %s\n", sfa_path); return 1; }

    int N = sfa->Nx;
    double L = sfa->Lx;
    long N3 = (long)N * N * N;
    int last_frame = sfa->total_frames - 1;

    printf("SFA: %s\n", sfa_path);
    printf("  Grid: %d^3, L=%.1f, frames=%d, columns=%d\n",
           N, L, sfa->total_frames, sfa->n_columns);
    printf("  Columns:");
    for (uint32_t c = 0; c < sfa->n_columns; c++)
        printf(" %s(%s)", sfa->columns[c].name, sfa_dtype_name(sfa->columns[c].dtype));
    printf("\n");
    printf("  Last frame: %d (t=%.2f)\n", last_frame, sfa_frame_time(sfa, last_frame));
    printf("  Frame bytes: %lu\n", (unsigned long)sfa->frame_bytes);

    /* ---- Read frame 0 (initial) ---- */
    printf("\nReading frame 0 (initial)...\n");
    void *buf0 = malloc(sfa->frame_bytes);
    if (sfa_read_frame(sfa, 0, buf0) < 0) {
        fprintf(stderr, "Failed to read frame 0\n");
        free(buf0); sfa_close(sfa); return 1;
    }

    double *phi0[3] = {NULL}, *theta0[3] = {NULL};
    double *phi_vel0[3] = {NULL}, *theta_vel0[3] = {NULL};
    extract_fields(sfa, buf0, N3, phi0, theta0, phi_vel0, theta_vel0);
    free(buf0);

    /* Centroid of initial state */
    double cx0, cy0, cz0;
    find_centroid(phi0, N, L, &cx0, &cy0, &cz0);
    printf("  Initial centroid: (%.3f, %.3f, %.3f)\n", cx0, cy0, cz0);

    ShellBin bins0[NR_BINS];
    compute_radial_profile(phi0, theta0, phi_vel0, theta_vel0, N, L, cx0, cy0, cz0, bins0);

    /* ---- Read last frame (converged) ---- */
    printf("Reading frame %d (converged)...\n", last_frame);
    void *buf = malloc(sfa->frame_bytes);
    if (sfa_read_frame(sfa, last_frame, buf) < 0) {
        fprintf(stderr, "Failed to read frame %d\n", last_frame);
        free(buf); sfa_close(sfa); return 1;
    }

    double *phi[3] = {NULL}, *theta[3] = {NULL};
    double *phi_vel[3] = {NULL}, *theta_vel[3] = {NULL};
    extract_fields(sfa, buf, N3, phi, theta, phi_vel, theta_vel);

    /* Centroid of converged state */
    double cx, cy, cz;
    find_centroid(phi, N, L, &cx, &cy, &cz);
    printf("  Converged centroid: (%.3f, %.3f, %.3f)\n", cx, cy, cz);
    printf("  Centroid drift: (%.3f, %.3f, %.3f)\n", cx - cx0, cy - cy0, cz - cz0);

    /* ---- Radial profile (converged) ---- */
    ShellBin bins[NR_BINS];
    compute_radial_profile(phi, theta, phi_vel, theta_vel, N, L, cx, cy, cz, bins);

    /* Write profile TSV */
    FILE *tsv = fopen("proton_profile.tsv", "w");
    fprintf(tsv, "r\tphi_rms\ttheta_rms\tP_abs\tphi_v_rms\ttheta_v_rms\tv_radial\tphi_max\ttheta_max\tcount"
                 "\tphi_rms_init\ttheta_rms_init\tP_abs_init\n");
    for (int b = 0; b < NR_BINS; b++) {
        double r = (b + 0.5) * DR;
        ShellBin *s = &bins[b];
        ShellBin *s0 = &bins0[b];
        double n_inv = (s->count > 0) ? 1.0 / s->count : 0;
        double n0_inv = (s0->count > 0) ? 1.0 / s0->count : 0;
        fprintf(tsv, "%.4f\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%d\t%.6e\t%.6e\t%.6e\n",
                r,
                s->phi_rms_sum * n_inv,
                s->theta_rms_sum * n_inv,
                s->P_abs_sum * n_inv,
                s->phi_v_rms_sum * n_inv,
                s->theta_v_rms_sum * n_inv,
                s->v_radial_sum * n_inv,
                s->phi_max, s->theta_max,
                s->count,
                s0->phi_rms_sum * n0_inv,
                s0->theta_rms_sum * n0_inv,
                s0->P_abs_sum * n0_inv);
    }
    fclose(tsv);
    printf("\nWrote proton_profile.tsv (%d radial bins, dr=%.2f)\n", NR_BINS, DR);

    /* Print key profile info */
    printf("\n--- Converged Radial Profile ---\n");
    printf("%6s  %10s  %10s  %10s  %10s  %10s  %6s\n",
           "r", "phi_rms", "theta_rms", "|P|", "phi_v", "v_rad", "count");
    for (int b = 0; b < NR_BINS && b < 32; b++) {
        double r = (b + 0.5) * DR;
        ShellBin *s = &bins[b];
        double n_inv = (s->count > 0) ? 1.0 / s->count : 0;
        printf("%6.2f  %10.4e  %10.4e  %10.4e  %10.4e  %10.4e  %6d\n",
               r, s->phi_rms_sum * n_inv, s->theta_rms_sum * n_inv,
               s->P_abs_sum * n_inv, s->phi_v_rms_sum * n_inv,
               s->v_radial_sum * n_inv, s->count);
    }

    /* ---- Extract 3D template ---- */
    printf("\nExtracting %d^3 template...\n", TEMPLATE_N);
    /* Re-read last frame for template extraction (buf already has raw data) */
    extract_template(sfa, buf, N, L, cx, cy, cz);

    free(buf);
    free_fields(phi, theta, phi_vel, theta_vel);
    free_fields(phi0, theta0, phi_vel0, theta_vel0);
    sfa_close(sfa);

    /* ---- Diagnostics analysis ---- */
    printf("\nAnalyzing diagnostics...\n");
    DiagRow *diag_rows = NULL;
    int n_diag = read_diag(diag_path, &diag_rows);
    printf("  Read %d diagnostic rows\n", n_diag);

    FILE *md = fopen("formation_analysis.md", "w");

    if (n_diag > 0) {
        analyze_formation(diag_rows, n_diag, md);
        free(diag_rows);
    }

    /* Add profile comparison to the report */
    fprintf(md, "## Radial Profile Comparison (Initial vs Converged)\n\n");
    fprintf(md, "| r | phi_rms (init) | phi_rms (conv) | ratio | P_abs (init) | P_abs (conv) | ratio |\n");
    fprintf(md, "|---|---|---|---|---|---|---|\n");
    for (int b = 0; b < NR_BINS && b < 24; b++) {
        double r = (b + 0.5) * DR;
        double n_inv  = (bins[b].count > 0)  ? 1.0 / bins[b].count  : 0;
        double n0_inv = (bins0[b].count > 0) ? 1.0 / bins0[b].count : 0;
        double phi_conv = bins[b].phi_rms_sum * n_inv;
        double phi_init = bins0[b].phi_rms_sum * n0_inv;
        double P_conv = bins[b].P_abs_sum * n_inv;
        double P_init = bins0[b].P_abs_sum * n0_inv;
        double phi_ratio = (phi_init > 1e-12) ? phi_conv / phi_init : 0;
        double P_ratio = (P_init > 1e-12) ? P_conv / P_init : 0;
        fprintf(md, "| %.2f | %.4e | %.4e | %.3f | %.4e | %.4e | %.3f |\n",
                r, phi_init, phi_conv, phi_ratio, P_init, P_conv, P_ratio);
    }
    fprintf(md, "\n");

    /* Template info */
    fprintf(md, "## Proton Template\n\n");
    fprintf(md, "- File: `proton_template.sfa`\n");
    fprintf(md, "- Grid: %d^3, 12 columns (f32)\n", TEMPLATE_N);
    fprintf(md, "- Centroid at grid center\n");
    fprintf(md, "- Converged centroid position: (%.3f, %.3f, %.3f)\n", cx, cy, cz);
    fprintf(md, "- Initial centroid position: (%.3f, %.3f, %.3f)\n", cx0, cy0, cz0);
    fprintf(md, "\nTo use: read this SFA, place the 64^3 block at desired (cx,cy,cz) in a larger grid,\n");
    fprintf(md, "add A_bg=0.1 background outside the proton, and write as a new seed SFA.\n");

    fclose(md);
    printf("\nWrote formation_analysis.md\n");
    printf("Done.\n");
    return 0;
}
