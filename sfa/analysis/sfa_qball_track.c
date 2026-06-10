/*  sfa_qball_track.c -- Q-ball cluster tracker for complexified SFA runs
 *
 *  Detects clusters per frame via 26-connected flood fill on the
 *  phase-invariant density
 *
 *      rho2 = sum_a (u_a^2 + v_a^2),   a = x,y,z
 *
 *  where u = phi (real part) and v = phiim (imaginary part), and reports
 *  per cluster:
 *
 *    - mass      = int rho2 dV
 *    - Q         = int sum_a (u_a vdot_a - v_a udot_a) dV   (Noether charge)
 *    - centroid  (rho2-weighted, periodic min-image, wrapped to [-L, L))
 *    - rms size  (rho2-weighted rms radius about centroid)
 *    - rho2_peak (max rho2 in cluster)
 *
 *  Requires a 24-column complex SFA (complex_phi=1 kernel output) with
 *  columns phi_{x,y,z}, phiim_{x,y,z}, phi_v{x,y,z}, phiim_v{x,y,z}.
 *  Exits with a clear error on 12-column real-field files.
 *
 *  Threshold: voxels with rho2 >= threshold_frac * rho2_max(frame) are
 *  cluster material (default 0.30). Clusters smaller than
 *  MIN_CLUSTER_VOXELS voxels are ignored as specks.
 *
 *  Build:
 *    gcc -O3 -fopenmp -o sfa_qball_track sfa_qball_track.c -lzstd -lm
 *
 *  Usage:
 *    ./sfa_qball_track input.sfa [threshold_frac] [--tsv out.tsv]
 *
 *  Human-readable summary goes to stdout; optional machine-readable TSV
 *  (one row per cluster per frame) to the --tsv file.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define SFA_IMPLEMENTATION
#include "../format/sfa.h"

#define MAX_CLUSTERS       4096
#define MAX_REPORT         32     /* clusters printed per frame (by mass) */
#define MIN_CLUSTER_VOXELS 20     /* speck filter */
#define DEFAULT_THR_FRAC   0.30

/* ---- column lookup by name ---- */
static int find_column(const SFA *sfa, const char *name) {
    for (uint32_t c = 0; c < sfa->n_columns; c++)
        if (strncmp(sfa->columns[c].name, name, sizeof(sfa->columns[c].name)) == 0)
            return (int)c;
    return -1;
}

/* ---- extract one column from an uncompressed frame buffer as f32 ---- */
static float *extract_column_f32(const void *buf, const SFA *sfa, int col_idx) {
    uint64_t N3 = sfa->N_total, off = 0;
    for (int c = 0; c < col_idx; c++)
        off += N3 * sfa_dtype_size[sfa->columns[c].dtype];
    int dtype = sfa->columns[col_idx].dtype;
    const uint8_t *src = (const uint8_t *)buf + off;
    float *arr = (float *)malloc(N3 * sizeof(float));
    if (!arr) { fprintf(stderr, "error: malloc failed (column extract)\n"); exit(1); }
    if (dtype == SFA_F32) {
        memcpy(arr, src, N3 * sizeof(float));
    } else if (dtype == SFA_F64) {
        const double *d = (const double *)src;
        #pragma omp parallel for
        for (uint64_t i = 0; i < N3; i++) arr[i] = (float)d[i];
    } else if (dtype == SFA_F16) {
        const uint16_t *h = (const uint16_t *)src;
        #pragma omp parallel for
        for (uint64_t i = 0; i < N3; i++) arr[i] = sfa_f16_to_f32(h[i]);
    } else {
        fprintf(stderr, "error: unsupported dtype %d in column %d\n", dtype, col_idx);
        exit(1);
    }
    return arr;
}

static void usage(const char *prog) {
    fprintf(stderr,
        "usage: %s input.sfa [threshold_frac] [--tsv out.tsv]\n"
        "  threshold_frac : cluster threshold as fraction of per-frame rho2_max\n"
        "                   (default %.2f)\n"
        "  --tsv out.tsv  : write machine-readable per-cluster rows to out.tsv\n",
        prog, DEFAULT_THR_FRAC);
}

int main(int argc, char **argv) {
    const char *in_path = NULL, *tsv_path = NULL;
    double thr_frac = DEFAULT_THR_FRAC;

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "--tsv") == 0) {
            if (++i >= argc) { usage(argv[0]); return 1; }
            tsv_path = argv[i];
        } else if (!in_path) {
            in_path = argv[i];
        } else {
            thr_frac = atof(argv[i]);
            if (thr_frac <= 0.0 || thr_frac >= 1.0) {
                fprintf(stderr, "error: threshold_frac must be in (0,1), got %s\n", argv[i]);
                return 1;
            }
        }
    }
    if (!in_path) { usage(argv[0]); return 1; }

    SFA *sfa = sfa_open(in_path);
    if (!sfa) { fprintf(stderr, "error: cannot open %s\n", in_path); return 1; }

    /* ---- resolve complex field columns by name ---- */
    static const char *re_names[3]  = {"phi_x",     "phi_y",     "phi_z"};
    static const char *im_names[3]  = {"phiim_x",   "phiim_y",   "phiim_z"};
    static const char *rev_names[3] = {"phi_vx",    "phi_vy",    "phi_vz"};
    static const char *imv_names[3] = {"phiim_vx",  "phiim_vy",  "phiim_vz"};
    int re_cols[3], im_cols[3], rev_cols[3], imv_cols[3];
    int missing = 0;
    for (int a = 0; a < 3; a++) {
        re_cols[a]  = find_column(sfa, re_names[a]);
        im_cols[a]  = find_column(sfa, im_names[a]);
        rev_cols[a] = find_column(sfa, rev_names[a]);
        imv_cols[a] = find_column(sfa, imv_names[a]);
        if (re_cols[a] < 0 || im_cols[a] < 0 || rev_cols[a] < 0 || imv_cols[a] < 0)
            missing = 1;
    }
    if (missing) {
        fprintf(stderr,
            "error: %s lacks complex Q-ball columns (need phi_*, phiim_*, "
            "phi_v*, phiim_v*; found %u columns).\n"
            "This tool requires 24-column output from the complexified kernel "
            "(complex_phi=1). 12-column real-field SFAs are not supported -- "
            "use cluster_profile for those.\n",
            in_path, sfa->n_columns);
        sfa_close(sfa);
        return 1;
    }

    int N = (int)sfa->Nx;
    double L = sfa->Lx;
    double dx = 2.0 * L / (N - 1);
    double dV = dx * dx * dx;
    long N3 = (long)N * N * N;
    long NN = (long)N * N;
    uint32_t total = sfa->total_frames;

    printf("# sfa_qball_track: %s\n", in_path);
    printf("# N=%d L=%g dx=%g frames=%u threshold_frac=%g min_voxels=%d\n",
           N, L, dx, total, thr_frac, MIN_CLUSTER_VOXELS);

    FILE *tsv = NULL;
    if (tsv_path) {
        tsv = fopen(tsv_path, "w");
        if (!tsv) { fprintf(stderr, "error: cannot write %s\n", tsv_path); sfa_close(sfa); return 1; }
        fprintf(tsv, "frame\tt\tn_clusters\tcid\tnvox\tmass\tQ\tcx\tcy\tcz\trho2_peak\trms_size\n");
    }

    void *buf     = malloc(sfa->frame_bytes);
    float *rho2   = (float *)malloc(N3 * sizeof(float));
    float *qd     = (float *)malloc(N3 * sizeof(float));
    int   *labels = (int *)malloc(N3 * sizeof(int));
    long  *queue  = (long *)malloc(N3 * sizeof(long));
    if (!buf || !rho2 || !qd || !labels || !queue) {
        fprintf(stderr, "error: allocation failed (N=%d)\n", N);
        return 1;
    }

    /* per-cluster accumulators */
    static long   c_nvox[MAX_CLUSTERS], c_seed[MAX_CLUSTERS];
    static double c_mass[MAX_CLUSTERS], c_q[MAX_CLUSTERS], c_peak[MAX_CLUSTERS];
    static double c_sx[MAX_CLUSTERS], c_sy[MAX_CLUSTERS], c_sz[MAX_CLUSTERS], c_srr[MAX_CLUSTERS];

    for (uint32_t fi = 0; fi < total; fi++) {
        double t = sfa_frame_time(sfa, fi);
        if (sfa_read_frame(sfa, fi, buf) != 0) {
            fprintf(stderr, "error: read frame %u failed\n", fi);
            break;
        }

        /* build rho2 and Noether charge density */
        memset(rho2, 0, N3 * sizeof(float));
        memset(qd,   0, N3 * sizeof(float));
        for (int a = 0; a < 3; a++) {
            float *u  = extract_column_f32(buf, sfa, re_cols[a]);
            float *v  = extract_column_f32(buf, sfa, im_cols[a]);
            float *du = extract_column_f32(buf, sfa, rev_cols[a]);
            float *dv = extract_column_f32(buf, sfa, imv_cols[a]);
            #pragma omp parallel for
            for (long i = 0; i < N3; i++) {
                rho2[i] += u[i]*u[i] + v[i]*v[i];
                qd[i]   += u[i]*dv[i] - v[i]*du[i];
            }
            free(u); free(v); free(du); free(dv);
        }

        /* global stats */
        double rmax = 0, rho_tot = 0, q_tot = 0;
        #pragma omp parallel for reduction(max:rmax) reduction(+:rho_tot,q_tot)
        for (long i = 0; i < N3; i++) {
            if (rho2[i] > rmax) rmax = rho2[i];
            rho_tot += rho2[i];
            q_tot   += qd[i];
        }

        /* flood fill (26-connected, periodic) */
        double thr = thr_frac * rmax;
        memset(labels, 0, N3 * sizeof(int));
        int nc = 0;
        for (long idx = 0; idx < N3; idx++) {
            if (rho2[idx] < thr || labels[idx]) continue;
            if (nc >= MAX_CLUSTERS) {
                fprintf(stderr, "warning: frame %u exceeds %d clusters, truncating\n",
                        fi, MAX_CLUSTERS);
                break;
            }
            int lab = ++nc;
            int c = lab - 1;
            c_nvox[c] = 0; c_mass[c] = 0; c_q[c] = 0; c_peak[c] = 0;
            c_sx[c] = 0; c_sy[c] = 0; c_sz[c] = 0; c_srr[c] = 0;
            c_seed[c] = idx;
            long front = 0, back = 0;
            queue[back++] = idx;
            labels[idx] = lab;
            int si = (int)(idx / NN), sj = (int)((idx / N) % N), sk = (int)(idx % N);
            while (front < back) {
                long cur = queue[front++];
                int ci = (int)(cur / NN), cj = (int)((cur / N) % N), ck = (int)(cur % N);
                double w = rho2[cur];
                /* min-image offsets relative to seed voxel */
                int dii = ci - si; if (dii >  N/2) dii -= N; if (dii < -N/2) dii += N;
                int djj = cj - sj; if (djj >  N/2) djj -= N; if (djj < -N/2) djj += N;
                int dkk = ck - sk; if (dkk >  N/2) dkk -= N; if (dkk < -N/2) dkk += N;
                c_nvox[c]++;
                c_mass[c] += w;
                c_q[c]    += qd[cur];
                if (w > c_peak[c]) c_peak[c] = w;
                c_sx[c]  += w * dii;
                c_sy[c]  += w * djj;
                c_sz[c]  += w * dkk;
                c_srr[c] += w * (double)(dii*dii + djj*djj + dkk*dkk);
                for (int di = -1; di <= 1; di++)
                for (int dj = -1; dj <= 1; dj++)
                for (int dk = -1; dk <= 1; dk++) {
                    if (!di && !dj && !dk) continue;
                    int ni = (ci + di + N) % N, nj = (cj + dj + N) % N, nk = (ck + dk + N) % N;
                    long nidx = (long)ni * NN + (long)nj * N + nk;
                    if (!labels[nidx] && rho2[nidx] >= thr) {
                        labels[nidx] = lab;
                        queue[back++] = nidx;
                    }
                }
            }
        }

        /* sort by mass descending */
        int order[MAX_CLUSTERS];
        for (int c = 0; c < nc; c++) order[c] = c;
        for (int a = 0; a < nc; a++)
            for (int b = a + 1; b < nc; b++)
                if (c_mass[order[b]] > c_mass[order[a]]) {
                    int tmp = order[a]; order[a] = order[b]; order[b] = tmp;
                }

        /* count significant clusters */
        int nsig = 0;
        for (int c = 0; c < nc; c++)
            if (c_nvox[c] >= MIN_CLUSTER_VOXELS) nsig++;

        printf("\nframe %u  t=%.2f  rho2_max=%.5g  int_rho2=%.5g  Q_total=%.4f  clusters=%d\n",
               fi, t, rmax, rho_tot * dV, q_tot * dV, nsig);
        printf("  %-4s %-8s %-10s %-10s %-8s %-8s %-8s %-10s %-8s\n",
               "cid", "nvox", "mass", "Q", "cx", "cy", "cz", "rho2_peak", "rms");

        int printed = 0;
        for (int oi = 0; oi < nc; oi++) {
            int c = order[oi];
            if (c_nvox[c] < MIN_CLUSTER_VOXELS) continue;
            double m = c_mass[c];
            long seed = c_seed[c];
            int si = (int)(seed / NN), sj = (int)((seed / N) % N), sk = (int)(seed % N);
            double mx = c_sx[c] / m, my = c_sy[c] / m, mz = c_sz[c] / m;
            double cx = -L + (si + mx) * dx;
            double cy = -L + (sj + my) * dx;
            double cz = -L + (sk + mz) * dx;
            double span = 2.0 * L;
            while (cx >= L) cx -= span;
            while (cx < -L) cx += span;
            while (cy >= L) cy -= span;
            while (cy < -L) cy += span;
            while (cz >= L) cz -= span;
            while (cz < -L) cz += span;
            double var = c_srr[c] / m - (mx*mx + my*my + mz*mz);
            if (var < 0) var = 0;
            double rms = sqrt(var) * dx;
            int cid = printed + 1;
            if (printed < MAX_REPORT) {
                printf("  %-4d %-8ld %-10.4f %-10.4f %-8.3f %-8.3f %-8.3f %-10.5f %-8.3f\n",
                       cid, c_nvox[c], m * dV, c_q[c] * dV, cx, cy, cz, c_peak[c], rms);
            }
            if (tsv) {
                fprintf(tsv, "%u\t%.4f\t%d\t%d\t%ld\t%.6g\t%.6g\t%.4f\t%.4f\t%.4f\t%.6g\t%.4f\n",
                        fi, t, nsig, cid, c_nvox[c], m * dV, c_q[c] * dV,
                        cx, cy, cz, c_peak[c], rms);
            }
            printed++;
        }
        if (printed > MAX_REPORT)
            printf("  ... (%d more clusters, see TSV)\n", printed - MAX_REPORT);
        fflush(stdout);
        if (tsv) fflush(tsv);
    }

    if (tsv) fclose(tsv);
    free(buf); free(rho2); free(qd); free(labels); free(queue);
    sfa_close(sfa);
    return 0;
}
