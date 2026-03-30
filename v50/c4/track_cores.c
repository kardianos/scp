/* track_cores.c — Track two P-maxima across multiple SFA frames
 *
 * For each frame, finds the two strongest |P| peaks (P = phi0*phi1*phi2),
 * separated by at least N/6 grid cells, and reports positions, P values,
 * separation, gap minimum P, and gap maximum theta_rms.
 *
 * Build: gcc -O3 -o track_cores track_cores.c -lzstd -lm
 * Usage: ./track_cores input.sfa first_frame last_frame [step]
 */

#define SFA_IMPLEMENTATION
#include "../../sfa/format/sfa.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

/* 3D index to flat index (row-major: i*N*N + j*N + k) */
static inline long idx3(int i, int j, int k, int N) {
    return (long)i * N * N + (long)j * N + k;
}

int main(int argc, char **argv) {
    if (argc < 4) {
        fprintf(stderr, "Usage: %s input.sfa first_frame last_frame [step]\n", argv[0]);
        return 1;
    }

    const char *path = argv[1];
    int first = atoi(argv[2]);
    int last  = atoi(argv[3]);
    int step  = (argc > 4) ? atoi(argv[4]) : 1;

    SFA *sfa = sfa_open(path);
    if (!sfa) { fprintf(stderr, "Cannot open %s\n", path); return 1; }

    int N = sfa->Nx;
    double L = sfa->Lx;
    double dx = 2.0 * L / (N - 1);
    long N3 = (long)N * N * N;
    int NN = N * N;
    int min_sep = N / 6;  /* minimum separation between peaks in grid cells */

    fprintf(stderr, "File: %s  N=%d L=%.1f dx=%.4f  frames %d..%d step %d\n",
            path, N, L, dx, first, last, step);
    fprintf(stderr, "Columns: %d  frame_bytes: %lu  min_sep: %d cells\n",
            sfa->n_columns, (unsigned long)sfa->frame_bytes, min_sep);

    /* Allocate frame buffer */
    size_t frame_bytes = 0;
    for (uint32_t c = 0; c < sfa->n_columns; c++)
        frame_bytes += N3 * sfa_dtype_size[sfa->columns[c].dtype];
    void *buf = malloc(frame_bytes);
    if (!buf) { fprintf(stderr, "Cannot allocate %zu bytes\n", frame_bytes); return 1; }

    /* Header */
    printf("frame\ttime\tpk1_x\tpk1_y\tpk1_z\tpk1_P\tpk2_x\tpk2_y\tpk2_z\tpk2_P\tseparation\tgap_min_P\tgap_max_theta\n");

    for (int fr = first; fr <= last; fr += step) {
        if (sfa_read_frame(sfa, fr, buf) < 0) {
            fprintf(stderr, "Cannot read frame %d, skipping\n", fr);
            continue;
        }
        double t = sfa_frame_time(sfa, fr);

        /* Pointer to each column (all float32, 6 columns: phi0,phi1,phi2,theta0,theta1,theta2) */
        float *phi[3], *theta[3];
        size_t col_bytes = N3 * sizeof(float);
        for (int a = 0; a < 3; a++) {
            phi[a]   = (float*)((char*)buf + a * col_bytes);
            theta[a] = (float*)((char*)buf + (3 + a) * col_bytes);
        }

        /* Find peak 1: global maximum of |P| */
        typedef struct { int i, j, k; double P; } Peak;
        Peak pk1 = {0, 0, 0, 0.0};
        for (long idx = 0; idx < N3; idx++) {
            double P = fabs((double)phi[0][idx] * phi[1][idx] * phi[2][idx]);
            if (P > pk1.P) {
                pk1.P = P;
                pk1.i = idx / NN;
                pk1.j = (idx / N) % N;
                pk1.k = idx % N;
            }
        }

        /* Find peak 2: maximum |P| at least min_sep cells from peak 1 */
        Peak pk2 = {0, 0, 0, 0.0};
        for (long idx = 0; idx < N3; idx++) {
            int i = idx / NN, j = (idx / N) % N, k = idx % N;
            int di = abs(i - pk1.i), dj = abs(j - pk1.j), dk = abs(k - pk1.k);
            if (di > N/2) di = N - di;
            if (dj > N/2) dj = N - dj;
            if (dk > N/2) dk = N - dk;
            if (di + dj + dk < min_sep) continue;
            double P = fabs((double)phi[0][idx] * phi[1][idx] * phi[2][idx]);
            if (P > pk2.P) {
                pk2.P = P;
                pk2.i = i; pk2.j = j; pk2.k = k;
            }
        }

        /* World coordinates */
        double x1 = -L + pk1.i * dx, y1 = -L + pk1.j * dx, z1 = -L + pk1.k * dx;
        double x2 = -L + pk2.i * dx, y2 = -L + pk2.j * dx, z2 = -L + pk2.k * dx;

        /* Separation in world coordinates */
        double ddx = (pk1.i - pk2.i) * dx;
        double ddy = (pk1.j - pk2.j) * dx;
        double ddz = (pk1.k - pk2.k) * dx;
        /* Handle periodic wrapping */
        if (fabs(ddx) > L) ddx = (ddx > 0) ? ddx - 2*L : ddx + 2*L;
        if (fabs(ddy) > L) ddy = (ddy > 0) ? ddy - 2*L : ddy + 2*L;
        if (fabs(ddz) > L) ddz = (ddz > 0) ? ddz - 2*L : ddz + 2*L;
        double sep = sqrt(ddx*ddx + ddy*ddy + ddz*ddz);

        /* Trace the line between pk1 and pk2, sample points along it
         * to find the minimum |P| in the gap and maximum theta_rms in the gap.
         * Use 200 sample points along the line segment. */
        double gap_min_P = 1e30;
        double gap_max_theta = 0.0;
        int n_samples = 200;
        for (int s = 1; s < n_samples - 1; s++) {
            /* Parameter along the line from pk1 to pk2 */
            double frac = (double)s / (n_samples - 1);
            /* Interpolate in grid coords (handling wrapping) */
            double gi = pk1.i + frac * (pk2.i - pk1.i);
            double gj = pk1.j + frac * (pk2.j - pk1.j);
            double gk = pk1.k + frac * (pk2.k - pk1.k);
            /* Nearest grid point */
            int ii = ((int)round(gi) % N + N) % N;
            int jj = ((int)round(gj) % N + N) % N;
            int kk = ((int)round(gk) % N + N) % N;

            long gidx = idx3(ii, jj, kk, N);
            double P = fabs((double)phi[0][gidx] * phi[1][gidx] * phi[2][gidx]);
            if (P < gap_min_P) gap_min_P = P;

            double th_rms = 0.0;
            for (int a = 0; a < 3; a++)
                th_rms += (double)theta[a][gidx] * theta[a][gidx];
            th_rms = sqrt(th_rms);
            if (th_rms > gap_max_theta) gap_max_theta = th_rms;
        }

        printf("%d\t%.4f\t%.3f\t%.3f\t%.3f\t%.6f\t%.3f\t%.3f\t%.3f\t%.6f\t%.4f\t%.6f\t%.6f\n",
               fr, t, x1, y1, z1, pk1.P, x2, y2, z2, pk2.P, sep, gap_min_P, gap_max_theta);

        fprintf(stderr, "  frame %d  t=%.2f  sep=%.3f  pk1=%.5f  pk2=%.5f\n",
                fr, t, sep, pk1.P, pk2.P);
    }

    free(buf);
    sfa_close(sfa);
    fprintf(stderr, "Done.\n");
    return 0;
}
