/*  mismatch_profile.c — Extract Cosserat mismatch |M|² along collision axis
 *
 *  M_a = curl(φ)_a/2 - θ_a  (Cosserat geometric mismatch)
 *  |M|² = M_0² + M_1² + M_2²
 *
 *  Extracts 1D profiles along the x-axis (y=z=center) for specified frames.
 *  Also computes |curl(φ)|², |θ|², and P for context.
 *
 *  Build: gcc -O3 -o mismatch_profile mismatch_profile.c -lzstd -lm
 *  Usage: ./mismatch_profile input.sfa [frame1,frame2,...] > profile.tsv
 */

#define SFA_IMPLEMENTATION
#include "../sfa/format/sfa.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

static float f16f(uint16_t h) {
    uint16_t s = h & 0x8000; int e = (h >> 10) & 0x1F; uint16_t m = h & 0x3FF;
    if (e == 0) return 0; if (e == 31) return s ? -1e30f : 1e30f;
    float f; uint32_t x = ((uint32_t)s << 16) | ((uint32_t)(e-15+127) << 23) | ((uint32_t)m << 13);
    memcpy(&f, &x, 4); return f;
}

/* Extract column c, voxel i as float from frame buffer */
static float col_f(void *buf, SFA *sfa, int c, long i) {
    long N3 = (long)sfa->Nx * sfa->Ny * sfa->Nz;
    uint64_t off = 0;
    for (int cc = 0; cc < c; cc++)
        off += (uint64_t)N3 * sfa_dtype_size[sfa->columns[cc].dtype];
    int dt = sfa->columns[c].dtype;
    uint8_t *src = (uint8_t*)buf + off;
    if (dt == SFA_F16) return f16f(((uint16_t*)src)[i]);
    if (dt == SFA_F32) return ((float*)src)[i];
    return (float)((double*)src)[i];
}

int main(int argc, char **argv) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s input.sfa [frame1,frame2,...]\n", argv[0]);
        return 1;
    }

    SFA *sfa = sfa_open(argv[1]);
    if (!sfa) { fprintf(stderr, "Cannot open %s\n", argv[1]); return 1; }

    int N = sfa->Nx;
    long N3 = (long)N * N * N;
    int NN = N * N;
    double L = sfa->Lx;
    double dx = 2.0 * L / (N - 1);
    double idx1 = 1.0 / (2.0 * dx);

    /* Parse frame list */
    int frames[256], nf = 0;
    if (argc > 2) {
        char *s = argv[2];
        while (*s && nf < 256) {
            frames[nf++] = atoi(s);
            while (*s && *s != ',') s++;
            if (*s == ',') s++;
        }
    } else {
        /* Default: every 5th frame */
        for (uint32_t f = 0; f < sfa->total_frames && nf < 256; f += 5)
            frames[nf++] = f;
    }

    void *buf = malloc(sfa->frame_bytes);
    if (!buf) { fprintf(stderr, "Cannot allocate frame buffer\n"); return 1; }

    /* Need columns: phi_x(0), phi_y(1), phi_z(2), theta_x(3), theta_y(4), theta_z(5) */
    if (sfa->n_columns < 6) { fprintf(stderr, "Need 6+ columns\n"); return 1; }

    /* Header */
    printf("frame\ttime\tx\tM_sq\tcurl_phi_sq\ttheta_sq\tP\tphi_sq\n");

    int cj = N / 2, ck = N / 2;  /* Central y,z line */

    for (int fi = 0; fi < nf; fi++) {
        int frame = frames[fi];
        if (frame < 0 || frame >= (int)sfa->total_frames) continue;

        if (sfa_read_frame(sfa, frame, buf) != 0) {
            fprintf(stderr, "Cannot read frame %d\n", frame);
            continue;
        }
        double t = sfa_frame_time(sfa, frame);

        /* Walk along x-axis (varying i, fixed j=N/2, k=N/2) */
        for (int i = 2; i < N - 2; i++) {
            long idx = (long)i * NN + cj * N + ck;
            double x_pos = -L + i * dx;

            /* Neighbor indices for curl computation */
            int ip = i + 1, im = i - 1;
            int jp = cj + 1, jm = cj - 1;
            int kp = ck + 1, km = ck - 1;

            long n_ip = (long)ip * NN + cj * N + ck;
            long n_im = (long)im * NN + cj * N + ck;
            long n_jp = (long)i * NN + jp * N + ck;
            long n_jm = (long)i * NN + jm * N + ck;
            long n_kp = (long)i * NN + cj * N + kp;
            long n_km = (long)i * NN + cj * N + km;

            /* curl(phi) = (dφ2/dy - dφ1/dz, dφ0/dz - dφ2/dx, dφ1/dx - dφ0/dy) */
            float cp0 = (col_f(buf, sfa, 2, n_jp) - col_f(buf, sfa, 2, n_jm)
                       - col_f(buf, sfa, 1, n_kp) + col_f(buf, sfa, 1, n_km)) * idx1;
            float cp1 = (col_f(buf, sfa, 0, n_kp) - col_f(buf, sfa, 0, n_km)
                       - col_f(buf, sfa, 2, n_ip) + col_f(buf, sfa, 2, n_im)) * idx1;
            float cp2 = (col_f(buf, sfa, 1, n_ip) - col_f(buf, sfa, 1, n_im)
                       - col_f(buf, sfa, 0, n_jp) + col_f(buf, sfa, 0, n_jm)) * idx1;

            /* theta at center */
            float t0 = col_f(buf, sfa, 3, idx);
            float t1 = col_f(buf, sfa, 4, idx);
            float t2 = col_f(buf, sfa, 5, idx);

            /* Mismatch: M_a = curl(phi)_a/2 - theta_a */
            float M0 = cp0 * 0.5f - t0;
            float M1 = cp1 * 0.5f - t1;
            float M2 = cp2 * 0.5f - t2;

            float M_sq = M0 * M0 + M1 * M1 + M2 * M2;
            float curl_sq = cp0 * cp0 + cp1 * cp1 + cp2 * cp2;
            float theta_sq = t0 * t0 + t1 * t1 + t2 * t2;

            /* phi and P */
            float p0 = col_f(buf, sfa, 0, idx);
            float p1 = col_f(buf, sfa, 1, idx);
            float p2 = col_f(buf, sfa, 2, idx);
            float P = p0 * p1 * p2;
            float phi_sq = p0 * p0 + p1 * p1 + p2 * p2;

            printf("%d\t%.4f\t%.3f\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\n",
                   frame, t, x_pos, M_sq, curl_sq, theta_sq, P, phi_sq);
        }
    }

    free(buf);
    sfa_close(sfa);
    return 0;
}
