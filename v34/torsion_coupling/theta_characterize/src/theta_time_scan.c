/*  theta_time_scan.c — Scan theta_phi at fixed r across multiple frames
 *
 *  Checks whether the azimuthal circulation is oscillating in time.
 *
 *  Build: gcc -O2 -o theta_time_scan src/theta_time_scan.c -lzstd -lm
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define SFA_IMPLEMENTATION
#include "../../viewer/sfa.h"

int main(int argc, char **argv) {
    int f_start = 50, f_end = 260, f_step = 10;
    char *sfapath = NULL;

    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "-fs") && i+1 < argc) f_start = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-fe") && i+1 < argc) f_end = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-step") && i+1 < argc) f_step = atoi(argv[++i]);
        else sfapath = argv[i];
    }
    if (!sfapath) {
        fprintf(stderr, "Usage: %s [-fs start] [-fe end] [-step step] <archive.sfa>\n", argv[0]);
        return 1;
    }

    SFA *s = sfa_open(sfapath);
    if (!s) { fprintf(stderr, "Cannot open %s\n", sfapath); return 1; }

    uint64_t N_total = (uint64_t)s->Nx * s->Ny * s->Nz;
    void *buf = malloc(s->frame_bytes);
    int Nx = s->Nx, Ny = s->Ny, Nz = s->Nz;
    double Lx = s->Lx, Ly = s->Ly;
    double dx = 2.0 * Lx / (Nx - 1);
    double dy = 2.0 * Ly / (Ny - 1);

    /* Probe radii */
    double r_probe[] = {2.0, 4.0, 6.0, 8.0};
    int n_probe = 4;
    double band_w = 1.0;

    FILE *fp = fopen("data/theta_time.tsv", "w");
    fprintf(fp, "frame\tt\ttheta_phi_r2\ttheta_phi_r4\ttheta_phi_r6\ttheta_phi_r8\n");

    printf("%6s %8s %12s %12s %12s %12s\n", "frame", "t", "tp_r2", "tp_r4", "tp_r6", "tp_r8");
    printf("----------------------------------------------------------------------\n");

    for (int frame = f_start; frame <= f_end && (uint32_t)frame < s->total_frames; frame += f_step) {
        if (sfa_read_frame(s, frame, buf) < 0) continue;
        double t = sfa_frame_time(s, frame);

        double *data = (double *)buf;
        double *phi_x   = data + 0 * N_total;
        double *phi_y   = data + 1 * N_total;
        double *phi_z   = data + 2 * N_total;
        double *theta_x = data + 3 * N_total;
        double *theta_y = data + 4 * N_total;

        /* Find braid center */
        double avg = 0;
        for (uint64_t idx = 0; idx < N_total; idx++)
            avg += phi_x[idx]*phi_x[idx] + phi_y[idx]*phi_y[idx] + phi_z[idx]*phi_z[idx];
        avg /= N_total;
        double thresh = 5.0 * avg;
        double wx = 0, wy = 0, wt = 0;
        for (int i = 0; i < Nx; i++) {
            double x = -Lx + i * dx;
            for (int j = 0; j < Ny; j++) {
                double y = -Ly + j * dy;
                for (int k = 0; k < Nz; k++) {
                    long idx = (long)i * (Ny * Nz) + j * Nz + k;
                    double p2 = phi_x[idx]*phi_x[idx] + phi_y[idx]*phi_y[idx] + phi_z[idx]*phi_z[idx];
                    if (p2 < thresh) continue;
                    wx += x*p2; wy += y*p2; wt += p2;
                }
            }
        }
        double cx = wx / wt, cy = wy / wt;

        /* Average signed theta_phi in each radial band */
        double sum[4] = {0,0,0,0};
        long cnt[4] = {0,0,0,0};

        for (int i = 0; i < Nx; i++) {
            double x = -Lx + i * dx;
            for (int j = 0; j < Ny; j++) {
                double y = -Ly + j * dy;
                double xc = x - cx, yc = y - cy;
                double r_perp = sqrt(xc*xc + yc*yc);
                double angle = atan2(yc, xc);
                double sa = sin(angle), ca = cos(angle);

                for (int p = 0; p < n_probe; p++) {
                    if (fabs(r_perp - r_probe[p]) > band_w / 2) continue;
                    for (int k = 0; k < Nz; k++) {
                        long idx = (long)i * (Ny * Nz) + j * Nz + k;
                        double t_phi = -theta_x[idx] * sa + theta_y[idx] * ca;
                        sum[p] += t_phi;
                        cnt[p]++;
                    }
                }
            }
        }

        fprintf(fp, "%d\t%.3f", frame, t);
        printf("%6d %8.2f", frame, t);
        for (int p = 0; p < n_probe; p++) {
            double val = (cnt[p] > 0) ? sum[p] / cnt[p] : 0;
            fprintf(fp, "\t%.6e", val);
            printf(" %12.4e", val);
        }
        fprintf(fp, "\n");
        printf("\n");
    }

    fclose(fp);
    printf("\nWritten: data/theta_time.tsv\n");
    free(buf);
    sfa_close(s);
    return 0;
}
