/*  theta_avg_compare.c — Time-averaged signed theta_phi for winding comparison
 *
 *  Averages the signed theta_phi over many frames to cancel oscillatory part,
 *  leaving only the DC (steady-state) component.
 *
 *  Build: gcc -O2 -o theta_avg_compare src/theta_avg_compare.c -lzstd -lm
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define SFA_IMPLEMENTATION
#include "../../viewer/sfa.h"

#define MAX_SHELLS 100

int main(int argc, char **argv) {
    int f_start = 10, f_end = -1, f_step = 1;
    char *sfapath = NULL;
    char *outpath = NULL;

    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "-fs") && i+1 < argc) f_start = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-fe") && i+1 < argc) f_end = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-step") && i+1 < argc) f_step = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-o") && i+1 < argc) outpath = argv[++i];
        else sfapath = argv[i];
    }
    if (!sfapath) {
        fprintf(stderr, "Usage: %s [-fs start] [-fe end] [-step step] [-o outfile] <archive.sfa>\n", argv[0]);
        return 1;
    }

    SFA *s = sfa_open(sfapath);
    if (!s) { fprintf(stderr, "Cannot open %s\n", sfapath); return 1; }

    if (f_end < 0) f_end = s->total_frames - 1;
    if ((uint32_t)f_end >= s->total_frames) f_end = s->total_frames - 1;

    uint64_t N_total = (uint64_t)s->Nx * s->Ny * s->Nz;
    void *buf = malloc(s->frame_bytes);
    int Nx = s->Nx, Ny = s->Ny, Nz = s->Nz;
    double Lx = s->Lx, Ly = s->Ly;
    double dx = 2.0 * Lx / (Nx - 1);
    double dy = 2.0 * Ly / (Ny - 1);

    double dr = 0.5;
    int n_shells = (int)(20.0 / dr) + 1;
    if (n_shells > MAX_SHELLS) n_shells = MAX_SHELLS;

    /* Accumulators for time average */
    double *sum_tphi = calloc(n_shells, sizeof(double));
    double *sum_tphi2 = calloc(n_shells, sizeof(double));
    long *total_count = calloc(n_shells, sizeof(long));
    int n_frames = 0;

    printf("Archive: %s (%u frames)\n", sfapath, s->total_frames);
    printf("Averaging frames %d to %d step %d\n", f_start, f_end, f_step);

    for (int frame = f_start; frame <= f_end; frame += f_step) {
        if (sfa_read_frame(s, frame, buf) < 0) continue;
        n_frames++;

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

        /* Accumulate in shells */
        for (int i = 0; i < Nx; i++) {
            double x = -Lx + i * dx;
            for (int j = 0; j < Ny; j++) {
                double y = -Ly + j * dy;
                double xc = x - cx, yc = y - cy;
                double r_perp = sqrt(xc*xc + yc*yc);
                int shell = (int)(r_perp / dr);
                if (shell >= n_shells) continue;

                double angle = atan2(yc, xc);
                double sa = sin(angle), ca = cos(angle);

                for (int k = 0; k < Nz; k++) {
                    long idx = (long)i * (Ny * Nz) + j * Nz + k;
                    double t_phi = -theta_x[idx] * sa + theta_y[idx] * ca;
                    sum_tphi[shell] += t_phi;
                    sum_tphi2[shell] += t_phi * t_phi;
                    total_count[shell]++;
                }
            }
        }
    }

    printf("Processed %d frames\n\n", n_frames);

    /* Output */
    if (!outpath) outpath = "data/theta_avg.tsv";
    FILE *fp = fopen(outpath, "w");
    fprintf(fp, "r\tmean_theta_phi\trms_theta_phi\tratio_mean_rms\n");

    printf("%6s %14s %14s %10s\n", "r", "mean(tphi)", "rms(tphi)", "|m|/rms");
    printf("--------------------------------------------------\n");
    for (int si = 0; si < n_shells; si++) {
        if (total_count[si] == 0) continue;
        double r = (si + 0.5) * dr;
        double n = (double)total_count[si];
        double mean = sum_tphi[si] / n;
        double rms = sqrt(sum_tphi2[si] / n);
        double ratio = (rms > 0) ? fabs(mean) / rms : 0;
        fprintf(fp, "%.2f\t%.6e\t%.6e\t%.4f\n", r, mean, rms, ratio);
        if (r <= 15.25)
            printf("%6.2f %14.4e %14.4e %10.4f\n", r, mean, rms, ratio);
    }

    fclose(fp);
    printf("\nWritten: %s\n", outpath);

    free(sum_tphi); free(sum_tphi2); free(total_count);
    free(buf);
    sfa_close(s);
    return 0;
}
