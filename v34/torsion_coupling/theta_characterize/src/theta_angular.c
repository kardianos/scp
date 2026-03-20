/*  theta_angular.c — Angular structure of theta_phi around braid axis
 *
 *  Checks whether theta_phi has a consistent sign (DC circulation)
 *  or is oscillatory (wave pattern) around the braid.
 *
 *  Build: gcc -O2 -o theta_angular src/theta_angular.c -lzstd -lm
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define SFA_IMPLEMENTATION
#include "../../viewer/sfa.h"

#define N_ANG 24  /* angular bins */
#define PI 3.14159265358979323846

int main(int argc, char **argv) {
    int frame_idx = 200;
    char *sfapath = NULL;

    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "-f") && i+1 < argc) frame_idx = atoi(argv[++i]);
        else sfapath = argv[i];
    }
    if (!sfapath) {
        fprintf(stderr, "Usage: %s [-f frame] <archive.sfa>\n", argv[0]);
        return 1;
    }

    SFA *s = sfa_open(sfapath);
    if (!s) { fprintf(stderr, "Cannot open %s\n", sfapath); return 1; }

    uint64_t N_total = (uint64_t)s->Nx * s->Ny * s->Nz;
    void *buf = malloc(s->frame_bytes);
    if (sfa_read_frame(s, frame_idx, buf) < 0) {
        fprintf(stderr, "Failed to read frame %d\n", frame_idx);
        free(buf); sfa_close(s); return 1;
    }
    double t = sfa_frame_time(s, frame_idx);
    printf("Frame %d: t = %.2f\n", frame_idx, t);

    double *data = (double *)buf;
    double *phi_x   = data + 0 * N_total;
    double *phi_y   = data + 1 * N_total;
    double *phi_z   = data + 2 * N_total;
    double *theta_x = data + 3 * N_total;
    double *theta_y = data + 4 * N_total;
    double *theta_z = data + 5 * N_total;

    int Nx = s->Nx, Ny = s->Ny, Nz = s->Nz;
    double Lx = s->Lx, Ly = s->Ly;
    double dx = 2.0 * Lx / (Nx - 1);
    double dy = 2.0 * Ly / (Ny - 1);

    /* Find braid center */
    double avg_phi2 = 0;
    for (uint64_t idx = 0; idx < N_total; idx++)
        avg_phi2 += phi_x[idx]*phi_x[idx] + phi_y[idx]*phi_y[idx] + phi_z[idx]*phi_z[idx];
    avg_phi2 /= N_total;

    double thresh = 5.0 * avg_phi2;
    double wx = 0, wy = 0, wt = 0;
    for (int i = 0; i < Nx; i++) {
        double x = -Lx + i * dx;
        for (int j = 0; j < Ny; j++) {
            double y = -Ly + j * dy;
            for (int k = 0; k < Nz; k++) {
                long idx = (long)i * (Ny * Nz) + j * Nz + k;
                double p2 = phi_x[idx]*phi_x[idx] + phi_y[idx]*phi_y[idx] + phi_z[idx]*phi_z[idx];
                if (p2 < thresh) continue;
                wx += x * p2; wy += y * p2; wt += p2;
            }
        }
    }
    double cx = wx / wt, cy = wy / wt;
    printf("Braid center: (%.3f, %.3f)\n", cx, cy);

    /* Angular bins of signed theta_phi at various radial distances */
    double r_bands[] = {2.0, 4.0, 6.0, 8.0, 10.0};
    int n_bands = 5;
    double band_width = 1.0;

    FILE *fp = fopen("data/theta_angular.tsv", "w");
    fprintf(fp, "angle_deg");
    for (int b = 0; b < n_bands; b++) fprintf(fp, "\ttheta_phi_r%.0f", r_bands[b]);
    fprintf(fp, "\n");

    double sum_phi[5][N_ANG];
    long cnt[5][N_ANG];
    memset(sum_phi, 0, sizeof(sum_phi));
    memset(cnt, 0, sizeof(cnt));

    for (int i = 0; i < Nx; i++) {
        double x = -Lx + i * dx;
        for (int j = 0; j < Ny; j++) {
            double y = -Ly + j * dy;
            double xc = x - cx, yc = y - cy;
            double r_perp = sqrt(xc*xc + yc*yc);
            double angle = atan2(yc, xc);
            double ca = cos(angle), sa = sin(angle);
            int ang_bin = (int)((angle + PI) / (2*PI) * N_ANG);
            if (ang_bin >= N_ANG) ang_bin = N_ANG - 1;

            for (int b = 0; b < n_bands; b++) {
                if (fabs(r_perp - r_bands[b]) > band_width / 2) continue;
                for (int k = 0; k < Nz; k++) {
                    long idx = (long)i * (Ny * Nz) + j * Nz + k;
                    double tx = theta_x[idx], ty = theta_y[idx];
                    double t_phi = -tx * sa + ty * ca;
                    sum_phi[b][ang_bin] += t_phi;
                    cnt[b][ang_bin]++;
                }
            }
        }
    }

    for (int a = 0; a < N_ANG; a++) {
        double deg = -180.0 + (a + 0.5) * 360.0 / N_ANG;
        fprintf(fp, "%.1f", deg);
        for (int b = 0; b < n_bands; b++) {
            double val = (cnt[b][a] > 0) ? sum_phi[b][a] / cnt[b][a] : 0;
            fprintf(fp, "\t%.6e", val);
        }
        fprintf(fp, "\n");
    }
    fclose(fp);

    /* Print summary: does theta_phi have consistent sign? */
    printf("\n=== Angular consistency of theta_phi ===\n");
    for (int b = 0; b < n_bands; b++) {
        int n_pos = 0, n_neg = 0;
        double mean = 0, rms = 0;
        int n_valid = 0;
        for (int a = 0; a < N_ANG; a++) {
            if (cnt[b][a] == 0) continue;
            double val = sum_phi[b][a] / cnt[b][a];
            mean += val;
            rms += val * val;
            n_valid++;
            if (val > 0) n_pos++; else n_neg++;
        }
        mean /= n_valid;
        rms = sqrt(rms / n_valid);
        printf("r=%.0f: mean(theta_phi) = %+.4e, rms = %.4e, ratio |mean|/rms = %.3f, "
               "sign: %d+/%d- of %d bins\n",
               r_bands[b], mean, rms, fabs(mean)/rms, n_pos, n_neg, n_valid);
    }

    free(buf);
    sfa_close(s);
    return 0;
}
