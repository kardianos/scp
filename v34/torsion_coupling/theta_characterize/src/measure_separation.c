/*  measure_separation.c — Track braid separation D(t) from SFA archive
 *
 *  For two-braid simulations, finds both braid centroids and measures
 *  the separation D(t) over all frames.
 *
 *  Build: gcc -O2 -o measure_separation src/measure_separation.c -lzstd -lm
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define SFA_IMPLEMENTATION
#include "../../viewer/sfa.h"

int main(int argc, char **argv) {
    char *sfapath = NULL;
    char *outpath = NULL;

    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "-o") && i+1 < argc) outpath = argv[++i];
        else sfapath = argv[i];
    }
    if (!sfapath) {
        fprintf(stderr, "Usage: %s [-o outfile] <archive.sfa>\n", argv[0]);
        return 1;
    }

    SFA *s = sfa_open(sfapath);
    if (!s) { fprintf(stderr, "Cannot open %s\n", sfapath); return 1; }

    printf("Archive: %s\n", sfapath);
    printf("Grid: %u x %u x %u, L=(%.1f, %.1f, %.1f)\n",
           s->Nx, s->Ny, s->Nz, s->Lx, s->Ly, s->Lz);
    printf("Frames: %u\n\n", s->total_frames);

    uint64_t N_total = (uint64_t)s->Nx * s->Ny * s->Nz;
    void *buf = malloc(s->frame_bytes);
    int Nx = s->Nx, Ny = s->Ny, Nz = s->Nz;
    double Lx = s->Lx, Ly = s->Ly;
    double dx = 2.0 * Lx / (Nx - 1);
    double dy = 2.0 * Ly / (Ny - 1);

    if (!outpath) outpath = "data/separation.tsv";
    FILE *fp = fopen(outpath, "w");
    fprintf(fp, "frame\tt\tx1\ty1\tx2\ty2\tD\n");

    printf("%6s %8s %8s %8s %8s %8s %8s\n", "frame", "t", "x1", "y1", "x2", "y2", "D");
    printf("---------------------------------------------------------------\n");

    for (uint32_t frame = 0; frame < s->total_frames; frame++) {
        if (sfa_read_frame(s, frame, buf) < 0) continue;
        double t = sfa_frame_time(s, frame);

        double *data = (double *)buf;
        double *phi_x = data + 0 * N_total;
        double *phi_y = data + 1 * N_total;
        double *phi_z = data + 2 * N_total;

        /* Compute average phi^2 */
        double avg = 0;
        for (uint64_t idx = 0; idx < N_total; idx++)
            avg += phi_x[idx]*phi_x[idx] + phi_y[idx]*phi_y[idx] + phi_z[idx]*phi_z[idx];
        avg /= N_total;
        double thresh = 5.0 * avg;

        /* Find overall centroid first */
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
        /* don't need overall centroid, use x=0 as split */

        /* Split into left (x<0) and right (x>0) half */
        double wx1 = 0, wy1 = 0, wt1 = 0;  /* x < 0 */
        double wx2 = 0, wy2 = 0, wt2 = 0;  /* x > 0 */

        for (int i = 0; i < Nx; i++) {
            double x = -Lx + i * dx;
            for (int j = 0; j < Ny; j++) {
                double y = -Ly + j * dy;
                for (int k = 0; k < Nz; k++) {
                    long idx = (long)i * (Ny * Nz) + j * Nz + k;
                    double p2 = phi_x[idx]*phi_x[idx] + phi_y[idx]*phi_y[idx] + phi_z[idx]*phi_z[idx];
                    if (p2 < thresh) continue;
                    if (x < 0) {
                        wx1 += x * p2; wy1 += y * p2; wt1 += p2;
                    } else {
                        wx2 += x * p2; wy2 += y * p2; wt2 += p2;
                    }
                }
            }
        }

        double x1 = (wt1 > 0) ? wx1/wt1 : -7.5;
        double y1 = (wt1 > 0) ? wy1/wt1 : 0;
        double x2 = (wt2 > 0) ? wx2/wt2 : 7.5;
        double y2 = (wt2 > 0) ? wy2/wt2 : 0;
        double D = sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1));

        fprintf(fp, "%u\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n",
                frame, t, x1, y1, x2, y2, D);
        printf("%6u %8.2f %8.3f %8.3f %8.3f %8.3f %8.3f\n",
               frame, t, x1, y1, x2, y2, D);
    }

    fclose(fp);
    printf("\nWritten: %s\n", outpath);
    free(buf);
    sfa_close(s);
    return 0;
}
