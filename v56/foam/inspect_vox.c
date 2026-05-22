/*  inspect_vox.c — dump per-channel stats from a voxel SFA file (FRMD frames). */
#define SFA_IMPLEMENTATION
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../../sfa/format/sfa.h"

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage: %s <voxel.sfa>\n", argv[0]); return 1; }
    SFA *s = sfa_open(argv[1]);
    if (!s) return 1;
    printf("file: %s  Nx=%u Ny=%u Nz=%u  frames=%u  cols=%u\n",
           argv[1], s->Nx, s->Ny, s->Nz, s->total_frames, s->n_columns);
    for (uint32_t c = 0; c < s->n_columns; c++) {
        printf("  col[%u]='%s' sem=%u comp=%u dtype=%u\n",
               c, s->columns[c].name, s->columns[c].semantic,
               s->columns[c].component, s->columns[c].dtype);
    }
    uint64_t Ntot = (uint64_t)s->Nx * s->Ny * s->Nz;
    uint8_t *buf = malloc(s->frame_bytes);
    if (!buf) return 1;
    int show[] = {0, 5, 10, 15, 20, -1};
    for (int si = 0; show[si] >= 0; si++) {
        int f = show[si];
        if ((uint32_t)f >= s->total_frames) break;
        if (sfa_read_frame(s, f, buf) != 0) continue;
        double t = (double)f * s->dt;
        printf("\n[frame %2d] t=%.3f\n", f, t);
        for (uint32_t c = 0; c < s->n_columns; c++) {
            double sum=0, sumsq=0, mx=0;
            for (uint64_t i = 0; i < Ntot; i++) {
                float v = sfa_read_voxel_f32(buf, s, c, i);
                sum += v; sumsq += v*v;
                if (fabsf(v) > mx) mx = fabsf(v);
            }
            printf("   col[%u] %-10s  mean=%+.4f  RMS=%.4f  |max|=%.4f\n",
                   c, s->columns[c].name, sum/Ntot, sqrt(sumsq/Ntot), mx);
        }
    }
    free(buf);
    sfa_close(s);
    return 0;
}
