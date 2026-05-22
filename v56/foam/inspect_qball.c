/*  inspect_qball.c — dump central-cell values from cell-native SFA */
#define SFA_IMPLEMENTATION
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../../sfa/format/sfa.h"

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage: %s <file.sfa>\n", argv[0]); return 1; }
    SFA *s = sfa_open(argv[1]);
    if (!s) { fprintf(stderr, "cannot open %s\n", argv[1]); return 1; }

    printf("file: %s\n", argv[1]);
    printf("n_columns=%u  total_frames=%u\n", s->n_columns, s->total_frames);
    for (uint32_t c = 0; c < s->n_columns; c++)
        printf("  col[%u] name='%s' semantic=%u component=%u dtype=%u\n",
               c, s->columns[c].name, s->columns[c].semantic,
               s->columns[c].component, s->columns[c].dtype);

    /* Walk frames sequentially via the FCEL reader. */
    int snap_idx = 0;
    for (uint32_t f = 0; f < s->total_frames; f++) {
        uint32_t N_cells = 0, n_cols = 0;
        uint8_t dtype = 0;
        void *data = NULL;
        int rc = sfa_read_cell_frame(s, f, &N_cells, &n_cols, &dtype, &data);
        if (rc != 0 || data == NULL) continue;
        if (snap_idx == 0 || snap_idx == 5 || snap_idx == 10 || snap_idx == 15 || snap_idx == 20) {
            printf("\n[snap %2d] frame=%u  N_cells=%u  n_cols=%u  dtype=%u\n",
                   snap_idx, f, N_cells, n_cols, dtype);
            /* Buffer layout: n_cols × N_cells, column-major, dtype elements. */
            for (uint32_t c = 0; c < n_cols; c++) {
                double sum = 0, sumsq = 0, mx = 0;
                for (uint32_t i = 0; i < N_cells; i++) {
                    float v = sfa_read_voxel_f32(data, s, c, i);
                    sum   += v;
                    sumsq += v*v;
                    if (fabsf(v) > mx) mx = fabsf(v);
                }
                printf("    col[%u] %-12s  mean=%+.4f  RMS=%.4f  |max|=%.4f\n",
                       c, s->columns[c].name, sum/N_cells, sqrt(sumsq/N_cells), mx);
            }
        }
        snap_idx++;
        free(data);
    }
    sfa_close(s);
    return 0;
}
