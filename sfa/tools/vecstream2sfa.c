/*  vecstream2sfa.c — Convert vecstream files back to SFA voxel format
 *
 *  Usage: vecstream2sfa input.vecstream -o output.sfa [-field 0] [-start 0] [-end -1]
 *
 *  Build: gcc -O3 -march=native -fopenmp -o vecstream2sfa vecstream2sfa.c -lzstd -lm
 */

#define SFA_IMPLEMENTATION
#include "../format/sfa.h"

#define VECSTREAM_IMPLEMENTATION
#include "../format/vecstream.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

int main(int argc, char **argv) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s input.vecstream -o output.sfa [-field 0] [-start 0] [-end -1]\n",
                argv[0]);
        return 1;
    }

    const char *input_path = argv[1];
    const char *output_path = NULL;
    int target_field = 0;
    int start_frame = 0;
    int end_frame = -1;

    for (int i = 2; i < argc; i++) {
        if (!strcmp(argv[i], "-o") && i+1 < argc) output_path = argv[++i];
        else if (!strcmp(argv[i], "-field") && i+1 < argc) target_field = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-start") && i+1 < argc) start_frame = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-end") && i+1 < argc) end_frame = atoi(argv[++i]);
    }

    if (!output_path) {
        fprintf(stderr, "Error: -o output.sfa required\n");
        return 1;
    }

    /* Open vecstream */
    VecStream *vs = vecstream_open(input_path);
    if (!vs) { fprintf(stderr, "Cannot open %s\n", input_path); return 1; }

    uint32_t Nx = vs->Nx, Ny = vs->Ny, Nz = vs->Nz;
    uint64_t N3 = (uint64_t)Nx * Ny * Nz;
    uint32_t total_frames = vecstream_frame_count(vs);

    printf("vecstream2sfa: %s -> %s\n", input_path, output_path);
    printf("  Grid: %u x %u x %u (%lu voxels)\n", Nx, Ny, Nz, (unsigned long)N3);
    printf("  Frames in vecstream: %u\n", total_frames);
    printf("  Target field: %d\n", target_field);

    /* Gather unique simulation times for the target field */
    /* We need to identify the logical frames: sets of entries that form
       an I->P->P->... chain for this field */
    double *times = (double *)malloc(total_frames * sizeof(double));
    int *frame_map = (int *)malloc(total_frames * sizeof(int)); /* vecstream index -> logical */
    int n_logical = 0;

    /* Find all frames for this field, in order */
    int *field_frames = (int *)malloc(total_frames * sizeof(int));
    int n_field_frames = 0;
    for (uint32_t i = 0; i < total_frames; i++) {
        if (vecstream_frame_field(vs, i) == target_field) {
            uint8_t ft = vecstream_frame_type(vs, i);
            if (ft == VS_FRAME_I || ft == VS_FRAME_P) {
                field_frames[n_field_frames] = (int)i;
                times[n_field_frames] = vecstream_frame_time(vs, i);
                n_field_frames++;
            }
        }
    }

    if (end_frame < 0 || end_frame >= n_field_frames)
        end_frame = n_field_frames - 1;
    if (start_frame < 0) start_frame = 0;
    if (start_frame > end_frame) {
        fprintf(stderr, "No frames in range [%d, %d]\n", start_frame, end_frame);
        vecstream_close(vs);
        return 1;
    }

    int out_frames = end_frame - start_frame + 1;
    printf("  Reconstructing frames %d..%d (%d frames)\n", start_frame, end_frame, out_frames);

    /* Create output SFA */
    SFA *sfa = sfa_create(output_path, Nx, Ny, Nz,
                          vs->Lx, vs->Ly, vs->Lz,
                          vs->dt > 0 ? vs->dt : 1.0);
    if (!sfa) { fprintf(stderr, "Cannot create %s\n", output_path); vecstream_close(vs); return 1; }

    char col_name[16];
    snprintf(col_name, sizeof(col_name), "field_%d", target_field);
    sfa_add_column(sfa, col_name, SFA_F32, SFA_POSITION, (uint8_t)target_field);
    sfa_finalize_header(sfa);

    /* Reconstruct each frame */
    float *voxels = (float *)malloc(N3 * sizeof(float));

    for (int li = start_frame; li <= end_frame; li++) {
        int vs_idx = field_frames[li];
        double t = vecstream_frame_time(vs, vs_idx);

        if (vecstream_reconstruct(vs, (uint32_t)vs_idx, (uint8_t)target_field, voxels) != 0) {
            fprintf(stderr, "  Failed to reconstruct frame %d (vs_idx=%d)\n", li, vs_idx);
            continue;
        }

        void *col_ptrs[1] = { voxels };
        sfa_write_frame(sfa, t, col_ptrs);

        if (li == start_frame || li == end_frame || (li - start_frame) % 10 == 0)
            printf("  Wrote frame %d/%d (t=%.4f)\n", li - start_frame + 1, out_frames, t);
    }

    sfa_close(sfa);
    vecstream_close(vs);
    free(voxels);
    free(times);
    free(frame_map);
    free(field_frames);

    printf("Done. Output: %s\n", output_path);
    return 0;
}
