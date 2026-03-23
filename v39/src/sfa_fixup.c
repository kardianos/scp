/*  sfa_fixup.c — Fix streaming SFA files that have total_frames=0
 *
 *  Build: gcc -O3 -o sfa_fixup src/sfa_fixup.c -lzstd
 *  Usage: ./sfa_fixup <file.sfa>
 */

#define SFA_IMPLEMENTATION
#include "../../sfa/format/sfa.h"

#include <stdio.h>

int main(int argc, char **argv) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <file.sfa>\n", argv[0]);
        return 1;
    }

    printf("Fixing index for: %s\n", argv[1]);
    int rc = sfa_fixup_index(argv[1]);
    if (rc < 0) {
        fprintf(stderr, "fixup failed with code %d\n", rc);
        return 1;
    }
    printf("Success: fixed %d frames\n", rc);

    /* Verify by opening */
    SFA *s = sfa_open(argv[1]);
    if (s) {
        printf("Verification: %d frames, grid %dx%dx%d\n", s->total_frames, s->Nx, s->Ny, s->Nz);
        sfa_close(s);
    }
    return 0;
}
