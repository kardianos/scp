/* sfa_test.c — Test SFA write + read roundtrip */
#define SFA_IMPLEMENTATION
#include "sfa.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

int main(void) {
    int Nx=32, Ny=32, Nz=32;
    int N_total = Nx * Ny * Nz;

    printf("=== SFA Write/Read Test ===\n");
    printf("Grid: %d×%d×%d = %d points\n\n", Nx, Ny, Nz, N_total);

    /* Create test data: 3 columns of f64 */
    double *phi_x = malloc(N_total * sizeof(double));
    double *phi_y = malloc(N_total * sizeof(double));
    double *phi_z = malloc(N_total * sizeof(double));

    /* ---- Write ---- */
    printf("Writing test.sfa...\n");
    SFA *w = sfa_create("test.sfa", Nx, Ny, Nz, 10.0, 10.0, 10.0, 0.01);
    sfa_add_column(w, "phi_x", SFA_F64, SFA_POSITION, 0);
    sfa_add_column(w, "phi_y", SFA_F64, SFA_POSITION, 1);
    sfa_add_column(w, "phi_z", SFA_F64, SFA_POSITION, 2);
    sfa_finalize_header(w);

    int n_frames = 5;
    for (int f = 0; f < n_frames; f++) {
        double t = f * 10.0;
        /* Fill with test pattern */
        for (int ix = 0; ix < Nx; ix++)
            for (int iy = 0; iy < Ny; iy++)
                for (int iz = 0; iz < Nz; iz++) {
                    int idx = ix*(Ny*Nz) + iy*Nz + iz;
                    double x = -10.0 + ix * 20.0/(Nx-1);
                    double y = -10.0 + iy * 20.0/(Ny-1);
                    double z = -10.0 + iz * 20.0/(Nz-1);
                    phi_x[idx] = sin(x + t*0.1) * exp(-(y*y+z*z)/25.0);
                    phi_y[idx] = cos(y + t*0.1) * exp(-(x*x+z*z)/25.0);
                    phi_z[idx] = sin(z + t*0.1) * exp(-(x*x+y*y)/25.0);
                }
        void *cols[] = {phi_x, phi_y, phi_z};
        sfa_write_frame(w, t, cols);
        printf("  Frame %d: t=%.1f written\n", f, t);
    }
    sfa_close(w);

    /* ---- Read ---- */
    printf("\nReading test.sfa...\n");
    SFA *r = sfa_open("test.sfa");
    if (!r) { fprintf(stderr, "Failed to open\n"); return 1; }

    printf("  Version: %d\n", r->version);
    printf("  Grid: %d×%d×%d\n", r->Nx, r->Ny, r->Nz);
    printf("  Domain: [%.1f]×[%.1f]×[%.1f]\n", r->Lx, r->Ly, r->Lz);
    printf("  Columns: %d\n", r->n_columns);
    for (uint32_t c = 0; c < r->n_columns; c++)
        printf("    %d: \"%s\" %s (sem=%d, comp=%d)\n",
               c, r->columns[c].name, sfa_dtype_name(r->columns[c].dtype),
               r->columns[c].semantic, r->columns[c].component);
    printf("  Frames: %d\n", r->total_frames);
    printf("  Frame bytes: %lu\n", (unsigned long)r->frame_bytes);

    /* Read each frame and verify */
    void *buf = malloc(r->frame_bytes);
    int errors = 0;
    for (uint32_t f = 0; f < r->total_frames; f++) {
        double t = sfa_frame_time(r, f);
        int rc = sfa_read_frame(r, f, buf);
        if (rc < 0) { printf("  Frame %d: READ ERROR %d\n", f, rc); errors++; continue; }

        /* Verify a few values */
        double *rd_phi_x = (double *)buf;
        /* Regenerate expected value at (ix=5, iy=5, iz=5) */
        int idx = 5*(Ny*Nz) + 5*Nz + 5;
        double x = -10.0 + 5 * 20.0/(Nx-1);
        double y = -10.0 + 5 * 20.0/(Ny-1);
        double z = -10.0 + 5 * 20.0/(Nz-1);
        double expected = sin(x + t*0.1) * exp(-(y*y+z*z)/25.0);
        double got = rd_phi_x[idx];
        double err = fabs(got - expected);
        printf("  Frame %d: t=%.1f, phi_x[5,5,5] = %.8f (expected %.8f, err=%.1e) %s\n",
               f, t, got, expected, err, (err < 1e-12) ? "OK" : "MISMATCH");
        if (err > 1e-12) errors++;
    }

    free(buf);
    sfa_close(r);

    /* File size comparison */
    FILE *fp = fopen("test.sfa", "rb");
    fseek(fp, 0, SEEK_END);
    long fsize = ftell(fp);
    fclose(fp);
    long raw_size = (long)n_frames * N_total * 3 * 8;
    printf("\nFile size: %ld bytes (%.1f KB)\n", fsize, fsize/1024.0);
    printf("Raw size:  %ld bytes (%.1f KB)\n", raw_size, raw_size/1024.0);
    printf("Ratio:     %.1f×\n", (double)raw_size / fsize);
    printf("\n%s\n", errors ? "ERRORS DETECTED" : "ALL TESTS PASSED");

    free(phi_x); free(phi_y); free(phi_z);
    return errors;
}
