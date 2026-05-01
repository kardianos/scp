/*  strip_background.c — Remove analytical background from an SFA file.
 *
 *  Reads the last frame, subtracts A_bg*cos(k*z+delta) from each phi component
 *  and omega*A_bg*sin(k*z+delta) from each phi velocity. Writes as new seed.
 *
 *  Build: gcc -O3 -fopenmp -o strip_background strip_background.c -I../../sfa/format -lzstd -lm
 *  Usage: strip_background input.sfa output.sfa [A_bg] [m]
 */
#define SFA_IMPLEMENTATION
#include "../../sfa/format/sfa.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define PI 3.14159265358979323846

int main(int argc, char **argv) {
    if (argc < 3) {
        fprintf(stderr, "Usage: %s input.sfa output.sfa [A_bg] [m]\n", argv[0]);
        return 1;
    }
    double A_bg = argc > 3 ? atof(argv[3]) : 0.1;
    double m = argc > 4 ? atof(argv[4]) : 3.0;
    double m2 = m * m;

    SFA *in = sfa_open(argv[1]);
    if (!in) { fprintf(stderr, "Cannot open %s\n", argv[1]); return 1; }
    int N = in->Nx;
    long N3 = (long)N*N*N;
    int NN = N*N;
    double L = in->Lx;
    double dx = 2.0 * L / (N - 1);
    double k_bg = PI / L;
    double omega_bg = sqrt(k_bg*k_bg + m2);
    double delta[3] = {0, 3.0005, 4.4325};

    /* Read last frame */
    uint32_t last = in->total_frames - 1;
    /* Fix total_frames from JTOP */
    {
        uint64_t jt = in->first_jtop_offset;
        uint32_t total = 0;
        while (jt) {
            fseek(in->fp, (long)jt + 12, SEEK_SET);
            uint32_t mx, cur; uint64_t nxt;
            fread(&mx, 4, 1, in->fp); fread(&cur, 4, 1, in->fp); fread(&nxt, 8, 1, in->fp);
            for (uint32_t i = 0; i < cur; i++) {
                uint64_t jo; uint32_t ff, fc;
                fread(&jo, 8, 1, in->fp); fread(&ff, 4, 1, in->fp); fread(&fc, 4, 1, in->fp);
                total += fc;
            }
            jt = nxt;
        }
        if (total > in->total_frames) in->total_frames = total;
    }
    last = in->total_frames - 1;

    /* Find last voxel frame */
    for (int fi = last; fi >= 0; fi--) {
        SFA_L2Entry e;
        if (sfa_find_frame(in, fi, &e) == 0 && e.frame_type == 0) {
            last = fi;
            break;
        }
    }

    void *buf = malloc(in->frame_bytes);
    printf("Reading frame %d from %s\n", last, argv[1]);
    sfa_read_frame(in, last, buf);

    /* Extract fields */
    double *phi[3], *theta[3], *vel[3], *tvel[3];
    for (int a = 0; a < 3; a++) {
        phi[a] = (double*)calloc(N3, sizeof(double));
        theta[a] = (double*)calloc(N3, sizeof(double));
        vel[a] = (double*)calloc(N3, sizeof(double));
        tvel[a] = (double*)calloc(N3, sizeof(double));
    }

    uint64_t off = 0;
    for (uint32_t col = 0; col < in->n_columns && col < 12; col++) {
        int dtype = in->columns[col].dtype;
        int sem = in->columns[col].semantic;
        int comp = in->columns[col].component;
        int es = sfa_dtype_size[dtype];
        uint8_t *src = (uint8_t*)buf + off;
        double *target = NULL;
        if (sem == 0 && comp < 3) target = phi[comp];
        else if (sem == 1 && comp < 3) target = theta[comp];
        else if (sem == 2 && comp < 3) target = vel[comp];
        else if (sem == 2 && comp >= 3 && comp < 6) target = tvel[comp-3];
        if (target) {
            for (long i = 0; i < N3; i++) {
                if (dtype == 0) { /* f16 */
                    uint16_t h; memcpy(&h, src+i*2, 2);
                    target[i] = sfa_f16_to_f32(h);
                } else if (dtype == 1) { /* f32 */
                    float f; memcpy(&f, src+i*4, 4);
                    target[i] = f;
                } else target[i] = ((double*)src)[i];
            }
        }
        off += (uint64_t)N3 * es;
    }
    free(buf); sfa_close(in);

    /* Stats before stripping */
    double phi_max_before = 0, P_int_before = 0;
    for (long i = 0; i < N3; i++) {
        double ps = phi[0][i]*phi[0][i]+phi[1][i]*phi[1][i]+phi[2][i]*phi[2][i];
        if (sqrt(ps) > phi_max_before) phi_max_before = sqrt(ps);
        P_int_before += fabs(phi[0][i]*phi[1][i]*phi[2][i]);
    }
    P_int_before *= dx*dx*dx;
    printf("Before: phi_max=%.4f P_int=%.2f\n", phi_max_before, P_int_before);

    /* Subtract background */
    printf("Stripping background: A_bg=%.3f m=%.3f k_bg=%.4f omega_bg=%.4f\n",
           A_bg, m, k_bg, omega_bg);

    for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
    for (int k = 0; k < N; k++) {
        long idx = (long)i*NN + j*N + k;
        double z = -L + k*dx;
        for (int a = 0; a < 3; a++) {
            phi[a][idx] -= A_bg * cos(k_bg*z + delta[a]);
            vel[a][idx] -= omega_bg * A_bg * sin(k_bg*z + delta[a]);
        }
    }

    /* Stats after stripping */
    double phi_max_after = 0, P_int_after = 0;
    for (long i = 0; i < N3; i++) {
        double ps = phi[0][i]*phi[0][i]+phi[1][i]*phi[1][i]+phi[2][i]*phi[2][i];
        if (sqrt(ps) > phi_max_after) phi_max_after = sqrt(ps);
        P_int_after += fabs(phi[0][i]*phi[1][i]*phi[2][i]);
    }
    P_int_after *= dx*dx*dx;
    printf("After:  phi_max=%.4f P_int=%.2f\n", phi_max_after, P_int_after);

    /* Write output */
    SFA *out = sfa_create(argv[2], N, N, N, L, L, L, 1.0);
    out->flags = 4; /* BSS */
    sfa_add_column(out, "phi_x", 1, 0, 0);
    sfa_add_column(out, "phi_y", 1, 0, 1);
    sfa_add_column(out, "phi_z", 1, 0, 2);
    sfa_add_column(out, "theta_x", 1, 1, 0);
    sfa_add_column(out, "theta_y", 1, 1, 1);
    sfa_add_column(out, "theta_z", 1, 1, 2);
    sfa_add_column(out, "phi_vx", 1, 2, 0);
    sfa_add_column(out, "phi_vy", 1, 2, 1);
    sfa_add_column(out, "phi_vz", 1, 2, 2);
    sfa_add_column(out, "theta_vx", 1, 2, 3);
    sfa_add_column(out, "theta_vy", 1, 2, 4);
    sfa_add_column(out, "theta_vz", 1, 2, 5);
    sfa_finalize_header(out);

    void *cols[12];
    float *cbufs[12];
    for (int c = 0; c < 12; c++) cbufs[c] = (float*)malloc(N3 * sizeof(float));
    for (long i = 0; i < N3; i++) {
        for (int a = 0; a < 3; a++) {
            cbufs[a][i] = (float)phi[a][i];
            cbufs[3+a][i] = (float)theta[a][i];
            cbufs[6+a][i] = (float)vel[a][i];
            cbufs[9+a][i] = (float)tvel[a][i];
        }
    }
    for (int c = 0; c < 12; c++) cols[c] = cbufs[c];
    sfa_write_frame(out, 0.0, cols);
    sfa_close(out);

    printf("Wrote: %s (background stripped)\n", argv[2]);

    for (int a = 0; a < 3; a++) { free(phi[a]); free(theta[a]); free(vel[a]); free(tvel[a]); }
    for (int c = 0; c < 12; c++) free(cbufs[c]);
    return 0;
}
