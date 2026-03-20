/*  theta_radial.c — Radial decay of cylindrical theta components
 *
 *  Opens an SFA archive, reads a frame, decomposes theta into
 *  cylindrical (r, phi, z) components around the braid axis,
 *  and averages in radial shells.
 *
 *  Build: gcc -O2 -o theta_radial src/theta_radial.c -lzstd -lm
 *  Usage: ./theta_radial [-f frame] [-signed] <archive.sfa>
 *         -signed: output signed theta_phi average (not squared)
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define SFA_IMPLEMENTATION
#include "../../viewer/sfa.h"

#define MAX_SHELLS 200

int main(int argc, char **argv) {
    int frame_idx = 200;
    int signed_mode = 0;
    char *sfapath = NULL;

    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "-f") && i+1 < argc) frame_idx = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-signed")) signed_mode = 1;
        else sfapath = argv[i];
    }
    if (!sfapath) {
        fprintf(stderr, "Usage: %s [-f frame] [-signed] <archive.sfa>\n", argv[0]);
        return 1;
    }

    SFA *s = sfa_open(sfapath);
    if (!s) { fprintf(stderr, "Cannot open %s\n", sfapath); return 1; }

    printf("Archive: %s\n", sfapath);
    printf("Grid: %u x %u x %u, L=(%.1f, %.1f, %.1f)\n",
           s->Nx, s->Ny, s->Nz, s->Lx, s->Ly, s->Lz);
    printf("Frames: %u, columns: %u\n", s->total_frames, s->n_columns);

    if ((uint32_t)frame_idx >= s->total_frames) {
        fprintf(stderr, "Frame %d out of range (max %u)\n", frame_idx, s->total_frames-1);
        sfa_close(s);
        return 1;
    }

    uint64_t N_total = (uint64_t)s->Nx * s->Ny * s->Nz;
    void *buf = malloc(s->frame_bytes);
    if (!buf) { fprintf(stderr, "malloc failed\n"); return 1; }

    printf("Reading frame %d...\n", frame_idx);
    if (sfa_read_frame(s, frame_idx, buf) < 0) {
        fprintf(stderr, "Failed to read frame %d\n", frame_idx);
        free(buf); sfa_close(s); return 1;
    }
    double t = sfa_frame_time(s, frame_idx);
    printf("Frame %d: t = %.2f\n", frame_idx, t);

    /* Column layout: phi_x, phi_y, phi_z, theta_x, theta_y, theta_z */
    double *data = (double *)buf;
    double *phi_x   = data + 0 * N_total;
    double *phi_y   = data + 1 * N_total;
    double *phi_z   = data + 2 * N_total;
    double *theta_x = data + 3 * N_total;
    double *theta_y = data + 4 * N_total;
    double *theta_z = data + 5 * N_total;

    int Nx = s->Nx, Ny = s->Ny, Nz = s->Nz;
    double Lx = s->Lx, Ly = s->Ly, Lz = s->Lz;
    double dx = 2.0 * Lx / (Nx - 1);
    double dy = 2.0 * Ly / (Ny - 1);
    double dz = 2.0 * Lz / (Nz - 1);

    /* Step 1: Find braid center (energy-weighted centroid, phi^2 > 5*avg) */
    double avg_phi2 = 0;
    for (uint64_t idx = 0; idx < N_total; idx++) {
        avg_phi2 += phi_x[idx]*phi_x[idx] + phi_y[idx]*phi_y[idx] + phi_z[idx]*phi_z[idx];
    }
    avg_phi2 /= N_total;
    printf("Average phi^2 = %.6e\n", avg_phi2);

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
                wx += x * p2;
                wy += y * p2;
                wt += p2;
            }
        }
    }
    double cx = wx / wt;
    double cy = wy / wt;
    printf("Braid center: (%.3f, %.3f)\n", cx, cy);

    /* Step 2: Accumulate in radial shells */
    double dr = 0.5;
    int n_shells = (int)(20.0 / dr) + 1;
    if (n_shells > MAX_SHELLS) n_shells = MAX_SHELLS;

    double *sum_tphi2  = calloc(n_shells, sizeof(double));  /* theta_phi^2 */
    double *sum_tr2    = calloc(n_shells, sizeof(double));  /* theta_r^2 */
    double *sum_tz2    = calloc(n_shells, sizeof(double));  /* theta_z^2 */
    double *sum_ttot2  = calloc(n_shells, sizeof(double));  /* total theta^2 */
    double *sum_phi2   = calloc(n_shells, sizeof(double));  /* phi^2 */
    double *sum_absP   = calloc(n_shells, sizeof(double));  /* |P| = |phi_x*phi_y*phi_z| */
    double *sum_tphi   = calloc(n_shells, sizeof(double));  /* signed theta_phi */
    long   *count      = calloc(n_shells, sizeof(long));

    for (int i = 0; i < Nx; i++) {
        double x = -Lx + i * dx;
        for (int j = 0; j < Ny; j++) {
            double y = -Ly + j * dy;
            double xc = x - cx;
            double yc = y - cy;
            double r_perp = sqrt(xc * xc + yc * yc);
            int shell = (int)(r_perp / dr);
            if (shell >= n_shells) continue;

            double angle = atan2(yc, xc);
            double ca = cos(angle), sa = sin(angle);

            for (int k = 0; k < Nz; k++) {
                long idx = (long)i * (Ny * Nz) + j * Nz + k;

                /* Cylindrical decomposition of theta */
                double tx = theta_x[idx], ty = theta_y[idx], tz = theta_z[idx];
                double t_r   =  tx * ca + ty * sa;    /* radial */
                double t_phi = -tx * sa + ty * ca;    /* azimuthal (circular) */
                double t_z   =  tz;                   /* axial */

                sum_tphi2[shell] += t_phi * t_phi;
                sum_tr2[shell]   += t_r * t_r;
                sum_tz2[shell]   += t_z * t_z;
                sum_ttot2[shell] += tx*tx + ty*ty + tz*tz;
                sum_tphi[shell]  += t_phi;

                double px = phi_x[idx], py = phi_y[idx], pz = phi_z[idx];
                sum_phi2[shell]  += px*px + py*py + pz*pz;
                sum_absP[shell]  += fabs(px * py * pz);
                count[shell]++;
            }
        }
    }

    /* Output */
    const char *outpath = signed_mode ? "data/theta_radial_signed.tsv" : "data/theta_radial.tsv";
    FILE *fp = fopen(outpath, "w");
    if (signed_mode) {
        fprintf(fp, "r\ttheta_phi\ttheta_phi2\ttheta_r2\ttheta_z2\ttheta_total2\tphi2\tabsP\n");
    } else {
        fprintf(fp, "r\ttheta_phi2\ttheta_r2\ttheta_z2\ttheta_total2\tphi2\tabsP\n");
    }

    printf("\n%5s %12s %12s %12s %12s %12s %12s %8s\n",
           "r", "theta_phi2", "theta_r2", "theta_z2", "theta_tot2", "phi2", "|P|", "count");
    printf("------------------------------------------------------------------------"
           "--------------------\n");

    for (int s_idx = 0; s_idx < n_shells; s_idx++) {
        if (count[s_idx] == 0) continue;
        double r = (s_idx + 0.5) * dr;
        double n = (double)count[s_idx];
        double tp2 = sum_tphi2[s_idx] / n;
        double tr2 = sum_tr2[s_idx] / n;
        double tz2 = sum_tz2[s_idx] / n;
        double tt2 = sum_ttot2[s_idx] / n;
        double p2  = sum_phi2[s_idx] / n;
        double aP  = sum_absP[s_idx] / n;
        double tp  = sum_tphi[s_idx] / n;

        if (signed_mode) {
            fprintf(fp, "%.2f\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\n",
                    r, tp, tp2, tr2, tz2, tt2, p2, aP);
        } else {
            fprintf(fp, "%.2f\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\n",
                    r, tp2, tr2, tz2, tt2, p2, aP);
        }
        printf("%5.2f %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e %8ld\n",
               r, tp2, tr2, tz2, tt2, p2, aP, count[s_idx]);
    }

    fclose(fp);
    printf("\nWritten: %s\n", outpath);

    free(sum_tphi2); free(sum_tr2); free(sum_tz2); free(sum_ttot2);
    free(sum_phi2); free(sum_absP); free(sum_tphi); free(count);
    free(buf);
    sfa_close(s);
    return 0;
}
