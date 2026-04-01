/*  gen_collision.c — Two-proton collision seed
 *
 *  Places two copies of a proton template at ±D/2 along x-axis,
 *  with opposite x-velocities (Galilean boost of the velocity field).
 *  Background generated analytically.
 *
 *  Build: gcc -O3 -o gen_collision gen_collision.c -lzstd -lm
 *  Usage: ./gen_collision -template /space/scp/v44/uud_seed.sfa
 *         -o /space/scp/v51/seed_collision.sfa
 *         [-N 384] [-L 50] [-D 25] [-v 0.1]
 */

#define SFA_IMPLEMENTATION
#include "../../sfa/format/sfa.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define NFIELDS 3
#define PI 3.14159265358979323846

static double f16_to_f64(uint16_t h) {
    uint16_t sign = h & 0x8000; int exp = (h >> 10) & 0x1F; uint16_t mant = h & 0x3FF;
    if (exp == 0) return sign ? -0.0 : 0.0;
    if (exp == 31) return sign ? -INFINITY : INFINITY;
    float f; uint32_t x = ((uint32_t)sign << 16) | ((uint32_t)(exp-15+127) << 23) | ((uint32_t)mant << 13);
    memcpy(&f, &x, 4); return (double)f;
}

int main(int argc, char **argv) {
    int N = 384;
    double L = 50.0, D = 25.0, v_kick = 0.1;
    double A_bg = 0.1, m2 = 2.25;
    double dt_factor = 0.025;
    char template_path[512] = "";
    char outpath[512] = "seed_collision.sfa";

    for (int i = 1; i < argc - 1; i += 2) {
        const char *k = argv[i], *v = argv[i+1];
        if      (!strcmp(k,"-N"))        N = atoi(v);
        else if (!strcmp(k,"-L"))        L = atof(v);
        else if (!strcmp(k,"-D"))        D = atof(v);
        else if (!strcmp(k,"-v"))        v_kick = atof(v);
        else if (!strcmp(k,"-A_bg"))     A_bg = atof(v);
        else if (!strcmp(k,"-template")) strncpy(template_path, v, 511);
        else if (!strcmp(k,"-o"))        strncpy(outpath, v, 511);
    }

    if (!template_path[0]) {
        fprintf(stderr, "Usage: %s -template proton.sfa -o output.sfa\n", argv[0]);
        return 1;
    }

    long N3 = (long)N * N * N;
    double dx = 2.0 * L / (N - 1);
    double dt = dt_factor * dx;
    int NN = N * N;

    fprintf(stderr, "gen_collision: N=%d L=%.1f D=%.1f v_kick=%.2f\n", N, L, D, v_kick);
    fprintf(stderr, "  Template: %s\n", template_path);

    /* Allocate main grid */
    double *phi[3], *theta[3], *phi_vel[3], *theta_vel[3];
    for (int a = 0; a < 3; a++) {
        phi[a]       = (double*)calloc(N3, sizeof(double));
        phi_vel[a]   = (double*)calloc(N3, sizeof(double));
        theta[a]     = (double*)calloc(N3, sizeof(double));
        theta_vel[a] = (double*)calloc(N3, sizeof(double));
    }

    /* Generate background */
    double k_bg = PI / L, omega_bg = sqrt(k_bg * k_bg + m2);
    for (long idx = 0; idx < N3; idx++) {
        int kk = (int)(idx % N);
        double z = -L + kk * dx;
        for (int a = 0; a < NFIELDS; a++) {
            double ph = k_bg * z + 2.0 * PI * a / 3.0;
            phi[a][idx] = A_bg * cos(ph);
            phi_vel[a][idx] = omega_bg * A_bg * sin(ph);
        }
    }

    /* Load template */
    SFA *tmpl = sfa_open(template_path);
    if (!tmpl) { fprintf(stderr, "FATAL: cannot open '%s'\n", template_path); return 1; }
    int TN = tmpl->Nx; double TL = tmpl->Lx, Tdx = 2.0 * TL / (TN - 1);
    long TN3 = (long)TN * TN * TN; int TNN = TN * TN;
    fprintf(stderr, "  Template: %d^3, L=%.1f\n", TN, TL);

    int frame = tmpl->total_frames - 1;
    void *buf = malloc(tmpl->frame_bytes);
    sfa_read_frame(tmpl, frame, buf);

    double *tphi[3], *tvel[3], *ttheta[3], *ttvel[3];
    for (int a = 0; a < 3; a++) {
        tphi[a] = (double*)calloc(TN3, sizeof(double));
        tvel[a] = (double*)calloc(TN3, sizeof(double));
        ttheta[a] = (double*)calloc(TN3, sizeof(double));
        ttvel[a] = (double*)calloc(TN3, sizeof(double));
    }
    uint64_t off = 0;
    for (int col = 0; col < tmpl->n_columns; col++) {
        int dtype = tmpl->columns[col].dtype, sem = tmpl->columns[col].semantic;
        int comp = tmpl->columns[col].component;
        int es = sfa_dtype_size[dtype]; uint8_t *src = (uint8_t*)buf + off;
        double *target = NULL;
        if (sem == SFA_POSITION && comp < 3) target = tphi[comp];
        else if (sem == SFA_ANGLE && comp < 3) target = ttheta[comp];
        else if (sem == SFA_VELOCITY && comp < 3) target = tvel[comp];
        else if (sem == SFA_VELOCITY && comp >= 3 && comp < 6) target = ttvel[comp - 3];
        if (target) {
            if (dtype == SFA_F64) for (long i = 0; i < TN3; i++) target[i] = ((double*)src)[i];
            else if (dtype == SFA_F32) for (long i = 0; i < TN3; i++) target[i] = (double)((float*)src)[i];
            else if (dtype == SFA_F16) for (long i = 0; i < TN3; i++) target[i] = f16_to_f64(((uint16_t*)src)[i]);
        }
        off += (uint64_t)TN3 * es;
    }
    free(buf); sfa_close(tmpl);

    /* Stamp two copies: proton A at x = -D/2, proton B at x = +D/2 */
    int half = TN / 2;
    double positions[2] = {-D / 2.0, +D / 2.0};
    double velocities[2] = {+v_kick, -v_kick};  /* moving toward each other */

    for (int p = 0; p < 2; p++) {
        /* Center of this proton in grid coords */
        int ci = (int)((positions[p] + L) / dx);
        int cj = N / 2, ck = N / 2;

        fprintf(stderr, "  Proton %d: x=%.1f (grid i=%d), v_x=%+.2f\n",
                p + 1, positions[p], ci, velocities[p]);

        int placed = 0;
        for (int ti = 0; ti < TN; ti++) {
            int gi = ci + ti - half; if (gi < 0 || gi >= N) continue;
            for (int tj = 0; tj < TN; tj++) {
                int gj = cj + tj - half; if (gj < 0 || gj >= N) continue;
                for (int tk = 0; tk < TN; tk++) {
                    int gk = ck + tk - half; if (gk < 0 || gk >= N) continue;
                    long tidx = (long)ti * TNN + tj * TN + tk;
                    long gidx = (long)gi * NN + gj * N + gk;

                    /* Template z for background subtraction */
                    double tz = -TL + tk * Tdx;
                    for (int a = 0; a < NFIELDS; a++) {
                        double ph_bg = k_bg * tz + 2.0 * PI * a / 3.0;
                        double bg_phi = A_bg * cos(ph_bg);
                        double bg_vel = omega_bg * A_bg * sin(ph_bg);

                        /* Perturbation = template - background */
                        phi[a][gidx] += tphi[a][tidx] - bg_phi;
                        phi_vel[a][gidx] += tvel[a][tidx] - bg_vel;
                        theta[a][gidx] += ttheta[a][tidx];
                        theta_vel[a][gidx] += ttvel[a][tidx];
                    }
                    /* Galilean velocity boost in x-direction:
                     * Add v_kick to phi_vel[0] (x-velocity of phi_x)
                     * This is approximate — proper boost would be Lorentz, but v << c */
                    phi_vel[0][gidx] += velocities[p] * phi[0][gidx];
                    phi_vel[1][gidx] += velocities[p] * phi[1][gidx];
                    phi_vel[2][gidx] += velocities[p] * phi[2][gidx];
                    placed++;
                }
            }
        }
        fprintf(stderr, "  Placed %d voxels\n", placed);
    }

    for (int a = 0; a < 3; a++) {
        free(tphi[a]); free(tvel[a]); free(ttheta[a]); free(ttvel[a]);
    }

    /* Write SFA */
    SFA *sfa = sfa_create(outpath, N, N, N, L, L, L, dt);
    uint8_t sfa_dtype = SFA_F64;
    sfa_add_column(sfa, "phi_x",    sfa_dtype, SFA_POSITION, 0);
    sfa_add_column(sfa, "phi_y",    sfa_dtype, SFA_POSITION, 1);
    sfa_add_column(sfa, "phi_z",    sfa_dtype, SFA_POSITION, 2);
    sfa_add_column(sfa, "theta_x",  sfa_dtype, SFA_ANGLE,    0);
    sfa_add_column(sfa, "theta_y",  sfa_dtype, SFA_ANGLE,    1);
    sfa_add_column(sfa, "theta_z",  sfa_dtype, SFA_ANGLE,    2);
    sfa_add_column(sfa, "phi_vx",   sfa_dtype, SFA_VELOCITY, 0);
    sfa_add_column(sfa, "phi_vy",   sfa_dtype, SFA_VELOCITY, 1);
    sfa_add_column(sfa, "phi_vz",   sfa_dtype, SFA_VELOCITY, 2);
    sfa_add_column(sfa, "theta_vx", sfa_dtype, SFA_VELOCITY, 3);
    sfa_add_column(sfa, "theta_vy", sfa_dtype, SFA_VELOCITY, 4);
    sfa_add_column(sfa, "theta_vz", sfa_dtype, SFA_VELOCITY, 5);
    sfa_finalize_header(sfa);

    void *cols[12] = {
        phi[0], phi[1], phi[2], theta[0], theta[1], theta[2],
        phi_vel[0], phi_vel[1], phi_vel[2],
        theta_vel[0], theta_vel[1], theta_vel[2]
    };
    sfa_write_frame(sfa, 0.0, cols);
    sfa_close(sfa);

    fprintf(stderr, "  Wrote: %s (%.0f MB)\n", outpath, (double)N3 * 12 * 8 / 1e6);
    for (int a = 0; a < 3; a++) {
        free(phi[a]); free(phi_vel[a]); free(theta[a]); free(theta_vel[a]);
    }
    return 0;
}
