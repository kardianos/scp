/* gen_qball_pair.c — stamp TWO Q-balls into one SFA seed for the v66
 * complexified kernel (complex_phi=1, 12 fields + 12 velocities = 24 columns).
 *
 * Ansatz (v66/THEORY.md §1,§3): a ball at rest with radial profile f(r),
 * internal phase delta and rotation frequency omega has, for all a = 0,1,2:
 *     u_a    =  f(r) cos(delta)          (phi real part)
 *     v_a    =  f(r) sin(delta)          (phi imaginary part)
 *     udot_a = -omega f(r) sin(delta)
 *     vdot_a = +omega f(r) cos(delta)
 * Theta sector (and velocities) zero. Negative omega = counter-rotating
 * (anti-ball, charge -Q); the sign flows through the velocity formulas.
 * Two balls are superposed linearly.
 *
 * Output: ONE SFA frame, f32, 24 columns named/typed exactly as the kernel
 * registers them (scp_sim.c sfa_add_column block / v66 SPEC §7.1), so
 * init=sfa loads every slot.
 *
 * Build:
 *   gcc -O3 -march=native -o bin/gen_qball_pair sfa/seed/gen_qball_pair.c -lzstd -lm
 *
 * Usage:
 *   gen_qball_pair N L profile1 omega1 delta1 x1 y1 z1 \
 *                      profile2 omega2 delta2 x2 y2 z2 out.sfa
 *   (delta in radians; profile = radial_qball output: '#' headers then "r f")
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define SFA_IMPLEMENTATION
#include "../format/sfa.h"

#define NCOLS 24

typedef struct {
    double *r, *f;
    size_t  n;
} Profile;

static void load_profile(const char *path, Profile *p) {
    FILE *fp = fopen(path, "r");
    if (!fp) { fprintf(stderr, "FATAL: cannot open profile '%s'\n", path); exit(1); }
    size_t cap = 1024;
    p->r = (double *)malloc(cap * sizeof(double));
    p->f = (double *)malloc(cap * sizeof(double));
    if (!p->r || !p->f) { fprintf(stderr, "FATAL: alloc\n"); exit(1); }
    p->n = 0;
    char line[512];
    while (fgets(line, sizeof(line), fp)) {
        const char *s = line;
        while (*s == ' ' || *s == '\t') s++;
        if (*s == '#' || *s == '\n' || *s == '\0' || *s == '\r') continue;
        double rv, fv;
        if (sscanf(s, "%lf %lf", &rv, &fv) != 2) {
            fprintf(stderr, "FATAL: parse error in '%s': %s", path, line);
            exit(1);
        }
        if (p->n > 0 && rv <= p->r[p->n - 1]) {
            fprintf(stderr, "FATAL: r not strictly increasing in '%s' at r=%g\n",
                    path, rv);
            exit(1);
        }
        if (p->n == cap) {
            cap *= 2;
            p->r = (double *)realloc(p->r, cap * sizeof(double));
            p->f = (double *)realloc(p->f, cap * sizeof(double));
            if (!p->r || !p->f) { fprintf(stderr, "FATAL: realloc\n"); exit(1); }
        }
        p->r[p->n] = rv;
        p->f[p->n] = fv;
        p->n++;
    }
    fclose(fp);
    if (p->n < 2) {
        fprintf(stderr, "FATAL: profile '%s' has <2 data points\n", path);
        exit(1);
    }
    printf("  profile '%s': %zu points, r=[%.3f,%.3f], f0=%.6f, f_last=%.3g\n",
           path, p->n, p->r[0], p->r[p->n - 1], p->f[0], p->f[p->n - 1]);
    if (fabs(p->f[p->n - 1]) > 1e-8)
        printf("  WARNING: profile '%s' truncated while still sizable "
               "(f_last=%.3g)\n", path, p->f[p->n - 1]);
}

/* piecewise-linear interpolation; flat below first r, zero above last */
static double interp(const Profile *p, double r) {
    if (r <= p->r[0]) return p->f[0];
    if (r >= p->r[p->n - 1]) return 0.0;
    size_t lo = 0, hi = p->n - 1;
    while (hi - lo > 1) {
        size_t mid = (lo + hi) / 2;
        if (p->r[mid] <= r) lo = mid; else hi = mid;
    }
    double t = (r - p->r[lo]) / (p->r[hi] - p->r[lo]);
    return p->f[lo] + t * (p->f[hi] - p->f[lo]);
}

int main(int argc, char **argv) {
    if (argc != 16) {
        fprintf(stderr,
            "Usage: %s N L profile1 omega1 delta1 x1 y1 z1 "
            "profile2 omega2 delta2 x2 y2 z2 out.sfa\n"
            "  delta in radians; negative omega = counter-rotating (charge -Q)\n",
            argv[0]);
        return 1;
    }
    int    N      = atoi(argv[1]);
    double L      = atof(argv[2]);
    const char *ppath[2] = { argv[3], argv[9] };
    double omega[2] = { atof(argv[4]),  atof(argv[10]) };
    double delta[2] = { atof(argv[5]),  atof(argv[11]) };
    double cx[2]    = { atof(argv[6]),  atof(argv[12]) };
    double cy[2]    = { atof(argv[7]),  atof(argv[13]) };
    double cz[2]    = { atof(argv[8]),  atof(argv[14]) };
    const char *outpath = argv[15];

    if (N < 2 || L <= 0.0) {
        fprintf(stderr, "FATAL: bad N=%d or L=%g\n", N, L);
        return 1;
    }

    printf("gen_qball_pair: N=%d L=%.4f -> %s\n", N, L, outpath);
    Profile prof[2];
    for (int b = 0; b < 2; b++) {
        load_profile(ppath[b], &prof[b]);
        printf("  ball %d: omega=%.6f delta=%.6f rad, center=(%.3f,%.3f,%.3f)\n",
               b + 1, omega[b], delta[b], cx[b], cy[b], cz[b]);
    }

    long N3 = (long)N * N * N;
    double dx = 2.0 * L / (N - 1);   /* grid: x = -L + i*dx (scp_init.h) */

    /* column order MUST match kernel registration (scp_sim.c / SPEC §7.1):
     *  0-2  phi_x/y/z       (u_a)        12-14 phiim_x/y/z     (v_a)
     *  3-5  theta_x/y/z     (tu_a)       15-17 thetaim_x/y/z   (tv_a)
     *  6-8  phi_vx/vy/vz    (udot_a)     18-20 phiim_vx/vy/vz  (vdot_a)
     *  9-11 theta_vx/vy/vz  (tudot_a)    21-23 thetaim_vx/...  (tvdot_a)  */
    float *cols[NCOLS];
    for (int c = 0; c < NCOLS; c++) {
        cols[c] = (float *)calloc((size_t)N3, sizeof(float));
        if (!cols[c]) { fprintf(stderr, "FATAL: alloc column %d\n", c); return 1; }
    }

    double cu[2], su[2];
    for (int b = 0; b < 2; b++) { cu[b] = cos(delta[b]); su[b] = sin(delta[b]); }

    for (int i = 0; i < N; i++) {
        double x = -L + i * dx;
        for (int j = 0; j < N; j++) {
            double y = -L + j * dx;
            for (int k = 0; k < N; k++) {
                double z = -L + k * dx;
                long idx = (long)i * N * N + (long)j * N + k;
                double u = 0.0, v = 0.0, ud = 0.0, vd = 0.0;
                for (int b = 0; b < 2; b++) {
                    double r = sqrt((x - cx[b]) * (x - cx[b]) +
                                    (y - cy[b]) * (y - cy[b]) +
                                    (z - cz[b]) * (z - cz[b]));
                    double f = interp(&prof[b], r);
                    u  += f * cu[b];
                    v  += f * su[b];
                    ud += -omega[b] * f * su[b];
                    vd += +omega[b] * f * cu[b];
                }
                for (int a = 0; a < 3; a++) {
                    cols[0 + a][idx]  = (float)u;   /* phi   (u_a)    */
                    cols[6 + a][idx]  = (float)ud;  /* phi_v (udot_a) */
                    cols[12 + a][idx] = (float)v;   /* phiim (v_a)    */
                    cols[18 + a][idx] = (float)vd;  /* phiim_v (vdot) */
                    /* theta sector (3-5, 9-11, 15-17, 21-23) stays zero */
                }
            }
        }
    }

    double dt = 0.025 * dx;  /* standard dt factor (metadata only) */
    SFA *sfa = sfa_create(outpath, N, N, N, L, L, L, dt);
    if (!sfa) { fprintf(stderr, "FATAL: cannot create '%s'\n", outpath); return 1; }

    /* metadata */
    char bN[32], bL[32], bo1[32], bo2[32], bd1[32], bd2[32], bp1[96], bp2[96];
    snprintf(bN, sizeof(bN), "%d", N);
    snprintf(bL, sizeof(bL), "%.6f", L);
    snprintf(bo1, sizeof(bo1), "%.6f", omega[0]);
    snprintf(bo2, sizeof(bo2), "%.6f", omega[1]);
    snprintf(bd1, sizeof(bd1), "%.6f", delta[0]);
    snprintf(bd2, sizeof(bd2), "%.6f", delta[1]);
    snprintf(bp1, sizeof(bp1), "%.6f,%.6f,%.6f", cx[0], cy[0], cz[0]);
    snprintf(bp2, sizeof(bp2), "%.6f,%.6f,%.6f", cx[1], cy[1], cz[1]);
    const char *keys[] = { "generator", "N", "L", "omega_1", "delta_1",
                           "position_1", "omega_2", "delta_2", "position_2" };
    const char *vals[] = { "gen_qball_pair", bN, bL, bo1, bd1, bp1, bo2, bd2, bp2 };
    sfa_add_kvmd(sfa, 0, 0xFFFFFFFF, 0xFFFFFFFF, keys, vals, 9);

    /* exact names/semantics/components as the kernel registers (24 cols) */
    sfa_add_column(sfa, "phi_x",      SFA_F32, SFA_POSITION, 0);
    sfa_add_column(sfa, "phi_y",      SFA_F32, SFA_POSITION, 1);
    sfa_add_column(sfa, "phi_z",      SFA_F32, SFA_POSITION, 2);
    sfa_add_column(sfa, "theta_x",    SFA_F32, SFA_ANGLE,    0);
    sfa_add_column(sfa, "theta_y",    SFA_F32, SFA_ANGLE,    1);
    sfa_add_column(sfa, "theta_z",    SFA_F32, SFA_ANGLE,    2);
    sfa_add_column(sfa, "phi_vx",     SFA_F32, SFA_VELOCITY, 0);
    sfa_add_column(sfa, "phi_vy",     SFA_F32, SFA_VELOCITY, 1);
    sfa_add_column(sfa, "phi_vz",     SFA_F32, SFA_VELOCITY, 2);
    sfa_add_column(sfa, "theta_vx",   SFA_F32, SFA_VELOCITY, 3);
    sfa_add_column(sfa, "theta_vy",   SFA_F32, SFA_VELOCITY, 4);
    sfa_add_column(sfa, "theta_vz",   SFA_F32, SFA_VELOCITY, 5);
    sfa_add_column(sfa, "phiim_x",    SFA_F32, SFA_POSITION, 3);
    sfa_add_column(sfa, "phiim_y",    SFA_F32, SFA_POSITION, 4);
    sfa_add_column(sfa, "phiim_z",    SFA_F32, SFA_POSITION, 5);
    sfa_add_column(sfa, "thetaim_x",  SFA_F32, SFA_ANGLE,    3);
    sfa_add_column(sfa, "thetaim_y",  SFA_F32, SFA_ANGLE,    4);
    sfa_add_column(sfa, "thetaim_z",  SFA_F32, SFA_ANGLE,    5);
    sfa_add_column(sfa, "phiim_vx",   SFA_F32, SFA_VELOCITY, 6);
    sfa_add_column(sfa, "phiim_vy",   SFA_F32, SFA_VELOCITY, 7);
    sfa_add_column(sfa, "phiim_vz",   SFA_F32, SFA_VELOCITY, 8);
    sfa_add_column(sfa, "thetaim_vx", SFA_F32, SFA_VELOCITY, 9);
    sfa_add_column(sfa, "thetaim_vy", SFA_F32, SFA_VELOCITY, 10);
    sfa_add_column(sfa, "thetaim_vz", SFA_F32, SFA_VELOCITY, 11);
    sfa_finalize_header(sfa);

    void *frame_cols[NCOLS];
    for (int c = 0; c < NCOLS; c++) frame_cols[c] = cols[c];
    sfa_write_frame(sfa, 0.0, frame_cols);
    sfa_close(sfa);

    printf("gen_qball_pair: wrote %s (%d^3, f32, 24 columns, 1 frame)\n",
           outpath, N);

    for (int c = 0; c < NCOLS; c++) free(cols[c]);
    for (int b = 0; b < 2; b++) { free(prof[b].r); free(prof[b].f); }
    return 0;
}
