/*  gen_qball_quark.c — per-COMPONENT Q-ball lump seeding ("quark" seeds)
 *
 *  The complexified theory binds only where all three component fields
 *  Phi_a overlap (product potential s = prod_a |Phi_a|^2). This generator
 *  stamps each component as its OWN displaced lump:
 *
 *      Phi_a(x, t=0) = f(|x - c_a|) e^{i omega t}   (component a at center c_a)
 *      u_a = f, v_a = 0, udot_a = 0, vdot_a = omega f
 *
 *  enabling single-component ("free quark"), two-component ("meson"), and
 *  three-displaced-component ("proton assembly") experiments. A component
 *  is disabled by mask digit '0' (its fields stay identically zero).
 *
 *  Output: 24-column complex SFA (matter only; gauge E is built by the
 *  kernel's init Gauss projection).
 *
 *  Build: gcc -O3 -fopenmp -o gen_qball_quark gen_qball_quark.c -lzstd -lm
 *  Usage:
 *    gen_qball_quark N L profile omega mask x0 y0 z0 x1 y1 z1 x2 y2 z2 out.sfa
 *    (mask: three digits, e.g. 111, 100, 110; profile: radial 'r f' table,
 *     extra columns ignored; negative omega = opposite charge, all comps)
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define SFA_IMPLEMENTATION
#include "../format/sfa.h"

#define NCOLS 24

typedef struct { double *r, *f; size_t n; } Profile;

static void load_profile(const char *path, Profile *p) {
    FILE *fp = fopen(path, "r");
    if (!fp) { fprintf(stderr, "FATAL: cannot open profile '%s'\n", path); exit(1); }
    size_t cap = 1024;
    p->r = (double *)malloc(cap * sizeof(double));
    p->f = (double *)malloc(cap * sizeof(double));
    p->n = 0;
    char line[512];
    while (fgets(line, sizeof(line), fp)) {
        const char *s = line;
        while (*s == ' ' || *s == '\t') s++;
        if (*s == '#' || *s == '\n' || *s == '\0' || *s == '\r') continue;
        double rv, fv;
        if (sscanf(s, "%lf %lf", &rv, &fv) != 2) continue;
        if (p->n == cap) {
            cap *= 2;
            p->r = (double *)realloc(p->r, cap * sizeof(double));
            p->f = (double *)realloc(p->f, cap * sizeof(double));
        }
        p->r[p->n] = rv; p->f[p->n] = fv; p->n++;
    }
    fclose(fp);
    if (p->n < 2) { fprintf(stderr, "FATAL: profile too short\n"); exit(1); }
    printf("  profile '%s': %zu points, f0=%.6f\n", path, p->n, p->f[0]);
}

static double interp(const Profile *p, double r) {
    if (r <= p->r[0]) return p->f[0];
    if (r >= p->r[p->n - 1]) return 0.0;
    size_t lo = 0, hi = p->n - 1;
    while (hi - lo > 1) { size_t mid = (lo + hi) / 2; if (p->r[mid] <= r) lo = mid; else hi = mid; }
    double t = (r - p->r[lo]) / (p->r[hi] - p->r[lo]);
    return p->f[lo] + t * (p->f[hi] - p->f[lo]);
}

int main(int argc, char **argv) {
    if (argc != 16) {
        fprintf(stderr,
            "Usage: %s N L profile omega mask x0 y0 z0 x1 y1 z1 x2 y2 z2 out.sfa\n"
            "  mask: 3 digits (111=all comps, 100=only comp0, ...)\n", argv[0]);
        return 1;
    }
    int N = atoi(argv[1]);
    double L = atof(argv[2]);
    const char *ppath = argv[3];
    double omega = atof(argv[4]);
    const char *mask = argv[5];
    double cx[3], cy[3], cz[3];
    for (int a = 0; a < 3; a++) {
        cx[a] = atof(argv[6 + 3*a]);
        cy[a] = atof(argv[7 + 3*a]);
        cz[a] = atof(argv[8 + 3*a]);
    }
    const char *outpath = argv[15];
    int on[3];
    for (int a = 0; a < 3; a++) on[a] = (strlen(mask) > (size_t)a && mask[a] == '1');

    Profile prof;
    load_profile(ppath, &prof);

    long N3 = (long)N * N * N, NN = (long)N * N;
    double dx = 2.0 * L / (N - 1);
    printf("gen_qball_quark: N=%d L=%g mask=%s omega=%g -> %s\n", N, L, mask, omega, outpath);
    for (int a = 0; a < 3; a++)
        printf("  comp %d: %s center=(%.2f,%.2f,%.2f)\n",
               a, on[a] ? "ON " : "off", cx[a], cy[a], cz[a]);

    float **col = (float **)malloc(NCOLS * sizeof(float *));
    for (int c = 0; c < NCOLS; c++) col[c] = (float *)calloc(N3, sizeof(float));
    /* column order mirrors gen_qball_pair / kernel registration:
       0-2 phi_xyz  3-5 theta_xyz  6-8 phi_vxyz  9-11 theta_vxyz
       12-14 phiim_xyz  15-17 thetaim_xyz  18-20 phiim_vxyz  21-23 thetaim_vxyz */
    #pragma omp parallel for
    for (long i = 0; i < N3; i++) {
        int ix = (int)(i / NN), iy = (int)((i / N) % N), iz = (int)(i % N);
        double x = -L + ix * dx, y = -L + iy * dx, z = -L + iz * dx;
        for (int a = 0; a < 3; a++) {
            if (!on[a]) continue;
            double dxa = x - cx[a], dya = y - cy[a], dza = z - cz[a];
            double r = sqrt(dxa*dxa + dya*dya + dza*dza);
            double f = interp(&prof, r);
            if (f == 0.0) continue;
            col[0 + a][i]  += (float)f;                 /* u_a   */
            col[18 + a][i] += (float)(omega * f);       /* vdot_a */
        }
    }

    SFA *sfa = sfa_create(outpath, N, N, N, L, L, L, 0.0);
    if (!sfa) { fprintf(stderr, "FATAL: cannot create %s\n", outpath); return 1; }
    static const char *names[NCOLS] = {
        "phi_x","phi_y","phi_z","theta_x","theta_y","theta_z",
        "phi_vx","phi_vy","phi_vz","theta_vx","theta_vy","theta_vz",
        "phiim_x","phiim_y","phiim_z","thetaim_x","thetaim_y","thetaim_z",
        "phiim_vx","phiim_vy","phiim_vz","thetaim_vx","thetaim_vy","thetaim_vz"};
    static const int sem[NCOLS] = {
        SFA_POSITION,SFA_POSITION,SFA_POSITION, SFA_ANGLE,SFA_ANGLE,SFA_ANGLE,
        SFA_VELOCITY,SFA_VELOCITY,SFA_VELOCITY, SFA_VELOCITY,SFA_VELOCITY,SFA_VELOCITY,
        SFA_POSITION,SFA_POSITION,SFA_POSITION, SFA_ANGLE,SFA_ANGLE,SFA_ANGLE,
        SFA_VELOCITY,SFA_VELOCITY,SFA_VELOCITY, SFA_VELOCITY,SFA_VELOCITY,SFA_VELOCITY};
    static const int comp[NCOLS] = {0,1,2, 0,1,2, 0,1,2, 3,4,5, 3,4,5, 3,4,5, 6,7,8, 9,10,11};
    for (int c = 0; c < NCOLS; c++)
        sfa_add_column(sfa, names[c], SFA_F32, sem[c], comp[c]);
    sfa_finalize_header(sfa);
    void **arrays = (void **)malloc(NCOLS * sizeof(void *));
    for (int c = 0; c < NCOLS; c++) arrays[c] = col[c];
    if (sfa_write_frame(sfa, 0.0, arrays) != 0) {
        fprintf(stderr, "FATAL: write failed\n"); return 1;
    }
    sfa_close(sfa);
    printf("gen_qball_quark: wrote %s (%d^3, f32, 24 columns, 1 frame)\n", outpath, N);
    return 0;
}
