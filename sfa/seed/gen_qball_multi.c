/*  gen_qball_multi.c — stamp N Q-balls into one complex SFA seed (nuclei seeds)
 *
 *  Generalization of gen_qball_pair to arbitrary ball count: each ball gets
 *  its own profile, frequency (negative = opposite charge), phase, and center.
 *
 *      Phi_a += f(|x-c|) e^{i(omega t + delta)}  for all three components a
 *      u_a += f cos(delta), v_a += f sin(delta),
 *      udot_a += -omega f sin(delta), vdot_a += +omega f cos(delta)
 *
 *  Output: 24-column complex SFA (matter only; gauge E built by the kernel's
 *  init Gauss projection).
 *
 *  Build: gcc -O3 -fopenmp -o gen_qball_multi gen_qball_multi.c -lzstd -lm
 *  Usage:
 *    gen_qball_multi N L out.sfa  profile omega delta x y z  [profile omega delta x y z ...]
 *    (one 6-argument group per ball, up to 16 balls)
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define SFA_IMPLEMENTATION
#include "../format/sfa.h"

#define NCOLS    24
#define MAXBALLS 16

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
    if (p->n < 2) { fprintf(stderr, "FATAL: profile '%s' too short\n", path); exit(1); }
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
    if (argc < 10 || (argc - 4) % 6 != 0) {
        fprintf(stderr,
            "Usage: %s N L out.sfa  profile omega delta x y z  [.. more balls ..]\n"
            "  one 6-arg group per ball (max %d); negative omega = opposite charge\n",
            argv[0], MAXBALLS);
        return 1;
    }
    int N = atoi(argv[1]);
    double L = atof(argv[2]);
    const char *outpath = argv[3];
    int nballs = (argc - 4) / 6;
    if (nballs > MAXBALLS) { fprintf(stderr, "FATAL: too many balls\n"); return 1; }

    Profile prof[MAXBALLS];
    double omega[MAXBALLS], delta[MAXBALLS], cx[MAXBALLS], cy[MAXBALLS], cz[MAXBALLS];
    printf("gen_qball_multi: N=%d L=%g, %d balls -> %s\n", N, L, nballs, outpath);
    for (int b = 0; b < nballs; b++) {
        const char **g = (const char **)&argv[4 + 6 * b];
        load_profile(g[0], &prof[b]);
        omega[b] = atof(g[1]); delta[b] = atof(g[2]);
        cx[b] = atof(g[3]); cy[b] = atof(g[4]); cz[b] = atof(g[5]);
        printf("  ball %d: omega=%+.4f delta=%.3f center=(%.2f,%.2f,%.2f) f0=%.4f\n",
               b, omega[b], delta[b], cx[b], cy[b], cz[b], prof[b].f[0]);
    }

    long N3 = (long)N * N * N, NN = (long)N * N;
    double dx = 2.0 * L / (N - 1);

    float **col = (float **)malloc(NCOLS * sizeof(float *));
    for (int c = 0; c < NCOLS; c++) col[c] = (float *)calloc(N3, sizeof(float));
    /* cols: 0-2 phi  3-5 theta  6-8 phi_v  9-11 theta_v
             12-14 phiim  15-17 thetaim  18-20 phiim_v  21-23 thetaim_v */
    #pragma omp parallel for
    for (long i = 0; i < N3; i++) {
        int ix = (int)(i / NN), iy = (int)((i / N) % N), iz = (int)(i % N);
        double x = -L + ix * dx, y = -L + iy * dx, z = -L + iz * dx;
        for (int b = 0; b < nballs; b++) {
            double dxa = x - cx[b], dya = y - cy[b], dza = z - cz[b];
            double r = sqrt(dxa*dxa + dya*dya + dza*dza);
            double f = interp(&prof[b], r);
            if (f == 0.0) continue;
            double cd = cos(delta[b]), sd = sin(delta[b]);
            for (int a = 0; a < 3; a++) {
                col[0 + a][i]  += (float)(f * cd);              /* u_a   */
                col[12 + a][i] += (float)(f * sd);              /* v_a   */
                col[6 + a][i]  += (float)(-omega[b] * f * sd);  /* udot  */
                col[18 + a][i] += (float)( omega[b] * f * cd);  /* vdot  */
            }
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
    printf("gen_qball_multi: wrote %s (%d^3, f32, 24 columns, %d balls)\n",
           outpath, N, nballs);
    return 0;
}
