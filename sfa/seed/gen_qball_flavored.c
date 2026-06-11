/*  gen_qball_flavored.c — seed FLAVORED (asymmetric) 3-component Q-balls
 *
 *  Stamps one or more balls, each with per-component radial profiles and
 *  per-component frequencies (the in-model flavor structure):
 *
 *      Phi_a(x, t=0) += f_a(|x-c|) ,  vdot_a += w_a f_a
 *
 *  Profile files: 4 columns "r f0 f1 f2" (flavored, from
 *  v71/analysis/flavored_qball.py) or 2 columns "r f" (symmetric — replicated
 *  to all three components). '#' comments ignored.
 *
 *  Output: 24-column complex SFA (matter only).
 *
 *  Build: gcc -O3 -fopenmp -o gen_qball_flavored gen_qball_flavored.c -lzstd -lm
 *  Usage:
 *    gen_qball_flavored N L out.sfa  profile w0 w1 w2 x y z  [.. more balls ..]
 *    (one 7-arg group per ball, max 8)
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define SFA_IMPLEMENTATION
#include "../format/sfa.h"

#define NCOLS    24
#define MAXBALLS 8

typedef struct {
    size_t n;
    double *r, *f[3];
} Prof;

static void load_prof(const char *path, Prof *p) {
    FILE *fp = fopen(path, "r");
    if (!fp) { fprintf(stderr, "FATAL: cannot open '%s'\n", path); exit(1); }
    size_t cap = 1024;
    p->r = malloc(cap * sizeof(double));
    for (int a = 0; a < 3; a++) p->f[a] = malloc(cap * sizeof(double));
    p->n = 0;
    int ncol = 0;
    char line[1024];
    while (fgets(line, sizeof(line), fp)) {
        const char *s = line;
        while (*s == ' ' || *s == '\t') s++;
        if (*s == '#' || *s == '\n' || *s == '\0' || *s == '\r') continue;
        double rv, a0, a1, a2;
        int got = sscanf(s, "%lf %lf %lf %lf", &rv, &a0, &a1, &a2);
        if (got < 2) continue;
        if (ncol == 0) ncol = got;
        if (p->n == cap) {
            cap *= 2;
            p->r = realloc(p->r, cap * sizeof(double));
            for (int a = 0; a < 3; a++) p->f[a] = realloc(p->f[a], cap * sizeof(double));
        }
        p->r[p->n] = rv;
        if (got >= 4) {
            p->f[0][p->n] = a0; p->f[1][p->n] = a1; p->f[2][p->n] = a2;
        } else {
            p->f[0][p->n] = p->f[1][p->n] = p->f[2][p->n] = a0;  /* symmetric */
        }
        p->n++;
    }
    fclose(fp);
    if (p->n < 2) { fprintf(stderr, "FATAL: profile '%s' too short\n", path); exit(1); }
    printf("  profile '%s': %zu rows (%s), f(0) = %.5f %.5f %.5f\n",
           path, p->n, ncol >= 4 ? "flavored 4-col" : "symmetric 2-col",
           p->f[0][0], p->f[1][0], p->f[2][0]);
}

static double interp(const Prof *p, int a, double r) {
    if (r <= p->r[0]) return p->f[a][0];
    if (r >= p->r[p->n - 1]) return 0.0;
    size_t lo = 0, hi = p->n - 1;
    while (hi - lo > 1) { size_t mid = (lo + hi) / 2; if (p->r[mid] <= r) lo = mid; else hi = mid; }
    double t = (r - p->r[lo]) / (p->r[hi] - p->r[lo]);
    return p->f[a][lo] + t * (p->f[a][hi] - p->f[a][lo]);
}

int main(int argc, char **argv) {
    if (argc < 11 || (argc - 4) % 7 != 0) {
        fprintf(stderr,
            "Usage: %s N L out.sfa  profile w0 w1 w2 x y z  [.. more balls ..]\n"
            "  one 7-arg group per ball (max %d); profile = 2-col (symmetric) or 4-col (flavored)\n",
            argv[0], MAXBALLS);
        return 1;
    }
    int N = atoi(argv[1]);
    double L = atof(argv[2]);
    const char *outpath = argv[3];
    int nballs = (argc - 4) / 7;
    if (nballs > MAXBALLS) { fprintf(stderr, "FATAL: too many balls\n"); return 1; }

    Prof prof[MAXBALLS];
    double w[MAXBALLS][3], cx[MAXBALLS], cy[MAXBALLS], cz[MAXBALLS];
    printf("gen_qball_flavored: N=%d L=%g, %d ball(s) -> %s\n", N, L, nballs, outpath);
    for (int b = 0; b < nballs; b++) {
        char **g = &argv[4 + 7 * b];
        load_prof(g[0], &prof[b]);
        for (int a = 0; a < 3; a++) w[b][a] = atof(g[1 + a]);
        cx[b] = atof(g[4]); cy[b] = atof(g[5]); cz[b] = atof(g[6]);
        printf("  ball %d: w=(%.4f,%.4f,%.4f) center=(%.2f,%.2f,%.2f)\n",
               b, w[b][0], w[b][1], w[b][2], cx[b], cy[b], cz[b]);
    }

    long N3 = (long)N * N * N, NN = (long)N * N;
    double dx = 2.0 * L / (N - 1);

    float **col = malloc(NCOLS * sizeof(float *));
    for (int c = 0; c < NCOLS; c++) col[c] = calloc(N3, sizeof(float));
    #pragma omp parallel for
    for (long i = 0; i < N3; i++) {
        int ix = (int)(i / NN), iy = (int)((i / N) % N), iz = (int)(i % N);
        double x = -L + ix * dx, y = -L + iy * dx, z = -L + iz * dx;
        for (int b = 0; b < nballs; b++) {
            double dxa = x - cx[b], dya = y - cy[b], dza = z - cz[b];
            double r = sqrt(dxa*dxa + dya*dya + dza*dza);
            for (int a = 0; a < 3; a++) {
                double f = interp(&prof[b], a, r);
                if (f == 0.0) continue;
                col[0 + a][i]  += (float)f;               /* u_a   */
                col[18 + a][i] += (float)(w[b][a] * f);   /* vdot_a */
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
    void **arrays = malloc(NCOLS * sizeof(void *));
    for (int c = 0; c < NCOLS; c++) arrays[c] = col[c];
    if (sfa_write_frame(sfa, 0.0, arrays) != 0) { fprintf(stderr, "FATAL: write failed\n"); return 1; }
    sfa_close(sfa);
    printf("gen_qball_flavored: wrote %s (%d balls)\n", outpath, nballs);
    return 0;
}
