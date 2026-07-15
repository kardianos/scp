/* gen_qball_pair_boost.c — two boosted Q-balls (orbit / collision seeds)
 *
 * Each ball is an exact eta=0 Lorentz boost of the radial rotating ansatz,
 * centered at (cx,cy,cz) with velocity (vx,vy,vz). Superposition is linear.
 * Negative omega ⇒ opposite charge (anti-ball).
 *
 * Usage:
 *   gen_qball_pair_boost N L out.sfa \
 *     profile1 omega1 delta1 x1 y1 z1 vx1 vy1 vz1 \
 *     profile2 omega2 delta2 x2 y2 z2 vx2 vy2 vz2
 *
 * Build:
 *   gcc -O3 -march=native -fopenmp -o bin/gen_qball_pair_boost \
 *       sfa/seed/gen_qball_pair_boost.c -lzstd -lm
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
    if (!fp) { fprintf(stderr, "FATAL: cannot open '%s'\n", path); exit(1); }
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
}

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

static double interp_deriv(const Profile *p, double r) {
    if (r <= p->r[0] || r >= p->r[p->n - 1]) return 0.0;
    size_t lo = 0, hi = p->n - 1;
    while (hi - lo > 1) {
        size_t mid = (lo + hi) / 2;
        if (p->r[mid] <= r) lo = mid; else hi = mid;
    }
    return (p->f[hi] - p->f[lo]) / (p->r[hi] - p->r[lo]);
}

/* Stamp one boosted ball (additive) into 24 columns. */
static void stamp_boosted(float **cols, int N, double L,
                          const Profile *prof, double omega, double delta,
                          double cx, double cy, double cz,
                          double vx, double vy, double vz) {
    double speed = sqrt(vx * vx + vy * vy + vz * vz);
    double dxg = 2.0 * L / (N - 1);
    double cdel = cos(delta), sdel = sin(delta);
    long N3 = (long)N * N * N;

    if (speed < 1e-12) {
        #pragma omp parallel for
        for (long idx = 0; idx < N3; idx++) {
            int i = (int)(idx / ((long)N * N));
            int j = (int)((idx / N) % N);
            int k = (int)(idx % N);
            double x = -L + i * dxg - cx;
            double y = -L + j * dxg - cy;
            double z = -L + k * dxg - cz;
            double r = sqrt(x * x + y * y + z * z);
            double f = interp(prof, r);
            if (f == 0.0) continue;
            double u = f * cdel, v = f * sdel;
            double ud = -omega * f * sdel, vd = omega * f * cdel;
            for (int a = 0; a < 3; a++) {
                cols[0 + a][idx]  += (float)u;
                cols[12 + a][idx] += (float)v;
                cols[6 + a][idx]  += (float)ud;
                cols[18 + a][idx] += (float)vd;
            }
        }
        return;
    }
    if (speed >= 0.9) {
        fprintf(stderr, "FATAL: |v|=%g >= 0.9\n", speed);
        exit(1);
    }
    double hx = vx / speed, hy = vy / speed, hz = vz / speed;
    double gamma = 1.0 / sqrt(1.0 - speed * speed);
    double gw = gamma * omega;
    double kdB = gamma * omega * speed;

    #pragma omp parallel for
    for (long idx = 0; idx < N3; idx++) {
        int i = (int)(idx / ((long)N * N));
        int j = (int)((idx / N) % N);
        int k = (int)(idx % N);
        double x = -L + i * dxg - cx;
        double y = -L + j * dxg - cy;
        double z = -L + k * dxg - cz;

        double xpar = x * hx + y * hy + z * hz;
        double px = x - xpar * hx;
        double py = y - xpar * hy;
        double pz = z - xpar * hz;
        double rp = sqrt(gamma * gamma * xpar * xpar + px * px + py * py + pz * pz);
        double f = interp(prof, rp);
        if (f == 0.0) continue;
        double fp = interp_deriv(prof, rp);
        double drdt = (rp > 1e-12) ? (-gamma * gamma * speed * xpar / rp) : 0.0;
        double ph = -kdB * xpar + delta;
        double cp = cos(ph), sp = sin(ph);
        double u = f * cp;
        double v = f * sp;
        double ud = fp * drdt * cp - f * gw * sp;
        double vd = fp * drdt * sp + f * gw * cp;
        for (int a = 0; a < 3; a++) {
            cols[0 + a][idx]  += (float)u;
            cols[12 + a][idx] += (float)v;
            cols[6 + a][idx]  += (float)ud;
            cols[18 + a][idx] += (float)vd;
        }
    }
}

int main(int argc, char **argv) {
    /* argc: prog + N L out + 9*2 ball args = 22 */
    if (argc != 22) {
        fprintf(stderr,
            "Usage: %s N L out.sfa \\\n"
            "  prof1 omega1 delta1 x1 y1 z1 vx1 vy1 vz1 \\\n"
            "  prof2 omega2 delta2 x2 y2 z2 vx2 vy2 vz2\n"
            "  (delta rad; neg omega = opposite charge; |v|<0.9)\n"
            "  (got argc=%d, want 22)\n",
            argv[0], argc);
        return 1;
    }
    int N = atoi(argv[1]);
    double L = atof(argv[2]);
    const char *outpath = argv[3];

    Profile p1, p2;
    load_profile(argv[4], &p1);
    double o1 = atof(argv[5]), d1 = atof(argv[6]);
    double c1x = atof(argv[7]), c1y = atof(argv[8]), c1z = atof(argv[9]);
    double v1x = atof(argv[10]), v1y = atof(argv[11]), v1z = atof(argv[12]);

    load_profile(argv[13], &p2);
    double o2 = atof(argv[14]), d2 = atof(argv[15]);
    double c2x = atof(argv[16]), c2y = atof(argv[17]), c2z = atof(argv[18]);
    double v2x = atof(argv[19]), v2y = atof(argv[20]), v2z = atof(argv[21]);

    printf("gen_qball_pair_boost: N=%d L=%g -> %s\n", N, L, outpath);
    printf("  ball1: ω=%+.4f δ=%.3f c=(%.2f,%.2f,%.2f) v=(%.4f,%.4f,%.4f)\n",
           o1, d1, c1x, c1y, c1z, v1x, v1y, v1z);
    printf("  ball2: ω=%+.4f δ=%.3f c=(%.2f,%.2f,%.2f) v=(%.4f,%.4f,%.4f)\n",
           o2, d2, c2x, c2y, c2z, v2x, v2y, v2z);

    long N3 = (long)N * N * N;
    double dx = 2.0 * L / (N - 1);
    float *cols[NCOLS];
    for (int c = 0; c < NCOLS; c++) {
        cols[c] = (float *)calloc((size_t)N3, sizeof(float));
        if (!cols[c]) { fprintf(stderr, "FATAL: OOM\n"); return 1; }
    }

    stamp_boosted(cols, N, L, &p1, o1, d1, c1x, c1y, c1z, v1x, v1y, v1z);
    stamp_boosted(cols, N, L, &p2, o2, d2, c2x, c2y, c2z, v2x, v2y, v2z);

    double dt = 0.025 * dx;
    SFA *sfa = sfa_create(outpath, N, N, N, L, L, L, dt);
    if (!sfa) { fprintf(stderr, "FATAL: sfa_create\n"); return 1; }

    char bN[32], bL[32];
    snprintf(bN, sizeof(bN), "%d", N);
    snprintf(bL, sizeof(bL), "%.6f", L);
    const char *keys[] = { "generator", "N", "L" };
    const char *vals[] = { "gen_qball_pair_boost", bN, bL };
    sfa_add_kvmd(sfa, 0, 0xFFFFFFFF, 0xFFFFFFFF, keys, vals, 3);

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
    printf("gen_qball_pair_boost: wrote %s\n", outpath);

    for (int c = 0; c < NCOLS; c++) free(cols[c]);
    free(p1.r); free(p1.f); free(p2.r); free(p2.f);
    return 0;
}
