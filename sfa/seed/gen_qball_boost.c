/* gen_qball_boost.c — stamp ONE boosted Q-ball into an SFA seed for the v66
 * complexified kernel (complex_phi=1, 12 fields + 12 velocities = 24 columns).
 *
 * Physics (v67/theta_dynamics/DEBROGLIE.md §1, theorem D1): at eta = 0 the
 * phi sector is exactly Lorentz invariant, so if Phi_a = f(r) e^{i omega t}
 * solves the EOM, the boost along x
 *
 *     Phi_a(t,x) = f(r'(t)) e^{i gamma omega (t - v x)},
 *     r'(t) = sqrt(gamma^2 (x - v t)^2 + y^2 + z^2),  gamma = 1/sqrt(1-v^2)
 *
 * is an EXACT moving solution (c = 1). This generator stamps the t = 0 slice
 * and its exact time derivative, for all a = 0,1,2 (symmetric ansatz):
 *
 *     r'      = sqrt(gamma^2 x^2 + y^2 + z^2)        (contracted envelope)
 *     phi_B   = -gamma omega v x                     (de Broglie phase tilt)
 *     u_a     =  f(r') cos(phi_B)
 *     v_a     =  f(r') sin(phi_B)
 *
 *     dr'/dt |_{t=0} = -gamma^2 v x / r'             (chain rule on r'(t))
 *     d(phase)/dt at fixed x = +gamma omega
 *     udot_a  = f'(r') (dr'/dt) cos(phi_B) - f(r') gamma omega sin(phi_B)
 *     vdot_a  = f'(r') (dr'/dt) sin(phi_B) + f(r') gamma omega cos(phi_B)
 *
 * f'(r') is the derivative of the piecewise-linear profile interpolant
 * (exact segment slope = finite difference on the profile grid).
 * Theta sector (and velocities) zero. Ball centered at the origin.
 *
 * Expected fingerprints (DEBROGLIE.md): phase stripes k_dB = gamma omega v,
 * fixed-voxel omega_core = gamma omega (blueshift), comoving beat omega/gamma,
 * E = gamma E0, p = gamma E0 v, translation at speed v.
 *
 * Output: ONE SFA frame, f32, 24 columns named/typed exactly as the kernel
 * registers them (scp_sim.c sfa_add_column block / v66 SPEC §7.1).
 *
 * Build:
 *   gcc -O3 -march=native -o bin/gen_qball_boost sfa/seed/gen_qball_boost.c -lzstd -lm
 *
 * Usage:
 *   gen_qball_boost N L profile omega vx out.sfa
 *   (boost along x only; |vx| < 0.9; profile = radial_qball output:
 *    '#' headers then "r f" rows)
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

/* locate segment lo with r[lo] <= r < r[lo+1] (caller guarantees in range) */
static size_t seg_lo(const Profile *p, double r) {
    size_t lo = 0, hi = p->n - 1;
    while (hi - lo > 1) {
        size_t mid = (lo + hi) / 2;
        if (p->r[mid] <= r) lo = mid; else hi = mid;
    }
    return lo;
}

/* piecewise-linear interpolation; flat below first r, zero above last */
static double interp(const Profile *p, double r) {
    if (r <= p->r[0]) return p->f[0];
    if (r >= p->r[p->n - 1]) return 0.0;
    size_t lo = seg_lo(p, r);
    double t = (r - p->r[lo]) / (p->r[lo + 1] - p->r[lo]);
    return p->f[lo] + t * (p->f[lo + 1] - p->f[lo]);
}

/* derivative of the interpolant: exact segment slope (finite difference on
 * the profile grid); 0 outside the tabulated range (flat/zero extensions) */
static double interp_deriv(const Profile *p, double r) {
    if (r <= p->r[0] || r >= p->r[p->n - 1]) return 0.0;
    size_t lo = seg_lo(p, r);
    return (p->f[lo + 1] - p->f[lo]) / (p->r[lo + 1] - p->r[lo]);
}

int main(int argc, char **argv) {
    if (argc != 7) {
        fprintf(stderr,
            "Usage: %s N L profile omega vx out.sfa\n"
            "  boost along x only, |vx| < 0.9; ball centered at origin;\n"
            "  profile = radial_qball output ('#' headers then \"r f\")\n",
            argv[0]);
        return 1;
    }
    int    N     = atoi(argv[1]);
    double L     = atof(argv[2]);
    const char *ppath = argv[3];
    double omega = atof(argv[4]);
    double vx    = atof(argv[5]);
    const char *outpath = argv[6];

    if (N < 2 || L <= 0.0) {
        fprintf(stderr, "FATAL: bad N=%d or L=%g\n", N, L);
        return 1;
    }
    if (fabs(vx) >= 0.9) {
        fprintf(stderr, "FATAL: |vx|=%g >= 0.9 (envelope under-resolved; "
                "use |vx| < 0.9)\n", fabs(vx));
        return 1;
    }

    double gamma = 1.0 / sqrt(1.0 - vx * vx);
    double gw    = gamma * omega;        /* lab (fixed-point) frequency  */
    double kdB   = gamma * omega * vx;   /* de Broglie tilt wavenumber   */

    printf("gen_qball_boost: N=%d L=%.4f -> %s\n", N, L, outpath);
    Profile prof;
    load_profile(ppath, &prof);
    printf("  omega=%.6f vx=%.6f gamma=%.6f\n", omega, vx, gamma);
    printf("  lab frequency gamma*omega=%.6f, comoving beat omega/gamma=%.6f\n",
           gw, omega / gamma);
    printf("  k_dB=gamma*omega*vx=%.6f (stripe wavelength %.3f)\n",
           kdB, (kdB != 0.0) ? 2.0 * M_PI / fabs(kdB) : INFINITY);

    long N3 = (long)N * N * N;
    double dx = 2.0 * L / (N - 1);   /* grid: x = -L + i*dx (scp_init.h) */
    if (kdB != 0.0 && 2.0 * M_PI / fabs(kdB) < 8.0 * dx)
        printf("  WARNING: de Broglie wavelength %.3f < 8 dx=%.3f "
               "(phase tilt under-resolved)\n", 2.0 * M_PI / fabs(kdB), 8.0 * dx);

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

    for (int i = 0; i < N; i++) {
        double x = -L + i * dx;
        double gx = gamma * x;
        for (int j = 0; j < N; j++) {
            double y = -L + j * dx;
            for (int k = 0; k < N; k++) {
                double z = -L + k * dx;
                long idx = (long)i * N * N + (long)j * N + k;

                double rp = sqrt(gx * gx + y * y + z * z); /* r' contracted */
                double f  = interp(&prof, rp);
                double fp = interp_deriv(&prof, rp);
                /* dr'/dt|_{t=0} = -gamma^2 vx x / r' (0 at the center where
                 * x = 0 forces the numerator to 0 first) */
                double drdt = (rp > 1e-12) ? (-gamma * gamma * vx * x / rp) : 0.0;

                double ph = -kdB * x;            /* phi_B = -gamma omega vx x */
                double cp = cos(ph), sp = sin(ph);

                double u  = f * cp;
                double v  = f * sp;
                double ud = fp * drdt * cp - f * gw * sp;
                double vd = fp * drdt * sp + f * gw * cp;

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
    char bN[32], bL[32], bo[32], bv[32], bg[32], bk[32];
    snprintf(bN, sizeof(bN), "%d", N);
    snprintf(bL, sizeof(bL), "%.6f", L);
    snprintf(bo, sizeof(bo), "%.6f", omega);
    snprintf(bv, sizeof(bv), "%.6f", vx);
    snprintf(bg, sizeof(bg), "%.6f", gamma);
    snprintf(bk, sizeof(bk), "%.6f", kdB);
    const char *keys[] = { "generator", "N", "L", "omega", "vx", "gamma", "k_dB" };
    const char *vals[] = { "gen_qball_boost", bN, bL, bo, bv, bg, bk };
    sfa_add_kvmd(sfa, 0, 0xFFFFFFFF, 0xFFFFFFFF, keys, vals, 7);

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

    printf("gen_qball_boost: wrote %s (%d^3, f32, 24 columns, 1 frame)\n",
           outpath, N);

    for (int c = 0; c < NCOLS; c++) free(cols[c]);
    free(prof.r);
    free(prof.f);
    return 0;
}
