/* gen_mf_shell_orbit.c — B4-style MF seeds: heavy C at origin + L shell with
 * optional tangential boosts for orbit kinematics.
 *
 * Writes two 24-col complex SFAs (no gauge):
 *   out_C.sfa — n_C balls (default 1 at origin) with boosts
 *   out_L.sfa — n_L balls on a sphere/shell with vt × ê_φ (tangential)
 *
 * Usage:
 *   gen_mf_shell_orbit N L R profC omegaC profL omegaL n_L vt out_C.sfa out_L.sfa
 *
 *   n_L ∈ {1,2,4,6}; shell geometry:
 *     1: (R,0,0)
 *     2: (±R,0,0)
 *     4: tetra on sphere radius R
 *     6: axial octahedron at ±R on each axis
 *   vt = tangential speed (c=1). vt=0 → rest shell (B4 rest).
 *   C is always a single rest ball at origin with profC/omegaC.
 *
 * B1 charge: same-sign omega C and L; q_Q=+1, q_L=-1 in cfg.
 *
 * Build:
 *   gcc -O3 -o bin/gen_mf_shell_orbit sfa/seed/gen_mf_shell_orbit.c -lzstd -lm
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define SFA_IMPLEMENTATION
#include "../format/sfa.h"
#define NCOLS 24
#define MAXB 16

typedef struct { double *r, *f; size_t n; } Profile;

static void load_profile(const char *path, Profile *p) {
    FILE *fp = fopen(path, "r");
    if (!fp) { fprintf(stderr, "FATAL: open %s\n", path); exit(1); }
    size_t cap = 1024; p->n = 0;
    p->r = malloc(cap * sizeof(double)); p->f = malloc(cap * sizeof(double));
    char line[512];
    while (fgets(line, sizeof(line), fp)) {
        const char *s = line;
        while (*s==' '||*s=='\t') s++;
        if (*s=='#'||*s=='\n'||!*s||*s=='\r') continue;
        double rv, fv;
        if (sscanf(s, "%lf %lf", &rv, &fv) != 2) continue;
        if (p->n == cap) { cap *= 2; p->r = realloc(p->r, cap*sizeof(double)); p->f = realloc(p->f, cap*sizeof(double)); }
        p->r[p->n] = rv; p->f[p->n] = fv; p->n++;
    }
    fclose(fp);
    if (p->n < 2) { fprintf(stderr, "FATAL: short profile\n"); exit(1); }
}
static double interp(const Profile *p, double r) {
    if (r <= p->r[0]) return p->f[0];
    if (r >= p->r[p->n-1]) return 0.0;
    size_t lo=0, hi=p->n-1;
    while (hi-lo>1) { size_t mid=(lo+hi)/2; if (p->r[mid]<=r) lo=mid; else hi=mid; }
    double t = (r - p->r[lo]) / (p->r[hi] - p->r[lo]);
    return p->f[lo] + t * (p->f[hi] - p->f[lo]);
}
static double interp_deriv(const Profile *p, double r) {
    if (r <= p->r[0] || r >= p->r[p->n-1]) return 0.0;
    size_t lo=0, hi=p->n-1;
    while (hi-lo>1) { size_t mid=(lo+hi)/2; if (p->r[mid]<=r) lo=mid; else hi=mid; }
    return (p->f[hi] - p->f[lo]) / (p->r[hi] - p->r[lo]);
}

/* Superpose one boosted ball (+=) into cols */
static void add_ball(float **cols, int N, double L, const Profile *prof,
                     double omega, double cx, double cy, double cz,
                     double vx, double vy, double vz) {
    double v2 = vx*vx+vy*vy+vz*vz;
    if (v2 >= 0.81) { fprintf(stderr, "FATAL: |v| too large\n"); exit(1); }
    double gamma = (v2 > 1e-30) ? 1.0/sqrt(1.0-v2) : 1.0;
    double dx = 2.0*L/(N-1);
    for (int i=0;i<N;i++) {
        double x = -L + i*dx - cx;
        for (int j=0;j<N;j++) {
            double y = -L + j*dx - cy;
            for (int k=0;k<N;k++) {
                double z = -L + k*dx - cz;
                long idx = (long)i*N*N + (long)j*N + k;
                double vdotx = vx*x + vy*y + vz*z;
                double r2 = x*x+y*y+z*z;
                double r_par2 = 0, r_perp2 = r2;
                if (v2 > 1e-30) {
                    double invv = 1.0/sqrt(v2);
                    double xpar = (x*vx+y*vy+z*vz)*invv;
                    r_par2 = xpar*xpar;
                    r_perp2 = r2 - r_par2;
                    if (r_perp2 < 0) r_perp2 = 0;
                }
                double rp = sqrt(gamma*gamma*r_par2 + r_perp2);
                double f = interp(prof, rp);
                if (f == 0.0 && rp > 0.1) continue;
                double fp = interp_deriv(prof, rp);
                double phiB = -gamma * omega * vdotx;
                double cB = cos(phiB), sB = sin(phiB);
                double drdt = 0;
                if (rp > 1e-14 && v2 > 1e-30)
                    drdt = -gamma*gamma * vdotx / rp;
                double udot = fp*drdt*cB - f*gamma*omega*sB;
                double vdot = fp*drdt*sB + f*gamma*omega*cB;
                for (int a=0;a<3;a++) {
                    cols[a][idx]      += (float)(f * cB);
                    cols[6+a][idx]    += (float)udot;
                    cols[12+a][idx]   += (float)(f * sB);
                    cols[18+a][idx]   += (float)vdot;
                }
            }
        }
    }
}

static void write_sfa(const char *path, int N, double L, float **cols) {
    SFA *sfa = sfa_create(path, N, N, N, L, L, L, 0.01);
    if (!sfa) { fprintf(stderr, "FATAL: create %s\n", path); exit(1); }
    sfa->flags = SFA_CODEC_COLZSTD | SFA_FLAG_STREAMING;
    const char *names[NCOLS] = {
        "phi_x","phi_y","phi_z","theta_x","theta_y","theta_z",
        "phi_vx","phi_vy","phi_vz","theta_vx","theta_vy","theta_vz",
        "phiim_x","phiim_y","phiim_z","thetaim_x","thetaim_y","thetaim_z",
        "phiim_vx","phiim_vy","phiim_vz","thetaim_vx","thetaim_vy","thetaim_vz"
    };
    int sem[NCOLS] = {
        SFA_POSITION,SFA_POSITION,SFA_POSITION,SFA_ANGLE,SFA_ANGLE,SFA_ANGLE,
        SFA_VELOCITY,SFA_VELOCITY,SFA_VELOCITY,SFA_VELOCITY,SFA_VELOCITY,SFA_VELOCITY,
        SFA_POSITION,SFA_POSITION,SFA_POSITION,SFA_ANGLE,SFA_ANGLE,SFA_ANGLE,
        SFA_VELOCITY,SFA_VELOCITY,SFA_VELOCITY,SFA_VELOCITY,SFA_VELOCITY,SFA_VELOCITY
    };
    int comp[NCOLS] = {0,1,2,0,1,2, 0,1,2,3,4,5, 3,4,5,3,4,5, 6,7,8,9,10,11};
    for (int c=0;c<NCOLS;c++) sfa_add_column(sfa, names[c], SFA_F32, sem[c], comp[c]);
    sfa_finalize_header(sfa);
    sfa_write_frame(sfa, 0.0, (void**)cols);
    sfa_close(sfa);
    printf("Wrote %s\n", path);
}

static int shell_centers(int n, double R, double c[][3], double vtang[][3], double vt) {
    /* positions + tangential velocity ê_φ = ê_z × r̂  (or ê_x × r̂ if polar) */
    if (n == 1) {
        c[0][0]=R; c[0][1]=0; c[0][2]=0;
        vtang[0][0]=0; vtang[0][1]=vt; vtang[0][2]=0;
        return 1;
    }
    if (n == 2) {
        c[0][0]=R; c[0][1]=0; c[0][2]=0; vtang[0][0]=0; vtang[0][1]=vt; vtang[0][2]=0;
        c[1][0]=-R; c[1][1]=0; c[1][2]=0; vtang[1][0]=0; vtang[1][1]=-vt; vtang[1][2]=0;
        return 2;
    }
    if (n == 4) {
        double raw[4][3] = {{1,1,1},{1,-1,-1},{-1,1,-1},{-1,-1,1}};
        for (int i=0;i<4;i++) {
            double nr = sqrt(raw[i][0]*raw[i][0]+raw[i][1]*raw[i][1]+raw[i][2]*raw[i][2]);
            c[i][0]=R*raw[i][0]/nr; c[i][1]=R*raw[i][1]/nr; c[i][2]=R*raw[i][2]/nr;
            /* ê_φ ≈ (-y, x, 0) normalized */
            double tx=-c[i][1], ty=c[i][0], tz=0;
            double tn=sqrt(tx*tx+ty*ty+tz*tz);
            if (tn < 1e-12) { tx=0; ty=-c[i][2]; tz=c[i][1]; tn=sqrt(tx*tx+ty*ty+tz*tz); }
            vtang[i][0]=vt*tx/tn; vtang[i][1]=vt*ty/tn; vtang[i][2]=vt*tz/tn;
        }
        return 4;
    }
    if (n == 6) {
        /* Axial octahedron; co-rotate about z: v = vt ê_z × r̂ (poles rest). */
        double ax[6][3] = {{R,0,0},{-R,0,0},{0,R,0},{0,-R,0},{0,0,R},{0,0,-R}};
        for (int i=0;i<6;i++) {
            c[i][0]=ax[i][0]; c[i][1]=ax[i][1]; c[i][2]=ax[i][2];
            /* ê_z × r = (-y, x, 0); |ê_z × r| = R_xy */
            double rxy = sqrt(c[i][0]*c[i][0] + c[i][1]*c[i][1]);
            if (rxy < 1e-12) {
                vtang[i][0]=0; vtang[i][1]=0; vtang[i][2]=0; /* poles */
            } else {
                vtang[i][0] = vt * (-c[i][1] / rxy);
                vtang[i][1] = vt * ( c[i][0] / rxy);
                vtang[i][2] = 0;
            }
        }
        return 6;
    }
    fprintf(stderr, "FATAL: n_L=%d not in {1,2,4,6}\n", n);
    exit(1);
}

int main(int argc, char **argv) {
    if (argc != 12) {
        fprintf(stderr,
            "Usage: %s N L R profC omegaC profL omegaL n_L vt out_C.sfa out_L.sfa\n",
            argv[0]);
        return 1;
    }
    int N = atoi(argv[1]);
    double Lbox = atof(argv[2]), R = atof(argv[3]);
    Profile profC, profL;
    load_profile(argv[4], &profC);
    double oC = atof(argv[5]);
    load_profile(argv[6], &profL);
    double oL = atof(argv[7]);
    int nL = atoi(argv[8]);
    double vt = atof(argv[9]);
    const char *outC = argv[10], *outL = argv[11];

    long N3 = (long)N*N*N;
    float *colsC[NCOLS], *colsL[NCOLS];
    for (int c=0;c<NCOLS;c++) {
        colsC[c]=calloc(N3,sizeof(float));
        colsL[c]=calloc(N3,sizeof(float));
        if (!colsC[c]||!colsL[c]) exit(1);
    }

    /* C: single rest ball at origin */
    add_ball(colsC, N, Lbox, &profC, oC, 0, 0, 0, 0, 0, 0);
    write_sfa(outC, N, Lbox, colsC);

    double cen[MAXB][3], vel[MAXB][3];
    int nb = shell_centers(nL, R, cen, vel, vt);
    printf("L shell: n=%d R=%.2f vt=%.4f\n", nb, R, vt);
    for (int b=0;b<nb;b++) {
        printf("  L[%d] center=(%.2f,%.2f,%.2f) v=(%.4f,%.4f,%.4f)\n",
               b, cen[b][0], cen[b][1], cen[b][2], vel[b][0], vel[b][1], vel[b][2]);
        add_ball(colsL, N, Lbox, &profL, oL,
                 cen[b][0], cen[b][1], cen[b][2],
                 vel[b][0], vel[b][1], vel[b][2]);
    }
    write_sfa(outL, N, Lbox, colsL);

    for (int c=0;c<NCOLS;c++) { free(colsC[c]); free(colsL[c]); }
    free(profC.r); free(profC.f); free(profL.r); free(profL.f);
    return 0;
}
