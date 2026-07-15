/* gen_mf_pair_boost.c — multi-fabric B1 seeds: separate heavy (C) and light (L) SFAs.
 *
 * Writes two 24-col complex seeds (no gauge; kernel Gauss-projects E):
 *   out_C.sfa — ball at (+D/2,0,0), omega_C, boost (vxC,vyC,vzC)
 *   out_L.sfa — ball at (-D/2,0,0), omega_L, boost (vxL,vyL,vzL)
 *
 * Usage (single profile — both fabrics):
 *   gen_mf_pair_boost N L D profile omegaC omegaL \
 *     vxC vyC vzC  vxL vyL vzL  out_C.sfa out_L.sfa
 *
 * Usage (dual profile — hierarchy / lighter L):
 *   gen_mf_pair_boost N L D profC omegaC profL omegaL \
 *     vxC vyC vzC  vxL vyL vzL  out_C.sfa out_L.sfa
 *
 * IMPORTANT (B1 charge assignment q_Q=+1, q_L=-1):
 *   Use the SAME sign of omega for C and L. Fabric charges already provide
 *   opposite EM (ρ_em = +ρ_C − ρ_L). Opposite omega *and* opposite q double-
 *   flips L → same-sign sources + opposite force charges → null relative
 *   force with COM co-drift (v75 G2 wrong_opp_omega).
 *
 * Typical equal-mass force:  omegaC=omegaL=+1.46, all v=0
 * Typical hierarchy force:   profC=f_w142, ωC=1.42; profL=f_w146, ωL=1.46
 * Typical head-on:           same-sign ω, vxC=-vr, vxL=+vr
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

static void write_one(int N, double L, const Profile *prof, double omega,
                      double cx, double cy, double cz,
                      double vx, double vy, double vz, const char *path) {
    double v2 = vx*vx+vy*vy+vz*vz;
    if (v2 >= 0.81) { fprintf(stderr, "FATAL: |v| too large\n"); exit(1); }
    double gamma = (v2 > 1e-30) ? 1.0/sqrt(1.0-v2) : 1.0;
    long N3 = (long)N*N*N;
    float *cols[NCOLS];
    for (int c=0;c<NCOLS;c++) { cols[c]=calloc(N3,sizeof(float)); if(!cols[c]) exit(1); }
    double dx = 2.0*L/(N-1);
    for (int i=0;i<N;i++) {
        double x = -L + i*dx - cx;
        for (int j=0;j<N;j++) {
            double y = -L + j*dx - cy;
            for (int k=0;k<N;k++) {
                double z = -L + k*dx - cz;
                long idx = (long)i*N*N + (long)j*N + k;
                /* boost along arbitrary v: longitudinal contract */
                double vdotx = vx*x + vy*y + vz*z;
                double r2_perp = x*x+y*y+z*z;
                double r_par2 = 0, r_perp2 = r2_perp;
                if (v2 > 1e-30) {
                    double vhat_x=vx/sqrt(v2), vhat_y=vy/sqrt(v2), vhat_z=vz/sqrt(v2);
                    double xpar = x*vhat_x+y*vhat_y+z*vhat_z;
                    r_par2 = xpar*xpar;
                    r_perp2 = r2_perp - r_par2;
                    if (r_perp2 < 0) r_perp2 = 0;
                } else {
                    r_par2 = 0; r_perp2 = r2_perp;
                }
                double rp = sqrt(gamma*gamma*r_par2 + r_perp2);
                double f = interp(prof, rp);
                double fp = interp_deriv(prof, rp);
                double phiB = -gamma * omega * vdotx; /* de Broglie */
                double cB = cos(phiB), sB = sin(phiB);
                double drdt = 0;
                if (rp > 1e-14 && v2 > 1e-30) {
                    /* d(r')/dt at t=0 ~ -gamma^2 (v·x)/r' for pure-x boost; general: */
                    drdt = -gamma*gamma * vdotx / rp;
                }
                double udot = fp*drdt*cB - f*gamma*omega*sB;
                double vdot = fp*drdt*sB + f*gamma*omega*cB;
                for (int a=0;a<3;a++) {
                    cols[a][idx] = (float)(f * cB);           /* phi */
                    cols[3+a][idx] = 0;                       /* theta */
                    cols[6+a][idx] = (float)udot;             /* phi_vel real */
                    cols[9+a][idx] = 0;                       /* theta_vel */
                    cols[12+a][idx] = (float)(f * sB);        /* phi_im */
                    cols[15+a][idx] = 0;
                    cols[18+a][idx] = (float)vdot;            /* phi_im_vel */
                    cols[21+a][idx] = 0;
                }
            }
        }
    }
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
    for (int c=0;c<NCOLS;c++) free(cols[c]);
    printf("Wrote %s (N=%d L=%.1f omega=%.4f center=(%.2f,%.2f,%.2f) v=(%.4f,%.4f,%.4f))\n",
           path, N, L, omega, cx, cy, cz, vx, vy, vz);
}

int main(int argc, char **argv) {
    /* argc 15: single profile; argc 16: dual profiles */
    if (argc != 15 && argc != 16) {
        fprintf(stderr,
            "Usage: %s N L D profile omegaC omegaL vxC vyC vzC vxL vyL vzL out_C.sfa out_L.sfa\n"
            "   or: %s N L D profC omegaC profL omegaL vxC vyC vzC vxL vyL vzL out_C.sfa out_L.sfa\n",
            argv[0], argv[0]);
        return 1;
    }
    int N = atoi(argv[1]);
    double L = atof(argv[2]), D = atof(argv[3]);
    Profile profC, profL;
    double oC, oL, vxC, vyC, vzC, vxL, vyL, vzL;
    const char *outC, *outL;
    if (argc == 15) {
        load_profile(argv[4], &profC);
        profL = profC; /* same arrays — free only once */
        oC = atof(argv[5]); oL = atof(argv[6]);
        vxC=atof(argv[7]); vyC=atof(argv[8]); vzC=atof(argv[9]);
        vxL=atof(argv[10]); vyL=atof(argv[11]); vzL=atof(argv[12]);
        outC = argv[13]; outL = argv[14];
        write_one(N, L, &profC, oC, +D/2, 0, 0, vxC, vyC, vzC, outC);
        write_one(N, L, &profL, oL, -D/2, 0, 0, vxL, vyL, vzL, outL);
        free(profC.r); free(profC.f);
    } else {
        load_profile(argv[4], &profC);
        oC = atof(argv[5]);
        load_profile(argv[6], &profL);
        oL = atof(argv[7]);
        vxC=atof(argv[8]); vyC=atof(argv[9]); vzC=atof(argv[10]);
        vxL=atof(argv[11]); vyL=atof(argv[12]); vzL=atof(argv[13]);
        outC = argv[14]; outL = argv[15];
        write_one(N, L, &profC, oC, +D/2, 0, 0, vxC, vyC, vzC, outC);
        write_one(N, L, &profL, oL, -D/2, 0, 0, vxL, vyL, vzL, outL);
        free(profC.r); free(profC.f);
        free(profL.r); free(profL.f);
    }
    return 0;
}
