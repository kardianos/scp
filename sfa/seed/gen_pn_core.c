/* gen_pn_core.c — proton/neutron multi-fabric nuclear seeds (P/N, Z/N).
 *
 * Writes:
 *   out_C.sfa  — bag fabric: ALL nuclear balls (Z protons + N neutrons)
 *   out_Q.sfa  — charge fabric: ONLY Z proton balls (neutrons absent → Q_em∝Z)
 *   out_L.sfa  — optional light fabric (zero if nL=0)
 *
 * Physics contract (B2 unlock, q_C=0, q_Q=+1, q_L=−1):
 *   ρ_em = q_Q ρ_Q + q_L ρ_L
 *   Proton analog  = co-located C+Q lump (bag + charge)
 *   Neutron analog = C-only lump (bag mass, no EM source)
 *
 * Usage:
 *   gen_pn_core N L profile omega \
 *     out_C.sfa out_Q.sfa out_L.sfa \
 *     nZ  x0 y0 z0  [x1 y1 z1 ...] \
 *     nN  x0 y0 z0  [x1 y1 z1 ...] \
 *     nL  x0 y0 z0  [..]    # nL may be 0 (writes empty L)
 *
 * Build:
 *   gcc -O3 -fopenmp -o bin/gen_pn_core sfa/seed/gen_pn_core.c -lzstd -lm
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define SFA_IMPLEMENTATION
#include "../format/sfa.h"

#define NCOLS 24
#define MAXB  32

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
        if (p->n == cap) {
            cap *= 2;
            p->r = realloc(p->r, cap*sizeof(double));
            p->f = realloc(p->f, cap*sizeof(double));
        }
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

static void stamp_balls(int N, double L, const Profile *prof, double omega,
                        int nb, const double *cx, const double *cy, const double *cz,
                        float *cols[NCOLS]) {
    long N3 = (long)N*N*N, NN = (long)N*N;
    double dx = 2.0*L/(N-1);
    #pragma omp parallel for
    for (long i = 0; i < N3; i++) {
        int ix = (int)(i / NN), iy = (int)((i / N) % N), iz = (int)(i % N);
        double x = -L + ix*dx, y = -L + iy*dx, z = -L + iz*dx;
        for (int b = 0; b < nb; b++) {
            double dxa = x-cx[b], dya = y-cy[b], dza = z-cz[b];
            double r = sqrt(dxa*dxa + dya*dya + dza*dza);
            double f = interp(prof, r);
            if (f == 0.0) continue;
            for (int a = 0; a < 3; a++) {
                cols[0+a][i]  += (float)f;                 /* u_a */
                cols[12+a][i] += 0.0f;                     /* v_a */
                cols[6+a][i]  += 0.0f;                     /* udot */
                cols[18+a][i] += (float)(omega * f);       /* vdot: Phi ~ e^{iωt} */
            }
        }
    }
}

static void write_sfa(const char *path, int N, double L, float *cols[NCOLS]) {
    SFA *sfa = sfa_create(path, N, N, N, L, L, L, 0.0);
    if (!sfa) { fprintf(stderr, "FATAL: create %s\n", path); exit(1); }
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
    static const int comp[NCOLS] = {
        0,1,2, 0,1,2, 0,1,2, 3,4,5,
        3,4,5, 3,4,5, 6,7,8, 9,10,11};
    for (int c = 0; c < NCOLS; c++)
        sfa_add_column(sfa, names[c], SFA_F32, sem[c], comp[c]);
    sfa_finalize_header(sfa);
    void *arrays[NCOLS];
    for (int c = 0; c < NCOLS; c++) arrays[c] = cols[c];
    if (sfa_write_frame(sfa, 0.0, arrays) != 0) {
        fprintf(stderr, "FATAL: write %s\n", path); exit(1);
    }
    sfa_close(sfa);
    printf("  wrote %s\n", path);
}

static int read_centers(int argc, char **argv, int *argi, int n,
                        double *cx, double *cy, double *cz, const char *label) {
    if (n < 0 || n > MAXB) {
        fprintf(stderr, "FATAL: %s count %d out of range\n", label, n);
        return 1;
    }
    for (int b = 0; b < n; b++) {
        if (*argi + 2 >= argc) {
            fprintf(stderr, "FATAL: need 3 coords for %s ball %d\n", label, b);
            return 1;
        }
        cx[b] = atof(argv[(*argi)++]);
        cy[b] = atof(argv[(*argi)++]);
        cz[b] = atof(argv[(*argi)++]);
    }
    return 0;
}

int main(int argc, char **argv) {
    if (argc < 12) {
        fprintf(stderr,
            "Usage: %s N L profile omega out_C out_Q out_L \\\n"
            "         nZ xz0 yz0 zz0 [...]  nN xn0 yn0 zn0 [...]  nL xl0 yl0 zl0 [...]\n"
            "  nL may be 0 (empty light fabric).\n", argv[0]);
        return 1;
    }
    int N = atoi(argv[1]);
    double L = atof(argv[2]);
    const char *prof_path = argv[3];
    double omega = atof(argv[4]);
    const char *outC = argv[5], *outQ = argv[6], *outL = argv[7];

    Profile prof;
    load_profile(prof_path, &prof);
    printf("gen_pn_core: N=%d L=%g omega=%.4f profile=%s (f0=%.4f)\n",
           N, L, omega, prof_path, prof.f[0]);

    int argi = 8;
    int nZ = atoi(argv[argi++]);
    double Zx[MAXB], Zy[MAXB], Zz[MAXB];
    if (read_centers(argc, argv, &argi, nZ, Zx, Zy, Zz, "Z")) return 1;

    if (argi >= argc) { fprintf(stderr, "FATAL: missing nN\n"); return 1; }
    int nN = atoi(argv[argi++]);
    double Nx_[MAXB], Ny_[MAXB], Nz_[MAXB];
    if (read_centers(argc, argv, &argi, nN, Nx_, Ny_, Nz_, "N")) return 1;

    if (argi >= argc) { fprintf(stderr, "FATAL: missing nL\n"); return 1; }
    int nL = atoi(argv[argi++]);
    double Lx[MAXB], Ly[MAXB], Lz[MAXB];
    if (read_centers(argc, argv, &argi, nL, Lx, Ly, Lz, "L")) return 1;

    printf("  Z=%d (protons on C+Q)  N=%d (neutrons on C only)  L=%d\n", nZ, nN, nL);
    for (int b = 0; b < nZ; b++)
        printf("    p%d center=(%.2f,%.2f,%.2f)\n", b, Zx[b], Zy[b], Zz[b]);
    for (int b = 0; b < nN; b++)
        printf("    n%d center=(%.2f,%.2f,%.2f)\n", b, Nx_[b], Ny_[b], Nz_[b]);

    long N3 = (long)N*N*N;
    float *cC[NCOLS], *cQ[NCOLS], *cL[NCOLS];
    for (int c = 0; c < NCOLS; c++) {
        cC[c] = calloc(N3, sizeof(float));
        cQ[c] = calloc(N3, sizeof(float));
        cL[c] = calloc(N3, sizeof(float));
        if (!cC[c] || !cQ[c] || !cL[c]) { fprintf(stderr, "FATAL: OOM\n"); return 1; }
    }

    /* C bag: all nuclear (Z + N) */
    double Cx[MAXB*2], Cy[MAXB*2], Cz[MAXB*2];
    int nC = 0;
    for (int b = 0; b < nZ; b++) {
        Cx[nC]=Zx[b]; Cy[nC]=Zy[b]; Cz[nC]=Zz[b]; nC++;
    }
    for (int b = 0; b < nN; b++) {
        Cx[nC]=Nx_[b]; Cy[nC]=Ny_[b]; Cz[nC]=Nz_[b]; nC++;
    }
    stamp_balls(N, L, &prof, omega, nC, Cx, Cy, Cz, cC);

    /* Q charge: protons only */
    stamp_balls(N, L, &prof, omega, nZ, Zx, Zy, Zz, cQ);

    /* L light: optional */
    if (nL > 0)
        stamp_balls(N, L, &prof, omega, nL, Lx, Ly, Lz, cL);

    write_sfa(outC, N, L, cC);
    write_sfa(outQ, N, L, cQ);
    write_sfa(outL, N, L, cL);

    printf("gen_pn_core: done  A=%d Z=%d N=%d  (EM charge sources = Z balls on Q)\n",
           nZ + nN, nZ, nN);

    for (int c = 0; c < NCOLS; c++) {
        free(cC[c]); free(cQ[c]); free(cL[c]);
    }
    free(prof.r); free(prof.f);
    return 0;
}
