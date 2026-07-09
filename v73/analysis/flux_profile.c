/*  flux_profile.c — the fabric LEDGER: local energy/charge flux from SFA frames.
 *
 *  Process-form stability (v73/PROCESS.md): a particle persists iff the fabric
 *  it takes up equals the fabric it lays down — the time-averaged flux through
 *  every closed surface vanishes. This tool measures that ledger.
 *
 *  Energy flux (derived from the continuity of the coupled system; the eta
 *  coupling L_int = eta[tu.(curl u) + tv.(curl v)] contributes a cross term):
 *      S_i = - sum_f  fdot d_i f   +  eta [ (udot x tu)_i + (vdot x tv)_i ]
 *  (sum over all 12 real fields; the cross product treats the component index
 *  as the spatial index, matching the kernel's curl convention.)
 *  U(1) charge current (continuity-consistent sign, J = v rho for a boost):
 *      J_i = sum_a [ v_a d_i u_a - u_a d_i v_a ] + (theta terms)
 *
 *  Mode 1 (default): shell-averaged radial profiles, time-averaged over frames
 *  with t in [t0, t1]:  r, <e>, <S_r>free, <S_r>eta, <S_r>tot, <J_r>,
 *  P(r) = 4 pi r^2 <S_r>tot (outward power through the shell).
 *  Mode 2 (--probe x): per-frame time series at grid point (x, 0, 0):
 *  t, u0, v0, |Phi|^2, e, rhoQ  (fabric-element history for uptake/layment).
 *
 *  Build: gcc -O3 -fopenmp -o flux_profile flux_profile.c -lzstd -lm
 *  Usage: flux_profile file.sfa eta [--t0 A --t1 B] [--out f.tsv]
 *                      [--probe x]... (up to 8 probes)
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define SFA_IMPLEMENTATION
#include "../../sfa/format/sfa.h"

static int N; static long N3, NN;
static inline long IX(int i,int j,int k){ return ((long)i*N+j)*N+k; }

int main(int argc, char **argv) {
    if (argc < 3) { fprintf(stderr, "usage: %s file.sfa eta [--t0 A --t1 B] [--out f.tsv] [--probe x]...\n", argv[0]); return 1; }
    const char *path = argv[1];
    double eta = atof(argv[2]);
    double t0 = 0, t1 = 1e30; const char *outp = NULL;
    double probes[8]; int nprobe = 0;
    for (int a = 3; a < argc; a++) {
        if (!strcmp(argv[a], "--t0") && a+1 < argc) t0 = atof(argv[++a]);
        else if (!strcmp(argv[a], "--t1") && a+1 < argc) t1 = atof(argv[++a]);
        else if (!strcmp(argv[a], "--out") && a+1 < argc) outp = argv[++a];
        else if (!strcmp(argv[a], "--probe") && a+1 < argc && nprobe < 8) probes[nprobe++] = atof(argv[++a]);
    }

    SFA *s = sfa_open(path);
    if (!s) { fprintf(stderr, "cannot open %s\n", path); return 1; }
    N = (int)s->Nx; N3 = (long)N*N*N; NN = (long)N*N;
    double L = s->Lx, dx = 2.0*L/(N-1), inv2dx = 1.0/(2.0*dx), dV = dx*dx*dx;

    /* column map: 12 fields + 12 velocities by semantic/component */
    int cu[3],cv[3],ctu[3],ctv[3], cud[3],cvd[3],ctud[3],ctvd[3];
    for (int a=0;a<3;a++) cu[a]=cv[a]=ctu[a]=ctv[a]=cud[a]=cvd[a]=ctud[a]=ctvd[a]=-1;
    for (uint32_t c = 0; c < s->n_columns; c++) {
        int sem = s->columns[c].semantic, cmp = s->columns[c].component;
        if (sem == SFA_POSITION && cmp < 3) cu[cmp] = c;
        else if (sem == SFA_POSITION && cmp < 6) cv[cmp-3] = c;
        else if (sem == SFA_ANGLE && cmp < 3) ctu[cmp] = c;
        else if (sem == SFA_ANGLE && cmp < 6) ctv[cmp-3] = c;
        else if (sem == SFA_VELOCITY && cmp < 3) cud[cmp] = c;
        else if (sem == SFA_VELOCITY && cmp < 6) ctud[cmp-3] = c;
        else if (sem == SFA_VELOCITY && cmp < 9) cvd[cmp-6] = c;
        else if (sem == SFA_VELOCITY && cmp < 12) ctvd[cmp-9] = c;
    }
    for (int a=0;a<3;a++)
        if (cu[a]<0||cv[a]<0||ctu[a]<0||ctv[a]<0||cud[a]<0||cvd[a]<0||ctud[a]<0||ctvd[a]<0) {
            fprintf(stderr, "missing columns (need 24-col complex file)\n"); return 1; }

    void *buf = malloc(s->frame_bytes);
    /* per-column base offsets into buf (all dtypes must be F32 here) */
    uint64_t *coff = malloc(s->n_columns * sizeof(uint64_t));
    { uint64_t off = 0;
      for (uint32_t c = 0; c < s->n_columns; c++) {
          if (s->columns[c].dtype != SFA_F32) { fprintf(stderr, "non-f32 column\n"); return 1; }
          coff[c] = off; off += (uint64_t)N3 * 4; } }
    #define F(c, ii) ((float*)((uint8_t*)buf + coff[(c)]))[(ii)]

    int nbin = N/2;
    double *se = calloc(nbin,sizeof(double)), *sfree = calloc(nbin,sizeof(double));
    double *seta = calloc(nbin,sizeof(double)), *sj = calloc(nbin,sizeof(double));
    long *cnt = calloc(nbin,sizeof(long));
    int navg = 0;

    FILE *pf[8]; char pnm[8][256];
    for (int p = 0; p < nprobe; p++) {
        snprintf(pnm[p],256,"%s.probe%.1f.tsv", outp?outp:"flux", probes[p]);
        pf[p] = fopen(pnm[p],"w");
        fprintf(pf[p],"t\tu0\tv0\tphisq\te\trhoQ\n");
    }

    for (uint32_t fr = 0; fr < s->total_frames; fr++) {
        double t = sfa_frame_time(s, fr);
        if (t < t0 - 1e-9 || t > t1 + 1e-9) { if (!nprobe) continue; }
        if (sfa_read_frame(s, fr, buf) != 0) break;

        /* probes: fabric-element history at (x,0,0) */
        for (int p = 0; p < nprobe; p++) {
            int i = (int)lround((probes[p] + L)/dx), j = N/2, k = N/2;
            if (i < 0 || i >= N) continue;
            long ii = IX(i,j,k);
            double u0=F(cu[0],ii), v0=F(cv[0],ii);
            double phisq=0, e=0, rq=0;
            for (int a=0;a<3;a++){
                double u=F(cu[a],ii), v=F(cv[a],ii), ud=F(cud[a],ii), vd=F(cvd[a],ii);
                double tu=F(ctu[a],ii), tv=F(ctv[a],ii), tud=F(ctud[a],ii), tvd=F(ctvd[a],ii);
                phisq += u*u+v*v;
                e += 0.5*(ud*ud+vd*vd+tud*tud+tvd*tvd);   /* kinetic only at probe (cheap) */
                rq += u*vd - v*ud;
            }
            fprintf(pf[p],"%.3f\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\n",t,u0,v0,phisq,e,rq);
        }
        if (t < t0 - 1e-9 || t > t1 + 1e-9) continue;

        /* shell-averaged radial fluxes */
        #pragma omp parallel
        {
            double *le = calloc(nbin,sizeof(double)), *lf = calloc(nbin,sizeof(double));
            double *lt = calloc(nbin,sizeof(double)), *lj = calloc(nbin,sizeof(double));
            long *lc = calloc(nbin,sizeof(long));
            #pragma omp for schedule(static)
            for (long ii = 0; ii < N3; ii++) {
                int i=(int)(ii/NN), j=(int)((ii/N)%N), k=(int)(ii%N);
                if (i<1||i>=N-1||j<1||j>=N-1||k<1||k>=N-1) continue;
                double x=-L+i*dx, y=-L+j*dx, z=-L+k*dx, r=sqrt(x*x+y*y+z*z);
                int b = (int)(r/dx); if (b >= nbin || r < 1e-9) continue;
                double rx=x/r, ry=y/r, rz=z/r;
                long Xp=IX(i+1,j,k),Xm=IX(i-1,j,k),Yp=IX(i,j+1,k),Ym=IX(i,j-1,k),Zp=IX(i,j,k+1),Zm=IX(i,j,k-1);
                double Sfx=0,Sfy=0,Sfz=0, Jx=0,Jy=0,Jz=0, e=0;
                double udv[3],vdv[3],tuv[3],tvv[3];   /* velocities as component-vectors */
                double tuvec[3],tvvec[3];
                for (int a=0;a<3;a++){
                    double u=F(cu[a],ii), v=F(cv[a],ii), tu=F(ctu[a],ii), tv=F(ctv[a],ii);
                    double ud=F(cud[a],ii), vd=F(cvd[a],ii), tud=F(ctud[a],ii), tvd=F(ctvd[a],ii);
                    double gux=(F(cu[a],Xp)-F(cu[a],Xm))*inv2dx, guy=(F(cu[a],Yp)-F(cu[a],Ym))*inv2dx, guz=(F(cu[a],Zp)-F(cu[a],Zm))*inv2dx;
                    double gvx=(F(cv[a],Xp)-F(cv[a],Xm))*inv2dx, gvy=(F(cv[a],Yp)-F(cv[a],Ym))*inv2dx, gvz=(F(cv[a],Zp)-F(cv[a],Zm))*inv2dx;
                    double gtux=(F(ctu[a],Xp)-F(ctu[a],Xm))*inv2dx, gtuy=(F(ctu[a],Yp)-F(ctu[a],Ym))*inv2dx, gtuz=(F(ctu[a],Zp)-F(ctu[a],Zm))*inv2dx;
                    double gtvx=(F(ctv[a],Xp)-F(ctv[a],Xm))*inv2dx, gtvy=(F(ctv[a],Yp)-F(ctv[a],Ym))*inv2dx, gtvz=(F(ctv[a],Zp)-F(ctv[a],Zm))*inv2dx;
                    Sfx -= ud*gux + vd*gvx + tud*gtux + tvd*gtvx;
                    Sfy -= ud*guy + vd*gvy + tud*gtuy + tvd*gtvy;
                    Sfz -= ud*guz + vd*gvz + tud*gtuz + tvd*gtvz;
                    Jx += v*gux - u*gvx + tv*gtux - tu*gtvx;
                    Jy += v*guy - u*gvy + tv*gtuy - tu*gtvy;
                    Jz += v*guz - u*gvz + tv*gtuz - tu*gtvz;
                    e  += 0.5*(ud*ud+vd*vd+tud*tud+tvd*tvd)
                        + 0.5*(gux*gux+guy*guy+guz*guz + gvx*gvx+gvy*gvy+gvz*gvz
                              +gtux*gtux+gtuy*gtuy+gtuz*gtuz + gtvx*gtvx+gtvy*gtvy+gtvz*gtvz);
                    udv[a]=ud; vdv[a]=vd; tuvec[a]=tu; tvvec[a]=tv;
                    (void)tuv;(void)tvv;
                }
                /* eta cross-flux: + eta (udot x tu + vdot x tv), component index = spatial */
                double Sex = eta*(udv[1]*tuvec[2]-udv[2]*tuvec[1] + vdv[1]*tvvec[2]-vdv[2]*tvvec[1]);
                double Sey = eta*(udv[2]*tuvec[0]-udv[0]*tuvec[2] + vdv[2]*tvvec[0]-vdv[0]*tvvec[2]);
                double Sez = eta*(udv[0]*tuvec[1]-udv[1]*tuvec[0] + vdv[0]*tvvec[1]-vdv[1]*tvvec[0]);
                double srf = Sfx*rx+Sfy*ry+Sfz*rz, sre = Sex*rx+Sey*ry+Sez*rz;
                double jr = Jx*rx+Jy*ry+Jz*rz;
                le[b]+=e; lf[b]+=srf; lt[b]+=sre; lj[b]+=jr; lc[b]++;
            }
            #pragma omp critical
            { for (int b=0;b<nbin;b++){ se[b]+=le[b]; sfree[b]+=lf[b]; seta[b]+=lt[b]; sj[b]+=lj[b]; cnt[b]+=lc[b]; } }
            free(le);free(lf);free(lt);free(lj);free(lc);
        }
        navg++;
    }

    FILE *of = outp ? fopen(outp,"w") : stdout;
    fprintf(of,"# %s eta=%.3f frames_averaged=%d t=[%g,%g]\n", path, eta, navg, t0, t1);
    fprintf(of,"r\te\tSr_free\tSr_eta\tSr_tot\tJr\tP_out\n");
    for (int b = 1; b < nbin; b++) {
        if (!cnt[b]) continue;
        double r = (b+0.5)*dx, n = (double)cnt[b];
        double sr = (sfree[b]+seta[b])/n;
        fprintf(of,"%.3f\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\n",
                r, se[b]/n, sfree[b]/n, seta[b]/n, sr, sj[b]/n, 4.0*M_PI*r*r*sr);
    }
    if (outp) fclose(of);
    for (int p = 0; p < nprobe; p++) { fclose(pf[p]); fprintf(stderr,"probe -> %s\n",pnm[p]); }
    (void)dV;
    sfa_close(s);
    return 0;
}
