/*  radial_eta_soliton.c — consistent (Phi, Theta) seed for the stationary
 *  eta-coupled Q-ball (the radiation-free geometry; FUTURE.md drift-floor work).
 *
 *  CONVENTION (fixed 2026-07-08, v72): the radial_qball profile f(r) is the
 *  PER-COMPONENT amplitude (gen_qball_boost convention): Phi_a = f(r) for each
 *  a, i.e. the component-space vector is sqrt(3) f(r) uhat, uhat=(1,1,1)/sqrt3,
 *  and s = prod|Phi_a|^2 = f^6. (This tool previously seeded Phi_a = f/sqrt3 —
 *  off-shell by 27x in s; all June-26 drift baselines carried that bug.) Curl:
 *      curl Phi = sqrt(3) f'(r) (rhat x uhat),
 *  which sources a TRANSVERSE torsion partner Theta = g(r) (rhat x uhat).
 *  Seeding Theta=0 leaves a ~1.4% drift over T=60 (the "eta-drain", actually an
 *  inconsistent-IC transient). The EXACT stationary g(r) solves the linear
 *  radial BVP (vector Laplacian of an l=1 transverse field gives the -2/r^2
 *  centrifugal term):
 *
 *      g'' + (2/r) g' - (2/r^2) g - (m_theta^2 - omega^2) g = -sqrt(3) eta f'
 *      g(0)=0  (g ~ r near origin),   g(r_max)=0  (Yukawa decay, kappa^2>0)
 *
 *  This tool reads a radial_qball profile (r f), solves g(r) by the Thomas
 *  algorithm, and writes a 24-column complex seed with Phi (symmetric ball) and
 *  the consistent Theta, both rotating at omega. Pass mode "zero" to seed
 *  Theta=0 (control) or "half" for g=f'/2, else the BVP solution.
 *
 *  NOTE (scope): the eta back-reaction on Phi (eta*curl Theta) has a component
 *  along rhat, so the fully self-consistent stationary solution is axisymmetric
 *  f(r,theta), not purely radial. This tool builds the EXACT transverse Theta for
 *  a given radial f and tests how much of the drain that removes; the residual
 *  measures the Phi back-reaction (the next, 2D, BVP step).
 *
 *  Build: gcc -O3 -fopenmp -o radial_eta_soliton radial_eta_soliton.c -lzstd -lm
 *  Usage: radial_eta_soliton N L profile.txt omega eta m_theta out.sfa [bvp|zero|half]
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define SFA_IMPLEMENTATION
#include "../format/sfa.h"
#define NCOLS 24

int main(int argc, char **argv) {
    if (argc < 8) { fprintf(stderr,"usage: %s N L profile omega eta m_theta out.sfa [bvp|zero|half]\n",argv[0]); return 1; }
    int N=atoi(argv[1]); double L=atof(argv[2]); const char *prof=argv[3];
    double omega=atof(argv[4]), eta=atof(argv[5]), mtheta=atof(argv[6]);
    const char *out=argv[7];
    const char *mode = (argc>8)? argv[8] : "bvp";

    /* read profile (r f), skipping '#' headers */
    int cap=8192, np=0; double *pr=malloc(cap*sizeof(double)), *pf=malloc(cap*sizeof(double));
    FILE *fp=fopen(prof,"r"); if(!fp){fprintf(stderr,"cannot open %s\n",prof);return 1;}
    char line[512];
    while(fgets(line,sizeof(line),fp)){ if(line[0]=='#'||line[0]=='\n')continue;
        double a,b; if(sscanf(line,"%lf %lf",&a,&b)==2){ if(np>=cap){cap*=2;pr=realloc(pr,cap*sizeof(double));pf=realloc(pf,cap*sizeof(double));} pr[np]=a;pf[np]=b;np++; } }
    fclose(fp);
    if(np<4){fprintf(stderr,"profile too short\n");return 1;}
    double dr=pr[1]-pr[0];

    /* f'(r) by central differences */
    double *df=malloc(np*sizeof(double));
    for(int i=1;i<np-1;i++) df[i]=(pf[i+1]-pf[i-1])/(2*dr);
    df[0]=(pf[1]-pf[0])/dr; df[np-1]=0;

    /* solve g(r): A_i g_{i-1} + B_i g_i + C_i g_{i+1} = D_i, g_0=g_{np-1}=0 (Thomas) */
    double *g=calloc(np,sizeof(double));
    double k2 = mtheta*mtheta - omega*omega;
    if (!strcmp(mode,"zero")) {
        /* g stays 0 */
    } else if (!strcmp(mode,"half")) {
        for(int i=0;i<np;i++) g[i]=0.5*sqrt(3.0)*df[i];   /* |(1/2)curl Phi| */
    } else { /* bvp */
        double *A=malloc(np*sizeof(double)),*B=malloc(np*sizeof(double)),*C=malloc(np*sizeof(double)),*D=malloc(np*sizeof(double));
        double *cp=malloc(np*sizeof(double)),*dp=malloc(np*sizeof(double));
        for(int i=1;i<np-1;i++){ double r=pr[i];
            A[i]=1.0/(dr*dr) - 1.0/(r*dr);
            B[i]=-2.0/(dr*dr) - 2.0/(r*r) - k2;
            C[i]=1.0/(dr*dr) + 1.0/(r*dr);
            D[i]=-sqrt(3.0)*eta*df[i];
        }
        /* Dirichlet g_0=0, g_{np-1}=0 : fold into ends */
        int lo=1, hi=np-2;
        cp[lo]=C[lo]/B[lo]; dp[lo]=D[lo]/B[lo];
        for(int i=lo+1;i<=hi;i++){ double m=B[i]-A[i]*cp[i-1]; cp[i]=C[i]/m; dp[i]=(D[i]-A[i]*dp[i-1])/m; }
        g[hi]=dp[hi];
        for(int i=hi-1;i>=lo;i--) g[i]=dp[i]-cp[i]*g[i+1];
        free(A);free(B);free(C);free(D);free(cp);free(dp);
    }
    double gmax=0; for(int i=0;i<np;i++) if(fabs(g[i])>gmax)gmax=fabs(g[i]);

    /* dump (r,f,g) table next to the seed */
    char tab[600]; snprintf(tab,sizeof(tab),"%s.fg.tsv",out);
    FILE *tf=fopen(tab,"w"); if(tf){ fprintf(tf,"r\tf\tdf\tg\n");
        for(int i=0;i<np;i++) fprintf(tf,"%.4f\t%.6e\t%.6e\t%.6e\n",pr[i],pf[i],df[i],g[i]); fclose(tf); }

    /* write seed: Phi_a=f(r) per component (real), Theta_a=g(r)(rhat x uhat)_a
     * (real), both rotating at omega */
    long N3=(long)N*N*N; double dx=2.0*L/(N-1), inv=1.0/sqrt(3.0);
    double ux=inv,uy=inv,uz=inv;
    float *cols[NCOLS]; for(int c=0;c<NCOLS;c++) cols[c]=calloc((size_t)N3,sizeof(float));
    double rmax=pr[np-1];
    #pragma omp parallel for
    for(long ii=0;ii<N3;ii++){
        int i=(int)(ii/((long)N*N)),j=(int)((ii/N)%N),k=(int)(ii%N);
        double x=-L+i*dx,y=-L+j*dx,z=-L+k*dx, r=sqrt(x*x+y*y+z*z);
        double f,gg;
        if(r>=rmax){ f=0; gg=0; }
        else { double t=r/dr; int idx=(int)t; if(idx>=np-1)idx=np-2; double frac=t-idx;
               f=pf[idx]+frac*(pf[idx+1]-pf[idx]); gg=g[idx]+frac*(g[idx+1]-g[idx]); }
        /* Phi_a = f per component (f is the per-component amplitude) */
        double pa[3]={f,f,f};
        /* rhat x u */
        double rx=(r>1e-9)?x/r:0, ry=(r>1e-9)?y/r:0, rz=(r>1e-9)?z/r:0;
        double cx=ry*uz-rz*uy, cy=rz*ux-rx*uz, cz=rx*uy-ry*ux;
        double ta[3]={gg*cx,gg*cy,gg*cz};
        for(int a=0;a<3;a++){
            cols[0 +a][ii]=(float)pa[a];          /* phi u */
            cols[18+a][ii]=(float)(omega*pa[a]);  /* phiim_v = vdot = omega*u */
            cols[3 +a][ii]=(float)ta[a];          /* theta re */
            cols[21+a][ii]=(float)(omega*ta[a]);  /* thetaim_v = omega*Theta_re */
        }
    }
    double dt=0.025*dx;
    SFA *sfa=sfa_create(out,N,N,N,L,L,L,dt);
    sfa_add_column(sfa,"phi_x",SFA_F32,SFA_POSITION,0);sfa_add_column(sfa,"phi_y",SFA_F32,SFA_POSITION,1);sfa_add_column(sfa,"phi_z",SFA_F32,SFA_POSITION,2);
    sfa_add_column(sfa,"theta_x",SFA_F32,SFA_ANGLE,0);sfa_add_column(sfa,"theta_y",SFA_F32,SFA_ANGLE,1);sfa_add_column(sfa,"theta_z",SFA_F32,SFA_ANGLE,2);
    sfa_add_column(sfa,"phi_vx",SFA_F32,SFA_VELOCITY,0);sfa_add_column(sfa,"phi_vy",SFA_F32,SFA_VELOCITY,1);sfa_add_column(sfa,"phi_vz",SFA_F32,SFA_VELOCITY,2);
    sfa_add_column(sfa,"theta_vx",SFA_F32,SFA_VELOCITY,3);sfa_add_column(sfa,"theta_vy",SFA_F32,SFA_VELOCITY,4);sfa_add_column(sfa,"theta_vz",SFA_F32,SFA_VELOCITY,5);
    sfa_add_column(sfa,"phiim_x",SFA_F32,SFA_POSITION,3);sfa_add_column(sfa,"phiim_y",SFA_F32,SFA_POSITION,4);sfa_add_column(sfa,"phiim_z",SFA_F32,SFA_POSITION,5);
    sfa_add_column(sfa,"thetaim_x",SFA_F32,SFA_ANGLE,3);sfa_add_column(sfa,"thetaim_y",SFA_F32,SFA_ANGLE,4);sfa_add_column(sfa,"thetaim_z",SFA_F32,SFA_ANGLE,5);
    sfa_add_column(sfa,"phiim_vx",SFA_F32,SFA_VELOCITY,6);sfa_add_column(sfa,"phiim_vy",SFA_F32,SFA_VELOCITY,7);sfa_add_column(sfa,"phiim_vz",SFA_F32,SFA_VELOCITY,8);
    sfa_add_column(sfa,"thetaim_vx",SFA_F32,SFA_VELOCITY,9);sfa_add_column(sfa,"thetaim_vy",SFA_F32,SFA_VELOCITY,10);sfa_add_column(sfa,"thetaim_vz",SFA_F32,SFA_VELOCITY,11);
    sfa_finalize_header(sfa);
    void *fc[NCOLS]; for(int c=0;c<NCOLS;c++) fc[c]=cols[c];
    sfa_write_frame(sfa,0.0,fc); sfa_close(sfa);
    printf("radial_eta_soliton: %s mode=%s  omega=%.3f eta=%.3f m_theta=%.3f kappa^2=%.4f g_max=%.4e\n",
           out, mode, omega, eta, mtheta, k2, gmax);
    for(int c=0;c<NCOLS;c++) free(cols[c]);
    free(pr);free(pf);free(df);free(g);
    return 0;
}
