/*  eta_relaxer.c — 3D gradient-flow relaxer for the stationary eta-coupled
 *  Q-ball (the radiation-free geometry; FUTURE.md / stationary eta-soliton).
 *
 *  Minimizes the EFFECTIVE static energy of the rotating ansatz
 *  Phi_a = phi_a(x) e^{i w t}, Theta_a = theta_a(x) e^{i w t}:
 *
 *    E = sum_a [ 1/2|grad phi_a|^2 + 1/2(m^2-w^2)|phi_a|^2
 *              + 1/2|grad theta_a|^2 + 1/2(mth^2-w^2)|theta_a|^2 ]
 *        + Vt(s)  - eta * Re[ theta^bar . (curl phi) ]   + V_ext
 *    s = prod_a |phi_a|^2,   Vt(s) = (mu/2) s/(1+kappa s)
 *
 *  Forces (F = -dE/dfield), gradient flow field += dtau*F:
 *    F_phi   = lap phi   - (m^2 -w^2) phi   - Vt'(s) ds/dphi   + eta*curl(theta) - Vext'
 *    F_theta = lap theta - (mth^2-w^2) theta                   + eta*curl(phi)   - Vext'
 *  (phi,theta each feel eta x curl of the OTHER's same re/im part; the relaxer
 *  finds the self-consistent Theta AND the 2D deformation of phi automatically.)
 *
 *  EXTERNAL PRESSURE (key to convergence): V_ext = P * max(0, r-R)^2 confines the
 *  field to r<~R so a definite compact configuration forms instead of dispersing
 *  or collapsing during the flow. P is ANNEALED from P0 down to 0; once the caged
 *  configuration has converged, removing the pressure lets the true free soliton
 *  emerge. (force coeff = -2 P max(0,r-R) along the field.)
 *
 *  Ungauged (g=0): the eta-drain lives in the matter/Theta sector; the gauge adds
 *  Coulomb without a drain, and the kernel rebuilds E by Gauss projection when the
 *  relaxed seed is run. phi initialized from a radial_qball profile (symmetric
 *  ball, u=(1,1,1)/sqrt3), theta from 0.
 *
 *  Build: gcc -O3 -fopenmp -o eta_relaxer eta_relaxer.c -lzstd -lm
 *  Usage: eta_relaxer N L profile omega eta m mth mu kappa out.sfa
 *                     [P0 R nlevel iters dtau]
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define SFA_IMPLEMENTATION
#include "../format/sfa.h"
#define NCOLS 24

static int N; static long N3, NN;
static inline long IX(int i,int j,int k){ return ((long)i*N+j)*N+k; }

/* periodic-safe neighbor (fields ~0 at boundary anyway) */
static inline int wrp(int i){ return (i+N)%N; }

int main(int argc,char**argv){
    if(argc<11){fprintf(stderr,"usage: %s N L profile omega eta m mth mu kappa out.sfa [P0 R nlevel iters dtau]\n",argv[0]);return 1;}
    N=atoi(argv[1]); double L=atof(argv[2]); const char*prof=argv[3];
    double omega=atof(argv[4]), eta=atof(argv[5]), m=atof(argv[6]), mth=atof(argv[7]);
    double mu=atof(argv[8]), kappa=atof(argv[9]); const char*out=argv[10];
    double P0   = (argc>11)?atof(argv[11]):2.0;
    double Rpr  = (argc>12)?atof(argv[12]):6.0;
    int nlevel  = (argc>13)?atoi(argv[13]):6;
    int iters   = (argc>14)?atoi(argv[14]):3000;
    double dtau = (argc>15)?atof(argv[15]):0.0;
    /* Higgs "fabric-pull" condensate (self-compression bag). vH=0 -> off. */
    double vH   = (argc>16)?atof(argv[16]):0.0;   /* Higgs VEV (fabric vacuum amplitude) */
    double lamH = (argc>17)?atof(argv[17]):1.0;   /* quartic: V_H=(lamH/4)(H^2-vH^2)^2  (gap m_H=sqrt(2 lamH) vH) */
    double kapH = (argc>18)?atof(argv[18]):2.0;   /* matter-Higgs coupling (bag): +kapH/2 * s * H^2 */
    double Rvoid= (argc>19)?atof(argv[19]):0.0;   /* pre-dig the condensate cavity of this radius (0=flat) */
    /* 2nd gauge: relative (color) U(1) Coulomb potential B, charges q_a.
     * Component a's rotation phase is gauged: omega_a(x) = omega + g2 q_a B(x), so
     * the effective mass becomes m^2 - omega_a^2 (color-dependent binding). B is
     * sourced by the relative charge rho_B = sum_a q_a omega_a |phi_a|^2 (Gauss). */
    double g2  = (argc>20)?atof(argv[20]):0.0;    /* relative-gauge coupling (0 = off) */
    double q0  = (argc>21)?atof(argv[21]):1.0;
    double q1  = (argc>22)?atof(argv[22]):1.0;
    double q2  = (argc>23)?atof(argv[23]):-2.0;   /* default T8 = (1,1,-2), traceless */
    double qch[3]={q0,q1,q2};
    double m2w=m*m-omega*omega, mt2w=mth*mth-omega*omega;

    N3=(long)N*N*N; NN=(long)N*N; double dx=2.0*L/(N-1), dx2=dx*dx;
    if(dtau<=0) dtau=0.04*dx2;        /* diffusion-stable default (reaction-safe) */

    /* fields: u[3],v[3],tu[3],tv[3] + new buffers (Jacobi double-buffering) */
    double *u[3],*v[3],*tu[3],*tv[3], *un[3],*vn[3],*tun[3],*tvn[3];
    for(int a=0;a<3;a++){ u[a]=calloc(N3,sizeof(double));v[a]=calloc(N3,sizeof(double));
                          tu[a]=calloc(N3,sizeof(double));tv[a]=calloc(N3,sizeof(double));
                          un[a]=calloc(N3,sizeof(double));vn[a]=calloc(N3,sizeof(double));
                          tun[a]=calloc(N3,sizeof(double));tvn[a]=calloc(N3,sizeof(double)); }
    double inv=1.0/sqrt(3.0);
    /* init Phi: a ".gen" file is a blob genome (e.g. the CMA 4-node winner);
     * otherwise a radial profile (r f) -> symmetric ball Phi_a=f(r)/sqrt3. */
    int slen=(int)strlen(prof);
    int is_genome = (slen>4 && !strcmp(prof+slen-4,".gen"));
    if (is_genome) {
        enum{MAXB=64}; double bc[MAXB][3],bs[MAXB],ba[MAXB][6]; int nb=0; double gom=omega;
        FILE*fp=fopen(prof,"r"); if(!fp){fprintf(stderr,"cannot open %s\n",prof);return 1;}
        char ln[512];
        while(fgets(ln,sizeof(ln),fp)){
            if(!strncmp(ln,"omega",5)) sscanf(ln+5,"%lf",&gom);
            else if(!strncmp(ln,"blob",4)&&nb<MAXB){ double*cc=bc[nb],*A=ba[nb];
                if(sscanf(ln+4,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
                   &cc[0],&cc[1],&cc[2],&bs[nb],&A[0],&A[1],&A[2],&A[3],&A[4],&A[5])==10) nb++; } }
        fclose(fp);
        printf("# genome %s: %d blobs (genome omega=%.3f, relaxing at omega=%.3f)\n",prof,nb,gom,omega);
        #pragma omp parallel for
        for(long ii=0;ii<N3;ii++){ int i=(int)(ii/NN),j=(int)((ii/N)%N),k=(int)(ii%N);
            double x=-L+i*dx,y=-L+j*dx,z=-L+k*dx;
            for(int b=0;b<nb;b++){ double dxx=x-bc[b][0],dyy=y-bc[b][1],dzz=z-bc[b][2];
                double gg=exp(-(dxx*dxx+dyy*dyy+dzz*dzz)/(2.0*bs[b]*bs[b]));
                for(int a=0;a<3;a++){ u[a][ii]+=ba[b][2*a]*gg; v[a][ii]+=ba[b][2*a+1]*gg; } } }
    } else {
        int cap=8192,np=0; double *pr=malloc(cap*sizeof(double)),*pf=malloc(cap*sizeof(double));
        FILE*fp=fopen(prof,"r"); if(!fp){fprintf(stderr,"cannot open %s\n",prof);return 1;}
        char ln[512]; while(fgets(ln,sizeof(ln),fp)){ if(ln[0]=='#'||ln[0]=='\n')continue; double a,b;
            if(sscanf(ln,"%lf %lf",&a,&b)==2){ if(np>=cap){cap*=2;pr=realloc(pr,cap*sizeof(double));pf=realloc(pf,cap*sizeof(double));} pr[np]=a;pf[np]=b;np++; } }
        fclose(fp); double prdr=pr[1]-pr[0], rmax=pr[np-1];
        #pragma omp parallel for
        for(long ii=0;ii<N3;ii++){ int i=(int)(ii/NN),j=(int)((ii/N)%N),k=(int)(ii%N);
            double x=-L+i*dx,y=-L+j*dx,z=-L+k*dx,r=sqrt(x*x+y*y+z*z); double f=0;
            if(r<rmax){ double t=r/prdr; int id=(int)t; if(id>=np-1)id=np-2; double fr=t-id; f=pf[id]+fr*(pf[id+1]-pf[id]); }
            for(int a=0;a<3;a++) u[a][ii]=f*inv;
        }
        free(pr); free(pf);
    }
    /* 2nd gauge relative-Coulomb potential B (real, init 0) */
    double *B=calloc(N3,sizeof(double)),*Bn=calloc(N3,sizeof(double));
    /* Higgs condensate: VEV outside, optionally a pre-dug cavity (radius Rvoid)
     * so the bag holds the matter from t=0 (breaks the dig-vs-hold bootstrap). */
    double *H=calloc(N3,sizeof(double)),*Hn=calloc(N3,sizeof(double));
    #pragma omp parallel for
    for(long ii=0;ii<N3;ii++){ int i=(int)(ii/NN),j=(int)((ii/N)%N),k=(int)(ii%N);
        double x=-L+i*dx,y=-L+j*dx,z=-L+k*dx,r=sqrt(x*x+y*y+z*z);
        H[ii]= (Rvoid>0)? vH*0.5*(1.0+tanh((r-Rvoid)/0.6)) : vH; }

    /* precompute pressure coefficient pc(r)=2*max(0,r-R) (times P each level) */
    float *pc=malloc(N3*sizeof(float));
    #pragma omp parallel for
    for(long ii=0;ii<N3;ii++){ int i=(int)(ii/NN),j=(int)((ii/N)%N),k=(int)(ii%N);
        double x=-L+i*dx,y=-L+j*dx,z=-L+k*dx,r=sqrt(x*x+y*y+z*z); double e=r-Rpr; pc[ii]=(float)(e>0?2.0*e:0.0); }

    double invdx2=1.0/dx2, inv2dx=1.0/(2.0*dx);
    /* target matter norm M0 = int sum_a|phi_a|^2 : held fixed each step so the
     * flow converges to the soliton (constrained min) instead of the vacuum. */
    double M0=0;
    #pragma omp parallel for reduction(+:M0)
    for(long ii=0;ii<N3;ii++) for(int a=0;a<3;a++) M0+=u[a][ii]*u[a][ii]+v[a][ii]*v[a][ii];
    printf("# eta_relaxer N=%d L=%g w=%.3f eta=%.3f m2w=%.4f mt2w=%.4f mu=%g kappa=%g\n",N,L,omega,eta,m2w,mt2w,mu,kappa);
    printf("# anneal P0=%.3g R=%.2f levels=%d iters=%d dtau=%.4g\n",P0,Rpr,nlevel,iters,dtau);
    printf("# %-6s %-10s %-12s %-12s %-12s\n","level","P","E_eff","maxF","s_max");

    for(int lev=0; lev<nlevel; lev++){
        double P = (nlevel>1)? P0*(1.0 - (double)lev/(nlevel-1)) : P0;  /* P0 -> 0 linearly */
        for(int it=0; it<iters; it++){
            double maxF=0;
            /* one gradient-flow sweep (Jacobi: read old, write new via temp not needed if we accept SOR-like; use in-place is fine for relaxation) */
            #pragma omp parallel for reduction(max:maxF)
            for(long ii=0;ii<N3;ii++){
                int i=(int)(ii/NN),j=(int)((ii/N)%N),k=(int)(ii%N);
                int ip=wrp(i+1),im=wrp(i-1),jp=wrp(j+1),jm=wrp(j-1),kp=wrp(k+1),km=wrp(k-1);
                long Xp=IX(ip,j,k),Xm=IX(im,j,k),Yp=IX(i,jp,k),Ym=IX(i,jm,k),Zp=IX(i,j,kp),Zm=IX(i,j,km);
                double Pcc=P*pc[ii];
                /* s and P_a for potential */
                double r0=u[0][ii]*u[0][ii]+v[0][ii]*v[0][ii];
                double r1=u[1][ii]*u[1][ii]+v[1][ii]*v[1][ii];
                double r2=u[2][ii]*u[2][ii]+v[2][ii]*v[2][ii];
                double s=r0*r1*r2; double Vtp=0.5*mu/((1+kappa*s)*(1+kappa*s));
                double Pa[3]={r1*r2,r0*r2,r0*r1};
                double hh=(vH>0)? H[ii]*H[ii] : 0.0;   /* Higgs density at this cell */
                for(int a=0;a<3;a++){
                    /* laplacians */
                    double lu =(u[a][Xp]+u[a][Xm]+u[a][Yp]+u[a][Ym]+u[a][Zp]+u[a][Zm]-6.0*u[a][ii])*invdx2;
                    double lv =(v[a][Xp]+v[a][Xm]+v[a][Yp]+v[a][Ym]+v[a][Zp]+v[a][Zm]-6.0*v[a][ii])*invdx2;
                    double ltu=(tu[a][Xp]+tu[a][Xm]+tu[a][Yp]+tu[a][Ym]+tu[a][Zp]+tu[a][Zm]-6.0*tu[a][ii])*invdx2;
                    double ltv=(tv[a][Xp]+tv[a][Xm]+tv[a][Yp]+tv[a][Ym]+tv[a][Zp]+tv[a][Zm]-6.0*tv[a][ii])*invdx2;
                    /* curl(theta)_a and curl(phi)_a : (curl W)_a = eps_abc d_b W_c */
                    int b=(a+1)%3, c=(a+2)%3;
                    long bp,bm,cp,cm;
                    /* d_b W_c - d_c W_b along axes b,c */
                    #define DERIV(arr,ax,comp) ( (ax==0)?(arr[comp][Xp]-arr[comp][Xm]):(ax==1)?(arr[comp][Yp]-arr[comp][Ym]):(arr[comp][Zp]-arr[comp][Zm]) )
                    (void)bp;(void)bm;(void)cp;(void)cm;
                    double curl_tu = (DERIV(tu,b,c) - DERIV(tu,c,b))*inv2dx;
                    double curl_tv = (DERIV(tv,b,c) - DERIV(tv,c,b))*inv2dx;
                    double curl_u  = (DERIV(u, b,c) - DERIV(u, c,b))*inv2dx;
                    double curl_v  = (DERIV(v, b,c) - DERIV(v, c,b))*inv2dx;
                    #undef DERIV
                    /* Higgs bag back-reaction: matter is expelled from the condensate
                     * (force -kapH H^2 u_a P_a), so it digs a cavity the condensate squeezes. */
                    double hbag = kapH*hh*Pa[a];
                    /* 2nd gauge: color-dependent effective mass m^2 - omega_a^2 */
                    double oma = omega + g2*qch[a]*B[ii];
                    double m2wa = (g2!=0.0) ? (m*m - oma*oma) : m2w;
                    double Fu = lu  - m2wa*u[a][ii]  - Vtp*2.0*u[a][ii]*Pa[a] + eta*curl_tu - Pcc*u[a][ii] - hbag*u[a][ii];
                    double Fv = lv  - m2wa*v[a][ii]  - Vtp*2.0*v[a][ii]*Pa[a] + eta*curl_tv - Pcc*v[a][ii] - hbag*v[a][ii];
                    double Ftu= ltu - mt2w*tu[a][ii] + eta*curl_u  - Pcc*tu[a][ii];
                    double Ftv= ltv - mt2w*tv[a][ii] + eta*curl_v  - Pcc*tv[a][ii];
                    un[a][ii]=u[a][ii]+dtau*Fu; vn[a][ii]=v[a][ii]+dtau*Fv;
                    tun[a][ii]=tu[a][ii]+dtau*Ftu; tvn[a][ii]=tv[a][ii]+dtau*Ftv;
                    double af=fabs(Fu); if(af>maxF)maxF=af; af=fabs(Ftu); if(af>maxF)maxF=af;
                }
                /* Higgs condensate force: F_H = lap H - lamH H(H^2-vH^2) - kapH s H */
                if(vH>0){
                    double lH=(H[Xp]+H[Xm]+H[Yp]+H[Ym]+H[Zp]+H[Zm]-6.0*H[ii])*invdx2;
                    double h=H[ii];
                    double FH = lH - lamH*h*(h*h-vH*vH) - kapH*s*h;
                    Hn[ii]=h+dtau*FH;
                }
                /* 2nd gauge Gauss/Poisson: F_B = lap B + g2 sum_a q_a omega_a |phi_a|^2 */
                if(g2!=0.0){
                    double lB=(B[Xp]+B[Xm]+B[Yp]+B[Ym]+B[Zp]+B[Zm]-6.0*B[ii])*invdx2;
                    double src=0;
                    for(int a=0;a<3;a++){ double oma=omega+g2*qch[a]*B[ii];
                        src += qch[a]*oma*(u[a][ii]*u[a][ii]+v[a][ii]*v[a][ii]); }
                    Bn[ii]=B[ii]+dtau*(lB + g2*src);
                }
            }
            /* swap old<->new (Jacobi) */
            for(int a=0;a<3;a++){ double*t;
                t=u[a];u[a]=un[a];un[a]=t;  t=v[a];v[a]=vn[a];vn[a]=t;
                t=tu[a];tu[a]=tun[a];tun[a]=t; t=tv[a];tv[a]=tvn[a];tvn[a]=t; }
            if(vH>0){ double*t=H;H=Hn;Hn=t; }
            if(g2!=0.0){ double*t=B;B=Bn;Bn=t; }
            /* norm-preserving projection (hold charge ~ fixed) */
            double M=0;
            #pragma omp parallel for reduction(+:M)
            for(long ii=0;ii<N3;ii++) for(int a=0;a<3;a++) M+=u[a][ii]*u[a][ii]+v[a][ii]*v[a][ii];
            if(M>0){ double sc=sqrt(M0/M);
                #pragma omp parallel for
                for(long ii=0;ii<N3;ii++) for(int a=0;a<3;a++){ u[a][ii]*=sc; v[a][ii]*=sc; } }
            if(it==iters-1){
                /* report energy + s_max at end of level */
                double E=0,smax=0,Hmin=vH>0?vH:0;
                #pragma omp parallel for reduction(+:E) reduction(max:smax) reduction(min:Hmin)
                for(long ii=0;ii<N3;ii++){
                    double r0=u[0][ii]*u[0][ii]+v[0][ii]*v[0][ii],r1=u[1][ii]*u[1][ii]+v[1][ii]*v[1][ii],r2=u[2][ii]*u[2][ii]+v[2][ii]*v[2][ii];
                    double s=r0*r1*r2; if(s>smax)smax=s;
                    double e=0;
                    for(int a=0;a<3;a++){ e+=0.5*m2w*(u[a][ii]*u[a][ii]+v[a][ii]*v[a][ii])+0.5*mt2w*(tu[a][ii]*tu[a][ii]+tv[a][ii]*tv[a][ii]); }
                    e+=0.5*mu*s/(1+kappa*s);
                    if(vH>0){ double h=H[ii]; e+=0.25*lamH*(h*h-vH*vH)*(h*h-vH*vH)+0.5*kapH*s*h*h; if(h<Hmin)Hmin=h; }
                    E+=e;
                }
                printf("  %-6d %-10.4g %-12.5e %-12.4e %-12.4e  Hmin=%.4f\n",lev,P,E*dx*dx*dx,maxF,smax,Hmin);
                fflush(stdout);
            }
        }
    }

    /* write relaxed seed (rotating at omega) */
    float *cols[NCOLS]; for(int c=0;c<NCOLS;c++) cols[c]=calloc((size_t)N3,sizeof(float));
    #pragma omp parallel for
    for(long ii=0;ii<N3;ii++) for(int a=0;a<3;a++){
        cols[0+a][ii]=(float)u[a][ii];  cols[12+a][ii]=(float)v[a][ii];
        cols[6+a][ii]=(float)(-omega*v[a][ii]); cols[18+a][ii]=(float)(omega*u[a][ii]);
        cols[3+a][ii]=(float)tu[a][ii]; cols[15+a][ii]=(float)tv[a][ii];
        cols[9+a][ii]=(float)(-omega*tv[a][ii]); cols[21+a][ii]=(float)(omega*tu[a][ii]);
    }
    double dt=0.025*dx; SFA*sfa=sfa_create(out,N,N,N,L,L,L,dt);
    const char*nm[24]={"phi_x","phi_y","phi_z","theta_x","theta_y","theta_z","phi_vx","phi_vy","phi_vz","theta_vx","theta_vy","theta_vz","phiim_x","phiim_y","phiim_z","thetaim_x","thetaim_y","thetaim_z","phiim_vx","phiim_vy","phiim_vz","thetaim_vx","thetaim_vy","thetaim_vz"};
    int sem[24]={1,1,1,3,3,3,2,2,2,2,2,2,1,1,1,3,3,3,2,2,2,2,2,2};
    int cmp[24]={0,1,2,0,1,2,0,1,2,3,4,5,3,4,5,3,4,5,6,7,8,9,10,11};
    for(int c=0;c<24;c++) sfa_add_column(sfa,nm[c],SFA_F32,sem[c],cmp[c]);
    sfa_finalize_header(sfa); void*fc[NCOLS]; for(int c=0;c<NCOLS;c++) fc[c]=cols[c];
    sfa_write_frame(sfa,0.0,fc); sfa_close(sfa);
    printf("eta_relaxer: wrote %s\n",out);
    return 0;
}
