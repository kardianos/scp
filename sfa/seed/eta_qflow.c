/*  eta_qflow.c — FIXED-CHARGE gradient-flow relaxer for the stationary
 *  eta-coupled Q-ball (the correct variational principle behind the
 *  stability heuristics, CONCEPT.md §2: the conserved charge supplies the
 *  Q^2/N pressure that turns the Derrick saddle into a minimum).
 *
 *  Rotating ansatz  Phi_a = phi_a(x) e^{iwt},  Theta_a = theta_a(x) e^{iwt}.
 *  Conserved U(1) charge  Q = w * N,   N = int sum_a (|phi_a|^2 + |theta_a|^2).
 *  Eliminating w = Q/N gives the fixed-Q energy functional
 *
 *    E_Q = int sum_a [ 1/2|grad phi_a|^2 + 1/2 m^2|phi_a|^2
 *                    + 1/2|grad theta_a|^2 + 1/2 mth^2|theta_a|^2 ]
 *          + Vt(s) - eta * Re[ theta^bar . (curl phi) ]   +  Q^2/(2N)
 *
 *  Plain gradient flow on E_Q at fixed Q (w is NOT fixed, the norm is NOT
 *  rescaled — w = Q/N emerges as the Lagrange multiplier):
 *
 *    F_phi   = lap phi   - (m^2  -w^2) phi   - 2 Vt'(s) phi_a prod_rest
 *                                            + eta*curl(theta)
 *    F_theta = lap theta - (mth^2-w^2) theta + eta*curl(phi)
 *
 *  with w^2 = (Q/N)^2 recomputed every step. This is monotone descent on a
 *  single functional; its fixed points are EXACTLY the kernel's stationary
 *  rotating solutions on the same stencil (7-pt Laplacian, central curl —
 *  copied verbatim from scp_sim.c compute_forces_complex).
 *
 *  Predecessor eta_relaxer fixed BOTH w and the norm (over-constrained; the
 *  flow leaked charge and needed an external pressure cage — FUTURE.md
 *  2026-06-26). It also seeded Phi_a = f/sqrt(3): the radial_qball profile f
 *  is the PER-COMPONENT amplitude (gen_qball_boost convention), so that seed
 *  was off-shell by 27x in s. This tool seeds Phi_a = f directly.
 *
 *  Theta init "bvp": the transverse partner g(r)(rhat x u)_a, u=(1,1,1)/sqrt3,
 *  from  g'' + (2/r)g' - (2/r^2)g - (mth^2-w0^2)g = -sqrt(3)*eta*f'
 *  (curl of Phi = f*(1,1,1) is sqrt(3) f' (rhat x u); radial_eta_soliton's
 *  source lacked the sqrt(3) because of its f/sqrt3 normalization).
 *
 *  WINDING MODE (v73, the real-space-circulation "Q-ring"): wind = n > 0 seeds
 *  Phi_a = f(r) * (rho/sqrt(rho^2+wtor^2))^n * e^{i n phi_az}  (all components,
 *  same winding; s is phase-blind so binding survives; f -> 0 on the axis).
 *  The state carries a REAL circulating charge current J_phi = n|f|^2/rho and
 *  angular momentum L_z = n Q exactly (spin = winding of the circulation).
 *  The flow may preserve or unwind it — either outcome is the experiment.
 *
 *  FIXED-(Q,J) MODE (v73 PROCESS.md §5.2, the two-multiplier principle):
 *  fixJ=1 constrains total angular momentum J alongside Q. Rigidly spinning
 *  ansatz  Phi(x,t) = R(Om t) phi(R(-Om t)x) e^{iwt}  (R rotates space AND
 *  the component index — the Cosserat vector structure needs both), so
 *  Phidot = iw phi - Om Jhat phi with  Jhat = d_phi - zhat x  (orbital+spin).
 *  Two multipliers from the 2x2 solve each step:
 *      Q = w N - Om l ,   J = w l - Om K ,
 *      N = int sum|f|^2,  l = int [u.Jv - v.Ju] (+Theta),  K = int |Jf|^2,
 *  E_kin = (wQ - Om J)/2, and the flow force gains Coriolis + centrifugal:
 *      F_u += w^2 u - 2 w Om (Jhat v) - Om^2 Jhat(Jhat u)   (v: +2wOm Jhat u).
 *  On an exact circulation eigenstate (Jhat phi = in phi) the pair (w,Om) is
 *  gauge-degenerate (only w_eff = w - n Om = Q/N physical); the solve is
 *  Tikhonov-regularized there — the force stays well-defined because it
 *  depends only on w_eff on that manifold. This closes the standing-wave
 *  fission channel BY CONSTRUCTION: rotating->standing changes J.
 *  Stability along a branch: the 2x2 VK matrix d(Q,J)/d(w,Om).
 *
 *  Build: gcc -O3 -march=native -fopenmp -o eta_qflow eta_qflow.c -lzstd -lm
 *  Usage: eta_qflow N L profile omega0 eta m mth mu kappa out.sfa
 *                   [Q iters dtau thinit printev wind wtor twist fixJ Jtarget]
 *    Q       target total charge (<=0: use omega0 * N_init)     [0]
 *    iters   flow iterations                                    [20000]
 *    dtau    flow step (<=0: 0.04*dx^2)                         [0]
 *    thinit  bvp | zero                                         [bvp]
 *    printev diagnostics every this many iterations             [500]
 *    wind    azimuthal winding number n (0 = ball)              [0]
 *    wtor    axis-regularization width for the torus seed       [2.0]
 *    twist   poloidal winding m                                 [0]
 *    fixJ    1 = constrain (Q, J) with two multipliers          [0]
 *    Jtarget target J; 0 with wind>0 -> wind*Q (0 with wind=0 -> J=0) [0]
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
static inline int wrp(int i){ return (i+N)%N; }

int main(int argc,char**argv){
    if(argc<11){fprintf(stderr,"usage: %s N L profile omega0 eta m mth mu kappa out.sfa [Q iters dtau thinit printev]\n",argv[0]);return 1;}
    N=atoi(argv[1]); double L=atof(argv[2]); const char*prof=argv[3];
    double omega0=atof(argv[4]), eta=atof(argv[5]), m=atof(argv[6]), mth=atof(argv[7]);
    double mu=atof(argv[8]), kappa=atof(argv[9]); const char*out=argv[10];
    double Qarg  = (argc>11)?atof(argv[11]):0.0;
    int    iters = (argc>12)?atoi(argv[12]):20000;
    double dtau  = (argc>13)?atof(argv[13]):0.0;
    const char*thinit = (argc>14)?argv[14]:"bvp";
    int printev  = (argc>15)?atoi(argv[15]):500;
    int wind     = (argc>16)?atoi(argv[16]):0;
    double wtor  = (argc>17)?atof(argv[17]):2.0;
    int twist    = (argc>18)?atoi(argv[18]):0;   /* poloidal winding m: phase = n*phi + m*chi */
    int fixJ     = (argc>19)?atoi(argv[19]):0;   /* two-multiplier fixed-(Q,J) mode */
    double Jarg  = (argc>20)?atof(argv[20]):0.0; /* J target (0 with wind>0 -> wind*Q) */

    N3=(long)N*N*N; NN=(long)N*N;
    double dx=2.0*L/(N-1), dx2=dx*dx, dV=dx*dx*dx;
    if(dtau<=0) dtau=0.04*dx2;

    /* ---- load radial profile (r f): f = PER-COMPONENT amplitude ---- */
    int cap=8192,np=0; double *pr=malloc(cap*sizeof(double)),*pf=malloc(cap*sizeof(double));
    FILE*fp=fopen(prof,"r"); if(!fp){fprintf(stderr,"cannot open %s\n",prof);return 1;}
    char ln[512];
    while(fgets(ln,sizeof(ln),fp)){ if(ln[0]=='#'||ln[0]=='\n')continue; double a,b;
        if(sscanf(ln,"%lf %lf",&a,&b)==2){ if(np>=cap){cap*=2;pr=realloc(pr,cap*sizeof(double));pf=realloc(pf,cap*sizeof(double));} pr[np]=a;pf[np]=b;np++; } }
    fclose(fp);
    if(np<4){fprintf(stderr,"profile too short\n");return 1;}
    double dr=pr[1]-pr[0], rmax=pr[np-1];

    /* f'(r) and the BVP Theta partner g(r) (Thomas), source -sqrt(3) eta f' */
    double *df=malloc(np*sizeof(double)), *gg=calloc(np,sizeof(double));
    for(int i=1;i<np-1;i++) df[i]=(pf[i+1]-pf[i-1])/(2*dr);
    df[0]=(pf[1]-pf[0])/dr; df[np-1]=0;
    if(!strcmp(thinit,"bvp") && eta!=0.0){
        double k2=mth*mth-omega0*omega0;
        double *A=malloc(np*sizeof(double)),*B=malloc(np*sizeof(double)),*C=malloc(np*sizeof(double)),*D=malloc(np*sizeof(double));
        double *cp=malloc(np*sizeof(double)),*dp=malloc(np*sizeof(double));
        for(int i=1;i<np-1;i++){ double r=pr[i];
            A[i]=1.0/(dr*dr)-1.0/(r*dr);
            B[i]=-2.0/(dr*dr)-2.0/(r*r)-k2;
            C[i]=1.0/(dr*dr)+1.0/(r*dr);
            D[i]=-sqrt(3.0)*eta*df[i]; }
        int lo=1,hi=np-2;
        cp[lo]=C[lo]/B[lo]; dp[lo]=D[lo]/B[lo];
        for(int i=lo+1;i<=hi;i++){ double mm=B[i]-A[i]*cp[i-1]; cp[i]=C[i]/mm; dp[i]=(D[i]-A[i]*dp[i-1])/mm; }
        gg[hi]=dp[hi];
        for(int i=hi-1;i>=lo;i--) gg[i]=dp[i]-cp[i]*gg[i+1];
        free(A);free(B);free(C);free(D);free(cp);free(dp);
    }
    double gmax=0; for(int i=0;i<np;i++) if(fabs(gg[i])>gmax)gmax=fabs(gg[i]);

    /* ---- fields (12): u,v = Re,Im Phi; tu,tv = Re,Im Theta; Jacobi buffers */
    double *u[3],*v[3],*tu[3],*tv[3], *un[3],*vn[3],*tun[3],*tvn[3];
    for(int a=0;a<3;a++){ u[a]=calloc(N3,sizeof(double));v[a]=calloc(N3,sizeof(double));
                          tu[a]=calloc(N3,sizeof(double));tv[a]=calloc(N3,sizeof(double));
                          un[a]=calloc(N3,sizeof(double));vn[a]=calloc(N3,sizeof(double));
                          tun[a]=calloc(N3,sizeof(double));tvn[a]=calloc(N3,sizeof(double)); }
    double inv=1.0/sqrt(3.0), ux=inv,uy=inv,uz=inv;
    #pragma omp parallel for
    for(long ii=0;ii<N3;ii++){
        int i=(int)(ii/NN),j=(int)((ii/N)%N),k=(int)(ii%N);
        double x=-L+i*dx,y=-L+j*dx,z=-L+k*dx,r=sqrt(x*x+y*y+z*z);
        double f=0,g=0;
        if(r<rmax){ double t=r/dr; int id=(int)t; if(id>=np-1)id=np-2; double fr=t-id;
            f=pf[id]+fr*(pf[id+1]-pf[id]); g=gg[id]+fr*(gg[id+1]-gg[id]); }
        if(wind>0 || twist!=0){
            /* Q-ring ("smoke ring"): a tube whose cross-section is the BALL
             * profile f(d), d = dist to the circle rho=wtor (major radius),
             * so the core sits at full binding amplitude; axis regularizer
             * (rho/sqrt(rho^2+1))^n keeps the winding phase regular at rho=0.
             * TWISTED ring (m = twist): phase gains m*chi (chi = poloidal
             * angle around the tube core); the tube must then be HOLLOW —
             * regularizer (d/sqrt(d^2+1))^|m| zeroes the amplitude on the
             * core circle. Zero set = z-axis (n) + core circle (m): a Hopf
             * link; circulation lines are (n,m) torus knots. */
            double rho=sqrt(x*x+y*y), phi_az=atan2(y,x);
            double drho=rho-wtor, d=sqrt(drho*drho+z*z);
            double chi=atan2(z,drho);
            double fd=0;
            if(d<rmax){ double t=d/dr; int id=(int)t; if(id>=np-1)id=np-2; double frq=t-id;
                fd=pf[id]+frq*(pf[id+1]-pf[id]); }
            double amp=fd;
            for(int q=0;q<wind;q++) amp*=rho/sqrt(rho*rho+1.0);
            for(int q=0;q<abs(twist);q++) amp*=d/sqrt(d*d+1.0);
            double ph=wind*phi_az + twist*chi;
            double cw=cos(ph), sw=sin(ph);
            for(int a=0;a<3;a++){ u[a][ii]=amp*cw; v[a][ii]=amp*sw; }
            continue;                                  /* Theta seeded 0 in winding mode */
        }
        for(int a=0;a<3;a++) u[a][ii]=f;              /* per-component amplitude */
        double rx=(r>1e-9)?x/r:0, ry=(r>1e-9)?y/r:0, rz=(r>1e-9)?z/r:0;
        tu[0][ii]=g*(ry*uz-rz*uy); tu[1][ii]=g*(rz*ux-rx*uz); tu[2][ii]=g*(rx*uy-ry*ux);
    }

    /* ---- fixed charge ---- */
    double Nnorm=0;
    #pragma omp parallel for reduction(+:Nnorm)
    for(long ii=0;ii<N3;ii++) for(int a=0;a<3;a++)
        Nnorm+=u[a][ii]*u[a][ii]+v[a][ii]*v[a][ii]+tu[a][ii]*tu[a][ii]+tv[a][ii]*tv[a][ii];
    Nnorm*=dV;
    double Q = (Qarg>0)? Qarg : omega0*Nnorm;
    double Jt = (fixJ && Jarg==0.0 && wind>0)? (double)wind*Q : Jarg;
    printf("# eta_qflow N=%d L=%g dx=%.4f eta=%.3f m=%.3f mth=%.3f mu=%g kappa=%g\n",N,L,dx,eta,m,mth,mu,kappa);
    printf("# thinit=%s g_max=%.4e  N_init=%.4f  Q=%.4f  omega_init=Q/N=%.6f  dtau=%.5g iters=%d\n",
           thinit,gmax,Nnorm,Q,Q/Nnorm,dtau,iters);
    if(fixJ){
        printf("# fixJ: J_target=%.4f  (two multipliers w,Om; Jhat = d_phi - zhat x)\n",Jt);
        printf("# %-8s %-10s %-10s %-12s %-12s %-12s %-12s %-10s %-10s\n",
               "iter","omega","Omega","E_QJ","maxF","s_max","N","J_act","cond");
    } else
    printf("# %-8s %-10s %-12s %-12s %-12s %-12s\n","iter","omega","E_Q","maxF","s_max","N");

    /* fixed-(Q,J) work arrays: Jhat applied to all 12 fields.
     * wmask: smooth cylindrical cutoff on the rotational force terms — a
     * rigidly rotating frame is only defined inside the light cylinder
     * rho < 1/|Om|; beyond it the Om^2 Jhat^2 term amplifies vacuum modes.
     * Fields decay as e^{-m rho} so masking at rho_c=11 biases nothing
     * (|phi| ~ 3e-5 there), it only kills the far-corner instability.
     * |Omega| is additionally capped at 0.15 (light cylinder at rho=6.7 >
     * any tested structure would be unphysical anyway). */
    double *ju[3],*jv[3],*jtu[3],*jtv[3], *wmask=NULL;
    const double OMEGA_CAP=0.15;
    if(fixJ){
        for(int a=0;a<3;a++){
            ju[a]=calloc(N3,sizeof(double));  jv[a]=calloc(N3,sizeof(double));
            jtu[a]=calloc(N3,sizeof(double)); jtv[a]=calloc(N3,sizeof(double)); }
        wmask=malloc(N3*sizeof(double));
        for(long ii=0;ii<N3;ii++){
            int i=(int)(ii/NN),j=(int)((ii/N)%N);
            double x=-L+i*dx, y=-L+j*dx, rho=sqrt(x*x+y*y);
            wmask[ii]=1.0/(1.0+exp((rho-11.0)/0.5));
        }
    }

    double invdx2=1.0/dx2, inv2dx=1.0/(2.0*dx);
    double omega=Q/Nnorm, Omega=0.0, ell=0.0, Kk=0.0, cond=1.0;
    if(fixJ && (argc<=13 || atof(argv[13])<=0)) dtau=0.02*dx2;  /* Om^2 Jhat^2 stiffness margin */

    /* (Jhat w)_a = x d_y w_a - y d_x w_a - (zhat x w)_a ; spin part:
     * (zhat x w)_0 = -w_1, (zhat x w)_1 = +w_0, (zhat x w)_2 = 0 */
    #define JSPIN(w,a,ii) ((a)==0? (w)[1][ii] : (a)==1? -(w)[0][ii] : 0.0)

    for(int it=0; it<=iters; it++){
        double Nn=0;
        if(fixJ){
            /* two-multiplier solve: compute Jhat of all 12 fields + moments */
            double el=0, kk=0;
            #pragma omp parallel for reduction(+:Nn,el,kk)
            for(long ii=0;ii<N3;ii++){
                int i=(int)(ii/NN),j=(int)((ii/N)%N),k=(int)(ii%N);
                int ip=wrp(i+1),im=wrp(i-1),jp=wrp(j+1),jm=wrp(j-1);
                long Xp=IX(ip,j,k),Xm=IX(im,j,k),Yp=IX(i,jp,k),Ym=IX(i,jm,k);
                double x=-L+i*dx, y=-L+j*dx;
                for(int a=0;a<3;a++){
                    double jua = x*(u[a][Yp]-u[a][Ym])*inv2dx - y*(u[a][Xp]-u[a][Xm])*inv2dx + JSPIN(u,a,ii);
                    double jva = x*(v[a][Yp]-v[a][Ym])*inv2dx - y*(v[a][Xp]-v[a][Xm])*inv2dx + JSPIN(v,a,ii);
                    double jta = x*(tu[a][Yp]-tu[a][Ym])*inv2dx - y*(tu[a][Xp]-tu[a][Xm])*inv2dx + JSPIN(tu,a,ii);
                    double jsa = x*(tv[a][Yp]-tv[a][Ym])*inv2dx - y*(tv[a][Xp]-tv[a][Xm])*inv2dx + JSPIN(tv,a,ii);
                    ju[a][ii]=jua; jv[a][ii]=jva; jtu[a][ii]=jta; jtv[a][ii]=jsa;
                    Nn += u[a][ii]*u[a][ii]+v[a][ii]*v[a][ii]+tu[a][ii]*tu[a][ii]+tv[a][ii]*tv[a][ii];
                    el += u[a][ii]*jva - v[a][ii]*jua + tu[a][ii]*jsa - tv[a][ii]*jta;
                    kk += jua*jua + jva*jva + jta*jta + jsa*jsa;
                }
            }
            Nn*=dV; ell=el*dV; Kk=kk*dV;
            double D = Nn*Kk - ell*ell, Dr = D + 1e-10*Nn*Kk;
            cond = D/(Nn*Kk);
            omega = (Kk*Q - ell*Jt)/Dr;
            Omega = (ell*Q - Nn*Jt)/Dr;
            if(Omega> OMEGA_CAP){Omega= OMEGA_CAP; omega=(Q+Omega*ell)/Nn;}
            if(Omega<-OMEGA_CAP){Omega=-OMEGA_CAP; omega=(Q+Omega*ell)/Nn;}
        } else {
            /* Lagrange multiplier from the charge constraint */
            #pragma omp parallel for reduction(+:Nn)
            for(long ii=0;ii<N3;ii++) for(int a=0;a<3;a++)
                Nn+=u[a][ii]*u[a][ii]+v[a][ii]*v[a][ii]+tu[a][ii]*tu[a][ii]+tv[a][ii]*tv[a][ii];
            Nn*=dV; omega=Q/Nn;
        }
        double m2w=m*m-omega*omega, mt2w=mth*mth-omega*omega;

        double maxF=0;
        #pragma omp parallel for reduction(max:maxF)
        for(long ii=0;ii<N3;ii++){
            int i=(int)(ii/NN),j=(int)((ii/N)%N),k=(int)(ii%N);
            int ip=wrp(i+1),im=wrp(i-1),jp=wrp(j+1),jm=wrp(j-1),kp=wrp(k+1),km=wrp(k-1);
            long Xp=IX(ip,j,k),Xm=IX(im,j,k),Yp=IX(i,jp,k),Ym=IX(i,jm,k),Zp=IX(i,j,kp),Zm=IX(i,j,km);
            double r0=u[0][ii]*u[0][ii]+v[0][ii]*v[0][ii];
            double r1=u[1][ii]*u[1][ii]+v[1][ii]*v[1][ii];
            double r2=u[2][ii]*u[2][ii]+v[2][ii]*v[2][ii];
            double s=r0*r1*r2, den=1.0+kappa*s;
            double Vtp=0.5*mu/(den*den);
            double Pa[3]={r1*r2,r0*r2,r0*r1};
            for(int a=0;a<3;a++){
                double lu =(u[a][Xp]+u[a][Xm]+u[a][Yp]+u[a][Ym]+u[a][Zp]+u[a][Zm]-6.0*u[a][ii])*invdx2;
                double lv =(v[a][Xp]+v[a][Xm]+v[a][Yp]+v[a][Ym]+v[a][Zp]+v[a][Zm]-6.0*v[a][ii])*invdx2;
                double ltu=(tu[a][Xp]+tu[a][Xm]+tu[a][Yp]+tu[a][Ym]+tu[a][Zp]+tu[a][Zm]-6.0*tu[a][ii])*invdx2;
                double ltv=(tv[a][Xp]+tv[a][Xm]+tv[a][Yp]+tv[a][Ym]+tv[a][Zp]+tv[a][Zm]-6.0*tv[a][ii])*invdx2;
                int b=(a+1)%3, c=(a+2)%3;
                #define DERIV(arr,ax,comp) ( (ax==0)?(arr[comp][Xp]-arr[comp][Xm]):(ax==1)?(arr[comp][Yp]-arr[comp][Ym]):(arr[comp][Zp]-arr[comp][Zm]) )
                double curl_tu = (DERIV(tu,b,c) - DERIV(tu,c,b))*inv2dx;
                double curl_tv = (DERIV(tv,b,c) - DERIV(tv,c,b))*inv2dx;
                double curl_u  = (DERIV(u, b,c) - DERIV(u, c,b))*inv2dx;
                double curl_v  = (DERIV(v, b,c) - DERIV(v, c,b))*inv2dx;
                #undef DERIV
                double Fu = lu  - m2w*u[a][ii]  - 2.0*Vtp*u[a][ii]*Pa[a] + eta*curl_tu;
                double Fv = lv  - m2w*v[a][ii]  - 2.0*Vtp*v[a][ii]*Pa[a] + eta*curl_tv;
                double Ftu= ltu - mt2w*tu[a][ii] + eta*curl_u;
                double Ftv= ltv - mt2w*tv[a][ii] + eta*curl_v;
                if(fixJ){
                    /* Coriolis + centrifugal (PROCESS.md §5.2), light-cylinder masked:
                     * F_u += -2wOm (Jhat v) - Om^2 Jhat(Jhat u), v: +2wOm (Jhat u) */
                    double x=-L+i*dx, y=-L+j*dx, wm=wmask[ii];
                    double wO=2.0*omega*Omega*wm, O2=Omega*Omega*wm;
                    double j2u = x*(ju[a][Yp]-ju[a][Ym])*inv2dx - y*(ju[a][Xp]-ju[a][Xm])*inv2dx + JSPIN(ju,a,ii);
                    double j2v = x*(jv[a][Yp]-jv[a][Ym])*inv2dx - y*(jv[a][Xp]-jv[a][Xm])*inv2dx + JSPIN(jv,a,ii);
                    double j2t = x*(jtu[a][Yp]-jtu[a][Ym])*inv2dx - y*(jtu[a][Xp]-jtu[a][Xm])*inv2dx + JSPIN(jtu,a,ii);
                    double j2s = x*(jtv[a][Yp]-jtv[a][Ym])*inv2dx - y*(jtv[a][Xp]-jtv[a][Xm])*inv2dx + JSPIN(jtv,a,ii);
                    Fu  += -wO*jv[a][ii]  - O2*j2u;
                    Fv  += +wO*ju[a][ii]  - O2*j2v;
                    Ftu += -wO*jtv[a][ii] - O2*j2t;
                    Ftv += +wO*jtu[a][ii] - O2*j2s;
                }
                un[a][ii]=u[a][ii]+dtau*Fu;   vn[a][ii]=v[a][ii]+dtau*Fv;
                tun[a][ii]=tu[a][ii]+dtau*Ftu; tvn[a][ii]=tv[a][ii]+dtau*Ftv;
                double af=fabs(Fu); if(af>maxF)maxF=af;
                af=fabs(Ftu); if(af>maxF)maxF=af;
            }
        }
        for(int a=0;a<3;a++){ double*t;
            t=u[a];u[a]=un[a];un[a]=t;  t=v[a];v[a]=vn[a];vn[a]=t;
            t=tu[a];tu[a]=tun[a];tun[a]=t; t=tv[a];tv[a]=tvn[a];tvn[a]=t; }

        if(it%printev==0 || it==iters){
            /* full E_Q: grad + mass + Vt + eta-cross + Q^2/(2N) */
            double Eg=0,Em=0,Ep=0,Ex=0,smax=0;
            #pragma omp parallel for reduction(+:Eg,Em,Ep,Ex) reduction(max:smax)
            for(long ii=0;ii<N3;ii++){
                int i=(int)(ii/NN),j=(int)((ii/N)%N),k=(int)(ii%N);
                int ip=wrp(i+1),im=wrp(i-1),jp=wrp(j+1),jm=wrp(j-1),kp=wrp(k+1),km=wrp(k-1);
                long Xp=IX(ip,j,k),Xm=IX(im,j,k),Yp=IX(i,jp,k),Ym=IX(i,jm,k),Zp=IX(i,j,kp),Zm=IX(i,j,km);
                double r0=u[0][ii]*u[0][ii]+v[0][ii]*v[0][ii];
                double r1=u[1][ii]*u[1][ii]+v[1][ii]*v[1][ii];
                double r2=u[2][ii]*u[2][ii]+v[2][ii]*v[2][ii];
                double s=r0*r1*r2; if(s>smax)smax=s;
                Ep += 0.5*mu*s/(1.0+kappa*s);
                for(int a=0;a<3;a++){
                    double gx,gy,gz;
                    gx=(u[a][Xp]-u[a][Xm])*inv2dx; gy=(u[a][Yp]-u[a][Ym])*inv2dx; gz=(u[a][Zp]-u[a][Zm])*inv2dx;
                    Eg+=0.5*(gx*gx+gy*gy+gz*gz);
                    gx=(v[a][Xp]-v[a][Xm])*inv2dx; gy=(v[a][Yp]-v[a][Ym])*inv2dx; gz=(v[a][Zp]-v[a][Zm])*inv2dx;
                    Eg+=0.5*(gx*gx+gy*gy+gz*gz);
                    gx=(tu[a][Xp]-tu[a][Xm])*inv2dx; gy=(tu[a][Yp]-tu[a][Ym])*inv2dx; gz=(tu[a][Zp]-tu[a][Zm])*inv2dx;
                    Eg+=0.5*(gx*gx+gy*gy+gz*gz);
                    gx=(tv[a][Xp]-tv[a][Xm])*inv2dx; gy=(tv[a][Yp]-tv[a][Ym])*inv2dx; gz=(tv[a][Zp]-tv[a][Zm])*inv2dx;
                    Eg+=0.5*(gx*gx+gy*gy+gz*gz);
                    Em+=0.5*m*m*(u[a][ii]*u[a][ii]+v[a][ii]*v[a][ii])
                       +0.5*mth*mth*(tu[a][ii]*tu[a][ii]+tv[a][ii]*tv[a][ii]);
                    int b=(a+1)%3, c=(a+2)%3;
                    #define DERIV(arr,ax,comp) ( (ax==0)?(arr[comp][Xp]-arr[comp][Xm]):(ax==1)?(arr[comp][Yp]-arr[comp][Ym]):(arr[comp][Zp]-arr[comp][Zm]) )
                    double curl_u  = (DERIV(u, b,c) - DERIV(u, c,b))*inv2dx;
                    double curl_v  = (DERIV(v, b,c) - DERIV(v, c,b))*inv2dx;
                    #undef DERIV
                    Ex += -eta*(tu[a][ii]*curl_u + tv[a][ii]*curl_v);
                }
            }
            if(fixJ){
                double EQJ=(Eg+Em+Ep+Ex)*dV + 0.5*(omega*Q - Omega*Jt);
                double Jact=omega*ell - Omega*Kk;
                printf("  %-8d %-10.6f %-10.6f %-12.6e %-12.4e %-12.4e %-12.4f %-10.4f %-10.3e\n",
                       it,omega,Omega,EQJ,maxF,smax,Nn,Jact,cond);
            } else {
            double EQ=(Eg+Em+Ep+Ex)*dV + Q*Q/(2.0*Nn);
            printf("  %-8d %-10.6f %-12.6e %-12.4e %-12.4e %-12.4f\n",it,omega,EQ,maxF,smax,Nn);
            }
            fflush(stdout);
        }
    }

    /* fixed-(Q,J): refresh Jhat arrays for the FINAL fields (the last update
     * ran after the last Jhat pass) — needed for velocities and the report */
    if(fixJ){
        double el=0, kk=0;
        #pragma omp parallel for reduction(+:el,kk)
        for(long ii=0;ii<N3;ii++){
            int i=(int)(ii/NN),j=(int)((ii/N)%N),k=(int)(ii%N);
            int ip=wrp(i+1),im=wrp(i-1),jp=wrp(j+1),jm=wrp(j-1);
            long Xp=IX(ip,j,k),Xm=IX(im,j,k),Yp=IX(i,jp,k),Ym=IX(i,jm,k);
            double x=-L+i*dx, y=-L+j*dx;
            for(int a=0;a<3;a++){
                double jua = x*(u[a][Yp]-u[a][Ym])*inv2dx - y*(u[a][Xp]-u[a][Xm])*inv2dx + JSPIN(u,a,ii);
                double jva = x*(v[a][Yp]-v[a][Ym])*inv2dx - y*(v[a][Xp]-v[a][Xm])*inv2dx + JSPIN(v,a,ii);
                double jta = x*(tu[a][Yp]-tu[a][Ym])*inv2dx - y*(tu[a][Xp]-tu[a][Xm])*inv2dx + JSPIN(tu,a,ii);
                double jsa = x*(tv[a][Yp]-tv[a][Ym])*inv2dx - y*(tv[a][Xp]-tv[a][Xm])*inv2dx + JSPIN(tv,a,ii);
                ju[a][ii]=jua; jv[a][ii]=jva; jtu[a][ii]=jta; jtv[a][ii]=jsa;
                el += u[a][ii]*jva - v[a][ii]*jua + tu[a][ii]*jsa - tv[a][ii]*jta;
                kk += jua*jua + jva*jva + jta*jta + jsa*jsa;
            }
        }
        ell=el*dV; Kk=kk*dV;
    }

    /* final charges */
    double Nphi=0,Nth=0;
    #pragma omp parallel for reduction(+:Nphi,Nth)
    for(long ii=0;ii<N3;ii++) for(int a=0;a<3;a++){
        Nphi+=u[a][ii]*u[a][ii]+v[a][ii]*v[a][ii];
        Nth +=tu[a][ii]*tu[a][ii]+tv[a][ii]*tv[a][ii]; }
    Nphi*=dV; Nth*=dV;
    printf("# final: omega=%.6f  Q_phi=%.4f  Q_theta=%.4f  Q_total=%.4f (target %.4f)\n",
           omega, omega*Nphi, omega*Nth, omega*(Nphi+Nth), Q);
    if(fixJ){
        double Nfin=Nphi+Nth;
        printf("# fixJ final: Omega=%.6f  J_actual=%.4f (target %.4f)  Q_actual=%.4f  cond=%.3e  w_eff(Q/N)=%.6f\n",
               Omega, omega*ell-Omega*Kk, Jt, omega*Nfin-Omega*ell, cond, Q/Nfin);
    }
    if(omega>=m)      printf("# WARNING: omega >= m (outside existence window)\n");
    if(omega*omega>=mth*mth) printf("# WARNING: omega >= m_theta (Theta unbound)\n");
    if(wind>0 || twist!=0){
        /* L_z = omega * int sum_a Im[phibar d_phi phi]; d_phi = x d_y - y d_x.
         * Also: winding number around the circle of peak amplitude (z=0), and
         * the peak radius rho_pk (ring geometry check). */
        double Limz=0; double rho_pk=0, amp_pk=0;
        int kz=N/2;
        for(int i=1;i<N-1;i++){ double x=-L+i*dx;
            for(int j=1;j<N-1;j++){ double y=-L+j*dx;
                long ii=IX(i,j,kz);
                double a2=0; for(int a=0;a<3;a++) a2+=u[a][ii]*u[a][ii]+v[a][ii]*v[a][ii];
                double rho=sqrt(x*x+y*y);
                if(a2>amp_pk){ amp_pk=a2; rho_pk=rho; } } }
        #pragma omp parallel for reduction(+:Limz)
        for(long ii=0;ii<N3;ii++){
            int i=(int)(ii/NN),j=(int)((ii/N)%N),k=(int)(ii%N);
            if(i<1||i>=N-1||j<1||j>=N-1||k<1||k>=N-1) continue;
            double x=-L+i*dx,y=-L+j*dx;
            long Xp=IX(i+1,j,k),Xm=IX(i-1,j,k),Yp=IX(i,j+1,k),Ym=IX(i,j-1,k);
            for(int a=0;a<3;a++){
                double dyu=(u[a][Yp]-u[a][Ym])/(2*dx), dyv=(v[a][Yp]-v[a][Ym])/(2*dx);
                double dxu=(u[a][Xp]-u[a][Xm])/(2*dx), dxv=(v[a][Xp]-v[a][Xm])/(2*dx);
                double dphiu=x*dyu-y*dxu, dphiv=x*dyv-y*dxv;
                Limz += u[a][ii]*dphiv - v[a][ii]*dphiu;
            }
        }
        double Lz=omega*Limz*dV;
        /* winding of component 0's phase around the peak circle */
        int M=256; double wsum=0, prev=0;
        for(int mth_i=0;mth_i<=M;mth_i++){
            double ph=2.0*M_PI*mth_i/M;
            double x=rho_pk*cos(ph), y=rho_pk*sin(ph);
            int i=(int)lround((x+L)/dx), j=(int)lround((y+L)/dx);
            if(i<0||i>=N||j<0||j>=N) continue;
            long ii=IX(i,j,kz);
            double arg=atan2(v[0][ii],u[0][ii]);
            if(mth_i>0){ double d=arg-prev;
                while(d>M_PI)d-=2*M_PI; while(d<-M_PI)d+=2*M_PI; wsum+=d; }
            prev=arg;
        }
        printf("# ring: L_z=%.4f  L_z/Q=%.4f (target wind=%d)  rho_peak=%.3f  winding@peak=%.2f\n",
               Lz, Lz/Q, wind, rho_pk, wsum/(2*M_PI));
    }

    /* write rotating seed, gen_qball_boost column order.
     * fixJ: Phidot = iw phi - Om Jhat phi -> udot = -w v - Om (Jhat u), etc. */
    float *cols[NCOLS]; for(int c=0;c<NCOLS;c++) cols[c]=calloc((size_t)N3,sizeof(float));
    #pragma omp parallel for
    for(long ii=0;ii<N3;ii++) for(int a=0;a<3;a++){
        double wm =fixJ? wmask[ii]:0;
        double Ju =fixJ? wm*ju[a][ii] :0, Jvv=fixJ? wm*jv[a][ii] :0;
        double Jtu=fixJ? wm*jtu[a][ii]:0, Jtv=fixJ? wm*jtv[a][ii]:0;
        cols[0+a][ii]=(float)u[a][ii];   cols[12+a][ii]=(float)v[a][ii];
        cols[6+a][ii]=(float)(-omega*v[a][ii]-Omega*Ju);  cols[18+a][ii]=(float)(omega*u[a][ii]-Omega*Jvv);
        cols[3+a][ii]=(float)tu[a][ii];  cols[15+a][ii]=(float)tv[a][ii];
        cols[9+a][ii]=(float)(-omega*tv[a][ii]-Omega*Jtu); cols[21+a][ii]=(float)(omega*tu[a][ii]-Omega*Jtv);
    }
    double dt=0.025*dx; SFA*sfa=sfa_create(out,N,N,N,L,L,L,dt);
    const char*nm[24]={"phi_x","phi_y","phi_z","theta_x","theta_y","theta_z",
        "phi_vx","phi_vy","phi_vz","theta_vx","theta_vy","theta_vz",
        "phiim_x","phiim_y","phiim_z","thetaim_x","thetaim_y","thetaim_z",
        "phiim_vx","phiim_vy","phiim_vz","thetaim_vx","thetaim_vy","thetaim_vz"};
    /* semantics: SFA_POSITION=0 (phi), SFA_ANGLE=1 (theta), SFA_VELOCITY=2.
     * NOTE: eta_relaxer shipped with {1,3,2} here — phi landed in the kernel's
     * theta and theta was tagged ACCELERATION; every relaxer seed was scrambled
     * on load (the June-26 "-8.1% and dies" release run). */
    int sem[24]={0,0,0,1,1,1,2,2,2,2,2,2,0,0,0,1,1,1,2,2,2,2,2,2};
    int cmp[24]={0,1,2,0,1,2,0,1,2,3,4,5,3,4,5,3,4,5,6,7,8,9,10,11};
    for(int c=0;c<24;c++) sfa_add_column(sfa,nm[c],SFA_F32,sem[c],cmp[c]);
    sfa_finalize_header(sfa); void*fc[NCOLS]; for(int c=0;c<NCOLS;c++) fc[c]=cols[c];
    sfa_write_frame(sfa,0.0,fc); sfa_close(sfa);
    printf("eta_qflow: wrote %s\n",out);
    return 0;
}
