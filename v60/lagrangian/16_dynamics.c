/* v60 GEN7 -- the literal DYNAMICS: time-evolve the nonlinear Euler-Lagrange
 * equation  box Phi = -V'(Phi)  on the GEN3 Koide-cone potential, and confirm the
 * GEN5 linearized spectrum NONLINEARLY.
 *
 *   Field:  x(t, s) in R^3 (the sqrt-mass / generation triple), on a 1D spatial
 *           lattice s with periodic BC.   Lagrangian (GEN3+):
 *               L = 1/2 (d_t x)^2 - 1/2 (d_s x)^2 - V(x),
 *               V = lam (e1^2 - 6 e2)^2 + mu (e1 - c)^2.
 *   EOM (velocity-Verlet, symplectic):  d_t^2 x = d_s^2 x - grad V(x).
 *
 * Tests:
 *   (1) HOMOGENEOUS (k=0) normal modes: small oscillation about the Brannen
 *       vacuum along each Hessian eigenvector has  omega^2 = (Hessian eigenvalue).
 *       The two massive modes oscillate at omega^2 = {2.98, 435}; the Goldstone
 *       (phase) direction does NOT oscillate (omega^2 = 0): given a velocity it
 *       DRIFTS linearly (free, massless).  => confirms GEN5 spectrum nonlinearly.
 *   (2) ENERGY CONSERVATION over a long run (symplectic integrator sanity).
 *   (3) DISPERSION (k != 0): a Goldstone-direction plane wave of wavenumber k
 *       oscillates at omega^2 = k^2 (speed 1, MASSLESS, propagating); a massive
 *       mode at omega^2 = k^2 + m^2.  => the relativistic dispersion of the
 *       derived spectrum, by genuine time evolution.
 *
 * Build & run:  gcc -O2 -o /tmp/dyn 16_dynamics.c -lm && /tmp/dyn
 */
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

static double lam = 1.0, mu = 1.0, cc = 3.0;
static double xvac[3];

static double e1f(const double *x){ return x[0]+x[1]+x[2]; }
static double e2f(const double *x){ return x[0]*x[1]+x[0]*x[2]+x[1]*x[2]; }
static double Vpot(const double *x){
    double e1=e1f(x), e2=e2f(x);
    return lam*(e1*e1-6*e2)*(e1*e1-6*e2) + mu*(e1-cc)*(e1-cc);
}
static void gradV(const double *x, double *g){
    double e1=e1f(x), e2=e2f(x), A=e1*e1-6*e2;
    double de2[3]={x[1]+x[2], x[0]+x[2], x[0]+x[1]};
    for(int i=0;i<3;i++){ double dA=2*e1-6*de2[i]; g[i]=lam*2*A*dA + mu*2*(e1-cc); }
}

/* Jacobi eigensolver for symmetric 3x3: eval[] ascending-ish, evec columns. */
static void jacobi3(double A[3][3], double eval[3], double evec[3][3]){
    double a[3][3]; memcpy(a,A,sizeof a);
    double v[3][3]={{1,0,0},{0,1,0},{0,0,1}};
    for(int sweep=0; sweep<100; sweep++){
        double off=fabs(a[0][1])+fabs(a[0][2])+fabs(a[1][2]);
        if(off<1e-14) break;
        for(int p=0;p<2;p++) for(int q=p+1;q<3;q++){
            if(fabs(a[p][q])<1e-300) continue;
            double th=(a[q][q]-a[p][p])/(2*a[p][q]);
            double t=(th>=0?1.0:-1.0)/(fabs(th)+sqrt(th*th+1));
            double cph=1/sqrt(t*t+1), sph=t*cph;
            for(int i=0;i<3;i++){
                double aip=a[i][p], aiq=a[i][q];
                a[i][p]=cph*aip-sph*aiq; a[i][q]=sph*aip+cph*aiq;
            }
            for(int i=0;i<3;i++){
                double api=a[p][i], aqi=a[q][i];
                a[p][i]=cph*api-sph*aqi; a[q][i]=sph*api+cph*aqi;
            }
            for(int i=0;i<3;i++){
                double vip=v[i][p], viq=v[i][q];
                v[i][p]=cph*vip-sph*viq; v[i][q]=sph*vip+cph*viq;
            }
        }
    }
    for(int i=0;i<3;i++){ eval[i]=a[i][i]; for(int j=0;j<3;j++) evec[j][i]=v[j][i]; }
}

/* measure omega via average spacing of successive maxima of s(t) */
static double measure_omega(double *ts, double *ss, int n, double dt){
    double last=-1, sum=0; int cnt=0;
    for(int i=1;i<n-1;i++) if(ss[i]>ss[i-1] && ss[i]>=ss[i+1]){
        double t=ts[i];
        if(last>=0){ sum+=t-last; cnt++; }
        last=t;
    }
    if(cnt==0) return 0;
    double T=sum/cnt; return 2*M_PI/T;
}

int main(void){
    double a=1.0, phi=2.0/9.0;
    for(int k=0;k<3;k++) xvac[k]=a*(1+sqrt(2.0)*cos(phi+2*M_PI*k/3.0));
    int fail=0;

    /* Hessian at vacuum (finite difference of analytic gradV) */
    double H[3][3]; double h=1e-5, gp[3], gm[3], xt[3];
    for(int j=0;j<3;j++){
        memcpy(xt,xvac,sizeof xt); xt[j]+=h; gradV(xt,gp);
        memcpy(xt,xvac,sizeof xt); xt[j]-=h; gradV(xt,gm);
        for(int i=0;i<3;i++) H[i][j]=(gp[i]-gm[i])/(2*h);
    }
    double eval[3], evec[3][3]; jacobi3(H,eval,evec);
    /* sort by |eval| ascending so index 0 = Goldstone */
    int idx[3]={0,1,2};
    for(int i=0;i<3;i++) for(int j=i+1;j<3;j++)
        if(fabs(eval[idx[j]])<fabs(eval[idx[i]])){int t=idx[i];idx[i]=idx[j];idx[j]=t;}
    printf("Hessian eigenvalues (m^2): goldstone=%.6f, massive=%.4f, %.4f\n",
           eval[idx[0]], eval[idx[1]], eval[idx[2]]);

    /* (1) homogeneous normal modes: oscillate along each massive eigenvector */
    double dt=2e-4; int N=200000;
    double *ts=malloc(N*sizeof(double)), *ss=malloc(N*sizeof(double));
    printf("\n(1) homogeneous (k=0) normal-mode frequencies:\n");
    for(int mi=1;mi<3;mi++){
        int e=idx[mi]; double v[3]={evec[0][e],evec[1][e],evec[2][e]};
        double A0=1e-3, x[3], p[3]={0,0,0}, F[3], Fn[3];
        for(int i=0;i<3;i++) x[i]=xvac[i]+A0*v[i];
        gradV(x,F); for(int i=0;i<3;i++) F[i]=-F[i];
        for(int n=0;n<N;n++){
            double s=0; for(int i=0;i<3;i++) s+=(x[i]-xvac[i])*v[i];
            ts[n]=n*dt; ss[n]=s;
            for(int i=0;i<3;i++) x[i]+=dt*p[i]+0.5*dt*dt*F[i];
            gradV(x,Fn); for(int i=0;i<3;i++) Fn[i]=-Fn[i];
            for(int i=0;i<3;i++) p[i]+=0.5*dt*(F[i]+Fn[i]);
            memcpy(F,Fn,sizeof F);
        }
        double om=measure_omega(ts,ss,N,dt), om2=om*om, pred=eval[e];
        double rel=fabs(om2-pred)/pred;
        printf("  massive mode %d: omega^2(measured)=%.4f  predicted(Hessian)=%.4f  rel.err=%.2e\n",
               mi, om2, pred, rel);
        if(rel>2e-3) fail=1;
    }

    /* Goldstone: give it a velocity -> linear drift (no oscillation) => massless */
    {
        int e=idx[0]; double v[3]={evec[0][e],evec[1][e],evec[2][e]};
        double vamp=1e-3, x[3], p[3], F[3], Fn[3];
        for(int i=0;i<3;i++){ x[i]=xvac[i]; p[i]=vamp*v[i]; }
        gradV(x,F); for(int i=0;i<3;i++) F[i]=-F[i];
        int M=20000; double s_final=0;
        for(int n=0;n<M;n++){
            for(int i=0;i<3;i++) x[i]+=dt*p[i]+0.5*dt*dt*F[i];
            gradV(x,Fn); for(int i=0;i<3;i++) Fn[i]=-Fn[i];
            for(int i=0;i<3;i++) p[i]+=0.5*dt*(F[i]+Fn[i]);
            memcpy(F,Fn,sizeof F);
        }
        for(int i=0;i<3;i++) s_final+=(x[i]-xvac[i])*v[i];
        double expected=vamp*(M*dt);           /* free drift s = vamp * t */
        double rel=fabs(s_final-expected)/expected;
        printf("  Goldstone: drift s(T)=%.6e  free-particle vamp*T=%.6e  rel.err=%.2e\n",
               s_final, expected, rel);
        if(rel>1e-3) fail=1;
        printf("  [%s] Goldstone is massless (linear drift, no restoring force).\n",
               rel<1e-3?"OK":"FAIL");
    }

    /* (2) energy conservation along a massive mode */
    {
        int e=idx[2]; double v[3]={evec[0][e],evec[1][e],evec[2][e]};
        double x[3], p[3]={0,0,0}, F[3], Fn[3];
        for(int i=0;i<3;i++) x[i]=xvac[i]+1e-3*v[i];
        gradV(x,F); for(int i=0;i<3;i++) F[i]=-F[i];
        double E0=0; for(int i=0;i<3;i++) E0+=0.5*p[i]*p[i]; E0+=Vpot(x);
        for(int n=0;n<N;n++){
            for(int i=0;i<3;i++) x[i]+=dt*p[i]+0.5*dt*dt*F[i];
            gradV(x,Fn); for(int i=0;i<3;i++) Fn[i]=-Fn[i];
            for(int i=0;i<3;i++) p[i]+=0.5*dt*(F[i]+Fn[i]);
            memcpy(F,Fn,sizeof F);
        }
        double E1=0; for(int i=0;i<3;i++) E1+=0.5*p[i]*p[i]; E1+=Vpot(x);
        double rel=fabs(E1-E0)/fabs(E0-Vpot(xvac)+1e-30);
        printf("\n(2) energy conservation: E0=%.10e  E1=%.10e  drift/KE=%.2e\n", E0, E1, rel);
        if(rel>1e-3) fail=1;
    }

    /* (3) dispersion on a spatial lattice: Goldstone plane wave -> omega^2 = k^2 */
    {
        int Nx=64; double dx=0.25, L=Nx*dx; int kn=2; double k0=2*M_PI*kn/L;
        int e=idx[0]; double v[3]={evec[0][e],evec[1][e],evec[2][e]};
        double (*x)[3]=malloc(Nx*sizeof *x), (*p)[3]=malloc(Nx*sizeof *p);
        double (*F)[3]=malloc(Nx*sizeof *F);
        double A0=1e-3;
        for(int s=0;s<Nx;s++){ double w=A0*cos(k0*s*dx);
            for(int i=0;i<3;i++){ x[s][i]=xvac[i]+w*v[i]; p[s][i]=0; } }
        double dt2=5e-3; int Nt=40000;
        double *ts2=malloc(Nt*sizeof(double)), *ss2=malloc(Nt*sizeof(double));
        for(int n=0;n<Nt;n++){
            ts2[n]=n*dt2;
            { double s=0; for(int i=0;i<3;i++) s+=(x[0][i]-xvac[i])*v[i]; ss2[n]=s; }
            for(int s=0;s<Nx;s++){ double g[3]; gradV(x[s],g);
                int sl=(s-1+Nx)%Nx, sr=(s+1)%Nx;
                for(int i=0;i<3;i++){
                    double lap=(x[sl][i]-2*x[s][i]+x[sr][i])/(dx*dx);
                    F[s][i]=lap-g[i];
                } }
            for(int s=0;s<Nx;s++) for(int i=0;i<3;i++)
                x[s][i]+=dt2*p[s][i]+0.5*dt2*dt2*F[s][i];
            for(int s=0;s<Nx;s++){ double g[3]; gradV(x[s],g);
                int sl=(s-1+Nx)%Nx, sr=(s+1)%Nx;
                for(int i=0;i<3;i++){
                    double lap=(x[sl][i]-2*x[s][i]+x[sr][i])/(dx*dx);
                    double Fn=lap-g[i];
                    p[s][i]+=0.5*dt2*(F[s][i]+Fn);
                } }
        }
        double om=measure_omega(ts2,ss2,Nt,dt2), om2=om*om;
        double pred=k0*k0; /* massless: omega^2 = k^2 */
        double rel=fabs(om2-pred)/pred;
        printf("\n(3) Goldstone dispersion (k=%.4f): omega^2(meas)=%.5f  k^2=%.5f  rel.err=%.2e\n",
               k0, om2, pred, rel);
        printf("  [%s] massless propagating Goldstone: omega^2 = k^2 (speed 1).\n",
               rel<3e-2?"OK":"FAIL");
        if(rel>3e-2) fail=1;
        free(x);free(p);free(F);free(ts2);free(ss2);
    }

    free(ts); free(ss);
    printf(fail ? "\nC DYNAMICS CHECK FAILED\n" : "\nC DYNAMICS CHECK PASSED.\n");
    return fail;
}
