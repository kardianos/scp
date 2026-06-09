/* v65 radial_cascade.c — test the user's cascade: phi activity -> density (rho), density
 * back-reacts on phi as pressure. Does the self-regulating cascade stabilize the phi-core
 * (a dissipative soliton), as Maxima predicts (charge K/lambda^3 -> stable minimum)?
 *
 * phi_a:  ddphi = lap(phi_a) - m2 phi_a - F_a(V) - g_back*rho*phi_a   (rho repels phi: pressure)
 * rho:    drho/dt = g_conv*Act - g_leak*rho + D*lap(rho)             (1st order; cascade source)
 *         Act = sum_a (1/2 dphi_a^2 + 1/2 (phi_a')^2)  (phi activity = the theta/density source)
 * Self-regulation: phi active -> rho grows -> g_back*rho suppresses phi where rho is high
 * -> phi forms a shell around the rho core; the conserved-ish rho resists compression.
 *
 * Build: gcc -O3 -march=native -o radial_cascade radial_cascade.c -lm
 * Usage: ./radial_cascade -A 1.2 -sigma 1.3 -gconv 0.5 -gback 2.0 -D 0.5 -leak 0.0 -T 200 ...
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#define NF 3
static double DELTA[3] = {0.0, 3.0005, 4.4325};

int main(int argc,char**argv){
    int N=1000; double R=30.0,A=1.2,sigma=1.3,T=200.0;
    double kappa=50.0,mu=-41.345,m2=2.25;
    double gconv=0.5, gback=2.0, D=0.5, leak=0.0;
    double rho0=0.0, rho_sig=1.5;   /* seed an initial density core (test Maxima mechanism) */
    double bsign=1.0, gdrift=1.0;   /* bsign:+1 bind phi to rho (attractive); gdrift unused (TF model) */
    double c_rho=1.0, Qd0=100.0;    /* c_rho: density pressure; Qd0: conserved charge (the cascade product) */
    double cool=0.0;                /* global phi cooling (gradient-flow toward the static minimum) */
    double damp=0.05, damp_width=6.0, dt_factor=0.2, diag_dt=5.0;
    for(int i=1;i<argc-1;i+=2){const char*k=argv[i],*v=argv[i+1];
        if(!strcmp(k,"-N"))N=atoi(v);else if(!strcmp(k,"-R"))R=atof(v);
        else if(!strcmp(k,"-A"))A=atof(v);else if(!strcmp(k,"-sigma"))sigma=atof(v);
        else if(!strcmp(k,"-T"))T=atof(v);else if(!strcmp(k,"-kappa"))kappa=atof(v);
        else if(!strcmp(k,"-mu"))mu=atof(v);else if(!strcmp(k,"-m2"))m2=atof(v);
        else if(!strcmp(k,"-gconv"))gconv=atof(v);else if(!strcmp(k,"-gback"))gback=atof(v);
        else if(!strcmp(k,"-D"))D=atof(v);else if(!strcmp(k,"-leak"))leak=atof(v);
        else if(!strcmp(k,"-rho0"))rho0=atof(v);else if(!strcmp(k,"-rho_sig"))rho_sig=atof(v);
        else if(!strcmp(k,"-bsign"))bsign=atof(v);else if(!strcmp(k,"-gdrift"))gdrift=atof(v);
        else if(!strcmp(k,"-c_rho"))c_rho=atof(v);else if(!strcmp(k,"-Qd"))Qd0=atof(v);
        else if(!strcmp(k,"-cool"))cool=atof(v);
        else if(!strcmp(k,"-diag_dt"))diag_dt=atof(v);
        else if(!strcmp(k,"-delta"))sscanf(v,"%lf,%lf,%lf",&DELTA[0],&DELTA[1],&DELTA[2]);
    }
    double dr=R/(N-1),dt=dt_factor*dr; int nsteps=(int)(T/dt),de=(int)(diag_dt/dt); if(de<1)de=1;
    double Dmax=0.4*dr*dr/dt; if(D>Dmax){fprintf(stderr,"  clamping D %.3f->%.3f (CFL)\n",D,Dmax);D=Dmax;}
    double *phi[NF],*vel[NF],*acc[NF],*rho=calloc(N,sizeof(double)),*r=malloc(N*sizeof(double));
    for(int a=0;a<NF;a++){phi[a]=calloc(N,sizeof(double));vel[a]=calloc(N,sizeof(double));acc[a]=calloc(N,sizeof(double));}
    for(int i=0;i<N;i++)r[i]=i*dr;
    for(int i=0;i<N;i++){double env=A*exp(-r[i]*r[i]/(2*sigma*sigma));for(int a=0;a<NF;a++)phi[a][i]=env*cos(DELTA[a]);
        rho[i]=rho0*exp(-r[i]*r[i]/(2*rho_sig*rho_sig));}   /* optional seeded density core */
    fprintf(stderr,"radial_cascade: A=%.2f sig=%.2f gconv=%.3f gback=%.3f D=%.2f leak=%.3f\n",A,sigma,gconv,gback,D,leak);
    #define LAP(f,i) ((i==0)?6.0*(f[1]-f[0])/(dr*dr): (i==N-1)?(f[N-2]-f[N-1])/(dr*dr): \
        (f[i+1]-2*f[i]+f[i-1])/(dr*dr)+(2.0/r[i])*(f[i+1]-f[i-1])/(2*dr))
    #define ACC() do{ for(int i=0;i<N;i++){ \
        double p0=phi[0][i],p1=phi[1][i],p2=phi[2][i],P=p0*p1*p2,den=1.0+kappa*P*P,inv=1.0/(den*den); \
        double pf[NF]; pf[0]=mu*p0*p1*p1*p2*p2*inv; pf[1]=mu*p1*p0*p0*p2*p2*inv; pf[2]=mu*p2*p0*p0*p1*p1*inv; \
        double dmp=(r[i]>R-damp_width)?damp*((r[i]-(R-damp_width))/damp_width):0.0; \
        for(int a=0;a<NF;a++){ double lap=LAP(phi[a],i); \
            acc[a][i]=lap - m2*phi[a][i] - pf[a] + bsign*gback*rho[i]*phi[a][i] - (dmp+cool)*vel[a][i]; } } }while(0)
    ACC();
    double Pmax0=0; for(int i=0;i<N;i++){double P=fabs(phi[0][i]*phi[1][i]*phi[2][i]);if(P>Pmax0)Pmax0=P;}
    printf("# t  Pmax  Qd(int rho)  rho_max  Ecore\n");
    for(int s=0;s<=nsteps;s++){
        if(s%de==0){double Pmax=0,Qd=0,rmax=0,Ec=0;
            for(int i=0;i<N;i++){double P=fabs(phi[0][i]*phi[1][i]*phi[2][i]);if(P>Pmax)Pmax=P;
                Qd+=rho[i]*r[i]*r[i]*dr*4*M_PI; if(rho[i]>rmax)rmax=rho[i];
                if(r[i]<5){double e=0;for(int a=0;a<NF;a++)e+=0.5*vel[a][i]*vel[a][i]+0.5*m2*phi[a][i]*phi[a][i];Ec+=e*r[i]*r[i]*dr;}}
            printf("%.1f %.5f %.4f %.5f %.4f\n",s*dt,Pmax,Qd,rmax,Ec);}
        /* leapfrog phi */
        for(int a=0;a<NF;a++)for(int i=0;i<N;i++){vel[a][i]+=0.5*dt*acc[a][i];phi[a][i]+=dt*vel[a][i];}
        /* Thomas-Fermi density: rho = max(0,(gback*phi2 + mu)/c_rho), mu set so int rho = Qd
           (conserved charge). gback binds rho to phi^2; c_rho is the pressure (caps the peak
           -> compressing the bound pair raises int rho^2 ~ Qd^2/lambda^3 = the Maxima barrier). */
        { double phi2[N]; for(int i=0;i<N;i++){phi2[i]=phi[0][i]*phi[0][i]+phi[1][i]*phi[1][i]+phi[2][i]*phi[2][i];}
          double mlo=-gback*phi2[0]-1.0, mhi=1.0+gback*phi2[0];   /* bracket mu */
          for(int it=0;it<60;it++){ double mm=0.5*(mlo+mhi),Q=0;
            for(int i=0;i<N;i++){double rr=(gback*phi2[i]+mm)/c_rho; if(rr>0)Q+=rr*r[i]*r[i]*dr*4*M_PI;}
            if(Q>Qd0)mhi=mm; else mlo=mm; }
          double mm=0.5*(mlo+mhi);
          for(int i=0;i<N;i++){double rr=(gback*phi2[i]+mm)/c_rho; rho[i]=(rr>0)?rr:0.0;} }
        ACC();
        for(int a=0;a<NF;a++)for(int i=0;i<N;i++)vel[a][i]+=0.5*dt*acc[a][i];
    }
    double Pmaxf=0,Qdf=0; for(int i=0;i<N;i++){double P=fabs(phi[0][i]*phi[1][i]*phi[2][i]);if(P>Pmaxf)Pmaxf=P;Qdf+=rho[i]*r[i]*r[i]*dr*4*M_PI;}
    printf("# RESULT gconv=%.3f gback=%.3f D=%.2f : Pmax0=%.4f Pmax_final=%.4f ratio=%.3f Qd_final=%.3f %s\n",
        gconv,gback,D,Pmax0,Pmaxf,Pmaxf/Pmax0,Qdf,(Pmaxf>0.3*Pmax0)?"STABILIZED":"decayed");
    return 0;
}
