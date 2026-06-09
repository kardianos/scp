/* v65 radial_oscillon.c — 1+1D spherically-symmetric solver for the 6-field Cosserat
 * phi-sector, to search for a coherent OSCILLON (the back-calculated stable structure).
 *
 * EOM (phi-only, eta=0; theta added later if a phi-oscillon is found):
 *   phi_a'' (time) = phi_a'' + (2/r)phi_a' - m^2 phi_a - mu*phi_a*PROD_{b!=a}phi_b^2/(1+kappa P^2)^2
 *   P = phi0*phi1*phi2.
 * Methods: leapfrog; r=0 regularity (lap0 = 6(phi[1]-phi[0])/dr^2); boundary sponge;
 * optional bulk cooling (to relax onto the oscillon attractor and shed radiation).
 *
 * Build: gcc -O3 -march=native -o radial_oscillon radial_oscillon.c -lm
 * Usage: ./radial_oscillon -N 1000 -R 30 -A 0.7 -sigma 1.3 -omega 1.3 -T 200 \
 *                          -kappa 50 -mu -41.345 -m2 2.25 -cool 0 -damp 0.05 -phase 1
 *   -phase 0: launch from rest, phi_a=A*env*cos(delta_a)  (the old seed)
 *   -phase 1: time-phase launch, vel_a = -A*env*omega*sin(delta_a)  (on the limit cycle)
 * Prints P_max(t) and core energy; final line = lifetime / coherence verdict.
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define NF 3
static double DELTA[3] = {0.0, 3.0005, 4.4325};

int main(int argc, char **argv){
    int N=1000; double R=30.0, A=0.7, sigma=1.3, omega=1.3, T=200.0;
    double kappa=50.0, mu=-41.345, m2=2.25;
    double cool=0.0, damp=0.05, damp_width=6.0;
    int phase=1; double dt_factor=0.2; double diag_dt=1.0;
    for(int i=1;i<argc-1;i+=2){ const char*k=argv[i],*v=argv[i+1];
        if(!strcmp(k,"-N"))N=atoi(v); else if(!strcmp(k,"-R"))R=atof(v);
        else if(!strcmp(k,"-A"))A=atof(v); else if(!strcmp(k,"-sigma"))sigma=atof(v);
        else if(!strcmp(k,"-omega"))omega=atof(v); else if(!strcmp(k,"-T"))T=atof(v);
        else if(!strcmp(k,"-kappa"))kappa=atof(v); else if(!strcmp(k,"-mu"))mu=atof(v);
        else if(!strcmp(k,"-m2"))m2=atof(v); else if(!strcmp(k,"-cool"))cool=atof(v);
        else if(!strcmp(k,"-damp"))damp=atof(v); else if(!strcmp(k,"-damp_width"))damp_width=atof(v);
        else if(!strcmp(k,"-phase"))phase=atoi(v); else if(!strcmp(k,"-diag_dt"))diag_dt=atof(v);
        else if(!strcmp(k,"-delta"))sscanf(v,"%lf,%lf,%lf",&DELTA[0],&DELTA[1],&DELTA[2]);
    }
    double dr=R/(N-1), dt=dt_factor*dr;
    int nsteps=(int)(T/dt), diag_every=(int)(diag_dt/dt); if(diag_every<1)diag_every=1;
    double *phi[NF],*vel[NF],*acc[NF],*r=malloc(N*sizeof(double));
    for(int a=0;a<NF;a++){phi[a]=calloc(N,sizeof(double));vel[a]=calloc(N,sizeof(double));acc[a]=calloc(N,sizeof(double));}
    for(int i=0;i<N;i++) r[i]=i*dr;
    /* init */
    for(int i=0;i<N;i++){ double env=A*exp(-r[i]*r[i]/(2*sigma*sigma));
        for(int a=0;a<NF;a++){ phi[a][i]=env*cos(DELTA[a]);
            vel[a][i]= (phase==1)? -env*omega*sin(DELTA[a]) : 0.0; }
    }
    fprintf(stderr,"radial_oscillon: N=%d R=%.1f dr=%.3f dt=%.4f A=%.3f sigma=%.2f omega=%.2f phase=%d cool=%.3f\n",
            N,R,dr,dt,A,sigma,omega,phase,cool);
    /* force */
    #define COMPUTE_ACC() do{ for(int i=0;i<N;i++){ \
        double p0=phi[0][i],p1=phi[1][i],p2=phi[2][i]; double P=p0*p1*p2; double den=1.0+kappa*P*P; double inv=1.0/(den*den); \
        double pf[NF]; pf[0]=mu*p0*(p1*p1)*(p2*p2)*inv; pf[1]=mu*p1*(p0*p0)*(p2*p2)*inv; pf[2]=mu*p2*(p0*p0)*(p1*p1)*inv; \
        for(int a=0;a<NF;a++){ double lap; \
            if(i==0) lap=6.0*(phi[a][1]-phi[a][0])/(dr*dr); \
            else if(i==N-1) lap=(phi[a][N-2]-phi[a][N-1])/(dr*dr); /* one-sided, sponge handles it */ \
            else lap=(phi[a][i+1]-2*phi[a][i]+phi[a][i-1])/(dr*dr) + (2.0/r[i])*(phi[a][i+1]-phi[a][i-1])/(2*dr); \
            double dmp = (r[i]>R-damp_width)? damp*( (r[i]-(R-damp_width))/damp_width ) : 0.0; \
            acc[a][i]=lap - m2*phi[a][i] - pf[a] - (dmp+cool)*vel[a][i]; \
        } } }while(0)
    COMPUTE_ACC();
    double Pmax0=0; for(int i=0;i<N;i++){double P=fabs(phi[0][i]*phi[1][i]*phi[2][i]); if(P>Pmax0)Pmax0=P;}
    printf("# t  P_max  E_core(r<5)\n");
    double life=0, Pmax_peak=Pmax0;
    for(int s=0;s<=nsteps;s++){
        if(s%diag_every==0){
            double Pmax=0,Ecore=0;
            for(int i=0;i<N;i++){ double P=fabs(phi[0][i]*phi[1][i]*phi[2][i]); if(P>Pmax)Pmax=P;
                if(r[i]<5.0){ double e=0; for(int a=0;a<NF;a++){ e+=0.5*vel[a][i]*vel[a][i]+0.5*m2*phi[a][i]*phi[a][i]; } Ecore+=e*r[i]*r[i]*dr; } }
            printf("%.2f %.6f %.6f\n", s*dt, Pmax, Ecore);
            if(Pmax>0.1*Pmax0) life=s*dt;       /* lifetime: P_max still >10% of initial peak */
        }
        /* leapfrog (velocity Verlet) */
        for(int a=0;a<NF;a++) for(int i=0;i<N;i++){ vel[a][i]+=0.5*dt*acc[a][i]; phi[a][i]+=dt*vel[a][i]; }
        COMPUTE_ACC();
        for(int a=0;a<NF;a++) for(int i=0;i<N;i++) vel[a][i]+=0.5*dt*acc[a][i];
    }
    /* coherence verdict */
    double Pmax_final=0; for(int i=0;i<N;i++){double P=fabs(phi[0][i]*phi[1][i]*phi[2][i]); if(P>Pmax_final)Pmax_final=P;}
    printf("# RESULT A=%.3f sigma=%.2f omega=%.2f phase=%d cool=%.3f : Pmax0=%.4f Pmax_final=%.4f ratio=%.3f lifetime=%.1f %s\n",
           A,sigma,omega,phase,cool,Pmax0,Pmax_final,Pmax_final/Pmax0,life,
           (Pmax_final>0.3*Pmax0)?"COHERENT":"decayed");
    return 0;
}
