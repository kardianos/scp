/*  t12_m7.c — Mechanism 7: Two-component field (S=braid, B=background)
 *  Coupling: acc_S += -g*B^2*S, acc_B += -g*S^2*B. g=0.01.
 *  Build: gcc -O3 -fopenmp -o t12_m7 src/t12_m7.c -lm
 */

#include "../../src/braid_core.h"

#define NBINS 60
#define A_BG  0.1
#define T_TOTAL 300.0
#define SNAP_DT 50.0

static double rho0_bg, mass2_phys;

/* Extra arrays for B fields */
static double *B_phi[NFIELDS], *B_vel[NFIELDS], *B_acc[NFIELDS];

static void compute_radial_rho(Grid *g,double *r_bins,double *rho_bins,int *counts,double dr){
    /* Measure TOTAL (S+B) energy density */
    int N=g->N,NN=N*N;double dx=g->dx,L=g->L;
    for(int b=0;b<NBINS;b++){r_bins[b]=(b+0.5)*dr;rho_bins[b]=0;counts[b]=0;}
    for(int i=1;i<N-1;i++){double x=-L+i*dx;for(int j=1;j<N-1;j++){double y=-L+j*dx;
        double rp=sqrt(x*x+y*y);int b=(int)(rp/dr);if(b>=NBINS)continue;
        for(int k=0;k<N;k++){int idx=i*NN+j*N+k;int kp=(k+1)%N,km=(k-1+N)%N;double e=0;
            /* S fields */
            for(int a=0;a<NFIELDS;a++){e+=0.5*g->vel[a][idx]*g->vel[a][idx];
                double gx=(g->phi[a][idx+NN]-g->phi[a][idx-NN])/(2*dx);
                double gy=(g->phi[a][idx+N]-g->phi[a][idx-N])/(2*dx);
                double gz=(g->phi[a][i*NN+j*N+kp]-g->phi[a][i*NN+j*N+km])/(2*dx);
                e+=0.5*(gx*gx+gy*gy+gz*gz);e+=0.5*mass2_phys*g->phi[a][idx]*g->phi[a][idx];}
            /* B fields */
            for(int a=0;a<NFIELDS;a++){e+=0.5*B_vel[a][idx]*B_vel[a][idx];
                double gx=(B_phi[a][idx+NN]-B_phi[a][idx-NN])/(2*dx);
                double gy=(B_phi[a][idx+N]-B_phi[a][idx-N])/(2*dx);
                double gz=(B_phi[a][i*NN+j*N+kp]-B_phi[a][i*NN+j*N+km])/(2*dx);
                e+=0.5*(gx*gx+gy*gy+gz*gz);e+=0.5*mass2_phys*B_phi[a][idx]*B_phi[a][idx];}
            rho_bins[b]+=e;counts[b]++;}}}
    for(int b=0;b<NBINS;b++)if(counts[b]>0)rho_bins[b]/=counts[b];}
static double get_rho_at_r(double *r,double *rho,int *c,double rt){
    double best=0,bd=1e30;for(int b=0;b<NBINS;b++){if(c[b]<5)continue;
    double d=fabs(r[b]-rt);if(d<bd){bd=d;best=rho[b];}}return best;}
static void apply_edge_damping(Grid *g){int N=g->N,NN=N*N;double dx=g->dx,L=g->L;
    double rs=0.85*L,re=0.98*L,idr=1.0/(re-rs+1e-30);
    for(int i=0;i<N;i++){double x=-L+i*dx;for(int j=0;j<N;j++){double y=-L+j*dx;
        double rp=sqrt(x*x+y*y);if(rp<=rs)continue;double f=(rp-rs)*idr;if(f>1)f=1;
        double d=1.0-0.95*f*f;for(int kk=0;kk<N;kk++){int idx=i*NN+j*N+kk;
            for(int a=0;a<NFIELDS;a++){g->vel[a][idx]*=d;B_vel[a][idx]*=d;}}}}}
static double compute_fc(Grid *g){int N=g->N,NN=N*N;double dx=g->dx,L=g->L;
    double Rc2=64.0,tot=0,cor=0;for(int i=1;i<N-1;i++){double x=-L+i*dx;
        for(int j=1;j<N-1;j++){double y=-L+j*dx;double rp2=x*x+y*y;
            for(int k=0;k<N;k++){int idx=i*NN+j*N+k;double r=0;
                for(int a=0;a<NFIELDS;a++)r+=g->phi[a][idx]*g->phi[a][idx];
                tot+=r;if(rp2<Rc2)cor+=r;}}}return cor/(tot+1e-30);}
static double compute_energy_m7(Grid *g,const double *phys,double gc){
    int N=g->N,NN=N*N;double dx=g->dx,dV=dx*dx*dx;
    double mu=phys[12],kappa=phys[13],lpw=phys[15];double E=0;
    for(int i=1;i<N-1;i++)for(int j=1;j<N-1;j++)for(int k=0;k<N;k++){
        int idx=i*NN+j*N+k;int kp=(k+1)%N,km=(k-1+N)%N;double ek=0,eg=0,em=0;
        for(int a=0;a<NFIELDS;a++){
            ek+=0.5*g->vel[a][idx]*g->vel[a][idx]+0.5*B_vel[a][idx]*B_vel[a][idx];
            double gx=(g->phi[a][idx+NN]-g->phi[a][idx-NN])/(2*dx);
            double gy=(g->phi[a][idx+N]-g->phi[a][idx-N])/(2*dx);
            double gz=(g->phi[a][i*NN+j*N+kp]-g->phi[a][i*NN+j*N+km])/(2*dx);
            eg+=0.5*(gx*gx+gy*gy+gz*gz);
            double bx=(B_phi[a][idx+NN]-B_phi[a][idx-NN])/(2*dx);
            double by=(B_phi[a][idx+N]-B_phi[a][idx-N])/(2*dx);
            double bz=(B_phi[a][i*NN+j*N+kp]-B_phi[a][i*NN+j*N+km])/(2*dx);
            eg+=0.5*(bx*bx+by*by+bz*bz);
            em+=0.5*mass2_phys*(g->phi[a][idx]*g->phi[a][idx]+B_phi[a][idx]*B_phi[a][idx]);}
        double p0=g->phi[0][idx],p1=g->phi[1][idx],p2=g->phi[2][idx],P=p0*p1*p2;
        double ep=(mu/2.0)*P*P/(1.0+kappa*P*P),epw=lpw*(p0*p1+p1*p2+p2*p0);
        /* Coupling energy */
        double S2=p0*p0+p1*p1+p2*p2;
        double B2=B_phi[0][idx]*B_phi[0][idx]+B_phi[1][idx]*B_phi[1][idx]+B_phi[2][idx]*B_phi[2][idx];
        double ecoup=gc*S2*B2;
        E+=(ek+eg+em+ep+epw+ecoup)*dV;}return E;}

/* ================================================================
   M7-specific: two-component forces
   S gets standard braid forces + coupling
   B gets free wave forces + coupling
   ================================================================ */

static void compute_forces_m7(Grid *g, const double *phys, double gc) {
    int N=g->N,NN=N*N,N3=N*N*N;
    double idx2=1.0/(g->dx*g->dx);
    double mu=phys[12],kappa=phys[13],lpw=phys[15];

    #pragma omp parallel for schedule(static)
    for(int idx=0;idx<N3;idx++){
        int i=idx/NN,j=(idx/N)%N,k=idx%N;
        if(i<1||i>=N-1||j<1||j>=N-1){
            for(int a=0;a<NFIELDS;a++){g->acc[a][idx]=0;B_acc[a][idx]=0;}continue;}
        int kp=(k+1)%N,km=(k-1+N)%N;
        int idx_kp=i*NN+j*N+kp,idx_km=i*NN+j*N+km;

        double p0=g->phi[0][idx],p1=g->phi[1][idx],p2=g->phi[2][idx];
        double P=p0*p1*p2;
        double denom=1.0+kappa*P*P;
        double mu_P_d2=mu*P/(denom*denom);

        double S2=p0*p0+p1*p1+p2*p2;
        double B2=B_phi[0][idx]*B_phi[0][idx]+B_phi[1][idx]*B_phi[1][idx]+B_phi[2][idx]*B_phi[2][idx];

        for(int a=0;a<NFIELDS;a++){
            /* S forces: standard + coupling */
            double lap_s=(g->phi[a][idx+NN]+g->phi[a][idx-NN]
                         +g->phi[a][idx+N]+g->phi[a][idx-N]
                         +g->phi[a][idx_kp]+g->phi[a][idx_km]
                         -6.0*g->phi[a][idx])*idx2;
            double dPda=(a==0)?p1*p2:(a==1)?p0*p2:p0*p1;
            double f_triple=mu_P_d2*dPda;
            double f_pw=lpw*(g->phi[(a+1)%3][idx]+g->phi[(a+2)%3][idx]);
            g->acc[a][idx]=lap_s-mass2_phys*g->phi[a][idx]-f_triple-f_pw
                           -gc*B2*g->phi[a][idx];

            /* B forces: free wave + coupling */
            double lap_b=(B_phi[a][idx+NN]+B_phi[a][idx-NN]
                         +B_phi[a][idx+N]+B_phi[a][idx-N]
                         +B_phi[a][idx_kp]+B_phi[a][idx_km]
                         -6.0*B_phi[a][idx])*idx2;
            B_acc[a][idx]=lap_b-mass2_phys*B_phi[a][idx]
                          -gc*S2*B_phi[a][idx];
        }
    }
}

int main(void) {
    bimodal_init_params();
    mass2_phys=2.25;
    double k_bg=PI/30.0,ob=sqrt(k_bg*k_bg+mass2_phys);
    rho0_bg=3.0*A_BG*A_BG*ob*ob;

    printf("=== T12 M7: Two-Component Field ===\n");
    printf("mass2=%.4f, rho0_bg=%.6f\n",mass2_phys,rho0_bg);

    int N=96; double L=30.0;
    Grid *g=grid_alloc(N,L);
    int N3=N*N*N;

    /* Allocate B fields */
    for(int a=0;a<NFIELDS;a++){
        B_phi[a]=calloc(N3,sizeof(double));
        B_vel[a]=calloc(N3,sizeof(double));
        B_acc[a]=calloc(N3,sizeof(double));
    }

    /* Initialize S = braid (no background) */
    init_braid(g,BIMODAL,-1);

    /* Initialize B = background traveling wave */
    int NN=N*N;double dx=g->dx;
    for(int i=0;i<N;i++)for(int j=0;j<N;j++)for(int k=0;k<N;k++){
        int idx=i*NN+j*N+k;double z=-L+k*dx;
        for(int a=0;a<NFIELDS;a++){
            double phase=PI/L*z+2.0*PI*a/3.0;
            B_phi[a][idx]=A_BG*cos(phase);
            B_vel[a][idx]=ob*A_BG*sin(phase);
        }
    }

    double gc=0.01;
    compute_forces_m7(g,BIMODAL,gc);

    double dr=L/NBINS;
    double r_bins[NBINS],rho_bins[NBINS];int counts[NBINS];
    double prev_drho15=0,prev_t=0;

    FILE *fout=fopen("data/t12_m7_timeseries.dat","w");
    fprintf(fout,"# t fc energy drho_r10 drho_r15 drho_r20\n");

    int n_total=(int)(T_TOTAL/g->dt),snap_every=(int)(SNAP_DT/g->dt);

    for(int step=0;step<=n_total;step++){
        if(step>0){
            double hdt=0.5*g->dt;
            /* Kick both S and B */
            for(int a=0;a<NFIELDS;a++){
                double *vs=g->vel[a],*as=g->acc[a],*vb=B_vel[a],*ab=B_acc[a];
                for(int idx=0;idx<N3;idx++){vs[idx]+=hdt*as[idx];vb[idx]+=hdt*ab[idx];}
            }
            /* Drift both */
            double dt=g->dt;
            for(int a=0;a<NFIELDS;a++){
                double *ps=g->phi[a],*vs=g->vel[a],*pb=B_phi[a],*vb=B_vel[a];
                for(int idx=0;idx<N3;idx++){ps[idx]+=dt*vs[idx];pb[idx]+=dt*vb[idx];}
            }
            /* Force */
            compute_forces_m7(g,BIMODAL,gc);
            /* Kick again */
            for(int a=0;a<NFIELDS;a++){
                double *vs=g->vel[a],*as=g->acc[a],*vb=B_vel[a],*ab=B_acc[a];
                for(int idx=0;idx<N3;idx++){vs[idx]+=hdt*as[idx];vb[idx]+=hdt*ab[idx];}
            }
            apply_edge_damping(g);
        }
        if(step%snap_every==0){
            double t=step*g->dt;
            compute_radial_rho(g,r_bins,rho_bins,counts,dr);
            double rho10=get_rho_at_r(r_bins,rho_bins,counts,10.0);
            double rho15=get_rho_at_r(r_bins,rho_bins,counts,15.0);
            double rho20=get_rho_at_r(r_bins,rho_bins,counts,20.0);
            double d10=rho10-rho0_bg,d15=rho15-rho0_bg,d20=rho20-rho0_bg;
            double fc=compute_fc(g),E=compute_energy_m7(g,BIMODAL,gc);
            double rate=(t>prev_t+1.0)?(d15-prev_drho15)/(t-prev_t):0;
            printf("t=%6.1f  fc=%.4f  E=%.1f  drho(10)=%+.4e  drho(15)=%+.4e  drho(20)=%+.4e  d/dt=%.2e\n",
                   t,fc,E,d10,d15,d20,rate);
            fprintf(fout,"%.2f %.6f %.2f %.8e %.8e %.8e\n",t,fc,E,d10,d15,d20);
            prev_drho15=d15;prev_t=t;
            if(check_blowup(g)){printf("BLOWUP at t=%.1f\n",t);break;}
        }
    }
    fclose(fout);
    for(int a=0;a<NFIELDS;a++){free(B_phi[a]);free(B_vel[a]);free(B_acc[a]);}
    grid_free(g);
    printf("=== M7 Complete ===\n");
    return 0;
}
