/*  v31_signfix.c — Test 4 candidate fixes for the wrong-sign gravity
 *
 *  The problem: M7+c(rho_B) creates B-accretion near the braid (c>1),
 *  causing repulsion between braids instead of attraction.
 *
 *  Fix 1: Large separation (D=60, in the depletion zone at r>16)
 *  Fix 2: Asymmetric coupling (g_absorb=0.02, g_radiate=0.002)
 *  Fix 3: Inverted c formula: c² = 1 + alpha_c*(1 - rho_B/rho0)
 *         → c LOWER where B enriched (core), c HIGHER where B depleted
 *  Fix 4: c depends on differential: c² = 1 - alpha_c*(rho_B - rho0)/rho0
 *         Same math as Fix 3 but conceptually: depletion raises c, accretion lowers c
 *
 *  For each fix: two braids, track D(t). Compare to control (alpha_c=0).
 *
 *  Build: gcc -O3 -fopenmp -o v31_signfix src/v31_signfix.c -lm
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>
#include <sys/stat.h>

#define NFIELDS 3
#define PI 3.14159265358979323846
#define NDIM 16

/* Bimodal params */
static double BIMODAL[NDIM];
static const double PA[NDIM] = {
    0.8,0.8,0.8, 0.00,1.67, 3.0,0.80,0.0, 1.0,0.0, 0.0,0.0, -29.7,50.0,1.50,0.0};
static const double PB[NDIM] = {
    0.8,0.8,0.8, 3.53,4.92, 3.0,0.25,0.0, 1.0,0.0, 0.0,0.0, -43.4,50.0,1.50,0.0};

static void bimodal_init(void) {
    for (int d = 0; d < NDIM; d++) BIMODAL[d] = 0.15*PA[d] + 0.85*PB[d];
}

/* Config for each fix */
typedef struct {
    const char *name;
    double D_init;       /* initial separation */
    double g_StoB;       /* S→B coupling (radiation) */
    double g_BtoS;       /* B→S coupling (absorption) */
    double alpha_c;
    int invert_c;        /* 0: standard, 1: inverted */
} FixConfig;

/* Grid */
typedef struct {
    double *S_phi[NFIELDS], *S_vel[NFIELDS], *S_acc[NFIELDS];
    double *B_phi[NFIELDS], *B_vel[NFIELDS], *B_acc[NFIELDS];
    double *rho_B, *c2_eff;
    int N; double L, dx, dt;
} Grid;

static Grid *grid_alloc(int N, double L) {
    Grid *g = calloc(1, sizeof(Grid));
    int N3 = N*N*N;
    for (int a = 0; a < NFIELDS; a++) {
        g->S_phi[a]=calloc(N3,sizeof(double)); g->S_vel[a]=calloc(N3,sizeof(double));
        g->S_acc[a]=calloc(N3,sizeof(double)); g->B_phi[a]=calloc(N3,sizeof(double));
        g->B_vel[a]=calloc(N3,sizeof(double)); g->B_acc[a]=calloc(N3,sizeof(double));
    }
    g->rho_B=calloc(N3,sizeof(double)); g->c2_eff=calloc(N3,sizeof(double));
    g->N=N; g->L=L; g->dx=2.0*L/(N-1); g->dt=0.15*g->dx;
    return g;
}

static void grid_free(Grid *g) {
    for(int a=0;a<NFIELDS;a++){
        free(g->S_phi[a]);free(g->S_vel[a]);free(g->S_acc[a]);
        free(g->B_phi[a]);free(g->B_vel[a]);free(g->B_acc[a]);}
    free(g->rho_B);free(g->c2_eff);free(g);
}

static double RHO0_BG = 0;

static void compute_rhoB_ceff(Grid *g, double alpha_c, int invert) {
    int N=g->N, NN=N*N, N3=N*N*N;
    double mass2=BIMODAL[14]*BIMODAL[14];

    #pragma omp parallel for schedule(static)
    for(int idx=0;idx<N3;idx++){
        double e=0;
        for(int a=0;a<NFIELDS;a++){
            e+=0.5*g->B_vel[a][idx]*g->B_vel[a][idx];
            e+=0.5*mass2*g->B_phi[a][idx]*g->B_phi[a][idx];
        }
        g->rho_B[idx]=e;
    }

    #pragma omp parallel for schedule(static)
    for(int idx=0;idx<N3;idx++){
        double ratio=g->rho_B[idx]/(RHO0_BG+1e-30);
        double c2;
        if(invert){
            /* Fix 3/4: c LOWER where B enriched (ratio>1), HIGHER where depleted (ratio<1) */
            c2 = 1.0 + alpha_c*(1.0 - ratio);
        } else {
            /* Standard: c LOWER where B depleted (ratio<1) */
            c2 = 1.0 - alpha_c*(1.0 - ratio);
        }
        if(c2<0.3) c2=0.3;
        if(c2>2.0) c2=2.0;
        g->c2_eff[idx]=c2;
    }
}

static void compute_forces(Grid *g, double g_BtoS, double g_StoB) {
    int N=g->N, NN=N*N, N3=N*N*N;
    double idx2=1.0/(g->dx*g->dx);
    double mu=BIMODAL[12], kappa=BIMODAL[13], mass2=BIMODAL[14]*BIMODAL[14];

    #pragma omp parallel for schedule(static)
    for(int idx=0;idx<N3;idx++){
        int i=idx/NN,j=(idx/N)%N,k=idx%N;
        if(i<1||i>=N-1||j<1||j>=N-1){
            for(int a=0;a<NFIELDS;a++){g->S_acc[a][idx]=0;g->B_acc[a][idx]=0;}
            continue;
        }
        int kp=(k+1)%N,km=(k-1+N)%N;
        int ikp=i*NN+j*N+kp, ikm=i*NN+j*N+km;

        double s0=g->S_phi[0][idx],s1=g->S_phi[1][idx],s2=g->S_phi[2][idx];
        double S2=s0*s0+s1*s1+s2*s2;
        double P=s0*s1*s2;
        double den=1.0+kappa*P*P;
        double mPd2=mu*P/(den*den);
        double B2=0;
        for(int a=0;a<NFIELDS;a++) B2+=g->B_phi[a][idx]*g->B_phi[a][idx];
        double c2=g->c2_eff[idx];

        for(int a=0;a<NFIELDS;a++){
            double lap_s=(g->S_phi[a][idx+NN]+g->S_phi[a][idx-NN]
                         +g->S_phi[a][idx+N]+g->S_phi[a][idx-N]
                         +g->S_phi[a][ikp]+g->S_phi[a][ikm]
                         -6.0*g->S_phi[a][idx])*idx2;
            double dPda=(a==0)?s1*s2:(a==1)?s0*s2:s0*s1;
            g->S_acc[a][idx]=c2*lap_s - mass2*g->S_phi[a][idx]
                             - mPd2*dPda - g_BtoS*B2*g->S_phi[a][idx];

            double lap_b=(g->B_phi[a][idx+NN]+g->B_phi[a][idx-NN]
                         +g->B_phi[a][idx+N]+g->B_phi[a][idx-N]
                         +g->B_phi[a][ikp]+g->B_phi[a][ikm]
                         -6.0*g->B_phi[a][idx])*idx2;
            g->B_acc[a][idx]=lap_b - mass2*g->B_phi[a][idx]
                             - g_StoB*S2*g->B_phi[a][idx];
        }
    }
}

static void verlet_step(Grid *g, FixConfig *cfg) {
    int N3=g->N*g->N*g->N;
    double hdt=0.5*g->dt, dt=g->dt;
    for(int a=0;a<NFIELDS;a++) for(int idx=0;idx<N3;idx++){
        g->S_vel[a][idx]+=hdt*g->S_acc[a][idx];
        g->B_vel[a][idx]+=hdt*g->B_acc[a][idx];
    }
    for(int a=0;a<NFIELDS;a++) for(int idx=0;idx<N3;idx++){
        g->S_phi[a][idx]+=dt*g->S_vel[a][idx];
        g->B_phi[a][idx]+=dt*g->B_vel[a][idx];
    }
    compute_rhoB_ceff(g, cfg->alpha_c, cfg->invert_c);
    compute_forces(g, cfg->g_BtoS, cfg->g_StoB);
    for(int a=0;a<NFIELDS;a++) for(int idx=0;idx<N3;idx++){
        g->S_vel[a][idx]+=hdt*g->S_acc[a][idx];
        g->B_vel[a][idx]+=hdt*g->B_acc[a][idx];
    }
}

static void apply_damping(Grid *g) {
    int N=g->N,NN=N*N;
    double dx=g->dx,L=g->L,rs=0.70*L,re=0.95*L,idr=1.0/(re-rs+1e-30);
    for(int i=0;i<N;i++){double x=-L+i*dx;
    for(int j=0;j<N;j++){double y=-L+j*dx;
        double rp=sqrt(x*x+y*y); if(rp<=rs)continue;
        double f=(rp-rs)*idr; if(f>1)f=1; double d=1.0-0.98*f*f;
        for(int kk=0;kk<N;kk++){int idx=i*NN+j*N+kk;
            for(int a=0;a<NFIELDS;a++){
                g->S_phi[a][idx]*=d;g->S_vel[a][idx]*=d;
                g->B_phi[a][idx]*=d;g->B_vel[a][idx]*=d;
            }}}}
}

static void init_two_braids(Grid *g, double D, double A_bg) {
    int N=g->N, NN=N*N;
    double dx=g->dx, L=g->L;
    double A[3]={BIMODAL[0],BIMODAL[1],BIMODAL[2]};
    double delta[3]={0,BIMODAL[3],BIMODAL[4]};
    double R_tube=BIMODAL[5],ellip=BIMODAL[6],ell_ang=BIMODAL[7];
    double k_fac=BIMODAL[8],mass=BIMODAL[14];
    double kw=k_fac*PI/L, omega=sqrt(kw*kw+mass*mass);
    double sx=1+ellip,sy=1-ellip,inv2R2=1.0/(2*R_tube*R_tube);
    double ca=cos(ell_ang),sa=sin(ell_ang);
    double k_bg=PI/L, omega_bg=sqrt(k_bg*k_bg+mass*mass);
    double x1=-D/2, x2=D/2;

    int N3=N*N*N;
    for(int a=0;a<NFIELDS;a++){
        memset(g->S_phi[a],0,N3*sizeof(double)); memset(g->S_vel[a],0,N3*sizeof(double));
        memset(g->S_acc[a],0,N3*sizeof(double)); memset(g->B_acc[a],0,N3*sizeof(double));
    }

    for(int i=0;i<N;i++){double x=-L+i*dx;
    for(int j=0;j<N;j++){double y=-L+j*dx;
    for(int kk=0;kk<N;kk++){double z=-L+kk*dx;
        int idx=i*NN+j*N+kk;
        /* Two S braids */
        for(int b=0;b<2;b++){
            double xc=(b==0)?x-x1:x-x2;
            double xr=xc*ca+y*sa, yr=-xc*sa+y*ca;
            double r2e=xr*xr/(sx*sx)+yr*yr/(sy*sy);
            double env=exp(-r2e*inv2R2);
            for(int a=0;a<NFIELDS;a++){
                double ph=kw*z+delta[a];
                g->S_phi[a][idx]+=A[a]*env*cos(ph);
                g->S_vel[a][idx]+=omega*A[a]*env*sin(ph);
            }
        }
        /* B background */
        for(int a=0;a<NFIELDS;a++){
            double ph=k_bg*z+2*PI*a/3.0;
            g->B_phi[a][idx]=A_bg*cos(ph);
            g->B_vel[a][idx]=omega_bg*A_bg*sin(ph);
        }
    }}}
}

static double measure_separation(Grid *g) {
    int N=g->N, NN=N*N, N3=N*N*N;
    double dx=g->dx, L=g->L;
    double wL=0,wR=0,xL=0,xR=0;
    for(int idx=0;idx<N3;idx++){
        int i=idx/NN; double x=-L+i*dx;
        double p2=0;
        for(int a=0;a<NFIELDS;a++) p2+=g->S_phi[a][idx]*g->S_phi[a][idx];
        if(x<0){wL+=p2; xL+=x*p2;}
        else   {wR+=p2; xR+=x*p2;}
    }
    if(wL>0)xL/=wL; if(wR>0)xR/=wR;
    return xR-xL;
}

/* ================================================================ */

int main(int argc, char **argv) {
    bimodal_init();
    omp_set_num_threads(16);

    int N=128; double L=40.0, T=400.0, A_bg=0.1;

    FixConfig fixes[] = {
        /* Control */
        {"control",  20.0, 0.01, 0.01, 0.0,  0},
        /* Fix 1: Large separation (D=60) */
        {"fix1_D60", 60.0, 0.01, 0.01, 0.2,  0},
        /* Fix 2: Asymmetric coupling (absorb 10× stronger than radiate) */
        {"fix2_asym",20.0, 0.002,0.02, 0.2,  0},
        /* Fix 3: Inverted c (c LOWER where B enriched) */
        {"fix3_inv", 20.0, 0.01, 0.01, 0.2,  1},
        /* Fix 4: Strong inverted c */
        {"fix4_inv5",20.0, 0.01, 0.01, 0.5,  1},
    };
    int n_fixes = 5;

    /* Fix 1 needs larger domain */
    double L_fix1 = 80.0;

    mkdir("data",0755);

    FILE *fp = fopen("data/signfix_results.tsv","w");
    fprintf(fp,"fix\tt\tD\tfc\tmin_c2\tavg_c2_mid\n");

    printf("=== V31 Sign-Fix Test: 5 configurations ===\n\n");

    for(int fi=0;fi<n_fixes;fi++){
        FixConfig *cfg = &fixes[fi];
        double this_L = (fi==1) ? L_fix1 : L;  /* larger L for fix1 */
        int this_N = (fi==1) ? 192 : N;        /* more points for larger L */

        printf("--- %s: D=%.0f, g_BtoS=%.3f, g_StoB=%.3f, ac=%.2f, inv=%d, N=%d, L=%.0f ---\n",
               cfg->name, cfg->D_init, cfg->g_BtoS, cfg->g_StoB,
               cfg->alpha_c, cfg->invert_c, this_N, this_L);

        Grid *g = grid_alloc(this_N, this_L);
        init_two_braids(g, cfg->D_init, A_bg);

        /* Set rho0 from far field */
        compute_rhoB_ceff(g, 0, 0);
        double sum=0; int cnt=0;
        int NN=this_N*this_N;
        for(int i=0;i<this_N;i++){double x=-this_L+i*g->dx;
        for(int j=0;j<this_N;j++){double y=-this_L+j*g->dx;
            if(sqrt(x*x+y*y)<this_L*0.5)continue;
            for(int kk=0;kk<this_N;kk++){
                sum+=g->rho_B[i*NN+j*this_N+kk]; cnt++;}}}
        RHO0_BG=sum/cnt;
        printf("  rho0_BG=%.4e\n",RHO0_BG);

        compute_rhoB_ceff(g, cfg->alpha_c, cfg->invert_c);
        compute_forces(g, cfg->g_BtoS, cfg->g_StoB);

        int n_steps=(int)(T/g->dt);
        int diag_every=n_steps/40; if(diag_every<1)diag_every=1;

        for(int step=0;step<=n_steps;step++){
            if(step>0){
                verlet_step(g, cfg);
                apply_damping(g);
            }
            if(step%diag_every==0){
                double t=step*g->dt;
                double D=measure_separation(g);
                /* fc */
                double p2c=0,p2t=0;
                int N3=this_N*this_N*this_N;
                for(int idx=0;idx<N3;idx++){
                    int i=idx/NN; double x=-this_L+i*g->dx;
                    int j_=(idx/this_N)%this_N; double y=-this_L+j_*g->dx;
                    double p2=0;
                    for(int a=0;a<NFIELDS;a++)p2+=g->S_phi[a][idx]*g->S_phi[a][idx];
                    p2t+=p2; if(x*x+y*y<100)p2c+=p2;
                }
                double fc=(p2t>0)?p2c/p2t:0;
                /* min c2 and avg c2 at midpoint */
                double min_c2=2,c2_mid_sum=0; int c2_mid_n=0;
                for(int idx=0;idx<N3;idx++){
                    if(g->c2_eff[idx]<min_c2)min_c2=g->c2_eff[idx];
                    int i=idx/NN; double x=-this_L+i*g->dx;
                    if(fabs(x)<5){c2_mid_sum+=g->c2_eff[idx];c2_mid_n++;}
                }
                double avg_c2_mid=(c2_mid_n>0)?c2_mid_sum/c2_mid_n:1;

                fprintf(fp,"%s\t%.1f\t%.4f\t%.4f\t%.4f\t%.4f\n",
                        cfg->name,t,D,fc,min_c2,avg_c2_mid);
                fflush(fp);

                if(step%(diag_every*5)==0)
                    printf("  t=%6.0f D=%7.2f fc=%.3f min_c2=%.3f c2_mid=%.3f\n",
                           t,D,fc,min_c2,avg_c2_mid);

                /* Blowup check */
                if(D!=D||fc!=fc){printf("  BLOWUP!\n");break;}
            }
        }
        printf("  Final: D=%.2f (started %.2f, delta=%+.2f)\n\n",
               measure_separation(g), cfg->D_init,
               measure_separation(g)-cfg->D_init);
        grid_free(g);
    }

    fclose(fp);

    /* Summary */
    printf("=== SUMMARY ===\n");
    printf("Fix            D_init  D_final  deltaD   Verdict\n");
    /* Re-read the file for final D values */
    FILE *fr=fopen("data/signfix_results.tsv","r");
    char line[512]; double last_D[5]={0}; int ci=-1;
    char last_name[64]="";
    while(fgets(line,sizeof(line),fr)){
        char name[64]; double t,D;
        if(sscanf(line,"%s %lf %lf",name,&t,&D)>=3){
            if(strcmp(name,last_name)!=0){ci++;strncpy(last_name,name,63);}
            if(ci>=0&&ci<5) last_D[ci]=D;
        }
    }
    fclose(fr);
    for(int i=0;i<n_fixes;i++){
        double dD=last_D[i]-fixes[i].D_init;
        const char *v=(dD<-1)?"ATTRACTION":(dD>1)?"REPULSION":"NEUTRAL";
        printf("%-14s  %5.0f   %7.2f  %+6.2f   %s\n",
               fixes[i].name,fixes[i].D_init,last_D[i],dD,v);
    }

    printf("\n=== Complete ===\n");
    return 0;
}
