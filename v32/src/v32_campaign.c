/*  v32_campaign.c — Extended campaign: steady state, multi-braid, force law
 *
 *  The CORRECT interaction: the braid talks to the fabric ONLY where
 *  binding is weak. The tightly-bound core is self-contained.
 *
 *  acc_a = ∇²φ_a + w(P)×(∇ρ/ρ)·∇φ_a - m²φ_a - V'(φ_a)
 *  w(P) = 1/(1 + |P|/P_thresh)
 *
 *  At core: w≈0.22 → gradient coupling suppressed (self-bound)
 *  At surface (r≈4-6): w≈0.5 → half coupling (interaction zone)
 *  At fabric: w≈1.0 → full coupling (pure ambient)
 *
 *  Single field. No split. No smoothing. Periodic BC.
 *
 *  Build: gcc -O3 -march=native -fopenmp -o v32_wgrad src/v32_weighted_grad.c -lm
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>
#include <sys/stat.h>

#define NFIELDS 3
#define PI 3.14159265358979323846

static double MU = -41.345, KAPPA = 50.0, MASS2 = 2.25;
static double A_BG = 0.1;
static double GRAD_ALPHA = 0.5;  /* overall gradient strength (0.5 = sweet spot) */
static double P_THRESH = 0.0;    /* auto-set from initial data */
static int    PINNED_BC = 0;     /* 1 = pin borders to background, 0 = periodic */
static int    N_BRAIDS  = 1;
static double BRAID_X[10], BRAID_Y[10]; /* braid center positions (up to 10) */

typedef struct {
    double *phi[NFIELDS], *vel[NFIELDS], *acc[NFIELDS];
    double *rho;
    double *grad_rho[3];
    double *binding_w;  /* w(P) at each point */
    int N; double L, dx, dt;
} Grid;

static Grid *grid_alloc(int N, double L) {
    Grid *g = calloc(1, sizeof(Grid));
    long N3 = (long)N*N*N;
    for (int a = 0; a < NFIELDS; a++) {
        g->phi[a] = calloc(N3, sizeof(double));
        g->vel[a] = calloc(N3, sizeof(double));
        g->acc[a] = calloc(N3, sizeof(double));
    }
    g->rho = calloc(N3, sizeof(double));
    g->binding_w = calloc(N3, sizeof(double));
    for (int d = 0; d < 3; d++) g->grad_rho[d] = calloc(N3, sizeof(double));
    g->N = N; g->L = L;
    g->dx = 2.0*L/(N-1);
    g->dt = 0.12 * g->dx;
    return g;
}

static void grid_free(Grid *g) {
    for (int a = 0; a < NFIELDS; a++) {
        free(g->phi[a]); free(g->vel[a]); free(g->acc[a]);
    }
    free(g->rho); free(g->binding_w);
    for (int d=0;d<3;d++) free(g->grad_rho[d]);
    free(g);
}

/* Compute ρ, ∇ρ, and binding weight w(P). All periodic. */
static void compute_fields(Grid *g) {
    int N = g->N, NN = N*N;
    long N3 = (long)N*N*N;
    double dx = g->dx;

    /* Pass 1: ρ and w(P) */
    #pragma omp parallel for schedule(static)
    for (long idx = 0; idx < N3; idx++) {
        int i=(int)(idx/NN), j=(int)((idx/N)%N), k=(int)(idx%N);
        int ip=(i+1)%N, im=(i-1+N)%N, jp=(j+1)%N, jm=(j-1+N)%N;
        int kp=(k+1)%N, km=(k-1+N)%N;

        double e = 0;
        for (int a = 0; a < NFIELDS; a++) {
            e += 0.5 * g->vel[a][idx] * g->vel[a][idx];
            double gx=(g->phi[a][(long)ip*NN+j*N+k]-g->phi[a][(long)im*NN+j*N+k])/(2*dx);
            double gy=(g->phi[a][(long)i*NN+jp*N+k]-g->phi[a][(long)i*NN+jm*N+k])/(2*dx);
            double gz=(g->phi[a][(long)i*NN+j*N+kp]-g->phi[a][(long)i*NN+j*N+km])/(2*dx);
            e += 0.5*(gx*gx + gy*gy + gz*gz);
            e += 0.5 * MASS2 * g->phi[a][idx] * g->phi[a][idx];
        }
        double P = g->phi[0][idx]*g->phi[1][idx]*g->phi[2][idx];
        e += (MU/2.0)*P*P/(1.0+KAPPA*P*P);
        if (e < 1e-15) e = 1e-15;
        g->rho[idx] = e;

        /* Binding weight: w = 1/(1 + |P|/P_thresh) */
        double absP = fabs(P);
        g->binding_w[idx] = 1.0 / (1.0 + absP / (P_THRESH + 1e-15));
    }

    /* Pass 2: ∇ρ */
    #pragma omp parallel for schedule(static)
    for (long idx = 0; idx < N3; idx++) {
        int i=(int)(idx/NN), j=(int)((idx/N)%N), k=(int)(idx%N);
        int ip=(i+1)%N, im=(i-1+N)%N, jp=(j+1)%N, jm=(j-1+N)%N;
        int kp=(k+1)%N, km=(k-1+N)%N;
        g->grad_rho[0][idx]=(g->rho[(long)ip*NN+j*N+k]-g->rho[(long)im*NN+j*N+k])/(2*dx);
        g->grad_rho[1][idx]=(g->rho[(long)i*NN+jp*N+k]-g->rho[(long)i*NN+jm*N+k])/(2*dx);
        g->grad_rho[2][idx]=(g->rho[(long)i*NN+j*N+kp]-g->rho[(long)i*NN+j*N+km])/(2*dx);
    }
}

/* Forces with binding-weighted gradient coupling */
static void compute_forces(Grid *g) {
    int N = g->N, NN = N*N;
    long N3 = (long)N*N*N;
    double dx = g->dx, idx2 = 1.0/(dx*dx);

    #pragma omp parallel for schedule(static)
    for (long idx = 0; idx < N3; idx++) {
        int i=(int)(idx/NN), j=(int)((idx/N)%N), k=(int)(idx%N);
        int ip=(i+1)%N, im=(i-1+N)%N, jp=(j+1)%N, jm=(j-1+N)%N;
        int kp=(k+1)%N, km=(k-1+N)%N;

        /* Binding weight at this point */
        double w = g->binding_w[idx];

        /* Weighted gradient coupling strength */
        double rho_inv = GRAD_ALPHA * w / (g->rho[idx] + 1e-15);
        double grx = g->grad_rho[0][idx];
        double gry = g->grad_rho[1][idx];
        double grz = g->grad_rho[2][idx];

        double p0=g->phi[0][idx], p1=g->phi[1][idx], p2=g->phi[2][idx];
        double P = p0*p1*p2;
        double den = 1.0+KAPPA*P*P;
        double mPd2 = MU*P/(den*den);

        for (int a = 0; a < NFIELDS; a++) {
            /* Standard Laplacian */
            double lap = (g->phi[a][(long)ip*NN+j*N+k] + g->phi[a][(long)im*NN+j*N+k]
                        + g->phi[a][(long)i*NN+jp*N+k] + g->phi[a][(long)i*NN+jm*N+k]
                        + g->phi[a][(long)i*NN+j*N+kp] + g->phi[a][(long)i*NN+j*N+km]
                        - 6.0*g->phi[a][idx]) * idx2;

            /* Gradient of φ_a */
            double dpx=(g->phi[a][(long)ip*NN+j*N+k]-g->phi[a][(long)im*NN+j*N+k])/(2*dx);
            double dpy=(g->phi[a][(long)i*NN+jp*N+k]-g->phi[a][(long)i*NN+jm*N+k])/(2*dx);
            double dpz=(g->phi[a][(long)i*NN+j*N+kp]-g->phi[a][(long)i*NN+j*N+km])/(2*dx);

            /* BINDING-WEIGHTED gradient coupling: w(P) × (∇ρ/ρ) · ∇φ_a */
            double grad_coupling = rho_inv * (grx*dpx + gry*dpy + grz*dpz);

            /* Triple product force */
            double dPda = (a==0)?p1*p2:(a==1)?p0*p2:p0*p1;

            g->acc[a][idx] = lap + grad_coupling - MASS2*g->phi[a][idx] - mPd2*dPda;
        }
    }
}

static void verlet_step(Grid *g) {
    long N3 = (long)g->N*g->N*g->N;
    double hdt = 0.5*g->dt, dt = g->dt;
    for (int a=0;a<NFIELDS;a++) for (long idx=0;idx<N3;idx++)
        g->vel[a][idx] += hdt * g->acc[a][idx];
    for (int a=0;a<NFIELDS;a++) for (long idx=0;idx<N3;idx++)
        g->phi[a][idx] += dt * g->vel[a][idx];
    compute_fields(g);
    compute_forces(g);
    for (int a=0;a<NFIELDS;a++) for (long idx=0;idx<N3;idx++)
        g->vel[a][idx] += hdt * g->acc[a][idx];
}

static void init_braid(Grid *g, double x_cen) {
    int N=g->N, NN=N*N;
    double dx=g->dx, L=g->L;
    double A[3]={0.8,0.8,0.8}, delta[3]={0,3.0005,4.4325};
    double R=3.0, el=0.3325, kf=1.0;
    double kw=kf*PI/L, om=sqrt(kw*kw+MASS2);
    double sx=1+el, sy=1-el, inv2R2=1.0/(2*R*R);
    double k_bg=PI/L, om_bg=sqrt(k_bg*k_bg+MASS2);
    for(int i=0;i<N;i++){double x=-L+i*dx;
    for(int j=0;j<N;j++){double y=-L+j*dx;
    for(int kk=0;kk<N;kk++){double z=-L+kk*dx;
        long idx=(long)i*NN+j*N+kk;
        double xc=x-x_cen, xr=xc, yr=y;
        double r2e=xr*xr/(sx*sx)+yr*yr/(sy*sy);
        double env=exp(-r2e*inv2R2);
        for(int a=0;a<NFIELDS;a++){
            double ph=kw*z+delta[a];
            double ph_bg=k_bg*z+2*PI*a/3.0;
            g->phi[a][idx]+=A[a]*env*cos(ph)+A_BG*cos(ph_bg);
            g->vel[a][idx]+=om*A[a]*env*sin(ph)+om_bg*A_BG*sin(ph_bg);
        }}}}
}

static double measure_sep(Grid *g) {
    long N3=(long)g->N*g->N*g->N;
    int NN=g->N*g->N, N=g->N;
    double dx=g->dx, L=g->L;
    double avg_rho = 0;
    for (long idx=0;idx<N3;idx++) avg_rho += g->rho[idx];
    avg_rho /= N3;
    double thresh = 3.0 * avg_rho;
    double wL=0,wR=0,xL=0,xR=0;
    for(long idx=0;idx<N3;idx++){
        if (g->rho[idx] < thresh) continue;
        int i=(int)(idx/NN); double x=-L+i*dx;
        double w = g->rho[idx];
        if(x<0){wL+=w;xL+=x*w;} else {wR+=w;xR+=x*w;}
    }
    if(wL>0)xL/=wL; if(wR>0)xR/=wR;
    return xR-xL;
}

static void save_field(Grid *g, double t, const char *dir) {
    long N3=(long)g->N*g->N*g->N;
    char fn[512]; snprintf(fn,sizeof(fn),"%s/field_t%04d.bin",dir,(int)(t+0.5));
    FILE *fp=fopen(fn,"wb"); if(!fp)return;
    int n=g->N; double l=g->L;
    fwrite(&n,sizeof(int),1,fp); fwrite(&l,sizeof(double),1,fp); fwrite(&t,sizeof(double),1,fp);
    for(int a=0;a<NFIELDS;a++) fwrite(g->phi[a],sizeof(double),N3,fp);
    fwrite(g->rho,sizeof(double),N3,fp);
    fwrite(g->binding_w,sizeof(double),N3,fp);
    fclose(fp);
}

int main(int argc, char **argv) {
    int N=128; double L=30, T=300;
    int n_braids=2; double D=20;
    char outdir[256]="data/wgrad";

    int single_alpha = 0;  /* 1 = run only the specified alpha, skip sweep */
    for(int i=1;i<argc;i++){
        if(!strcmp(argv[i],"-N"))N=atoi(argv[++i]);
        else if(!strcmp(argv[i],"-L"))L=atof(argv[++i]);
        else if(!strcmp(argv[i],"-T"))T=atof(argv[++i]);
        else if(!strcmp(argv[i],"-ga")){GRAD_ALPHA=atof(argv[++i]); single_alpha=1;}
        else if(!strcmp(argv[i],"-braids")){N_BRAIDS=atoi(argv[++i]);n_braids=N_BRAIDS;}
        else if(!strcmp(argv[i],"-D"))D=atof(argv[++i]);
        else if(!strcmp(argv[i],"-o"))strncpy(outdir,argv[++i],255);
        else if(!strcmp(argv[i],"-pinned"))PINNED_BC=1;
        else if(!strcmp(argv[i],"-bg"))A_BG=atof(argv[++i]);
        else if(!strcmp(argv[i],"-pt"))P_THRESH=atof(argv[++i]);
        else if(!strcmp(argv[i],"-snap")){/* handled below */}
        else if(!strcmp(argv[i],"-diag")){/* handled below */}
    }
    n_braids = N_BRAIDS;

    omp_set_num_threads(16);
    printf("=== V32 Binding-Weighted Gradient Coupling ===\n");
    printf("acc = ∇²φ + w(P)×(∇ρ/ρ)·∇φ - m²φ - V'(φ)\n");
    printf("w(P) = 1/(1+|P|/P_thresh): core≈0.22, surface≈0.5, fabric≈1.0\n");
    printf("N=%d L=%.0f T=%.0f alpha=%.2f braids=%d D=%.0f\n\n",
           N,L,T,GRAD_ALPHA,n_braids,D);

    mkdir("data",0755); mkdir(outdir,0755);
    Grid *g = grid_alloc(N, L);
    long N3=(long)N*N*N;
    printf("dx=%.4f dt=%.5f\n",g->dx,g->dt);

    for(int a=0;a<NFIELDS;a++){memset(g->phi[a],0,N3*sizeof(double));memset(g->vel[a],0,N3*sizeof(double));}

    /* Initialize braids at specified positions */
    if (N_BRAIDS == 1) {
        BRAID_X[0] = 0; BRAID_Y[0] = 0;
    } else if (N_BRAIDS == 2 && BRAID_X[0] == 0 && BRAID_X[1] == 0) {
        BRAID_X[0] = -D/2; BRAID_X[1] = D/2;
        BRAID_Y[0] = 0; BRAID_Y[1] = 0;
    }
    /* For 5 braids, positions should be set via -bx/-by args or defaults */
    if (N_BRAIDS == 5 && BRAID_X[0] == 0 && BRAID_X[1] == 0) {
        /* Pentagon arrangement */
        for (int b = 0; b < 5; b++) {
            double angle = 2*PI*b/5;
            BRAID_X[b] = D/2 * cos(angle);
            BRAID_Y[b] = D/2 * sin(angle);
        }
    }
    for (int b = 0; b < N_BRAIDS; b++) {
        printf("  Braid %d at (%.1f, %.1f)\n", b, BRAID_X[b], BRAID_Y[b]);
        init_braid(g, BRAID_X[b]);
    }

    /* Auto-set P_THRESH from initial |P| at braid core */
    {
        double P_max = 0;
        for (long idx=0;idx<N3;idx++){
            double P=fabs(g->phi[0][idx]*g->phi[1][idx]*g->phi[2][idx]);
            if(P>P_max)P_max=P;
        }
        P_THRESH = P_max * 0.1;  /* 10% of peak */
        printf("P_max=%.6f, P_THRESH=%.6f (10%% of peak)\n", P_max, P_THRESH);
    }

    compute_fields(g);
    compute_forces(g);

    /* Scan: alpha=0 (control), 0.5, 1.0, 2.0  OR single alpha from -ga */
    double alphas_default[] = {0.0, 0.5, 1.0, 2.0};
    double alphas_single[1];
    double *alphas;
    int n_alphas;
    if (single_alpha) {
        alphas_single[0] = GRAD_ALPHA;
        alphas = alphas_single;
        n_alphas = 1;
    } else {
        alphas = alphas_default;
        n_alphas = 4;
    }

    char sumpath[512]; snprintf(sumpath,sizeof(sumpath),"%s/summary.tsv",outdir);
    FILE *fp_sum = fopen(sumpath,"w");
    fprintf(fp_sum,"alpha\tD_init\tD_final\tdeltaD\tE_init\tE_final\tstable\n");

    for (int ai = 0; ai < n_alphas; ai++) {
        GRAD_ALPHA = alphas[ai];
        char odir[256];
        if (single_alpha) {
            strncpy(odir, outdir, 255);
        } else {
            snprintf(odir,sizeof(odir),"data/wgrad_a%.1f",GRAD_ALPHA);
        }
        mkdir(odir,0755);

        printf("\n--- grad_alpha = %.1f ---\n", GRAD_ALPHA);

        /* Re-initialize */
        for(int a=0;a<NFIELDS;a++){memset(g->phi[a],0,N3*sizeof(double));memset(g->vel[a],0,N3*sizeof(double));}
        for (int b = 0; b < N_BRAIDS; b++) init_braid(g, BRAID_X[b]);
        compute_fields(g);
        compute_forces(g);

        char tsp[512]; snprintf(tsp,sizeof(tsp),"%s/timeseries.tsv",odir);
        FILE *fp=fopen(tsp,"w");
        fprintf(fp,"t\tE\tfc\tmax_rho\tD\tmax_grad_rho\tavg_w_core\n");

        int nsteps=(int)(T/g->dt);
        int diag_every=nsteps/60; if(diag_every<1)diag_every=1;
        int snap_every=nsteps/5; if(snap_every<1)snap_every=1;

        double E_init=0, D_init=0;
        int stable = 1;

        for(int step=0;step<=nsteps;step++){
            if(step>0) verlet_step(g);
            double t=step*g->dt;
            if(step%diag_every==0){
                double E=0,p2c=0,p2t=0,mr=0,mgr=0;
                double w_core_sum=0; int w_core_n=0;
                int NN=N*N;
                for(long idx=0;idx<N3;idx++){
                    int i=(int)(idx/NN),j=(int)((idx/N)%N);
                    double x=-L+i*g->dx, y=-L+j*g->dx;
                    double p2=0;
                    for(int a=0;a<NFIELDS;a++){
                        p2+=g->phi[a][idx]*g->phi[a][idx];
                        E+=0.5*g->vel[a][idx]*g->vel[a][idx]*g->dx*g->dx*g->dx;
                    }
                    p2t+=p2; if(x*x+y*y<64)p2c+=p2;
                    if(g->rho[idx]>mr)mr=g->rho[idx];
                    double gr=sqrt(g->grad_rho[0][idx]*g->grad_rho[0][idx]
                                  +g->grad_rho[1][idx]*g->grad_rho[1][idx]
                                  +g->grad_rho[2][idx]*g->grad_rho[2][idx]);
                    if(gr>mgr)mgr=gr;
                    if(x*x+y*y<25){w_core_sum+=g->binding_w[idx];w_core_n++;}
                }
                double fc=(p2t>0)?p2c/p2t:0;
                double Dsep=(n_braids>1)?measure_sep(g):0;
                double avg_w=(w_core_n>0)?w_core_sum/w_core_n:0;

                if(step==0){E_init=E;D_init=Dsep;}

                fprintf(fp,"%.1f\t%.2e\t%.4f\t%.4e\t%.2f\t%.4e\t%.4f\n",
                        t,E,fc,mr,Dsep,mgr,avg_w);
                fflush(fp);

                if(step%(diag_every*5)==0)
                    printf("t=%7.1f E=%.2e fc=%.3f D=%.1f |∇ρ|=%.2e w_core=%.3f\n",
                           t,E,fc,Dsep,mgr,avg_w);

                if(E!=E||mr>1e6){stable=0;printf("  BLOWUP!\n");break;}
            }
            if(step>0 && step%snap_every==0) save_field(g,t,odir);
        }
        fclose(fp);

        double D_final = measure_sep(g);
        double E_final = 0;
        for(long idx=0;idx<N3;idx++) for(int a=0;a<NFIELDS;a++)
            E_final+=0.5*g->vel[a][idx]*g->vel[a][idx]*g->dx*g->dx*g->dx;

        printf("  RESULT: D: %.2f → %.2f (ΔD=%+.2f) E: %.0f → %.0f %s\n",
               D_init, D_final, D_final-D_init, E_init, E_final,
               stable?"STABLE":"BLOWUP");

        fprintf(fp_sum,"%.1f\t%.2f\t%.2f\t%+.2f\t%.0f\t%.0f\t%d\n",
                GRAD_ALPHA, D_init, D_final, D_final-D_init, E_init, E_final, stable);
        fflush(fp_sum);
    }
    fclose(fp_sum);

    printf("\n=== SUMMARY ===\n");
    printf("alpha  D_init  D_final  deltaD   Verdict\n");
    /* Re-read summary */
    FILE *fr=fopen(sumpath,"r");
    char line[512]; fgets(line,sizeof(line),fr); /* skip header */
    while(fgets(line,sizeof(line),fr)){
        double a,di,df,dd,ei,ef; int s;
        sscanf(line,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%d",&a,&di,&df,&dd,&ei,&ef,&s);
        const char *v = !s?"BLOWUP":(dd<-1?"ATTRACTION":(dd>1?"REPULSION":"NEUTRAL"));
        printf("%.1f    %.1f    %.1f    %+.1f     %s\n",a,di,df,dd,v);
    }
    fclose(fr);

    grid_free(g);
    printf("\n=== Complete ===\n");
    return 0;
}
