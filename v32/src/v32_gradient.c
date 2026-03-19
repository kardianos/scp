/*  v32_gradient.c — Single-field with gradient coupling term
 *
 *  The CORRECT wave equation on a non-uniform medium:
 *     ∂²φ/∂t² = (1/ρ)∇·(ρ∇φ) - m²φ - ∂V/∂φ
 *             = ∇²φ + (∇ρ/ρ)·∇φ - m²φ - ∂V/∂φ
 *
 *  The gradient coupling (∇ρ/ρ)·∇φ makes the braid drift along ∇ρ.
 *  This is the missing term that creates gravity.
 *
 *  Single field. No split. No smoothing. No c(ρ).
 *  Explicit symplectic Verlet. Periodic BC.
 *
 *  Build: gcc -O3 -march=native -fopenmp -o v32_grad src/v32_gradient.c -lm
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
static double GRAD_ALPHA = 1.0;  /* strength of gradient coupling (1.0 = physical) */

typedef struct {
    double *phi[NFIELDS], *vel[NFIELDS], *acc[NFIELDS];
    double *rho;      /* energy density */
    double *grad_rho[3]; /* ∇ρ components */
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
    for (int d = 0; d < 3; d++) g->grad_rho[d] = calloc(N3, sizeof(double));
    g->N = N; g->L = L;
    g->dx = 2.0*L/(N-1);
    g->dt = 0.12 * g->dx;  /* slightly conservative CFL */
    return g;
}

static void grid_free(Grid *g) {
    for (int a = 0; a < NFIELDS; a++) {
        free(g->phi[a]); free(g->vel[a]); free(g->acc[a]);
    }
    free(g->rho); for (int d=0;d<3;d++) free(g->grad_rho[d]);
    free(g);
}

/* Compute ρ(x) and ∇ρ(x). All periodic. */
static void compute_rho_and_grad(Grid *g) {
    int N = g->N, NN = N*N;
    long N3 = (long)N*N*N;
    double dx = g->dx;

    /* Step 1: compute ρ */
    #pragma omp parallel for schedule(static)
    for (long idx = 0; idx < N3; idx++) {
        int i=(int)(idx/NN), j=(int)((idx/N)%N), k=(int)(idx%N);
        int ip=(i+1)%N, im=(i-1+N)%N;
        int jp=(j+1)%N, jm=(j-1+N)%N;
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
    }

    /* Step 2: compute ∇ρ (central differences, periodic) */
    #pragma omp parallel for schedule(static)
    for (long idx = 0; idx < N3; idx++) {
        int i=(int)(idx/NN), j=(int)((idx/N)%N), k=(int)(idx%N);
        int ip=(i+1)%N, im=(i-1+N)%N;
        int jp=(j+1)%N, jm=(j-1+N)%N;
        int kp=(k+1)%N, km=(k-1+N)%N;

        g->grad_rho[0][idx] = (g->rho[(long)ip*NN+j*N+k] - g->rho[(long)im*NN+j*N+k]) / (2*dx);
        g->grad_rho[1][idx] = (g->rho[(long)i*NN+jp*N+k] - g->rho[(long)i*NN+jm*N+k]) / (2*dx);
        g->grad_rho[2][idx] = (g->rho[(long)i*NN+j*N+kp] - g->rho[(long)i*NN+j*N+km]) / (2*dx);
    }
}

/* Forces: ∇²φ + α(∇ρ/ρ)·∇φ - m²φ - V'(φ). All periodic. */
static void compute_forces(Grid *g) {
    int N = g->N, NN = N*N;
    long N3 = (long)N*N*N;
    double dx = g->dx, idx2 = 1.0/(dx*dx);

    #pragma omp parallel for schedule(static)
    for (long idx = 0; idx < N3; idx++) {
        int i=(int)(idx/NN), j=(int)((idx/N)%N), k=(int)(idx%N);
        int ip=(i+1)%N, im=(i-1+N)%N;
        int jp=(j+1)%N, jm=(j-1+N)%N;
        int kp=(k+1)%N, km=(k-1+N)%N;

        double rho_inv = GRAD_ALPHA / (g->rho[idx] + 1e-15);
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
            double dpx = (g->phi[a][(long)ip*NN+j*N+k] - g->phi[a][(long)im*NN+j*N+k]) / (2*dx);
            double dpy = (g->phi[a][(long)i*NN+jp*N+k] - g->phi[a][(long)i*NN+jm*N+k]) / (2*dx);
            double dpz = (g->phi[a][(long)i*NN+j*N+kp] - g->phi[a][(long)i*NN+j*N+km]) / (2*dx);

            /* GRADIENT COUPLING: (∇ρ/ρ) · ∇φ_a */
            double grad_coupling = rho_inv * (grx*dpx + gry*dpy + grz*dpz);

            /* Triple product */
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
    compute_rho_and_grad(g);
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
    /* Track braid peaks, not center of mass */
    long N3=(long)g->N*g->N*g->N;
    int NN=g->N*g->N, N=g->N;
    double dx=g->dx, L=g->L;

    /* Find high-energy particles in each half */
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

int main(int argc, char **argv) {
    int N=128; double L=30, T=300;
    int n_braids=2; double D=20;
    char outdir[256]="data/grad";

    for(int i=1;i<argc;i++){
        if(!strcmp(argv[i],"-N"))N=atoi(argv[++i]);
        else if(!strcmp(argv[i],"-L"))L=atof(argv[++i]);
        else if(!strcmp(argv[i],"-T"))T=atof(argv[++i]);
        else if(!strcmp(argv[i],"-ga"))GRAD_ALPHA=atof(argv[++i]);
        else if(!strcmp(argv[i],"-braids"))n_braids=atoi(argv[++i]);
        else if(!strcmp(argv[i],"-D"))D=atof(argv[++i]);
        else if(!strcmp(argv[i],"-o"))strncpy(outdir,argv[++i],255);
    }

    omp_set_num_threads(16);
    printf("=== V32 Gradient Coupling: (∇ρ/ρ)·∇φ ===\n");
    printf("N=%d L=%.0f T=%.0f grad_alpha=%.2f braids=%d D=%.0f\n\n",
           N,L,T,GRAD_ALPHA,n_braids,D);
    printf("Equation: ∂²φ/∂t² = ∇²φ + α(∇ρ/ρ)·∇φ - m²φ - V'(φ)\n");
    printf("Single field. No split. No smoothing. No c(ρ). Periodic BC.\n\n");

    mkdir("data",0755); mkdir(outdir,0755);
    Grid *g = grid_alloc(N, L);
    long N3=(long)N*N*N;
    printf("dx=%.4f dt=%.5f\n",g->dx,g->dt);

    for(int a=0;a<NFIELDS;a++){memset(g->phi[a],0,N3*sizeof(double));memset(g->vel[a],0,N3*sizeof(double));}
    if(n_braids==1) init_braid(g,0);
    else { init_braid(g,-D/2); init_braid(g,+D/2); }

    compute_rho_and_grad(g);
    compute_forces(g);

    char tsp[512]; snprintf(tsp,sizeof(tsp),"%s/timeseries.tsv",outdir);
    FILE *fp=fopen(tsp,"w");
    fprintf(fp,"t\tE\tfc\tmax_rho\tD\tmax_grad_rho\n");

    int nsteps=(int)(T/g->dt);
    int diag_every=nsteps/60; if(diag_every<1)diag_every=1;
    printf("Steps=%d\n\n",nsteps);

    double wall0=omp_get_wtime();

    for(int step=0;step<=nsteps;step++){
        if(step>0) verlet_step(g);
        double t=step*g->dt;
        if(step%diag_every==0){
            double E=0,p2c=0,p2t=0,mr=0,mgr=0;
            for(long idx=0;idx<N3;idx++){
                int i=(int)(idx/(N*N)),j=(int)((idx/N)%N);
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
            }
            double fc=(p2t>0)?p2c/p2t:0;
            double Dsep=(n_braids>1)?measure_sep(g):0;
            fprintf(fp,"%.1f\t%.2e\t%.4f\t%.4e\t%.2f\t%.4e\n",t,E,fc,mr,Dsep,mgr);
            fflush(fp);
            printf("t=%7.1f E=%.2e fc=%.3f rho_max=%.2e D=%.1f |∇ρ|_max=%.2e\n",
                   t,E,fc,mr,Dsep,mgr);
            fflush(stdout);
        }
    }
    fclose(fp);
    printf("\n=== Complete: %.0fs ===\n",omp_get_wtime()-wall0);
    grid_free(g);
    return 0;
}
