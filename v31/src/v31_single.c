/*  v31_single.c — Single-field explicit Verlet with rational c(rho)
 *
 *  ONE field phi_a. No S/B split. No smoothing. No implicit solver.
 *  Symplectic Velocity Verlet (energy-conserving).
 *  Fully periodic BC (background sustains itself).
 *
 *  c²(x) = 1 + alpha_c * (rho0 - rho(x)) / (rho0 + rho(x))
 *  Bounded: c² ∈ [1-alpha_c, 1+alpha_c]. No extremes.
 *  At rho=rho0: c²=1. Dense→slow. Dilute→fast.
 *
 *  Build: gcc -O3 -march=native -fopenmp -o v31_single src/v31_single.c -lm
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
static double ALPHA_C = 0.3;
static double A_BG = 0.1;
static double RHO0 = 0.0;

typedef struct {
    double *phi[NFIELDS], *vel[NFIELDS], *acc[NFIELDS];
    double *c2;
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
    g->c2 = calloc(N3, sizeof(double));
    g->N = N; g->L = L;
    g->dx = 2.0*L/(N-1);
    /* CFL: dt < dx / sqrt(c2_max). c2_max = 1+alpha_c. */
    g->dt = 0.15 * g->dx / sqrt(1.0 + ALPHA_C);
    return g;
}

static void grid_free(Grid *g) {
    for (int a = 0; a < NFIELDS; a++) {
        free(g->phi[a]); free(g->vel[a]); free(g->acc[a]);
    }
    free(g->c2); free(g);
}

/* Compute c²(x) from current phi, vel — NO smoothing */
static void compute_c2(Grid *g) {
    int N = g->N, NN = N*N;
    long N3 = (long)N*N*N;
    double dx = g->dx;

    #pragma omp parallel for schedule(static)
    for (long idx = 0; idx < N3; idx++) {
        int i = (int)(idx/NN), j = (int)((idx/N)%N), k = (int)(idx%N);
        int ip=(i+1)%N, im=(i-1+N)%N;
        int jp=(j+1)%N, jm=(j-1+N)%N;
        int kp=(k+1)%N, km=(k-1+N)%N;

        double e = 0;
        for (int a = 0; a < NFIELDS; a++) {
            e += 0.5 * g->vel[a][idx] * g->vel[a][idx];
            double gx = (g->phi[a][(long)ip*NN+j*N+k] - g->phi[a][(long)im*NN+j*N+k]) / (2*dx);
            double gy = (g->phi[a][(long)i*NN+jp*N+k] - g->phi[a][(long)i*NN+jm*N+k]) / (2*dx);
            double gz = (g->phi[a][(long)i*NN+j*N+kp] - g->phi[a][(long)i*NN+j*N+km]) / (2*dx);
            e += 0.5*(gx*gx + gy*gy + gz*gz);
            e += 0.5 * MASS2 * g->phi[a][idx] * g->phi[a][idx];
        }
        double P = g->phi[0][idx]*g->phi[1][idx]*g->phi[2][idx];
        e += (MU/2.0)*P*P/(1.0+KAPPA*P*P);
        if (e < 1e-15) e = 1e-15;

        /* Rational c²: bounded [1-alpha, 1+alpha] */
        g->c2[idx] = 1.0 + ALPHA_C * (RHO0 - e) / (RHO0 + e);
    }
}

/* Forces: c²(x) × Laplacian + mass + triple product. All periodic. */
static void compute_forces(Grid *g) {
    int N = g->N, NN = N*N;
    long N3 = (long)N*N*N;
    double idx2 = 1.0 / (g->dx * g->dx);

    #pragma omp parallel for schedule(static)
    for (long idx = 0; idx < N3; idx++) {
        int i=(int)(idx/NN), j=(int)((idx/N)%N), k=(int)(idx%N);
        int ip=(i+1)%N, im=(i-1+N)%N;
        int jp=(j+1)%N, jm=(j-1+N)%N;
        int kp=(k+1)%N, km=(k-1+N)%N;

        double c2 = g->c2[idx];
        double p0=g->phi[0][idx], p1=g->phi[1][idx], p2=g->phi[2][idx];
        double P = p0*p1*p2;
        double den = 1.0 + KAPPA*P*P;
        double mPd2 = MU*P/(den*den);

        for (int a = 0; a < NFIELDS; a++) {
            double lap = (g->phi[a][(long)ip*NN+j*N+k] + g->phi[a][(long)im*NN+j*N+k]
                        + g->phi[a][(long)i*NN+jp*N+k] + g->phi[a][(long)i*NN+jm*N+k]
                        + g->phi[a][(long)i*NN+j*N+kp] + g->phi[a][(long)i*NN+j*N+km]
                        - 6.0*g->phi[a][idx]) * idx2;

            double dPda = (a==0)?p1*p2:(a==1)?p0*p2:p0*p1;

            g->acc[a][idx] = c2 * lap - MASS2 * g->phi[a][idx] - mPd2 * dPda;
        }
    }
}

/* Symplectic Velocity Verlet — energy-conserving */
static void verlet_step(Grid *g) {
    long N3 = (long)g->N * g->N * g->N;
    double hdt = 0.5*g->dt, dt = g->dt;

    for (int a = 0; a < NFIELDS; a++)
        for (long idx = 0; idx < N3; idx++)
            g->vel[a][idx] += hdt * g->acc[a][idx];
    for (int a = 0; a < NFIELDS; a++)
        for (long idx = 0; idx < N3; idx++)
            g->phi[a][idx] += dt * g->vel[a][idx];

    compute_c2(g);
    compute_forces(g);

    for (int a = 0; a < NFIELDS; a++)
        for (long idx = 0; idx < N3; idx++)
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

    for (int i=0;i<N;i++){double x=-L+i*dx;
    for (int j=0;j<N;j++){double y=-L+j*dx;
    for (int kk=0;kk<N;kk++){double z=-L+kk*dx;
        long idx=(long)i*NN+j*N+kk;
        double xc=x-x_cen;
        double xr=xc, yr=y;
        double r2e=xr*xr/(sx*sx)+yr*yr/(sy*sy);
        double env=exp(-r2e*inv2R2);
        for (int a=0;a<NFIELDS;a++){
            double ph=kw*z+delta[a];
            double ph_bg=k_bg*z+2*PI*a/3.0;
            g->phi[a][idx] += A[a]*env*cos(ph) + A_BG*cos(ph_bg);
            g->vel[a][idx] += om*A[a]*env*sin(ph) + om_bg*A_BG*sin(ph_bg);
        }
    }}}
}

static double measure_sep(Grid *g) {
    long N3=(long)g->N*g->N*g->N;
    int NN=g->N*g->N;
    double dx=g->dx, L=g->L;
    double wL=0,wR=0,xL=0,xR=0;
    for(long idx=0;idx<N3;idx++){
        int i=(int)(idx/NN); double x=-L+i*dx;
        double p2=0;
        for(int a=0;a<NFIELDS;a++) p2+=g->phi[a][idx]*g->phi[a][idx];
        if(x<0){wL+=p2;xL+=x*p2;} else {wR+=p2;xR+=x*p2;}
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
    fwrite(g->c2,sizeof(double),N3,fp);
    fclose(fp);
    printf("  Saved: %s (%.0f MB)\n",fn,(4)*N3*8.0/1e6);
}

int main(int argc, char **argv) {
    int N=128; double L=20, T=500;
    int n_braids=1; double D=20;
    char outdir[256]="data/single";

    for(int i=1;i<argc;i++){
        if(!strcmp(argv[i],"-N"))N=atoi(argv[++i]);
        else if(!strcmp(argv[i],"-L"))L=atof(argv[++i]);
        else if(!strcmp(argv[i],"-T"))T=atof(argv[++i]);
        else if(!strcmp(argv[i],"-ac"))ALPHA_C=atof(argv[++i]);
        else if(!strcmp(argv[i],"-braids"))n_braids=atoi(argv[++i]);
        else if(!strcmp(argv[i],"-D"))D=atof(argv[++i]);
        else if(!strcmp(argv[i],"-bg"))A_BG=atof(argv[++i]);
        else if(!strcmp(argv[i],"-o"))strncpy(outdir,argv[++i],255);
    }

    omp_set_num_threads(16);
    printf("=== V31 Single-Field Explicit Verlet + Rational c(rho) ===\n");
    printf("N=%d L=%.0f T=%.0f alpha_c=%.2f braids=%d D=%.0f bg=%.2f\n",
           N,L,T,ALPHA_C,n_braids,D,A_BG);
    printf("No split. No smoothing. Symplectic. Periodic BC.\n\n");

    mkdir("data",0755); mkdir(outdir,0755);
    Grid *g = grid_alloc(N, L);
    long N3=(long)N*N*N;
    printf("dx=%.4f dt=%.5f c2_range=[%.3f,%.3f]\n",g->dx,g->dt,1-ALPHA_C,1+ALPHA_C);

    for(int a=0;a<NFIELDS;a++){memset(g->phi[a],0,N3*sizeof(double));memset(g->vel[a],0,N3*sizeof(double));}
    if(n_braids==1) init_braid(g,0);
    else { init_braid(g,-D/2); init_braid(g,+D/2); }

    /* Set RHO0 from far field */
    compute_c2(g); /* temp, just to get rho */
    { double s=0;int c=0;int NN=N*N;
      for(int i=0;i<N;i++){double x=-L+i*g->dx;
      for(int j=0;j<N;j++){double y=-L+j*g->dx;
        if(sqrt(x*x+y*y)<L*0.4)continue;
        for(int kk=0;kk<N;kk++){long idx=(long)i*NN+j*N+kk;
          double e=0;
          for(int a=0;a<NFIELDS;a++){e+=0.5*g->vel[a][idx]*g->vel[a][idx]+0.5*MASS2*g->phi[a][idx]*g->phi[a][idx];}
          s+=e;c++;}}}
      RHO0=s/c;
    }
    printf("rho0=%.6e\n",RHO0);
    compute_c2(g); compute_forces(g);

    char tsp[512]; snprintf(tsp,sizeof(tsp),"%s/timeseries.tsv",outdir);
    FILE *fp=fopen(tsp,"w");
    fprintf(fp,"t\tE\tfc\tmax_rho\tmin_c2\tmax_c2\tD\n");

    int nsteps=(int)(T/g->dt);
    int diag_every=nsteps/50; if(diag_every<1)diag_every=1;
    int snap_every=nsteps/5; if(snap_every<1)snap_every=1;
    printf("Steps=%d diag_every=%d snap_every=%d\n\n",nsteps,diag_every,snap_every);

    double wall0=omp_get_wtime();
    save_field(g,0,outdir);

    for(int step=0;step<=nsteps;step++){
        if(step>0) verlet_step(g);
        double t=step*g->dt;

        if(step%diag_every==0){
            double E=0,p2c=0,p2t=0,mr=0,mnc=2,mxc=0;
            #pragma omp parallel for reduction(+:E,p2c,p2t) reduction(max:mr,mxc) reduction(min:mnc)
            for(long idx=0;idx<N3;idx++){
                int i=(int)(idx/(N*N)),j=(int)((idx/N)%N);
                double x=-L+i*g->dx, y=-L+j*g->dx;
                double p2=0;
                for(int a=0;a<NFIELDS;a++){
                    p2+=g->phi[a][idx]*g->phi[a][idx];
                    E+=0.5*g->vel[a][idx]*g->vel[a][idx]*g->dx*g->dx*g->dx;
                }
                p2t+=p2; if(x*x+y*y<64)p2c+=p2;
                double e=0; for(int a=0;a<NFIELDS;a++) e+=g->phi[a][idx]*g->phi[a][idx];
                if(e>mr)mr=e;
                if(g->c2[idx]<mnc)mnc=g->c2[idx];
                if(g->c2[idx]>mxc)mxc=g->c2[idx];
            }
            double fc=(p2t>0)?p2c/p2t:0;
            double Dsep=(n_braids>1)?measure_sep(g):0;

            fprintf(fp,"%.1f\t%.2e\t%.4f\t%.4e\t%.4f\t%.4f\t%.2f\n",
                    t,E,fc,mr,mnc,mxc,Dsep);
            fflush(fp);
            printf("t=%7.1f E=%.2e fc=%.3f max_rho=%.2e c2=[%.3f,%.3f] D=%.1f\n",
                   t,E,fc,mr,mnc,mxc,Dsep);
            fflush(stdout);
        }
        if(step>0 && step%snap_every==0) save_field(g,t,outdir);
    }
    save_field(g,T,outdir);
    fclose(fp);

    printf("\n=== Complete: %.0fs (%.1f min) ===\n",omp_get_wtime()-wall0,(omp_get_wtime()-wall0)/60);
    printf("Avg step: %.4fs\n",(omp_get_wtime()-wall0)/nsteps);
    grid_free(g);
    return 0;
}
