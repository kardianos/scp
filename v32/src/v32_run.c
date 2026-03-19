/*  v32_run.c — Single-run binding-weighted gradient coupling
 *
 *  For long E1/E2/E3 runs at high resolution (N=512).
 *  Same physics as v32_weighted_grad.c but:
 *  - No sweep loop (single run with CLI params)
 *  - More frequent snapshots for long runs
 *  - Radial profile output around each braid
 *  - Better diagnostics (potential + kinetic energy, peak tracking)
 *
 *  acc_a = nabla^2 phi_a + w(P)*alpha*(grad_rho/rho).grad_phi_a - m^2 phi_a - V'(phi_a)
 *  w(P) = 1/(1 + |P|/P_thresh)
 *
 *  Build: gcc -O3 -march=native -fopenmp -o v32_run src/v32_run.c -lm
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>
#include <sys/stat.h>
#include <time.h>

#define NFIELDS 3
#define PI 3.14159265358979323846

static double MU = -41.345, KAPPA = 50.0, MASS2 = 2.25;
static double A_BG = 0.1;
static double GRAD_ALPHA = 0.5;  /* best from sweep */
static double P_THRESH = 0.0;    /* auto-set from initial data */
static int    N_BRAIDS  = 1;

typedef struct {
    double *phi[NFIELDS], *vel[NFIELDS], *acc[NFIELDS];
    double *rho;
    double *grad_rho[3];
    double *binding_w;
    int N; double L, dx, dt;
} Grid;

static Grid *grid_alloc(int N, double L) {
    Grid *g = calloc(1, sizeof(Grid));
    if (!g) { fprintf(stderr, "Failed to alloc Grid struct\n"); exit(1); }
    long N3 = (long)N*N*N;
    printf("Allocating: %ld points x %d arrays = %.1f GB\n",
           N3, 3*NFIELDS+1+1+3, N3*(3.0*NFIELDS+1+1+3)*8.0/1e9);
    for (int a = 0; a < NFIELDS; a++) {
        g->phi[a] = calloc(N3, sizeof(double));
        g->vel[a] = calloc(N3, sizeof(double));
        g->acc[a] = calloc(N3, sizeof(double));
        if (!g->phi[a] || !g->vel[a] || !g->acc[a]) {
            fprintf(stderr, "ALLOC FAILED at field %d\n", a); exit(1);
        }
    }
    g->rho = calloc(N3, sizeof(double));
    g->binding_w = calloc(N3, sizeof(double));
    for (int d = 0; d < 3; d++) g->grad_rho[d] = calloc(N3, sizeof(double));
    if (!g->rho || !g->binding_w) {
        fprintf(stderr, "ALLOC FAILED for rho/w\n"); exit(1);
    }
    g->N = N; g->L = L;
    g->dx = 2.0*L/(N-1);
    g->dt = 0.25 * g->dx;  /* default; overridden by main after alloc */
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

/* Compute rho, grad_rho, and binding weight w(P). All periodic. */
static void compute_fields(Grid *g) {
    int N = g->N, NN = N*N;
    long N3 = (long)N*N*N;
    double dx = g->dx;

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

        double absP = fabs(P);
        g->binding_w[idx] = 1.0 / (1.0 + absP / (P_THRESH + 1e-15));
    }

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

static void compute_forces(Grid *g) {
    int N = g->N, NN = N*N;
    long N3 = (long)N*N*N;
    double dx = g->dx, idx2 = 1.0/(dx*dx);

    #pragma omp parallel for schedule(static)
    for (long idx = 0; idx < N3; idx++) {
        int i=(int)(idx/NN), j=(int)((idx/N)%N), k=(int)(idx%N);
        int ip=(i+1)%N, im=(i-1+N)%N, jp=(j+1)%N, jm=(j-1+N)%N;
        int kp=(k+1)%N, km=(k-1+N)%N;

        double w = g->binding_w[idx];
        double rho_inv = GRAD_ALPHA * w / (g->rho[idx] + 1e-15);
        double grx = g->grad_rho[0][idx];
        double gry = g->grad_rho[1][idx];
        double grz = g->grad_rho[2][idx];

        double p0=g->phi[0][idx], p1=g->phi[1][idx], p2=g->phi[2][idx];
        double P = p0*p1*p2;
        double den = 1.0+KAPPA*P*P;
        double mPd2 = MU*P/(den*den);

        for (int a = 0; a < NFIELDS; a++) {
            double lap = (g->phi[a][(long)ip*NN+j*N+k] + g->phi[a][(long)im*NN+j*N+k]
                        + g->phi[a][(long)i*NN+jp*N+k] + g->phi[a][(long)i*NN+jm*N+k]
                        + g->phi[a][(long)i*NN+j*N+kp] + g->phi[a][(long)i*NN+j*N+km]
                        - 6.0*g->phi[a][idx]) * idx2;

            double dpx=(g->phi[a][(long)ip*NN+j*N+k]-g->phi[a][(long)im*NN+j*N+k])/(2*dx);
            double dpy=(g->phi[a][(long)i*NN+jp*N+k]-g->phi[a][(long)i*NN+jm*N+k])/(2*dx);
            double dpz=(g->phi[a][(long)i*NN+j*N+kp]-g->phi[a][(long)i*NN+j*N+km])/(2*dx);

            double grad_coupling = rho_inv * (grx*dpx + gry*dpy + grz*dpz);
            double dPda = (a==0)?p1*p2:(a==1)?p0*p2:p0*p1;

            g->acc[a][idx] = lap + grad_coupling - MASS2*g->phi[a][idx] - mPd2*dPda;
        }
    }
}

static void verlet_step(Grid *g) {
    long N3 = (long)g->N*g->N*g->N;
    double hdt = 0.5*g->dt, dt = g->dt;
    for (int a=0;a<NFIELDS;a++) {
        #pragma omp parallel for schedule(static)
        for (long idx=0;idx<N3;idx++)
            g->vel[a][idx] += hdt * g->acc[a][idx];
    }
    for (int a=0;a<NFIELDS;a++) {
        #pragma omp parallel for schedule(static)
        for (long idx=0;idx<N3;idx++)
            g->phi[a][idx] += dt * g->vel[a][idx];
    }
    compute_fields(g);
    compute_forces(g);
    for (int a=0;a<NFIELDS;a++) {
        #pragma omp parallel for schedule(static)
        for (long idx=0;idx<N3;idx++)
            g->vel[a][idx] += hdt * g->acc[a][idx];
    }
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

/* Find braid center by energy-weighted centroid */
static void find_braid_centers(Grid *g, int n_braids, double *cx, double *cy) {
    long N3=(long)g->N*g->N*g->N;
    int NN=g->N*g->N, N=g->N;
    double dx=g->dx, L=g->L;
    double avg_rho = 0;
    for (long idx=0;idx<N3;idx++) avg_rho += g->rho[idx];
    avg_rho /= N3;
    double thresh = 3.0 * avg_rho;

    if (n_braids == 1) {
        double wx=0, wy=0, wt=0;
        for(long idx=0;idx<N3;idx++){
            if (g->rho[idx] < thresh) continue;
            int i=(int)(idx/NN), j=(int)((idx/N)%N);
            double x=-L+i*dx, y=-L+j*dx;
            double w = g->rho[idx];
            wx+=x*w; wy+=y*w; wt+=w;
        }
        if(wt>0){cx[0]=wx/wt; cy[0]=wy/wt;}
    } else {
        /* Two braids: split by x<0 and x>0 */
        double wL=0,wR=0,xL=0,yL=0,xR=0,yR=0;
        for(long idx=0;idx<N3;idx++){
            if (g->rho[idx] < thresh) continue;
            int i=(int)(idx/NN), j=(int)((idx/N)%N);
            double x=-L+i*dx, y=-L+j*dx;
            double w = g->rho[idx];
            if(x<0){wL+=w;xL+=x*w;yL+=y*w;}
            else   {wR+=w;xR+=x*w;yR+=y*w;}
        }
        if(wL>0){cx[0]=xL/wL;cy[0]=yL/wL;}
        if(wR>0){cx[1]=xR/wR;cy[1]=yR/wR;}
    }
}

static double measure_sep(Grid *g) {
    double cx[2]={0}, cy[2]={0};
    find_braid_centers(g, 2, cx, cy);
    double ddx=cx[1]-cx[0], ddy=cy[1]-cy[0];
    return sqrt(ddx*ddx+ddy*ddy);
}

static int COMPACT_SNAP = 0;  /* 1 = save only phi (60% smaller) */

static void save_field(Grid *g, double t, const char *dir) {
    long N3=(long)g->N*g->N*g->N;
    char fn[512]; snprintf(fn,sizeof(fn),"%s/field_t%04d.bin",dir,(int)(t+0.5));
    FILE *fp=fopen(fn,"wb"); if(!fp){fprintf(stderr,"Cannot write %s\n",fn);return;}
    int n=g->N; double l=g->L;
    fwrite(&n,sizeof(int),1,fp); fwrite(&l,sizeof(double),1,fp); fwrite(&t,sizeof(double),1,fp);
    for(int a=0;a<NFIELDS;a++) fwrite(g->phi[a],sizeof(double),N3,fp);
    if (!COMPACT_SNAP) {
        fwrite(g->rho,sizeof(double),N3,fp);
        fwrite(g->binding_w,sizeof(double),N3,fp);
    }
    fclose(fp);
    double nfields_saved = COMPACT_SNAP ? 3.0 : 5.0;
    printf("  [SNAP] Saved %s (%.1f GB)%s\n", fn, N3*nfields_saved*8.0/1e9,
           COMPACT_SNAP?" [compact]":"");
}

/* Save radial profile around a braid center */
static void save_radial_profile(Grid *g, double cx, double cy, double t, const char *dir, int braid_id) {
    int N=g->N, NN=N*N;
    double dx=g->dx, L=g->L;
    int nbins = 100;
    double rmax = 15.0;
    double dr = rmax / nbins;
    double *rho_sum = calloc(nbins, sizeof(double));
    double *w_sum   = calloc(nbins, sizeof(double));
    double *phi2_sum = calloc(nbins, sizeof(double));
    int *count = calloc(nbins, sizeof(int));

    /* Average over z, bin by (x-cx, y-cy) distance */
    for (int i=0;i<N;i++){double x=-L+i*dx;
    for (int j=0;j<N;j++){double y=-L+j*dx;
        double rx=x-cx, ry=y-cy;
        double r=sqrt(rx*rx+ry*ry);
        int bin=(int)(r/dr);
        if(bin>=nbins) continue;
        for(int kk=0;kk<N;kk++){
            long idx=(long)i*NN+j*N+kk;
            rho_sum[bin] += g->rho[idx];
            w_sum[bin]   += g->binding_w[idx];
            double p2=0;
            for(int a=0;a<NFIELDS;a++) p2+=g->phi[a][idx]*g->phi[a][idx];
            phi2_sum[bin] += p2;
            count[bin]++;
        }
    }}

    char fn[512]; snprintf(fn,sizeof(fn),"%s/profile_b%d_t%04d.tsv",dir,braid_id,(int)(t+0.5));
    FILE *fp=fopen(fn,"w");
    fprintf(fp,"r\trho\tw\tphi2\n");
    for(int b=0;b<nbins;b++){
        if(count[b]==0) continue;
        double n=(double)count[b];
        fprintf(fp,"%.3f\t%.6e\t%.6f\t%.6e\n",
                (b+0.5)*dr, rho_sum[b]/n, w_sum[b]/n, phi2_sum[b]/n);
    }
    fclose(fp);
    free(rho_sum); free(w_sum); free(phi2_sum); free(count);
}

/* Compute total energy (kinetic + gradient + mass + potential) */
static void compute_total_energy(Grid *g, double *E_kin, double *E_grad, double *E_mass, double *E_pot) {
    long N3=(long)g->N*g->N*g->N;
    int N=g->N, NN=N*N;
    double dx=g->dx, dV=dx*dx*dx;
    double ek=0, eg=0, em=0, ep=0;

    #pragma omp parallel for schedule(static) reduction(+:ek,eg,em,ep)
    for(long idx=0;idx<N3;idx++){
        int i=(int)(idx/NN), j=(int)((idx/N)%N), k=(int)(idx%N);
        int ip=(i+1)%N, im=(i-1+N)%N, jp=(j+1)%N, jm=(j-1+N)%N;
        int kp=(k+1)%N, km=(k-1+N)%N;

        for(int a=0;a<NFIELDS;a++){
            ek += 0.5*g->vel[a][idx]*g->vel[a][idx]*dV;
            double gx=(g->phi[a][(long)ip*NN+j*N+k]-g->phi[a][(long)im*NN+j*N+k])/(2*dx);
            double gy=(g->phi[a][(long)i*NN+jp*N+k]-g->phi[a][(long)i*NN+jm*N+k])/(2*dx);
            double gz=(g->phi[a][(long)i*NN+j*N+kp]-g->phi[a][(long)i*NN+j*N+km])/(2*dx);
            eg += 0.5*(gx*gx+gy*gy+gz*gz)*dV;
            em += 0.5*MASS2*g->phi[a][idx]*g->phi[a][idx]*dV;
        }
        double P=g->phi[0][idx]*g->phi[1][idx]*g->phi[2][idx];
        ep += (MU/2.0)*P*P/(1.0+KAPPA*P*P)*dV;
    }
    *E_kin=ek; *E_grad=eg; *E_mass=em; *E_pot=ep;
}

int main(int argc, char **argv) {
    int N=512; double L=20, T=2000, D=20;
    double dt_factor = 0.25;  /* dt = dt_factor * dx; CFL limit 0.577 */
    char outdir[256]="data/E1";
    double snap_interval = 50.0;  /* save snapshot every 50 time units */
    double diag_interval = 5.0;   /* diagnostics every 5 time units */
    double profile_interval = 100.0; /* radial profile every 100 time units */
    int nthreads = 16;

    for(int i=1;i<argc;i++){
        if(!strcmp(argv[i],"-N"))N=atoi(argv[++i]);
        else if(!strcmp(argv[i],"-L"))L=atof(argv[++i]);
        else if(!strcmp(argv[i],"-T"))T=atof(argv[++i]);
        else if(!strcmp(argv[i],"-ga"))GRAD_ALPHA=atof(argv[++i]);
        else if(!strcmp(argv[i],"-braids"))N_BRAIDS=atoi(argv[++i]);
        else if(!strcmp(argv[i],"-D"))D=atof(argv[++i]);
        else if(!strcmp(argv[i],"-o"))strncpy(outdir,argv[++i],255);
        else if(!strcmp(argv[i],"-bg"))A_BG=atof(argv[++i]);
        else if(!strcmp(argv[i],"-pt"))P_THRESH=atof(argv[++i]);
        else if(!strcmp(argv[i],"-snap"))snap_interval=atof(argv[++i]);
        else if(!strcmp(argv[i],"-diag"))diag_interval=atof(argv[++i]);
        else if(!strcmp(argv[i],"-prof"))profile_interval=atof(argv[++i]);
        else if(!strcmp(argv[i],"-threads"))nthreads=atoi(argv[++i]);
        else if(!strcmp(argv[i],"-dtf"))dt_factor=atof(argv[++i]);
        else if(!strcmp(argv[i],"-compact"))COMPACT_SNAP=1;
    }

    omp_set_num_threads(nthreads);
    setvbuf(stdout, NULL, _IOLBF, 0);  /* line-buffered for log files */

    printf("=== V32 Single Run: Binding-Weighted Gradient Coupling ===\n");
    printf("acc = nabla^2 phi + w(P)*alpha*(grad_rho/rho).grad_phi - m^2 phi - V'(phi)\n");
    printf("w(P) = 1/(1+|P|/P_thresh)\n");
    printf("N=%d L=%.1f T=%.0f alpha=%.2f braids=%d D=%.1f threads=%d\n",
           N,L,T,GRAD_ALPHA,N_BRAIDS,D,nthreads);
    printf("snap=%.0f diag=%.1f prof=%.0f\n\n", snap_interval, diag_interval, profile_interval);

    mkdir("data",0755); mkdir(outdir,0755);
    Grid *g = grid_alloc(N, L);
    g->dt = dt_factor * g->dx;
    long N3=(long)N*N*N;
    printf("dx=%.4f dt=%.6f (%.2f*dx, CFL=%.3f*dx)  steps=%d\n",
           g->dx, g->dt, dt_factor, 1.0/sqrt(3.0), (int)(T/g->dt));

    /* Initialize */
    for(int a=0;a<NFIELDS;a++){memset(g->phi[a],0,N3*sizeof(double));memset(g->vel[a],0,N3*sizeof(double));}
    if(N_BRAIDS==1){
        init_braid(g,0);
    } else if(N_BRAIDS==2){
        init_braid(g,-D/2);
        init_braid(g,+D/2);
    } else if(N_BRAIDS==5){
        for(int b=0;b<5;b++){
            double angle=2*PI*b/5;
            init_braid(g, D/2*cos(angle));
            /* Note: init_braid only offsets in x. For 5-braid need y offset too.
               For now, place them along x axis with spacing D */
        }
    }

    /* Auto-set P_THRESH */
    {
        double P_max = 0;
        for (long idx=0;idx<N3;idx++){
            double P=fabs(g->phi[0][idx]*g->phi[1][idx]*g->phi[2][idx]);
            if(P>P_max)P_max=P;
        }
        if (P_THRESH <= 0) P_THRESH = P_max * 0.1;
        printf("P_max=%.6f, P_THRESH=%.6f (%.0f%% of peak)\n\n",
               P_max, P_THRESH, 100.0*P_THRESH/P_max);
    }

    compute_fields(g);
    compute_forces(g);

    /* Open timeseries file */
    char tsp[512]; snprintf(tsp,sizeof(tsp),"%s/timeseries.tsv",outdir);
    FILE *fp=fopen(tsp,"w");
    fprintf(fp,"t\tE_kin\tE_grad\tE_mass\tE_pot\tE_total\tfc\tmax_rho\tD\tavg_w_core\tcx0\tcy0\n");

    int nsteps=(int)(T/g->dt);
    int diag_every=(int)(diag_interval/g->dt); if(diag_every<1) diag_every=1;
    int snap_every=(int)(snap_interval/g->dt); if(snap_every<1) snap_every=1;
    int prof_every=(int)(profile_interval/g->dt); if(prof_every<1) prof_every=1;

    printf("Total steps: %d, diag every %d, snap every %d, prof every %d\n\n",
           nsteps, diag_every, snap_every, prof_every);

    time_t t_start = time(NULL);
    double E_init_total = 0;

    /* Save initial snapshot */
    save_field(g, 0, outdir);

    for(int step=0;step<=nsteps;step++){
        if(step>0) verlet_step(g);
        double t=step*g->dt;

        if(step%diag_every==0){
            double E_kin, E_grad, E_mass, E_pot;
            compute_total_energy(g, &E_kin, &E_grad, &E_mass, &E_pot);
            double E_total = E_kin + E_grad + E_mass + E_pot;

            if(step==0) E_init_total = E_total;

            /* Braid tracking */
            double cx[2]={0}, cy[2]={0};
            find_braid_centers(g, N_BRAIDS, cx, cy);
            double Dsep = (N_BRAIDS>1) ? measure_sep(g) : 0;

            /* Confinement fraction */
            double p2c=0, p2t=0;
            int NN=N*N;
            double max_rho=0;
            double w_core_sum=0; long w_core_n=0;
            #pragma omp parallel for schedule(static) reduction(+:p2c,p2t,w_core_sum,w_core_n) reduction(max:max_rho)
            for(long idx=0;idx<N3;idx++){
                int i=(int)(idx/NN), j=(int)((idx/N)%N);
                double x=-L+i*g->dx, y=-L+j*g->dx;
                double p2=0;
                for(int a=0;a<NFIELDS;a++) p2+=g->phi[a][idx]*g->phi[a][idx];
                p2t+=p2;
                /* For single braid: center at origin; for two: each within R=8 */
                if(N_BRAIDS==1){
                    if(x*x+y*y<64) p2c+=p2;
                    if(x*x+y*y<25){w_core_sum+=g->binding_w[idx];w_core_n++;}
                } else {
                    double dx0=x-cx[0], dy0=y-cy[0], dx1=x-cx[1], dy1=y-cy[1];
                    double r0=dx0*dx0+dy0*dy0, r1=dx1*dx1+dy1*dy1;
                    if(r0<64||r1<64) p2c+=p2;
                    if(r0<25||r1<25){w_core_sum+=g->binding_w[idx];w_core_n++;}
                }
                if(g->rho[idx]>max_rho) max_rho=g->rho[idx];
            }
            double fc=(p2t>0)?p2c/p2t:0;
            double avg_w=(w_core_n>0)?w_core_sum/(double)w_core_n:0;

            fprintf(fp,"%.1f\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%.4f\t%.4e\t%.2f\t%.4f\t%.2f\t%.2f\n",
                    t,E_kin,E_grad,E_mass,E_pot,E_total,fc,max_rho,Dsep,avg_w,cx[0],cy[0]);
            fflush(fp);

            /* Print progress */
            if(step%(diag_every*10)==0 || step==0){
                double elapsed = difftime(time(NULL), t_start);
                double frac = (double)step/nsteps;
                double eta = (frac>0.001) ? elapsed*(1-frac)/frac : 0;
                printf("t=%7.1f [%.1f%%] E=%.3e (%.2f×E0) fc=%.3f D=%.1f w=%.3f max_rho=%.2e [%.0fs elapsed, ETA %.0fs]\n",
                       t, 100*frac, E_total, (E_init_total!=0)?E_total/E_init_total:1.0,
                       fc, Dsep, avg_w, max_rho, elapsed, eta);
            }

            if(E_total!=E_total || max_rho>1e8){
                printf("  BLOWUP at t=%.1f! Saving final state and aborting.\n", t);
                save_field(g, t, outdir);
                break;
            }
        }

        /* Snapshots */
        if(step>0 && step%snap_every==0){
            save_field(g, t, outdir);
        }

        /* Radial profiles */
        if(step%prof_every==0){
            double cx[2]={0}, cy[2]={0};
            find_braid_centers(g, N_BRAIDS, cx, cy);
            for(int b=0;b<N_BRAIDS && b<2;b++){
                save_radial_profile(g, cx[b], cy[b], t, outdir, b);
            }
        }
    }
    fclose(fp);

    /* Final summary */
    double E_kin, E_grad, E_mass, E_pot;
    compute_total_energy(g, &E_kin, &E_grad, &E_mass, &E_pot);
    double E_total = E_kin + E_grad + E_mass + E_pot;
    double elapsed = difftime(time(NULL), t_start);

    printf("\n=== FINAL RESULT ===\n");
    printf("T=%.0f, E_total=%.4e (%.2f x initial)\n", T, E_total,
           (E_init_total!=0)?E_total/E_init_total:1.0);
    printf("E_kin=%.4e  E_grad=%.4e  E_mass=%.4e  E_pot=%.4e\n",
           E_kin, E_grad, E_mass, E_pot);
    if(N_BRAIDS>1){
        printf("D_final=%.2f\n", measure_sep(g));
    }
    printf("Wall time: %.0f seconds (%.1f hours)\n", elapsed, elapsed/3600);

    /* Save final snapshot */
    save_field(g, T, outdir);

    grid_free(g);
    printf("=== Complete ===\n");
    return 0;
}
