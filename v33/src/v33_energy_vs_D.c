/*  v33_energy_vs_D.c — Measure total system energy as a function of braid separation D
 *
 *  For the energy minimization test: if F = -dE/dD, then gravity is
 *  the field relaxing to a lower energy configuration.
 *
 *  Method: Initialize two braids at separation D, let settle T=30,
 *  measure E_total averaged over t=20-30 (after transient dies down).
 *  Repeat for many D values. Compare -dE/dD to measured force from C1.
 *
 *  Build: gcc -O3 -march=native -fopenmp -o energy_vs_D src/v33_energy_vs_D.c -lm
 *  Usage: ./energy_vs_D [-single] [-D 15] [-N 128] [-L 40]
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>

#define NFIELDS 3
#define PI 3.14159265358979323846

static double MU = -41.345, KAPPA = 50.0, MASS2 = 2.25, A_BG = 0.1;

typedef struct {
    double *mem;
    double *phi[NFIELDS], *vel[NFIELDS], *acc[NFIELDS];
    int N; long N3;
    double L, dx, dt;
} Grid;

static Grid *grid_alloc(int N, double L) {
    Grid *g = calloc(1, sizeof(Grid));
    g->N = N; g->N3 = (long)N*N*N;
    g->L = L; g->dx = 2.0*L/(N-1);
    g->dt = 0.12 * g->dx;
    long total = 9 * g->N3;
    g->mem = malloc(total * sizeof(double));
    if (!g->mem) { fprintf(stderr, "FATAL: malloc\n"); exit(1); }
    memset(g->mem, 0, total * sizeof(double));
    for (int a = 0; a < NFIELDS; a++) {
        g->phi[a] = g->mem + (0+a)*g->N3;
        g->vel[a] = g->mem + (3+a)*g->N3;
        g->acc[a] = g->mem + (6+a)*g->N3;
    }
    return g;
}
static void grid_free(Grid *g) { free(g->mem); free(g); }

static void compute_forces(Grid *g) {
    const int N=g->N, NN=N*N; const long N3=g->N3;
    const double idx2 = 1.0/(g->dx*g->dx);
    #pragma omp parallel for schedule(static)
    for (long idx=0; idx<N3; idx++) {
        int i=(int)(idx/NN), j=(int)((idx/N)%N), k=(int)(idx%N);
        int ip=(i+1)%N,im=(i-1+N)%N,jp=(j+1)%N,jm=(j-1+N)%N,kp=(k+1)%N,km=(k-1+N)%N;
        double p0=g->phi[0][idx],p1=g->phi[1][idx],p2=g->phi[2][idx];
        double P=p0*p1*p2, den=1.0+KAPPA*P*P, mPd2=MU*P/(den*den);
        for (int a=0; a<NFIELDS; a++) {
            double lap = (g->phi[a][(long)ip*NN+j*N+k]+g->phi[a][(long)im*NN+j*N+k]
                         +g->phi[a][(long)i*NN+jp*N+k]+g->phi[a][(long)i*NN+jm*N+k]
                         +g->phi[a][(long)i*NN+j*N+kp]+g->phi[a][(long)i*NN+j*N+km]
                         -6.0*g->phi[a][idx])*idx2;
            double dPda=(a==0)?p1*p2:(a==1)?p0*p2:p0*p1;
            g->acc[a][idx] = lap - MASS2*g->phi[a][idx] - mPd2*dPda;
        }
    }
}

static void verlet_step(Grid *g) {
    const long N3=g->N3; const double hdt=0.5*g->dt, dt=g->dt;
    for (int a=0;a<NFIELDS;a++) { double *v=g->vel[a],*ac=g->acc[a]; for(long i=0;i<N3;i++) v[i]+=hdt*ac[i]; }
    for (int a=0;a<NFIELDS;a++) { double *p=g->phi[a],*v=g->vel[a]; for(long i=0;i<N3;i++) p[i]+=dt*v[i]; }
    compute_forces(g);
    for (int a=0;a<NFIELDS;a++) { double *v=g->vel[a],*ac=g->acc[a]; for(long i=0;i<N3;i++) v[i]+=hdt*ac[i]; }
}

static void init_braid(Grid *g, double x_cen, double y_cen) {
    const int N=g->N, NN=N*N; const double dx=g->dx, L=g->L;
    const double A[3]={0.8,0.8,0.8}, delta[3]={0,3.0005,4.4325};
    const double R_tube=3.0, ellip=0.3325, kw=PI/L;
    const double omega=sqrt(kw*kw+MASS2), sx=1+ellip, sy=1-ellip;
    const double inv2R2=1.0/(2*R_tube*R_tube);
    const double k_bg=PI/L, omega_bg=sqrt(k_bg*k_bg+MASS2);
    for (int i=0;i<N;i++) { double x=-L+i*dx;
    for (int j=0;j<N;j++) { double y=-L+j*dx;
    for (int kk=0;kk<N;kk++) { double z=-L+kk*dx;
        long idx=(long)i*NN+j*N+kk;
        double xc=x-x_cen, yc=y-y_cen;
        double r2e=xc*xc/(sx*sx)+yc*yc/(sy*sy);
        double env=exp(-r2e*inv2R2);
        for (int a=0;a<NFIELDS;a++) {
            double ph=kw*z+delta[a], ph_bg=k_bg*z+2*PI*a/3.0;
            g->phi[a][idx] += A[a]*env*cos(ph)+A_BG*cos(ph_bg);
            g->vel[a][idx] += omega*A[a]*env*sin(ph)+omega_bg*A_BG*sin(ph_bg);
        }
    }}}
}

static double compute_energy(Grid *g) {
    const int N=g->N, NN=N*N; const long N3=g->N3;
    const double dx=g->dx, dV=dx*dx*dx;
    double et=0;
    #pragma omp parallel for reduction(+:et) schedule(static)
    for (long idx=0; idx<N3; idx++) {
        int i=(int)(idx/NN),j=(int)((idx/N)%N),k=(int)(idx%N);
        int ip=(i+1)%N,im=(i-1+N)%N,jp=(j+1)%N,jm=(j-1+N)%N,kp=(k+1)%N,km=(k-1+N)%N;
        for (int a=0;a<NFIELDS;a++) {
            et += 0.5*g->vel[a][idx]*g->vel[a][idx]*dV;
            double gx=(g->phi[a][(long)ip*NN+j*N+k]-g->phi[a][(long)im*NN+j*N+k])/(2*dx);
            double gy=(g->phi[a][(long)i*NN+jp*N+k]-g->phi[a][(long)i*NN+jm*N+k])/(2*dx);
            double gz=(g->phi[a][(long)i*NN+j*N+kp]-g->phi[a][(long)i*NN+j*N+km])/(2*dx);
            et += 0.5*(gx*gx+gy*gy+gz*gz)*dV;
            et += 0.5*MASS2*g->phi[a][idx]*g->phi[a][idx]*dV;
        }
        double P=g->phi[0][idx]*g->phi[1][idx]*g->phi[2][idx];
        et += (MU/2.0)*P*P/(1.0+KAPPA*P*P)*dV;
    }
    return et;
}

int main(int argc, char **argv) {
    int N=128; double L=40.0, D=20.0, T_settle=30.0, T_avg_start=20.0;
    int single=0;

    for (int i=1;i<argc;i++) {
        if (!strcmp(argv[i],"-N")) N=atoi(argv[++i]);
        else if (!strcmp(argv[i],"-L")) L=atof(argv[++i]);
        else if (!strcmp(argv[i],"-D")) D=atof(argv[++i]);
        else if (!strcmp(argv[i],"-T")) T_settle=atof(argv[++i]);
        else if (!strcmp(argv[i],"-single")) single=1;
    }

    omp_set_num_threads(4);
    Grid *g = grid_alloc(N, L);

    if (single) {
        init_braid(g, 0, 0);
    } else {
        init_braid(g, -D/2, 0);
        init_braid(g, D/2, 0);
    }
    compute_forces(g);

    /* Evolve to T_settle, collecting energy averages in [T_avg_start, T_settle] */
    int n_steps = (int)(T_settle/g->dt);
    int avg_start_step = (int)(T_avg_start/g->dt);
    double E_sum = 0; int E_count = 0;
    int sample_every = (int)(1.0/g->dt);  /* every 1 time unit */
    if (sample_every < 1) sample_every = 1;

    for (int step = 0; step <= n_steps; step++) {
        if (step > 0) verlet_step(g);
        if (step >= avg_start_step && step % sample_every == 0) {
            E_sum += compute_energy(g);
            E_count++;
        }
    }

    double E_avg = E_sum / E_count;

    if (single)
        printf("SINGLE\t%.4f\t%d\t%.8e\t%d\n", L, N, E_avg, E_count);
    else
        printf("D=%.1f\t%.4f\t%d\t%.8e\t%d\n", D, L, N, E_avg, E_count);

    grid_free(g);
    return 0;
}
