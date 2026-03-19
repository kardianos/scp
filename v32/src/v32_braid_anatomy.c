/*  v32_braid_anatomy.c — Geometric analysis of the braid's binding structure
 *
 *  Analyze: WHERE is the braid bound? WHERE does it interact with fabric?
 *  Map |P(x)|, ∇ρ, energy flux across the binding isosurface.
 *
 *  Output: 2D slices through the braid showing:
 *  - |P(x)|: binding strength
 *  - ρ(x): energy density
 *  - ∇ρ(x): gradient vectors
 *  - Energy flux: v_a × ∇φ_a (direction of energy flow)
 *  - Binding weakness: w(P) = 1/(1 + |P|/P_thresh)
 *
 *  Build: gcc -O3 -march=native -fopenmp -o v32_anatomy src/v32_braid_anatomy.c -lm
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

typedef struct {
    double *phi[NFIELDS], *vel[NFIELDS], *acc[NFIELDS];
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
    g->N = N; g->L = L;
    g->dx = 2.0*L/(N-1); g->dt = 0.12*g->dx;
    return g;
}

/* Standard forces (for equilibration) */
static void compute_forces(Grid *g) {
    int N=g->N, NN=N*N; double idx2=1.0/(g->dx*g->dx);
    long N3=(long)N*N*N;
    #pragma omp parallel for schedule(static)
    for(long idx=0;idx<N3;idx++){
        int i=(int)(idx/NN),j=(int)((idx/N)%N),k=(int)(idx%N);
        int ip=(i+1)%N,im=(i-1+N)%N,jp=(j+1)%N,jm=(j-1+N)%N,kp=(k+1)%N,km=(k-1+N)%N;
        double p0=g->phi[0][idx],p1=g->phi[1][idx],p2=g->phi[2][idx];
        double P=p0*p1*p2,den=1.0+KAPPA*P*P,mPd2=MU*P/(den*den);
        for(int a=0;a<NFIELDS;a++){
            double lap=(g->phi[a][(long)ip*NN+j*N+k]+g->phi[a][(long)im*NN+j*N+k]
                       +g->phi[a][(long)i*NN+jp*N+k]+g->phi[a][(long)i*NN+jm*N+k]
                       +g->phi[a][(long)i*NN+j*N+kp]+g->phi[a][(long)i*NN+j*N+km]
                       -6.0*g->phi[a][idx])*idx2;
            double dPda=(a==0)?p1*p2:(a==1)?p0*p2:p0*p1;
            g->acc[a][idx]=lap-MASS2*g->phi[a][idx]-mPd2*dPda;
        }
    }
}

static void verlet_step(Grid *g) {
    long N3=(long)g->N*g->N*g->N;
    double hdt=0.5*g->dt,dt=g->dt;
    for(int a=0;a<NFIELDS;a++) for(long idx=0;idx<N3;idx++)
        g->vel[a][idx]+=hdt*g->acc[a][idx];
    for(int a=0;a<NFIELDS;a++) for(long idx=0;idx<N3;idx++)
        g->phi[a][idx]+=dt*g->vel[a][idx];
    compute_forces(g);
    for(int a=0;a<NFIELDS;a++) for(long idx=0;idx<N3;idx++)
        g->vel[a][idx]+=hdt*g->acc[a][idx];
}

static void init_braid(Grid *g) {
    int N=g->N,NN=N*N; double dx=g->dx,L=g->L;
    double A[3]={0.8,0.8,0.8},delta[3]={0,3.0005,4.4325};
    double R=3.0,el=0.3325,kw=PI/L,om=sqrt(kw*kw+MASS2);
    double sx=1+el,sy=1-el,inv2R2=1.0/(2*R*R);
    double k_bg=PI/L,om_bg=sqrt(k_bg*k_bg+MASS2);
    long N3=(long)N*N*N;
    for(int a=0;a<NFIELDS;a++){memset(g->phi[a],0,N3*sizeof(double));memset(g->vel[a],0,N3*sizeof(double));}
    for(int i=0;i<N;i++){double x=-L+i*dx;
    for(int j=0;j<N;j++){double y=-L+j*dx;
    for(int kk=0;kk<N;kk++){double z=-L+kk*dx;
        long idx=(long)i*NN+j*N+kk;
        double xr=x,yr=y;
        double r2e=xr*xr/(sx*sx)+yr*yr/(sy*sy);
        double env=exp(-r2e*inv2R2);
        for(int a=0;a<NFIELDS;a++){
            double ph=kw*z+delta[a];
            double ph_bg=k_bg*z+2*PI*a/3.0;
            g->phi[a][idx]=A[a]*env*cos(ph)+A_BG*cos(ph_bg);
            g->vel[a][idx]=om*A[a]*env*sin(ph)+om_bg*A_BG*sin(ph_bg);
        }}}}
}

/* ================================================================ */

int main(int argc, char **argv) {
    omp_set_num_threads(16);
    int N = 128; double L = 20.0;

    printf("=== V32 Braid Anatomy: Binding Structure Analysis ===\n\n");

    mkdir("data", 0755); mkdir("data/anatomy", 0755);

    Grid *g = grid_alloc(N, L);
    init_braid(g);
    compute_forces(g);

    /* Equilibrate for T=200 (let braid settle) */
    printf("Equilibrating for T=200...\n");
    int eq_steps = (int)(200.0 / g->dt);
    for (int s = 0; s < eq_steps; s++) verlet_step(g);
    printf("Done.\n\n");

    /* Now analyze the equilibrated braid at multiple time snapshots */
    int N_snaps = 10;
    double snap_interval = 5.0;
    int snap_steps = (int)(snap_interval / g->dt);

    int NN = N*N;
    double dx = g->dx;
    int k_mid = N/2;

    /* Accumulate time-averaged maps */
    int slice_N = N;
    double *avg_P    = calloc(slice_N * slice_N, sizeof(double));
    double *avg_rho  = calloc(slice_N * slice_N, sizeof(double));
    double *avg_grx  = calloc(slice_N * slice_N, sizeof(double));
    double *avg_gry  = calloc(slice_N * slice_N, sizeof(double));
    double *avg_flux_r = calloc(slice_N * slice_N, sizeof(double));
    double *avg_w    = calloc(slice_N * slice_N, sizeof(double));

    double P_max_global = 0;

    for (int snap = 0; snap < N_snaps; snap++) {
        /* Advance */
        for (int s = 0; s < snap_steps; s++) verlet_step(g);
        printf("Snap %d/%d (t=%.0f)\n", snap+1, N_snaps, 200.0 + (snap+1)*snap_interval);

        /* Analyze z=z_mid slice */
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                long idx = (long)i*NN + j*N + k_mid;
                int si = i * slice_N / N;
                int sj = j * slice_N / N;
                int sidx = si * slice_N + sj;

                double x = -L + i*dx, y = -L + j*dx;
                int ip=(i+1)%N,im=(i-1+N)%N,jp=(j+1)%N,jm=(j-1+N)%N;
                int kp=(k_mid+1)%N,km=(k_mid-1+N)%N;

                /* |P| = triple product magnitude */
                double p0=g->phi[0][idx],p1=g->phi[1][idx],p2=g->phi[2][idx];
                double P = fabs(p0*p1*p2);
                if (P > P_max_global) P_max_global = P;
                avg_P[sidx] += P;

                /* Energy density */
                double e = 0;
                for (int a=0;a<NFIELDS;a++) {
                    e += 0.5*g->vel[a][idx]*g->vel[a][idx];
                    double gx=(g->phi[a][(long)ip*NN+j*N+k_mid]-g->phi[a][(long)im*NN+j*N+k_mid])/(2*dx);
                    double gy=(g->phi[a][(long)i*NN+jp*N+k_mid]-g->phi[a][(long)i*NN+jm*N+k_mid])/(2*dx);
                    double gz=(g->phi[a][(long)i*NN+j*N+kp]-g->phi[a][(long)i*NN+j*N+km])/(2*dx);
                    e += 0.5*(gx*gx+gy*gy+gz*gz);
                    e += 0.5*MASS2*g->phi[a][idx]*g->phi[a][idx];
                }
                avg_rho[sidx] += e;

                /* Gradient of ρ (compute inline) */
                /* Need ρ at neighbors — approximate from |phi|² */
                double rho_ip=0,rho_im=0,rho_jp=0,rho_jm=0;
                for(int a=0;a<NFIELDS;a++){
                    rho_ip+=g->phi[a][(long)ip*NN+j*N+k_mid]*g->phi[a][(long)ip*NN+j*N+k_mid];
                    rho_im+=g->phi[a][(long)im*NN+j*N+k_mid]*g->phi[a][(long)im*NN+j*N+k_mid];
                    rho_jp+=g->phi[a][(long)i*NN+jp*N+k_mid]*g->phi[a][(long)i*NN+jp*N+k_mid];
                    rho_jm+=g->phi[a][(long)i*NN+jm*N+k_mid]*g->phi[a][(long)i*NN+jm*N+k_mid];
                }
                double drho_dx = (rho_ip - rho_im) / (2*dx);
                double drho_dy = (rho_jp - rho_jm) / (2*dx);
                avg_grx[sidx] += drho_dx;
                avg_gry[sidx] += drho_dy;

                /* Radial energy flux: Σ_a v_a × (r̂ · ∇φ_a) */
                double r = sqrt(x*x + y*y) + 1e-10;
                double rx = x/r, ry = y/r;
                double flux = 0;
                for (int a=0;a<NFIELDS;a++) {
                    double dpx=(g->phi[a][(long)ip*NN+j*N+k_mid]-g->phi[a][(long)im*NN+j*N+k_mid])/(2*dx);
                    double dpy=(g->phi[a][(long)i*NN+jp*N+k_mid]-g->phi[a][(long)i*NN+jm*N+k_mid])/(2*dx);
                    double radial_grad = rx*dpx + ry*dpy;
                    flux += g->vel[a][idx] * radial_grad;
                }
                avg_flux_r[sidx] += flux;
            }
        }
    }

    /* Normalize averages */
    for (int s = 0; s < slice_N*slice_N; s++) {
        avg_P[s] /= N_snaps;
        avg_rho[s] /= N_snaps;
        avg_grx[s] /= N_snaps;
        avg_gry[s] /= N_snaps;
        avg_flux_r[s] /= N_snaps;
    }

    /* Compute binding weakness: w = 1/(1 + |P|/P_thresh) */
    double P_thresh = P_max_global * 0.1;  /* 10% of max */
    printf("\nP_max = %.6f, P_threshold = %.6f\n", P_max_global, P_thresh);

    for (int s = 0; s < slice_N*slice_N; s++)
        avg_w[s] = 1.0 / (1.0 + avg_P[s] / P_thresh);

    /* Write all maps as TSV */
    FILE *fp = fopen("data/anatomy/braid_maps.tsv", "w");
    fprintf(fp, "x\ty\tP\trho\tgrad_rho_x\tgrad_rho_y\tflux_r\tbinding_weakness\n");
    for (int i = 0; i < slice_N; i++) {
        double x = -L + i * (2*L/(slice_N-1));
        for (int j = 0; j < slice_N; j++) {
            double y = -L + j * (2*L/(slice_N-1));
            int sidx = i*slice_N + j;
            fprintf(fp, "%.2f\t%.2f\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6f\n",
                    x, y, avg_P[sidx], avg_rho[sidx],
                    avg_grx[sidx], avg_gry[sidx],
                    avg_flux_r[sidx], avg_w[sidx]);
        }
    }
    fclose(fp);

    /* Radial profile of each quantity */
    printf("\n=== Radial Binding Structure ===\n");
    printf("   r     |P|      ρ        |∇ρ|      flux_r    w(binding)\n");

    FILE *fp2 = fopen("data/anatomy/radial_binding.tsv", "w");
    fprintf(fp2, "r\tP\trho\tgrad_rho\tflux_r\tw\n");

    int nbins = 40;
    double dr = L / nbins;
    for (int b = 0; b < nbins; b++) {
        double r_lo = b * dr, r_hi = (b+1) * dr;
        double sum_P=0, sum_rho=0, sum_gr=0, sum_flux=0, sum_w=0;
        int cnt = 0;
        for (int i = 0; i < slice_N; i++) {
            double x = -L + i * (2*L/(slice_N-1));
            for (int j = 0; j < slice_N; j++) {
                double y = -L + j * (2*L/(slice_N-1));
                double r = sqrt(x*x + y*y);
                if (r < r_lo || r >= r_hi) continue;
                int sidx = i*slice_N + j;
                sum_P += avg_P[sidx];
                sum_rho += avg_rho[sidx];
                sum_gr += sqrt(avg_grx[sidx]*avg_grx[sidx]+avg_gry[sidx]*avg_gry[sidx]);
                sum_flux += avg_flux_r[sidx];
                sum_w += avg_w[sidx];
                cnt++;
            }
        }
        if (cnt > 0) {
            double r = (r_lo + r_hi) / 2;
            sum_P/=cnt; sum_rho/=cnt; sum_gr/=cnt; sum_flux/=cnt; sum_w/=cnt;
            printf("  %5.1f  %.4e  %.4e  %.4e  %+.4e  %.4f\n",
                   r, sum_P, sum_rho, sum_gr, sum_flux, sum_w);
            fprintf(fp2, "%.2f\t%.6e\t%.6e\t%.6e\t%.6e\t%.6f\n",
                    r, sum_P, sum_rho, sum_gr, sum_flux, sum_w);
        }
    }
    fclose(fp2);

    printf("\n=== KEY: Where w > 0.5 is the interaction surface ===\n");
    printf("=== Where flux_r > 0 is outtake, < 0 is intake ===\n");

    free(avg_P); free(avg_rho); free(avg_grx); free(avg_gry);
    free(avg_flux_r); free(avg_w);
    printf("\n=== Complete ===\n");
    return 0;
}
