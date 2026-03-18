/*  t12_m2.c — Mechanism 2: Back-pressure from density gradient
 *  Add beta * d(rho)/dx_i to acceleration. Compute rho(x) as local energy density
 *  smoothed over 3 cells, take gradient. beta=0.01.
 *
 *  Build: gcc -O3 -fopenmp -o t12_m2 src/t12_m2.c -lm
 */

#include "../../src/braid_core.h"

#define NBINS 60
#define A_BG  0.1
#define T_TOTAL 300.0
#define SNAP_DT 50.0

static double rho0_bg, mass2_phys;

/* --- shared measurement (identical across all mN) --- */
static void compute_radial_rho(Grid *g, double *r_bins, double *rho_bins,
                                int *counts, double dr) {
    int N = g->N, NN = N*N;
    double dx = g->dx, L = g->L;
    for (int b = 0; b < NBINS; b++) { r_bins[b]=(b+0.5)*dr; rho_bins[b]=0; counts[b]=0; }
    for (int i = 1; i < N-1; i++) {
        double x = -L + i*dx;
        for (int j = 1; j < N-1; j++) {
            double y = -L + j*dx;
            double rp = sqrt(x*x+y*y);
            int b = (int)(rp/dr); if (b >= NBINS) continue;
            for (int k = 0; k < N; k++) {
                int idx = i*NN+j*N+k;
                int kp=(k+1)%N, km=(k-1+N)%N;
                double e = 0;
                for (int a = 0; a < NFIELDS; a++) {
                    e += 0.5*g->vel[a][idx]*g->vel[a][idx];
                    double gx=(g->phi[a][idx+NN]-g->phi[a][idx-NN])/(2*dx);
                    double gy=(g->phi[a][idx+N]-g->phi[a][idx-N])/(2*dx);
                    double gz=(g->phi[a][i*NN+j*N+kp]-g->phi[a][i*NN+j*N+km])/(2*dx);
                    e += 0.5*(gx*gx+gy*gy+gz*gz);
                    e += 0.5*mass2_phys*g->phi[a][idx]*g->phi[a][idx];
                }
                rho_bins[b] += e; counts[b]++;
            }
        }
    }
    for (int b = 0; b < NBINS; b++) if (counts[b]>0) rho_bins[b]/=counts[b];
}

static double get_rho_at_r(double *r_bins, double *rho_bins, int *counts, double r_target) {
    double best=0, best_dist=1e30;
    for (int b=0; b<NBINS; b++) {
        if (counts[b]<5) continue;
        double d=fabs(r_bins[b]-r_target);
        if (d<best_dist) { best_dist=d; best=rho_bins[b]; }
    }
    return best;
}

static void add_background(Grid *g) {
    int N=g->N, NN=N*N; double dx=g->dx, L=g->L;
    double k_bg=PI/L, omega_bg=sqrt(k_bg*k_bg+mass2_phys);
    for (int i=0;i<N;i++) for (int j=0;j<N;j++) for (int k=0;k<N;k++) {
        int idx=i*NN+j*N+k; double z=-L+k*dx;
        for (int a=0;a<NFIELDS;a++) {
            double phase=k_bg*z+2.0*PI*a/3.0;
            g->phi[a][idx]+=A_BG*cos(phase);
            g->vel[a][idx]+=omega_bg*A_BG*sin(phase);
        }
    }
}

static void apply_edge_damping(Grid *g) {
    int N=g->N, NN=N*N; double dx=g->dx, L=g->L;
    double r_start=0.85*L, r_end=0.98*L, inv_dr=1.0/(r_end-r_start+1e-30);
    for (int i=0;i<N;i++) { double x=-L+i*dx;
        for (int j=0;j<N;j++) { double y=-L+j*dx;
            double rp=sqrt(x*x+y*y); if (rp<=r_start) continue;
            double f=(rp-r_start)*inv_dr; if(f>1.0) f=1.0;
            double damp=1.0-0.95*f*f;
            for (int kk=0;kk<N;kk++) { int idx=i*NN+j*N+kk;
                for (int a=0;a<NFIELDS;a++) g->vel[a][idx]*=damp;
            }
        }
    }
}

static double compute_fc(Grid *g) {
    int N=g->N, NN=N*N; double dx=g->dx, L=g->L;
    double Rc2=64.0, total=0, core=0;
    for (int i=1;i<N-1;i++) { double x=-L+i*dx;
        for (int j=1;j<N-1;j++) { double y=-L+j*dx;
            double rp2=x*x+y*y;
            for (int k=0;k<N;k++) { int idx=i*NN+j*N+k;
                double rho=0;
                for (int a=0;a<NFIELDS;a++) rho+=g->phi[a][idx]*g->phi[a][idx];
                total+=rho; if(rp2<Rc2) core+=rho;
            }
        }
    }
    return core/(total+1e-30);
}

static double compute_energy(Grid *g, const double *phys) {
    int N=g->N, NN=N*N; double dx=g->dx, L=g->L, dV=dx*dx*dx;
    double mu=phys[12],kappa=phys[13],lpw=phys[15]; double E=0;
    for (int i=1;i<N-1;i++) for (int j=1;j<N-1;j++) for (int k=0;k<N;k++) {
        int idx=i*NN+j*N+k; int kp=(k+1)%N,km=(k-1+N)%N;
        double ek=0,eg=0,em=0;
        for (int a=0;a<NFIELDS;a++) {
            ek+=0.5*g->vel[a][idx]*g->vel[a][idx];
            double gx=(g->phi[a][idx+NN]-g->phi[a][idx-NN])/(2*dx);
            double gy=(g->phi[a][idx+N]-g->phi[a][idx-N])/(2*dx);
            double gz=(g->phi[a][i*NN+j*N+kp]-g->phi[a][i*NN+j*N+km])/(2*dx);
            eg+=0.5*(gx*gx+gy*gy+gz*gz);
            em+=0.5*mass2_phys*g->phi[a][idx]*g->phi[a][idx];
        }
        double p0=g->phi[0][idx],p1=g->phi[1][idx],p2=g->phi[2][idx];
        double P=p0*p1*p2; double ep=(mu/2.0)*P*P/(1.0+kappa*P*P);
        double epw=lpw*(p0*p1+p1*p2+p2*p0);
        E+=(ek+eg+em+ep+epw)*dV;
    }
    return E;
}

/* ================================================================
   M2-specific: back-pressure from density gradient
   Need: rho field (smoothed), then add beta*grad(rho) to acc
   ================================================================ */

static double *rho_field = NULL;

static void compute_rho_field(Grid *g) {
    int N = g->N, NN = N*N, N3 = N*N*N;
    if (!rho_field) rho_field = calloc(N3, sizeof(double));

    /* Raw rho = sum_a [v_a^2/2 + m^2 phi_a^2/2] */
    #pragma omp parallel for schedule(static)
    for (int idx = 0; idx < N3; idx++) {
        double r = 0;
        for (int a = 0; a < NFIELDS; a++)
            r += 0.5*g->vel[a][idx]*g->vel[a][idx]
               + 0.5*mass2_phys*g->phi[a][idx]*g->phi[a][idx];
        rho_field[idx] = r;
    }

    /* Simple smoothing: replace with average of 7 neighbors (self + 6 nearest)
       Do one pass. For boundary safety, skip edges. */
    double *tmp = calloc(N3, sizeof(double));
    memcpy(tmp, rho_field, N3*sizeof(double));
    #pragma omp parallel for schedule(static)
    for (int idx = 0; idx < N3; idx++) {
        int i = idx/(NN), j = (idx/N)%N, k = idx%N;
        if (i<1||i>=N-1||j<1||j>=N-1) { rho_field[idx]=tmp[idx]; continue; }
        int kp=(k+1)%N, km=(k-1+N)%N;
        rho_field[idx] = (tmp[idx] + tmp[idx+NN] + tmp[idx-NN]
                        + tmp[idx+N] + tmp[idx-N]
                        + tmp[i*NN+j*N+kp] + tmp[i*NN+j*N+km]) / 7.0;
    }
    free(tmp);
}

static void apply_m2_backpressure(Grid *g, double beta) {
    int N = g->N, NN = N*N;
    double dx = g->dx, idx2 = 1.0/(2.0*dx);

    compute_rho_field(g);

    #pragma omp parallel for schedule(static)
    for (int i = 2; i < N-2; i++)
        for (int j = 2; j < N-2; j++)
            for (int k = 0; k < N; k++) {
                int idx = i*NN + j*N + k;
                int kp=(k+1)%N, km=(k-1+N)%N;
                double drho_dx = (rho_field[idx+NN] - rho_field[idx-NN]) * idx2;
                double drho_dy = (rho_field[idx+N]  - rho_field[idx-N])  * idx2;
                double drho_dz = (rho_field[i*NN+j*N+kp] - rho_field[i*NN+j*N+km]) * idx2;
                /* Add gradient of rho as extra acceleration to each field component.
                   This pushes field toward higher density (away from depletion). */
                for (int a = 0; a < NFIELDS; a++) {
                    /* Approximate: acc += beta * grad(rho) dot e_a ???
                       Actually we want: the back-pressure pushes velocity,
                       so add beta*grad(rho) as force on each field.
                       But grad(rho) is a spatial vector; we add it times phi_a direction. */
                    /* Simplest: acc_a += beta * (d rho/dx * d phi_a/dx + ...) / phi_a
                       No, just: acc_a += beta * laplacian(rho) * phi_a ...
                       The proposal says: add beta * d(rho)/dx_i to acceleration.
                       Since acc is scalar per field, interpret as:
                       acc_a += beta * (grad rho . grad phi_a) / (rho + eps) */
                    double dphi_dx = (g->phi[a][idx+NN]-g->phi[a][idx-NN])*idx2;
                    double dphi_dy = (g->phi[a][idx+N]-g->phi[a][idx-N])*idx2;
                    double dphi_dz = (g->phi[a][i*NN+j*N+kp]-g->phi[a][i*NN+j*N+km])*idx2;
                    g->acc[a][idx] += beta * (drho_dx*dphi_dx + drho_dy*dphi_dy + drho_dz*dphi_dz)
                                      / (rho_field[idx] + 1e-10);
                }
            }
}

int main(void) {
    bimodal_init_params();
    mass2_phys = 2.25;
    double k_bg = PI/30.0, omega_bg = sqrt(k_bg*k_bg+mass2_phys);
    rho0_bg = 3.0*A_BG*A_BG*omega_bg*omega_bg;

    printf("=== T12 M2: Back-Pressure ===\n");
    printf("mass2=%.4f, rho0_bg=%.6f\n", mass2_phys, rho0_bg);

    Grid *g = grid_alloc(96, 30.0);
    init_braid(g, BIMODAL, -1);
    add_background(g);
    compute_forces(g, BIMODAL, mass2_phys);

    double dr = 30.0/NBINS;
    double r_bins[NBINS], rho_bins[NBINS]; int counts[NBINS];
    double prev_drho15=0, prev_t=0;

    FILE *fout = fopen("data/t12_m2_timeseries.dat","w");
    fprintf(fout,"# t fc energy drho_r10 drho_r15 drho_r20\n");

    int n_total=(int)(T_TOTAL/g->dt), snap_every=(int)(SNAP_DT/g->dt);
    double beta = 0.01;

    printf("dt=%.5f, n_total=%d\n", g->dt, n_total);

    for (int step=0; step<=n_total; step++) {
        if (step>0) {
            /* Modified Verlet: kick-drift-force-backpressure-kick */
            double hdt=0.5*g->dt;
            verlet_kick(g,hdt);
            verlet_drift(g);
            compute_forces(g, BIMODAL, mass2_phys);
            apply_m2_backpressure(g, beta);
            verlet_kick(g,hdt);
            apply_edge_damping(g);
        }

        if (step%snap_every==0) {
            double t=step*g->dt;
            compute_radial_rho(g, r_bins, rho_bins, counts, dr);
            double rho10=get_rho_at_r(r_bins,rho_bins,counts,10.0);
            double rho15=get_rho_at_r(r_bins,rho_bins,counts,15.0);
            double rho20=get_rho_at_r(r_bins,rho_bins,counts,20.0);
            double drho10=rho10-rho0_bg, drho15=rho15-rho0_bg, drho20=rho20-rho0_bg;
            double fc=compute_fc(g), E=compute_energy(g,BIMODAL);
            double rate15=(t>prev_t+1.0)?(drho15-prev_drho15)/(t-prev_t):0;

            printf("t=%6.1f  fc=%.4f  E=%.1f  drho(10)=%+.4e  drho(15)=%+.4e  drho(20)=%+.4e  d/dt=%.2e\n",
                   t,fc,E,drho10,drho15,drho20,rate15);
            fprintf(fout,"%.2f %.6f %.2f %.8e %.8e %.8e\n",t,fc,E,drho10,drho15,drho20);
            prev_drho15=drho15; prev_t=t;
            if (check_blowup(g)) { printf("BLOWUP at t=%.1f\n",t); break; }
        }
    }
    fclose(fout);
    if (rho_field) free(rho_field);
    grid_free(g);
    printf("=== M2 Complete ===\n");
    return 0;
}
