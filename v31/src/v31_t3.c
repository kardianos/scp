/*  v31_t3.c — Pinned boundaries steady-state accretion test
 *
 *  Single braid, M7+c(rho_B), B fields PINNED at domain edges.
 *  Edges act as infinite reservoir -> steady-state depletion profile.
 *
 *  Build: gcc -O3 -fopenmp -o v31_t3 src/v31_t3.c -lm
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

static double BIMODAL[NDIM];
static const double PATH_A[NDIM] = {
    0.8,0.8,0.8, 0.00,1.67, 3.0,0.80,0.0, 1.0,0.0, 0.0,0.0, -29.7,50.0,1.50,0.0};
static const double PATH_B_[NDIM] = {
    0.8,0.8,0.8, 3.53,4.92, 3.0,0.25,0.0, 1.0,0.0, 0.0,0.0, -43.4,50.0,1.50,0.0};

static void bimodal_init(void) {
    for (int d = 0; d < NDIM; d++)
        BIMODAL[d] = 0.15*PATH_A[d] + 0.85*PATH_B_[d];
}

static double G_COUP   = 0.01;
static double ALPHA_C  = 0.2;
static double A_BG     = 0.1;
static double RHO0_BG  = 0.0;

typedef struct {
    double *S_phi[NFIELDS], *S_vel[NFIELDS], *S_acc[NFIELDS];
    double *B_phi[NFIELDS], *B_vel[NFIELDS], *B_acc[NFIELDS];
    double *rho_B;
    double *c2_eff;
    /* Pinned B initial values */
    double *B_phi_init[NFIELDS];
    double *B_vel_init[NFIELDS];
    int N; double L, dx, dt;
} Grid;

static Grid *grid_alloc(int N, double L) {
    Grid *g = calloc(1, sizeof(Grid));
    int N3 = N*N*N;
    for (int a = 0; a < NFIELDS; a++) {
        g->S_phi[a] = calloc(N3, sizeof(double));
        g->S_vel[a] = calloc(N3, sizeof(double));
        g->S_acc[a] = calloc(N3, sizeof(double));
        g->B_phi[a] = calloc(N3, sizeof(double));
        g->B_vel[a] = calloc(N3, sizeof(double));
        g->B_acc[a] = calloc(N3, sizeof(double));
        g->B_phi_init[a] = calloc(N3, sizeof(double));
        g->B_vel_init[a] = calloc(N3, sizeof(double));
    }
    g->rho_B  = calloc(N3, sizeof(double));
    g->c2_eff = calloc(N3, sizeof(double));
    g->N = N; g->L = L;
    g->dx = 2.0*L/(N-1);
    g->dt = 0.15 * g->dx;
    return g;
}

static void grid_free(Grid *g) {
    for (int a = 0; a < NFIELDS; a++) {
        free(g->S_phi[a]); free(g->S_vel[a]); free(g->S_acc[a]);
        free(g->B_phi[a]); free(g->B_vel[a]); free(g->B_acc[a]);
        free(g->B_phi_init[a]); free(g->B_vel_init[a]);
    }
    free(g->rho_B); free(g->c2_eff); free(g);
}

static void compute_rhoB_and_ceff(Grid *g) {
    int N3 = g->N*g->N*g->N;
    double mass2 = BIMODAL[14]*BIMODAL[14];

    #pragma omp parallel for schedule(static)
    for (int idx = 0; idx < N3; idx++) {
        double e = 0;
        for (int a = 0; a < NFIELDS; a++) {
            e += 0.5 * g->B_vel[a][idx] * g->B_vel[a][idx];
            e += 0.5 * mass2 * g->B_phi[a][idx] * g->B_phi[a][idx];
        }
        g->rho_B[idx] = e;
    }

    #pragma omp parallel for schedule(static)
    for (int idx = 0; idx < N3; idx++) {
        double ratio = g->rho_B[idx] / (RHO0_BG + 1e-30);
        double c2 = 1.0 - ALPHA_C * (1.0 - ratio);
        if (c2 < 0.5) c2 = 0.5;
        if (c2 > 1.5) c2 = 1.5;
        g->c2_eff[idx] = c2;
    }
}

static void compute_forces(Grid *g) {
    int N = g->N, NN = N*N, N3 = N*N*N;
    double idx2 = 1.0 / (g->dx * g->dx);
    double mu = BIMODAL[12], kappa = BIMODAL[13];
    double mass2 = BIMODAL[14] * BIMODAL[14];

    #pragma omp parallel for schedule(static)
    for (int idx = 0; idx < N3; idx++) {
        int i = idx/NN, j = (idx/N)%N, k = idx%N;

        if (i < 1 || i >= N-1 || j < 1 || j >= N-1) {
            for (int a = 0; a < NFIELDS; a++) {
                g->S_acc[a][idx] = 0; g->B_acc[a][idx] = 0;
            }
            continue;
        }

        int kp = (k+1)%N, km = (k-1+N)%N;
        int idx_kp = i*NN+j*N+kp, idx_km = i*NN+j*N+km;

        double s0 = g->S_phi[0][idx], s1 = g->S_phi[1][idx], s2 = g->S_phi[2][idx];
        double S2 = s0*s0 + s1*s1 + s2*s2;
        double P = s0*s1*s2;
        double denom = 1.0 + kappa*P*P;
        double mu_P_d2 = mu*P/(denom*denom);

        double B2 = 0;
        for (int a = 0; a < NFIELDS; a++)
            B2 += g->B_phi[a][idx] * g->B_phi[a][idx];

        double c2_s = g->c2_eff[idx];

        for (int a = 0; a < NFIELDS; a++) {
            double lap_s = (g->S_phi[a][idx+NN] + g->S_phi[a][idx-NN]
                          + g->S_phi[a][idx+N]  + g->S_phi[a][idx-N]
                          + g->S_phi[a][idx_kp] + g->S_phi[a][idx_km]
                          - 6.0*g->S_phi[a][idx]) * idx2;

            double dPda = (a==0)?s1*s2:(a==1)?s0*s2:s0*s1;

            g->S_acc[a][idx] = c2_s * lap_s
                             - mass2 * g->S_phi[a][idx]
                             - mu_P_d2 * dPda
                             - G_COUP * B2 * g->S_phi[a][idx];

            double lap_b = (g->B_phi[a][idx+NN] + g->B_phi[a][idx-NN]
                          + g->B_phi[a][idx+N]  + g->B_phi[a][idx-N]
                          + g->B_phi[a][idx_kp] + g->B_phi[a][idx_km]
                          - 6.0*g->B_phi[a][idx]) * idx2;

            g->B_acc[a][idx] = lap_b
                             - mass2 * g->B_phi[a][idx]
                             - G_COUP * S2 * g->B_phi[a][idx];
        }
    }
}

static void verlet_step(Grid *g) {
    int N3 = g->N*g->N*g->N;
    double hdt = 0.5*g->dt, dt = g->dt;

    for (int a = 0; a < NFIELDS; a++) {
        double *sv=g->S_vel[a], *sa=g->S_acc[a], *bv=g->B_vel[a], *ba=g->B_acc[a];
        #pragma omp parallel for schedule(static)
        for (int idx = 0; idx < N3; idx++) {
            sv[idx] += hdt * sa[idx];
            bv[idx] += hdt * ba[idx];
        }
    }
    for (int a = 0; a < NFIELDS; a++) {
        double *sp=g->S_phi[a], *sv=g->S_vel[a], *bp=g->B_phi[a], *bv=g->B_vel[a];
        #pragma omp parallel for schedule(static)
        for (int idx = 0; idx < N3; idx++) {
            sp[idx] += dt * sv[idx];
            bp[idx] += dt * bv[idx];
        }
    }

    compute_rhoB_and_ceff(g);
    compute_forces(g);

    for (int a = 0; a < NFIELDS; a++) {
        double *sv=g->S_vel[a], *sa=g->S_acc[a], *bv=g->B_vel[a], *ba=g->B_acc[a];
        #pragma omp parallel for schedule(static)
        for (int idx = 0; idx < N3; idx++) {
            sv[idx] += hdt * sa[idx];
            bv[idx] += hdt * ba[idx];
        }
    }
}

/* Apply damping to S at edges, and PIN B at edges */
static void apply_bc(Grid *g, double t) {
    int N = g->N, NN = N*N;
    double dx = g->dx, L = g->L;

    /* S damping at far field */
    double rs = 0.70*L, re = 0.95*L, idr = 1.0/(re-rs+1e-30);
    for (int i = 0; i < N; i++) {
        double x = -L + i*dx;
        for (int j = 0; j < N; j++) {
            double y = -L + j*dx;
            double rp = sqrt(x*x+y*y);
            if (rp <= rs) continue;
            double f = (rp-rs)*idr; if(f>1)f=1;
            double d = 1.0 - 0.98*f*f;
            for (int kk = 0; kk < N; kk++) {
                int idx = i*NN+j*N+kk;
                for (int a = 0; a < NFIELDS; a++) {
                    g->S_phi[a][idx]*=d; g->S_vel[a][idx]*=d;
                }
            }
        }
    }

    /* PIN B at edges: i<2 or i>=N-2 or j<2 or j>=N-2 */
    /* Reset B_phi and B_vel to initial background values (time-evolved analytically) */
    double k_bg = PI/L, omega_bg = sqrt(k_bg*k_bg+BIMODAL[14]*BIMODAL[14]);
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            if (i >= 2 && i < N-2 && j >= 2 && j < N-2) continue;
            for (int kk = 0; kk < N; kk++) {
                int idx = i*NN+j*N+kk;
                double z = -L + kk*dx;
                for (int a = 0; a < NFIELDS; a++) {
                    double ph = k_bg*z + 2*PI*a/3.0;
                    /* Analytic time evolution of free massive wave */
                    g->B_phi[a][idx] = A_BG * cos(ph) * cos(omega_bg*t)
                                     + A_BG * sin(ph) * sin(omega_bg*t);
                    g->B_vel[a][idx] = omega_bg * A_BG * sin(ph) * cos(omega_bg*t)
                                     - omega_bg * A_BG * cos(ph) * sin(omega_bg*t);
                    /* Simplification: cos(ph-wt) = cos(ph)cos(wt)+sin(ph)sin(wt) */
                    /* vel = omega*sin(ph-wt) = omega*(sin(ph)cos(wt)-cos(ph)sin(wt)) */
                }
            }
        }
    }
}

static void init_braid(Grid *g) {
    int N = g->N, NN = N*N;
    double dx = g->dx, L = g->L;
    double A[3]={BIMODAL[0],BIMODAL[1],BIMODAL[2]};
    double delta[3]={0,BIMODAL[3],BIMODAL[4]};
    double R_tube=BIMODAL[5], ellip=BIMODAL[6], ell_ang=BIMODAL[7];
    double k_fac=BIMODAL[8], mass=BIMODAL[14];
    double k=k_fac*PI/L, omega=sqrt(k*k+mass*mass);
    double sx=1+ellip, sy=1-ellip;
    double inv2R2=1.0/(2*R_tube*R_tube);

    for(int i=0;i<N;i++){double x=-L+i*dx;
    for(int j=0;j<N;j++){double y=-L+j*dx;
        double ca=cos(ell_ang),sa=sin(ell_ang);
        double xr=x*ca+y*sa, yr=-x*sa+y*ca;
        double r2e=xr*xr/(sx*sx)+yr*yr/(sy*sy);
        double env=exp(-r2e*inv2R2);
    for(int kk=0;kk<N;kk++){double z=-L+kk*dx;
        int idx=i*NN+j*N+kk;
        for(int a=0;a<NFIELDS;a++){
            double ph=k*z+delta[a];
            g->S_phi[a][idx]=A[a]*env*cos(ph);
            g->S_vel[a][idx]=omega*A[a]*env*sin(ph);
        }}}}
}

static void init_background(Grid *g) {
    int N = g->N, NN = N*N;
    double dx = g->dx, L = g->L;
    double k_bg = PI/L, omega_bg = sqrt(k_bg*k_bg+BIMODAL[14]*BIMODAL[14]);
    for(int i=0;i<N;i++)for(int j=0;j<N;j++)for(int kk=0;kk<N;kk++){
        int idx=i*NN+j*N+kk; double z=-L+kk*dx;
        for(int a=0;a<NFIELDS;a++){
            double ph=k_bg*z+2*PI*a/3.0;
            g->B_phi[a][idx]=A_BG*cos(ph);
            g->B_vel[a][idx]=omega_bg*A_BG*sin(ph);
            g->B_phi_init[a][idx]=g->B_phi[a][idx];
            g->B_vel_init[a][idx]=g->B_vel[a][idx];
        }}
}

/* Radial profile measurement */
#define NBIN 80
#define DR 0.5

int main(int argc, char **argv) {
    setbuf(stdout, NULL);
    bimodal_init();
    omp_set_num_threads(16);

    int N = 128; double L = 40.0, T = 500.0;
    ALPHA_C = 0.2;

    for (int i=1;i<argc;i++) {
        if(!strcmp(argv[i],"-N")&&i+1<argc) N=atoi(argv[++i]);
        else if(!strcmp(argv[i],"-L")&&i+1<argc) L=atof(argv[++i]);
        else if(!strcmp(argv[i],"-T")&&i+1<argc) T=atof(argv[++i]);
        else if(!strcmp(argv[i],"-ac")&&i+1<argc) ALPHA_C=atof(argv[++i]);
    }

    printf("=== V31 T3: Pinned Boundaries Steady State ===\n");
    printf("N=%d L=%.0f T=%.0f alpha_c=%.3f\n\n", N, L, T, ALPHA_C);

    mkdir("data", 0755);

    Grid *g = grid_alloc(N, L);
    init_braid(g);
    init_background(g);

    compute_rhoB_and_ceff(g);
    /* measure rho0 */
    {
        int NN=N*N;
        double dx=g->dx, sum=0; int cnt=0;
        for(int i=0;i<N;i++){double x=-L+i*dx;
        for(int j=0;j<N;j++){double y=-L+j*dx;
            if(sqrt(x*x+y*y)<20) continue;
            for(int kk=0;kk<N;kk++){
                sum+=g->rho_B[i*NN+j*N+kk]; cnt++;
            }}}
        RHO0_BG = sum/cnt;
    }
    printf("rho0_BG = %.6e\n", RHO0_BG);
    compute_rhoB_and_ceff(g);
    compute_forces(g);

    int n_steps = (int)(T/g->dt);
    int profile_every = (int)(100.0/g->dt); /* every T=100 */
    int diag_every = n_steps/100;
    int NN = N*N, N3=N*N*N;
    double dx = g->dx;

    FILE *fp_prof = fopen("data/t3_profiles.tsv", "w");
    fprintf(fp_prof, "snapshot\tr\trhoB\tc2_eff\trhoS\n");

    FILE *fp_acc = fopen("data/t3_accretion.tsv", "w");
    fprintf(fp_acc, "t\tE_B_total\tE_B_r15\tE_S\tfc\n");

    /* Store previous E_B inside r=15 for accretion rate */
    double prev_EB_r15 = 0;
    int snap = 0;

    for (int step = 0; step <= n_steps; step++) {
        double t = step * g->dt;

        if (step > 0) {
            verlet_step(g);
            apply_bc(g, t);
        }

        /* Radial profile at T=100,200,...,1000 */
        if (step % profile_every == 0 && step > 0) {
            snap++;
            double rhoB_bin[NBIN]={0}, c2_bin[NBIN]={0}, rhoS_bin[NBIN]={0};
            long cnt_bin[NBIN]={0};

            for(int i=0;i<N;i++){double x=-L+i*dx;
            for(int j=0;j<N;j++){double y=-L+j*dx;
                double rp = sqrt(x*x+y*y);
                int bin = (int)(rp/DR);
                if (bin >= NBIN) continue;
                for(int kk=0;kk<N;kk++){
                    int idx=i*NN+j*N+kk;
                    rhoB_bin[bin] += g->rho_B[idx];
                    c2_bin[bin]   += g->c2_eff[idx];
                    double p2=0;
                    for(int a=0;a<NFIELDS;a++)
                        p2+=g->S_phi[a][idx]*g->S_phi[a][idx];
                    rhoS_bin[bin] += p2;
                    cnt_bin[bin]++;
                }}}

            for (int b = 0; b < NBIN; b++) {
                if (cnt_bin[b] == 0) continue;
                double r = (b+0.5)*DR;
                fprintf(fp_prof, "%d\t%.2f\t%.6e\t%.6f\t%.6e\n",
                        snap, r, rhoB_bin[b]/cnt_bin[b],
                        c2_bin[b]/cnt_bin[b], rhoS_bin[b]/cnt_bin[b]);
            }
            fflush(fp_prof);
            printf("  Snapshot %d at t=%.0f\n", snap, t);
        }

        /* Accretion diagnostic */
        if (step % diag_every == 0) {
            double E_B_total=0, E_B_r15=0, E_S=0;
            double phi2_c=0, phi2_t=0;
            double mass2 = BIMODAL[14]*BIMODAL[14];

            for(int idx=0;idx<N3;idx++){
                int ii=idx/NN,jj=(idx/N)%N;
                double x=-L+ii*dx, y=-L+jj*dx;
                double rp2=x*x+y*y;

                double eB=0, eS=0, p2=0;
                for(int a=0;a<NFIELDS;a++){
                    eB += 0.5*g->B_vel[a][idx]*g->B_vel[a][idx]
                        + 0.5*mass2*g->B_phi[a][idx]*g->B_phi[a][idx];
                    eS += 0.5*g->S_vel[a][idx]*g->S_vel[a][idx];
                    p2 += g->S_phi[a][idx]*g->S_phi[a][idx];
                }
                eB *= dx*dx*dx;
                eS *= dx*dx*dx;
                E_B_total += eB;
                if(rp2 < 225) E_B_r15 += eB;  /* r<15 */
                E_S += eS;
                phi2_t += p2;
                double r1 = rp2;
                if(r1<100) phi2_c+=p2;
            }
            double fc = (phi2_t>0)?phi2_c/phi2_t:0;

            fprintf(fp_acc, "%.1f\t%.4e\t%.4e\t%.1f\t%.4f\n",
                    t, E_B_total, E_B_r15, E_S, fc);
            fflush(fp_acc);

            if (step % (n_steps/10) == 0) {
                double accr_rate = (step>0 && prev_EB_r15>0) ?
                    (E_B_r15-prev_EB_r15)/(diag_every*g->dt) : 0;
                printf("  t=%6.0f  fc=%.3f  E_B(r<15)=%.2e  accr_rate=%.2e\n",
                       t, fc, E_B_r15, accr_rate);
            }
            prev_EB_r15 = E_B_r15;

            if(E_S!=E_S || E_S>1e10) {
                printf("  BLOWUP!\n"); break;
            }
        }
    }

    fclose(fp_prof);
    fclose(fp_acc);

    printf("\n=== T3 Complete ===\n");
    grid_free(g);
    return 0;
}
