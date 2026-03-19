/*  v31_t2.c — Two-braid attraction test (THE GRAVITY TEST)
 *
 *  Two braids at (+-10, 0, 0) in M7+c(rho_B) model.
 *  Three configs: alpha_c=0.2, 0.0 (control), 0.5 (strong)
 *  Track separation D(t) via center-of-mass of |S|^2 in left/right halves.
 *
 *  Build: gcc -O3 -fopenmp -o v31_t2 src/v31_t2.c -lm
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

static void apply_damping(Grid *g) {
    int N = g->N, NN = N*N;
    double dx = g->dx, L = g->L;
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
                    g->B_phi[a][idx]*=d; g->B_vel[a][idx]*=d;
                }
            }
        }
    }
}

/* Initialize TWO braids at (cx1,0,0) and (cx2,0,0) */
static void init_two_braids(Grid *g, double cx1, double cx2) {
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
    for(int kk=0;kk<N;kk++){double z=-L+kk*dx;
        int idx=i*NN+j*N+kk;

        /* Braid 1 centered at (cx1, 0) */
        double x1=x-cx1, y1=y;
        double ca=cos(ell_ang),sa=sin(ell_ang);
        double xr1=x1*ca+y1*sa, yr1=-x1*sa+y1*ca;
        double r2e1=xr1*xr1/(sx*sx)+yr1*yr1/(sy*sy);
        double env1=exp(-r2e1*inv2R2);

        /* Braid 2 centered at (cx2, 0) */
        double x2=x-cx2, y2=y;
        double xr2=x2*ca+y2*sa, yr2=-x2*sa+y2*ca;
        double r2e2=xr2*xr2/(sx*sx)+yr2*yr2/(sy*sy);
        double env2=exp(-r2e2*inv2R2);

        for(int a=0;a<NFIELDS;a++){
            double ph=k*z+delta[a];
            /* Sum the two envelopes (they are well-separated) */
            g->S_phi[a][idx]=(env1+env2)*A[a]*cos(ph);
            g->S_vel[a][idx]=(env1+env2)*omega*A[a]*sin(ph);
        }
    }}}
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
        }}
}

/* Track separation: center-of-mass of |S|^2 in left (x<0) and right (x>0) halves */
static double measure_separation(Grid *g) {
    int N = g->N, NN = N*N;
    double dx = g->dx, L = g->L;
    double xL_sum=0, wL_sum=0, xR_sum=0, wR_sum=0;

    for(int i=0;i<N;i++){double x=-L+i*dx;
    for(int j=0;j<N;j++){double y=-L+j*dx;
    for(int kk=0;kk<N;kk++){
        int idx=i*NN+j*N+kk;
        double p2=0;
        for(int a=0;a<NFIELDS;a++)
            p2+=g->S_phi[a][idx]*g->S_phi[a][idx];
        if(x<0) { xL_sum+=x*p2; wL_sum+=p2; }
        else     { xR_sum+=x*p2; wR_sum+=p2; }
    }}}

    double xL = (wL_sum>0) ? xL_sum/wL_sum : -10.0;
    double xR = (wR_sum>0) ? xR_sum/wR_sum :  10.0;
    return xR - xL;
}

int main(int argc, char **argv) {
    setbuf(stdout, NULL);
    bimodal_init();
    omp_set_num_threads(16);

    int N = 128; double L = 40.0, T = 400.0;
    double sep = 10.0;  /* half-separation: braids at +-sep */

    for (int i=1;i<argc;i++) {
        if(!strcmp(argv[i],"-N")&&i+1<argc) N=atoi(argv[++i]);
        else if(!strcmp(argv[i],"-L")&&i+1<argc) L=atof(argv[++i]);
        else if(!strcmp(argv[i],"-T")&&i+1<argc) T=atof(argv[++i]);
    }

    printf("=== V31 T2: Two-Braid Attraction Test ===\n");
    printf("N=%d L=%.0f T=%.0f sep=+-%.0f\n\n", N, L, T, sep);

    mkdir("data", 0755);

    double alphas[] = {0.0, 0.2, 0.5};
    const char *labels[] = {"ac0.0", "ac0.2", "ac0.5"};
    int n_configs = 3;

    FILE *fp_sep = fopen("data/t2_separation.tsv", "w");
    fprintf(fp_sep, "config\tt\tD\txL\txR\tfc\n");

    FILE *fp_ts = fopen("data/t2_timeseries.tsv", "w");
    fprintf(fp_ts, "config\tt\tD\tE_S\tmin_c2\tavg_c2_mid\n");

    for (int ci = 0; ci < n_configs; ci++) {
        ALPHA_C = alphas[ci];
        printf("\n--- Config: %s (alpha_c=%.2f) ---\n", labels[ci], ALPHA_C);

        Grid *g = grid_alloc(N, L);
        init_two_braids(g, -sep, +sep);
        init_background(g);

        compute_rhoB_and_ceff(g);
        /* Measure rho0 from far-field */
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
        int diag_every = n_steps/50;  /* 50 snapshots */
        int NN = N*N, N3=N*N*N;
        double dx = g->dx;

        for (int step = 0; step <= n_steps; step++) {
            if (step > 0) {
                verlet_step(g);
                apply_damping(g);
            }

            if (step % diag_every == 0) {
                double t = step*g->dt;

                /* Separation */
                double xL_sum=0, wL_sum=0, xR_sum=0, wR_sum=0;
                double E_S=0, min_c2=2, c2_mid_sum=0; int c2_mid_n=0;

                for(int idx=0;idx<N3;idx++){
                    int ii=idx/NN,jj=(idx/N)%N;
                    double x=-L+ii*dx, y=-L+jj*dx;
                    double p2=0;
                    for(int a=0;a<NFIELDS;a++){
                        p2+=g->S_phi[a][idx]*g->S_phi[a][idx];
                        E_S+=0.5*g->S_vel[a][idx]*g->S_vel[a][idx]*dx*dx*dx;
                    }
                    if(x<0) { xL_sum+=x*p2; wL_sum+=p2; }
                    else     { xR_sum+=x*p2; wR_sum+=p2; }
                    if(g->c2_eff[idx]<min_c2) min_c2=g->c2_eff[idx];
                    /* midpoint region: |x|<5 */
                    if(fabs(x)<5 && fabs(y)<5) {
                        c2_mid_sum+=g->c2_eff[idx]; c2_mid_n++;
                    }
                }

                double xL = (wL_sum>0) ? xL_sum/wL_sum : -sep;
                double xR = (wR_sum>0) ? xR_sum/wR_sum :  sep;
                double D = xR - xL;
                double fc_L = wL_sum / (wL_sum+wR_sum+1e-30);
                double fc_total = 0;
                {
                    double phi2_c=0, phi2_t=0;
                    for(int idx=0;idx<N3;idx++){
                        int ii=idx/NN,jj=(idx/N)%N;
                        double x=-L+ii*dx, y=-L+jj*dx;
                        double p2=0;
                        for(int a=0;a<NFIELDS;a++)
                            p2+=g->S_phi[a][idx]*g->S_phi[a][idx];
                        phi2_t+=p2;
                        /* "core" = near either braid center */
                        double r1=(x+sep)*(x+sep)+y*y;
                        double r2=(x-sep)*(x-sep)+y*y;
                        if(r1<100||r2<100) phi2_c+=p2;
                    }
                    fc_total=(phi2_t>0)?phi2_c/phi2_t:0;
                }
                double avg_c2_mid = (c2_mid_n>0)?c2_mid_sum/c2_mid_n:1;

                fprintf(fp_sep, "%s\t%.1f\t%.4f\t%.4f\t%.4f\t%.4f\n",
                        labels[ci], t, D, xL, xR, fc_total);
                fprintf(fp_ts, "%s\t%.1f\t%.4f\t%.1f\t%.4f\t%.4f\n",
                        labels[ci], t, D, E_S, min_c2, avg_c2_mid);
                fflush(fp_sep); fflush(fp_ts);

                if (step % (n_steps/10) == 0) {
                    printf("  t=%6.0f  D=%.3f  xL=%.2f xR=%.2f  fc=%.3f  min_c2=%.3f\n",
                           t, D, xL, xR, fc_total, min_c2);
                }

                /* Blowup check */
                if(E_S!=E_S || E_S>1e10) {
                    printf("  BLOWUP at t=%.0f!\n", t);
                    break;
                }
            }
        }

        grid_free(g);
    }

    fclose(fp_sep);
    fclose(fp_ts);
    printf("\n=== T2 Complete ===\n");
    return 0;
}
