/*  v31_t1.c — Radial c_eff profile: braid+bg vs bg-only control
 *
 *  Run A: M7+c(rho_B) with braid + background (alpha_c=0.2)
 *  Run B: M7+c(rho_B) with background only (no braid, control)
 *  Differential: delta_c_eff2(r) = c_eff2_A(r) - c_eff2_B(r)
 *
 *  Build: gcc -O3 -fopenmp -o v31_t1 src/v31_t1.c -lm
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
static const double PATH_A[NDIM] = {
    0.8,0.8,0.8, 0.00,1.67, 3.0,0.80,0.0, 1.0,0.0, 0.0,0.0, -29.7,50.0,1.50,0.0};
static const double PATH_B_[NDIM] = {
    0.8,0.8,0.8, 3.53,4.92, 3.0,0.25,0.0, 1.0,0.0, 0.0,0.0, -43.4,50.0,1.50,0.0};

static void bimodal_init(void) {
    for (int d = 0; d < NDIM; d++)
        BIMODAL[d] = 0.15*PATH_A[d] + 0.85*PATH_B_[d];
}

/* Parameters */
static double G_COUP   = 0.01;
static double ALPHA_C  = 0.2;
static double A_BG     = 0.1;
static double RHO0_BG  = 0.0;

/* Grid */
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

/* Initialize S braid centered at (cx, cy) */
static void init_braid(Grid *g, double cx, double cy) {
    int N = g->N, NN = N*N;
    double dx = g->dx, L = g->L;
    double A[3]={BIMODAL[0],BIMODAL[1],BIMODAL[2]};
    double delta[3]={0,BIMODAL[3],BIMODAL[4]};
    double R_tube=BIMODAL[5], ellip=BIMODAL[6], ell_ang=BIMODAL[7];
    double k_fac=BIMODAL[8], mass=BIMODAL[14];
    double k=k_fac*PI/L, omega=sqrt(k*k+mass*mass);
    double sx=1+ellip, sy=1-ellip;
    double inv2R2=1.0/(2*R_tube*R_tube);

    for(int i=0;i<N;i++){double x=-L+i*dx-cx;
    for(int j=0;j<N;j++){double y=-L+j*dx-cy;
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

/* Initialize B background */
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

/* Measure reference rho0 from far-field */
static void measure_rho0(Grid *g) {
    int N = g->N, NN = N*N;
    double dx = g->dx, L = g->L;
    double sum=0; int cnt=0;
    for(int i=0;i<N;i++){double x=-L+i*dx;
    for(int j=0;j<N;j++){double y=-L+j*dx;
        if(sqrt(x*x+y*y)<15) continue;
        for(int kk=0;kk<N;kk++){
            sum+=g->rho_B[i*NN+j*N+kk]; cnt++;
        }}}
    RHO0_BG = sum/cnt;
}

/* Radial profile measurement */
#define NBIN 60
#define DR 0.5

int main(int argc, char **argv) {
    setbuf(stdout, NULL);
    bimodal_init();
    omp_set_num_threads(16);

    int N = 128; double L = 30.0, T = 400.0;
    ALPHA_C = 0.2;

    for (int i=1;i<argc;i++) {
        if(!strcmp(argv[i],"-N")&&i+1<argc) N=atoi(argv[++i]);
        else if(!strcmp(argv[i],"-L")&&i+1<argc) L=atof(argv[++i]);
        else if(!strcmp(argv[i],"-T")&&i+1<argc) T=atof(argv[++i]);
        else if(!strcmp(argv[i],"-ac")&&i+1<argc) ALPHA_C=atof(argv[++i]);
    }

    printf("=== V31 T1: Radial c_eff Profile ===\n");
    printf("N=%d L=%.0f T=%.0f alpha_c=%.3f\n\n", N, L, T, ALPHA_C);

    mkdir("data", 0755);

    /* Accumulators for time-averaged radial profiles */
    double rhoB_A[NBIN]={0}, c2_A[NBIN]={0}, rhoS_A[NBIN]={0};
    double rhoB_B[NBIN]={0}, c2_B[NBIN]={0};
    long   cnt_A[NBIN]={0}, cnt_B[NBIN]={0};
    int    n_avg = 0;

    /* ===== RUN A: Braid + Background ===== */
    printf("--- Run A: Braid + Background ---\n");
    {
        Grid *g = grid_alloc(N, L);
        init_braid(g, 0.0, 0.0);
        init_background(g);
        compute_rhoB_and_ceff(g);
        measure_rho0(g);
        printf("rho0_BG = %.6e\n", RHO0_BG);
        compute_rhoB_and_ceff(g);
        compute_forces(g);

        int n_steps = (int)(T/g->dt);
        int diag_every = n_steps/50;
        double T_avg_start = 300.0;
        int avg_start_step = (int)(T_avg_start/g->dt);
        double dx = g->dx;
        int NN = N*N;

        for (int step = 0; step <= n_steps; step++) {
            if (step > 0) {
                verlet_step(g);
                apply_damping(g);
            }

            /* Print progress */
            if (step % (n_steps/10) == 0) {
                double t = step*g->dt;
                double fc=0, phi2_c=0, phi2_t=0;
                for(int idx=0;idx<N*N*N;idx++){
                    int ii=idx/NN,jj=(idx/N)%N;
                    double x=-L+ii*dx, y=-L+jj*dx;
                    double p2=0;
                    for(int a=0;a<NFIELDS;a++)
                        p2+=g->S_phi[a][idx]*g->S_phi[a][idx];
                    phi2_t+=p2;
                    if(x*x+y*y<64) phi2_c+=p2;
                }
                fc=(phi2_t>0)?phi2_c/phi2_t:0;
                printf("  t=%6.0f  fc=%.3f\n", t, fc);
            }

            /* Accumulate radial profiles for T > T_avg_start */
            if (step >= avg_start_step && step % diag_every == 0) {
                n_avg++;
                for(int i=0;i<N;i++){double x=-L+i*dx;
                for(int j=0;j<N;j++){double y=-L+j*dx;
                    double rp = sqrt(x*x+y*y);
                    int bin = (int)(rp/DR);
                    if (bin >= NBIN) continue;
                    for(int kk=0;kk<N;kk++){
                        int idx=i*NN+j*N+kk;
                        rhoB_A[bin] += g->rho_B[idx];
                        c2_A[bin]   += g->c2_eff[idx];
                        double p2=0;
                        for(int a=0;a<NFIELDS;a++)
                            p2+=g->S_phi[a][idx]*g->S_phi[a][idx];
                        rhoS_A[bin] += p2;
                        cnt_A[bin]++;
                    }}}
            }
        }
        printf("Run A: %d time-averaging snapshots\n", n_avg);
        grid_free(g);
    }

    /* ===== RUN B: Background only (no braid) ===== */
    printf("\n--- Run B: Background only (control) ---\n");
    int n_avg_B = 0;
    {
        Grid *g = grid_alloc(N, L);
        /* NO braid — S fields stay zero */
        init_background(g);
        compute_rhoB_and_ceff(g);
        /* Use same rho0 as run A for fair comparison */
        compute_rhoB_and_ceff(g);
        compute_forces(g);

        int n_steps = (int)(T/g->dt);
        int diag_every = n_steps/50;
        double T_avg_start = 300.0;
        int avg_start_step = (int)(T_avg_start/g->dt);
        double dx = g->dx;
        int NN = N*N;

        for (int step = 0; step <= n_steps; step++) {
            if (step > 0) {
                verlet_step(g);
                apply_damping(g);
            }

            if (step % (n_steps/10) == 0) {
                double t = step*g->dt;
                printf("  t=%6.0f\n", t);
            }

            if (step >= avg_start_step && step % diag_every == 0) {
                n_avg_B++;
                for(int i=0;i<N;i++){double x=-L+i*dx;
                for(int j=0;j<N;j++){double y=-L+j*dx;
                    double rp = sqrt(x*x+y*y);
                    int bin = (int)(rp/DR);
                    if (bin >= NBIN) continue;
                    for(int kk=0;kk<N;kk++){
                        int idx=i*NN+j*N+kk;
                        rhoB_B[bin] += g->rho_B[idx];
                        c2_B[bin]   += g->c2_eff[idx];
                        cnt_B[bin]++;
                    }}}
            }
        }
        printf("Run B: %d time-averaging snapshots\n", n_avg_B);
        grid_free(g);
    }

    /* ===== Compute profiles and write output ===== */
    printf("\n--- Computing profiles ---\n");

    /* Average */
    for (int b = 0; b < NBIN; b++) {
        if (cnt_A[b] > 0) { rhoB_A[b]/=cnt_A[b]; c2_A[b]/=cnt_A[b]; rhoS_A[b]/=cnt_A[b]; }
        if (cnt_B[b] > 0) { rhoB_B[b]/=cnt_B[b]; c2_B[b]/=cnt_B[b]; }
    }

    FILE *fp = fopen("data/t1_profiles.tsv", "w");
    fprintf(fp, "r\trhoB_A\trhoB_B\tdelta_rhoB\tc2_A\tc2_B\tdelta_c2\trhoS\n");
    for (int b = 0; b < NBIN; b++) {
        double r = (b+0.5)*DR;
        double drho = rhoB_A[b] - rhoB_B[b];
        double dc2  = c2_A[b]   - c2_B[b];
        fprintf(fp, "%.2f\t%.6e\t%.6e\t%.6e\t%.6f\t%.6f\t%.6e\t%.6e\n",
                r, rhoB_A[b], rhoB_B[b], drho, c2_A[b], c2_B[b], dc2, rhoS_A[b]);
    }
    fclose(fp);
    printf("Wrote data/t1_profiles.tsv\n");

    /* Power-law fit for r=5..20: delta ~ A/r^alpha */
    /* Use log-log linear regression */
    double sum_lnr=0, sum_lny=0, sum_lnr2=0, sum_lnrlny=0;
    int nfit=0;
    double sum_lnr_c=0, sum_lny_c=0, sum_lnr2_c=0, sum_lnrlny_c=0;
    int nfit_c=0;

    for (int b = 0; b < NBIN; b++) {
        double r = (b+0.5)*DR;
        if (r < 5.0 || r > 20.0) continue;
        double drho = rhoB_A[b] - rhoB_B[b];
        double dc2  = c2_A[b]   - c2_B[b];

        /* Fit |delta_rhoB| */
        if (fabs(drho) > 1e-20) {
            double lr = log(r), ly = log(fabs(drho));
            sum_lnr += lr; sum_lny += ly;
            sum_lnr2 += lr*lr; sum_lnrlny += lr*ly;
            nfit++;
        }
        /* Fit |delta_c2| */
        if (fabs(dc2) > 1e-20) {
            double lr = log(r), ly = log(fabs(dc2));
            sum_lnr_c += lr; sum_lny_c += ly;
            sum_lnr2_c += lr*lr; sum_lnrlny_c += lr*ly;
            nfit_c++;
        }
    }

    FILE *fp2 = fopen("data/t1_fit.tsv", "w");
    fprintf(fp2, "quantity\texponent\tamplitude\tnpoints\tsign\n");

    if (nfit > 2) {
        double alpha_rho = -(nfit*sum_lnrlny - sum_lnr*sum_lny) /
                            (nfit*sum_lnr2 - sum_lnr*sum_lnr);
        double lnA_rho   = (sum_lny + alpha_rho*sum_lnr) / nfit;
        /* Check sign: is the delta mostly positive or negative? */
        double drho_mid = rhoB_A[(int)(10.0/DR)] - rhoB_B[(int)(10.0/DR)];
        fprintf(fp2, "delta_rhoB\t%.4f\t%.6e\t%d\t%s\n",
                alpha_rho, exp(lnA_rho), nfit, drho_mid>0?"positive":"negative");
        printf("delta_rhoB ~ %.3e / r^%.3f  (sign=%s, %d pts)\n",
               exp(lnA_rho), alpha_rho, drho_mid>0?"positive":"negative", nfit);
    }

    if (nfit_c > 2) {
        double alpha_c2 = -(nfit_c*sum_lnrlny_c - sum_lnr_c*sum_lny_c) /
                           (nfit_c*sum_lnr2_c - sum_lnr_c*sum_lnr_c);
        double lnA_c2   = (sum_lny_c + alpha_c2*sum_lnr_c) / nfit_c;
        double dc2_mid = c2_A[(int)(10.0/DR)] - c2_B[(int)(10.0/DR)];
        fprintf(fp2, "delta_c2\t%.4f\t%.6e\t%d\t%s\n",
                alpha_c2, exp(lnA_c2), nfit_c, dc2_mid>0?"positive":"negative");
        printf("delta_c2   ~ %.3e / r^%.3f  (sign=%s, %d pts)\n",
               exp(lnA_c2), alpha_c2, dc2_mid>0?"positive":"negative", nfit_c);
    }
    fclose(fp2);

    printf("\n=== T1 Complete ===\n");
    return 0;
}
