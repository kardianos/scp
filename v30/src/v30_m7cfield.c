/*  v30_m7cfield.c — M7 two-component + c(rho_B) at braid scale
 *
 *  Best model: S (braid) propagates at c_eff set by B (background).
 *  B propagates at c=1. Coupling g×S²×B converts B→S.
 *  c_eff² = 1 - alpha_c × (1 - rho_B/rho0)
 *
 *  This separates source (S) from medium (B), avoiding the freezing
 *  problem of V30's uniform c(rho).
 *
 *  Build: gcc -O3 -fopenmp -o v30_m7c src/v30_m7cfield.c -lm
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
static double G_COUP   = 0.01;    /* S-B coupling */
static double ALPHA_C  = 0.1;     /* c(rho) perturbation strength */
static double A_BG     = 0.1;     /* background amplitude */
static double RHO0_BG  = 0.0;     /* reference background density (auto-set) */

/* Grid */
typedef struct {
    double *S_phi[NFIELDS], *S_vel[NFIELDS], *S_acc[NFIELDS];
    double *B_phi[NFIELDS], *B_vel[NFIELDS], *B_acc[NFIELDS];
    double *rho_B;    /* B-field energy density */
    double *c2_eff;   /* effective c² for S fields */
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

/* ================================================================ */

static void compute_rhoB_and_ceff(Grid *g) {
    int N = g->N, NN = N*N, N3 = N*N*N;
    double dx = g->dx, mass2 = BIMODAL[14]*BIMODAL[14];

    /* B-field energy density (skip gradient for speed, use KE+mass) */
    #pragma omp parallel for schedule(static)
    for (int idx = 0; idx < N3; idx++) {
        double e = 0;
        for (int a = 0; a < NFIELDS; a++) {
            e += 0.5 * g->B_vel[a][idx] * g->B_vel[a][idx];
            e += 0.5 * mass2 * g->B_phi[a][idx] * g->B_phi[a][idx];
        }
        g->rho_B[idx] = e;
    }

    /* c_eff² = 1 - alpha_c × (1 - rho_B/rho0)
       At rho_B = rho0: c² = 1
       At rho_B = 0:    c² = 1 - alpha_c
       At rho_B > rho0: c² > 1 slightly (denser → faster) */
    #pragma omp parallel for schedule(static)
    for (int idx = 0; idx < N3; idx++) {
        double ratio = g->rho_B[idx] / (RHO0_BG + 1e-30);
        double c2 = 1.0 - ALPHA_C * (1.0 - ratio);
        if (c2 < 0.5) c2 = 0.5;   /* floor */
        if (c2 > 1.5) c2 = 1.5;   /* cap */
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

        /* S fields */
        double s0 = g->S_phi[0][idx], s1 = g->S_phi[1][idx], s2 = g->S_phi[2][idx];
        double S2 = s0*s0 + s1*s1 + s2*s2;
        double P = s0*s1*s2;
        double denom = 1.0 + kappa*P*P;
        double mu_P_d2 = mu*P/(denom*denom);

        /* B fields */
        double B2 = 0;
        for (int a = 0; a < NFIELDS; a++)
            B2 += g->B_phi[a][idx] * g->B_phi[a][idx];

        /* c_eff² for S (from B density) */
        double c2_s = g->c2_eff[idx];

        for (int a = 0; a < NFIELDS; a++) {
            /* S: Laplacian with c_eff² */
            double lap_s = (g->S_phi[a][idx+NN] + g->S_phi[a][idx-NN]
                          + g->S_phi[a][idx+N]  + g->S_phi[a][idx-N]
                          + g->S_phi[a][idx_kp] + g->S_phi[a][idx_km]
                          - 6.0*g->S_phi[a][idx]) * idx2;

            double dPda = (a==0)?s1*s2:(a==1)?s0*s2:s0*s1;

            g->S_acc[a][idx] = c2_s * lap_s       /* c_eff from B! */
                             - mass2 * g->S_phi[a][idx]
                             - mu_P_d2 * dPda
                             - G_COUP * B2 * g->S_phi[a][idx];

            /* B: Laplacian at c=1 (B is the absolute medium) */
            double lap_b = (g->B_phi[a][idx+NN] + g->B_phi[a][idx-NN]
                          + g->B_phi[a][idx+N]  + g->B_phi[a][idx-N]
                          + g->B_phi[a][idx_kp] + g->B_phi[a][idx_km]
                          - 6.0*g->B_phi[a][idx]) * idx2;

            g->B_acc[a][idx] = lap_b              /* c=1 for B always */
                             - mass2 * g->B_phi[a][idx]
                             - G_COUP * S2 * g->B_phi[a][idx];
        }
    }
}

static void verlet_step(Grid *g) {
    int N3 = g->N*g->N*g->N;
    double hdt = 0.5*g->dt, dt = g->dt;

    for (int a = 0; a < NFIELDS; a++)
        for (int idx = 0; idx < N3; idx++) {
            g->S_vel[a][idx] += hdt * g->S_acc[a][idx];
            g->B_vel[a][idx] += hdt * g->B_acc[a][idx];
        }
    for (int a = 0; a < NFIELDS; a++)
        for (int idx = 0; idx < N3; idx++) {
            g->S_phi[a][idx] += dt * g->S_vel[a][idx];
            g->B_phi[a][idx] += dt * g->B_vel[a][idx];
        }

    compute_rhoB_and_ceff(g);
    compute_forces(g);

    for (int a = 0; a < NFIELDS; a++)
        for (int idx = 0; idx < N3; idx++) {
            g->S_vel[a][idx] += hdt * g->S_acc[a][idx];
            g->B_vel[a][idx] += hdt * g->B_acc[a][idx];
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

/* ================================================================ */

int main(int argc, char **argv) {
    bimodal_init();
    omp_set_num_threads(16);

    int N = 128; double L = 20.0, T = 500.0;
    for (int i=1;i<argc;i++) {
        if(!strcmp(argv[i],"-N")&&i+1<argc) N=atoi(argv[++i]);
        else if(!strcmp(argv[i],"-L")&&i+1<argc) L=atof(argv[++i]);
        else if(!strcmp(argv[i],"-T")&&i+1<argc) T=atof(argv[++i]);
        else if(!strcmp(argv[i],"-ac")&&i+1<argc) ALPHA_C=atof(argv[++i]);
        else if(!strcmp(argv[i],"-g")&&i+1<argc) G_COUP=atof(argv[++i]);
        else if(!strcmp(argv[i],"-bg")&&i+1<argc) A_BG=atof(argv[++i]);
    }

    printf("=== V30 M7+c(rho_B): Best Model Test ===\n");
    printf("N=%d L=%.0f T=%.0f alpha_c=%.3f g=%.3f A_bg=%.2f\n\n",
           N, L, T, ALPHA_C, G_COUP, A_BG);

    Grid *g = grid_alloc(N, L);
    int NN = N*N, N3 = N*N*N;
    double dx = g->dx;
    printf("dx=%.4f dt=%.5f\n", dx, g->dt);

    /* Init S = bimodal braid */
    {
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

    /* Init B = uniform background */
    {
        double k_bg = PI/L, omega_bg = sqrt(k_bg*k_bg+BIMODAL[14]*BIMODAL[14]);
        for(int i=0;i<N;i++)for(int j=0;j<N;j++)for(int kk=0;kk<N;kk++){
            int idx=i*NN+j*N+kk; double z=-L+kk*dx;
            for(int a=0;a<NFIELDS;a++){
                double ph=k_bg*z+2*PI*a/3.0;
                g->B_phi[a][idx]=A_BG*cos(ph);
                g->B_vel[a][idx]=omega_bg*A_BG*sin(ph);
            }}
    }

    /* Set reference rho_B from initial background */
    compute_rhoB_and_ceff(g);
    {
        double sum=0; int cnt=0;
        /* Sample far-field (r>15) */
        for(int i=0;i<N;i++){double x=-L+i*dx;
        for(int j=0;j<N;j++){double y=-L+j*dx;
            if(sqrt(x*x+y*y)<15) continue;
            for(int kk=0;kk<N;kk++){
                sum+=g->rho_B[i*NN+j*N+kk]; cnt++;
            }}}
        RHO0_BG = sum/cnt;
    }
    printf("rho0_BG = %.6e (far-field B-energy density)\n", RHO0_BG);
    compute_rhoB_and_ceff(g);  /* recompute with correct rho0 */
    compute_forces(g);

    /* Scan alpha_c values: 0 (control), 0.05, 0.1, 0.2, 0.5 */
    double alphas[] = {0.0, 0.05, 0.1, 0.2, 0.5};
    int n_alphas = 5;

    mkdir("data", 0755);
    mkdir("data/m7c", 0755);

    FILE *fp_summary = fopen("data/m7c/summary.tsv", "w");
    fprintf(fp_summary, "alpha_c\tfc_final\tE_final\tmax_rho_S\tmin_c2\tavg_c2_core\tstable\n");

    for (int ai = 0; ai < n_alphas; ai++) {
        ALPHA_C = alphas[ai];
        printf("\n--- alpha_c = %.3f ---\n", ALPHA_C);

        /* Re-initialize for each alpha */
        {
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
                double xr=x*ca+y*sa,yr=-x*sa+y*ca;
                double r2e=xr*xr/(sx*sx)+yr*yr/(sy*sy);
                double env=exp(-r2e*inv2R2);
            for(int kk=0;kk<N;kk++){double z=-L+kk*dx;
                int idx=i*NN+j*N+kk;
                for(int a=0;a<NFIELDS;a++){
                    double ph=k*z+delta[a];
                    g->S_phi[a][idx]=A[a]*env*cos(ph);
                    g->S_vel[a][idx]=omega*A[a]*env*sin(ph);
                }}}}
            double k_bg=PI/L, omega_bg=sqrt(k_bg*k_bg+BIMODAL[14]*BIMODAL[14]);
            for(int i=0;i<N;i++)for(int j=0;j<N;j++)for(int kk=0;kk<N;kk++){
                int idx=i*NN+j*N+kk; double z=-L+kk*dx;
                for(int a=0;a<NFIELDS;a++){
                    double ph=k_bg*z+2*PI*a/3.0;
                    g->B_phi[a][idx]=A_BG*cos(ph);
                    g->B_vel[a][idx]=omega_bg*A_BG*sin(ph);
                    g->S_acc[a][idx]=0; g->B_acc[a][idx]=0;
                }}
        }

        compute_rhoB_and_ceff(g);
        compute_forces(g);

        int n_steps = (int)(T/g->dt);
        int diag_every = n_steps/10;
        int stable = 1;

        char tspath[256];
        snprintf(tspath,sizeof(tspath),"data/m7c/ts_ac%.3f.tsv",ALPHA_C);
        FILE *fp=fopen(tspath,"w");
        fprintf(fp,"t\tE_S\tfc_S\tmax_rhoS\tmin_c2\tavg_c2_r3\n");

        for (int step = 0; step <= n_steps; step++) {
            if (step > 0) {
                verlet_step(g);
                apply_damping(g);
            }

            if (step % diag_every == 0) {
                double t = step*g->dt;
                /* Compute S energy, fc */
                double E_S=0, phi2_core=0, phi2_total=0, max_rhoS=0;
                double min_c2=2, c2_core_sum=0; int c2_core_n=0;

                for(int idx=0;idx<N3;idx++){
                    int i=idx/NN,j=(idx/N)%N;
                    double x=-L+i*dx, y=-L+j*dx;
                    double rp2=x*x+y*y;
                    double p2=0;
                    for(int a=0;a<NFIELDS;a++){
                        p2+=g->S_phi[a][idx]*g->S_phi[a][idx];
                        E_S+=0.5*g->S_vel[a][idx]*g->S_vel[a][idx]*dx*dx*dx;
                    }
                    phi2_total+=p2;
                    if(rp2<64) phi2_core+=p2;
                    if(p2>max_rhoS) max_rhoS=p2;
                    if(g->c2_eff[idx]<min_c2) min_c2=g->c2_eff[idx];
                    if(rp2<9){c2_core_sum+=g->c2_eff[idx];c2_core_n++;}
                }
                double fc=(phi2_total>0)?phi2_core/phi2_total:0;
                double avg_c2_core=(c2_core_n>0)?c2_core_sum/c2_core_n:1;

                fprintf(fp,"%.1f\t%.1f\t%.4f\t%.4e\t%.4f\t%.4f\n",
                        t,E_S,fc,max_rhoS,min_c2,avg_c2_core);
                fflush(fp);
                printf("  t=%6.0f  fc=%.3f  E_S=%.0f  min_c2=%.3f  avg_c2(core)=%.3f\n",
                       t,fc,E_S,min_c2,avg_c2_core);

                /* Check for blowup */
                if(max_rhoS>100 || E_S!=E_S){ stable=0; printf("  BLOWUP!\n"); break; }
            }
        }
        fclose(fp);

        /* Final state */
        double fc_final=0,E_final=0,max_rhoS_final=0,min_c2_final=2,avg_c2_final=1;
        {
            double phi2_c=0,phi2_t=0;
            double c2s=0; int c2n=0;
            for(int idx=0;idx<N3;idx++){
                int i=idx/NN,j=(idx/N)%N;
                double x=-L+i*dx,y=-L+j*dx;
                double p2=0;
                for(int a=0;a<NFIELDS;a++){
                    p2+=g->S_phi[a][idx]*g->S_phi[a][idx];
                    E_final+=0.5*g->S_vel[a][idx]*g->S_vel[a][idx]*dx*dx*dx;
                }
                phi2_t+=p2;
                if(x*x+y*y<64)phi2_c+=p2;
                if(p2>max_rhoS_final)max_rhoS_final=p2;
                if(g->c2_eff[idx]<min_c2_final)min_c2_final=g->c2_eff[idx];
                if(x*x+y*y<9){c2s+=g->c2_eff[idx];c2n++;}
            }
            fc_final=(phi2_t>0)?phi2_c/phi2_t:0;
            avg_c2_final=(c2n>0)?c2s/c2n:1;
        }

        fprintf(fp_summary,"%.3f\t%.4f\t%.1f\t%.4e\t%.4f\t%.4f\t%d\n",
                ALPHA_C,fc_final,E_final,max_rhoS_final,min_c2_final,avg_c2_final,stable);
        fflush(fp_summary);
    }

    fclose(fp_summary);
    grid_free(g);
    printf("\n=== Complete ===\n");
    return 0;
}
