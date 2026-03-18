/*  t12_largescale.c — Test if depletion exponent alpha -> 1.0 at large r
 *
 *  The spectral analysis found alpha=1.19 at L=30 (r=4-20).
 *  If the deviation from 1/r is geometric (cylindrical braid, finite domain),
 *  then alpha should approach 1.0 at larger r/L.
 *
 *  Run at THREE domain sizes to extrapolate:
 *    L=30  (N=128)  — baseline, fit over r=8-25
 *    L=60  (N=128)  — medium, fit over r=15-50
 *    L=100 (N=128)  — large, fit over r=25-80
 *
 *  Same dx at all sizes (dx~0.47-1.57) — resolution varies but the
 *  FAR FIELD is what matters, where the field is smooth.
 *
 *  Method: M7 two-component (best mechanism), braid+bg vs control,
 *  spectral power measurement at probe points.
 *
 *  Build: gcc -O3 -fopenmp -o t12_large src/t12_largescale.c -lm
 */

#include "../../src/braid_core.h"
#include <float.h>

/* M7 two-component arrays */
static double *B_phi[NFIELDS], *B_vel[NFIELDS], *B_acc[NFIELDS];
static double G_COUP = 0.01;

static void compute_forces_m7(Grid *g, const double *phys) {
    int N = g->N, NN = N*N, N3 = N*N*N;
    double idx2 = 1.0 / (g->dx * g->dx);
    double mu = phys[12], kappa = phys[13];
    double mass2 = phys[14] * phys[14];
    double lpw = phys[15];

    #pragma omp parallel for schedule(static)
    for (int idx = 0; idx < N3; idx++) {
        int i = idx / NN, j = (idx / N) % N, k = idx % N;
        if (i < 1 || i >= N-1 || j < 1 || j >= N-1) {
            for (int a = 0; a < NFIELDS; a++) {
                g->acc[a][idx] = 0; B_acc[a][idx] = 0;
            }
            continue;
        }
        int kp = (k+1)%N, km = (k-1+N)%N;
        int idx_kp = i*NN+j*N+kp, idx_km = i*NN+j*N+km;

        double s0=g->phi[0][idx], s1=g->phi[1][idx], s2=g->phi[2][idx];
        double S2 = s0*s0+s1*s1+s2*s2;
        double P = s0*s1*s2;
        double denom = 1.0+kappa*P*P;
        double mu_P_d2 = mu*P/(denom*denom);
        double b0=B_phi[0][idx], b1=B_phi[1][idx], b2=B_phi[2][idx];
        double B2 = b0*b0+b1*b1+b2*b2;

        for (int a = 0; a < NFIELDS; a++) {
            double lap_s = (g->phi[a][idx+NN]+g->phi[a][idx-NN]
                          +g->phi[a][idx+N]+g->phi[a][idx-N]
                          +g->phi[a][idx_kp]+g->phi[a][idx_km]
                          -6.0*g->phi[a][idx])*idx2;
            double dPda = (a==0)?s1*s2:(a==1)?s0*s2:s0*s1;
            g->acc[a][idx] = lap_s - mass2*g->phi[a][idx]
                           - mu_P_d2*dPda
                           - lpw*(g->phi[(a+1)%3][idx]+g->phi[(a+2)%3][idx])
                           - G_COUP*B2*g->phi[a][idx];
            double lap_b = (B_phi[a][idx+NN]+B_phi[a][idx-NN]
                          +B_phi[a][idx+N]+B_phi[a][idx-N]
                          +B_phi[a][idx_kp]+B_phi[a][idx_km]
                          -6.0*B_phi[a][idx])*idx2;
            B_acc[a][idx] = lap_b - mass2*B_phi[a][idx] - G_COUP*S2*B_phi[a][idx];
        }
    }
}

static void verlet_m7(Grid *g, const double *phys) {
    int N3 = g->N*g->N*g->N;
    double hdt = 0.5*g->dt, dt = g->dt;
    for (int a = 0; a < NFIELDS; a++)
        for (int idx = 0; idx < N3; idx++) {
            g->vel[a][idx] += hdt*g->acc[a][idx];
            B_vel[a][idx]  += hdt*B_acc[a][idx];
        }
    for (int a = 0; a < NFIELDS; a++)
        for (int idx = 0; idx < N3; idx++) {
            g->phi[a][idx] += dt*g->vel[a][idx];
            B_phi[a][idx]  += dt*B_vel[a][idx];
        }
    compute_forces_m7(g, phys);
    for (int a = 0; a < NFIELDS; a++)
        for (int idx = 0; idx < N3; idx++) {
            g->vel[a][idx] += hdt*g->acc[a][idx];
            B_vel[a][idx]  += hdt*B_acc[a][idx];
        }
}

static void damp_edges_m7(Grid *g) {
    int N = g->N, NN = N*N;
    double dx = g->dx, L = g->L;
    double rs = 0.85*L, re = 0.98*L, idr = 1.0/(re-rs+1e-30);
    for (int i = 0; i < N; i++) {
        double x = -L + i*dx;
        for (int j = 0; j < N; j++) {
            double y = -L + j*dx;
            double rp = sqrt(x*x+y*y);
            if (rp <= rs) continue;
            double f = (rp-rs)*idr; if (f>1)f=1;
            double d = 1.0-0.5*f*f;
            for (int kk = 0; kk < N; kk++) {
                int idx = i*NN+j*N+kk;
                for (int a = 0; a < NFIELDS; a++) {
                    g->phi[a][idx]*=d; g->vel[a][idx]*=d;
                    B_phi[a][idx]*=d;  B_vel[a][idx]*=d;
                }
            }
        }
    }
}

/* Radial B-field energy profile (time-averaged over late time) */
typedef struct { double rho_sum; int count; } Bin;

static void measure_B_profile(Grid *g, Bin *bins, int nbins, double dr) {
    int N = g->N, NN = N*N;
    double dx = g->dx, L = g->L;
    double mass2 = BIMODAL[14]*BIMODAL[14];
    for (int b = 0; b < nbins; b++) { bins[b].rho_sum = 0; bins[b].count = 0; }
    for (int i = 1; i < N-1; i++) {
        double x = -L+i*dx;
        for (int j = 1; j < N-1; j++) {
            double y = -L+j*dx;
            double rp = sqrt(x*x+y*y);
            int b = (int)(rp/dr); if (b >= nbins) continue;
            for (int k = 0; k < N; k++) {
                int idx = i*NN+j*N+k;
                int kp=(k+1)%N, km=(k-1+N)%N;
                double e = 0;
                for (int a = 0; a < NFIELDS; a++) {
                    e += 0.5*B_vel[a][idx]*B_vel[a][idx];
                    double gx = (B_phi[a][idx+NN]-B_phi[a][idx-NN])/(2*dx);
                    double gy = (B_phi[a][idx+N]-B_phi[a][idx-N])/(2*dx);
                    double gz = (B_phi[a][i*NN+j*N+kp]-B_phi[a][i*NN+j*N+km])/(2*dx);
                    e += 0.5*(gx*gx+gy*gy+gz*gz);
                    e += 0.5*mass2*B_phi[a][idx]*B_phi[a][idx];
                }
                bins[b].rho_sum += e; bins[b].count++;
            }
        }
    }
}

/* Run one (L, N) configuration. Returns alpha from power-law fit. */
static double run_one_scale(int N, double L, double T_total, double A_bg,
                            const char *tag) {
    printf("\n============================================================\n");
    printf("  Scale: L=%.0f, N=%d, dx=%.4f, T=%.0f  [%s]\n", L, N, 2*L/(N-1), T_total, tag);
    printf("============================================================\n");

    double dx = 2.0*L/(N-1);
    double dt = 0.20*dx;
    double mass2 = BIMODAL[14]*BIMODAL[14];
    double k_bg = PI/L;
    double omega_bg = sqrt(k_bg*k_bg + mass2);
    int NN = N*N, N3 = N*N*N;

    /* ---- BRAID RUN ---- */
    Grid *g = grid_alloc(N, L);
    for (int a = 0; a < NFIELDS; a++) {
        B_phi[a] = realloc(B_phi[a], N3*sizeof(double));
        B_vel[a] = realloc(B_vel[a], N3*sizeof(double));
        B_acc[a] = realloc(B_acc[a], N3*sizeof(double));
        memset(B_acc[a], 0, N3*sizeof(double));
    }
    init_braid(g, BIMODAL, -1);
    for (int i=0;i<N;i++) for(int j=0;j<N;j++) for(int k=0;k<N;k++) {
        int idx=i*NN+j*N+k; double z=-L+k*dx;
        for(int a=0;a<NFIELDS;a++){
            double ph=k_bg*z+2.0*PI*a/3.0;
            B_phi[a][idx]=A_bg*cos(ph); B_vel[a][idx]=omega_bg*A_bg*sin(ph);
        }
    }
    compute_forces_m7(g, BIMODAL);

    int n_steps = (int)(T_total/dt);
    int avg_start = n_steps/2;
    int nbins = (int)(0.85*L/1.0)+1;  /* 1 code unit bins */
    double dr = 1.0;
    Bin *bins_braid = calloc(nbins, sizeof(Bin));
    Bin *avg_braid  = calloc(nbins, sizeof(Bin));
    int n_avg = 0;

    printf("  Braid run: %d steps...\n", n_steps); fflush(stdout);
    double t0 = omp_get_wtime();
    for (int s = 0; s < n_steps; s++) {
        verlet_m7(g, BIMODAL);
        damp_edges_m7(g);
        if (s >= avg_start && s % 50 == 0) {
            measure_B_profile(g, bins_braid, nbins, dr);
            for (int b=0;b<nbins;b++) {
                avg_braid[b].rho_sum += (bins_braid[b].count>0) ?
                    bins_braid[b].rho_sum/bins_braid[b].count : 0;
                avg_braid[b].count++;
            }
            n_avg++;
        }
        if (s%(n_steps/5)==0) {
            printf("    step %d/%d (%.0f%%, %.0fs)\n", s, n_steps,
                   100.0*s/n_steps, omp_get_wtime()-t0);
            fflush(stdout);
        }
    }
    printf("  Braid: %d avg samples, %.0fs\n", n_avg, omp_get_wtime()-t0);
    grid_free(g);

    /* ---- CONTROL RUN ---- */
    g = grid_alloc(N, L);
    grid_zero(g);
    for (int i=0;i<N;i++) for(int j=0;j<N;j++) for(int k=0;k<N;k++) {
        int idx=i*NN+j*N+k; double z=-L+k*dx;
        for(int a=0;a<NFIELDS;a++){
            double ph=k_bg*z+2.0*PI*a/3.0;
            B_phi[a][idx]=A_bg*cos(ph); B_vel[a][idx]=omega_bg*A_bg*sin(ph);
        }
    }
    memset(B_acc[0],0,N3*sizeof(double));
    memset(B_acc[1],0,N3*sizeof(double));
    memset(B_acc[2],0,N3*sizeof(double));
    compute_forces_m7(g, BIMODAL);

    Bin *bins_ctrl = calloc(nbins, sizeof(Bin));
    Bin *avg_ctrl  = calloc(nbins, sizeof(Bin));
    int n_avg_c = 0;

    printf("  Control run: %d steps...\n", n_steps); fflush(stdout);
    t0 = omp_get_wtime();
    for (int s = 0; s < n_steps; s++) {
        verlet_m7(g, BIMODAL);
        damp_edges_m7(g);
        if (s >= avg_start && s % 50 == 0) {
            measure_B_profile(g, bins_ctrl, nbins, dr);
            for (int b=0;b<nbins;b++) {
                avg_ctrl[b].rho_sum += (bins_ctrl[b].count>0) ?
                    bins_ctrl[b].rho_sum/bins_ctrl[b].count : 0;
                avg_ctrl[b].count++;
            }
            n_avg_c++;
        }
        if (s%(n_steps/5)==0) {
            printf("    step %d/%d (%.0f%%, %.0fs)\n", s, n_steps,
                   100.0*s/n_steps, omp_get_wtime()-t0);
            fflush(stdout);
        }
    }
    printf("  Control: %d avg samples, %.0fs\n", n_avg_c, omp_get_wtime()-t0);
    grid_free(g);

    /* ---- DIFFERENTIAL PROFILE AND FIT ---- */
    char fname[256];
    snprintf(fname, sizeof(fname), "data/largescale_%s_profile.tsv", tag);
    FILE *fp = fopen(fname, "w");
    fprintf(fp, "r\trho_B_braid\trho_B_ctrl\tdelta_rho\n");

    /* Fit region: r_min to r_max (skip core and edges) */
    double r_min = 0.2 * L;  /* well outside braid core */
    double r_max = 0.75 * L; /* well inside damping zone */

    double sum_lr=0, sum_lp=0, sum_lr2=0, sum_lrlp=0;
    int nfit = 0;

    printf("\n  r       rho_B_braid    rho_B_ctrl     delta_rho\n");
    for (int b = 0; b < nbins; b++) {
        double r = (b+0.5)*dr;
        if (avg_braid[b].count < 1 || avg_ctrl[b].count < 1) continue;
        double rb = avg_braid[b].rho_sum / avg_braid[b].count;
        double rc = avg_ctrl[b].rho_sum / avg_ctrl[b].count;
        double delta = rb - rc;

        fprintf(fp, "%.2f\t%.8e\t%.8e\t%+.8e\n", r, rb, rc, delta);

        if (b % (int)(5.0/dr) == 0 && r < 0.8*L) {
            printf("  %6.1f   %.6e   %.6e   %+.6e\n", r, rb, rc, delta);
        }

        /* Fit region */
        if (r >= r_min && r <= r_max && fabs(delta) > 1e-20) {
            double lr = log(r), lp = log(fabs(delta));
            sum_lr += lr; sum_lp += lp;
            sum_lr2 += lr*lr; sum_lrlp += lr*lp;
            nfit++;
        }
    }
    fclose(fp);

    double alpha = 0;
    if (nfit >= 3) {
        alpha = -(nfit*sum_lrlp - sum_lr*sum_lp) / (nfit*sum_lr2 - sum_lr*sum_lr);
        double lnA = (sum_lp + alpha*sum_lr) / nfit;
        printf("\n  FIT: |delta_rho| ~ %.4e / r^%.4f  (from %d points, r=%.0f..%.0f)\n",
               exp(lnA), alpha, nfit, r_min, r_max);
    } else {
        printf("\n  FIT: insufficient points (%d)\n", nfit);
    }

    free(bins_braid); free(avg_braid);
    free(bins_ctrl);  free(avg_ctrl);

    return alpha;
}

int main(int argc, char **argv) {
    omp_set_num_threads(16);
    bimodal_init_params();

    printf("=== T12 Large-Scale Depletion: Does alpha -> 1.0? ===\n\n");

    double A_bg = 0.1;

    /* Allocate B arrays (will be realloc'd per scale) */
    for (int a = 0; a < NFIELDS; a++) {
        B_phi[a] = NULL; B_vel[a] = NULL; B_acc[a] = NULL;
    }

    /* Three scales */
    double alpha1 = run_one_scale(128,  30.0,  500.0, A_bg, "L30");
    double alpha2 = run_one_scale(128,  60.0,  600.0, A_bg, "L60");
    double alpha3 = run_one_scale(128, 100.0,  700.0, A_bg, "L100");

    /* Summary and extrapolation */
    printf("\n============================================================\n");
    printf("  SUMMARY: Alpha vs Domain Size\n");
    printf("============================================================\n");
    printf("  L=30:   alpha = %.4f  (fit r=6..22)\n", alpha1);
    printf("  L=60:   alpha = %.4f  (fit r=12..45)\n", alpha2);
    printf("  L=100:  alpha = %.4f  (fit r=20..75)\n", alpha3);

    /* Linear extrapolation: alpha vs 1/L */
    /* If alpha = 1 + C/L, then at L->inf, alpha->1 */
    double invL1 = 1.0/30, invL2 = 1.0/60, invL3 = 1.0/100;
    double slope = (alpha3 - alpha1) / (invL3 - invL1);
    double alpha_inf = alpha1 - slope * invL1;

    printf("\n  Extrapolation: alpha(1/L) linear fit\n");
    printf("  slope = %.4f\n", slope);
    printf("  alpha(L=inf) = %.4f\n", alpha_inf);

    if (alpha_inf > 0.8 && alpha_inf < 1.2) {
        printf("\n  >>> CONSISTENT WITH 1/r (alpha->1) AT LARGE SCALE <<<\n");
    } else if (alpha_inf > 1.5) {
        printf("\n  >>> NOT approaching 1/r — depletion is steeper <<<\n");
    } else if (alpha_inf < 0.5) {
        printf("\n  >>> Shallower than 1/r — possible artifact <<<\n");
    }

    /* Write summary */
    FILE *fp = fopen("data/largescale_summary.tsv", "w");
    fprintf(fp, "L\tN\talpha\tinvL\n");
    fprintf(fp, "30\t128\t%.6f\t%.6f\n", alpha1, invL1);
    fprintf(fp, "60\t128\t%.6f\t%.6f\n", alpha2, invL2);
    fprintf(fp, "100\t128\t%.6f\t%.6f\n", alpha3, invL3);
    fprintf(fp, "# alpha_extrapolated_inf = %.6f\n", alpha_inf);
    fclose(fp);

    for (int a = 0; a < NFIELDS; a++) {
        free(B_phi[a]); free(B_vel[a]); free(B_acc[a]);
    }

    printf("\n=== T12 Large-Scale Complete ===\n");
    return 0;
}
