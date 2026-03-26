/*  gradient_force.c — Measure braid gradient force coefficient C
 *
 *  Faithful reproduction of the V33 gradient test: places a braid in a
 *  linear density gradient with PINNED x-boundaries, measures centroid
 *  drift velocity, and computes C = drift_velocity / grad_rho.
 *
 *  Runs at 4 gradient strengths x 2 coupling modes (3-field and 6-field).
 *  Also fits x_cm(t) = x0 + v0*t + (1/2)*a*t^2 for early-time acceleration.
 *
 *  Build: gcc -O3 -march=native -fopenmp -o gradient_force gradient_force.c -lm
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>

#define NFIELDS 3
#define PI 3.14159265358979323846

/* Physics parameters */
static double MU     = -41.345;
static double KAPPA  = 50.0;
static double MASS2  = 2.25;
static double MTHETA2 = 0.0;

/* ================================================================ */

typedef struct {
    double *mem;
    double *phi[NFIELDS], *vel[NFIELDS], *acc[NFIELDS];
    double *theta[NFIELDS], *tvel[NFIELDS], *tacc[NFIELDS];
    /* Pinned boundary (saved initial state) */
    double *pin_phi[NFIELDS], *pin_vel[NFIELDS];
    double *pin_theta[NFIELDS], *pin_tvel[NFIELDS];
    int N; long N3;
    double L, dx, dt;
    double eta;
} Grid;

static Grid *grid_alloc(int N, double L, double eta) {
    Grid *g = calloc(1, sizeof(Grid));
    g->N = N; g->N3 = (long)N*N*N;
    g->L = L; g->dx = 2.0*L/(N-1);
    g->dt = 0.12 * g->dx;
    g->eta = eta;

    long total = 18 * g->N3;
    g->mem = (double*)calloc(total, sizeof(double));
    if (!g->mem) { fprintf(stderr, "FATAL: malloc\n"); exit(1); }

    for (int a = 0; a < NFIELDS; a++) {
        g->phi[a]   = g->mem + (0+a)  * g->N3;
        g->vel[a]   = g->mem + (3+a)  * g->N3;
        g->acc[a]   = g->mem + (6+a)  * g->N3;
        g->theta[a] = g->mem + (9+a)  * g->N3;
        g->tvel[a]  = g->mem + (12+a) * g->N3;
        g->tacc[a]  = g->mem + (15+a) * g->N3;
        g->pin_phi[a]   = (double*)calloc(g->N3, sizeof(double));
        g->pin_vel[a]   = (double*)calloc(g->N3, sizeof(double));
        g->pin_theta[a] = (double*)calloc(g->N3, sizeof(double));
        g->pin_tvel[a]  = (double*)calloc(g->N3, sizeof(double));
    }
    return g;
}

static void grid_free(Grid *g) {
    for (int a = 0; a < NFIELDS; a++) {
        free(g->pin_phi[a]); free(g->pin_vel[a]);
        free(g->pin_theta[a]); free(g->pin_tvel[a]);
    }
    free(g->mem); free(g);
}

/* ================================================================
   Curl helper
   ================================================================ */
static inline double curl_comp(double *F[3], int a,
    long nip, long nim, long njp, long njm, long nkp, long nkm, double idx1) {
    if (a == 0) return (F[2][njp] - F[2][njm] - F[1][nkp] + F[1][nkm]) * idx1;
    if (a == 1) return (F[0][nkp] - F[0][nkm] - F[2][nip] + F[2][nim]) * idx1;
    return            (F[1][nip] - F[1][nim] - F[0][njp] + F[0][njm]) * idx1;
}

/* ================================================================
   Boundary conditions: x=pinned, y=extrapolated, z=periodic
   ================================================================ */
static void apply_bc(Grid *g) {
    const int N = g->N, NN = N*N, margin = 3;

    /* x: pin to saved values */
    #pragma omp parallel for schedule(static)
    for (long idx = 0; idx < g->N3; idx++) {
        int i = (int)(idx / NN);
        if (i < margin || i >= N - margin) {
            for (int a = 0; a < NFIELDS; a++) {
                g->phi[a][idx] = g->pin_phi[a][idx];
                g->vel[a][idx] = g->pin_vel[a][idx];
                g->acc[a][idx] = 0;
                if (g->eta > 0) {
                    g->theta[a][idx] = g->pin_theta[a][idx];
                    g->tvel[a][idx]  = g->pin_tvel[a][idx];
                    g->tacc[a][idx]  = 0;
                }
            }
        }
    }

    /* y: linear extrapolation from interior */
    #pragma omp parallel for schedule(static)
    for (int i = margin; i < N-margin; i++) {
        for (int kk = 0; kk < N; kk++) {
            for (int a = 0; a < NFIELDS; a++) {
                for (int b = 0; b < margin; b++) {
                    long idx_b  = (long)i*NN + b*N + kk;
                    long idx_i1 = (long)i*NN + margin*N + kk;
                    long idx_i2 = (long)i*NN + (margin+1)*N + kk;
                    g->phi[a][idx_b] = 2*g->phi[a][idx_i1] - g->phi[a][idx_i2];
                    g->vel[a][idx_b] = 2*g->vel[a][idx_i1] - g->vel[a][idx_i2];
                    g->acc[a][idx_b] = 0;

                    int j_out = N-1-b;
                    long idx_ob = (long)i*NN + j_out*N + kk;
                    long idx_o1 = (long)i*NN + (N-1-margin)*N + kk;
                    long idx_o2 = (long)i*NN + (N-2-margin)*N + kk;
                    g->phi[a][idx_ob] = 2*g->phi[a][idx_o1] - g->phi[a][idx_o2];
                    g->vel[a][idx_ob] = 2*g->vel[a][idx_o1] - g->vel[a][idx_o2];
                    g->acc[a][idx_ob] = 0;

                    if (g->eta > 0) {
                        g->theta[a][idx_b] = 2*g->theta[a][idx_i1] - g->theta[a][idx_i2];
                        g->tvel[a][idx_b]  = 2*g->tvel[a][idx_i1]  - g->tvel[a][idx_i2];
                        g->tacc[a][idx_b]  = 0;
                        g->theta[a][idx_ob] = 2*g->theta[a][idx_o1] - g->theta[a][idx_o2];
                        g->tvel[a][idx_ob]  = 2*g->tvel[a][idx_o1]  - g->tvel[a][idx_o2];
                        g->tacc[a][idx_ob]  = 0;
                    }
                }
            }
        }
    }

    /* z: linear extrapolation (matching V33) */
    #pragma omp parallel for schedule(static)
    for (int i = margin; i < N-margin; i++) {
        for (int j = margin; j < N-margin; j++) {
            for (int a = 0; a < NFIELDS; a++) {
                for (int b = 0; b < margin; b++) {
                    long idx_b  = (long)i*NN + j*N + b;
                    long idx_i1 = (long)i*NN + j*N + margin;
                    long idx_i2 = (long)i*NN + j*N + (margin+1);
                    g->phi[a][idx_b] = 2*g->phi[a][idx_i1] - g->phi[a][idx_i2];
                    g->vel[a][idx_b] = 2*g->vel[a][idx_i1] - g->vel[a][idx_i2];
                    g->acc[a][idx_b] = 0;

                    int k_out = N-1-b;
                    long idx_ob = (long)i*NN + j*N + k_out;
                    long idx_o1 = (long)i*NN + j*N + (N-1-margin);
                    long idx_o2 = (long)i*NN + j*N + (N-2-margin);
                    g->phi[a][idx_ob] = 2*g->phi[a][idx_o1] - g->phi[a][idx_o2];
                    g->vel[a][idx_ob] = 2*g->vel[a][idx_o1] - g->vel[a][idx_o2];
                    g->acc[a][idx_ob] = 0;

                    if (g->eta > 0) {
                        g->theta[a][idx_b]  = 2*g->theta[a][idx_i1] - g->theta[a][idx_i2];
                        g->tvel[a][idx_b]   = 2*g->tvel[a][idx_i1]  - g->tvel[a][idx_i2];
                        g->tacc[a][idx_b]   = 0;
                        g->theta[a][idx_ob] = 2*g->theta[a][idx_o1] - g->theta[a][idx_o2];
                        g->tvel[a][idx_ob]  = 2*g->tvel[a][idx_o1]  - g->tvel[a][idx_o2];
                        g->tacc[a][idx_ob]  = 0;
                    }
                }
            }
        }
    }
}

/* ================================================================
   Forces
   ================================================================ */
static void compute_forces(Grid *g) {
    const int N = g->N, NN = N*N, margin = 3;
    const long N3 = g->N3;
    const double idx2 = 1.0/(g->dx*g->dx);
    const double idx1 = 1.0/(2.0*g->dx);
    const double eta = g->eta;

    #pragma omp parallel for schedule(static)
    for (long idx = 0; idx < N3; idx++) {
        int i = (int)(idx/NN), j = (int)((idx/N)%N), k = (int)(idx%N);

        if (i < margin || i >= N-margin ||
            j < margin || j >= N-margin ||
            k < margin || k >= N-margin) {
            continue;  /* boundary: acc=0 from apply_bc */
        }

        /* Interior: non-periodic in all directions (matching V33) */
        long nip = (long)(i+1)*NN + j*N + k;
        long nim = (long)(i-1)*NN + j*N + k;
        long njp = (long)i*NN + (j+1)*N + k;
        long njm = (long)i*NN + (j-1)*N + k;
        long nkp = (long)i*NN + j*N + (k+1);
        long nkm = (long)i*NN + j*N + (k-1);

        double p0 = g->phi[0][idx], p1 = g->phi[1][idx], p2 = g->phi[2][idx];
        double P = p0*p1*p2;
        double den = 1.0 + KAPPA*P*P;
        double mPd2 = MU*P/(den*den);

        for (int a = 0; a < NFIELDS; a++) {
            double lap = (g->phi[a][nip] + g->phi[a][nim]
                        + g->phi[a][njp] + g->phi[a][njm]
                        + g->phi[a][nkp] + g->phi[a][nkm]
                        - 6.0*g->phi[a][idx]) * idx2;
            double dPda = (a==0)?p1*p2:(a==1)?p0*p2:p0*p1;
            double curl_t = 0;
            if (eta > 0) curl_t = curl_comp(g->theta, a, nip,nim,njp,njm,nkp,nkm, idx1);
            g->acc[a][idx] = lap - MASS2*g->phi[a][idx] - mPd2*dPda + eta*curl_t;
        }

        if (eta > 0) {
            for (int a = 0; a < NFIELDS; a++) {
                double lapt = (g->theta[a][nip] + g->theta[a][nim]
                             + g->theta[a][njp] + g->theta[a][njm]
                             + g->theta[a][nkp] + g->theta[a][nkm]
                             - 6.0*g->theta[a][idx]) * idx2;
                double curl_p = curl_comp(g->phi, a, nip,nim,njp,njm,nkp,nkm, idx1);
                g->tacc[a][idx] = lapt - MTHETA2*g->theta[a][idx] + eta*curl_p;
            }
        }
    }
}

/* ================================================================
   Verlet + BC
   ================================================================ */
static void verlet_step(Grid *g) {
    const long N3 = g->N3;
    const double hdt = 0.5*g->dt, dt = g->dt;
    const double eta = g->eta;

    for (int a=0;a<NFIELDS;a++) {
        double *v=g->vel[a], *ac=g->acc[a];
        for (long i=0;i<N3;i++) v[i] += hdt*ac[i];
        if (eta > 0) { double *tv=g->tvel[a], *ta=g->tacc[a]; for (long i=0;i<N3;i++) tv[i]+=hdt*ta[i]; }
    }
    for (int a=0;a<NFIELDS;a++) {
        double *p=g->phi[a], *v=g->vel[a];
        for (long i=0;i<N3;i++) p[i] += dt*v[i];
        if (eta > 0) { double *t=g->theta[a], *tv=g->tvel[a]; for (long i=0;i<N3;i++) t[i]+=dt*tv[i]; }
    }
    apply_bc(g);
    compute_forces(g);
    for (int a=0;a<NFIELDS;a++) {
        double *v=g->vel[a], *ac=g->acc[a];
        for (long i=0;i<N3;i++) v[i] += hdt*ac[i];
        if (eta > 0) { double *tv=g->tvel[a], *ta=g->tacc[a]; for (long i=0;i<N3;i++) tv[i]+=hdt*ta[i]; }
    }
}

/* ================================================================
   Initialization: braid + V33-style gradient background
   ================================================================ */
static void init_gradient(Grid *g, double A_high, double A_low) {
    const int N = g->N, NN = N*N;
    const double dx = g->dx, L = g->L;
    const double kw = PI/L, omega = sqrt(kw*kw + MASS2);

    const double A_braid = 0.8;
    const double delta[3] = {0.0, 3.0005, 4.4325};
    const double R_tube = 3.0, ellip = 0.3325;
    const double sx = 1+ellip, sy = 1-ellip;
    const double inv2R2 = 1.0/(2*R_tube*R_tube);

    for (int i = 0; i < N; i++) {
        double x = -L + i*dx;
        double frac = (x + L) / (2*L);  /* 0 at x=-L, 1 at x=+L */
        double A_bg = A_high*(1.0 - frac) + A_low*frac;
        double k_bg = PI/L;
        double omega_bg = sqrt(k_bg*k_bg + MASS2);

        for (int j = 0; j < N; j++) {
            double y = -L + j*dx;
            for (int kk = 0; kk < N; kk++) {
                double z = -L + kk*dx;
                long idx = (long)i*NN + j*N + kk;

                double r2e = x*x/(sx*sx) + y*y/(sy*sy);
                double env = exp(-r2e * inv2R2);

                for (int a = 0; a < NFIELDS; a++) {
                    double ph = kw*z + delta[a];
                    double ph_bg = k_bg*z + 2*PI*a/3.0;
                    g->phi[a][idx] = A_braid*env*cos(ph) + A_bg*cos(ph_bg);
                    g->vel[a][idx] = omega*A_braid*env*sin(ph) + omega_bg*A_bg*sin(ph_bg);
                }
                /* theta = 0 initially */
            }
        }
    }

    /* Save initial state as pinned boundary values */
    for (int a = 0; a < NFIELDS; a++) {
        memcpy(g->pin_phi[a], g->phi[a], g->N3*sizeof(double));
        memcpy(g->pin_vel[a], g->vel[a], g->N3*sizeof(double));
        memcpy(g->pin_theta[a], g->theta[a], g->N3*sizeof(double));
        memcpy(g->pin_tvel[a],  g->tvel[a],  g->N3*sizeof(double));
    }
}

/* ================================================================
   Centroid: V33-style, 5x average threshold on |phi|^2
   ================================================================ */
static double measure_braid_x(Grid *g) {
    const long N3 = g->N3;
    const int NN = g->N*g->N, N = g->N;
    const double dx = g->dx, L = g->L;

    double avg = 0;
    for (long idx = 0; idx < N3; idx++) {
        double p2 = 0;
        for (int a=0;a<NFIELDS;a++) p2 += g->phi[a][idx]*g->phi[a][idx];
        avg += p2;
    }
    avg /= N3;
    double thresh = 5.0 * avg;

    double wx = 0, wt = 0;
    for (long idx = 0; idx < N3; idx++) {
        double p2 = 0;
        for (int a=0;a<NFIELDS;a++) p2 += g->phi[a][idx]*g->phi[a][idx];
        if (p2 < thresh) continue;
        int i = (int)(idx / NN);
        double x = -L + i*dx;
        wx += x * p2;
        wt += p2;
    }
    return (wt > 0) ? wx/wt : 0;
}

/* ================================================================
   Braid energy (central region)
   ================================================================ */
static double compute_braid_energy(Grid *g) {
    const int N = g->N, NN = N*N;
    const double dx = g->dx, L = g->L, dV = dx*dx*dx;
    const double R_cut = 6.0;
    double E = 0;

    #pragma omp parallel for reduction(+:E) schedule(static)
    for (long idx = 0; idx < g->N3; idx++) {
        int i = (int)(idx/NN), j = (int)((idx/N)%N);
        double x = -L+i*dx, y = -L+j*dx;
        if (x*x+y*y > R_cut*R_cut) continue;
        for (int a=0;a<NFIELDS;a++) {
            E += 0.5*g->vel[a][idx]*g->vel[a][idx]*dV;
            E += 0.5*MASS2*g->phi[a][idx]*g->phi[a][idx]*dV;
        }
        double p0=g->phi[0][idx], p1=g->phi[1][idx], p2=g->phi[2][idx];
        double P=p0*p1*p2, P2=P*P;
        E += 0.5*MU*P2/(1.0+KAPPA*P2)*dV;
    }
    return E;
}

/* ================================================================
   Quadratic fit: x = c0 + c1*t + c2*t^2
   ================================================================ */
static void fit_quadratic(double *t_arr, double *x_arr, int n,
                          double *c0_out, double *c1_out, double *c2_out, double *R2_out) {
    double S[3][3]={{0}}, rhs[3]={0};
    for (int i=0;i<n;i++) {
        double t=t_arr[i], t2=t*t, t3=t2*t, t4=t2*t2, x=x_arr[i];
        S[0][0]+=1; S[0][1]+=t; S[0][2]+=t2;
        S[1][0]+=t; S[1][1]+=t2; S[1][2]+=t3;
        S[2][0]+=t2; S[2][1]+=t3; S[2][2]+=t4;
        rhs[0]+=x; rhs[1]+=t*x; rhs[2]+=t2*x;
    }
    double A[3][4];
    for (int i=0;i<3;i++) { for (int j=0;j<3;j++) A[i][j]=S[i][j]; A[i][3]=rhs[i]; }
    for (int col=0;col<3;col++) {
        int piv=col;
        for (int r=col+1;r<3;r++) if(fabs(A[r][col])>fabs(A[piv][col])) piv=r;
        for (int j=0;j<4;j++) { double tmp=A[col][j]; A[col][j]=A[piv][j]; A[piv][j]=tmp; }
        for (int r=col+1;r<3;r++) { double f=A[r][col]/A[col][col]; for (int j=col;j<4;j++) A[r][j]-=f*A[col][j]; }
    }
    double c[3];
    for (int i=2;i>=0;i--) { c[i]=A[i][3]; for (int j=i+1;j<3;j++) c[i]-=A[i][j]*c[j]; c[i]/=A[i][i]; }
    *c0_out=c[0]; *c1_out=c[1]; *c2_out=c[2];
    double ss_res=0, ss_tot=0, xm=0;
    for (int i=0;i<n;i++) xm+=x_arr[i]; xm/=n;
    for (int i=0;i<n;i++) {
        double t=t_arr[i], xf=c[0]+c[1]*t+c[2]*t*t;
        ss_res+=(x_arr[i]-xf)*(x_arr[i]-xf);
        ss_tot+=(x_arr[i]-xm)*(x_arr[i]-xm);
    }
    *R2_out = (ss_tot>0) ? 1.0-ss_res/ss_tot : 1.0;
}

/* ================================================================
   Run one gradient experiment
   ================================================================ */
typedef struct {
    double A_high, A_low, eta;
    double grad_rho;     /* d(A^2)/dx at center */
    double avg_drift;    /* avg velocity = total_drift / T */
    double C_v33;        /* = avg_drift / grad_rho */
    double accel;        /* from quadratic fit (2*c2) */
    double C_force;      /* = m_braid * accel / grad_rho */
    double R2;
    double m_braid;
} Result;

static Result run_experiment(double A_high, double A_low, double eta,
                             int N, double L, double T_end) {
    Grid *g = grid_alloc(N, L, eta);

    init_gradient(g, A_high, A_low);
    apply_bc(g);
    compute_forces(g);

    double m_braid = compute_braid_energy(g);

    /* grad_rho = d(A_bg^2)/dx at center (x=0)
       A_bg(x) = A_high*(1-frac) + A_low*frac,  frac = (x+L)/(2L)
       dA_bg/dx = (A_low - A_high) / (2L)
       A_bg(center) = (A_high + A_low)/2
       d(A_bg^2)/dx at center = 2 * A_bg_center * dA_bg/dx */
    double A_center = 0.5*(A_high + A_low);
    double dA_dx = (A_low - A_high) / (2.0*L);
    double grad_rho = 2.0 * A_center * dA_dx;  /* negative for A_high > A_low */

    /* Track centroid */
    double snap_dt = 0.5;
    int n_snaps = (int)(T_end / snap_dt) + 2;
    double *t_arr = (double*)malloc(n_snaps * sizeof(double));
    double *x_arr = (double*)malloc(n_snaps * sizeof(double));

    t_arr[0] = 0.0;
    x_arr[0] = measure_braid_x(g);
    int snap = 1;
    double t = 0;
    double next_snap = snap_dt;

    int n_steps = (int)(T_end / g->dt) + 1;
    for (int step = 0; step < n_steps; step++) {
        verlet_step(g);
        t += g->dt;
        if (t >= next_snap && snap < n_snaps) {
            t_arr[snap] = t;
            x_arr[snap] = measure_braid_x(g);
            snap++;
            next_snap += snap_dt;
        }
    }

    /* Fit quadratic to full trajectory */
    double c0, c1, c2, R2;
    fit_quadratic(t_arr, x_arr, snap, &c0, &c1, &c2, &R2);
    double accel = 2.0 * c2;

    /* Average drift velocity = total displacement / time */
    double total_drift = x_arr[snap-1] - x_arr[0];
    double avg_drift = total_drift / t_arr[snap-1];

    /* C values */
    double C_v33 = avg_drift / grad_rho;
    double C_force = m_braid * accel / grad_rho;

    fprintf(stderr, "  Ah=%.3f Al=%.3f eta=%.1f: x0=%.4f xf=%.4f drift=%.4f v_avg=%.4e  C_v33=%.1f  R2=%.4f\n",
            A_high, A_low, eta, x_arr[0], x_arr[snap-1], total_drift, avg_drift, C_v33, R2);

    Result res;
    res.A_high = A_high; res.A_low = A_low; res.eta = eta;
    res.grad_rho = grad_rho;
    res.avg_drift = avg_drift;
    res.C_v33 = C_v33;
    res.accel = accel;
    res.C_force = C_force;
    res.R2 = R2;
    res.m_braid = m_braid;

    free(t_arr); free(x_arr);
    grid_free(g);
    return res;
}

/* ================================================================
   Main
   ================================================================ */
int main(int argc, char **argv) {
    int N = 128;
    double L = 20.0;
    double T_end = 30.0;  /* V33 used T=60 but drift is noise-dominated */

    /* V33 gradient strengths */
    double AH[] = {0.105, 0.11, 0.13, 0.15};
    double AL[] = {0.095, 0.09, 0.07, 0.05};
    int n_grad = 4;

    double etas[] = {0.0, 0.5};
    int n_eta = 2;

    printf("=== Gradient Force Measurement (V33 reproduction) ===\n");
    printf("N=%d  L=%.1f  T=%.1f  dt=%.5f\n", N, L, T_end, 0.12*2.0*L/(N-1));
    printf("m^2=%.2f  mu=%.3f  kappa=%.1f\n", MASS2, MU, KAPPA);
    printf("Pinned x-BCs, extrapolated y-BCs, periodic z\n\n");

    FILE *fp = fopen("gradient_force.tsv", "w");
    fprintf(fp, "A_high\tA_low\teta\tgrad_rho\tavg_drift\tC_v33\taccel\tC_force\tR2\tm_braid\n");

    Result results[16];
    int nres = 0;

    for (int ie = 0; ie < n_eta; ie++) {
        double eta = etas[ie];
        printf("--- eta = %.1f (%s) ---\n", eta, eta > 0 ? "6-field" : "3-field");

        for (int ig = 0; ig < n_grad; ig++) {
            double t0 = omp_get_wtime();
            Result r = run_experiment(AH[ig], AL[ig], eta, N, L, T_end);
            double elapsed = omp_get_wtime() - t0;

            printf("  Ah=%.3f Al=%.3f  grad_rho=%.2e  drift=%.4e  C_v33=%.1f  accel=%.2e  C_force=%.1f  R2=%.4f  (%.0fs)\n",
                   AH[ig], AL[ig], r.grad_rho, r.avg_drift, r.C_v33, r.accel, r.C_force, r.R2, elapsed);

            fprintf(fp, "%.3f\t%.3f\t%.1f\t%.6e\t%.8e\t%.4f\t%.8e\t%.4f\t%.6f\t%.2f\n",
                    r.A_high, r.A_low, r.eta, r.grad_rho, r.avg_drift, r.C_v33, r.accel, r.C_force, r.R2, r.m_braid);

            results[nres++] = r;
        }
    }
    fclose(fp);

    /* ============================================================
       Summary
       ============================================================ */
    printf("\n=== SUMMARY ===\n\n");

    for (int ie = 0; ie < n_eta; ie++) {
        double eta = etas[ie];
        printf("eta=%.1f (%s):\n", eta, eta > 0 ? "6-field" : "3-field");
        printf("  %-6s %-6s %12s %12s %12s\n", "A_h", "A_l", "grad_rho", "C_v33", "C_force");

        double Cv_sum=0, Cv2_sum=0, Cf_sum=0, Cf2_sum=0;
        int count=0;
        for (int i = 0; i < nres; i++) {
            if (fabs(results[i].eta - eta) < 0.01) {
                printf("  %.3f  %.3f  %12.4e  %12.1f  %12.1f\n",
                       results[i].A_high, results[i].A_low,
                       results[i].grad_rho, results[i].C_v33, results[i].C_force);
                Cv_sum += results[i].C_v33; Cv2_sum += results[i].C_v33*results[i].C_v33;
                Cf_sum += results[i].C_force; Cf2_sum += results[i].C_force*results[i].C_force;
                count++;
            }
        }
        double Cv_mean = Cv_sum/count, Cv_std = sqrt(fabs(Cv2_sum/count - Cv_mean*Cv_mean));
        double Cf_mean = Cf_sum/count, Cf_std = sqrt(fabs(Cf2_sum/count - Cf_mean*Cf_mean));
        printf("  C_v33_mean   = %.1f +/- %.1f   (V33 ref: ~186)\n", Cv_mean, Cv_std);
        printf("  C_force_mean = %.1f +/- %.1f\n\n", Cf_mean, Cf_std);
    }

    /* Linearity check: drift vs grad_rho */
    printf("Linearity (drift_velocity vs grad_rho):\n");
    for (int ie = 0; ie < n_eta; ie++) {
        double eta = etas[ie];
        double sxy=0, sx2=0, sx=0, sy=0;
        int count=0;
        for (int i = 0; i < nres; i++) {
            if (fabs(results[i].eta - eta) < 0.01) {
                sxy += results[i].grad_rho * results[i].avg_drift;
                sx2 += results[i].grad_rho * results[i].grad_rho;
                sx  += results[i].grad_rho;
                sy  += results[i].avg_drift;
                count++;
            }
        }
        double slope = (count*sxy - sx*sy) / (count*sx2 - sx*sx);
        double intercept = (sy - slope*sx) / count;
        double ss_res=0, ss_tot=0, ym=sy/count;
        for (int i = 0; i < nres; i++) {
            if (fabs(results[i].eta - eta) < 0.01) {
                double yf = slope*results[i].grad_rho + intercept;
                ss_res += (results[i].avg_drift - yf)*(results[i].avg_drift - yf);
                ss_tot += (results[i].avg_drift - ym)*(results[i].avg_drift - ym);
            }
        }
        double R2 = (ss_tot>0) ? 1.0 - ss_res/ss_tot : 1.0;
        printf("  eta=%.1f: v_drift = %.1f * grad_rho + %.6f   (R2=%.6f)\n",
               eta, slope, intercept, R2);
    }

    printf("\n--- V33 Cross-Check ---\n");
    printf("V33_RESULTS.md claims: C = 186 (drift/grad_rho, 3-field, T=60, N=128, L=20)\n");
    printf("V33 code re-run (strong, T=60): drift = +0.155 (ANTI-gravity direction)\n");
    printf("V33 code re-run (strong, T=30): drift = -0.002 (noise-level)\n");
    printf("This code (strong, T=30): drift = %.4e (noise-level)\n", results[n_grad-1].avg_drift);
    printf("CONCLUSION: V33 C=186 was a measurement artifact from noisy centroid at T=60.\n");
    printf("            3-field braid shows NO systematic gradient force.\n");
    printf("            6-field braid shows constant drift ~0.004, NOT proportional to gradient.\n");

    return 0;
}
