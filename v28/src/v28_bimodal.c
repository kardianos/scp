/*  v28_bimodal.c — Test bimodal superposition of Path A (gravity) + Path B (EM)
 *
 *  Hypothesis: combining the two modes creates a lower-energy, higher-binding
 *  state where both trans_l2 AND torsion flux are present.
 *
 *  Tests:
 *    1. Pure A, Pure B (controls)
 *    2. Parameter interpolation A→B at t=0.25, 0.5, 0.75
 *    3. Field superposition α*A + β*B at various (α,β)
 *    4. Hybrid: A's envelope + B's phases, B's envelope + A's phases
 *    5. Energy-minimizing scan: fine grid near best interpolation point
 *
 *  Build:  gcc -O3 -fopenmp -o v28_bimodal v28_bimodal.c -lm
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include <omp.h>

#define NFIELDS 3
#define PI      3.14159265358979323846

/* Grid parameters */
static int    SIM_N    = 80;
static double SIM_T    = 200.0;
static double DOMAIN_L = 20.0;
static int    NTHREADS = 16;

/* ================================================================
   Path A: best trans_l2 (gravity) — from v28 Tier 3 #1
   Path B: best torsion (EM+stable) — from v28 Tier 3 #2
   ================================================================ */

/* Parameter order: A1, A2, A3, delta2, delta3, R_tube, ellip, ell_ang,
                    k_fac, A_bg, R_disp, ell_rot, mu, kappa, mass, lam_pw */
#define NDIM 16

static double PATH_A[NDIM] = {
    /* Best trans_l2 = 0.344, torsion = 0.42, fc = 0.49 */
    /* From Tier2 #4: m=1.50, mu=-29.7, delta=(0.00, 1.67), ellip=0.80 */
    0.8, 0.8, 0.8,   0.00, 1.67,
    3.0, 0.80, 0.0,  1.0, 0.0,
    0.0, 0.0,        -29.7, 50.0, 1.50,
    0.0
};

static double PATH_B[NDIM] = {
    /* Best torsion = 2.21, fc = 0.87, trans_l2 = 0.13 */
    /* From Tier2 #2: m=1.50, mu=-43.4, delta=(3.53, 4.92), ellip=0.25 */
    0.8, 0.8, 0.8,   3.53, 4.92,
    3.0, 0.25, 0.0,  1.0, 0.0,
    0.0, 0.0,        -43.4, 50.0, 1.50,
    0.0
};

static const char *PNAME[NDIM] = {
    "A1","A2","A3","delta2","delta3",
    "R_tube","ellip","ell_ang","k_fac","A_bg",
    "R_disp","ell_rot","mu","kappa","mass","lam_pw"
};

/* ================================================================
   Types
   ================================================================ */

typedef struct {
    double *phi[NFIELDS];
    double *vel[NFIELDS];
    double *acc[NFIELDS];
    int N;
    double L, dx, dt;
} Grid;

typedef struct {
    char label[64];
    double fitness;
    double l2_frac;
    double transverse_l2;
    double torsion_flux;
    double fc;
    double peak_P;
    double winding;
    double energy;
    int stable;
} Result;

/* ================================================================
   Grid
   ================================================================ */

static Grid *alloc_grid(void) {
    Grid *g = calloc(1, sizeof(Grid));
    int N3 = SIM_N * SIM_N * SIM_N;
    for (int a = 0; a < NFIELDS; a++) {
        g->phi[a] = calloc(N3, sizeof(double));
        g->vel[a] = calloc(N3, sizeof(double));
        g->acc[a] = calloc(N3, sizeof(double));
    }
    g->N = SIM_N;
    g->L = DOMAIN_L;
    g->dx = 2.0 * DOMAIN_L / (SIM_N - 1);
    g->dt = 0.20 * g->dx;
    return g;
}

static void free_grid(Grid *g) {
    for (int a = 0; a < NFIELDS; a++) {
        free(g->phi[a]); free(g->vel[a]); free(g->acc[a]);
    }
    free(g);
}

/* ================================================================
   Single-config initialization (same as v28_search)
   ================================================================ */

static void init_single(Grid *g, const double *phys,
                         double *phi_out[3], double *vel_out[3]) {
    double A[3]     = {phys[0], phys[1], phys[2]};
    double delta[3] = {0.0, phys[3], phys[4]};
    double R_tube   = phys[5];
    double ellip    = phys[6];
    double ell_ang  = phys[7];
    double k_fac    = phys[8];
    double A_bg     = phys[9];
    double R_disp   = phys[10];
    double ell_rot  = phys[11];
    double mass     = phys[14];

    int N = g->N, NN = N * N;
    double dx = g->dx, L = g->L;
    double k = k_fac * PI / L;
    double omega = sqrt(k * k + mass * mass);
    double inv_2R2 = 1.0 / (2.0 * R_tube * R_tube);
    double sx = 1.0 + ellip, sy = 1.0 - ellip;

    double cx[3], cy[3], ea[3];
    for (int a = 0; a < 3; a++) {
        double ang = 2.0 * PI * a / 3.0;
        cx[a] = R_disp * cos(ang);
        cy[a] = R_disp * sin(ang);
        ea[a] = (ell_rot > 0.5) ? ell_ang + 2.0 * PI * a / 3.0 : ell_ang;
    }

    for (int i = 0; i < N; i++) {
        double x = -L + i * dx;
        for (int j = 0; j < N; j++) {
            double y = -L + j * dx;
            for (int kk = 0; kk < N; kk++) {
                double z = -L + kk * dx;
                int idx = i * NN + j * N + kk;
                for (int a = 0; a < NFIELDS; a++) {
                    double xc = x - cx[a], yc = y - cy[a];
                    double ca = cos(ea[a]), sa = sin(ea[a]);
                    double xr = xc * ca + yc * sa;
                    double yr = -xc * sa + yc * ca;
                    double r2e = xr*xr/(sx*sx) + yr*yr/(sy*sy);
                    double env = exp(-r2e * inv_2R2);
                    double ph = k * z + delta[a];
                    double amp = A[a] * env + A_bg;
                    phi_out[a][idx] = amp * cos(ph);
                    vel_out[a][idx] = omega * amp * sin(ph);
                }
            }
        }
    }
}

/* ================================================================
   Initialization modes
   ================================================================ */

/* Mode 1: Single parameter set */
static void init_params(Grid *g, const double *phys) {
    int N3 = g->N * g->N * g->N;
    for (int a = 0; a < NFIELDS; a++) {
        memset(g->phi[a], 0, N3 * sizeof(double));
        memset(g->vel[a], 0, N3 * sizeof(double));
        memset(g->acc[a], 0, N3 * sizeof(double));
    }
    init_single(g, phys, g->phi, g->vel);
}

/* Mode 2: Parameter interpolation — linearly interpolate all params */
static void init_interp(Grid *g, const double *pA, const double *pB, double t) {
    double phys[NDIM];
    for (int d = 0; d < NDIM; d++)
        phys[d] = (1.0 - t) * pA[d] + t * pB[d];
    init_params(g, phys);
}

/* Mode 3: Field superposition — α*fields_A + β*fields_B
   Uses the DYNAMICS (mu, kappa, mass) from whichever has more weight,
   or from a specified blend */
static void init_superpose(Grid *g, const double *pA, const double *pB,
                           double alpha, double beta, const double *dyn_params) {
    int N3 = g->N * g->N * g->N;
    double *phiA[3], *velA[3], *phiB[3], *velB[3];
    for (int a = 0; a < 3; a++) {
        phiA[a] = calloc(N3, sizeof(double));
        velA[a] = calloc(N3, sizeof(double));
        phiB[a] = calloc(N3, sizeof(double));
        velB[a] = calloc(N3, sizeof(double));
    }

    init_single(g, pA, phiA, velA);
    init_single(g, pB, phiB, velB);

    for (int a = 0; a < 3; a++) {
        memset(g->acc[a], 0, N3 * sizeof(double));
        for (int idx = 0; idx < N3; idx++) {
            g->phi[a][idx] = alpha * phiA[a][idx] + beta * phiB[a][idx];
            g->vel[a][idx] = alpha * velA[a][idx] + beta * velB[a][idx];
        }
    }

    for (int a = 0; a < 3; a++) {
        free(phiA[a]); free(velA[a]);
        free(phiB[a]); free(velB[a]);
    }
}

/* Mode 4: Hybrid — A's geometric envelope with B's phase structure */
static void init_hybrid(Grid *g, const double *p_env, const double *p_phase) {
    /* Take envelope params (A, R_tube, ellip, etc.) from p_env,
       phase params (delta, k_fac) from p_phase,
       dynamics (mu, kappa, mass) interpolated */
    double phys[NDIM];
    /* Envelope from p_env: A1-A3, R_tube, ellip, ell_ang, R_disp, ell_rot, A_bg */
    phys[0] = p_env[0]; phys[1] = p_env[1]; phys[2] = p_env[2];
    phys[5] = p_env[5]; phys[6] = p_env[6]; phys[7] = p_env[7];
    phys[9] = p_env[9]; phys[10] = p_env[10]; phys[11] = p_env[11];
    /* Phase from p_phase: delta2, delta3, k_fac */
    phys[3] = p_phase[3]; phys[4] = p_phase[4]; phys[8] = p_phase[8];
    /* Dynamics: average of both */
    phys[12] = 0.5 * (p_env[12] + p_phase[12]);  /* mu */
    phys[13] = 0.5 * (p_env[13] + p_phase[13]);  /* kappa */
    phys[14] = 0.5 * (p_env[14] + p_phase[14]);  /* mass */
    phys[15] = 0.5 * (p_env[15] + p_phase[15]);  /* lam_pw */

    init_params(g, phys);
}

/* ================================================================
   Force computation (copied from v28_search, with OMP)
   ================================================================ */

static void compute_forces(Grid *g, const double *phys) {
    int N = g->N, NN = N * N, N3 = N * N * N;
    double idx2 = 1.0 / (g->dx * g->dx);
    double mu    = phys[12];
    double kappa = phys[13];
    double mass2 = phys[14] * phys[14];
    double lpw   = phys[15];

    #pragma omp parallel for schedule(static)
    for (int idx = 0; idx < N3; idx++) {
        int i = idx / NN, j = (idx / N) % N, k = idx % N;
        if (i < 1 || i >= N-1 || j < 1 || j >= N-1) {
            g->acc[0][idx] = g->acc[1][idx] = g->acc[2][idx] = 0;
            continue;
        }
        int kp = (k+1)%N, km = (k-1+N)%N;
        int idx_kp = i*NN + j*N + kp;
        int idx_km = i*NN + j*N + km;

        double p0 = g->phi[0][idx], p1 = g->phi[1][idx], p2 = g->phi[2][idx];
        double P = p0 * p1 * p2;
        double denom = 1.0 + kappa * P * P;
        double mu_P_d2 = mu * P / (denom * denom);

        for (int a = 0; a < NFIELDS; a++) {
            double lap = (g->phi[a][idx+NN] + g->phi[a][idx-NN]
                        + g->phi[a][idx+N]  + g->phi[a][idx-N]
                        + g->phi[a][idx_kp] + g->phi[a][idx_km]
                        - 6.0 * g->phi[a][idx]) * idx2;
            double dPda = (a==0) ? p1*p2 : (a==1) ? p0*p2 : p0*p1;
            double f_triple = mu_P_d2 * dPda;
            double f_pw = lpw * (g->phi[(a+1)%3][idx] + g->phi[(a+2)%3][idx]);
            g->acc[a][idx] = lap - mass2 * g->phi[a][idx] - f_triple - f_pw;
        }
    }
}

/* ================================================================
   Time stepping
   ================================================================ */

static void verlet_step(Grid *g) {
    int N3 = g->N * g->N * g->N;
    double dt = g->dt, hdt = 0.5 * dt;
    for (int a = 0; a < NFIELDS; a++) {
        double *p = g->phi[a], *v = g->vel[a], *ac = g->acc[a];
        for (int idx = 0; idx < N3; idx++) {
            v[idx] += hdt * ac[idx];
            p[idx] += dt * v[idx];
        }
    }
}

static void verlet_finish(Grid *g) {
    int N3 = g->N * g->N * g->N;
    double hdt = 0.5 * g->dt;
    for (int a = 0; a < NFIELDS; a++) {
        double *v = g->vel[a], *ac = g->acc[a];
        for (int idx = 0; idx < N3; idx++)
            v[idx] += hdt * ac[idx];
    }
}

static void apply_damping(Grid *g) {
    int N = g->N, NN = N * N;
    double dx = g->dx, L = g->L;
    double r_start = 0.70 * L, r_end = 0.95 * L;
    double inv_dr = 1.0 / (r_end - r_start + 1e-30);

    for (int i = 0; i < N; i++) {
        double x = -L + i * dx;
        for (int j = 0; j < N; j++) {
            double y = -L + j * dx;
            double rp = sqrt(x*x + y*y);
            if (rp <= r_start) continue;
            double f = (rp - r_start) * inv_dr;
            if (f > 1.0) f = 1.0;
            double damp = 1.0 - 0.98 * f * f;
            for (int kk = 0; kk < N; kk++) {
                int idx = i * NN + j * N + kk;
                for (int a = 0; a < NFIELDS; a++) {
                    g->phi[a][idx] *= damp;
                    g->vel[a][idx] *= damp;
                }
            }
        }
    }
}

/* ================================================================
   Diagnostics
   ================================================================ */

static double compute_winding(Grid *g) {
    int N = g->N, NN = N * N, ic = N/2, jc = N/2;
    double total = 0;
    for (int k = 0; k < N; k++) {
        int k1 = (k + 1) % N;
        int idx0 = ic*NN + jc*N + k;
        int idx1 = ic*NN + jc*N + k1;
        double re0 = g->phi[0][idx0], im0 = g->phi[1][idx0];
        double re1 = g->phi[0][idx1], im1 = g->phi[1][idx1];
        total += atan2(im1*re0 - re1*im0, re1*re0 + im1*im0);
    }
    return total / (2.0 * PI);
}

static void compute_diagnostics(Grid *g, const double *phys, Result *res) {
    int N = g->N, NN = N*N;
    double dx = g->dx, L = g->L, dV = dx*dx*dx;
    double mu = phys[12], kappa = phys[13], mass2 = phys[14]*phys[14];
    double lpw = phys[15];
    double R_core = 8.0;
    if (R_core > L/3.0) R_core = L/3.0;
    double Rc2 = R_core * R_core;

    double phi2_total = 0, phi2_core = 0, E_total = 0, peak_P = 0;
    double Ixx = 0, Iyy = 0, Izz = 0, Ixy = 0, M0 = 0, R2sum = 0;

    for (int i = 1; i < N-1; i++) {
        double x = -L + i*dx;
        for (int j = 1; j < N-1; j++) {
            double y = -L + j*dx;
            double rp2 = x*x + y*y;
            for (int k = 0; k < N; k++) {
                double z = -L + k*dx;
                int idx = i*NN + j*N + k;
                double p0 = g->phi[0][idx], p1 = g->phi[1][idx], p2 = g->phi[2][idx];
                double rho = p0*p0 + p1*p1 + p2*p2;
                phi2_total += rho;
                if (rp2 < Rc2) phi2_core += rho;

                double P = p0*p1*p2;
                double absP = fabs(P);
                if (absP > peak_P) peak_P = absP;

                double ek = 0, eg = 0;
                for (int a = 0; a < NFIELDS; a++) {
                    ek += 0.5 * g->vel[a][idx] * g->vel[a][idx];
                    int kp = (k+1)%N, km = (k-1+N)%N;
                    double gx = (g->phi[a][idx+NN]-g->phi[a][idx-NN])/(2*dx);
                    double gy = (g->phi[a][idx+N] -g->phi[a][idx-N]) /(2*dx);
                    double gz = (g->phi[a][i*NN+j*N+kp]-g->phi[a][i*NN+j*N+km])/(2*dx);
                    eg += 0.5*(gx*gx + gy*gy + gz*gz);
                }
                double em = 0.5 * mass2 * rho;
                double ep = (mu/2.0)*P*P/(1.0+kappa*P*P);
                double epw = lpw*(p0*p1 + p1*p2 + p2*p0);
                E_total += (ek+eg+em+ep+epw) * dV;

                double r2 = rp2 + z*z;
                M0 += rho; R2sum += rho*r2;
                Ixx += rho*x*x; Iyy += rho*y*y; Izz += rho*z*z;
                Ixy += rho*x*y;
            }
        }
    }

    /* Total quadrupole */
    double Qxx = 3*Ixx-R2sum, Qyy = 3*Iyy-R2sum, Qzz = 3*Izz-R2sum, Qxy = 3*Ixy;
    double Q2 = Qxx*Qxx + Qyy*Qyy + Qzz*Qzz + 2*Qxy*Qxy;
    double R2_mean = R2sum / (M0 + 1e-30);
    res->l2_frac = sqrt(Q2) / (M0 * R2_mean + 1e-30);
    res->transverse_l2 = fabs(Ixx - Iyy) / (Ixx + Iyy + 1e-30);
    res->fc = phi2_core / (phi2_total + 1e-30);
    res->peak_P = peak_P;
    res->energy = E_total;

    /* Torsion flux (all 3 components) */
    double Phi_xy = 0, Phi_xz = 0, Phi_yz = 0;
    int k_mid = N/2, j_mid = N/2, i_mid = N/2;
    for (int i = 2; i < N-2; i++) {
        for (int j = 2; j < N-2; j++) {
            {
                int idx = i*NN + j*N + k_mid;
                double dx_p1 = (g->phi[1][idx+NN]-g->phi[1][idx-NN])/(2*dx);
                double dy_p0 = (g->phi[0][idx+N] -g->phi[0][idx-N]) /(2*dx);
                Phi_xy += 0.5*(dx_p1 - dy_p0)*dx*dx;
            }
            {
                int idx = i*NN + j_mid*N + j;
                int kp = (j+1)%N, km = (j-1+N)%N;
                double dx_p2 = (g->phi[2][idx+NN]-g->phi[2][idx-NN])/(2*dx);
                double dz_p0 = (g->phi[0][i*NN+j_mid*N+kp]-g->phi[0][i*NN+j_mid*N+km])/(2*dx);
                Phi_xz += 0.5*(dx_p2 - dz_p0)*dx*dx;
            }
            {
                int idx = i_mid*NN + i*N + j;
                int kp = (j+1)%N, km = (j-1+N)%N;
                double dy_p2 = (g->phi[2][idx+N] -g->phi[2][idx-N]) /(2*dx);
                double dz_p1 = (g->phi[1][i_mid*NN+i*N+kp]-g->phi[1][i_mid*NN+i*N+km])/(2*dx);
                Phi_yz += 0.5*(dy_p2 - dz_p1)*dx*dx;
            }
        }
    }
    res->torsion_flux = fmax(fabs(Phi_xy), fmax(fabs(Phi_xz), fabs(Phi_yz)));
    res->winding = compute_winding(g);
}

/* ================================================================
   Run simulation and measure
   ================================================================ */

static int check_blowup(Grid *g) {
    int N3 = g->N * g->N * g->N;
    for (int a = 0; a < NFIELDS; a++)
        for (int idx = 0; idx < N3; idx += 37)
            if (fabs(g->phi[a][idx]) > 50.0) return 1;
    return 0;
}

static void run_sim(Grid *g, const double *dyn_params, Result *res) {
    int Nsteps = (int)(SIM_T / g->dt);
    int check_every = Nsteps / 4;
    if (check_every < 1) check_every = 1;

    for (int step = 0; step < Nsteps; step++) {
        compute_forces(g, dyn_params);
        verlet_step(g);
        compute_forces(g, dyn_params);
        verlet_finish(g);
        apply_damping(g);

        if (step > 0 && step % check_every == 0) {
            if (check_blowup(g)) {
                res->stable = 0;
                printf("    [BLOWUP at step %d]\n", step);
                return;
            }
        }
    }
    res->stable = 1;
    compute_diagnostics(g, dyn_params, res);
}

static void print_result(const Result *r) {
    printf("  %-24s  E=%10.1f  fc=%.4f  l2=%.4f  trans=%.4f  tor=%.4f  |P|=%.4f  w=%.3f  %s\n",
           r->label, r->energy, r->fc, r->l2_frac, r->transverse_l2,
           r->torsion_flux, r->peak_P, r->winding,
           r->stable ? "OK" : "BLOWUP");
}

/* ================================================================
   Main
   ================================================================ */

int main(int argc, char **argv) {
    /* Parse optional args */
    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "-N") && i+1 < argc) SIM_N = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-T") && i+1 < argc) SIM_T = atof(argv[++i]);
        else if (!strcmp(argv[i], "-threads") && i+1 < argc) NTHREADS = atoi(argv[++i]);
    }
    omp_set_num_threads(NTHREADS);

    printf("V28 Bimodal Superposition Test\n");
    printf("Grid: N=%d, T=%.0f, L=%.0f, threads=%d\n\n", SIM_N, SIM_T, DOMAIN_L, NTHREADS);

    printf("Path A (gravity):  mu=%.1f  delta=(%.2f, %.2f)  ellip=%.2f  m=%.2f\n",
           PATH_A[12], PATH_A[3], PATH_A[4], PATH_A[6], PATH_A[14]);
    printf("Path B (EM):       mu=%.1f  delta=(%.2f, %.2f)  ellip=%.2f  m=%.2f\n\n",
           PATH_B[12], PATH_B[3], PATH_B[4], PATH_B[6], PATH_B[14]);

    Grid *g = alloc_grid();
    Result results[64];
    int nr = 0;

    /* ==================== SECTION 1: Controls ==================== */
    printf("=== Section 1: Pure Controls ===\n");

    /* Pure A */
    init_params(g, PATH_A);
    strcpy(results[nr].label, "Pure A (gravity)");
    run_sim(g, PATH_A, &results[nr]);
    print_result(&results[nr]); nr++;

    /* Pure B */
    init_params(g, PATH_B);
    strcpy(results[nr].label, "Pure B (EM)");
    run_sim(g, PATH_B, &results[nr]);
    print_result(&results[nr]); nr++;

    /* ==================== SECTION 2: Parameter Interpolation ==================== */
    printf("\n=== Section 2: Parameter Interpolation A→B ===\n");

    double t_values[] = {0.15, 0.25, 0.35, 0.50, 0.65, 0.75, 0.85};
    for (int ti = 0; ti < 7; ti++) {
        double t = t_values[ti];
        double phys[NDIM];
        for (int d = 0; d < NDIM; d++)
            phys[d] = (1.0 - t) * PATH_A[d] + t * PATH_B[d];

        init_params(g, phys);
        snprintf(results[nr].label, 64, "Interp t=%.2f", t);
        run_sim(g, phys, &results[nr]);
        print_result(&results[nr]); nr++;
    }

    /* ==================== SECTION 3: Field Superposition ==================== */
    printf("\n=== Section 3: Field Superposition α*A + β*B ===\n");

    /* Dynamics: use average of A and B */
    double dyn_avg[NDIM];
    for (int d = 0; d < NDIM; d++)
        dyn_avg[d] = 0.5 * (PATH_A[d] + PATH_B[d]);

    struct { double alpha, beta; const char *label; } super[] = {
        {0.7, 0.7, "0.7A + 0.7B (equal)"},
        {0.8, 0.4, "0.8A + 0.4B (A-heavy)"},
        {0.4, 0.8, "0.4A + 0.8B (B-heavy)"},
        {0.6, 0.6, "0.6A + 0.6B (mild)"},
        {0.9, 0.3, "0.9A + 0.3B"},
        {0.3, 0.9, "0.3A + 0.9B"},
        {1.0, 0.5, "1.0A + 0.5B (A+half B)"},
        {0.5, 1.0, "0.5A + 1.0B (B+half A)"},
    };
    int n_super = sizeof(super) / sizeof(super[0]);

    for (int si = 0; si < n_super; si++) {
        init_superpose(g, PATH_A, PATH_B, super[si].alpha, super[si].beta, dyn_avg);
        snprintf(results[nr].label, 64, "Super %s", super[si].label);
        run_sim(g, dyn_avg, &results[nr]);
        print_result(&results[nr]); nr++;
    }

    /* ==================== SECTION 4: Hybrid (envelope/phase swap) ==================== */
    printf("\n=== Section 4: Hybrid Configurations ===\n");

    /* A's envelope (ellip=0.80) + B's phases (delta=3.53, 4.92) */
    init_hybrid(g, PATH_A, PATH_B);
    snprintf(results[nr].label, 64, "A-env + B-phase");
    double phys_hybrid1[NDIM];
    for (int d = 0; d < NDIM; d++) phys_hybrid1[d] = dyn_avg[d];
    phys_hybrid1[14] = 1.5; /* keep mass=1.5 */
    run_sim(g, phys_hybrid1, &results[nr]);
    print_result(&results[nr]); nr++;

    /* B's envelope (ellip=0.25) + A's phases (delta=0.00, 1.67) */
    init_hybrid(g, PATH_B, PATH_A);
    snprintf(results[nr].label, 64, "B-env + A-phase");
    run_sim(g, phys_hybrid1, &results[nr]);
    print_result(&results[nr]); nr++;

    /* A's envelope + average phases */
    {
        double phys[NDIM];
        memcpy(phys, PATH_A, sizeof(phys));
        phys[3] = 0.5*(PATH_A[3]+PATH_B[3]);  /* avg delta2 */
        phys[4] = 0.5*(PATH_A[4]+PATH_B[4]);  /* avg delta3 */
        phys[12] = 0.5*(PATH_A[12]+PATH_B[12]); /* avg mu */
        init_params(g, phys);
        snprintf(results[nr].label, 64, "A-env + avg-phase");
        run_sim(g, phys, &results[nr]);
        print_result(&results[nr]); nr++;
    }

    /* ==================== SECTION 5: Fine-scan t ∈ [0.75, 0.95] ==================== */
    printf("\n=== Section 5: Fine-scan around t=0.85 ===\n");

    for (int ti = 0; ti <= 20; ti++) {
        double t = 0.75 + ti * 0.01;
        double phys[NDIM];
        for (int d = 0; d < NDIM; d++)
            phys[d] = (1.0 - t) * PATH_A[d] + t * PATH_B[d];

        init_params(g, phys);
        snprintf(results[nr].label, 64, "Fine t=%.2f", t);
        run_sim(g, phys, &results[nr]);
        print_result(&results[nr]); nr++;
    }

    /* ==================== SECTION 6: N=128 T=500 Validation of best ==================== */
    printf("\n=== Section 6: High-res validation (N=128, T=500) ===\n");
    {
        /* Find best from fine-scan (sections 5) based on combined score */
        int best_fine = -1;
        double best_score = -1;
        /* Fine-scan results start after sections 1-4 (2+7+8+3 = 20 results) */
        int fine_start = 20;  /* index where fine-scan begins */
        for (int i = fine_start; i < nr; i++) {
            if (!results[i].stable) continue;
            double sc = 2.0 * results[i].transverse_l2
                      + 1.5 * fmin(results[i].torsion_flux / 1.0, 1.0)
                      + results[i].fc
                      + 0.5 * fmin(results[i].peak_P / 2.0, 1.0);
            if (sc > best_score) { best_score = sc; best_fine = i; }
        }

        /* Also validate pure A, pure B, and t=0.85 at high res */
        typedef struct { const char *label; double t; } ValConfig;
        ValConfig val_configs[] = {
            {"HiRes Pure A", 0.0},
            {"HiRes Pure B", 1.0},
            {"HiRes t=0.85", 0.85},
        };
        int n_val = 3;

        /* Find best_fine's t value from label */
        double best_t = 0.85;
        if (best_fine >= 0) {
            sscanf(results[best_fine].label, "Fine t=%lf", &best_t);
            if (fabs(best_t - 0.85) > 0.005) n_val = 4; /* add if different from 0.85 */
        }

        /* Reallocate grid for N=128 */
        free_grid(g);
        int save_N = SIM_N;
        double save_T = SIM_T;
        SIM_N = 128;
        SIM_T = 500.0;
        g = alloc_grid();

        for (int vi = 0; vi < n_val; vi++) {
            double t;
            const char *lbl;
            if (vi < 3) {
                t = val_configs[vi].t;
                lbl = val_configs[vi].label;
            } else {
                t = best_t;
                lbl = "HiRes Best-Fine";
            }

            double phys[NDIM];
            for (int d = 0; d < NDIM; d++)
                phys[d] = (1.0 - t) * PATH_A[d] + t * PATH_B[d];

            init_params(g, phys);
            snprintf(results[nr].label, 64, "%s (t=%.2f)", lbl, t);
            printf("  Running %s ...\n", results[nr].label); fflush(stdout);
            run_sim(g, phys, &results[nr]);
            print_result(&results[nr]); nr++;
        }

        /* Restore */
        free_grid(g);
        SIM_N = save_N;
        SIM_T = save_T;
        g = alloc_grid();
    }

    /* ==================== SECTION 7: Summary ==================== */
    printf("\n=== SUMMARY: Sorted by combined score ===\n");
    printf("  (score = 2*trans_l2 + 1.5*min(torsion,1) + fc + 0.5*|P|_norm)\n\n");

    /* Compute combined score */
    double scores[64];
    for (int i = 0; i < nr; i++) {
        if (!results[i].stable) { scores[i] = -1; continue; }
        scores[i] = 2.0 * results[i].transverse_l2
                  + 1.5 * fmin(results[i].torsion_flux / 1.0, 1.0)
                  + results[i].fc
                  + 0.5 * fmin(results[i].peak_P / 2.0, 1.0);
    }

    /* Sort by score (simple insertion sort) */
    int order[64];
    for (int i = 0; i < nr; i++) order[i] = i;
    for (int i = 1; i < nr; i++) {
        int j = i;
        while (j > 0 && scores[order[j-1]] < scores[order[j]]) {
            int tmp = order[j]; order[j] = order[j-1]; order[j-1] = tmp;
            j--;
        }
    }

    printf("  Rank  Score   trans_l2  torsion  fc      |P|     wind    Energy        Label\n");
    printf("  ----  ------  --------  -------  ------  ------  ------  ----------    -----\n");
    for (int i = 0; i < nr; i++) {
        int idx = order[i];
        Result *r = &results[idx];
        printf("  %3d   %5.3f   %.4f    %.4f   %.4f  %.4f  %+.3f  %10.1f    %s\n",
               i+1, scores[idx], r->transverse_l2, r->torsion_flux,
               r->fc, r->peak_P, r->winding, r->energy, r->label);
    }

    /* Check if any combined beats both pure A and B */
    printf("\n=== KEY QUESTION: Does any combination beat BOTH controls? ===\n");
    double A_trans = results[0].transverse_l2;
    double B_torsion = results[1].torsion_flux;
    printf("  Pure A trans_l2 = %.4f, Pure B torsion = %.4f\n", A_trans, B_torsion);

    int found_better = 0;
    for (int i = 2; i < nr; i++) {
        if (!results[i].stable) continue;
        if (results[i].transverse_l2 > 0.8 * A_trans &&
            results[i].torsion_flux > 0.8 * B_torsion) {
            printf("  ** %s: trans=%.4f (%.0f%% of A), tor=%.4f (%.0f%% of B) **\n",
                   results[i].label, results[i].transverse_l2,
                   100*results[i].transverse_l2/A_trans,
                   results[i].torsion_flux,
                   100*results[i].torsion_flux/B_torsion);
            found_better = 1;
        }
    }
    if (!found_better) {
        printf("  No combination achieved >80%% of both A's trans_l2 AND B's torsion.\n");
        printf("  Closest:\n");
        /* Find best combined (trans + torsion) */
        double best_combined = -1;
        int best_idx = -1;
        for (int i = 2; i < nr; i++) {
            if (!results[i].stable) continue;
            double comb = results[i].transverse_l2 / (A_trans + 1e-30)
                        + results[i].torsion_flux / (B_torsion + 1e-30);
            if (comb > best_combined) { best_combined = comb; best_idx = i; }
        }
        if (best_idx >= 0) {
            printf("    %s: trans=%.4f (%.0f%%), tor=%.4f (%.0f%%)\n",
                   results[best_idx].label,
                   results[best_idx].transverse_l2,
                   100*results[best_idx].transverse_l2/A_trans,
                   results[best_idx].torsion_flux,
                   100*results[best_idx].torsion_flux/B_torsion);
        }
    }

    free_grid(g);
    return 0;
}
