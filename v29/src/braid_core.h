/*  braid_core.h — Shared PDE solver for V29 bimodal braid tests
 *
 *  Include this header in each test's .c file. It provides:
 *  - Grid allocation and setup
 *  - Bimodal sweet spot parameters
 *  - Braid initialization (with configurable m_init)
 *  - Force computation (OMP-parallel)
 *  - Velocity Verlet time stepping
 *  - Absorbing damping (optional, caller decides)
 *  - Full diagnostics (fc, l2, torsion, winding, energy, |P|)
 *
 *  Build: gcc -O3 -fopenmp -o test test.c -lm
 */

#ifndef BRAID_CORE_H
#define BRAID_CORE_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>

#define NFIELDS 3
#define NDIM    16
#define PI      3.14159265358979323846

/* ================================================================
   Bimodal sweet spot parameters
   ================================================================ */

/* Parameter order: A1, A2, A3, delta2, delta3, R_tube, ellip, ell_ang,
                    k_fac, A_bg, R_disp, ell_rot, mu, kappa, mass, lam_pw */

static const double PATH_A[NDIM] = {
    0.8, 0.8, 0.8,   0.00, 1.67,
    3.0, 0.80, 0.0,  1.0, 0.0,
    0.0, 0.0,        -29.7, 50.0, 1.50, 0.0
};

static const double PATH_B[NDIM] = {
    0.8, 0.8, 0.8,   3.53, 4.92,
    3.0, 0.25, 0.0,  1.0, 0.0,
    0.0, 0.0,        -43.4, 50.0, 1.50, 0.0
};

/* t=0.85 interpolation (the validated sweet spot) */
static double BIMODAL[NDIM];

static void bimodal_init_params(void) {
    for (int d = 0; d < NDIM; d++)
        BIMODAL[d] = 0.15 * PATH_A[d] + 0.85 * PATH_B[d];
}

/* Convenience: interpolate at arbitrary t */
static void interp_params(double *out, double t) {
    for (int d = 0; d < NDIM; d++)
        out[d] = (1.0 - t) * PATH_A[d] + t * PATH_B[d];
}

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
   Grid management
   ================================================================ */

static Grid *grid_alloc(int N, double L) {
    Grid *g = calloc(1, sizeof(Grid));
    int N3 = N * N * N;
    for (int a = 0; a < NFIELDS; a++) {
        g->phi[a] = calloc(N3, sizeof(double));
        g->vel[a] = calloc(N3, sizeof(double));
        g->acc[a] = calloc(N3, sizeof(double));
    }
    g->N = N; g->L = L;
    g->dx = 2.0 * L / (N - 1);
    g->dt = 0.20 * g->dx;
    return g;
}

static void grid_free(Grid *g) {
    for (int a = 0; a < NFIELDS; a++) {
        free(g->phi[a]); free(g->vel[a]); free(g->acc[a]);
    }
    free(g);
}

static void grid_zero(Grid *g) {
    int N3 = g->N * g->N * g->N;
    for (int a = 0; a < NFIELDS; a++) {
        memset(g->phi[a], 0, N3 * sizeof(double));
        memset(g->vel[a], 0, N3 * sizeof(double));
        memset(g->acc[a], 0, N3 * sizeof(double));
    }
}

/* ================================================================
   Braid initialization
   m_init_override: if >= 0, use this instead of phys[14] for omega.
                    Set to -1 to use phys[14] (default behavior).
   ================================================================ */

static void init_braid(Grid *g, const double *phys, double m_init_override) {
    double A[3]     = {phys[0], phys[1], phys[2]};
    double delta[3] = {0.0, phys[3], phys[4]};
    double R_tube   = phys[5];
    double ellip    = phys[6];
    double ell_ang  = phys[7];
    double k_fac    = phys[8];
    double A_bg     = phys[9];
    double R_disp   = phys[10];
    double ell_rot  = phys[11];
    double m_init   = (m_init_override >= 0) ? m_init_override : phys[14];

    int N = g->N, NN = N * N;
    double dx = g->dx, L = g->L;
    double k = k_fac * PI / L;
    double omega = sqrt(k * k + m_init * m_init);
    double inv_2R2 = 1.0 / (2.0 * R_tube * R_tube);
    double sx = 1.0 + ellip, sy = 1.0 - ellip;

    double cx[3], cy[3], ea[3];
    for (int a = 0; a < 3; a++) {
        double ang = 2.0 * PI * a / 3.0;
        cx[a] = R_disp * cos(ang);
        cy[a] = R_disp * sin(ang);
        ea[a] = (ell_rot > 0.5) ? ell_ang + 2.0*PI*a/3.0 : ell_ang;
    }

    grid_zero(g);

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
                    double xr = xc*ca + yc*sa;
                    double yr = -xc*sa + yc*ca;
                    double r2e = xr*xr/(sx*sx) + yr*yr/(sy*sy);
                    double env = exp(-r2e * inv_2R2);
                    double ph = k * z + delta[a];
                    double amp = A[a] * env + A_bg;
                    g->phi[a][idx] = amp * cos(ph);
                    g->vel[a][idx] = omega * amp * sin(ph);
                }
            }
        }
    }
}

/* ================================================================
   Force computation
   mass2_override: if >= 0, use this as m² in the EOM (NOT init).
                   Set to -1 to use phys[14]².
   ================================================================ */

static void compute_forces(Grid *g, const double *phys, double mass2_override) {
    int N = g->N, NN = N*N, N3 = N*N*N;
    double idx2 = 1.0 / (g->dx * g->dx);
    double mu    = phys[12];
    double kappa = phys[13];
    double mass2 = (mass2_override >= 0) ? mass2_override : phys[14]*phys[14];
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
        double P  = p0 * p1 * p2;
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
   Time stepping (Velocity Verlet)
   ================================================================ */

static void verlet_kick(Grid *g, double half_dt) {
    int N3 = g->N * g->N * g->N;
    for (int a = 0; a < NFIELDS; a++) {
        double *v = g->vel[a], *ac = g->acc[a];
        for (int idx = 0; idx < N3; idx++)
            v[idx] += half_dt * ac[idx];
    }
}

static void verlet_drift(Grid *g) {
    int N3 = g->N * g->N * g->N;
    double dt = g->dt;
    for (int a = 0; a < NFIELDS; a++) {
        double *p = g->phi[a], *v = g->vel[a];
        for (int idx = 0; idx < N3; idx++)
            p[idx] += dt * v[idx];
    }
}

/* Full Verlet step: kick-drift-force-kick */
static void verlet_full_step(Grid *g, const double *phys, double mass2_ov) {
    double hdt = 0.5 * g->dt;
    verlet_kick(g, hdt);
    verlet_drift(g);
    compute_forces(g, phys, mass2_ov);
    verlet_kick(g, hdt);
}

/* ================================================================
   Absorbing damping layer (xy only, cylindrical)
   Call this AFTER verlet_full_step if you want absorbing BC.
   ================================================================ */

static void apply_damping_xy(Grid *g) {
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

static int check_blowup(Grid *g) {
    int N3 = g->N * g->N * g->N;
    for (int a = 0; a < NFIELDS; a++)
        for (int idx = 0; idx < N3; idx += 37)
            if (fabs(g->phi[a][idx]) > 50.0) return 1;
    return 0;
}

static double compute_winding(Grid *g) {
    int N = g->N, NN = N*N, ic = N/2, jc = N/2;
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
                if (fabs(P) > peak_P) peak_P = fabs(P);

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

    double Qxx = 3*Ixx-R2sum, Qyy = 3*Iyy-R2sum, Qzz = 3*Izz-R2sum, Qxy = 3*Ixy;
    double Q2 = Qxx*Qxx + Qyy*Qyy + Qzz*Qzz + 2*Qxy*Qxy;
    double R2m = R2sum / (M0 + 1e-30);
    res->l2_frac = sqrt(Q2) / (M0 * R2m + 1e-30);
    res->transverse_l2 = fabs(Ixx - Iyy) / (Ixx + Iyy + 1e-30);
    res->fc = phi2_core / (phi2_total + 1e-30);
    res->peak_P = peak_P;
    res->energy = E_total;

    /* Torsion flux — all 3 components */
    double Phi_xy = 0, Phi_xz = 0, Phi_yz = 0;
    int k_mid = N/2, j_mid = N/2, i_mid = N/2;
    for (int i = 2; i < N-2; i++) {
        for (int j = 2; j < N-2; j++) {
            {
                int idx = i*NN + j*N + k_mid;
                double dxp1 = (g->phi[1][idx+NN]-g->phi[1][idx-NN])/(2*dx);
                double dyp0 = (g->phi[0][idx+N] -g->phi[0][idx-N]) /(2*dx);
                Phi_xy += 0.5*(dxp1 - dyp0)*dx*dx;
            }
            {
                int idx = i*NN + j_mid*N + j;
                int kp=(j+1)%N, km=(j-1+N)%N;
                double dxp2 = (g->phi[2][idx+NN]-g->phi[2][idx-NN])/(2*dx);
                double dzp0 = (g->phi[0][i*NN+j_mid*N+kp]-g->phi[0][i*NN+j_mid*N+km])/(2*dx);
                Phi_xz += 0.5*(dxp2 - dzp0)*dx*dx;
            }
            {
                int idx = i_mid*NN + i*N + j;
                int kp=(j+1)%N, km=(j-1+N)%N;
                double dyp2 = (g->phi[2][idx+N] -g->phi[2][idx-N]) /(2*dx);
                double dzp1 = (g->phi[1][i_mid*NN+i*N+kp]-g->phi[1][i_mid*NN+i*N+km])/(2*dx);
                Phi_yz += 0.5*(dyp2 - dzp1)*dx*dx;
            }
        }
    }
    res->torsion_flux = fmax(fabs(Phi_xy), fmax(fabs(Phi_xz), fabs(Phi_yz)));
    res->winding = compute_winding(g);
}

static void print_result(const Result *r) {
    printf("  %-28s  E=%10.1f  fc=%.4f  l2=%.4f  trans=%.4f  tor=%.4f  |P|=%.4f  w=%+.3f  %s\n",
           r->label, r->energy, r->fc, r->l2_frac, r->transverse_l2,
           r->torsion_flux, r->peak_P, r->winding,
           r->stable ? "OK" : "BLOWUP");
}

#endif /* BRAID_CORE_H */
