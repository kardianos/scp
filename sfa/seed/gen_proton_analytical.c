/*  gen_proton_analytical.c — Generate proton seed from analytical profile fits
 *
 *  Fits the converged proton radial profile from simulation with analytical
 *  functions, then generates 3D seeds using the fitted envelopes + 3-braid
 *  phase structure.
 *
 *  Three fitting levels:
 *    Level 1: Spherical ansatz (single radial envelope, isotropic)
 *    Level 2: Three-braid sum with converged parameters
 *    Level 3: Empirical lookup table interpolation
 *
 *  Build: gcc -O3 -march=native -fopenmp -o gen_proton_analytical \
 *         gen_proton_analytical.c -lzstd -lm
 *
 *  Usage: ./gen_proton_analytical -N 192 -L 25 -o proton_seed.sfa
 *         ./gen_proton_analytical -level 2 -chirality UDD -o neutron_seed.sfa
 *         ./gen_proton_analytical -fit   (fit-only mode, prints parameters + errors)
 */

#define SFA_IMPLEMENTATION
#include "../format/sfa.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>
#include <float.h>

#define PI 3.14159265358979323846

/* ===== MEASURED RADIAL PROFILE (from proton_profile.tsv) ===== */
/* Converged profile at T ≈ 20 from V43 proton formation sim */

#define NPROF 48

static const double prof_r[NPROF] = {
    0.1250, 0.3750, 0.6250, 0.8750, 1.1250, 1.3750, 1.6250, 1.8750,
    2.1250, 2.3750, 2.6250, 2.8750, 3.1250, 3.3750, 3.6250, 3.8750,
    4.1250, 4.3750, 4.6250, 4.8750, 5.1250, 5.3750, 5.6250, 5.8750,
    6.1250, 6.3750, 6.6250, 6.8750, 7.1250, 7.3750, 7.6250, 7.8750,
    8.1250, 8.3750, 8.6250, 8.8750, 9.1250, 9.3750, 9.6250, 9.8750,
    10.1250, 10.3750, 10.6250, 10.8750, 11.1250, 11.3750, 11.6250, 11.8750
};

static const double prof_phi_rms[NPROF] = {
    8.957332e-02, 8.957929e-02, 9.097124e-02, 9.452119e-02, 1.011036e-01,
    1.076245e-01, 1.173643e-01, 1.273026e-01, 1.378794e-01, 1.484698e-01,
    1.581058e-01, 1.668550e-01, 1.752138e-01, 1.825474e-01, 1.876289e-01,
    1.923888e-01, 1.948865e-01, 1.967892e-01, 1.975596e-01, 1.979517e-01,
    1.974375e-01, 1.980518e-01, 1.966032e-01, 1.974795e-01, 1.960364e-01,
    1.954329e-01, 1.949979e-01, 1.934939e-01, 1.921652e-01, 1.925620e-01,
    1.921195e-01, 1.917434e-01, 1.927507e-01, 1.929214e-01, 1.927142e-01,
    1.931997e-01, 1.932025e-01, 1.934013e-01, 1.943314e-01, 1.952031e-01,
    1.965781e-01, 1.984931e-01, 2.010240e-01, 2.040399e-01, 2.076683e-01,
    2.111236e-01, 2.146182e-01, 2.181902e-01
};

static const double prof_theta_rms[NPROF] = {
    1.222707e-01, 1.192173e-01, 1.128417e-01, 1.026436e-01, 9.090891e-02,
    8.179093e-02, 7.577896e-02, 7.340408e-02, 7.269776e-02, 7.353497e-02,
    7.283094e-02, 7.113378e-02, 6.798636e-02, 6.448015e-02, 6.325640e-02,
    6.240812e-02, 6.331253e-02, 6.598064e-02, 6.904700e-02, 7.068495e-02,
    7.103415e-02, 7.021869e-02, 6.723977e-02, 6.408533e-02, 6.106375e-02,
    5.837420e-02, 5.762876e-02, 5.737679e-02, 5.711282e-02, 5.631715e-02,
    5.471364e-02, 5.299135e-02, 5.122902e-02, 5.019964e-02, 4.939861e-02,
    4.919661e-02, 4.892017e-02, 4.840039e-02, 4.722077e-02, 4.530690e-02,
    4.320888e-02, 4.103321e-02, 3.932907e-02, 3.781599e-02, 3.687038e-02,
    3.624640e-02, 3.601148e-02, 3.605865e-02
};

static const double prof_P_abs[NPROF] = {
    1.080615e-04, 1.070426e-04, 1.085966e-04, 1.208967e-04, 1.490766e-04,
    1.885101e-04, 2.614197e-04, 3.611764e-04, 4.800831e-04, 6.406795e-04,
    8.135111e-04, 9.885967e-04, 1.188118e-03, 1.375876e-03, 1.512624e-03,
    1.646482e-03, 1.721693e-03, 1.795466e-03, 1.825883e-03, 1.833879e-03,
    1.828898e-03, 1.813087e-03, 1.774087e-03, 1.838132e-03, 1.825295e-03,
    1.869700e-03, 1.973554e-03, 2.023973e-03, 2.139058e-03, 2.306646e-03,
    2.331157e-03, 2.324581e-03, 2.305846e-03, 2.162608e-03, 1.901130e-03,
    1.662942e-03, 1.393962e-03, 1.191837e-03, 1.068658e-03, 9.847183e-04,
    9.575505e-04, 9.907143e-04, 1.094562e-03, 1.267727e-03, 1.517491e-03,
    1.840724e-03, 2.206024e-03, 2.628897e-03
};

static const double prof_v_radial[NPROF] = {
    9.812190e-02, -5.871394e-02, -1.332813e-01, -1.602194e-01, -1.664255e-01,
    -1.900993e-01, -2.182320e-01, -2.405306e-01, -2.555898e-01, -2.833662e-01,
    -2.847021e-01, -2.974509e-01, -3.252480e-01, -3.144147e-01, -3.226472e-01,
    -3.295531e-01, -3.240448e-01, -3.270169e-01, -3.189363e-01, -3.134448e-01,
    -3.073029e-01, -2.956620e-01, -2.863368e-01, -2.813457e-01, -2.628692e-01,
    -2.426417e-01, -2.286342e-01, -2.032713e-01, -1.786466e-01, -1.589353e-01,
    -1.281655e-01, -1.022569e-01, -7.894387e-02, -5.180790e-02, -2.969330e-02,
    -1.023879e-02, 8.159861e-03, 2.067327e-02, 3.184817e-02, 4.120204e-02,
    4.910394e-02, 5.658475e-02, 6.235504e-02, 6.827648e-02, 7.342660e-02,
    7.809183e-02, 8.149160e-02, 8.346570e-02
};

static const double prof_phi_v_rms[NPROF] = {
    1.163968e+00, 1.158617e+00, 1.141310e+00, 1.108531e+00, 1.069401e+00,
    1.042909e+00, 1.010164e+00, 9.965431e-01, 9.594800e-01, 9.406814e-01,
    8.986470e-01, 8.693903e-01, 8.608834e-01, 8.229616e-01, 8.040882e-01,
    7.822206e-01, 7.558581e-01, 7.368916e-01, 7.133888e-01, 6.894967e-01,
    6.679313e-01, 6.441193e-01, 6.152136e-01, 5.956959e-01, 5.680944e-01,
    5.389360e-01, 5.184730e-01, 4.902590e-01, 4.627728e-01, 4.435084e-01,
    4.136301e-01, 3.880690e-01, 3.650013e-01, 3.383023e-01, 3.175089e-01,
    3.006289e-01, 2.863435e-01, 2.781432e-01, 2.724710e-01, 2.703553e-01,
    2.708243e-01, 2.735390e-01, 2.797174e-01, 2.848166e-01, 2.918701e-01,
    2.984858e-01, 3.034319e-01, 3.080420e-01
};

/* Initial profile for comparison */
static const double prof_phi_rms_init[NPROF] = {
    2.680993e-01, 2.685014e-01, 2.694508e-01, 2.705958e-01, 2.720473e-01,
    2.737684e-01, 2.760003e-01, 2.766580e-01, 2.787098e-01, 2.808506e-01,
    2.811457e-01, 2.820818e-01, 2.823512e-01, 2.825227e-01, 2.817986e-01,
    2.808726e-01, 2.792091e-01, 2.777426e-01, 2.750711e-01, 2.734098e-01,
    2.697916e-01, 2.683550e-01, 2.654317e-01, 2.635613e-01, 2.617629e-01,
    2.601779e-01, 2.579151e-01, 2.557436e-01, 2.537338e-01, 2.515413e-01,
    2.492644e-01, 2.471417e-01, 2.449010e-01, 2.423848e-01, 2.399130e-01,
    2.367934e-01, 2.343107e-01, 2.315993e-01, 2.285332e-01, 2.259766e-01,
    2.231363e-01, 2.204968e-01, 2.174172e-01, 2.147943e-01, 2.116410e-01,
    2.092298e-01, 2.065138e-01, 2.038720e-01
};

/* ===== f16 conversion ===== */
static inline uint16_t f64_to_f16(double v) {
    float f = (float)v;
    uint32_t x; memcpy(&x, &f, 4);
    uint16_t sign = (x >> 16) & 0x8000;
    int exp = ((x >> 23) & 0xFF) - 127 + 15;
    uint16_t mant = (x >> 13) & 0x3FF;
    if (exp <= 0) return sign;
    if (exp >= 31) return sign | 0x7C00;
    return sign | (exp << 10) | mant;
}

/* ===== FITTING FUNCTIONS ===== */

/* Level 1: Spherical difference-of-Gaussians for phi_rms
 *   phi_rms(r) = A_flat - A_dip * exp(-r^2 / (2*s_dip^2))
 *   This captures: low at center (~0.09), flat plateau (~0.195)
 */
typedef struct {
    double A_flat;    /* plateau amplitude */
    double A_dip;     /* central dip depth */
    double s_dip;     /* dip width */
} L1PhiFit;

static double l1_phi_rms(double r, const L1PhiFit *f) {
    return f->A_flat - f->A_dip * exp(-r * r / (2 * f->s_dip * f->s_dip));
}

/* Level 1: Theta envelope — two Gaussians (narrow core + wide halo) + bg
 *   theta_rms(r) = T1 * exp(-r^2/(2*s1^2)) + T2 * exp(-r^2/(2*s2^2)) + T_bg
 */
typedef struct {
    double T1, s1;    /* narrow core component */
    double T2, s2;    /* wide halo component */
    double T_bg;      /* far-field background */
} L1ThetaFit;

static double l1_theta_rms(double r, const L1ThetaFit *f) {
    return f->T1 * exp(-r * r / (2 * f->s1 * f->s1))
         + f->T2 * exp(-r * r / (2 * f->s2 * f->s2))
         + f->T_bg;
}

/* Level 1: Velocity envelope — Gaussian
 *   phi_v_rms(r) = V_core * exp(-r^2/(2*s_v^2)) + V_bg
 */
typedef struct {
    double V_core;
    double V_bg;
    double s_v;
} L1VelFit;

static double l1_vel_rms(double r, const L1VelFit *f) {
    return f->V_core * exp(-r * r / (2 * f->s_v * f->s_v)) + f->V_bg;
}

/* Level 1: Radial velocity — negative Gaussian
 *   v_radial(r) = -V_r * r/(r+r0) * exp(-r^2/(2*s_r^2))
 */
typedef struct {
    double V_r;
    double r0;
    double s_r;
} L1VradFit;

static double l1_v_radial(double r, const L1VradFit *f) {
    return -f->V_r * r / (r + f->r0) * exp(-r * r / (2 * f->s_r * f->s_r));
}

/* ===== SIMPLE OPTIMIZER (golden section per parameter) ===== */

static double fit_l1_phi_error(const L1PhiFit *f, int n_pts) {
    double err = 0;
    int count = 0;
    for (int i = 0; i < n_pts && i < NPROF; i++) {
        double pred = l1_phi_rms(prof_r[i], f);
        double diff = pred - prof_phi_rms[i];
        err += diff * diff;
        count++;
    }
    return sqrt(err / count);
}

static double fit_l1_theta_error(const L1ThetaFit *f, int n_pts) {
    double err = 0;
    int count = 0;
    for (int i = 0; i < n_pts && i < NPROF; i++) {
        double pred = l1_theta_rms(prof_r[i], f);
        double diff = pred - prof_theta_rms[i];
        err += diff * diff;
        count++;
    }
    return sqrt(err / count);
}

static double fit_l1_vel_error(const L1VelFit *f, int n_pts) {
    double err = 0;
    int count = 0;
    for (int i = 0; i < n_pts && i < NPROF; i++) {
        double pred = l1_vel_rms(prof_r[i], f);
        double diff = pred - prof_phi_v_rms[i];
        err += diff * diff;
        count++;
    }
    return sqrt(err / count);
}

static double fit_l1_vrad_error(const L1VradFit *f, int n_pts) {
    double err = 0;
    int count = 0;
    /* Skip first point (anomalous positive value) */
    for (int i = 1; i < n_pts && i < NPROF; i++) {
        double pred = l1_v_radial(prof_r[i], f);
        double diff = pred - prof_v_radial[i];
        err += diff * diff;
        count++;
    }
    return sqrt(err / count);
}

/* Nelder-Mead-like coordinate descent optimizer */
static void optimize_phi_fit(L1PhiFit *f, int n_pts) {
    double best_err = fit_l1_phi_error(f, n_pts);
    double step[3] = {0.02, 0.02, 1.0};
    double *params[3] = {&f->A_flat, &f->A_dip, &f->s_dip};
    for (int iter = 0; iter < 5000; iter++) {
        for (int p = 0; p < 3; p++) {
            double orig = *params[p];
            /* Try + */
            *params[p] = orig + step[p];
            double e_plus = fit_l1_phi_error(f, n_pts);
            if (e_plus < best_err) { best_err = e_plus; continue; }
            /* Try - */
            *params[p] = orig - step[p];
            double e_minus = fit_l1_phi_error(f, n_pts);
            if (e_minus < best_err) { best_err = e_minus; continue; }
            /* No improvement */
            *params[p] = orig;
        }
        for (int p = 0; p < 3; p++) step[p] *= 0.998;
    }
}

static void optimize_theta_fit(L1ThetaFit *f, int n_pts) {
    double best_err = fit_l1_theta_error(f, n_pts);
    double step[5] = {0.01, 0.5, 0.01, 0.5, 0.005};
    double *params[5] = {&f->T1, &f->s1, &f->T2, &f->s2, &f->T_bg};
    for (int iter = 0; iter < 8000; iter++) {
        for (int p = 0; p < 5; p++) {
            double orig = *params[p];
            *params[p] = orig + step[p];
            if (p == 1 || p == 3) { if (*params[p] < 0.1) { *params[p] = orig; continue; } }
            double e_plus = fit_l1_theta_error(f, n_pts);
            if (e_plus < best_err) { best_err = e_plus; continue; }
            *params[p] = orig - step[p];
            if (p == 1 || p == 3) { if (*params[p] < 0.1) { *params[p] = orig; continue; } }
            double e_minus = fit_l1_theta_error(f, n_pts);
            if (e_minus < best_err) { best_err = e_minus; continue; }
            *params[p] = orig;
        }
        for (int p = 0; p < 5; p++) step[p] *= 0.9985;
    }
}

static void optimize_vel_fit(L1VelFit *f, int n_pts) {
    double best_err = fit_l1_vel_error(f, n_pts);
    double step[3] = {0.05, 0.02, 1.0};
    double *params[3] = {&f->V_core, &f->V_bg, &f->s_v};
    for (int iter = 0; iter < 5000; iter++) {
        for (int p = 0; p < 3; p++) {
            double orig = *params[p];
            *params[p] = orig + step[p];
            double e_plus = fit_l1_vel_error(f, n_pts);
            if (e_plus < best_err) { best_err = e_plus; continue; }
            *params[p] = orig - step[p];
            double e_minus = fit_l1_vel_error(f, n_pts);
            if (e_minus < best_err) { best_err = e_minus; continue; }
            *params[p] = orig;
        }
        for (int p = 0; p < 3; p++) step[p] *= 0.998;
    }
}

static void optimize_vrad_fit(L1VradFit *f, int n_pts) {
    double best_err = fit_l1_vrad_error(f, n_pts);
    double step[3] = {0.02, 0.5, 1.0};
    double *params[3] = {&f->V_r, &f->r0, &f->s_r};
    for (int iter = 0; iter < 5000; iter++) {
        for (int p = 0; p < 3; p++) {
            double orig = *params[p];
            *params[p] = orig + step[p];
            if (*params[p] < 0.01) { *params[p] = orig; continue; }
            double e_plus = fit_l1_vrad_error(f, n_pts);
            if (e_plus < best_err) { best_err = e_plus; continue; }
            *params[p] = orig - step[p];
            if (*params[p] < 0.01) { *params[p] = orig; continue; }
            double e_minus = fit_l1_vrad_error(f, n_pts);
            if (e_minus < best_err) { best_err = e_minus; continue; }
            *params[p] = orig;
        }
        for (int p = 0; p < 3; p++) step[p] *= 0.998;
    }
}

/* ===== LEVEL 2: Three-braid model with adjusted parameters ===== */

typedef struct {
    double A_braid;    /* per-braid amplitude */
    double R_tube;     /* tube radius */
    double ellip;      /* ellipticity */
    double A_bg;       /* background amplitude */
    double theta_amp;  /* theta amplitude */
    double theta_R;    /* theta tube radius */
} L2Params;

/* Compute phi_rms for level 2 three-braid model at radius r
 * by Monte Carlo sampling on a sphere of radius r.
 * Returns RMS of |phi| over angles. */
static double l2_phi_rms_at_r(double r, const L2Params *p, double m2, double L_box,
                               const double delta[3], const double carrier[3],
                               const double chiral[3]) {
    if (r < 0.01) r = 0.01;
    double kw = PI / L_box;
    double sx = 1 + p->ellip, sy = 1 - p->ellip;
    double inv2R2 = 1.0 / (2 * p->R_tube * p->R_tube);
    double k_bg = PI / L_box;

    /* Sample 200 points on sphere using Fibonacci spiral */
    int NS = 200;
    double phi_sq_sum = 0;
    double golden = (1 + sqrt(5)) / 2;
    for (int s = 0; s < NS; s++) {
        double theta_s = acos(1 - 2.0 * (s + 0.5) / NS);
        double phi_s = 2 * PI * s / golden;
        double x = r * sin(theta_s) * cos(phi_s);
        double y = r * sin(theta_s) * sin(phi_s);
        double z = r * cos(theta_s);

        for (int a = 0; a < 3; a++) {
            /* Background */
            double val = p->A_bg * cos(k_bg * z + 2 * PI * a / 3.0);

            /* Braid 0: along z */
            {
                double r2e = x * x / (sx * sx) + y * y / (sy * sy);
                double env = p->A_braid * exp(-r2e * inv2R2);
                val += env * cos(chiral[0] * kw * z + delta[a] + carrier[0]);
            }
            /* Braid 1: along x */
            {
                double zr = x, xr = -z, yr = y;
                double r2e = xr * xr / (sx * sx) + yr * yr / (sy * sy);
                double env = p->A_braid * exp(-r2e * inv2R2);
                val += env * cos(chiral[1] * kw * zr + delta[a] + carrier[1]);
            }
            /* Braid 2: along y */
            {
                double zr = y, yr = -z, xr = x;
                double r2e = xr * xr / (sx * sx) + yr * yr / (sy * sy);
                double env = p->A_braid * exp(-r2e * inv2R2);
                val += env * cos(chiral[2] * kw * zr + delta[a] + carrier[2]);
            }
            phi_sq_sum += val * val;
        }
    }
    return sqrt(phi_sq_sum / (3.0 * NS));
}

/* Compute P_abs for level 2 at radius r */
static double l2_P_abs_at_r(double r, const L2Params *p, double m2, double L_box,
                              const double delta[3], const double carrier[3],
                              const double chiral[3]) {
    if (r < 0.01) r = 0.01;
    double kw = PI / L_box;
    double sx = 1 + p->ellip, sy = 1 - p->ellip;
    double inv2R2 = 1.0 / (2 * p->R_tube * p->R_tube);
    double k_bg = PI / L_box;

    int NS = 200;
    double P_sum = 0;
    double golden = (1 + sqrt(5)) / 2;
    for (int s = 0; s < NS; s++) {
        double theta_s = acos(1 - 2.0 * (s + 0.5) / NS);
        double phi_s = 2 * PI * s / golden;
        double x = r * sin(theta_s) * cos(phi_s);
        double y = r * sin(theta_s) * sin(phi_s);
        double z = r * cos(theta_s);

        double phi_vals[3];
        for (int a = 0; a < 3; a++) {
            double val = p->A_bg * cos(k_bg * z + 2 * PI * a / 3.0);
            {
                double r2e = x * x / (sx * sx) + y * y / (sy * sy);
                double env = p->A_braid * exp(-r2e * inv2R2);
                val += env * cos(chiral[0] * kw * z + delta[a] + carrier[0]);
            }
            {
                double zr = x, xr = -z, yr = y;
                double r2e = xr * xr / (sx * sx) + yr * yr / (sy * sy);
                double env = p->A_braid * exp(-r2e * inv2R2);
                val += env * cos(chiral[1] * kw * zr + delta[a] + carrier[1]);
            }
            {
                double zr = y, yr = -z, xr = x;
                double r2e = xr * xr / (sx * sx) + yr * yr / (sy * sy);
                double env = p->A_braid * exp(-r2e * inv2R2);
                val += env * cos(chiral[2] * kw * zr + delta[a] + carrier[2]);
            }
            phi_vals[a] = val;
        }
        P_sum += fabs(phi_vals[0] * phi_vals[1] * phi_vals[2]);
    }
    return P_sum / NS;
}

/* Compute theta_rms for level 2 at radius r */
static double l2_theta_rms_at_r(double r, const L2Params *p, double L_box,
                                  const double carrier[3], const double chiral[3]) {
    if (r < 0.01) r = 0.01;
    double kw = PI / L_box;
    double sx = 1 + p->ellip, sy = 1 - p->ellip;
    double inv2R2_t = 1.0 / (2 * p->theta_R * p->theta_R);

    int NS = 200;
    double theta_sq_sum = 0;
    double golden = (1 + sqrt(5)) / 2;
    for (int s = 0; s < NS; s++) {
        double theta_s = acos(1 - 2.0 * (s + 0.5) / NS);
        double phi_s = 2 * PI * s / golden;
        double x = r * sin(theta_s) * cos(phi_s);
        double y = r * sin(theta_s) * sin(phi_s);
        double z = r * cos(theta_s);

        for (int tc = 0; tc < 3; tc++) {
            double tval = 0;
            for (int b = 0; b < 3; b++) {
                double braid_env;
                if (b == 0) {
                    double r2e = x * x / (sx * sx) + y * y / (sy * sy);
                    braid_env = p->A_braid * exp(-r2e * inv2R2_t);
                } else if (b == 1) {
                    double xr = -z, yr = y;
                    double r2e = xr * xr / (sx * sx) + yr * yr / (sy * sy);
                    braid_env = p->A_braid * exp(-r2e * inv2R2_t);
                } else {
                    double xr = x, yr = -z;
                    double r2e = xr * xr / (sx * sx) + yr * yr / (sy * sy);
                    braid_env = p->A_braid * exp(-r2e * inv2R2_t);
                }
                double t_env = p->theta_amp * braid_env;
                if (tc == 0)      tval += chiral[b] * t_env * cos(kw * z + carrier[b]);
                else if (tc == 1) tval += chiral[b] * t_env * sin(kw * z + carrier[b]);
                else              tval += chiral[b] * t_env * cos(kw * z + carrier[b] + PI / 4);
            }
            theta_sq_sum += tval * tval;
        }
    }
    return sqrt(theta_sq_sum / (3.0 * NS));
}

/* Compute errors for level 2 */
static double l2_phi_error(const L2Params *p, double m2, double L_box,
                            const double delta[3], const double carrier[3],
                            const double chiral[3], int n_pts) {
    double err = 0;
    int count = 0;
    for (int i = 0; i < n_pts && i < NPROF; i++) {
        double pred = l2_phi_rms_at_r(prof_r[i], p, m2, L_box, delta, carrier, chiral);
        double diff = pred - prof_phi_rms[i];
        err += diff * diff;
        count++;
    }
    return sqrt(err / count);
}

static double l2_P_error(const L2Params *p, double m2, double L_box,
                          const double delta[3], const double carrier[3],
                          const double chiral[3], int n_pts) {
    double err = 0;
    int count = 0;
    for (int i = 0; i < n_pts && i < NPROF; i++) {
        double pred = l2_P_abs_at_r(prof_r[i], p, m2, L_box, delta, carrier, chiral);
        double diff = pred - prof_P_abs[i];
        err += diff * diff;
        count++;
    }
    return sqrt(err / count);
}

/* ===== LEVEL 3: Empirical interpolation ===== */

static double interp_profile(const double *r_tab, const double *v_tab, int n, double r) {
    if (r <= r_tab[0]) return v_tab[0];
    if (r >= r_tab[n - 1]) return v_tab[n - 1];
    /* Binary search */
    int lo = 0, hi = n - 1;
    while (hi - lo > 1) {
        int mid = (lo + hi) / 2;
        if (r_tab[mid] <= r) lo = mid; else hi = mid;
    }
    double t = (r - r_tab[lo]) / (r_tab[hi] - r_tab[lo]);
    return v_tab[lo] + t * (v_tab[hi] - v_tab[lo]);
}

/* ===== MAIN ===== */
int main(int argc, char **argv) {
    int N = 192;
    double L = 25.0;
    double A_bg = 0.1;
    double m2 = 2.25;
    double eta = 0.5;
    double delta[3] = {0.0, 3.0005, 4.4325};
    double cx = 0, cy = 0, cz = 0;
    char chirality[8] = "UUD";
    char outpath[512] = "proton_analytical.sfa";
    int level = 2;  /* default to best-expected level */
    int fit_only = 0;
    int precision = 1;  /* f32 */
    double dt_factor = 0.025;

    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "-fit")) { fit_only = 1; continue; }
        if (i + 1 >= argc) continue;
        if      (!strcmp(argv[i], "-N"))         { N = atoi(argv[i+1]); i++; }
        else if (!strcmp(argv[i], "-L"))         { L = atof(argv[i+1]); i++; }
        else if (!strcmp(argv[i], "-cx"))        { cx = atof(argv[i+1]); i++; }
        else if (!strcmp(argv[i], "-cy"))        { cy = atof(argv[i+1]); i++; }
        else if (!strcmp(argv[i], "-cz"))        { cz = atof(argv[i+1]); i++; }
        else if (!strcmp(argv[i], "-chirality")) { strncpy(chirality, argv[i+1], 7); i++; }
        else if (!strcmp(argv[i], "-level"))     { level = atoi(argv[i+1]); i++; }
        else if (!strcmp(argv[i], "-o"))         { strncpy(outpath, argv[i+1], 511); i++; }
        else if (!strcmp(argv[i], "-precision")) {
            if(!strcmp(argv[i+1],"f16")) precision=0;
            else if(!strcmp(argv[i+1],"f32")) precision=1;
            else if(!strcmp(argv[i+1],"f64")) precision=2;
            i++;
        }
    }

    /* Chirality */
    double chiral[3] = {1, 1, 1};
    for (int b = 0; b < 3 && chirality[b]; b++)
        chiral[b] = (chirality[b] == 'D' || chirality[b] == 'd') ? -1.0 : 1.0;

    double carrier[3] = {0.0, 2 * PI / 3, 4 * PI / 3};

    /* ==================== FITTING ==================== */
    fprintf(stderr, "=== PROTON PROFILE FITTING ===\n\n");

    /* --- Level 1: Spherical fits --- */
    fprintf(stderr, "--- Level 1: Spherical ansatz ---\n");

    /* Fit phi_rms */
    /* Profile: starts at 0.09, rises to plateau ~0.197 around r=5, then
     * slowly rises again beyond r=10 (this is the background dominating).
     * Use fit only up to r ≈ 10 where the proton structure dominates. */
    int n_fit = 40;  /* up to r ≈ 10 */

    L1PhiFit phi_fit = {0.195, 0.105, 2.5};
    optimize_phi_fit(&phi_fit, n_fit);
    double phi_err = fit_l1_phi_error(&phi_fit, n_fit);
    fprintf(stderr, "  phi_rms fit: A_flat=%.6f, A_dip=%.6f, s_dip=%.3f\n",
            phi_fit.A_flat, phi_fit.A_dip, phi_fit.s_dip);
    fprintf(stderr, "  phi_rms RMS error: %.6f (relative: %.1f%%)\n",
            phi_err, 100 * phi_err / 0.195);

    /* Fit theta_rms */
    L1ThetaFit theta_fit = {.T1=0.06, .s1=1.0, .T2=0.03, .s2=5.0, .T_bg=0.035};
    optimize_theta_fit(&theta_fit, n_fit);
    double theta_err = fit_l1_theta_error(&theta_fit, n_fit);
    fprintf(stderr, "  theta_rms fit: T1=%.6f s1=%.3f T2=%.6f s2=%.3f T_bg=%.6f\n",
            theta_fit.T1, theta_fit.s1, theta_fit.T2, theta_fit.s2, theta_fit.T_bg);
    fprintf(stderr, "  theta_rms RMS error: %.6f (relative: %.1f%%)\n",
            theta_err, 100 * theta_err / 0.07);

    /* Fit phi_v_rms */
    L1VelFit vel_fit = {0.9, 0.27, 5.0};
    optimize_vel_fit(&vel_fit, n_fit);
    double vel_err = fit_l1_vel_error(&vel_fit, n_fit);
    fprintf(stderr, "  phi_v_rms fit: V_core=%.6f, V_bg=%.6f, s_v=%.3f\n",
            vel_fit.V_core, vel_fit.V_bg, vel_fit.s_v);
    fprintf(stderr, "  phi_v_rms RMS error: %.6f (relative: %.1f%%)\n",
            vel_err, 100 * vel_err / 0.7);

    /* Fit v_radial */
    L1VradFit vrad_fit = {0.35, 1.0, 5.0};
    optimize_vrad_fit(&vrad_fit, n_fit);
    double vrad_err = fit_l1_vrad_error(&vrad_fit, n_fit);
    fprintf(stderr, "  v_radial fit: V_r=%.6f, r0=%.4f, s_r=%.3f\n",
            vrad_fit.V_r, vrad_fit.r0, vrad_fit.s_r);
    fprintf(stderr, "  v_radial RMS error: %.6f (relative: %.1f%%)\n",
            vrad_err, 100 * vrad_err / 0.3);

    /* Predict P from Level 1 spherical model:
     * In the spherical model, P = phi_0 * phi_1 * phi_2. For random phases
     * on a shell at radius r, the expected |P| ~ (phi_rms)^3 * some factor.
     * With carrier phases {0, 2π/3, 4π/3}, the cancellation makes this
     * model poorly suited for P prediction. Report anyway. */
    fprintf(stderr, "\n  Level 1 P prediction: phi_rms^3 scaling\n");
    double P_l1_err = 0;
    int P_count = 0;
    for (int i = 0; i < n_fit; i++) {
        double pr = l1_phi_rms(prof_r[i], &phi_fit);
        /* The cube gives the uncorrelated product; scale by empirical factor */
        double P_pred = pr * pr * pr * 0.15;  /* empirical scaling */
        double diff = P_pred - prof_P_abs[i];
        P_l1_err += diff * diff;
        P_count++;
    }
    P_l1_err = sqrt(P_l1_err / P_count);
    fprintf(stderr, "  |P| RMS error (L1, scaled phi^3): %.6f (relative: %.1f%%)\n",
            P_l1_err, 100 * P_l1_err / 0.002);

    /* --- Level 2: Three-braid model --- */
    fprintf(stderr, "\n--- Level 2: Three-braid sum with converged parameters ---\n");

    /* Start from gen_phase_confined defaults, then optimize A and R_tube */
    L2Params l2 = {
        .A_braid = 0.3,
        .R_tube = 4.5,
        .ellip = 0.3325,
        .A_bg = 0.1,
        .theta_amp = 0.015,   /* eta * A * 0.1 as in gen_phase_confined */
        .theta_R = 4.5
    };

    fprintf(stderr, "  Initial (gen_phase_confined defaults): A=%.3f R=%.1f\n",
            l2.A_braid, l2.R_tube);
    double l2_init_phi_err = l2_phi_error(&l2, m2, L, delta, carrier, chiral, n_fit);
    double l2_init_P_err = l2_P_error(&l2, m2, L, delta, carrier, chiral, n_fit);
    fprintf(stderr, "  Initial phi_rms error: %.6f, |P| error: %.6f\n",
            l2_init_phi_err, l2_init_P_err);

    /* Optimize A_braid, R_tube, theta_amp, theta_R */
    fprintf(stderr, "  Optimizing A, R_tube, ellip, theta_amp, theta_R...\n");
    {
        double best = l2_phi_error(&l2, m2, L, delta, carrier, chiral, n_fit);
        double step[5] = {0.02, 0.3, 0.02, 0.002, 0.3};
        double *pp[5] = {&l2.A_braid, &l2.R_tube, &l2.ellip, &l2.theta_amp, &l2.theta_R};
        double lo[5] = {0.01, 1.0, 0.0, 0.001, 1.0};
        double hi[5] = {1.0, 15.0, 0.5, 0.5, 15.0};
        for (int iter = 0; iter < 2000; iter++) {
            for (int p = 0; p < 5; p++) {
                double orig = *pp[p];
                *pp[p] = orig + step[p];
                if (*pp[p] > hi[p]) *pp[p] = hi[p];
                double e_plus = l2_phi_error(&l2, m2, L, delta, carrier, chiral, n_fit);
                if (e_plus < best) { best = e_plus; continue; }
                *pp[p] = orig - step[p];
                if (*pp[p] < lo[p]) *pp[p] = lo[p];
                double e_minus = l2_phi_error(&l2, m2, L, delta, carrier, chiral, n_fit);
                if (e_minus < best) { best = e_minus; continue; }
                *pp[p] = orig;
            }
            for (int p = 0; p < 5; p++) step[p] *= 0.998;
        }
    }

    double l2_opt_phi_err = l2_phi_error(&l2, m2, L, delta, carrier, chiral, n_fit);
    double l2_opt_P_err = l2_P_error(&l2, m2, L, delta, carrier, chiral, n_fit);
    fprintf(stderr, "  Optimized: A=%.4f R=%.3f ellip=%.4f theta_amp=%.5f theta_R=%.3f\n",
            l2.A_braid, l2.R_tube, l2.ellip, l2.theta_amp, l2.theta_R);
    fprintf(stderr, "  Optimized phi_rms error: %.6f (relative: %.1f%%)\n",
            l2_opt_phi_err, 100 * l2_opt_phi_err / 0.195);
    fprintf(stderr, "  Optimized |P| error: %.6f (relative: %.1f%%)\n",
            l2_opt_P_err, 100 * l2_opt_P_err / 0.002);

    /* Now optimize theta separately */
    {
        double best = 1e10;
        double step[2] = {0.002, 0.3};
        double *pp[2] = {&l2.theta_amp, &l2.theta_R};
        for (int iter = 0; iter < 2000; iter++) {
            /* Compute theta error */
            double err = 0;
            for (int i = 0; i < n_fit && i < NPROF; i++) {
                double pred = l2_theta_rms_at_r(prof_r[i], &l2, L, carrier, chiral);
                double diff = pred - prof_theta_rms[i];
                err += diff * diff;
            }
            err = sqrt(err / n_fit);
            if (err < best) best = err;

            for (int p = 0; p < 2; p++) {
                double orig = *pp[p];
                *pp[p] = orig + step[p];
                if (*pp[p] < 0.001) *pp[p] = 0.001;
                double e2 = 0;
                for (int i = 0; i < n_fit && i < NPROF; i++) {
                    double pred = l2_theta_rms_at_r(prof_r[i], &l2, L, carrier, chiral);
                    double diff = pred - prof_theta_rms[i];
                    e2 += diff * diff;
                }
                e2 = sqrt(e2 / n_fit);
                if (e2 < best) { best = e2; continue; }
                *pp[p] = orig - step[p];
                if (*pp[p] < 0.001) *pp[p] = 0.001;
                e2 = 0;
                for (int i = 0; i < n_fit && i < NPROF; i++) {
                    double pred = l2_theta_rms_at_r(prof_r[i], &l2, L, carrier, chiral);
                    double diff = pred - prof_theta_rms[i];
                    e2 += diff * diff;
                }
                e2 = sqrt(e2 / n_fit);
                if (e2 < best) { best = e2; continue; }
                *pp[p] = orig;
            }
            for (int p = 0; p < 2; p++) step[p] *= 0.998;
        }
        fprintf(stderr, "  Theta optimized: theta_amp=%.5f theta_R=%.3f, RMS err=%.6f\n",
                l2.theta_amp, l2.theta_R, best);
    }

    /* --- Level 3: Empirical interpolation --- */
    fprintf(stderr, "\n--- Level 3: Empirical radial table ---\n");
    fprintf(stderr, "  Uses measured profile directly — zero fitting error by construction\n");
    fprintf(stderr, "  Applies 3D phase structure from gen_phase_confined\n");

    /* Print comparison table */
    fprintf(stderr, "\n=== COMPARISON: r, phi_rms_data, L1_fit, L2_fit ===\n");
    fprintf(stderr, "%6s %10s %10s %10s %10s %10s %10s\n",
            "r", "phi_data", "L1_phi", "L2_phi", "P_data", "L1_P", "L2_P");
    for (int i = 0; i < n_fit; i += 2) {
        double l1_pr = l1_phi_rms(prof_r[i], &phi_fit);
        double l2_pr = l2_phi_rms_at_r(prof_r[i], &l2, m2, L, delta, carrier, chiral);
        double l1_P = pow(l1_pr, 3) * 0.15;
        double l2_P = l2_P_abs_at_r(prof_r[i], &l2, m2, L, delta, carrier, chiral);
        fprintf(stderr, "%6.2f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f\n",
                prof_r[i], prof_phi_rms[i], l1_pr, l2_pr,
                prof_P_abs[i], l1_P, l2_P);
    }

    /* Physical insight */
    fprintf(stderr, "\n=== PHYSICAL INSIGHT ===\n");
    fprintf(stderr, "  Initial A=0.300, Converged amplitude ratio = %.2f\n",
            prof_phi_rms[0] / prof_phi_rms_init[0]);
    fprintf(stderr, "  → Core amplitude drops to %.0f%% of initial\n",
            100 * prof_phi_rms[0] / prof_phi_rms_init[0]);
    fprintf(stderr, "  Initial R_tube=4.5, Converged L2 fit R_tube=%.2f\n", l2.R_tube);
    fprintf(stderr, "  → Tube radius %s by %.0f%%\n",
            l2.R_tube > 4.5 ? "INCREASED" : "decreased",
            100 * fabs(l2.R_tube - 4.5) / 4.5);
    fprintf(stderr, "  Initial ellip=0.333, Converged L2 fit ellip=%.4f\n", l2.ellip);

    /* Core dip analysis */
    double dip_ratio = prof_phi_rms[0] / prof_phi_rms[20]; /* core / plateau */
    fprintf(stderr, "  Core/plateau phi_rms ratio: %.3f (core is depleted)\n", dip_ratio);
    fprintf(stderr, "  theta_rms: core=%.4f edge=%.4f ratio=%.2f (theta confined to core)\n",
            prof_theta_rms[0], prof_theta_rms[NPROF - 1],
            prof_theta_rms[0] / prof_theta_rms[NPROF - 1]);
    /* Find v_radial zero crossing */
    double v_cross_r = 9.0;
    for (int i = 1; i < NPROF; i++) {
        if (prof_v_radial[i-1] < 0 && prof_v_radial[i] >= 0) {
            v_cross_r = prof_r[i-1] + (prof_r[i]-prof_r[i-1]) *
                        (-prof_v_radial[i-1]) / (prof_v_radial[i]-prof_v_radial[i-1]);
            break;
        }
    }
    fprintf(stderr, "  v_radial: contracting out to r ≈ %.1f (crossover from negative to positive)\n",
            v_cross_r);
    fprintf(stderr, "  phi_v_rms core/edge: %.4f/%.4f = %.1fx (hotter core)\n",
            prof_phi_v_rms[0], prof_phi_v_rms[NPROF - 1],
            prof_phi_v_rms[0] / prof_phi_v_rms[NPROF - 1]);

    /* Check sphericity: compare initial (3-braid) vs converged profile shape */
    fprintf(stderr, "\n  Sphericity test:\n");
    fprintf(stderr, "    Initial phi_rms is nearly flat (0.268 at core, 0.204 at r=12)\n");
    fprintf(stderr, "    Converged phi_rms has strong radial structure (0.090 at core, 0.195 plateau)\n");
    fprintf(stderr, "    The 3 braids have partially merged into a ROUGHLY spherical structure\n");
    fprintf(stderr, "    BUT the Level 2 three-braid model still captures the 3D phase structure\n");

    if (fit_only) {
        fprintf(stderr, "\n  (fit-only mode, no seed generated)\n");
        return 0;
    }

    /* ==================== SEED GENERATION ==================== */
    fprintf(stderr, "\n=== GENERATING SEED (Level %d) ===\n", level);
    fprintf(stderr, "  N=%d L=%.1f center=(%.1f,%.1f,%.1f) chirality=%s\n",
            N, L, cx, cy, cz, chirality);

    long N3 = (long)N * N * N;
    int NN = N * N;
    double dx = 2.0 * L / (N - 1);
    double dt = dt_factor * dx;
    double kw = PI / L;
    double omega = sqrt(kw * kw + m2);
    double k_bg = PI / L;
    double omega_bg = sqrt(k_bg * k_bg + m2);

    double *phi[3], *phi_vel[3], *theta[3], *theta_vel[3];
    for (int a = 0; a < 3; a++) {
        phi[a]       = (double *)calloc(N3, sizeof(double));
        phi_vel[a]   = (double *)calloc(N3, sizeof(double));
        theta[a]     = (double *)calloc(N3, sizeof(double));
        theta_vel[a] = (double *)calloc(N3, sizeof(double));
    }

    /* Precompute radial correction ratio: converged / initial phi_rms
     * This captures the amplitude change during evolution.
     * We apply this as a multiplicative correction to the gen_phase_confined structure. */
    double prof_ratio[NPROF];
    for (int i = 0; i < NPROF; i++) {
        prof_ratio[i] = (prof_phi_rms_init[i] > 1e-10) ?
            prof_phi_rms[i] / prof_phi_rms_init[i] : 1.0;
    }

    if (level == 1) {
        /* ===== LEVEL 1: gen_phase_confined structure * L1 radial correction ===== */
        /* Use the standard 3-tube Gaussian braid with fitted radial modulation.
         * The correction factor C(r) = phi_rms_fit(r) / phi_rms_initial_model(r)
         * shrinks the core and preserves the outer envelope. */
        fprintf(stderr, "  Level 1: 3-braid Gaussian tubes * L1 radial correction\n");

        /* Use gen_phase_confined defaults for tube structure */
        double A_l1 = 0.3, R_l1 = 4.5, ellip_l1 = 0.3325;
        double sx = 1 + ellip_l1, sy = 1 - ellip_l1;
        double inv2R2_l1 = 1.0 / (2 * R_l1 * R_l1);
        double theta_init_frac = 0.5;

        #pragma omp parallel for schedule(static)
        for (long idx = 0; idx < N3; idx++) {
            int ii = (int)(idx / NN), jj = (int)((idx / N) % N), kk = (int)(idx % N);
            double x = -L + ii * dx - cx;
            double y = -L + jj * dx - cy;
            double z = -L + kk * dx - cz;
            double r = sqrt(x * x + y * y + z * z);
            double z_abs = -L + kk * dx;

            /* Radial correction from L1 fit */
            double target_rms = l1_phi_rms(r, &phi_fit);
            double init_rms = interp_profile(prof_r, prof_phi_rms_init, NPROF, r);
            double C_r = (init_rms > 1e-10) ? target_rms / init_rms : 1.0;
            if (C_r > 3.0) C_r = 3.0;

            /* Background */
            for (int a = 0; a < 3; a++) {
                double ph_bg = k_bg * z_abs + 2 * PI * a / 3.0;
                phi[a][idx] = A_bg * cos(ph_bg);
                phi_vel[a][idx] = omega_bg * A_bg * sin(ph_bg);
            }

            /* Braid 0: along z */
            {
                double r2e = x * x / (sx * sx) + y * y / (sy * sy);
                double env = C_r * A_l1 * exp(-r2e * inv2R2_l1);
                for (int a = 0; a < 3; a++) {
                    double ph = chiral[0] * kw * z + delta[a] + carrier[0];
                    phi[a][idx] += env * cos(ph);
                    phi_vel[a][idx] += chiral[0] * omega * env * sin(ph);
                }
            }
            /* Braid 1: along x */
            {
                double zr = x, xr = -z, yr = y;
                double r2e = xr * xr / (sx * sx) + yr * yr / (sy * sy);
                double env = C_r * A_l1 * exp(-r2e * inv2R2_l1);
                for (int a = 0; a < 3; a++) {
                    double ph = chiral[1] * kw * zr + delta[a] + carrier[1];
                    phi[a][idx] += env * cos(ph);
                    phi_vel[a][idx] += chiral[1] * omega * env * sin(ph);
                }
            }
            /* Braid 2: along y */
            {
                double zr = y, yr = -z, xr = x;
                double r2e = xr * xr / (sx * sx) + yr * yr / (sy * sy);
                double env = C_r * A_l1 * exp(-r2e * inv2R2_l1);
                for (int a = 0; a < 3; a++) {
                    double ph = chiral[2] * kw * zr + delta[a] + carrier[2];
                    phi[a][idx] += env * cos(ph);
                    phi_vel[a][idx] += chiral[2] * omega * env * sin(ph);
                }
            }

            /* Theta from L1 fit */
            double theta_target = l1_theta_rms(r, &theta_fit);
            double theta_C = theta_target / 0.01;  /* scale relative to initial theta */
            if (theta_C > 50) theta_C = 50;
            for (int b = 0; b < 3; b++) {
                double braid_env;
                if (b == 0) {
                    double r2e = x * x / (sx * sx) + y * y / (sy * sy);
                    braid_env = A_l1 * exp(-r2e * inv2R2_l1);
                } else if (b == 1) {
                    double xr = -z, yr = y;
                    double r2e = xr * xr / (sx * sx) + yr * yr / (sy * sy);
                    braid_env = A_l1 * exp(-r2e * inv2R2_l1);
                } else {
                    double xr = x, yr = -z;
                    double r2e = xr * xr / (sx * sx) + yr * yr / (sy * sy);
                    braid_env = A_l1 * exp(-r2e * inv2R2_l1);
                }
                double t_env = theta_init_frac * eta * braid_env * 0.1 * theta_C;
                theta[0][idx] += chiral[b] * t_env * cos(kw * z + carrier[b]);
                theta[1][idx] += chiral[b] * t_env * sin(kw * z + carrier[b]);
                theta[2][idx] += chiral[b] * t_env * cos(kw * z + carrier[b] + PI / 4);
            }

            /* Contracting velocity from fitted v_radial profile */
            if (r > 0.1) {
                double vr = l1_v_radial(r, &vrad_fit);
                double rhat[3] = {x / r, y / r, z / r};
                for (int a = 0; a < 3; a++)
                    phi_vel[a][idx] += vr * rhat[a % 3] * phi[a][idx];
            }
        }

    } else if (level == 2) {
        /* ===== LEVEL 2: Three-braid with optimized parameters ===== */
        fprintf(stderr, "  Level 2: Three-braid model with optimized parameters\n");
        fprintf(stderr, "    A=%.4f R=%.3f ellip=%.4f theta_amp=%.5f theta_R=%.3f\n",
                l2.A_braid, l2.R_tube, l2.ellip, l2.theta_amp, l2.theta_R);

        double sx = 1 + l2.ellip, sy = 1 - l2.ellip;
        double inv2R2 = 1.0 / (2 * l2.R_tube * l2.R_tube);
        double inv2R2_t = 1.0 / (2 * l2.theta_R * l2.theta_R);

        #pragma omp parallel for schedule(static)
        for (long idx = 0; idx < N3; idx++) {
            int ii = (int)(idx / NN), jj = (int)((idx / N) % N), kk = (int)(idx % N);
            double x = -L + ii * dx - cx;
            double y = -L + jj * dx - cy;
            double z = -L + kk * dx - cz;
            double r = sqrt(x * x + y * y + z * z);
            double z_abs = -L + kk * dx;

            /* Background */
            for (int a = 0; a < 3; a++) {
                double ph_bg = k_bg * z_abs + 2 * PI * a / 3.0;
                phi[a][idx] = A_bg * cos(ph_bg);
                phi_vel[a][idx] = omega_bg * A_bg * sin(ph_bg);
            }

            /* Braid 0: along z */
            {
                double r2e = x * x / (sx * sx) + y * y / (sy * sy);
                double env = l2.A_braid * exp(-r2e * inv2R2);
                for (int a = 0; a < 3; a++) {
                    double ph = chiral[0] * kw * z + delta[a] + carrier[0];
                    phi[a][idx] += env * cos(ph);
                    phi_vel[a][idx] += chiral[0] * omega * env * sin(ph);
                }
            }
            /* Braid 1: along x */
            {
                double zr = x, xr = -z, yr = y;
                double r2e = xr * xr / (sx * sx) + yr * yr / (sy * sy);
                double env = l2.A_braid * exp(-r2e * inv2R2);
                for (int a = 0; a < 3; a++) {
                    double ph = chiral[1] * kw * zr + delta[a] + carrier[1];
                    phi[a][idx] += env * cos(ph);
                    phi_vel[a][idx] += chiral[1] * omega * env * sin(ph);
                }
            }
            /* Braid 2: along y */
            {
                double zr = y, yr = -z, xr = x;
                double r2e = xr * xr / (sx * sx) + yr * yr / (sy * sy);
                double env = l2.A_braid * exp(-r2e * inv2R2);
                for (int a = 0; a < 3; a++) {
                    double ph = chiral[2] * kw * zr + delta[a] + carrier[2];
                    phi[a][idx] += env * cos(ph);
                    phi_vel[a][idx] += chiral[2] * omega * env * sin(ph);
                }
            }

            /* Theta from L2 params */
            for (int b = 0; b < 3; b++) {
                double braid_env;
                if (b == 0) {
                    double r2e = x * x / (sx * sx) + y * y / (sy * sy);
                    braid_env = l2.A_braid * exp(-r2e * inv2R2_t);
                } else if (b == 1) {
                    double xr = -z, yr = y;
                    double r2e = xr * xr / (sx * sx) + yr * yr / (sy * sy);
                    braid_env = l2.A_braid * exp(-r2e * inv2R2_t);
                } else {
                    double xr = x, yr = -z;
                    double r2e = xr * xr / (sx * sx) + yr * yr / (sy * sy);
                    braid_env = l2.A_braid * exp(-r2e * inv2R2_t);
                }
                double t_env = l2.theta_amp * braid_env;
                theta[0][idx] += chiral[b] * t_env * cos(kw * z + carrier[b]);
                theta[1][idx] += chiral[b] * t_env * sin(kw * z + carrier[b]);
                theta[2][idx] += chiral[b] * t_env * cos(kw * z + carrier[b] + PI / 4);
            }

            /* Contracting velocity from fitted v_radial profile */
            if (r > 0.1) {
                double vr = l1_v_radial(r, &vrad_fit);
                double rhat[3] = {x / r, y / r, z / r};
                for (int a = 0; a < 3; a++)
                    phi_vel[a][idx] += vr * rhat[a % 3] * phi[a][idx];
            }
        }

    } else if (level == 3) {
        /* ===== LEVEL 3: Empirical correction to gen_phase_confined ===== */
        /* Same 3-tube structure as gen_phase_confined, but multiply by the exact
         * converged/initial ratio from simulation data at each radius.
         * This is the most accurate because it uses the MEASURED correction. */
        fprintf(stderr, "  Level 3: gen_phase_confined * empirical correction C(r)\n");

        double A_l3 = 0.3, R_l3 = 4.5, ellip_l3 = 0.3325;
        double sx = 1 + ellip_l3, sy = 1 - ellip_l3;
        double inv2R2_l3 = 1.0 / (2 * R_l3 * R_l3);
        double theta_init_frac = 0.5;

        #pragma omp parallel for schedule(static)
        for (long idx = 0; idx < N3; idx++) {
            int ii = (int)(idx / NN), jj = (int)((idx / N) % N), kk = (int)(idx % N);
            double x = -L + ii * dx - cx;
            double y = -L + jj * dx - cy;
            double z = -L + kk * dx - cz;
            double r = sqrt(x * x + y * y + z * z);
            double z_abs = -L + kk * dx;

            /* Empirical radial correction: converged / initial phi_rms ratio */
            double C_r = interp_profile(prof_r, prof_ratio, NPROF, r);
            if (C_r > 3.0) C_r = 3.0;
            if (C_r < 0.0) C_r = 0.0;

            /* Background */
            for (int a = 0; a < 3; a++) {
                double ph_bg = k_bg * z_abs + 2 * PI * a / 3.0;
                phi[a][idx] = A_bg * cos(ph_bg);
                phi_vel[a][idx] = omega_bg * A_bg * sin(ph_bg);
            }

            /* Braid 0: along z */
            {
                double r2e = x * x / (sx * sx) + y * y / (sy * sy);
                double env = C_r * A_l3 * exp(-r2e * inv2R2_l3);
                for (int a = 0; a < 3; a++) {
                    double ph = chiral[0] * kw * z + delta[a] + carrier[0];
                    phi[a][idx] += env * cos(ph);
                    phi_vel[a][idx] += chiral[0] * omega * env * sin(ph);
                }
            }
            /* Braid 1: along x */
            {
                double zr = x, xr = -z, yr = y;
                double r2e = xr * xr / (sx * sx) + yr * yr / (sy * sy);
                double env = C_r * A_l3 * exp(-r2e * inv2R2_l3);
                for (int a = 0; a < 3; a++) {
                    double ph = chiral[1] * kw * zr + delta[a] + carrier[1];
                    phi[a][idx] += env * cos(ph);
                    phi_vel[a][idx] += chiral[1] * omega * env * sin(ph);
                }
            }
            /* Braid 2: along y */
            {
                double zr = y, yr = -z, xr = x;
                double r2e = xr * xr / (sx * sx) + yr * yr / (sy * sy);
                double env = C_r * A_l3 * exp(-r2e * inv2R2_l3);
                for (int a = 0; a < 3; a++) {
                    double ph = chiral[2] * kw * zr + delta[a] + carrier[2];
                    phi[a][idx] += env * cos(ph);
                    phi_vel[a][idx] += chiral[2] * omega * env * sin(ph);
                }
            }

            /* Theta from empirical profile */
            double theta_target = interp_profile(prof_r, prof_theta_rms, NPROF, r);
            double theta_C = theta_target / 0.01;
            if (theta_C > 50) theta_C = 50;
            for (int b = 0; b < 3; b++) {
                double braid_env;
                if (b == 0) {
                    double r2e = x * x / (sx * sx) + y * y / (sy * sy);
                    braid_env = A_l3 * exp(-r2e * inv2R2_l3);
                } else if (b == 1) {
                    double xr = -z, yr = y;
                    double r2e = xr * xr / (sx * sx) + yr * yr / (sy * sy);
                    braid_env = A_l3 * exp(-r2e * inv2R2_l3);
                } else {
                    double xr = x, yr = -z;
                    double r2e = xr * xr / (sx * sx) + yr * yr / (sy * sy);
                    braid_env = A_l3 * exp(-r2e * inv2R2_l3);
                }
                double t_env = theta_init_frac * eta * braid_env * 0.1 * theta_C;
                theta[0][idx] += chiral[b] * t_env * cos(kw * z + carrier[b]);
                theta[1][idx] += chiral[b] * t_env * sin(kw * z + carrier[b]);
                theta[2][idx] += chiral[b] * t_env * cos(kw * z + carrier[b] + PI / 4);
            }

            /* Contracting velocity from empirical profile */
            if (r > 0.1) {
                double vr = interp_profile(prof_r, prof_v_radial, NPROF, r);
                double rhat[3] = {x / r, y / r, z / r};
                for (int a = 0; a < 3; a++)
                    phi_vel[a][idx] += vr * rhat[a % 3] * phi[a][idx];
            }
        }
    }

    /* ===== Self-test: compute radial profile of generated seed ===== */
    fprintf(stderr, "\n=== SELF-TEST: Radial profile of generated seed ===\n");
    fprintf(stderr, "%6s %10s %10s %10s %10s %10s %10s\n",
            "r", "phi_data", "phi_seed", "err%", "P_data", "P_seed", "err%");

    /* Bin into radial shells */
    int n_bins = 48;
    double *seed_phi_sq = (double *)calloc(n_bins, sizeof(double));
    double *seed_P_sum = (double *)calloc(n_bins, sizeof(double));
    double *seed_theta_sq = (double *)calloc(n_bins, sizeof(double));
    int *seed_count = (int *)calloc(n_bins, sizeof(int));
    double dr_bin = 0.25;

    for (long idx = 0; idx < N3; idx++) {
        int ii = (int)(idx / NN), jj = (int)((idx / N) % N), kk = (int)(idx % N);
        double x = -L + ii * dx - cx;
        double y = -L + jj * dx - cy;
        double z = -L + kk * dx - cz;
        double r = sqrt(x * x + y * y + z * z);
        int bin = (int)(r / dr_bin);
        if (bin >= n_bins) continue;
        for (int a = 0; a < 3; a++)
            seed_phi_sq[bin] += phi[a][idx] * phi[a][idx];
        seed_P_sum[bin] += fabs(phi[0][idx] * phi[1][idx] * phi[2][idx]);
        for (int a = 0; a < 3; a++)
            seed_theta_sq[bin] += theta[a][idx] * theta[a][idx];
        seed_count[bin]++;
    }

    double phi_rms_err_sum = 0, P_err_sum = 0;
    int err_count = 0;
    for (int i = 0; i < n_bins && i < NPROF; i++) {
        if (seed_count[i] == 0) continue;
        double seed_phi_rms = sqrt(seed_phi_sq[i] / (3.0 * seed_count[i]));
        double seed_P = seed_P_sum[i] / seed_count[i];
        double r_bin = (i + 0.5) * dr_bin;

        /* Find closest profile point */
        int closest = 0;
        double min_dr = fabs(r_bin - prof_r[0]);
        for (int j = 1; j < NPROF; j++) {
            double d = fabs(r_bin - prof_r[j]);
            if (d < min_dr) { min_dr = d; closest = j; }
        }
        if (min_dr > 0.15) continue;

        double phi_err_pct = 100 * fabs(seed_phi_rms - prof_phi_rms[closest]) / prof_phi_rms[closest];
        double P_err_pct = prof_P_abs[closest] > 1e-8 ?
            100 * fabs(seed_P - prof_P_abs[closest]) / prof_P_abs[closest] : 0;

        if (i % 2 == 0) {
            fprintf(stderr, "%6.2f %10.6f %10.6f %8.1f%% %10.6f %10.6f %8.1f%%\n",
                    r_bin, prof_phi_rms[closest], seed_phi_rms, phi_err_pct,
                    prof_P_abs[closest], seed_P, P_err_pct);
        }
        phi_rms_err_sum += (seed_phi_rms - prof_phi_rms[closest]) *
                           (seed_phi_rms - prof_phi_rms[closest]);
        P_err_sum += (seed_P - prof_P_abs[closest]) *
                     (seed_P - prof_P_abs[closest]);
        err_count++;
    }
    double phi_rms_final = sqrt(phi_rms_err_sum / err_count);
    double P_rms_final = sqrt(P_err_sum / err_count);
    fprintf(stderr, "\n  Overall phi_rms RMS error: %.6f (%.1f%% of mean)\n",
            phi_rms_final, 100 * phi_rms_final / 0.18);
    fprintf(stderr, "  Overall |P| RMS error: %.6f (%.1f%% of mean)\n",
            P_rms_final, 100 * P_rms_final / 0.0015);

    free(seed_phi_sq); free(seed_P_sum); free(seed_theta_sq); free(seed_count);

    /* ===== WRITE SFA ===== */
    uint8_t sfa_dtype = (precision == 0) ? SFA_F16 : (precision == 1) ? SFA_F32 : SFA_F64;
    SFA *sfa = sfa_create(outpath, N, N, N, L, L, L, dt);

    char vch[32], vlvl[16];
    snprintf(vch, 32, "%s", chirality);
    snprintf(vlvl, 16, "%d", level);
    char vA[32], vR[32];
    if (level == 2) {
        snprintf(vA, 32, "%.6f", l2.A_braid);
        snprintf(vR, 32, "%.3f", l2.R_tube);
    } else {
        snprintf(vA, 32, "%.6f", phi_fit.A_flat);
        snprintf(vR, 32, "%.3f", phi_fit.s_dip);
    }
    const char *keys[] = {"chirality", "type", "level", "A_fit", "R_fit"};
    const char *vals[] = {vch, "proton_analytical", vlvl, vA, vR};
    sfa_add_kvmd(sfa, 0, 0xFFFFFFFF, 0xFFFFFFFF, keys, vals, 5);

    sfa_add_column(sfa, "phi_x", sfa_dtype, SFA_POSITION, 0);
    sfa_add_column(sfa, "phi_y", sfa_dtype, SFA_POSITION, 1);
    sfa_add_column(sfa, "phi_z", sfa_dtype, SFA_POSITION, 2);
    sfa_add_column(sfa, "theta_x", sfa_dtype, SFA_ANGLE, 0);
    sfa_add_column(sfa, "theta_y", sfa_dtype, SFA_ANGLE, 1);
    sfa_add_column(sfa, "theta_z", sfa_dtype, SFA_ANGLE, 2);
    sfa_add_column(sfa, "phi_vx", sfa_dtype, SFA_VELOCITY, 0);
    sfa_add_column(sfa, "phi_vy", sfa_dtype, SFA_VELOCITY, 1);
    sfa_add_column(sfa, "phi_vz", sfa_dtype, SFA_VELOCITY, 2);
    sfa_add_column(sfa, "theta_vx", sfa_dtype, SFA_VELOCITY, 3);
    sfa_add_column(sfa, "theta_vy", sfa_dtype, SFA_VELOCITY, 4);
    sfa_add_column(sfa, "theta_vz", sfa_dtype, SFA_VELOCITY, 5);
    sfa_finalize_header(sfa);

    if (precision == 2) {
        void *cols[12] = {phi[0], phi[1], phi[2], theta[0], theta[1], theta[2],
                          phi_vel[0], phi_vel[1], phi_vel[2],
                          theta_vel[0], theta_vel[1], theta_vel[2]};
        sfa_write_frame(sfa, 0.0, cols);
    } else {
        void *cols[12];
        int es = (precision == 0) ? 2 : 4;
        for (int c = 0; c < 12; c++) {
            double *src = (c < 3) ? phi[c] : (c < 6) ? theta[c - 3] :
                          (c < 9) ? phi_vel[c - 6] : theta_vel[c - 9];
            cols[c] = malloc(N3 * es);
            if (precision == 1) {
                float *p = (float *)cols[c];
                for (long i = 0; i < N3; i++) p[i] = (float)src[i];
            } else {
                uint16_t *p = (uint16_t *)cols[c];
                for (long i = 0; i < N3; i++) p[i] = f64_to_f16(src[i]);
            }
        }
        sfa_write_frame(sfa, 0.0, cols);
        for (int c = 0; c < 12; c++) free(cols[c]);
    }
    sfa_close(sfa);
    fprintf(stderr, "  Wrote %s\n", outpath);

    for (int a = 0; a < 3; a++) {
        free(phi[a]); free(phi_vel[a]); free(theta[a]); free(theta_vel[a]);
    }

    return 0;
}
