/*  scp_sim_fixed.c — Fixed-point int32 6-field Cosserat simulation kernel
 *
 *  Same PDE as scp_sim.c, same leapfrog integrator, same stencil.
 *  ALL physics computation in int32/int64 fixed-point (Q16.16 format).
 *  Floating point is used ONLY at I/O boundaries (config parse, SFA output).
 *
 *  Physics: Variant B v2 unified transfer potential
 *    V_base(P) = (mu/2) * P^2 / (1 + kappa*P^2)
 *    f_transfer = epsilon + (1-epsilon) * Theta/(Theta+Theta_c)
 *    W_confine = |V_base| * gamma * ln(1 + Theta/Theta_c)
 *
 *  Build: gcc -O3 -march=native -fopenmp -o scp_sim_fixed scp_sim_fixed.c -lzstd -lm
 *  Run:   OMP_NUM_THREADS=8 ./scp_sim_fixed config.cfg [-key value ...]
 */

#define SFA_IMPLEMENTATION
#include "../../sfa/format/sfa.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>
#include <sys/stat.h>
#include <unistd.h>

#define NFIELDS 3
#define PI 3.14159265358979323846

/* ================================================================
   Fixed-point Q16.16 arithmetic
   ================================================================ */

#define FP_SHIFT  16
#define FP_ONE    (1 << FP_SHIFT)       /* 65536 = 1.0 */
#define FP_HALF   (1 << (FP_SHIFT - 1)) /* 32768 = 0.5 */

/* Convert float <-> fixed-point */
static inline int32_t fp_from_f64(double v) {
    double scaled = v * (double)FP_ONE;
    if (scaled > 2147483647.0) return 2147483647;
    if (scaled < -2147483648.0) return -2147483647 - 1;
    return (int32_t)lround(scaled);
}

static inline double fp_to_f64(int32_t v) {
    return (double)v / (double)FP_ONE;
}

/* Fixed-point multiply: (a * b) >> 16 using int64 intermediate */
static inline int32_t fp_mul(int32_t a, int32_t b) {
    return (int32_t)(((int64_t)a * (int64_t)b) >> FP_SHIFT);
}

/* Fixed-point divide: (a << 16) / b using int64 intermediate */
static inline int32_t fp_div(int32_t a, int32_t b) {
    if (b == 0) return (a >= 0) ? 2147483647 : -2147483647;
    return (int32_t)(((int64_t)a << FP_SHIFT) / (int64_t)b);
}

/* Fixed-point square */
static inline int32_t fp_sq(int32_t a) {
    return (int32_t)(((int64_t)a * (int64_t)a) >> FP_SHIFT);
}

/* Fixed-point absolute value */
static inline int32_t fp_abs(int32_t a) {
    return (a >= 0) ? a : -a;
}

/* Fixed-point sqrt via Newton's method.
 * Input: Q16.16 non-negative. Output: Q16.16.
 * sqrt(x_fp) in Q16.16: result = sqrt(x * 65536) = sqrt(x) * 256
 * We compute sqrt(x_fp << 16) to get Q16.16 result. */
static inline int32_t fp_sqrt(int32_t x) {
    if (x <= 0) return 0;
    /* We need sqrt(x) where x is in Q16.16.
     * sqrt_fp(x_fp) = sqrt(x_fp * FP_ONE) = sqrt(x_fp) * 256
     * Use int64 to avoid overflow: val = x_fp << 16 */
    uint64_t val = (uint64_t)(uint32_t)x << FP_SHIFT;
    uint64_t r = val;
    /* Initial guess: half the bit length */
    uint64_t guess = 1;
    uint64_t tmp = val;
    while (tmp > 1) { tmp >>= 2; guess <<= 1; }
    r = guess;
    for (int i = 0; i < 20; i++) {
        if (r == 0) break;
        uint64_t next = (r + val / r) >> 1;
        if (next >= r) break;
        r = next;
    }
    return (int32_t)r;
}

/* Fixed-point natural log approximation.
 * ln(x) for x in Q16.16, x > 0. Result in Q16.16.
 *
 * Method: reduce x to [1,2) by extracting integer part of log2,
 * then use polynomial on the reduced range.
 *   ln(x) = ln(2^n * m) = n*ln(2) + ln(m)  where m in [1,2)
 *   ln(m) ≈ (m-1) - (m-1)^2/2 + (m-1)^3/3  for m in [1,2)
 */
static int32_t fp_ln(int32_t x) {
    if (x <= 0) return -2147483647;  /* -inf approx */
    if (x == FP_ONE) return 0;

    /* Find n such that x = 2^n * m, m in [FP_ONE, 2*FP_ONE) */
    int n = 0;
    int32_t m = x;
    while (m >= 2 * FP_ONE) { m >>= 1; n++; }
    while (m < FP_ONE) { m <<= 1; n--; }
    /* Now m is in [FP_ONE, 2*FP_ONE), n is the exponent */

    /* u = m - 1.0 in Q16.16 */
    int32_t u = m - FP_ONE;  /* u in [0, FP_ONE) */

    /* ln(m) = u - u^2/2 + u^3/3 - u^4/4 + u^5/5 - u^6/6 + u^7/7
     * (Maclaurin series for ln(1+u), |u|<1)
     *
     * TODO: PRECISION — 4-term series has ~16% relative error at u->1
     * (the boundary where m is just below 2*FP_ONE). Using 7 terms
     * reduces worst-case error to ~3%. For production, consider a
     * minimax polynomial fit over [0,1) which can achieve <0.1% with
     * the same number of terms. */
    int32_t u2 = fp_mul(u, u);
    int32_t u3 = fp_mul(u2, u);
    int32_t u4 = fp_mul(u3, u);
    int32_t u5 = fp_mul(u4, u);
    int32_t u6 = fp_mul(u5, u);
    int32_t u7 = fp_mul(u6, u);

    /* ln(m) = u - u2/2 + u3/3 - u4/4 + u5/5 - u6/6 + u7/7 */
    int32_t ln_m = u - (u2 >> 1)
                 + fp_div(u3, 3 * FP_ONE) - (u4 >> 2)
                 + fp_div(u5, 5 * FP_ONE) - fp_div(u6, 6 * FP_ONE)
                 + fp_div(u7, 7 * FP_ONE);

    /* ln(2) in Q16.16 = 0.693147... * 65536 = 45426 */
    const int32_t LN2_FP = 45426;
    int32_t result = (int32_t)((int64_t)n * LN2_FP) + ln_m;
    return result;
}

/* Fixed-point ln(1 + x/y) = ln((y+x)/y) = ln(y+x) - ln(y)
 * More numerically stable version for the confinement term */
static inline int32_t fp_ln1p_ratio(int32_t x, int32_t y) {
    /* ln(1 + x/y) where x >= 0, y > 0 */
    if (x == 0) return 0;
    /* For small x/y, use the series: ln(1+r) ~ r - r^2/2 + r^3/3
     * where r = x/y */
    int32_t sum = x + y;
    if (sum <= 0) return 0;  /* overflow protection */
    return fp_ln(sum) - fp_ln(y);
}

/* ================================================================
   Configuration — same as scp_sim.c plus variant B params
   ================================================================ */

typedef struct {
    /* Grid */
    int N;
    double L, T, dt_factor;

    /* Physics */
    double m2, mtheta2, eta, mu, kappa;

    /* Variant B v2 params */
    double theta_c;     /* Transfer midpoint */
    double epsilon;     /* Pilot light fraction */
    double gamma_conf;  /* Confinement strength */

    /* Boundary */
    int bc_type;
    double damp_width, damp_rate;
    double bc_switch_time;
    double gradient_A_high, gradient_A_low;
    int gradient_margin;

    /* Init */
    char init[32];
    double A, sigma, A_bg, ellip, R_tube;
    double delta[3];
    char init_sfa[512];
    int init_frame;
    char init_exec[1024];

    /* Output */
    char output[512];
    char diag_file[512];
    int precision;
    int output_split;
    double snap_dt, diag_dt;
} Config;

static Config cfg_defaults(void) {
    Config c = {0};
    c.N = 128;  c.L = 10.0;  c.T = 200.0;  c.dt_factor = 0.025;
    c.m2 = 2.25;  c.mtheta2 = 0.0;  c.eta = 0.5;
    c.mu = -41.345;  c.kappa = 50.0;
    c.theta_c = 0.05;  c.epsilon = 0.1;  c.gamma_conf = 0.1;
    c.bc_type = 0;
    c.damp_width = 3.0;  c.damp_rate = 0.01;  c.bc_switch_time = 0.0;
    c.gradient_A_high = 0.15;  c.gradient_A_low = 0.05;  c.gradient_margin = 3;
    strcpy(c.init, "oscillon");
    c.A = 0.8;  c.sigma = 3.0;  c.A_bg = 0.1;
    c.ellip = 0.3325;  c.R_tube = 3.0;
    c.delta[0] = 0.0;  c.delta[1] = 3.0005;  c.delta[2] = 4.4325;
    c.init_frame = -1;
    strcpy(c.output, "output.sfa");
    strcpy(c.diag_file, "diag.tsv");
    c.precision = 1;
    c.output_split = 0;
    c.snap_dt = 5.0;  c.diag_dt = 2.0;
    return c;
}

static void cfg_set(Config *c, const char *key, const char *val) {
    if      (!strcmp(key,"N"))           c->N = atoi(val);
    else if (!strcmp(key,"L"))           c->L = atof(val);
    else if (!strcmp(key,"T"))           c->T = atof(val);
    else if (!strcmp(key,"dt_factor"))   c->dt_factor = atof(val);
    else if (!strcmp(key,"m"))         { double m = atof(val); c->m2 = m*m; }
    else if (!strcmp(key,"m2"))          c->m2 = atof(val);
    else if (!strcmp(key,"m_theta"))   { double m = atof(val); c->mtheta2 = m*m; }
    else if (!strcmp(key,"mtheta2"))     c->mtheta2 = atof(val);
    else if (!strcmp(key,"eta"))         c->eta = atof(val);
    else if (!strcmp(key,"mu"))          c->mu = atof(val);
    else if (!strcmp(key,"kappa"))       c->kappa = atof(val);
    else if (!strcmp(key,"theta_c"))     c->theta_c = atof(val);
    else if (!strcmp(key,"epsilon"))     c->epsilon = atof(val);
    else if (!strcmp(key,"gamma"))       c->gamma_conf = atof(val);
    else if (!strcmp(key,"bc_type"))     c->bc_type = atoi(val);
    else if (!strcmp(key,"damp_width"))  c->damp_width = atof(val);
    else if (!strcmp(key,"damp_rate"))   c->damp_rate = atof(val);
    else if (!strcmp(key,"bc_switch_time")) c->bc_switch_time = atof(val);
    else if (!strcmp(key,"gradient_A_high")) c->gradient_A_high = atof(val);
    else if (!strcmp(key,"gradient_A_low"))  c->gradient_A_low = atof(val);
    else if (!strcmp(key,"gradient_margin")) c->gradient_margin = atoi(val);
    else if (!strcmp(key,"init"))        strncpy(c->init, val, 31);
    else if (!strcmp(key,"A"))           c->A = atof(val);
    else if (!strcmp(key,"sigma"))       c->sigma = atof(val);
    else if (!strcmp(key,"A_bg"))        c->A_bg = atof(val);
    else if (!strcmp(key,"ellip"))       c->ellip = atof(val);
    else if (!strcmp(key,"R_tube"))      c->R_tube = atof(val);
    else if (!strcmp(key,"delta"))       sscanf(val, "%lf,%lf,%lf", &c->delta[0], &c->delta[1], &c->delta[2]);
    else if (!strcmp(key,"init_sfa"))    strncpy(c->init_sfa, val, 511);
    else if (!strcmp(key,"init_frame"))  c->init_frame = atoi(val);
    else if (!strcmp(key,"init_exec"))   strncpy(c->init_exec, val, 1023);
    else if (!strcmp(key,"output"))      strncpy(c->output, val, 511);
    else if (!strcmp(key,"diag_file"))   strncpy(c->diag_file, val, 511);
    else if (!strcmp(key,"snap_dt"))     c->snap_dt = atof(val);
    else if (!strcmp(key,"diag_dt"))     c->diag_dt = atof(val);
    else if (!strcmp(key,"precision")) {
        if      (!strcmp(val,"f16")) c->precision = 0;
        else if (!strcmp(val,"f32")) c->precision = 1;
        else if (!strcmp(val,"f64")) c->precision = 2;
    }
    else if (!strcmp(key,"output_split")) c->output_split = atoi(val);
    /* Silently ignore old params that don't apply to fixed-point variant B */
    else if (!strcmp(key,"lambda_theta") || !strcmp(key,"eta1") ||
             !strcmp(key,"mode") || !strcmp(key,"inv_alpha") ||
             !strcmp(key,"inv_beta") || !strcmp(key,"kappa_gamma")) {
        /* ignored */
    }
    else fprintf(stderr, "WARNING: unknown config key '%s'\n", key);
}

static void cfg_load(Config *c, const char *path) {
    FILE *fp = fopen(path, "r");
    if (!fp) { fprintf(stderr, "Cannot open config: %s\n", path); exit(1); }
    char line[2048];
    while (fgets(line, sizeof(line), fp)) {
        char *p = line;
        while (*p == ' ' || *p == '\t') p++;
        if (*p == '#' || *p == '\n' || *p == '\0') continue;
        char *eq = strchr(p, '=');
        if (!eq) continue;
        *eq = '\0';
        char *key = p, *val = eq + 1;
        char *ke = eq - 1;
        while (ke > key && (*ke == ' ' || *ke == '\t')) *ke-- = '\0';
        while (*val == ' ' || *val == '\t') val++;
        char *ve = val + strlen(val) - 1;
        while (ve > val && (*ve == '\n' || *ve == '\r' || *ve == ' ')) *ve-- = '\0';
        char *hash = strchr(val, '#');
        if (hash) { *hash = '\0'; ve = hash - 1; while (ve > val && *ve == ' ') *ve-- = '\0'; }
        cfg_set(c, key, val);
    }
    fclose(fp);
}

static void cfg_print(const Config *c) {
    const char *prec_names[] = {"f16", "f32", "f64"};
    printf("=== scp_sim_fixed: Fixed-Point int32 6-field Cosserat (Variant B v2) ===\n");
    printf("d^2phi/dt^2 = lap(phi) - m^2*phi - dVdP*dPda*f_phi + eta*curl(theta)\n");
    printf("d^2theta/dt^2 = lap(theta) - |V_base|*mass_coeff*theta + eta*curl(phi)\n\n");
    printf("Grid:    N=%d L=%.1f T=%.0f dt_factor=%.4f\n", c->N, c->L, c->T, c->dt_factor);
    printf("Physics: m^2=%.4f m_theta^2=%.4f eta=%.3f mu=%.3f kappa=%.1f\n",
           c->m2, c->mtheta2, c->eta, c->mu, c->kappa);
    printf("TransferB: theta_c=%.4f epsilon=%.3f gamma=%.3f\n",
           c->theta_c, c->epsilon, c->gamma_conf);
    if (c->bc_type == 0)
        printf("BC:      absorbing sphere (width=%.1f rate=%.4f)\n", c->damp_width, c->damp_rate);
    else if (c->bc_type == 1)
        printf("BC:      gradient pinned (A_high=%.3f A_low=%.3f margin=%d)\n",
               c->gradient_A_high, c->gradient_A_low, c->gradient_margin);
    printf("Init:    %s", c->init);
    if (!strcmp(c->init, "sfa")) printf(" (%s frame=%d)", c->init_sfa, c->init_frame);
    printf("\nOutput:  %s (%s, snap=%.1f diag=%.1f)\n\n",
           c->output, prec_names[c->precision], c->snap_dt, c->diag_dt);
}

/* ================================================================
   f16 conversion helpers (for SFA output)
   ================================================================ */

static inline uint16_t f64_to_f16(double v) {
    float f = (float)v;
    uint32_t x;
    memcpy(&x, &f, 4);
    uint16_t sign = (x >> 16) & 0x8000;
    int exp = ((x >> 23) & 0xFF) - 127 + 15;
    uint16_t mant = (x >> 13) & 0x3FF;
    if (exp <= 0) return sign;
    if (exp >= 31) return sign | 0x7C00;
    return sign | (exp << 10) | mant;
}

static inline double f16_to_f64(uint16_t h) {
    uint16_t sign = h & 0x8000;
    int exp = (h >> 10) & 0x1F;
    uint16_t mant = h & 0x3FF;
    if (exp == 0) return sign ? -0.0 : 0.0;
    if (exp == 31) return sign ? -INFINITY : INFINITY;
    float f;
    uint32_t x = ((uint32_t)sign << 16) | ((uint32_t)(exp - 15 + 127) << 23) | ((uint32_t)mant << 13);
    memcpy(&f, &x, 4);
    return (double)f;
}

/* ================================================================
   Fixed-point pre-computed constants
   ================================================================ */

typedef struct {
    /* Physics constants in Q16.16 */
    int32_t MASS2;       /* m^2 */
    int32_t MTHETA2;     /* m_theta^2 */
    int32_t ETA;         /* eta */
    int32_t MU;          /* mu (NEGATIVE) */
    int32_t MU_HALF;     /* mu/2 */
    int32_t KAPPA;       /* kappa */
    int32_t THETA_C;     /* theta_c */
    int32_t EPSILON;     /* epsilon */
    int32_t ONE_MINUS_EPS; /* 1 - epsilon */
    int32_t GAMMA;       /* gamma (confinement) */

    /* Grid constants in Q16.16 */
    int32_t DX;          /* dx */
    int32_t DT;          /* dt */
    int32_t HALF_DT;     /* dt/2 */
    int32_t IDX2;        /* 1/dx^2 */
    int32_t IDX1;        /* 1/(2*dx) */

    /* Damping */
    int32_t DAMP_RATE;   /* damp_rate */
    int32_t R_DAMP;      /* L - damp_width */
} FPConst;

static FPConst fp_constants(const Config *c, double dx, double dt) {
    FPConst f;
    f.MASS2       = fp_from_f64(c->m2);
    f.MTHETA2     = fp_from_f64(c->mtheta2);
    f.ETA         = fp_from_f64(c->eta);
    f.MU          = fp_from_f64(c->mu);
    f.MU_HALF     = fp_from_f64(c->mu / 2.0);
    f.KAPPA       = fp_from_f64(c->kappa);
    f.THETA_C     = fp_from_f64(c->theta_c);
    f.EPSILON     = fp_from_f64(c->epsilon);
    f.ONE_MINUS_EPS = fp_from_f64(1.0 - c->epsilon);
    f.GAMMA       = fp_from_f64(c->gamma_conf);
    f.DX          = fp_from_f64(dx);
    f.DT          = fp_from_f64(dt);
    f.HALF_DT     = fp_from_f64(dt / 2.0);
    f.IDX2        = fp_from_f64(1.0 / (dx * dx));
    f.IDX1        = fp_from_f64(1.0 / (2.0 * dx));
    f.DAMP_RATE   = fp_from_f64(c->damp_rate);
    f.R_DAMP      = fp_from_f64(c->L - c->damp_width);

    printf("Fixed-point constants (Q16.16):\n");
    printf("  MASS2=%d (%.4f) MU=%d (%.4f) KAPPA=%d (%.4f)\n",
           f.MASS2, fp_to_f64(f.MASS2), f.MU, fp_to_f64(f.MU), f.KAPPA, fp_to_f64(f.KAPPA));
    printf("  ETA=%d (%.4f) THETA_C=%d (%.6f) EPSILON=%d (%.4f) GAMMA=%d (%.4f)\n",
           f.ETA, fp_to_f64(f.ETA), f.THETA_C, fp_to_f64(f.THETA_C),
           f.EPSILON, fp_to_f64(f.EPSILON), f.GAMMA, fp_to_f64(f.GAMMA));
    printf("  DX=%d (%.4f) DT=%d (%.6f) IDX2=%d (%.4f) IDX1=%d (%.4f)\n",
           f.DX, fp_to_f64(f.DX), f.DT, fp_to_f64(f.DT),
           f.IDX2, fp_to_f64(f.IDX2), f.IDX1, fp_to_f64(f.IDX1));
    printf("\n");
    return f;
}

/* ================================================================
   Grid: 18 int32 arrays (6 fields x {val, vel, acc})
   ================================================================ */

typedef struct {
    int32_t *mem;
    int32_t *phi[NFIELDS], *phi_vel[NFIELDS], *phi_acc[NFIELDS];
    int32_t *theta[NFIELDS], *theta_vel[NFIELDS], *theta_acc[NFIELDS];
    int32_t *pin_phi[NFIELDS], *pin_vel[NFIELDS];
    int32_t *pin_theta[NFIELDS], *pin_tvel[NFIELDS];
    int N; long N3;
    double L, dx, dt;  /* Keep fp64 copies for coordinate calcs and output */
} Grid;

static Grid *grid_alloc(const Config *c) {
    Grid *g = calloc(1, sizeof(Grid));
    g->N  = c->N;
    g->N3 = (long)c->N * c->N * c->N;
    g->L  = c->L;
    g->dx = 2.0 * c->L / (c->N - 1);
    g->dt = c->dt_factor * g->dx;

    long total = 18 * g->N3;
    printf("Allocating %.2f GB (%ld int32s, N=%d, 6 fields, fixed-point)\n",
           total * 4.0 / 1e9, total, c->N);
    g->mem = malloc(total * sizeof(int32_t));
    if (!g->mem) { fprintf(stderr, "FATAL: malloc failed\n"); exit(1); }
    memset(g->mem, 0, total * sizeof(int32_t));

    for (int a = 0; a < NFIELDS; a++) {
        g->phi[a]       = g->mem + (0  + a) * g->N3;
        g->phi_vel[a]   = g->mem + (3  + a) * g->N3;
        g->phi_acc[a]   = g->mem + (6  + a) * g->N3;
        g->theta[a]     = g->mem + (9  + a) * g->N3;
        g->theta_vel[a] = g->mem + (12 + a) * g->N3;
        g->theta_acc[a] = g->mem + (15 + a) * g->N3;
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

static void grid_save_pinned(Grid *g) {
    for (int a = 0; a < NFIELDS; a++) {
        g->pin_phi[a]   = malloc(g->N3 * sizeof(int32_t));
        g->pin_vel[a]   = malloc(g->N3 * sizeof(int32_t));
        g->pin_theta[a] = malloc(g->N3 * sizeof(int32_t));
        g->pin_tvel[a]  = malloc(g->N3 * sizeof(int32_t));
        memcpy(g->pin_phi[a],   g->phi[a],       g->N3 * sizeof(int32_t));
        memcpy(g->pin_vel[a],   g->phi_vel[a],   g->N3 * sizeof(int32_t));
        memcpy(g->pin_theta[a], g->theta[a],     g->N3 * sizeof(int32_t));
        memcpy(g->pin_tvel[a],  g->theta_vel[a], g->N3 * sizeof(int32_t));
    }
}

/* ================================================================
   Curl helper (fixed-point)
   curl(F)_0 = dF2/dy - dF1/dz
   curl(F)_1 = dF0/dz - dF2/dx
   curl(F)_2 = dF1/dx - dF0/dy
   Each derivative is central difference: (F[+1] - F[-1]) / (2*dx)
   In fixed-point: (F[+1] - F[-1]) * IDX1, where IDX1 = fp(1/(2*dx))
   ================================================================ */

static inline int32_t curl_component_fp(int32_t *F[3], int a,
    long n_ip, long n_im, long n_jp, long n_jm, long n_kp, long n_km,
    int32_t idx1) {
    int32_t d1, d2;
    if (a == 0) {
        d1 = F[2][n_jp] - F[2][n_jm];  /* dF2/dy */
        d2 = F[1][n_kp] - F[1][n_km];  /* dF1/dz */
    } else if (a == 1) {
        d1 = F[0][n_kp] - F[0][n_km];  /* dF0/dz */
        d2 = F[2][n_ip] - F[2][n_im];  /* dF2/dx */
    } else {
        d1 = F[1][n_ip] - F[1][n_im];  /* dF1/dx */
        d2 = F[0][n_jp] - F[0][n_jm];  /* dF0/dy */
    }
    /* (d1 - d2) * idx1, using int64 for the multiply */
    return (int32_t)(((int64_t)(d1 - d2) * (int64_t)idx1) >> FP_SHIFT);
}

/* ================================================================
   Initialization — generate in float, convert to fixed-point
   ================================================================ */

static void init_oscillon(Grid *g, const Config *c) {
    const int N = g->N, NN = N * N;
    const double L = g->L, dx = g->dx;
    printf("Init: oscillon (A=%.3f sigma=%.3f)\n", c->A, c->sigma);
    for (int i = 0; i < N; i++) { double x = -L + i * dx;
    for (int j = 0; j < N; j++) { double y = -L + j * dx;
    for (int k = 0; k < N; k++) { double z = -L + k * dx;
        long idx = (long)i * NN + j * N + k;
        double r2 = x * x + y * y + z * z;
        double env = c->A * exp(-r2 / (2.0 * c->sigma * c->sigma));
        for (int a = 0; a < NFIELDS; a++)
            g->phi[a][idx] = fp_from_f64(env * cos(c->delta[a]));
    }}}
}

static void init_braid(Grid *g, const Config *c) {
    const int N = g->N, NN = N * N;
    const double L = g->L, dx = g->dx;
    const double kw = PI / L, omega = sqrt(kw * kw + c->m2);
    const double sx = 1 + c->ellip, sy = 1 - c->ellip;
    const double inv2R2 = 1.0 / (2 * c->R_tube * c->R_tube);
    const double k_bg = PI / L, omega_bg = sqrt(k_bg * k_bg + c->m2);
    printf("Init: braid (R=%.1f ellip=%.4f A=%.2f A_bg=%.2f)\n",
           c->R_tube, c->ellip, c->A, c->A_bg);
    for (int i = 0; i < N; i++) { double x = -L + i * dx;
    for (int j = 0; j < N; j++) { double y = -L + j * dx;
    for (int k = 0; k < N; k++) { double z = -L + k * dx;
        long idx = (long)i * NN + j * N + k;
        double r2e = x * x / (sx * sx) + y * y / (sy * sy);
        double env = exp(-r2e * inv2R2);
        for (int a = 0; a < NFIELDS; a++) {
            double ph = kw * z + c->delta[a];
            double ph_bg = k_bg * z + 2 * PI * a / 3.0;
            g->phi[a][idx]     = fp_from_f64(c->A * env * cos(ph) + c->A_bg * cos(ph_bg));
            g->phi_vel[a][idx] = fp_from_f64(omega * c->A * env * sin(ph) + omega_bg * c->A_bg * sin(ph_bg));
        }
    }}}
}

static void init_from_sfa(Grid *g, const Config *c) {
    printf("Init: SFA file '%s' frame=%d\n", c->init_sfa, c->init_frame);
    SFA *sfa = sfa_open(c->init_sfa);
    if (!sfa) { fprintf(stderr, "FATAL: cannot open SFA '%s'\n", c->init_sfa); exit(1); }
    if ((int)sfa->Nx != g->N || (int)sfa->Ny != g->N || (int)sfa->Nz != g->N) {
        fprintf(stderr, "FATAL: SFA grid %ux%ux%u != sim grid %d^3\n",
                sfa->Nx, sfa->Ny, sfa->Nz, g->N);
        exit(1);
    }
    int frame = c->init_frame;
    if (frame < 0) frame = sfa->total_frames + frame;
    if (frame < 0 || frame >= (int)sfa->total_frames) {
        fprintf(stderr, "FATAL: frame %d out of range [0,%u)\n", frame, sfa->total_frames);
        exit(1);
    }
    printf("  Grid: %ux%ux%u, %d columns, %u frames, reading frame %d\n",
           sfa->Nx, sfa->Ny, sfa->Nz, sfa->n_columns, sfa->total_frames, frame);

    void *buf = malloc(sfa->frame_bytes);
    if (!buf) { fprintf(stderr, "FATAL: frame buffer alloc\n"); exit(1); }
    sfa_read_frame(sfa, frame, buf);

    /* Walk columns: read as float, convert to fixed-point */
    int loaded[12] = {0};
    uint64_t off = 0;
    for (int col = 0; col < (int)sfa->n_columns; col++) {
        int dtype = sfa->columns[col].dtype;
        int sem   = sfa->columns[col].semantic;
        int comp  = sfa->columns[col].component;
        int es    = sfa_dtype_size[dtype];
        uint8_t *src = (uint8_t*)buf + off;

        int32_t *target = NULL;
        int slot = -1;
        if (sem == SFA_POSITION && comp < 3) { target = g->phi[comp]; slot = comp; }
        else if (sem == SFA_ANGLE && comp < 3) { target = g->theta[comp]; slot = 3 + comp; }
        else if (sem == SFA_VELOCITY && comp < 3) { target = g->phi_vel[comp]; slot = 6 + comp; }
        else if (sem == SFA_VELOCITY && comp >= 3 && comp < 6) { target = g->theta_vel[comp - 3]; slot = 9 + comp - 3; }

        if (target) {
            long N3 = g->N3;
            for (long i = 0; i < N3; i++) {
                double val;
                if (dtype == SFA_F64) val = ((double*)src)[i];
                else if (dtype == SFA_F32) val = (double)((float*)src)[i];
                else if (dtype == SFA_F16) val = f16_to_f64(((uint16_t*)src)[i]);
                else val = 0.0;
                target[i] = fp_from_f64(val);
            }
            loaded[slot] = 1;
            printf("  Loaded col %d '%s' (sem=%d comp=%d) -> slot %d (fixed-point)\n",
                   col, sfa->columns[col].name, sem, comp, slot);
        }
        off += (uint64_t)g->N3 * es;
    }

    int n_fields = 0, n_vels = 0;
    for (int i = 0; i < 6; i++) n_fields += loaded[i];
    for (int i = 6; i < 12; i++) n_vels += loaded[i];
    printf("  Loaded: %d/6 field arrays, %d/6 velocity arrays\n", n_fields, n_vels);
    if (n_vels == 0)
        printf("  WARNING: no velocity data -- starting from rest (cold restart)\n");

    free(buf);
    sfa_close(sfa);
}

static void init_from_exec(Grid *g, const Config *c) {
    printf("Init: exec '%s'\n", c->init_exec);
    char tmppath[512];
    snprintf(tmppath, sizeof(tmppath), "/tmp/scp_sim_fixed_init_%d.sfa", (int)getpid());
    char cmd[2048];
    snprintf(cmd, sizeof(cmd), "%s > %s", c->init_exec, tmppath);
    printf("  Running: %s\n", cmd);
    int ret = system(cmd);
    if (ret != 0) { fprintf(stderr, "FATAL: exec init failed (exit %d)\n", ret); exit(1); }
    Config tmp = *c;
    strncpy(tmp.init_sfa, tmppath, 511);
    tmp.init_frame = -1;
    init_from_sfa(g, &tmp);
    remove(tmppath);
}

static void init_template(Grid *g, const Config *c) {
    printf("Init: template '%s'\n", c->init_sfa);
    const int N = g->N, NN = N * N;
    const double L = g->L, dx = g->dx;
    const double k_bg = PI / L, omega_bg = sqrt(k_bg * k_bg + c->m2);

    /* Generate background */
    #pragma omp parallel for schedule(static)
    for (long idx = 0; idx < g->N3; idx++) {
        int i = (int)(idx / NN), j = (int)((idx / N) % N), k = (int)(idx % N);
        double x = -L + i * dx, z = -L + k * dx;
        double A_bg_x = c->A_bg;
        if (c->bc_type == 1) {
            double frac = (x + L) / (2.0 * L);
            A_bg_x = c->gradient_A_high * (1.0 - frac) + c->gradient_A_low * frac;
        }
        for (int a = 0; a < NFIELDS; a++) {
            double ph_bg = k_bg * z + 2 * PI * a / 3.0;
            g->phi[a][idx] = fp_from_f64(A_bg_x * cos(ph_bg));
            g->phi_vel[a][idx] = fp_from_f64(omega_bg * A_bg_x * sin(ph_bg));
        }
    }

    /* Load template and stamp at center */
    SFA *tmpl = sfa_open(c->init_sfa);
    if (!tmpl) { fprintf(stderr, "FATAL: cannot open template '%s'\n", c->init_sfa); exit(1); }
    int TN = tmpl->Nx; double TL = tmpl->Lx, Tdx = 2.0 * TL / (TN - 1);
    long TN3 = (long)TN * TN * TN; int TNN = TN * TN;
    printf("  Template: %d^3, L=%.2f, dx=%.4f\n", TN, TL, Tdx);

    int frame = tmpl->total_frames - 1;
    void *buf = malloc(tmpl->frame_bytes);
    sfa_read_frame(tmpl, frame, buf);

    /* Read template into temporary float arrays, then convert */
    double *tphi[3], *tvel[3], *ttheta[3], *ttvel[3];
    for (int a = 0; a < 3; a++) {
        tphi[a] = calloc(TN3, sizeof(double)); tvel[a] = calloc(TN3, sizeof(double));
        ttheta[a] = calloc(TN3, sizeof(double)); ttvel[a] = calloc(TN3, sizeof(double));
    }
    uint64_t off = 0;
    for (int col = 0; col < (int)tmpl->n_columns; col++) {
        int dtype = tmpl->columns[col].dtype, sem = tmpl->columns[col].semantic, comp = tmpl->columns[col].component;
        int es = sfa_dtype_size[dtype]; uint8_t *src = (uint8_t*)buf + off;
        double *target = NULL;
        if (sem == SFA_POSITION && comp < 3) target = tphi[comp];
        else if (sem == SFA_ANGLE && comp < 3) target = ttheta[comp];
        else if (sem == SFA_VELOCITY && comp < 3) target = tvel[comp];
        else if (sem == SFA_VELOCITY && comp >= 3 && comp < 6) target = ttvel[comp - 3];
        if (target) {
            if (dtype == SFA_F64) for (long i = 0; i < TN3; i++) target[i] = ((double*)src)[i];
            else if (dtype == SFA_F32) for (long i = 0; i < TN3; i++) target[i] = (double)((float*)src)[i];
            else if (dtype == SFA_F16) for (long i = 0; i < TN3; i++) target[i] = f16_to_f64(((uint16_t*)src)[i]);
        }
        off += (uint64_t)TN3 * es;
    }
    free(buf); sfa_close(tmpl);

    int ci = N / 2, cj = N / 2, ck = N / 2, half = TN / 2;
    int placed = 0;
    for (int ti = 0; ti < TN; ti++) {
        int gi = ci + ti - half; if (gi < 0 || gi >= N) continue;
        for (int tj = 0; tj < TN; tj++) {
            int gj = cj + tj - half; if (gj < 0 || gj >= N) continue;
            for (int tk = 0; tk < TN; tk++) {
                int gk = ck + tk - half; if (gk < 0 || gk >= N) continue;
                long tidx = (long)ti * TNN + tj * TN + tk;
                long gidx = (long)gi * NN + gj * N + gk;
                double tz = -TL + tk * Tdx;
                for (int a = 0; a < NFIELDS; a++) {
                    double ph_bg_t = k_bg * tz + 2 * PI * a / 3.0;
                    double bg_phi = c->A_bg * cos(ph_bg_t);
                    double bg_vel = omega_bg * c->A_bg * sin(ph_bg_t);
                    double new_phi = fp_to_f64(g->phi[a][gidx]) + tphi[a][tidx] - bg_phi;
                    double new_vel = fp_to_f64(g->phi_vel[a][gidx]) + tvel[a][tidx] - bg_vel;
                    g->phi[a][gidx] = fp_from_f64(new_phi);
                    g->phi_vel[a][gidx] = fp_from_f64(new_vel);
                    g->theta[a][gidx] = fp_from_f64(fp_to_f64(g->theta[a][gidx]) + ttheta[a][tidx]);
                    g->theta_vel[a][gidx] = fp_from_f64(fp_to_f64(g->theta_vel[a][gidx]) + ttvel[a][tidx]);
                }
                placed++;
            }
        }
    }
    for (int a = 0; a < 3; a++) { free(tphi[a]); free(tvel[a]); free(ttheta[a]); free(ttvel[a]); }
    printf("  Placed %d voxels from template\n", placed);
}

static void do_init(Grid *g, const Config *c) {
    if      (!strcmp(c->init, "oscillon")) init_oscillon(g, c);
    else if (!strcmp(c->init, "braid"))    init_braid(g, c);
    else if (!strcmp(c->init, "sfa"))      init_from_sfa(g, c);
    else if (!strcmp(c->init, "exec"))     init_from_exec(g, c);
    else if (!strcmp(c->init, "template")) init_template(g, c);
    else { fprintf(stderr, "FATAL: unknown init mode '%s'\n", c->init); exit(1); }
}

/* ================================================================
   Forces: 6-field Cosserat with Variant B v2 unified transfer
   ALL computation in int32/int64 fixed-point.
   ================================================================ */

static void compute_forces(Grid *g, const FPConst *fp, const Config *c) {
    const int N = g->N, NN = N * N;
    const long N3 = g->N3;
    const int32_t IDX2 = fp->IDX2;
    const int32_t IDX1 = fp->IDX1;
    const int32_t MASS2 = fp->MASS2;
    const int32_t MTHETA2 = fp->MTHETA2;
    const int32_t ETA = fp->ETA;
    const int32_t MU = fp->MU;
    const int32_t MU_HALF = fp->MU_HALF;
    const int32_t KAPPA = fp->KAPPA;
    const int32_t THETA_C = fp->THETA_C;
    const int32_t EPSILON = fp->EPSILON;
    const int32_t ONE_MINUS_EPS = fp->ONE_MINUS_EPS;
    const int32_t GAMMA = fp->GAMMA;

    #pragma omp parallel for schedule(static)
    for (long idx = 0; idx < N3; idx++) {
        int i = (int)(idx / NN), j = (int)((idx / N) % N), k = (int)(idx % N);
        int ip = (i + 1) % N, im = (i - 1 + N) % N;
        int jp = (j + 1) % N, jm = (j - 1 + N) % N;
        int kp = (k + 1) % N, km = (k - 1 + N) % N;
        long n_ip = (long)ip * NN + j * N + k, n_im = (long)im * NN + j * N + k;
        long n_jp = (long)i * NN + jp * N + k, n_jm = (long)i * NN + jm * N + k;
        long n_kp = (long)i * NN + j * N + kp, n_km = (long)i * NN + j * N + km;

        int32_t p0 = g->phi[0][idx], p1 = g->phi[1][idx], p2 = g->phi[2][idx];

        /* Triple product P = phi[0]*phi[1]*phi[2] in Q16.16
         *
         * TODO: PRECISION — At background phi=0.1 (int ~6554), the chained
         * fp_mul gives P~65 (representing 0.001). Then P² = fp_sq(65) =
         * 4225>>16 = 0 — underflows to exactly zero. This means V_base,
         * dVdP, and all potential forces are hard-zero in the background.
         * Physically acceptable (background P should be negligible) but
         * creates a hard step rather than smooth rolloff at the braid edge.
         *
         * At braid peak (phi~0.3, int ~19661): P~1769 (0.027), P²~47
         * (0.0007) — only ~6 bits of precision. Marginal for force accuracy.
         *
         * For production: compute P and P² in a wider intermediate format
         * (e.g., Q8.24 for the triple product chain, then convert back to
         * Q16.16 for the force), or use int64 for the P² computation:
         *   int64_t P2_wide = ((int64_t)P * P);  // Q32.32
         *   then scale kappa*P2 and den in int64 before truncating.
         */
        int32_t P = fp_mul(fp_mul(p0, p1), p2);
        int32_t P2 = fp_sq(P);

        /* V_base = (mu/2) * P^2 / (1 + kappa*P^2)
         * den = FP_ONE + kappa*P^2 */
        int32_t kP2 = fp_mul(KAPPA, P2);
        int32_t den = FP_ONE + kP2;

        /* V_base = MU_HALF * P2 / den  (NEGATIVE, since mu < 0) */
        int32_t V_base = fp_mul(MU_HALF, fp_div(P2, den));
        int32_t V_abs = fp_abs(V_base);  /* |V_base|, POSITIVE */

        /* dV_base/dP = mu * P / den^2 */
        int32_t den2 = fp_mul(den, den);
        int32_t dVdP = fp_div(fp_mul(MU, P), den2);

        /* Theta magnitude: Theta = t0^2 + t1^2 + t2^2 */
        int32_t t0 = g->theta[0][idx], t1 = g->theta[1][idx], t2 = g->theta[2][idx];
        int32_t Theta = fp_sq(t0) + fp_sq(t1) + fp_sq(t2);

        /* Transfer function: f = epsilon + (1-epsilon) * Theta / (Theta + Theta_c) */
        int32_t Theta_sum = Theta + THETA_C;
        int32_t f_transfer;
        if (Theta_sum > 0) {
            f_transfer = EPSILON + fp_mul(ONE_MINUS_EPS, fp_div(Theta, Theta_sum));
        } else {
            f_transfer = EPSILON;
        }

        /* Confinement: f_confine = gamma * ln(1 + Theta/Theta_c)
         * = gamma * (ln(Theta + Theta_c) - ln(Theta_c)) */
        int32_t f_confine;
        if (Theta > 0 && GAMMA > 0 && THETA_C > 0) {
            f_confine = fp_mul(GAMMA, fp_ln1p_ratio(Theta, THETA_C));
        } else {
            f_confine = 0;
        }

        /* Phi force modulation: transfer ATTRACTS, confine REPELS.
         * W = V_base * f_transfer  →  dW/dP = dVdP * f_transfer
         * W_confine = |V_base| * f_confine  →  dWc/dP = -dVdP * f_confine
         * Total: dU/dP = dVdP * (f_transfer - f_confine)
         * Force: -dVdP * (f_transfer - f_confine) * dPda
         * The MINUS on f_confine provides negative feedback. */
        int32_t f_phi = f_transfer - f_confine;

        /* Theta force: F = +V_abs × (df/dΘ - confine_deriv) × 2 × theta_a
         * Small Θ: df > confine → DRIVES theta growth (transfer)
         * Large Θ: confine > df → CONFINES theta (pulls back)
         * Equilibrium where df_dTheta = confine_deriv. */
        int32_t Theta_sum2 = fp_mul(Theta_sum, Theta_sum);
        int32_t df_dTheta = 0;
        int32_t confine_deriv = 0;
        if (Theta_sum > 0 && Theta_sum2 > 0) {
            df_dTheta = fp_div(fp_mul(ONE_MINUS_EPS, THETA_C), Theta_sum2);
            confine_deriv = fp_div(GAMMA, Theta_sum);
        }
        /* drive_coeff = (df_dTheta - confine_deriv) * 2
         * Positive at small Θ (drives growth), negative at large Θ (confines) */
        int32_t drive_coeff = (df_dTheta - confine_deriv) << 1;
        /* theta_drive = V_abs * drive_coeff */
        int32_t theta_drive = fp_mul(V_abs, drive_coeff);

        /* --- Phi forces --- */
        for (int a = 0; a < NFIELDS; a++) {
            /* Laplacian: (sum_neighbors - 6*center) * idx2 */
            int32_t center = g->phi[a][idx];
            int64_t lap64 = (int64_t)g->phi[a][n_ip] + g->phi[a][n_im]
                          + g->phi[a][n_jp] + g->phi[a][n_jm]
                          + g->phi[a][n_kp] + g->phi[a][n_km]
                          - 6 * (int64_t)center;
            /* lap = lap64 * IDX2, but lap64 is in Q16.16 and IDX2 is in Q16.16
             * so the product is Q32.32, shift right by 16 to get Q16.16 */
            int32_t lap = (int32_t)((lap64 * (int64_t)IDX2) >> FP_SHIFT);

            /* dP/dphi_a */
            int32_t dPda;
            if      (a == 0) dPda = fp_mul(p1, p2);
            else if (a == 1) dPda = fp_mul(p0, p2);
            else             dPda = fp_mul(p0, p1);

            /* Potential force: -dVdP * dPda * f_phi */
            int32_t pot = fp_mul(fp_mul(dVdP, dPda), f_phi);
            /* Note: pot = dVdP*dPda*f_phi. Since dVdP < 0 (mu<0):
             * -pot = -dVdP*dPda*f_phi = |dVdP|*dPda*f_phi
             * If f_phi > 0: attractive. If f_phi < 0: repulsive (confine dominates). */

            /* Mass term: -m^2 * phi */
            int32_t mass_term = fp_mul(MASS2, center);

            /* Curl coupling: eta * curl(theta)_a */
            int32_t ct = curl_component_fp(g->theta, a, n_ip, n_im, n_jp, n_jm, n_kp, n_km, IDX1);
            int32_t curl_term = fp_mul(ETA, ct);

            /* Total: lap - m^2*phi - dVdP*dPda*f_phi + eta*curl(theta)
             * = lap - mass_term - pot + curl_term */
            g->phi_acc[a][idx] = lap - mass_term - pot + curl_term;
        }

        /* --- Theta forces --- */
        for (int a = 0; a < NFIELDS; a++) {
            /* Laplacian */
            int32_t center_t = g->theta[a][idx];
            int64_t lapt64 = (int64_t)g->theta[a][n_ip] + g->theta[a][n_im]
                           + g->theta[a][n_jp] + g->theta[a][n_jm]
                           + g->theta[a][n_kp] + g->theta[a][n_km]
                           - 6 * (int64_t)center_t;
            int32_t lapt = (int32_t)((lapt64 * (int64_t)IDX2) >> FP_SHIFT);

            /* Curl coupling: eta * curl(phi)_a */
            int32_t cp = curl_component_fp(g->phi, a, n_ip, n_im, n_jp, n_jm, n_kp, n_km, IDX1);
            int32_t curl_p = fp_mul(ETA, cp);

            /* Unified transfer force: +theta_drive * theta (drives or confines)
             * Plus bare theta mass: -MTHETA2 * theta (always confining) */
            int32_t drive_t = fp_mul(theta_drive, center_t);
            int32_t bare_mass_t = fp_mul(MTHETA2, center_t);

            /* Total: lap + theta_drive*theta - mtheta2*theta + eta*curl(phi) */
            g->theta_acc[a][idx] = lapt + drive_t - bare_mass_t + curl_p;
        }
    }
}

/* ================================================================
   Boundary conditions (in fixed-point)
   ================================================================ */

/* bc_type=0: Spherical absorbing boundary */
static void apply_damping(Grid *g, const FPConst *fp, const Config *c) {
    const int N = g->N, NN = N * N;
    const double L = g->L, dx = g->dx, DW = c->damp_width;
    if (DW <= 0 || c->damp_rate <= 0) return;
    const double R_damp = L - DW;

    #pragma omp parallel for schedule(static)
    for (long idx = 0; idx < g->N3; idx++) {
        int i = (int)(idx / NN), j = (int)((idx / N) % N), k = (int)(idx % N);
        double x = -L + i * dx, y = -L + j * dx, z = -L + k * dx;
        double r = sqrt(x * x + y * y + z * z);
        if (r > R_damp) {
            double s = (r - R_damp) / DW; if (s > 1) s = 1;
            /* Damping factor: d = 1 - rate*s^2 */
            int32_t d = fp_from_f64(1.0 - c->damp_rate * s * s);
            for (int a = 0; a < NFIELDS; a++) {
                g->phi_vel[a][idx] = fp_mul(g->phi_vel[a][idx], d);
                g->theta_vel[a][idx] = fp_mul(g->theta_vel[a][idx], d);
            }
        }
    }
}

/* bc_type=1: Gradient pinned boundary */
static void apply_gradient_bc(Grid *g, const Config *c) {
    const int N = g->N, NN = N * N;
    const int M = c->gradient_margin;

    /* x-direction: pin boundary slabs */
    #pragma omp parallel for schedule(static)
    for (long idx = 0; idx < g->N3; idx++) {
        int i = (int)(idx / NN);
        if (i < M || i >= N - M) {
            for (int a = 0; a < NFIELDS; a++) {
                g->phi[a][idx]       = g->pin_phi[a][idx];
                g->phi_vel[a][idx]   = g->pin_vel[a][idx];
                g->phi_acc[a][idx]   = 0;
                g->theta[a][idx]     = g->pin_theta[a][idx];
                g->theta_vel[a][idx] = g->pin_tvel[a][idx];
                g->theta_acc[a][idx] = 0;
            }
        }
    }

    /* y-direction: linear extrapolation */
    #pragma omp parallel for schedule(static)
    for (int i = M; i < N - M; i++) {
        for (int k = 0; k < N; k++) {
            for (int a = 0; a < NFIELDS; a++) {
                long idx0 = (long)i * NN + 0 * N + k;
                long idx1 = (long)i * NN + 1 * N + k;
                long idx2 = (long)i * NN + 2 * N + k;
                g->phi[a][idx0] = 2 * g->phi[a][idx1] - g->phi[a][idx2];
                g->phi_vel[a][idx0] = 2 * g->phi_vel[a][idx1] - g->phi_vel[a][idx2];
                g->theta[a][idx0] = 2 * g->theta[a][idx1] - g->theta[a][idx2];
                g->theta_vel[a][idx0] = 2 * g->theta_vel[a][idx1] - g->theta_vel[a][idx2];

                idx0 = (long)i * NN + (N - 1) * N + k;
                idx1 = (long)i * NN + (N - 2) * N + k;
                idx2 = (long)i * NN + (N - 3) * N + k;
                g->phi[a][idx0] = 2 * g->phi[a][idx1] - g->phi[a][idx2];
                g->phi_vel[a][idx0] = 2 * g->phi_vel[a][idx1] - g->phi_vel[a][idx2];
                g->theta[a][idx0] = 2 * g->theta[a][idx1] - g->theta[a][idx2];
                g->theta_vel[a][idx0] = 2 * g->theta_vel[a][idx1] - g->theta_vel[a][idx2];
            }
        }
    }
}

/* ================================================================
   Verlet integrator (fixed-point)
   ================================================================ */

static void verlet_step(Grid *g, const FPConst *fp, const Config *c) {
    const long N3 = g->N3;
    const int32_t HALF_DT = fp->HALF_DT;
    const int32_t DT = fp->DT;

    /* Half-kick: v += (dt/2) * a */
    for (int a = 0; a < NFIELDS; a++) {
        int32_t *vp = g->phi_vel[a], *ap = g->phi_acc[a];
        int32_t *vt = g->theta_vel[a], *at = g->theta_acc[a];
        #pragma omp parallel for schedule(static)
        for (long i = 0; i < N3; i++) {
            vp[i] += fp_mul(HALF_DT, ap[i]);
            vt[i] += fp_mul(HALF_DT, at[i]);
        }
    }

    /* Drift: x += dt * v */
    for (int a = 0; a < NFIELDS; a++) {
        int32_t *pp = g->phi[a], *vp = g->phi_vel[a];
        int32_t *pt = g->theta[a], *vt = g->theta_vel[a];
        #pragma omp parallel for schedule(static)
        for (long i = 0; i < N3; i++) {
            pp[i] += fp_mul(DT, vp[i]);
            pt[i] += fp_mul(DT, vt[i]);
        }
    }

    /* Compute new forces */
    compute_forces(g, fp, c);

    /* Half-kick again: v += (dt/2) * a */
    for (int a = 0; a < NFIELDS; a++) {
        int32_t *vp = g->phi_vel[a], *ap = g->phi_acc[a];
        int32_t *vt = g->theta_vel[a], *at = g->theta_acc[a];
        #pragma omp parallel for schedule(static)
        for (long i = 0; i < N3; i++) {
            vp[i] += fp_mul(HALF_DT, ap[i]);
            vt[i] += fp_mul(HALF_DT, at[i]);
        }
    }

    /* Boundary conditions */
    if (c->bc_type == 0)
        apply_damping(g, fp, c);
    else if (c->bc_type == 1)
        apply_gradient_bc(g, c);
}

/* ================================================================
   Diagnostics — compute energies using int64 accumulation,
   convert to double for output
   ================================================================ */

static void compute_energy(Grid *g, const FPConst *fp, const Config *c,
    double *epk, double *etk, double *eg, double *em, double *ep,
    double *etg, double *etm, double *ec, double *et,
    double *phi_max_out, double *P_max_out,
    double *e_transfer_out, double *e_confine_out, double *f_avg_out) {

    const int N = g->N, NN = N * N;
    const long N3 = g->N3;
    const double dx = g->dx, dV = dx * dx * dx;
    const double idx1_f = 1.0 / (2.0 * dx);

    /* We accumulate in double since energy is a diagnostic (not physics).
     * The physics loop is all fixed-point; energy is computed at diagnostic
     * intervals only, so float here is acceptable. */
    double s_epk = 0, s_etk = 0, s_eg = 0, s_em = 0, s_ep = 0;
    double s_etg = 0, s_etm = 0, s_ec = 0;
    double s_pm = 0, s_Pm = 0;
    double s_etr = 0, s_ecf = 0, s_favg = 0;

    #pragma omp parallel for reduction(+:s_epk,s_etk,s_eg,s_em,s_ep,s_etg,s_etm,s_ec,s_etr,s_ecf,s_favg) \
        reduction(max:s_pm,s_Pm) schedule(static)
    for (long idx = 0; idx < N3; idx++) {
        int i = (int)(idx / NN), j = (int)((idx / N) % N), k = (int)(idx % N);
        int ip = (i + 1) % N, im = (i - 1 + N) % N;
        int jp = (j + 1) % N, jm = (j - 1 + N) % N;
        int kp = (k + 1) % N, km = (k - 1 + N) % N;
        long n_ip = (long)ip * NN + j * N + k, n_im = (long)im * NN + j * N + k;
        long n_jp = (long)i * NN + jp * N + k, n_jm = (long)i * NN + jm * N + k;
        long n_kp = (long)i * NN + j * N + kp, n_km = (long)i * NN + j * N + km;

        /* Convert fields to double for energy computation */
        double p0 = fp_to_f64(g->phi[0][idx]);
        double p1 = fp_to_f64(g->phi[1][idx]);
        double p2 = fp_to_f64(g->phi[2][idx]);

        double P = p0 * p1 * p2, P2 = P * P;
        double KAPPA_f = fp_to_f64(fp->KAPPA);
        double MU_f = fp_to_f64(fp->MU);
        double den = 1.0 + KAPPA_f * P2;

        double V_base = (MU_f / 2.0) * P2 / den;
        double V_abs = fabs(V_base);

        double tt0 = fp_to_f64(g->theta[0][idx]);
        double tt1 = fp_to_f64(g->theta[1][idx]);
        double tt2 = fp_to_f64(g->theta[2][idx]);
        double Theta = tt0 * tt0 + tt1 * tt1 + tt2 * tt2;

        double THETA_C_f = fp_to_f64(fp->THETA_C);
        double EPS_f = fp_to_f64(fp->EPSILON);
        double GAMMA_f = fp_to_f64(fp->GAMMA);
        double Theta_sum = Theta + THETA_C_f;

        double f_transfer = EPS_f + (1.0 - EPS_f) * Theta / (Theta_sum > 0 ? Theta_sum : 1e-30);
        double f_confine = (GAMMA_f > 0 && THETA_C_f > 0) ? GAMMA_f * log(1.0 + Theta / THETA_C_f) : 0;

        s_etr += V_base * f_transfer * dV;
        s_ecf += V_abs * f_confine * dV;
        s_favg += f_transfer;

        for (int a = 0; a < NFIELDS; a++) {
            double pv = fp_to_f64(g->phi_vel[a][idx]);
            double tv = fp_to_f64(g->theta_vel[a][idx]);
            s_epk += 0.5 * pv * pv * dV;
            s_etk += 0.5 * tv * tv * dV;

            double phi_here = fp_to_f64(g->phi[a][idx]);
            double gx = (fp_to_f64(g->phi[a][n_ip]) - fp_to_f64(g->phi[a][n_im])) * idx1_f;
            double gy = (fp_to_f64(g->phi[a][n_jp]) - fp_to_f64(g->phi[a][n_jm])) * idx1_f;
            double gz = (fp_to_f64(g->phi[a][n_kp]) - fp_to_f64(g->phi[a][n_km])) * idx1_f;
            s_eg += 0.5 * (gx * gx + gy * gy + gz * gz) * dV;
            s_em += 0.5 * fp_to_f64(fp->MASS2) * phi_here * phi_here * dV;

            double theta_here = fp_to_f64(g->theta[a][idx]);
            double tgx = (fp_to_f64(g->theta[a][n_ip]) - fp_to_f64(g->theta[a][n_im])) * idx1_f;
            double tgy = (fp_to_f64(g->theta[a][n_jp]) - fp_to_f64(g->theta[a][n_jm])) * idx1_f;
            double tgz = (fp_to_f64(g->theta[a][n_kp]) - fp_to_f64(g->theta[a][n_km])) * idx1_f;
            s_etg += 0.5 * (tgx * tgx + tgy * tgy + tgz * tgz) * dV;

            /* Theta mass from transfer + confinement */
            double theta_coeff_f = V_abs * ((1.0 - EPS_f) * THETA_C_f / (Theta_sum * Theta_sum) + GAMMA_f / Theta_sum) * 2.0;
            s_etm += 0.5 * (fp_to_f64(fp->MTHETA2) + theta_coeff_f) * theta_here * theta_here * dV;

            double ap = fabs(phi_here); if (ap > s_pm) s_pm = ap;
        }

        s_ep += V_base * dV;  /* Raw potential (without transfer modulation, for reference) */
        double Pa = fabs(P); if (Pa > s_Pm) s_Pm = Pa;

        /* Curl coupling energy */
        double eta_f = fp_to_f64(fp->ETA);
        for (int a = 0; a < NFIELDS; a++) {
            double phi_here = fp_to_f64(g->phi[a][idx]);
            double d1, d2;
            if (a == 0) {
                d1 = fp_to_f64(g->theta[2][n_jp]) - fp_to_f64(g->theta[2][n_jm]);
                d2 = fp_to_f64(g->theta[1][n_kp]) - fp_to_f64(g->theta[1][n_km]);
            } else if (a == 1) {
                d1 = fp_to_f64(g->theta[0][n_kp]) - fp_to_f64(g->theta[0][n_km]);
                d2 = fp_to_f64(g->theta[2][n_ip]) - fp_to_f64(g->theta[2][n_im]);
            } else {
                d1 = fp_to_f64(g->theta[1][n_ip]) - fp_to_f64(g->theta[1][n_im]);
                d2 = fp_to_f64(g->theta[0][n_jp]) - fp_to_f64(g->theta[0][n_jm]);
            }
            double ct = (d1 - d2) * idx1_f;
            s_ec -= eta_f * phi_here * ct * dV;
        }
    }

    *epk = s_epk; *etk = s_etk; *eg = s_eg; *em = s_em; *ep = s_ep;
    *etg = s_etg; *etm = s_etm; *ec = s_ec;
    *et = s_epk + s_etk + s_eg + s_em + s_etr + s_ecf + s_etg + s_etm + s_ec;
    *phi_max_out = s_pm; *P_max_out = s_Pm;
    *e_transfer_out = s_etr;
    *e_confine_out = s_ecf;
    *f_avg_out = s_favg / N3;
}

static double theta_rms(Grid *g) {
    double sum = 0;
    #pragma omp parallel for reduction(+:sum)
    for (long i = 0; i < g->N3; i++)
        for (int a = 0; a < NFIELDS; a++) {
            double v = fp_to_f64(g->theta[a][i]);
            sum += v * v;
        }
    return sqrt(sum / (3.0 * g->N3));
}

static double P_integrated(Grid *g) {
    double t = 0;
    const double dV = g->dx * g->dx * g->dx;
    #pragma omp parallel for reduction(+:t)
    for (long i = 0; i < g->N3; i++) {
        double p = fp_to_f64(g->phi[0][i]) * fp_to_f64(g->phi[1][i]) * fp_to_f64(g->phi[2][i]);
        t += fabs(p) * dV;
    }
    return t;
}

/* ================================================================
   SFA output: convert int32 fixed-point to float for SFA
   ================================================================ */

static void sfa_snap(SFA *sfa, Grid *g, double t, int precision) {
    long n = g->N3;

    /* Convert all 12 channels from fixed-point to float */
    int32_t *src_arrays[12] = {
        g->phi[0], g->phi[1], g->phi[2],
        g->theta[0], g->theta[1], g->theta[2],
        g->phi_vel[0], g->phi_vel[1], g->phi_vel[2],
        g->theta_vel[0], g->theta_vel[1], g->theta_vel[2]
    };

    if (precision == 2) {
        /* f64 output */
        double *cols_f64[12];
        for (int c = 0; c < 12; c++) {
            cols_f64[c] = malloc(n * sizeof(double));
            for (long i = 0; i < n; i++)
                cols_f64[c][i] = fp_to_f64(src_arrays[c][i]);
        }
        sfa_write_frame(sfa, t, (void**)cols_f64);
        for (int c = 0; c < 12; c++) free(cols_f64[c]);
    } else if (precision == 1) {
        /* f32 output */
        float *cols_f32[12];
        for (int c = 0; c < 12; c++) {
            cols_f32[c] = malloc(n * sizeof(float));
            for (long i = 0; i < n; i++)
                cols_f32[c][i] = (float)fp_to_f64(src_arrays[c][i]);
        }
        sfa_write_frame(sfa, t, (void**)cols_f32);
        for (int c = 0; c < 12; c++) free(cols_f32[c]);
    } else {
        /* f16 output */
        uint16_t *cols_f16[12];
        for (int c = 0; c < 12; c++) {
            cols_f16[c] = malloc(n * sizeof(uint16_t));
            for (long i = 0; i < n; i++)
                cols_f16[c][i] = f64_to_f16(fp_to_f64(src_arrays[c][i]));
        }
        sfa_write_frame(sfa, t, (void**)cols_f16);
        for (int c = 0; c < 12; c++) free(cols_f16[c]);
    }
}

/* ================================================================
   Main
   ================================================================ */

int main(int argc, char **argv) {
    Config c = cfg_defaults();

    if (argc < 2) {
        fprintf(stderr, "Usage: %s config.cfg [-key value ...]\n", argv[0]);
        fprintf(stderr, "       %s input.sfa [-key value ...]   (restart from SFA)\n", argv[0]);
        return 1;
    }

    int arg_start = 2;
    const char *arg1 = argv[1];
    int len1 = strlen(arg1);

    if (len1 > 4 && !strcmp(arg1 + len1 - 4, ".sfa")) {
        /* First arg is an SFA file: read KVMD as base config */
        printf("Loading parameters from SFA: %s\n", arg1);
        SFA *init_sfa = sfa_open(arg1);
        if (!init_sfa) { fprintf(stderr, "Cannot open SFA: %s\n", arg1); return 1; }

        SFA_KVMDSet kv[SFA_MAX_KVMD_SETS];
        int n_kv = sfa_read_kvmd(init_sfa, kv, SFA_MAX_KVMD_SETS);
        if (n_kv == 0) {
            printf("  WARNING: no KVMD metadata found -- using defaults\n");
        } else {
            SFA_KVMDSet *use = &kv[0];
            for (int i = 0; i < n_kv; i++) {
                if (kv[i].first_frame == 0xFFFFFFFF) { use = &kv[i]; break; }
            }
            printf("  KVMD set %d: %d parameters\n", use->set_id, use->n_pairs);
            for (int i = 0; i < use->n_pairs; i++)
                cfg_set(&c, use->keys[i], use->values[i]);
        }

        char tmp[64];
        snprintf(tmp, 64, "%u", init_sfa->Nx); cfg_set(&c, "N", tmp);
        snprintf(tmp, 64, "%.6f", init_sfa->Lx); cfg_set(&c, "L", tmp);

        cfg_set(&c, "init", "sfa");
        strncpy(c.init_sfa, arg1, 511);
        c.init_frame = -1;

        sfa_close(init_sfa);
    } else {
        cfg_load(&c, arg1);
    }

    /* CLI overrides */
    for (int i = arg_start; i < argc - 1; i += 2) {
        const char *key = argv[i];
        if (key[0] == '-') key++;
        cfg_set(&c, key, argv[i + 1]);
    }

    /* OpenMP setup */
    int nth = 4;
    char *env_t = getenv("OMP_NUM_THREADS");
    if (env_t) nth = atoi(env_t);
    omp_set_num_threads(nth);

    cfg_print(&c);

    /* Grid */
    Grid *g = grid_alloc(&c);
    printf("dx=%.4f dt=%.6f threads=%d\n", g->dx, g->dt, nth);

    /* Fixed-point constants */
    FPConst fp = fp_constants(&c, g->dx, g->dt);

    /* Verify fixed-point ranges */
    printf("Range check:\n");
    printf("  mu=%.3f -> FP=%d (%.3f) %s\n",
           c.mu, fp.MU, fp_to_f64(fp.MU), (fp_abs(fp.MU) > FP_ONE * 100) ? "WARNING: large" : "OK");
    printf("  kappa=%.1f -> FP=%d (%.3f) %s\n",
           c.kappa, fp.KAPPA, fp_to_f64(fp.KAPPA), (fp.KAPPA > FP_ONE * 100) ? "WARNING: large" : "OK");
    printf("  1/dx^2=%.3f -> FP=%d (%.3f) %s\n",
           1.0/(g->dx*g->dx), fp.IDX2, fp_to_f64(fp.IDX2),
           (fp.IDX2 > FP_ONE * 1000) ? "WARNING: large" : "OK");
    printf("\n");

    /* Initialize */
    do_init(g, &c);
    compute_forces(g, &fp, &c);

    /* For gradient_pinned BC: save initial state as pinned boundary values */
    if (c.bc_type == 1) {
        grid_save_pinned(g);
        printf("Gradient BC: pinned %d slabs on each x-face\n\n", c.gradient_margin);
    }

    /* SFA archive */
    uint8_t sfa_dtype = (c.precision == 0) ? SFA_F16 : (c.precision == 1) ? SFA_F32 : SFA_F64;
    SFA *sfa = sfa_create(c.output, c.N, c.N, c.N, c.L, c.L, c.L, g->dt);
    sfa->flags = SFA_CODEC_COLZSTD | SFA_FLAG_STREAMING;
    if (c.output_split) {
        sfa_enable_split(sfa);
        strncpy(sfa->split_base, c.output, 507);
        char *ext = strrchr(sfa->split_base, '.');
        if (ext && !strcmp(ext, ".sfa")) *ext = '\0';
    }

    /* Embed physics parameters as KVMD metadata */
    {
        char vN[32], vL[32], vT[32], vdt[32], vm[32], vmt[32], veta[32], vmu[32], vkappa[32];
        char vtc[32], veps[32], vgam[32], vdw[32], vdr[32], vprec[32], vdelta[64], vbcsw[32];
        snprintf(vN, 32, "%d", c.N); snprintf(vL, 32, "%.6f", c.L); snprintf(vT, 32, "%.6f", c.T);
        snprintf(vdt, 32, "%.6f", c.dt_factor);
        snprintf(vm, 32, "%.6f", sqrt(c.m2)); snprintf(vmt, 32, "%.6f", sqrt(c.mtheta2));
        snprintf(veta, 32, "%.6f", c.eta);
        snprintf(vmu, 32, "%.6f", c.mu); snprintf(vkappa, 32, "%.6f", c.kappa);
        snprintf(vtc, 32, "%.6f", c.theta_c); snprintf(veps, 32, "%.6f", c.epsilon);
        snprintf(vgam, 32, "%.6f", c.gamma_conf);
        snprintf(vdw, 32, "%.6f", c.damp_width); snprintf(vdr, 32, "%.6f", c.damp_rate);
        snprintf(vprec, 32, "%s", (const char*[]){"f16", "f32", "f64"}[c.precision]);
        snprintf(vdelta, 64, "%.6f,%.6f,%.6f", c.delta[0], c.delta[1], c.delta[2]);
        snprintf(vbcsw, 32, "%.6f", c.bc_switch_time);
        char vbc[32], vgah[32], vgal[32], vgm[32];
        snprintf(vbc, 32, "%d", c.bc_type);
        snprintf(vgah, 32, "%.6f", c.gradient_A_high);
        snprintf(vgal, 32, "%.6f", c.gradient_A_low);
        snprintf(vgm, 32, "%d", c.gradient_margin);
        const char *keys[] = {"N", "L", "T", "dt_factor", "m", "m_theta", "eta", "mu", "kappa",
                              "theta_c", "epsilon", "gamma",
                              "damp_width", "damp_rate", "precision", "delta",
                              "bc_switch_time", "bc_type",
                              "gradient_A_high", "gradient_A_low", "gradient_margin",
                              "kernel"};
        const char *vals[] = {vN, vL, vT, vdt, vm, vmt, veta, vmu, vkappa,
                              vtc, veps, vgam,
                              vdw, vdr, vprec, vdelta,
                              vbcsw, vbc,
                              vgah, vgal, vgm,
                              "fixed_point_q16.16"};
        sfa_add_kvmd(sfa, 0, 0xFFFFFFFF, 0xFFFFFFFF, keys, vals, 22);
    }
    sfa_add_column(sfa, "phi_x",    sfa_dtype, SFA_POSITION, 0);
    sfa_add_column(sfa, "phi_y",    sfa_dtype, SFA_POSITION, 1);
    sfa_add_column(sfa, "phi_z",    sfa_dtype, SFA_POSITION, 2);
    sfa_add_column(sfa, "theta_x",  sfa_dtype, SFA_ANGLE,    0);
    sfa_add_column(sfa, "theta_y",  sfa_dtype, SFA_ANGLE,    1);
    sfa_add_column(sfa, "theta_z",  sfa_dtype, SFA_ANGLE,    2);
    sfa_add_column(sfa, "phi_vx",   sfa_dtype, SFA_VELOCITY, 0);
    sfa_add_column(sfa, "phi_vy",   sfa_dtype, SFA_VELOCITY, 1);
    sfa_add_column(sfa, "phi_vz",   sfa_dtype, SFA_VELOCITY, 2);
    sfa_add_column(sfa, "theta_vx", sfa_dtype, SFA_VELOCITY, 3);
    sfa_add_column(sfa, "theta_vy", sfa_dtype, SFA_VELOCITY, 4);
    sfa_add_column(sfa, "theta_vz", sfa_dtype, SFA_VELOCITY, 5);
    sfa_finalize_header(sfa);
    const char *pn[] = {"f16", "f32", "f64"};
    printf("SFA: %s (12 cols, %s, BSS+zstd, fixed-point kernel)\n\n", c.output, pn[c.precision]);

    sfa_snap(sfa, g, 0.0, c.precision);

    /* Diagnostics file */
    FILE *diag_fp = fopen(c.diag_file, "w");
    fprintf(diag_fp, "t\tE_phi_kin\tE_theta_kin\tE_grad\tE_mass\tE_pot\tE_tgrad\tE_tmass\t"
                     "E_coupling\tE_total\tphi_max\tP_max\tP_int\ttheta_rms\tE_transfer\tE_confine\tf_avg\n");

    int n_steps = (int)lround(c.T / g->dt);
    int diag_every = (int)lround(c.diag_dt / g->dt); if (diag_every < 1) diag_every = 1;
    int snap_every = (int)lround(c.snap_dt / g->dt); if (snap_every < 1) snap_every = 1;

    /* Initial diagnostic */
    double epk, etk, eg, em, ep, etg, etm, ec, et, pm, Pm, etr, ecf, favg;
    compute_energy(g, &fp, &c, &epk, &etk, &eg, &em, &ep, &etg, &etm, &ec, &et, &pm, &Pm, &etr, &ecf, &favg);
    double Pint0 = P_integrated(g), trms0 = theta_rms(g);
    double E0 = et;
    printf("INIT: E_total=%.4e E_pot=%.4f phi_max=%.4f P_int=%.4e theta_rms=%.3e f_avg=%.4f\n\n",
           et, ep, pm, Pint0, trms0, favg);
    fprintf(diag_fp, "%.2f\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\n",
            0.0, epk, etk, eg, em, ep, etg, etm, ec, et, pm, Pm, Pint0, trms0, etr, ecf, favg);

    double wall0 = omp_get_wtime();
    int major = diag_every * 25; if (major < 1) major = 1;
    if (major > n_steps / 10) {
        major = (n_steps / 10 / diag_every) * diag_every;
        if (major < diag_every) major = diag_every;
    }

    int bc_switched = 0;
    int bc_switch_step = (c.bc_switch_time > 0) ? (int)lround(c.bc_switch_time / g->dt) : 0;

    for (int step = 1; step <= n_steps; step++) {
        if (bc_switch_step > 0 && step == bc_switch_step && !bc_switched) {
            c.damp_rate = 0.0;
            fp.DAMP_RATE = 0;
            printf("\n*** BC SWITCH at t=%.1f: absorbing -> periodic ***\n\n", step * g->dt);
            bc_switched = 1;
        }
        verlet_step(g, &fp, &c);
        double t = step * g->dt;

        if (step % snap_every == 0)
            sfa_snap(sfa, g, t, c.precision);

        if (step % diag_every == 0) {
            compute_energy(g, &fp, &c, &epk, &etk, &eg, &em, &ep, &etg, &etm, &ec, &et, &pm, &Pm, &etr, &ecf, &favg);
            double Pint = P_integrated(g), trms = theta_rms(g);
            fprintf(diag_fp, "%.2f\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\n",
                    t, epk, etk, eg, em, ep, etg, etm, ec, et, pm, Pm, Pint, trms, etr, ecf, favg);
            fflush(diag_fp);

            if (step % major == 0) {
                double wall = omp_get_wtime() - wall0;
                double drift = 100 * (et - E0) / (fabs(E0) + 1e-30);
                printf("t=%7.1f E=%.3e (drift %+.3f%%) Ep=%.1f phi=%.3f theta_rms=%.2e f_avg=%.4f "
                       "[%.0f%% %.1fs %.2fms/step]\n",
                       t, et, drift, ep, pm, trms, favg, 100.0 * step / n_steps, wall, 1000 * wall / step);
                fflush(stdout);
            }
        }
    }

    /* Final frame */
    {
        int last_snapped = (n_steps / snap_every) * snap_every;
        int gap = n_steps - last_snapped;
        if (gap > 0)
            sfa_snap(sfa, g, n_steps * g->dt, c.precision);
    }
    uint32_t nf = sfa->total_frames;
    sfa_close(sfa);

    /* Final summary */
    compute_energy(g, &fp, &c, &epk, &etk, &eg, &em, &ep, &etg, &etm, &ec, &et, &pm, &Pm, &etr, &ecf, &favg);
    double trms = theta_rms(g), Pint = P_integrated(g);
    double wall = omp_get_wtime() - wall0;

    printf("\n=== COMPLETE (fixed-point Q16.16) ===\n");
    printf("E_total=%.4e (drift %.3f%%) E_pot=%.4f\n", et, 100 * (et - E0) / (fabs(E0) + 1e-30), ep);
    printf("E_transfer=%.4e E_confine=%.4e f_avg=%.4f\n", etr, ecf, favg);
    printf("phi_max=%.4f P_int=%.4e theta_rms=%.3e\n", pm, Pint, trms);
    printf("SFA: %s (%u frames)\n", c.output, nf);
    printf("[%s] theta_rms grew: %.2e -> %.2e\n", (trms > trms0 + 1e-10) ? "OK" : "WARN", trms0, trms);
    printf("Wall: %.1fs (%.1f min) %.2fms/step\n", wall, wall / 60, 1000 * wall / n_steps);

    fclose(diag_fp);
    grid_free(g);
    return 0;
}
