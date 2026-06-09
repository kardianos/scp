/*  scp_config.h — Shared configuration for CPU and CUDA simulation kernels
 *
 *  Includes: Config struct, defaults, parsing, loading, printing.
 *  Used by both scp_sim.c (CPU) and scp_sim.cu (CUDA).
 */

#ifndef SCP_CONFIG_H
#define SCP_CONFIG_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define NFIELDS 3
#define PI 3.14159265358979323846

/* ================================================================
   Configuration — all parameters with defaults
   ================================================================ */

typedef struct {
    /* Grid */
    int N;
    double L, T, dt_factor;

    /* Physics */
    double m2, mtheta2, eta, mu, kappa;
    double kappa_h;             /* chiral helicity coupling: κ_h P² φ·curl(φ) */
    double alpha_cs;            /* Cosserat strain: α|curl(φ)/2 - θ|² */
    double beta_h;              /* curl-squared hardening: (β/2)|θ|²|∇×φ|² */

    /* Mass coupling mode */
    int mode;
    double inv_alpha, inv_beta, kappa_gamma;
    double gamma_conf;          /* mismatch confining potential: γ|M|²|φ|² (Path 7) */
    double theta_sat;           /* theta saturation threshold (0=disabled) */
    double gamma_conv;          /* theta→phi conversion rate at saturation */

    /* V54 theta self-interaction modes (0=off) */
    double sigma_grad;          /* Test A: gradient-based |∇θ|² conversion */
    double chi_chiral;          /* Test A: chirality bias for gradient mode */
    double sigma_cubic;         /* Test B: cubic self-interaction |θ|²·θ */
    double sigma_freq;          /* Test C: frequency-mismatch conversion */

    /* V54 Lagrangian-derived couplings (energy-conserving) */
    double sigma_cross;         /* Cross-potential: -(σ/2)|θ|²|φ|² */
    double lambda_self;         /* Theta self-potential: -(λ/4)|θ|⁴ */
    double theta_vev;           /* Theta VEV: V = (λ/4)(|θ|²-v²)², 0=standard */
    double tune_dt;             /* Auto-tune interval (0=disabled) */
    int sweep;                  /* Sweep mode: 0=off, 1=parameter sweep */
    double sweep_T;             /* Time per sweep trial (default 50) */

    /* v65 self-tuning (mechanism A): kappa self-organizes to the stability edge at fixed
       conserved charge Q=int|phi|^2. Action pulls kappa down (dE/dkappa>0); a collapse
       sensor pushes it up; the balance is the self-organized-critical edge kappa*(Q). */
    int    self_tune;           /* 0=off, 1=on (kappa becomes dynamical) */
    double st_dt;               /* self-tune update interval (sim time units) */
    double st_eps;              /* stable: drift kappa DOWN by this fraction (action pull) */
    double st_gamma;            /* collapse: push kappa UP by this fraction (stability feedback) */
    double st_pcrit;            /* P_max collapse threshold (settled-state, SOC mode) */
    double st_charge;           /* target charge Q (<=0: use value at first tune step) */
    int    st_project;          /* 1=enforce fixed charge by rescaling phi (mechanism A) */
    /* Bidirectional density homeostasis: regulate core density P_max to st_ptarget by
       adjusting kappa (higher kappa -> lower density). Corrects BOTH dissolution (density
       too low -> kappa down -> re-concentrate) and collapse (too high -> kappa up ->
       spread), so the particle persists. st_ptarget>0 enables this instead of SOC. */
    double st_ptarget;          /* target core density P_max (>0 enables density feedback) */
    double st_gain;             /* proportional gain for the density feedback */

    /* Boundary */
    int bc_type;                /* 0=absorb_sphere, 1=gradient_pinned, 2=periodic */
    double damp_width, damp_rate;
    double bc_switch_time;      /* switch absorb->periodic at this sim time (0=never) */
    double gradient_A_high, gradient_A_low;
    int gradient_margin;

    /* Init */
    char init[32];          /* "oscillon", "braid", "sfa", "exec", "template" */
    double A, sigma, A_bg, ellip, R_tube;
    double delta[3];
    char init_sfa[512];
    int init_frame;
    char init_exec[1024];

    /* Output */
    char output[512];
    char diag_file[512];
    int precision;          /* 0=f16, 1=f32, 2=f64 */
    int output_split;       /* CUDA only: split output files */
    double snap_dt, diag_dt;
    double burst_start, burst_end;
    double burst_every;

    /* Vector frame output (FRVD in same SFA) */
    double vec_snap_dt;
    int vec_iframe_interval;
    double vec_delta_tol;
    int vec_block_size;
} Config;

static Config cfg_defaults(void) {
    Config c;
    memset(&c, 0, sizeof(c));
    c.N = 128;  c.L = 10.0;  c.T = 200.0;  c.dt_factor = 0.025;
    c.m2 = 2.25;  c.mtheta2 = 0.0;  c.eta = 0.5;
    c.mu = -41.345;  c.kappa = 50.0;  c.kappa_h = 0.0;  c.alpha_cs = 0.0;  c.beta_h = 0.0;
    c.mode = 0;  c.inv_alpha = 2.25;  c.inv_beta = 5.0;  c.kappa_gamma = 2.0;  c.gamma_conf = 0;
    c.theta_sat = 0;  c.gamma_conv = 0;
    c.sigma_grad = 0;  c.chi_chiral = 0;  c.sigma_cubic = 0;  c.sigma_freq = 0;
    c.sigma_cross = 0;  c.lambda_self = 0;  c.theta_vev = 0;  c.tune_dt = 0;
    c.sweep = 0;  c.sweep_T = 50;
    c.self_tune = 0;  c.st_dt = 2.0;  c.st_eps = 0.04;  c.st_gamma = 0.12;
    c.st_pcrit = 2.0;  c.st_charge = 0;  c.st_project = 1;
    c.st_ptarget = 0;  c.st_gain = 0.10;
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
    c.precision = 1;  /* f32 */
    c.output_split = 0;
    c.snap_dt = 5.0;  c.diag_dt = 2.0;
    c.burst_start = 0;  c.burst_end = 0;  c.burst_every = 0;
    c.vec_snap_dt = 0;
    c.vec_iframe_interval = 100;
    c.vec_delta_tol = 0.01;  c.vec_block_size = 8;
    return c;
}

static void cfg_set(Config *c, const char *key, const char *val) {
    if      (!strcmp(key,"N"))           c->N = atoi(val);
    else if (!strcmp(key,"L"))           c->L = atof(val);
    else if (!strcmp(key,"T"))           c->T = atof(val);
    else if (!strcmp(key,"dt_factor"))   c->dt_factor = atof(val);
    else if (!strcmp(key,"m"))         { double m = atof(val); c->m2 = m*m; }
    else if (!strcmp(key,"m_theta"))   { double m = atof(val); c->mtheta2 = m*m; }
    else if (!strcmp(key,"eta"))         c->eta = atof(val);
    else if (!strcmp(key,"mu"))          c->mu = atof(val);
    else if (!strcmp(key,"kappa"))       c->kappa = atof(val);
    else if (!strcmp(key,"kappa_h"))    c->kappa_h = atof(val);
    else if (!strcmp(key,"alpha_cs"))   c->alpha_cs = atof(val);
    else if (!strcmp(key,"beta_h"))    c->beta_h = atof(val);
    else if (!strcmp(key,"mode"))        c->mode = atoi(val);
    else if (!strcmp(key,"inv_alpha"))   c->inv_alpha = atof(val);
    else if (!strcmp(key,"inv_beta"))    c->inv_beta = atof(val);
    else if (!strcmp(key,"kappa_gamma")) c->kappa_gamma = atof(val);
    else if (!strcmp(key,"gamma_conf")) c->gamma_conf = atof(val);
    else if (!strcmp(key,"theta_sat"))  c->theta_sat = atof(val);
    else if (!strcmp(key,"gamma_conv")) c->gamma_conv = atof(val);
    else if (!strcmp(key,"sigma_grad"))  c->sigma_grad = atof(val);
    else if (!strcmp(key,"chi_chiral")) c->chi_chiral = atof(val);
    else if (!strcmp(key,"sigma_cubic")) c->sigma_cubic = atof(val);
    else if (!strcmp(key,"sigma_freq")) c->sigma_freq = atof(val);
    else if (!strcmp(key,"sigma_cross")) c->sigma_cross = atof(val);
    else if (!strcmp(key,"lambda_self")) c->lambda_self = atof(val);
    else if (!strcmp(key,"theta_vev"))  c->theta_vev = atof(val);
    else if (!strcmp(key,"tune_dt"))    c->tune_dt = atof(val);
    else if (!strcmp(key,"self_tune"))  c->self_tune = atoi(val);
    else if (!strcmp(key,"st_dt"))      c->st_dt = atof(val);
    else if (!strcmp(key,"st_eps"))     c->st_eps = atof(val);
    else if (!strcmp(key,"st_gamma"))   c->st_gamma = atof(val);
    else if (!strcmp(key,"st_pcrit"))   c->st_pcrit = atof(val);
    else if (!strcmp(key,"st_charge"))  c->st_charge = atof(val);
    else if (!strcmp(key,"st_project")) c->st_project = atoi(val);
    else if (!strcmp(key,"st_ptarget")) c->st_ptarget = atof(val);
    else if (!strcmp(key,"st_gain"))    c->st_gain = atof(val);
    else if (!strcmp(key,"sweep"))      c->sweep = atoi(val);
    else if (!strcmp(key,"sweep_T"))    c->sweep_T = atof(val);
    else if (!strcmp(key,"bc_type"))      c->bc_type = atoi(val);
    else if (!strcmp(key,"damp_width"))  c->damp_width = atof(val);
    else if (!strcmp(key,"bc_switch_time")) c->bc_switch_time = atof(val);
    else if (!strcmp(key,"damp_rate"))   c->damp_rate = atof(val);
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
    else if (!strcmp(key,"output_split")) c->output_split = atoi(val);
    else if (!strcmp(key,"vec_snap_dt"))        c->vec_snap_dt = atof(val);
    else if (!strcmp(key,"vec_iframe_interval")) c->vec_iframe_interval = atoi(val);
    else if (!strcmp(key,"vec_delta_tol"))      c->vec_delta_tol = atof(val);
    else if (!strcmp(key,"vec_block_size"))     c->vec_block_size = atoi(val);
    else if (!strcmp(key,"burst_start")) c->burst_start = atof(val);
    else if (!strcmp(key,"burst_end"))   c->burst_end = atof(val);
    else if (!strcmp(key,"burst_every")) c->burst_every = atof(val);
    else if (!strcmp(key,"precision")) {
        if      (!strcmp(val,"f16")) c->precision = 0;
        else if (!strcmp(val,"f32")) c->precision = 1;
        else if (!strcmp(val,"f64")) c->precision = 2;
    }
    else if (!strcmp(key,"vec_output")) {}           /* deprecated */
    else if (!strcmp(key,"vec_kframe_interval")) {}  /* deprecated */
    else if (!strcmp(key,"vec_n_fields")) {}         /* deprecated */
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
    printf("=== scp_sim_c4: V44 + Cosserat + Curl²-Hardening (+ optional chiral) ===\n");
    printf("d²φ/dt² = ∇²φ - m²φ - V'(P) + η·curl(θ) - α·curl(M) - 2curl(Q)\n");
    printf("d²θ/dt² = ∇²θ - m_θ²θ + η·curl(φ) + 2α·M - β|∇×φ|²θ\n");
    printf("  M = curl(φ)/2 - θ,  Q = (β/2)|θ|²curl(φ)\n\n");
    printf("Grid:    N=%d L=%.1f T=%.0f dt_factor=%.4f\n", c->N, c->L, c->T, c->dt_factor);
    printf("Physics: m²=%.4f m_θ²=%.4f η=%.3f μ=%.3f κ=%.1f κ_h=%.3f α=%.3f β=%.3f\n",
           c->m2, c->mtheta2, c->eta, c->mu, c->kappa, c->kappa_h, c->alpha_cs, c->beta_h);
    printf("Mode:    %d", c->mode);
    if (c->mode == 1) printf(" (inverse: α=%.3f β=%.3f)", c->inv_alpha, c->inv_beta);
    if (c->mode == 3) printf(" (density-κ: γ=%.3f)", c->kappa_gamma);
    if (c->bc_type == 0) {
        printf("\nBC:      absorbing sphere (width=%.1f rate=%.4f)", c->damp_width, c->damp_rate);
        if (c->bc_switch_time > 0) printf(" -> periodic at t=%.0f", c->bc_switch_time);
        printf("\n");
    } else if (c->bc_type == 1)
        printf("\nBC:      gradient pinned (A_high=%.3f A_low=%.3f margin=%d)\n",
               c->gradient_A_high, c->gradient_A_low, c->gradient_margin);
    else if (c->bc_type == 2)
        printf("\nBC:      periodic\n");
    printf("Init:    %s", c->init);
    if (!strcmp(c->init, "sfa")) printf(" (%s frame=%d)", c->init_sfa, c->init_frame);
    if (!strcmp(c->init, "exec")) printf(" (%s)", c->init_exec);
    if (!strcmp(c->init, "template")) printf(" (%s)", c->init_sfa);
    printf("\nOutput:  %s (%s, snap=%.1f diag=%.1f)\n\n",
           c->output, prec_names[c->precision], c->snap_dt, c->diag_dt);
}

/* ================================================================
   f16 conversion helpers
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

#endif /* SCP_CONFIG_H */
