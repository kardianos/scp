/*  scp_sim_unified.c — V49 Variant B: Unified Transfer Potential
 *
 *  Modified copy of sfa/sim/scp_sim.c with the unified transfer potential
 *  replacing V(P) + lambda_theta.
 *
 *  Physics:
 *    W(P, Theta) = V_base(P) * f(Theta)           (transfer, negative)
 *    W_confine   = |V_base(P)| * gamma * ln(1 + Theta/theta_c)  (positive)
 *
 *    f(Theta) = epsilon + (1-epsilon) * Theta / (Theta + theta_c)
 *
 *    Phi force:   -dVdP * dPda * f_phi + eta * curl(theta)
 *                 where f_phi = f_transfer - f_confine (confine OPPOSES binding)
 *    Theta force: -V_abs * (df_dTheta + confine_deriv) * 2 * theta_a + eta * curl(phi)
 *                 where V_abs = -V_base (positive), df_dTheta and confine_deriv as in proposal
 *
 *  New params: theta_c (default 0.05), epsilon (default 0.1), gamma (default 0.1)
 *  Removed:    lambda_theta, eta1
 *
 *  Build: gcc -O3 -march=native -fopenmp -o scp_sim_unified scp_sim_unified.c -lzstd -lm
 *  Run:   OMP_NUM_THREADS=8 ./scp_sim_unified config.cfg [-key value ...]
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
   Configuration — all parameters with defaults
   ================================================================ */

typedef struct {
    /* Grid */
    int N;
    double L, T, dt_factor;

    /* Physics */
    double m2, mtheta2, eta, mu, kappa;

    /* Unified transfer potential */
    double theta_c, epsilon, gamma_conf;

    /* Mass coupling mode */
    int mode;
    double inv_alpha, inv_beta, kappa_gamma;

    /* Boundary */
    int bc_type;                /* 0=absorb_sphere (default), 1=gradient_pinned */
    double damp_width, damp_rate;
    double bc_switch_time;      /* switch absorb->periodic at this sim time (0=never) */
    double gradient_A_high, gradient_A_low;
    int gradient_margin;

    /* Init */
    char init[32];          /* "oscillon", "braid", "sfa", "exec" */
    double A, sigma, A_bg, ellip, R_tube;
    double delta[3];
    char init_sfa[512];
    int init_frame;
    char init_exec[1024];

    /* Output */
    char output[512];
    char diag_file[512];
    int precision;          /* 0=f16, 1=f32, 2=f64 */
    int output_split;       /* 0=single .sfa, 1=header .sfa + per-frame .sfp */
    double snap_dt, diag_dt;
} Config;

static Config cfg_defaults(void) {
    Config c = {0};
    c.N = 128;  c.L = 10.0;  c.T = 200.0;  c.dt_factor = 0.025;
    c.m2 = 2.25;  c.mtheta2 = 0.0;  c.eta = 0.5;
    c.mu = -41.345;  c.kappa = 50.0;
    c.theta_c = 0.05;  c.epsilon = 0.1;  c.gamma_conf = 0.1;
    c.mode = 0;  c.inv_alpha = 2.25;  c.inv_beta = 5.0;  c.kappa_gamma = 2.0;
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
    else if (!strcmp(key,"theta_c"))     c->theta_c = atof(val);
    else if (!strcmp(key,"epsilon"))     c->epsilon = atof(val);
    else if (!strcmp(key,"gamma"))       c->gamma_conf = atof(val);
    else if (!strcmp(key,"mode"))        c->mode = atoi(val);
    else if (!strcmp(key,"inv_alpha"))   c->inv_alpha = atof(val);
    else if (!strcmp(key,"inv_beta"))    c->inv_beta = atof(val);
    else if (!strcmp(key,"kappa_gamma")) c->kappa_gamma = atof(val);
    else if (!strcmp(key,"bc_type"))      c->bc_type = atoi(val);
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
    /* Silently ignore old params that may appear in SFA KVMD */
    else if (!strcmp(key,"lambda_theta") || !strcmp(key,"eta1")) { /* removed */ }
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
        /* trim key */
        char *ke = eq - 1;
        while (ke > key && (*ke == ' ' || *ke == '\t')) *ke-- = '\0';
        /* trim val */
        while (*val == ' ' || *val == '\t') val++;
        char *ve = val + strlen(val) - 1;
        while (ve > val && (*ve == '\n' || *ve == '\r' || *ve == ' ')) *ve-- = '\0';
        /* strip inline comments */
        char *hash = strchr(val, '#');
        if (hash) { *hash = '\0'; ve = hash - 1; while (ve > val && *ve == ' ') *ve-- = '\0'; }
        cfg_set(c, key, val);
    }
    fclose(fp);
}

static void cfg_print(const Config *c) {
    const char *prec_names[] = {"f16", "f32", "f64"};
    printf("=== scp_sim_unified: V49 Variant B — Unified Transfer Potential ===\n");
    printf("d²phi/dt² = nabla²phi - m²phi - dVdP*dPda*f_phi + eta*curl(theta)\n");
    printf("d²theta/dt² = nabla²theta - |V_base|*(df/dTheta + gamma/(Theta+theta_c))*2*theta + eta*curl(phi)\n");
    printf("  V_base(P) = (mu/2)*P²/(1+kappa*P²)\n");
    printf("  f(Theta) = eps + (1-eps)*Theta/(Theta+theta_c)\n");
    printf("  f_phi = f_transfer - gamma*ln(1+Theta/theta_c)  [confine opposes binding]\n\n");
    printf("Grid:    N=%d L=%.1f T=%.0f dt_factor=%.4f\n", c->N, c->L, c->T, c->dt_factor);
    printf("Physics: m²=%.4f m_theta²=%.4f eta=%.3f mu=%.3f kappa=%.1f\n",
           c->m2, c->mtheta2, c->eta, c->mu, c->kappa);
    printf("Transfer: theta_c=%.4f epsilon=%.4f gamma=%.4f\n",
           c->theta_c, c->epsilon, c->gamma_conf);
    printf("Mode:    %d", c->mode);
    if (c->mode == 1) printf(" (inverse: alpha=%.3f beta=%.3f)", c->inv_alpha, c->inv_beta);
    if (c->mode == 3) printf(" (density-kappa: gamma=%.3f)", c->kappa_gamma);
    if (c->bc_type == 0)
        printf("\nBC:      absorbing sphere (width=%.1f rate=%.4f)\n", c->damp_width, c->damp_rate);
    else if (c->bc_type == 1)
        printf("\nBC:      gradient pinned (A_high=%.3f A_low=%.3f margin=%d)\n",
               c->gradient_A_high, c->gradient_A_low, c->gradient_margin);
    printf("Init:    %s", c->init);
    if (!strcmp(c->init, "sfa")) printf(" (%s frame=%d)", c->init_sfa, c->init_frame);
    if (!strcmp(c->init, "exec")) printf(" (%s)", c->init_exec);
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

/* ================================================================
   Grid: 18 arrays (6 fields x {val, vel, acc})
   ================================================================ */

typedef struct {
    double *mem;
    double *phi[NFIELDS], *phi_vel[NFIELDS], *phi_acc[NFIELDS];
    double *theta[NFIELDS], *theta_vel[NFIELDS], *theta_acc[NFIELDS];
    double *pin_phi[NFIELDS], *pin_vel[NFIELDS];  /* pinned BC storage (bc_type=1) */
    double *pin_theta[NFIELDS], *pin_tvel[NFIELDS];
    int N; long N3;
    double L, dx, dt;
} Grid;

static Grid *grid_alloc(const Config *c) {
    Grid *g = calloc(1, sizeof(Grid));
    g->N  = c->N;
    g->N3 = (long)c->N * c->N * c->N;
    g->L  = c->L;
    g->dx = 2.0 * c->L / (c->N - 1);
    g->dt = c->dt_factor * g->dx;

    long total = 18 * g->N3;
    printf("Allocating %.2f GB (%ld doubles, N=%d, 6 fields)\n", total*8.0/1e9, total, c->N);
    g->mem = malloc(total * sizeof(double));
    if (!g->mem) { fprintf(stderr, "FATAL: malloc failed\n"); exit(1); }
    memset(g->mem, 0, total * sizeof(double));

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

/* Save current state as pinned boundary values (for bc_type=1) */
static void grid_save_pinned(Grid *g) {
    for (int a = 0; a < NFIELDS; a++) {
        g->pin_phi[a]   = malloc(g->N3 * sizeof(double));
        g->pin_vel[a]   = malloc(g->N3 * sizeof(double));
        g->pin_theta[a] = malloc(g->N3 * sizeof(double));
        g->pin_tvel[a]  = malloc(g->N3 * sizeof(double));
        memcpy(g->pin_phi[a],   g->phi[a],       g->N3 * sizeof(double));
        memcpy(g->pin_vel[a],   g->phi_vel[a],   g->N3 * sizeof(double));
        memcpy(g->pin_theta[a], g->theta[a],     g->N3 * sizeof(double));
        memcpy(g->pin_tvel[a],  g->theta_vel[a], g->N3 * sizeof(double));
    }
}

/* ================================================================
   Curl helper
   ================================================================ */

static inline double curl_component(double *F[3], int a,
    long n_ip, long n_im, long n_jp, long n_jm, long n_kp, long n_km,
    double idx1) {
    if (a == 0) return (F[2][n_jp] - F[2][n_jm] - F[1][n_kp] + F[1][n_km]) * idx1;
    if (a == 1) return (F[0][n_kp] - F[0][n_km] - F[2][n_ip] + F[2][n_im]) * idx1;
    return            (F[1][n_ip] - F[1][n_im] - F[0][n_jp] + F[0][n_jm]) * idx1;
}

/* ================================================================
   Initialization
   ================================================================ */

static void init_oscillon(Grid *g, const Config *c) {
    const int N = g->N, NN = N * N;
    const double L = g->L, dx = g->dx;
    printf("Init: oscillon (A=%.3f sigma=%.3f)\n", c->A, c->sigma);
    for (int i = 0; i < N; i++) { double x = -L + i*dx;
    for (int j = 0; j < N; j++) { double y = -L + j*dx;
    for (int k = 0; k < N; k++) { double z = -L + k*dx;
        long idx = (long)i*NN + j*N + k;
        double r2 = x*x + y*y + z*z;
        double env = c->A * exp(-r2 / (2.0 * c->sigma * c->sigma));
        for (int a = 0; a < NFIELDS; a++)
            g->phi[a][idx] = env * cos(c->delta[a]);
    }}}
}

static void init_braid(Grid *g, const Config *c) {
    const int N = g->N, NN = N * N;
    const double L = g->L, dx = g->dx;
    const double kw = PI/L, omega = sqrt(kw*kw + c->m2);
    const double sx = 1+c->ellip, sy = 1-c->ellip;
    const double inv2R2 = 1.0/(2*c->R_tube*c->R_tube);
    const double k_bg = PI/L, omega_bg = sqrt(k_bg*k_bg + c->m2);
    printf("Init: braid (R=%.1f ellip=%.4f A=%.2f A_bg=%.2f)\n",
           c->R_tube, c->ellip, c->A, c->A_bg);
    for (int i = 0; i < N; i++) { double x = -L + i*dx;
    for (int j = 0; j < N; j++) { double y = -L + j*dx;
    for (int k = 0; k < N; k++) { double z = -L + k*dx;
        long idx = (long)i*NN + j*N + k;
        double r2e = x*x/(sx*sx) + y*y/(sy*sy);
        double env = exp(-r2e * inv2R2);
        for (int a = 0; a < NFIELDS; a++) {
            double ph = kw*z + c->delta[a];
            double ph_bg = k_bg*z + 2*PI*a/3.0;
            g->phi[a][idx] = c->A*env*cos(ph) + c->A_bg*cos(ph_bg);
            g->phi_vel[a][idx] = omega*c->A*env*sin(ph) + omega_bg*c->A_bg*sin(ph_bg);
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

    /* Walk columns and load by semantic+component */
    int loaded[12] = {0};  /* track which of our 12 fields got loaded */
    uint64_t off = 0;
    for (int col = 0; col < sfa->n_columns; col++) {
        int dtype = sfa->columns[col].dtype;
        int sem   = sfa->columns[col].semantic;
        int comp  = sfa->columns[col].component;
        int es    = sfa_dtype_size[dtype];
        uint8_t *src = (uint8_t*)buf + off;

        double *target = NULL;
        int slot = -1;
        if (sem == SFA_POSITION && comp < 3) { target = g->phi[comp]; slot = comp; }
        else if (sem == SFA_ANGLE && comp < 3) { target = g->theta[comp]; slot = 3+comp; }
        else if (sem == SFA_VELOCITY && comp < 3) { target = g->phi_vel[comp]; slot = 6+comp; }
        else if (sem == SFA_VELOCITY && comp >= 3 && comp < 6) { target = g->theta_vel[comp-3]; slot = 9+comp-3; }

        if (target) {
            long N3 = g->N3;
            if (dtype == SFA_F64)
                for (long i = 0; i < N3; i++) target[i] = ((double*)src)[i];
            else if (dtype == SFA_F32)
                for (long i = 0; i < N3; i++) target[i] = (double)((float*)src)[i];
            else if (dtype == SFA_F16)
                for (long i = 0; i < N3; i++) target[i] = f16_to_f64(((uint16_t*)src)[i]);
            loaded[slot] = 1;
            printf("  Loaded col %d '%s' (sem=%d comp=%d dtype=%d) -> slot %d\n",
                   col, sfa->columns[col].name, sem, comp, dtype, slot);
        }
        off += (uint64_t)g->N3 * es;
    }

    /* Check what was loaded */
    int n_fields = 0, n_vels = 0;
    for (int i = 0; i < 6; i++) n_fields += loaded[i];
    for (int i = 6; i < 12; i++) n_vels += loaded[i];
    printf("  Loaded: %d/6 field arrays, %d/6 velocity arrays\n", n_fields, n_vels);
    if (n_vels == 0)
        printf("  WARNING: no velocity data -- starting from rest (cold restart)\n");
    if (fabs(sfa->Lx - g->L) > 0.01)
        printf("  WARNING: SFA L=%.2f != config L=%.2f\n", sfa->Lx, g->L);

    free(buf);
    sfa_close(sfa);
}

static void init_from_exec(Grid *g, const Config *c) {
    printf("Init: exec '%s'\n", c->init_exec);
    /* Write to temp file, then load as SFA */
    char tmppath[512];
    snprintf(tmppath, sizeof(tmppath), "/tmp/scp_sim_init_%d.sfa", (int)getpid());
    char cmd[2048];
    snprintf(cmd, sizeof(cmd), "%s > %s", c->init_exec, tmppath);
    printf("  Running: %s\n", cmd);
    int ret = system(cmd);
    if (ret != 0) { fprintf(stderr, "FATAL: exec init failed (exit %d)\n", ret); exit(1); }

    /* Temporarily redirect to SFA loader */
    Config tmp = *c;
    strncpy(tmp.init_sfa, tmppath, 511);
    tmp.init_frame = -1;
    init_from_sfa(g, &tmp);
    remove(tmppath);
}

/* init=template: Load a small template SFA and stamp it into the grid.
 * Background generated analytically (with gradient if bc_type=1).
 * Uses init_sfa for the template path. */
static void init_template(Grid *g, const Config *c) {
    printf("Init: template '%s'\n", c->init_sfa);

    const int N=g->N, NN=N*N;
    const double L=g->L, dx=g->dx;
    const double k_bg=PI/L, omega_bg=sqrt(k_bg*k_bg+c->m2);

    /* Generate background (with optional gradient) */
    #pragma omp parallel for schedule(static)
    for (long idx=0; idx<g->N3; idx++) {
        int i=(int)(idx/NN), j=(int)((idx/N)%N), k=(int)(idx%N);
        double x=-L+i*dx, z=-L+k*dx;
        double A_bg_x = c->A_bg;
        if (c->bc_type == 1) {
            double frac = (x + L) / (2.0 * L);
            A_bg_x = c->gradient_A_high * (1.0 - frac) + c->gradient_A_low * frac;
        }
        for (int a=0; a<NFIELDS; a++) {
            double ph_bg = k_bg*z + 2*PI*a/3.0;
            g->phi[a][idx] = A_bg_x * cos(ph_bg);
            g->phi_vel[a][idx] = omega_bg * A_bg_x * sin(ph_bg);
        }
    }

    /* Load template and stamp at center */
    SFA *tmpl = sfa_open(c->init_sfa);
    if (!tmpl) { fprintf(stderr, "FATAL: cannot open template '%s'\n", c->init_sfa); exit(1); }
    int TN=tmpl->Nx; double TL=tmpl->Lx, Tdx=2.0*TL/(TN-1);
    long TN3=(long)TN*TN*TN; int TNN=TN*TN;
    printf("  Template: %d^3, L=%.2f, dx=%.4f\n", TN, TL, Tdx);

    int frame = tmpl->total_frames - 1;
    void *buf = malloc(tmpl->frame_bytes);
    sfa_read_frame(tmpl, frame, buf);

    double *tphi[3], *tvel[3], *ttheta[3], *ttvel[3];
    for (int a=0; a<3; a++) {
        tphi[a]=calloc(TN3,sizeof(double)); tvel[a]=calloc(TN3,sizeof(double));
        ttheta[a]=calloc(TN3,sizeof(double)); ttvel[a]=calloc(TN3,sizeof(double));
    }
    uint64_t off=0;
    for (int col=0; col<tmpl->n_columns; col++) {
        int dtype=tmpl->columns[col].dtype, sem=tmpl->columns[col].semantic, comp=tmpl->columns[col].component;
        int es=sfa_dtype_size[dtype]; uint8_t *src=(uint8_t*)buf+off;
        double *target=NULL;
        if (sem==SFA_POSITION&&comp<3) target=tphi[comp];
        else if (sem==SFA_ANGLE&&comp<3) target=ttheta[comp];
        else if (sem==SFA_VELOCITY&&comp<3) target=tvel[comp];
        else if (sem==SFA_VELOCITY&&comp>=3&&comp<6) target=ttvel[comp-3];
        if (target) {
            if (dtype==SFA_F64) for(long i=0;i<TN3;i++) target[i]=((double*)src)[i];
            else if (dtype==SFA_F32) for(long i=0;i<TN3;i++) target[i]=(double)((float*)src)[i];
            else if (dtype==SFA_F16) for(long i=0;i<TN3;i++) target[i]=f16_to_f64(((uint16_t*)src)[i]);
        }
        off += (uint64_t)TN3 * es;
    }
    free(buf); sfa_close(tmpl);

    int ci=N/2, cj=N/2, ck=N/2, half=TN/2;
    int placed=0;
    for (int ti=0; ti<TN; ti++) {
        int gi=ci+ti-half; if(gi<0||gi>=N) continue;
        for (int tj=0; tj<TN; tj++) {
            int gj=cj+tj-half; if(gj<0||gj>=N) continue;
            for (int tk=0; tk<TN; tk++) {
                int gk=ck+tk-half; if(gk<0||gk>=N) continue;
                long tidx=(long)ti*TNN+tj*TN+tk, gidx=(long)gi*NN+gj*N+gk;
                double tz=-TL+tk*Tdx;
                for (int a=0; a<NFIELDS; a++) {
                    double ph_bg_t=k_bg*tz+2*PI*a/3.0;
                    double bg_phi=c->A_bg*cos(ph_bg_t), bg_vel=omega_bg*c->A_bg*sin(ph_bg_t);
                    g->phi[a][gidx] += tphi[a][tidx] - bg_phi;
                    g->phi_vel[a][gidx] += tvel[a][tidx] - bg_vel;
                    g->theta[a][gidx] += ttheta[a][tidx];
                    g->theta_vel[a][gidx] += ttvel[a][tidx];
                }
                placed++;
            }
        }
    }
    for (int a=0;a<3;a++) { free(tphi[a]);free(tvel[a]);free(ttheta[a]);free(ttvel[a]); }
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
   Forces: Unified Transfer Potential (V49 Variant B)

   V_base(P) = (mu/2) * P^2 / (1 + kappa * P^2)       (NEGATIVE, mu < 0)
   V_abs     = -V_base                                  (POSITIVE)
   f_transfer = epsilon + (1-epsilon) * Theta / (Theta + theta_c)
   f_confine  = gamma * ln(1 + Theta / theta_c)
   f_phi      = f_transfer - f_confine  (confine OPPOSES binding in phi)

   Phi force:   lap - m^2*phi - dVdP*dPda*f_phi + eta*curl(theta)
   Theta force: lap - V_abs*(df_dTheta + confine_deriv)*2*theta + eta*curl(phi)
   ================================================================ */

static void compute_forces(Grid *g, const Config *c) {
    const int N = g->N, NN = N * N;
    const long N3 = g->N3;
    const double idx2 = 1.0 / (g->dx * g->dx);
    const double idx1 = 1.0 / (2.0 * g->dx);
    const double MU = c->mu, KAPPA = c->kappa, MASS2 = c->m2;
    const double MTHETA2 = c->mtheta2, ETA = c->eta;
    const double THETA_C = c->theta_c, EPSILON = c->epsilon, GAMMA = c->gamma_conf;
    const int MODE = c->mode;
    const double KG = c->kappa_gamma, IA = c->inv_alpha, IB = c->inv_beta;

    #pragma omp parallel for schedule(static)
    for (long idx = 0; idx < N3; idx++) {
        int i = (int)(idx / NN), j = (int)((idx / N) % N), k = (int)(idx % N);
        int ip=(i+1)%N, im=(i-1+N)%N, jp=(j+1)%N, jm=(j-1+N)%N;
        int kp=(k+1)%N, km=(k-1+N)%N;
        long n_ip=(long)ip*NN+j*N+k, n_im=(long)im*NN+j*N+k;
        long n_jp=(long)i*NN+jp*N+k, n_jm=(long)i*NN+jm*N+k;
        long n_kp=(long)i*NN+j*N+kp, n_km=(long)i*NN+j*N+km;

        double p0=g->phi[0][idx], p1=g->phi[1][idx], p2=g->phi[2][idx];
        double sig = p0*p0 + p1*p1 + p2*p2;
        double me2 = (MODE==1) ? IA/(1.0+IB*sig) : MASS2;
        double keff = (MODE==3) ? KAPPA/(1.0+KG*sig) : KAPPA;

        double P = p0*p1*p2, P2 = P*P;
        double den = 1.0 + keff*P2;
        double V_base = (MU/2.0) * P2 / den;           /* NEGATIVE (mu < 0) */
        double V_abs = -V_base;                          /* POSITIVE */
        double dVdP = MU * P / (den*den);               /* NEGATIVE */

        double t2c = 0;
        if (MODE == 3 && KG > 0) {
            double D = 1.0 + KG*sig + KAPPA*P2;
            t2c = MU * KG * KAPPA * P2 * P2 / (D*D);
        }

        /* Theta energy density |theta|^2 */
        double t0=g->theta[0][idx], t1=g->theta[1][idx], t2=g->theta[2][idx];
        double Theta = t0*t0 + t1*t1 + t2*t2;
        double Theta_sum = Theta + THETA_C;

        /* Transfer function: f(Theta) = eps + (1-eps)*Theta/(Theta+theta_c) */
        double f_transfer = EPSILON + (1.0-EPSILON) * Theta / Theta_sum;

        /* Confinement function: gamma * ln(1 + Theta/theta_c) */
        double f_confine = GAMMA * log(1.0 + Theta / THETA_C);

        /* Phi force modulation: transfer ATTRACTS, confine REPELS.
         * W = V_base * f_transfer  →  dW/dP = dVdP * f_transfer
         * W_confine = |V_base| * f_confine  →  dWc/dP = -dVdP * f_confine
         * Total dU/dP = dVdP * (f_transfer - f_confine)
         * Force: -dVdP * (f_transfer - f_confine) * dPda
         *
         * The MINUS on f_confine is because |V_base| = -V_base has
         * the opposite P-derivative. Confinement COSTS energy when P
         * grows, providing negative feedback that prevents runaway. */
        double f_phi = f_transfer - f_confine;

        for (int a = 0; a < NFIELDS; a++) {
            double lap = (g->phi[a][n_ip]+g->phi[a][n_im]+g->phi[a][n_jp]
                        +g->phi[a][n_jm]+g->phi[a][n_kp]+g->phi[a][n_km]
                        -6.0*g->phi[a][idx]) * idx2;
            double dPda = (a==0)?p1*p2:(a==1)?p0*p2:p0*p1;
            double ct = curl_component(g->theta, a, n_ip,n_im,n_jp,n_jm,n_kp,n_km, idx1);
            g->phi_acc[a][idx] = lap - me2*g->phi[a][idx] - dVdP*dPda*f_phi - t2c*g->phi[a][idx]
                               + ETA*ct;
        }

        /* Theta force: F = +V_abs × (df/dΘ - confine_deriv) × 2 × theta_a
         * Small Θ: df > confine → DRIVES theta growth (energy transfer)
         * Large Θ: confine > df → CONFINES theta (pulls back)
         * Equilibrium where df_dTheta = confine_deriv. */
        double df_dTheta = (1.0-EPSILON) * THETA_C / (Theta_sum*Theta_sum);
        double confine_deriv = GAMMA / Theta_sum;
        double theta_drive = V_abs * (df_dTheta - confine_deriv) * 2.0;
        /* theta_drive > 0 at small Θ (growth), < 0 at large Θ (confinement) */

        for (int a = 0; a < NFIELDS; a++) {
            double lapt = (g->theta[a][n_ip]+g->theta[a][n_im]+g->theta[a][n_jp]
                         +g->theta[a][n_jm]+g->theta[a][n_kp]+g->theta[a][n_km]
                         -6.0*g->theta[a][idx]) * idx2;
            double cp = curl_component(g->phi, a, n_ip,n_im,n_jp,n_jm,n_kp,n_km, idx1);
            g->theta_acc[a][idx] = lapt - MTHETA2*g->theta[a][idx]
                                 + theta_drive * g->theta[a][idx] + ETA*cp;
        }
    }
}

/* ================================================================
   Boundary conditions
   ================================================================ */

/* bc_type=0: Spherical absorbing boundary */
static void apply_damping(Grid *g, const Config *c) {
    const int N=g->N, NN=N*N;
    const double L=g->L, dx=g->dx, DW=c->damp_width, DR=c->damp_rate;
    if (DW <= 0 || DR <= 0) return;
    const double R_damp = L - DW;

    #pragma omp parallel for schedule(static)
    for (long idx = 0; idx < g->N3; idx++) {
        int i=(int)(idx/NN), j=(int)((idx/N)%N), k=(int)(idx%N);
        double x=-L+i*dx, y=-L+j*dx, z=-L+k*dx;
        double r = sqrt(x*x+y*y+z*z);
        if (r > R_damp) {
            double s = (r-R_damp)/DW; if (s>1) s=1;
            double d = 1.0 - DR*s*s;
            for (int a=0;a<NFIELDS;a++) { g->phi_vel[a][idx]*=d; g->theta_vel[a][idx]*=d; }
        }
    }
}

/* bc_type=1: Gradient pinned boundary
 * x-direction: pin first/last margin slabs to initial values (maintains gradient)
 * y,z-directions: linear extrapolation from interior (free-floating outflow)
 * z-direction: periodic (handled by Laplacian stencil's modular indexing)
 */
static void apply_gradient_bc(Grid *g, const Config *c) {
    const int N=g->N, NN=N*N;
    const int M = c->gradient_margin;

    /* x-direction: pin boundary slabs to saved initial state */
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

    /* y-direction: linear extrapolation from interior (free-floating) */
    #pragma omp parallel for schedule(static)
    for (int i = M; i < N - M; i++) {
        for (int k = 0; k < N; k++) {
            for (int a = 0; a < NFIELDS; a++) {
                /* j=0 edge: extrapolate from j=1,2 */
                long idx0 = (long)i*NN + 0*N + k;
                long idx1 = (long)i*NN + 1*N + k;
                long idx2 = (long)i*NN + 2*N + k;
                g->phi[a][idx0] = 2*g->phi[a][idx1] - g->phi[a][idx2];
                g->phi_vel[a][idx0] = 2*g->phi_vel[a][idx1] - g->phi_vel[a][idx2];
                g->theta[a][idx0] = 2*g->theta[a][idx1] - g->theta[a][idx2];
                g->theta_vel[a][idx0] = 2*g->theta_vel[a][idx1] - g->theta_vel[a][idx2];

                /* j=N-1 edge: extrapolate from j=N-2,N-3 */
                idx0 = (long)i*NN + (N-1)*N + k;
                idx1 = (long)i*NN + (N-2)*N + k;
                idx2 = (long)i*NN + (N-3)*N + k;
                g->phi[a][idx0] = 2*g->phi[a][idx1] - g->phi[a][idx2];
                g->phi_vel[a][idx0] = 2*g->phi_vel[a][idx1] - g->phi_vel[a][idx2];
                g->theta[a][idx0] = 2*g->theta[a][idx1] - g->theta[a][idx2];
                g->theta_vel[a][idx0] = 2*g->theta_vel[a][idx1] - g->theta_vel[a][idx2];
            }
        }
    }
}

/* ================================================================
   Verlet integrator
   ================================================================ */

static void verlet_step(Grid *g, const Config *c) {
    const long N3=g->N3; const double hdt=0.5*g->dt, dt=g->dt;
    for (int a=0;a<NFIELDS;a++) {
        double *vp=g->phi_vel[a],*ap=g->phi_acc[a],*vt=g->theta_vel[a],*at=g->theta_acc[a];
        #pragma omp parallel for schedule(static)
        for (long i=0;i<N3;i++) { vp[i]+=hdt*ap[i]; vt[i]+=hdt*at[i]; }
    }
    for (int a=0;a<NFIELDS;a++) {
        double *pp=g->phi[a],*vp=g->phi_vel[a],*pt=g->theta[a],*vt=g->theta_vel[a];
        #pragma omp parallel for schedule(static)
        for (long i=0;i<N3;i++) { pp[i]+=dt*vp[i]; pt[i]+=dt*vt[i]; }
    }
    compute_forces(g, c);
    for (int a=0;a<NFIELDS;a++) {
        double *vp=g->phi_vel[a],*ap=g->phi_acc[a],*vt=g->theta_vel[a],*at=g->theta_acc[a];
        #pragma omp parallel for schedule(static)
        for (long i=0;i<N3;i++) { vp[i]+=hdt*ap[i]; vt[i]+=hdt*at[i]; }
    }
    if (c->bc_type == 0)
        apply_damping(g, c);
    else if (c->bc_type == 1)
        apply_gradient_bc(g, c);
}

/* ================================================================
   Diagnostics — unified transfer energy decomposition

   E_transfer = V_base(P) * f_transfer(Theta)       (negative, binding)
   E_confine  = |V_base(P)| * gamma * ln(1 + Theta/theta_c)  (positive, theta trap)

   TSV columns reuse existing names:
     E_pot    -> E_transfer
     E_tmass  -> E_confine
   ================================================================ */

static void compute_energy(Grid *g, const Config *c,
    double *epk, double *etk, double *eg, double *em, double *e_transfer,
    double *etg, double *e_confine, double *ec, double *et,
    double *phi_max, double *P_max) {
    const int N=g->N,NN=N*N; const long N3=g->N3;
    const double dx=g->dx, dV=dx*dx*dx, idx1=1.0/(2.0*dx);
    const double MU=c->mu, KAPPA=c->kappa, MASS2=c->m2, MTHETA2=c->mtheta2;
    const double ETA=c->eta;
    const double THETA_C=c->theta_c, EPSILON=c->epsilon, GAMMA=c->gamma_conf;
    const int MODE=c->mode; const double KG=c->kappa_gamma;
    double s_epk=0,s_etk=0,s_eg=0,s_em=0,s_etr=0,s_etg=0,s_eco=0,s_ec=0,s_pm=0,s_Pm=0;

    #pragma omp parallel for reduction(+:s_epk,s_etk,s_eg,s_em,s_etr,s_etg,s_eco,s_ec) \
        reduction(max:s_pm,s_Pm) schedule(static)
    for (long idx=0;idx<N3;idx++) {
        int i=(int)(idx/NN),j=(int)((idx/N)%N),k=(int)(idx%N);
        int ip=(i+1)%N,im=(i-1+N)%N,jp=(j+1)%N,jm=(j-1+N)%N,kp=(k+1)%N,km=(k-1+N)%N;
        long n_ip=(long)ip*NN+j*N+k,n_im=(long)im*NN+j*N+k;
        long n_jp=(long)i*NN+jp*N+k,n_jm=(long)i*NN+jm*N+k;
        long n_kp=(long)i*NN+j*N+kp,n_km=(long)i*NN+j*N+km;
        double p0=g->phi[0][idx],p1=g->phi[1][idx],p2=g->phi[2][idx];
        double sig=p0*p0+p1*p1+p2*p2;
        double me2=(MODE==1)?c->inv_alpha/(1.0+c->inv_beta*sig):MASS2;
        double keff=(MODE==3)?KAPPA/(1.0+KG*sig):KAPPA;
        double P=p0*p1*p2;
        double P2=P*P;

        /* Theta energy density */
        double t0=g->theta[0][idx], t1=g->theta[1][idx], t2=g->theta[2][idx];
        double Theta = t0*t0 + t1*t1 + t2*t2;
        double Theta_sum = Theta + THETA_C;

        /* Transfer function */
        double f_transfer = EPSILON + (1.0-EPSILON) * Theta / Theta_sum;

        for (int a=0;a<NFIELDS;a++) {
            s_epk+=0.5*g->phi_vel[a][idx]*g->phi_vel[a][idx]*dV;
            s_etk+=0.5*g->theta_vel[a][idx]*g->theta_vel[a][idx]*dV;
            double gx=(g->phi[a][n_ip]-g->phi[a][n_im])*idx1;
            double gy=(g->phi[a][n_jp]-g->phi[a][n_jm])*idx1;
            double gz=(g->phi[a][n_kp]-g->phi[a][n_km])*idx1;
            s_eg+=0.5*(gx*gx+gy*gy+gz*gz)*dV;
            s_em+=0.5*me2*g->phi[a][idx]*g->phi[a][idx]*dV;
            double tgx=(g->theta[a][n_ip]-g->theta[a][n_im])*idx1;
            double tgy=(g->theta[a][n_jp]-g->theta[a][n_jm])*idx1;
            double tgz=(g->theta[a][n_kp]-g->theta[a][n_km])*idx1;
            s_etg+=0.5*(tgx*tgx+tgy*tgy+tgz*tgz)*dV;
            /* Bare theta mass energy (if m_theta != 0) included in confine column */
            s_eco+=0.5*MTHETA2*g->theta[a][idx]*g->theta[a][idx]*dV;
            double ap=fabs(g->phi[a][idx]); if(ap>s_pm) s_pm=ap;
        }

        /* E_transfer = V_base(P) * f_transfer  (negative, binding) */
        double V_base;
        if (MODE==3) {
            double D=1.0+KG*sig+KAPPA*P2;
            V_base = (MU/2.0)*P2*(1.0+KG*sig)/D;
        } else {
            V_base = (MU/2.0)*P2/(1.0+keff*P2);
        }
        s_etr += V_base * f_transfer * dV;

        /* E_confine = |V_base| * gamma * ln(1 + Theta/theta_c)  (positive) */
        double V_abs = -V_base;  /* positive */
        s_eco += V_abs * GAMMA * log(1.0 + Theta / THETA_C) * dV;

        double Pa=fabs(P); if(Pa>s_Pm) s_Pm=Pa;

        /* Curl coupling energy */
        for (int a=0;a<NFIELDS;a++) {
            double ct=curl_component(g->theta,a,n_ip,n_im,n_jp,n_jm,n_kp,n_km,idx1);
            s_ec-=ETA*g->phi[a][idx]*ct*dV;
        }
    }
    *epk=s_epk;*etk=s_etk;*eg=s_eg;*em=s_em;*e_transfer=s_etr;
    *etg=s_etg;*e_confine=s_eco;*ec=s_ec;
    *et=s_epk+s_etk+s_eg+s_em+s_etr+s_etg+s_eco+s_ec;
    *phi_max=s_pm;*P_max=s_Pm;
}

static double theta_rms(Grid *g) {
    double sum=0;
    #pragma omp parallel for reduction(+:sum)
    for (long i=0;i<g->N3;i++)
        for (int a=0;a<NFIELDS;a++) sum+=g->theta[a][i]*g->theta[a][i];
    return sqrt(sum/(3.0*g->N3));
}

static double P_integrated(Grid *g) {
    double t=0; const double dV=g->dx*g->dx*g->dx;
    #pragma omp parallel for reduction(+:t)
    for (long i=0;i<g->N3;i++) t+=fabs(g->phi[0][i]*g->phi[1][i]*g->phi[2][i])*dV;
    return t;
}

/* Compute average transfer fraction f(Theta) over grid */
static double f_avg(Grid *g, const Config *c) {
    double sum=0;
    const double THETA_C=c->theta_c, EPSILON=c->epsilon;
    #pragma omp parallel for reduction(+:sum)
    for (long i=0;i<g->N3;i++) {
        double t0=g->theta[0][i], t1=g->theta[1][i], t2=g->theta[2][i];
        double Theta = t0*t0 + t1*t1 + t2*t2;
        sum += EPSILON + (1.0-EPSILON) * Theta / (Theta + THETA_C);
    }
    return sum / (double)g->N3;
}

/* ================================================================
   SFA output: 12 columns, configurable precision
   ================================================================ */

static void *cast_buf = NULL;
static long cast_buf_n = 0;

static void *downcast(double *src, long n, int precision) {
    if (precision == 2) return src;  /* f64: no conversion */
    long need = n * ((precision == 0) ? 2 : 4);
    if (need > cast_buf_n) { cast_buf_n = need; cast_buf = realloc(cast_buf, need); }
    if (precision == 1) { float *p=(float*)cast_buf; for(long i=0;i<n;i++) p[i]=(float)src[i]; }
    else { uint16_t *p=(uint16_t*)cast_buf; for(long i=0;i<n;i++) p[i]=f64_to_f16(src[i]); }
    return cast_buf;
}

static void sfa_snap(SFA *sfa, Grid *g, double t, int precision) {
    long n = g->N3;
    /* 12 arrays: phi[3], theta[3], phi_vel[3], theta_vel[3] */
    double *arrays[12] = {
        g->phi[0], g->phi[1], g->phi[2],
        g->theta[0], g->theta[1], g->theta[2],
        g->phi_vel[0], g->phi_vel[1], g->phi_vel[2],
        g->theta_vel[0], g->theta_vel[1], g->theta_vel[2]
    };

    /* For f64, we can pass pointers directly. For f16/f32, we need to
       downcast each column into separate buffers since sfa_write_frame
       uses all pointers simultaneously. */
    if (precision == 2) {
        sfa_write_frame(sfa, t, (void**)arrays);
    } else {
        void *cols[12];
        int es = (precision == 0) ? 2 : 4;
        for (int c = 0; c < 12; c++) {
            cols[c] = malloc(n * es);
            if (precision == 1) { float *p=(float*)cols[c]; for(long i=0;i<n;i++) p[i]=(float)arrays[c][i]; }
            else { uint16_t *p=(uint16_t*)cols[c]; for(long i=0;i<n;i++) p[i]=f64_to_f16(arrays[c][i]); }
        }
        sfa_write_frame(sfa, t, cols);
        for (int c = 0; c < 12; c++) free(cols[c]);
    }
}

/* ================================================================
   Main
   ================================================================ */

int main(int argc, char **argv) {
    Config c = cfg_defaults();

    /* Load config: from .cfg file, from .sfa file (KVMD), or CLI only */
    if (argc < 2) {
        fprintf(stderr, "Usage: %s config.cfg [-key value ...]\n", argv[0]);
        fprintf(stderr, "       %s input.sfa [-key value ...]   (restart from SFA)\n", argv[0]);
        return 1;
    }

    int arg_start = 2;
    const char *arg1 = argv[1];
    int len1 = strlen(arg1);

    if (len1 > 4 && !strcmp(arg1 + len1 - 4, ".sfa")) {
        /* First arg is an SFA file: read KVMD as base config, set init=sfa */
        printf("Loading parameters from SFA: %s\n", arg1);
        SFA *init_sfa = sfa_open(arg1);
        if (!init_sfa) { fprintf(stderr, "Cannot open SFA: %s\n", arg1); return 1; }

        SFA_KVMDSet kv[SFA_MAX_KVMD_SETS];
        int n_kv = sfa_read_kvmd(init_sfa, kv, SFA_MAX_KVMD_SETS);
        if (n_kv == 0) {
            printf("  WARNING: no KVMD metadata found -- using defaults\n");
        } else {
            /* Apply KVMD set 0 (or the one covering the last frame) as base config */
            SFA_KVMDSet *use = &kv[0];
            for (int i = 0; i < n_kv; i++) {
                if (kv[i].first_frame == 0xFFFFFFFF) { use = &kv[i]; break; }
            }
            printf("  KVMD set %d: %d parameters\n", use->set_id, use->n_pairs);
            for (int i = 0; i < use->n_pairs; i++)
                cfg_set(&c, use->keys[i], use->values[i]);
        }

        /* Override grid from SFA header (authoritative) */
        char tmp[64];
        snprintf(tmp, 64, "%u", init_sfa->Nx); cfg_set(&c, "N", tmp);
        snprintf(tmp, 64, "%.6f", init_sfa->Lx); cfg_set(&c, "L", tmp);

        /* Set init=sfa pointing to this file */
        cfg_set(&c, "init", "sfa");
        strncpy(c.init_sfa, arg1, 511);
        c.init_frame = -1;

        sfa_close(init_sfa);
    } else {
        /* First arg is a config file */
        cfg_load(&c, arg1);
    }

    /* CLI overrides (highest priority) */
    for (int i = arg_start; i < argc - 1; i += 2) {
        const char *key = argv[i];
        if (key[0] == '-') key++;
        cfg_set(&c, key, argv[i+1]);
    }

    /* OpenMP setup */
    int nth = 4;
    char *env_t = getenv("OMP_NUM_THREADS");
    if (env_t) nth = atoi(env_t);
    omp_set_num_threads(nth);

    cfg_print(&c);

    /* Grid */
    Grid *g = grid_alloc(&c);
    printf("dx=%.4f dt=%.6f threads=%d\n\n", g->dx, g->dt, nth);

    /* Initialize */
    do_init(g, &c);
    compute_forces(g, &c);

    /* For gradient_pinned BC: save initial state as pinned boundary values */
    if (c.bc_type == 1) {
        grid_save_pinned(g);
        printf("Gradient BC: pinned %d slabs on each x-face (A_high=%.3f, A_low=%.3f)\n\n",
               c.gradient_margin, c.gradient_A_high, c.gradient_A_low);
    }

    /* SFA archive: 12 columns + KVMD metadata */
    uint8_t sfa_dtype = (c.precision == 0) ? SFA_F16 : (c.precision == 1) ? SFA_F32 : SFA_F64;
    SFA *sfa = sfa_create(c.output, c.N, c.N, c.N, c.L, c.L, c.L, g->dt);
    sfa->flags = SFA_CODEC_COLZSTD | SFA_FLAG_STREAMING;  /* per-column parallel compression */
    if (c.output_split) {
        sfa_enable_split(sfa);
        /* Derive split base: strip .sfa extension */
        strncpy(sfa->split_base, c.output, 507);
        char *ext = strrchr(sfa->split_base, '.');
        if (ext && !strcmp(ext, ".sfa")) *ext = '\0';
    }

    /* Embed physics parameters as KVMD metadata */
    {
        char vN[32],vL[32],vT[32],vdt[32],vm[32],vmt[32],veta[32],vmu[32],vkappa[32];
        char vmode[32],via[32],vib[32],vkg[32],vdw[32],vdr[32],vprec[32],vdelta[64],vbcsw[32];
        char vtc[32],veps[32],vgam[32];
        snprintf(vN,32,"%d",c.N); snprintf(vL,32,"%.6f",c.L); snprintf(vT,32,"%.6f",c.T);
        snprintf(vdt,32,"%.6f",c.dt_factor);
        snprintf(vm,32,"%.6f",sqrt(c.m2)); snprintf(vmt,32,"%.6f",sqrt(c.mtheta2));
        snprintf(veta,32,"%.6f",c.eta);
        snprintf(vmu,32,"%.6f",c.mu);
        snprintf(vkappa,32,"%.6f",c.kappa); snprintf(vmode,32,"%d",c.mode);
        snprintf(via,32,"%.6f",c.inv_alpha); snprintf(vib,32,"%.6f",c.inv_beta);
        snprintf(vkg,32,"%.6f",c.kappa_gamma);
        snprintf(vdw,32,"%.6f",c.damp_width); snprintf(vdr,32,"%.6f",c.damp_rate);
        snprintf(vprec,32,"%s", (const char*[]){"f16","f32","f64"}[c.precision]);
        snprintf(vdelta,64,"%.6f,%.6f,%.6f",c.delta[0],c.delta[1],c.delta[2]);
        snprintf(vbcsw,32,"%.6f",c.bc_switch_time);
        snprintf(vtc,32,"%.6f",c.theta_c);
        snprintf(veps,32,"%.6f",c.epsilon);
        snprintf(vgam,32,"%.6f",c.gamma_conf);
        char vbc[32], vgah[32], vgal[32], vgm[32];
        snprintf(vbc,32,"%d",c.bc_type);
        snprintf(vgah,32,"%.6f",c.gradient_A_high);
        snprintf(vgal,32,"%.6f",c.gradient_A_low);
        snprintf(vgm,32,"%d",c.gradient_margin);
        const char *keys[] = {"N","L","T","dt_factor","m","m_theta","eta","mu","kappa",
                              "theta_c","epsilon","gamma","bc_switch_time",
                              "mode","inv_alpha","inv_beta","kappa_gamma",
                              "damp_width","damp_rate","precision","delta",
                              "bc_type","gradient_A_high","gradient_A_low","gradient_margin"};
        const char *vals[] = {vN,vL,vT,vdt,vm,vmt,veta,vmu,vkappa,
                              vtc,veps,vgam,vbcsw,
                              vmode,via,vib,vkg,vdw,vdr,vprec,vdelta,
                              vbc,vgah,vgal,vgm};
        sfa_add_kvmd(sfa, 0, 0xFFFFFFFF, 0xFFFFFFFF, keys, vals, 25);
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
    const char *pn[] = {"f16","f32","f64"};
    printf("SFA: %s (12 cols, %s, BSS+zstd)\n\n", c.output, pn[c.precision]);

    sfa_snap(sfa, g, 0.0, c.precision);

    /* Diagnostics file — column names updated for unified transfer */
    FILE *fp = fopen(c.diag_file, "w");
    fprintf(fp, "t\tE_phi_kin\tE_theta_kin\tE_grad\tE_mass\tE_transfer\tE_tgrad\tE_confine\t"
                "E_coupling\tE_total\tphi_max\tP_max\tP_int\ttheta_rms\tf_avg\n");

    int n_steps = (int)lround(c.T / g->dt);
    int diag_every = (int)lround(c.diag_dt / g->dt); if (diag_every<1) diag_every=1;
    int snap_every = (int)lround(c.snap_dt / g->dt); if (snap_every<1) snap_every=1;

    /* Initial diagnostic */
    double epk,etk,eg,em,etr,etg,eco,ec,et,pm,Pm;
    compute_energy(g,&c,&epk,&etk,&eg,&em,&etr,&etg,&eco,&ec,&et,&pm,&Pm);
    double Pint0=P_integrated(g), trms0=theta_rms(g), favg0=f_avg(g,&c);
    double E0 = et;
    printf("INIT: E_total=%.4e E_transfer=%.4f E_confine=%.4f phi_max=%.4f P_int=%.4e theta_rms=%.3e f_avg=%.4f\n\n",
           et, etr, eco, pm, Pint0, trms0, favg0);
    fprintf(fp,"%.2f\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6f\n",
            0.0,epk,etk,eg,em,etr,etg,eco,ec,et,pm,Pm,Pint0,trms0,favg0);

    double wall0 = omp_get_wtime();
    int major = diag_every * 25; if (major<1) major=1;
    /* Cap so we get ~10 progress lines even for short runs */
    if (major > n_steps/10) { major = (n_steps/10/diag_every)*diag_every; if (major<diag_every) major=diag_every; }

    int bc_switched = 0;
    int bc_switch_step = (c.bc_switch_time > 0) ? (int)lround(c.bc_switch_time / g->dt) : 0;

    for (int step = 1; step <= n_steps; step++) {
        if (bc_switch_step > 0 && step == bc_switch_step && !bc_switched) {
            c.damp_rate = 0.0;
            printf("\n*** BC SWITCH at t=%.1f: absorbing -> periodic ***\n\n", step * g->dt);
            bc_switched = 1;
        }
        verlet_step(g, &c);
        double t = step * g->dt;

        if (step % snap_every == 0)
            sfa_snap(sfa, g, t, c.precision);

        if (step % diag_every == 0) {
            compute_energy(g,&c,&epk,&etk,&eg,&em,&etr,&etg,&eco,&ec,&et,&pm,&Pm);
            double Pint=P_integrated(g), trms=theta_rms(g), fa=f_avg(g,&c);
            fprintf(fp,"%.2f\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6f\n",
                    t,epk,etk,eg,em,etr,etg,eco,ec,et,pm,Pm,Pint,trms,fa);
            fflush(fp);

            if (step % major == 0) {
                double wall=omp_get_wtime()-wall0;
                double drift=100*(et-E0)/(fabs(E0)+1e-30);
                printf("t=%7.1f E=%.3e (drift %+.3f%%) Etr=%.1f Eco=%.1f phi=%.3f theta_rms=%.2e f_avg=%.4f "
                       "[%.0f%% %.1fs %.2fms/step]\n",
                       t,et,drift,etr,eco,pm,trms,fa,100.0*step/n_steps,wall,1000*wall/step);
                fflush(stdout);
            }
        }
    }

    /* Final frame -- write if the last step wasn't exactly a snap point. */
    {
        int last_snapped = (n_steps / snap_every) * snap_every;
        int gap = n_steps - last_snapped;
        if (gap > 0)
            sfa_snap(sfa, g, n_steps*g->dt, c.precision);
    }
    uint32_t nf = sfa->total_frames;
    sfa_close(sfa);

    /* Final summary */
    compute_energy(g,&c,&epk,&etk,&eg,&em,&etr,&etg,&eco,&ec,&et,&pm,&Pm);
    double trms=theta_rms(g), Pint=P_integrated(g), fa=f_avg(g,&c);
    double wall=omp_get_wtime()-wall0;

    printf("\n=== COMPLETE ===\n");
    printf("E_total=%.4e (drift %.3f%%) E_transfer=%.4f E_confine=%.4f\n",
           et, 100*(et-E0)/(fabs(E0)+1e-30), etr, eco);
    printf("phi_max=%.4f P_int=%.4e theta_rms=%.3e f_avg=%.4f\n", pm, Pint, trms, fa);
    printf("SFA: %s (%u frames)\n", c.output, nf);
    printf("[%s] theta_rms grew: %.2e -> %.2e\n", (trms>trms0+1e-10)?"OK":"WARN", trms0, trms);
    printf("Wall: %.1fs (%.1f min) %.2fms/step\n", wall, wall/60, 1000*wall/n_steps);

    fclose(fp);
    if (cast_buf) free(cast_buf);
    grid_free(g);
    return 0;
}
