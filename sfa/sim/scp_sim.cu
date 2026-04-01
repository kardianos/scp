/*  scp_sim.cu — Unified 6-field Cosserat simulation kernel (CUDA)
 *
 *  GPU version of scp_sim.c — same config, same SFA output, same physics.
 *  Config parsing, SFA I/O, and init run on CPU. Physics runs on GPU.
 *
 *  Build: nvcc -O3 -arch=sm_70 -o scp_sim_cuda scp_sim.cu -lzstd -lm -lpthread
 *  Run:   ./scp_sim_cuda config.cfg [-key value ...]
 */

#define SFA_IMPLEMENTATION
#include "../format/sfa.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>
#include <pthread.h>

#define NFIELDS 3
#define PI 3.14159265358979323846
#define THREADS_PER_BLOCK 256

/* ================================================================
   Configuration — identical to scp_sim.c
   ================================================================ */

typedef struct {
    int N;
    double L, T, dt_factor;
    double m2, mtheta2, eta, mu, kappa;
    double kappa_h;             /* chiral helicity coupling */
    double alpha_cs;            /* Cosserat strain */
    double beta_h;              /* curl-squared hardening */
    int mode;
    double inv_alpha, inv_beta, kappa_gamma;
    int bc_type;                /* 0=absorb_sphere, 1=gradient_pinned */
    double damp_width, damp_rate;
    double bc_switch_time;      /* switch absorb->periodic at this sim time (0=never) */
    double gradient_A_high, gradient_A_low;
    int gradient_margin;
    char init[32];
    double A, sigma, A_bg, ellip, R_tube;
    double delta[3];
    char init_sfa[512];
    int init_frame;
    char init_exec[1024];
    char output[512];
    char diag_file[512];
    int precision;
    int output_split;
    double snap_dt, diag_dt;
    double burst_start, burst_end, burst_every;
} Config;

static Config cfg_defaults(void) {
    Config c = {};
    c.N = 128;  c.L = 10.0;  c.T = 200.0;  c.dt_factor = 0.025;
    c.m2 = 2.25;  c.mtheta2 = 0.0;  c.eta = 0.5;
    c.mu = -41.345;  c.kappa = 50.0;  c.kappa_h = 0.0;  c.alpha_cs = 0.0;  c.beta_h = 0.0;
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
    c.precision = 1;
    c.output_split = 0;
    c.snap_dt = 5.0;  c.diag_dt = 2.0;
    c.burst_start = 0;  c.burst_end = 0;  c.burst_every = 0;
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
    else if (!strcmp(key,"burst_start")) c->burst_start = atof(val);
    else if (!strcmp(key,"burst_end"))   c->burst_end = atof(val);
    else if (!strcmp(key,"burst_every")) c->burst_every = atof(val);
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
    printf("=== scp_sim CUDA: V50/C4 Cosserat + Curl²-Hardening (GPU) ===\n");
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
    }
    else if (c->bc_type == 1)
        printf("\nBC:      gradient pinned (A_high=%.3f A_low=%.3f margin=%d)\n",
               c->gradient_A_high, c->gradient_A_low, c->gradient_margin);
    printf("Init:    %s", c->init);
    if (!strcmp(c->init, "sfa")) printf(" (%s frame=%d)", c->init_sfa, c->init_frame);
    printf("\nOutput:  %s (%s, snap=%.1f diag=%.1f)\n\n",
           c->output, prec_names[c->precision], c->snap_dt, c->diag_dt);
}

/* ================================================================
   f16 helpers (host-side)
   ================================================================ */

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

static inline double f16_to_f64(uint16_t h) {
    uint16_t sign = h & 0x8000;
    int exp = (h >> 10) & 0x1F;
    uint16_t mant = h & 0x3FF;
    if (exp == 0) return sign ? -0.0 : 0.0;
    if (exp == 31) return sign ? -INFINITY : INFINITY;
    float f; uint32_t x = ((uint32_t)sign << 16) | ((uint32_t)(exp-15+127) << 23) | ((uint32_t)mant << 13);
    memcpy(&f, &x, 4); return (double)f;
}

/* ================================================================
   Hook-based async pipeline — FieldState, hooks, contexts
   ================================================================ */

typedef struct {
    double *phi[3], *vel_phi[3], *acc_phi[3];
    double *theta[3], *vel_theta[3], *acc_theta[3];
    long N3;
    int N;
    double L, dx, dt;
} FieldState;

#define MAX_HOOKS 4

typedef void (*HookFn)(int step, double t, const FieldState *state, void *ctx);

typedef struct {
    HookFn fn;
    void *ctx;
} StepHook;

static StepHook hooks[MAX_HOOKS];
static int n_hooks = 0;

static void register_hook(HookFn fn, void *ctx) {
    if (n_hooks < MAX_HOOKS) {
        hooks[n_hooks].fn = fn;
        hooks[n_hooks].ctx = ctx;
        n_hooks++;
    }
}

/* Snapshot hook context — async f16 conversion + DMA + compress/write */
typedef struct {
    SFA *sfa;
    int precision;        /* 0=f16, 1=f32, 2=f64 */
    int snap_every;

    /* GPU resources */
    cudaStream_t stream;
    uint16_t *d_f16_buf;  /* GPU-side f16 staging (12 * N3 * 2 bytes) */
    float *d_f32_buf;     /* GPU-side f32 staging (12 * N3 * 4 bytes) — for precision=1 */

    /* Host resources */
    void *h_pin_buf;      /* PINNED host memory for DMA target */
    long buf_bytes;       /* size of staging buffer */
    long N3;

    /* Async writer thread */
    pthread_t writer_thread;
    pthread_mutex_t mutex;
    pthread_cond_t cond;
    int writer_busy;
    int writer_has_data;
    int writer_shutdown;
    double frame_time;
} SnapHookCtx;

/* Diagnostics hook context — GPU-side reduction */
typedef struct {
    FILE *fp;
    int diag_every;
    int major_every;      /* major print interval */
    long N3;
    double dV;            /* dx^3 */
    double dx;
    double E0;            /* initial energy for drift */
    int n_steps;          /* total steps for % progress */
    cudaEvent_t t_start;  /* for wall-clock timing */

    /* GPU reduction output */
    double *d_results;    /* 12 doubles on GPU */
    double *h_results;    /* 12 doubles on host (pinned) */

    /* Config params for potential energy */
    double m2, mtheta2, eta, mu, kappa;
    double kappa_h, alpha_cs, beta_h;
    int mode;
    double inv_alpha, inv_beta, kappa_gamma;
} DiagHookCtx;

/* ================================================================
   Host grid (for init, diagnostics, SFA I/O)
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
    Grid *g = (Grid*)calloc(1, sizeof(Grid));
    g->N = c->N; g->N3 = (long)c->N * c->N * c->N;
    g->L = c->L; g->dx = 2.0 * c->L / (c->N - 1);
    g->dt = c->dt_factor * g->dx;
    long total = 18 * g->N3;
    g->mem = (double*)malloc(total * sizeof(double));
    if (!g->mem) { fprintf(stderr, "FATAL: host malloc\n"); exit(1); }
    memset(g->mem, 0, total * sizeof(double));
    for (int a = 0; a < NFIELDS; a++) {
        g->phi[a]       = g->mem + (0+a)*g->N3;
        g->phi_vel[a]   = g->mem + (3+a)*g->N3;
        g->phi_acc[a]   = g->mem + (6+a)*g->N3;
        g->theta[a]     = g->mem + (9+a)*g->N3;
        g->theta_vel[a] = g->mem + (12+a)*g->N3;
        g->theta_acc[a] = g->mem + (15+a)*g->N3;
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
   Initialization (host-side, identical to scp_sim.c)
   ================================================================ */

static void init_oscillon(Grid *g, const Config *c) {
    const int N=g->N, NN=N*N; const double L=g->L, dx=g->dx;
    printf("Init: oscillon (A=%.3f sigma=%.3f)\n", c->A, c->sigma);
    for (int i=0;i<N;i++) { double x=-L+i*dx;
    for (int j=0;j<N;j++) { double y=-L+j*dx;
    for (int k=0;k<N;k++) { double z=-L+k*dx;
        long idx=(long)i*NN+j*N+k;
        double r2=x*x+y*y+z*z;
        double env=c->A*exp(-r2/(2.0*c->sigma*c->sigma));
        for (int a=0;a<NFIELDS;a++) g->phi[a][idx]=env*cos(c->delta[a]);
    }}}
}

static void init_braid(Grid *g, const Config *c) {
    const int N=g->N, NN=N*N; const double L=g->L, dx=g->dx;
    const double kw=PI/L, omega=sqrt(kw*kw+c->m2);
    const double sx=1+c->ellip, sy=1-c->ellip;
    const double inv2R2=1.0/(2*c->R_tube*c->R_tube);
    const double k_bg=PI/L, omega_bg=sqrt(k_bg*k_bg+c->m2);
    printf("Init: braid (R=%.1f ellip=%.4f A=%.2f A_bg=%.2f)\n",
           c->R_tube, c->ellip, c->A, c->A_bg);
    for (int i=0;i<N;i++) { double x=-L+i*dx;
    for (int j=0;j<N;j++) { double y=-L+j*dx;
    for (int k=0;k<N;k++) { double z=-L+k*dx;
        long idx=(long)i*NN+j*N+k;
        double r2e=x*x/(sx*sx)+y*y/(sy*sy);
        double env=exp(-r2e*inv2R2);
        for (int a=0;a<NFIELDS;a++) {
            double ph=kw*z+c->delta[a], ph_bg=k_bg*z+2*PI*a/3.0;
            g->phi[a][idx]=c->A*env*cos(ph)+c->A_bg*cos(ph_bg);
            g->phi_vel[a][idx]=omega*c->A*env*sin(ph)+omega_bg*c->A_bg*sin(ph_bg);
        }
    }}}
}

static void init_from_sfa(Grid *g, const Config *c) {
    printf("Init: SFA file '%s' frame=%d\n", c->init_sfa, c->init_frame);
    SFA *sfa = sfa_open(c->init_sfa);
    if (!sfa) { fprintf(stderr, "FATAL: cannot open SFA '%s'\n", c->init_sfa); exit(1); }
    if ((int)sfa->Nx != g->N || (int)sfa->Ny != g->N || (int)sfa->Nz != g->N) {
        fprintf(stderr, "FATAL: SFA grid %ux%ux%u != sim grid %d^3\n",
                sfa->Nx, sfa->Ny, sfa->Nz, g->N); exit(1);
    }
    int frame = c->init_frame;
    if (frame < 0) frame = sfa->total_frames + frame;
    printf("  Grid: %ux%ux%u, %d columns, %u frames, reading frame %d\n",
           sfa->Nx, sfa->Ny, sfa->Nz, sfa->n_columns, sfa->total_frames, frame);
    void *buf = malloc(sfa->frame_bytes);
    sfa_read_frame(sfa, frame, buf);
    int loaded[12] = {0};
    uint64_t off = 0;
    for (int col = 0; col < sfa->n_columns; col++) {
        int dtype=sfa->columns[col].dtype, sem=sfa->columns[col].semantic, comp=sfa->columns[col].component;
        int es=sfa_dtype_size[dtype]; uint8_t *src=(uint8_t*)buf+off;
        double *target=NULL; int slot=-1;
        if (sem==SFA_POSITION&&comp<3) {target=g->phi[comp];slot=comp;}
        else if (sem==SFA_ANGLE&&comp<3) {target=g->theta[comp];slot=3+comp;}
        else if (sem==SFA_VELOCITY&&comp<3) {target=g->phi_vel[comp];slot=6+comp;}
        else if (sem==SFA_VELOCITY&&comp>=3&&comp<6) {target=g->theta_vel[comp-3];slot=9+comp-3;}
        if (target) {
            long N3=g->N3;
            if (dtype==SFA_F64) for(long i=0;i<N3;i++) target[i]=((double*)src)[i];
            else if (dtype==SFA_F32) for(long i=0;i<N3;i++) target[i]=(double)((float*)src)[i];
            else if (dtype==SFA_F16) for(long i=0;i<N3;i++) target[i]=f16_to_f64(((uint16_t*)src)[i]);
            loaded[slot]=1;
        }
        off += (uint64_t)g->N3 * es;
    }
    int nf=0,nv=0; for(int i=0;i<6;i++) nf+=loaded[i]; for(int i=6;i<12;i++) nv+=loaded[i];
    printf("  Loaded: %d/6 fields, %d/6 velocities\n", nf, nv);
    if (nv==0) printf("  WARNING: no velocity data — cold restart\n");
    free(buf); sfa_close(sfa);
}

/* init=template: Load a small template SFA and stamp it into the grid.
 * Background generated analytically (with gradient if bc_type=1).
 * Template is centered at (cx, cy, cz) in world coords. */
static void init_template(Grid *g, const Config *c) {
    printf("Init: template '%s' at (%.1f, %.1f, %.1f)\n",
           c->init_sfa, 0.0, 0.0, 0.0);  /* center for now */

    /* Step 1: Generate background (with optional gradient) */
    const int N=g->N, NN=N*N;
    const double L=g->L, dx=g->dx;
    const double k_bg=PI/L, omega_bg=sqrt(k_bg*k_bg+c->m2);

    for (int i=0; i<N; i++) { double x=-L+i*dx;
        double A_bg_x = c->A_bg;
        if (c->bc_type == 1) {
            double frac = (x + L) / (2.0 * L);
            A_bg_x = c->gradient_A_high * (1.0 - frac) + c->gradient_A_low * frac;
        }
        for (int j=0; j<N; j++) {
        for (int k=0; k<N; k++) { double z=-L+k*dx;
            long idx = (long)i*NN + j*N + k;
            for (int a=0; a<NFIELDS; a++) {
                double ph_bg = k_bg*z + 2*PI*a/3.0;
                g->phi[a][idx] = A_bg_x * cos(ph_bg);
                g->phi_vel[a][idx] = omega_bg * A_bg_x * sin(ph_bg);
                g->theta[a][idx] = 0;
                g->theta_vel[a][idx] = 0;
            }
        }}
    }

    /* Step 2: Load template SFA (small, e.g. 64^3) */
    SFA *tmpl = sfa_open(c->init_sfa);
    if (!tmpl) { fprintf(stderr, "FATAL: cannot open template '%s'\n", c->init_sfa); exit(1); }
    int TN = tmpl->Nx;
    double TL = tmpl->Lx;
    double Tdx = 2.0 * TL / (TN - 1);
    long TN3 = (long)TN * TN * TN;
    int TNN = TN * TN;

    printf("  Template: %dx%dx%d, L=%.2f, dx=%.4f\n", TN, TN, TN, TL, Tdx);

    /* Read last frame */
    int frame = tmpl->total_frames - 1;
    void *buf = malloc(tmpl->frame_bytes);
    sfa_read_frame(tmpl, frame, buf);

    /* Extract template fields */
    double *tphi[3]={0}, *tvel[3]={0}, *ttheta[3]={0}, *ttvel[3]={0};
    for (int a=0; a<3; a++) {
        tphi[a] = (double*)calloc(TN3, sizeof(double));
        tvel[a] = (double*)calloc(TN3, sizeof(double));
        ttheta[a] = (double*)calloc(TN3, sizeof(double));
        ttvel[a] = (double*)calloc(TN3, sizeof(double));
    }

    uint64_t off = 0;
    for (int col = 0; col < tmpl->n_columns; col++) {
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

    /* Step 3: Subtract template's background, add perturbation to main grid */
    /* Template center maps to (cx, cy, cz) in the main grid (default: 0,0,0) */
    double cx=0, cy=0, cz=0;
    int ci = (int)((cx + L) / dx + 0.5);
    int cj = (int)((cy + L) / dx + 0.5);
    int ck = (int)((cz + L) / dx + 0.5);
    int half = TN / 2;

    int placed = 0;
    for (int ti=0; ti<TN; ti++) {
        int gi = ci + ti - half;
        if (gi < 0 || gi >= N) continue;
        double tx = -TL + ti * Tdx;
        double gx = -L + gi * dx;

        /* Template background at this x */
        double T_A_bg = c->A_bg;  /* template was generated with uniform A_bg */

        for (int tj=0; tj<TN; tj++) {
            int gj = cj + tj - half;
            if (gj < 0 || gj >= N) continue;
            for (int tk=0; tk<TN; tk++) {
                int gk = ck + tk - half;
                if (gk < 0 || gk >= N) continue;

                long tidx = (long)ti*TNN + tj*TN + tk;
                long gidx = (long)gi*NN + gj*N + gk;
                double tz = -TL + tk * Tdx;
                double gz = -L + gk * dx;

                for (int a=0; a<NFIELDS; a++) {
                    /* Template background value */
                    double ph_bg_t = k_bg * tz + 2*PI*a/3.0;  /* use template z for bg subtraction */
                    double bg_phi = T_A_bg * cos(ph_bg_t);
                    double bg_vel = omega_bg * T_A_bg * sin(ph_bg_t);

                    /* Perturbation = template - background */
                    double dphi = tphi[a][tidx] - bg_phi;
                    double dvel = tvel[a][tidx] - bg_vel;

                    /* Add perturbation to main grid */
                    g->phi[a][gidx] += dphi;
                    g->phi_vel[a][gidx] += dvel;
                    g->theta[a][gidx] += ttheta[a][tidx];
                    g->theta_vel[a][gidx] += ttvel[a][tidx];
                }
                placed++;
            }
        }
    }

    for (int a=0; a<3; a++) { free(tphi[a]); free(tvel[a]); free(ttheta[a]); free(ttvel[a]); }
    printf("  Placed %d voxels from template (%d^3 = %ld)\n", placed, TN, TN3);
}

static void do_init(Grid *g, const Config *c) {
    if      (!strcmp(c->init, "oscillon")) init_oscillon(g, c);
    else if (!strcmp(c->init, "braid"))    init_braid(g, c);
    else if (!strcmp(c->init, "sfa"))      init_from_sfa(g, c);
    else if (!strcmp(c->init, "template")) init_template(g, c);
    else { fprintf(stderr, "FATAL: unknown init '%s'\n", c->init); exit(1); }
}

/* ================================================================
   GPU constant memory
   ================================================================ */

__constant__ double d_MU, d_KAPPA, d_MASS2, d_MTHETA2, d_ETA;
__constant__ double d_KAPPA_H, d_ALPHA_CS, d_BETA_H;
__constant__ double d_inv_alpha, d_inv_beta, d_kappa_gamma;
__constant__ double d_idx2, d_idx1;
__constant__ double d_L, d_dx, d_DAMP_WIDTH, d_DAMP_RATE;
__constant__ int d_N, d_NN, d_MODE;
__constant__ long d_N3;
__constant__ int d_BC_TYPE, d_GRAD_MARGIN;
__constant__ double d_GRAD_A_HIGH, d_GRAD_A_LOW;

/* ================================================================
   GPU kernels
   ================================================================ */

/* Pass 1: compute intermediate fields (mismatch M, hardening Q) at all voxels.
 * Only launched when alpha_cs > 0 or beta_h > 0. */
__global__ void compute_intermediates_kernel(
    const double *phi0, const double *phi1, const double *phi2,
    const double *theta0, const double *theta1, const double *theta2,
    double *mis0, double *mis1, double *mis2,
    double *Q0, double *Q1, double *Q2)
{
    long idx = (long)blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= d_N3) return;

    int N = d_N, NN = d_NN;
    int i = (int)(idx / NN), j = (int)((idx / N) % N), k = (int)(idx % N);
    int ip=(i+1)%N, im=(i-1+N)%N, jp=(j+1)%N, jm=(j-1+N)%N, kp=(k+1)%N, km=(k-1+N)%N;
    long n_ip=(long)ip*NN+j*N+k, n_im=(long)im*NN+j*N+k;
    long n_jp=(long)i*NN+jp*N+k, n_jm=(long)i*NN+jm*N+k;
    long n_kp=(long)i*NN+j*N+kp, n_km=(long)i*NN+j*N+km;

    const double *phi[3] = {phi0, phi1, phi2};
    /* curl(phi) at this voxel */
    double cp0 = (phi[2][n_jp]-phi[2][n_jm]-phi[1][n_kp]+phi[1][n_km])*d_idx1;
    double cp1 = (phi[0][n_kp]-phi[0][n_km]-phi[2][n_ip]+phi[2][n_im])*d_idx1;
    double cp2 = (phi[1][n_ip]-phi[1][n_im]-phi[0][n_jp]+phi[0][n_jm])*d_idx1;

    /* Cosserat mismatch: M_a = curl(φ)_a/2 - θ_a */
    if (d_ALPHA_CS != 0) {
        mis0[idx] = cp0 * 0.5 - theta0[idx];
        mis1[idx] = cp1 * 0.5 - theta1[idx];
        mis2[idx] = cp2 * 0.5 - theta2[idx];
    }

    /* Hardening vector: Q_a = (β/2)|θ|²·curl(φ)_a */
    if (d_BETA_H != 0) {
        double t0=theta0[idx], t1=theta1[idx], t2=theta2[idx];
        double T2 = t0*t0 + t1*t1 + t2*t2;
        double coeff = 0.5 * d_BETA_H * T2;
        Q0[idx] = coeff * cp0;
        Q1[idx] = coeff * cp1;
        Q2[idx] = coeff * cp2;
    }
}

/* Pass 2: compute all forces using pre-computed intermediates */
__global__ void compute_forces_kernel(
    const double *phi0, const double *phi1, const double *phi2,
    const double *theta0, const double *theta1, const double *theta2,
    double *acc_phi0, double *acc_phi1, double *acc_phi2,
    double *acc_theta0, double *acc_theta1, double *acc_theta2,
    const double *mis0, const double *mis1, const double *mis2,
    const double *Q0, const double *Q1, const double *Q2)
{
    long idx = (long)blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= d_N3) return;

    int N = d_N, NN = d_NN;
    int i = (int)(idx / NN), j = (int)((idx / N) % N), k = (int)(idx % N);
    int ip=(i+1)%N, im=(i-1+N)%N, jp=(j+1)%N, jm=(j-1+N)%N, kp=(k+1)%N, km=(k-1+N)%N;
    long n_ip=(long)ip*NN+j*N+k, n_im=(long)im*NN+j*N+k;
    long n_jp=(long)i*NN+jp*N+k, n_jm=(long)i*NN+jm*N+k;
    long n_kp=(long)i*NN+j*N+kp, n_km=(long)i*NN+j*N+km;

    /* Load phi at center */
    double p0=phi0[idx], p1=phi1[idx], p2=phi2[idx];
    double sig = p0*p0 + p1*p1 + p2*p2;

    /* Mode-dependent effective mass and kappa */
    double me2 = (d_MODE==1) ? d_inv_alpha/(1.0+d_inv_beta*sig) : d_MASS2;
    double keff = (d_MODE==3) ? d_KAPPA/(1.0+d_kappa_gamma*sig) : d_KAPPA;

    double P = p0*p1*p2, P2 = P*P;
    double den = 1.0 + keff*P2;
    double dVdP = d_MU * P / (den*den);

    /* Mode 3 Term II */
    double t2c = 0;
    if (d_MODE == 3 && d_kappa_gamma > 0) {
        double D = 1.0 + d_kappa_gamma*sig + d_KAPPA*P2;
        t2c = d_MU * d_kappa_gamma * d_KAPPA * P2 * P2 / (D*D);
    }

    /* Phi field arrays for curl */
    const double *phi[3] = {phi0, phi1, phi2};
    const double *th[3] = {theta0, theta1, theta2};

    /* curl(φ) at this voxel — needed for chiral and theta hardening */
    double curl_phi[3];
    curl_phi[0] = (phi[2][n_jp]-phi[2][n_jm]-phi[1][n_kp]+phi[1][n_km])*d_idx1;
    curl_phi[1] = (phi[0][n_kp]-phi[0][n_km]-phi[2][n_ip]+phi[2][n_im])*d_idx1;
    curl_phi[2] = (phi[1][n_ip]-phi[1][n_im]-phi[0][n_jp]+phi[0][n_jm])*d_idx1;

    /* Chiral helicity: κ_h P² φ·curl(φ) */
    double h_chiral = p0*curl_phi[0] + p1*curl_phi[1] + p2*curl_phi[2];
    double dP2dx=0, dP2dy=0, dP2dz=0;
    if (d_KAPPA_H != 0) {
        double Pip=phi0[n_ip]*phi1[n_ip]*phi2[n_ip], Pim=phi0[n_im]*phi1[n_im]*phi2[n_im];
        double Pjp=phi0[n_jp]*phi1[n_jp]*phi2[n_jp], Pjm=phi0[n_jm]*phi1[n_jm]*phi2[n_jm];
        double Pkp=phi0[n_kp]*phi1[n_kp]*phi2[n_kp], Pkm=phi0[n_km]*phi1[n_km]*phi2[n_km];
        dP2dx=(Pip*Pip-Pim*Pim)*d_idx1; dP2dy=(Pjp*Pjp-Pjm*Pjm)*d_idx1; dP2dz=(Pkp*Pkp-Pkm*Pkm)*d_idx1;
    }

    /* Cosserat strain: curl(M) for phi force */
    double curl_M[3] = {0, 0, 0};
    if (d_ALPHA_CS != 0) {
        const double *M[3] = {mis0, mis1, mis2};
        curl_M[0] = (M[2][n_jp]-M[2][n_jm]-M[1][n_kp]+M[1][n_km])*d_idx1;
        curl_M[1] = (M[0][n_kp]-M[0][n_km]-M[2][n_ip]+M[2][n_im])*d_idx1;
        curl_M[2] = (M[1][n_ip]-M[1][n_im]-M[0][n_jp]+M[0][n_jm])*d_idx1;
    }

    /* Curl-squared hardening: curl(Q) for phi force */
    double curl_Q[3] = {0, 0, 0};
    if (d_BETA_H != 0) {
        const double *Qv[3] = {Q0, Q1, Q2};
        curl_Q[0] = (Qv[2][n_jp]-Qv[2][n_jm]-Qv[1][n_kp]+Qv[1][n_km])*d_idx1;
        curl_Q[1] = (Qv[0][n_kp]-Qv[0][n_km]-Qv[2][n_ip]+Qv[2][n_im])*d_idx1;
        curl_Q[2] = (Qv[1][n_ip]-Qv[1][n_im]-Qv[0][n_jp]+Qv[0][n_jm])*d_idx1;
    }

    /* |curl(φ)|² for theta hardening */
    double curl_sq = curl_phi[0]*curl_phi[0] + curl_phi[1]*curl_phi[1] + curl_phi[2]*curl_phi[2];

    double acc_p[3], acc_t[3];

    /* --- Phi forces --- */
    for (int a = 0; a < 3; a++) {
        double lap = (phi[a][n_ip]+phi[a][n_im]+phi[a][n_jp]+phi[a][n_jm]
                     +phi[a][n_kp]+phi[a][n_km]-6.0*phi[a][idx]) * d_idx2;
        double dPda = (a==0)?p1*p2:(a==1)?p0*p2:p0*p1;

        /* curl(theta)_a */
        double ct;
        if (a==0) ct = (th[2][n_jp]-th[2][n_jm]-th[1][n_kp]+th[1][n_km])*d_idx1;
        else if (a==1) ct = (th[0][n_kp]-th[0][n_km]-th[2][n_ip]+th[2][n_im])*d_idx1;
        else ct = (th[1][n_ip]-th[1][n_im]-th[0][n_jp]+th[0][n_jm])*d_idx1;

        /* Chiral force */
        double chiral = 0;
        if (d_KAPPA_H != 0) {
            double t1c = 2.0 * P * dPda * h_chiral;
            double t2c_ = 2.0 * P2 * curl_phi[a];
            double t3c;
            if (a == 0)      t3c = -(dP2dz * p1 - dP2dy * p2);
            else if (a == 1) t3c = -(dP2dx * p2 - dP2dz * p0);
            else             t3c = -(dP2dy * p0 - dP2dx * p1);
            chiral = d_KAPPA_H * (t1c + t2c_ + t3c);
        }

        acc_p[a] = lap - me2*phi[a][idx] - dVdP*dPda - t2c*phi[a][idx]
                 + d_ETA*ct + chiral
                 - d_ALPHA_CS * curl_M[a]    /* Cosserat strain */
                 - 2.0 * curl_Q[a];          /* Curl² hardening */
    }

    /* --- Theta forces --- */
    for (int a = 0; a < 3; a++) {
        double lapt = (th[a][n_ip]+th[a][n_im]+th[a][n_jp]+th[a][n_jm]
                      +th[a][n_kp]+th[a][n_km]-6.0*th[a][idx]) * d_idx2;
        double cp;
        if (a==0) cp = (phi[2][n_jp]-phi[2][n_jm]-phi[1][n_kp]+phi[1][n_km])*d_idx1;
        else if (a==1) cp = (phi[0][n_kp]-phi[0][n_km]-phi[2][n_ip]+phi[2][n_im])*d_idx1;
        else cp = (phi[1][n_ip]-phi[1][n_im]-phi[0][n_jp]+phi[0][n_jm])*d_idx1;

        /* Cosserat strain: +2α·M_a (computed inline from M at center) */
        double cosserat_theta = 0;
        if (d_ALPHA_CS != 0) {
            const double *M[3] = {mis0, mis1, mis2};
            cosserat_theta = 2.0 * d_ALPHA_CS * M[a][idx];
        }

        /* Hardening: -β|∇×φ|²·θ_a */
        double harden_theta = -d_BETA_H * curl_sq * th[a][idx];

        acc_t[a] = lapt - d_MTHETA2*th[a][idx] + d_ETA*cp
                 + cosserat_theta + harden_theta;
    }

    acc_phi0[idx]=acc_p[0]; acc_phi1[idx]=acc_p[1]; acc_phi2[idx]=acc_p[2];
    acc_theta0[idx]=acc_t[0]; acc_theta1[idx]=acc_t[1]; acc_theta2[idx]=acc_t[2];
}

__global__ void verlet_halfkick_kernel(
    double *vp0, double *vp1, double *vp2,
    double *vt0, double *vt1, double *vt2,
    const double *ap0, const double *ap1, const double *ap2,
    const double *at0, const double *at1, const double *at2,
    double hdt)
{
    long idx = (long)blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= d_N3) return;
    vp0[idx]+=hdt*ap0[idx]; vp1[idx]+=hdt*ap1[idx]; vp2[idx]+=hdt*ap2[idx];
    vt0[idx]+=hdt*at0[idx]; vt1[idx]+=hdt*at1[idx]; vt2[idx]+=hdt*at2[idx];
}

__global__ void verlet_drift_kernel(
    double *p0, double *p1, double *p2,
    double *t0, double *t1, double *t2,
    const double *vp0, const double *vp1, const double *vp2,
    const double *vt0, const double *vt1, const double *vt2,
    double dt)
{
    long idx = (long)blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= d_N3) return;
    p0[idx]+=dt*vp0[idx]; p1[idx]+=dt*vp1[idx]; p2[idx]+=dt*vp2[idx];
    t0[idx]+=dt*vt0[idx]; t1[idx]+=dt*vt1[idx]; t2[idx]+=dt*vt2[idx];
}

__global__ void absorbing_boundary_kernel(
    double *vp0, double *vp1, double *vp2,
    double *vt0, double *vt1, double *vt2)
{
    long idx = (long)blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= d_N3) return;
    if (d_DAMP_WIDTH <= 0 || d_DAMP_RATE <= 0) return;

    int N=d_N, NN=d_NN;
    int i=(int)(idx/NN), j=(int)((idx/N)%N), k=(int)(idx%N);
    double x = -d_L + i*d_dx, y = -d_L + j*d_dx, z = -d_L + k*d_dx;
    double r = sqrt(x*x + y*y + z*z);
    double R_damp = d_L - d_DAMP_WIDTH;

    if (r > R_damp) {
        double s = (r - R_damp) / d_DAMP_WIDTH;
        if (s > 1.0) s = 1.0;
        double damp = 1.0 - d_DAMP_RATE * s * s;
        vp0[idx]*=damp; vp1[idx]*=damp; vp2[idx]*=damp;
        vt0[idx]*=damp; vt1[idx]*=damp; vt2[idx]*=damp;
    }
}

/* bc_type=1: GPU gradient pinned boundary kernel
 * Computes boundary values ANALYTICALLY from the gradient parameters.
 * x-boundary: phi_a = A_bg(x) * cos(k*z + 2*pi*a/3), theta=0
 * y-boundary: linear extrapolation from interior
 * No pinned arrays needed — zero extra GPU memory. */
__global__ void gradient_bc_kernel(
    double *p0, double *p1, double *p2,
    double *vp0, double *vp1, double *vp2,
    double *ap0, double *ap1, double *ap2,
    double *t0, double *t1, double *t2,
    double *vt0, double *vt1, double *vt2,
    double *at0, double *at1, double *at2)
{
    long idx = (long)blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= d_N3) return;

    int N = d_N, NN = d_NN;
    int M = d_GRAD_MARGIN;
    int i = (int)(idx / NN), j = (int)((idx / N) % N), k = (int)(idx % N);

    /* x-direction: pin boundary slabs to analytical gradient background */
    if (i < M || i >= N - M) {
        double x = -d_L + i * d_dx;
        double z = -d_L + k * d_dx;
        double frac = (x + d_L) / (2.0 * d_L);
        double A_bg = d_GRAD_A_HIGH * (1.0 - frac) + d_GRAD_A_LOW * frac;
        double k_bg = 3.14159265358979323846 / d_L;
        double omega_bg = sqrt(k_bg * k_bg + d_MASS2);
        double delta[3] = {0.0, 2.0943951023931953, 4.1887902047863905}; /* 0, 2pi/3, 4pi/3 */

        p0[idx] = A_bg * cos(k_bg * z + delta[0]);
        p1[idx] = A_bg * cos(k_bg * z + delta[1]);
        p2[idx] = A_bg * cos(k_bg * z + delta[2]);
        vp0[idx] = omega_bg * A_bg * sin(k_bg * z + delta[0]);
        vp1[idx] = omega_bg * A_bg * sin(k_bg * z + delta[1]);
        vp2[idx] = omega_bg * A_bg * sin(k_bg * z + delta[2]);
        ap0[idx] = 0; ap1[idx] = 0; ap2[idx] = 0;
        t0[idx] = 0; t1[idx] = 0; t2[idx] = 0;
        vt0[idx] = 0; vt1[idx] = 0; vt2[idx] = 0;
        at0[idx] = 0; at1[idx] = 0; at2[idx] = 0;
        return;
    }

    /* y-direction: linear extrapolation at j=0 and j=N-1 */
    if (j == 0) {
        long idx1 = (long)i*NN + 1*N + k;
        long idx2 = (long)i*NN + 2*N + k;
        p0[idx]=2*p0[idx1]-p0[idx2]; p1[idx]=2*p1[idx1]-p1[idx2]; p2[idx]=2*p2[idx1]-p2[idx2];
        vp0[idx]=2*vp0[idx1]-vp0[idx2]; vp1[idx]=2*vp1[idx1]-vp1[idx2]; vp2[idx]=2*vp2[idx1]-vp2[idx2];
        t0[idx]=2*t0[idx1]-t0[idx2]; t1[idx]=2*t1[idx1]-t1[idx2]; t2[idx]=2*t2[idx1]-t2[idx2];
        vt0[idx]=2*vt0[idx1]-vt0[idx2]; vt1[idx]=2*vt1[idx1]-vt1[idx2]; vt2[idx]=2*vt2[idx1]-vt2[idx2];
    } else if (j == N - 1) {
        long idx1 = (long)i*NN + (N-2)*N + k;
        long idx2 = (long)i*NN + (N-3)*N + k;
        p0[idx]=2*p0[idx1]-p0[idx2]; p1[idx]=2*p1[idx1]-p1[idx2]; p2[idx]=2*p2[idx1]-p2[idx2];
        vp0[idx]=2*vp0[idx1]-vp0[idx2]; vp1[idx]=2*vp1[idx1]-vp1[idx2]; vp2[idx]=2*vp2[idx1]-vp2[idx2];
        t0[idx]=2*t0[idx1]-t0[idx2]; t1[idx]=2*t1[idx1]-t1[idx2]; t2[idx]=2*t2[idx1]-t2[idx2];
        vt0[idx]=2*vt0[idx1]-vt0[idx2]; vt1[idx]=2*vt1[idx1]-vt1[idx2]; vt2[idx]=2*vt2[idx1]-vt2[idx2];
    }
}

__global__ void downcast_f64_to_f32_kernel(const double *src, float *dst, long n) {
    long idx = (long)blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n) return;
    dst[idx] = (float)src[idx];
}

__global__ void downcast_f64_to_f16_kernel(const double *src, uint16_t *dst, long n) {
    long idx = (long)blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n) return;
    float f = (float)src[idx];
    uint32_t x = __float_as_uint(f);
    uint16_t sign = (x >> 16) & 0x8000;
    int exp = ((x >> 23) & 0xFF) - 127 + 15;
    uint16_t mant = (x >> 13) & 0x3FF;
    if (exp <= 0) { dst[idx] = sign; return; }
    if (exp >= 31) { dst[idx] = sign | 0x7C00; return; }
    dst[idx] = sign | (exp << 10) | mant;
}

/* ================================================================
   GPU diagnostics reduction kernel
   Computes 12 values in a single pass via block-level reduction:
   [0] E_phi_kin, [1] E_theta_kin, [2] E_grad, [3] E_mass,
   [4] E_pot, [5] E_tgrad, [6] E_tmass, [7] E_coupling,
   [8] phi_max, [9] P_max, [10] P_int, [11] theta_rms_sum
   ================================================================ */
#define DIAG_NVALS 12

__global__ void reduce_diagnostics_kernel(
    const double *p0, const double *p1, const double *p2,
    const double *t0, const double *t1, const double *t2,
    const double *vp0, const double *vp1, const double *vp2,
    const double *vt0, const double *vt1, const double *vt2,
    double dV, double idx1, double idx2_val,
    double mass2, double mtheta2, double eta, double mu, double kappa,
    double alpha_cs, double beta_h, double kappa_h,
    int mode, double inv_alpha, double inv_beta, double kappa_gamma,
    double *d_results)
{
    __shared__ double sdata[DIAG_NVALS * 256];  /* THREADS_PER_BLOCK * DIAG_NVALS */

    long idx = (long)blockIdx.x * blockDim.x + threadIdx.x;
    int tid = threadIdx.x;

    double local[DIAG_NVALS];
    for (int v = 0; v < DIAG_NVALS; v++) local[v] = 0.0;

    if (idx < d_N3) {
        int N = d_N, NN = d_NN;
        int i = (int)(idx / NN), j = (int)((idx / N) % N), k = (int)(idx % N);
        int ip=(i+1)%N, im=(i-1+N)%N, jp=(j+1)%N, jm=(j-1+N)%N, kp=(k+1)%N, km=(k-1+N)%N;
        long n_ip=(long)ip*NN+j*N+k, n_im=(long)im*NN+j*N+k;
        long n_jp=(long)i*NN+jp*N+k, n_jm=(long)i*NN+jm*N+k;
        long n_kp=(long)i*NN+j*N+kp, n_km=(long)i*NN+j*N+km;

        double pp0=p0[idx], pp1=p1[idx], pp2=p2[idx];
        double sig = pp0*pp0 + pp1*pp1 + pp2*pp2;
        double me2 = (mode==1) ? inv_alpha/(1.0+inv_beta*sig) : mass2;
        double keff = (mode==3) ? kappa/(1.0+kappa_gamma*sig) : kappa;

        const double *phi[3] = {p0, p1, p2};
        const double *th[3] = {t0, t1, t2};
        const double *vphi[3] = {vp0, vp1, vp2};
        const double *vth[3] = {vt0, vt1, vt2};

        double P = pp0*pp1*pp2, P2 = P*P;

        double s_epk=0, s_etk=0, s_eg=0, s_em=0, s_etg=0, s_etm=0, s_ec=0;
        double s_pm=0, s_trms=0;

        for (int a = 0; a < 3; a++) {
            s_epk += 0.5*vphi[a][idx]*vphi[a][idx]*dV;
            s_etk += 0.5*vth[a][idx]*vth[a][idx]*dV;
            double gx=(phi[a][n_ip]-phi[a][n_im])*idx1;
            double gy=(phi[a][n_jp]-phi[a][n_jm])*idx1;
            double gz=(phi[a][n_kp]-phi[a][n_km])*idx1;
            s_eg += 0.5*(gx*gx+gy*gy+gz*gz)*dV;
            s_em += 0.5*me2*phi[a][idx]*phi[a][idx]*dV;
            double tgx=(th[a][n_ip]-th[a][n_im])*idx1;
            double tgy=(th[a][n_jp]-th[a][n_jm])*idx1;
            double tgz=(th[a][n_kp]-th[a][n_km])*idx1;
            s_etg += 0.5*(tgx*tgx+tgy*tgy+tgz*tgz)*dV;
            s_etm += 0.5*mtheta2*th[a][idx]*th[a][idx]*dV;
            double ap = fabs(phi[a][idx]); if (ap > s_pm) s_pm = ap;
            s_trms += th[a][idx]*th[a][idx];
        }

        double s_ep;
        if (mode==3) {
            double D = 1.0+kappa_gamma*sig+kappa*P2;
            s_ep = (mu/2.0)*P2*(1.0+kappa_gamma*sig)/D*dV;
        } else {
            s_ep = (mu/2.0)*P2/(1.0+keff*P2)*dV;
        }

        /* Cosserat/hardening/chiral energy contributions (added to E_pot) */
        double cx0=(phi[2][n_jp]-phi[2][n_jm]-phi[1][n_kp]+phi[1][n_km])*idx1;
        double cx1=(phi[0][n_kp]-phi[0][n_km]-phi[2][n_ip]+phi[2][n_im])*idx1;
        double cx2=(phi[1][n_ip]-phi[1][n_im]-phi[0][n_jp]+phi[0][n_jm])*idx1;
        if (alpha_cs != 0) {
            double m0=cx0*0.5-th[0][idx], m1=cx1*0.5-th[1][idx], m2c=cx2*0.5-th[2][idx];
            s_ep += alpha_cs * (m0*m0 + m1*m1 + m2c*m2c) * dV;
        }
        if (beta_h != 0) {
            double csq = cx0*cx0 + cx1*cx1 + cx2*cx2;
            double T2 = th[0][idx]*th[0][idx]+th[1][idx]*th[1][idx]+th[2][idx]*th[2][idx];
            s_ep += 0.5 * beta_h * T2 * csq * dV;
        }
        if (kappa_h != 0) {
            double hel = pp0*cx0 + pp1*cx1 + pp2*cx2;
            s_ep += kappa_h * P2 * hel * dV;
        }

        /* Curl coupling energy */
        double curl_th[3];
        curl_th[0]=(th[2][n_jp]-th[2][n_jm]-th[1][n_kp]+th[1][n_km])*idx1;
        curl_th[1]=(th[0][n_kp]-th[0][n_km]-th[2][n_ip]+th[2][n_im])*idx1;
        curl_th[2]=(th[1][n_ip]-th[1][n_im]-th[0][n_jp]+th[0][n_jm])*idx1;
        for (int a=0;a<3;a++) s_ec -= eta * phi[a][idx] * curl_th[a] * dV;

        local[0] = s_epk;
        local[1] = s_etk;
        local[2] = s_eg;
        local[3] = s_em;
        local[4] = s_ep;
        local[5] = s_etg;
        local[6] = s_etm;
        local[7] = s_ec;
        local[8] = s_pm;        /* phi_max — will use max reduction */
        local[9] = fabs(P);     /* P_max — will use max reduction */
        local[10] = fabs(P)*dV; /* P_int */
        local[11] = s_trms;     /* theta_rms_sum (sum of squares) */
    }

    /* Store to shared memory */
    for (int v = 0; v < DIAG_NVALS; v++)
        sdata[v * blockDim.x + tid] = local[v];
    __syncthreads();

    /* Block-level reduction */
    for (int s = blockDim.x/2; s > 0; s >>= 1) {
        if (tid < s) {
            for (int v = 0; v < DIAG_NVALS; v++) {
                int si = v * blockDim.x;
                if (v == 8 || v == 9) {
                    /* max reduction for phi_max, P_max */
                    double other = sdata[si + tid + s];
                    if (other > sdata[si + tid]) sdata[si + tid] = other;
                } else {
                    sdata[si + tid] += sdata[si + tid + s];
                }
            }
        }
        __syncthreads();
    }

    /* Write block results to global memory via atomicAdd / atomicMax */
    if (tid == 0) {
        for (int v = 0; v < DIAG_NVALS; v++) {
            if (v == 8 || v == 9) {
                /* Atomic max for doubles — use unsigned long long CAS loop */
                unsigned long long *addr = (unsigned long long *)&d_results[v];
                unsigned long long old_val = *addr, assumed;
                double new_val = sdata[v * blockDim.x];
                do {
                    assumed = old_val;
                    double cur = __longlong_as_double(assumed);
                    if (new_val <= cur) break;
                    unsigned long long new_bits = __double_as_longlong(new_val);
                    old_val = atomicCAS(addr, assumed, new_bits);
                } while (assumed != old_val);
            } else {
                atomicAdd(&d_results[v], sdata[v * blockDim.x]);
            }
        }
    }
}

/* Zero the diagnostics results buffer */
__global__ void zero_diag_kernel(double *d_results, int n) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < n) d_results[idx] = 0.0;
}

/* ================================================================
   GPU memory management
   ================================================================ */

static double *d_phi[3], *d_vel_phi[3], *d_acc_phi[3];
static double *d_theta[3], *d_vel_theta[3], *d_acc_theta[3];
static double *d_mismatch[3], *d_harden_Q[3];  /* intermediates for Cosserat/hardening */
static int gpu_blocks;
static int gpu_has_intermediates = 0;  /* nonzero if intermediate arrays allocated */

static void gpu_alloc(long N3, int need_intermediates) {
    size_t bytes = N3 * sizeof(double);
    for (int a = 0; a < 3; a++) {
        cudaMalloc(&d_phi[a], bytes);       cudaMalloc(&d_vel_phi[a], bytes);   cudaMalloc(&d_acc_phi[a], bytes);
        cudaMalloc(&d_theta[a], bytes);     cudaMalloc(&d_vel_theta[a], bytes); cudaMalloc(&d_acc_theta[a], bytes);
    }
    double total_gb = 18.0*bytes/1e9;
    if (need_intermediates) {
        for (int a = 0; a < 3; a++) {
            cudaMalloc(&d_mismatch[a], bytes);
            cudaMalloc(&d_harden_Q[a], bytes);
            cudaMemset(d_mismatch[a], 0, bytes);
            cudaMemset(d_harden_Q[a], 0, bytes);
        }
        gpu_has_intermediates = 1;
        total_gb += 6.0*bytes/1e9;
    }
    gpu_blocks = (int)((N3 + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK);
    printf("GPU: allocated %.2f GB (physics%s), %d blocks × %d threads\n",
           total_gb, need_intermediates ? "+Cosserat" : "", gpu_blocks, THREADS_PER_BLOCK);
}

static void gpu_upload(Grid *g) {
    size_t bytes = g->N3 * sizeof(double);
    for (int a = 0; a < 3; a++) {
        cudaMemcpy(d_phi[a], g->phi[a], bytes, cudaMemcpyHostToDevice);
        cudaMemcpy(d_vel_phi[a], g->phi_vel[a], bytes, cudaMemcpyHostToDevice);
        cudaMemcpy(d_theta[a], g->theta[a], bytes, cudaMemcpyHostToDevice);
        cudaMemcpy(d_vel_theta[a], g->theta_vel[a], bytes, cudaMemcpyHostToDevice);
        cudaMemset(d_acc_phi[a], 0, bytes);
        cudaMemset(d_acc_theta[a], 0, bytes);
    }
}

static void gpu_download(Grid *g) {
    size_t bytes = g->N3 * sizeof(double);
    for (int a = 0; a < 3; a++) {
        cudaMemcpy(g->phi[a], d_phi[a], bytes, cudaMemcpyDeviceToHost);
        cudaMemcpy(g->phi_vel[a], d_vel_phi[a], bytes, cudaMemcpyDeviceToHost);
        cudaMemcpy(g->theta[a], d_theta[a], bytes, cudaMemcpyDeviceToHost);
        cudaMemcpy(g->theta_vel[a], d_vel_theta[a], bytes, cudaMemcpyDeviceToHost);
    }
}

static void gpu_free(void) {
    for (int a = 0; a < 3; a++) {
        cudaFree(d_phi[a]); cudaFree(d_vel_phi[a]); cudaFree(d_acc_phi[a]);
        cudaFree(d_theta[a]); cudaFree(d_vel_theta[a]); cudaFree(d_acc_theta[a]);
    }
    if (gpu_has_intermediates) {
        for (int a = 0; a < 3; a++) {
            cudaFree(d_mismatch[a]); cudaFree(d_harden_Q[a]);
        }
    }
}

static void gpu_set_constants(const Config *c, double dx) {
    double idx2 = 1.0/(dx*dx), idx1 = 1.0/(2.0*dx);
    int N = c->N, NN = N*N;
    long N3 = (long)N*N*N;
    cudaMemcpyToSymbol(d_MU, &c->mu, sizeof(double));
    cudaMemcpyToSymbol(d_KAPPA, &c->kappa, sizeof(double));
    cudaMemcpyToSymbol(d_MASS2, &c->m2, sizeof(double));
    cudaMemcpyToSymbol(d_MTHETA2, &c->mtheta2, sizeof(double));
    cudaMemcpyToSymbol(d_ETA, &c->eta, sizeof(double));
    cudaMemcpyToSymbol(d_KAPPA_H, &c->kappa_h, sizeof(double));
    cudaMemcpyToSymbol(d_ALPHA_CS, &c->alpha_cs, sizeof(double));
    cudaMemcpyToSymbol(d_BETA_H, &c->beta_h, sizeof(double));
    cudaMemcpyToSymbol(d_inv_alpha, &c->inv_alpha, sizeof(double));
    cudaMemcpyToSymbol(d_inv_beta, &c->inv_beta, sizeof(double));
    cudaMemcpyToSymbol(d_kappa_gamma, &c->kappa_gamma, sizeof(double));
    cudaMemcpyToSymbol(d_idx2, &idx2, sizeof(double));
    cudaMemcpyToSymbol(d_idx1, &idx1, sizeof(double));
    cudaMemcpyToSymbol(d_L, &c->L, sizeof(double));
    cudaMemcpyToSymbol(d_dx, &dx, sizeof(double));
    cudaMemcpyToSymbol(d_DAMP_WIDTH, &c->damp_width, sizeof(double));
    cudaMemcpyToSymbol(d_DAMP_RATE, &c->damp_rate, sizeof(double));
    cudaMemcpyToSymbol(d_N, &N, sizeof(int));
    cudaMemcpyToSymbol(d_NN, &NN, sizeof(int));
    cudaMemcpyToSymbol(d_MODE, &c->mode, sizeof(int));
    cudaMemcpyToSymbol(d_N3, &N3, sizeof(long));
}

/* GPU Verlet step */
static void gpu_verlet_step(double dt) {
    double hdt = 0.5 * dt;
    /* Half-kick */
    verlet_halfkick_kernel<<<gpu_blocks, THREADS_PER_BLOCK>>>(
        d_vel_phi[0], d_vel_phi[1], d_vel_phi[2],
        d_vel_theta[0], d_vel_theta[1], d_vel_theta[2],
        d_acc_phi[0], d_acc_phi[1], d_acc_phi[2],
        d_acc_theta[0], d_acc_theta[1], d_acc_theta[2], hdt);
    /* Drift */
    verlet_drift_kernel<<<gpu_blocks, THREADS_PER_BLOCK>>>(
        d_phi[0], d_phi[1], d_phi[2],
        d_theta[0], d_theta[1], d_theta[2],
        d_vel_phi[0], d_vel_phi[1], d_vel_phi[2],
        d_vel_theta[0], d_vel_theta[1], d_vel_theta[2], dt);
    /* Intermediates (Cosserat mismatch + hardening Q) */
    if (gpu_has_intermediates) {
        compute_intermediates_kernel<<<gpu_blocks, THREADS_PER_BLOCK>>>(
            d_phi[0], d_phi[1], d_phi[2],
            d_theta[0], d_theta[1], d_theta[2],
            d_mismatch[0], d_mismatch[1], d_mismatch[2],
            d_harden_Q[0], d_harden_Q[1], d_harden_Q[2]);
    }
    /* Forces */
    compute_forces_kernel<<<gpu_blocks, THREADS_PER_BLOCK>>>(
        d_phi[0], d_phi[1], d_phi[2],
        d_theta[0], d_theta[1], d_theta[2],
        d_acc_phi[0], d_acc_phi[1], d_acc_phi[2],
        d_acc_theta[0], d_acc_theta[1], d_acc_theta[2],
        d_mismatch[0], d_mismatch[1], d_mismatch[2],
        d_harden_Q[0], d_harden_Q[1], d_harden_Q[2]);
    /* Half-kick */
    verlet_halfkick_kernel<<<gpu_blocks, THREADS_PER_BLOCK>>>(
        d_vel_phi[0], d_vel_phi[1], d_vel_phi[2],
        d_vel_theta[0], d_vel_theta[1], d_vel_theta[2],
        d_acc_phi[0], d_acc_phi[1], d_acc_phi[2],
        d_acc_theta[0], d_acc_theta[1], d_acc_theta[2], hdt);
    /* Boundary — type set externally before calling gpu_verlet_step */
    absorbing_boundary_kernel<<<gpu_blocks, THREADS_PER_BLOCK>>>(
        d_vel_phi[0], d_vel_phi[1], d_vel_phi[2],
        d_vel_theta[0], d_vel_theta[1], d_vel_theta[2]);
}

/* Save current host grid state as pinned boundary values (for bc_type=1) */
static void grid_save_pinned(Grid *g) {
    for (int a = 0; a < NFIELDS; a++) {
        g->pin_phi[a]   = (double*)malloc(g->N3 * sizeof(double));
        g->pin_vel[a]   = (double*)malloc(g->N3 * sizeof(double));
        g->pin_theta[a] = (double*)malloc(g->N3 * sizeof(double));
        g->pin_tvel[a]  = (double*)malloc(g->N3 * sizeof(double));
        memcpy(g->pin_phi[a],   g->phi[a],       g->N3 * sizeof(double));
        memcpy(g->pin_vel[a],   g->phi_vel[a],   g->N3 * sizeof(double));
        memcpy(g->pin_theta[a], g->theta[a],     g->N3 * sizeof(double));
        memcpy(g->pin_tvel[a],  g->theta_vel[a], g->N3 * sizeof(double));
    }
}

/* Apply gradient pinned BC on GPU: download boundary slabs, overwrite, upload.
 * This is done host-side because pinned data lives on the host. The boundary
 * slabs are small (margin/N fraction of the grid), so the transfer cost is low
 * relative to the GPU compute time for large grids. */
static void gpu_apply_gradient_bc(Grid *g, const Config *c) {
    const int N = g->N, NN = N*N;
    const int M = c->gradient_margin;
    const size_t bytes = g->N3 * sizeof(double);

    /* Download full fields from GPU to host */
    for (int a = 0; a < NFIELDS; a++) {
        cudaMemcpy(g->phi[a],       d_phi[a],       bytes, cudaMemcpyDeviceToHost);
        cudaMemcpy(g->phi_vel[a],   d_vel_phi[a],   bytes, cudaMemcpyDeviceToHost);
        cudaMemcpy(g->theta[a],     d_theta[a],     bytes, cudaMemcpyDeviceToHost);
        cudaMemcpy(g->theta_vel[a], d_vel_theta[a], bytes, cudaMemcpyDeviceToHost);
    }

    /* x-direction: pin boundary slabs */
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

    /* y-direction: linear extrapolation from interior */
    for (int i = M; i < N - M; i++) {
        for (int k = 0; k < N; k++) {
            for (int a = 0; a < NFIELDS; a++) {
                long idx0 = (long)i*NN + 0*N + k;
                long idx1 = (long)i*NN + 1*N + k;
                long idx2 = (long)i*NN + 2*N + k;
                g->phi[a][idx0] = 2*g->phi[a][idx1] - g->phi[a][idx2];
                g->phi_vel[a][idx0] = 2*g->phi_vel[a][idx1] - g->phi_vel[a][idx2];
                g->theta[a][idx0] = 2*g->theta[a][idx1] - g->theta[a][idx2];
                g->theta_vel[a][idx0] = 2*g->theta_vel[a][idx1] - g->theta_vel[a][idx2];

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

    /* Upload modified fields back to GPU */
    for (int a = 0; a < NFIELDS; a++) {
        cudaMemcpy(d_phi[a],       g->phi[a],       bytes, cudaMemcpyHostToDevice);
        cudaMemcpy(d_vel_phi[a],   g->phi_vel[a],   bytes, cudaMemcpyHostToDevice);
        cudaMemcpy(d_theta[a],     g->theta[a],     bytes, cudaMemcpyHostToDevice);
        cudaMemcpy(d_vel_theta[a], g->theta_vel[a], bytes, cudaMemcpyHostToDevice);
    }
}

/* ================================================================
   Snapshot hook — async f16/f32 conversion + DMA + compress/write
   ================================================================ */

static void *snap_writer_thread(void *arg) {
    SnapHookCtx *ctx = (SnapHookCtx *)arg;
    while (1) {
        pthread_mutex_lock(&ctx->mutex);
        while (!ctx->writer_has_data && !ctx->writer_shutdown)
            pthread_cond_wait(&ctx->cond, &ctx->mutex);
        if (ctx->writer_shutdown) { pthread_mutex_unlock(&ctx->mutex); break; }
        ctx->writer_has_data = 0;
        ctx->writer_busy = 1;
        pthread_mutex_unlock(&ctx->mutex);

        /* Compress and write frame to SFA */
        void *cols[12];
        long N3 = ctx->N3;
        if (ctx->precision == 0) {
            /* f16: h_pin_buf is uint16_t[12*N3] */
            uint16_t *base = (uint16_t *)ctx->h_pin_buf;
            for (int c = 0; c < 12; c++)
                cols[c] = base + c * N3;
        } else if (ctx->precision == 1) {
            /* f32: h_pin_buf is float[12*N3] */
            float *base = (float *)ctx->h_pin_buf;
            for (int c = 0; c < 12; c++)
                cols[c] = base + c * N3;
        } else {
            /* f64: h_pin_buf is double[12*N3] */
            double *base = (double *)ctx->h_pin_buf;
            for (int c = 0; c < 12; c++)
                cols[c] = base + c * N3;
        }
        sfa_write_frame(ctx->sfa, ctx->frame_time, cols);

        pthread_mutex_lock(&ctx->mutex);
        ctx->writer_busy = 0;
        pthread_cond_signal(&ctx->cond);
        pthread_mutex_unlock(&ctx->mutex);
    }
    return NULL;
}

static SnapHookCtx create_snap_hook(SFA *sfa, int precision, int snap_every, long N3) {
    SnapHookCtx ctx;
    memset(&ctx, 0, sizeof(ctx));
    ctx.sfa = sfa;
    ctx.precision = precision;
    ctx.snap_every = snap_every;
    ctx.N3 = N3;

    cudaStreamCreate(&ctx.stream);

    /* Allocate GPU staging buffer and pinned host buffer */
    if (precision == 0) {
        ctx.buf_bytes = 12L * N3 * sizeof(uint16_t);
        cudaMalloc(&ctx.d_f16_buf, ctx.buf_bytes);
        ctx.d_f32_buf = NULL;
    } else if (precision == 1) {
        ctx.buf_bytes = 12L * N3 * sizeof(float);
        ctx.d_f16_buf = NULL;
        cudaMalloc(&ctx.d_f32_buf, ctx.buf_bytes);
    } else {
        /* f64: DMA directly from physics arrays, staging = pinned host only */
        ctx.buf_bytes = 12L * N3 * sizeof(double);
        ctx.d_f16_buf = NULL;
        ctx.d_f32_buf = NULL;
    }
    cudaMallocHost(&ctx.h_pin_buf, ctx.buf_bytes);

    printf("Snap hook: %.2f GB GPU staging + %.2f GB pinned host\n",
           (precision < 2 ? ctx.buf_bytes : 0) / 1e9, ctx.buf_bytes / 1e9);

    pthread_mutex_init(&ctx.mutex, NULL);
    pthread_cond_init(&ctx.cond, NULL);
    ctx.writer_busy = 0;
    ctx.writer_has_data = 0;
    ctx.writer_shutdown = 0;

    /* NOTE: writer thread is NOT started here — must call start_snap_writer()
     * after the ctx is in its final memory location (returned by value). */
    return ctx;
}

static void start_snap_writer(SnapHookCtx *ctx) {
    pthread_create(&ctx->writer_thread, NULL, snap_writer_thread, ctx);
}

static void destroy_snap_hook(SnapHookCtx *ctx) {
    /* Wait for any pending write AND any queued data to be consumed.
     * Signal the cond to wake the writer if it's waiting for data —
     * otherwise both threads deadlock on the same cond_wait. */
    pthread_mutex_lock(&ctx->mutex);
    while (ctx->writer_busy || ctx->writer_has_data) {
        pthread_cond_signal(&ctx->cond);  /* wake writer to consume data */
        pthread_cond_wait(&ctx->cond, &ctx->mutex);
    }
    ctx->writer_shutdown = 1;
    pthread_cond_signal(&ctx->cond);
    pthread_mutex_unlock(&ctx->mutex);
    pthread_join(ctx->writer_thread, NULL);

    pthread_mutex_destroy(&ctx->mutex);
    pthread_cond_destroy(&ctx->cond);

    cudaStreamDestroy(ctx->stream);
    if (ctx->d_f16_buf) cudaFree(ctx->d_f16_buf);
    if (ctx->d_f32_buf) cudaFree(ctx->d_f32_buf);
    cudaFreeHost(ctx->h_pin_buf);
}

static void snap_hook(int step, double t, const FieldState *state, void *vctx) {
    SnapHookCtx *ctx = (SnapHookCtx *)vctx;
    if (step % ctx->snap_every != 0) return;

    /* Wait for previous write to finish (the only possible stall) */
    pthread_mutex_lock(&ctx->mutex);
    while (ctx->writer_busy)
        pthread_cond_wait(&ctx->cond, &ctx->mutex);
    pthread_mutex_unlock(&ctx->mutex);

    long N3 = state->N3;
    int blocks = (int)((N3 + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK);

    const double *src[12] = {
        state->phi[0], state->phi[1], state->phi[2],
        state->theta[0], state->theta[1], state->theta[2],
        state->vel_phi[0], state->vel_phi[1], state->vel_phi[2],
        state->vel_theta[0], state->vel_theta[1], state->vel_theta[2]
    };

    if (ctx->precision == 0) {
        /* f16 path: convert on GPU, DMA f16 */
        for (int c = 0; c < 12; c++) {
            downcast_f64_to_f16_kernel<<<blocks, THREADS_PER_BLOCK, 0, ctx->stream>>>(
                src[c], ctx->d_f16_buf + c * N3, N3);
        }
        cudaStreamSynchronize(ctx->stream);
        cudaMemcpyAsync(ctx->h_pin_buf, ctx->d_f16_buf, ctx->buf_bytes,
                        cudaMemcpyDeviceToHost, ctx->stream);
    } else if (ctx->precision == 1) {
        /* f32 path: convert on GPU, DMA f32 */
        for (int c = 0; c < 12; c++) {
            downcast_f64_to_f32_kernel<<<blocks, THREADS_PER_BLOCK, 0, ctx->stream>>>(
                src[c], ctx->d_f32_buf + c * N3, N3);
        }
        cudaStreamSynchronize(ctx->stream);
        cudaMemcpyAsync(ctx->h_pin_buf, ctx->d_f32_buf, ctx->buf_bytes,
                        cudaMemcpyDeviceToHost, ctx->stream);
    } else {
        /* f64 path: DMA raw doubles (12 separate copies) */
        double *dst = (double *)ctx->h_pin_buf;
        for (int c = 0; c < 12; c++) {
            cudaMemcpyAsync(dst + c * N3, src[c], N3 * sizeof(double),
                            cudaMemcpyDeviceToHost, ctx->stream);
        }
    }

    /* Wait for DMA to complete before signaling writer */
    cudaStreamSynchronize(ctx->stream);

    /* Signal the writer thread */
    ctx->frame_time = t;
    pthread_mutex_lock(&ctx->mutex);
    ctx->writer_has_data = 1;
    pthread_cond_signal(&ctx->cond);
    pthread_mutex_unlock(&ctx->mutex);
    /* Return immediately — physics continues on stream 0 */
}

/* ================================================================
   Diagnostics hook — GPU-side reduction, no full download
   ================================================================ */

static DiagHookCtx create_diag_hook(FILE *fp, int diag_every, int major_every,
                                     long N3, double dx, int n_steps,
                                     cudaEvent_t t_start, const Config *c) {
    DiagHookCtx ctx;
    memset(&ctx, 0, sizeof(ctx));
    ctx.fp = fp;
    ctx.diag_every = diag_every;
    ctx.major_every = major_every;
    ctx.N3 = N3;
    ctx.dV = dx * dx * dx;
    ctx.dx = dx;
    ctx.E0 = 0;
    ctx.n_steps = n_steps;
    ctx.t_start = t_start;
    ctx.m2 = c->m2;
    ctx.mtheta2 = c->mtheta2;
    ctx.eta = c->eta;
    ctx.mu = c->mu;
    ctx.kappa = c->kappa;
    ctx.kappa_h = c->kappa_h;
    ctx.alpha_cs = c->alpha_cs;
    ctx.beta_h = c->beta_h;
    ctx.mode = c->mode;
    ctx.inv_alpha = c->inv_alpha;
    ctx.inv_beta = c->inv_beta;
    ctx.kappa_gamma = c->kappa_gamma;

    cudaMalloc(&ctx.d_results, DIAG_NVALS * sizeof(double));
    cudaMallocHost(&ctx.h_results, DIAG_NVALS * sizeof(double));
    return ctx;
}

static void destroy_diag_hook(DiagHookCtx *ctx) {
    cudaFree(ctx->d_results);
    cudaFreeHost(ctx->h_results);
}

/* Run the GPU reduction and populate h_results. Used by diag_hook and for initial/final diag. */
static void run_gpu_diagnostics(const FieldState *state, DiagHookCtx *ctx) {
    /* Zero the results buffer */
    zero_diag_kernel<<<1, DIAG_NVALS>>>(ctx->d_results, DIAG_NVALS);

    double idx1 = 1.0 / (2.0 * ctx->dx);
    double idx2_val = 1.0 / (ctx->dx * ctx->dx);

    reduce_diagnostics_kernel<<<gpu_blocks, THREADS_PER_BLOCK>>>(
        state->phi[0], state->phi[1], state->phi[2],
        state->theta[0], state->theta[1], state->theta[2],
        state->vel_phi[0], state->vel_phi[1], state->vel_phi[2],
        state->vel_theta[0], state->vel_theta[1], state->vel_theta[2],
        ctx->dV, idx1, idx2_val,
        ctx->m2, ctx->mtheta2, ctx->eta, ctx->mu, ctx->kappa,
        ctx->alpha_cs, ctx->beta_h, ctx->kappa_h,
        ctx->mode, ctx->inv_alpha, ctx->inv_beta, ctx->kappa_gamma,
        ctx->d_results);

    cudaMemcpy(ctx->h_results, ctx->d_results, DIAG_NVALS * sizeof(double),
               cudaMemcpyDeviceToHost);
}

static void diag_hook(int step, double t, const FieldState *state, void *vctx) {
    DiagHookCtx *ctx = (DiagHookCtx *)vctx;
    if (step % ctx->diag_every != 0) return;

    run_gpu_diagnostics(state, ctx);

    double *r = ctx->h_results;
    /* r: [0]epk [1]etk [2]eg [3]em [4]ep [5]etg [6]etm [7]ec [8]pm [9]Pm [10]Pint [11]trms_sum */
    double et = r[0]+r[1]+r[2]+r[3]+r[4]+r[5]+r[6]+r[7];
    double trms = sqrt(r[11] / (3.0 * ctx->N3));

    fprintf(ctx->fp, "%.2f\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\n",
            t, r[0], r[1], r[2], r[3], r[4], r[5], r[6], r[7], et, r[8], r[9], r[10], trms);
    fflush(ctx->fp);

    if (step % ctx->major_every == 0) {
        cudaEvent_t t_stop;
        cudaEventCreate(&t_stop);
        cudaEventRecord(t_stop);
        cudaEventSynchronize(t_stop);
        float ms; cudaEventElapsedTime(&ms, ctx->t_start, t_stop);
        cudaEventDestroy(t_stop);
        double drift = 100.0 * (et - ctx->E0) / (fabs(ctx->E0) + 1e-30);
        printf("t=%7.1f E=%.3e (drift %+.3f%%) Ep=%.1f phi=%.3f θ_rms=%.2e "
               "[%.0f%% %.1fs %.2fms/step]\n",
               t, et, drift, r[4], r[8], trms,
               100.0 * step / ctx->n_steps, ms / 1000, ms / step);
        fflush(stdout);
    }
}

/* ================================================================
   Diagnostics (host-side, same as CPU version) — kept for reference/fallback
   ================================================================ */

static void compute_energy(Grid *g, const Config *c,
    double *epk, double *etk, double *eg, double *em, double *ep,
    double *etg, double *etm, double *ec, double *et,
    double *phi_max, double *P_max) {
    const int N=g->N,NN=N*N; const long N3=g->N3;
    const double dx=g->dx, dV=dx*dx*dx, idx1=1.0/(2.0*dx);
    double s_epk=0,s_etk=0,s_eg=0,s_em=0,s_ep=0,s_etg=0,s_etm=0,s_ec=0,s_pm=0,s_Pm=0;
    for (long idx=0;idx<N3;idx++) {
        int i=(int)(idx/NN),j=(int)((idx/N)%N),k=(int)(idx%N);
        int ip=(i+1)%N,im=(i-1+N)%N,jp=(j+1)%N,jm=(j-1+N)%N,kp=(k+1)%N,km=(k-1+N)%N;
        long n_ip=(long)ip*NN+j*N+k,n_im=(long)im*NN+j*N+k;
        long n_jp=(long)i*NN+jp*N+k,n_jm=(long)i*NN+jm*N+k;
        long n_kp=(long)i*NN+j*N+kp,n_km=(long)i*NN+j*N+km;
        double p0=g->phi[0][idx],p1=g->phi[1][idx],p2=g->phi[2][idx];
        double sig=p0*p0+p1*p1+p2*p2;
        double me2=(c->mode==1)?c->inv_alpha/(1.0+c->inv_beta*sig):c->m2;
        double keff=(c->mode==3)?c->kappa/(1.0+c->kappa_gamma*sig):c->kappa;
        double P=p0*p1*p2, P2=P*P;
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
            s_etm+=0.5*c->mtheta2*g->theta[a][idx]*g->theta[a][idx]*dV;
            double ap=fabs(g->phi[a][idx]); if(ap>s_pm) s_pm=ap;
        }
        if (c->mode==3) {
            double D=1.0+c->kappa_gamma*sig+c->kappa*P2;
            s_ep+=(c->mu/2.0)*P2*(1.0+c->kappa_gamma*sig)/D*dV;
        } else s_ep+=(c->mu/2.0)*P2/(1.0+keff*P2)*dV;
        double Pa=fabs(P); if(Pa>s_Pm) s_Pm=Pa;
        /* Cosserat/hardening/chiral energy */
        double cx0=(g->phi[2][n_jp]-g->phi[2][n_jm]-g->phi[1][n_kp]+g->phi[1][n_km])*idx1;
        double cx1=(g->phi[0][n_kp]-g->phi[0][n_km]-g->phi[2][n_ip]+g->phi[2][n_im])*idx1;
        double cx2=(g->phi[1][n_ip]-g->phi[1][n_im]-g->phi[0][n_jp]+g->phi[0][n_jm])*idx1;
        if (c->alpha_cs != 0) {
            double m0=cx0*0.5-g->theta[0][idx], m1=cx1*0.5-g->theta[1][idx], m2c=cx2*0.5-g->theta[2][idx];
            s_ep += c->alpha_cs * (m0*m0 + m1*m1 + m2c*m2c) * dV;
        }
        if (c->beta_h != 0) {
            double csq = cx0*cx0 + cx1*cx1 + cx2*cx2;
            double T2 = g->theta[0][idx]*g->theta[0][idx]+g->theta[1][idx]*g->theta[1][idx]+g->theta[2][idx]*g->theta[2][idx];
            s_ep += 0.5 * c->beta_h * T2 * csq * dV;
        }
        if (c->kappa_h != 0) {
            double hel = p0*cx0 + p1*cx1 + p2*cx2;
            s_ep += c->kappa_h * P2 * hel * dV;
        }
        /* Curl coupling energy */
        double curl_th[3];
        curl_th[0]=(g->theta[2][n_jp]-g->theta[2][n_jm]-g->theta[1][n_kp]+g->theta[1][n_km])*idx1;
        curl_th[1]=(g->theta[0][n_kp]-g->theta[0][n_km]-g->theta[2][n_ip]+g->theta[2][n_im])*idx1;
        curl_th[2]=(g->theta[1][n_ip]-g->theta[1][n_im]-g->theta[0][n_jp]+g->theta[0][n_jm])*idx1;
        for (int a=0;a<3;a++) s_ec -= c->eta * g->phi[a][idx] * curl_th[a] * dV;
    }
    *epk=s_epk;*etk=s_etk;*eg=s_eg;*em=s_em;*ep=s_ep;
    *etg=s_etg;*etm=s_etm;*ec=s_ec;
    *et=s_epk+s_etk+s_eg+s_em+s_ep+s_etg+s_etm+s_ec;
    *phi_max=s_pm;*P_max=s_Pm;
}

static double theta_rms(Grid *g) {
    double sum=0;
    for (long i=0;i<g->N3;i++)
        for (int a=0;a<NFIELDS;a++) sum+=g->theta[a][i]*g->theta[a][i];
    return sqrt(sum/(3.0*g->N3));
}

static double P_integrated(Grid *g) {
    double t=0; const double dV=g->dx*g->dx*g->dx;
    for (long i=0;i<g->N3;i++) t+=fabs(g->phi[0][i]*g->phi[1][i]*g->phi[2][i])*dV;
    return t;
}

/* ================================================================
   SFA output (host-side, using GPU downcast for f32)
   ================================================================ */

static void sfa_snap_gpu(SFA *sfa, Grid *g, double t, int precision) {
    long n = g->N3;
    /* Download all fields+velocities from GPU */
    gpu_download(g);

    if (precision == 2) { /* f64 */
        double *arrays[12] = {
            g->phi[0],g->phi[1],g->phi[2],g->theta[0],g->theta[1],g->theta[2],
            g->phi_vel[0],g->phi_vel[1],g->phi_vel[2],
            g->theta_vel[0],g->theta_vel[1],g->theta_vel[2]
        };
        sfa_write_frame(sfa, t, (void**)arrays);
    } else {
        void *cols[12];
        double *arrays[12] = {
            g->phi[0],g->phi[1],g->phi[2],g->theta[0],g->theta[1],g->theta[2],
            g->phi_vel[0],g->phi_vel[1],g->phi_vel[2],
            g->theta_vel[0],g->theta_vel[1],g->theta_vel[2]
        };
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

    if (argc < 2) {
        fprintf(stderr, "Usage: %s config.cfg [-key value ...]\n", argv[0]);
        fprintf(stderr, "       %s input.sfa [-key value ...]   (restart from SFA)\n", argv[0]);
        return 1;
    }

    const char *arg1 = argv[1];
    int len1 = strlen(arg1);
    if (len1 > 4 && !strcmp(arg1 + len1 - 4, ".sfa")) {
        printf("Loading parameters from SFA: %s\n", arg1);
        SFA *init_sfa = sfa_open(arg1);
        if (!init_sfa) { fprintf(stderr, "Cannot open SFA: %s\n", arg1); return 1; }
        SFA_KVMDSet kv[SFA_MAX_KVMD_SETS];
        int n_kv = sfa_read_kvmd(init_sfa, kv, SFA_MAX_KVMD_SETS);
        if (n_kv > 0) {
            SFA_KVMDSet *use = &kv[0];
            for (int i=0;i<n_kv;i++) if(kv[i].first_frame==0xFFFFFFFF){use=&kv[i];break;}
            printf("  KVMD set %d: %d parameters\n", use->set_id, use->n_pairs);
            for (int i=0;i<use->n_pairs;i++) cfg_set(&c, use->keys[i], use->values[i]);
        }
        char tmp[64];
        snprintf(tmp,64,"%u",init_sfa->Nx); cfg_set(&c,"N",tmp);
        snprintf(tmp,64,"%.6f",init_sfa->Lx); cfg_set(&c,"L",tmp);
        cfg_set(&c,"init","sfa");
        strncpy(c.init_sfa, arg1, 511);
        c.init_frame = -1;
        sfa_close(init_sfa);
    } else {
        cfg_load(&c, arg1);
    }

    for (int i=2; i<argc-1; i+=2) {
        const char *key=argv[i]; if(key[0]=='-') key++;
        cfg_set(&c, key, argv[i+1]);
    }

    /* GPU info */
    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, 0);
    printf("GPU: %s (%.1f GB, SM %d.%d)\n\n", prop.name,
           prop.totalGlobalMem/1e9, prop.major, prop.minor);

    cfg_print(&c);

    Grid *g = grid_alloc(&c);
    printf("dx=%.4f dt=%.6f\n\n", g->dx, g->dt);

    do_init(g, &c);

    /* For gradient_pinned BC: save initial state, disable spherical damping */
    if (c.bc_type == 1) {
        grid_save_pinned(g);
        c.damp_width = 0;  c.damp_rate = 0;  /* disable spherical absorb on GPU */
        printf("Gradient BC: pinned %d slabs on each x-face (A_high=%.3f, A_low=%.3f)\n\n",
               c.gradient_margin, c.gradient_A_high, c.gradient_A_low);
    }

    /* GPU setup */
    gpu_alloc(g->N3, (c.alpha_cs != 0 || c.beta_h != 0));
    gpu_set_constants(&c, g->dx);
    cudaMemcpyToSymbol(d_BC_TYPE, &c.bc_type, sizeof(int));
    cudaMemcpyToSymbol(d_GRAD_MARGIN, &c.gradient_margin, sizeof(int));
    cudaMemcpyToSymbol(d_GRAD_A_HIGH, &c.gradient_A_high, sizeof(double));
    cudaMemcpyToSymbol(d_GRAD_A_LOW, &c.gradient_A_low, sizeof(double));
    gpu_upload(g);

    /* Initial force computation on GPU */
    if (gpu_has_intermediates) {
        compute_intermediates_kernel<<<gpu_blocks, THREADS_PER_BLOCK>>>(
            d_phi[0],d_phi[1],d_phi[2], d_theta[0],d_theta[1],d_theta[2],
            d_mismatch[0],d_mismatch[1],d_mismatch[2],
            d_harden_Q[0],d_harden_Q[1],d_harden_Q[2]);
    }
    compute_forces_kernel<<<gpu_blocks, THREADS_PER_BLOCK>>>(
        d_phi[0],d_phi[1],d_phi[2], d_theta[0],d_theta[1],d_theta[2],
        d_acc_phi[0],d_acc_phi[1],d_acc_phi[2],
        d_acc_theta[0],d_acc_theta[1],d_acc_theta[2],
        d_mismatch[0],d_mismatch[1],d_mismatch[2],
        d_harden_Q[0],d_harden_Q[1],d_harden_Q[2]);
    cudaDeviceSynchronize();

    /* SFA archive with KVMD */
    uint8_t sfa_dtype = (c.precision==0)?SFA_F16:(c.precision==1)?SFA_F32:SFA_F64;
    SFA *sfa = sfa_create(c.output, c.N, c.N, c.N, c.L, c.L, c.L, g->dt);
    sfa->flags = SFA_CODEC_COLZSTD | SFA_FLAG_STREAMING;  /* per-column parallel compression */
    if (c.output_split) {
        sfa_enable_split(sfa);
        strncpy(sfa->split_base, c.output, 507);
        char *ext = strrchr(sfa->split_base, '.');
        if (ext && !strcmp(ext, ".sfa")) *ext = '\0';
    }
    {
        char vN[32],vL[32],vT[32],vdt[32],vm[32],vmt[32],veta[32],vmu[32],vkappa[32],vkh[32],vacs[32],vbh[32];
        char vmode[32],via[32],vib[32],vkg[32],vdw[32],vdr[32],vprec[32],vdelta[64],vbcsw[32];
        snprintf(vN,32,"%d",c.N); snprintf(vL,32,"%.6f",c.L); snprintf(vT,32,"%.6f",c.T);
        snprintf(vdt,32,"%.6f",c.dt_factor);
        snprintf(vm,32,"%.6f",sqrt(c.m2)); snprintf(vmt,32,"%.6f",sqrt(c.mtheta2));
        snprintf(veta,32,"%.6f",c.eta); snprintf(vmu,32,"%.6f",c.mu);
        snprintf(vkappa,32,"%.6f",c.kappa); snprintf(vkh,32,"%.6f",c.kappa_h);
        snprintf(vacs,32,"%.6f",c.alpha_cs); snprintf(vbh,32,"%.6f",c.beta_h);
        snprintf(vmode,32,"%d",c.mode);
        snprintf(via,32,"%.6f",c.inv_alpha); snprintf(vib,32,"%.6f",c.inv_beta);
        snprintf(vkg,32,"%.6f",c.kappa_gamma);
        snprintf(vdw,32,"%.6f",c.damp_width); snprintf(vdr,32,"%.6f",c.damp_rate);
        snprintf(vprec,32,"%s",(const char*[]){"f16","f32","f64"}[c.precision]);
        snprintf(vdelta,64,"%.6f,%.6f,%.6f",c.delta[0],c.delta[1],c.delta[2]);
        snprintf(vbcsw,32,"%.6f",c.bc_switch_time);
        char vbc[32], vgah[32], vgal[32], vgm[32];
        snprintf(vbc,32,"%d",c.bc_type);
        snprintf(vgah,32,"%.6f",c.gradient_A_high);
        snprintf(vgal,32,"%.6f",c.gradient_A_low);
        snprintf(vgm,32,"%d",c.gradient_margin);
        const char *keys[]={"N","L","T","dt_factor","m","m_theta","eta","mu","kappa",
                            "kappa_h","alpha_cs","beta_h","bc_switch_time",
                            "mode","inv_alpha","inv_beta","kappa_gamma",
                            "damp_width","damp_rate","precision","delta",
                            "bc_type","gradient_A_high","gradient_A_low","gradient_margin"};
        const char *vals[]={vN,vL,vT,vdt,vm,vmt,veta,vmu,vkappa,
                            vkh,vacs,vbh,vbcsw,
                            vmode,via,vib,vkg,vdw,vdr,vprec,vdelta,
                            vbc,vgah,vgal,vgm};
        sfa_add_kvmd(sfa, 0, 0xFFFFFFFF, 0xFFFFFFFF, keys, vals, 25);
    }
    sfa_add_column(sfa,"phi_x",sfa_dtype,SFA_POSITION,0);
    sfa_add_column(sfa,"phi_y",sfa_dtype,SFA_POSITION,1);
    sfa_add_column(sfa,"phi_z",sfa_dtype,SFA_POSITION,2);
    sfa_add_column(sfa,"theta_x",sfa_dtype,SFA_ANGLE,0);
    sfa_add_column(sfa,"theta_y",sfa_dtype,SFA_ANGLE,1);
    sfa_add_column(sfa,"theta_z",sfa_dtype,SFA_ANGLE,2);
    sfa_add_column(sfa,"phi_vx",sfa_dtype,SFA_VELOCITY,0);
    sfa_add_column(sfa,"phi_vy",sfa_dtype,SFA_VELOCITY,1);
    sfa_add_column(sfa,"phi_vz",sfa_dtype,SFA_VELOCITY,2);
    sfa_add_column(sfa,"theta_vx",sfa_dtype,SFA_VELOCITY,3);
    sfa_add_column(sfa,"theta_vy",sfa_dtype,SFA_VELOCITY,4);
    sfa_add_column(sfa,"theta_vz",sfa_dtype,SFA_VELOCITY,5);
    sfa_finalize_header(sfa);
    const char *pn[]={"f16","f32","f64"};
    printf("SFA: %s (12 cols, %s, colzstd)\n\n", c.output, pn[c.precision]);

    /* Timing, step counts — use lround to avoid truncation errors */
    int n_steps=(int)lround(c.T/g->dt);
    int diag_every=(int)lround(c.diag_dt/g->dt); if(diag_every<1) diag_every=1;
    int snap_every=(int)lround(c.snap_dt/g->dt); if(snap_every<1) snap_every=1;
    int major = diag_every * 25; if(major<1) major=1;
    /* Cap so we get ~10 progress lines even for short runs */
    if(major > n_steps/10) { major = (n_steps/10/diag_every)*diag_every; if(major<diag_every) major=diag_every; }

    /* Build FieldState from device pointers */
    FieldState fstate;
    for (int a = 0; a < 3; a++) {
        fstate.phi[a] = d_phi[a]; fstate.vel_phi[a] = d_vel_phi[a]; fstate.acc_phi[a] = d_acc_phi[a];
        fstate.theta[a] = d_theta[a]; fstate.vel_theta[a] = d_vel_theta[a]; fstate.acc_theta[a] = d_acc_theta[a];
    }
    fstate.N3 = g->N3; fstate.N = g->N;
    fstate.L = g->L; fstate.dx = g->dx; fstate.dt = g->dt;

    /* Create hooks */
    SnapHookCtx snap_ctx = create_snap_hook(sfa, c.precision, snap_every, g->N3);
    start_snap_writer(&snap_ctx);
    register_hook(snap_hook, &snap_ctx);

    /* Diagnostics file */
    FILE *fp = fopen(c.diag_file, "w");
    fprintf(fp, "t\tE_phi_kin\tE_theta_kin\tE_grad\tE_mass\tE_pot\tE_tgrad\tE_tmass\t"
                "E_coupling\tE_total\tphi_max\tP_max\tP_int\ttheta_rms\n");

    cudaEvent_t t_start;
    cudaEventCreate(&t_start);
    cudaEventRecord(t_start);

    DiagHookCtx diag_ctx = create_diag_hook(fp, diag_every, major, g->N3, g->dx, n_steps, t_start, &c);
    register_hook(diag_hook, &diag_ctx);

    /* Initial diagnostics (GPU-side reduction, no full download) */
    run_gpu_diagnostics(&fstate, &diag_ctx);
    {
        double *r = diag_ctx.h_results;
        double et0 = r[0]+r[1]+r[2]+r[3]+r[4]+r[5]+r[6]+r[7];
        double trms0 = sqrt(r[11] / (3.0 * g->N3));
        diag_ctx.E0 = et0;
        printf("INIT: E_total=%.4e E_pot=%.4f phi_max=%.4f theta_rms=%.3e\n\n", et0, r[4], r[8], trms0);
        fprintf(fp, "%.2f\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\n",
                0.0, r[0], r[1], r[2], r[3], r[4], r[5], r[6], r[7], et0, r[8], r[9], r[10], trms0);
        fflush(fp);
    }

    /* t=0 snapshot via hook */
    snap_hook(0, 0.0, &fstate, &snap_ctx);

    printf("Async pipeline: snap every %d steps, diag every %d steps, %d total\n\n",
           snap_every, diag_every, n_steps);

    /* BC switch: absorb -> periodic at bc_switch_time */
    int bc_switched = 0;
    int bc_switch_step = (c.bc_switch_time > 0) ? (int)lround(c.bc_switch_time / g->dt) : 0;

    /* ===== Main loop ===== */
    for (int step = 1; step <= n_steps; step++) {
        /* Check for BC switch */
        if (bc_switch_step > 0 && step == bc_switch_step && !bc_switched) {
            double zero = 0.0;
            cudaMemcpyToSymbol(d_DAMP_RATE, &zero, sizeof(double));
            printf("\n*** BC SWITCH at t=%.1f: absorbing -> periodic ***\n\n", step * g->dt);
            bc_switched = 1;
        }

        gpu_verlet_step(g->dt);

        if (c.bc_type == 1) {
            gradient_bc_kernel<<<gpu_blocks, THREADS_PER_BLOCK>>>(
                d_phi[0], d_phi[1], d_phi[2],
                d_vel_phi[0], d_vel_phi[1], d_vel_phi[2],
                d_acc_phi[0], d_acc_phi[1], d_acc_phi[2],
                d_theta[0], d_theta[1], d_theta[2],
                d_vel_theta[0], d_vel_theta[1], d_vel_theta[2],
                d_acc_theta[0], d_acc_theta[1], d_acc_theta[2]);
        }

        /* Dispatch hooks — with burst mode snap override */
        double t = step * g->dt;
        double burst_dur = c.burst_end - c.burst_start;
        int in_burst = 0;
        if (burst_dur > 0) {
            double t_rel = t - c.burst_start;
            if (t_rel >= 0) {
                if (c.burst_every > 0) t_rel = fmod(t_rel, c.burst_every);
                in_burst = (t_rel >= 0 && t_rel < burst_dur);
            }
        }
        if (in_burst) {
            /* Burst: snap every timestep. Call snap directly, skip it in hook dispatch. */
            int saved = snap_ctx.snap_every;
            snap_ctx.snap_every = 1;
            snap_hook(1, t, &fstate, &snap_ctx);
            snap_ctx.snap_every = saved;
            /* Dispatch only non-snap hooks */
            for (int h = 0; h < n_hooks; h++)
                if (hooks[h].fn != snap_hook)
                    hooks[h].fn(step, t, &fstate, hooks[h].ctx);
        } else {
            for (int h = 0; h < n_hooks; h++)
                hooks[h].fn(step, t, &fstate, hooks[h].ctx);
        }
    }

    /* Final frame — write if the last step wasn't exactly a snap point. */
    {
        int last_snapped = (n_steps / snap_every) * snap_every;
        int gap = n_steps - last_snapped;
        if (gap > 0) {
            int saved_every = snap_ctx.snap_every;
            snap_ctx.snap_every = 1;
            snap_hook(1, n_steps * g->dt, &fstate, &snap_ctx);
            snap_ctx.snap_every = saved_every;
        }
    }

    /* Wait for writer to finish, then destroy snap hook */
    destroy_snap_hook(&snap_ctx);
    uint32_t nf = sfa->total_frames;
    sfa_close(sfa);

    /* Final diagnostics (GPU-side, no download) */
    run_gpu_diagnostics(&fstate, &diag_ctx);
    {
        double *r = diag_ctx.h_results;
        double et = r[0]+r[1]+r[2]+r[3]+r[4]+r[5]+r[6]+r[7];
        double trms = sqrt(r[11] / (3.0 * g->N3));

        cudaEvent_t t_stop;
        cudaEventCreate(&t_stop);
        cudaEventRecord(t_stop);
        cudaEventSynchronize(t_stop);
        float total_ms; cudaEventElapsedTime(&total_ms, t_start, t_stop);
        cudaEventDestroy(t_stop);

        printf("\n=== COMPLETE (GPU, async pipeline) ===\n");
        printf("E_total=%.4e (drift %.3f%%) E_pot=%.4f\n", et, 100*(et-diag_ctx.E0)/(fabs(diag_ctx.E0)+1e-30), r[4]);
        printf("phi_max=%.4f theta_rms=%.3e\n", r[8], trms);
        printf("SFA: %s (%u frames)\n", c.output, nf);
        printf("Wall: %.1fs (%.1f min) %.2fms/step\n", total_ms/1000, total_ms/60000, total_ms/n_steps);
    }

    destroy_diag_hook(&diag_ctx);
    fclose(fp);
    gpu_free();
    grid_free(g);
    cudaEventDestroy(t_start);
    return 0;
}
