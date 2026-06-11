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
#include "scp_config.h"

#include <sys/stat.h>
#include <pthread.h>
#define THREADS_PER_BLOCK 256

/* ================================================================
   Configuration — identical to scp_sim.c
   ================================================================ */

/* Config struct defined in scp_config.h */

/* cfg_defaults through f16 helpers now in scp_config.h */
#if 0  /* deleted — see scp_config.h */
static Config DELETED_cfg_defaults(void) {
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
    else if (!strcmp(key,"vec_snap_dt"))        c->vec_snap_dt = atof(val);
    else if (!strcmp(key,"vec_iframe_interval")) c->vec_iframe_interval = atoi(val);
    else if (!strcmp(key,"vec_delta_tol"))      c->vec_delta_tol = atof(val);
    else if (!strcmp(key,"vec_block_size"))     c->vec_block_size = atoi(val);
    else if (!strcmp(key,"vec_output")) {} /* deprecated, ignored */
    else if (!strcmp(key,"vec_kframe_interval")) {} /* deprecated */
    else if (!strcmp(key,"vec_n_fields")) {} /* deprecated */
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
#endif /* deleted config/f16 block */

/* ================================================================
   Hook-based async pipeline — FieldState, hooks, contexts
   ================================================================ */

typedef struct {
    double *phi[3], *vel_phi[3], *acc_phi[3];
    double *theta[3], *vel_theta[3], *acc_theta[3];
    /* v66 complex sector (device pointers; NULL when complex_mode=0) */
    int complex_mode;
    double *phi_im[3], *vel_phi_im[3], *acc_phi_im[3];
    double *theta_im[3], *vel_theta_im[3], *acc_theta_im[3];
    /* v69 gauge sector (device pointers; NULL when gauge_mode=0) */
    int gauge_mode;
    double *th[3], *Efield[3], *E_acc[3];
    long N3;
    int N;
    double L, dx, dt;
} FieldState;

#define MAX_HOOKS 8

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

/* Forward declare — defined before snap_hook implementation */
typedef struct FrameWriter_ FrameWriter;

/* Snapshot hook context — async f16 conversion + DMA + compress/write */
typedef struct {
    FrameWriter *fw;  /* shared frame writer (thread-safe) */
    int precision;        /* 0=f16, 1=f32, 2=f64 */
    int snap_every;
    int nf;               /* columns per frame: 12 real, 24 complex (v66), 30 gauged (v69) */

    /* GPU resources */
    cudaStream_t stream;
    uint16_t *d_f16_buf;  /* GPU-side f16 staging (nf * N3 * 2 bytes) */
    float *d_f32_buf;     /* GPU-side f32 staging (nf * N3 * 4 bytes) — for precision=1 */

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

    /* GPU reduction output (12 doubles real mode; 21 in complex mode, v66/v67) */
    int complex_mode;
    double *d_results;
    double *h_results;

    /* v67 frequency-domain diagnostics (allocated/used only in complex mode) */
    unsigned long long *d_sidx;  /* argmax-s voxel index (device, init ~0ULL) */
    double qdiag_R2;             /* Q_core integration radius^2 */
    long probe1_idx, probe2_idx; /* theta probe voxel indices, (x=+r_p,0,0) */
    double omega_core;           /* filled by run_gpu_diagnostics */
    double probes[4];            /* thp1u, thp1v, thp2u, thp2v */

    /* v69 gauge diagnostics (SPEC §4) */
    int gauge_mode;              /* c->complex_gauge */
    int gauged;                  /* complex_gauge && g_gauge != 0 */
    double g_gauge;
    double G_offset;             /* frozen jellium offset (set after init_gauss_project) */
    double Rint2;                /* interior Omega radius^2 (sponge excluded, bc_type=0) */
    double Rcube;                /* Q_flux cube half-width = qdiag_radius */
    double gauss_max, gauss_l2, e_em, q_flux;  /* filled by run_gpu_diagnostics */

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
    /* v66 complex sector: u=phi (Re Φ), v=phi_im (Im Φ), tu=theta, tv=theta_im */
    int complex_mode;                       /* copied from c->complex_phi at alloc */
    double *phi_im[NFIELDS],   *phi_im_vel[NFIELDS],   *phi_im_acc[NFIELDS];
    double *theta_im[NFIELDS], *theta_im_vel[NFIELDS], *theta_im_acc[NFIELDS];
    /* v69 gauged U(1) sector (SPEC §2.2) — host copies for init/projection;
     * NULL when complex_gauge=0. (No host E_acc/scratch: forces run on GPU.) */
    int    gauge_mode;            /* copied from c->complex_gauge at alloc */
    double g_gauge;               /* cached coupling */
    double G_offset;              /* frozen uniform Gauss offset, measured post-init (§4.1) */
    double *th[3];                /* link angles theta_i(x) on link (x, x+i), wrapped (-pi,pi] */
    double *Efield[3];            /* E_i(x) on the same link (noncompact) */
    int N; long N3;
    double L, dx, dt;
} Grid;

static Grid *grid_alloc(const Config *c) {
    Grid *g = (Grid*)calloc(1, sizeof(Grid));
    g->N = c->N; g->N3 = (long)c->N * c->N * c->N;
    g->L = c->L; g->dx = 2.0 * c->L / (c->N - 1);
    g->dt = c->dt_factor * g->dx;
    g->complex_mode = c->complex_phi;
    g->gauge_mode   = c->complex_gauge;
    g->g_gauge      = c->g_gauge;
    g->G_offset     = 0.0;
    /* gauge adds 6 host blocks (th 36-38, E 39-41) after the 36 complex blocks */
    long nblocks = g->gauge_mode ? 42 : (g->complex_mode ? 36 : 18);
    long total = nblocks * g->N3;
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
    /* v66: imaginary copy appended after the 18 real blocks (real offsets unchanged) */
    for (int a = 0; a < NFIELDS; a++) {
        if (g->complex_mode) {
            g->phi_im[a]       = g->mem + (18+a)*g->N3;
            g->phi_im_vel[a]   = g->mem + (21+a)*g->N3;
            g->phi_im_acc[a]   = g->mem + (24+a)*g->N3;
            g->theta_im[a]     = g->mem + (27+a)*g->N3;
            g->theta_im_vel[a] = g->mem + (30+a)*g->N3;
            g->theta_im_acc[a] = g->mem + (33+a)*g->N3;
        } else {
            g->phi_im[a] = g->phi_im_vel[a] = g->phi_im_acc[a] = NULL;
            g->theta_im[a] = g->theta_im_vel[a] = g->theta_im_acc[a] = NULL;
        }
    }
    /* v69 gauge sector: host blocks 36-41 (th 36-38, E 39-41), zero-initialized */
    for (int i = 0; i < 3; i++) {
        if (g->gauge_mode) {
            g->th[i]     = g->mem + (36 + i) * g->N3;
            g->Efield[i] = g->mem + (39 + i) * g->N3;
        } else {
            g->th[i] = g->Efield[i] = NULL;
        }
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

/* Initialization from shared header.
 * SCP_COMPLEX_FIELDS enables the v66 complex-sector init code paths in
 * scp_init.h (init=qball ansatz fill + imaginary-column SFA loading);
 * the Grid above provides the phi_im/theta_im members they reference. */
#define SCP_COMPLEX_FIELDS 1
#include "scp_init.h"

#if 0 /* deleted — now in scp_init.h */
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
#endif /* deleted init block */

/* ================================================================
   v69 host-side gauged init (mirrors scp_sim.c): Gauss residual,
   net-charge refusal (§1.2), one-time CG Gauss projection (§5.4).
   These run on the HOST grid after do_init and BEFORE gpu_upload,
   so the device starts from the exact projected state.
   ================================================================ */

/* G(x) = (1/a) sum_i [E_i(x) - E_i(x-i)] - g*rho_Q(x) at voxel idx */
static inline double gauss_residual_at(Grid *g, long idx, const long nm[3]) {
    double divE = (g->Efield[0][idx]-g->Efield[0][nm[0]]
                  +g->Efield[1][idx]-g->Efield[1][nm[1]]
                  +g->Efield[2][idx]-g->Efield[2][nm[2]]) / g->dx;
    double rho = 0;
    for (int a=0;a<NFIELDS;a++)
        rho += g->phi[a][idx]  *g->phi_im_vel[a][idx] - g->phi_im[a][idx]  *g->phi_vel[a][idx]
             + g->theta[a][idx]*g->theta_im_vel[a][idx] - g->theta_im[a][idx]*g->theta_vel[a][idx];
    return divE - g->g_gauge*rho;
}

/* §1.2 runtime net-charge check (bc_type=2 only) */
static void net_charge_check(Grid *g, const Config *c) {
    (void)c;
    const long N3=g->N3; const double dV=g->dx*g->dx*g->dx;
    double qn=0, qa=0;
    for (long idx=0;idx<N3;idx++) {
        double rho=0;
        for (int a=0;a<NFIELDS;a++)
            rho += g->phi[a][idx]  *g->phi_im_vel[a][idx] - g->phi_im[a][idx]  *g->phi_vel[a][idx]
                 + g->theta[a][idx]*g->theta_im_vel[a][idx] - g->theta_im[a][idx]*g->theta_vel[a][idx];
        qn += rho*dV; qa += fabs(rho)*dV;
    }
    if (fabs(qn) > 1e-6 * fmax(qa, 1.0)) {
        fprintf(stderr, "ERROR: bc_type=2 (periodic) requires net-neutral seed "
                "(Q_net=%.3e, Q_abs=%.3e); use +/- pairs or bc_type=0. "
                "jellium=1 is deferred (not implemented in v1).\n", qn, qa);
        exit(1);
    }
}

/* §5.4 one-time Gauss projection: CG-solve the periodic 7-point lattice
 * Poisson lap(chi) = -(G - Gbar), correct E by the forward gradient, then
 * measure and freeze G_offset. Idempotent; runs for EVERY init mode.
 * (Host-side, serial — one-time cost at init.) */
static void init_gauss_project(Grid *g, const Config *c) {
    (void)c;
    const int N=g->N, NN=N*N; const long N3=g->N3;
    const double dx=g->dx, idx2=1.0/(dx*dx);
    double *b  = (double*)malloc(N3*sizeof(double));
    double *xv = (double*)calloc(N3, sizeof(double));
    double *r  = (double*)malloc(N3*sizeof(double));
    double *p  = (double*)malloc(N3*sizeof(double));
    double *Ap = (double*)malloc(N3*sizeof(double));
    if (!b||!xv||!r||!p||!Ap) { fprintf(stderr,"FATAL: gauss projection alloc\n"); exit(1); }

    /* b = G - Gbar (zero-mean RHS; we solve (-lap) chi = b) */
    double gsum=0, gmax0=0;
    for (long idx=0;idx<N3;idx++) {
        int i=(int)(idx/NN), j=(int)((idx/N)%N), k=(int)(idx%N);
        long nm[3];
        nm[0]=(long)((i-1+N)%N)*NN+(long)j*N+k;
        nm[1]=(long)i*NN+(long)((j-1+N)%N)*N+k;
        nm[2]=(long)i*NN+(long)j*N+(k-1+N)%N;
        b[idx] = gauss_residual_at(g, idx, nm);
        gsum += b[idx];
    }
    double gbar = gsum/(double)N3;
    for (long idx=0;idx<N3;idx++) {
        b[idx] -= gbar;
        double d=fabs(b[idx]); if (d>gmax0) gmax0=d;
    }

    const double tol = 1e-13 * fmax(1.0, gmax0);
    memcpy(r, b, N3*sizeof(double));
    memcpy(p, b, N3*sizeof(double));
    double rz=0;
    for (long idx=0;idx<N3;idx++) rz += r[idx]*r[idx];

    int it=0;
    const int itmax=20000;
    double rmax=gmax0;
    while (rmax > tol && it < itmax) {
        double pAp=0;
        for (long idx=0;idx<N3;idx++) {
            int i=(int)(idx/NN), j=(int)((idx/N)%N), k=(int)(idx%N);
            long nip=(long)((i+1)%N)*NN+(long)j*N+k, nim=(long)((i-1+N)%N)*NN+(long)j*N+k;
            long njp=(long)i*NN+(long)((j+1)%N)*N+k, njm=(long)i*NN+(long)((j-1+N)%N)*N+k;
            long nkp=(long)i*NN+(long)j*N+(k+1)%N,   nkm=(long)i*NN+(long)j*N+(k-1+N)%N;
            double lap = (p[nip]+p[nim]+p[njp]+p[njm]+p[nkp]+p[nkm]-6.0*p[idx])*idx2;
            Ap[idx] = -lap;                       /* A = -lap (PSD on zero-mean) */
            pAp += p[idx]*Ap[idx];
        }
        double alpha = rz/pAp;
        double rz2=0, xsum=0, rm=0;
        for (long idx=0;idx<N3;idx++) {
            xv[idx] += alpha*p[idx];
            r[idx]  -= alpha*Ap[idx];
            rz2 += r[idx]*r[idx];
            xsum += xv[idx];
            double d=fabs(r[idx]); if (d>rm) rm=d;
        }
        /* keep chi zero-mean per iteration (torus null space) */
        double xbar = xsum/(double)N3;
        double beta = rz2/rz;
        for (long idx=0;idx<N3;idx++) {
            xv[idx] -= xbar;
            p[idx] = r[idx] + beta*p[idx];
        }
        rz = rz2; rmax = rm; it++;
    }
    if (rmax > tol) {
        fprintf(stderr, "FATAL: Gauss projection CG failed to converge "
                "(%d iters, max|res|=%.3e > tol=%.3e)\n", it, rmax, tol);
        exit(1);
    }

    /* E correction by the forward gradient: div-correction = lap(chi) exactly */
    for (long idx=0;idx<N3;idx++) {
        int i=(int)(idx/NN), j=(int)((idx/N)%N), k=(int)(idx%N);
        long np[3];
        np[0]=(long)((i+1)%N)*NN+(long)j*N+k;
        np[1]=(long)i*NN+(long)((j+1)%N)*N+k;
        np[2]=(long)i*NN+(long)j*N+(k+1)%N;
        for (int d=0;d<3;d++)
            g->Efield[d][idx] += (xv[np[d]] - xv[idx])/dx;
    }

    /* measure + freeze G_offset; report post-projection residual */
    double osum=0;
    for (long idx=0;idx<N3;idx++) {
        int i=(int)(idx/NN), j=(int)((idx/N)%N), k=(int)(idx%N);
        long nm[3];
        nm[0]=(long)((i-1+N)%N)*NN+(long)j*N+k;
        nm[1]=(long)i*NN+(long)((j-1+N)%N)*N+k;
        nm[2]=(long)i*NN+(long)j*N+(k-1+N)%N;
        osum += gauss_residual_at(g, idx, nm);
    }
    g->G_offset = osum/(double)N3;
    double gres=0;
    for (long idx=0;idx<N3;idx++) {
        int i=(int)(idx/NN), j=(int)((idx/N)%N), k=(int)(idx%N);
        long nm[3];
        nm[0]=(long)((i-1+N)%N)*NN+(long)j*N+k;
        nm[1]=(long)i*NN+(long)((j-1+N)%N)*N+k;
        nm[2]=(long)i*NN+(long)j*N+(k-1+N)%N;
        double d=fabs(gauss_residual_at(g, idx, nm) - g->G_offset);
        if (d>gres) gres=d;
    }
    printf("Gauss projection: %d CG iters, gauss_max(0)=%.2e (offset %.3e)\n",
           it, gres, g->G_offset);
    free(b); free(xv); free(r); free(p); free(Ap);
}

/* ================================================================
   GPU constant memory
   ================================================================ */

__constant__ double d_MU, d_KAPPA, d_MASS2, d_MTHETA2, d_ETA;
__constant__ double d_KAPPA_H, d_ALPHA_CS, d_BETA_H;
__constant__ double d_inv_alpha, d_inv_beta, d_kappa_gamma;
__constant__ double d_THETA_SAT, d_GAMMA_CONV;
__constant__ double d_SIGMA_GRAD, d_CHI_CHIRAL, d_SIGMA_CUBIC, d_SIGMA_FREQ;
__constant__ double d_SIGMA_CROSS, d_LAMBDA_SELF;
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
    const double *tvel0, const double *tvel1, const double *tvel2,
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

    /* --- Theta saturation → phi conversion --- */
    if (d_THETA_SAT > 0 && d_GAMMA_CONV > 0) {
        double th2 = th[0][idx]*th[0][idx] + th[1][idx]*th[1][idx] + th[2][idx]*th[2][idx];
        double excess = th2 - d_THETA_SAT * d_THETA_SAT;
        if (excess > 0) {
            /* curl(theta) at this voxel */
            double ct_sat[3];
            ct_sat[0] = (th[2][n_jp]-th[2][n_jm]-th[1][n_kp]+th[1][n_km])*d_idx1;
            ct_sat[1] = (th[0][n_kp]-th[0][n_km]-th[2][n_ip]+th[2][n_im])*d_idx1;
            ct_sat[2] = (th[1][n_ip]-th[1][n_im]-th[0][n_jp]+th[0][n_jm])*d_idx1;
            for (int a = 0; a < 3; a++) {
                /* Drain theta: -γ·excess·θ_a */
                acc_t[a] -= d_GAMMA_CONV * excess * th[a][idx];
                /* Source phi from curl(theta): +γ·excess·curl(θ)_a */
                acc_p[a] += d_GAMMA_CONV * excess * ct_sat[a];
            }
        }
    }

    /* --- Test A: Gradient-based theta self-interaction --- */
    if (d_SIGMA_GRAD > 0) {
        /* |∇θ|² = sum over components of |grad(theta_a)|² */
        double grad_th_sq = 0;
        double ct_g[3];  /* curl(theta) for conversion */
        ct_g[0] = (th[2][n_jp]-th[2][n_jm]-th[1][n_kp]+th[1][n_km])*d_idx1;
        ct_g[1] = (th[0][n_kp]-th[0][n_km]-th[2][n_ip]+th[2][n_im])*d_idx1;
        ct_g[2] = (th[1][n_ip]-th[1][n_im]-th[0][n_jp]+th[0][n_jm])*d_idx1;
        for (int a = 0; a < 3; a++) {
            double dtdx = (th[a][n_ip]-th[a][n_im])*d_idx1;
            double dtdy = (th[a][n_jp]-th[a][n_jm])*d_idx1;
            double dtdz = (th[a][n_kp]-th[a][n_km])*d_idx1;
            grad_th_sq += dtdx*dtdx + dtdy*dtdy + dtdz*dtdz;
        }
        /* Chirality bias: h = θ·curl(θ) */
        double th2_g = th[0][idx]*th[0][idx]+th[1][idx]*th[1][idx]+th[2][idx]*th[2][idx];
        double h = th[0][idx]*ct_g[0]+th[1][idx]*ct_g[1]+th[2][idx]*ct_g[2];
        double sigma_eff = d_SIGMA_GRAD * (1.0 + d_CHI_CHIRAL * h / (th2_g + 0.01));
        if (sigma_eff < 0) sigma_eff = 0;
        for (int a = 0; a < 3; a++) {
            acc_t[a] -= sigma_eff * grad_th_sq * th[a][idx];
            acc_p[a] += sigma_eff * grad_th_sq * ct_g[a];
        }
    }

    /* --- Test B: Cubic theta self-interaction --- */
    if (d_SIGMA_CUBIC > 0) {
        double th2_c = th[0][idx]*th[0][idx]+th[1][idx]*th[1][idx]+th[2][idx]*th[2][idx];
        double ct_c[3];
        ct_c[0] = (th[2][n_jp]-th[2][n_jm]-th[1][n_kp]+th[1][n_km])*d_idx1;
        ct_c[1] = (th[0][n_kp]-th[0][n_km]-th[2][n_ip]+th[2][n_im])*d_idx1;
        ct_c[2] = (th[1][n_ip]-th[1][n_im]-th[0][n_jp]+th[0][n_jm])*d_idx1;
        for (int a = 0; a < 3; a++) {
            acc_t[a] -= d_SIGMA_CUBIC * th2_c * th[a][idx];
            acc_p[a] += d_SIGMA_CUBIC * th2_c * ct_c[a];
        }
    }

    /* --- Test C: Frequency-mismatch conversion --- */
    if (d_SIGMA_FREQ > 0) {
        /* Local omega at center: |theta_vel| / (|theta| + eps) */
        double tv_sq = tvel0[idx]*tvel0[idx]+tvel1[idx]*tvel1[idx]+tvel2[idx]*tvel2[idx];
        double th2_f = th[0][idx]*th[0][idx]+th[1][idx]*th[1][idx]+th[2][idx]*th[2][idx];
        double omega_c = sqrt(tv_sq) / (sqrt(th2_f) + 0.01);

        /* Average omega of 6 neighbors */
        double dw_sum = 0;
        int nn = 0;
        long nbrs[6] = {n_ip, n_im, n_jp, n_jm, n_kp, n_km};
        for (int nb = 0; nb < 6; nb++) {
            long nidx = nbrs[nb];
            double tv_n = tvel0[nidx]*tvel0[nidx]+tvel1[nidx]*tvel1[nidx]+tvel2[nidx]*tvel2[nidx];
            double th_n = th[0][nidx]*th[0][nidx]+th[1][nidx]*th[1][nidx]+th[2][nidx]*th[2][nidx];
            double omega_n = sqrt(tv_n) / (sqrt(th_n) + 0.01);
            dw_sum += fabs(omega_c - omega_n);
            nn++;
        }
        double dw_avg = dw_sum / nn;

        double ct_f[3];
        ct_f[0] = (th[2][n_jp]-th[2][n_jm]-th[1][n_kp]+th[1][n_km])*d_idx1;
        ct_f[1] = (th[0][n_kp]-th[0][n_km]-th[2][n_ip]+th[2][n_im])*d_idx1;
        ct_f[2] = (th[1][n_ip]-th[1][n_im]-th[0][n_jp]+th[0][n_jm])*d_idx1;
        for (int a = 0; a < 3; a++) {
            acc_t[a] -= d_SIGMA_FREQ * dw_avg * th[a][idx];
            acc_p[a] += d_SIGMA_FREQ * dw_avg * ct_f[a];
        }
    }

    /* --- Lagrangian cross-potential: -(sigma/2)|θ|²|φ|² --- */
    if (d_SIGMA_CROSS > 0) {
        double th2_x = th[0][idx]*th[0][idx]+th[1][idx]*th[1][idx]+th[2][idx]*th[2][idx];
        double ph2_x = phi[0][idx]*phi[0][idx]+phi[1][idx]*phi[1][idx]+phi[2][idx]*phi[2][idx];
        for (int a = 0; a < 3; a++) {
            acc_p[a] -= d_SIGMA_CROSS * th2_x * phi[a][idx];
            acc_t[a] -= d_SIGMA_CROSS * ph2_x * th[a][idx];
        }
    }

    /* --- Lagrangian theta self-potential: -(lambda/4)|θ|⁴ --- */
    if (d_LAMBDA_SELF > 0) {
        double th2_l = th[0][idx]*th[0][idx]+th[1][idx]*th[1][idx]+th[2][idx]*th[2][idx];
        for (int a = 0; a < 3; a++) {
            acc_t[a] -= d_LAMBDA_SELF * th2_l * th[a][idx];
        }
    }

    acc_phi0[idx]=acc_p[0]; acc_phi1[idx]=acc_p[1]; acc_phi2[idx]=acc_p[2];
    acc_theta0[idx]=acc_t[0]; acc_theta1[idx]=acc_t[1]; acc_theta2[idx]=acc_t[2];
}

/* v66 complex forces (SPEC §3, mirrors CPU compute_forces_complex):
 * minimal 12-field loop, THEORY.md §1 EOM verbatim. Curl pairs (u,tu)
 * and (v,tv) ONLY — no u<->tv or v<->tu coupling. cfg_validate guarantees
 * alpha_cs=beta_h=kappa_h=0, mode=0, sigma_*=theta_sat=0 under complex,
 * so no intermediates pass, no chiral block, no mode/conversion terms. */
__global__ void compute_forces_complex_kernel(
    const double *u0, const double *u1, const double *u2,
    const double *v0, const double *v1, const double *v2,
    const double *tu0, const double *tu1, const double *tu2,
    const double *tv0, const double *tv1, const double *tv2,
    double *acc_u0, double *acc_u1, double *acc_u2,
    double *acc_v0, double *acc_v1, double *acc_v2,
    double *acc_tu0, double *acc_tu1, double *acc_tu2,
    double *acc_tv0, double *acc_tv1, double *acc_tv2)
{
    long idx = (long)blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= d_N3) return;

    int N = d_N, NN = d_NN;
    int i = (int)(idx / NN), j = (int)((idx / N) % N), k = (int)(idx % N);
    int ip=(i+1)%N, im=(i-1+N)%N, jp=(j+1)%N, jm=(j-1+N)%N, kp=(k+1)%N, km=(k-1+N)%N;
    long n_ip=(long)ip*NN+j*N+k, n_im=(long)im*NN+j*N+k;
    long n_jp=(long)i*NN+jp*N+k, n_jm=(long)i*NN+jm*N+k;
    long n_kp=(long)i*NN+j*N+kp, n_km=(long)i*NN+j*N+km;

    const double *u[3]  = {u0, u1, u2};
    const double *v[3]  = {v0, v1, v2};
    const double *tu[3] = {tu0, tu1, tu2};
    const double *tv[3] = {tv0, tv1, tv2};

    /* amplitudes */
    double uu0=u0[idx], uu1=u1[idx], uu2=u2[idx];
    double vv0=v0[idx], vv1=v1[idx], vv2=v2[idx];
    double s2_0 = uu0*uu0 + vv0*vv0;          /* |Phi_0|^2 */
    double s2_1 = uu1*uu1 + vv1*vv1;
    double s2_2 = uu2*uu2 + vv2*vv2;
    double sv   = s2_0 * s2_1 * s2_2;         /* s = prod |Phi_a|^2 */
    double den  = 1.0 + d_KAPPA*sv;
    double Vp   = 0.5*d_MU / (den*den);       /* Vt'(s) = (mu/2)/(1+kappa s)^2 */
    /* product over b != a — explicit pair products, NEVER s/s2_a */
    double prod_rest[3] = { s2_1*s2_2, s2_0*s2_2, s2_0*s2_1 };

    double acc_u[3], acc_v[3], acc_tu[3], acc_tv[3];
    for (int a = 0; a < 3; a++) {
        /* u sector: pairs with tu via curl */
        double lap_u = (u[a][n_ip]+u[a][n_im]+u[a][n_jp]+u[a][n_jm]
                       +u[a][n_kp]+u[a][n_km]-6.0*u[a][idx]) * d_idx2;
        double ctu;
        if (a==0)      ctu = (tu[2][n_jp]-tu[2][n_jm]-tu[1][n_kp]+tu[1][n_km])*d_idx1;
        else if (a==1) ctu = (tu[0][n_kp]-tu[0][n_km]-tu[2][n_ip]+tu[2][n_im])*d_idx1;
        else           ctu = (tu[1][n_ip]-tu[1][n_im]-tu[0][n_jp]+tu[0][n_jm])*d_idx1;
        acc_u[a] = lap_u - d_MASS2*u[a][idx]
                 - 2.0*Vp*u[a][idx]*prod_rest[a] + d_ETA*ctu;

        /* v sector: pairs with tv via curl */
        double lap_v = (v[a][n_ip]+v[a][n_im]+v[a][n_jp]+v[a][n_jm]
                       +v[a][n_kp]+v[a][n_km]-6.0*v[a][idx]) * d_idx2;
        double ctv;
        if (a==0)      ctv = (tv[2][n_jp]-tv[2][n_jm]-tv[1][n_kp]+tv[1][n_km])*d_idx1;
        else if (a==1) ctv = (tv[0][n_kp]-tv[0][n_km]-tv[2][n_ip]+tv[2][n_im])*d_idx1;
        else           ctv = (tv[1][n_ip]-tv[1][n_im]-tv[0][n_jp]+tv[0][n_jm])*d_idx1;
        acc_v[a] = lap_v - d_MASS2*v[a][idx]
                 - 2.0*Vp*v[a][idx]*prod_rest[a] + d_ETA*ctv;

        /* tu sector */
        double lap_tu = (tu[a][n_ip]+tu[a][n_im]+tu[a][n_jp]+tu[a][n_jm]
                        +tu[a][n_kp]+tu[a][n_km]-6.0*tu[a][idx]) * d_idx2;
        double cu;
        if (a==0)      cu = (u[2][n_jp]-u[2][n_jm]-u[1][n_kp]+u[1][n_km])*d_idx1;
        else if (a==1) cu = (u[0][n_kp]-u[0][n_km]-u[2][n_ip]+u[2][n_im])*d_idx1;
        else           cu = (u[1][n_ip]-u[1][n_im]-u[0][n_jp]+u[0][n_jm])*d_idx1;
        acc_tu[a] = lap_tu - d_MTHETA2*tu[a][idx] + d_ETA*cu;

        /* tv sector */
        double lap_tv = (tv[a][n_ip]+tv[a][n_im]+tv[a][n_jp]+tv[a][n_jm]
                        +tv[a][n_kp]+tv[a][n_km]-6.0*tv[a][idx]) * d_idx2;
        double cv;
        if (a==0)      cv = (v[2][n_jp]-v[2][n_jm]-v[1][n_kp]+v[1][n_km])*d_idx1;
        else if (a==1) cv = (v[0][n_kp]-v[0][n_km]-v[2][n_ip]+v[2][n_im])*d_idx1;
        else           cv = (v[1][n_ip]-v[1][n_im]-v[0][n_jp]+v[0][n_jm])*d_idx1;
        acc_tv[a] = lap_tv - d_MTHETA2*tv[a][idx] + d_ETA*cv;
    }

    acc_u0[idx]=acc_u[0]; acc_u1[idx]=acc_u[1]; acc_u2[idx]=acc_u[2];
    acc_v0[idx]=acc_v[0]; acc_v1[idx]=acc_v[1]; acc_v2[idx]=acc_v[2];
    acc_tu0[idx]=acc_tu[0]; acc_tu1[idx]=acc_tu[1]; acc_tu2[idx]=acc_tu[2];
    acc_tv0[idx]=acc_tv[0]; acc_tv1[idx]=acc_tv[1]; acc_tv2[idx]=acc_tv[2];
}

/* ================================================================
   v69 gauged complex forces (SPEC §3, mirrors CPU
   compute_forces_complex_gauge): 12 matter fields + compact U(1)
   links. ALL spatial matter differences are link-covariant
   (forward+backward transported neighbors, §3.1-3.3); the E-kick is
   plaquette staples + the unique Gauss-conserving lattice current
   (§3.4). Never launched when g_gauge==0 (§1.3 dispatch).
   The CPU's link cos/sin (pass A) and plaquette-sine (pass B)
   scratch arrays are recomputed INLINE here — only th/E/E_acc
   (9 blocks) live on the device.
   ================================================================ */

/* sin(theta_P) of the plane-(a,b) plaquette anchored at lattice site
 * (i,j,k): ang = th_a(x) + th_b(x+a) - th_a(x+b) - th_b(x) */
__device__ inline double gauge_plaq_sin(const double *tha, const double *thb,
    int a, int b, int i, int j, int k, int N, int NN)
{
    long s = (long)i*NN + (long)j*N + k;
    int ia=i, ja=j, ka=k;
    if (a==0) ia=(i+1)%N; else if (a==1) ja=(j+1)%N; else ka=(k+1)%N;
    long sa = (long)ia*NN + (long)ja*N + ka;
    int ib=i, jb=j, kb=k;
    if (b==0) ib=(i+1)%N; else if (b==1) jb=(j+1)%N; else kb=(k+1)%N;
    long sb = (long)ib*NN + (long)jb*N + kb;
    return sin(tha[s] + thb[sa] - tha[sb] - thb[s]);
}

__global__ void compute_forces_complex_gauge_kernel(
    const double *u0, const double *u1, const double *u2,
    const double *v0, const double *v1, const double *v2,
    const double *tu0, const double *tu1, const double *tu2,
    const double *tv0, const double *tv1, const double *tv2,
    const double *th0, const double *th1, const double *th2,
    double *acc_u0, double *acc_u1, double *acc_u2,
    double *acc_v0, double *acc_v1, double *acc_v2,
    double *acc_tu0, double *acc_tu1, double *acc_tu2,
    double *acc_tv0, double *acc_tv1, double *acc_tv2,
    double *Eacc0, double *Eacc1, double *Eacc2,
    double GG, double STP, double inva)
{
    long idx = (long)blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= d_N3) return;

    int N = d_N, NN = d_NN;
    int i = (int)(idx / NN), j = (int)((idx / N) % N), k = (int)(idx % N);
    long np[3], nm[3];
    np[0] = (long)((i+1)%N)*NN + (long)j*N + k;
    nm[0] = (long)((i-1+N)%N)*NN + (long)j*N + k;
    np[1] = (long)i*NN + (long)((j+1)%N)*N + k;
    nm[1] = (long)i*NN + (long)((j-1+N)%N)*N + k;
    np[2] = (long)i*NN + (long)j*N + (k+1)%N;
    nm[2] = (long)i*NN + (long)j*N + (k-1+N)%N;

    const double *u[3]  = {u0, u1, u2};
    const double *v[3]  = {v0, v1, v2};
    const double *tu[3] = {tu0, tu1, tu2};
    const double *tv[3] = {tv0, tv1, tv2};
    const double *th[3] = {th0, th1, th2};
    double *Eacc[3] = {Eacc0, Eacc1, Eacc2};

    /* fields at x: f=0..2 Phi_a (u,v), f=3..5 Theta_a (tu,tv) */
    double fu[6], fv[6];
    for (int a = 0; a < 3; a++) {
        fu[a]   = u[a][idx];   fv[a]   = v[a][idx];
        fu[3+a] = tu[a][idx];  fv[3+a] = tv[a][idx];
    }

    /* transported neighbors (SPEC §3.1) per direction d, field f */
    double TRp[3][6], TIp[3][6], TRm[3][6], TIm[3][6];
    /* raw +d theta neighbors needed for the symmetrized eta current */
    double tup[3][3], tvp[3][3];
    /* forward-link cos/sin at x (CPU link_c/link_s scratch, inline) */
    double cPd[3], sPd[3];
    for (int d = 0; d < 3; d++) {
        double cP = cos(th[d][idx]),   sP = sin(th[d][idx]);
        double cM = cos(th[d][nm[d]]), sM = sin(th[d][nm[d]]);
        cPd[d] = cP;  sPd[d] = sP;
        for (int a = 0; a < 3; a++) {
            double upn = u[a][np[d]],   vpn = v[a][np[d]];
            double umn = u[a][nm[d]],   vmn = v[a][nm[d]];
            double tpn = tu[a][np[d]],  wpn = tv[a][np[d]];
            double tmn = tu[a][nm[d]],  wmn = tv[a][nm[d]];
            TRp[d][a]   = cP*upn - sP*vpn;   TIp[d][a]   = cP*vpn + sP*upn;
            TRm[d][a]   = cM*umn + sM*vmn;   TIm[d][a]   = cM*vmn - sM*umn;
            TRp[d][3+a] = cP*tpn - sP*wpn;   TIp[d][3+a] = cP*wpn + sP*tpn;
            TRm[d][3+a] = cM*tmn + sM*wmn;   TIm[d][3+a] = cM*wmn - sM*tmn;
            tup[d][a] = tpn;  tvp[d][a] = wpn;
        }
    }

    /* potential pieces (identical to compute_forces_complex) */
    double s2_0 = fu[0]*fu[0] + fv[0]*fv[0];
    double s2_1 = fu[1]*fu[1] + fv[1]*fv[1];
    double s2_2 = fu[2]*fu[2] + fv[2]*fv[2];
    double s    = s2_0 * s2_1 * s2_2;
    double den  = 1.0 + d_KAPPA*s;
    double Vp   = 0.5*d_MU / (den*den);
    double prod_rest[3] = { s2_1*s2_2, s2_0*s2_2, s2_0*s2_1 };

    /* covariant central differences (SPEC §3.3): Dc[d][f] */
    double DcU[3][6], DcV[3][6];
    for (int d = 0; d < 3; d++)
        for (int f = 0; f < 6; f++) {
            DcU[d][f] = (TRp[d][f] - TRm[d][f]) * d_idx1;
            DcV[d][f] = (TIp[d][f] - TIm[d][f]) * d_idx1;
        }

    /* covariant curls: (DxF)_0 = Dc_1 F_2 - Dc_2 F_1, cyclic */
    const int ci1[3] = {1,2,0}, ci2[3] = {2,0,1};
    for (int a = 0; a < 3; a++) {
        int d1 = ci1[a], d2 = ci2[a];
        double reDxT = DcU[d1][3+ci2[a]] - DcU[d2][3+ci1[a]];
        double imDxT = DcV[d1][3+ci2[a]] - DcV[d2][3+ci1[a]];
        double reDxP = DcU[d1][ci2[a]]   - DcU[d2][ci1[a]];
        double imDxP = DcV[d1][ci2[a]]   - DcV[d2][ci1[a]];

        /* covariant Laplacians (SPEC §3.2) */
        double lapU_P = (TRp[0][a]+TRm[0][a] + TRp[1][a]+TRm[1][a]
                       + TRp[2][a]+TRm[2][a] - 6.0*fu[a]) * d_idx2;
        double lapV_P = (TIp[0][a]+TIm[0][a] + TIp[1][a]+TIm[1][a]
                       + TIp[2][a]+TIm[2][a] - 6.0*fv[a]) * d_idx2;
        double lapU_T = (TRp[0][3+a]+TRm[0][3+a] + TRp[1][3+a]+TRm[1][3+a]
                       + TRp[2][3+a]+TRm[2][3+a] - 6.0*fu[3+a]) * d_idx2;
        double lapV_T = (TIp[0][3+a]+TIm[0][3+a] + TIp[1][3+a]+TIm[1][3+a]
                       + TIp[2][3+a]+TIm[2][3+a] - 6.0*fv[3+a]) * d_idx2;

        double au = lapU_P - d_MASS2*fu[a]
                  - 2.0*Vp*fu[a]*prod_rest[a] + d_ETA*reDxT;
        double av = lapV_P - d_MASS2*fv[a]
                  - 2.0*Vp*fv[a]*prod_rest[a] + d_ETA*imDxT;
        double atu = lapU_T - d_MTHETA2*fu[3+a] + d_ETA*reDxP;
        double atv = lapV_T - d_MTHETA2*fv[3+a] + d_ETA*imDxP;
        if (a == 0)      { acc_u0[idx]=au; acc_v0[idx]=av; acc_tu0[idx]=atu; acc_tv0[idx]=atv; }
        else if (a == 1) { acc_u1[idx]=au; acc_v1[idx]=av; acc_tu1[idx]=atu; acc_tv1[idx]=atv; }
        else             { acc_u2[idx]=au; acc_v2[idx]=av; acc_tu2[idx]=atu; acc_tv2[idx]=atv; }
    }

    /* staple sums (SPEC §3.4 table): plaquette sines recomputed inline,
     * planes p0=(0,1), p1=(1,2), p2=(2,0), at x and x-1/x-j/x-k */
    int im = (i-1+N)%N, jm = (j-1+N)%N, km = (k-1+N)%N;
    double P0_x  = gauge_plaq_sin(th[0], th[1], 0, 1, i,  j,  k,  N, NN);
    double P0_im = gauge_plaq_sin(th[0], th[1], 0, 1, im, j,  k,  N, NN);
    double P0_jm = gauge_plaq_sin(th[0], th[1], 0, 1, i,  jm, k,  N, NN);
    double P1_x  = gauge_plaq_sin(th[1], th[2], 1, 2, i,  j,  k,  N, NN);
    double P1_jm = gauge_plaq_sin(th[1], th[2], 1, 2, i,  jm, k,  N, NN);
    double P1_km = gauge_plaq_sin(th[1], th[2], 1, 2, i,  j,  km, N, NN);
    double P2_x  = gauge_plaq_sin(th[2], th[0], 2, 0, i,  j,  k,  N, NN);
    double P2_im = gauge_plaq_sin(th[2], th[0], 2, 0, im, j,  k,  N, NN);
    double P2_km = gauge_plaq_sin(th[2], th[0], 2, 0, i,  j,  km, N, NN);
    double st[3];
    st[0] =  (P0_x - P0_jm) - (P2_x - P2_km);
    st[1] = -(P0_x - P0_im) + (P1_x - P1_km);
    st[2] =  (P2_x - P2_im) - (P1_x - P1_jm);

    /* lattice current J_i^lat (SPEC §3.4) + E-kick */
    for (int d = 0; d < 3; d++) {
        double cP = cPd[d], sP = sPd[d];
        double Jg = 0.0;
        for (int a = 0; a < 3; a++) {
            Jg += fu[a]  *TIp[d][a]   - fv[a]  *TRp[d][a];     /* Phi link current */
            Jg += fu[3+a]*TIp[d][3+a] - fv[3+a]*TRp[d][3+a];   /* Theta ditto */
        }
        /* eta seagull, symmetrized over link ends:
         * T(j,k) = Im[Th_j(x)~ U Phi_k(x+d)] + Im[Th_j(x+d)~ U+ Phi_k(x)],
         * J_eta_d = -(eta/2) * (T(d1,d2) - T(d2,d1)), (d,d1,d2) cyclic */
        int d1 = ci1[d], d2 = ci2[d];
        double WpR1 = cP*fu[d2] + sP*fv[d2], WpI1 = cP*fv[d2] - sP*fu[d2];
        double WpR2 = cP*fu[d1] + sP*fv[d1], WpI2 = cP*fv[d1] - sP*fu[d1];
        double T12 = (fu[3+d1]*TIp[d][d2] - fv[3+d1]*TRp[d][d2])
                   + (tup[d][d1]*WpI1     - tvp[d][d1]*WpR1);
        double T21 = (fu[3+d2]*TIp[d][d1] - fv[3+d2]*TRp[d][d1])
                   + (tup[d][d2]*WpI2     - tvp[d][d2]*WpR2);
        double Jlat = inva*Jg - 0.5*d_ETA*(T12 - T21);
        Eacc[d][idx] = STP*st[d] + GG*Jlat;
    }
}

/* v69 E half-kick: E_i += hdt * K_i (3 link-field arrays) */
__global__ void gauge_halfkick_kernel(
    double *E0, double *E1, double *E2,
    const double *K0, const double *K1, const double *K2, double hdt)
{
    long idx = (long)blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= d_N3) return;
    E0[idx]+=hdt*K0[idx]; E1[idx]+=hdt*K1[idx]; E2[idx]+=hdt*K2[idx];
}

/* v69 link drift: th_i += gad*E_i (gad = -g*a*dt), wrapped to (-pi,pi] */
__global__ void gauge_link_drift_kernel(
    double *t0, double *t1, double *t2,
    const double *E0, const double *E1, const double *E2, double gad)
{
    long idx = (long)blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= d_N3) return;
    double w;
    w = t0[idx] + gad*E0[idx];  t0[idx] = w - 2.0*PI*rint(w/(2.0*PI));
    w = t1[idx] + gad*E1[idx];  t1[idx] = w - 2.0*PI*rint(w/(2.0*PI));
    w = t2[idx] + gad*E2[idx];  t2[idx] = w - 2.0*PI*rint(w/(2.0*PI));
}

/* v69 §3.6: damp E (velocity-like) in the sponge; th links are
 * positions and are NOT damped */
__global__ void absorbing_boundary_gauge_kernel(
    double *E0, double *E1, double *E2)
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
        E0[idx]*=damp; E1[idx]*=damp; E2[idx]*=damp;
    }
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

/* v66: bc_type=1 imaginary sector (complex_phi=1). The analytical gradient
 * background is purely real, so the x-slabs pin the imaginary sector to
 * ZERO (the CPU kernel pins to the saved initial imaginary state, which is
 * zero for every init mode legal under bc_type=1 + complex). y edges:
 * linear extrapolation, same stencil as the real kernel. Launched right
 * after gradient_bc_kernel with the *_im device pointers. */
__global__ void gradient_bc_im_kernel(
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

    /* x-direction: pin boundary slabs (imaginary background = 0) */
    if (i < M || i >= N - M) {
        p0[idx] = 0; p1[idx] = 0; p2[idx] = 0;
        vp0[idx] = 0; vp1[idx] = 0; vp2[idx] = 0;
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

/* ================================================================
   v66 complex diagnostics reduction (mirrors CPU compute_energy_complex,
   theta_rms_complex, P_integrated_complex and pass 1 of compute_charges).
   19 values in one pass:
   [0] E_phi_kin, [1] E_theta_kin, [2] E_grad, [3] E_mass,
   [4] E_pot = int Vt(s) dV, [5] E_tgrad, [6] E_tmass, [7] E_coupling,
   [8] phi_max = max_a |Phi_a| (max), [9] P_max = max sqrt(s) (max),
   [10] P_int = int sqrt(s) dV, [11] theta_rms_sum (tu^2+tv^2, /(6 N^3) on host),
   [12] Q_phi_sum  = sum_a (u*vdot - v*udot)   (raw, *dV on host),
   [13] Q_theta_sum= sum_a (tu*tvdot - tv*tudot) (raw, *dV on host),
   [14] s_max (max), [15] M0_sum = sum rho2 (raw),
   [16] cx_sum = x*rho2, [17] cy_sum, [18] cz_sum,
   [19] Q_core_sum = (qp+qt) masked to x^2+y^2+z^2 <= qdiag_radius^2
        (raw, *dV on host; v67).
   [20] (filled by reduce_rcore_kernel, second pass): rsum for r_core.
   ================================================================ */
#define CDIAG_NVALS 20
#define CDIAG_TOTAL 21
/* v69: gauss reduction appends 6 slots after the complex 21:
 * [21] gauss_max (max), [22] gauss_l2 sum, [23] nOmega,
 * [24] E_em, [25] divE cube sum, [26] nC */
#define GAUSS_BASE  21
#define GDIAG_TOTAL 27

__global__ void reduce_diagnostics_complex_kernel(
    const double *u0, const double *u1, const double *u2,
    const double *v0, const double *v1, const double *v2,
    const double *tu0, const double *tu1, const double *tu2,
    const double *tv0, const double *tv1, const double *tv2,
    const double *vu0, const double *vu1, const double *vu2,
    const double *vv0, const double *vv1, const double *vv2,
    const double *vtu0, const double *vtu1, const double *vtu2,
    const double *vtv0, const double *vtv1, const double *vtv2,
    double dV, double idx1,
    double mass2, double mtheta2, double eta, double mu, double kappa,
    double qr2, double *d_results)
{
    __shared__ double sdata[CDIAG_NVALS * 256];  /* THREADS_PER_BLOCK * CDIAG_NVALS */

    long idx = (long)blockIdx.x * blockDim.x + threadIdx.x;
    int tid = threadIdx.x;

    double local[CDIAG_NVALS];
    for (int vv = 0; vv < CDIAG_NVALS; vv++) local[vv] = 0.0;

    if (idx < d_N3) {
        int N = d_N, NN = d_NN;
        int i = (int)(idx / NN), j = (int)((idx / N) % N), k = (int)(idx % N);
        int ip=(i+1)%N, im=(i-1+N)%N, jp=(j+1)%N, jm=(j-1+N)%N, kp=(k+1)%N, km=(k-1+N)%N;
        long n_ip=(long)ip*NN+j*N+k, n_im=(long)im*NN+j*N+k;
        long n_jp=(long)i*NN+jp*N+k, n_jm=(long)i*NN+jm*N+k;
        long n_kp=(long)i*NN+j*N+kp, n_km=(long)i*NN+j*N+km;

        const double *u[3]   = {u0, u1, u2};
        const double *v[3]   = {v0, v1, v2};
        const double *tu[3]  = {tu0, tu1, tu2};
        const double *tv[3]  = {tv0, tv1, tv2};
        const double *vu[3]  = {vu0, vu1, vu2};
        const double *vv_[3] = {vv0, vv1, vv2};
        const double *vtu[3] = {vtu0, vtu1, vtu2};
        const double *vtv[3] = {vtv0, vtv1, vtv2};

        double s2[3];
        double s_epk=0, s_etk=0, s_eg=0, s_em=0, s_etg=0, s_etm=0, s_ec=0;
        double s_pm=0, s_trms=0, qp=0, qt=0, rho2=0;

        for (int a = 0; a < 3; a++) {
            double ua=u[a][idx], va=v[a][idx];
            double tua=tu[a][idx], tva=tv[a][idx];
            s2[a] = ua*ua + va*va;
            s_epk += 0.5*(vu[a][idx]*vu[a][idx] + vv_[a][idx]*vv_[a][idx])*dV;
            s_etk += 0.5*(vtu[a][idx]*vtu[a][idx] + vtv[a][idx]*vtv[a][idx])*dV;
            double gx=(u[a][n_ip]-u[a][n_im])*idx1;
            double gy=(u[a][n_jp]-u[a][n_jm])*idx1;
            double gz=(u[a][n_kp]-u[a][n_km])*idx1;
            double hx=(v[a][n_ip]-v[a][n_im])*idx1;
            double hy=(v[a][n_jp]-v[a][n_jm])*idx1;
            double hz=(v[a][n_kp]-v[a][n_km])*idx1;
            s_eg += 0.5*(gx*gx+gy*gy+gz*gz+hx*hx+hy*hy+hz*hz)*dV;
            s_em += 0.5*mass2*(ua*ua+va*va)*dV;
            double tgx=(tu[a][n_ip]-tu[a][n_im])*idx1;
            double tgy=(tu[a][n_jp]-tu[a][n_jm])*idx1;
            double tgz=(tu[a][n_kp]-tu[a][n_km])*idx1;
            double thx=(tv[a][n_ip]-tv[a][n_im])*idx1;
            double thy=(tv[a][n_jp]-tv[a][n_jm])*idx1;
            double thz=(tv[a][n_kp]-tv[a][n_km])*idx1;
            s_etg += 0.5*(tgx*tgx+tgy*tgy+tgz*tgz+thx*thx+thy*thy+thz*thz)*dV;
            s_etm += 0.5*mtheta2*(tua*tua+tva*tva)*dV;
            /* phi_max as complex modulus |Phi_a| (time-stationary on a Q-ball) */
            double mod = sqrt(s2[a]); if (mod > s_pm) s_pm = mod;
            /* E_coupling = -eta * (u.curl(tu) + v.curl(tv)) */
            double ctu, ctv;
            if (a==0) {
                ctu = (tu[2][n_jp]-tu[2][n_jm]-tu[1][n_kp]+tu[1][n_km])*idx1;
                ctv = (tv[2][n_jp]-tv[2][n_jm]-tv[1][n_kp]+tv[1][n_km])*idx1;
            } else if (a==1) {
                ctu = (tu[0][n_kp]-tu[0][n_km]-tu[2][n_ip]+tu[2][n_im])*idx1;
                ctv = (tv[0][n_kp]-tv[0][n_km]-tv[2][n_ip]+tv[2][n_im])*idx1;
            } else {
                ctu = (tu[1][n_ip]-tu[1][n_im]-tu[0][n_jp]+tu[0][n_jm])*idx1;
                ctv = (tv[1][n_ip]-tv[1][n_im]-tv[0][n_jp]+tv[0][n_jm])*idx1;
            }
            s_ec -= eta*(ua*ctu + va*ctv)*dV;
            s_trms += tua*tua + tva*tva;
            /* Noether charges (THEORY §2): Q > 0 for Phi ~ e^{+i omega t} */
            qp += ua*vv_[a][idx] - va*vu[a][idx];
            qt += tua*vtv[a][idx] - tva*vtu[a][idx];
            rho2 += s2[a];
        }
        double sv = s2[0]*s2[1]*s2[2];
        double x = -d_L + i*d_dx, y = -d_L + j*d_dx, z = -d_L + k*d_dx;

        local[0]  = s_epk;
        local[1]  = s_etk;
        local[2]  = s_eg;
        local[3]  = s_em;
        local[4]  = (mu/2.0)*sv/(1.0+kappa*sv)*dV;  /* Vt(s) */
        local[5]  = s_etg;
        local[6]  = s_etm;
        local[7]  = s_ec;
        local[8]  = s_pm;          /* phi_max — max reduction */
        local[9]  = sqrt(sv);      /* P_max = sqrt(s) — max reduction */
        local[10] = sqrt(sv)*dV;   /* P_int */
        local[11] = s_trms;        /* theta_rms sum (6 components) */
        local[12] = qp;            /* Q_phi (raw) */
        local[13] = qt;            /* Q_theta (raw) */
        local[14] = sv;            /* s_max — max reduction */
        local[15] = rho2;          /* M0 (raw) */
        local[16] = x*rho2;
        local[17] = y*rho2;
        local[18] = z*rho2;
        /* v67 Q_core: total Noether charge density masked to the interior ball */
        local[19] = (x*x + y*y + z*z <= qr2) ? (qp + qt) : 0.0;
    }

    for (int vv = 0; vv < CDIAG_NVALS; vv++)
        sdata[vv * blockDim.x + tid] = local[vv];
    __syncthreads();

    for (int s = blockDim.x/2; s > 0; s >>= 1) {
        if (tid < s) {
            for (int vv = 0; vv < CDIAG_NVALS; vv++) {
                int si = vv * blockDim.x;
                if (vv == 8 || vv == 9 || vv == 14) {
                    double other = sdata[si + tid + s];
                    if (other > sdata[si + tid]) sdata[si + tid] = other;
                } else {
                    sdata[si + tid] += sdata[si + tid + s];
                }
            }
        }
        __syncthreads();
    }

    if (tid == 0) {
        for (int vv = 0; vv < CDIAG_NVALS; vv++) {
            if (vv == 8 || vv == 9 || vv == 14) {
                unsigned long long *addr = (unsigned long long *)&d_results[vv];
                unsigned long long old_val = *addr, assumed;
                double new_val = sdata[vv * blockDim.x];
                do {
                    assumed = old_val;
                    double cur = __longlong_as_double(assumed);
                    if (new_val <= cur) break;
                    unsigned long long new_bits = __double_as_longlong(new_val);
                    old_val = atomicCAS(addr, assumed, new_bits);
                } while (assumed != old_val);
            } else {
                atomicAdd(&d_results[vv], sdata[vv * blockDim.x]);
            }
        }
    }
}

/* ================================================================
   v69 gauged complex diagnostics reduction (mirrors CPU
   compute_energy_complex_gauge): identical 20-slot layout to
   reduce_diagnostics_complex_kernel, but E_grad [2], E_tgrad [5] and
   E_coupling [7] use covariant CENTRAL differences (link-transported
   neighbors). E_em is NOT accumulated here — it comes from
   reduce_gauss_kernel (identical formula, computed once) and is added
   to E_total on the host. Never launched when g_gauge==0.
   ================================================================ */
__global__ void reduce_diagnostics_complex_gauge_kernel(
    const double *u0, const double *u1, const double *u2,
    const double *v0, const double *v1, const double *v2,
    const double *tu0, const double *tu1, const double *tu2,
    const double *tv0, const double *tv1, const double *tv2,
    const double *vu0, const double *vu1, const double *vu2,
    const double *vv0, const double *vv1, const double *vv2,
    const double *vtu0, const double *vtu1, const double *vtu2,
    const double *vtv0, const double *vtv1, const double *vtv2,
    const double *th0, const double *th1, const double *th2,
    double dV, double idx1,
    double mass2, double mtheta2, double eta, double mu, double kappa,
    double qr2, double *d_results)
{
    __shared__ double sdata[CDIAG_NVALS * 256];

    long idx = (long)blockIdx.x * blockDim.x + threadIdx.x;
    int tid = threadIdx.x;

    double local[CDIAG_NVALS];
    for (int vv = 0; vv < CDIAG_NVALS; vv++) local[vv] = 0.0;

    if (idx < d_N3) {
        int N = d_N, NN = d_NN;
        int i = (int)(idx / NN), j = (int)((idx / N) % N), k = (int)(idx % N);
        long np[3], nm[3];
        np[0]=(long)((i+1)%N)*NN+(long)j*N+k;   nm[0]=(long)((i-1+N)%N)*NN+(long)j*N+k;
        np[1]=(long)i*NN+(long)((j+1)%N)*N+k;   nm[1]=(long)i*NN+(long)((j-1+N)%N)*N+k;
        np[2]=(long)i*NN+(long)j*N+(k+1)%N;     nm[2]=(long)i*NN+(long)j*N+(k-1+N)%N;

        const double *u[3]   = {u0, u1, u2};
        const double *v[3]   = {v0, v1, v2};
        const double *tu[3]  = {tu0, tu1, tu2};
        const double *tv[3]  = {tv0, tv1, tv2};
        const double *vu[3]  = {vu0, vu1, vu2};
        const double *vv_[3] = {vv0, vv1, vv2};
        const double *vtu[3] = {vtu0, vtu1, vtu2};
        const double *vtv[3] = {vtv0, vtv1, vtv2};
        const double *th[3]  = {th0, th1, th2};

        double fu[6], fv[6];
        for (int a=0;a<3;a++) {
            fu[a]   = u[a][idx];   fv[a]   = v[a][idx];
            fu[3+a] = tu[a][idx];  fv[3+a] = tv[a][idx];
        }
        /* covariant central differences (link cos/sin inline) */
        double DcU[3][6], DcV[3][6];
        for (int d=0; d<3; d++) {
            double cP=cos(th[d][idx]),   sP=sin(th[d][idx]);
            double cM=cos(th[d][nm[d]]), sM=sin(th[d][nm[d]]);
            for (int a=0;a<3;a++) {
                double upn=u[a][np[d]],   vpn=v[a][np[d]];
                double umn=u[a][nm[d]],   vmn=v[a][nm[d]];
                double tpn=tu[a][np[d]],  wpn=tv[a][np[d]];
                double tmn=tu[a][nm[d]],  wmn=tv[a][nm[d]];
                DcU[d][a]   = ((cP*upn - sP*vpn) - (cM*umn + sM*vmn))*idx1;
                DcV[d][a]   = ((cP*vpn + sP*upn) - (cM*vmn - sM*umn))*idx1;
                DcU[d][3+a] = ((cP*tpn - sP*wpn) - (cM*tmn + sM*wmn))*idx1;
                DcV[d][3+a] = ((cP*wpn + sP*tpn) - (cM*wmn - sM*tmn))*idx1;
            }
        }

        double s2[3];
        double s_epk=0, s_etk=0, s_eg=0, s_em=0, s_etg=0, s_etm=0, s_ec=0;
        double s_pm=0, s_trms=0, qp=0, qt=0, rho2=0;
        const int ci1[3]={1,2,0}, ci2[3]={2,0,1};

        for (int a = 0; a < 3; a++) {
            s2[a] = fu[a]*fu[a] + fv[a]*fv[a];
            s_epk += 0.5*(vu[a][idx]*vu[a][idx] + vv_[a][idx]*vv_[a][idx])*dV;
            s_etk += 0.5*(vtu[a][idx]*vtu[a][idx] + vtv[a][idx]*vtv[a][idx])*dV;
            s_eg += 0.5*(DcU[0][a]*DcU[0][a]+DcU[1][a]*DcU[1][a]+DcU[2][a]*DcU[2][a]
                        +DcV[0][a]*DcV[0][a]+DcV[1][a]*DcV[1][a]+DcV[2][a]*DcV[2][a])*dV;
            s_em += 0.5*mass2*s2[a]*dV;
            s_etg += 0.5*(DcU[0][3+a]*DcU[0][3+a]+DcU[1][3+a]*DcU[1][3+a]+DcU[2][3+a]*DcU[2][3+a]
                         +DcV[0][3+a]*DcV[0][3+a]+DcV[1][3+a]*DcV[1][3+a]+DcV[2][3+a]*DcV[2][3+a])*dV;
            s_etm += 0.5*mtheta2*(fu[3+a]*fu[3+a]+fv[3+a]*fv[3+a])*dV;
            /* phi_max as complex modulus |Phi_a| */
            double mod = sqrt(s2[a]); if (mod > s_pm) s_pm = mod;
            /* E_coupling = -eta*(u.Re(DxTheta) + v.Im(DxTheta)) (transported curls) */
            double reDxT = DcU[ci1[a]][3+ci2[a]] - DcU[ci2[a]][3+ci1[a]];
            double imDxT = DcV[ci1[a]][3+ci2[a]] - DcV[ci2[a]][3+ci1[a]];
            s_ec -= eta*(fu[a]*reDxT + fv[a]*imDxT)*dV;
            s_trms += fu[3+a]*fu[3+a] + fv[3+a]*fv[3+a];
            /* Noether charges (THEORY §2): Q > 0 for Phi ~ e^{+i omega t} */
            qp += fu[a]*vv_[a][idx] - fv[a]*vu[a][idx];
            qt += fu[3+a]*vtv[a][idx] - fv[3+a]*vtu[a][idx];
            rho2 += s2[a];
        }
        double sv = s2[0]*s2[1]*s2[2];
        double x = -d_L + i*d_dx, y = -d_L + j*d_dx, z = -d_L + k*d_dx;

        local[0]  = s_epk;
        local[1]  = s_etk;
        local[2]  = s_eg;
        local[3]  = s_em;
        local[4]  = (mu/2.0)*sv/(1.0+kappa*sv)*dV;  /* Vt(s) */
        local[5]  = s_etg;
        local[6]  = s_etm;
        local[7]  = s_ec;
        local[8]  = s_pm;          /* phi_max — max reduction */
        local[9]  = sqrt(sv);      /* P_max = sqrt(s) — max reduction */
        local[10] = sqrt(sv)*dV;   /* P_int */
        local[11] = s_trms;        /* theta_rms sum (6 components) */
        local[12] = qp;            /* Q_phi (raw) */
        local[13] = qt;            /* Q_theta (raw) */
        local[14] = sv;            /* s_max — max reduction */
        local[15] = rho2;          /* M0 (raw) */
        local[16] = x*rho2;
        local[17] = y*rho2;
        local[18] = z*rho2;
        local[19] = (x*x + y*y + z*z <= qr2) ? (qp + qt) : 0.0;
    }

    for (int vv = 0; vv < CDIAG_NVALS; vv++)
        sdata[vv * blockDim.x + tid] = local[vv];
    __syncthreads();

    for (int s = blockDim.x/2; s > 0; s >>= 1) {
        if (tid < s) {
            for (int vv = 0; vv < CDIAG_NVALS; vv++) {
                int si = vv * blockDim.x;
                if (vv == 8 || vv == 9 || vv == 14) {
                    double other = sdata[si + tid + s];
                    if (other > sdata[si + tid]) sdata[si + tid] = other;
                } else {
                    sdata[si + tid] += sdata[si + tid + s];
                }
            }
        }
        __syncthreads();
    }

    if (tid == 0) {
        for (int vv = 0; vv < CDIAG_NVALS; vv++) {
            if (vv == 8 || vv == 9 || vv == 14) {
                unsigned long long *addr = (unsigned long long *)&d_results[vv];
                unsigned long long old_val = *addr, assumed;
                double new_val = sdata[vv * blockDim.x];
                do {
                    assumed = old_val;
                    double cur = __longlong_as_double(assumed);
                    if (new_val <= cur) break;
                    unsigned long long new_bits = __double_as_longlong(new_val);
                    old_val = atomicCAS(addr, assumed, new_bits);
                } while (assumed != old_val);
            } else {
                atomicAdd(&d_results[vv], sdata[vv * blockDim.x]);
            }
        }
    }
}

/* ================================================================
   v69 Gauss / gauge diagnostics reduction (mirrors CPU compute_gauss):
   6 values into d_out[0..5]:
   [0] gauss_max = max |G - OFF| over interior Omega (max reduction)
   [1] gauss_l2 sum = sum (G-OFF)^2 over Omega (host: sqrt(sum/nOmega))
   [2] nOmega   (voxel count, summed as double)
   [3] E_em     = int [E^2/2 + (1-cos th_P)/(g^2 a^4)] dV
   [4] fsum     = sum divE*dV over the centered Q_flux cube
   [5] nC       (cube voxel count)
   Host post-processing: q_flux = (fsum - OFF*nC*dV)/g  (g != 0).
   Launched whenever complex_gauge=1 (handles GG==0: no magnetic term).
   ================================================================ */
#define GAUSS_NVALS 6

__global__ void reduce_gauss_kernel(
    const double *u0, const double *u1, const double *u2,
    const double *v0, const double *v1, const double *v2,
    const double *tu0, const double *tu1, const double *tu2,
    const double *tv0, const double *tv1, const double *tv2,
    const double *vu0, const double *vu1, const double *vu2,
    const double *vv0, const double *vv1, const double *vv2,
    const double *vtu0, const double *vtu1, const double *vtu2,
    const double *vtv0, const double *vtv1, const double *vtv2,
    const double *E0, const double *E1, const double *E2,
    const double *th0, const double *th1, const double *th2,
    double GG, double OFF, double Rint2, double Rcube, double dV,
    double *d_out)
{
    __shared__ double sdata[GAUSS_NVALS * 256];

    long idx = (long)blockIdx.x * blockDim.x + threadIdx.x;
    int tid = threadIdx.x;

    double local[GAUSS_NVALS];
    for (int vv = 0; vv < GAUSS_NVALS; vv++) local[vv] = 0.0;

    if (idx < d_N3) {
        int N = d_N, NN = d_NN;
        int i = (int)(idx / NN), j = (int)((idx / N) % N), k = (int)(idx % N);
        long np[3], nm[3];
        np[0]=(long)((i+1)%N)*NN+(long)j*N+k;   nm[0]=(long)((i-1+N)%N)*NN+(long)j*N+k;
        np[1]=(long)i*NN+(long)((j+1)%N)*N+k;   nm[1]=(long)i*NN+(long)((j-1+N)%N)*N+k;
        np[2]=(long)i*NN+(long)j*N+(k+1)%N;     nm[2]=(long)i*NN+(long)j*N+(k-1+N)%N;
        double x = -d_L + i*d_dx, y = -d_L + j*d_dx, z = -d_L + k*d_dx;

        double divE = (E0[idx]-E0[nm[0]]
                      +E1[idx]-E1[nm[1]]
                      +E2[idx]-E2[nm[2]]) / d_dx;
        /* rho_Q at idx (Noether density) */
        double rho = (u0[idx]*vv0[idx] - v0[idx]*vu0[idx])
                   + (u1[idx]*vv1[idx] - v1[idx]*vu1[idx])
                   + (u2[idx]*vv2[idx] - v2[idx]*vu2[idx])
                   + (tu0[idx]*vtv0[idx] - tv0[idx]*vtu0[idx])
                   + (tu1[idx]*vtv1[idx] - tv1[idx]*vtu1[idx])
                   + (tu2[idx]*vtv2[idx] - tv2[idx]*vtu2[idx]);
        double G = divE - GG*rho;

        if (x*x + y*y + z*z < Rint2) {
            double d = fabs(G - OFF);
            local[0] = d;                /* gauss_max — max reduction */
            local[1] = (G-OFF)*(G-OFF);  /* gauss_l2 sum */
            local[2] = 1.0;              /* nOmega */
        }
        double ee = 0.5*(E0[idx]*E0[idx] + E1[idx]*E1[idx] + E2[idx]*E2[idx]);
        if (GG != 0) {
            const double *th[3] = {th0, th1, th2};
            const int pa[3]={0,1,2}, pb[3]={1,2,0};
            double bb=0;
            for (int p=0;p<3;p++) {
                int a=pa[p], b=pb[p];
                double ang = th[a][idx]+th[b][np[a]]-th[a][np[b]]-th[b][idx];
                bb += 1.0 - cos(ang);
            }
            ee += bb/(GG*GG*d_dx*d_dx*d_dx*d_dx);
        }
        local[3] = ee*dV;
        if (fabs(x)<=Rcube && fabs(y)<=Rcube && fabs(z)<=Rcube) {
            local[4] = divE*dV;          /* a^3 * (1/a) sum_i (E_i(x)-E_i(x-i)) */
            local[5] = 1.0;              /* nC */
        }
    }

    for (int vv = 0; vv < GAUSS_NVALS; vv++)
        sdata[vv * blockDim.x + tid] = local[vv];
    __syncthreads();

    for (int s = blockDim.x/2; s > 0; s >>= 1) {
        if (tid < s) {
            for (int vv = 0; vv < GAUSS_NVALS; vv++) {
                int si = vv * blockDim.x;
                if (vv == 0) {
                    double other = sdata[si + tid + s];
                    if (other > sdata[si + tid]) sdata[si + tid] = other;
                } else {
                    sdata[si + tid] += sdata[si + tid + s];
                }
            }
        }
        __syncthreads();
    }

    if (tid == 0) {
        for (int vv = 0; vv < GAUSS_NVALS; vv++) {
            if (vv == 0) {
                unsigned long long *addr = (unsigned long long *)&d_out[vv];
                unsigned long long old_val = *addr, assumed;
                double new_val = sdata[vv * blockDim.x];
                do {
                    assumed = old_val;
                    double cur = __longlong_as_double(assumed);
                    if (new_val <= cur) break;
                    unsigned long long new_bits = __double_as_longlong(new_val);
                    old_val = atomicCAS(addr, assumed, new_bits);
                } while (assumed != old_val);
            } else {
                atomicAdd(&d_out[vv], sdata[vv * blockDim.x]);
            }
        }
    }
}

/* v66 second pass of compute_charges: centroid-based rms radius of
 * rho2 = sum_a (u_a^2 + v_a^2). Accumulates the RAW rsum into d_rsum
 * (host applies *dV/M0 and sqrt, exactly like the CPU). */
__global__ void reduce_rcore_kernel(
    const double *u0, const double *u1, const double *u2,
    const double *v0, const double *v1, const double *v2,
    double cx, double cy, double cz, double *d_rsum)
{
    __shared__ double sdata[256];
    long idx = (long)blockIdx.x * blockDim.x + threadIdx.x;
    int tid = threadIdx.x;
    double val = 0.0;
    if (idx < d_N3) {
        int N = d_N, NN = d_NN;
        int i = (int)(idx / NN), j = (int)((idx / N) % N), k = (int)(idx % N);
        double x = -d_L + i*d_dx, y = -d_L + j*d_dx, z = -d_L + k*d_dx;
        double rho2 = u0[idx]*u0[idx] + v0[idx]*v0[idx]
                    + u1[idx]*u1[idx] + v1[idx]*v1[idx]
                    + u2[idx]*u2[idx] + v2[idx]*v2[idx];
        val = ((x-cx)*(x-cx)+(y-cy)*(y-cy)+(z-cz)*(z-cz))*rho2;
    }
    sdata[tid] = val;
    __syncthreads();
    for (int s = blockDim.x/2; s > 0; s >>= 1) {
        if (tid < s) sdata[tid] += sdata[tid + s];
        __syncthreads();
    }
    if (tid == 0) atomicAdd(d_rsum, sdata[0]);
}

/* v67: find the voxel index where s = prod_a(u_a^2+v_a^2) attains the global
 * max found by the main complex reduction (omega_core sample point). The
 * relative tolerance guards 1-ulp cross-kernel rounding differences; ties
 * resolve to the LOWEST index via atomicMin (matches CPU compute_charges).
 * d_idx must be initialized to ~0ULL before launch. */
__global__ void find_smax_idx_kernel(
    const double *u0, const double *u1, const double *u2,
    const double *v0, const double *v1, const double *v2,
    double smax, unsigned long long *d_idx)
{
    long idx = (long)blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= d_N3) return;
    double s2_0 = u0[idx]*u0[idx] + v0[idx]*v0[idx];
    double s2_1 = u1[idx]*u1[idx] + v1[idx]*v1[idx];
    double s2_2 = u2[idx]*u2[idx] + v2[idx]*v2[idx];
    double s = s2_0*s2_1*s2_2;
    if (s >= smax * (1.0 - 1e-12))
        atomicMin(d_idx, (unsigned long long)idx);
}

/* Zero the diagnostics results buffer */
__global__ void zero_diag_kernel(double *d_results, int n) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < n) d_results[idx] = 0.0;
}

/* v65 self-tuning: rescale phi (and its velocity) by f to enforce a fixed charge
   Q=int|phi|^2 (the mechanism-A Lagrange-multiplier projection). */
__global__ void rescale_phi_kernel(double *p0, double *p1, double *p2,
                                   double *v0, double *v1, double *v2,
                                   double f, long N3) {
    long i = (long)blockIdx.x * blockDim.x + threadIdx.x;
    if (i < N3) { p0[i]*=f; p1[i]*=f; p2[i]*=f; v0[i]*=f; v1[i]*=f; v2[i]*=f; }
}

/* ================================================================
   GPU memory management
   ================================================================ */

static double *d_phi[3], *d_vel_phi[3], *d_acc_phi[3];
static double *d_theta[3], *d_vel_theta[3], *d_acc_theta[3];
static double *d_mismatch[3], *d_harden_Q[3];  /* intermediates for Cosserat/hardening */
/* v66 complex sector device arrays (allocated only when complex_phi=1) */
static double *d_phi_im[3], *d_vel_phi_im[3], *d_acc_phi_im[3];
static double *d_theta_im[3], *d_vel_theta_im[3], *d_acc_theta_im[3];
/* v69 gauge sector device arrays (allocated only when complex_gauge=1):
 * +9 f64 blocks — links th_i, E_i, E-kick E_acc_i (CPU layout, no scratch) */
static double *d_th[3], *d_E[3], *d_E_acc[3];
static int gpu_blocks;
static int gpu_has_intermediates = 0;  /* nonzero if intermediate arrays allocated */
static int gpu_complex_mode = 0;       /* nonzero if complex arrays allocated (v66) */
static int gpu_gauge_mode = 0;         /* nonzero if gauge arrays allocated (v69) */
static double gpu_g_gauge = 0.0;       /* cached coupling for verlet dispatch */
static double gpu_dx_cached = 0.0;     /* cached dx for the link-drift prefactor */

static void gpu_alloc(long N3, int need_intermediates, int complex_mode, int gauge_mode) {
    size_t bytes = N3 * sizeof(double);
    gpu_complex_mode = complex_mode;
    gpu_gauge_mode = gauge_mode;
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
    if (complex_mode) {
        /* 18 more arrays: the imaginary copy of the 6-field system */
        for (int a = 0; a < 3; a++) {
            cudaMalloc(&d_phi_im[a], bytes);   cudaMalloc(&d_vel_phi_im[a], bytes);   cudaMalloc(&d_acc_phi_im[a], bytes);
            cudaMalloc(&d_theta_im[a], bytes); cudaMalloc(&d_vel_theta_im[a], bytes); cudaMalloc(&d_acc_theta_im[a], bytes);
        }
        total_gb += 18.0*bytes/1e9;
    } else {
        for (int a = 0; a < 3; a++) {
            d_phi_im[a] = d_vel_phi_im[a] = d_acc_phi_im[a] = NULL;
            d_theta_im[a] = d_vel_theta_im[a] = d_acc_theta_im[a] = NULL;
        }
    }
    if (gauge_mode) {
        /* v69: +9 f64 blocks (th_i, E_i, E_acc_i), zero-initialized */
        for (int a = 0; a < 3; a++) {
            cudaMalloc(&d_th[a], bytes);
            cudaMalloc(&d_E[a], bytes);
            cudaMalloc(&d_E_acc[a], bytes);
            cudaMemset(d_th[a], 0, bytes);
            cudaMemset(d_E[a], 0, bytes);
            cudaMemset(d_E_acc[a], 0, bytes);
        }
        total_gb += 9.0*bytes/1e9;
    } else {
        for (int a = 0; a < 3; a++) d_th[a] = d_E[a] = d_E_acc[a] = NULL;
    }
    gpu_blocks = (int)((N3 + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK);
    printf("GPU: allocated %.2f GB (physics%s%s%s), %d blocks × %d threads\n",
           total_gb, need_intermediates ? "+Cosserat" : "",
           complex_mode ? ", complex 36 arrays" : "",
           gauge_mode ? " + U(1) links 9 arrays" : "",
           gpu_blocks, THREADS_PER_BLOCK);
    if (complex_mode) {
        /* v66 memory guard: 36 N^3 doubles of physics, plus snapshot staging
         * (24 N^3 f32/f16). Physics alone: N=320 -> 9.4 GB, N=384 -> 16.3 GB,
         * N=400 -> 18.4 GB. Practical limits: N<=320 on 16 GB, N<=400 on 32 GB.
         * v69 gauged: +9 blocks -> 45 N^3 doubles (N=320 -> 11.8 GB,
         * N=384 -> 20.4 GB) plus 30-column staging; N<=288 on 16 GB,
         * N<=384 on 32 GB. */
        size_t free_b = 0, total_b = 0;
        cudaMemGetInfo(&free_b, &total_b);
        int nphys = gauge_mode ? 45 : 36;
        double need_gb = (double)nphys*bytes/1e9;
        if (gauge_mode)
            printf("Gauged mode: %.2f GB physics (45 N^3 arrays); max N<=288 on 16 GB, N<=384 on 32 GB\n",
                   need_gb);
        else
            printf("Complex mode: %.2f GB physics (36 N^3 arrays); max N<=320 on 16 GB, N<=400 on 32 GB\n",
                   need_gb);
        if (need_gb*1e9 > (double)total_b)
            printf("WARNING: %s physics arrays (%.2f GB) exceed GPU memory (%.2f GB)\n",
                   gauge_mode ? "gauged" : "complex", need_gb, total_b/1e9);
    }
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
    if (gpu_complex_mode) {
        for (int a = 0; a < 3; a++) {
            cudaMemcpy(d_phi_im[a], g->phi_im[a], bytes, cudaMemcpyHostToDevice);
            cudaMemcpy(d_vel_phi_im[a], g->phi_im_vel[a], bytes, cudaMemcpyHostToDevice);
            cudaMemcpy(d_theta_im[a], g->theta_im[a], bytes, cudaMemcpyHostToDevice);
            cudaMemcpy(d_vel_theta_im[a], g->theta_im_vel[a], bytes, cudaMemcpyHostToDevice);
            cudaMemset(d_acc_phi_im[a], 0, bytes);
            cudaMemset(d_acc_theta_im[a], 0, bytes);
        }
    }
    if (gpu_gauge_mode) {
        /* v69: upload host-side init (post Gauss-projection) links + E */
        for (int a = 0; a < 3; a++) {
            cudaMemcpy(d_th[a], g->th[a], bytes, cudaMemcpyHostToDevice);
            cudaMemcpy(d_E[a], g->Efield[a], bytes, cudaMemcpyHostToDevice);
            cudaMemset(d_E_acc[a], 0, bytes);
        }
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
    if (gpu_complex_mode) {
        for (int a = 0; a < 3; a++) {
            cudaMemcpy(g->phi_im[a], d_phi_im[a], bytes, cudaMemcpyDeviceToHost);
            cudaMemcpy(g->phi_im_vel[a], d_vel_phi_im[a], bytes, cudaMemcpyDeviceToHost);
            cudaMemcpy(g->theta_im[a], d_theta_im[a], bytes, cudaMemcpyDeviceToHost);
            cudaMemcpy(g->theta_im_vel[a], d_vel_theta_im[a], bytes, cudaMemcpyDeviceToHost);
        }
    }
    if (gpu_gauge_mode) {
        for (int a = 0; a < 3; a++) {
            cudaMemcpy(g->th[a], d_th[a], bytes, cudaMemcpyDeviceToHost);
            cudaMemcpy(g->Efield[a], d_E[a], bytes, cudaMemcpyDeviceToHost);
        }
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
    if (gpu_complex_mode) {
        for (int a = 0; a < 3; a++) {
            cudaFree(d_phi_im[a]); cudaFree(d_vel_phi_im[a]); cudaFree(d_acc_phi_im[a]);
            cudaFree(d_theta_im[a]); cudaFree(d_vel_theta_im[a]); cudaFree(d_acc_theta_im[a]);
        }
    }
    if (gpu_gauge_mode) {
        for (int a = 0; a < 3; a++) {
            cudaFree(d_th[a]); cudaFree(d_E[a]); cudaFree(d_E_acc[a]);
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
    cudaMemcpyToSymbol(d_THETA_SAT, &c->theta_sat, sizeof(double));
    cudaMemcpyToSymbol(d_GAMMA_CONV, &c->gamma_conv, sizeof(double));
    cudaMemcpyToSymbol(d_SIGMA_GRAD, &c->sigma_grad, sizeof(double));
    cudaMemcpyToSymbol(d_CHI_CHIRAL, &c->chi_chiral, sizeof(double));
    cudaMemcpyToSymbol(d_SIGMA_CUBIC, &c->sigma_cubic, sizeof(double));
    cudaMemcpyToSymbol(d_SIGMA_FREQ, &c->sigma_freq, sizeof(double));
    cudaMemcpyToSymbol(d_SIGMA_CROSS, &c->sigma_cross, sizeof(double));
    cudaMemcpyToSymbol(d_LAMBDA_SELF, &c->lambda_self, sizeof(double));
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
    /* v69 §1.3/§3.5: gauge sector active only when complex_gauge && g != 0 */
    const int gauged = gpu_gauge_mode && gpu_g_gauge != 0.0;
    /* stage 1: half-kick (matter velocities AND E) */
    verlet_halfkick_kernel<<<gpu_blocks, THREADS_PER_BLOCK>>>(
        d_vel_phi[0], d_vel_phi[1], d_vel_phi[2],
        d_vel_theta[0], d_vel_theta[1], d_vel_theta[2],
        d_acc_phi[0], d_acc_phi[1], d_acc_phi[2],
        d_acc_theta[0], d_acc_theta[1], d_acc_theta[2], hdt);
    if (gpu_complex_mode)
        verlet_halfkick_kernel<<<gpu_blocks, THREADS_PER_BLOCK>>>(
            d_vel_phi_im[0], d_vel_phi_im[1], d_vel_phi_im[2],
            d_vel_theta_im[0], d_vel_theta_im[1], d_vel_theta_im[2],
            d_acc_phi_im[0], d_acc_phi_im[1], d_acc_phi_im[2],
            d_acc_theta_im[0], d_acc_theta_im[1], d_acc_theta_im[2], hdt);
    if (gauged)
        gauge_halfkick_kernel<<<gpu_blocks, THREADS_PER_BLOCK>>>(
            d_E[0], d_E[1], d_E[2],
            d_E_acc[0], d_E_acc[1], d_E_acc[2], hdt);
    /* stage 2: drift (fields AND link angles, wrapped to (-pi,pi]) */
    verlet_drift_kernel<<<gpu_blocks, THREADS_PER_BLOCK>>>(
        d_phi[0], d_phi[1], d_phi[2],
        d_theta[0], d_theta[1], d_theta[2],
        d_vel_phi[0], d_vel_phi[1], d_vel_phi[2],
        d_vel_theta[0], d_vel_theta[1], d_vel_theta[2], dt);
    if (gpu_complex_mode)
        verlet_drift_kernel<<<gpu_blocks, THREADS_PER_BLOCK>>>(
            d_phi_im[0], d_phi_im[1], d_phi_im[2],
            d_theta_im[0], d_theta_im[1], d_theta_im[2],
            d_vel_phi_im[0], d_vel_phi_im[1], d_vel_phi_im[2],
            d_vel_theta_im[0], d_vel_theta_im[1], d_vel_theta_im[2], dt);
    if (gauged) {
        const double gad = -gpu_g_gauge * gpu_dx_cached * dt;  /* th_dot = -g*a*E */
        gauge_link_drift_kernel<<<gpu_blocks, THREADS_PER_BLOCK>>>(
            d_th[0], d_th[1], d_th[2],
            d_E[0], d_E[1], d_E[2], gad);
    }
    /* stage 3: forces (matter accs + E_acc) */
    if (gauged) {
        const double GG  = gpu_g_gauge;
        const double STP = 1.0 / (GG * gpu_dx_cached * gpu_dx_cached * gpu_dx_cached);
        const double inva = 1.0 / gpu_dx_cached;
        compute_forces_complex_gauge_kernel<<<gpu_blocks, THREADS_PER_BLOCK>>>(
            d_phi[0], d_phi[1], d_phi[2],
            d_phi_im[0], d_phi_im[1], d_phi_im[2],
            d_theta[0], d_theta[1], d_theta[2],
            d_theta_im[0], d_theta_im[1], d_theta_im[2],
            d_th[0], d_th[1], d_th[2],
            d_acc_phi[0], d_acc_phi[1], d_acc_phi[2],
            d_acc_phi_im[0], d_acc_phi_im[1], d_acc_phi_im[2],
            d_acc_theta[0], d_acc_theta[1], d_acc_theta[2],
            d_acc_theta_im[0], d_acc_theta_im[1], d_acc_theta_im[2],
            d_E_acc[0], d_E_acc[1], d_E_acc[2],
            GG, STP, inva);
    } else if (gpu_complex_mode) {
        /* v66: 12-field minimal loop; intermediates guarded off by cfg_validate */
        compute_forces_complex_kernel<<<gpu_blocks, THREADS_PER_BLOCK>>>(
            d_phi[0], d_phi[1], d_phi[2],
            d_phi_im[0], d_phi_im[1], d_phi_im[2],
            d_theta[0], d_theta[1], d_theta[2],
            d_theta_im[0], d_theta_im[1], d_theta_im[2],
            d_acc_phi[0], d_acc_phi[1], d_acc_phi[2],
            d_acc_phi_im[0], d_acc_phi_im[1], d_acc_phi_im[2],
            d_acc_theta[0], d_acc_theta[1], d_acc_theta[2],
            d_acc_theta_im[0], d_acc_theta_im[1], d_acc_theta_im[2]);
    } else {
        /* Intermediates (Cosserat mismatch + hardening Q) */
        if (gpu_has_intermediates) {
            compute_intermediates_kernel<<<gpu_blocks, THREADS_PER_BLOCK>>>(
                d_phi[0], d_phi[1], d_phi[2],
                d_theta[0], d_theta[1], d_theta[2],
                d_mismatch[0], d_mismatch[1], d_mismatch[2],
                d_harden_Q[0], d_harden_Q[1], d_harden_Q[2]);
        }
        compute_forces_kernel<<<gpu_blocks, THREADS_PER_BLOCK>>>(
            d_phi[0], d_phi[1], d_phi[2],
            d_theta[0], d_theta[1], d_theta[2],
            d_vel_theta[0], d_vel_theta[1], d_vel_theta[2],
            d_acc_phi[0], d_acc_phi[1], d_acc_phi[2],
            d_acc_theta[0], d_acc_theta[1], d_acc_theta[2],
            d_mismatch[0], d_mismatch[1], d_mismatch[2],
            d_harden_Q[0], d_harden_Q[1], d_harden_Q[2]);
    }
    /* stage 4: half-kick (matter + E) */
    verlet_halfkick_kernel<<<gpu_blocks, THREADS_PER_BLOCK>>>(
        d_vel_phi[0], d_vel_phi[1], d_vel_phi[2],
        d_vel_theta[0], d_vel_theta[1], d_vel_theta[2],
        d_acc_phi[0], d_acc_phi[1], d_acc_phi[2],
        d_acc_theta[0], d_acc_theta[1], d_acc_theta[2], hdt);
    if (gpu_complex_mode)
        verlet_halfkick_kernel<<<gpu_blocks, THREADS_PER_BLOCK>>>(
            d_vel_phi_im[0], d_vel_phi_im[1], d_vel_phi_im[2],
            d_vel_theta_im[0], d_vel_theta_im[1], d_vel_theta_im[2],
            d_acc_phi_im[0], d_acc_phi_im[1], d_acc_phi_im[2],
            d_acc_theta_im[0], d_acc_theta_im[1], d_acc_theta_im[2], hdt);
    if (gauged)
        gauge_halfkick_kernel<<<gpu_blocks, THREADS_PER_BLOCK>>>(
            d_E[0], d_E[1], d_E[2],
            d_E_acc[0], d_E_acc[1], d_E_acc[2], hdt);
    /* stage 5: boundary — type set externally before calling gpu_verlet_step */
    absorbing_boundary_kernel<<<gpu_blocks, THREADS_PER_BLOCK>>>(
        d_vel_phi[0], d_vel_phi[1], d_vel_phi[2],
        d_vel_theta[0], d_vel_theta[1], d_vel_theta[2]);
    if (gpu_complex_mode) {
        absorbing_boundary_kernel<<<gpu_blocks, THREADS_PER_BLOCK>>>(
            d_vel_phi_im[0], d_vel_phi_im[1], d_vel_phi_im[2],
            d_vel_theta_im[0], d_vel_theta_im[1], d_vel_theta_im[2]);
        /* v69 §3.6: damp E in the sponge (links NOT damped) */
        if (gpu_gauge_mode)
            absorbing_boundary_gauge_kernel<<<gpu_blocks, THREADS_PER_BLOCK>>>(
                d_E[0], d_E[1], d_E[2]);
    }
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
   FrameWriter — thread-safe wrapper for SFA writes.
   Both snap_hook (async writer thread) and vec_hook (main thread)
   write through this. The mutex serializes all writes to one file.
   ================================================================ */

typedef struct FrameWriter_ {
    SFA *sfa;
    pthread_mutex_t mutex;  /* only held during fseek+fwrite+index — NOT during compression */
} FrameWriter;

static void fw_init(FrameWriter *fw, SFA *sfa) {
    fw->sfa = sfa;
    pthread_mutex_init(&fw->mutex, NULL);
}

static void fw_destroy(FrameWriter *fw) {
    pthread_mutex_destroy(&fw->mutex);
}

/* Write a pre-compressed chunk to the SFA file under lock.
 * chunk_type: SFA_CHUNK_FRMD or SFA_CHUNK_FRVD
 * frame_type: SFA_FRAME_VOXEL, SFA_FRAME_VEC_I, SFA_FRAME_VEC_P
 * The caller is responsible for compression — this only does file I/O. */
static int fw_write_raw(FrameWriter *fw, uint32_t chunk_type, uint32_t frame_type,
                         double time, const void *compressed, uint64_t comp_size,
                         uint32_t checksum) {
    pthread_mutex_lock(&fw->mutex);

    SFA *s = fw->sfa;
    FILE *fp = s->fp;

    /* Append chunk: FRMD = type(4)+size(8)+data; FRVD = type(4)+size(8)+time(8)+data */
    fseek(fp, 0, SEEK_END);
    uint64_t frame_offset = ftell(fp);
    if (chunk_type == SFA_CHUNK_FRVD) {
        uint64_t chunk_size = 8 + comp_size;  /* time + data */
        fwrite(&chunk_type, 4, 1, fp);
        fwrite(&chunk_size, 8, 1, fp);
        fwrite(&time, 8, 1, fp);
        fwrite(compressed, 1, comp_size, fp);
    } else {
        uint64_t chunk_size = 12 + comp_size;  /* FRMD includes header in size */
        fwrite(&chunk_type, 4, 1, fp);
        fwrite(&chunk_size, 8, 1, fp);
        fwrite(compressed, 1, comp_size, fp);
    }

    /* Update JMPF index */
    if (s->cur_jmpf_offset) {
        long jmpf_entry = (long)s->cur_jmpf_offset + 12 + 8 + s->cur_jmpf_entries * 32;
        fseek(fp, jmpf_entry, SEEK_SET);
        fwrite(&time, 8, 1, fp);
        fwrite(&frame_offset, 8, 1, fp);
        fwrite(&comp_size, 8, 1, fp);
        fwrite(&checksum, 4, 1, fp);
        fwrite(&frame_type, 4, 1, fp);

        s->cur_jmpf_entries++;
        fseek(fp, (long)s->cur_jmpf_offset + 12 + 4, SEEK_SET);
        fwrite(&s->cur_jmpf_entries, 4, 1, fp);

        long jtop_fc = (long)s->cur_jtop_offset + 12 + 16 + s->cur_l1_index * 16 + 12;
        fseek(fp, jtop_fc, SEEK_SET);
        fwrite(&s->cur_jmpf_entries, 4, 1, fp);

        s->total_frames++;
        fseek(fp, 0, SEEK_END);
    }

    pthread_mutex_unlock(&fw->mutex);
    return 0;
}

/* High-level: compress + write a voxel frame.
 * Compression happens OUTSIDE the lock — only fwrite under lock. */
static int fw_write_voxel(FrameWriter *fw, double time, void **cols) {
    SFA *s = fw->sfa;

    /* Phase 1: Assemble + compress OUTSIDE the lock */
    uint8_t *raw = (uint8_t*)malloc(s->frame_bytes);
    uint64_t off = 0;
    for (uint32_t c = 0; c < s->n_columns; c++) {
        uint64_t col_bytes = s->N_total * sfa_dtype_size[s->columns[c].dtype];
        memcpy(raw + off, cols[c], col_bytes);
        off += col_bytes;
    }
    uint32_t checksum = sfa_crc32(raw, s->frame_bytes);

    /* Per-column BSS+zstd compression (parallel, no lock) */
    uint32_t nc = s->n_columns;
    uint64_t *col_comp_sizes = (uint64_t*)calloc(nc, sizeof(uint64_t));
    uint8_t **col_comp_data = (uint8_t**)calloc(nc, sizeof(uint8_t*));
    uint64_t *col_offsets = (uint64_t*)calloc(nc, sizeof(uint64_t));
    uint64_t *col_bytesv = (uint64_t*)calloc(nc, sizeof(uint64_t));
    { uint64_t o2 = 0; for (uint32_t c = 0; c < nc; c++) {
        col_offsets[c] = o2; col_bytesv[c] = s->N_total * sfa_dtype_size[s->columns[c].dtype]; o2 += col_bytesv[c];
    }}

    #ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic)
    #endif
    for (uint32_t c = 0; c < nc; c++) {
        int es = sfa_dtype_size[s->columns[c].dtype];
        uint64_t col_bytes = col_bytesv[c];
        uint8_t *bss_col = (uint8_t*)malloc(col_bytes);
        uint8_t *src_col = raw + col_offsets[c];
        for (int b = 0; b < es; b++)
            for (uint64_t v = 0; v < s->N_total; v++)
                bss_col[b * s->N_total + v] = src_col[v * es + b];
        size_t bound = ZSTD_compressBound(col_bytes);
        col_comp_data[c] = (uint8_t*)malloc(bound);
        col_comp_sizes[c] = ZSTD_compress(col_comp_data[c], bound, bss_col, col_bytes, 3);
        free(bss_col);
    }
    free(col_offsets); free(col_bytesv);

    /* Assemble compressed payload */
    uint64_t total_comp = 4 + nc * 8;
    for (uint32_t c = 0; c < nc; c++) total_comp += col_comp_sizes[c];
    void *comp = malloc(total_comp);
    uint8_t *p = (uint8_t*)comp;
    memcpy(p, &nc, 4); p += 4;
    memcpy(p, col_comp_sizes, nc * 8); p += nc * 8;
    for (uint32_t c = 0; c < nc; c++) {
        memcpy(p, col_comp_data[c], col_comp_sizes[c]); p += col_comp_sizes[c];
        free(col_comp_data[c]);
    }
    free(col_comp_sizes); free(col_comp_data); free(raw);

    /* Phase 2: Write under lock — just fseek + fwrite + index update */
    int ret = fw_write_raw(fw, SFA_CHUNK_FRMD, SFA_FRAME_VOXEL, time, comp, total_comp, checksum);
    free(comp);
    return ret;
}

/* High-level: compress + write a vector I-frame.
 * Compression happens OUTSIDE the lock, only fwrite under lock. */
static int fw_write_vec_iframe(FrameWriter *fw, double time,
                                uint32_t n_patches, uint8_t block_size, uint16_t n_coeffs,
                                const int16_t *origins, const float *coeffs) {
    /* Build payload */
    uint8_t pad = 0;
    uint64_t origin_bytes = (uint64_t)n_patches * 6;
    uint64_t coeff_bytes = (uint64_t)n_patches * n_coeffs * 4;
    uint64_t payload_size = 8 + origin_bytes + coeff_bytes;

    uint8_t *payload = (uint8_t*)malloc(payload_size);
    uint64_t off = 0;
    memcpy(payload + off, &n_patches, 4); off += 4;
    memcpy(payload + off, &block_size, 1); off += 1;
    memcpy(payload + off, &n_coeffs, 2); off += 2;
    memcpy(payload + off, &pad, 1); off += 1;
    memcpy(payload + off, origins, origin_bytes); off += origin_bytes;
    memcpy(payload + off, coeffs, coeff_bytes);

    /* Compress OUTSIDE the lock */
    uint32_t checksum = sfa_crc32(payload, payload_size);
    size_t bound = ZSTD_compressBound(payload_size);
    void *comp = malloc(bound);
    size_t comp_size = ZSTD_compress(comp, bound, payload, payload_size, 3);
    free(payload);

    if (ZSTD_isError(comp_size)) { free(comp); return -1; }

    /* Only the write is under lock */
    int ret = fw_write_raw(fw, SFA_CHUNK_FRVD, SFA_FRAME_VEC_I, time, comp, comp_size, checksum);
    free(comp);
    return ret;
}

/* High-level: compress + write a vector P-frame.
 * Compression OUTSIDE lock. */
static int fw_write_vec_pframe(FrameWriter *fw, double time,
                                uint32_t n_deltas, uint16_t n_coeffs,
                                const uint32_t *indices, const float *delta_coeffs) {
    uint16_t pad = 0;
    uint64_t idx_bytes = (uint64_t)n_deltas * 4;
    uint64_t coeff_bytes = (uint64_t)n_deltas * n_coeffs * 4;
    uint64_t payload_size = 8 + idx_bytes + coeff_bytes;

    uint8_t *payload = (uint8_t*)malloc(payload_size);
    uint64_t off = 0;
    memcpy(payload + off, &n_deltas, 4); off += 4;
    memcpy(payload + off, &n_coeffs, 2); off += 2;
    memcpy(payload + off, &pad, 2); off += 2;
    if (n_deltas > 0 && indices && delta_coeffs) {
        memcpy(payload + off, indices, idx_bytes); off += idx_bytes;
        memcpy(payload + off, delta_coeffs, coeff_bytes);
    }

    uint32_t checksum = sfa_crc32(payload, payload_size);
    size_t bound = ZSTD_compressBound(payload_size);
    void *comp = malloc(bound);
    size_t comp_size = ZSTD_compress(comp, bound, payload, payload_size, 3);
    free(payload);

    if (ZSTD_isError(comp_size)) { free(comp); return -1; }

    int ret = fw_write_raw(fw, SFA_CHUNK_FRVD, SFA_FRAME_VEC_P, time, comp, comp_size, checksum);
    free(comp);
    return ret;
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

        /* Compress and write frame to SFA (nf = 12 real, 24 complex, 30 gauged) */
        void *cols[30];
        long N3 = ctx->N3;
        int nf = ctx->nf;
        if (ctx->precision == 0) {
            /* f16: h_pin_buf is uint16_t[nf*N3] */
            uint16_t *base = (uint16_t *)ctx->h_pin_buf;
            for (int c = 0; c < nf; c++)
                cols[c] = base + c * N3;
        } else if (ctx->precision == 1) {
            /* f32: h_pin_buf is float[nf*N3] */
            float *base = (float *)ctx->h_pin_buf;
            for (int c = 0; c < nf; c++)
                cols[c] = base + c * N3;
        } else {
            /* f64: h_pin_buf is double[nf*N3] */
            double *base = (double *)ctx->h_pin_buf;
            for (int c = 0; c < nf; c++)
                cols[c] = base + c * N3;
        }
        fw_write_voxel(ctx->fw, ctx->frame_time, cols);

        pthread_mutex_lock(&ctx->mutex);
        ctx->writer_busy = 0;
        pthread_cond_signal(&ctx->cond);
        pthread_mutex_unlock(&ctx->mutex);
    }
    return NULL;
}

static SnapHookCtx create_snap_hook(FrameWriter *fw, int precision, int snap_every, long N3, int nf) {
    SnapHookCtx ctx;
    memset(&ctx, 0, sizeof(ctx));
    ctx.fw = fw;
    ctx.precision = precision;
    ctx.snap_every = snap_every;
    ctx.N3 = N3;
    ctx.nf = nf;   /* 12 real-mode columns, 24 complex (v66), 30 gauged (v69) */

    cudaStreamCreate(&ctx.stream);

    /* Allocate GPU staging buffer and pinned host buffer */
    if (precision == 0) {
        ctx.buf_bytes = (long)nf * N3 * sizeof(uint16_t);
        cudaMalloc(&ctx.d_f16_buf, ctx.buf_bytes);
        ctx.d_f32_buf = NULL;
    } else if (precision == 1) {
        ctx.buf_bytes = (long)nf * N3 * sizeof(float);
        ctx.d_f16_buf = NULL;
        cudaMalloc(&ctx.d_f32_buf, ctx.buf_bytes);
    } else {
        /* f64: DMA directly from physics arrays, staging = pinned host only */
        ctx.buf_bytes = (long)nf * N3 * sizeof(double);
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
    int nf = ctx->nf;

    /* Order MUST match column registration: 12 real arrays, then the
       imaginary copy (phi_im, theta_im, phi_im_vel, theta_im_vel) — v66,
       then the v69 gauge sector (th links, E). */
    const double *src[30] = {
        state->phi[0], state->phi[1], state->phi[2],
        state->theta[0], state->theta[1], state->theta[2],
        state->vel_phi[0], state->vel_phi[1], state->vel_phi[2],
        state->vel_theta[0], state->vel_theta[1], state->vel_theta[2],
        state->phi_im[0], state->phi_im[1], state->phi_im[2],
        state->theta_im[0], state->theta_im[1], state->theta_im[2],
        state->vel_phi_im[0], state->vel_phi_im[1], state->vel_phi_im[2],
        state->vel_theta_im[0], state->vel_theta_im[1], state->vel_theta_im[2],
        state->th[0], state->th[1], state->th[2],
        state->Efield[0], state->Efield[1], state->Efield[2]
    };

    if (ctx->precision == 0) {
        /* f16 path: convert on GPU, DMA f16 */
        for (int c = 0; c < nf; c++) {
            downcast_f64_to_f16_kernel<<<blocks, THREADS_PER_BLOCK, 0, ctx->stream>>>(
                src[c], ctx->d_f16_buf + c * N3, N3);
        }
        cudaStreamSynchronize(ctx->stream);
        cudaMemcpyAsync(ctx->h_pin_buf, ctx->d_f16_buf, ctx->buf_bytes,
                        cudaMemcpyDeviceToHost, ctx->stream);
    } else if (ctx->precision == 1) {
        /* f32 path: convert on GPU, DMA f32 */
        for (int c = 0; c < nf; c++) {
            downcast_f64_to_f32_kernel<<<blocks, THREADS_PER_BLOCK, 0, ctx->stream>>>(
                src[c], ctx->d_f32_buf + c * N3, N3);
        }
        cudaStreamSynchronize(ctx->stream);
        cudaMemcpyAsync(ctx->h_pin_buf, ctx->d_f32_buf, ctx->buf_bytes,
                        cudaMemcpyDeviceToHost, ctx->stream);
    } else {
        /* f64 path: DMA raw doubles (nf separate copies) */
        double *dst = (double *)ctx->h_pin_buf;
        for (int c = 0; c < nf; c++) {
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
    ctx.complex_mode = c->complex_phi;

    /* v69 gauge diagnostics config (G_offset set by the caller after
     * init_gauss_project — the host Grid owns the measured offset) */
    ctx.gauge_mode = c->complex_gauge;
    ctx.gauged = (c->complex_gauge && c->g_gauge != 0.0);
    ctx.g_gauge = c->g_gauge;
    ctx.G_offset = 0.0;
    /* Interior Omega: one extra dx inside the sponge edge (CPU compute_gauss) */
    double Rint = (c->bc_type == 0) ? (c->L - c->damp_width - dx) : 1e300;
    ctx.Rint2 = Rint * Rint;
    ctx.Rcube = c->qdiag_radius;
    ctx.gauss_max = ctx.gauss_l2 = ctx.e_em = ctx.q_flux = 0.0;

    int nvals = ctx.gauge_mode ? GDIAG_TOTAL : (ctx.complex_mode ? CDIAG_TOTAL : DIAG_NVALS);
    cudaMalloc(&ctx.d_results, nvals * sizeof(double));
    cudaMallocHost(&ctx.h_results, nvals * sizeof(double));

    /* v67: argmax-s index buffer + Q_core radius + theta probe voxel indices
     * (complex mode only; real mode allocates and launches nothing new) */
    if (ctx.complex_mode) {
        cudaMalloc(&ctx.d_sidx, sizeof(unsigned long long));
        ctx.qdiag_R2 = c->qdiag_radius * c->qdiag_radius;
        int N = c->N;
        double L = c->L;
        int jc = (int)lround(L/dx);                    if (jc<0) jc=0; if (jc>N-1) jc=N-1;
        int i1 = (int)lround((c->qdiag_probe1+L)/dx);  if (i1<0) i1=0; if (i1>N-1) i1=N-1;
        int i2 = (int)lround((c->qdiag_probe2+L)/dx);  if (i2<0) i2=0; if (i2>N-1) i2=N-1;
        ctx.probe1_idx = (long)i1*N*N + (long)jc*N + jc;
        ctx.probe2_idx = (long)i2*N*N + (long)jc*N + jc;
    }
    return ctx;
}

static void destroy_diag_hook(DiagHookCtx *ctx) {
    cudaFree(ctx->d_results);
    cudaFreeHost(ctx->h_results);
    if (ctx->d_sidx) cudaFree(ctx->d_sidx);
}

/* Run the GPU reduction and populate h_results. Used by diag_hook and for initial/final diag. */
static void run_gpu_diagnostics(const FieldState *state, DiagHookCtx *ctx) {
    double idx1 = 1.0 / (2.0 * ctx->dx);
    double idx2_val = 1.0 / (ctx->dx * ctx->dx);

    if (ctx->complex_mode) {
        /* v66 complex path: 20-value reduction + two-pass centroid r_core
         * (mirrors CPU compute_energy_complex + compute_charges).
         * v69 gauged: covariant-difference energy variant + gauss reduction. */
        int nvals = ctx->gauge_mode ? GDIAG_TOTAL : CDIAG_TOTAL;
        zero_diag_kernel<<<1, nvals>>>(ctx->d_results, nvals);
        if (ctx->gauged) {
            reduce_diagnostics_complex_gauge_kernel<<<gpu_blocks, THREADS_PER_BLOCK>>>(
                state->phi[0], state->phi[1], state->phi[2],
                state->phi_im[0], state->phi_im[1], state->phi_im[2],
                state->theta[0], state->theta[1], state->theta[2],
                state->theta_im[0], state->theta_im[1], state->theta_im[2],
                state->vel_phi[0], state->vel_phi[1], state->vel_phi[2],
                state->vel_phi_im[0], state->vel_phi_im[1], state->vel_phi_im[2],
                state->vel_theta[0], state->vel_theta[1], state->vel_theta[2],
                state->vel_theta_im[0], state->vel_theta_im[1], state->vel_theta_im[2],
                state->th[0], state->th[1], state->th[2],
                ctx->dV, idx1,
                ctx->m2, ctx->mtheta2, ctx->eta, ctx->mu, ctx->kappa,
                ctx->qdiag_R2, ctx->d_results);
        } else {
            reduce_diagnostics_complex_kernel<<<gpu_blocks, THREADS_PER_BLOCK>>>(
                state->phi[0], state->phi[1], state->phi[2],
                state->phi_im[0], state->phi_im[1], state->phi_im[2],
                state->theta[0], state->theta[1], state->theta[2],
                state->theta_im[0], state->theta_im[1], state->theta_im[2],
                state->vel_phi[0], state->vel_phi[1], state->vel_phi[2],
                state->vel_phi_im[0], state->vel_phi_im[1], state->vel_phi_im[2],
                state->vel_theta[0], state->vel_theta[1], state->vel_theta[2],
                state->vel_theta_im[0], state->vel_theta_im[1], state->vel_theta_im[2],
                ctx->dV, idx1,
                ctx->m2, ctx->mtheta2, ctx->eta, ctx->mu, ctx->kappa,
                ctx->qdiag_R2, ctx->d_results);
        }
        if (ctx->gauge_mode) {
            /* v69 gauss/Q_flux reduction (also when g==0: gauss cols still written) */
            reduce_gauss_kernel<<<gpu_blocks, THREADS_PER_BLOCK>>>(
                state->phi[0], state->phi[1], state->phi[2],
                state->phi_im[0], state->phi_im[1], state->phi_im[2],
                state->theta[0], state->theta[1], state->theta[2],
                state->theta_im[0], state->theta_im[1], state->theta_im[2],
                state->vel_phi[0], state->vel_phi[1], state->vel_phi[2],
                state->vel_phi_im[0], state->vel_phi_im[1], state->vel_phi_im[2],
                state->vel_theta[0], state->vel_theta[1], state->vel_theta[2],
                state->vel_theta_im[0], state->vel_theta_im[1], state->vel_theta_im[2],
                state->Efield[0], state->Efield[1], state->Efield[2],
                state->th[0], state->th[1], state->th[2],
                ctx->g_gauge, ctx->G_offset, ctx->Rint2, ctx->Rcube, ctx->dV,
                ctx->d_results + GAUSS_BASE);
        }
        cudaMemcpy(ctx->h_results, ctx->d_results, nvals * sizeof(double),
                   cudaMemcpyDeviceToHost);
        if (ctx->gauge_mode) {
            /* host post-processing (mirrors CPU compute_gauss outputs) */
            const double *h = ctx->h_results;
            double nOmega = h[GAUSS_BASE+2], nC = h[GAUSS_BASE+5];
            ctx->gauss_max = h[GAUSS_BASE+0];
            ctx->gauss_l2  = (nOmega > 0) ? sqrt(h[GAUSS_BASE+1]/nOmega) : 0.0;
            ctx->e_em      = h[GAUSS_BASE+3];
            ctx->q_flux    = (ctx->g_gauge != 0)
                           ? (h[GAUSS_BASE+4] - ctx->G_offset*nC*ctx->dV)/ctx->g_gauge : 0.0;
        }
        /* pass 2: rms charge radius about the rho2 centroid (dV cancels in centroid) */
        double M0_raw = ctx->h_results[15];
        if (M0_raw * ctx->dV >= 1e-30) {
            double cx = ctx->h_results[16] / M0_raw;
            double cy = ctx->h_results[17] / M0_raw;
            double cz = ctx->h_results[18] / M0_raw;
            reduce_rcore_kernel<<<gpu_blocks, THREADS_PER_BLOCK>>>(
                state->phi[0], state->phi[1], state->phi[2],
                state->phi_im[0], state->phi_im[1], state->phi_im[2],
                cx, cy, cz, ctx->d_results + 20);
            cudaMemcpy(&ctx->h_results[20], ctx->d_results + 20, sizeof(double),
                       cudaMemcpyDeviceToHost);
        } else {
            ctx->h_results[20] = 0.0;
        }

        /* v67 pass 3: omega_core at the argmax-s voxel (tiny kernel reading the
         * s_max value found by pass 1, then 4 single-voxel DtoH copies) */
        double smax = ctx->h_results[14];
        cudaMemset(ctx->d_sidx, 0xFF, sizeof(unsigned long long));
        find_smax_idx_kernel<<<gpu_blocks, THREADS_PER_BLOCK>>>(
            state->phi[0], state->phi[1], state->phi[2],
            state->phi_im[0], state->phi_im[1], state->phi_im[2],
            smax, ctx->d_sidx);
        unsigned long long h_sidx = ~0ULL;
        cudaMemcpy(&h_sidx, ctx->d_sidx, sizeof(unsigned long long),
                   cudaMemcpyDeviceToHost);
        ctx->omega_core = 0.0;
        if (h_sidx != ~0ULL && (long)h_sidx < ctx->N3) {
            long si = (long)h_sidx;
            double u0v, v0v, ud0, vd0;
            cudaMemcpy(&u0v, state->phi[0]        + si, sizeof(double), cudaMemcpyDeviceToHost);
            cudaMemcpy(&v0v, state->phi_im[0]     + si, sizeof(double), cudaMemcpyDeviceToHost);
            cudaMemcpy(&ud0, state->vel_phi[0]    + si, sizeof(double), cudaMemcpyDeviceToHost);
            cudaMemcpy(&vd0, state->vel_phi_im[0] + si, sizeof(double), cudaMemcpyDeviceToHost);
            double den = u0v*u0v + v0v*v0v;
            if (den >= 1e-20) ctx->omega_core = (u0v*vd0 - v0v*ud0)/den;
        }

        /* v67 theta probes: single-voxel point samples for offline FFT.
         * Component 1 (not 0): for a centered symmetric ball the theta source
         * eta*curl(phi) has component 0 ~ (y-z), identically zero on the y=z
         * probe diagonal; component 1 ~ (z-x) is maximal on the +x axis. */
        cudaMemcpy(&ctx->probes[0], state->theta[1]    + ctx->probe1_idx, sizeof(double), cudaMemcpyDeviceToHost);
        cudaMemcpy(&ctx->probes[1], state->theta_im[1] + ctx->probe1_idx, sizeof(double), cudaMemcpyDeviceToHost);
        cudaMemcpy(&ctx->probes[2], state->theta[1]    + ctx->probe2_idx, sizeof(double), cudaMemcpyDeviceToHost);
        cudaMemcpy(&ctx->probes[3], state->theta_im[1] + ctx->probe2_idx, sizeof(double), cudaMemcpyDeviceToHost);
        return;
    }

    /* Zero the results buffer */
    zero_diag_kernel<<<1, DIAG_NVALS>>>(ctx->d_results, DIAG_NVALS);

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

/* v66: extract Noether charges + core diagnostics from a completed complex
 * reduction (raw sums -> physical values, mirroring CPU compute_charges). */
static void complex_charges_from_results(const DiagHookCtx *ctx,
    double *Qp, double *Qt, double *smax, double *rcore,
    double *omega_core, double *Q_core) {
    const double *r = ctx->h_results;
    *Qp = r[12] * ctx->dV;
    *Qt = r[13] * ctx->dV;
    *smax = r[14];
    double M0 = r[15] * ctx->dV;
    *rcore = (M0 < 1e-30) ? 0.0 : sqrt(r[20] * ctx->dV / M0);
    *omega_core = ctx->omega_core;   /* computed in run_gpu_diagnostics pass 3 */
    *Q_core = r[19] * ctx->dV;
}

static void diag_hook(int step, double t, const FieldState *state, void *vctx) {
    DiagHookCtx *ctx = (DiagHookCtx *)vctx;
    if (step % ctx->diag_every != 0) return;

    run_gpu_diagnostics(state, ctx);

    double *r = ctx->h_results;
    /* r: [0]epk [1]etk [2]eg [3]em [4]ep [5]etg [6]etm [7]ec [8]pm [9]Pm [10]Pint [11]trms_sum */
    /* v69: gauged E_total includes E_em (CPU compute_energy_complex_gauge) */
    double et = r[0]+r[1]+r[2]+r[3]+r[4]+r[5]+r[6]+r[7]
              + (ctx->gauged ? ctx->e_em : 0.0);
    /* theta_rms: complex mode divides by the actual 6 theta components (v66) */
    double trms = sqrt(r[11] / ((ctx->complex_mode ? 6.0 : 3.0) * ctx->N3));
    double Qp = 0, Qt = 0, smax = 0, rcore = 0, ocore = 0, qcore = 0;
    if (ctx->complex_mode)
        complex_charges_from_results(ctx, &Qp, &Qt, &smax, &rcore, &ocore, &qcore);

    fprintf(ctx->fp, "%.2f\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e",
            t, r[0], r[1], r[2], r[3], r[4], r[5], r[6], r[7], et, r[8], r[9], r[10], trms);
    if (ctx->complex_mode)
        fprintf(ctx->fp, "\t%.12e\t%.12e\t%.12e\t%.6e\t%.6e\t%.12e\t%.12e\t%.6e\t%.6e\t%.6e\t%.6e",
                Qp, Qt, Qp+Qt, smax, rcore, ocore, qcore,
                ctx->probes[0], ctx->probes[1], ctx->probes[2], ctx->probes[3]);
    if (ctx->gauge_mode)
        fprintf(ctx->fp, "\t%.6e\t%.6e\t%.12e\t%.12e",
                ctx->gauss_max, ctx->gauss_l2, ctx->e_em, ctx->q_flux);
    fprintf(ctx->fp, "\n");
    fflush(ctx->fp);

    if (step % ctx->major_every == 0) {
        cudaEvent_t t_stop;
        cudaEventCreate(&t_stop);
        cudaEventRecord(t_stop);
        cudaEventSynchronize(t_stop);
        float ms; cudaEventElapsedTime(&ms, ctx->t_start, t_stop);
        cudaEventDestroy(t_stop);
        double drift = 100.0 * (et - ctx->E0) / (fabs(ctx->E0) + 1e-30);
        printf("t=%7.1f E=%.3e (drift %+.3f%%) Ep=%.1f phi=%.3f θ_rms=%.2e ",
               t, et, drift, r[4], r[8], trms);
        if (ctx->complex_mode) printf("Q=%.6e s_max=%.4e ", Qp+Qt, smax);
        if (ctx->gauge_mode) printf("gauss=%.1e E_em=%.3e ", ctx->gauss_max, ctx->e_em);
        printf("[%.0f%% %.1fs %.2fms/step]\n",
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
   Vecstream output hook
   ================================================================ */

/* Vector frame output uses SFA's native FRVD chunks — no separate vecstream format */
#include "vecstream_gpu_v2.cuh"

typedef struct {
    /* Config */
    float delta_tol;
    int iframe_interval;
    int block_size;
    int refit_interval;   /* refit temporal model every N frames */
    float omega;          /* breathing frequency for temporal model */
    /* State */
    FrameWriter  *fw;    /* shared frame writer (same file as voxel frames) */
    int           enabled;
    int           N;
    long          N3;
    int           BN, n_patches;
    int           vec_frame;
    int           vec_snap_every;
    int           cpp;            /* coefficients per patch (VS2_MULTI_TOTAL) */
    long          n_total;        /* n_patches × cpp */
    /* GPU buffers */
    float        *d_coeffs;       /* n_patches × cpp (current frame) */
    float        *d_prev_coeffs;  /* n_patches × cpp (previous frame for fallback) */
    float        *d_predicted;    /* n_patches × cpp (temporal prediction) */
    float        *d_residual;     /* n_patches × cpp (actual - predicted) */
    uint32_t     *d_nz_patches;   /* nonzero patch indices */
    uint32_t     *d_n_nonzero;    /* atomic counter */
    /* Temporal model on GPU */
    float        *d_temp_mean;    /* n_total */
    float        *d_temp_amp;     /* n_total */
    float        *d_temp_phase;   /* n_total */
    float        *d_sum_cos;      /* n_total accumulators */
    float        *d_sum_sin;
    float        *d_sum_mean;
    int           temporal_count; /* frames since last refit */
    /* Host pinned buffers */
    float        *h_coeffs;       /* n_total pinned */
    float        *h_residual;     /* n_total pinned */
    uint32_t     *h_nz_patches;   /* n_patches pinned */
    float        *h_gather_buf;   /* n_total pinned — pre-allocated gather buffer */
    cudaStream_t  stream;
} VecHookCtx;

static VecHookCtx create_vec_hook(const Config *c, double dt, int N3, FrameWriter *fw) {
    VecHookCtx vh = {};
    if (c->vec_snap_dt <= 0) { vh.enabled = 0; return vh; }

    vh.enabled = 1;
    vh.N = c->N;
    vh.N3 = N3;
    int BS = c->vec_block_size > 0 ? c->vec_block_size : 8;
    vh.BN = c->N / BS;
    vh.n_patches = vh.BN * vh.BN * vh.BN;
    vh.delta_tol = c->vec_delta_tol > 0 ? c->vec_delta_tol : 0.01f;
    vh.iframe_interval = c->vec_iframe_interval > 0 ? c->vec_iframe_interval : 100;
    vh.block_size = BS;
    vh.cpp = VS2_MULTI_TOTAL;  /* 6 × 64 = 384 coefficients per patch */
    vh.n_total = (long)vh.n_patches * vh.cpp;
    vh.omega = 2.0 * 3.14159265358979 / 2.2;  /* breathing frequency */
    vh.refit_interval = 50;  /* refit temporal model every 50 vec frames */

    double snap_dt = c->vec_snap_dt > 0 ? c->vec_snap_dt : c->snap_dt;
    vh.vec_snap_every = (int)lround(snap_dt / dt);
    if (vh.vec_snap_every < 1) vh.vec_snap_every = 1;

    /* Use the shared FrameWriter — vector frames interleave with voxel frames in one file */
    vh.fw = fw;
    if (!vh.fw) {
        fprintf(stderr, "vecstream: no frame writer provided\n");
        vh.enabled = 0;
        return vh;
    }

    /* GPU buffers — check each allocation */
    long total_bytes = vh.n_total * sizeof(float);
    cudaError_t err;
    size_t free_mem, total_mem;
    cudaMemGetInfo(&free_mem, &total_mem);
    fprintf(stderr, "Vecstream: need %.0f MB GPU, %.0f MB free / %.0f MB total\n",
            total_bytes * 3.0 / 1e6, free_mem / 1e6, total_mem / 1e6);

    err = cudaMalloc(&vh.d_coeffs, total_bytes);       if (err) { fprintf(stderr, "vec cudaMalloc: %s\n", cudaGetErrorString(err)); vh.enabled=0; return vh; }
    err = cudaMalloc(&vh.d_prev_coeffs, total_bytes);   if (err) { vh.enabled=0; return vh; }
    err = cudaMalloc(&vh.d_predicted, total_bytes);     if (err) { vh.enabled=0; return vh; }
    err = cudaMalloc(&vh.d_residual, total_bytes);      if (err) { vh.enabled=0; return vh; }
    err = cudaMalloc(&vh.d_nz_patches, vh.n_patches * sizeof(uint32_t)); if (err) { vh.enabled=0; return vh; }
    err = cudaMalloc(&vh.d_n_nonzero, sizeof(uint32_t)); if (err) { vh.enabled=0; return vh; }
    err = cudaMalloc(&vh.d_temp_mean, total_bytes);     if (err) { vh.enabled=0; return vh; }
    err = cudaMalloc(&vh.d_temp_amp, total_bytes);      if (err) { vh.enabled=0; return vh; }
    err = cudaMalloc(&vh.d_temp_phase, total_bytes);    if (err) { vh.enabled=0; return vh; }
    err = cudaMalloc(&vh.d_sum_cos, total_bytes);       if (err) { vh.enabled=0; return vh; }
    err = cudaMalloc(&vh.d_sum_sin, total_bytes);       if (err) { vh.enabled=0; return vh; }
    err = cudaMalloc(&vh.d_sum_mean, total_bytes);      if (err) { vh.enabled=0; return vh; }
    cudaMemset(vh.d_prev_coeffs, 0, total_bytes);
    cudaMemset(vh.d_temp_mean, 0, total_bytes);
    cudaMemset(vh.d_temp_amp, 0, total_bytes);
    cudaMemset(vh.d_temp_phase, 0, total_bytes);
    cudaMemset(vh.d_sum_cos, 0, total_bytes);
    cudaMemset(vh.d_sum_sin, 0, total_bytes);
    cudaMemset(vh.d_sum_mean, 0, total_bytes);
    vh.temporal_count = 0;

    cudaMallocHost(&vh.h_coeffs, total_bytes);
    cudaMallocHost(&vh.h_residual, total_bytes);
    cudaMallocHost(&vh.h_nz_patches, vh.n_patches * sizeof(uint32_t));
    cudaMallocHost(&vh.h_gather_buf, total_bytes);
    cudaStreamCreate(&vh.stream);

    vs2_gpu_init_projection();

    vh.vec_frame = 0;
    printf("Vec frames: BS=%d, %d patches, %d coeffs/patch, snap every %d steps\n",
           BS, vh.n_patches, vh.cpp, vh.vec_snap_every);
    return vh;
}

/* Vector frame hook: fit patches, temporal model, I/P frames */
static void vec_hook(int step, double t, const FieldState *state, void *vctx) {
    VecHookCtx *vh = (VecHookCtx*)vctx;
    if (!vh->enabled) return;
    if (step % vh->vec_snap_every != 0) return;

    int np = vh->n_patches;
    int BN = vh->BN;
    long n_total = vh->n_total;
    int blk256 = (int)((n_total + 255) / 256);

    /* Ensure physics kernels on default stream are complete before reading fields */
    cudaDeviceSynchronize();

    /* Fit all 6 fields into d_coeffs */
    vs2_fit_multi_patches<<<np, 512, 0, vh->stream>>>(
        state->phi[0], state->phi[1], state->phi[2],
        state->theta[0], state->theta[1], state->theta[2],
        vh->N, BN, vh->d_coeffs);

    /* Update temporal accumulators */
    float cos_wt = cosf(vh->omega * (float)t);
    float sin_wt = sinf(vh->omega * (float)t);
    vs2_temporal_update<<<blk256, 256, 0, vh->stream>>>(
        vh->d_coeffs, vh->d_sum_cos, vh->d_sum_sin, vh->d_sum_mean,
        cos_wt, sin_wt, (int)n_total);
    vh->temporal_count++;

    int is_iframe = (vh->vec_frame == 0) || (vh->vec_frame % vh->iframe_interval == 0);

    /* Bootstrap: on first frame, set mean = actual coefficients */
    if (vh->vec_frame == 0) {
        cudaMemcpyAsync(vh->d_temp_mean, vh->d_coeffs, n_total * sizeof(float),
                        cudaMemcpyDeviceToDevice, vh->stream);
        cudaMemsetAsync(vh->d_temp_amp, 0, n_total * sizeof(float), vh->stream);
        cudaMemsetAsync(vh->d_temp_phase, 0, n_total * sizeof(float), vh->stream);
    }

    /* Refit at every I-frame boundary or periodically */
    if ((is_iframe || vh->temporal_count >= vh->refit_interval) && vh->temporal_count > 2) {
        vs2_temporal_refit<<<blk256, 256, 0, vh->stream>>>(
            vh->d_sum_cos, vh->d_sum_sin, vh->d_sum_mean,
            vh->temporal_count,
            vh->d_temp_mean, vh->d_temp_amp, vh->d_temp_phase, (int)n_total);
        cudaMemsetAsync(vh->d_sum_cos, 0, n_total * sizeof(float), vh->stream);
        cudaMemsetAsync(vh->d_sum_sin, 0, n_total * sizeof(float), vh->stream);
        cudaMemsetAsync(vh->d_sum_mean, 0, n_total * sizeof(float), vh->stream);
        vh->temporal_count = 0;
    }

    if (is_iframe) {
        /* Download coefficients + temporal model */
        cudaMemcpyAsync(vh->h_coeffs, vh->d_coeffs, n_total * sizeof(float),
                        cudaMemcpyDeviceToHost, vh->stream);
        cudaStreamSynchronize(vh->stream);

        int16_t *origins = (int16_t*)vh->h_gather_buf;
        for (int p = 0; p < np; p++) {
            int bk = p % BN, bj = (p / BN) % BN, bi = p / (BN * BN);
            origins[p*3 + 0] = (int16_t)(bi * vh->block_size);
            origins[p*3 + 1] = (int16_t)(bj * vh->block_size);
            origins[p*3 + 2] = (int16_t)(bk * vh->block_size);
        }

        /* Download temporal model for embedding */
        float *h_mean = (float*)malloc(n_total * sizeof(float));
        float *h_amp  = (float*)malloc(n_total * sizeof(float));
        float *h_phase = (float*)malloc(n_total * sizeof(float));
        cudaMemcpy(h_mean, vh->d_temp_mean, n_total * sizeof(float), cudaMemcpyDeviceToHost);
        cudaMemcpy(h_amp, vh->d_temp_amp, n_total * sizeof(float), cudaMemcpyDeviceToHost);
        cudaMemcpy(h_phase, vh->d_temp_phase, n_total * sizeof(float), cudaMemcpyDeviceToHost);

        /* Write I-frame with temporal model via FrameWriter */
        {
            uint8_t flags = 0x01;
            uint64_t origin_bytes = (uint64_t)np * 6;
            uint64_t coeff_bytes = n_total * 4;
            uint64_t model_bytes = 4 + n_total * 4 * 3;
            uint64_t payload_size = 8 + origin_bytes + coeff_bytes + model_bytes;
            uint8_t *pl = (uint8_t*)malloc(payload_size);
            uint64_t off = 0;
            uint32_t np32 = np; uint8_t bs8 = vh->block_size; uint16_t nc16 = vh->cpp;
            memcpy(pl+off, &np32, 4); off += 4;
            memcpy(pl+off, &bs8, 1); off += 1;
            memcpy(pl+off, &nc16, 2); off += 2;
            memcpy(pl+off, &flags, 1); off += 1;
            memcpy(pl+off, origins, origin_bytes); off += origin_bytes;
            memcpy(pl+off, vh->h_coeffs, coeff_bytes); off += coeff_bytes;
            memcpy(pl+off, &vh->omega, 4); off += 4;
            memcpy(pl+off, h_mean, n_total * 4); off += n_total * 4;
            memcpy(pl+off, h_amp, n_total * 4); off += n_total * 4;
            memcpy(pl+off, h_phase, n_total * 4); off += n_total * 4;

            size_t bound = ZSTD_compressBound(payload_size);
            void *comp = malloc(bound);
            size_t comp_size = ZSTD_compress(comp, bound, pl, payload_size, 3);
            uint32_t checksum = sfa_crc32(pl, payload_size);
            fw_write_raw(vh->fw, SFA_CHUNK_FRVD, SFA_FRAME_VEC_I, t, comp, comp_size, checksum);
            free(comp); free(pl);
        }
        free(h_mean); free(h_amp); free(h_phase);

        cudaMemcpyAsync(vh->d_prev_coeffs, vh->d_coeffs, n_total * sizeof(float),
                        cudaMemcpyDeviceToDevice, vh->stream);
    } else {
        /* P-frame: delta = actual - predicted(t) */
        vs2_temporal_predict<<<blk256, 256, 0, vh->stream>>>(
            vh->d_temp_mean, vh->d_temp_amp, vh->d_temp_phase,
            vh->omega, (float)t, vh->d_predicted, (int)n_total);

        uint32_t zero = 0;
        cudaMemcpyAsync(vh->d_n_nonzero, &zero, sizeof(uint32_t),
                        cudaMemcpyHostToDevice, vh->stream);
        vs2_temporal_residual<<<(np+255)/256, 256, 0, vh->stream>>>(
            vh->d_coeffs, vh->d_predicted, vh->d_residual,
            vh->delta_tol, vh->d_nz_patches, vh->d_n_nonzero,
            np, vh->cpp);
        uint32_t h_nnz;
        cudaMemcpyAsync(&h_nnz, vh->d_n_nonzero, sizeof(uint32_t),
                        cudaMemcpyDeviceToHost, vh->stream);
        cudaStreamSynchronize(vh->stream);

        if (h_nnz > 0 && h_nnz <= (uint32_t)np) {
            cudaMemcpyAsync(vh->h_nz_patches, vh->d_nz_patches,
                            h_nnz * sizeof(uint32_t), cudaMemcpyDeviceToHost, vh->stream);
            cudaMemcpyAsync(vh->h_residual, vh->d_residual, n_total * sizeof(float),
                            cudaMemcpyDeviceToHost, vh->stream);
            cudaStreamSynchronize(vh->stream);

            for (uint32_t i = 0; i < h_nnz; i++)
                memcpy(vh->h_gather_buf + (long)i * vh->cpp,
                       vh->h_residual + (long)vh->h_nz_patches[i] * vh->cpp,
                       vh->cpp * sizeof(float));

            fw_write_vec_pframe(vh->fw, t, h_nnz, vh->cpp,
                                vh->h_nz_patches, vh->h_gather_buf);
        } else {
            fw_write_vec_pframe(vh->fw, t, 0, vh->cpp, NULL, NULL);
        }
        cudaMemcpyAsync(vh->d_prev_coeffs, vh->d_coeffs, n_total * sizeof(float),
                        cudaMemcpyDeviceToDevice, vh->stream);
    }
    vh->vec_frame++;
}

static void destroy_vec_hook(VecHookCtx *vh) {
    if (!vh->enabled || !vh->fw) return;
    vh->fw = NULL;
    cudaFree(vh->d_coeffs); cudaFree(vh->d_prev_coeffs);
    cudaFree(vh->d_predicted); cudaFree(vh->d_residual);
    cudaFree(vh->d_nz_patches); cudaFree(vh->d_n_nonzero);
    cudaFree(vh->d_temp_mean); cudaFree(vh->d_temp_amp); cudaFree(vh->d_temp_phase);
    cudaFree(vh->d_sum_cos); cudaFree(vh->d_sum_sin); cudaFree(vh->d_sum_mean);
    cudaFreeHost(vh->h_coeffs); cudaFreeHost(vh->h_residual);
    cudaFreeHost(vh->h_nz_patches); cudaFreeHost(vh->h_gather_buf);
    cudaStreamDestroy(vh->stream);
    fprintf(stderr, "Vecstream v2: %d frames written\n", vh->vec_frame);
}

/* ================================================================
   Inline cluster analysis + parameter tuning
   ================================================================ */

static void inline_cluster_analysis(Config *c, int N, long N3, double L, double dx,
                                     double t, int step) {
    /* Download phi fields from GPU to temporary host buffer */
    double *h_phi[3];
    for (int a = 0; a < 3; a++) {
        h_phi[a] = (double*)malloc(N3 * sizeof(double));
        cudaMemcpy(h_phi[a], d_phi[a], N3*sizeof(double), cudaMemcpyDeviceToHost);
    }

    /* Compute |phi|^2 and P = phi0*phi1*phi2 at each voxel */
    double phi2_max = 0, phi2_mean = 0, P_max = 0, P_mean = 0;
    int NN = N*N;
    for (long idx = 0; idx < N3; idx++) {
        double p2 = h_phi[0][idx]*h_phi[0][idx] + h_phi[1][idx]*h_phi[1][idx]
                   + h_phi[2][idx]*h_phi[2][idx];
        double P = fabs(h_phi[0][idx] * h_phi[1][idx] * h_phi[2][idx]);
        if (p2 > phi2_max) phi2_max = p2;
        phi2_mean += p2;
        if (P > P_max) P_max = P;
        P_mean += P;
    }
    phi2_mean /= N3;
    P_mean /= N3;

    /* Flood-fill cluster detection on |phi|^2 > threshold */
    /* Threshold: 3x the mean — anything above is "concentrated" */
    double threshold = 3.0 * phi2_mean;
    if (threshold < 0.01) threshold = 0.01;

    int *label = (int*)calloc(N3, sizeof(int));  /* 0 = unlabeled */
    int *stack = (int*)malloc(N3 * sizeof(int));
    int n_clusters = 0;
    double best_mass = 0, best_rms = 0, best_cx = 0, best_cy = 0, best_cz = 0;
    int best_nvox = 0;

    for (long seed = 0; seed < N3; seed++) {
        double p2 = h_phi[0][seed]*h_phi[0][seed] + h_phi[1][seed]*h_phi[1][seed]
                   + h_phi[2][seed]*h_phi[2][seed];
        if (p2 < threshold || label[seed] != 0) continue;

        /* New cluster */
        n_clusters++;
        int cid = n_clusters;
        int sp = 0;
        stack[sp++] = (int)seed;
        label[seed] = cid;

        double mass = 0, cx = 0, cy = 0, cz = 0;
        int nvox = 0;

        while (sp > 0) {
            int idx = stack[--sp];
            int i = idx / NN, j = (idx / N) % N, k = idx % N;
            double p2v = h_phi[0][idx]*h_phi[0][idx] + h_phi[1][idx]*h_phi[1][idx]
                       + h_phi[2][idx]*h_phi[2][idx];
            double w = p2v * dx*dx*dx;
            mass += w;
            cx += (-L + i*dx) * w;
            cy += (-L + j*dx) * w;
            cz += (-L + k*dx) * w;
            nvox++;

            /* 6-connected neighbors (periodic) */
            int nbrs[6] = {
                ((i+1)%N)*NN+j*N+k, ((i-1+N)%N)*NN+j*N+k,
                i*NN+((j+1)%N)*N+k, i*NN+((j-1+N)%N)*N+k,
                i*NN+j*N+((k+1)%N), i*NN+j*N+((k-1+N)%N)
            };
            for (int n = 0; n < 6; n++) {
                int ni = nbrs[n];
                if (label[ni] != 0) continue;
                double np2 = h_phi[0][ni]*h_phi[0][ni] + h_phi[1][ni]*h_phi[1][ni]
                            + h_phi[2][ni]*h_phi[2][ni];
                if (np2 >= threshold) {
                    label[ni] = cid;
                    stack[sp++] = ni;
                }
            }
        }

        if (mass > 0) { cx /= mass; cy /= mass; cz /= mass; }

        /* Compute rms_r for this cluster */
        double r2_sum = 0;
        for (long idx2 = 0; idx2 < N3; idx2++) {
            if (label[idx2] != cid) continue;
            int i2 = (int)(idx2 / NN), j2 = (int)((idx2 / N) % N), k2 = (int)(idx2 % N);
            double xi = -L+i2*dx-cx, yi = -L+j2*dx-cy, zi = -L+k2*dx-cz;
            double p2v = h_phi[0][idx2]*h_phi[0][idx2]+h_phi[1][idx2]*h_phi[1][idx2]
                        +h_phi[2][idx2]*h_phi[2][idx2];
            r2_sum += (xi*xi+yi*yi+zi*zi) * p2v * dx*dx*dx;
        }
        double rms_r = (mass > 0) ? sqrt(r2_sum / mass) : L;

        if (mass > best_mass) {
            best_mass = mass; best_rms = rms_r; best_nvox = nvox;
            best_cx = cx; best_cy = cy; best_cz = cz;
        }
    }

    /* Metrics */
    double box_rms = L * sqrt(3.0/5.0);
    double compactness = (best_rms < box_rms) ? 1.0 - best_rms/box_rms : 0.0;
    double density_contrast = (phi2_mean > 0) ? phi2_max / phi2_mean : 0;
    double vox_frac = (double)best_nvox / N3;
    int has_particle = (best_nvox > 0 && best_nvox < N3/2 && best_rms < L*0.7);

    printf("\n[TUNE t=%.0f] clusters=%d top_mass=%.1f top_rms=%.2f top_nvox=%d(%.1f%%) "
           "compact=%.3f contrast=%.1f particle=%s\n",
           t, n_clusters, best_mass, best_rms, best_nvox, vox_frac*100,
           compactness, density_contrast, has_particle ? "YES" : "NO");

    /* === Parameter tuning === */
    int changed = 0;
    static int tune_count = 0;
    tune_count++;

    if (!has_particle) {
        /* No particle detected — nudge toward binding */
        if (fabs(c->mu) < 300) {
            c->mu -= 1.0;  /* more negative = stronger binding */
            printf("  TUNE: mu -> %.2f (no particle)\n", c->mu);
            changed = 1;
        }
    } else {
        /* Have particle — enlarge and harden */

        /* Tighten: if compactness < 0.5, increase m2 */
        if (compactness < 0.5) {
            c->m2 += 0.1;
            printf("  TUNE: m2 -> %.2f (compact=%.3f, tighten)\n", c->m2, compactness);
            changed = 1;
        }

        /* Harden: if density contrast < 50, strengthen binding */
        if (density_contrast < 50) {
            c->mu -= 0.5;
            printf("  TUNE: mu -> %.2f (contrast=%.1f, harden)\n", c->mu, density_contrast);
            changed = 1;
        }

        /* Enlarge: if top cluster rms_r < 3, boost eta to spread structure */
        if (best_rms < 3.0) {
            c->eta += 0.01;
            printf("  TUNE: eta -> %.3f (rms=%.2f, enlarge)\n", c->eta, best_rms);
            changed = 1;
        }

        /* If too many clusters, strengthen kappa to merge */
        if (n_clusters > 10) {
            c->kappa += 1.0;
            printf("  TUNE: kappa -> %.1f (%d clusters, merge)\n", c->kappa, n_clusters);
            changed = 1;
        }

        /* Report particle health */
        printf("  PARTICLE: mass=%.1f rms=%.2f nvox=%d contrast=%.0f compact=%.3f\n",
               best_mass, best_rms, best_nvox, density_contrast, compactness);
    }

    if (changed) {
        gpu_set_constants(c, dx);
        printf("  -> constants uploaded\n");
    }

    fflush(stdout);

    free(label); free(stack);
    for (int a = 0; a < 3; a++) free(h_phi[a]);
}

/* ================================================================
   Parameter sweep — runs entirely on GPU
   ================================================================ */

typedef struct {
    double mu, eta, m2, kappa, sigma_cross, lambda_self;
    double score;         /* higher = better particle */
    double phi_max, P_max, compactness, density_contrast;
    int n_clusters, best_nvox;
    double best_rms, best_mass;
} SweepResult;

static void run_sweep(Config *c, Grid *g) {
    int N = g->N;
    long N3 = g->N3;
    double dx = g->dx, dt = g->dt;

    /* Save initial field state */
    size_t field_bytes = N3 * sizeof(double);
    double *save_phi[3], *save_vel[3], *save_theta[3], *save_tvel[3];
    for (int a = 0; a < 3; a++) {
        save_phi[a] = (double*)malloc(field_bytes);
        save_vel[a] = (double*)malloc(field_bytes);
        save_theta[a] = (double*)malloc(field_bytes);
        save_tvel[a] = (double*)malloc(field_bytes);
        cudaMemcpy(save_phi[a], d_phi[a], field_bytes, cudaMemcpyDeviceToHost);
        cudaMemcpy(save_vel[a], d_vel_phi[a], field_bytes, cudaMemcpyDeviceToHost);
        cudaMemcpy(save_theta[a], d_theta[a], field_bytes, cudaMemcpyDeviceToHost);
        cudaMemcpy(save_tvel[a], d_vel_theta[a], field_bytes, cudaMemcpyDeviceToHost);
    }

    /* Parameter grid */
    double mu_vals[]    = {-20, -30, -41.345, -60, -80, -100};
    double eta_vals[]   = {0.2, 0.35, 0.5, 0.7, 1.0};
    double m2_vals[]    = {1.5, 2.25, 4.0, 6.25, 9.0};
    double kappa_vals[] = {20, 35, 50, 75, 100};
    int n_mu = sizeof(mu_vals)/sizeof(double);
    int n_eta = sizeof(eta_vals)/sizeof(double);
    int n_m2 = sizeof(m2_vals)/sizeof(double);
    int n_kappa = sizeof(kappa_vals)/sizeof(double);
    int total = n_mu * n_eta * n_m2 * n_kappa;

    printf("\n=== PARAMETER SWEEP: %d combinations, T=%.0f each ===\n", total, c->sweep_T);
    printf("  mu: "); for(int i=0;i<n_mu;i++) printf("%.1f ",mu_vals[i]); printf("\n");
    printf("  eta: "); for(int i=0;i<n_eta;i++) printf("%.2f ",eta_vals[i]); printf("\n");
    printf("  m2: "); for(int i=0;i<n_m2;i++) printf("%.2f ",m2_vals[i]); printf("\n");
    printf("  kappa: "); for(int i=0;i<n_kappa;i++) printf("%.0f ",kappa_vals[i]); printf("\n\n");
    fflush(stdout);

    int sweep_steps = (int)lround(c->sweep_T / dt);

    SweepResult *results = (SweepResult*)calloc(total, sizeof(SweepResult));
    int trial = 0;
    double best_score = -1;
    int best_trial = -1;

    for (int im = 0; im < n_mu; im++)
    for (int ie = 0; ie < n_eta; ie++)
    for (int i2 = 0; i2 < n_m2; i2++)
    for (int ik = 0; ik < n_kappa; ik++) {
        /* Set parameters */
        c->mu = mu_vals[im];
        c->eta = eta_vals[ie];
        c->m2 = m2_vals[i2];
        c->kappa = kappa_vals[ik];

        /* Restore initial field */
        for (int a = 0; a < 3; a++) {
            cudaMemcpy(d_phi[a], save_phi[a], field_bytes, cudaMemcpyHostToDevice);
            cudaMemcpy(d_vel_phi[a], save_vel[a], field_bytes, cudaMemcpyHostToDevice);
            cudaMemcpy(d_theta[a], save_theta[a], field_bytes, cudaMemcpyHostToDevice);
            cudaMemcpy(d_vel_theta[a], save_tvel[a], field_bytes, cudaMemcpyHostToDevice);
        }

        /* Upload new constants */
        gpu_set_constants(c, dx);

        /* Run for sweep_steps */
        /* Need initial force computation */
        if (gpu_has_intermediates) {
            compute_intermediates_kernel<<<gpu_blocks, THREADS_PER_BLOCK>>>(
                d_phi[0],d_phi[1],d_phi[2], d_theta[0],d_theta[1],d_theta[2],
                d_mismatch[0],d_mismatch[1],d_mismatch[2],
                d_harden_Q[0],d_harden_Q[1],d_harden_Q[2]);
        }
        compute_forces_kernel<<<gpu_blocks, THREADS_PER_BLOCK>>>(
            d_phi[0],d_phi[1],d_phi[2], d_theta[0],d_theta[1],d_theta[2],
            d_vel_theta[0],d_vel_theta[1],d_vel_theta[2],
            d_acc_phi[0],d_acc_phi[1],d_acc_phi[2],
            d_acc_theta[0],d_acc_theta[1],d_acc_theta[2],
            d_mismatch[0],d_mismatch[1],d_mismatch[2],
            d_harden_Q[0],d_harden_Q[1],d_harden_Q[2]);

        for (int s = 0; s < sweep_steps; s++)
            gpu_verlet_step(dt);

        cudaDeviceSynchronize();

        /* Analyze — download phi, flood fill */
        double *h_phi[3];
        for (int a = 0; a < 3; a++) {
            h_phi[a] = (double*)malloc(field_bytes);
            cudaMemcpy(h_phi[a], d_phi[a], field_bytes, cudaMemcpyDeviceToHost);
        }

        int NN = N*N;
        double phi2_max = 0, phi2_mean = 0, P_abs_max = 0;
        for (long idx = 0; idx < N3; idx++) {
            double p2 = h_phi[0][idx]*h_phi[0][idx]+h_phi[1][idx]*h_phi[1][idx]+h_phi[2][idx]*h_phi[2][idx];
            if (p2 > phi2_max) phi2_max = p2;
            phi2_mean += p2;
            double P = fabs(h_phi[0][idx]*h_phi[1][idx]*h_phi[2][idx]);
            if (P > P_abs_max) P_abs_max = P;
        }
        phi2_mean /= N3;

        double threshold = 3.0 * phi2_mean;
        if (threshold < 0.001) threshold = 0.001;

        int *label = (int*)calloc(N3, sizeof(int));
        int *stack = (int*)malloc(N3 * sizeof(int));
        int n_clusters = 0;
        double best_mass = 0, best_rms = 999;
        int best_nvox = 0;

        for (long seed = 0; seed < N3; seed++) {
            double p2 = h_phi[0][seed]*h_phi[0][seed]+h_phi[1][seed]*h_phi[1][seed]+h_phi[2][seed]*h_phi[2][seed];
            if (p2 < threshold || label[seed] != 0) continue;
            n_clusters++;
            int cid = n_clusters;
            int sp = 0;
            stack[sp++] = (int)seed;
            label[seed] = cid;
            double mass = 0, cx=0,cy=0,cz=0;
            int nvox = 0;
            while (sp > 0) {
                int idx = stack[--sp];
                int i=idx/NN, j=(idx/N)%N, k=idx%N;
                double p2v = h_phi[0][idx]*h_phi[0][idx]+h_phi[1][idx]*h_phi[1][idx]+h_phi[2][idx]*h_phi[2][idx];
                double w = p2v*dx*dx*dx;
                mass += w; cx+=(-g->L+i*dx)*w; cy+=(-g->L+j*dx)*w; cz+=(-g->L+k*dx)*w;
                nvox++;
                int nbrs[6]={((i+1)%N)*NN+j*N+k,((i-1+N)%N)*NN+j*N+k,
                             i*NN+((j+1)%N)*N+k,i*NN+((j-1+N)%N)*N+k,
                             i*NN+j*N+((k+1)%N),i*NN+j*N+((k-1+N)%N)};
                for(int n=0;n<6;n++){
                    int ni=nbrs[n];
                    if(label[ni]!=0) continue;
                    double np2=h_phi[0][ni]*h_phi[0][ni]+h_phi[1][ni]*h_phi[1][ni]+h_phi[2][ni]*h_phi[2][ni];
                    if(np2>=threshold){label[ni]=cid;stack[sp++]=ni;}
                }
            }
            if(mass>0){cx/=mass;cy/=mass;cz/=mass;}
            double r2s=0;
            for(long ii=0;ii<N3;ii++){
                if(label[ii]!=cid) continue;
                int i2=ii/NN,j2=(ii/N)%N,k2=ii%N;
                double xi=-g->L+i2*dx-cx,yi=-g->L+j2*dx-cy,zi=-g->L+k2*dx-cz;
                double p2v=h_phi[0][ii]*h_phi[0][ii]+h_phi[1][ii]*h_phi[1][ii]+h_phi[2][ii]*h_phi[2][ii];
                r2s+=(xi*xi+yi*yi+zi*zi)*p2v*dx*dx*dx;
            }
            double rms=(mass>0)?sqrt(r2s/mass):g->L;
            if(mass>best_mass){best_mass=mass;best_rms=rms;best_nvox=nvox;}
        }

        double box_rms = g->L*sqrt(3.0/5.0);
        double compactness = (best_rms<box_rms)?1.0-best_rms/box_rms:0.0;
        double contrast = (phi2_mean>0)?phi2_max/phi2_mean:0;
        int has_particle = (best_nvox>0 && best_nvox<N3/2 && best_rms<g->L*0.7);

        /* Score: compactness × contrast × P_max — rewards tight, dense, bound clusters */
        double score = compactness * contrast * sqrt(P_abs_max + 0.001);
        if (!has_particle) score = 0;

        results[trial].mu = c->mu;
        results[trial].eta = c->eta;
        results[trial].m2 = c->m2;
        results[trial].kappa = c->kappa;
        results[trial].score = score;
        results[trial].phi_max = sqrt(phi2_max);
        results[trial].P_max = P_abs_max;
        results[trial].compactness = compactness;
        results[trial].density_contrast = contrast;
        results[trial].n_clusters = n_clusters;
        results[trial].best_nvox = best_nvox;
        results[trial].best_rms = best_rms;
        results[trial].best_mass = best_mass;

        if (score > best_score) { best_score = score; best_trial = trial; }

        if (trial % 25 == 0 || score > 0) {
            printf("[%4d/%d] mu=%6.1f eta=%.2f m2=%.2f k=%3.0f | clust=%2d mass=%6.1f rms=%5.2f "
                   "compact=%.3f contrast=%5.1f P=%.4f %s score=%.2f\n",
                   trial+1, total, c->mu, c->eta, c->m2, c->kappa,
                   n_clusters, best_mass, best_rms, compactness, contrast,
                   P_abs_max, has_particle?"YES":"no ", score);
            fflush(stdout);
        }

        free(label); free(stack);
        for (int a=0;a<3;a++) free(h_phi[a]);
        trial++;
    }

    /* Report top 20 */
    printf("\n=== TOP 20 RESULTS ===\n");
    printf("rank  score    mu     eta    m2    kappa  clust  mass    rms   compact contrast P_max\n");
    for (int rank = 0; rank < 20 && rank < total; rank++) {
        /* Find rank-th best */
        int bi = -1; double bs = -1;
        for (int i = 0; i < total; i++) {
            if (results[i].score > bs) {
                /* Check not already printed */
                int used = 0;
                for (int r2 = 0; r2 < rank; r2++) {
                    /* Simple: just zero out score after printing */
                }
                if (!used) { bs = results[i].score; bi = i; }
            }
        }
        if (bi < 0) break;
        SweepResult *r = &results[bi];
        printf("%4d  %6.2f  %6.1f  %.2f  %5.2f  %5.0f  %4d  %6.1f  %5.2f  %.3f   %5.1f   %.4f\n",
               rank+1, r->score, r->mu, r->eta, r->m2, r->kappa,
               r->n_clusters, r->best_mass, r->best_rms, r->compactness,
               r->density_contrast, r->P_max);
        r->score = -999;  /* mark as printed */
        fflush(stdout);
    }

    if (best_trial >= 0) {
        SweepResult *b = &results[best_trial];
        printf("\n=== BEST: mu=%.1f eta=%.2f m2=%.2f kappa=%.0f (score=%.2f) ===\n",
               b->mu, b->eta, b->m2, b->kappa, best_score);
        /* Set config to best parameters for subsequent run */
        c->mu = b->mu; c->eta = b->eta; c->m2 = b->m2; c->kappa = b->kappa;
    }

    /* Cleanup */
    for (int a=0;a<3;a++){free(save_phi[a]);free(save_vel[a]);free(save_theta[a]);free(save_tvel[a]);}
    free(results);
}

/* ================================================================
   Chirality pair sweep
   ================================================================ */

static void init_braid_at(Grid *g, const Config *c,
                           double cx, double cy, double cz,
                           double d0, double d1, double d2,
                           double theta_sign, int add) {
    /* Place a braid centered at (cx,cy,cz) with phase offsets (d0,d1,d2).
     * theta_sign: +1 or -1 for chirality.
     * add=0: overwrite, add=1: add to existing field. */
    int N = g->N, NN = N*N;
    double L = g->L, dx = g->dx;
    double kw = PI/L, omega = sqrt(kw*kw + c->m2);
    double sx = 1+c->ellip, sy = 1-c->ellip;
    double inv2R2 = 1.0/(2*c->R_tube*c->R_tube);
    double k_bg = PI/L;
    double delta[3] = {d0, d1, d2};

    for (int i = 0; i < N; i++) { double x = -L + i*dx - cx;
    for (int j = 0; j < N; j++) { double y = -L + j*dx - cy;
    for (int k = 0; k < N; k++) { double z = -L + k*dx - cz;
        long idx = (long)i*NN + j*N + k;
        double r2e = x*x/(sx*sx) + y*y/(sy*sy);
        double env = exp(-r2e * inv2R2);
        for (int a = 0; a < 3; a++) {
            double ph = kw*z + delta[a];
            double val = c->A * env * cos(ph);
            if (add) g->phi[a][idx] += val;
            else     g->phi[a][idx] = val;
            double vval = omega * c->A * env * sin(ph);
            if (add) g->phi_vel[a][idx] += vval;
            else     g->phi_vel[a][idx] = vval;
        }
    }}}

    /* Set theta = theta_sign * gain * curl(phi) */
    double gain = 2.5 * theta_sign;
    double idx2 = 1.0 / (2.0 * dx);
    for (long vi = 0; vi < g->N3; vi++) {
        int i = (int)(vi/NN), j = (int)((vi/N)%N), k = (int)(vi%N);
        int ip=(i+1)%N, im=(i-1+N)%N, jp=(j+1)%N, jm=(j-1+N)%N, kp=(k+1)%N, km=(k-1+N)%N;
        double cx_ = (g->phi[2][(long)i*NN+jp*N+k]-g->phi[2][(long)i*NN+jm*N+k])*idx2
                   - (g->phi[1][(long)i*NN+j*N+kp]-g->phi[1][(long)i*NN+j*N+km])*idx2;
        double cy_ = (g->phi[0][(long)i*NN+j*N+kp]-g->phi[0][(long)i*NN+j*N+km])*idx2
                   - (g->phi[2][(long)ip*NN+j*N+k]-g->phi[2][(long)im*NN+j*N+k])*idx2;
        double cz_ = (g->phi[1][(long)ip*NN+j*N+k]-g->phi[1][(long)im*NN+j*N+k])*idx2
                   - (g->phi[0][(long)i*NN+jp*N+k]-g->phi[0][(long)i*NN+jm*N+k])*idx2;
        if (add) {
            g->theta[0][vi] += gain * cx_;
            g->theta[1][vi] += gain * cy_;
            g->theta[2][vi] += gain * cz_;
        } else {
            g->theta[0][vi] = gain * cx_;
            g->theta[1][vi] = gain * cy_;
            g->theta[2][vi] = gain * cz_;
        }
    }
}

static void run_chiral_sweep(Config *c, Grid *g) {
    int N = g->N;
    long N3 = g->N3;
    double dx = g->dx, dt = g->dt, L = g->L;
    int NN = N*N;
    size_t fb = N3 * sizeof(double);

    int sweep_steps = (int)lround(c->sweep_T / dt);

    /* Chirality methods for particle B (particle A is always standard) */
    /* Standard: delta=(0, 3.0005, 4.4325), theta_sign=-1 */
    double std_d[3] = {0, 3.0005, 4.4325};

    typedef struct {
        const char *name;
        double d[3];          /* phase offsets */
        double theta_sign;    /* +1 or -1 */
        double ellip_mult;    /* multiply ellip by this */
    } ChiralMethod;

    ChiralMethod methods[] = {
        {"same (control)",        {0, 3.0005, 4.4325},       -1,  1},
        {"flip_theta",            {0, 3.0005, 4.4325},       +1,  1},
        {"reverse_delta",         {4.4325, 3.0005, 0},       -1,  1},
        {"negate_delta",          {0, -3.0005, -4.4325},     -1,  1},
        {"flip_theta+rev_delta",  {4.4325, 3.0005, 0},       +1,  1},
        {"flip_theta+neg_delta",  {0, -3.0005, -4.4325},     +1,  1},
        {"flip_ellip",            {0, 3.0005, 4.4325},       -1, -1},
        {"flip_all",              {0, -3.0005, -4.4325},     +1, -1},
        {"pi_shift",              {PI, 3.0005+PI, 4.4325+PI},-1,  1},
        {"pi_shift+flip_theta",   {PI, 3.0005+PI, 4.4325+PI},+1,  1},
    };
    int n_methods = sizeof(methods) / sizeof(methods[0]);

    /* Separation distances to test */
    double separations[] = {6, 8, 10};
    int n_sep = sizeof(separations) / sizeof(double);

    int total = n_methods * n_sep;
    printf("\n=== CHIRALITY PAIR SWEEP: %d methods × %d separations = %d trials ===\n\n",
           n_methods, n_sep, total);
    printf("Particle A (standard): delta=(0, 3.0005, 4.4325), theta_sign=-1\n");
    printf("Sweep T=%.0f, params: mu=%.1f eta=%.2f m2=%.2f kappa=%.0f\n\n",
           c->sweep_T, c->mu, c->eta, c->m2, c->kappa);
    printf("%-25s  sep  clustA  clustB  H_A      H_B      same?  both_alive?\n",
           "method");
    printf("---------------------------------------------------------------------------------\n");
    fflush(stdout);

    double orig_ellip = c->ellip;

    for (int im = 0; im < n_methods; im++)
    for (int is = 0; is < n_sep; is++) {
        double D = separations[is];
        ChiralMethod *m = &methods[im];

        /* Clear grid */
        for (int a = 0; a < 3; a++) {
            memset(g->phi[a], 0, fb);
            memset(g->phi_vel[a], 0, fb);
            memset(g->theta[a], 0, fb);
            memset(g->theta_vel[a], 0, fb);
        }

        /* Add background */
        double k_bg = PI/L, omega_bg = sqrt(k_bg*k_bg + c->m2);
        for (int i=0;i<N;i++) for (int j=0;j<N;j++) for (int k_=0;k_<N;k_++) {
            long idx = (long)i*NN+j*N+k_;
            double z = -L + k_*dx;
            for (int a=0;a<3;a++) {
                double ph = k_bg*z + 2*PI*a/3.0;
                g->phi[a][idx] += c->A_bg * cos(ph);
                g->phi_vel[a][idx] += omega_bg * c->A_bg * sin(ph);
            }
        }

        /* Place braid A (standard) at (-D/2, 0, 0) */
        c->ellip = orig_ellip;
        init_braid_at(g, c, -D/2, 0, 0, std_d[0], std_d[1], std_d[2], -1, 1);

        /* Place braid B (test chirality) at (+D/2, 0, 0) */
        c->ellip = orig_ellip * m->ellip_mult;
        init_braid_at(g, c, D/2, 0, 0, m->d[0], m->d[1], m->d[2], m->theta_sign, 1);
        c->ellip = orig_ellip;

        /* Upload and run */
        gpu_upload(g);
        gpu_set_constants(c, dx);
        if (gpu_has_intermediates) {
            compute_intermediates_kernel<<<gpu_blocks, THREADS_PER_BLOCK>>>(
                d_phi[0],d_phi[1],d_phi[2], d_theta[0],d_theta[1],d_theta[2],
                d_mismatch[0],d_mismatch[1],d_mismatch[2],
                d_harden_Q[0],d_harden_Q[1],d_harden_Q[2]);
        }
        compute_forces_kernel<<<gpu_blocks, THREADS_PER_BLOCK>>>(
            d_phi[0],d_phi[1],d_phi[2], d_theta[0],d_theta[1],d_theta[2],
            d_vel_theta[0],d_vel_theta[1],d_vel_theta[2],
            d_acc_phi[0],d_acc_phi[1],d_acc_phi[2],
            d_acc_theta[0],d_acc_theta[1],d_acc_theta[2],
            d_mismatch[0],d_mismatch[1],d_mismatch[2],
            d_harden_Q[0],d_harden_Q[1],d_harden_Q[2]);

        for (int s = 0; s < sweep_steps; s++)
            gpu_verlet_step(dt);
        cudaDeviceSynchronize();

        /* Download phi and theta for analysis */
        double *h_phi[3], *h_theta[3];
        for (int a=0;a<3;a++) {
            h_phi[a] = (double*)malloc(fb);
            h_theta[a] = (double*)malloc(fb);
            cudaMemcpy(h_phi[a], d_phi[a], fb, cudaMemcpyDeviceToHost);
            cudaMemcpy(h_theta[a], d_theta[a], fb, cudaMemcpyDeviceToHost);
        }

        /* Find clusters and compute H_cross for left (x<0) and right (x>0) halves */
        double H_left = 0, H_right = 0;
        double mass_left = 0, mass_right = 0;
        double phi_max_left = 0, phi_max_right = 0;
        double idx1 = 1.0 / (2.0*dx);

        for (long vi = 0; vi < N3; vi++) {
            int i=(int)(vi/NN), j=(int)((vi/N)%N), k_=(int)(vi%N);
            double x = -L + i*dx;
            int ip=(i+1)%N, im=(i-1+N)%N, jp=(j+1)%N, jm=(j-1+N)%N, kp=(k_+1)%N, km=(k_-1+N)%N;

            double p2 = h_phi[0][vi]*h_phi[0][vi]+h_phi[1][vi]*h_phi[1][vi]+h_phi[2][vi]*h_phi[2][vi];
            double dV = dx*dx*dx;

            /* curl(phi) */
            double cp0 = (h_phi[2][(long)i*NN+jp*N+k_]-h_phi[2][(long)i*NN+jm*N+k_])*idx1
                       - (h_phi[1][(long)i*NN+j*N+kp]-h_phi[1][(long)i*NN+j*N+km])*idx1;
            double cp1 = (h_phi[0][(long)i*NN+j*N+kp]-h_phi[0][(long)i*NN+j*N+km])*idx1
                       - (h_phi[2][(long)ip*NN+j*N+k_]-h_phi[2][(long)im*NN+j*N+k_])*idx1;
            double cp2 = (h_phi[1][(long)ip*NN+j*N+k_]-h_phi[1][(long)im*NN+j*N+k_])*idx1
                       - (h_phi[0][(long)i*NN+jp*N+k_]-h_phi[0][(long)i*NN+jm*N+k_])*idx1;

            /* H_cross = theta · curl(phi) */
            double hc = (h_theta[0][vi]*cp0 + h_theta[1][vi]*cp1 + h_theta[2][vi]*cp2) * dV;

            if (x < 0) {
                H_left += hc;
                mass_left += p2 * dV;
                if (sqrt(p2) > phi_max_left) phi_max_left = sqrt(p2);
            } else {
                H_right += hc;
                mass_right += p2 * dV;
                if (sqrt(p2) > phi_max_right) phi_max_right = sqrt(p2);
            }
        }

        int A_alive = (phi_max_left > 0.3);
        int B_alive = (phi_max_right > 0.3);
        int same_sign = ((H_left < 0) == (H_right < 0));

        printf("%-25s  %3.0f   %5.1f/%c  %5.1f/%c  %+8.1f %+8.1f  %s    %s\n",
               m->name, D,
               phi_max_left, A_alive?'Y':'n',
               phi_max_right, B_alive?'Y':'n',
               H_left, H_right,
               same_sign ? "SAME" : "DIFF",
               (A_alive && B_alive) ? "BOTH" : (A_alive || B_alive) ? "ONE " : "NONE");
        fflush(stdout);

        for (int a=0;a<3;a++) { free(h_phi[a]); free(h_theta[a]); }
    }

    c->ellip = orig_ellip;
    printf("\nChirality sweep complete.\n");
    fflush(stdout);
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
    cfg_validate(&c);

    Grid *g = grid_alloc(&c);
    printf("dx=%.4f dt=%.6f\n\n", g->dx, g->dt);

    do_init(g, &c);

    /* v69: net-charge refusal (§1.2) + mandatory Gauss projection (§5.4),
     * HOST-side on the freshly initialized grid, BEFORE gpu_upload */
    const int gauged = (c.complex_gauge && c.g_gauge != 0.0);
    if (gauged) {
        if (c.bc_type == 2) net_charge_check(g, &c);
        init_gauss_project(g, &c);
    }
    if (c.test_gauge_xform) {
        fprintf(stderr, "ERROR: test_gauge_xform=1 is CPU-only (use scp_sim, not scp_sim_cuda)\n");
        exit(1);
    }

    /* For gradient_pinned BC: save initial state, disable spherical damping */
    if (c.bc_type == 1) {
        grid_save_pinned(g);
        c.damp_width = 0;  c.damp_rate = 0;  /* disable spherical absorb on GPU */
        printf("Gradient BC: pinned %d slabs on each x-face (A_high=%.3f, A_low=%.3f)\n\n",
               c.gradient_margin, c.gradient_A_high, c.gradient_A_low);
    }
    /* bc_type=2 (periodic): CPU applies NO damping; the GPU absorbing kernel
     * runs unconditionally, so zero its parameters (no-op in the kernel). */
    if (c.bc_type == 2) {
        c.damp_width = 0;  c.damp_rate = 0;
    }

    /* GPU setup */
    gpu_alloc(g->N3, (c.alpha_cs != 0 || c.beta_h != 0), c.complex_phi, c.complex_gauge);
    gpu_g_gauge = c.g_gauge;
    gpu_dx_cached = g->dx;
    gpu_set_constants(&c, g->dx);
    cudaMemcpyToSymbol(d_BC_TYPE, &c.bc_type, sizeof(int));
    cudaMemcpyToSymbol(d_GRAD_MARGIN, &c.gradient_margin, sizeof(int));
    cudaMemcpyToSymbol(d_GRAD_A_HIGH, &c.gradient_A_high, sizeof(double));
    cudaMemcpyToSymbol(d_GRAD_A_LOW, &c.gradient_A_low, sizeof(double));
    gpu_upload(g);

    /* Initial force computation on GPU */
    if (gauged) {
        const double STP = 1.0 / (c.g_gauge * g->dx * g->dx * g->dx);
        compute_forces_complex_gauge_kernel<<<gpu_blocks, THREADS_PER_BLOCK>>>(
            d_phi[0],d_phi[1],d_phi[2], d_phi_im[0],d_phi_im[1],d_phi_im[2],
            d_theta[0],d_theta[1],d_theta[2], d_theta_im[0],d_theta_im[1],d_theta_im[2],
            d_th[0],d_th[1],d_th[2],
            d_acc_phi[0],d_acc_phi[1],d_acc_phi[2],
            d_acc_phi_im[0],d_acc_phi_im[1],d_acc_phi_im[2],
            d_acc_theta[0],d_acc_theta[1],d_acc_theta[2],
            d_acc_theta_im[0],d_acc_theta_im[1],d_acc_theta_im[2],
            d_E_acc[0],d_E_acc[1],d_E_acc[2],
            c.g_gauge, STP, 1.0/g->dx);
    } else if (c.complex_phi) {
        compute_forces_complex_kernel<<<gpu_blocks, THREADS_PER_BLOCK>>>(
            d_phi[0],d_phi[1],d_phi[2], d_phi_im[0],d_phi_im[1],d_phi_im[2],
            d_theta[0],d_theta[1],d_theta[2], d_theta_im[0],d_theta_im[1],d_theta_im[2],
            d_acc_phi[0],d_acc_phi[1],d_acc_phi[2],
            d_acc_phi_im[0],d_acc_phi_im[1],d_acc_phi_im[2],
            d_acc_theta[0],d_acc_theta[1],d_acc_theta[2],
            d_acc_theta_im[0],d_acc_theta_im[1],d_acc_theta_im[2]);
    } else {
        if (gpu_has_intermediates) {
            compute_intermediates_kernel<<<gpu_blocks, THREADS_PER_BLOCK>>>(
                d_phi[0],d_phi[1],d_phi[2], d_theta[0],d_theta[1],d_theta[2],
                d_mismatch[0],d_mismatch[1],d_mismatch[2],
                d_harden_Q[0],d_harden_Q[1],d_harden_Q[2]);
        }
        compute_forces_kernel<<<gpu_blocks, THREADS_PER_BLOCK>>>(
            d_phi[0],d_phi[1],d_phi[2], d_theta[0],d_theta[1],d_theta[2],
            d_vel_theta[0],d_vel_theta[1],d_vel_theta[2],
            d_acc_phi[0],d_acc_phi[1],d_acc_phi[2],
            d_acc_theta[0],d_acc_theta[1],d_acc_theta[2],
            d_mismatch[0],d_mismatch[1],d_mismatch[2],
            d_harden_Q[0],d_harden_Q[1],d_harden_Q[2]);
    }
    cudaDeviceSynchronize();

    /* Parameter sweep mode */
    if (c.sweep == 2) {
        run_chiral_sweep(&c, g);
        printf("\nChiral sweep done. Exiting.\n");
        return 0;
    }
    if (c.sweep == 1) {
        run_sweep(&c, g);
        printf("\nSweep complete. Best params set in config. Continuing with normal run...\n\n");
        /* Re-init with best params and re-upload */
        do_init(g, &c);
        gpu_upload(g);
        gpu_set_constants(&c, g->dx);
        if (gpu_has_intermediates) {
            compute_intermediates_kernel<<<gpu_blocks, THREADS_PER_BLOCK>>>(
                d_phi[0],d_phi[1],d_phi[2], d_theta[0],d_theta[1],d_theta[2],
                d_mismatch[0],d_mismatch[1],d_mismatch[2],
                d_harden_Q[0],d_harden_Q[1],d_harden_Q[2]);
        }
        compute_forces_kernel<<<gpu_blocks, THREADS_PER_BLOCK>>>(
            d_phi[0],d_phi[1],d_phi[2], d_theta[0],d_theta[1],d_theta[2],
            d_vel_theta[0],d_vel_theta[1],d_vel_theta[2],
            d_acc_phi[0],d_acc_phi[1],d_acc_phi[2],
            d_acc_theta[0],d_acc_theta[1],d_acc_theta[2],
            d_mismatch[0],d_mismatch[1],d_mismatch[2],
            d_harden_Q[0],d_harden_Q[1],d_harden_Q[2]);
        cudaDeviceSynchronize();
    }

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
    sfa_embed_kvmd(sfa, &c);
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
    if (c.complex_phi) {
        /* v66 imaginary sector (SPEC §7.1): component offset marks Im parts */
        sfa_add_column(sfa,"phiim_x",   sfa_dtype,SFA_POSITION,3);
        sfa_add_column(sfa,"phiim_y",   sfa_dtype,SFA_POSITION,4);
        sfa_add_column(sfa,"phiim_z",   sfa_dtype,SFA_POSITION,5);
        sfa_add_column(sfa,"thetaim_x", sfa_dtype,SFA_ANGLE,3);
        sfa_add_column(sfa,"thetaim_y", sfa_dtype,SFA_ANGLE,4);
        sfa_add_column(sfa,"thetaim_z", sfa_dtype,SFA_ANGLE,5);
        sfa_add_column(sfa,"phiim_vx",  sfa_dtype,SFA_VELOCITY,6);
        sfa_add_column(sfa,"phiim_vy",  sfa_dtype,SFA_VELOCITY,7);
        sfa_add_column(sfa,"phiim_vz",  sfa_dtype,SFA_VELOCITY,8);
        sfa_add_column(sfa,"thetaim_vx",sfa_dtype,SFA_VELOCITY,9);
        sfa_add_column(sfa,"thetaim_vy",sfa_dtype,SFA_VELOCITY,10);
        sfa_add_column(sfa,"thetaim_vz",sfa_dtype,SFA_VELOCITY,11);
    }
    if (c.complex_gauge) {
        /* v69 gauge sector (SPEC §6): links as ANGLE 6-8, E as VELOCITY 12-14 */
        sfa_add_column(sfa,"th_x",sfa_dtype,SFA_ANGLE,   6);
        sfa_add_column(sfa,"th_y",sfa_dtype,SFA_ANGLE,   7);
        sfa_add_column(sfa,"th_z",sfa_dtype,SFA_ANGLE,   8);
        sfa_add_column(sfa,"E_x", sfa_dtype,SFA_VELOCITY,12);
        sfa_add_column(sfa,"E_y", sfa_dtype,SFA_VELOCITY,13);
        sfa_add_column(sfa,"E_z", sfa_dtype,SFA_VELOCITY,14);
    }
    sfa_finalize_header(sfa);
    const char *pn[]={"f16","f32","f64"};
    printf("SFA: %s (%d cols, %s, colzstd)\n\n", c.output,
           c.complex_gauge ? 30 : (c.complex_phi ? 24 : 12), pn[c.precision]);

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
        fstate.phi_im[a] = d_phi_im[a]; fstate.vel_phi_im[a] = d_vel_phi_im[a]; fstate.acc_phi_im[a] = d_acc_phi_im[a];
        fstate.theta_im[a] = d_theta_im[a]; fstate.vel_theta_im[a] = d_vel_theta_im[a]; fstate.acc_theta_im[a] = d_acc_theta_im[a];
        fstate.th[a] = d_th[a]; fstate.Efield[a] = d_E[a]; fstate.E_acc[a] = d_E_acc[a];
    }
    fstate.complex_mode = c.complex_phi;
    fstate.gauge_mode = c.complex_gauge;
    fstate.N3 = g->N3; fstate.N = g->N;
    fstate.L = g->L; fstate.dx = g->dx; fstate.dt = g->dt;

    /* Create shared frame writer — one file, both hooks write through this */
    FrameWriter frame_writer;
    fw_init(&frame_writer, sfa);

    /* Create hooks */
    SnapHookCtx snap_ctx = create_snap_hook(&frame_writer, c.precision, snap_every, g->N3,
                                            c.complex_gauge ? 30 : (c.complex_phi ? 24 : 12));
    start_snap_writer(&snap_ctx);
    register_hook(snap_hook, &snap_ctx);

    /* Diagnostics file */
    FILE *fp = fopen(c.diag_file, "w");
    fprintf(fp, "t\tE_phi_kin\tE_theta_kin\tE_grad\tE_mass\tE_pot\tE_tgrad\tE_tmass\t"
                "E_coupling\tE_total\tphi_max\tP_max\tP_int\ttheta_rms");
    if (c.complex_phi) fprintf(fp, "\tQ_phi\tQ_theta\tQ_total\ts_max\tr_core"
                                   "\tomega_core\tQ_core\tthp1u\tthp1v\tthp2u\tthp2v");
    if (c.complex_gauge) fprintf(fp, "\tgauss_max\tgauss_l2\tE_em\tQ_flux");
    fprintf(fp, "\n");

    cudaEvent_t t_start;
    cudaEventCreate(&t_start);
    cudaEventRecord(t_start);

    DiagHookCtx diag_ctx = create_diag_hook(fp, diag_every, major, g->N3, g->dx, n_steps, t_start, &c);
    diag_ctx.G_offset = g->G_offset;   /* frozen jellium offset from init_gauss_project */
    register_hook(diag_hook, &diag_ctx);

    /* Vector output hook (optional) — uses same FrameWriter */
    VecHookCtx vec_ctx = create_vec_hook(&c, g->dt, g->N3, &frame_writer);
    if (vec_ctx.enabled)
        register_hook(vec_hook, &vec_ctx);

    /* Initial diagnostics (GPU-side reduction, no full download) */
    run_gpu_diagnostics(&fstate, &diag_ctx);
    {
        double *r = diag_ctx.h_results;
        double et0 = r[0]+r[1]+r[2]+r[3]+r[4]+r[5]+r[6]+r[7]
                   + (diag_ctx.gauged ? diag_ctx.e_em : 0.0);
        double trms0 = sqrt(r[11] / ((c.complex_phi ? 6.0 : 3.0) * g->N3));
        diag_ctx.E0 = et0;
        printf("INIT: E_total=%.4e E_pot=%.4f phi_max=%.4f theta_rms=%.3e\n\n", et0, r[4], r[8], trms0);
        fprintf(fp, "%.2f\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e",
                0.0, r[0], r[1], r[2], r[3], r[4], r[5], r[6], r[7], et0, r[8], r[9], r[10], trms0);
        if (c.complex_phi) {
            double Qp, Qt, smax, rcore, ocore, qcore;
            complex_charges_from_results(&diag_ctx, &Qp, &Qt, &smax, &rcore, &ocore, &qcore);
            fprintf(fp, "\t%.12e\t%.12e\t%.12e\t%.6e\t%.6e\t%.12e\t%.12e\t%.6e\t%.6e\t%.6e\t%.6e",
                    Qp, Qt, Qp+Qt, smax, rcore, ocore, qcore,
                    diag_ctx.probes[0], diag_ctx.probes[1], diag_ctx.probes[2], diag_ctx.probes[3]);
            printf("INIT charges: Q_phi=%.6e Q_theta=%.6e s_max=%.4e r_core=%.3f "
                   "omega_core=%.4f Q_core=%.6e\n\n",
                   Qp, Qt, smax, rcore, ocore, qcore);
        }
        if (c.complex_gauge) {
            fprintf(fp, "\t%.6e\t%.6e\t%.12e\t%.12e",
                    diag_ctx.gauss_max, diag_ctx.gauss_l2, diag_ctx.e_em, diag_ctx.q_flux);
            printf("INIT gauge: gauss_max=%.3e gauss_l2=%.3e E_em=%.6e Q_flux=%.6e\n\n",
                   diag_ctx.gauss_max, diag_ctx.gauss_l2, diag_ctx.e_em, diag_ctx.q_flux);
        }
        fprintf(fp, "\n");
        fflush(fp);
    }

    /* t=0 snapshot — wait for voxel write to complete before vec frame */
    snap_hook(0, 0.0, &fstate, &snap_ctx);
    /* Block until the snap writer thread finishes the voxel frame */
    {
        pthread_mutex_lock(&snap_ctx.mutex);
        while (snap_ctx.writer_busy || snap_ctx.writer_has_data)
            pthread_cond_wait(&snap_ctx.cond, &snap_ctx.mutex);
        pthread_mutex_unlock(&snap_ctx.mutex);
    }
    if (vec_ctx.enabled)
        vec_hook(0, 0.0, &fstate, &vec_ctx);

    printf("Async pipeline: snap every %d steps, diag every %d steps, %d total\n\n",
           snap_every, diag_every, n_steps);

    /* BC switch: absorb -> periodic at bc_switch_time */
    int bc_switched = 0;
    int bc_switch_step = (c.bc_switch_time > 0) ? (int)lround(c.bc_switch_time / g->dt) : 0;

    /* Auto-tune interval */
    int tune_every = (c.tune_dt > 0) ? (int)lround(c.tune_dt / g->dt) : 0;
    if (tune_every > 0)
        printf("Auto-tune: every %d steps (%.1f time units)\n\n", tune_every, c.tune_dt);

    /* v65 self-tuning: kappa becomes a dynamical variable that self-organizes to the
       stability edge at fixed conserved charge Q=int|phi|^2 (mechanism A). */
    int st_every = (c.self_tune && c.st_dt > 0) ? (int)lround(c.st_dt / g->dt) : 0;
    double st_Q0 = -1.0;   /* target charge, set at first tune step */
    if (st_every > 0)
        printf("Self-tune: kappa dynamical, every %d steps (%.1f t.u.); "
               "eps=%.3f gamma=%.3f pcrit=%.2f project=%d  kappa0=%.3f\n\n",
               st_every, c.st_dt, c.st_eps, c.st_gamma, c.st_pcrit, c.st_project, c.kappa);

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

        /* Auto-tune: cluster analysis + parameter adjustment */
        if (tune_every > 0 && step % tune_every == 0) {
            cudaDeviceSynchronize();
            inline_cluster_analysis(&c, g->N, g->N3, g->L, g->dx, step * g->dt, step);
        }

        /* v65 self-tuning (mechanism A): kappa -> stability edge at fixed charge.
           The action pull (dE/dkappa>0) drifts kappa down when the core is stable; a
           collapse sensor (P_max) pushes it up; the balance is the SOC critical edge.
           The charge Q=int|phi|^2 (= 2*E_mass/m2) is held fixed by rescaling phi, so the
           self-tuned kappa*(Q) is a function of the conserved charge. */
        if (st_every > 0 && step % st_every == 0) {
            cudaDeviceSynchronize();
            run_gpu_diagnostics(&fstate, &diag_ctx);
            double Pmax  = diag_ctx.h_results[9];
            double Emass = diag_ctx.h_results[3];
            double Q = (c.m2 > 1e-12) ? 2.0 * Emass / c.m2 : Emass;
            if (st_Q0 < 0) st_Q0 = (c.st_charge > 0 ? c.st_charge : Q);  /* set target once */
            if (c.st_project && c.st_ptarget <= 0 && Q > 1e-12) {        /* charge projection (SOC mode only) */
                double f = sqrt(st_Q0 / Q);
                rescale_phi_kernel<<<gpu_blocks, THREADS_PER_BLOCK>>>(
                    d_phi[0], d_phi[1], d_phi[2],
                    d_vel_phi[0], d_vel_phi[1], d_vel_phi[2], f, g->N3);
            }
            int nanflag = (Pmax != Pmax) || (Pmax > 1e30);
            if (c.st_ptarget > 0.0) {
                /* Bidirectional density homeostasis: regulate P_max -> st_ptarget. Higher
                   kappa lowers the equilibrium density, so kappa *= (1 + gain*err) with
                   err=(P_max-target)/target is stable NEGATIVE feedback. Dissolution
                   (P_max<target) -> kappa DOWN -> attraction re-concentrates the particle;
                   collapse (P_max>target) -> kappa UP -> saturation spreads it. */
                double err = nanflag ? 1.0 : (Pmax - c.st_ptarget) / c.st_ptarget;
                if (err >  1.0) err =  1.0;
                if (err < -1.0) err = -1.0;
                c.kappa *= (1.0 + c.st_gain * err);
                if (c.kappa < 1e-3) c.kappa = 1e-3;
                cudaMemcpyToSymbol(d_KAPPA, &c.kappa, sizeof(double));
                printf("[self-tune] t=%6.1f kappa=%9.4f Pmax=%10.3e target=%.4f err=%+.3f Q=%9.3e\n",
                       step * g->dt, c.kappa, Pmax, c.st_ptarget, err, Q);
            } else {
                int collapse = nanflag || (Pmax > c.st_pcrit);
                c.kappa *= collapse ? (1.0 + c.st_gamma) : (1.0 - c.st_eps);  /* SOC update */
                if (c.kappa < 1e-3) c.kappa = 1e-3;
                cudaMemcpyToSymbol(d_KAPPA, &c.kappa, sizeof(double));
                printf("[self-tune] t=%6.1f kappa=%9.4f Pmax=%10.3e Q=%9.3e %s\n",
                       step * g->dt, c.kappa, Pmax, Q, collapse ? "COLLAPSE" : "stable");
            }
            fflush(stdout);
        }

        if (c.bc_type == 1) {
            gradient_bc_kernel<<<gpu_blocks, THREADS_PER_BLOCK>>>(
                d_phi[0], d_phi[1], d_phi[2],
                d_vel_phi[0], d_vel_phi[1], d_vel_phi[2],
                d_acc_phi[0], d_acc_phi[1], d_acc_phi[2],
                d_theta[0], d_theta[1], d_theta[2],
                d_vel_theta[0], d_vel_theta[1], d_vel_theta[2],
                d_acc_theta[0], d_acc_theta[1], d_acc_theta[2]);
            if (c.complex_phi)
                gradient_bc_im_kernel<<<gpu_blocks, THREADS_PER_BLOCK>>>(
                    d_phi_im[0], d_phi_im[1], d_phi_im[2],
                    d_vel_phi_im[0], d_vel_phi_im[1], d_vel_phi_im[2],
                    d_acc_phi_im[0], d_acc_phi_im[1], d_acc_phi_im[2],
                    d_theta_im[0], d_theta_im[1], d_theta_im[2],
                    d_vel_theta_im[0], d_vel_theta_im[1], d_vel_theta_im[2],
                    d_acc_theta_im[0], d_acc_theta_im[1], d_acc_theta_im[2]);
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
            /* Dispatch snap hook first */
            snap_hook(step, t, &fstate, &snap_ctx);
            /* If snap fired, wait for writer to finish before vec writes */
            if (step % snap_ctx.snap_every == 0) {
                pthread_mutex_lock(&snap_ctx.mutex);
                while (snap_ctx.writer_busy || snap_ctx.writer_has_data)
                    pthread_cond_wait(&snap_ctx.cond, &snap_ctx.mutex);
                pthread_mutex_unlock(&snap_ctx.mutex);
            }
            /* Now dispatch remaining hooks (diag, vec) */
            for (int h = 0; h < n_hooks; h++)
                if (hooks[h].fn != snap_hook)
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

    /* Wait for writer to finish, then destroy hooks and frame writer */
    destroy_snap_hook(&snap_ctx);
    destroy_vec_hook(&vec_ctx);
    fw_destroy(&frame_writer);
    uint32_t nf = sfa->total_frames;
    sfa_close(sfa);

    /* Final diagnostics (GPU-side, no download) */
    run_gpu_diagnostics(&fstate, &diag_ctx);
    {
        double *r = diag_ctx.h_results;
        double et = r[0]+r[1]+r[2]+r[3]+r[4]+r[5]+r[6]+r[7]
                  + (diag_ctx.gauged ? diag_ctx.e_em : 0.0);
        double trms = sqrt(r[11] / ((c.complex_phi ? 6.0 : 3.0) * g->N3));

        cudaEvent_t t_stop;
        cudaEventCreate(&t_stop);
        cudaEventRecord(t_stop);
        cudaEventSynchronize(t_stop);
        float total_ms; cudaEventElapsedTime(&total_ms, t_start, t_stop);
        cudaEventDestroy(t_stop);

        printf("\n=== COMPLETE (GPU, async pipeline) ===\n");
        printf("E_total=%.4e (drift %.3f%%) E_pot=%.4f\n", et, 100*(et-diag_ctx.E0)/(fabs(diag_ctx.E0)+1e-30), r[4]);
        printf("phi_max=%.4f theta_rms=%.3e\n", r[8], trms);
        if (c.complex_phi) {
            double Qp, Qt, smax, rcore, ocore, qcore;
            complex_charges_from_results(&diag_ctx, &Qp, &Qt, &smax, &rcore, &ocore, &qcore);
            printf("Q_phi=%.6e Q_theta=%.6e Q_total=%.6e s_max=%.4e r_core=%.3f "
                   "omega_core=%.4f Q_core=%.6e\n",
                   Qp, Qt, Qp+Qt, smax, rcore, ocore, qcore);
        }
        if (c.complex_gauge)
            printf("gauss_max=%.3e gauss_l2=%.3e E_em=%.6e Q_flux=%.6e\n",
                   diag_ctx.gauss_max, diag_ctx.gauss_l2, diag_ctx.e_em, diag_ctx.q_flux);
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
