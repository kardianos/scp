/*  scp_sim.cu — Unified 6-field Cosserat simulation kernel (CUDA)
 *
 *  GPU version of scp_sim.c — same config, same SFA output, same physics.
 *  Config parsing, SFA I/O, and init run on CPU. Physics runs on GPU.
 *
 *  Build: nvcc -O3 -arch=sm_70 -o scp_sim_cuda scp_sim.cu -lzstd -lm
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
    int mode;
    double inv_alpha, inv_beta, kappa_gamma;
    double damp_width, damp_rate;
    char init[32];
    double A, sigma, A_bg, ellip, R_tube;
    double delta[3];
    char init_sfa[512];
    int init_frame;
    char init_exec[1024];
    char output[512];
    char diag_file[512];
    int precision;
    double snap_dt, diag_dt;
} Config;

static Config cfg_defaults(void) {
    Config c = {};
    c.N = 128;  c.L = 10.0;  c.T = 200.0;  c.dt_factor = 0.025;
    c.m2 = 2.25;  c.mtheta2 = 0.0;  c.eta = 0.5;
    c.mu = -41.345;  c.kappa = 50.0;
    c.mode = 0;  c.inv_alpha = 2.25;  c.inv_beta = 5.0;  c.kappa_gamma = 2.0;
    c.damp_width = 3.0;  c.damp_rate = 0.01;
    strcpy(c.init, "oscillon");
    c.A = 0.8;  c.sigma = 3.0;  c.A_bg = 0.1;
    c.ellip = 0.3325;  c.R_tube = 3.0;
    c.delta[0] = 0.0;  c.delta[1] = 3.0005;  c.delta[2] = 4.4325;
    c.init_frame = -1;
    strcpy(c.output, "output.sfa");
    strcpy(c.diag_file, "diag.tsv");
    c.precision = 1;
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
    else if (!strcmp(key,"mode"))        c->mode = atoi(val);
    else if (!strcmp(key,"inv_alpha"))   c->inv_alpha = atof(val);
    else if (!strcmp(key,"inv_beta"))    c->inv_beta = atof(val);
    else if (!strcmp(key,"kappa_gamma")) c->kappa_gamma = atof(val);
    else if (!strcmp(key,"damp_width"))  c->damp_width = atof(val);
    else if (!strcmp(key,"damp_rate"))   c->damp_rate = atof(val);
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
    printf("=== scp_sim CUDA: Unified 6-field Cosserat Kernel (GPU) ===\n");
    printf("d²φ/dt² = ∇²φ - m²φ - V'(P) + η·curl(θ)\n");
    printf("d²θ/dt² = ∇²θ - m_θ²θ + η·curl(φ)\n\n");
    printf("Grid:    N=%d L=%.1f T=%.0f dt_factor=%.4f\n", c->N, c->L, c->T, c->dt_factor);
    printf("Physics: m²=%.4f m_θ²=%.4f η=%.3f μ=%.3f κ=%.1f\n",
           c->m2, c->mtheta2, c->eta, c->mu, c->kappa);
    printf("Mode:    %d", c->mode);
    if (c->mode == 1) printf(" (inverse: α=%.3f β=%.3f)", c->inv_alpha, c->inv_beta);
    if (c->mode == 3) printf(" (density-κ: γ=%.3f)", c->kappa_gamma);
    printf("\nBC:      absorbing sphere (width=%.1f rate=%.4f)\n", c->damp_width, c->damp_rate);
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
   Host grid (for init, diagnostics, SFA I/O)
   ================================================================ */

typedef struct {
    double *mem;
    double *phi[NFIELDS], *phi_vel[NFIELDS], *phi_acc[NFIELDS];
    double *theta[NFIELDS], *theta_vel[NFIELDS], *theta_acc[NFIELDS];
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
static void grid_free(Grid *g) { free(g->mem); free(g); }

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

static void do_init(Grid *g, const Config *c) {
    if      (!strcmp(c->init, "oscillon")) init_oscillon(g, c);
    else if (!strcmp(c->init, "braid"))    init_braid(g, c);
    else if (!strcmp(c->init, "sfa"))      init_from_sfa(g, c);
    else { fprintf(stderr, "FATAL: unknown init '%s'\n", c->init); exit(1); }
}

/* ================================================================
   GPU constant memory
   ================================================================ */

__constant__ double d_MU, d_KAPPA, d_MASS2, d_MTHETA2, d_ETA;
__constant__ double d_inv_alpha, d_inv_beta, d_kappa_gamma;
__constant__ double d_idx2, d_idx1;
__constant__ double d_L, d_dx, d_DAMP_WIDTH, d_DAMP_RATE;
__constant__ int d_N, d_NN, d_MODE;
__constant__ long d_N3;

/* ================================================================
   GPU kernels
   ================================================================ */

__global__ void compute_forces_kernel(
    const double *phi0, const double *phi1, const double *phi2,
    const double *theta0, const double *theta1, const double *theta2,
    double *acc_phi0, double *acc_phi1, double *acc_phi2,
    double *acc_theta0, double *acc_theta1, double *acc_theta2)
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
    double acc_p[3], acc_t[3];

    for (int a = 0; a < 3; a++) {
        /* Laplacian */
        double lap = (phi[a][n_ip]+phi[a][n_im]+phi[a][n_jp]+phi[a][n_jm]
                     +phi[a][n_kp]+phi[a][n_km]-6.0*phi[a][idx]) * d_idx2;
        double dPda = (a==0)?p1*p2:(a==1)?p0*p2:p0*p1;

        /* curl(theta)_a */
        double ct;
        if (a==0) ct = (th[2][n_jp]-th[2][n_jm]-th[1][n_kp]+th[1][n_km])*d_idx1;
        else if (a==1) ct = (th[0][n_kp]-th[0][n_km]-th[2][n_ip]+th[2][n_im])*d_idx1;
        else ct = (th[1][n_ip]-th[1][n_im]-th[0][n_jp]+th[0][n_jm])*d_idx1;

        acc_p[a] = lap - me2*phi[a][idx] - dVdP*dPda - t2c*phi[a][idx] + d_ETA*ct;
    }

    for (int a = 0; a < 3; a++) {
        double lapt = (th[a][n_ip]+th[a][n_im]+th[a][n_jp]+th[a][n_jm]
                      +th[a][n_kp]+th[a][n_km]-6.0*th[a][idx]) * d_idx2;
        double cp;
        if (a==0) cp = (phi[2][n_jp]-phi[2][n_jm]-phi[1][n_kp]+phi[1][n_km])*d_idx1;
        else if (a==1) cp = (phi[0][n_kp]-phi[0][n_km]-phi[2][n_ip]+phi[2][n_im])*d_idx1;
        else cp = (phi[1][n_ip]-phi[1][n_im]-phi[0][n_jp]+phi[0][n_jm])*d_idx1;

        acc_t[a] = lapt - d_MTHETA2*th[a][idx] + d_ETA*cp;
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

__global__ void downcast_f64_to_f32_kernel(const double *src, float *dst, long n) {
    long idx = (long)blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n) return;
    dst[idx] = (float)src[idx];
}

/* ================================================================
   GPU memory management
   ================================================================ */

static double *d_phi[3], *d_vel_phi[3], *d_acc_phi[3];
static double *d_theta[3], *d_vel_theta[3], *d_acc_theta[3];
static float *d_f32_buf;
static int gpu_blocks;

static void gpu_alloc(long N3) {
    size_t bytes = N3 * sizeof(double);
    for (int a = 0; a < 3; a++) {
        cudaMalloc(&d_phi[a], bytes);       cudaMalloc(&d_vel_phi[a], bytes);   cudaMalloc(&d_acc_phi[a], bytes);
        cudaMalloc(&d_theta[a], bytes);     cudaMalloc(&d_vel_theta[a], bytes); cudaMalloc(&d_acc_theta[a], bytes);
    }
    cudaMalloc(&d_f32_buf, N3 * sizeof(float));
    gpu_blocks = (int)((N3 + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK);
    printf("GPU: allocated %.2f GB, %d blocks × %d threads\n",
           (18.0*bytes + N3*sizeof(float))/1e9, gpu_blocks, THREADS_PER_BLOCK);
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
    cudaFree(d_f32_buf);
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
    /* Forces */
    compute_forces_kernel<<<gpu_blocks, THREADS_PER_BLOCK>>>(
        d_phi[0], d_phi[1], d_phi[2],
        d_theta[0], d_theta[1], d_theta[2],
        d_acc_phi[0], d_acc_phi[1], d_acc_phi[2],
        d_acc_theta[0], d_acc_theta[1], d_acc_theta[2]);
    /* Half-kick */
    verlet_halfkick_kernel<<<gpu_blocks, THREADS_PER_BLOCK>>>(
        d_vel_phi[0], d_vel_phi[1], d_vel_phi[2],
        d_vel_theta[0], d_vel_theta[1], d_vel_theta[2],
        d_acc_phi[0], d_acc_phi[1], d_acc_phi[2],
        d_acc_theta[0], d_acc_theta[1], d_acc_theta[2], hdt);
    /* Boundary damping */
    absorbing_boundary_kernel<<<gpu_blocks, THREADS_PER_BLOCK>>>(
        d_vel_phi[0], d_vel_phi[1], d_vel_phi[2],
        d_vel_theta[0], d_vel_theta[1], d_vel_theta[2]);
}

/* ================================================================
   Diagnostics (host-side, same as CPU version)
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
        double P=p0*p1*p2, P2=P*P;
        if (c->mode==3) {
            double D=1.0+c->kappa_gamma*sig+c->kappa*P2;
            s_ep+=(c->mu/2.0)*P2*(1.0+c->kappa_gamma*sig)/D*dV;
        } else s_ep+=(c->mu/2.0)*P2/(1.0+keff*P2)*dV;
        double Pa=fabs(P); if(Pa>s_Pm) s_Pm=Pa;
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

    /* GPU setup */
    gpu_alloc(g->N3);
    gpu_set_constants(&c, g->dx);
    gpu_upload(g);

    /* Initial force computation on GPU */
    compute_forces_kernel<<<gpu_blocks, THREADS_PER_BLOCK>>>(
        d_phi[0],d_phi[1],d_phi[2], d_theta[0],d_theta[1],d_theta[2],
        d_acc_phi[0],d_acc_phi[1],d_acc_phi[2],
        d_acc_theta[0],d_acc_theta[1],d_acc_theta[2]);
    cudaDeviceSynchronize();

    /* SFA archive with KVMD */
    uint8_t sfa_dtype = (c.precision==0)?SFA_F16:(c.precision==1)?SFA_F32:SFA_F64;
    SFA *sfa = sfa_create(c.output, c.N, c.N, c.N, c.L, c.L, c.L, g->dt);
    {
        char vN[32],vL[32],vT[32],vdt[32],vm[32],vmt[32],veta[32],vmu[32],vkappa[32];
        char vmode[32],via[32],vib[32],vkg[32],vdw[32],vdr[32],vprec[32],vdelta[64];
        snprintf(vN,32,"%d",c.N); snprintf(vL,32,"%.6f",c.L); snprintf(vT,32,"%.6f",c.T);
        snprintf(vdt,32,"%.6f",c.dt_factor);
        snprintf(vm,32,"%.6f",sqrt(c.m2)); snprintf(vmt,32,"%.6f",sqrt(c.mtheta2));
        snprintf(veta,32,"%.6f",c.eta); snprintf(vmu,32,"%.6f",c.mu);
        snprintf(vkappa,32,"%.6f",c.kappa); snprintf(vmode,32,"%d",c.mode);
        snprintf(via,32,"%.6f",c.inv_alpha); snprintf(vib,32,"%.6f",c.inv_beta);
        snprintf(vkg,32,"%.6f",c.kappa_gamma);
        snprintf(vdw,32,"%.6f",c.damp_width); snprintf(vdr,32,"%.6f",c.damp_rate);
        snprintf(vprec,32,"%s",(const char*[]){"f16","f32","f64"}[c.precision]);
        snprintf(vdelta,64,"%.6f,%.6f,%.6f",c.delta[0],c.delta[1],c.delta[2]);
        const char *keys[]={"N","L","T","dt_factor","m","m_theta","eta","mu","kappa",
                            "mode","inv_alpha","inv_beta","kappa_gamma",
                            "damp_width","damp_rate","precision","delta"};
        const char *vals[]={vN,vL,vT,vdt,vm,vmt,veta,vmu,vkappa,
                            vmode,via,vib,vkg,vdw,vdr,vprec,vdelta};
        sfa_add_kvmd(sfa, 0, 0xFFFFFFFF, 0xFFFFFFFF, keys, vals, 17);
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
    printf("SFA: %s (12 cols, %s, BSS+zstd)\n\n", c.output, pn[c.precision]);

    /* t=0 snapshot */
    sfa_snap_gpu(sfa, g, 0.0, c.precision);

    /* Diagnostics */
    FILE *fp = fopen(c.diag_file, "w");
    fprintf(fp, "t\tE_phi_kin\tE_theta_kin\tE_grad\tE_mass\tE_pot\tE_tgrad\tE_tmass\t"
                "E_coupling\tE_total\tphi_max\tP_max\tP_int\ttheta_rms\n");

    int n_steps=(int)(c.T/g->dt);
    int diag_every=(int)(c.diag_dt/g->dt); if(diag_every<1) diag_every=1;
    int snap_every=(int)(c.snap_dt/g->dt); if(snap_every<1) snap_every=1;

    /* Download for initial diagnostics */
    gpu_download(g);
    double epk,etk,eg,em,ep,etg,etm,ec,et,pm,Pm;
    compute_energy(g,&c,&epk,&etk,&eg,&em,&ep,&etg,&etm,&ec,&et,&pm,&Pm);
    double Pint0=P_integrated(g), trms0=theta_rms(g), E0=et;
    printf("INIT: E_total=%.4e E_pot=%.4f phi_max=%.4f theta_rms=%.3e\n\n", et,ep,pm,trms0);
    fprintf(fp,"%.2f\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\n",
            0.0,epk,etk,eg,em,ep,etg,etm,ec,et,pm,Pm,Pint0,trms0);

    cudaEvent_t t_start, t_stop;
    cudaEventCreate(&t_start); cudaEventCreate(&t_stop);
    cudaEventRecord(t_start);

    int major = diag_every * 25; if(major<1) major=1;

    for (int step = 1; step <= n_steps; step++) {
        gpu_verlet_step(g->dt);

        int need_diag = (step % diag_every == 0);
        int need_snap = (step % snap_every == 0);

        if (need_snap) {
            double t = step * g->dt;
            sfa_snap_gpu(sfa, g, t, c.precision);
        }

        if (need_diag) {
            double t = step * g->dt;
            if (!need_snap) { cudaDeviceSynchronize(); gpu_download(g); }
            compute_energy(g,&c,&epk,&etk,&eg,&em,&ep,&etg,&etm,&ec,&et,&pm,&Pm);
            double Pint=P_integrated(g), trms=theta_rms(g);
            fprintf(fp,"%.2f\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\n",
                    t,epk,etk,eg,em,ep,etg,etm,ec,et,pm,Pm,Pint,trms);
            fflush(fp);

            if (step % major == 0) {
                cudaEventRecord(t_stop); cudaEventSynchronize(t_stop);
                float ms; cudaEventElapsedTime(&ms, t_start, t_stop);
                double drift=100*(et-E0)/(fabs(E0)+1e-30);
                printf("t=%7.1f E=%.3e (drift %+.3f%%) Ep=%.1f phi=%.3f θ_rms=%.2e "
                       "[%.0f%% %.1fs %.2fms/step]\n",
                       t,et,drift,ep,pm,trms,100.0*step/n_steps,ms/1000,ms/step);
                fflush(stdout);
            }
        }
    }

    /* Final frame */
    sfa_snap_gpu(sfa, g, n_steps*g->dt, c.precision);
    uint32_t nf = sfa->total_frames;
    sfa_close(sfa);

    /* Final summary */
    compute_energy(g,&c,&epk,&etk,&eg,&em,&ep,&etg,&etm,&ec,&et,&pm,&Pm);
    double trms=theta_rms(g);
    cudaEventRecord(t_stop); cudaEventSynchronize(t_stop);
    float total_ms; cudaEventElapsedTime(&total_ms, t_start, t_stop);

    printf("\n=== COMPLETE (GPU) ===\n");
    printf("E_total=%.4e (drift %.3f%%) E_pot=%.4f\n", et, 100*(et-E0)/(fabs(E0)+1e-30), ep);
    printf("phi_max=%.4f theta_rms=%.3e\n", pm, trms);
    printf("SFA: %s (%u frames)\n", c.output, nf);
    printf("[%s] theta_rms grew: %.2e -> %.2e\n", (trms>trms0+1e-10)?"OK":"WARN", trms0, trms);
    printf("Wall: %.1fs (%.1f min) %.2fms/step\n", total_ms/1000, total_ms/60000, total_ms/n_steps);

    fclose(fp);
    gpu_free();
    grid_free(g);
    cudaEventDestroy(t_start); cudaEventDestroy(t_stop);
    return 0;
}
