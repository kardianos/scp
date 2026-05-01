/*  scp_sim.c — Unified 6-field Cosserat simulation kernel
 *
 *  Reads ALL parameters from a config file. No hardcoded physics.
 *  Full 6-field Cosserat equation: 3 phi + 3 theta with curl coupling.
 *
 *  Features:
 *    - 12-column SFA output (fields + velocities) for exact restart
 *    - Configurable output precision: f16, f32, f64
 *    - Init from: built-in (braid/oscillon), SFA file (any frame/dtype), exec
 *    - Mass coupling modes: constant (0), inverse (1), density-kappa (3)
 *    - Spherical absorbing boundary conditions
 *    - Diagnostics: energy decomposition, theta_rms, aspect ratio
 *
 *  Build: gcc -O3 -march=native -fopenmp -o scp_sim scp_sim.c -lzstd -lm
 *  Run:   OMP_NUM_THREADS=8 ./scp_sim config.cfg [-key value ...]
 */

#define SFA_IMPLEMENTATION
#include "../format/sfa.h"
#include "scp_config.h"

#include <omp.h>
#include <sys/stat.h>

/* Config, f16 helpers, and constants defined in scp_config.h */


/* ================================================================
   Grid: 18 arrays (6 fields × {val, vel, acc})
   ================================================================ */

typedef struct {
    double *mem;
    double *phi[NFIELDS], *phi_vel[NFIELDS], *phi_acc[NFIELDS];
    double *theta[NFIELDS], *theta_vel[NFIELDS], *theta_acc[NFIELDS];
    double *mismatch[NFIELDS];  /* Cosserat mismatch: M_a = curl(φ)_a/2 - θ_a */
    double *harden_Q[NFIELDS]; /* Hardening vector: Q_a = (β/2)|θ|²curl(φ)_a */
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
    /* Temporary arrays for two-pass force computations */
    for (int a = 0; a < NFIELDS; a++) {
        g->mismatch[a] = calloc(g->N3, sizeof(double));
        g->harden_Q[a] = calloc(g->N3, sizeof(double));
    }
    return g;
}

static void grid_free(Grid *g) {
    for (int a = 0; a < NFIELDS; a++) { free(g->mismatch[a]); free(g->harden_Q[a]); }
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

/* Initialization and KVMD embedding from shared header */
#include "scp_init.h"

/* Forces start below — init functions are in scp_init.h */

/*  init_oscillon through do_init deleted — now in scp_init.h */

#if 0 /* deleted init block */
static void DELETED_init_oscillon(Grid *g, const Config *c) {
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
    for (uint32_t col = 0; col < sfa->n_columns; col++) {
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
        printf("  WARNING: no velocity data — starting from rest (cold restart)\n");
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
    snprintf(tmp.init_sfa, sizeof(tmp.init_sfa), "%s", tmppath);
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
        int i=(int)(idx/NN), k=(int)(idx%N);
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
    for (uint32_t col=0; col<tmpl->n_columns; col++) {
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
#endif /* deleted init block */

/* ================================================================
   Forces: V44 + chiral helicity (C2) + Cosserat strain (C3)

   Two-pass when alpha_cs > 0:
   Pass 1: compute mismatch M_a = curl(φ)_a/2 - θ_a at all voxels
   Pass 2: compute forces using M (for theta) and curl(M) (for phi)
   ================================================================ */

static void compute_forces(Grid *g, const Config *c) {
    const int N = g->N, NN = N * N;
    const long N3 = g->N3;
    const double idx2 = 1.0 / (g->dx * g->dx);
    const double idx1 = 1.0 / (2.0 * g->dx);
    const double MU = c->mu, KAPPA = c->kappa, MASS2 = c->m2;
    const double MTHETA2 = c->mtheta2, ETA = c->eta;
    const double KH = c->kappa_h, ACS = c->alpha_cs, BH = c->beta_h;
    const int MODE = c->mode;
    const double KG = c->kappa_gamma, IA = c->inv_alpha, IB = c->inv_beta;

    /* === Pass 1: Compute temporary fields at all voxels === */
    if (ACS != 0 || BH != 0) {
        #pragma omp parallel for schedule(static)
        for (long idx = 0; idx < N3; idx++) {
            int i=(int)(idx/NN), j=(int)((idx/N)%N), k=(int)(idx%N);
            int ip=(i+1)%N,im=(i-1+N)%N,jp=(j+1)%N,jm=(j-1+N)%N,kp=(k+1)%N,km=(k-1+N)%N;
            long n_ip=(long)ip*NN+j*N+k, n_im=(long)im*NN+j*N+k;
            long n_jp=(long)i*NN+jp*N+k, n_jm=(long)i*NN+jm*N+k;
            long n_kp=(long)i*NN+j*N+kp, n_km=(long)i*NN+j*N+km;

            double cp0 = curl_component(g->phi, 0, n_ip,n_im,n_jp,n_jm,n_kp,n_km, idx1);
            double cp1 = curl_component(g->phi, 1, n_ip,n_im,n_jp,n_jm,n_kp,n_km, idx1);
            double cp2 = curl_component(g->phi, 2, n_ip,n_im,n_jp,n_jm,n_kp,n_km, idx1);

            /* Cosserat mismatch: M_a = curl(φ)_a/2 - θ_a */
            if (ACS != 0) {
                g->mismatch[0][idx] = cp0 * 0.5 - g->theta[0][idx];
                g->mismatch[1][idx] = cp1 * 0.5 - g->theta[1][idx];
                g->mismatch[2][idx] = cp2 * 0.5 - g->theta[2][idx];
            }

            /* Hardening vector: Q_a = (β/2)|θ|² × curl(φ)_a
             * (Maxima-verified: F_phi = -2 curl(Q), so Q is the intermediate) */
            if (BH != 0) {
                double t0=g->theta[0][idx], t1=g->theta[1][idx], t2=g->theta[2][idx];
                double T2 = t0*t0 + t1*t1 + t2*t2;
                double coeff = 0.5 * BH * T2;
                g->harden_Q[0][idx] = coeff * cp0;
                g->harden_Q[1][idx] = coeff * cp1;
                g->harden_Q[2][idx] = coeff * cp2;
            }
        }
    }

    /* === Pass 2: Compute all forces === */
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
        double dVdP = MU * P / (den*den);

        double t2c = 0;
        if (MODE == 3 && KG > 0) {
            double D = 1.0 + KG*sig + KAPPA*P2;
            t2c = MU * KG * KAPPA * P2 * P2 / (D*D);
        }

        /* curl(φ) at this voxel — needed for chiral and Cosserat */
        double curl_phi[3];
        curl_phi[0] = (g->phi[2][n_jp]-g->phi[2][n_jm]-g->phi[1][n_kp]+g->phi[1][n_km])*idx1;
        curl_phi[1] = (g->phi[0][n_kp]-g->phi[0][n_km]-g->phi[2][n_ip]+g->phi[2][n_im])*idx1;
        curl_phi[2] = (g->phi[1][n_ip]-g->phi[1][n_im]-g->phi[0][n_jp]+g->phi[0][n_jm])*idx1;

        /* Chiral helicity: κ_h P² φ·curl(φ) */
        double h = p0*curl_phi[0] + p1*curl_phi[1] + p2*curl_phi[2];
        double P2_ip=0,P2_im=0,P2_jp=0,P2_jm=0,P2_kp=0,P2_km=0;
        double dP2dx=0, dP2dy=0, dP2dz=0;
        if (KH != 0) {
            double Pip=g->phi[0][n_ip]*g->phi[1][n_ip]*g->phi[2][n_ip];
            double Pim=g->phi[0][n_im]*g->phi[1][n_im]*g->phi[2][n_im];
            double Pjp=g->phi[0][n_jp]*g->phi[1][n_jp]*g->phi[2][n_jp];
            double Pjm=g->phi[0][n_jm]*g->phi[1][n_jm]*g->phi[2][n_jm];
            double Pkp=g->phi[0][n_kp]*g->phi[1][n_kp]*g->phi[2][n_kp];
            double Pkm=g->phi[0][n_km]*g->phi[1][n_km]*g->phi[2][n_km];
            P2_ip=Pip*Pip; P2_im=Pim*Pim; P2_jp=Pjp*Pjp;
            P2_jm=Pjm*Pjm; P2_kp=Pkp*Pkp; P2_km=Pkm*Pkm;
            dP2dx=(P2_ip-P2_im)*idx1; dP2dy=(P2_jp-P2_jm)*idx1; dP2dz=(P2_kp-P2_km)*idx1;
        }

        /* Cosserat strain: curl(M) for phi force, M for theta force */
        double curl_M[3] = {0, 0, 0};
        if (ACS != 0) {
            curl_M[0] = (g->mismatch[2][n_jp]-g->mismatch[2][n_jm]
                        -g->mismatch[1][n_kp]+g->mismatch[1][n_km])*idx1;
            curl_M[1] = (g->mismatch[0][n_kp]-g->mismatch[0][n_km]
                        -g->mismatch[2][n_ip]+g->mismatch[2][n_im])*idx1;
            curl_M[2] = (g->mismatch[1][n_ip]-g->mismatch[1][n_im]
                        -g->mismatch[0][n_jp]+g->mismatch[0][n_jm])*idx1;
        }

        /* Curl-squared hardening: curl(Q) for phi force
         * Q_a = (β/2)|θ|²curl(φ)_a (computed in pass 1)
         * F_phi_a = -2 curl(Q)_a  (Maxima-verified sign) */
        double curl_Q[3] = {0, 0, 0};
        if (BH != 0) {
            curl_Q[0] = (g->harden_Q[2][n_jp]-g->harden_Q[2][n_jm]
                        -g->harden_Q[1][n_kp]+g->harden_Q[1][n_km])*idx1;
            curl_Q[1] = (g->harden_Q[0][n_kp]-g->harden_Q[0][n_km]
                        -g->harden_Q[2][n_ip]+g->harden_Q[2][n_im])*idx1;
            curl_Q[2] = (g->harden_Q[1][n_ip]-g->harden_Q[1][n_im]
                        -g->harden_Q[0][n_jp]+g->harden_Q[0][n_jm])*idx1;
        }

        /* |curl(φ)|² at this voxel — needed for theta hardening force */
        double curl_sq = curl_phi[0]*curl_phi[0] + curl_phi[1]*curl_phi[1] + curl_phi[2]*curl_phi[2];

        /* --- Phi forces --- */
        for (int a = 0; a < NFIELDS; a++) {
            double lap = (g->phi[a][n_ip]+g->phi[a][n_im]+g->phi[a][n_jp]
                        +g->phi[a][n_jm]+g->phi[a][n_kp]+g->phi[a][n_km]
                        -6.0*g->phi[a][idx]) * idx2;
            double dPda = (a==0)?p1*p2:(a==1)?p0*p2:p0*p1;
            double ct = curl_component(g->theta, a, n_ip,n_im,n_jp,n_jm,n_kp,n_km, idx1);

            /* Chiral force: κ_h × [2P dPda h + 2P² curl(φ)_a - ∇P² cross terms]
             * The t3 sign is MINUS (Maxima-verified EL derivation).
             * F_a = dL/dφ_a - d_j(dL/d(d_j φ_a))
             *      = kh*(2P dPda h + 2P² curl_a - dP²cross_a) */
            double chiral = 0;
            if (KH != 0) {
                double t1 = 2.0 * P * dPda * h;
                double t2 = 2.0 * P2 * curl_phi[a];
                double t3;
                if (a == 0)      t3 = -(dP2dz * p1 - dP2dy * p2);
                else if (a == 1) t3 = -(dP2dx * p2 - dP2dz * p0);
                else             t3 = -(dP2dy * p0 - dP2dx * p1);
                chiral = KH * (t1 + t2 + t3);
            }

            /* Cosserat strain force on phi: -α × curl(M)_a
             * (Maxima-verified: EL derivative gives MINUS curl(M)) */
            double cosserat_phi = -ACS * curl_M[a];

            /* Curl-squared hardening force on phi: -2 × curl(Q)_a
             * (Maxima-verified: F = -2 curl(Q) where Q=(β/2)|θ|²curl(φ)) */
            double harden_phi = -2.0 * curl_Q[a];

            g->phi_acc[a][idx] = lap - me2*g->phi[a][idx] - dVdP*dPda
                               - t2c*g->phi[a][idx] + ETA*ct + chiral
                               + cosserat_phi + harden_phi;
        }

        /* --- Theta forces --- */
        for (int a = 0; a < NFIELDS; a++) {
            double lapt = (g->theta[a][n_ip]+g->theta[a][n_im]+g->theta[a][n_jp]
                         +g->theta[a][n_jm]+g->theta[a][n_kp]+g->theta[a][n_km]
                         -6.0*g->theta[a][idx]) * idx2;
            double cp = curl_component(g->phi, a, n_ip,n_im,n_jp,n_jm,n_kp,n_km, idx1);

            /* Cosserat strain force on theta: +2α × M_a = 2α(curl(φ)_a/2 - θ_a)
             * = α curl(φ)_a - 2α θ_a
             * This pulls theta toward the geometrically correct value curl(φ)/2.
             * At equilibrium (θ = curl(φ)/2): force is zero. */
            double cosserat_theta = 0;
            if (ACS != 0)
                cosserat_theta = 2.0 * ACS * g->mismatch[a][idx];

            /* Curl-squared hardening force on theta: -β|∇×φ|²θ_a
             * Theta becomes heavy where twist is strong → confines to shell.
             * Massless in vacuum (|∇×φ|=0). */
            double harden_theta = -BH * curl_sq * g->theta[a][idx];

            g->theta_acc[a][idx] = lapt - MTHETA2*g->theta[a][idx] + ETA*cp
                                 + cosserat_theta + harden_theta;
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
   Diagnostics
   ================================================================ */

static void compute_energy(Grid *g, const Config *c,
    double *epk, double *etk, double *eg, double *em, double *ep,
    double *etg, double *etm, double *ec, double *et,
    double *phi_max, double *P_max) {
    const int N=g->N,NN=N*N; const long N3=g->N3;
    const double dx=g->dx, dV=dx*dx*dx, idx1=1.0/(2.0*dx);
    const double MU=c->mu, KAPPA=c->kappa, MASS2=c->m2, MTHETA2=c->mtheta2, ETA=c->eta;
    const int MODE=c->mode; const double KG=c->kappa_gamma;
    double s_epk=0,s_etk=0,s_eg=0,s_em=0,s_ep=0,s_etg=0,s_etm=0,s_ec=0,s_pm=0,s_Pm=0;

    #pragma omp parallel for reduction(+:s_epk,s_etk,s_eg,s_em,s_ep,s_etg,s_etm,s_ec) \
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
            s_etm+=0.5*MTHETA2*g->theta[a][idx]*g->theta[a][idx]*dV;
            double ap=fabs(g->phi[a][idx]); if(ap>s_pm) s_pm=ap;
        }
        double P=p0*p1*p2, P2=P*P;
        if (MODE==3) {
            double D=1.0+KG*sig+KAPPA*P2;
            s_ep+=(MU/2.0)*P2*(1.0+KG*sig)/D*dV;
        } else {
            s_ep+=(MU/2.0)*P2/(1.0+keff*P2)*dV;
        }
        double Pa=fabs(P); if(Pa>s_Pm) s_Pm=Pa;
        /* Cosserat strain energy: α|curl(φ)/2 - θ|² */
        if (c->alpha_cs != 0) {
            double cx0=(g->phi[2][n_jp]-g->phi[2][n_jm]-g->phi[1][n_kp]+g->phi[1][n_km])*idx1;
            double cx1=(g->phi[0][n_kp]-g->phi[0][n_km]-g->phi[2][n_ip]+g->phi[2][n_im])*idx1;
            double cx2=(g->phi[1][n_ip]-g->phi[1][n_im]-g->phi[0][n_jp]+g->phi[0][n_jm])*idx1;
            double m0=cx0*0.5-g->theta[0][idx], m1=cx1*0.5-g->theta[1][idx], m2c=cx2*0.5-g->theta[2][idx];
            s_ep += c->alpha_cs * (m0*m0 + m1*m1 + m2c*m2c) * dV;
        }
        /* Curl-squared hardening energy: (β/2)|θ|²|∇×φ|² */
        if (c->beta_h != 0) {
            double cx0=(g->phi[2][n_jp]-g->phi[2][n_jm]-g->phi[1][n_kp]+g->phi[1][n_km])*idx1;
            double cx1=(g->phi[0][n_kp]-g->phi[0][n_km]-g->phi[2][n_ip]+g->phi[2][n_im])*idx1;
            double cx2=(g->phi[1][n_ip]-g->phi[1][n_im]-g->phi[0][n_jp]+g->phi[0][n_jm])*idx1;
            double csq = cx0*cx0 + cx1*cx1 + cx2*cx2;
            double T2 = g->theta[0][idx]*g->theta[0][idx]+g->theta[1][idx]*g->theta[1][idx]+g->theta[2][idx]*g->theta[2][idx];
            s_ep += 0.5 * c->beta_h * T2 * csq * dV;
        }
        /* Chiral helicity energy: E_h = κ_h P² φ·curl(φ) */
        if (c->kappa_h != 0) {
            double cx0=(g->phi[2][n_jp]-g->phi[2][n_jm]-g->phi[1][n_kp]+g->phi[1][n_km])*idx1;
            double cx1=(g->phi[0][n_kp]-g->phi[0][n_km]-g->phi[2][n_ip]+g->phi[2][n_im])*idx1;
            double cx2=(g->phi[1][n_ip]-g->phi[1][n_im]-g->phi[0][n_jp]+g->phi[0][n_jm])*idx1;
            double hel = p0*cx0 + p1*cx1 + p2*cx2;
            s_ep += c->kappa_h * P2 * hel * dV;
        }
        for (int a=0;a<NFIELDS;a++) {
            double ct=curl_component(g->theta,a,n_ip,n_im,n_jp,n_jm,n_kp,n_km,idx1);
            s_ec-=ETA*g->phi[a][idx]*ct*dV;
        }
    }
    *epk=s_epk;*etk=s_etk;*eg=s_eg;*em=s_em;*ep=s_ep;
    *etg=s_etg;*etm=s_etm;*ec=s_ec;
    *et=s_epk+s_etk+s_eg+s_em+s_ep+s_etg+s_etm+s_ec;
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

/* ================================================================
   SFA output: 12 columns, configurable precision
   ================================================================ */

/* ---- Vector frame fitting ---- */
#define VS_ORDER 3
#define VS_NC    ((VS_ORDER+1)*(VS_ORDER+1)*(VS_ORDER+1))  /* 64 */

/* Fit a tricubic polynomial to an 8³ block via least-squares (normal equations).
 * Writes VS_NC coefficients to out_coeffs. */
static void fit_patch(const float *field, int N, int ox, int oy, int oz,
                      int bs, float *out_coeffs) {
    int o1 = VS_ORDER + 1;  /* 4 */
    int nc = o1 * o1 * o1;  /* 64 */
    /* Use AᵀA solve: accumulate in double for precision */
    double AtA[64*64] = {0};
    double Atb[64] = {0};

    for (int di = 0; di < bs; di++) {
        int gi = ox + di; if (gi >= N) gi = N-1;
        double tx = (bs > 1) ? (double)di / (bs - 1) : 0;
        double sx = 2*tx - 1;  /* Chebyshev: map [0,1] -> [-1,1] */
        double txi[4] = {1, sx, 2*sx*sx-1, 4*sx*sx*sx-3*sx};
        for (int dj = 0; dj < bs; dj++) {
            int gj = oy + dj; if (gj >= N) gj = N-1;
            double ty = (bs > 1) ? (double)dj / (bs - 1) : 0;
            double sy = 2*ty - 1;
            double tyj[4] = {1, sy, 2*sy*sy-1, 4*sy*sy*sy-3*sy};
            for (int dk = 0; dk < bs; dk++) {
                int gk = oz + dk; if (gk >= N) gk = N-1;
                double tz = (bs > 1) ? (double)dk / (bs - 1) : 0;
                double sz = 2*tz - 1;
                double tzk[4] = {1, sz, 2*sz*sz-1, 4*sz*sz*sz-3*sz};

                /* Compute basis values */
                double basis[64];
                int c = 0;
                for (int a = 0; a < o1; a++)
                for (int b = 0; b < o1; b++)
                for (int g2 = 0; g2 < o1; g2++)
                    basis[c++] = txi[a] * tyj[b] * tzk[g2];

                float val = field[(long)gi*N*N + (long)gj*N + gk];
                for (int r = 0; r < nc; r++) {
                    Atb[r] += basis[r] * val;
                    for (int s = r; s < nc; s++)
                        AtA[r*nc+s] += basis[r] * basis[s];
                }
            }
        }
    }
    /* Symmetrize */
    for (int r = 0; r < nc; r++)
        for (int s = 0; s < r; s++)
            AtA[r*nc+s] = AtA[s*nc+r];

    /* Solve via Cholesky (AtA is positive definite for non-degenerate data) */
    /* Simple in-place Cholesky: L * L^T = AtA, then solve L*y=Atb, L^T*x=y */
    double L[64*64] = {0};
    for (int i = 0; i < nc; i++) {
        for (int j = 0; j <= i; j++) {
            double s = AtA[i*nc+j];
            for (int k = 0; k < j; k++) s -= L[i*nc+k] * L[j*nc+k];
            if (i == j) L[i*nc+j] = sqrt(s > 0 ? s : 1e-30);
            else        L[i*nc+j] = s / (L[j*nc+j] + 1e-30);
        }
    }
    /* Forward: L*y = Atb */
    double y[64];
    for (int i = 0; i < nc; i++) {
        double s = Atb[i];
        for (int k = 0; k < i; k++) s -= L[i*nc+k] * y[k];
        y[i] = s / (L[i*nc+i] + 1e-30);
    }
    /* Back: L^T*x = y */
    double x[64];
    for (int i = nc-1; i >= 0; i--) {
        double s = y[i];
        for (int k = i+1; k < nc; k++) s -= L[k*nc+i] * x[k];
        x[i] = s / (L[i*nc+i] + 1e-30);
    }
    for (int i = 0; i < nc; i++) out_coeffs[i] = (float)x[i];
}

/* ---- Voxel snapshot ---- */

static void *cast_buf = NULL;

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
            printf("  WARNING: no KVMD metadata found — using defaults\n");
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

    sfa_embed_kvmd(sfa, &c);
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

    /* Vector frame output (FRVD into same SFA file) */
    int vec_snap_every = 0;
    int vec_frame = 0;
    int vec_enabled = (c.vec_snap_dt > 0);
    int vec_BS = c.vec_block_size > 0 ? c.vec_block_size : 8;
    int vec_BN = c.N / vec_BS;
    int vec_n_patches = vec_BN * vec_BN * vec_BN;
    int vec_nc = (VS_ORDER+1)*(VS_ORDER+1)*(VS_ORDER+1);  /* coeffs per patch per field */
    int vec_nc_total = 6 * vec_nc;  /* all 6 fields per patch */
    int16_t *vec_origins = NULL;
    float *vec_coeffs = NULL;
    float *vec_prev_coeffs = NULL;
    float *vec_buf = NULL;
    /* Temporal model: mean + amp*cos(omega*t + phase) per coefficient */
    float vec_omega = 2.0f * 3.14159265f / 2.2f;  /* breathing frequency */
    float *vec_temp_mean = NULL, *vec_temp_amp = NULL, *vec_temp_phase = NULL;
    float *vec_sum_cos = NULL, *vec_sum_sin = NULL, *vec_sum_mean = NULL;
    int vec_temporal_count = 0;
    int vec_refit_interval = 50;

    if (vec_enabled) {
        vec_snap_every = (int)lround(c.vec_snap_dt / g->dt);
        if (vec_snap_every < 1) vec_snap_every = 1;

        /* Pre-compute patch origins */
        vec_origins = (int16_t*)malloc(vec_n_patches * 3 * sizeof(int16_t));
        int pi = 0;
        for (int bi = 0; bi < vec_BN; bi++)
        for (int bj = 0; bj < vec_BN; bj++)
        for (int bk = 0; bk < vec_BN; bk++) {
            vec_origins[pi*3+0] = (int16_t)(bi * vec_BS);
            vec_origins[pi*3+1] = (int16_t)(bj * vec_BS);
            vec_origins[pi*3+2] = (int16_t)(bk * vec_BS);
            pi++;
        }
        vec_coeffs = (float*)calloc((long)vec_n_patches * vec_nc_total, sizeof(float));
        vec_prev_coeffs = (float*)calloc((long)vec_n_patches * vec_nc_total, sizeof(float));
        vec_buf = (float*)malloc(g->N3 * sizeof(float));
        long vec_n_total = (long)vec_n_patches * vec_nc_total;
        vec_temp_mean  = (float*)calloc(vec_n_total, sizeof(float));
        vec_temp_amp   = (float*)calloc(vec_n_total, sizeof(float));
        vec_temp_phase = (float*)calloc(vec_n_total, sizeof(float));
        vec_sum_cos    = (float*)calloc(vec_n_total, sizeof(float));
        vec_sum_sin    = (float*)calloc(vec_n_total, sizeof(float));
        vec_sum_mean   = (float*)calloc(vec_n_total, sizeof(float));
        printf("Vec frames: BS=%d, %d patches, %d coeffs/patch, snap every %d steps (temporal)\n",
               vec_BS, vec_n_patches, vec_nc_total, vec_snap_every);
    }

    /* Diagnostics file */
    FILE *fp = fopen(c.diag_file, "w");
    fprintf(fp, "t\tE_phi_kin\tE_theta_kin\tE_grad\tE_mass\tE_pot\tE_tgrad\tE_tmass\t"
                "E_coupling\tE_total\tphi_max\tP_max\tP_int\ttheta_rms\n");

    int n_steps = (int)lround(c.T / g->dt);
    int diag_every = (int)lround(c.diag_dt / g->dt); if (diag_every<1) diag_every=1;
    int snap_every = (int)lround(c.snap_dt / g->dt); if (snap_every<1) snap_every=1;
    double burst_dur = c.burst_end - c.burst_start;
    if (c.burst_start > 0 && burst_dur > 0) {
        if (c.burst_every > 0)
            printf("Burst snap: %.1f duration every %.0f time units starting t=%.1f\n",
                   burst_dur, c.burst_every, c.burst_start);
        else
            printf("Burst snap: every timestep from t=%.1f to t=%.1f\n",
                   c.burst_start, c.burst_end);
    }

    /* Initial diagnostic */
    double epk,etk,eg,em,ep,etg,etm,ec,et,pm,Pm;
    compute_energy(g,&c,&epk,&etk,&eg,&em,&ep,&etg,&etm,&ec,&et,&pm,&Pm);
    double Pint0=P_integrated(g), trms0=theta_rms(g);
    double E0 = et;
    printf("INIT: E_total=%.4e E_pot=%.4f phi_max=%.4f P_int=%.4e theta_rms=%.3e\n\n",
           et, ep, pm, Pint0, trms0);
    fprintf(fp,"%.2f\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\n",
            0.0,epk,etk,eg,em,ep,etg,etm,ec,et,pm,Pm,Pint0,trms0);

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

        /* Snap: normal interval OR every step during burst window(s) */
        int in_burst = 0;
        if (burst_dur > 0) {
            double t_rel = t - c.burst_start;
            if (t_rel >= 0) {
                if (c.burst_every > 0)
                    t_rel = fmod(t_rel, c.burst_every);
                in_burst = (t_rel >= 0 && t_rel < burst_dur);
            }
        }
        if (in_burst || step % snap_every == 0)
            sfa_snap(sfa, g, t, c.precision);

        /* Vector frame output (FRVD into same SFA) */
        if (vec_enabled && step % vec_snap_every == 0) {
            int is_iframe = (vec_frame == 0) || (vec_frame % c.vec_iframe_interval == 0);

            /* Fit all 6 fields into one coefficient array: phi[3] + theta[3]
             * Each patch gets 6 * vec_nc coefficients, interleaved per patch */
            int nc_total = vec_nc_total;
            float *all_coeffs = (float*)malloc((long)vec_n_patches * nc_total * sizeof(float));

            for (int f = 0; f < 6; f++) {
                double *src = (f < 3) ? g->phi[f] : g->theta[f-3];
                #pragma omp parallel for schedule(static)
                for (long i = 0; i < g->N3; i++) vec_buf[i] = (float)src[i];

                #pragma omp parallel for schedule(dynamic)
                for (int pi = 0; pi < vec_n_patches; pi++) {
                    fit_patch(vec_buf, g->N,
                              vec_origins[pi*3], vec_origins[pi*3+1], vec_origins[pi*3+2],
                              vec_BS, &all_coeffs[(long)pi * nc_total + f * vec_nc]);
                }
            }

            /* Update temporal accumulators */
            {
                long vnt = (long)vec_n_patches * nc_total;
                float cos_wt = cosf(vec_omega * (float)t);
                float sin_wt = sinf(vec_omega * (float)t);
                for (long i = 0; i < vnt; i++) {
                    vec_sum_mean[i] += all_coeffs[i];
                    vec_sum_cos[i]  += all_coeffs[i] * cos_wt;
                    vec_sum_sin[i]  += all_coeffs[i] * sin_wt;
                }
                vec_temporal_count++;

                /* Bootstrap: on first frame, set mean = actual so P-frames
                 * predict the real initial state instead of zero */
                if (vec_frame == 0) {
                    for (long i = 0; i < vnt; i++) {
                        vec_temp_mean[i] = all_coeffs[i];
                        vec_temp_amp[i] = 0;
                        vec_temp_phase[i] = 0;
                    }
                }

                /* Refit at every I-frame boundary (not just every N frames) */
                if ((is_iframe || vec_temporal_count >= vec_refit_interval) && vec_temporal_count > 2) {
                    float inv_n = 1.0f / vec_temporal_count;
                    for (long i = 0; i < vnt; i++) {
                        vec_temp_mean[i] = vec_sum_mean[i] * inv_n;
                        float sc = vec_sum_cos[i] * inv_n - vec_temp_mean[i] * cos_wt;
                        float ss = vec_sum_sin[i] * inv_n - vec_temp_mean[i] * sin_wt;
                        vec_temp_amp[i] = 2.0f * sqrtf(sc*sc + ss*ss);
                        vec_temp_phase[i] = atan2f(-ss, sc);
                    }
                    memset(vec_sum_cos, 0, vnt * sizeof(float));
                    memset(vec_sum_sin, 0, vnt * sizeof(float));
                    memset(vec_sum_mean, 0, vnt * sizeof(float));
                    vec_temporal_count = 0;
                }
            }

            if (is_iframe) {
                long vnt = (long)vec_n_patches * nc_total;
                sfa_write_vec_iframe_temporal(sfa, t, vec_n_patches, vec_BS, nc_total,
                                     vec_origins, all_coeffs,
                                     vec_omega, vec_temp_mean, vec_temp_amp, vec_temp_phase);
                memcpy(vec_prev_coeffs, all_coeffs, vnt * sizeof(float));
            } else {
                /* P-frame: delta = actual - predicted(t) */
                long vnt = (long)vec_n_patches * nc_total;
                uint32_t *didx = (uint32_t*)malloc(vec_n_patches * sizeof(uint32_t));
                float *dcoeffs = (float*)malloc(vnt * sizeof(float));
                int nd = 0;
                for (int p = 0; p < vec_n_patches; p++) {
                    float mx = 0;
                    long base = (long)p * nc_total;
                    for (int ci = 0; ci < nc_total; ci++) {
                        float predicted = vec_temp_mean[base+ci]
                            + vec_temp_amp[base+ci] * cosf(vec_omega * (float)t + vec_temp_phase[base+ci]);
                        float d = all_coeffs[base+ci] - predicted;
                        dcoeffs[(long)nd*nc_total+ci] = d;  /* tentative */
                        if (fabsf(d) > mx) mx = fabsf(d);
                    }
                    if (mx > c.vec_delta_tol) {
                        didx[nd] = p;
                        /* dcoeffs already written above at nd offset */
                        nd++;
                    }
                }
                sfa_write_vec_pframe(sfa, t, nd, nc_total, didx, dcoeffs);
                free(didx); free(dcoeffs);
                /* Update prev (for fallback if reader lacks temporal) */
                memcpy(vec_prev_coeffs, all_coeffs, vnt * sizeof(float));
            }
            free(all_coeffs);
            vec_frame++;
        }

        if (step % diag_every == 0) {
            compute_energy(g,&c,&epk,&etk,&eg,&em,&ep,&etg,&etm,&ec,&et,&pm,&Pm);
            double Pint=P_integrated(g), trms=theta_rms(g);
            fprintf(fp,"%.2f\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\n",
                    t,epk,etk,eg,em,ep,etg,etm,ec,et,pm,Pm,Pint,trms);
            fflush(fp);

            if (step % major == 0) {
                double wall=omp_get_wtime()-wall0;
                double drift=100*(et-E0)/(fabs(E0)+1e-30);
                printf("t=%7.1f E=%.3e (drift %+.3f%%) Ep=%.1f phi=%.3f θ_rms=%.2e "
                       "[%.0f%% %.1fs %.2fms/step]\n",
                       t,et,drift,ep,pm,trms,100.0*step/n_steps,wall,1000*wall/step);
                fflush(stdout);
            }
        }
    }

    /* Final frame — only if the last step wasn't already (or nearly) a snap point.
     * Skip if the gap is less than half a snap interval to avoid near-duplicate frames. */
    {
        int last_snapped = (n_steps / snap_every) * snap_every;
        int gap = n_steps - last_snapped;
        if (gap > snap_every / 2)
            sfa_snap(sfa, g, n_steps*g->dt, c.precision);
    }
    uint32_t nf = sfa->total_frames;
    sfa_close(sfa);

    /* Cleanup vec frame buffers */
    if (vec_enabled) {
        free(vec_origins); free(vec_coeffs); free(vec_prev_coeffs); free(vec_buf);
        free(vec_temp_mean); free(vec_temp_amp); free(vec_temp_phase);
        free(vec_sum_cos); free(vec_sum_sin); free(vec_sum_mean);
        printf("Vec frames: %d written to %s (temporal)\n", vec_frame, c.output);
    }

    /* Final summary */
    compute_energy(g,&c,&epk,&etk,&eg,&em,&ep,&etg,&etm,&ec,&et,&pm,&Pm);
    double trms=theta_rms(g), Pint=P_integrated(g);
    double wall=omp_get_wtime()-wall0;

    printf("\n=== COMPLETE ===\n");
    printf("E_total=%.4e (drift %.3f%%) E_pot=%.4f\n", et, 100*(et-E0)/(fabs(E0)+1e-30), ep);
    printf("phi_max=%.4f P_int=%.4e theta_rms=%.3e\n", pm, Pint, trms);
    printf("SFA: %s (%u frames)\n", c.output, nf);
    printf("[%s] theta_rms grew: %.2e -> %.2e\n", (trms>trms0+1e-10)?"OK":"WARN", trms0, trms);
    printf("Wall: %.1fs (%.1f min) %.2fms/step\n", wall, wall/60, 1000*wall/n_steps);

    fclose(fp);
    if (cast_buf) free(cast_buf);
    grid_free(g);
    return 0;
}
