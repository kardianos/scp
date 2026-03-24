/*  scp_multi.c — Multi-resolution nested grid Cosserat simulation
 *
 *  Extends scp_sim with Berger-Oliger subcycling for nested grids:
 *  - Root domain (coarsest) with periodic or absorbing BC
 *  - Child domains (finer) with ghost zones filled by trilinear prolongation
 *  - Restriction (fine → coarse) via volume averaging
 *  - Temporal interpolation at coarse/fine boundaries during subcycling
 *
 *  Build: gcc -O3 -march=native -fopenmp -o scp_multi scp_multi.c -lzstd -lm
 *  Run:   OMP_NUM_THREADS=8 ./scp_multi config.cfg
 */

#define SFA_IMPLEMENTATION
#include "../format/sfa.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>
#include <sys/stat.h>
#include <unistd.h>

#define NFIELDS 3
#define PI 3.14159265358979323846
#define MAX_DOMAINS 16

/* ================================================================
   Configuration
   ================================================================ */

typedef struct {
    int N;
    double L, cx, cy, cz;
    int parent_id;      /* -1 for root */
    int ghost;
    double snap_dt;
    /* init */
    char init_mode[32]; /* "braid", "oscillon", "empty" */
    double A, sigma, R_tube, ellip, A_bg;
    double delta[3];
} DomainConfig;

typedef struct {
    /* Physics (shared) */
    double m2, mtheta2, eta, mu, kappa;
    int mode;
    double inv_alpha, inv_beta, kappa_gamma;
    double damp_width, damp_rate;
    double T, dt_factor;
    /* Domains */
    int n_domains;
    DomainConfig dom[MAX_DOMAINS];
    /* Output */
    char output[512];
    char diag_file[512];
    int precision;
    double diag_dt;
} MultiConfig;

static MultiConfig mcfg_defaults(void) {
    MultiConfig c = {0};
    c.m2 = 2.25; c.mtheta2 = 0; c.eta = 0.5;
    c.mu = -41.345; c.kappa = 50.0;
    c.mode = 0; c.inv_alpha = 2.25; c.inv_beta = 5.0; c.kappa_gamma = 2.0;
    c.damp_width = 3.0; c.damp_rate = 0.01;
    c.T = 200; c.dt_factor = 0.025;
    strcpy(c.output, "multi_output.sfa");
    strcpy(c.diag_file, "multi_diag.tsv");
    c.precision = 1; c.diag_dt = 2.0;
    return c;
}

static void mcfg_set(MultiConfig *c, const char *key, const char *val) {
    /* Physics keys */
    if      (!strcmp(key,"m"))           { double m=atof(val); c->m2=m*m; }
    else if (!strcmp(key,"m_theta"))     { double m=atof(val); c->mtheta2=m*m; }
    else if (!strcmp(key,"eta"))         c->eta = atof(val);
    else if (!strcmp(key,"mu"))          c->mu = atof(val);
    else if (!strcmp(key,"kappa"))       c->kappa = atof(val);
    else if (!strcmp(key,"mode"))        c->mode = atoi(val);
    else if (!strcmp(key,"inv_alpha"))   c->inv_alpha = atof(val);
    else if (!strcmp(key,"inv_beta"))    c->inv_beta = atof(val);
    else if (!strcmp(key,"kappa_gamma")) c->kappa_gamma = atof(val);
    else if (!strcmp(key,"damp_width"))  c->damp_width = atof(val);
    else if (!strcmp(key,"damp_rate"))   c->damp_rate = atof(val);
    else if (!strcmp(key,"T"))           c->T = atof(val);
    else if (!strcmp(key,"dt_factor"))   c->dt_factor = atof(val);
    else if (!strcmp(key,"output"))      strncpy(c->output, val, 511);
    else if (!strcmp(key,"diag_file"))   strncpy(c->diag_file, val, 511);
    else if (!strcmp(key,"diag_dt"))     c->diag_dt = atof(val);
    else if (!strcmp(key,"n_grids"))     c->n_domains = atoi(val);
    else if (!strcmp(key,"precision")) {
        if (!strcmp(val,"f16")) c->precision=0;
        else if (!strcmp(val,"f32")) c->precision=1;
        else if (!strcmp(val,"f64")) c->precision=2;
    }
    /* Per-domain keys: grid.X.Y */
    else if (!strncmp(key, "grid.", 5)) {
        int gid = atoi(key + 5);
        if (gid < 0 || gid >= MAX_DOMAINS) return;
        DomainConfig *d = &c->dom[gid];
        const char *sub = strchr(key + 5, '.');
        if (!sub) return; sub++;
        if      (!strcmp(sub,"N"))       d->N = atoi(val);
        else if (!strcmp(sub,"L"))       d->L = atof(val);
        else if (!strcmp(sub,"center"))  sscanf(val,"%lf,%lf,%lf",&d->cx,&d->cy,&d->cz);
        else if (!strcmp(sub,"parent"))  d->parent_id = atoi(val);
        else if (!strcmp(sub,"ghost"))   d->ghost = atoi(val);
        else if (!strcmp(sub,"snap_dt")) d->snap_dt = atof(val);
        else if (!strcmp(sub,"init"))    strncpy(d->init_mode, val, 31);
        else if (!strcmp(sub,"A"))       d->A = atof(val);
        else if (!strcmp(sub,"sigma"))   d->sigma = atof(val);
        else if (!strcmp(sub,"R_tube"))  d->R_tube = atof(val);
        else if (!strcmp(sub,"ellip"))   d->ellip = atof(val);
        else if (!strcmp(sub,"A_bg"))    d->A_bg = atof(val);
        else if (!strcmp(sub,"delta"))   sscanf(val,"%lf,%lf,%lf",&d->delta[0],&d->delta[1],&d->delta[2]);
    }
}

static void mcfg_load(MultiConfig *c, const char *path) {
    FILE *fp = fopen(path, "r");
    if (!fp) { fprintf(stderr, "Cannot open config: %s\n", path); exit(1); }
    char line[2048];
    while (fgets(line, sizeof(line), fp)) {
        char *p = line;
        while (*p==' '||*p=='\t') p++;
        if (*p=='#'||*p=='\n'||*p=='\0') continue;
        char *eq = strchr(p, '=');
        if (!eq) continue;
        *eq = '\0';
        char *key=p, *val=eq+1;
        char *ke=eq-1; while(ke>key&&(*ke==' '||*ke=='\t')) *ke--='\0';
        while(*val==' '||*val=='\t') val++;
        char *ve=val+strlen(val)-1; while(ve>val&&(*ve=='\n'||*ve=='\r'||*ve==' ')) *ve--='\0';
        char *hash=strchr(val,'#'); if(hash){*hash='\0';ve=hash-1;while(ve>val&&*ve==' ')*ve--='\0';}
        mcfg_set(c, key, val);
    }
    fclose(fp);
    /* Fill defaults for domains */
    for (int i = 0; i < c->n_domains; i++) {
        DomainConfig *d = &c->dom[i];
        if (d->parent_id == 0 && i == 0) d->parent_id = -1;  /* root */
        if (d->ghost == 0 && d->parent_id >= 0) d->ghost = 2;
        if (d->snap_dt <= 0) d->snap_dt = 5.0;
        if (d->A <= 0) d->A = 0.8;
        if (d->A_bg <= 0) d->A_bg = 0.1;
        if (d->R_tube <= 0) d->R_tube = 3.0;
        if (d->sigma <= 0) d->sigma = 3.0;
        if (d->delta[1] == 0 && d->delta[2] == 0) {
            d->delta[0]=0; d->delta[1]=3.0005; d->delta[2]=4.4325;
        }
        if (d->init_mode[0] == '\0') strcpy(d->init_mode, "empty");
    }
}

/* ================================================================
   Grid with ghost zones
   ================================================================ */

typedef struct {
    double *mem;
    double *phi[NFIELDS], *phi_vel[NFIELDS], *phi_acc[NFIELDS];
    double *theta[NFIELDS], *theta_vel[NFIELDS], *theta_acc[NFIELDS];
    int N_phys;    /* physical cells per axis */
    int ghost;     /* ghost cells per side */
    int N_total;   /* N_phys + 2*ghost */
    long N3;       /* N_total^3 */
    double L, dx, dt;
} GGrid;

static GGrid *ggrid_alloc(int N_phys, double L, int ghost, double dt_factor) {
    GGrid *g = (GGrid*)calloc(1, sizeof(GGrid));
    g->N_phys = N_phys;
    g->ghost = ghost;
    g->N_total = N_phys + 2 * ghost;
    g->N3 = (long)g->N_total * g->N_total * g->N_total;
    g->L = L;
    g->dx = 2.0 * L / (N_phys - 1);
    g->dt = dt_factor * g->dx;

    long total = 18 * g->N3;
    g->mem = (double*)malloc(total * sizeof(double));
    if (!g->mem) { fprintf(stderr, "FATAL: malloc %ld doubles\n", total); exit(1); }
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

static void ggrid_free(GGrid *g) { free(g->mem); free(g); }

/* Index into GGrid: (i,j,k) are absolute indices in N_total grid */
static inline long gidx(GGrid *g, int i, int j, int k) {
    return (long)i * g->N_total * g->N_total + (long)j * g->N_total + k;
}

/* World coordinate from grid index (in physical frame) */
static inline double gworld(GGrid *g, double center, int idx) {
    /* idx is in N_total frame; physical cell 0 is at index ghost */
    return center - g->L + (idx - g->ghost) * g->dx;
}

/* ================================================================
   Domain
   ================================================================ */

typedef struct Domain {
    int id;
    GGrid *grid;
    double cx, cy, cz;        /* center in world coords */
    int parent_idx;            /* index into domain array, -1 for root */
    int substeps;              /* dt_parent / dt_self */
    double snap_dt;
    double next_snap;

    /* Temporal interpolation buffers (parent overlap region) */
    double *pbuf_t0[12];      /* parent state at t_n */
    double *pbuf_t1[12];      /* parent state at t_{n+1} */
    int pb_i0, pb_j0, pb_k0;  /* start indices in parent grid */
    int pb_ni, pb_nj, pb_nk;  /* size of parent buffer subgrid */
} Domain;

/* ================================================================
   Curl helper (same as scp_sim.c)
   ================================================================ */

static inline double curl_comp(double *F[3], int a,
    long n_ip, long n_im, long n_jp, long n_jm, long n_kp, long n_km, double idx1) {
    if (a==0) return (F[2][n_jp]-F[2][n_jm]-F[1][n_kp]+F[1][n_km])*idx1;
    if (a==1) return (F[0][n_kp]-F[0][n_km]-F[2][n_ip]+F[2][n_im])*idx1;
    return (F[1][n_ip]-F[1][n_im]-F[0][n_jp]+F[0][n_jm])*idx1;
}

/* ================================================================
   Force computation — with ghost zones (no wrapping)
   Forces only computed for physical cells [ghost .. ghost+N_phys-1]
   Ghost cells provide boundary stencil data.
   ================================================================ */

static void compute_forces_ghost(GGrid *g, const MultiConfig *c) {
    const int Nt = g->N_total, NN = Nt*Nt;
    const int g0 = g->ghost, g1 = g->ghost + g->N_phys;
    const double idx2 = 1.0/(g->dx*g->dx), idx1 = 1.0/(2.0*g->dx);
    const double MU=c->mu, KAPPA=c->kappa, MASS2=c->m2;
    const double MTHETA2=c->mtheta2, ETA=c->eta;
    const int MODE=c->mode;
    const double KG=c->kappa_gamma, IA=c->inv_alpha, IB=c->inv_beta;

    #pragma omp parallel for schedule(static) collapse(3)
    for (int i = g0; i < g1; i++)
    for (int j = g0; j < g1; j++)
    for (int k = g0; k < g1; k++) {
        long idx = (long)i*NN + j*Nt + k;
        /* Neighbors: no wrapping, ghost zones provide boundary data */
        long n_ip = (long)(i+1)*NN+j*Nt+k, n_im = (long)(i-1)*NN+j*Nt+k;
        long n_jp = (long)i*NN+(j+1)*Nt+k, n_jm = (long)i*NN+(j-1)*Nt+k;
        long n_kp = (long)i*NN+j*Nt+(k+1), n_km = (long)i*NN+j*Nt+(k-1);

        double p0=g->phi[0][idx], p1=g->phi[1][idx], p2=g->phi[2][idx];
        double sig = p0*p0+p1*p1+p2*p2;
        double me2 = (MODE==1)?IA/(1.0+IB*sig):MASS2;
        double keff = (MODE==3)?KAPPA/(1.0+KG*sig):KAPPA;

        double P=p0*p1*p2, P2=P*P;
        double den=1.0+keff*P2;
        double dVdP = MU*P/(den*den);
        double t2c = 0;
        if (MODE==3 && KG>0) { double D=1.0+KG*sig+KAPPA*P2; t2c=MU*KG*KAPPA*P2*P2/(D*D); }

        for (int a=0;a<NFIELDS;a++) {
            double lap = (g->phi[a][n_ip]+g->phi[a][n_im]+g->phi[a][n_jp]
                        +g->phi[a][n_jm]+g->phi[a][n_kp]+g->phi[a][n_km]
                        -6.0*g->phi[a][idx]) * idx2;
            double dPda = (a==0)?p1*p2:(a==1)?p0*p2:p0*p1;
            double ct = curl_comp(g->theta, a, n_ip,n_im,n_jp,n_jm,n_kp,n_km, idx1);
            g->phi_acc[a][idx] = lap - me2*g->phi[a][idx] - dVdP*dPda - t2c*g->phi[a][idx] + ETA*ct;
        }
        for (int a=0;a<NFIELDS;a++) {
            double lapt = (g->theta[a][n_ip]+g->theta[a][n_im]+g->theta[a][n_jp]
                         +g->theta[a][n_jm]+g->theta[a][n_kp]+g->theta[a][n_km]
                         -6.0*g->theta[a][idx]) * idx2;
            double cp = curl_comp(g->phi, a, n_ip,n_im,n_jp,n_jm,n_kp,n_km, idx1);
            g->theta_acc[a][idx] = lapt - MTHETA2*g->theta[a][idx] + ETA*cp;
        }
    }
}

/* Root domain: periodic BC force computation (same as scp_sim.c) */
static void compute_forces_periodic(GGrid *g, const MultiConfig *c) {
    const int N=g->N_phys, NN=N*N;  /* root has ghost=0, N_total=N_phys */
    const long N3 = (long)N*N*N;
    const double idx2=1.0/(g->dx*g->dx), idx1=1.0/(2.0*g->dx);
    const double MU=c->mu, KAPPA=c->kappa, MASS2=c->m2;
    const double MTHETA2=c->mtheta2, ETA=c->eta;
    const int MODE=c->mode;
    const double KG=c->kappa_gamma, IA=c->inv_alpha, IB=c->inv_beta;

    #pragma omp parallel for schedule(static)
    for (long idx=0; idx<N3; idx++) {
        int i=(int)(idx/NN),j=(int)((idx/N)%N),k=(int)(idx%N);
        int ip=(i+1)%N,im=(i-1+N)%N,jp=(j+1)%N,jm=(j-1+N)%N,kp=(k+1)%N,km=(k-1+N)%N;
        long n_ip=(long)ip*NN+j*N+k, n_im=(long)im*NN+j*N+k;
        long n_jp=(long)i*NN+jp*N+k, n_jm=(long)i*NN+jm*N+k;
        long n_kp=(long)i*NN+j*N+kp, n_km=(long)i*NN+j*N+km;

        double p0=g->phi[0][idx],p1=g->phi[1][idx],p2=g->phi[2][idx];
        double sig=p0*p0+p1*p1+p2*p2;
        double me2=(MODE==1)?IA/(1.0+IB*sig):MASS2;
        double keff=(MODE==3)?KAPPA/(1.0+KG*sig):KAPPA;
        double P=p0*p1*p2,P2=P*P;
        double den=1.0+keff*P2;
        double dVdP=MU*P/(den*den);
        double t2c=0;
        if(MODE==3&&KG>0){double D=1.0+KG*sig+KAPPA*P2;t2c=MU*KG*KAPPA*P2*P2/(D*D);}

        for (int a=0;a<NFIELDS;a++) {
            double lap=(g->phi[a][n_ip]+g->phi[a][n_im]+g->phi[a][n_jp]
                       +g->phi[a][n_jm]+g->phi[a][n_kp]+g->phi[a][n_km]
                       -6.0*g->phi[a][idx])*idx2;
            double dPda=(a==0)?p1*p2:(a==1)?p0*p2:p0*p1;
            double ct=curl_comp(g->theta,a,n_ip,n_im,n_jp,n_jm,n_kp,n_km,idx1);
            g->phi_acc[a][idx]=lap-me2*g->phi[a][idx]-dVdP*dPda-t2c*g->phi[a][idx]+ETA*ct;
        }
        for (int a=0;a<NFIELDS;a++) {
            double lapt=(g->theta[a][n_ip]+g->theta[a][n_im]+g->theta[a][n_jp]
                        +g->theta[a][n_jm]+g->theta[a][n_kp]+g->theta[a][n_km]
                        -6.0*g->theta[a][idx])*idx2;
            double cp=curl_comp(g->phi,a,n_ip,n_im,n_jp,n_jm,n_kp,n_km,idx1);
            g->theta_acc[a][idx]=lapt-MTHETA2*g->theta[a][idx]+ETA*cp;
        }
    }
}

/* ================================================================
   Absorbing boundary (root domain only)
   ================================================================ */

static void apply_damping(GGrid *g, const MultiConfig *c) {
    const int N=g->N_phys, NN=N*N;
    const double L=g->L, dx=g->dx, DW=c->damp_width, DR=c->damp_rate;
    if (DW<=0||DR<=0) return;
    const double R_damp = L - DW;
    #pragma omp parallel for schedule(static)
    for (long idx=0; idx<(long)N*N*N; idx++) {
        int i=(int)(idx/NN),j=(int)((idx/N)%N),k=(int)(idx%N);
        double x=-L+i*dx, y=-L+j*dx, z=-L+k*dx;
        double r=sqrt(x*x+y*y+z*z);
        if (r>R_damp) {
            double s=(r-R_damp)/DW; if(s>1)s=1;
            double d=1.0-DR*s*s;
            for(int a=0;a<NFIELDS;a++){g->phi_vel[a][idx]*=d;g->theta_vel[a][idx]*=d;}
        }
    }
}

/* ================================================================
   Split Verlet operations (for orchestrating subcycling)
   ================================================================ */

static void kick_half(GGrid *g, double hdt, int use_ghost) {
    const int Nt=g->N_total, NN=Nt*Nt;
    const int g0=use_ghost?g->ghost:0, g1=use_ghost?(g->ghost+g->N_phys):g->N_phys;
    #pragma omp parallel for schedule(static) collapse(3)
    for (int i=g0;i<g1;i++) for(int j=g0;j<g1;j++) for(int k=g0;k<g1;k++) {
        long idx=(long)i*NN+j*Nt+k;
        for(int a=0;a<NFIELDS;a++){
            g->phi_vel[a][idx]+=hdt*g->phi_acc[a][idx];
            g->theta_vel[a][idx]+=hdt*g->theta_acc[a][idx];
        }
    }
}

static void drift(GGrid *g, double dt, int use_ghost) {
    const int Nt=g->N_total, NN=Nt*Nt;
    const int g0=use_ghost?g->ghost:0, g1=use_ghost?(g->ghost+g->N_phys):g->N_phys;
    #pragma omp parallel for schedule(static) collapse(3)
    for (int i=g0;i<g1;i++) for(int j=g0;j<g1;j++) for(int k=g0;k<g1;k++) {
        long idx=(long)i*NN+j*Nt+k;
        for(int a=0;a<NFIELDS;a++){
            g->phi[a][idx]+=dt*g->phi_vel[a][idx];
            g->theta[a][idx]+=dt*g->theta_vel[a][idx];
        }
    }
}

/* ================================================================
   Prolongation: trilinear interpolation from parent to child ghost zones
   ================================================================ */

/* Sample a single value from the parent buffer via trilinear interpolation.
   (fi, fj, fk) are fractional indices in the parent buffer. */
static inline double trilinear(const double *buf, int ni, int nj, int nk,
                               double fi, double fj, double fk) {
    int i0=(int)floor(fi), j0=(int)floor(fj), k0=(int)floor(fk);
    double di=fi-i0, dj=fj-j0, dk=fk-k0;
    /* Clamp to buffer bounds */
    int i1=i0+1, j1=j0+1, k1=k0+1;
    if(i0<0)i0=0; if(j0<0)j0=0; if(k0<0)k0=0;
    if(i1>=ni)i1=ni-1; if(j1>=nj)j1=nj-1; if(k1>=nk)k1=nk-1;
    if(i0>=ni)i0=ni-1; if(j0>=nj)j0=nj-1; if(k0>=nk)k0=nk-1;

    double c000=buf[(long)i0*nj*nk+j0*nk+k0], c100=buf[(long)i1*nj*nk+j0*nk+k0];
    double c010=buf[(long)i0*nj*nk+j1*nk+k0], c110=buf[(long)i1*nj*nk+j1*nk+k0];
    double c001=buf[(long)i0*nj*nk+j0*nk+k1], c101=buf[(long)i1*nj*nk+j0*nk+k1];
    double c011=buf[(long)i0*nj*nk+j1*nk+k1], c111=buf[(long)i1*nj*nk+j1*nk+k1];

    double c00=c000*(1-di)+c100*di, c01=c001*(1-di)+c101*di;
    double c10=c010*(1-di)+c110*di, c11=c011*(1-di)+c111*di;
    double c0=c00*(1-dj)+c10*dj, c1=c01*(1-dj)+c11*dj;
    return c0*(1-dk)+c1*dk;
}

/* Fill parent buffer from parent grid (subset covering child extent) */
static void fill_parent_buffer(Domain *child, Domain *parent, double **buf) {
    GGrid *pg = parent->grid, *cg = child->grid;
    double child_xlo = child->cx - child->grid->L - child->grid->ghost * child->grid->dx;
    double child_ylo = child->cy - child->grid->L - child->grid->ghost * child->grid->dx;
    double child_zlo = child->cz - child->grid->L - child->grid->ghost * child->grid->dx;

    /* Map child extent to parent grid indices */
    int pi0 = (int)floor((child_xlo - (parent->cx - pg->L)) / pg->dx) - 1;
    int pj0 = (int)floor((child_ylo - (parent->cy - pg->L)) / pg->dx) - 1;
    int pk0 = (int)floor((child_zlo - (parent->cz - pg->L)) / pg->dx) - 1;
    if(pi0<0)pi0=0; if(pj0<0)pj0=0; if(pk0<0)pk0=0;

    int Np = pg->N_phys;
    int pi1 = pi0 + (int)ceil(2.0*cg->L/pg->dx) + 4;
    int pj1 = pj0 + (int)ceil(2.0*cg->L/pg->dx) + 4;
    int pk1 = pk0 + (int)ceil(2.0*cg->L/pg->dx) + 4;
    if(pi1>Np)pi1=Np; if(pj1>Np)pj1=Np; if(pk1>Np)pk1=Np;

    child->pb_i0=pi0; child->pb_j0=pj0; child->pb_k0=pk0;
    child->pb_ni=pi1-pi0; child->pb_nj=pj1-pj0; child->pb_nk=pk1-pk0;
    long pb_n3 = (long)child->pb_ni * child->pb_nj * child->pb_nk;

    double *arrays[12] = {
        pg->phi[0],pg->phi[1],pg->phi[2],pg->theta[0],pg->theta[1],pg->theta[2],
        pg->phi_vel[0],pg->phi_vel[1],pg->phi_vel[2],
        pg->theta_vel[0],pg->theta_vel[1],pg->theta_vel[2]
    };

    for (int f = 0; f < 12; f++) {
        if (!buf[f]) buf[f] = (double*)malloc(pb_n3 * sizeof(double));
        for (int i=0;i<child->pb_ni;i++)
        for (int j=0;j<child->pb_nj;j++)
        for (int k=0;k<child->pb_nk;k++) {
            int pi=pi0+i, pj=pj0+j, pk=pk0+k;
            /* Wrap for periodic parent */
            pi=((pi%Np)+Np)%Np; pj=((pj%Np)+Np)%Np; pk=((pk%Np)+Np)%Np;
            long pidx = (long)pi*Np*Np + pj*Np + pk;
            buf[f][(long)i*child->pb_nj*child->pb_nk + j*child->pb_nk + k] = arrays[f][pidx];
        }
    }
}

/* Prolongate: fill child ghost zones from parent buffers with temporal interpolation */
static void prolongate(Domain *child, Domain *parent, double alpha) {
    GGrid *cg = child->grid, *pg = parent->grid;
    int Nt=cg->N_total, NN=Nt*Nt;
    int g0=cg->ghost, g1=cg->ghost+cg->N_phys;

    double *arrays[12] = {
        cg->phi[0],cg->phi[1],cg->phi[2],cg->theta[0],cg->theta[1],cg->theta[2],
        cg->phi_vel[0],cg->phi_vel[1],cg->phi_vel[2],
        cg->theta_vel[0],cg->theta_vel[1],cg->theta_vel[2]
    };

    /* Only fill ghost cells (outside [g0, g1) in each dimension) */
    for (int i=0;i<Nt;i++) for(int j=0;j<Nt;j++) for(int k=0;k<Nt;k++) {
        int is_ghost = (i<g0||i>=g1||j<g0||j>=g1||k<g0||k>=g1);
        if (!is_ghost) continue;

        long idx = (long)i*NN + j*Nt + k;
        /* World position of this ghost cell */
        double wx = gworld(cg, child->cx, i);
        double wy = gworld(cg, child->cy, j);
        double wz = gworld(cg, child->cz, k);
        /* Fractional index in parent buffer */
        double fi = (wx - (parent->cx - pg->L))/pg->dx - child->pb_i0;
        double fj = (wy - (parent->cy - pg->L))/pg->dx - child->pb_j0;
        double fk = (wz - (parent->cz - pg->L))/pg->dx - child->pb_k0;

        for (int f = 0; f < 12; f++) {
            double v0 = trilinear(child->pbuf_t0[f], child->pb_ni, child->pb_nj, child->pb_nk, fi,fj,fk);
            double v1 = trilinear(child->pbuf_t1[f], child->pb_ni, child->pb_nj, child->pb_nk, fi,fj,fk);
            arrays[f][idx] = (1.0 - alpha)*v0 + alpha*v1;
        }
    }
}

/* ================================================================
   Restriction: volume-average fine interior → coarse overlap cells
   ================================================================ */

static void restrict_to_parent(Domain *child, Domain *parent) {
    GGrid *cg = child->grid, *pg = parent->grid;
    int Np = pg->N_phys;
    double ratio = pg->dx / cg->dx;  /* cells per parent cell */
    int iratio = (int)round(ratio);
    if (iratio < 1) iratio = 1;

    double *p_arrays[12] = {
        pg->phi[0],pg->phi[1],pg->phi[2],pg->theta[0],pg->theta[1],pg->theta[2],
        pg->phi_vel[0],pg->phi_vel[1],pg->phi_vel[2],
        pg->theta_vel[0],pg->theta_vel[1],pg->theta_vel[2]
    };
    double *c_arrays[12] = {
        cg->phi[0],cg->phi[1],cg->phi[2],cg->theta[0],cg->theta[1],cg->theta[2],
        cg->phi_vel[0],cg->phi_vel[1],cg->phi_vel[2],
        cg->theta_vel[0],cg->theta_vel[1],cg->theta_vel[2]
    };

    int cg0 = cg->ghost;
    int Nc = cg->N_phys;
    int Ct = cg->N_total;

    /* For each parent cell overlapping the child's physical region */
    for (int pi=0;pi<Np;pi++) for(int pj=0;pj<Np;pj++) for(int pk=0;pk<Np;pk++) {
        double px = parent->cx - pg->L + pi*pg->dx;
        double py = parent->cy - pg->L + pj*pg->dx;
        double pz = parent->cz - pg->L + pk*pg->dx;

        /* Check if this parent cell is inside the child's physical region */
        if (px < child->cx-cg->L || px > child->cx+cg->L) continue;
        if (py < child->cy-cg->L || py > child->cy+cg->L) continue;
        if (pz < child->cz-cg->L || pz > child->cz+cg->L) continue;

        /* Find child cells covering this parent cell */
        int ci0 = cg0 + (int)round((px - (child->cx-cg->L))/cg->dx);
        int cj0 = cg0 + (int)round((py - (child->cy-cg->L))/cg->dx);
        int ck0 = cg0 + (int)round((pz - (child->cz-cg->L))/cg->dx);

        long pidx = (long)pi*Np*Np + pj*Np + pk;

        for (int f = 0; f < 12; f++) {
            double sum = 0; int cnt = 0;
            for (int di=0;di<iratio;di++) for(int dj=0;dj<iratio;dj++) for(int dk=0;dk<iratio;dk++) {
                int ci=ci0+di, cj=cj0+dj, ck=ck0+dk;
                if (ci>=cg0 && ci<cg0+Nc && cj>=cg0 && cj<cg0+Nc && ck>=cg0 && ck<cg0+Nc) {
                    sum += c_arrays[f][(long)ci*Ct*Ct + cj*Ct + ck];
                    cnt++;
                }
            }
            if (cnt > 0) p_arrays[f][pidx] = sum / cnt;
        }
    }
}

/* ================================================================
   Initialization
   ================================================================ */

static void init_domain(Domain *dom, const DomainConfig *dc, const MultiConfig *c) {
    GGrid *g = dom->grid;
    int Nt=g->N_total, NN=Nt*Nt;
    int g0=g->ghost;
    double kw=PI/g->L, omega=sqrt(kw*kw+c->m2);
    double k_bg=PI/g->L, omega_bg=sqrt(k_bg*k_bg+c->m2);

    if (!strcmp(dc->init_mode, "braid")) {
        double sx=1+dc->ellip, sy=1-dc->ellip;
        double inv2R2=1.0/(2*dc->R_tube*dc->R_tube);
        for(int i=g0;i<g0+g->N_phys;i++) for(int j=g0;j<g0+g->N_phys;j++)
        for(int k=g0;k<g0+g->N_phys;k++) {
            long idx=(long)i*NN+j*Nt+k;
            double x=gworld(g,dom->cx,i), y=gworld(g,dom->cy,j), z=gworld(g,dom->cz,k);
            double r2e=x*x/(sx*sx)+y*y/(sy*sy);
            double env=exp(-r2e*inv2R2);
            for(int a=0;a<NFIELDS;a++){
                double ph=kw*z+dc->delta[a], ph_bg=k_bg*z+2*PI*a/3.0;
                g->phi[a][idx]=dc->A*env*cos(ph)+dc->A_bg*cos(ph_bg);
                g->phi_vel[a][idx]=omega*dc->A*env*sin(ph)+omega_bg*dc->A_bg*sin(ph_bg);
            }
        }
    } else if (!strcmp(dc->init_mode, "oscillon")) {
        for(int i=g0;i<g0+g->N_phys;i++) for(int j=g0;j<g0+g->N_phys;j++)
        for(int k=g0;k<g0+g->N_phys;k++) {
            long idx=(long)i*NN+j*Nt+k;
            double x=gworld(g,dom->cx,i), y=gworld(g,dom->cy,j), z=gworld(g,dom->cz,k);
            double r2=x*x+y*y+z*z;
            double env=dc->A*exp(-r2/(2.0*dc->sigma*dc->sigma));
            for(int a=0;a<NFIELDS;a++) g->phi[a][idx]=env*cos(dc->delta[a]);
        }
    } else {
        /* "empty" — background only */
        for(int i=g0;i<g0+g->N_phys;i++) for(int j=g0;j<g0+g->N_phys;j++)
        for(int k=g0;k<g0+g->N_phys;k++) {
            long idx=(long)i*NN+j*Nt+k;
            double z=gworld(g,dom->cz,k);
            for(int a=0;a<NFIELDS;a++){
                double ph_bg=k_bg*z+2*PI*a/3.0;
                g->phi[a][idx]=dc->A_bg*cos(ph_bg);
                g->phi_vel[a][idx]=omega_bg*dc->A_bg*sin(ph_bg);
            }
        }
    }
}

/* ================================================================
   Berger-Oliger stepping
   ================================================================ */

static void step_multigrid(Domain *domains, int n_dom, const MultiConfig *c) {
    Domain *root = &domains[0];
    GGrid *rg = root->grid;
    double dt_coarse = rg->dt;

    /* 1. Store parent buffers at t_n for each child */
    for (int d = 1; d < n_dom; d++)
        fill_parent_buffer(&domains[d], &domains[domains[d].parent_idx], domains[d].pbuf_t0);

    /* 2. Root: half-kick + drift */
    kick_half(rg, 0.5*dt_coarse, 0);
    drift(rg, dt_coarse, 0);

    /* 3. Root: compute forces */
    compute_forces_periodic(rg, c);

    /* 4. Root: half-kick (second half) */
    kick_half(rg, 0.5*dt_coarse, 0);

    /* 5. Root: absorbing BC */
    apply_damping(rg, c);

    /* 6. Store parent buffers at t_{n+1} */
    for (int d = 1; d < n_dom; d++)
        fill_parent_buffer(&domains[d], &domains[domains[d].parent_idx], domains[d].pbuf_t1);

    /* 7. Subcycle each child domain */
    for (int d = 1; d < n_dom; d++) {
        Domain *child = &domains[d];
        GGrid *cg = child->grid;
        double dt_fine = cg->dt;
        int substeps = child->substeps;

        for (int sub = 0; sub < substeps; sub++) {
            double alpha = ((double)sub + 0.5) / substeps;
            /* Fill ghost zones */
            prolongate(child, &domains[child->parent_idx], alpha);
            /* Verlet step */
            kick_half(cg, 0.5*dt_fine, 1);
            drift(cg, dt_fine, 1);
            compute_forces_ghost(cg, c);
            kick_half(cg, 0.5*dt_fine, 1);
        }

        /* 8. Restrict child back to parent */
        restrict_to_parent(child, &domains[child->parent_idx]);
    }
}

/* ================================================================
   Diagnostics (per-domain energy, computed on physical cells only)
   ================================================================ */

static void domain_energy(GGrid *g, const MultiConfig *c, int is_root,
    double *E_tot, double *E_pot, double *phi_max, double *theta_rms_out) {
    int Nt=g->N_total, NN=Nt*Nt;
    int g0=is_root?0:g->ghost, g1=is_root?g->N_phys:(g->ghost+g->N_phys);
    double dx=g->dx, dV=dx*dx*dx;
    double idx1=1.0/(2.0*dx);
    double epk=0,etk=0,eg=0,em=0,ep=0,etg=0,etm=0,ec=0,pm=0,trms=0;

    for(int i=g0;i<g1;i++) for(int j=g0;j<g1;j++) for(int k=g0;k<g1;k++) {
        long idx=(long)i*NN+j*Nt+k;
        long n_ip=(long)(i+1)*NN+j*Nt+k, n_im=(long)(i-1)*NN+j*Nt+k;
        long n_jp=(long)i*NN+(j+1)*Nt+k, n_jm=(long)i*NN+(j-1)*Nt+k;
        long n_kp=(long)i*NN+j*Nt+(k+1), n_km=(long)i*NN+j*Nt+(k-1);
        /* Wrap for root (ghost=0) */
        if (is_root) {
            int N=g->N_phys;
            int ii=i,jj=j,kk=k;
            n_ip=(long)((ii+1)%N)*N*N+jj*N+kk; n_im=(long)((ii-1+N)%N)*N*N+jj*N+kk;
            n_jp=(long)ii*N*N+((jj+1)%N)*N+kk; n_jm=(long)ii*N*N+((jj-1+N)%N)*N+kk;
            n_kp=(long)ii*N*N+jj*N+((kk+1)%N); n_km=(long)ii*N*N+jj*N+((kk-1+N)%N);
        }

        double p0=g->phi[0][idx],p1=g->phi[1][idx],p2=g->phi[2][idx];
        for(int a=0;a<NFIELDS;a++){
            epk+=0.5*g->phi_vel[a][idx]*g->phi_vel[a][idx]*dV;
            etk+=0.5*g->theta_vel[a][idx]*g->theta_vel[a][idx]*dV;
            double gx=(g->phi[a][n_ip]-g->phi[a][n_im])*idx1;
            double gy=(g->phi[a][n_jp]-g->phi[a][n_jm])*idx1;
            double gz=(g->phi[a][n_kp]-g->phi[a][n_km])*idx1;
            eg+=0.5*(gx*gx+gy*gy+gz*gz)*dV;
            em+=0.5*c->m2*g->phi[a][idx]*g->phi[a][idx]*dV;
            double ap=fabs(g->phi[a][idx]); if(ap>pm) pm=ap;
            trms+=g->theta[a][idx]*g->theta[a][idx];
        }
        double P=p0*p1*p2,P2=P*P;
        ep+=(c->mu/2.0)*P2/(1.0+c->kappa*P2)*dV;
    }
    long nphys3=(long)(g1-g0)*(g1-g0)*(g1-g0);
    *E_tot=epk+etk+eg+em+ep;
    *E_pot=ep; *phi_max=pm;
    *theta_rms_out=sqrt(trms/(3.0*nphys3));
}

/* ================================================================
   f16 helpers
   ================================================================ */

static inline uint16_t f64_to_f16(double v) {
    float f=(float)v; uint32_t x; memcpy(&x,&f,4);
    uint16_t sign=(x>>16)&0x8000; int exp=((x>>23)&0xFF)-127+15; uint16_t mant=(x>>13)&0x3FF;
    if(exp<=0)return sign; if(exp>=31)return sign|0x7C00;
    return sign|(exp<<10)|mant;
}

/* ================================================================
   SFA output: write one domain's frame
   ================================================================ */

static void snap_domain(SFA *sfa, Domain *dom, double t, int precision) {
    GGrid *g = dom->grid;
    int Np=g->N_phys, Nt=g->N_total, g0=g->ghost;
    long Np3 = (long)Np*Np*Np;

    /* Extract physical cells (strip ghost zones) into contiguous buffer */
    double *phys[12];
    double *arrays[12] = {
        g->phi[0],g->phi[1],g->phi[2],g->theta[0],g->theta[1],g->theta[2],
        g->phi_vel[0],g->phi_vel[1],g->phi_vel[2],
        g->theta_vel[0],g->theta_vel[1],g->theta_vel[2]
    };
    for (int f=0;f<12;f++) {
        phys[f] = (double*)malloc(Np3*sizeof(double));
        long out=0;
        for(int i=g0;i<g0+Np;i++) for(int j=g0;j<g0+Np;j++) for(int k=g0;k<g0+Np;k++)
            phys[f][out++] = arrays[f][(long)i*Nt*Nt+j*Nt+k];
    }

    /* Compute frame_bytes for this grid */
    uint64_t bytes_per_point = 0;
    for (uint32_t c=0;c<sfa->n_columns;c++) bytes_per_point+=sfa_dtype_size[sfa->columns[c].dtype];
    uint64_t fb = (uint64_t)Np3 * bytes_per_point;

    if (precision == 2) {
        sfa_write_frame_ex(sfa, t, (void**)phys, (uint32_t)dom->id, fb);
    } else {
        void *cols[12];
        int es = (precision==0)?2:4;
        for (int f=0;f<12;f++) {
            cols[f]=malloc(Np3*es);
            if(precision==1){float *p=(float*)cols[f]; for(long i=0;i<Np3;i++) p[i]=(float)phys[f][i];}
            else{uint16_t *p=(uint16_t*)cols[f]; for(long i=0;i<Np3;i++) p[i]=f64_to_f16(phys[f][i]);}
        }
        sfa_write_frame_ex(sfa, t, cols, (uint32_t)dom->id, fb);
        for(int f=0;f<12;f++) free(cols[f]);
    }
    for(int f=0;f<12;f++) free(phys[f]);
}

/* ================================================================
   Main
   ================================================================ */

int main(int argc, char **argv) {
    MultiConfig c = mcfg_defaults();
    if (argc < 2) {
        fprintf(stderr, "Usage: %s config.cfg [-key value ...]\n", argv[0]);
        return 1;
    }
    mcfg_load(&c, argv[1]);
    for (int i=2;i<argc-1;i+=2) { const char *k=argv[i]; if(k[0]=='-')k++; mcfg_set(&c,k,argv[i+1]); }

    int nth=4; char *et=getenv("OMP_NUM_THREADS"); if(et) nth=atoi(et);
    omp_set_num_threads(nth);

    printf("=== scp_multi: Multi-Resolution Cosserat Kernel ===\n");
    printf("Physics: m²=%.4f η=%.3f μ=%.3f κ=%.1f\n", c.m2, c.eta, c.mu, c.kappa);
    printf("Domains: %d\n\n", c.n_domains);

    /* Allocate domains */
    Domain domains[MAX_DOMAINS] = {0};
    for (int d=0; d<c.n_domains; d++) {
        DomainConfig *dc = &c.dom[d];
        int ghost = (dc->parent_id < 0) ? 0 : dc->ghost;
        domains[d].id = d;
        domains[d].grid = ggrid_alloc(dc->N, dc->L, ghost, c.dt_factor);
        domains[d].cx = dc->cx; domains[d].cy = dc->cy; domains[d].cz = dc->cz;
        domains[d].parent_idx = dc->parent_id;
        domains[d].snap_dt = dc->snap_dt;
        domains[d].next_snap = 0;

        if (dc->parent_id >= 0) {
            double ratio = domains[dc->parent_id].grid->dx / domains[d].grid->dx;
            domains[d].substeps = (int)round(ratio);
            if (domains[d].substeps < 1) domains[d].substeps = 1;
        } else {
            domains[d].substeps = 1;
        }

        GGrid *g = domains[d].grid;
        printf("  Domain %d: N=%d L=%.1f center=(%.1f,%.1f,%.1f) ghost=%d dx=%.4f dt=%.6f",
               d, dc->N, dc->L, dc->cx, dc->cy, dc->cz, ghost, g->dx, g->dt);
        if (dc->parent_id >= 0)
            printf(" parent=%d substeps=%d", dc->parent_id, domains[d].substeps);
        printf(" init=%s\n", dc->init_mode);

        init_domain(&domains[d], dc, &c);

        /* Allocate temporal interpolation buffers for children */
        if (dc->parent_id >= 0) {
            for (int f=0;f<12;f++) { domains[d].pbuf_t0[f]=NULL; domains[d].pbuf_t1[f]=NULL; }
        }
    }

    /* Initial force computation */
    for (int d=0; d<c.n_domains; d++) {
        if (c.dom[d].parent_id < 0)
            compute_forces_periodic(domains[d].grid, &c);
        else
            compute_forces_ghost(domains[d].grid, &c);
    }

    /* SFA output with GDEF */
    uint8_t sfa_dtype = (c.precision==0)?SFA_F16:(c.precision==1)?SFA_F32:SFA_F64;
    SFA *sfa = sfa_create(c.output, c.dom[0].N, c.dom[0].N, c.dom[0].N,
                          c.dom[0].L, c.dom[0].L, c.dom[0].L, domains[0].grid->dt);
    for (int d=0;d<c.n_domains;d++) {
        DomainConfig *dc = &c.dom[d];
        sfa_add_grid(sfa, d, (dc->parent_id<0)?0xFFFF:dc->parent_id,
                     dc->N, dc->N, dc->N, dc->L, dc->L, dc->L,
                     dc->cx, dc->cy, dc->cz, dc->ghost, 0);
    }
    /* KVMD: shared physics */
    {
        char vm[32],veta[32],vmu[32],vkappa[32],vmode[32],vnd[32];
        snprintf(vm,32,"%.6f",sqrt(c.m2)); snprintf(veta,32,"%.6f",c.eta);
        snprintf(vmu,32,"%.6f",c.mu); snprintf(vkappa,32,"%.6f",c.kappa);
        snprintf(vmode,32,"%d",c.mode); snprintf(vnd,32,"%d",c.n_domains);
        const char *keys[]={"m","eta","mu","kappa","mode","n_grids"};
        const char *vals[]={vm,veta,vmu,vkappa,vmode,vnd};
        sfa_add_kvmd(sfa,0,0xFFFFFFFF,0xFFFFFFFF,keys,vals,6);
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
    printf("\nSFA: %s (GDEF: %d grids)\n\n", c.output, c.n_domains);

    /* t=0 snapshots */
    for (int d=0;d<c.n_domains;d++) snap_domain(sfa, &domains[d], 0.0, c.precision);

    /* Diagnostics file */
    FILE *fp = fopen(c.diag_file, "w");
    fprintf(fp,"t\tdomain\tE_total\tE_pot\tphi_max\ttheta_rms\n");

    double dt_coarse = domains[0].grid->dt;
    int n_coarse_steps = (int)(c.T / dt_coarse);
    int diag_every = (int)(c.diag_dt / dt_coarse); if(diag_every<1) diag_every=1;

    /* Initial diagnostics */
    double E0_total = 0;
    for (int d=0;d<c.n_domains;d++) {
        double et,ep,pm,trms;
        domain_energy(domains[d].grid, &c, d==0, &et, &ep, &pm, &trms);
        E0_total += et;
        fprintf(fp,"%.2f\t%d\t%.6e\t%.6e\t%.6e\t%.6e\n", 0.0, d, et, ep, pm, trms);
    }
    printf("INIT: E_total=%.4e (all domains)\n\n", E0_total);

    double wall0 = omp_get_wtime();

    for (int step=1; step<=n_coarse_steps; step++) {
        double t = step * dt_coarse;

        step_multigrid(domains, c.n_domains, &c);

        /* Per-domain snapshots */
        for (int d=0;d<c.n_domains;d++) {
            if (t >= domains[d].next_snap - 0.5*dt_coarse) {
                snap_domain(sfa, &domains[d], t, c.precision);
                domains[d].next_snap = t + c.dom[d].snap_dt;
            }
        }

        /* Diagnostics */
        if (step % diag_every == 0) {
            double et_sum = 0;
            for (int d=0;d<c.n_domains;d++) {
                double et,ep,pm,trms;
                domain_energy(domains[d].grid, &c, d==0, &et, &ep, &pm, &trms);
                et_sum += et;
                fprintf(fp,"%.2f\t%d\t%.6e\t%.6e\t%.6e\t%.6e\n", t, d, et, ep, pm, trms);
            }
            fflush(fp);

            if (step % (diag_every*10) == 0) {
                double wall=omp_get_wtime()-wall0;
                double drift=100*(et_sum-E0_total)/(fabs(E0_total)+1e-30);
                printf("t=%7.1f E=%.3e (drift %+.3f%%) [%.0f%% %.1fs]\n",
                       t, et_sum, drift, 100.0*step/n_coarse_steps, wall);
            }
        }
    }

    /* Final snapshots */
    double t_final = n_coarse_steps * dt_coarse;
    for (int d=0;d<c.n_domains;d++) snap_domain(sfa, &domains[d], t_final, c.precision);

    uint32_t nf = sfa->total_frames;
    sfa_close(sfa);

    /* Summary */
    double wall = omp_get_wtime() - wall0;
    printf("\n=== COMPLETE ===\n");
    printf("SFA: %s (%u frames)\n", c.output, nf);
    for (int d=0;d<c.n_domains;d++) {
        double et,ep,pm,trms;
        domain_energy(domains[d].grid, &c, d==0, &et, &ep, &pm, &trms);
        printf("  Domain %d: E=%.3e Ep=%.1f phi_max=%.3f theta_rms=%.3e\n", d, et, ep, pm, trms);
    }
    printf("Wall: %.1fs (%.1f min)\n", wall, wall/60);

    fclose(fp);
    for (int d=0;d<c.n_domains;d++) {
        ggrid_free(domains[d].grid);
        for (int f=0;f<12;f++) { free(domains[d].pbuf_t0[f]); free(domains[d].pbuf_t1[f]); }
    }
    return 0;
}
