/*  scp_multi.cu — Multi-resolution nested grid Cosserat simulation (CUDA)
 *
 *  GPU version of scp_multi.c — same config, same physics, same SFA output.
 *  Each domain runs its force kernel on GPU. Prolongation/restriction on CPU.
 *  Inter-domain data transfer via D2H/H2D around prolongation/restriction.
 *
 *  Build: nvcc -O3 -arch=sm_70 -o scp_multi_cuda scp_multi.cu -lzstd -lm
 *  Run:   ./scp_multi_cuda config.cfg
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
#define MAX_DOMAINS 16
#define TPB 256

/* ================================================================
   Configuration (identical to scp_multi.c)
   ================================================================ */

typedef struct {
    int N; double L, cx, cy, cz;
    int parent_id, ghost;
    double snap_dt;
    char init_mode[32];
    double A, sigma, R_tube, ellip, A_bg, delta[3];
} DomainConfig;

typedef struct {
    double m2, mtheta2, eta, mu, kappa;
    int mode; double inv_alpha, inv_beta, kappa_gamma;
    double damp_width, damp_rate, T, dt_factor;
    int n_domains; DomainConfig dom[MAX_DOMAINS];
    char output[512], diag_file[512];
    int precision; double diag_dt;
} MultiConfig;

static MultiConfig mcfg_defaults(void) {
    MultiConfig c = {};
    c.m2=2.25;c.mtheta2=0;c.eta=0.5;c.mu=-41.345;c.kappa=50;
    c.mode=0;c.inv_alpha=2.25;c.inv_beta=5;c.kappa_gamma=2;
    c.damp_width=3;c.damp_rate=0.01;c.T=200;c.dt_factor=0.025;
    strcpy(c.output,"multi_output.sfa");strcpy(c.diag_file,"multi_diag.tsv");
    c.precision=1;c.diag_dt=2;
    return c;
}

static void mcfg_set(MultiConfig *c, const char *key, const char *val) {
    if(!strcmp(key,"m")){double m=atof(val);c->m2=m*m;}
    else if(!strcmp(key,"m_theta")){double m=atof(val);c->mtheta2=m*m;}
    else if(!strcmp(key,"eta"))c->eta=atof(val);
    else if(!strcmp(key,"mu"))c->mu=atof(val);
    else if(!strcmp(key,"kappa"))c->kappa=atof(val);
    else if(!strcmp(key,"mode"))c->mode=atoi(val);
    else if(!strcmp(key,"kappa_gamma"))c->kappa_gamma=atof(val);
    else if(!strcmp(key,"damp_width"))c->damp_width=atof(val);
    else if(!strcmp(key,"damp_rate"))c->damp_rate=atof(val);
    else if(!strcmp(key,"T"))c->T=atof(val);
    else if(!strcmp(key,"dt_factor"))c->dt_factor=atof(val);
    else if(!strcmp(key,"output"))strncpy(c->output,val,511);
    else if(!strcmp(key,"diag_file"))strncpy(c->diag_file,val,511);
    else if(!strcmp(key,"diag_dt"))c->diag_dt=atof(val);
    else if(!strcmp(key,"n_grids"))c->n_domains=atoi(val);
    else if(!strcmp(key,"precision")){
        if(!strcmp(val,"f16"))c->precision=0;
        else if(!strcmp(val,"f32"))c->precision=1;
        else if(!strcmp(val,"f64"))c->precision=2;
    }
    else if(!strncmp(key,"grid.",5)){
        int gid=atoi(key+5); if(gid<0||gid>=MAX_DOMAINS)return;
        DomainConfig *d=&c->dom[gid];
        const char *sub=strchr(key+5,'.'); if(!sub)return; sub++;
        if(!strcmp(sub,"N"))d->N=atoi(val);
        else if(!strcmp(sub,"L"))d->L=atof(val);
        else if(!strcmp(sub,"center"))sscanf(val,"%lf,%lf,%lf",&d->cx,&d->cy,&d->cz);
        else if(!strcmp(sub,"parent"))d->parent_id=atoi(val);
        else if(!strcmp(sub,"ghost"))d->ghost=atoi(val);
        else if(!strcmp(sub,"snap_dt"))d->snap_dt=atof(val);
        else if(!strcmp(sub,"init"))strncpy(d->init_mode,val,31);
        else if(!strcmp(sub,"A"))d->A=atof(val);
        else if(!strcmp(sub,"sigma"))d->sigma=atof(val);
        else if(!strcmp(sub,"R_tube"))d->R_tube=atof(val);
        else if(!strcmp(sub,"ellip"))d->ellip=atof(val);
        else if(!strcmp(sub,"A_bg"))d->A_bg=atof(val);
        else if(!strcmp(sub,"delta"))sscanf(val,"%lf,%lf,%lf",&d->delta[0],&d->delta[1],&d->delta[2]);
    }
}

static void mcfg_load(MultiConfig *c, const char *path) {
    FILE *fp=fopen(path,"r"); if(!fp){fprintf(stderr,"Cannot open %s\n",path);exit(1);}
    char line[2048];
    while(fgets(line,sizeof(line),fp)){
        char *p=line; while(*p==' '||*p=='\t')p++;
        if(*p=='#'||*p=='\n'||*p=='\0')continue;
        char *eq=strchr(p,'='); if(!eq)continue; *eq='\0';
        char *key=p,*val=eq+1;
        char *ke=eq-1; while(ke>key&&(*ke==' '||*ke=='\t'))*ke--='\0';
        while(*val==' '||*val=='\t')val++;
        char *ve=val+strlen(val)-1; while(ve>val&&(*ve=='\n'||*ve=='\r'||*ve==' '))*ve--='\0';
        char *hash=strchr(val,'#'); if(hash){*hash='\0';ve=hash-1;while(ve>val&&*ve==' ')*ve--='\0';}
        mcfg_set(c,key,val);
    }
    fclose(fp);
    for(int i=0;i<c->n_domains;i++){
        DomainConfig *d=&c->dom[i];
        if(d->parent_id==0&&i==0)d->parent_id=-1;
        if(d->ghost==0&&d->parent_id>=0)d->ghost=2;
        if(d->snap_dt<=0)d->snap_dt=5;
        if(d->A<=0)d->A=0.8;if(d->A_bg<=0)d->A_bg=0.1;
        if(d->R_tube<=0)d->R_tube=3;if(d->sigma<=0)d->sigma=3;
        if(d->delta[1]==0&&d->delta[2]==0){d->delta[0]=0;d->delta[1]=3.0005;d->delta[2]=4.4325;}
        if(d->init_mode[0]=='\0')strcpy(d->init_mode,"empty");
    }
}

/* ================================================================
   Host grid (same as CPU version)
   ================================================================ */

typedef struct {
    double *mem;
    double *phi[NFIELDS],*phi_vel[NFIELDS],*phi_acc[NFIELDS];
    double *theta[NFIELDS],*theta_vel[NFIELDS],*theta_acc[NFIELDS];
    int N_phys,ghost,N_total; long N3; double L,dx,dt;
} GGrid;

static GGrid *ggrid_alloc(int Np,double L,int gh,double dtf){
    GGrid *g=(GGrid*)calloc(1,sizeof(GGrid));
    g->N_phys=Np;g->ghost=gh;g->N_total=Np+2*gh;
    g->N3=(long)g->N_total*g->N_total*g->N_total;
    g->L=L;g->dx=2.0*L/(Np-1);g->dt=dtf*g->dx;
    long total=18*g->N3;
    g->mem=(double*)malloc(total*sizeof(double));
    if(!g->mem){fprintf(stderr,"FATAL: malloc\n");exit(1);}
    memset(g->mem,0,total*sizeof(double));
    for(int a=0;a<NFIELDS;a++){
        g->phi[a]=g->mem+(0+a)*g->N3; g->phi_vel[a]=g->mem+(3+a)*g->N3;
        g->phi_acc[a]=g->mem+(6+a)*g->N3; g->theta[a]=g->mem+(9+a)*g->N3;
        g->theta_vel[a]=g->mem+(12+a)*g->N3; g->theta_acc[a]=g->mem+(15+a)*g->N3;
    }
    return g;
}
static void ggrid_free(GGrid *g){free(g->mem);free(g);}
static inline long gidx(GGrid *g,int i,int j,int k){return(long)i*g->N_total*g->N_total+j*g->N_total+k;}
static inline double gworld(GGrid *g,double center,int idx){return center-g->L+(idx-g->ghost)*g->dx;}

/* ================================================================
   Domain
   ================================================================ */

typedef struct {
    int id; GGrid *grid;
    double cx,cy,cz;
    int parent_idx,substeps;
    double snap_dt,next_snap;
    /* GPU arrays */
    double *d_phi[3],*d_vel_phi[3],*d_acc_phi[3];
    double *d_theta[3],*d_vel_theta[3],*d_acc_theta[3];
    int gpu_blocks;
    /* Prolongation buffers (host) */
    double *pbuf_t0[12],*pbuf_t1[12];
    int pb_i0,pb_j0,pb_k0,pb_ni,pb_nj,pb_nk;
} Domain;

/* ================================================================
   GPU constant memory (per-kernel-call, set before each domain's kernel)
   ================================================================ */

__constant__ double d_MU,d_KAPPA,d_MASS2,d_MTHETA2,d_ETA;
__constant__ double d_idx2,d_idx1,d_L,d_dx,d_DAMP_WIDTH,d_DAMP_RATE;
__constant__ double d_inv_alpha,d_inv_beta,d_kappa_gamma;
__constant__ int d_N,d_NN,d_MODE,d_ghost,d_Nphys;
__constant__ long d_N3;

static void gpu_set_constants(const MultiConfig *c, GGrid *g) {
    double idx2=1.0/(g->dx*g->dx),idx1=1.0/(2.0*g->dx);
    int N=g->N_total,NN=N*N,gh=g->ghost,Np=g->N_phys;
    long N3=g->N3;
    cudaMemcpyToSymbol(d_MU,&c->mu,8);cudaMemcpyToSymbol(d_KAPPA,&c->kappa,8);
    cudaMemcpyToSymbol(d_MASS2,&c->m2,8);cudaMemcpyToSymbol(d_MTHETA2,&c->mtheta2,8);
    cudaMemcpyToSymbol(d_ETA,&c->eta,8);cudaMemcpyToSymbol(d_idx2,&idx2,8);
    cudaMemcpyToSymbol(d_idx1,&idx1,8);cudaMemcpyToSymbol(d_L,&g->L,8);
    cudaMemcpyToSymbol(d_dx,&g->dx,8);cudaMemcpyToSymbol(d_DAMP_WIDTH,&c->damp_width,8);
    cudaMemcpyToSymbol(d_DAMP_RATE,&c->damp_rate,8);
    cudaMemcpyToSymbol(d_inv_alpha,&c->inv_alpha,8);cudaMemcpyToSymbol(d_inv_beta,&c->inv_beta,8);
    cudaMemcpyToSymbol(d_kappa_gamma,&c->kappa_gamma,8);
    cudaMemcpyToSymbol(d_N,&N,4);cudaMemcpyToSymbol(d_NN,&NN,4);
    cudaMemcpyToSymbol(d_MODE,&c->mode,4);cudaMemcpyToSymbol(d_ghost,&gh,4);
    cudaMemcpyToSymbol(d_Nphys,&Np,4);cudaMemcpyToSymbol(d_N3,&N3,8);
}

/* ================================================================
   GPU kernels
   ================================================================ */

/* Force kernel for ghost-zone domains (no wrapping) */
__global__ void forces_ghost_kernel(
    const double *p0,const double *p1,const double *p2,
    const double *t0,const double *t1,const double *t2,
    double *ap0,double *ap1,double *ap2,
    double *at0,double *at1,double *at2)
{
    long idx=(long)blockIdx.x*blockDim.x+threadIdx.x;
    if(idx>=d_N3)return;
    int N=d_N,NN=d_NN,g0=d_ghost,g1=d_ghost+d_Nphys;
    int i=(int)(idx/NN),j=(int)((idx/N)%N),k=(int)(idx%N);
    if(i<g0||i>=g1||j<g0||j>=g1||k<g0||k>=g1)return;

    long nip=(long)(i+1)*NN+j*N+k,nim=(long)(i-1)*NN+j*N+k;
    long njp=(long)i*NN+(j+1)*N+k,njm=(long)i*NN+(j-1)*N+k;
    long nkp=(long)i*NN+j*N+(k+1),nkm=(long)i*NN+j*N+(k-1);

    double pp0=p0[idx],pp1=p1[idx],pp2=p2[idx];
    double sig=pp0*pp0+pp1*pp1+pp2*pp2;
    double me2=(d_MODE==1)?d_inv_alpha/(1.0+d_inv_beta*sig):d_MASS2;
    double keff=(d_MODE==3)?d_KAPPA/(1.0+d_kappa_gamma*sig):d_KAPPA;
    double P=pp0*pp1*pp2,P2=P*P;
    double den=1.0+keff*P2;
    double dVdP=d_MU*P/(den*den);
    double t2c=0;
    if(d_MODE==3&&d_kappa_gamma>0){double D=1.0+d_kappa_gamma*sig+d_KAPPA*P2;t2c=d_MU*d_kappa_gamma*d_KAPPA*P2*P2/(D*D);}

    const double *phi[3]={p0,p1,p2},*th[3]={t0,t1,t2};
    for(int a=0;a<3;a++){
        double lap=(phi[a][nip]+phi[a][nim]+phi[a][njp]+phi[a][njm]+phi[a][nkp]+phi[a][nkm]-6.0*phi[a][idx])*d_idx2;
        double dPda=(a==0)?pp1*pp2:(a==1)?pp0*pp2:pp0*pp1;
        double ct;
        if(a==0)ct=(th[2][njp]-th[2][njm]-th[1][nkp]+th[1][nkm])*d_idx1;
        else if(a==1)ct=(th[0][nkp]-th[0][nkm]-th[2][nip]+th[2][nim])*d_idx1;
        else ct=(th[1][nip]-th[1][nim]-th[0][njp]+th[0][njm])*d_idx1;
        double *acc=(a==0)?ap0:(a==1)?ap1:ap2;
        acc[idx]=lap-me2*phi[a][idx]-dVdP*dPda-t2c*phi[a][idx]+d_ETA*ct;
    }
    for(int a=0;a<3;a++){
        double lapt=(th[a][nip]+th[a][nim]+th[a][njp]+th[a][njm]+th[a][nkp]+th[a][nkm]-6.0*th[a][idx])*d_idx2;
        double cp;
        if(a==0)cp=(phi[2][njp]-phi[2][njm]-phi[1][nkp]+phi[1][nkm])*d_idx1;
        else if(a==1)cp=(phi[0][nkp]-phi[0][nkm]-phi[2][nip]+phi[2][nim])*d_idx1;
        else cp=(phi[1][nip]-phi[1][nim]-phi[0][njp]+phi[0][njm])*d_idx1;
        double *acc=(a==0)?at0:(a==1)?at1:at2;
        acc[idx]=lapt-d_MTHETA2*th[a][idx]+d_ETA*cp;
    }
}

/* Force kernel for root domain (periodic wrapping) */
__global__ void forces_periodic_kernel(
    const double *p0,const double *p1,const double *p2,
    const double *t0,const double *t1,const double *t2,
    double *ap0,double *ap1,double *ap2,
    double *at0,double *at1,double *at2)
{
    long idx=(long)blockIdx.x*blockDim.x+threadIdx.x;
    if(idx>=d_N3)return;
    int N=d_N,NN=d_NN;
    int i=(int)(idx/NN),j=(int)((idx/N)%N),k=(int)(idx%N);
    int ip=(i+1)%N,im=(i-1+N)%N,jp=(j+1)%N,jm=(j-1+N)%N,kp=(k+1)%N,km=(k-1+N)%N;
    long nip=(long)ip*NN+j*N+k,nim=(long)im*NN+j*N+k;
    long njp=(long)i*NN+jp*N+k,njm=(long)i*NN+jm*N+k;
    long nkp=(long)i*NN+j*N+kp,nkm=(long)i*NN+j*N+km;

    double pp0=p0[idx],pp1=p1[idx],pp2=p2[idx];
    double sig=pp0*pp0+pp1*pp1+pp2*pp2;
    double me2=(d_MODE==1)?d_inv_alpha/(1.0+d_inv_beta*sig):d_MASS2;
    double keff=(d_MODE==3)?d_KAPPA/(1.0+d_kappa_gamma*sig):d_KAPPA;
    double P=pp0*pp1*pp2,P2=P*P;
    double den=1.0+keff*P2;
    double dVdP=d_MU*P/(den*den);
    double t2c=0;
    if(d_MODE==3&&d_kappa_gamma>0){double D=1.0+d_kappa_gamma*sig+d_KAPPA*P2;t2c=d_MU*d_kappa_gamma*d_KAPPA*P2*P2/(D*D);}

    const double *phi[3]={p0,p1,p2},*th[3]={t0,t1,t2};
    for(int a=0;a<3;a++){
        double lap=(phi[a][nip]+phi[a][nim]+phi[a][njp]+phi[a][njm]+phi[a][nkp]+phi[a][nkm]-6.0*phi[a][idx])*d_idx2;
        double dPda=(a==0)?pp1*pp2:(a==1)?pp0*pp2:pp0*pp1;
        double ct;
        if(a==0)ct=(th[2][njp]-th[2][njm]-th[1][nkp]+th[1][nkm])*d_idx1;
        else if(a==1)ct=(th[0][nkp]-th[0][nkm]-th[2][nip]+th[2][nim])*d_idx1;
        else ct=(th[1][nip]-th[1][nim]-th[0][njp]+th[0][njm])*d_idx1;
        double *acc=(a==0)?ap0:(a==1)?ap1:ap2;
        acc[idx]=lap-me2*phi[a][idx]-dVdP*dPda-t2c*phi[a][idx]+d_ETA*ct;
    }
    for(int a=0;a<3;a++){
        double lapt=(th[a][nip]+th[a][nim]+th[a][njp]+th[a][njm]+th[a][nkp]+th[a][nkm]-6.0*th[a][idx])*d_idx2;
        double cp;
        if(a==0)cp=(phi[2][njp]-phi[2][njm]-phi[1][nkp]+phi[1][nkm])*d_idx1;
        else if(a==1)cp=(phi[0][nkp]-phi[0][nkm]-phi[2][nip]+phi[2][nim])*d_idx1;
        else cp=(phi[1][nip]-phi[1][nim]-phi[0][njp]+phi[0][njm])*d_idx1;
        double *acc=(a==0)?at0:(a==1)?at1:at2;
        acc[idx]=lapt-d_MTHETA2*th[a][idx]+d_ETA*cp;
    }
}

__global__ void halfkick_kernel(double *vp0,double *vp1,double *vp2,
    double *vt0,double *vt1,double *vt2,
    const double *ap0,const double *ap1,const double *ap2,
    const double *at0,const double *at1,const double *at2,double hdt){
    long idx=(long)blockIdx.x*blockDim.x+threadIdx.x;
    if(idx>=d_N3)return;
    vp0[idx]+=hdt*ap0[idx];vp1[idx]+=hdt*ap1[idx];vp2[idx]+=hdt*ap2[idx];
    vt0[idx]+=hdt*at0[idx];vt1[idx]+=hdt*at1[idx];vt2[idx]+=hdt*at2[idx];
}

__global__ void drift_kernel(double *p0,double *p1,double *p2,
    double *t0,double *t1,double *t2,
    const double *vp0,const double *vp1,const double *vp2,
    const double *vt0,const double *vt1,const double *vt2,double dt){
    long idx=(long)blockIdx.x*blockDim.x+threadIdx.x;
    if(idx>=d_N3)return;
    p0[idx]+=dt*vp0[idx];p1[idx]+=dt*vp1[idx];p2[idx]+=dt*vp2[idx];
    t0[idx]+=dt*vt0[idx];t1[idx]+=dt*vt1[idx];t2[idx]+=dt*vt2[idx];
}

__global__ void absorb_kernel(double *vp0,double *vp1,double *vp2,
    double *vt0,double *vt1,double *vt2){
    long idx=(long)blockIdx.x*blockDim.x+threadIdx.x;
    if(idx>=d_N3)return;
    if(d_DAMP_WIDTH<=0||d_DAMP_RATE<=0)return;
    int N=d_N,NN=d_NN;
    int i=(int)(idx/NN),j=(int)((idx/N)%N),k=(int)(idx%N);
    double x=-d_L+i*d_dx,y=-d_L+j*d_dx,z=-d_L+k*d_dx;
    double r=sqrt(x*x+y*y+z*z);
    double Rd=d_L-d_DAMP_WIDTH;
    if(r>Rd){double s=(r-Rd)/d_DAMP_WIDTH;if(s>1)s=1;double d=1.0-d_DAMP_RATE*s*s;
        vp0[idx]*=d;vp1[idx]*=d;vp2[idx]*=d;vt0[idx]*=d;vt1[idx]*=d;vt2[idx]*=d;}
}

/* ================================================================
   GPU memory management per domain
   ================================================================ */

static void dom_gpu_alloc(Domain *dom) {
    size_t bytes=dom->grid->N3*sizeof(double);
    for(int a=0;a<3;a++){
        cudaMalloc(&dom->d_phi[a],bytes);cudaMalloc(&dom->d_vel_phi[a],bytes);cudaMalloc(&dom->d_acc_phi[a],bytes);
        cudaMalloc(&dom->d_theta[a],bytes);cudaMalloc(&dom->d_vel_theta[a],bytes);cudaMalloc(&dom->d_acc_theta[a],bytes);
    }
    dom->gpu_blocks=(int)((dom->grid->N3+TPB-1)/TPB);
}

static void dom_gpu_upload(Domain *dom) {
    size_t b=dom->grid->N3*sizeof(double);
    for(int a=0;a<3;a++){
        cudaMemcpy(dom->d_phi[a],dom->grid->phi[a],b,cudaMemcpyHostToDevice);
        cudaMemcpy(dom->d_vel_phi[a],dom->grid->phi_vel[a],b,cudaMemcpyHostToDevice);
        cudaMemcpy(dom->d_theta[a],dom->grid->theta[a],b,cudaMemcpyHostToDevice);
        cudaMemcpy(dom->d_vel_theta[a],dom->grid->theta_vel[a],b,cudaMemcpyHostToDevice);
        /* Do NOT zero accelerations — Berger-Oliger needs them for next half-kick */
        cudaMemcpy(dom->d_acc_phi[a],dom->grid->phi_acc[a],b,cudaMemcpyHostToDevice);
        cudaMemcpy(dom->d_acc_theta[a],dom->grid->theta_acc[a],b,cudaMemcpyHostToDevice);
    }
}

static void dom_gpu_download(Domain *dom) {
    size_t b=dom->grid->N3*sizeof(double);
    for(int a=0;a<3;a++){
        cudaMemcpy(dom->grid->phi[a],dom->d_phi[a],b,cudaMemcpyDeviceToHost);
        cudaMemcpy(dom->grid->phi_vel[a],dom->d_vel_phi[a],b,cudaMemcpyDeviceToHost);
        cudaMemcpy(dom->grid->phi_acc[a],dom->d_acc_phi[a],b,cudaMemcpyDeviceToHost);
        cudaMemcpy(dom->grid->theta[a],dom->d_theta[a],b,cudaMemcpyDeviceToHost);
        cudaMemcpy(dom->grid->theta_vel[a],dom->d_vel_theta[a],b,cudaMemcpyDeviceToHost);
        cudaMemcpy(dom->grid->theta_acc[a],dom->d_acc_theta[a],b,cudaMemcpyDeviceToHost);
    }
}

static void dom_gpu_free(Domain *dom) {
    for(int a=0;a<3;a++){
        cudaFree(dom->d_phi[a]);cudaFree(dom->d_vel_phi[a]);cudaFree(dom->d_acc_phi[a]);
        cudaFree(dom->d_theta[a]);cudaFree(dom->d_vel_theta[a]);cudaFree(dom->d_acc_theta[a]);
    }
}

/* GPU step helpers */
static void gpu_kick(Domain *dom, double hdt) {
    halfkick_kernel<<<dom->gpu_blocks,TPB>>>(
        dom->d_vel_phi[0],dom->d_vel_phi[1],dom->d_vel_phi[2],
        dom->d_vel_theta[0],dom->d_vel_theta[1],dom->d_vel_theta[2],
        dom->d_acc_phi[0],dom->d_acc_phi[1],dom->d_acc_phi[2],
        dom->d_acc_theta[0],dom->d_acc_theta[1],dom->d_acc_theta[2],hdt);
}
static void gpu_drift(Domain *dom, double dt) {
    drift_kernel<<<dom->gpu_blocks,TPB>>>(
        dom->d_phi[0],dom->d_phi[1],dom->d_phi[2],
        dom->d_theta[0],dom->d_theta[1],dom->d_theta[2],
        dom->d_vel_phi[0],dom->d_vel_phi[1],dom->d_vel_phi[2],
        dom->d_vel_theta[0],dom->d_vel_theta[1],dom->d_vel_theta[2],dt);
}
static void gpu_forces(Domain *dom, const MultiConfig *c, int is_root) {
    gpu_set_constants(c, dom->grid);
    if (is_root)
        forces_periodic_kernel<<<dom->gpu_blocks,TPB>>>(
            dom->d_phi[0],dom->d_phi[1],dom->d_phi[2],
            dom->d_theta[0],dom->d_theta[1],dom->d_theta[2],
            dom->d_acc_phi[0],dom->d_acc_phi[1],dom->d_acc_phi[2],
            dom->d_acc_theta[0],dom->d_acc_theta[1],dom->d_acc_theta[2]);
    else
        forces_ghost_kernel<<<dom->gpu_blocks,TPB>>>(
            dom->d_phi[0],dom->d_phi[1],dom->d_phi[2],
            dom->d_theta[0],dom->d_theta[1],dom->d_theta[2],
            dom->d_acc_phi[0],dom->d_acc_phi[1],dom->d_acc_phi[2],
            dom->d_acc_theta[0],dom->d_acc_theta[1],dom->d_acc_theta[2]);
}
static void gpu_absorb(Domain *dom, const MultiConfig *c) {
    gpu_set_constants(c, dom->grid);
    absorb_kernel<<<dom->gpu_blocks,TPB>>>(
        dom->d_vel_phi[0],dom->d_vel_phi[1],dom->d_vel_phi[2],
        dom->d_vel_theta[0],dom->d_vel_theta[1],dom->d_vel_theta[2]);
}

/* ================================================================
   Prolongation / restriction (CPU, same as scp_multi.c)
   Requires D2H before, H2D after.
   ================================================================ */

static inline double trilinear(const double *buf,int ni,int nj,int nk,double fi,double fj,double fk){
    int i0=(int)floor(fi),j0=(int)floor(fj),k0=(int)floor(fk);
    double di=fi-i0,dj=fj-j0,dk=fk-k0;
    int i1=i0+1,j1=j0+1,k1=k0+1;
    if(i0<0)i0=0;if(j0<0)j0=0;if(k0<0)k0=0;
    if(i1>=ni)i1=ni-1;if(j1>=nj)j1=nj-1;if(k1>=nk)k1=nk-1;
    if(i0>=ni)i0=ni-1;if(j0>=nj)j0=nj-1;if(k0>=nk)k0=nk-1;
    double c000=buf[(long)i0*nj*nk+j0*nk+k0],c100=buf[(long)i1*nj*nk+j0*nk+k0];
    double c010=buf[(long)i0*nj*nk+j1*nk+k0],c110=buf[(long)i1*nj*nk+j1*nk+k0];
    double c001=buf[(long)i0*nj*nk+j0*nk+k1],c101=buf[(long)i1*nj*nk+j0*nk+k1];
    double c011=buf[(long)i0*nj*nk+j1*nk+k1],c111=buf[(long)i1*nj*nk+j1*nk+k1];
    double c00=c000*(1-di)+c100*di,c01=c001*(1-di)+c101*di;
    double c10=c010*(1-di)+c110*di,c11=c011*(1-di)+c111*di;
    double c0=c00*(1-dj)+c10*dj,c1=c01*(1-dj)+c11*dj;
    return c0*(1-dk)+c1*dk;
}

static void fill_parent_buffer(Domain *child,Domain *parent,double **buf){
    GGrid *pg=parent->grid,*cg=child->grid;
    int Np=pg->N_phys;
    double cxlo=child->cx-cg->L-cg->ghost*cg->dx;
    double cylo=child->cy-cg->L-cg->ghost*cg->dx;
    double czlo=child->cz-cg->L-cg->ghost*cg->dx;
    int pi0=(int)floor((cxlo-(parent->cx-pg->L))/pg->dx)-1;
    int pj0=(int)floor((cylo-(parent->cy-pg->L))/pg->dx)-1;
    int pk0=(int)floor((czlo-(parent->cz-pg->L))/pg->dx)-1;
    if(pi0<0)pi0=0;if(pj0<0)pj0=0;if(pk0<0)pk0=0;
    int pi1=pi0+(int)ceil(2.0*cg->L/pg->dx)+4;
    int pj1=pj0+(int)ceil(2.0*cg->L/pg->dx)+4;
    int pk1=pk0+(int)ceil(2.0*cg->L/pg->dx)+4;
    if(pi1>Np)pi1=Np;if(pj1>Np)pj1=Np;if(pk1>Np)pk1=Np;
    child->pb_i0=pi0;child->pb_j0=pj0;child->pb_k0=pk0;
    child->pb_ni=pi1-pi0;child->pb_nj=pj1-pj0;child->pb_nk=pk1-pk0;
    long pb_n3=(long)child->pb_ni*child->pb_nj*child->pb_nk;
    double *arrays[12]={pg->phi[0],pg->phi[1],pg->phi[2],pg->theta[0],pg->theta[1],pg->theta[2],
        pg->phi_vel[0],pg->phi_vel[1],pg->phi_vel[2],pg->theta_vel[0],pg->theta_vel[1],pg->theta_vel[2]};
    for(int f=0;f<12;f++){
        if(!buf[f])buf[f]=(double*)malloc(pb_n3*sizeof(double));
        for(int i=0;i<child->pb_ni;i++)for(int j=0;j<child->pb_nj;j++)for(int k=0;k<child->pb_nk;k++){
            int pi=pi0+i,pj=pj0+j,pk=pk0+k;
            pi=((pi%Np)+Np)%Np;pj=((pj%Np)+Np)%Np;pk=((pk%Np)+Np)%Np;
            buf[f][(long)i*child->pb_nj*child->pb_nk+j*child->pb_nk+k]=arrays[f][(long)pi*Np*Np+pj*Np+pk];
        }
    }
}

static void prolongate(Domain *child,Domain *parent,double alpha){
    GGrid *cg=child->grid,*pg=parent->grid;
    int Nt=cg->N_total,NN=Nt*Nt,g0=cg->ghost,g1=cg->ghost+cg->N_phys;
    double *arrays[12]={cg->phi[0],cg->phi[1],cg->phi[2],cg->theta[0],cg->theta[1],cg->theta[2],
        cg->phi_vel[0],cg->phi_vel[1],cg->phi_vel[2],cg->theta_vel[0],cg->theta_vel[1],cg->theta_vel[2]};
    for(int i=0;i<Nt;i++)for(int j=0;j<Nt;j++)for(int k=0;k<Nt;k++){
        if(i>=g0&&i<g1&&j>=g0&&j<g1&&k>=g0&&k<g1)continue;
        long idx=(long)i*NN+j*Nt+k;
        double wx=gworld(cg,child->cx,i),wy=gworld(cg,child->cy,j),wz=gworld(cg,child->cz,k);
        double fi=(wx-(parent->cx-pg->L))/pg->dx-child->pb_i0;
        double fj=(wy-(parent->cy-pg->L))/pg->dx-child->pb_j0;
        double fk=(wz-(parent->cz-pg->L))/pg->dx-child->pb_k0;
        for(int f=0;f<12;f++){
            double v0=trilinear(child->pbuf_t0[f],child->pb_ni,child->pb_nj,child->pb_nk,fi,fj,fk);
            double v1=trilinear(child->pbuf_t1[f],child->pb_ni,child->pb_nj,child->pb_nk,fi,fj,fk);
            arrays[f][idx]=(1.0-alpha)*v0+alpha*v1;
        }
    }
}

static void restrict_to_parent(Domain *child,Domain *parent){
    GGrid *cg=child->grid,*pg=parent->grid;
    int Np=pg->N_phys,Ct=cg->N_total,cg0=cg->ghost,Nc=cg->N_phys;
    double ratio=pg->dx/cg->dx;int iratio=(int)round(ratio);if(iratio<1)iratio=1;
    double *pa[12]={pg->phi[0],pg->phi[1],pg->phi[2],pg->theta[0],pg->theta[1],pg->theta[2],
        pg->phi_vel[0],pg->phi_vel[1],pg->phi_vel[2],pg->theta_vel[0],pg->theta_vel[1],pg->theta_vel[2]};
    double *ca[12]={cg->phi[0],cg->phi[1],cg->phi[2],cg->theta[0],cg->theta[1],cg->theta[2],
        cg->phi_vel[0],cg->phi_vel[1],cg->phi_vel[2],cg->theta_vel[0],cg->theta_vel[1],cg->theta_vel[2]};
    for(int pi=0;pi<Np;pi++)for(int pj=0;pj<Np;pj++)for(int pk=0;pk<Np;pk++){
        double px=parent->cx-pg->L+pi*pg->dx,py=parent->cy-pg->L+pj*pg->dx,pz=parent->cz-pg->L+pk*pg->dx;
        if(px<child->cx-cg->L||px>child->cx+cg->L)continue;
        if(py<child->cy-cg->L||py>child->cy+cg->L)continue;
        if(pz<child->cz-cg->L||pz>child->cz+cg->L)continue;
        int ci0=cg0+(int)round((px-(child->cx-cg->L))/cg->dx);
        int cj0=cg0+(int)round((py-(child->cy-cg->L))/cg->dx);
        int ck0=cg0+(int)round((pz-(child->cz-cg->L))/cg->dx);
        long pidx=(long)pi*Np*Np+pj*Np+pk;
        for(int f=0;f<12;f++){
            double sum=0;int cnt=0;
            for(int di=0;di<iratio;di++)for(int dj=0;dj<iratio;dj++)for(int dk=0;dk<iratio;dk++){
                int ci=ci0+di,cj=cj0+dj,ck=ck0+dk;
                if(ci>=cg0&&ci<cg0+Nc&&cj>=cg0&&cj<cg0+Nc&&ck>=cg0&&ck<cg0+Nc){
                    sum+=ca[f][(long)ci*Ct*Ct+cj*Ct+ck];cnt++;
                }
            }
            if(cnt>0)pa[f][pidx]=sum/cnt;
        }
    }
}

/* ================================================================
   Initialization (same as CPU scp_multi.c)
   ================================================================ */

static void init_domain(Domain *dom,const DomainConfig *dc,const MultiConfig *c){
    GGrid *g=dom->grid;int Nt=g->N_total,NN=Nt*Nt,g0=g->ghost;
    double kw=PI/g->L,omega=sqrt(kw*kw+c->m2);
    double k_bg=PI/g->L,omega_bg=sqrt(k_bg*k_bg+c->m2);
    if(!strcmp(dc->init_mode,"braid")){
        double sx=1+dc->ellip,sy=1-dc->ellip,inv2R2=1.0/(2*dc->R_tube*dc->R_tube);
        for(int i=g0;i<g0+g->N_phys;i++)for(int j=g0;j<g0+g->N_phys;j++)
        for(int k=g0;k<g0+g->N_phys;k++){
            long idx=(long)i*NN+j*Nt+k;
            double x=gworld(g,dom->cx,i),y=gworld(g,dom->cy,j),z=gworld(g,dom->cz,k);
            double r2e=x*x/(sx*sx)+y*y/(sy*sy);double env=exp(-r2e*inv2R2);
            for(int a=0;a<NFIELDS;a++){
                double ph=kw*z+dc->delta[a],ph_bg=k_bg*z+2*PI*a/3.0;
                g->phi[a][idx]=dc->A*env*cos(ph)+dc->A_bg*cos(ph_bg);
                g->phi_vel[a][idx]=omega*dc->A*env*sin(ph)+omega_bg*dc->A_bg*sin(ph_bg);
            }
        }
    } else if(!strcmp(dc->init_mode,"oscillon")){
        for(int i=g0;i<g0+g->N_phys;i++)for(int j=g0;j<g0+g->N_phys;j++)
        for(int k=g0;k<g0+g->N_phys;k++){
            long idx=(long)i*NN+j*Nt+k;
            double x=gworld(g,dom->cx,i),y=gworld(g,dom->cy,j),z=gworld(g,dom->cz,k);
            double r2=x*x+y*y+z*z;double env=dc->A*exp(-r2/(2.0*dc->sigma*dc->sigma));
            for(int a=0;a<NFIELDS;a++)g->phi[a][idx]=env*cos(dc->delta[a]);
        }
    } else { /* empty: background only */
        for(int i=g0;i<g0+g->N_phys;i++)for(int j=g0;j<g0+g->N_phys;j++)
        for(int k=g0;k<g0+g->N_phys;k++){
            long idx=(long)i*NN+j*Nt+k;double z=gworld(g,dom->cz,k);
            for(int a=0;a<NFIELDS;a++){
                double ph_bg=k_bg*z+2*PI*a/3.0;
                g->phi[a][idx]=dc->A_bg*cos(ph_bg);
                g->phi_vel[a][idx]=omega_bg*dc->A_bg*sin(ph_bg);
            }
        }
    }
}

/* ================================================================
   Diagnostics + SFA output (host-side)
   ================================================================ */

static void domain_energy(GGrid *g,const MultiConfig *c,int is_root,
    double *E_tot,double *E_pot,double *phi_max,double *theta_rms_out){
    int Nt=g->N_total,NN=Nt*Nt;
    int g0=is_root?0:g->ghost,g1=is_root?g->N_phys:(g->ghost+g->N_phys);
    double dx=g->dx,dV=dx*dx*dx,idx1=1.0/(2.0*dx);
    double epk=0,etk=0,eg=0,em=0,ep=0,pm=0,trms=0;
    for(int i=g0;i<g1;i++)for(int j=g0;j<g1;j++)for(int k=g0;k<g1;k++){
        long idx=(long)i*NN+j*Nt+k;
        long nip,nim,njp,njm,nkp,nkm;
        if(is_root){int N=g->N_phys;
            nip=(long)((i+1)%N)*N*N+j*N+k;nim=(long)((i-1+N)%N)*N*N+j*N+k;
            njp=(long)i*N*N+((j+1)%N)*N+k;njm=(long)i*N*N+((j-1+N)%N)*N+k;
            nkp=(long)i*N*N+j*N+((k+1)%N);nkm=(long)i*N*N+j*N+((k-1+N)%N);
        } else {
            nip=(long)(i+1)*NN+j*Nt+k;nim=(long)(i-1)*NN+j*Nt+k;
            njp=(long)i*NN+(j+1)*Nt+k;njm=(long)i*NN+(j-1)*Nt+k;
            nkp=(long)i*NN+j*Nt+(k+1);nkm=(long)i*NN+j*Nt+(k-1);
        }
        double p0=g->phi[0][idx],p1=g->phi[1][idx],p2=g->phi[2][idx];
        for(int a=0;a<NFIELDS;a++){
            epk+=0.5*g->phi_vel[a][idx]*g->phi_vel[a][idx]*dV;
            etk+=0.5*g->theta_vel[a][idx]*g->theta_vel[a][idx]*dV;
            double gx=(g->phi[a][nip]-g->phi[a][nim])*idx1;
            double gy=(g->phi[a][njp]-g->phi[a][njm])*idx1;
            double gz=(g->phi[a][nkp]-g->phi[a][nkm])*idx1;
            eg+=0.5*(gx*gx+gy*gy+gz*gz)*dV;
            em+=0.5*c->m2*g->phi[a][idx]*g->phi[a][idx]*dV;
            double ap=fabs(g->phi[a][idx]);if(ap>pm)pm=ap;
            trms+=g->theta[a][idx]*g->theta[a][idx];
        }
        double P=p0*p1*p2,P2=P*P;
        ep+=(c->mu/2.0)*P2/(1.0+c->kappa*P2)*dV;
    }
    long nphys3=(long)(g1-g0)*(g1-g0)*(g1-g0);
    *E_tot=epk+etk+eg+em+ep;*E_pot=ep;*phi_max=pm;
    *theta_rms_out=sqrt(trms/(3.0*nphys3));
}

static inline uint16_t f64_to_f16(double v){
    float f=(float)v;uint32_t x;memcpy(&x,&f,4);
    uint16_t sign=(x>>16)&0x8000;int exp=((x>>23)&0xFF)-127+15;uint16_t mant=(x>>13)&0x3FF;
    if(exp<=0)return sign;if(exp>=31)return sign|0x7C00;return sign|(exp<<10)|mant;
}

static void snap_domain(SFA *sfa,Domain *dom,double t,int precision){
    GGrid *g=dom->grid;int Np=g->N_phys,Nt=g->N_total,g0=g->ghost;
    long Np3=(long)Np*Np*Np;
    double *phys[12];
    double *arrays[12]={g->phi[0],g->phi[1],g->phi[2],g->theta[0],g->theta[1],g->theta[2],
        g->phi_vel[0],g->phi_vel[1],g->phi_vel[2],g->theta_vel[0],g->theta_vel[1],g->theta_vel[2]};
    for(int f=0;f<12;f++){
        phys[f]=(double*)malloc(Np3*sizeof(double));long out=0;
        for(int i=g0;i<g0+Np;i++)for(int j=g0;j<g0+Np;j++)for(int k=g0;k<g0+Np;k++)
            phys[f][out++]=arrays[f][(long)i*Nt*Nt+j*Nt+k];
    }
    uint64_t bpp=0;
    for(uint32_t c=0;c<sfa->n_columns;c++)bpp+=sfa_dtype_size[sfa->columns[c].dtype];
    uint64_t fb=(uint64_t)Np3*bpp;
    if(precision==2){sfa_write_frame_ex(sfa,t,(void**)phys,(uint32_t)dom->id,fb);}
    else{
        void *cols[12];int es=(precision==0)?2:4;
        for(int f=0;f<12;f++){
            cols[f]=malloc(Np3*es);
            if(precision==1){float *p=(float*)cols[f];for(long i=0;i<Np3;i++)p[i]=(float)phys[f][i];}
            else{uint16_t *p=(uint16_t*)cols[f];for(long i=0;i<Np3;i++)p[i]=f64_to_f16(phys[f][i]);}
        }
        sfa_write_frame_ex(sfa,t,cols,(uint32_t)dom->id,fb);
        for(int f=0;f<12;f++)free(cols[f]);
    }
    for(int f=0;f<12;f++)free(phys[f]);
}

/* ================================================================
   Berger-Oliger step (GPU)
   ================================================================ */

static void step_multigrid_gpu(Domain *doms,int n_dom,const MultiConfig *c){
    Domain *root=&doms[0];
    double dt_c=root->grid->dt;

    /* 1. Download root, fill parent buffers at t_n */
    dom_gpu_download(root);
    for(int d=1;d<n_dom;d++)
        fill_parent_buffer(&doms[d],&doms[doms[d].parent_idx],doms[d].pbuf_t0);

    /* 2. Root: full Verlet step on GPU */
    gpu_set_constants(c,root->grid);
    gpu_kick(root,0.5*dt_c);
    gpu_drift(root,dt_c);
    gpu_forces(root,c,1);
    gpu_kick(root,0.5*dt_c);
    gpu_absorb(root, c);
    cudaDeviceSynchronize();

    /* 3. Download root at t_{n+1}, fill parent buffers */
    dom_gpu_download(root);
    for(int d=1;d<n_dom;d++)
        fill_parent_buffer(&doms[d],&doms[doms[d].parent_idx],doms[d].pbuf_t1);

    /* 4. Subcycle each child */
    for(int d=1;d<n_dom;d++){
        Domain *child=&doms[d];
        double dt_f=child->grid->dt;
        int substeps=child->substeps;

        for(int sub=0;sub<substeps;sub++){
            double alpha=((double)sub+0.5)/substeps;
            /* Prolongate on CPU, upload ghost zones */
            dom_gpu_download(child); /* get current state */
            prolongate(child,&doms[child->parent_idx],alpha);
            dom_gpu_upload(child); /* push with filled ghost zones */
            /* Verlet step on GPU */
            gpu_set_constants(c,child->grid);
            gpu_kick(child,0.5*dt_f);
            gpu_drift(child,dt_f);
            gpu_forces(child,c,0);
            gpu_kick(child,0.5*dt_f);
            cudaDeviceSynchronize();
        }

        /* 5. Restrict child back to parent */
        dom_gpu_download(child);
        restrict_to_parent(child,&doms[child->parent_idx]);
        dom_gpu_upload(root); /* push restricted data back to GPU */
    }
}

/* ================================================================
   Main
   ================================================================ */

int main(int argc,char **argv){
    MultiConfig c=mcfg_defaults();
    if(argc<2){fprintf(stderr,"Usage: %s config.cfg\n",argv[0]);return 1;}
    mcfg_load(&c,argv[1]);
    for(int i=2;i<argc-1;i+=2){const char *k=argv[i];if(k[0]=='-')k++;mcfg_set(&c,k,argv[i+1]);}

    cudaDeviceProp prop;cudaGetDeviceProperties(&prop,0);
    printf("GPU: %s (%.1f GB)\n\n",prop.name,prop.totalGlobalMem/1e9);
    printf("=== scp_multi CUDA: Multi-Resolution Cosserat (GPU) ===\n");
    printf("Physics: m²=%.4f η=%.3f μ=%.3f κ=%.1f\n",c.m2,c.eta,c.mu,c.kappa);
    printf("Domains: %d\n\n",c.n_domains);

    Domain doms[MAX_DOMAINS]={};
    for(int d=0;d<c.n_domains;d++){
        DomainConfig *dc=&c.dom[d];
        int ghost=(dc->parent_id<0)?0:dc->ghost;
        doms[d].id=d;doms[d].grid=ggrid_alloc(dc->N,dc->L,ghost,c.dt_factor);
        doms[d].cx=dc->cx;doms[d].cy=dc->cy;doms[d].cz=dc->cz;
        doms[d].parent_idx=dc->parent_id;doms[d].snap_dt=dc->snap_dt;doms[d].next_snap=0;
        if(dc->parent_id>=0){
            double ratio=doms[dc->parent_id].grid->dx/doms[d].grid->dx;
            doms[d].substeps=(int)round(ratio);if(doms[d].substeps<1)doms[d].substeps=1;
        } else doms[d].substeps=1;
        GGrid *g=doms[d].grid;
        printf("  Domain %d: N=%d L=%.1f ghost=%d dx=%.4f dt=%.6f",d,dc->N,dc->L,ghost,g->dx,g->dt);
        if(dc->parent_id>=0)printf(" parent=%d sub=%d",dc->parent_id,doms[d].substeps);
        printf(" init=%s\n",dc->init_mode);
        init_domain(&doms[d],dc,&c);
        dom_gpu_alloc(&doms[d]);dom_gpu_upload(&doms[d]);
        for(int f=0;f<12;f++){doms[d].pbuf_t0[f]=NULL;doms[d].pbuf_t1[f]=NULL;}
    }

    /* Initial forces */
    for(int d=0;d<c.n_domains;d++){
        gpu_set_constants(&c,doms[d].grid);
        gpu_forces(&doms[d],&c,d==0);
    }
    cudaDeviceSynchronize();

    /* SFA */
    uint8_t sfa_dtype=(c.precision==0)?SFA_F16:(c.precision==1)?SFA_F32:SFA_F64;
    SFA *sfa=sfa_create(c.output,c.dom[0].N,c.dom[0].N,c.dom[0].N,c.dom[0].L,c.dom[0].L,c.dom[0].L,doms[0].grid->dt);
    for(int d=0;d<c.n_domains;d++){
        DomainConfig *dc=&c.dom[d];
        sfa_add_grid(sfa,d,(dc->parent_id<0)?0xFFFF:dc->parent_id,
            dc->N,dc->N,dc->N,dc->L,dc->L,dc->L,dc->cx,dc->cy,dc->cz,dc->ghost,0);
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
    printf("\nSFA: %s (GDEF: %d grids)\n\n",c.output,c.n_domains);

    /* t=0 snapshots */
    for(int d=0;d<c.n_domains;d++){dom_gpu_download(&doms[d]);snap_domain(sfa,&doms[d],0,c.precision);}

    FILE *fp=fopen(c.diag_file,"w");
    fprintf(fp,"t\tdomain\tE_total\tE_pot\tphi_max\ttheta_rms\n");

    double dt_c=doms[0].grid->dt;
    int n_steps=(int)(c.T/dt_c);
    int diag_every=(int)(c.diag_dt/dt_c);if(diag_every<1)diag_every=1;

    double E0=0;
    for(int d=0;d<c.n_domains;d++){
        dom_gpu_download(&doms[d]);
        double et,ep,pm,trms;domain_energy(doms[d].grid,&c,d==0,&et,&ep,&pm,&trms);
        E0+=et;fprintf(fp,"%.2f\t%d\t%.6e\t%.6e\t%.6e\t%.6e\n",0.0,d,et,ep,pm,trms);
    }
    printf("INIT: E_total=%.4e\n\n",E0);

    cudaEvent_t t_start,t_stop;
    cudaEventCreate(&t_start);cudaEventCreate(&t_stop);cudaEventRecord(t_start);

    for(int step=1;step<=n_steps;step++){
        double t=step*dt_c;
        step_multigrid_gpu(doms,c.n_domains,&c);

        for(int d=0;d<c.n_domains;d++){
            if(t>=doms[d].next_snap-0.5*dt_c){
                dom_gpu_download(&doms[d]);snap_domain(sfa,&doms[d],t,c.precision);
                doms[d].next_snap=t+c.dom[d].snap_dt;
            }
        }
        if(step%diag_every==0){
            double et_sum=0;
            for(int d=0;d<c.n_domains;d++){
                dom_gpu_download(&doms[d]);
                double et,ep,pm,trms;domain_energy(doms[d].grid,&c,d==0,&et,&ep,&pm,&trms);
                et_sum+=et;fprintf(fp,"%.2f\t%d\t%.6e\t%.6e\t%.6e\t%.6e\n",t,d,et,ep,pm,trms);
            }
            fflush(fp);
            if(step%(diag_every*10)==0){
                cudaEventRecord(t_stop);cudaEventSynchronize(t_stop);
                float ms;cudaEventElapsedTime(&ms,t_start,t_stop);
                printf("t=%7.1f E=%.3e (drift %+.3f%%) [%.0f%% %.1fs]\n",
                    t,et_sum,100*(et_sum-E0)/(fabs(E0)+1e-30),100.0*step/n_steps,ms/1000);
            }
        }
    }

    double t_final=n_steps*dt_c;
    for(int d=0;d<c.n_domains;d++){dom_gpu_download(&doms[d]);snap_domain(sfa,&doms[d],t_final,c.precision);}
    uint32_t nf=sfa->total_frames;sfa_close(sfa);

    cudaEventRecord(t_stop);cudaEventSynchronize(t_stop);
    float total_ms;cudaEventElapsedTime(&total_ms,t_start,t_stop);
    printf("\n=== COMPLETE (GPU) ===\nSFA: %s (%u frames)\n",c.output,nf);
    for(int d=0;d<c.n_domains;d++){
        dom_gpu_download(&doms[d]);
        double et,ep,pm,trms;domain_energy(doms[d].grid,&c,d==0,&et,&ep,&pm,&trms);
        printf("  Domain %d: E=%.3e Ep=%.1f phi=%.3f trms=%.3e\n",d,et,ep,pm,trms);
    }
    printf("Wall: %.1fs (%.1f min)\n",total_ms/1000,total_ms/60000);

    fclose(fp);
    for(int d=0;d<c.n_domains;d++){dom_gpu_free(&doms[d]);ggrid_free(doms[d].grid);
        for(int f=0;f<12;f++){free(doms[d].pbuf_t0[f]);free(doms[d].pbuf_t1[f]);}}
    cudaEventDestroy(t_start);cudaEventDestroy(t_stop);
    return 0;
}
