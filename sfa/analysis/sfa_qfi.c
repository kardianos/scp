/*  sfa_qfi.c — Quantum Fisher Information entanglement witness for complex SFA runs
 *
 *  CONCEPT.md §9 / FUTURE.md "Classical->Quantum Bridge". Read-only check that
 *  tests whether the classical SCP field reproduces a *quantum* many-body
 *  signature: the quantum Fisher information of Mazza et al. (Nat. Phys. 2026,
 *  doi:10.1038/s41567-026-03298-0), an entanglement witness from a dynamical
 *  structure factor S(q,omega).
 *
 *  KEY PHYSICS (why a single Q-ball is null, and what fixes it):
 *    A clean isolated matter ball factorizes -> static |Phi|^2 -> trivial witness.
 *    Genuine multipartite entanglement lives in the CONSTRAINT-INDUCED
 *    correlations between sectors:
 *      - gauge:     Gauss law  div E = g rho_Q  ties E to the matter charge.
 *      - geometric: Cosserat torsion lock  eta|curl(Phi)/2 - Theta|^2  ties the
 *                   Theta partner to the matter vorticity.
 *    So this tool reports, for each gauge-invariant scalar observable, the QFI
 *    witness f_Q, AND the matter<->gauge / matter<->torsion cross-coherence
 *    Gamma(q) (the correlation an entanglement witness detects).
 *
 *  Observables (gauge-invariant scalar densities, projected to q):
 *    rhoQ     : sum_a (u_a vdot_a - v_a udot_a)        charge density
 *    phiMod   : sum_a (u_a^2 + v_a^2)                  matter density
 *    thetaMod : sum_a (tu_a^2 + tv_a^2)                torsion-partner density   [if present]
 *    Emag     : E_x^2+E_y^2+E_z^2                      gauge electric density    [if present]
 *    Amag     : th_x^2+th_y^2+th_z^2                   gauge connection density  [if present]
 *    field0   : u_0 + i v_0  (complex)                 clock diagnostic (NOT a witness:
 *                                                      coherent carrier, gauge-dependent)
 *
 *  hbar is supplied by the soliton identity hbar_eff = Q (CONCEPT.md §7), measured.
 *  Operating point (remainder item 3): the kernel only discriminates when
 *  hbar*w ~ T. With --auto-T, T_eff is set per series to hbar*w_peak/2 so the
 *  kernel sits at O(1) at the dominant mode instead of saturating.
 *
 *  Build: gcc -O3 -fopenmp -o sfa_qfi sfa_qfi.c -lzstd -lm
 *  Usage: sfa_qfi input.sfa [options]
 *    --T <Teff>          fixed effective temperature (default 1.0)
 *    --auto-T            set T_eff = hbar*w_peak/2 per (obs,q)  [operating point]
 *    --hbar <H>          override hbar_eff (default: |Q| from frame 0)
 *    --hrange <d>        operator range (h_max-h_min) for nQFI (default 1)
 *    --q nx,ny,nz        add wavevector q = (pi/L)*(nx,ny,nz) (repeatable)
 *                        default scan: (0..4,0,0)
 *    --tsv out.tsv       write per-(obs,q) table
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define SFA_IMPLEMENTATION
#include "../format/sfa.h"

static int find_column(const SFA *sfa, const char *name) {
    for (uint32_t c = 0; c < sfa->n_columns; c++)
        if (strncmp(sfa->columns[c].name, name, sizeof(sfa->columns[c].name)) == 0)
            return (int)c;
    return -1;
}

static float *extract_column_f32(const void *buf, const SFA *sfa, int col_idx) {
    uint64_t N3 = sfa->N_total, off = 0;
    for (int c = 0; c < col_idx; c++)
        off += N3 * sfa_dtype_size[sfa->columns[c].dtype];
    int dtype = sfa->columns[col_idx].dtype;
    const uint8_t *src = (const uint8_t *)buf + off;
    float *arr = (float *)malloc(N3 * sizeof(float));
    if (!arr) { fprintf(stderr, "error: malloc\n"); exit(1); }
    if (dtype == SFA_F32) memcpy(arr, src, N3 * sizeof(float));
    else if (dtype == SFA_F64) {
        const double *d = (const double *)src;
        #pragma omp parallel for
        for (uint64_t i = 0; i < N3; i++) arr[i] = (float)d[i];
    } else if (dtype == SFA_F16) {
        const uint16_t *h = (const uint16_t *)src;
        #pragma omp parallel for
        for (uint64_t i = 0; i < N3; i++) arr[i] = sfa_f16_to_f32(h[i]);
    } else { fprintf(stderr, "error: dtype %d\n", dtype); exit(1); }
    return arr;
}

#define MAXQ 32
/* observable ids */
enum { O_RHOQ=0, O_PHIMOD, O_THETAMOD, O_EMAG, O_AMAG, O_FIELD0, NOBS };
static const char *OBS_NAME[NOBS] = {"rhoQ","phiMod","thetaMod","Emag","Amag","field0"};

int main(int argc, char **argv) {
    const char *in_path = NULL, *tsv_path = NULL;
    double Teff = 1.0, hbar_override = -1.0, hrange = 1.0;
    int autoT = 0;
    int qn[MAXQ][3]; int nq = 0;

    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "--T") && i+1 < argc) Teff = atof(argv[++i]);
        else if (!strcmp(argv[i], "--auto-T")) autoT = 1;
        else if (!strcmp(argv[i], "--hbar") && i+1 < argc) hbar_override = atof(argv[++i]);
        else if (!strcmp(argv[i], "--hrange") && i+1 < argc) hrange = atof(argv[++i]);
        else if (!strcmp(argv[i], "--tsv") && i+1 < argc) tsv_path = argv[++i];
        else if (!strcmp(argv[i], "--q") && i+1 < argc) {
            if (nq < MAXQ && sscanf(argv[++i], "%d,%d,%d",
                &qn[nq][0], &qn[nq][1], &qn[nq][2]) == 3) nq++;
        }
        else in_path = argv[i];
    }
    if (!in_path) { fprintf(stderr, "usage: %s input.sfa [--T t|--auto-T][--hbar H][--hrange d][--q nx,ny,nz][--tsv f]\n", argv[0]); return 1; }
    if (nq == 0) { for (int n=0;n<=4;n++){ qn[nq][0]=n;qn[nq][1]=0;qn[nq][2]=0;nq++; } }

    SFA *sfa = sfa_open(in_path);
    if (!sfa) { fprintf(stderr, "error: cannot open %s\n", in_path); return 1; }

    /* matter columns (required) */
    static const char *u_n[3]={"phi_x","phi_y","phi_z"}, *v_n[3]={"phiim_x","phiim_y","phiim_z"};
    static const char *ud_n[3]={"phi_vx","phi_vy","phi_vz"}, *vd_n[3]={"phiim_vx","phiim_vy","phiim_vz"};
    int cu[3],cv[3],cud[3],cvd[3];
    for (int a=0;a<3;a++){ cu[a]=find_column(sfa,u_n[a]);cv[a]=find_column(sfa,v_n[a]);
        cud[a]=find_column(sfa,ud_n[a]);cvd[a]=find_column(sfa,vd_n[a]);
        if(cu[a]<0||cv[a]<0||cud[a]<0||cvd[a]<0){fprintf(stderr,"error: needs complex_phi=1 output\n");return 1;} }
    /* optional torsion (theta) columns */
    static const char *tu_n[3]={"theta_x","theta_y","theta_z"}, *tv_n[3]={"thetaim_x","thetaim_y","thetaim_z"};
    int ctu[3],ctv[3], have_theta=1;
    for(int a=0;a<3;a++){ ctu[a]=find_column(sfa,tu_n[a]); ctv[a]=find_column(sfa,tv_n[a]); if(ctu[a]<0||ctv[a]<0)have_theta=0; }
    /* optional gauge columns */
    static const char *e_n[3]={"E_x","E_y","E_z"}, *a_n[3]={"th_x","th_y","th_z"};
    int ce[3],ca[3], have_gauge=1;
    for(int a=0;a<3;a++){ ce[a]=find_column(sfa,e_n[a]); ca[a]=find_column(sfa,a_n[a]); if(ce[a]<0||ca[a]<0)have_gauge=0; }

    int obs_on[NOBS]={1,1,have_theta,have_gauge,have_gauge,1};

    int N=(int)sfa->Nx; double L=sfa->Lx, dx=2.0*L/(N-1), dV=dx*dx*dx, Vbox=pow(2.0*L,3.0);
    long N3=(long)N*N*N, NN=(long)N*N; uint32_t M=sfa->total_frames;
    if (M<8){ fprintf(stderr,"error: need >=8 frames (have %u)\n",M); return 1; }

    double q0=M_PI/L;
    double *tc=malloc((size_t)nq*3*N*sizeof(double)), *ts=malloc((size_t)nq*3*N*sizeof(double));
    for(int qi=0;qi<nq;qi++)for(int ax=0;ax<3;ax++){ double qv=q0*qn[qi][ax];
        for(int idx=0;idx<N;idx++){ double c=-L+idx*dx; tc[((size_t)qi*3+ax)*N+idx]=cos(qv*c); ts[((size_t)qi*3+ax)*N+idx]=sin(qv*c);} }

    /* time series re/im[obs][qi][frame] */
    size_t SER=(size_t)NOBS*nq*M;
    double *re=calloc(SER,sizeof(double)), *im=calloc(SER,sizeof(double));
    double *tarr=malloc(M*sizeof(double)); double Qtot0=0;
    void *buf=malloc(sfa->frame_bytes);
    fprintf(stderr,"# reading %u frames N=%d L=%g  theta=%d gauge=%d\n",M,N,L,have_theta,have_gauge);

    for(uint32_t fi=0;fi<M;fi++){
        tarr[fi]=sfa_frame_time(sfa,fi);
        if(sfa_read_frame(sfa,fi,buf)!=0){ M=fi; break; }
        float *u[3],*v[3],*du[3],*dv[3]; for(int a=0;a<3;a++){u[a]=extract_column_f32(buf,sfa,cu[a]);v[a]=extract_column_f32(buf,sfa,cv[a]);du[a]=extract_column_f32(buf,sfa,cud[a]);dv[a]=extract_column_f32(buf,sfa,cvd[a]);}
        float *tu[3]={0},*tv[3]={0}; if(have_theta)for(int a=0;a<3;a++){tu[a]=extract_column_f32(buf,sfa,ctu[a]);tv[a]=extract_column_f32(buf,sfa,ctv[a]);}
        float *E[3]={0},*Av[3]={0}; if(have_gauge)for(int a=0;a<3;a++){E[a]=extract_column_f32(buf,sfa,ce[a]);Av[a]=extract_column_f32(buf,sfa,ca[a]);}
        if(fi==0){ double q=0;
            #pragma omp parallel for reduction(+:q)
            for(long i=0;i<N3;i++)for(int a=0;a<3;a++)q+=(double)u[a][i]*dv[a][i]-(double)v[a][i]*du[a][i];
            Qtot0=q*dV; }
        for(int qi=0;qi<nq;qi++){
            double *txc=&tc[((size_t)qi*3+0)*N],*txs=&ts[((size_t)qi*3+0)*N];
            double *tyc=&tc[((size_t)qi*3+1)*N],*tys=&ts[((size_t)qi*3+1)*N];
            double *tzc=&tc[((size_t)qi*3+2)*N],*tzs=&ts[((size_t)qi*3+2)*N];
            double rR=0,iR=0, rP=0,iP=0, rT=0,iT=0, rE=0,iE=0, rA=0,iA=0, rF=0,iF=0;
            #pragma omp parallel for reduction(+:rR,iR,rP,iP,rT,iT,rE,iE,rA,iA,rF,iF)
            for(long i=0;i<N3;i++){
                int ix=(int)(i/NN),iy=(int)((i/N)%N),iz=(int)(i%N);
                double ar=txc[ix],ai=-txs[ix], br=tyc[iy],bi=-tys[iy], cr2=tzc[iz],ci2=-tzs[iz];
                double abr=ar*br-ai*bi, abi=ar*bi+ai*br;
                double pr=abr*cr2-abi*ci2, pi=abr*ci2+abi*cr2;   /* e^{-iq.x} */
                double rho=0,pm=0;
                for(int a=0;a<3;a++){ double uu=u[a][i],vv=v[a][i];
                    rho+=uu*(double)dv[a][i]-vv*(double)du[a][i]; pm+=uu*uu+vv*vv; }
                rR+=rho*pr; iR+=rho*pi;  rP+=pm*pr; iP+=pm*pi;
                if(have_theta){ double tm=0; for(int a=0;a<3;a++){tm+=(double)tu[a][i]*tu[a][i]+(double)tv[a][i]*tv[a][i];} rT+=tm*pr; iT+=tm*pi; }
                if(have_gauge){ double em=(double)E[0][i]*E[0][i]+(double)E[1][i]*E[1][i]+(double)E[2][i]*E[2][i];
                                double am=(double)Av[0][i]*Av[0][i]+(double)Av[1][i]*Av[1][i]+(double)Av[2][i]*Av[2][i];
                                rE+=em*pr; iE+=em*pi; rA+=am*pr; iA+=am*pi; }
                double f0u=u[0][i],f0v=v[0][i]; rF+=f0u*pr-f0v*pi; iF+=f0u*pi+f0v*pr;
            }
            #define PUT(ob,RR,II) do{ size_t b=(((size_t)(ob))*nq+qi)*M+fi; re[b]=(RR)*dV; im[b]=(II)*dV; }while(0)
            PUT(O_RHOQ,rR,iR); PUT(O_PHIMOD,rP,iP);
            if(have_theta)PUT(O_THETAMOD,rT,iT);
            if(have_gauge){PUT(O_EMAG,rE,iE);PUT(O_AMAG,rA,iA);}
            PUT(O_FIELD0,rF,iF);
            #undef PUT
        }
        for(int a=0;a<3;a++){free(u[a]);free(v[a]);free(du[a]);free(dv[a]); if(have_theta){free(tu[a]);free(tv[a]);} if(have_gauge){free(E[a]);free(Av[a]);}}
    }
    free(buf); free(tc); free(ts);

    double dt=(tarr[M-1]-tarr[0])/(M-1), maxdev=0;
    for(uint32_t k=1;k<M;k++){double d=fabs((tarr[k]-tarr[k-1])-dt); if(d>maxdev)maxdev=d;}
    if(maxdev>1e-6*fabs(dt)) fprintf(stderr,"# WARNING: non-uniform spacing (dev %.3g vs dt=%.3g)\n",maxdev,dt);
    double hbar=(hbar_override>0)?hbar_override:fabs(Qtot0);
    double Tspan=M*dt, dom=2.0*M_PI/Tspan;

    /* spectra X[obs][qi][freq] (complex) */
    double *Xr=calloc(SER,sizeof(double)), *Xi=calloc(SER,sizeof(double));
    for(int ob=0;ob<NOBS;ob++){ if(!obs_on[ob])continue;
        for(int qi=0;qi<nq;qi++){
            double *xr=&re[((size_t)ob*nq+qi)*M], *xi=&im[((size_t)ob*nq+qi)*M];
            double mr=0,mi=0; for(uint32_t k=0;k<M;k++){mr+=xr[k];mi+=xi[k];} mr/=M;mi/=M;
            for(uint32_t j=0;j<M;j++){ int jj=(j<=M/2)?(int)j:(int)j-(int)M; double w=jj*dom; double Ar=0,Ai=0;
                for(uint32_t n=0;n<M;n++){ double ph=-w*(tarr[n]-tarr[0]); double drr=xr[n]-mr,dii=xi[n]-mi; double c=cos(ph),s=sin(ph); Ar+=drr*c-dii*s; Ai+=drr*s+dii*c; }
                size_t b=((size_t)ob*nq+qi)*M+j; Xr[b]=Ar*dt; Xi[b]=Ai*dt;
            }
        }
    }

    printf("# %s\n",in_path);
    printf("# M=%u dt=%.4f T_span=%.2f w_Nyq=%.4f d_w=%.4f  V=%.4g Q=%.4f hbar_eff=%.4f%s hrange=%.4g\n",
           M,dt,Tspan,M_PI/dt,dom,Vbox,Qtot0,hbar, autoT?" T=auto(hbar*w_pk/2)":"", hrange);
    printf("# f_Q = 4 sum_{w>0} tanh(hbar w/2T)(1-e^{-hbar w/T}) S dw ;  nQFI=f_Q/hrange^2\n#\n");
    printf("%-9s %-7s %9s %10s %9s %12s %12s\n","obs","q","w_peak","intS","T_eff","f_Q","nQFI");

    FILE *tsv=NULL; if(tsv_path){tsv=fopen(tsv_path,"w"); if(tsv)fprintf(tsv,"obs\tqx\tqy\tqz\tw_peak\tintS\tT_eff\tf_Q\tnQFI\n");}

    for(int ob=0;ob<NOBS;ob++){ if(!obs_on[ob])continue;
        for(int qi=0;qi<nq;qi++){
            double *xr=&Xr[((size_t)ob*nq+qi)*M], *xi=&Xi[((size_t)ob*nq+qi)*M];
            double intS=0,wpeak=0,Speak=-1;
            for(uint32_t j=0;j<M;j++){ int jj=(j<=M/2)?(int)j:(int)j-(int)M; double w=jj*dom; if(w<=0)continue;
                double S=(xr[j]*xr[j]+xi[j]*xi[j])/(Tspan*Vbox); intS+=S*dom; if(S>Speak){Speak=S;wpeak=w;} }
            double Tuse = autoT ? (hbar*wpeak/2.0>0?hbar*wpeak/2.0:Teff) : Teff;
            double fQ=0;
            for(uint32_t j=0;j<M;j++){ int jj=(j<=M/2)?(int)j:(int)j-(int)M; double w=jj*dom; if(w<=0)continue;
                double S=(xr[j]*xr[j]+xi[j]*xi[j])/(Tspan*Vbox);
                fQ+=4.0*tanh(hbar*w/(2.0*Tuse))*(1.0-exp(-hbar*w/Tuse))*S*dom; }
            double nQFI=fQ/(hrange*hrange);
            char ql[24]; snprintf(ql,sizeof(ql),"%d,%d,%d",qn[qi][0],qn[qi][1],qn[qi][2]);
            printf("%-9s %-7s %9.4f %10.3e %9.3f %12.5e %12.5e\n",OBS_NAME[ob],ql,wpeak,intS,Tuse,fQ,nQFI);
            if(tsv)fprintf(tsv,"%s\t%d\t%d\t%d\t%.6f\t%.6e\t%.4f\t%.6e\t%.6e\n",OBS_NAME[ob],qn[qi][0],qn[qi][1],qn[qi][2],wpeak,intS,Tuse,fQ,nQFI);
        }
    }

    /* cross-sector coherence: rhoQ <-> {phiMod,thetaMod,Emag} (constraint-induced correlation) */
    printf("#\n# cross-sector coherence Gamma(q) = sum|S_AB| / sqrt(sum S_AA sum S_BB), A=rhoQ  [0..1]\n");
    printf("# (matter<->gauge / matter<->torsion correlation; the structure a witness detects)\n");
    printf("%-18s %-7s %10s\n","pair","q","Gamma");
    int partners[3]={O_PHIMOD,O_THETAMOD,O_EMAG};
    for(int p=0;p<3;p++){ int ob=partners[p]; if(!obs_on[ob])continue;
        for(int qi=0;qi<nq;qi++){
            double *ar=&Xr[((size_t)O_RHOQ*nq+qi)*M], *ai=&Xi[((size_t)O_RHOQ*nq+qi)*M];
            double *br=&Xr[((size_t)ob*nq+qi)*M], *bi=&Xi[((size_t)ob*nq+qi)*M];
            double sAB=0,sAA=0,sBB=0;
            for(uint32_t j=0;j<M;j++){ int jj=(j<=M/2)?(int)j:(int)j-(int)M; double w=jj*dom; if(w<=0)continue;
                double cr=ar[j]*br[j]+ai[j]*bi[j], ci=ai[j]*br[j]-ar[j]*bi[j]; /* A conj? use A * conj(B) */
                /* S_AB = A * conj(B) : real=ar*br+ai*bi, imag=ai*br-ar*bi */
                sAB+=sqrt(cr*cr+ci*ci); sAA+=ar[j]*ar[j]+ai[j]*ai[j]; sBB+=br[j]*br[j]+bi[j]*bi[j]; }
            double G=(sAA>0&&sBB>0)?sAB/sqrt(sAA*sBB):0;
            char ql[24]; snprintf(ql,sizeof(ql),"%d,%d,%d",qn[qi][0],qn[qi][1],qn[qi][2]);
            char pr[18]; snprintf(pr,sizeof(pr),"rhoQ-%s",OBS_NAME[ob]);
            printf("%-18s %-7s %10.4f\n",pr,ql,G);
        }
    }

    if(tsv)fclose(tsv);
    free(re);free(im);free(Xr);free(Xi);free(tarr);
    sfa_close(sfa);
    return 0;
}
