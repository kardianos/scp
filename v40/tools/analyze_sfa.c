/*  analyze_sfa.c — Compute stability metrics from an SFA simulation file
 *
 *  Reads each frame, computes energy decomposition, triple product coherence,
 *  theta_rms, aspect ratio, and outputs a per-frame JSON analysis plus a
 *  summary score.
 *
 *  Build: gcc -O3 -fopenmp -o analyze_sfa analyze_sfa.c -lzstd -lm
 *  Usage: ./analyze_sfa input.sfa [--json output.json]
 */

#define SFA_IMPLEMENTATION
#include "../../sfa/format/sfa.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>

#define PI 3.14159265358979323846

/* Physics defaults (can be overridden from KVMD) */
static double MU = -41.345, KAPPA = 50.0, MASS2 = 2.25;
static double MTHETA2 = 0.0, ETA = 0.5;

typedef struct {
    double time;
    double E_kin, E_grad, E_mass, E_pot, E_theta_kin, E_theta_grad, E_coupling;
    double E_total;
    double phi_max, P_max, P_int, theta_rms;
    double aspect;  /* I_max / I_min */
    double score;   /* composite stability metric */
    int alive;      /* 1 = coherent, 0 = dead */
} FrameMetrics;

static void compute_metrics(double *phi[3], double *theta[3],
    double *phi_vel[3], double *theta_vel[3],
    int N, double L, FrameMetrics *m) {
    double dx = 2.0*L/(N-1), dV = dx*dx*dx, idx1 = 1.0/(2.0*dx);
    int NN = N*N;
    long N3 = (long)N*N*N;

    double epk=0,etk=0,eg=0,em=0,ep=0,etg=0,ec=0;
    double pm=0,Pm=0,Pint=0,trms=0;

    #pragma omp parallel for reduction(+:epk,etk,eg,em,ep,etg,ec,Pint,trms) \
        reduction(max:pm,Pm) schedule(static)
    for (long idx=0; idx<N3; idx++) {
        int i=(int)(idx/(NN)),j=(int)((idx/N)%N),k=(int)(idx%N);
        int ip=(i+1)%N,im=(i-1+N)%N,jp=(j+1)%N,jm=(j-1+N)%N,kp=(k+1)%N,km=(k-1+N)%N;
        long nip=(long)ip*NN+j*N+k,nim=(long)im*NN+j*N+k;
        long njp=(long)i*NN+jp*N+k,njm=(long)i*NN+jm*N+k;
        long nkp=(long)i*NN+j*N+kp,nkm=(long)i*NN+j*N+km;

        double p0=phi[0][idx],p1=phi[1][idx],p2=phi[2][idx];
        for (int a=0;a<3;a++) {
            if (phi_vel[a]) epk+=0.5*phi_vel[a][idx]*phi_vel[a][idx]*dV;
            if (theta_vel[a]) etk+=0.5*theta_vel[a][idx]*theta_vel[a][idx]*dV;
            double gx=(phi[a][nip]-phi[a][nim])*idx1;
            double gy=(phi[a][njp]-phi[a][njm])*idx1;
            double gz=(phi[a][nkp]-phi[a][nkm])*idx1;
            eg+=0.5*(gx*gx+gy*gy+gz*gz)*dV;
            em+=0.5*MASS2*phi[a][idx]*phi[a][idx]*dV;
            if (theta[a]) {
                double tgx=(theta[a][nip]-theta[a][nim])*idx1;
                double tgy=(theta[a][njp]-theta[a][njm])*idx1;
                double tgz=(theta[a][nkp]-theta[a][nkm])*idx1;
                etg+=0.5*(tgx*tgx+tgy*tgy+tgz*tgz)*dV;
                trms+=theta[a][idx]*theta[a][idx];
            }
            double ap=fabs(phi[a][idx]); if(ap>pm) pm=ap;
        }
        double P=p0*p1*p2, P2=P*P;
        ep+=(MU/2.0)*P2/(1.0+KAPPA*P2)*dV;
        Pint+=fabs(P)*dV;
        double Pa=fabs(P); if(Pa>Pm) Pm=Pa;
        /* Coupling energy (approximate — need curl) */
        if (theta[0] && theta[1] && theta[2]) {
            double ct0=(theta[2][njp]-theta[2][njm]-theta[1][nkp]+theta[1][nkm])*idx1;
            double ct1=(theta[0][nkp]-theta[0][nkm]-theta[2][nip]+theta[2][nim])*idx1;
            double ct2=(theta[1][nip]-theta[1][nim]-theta[0][njp]+theta[0][njm])*idx1;
            ec -= ETA*(phi[0][idx]*ct0+phi[1][idx]*ct1+phi[2][idx]*ct2)*dV;
        }
    }

    m->E_kin=epk; m->E_grad=eg; m->E_mass=em; m->E_pot=ep;
    m->E_theta_kin=etk; m->E_theta_grad=etg; m->E_coupling=ec;
    m->E_total=epk+etk+eg+em+ep+etg+ec;
    m->phi_max=pm; m->P_max=Pm; m->P_int=Pint;
    m->theta_rms=sqrt(trms/(3.0*N3));

    /* Moment of inertia for aspect ratio */
    double cm[3]={0}, E_sum=0;
    for (long idx=0;idx<N3;idx++) {
        int i=(int)(idx/(NN)),j=(int)((idx/N)%N),k=(int)(idx%N);
        double x=-L+i*dx, y=-L+j*dx, z=-L+k*dx;
        double w=0; for(int a=0;a<3;a++) w+=phi[a][idx]*phi[a][idx]; w*=dV;
        cm[0]+=x*w; cm[1]+=y*w; cm[2]+=z*w; E_sum+=w;
    }
    if(E_sum>0){cm[0]/=E_sum;cm[1]/=E_sum;cm[2]/=E_sum;}
    double I[3]={0}; /* diagonal approximation */
    for (long idx=0;idx<N3;idx++) {
        int i=(int)(idx/(NN)),j=(int)((idx/N)%N),k=(int)(idx%N);
        double x=-L+i*dx-cm[0], y=-L+j*dx-cm[1], z=-L+k*dx-cm[2];
        double w=0; for(int a=0;a<3;a++) w+=phi[a][idx]*phi[a][idx]; w*=dV;
        I[0]+=(y*y+z*z)*w; I[1]+=(x*x+z*z)*w; I[2]+=(x*x+y*y)*w;
    }
    double Imax=I[0],Imin=I[0];
    for(int i=1;i<3;i++){if(I[i]>Imax)Imax=I[i];if(I[i]<Imin)Imin=I[i];}
    m->aspect = (Imin>1e-30) ? Imax/Imin : 999;

    /* Alive check */
    m->alive = (m->P_int > 0.1 && m->E_pot < 0 && m->phi_max > 0.2);
}

static double compute_score(FrameMetrics *m, double E_pot_braid, double P_int_0) {
    double energy_score = (E_pot_braid != 0) ? fmax(0, -m->E_pot / fabs(E_pot_braid)) : 0;
    double coherence_score = (P_int_0 > 0) ? m->P_int / P_int_0 : 0;
    double theta_score = fmin(1.0, m->theta_rms / 0.05);  /* 0.05 is approximate equilibrium */
    double shape_score = 1.0 / (1.0 + fmax(0, m->aspect - 1));
    m->score = 0.4*energy_score + 0.3*coherence_score + 0.2*theta_score + 0.1*shape_score;
    return m->score;
}

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "Usage: %s input.sfa [--json out.json]\n", argv[0]); return 1; }
    const char *sfa_path = argv[1];
    const char *json_path = NULL;
    for (int i=2;i<argc;i++) if(!strcmp(argv[i],"--json")&&i+1<argc) json_path=argv[++i];

    SFA *sfa = sfa_open(sfa_path);
    if (!sfa) { fprintf(stderr, "Cannot open %s\n", sfa_path); return 1; }

    int N = sfa->Nx;
    double L = sfa->Lx;
    long N3 = (long)N*N*N;

    /* Try to read physics from KVMD */
    SFA_KVMDSet kv[16];
    int n_kv = sfa_read_kvmd(sfa, kv, 16);
    for (int s=0;s<n_kv;s++) for(int p=0;p<kv[s].n_pairs;p++) {
        if(!strcmp(kv[s].keys[p],"mu")) MU=atof(kv[s].values[p]);
        if(!strcmp(kv[s].keys[p],"kappa")) KAPPA=atof(kv[s].values[p]);
        if(!strcmp(kv[s].keys[p],"m")){double m=atof(kv[s].values[p]);MASS2=m*m;}
        if(!strcmp(kv[s].keys[p],"eta")) ETA=atof(kv[s].values[p]);
    }

    printf("Analyzing: %s (N=%d, L=%.1f, %u frames)\n", sfa_path, N, L, sfa->total_frames);
    printf("Physics: mu=%.3f kappa=%.1f m^2=%.4f eta=%.3f\n\n", MU, KAPPA, MASS2, ETA);

    void *buf = malloc(sfa->frame_bytes);
    FrameMetrics *metrics = (FrameMetrics*)calloc(sfa->total_frames, sizeof(FrameMetrics));

    double P_int_0 = 0;

    for (uint32_t f = 0; f < sfa->total_frames; f++) {
        sfa_read_frame(sfa, f, buf);
        double t = sfa_frame_time(sfa, f);

        /* Extract fields by semantic */
        double *phi[3]={NULL}, *theta[3]={NULL}, *phi_vel[3]={NULL}, *theta_vel[3]={NULL};
        uint64_t off = 0;
        for (int c=0;c<sfa->n_columns;c++) {
            int dtype=sfa->columns[c].dtype, sem=sfa->columns[c].semantic, comp=sfa->columns[c].component;
            int es=sfa_dtype_size[dtype];
            /* Upconvert to double */
            double *arr = (double*)malloc(N3*sizeof(double));
            uint8_t *src = (uint8_t*)buf + off;
            if(dtype==SFA_F64) memcpy(arr,src,N3*8);
            else if(dtype==SFA_F32) for(long i=0;i<N3;i++) arr[i]=(double)((float*)src)[i];
            else if(dtype==SFA_F16) for(long i=0;i<N3;i++){
                uint16_t h=((uint16_t*)src)[i];int e=(h>>10)&0x1F;uint16_t m=h&0x3FF;
                uint16_t s=h&0x8000;float fv;
                if(e==0){arr[i]=0;}else{uint32_t x=((uint32_t)s<<16)|((uint32_t)(e-15+127)<<23)|((uint32_t)m<<13);
                memcpy(&fv,&x,4);arr[i]=(double)fv;}
            }
            if(sem==SFA_POSITION&&comp<3) phi[comp]=arr;
            else if(sem==SFA_ANGLE&&comp<3) theta[comp]=arr;
            else if(sem==SFA_VELOCITY&&comp<3) phi_vel[comp]=arr;
            else if(sem==SFA_VELOCITY&&comp>=3&&comp<6) theta_vel[comp-3]=arr;
            else free(arr);
            off += (uint64_t)N3*es;
        }

        metrics[f].time = t;
        compute_metrics(phi, theta, phi_vel, theta_vel, N, L, &metrics[f]);
        if (f == 0) P_int_0 = metrics[f].P_int;
        compute_score(&metrics[f], -80.0, P_int_0);

        printf("Frame %2u t=%6.1f: E=%.2e Ep=%7.1f Pint=%6.1f phi=%.3f θ=%.3e asp=%.2f S=%.3f %s\n",
               f, t, metrics[f].E_total, metrics[f].E_pot, metrics[f].P_int,
               metrics[f].phi_max, metrics[f].theta_rms, metrics[f].aspect,
               metrics[f].score, metrics[f].alive?"ALIVE":"DEAD");

        for(int a=0;a<3;a++){free(phi[a]);if(theta[a])free(theta[a]);
            if(phi_vel[a])free(phi_vel[a]);if(theta_vel[a])free(theta_vel[a]);}
    }

    /* Summary */
    int last = sfa->total_frames - 1;
    double S_final = metrics[last].score;
    double S_mean = 0;
    int alive_count = 0;
    for (uint32_t f=0;f<sfa->total_frames;f++) { S_mean+=metrics[f].score; alive_count+=metrics[f].alive; }
    S_mean /= sfa->total_frames;

    printf("\n=== SUMMARY ===\n");
    printf("S_final = %.4f\n", S_final);
    printf("S_mean  = %.4f\n", S_mean);
    printf("Alive frames: %d / %u\n", alive_count, sfa->total_frames);
    printf("E_pot(final) = %.1f\n", metrics[last].E_pot);
    printf("P_int retention = %.1f%%\n", 100.0*metrics[last].P_int/P_int_0);

    /* Write JSON if requested */
    if (json_path) {
        FILE *jfp = fopen(json_path, "w");
        fprintf(jfp, "{\n  \"sfa\": \"%s\",\n  \"N\": %d, \"L\": %.1f,\n", sfa_path, N, L);
        fprintf(jfp, "  \"S_final\": %.6f,\n  \"S_mean\": %.6f,\n", S_final, S_mean);
        fprintf(jfp, "  \"alive_frames\": %d,\n  \"total_frames\": %u,\n", alive_count, sfa->total_frames);
        fprintf(jfp, "  \"E_pot_final\": %.6f,\n  \"P_int_retention\": %.6f,\n",
                metrics[last].E_pot, metrics[last].P_int/P_int_0);
        fprintf(jfp, "  \"frames\": [\n");
        for (uint32_t f=0;f<sfa->total_frames;f++) {
            FrameMetrics *fm = &metrics[f];
            fprintf(jfp, "    {\"t\":%.2f,\"E\":%.4e,\"Ep\":%.4f,\"Pint\":%.4f,"
                   "\"phi_max\":%.4f,\"theta_rms\":%.4e,\"aspect\":%.3f,\"S\":%.4f,\"alive\":%d}%s\n",
                   fm->time,fm->E_total,fm->E_pot,fm->P_int,fm->phi_max,fm->theta_rms,
                   fm->aspect,fm->score,fm->alive, (f<sfa->total_frames-1)?",":"");
        }
        fprintf(jfp, "  ]\n}\n");
        fclose(jfp);
        printf("JSON written to %s\n", json_path);
    }

    free(buf); free(metrics);
    sfa_close(sfa);
    return 0;
}
