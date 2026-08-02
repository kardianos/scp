/*
 * charge_pack.c — the packing question (user directive 2026-08-02):
 * "how much field can we densely pack into a single harmonic structure."
 *
 * FK realization (standing gate law, p=8, A=1): a cycle of N sites packed
 * with k units of INTEGER winding. Two regimes:
 *   - k small: the winding packs as k distinct kinks, E ~ k*E1, each core
 *     localized and the rest of the cycle at the vacuum — CHARGED packing;
 *   - k large: the cores overlap, no core is addressable, the state is a
 *     uniform phase climb — STRAIN packing.
 * The crossover k* is the charge-packing capacity of one cycle.
 *
 * 2026-08-02 REVISION — two defects repaired:
 *
 *  (D10) THE KINK COUNTER COUNTED BONDS, NOT KINKS. nkink counted bonds
 *        with step > 1 rad, which is 4 for a SINGLE kink (that is the core
 *        width) and 100 at k=25. The regime classifier built on nk >= k
 *        was therefore meaningless — it labelled k=25 "charged(kinks)"
 *        while the energy said the identity was already lost there. FIX:
 *        count core CLUSTERS (maximal runs of large-step bonds) and
 *        classify on the energy crossover, which is the real signal.
 *
 *  (D11) THE "10.8x CHEAPER" HEADLINE COMPARED AGAINST A NON-STATIONARY
 *        STATE. E_strain(k) was evaluated on the UNRELAXED uniform ramp.
 *        At small k that configuration is not a competing packing at all —
 *        it relaxes to the kink. Quoting E_strain(1)/E1 = 10.8 as "the
 *        cost of incompatible packing" has no content. What the uniform
 *        state IS good for is the crossover: it is the true relaxed state
 *        at large k, so the k where E_relax meets E_strain marks capacity.
 *        FIX: E_strain is kept only as the crossover reference and is
 *        labelled as such; the low-k ratio is not quoted as a result.
 *        The meaningful cost comparison is merged-jump vs relaxed lattice,
 *        both stationary, which is the modest factor printed below.
 *
 * Build: gcc -O2 -o charge_pack charge_pack.c -lm
 */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static int P=8;
static double A=1.0;
static double Vgate(double th){double b=0.5*(1.0+cos(th));return 1.0-pow(b,P);}
static double dVgate(double th){double b=0.5*(1.0+cos(th));return 0.5*P*sin(th)*pow(b,P-1.0);}

static int N; static double TWIST;

static double chain_E(const double*th){
    double E=0;
    for(int i=0;i<N;i++)E+=Vgate(th[i]);
    for(int i=0;i<N-1;i++){double d=th[i+1]-th[i];E+=0.5*A*d*d;}
    double d=th[0]+TWIST-th[N-1];E+=0.5*A*d*d;
    return E;
}
static double relax(double*th,int iters,double lr,double*maxg){
    double*g=malloc(N*sizeof(double));
    double mf=0;
    for(int it=0;it<iters;it++){
        for(int i=0;i<N;i++){
            g[i]=dVgate(th[i]);
            if(i>0)g[i]+=A*(th[i]-th[i-1]);
            if(i<N-1)g[i]+=A*(th[i]-th[i+1]);
            if(i==0)g[i]+=A*(th[0]+TWIST-th[N-1]);
            if(i==N-1)g[i]-=A*(th[0]+TWIST-th[N-1]);
        }
        mf=0;
        for(int i=0;i<N;i++){th[i]-=lr*g[i];double f=fabs(g[i]);if(f>mf)mf=f;}
        if(mf<1e-14)break;
    }
    if(maxg)*maxg=mf;
    free(g);
    return chain_E(th);
}

/* Count CORE CLUSTERS, not large bonds (D10). A core is a maximal run of
 * consecutive bonds whose step exceeds half the maximum step; the number of
 * such runs around the ring is the number of resolvable kinks. */
static void diag(const double*th,double*maxstep,int*ncore){
    double steps[4096];
    for(int i=0;i<N-1;i++) steps[i]=fabs(th[i+1]-th[i]);
    steps[N-1]=fabs(th[0]+TWIST-th[N-1]);
    *maxstep=0;
    for(int i=0;i<N;i++) if(steps[i]>*maxstep)*maxstep=steps[i];
    double thr=0.5*(*maxstep);
    int runs=0;
    for(int i=0;i<N;i++){
        int prev=(i+N-1)%N;
        if(steps[i]>thr && !(steps[prev]>thr)) runs++;
    }
    /* a fully uniform profile has every bond above threshold and no
     * boundary anywhere: that is zero resolvable cores, not one. */
    int all=1; for(int i=0;i<N;i++) if(!(steps[i]>thr)) all=0;
    *ncore = all ? 0 : runs;
}
/* the uniform-strain state: exact stationary state of the ring at any k */
static double E_uniform(int k){
    double d=2*M_PI*k/N,E=0;
    for(int i=0;i<N;i++)E+=Vgate(i*d);
    E+=0.5*A*N*d*d;
    return E;
}

int main(void){
    printf("# charge_pack — the packing capacity of one harmonic cycle (p=8,A=1)\n");
    N=100;
    TWIST=2*M_PI;
    double*th=calloc(N,sizeof(double));
    /* site-centred seed so the reference is the MINIMUM, not the PN saddle
     * (charge_fk.c D4) */
    for(int i=0;i<N;i++) th[i]=(i<N/2)?0.0:2*M_PI;
    th[N/2]=M_PI;
    double mg=0;
    double E1=relax(th,400000,0.02,&mg);
    printf("# N=%d, E_core(k=1)=%.6f (max|grad|=%.1e)\n\n",N,E1,mg);

    printf("== winding packing scan: k integer charges on one N=%d cycle ==\n",N);
    printf("%-5s %-7s %-12s %-10s %-9s %-7s %-12s %s\n",
           "k","N/k","E_relax","E/(k*E1)","maxstep","ncore","E_uniform","regime");
    int ks[]={1,2,3,4,5,6,8,10,12,15,18,20,22,25,28,30,35,40,45,50};
    int kstar=0;
    for(int ki=0;ki<20;ki++){
        int k=ks[ki];
        TWIST=2*M_PI*k;
        for(int i=0;i<N;i++)th[i]=2*M_PI*k*(double)i/N;
        double E=relax(th,400000,0.02,&mg);
        double ms;int nc;
        diag(th,&ms,&nc);
        double Eu=E_uniform(k);
        /* capacity: the relaxed lattice has become the uniform state */
        int uniform = fabs(E-Eu) < 1e-3*Eu;
        if(uniform && !kstar) kstar=k;
        const char*reg = uniform ? "STRAIN (= the uniform state)"
                       : (nc>=k   ? "charged (k resolvable cores)"
                       : (nc>0    ? "transitional (cores merging)"
                                  : "no resolvable cores"));
        printf("%-5d %-7.2f %-12.4f %-10.4f %-9.4f %-7d %-12.4f %s\n",
               k,(double)N/k,E,E/(k*E1),ms,nc,Eu,reg);
    }
    printf("\ncapacity k* (first k where the relaxed lattice IS the uniform\n");
    printf("state, to 0.1%%): k* = %d, i.e. N/k* = %.1f sites per charge.\n",
           kstar, kstar?(double)N/kstar:0.0);
    printf("Compare charge_fk.c's independent floor from the energy excess on\n");
    printf("small cycles: N/k >= 6 costs <0.3%%, N/k = 4 costs 6.6%%. The two\n");
    printf("agree in kind (a few sites per charge) but NOT in value; the\n");
    printf("capacity is a crossover, not a sharp edge, and no single 'xi'\n");
    printf("reproduces both. Do not quote one number for it.\n");

    printf("\n== merged jump vs relaxed lattice at fixed winding W=8 ==\n");
    printf("(both stationary states of the same ring — a fair comparison)\n");
    int W=8; TWIST=2*M_PI*W;
    for(int i=0;i<N;i++)th[i]=0;
    double E_merged=relax(th,400000,0.02,&mg);
    double ms_m;int nc_m; diag(th,&ms_m,&nc_m);
    for(int i=0;i<N;i++)th[i]=2*M_PI*W*(double)i/N;
    double E_lattice=relax(th,400000,0.02,&mg);
    double ms;int nc; diag(th,&ms,&nc);
    printf("merged (all twist on one bond): E=%.4f maxstep=%.4f ncore=%d\n",E_merged,ms_m,nc_m);
    printf("lattice (relaxed):              E=%.4f maxstep=%.4f ncore=%d\n",E_lattice,ms,nc);
    printf("reading: kink repulsion (charge_fk2 N1) spreads the winding into an\n");
    printf("even lattice; the merged defect costs %.2fx more. That factor — not\n",
           E_merged/E_lattice);
    printf("any comparison against an unrelaxed ramp — is the cost of packing\n");
    printf("the same winding incompatibly.\n");

    printf("\n# charge_pack DONE\n");
    return 0;
}
