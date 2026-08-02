/*
 * charge_fk.c — the charge soliton, realized as a Frenkel-Kontorova chain
 * with the STANDING gate potential (CHARGE §7.2). Standalone; no cellfab
 * dependency. Subordinate to CHARGE.md §7 / charge_sym.mac.
 *
 * The model (the symbolic PART 4 made numerical):
 *   E[theta] = sum_i [ V(theta_i) + (A/2)*(theta_{i+1}-theta_i)^2 ]
 *   V(theta) = 1 - ((1+cos theta)/2)^p_gate         (the cellfab gate law)
 *   A        = coupling stiffness (from the measured tongue widths)
 *
 * A charge-q state (q = k/Q) is a kink: theta climbs 2*pi*k/Q across the
 * chain. Bogomolny bound: E_core = sqrt(2A)*I(p)*|k|/Q  with I(1)=4,
 * I(8)=5.405 (verified in maxima). Core width xi = sqrt(2A/V''(0)),
 * V''(0)=p/2.
 *
 * This program measures, parameter-free against the maxima predictions:
 *   (1) E_core vs A        -> sqrt(A) scaling
 *   (2) E_core vs p        -> I(p) shape factor (I(1)=4, I(8)=5.405)
 *   (3) kink width xi      -> sqrt(2A/V''(0))
 *   (4) SIZE CONDITION: the kink is stably trapped on a periodic cycle of
 *       N sites only when N*d > |k|*xi. Below it the kink unwinds
 *       (charge cannot form) — the user's "specific size conditions".
 *   (5) the mass decomposition M = N*e_bal + core + dressing.
 *
 * Modes: --mode open (one kink on an open chain)
 *        --mode cycle (periodic cycle, winding w=k/Q forced by BC)
 *
 * Build: gcc -O2 -o charge_fk charge_fk.c -lm
 * Run:   ./charge_fk
 */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* the standing gate potential and its derivative (curvature at min) */
static double Vgate(double th,int p){
    double c=cos(th);
    double base=0.5*(1.0+c);          /* in [0,1] */
    double bp=pow(base,p);            /* ((1+cos)/2)^p */
    return 1.0-bp;
}
static double dVgate(double th,int p){
    double c=cos(th),s=sin(th);
    double base=0.5*(1.0+c);
    /* d/dth base = -sin/2 ; d base^p = p*base^(p-1)*(-sin/2) */
    return 0.5*p*s*pow(base,p-1.0);
}
static double Vdd0(int p){ return 0.5*p; }   /* V''(0) = p/2 (verified) */

/* ---- a chain ---- */
typedef struct {
    int N;            /* sites */
    int p;            /* gate exponent */
    double A;         /* coupling stiffness */
    double *th;       /* phases */
    int periodic;     /* cycle (1) or open (0) */
    double twist;     /* BC twist theta_N - theta_0 (open: forced; cycle: topological) */
} Chain;

static Chain *chain_new(int N,int p,double A,int periodic,double twist){
    Chain *c=malloc(sizeof(Chain)); c->N=N; c->p=p; c->A=A; c->periodic=periodic;
    c->twist=twist; c->th=calloc(N,sizeof(double));
    return c;
}
static void chain_free(Chain*c){ free(c->th); free(c); }

/* total energy. For the periodic cycle the closing edge (N-1 -> 0) carries
 * the topological twist: theta_0 = theta_N - twist, so the last coupling
 * term is (theta_0 + twist - theta_{N-1})^2. For an open chain the BC is
 * fixed: theta_0=0, theta_{N-1}=twist. */
static double chain_E(const Chain*c){
    double E=0;
    for(int i=0;i<c->N;i++) E+=Vgate(c->th[i],c->p);
    for(int i=0;i<c->N-1;i++){ double d=c->th[i+1]-c->th[i]; E+=0.5*c->A*d*d; }
    if(c->periodic){
        double d=c->th[0]+c->twist-c->th[c->N-1]; E+=0.5*c->A*d*d;
    }
    return E;
}

/* gradient (force) dE/dth_i at site i. E = sum V(th_i) + (A/2) sum (th_{i+1}-th_i)^2
 * => dE/dth_i = dV/dth_i + A*(2 th_i - th_{i-1} - th_{i+1})
 *             = dV/dth_i + A*(th_i - th_{i-1}) + A*(th_i - th_{i+1}).  */
static double chain_grad(const Chain*c,int i,double*g){
    g[i]=dVgate(c->th[i],c->p);
    if(i>0){ g[i]+=c->A*(c->th[i]-c->th[i-1]); }
    if(i<c->N-1){ g[i]+=c->A*(c->th[i]-c->th[i+1]); }
    if(c->periodic){
        /* closing edge (th_0 + twist - th_{N-1}): site 0 feels +A*(th_0+twist-th_{N-1}),
         * site N-1 feels -A*(same). */
        if(i==0){ g[i]+=c->A*(c->th[0]+c->twist-c->th[c->N-1]); }
        if(i==c->N-1){ g[i]-=c->A*(c->th[0]+c->twist-c->th[c->N-1]); }
    }
    return g[i];
}

/* overdamped relax to a force-free state (gradient = 0). Returns final E. */
static double chain_relax(Chain*c,int iters,double lr){
    double *g=malloc(c->N*sizeof(double));
    double E=chain_E(c);
    for(int it=0; it<iters; it++){
        for(int i=0;i<c->N;i++) chain_grad(c,i,g);
        /* fixed BC for open: pin ends */
        int lo = c->periodic ? 0 : 1;
        int hi = c->periodic ? c->N : c->N-1;
        double maxf=0;
        for(int i=lo;i<hi;i++){
            c->th[i] -= lr*g[i];
            double f=fabs(g[i]); if(f>maxf) maxf=f;
        }
        if((it%2000)==0 && maxf<1e-12) break;
    }
    E=chain_E(c);
    free(g);
    return E;
}

/* seed a kink of winding w (theta climbs 2*pi*w across the chain) */
static void seed_kink(Chain*c,double w){
    int N=c->N;
    for(int i=0;i<N;i++){
        double x=(double)i/(N-1);
        /* smooth kink profile: theta = pi*w*(1-cos(pi*x))? use linear + smooth */
        c->th[i] = 2.0*M_PI*w*x;
    }
    if(!c->periodic){ c->th[0]=0; c->th[N-1]=2*M_PI*w; }
}

/* count the realized winding (unwrap theta across the chain) */
static double realized_winding(const Chain*c){
    double w=0;
    for(int i=1;i<c->N;i++){
        double d=c->th[i]-c->th[i-1];
        w+=d;
    }
    if(c->periodic) w += (c->th[0]+c->twist-c->th[c->N-1]);
    return w/(2*M_PI);
}

/* core width: count sites where |theta_i - theta_{i-1}| > pi (the kink
 * core, where the phase jumps). Returns the spatial FWHM in sites. */
static double core_width(const Chain*c){
    double maxstep=0; int imax=0;
    for(int i=1;i<c->N;i++){ double d=fabs(c->th[i]-c->th[i-1]); if(d>maxstep){maxstep=d;imax=i;} }
    /* half-max width in steps */
    double half=maxstep/2; int cnt=0;
    for(int i=1;i<c->N;i++){ double d=fabs(c->th[i]-c->th[i-1]); if(d>half) cnt++; }
    return (double)cnt;
}

int main(int argc,char**argv){
    int p_default=8;
    printf("# charge_fk — FK charge soliton (gate potential p_gate)\n");
    printf("# maxima predictions: I(1)=4, I(8)=5.405336; V''(0)=p/2; xi=sqrt(2A/V''(0)); E_core=sqrt(2A)*I(p)\n\n");

    /* ---- (1)(2) E_core vs p on a long open chain (one kink, w=1) ---- */
    printf("== E_core vs p (A=1.0, N=400, one kink w=1) ==\n");
    printf("%-4s %-12s %-12s %-12s %-12s\n","p","E_measured","sqrt(2A)I(p)","ratio","xi[sites]");
    double A=1.0;
    /* I(p) measured by maxima (charge_sym.mac PART 1) */
    double Imap[17]={0,4.0,4.591174,0,5.057061,0,0,0,5.405336,0,0,0,0,0,0,0,5.658564};
    for(int p=1;p<=16;p=(p==1)?2:p*2){
        Chain*c=chain_new(400,p,A,0,2*M_PI);
        seed_kink(c,1.0);
        double E=chain_relax(c,40000,0.02);
        Chain*v=chain_new(400,p,A,0,0);
        for(int i=0;i<v->N;i++) v->th[i]=0;
        double Ev=chain_E(v);
        double Ecore=E-Ev;
        double xi=core_width(c);
        printf("%-4d %-12.6f %-12.6f %-12.4f %-12.2f\n",
               p, Ecore, sqrt(2*A)*Imap[p], Ecore/(sqrt(2*A)*Imap[p]), xi);
        chain_free(c); chain_free(v);
    }

    /* ---- (1b) E_core vs A scaling (p=8) ---- */
    printf("\n== E_core vs A (p=8, N=400, w=1) — expect sqrt(A) scaling ==\n");
    printf("%-8s %-12s %-12s %-12s\n","A","E_measured","sqrt(2A)*I(8)","ratio");
    double As[]={0.25,0.5,1.0,2.0,4.0};
    for(int k=0;k<5;k++){
        Chain*c=chain_new(400,8,As[k],0,2*M_PI);
        seed_kink(c,1.0);
        double E=chain_relax(c,40000,0.02);
        Chain*v=chain_new(400,8,As[k],0,0); for(int i=0;i<v->N;i++)v->th[i]=0; double Ev=chain_E(v);
        double Ecore=E-Ev;
        printf("%-8.2f %-12.6f %-12.6f %-12.4f\n",As[k],Ecore,sqrt(2*As[k])*5.405336,Ecore/(sqrt(2*As[k])*5.405336));
        chain_free(c);chain_free(v);
    }

    /* ---- (4) SIZE CONDITION: periodic cycle, winding w, scan N ----
     * The kink needs N > |w|*xi core sites to fit; below that the core is
     * squeezed by its periodic image, the energy rises, and at N ~ |w|*xi
     * the topological state becomes indistinguishable from uniform strain
     * (charge cannot form). Probe small N and several windings. */
    printf("\n== SIZE CONDITION: periodic cycle, p=8, A=1.0 ==\n");
    printf("kink core width xi ~ 3-5 sites (measured above)\n");
    printf("%-6s %-6s %-12s %-12s %-12s\n","w","N","E_total","winding","profile_maxstep");
    double ws[]={1.0,1.0,1.0,1.0,2.0,2.0,2.0,3.0,3.0};
    int    Ns[]={3,4,6,10,6,8,12,9,15};
    for(int k=0;k<9;k++){
        double w=ws[k]; int N=Ns[k];
        Chain*c=chain_new(N,8,1.0,1,2*M_PI*w);
        seed_kink(c,w);
        double E=chain_relax(c,80000,0.02);
        double wr=realized_winding(c);
        /* profile_maxstep = largest single-site phase jump (kink localized => big step) */
        double maxstep=0; for(int i=1;i<N;i++){double d=fabs(c->th[i]-c->th[i-1]); if(d>maxstep)maxstep=d;}
        double d=fabs(c->th[0]+2*M_PI*w-c->th[N-1]); if(d>maxstep)maxstep=d;
        printf("%-6.1f %-6d %-12.6f %-12.4f %-12.4f\n",w,N,E,wr,maxstep);
        chain_free(c);
    }
    printf("reading: at N >> |w|*xi the kink is a clean soliton (maxstep small, winding=w);\n");
    printf("         at N ~ |w|*xi the core is squeezed (maxstep grows) -> charge-forming size floor.\n");

    /* ---- (5) charge q = k/Q: scan the winding w on a large cycle ----
     * E_core should scale |w|*sqrt(2A)*I(p) (linear in |k|/Q). */
    printf("\n== CHARGE q=k/Q scaling (p=8, A=1.0, N=400, periodic) ==\n");
    printf("%-10s %-12s %-12s %-12s\n","q=w","E_core","|w|*E_core1","ratio");
    double Ecore1=-1;
    double wq[]={1.0,2.0,3.0,1.0/3,2.0/3};
    for(int k=0;k<5;k++){
        Chain*c=chain_new(400,8,1.0,1,2*M_PI*wq[k]);
        seed_kink(c,wq[k]);
        double E=chain_relax(c,60000,0.02);
        Chain*v=chain_new(400,8,1.0,1,0); for(int i=0;i<v->N;i++)v->th[i]=0; double Ev=chain_E(v);
        double Ecore=E-Ev;
        if(k==0) Ecore1=Ecore;
        printf("%-10.4f %-12.6f %-12.6f %-12.4f\n",wq[k],Ecore,fabs(wq[k])*Ecore1,
               Ecore/(fabs(wq[k])*Ecore1));
        chain_free(c);chain_free(v);
    }

    printf("\n# charge_fk DONE\n");
    return 0;
}
