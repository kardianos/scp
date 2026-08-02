/*
 * charge_fk.c — the charge soliton, realized as a Frenkel-Kontorova chain
 * with the STANDING gate potential (CHARGE §7.2). Standalone; no cellfab
 * dependency. Subordinate to CHARGE.md §7/§8 and charge_sym.mac.
 *
 * The model:
 *   E[theta] = sum_i [ V(theta_i) + (A/2)*(theta_{i+1}-theta_i)^2 ]
 *   V(theta) = 1 - ((1+cos theta)/2)^p_gate         (the cellfab gate law)
 *   A        = coupling stiffness (from the measured tongue widths)
 *
 * An INTEGER charge k is a kink: theta climbs 2*pi*k across the chain.
 * Continuum Bogomolny bound: E_core = sqrt(2A)*I(p)*|k|, I(1)=4,
 * I(8)=5.405336 (charge_sym.mac). Core width xi = sqrt(2A/V''(0)),
 * V''(0) = p/2.
 *
 * 2026-08-02 REVISION — two defects of the 2026-08-01 version repaired:
 *
 *  (D4) THE RELAXATION LANDED ON THE PEIERLS-NABARRO SADDLE, NOT THE
 *       MINIMUM. A single uniform-ramp seed on an even-length chain is
 *       symmetric about a BOND, and gradient descent preserves that
 *       symmetry: it converges to the bond-centred stationary state
 *       (E = 7.598210 at p=8, A=1, N=400; max|grad| ~ 9e-15, so it is a
 *       genuine stationary point -- just not the ground state). The
 *       minimum is the site-centred kink, E = 7.432761, 2.2% lower, and
 *       it is the value charge_fk2.c and charge_pack.c report for the
 *       same object. Every ratio in the old E_core tables was high for
 *       this reason (p=8 read 0.994 of the bound; it is really 0.972).
 *       FIX: relax from several seeds and keep the lowest (relax_best).
 *
 *  (D5) THE FRACTIONAL-CHARGE SCAN MEASURED AN ARTIFACT. A twist of
 *       2*pi*k/Q on a UNISON ring is not a topological state at all --
 *       2*pi/Q is not a period of V -- so the relaxed configuration is
 *       not a soliton but an elastic strain sitting on the ring's
 *       closing bond, with every other site at the vacuum. The old
 *       "integer linear / fractional sub-linear = the quark pattern"
 *       reading was comparing a soliton against a boundary condition.
 *       FIX: section (5) is now an explicit NEGATIVE CONTROL that prints
 *       the diagnostic (seam-bond load, non-vacuum site count) so the
 *       artifact is self-documenting. The 1/Q core scaling of §7.2 is
 *       untested and must be tested on a Z_Q INTERVAL chain instead.
 *
 * Measures, against the maxima predictions:
 *   (1)(2) E_core vs p, vs A -> I(p) shape factor and sqrt(A) scaling
 *   (3)    kink width xi
 *   (4)    SIZE CONDITION: winding k on a cycle of N voices
 *   (5)    NEGATIVE CONTROL: fractional winding on a unison ring
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

/* the standing gate potential and its derivative */
static double Vgate(double th,int p){
    double base=0.5*(1.0+cos(th));
    return 1.0-pow(base,p);
}
static double dVgate(double th,int p){
    double base=0.5*(1.0+cos(th));
    return 0.5*p*sin(th)*pow(base,p-1.0);
}
static double Vdd0(int p){ return 0.5*p; }   /* V''(0) = p/2 (verified) */

/* ---- a chain ---- */
typedef struct {
    int N;            /* sites */
    int p;            /* gate exponent */
    double A;         /* coupling stiffness */
    double *th;       /* phases */
    int periodic;     /* cycle (1) or open (0) */
    double twist;     /* BC twist */
} Chain;

static Chain *chain_new(int N,int p,double A,int periodic,double twist){
    Chain *c=malloc(sizeof(Chain)); c->N=N; c->p=p; c->A=A; c->periodic=periodic;
    c->twist=twist; c->th=calloc(N,sizeof(double));
    return c;
}
static void chain_free(Chain*c){ free(c->th); free(c); }

/* total energy. Periodic: the closing edge carries the topological twist,
 * theta_N == theta_0 + twist. Open: the ends are pinned at 0 and twist. */
static double chain_E(const Chain*c){
    double E=0;
    for(int i=0;i<c->N;i++) E+=Vgate(c->th[i],c->p);
    for(int i=0;i<c->N-1;i++){ double d=c->th[i+1]-c->th[i]; E+=0.5*c->A*d*d; }
    if(c->periodic){
        double d=c->th[0]+c->twist-c->th[c->N-1]; E+=0.5*c->A*d*d;
    }
    return E;
}

static void chain_grad(const Chain*c,double*g){
    int N=c->N;
    for(int i=0;i<N;i++){
        g[i]=dVgate(c->th[i],c->p);
        if(i>0)   g[i]+=c->A*(c->th[i]-c->th[i-1]);
        if(i<N-1) g[i]+=c->A*(c->th[i]-c->th[i+1]);
        if(c->periodic){
            if(i==0)     g[i]+=c->A*(c->th[0]+c->twist-c->th[N-1]);
            if(i==N-1)   g[i]-=c->A*(c->th[0]+c->twist-c->th[N-1]);
        }
    }
}
/* overdamped relax to a stationary state. Returns final E and, if maxg is
 * non-NULL, the residual max|grad| over the moving sites (the convergence
 * witness -- print it, do not trust a relaxation that does not report it). */
static double chain_relax(Chain*c,int iters,double lr,double*maxg){
    double *g=malloc(c->N*sizeof(double));
    int lo = c->periodic ? 0 : 1;
    int hi = c->periodic ? c->N : c->N-1;
    double mf=0;
    for(int it=0; it<iters; it++){
        chain_grad(c,g);
        mf=0;
        for(int i=lo;i<hi;i++){
            c->th[i] -= lr*g[i];
            double f=fabs(g[i]); if(f>mf) mf=f;
        }
        if(mf<1e-14) break;
    }
    if(maxg) *maxg=mf;
    double E=chain_E(c);
    free(g);
    return E;
}

/* --- seeds --- */
/* (a) uniform ramp: the whole twist spread evenly. */
static void seed_ramp(Chain*c,double w){
    int N=c->N;
    double denom = c->periodic ? (double)N : (double)(N-1);
    for(int i=0;i<N;i++) c->th[i]=2.0*M_PI*w*(double)i/denom;
    if(!c->periodic){ c->th[0]=0; c->th[N-1]=2*M_PI*w; }
}
/* (b) |k| sharp steps, each with one site parked at the potential maximum
 * (the site-centred configuration -- the PN MINIMUM for this potential). */
static void seed_steps(Chain*c,int k,int shift){
    int N=c->N;
    for(int i=0;i<N;i++){
        int units=0;
        for(int m=0;m<k;m++){
            int centre=(int)(((double)m+0.5)*N/k)+shift;
            if(i>centre) units++;
        }
        c->th[i]=2.0*M_PI*units;
        for(int m=0;m<k;m++){
            int centre=(int)(((double)m+0.5)*N/k)+shift;
            if(i==centre) c->th[i]+=M_PI;
        }
    }
    if(!c->periodic){ c->th[0]=0; c->th[N-1]=2*M_PI*k; }
}

/* Relax from every seed and keep the LOWEST stationary state. This is the
 * D4 fix: one seed selects one branch, and on a symmetric chain that
 * branch is the saddle. */
static double relax_best(Chain*c,double w,int iters,double lr,double*maxg){
    int N=c->N;
    double *best=malloc(N*sizeof(double));
    double Ebest=1e30, mgbest=0;
    int intw = (fabs(w-floor(w+0.5))<1e-9) ? (int)floor(fabs(w)+0.5) : 0;
    for(int s=0;s<3;s++){
        double mg=0;
        if(s==0)                seed_ramp(c,w);
        else if(intw>0)         seed_steps(c,intw,(s==1)?0:1);
        else                    continue;
        double E=chain_relax(c,iters,lr,&mg);
        if(E<Ebest){ Ebest=E; mgbest=mg; memcpy(best,c->th,N*sizeof(double)); }
    }
    memcpy(c->th,best,N*sizeof(double));
    free(best);
    if(maxg) *maxg=mgbest;
    return Ebest;
}

/* largest single-bond phase step (localised core => large step) */
static double max_step(const Chain*c){
    double mx=0;
    for(int i=1;i<c->N;i++){ double d=fabs(c->th[i]-c->th[i-1]); if(d>mx)mx=d; }
    if(c->periodic){
        double d=fabs(c->th[0]+c->twist-c->th[c->N-1]); if(d>mx)mx=d;
    }
    return mx;
}
/* FWHM of the step profile, in sites */
static double core_width(const Chain*c){
    double mx=max_step(c); int cnt=0;
    for(int i=1;i<c->N;i++){ double d=fabs(c->th[i]-c->th[i-1]); if(d>mx/2) cnt++; }
    return (double)cnt;
}
/* how many sites sit away from a vacuum minimum of V */
static int nonvacuum(const Chain*c){
    int n=0;
    for(int i=0;i<c->N;i++) if(Vgate(c->th[i],c->p)>0.01) n++;
    return n;
}

int main(void){
    printf("# charge_fk — FK charge soliton (gate potential p_gate)\n");
    printf("# maxima: I(1)=4, I(8)=5.405336; V''(0)=p/2; xi=sqrt(2A/V''(0));\n");
    printf("# continuum bound E_core = sqrt(2A)*I(p)*|k|\n");
    printf("# all states relaxed from multiple seeds; lowest kept (see D4)\n\n");

    double A=1.0;
    /* I(p) from charge_sym.mac PART 1 */
    double Imap[17]={0,4.0,4.591174,0,5.057061,0,0,0,5.405336,0,0,0,0,0,0,0,5.658564};

    /* ---- (1)(2) E_core vs p, open chain, one kink ---- */
    printf("== E_core vs p (A=1.0, N=400, one kink k=1) ==\n");
    printf("%-4s %-12s %-12s %-9s %-9s %-10s\n",
           "p","E_measured","sqrt(2A)I(p)","ratio","xi[sites]","max|grad|");
    for(int p=1;p<=16;p=(p==1)?2:p*2){
        Chain*c=chain_new(400,p,A,0,2*M_PI);
        double mg=0;
        double E=relax_best(c,1.0,400000,0.02,&mg);
        printf("%-4d %-12.6f %-12.6f %-9.4f %-9.0f %-10.1e\n",
               p, E, sqrt(2*A)*Imap[p], E/(sqrt(2*A)*Imap[p]), core_width(c), mg);
        chain_free(c);
    }
    printf("reading: the discrete minimum sits BELOW the continuum bound by the\n");
    printf("         lattice correction (2.8%% at p=8); the old 0.994 was the saddle.\n");

    /* ---- (1b) E_core vs A (p=8) ---- */
    printf("\n== E_core vs A (p=8, N=400, k=1) — expect sqrt(A) scaling ==\n");
    printf("%-8s %-12s %-14s %-9s %-10s\n","A","E_measured","sqrt(2A)*I(8)","ratio","max|grad|");
    double As[]={0.25,0.5,1.0,2.0,4.0};
    for(int k=0;k<5;k++){
        Chain*c=chain_new(400,8,As[k],0,2*M_PI);
        double mg=0;
        double E=relax_best(c,1.0,400000,0.02,&mg);
        printf("%-8.2f %-12.6f %-14.6f %-9.4f %-10.1e\n",
               As[k],E,sqrt(2*As[k])*5.405336,E/(sqrt(2*As[k])*5.405336),mg);
        chain_free(c);
    }
    printf("reading: the continuum bound is approached from below as A grows\n");
    printf("         (the core widens and the lattice correction shrinks).\n");

    /* ---- (4) SIZE CONDITION: periodic cycle, integer winding k, scan N ---- */
    printf("\n== SIZE CONDITION: periodic cycle, integer k, p=8, A=1.0 ==\n");
    printf("# NOTE on xi. The continuum core width sqrt(2A/V''(0)) = %.2f site at\n",
           sqrt(2*1.0/Vdd0(8)));
    printf("# the standing p=8, A=1 -- SUB-LATTICE. The kink is deep in the\n");
    printf("# discrete regime, so the continuum width is NOT the size scale, and\n");
    printf("# the '~4 sites' quoted elsewhere is the FWHM of the step profile (a\n");
    printf("# lattice quantity). They are not interchangeable. The size floor\n");
    printf("# below is measured, and it collapses onto N/k.\n");
    printf("%-4s %-6s %-7s %-12s %-12s %-12s %-10s\n",
           "k","N","N/k","E_total","E/(k*E1)","maxstep","max|grad|");
    double E1=0;
    {   Chain*c=chain_new(400,8,1.0,1,2*M_PI);
        double mg=0; E1=relax_best(c,1.0,400000,0.02,&mg); chain_free(c); }
    printf("# reference E1 (k=1, N=400 cycle) = %.6f\n",E1);
    int kk[]={1,1,1,1,2,2,2,3,3};
    int Ns[]={3,4,6,10,6,8,12,9,15};
    for(int j=0;j<9;j++){
        int k=kk[j],N=Ns[j];
        Chain*c=chain_new(N,8,1.0,1,2*M_PI*k);
        double mg=0;
        double E=relax_best(c,(double)k,400000,0.02,&mg);
        printf("%-4d %-6d %-7.2f %-12.6f %-12.4f %-12.4f %-10.1e\n",
               k,N,(double)N/k,E,E/(k*E1),max_step(c),mg);
        chain_free(c);
    }
    printf("reading: E/(k*E1) depends ONLY on N/k (the k=1,2,3 rows at equal N/k\n");
    printf("         agree to all printed digits) -- the cores are independent and\n");
    printf("         each needs its own room. Measured floor: N/k >= 6 costs <0.3%%,\n");
    printf("         N/k = 4 costs 6.6%%, N/k = 3 costs 15.4%%.\n");
    printf("         The SHARP criterion is the PN barrier (charge_fk2.c N2),\n");
    printf("         not this energy excess.\n");

    /* ---- (5) integer charge scaling: linear in |k| ---- */
    printf("\n== INTEGER charge scaling (p=8, A=1.0, N=400 cycle) ==\n");
    printf("%-6s %-12s %-12s %-9s\n","k","E_core","|k|*E1","ratio");
    for(int k=1;k<=3;k++){
        Chain*c=chain_new(400,8,1.0,1,2*M_PI*k);
        double mg=0;
        double E=relax_best(c,(double)k,400000,0.02,&mg);
        printf("%-6d %-12.6f %-12.6f %-9.4f\n",k,E,k*E1,E/(k*E1));
        chain_free(c);
    }
    printf("reading: integer charges are additive to <1%% — the cores are\n");
    printf("         contact-range (charge_sym.mac PART 4) so they do not\n");
    printf("         interact at N/k = 133 sites apart.\n");

    /* ---- (6) NEGATIVE CONTROL: fractional winding on a UNISON ring ---- */
    printf("\n== (6) NEGATIVE CONTROL: fractional twist on a unison ring ==\n");
    printf("2*pi/Q is NOT a period of V, so a 2*pi*k/Q twist on a unison ring\n");
    printf("is not a soliton. The diagnostic below shows what it actually is.\n");
    printf("%-8s %-12s %-12s %-12s %-14s\n",
           "w=k/Q","E_relax","seam bond","max interior","non-vacuum sites");
    double wf[]={1.0,1.0/3.0,2.0/3.0};
    for(int j=0;j<3;j++){
        double w=wf[j];
        Chain*c=chain_new(400,8,1.0,1,2*M_PI*w);
        double mg=0;
        double E=relax_best(c,w,400000,0.02,&mg);
        double seam=fabs(c->th[0]+c->twist-c->th[c->N-1]);
        double mxi=0;
        for(int i=1;i<c->N;i++){ double d=fabs(c->th[i]-c->th[i-1]); if(d>mxi)mxi=d; }
        printf("%-8.4f %-12.6f %-12.4f %-12.4f %-14d\n",w,E,seam,mxi,nonvacuum(c));
        chain_free(c);
    }
    printf("reading: at w=1 the twist is carried by a localised CORE (seam ~ 0,\n");
    printf("         a large interior step, a handful of non-vacuum sites).\n");
    printf("         At w=1/3 and 2/3 the twist sits on the SEAM BOND and the\n");
    printf("         chain is at the vacuum everywhere else: an elastic strain\n");
    printf("         imposed by the boundary condition, not a fractional charge.\n");
    printf("         => no conclusion about fractional charge can be drawn from\n");
    printf("         a unison ring. Test the 1/Q scaling on a Z_Q INTERVAL chain\n");
    printf("         whose gate variable is the comb phase q*th_i - p*th_j.\n");

    printf("\n# charge_fk DONE\n");
    return 0;
}
