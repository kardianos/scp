/*  fabric_measure.c — inverse programme for fabric_harmonic
 *
 *  Controlled measurements of the *actual* instrument dynamics:
 *    M1  instantaneous rates (RHS) at prescribed two-cell states
 *    M2  secular (time-averaged) transfer vs detuning
 *    M3  theta response to action (does density do work?)
 *    M4  two-lump interaction vs separation (attractive?)
 *    M5  free-spectrum tower: commensurate vs detuned (species control?)
 *    M6  single-cell multi-mode lock residual evolution
 *
 *  Same EOM as fabric_harmonic.c (duplicated deliberately so this file is
 *  self-contained and the measured model is frozen).
 *
 *  Build: gcc -O3 -march=native -fopenmp -o fabric_measure fabric_measure.c -lm
 *  Usage: fabric_measure [M1|M2|M3|M4|M5|M6|all] [outdir]
 *  Writes TSV tables under outdir (default v88/measure_out).
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>
#include <omp.h>

#define DMAX 4
#define MMAX 4

static int D = 1, L = 64, N = 0, M = 3;   /* default: 1D chain for two-cell clarity */
static int stride[DMAX];
static double WBAR[MMAX];
static double SIG = 0.40, THMAX = 0.80, GAM = 1.00, KTH = 1.00, MU = 0.50, ETA = 0.30;
static double EPS1 = 0.15, EPS2 = 0.05, GINT = 0.08, DT = 0.004;
static int TOWER_MODE = 0; /* 0=detuned (default instrument), 1=commensurate, 2=irrational */

static double *zr, *zi, *th, *Isum;

static inline int ZI(int i, int a) { return i * M + a; }
static inline int nbr(int idx, int a, int dir) {
    int c = (idx / stride[a]) % L;
    int nc = ((c + dir) % L + L) % L;
    return idx + (nc - c) * stride[a];
}

static void set_spectrum(void) {
    for (int a = 0; a < M; a++) {
        double al = (double)(a + 1);
        if (TOWER_MODE == 1)      WBAR[a] = al;                          /* commensurate */
        else if (TOWER_MODE == 2) WBAR[a] = al * (1.0 + 0.17 * sin(1.7*al)); /* strong detune */
        else                      WBAR[a] = al * (1.0 + 0.03 * al);      /* instrument default */
    }
}
static inline double w_alpha(int a, double theta) {
    return WBAR[a] * (1.0 - SIG * theta);
}
static inline double theta_eq_from_I(const double *Il) {
    double s = 0;
    for (int a = 0; a < M; a++) s += (a+1.0)*(a+1.0) * Il[a];
    return -THMAX * tanh(GAM * s);
}

static void alloc_fields(void) {
    N = 1; for (int a = 0; a < D; a++) N *= L;
    stride[0] = 1;
    for (int a = 1; a < D; a++) stride[a] = stride[a-1] * L;
    size_t nm = (size_t)N * M;
    zr = realloc(zr, sizeof(double)*nm);
    zi = realloc(zi, sizeof(double)*nm);
    th = realloc(th, sizeof(double)*N);
    Isum = realloc(Isum, sizeof(double)*N);
    memset(zr, 0, sizeof(double)*nm);
    memset(zi, 0, sizeof(double)*nm);
    memset(th, 0, sizeof(double)*N);
    set_spectrum();
}

/* ---- identical EOM core to fabric_harmonic (frozen) ---- */
static void deriv_z(const double *ar, const double *ai, const double *theta,
                    double *br, double *bi) {
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < N; i++) {
        for (int a = 0; a < M; a++) {
            int k = ZI(i, a);
            double wr = 0, wi = 0;
            double om = w_alpha(a, theta[i]);
            wr += om * ar[k]; wi += om * ai[k];

            if (EPS1 != 0.0) {
                double sr = 0, si = 0;
                for (int d = 0; d < D; d++) {
                    int p = nbr(i,d,+1), m = nbr(i,d,-1);
                    int kp = ZI(p,a), km = ZI(m,a);
                    sr += ar[kp]+ar[km]-2.0*ar[k];
                    si += ai[kp]+ai[km]-2.0*ai[k];
                }
                wr += EPS1*sr; wi += EPS1*si;
            }
            if (EPS2 != 0.0 && M >= 2) {
                int a1 = a + 1;
                if (2*a1 <= M) {
                    int b = 2*a1 - 1;
                    for (int d = 0; d < D; d++) {
                        int p = nbr(i,d,+1), m = nbr(i,d,-1);
                        for (int t = 0; t < 2; t++) {
                            int j = (t==0)?p:m;
                            int kb = ZI(j,b);
                            double car = ar[k], cai = -ai[k];
                            double zbr = ar[kb], zbi = ai[kb];
                            wr += 2.0*EPS2*(car*zbr - cai*zbi);
                            wi += 2.0*EPS2*(car*zbi + cai*zbr);
                        }
                    }
                }
                if (a1 % 2 == 0) {
                    int lo = a1/2 - 1;
                    for (int d = 0; d < D; d++) {
                        int p = nbr(i,d,+1), m = nbr(i,d,-1);
                        for (int t = 0; t < 2; t++) {
                            int j = (t==0)?p:m;
                            int kl = ZI(j,lo);
                            wr += EPS2*(ar[kl]*ar[kl]-ai[kl]*ai[kl]);
                            wi += EPS2*(2.0*ar[kl]*ai[kl]);
                        }
                    }
                }
            }
            if (GINT != 0.0 && M >= 3) {
                int k0 = ZI(i,0), k1 = ZI(i,1), k2 = ZI(i,2);
                if (a==0) {
                    double c2r=ar[k1], c2i=-ai[k1], t3r=ar[k2], t3i=ai[k2];
                    wr += GINT*(c2r*t3r - c2i*t3i);
                    wi += GINT*(c2r*t3i + c2i*t3r);
                } else if (a==1) {
                    double c1r=ar[k0], c1i=-ai[k0], t3r=ar[k2], t3i=ai[k2];
                    wr += GINT*(c1r*t3r - c1i*t3i);
                    wi += GINT*(c1r*t3i + c1i*t3r);
                } else if (a==2) {
                    wr += GINT*(ar[k0]*ar[k1]-ai[k0]*ai[k1]);
                    wi += GINT*(ar[k0]*ai[k1]+ai[k0]*ar[k1]);
                }
            }
            br[k] =  wi;
            bi[k] = -wr;
        }
    }
}

static double *nth_buf, *kr_buf, *ki_buf, *tr_buf, *ti_buf;
static size_t buf_n = 0, buf_nm = 0;

static void ensure_bufs(void) {
    size_t nm = (size_t)N * M;
    if (N > (int)buf_n) {
        nth_buf = realloc(nth_buf, sizeof(double) * N);
        buf_n = (size_t)N;
    }
    if (nm > buf_nm) {
        kr_buf = realloc(kr_buf, sizeof(double) * nm);
        ki_buf = realloc(ki_buf, sizeof(double) * nm);
        tr_buf = realloc(tr_buf, sizeof(double) * nm);
        ti_buf = realloc(ti_buf, sizeof(double) * nm);
        buf_nm = nm;
    }
}

static void step_theta(void) {
    ensure_bufs();
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < N; i++) {
        double Il[MMAX];
        for (int a = 0; a < M; a++) {
            int k = ZI(i,a);
            Il[a] = zr[k]*zr[k]+zi[k]*zi[k];
        }
        double teq = theta_eq_from_I(Il);
        double lap = 0;
        for (int d = 0; d < D; d++) {
            int p = nbr(i,d,+1), m = nbr(i,d,-1);
            lap += th[p]+th[m]-2.0*th[i];
        }
        double force = -KTH*(th[i]-teq) + MU*lap;
        nth_buf[i] = th[i] + DT*ETA*force;
        if (nth_buf[i] > 1.5) nth_buf[i]=1.5;
        if (nth_buf[i] <-1.5) nth_buf[i]=-1.5;
    }
    memcpy(th, nth_buf, sizeof(double)*N);
}

static void step_z(void) {
    ensure_bufs();
    size_t nm = (size_t)N*M;
    deriv_z(zr,zi,th,kr_buf,ki_buf);
    for (size_t q=0;q<nm;q++){ tr_buf[q]=zr[q]+0.5*DT*kr_buf[q]; ti_buf[q]=zi[q]+0.5*DT*ki_buf[q]; }
    deriv_z(tr_buf,ti_buf,th,kr_buf,ki_buf);
    for (size_t q=0;q<nm;q++){ zr[q]+=DT*kr_buf[q]; zi[q]+=DT*ki_buf[q]; }
}

static void step(void) { step_z(); step_theta(); }

static void zero_all(void) {
    size_t nm = (size_t)N*M;
    memset(zr,0,sizeof(double)*nm); memset(zi,0,sizeof(double)*nm);
    memset(th,0,sizeof(double)*N);
}

static void set_cell_mode(int i, int a, double amp, double phi) {
    int k = ZI(i,a);
    zr[k] = amp * cos(phi);
    zi[k] = amp * sin(phi);
}

static void get_I_phi(int i, int a, double *I, double *phi) {
    int k = ZI(i,a);
    *I = zr[k]*zr[k]+zi[k]*zi[k];
    *phi = atan2(zi[k], zr[k]);
}

/* dI/dt = 2(re*bre + im*bim), dphi/dt from Im/Re of ż/z */
static void rates_at_cell(int i, int a, const double *br, const double *bi,
                          double *dI, double *dphi) {
    int k = ZI(i,a);
    double re=zr[k], im=zi[k], I=re*re+im*im;
    *dI = 2.0*(re*br[k] + im*bi[k]);
    if (I < 1e-30) { *dphi = 0; return; }
    /* dphi/dt = (re bim - im bre)/I */
    *dphi = (re*bi[k] - im*br[k]) / I;
}

/* ================================================================
 * M1: instantaneous two-cell 1:1 transfer surface
 * cells 0 and 1 on a 1D ring; only mode 0 excited; vary I, dphi, dtheta
 * ================================================================ */
static void measure_M1(const char *outdir) {
    char path[512]; snprintf(path,sizeof path,"%s/M1_transfer.tsv",outdir);
    FILE *f = fopen(path,"w");
    if (!f) { perror(path); return; }
    fprintf(f, "# M1 instantaneous 1:1 transfer (mode 0)\n");
    fprintf(f, "# eps1=%.4f eps2=%.4f gint=%.4f M=%d D=%d L=%d\n", EPS1,EPS2,GINT,M,D,L);
    fprintf(f, "I0\tI1\tdphi\tth0\tth1\tdetune\tdI0\tdI1\tdphi0\tdphi1\tP_theory\n");

    D=1; L=32; M=3; alloc_fields();
    size_t nm=(size_t)N*M;
    double *br=malloc(sizeof(double)*nm), *bi=malloc(sizeof(double)*nm);

    double Is[] = {0.05, 0.1, 0.2, 0.4, 0.8};
    int nI = 5;
    for (int i0=0;i0<nI;i0++) for (int i1=0;i1<nI;i1++) {
        for (int ip=0; ip<16; ip++) {
            double dphi = ip * (2.0*M_PI/16.0);
            for (int it=0; it<5; it++) {
                double th0 = -0.2*it; /* 0, -0.2, -0.4, -0.6, -0.8 */
                double th1 = 0.0;
                zero_all();
                set_cell_mode(0, 0, sqrt(Is[i0]), 0.0);
                set_cell_mode(1, 0, sqrt(Is[i1]), dphi);
                th[0]=th0; th[1]=th1;
                /* isolate: kill all other couplings by zeroing other cells (already zero) */
                deriv_z(zr,zi,th,br,bi);
                double dI0,dI1,dp0,dp1;
                rates_at_cell(0,0,br,bi,&dI0,&dp0);
                rates_at_cell(1,0,br,bi,&dI1,&dp1);
                double om0=w_alpha(0,th0), om1=w_alpha(0,th1);
                double det = om0-om1;
                /* pure 1:1 theory for H=eps1*(z0* z1 + c.c.) without Laplacian onsite:
                   with Laplacian, hopping to empty nbrs adds -2D eps1 term on dphi, not dI for isolated pair:
                   for two sites on ring of L, each has one partner "on" and one "off"
                   dI0 = 2 eps1 sqrt(I0 I1) sin(phi0-phi1) * (connectivity factor)
                */
                double Pth = 2.0*EPS1*sqrt(Is[i0]*Is[i1])*sin(-dphi); /* expect form */
                fprintf(f, "%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\n",
                        Is[i0],Is[i1],dphi,th0,th1,det,dI0,dI1,dp0,dp1,Pth);
            }
        }
    }
    fclose(f);
    free(br); free(bi);
    printf("M1 wrote %s\n", path);
}

/* ================================================================
 * M2: secular transfer vs detuning (short time average)
 * Two cells, fixed amplitudes projected each step? No: free evolve short T
 * with relative phase scan; average dI over common period estimate.
 * Better: fix amplitudes by re-normalizing each step (gauge secular drive)
 * Actually: use instantaneous dI at many phases and average analytically
 * via sampling — secular coefficient = <dI * sin(dphi)> / ...
 * Also scan detuning by th0.
 * ================================================================ */
static void measure_M2(const char *outdir) {
    char path[512]; snprintf(path,sizeof path,"%s/M2_secular.tsv",outdir);
    FILE *f = fopen(path,"w");
    fprintf(f, "# M2 secular 1:1 transfer amplitude vs detuning\n");
    fprintf(f, "# A_sec := max_dphi |dI0| at fixed I; also mean dI0 over dphi (should ~0)\n");
    fprintf(f, "I\tdetune\tth0\tom0\tom1\tA_sec\tmean_dI\trms_dI\tA_theory\tphase_of_max\n");

    D=1; L=32; M=3; alloc_fields();
    size_t nm=(size_t)N*M;
    double *br=malloc(sizeof(double)*nm), *bi=malloc(sizeof(double)*nm);

    double I = 0.25;
    for (int it=0; it<=20; it++) {
        double th0 = -0.04 * it; /* 0 .. -0.8 */
        double th1 = 0.0;
        double om0=w_alpha(0,th0), om1=w_alpha(0,th1);
        double det = om0 - om1;
        double Asec=0, mean=0, m2=0, phmax=0;
        int nph=64;
        for (int ip=0; ip<nph; ip++) {
            double dphi = ip*(2.0*M_PI/nph);
            zero_all();
            set_cell_mode(0,0,sqrt(I),0.0);
            set_cell_mode(1,0,sqrt(I),dphi);
            th[0]=th0; th[1]=th1;
            deriv_z(zr,zi,th,br,bi);
            double dI0,dp0; rates_at_cell(0,0,br,bi,&dI0,&dp0);
            mean += dI0; m2 += dI0*dI0;
            if (fabs(dI0) > Asec) { Asec=fabs(dI0); phmax=dphi; }
        }
        mean /= nph; m2 = sqrt(m2/nph);
        /* theory for resonant 1:1 Laplacian pair: |dI|_max = 2 eps1 I  (same I)
           OFF resonance the INSTANTANEOUS amp is unchanged (detune enters dphi not dI
           for bilinear hopping!). Secular average of transfer over free evolution
           is what dies with detuning — measure that next. */
        double Ath = 2.0*EPS1*I;
        fprintf(f, "%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.8e\t%.8e\t%.8e\t%.8e\t%.6f\n",
                I, det, th0, om0, om1, Asec, mean, m2, Ath, phmax);
    }
    fclose(f);
    free(br); free(bi);

    /* M2b: true secular — free evolve STRICT two-site ring (L=2), net Delta I */
    snprintf(path,sizeof path,"%s/M2b_secular_evolve.tsv",outdir);
    f = fopen(path,"w");
    fprintf(f, "# M2b free evolve two-site ring: net Delta I0 over T_avg vs detuning\n");
    fprintf(f, "I\tdetune\tth0_fixed\tT\tDeltaI0\tDeltaI1\tI0_final\tI1_final\n");
    double Tavg = 40.0;
    int nsteps = (int)(Tavg/DT);
    D=1; L=2; M=3; alloc_fields(); ensure_bufs();
    for (int it=0; it<=20; it++) {
        double th0 = -0.04*it, th1=0.0;
        double om0=w_alpha(0,th0), om1=w_alpha(0,th1);
        double det = om0-om1;
        zero_all();
        set_cell_mode(0,0,sqrt(I),0.0);
        set_cell_mode(1,0,sqrt(I),0.5);
        th[0]=th0; th[1]=th1;
        double I0i,I1i,ph; get_I_phi(0,0,&I0i,&ph); get_I_phi(1,0,&I1i,&ph);
        for (int s=0;s<nsteps;s++) {
            step_z();
            th[0]=th0; th[1]=th1;
        }
        double I0f,I1f; get_I_phi(0,0,&I0f,&ph); get_I_phi(1,0,&I1f,&ph);
        fprintf(f, "%.6f\t%.6f\t%.6f\t%.4f\t%.8e\t%.8e\t%.8e\t%.8e\n",
                I, det, th0, Tavg, I0f-I0i, I1f-I1i, I0f, I1f);
    }
    fclose(f);
    printf("M2 wrote secular tables\n");
    fflush(stdout);
}

/* ================================================================
 * M3: theta response — is density decorative?
 * ================================================================ */
static void measure_M3(const char *outdir) {
    char path[512]; snprintf(path,sizeof path,"%s/M3_theta.tsv",outdir);
    FILE *f = fopen(path,"w");
    fprintf(f, "# M3 theta response: single cell forced I, relax theta, measure omega shift\n");
    fprintf(f, "I_tot\tI0\tI1\tI2\ttheta_eq_analytic\ttheta_relaxed\tdom0\tdom_rel\n");

    D=1; L=16; M=3; alloc_fields(); ensure_bufs();
    for (int k=0; k<=30; k++) {
        double Itot = 0.02 * k; /* 0 .. 0.6 */
        /* put all in mode 0 */
        zero_all();
        set_cell_mode(L/2, 0, sqrt(Itot), 0.0);
        double Il[MMAX]={0}; Il[0]=Itot;
        double teq = theta_eq_from_I(Il);
        /* relax theta only, hold z */
        th[L/2]=0;
        for (int s=0;s<2000;s++) {
            step_theta();
        }
        double thr = th[L/2];
        double om0 = w_alpha(0, 0.0);
        double omr = w_alpha(0, thr);
        fprintf(f, "%.6f\t%.6f\t0\t0\t%.8e\t%.8e\t%.8e\t%.8e\n",
                Itot, Itot, teq, thr, omr-om0, (omr-om0)/ (om0>0?om0:1));
    }
    /* also: required I to reach theta = -0.5 */
    fclose(f);

    snprintf(path,sizeof path,"%s/M3_scale.tsv",outdir);
    f = fopen(path,"w");
    fprintf(f, "# What GAM, THMAX give |theta|=0.5 at I=0.25 (mode0 only)? scan\n");
    fprintf(f, "GAM\tTHMAX\ttheta_eq_at_I0.25\n");
    double gams[]={0.5,1,2,4,8,16,32,64};
    double thms[]={0.4,0.8,1.2,1.6,2.0};
    for (int ig=0;ig<8;ig++) for (int it=0;it<5;it++) {
        double g=gams[ig], tm=thms[it];
        double s = 1.0*1.0*0.25; /* w_alpha weight alpha^2 * I */
        double teq = -tm * tanh(g * s);
        fprintf(f, "%.4f\t%.4f\t%.8e\n", g, tm, teq);
    }
    fclose(f);
    printf("M3 wrote theta tables\n");
}

/* ================================================================
 * M4: two-lump interaction vs separation
 * Seed two Gaussian-ish single-mode lumps on 1D chain, measure
 *   - total energy vs separation (static, after short relax)
 *   - d(sep)/dt from centroid motion over time
 * ================================================================ */
static double total_I(void) {
    double s=0;
    for (int i=0;i<N;i++) for (int a=0;a<M;a++) {
        int k=ZI(i,a); s+=zr[k]*zr[k]+zi[k]*zi[k];
    }
    return s;
}
static double free_energy(void) {
    /* E ≈ Σ ω(θ)I + ½K(θ-θeq)² + (μ/2)Σ(θi-θj)² + eps1 hop (real part) */
    double e=0;
    for (int i=0;i<N;i++) {
        double Il[MMAX];
        for (int a=0;a<M;a++) {
            int k=ZI(i,a);
            Il[a]=zr[k]*zr[k]+zi[k]*zi[k];
            e += w_alpha(a,th[i])*Il[a];
        }
        double teq=theta_eq_from_I(Il);
        double d=th[i]-teq; e += 0.5*KTH*d*d;
        for (int d2=0;d2<D;d2++) {
            int p=nbr(i,d2,+1);
            double dt=th[i]-th[p];
            e += 0.25*MU*dt*dt;
            for (int a=0;a<M;a++) {
                int k=ZI(i,a), kp=ZI(p,a);
                /* hop energy eps1 Re(z_i* z_p) with Laplacian form:
                   actually H_hop = (eps1/2) Σ |z_i-z_j|^2 style */
                double dr=zr[k]-zr[kp], di=zi[k]-zi[kp];
                e += 0.5*EPS1*(dr*dr+di*di);
            }
        }
    }
    return e;
}

static void place_lump(int c, double amp, double width, double phi) {
    for (int i=0;i<N;i++) {
        int dx=i-c; if (dx>L/2) dx-=L; if (dx<-L/2) dx+=L;
        double env = amp * exp(-0.5*(dx*dx)/(width*width));
        if (env < 1e-6) continue;
        zr[ZI(i,0)] += env * cos(phi);
        zi[ZI(i,0)] += env * sin(phi);
        th[i] += -0.15 * (env/amp); /* mild densification hint */
    }
}

static double centroid_mode0(void) {
    double m=0, w=0;
    for (int i=0;i<N;i++) {
        double I=zr[ZI(i,0)]*zr[ZI(i,0)]+zi[ZI(i,0)]*zi[ZI(i,0)];
        /* unwrap around max */
        m += I * i; w += I;
    }
    return (w>0)? m/w : 0;
}

/* two centroids via half-split — crude: find two peaks */
static void two_centroids(double *c1, double *c2, double *I1, double *I2) {
    int p1=-1; double pk=0;
    for (int i=0;i<N;i++) {
        double I=zr[ZI(i,0)]*zr[ZI(i,0)]+zi[ZI(i,0)]*zi[ZI(i,0)];
        if (I>pk){pk=I;p1=i;}
    }
    int p2=-1; pk=0;
    for (int i=0;i<N;i++) {
        int d=i-p1; if (d>L/2)d-=L; if(d<-L/2)d+=L;
        if (abs(d)<3) continue;
        double I=zr[ZI(i,0)]*zr[ZI(i,0)]+zi[ZI(i,0)]*zi[ZI(i,0)];
        if (I>pk){pk=I;p2=i;}
    }
    /* COM in neighbourhood of each peak */
    double m1=0,w1=0,m2=0,w2=0;
    for (int i=0;i<N;i++) {
        double I=zr[ZI(i,0)]*zr[ZI(i,0)]+zi[ZI(i,0)]*zi[ZI(i,0)];
        int d1=i-p1; if(d1>L/2)d1-=L; if(d1<-L/2)d1+=L;
        int d2=i-p2; if(d2>L/2)d2-=L; if(d2<-L/2)d2+=L;
        if (abs(d1)<=5){ m1+=I*(p1+d1); w1+=I; }
        if (p2>=0 && abs(d2)<=5){ m2+=I*(p2+d2); w2+=I; }
    }
    *c1 = (w1>0)?m1/w1:p1; *I1=w1;
    *c2 = (w2>0)?m2/w2:p2; *I2=w2;
}

static void measure_M4(const char *outdir) {
    char path[512]; snprintf(path,sizeof path,"%s/M4_interaction.tsv",outdir);
    FILE *f = fopen(path,"w");
    fprintf(f, "# M4 two-lump static energy vs separation (1D)\n");
    fprintf(f, "sep\tE\tE_minus_2E1\tI_tot\ttheta_min\n");

    D=1; L=64; M=3; alloc_fields(); ensure_bufs();
    /* single lump energy reference — short relax */
    zero_all();
    place_lump(L/2, 0.7, 2.0, 0.0);
    for (int s=0;s<800;s++) step();
    double E1 = free_energy();
    double Iref = total_I();

    for (int sep=2; sep<=24; sep++) {
        zero_all();
        int c1 = L/2 - sep/2;
        int c2 = L/2 + (sep - sep/2);
        place_lump(c1, 0.7, 2.0, 0.0);
        place_lump(c2, 0.7, 2.0, 0.3);
        for (int s=0;s<800;s++) step();
        double E = free_energy();
        double It = total_I();
        double tmin=0; for(int i=0;i<N;i++) if(th[i]<tmin)tmin=th[i];
        fprintf(f, "%d\t%.8e\t%.8e\t%.8e\t%.6f\n", sep, E, E-2.0*E1, It, tmin);
    }
    fclose(f);

    /* M4b: dynamical approach/separation rate */
    snprintf(path,sizeof path,"%s/M4b_force.tsv",outdir);
    f = fopen(path,"w");
    fprintf(f, "# M4b d(sep)/dt after short evolve from initial sep\n");
    fprintf(f, "sep0\tsep_T\tdsep_dt\tI1\tI2\tE\n");
    double T=20.0; int nst=(int)(T/DT);
    for (int sep=3; sep<=20; sep+=1) {
        zero_all();
        int c1 = L/2 - sep/2;
        int c2 = L/2 + (sep - sep/2);
        place_lump(c1, 0.7, 2.0, 0.0);
        place_lump(c2, 0.7, 2.0, 0.0); /* same phase */
        double ca,cb,ia,ib; two_centroids(&ca,&cb,&ia,&ib);
        double sep0 = fabs(cb-ca); if (sep0>L/2) sep0=L-sep0;
        for (int s=0;s<nst;s++) step();
        two_centroids(&ca,&cb,&ia,&ib);
        double sep1 = fabs(cb-ca); if (sep1>L/2) sep1=L-sep1;
        fprintf(f, "%d\t%.6f\t%.8e\t%.6f\t%.6f\t%.8e\n",
                sep, sep1, (sep1-sep0)/T, ia, ib, free_energy());
    }
    fclose(f);

    /* M4c: same with opposite phase */
    snprintf(path,sizeof path,"%s/M4c_force_antiphase.tsv",outdir);
    f = fopen(path,"w");
    fprintf(f, "# M4c d(sep)/dt antiphase (phi=pi)\n");
    fprintf(f, "sep0\tsep_T\tdsep_dt\tI1\tI2\tE\n");
    for (int sep=3; sep<=20; sep+=1) {
        zero_all();
        int c1 = L/2 - sep/2;
        int c2 = L/2 + (sep - sep/2);
        place_lump(c1, 0.7, 2.0, 0.0);
        place_lump(c2, 0.7, 2.0, M_PI);
        double ca,cb,ia,ib; two_centroids(&ca,&cb,&ia,&ib);
        double sep0 = fabs(cb-ca); if (sep0>L/2) sep0=L-sep0;
        for (int s=0;s<nst;s++) step();
        two_centroids(&ca,&cb,&ia,&ib);
        double sep1 = fabs(cb-ca); if (sep1>L/2) sep1=L-sep1;
        fprintf(f, "%d\t%.6f\t%.8e\t%.6f\t%.6f\t%.8e\n",
                sep, sep1, (sep1-sep0)/T, ia, ib, free_energy());
    }
    fclose(f);
    printf("M4 wrote interaction tables (E1=%.6f Iref=%.6f)\n", E1, Iref);
    fflush(stdout);
}

/* ================================================================
 * M5: tower mode — commensurate vs detuned, lock residual + transfer
 * ================================================================ */
static void measure_M5(const char *outdir) {
    char path[512]; snprintf(path,sizeof path,"%s/M5_tower.tsv",outdir);
    FILE *f = fopen(path,"w");
    fprintf(f, "# M5 tower: lock residual |n·ω| and 3-wave dI rates for tower modes\n");
    fprintf(f, "tower\tW0\tW1\tW2\tlock_210\tlock_12m1\tdI0_3wave\tdI1_3wave\tdI2_3wave\n");

    D=1; L=16; M=3;
    for (int tm=0; tm<3; tm++) {
        TOWER_MODE = tm;
        alloc_fields();
        /* equal amplitude multi-mode cell */
        zero_all();
        double amp=0.4;
        set_cell_mode(0,0,amp,0.0);
        set_cell_mode(0,1,amp,0.2);
        set_cell_mode(0,2,amp,0.4);
        th[0]=0;
        double om[3];
        for (int a=0;a<3;a++) om[a]=w_alpha(a,0);
        double lock210 = fabs(2*om[0]-om[1]);
        double lock12m1 = fabs(-om[0]+2*om[1]-om[2]);
        size_t nm=(size_t)N*M;
        double *br=malloc(sizeof(double)*nm), *bi=malloc(sizeof(double)*nm);
        /* only GINT: kill hop for pure intra-cell */
        double e1=EPS1,e2=EPS2; EPS1=0; EPS2=0;
        deriv_z(zr,zi,th,br,bi);
        EPS1=e1; EPS2=e2;
        double dI[3], dp;
        for (int a=0;a<3;a++) rates_at_cell(0,a,br,bi,&dI[a],&dp);
        fprintf(f, "%d\t%.6f\t%.6f\t%.6f\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\n",
                tm, om[0],om[1],om[2], lock210, lock12m1, dI[0],dI[1],dI[2]);
        free(br); free(bi);
    }
    TOWER_MODE=0;
    fclose(f);
    printf("M5 wrote tower table\n");
}

/* ================================================================
 * M6: single-cell lock evolution + what Nc the field produces
 * from a multi-cell seed of width W, after evolution, measure Nc
 * ================================================================ */
static void measure_M6(const char *outdir) {
    char path[512]; snprintf(path,sizeof path,"%s/M6_Nc_vs_width.tsv",outdir);
    FILE *f = fopen(path,"w");
    fprintf(f, "# M6: seed Gaussian width W, evolve T, measure IPR and peak-cluster Nc\n");
    fprintf(f, "W\tamp\tT\tIPR_cells\tpeakI\tn_above_frac\tmean_theta_core\n");

    D=1; L=64; M=3; alloc_fields(); ensure_bufs();
    double T=20.0; int nst=(int)(T/DT);
    for (int iw=1; iw<=12; iw++) {
        double W = 0.5 * iw; /* 0.5 .. 6 */
        for (int ia=0; ia<3; ia++) {
            double amp = 0.4 + 0.3*ia;
            zero_all();
            place_lump(L/2, amp, W, 0.0);
            /* also seed weak mode 2 for multi-mode */
            for (int i=0;i<N;i++) {
                int dx=i-L/2; if(dx>L/2)dx-=L; if(dx<-L/2)dx+=L;
                double env = 0.3*amp*exp(-0.5*(dx*dx)/(W*W));
                zr[ZI(i,2)] += env;
            }
            for (int s=0;s<nst;s++) step();
            double s2=0,s4=0,pk=0,thc=0; int ncore=0;
            for (int i=0;i<N;i++) {
                double I=0; for(int a=0;a<M;a++){int k=ZI(i,a);I+=zr[k]*zr[k]+zi[k]*zi[k];}
                s2+=I; s4+=I*I; if(I>pk)pk=I;
            }
            double ipr = (s4>0)? (s2*s2)/s4 : 0; /* cells participating */
            double thr=0.22*pk; int nabove=0;
            for (int i=0;i<N;i++) {
                double I=0; for(int a=0;a<M;a++){int k=ZI(i,a);I+=zr[k]*zr[k]+zi[k]*zi[k];}
                if (I>=thr){ nabove++; thc+=th[i]; ncore++; }
            }
            fprintf(f, "%.4f\t%.4f\t%.2f\t%.6f\t%.6f\t%d\t%.6f\n",
                    W, amp, T, ipr, pk, nabove, ncore?thc/ncore:0);
        }
    }
    fclose(f);

    /* M6b: lock residual over time for single multi-mode cell (no hop) */
    snprintf(path,sizeof path,"%s/M6b_lock_evolve.tsv",outdir);
    f = fopen(path,"w");
    fprintf(f, "# M6b isolated cell (EPS1=EPS2=0) multi-mode phase residual\n");
    fprintf(f, "t\tI0\tI1\tI2\tpsi\tdpsi_dt_est\ttheta\n");
    D=1; L=8; M=3; alloc_fields(); ensure_bufs();
    double e1=EPS1,e2=EPS2; EPS1=0; EPS2=0;
    zero_all();
    set_cell_mode(0,0,0.5,0.0);
    set_cell_mode(0,1,0.4,0.5);
    set_cell_mode(0,2,0.3,1.0);
    double T2=30; int ns2=(int)(T2/DT);
    double prev_psi=0; int has_prev=0;
    for (int s=0;s<=ns2;s++) {
        if (s%((int)(0.5/DT))==0) {
            double I0,I1,I2,p0,p1,p2;
            get_I_phi(0,0,&I0,&p0); get_I_phi(0,1,&I1,&p1); get_I_phi(0,2,&I2,&p2);
            /* lock phase for 1+2-3: psi = p0+p1-p2 */
            double psi = p0+p1-p2;
            /* wrap */
            while (psi>M_PI) psi-=2*M_PI;
            while (psi<-M_PI) psi+=2*M_PI;
            double dpsi=0;
            if (has_prev) {
                dpsi = psi - prev_psi;
                while (dpsi>M_PI) dpsi-=2*M_PI;
                while (dpsi<-M_PI) dpsi+=2*M_PI;
                dpsi /= 0.5;
            }
            fprintf(f, "%.4f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n",
                    s*DT, I0,I1,I2, psi, dpsi, th[0]);
            prev_psi=psi; has_prev=1;
        }
        step();
    }
    EPS1=e1; EPS2=e2;
    fclose(f);
    printf("M6 wrote Nc and lock tables\n");
}

/* M7: effective two-cell potential from measured rates
 * V_eff(sep) is already M4; also measure action-exchange kernel K(delta)
 * from M2b.
 */
static void measure_M7_summary(const char *outdir) {
    char path[512]; snprintf(path,sizeof path,"%s/M7_params.tsv",outdir);
    FILE *f = fopen(path,"w");
    fprintf(f, "# instrument parameters frozen for this measurement campaign\n");
    fprintf(f, "name\tvalue\n");
    fprintf(f, "SIG\t%.6f\n", SIG);
    fprintf(f, "THMAX\t%.6f\n", THMAX);
    fprintf(f, "GAM\t%.6f\n", GAM);
    fprintf(f, "KTH\t%.6f\n", KTH);
    fprintf(f, "MU\t%.6f\n", MU);
    fprintf(f, "ETA\t%.6f\n", ETA);
    fprintf(f, "EPS1\t%.6f\n", EPS1);
    fprintf(f, "EPS2\t%.6f\n", EPS2);
    fprintf(f, "GINT\t%.6f\n", GINT);
    fprintf(f, "DT\t%.6f\n", DT);
    for (int a=0;a<M;a++) fprintf(f, "WBAR_%d\t%.6f\n", a, WBAR[a]);
    fclose(f);
}

int main(int argc, char **argv) {
    setvbuf(stdout, NULL, _IONBF, 0);
    const char *which = (argc>1)?argv[1]:"all";
    const char *outdir = (argc>2)?argv[2]:"measure_out";
    mkdir(outdir, 0755);
    set_spectrum();
    printf("fabric_measure: which=%s outdir=%s\n", which, outdir);

    if (!strcmp(which,"all") || !strcmp(which,"M1")) measure_M1(outdir);
    if (!strcmp(which,"all") || !strcmp(which,"M2")) measure_M2(outdir);
    if (!strcmp(which,"all") || !strcmp(which,"M3")) measure_M3(outdir);
    if (!strcmp(which,"all") || !strcmp(which,"M4")) measure_M4(outdir);
    if (!strcmp(which,"all") || !strcmp(which,"M5")) measure_M5(outdir);
    if (!strcmp(which,"all") || !strcmp(which,"M6")) measure_M6(outdir);
    measure_M7_summary(outdir);
    printf("done.\n");
    return 0;
}
