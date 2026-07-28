/*  fabric_harmonic.c — v88 multi-dimensional harmonic fabric
 *
 *  NOT a discrete-breather / DNLS model. fabric_trap.c is the cautionary
 *  example: one complex mode, amplitude-shifted frequency, continuous dial.
 *
 *  HERE each cell carries m>=2 internal harmonics z_a plus densification
 *  theta. Net energy transfer between neighbours is resonant (complete-cycle
 *  secular average): 1:1 and 2:1 channels, plus an intra-cell 1+2<->3 lock
 *  that is the emergent interior configuration ("Higgs" = pattern, not a
 *  field). Integers enter as mode indices, resonance orders (p:q), and lock
 *  vectors -- not as amplitude thresholds.
 *
 *  Theory: v88/GROK_V88_DESIGN.md
 *  Ontology: v88/ONTOLOGY.md
 *
 *  Build: gcc -O3 -march=native -fopenmp -o fabric_harmonic fabric_harmonic.c -lm
 *  Usage: fabric_harmonic [schedule] [d] [seed]
 *    schedule: smoke | quench | thermal | cyclic | driven | staged |
 *              coalesce | detune
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>

#define DMAX 4
#define MMAX 4
#define MAX_CLU 8192

/* ---- parameters (defaults from GROK_V88_DESIGN.md) ---- */
static int D = 3, L = 20, N = 0, M = 3;
static int stride[DMAX];

static double WBAR[MMAX];       /* free frequencies at theta=0 */
static double SIG = 0.40;       /* omega shifts with theta */
static double THMAX = 0.80;     /* hard-core densification scale */
static double GAM = 1.00;       /* cyclic energy -> tightness */
static double KTH = 1.00;       /* pull of theta toward theta_eq */
static double MU = 0.50;        /* elastic stiffness */
static double ETA = 0.30;       /* overdamped rate for theta */
static double EPS1 = 0.15;      /* 1:1 neighbour coupling */
static double EPS2 = 0.05;      /* 2:1 neighbour coupling */
static double GINT = 0.08;      /* intra-cell 1+2 <-> 3 */
static double DT = 0.008;
static double I_THR_FRAC = 0.22; /* cluster if I_i >= frac * peak I */

static double *zr, *zi;         /* N*M complex amplitudes (split) */
static double *th;              /* N dilations */
static double *Isum;            /* N total action cache */

static unsigned long long rs;

static double urand(void) {
    rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17;
    return (double)(rs >> 11) * (1.0 / 9007199254740992.0);
}
static double grand(void) { /* Box-Muller */
    double u = urand() + 1e-300, v = urand();
    return sqrt(-2.0 * log(u)) * cos(2.0 * M_PI * v);
}

static inline int nbr(int idx, int a, int dir) {
    int c = (idx / stride[a]) % L;
    int nc = ((c + dir) % L + L) % L;
    return idx + (nc - c) * stride[a];
}
static inline int ZI(int i, int a) { return i * M + a; }

static void set_spectrum(void) {
    for (int a = 0; a < M; a++) {
        double al = (double)(a + 1);
        /* weakly detuned harmonic tower: resonance is special, not automatic */
        WBAR[a] = al * (1.0 + 0.03 * al);
    }
}

static inline double w_alpha(int a, double theta) {
    return WBAR[a] * (1.0 - SIG * theta);
}

static inline double theta_eq_from_I(const double *Iloc) {
    double s = 0;
    for (int a = 0; a < M; a++) {
        double al = (double)(a + 1);
        s += al * al * Iloc[a];
    }
    return -THMAX * tanh(GAM * s);
}

static void cache_I(void) {
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < N; i++) {
        double s = 0;
        for (int a = 0; a < M; a++) {
            int k = ZI(i, a);
            s += zr[k] * zr[k] + zi[k] * zi[k];
        }
        Isum[i] = s;
    }
}

/* Right-hand side for z: i dz/dt = dH/dz*
 * Written as dz/dt = -i * ( ... )  =>  dr/dt, di/dt form. */
static void deriv_z(const double *ar, const double *ai, const double *theta,
                    double *br, double *bi) {
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < N; i++) {
        double Il[MMAX];
        for (int a = 0; a < M; a++) {
            int k = ZI(i, a);
            Il[a] = ar[k] * ar[k] + ai[k] * ai[k];
        }
        for (int a = 0; a < M; a++) {
            int k = ZI(i, a);
            double wr = 0, wi = 0; /* complex field multiplying z contribution */
            /* free rotation: omega * z  =>  i dz = omega z  => dz = -i omega z */
            double om = w_alpha(a, theta[i]);
            wr += om * ar[k];
            wi += om * ai[k];

            /* 1:1 neighbour Laplacian: eps1 * sum_j (z_j - z_i)
             * preserves total action; band is [0, 4 d eps1] */
            if (EPS1 != 0.0) {
                double sr = 0, si = 0;
                for (int d = 0; d < D; d++) {
                    int p = nbr(i, d, +1), m = nbr(i, d, -1);
                    int kp = ZI(p, a), km = ZI(m, a);
                    sr += ar[kp] + ar[km] - 2.0 * ar[k];
                    si += ai[kp] + ai[km] - 2.0 * ai[k];
                }
                wr += EPS1 * sr;
                wi += EPS1 * si;
            }

            /* 2:1 channels: for beta = 2*alpha (1-based: a+1 and 2(a+1))
             * H ~ eps2 ( z_a^{*2} z_b + c.c. ) with b = 2a (1-based)
             * Contribution to mode a: 2 eps2 z_a* z_b
             * Contribution to mode b:   eps2 z_a^2
             */
            if (EPS2 != 0.0 && M >= 2) {
                int a1 = a + 1; /* 1-based index */
                /* as lower mode of a 2:1 pair: b = 2*a1 */
                if (2 * a1 <= M) {
                    int b = 2 * a1 - 1; /* 0-based */
                    for (int d = 0; d < D; d++) {
                        int p = nbr(i, d, +1), m = nbr(i, d, -1);
                        for (int t = 0; t < 2; t++) {
                            int j = (t == 0) ? p : m;
                            int kb = ZI(j, b);
                            /* 2 eps2 conj(z_a) * z_b  but z_a is local:
                             * use local z_a with neighbour z_b */
                            /* term in i dz_a = 2 eps2 z_a* z_b^{nbr}  is wrong
                             * canonical: H = eps2 (z_a*^2 z_b + z_a^2 z_b*)
                             * i dz_a = 2 eps2 z_a* z_b  (same cell form)
                             * For inter-cell we use
                             * H = eps2 (z_i,a*^2 z_j,b + z_j,a*^2 z_i,b)/2 + c.c. simplified:
                             * couple local lower to neighbour higher: */
                            double zbr = ar[kb], zbi = ai[kb];
                            /* i d z_a += 2 eps2 conj(z_a) * z_b  ? for H=eps2 z_a*^2 z_b
                             * ∂H/∂z_a* = 2 eps2 z_a* z_b  => i ż_a = 2 eps2 z_a* z_b
                             * Here z_b is neighbour's. Use local z_a. */
                            double car = ar[k], cai = -ai[k]; /* conj z_a */
                            /* 2 eps2 * conj(z_a) * z_b */
                            double tr = 2.0 * EPS2 * (car * zbr - cai * zbi);
                            double ti = 2.0 * EPS2 * (car * zbi + cai * zbr);
                            wr += tr;
                            wi += ti;
                        }
                    }
                }
                /* as higher mode: a1 even, lower = a1/2 */
                if (a1 % 2 == 0) {
                    int lo = a1 / 2 - 1; /* 0-based lower */
                    for (int d = 0; d < D; d++) {
                        int p = nbr(i, d, +1), m = nbr(i, d, -1);
                        for (int t = 0; t < 2; t++) {
                            int j = (t == 0) ? p : m;
                            int kl = ZI(j, lo);
                            /* i dz_b = eps2 z_a^2  with a=lo neighbour */
                            double zr2 = ar[kl]*ar[kl] - ai[kl]*ai[kl];
                            double zi2 = 2.0*ar[kl]*ai[kl];
                            wr += EPS2 * zr2;
                            wi += EPS2 * zi2;
                        }
                    }
                }
            }

            /* intra-cell 3-wave: H = g (z1 z2 z3* + c.c.)
             * i ż1 = g z2* z3?  Wait:
             * H = g (z1 z2 conj(z3) + conj(z1) conj(z2) z3)
             * ∂H/∂z1* = g conj(z2) z3
             * ∂H/∂z2* = g conj(z1) z3
             * ∂H/∂z3* = g z1 z2
             */
            if (GINT != 0.0 && M >= 3) {
                int k0 = ZI(i, 0), k1 = ZI(i, 1), k2 = ZI(i, 2);
                if (a == 0) {
                    /* g * conj(z2) * z3 */
                    double c2r = ar[k1], c2i = -ai[k1];
                    double t3r = ar[k2], t3i = ai[k2];
                    wr += GINT * (c2r * t3r - c2i * t3i);
                    wi += GINT * (c2r * t3i + c2i * t3r);
                } else if (a == 1) {
                    double c1r = ar[k0], c1i = -ai[k0];
                    double t3r = ar[k2], t3i = ai[k2];
                    wr += GINT * (c1r * t3r - c1i * t3i);
                    wi += GINT * (c1r * t3i + c1i * t3r);
                } else if (a == 2) {
                    /* g * z1 * z2 */
                    wr += GINT * (ar[k0]*ar[k1] - ai[k0]*ai[k1]);
                    wi += GINT * (ar[k0]*ai[k1] + ai[k0]*ar[k1]);
                }
            }

            /* dz/dt = -i * (wr + i wi) = wi - i wr */
            br[k] =  wi;
            bi[k] = -wr;
        }
        (void)Il;
    }
}

static void step_theta(double noise_amp) {
    cache_I();
    static double *nth; nth = realloc(nth, sizeof(double) * N);
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < N; i++) {
        double Il[MMAX];
        for (int a = 0; a < M; a++) {
            int k = ZI(i, a);
            Il[a] = zr[k]*zr[k] + zi[k]*zi[k];
        }
        double teq = theta_eq_from_I(Il);
        double lap = 0;
        for (int d = 0; d < D; d++) {
            int p = nbr(i, d, +1), m = nbr(i, d, -1);
            lap += th[p] + th[m] - 2.0 * th[i];
        }
        double force = -KTH * (th[i] - teq) + MU * lap;
        double nse = (noise_amp > 0) ? noise_amp * grand() : 0.0;
        nth[i] = th[i] + DT * (ETA * force + nse);
        /* soft clamp */
        if (nth[i] >  1.5) nth[i] =  1.5;
        if (nth[i] < -1.5) nth[i] = -1.5;
    }
    memcpy(th, nth, sizeof(double) * N);
}

/* RK2 on z, then theta */
static void step(double noise_z) {
    static double *kr, *ki, *tr, *ti;
    size_t nm = (size_t)N * M;
    kr = realloc(kr, sizeof(double) * nm);
    ki = realloc(ki, sizeof(double) * nm);
    tr = realloc(tr, sizeof(double) * nm);
    ti = realloc(ti, sizeof(double) * nm);

    deriv_z(zr, zi, th, kr, ki);
    #pragma omp parallel for schedule(static)
    for (size_t q = 0; q < nm; q++) {
        tr[q] = zr[q] + 0.5 * DT * kr[q];
        ti[q] = zi[q] + 0.5 * DT * ki[q];
    }
    deriv_z(tr, ti, th, kr, ki);
    #pragma omp parallel for schedule(static)
    for (size_t q = 0; q < nm; q++) {
        zr[q] += DT * kr[q];
        zi[q] += DT * ki[q];
        if (noise_z > 0) {
            zr[q] += noise_z * sqrt(DT) * grand();
            zi[q] += noise_z * sqrt(DT) * grand();
        }
    }
    step_theta(noise_z > 0 ? 0.3 * noise_z : 0.0);
}

static double total_action(void) {
    cache_I();
    double s = 0;
    for (int i = 0; i < N; i++) s += Isum[i];
    return s;
}

static double peak_action(void) {
    cache_I();
    double p = 0;
    for (int i = 0; i < N; i++) if (Isum[i] > p) p = Isum[i];
    return p;
}

/* free + geometric part of H (interaction energy omitted — O(eps) diagnostic) */
static double energy_diag(void) {
    cache_I();
    double e = 0;
    for (int i = 0; i < N; i++) {
        double Il[MMAX];
        for (int a = 0; a < M; a++) {
            int k = ZI(i, a);
            Il[a] = zr[k]*zr[k] + zi[k]*zi[k];
            e += w_alpha(a, th[i]) * Il[a];
        }
        double teq = theta_eq_from_I(Il);
        double dth = th[i] - teq;
        e += 0.5 * KTH * dth * dth;
        for (int d = 0; d < D; d++) {
            int p = nbr(i, d, +1);
            double dt = th[i] - th[p];
            e += 0.25 * MU * dt * dt; /* half-count each edge */
        }
    }
    return e;
}

/* ---- clustering + fingerprints ---- */
typedef struct {
    int n_c;
    double s[MMAX];     /* fractional action per mode */
    int lock[MMAX];     /* best small integer lock vector */
    double lock_res;    /* residual |n·ω| */
    double I_tot;
    double th_mean;
} Finger;

static int find_clusters(int *label, int *sizes, int maxc, Finger *fp) {
    cache_I();
    double pk = 0;
    for (int i = 0; i < N; i++) if (Isum[i] > pk) pk = Isum[i];
    double thr = I_THR_FRAC * pk;
    if (thr < 1e-8) thr = 1e-8;
    for (int i = 0; i < N; i++) label[i] = -1;
    int *stack = malloc(sizeof(int) * N);
    int nc = 0;
    for (int i = 0; i < N; i++) {
        if (Isum[i] < thr || label[i] >= 0) continue;
        int sp = 0, cnt = 0;
        stack[sp++] = i; label[i] = nc;
        double modeI[MMAX]; memset(modeI, 0, sizeof(modeI));
        double ths = 0;
        while (sp) {
            int c = stack[--sp]; cnt++;
            ths += th[c];
            for (int a = 0; a < M; a++) {
                int k = ZI(c, a);
                modeI[a] += zr[k]*zr[k] + zi[k]*zi[k];
            }
            for (int d = 0; d < D; d++)
                for (int s = -1; s <= 1; s += 2) {
                    int j = nbr(c, d, s);
                    if (Isum[j] >= thr && label[j] < 0) {
                        label[j] = nc; stack[sp++] = j;
                    }
                }
        }
        if (nc < maxc) {
            sizes[nc] = cnt;
            Finger *f = &fp[nc];
            f->n_c = cnt;
            f->I_tot = 0;
            for (int a = 0; a < M; a++) f->I_tot += modeI[a];
            for (int a = 0; a < M; a++)
                f->s[a] = (f->I_tot > 1e-15) ? modeI[a] / f->I_tot : 0;
            f->th_mean = ths / (cnt > 0 ? cnt : 1);

            /* best lock vector |n_a|<=2, not all zero: minimise |sum n_a w_a|
             * using mean theta of cluster */
            double best = 1e99;
            int bn[MMAX]; memset(bn, 0, sizeof(bn));
            double thm = f->th_mean;
            double om[MMAX];
            for (int a = 0; a < M; a++) om[a] = w_alpha(a, thm);
            /* exhaustive small search */
            int lim = 2;
            int ntot = 1; for (int a = 0; a < M; a++) ntot *= (2*lim+1);
            for (int code = 0; code < ntot; code++) {
                int n[MMAX], c = code, any = 0;
                for (int a = 0; a < M; a++) {
                    n[a] = (c % (2*lim+1)) - lim; c /= (2*lim+1);
                    if (n[a]) any = 1;
                }
                if (!any) continue;
                /* gcd-normalise skip; just score */
                double r = 0;
                for (int a = 0; a < M; a++) r += n[a] * om[a];
                r = fabs(r);
                if (r < best) { best = r; memcpy(bn, n, sizeof(bn)); }
            }
            memcpy(f->lock, bn, sizeof(int)*MMAX);
            f->lock_res = best;
        }
        nc++;
    }
    free(stack);
    return nc;
}

/* coarse species bin: N_c bin + mode occupancy pattern (3 bits each s_a) */
static unsigned species_key(const Finger *f) {
    int nb = f->n_c;
    if (nb > 63) nb = 63;
    unsigned key = (unsigned)nb;
    for (int a = 0; a < M && a < 3; a++) {
        int b = (int)floor(f->s[a] * 4.0); /* 0..3 */
        if (b < 0) b = 0; if (b > 3) b = 3;
        key |= ((unsigned)b) << (6 + 2*a);
    }
    /* lock signature: which modes participate in lock */
    unsigned lp = 0;
    for (int a = 0; a < M && a < 3; a++)
        if (f->lock[a]) lp |= (1u << a);
    key |= lp << 12;
    return key;
}

static void report_clusters(const char *tag, double t) {
    int *lab = malloc(sizeof(int) * N);
    int *sz = malloc(sizeof(int) * MAX_CLU);
    Finger *fp = calloc(MAX_CLU, sizeof(Finger));
    int nc = find_clusters(lab, sz, MAX_CLU, fp);
    int cap = nc < MAX_CLU ? nc : MAX_CLU;

    double mean_nc = 0, var_nc = 0;
    for (int i = 0; i < cap; i++) mean_nc += fp[i].n_c;
    if (cap) mean_nc /= cap;
    for (int i = 0; i < cap; i++) {
        double d = fp[i].n_c - mean_nc;
        var_nc += d * d;
    }
    if (cap) var_nc = sqrt(var_nc / cap);

    /* histogram top species keys */
    unsigned keys[MAX_CLU];
    int counts[MAX_CLU]; int nk = 0;
    for (int i = 0; i < cap; i++) {
        unsigned k = species_key(&fp[i]);
        int found = -1;
        for (int j = 0; j < nk; j++) if (keys[j] == k) { found = j; break; }
        if (found >= 0) counts[found]++;
        else if (nk < MAX_CLU) { keys[nk] = k; counts[nk] = 1; nk++; }
    }
    /* sort by count */
    for (int i = 0; i < nk; i++)
        for (int j = i+1; j < nk; j++)
            if (counts[j] > counts[i]) {
                int tc = counts[i]; counts[i] = counts[j]; counts[j] = tc;
                unsigned tk = keys[i]; keys[i] = keys[j]; keys[j] = tk;
            }

    double act = total_action(), pk = peak_action(), en = energy_diag();
    printf("  [%s t=%.1f] n_lumps=%d  mean_Nc=%.2f  sd/mean=%.3f  act=%.3f  E~%.3f  peakI=%.3f\n",
           tag, t, nc, mean_nc, mean_nc > 0 ? var_nc / mean_nc : 0.0, act, en, pk);
    int show = nk < 8 ? nk : 8;
    for (int i = 0; i < show; i++) {
        /* decode a representative */
        Finger *rep = NULL;
        for (int j = 0; j < cap; j++)
            if (species_key(&fp[j]) == keys[i]) { rep = &fp[j]; break; }
        if (!rep) continue;
        printf("    species#%d count=%d  Nc~%d  s=(", i, counts[i], rep->n_c);
        for (int a = 0; a < M; a++)
            printf("%s%.2f", a ? "," : "", rep->s[a]);
        printf(")  lock=(");
        for (int a = 0; a < M; a++)
            printf("%s%d", a ? "," : "", rep->lock[a]);
        printf(")  |n·ω|=%.4f  θ=%.3f\n", rep->lock_res, rep->th_mean);
    }
    fflush(stdout);
    free(lab); free(sz); free(fp);
}

static void init_random(double amp) {
    /* sparse random: only a fraction of cells start hot, so multiple
       seeds can compete instead of one percolating sea */
    double fill = 0.12;
    for (int i = 0; i < N; i++) {
        th[i] = 0.0;
        int hot = (urand() < fill);
        for (int a = 0; a < M; a++) {
            int k = ZI(i, a);
            if (hot) {
                zr[k] = amp * (urand() * 2 - 1);
                zi[k] = amp * (urand() * 2 - 1);
            } else {
                zr[k] = 0.05 * amp * (urand() * 2 - 1);
                zi[k] = 0.05 * amp * (urand() * 2 - 1);
            }
        }
    }
}

static void init_seeded_lumps(int nseeds, double sep) {
    /* zero field, place multi-mode Gaussian-ish lumps */
    size_t nm = (size_t)N * M;
    memset(zr, 0, sizeof(double) * nm);
    memset(zi, 0, sizeof(double) * nm);
    memset(th, 0, sizeof(double) * N);
    int c0 = L / 2;
    for (int s = 0; s < nseeds; s++) {
        int shift = (int)((s - (nseeds - 1) * 0.5) * sep);
        int cx = (c0 + shift + L) % L;
        int cy = c0, cz = c0;
        for (int i = 0; i < N; i++) {
            int x = i % L;
            int y = (i / stride[1 % D]) % L;
            int z = (D > 2) ? (i / stride[2]) % L : 0;
            int dx = x - cx; if (dx > L/2) dx -= L; if (dx < -L/2) dx += L;
            int dy = y - cy; if (dy > L/2) dy -= L; if (dy < -L/2) dy += L;
            int dz = z - cz; if (dz > L/2) dz -= L; if (dz < -L/2) dz += L;
            double r2 = dx*dx + dy*dy + (D > 2 ? dz*dz : 0);
            double env = exp(-0.5 * r2 / 2.25);
            if (env < 1e-4) continue;
            /* S3-like: modes 1,2,3 with phases for 1+2~3 lock */
            double ph = 0.3 * s;
            zr[ZI(i,0)] += env * 0.55 * cos(ph);
            zi[ZI(i,0)] += env * 0.55 * sin(ph);
            if (M > 1) {
                zr[ZI(i,1)] += env * 0.40 * cos(1.7*ph);
                zi[ZI(i,1)] += env * 0.40 * sin(1.7*ph);
            }
            if (M > 2) {
                zr[ZI(i,2)] += env * 0.30 * cos(2.9*ph);
                zi[ZI(i,2)] += env * 0.30 * sin(2.9*ph);
            }
            th[i] += -0.4 * env;
        }
    }
}

static void alloc_fields(void) {
    N = 1; for (int a = 0; a < D; a++) N *= L;
    stride[0] = 1;
    for (int a = 1; a < D; a++) stride[a] = stride[a-1] * L;
    size_t nm = (size_t)N * M;
    zr = realloc(zr, sizeof(double) * nm);
    zi = realloc(zi, sizeof(double) * nm);
    th = realloc(th, sizeof(double) * N);
    Isum = realloc(Isum, sizeof(double) * N);
    set_spectrum();
}

static void choose_L(void) {
    /* keep N roughly 8k-20k */
    int Ls[] = {0, 4096, 90, 22, 12};
    if (D >= 1 && D <= 4) L = Ls[D];
    if (D == 1) L = 4096;
    if (D == 2) L = 90;
    if (D == 3) L = 22;
    if (D == 4) L = 12;
}

static void run_steps(int nsteps, double noise, int report_every, const char *tag) {
    for (int t = 1; t <= nsteps; t++) {
        step(noise);
        if (report_every > 0 && (t % report_every == 0 || t == nsteps))
            report_clusters(tag, t * DT);
    }
}

static void banner(const char *sched) {
    printf("=====================================================================\n");
    printf("v88 multi-harmonic fabric  schedule=%s  d=%d L=%d N=%d M=%d\n",
           sched, D, L, N, M);
    printf("=====================================================================\n");
    printf("  omega_bar =");
    for (int a = 0; a < M; a++) printf(" %.4f", WBAR[a]);
    printf("\n  sig=%.2f thmax=%.2f gam=%.2f eps1=%.3f eps2=%.3f gint=%.3f mu=%.2f\n",
           SIG, THMAX, GAM, EPS1, EPS2, GINT, MU);
    printf("  Integers enter via mode indices, p:q channels, lock vectors.\n");
    printf("  Not a DNLS breather: m>=2 + geometry + multi-channel resonance.\n\n");
}

int main(int argc, char **argv) {
    const char *sched = (argc > 1) ? argv[1] : "smoke";
    D = (argc > 2) ? atoi(argv[2]) : 3;
    unsigned long long seed = (argc > 3) ? strtoull(argv[3], NULL, 10)
                                         : 0xC0FFEEULL;
    if (D < 1) D = 1; if (D > DMAX) D = DMAX;
    rs = seed ? seed : 0xC0FFEEULL;

    /* detune control: kill multi-mode channels */
    if (strcmp(sched, "detune") == 0) {
        EPS2 = 0.0;
        GINT = 0.0;
    }

    choose_L();
    alloc_fields();
    banner(sched);

    if (strcmp(sched, "smoke") == 0) {
        init_random(0.35);
        printf("  SMOKE: short evolution, expect finite action and some lumps\n");
        report_clusters("init", 0);
        run_steps(800, 0.0, 400, "smoke");
    } else if (strcmp(sched, "quench") == 0) {
        init_random(0.45);
        report_clusters("init", 0);
        run_steps((int)(400.0 / DT), 0.0, (int)(100.0 / DT), "quench");
    } else if (strcmp(sched, "thermal") == 0) {
        init_random(0.45);
        report_clusters("init", 0);
        /* cool noise: start 0.04 -> 0 */
        int total = (int)(600.0 / DT);
        int chunk = total / 6;
        for (int stage = 0; stage < 6; stage++) {
            double nse = 0.04 * (5 - stage) / 5.0;
            char tag[32]; snprintf(tag, sizeof tag, "therm_s%d", stage);
            run_steps(chunk, nse, chunk, tag);
        }
    } else if (strcmp(sched, "cyclic") == 0) {
        init_random(0.45);
        for (int cyc = 0; cyc < 3; cyc++) {
            char tag[32];
            snprintf(tag, sizeof tag, "cyc%d_hot", cyc);
            run_steps((int)(80.0 / DT), 0.035, (int)(80.0 / DT), tag);
            snprintf(tag, sizeof tag, "cyc%d_cool", cyc);
            run_steps((int)(150.0 / DT), 0.0, (int)(150.0 / DT), tag);
        }
    } else if (strcmp(sched, "driven") == 0) {
        init_random(0.25);
        /* continuous weak injection into mode 0 */
        int total = (int)(400.0 / DT);
        int rep = (int)(100.0 / DT);
        for (int t = 1; t <= total; t++) {
            step(0.0);
            /* pump */
            double amp = 0.002;
            for (int i = 0; i < N; i++) {
                int k = ZI(i, 0);
                zr[k] += amp * DT * cos(WBAR[0] * t * DT);
                zi[k] += amp * DT * sin(WBAR[0] * t * DT);
            }
            if (t % rep == 0 || t == total)
                report_clusters("driven", t * DT);
        }
    } else if (strcmp(sched, "staged") == 0) {
        /* raise action density in stages by rescaling */
        init_random(0.20);
        for (int st = 0; st < 4; st++) {
            double scale = 1.0 + 0.35 * st;
            size_t nm = (size_t)N * M;
            for (size_t q = 0; q < nm; q++) { zr[q] *= scale; zi[q] *= scale; }
            /* undo cumulative: only scale relative to previous once —
               simpler: re-init with increasing amp */
            (void)scale;
        }
        /* cleaner staged: re-init each stage with higher amp, anneal */
        double amps[] = {0.25, 0.35, 0.45, 0.55};
        for (int st = 0; st < 4; st++) {
            init_random(amps[st]);
            char tag[32]; snprintf(tag, sizeof tag, "stage%d", st);
            run_steps((int)(200.0 / DT), 0.01, (int)(200.0 / DT), tag);
            run_steps((int)(150.0 / DT), 0.0, (int)(150.0 / DT), tag);
        }
    } else if (strcmp(sched, "coalesce") == 0) {
        init_seeded_lumps(3, 5.0);
        report_clusters("seeds", 0);
        run_steps((int)(300.0 / DT), 0.0, (int)(100.0 / DT), "coalesce");
    } else if (strcmp(sched, "detune") == 0) {
        /* EPS2=GINT already zeroed */
        init_random(0.45);
        report_clusters("init", 0);
        run_steps((int)(400.0 / DT), 0.0, (int)(100.0 / DT), "detune");
    } else {
        fprintf(stderr, "unknown schedule '%s'\n", sched);
        fprintf(stderr, "use: smoke|quench|thermal|cyclic|driven|staged|coalesce|detune\n");
        return 1;
    }

    printf("\n  DONE schedule=%s seed=%llu action=%.6f peakI=%.6f\n",
           sched, (unsigned long long)seed, total_action(), peak_action());
    printf("  Species recurrence across seeds is the success criterion.\n");
    printf("  Compare quench/thermal multi-peak fingerprints to detune control.\n");
    return 0;
}
