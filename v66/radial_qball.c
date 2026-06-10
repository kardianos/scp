/* v66 radial_qball.c — static radial Q-ball profile solver (SPEC.md §8, THEORY.md §3-4)
 *
 * ODE (symmetric ansatz Phi_a = f(r) e^{i w t}, Theta = 0):
 *     f'' + (2/r) f' = U'(f)
 *     U'(f)  = (m^2 - w^2) f + 2 Vt'(f^6) f^5,   Vt'(s) = (mu/2)/(1 + kappa s)^2
 *     U(f)   = (1/2)(m^2 - w^2) f^2 + (mu/6) f^6/(1 + kappa f^6)
 *     f'(0) = 0, f(inf) = 0
 * Existence window: w^2 > wmin^2 = m^2 + (mu/9)(2/kappa)^(2/3), w < m.
 *
 * Method: long-double RK4 shooting on f0 in (f1, f_dip) of U_w (under/over
 * bisection, matching theory/qball_checks.py classification); series start at
 * r=0 for regularity; analytic tail C e^{-beta r}/r attached where f < 1e-3 f0;
 * composite Simpson integrals with 4 pi r^2 volume weight.
 *
 * Observables (THEORY C10, C11):
 *     Q = 3 w  int f^2 4 pi r^2 dr
 *     E = int [ (3/2)(w^2 f^2 + f'^2 + m^2 f^2) + (mu/2) f^6/(1+kappa f^6) ] 4 pi r^2 dr
 *
 * Build: gcc -O3 -o radial_qball radial_qball.c -lm
 * Usage: ./radial_qball -omega 1.39 [options]                     (single profile)
 *        ./radial_qball -omega_min 1.315 -omega_max 1.495 -omega_step 0.005   (scan)
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>

typedef long double ld;

/* ---- parameters (global; set from CLI) ---- */
static ld M2  = 2.25L;
static ld MU  = -41.345L;
static ld KAP = 50.0L;

static ld Vt (ld s) { return 0.5L*MU*s/(1.0L + KAP*s); }
static ld Vtp(ld s) { ld d = 1.0L + KAP*s; return 0.5L*MU/(d*d); }
static ld Uw (ld f, ld w2) { ld f2=f*f, f6=f2*f2*f2; return 0.5L*(M2-w2)*f2 + Vt(f6)/3.0L; }
static ld Uwp(ld f, ld w2) { ld f2=f*f, f4=f2*f2;    return (M2-w2)*f + 2.0L*Vtp(f4*f2)*f4*f; }

/* ---- RK4 step for (f, p=f'), p' = U'(f) - (2/r) p ---- */
static void rk4(ld *f, ld *p, ld r, ld dr, ld w2)
{
    ld k1f = *p,             k1p = Uwp(*f, w2) - 2.0L/r*(*p);
    ld r2  = r + 0.5L*dr;
    ld f2  = *f + 0.5L*dr*k1f, p2 = *p + 0.5L*dr*k1p;
    ld k2f = p2,             k2p = Uwp(f2, w2) - 2.0L/r2*p2;
    ld f3  = *f + 0.5L*dr*k2f, p3 = *p + 0.5L*dr*k2p;
    ld k3f = p3,             k3p = Uwp(f3, w2) - 2.0L/r2*p3;
    ld r4  = r + dr;
    ld f4  = *f + dr*k3f,    p4 = *p + dr*k3p;
    ld k4f = p4,             k4p = Uwp(f4, w2) - 2.0L/r4*p4;
    *f += dr/6.0L*(k1f + 2.0L*k2f + 2.0L*k3f + k4f);
    *p += dr/6.0L*(k1p + 2.0L*k2p + 2.0L*k3p + k4p);
}

/* ---- bracket auto-detection ---- */
/* f1: smallest positive zero of U (entering U<0), bisected */
static ld find_f1(ld w2)
{
    ld fprev = 1e-4L;
    for (ld fg = 2e-4L; fg <= 3.0L; fg += 1e-4L) {
        if (Uw(fg, w2) < 0) {
            ld a = fprev, b = fg;             /* U(a)>=0, U(b)<0 */
            for (int i = 0; i < 200; i++) {
                ld m = 0.5L*(a+b);
                if (Uw(m, w2) < 0) b = m; else a = m;
            }
            return b;                          /* just inside U<0 */
        }
        fprev = fg;
    }
    return -1.0L;
}
/* upper end of the U<0 region (for dip bracketing) */
static ld find_f2(ld f1, ld w2)
{
    for (ld fg = f1 + 1e-4L; fg <= 5.0L; fg += 1e-4L)
        if (Uw(fg, w2) >= 0) return fg;
    return 5.0L;
}
/* f_dip = argmin U: bisect U'(f)=0 in (f1, f2); return the U'<0 side (matches qball_checks.py) */
static ld find_fdip(ld f1, ld f2, ld w2)
{
    ld a = f1, b = f2;
    for (int i = 0; i < 200; i++) {
        ld m = 0.5L*(a+b);
        if (Uwp(m, w2) < 0) a = m; else b = m;
    }
    return a;
}

/* ---- shooting classification (matches theory/qball_checks.py) ----
 * UNDER: f' > 0 while f < 0.99 f0 (friction won)        -> raise f0
 * OVER : f crosses 0, or runs away above 1.05 f0, or reaches r_stop -> lower f0
 */
enum { TAG_UNDER = 0, TAG_OVER = 1 };

static int classify(ld f0, ld w2, ld dr, ld rstop)
{
    ld up0 = Uwp(f0, w2);
    ld f = f0 + dr*dr/6.0L*up0;     /* series start, f'(0)=0 regularity */
    ld p = dr/3.0L*up0;
    ld r = dr;
    long n = (long)(rstop/dr + 0.5);
    for (long i = 1; i < n; i++) {
        rk4(&f, &p, r, dr, w2);
        r += dr;
        if (f < 0)                    return TAG_OVER;
        if (f > 1.05L*f0)             return TAG_OVER;   /* 'right' */
        if (p > 0 && f < 0.99L*f0)    return TAG_UNDER;
    }
    return TAG_OVER;                                      /* 'end' */
}

/* ---- assembled profile (numerical trace + analytic tail) ---- */
typedef struct {
    double *tf, *tp;     /* trace on i*dr grid, valid 0..i_last */
    long    i_last;
    long    i_attach;    /* tail attach index */
    ld      dr, r_t, f_t, beta;
} Prof;

static void prof_free(Prof *P) { free(P->tf); free(P->tp); P->tf = P->tp = NULL; }

static void prof_eval(const Prof *P, ld r, ld *f, ld *p)
{
    if (r > P->r_t || P->i_attach < 1) {
        if (P->f_t <= 0) { *f = 0; *p = 0; return; }
        ld e = P->f_t * (P->r_t / r) * expl(-P->beta*(r - P->r_t));
        *f = e;
        *p = -e*(P->beta + 1.0L/r);
        return;
    }
    ld x = r / P->dr;
    long i = (long)x; if (i >= P->i_attach) i = P->i_attach - 1;
    ld fr = x - (ld)i;
    *f = (1.0L-fr)*P->tf[i] + fr*P->tf[i+1];
    *p = (1.0L-fr)*P->tp[i] + fr*P->tp[i+1];
}

typedef struct {
    double w, f0, Q, E, EomQ, rhalf, rQ, residual, e0, s0;
} Obs;

#define R_STOP_DEF 80.0L

/* full solve at one omega; fills obs; fills prof if non-NULL (caller frees) */
static int solve_qball(double w, double dr_in, int iters,
                       double f0_lo_in, double f0_hi_in,
                       Obs *out, Prof *prof_out)
{
    ld w2 = (ld)w*(ld)w;
    ld dr = (ld)dr_in;
    ld rstop = R_STOP_DEF;

    /* bracket */
    ld f1 = find_f1(w2);
    if (f1 < 0) { fprintf(stderr, "ERROR: no U<0 region at omega=%.6f (outside window)\n", w); return 1; }
    ld lo, hi;
    if (f0_lo_in > 0) lo = (ld)f0_lo_in;
    else              lo = f1*(1.0L + 1e-6L);
    if (f0_hi_in > 0) hi = (ld)f0_hi_in;
    else              hi = find_fdip(f1, find_f2(f1, w2), w2);
    if (classify(lo, w2, dr, rstop) != TAG_UNDER) {
        fprintf(stderr, "ERROR: lo bracket f0=%.9Lf not 'under' at omega=%.6f\n", lo, w);
        return 1;
    }
    /* bisection (long double) */
    for (int it = 0; it < iters; it++) {
        if (hi - lo < 1e-17L*lo) break;
        ld mid = 0.5L*(lo + hi);
        if (classify(mid, w2, dr, rstop) == TAG_UNDER) lo = mid; else hi = mid;
    }
    ld f0 = lo;

    /* final trace with storage */
    long n_stop = (long)(rstop/dr + 0.5);
    double *tf = (double*)malloc((size_t)(n_stop+1)*sizeof(double));
    double *tp = (double*)malloc((size_t)(n_stop+1)*sizeof(double));
    if (!tf || !tp) { fprintf(stderr, "ERROR: out of memory\n"); exit(1); }
    ld up0 = Uwp(f0, w2);
    ld f = f0 + dr*dr/6.0L*up0, p = dr/3.0L*up0, r = dr;
    tf[0] = (double)f0; tp[0] = 0.0;
    tf[1] = (double)f;  tp[1] = (double)p;
    long i_last = 1;
    for (long i = 2; i <= n_stop; i++) {
        rk4(&f, &p, r, dr, w2);
        r += dr;
        tf[i] = (double)f; tp[i] = (double)p;
        i_last = i;
        if (f <= 0 || (p >= 0 && f < 0.99L*f0)) break;   /* breakdown (bracket noise) */
    }

    /* tail attach: first index in the decaying segment with f < 1e-3 f0 */
    long i_attach = -1;
    for (long i = 1; i <= i_last; i++)
        if (tf[i] > 0 && tf[i] < 1e-3*(double)f0) { i_attach = i; break; }
    if (i_attach < 0) {
        i_attach = i_last;
        while (i_attach > 1 && tf[i_attach] <= 0) i_attach--;
    }
    Prof P;
    P.tf = tf; P.tp = tp; P.i_last = i_last; P.i_attach = i_attach;
    P.dr = dr;
    P.r_t  = (ld)i_attach*dr;
    P.f_t  = (ld)tf[i_attach]; if (P.f_t < 0) P.f_t = 0;
    P.beta = sqrtl(M2 - w2);

    /* relative ODE residual over the kept numerical segment (r>0.5, f>1e-6 f0) */
    double maxres = 0.0, maxU = 0.0;
    long i0 = (long)(0.5L/dr); if (i0 < 1) i0 = 1;
    for (long i = i0; i < i_attach; i++) {
        if (tf[i] <= 1e-6*(double)f0) break;
        double fpp = (tp[i+1] - tp[i-1]) / (2.0*(double)dr);
        double rr  = (double)i*(double)dr;
        double rhs = (double)Uwp((ld)tf[i], w2);
        double res = fpp + 2.0/rr*tp[i] - rhs;
        if (fabs(res) > maxres) maxres = fabs(res);
        if (fabs(rhs) > maxU)   maxU   = fabs(rhs);
    }
    double residual = (maxU > 0) ? maxres/maxU : 0.0;

    /* composite Simpson over assembled grid 0..rstop (n_stop even by construction) */
    if (n_stop % 2) n_stop--;
    ld IQ2 = 0, IE = 0, IR4 = 0;             /* int f^2 dV; int e dV; int r^2 f^2 dV */
    for (long i = 0; i <= n_stop; i++) {
        ld ri = (ld)i*dr, fi, pi;
        if (i <= i_attach) { fi = (ld)tf[i]; pi = (ld)tp[i]; }
        else               prof_eval(&P, ri, &fi, &pi);
        ld wgt = (i == 0 || i == n_stop) ? 1.0L : (i % 2 ? 4.0L : 2.0L);
        ld r2  = ri*ri, f2 = fi*fi, f6 = f2*f2*f2;
        ld vol = 4.0L*(ld)M_PI*r2;
        IQ2 += wgt * vol * f2;
        IE  += wgt * vol * (1.5L*(w2*f2 + pi*pi + M2*f2) + Vt(f6));
        IR4 += wgt * vol * r2 * f2;
    }
    IQ2 *= dr/3.0L; IE *= dr/3.0L; IR4 *= dr/3.0L;
    /* analytic remainder beyond rstop: ~e^{-2 beta rstop}, negligible (<1e-30) */

    ld Q = 3.0L*(ld)w*IQ2;
    ld E = IE;

    /* r_half: largest r with f >= f0/2 (monotone decreasing core; linear interp) */
    double rhalf = 0.0, h = 0.5*(double)f0;
    for (long i = 0; i < i_attach; i++) {
        if (tf[i] >= h && tf[i+1] < h) {
            rhalf = ((double)i + (tf[i]-h)/(tf[i]-tf[i+1])) * (double)dr;
            /* keep scanning: spec says LARGEST r (handles non-monotone noise) */
        }
    }
    ld f0_2 = f0*f0, f0_6 = f0_2*f0_2*f0_2;
    out->w   = w;
    out->f0  = (double)f0;
    out->Q   = (double)Q;
    out->E   = (double)E;
    out->EomQ= (double)(E/(sqrtl(M2)*Q));
    out->rhalf = rhalf;
    out->rQ  = (double)sqrtl(IR4/IQ2);
    out->residual = residual;
    out->e0  = (double)(1.5L*(w2+M2)*f0_2 + Vt(f0_6));
    out->s0  = (double)f0_6;

    if (prof_out) *prof_out = P;
    else          prof_free(&P);
    return 0;
}

/* ---- output helpers ---- */
static void mkdirs_for(const char *path)
{
    char tmp[1024];
    strncpy(tmp, path, sizeof(tmp)-1); tmp[sizeof(tmp)-1] = 0;
    for (char *q = tmp + 1; *q; q++)
        if (*q == '/') { *q = 0; mkdir(tmp, 0755); *q = '/'; }
}

static void write_profile(const char *path, const Obs *o, const Prof *P,
                          double rmax, double out_dr)
{
    mkdirs_for(path);
    FILE *fp = fopen(path, "w");
    if (!fp) { fprintf(stderr, "ERROR: cannot write %s: %s\n", path, strerror(errno)); exit(1); }
    fprintf(fp, "# radial_qball profile: omega=%.6f m2=%.6f mu=%.6f kappa=%.6f\n",
            o->w, (double)M2, (double)MU, (double)KAP);
    fprintf(fp, "# f0=%.6f  E=%.1f  Q=%.1f  E/(m*Q)=%.4f  residual=%.1e\n",
            o->f0, o->E, o->Q, o->EomQ, o->residual);
    fprintf(fp, "# r f\n");
    long n_out = (long)(rmax/out_dr + 0.5);
    for (long i = 0; i <= n_out; i++) {
        ld r = (ld)i*(ld)out_dr, f, p;
        if (i == 0) f = (ld)P->tf[0];
        else prof_eval(P, r, &f, &p);
        if (f < 0) f = 0;
        fprintf(fp, "%.6f %.9f\n", (double)r, (double)f);
    }
    fclose(fp);
    printf("wrote profile %s (rmax=%.1f, out_dr=%.3f)\n", path, rmax, out_dr);
}

static void scan_append(const char *path, const Obs *o, int n, const int *sign)
{
    mkdirs_for(path);
    struct stat st;
    int fresh = (stat(path, &st) != 0);
    FILE *fp = fopen(path, "a");
    if (!fp) { fprintf(stderr, "ERROR: cannot write %s: %s\n", path, strerror(errno)); exit(1); }
    if (fresh)
        fprintf(fp, "omega\tf0\tE\tQ\tE_over_mQ\tdQ_domega_sign\tr_half\tr_Q\tresidual\n");
    for (int i = 0; i < n; i++)
        fprintf(fp, "%.4f\t%.6f\t%.1f\t%.1f\t%.4f\t%d\t%.2f\t%.2f\t%.1e\n",
                o[i].w, o[i].f0, o[i].E, o[i].Q, o[i].EomQ, sign[i],
                o[i].rhalf, o[i].rQ, o[i].residual);
    fclose(fp);
    printf("appended %d row%s to %s\n", n, n == 1 ? "" : "s", path);
}

/* sanity gates vs THEORY.md §4 (omega = 1.39 reference row) */
static void sanity_check(const Obs *o)
{
    if (fabs(o->w - 1.39) > 1e-6) return;
    struct { const char *name; double got, want; } g[] = {
        { "f0",     o->f0,    0.6405 },
        { "Q",      o->Q,     482.2  },
        { "E",      o->E,     691.9  },
        { "r_half", o->rhalf, 4.41   },
    };
    for (int i = 0; i < 4; i++) {
        double rel = fabs(g[i].got - g[i].want)/g[i].want;
        if (rel > 0.05)
            fprintf(stderr, "WARN: omega=1.39 sanity gate: %s=%.4f vs THEORY %.4f (%.1f%% off)\n",
                    g[i].name, g[i].got, g[i].want, 100.0*rel);
    }
}

int main(int argc, char **argv)
{
    double omega = -1, omega_min = -1, omega_max = -1, omega_step = -1;
    double rmax = 60.0, dr = 0.001, out_dr = 0.02;
    double f0_lo = -1, f0_hi = -1;
    int iters = 120;
    char profile_path[1024] = "";
    char scan_path[1024] = "v66/results/scan.tsv";

    for (int i = 1; i < argc - 1; i += 2) {
        const char *k = argv[i], *v = argv[i+1];
        if      (!strcmp(k, "-omega"))      omega = atof(v);
        else if (!strcmp(k, "-omega_min"))  omega_min = atof(v);
        else if (!strcmp(k, "-omega_max"))  omega_max = atof(v);
        else if (!strcmp(k, "-omega_step")) omega_step = atof(v);
        else if (!strcmp(k, "-m2"))         M2 = (ld)atof(v);
        else if (!strcmp(k, "-mu"))         MU = (ld)atof(v);
        else if (!strcmp(k, "-kappa"))      KAP = (ld)atof(v);
        else if (!strcmp(k, "-rmax"))       rmax = atof(v);
        else if (!strcmp(k, "-dr"))         dr = atof(v);
        else if (!strcmp(k, "-out_dr"))     out_dr = atof(v);
        else if (!strcmp(k, "-f0_lo"))      f0_lo = atof(v);
        else if (!strcmp(k, "-f0_hi"))      f0_hi = atof(v);
        else if (!strcmp(k, "-iters"))      iters = atoi(v);
        else if (!strcmp(k, "-profile"))    { strncpy(profile_path, v, sizeof(profile_path)-1); }
        else if (!strcmp(k, "-scan"))       { strncpy(scan_path, v, sizeof(scan_path)-1); }
        else { fprintf(stderr, "ERROR: unknown flag %s\n", k); return 1; }
    }

    /* existence window */
    double m = (double)sqrtl(M2);
    double wmin2 = (double)(M2 + MU/9.0L*powl(2.0L/KAP, 2.0L/3.0L));
    if (wmin2 >= (double)M2 || wmin2 < 0) {
        fprintf(stderr, "ERROR: no existence window for m2=%.4f mu=%.4f kappa=%.4f\n",
                (double)M2, (double)MU, (double)KAP);
        return 1;
    }
    double wmin = sqrt(wmin2);
    printf("radial_qball: m2=%.4f mu=%.4f kappa=%.4f -> window omega in (%.6f, %.6f)\n",
           (double)M2, (double)MU, (double)KAP, wmin, m);

    int scan_mode = (omega_min > 0 && omega_max > 0 && omega_step > 0);
    if (!scan_mode && omega <= 0) {
        fprintf(stderr, "usage: radial_qball -omega W | -omega_min A -omega_max B -omega_step S  [options]\n");
        return 1;
    }

    if (scan_mode) {
        int n = (int)((omega_max - omega_min)/omega_step + 1e-9) + 1;
        Obs *obs = (Obs*)malloc((size_t)n*sizeof(Obs));
        int *sgn = (int*)malloc((size_t)n*sizeof(int));
        for (int i = 0; i < n; i++) {
            double w = omega_min + i*omega_step;
            if (w <= wmin || w >= m) {
                fprintf(stderr, "ERROR: omega=%.4f outside window (%.6f, %.6f)\n", w, wmin, m);
                return 1;
            }
            Prof P;
            if (solve_qball(w, dr, iters, f0_lo, f0_hi, &obs[i], &P)) return 1;
            printf("omega=%.4f  f0=%.6f  Q=%.2f  E=%.2f  E/Q=%.4f  r_half=%.3f  r_Q=%.3f  res=%.1e\n",
                   w, obs[i].f0, obs[i].Q, obs[i].E, obs[i].E/obs[i].Q, obs[i].rhalf,
                   obs[i].rQ, obs[i].residual);
            char ppath[1024];
            if (profile_path[0]) snprintf(ppath, sizeof(ppath), "%s", profile_path);
            else snprintf(ppath, sizeof(ppath), "v66/results/profile_omega%.4f.txt", w);
            write_profile(ppath, &obs[i], &P, rmax, out_dr);
            sanity_check(&obs[i]);
            prof_free(&P);
        }
        for (int i = 0; i < n; i++) {            /* centered dQ/domega sign */
            double dQ;
            if (i == 0)        dQ = obs[1].Q   - obs[0].Q;
            else if (i == n-1) dQ = obs[n-1].Q - obs[n-2].Q;
            else               dQ = obs[i+1].Q - obs[i-1].Q;
            sgn[i] = (dQ < 0) ? -1 : 1;
        }
        scan_append(scan_path, obs, n, sgn);
        free(obs); free(sgn);
    } else {
        if (omega <= wmin || omega >= m) {
            fprintf(stderr, "ERROR: omega=%.4f outside window (%.6f, %.6f)\n", omega, wmin, m);
            return 1;
        }
        Obs o; Prof P;
        if (solve_qball(omega, dr, iters, f0_lo, f0_hi, &o, &P)) return 1;
        printf("omega=%.4f  f0=%.6f  Q=%.2f  E=%.2f  E/Q=%.4f  r_half=%.3f  r_Q=%.3f  res=%.1e\n",
               omega, o.f0, o.Q, o.E, o.E/o.Q, o.rhalf, o.rQ, o.residual);
        /* dQ/domega sign from internal solves at omega +/- 0.002 (clamped to window) */
        double h = 0.002;
        double wlo = omega - h, whi = omega + h;
        if (wlo <= wmin) wlo = omega;
        if (whi >= m)    whi = omega;
        Obs olo = o, ohi = o;
        if (wlo != omega && solve_qball(wlo, dr, iters, -1, -1, &olo, NULL)) return 1;
        if (whi != omega && solve_qball(whi, dr, iters, -1, -1, &ohi, NULL)) return 1;
        int sgn = ((ohi.Q - olo.Q) < 0) ? -1 : 1;

        char ppath[1024];
        if (profile_path[0]) snprintf(ppath, sizeof(ppath), "%s", profile_path);
        else snprintf(ppath, sizeof(ppath), "v66/results/profile_omega%.4f.txt", omega);
        write_profile(ppath, &o, &P, rmax, out_dr);
        sanity_check(&o);
        scan_append(scan_path, &o, 1, &sgn);
        prof_free(&P);
    }
    return 0;
}
