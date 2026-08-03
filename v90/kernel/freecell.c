/* freecell.c — v90 KERNEL OF RECORD (carried from v89/freecell.c).
 *
 * v90 role: the production C kernel (GPU path preserved). The Go kernel
 * (v90/fab) is the side-by-side experiment, A/B-verified against this
 * file (v90/VERIFY.md). New apparatus lands here first, battery-gated
 * (go run ./cmd/battery). The v89 original stays frozen as evidence.
 *
 * ORIGINAL HEADER (v89, 2026-08-02) follows.
 *
 * freecell.c — v89 — FREE CELLS: the substance test.
 *
 * Standalone. Does not touch cellfab.c. Build:
 *   gcc -O2 -Wall -Wextra -o freecell freecell.c -lm
 * Run:
 *   ./freecell [config.cfg] [key=value ...]
 *
 * WHAT THIS IS
 * ------------
 * FREECELL.md's substance half: can a free-cell fabric hold a localised
 * object together? The execution half (local clocks) is solved and is
 * not re-tested here; this instrument tests COHESION on a substrate
 * whose geometry is state:
 *
 *   - cell positions are DYNAMICAL (overdamped motion under forces);
 *   - radii are live, r = r0*(Es/e_s0)^(1/3)  (as in cellfab pass 1);
 *   - the channel set is re-derived every step from current geometry
 *     (contact rule d < cfac*(ri+rj), cfac=1.15 — the LIVEFAB rule);
 *   - channels have IDENTITY: born free at (e_mid=0, phase=0), they die
 *     ONLY at e_mid=0 after leaving the candidate set (LIVEFAB rule
 *     alpha: flush-on-death — a dying channel finishes its in-flight
 *     cycles and delivers; unreceivable residue returns to its source,
 *     counted);
 *   - the LAW is laws_V2g VERBATIM. The kernel's law never reads
 *     positions (cellfab.c:664 "scaffold positions — init + diagnostics
 *     only"); it reads d, u_hat, and the live lens area A. Here those
 *     come from live geometry instead of a frozen scaffold. No law
 *     constant is changed; defaults below ARE laws_V2g.cfg (2026-07-28).
 *
 * THE TWO GEOMETRIC FORCES (both already in the programme)
 * --------------------------------------------------------
 *  repulsion  pure contact repulsion below d = ri+rj, relaxed to
 *             jamming (LIVEFAB: "a contact rule with energy-dependent
 *             radii, relaxed by pure repulsion to jamming — not a
 *             spring over a cutoff"; the spring variant densified
 *             degree 8.6 -> 47 and is rejected).
 *  the bond   P15 retardation plasticity (cellfab.c:3068-3103), the S2
 *             odd term on the link's geometric conjugate:
 *                 dd = -kappa * base * Sm * dV/dd,
 *                 V = 1 - G(psi_f) G(psi_b),
 *             applied here to the endpoint POSITIONS (+-u_hat*dd/2)
 *             instead of a frozen ld[] parameter. Its lock condition is
 *             the measured pair-separation ladder
 *                 (q w_i + p w_j) d / C = 2 pi m,
 *             so tuned voices acquire quantized bond lengths, while the
 *             vacuum is EXACTLY inert (Sm = 0 when either side is
 *             silent). kappa_bond = 0 removes the bond, changing
 *             nothing else.
 *
 * EXPERIMENTS (exp=...)
 * ---------------------
 *  FC0 exp=bath     vacuum sanity: jam-settles static, no densification,
 *                   ledger at the FP floor with noise churn (e1-class).
 *  FC1 exp=pair     two tuned voices, no bath: the bond. d must walk to
 *                   d* = 2 pi m C/(q w_i + p w_j) from BOTH sides and
 *                   stay. Controls: kappa_bond=0; unloaded pair (force
 *                   exactly 0); detuned pair.
 *  FC2 exp=ring/oct rung-locked truss objects (ring6 floppy vs
 *                   octahedron isostatic): shape held? shear response?
 *                   energy retained?
 *  FC3 exp=blob     e3a-class dense blob ON THE FREE SUBSTRATE at
 *                   density, with the packing/degree columns printed
 *                   beside every row (the FREECELL 5a test), plus the
 *                   dilute discriminator via bath_frac. kx>0 = the e3b
 *                   tilt (phase ramp + prealign), RESULT blob_drift.
 *  FC4              footprint/degree under s_k (g1-analog columns in
 *                   the blob rows).
 *  FC5 exp=pulse    Tier-A e2 analog: field packet on the live
 *                   substrate, RESULT front_speed (cellfab convention).
 *
 * RESULTS (2026-08-02): FREECELL.md §9 and FREECELL_LOG.md §6. The bond
 * is real (two-sided lock at the ladder, offset = pressure/K_b exactly);
 * species select by parity (even rings live, odd π-frustrated dead) and
 * collimation (chains/rings/tubes); yield 5–8% vs derived 7.9%; embedded
 * ring re-locks after bath crush and bleeds at the frozen foam's c0;
 * death = locked shrink along the ladder -> cage-tracking unlock at
 * x~0.19 -> husk; blob holds at density (live >= frozen), fragments
 * dilute; pulse v/C = 0.574; conservation at the FP floor throughout.
 *
 * Predictions pre-registered in FREECELL_LOG.md BEFORE first run.
 * Log convention: stdout -> v89/charge/runs/freecell_<exp>.log
 */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include "../../v89/int_ledger/trig_lut.inc"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#define TWO_PI (2.0 * M_PI)

/* ------------------------------------------------------------------ */
/* config — LAW KEYS default to laws_V2g.cfg (2026-07-28) VERBATIM     */
/* ------------------------------------------------------------------ */

typedef struct {
    unsigned long seed;
    /* --- laws_V2g (do not change defaults: they ARE the table) --- */
    double C, dmin, r0, rjit;
    double w1, w2, q_detune;
    double gamma_res, gamma_res_m, p_gate, lock_floor;
    double k_dep, k_dep_m, cap;
    double e_s0, es_floor, e_cond;
    double f_conv, f_evap, s_pull;
    double kappa_lock, kappa_align, kappa_freq, kappa_reac;
    double s_k, s_disp, sigma_tumble;
    int comb_limit;
    double rough_k, gamma_rough;
    int mob_sym;
    double mob_floor, field_J, quant_A0;
    int quant_mode;
    /* --- freecell geometry sector (apparatus, not law) --- */
    double cfac;          /* candidate rule d < cfac*(ri+rj)  (LIVEFAB) */
    double k_rep;         /* contact repulsion stiffness */
    double mob_geo;       /* overdamped position mobility */
    double kappa_bond;    /* P15 force onto positions (0 = no bond) */
    int    jam_sweeps;    /* init jam settle sweeps */
    double jam_k;         /* settle step (livefab used 0.3) */
    double jam_tol;       /* early-exit witness */
    int    freeze_geo;    /* 1 = positions frozen (frozen-scaffold control) */
    double bath_frac;     /* thin the dart bath (dilute control), 1=full */
    /* --- run + apparatus --- */
    double L, dt, T;
    int diag_every, snap_every;
    char exp[32];
    double noise_amp;
    double amp, sigma;              /* blob / pulse seed */
    double kx;                      /* pulse phase tilt */
    double px, py, pz;
    int    bath;                    /* 1 = dart bath present */
    double pair_x0, pair_x1, pair_doff;  /* pair loads (x1<0: =x0); initial d - d* */
    int    pair_m;                  /* rung index m */
    int    pair_pp, pair_qq;        /* interval p:q */
    int    seedlock;                /* 1 = phases seeded at the lock */
    int    ring_n;                  /* ring voices */
    double ring_x, ring_doff;
    /* FCQ — the fifth-triangle (incomplete harmonic; quark analog).
     * Voices: U at tri_xU, D at the exact fifth load
     * x_D = (0.5 + 1.8 x_U)/1.2 (w_D = (2/3) w_U). tri_kind 0 = UUD
     * (equilateral by the ladder identity d*_UU(m=1) = d*_UD(m=2)),
     * 1 = UDD (isoceles; the D-D edge exceeds contact -> open chain).
     * tri_branch = the Z3 phase branch k seeded on the D voice. */
    int    tri_kind, tri_branch;
    double tri_xU, tri_xD, tri_doff;
    double tri2_sep;                /* COM separation of the second triangle */
    int    tri2_k2;                 /* Z3 branch of the second triangle */
    double oct_x, oct_doff;
    double shear_eps, shear_t;      /* instantaneous deviatoric strain test */
    /* DS — double slit on the free substrate (v90 P1; DS.md). Wall is
     * CARVED VACUUM (no cells = no contact = no transport); 2D bath. */
    double slit_wallx, slit_th, slit_sep, slit_hw;
    double slit_screenx, slit_srcx, slit_sy;
    double slit_t0, slit_t1;        /* screen gate window */
    double slit_pinw;               /* fixture half-width about the wall plane */
    int    slit_mask;               /* 0 both, 1 A only, 2 B only, 3 no wall */
    char snap_dir[256];
    char snap_file[256];            /* FCS v3 single-stream output */
} Cfg;

static Cfg P;

static void cfg_defaults(void)
{
    memset(&P, 0, sizeof(P));
    P.seed = 20260802UL;
    /* laws_V2g.cfg verbatim */
    P.C = 1.0; P.dmin = 1.0; P.r0 = 0.85; P.rjit = 0.06;
    P.w1 = 1.65; P.w2 = 2.9; P.q_detune = 1.2;
    P.gamma_res = 0.25; P.gamma_res_m = 0.10; P.p_gate = 8; P.lock_floor = 0.005;
    P.k_dep = 1.2; P.k_dep_m = 2.0; P.cap = 2.5;
    P.e_s0 = 1.0; P.es_floor = 0.05; P.e_cond = 0;
    P.f_conv = 0.25; P.f_evap = 0.5; P.s_pull = 0.15;
    P.kappa_lock = 1.0; P.kappa_align = 0.5; P.kappa_freq = 0; P.kappa_reac = 0;
    P.s_k = 0.06; P.s_disp = 0.3; P.sigma_tumble = 0.01;
    P.comb_limit = 6; P.rough_k = 0.35; P.gamma_rough = 0.5;
    P.mob_sym = 1; P.mob_floor = 0.004; P.field_J = 1.8;
    P.quant_A0 = 1.15; P.quant_mode = 2;
    /* freecell geometry */
    P.cfac = 1.15; P.k_rep = 1.0; P.mob_geo = 1.0; P.kappa_bond = 1.0;
    P.jam_sweeps = 800; P.jam_k = 0.10; P.jam_tol = 1e-4;
    P.freeze_geo = 0; P.bath_frac = 1.0;
    /* run */
    P.L = 16.0; P.dt = 0.02; P.T = 40.0;
    P.diag_every = 100; P.snap_every = 0;
    strcpy(P.exp, "bath");
    P.noise_amp = 0;
    P.amp = 1.6; P.sigma = 2.2; P.kx = 0;
    P.px = -1; P.py = -1; P.pz = -1;
    P.bath = 1;
    P.pair_x0 = 0.325; P.pair_x1 = -1; P.pair_doff = 0.05; P.pair_m = 1;
    P.pair_pp = 1; P.pair_qq = 1; P.seedlock = 1;
    P.ring_n = 6; P.ring_x = 0.325; P.ring_doff = 0.0;
    P.tri_kind = 0; P.tri_branch = 0; P.tri_xU = 0.20; P.tri_xD = -1; P.tri_doff = 0;
    P.tri2_sep = 2.6; P.tri2_k2 = 0;
    P.oct_x = 0.325; P.oct_doff = 0.0;
    P.shear_eps = 0; P.shear_t = 0;
    P.slit_wallx = 16.0; P.slit_th = 1.6; P.slit_sep = 9.0; P.slit_hw = 2.0;
    P.slit_screenx = 28.0; P.slit_srcx = 8.0; P.slit_sy = 10.0;
    P.slit_t0 = 16.0; P.slit_t1 = 60.0; P.slit_mask = 0;
    P.slit_pinw = 3.0;
    strcpy(P.snap_dir, "");
    strcpy(P.snap_file, "");
}

static void set_kv(const char *k, const char *v)
{
    if (0) {}
    else if (!strcmp(k, "seed")) P.seed = strtoul(v, NULL, 10);
    else if (!strcmp(k, "C")) P.C = atof(v);
    else if (!strcmp(k, "dmin")) P.dmin = atof(v);
    else if (!strcmp(k, "r0")) P.r0 = atof(v);
    else if (!strcmp(k, "rjit")) P.rjit = atof(v);
    else if (!strcmp(k, "w1")) P.w1 = atof(v);
    else if (!strcmp(k, "w2")) P.w2 = atof(v);
    else if (!strcmp(k, "q_detune")) P.q_detune = atof(v);
    else if (!strcmp(k, "gamma_res")) P.gamma_res = atof(v);
    else if (!strcmp(k, "gamma_res_m")) P.gamma_res_m = atof(v);
    else if (!strcmp(k, "p_gate")) P.p_gate = atof(v);
    else if (!strcmp(k, "lock_floor")) P.lock_floor = atof(v);
    else if (!strcmp(k, "k_dep")) P.k_dep = atof(v);
    else if (!strcmp(k, "k_dep_m")) P.k_dep_m = atof(v);
    else if (!strcmp(k, "cap")) P.cap = atof(v);
    else if (!strcmp(k, "e_s0")) P.e_s0 = atof(v);
    else if (!strcmp(k, "es_floor")) P.es_floor = atof(v);
    else if (!strcmp(k, "e_cond")) P.e_cond = atof(v);
    else if (!strcmp(k, "f_conv")) P.f_conv = atof(v);
    else if (!strcmp(k, "f_evap")) P.f_evap = atof(v);
    else if (!strcmp(k, "s_pull")) P.s_pull = atof(v);
    else if (!strcmp(k, "kappa_lock")) P.kappa_lock = atof(v);
    else if (!strcmp(k, "kappa_align")) P.kappa_align = atof(v);
    else if (!strcmp(k, "kappa_freq")) P.kappa_freq = atof(v);
    else if (!strcmp(k, "kappa_reac")) P.kappa_reac = atof(v);
    else if (!strcmp(k, "s_k")) P.s_k = atof(v);
    else if (!strcmp(k, "s_disp")) P.s_disp = atof(v);
    else if (!strcmp(k, "sigma_tumble")) P.sigma_tumble = atof(v);
    else if (!strcmp(k, "comb_limit")) P.comb_limit = atoi(v);
    else if (!strcmp(k, "rough_k")) P.rough_k = atof(v);
    else if (!strcmp(k, "gamma_rough")) P.gamma_rough = atof(v);
    else if (!strcmp(k, "mob_sym")) P.mob_sym = atoi(v);
    else if (!strcmp(k, "mob_floor")) P.mob_floor = atof(v);
    else if (!strcmp(k, "field_J")) P.field_J = atof(v);
    else if (!strcmp(k, "quant_A0")) P.quant_A0 = atof(v);
    else if (!strcmp(k, "quant_mode")) P.quant_mode = atoi(v);
    else if (!strcmp(k, "cfac")) P.cfac = atof(v);
    else if (!strcmp(k, "k_rep")) P.k_rep = atof(v);
    else if (!strcmp(k, "mob_geo")) P.mob_geo = atof(v);
    else if (!strcmp(k, "kappa_bond")) P.kappa_bond = atof(v);
    else if (!strcmp(k, "jam_sweeps")) P.jam_sweeps = atoi(v);
    else if (!strcmp(k, "jam_k")) P.jam_k = atof(v);
    else if (!strcmp(k, "jam_tol")) P.jam_tol = atof(v);
    else if (!strcmp(k, "freeze_geo")) P.freeze_geo = atoi(v);
    else if (!strcmp(k, "bath_frac")) P.bath_frac = atof(v);
    else if (!strcmp(k, "L")) P.L = atof(v);
    else if (!strcmp(k, "dt")) P.dt = atof(v);
    else if (!strcmp(k, "T")) P.T = atof(v);
    else if (!strcmp(k, "diag_every")) P.diag_every = atoi(v);
    else if (!strcmp(k, "snap_every")) P.snap_every = atoi(v);
    else if (!strcmp(k, "exp")) { strncpy(P.exp, v, 31); P.exp[31] = 0; }
    else if (!strcmp(k, "noise_amp")) P.noise_amp = atof(v);
    else if (!strcmp(k, "amp")) P.amp = atof(v);
    else if (!strcmp(k, "sigma")) P.sigma = atof(v);
    else if (!strcmp(k, "kx")) P.kx = atof(v);
    else if (!strcmp(k, "px")) P.px = atof(v);
    else if (!strcmp(k, "py")) P.py = atof(v);
    else if (!strcmp(k, "pz")) P.pz = atof(v);
    else if (!strcmp(k, "bath")) P.bath = atoi(v);
    else if (!strcmp(k, "pair_x0")) P.pair_x0 = atof(v);
    else if (!strcmp(k, "pair_x1")) P.pair_x1 = atof(v);
    else if (!strcmp(k, "pair_doff")) P.pair_doff = atof(v);
    else if (!strcmp(k, "pair_m")) P.pair_m = atoi(v);
    else if (!strcmp(k, "pair_p")) P.pair_pp = atoi(v);
    else if (!strcmp(k, "pair_q")) P.pair_qq = atoi(v);
    else if (!strcmp(k, "seedlock")) P.seedlock = atoi(v);
    else if (!strcmp(k, "ring_n")) P.ring_n = atoi(v);
    else if (!strcmp(k, "tri_kind")) P.tri_kind = atoi(v);
    else if (!strcmp(k, "tri_branch")) P.tri_branch = atoi(v);
    else if (!strcmp(k, "tri_xU")) P.tri_xU = atof(v);
    else if (!strcmp(k, "tri_xD")) P.tri_xD = atof(v);
    else if (!strcmp(k, "tri_doff")) P.tri_doff = atof(v);
    else if (!strcmp(k, "tri2_sep")) P.tri2_sep = atof(v);
    else if (!strcmp(k, "tri2_k2")) P.tri2_k2 = atoi(v);
    else if (!strcmp(k, "ring_x")) P.ring_x = atof(v);
    else if (!strcmp(k, "ring_doff")) P.ring_doff = atof(v);
    else if (!strcmp(k, "oct_x")) P.oct_x = atof(v);
    else if (!strcmp(k, "oct_doff")) P.oct_doff = atof(v);
    else if (!strcmp(k, "shear_eps")) P.shear_eps = atof(v);
    else if (!strcmp(k, "shear_t")) P.shear_t = atof(v);
    else if (!strcmp(k, "slit_wallx")) P.slit_wallx = atof(v);
    else if (!strcmp(k, "slit_th")) P.slit_th = atof(v);
    else if (!strcmp(k, "slit_sep")) P.slit_sep = atof(v);
    else if (!strcmp(k, "slit_hw")) P.slit_hw = atof(v);
    else if (!strcmp(k, "slit_screenx")) P.slit_screenx = atof(v);
    else if (!strcmp(k, "slit_srcx")) P.slit_srcx = atof(v);
    else if (!strcmp(k, "slit_sy")) P.slit_sy = atof(v);
    else if (!strcmp(k, "slit_t0")) P.slit_t0 = atof(v);
    else if (!strcmp(k, "slit_t1")) P.slit_t1 = atof(v);
    else if (!strcmp(k, "slit_mask")) P.slit_mask = atoi(v);
    else if (!strcmp(k, "slit_pinw")) P.slit_pinw = atof(v);
    else if (!strcmp(k, "snap_dir")) { strncpy(P.snap_dir, v, 255); P.snap_dir[255] = 0; }
    else if (!strcmp(k, "snap_file")) { strncpy(P.snap_file, v, 255); P.snap_file[255] = 0; }
    else fprintf(stderr, "# WARN unknown key %s\n", k);
}

static void load_cfg(const char *path)
{
    FILE *f = fopen(path, "r");
    if (!f) { fprintf(stderr, "cannot open %s\n", path); exit(1); }
    char line[512];
    while (fgets(line, sizeof line, f)) {
        if (line[0] == '#' || line[0] == '\n') continue;
        char *eq = strchr(line, '=');
        if (!eq) continue;
        *eq = 0;
        char *k = line, *v = eq + 1;
        char *nl = strchr(v, '\n'); if (nl) *nl = 0;
        char *hs = strchr(v, '#'); if (hs) *hs = 0;
        while (*k == ' ') k++;
        char *ke = k + strlen(k); while (ke > k && ke[-1] == ' ') *--ke = 0;
        while (*v == ' ') v++;
        char *ve = v + strlen(v); while (ve > v && (ve[-1] == ' ' || ve[-1] == '\r')) *--ve = 0;
        if (*k) set_kv(k, v);
    }
    fclose(f);
}

/* ------------------------------------------------------------------ */
/* rng + small math (cellfab-identical xorshift64)                     */
/* ------------------------------------------------------------------ */

static uint64_t rng_s;
static uint64_t xrand(void)
{ uint64_t x = rng_s; x ^= x << 13; x ^= x >> 7; x ^= x << 17; return rng_s = x; }
static double frand(void) { return (double)(xrand() >> 11) * (1.0 / 9007199254740992.0); }
static double grand_(void)
{
    double u1 = frand() + 1e-18, u2 = frand();
    return sqrt(-2.0 * log(u1)) * cos(TWO_PI * u2);
}
static double wrap_pi(double a)
{
    a = fmod(a + M_PI, TWO_PI);
    if (a < 0) a += TWO_PI;
    return a - M_PI;
}
/* periodic box: minimum-image separation and position fold */
static double wr(double d)
{
    if (d >  0.5 * P.L) d -= P.L;
    if (d < -0.5 * P.L) d += P.L;
    return d;
}
static double fold(double x)
{
    x = fmod(x, P.L);
    if (x < 0) x += P.L;
    return x;
}
static double gate_of(double dphi)      /* cellfab.c:638 verbatim */
{
    double g = 0.5 * (1.0 + lut_cos(dphi));
    int ip = (int)P.p_gate;
    if ((double)ip == P.p_gate && ip >= 1 && ip <= 8) {
        double r = 1.0;
        for (int q = 0; q < ip; q++) r *= g;
        return r;
    }
    return pow(g, P.p_gate);
}
static void rand_unit(double *x, double *y, double *z)
{
    double a, b, c, n2;
    do { a = grand_(); b = grand_(); c = grand_(); n2 = a*a + b*b + c*c; } while (n2 < 1e-12);
    double inv = 1.0 / sqrt(n2);
    *x = a * inv; *y = b * inv; *z = c * inv;
}

/* the partial comb — cellfab.c:589 verbatim */
static int ncomb = 0;
static int combp[24], combq[24];
static int gcd_i(int a, int b) { while (b) { int t = a % b; a = b; b = t; } return a; }
static void comb_build(void)
{
    ncomb = 0;
    for (int pp = 1; pp <= 6; pp++)
        for (int qq = 1; qq <= 6; qq++) {
            if (pp * qq > P.comb_limit) continue;
            if (gcd_i(pp, qq) != 1) continue;
            combp[ncomb] = pp; combq[ncomb] = qq; ncomb++;
        }
    for (int a = 0; a < ncomb; a++)
        for (int b = a + 1; b < ncomb; b++)
            if (combp[b] * combq[b] < combp[a] * combq[a]) {
                int t;
                t = combp[a]; combp[a] = combp[b]; combp[b] = t;
                t = combq[a]; combq[a] = combq[b]; combq[b] = t;
            }
}

/* ------------------------------------------------------------------ */
/* state                                                               */
/* ------------------------------------------------------------------ */

static int NC = 0;
/* cells */
static double *px_, *py_, *pz_;      /* POSITIONS — dynamical here      */
static double *cr0, *cr;
static double *n1x, *n1y, *n1z, *n2x, *n2y, *n2z;
static double *th2, *cbeta;
static double *w1e, *w2e;
static double *Es, *Em, *Ee;
static double *fa1, *fa2;
static double *flload;
static double *qcnvD, *qcnvF;
static double *roughq;
static double *req1, *scl1;
static double *sprq, *sscl;
static double *fsum_;
static unsigned char *tag;
static unsigned char *pin;   /* apparatus fixture: pass D skips pinned cells */
static double *fxb, *fyb, *fzb;      /* geometric force gather buffers  */
static double *rngbuf, *nsnap, *th2s;

/* channels — PERSISTENT SLOTS with identity (the birth/death ledger).
 * A slot is FREE, ALIVE (in the candidate set), or DYING (out of the
 * candidate set, flushing in-flight energy under rule alpha). Slots
 * carry the DENSE sector only: field gated transport is retired in the
 * standing kernel (cellfab.c:2952-2954), the field moves in pass F. */
#define S_FREE  0
#define S_ALIVE 1
#define S_DYING 2
static int NLMAX = 0, NSLOT = 0;     /* high-water slot count           */
static int *sli, *slj;
static unsigned char *sst;
static double *sd, *sux, *suy, *suz; /* live d, u_hat (i -> j)          */
static double *sA;                   /* live lens area                  */
static double *slem, *slph;          /* [slot][dir]: in-flight, phase   */
static signed char *slp, *slq;       /* locked-in partial ratio p:q     */
static double *swant;                /* [slot][dir] wants               */
static double *sflux;                /* dense deposits this step        */
static double *sldd;                 /* bond misfit buffer (Jacobi)     */
static double *swl;                  /* space flux per slot             */
static double *sfluxd;               /* cumulative DIRECTED dense flux [slot][dir] */
static int *freelist, nfree = 0;
static long births = 0, deaths = 0, beta_returns = 0;
static double beta_energy = 0;

/* pair -> slot hash (open addressing, tombstones) */
static int64_t *hkey; static int *hval; static int hsize = 0;
#define H_EMPTY (-1)
#define H_TOMB  (-2)
static int64_t pkey(int i, int j) { return ((int64_t)i << 20) | j; }
static int hfind(int i, int j)
{
    int64_t k = pkey(i, j);
    int h = (int)((k * 0x9E3779B97F4A7C15ULL) >> 40) & (hsize - 1);
    for (;;) {
        int64_t kk = hkey[h];
        if (kk == k) return hval[h];
        if (kk == H_EMPTY) return -1;
        h = (h + 1) & (hsize - 1);
    }
}
static void hput(int i, int j, int slot)
{
    int64_t k = pkey(i, j);
    int h = (int)((k * 0x9E3779B97F4A7C15ULL) >> 40) & (hsize - 1);
    for (;;) {
        int64_t kk = hkey[h];
        if (kk == H_EMPTY || kk == H_TOMB || kk == k) { hkey[h] = k; hval[h] = slot; return; }
        h = (h + 1) & (hsize - 1);
    }
}
static void hdel(int i, int j)
{
    int64_t k = pkey(i, j);
    int h = (int)((k * 0x9E3779B97F4A7C15ULL) >> 40) & (hsize - 1);
    for (;;) {
        int64_t kk = hkey[h];
        if (kk == k) { hkey[h] = H_TOMB; return; }
        if (kk == H_EMPTY) return;
        h = (h + 1) & (hsize - 1);
    }
}

/* cell -> incident slots (CSR, canonical: sorted by neighbour id) */
static int *cls_, *clidx;

/* spatial hash for candidate discovery */
static int *bin_head, *bin_next, bin_g = 0;
static double bin_sz = 0;

/* conversion ledgers */
static double rough_total = 0, cond_total = 0, evap_total = 0, backs_total = 0;
static double A0eff = 0;
static long qfire_n = 0;
static double sim_t = 0;
static double E0_total = 0;
static double cenx, ceny, cenz;
static long ncand_last = 0;
static double fs_t[512], fs_r[512], fs_y[512], fs_z[512];
static int nfsamp = 0;
#define DS_NBIN 96
static double ds_I[DS_NBIN];
static double ds_expo = 0;
static double *ds_cellI = NULL;   /* per-cell time-integrated Ee on the strip */

/* ------------------------------------------------------------------ */
/* atoms at mode boundaries — cellfab.c:2739 verbatim                  */
/* ------------------------------------------------------------------ */

static double atoms_eps(double eps_src, double eps_dst)
{ return P.quant_mode >= 3 ? eps_dst : eps_src; }
static double atoms_w(double w_src, double w_dst)
{ return P.quant_mode >= 3 ? w_dst : w_src; }
static double atoms_fire(double demand, double eps_src, double eps_dst, double *cred)
{
    if (A0eff <= 0 || demand <= 0) return demand;
    double eps = atoms_eps(eps_src, eps_dst);
    if (eps <= 0) return demand;
    if (P.quant_mode == 2 || P.quant_mode == 4) {
        *cred += demand;
        if (*cred > 2.0 * eps) *cred = 2.0 * eps;
        double mv = floor(*cred / eps) * eps;
        *cred -= mv;
        return mv;
    }
    return floor(demand / eps) * eps;
}
static double atoms_clamp(double mv, double ceil_e, double eps_src, double eps_dst,
                          double *cred)
{
    if (mv <= ceil_e) return mv;
    if (A0eff <= 0) return ceil_e;
    double eps = atoms_eps(eps_src, eps_dst);
    double keep = eps > 0 ? floor(ceil_e / eps) * eps : ceil_e;
    if ((P.quant_mode == 2 || P.quant_mode == 4) && cred) *cred += mv - keep;
    return keep;
}
static void qatom_diag(int fd, double w, double e)
{
    if (A0eff <= 0 || e <= 0) return;
    if ((qfire_n++ % 200) == 0)
        printf("# QATOM t=%.2f dir=%s w=%.9g e=%.12g\n",
               sim_t, fd ? "FD" : "DF", w, e);
}

static void field_inject(int i, double dE)   /* cellfab.c:2711 verbatim */
{
    if (dE <= 0) return;
    double e = fa1[i]*fa1[i] + fa2[i]*fa2[i];
    if (e > 1e-20) {
        double fac = sqrt((e + dE) / e);
        fa1[i] *= fac; fa2[i] *= fac;
    } else {
        fa1[i] = sqrt(dE) * lut_cos(th2[i]);
        fa2[i] = sqrt(dE) * lut_sin(th2[i]);
    }
    Ee[i] = e + dE;
}

/* ------------------------------------------------------------------ */
/* allocation                                                          */
/* ------------------------------------------------------------------ */

static void alloc_all(int nc)
{
    NC = nc;
    NLMAX = 80 * nc + 1024;
    px_ = malloc(nc * sizeof(double)); py_ = malloc(nc * sizeof(double)); pz_ = malloc(nc * sizeof(double));
    cr0 = malloc(nc * sizeof(double)); cr = malloc(nc * sizeof(double));
    n1x = malloc(nc * sizeof(double)); n1y = malloc(nc * sizeof(double)); n1z = malloc(nc * sizeof(double));
    n2x = malloc(nc * sizeof(double)); n2y = malloc(nc * sizeof(double)); n2z = malloc(nc * sizeof(double));
    th2 = malloc(nc * sizeof(double)); cbeta = malloc(nc * sizeof(double));
    w1e = malloc(nc * sizeof(double)); w2e = malloc(nc * sizeof(double));
    Es = malloc(nc * sizeof(double)); Em = calloc(nc, sizeof(double)); Ee = calloc(nc, sizeof(double));
    fa1 = calloc(nc, sizeof(double)); fa2 = calloc(nc, sizeof(double));
    flload = calloc(nc, sizeof(double));
    qcnvD = calloc(nc, sizeof(double)); qcnvF = calloc(nc, sizeof(double));
    roughq = calloc(nc, sizeof(double));
    req1 = malloc(nc * sizeof(double)); scl1 = malloc(nc * sizeof(double));
    sprq = malloc(nc * sizeof(double)); sscl = malloc(nc * sizeof(double));
    fsum_ = malloc(nc * sizeof(double));
    tag = calloc(nc, 1);
    pin = calloc(nc, 1);
    fxb = malloc(nc * sizeof(double)); fyb = malloc(nc * sizeof(double)); fzb = malloc(nc * sizeof(double));
    rngbuf = malloc(6 * (size_t)nc * sizeof(double));
    nsnap = malloc(6 * (size_t)nc * sizeof(double));
    th2s = malloc(nc * sizeof(double));
    cls_ = malloc((nc + 1) * sizeof(int));
    clidx = malloc((size_t)2 * NLMAX * sizeof(int));

    sli = malloc(NLMAX * sizeof(int)); slj = malloc(NLMAX * sizeof(int));
    sst = calloc(NLMAX, 1);
    sd = malloc(NLMAX * sizeof(double));
    sux = malloc(NLMAX * sizeof(double)); suy = malloc(NLMAX * sizeof(double)); suz = malloc(NLMAX * sizeof(double));
    sA = calloc(NLMAX, sizeof(double));
    slem = calloc((size_t)2 * NLMAX, sizeof(double));
    slph = calloc((size_t)2 * NLMAX, sizeof(double));
    slp = malloc(NLMAX); slq = malloc(NLMAX);
    swant = calloc((size_t)2 * NLMAX, sizeof(double));
    sflux = calloc(NLMAX, sizeof(double));
    sldd = calloc(NLMAX, sizeof(double));
    swl = calloc(NLMAX, sizeof(double));
    sfluxd = calloc((size_t)2 * NLMAX, sizeof(double));
    freelist = malloc(NLMAX * sizeof(int));
    nfree = 0;
    for (int s = NLMAX - 1; s >= 0; s--) freelist[nfree++] = s;
    NSLOT = 0;

    hsize = 1; while (hsize < 4 * NLMAX) hsize <<= 1;
    hkey = malloc((size_t)hsize * sizeof(int64_t));
    hval = malloc((size_t)hsize * sizeof(int));
    for (int h = 0; h < hsize; h++) hkey[h] = H_EMPTY;

    int nb = 1;
    bin_head = NULL; bin_next = malloc(nc * sizeof(int));
    (void)nb;
}

/* ------------------------------------------------------------------ */
/* geometry: spatial hash + candidate refresh + birth/death (rule α)   */
/* ------------------------------------------------------------------ */

static void bins_build(void)
{
    double rmax = 0;
    for (int i = 0; i < NC; i++) if (cr[i] > rmax) rmax = cr[i];
    double cut = P.cfac * 2.0 * rmax;
    if (cut < 1e-6) cut = 1e-6;
    int g = (int)(P.L / cut); if (g < 1) g = 1; if (g > 128) g = 128;
    if (g != bin_g) {
        free(bin_head);
        bin_head = malloc((size_t)g * g * g * sizeof(int));
        bin_g = g;
    }
    bin_sz = P.L / g;
    for (int b = 0; b < g * g * g; b++) bin_head[b] = -1;
    for (int i = 0; i < NC; i++) {
        int bx = (int)(px_[i] / bin_sz), by = (int)(py_[i] / bin_sz), bz = (int)(pz_[i] / bin_sz);
        if (bx >= g) bx = g - 1; if (by >= g) by = g - 1; if (bz >= g) bz = g - 1;
        if (bx < 0) bx = 0; if (by < 0) by = 0; if (bz < 0) bz = 0;
        int b = (bx * g + by) * g + bz;
        bin_next[i] = bin_head[b]; bin_head[b] = i;
    }
}

static int slot_new(int i, int j)
{
    if (nfree == 0) {
        fprintf(stderr, "FATAL: slot table full at pair (%d,%d) NSLOT=%d NC=%d "
                "pos_i=(%.2f,%.2f,%.2f) pos_j=(%.2f,%.2f,%.2f)\n",
                i, j, NSLOT, NC, px_[i], py_[i], pz_[i], px_[j], py_[j], pz_[j]);
        exit(2);
    }
    int s = freelist[--nfree];
    sli[s] = i; slj[s] = j; sst[s] = S_ALIVE;
    slem[2*s] = slem[2*s+1] = 0; slph[2*s] = slph[2*s+1] = 0;
    slp[s] = 1; slq[s] = 1; sA[s] = 0; sldd[s] = 0; swl[s] = 0;
    swant[2*s] = swant[2*s+1] = 0; sflux[s] = 0;
    hput(i, j, s);
    if (s >= NSLOT) NSLOT = s + 1;
    births++;
    return s;
}

/* candidate refresh: births, dying transitions, rule-α frees, live d/u.
 * Canonical order everywhere (cells ascending, neighbours sorted). */
static int nbuf[512];
static void topo_refresh(void)
{
    bins_build();
    long ncand = 0;
    int g = bin_g;
    /* mark: any ALIVE slot no longer in the candidate set -> DYING.
     * First pass: flag all non-free slots, then re-affirm candidates. */
    for (int s = 0; s < NSLOT; s++)
        if (sst[s] == S_ALIVE) sst[s] = S_DYING;   /* provisional */
    for (int i = 0; i < NC; i++) {
        int bx = (int)(px_[i] / bin_sz), by = (int)(py_[i] / bin_sz), bz = (int)(pz_[i] / bin_sz);
        if (bx >= g) bx = g - 1; if (by >= g) by = g - 1; if (bz >= g) bz = g - 1;
        if (bx < 0) bx = 0; if (by < 0) by = 0; if (bz < 0) bz = 0;
        int nn = 0;
        for (int ax = bx - 1; ax <= bx + 1; ax++)
            for (int ay = by - 1; ay <= by + 1; ay++)
                for (int az = bz - 1; az <= bz + 1; az++) {
                    int wx = (ax + g) % g, wy2 = (ay + g) % g, wz = (az + g) % g;
                    for (int q = bin_head[(wx * g + wy2) * g + wz]; q >= 0; q = bin_next[q]) {
                        if (q <= i) continue;
                        double dx = wr(px_[q] - px_[i]), dy = wr(py_[q] - py_[i]), dz = wr(pz_[q] - pz_[i]);
                        double d2 = dx*dx + dy*dy + dz*dz;
                        double cut = P.cfac * (cr[i] + cr[q]);
                        if (d2 >= cut * cut || d2 <= 1e-12) continue;
                        int dup = 0;
                        for (int a2 = 0; a2 < nn; a2++) if (nbuf[a2] == q) dup = 1;
                        if (!dup && nn < 512) nbuf[nn++] = q;
                    }
                }
        /* canonical: sort the neighbour list ascending */
        for (int a = 1; a < nn; a++) {
            int v = nbuf[a], b = a - 1;
            while (b >= 0 && nbuf[b] > v) { nbuf[b + 1] = nbuf[b]; b--; }
            nbuf[b + 1] = v;
        }
        for (int a = 0; a < nn; a++) {
            int j = nbuf[a];
            ncand++;
            int s = hfind(i, j);
            if (s < 0) s = slot_new(i, j);
            else sst[s] = S_ALIVE;                 /* re-affirm / resurrect */
        }
    }
    ncand_last = ncand;
    /* rule α: a DYING slot with no in-flight energy is freed */
    for (int s = 0; s < NSLOT; s++) {
        if (sst[s] == S_DYING && slem[2*s] == 0 && slem[2*s+1] == 0) {
            hdel(sli[s], slj[s]);
            sst[s] = S_FREE;
            freelist[nfree++] = s;
            deaths++;
        }
    }
    /* live geometry for every non-free slot (+ lens area for the t=0
     * diagnostics; pass 2 recomputes it identically each step) */
    for (int s = 0; s < NSLOT; s++) {
        if (sst[s] == S_FREE) continue;
        int i = sli[s], j = slj[s];
        double dx = wr(px_[j] - px_[i]), dy = wr(py_[j] - py_[i]), dz = wr(pz_[j] - pz_[i]);
        double d = sqrt(dx*dx + dy*dy + dz*dz);
        if (d < 1e-9) d = 1e-9;
        sd[s] = d; sux[s] = dx / d; suy[s] = dy / d; suz[s] = dz / d;
        double ri = cr[i], rj = cr[j], A = 0;
        if (d < ri + rj) {
            double t = d * d - rj * rj + ri * ri;
            double a2 = (4.0 * d * d * ri * ri - t * t) / (4.0 * d * d);
            if (a2 > 0) A = M_PI * a2;
            else if (d < fabs(ri - rj)) { double rm = ri < rj ? ri : rj; A = M_PI * rm * rm; }
        }
        sA[s] = A;
    }
    /* CSR, canonical (slots visited in index order after per-cell sort) */
    for (int i = 0; i <= NC; i++) cls_[i] = 0;
    for (int s = 0; s < NSLOT; s++) {
        if (sst[s] == S_FREE) continue;
        cls_[sli[s] + 1]++; cls_[slj[s] + 1]++;
    }
    for (int i = 0; i < NC; i++) cls_[i + 1] += cls_[i];
    int *fill = malloc(NC * sizeof(int));
    memset(fill, 0, NC * sizeof(int));
    for (int s = 0; s < NSLOT; s++) {
        if (sst[s] == S_FREE) continue;
        clidx[cls_[sli[s]] + fill[sli[s]]++] = s;
        clidx[cls_[slj[s]] + fill[slj[s]]++] = s;
    }
    free(fill);
    /* sort each cell's incident list by (other endpoint) for canonical order */
    for (int i = 0; i < NC; i++) {
        int lo = cls_[i], hi = cls_[i + 1];
        for (int a = lo + 1; a < hi; a++) {
            int v = clidx[a];
            int vk = (sli[v] == i) ? slj[v] : sli[v];
            int b = a - 1;
            while (b >= lo) {
                int u = clidx[b];
                int uk = (sli[u] == i) ? slj[u] : sli[u];
                if (uk <= vk) break;
                clidx[b + 1] = u; b--;
            }
            clidx[b + 1] = v;
        }
    }
}

/* ------------------------------------------------------------------ */
/* the step — passes ported from cellfab.c step_field() in order       */
/* ------------------------------------------------------------------ */

static double gm_, G2m_, lockf_, Aref_, dref_;

static void step(void)
{
    double dt = P.dt;

    /* pass 0: pitch share of in-flight dense energy (cellfab.c:2795) */
    for (int i = 0; i < NC; i++) {
        double fl = 0;
        for (int q = cls_[i]; q < cls_[i + 1]; q++) {
            int s = clidx[q];
            fl += slem[2*s] + slem[2*s+1];
        }
        flload[i] = 0.5 * fl;
    }

    /* pass S: space-mode transport (cellfab.c:2811) — stale-by-one lA,
     * exactly as the kernel (its pass S precedes the lA recompute) */
    if (P.s_k > 0) {
        for (int s = 0; s < NSLOT; s++) {
            swl[s] = 0.0;
            if (sst[s] == S_FREE || sA[s] <= 0) continue;
            int i = sli[s], j = slj[s];
            double pi_ = Es[i] + P.s_disp * (Em[i] + Ee[i]);
            double pj_ = Es[j] + P.s_disp * (Em[j] + Ee[j]);
            double dp = pi_ - pj_;
            if (dp == 0) continue;
            double w = (sA[s] / Aref_) * (dref_ / sd[s]);
            swl[s] = P.s_k * dt * w * dp;
        }
        for (int i = 0; i < NC; i++) {
            double rq = 0;
            for (int q = cls_[i]; q < cls_[i + 1]; q++) {
                int s = clidx[q];
                double f = swl[s];
                if ((f > 0 && sli[s] == i) || (f < 0 && slj[s] == i)) rq += fabs(f);
            }
            sprq[i] = rq;
            double avail = Es[i] - P.es_floor;
            sscl[i] = (avail <= 0) ? 0.0 : (rq > avail ? avail / rq : 1.0);
        }
        for (int i = 0; i < NC; i++) {
            double de = 0;
            for (int q = cls_[i]; q < cls_[i + 1]; q++) {
                int s = clidx[q];
                double f = swl[s];
                if (f == 0) continue;
                int src = f > 0 ? sli[s] : slj[s];
                double mag = fabs(f) * sscl[src];
                de += (src == i) ? -mag : mag;
            }
            Es[i] += de;
        }
    }

    /* pass 1: live radii, effective frequencies (cellfab.c:2890) */
    for (int i = 0; i < NC; i++) {
        double ratio = Es[i] / P.e_s0;
        cr[i] = cr0[i] * cbrt(ratio > 0 ? ratio : 0);
        double x = (Em[i] + flload[i]) / P.cap;
        double det = 1.0 + P.q_detune * x;
        w1e[i] = P.w1 / det;
        w2e[i] = P.w2 / det;
    }

    /* G0: topology from live geometry (the free-cell replacement for
     * the init-time dart throw; LIVEFAB §5 item 3) */
    if (!P.freeze_geo) topo_refresh();

    /* pass 2: channel lens area + dense wants + bond misfit buffer
     * (cellfab.c:2911, c==1 branch; field gated transport retired) */
    for (int s = 0; s < NSLOT; s++) {
        if (sst[s] == S_FREE) { sA[s] = 0; continue; }
        int i = sli[s], j = slj[s];
        double d = sd[s], ri = cr[i], rj = cr[j];
        double A = 0;
        if (d < ri + rj) {
            double t = d * d - rj * rj + ri * ri;
            double a2 = (4.0 * d * d * ri * ri - t * t) / (4.0 * d * d);
            if (a2 > 0) A = M_PI * a2;
            else if (d < fabs(ri - rj)) { double rm = ri < rj ? ri : rj; A = M_PI * rm * rm; }
        }
        sA[s] = A;
        swant[2*s] = swant[2*s+1] = 0;
        sflux[s] = 0;
        sldd[s] = 0;
        if (A <= 0) continue;

        double geo = (A / Aref_) * (dref_ / d);
        double occ_i = (Em[i] + Ee[i]) / P.cap;
        double occ_j = (Em[j] + Ee[j]) / P.cap;
        double head_i = 1.0 - occ_i; if (head_i < 0) head_i = 0; if (head_i > 1) head_i = 1;
        double head_j = 1.0 - occ_j; if (head_j < 0) head_j = 0; if (head_j > 1) head_j = 1;

        double d1i = n1x[i]*sux[s] + n1y[i]*suy[s] + n1z[i]*suz[s];
        double d2i = n2x[i]*sux[s] + n2y[i]*suy[s] + n2z[i]*suz[s];
        double d1j = n1x[j]*sux[s] + n1y[j]*suy[s] + n1z[j]*suz[s];
        double d2j = n2x[j]*sux[s] + n2y[j]*suy[s] + n2z[j]*suz[s];
        double axi = (1.0 - d1i*d1i) * (1.0 - d2i*d2i); if (axi < 0) axi = 0;
        double axj = (1.0 - d1j*d1j) * (1.0 - d2j*d2j); if (axj < 0) axj = 0;
        double axij = axi * axj;
        if (axij < 1e-10) continue;

        double mob_i = Em[i], mob_j = Em[j];
        if (mob_i < 1e-15 && mob_j < 1e-15) continue;
        double dnn = n2x[i]*n2x[j] + n2y[i]*n2y[j] + n2z[i]*n2z[j];
        double gpl = axij * dnn * dnn;
        if (gpl < 1e-8) continue;

        double wi = w2e[i], wj = w2e[j];
        double thi = th2[i], thj = th2[j];

        /* C1 sympathetic joining through coincident partials */
        double best = -1.0;
        int bp = 1, bq = 1;
        for (int tt = 0; tt < ncomb; tt++) {
            double pp = combp[tt], qq = combq[tt];
            double gw = gm_ / (pp * qq);
            double det = qq * wi - pp * wj;
            double sc = (gw * gw / (gw * gw + det * det)) / (pp * qq);
            if (sc > best) { best = sc; bp = combp[tt]; bq = combq[tt]; }
        }
        double res = best;
        slp[s] = (signed char)bp;
        slq[s] = (signed char)bq;
        double g_ij = gate_of(bq * thi - bq * wi * d / P.C - bp * thj);
        double g_ji = gate_of(bp * thj - bp * wj * d / P.C - bq * thi);
        double mi_eff = mob_i, mj_eff = mob_j;
        if (P.mob_sym) {
            double fl0 = P.mob_floor * P.cap;
            double mjr = mob_j > fl0 ? mob_j : fl0;
            double mir = mob_i > fl0 ? mob_i : fl0;
            mi_eff = sqrt(mob_i * mjr);
            mj_eff = sqrt(mob_j * mir);
        }
        double kd = P.k_dep * P.k_dep_m;
        double base = kd * dt * geo * gpl * res;
        double w_ij = base * g_ij * head_j * mi_eff;
        double w_ji = base * g_ji * head_i * mj_eff;

        if (P.kappa_reac > 0) {
            /* S2 — the choir's correction, DERIVED (cellfab.c:3042
             * verbatim): the reactive (odd) component of sympathetic
             * exchange — the interference cross-flow the even-gate rate
             * compression discards. kappa_reac = 1 is the unitarity
             * point. Zero in laws_V2g (retired at rate level); carried
             * here for the S2-full acceptance probe on the fifth. */
            double ps_ij = wrap_pi(bq * thi - bq * wi * d / P.C - bp * thj);
            double ps_ji = wrap_pi(bp * thj - bp * wj * d / P.C - bq * thi);
            double Sm2 = sqrt(mi_eff * mj_eff);
            double hh = sqrt(head_i * head_j);
            double reac = P.kappa_reac * 0.5 * base * hh * Sm2;
            w_ij -= reac * g_ij * lut_sin(ps_ij);
            w_ji -= reac * g_ji * lut_sin(ps_ji);
            if (w_ij < 0) w_ij = 0;
            if (w_ji < 0) w_ji = 0;
        }

        /* THE BOND — P15 (cellfab.c:3068) on the geometric conjugate,
         * buffered Jacobi, applied to POSITIONS in pass D */
        if (P.kappa_bond > 0) {
            double ps_f = wrap_pi(bq * thi - bq * wi * d / P.C - bp * thj);
            double ps_b = wrap_pi(bp * thj - bp * wj * d / P.C - bq * thi);
            double Smp = sqrt(mi_eff * mj_eff);
            if (Smp > 0) {
                double hf = 0.5 * (1.0 + lut_cos(ps_f));
                double hb = 0.5 * (1.0 + lut_cos(ps_b));
                double Gf = pow(hf, P.p_gate), Gb = pow(hb, P.p_gate);
                double Gfp = hf > 1e-12 ? -0.5 * P.p_gate * (Gf / hf) * lut_sin(ps_f) : 0.0;
                double Gbp = hb > 1e-12 ? -0.5 * P.p_gate * (Gb / hb) * lut_sin(ps_b) : 0.0;
                double dVdd = (Gfp * Gb * bq * wi + Gf * Gbp * bp * wj) / P.C;
                double dd = -P.kappa_bond * base * Smp * dVdd;
                /* integrator guard: the plastic walk is slow (~1e-4/step
                 * in the anneal probe); cap the per-step kick so a deep-
                 * overlap encounter (geo ~ 1/d) cannot throw cells */
                if (dd > 0.05) dd = 0.05;
                if (dd < -0.05) dd = -0.05;
                sldd[s] = dd;
            }
        }

        if (w_ij > 0) swant[2*s]   = w_ij;
        if (w_ji > 0) swant[2*s+1] = w_ji;
    }

    /* pass 3: outflow limiter (cellfab.c:3109), dense only */
    for (int i = 0; i < NC; i++) {
        double r1 = 0;
        for (int q = cls_[i]; q < cls_[i + 1]; q++) {
            int s = clidx[q];
            int dir = sli[s] == i ? 0 : 1;
            r1 += swant[2*s + dir];
        }
        req1[i] = r1;
        double a1 = 0.98 * Em[i];
        scl1[i] = (r1 > a1 && r1 > 0) ? a1 / r1 : 1.0;
    }

    /* pass 4: resolve deposits, debit sources, entrain receivers
     * (cellfab.c:3126) */
    for (int i = 0; i < NC; i++) th2s[i] = th2[i];
    for (int s = 0; s < NSLOT; s++) {
        if (sst[s] == S_FREE || sA[s] <= 0) continue;
        for (int dir = 0; dir < 2; dir++) {
            double f = swant[2*s + dir];
            if (f <= 0) continue;
            int src = dir == 0 ? sli[s] : slj[s];
            f *= scl1[src];
            swant[2*s + dir] = f;
            if (slem[2*s + dir] <= 0) slph[2*s + dir] = 0;
            slem[2*s + dir] += f;
            sflux[s] += f;
            sfluxd[2*s + dir] += f;   /* directed ledger for circulation */
        }
    }
    for (int i = 0; i < NC; i++) {
        double d1 = 0;
        for (int q = cls_[i]; q < cls_[i + 1]; q++) {
            int s = clidx[q];
            if (sA[s] <= 0) continue;
            int dir = sli[s] == i ? 0 : 1;
            double f = swant[2*s + dir];
            if (f > 0) d1 += f;
        }
        Em[i] -= d1;
    }
    for (int i = 0; i < NC; i++) {
        for (int q = cls_[i]; q < cls_[i + 1]; q++) {
            int s = clidx[q];
            if (sA[s] <= 0) continue;
            int dir = sli[s] == i ? 1 : 0;   /* the direction ARRIVING at i */
            int src = dir == 0 ? sli[s] : slj[s];
            double f = swant[2*s + dir];
            if (f <= 0) continue;
            double mobr = Em[i];
            double wsrc = w2e[src];
            double thsrc = th2s[src];
            double ms = dir == 0 ? slq[s] : slp[s];
            double mr = dir == 0 ? slp[s] : slq[s];
            double err = wrap_pi(ms * (thsrc - wsrc * sd[s] / P.C) - mr * th2[i]) / mr;
            double mix = f / (f + mobr + lockf_);
            th2[i] += P.kappa_lock * mix * err;
        }
    }

    /* pass 5: advance cycles; deliver on completion (cellfab.c:3198).
     * ALIVE pinched channels hold flight (kernel semantics); DYING
     * channels advance and flush (rule α), returning unreceivable
     * residue to their source (counted). */
    for (int i2 = 0; i2 < NC; i2++) th2s[i2] = th2[i2];
    for (int recv = 0; recv < NC; recv++) {
        for (int q = cls_[recv]; q < cls_[recv + 1]; q++) {
            int s = clidx[q];
            int dying = sst[s] == S_DYING;
            if (sA[s] <= 0 && !dying) continue;
            int i = sli[s], j = slj[s];
            int dir = sli[s] == recv ? 1 : 0;
            int send = dir == 0 ? i : j;
            double adv = dt * P.C / sd[s];
            int sslot = 2*s + dir;
            if (slem[sslot] <= 0) continue;
            slph[sslot] += adv;
            if (slph[sslot] < 1.0) continue;
            double freec = P.cap - (Em[recv] + Ee[recv]);
            double take = slem[sslot];
            if (take > freec) take = freec > 0 ? freec : 0;
            if (take > 0) {
                double mobprev = Em[recv];
                slem[sslot] -= take;
                /* C2: dissonance radiates (cellfab.c:3228) */
                double det = slq[s] * w2e[i] - slp[s] * w2e[j];
                double R = 2.0 * fabs(det) * P.gamma_rough
                           / (det * det + P.gamma_rough * P.gamma_rough);
                double rough = take * P.rough_k * R;
                double reF = A0eff * w2e[recv] / TWO_PI;
                double reD = A0eff * w1e[recv] / TWO_PI;
                rough = atoms_fire(rough, reF, reD, &qcnvF[recv]);
                rough = atoms_clamp(rough, take, reF, reD, &qcnvF[recv]);
                if (rough > 0) roughq[recv] += rough;
                double back_s = rough * P.s_pull / (1.0 + P.s_pull);
                backs_total += back_s;
                Em[recv] += take - rough;
                field_inject(recv, rough - back_s);
                Es[recv] += back_s;
                rough_total += rough;
                double wsend = w2e[send];
                double thsend = th2s[send];
                double ms = dir == 0 ? slq[s] : slp[s];
                double mr = dir == 0 ? slp[s] : slq[s];
                double err = wrap_pi(ms * (thsend - wsend * sd[s] / P.C) - mr * th2[recv]) / mr;
                double mix = take / (take + mobprev + lockf_);
                th2[recv] += P.kappa_lock * mix * err;
            }
            if (dying && slem[sslot] > 1e-17) {
                /* rule α flush: the unreceivable residue returns to its
                 * source (the β fallback, applied only to dying slots) */
                Em[send] += slem[sslot];
                beta_energy += slem[sslot];
                beta_returns++;
                slem[sslot] = 0; slph[sslot] = 0;
            }
            if (slem[sslot] <= 1e-17) { slem[sslot] = 0; slph[sslot] = 0; }
            else if (take <= 0) slph[sslot] = 0;
            else slph[sslot] -= 1.0;
        }
    }

    /* pass D: MOTION — the free-cell core. Bond displacements (Jacobi
     * buffered, ±u_hat/2) + contact repulsion, overdamped, clamped
     * walls (cellfab anneal convention). Geometry is not a ledger
     * (kernel precedent: P15 moves ld with no energy entry). */
    if (!P.freeze_geo) {
        for (int i = 0; i < NC; i++) { fxb[i] = fyb[i] = fzb[i] = 0; }
        for (int i = 0; i < NC; i++) {
            double dxx = 0, dyy = 0, dzz = 0;
            if (pin[i]) { fxb[i] = fyb[i] = fzb[i] = 0; continue; }
            for (int q = cls_[i]; q < cls_[i + 1]; q++) {
                int s = clidx[q];
                double sgn = (sli[s] == i) ? -1.0 : 1.0;
                /* bond: d changes by sldd -> each endpoint ±sldd/2 */
                if (sldd[s] != 0) {
                    double h = 0.5 * sldd[s] * sgn;
                    dxx += h * sux[s]; dyy += h * suy[s]; dzz += h * suz[s];
                }
                /* repulsion below contact */
                double contact = cr[sli[s]] + cr[slj[s]];
                if (sd[s] < contact) {
                    double f = P.k_rep * (contact - sd[s]) * sgn;
                    dxx += P.mob_geo * P.dt * f * sux[s];
                    dyy += P.mob_geo * P.dt * f * suy[s];
                    dzz += P.mob_geo * P.dt * f * suz[s];
                }
            }
            /* per-step displacement cap (quasi-static motion guard) */
            double dm = sqrt(dxx*dxx + dyy*dyy + dzz*dzz);
            if (dm > 0.1) { double sc = 0.1 / dm; dxx *= sc; dyy *= sc; dzz *= sc; }
            fxb[i] = dxx; fyb[i] = dyy; fzb[i] = dzz;
        }
        for (int i = 0; i < NC; i++) {
            px_[i] = fold(px_[i] + fxb[i]);
            py_[i] = fold(py_[i] + fyb[i]);
            pz_[i] = fold(pz_[i] + fzb[i]);
            if (!(px_[i] == px_[i]) || !(py_[i] == py_[i]) || !(pz_[i] == pz_[i])) {
                fprintf(stderr, "FATAL NaN position cell %d at t=%.3f: f=(%g,%g,%g) "
                        "Em=%g Es=%g cr=%g\n", i, sim_t, fxb[i], fyb[i], fzb[i],
                        Em[i], Es[i], cr[i]);
                for (int q = cls_[i]; q < cls_[i + 1]; q++) {
                    int s = clidx[q];
                    fprintf(stderr, "  slot %d (%d-%d) st=%d d=%g A=%g ldd=%g lem=%g/%g\n",
                            s, sli[s], slj[s], sst[s], sd[s], sA[s], sldd[s],
                            slem[2*s], slem[2*s+1]);
                }
                exit(4);
            }
        }
    }

    /* pass F: unitary field hops (cellfab.c:3289), canonical link order */
    for (int i = 0; i < NC; i++) {
        double ang = w1e[i] * dt;
        double cc, ss;
        sincos(ang, &ss, &cc);
        double a1 = fa1[i], a2 = fa2[i];
        fa1[i] = cc * a1 + ss * a2;
        fa2[i] = -ss * a1 + cc * a2;
    }
    for (int i = 0; i < NC; i++) {
        double sacc = 0;
        for (int q = cls_[i]; q < cls_[i + 1]; q++) {
            int s = clidx[q];
            if (sA[s] <= 0) continue;
            sacc += (sA[s] / Aref_) * (dref_ / sd[s]);
        }
        fsum_[i] = sacc;
    }
    for (int i = 0; i < NC; i++) {
        for (int q = cls_[i]; q < cls_[i + 1]; q++) {
            int s = clidx[q];
            if (sli[s] != i) continue;       /* canonical: apply from the i side */
            if (sA[s] <= 0) continue;
            int j = slj[s];
            double si = fsum_[i], sj = fsum_[j];
            if (si <= 1e-12 || sj <= 1e-12) continue;
            double w = (sA[s] / Aref_) * (dref_ / sd[s]);
            double tau = P.field_J * w / sqrt(si * sj) * dt;
            double cc, ss;
            sincos(tau, &ss, &cc);
            double a1i = fa1[i], a2i = fa2[i], a1j = fa1[j], a2j = fa2[j];
            fa1[i] = cc * a1i + ss * a2j;
            fa2[i] = cc * a2i - ss * a1j;
            fa1[j] = cc * a1j + ss * a2i;
            fa2[j] = cc * a2j - ss * a1i;
        }
    }
    for (int i = 0; i < NC; i++)
        Ee[i] = fa1[i]*fa1[i] + fa2[i]*fa2[i];

    /* pass 6: dense clock + beat-gated conversion (cellfab.c:3448),
     * serial by design */
    for (int i = 0; i < NC; i++) {
        if (roughq[i] > 0) {
            qatom_diag(0, atoms_w(w2e[i], w1e[i]), roughq[i]);
            roughq[i] = 0;
        }
    }
    for (int i = 0; i < NC; i++) {
        th2[i] = fmod(th2[i] + w2e[i] * dt, TWO_PI);
        cbeta[i] += (w1e[i] - w2e[i]) * dt;
        int beat_fire = 0;
        if (cbeta[i] >= TWO_PI) { cbeta[i] -= TWO_PI; beat_fire = 1; }
        else if (cbeta[i] <= -TWO_PI) { cbeta[i] += TWO_PI; beat_fire = 1; }
        if (beat_fire) {
            if (Ee[i] > P.e_cond) {
                double d1 = P.f_conv * (Ee[i] - P.e_cond);
                double eF = A0eff * w1e[i] / TWO_PI;
                double eD = A0eff * w2e[i] / TWO_PI;
                d1 = atoms_fire(d1, eF, eD, &qcnvD[i]);
                d1 = atoms_clamp(d1, 0.98 * Ee[i], eF, eD, &qcnvD[i]);
                if (d1 > 0) {
                    cond_total += d1;
                    qatom_diag(1, atoms_w(w1e[i], w2e[i]), d1);
                    double dsp = P.s_pull * d1;
                    double avail = Es[i] - P.es_floor;
                    if (avail < 0) avail = 0;
                    if (dsp > avail) dsp = avail;
                    double fac = Ee[i] > 0 ? sqrt((Ee[i] - d1) / Ee[i]) : 0;
                    fa1[i] *= fac; fa2[i] *= fac;
                    Ee[i] -= d1;
                    Es[i] -= dsp;
                    Em[i] += d1 + dsp;
                }
            }
            double tot = Em[i] + Ee[i];
            if (tot > P.cap) {
                double d2 = P.f_evap * (tot - P.cap);
                double eF2 = A0eff * w2e[i] / TWO_PI;
                double eD2 = A0eff * w1e[i] / TWO_PI;
                d2 = atoms_fire(d2, eF2, eD2, &qcnvF[i]);
                d2 = atoms_clamp(d2, Em[i], eF2, eD2, &qcnvF[i]);
                if (d2 > 0) {
                    evap_total += d2;
                    qatom_diag(0, atoms_w(w2e[i], w1e[i]), d2);
                    double bs = d2 * P.s_pull / (1.0 + P.s_pull);
                    backs_total += bs;
                    Em[i] -= d2;
                    field_inject(i, d2 - bs);
                    Es[i] += bs;
                }
            }
        }
    }

    /* pass 7: plane re-alignment + tumble (cellfab.c:3515) */
    double sq = P.sigma_tumble * sqrt(dt);
    if (sq > 0)
        for (int i = 0; i < NC; i++)
            for (int k = 0; k < 6; k++) rngbuf[6 * (size_t)i + k] = grand_();
    for (int i = 0; i < NC; i++) {
        nsnap[6*(size_t)i+0] = n1x[i]; nsnap[6*(size_t)i+1] = n1y[i]; nsnap[6*(size_t)i+2] = n1z[i];
        nsnap[6*(size_t)i+3] = n2x[i]; nsnap[6*(size_t)i+4] = n2y[i]; nsnap[6*(size_t)i+5] = n2z[i];
    }
    for (int i = 0; i < NC; i++) {
        for (int c = 0; c < 2; c++) {
            double mx = c == 0 ? n1x[i] : n2x[i];
            double my = c == 0 ? n1y[i] : n2y[i];
            double mz = c == 0 ? n1z[i] : n2z[i];
            double ax = 0, ay = 0, az = 0, fs = 0;
            if (c == 1) {
                for (int q = cls_[i]; q < cls_[i + 1]; q++) {
                    int s = clidx[q];
                    double fl = sflux[s];
                    if (fl <= 0) continue;
                    int o = sli[s] == i ? slj[s] : sli[s];
                    const double *no = nsnap + 6 * (size_t)o + 3;
                    double sgn = (mx * no[0] + my * no[1] + mz * no[2]) >= 0 ? 1.0 : -1.0;
                    double du = mx * sux[s] + my * suy[s] + mz * suz[s];
                    ax += fl * (sgn * no[0] - du * sux[s]);
                    ay += fl * (sgn * no[1] - du * suy[s]);
                    az += fl * (sgn * no[2] - du * suz[s]);
                    fs += fl;
                }
            }
            double vx = mx, vy = my, vz = mz;
            if (fs > 0) {
                double w = P.kappa_align * dt / fs;
                vx += w * ax; vy += w * ay; vz += w * az;
            }
            if (sq > 0) {
                const double *g = rngbuf + 6 * (size_t)i + 3 * c;
                vx += sq * g[0]; vy += sq * g[1]; vz += sq * g[2];
            }
            double nn = sqrt(vx*vx + vy*vy + vz*vz);
            if (nn > 1e-12) {
                double *nx = c == 0 ? &n1x[i] : &n2x[i];
                double *ny = c == 0 ? &n1y[i] : &n2y[i];
                double *nz = c == 0 ? &n1z[i] : &n2z[i];
                *nx = vx / nn; *ny = vy / nn; *nz = vz / nn;
            }
        }
    }
}

/* ------------------------------------------------------------------ */
/* totals + diagnostics                                                */
/* ------------------------------------------------------------------ */

static double total_energy(void)
{
    double s = 0, comp = 0;
    for (int i = 0; i < NC; i++) {
        double v = Es[i] + Em[i] + Ee[i];
        double y = v - comp, t = s + y;
        comp = (t - s) - y; s = t;
    }
    for (int sl = 0; sl < NSLOT; sl++) {
        if (sst[sl] == S_FREE) continue;
        double v = slem[2*sl] + slem[2*sl+1];
        double y = v - comp, t = s + y;
        comp = (t - s) - y; s = t;
    }
    return s;
}

static void geo_stats(double *phi, double *zl, int *nla, int *nld,
                      double *dbar, double *sig_d, double *maxldd)
{
    double vol = 0;
    for (int i = 0; i < NC; i++) vol += (4.0/3.0) * M_PI * cr[i]*cr[i]*cr[i];
    *phi = vol / (P.L * P.L * P.L);
    long na = 0, nd = 0; double sm = 0, sm2 = 0; long nl = 0;
    double ml = 0;
    for (int s = 0; s < NSLOT; s++) {
        if (sst[s] == S_FREE) continue;
        if (sst[s] == S_DYING) { nd++; continue; }
        if (sA[s] > 0) {
            na++;
            sm += sd[s]; sm2 += sd[s]*sd[s]; nl++;
        }
        if (fabs(sldd[s]) > ml) ml = fabs(sldd[s]);
    }
    *zl = 2.0 * na / (NC > 0 ? NC : 1);
    *nla = (int)na; *nld = (int)nd;
    double m = nl ? sm / nl : 0;
    double var = nl ? sm2 / nl - m * m : 0; if (var < 0) var = 0;
    *dbar = m; *sig_d = m > 0 ? sqrt(var) / m : 0;
    *maxldd = ml;
}

static double max_net_force(void)
{
    /* geometric net displacement per step (bond+repulsion buffers) */
    double mx = 0;
    for (int i = 0; i < NC; i++) {
        double f = sqrt(fxb[i]*fxb[i] + fyb[i]*fyb[i] + fzb[i]*fzb[i]);
        if (f > mx) mx = f;
    }
    return mx;
}

/* tagged-object metrics: Em retention, COM, RMS radius, connectivity */
static void tag_stats(double *emt, double *comx, double *comy, double *comz,
                      double *rms, double *conn, double *zblob)
{
    double em = 0, cx2 = 0, cy2 = 0, cz2 = 0, wt = 0;
    int nt = 0;
    for (int i = 0; i < NC; i++) {
        if (!tag[i]) continue;
        nt++;
        /* dense inventory NOW: store + this cell's live half of in-flight
         * (flload[] is a start-of-step snapshot — stale at diag time) */
        double fl = 0;
        for (int q = cls_[i]; q < cls_[i + 1]; q++) {
            int s = clidx[q];
            fl += slem[2*s] + slem[2*s+1];
        }
        double w = Em[i] + 0.5 * fl;
        em += w;
        /* COM in minimum-image coordinates about the seed centre */
        cx2 += w * (cenx + wr(px_[i] - cenx));
        cy2 += w * (ceny + wr(py_[i] - ceny));
        cz2 += w * (cenz + wr(pz_[i] - cenz));
        wt += w;
    }
    *emt = em;
    if (wt > 1e-12) { cx2 /= wt; cy2 /= wt; cz2 /= wt; }
    *comx = cx2; *comy = cy2; *comz = cz2;
    double r2 = 0;
    if (wt > 1e-12)
        for (int i = 0; i < NC; i++) {
            if (!tag[i]) continue;
            double w = Em[i] + flload[i];
            double dx = wr(px_[i]-cx2), dy = wr(py_[i]-cy2), dz = wr(pz_[i]-cz2);
            r2 += w * (dx*dx + dy*dy + dz*dz);
        }
    *rms = wt > 1e-12 ? sqrt(r2 / wt) : 0;
    /* largest connected component among tagged via live (A>0) links */
    static int *comp = NULL, *stack = NULL;
    if (!comp) { comp = malloc(NC * sizeof(int)); stack = malloc(NC * sizeof(int)); }
    for (int i = 0; i < NC; i++) comp[i] = -1;
    int bigc = 0;
    long degsum = 0;
    for (int i = 0; i < NC; i++) {
        if (!tag[i]) continue;
        long dg = 0;
        for (int q = cls_[i]; q < cls_[i + 1]; q++)
            if (sA[clidx[q]] > 0) dg++;
        degsum += dg;
        if (comp[i] >= 0) continue;
        int sz = 0, sp = 0;
        stack[sp++] = i; comp[i] = i;
        while (sp) {
            int u = stack[--sp]; sz++;
            for (int q = cls_[u]; q < cls_[u + 1]; q++) {
                int s = clidx[q];
                if (sA[s] <= 0) continue;
                int v = sli[s] == u ? slj[s] : sli[s];
                if (tag[v] && comp[v] < 0) { comp[v] = i; stack[sp++] = v; }
            }
        }
        if (sz > bigc) bigc = sz;
    }
    *conn = nt ? (double)bigc / nt : 0;
    *zblob = nt ? (double)degsum / nt : 0;
}

static void es_shells(double sh[8])
{
    double cnt[8] = {0};
    for (int k = 0; k < 8; k++) sh[k] = 0;
    double drs = 0.5 * P.L / 8.0;
    for (int i = 0; i < NC; i++) {
        double dx = wr(px_[i]-cenx), dy = wr(py_[i]-ceny), dz = wr(pz_[i]-cenz);
        double r = sqrt(dx*dx + dy*dy + dz*dz);
        int k = (int)(r / drs);
        if (k < 8) { sh[k] += Es[i]; cnt[k] += 1; }
    }
    for (int k = 0; k < 8; k++) if (cnt[k] > 0) sh[k] /= cnt[k];
}

/* ------------------------------------------------------------------ */
/* init                                                                */
/* ------------------------------------------------------------------ */

static int build_bath(double *bx, double *by, double *bz, int nmax)
{
    /* blue-noise dart throw (cellfab.c:799 structure) + optional thinning */
    int gn = (int)ceil(P.L / P.dmin); if (gn < 1) gn = 1;
    int nbins = gn * gn * gn;
    int *head = malloc(nbins * sizeof(int));
    int *nxt = malloc(nmax * sizeof(int));
    for (int b = 0; b < nbins; b++) head[b] = -1;
    int n = 0, fails = 0;
    double d2min = P.dmin * P.dmin;
    while (fails < 30000 && n < nmax) {
        double x = frand() * P.L, y = frand() * P.L, z = frand() * P.L;
        int bx2 = (int)(x / P.L * gn), by2 = (int)(y / P.L * gn), bz2 = (int)(z / P.L * gn);
        if (bx2 >= gn) bx2 = gn - 1; if (by2 >= gn) by2 = gn - 1; if (bz2 >= gn) bz2 = gn - 1;
        int ok = 1;
        for (int ax = bx2 - 1; ax <= bx2 + 1 && ok; ax++)
            for (int ay = by2 - 1; ay <= by2 + 1 && ok; ay++)
                for (int az = bz2 - 1; az <= bz2 + 1 && ok; az++) {
                    int wx = (ax + gn) % gn, wy2 = (ay + gn) % gn, wz = (az + gn) % gn;
                    for (int q = head[(wx * gn + wy2) * gn + wz]; q >= 0; q = nxt[q]) {
                        double dx = wr(bx[q] - x), dy = wr(by[q] - y), dz = wr(bz[q] - z);
                        if (dx*dx + dy*dy + dz*dz < d2min) { ok = 0; break; }
                    }
                }
        if (!ok) { fails++; continue; }
        bx[n] = x; by[n] = y; bz[n] = z;
        int b = (bx2 * gn + by2) * gn + bz2;
        nxt[n] = head[b]; head[b] = n;
        n++; fails = 0;
    }
    free(head); free(nxt);
    if (P.bath_frac < 1.0) {
        int m = 0;
        for (int i = 0; i < n; i++)
            if (frand() < P.bath_frac) { bx[m] = bx[i]; by[m] = by[i]; bz[m] = bz[i]; m++; }
        n = m;
    }
    return n;
}

/* 2D dart throw in the z = L/2 plane (DS medium; coplanar dynamics is
 * exact — every force the kernel applies to coplanar cells is in-plane) */
static int build_bath2d(double *bx, double *by, double *bz, int nmax)
{
    int gn = (int)ceil(P.L / P.dmin); if (gn < 1) gn = 1;
    int nbins = gn * gn;
    int *head = malloc(nbins * sizeof(int));
    int *nxt = malloc(nmax * sizeof(int));
    for (int b = 0; b < nbins; b++) head[b] = -1;
    int n = 0, fails = 0;
    double d2min = P.dmin * P.dmin;
    while (fails < 30000 && n < nmax) {
        double x = frand() * P.L, y = frand() * P.L;
        int bx2 = (int)(x / P.L * gn), by2 = (int)(y / P.L * gn);
        if (bx2 >= gn) bx2 = gn - 1;
        if (by2 >= gn) by2 = gn - 1;
        int ok = 1;
        for (int ax = bx2 - 1; ax <= bx2 + 1 && ok; ax++)
            for (int ay = by2 - 1; ay <= by2 + 1 && ok; ay++) {
                int wx = (ax + gn) % gn, wy2 = (ay + gn) % gn;
                for (int q = head[wx * gn + wy2]; q >= 0; q = nxt[q]) {
                    double dx = wr(bx[q] - x), dy = wr(by[q] - y);
                    if (dx*dx + dy*dy < d2min) { ok = 0; break; }
                }
            }
        if (!ok) { fails++; continue; }
        bx[n] = x; by[n] = y; bz[n] = 0.5 * P.L;
        int b = bx2 * gn + by2;
        nxt[n] = head[b]; head[b] = n;
        n++; fails = 0;
    }
    free(head); free(nxt);
    if (P.bath_frac < 1.0) {
        int m = 0;
        for (int i = 0; i < n; i++)
            if (frand() < P.bath_frac) { bx[m] = bx[i]; by[m] = by[i]; bz[m] = bz[i]; m++; }
        n = m;
    }
    return n;
}

static double total_energy(void);
static void geo_stats(double *phi, double *zl, int *nla, int *nld,
                      double *dbar, double *sig_d, double *maxldd);

/* ---- FCS v3 chunked stream (FCS.md) — same bytes as the Go kernel.
 * One CELL+LINK frame pair per snapshot; ANLZ instrumentation frames
 * interleave freely. Pure output — consumes no RNG. ---- */
static FILE *fcs_stream = NULL;

static void fcs_wu32(FILE *f, uint32_t v) { fwrite(&v, 4, 1, f); }
static void fcs_wu64(FILE *f, uint64_t v) { fwrite(&v, 8, 1, f); }
static void fcs_wf64(FILE *f, double v) { fwrite(&v, 8, 1, f); }

static const char *FCS_CELL_COLS[] =
    {"x","y","z","r","es","em","ee","xload","tag","fa1","fa2","th2"};
#define FCS_NCELLCOLS 12
static const char *FCS_LINK_COLS[] = {"d","A","lem","gg"};
#define FCS_NLINKCOLS 4

static void fcs_begin(FILE *f)
{
    fwrite("FCS1", 1, 4, f);
    fcs_wu32(f, 3);
    char cfg[2048];
    int n = snprintf(cfg, sizeof cfg,
        "exp=%s seed=%lu L=%g dt=%g T=%g\n"
        "laws_V2g: C=%g w1=%g w2=%g q_detune=%g gamma_res=%g gamma_res_m=%g "
        "p_gate=%g lock_floor=%g k_dep=%g k_dep_m=%g cap=%g e_s0=%g es_floor=%g "
        "e_cond=%g f_conv=%g f_evap=%g s_pull=%g kappa_lock=%g kappa_align=%g "
        "kappa_reac=%g s_k=%g s_disp=%g sigma_tumble=%g comb_limit=%d rough_k=%g "
        "gamma_rough=%g mob_sym=%d mob_floor=%g field_J=%g quant_A0=%g quant_mode=%d\n"
        "geometry: cfac=%g k_rep=%g mob_geo=%g kappa_bond=%g freeze_geo=%d "
        "bath=%d bath_frac=%g\n",
        P.exp, P.seed, P.L, P.dt, P.T,
        P.C, P.w1, P.w2, P.q_detune, P.gamma_res, P.gamma_res_m,
        P.p_gate, P.lock_floor, P.k_dep, P.k_dep_m, P.cap, P.e_s0, P.es_floor,
        P.e_cond, P.f_conv, P.f_evap, P.s_pull, P.kappa_lock, P.kappa_align,
        P.kappa_reac, P.s_k, P.s_disp, P.sigma_tumble, P.comb_limit, P.rough_k,
        P.gamma_rough, P.mob_sym, P.mob_floor, P.field_J, P.quant_A0, P.quant_mode,
        P.cfac, P.k_rep, P.mob_geo, P.kappa_bond, P.freeze_geo,
        P.bath, P.bath_frac);
    fwrite("CFG ", 1, 4, f); fcs_wu64(f, (uint64_t)n); fwrite(cfg, 1, n, f);
    unsigned char sb[512];
    int o = 0;
    uint32_t ncc = FCS_NCELLCOLS;
    memcpy(sb + o, &ncc, 4); o += 4;
    for (int k = 0; k < FCS_NCELLCOLS; k++) {
        size_t l = strlen(FCS_CELL_COLS[k]);
        sb[o++] = (unsigned char)l;
        memcpy(sb + o, FCS_CELL_COLS[k], l); o += (int)l;
    }
    uint32_t nlc = FCS_NLINKCOLS;
    memcpy(sb + o, &nlc, 4); o += 4;
    for (int k = 0; k < FCS_NLINKCOLS; k++) {
        size_t l = strlen(FCS_LINK_COLS[k]);
        sb[o++] = (unsigned char)l;
        memcpy(sb + o, FCS_LINK_COLS[k], l); o += (int)l;
    }
    fwrite("SCHM", 1, 4, f); fcs_wu64(f, (uint64_t)o); fwrite(sb, 1, o, f);
}

static void fcs_cell_frame(FILE *f)
{
    uint64_t plen = 8 + 8 + 4 + (uint64_t)NC * FCS_NCELLCOLS * 4;
    fwrite("CELL", 1, 4, f); fcs_wu64(f, plen);
    fcs_wf64(f, sim_t); fcs_wf64(f, P.L); fcs_wu32(f, (uint32_t)NC);
    for (int i = 0; i < NC; i++) {
        float rec[FCS_NCELLCOLS];
        rec[0] = (float)px_[i]; rec[1] = (float)py_[i]; rec[2] = (float)pz_[i];
        rec[3] = (float)cr[i]; rec[4] = (float)Es[i]; rec[5] = (float)Em[i];
        rec[6] = (float)Ee[i]; rec[7] = (float)((Em[i] + flload[i]) / P.cap);
        rec[8] = (float)tag[i];
        rec[9] = (float)fa1[i]; rec[10] = (float)fa2[i]; rec[11] = (float)th2[i];
        fwrite(rec, 4, FCS_NCELLCOLS, f);
    }
    uint32_t nl = 0;
    for (int s = 0; s < NSLOT; s++) if (sst[s] != S_FREE) nl++;
    uint64_t llen = 8 + 4 + (uint64_t)nl * (8 + FCS_NLINKCOLS * 4);
    fwrite("LINK", 1, 4, f); fcs_wu64(f, llen);
    fcs_wf64(f, sim_t); fcs_wu32(f, nl);
    for (int s = 0; s < NSLOT; s++) {
        if (sst[s] == S_FREE) continue;
        int i = sli[s], j = slj[s];
        uint32_t ij[2] = { (uint32_t)i, (uint32_t)j };
        double gf = gate_of(slq[s]*th2[i] - slq[s]*w2e[i]*sd[s]/P.C - slp[s]*th2[j]);
        double gb = gate_of(slp[s]*th2[j] - slp[s]*w2e[j]*sd[s]/P.C - slq[s]*th2[i]);
        float lr[FCS_NLINKCOLS];
        lr[0] = (float)sd[s]; lr[1] = (float)sA[s];
        lr[2] = (float)(slem[2*s] + slem[2*s+1]);
        lr[3] = (float)(gf * gb);
        fwrite(ij, 4, 2, f);
        fwrite(lr, 4, FCS_NLINKCOLS, f);
    }
    fflush(f);
}

static void fcs_anlz_table(FILE *f, const char *name, const char **cols,
                           int ncols, const double *rows, int nrows)
{
    size_t nn = strlen(name);
    uint64_t plen = 1 + 1 + nn + 8 + 4 + 4 + (uint64_t)nrows * ncols * 8;
    for (int c = 0; c < ncols; c++) plen += 1 + strlen(cols[c]);
    fwrite("ANLZ", 1, 4, f); fcs_wu64(f, plen);
    unsigned char kind = 1, l = (unsigned char)nn;
    fwrite(&kind, 1, 1, f);
    fwrite(&l, 1, 1, f); fwrite(name, 1, nn, f);
    fcs_wf64(f, sim_t);
    fcs_wu32(f, (uint32_t)ncols);
    for (int c = 0; c < ncols; c++) {
        unsigned char cl = (unsigned char)strlen(cols[c]);
        fwrite(&cl, 1, 1, f); fwrite(cols[c], 1, cl, f);
    }
    fcs_wu32(f, (uint32_t)nrows);
    fwrite(rows, 8, (size_t)nrows * ncols, f);
    fflush(f);
}

static void fcs_anlz_text(FILE *f, const char *txt)
{
    uint64_t plen = 1 + strlen(txt);
    fwrite("ANLZ", 1, 4, f); fcs_wu64(f, plen);
    unsigned char kind = 0;
    fwrite(&kind, 1, 1, f);
    fwrite(txt, 1, strlen(txt), f);
    fflush(f);
}

/* per-snapshot instrumentation: scalar meters, and for exp=slit the
 * ACCUMULATING screen profile — the fringe build-up lives in the stream */
static void fcs_instrument(FILE *f)
{
    double Etot = total_energy();
    double phi, zl, dbar, sig, mld; int nla, nld;
    geo_stats(&phi, &zl, &nla, &nld, &dbar, &sig, &mld);
    static const char *MC[] = {"drift","phi","z_live","nla","nld","births","deaths"};
    double mr[7];
    mr[0] = (Etot - E0_total) / (E0_total != 0 ? E0_total : 1);
    mr[1] = phi; mr[2] = zl; mr[3] = (double)nla; mr[4] = (double)nld;
    mr[5] = (double)births; mr[6] = (double)deaths;
    fcs_anlz_table(f, "meters", MC, 7, mr, 1);
    if (!strcmp(P.exp, "slit")) {
        static const char *SC[] = {"y","I"};
        double rows[DS_NBIN * 2];
        for (int b = 0; b < DS_NBIN; b++) {
            rows[2*b] = (b + 0.5) * P.L / DS_NBIN;
            rows[2*b+1] = ds_I[b];
        }
        fcs_anlz_table(f, "ds_screen", SC, 2, rows, DS_NBIN);
    }
}

/* per-frame v3 file under snap_dir (a v3 stream with one frame) */
static void write_fcs(int idx)
{
    char path[512];
    snprintf(path, sizeof path, "%s/snap_%05d_t%.3f.fcs", P.snap_dir, idx, sim_t);
    FILE *f = fopen(path, "wb");
    if (!f) { fprintf(stderr, "# WARN snap open %s failed\n", path); return; }
    fcs_begin(f);
    fcs_cell_frame(f);
    fclose(f);
}

static void jam_settle(void)
{
    /* pure repulsion to jamming (LIVEFAB relax), geometry only.
     * Witness printed: this is the §5a "configuration announces itself". */
    int sw;
    double mf = 0, mf_prev = 1e300, jk = P.jam_k;
    for (sw = 0; sw < P.jam_sweeps; sw++) {
        bins_build();
        for (int i = 0; i < NC; i++) { fxb[i] = fyb[i] = fzb[i] = 0; }
        int g = bin_g;
        for (int i = 0; i < NC; i++) {
            int bx = (int)(px_[i] / bin_sz), by = (int)(py_[i] / bin_sz), bz = (int)(pz_[i] / bin_sz);
            if (bx >= g) bx = g - 1; if (by >= g) by = g - 1; if (bz >= g) bz = g - 1;
            if (bx < 0) bx = 0; if (by < 0) by = 0; if (bz < 0) bz = 0;
            for (int ax = bx - 1; ax <= bx + 1; ax++)
                for (int ay = by - 1; ay <= by + 1; ay++)
                    for (int az = bz - 1; az <= bz + 1; az++) {
                        int wx = (ax + g) % g, wy2 = (ay + g) % g, wz = (az + g) % g;
                        for (int q = bin_head[(wx * g + wy2) * g + wz]; q >= 0; q = bin_next[q]) {
                            if (q <= i) continue;
                            double dx = wr(px_[q]-px_[i]), dy = wr(py_[q]-py_[i]), dz = wr(pz_[q]-pz_[i]);
                            double d2 = dx*dx + dy*dy + dz*dz;
                            double contact = cr[i] + cr[q];
                            if (d2 >= contact * contact || d2 <= 0) continue;
                            double d = sqrt(d2);
                            double push = jk * (contact - d);
                            if (push > 0.3) push = 0.3;   /* bounded shove */
                            double f = d > 1e-6 ? push / d : 0.0;
                            fxb[i] -= f * dx; fyb[i] -= f * dy; fzb[i] -= f * dz;
                            fxb[q] += f * dx; fyb[q] += f * dy; fzb[q] += f * dz;
                        }
                    }
        }
        mf = 0;
        for (int i = 0; i < NC; i++) {
            px_[i] = fold(px_[i] + fxb[i]);
            py_[i] = fold(py_[i] + fyb[i]);
            pz_[i] = fold(pz_[i] + fzb[i]);
            double f = sqrt(fxb[i]*fxb[i] + fyb[i]*fyb[i] + fzb[i]*fzb[i]);
            if (f > mf) mf = f;
        }
        /* adaptive damping: a rising witness means overshoot — soften */
        if (mf > mf_prev && jk > 1e-3) jk *= 0.7;
        mf_prev = mf;
        if (mf < P.jam_tol) { sw++; break; }
    }
    printf("# JAM sweeps=%d max_step=%.3e (tol %.1e, k_end=%.4f) — witness for the settle\n",
           sw, mf, P.jam_tol, jk);
}

/* load a voice to occupancy x: Em = x*cap exactly, space pulled
 * (cellfab.c:1180 pair recipe verbatim) */
static void load_voice(int i, double x)
{
    double add = x * P.cap / (1.0 + P.s_pull);
    Em[i] += add;
    double pull = P.s_pull * add;
    double avail = Es[i] - P.es_floor;
    if (pull > avail) pull = avail > 0 ? avail : 0;
    Es[i] -= pull;
    Em[i] += pull;
}

static double dstar_of(double x, int pp, int qq, int m)
{
    /* the pair-separation ladder at equal load x:
     * (q*w_i + p*w_j)*d/C = 2 pi m with w = w2/(1+q_detune*x) */
    double w = P.w2 / (1.0 + P.q_detune * x);
    return TWO_PI * m * P.C / ((qq + pp) * w);
}

static void seed_cell_defaults(int i)
{
    cr0[i] = P.r0;
    cr[i] = cr0[i];
    rand_unit(&n1x[i], &n1y[i], &n1z[i]);
    rand_unit(&n2x[i], &n2y[i], &n2z[i]);
    th2[i] = frand() * TWO_PI;
    cbeta[i] = frand() * TWO_PI;
    Es[i] = P.e_s0;
}

static void normals_transverse(int i, double ux, double uy, double uz)
{
    /* both plane normals perpendicular to u (prealign for a link) */
    double ax = fabs(ux) < 0.9 ? 1 : 0, ay = fabs(ux) < 0.9 ? 0 : 1, az = 0;
    double v1x = uy * az - uz * ay, v1y = uz * ax - ux * az, v1z = ux * ay - uy * ax;
    double n = sqrt(v1x*v1x + v1y*v1y + v1z*v1z);
    v1x /= n; v1y /= n; v1z /= n;
    double v2x = uy * v1z - uz * v1y, v2y = uz * v1x - ux * v1z, v2z = ux * v1y - uy * v1x;
    n1x[i] = v1x; n1y[i] = v1y; n1z[i] = v1z;
    n2x[i] = v2x; n2y[i] = v2y; n2z[i] = v2z;
}

/* pair/truss registry for reporting */
static int rep_pair_i = -1, rep_pair_j = -1;
static int truss_n = 0;
static int truss_e[64][2];
static double truss_dstar[64];
static double truss_shape0[3];
/* FCQ triangle registry: declared chart (p,q) per cycle edge, cycle
 * direction, and the second triangle for exp=tri2 */
static int tri_on = 0, tri_v[2][3];
static signed char tri_p[2][3], tri_q[2][3];
static int ntri = 0;

/* the fixed-chart edge defect psi_e = wrap(q*th_src - q*w_src*d/C - p*th_dst)
 * and cumulative circulation around a declared 3-cycle */
static void tri_meters(int t, double psi[3], double *ggm, double *circ)
{
    *ggm = 1; *circ = 0;
    for (int e = 0; e < 3; e++) {
        int a = tri_v[t][e], b = tri_v[t][(e + 1) % 3];
        int s = hfind(a < b ? a : b, a < b ? b : a);
        if (s < 0 || sst[s] == S_FREE) { psi[e] = 99; *ggm = 0; continue; }
        double p = tri_p[t][e], q = tri_q[t][e];
        psi[e] = wrap_pi(q * th2[a] - q * w2e[a] * sd[s] / P.C - p * th2[b]);
        double gf = gate_of(q * th2[a] - q * w2e[a] * sd[s] / P.C - p * th2[b]);
        double gb = gate_of(p * th2[b] - p * w2e[b] * sd[s] / P.C - q * th2[a]);
        *ggm *= sqrt(gf * gb);
        /* cycle-direction flux minus counter-flux: a->b is dir0 iff a==sli */
        int fwd = (sli[s] == a) ? 0 : 1;
        *circ += sfluxd[2*s + fwd] - sfluxd[2*s + (1 - fwd)];
    }
}

static void gyration_eigs(double out[3])
{
    /* sorted eigenvalues of the tagged gyration tensor (shape metric) */
    double cx2 = 0, cy2 = 0, cz2 = 0; int nt = 0;
    for (int i = 0; i < NC; i++)
        if (tag[i]) { cx2 += cenx + wr(px_[i]-cenx); cy2 += ceny + wr(py_[i]-ceny);
                      cz2 += cenz + wr(pz_[i]-cenz); nt++; }
    if (!nt) { out[0] = out[1] = out[2] = 0; return; }
    cx2 /= nt; cy2 /= nt; cz2 /= nt;
    double m[6] = {0};
    for (int i = 0; i < NC; i++) {
        if (!tag[i]) continue;
        double dx = wr(px_[i]-cx2), dy = wr(py_[i]-cy2), dz = wr(pz_[i]-cz2);
        m[0] += dx*dx; m[1] += dy*dy; m[2] += dz*dz;
        m[3] += dx*dy; m[4] += dx*dz; m[5] += dy*dz;
    }
    for (int k = 0; k < 6; k++) m[k] /= nt;
    /* eigenvalues of symmetric 3x3 (analytic, Smith) */
    double p1 = m[3]*m[3] + m[4]*m[4] + m[5]*m[5];
    double q = (m[0] + m[1] + m[2]) / 3.0;
    double p2 = (m[0]-q)*(m[0]-q) + (m[1]-q)*(m[1]-q) + (m[2]-q)*(m[2]-q) + 2*p1;
    double p = sqrt(p2 / 6.0);
    if (p < 1e-15) { out[0] = out[1] = out[2] = q; return; }
    double B[6];
    B[0] = (m[0]-q)/p; B[1] = (m[1]-q)/p; B[2] = (m[2]-q)/p;
    B[3] = m[3]/p; B[4] = m[4]/p; B[5] = m[5]/p;
    double detB = B[0]*(B[1]*B[2]-B[5]*B[5]) - B[3]*(B[3]*B[2]-B[5]*B[4]) + B[4]*(B[3]*B[5]-B[1]*B[4]);
    double r = detB / 2.0;
    if (r < -1) r = -1; if (r > 1) r = 1;
    double phi = acos(r) / 3.0;
    double e1 = q + 2*p*cos(phi);
    double e3 = q + 2*p*cos(phi + 2.0*M_PI/3.0);
    double e2 = 3*q - e1 - e3;
    out[0] = e3; out[1] = e2; out[2] = e1;   /* ascending */
}

/* ------------------------------------------------------------------ */
/* main                                                                */
/* ------------------------------------------------------------------ */

int main(int argc, char **argv)
{
    setvbuf(stdout, NULL, _IOLBF, 0);
    cfg_defaults();
    int argi = 1;
    if (argc > 1 && !strchr(argv[1], '=')) { load_cfg(argv[1]); argi = 2; }
    for (; argi < argc; argi++) {
        char *eq = strchr(argv[argi], '=');
        if (!eq) continue;
        *eq = 0;
        set_kv(argv[argi], eq + 1);
    }

    comb_build();
    rng_s = P.seed ? P.seed : 1;
    for (int w = 0; w < 8; w++) xrand();   /* cellfab warm-up */

    gm_ = P.gamma_res_m < 0 ? P.gamma_res : P.gamma_res_m;
    G2m_ = gm_ * gm_;
    lockf_ = P.lock_floor * P.cap;
    Aref_ = M_PI * P.r0 * P.r0;
    dref_ = 2.0 * P.r0;
    A0eff = P.quant_A0 < 0 ? P.e_s0 * 1.5053 / P.C : P.quant_A0;

    printf("# freecell — v89 free-cell substance test (exp=%s)\n", P.exp);
    printf("# laws = laws_V2g VERBATIM (defaults): C=%g w1=%g w2=%g q_detune=%g\n",
           P.C, P.w1, P.w2, P.q_detune);
    printf("# gamma_res=%g gamma_res_m=%g p_gate=%g lock_floor=%g k_dep=%g k_dep_m=%g cap=%g\n",
           P.gamma_res, P.gamma_res_m, P.p_gate, P.lock_floor, P.k_dep, P.k_dep_m, P.cap);
    printf("# e_s0=%g es_floor=%g e_cond=%g f_conv=%g f_evap=%g s_pull=%g\n",
           P.e_s0, P.es_floor, P.e_cond, P.f_conv, P.f_evap, P.s_pull);
    printf("# kappa_lock=%g kappa_align=%g s_k=%g s_disp=%g sigma_tumble=%g\n",
           P.kappa_lock, P.kappa_align, P.s_k, P.s_disp, P.sigma_tumble);
    printf("# comb_limit=%d rough_k=%g gamma_rough=%g mob_sym=%d mob_floor=%g field_J=%g\n",
           P.comb_limit, P.rough_k, P.gamma_rough, P.mob_sym, P.mob_floor, P.field_J);
    printf("# quant_A0=%g quant_mode=%d (A0eff=%g)\n", P.quant_A0, P.quant_mode, A0eff);
    printf("# GEOMETRY (apparatus): cfac=%g k_rep=%g mob_geo=%g kappa_bond=%g freeze_geo=%d\n",
           P.cfac, P.k_rep, P.mob_geo, P.kappa_bond, P.freeze_geo);
    printf("# bath=%d bath_frac=%g jam_sweeps=%d jam_k=%g L=%g dt=%g T=%g seed=%lu\n",
           P.bath, P.bath_frac, P.jam_sweeps, P.jam_k, P.L, P.dt, P.T, P.seed);

    /* ---------------- build ---------------- */
    double V = P.L * P.L * P.L;
    int nmax = (int)(0.8 * V / (P.dmin * P.dmin * P.dmin)) + 4096;
    double *bx = malloc(nmax * sizeof(double));
    double *by = malloc(nmax * sizeof(double));
    double *bz = malloc(nmax * sizeof(double));
    int nb = 0;
    if (P.bath) nb = !strcmp(P.exp, "slit")
        ? build_bath2d(bx, by, bz, nmax)
        : build_bath(bx, by, bz, nmax);

    int extra = 0;
    if (!strcmp(P.exp, "pair") && !P.bath) extra = 2;
    if (!strcmp(P.exp, "ring")) extra = P.ring_n;
    if (!strcmp(P.exp, "oct")) extra = 6;
    if (!strcmp(P.exp, "tri")) extra = 3;
    if (!strcmp(P.exp, "tri2")) extra = 6;

    alloc_all(nb + extra);
    for (int i = 0; i < nb; i++) {
        px_[i] = bx[i]; py_[i] = by[i]; pz_[i] = bz[i];
        /* cellfab per-cell draw order: cr0 jitter, n1, n2, th (th1 skipped:
         * field phase lives in fa), th2, cbeta */
        cr0[i] = P.r0 * (1.0 + P.rjit * (2.0 * frand() - 1.0));
        cr[i] = cr0[i];
        rand_unit(&n1x[i], &n1y[i], &n1z[i]);
        rand_unit(&n2x[i], &n2y[i], &n2z[i]);
        (void)frand();               /* th1 draw consumed, stream aligned */
        th2[i] = frand() * TWO_PI;
        cbeta[i] = frand() * TWO_PI;
        Es[i] = P.e_s0;
    }
    free(bx); free(by); free(bz);
    NC = nb + extra;

    double cx0 = P.px < 0 ? 0.5 * P.L : P.px;
    double cy0 = P.py < 0 ? 0.5 * P.L : P.py;
    double cz0 = P.pz < 0 ? 0.5 * P.L : P.pz;
    cenx = cx0; ceny = cy0; cenz = cz0;

    /* jam-settle the bath BEFORE seeding (livefab: relaxed to jamming).
     * Only the nb bath cells exist yet — truss cells are appended after
     * the settle, so the jam must not see their uninitialized slots. */
    if (nb > 0) {
        int NCfull = NC;
        NC = nb;
        jam_settle();
        NC = NCfull;
    }
    else printf("# JAM skipped (no bath)\n");

    /* ---------------- seeds ---------------- */
    if (!strcmp(P.exp, "pair")) {
        int i, j;
        double x = P.pair_x0;
        double x1 = P.pair_x1 < 0 ? x : P.pair_x1;
        /* the ladder with per-voice pitches: (q w_i + p w_j) d / C = 2 pi m */
        double wA = P.w2 / (1.0 + P.q_detune * x);
        double wB = P.w2 / (1.0 + P.q_detune * x1);
        double ds = TWO_PI * P.pair_m * P.C / (P.pair_qq * wA + P.pair_pp * wB);
        double d0 = ds + P.pair_doff;
        if (!P.bath) {
            i = nb; j = nb + 1;
            seed_cell_defaults(i); seed_cell_defaults(j);
            px_[i] = cx0 - 0.5 * d0; py_[i] = cy0; pz_[i] = cz0;
            px_[j] = cx0 + 0.5 * d0; py_[j] = cy0; pz_[j] = cz0;
        } else {
            /* pick the two bath cells whose pair straddles the centre with
             * separation nearest d0 (e7-style: the pair is OF the bath) */
            int bi = -1, bj = -1; double bcost = 1e300;
            for (int a = 0; a < nb; a++)
                for (int b2 = a + 1; b2 < nb; b2++) {
                    double dx = wr(px_[b2]-px_[a]), dy = wr(py_[b2]-py_[a]), dz = wr(pz_[b2]-pz_[a]);
                    double d = sqrt(dx*dx + dy*dy + dz*dz);
                    double mx = wr(px_[a]-cx0) + 0.5*dx, my = wr(py_[a]-cy0) + 0.5*dy,
                           mz = wr(pz_[a]-cz0) + 0.5*dz;
                    double cost = fabs(d - d0) * 10 + sqrt(mx*mx+my*my+mz*mz);
                    if (d > 0.9 * d0 && d < 1.3 * d0 && cost < bcost)
                        { bcost = cost; bi = a; bj = b2; }
                }
            i = bi; j = bj;
            if (i < 0) { fprintf(stderr, "no bath pair found near d0\n"); return 1; }
        }
        double dx = wr(px_[j]-px_[i]), dy = wr(py_[j]-py_[i]), dz = wr(pz_[j]-pz_[i]);
        double d = sqrt(dx*dx + dy*dy + dz*dz);
        normals_transverse(i, dx/d, dy/d, dz/d);
        normals_transverse(j, dx/d, dy/d, dz/d);
        if (x > 0) load_voice(i, x);
        if (x1 > 0) load_voice(j, x1);
        if (P.seedlock) {
            th2[i] = frand() * TWO_PI;
            th2[j] = fmod((P.pair_qq * (th2[i] - wA * d / P.C)) / P.pair_pp
                          + 8.0 * TWO_PI, TWO_PI);
        }
        tag[i] = tag[j] = 1;
        rep_pair_i = i; rep_pair_j = j;
        printf("# SEED pair: i=%d j=%d x=%.4f/%.4f p:q=%d:%d m=%d d*=%.6f d0=%.6f (doff %+0.4f)\n",
               i, j, x, x1, P.pair_pp, P.pair_qq, P.pair_m, ds, d, P.pair_doff);
        printf("# SEED pair radii r=%.4f contact=%.4f  Em/voice=%.4f\n",
               cr0[i] * cbrt(Es[i]), cr0[i]*cbrt(Es[i]) + cr0[j]*cbrt(Es[j]), Em[i]);
    }
    else if (!strcmp(P.exp, "ring") || !strcmp(P.exp, "oct")) {
        int nvo = !strcmp(P.exp, "oct") ? 6 : P.ring_n;
        double x = !strcmp(P.exp, "oct") ? P.oct_x : P.ring_x;
        double doff = !strcmp(P.exp, "oct") ? P.oct_doff : P.ring_doff;
        double ds = dstar_of(x, 1, 1, P.pair_m);
        double de = ds + doff;
        /* bath embedding: carve a cavity — remove bath cells within dmin
         * of any truss site, compact, then append the truss */
        if (P.bath && nb > 0) {
            double R = de / (2.0 * sin(M_PI / (nvo > 2 ? nvo : 3)));
            double clear = R + P.dmin + 0.2;
            int m = 0;
            for (int i = 0; i < nb; i++) {
                double dx = wr(px_[i]-cx0), dy = wr(py_[i]-cy0), dz = wr(pz_[i]-cz0);
                double keep_r2 = dx*dx + dy*dy + dz*dz;
                int keep = 1;
                if (keep_r2 < clear * clear) {
                    /* fine test against actual truss sites below is
                     * shape-specific; the spherical clear is enough */
                    keep = 0;
                }
                if (keep) {
                    if (m != i) {
                        px_[m]=px_[i]; py_[m]=py_[i]; pz_[m]=pz_[i];
                        cr0[m]=cr0[i]; cr[m]=cr[i];
                        n1x[m]=n1x[i]; n1y[m]=n1y[i]; n1z[m]=n1z[i];
                        n2x[m]=n2x[i]; n2y[m]=n2y[i]; n2z[m]=n2z[i];
                        th2[m]=th2[i]; cbeta[m]=cbeta[i];
                        Es[m]=Es[i]; Em[m]=Em[i]; Ee[m]=Ee[i];
                        fa1[m]=fa1[i]; fa2[m]=fa2[i];
                    }
                    m++;
                }
            }
            printf("# CARVE removed %d bath cells (clear=%.2f)\n", nb - m, clear);
            nb = m;
            NC = nb + extra;
        }
        int base_i = nb;
        truss_n = 0;
        if (!strcmp(P.exp, "oct")) {
            double a = de / sqrt(2.0);
            double vx[6] = { a, -a, 0, 0, 0, 0 };
            double vy[6] = { 0, 0, a, -a, 0, 0 };
            double vz[6] = { 0, 0, 0, 0, a, -a };
            for (int k = 0; k < 6; k++) {
                int i = base_i + k;
                seed_cell_defaults(i);
                px_[i] = cx0 + vx[k]; py_[i] = cy0 + vy[k]; pz_[i] = cz0 + vz[k];
                /* global-aligned normals: dense plane n2 = z, n1 = x */
                n1x[i] = 1; n1y[i] = 0; n1z[i] = 0;
                n2x[i] = 0; n2y[i] = 0; n2z[i] = 1;
                load_voice(i, x);
                tag[i] = 1;
            }
            for (int a2 = 0; a2 < 6; a2++)
                for (int b2 = a2 + 1; b2 < 6; b2++) {
                    double dx = px_[base_i+b2]-px_[base_i+a2], dy = py_[base_i+b2]-py_[base_i+a2],
                           dz = pz_[base_i+b2]-pz_[base_i+a2];
                    double d = sqrt(dx*dx+dy*dy+dz*dz);
                    if (d < 1.3 * de) {
                        truss_e[truss_n][0] = base_i + a2;
                        truss_e[truss_n][1] = base_i + b2;
                        truss_dstar[truss_n] = ds;
                        truss_n++;
                    }
                }
        } else {
            double R = de / (2.0 * sin(M_PI / nvo));
            for (int k = 0; k < nvo; k++) {
                int i = base_i + k;
                seed_cell_defaults(i);
                double a = TWO_PI * k / nvo;
                px_[i] = cx0 + R * cos(a); py_[i] = cy0 + R * sin(a); pz_[i] = cz0;
                /* COPLANAR normals: both planes = the ring plane. The
                 * collimation axis (the planes' intersection) is then
                 * unconstrained in-plane: axi = 1 and dnn = 1 on every
                 * in-plane edge — gpl = 1 exactly. (A radial n1 gave
                 * gpl = 0.56 and the bond lost to repulsion by t~40.) */
                n1x[i] = 0; n1y[i] = 0; n1z[i] = 1;
                n2x[i] = 0; n2y[i] = 0; n2z[i] = 1;
                load_voice(i, x);
                tag[i] = 1;
            }
            for (int k = 0; k < nvo; k++) {
                truss_e[truss_n][0] = base_i + k;
                truss_e[truss_n][1] = base_i + (k + 1) % nvo;
                truss_dstar[truss_n] = ds;
                truss_n++;
            }
        }
        /* phases: lock around the edge list (unison chain) */
        if (P.seedlock) {
            double we = P.w2 / (1.0 + P.q_detune * x);
            th2[base_i] = frand() * TWO_PI;
            for (int e = 0; e < truss_n; e++) {
                int i = truss_e[e][0], j = truss_e[e][1];
                double dx = px_[j]-px_[i], dy = py_[j]-py_[i], dz = pz_[j]-pz_[i];
                double d = sqrt(dx*dx+dy*dy+dz*dz);
                th2[j] = fmod(th2[i] - we * d / P.C + 8.0 * TWO_PI, TWO_PI);
            }
        }
        printf("# SEED %s: n=%d x=%.4f d*=%.6f d_edge=%.6f edges=%d\n",
               P.exp, nvo, x, ds, de, truss_n);
    }
    else if (!strcmp(P.exp, "blob")) {
        /* e3a-class Gaussian dense blob (cellfab.c:1754 recipe) */
        int nload = 0;
        for (int i = 0; i < NC; i++) {
            double dx = wr(px_[i]-cx0), dy = wr(py_[i]-cy0), dz = wr(pz_[i]-cz0);
            double r2 = dx*dx + dy*dy + dz*dz;
            double add = P.amp * exp(-r2 / (2.0 * P.sigma * P.sigma));
            if (add < 1e-4) continue;
            if (Em[i] + add > 0.95 * P.cap) add = 0.95 * P.cap - Em[i];
            if (add <= 0) continue;
            Em[i] += add;
            double pull = P.s_pull * add;
            double avail = Es[i] - P.es_floor;
            if (pull > avail) pull = avail > 0 ? avail : 0;
            Es[i] -= pull;
            Em[i] += pull;
            /* e3b tilt: phase ramp along x (cellfab.c:1750 convention) */
            th2[i] = P.kx != 0 ? fmod(-(P.kx * dx) + 8.0 * TWO_PI, TWO_PI) : 0;
            if (add > 0.05 * P.amp) { tag[i] = 1; nload++; }
        }
        if (P.kx != 0)
            for (int i = 0; i < NC; i++) normals_transverse(i, 1, 0, 0);
        printf("# SEED blob: amp=%g sigma=%g kx=%g tagged=%d (add>5%% of amp)\n",
               P.amp, P.sigma, P.kx, nload);
    }
    else if (!strcmp(P.exp, "tri") || !strcmp(P.exp, "tri2")) {
        /* FCQ — the fifth-triangle: the incomplete harmonic as binding.
         * UUD: pitches (w, w, 2w/3); the ladder identity S_UD = 4 w_U
         * makes d*_UD(m=2) = d*_UU(m=1) = pi/w_U EXACTLY: equilateral.
         * UDD: isoceles; its D-D edge d* = pi/w_D = 1.5*pi/w_U exceeds
         * contact -> predicted open chain. The Z3 branch k seeds the
         * D-voice phase (th_D = (2(th_U - w_U d) + 2 pi k)/3). */
        double xU = P.tri_xU;
        double xD = P.tri_xD < 0 ? (0.5 + 1.8 * xU) / 1.2 : P.tri_xD;
        double wU = P.w2 / (1.0 + P.q_detune * xU);
        double wD = P.w2 / (1.0 + P.q_detune * xD);
        double dUU = M_PI / wU;                      /* m=1 unison rung */
        double dDD = M_PI / wD;                      /* m=1 for D-D     */
        double dUD = 2.0 * TWO_PI / (2.0 * wU + 3.0 * wD);   /* m=2 fifth */
        ntri = !strcmp(P.exp, "tri2") ? 2 : 1;
        tri_on = 1;
        int base_i = nb;
        truss_n = 0;
        for (int t = 0; t < ntri; t++) {
            double ox = t == 0 ? 0.0 : P.tri2_sep;
            int kbr = t == 0 ? P.tri_branch : P.tri2_k2;
            int v0 = base_i + 3 * t;
            for (int kk = 0; kk < 3; kk++) {
                int i = v0 + kk;
                seed_cell_defaults(i);
                n1x[i] = 0; n1y[i] = 0; n1z[i] = 1;  /* coplanar, gpl=1 */
                n2x[i] = 0; n2y[i] = 0; n2z[i] = 1;
                tag[i] = 1;
                tri_v[t][kk] = i;
            }
            if (P.tri_kind == 0) {
                /* UUD equilateral, edge dUU (= dUD to FP): vertices at
                 * 90/210/330 degrees, circumradius d/sqrt(3) */
                double de = dUU + P.tri_doff, R = de / sqrt(3.0);
                for (int kk = 0; kk < 3; kk++) {
                    double a = M_PI / 2.0 + TWO_PI * kk / 3.0;
                    px_[v0 + kk] = fold(cx0 + ox + R * cos(a));
                    py_[v0 + kk] = fold(cy0 + R * sin(a));
                    pz_[v0 + kk] = cz0;
                }
                load_voice(v0 + 0, xU);
                load_voice(v0 + 1, xU);
                load_voice(v0 + 2, xD);
                /* chart: 0->1 unison, 1->2 fifth (src U: p=3,q=2),
                 * 2->0 inverse fifth (src D: p=2,q=3) */
                tri_p[t][0] = 1; tri_q[t][0] = 1;
                tri_p[t][1] = 3; tri_q[t][1] = 2;
                tri_p[t][2] = 2; tri_q[t][2] = 3;
                if (P.seedlock) {
                    th2[v0] = frand() * TWO_PI;
                    th2[v0+1] = fmod(th2[v0] - wU * de + 8.0 * TWO_PI, TWO_PI);
                    th2[v0+2] = fmod((2.0 * (th2[v0+1] - wU * dUD)
                                      + TWO_PI * kbr) / 3.0
                                     + 8.0 * TWO_PI, TWO_PI);
                }
                truss_e[truss_n][0] = v0;     truss_e[truss_n][1] = v0 + 1;
                truss_dstar[truss_n++] = dUU;
                truss_e[truss_n][0] = v0 + 1; truss_e[truss_n][1] = v0 + 2;
                truss_dstar[truss_n++] = dUD;
                truss_e[truss_n][0] = v0 + 2; truss_e[truss_n][1] = v0;
                truss_dstar[truss_n++] = dUD;
            } else {
                /* UDD isoceles: U apex, D-D base at dDD, sides dUD */
                double hb = 0.5 * dDD;
                double hh = sqrt(dUD * dUD - hb * hb);
                px_[v0 + 0] = fold(cx0 + ox);       py_[v0 + 0] = fold(cy0 + hh); pz_[v0 + 0] = cz0;
                px_[v0 + 1] = fold(cx0 + ox - hb);  py_[v0 + 1] = cy0;            pz_[v0 + 1] = cz0;
                px_[v0 + 2] = fold(cx0 + ox + hb);  py_[v0 + 2] = cy0;            pz_[v0 + 2] = cz0;
                load_voice(v0 + 0, xU);
                load_voice(v0 + 1, xD);
                load_voice(v0 + 2, xD);
                tri_p[t][0] = 3; tri_q[t][0] = 2;   /* U->D fifth */
                tri_p[t][1] = 1; tri_q[t][1] = 1;   /* D->D unison */
                tri_p[t][2] = 2; tri_q[t][2] = 3;   /* D->U inverse */
                if (P.seedlock) {
                    th2[v0] = frand() * TWO_PI;
                    th2[v0+1] = fmod((2.0 * (th2[v0] - wU * dUD)
                                      + TWO_PI * kbr) / 3.0
                                     + 8.0 * TWO_PI, TWO_PI);
                    th2[v0+2] = fmod(th2[v0+1] - wD * dDD + 8.0 * TWO_PI, TWO_PI);
                }
                truss_e[truss_n][0] = v0;     truss_e[truss_n][1] = v0 + 1;
                truss_dstar[truss_n++] = dUD;
                truss_e[truss_n][0] = v0 + 1; truss_e[truss_n][1] = v0 + 2;
                truss_dstar[truss_n++] = dDD;
                truss_e[truss_n][0] = v0 + 2; truss_e[truss_n][1] = v0;
                truss_dstar[truss_n++] = dUD;
            }
        }
        printf("# SEED tri: kind=%s ntri=%d xU=%.4f xD=%.4f wU=%.4f wD=%.4f (wD/wU=%.6f)\n",
               P.tri_kind ? "UDD" : "UUD", ntri, xU, xD, wU, wD, wD / wU);
        printf("# SEED tri: d*_UU=%.6f d*_UD(m=2)=%.6f d*_DD=%.6f branch_k=%d doff=%+.3f\n",
               dUU, dUD, dDD, P.tri_branch, P.tri_doff);
        printf("# SEED tri: contacts UU=%.4f UD=%.4f DD=%.4f (DD edge %s)\n",
               2 * P.r0 * cbrt(Es[base_i]),
               P.r0 * (cbrt(Es[base_i]) + cbrt(1.0 - P.s_pull * xD * P.cap / (1 + P.s_pull))),
               2 * P.r0 * cbrt(1.0 - P.s_pull * xD * P.cap / (1 + P.s_pull)),
               dDD > 2 * P.r0 * cbrt(1.0 - P.s_pull * xD * P.cap / (1 + P.s_pull))
               ? "OPEN (no contact)" : "closed");
    }
    else if (!strcmp(P.exp, "pulse")) {
        /* Tier-A e2 analog: field packet with phase tilt on the LIVE
         * substrate (cellfab.c:1803 aux-packet recipe); prealigned
         * normals transverse to k̂ = x̂ */
        for (int i = 0; i < NC; i++) {
            double dx = wr(px_[i]-cx0), dy = wr(py_[i]-cy0), dz = wr(pz_[i]-cz0);
            double g = exp(-(dx*dx + dy*dy + dz*dz) / (2.0 * P.sigma * P.sigma));
            if (g < 1e-8) continue;
            double tilt = -(P.kx * dx);
            fa1[i] += sqrt(P.amp * g) * cos(tilt);
            fa2[i] += sqrt(P.amp * g) * sin(tilt);
            Ee[i] = fa1[i]*fa1[i] + fa2[i]*fa2[i];
        }
        for (int i = 0; i < NC; i++) normals_transverse(i, 1, 0, 0);
        printf("# SEED pulse: amp=%g sigma=%g kx=%g (prealigned transverse)\n",
               P.amp, P.sigma, P.kx);
    }
    else if (!strcmp(P.exp, "slit")) {
        /* DS tier-0 (DS.md): carve the vacuum wall out of the jammed 2D
         * bath (except the slit bridges), then seed a quasi-plane field
         * packet. Optics regime declared in DS.md: q_detune=0 e_cond=99
         * (given on the command line, printed in the header above). */
        double yA = cy0 - 0.5 * P.slit_sep, yB = cy0 + 0.5 * P.slit_sep;
        int m = 0;
        for (int i = 0; i < nb; i++) {
            int inwall = fabs(wr(px_[i] - P.slit_wallx)) < 0.5 * P.slit_th;
            int inA = fabs(py_[i] - yA) < P.slit_hw;
            int inB = fabs(py_[i] - yB) < P.slit_hw;
            int keep = 1;
            if (P.slit_mask != 3 && inwall) {
                keep = 0;
                if (P.slit_mask == 0 && (inA || inB)) keep = 1;
                if (P.slit_mask == 1 && inA) keep = 1;
                if (P.slit_mask == 2 && inB) keep = 1;
            }
            if (keep) {
                if (m != i) {
                    px_[m]=px_[i]; py_[m]=py_[i]; pz_[m]=pz_[i];
                    cr0[m]=cr0[i]; cr[m]=cr[i];
                    n1x[m]=n1x[i]; n1y[m]=n1y[i]; n1z[m]=n1z[i];
                    n2x[m]=n2x[i]; n2y[m]=n2y[i]; n2z[m]=n2z[i];
                    th2[m]=th2[i]; cbeta[m]=cbeta[i];
                    Es[m]=Es[i]; Em[m]=Em[i]; Ee[m]=Ee[i];
                    fa1[m]=fa1[i]; fa2[m]=fa2[i];
                }
                m++;
            }
        }
        printf("# CARVE slit wall: removed %d cells (mask=%d wall_x=%g th=%g slits y=%.2f/%.2f hw=%g)\n",
               nb - m, P.slit_mask, P.slit_wallx, P.slit_th, yA, yB, P.slit_hw);
        nb = m; NC = nb;
        /* THE WALL IS A FIXTURE (v89 DS precedent: the wall was imposed
         * apparatus — a detuned mute). A carved vacuum gap in the LIVE
         * foam heals in ~15 t.u. (measured, DS.md: the substrate defends
         * no shape without energy) — so the face cells are pinned: pass D
         * skips them; every other sector still runs on them. */
        int npin = 0;
        if (P.slit_mask != 3)
            for (int i = 0; i < NC; i++)
                if (fabs(wr(px_[i] - P.slit_wallx)) < P.slit_pinw) { pin[i] = 1; npin++; }
        printf("# PIN wall fixture: %d cells within +-%g of x=%g\n",
               npin, P.slit_pinw, P.slit_wallx);
        /* quasi-plane packet: Gaussian sigma_x = P.sigma, sigma_y = slit_sy */
        for (int i = 0; i < NC; i++) {
            double dx = wr(px_[i] - P.slit_srcx), dy = wr(py_[i] - cy0);
            double g = exp(-dx*dx / (2.0 * P.sigma * P.sigma)
                           - dy*dy / (2.0 * P.slit_sy * P.slit_sy));
            if (g < 1e-8) continue;
            double tilt = -(P.kx * dx);
            fa1[i] += sqrt(P.amp * g) * cos(tilt);
            fa2[i] += sqrt(P.amp * g) * sin(tilt);
            Ee[i] = fa1[i]*fa1[i] + fa2[i]*fa2[i];
        }
        for (int i = 0; i < NC; i++) normals_transverse(i, 1, 0, 0);
        double lam = P.kx != 0 ? TWO_PI / P.kx : 0;
        double D = P.slit_screenx - P.slit_wallx;
        printf("# SEED slit: src_x=%g sigma_x=%g sigma_y=%g kx=%g amp=%g -> lam_seed=%.4f\n",
               P.slit_srcx, P.sigma, P.slit_sy, P.kx, P.amp, lam);
        printf("# SEED slit loci: D=%.2f s=%.2f dy_smallangle=%.3f screen_x=%g gate=[%g,%g]\n",
               D, P.slit_sep, lam * D / P.slit_sep, P.slit_screenx, P.slit_t0, P.slit_t1);
    }
    /* exp=bath: nothing to seed */

    if (P.noise_amp > 0) {
        for (int i = 0; i < NC; i++) {
            double e0 = frand() * P.noise_amp;
            double ph = frand() * TWO_PI;
            fa1[i] += sqrt(e0) * cos(ph);
            fa2[i] += sqrt(e0) * sin(ph);
            Ee[i] = fa1[i]*fa1[i] + fa2[i]*fa2[i];
        }
        printf("# SEED noise: amp=%g (e1-class churn)\n", P.noise_amp);
    }

    /* initial radii + topology + ledger */
    for (int i = 0; i < NC; i++) {
        double ratio = Es[i] / P.e_s0;
        cr[i] = cr0[i] * cbrt(ratio > 0 ? ratio : 0);
    }
    topo_refresh();
    for (int i = 0; i < NC; i++) { fxb[i] = fyb[i] = fzb[i] = 0; }
    E0_total = total_energy();
    births = deaths = beta_returns = 0; beta_energy = 0;
    {
        double phi, zl, dbar, sig, mld; int nla, nld;
        geo_stats(&phi, &zl, &nla, &nld, &dbar, &sig, &mld);
        printf("# INIT NC=%d slots=%d live=%d phi_nom=%.4f z_live=%.2f dbar=%.4f sigma_d=%.2f%% E0=%.9f\n",
               NC, NSLOT, nla, phi, zl, dbar, 100 * sig, E0_total);
    }
    if (truss_n > 0) gyration_eigs(truss_shape0);

    if (P.snap_file[0]) {
        fcs_stream = fopen(P.snap_file, "wb");
        if (!fcs_stream) fprintf(stderr, "# WARN snap_file open failed\n");
        else fcs_begin(fcs_stream);
    }

    /* ---------------- run ---------------- */
    int steps = (int)(P.T / P.dt + 0.5);
    int sheared = 0;
    printf("#\n# t | drift_rel | phi_nom z_live NLa NLd births deaths | dbar sigma_d%% | maxstep maxbond | Em_tag ret com_dev rms conn z_tag | pair: d delta gf gb | core/far\n");
    double emt0 = -1;
    for (int st = 0; st <= steps; st++) {
        sim_t = st * P.dt;
        /* instantaneous deviatoric strain test (truss shear leg) */
        if (P.shear_eps > 0 && !sheared && sim_t >= P.shear_t && P.shear_t > 0) {
            double sx = 1.0 + P.shear_eps, sy = 1.0 / sqrt(1.0 + P.shear_eps),
                   sz = sy;
            for (int i = 0; i < NC; i++) {
                if (!tag[i]) continue;
                px_[i] = fold(cx0 + wr(px_[i] - cx0) * sx);
                py_[i] = fold(cy0 + wr(py_[i] - cy0) * sy);
                pz_[i] = fold(cz0 + wr(pz_[i] - cz0) * sz);
            }
            sheared = 1;
            printf("# SHEAR applied: eps=%g at t=%.2f (volume-preserving deviatoric)\n",
                   P.shear_eps, sim_t);
        }
        if (st % P.diag_every == 0 || st == steps) {
            double Etot = total_energy();
            double drift = (Etot - E0_total) / (E0_total != 0 ? E0_total : 1);
            double phi, zl, dbar, sig, mld; int nla, nld;
            geo_stats(&phi, &zl, &nla, &nld, &dbar, &sig, &mld);
            double emt = 0, comx = 0, comy = 0, comz = 0, rms = 0, conn = 0, zb = 0;
            tag_stats(&emt, &comx, &comy, &comz, &rms, &conn, &zb);
            if (emt0 < 0) emt0 = emt > 0 ? emt : -1;
            double ret = emt0 > 0 ? emt / emt0 : 0;
            double comdev = emt > 1e-12
                ? sqrt((comx-cenx)*(comx-cenx) + (comy-ceny)*(comy-ceny)
                       + (comz-cenz)*(comz-cenz)) : 0;
            double mstep = max_net_force();
            printf("t=%8.2f %+9.2e | %.4f %5.2f %6d %4d %6ld %6ld | %.4f %5.2f | %8.2e %8.2e | %8.4f %6.4f %7.4f %7.4f %5.3f %5.2f",
                   sim_t, drift, phi, zl, nla, nld, births, deaths,
                   dbar, 100 * sig, mstep, mld, emt, ret, comdev, rms, conn, zb);
            if (rep_pair_i >= 0) {
                int i = rep_pair_i, j = rep_pair_j;
                int s = hfind(i < j ? i : j, i < j ? j : i);
                if (s >= 0 && sst[s] != S_FREE) {
                    double delta = wrap_pi((slq[s] * w2e[i] + slp[s] * w2e[j]) * sd[s] / P.C);
                    double gf = gate_of(slq[s]*th2[i] - slq[s]*w2e[i]*sd[s]/P.C - slp[s]*th2[j]);
                    double gb = gate_of(slp[s]*th2[j] - slp[s]*w2e[j]*sd[s]/P.C - slq[s]*th2[i]);
                    printf(" | d=%.6f dl=%+.4f gf=%.4f gb=%.4f x=%.4f",
                           sd[s], delta, gf, gb, (Em[i]+flload[i]) / P.cap);
                } else printf(" | PAIR-CHANNEL-DEAD");
            }
            if (!strcmp(P.exp, "blob")) {
                double sh[8]; es_shells(sh);
                double far = (sh[5] + sh[6] + sh[7]) / 3.0;
                printf(" | es_core/far=%.4f", far > 0 ? sh[0] / far : 0);
                if (nfsamp < 512) {
                    fs_t[nfsamp] = sim_t;
                    fs_r[nfsamp] = comx; fs_y[nfsamp] = comy; fs_z[nfsamp] = comz;
                    nfsamp++;
                }
            }
            if (!strcmp(P.exp, "pulse")) {
                /* field front = max radius over cells with Ee > 5% amp
                 * (cellfab.c:3701 convention, verbatim) */
                double rf = 0, thr = 0.05 * (P.amp > 0 ? P.amp : 1.0);
                for (int i2 = 0; i2 < NC; i2++) {
                    if (Ee[i2] <= thr) continue;
                    double dx = wr(px_[i2]-cenx), dy = wr(py_[i2]-ceny), dz = wr(pz_[i2]-cenz);
                    double rr = sqrt(dx*dx + dy*dy + dz*dz);
                    if (rr > rf) rf = rr;
                }
                printf(" | Ee_front=%.3f", rf);
                if (nfsamp < 512) { fs_t[nfsamp] = sim_t; fs_r[nfsamp] = rf; nfsamp++; }
            }
            if (truss_n > 0) {
                double mx = 0, sm = 0, gmean = 0;
                for (int e = 0; e < truss_n; e++) {
                    int i = truss_e[e][0], j = truss_e[e][1];
                    double dx = px_[j]-px_[i], dy = py_[j]-py_[i], dz = pz_[j]-pz_[i];
                    double d = sqrt(dx*dx+dy*dy+dz*dz);
                    double dev = fabs(d - truss_dstar[e]);
                    if (dev > mx) mx = dev;
                    sm += dev;
                    int s = hfind(i < j ? i : j, i < j ? j : i);
                    if (s >= 0 && sst[s] != S_FREE) {
                        double gf = gate_of(slq[s]*th2[i] - slq[s]*w2e[i]*sd[s]/P.C - slp[s]*th2[j]);
                        double gb = gate_of(slp[s]*th2[j] - slp[s]*w2e[j]*sd[s]/P.C - slq[s]*th2[i]);
                        gmean += gf * gb;
                    }
                }
                double eg[3]; gyration_eigs(eg);
                double sdev = fabs(eg[0]-truss_shape0[0]) + fabs(eg[1]-truss_shape0[1])
                            + fabs(eg[2]-truss_shape0[2]);
                double xm = 0; int nxm = 0;
                for (int i2 = 0; i2 < NC; i2++)
                    if (tag[i2]) { xm += (Em[i2] + flload[i2]) / P.cap; nxm++; }
                printf(" | edge_dev mean=%.4f max=%.4f gg=%.4f shape_dev=%.4f xbar=%.4f",
                       sm / truss_n, mx, gmean / truss_n, sdev, nxm ? xm / nxm : 0);
            }
            if (tri_on) {
                for (int t2 = 0; t2 < ntri; t2++) {
                    double psi[3], ggm, circ;
                    tri_meters(t2, psi, &ggm, &circ);
                    printf(" | T%d psi=(%+.3f,%+.3f,%+.3f) ggm=%.3f circ=%+.4f xUDD=(%.3f,%.3f,%.3f)",
                           t2, psi[0], psi[1], psi[2], ggm, circ,
                           (Em[tri_v[t2][0]] + flload[tri_v[t2][0]]) / P.cap,
                           (Em[tri_v[t2][1]] + flload[tri_v[t2][1]]) / P.cap,
                           (Em[tri_v[t2][2]] + flload[tri_v[t2][2]]) / P.cap);
                }
                if (ntri == 2) {
                    double c1x = 0, c1y = 0, c2x = 0, c2y = 0;
                    for (int kk = 0; kk < 3; kk++) {
                        c1x += cenx + wr(px_[tri_v[0][kk]] - cenx);
                        c1y += ceny + wr(py_[tri_v[0][kk]] - ceny);
                        c2x += cenx + wr(px_[tri_v[1][kk]] - cenx);
                        c2y += ceny + wr(py_[tri_v[1][kk]] - ceny);
                    }
                    double ddx = (c2x - c1x) / 3, ddy = (c2y - c1y) / 3;
                    printf(" | sep=%.4f", sqrt(ddx*ddx + ddy*ddy));
                }
            }
            printf("\n");
        }
        if (P.snap_every > 0 && P.snap_dir[0] && st % P.snap_every == 0)
            write_fcs(st / P.snap_every);
        if (fcs_stream && P.snap_every > 0 && st % P.snap_every == 0) {
            fcs_cell_frame(fcs_stream);
            fcs_instrument(fcs_stream);
        }
        if (st == steps) break;
        step();
        if (!strcmp(P.exp, "slit")) {
            double tn = (st + 1) * P.dt;
            if (tn >= P.slit_t0 && tn <= P.slit_t1) {
                if (!ds_cellI) ds_cellI = calloc(NC, sizeof(double));
                for (int i = 0; i < NC; i++) {
                    if (fabs(wr(px_[i] - P.slit_screenx)) > 0.75) continue;
                    int b = (int)(py_[i] / P.L * DS_NBIN);
                    if (b < 0) b = 0;
                    if (b >= DS_NBIN) b = DS_NBIN - 1;
                    ds_I[b] += Ee[i] * P.dt;
                    ds_cellI[i] += Ee[i] * P.dt;
                    ds_expo += Ee[i] * P.dt;
                }
            }
        }
    }

    /* ---------------- final report ---------------- */
    {
        double Etot = total_energy();
        printf("#\n# RESULT drift_rel %.3e\n", (Etot - E0_total) / (E0_total != 0 ? E0_total : 1));
        printf("# RESULT births %ld deaths %ld beta_returns %ld beta_energy %.6e\n",
               births, deaths, beta_returns, beta_energy);
        printf("# RESULT conv rough=%.6f cond=%.6f evap=%.6f backs=%.6f\n",
               rough_total, cond_total, evap_total, backs_total);
        double phi, zl, dbar, sig, mld; int nla, nld;
        geo_stats(&phi, &zl, &nla, &nld, &dbar, &sig, &mld);
        printf("# RESULT geometry phi_nom=%.4f z_live=%.2f dbar=%.4f sigma_d=%.4f\n",
               phi, zl, dbar, sig);
        if (rep_pair_i >= 0) {
            int i = rep_pair_i, j = rep_pair_j;
            int s = hfind(i < j ? i : j, i < j ? j : i);
            if (s >= 0 && sst[s] != S_FREE) {
                double x1r = P.pair_x1 < 0 ? P.pair_x0 : P.pair_x1;
                double wA = P.w2 / (1.0 + P.q_detune * P.pair_x0);
                double wB = P.w2 / (1.0 + P.q_detune * x1r);
                double ds = TWO_PI * P.pair_m * P.C / (P.pair_qq * wA + P.pair_pp * wB);
                /* the LIVE rung from current pitches — the honest target */
                double dsl = TWO_PI * P.pair_m * P.C / (slq[s] * w2e[i] + slp[s] * w2e[j]);
                printf("# RESULT pair d_final=%.6f d_star_seed=%.6f d_star_live=%.6f off_live=%+.6f Em_i=%.5f Em_j=%.5f\n",
                       sd[s], ds, dsl, sd[s] - dsl, Em[i], Em[j]);
            } else printf("# RESULT pair CHANNEL DEAD\n");
        }
        if (truss_n > 0) {
            double mx = 0, sm = 0;
            for (int e = 0; e < truss_n; e++) {
                int i = truss_e[e][0], j = truss_e[e][1];
                double dx = px_[j]-px_[i], dy = py_[j]-py_[i], dz = pz_[j]-pz_[i];
                double d = sqrt(dx*dx+dy*dy+dz*dz);
                double dev = fabs(d - truss_dstar[e]);
                if (dev > mx) mx = dev;
                sm += dev;
            }
            double eg[3]; gyration_eigs(eg);
            printf("# RESULT truss edge_dev_mean=%.6f edge_dev_max=%.6f shape=(%.4f,%.4f,%.4f) shape0=(%.4f,%.4f,%.4f)\n",
                   sm / truss_n, mx, eg[0], eg[1], eg[2],
                   truss_shape0[0], truss_shape0[1], truss_shape0[2]);
            /* per-edge autopsy: d, dev, gates, and the collimation each
             * edge actually gets (the bond strength carries gpl) */
            for (int e = 0; e < truss_n; e++) {
                int i = truss_e[e][0], j = truss_e[e][1];
                int s = hfind(i < j ? i : j, i < j ? j : i);
                if (s < 0 || sst[s] == S_FREE) {
                    printf("# EDGE %d %d-%d CHANNEL DEAD\n", e, i, j);
                    continue;
                }
                double gf = gate_of(slq[s]*th2[i] - slq[s]*w2e[i]*sd[s]/P.C - slp[s]*th2[j]);
                double gb = gate_of(slp[s]*th2[j] - slp[s]*w2e[j]*sd[s]/P.C - slq[s]*th2[i]);
                double d1i = n1x[i]*sux[s] + n1y[i]*suy[s] + n1z[i]*suz[s];
                double d2i = n2x[i]*sux[s] + n2y[i]*suy[s] + n2z[i]*suz[s];
                double d1j = n1x[j]*sux[s] + n1y[j]*suy[s] + n1z[j]*suz[s];
                double d2j = n2x[j]*sux[s] + n2y[j]*suy[s] + n2z[j]*suz[s];
                double axi = (1.0 - d1i*d1i)*(1.0 - d2i*d2i); if (axi < 0) axi = 0;
                double axj = (1.0 - d1j*d1j)*(1.0 - d2j*d2j); if (axj < 0) axj = 0;
                double dnn = n2x[i]*n2x[j] + n2y[i]*n2y[j] + n2z[i]*n2z[j];
                printf("# EDGE %d %d-%d d=%.5f dev=%+.5f gg=%.4f gpl=%.4f\n",
                       e, i, j, sd[s], sd[s] - truss_dstar[e], gf * gb,
                       axi * axj * dnn * dnn);
            }
        }
        if (tri_on) {
            for (int t2 = 0; t2 < ntri; t2++) {
                double psi[3], ggm, circ;
                tri_meters(t2, psi, &ggm, &circ);
                printf("# RESULT tri%d psi=(%+.4f,%+.4f,%+.4f) ggm=%.4f circ=%+.5f Em=(%.4f,%.4f,%.4f)\n",
                       t2, psi[0], psi[1], psi[2], ggm, circ,
                       Em[tri_v[t2][0]], Em[tri_v[t2][1]], Em[tri_v[t2][2]]);
            }
        }
        if (!strcmp(P.exp, "pulse") && nfsamp >= 6) {
            int h = nfsamp / 2, n = nfsamp - h;
            double sx = 0, sy = 0, sxx = 0, sxy = 0;
            for (int k = h; k < nfsamp; k++)
                { sx += fs_t[k]; sy += fs_r[k]; sxx += fs_t[k]*fs_t[k]; sxy += fs_t[k]*fs_r[k]; }
            double den = n * sxx - sx * sx;
            double v = den != 0 ? (n * sxy - sx * sy) / den : 0;
            printf("# RESULT front_speed v=%.4f v_over_C=%.4f (second-half fit, n=%d)\n",
                   v, v / P.C, n);
        }
        if (!strcmp(P.exp, "blob")) {
            double sh[8]; es_shells(sh);
            printf("# RESULT es_shells %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f\n",
                   sh[0], sh[1], sh[2], sh[3], sh[4], sh[5], sh[6], sh[7]);
            double emt, comx, comy, comz, rms, conn, zb;
            tag_stats(&emt, &comx, &comy, &comz, &rms, &conn, &zb);
            printf("# RESULT blob Em_tag=%.5f com=(%.3f,%.3f,%.3f) rms=%.4f conn=%.4f z_tag=%.2f\n",
                   emt, comx, comy, comz, rms, conn, zb);
            if (nfsamp >= 6) {
                /* e3-class drift: linear fit of the dense centroid over
                 * the second half (cellfab RESULT blob_drift convention) */
                int h = nfsamp / 2, n = nfsamp - h;
                double sx = 0, sxx = 0, syx = 0, syy = 0, syz = 0, sxyx = 0, sxyy = 0, sxyz = 0;
                for (int k = h; k < nfsamp; k++) {
                    sx += fs_t[k]; sxx += fs_t[k]*fs_t[k];
                    syx += fs_r[k]; sxyx += fs_t[k]*fs_r[k];
                    syy += fs_y[k]; sxyy += fs_t[k]*fs_y[k];
                    syz += fs_z[k]; sxyz += fs_t[k]*fs_z[k];
                }
                double den = n * sxx - sx * sx;
                double vx = den != 0 ? (n * sxyx - sx * syx) / den : 0;
                double vy = den != 0 ? (n * sxyy - sx * syy) / den : 0;
                double vz = den != 0 ? (n * sxyz - sx * syz) / den : 0;
                double sp = sqrt(vx*vx + vy*vy + vz*vz);
                printf("# RESULT blob_drift speed=%.6f cos_to_kdir=%.4f v=(%.2e,%.2e,%.2e)\n",
                       sp, sp > 0 ? vx / sp : 0, vx, vy, vz);
            }
        }
        if (fcs_stream) {
            char fin[512];
            snprintf(fin, sizeof fin,
                "RESULT drift_rel %.3e\nRESULT births %ld deaths %ld beta_returns %ld\n"
                "RESULT conv rough=%.6f cond=%.6f evap=%.6f backs=%.6f\n",
                (Etot - E0_total) / (E0_total != 0 ? E0_total : 1),
                births, deaths, beta_returns,
                rough_total, cond_total, evap_total, backs_total);
            fcs_anlz_text(fcs_stream, fin);
            fclose(fcs_stream);
            fcs_stream = NULL;
        }
        if (!strcmp(P.exp, "slit")) {
            printf("# RESULT ds exposure=%.6f gate=[%g,%g] screen_x=%g nbin=%d\n",
                   ds_expo, P.slit_t0, P.slit_t1, P.slit_screenx, DS_NBIN);
            for (int b = 0; b < DS_NBIN; b++)
                printf("# SCREEN y=%.4f I=%.8f\n",
                       (b + 0.5) * P.L / DS_NBIN, ds_I[b]);
            /* per-cell record: the analyzer smooths these (foam speckle
             * is per-cell; the r(y) ratio cancels what smoothing misses) */
            if (ds_cellI)
                for (int i = 0; i < NC; i++)
                    if (ds_cellI[i] > 0)
                        printf("# SCREENCELL y=%.4f x=%.4f I=%.8f\n",
                               py_[i], px_[i], ds_cellI[i]);
        }
    }
    return 0;
}
