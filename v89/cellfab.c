/*
 * cellfab.c — v89 cell-fabric kernel. Written fresh for the v89 program;
 * shares no code, no lattice, no format with any pre-v89 instrument.
 *
 * Ontology (see CELLFAB.md, subordinate to PRINCIPLE.md):
 *   - The field is a set of cells; each cell IS a parcel of energy. Its
 *     space ledger E_s is the cell being space. No grid: cells have areas
 *     of effect and interact only through a list of overlap channels.
 *   - Each cell aligns internally on two planes (normals n1, n2) carrying
 *     two different harmonic frequencies (w1 fast: field sector, w2 slow:
 *     dense sector). Their beat is the conversion clock between sectors.
 *   - Energy transfers between cells only by resonant joining: geometric
 *     conductance (contact area / distance) x plane overlap x frequency
 *     resonance x tail-phase match, cycle-gated with latency d/C. Nothing
 *     moves; clocks and planes re-align and energy converts.
 *   - Saturation (capacity + load detuning) is why everything does not
 *     collapse into one cell: the harmonics refuse the energy, and energy
 *     is never destroyed. Conservation is THE law: every operation is a
 *     paired ledger move; total drift is reported every diagnostic row.
 *
 * Modes:
 *   mode=field  cell-fabric evolution (init = vacuum|noise|pulse|blob)
 *   mode=bell   event-level CHSH on the kernel's transfer law: joint
 *               unfinished harmonic (pair = one process) vs LHV control
 *
 * Output: '#' headers + TSV diagnostics on stdout, '# RESULT' lines at end,
 * optional per-cell TSV snapshots. No SFA (lattice instrument, not v89).
 *
 * Build:  gcc -O2 -o v89/cellfab v89/cellfab.c -lm
 * Run:    ./v89/cellfab [config.cfg] [key=value ...]
 */

#define _GNU_SOURCE   /* sincos */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <sys/stat.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#define TWO_PI (2.0 * M_PI)

/* ------------------------------------------------------------------ */
/* config                                                              */
/* ------------------------------------------------------------------ */

typedef struct {
    unsigned long seed;
    int mode;        /* 0 field, 1 bell, 2 ladder (analytic) */
    int init;        /* 0 vacuum, 1 noise, 2 pulse, 3 blob, 4 pairs */
    double L, dmin, r0, rjit;
    double C;                    /* conversion rate */
    double w1, w2;               /* the two plane frequencies */
    double q_detune;             /* load detuning strength */
    double gamma_res;            /* resonance acceptance width (field sector) */
    double gamma_res_m;          /* resonance width, dense sector (<0: same) */
    double p_gate;               /* tail-gate sharpness exponent */
    double lock_floor;           /* injection-lock inertia floor (x cap) */
    double k_dep, k_dep_m;       /* deposit rate scale; dense multiplier */
    double cap;                  /* internal bound on E_m + E_e */
    double e_s0, es_floor;       /* space ledger nominal / conversion floor */
    double e_cond, f_conv, f_evap, s_pull;
    double kappa_lock, kappa_align, sigma_tumble;
    double kappa_freq;           /* A1: dispersive exchange bias (frequency entrainment) */
    double kappa_reac;           /* S2: derived reactive exchange (choir's correction, no posit) */
    double s_k, s_disp;          /* G: space-mode transport (pressure pushes) + displacement */
    double es_gx;                /* apparatus: imposed initial space gradient along x */
    double dt, T;
    int diag_every, snap_every;
    char snap_dir[256];
    double noise_amp;
    double px, py, pz;           /* seed center (<0 -> box center) */
    double amp, sigma;           /* seed amplitude / width */
    double kx, ky, kz;           /* phase tilt (the ghost of momentum) */
    double aux_amp, aux_x, aux_y, aux_z;   /* auxiliary field pulse (apparatus) */
    double aux_kx, aux_ky, aux_kz, aux_sigma;
    int prealign;                /* align plane normals transverse at init */
    int rad_diag;                /* G4 apparatus: per-shell radial flux report */
    int lump_diag;               /* MASS apparatus: connected-lump census rows */
    double lump_thr;             /* census membership threshold on free Em */
    int ring_n, ring_wind;       /* MASS apparatus: ring seed (voices, winding) */
    int ring_m;                  /* target closure integer (<0 = per-link auto) */
    int shell_n;                 /* H: closed-shell seed (8=cube; 0=pocket only) */
    double shell_a, shell_x;     /* shell edge target; load (<0 = pi-rung auto) */
    double es_core, r_core;      /* H: the hole — core space store; core radius */
    double ring_d, ring_x;       /* ring spacing; load (<0 = tuning-curve auto) */
    long trials;                 /* bell trials per angle combo */
    double aA1, aA2, aB1, aB2;   /* analyzer angles, degrees */
    int debug;                   /* extra gate/flux statistics per diag row */
    int npairs;                  /* standing-pair seeds (init=pairs) */
    double pair_x0;              /* pair occupancy; <0: each at its x*(d) */
    double pair_x1;              /* second cell occupancy; <0: same as x0 */
    int pair_p, pair_q;          /* seeded interval p:q (1:1 unison) */
    int pair_seedlock;           /* 1: seed phases at lock; 0: random */
    double pair_dlo, pair_dhi;   /* link-length window for pair selection */
    double pair_gap;             /* min distance between distinct pairs */
    /* round 2 — music theory in the kernel */
    int comb_limit;              /* C1: max p*q in the partial comb (1 = legacy 1:1) */
    double rough_k;              /* C2: fraction of roughness radiated D->F (0 = legacy) */
    double gamma_rough;          /* C2: Plomp-Levelt roughness peak detuning */
    int mob_sym;                 /* C3: 1 = mutual (sqrt) dense coupling, 0 = legacy */
    double mob_floor;            /* C3: sympathetic readiness floor (x cap) */
    /* double slit (init=slit): wall as detuned medium, slits as vacuum
     * windows, screen/sinks as condensation records, which-path recorders */
    double wall_x, wall_th;      /* wall plane center / thickness */
    double slit_sep, slit_hw;    /* slit strip separation / half-width (y) */
    double screen_x, sink_m;     /* screen plane; edge sink margin */
    double edge_sink;            /* any-init absorbing edge margin (0=off) */
    double src_lo, src_hi;       /* source slab x-range */
    double absorb_rate;          /* record medium condensation rate */
    double wall_detune;          /* wall field-plane pitch offset */
    int which_path;              /* 0 none, 1 both slits, 2 slit A only */
    double det_tap, sigma_det;   /* recorder tap rate; recorder clock walk */
    double e_click;              /* record grain size for click logging (0 off) */
    int block_slit;              /* 0 none, 1 block A, 2 block B */
    double t_expose;             /* shutter: freeze the screen record at this time (0 off) */
    double field_J;              /* repaired field sector: hop coupling scale */
    int n_quanta;                /* tier 1: sample N single-quantum claims (0 off) */
    double tag_rate;             /* tier 2: unitary which-path tag rotation at slit A */
    double analyzer_deg;         /* tier 2: screen analyzer basis (deg; <-900 = native H/V) */
    double t_choice;             /* tier 3: flip analyzer basis at this time (0 = fixed) */
    /* tier 4 (init=hom): two waveguides + tunneling coupler; two quanta as
     * the two field components; coincidence from the exchange-symmetrized
     * joint amplitude (exact for non-interacting quanta) */
    double hom_c0, hom_c1;       /* coupler window in x */
    double hom_barrier;          /* coupler gap detune (smaller = stronger coupling) */
    double hom_dx;               /* launch offset of quantum B (the delay knob) */
    int hom_solo;                /* 1: launch only quantum A (coupler calibration) */
    double hom_det;              /* detector window length before the end sink */
    /* F-h route (HBAR.md candidate 1+3): the atom of transfer is action.
     * Dense deposits and beat conversions move in whole quanta
     * eps(w) = A0*w/2pi. 0 = off (legacy float-identical); -1 = auto:
     * A0 = e_s0 * dbar / C — the space-cell grain, no new constant. */
    double quant_A0;
    int quant_mode;              /* U2 variant: 1 src+floor, 2 src+credit, 3 dst+floor, 4 dst+credit */
} Cfg;

static Cfg P;

static void cfg_defaults(void)
{
    memset(&P, 0, sizeof(P));
    P.seed = 12345;
    P.mode = 0;
    P.init = 0;
    P.L = 24.0; P.dmin = 1.0; P.r0 = 0.85; P.rjit = 0.06;
    P.C = 1.0;
    /* w*d/C should sit near pi/2 for typical link length d: the backward
     * tail gate wraps to ~pi (sealed). Near w*d ~ pi it wraps to ~0 and
     * reopens — a foam-geometry resonance to keep clear of. */
    P.w1 = 1.5; P.w2 = 0.93;
    P.q_detune = 0.35;
    P.gamma_res = 0.25;
    P.gamma_res_m = -1.0;
    P.p_gate = 4.0;
    P.lock_floor = 0.005;
    P.k_dep = 1.2; P.k_dep_m = 1.0;
    P.cap = 2.5;
    P.e_s0 = 1.0; P.es_floor = 0.05;
    P.e_cond = 0.30; P.f_conv = 0.25; P.f_evap = 0.5; P.s_pull = 0.5;
    P.kappa_lock = 0.9; P.kappa_align = 0.5; P.sigma_tumble = 0.01;
    P.kappa_freq = 0.0;
    P.kappa_reac = 0.0;
    P.s_k = 0.0; P.s_disp = 0.0; P.es_gx = 0.0;
    P.dt = 0.02; P.T = 40.0;
    P.diag_every = 50; P.snap_every = 0;
    strcpy(P.snap_dir, "snaps");
    P.noise_amp = 0.0;
    P.px = -1; P.py = -1; P.pz = -1;
    P.amp = 1.0; P.sigma = 2.5;
    P.kx = 0; P.ky = 0; P.kz = 0;
    P.aux_amp = 0; P.aux_x = -1; P.aux_y = -1; P.aux_z = -1;
    P.aux_kx = 0; P.aux_ky = 0; P.aux_kz = 0; P.aux_sigma = -1;
    P.prealign = 0;
    P.rad_diag = 0;
    P.lump_diag = 0; P.lump_thr = 0.1;
    P.ring_n = 0; P.ring_wind = 0; P.ring_d = 1.25; P.ring_x = -1;
    P.ring_m = -1;
    P.shell_n = 8; P.shell_a = 1.5; P.shell_x = -1;
    P.es_core = -1; P.r_core = 0.8;
    P.trials = 200000;
    P.aA1 = 0.0; P.aA2 = 45.0; P.aB1 = 22.5; P.aB2 = 67.5;
    P.npairs = 48;
    P.pair_x0 = 0.4;
    P.pair_x1 = -1.0;
    P.pair_p = 1; P.pair_q = 1;
    P.pair_seedlock = 1;
    P.pair_dlo = 1.0; P.pair_dhi = 1.55;
    P.pair_gap = 5.0;
    P.comb_limit = 6;
    P.rough_k = 0.35;
    P.gamma_rough = 0.5;
    P.mob_sym = 1;
    P.mob_floor = 0.01;
    P.wall_x = 11.5; P.wall_th = 2.4;
    P.slit_sep = 9.0; P.slit_hw = 1.2;
    P.screen_x = 26.5; P.sink_m = 1.5;
    P.edge_sink = 0;
    P.src_lo = 1.5; P.src_hi = 8.5;
    P.absorb_rate = 2.0;
    P.wall_detune = 8.0;
    P.which_path = 0;
    P.det_tap = 0.06; P.sigma_det = 6.0;
    P.e_click = 0.0;
    P.block_slit = 0;
    P.t_expose = 0.0;
    P.field_J = 0.06;
    P.n_quanta = 0;
    P.tag_rate = 0.0;
    P.analyzer_deg = -999.0;
    P.t_choice = 0.0;
    P.hom_c0 = 15.0; P.hom_c1 = 21.0;
    P.hom_barrier = 1.2;
    P.hom_dx = 0.0;
    P.hom_solo = 0;
    P.hom_det = 7.0;
    P.quant_A0 = 0.0;
    P.quant_mode = 1;
}

static void set_kv(const char *k, const char *v)
{
    if (!strcmp(k, "seed")) P.seed = strtoul(v, NULL, 10);
    else if (!strcmp(k, "mode")) P.mode = !strcmp(v, "bell") ? 1 : (!strcmp(v, "ladder") ? 2 : 0);
    else if (!strcmp(k, "init")) {
        if (!strcmp(v, "vacuum")) P.init = 0;
        else if (!strcmp(v, "noise")) P.init = 1;
        else if (!strcmp(v, "pulse")) P.init = 2;
        else if (!strcmp(v, "blob")) P.init = 3;
        else if (!strcmp(v, "pairs")) P.init = 4;
        else if (!strcmp(v, "slit")) P.init = 5;
        else if (!strcmp(v, "hom")) P.init = 6;
        else if (!strcmp(v, "ring")) P.init = 7;
        else if (!strcmp(v, "shell")) P.init = 8;
        else fprintf(stderr, "# WARN unknown init '%s'\n", v);
    }
    else if (!strcmp(k, "L")) P.L = atof(v);
    else if (!strcmp(k, "dmin")) P.dmin = atof(v);
    else if (!strcmp(k, "r0")) P.r0 = atof(v);
    else if (!strcmp(k, "rjit")) P.rjit = atof(v);
    else if (!strcmp(k, "C")) P.C = atof(v);
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
    else if (!strcmp(k, "es_gx")) P.es_gx = atof(v);
    else if (!strcmp(k, "sigma_tumble")) P.sigma_tumble = atof(v);
    else if (!strcmp(k, "dt")) P.dt = atof(v);
    else if (!strcmp(k, "T")) P.T = atof(v);
    else if (!strcmp(k, "diag_every")) P.diag_every = atoi(v);
    else if (!strcmp(k, "snap_every")) P.snap_every = atoi(v);
    else if (!strcmp(k, "snap_dir")) { strncpy(P.snap_dir, v, 255); P.snap_dir[255] = 0; }
    else if (!strcmp(k, "noise_amp")) P.noise_amp = atof(v);
    else if (!strcmp(k, "px")) P.px = atof(v);
    else if (!strcmp(k, "py")) P.py = atof(v);
    else if (!strcmp(k, "pz")) P.pz = atof(v);
    else if (!strcmp(k, "amp")) P.amp = atof(v);
    else if (!strcmp(k, "sigma")) P.sigma = atof(v);
    else if (!strcmp(k, "kx")) P.kx = atof(v);
    else if (!strcmp(k, "ky")) P.ky = atof(v);
    else if (!strcmp(k, "kz")) P.kz = atof(v);
    else if (!strcmp(k, "aux_amp")) P.aux_amp = atof(v);
    else if (!strcmp(k, "aux_x")) P.aux_x = atof(v);
    else if (!strcmp(k, "aux_y")) P.aux_y = atof(v);
    else if (!strcmp(k, "aux_z")) P.aux_z = atof(v);
    else if (!strcmp(k, "aux_kx")) P.aux_kx = atof(v);
    else if (!strcmp(k, "aux_ky")) P.aux_ky = atof(v);
    else if (!strcmp(k, "aux_kz")) P.aux_kz = atof(v);
    else if (!strcmp(k, "aux_sigma")) P.aux_sigma = atof(v);
    else if (!strcmp(k, "prealign")) P.prealign = atoi(v);
    else if (!strcmp(k, "rad_diag")) P.rad_diag = atoi(v);
    else if (!strcmp(k, "lump_diag")) P.lump_diag = atoi(v);
    else if (!strcmp(k, "lump_thr")) P.lump_thr = atof(v);
    else if (!strcmp(k, "ring_n")) P.ring_n = atoi(v);
    else if (!strcmp(k, "ring_wind")) P.ring_wind = atoi(v);
    else if (!strcmp(k, "ring_m")) P.ring_m = atoi(v);
    else if (!strcmp(k, "shell_n")) P.shell_n = atoi(v);
    else if (!strcmp(k, "shell_a")) P.shell_a = atof(v);
    else if (!strcmp(k, "shell_x")) P.shell_x = atof(v);
    else if (!strcmp(k, "es_core")) P.es_core = atof(v);
    else if (!strcmp(k, "r_core")) P.r_core = atof(v);
    else if (!strcmp(k, "ring_d")) P.ring_d = atof(v);
    else if (!strcmp(k, "ring_x")) P.ring_x = atof(v);
    else if (!strcmp(k, "trials")) P.trials = atol(v);
    else if (!strcmp(k, "aA1")) P.aA1 = atof(v);
    else if (!strcmp(k, "aA2")) P.aA2 = atof(v);
    else if (!strcmp(k, "aB1")) P.aB1 = atof(v);
    else if (!strcmp(k, "aB2")) P.aB2 = atof(v);
    else if (!strcmp(k, "debug")) P.debug = atoi(v);
    else if (!strcmp(k, "npairs")) P.npairs = atoi(v);
    else if (!strcmp(k, "pair_x0")) P.pair_x0 = atof(v);
    else if (!strcmp(k, "pair_x1")) P.pair_x1 = atof(v);
    else if (!strcmp(k, "pair_p")) P.pair_p = atoi(v);
    else if (!strcmp(k, "pair_q")) P.pair_q = atoi(v);
    else if (!strcmp(k, "pair_seedlock")) P.pair_seedlock = atoi(v);
    else if (!strcmp(k, "pair_dlo")) P.pair_dlo = atof(v);
    else if (!strcmp(k, "pair_dhi")) P.pair_dhi = atof(v);
    else if (!strcmp(k, "pair_gap")) P.pair_gap = atof(v);
    else if (!strcmp(k, "comb_limit")) P.comb_limit = atoi(v);
    else if (!strcmp(k, "rough_k")) P.rough_k = atof(v);
    else if (!strcmp(k, "gamma_rough")) P.gamma_rough = atof(v);
    else if (!strcmp(k, "mob_sym")) P.mob_sym = atoi(v);
    else if (!strcmp(k, "mob_floor")) P.mob_floor = atof(v);
    else if (!strcmp(k, "wall_x")) P.wall_x = atof(v);
    else if (!strcmp(k, "wall_th")) P.wall_th = atof(v);
    else if (!strcmp(k, "slit_sep")) P.slit_sep = atof(v);
    else if (!strcmp(k, "slit_hw")) P.slit_hw = atof(v);
    else if (!strcmp(k, "screen_x")) P.screen_x = atof(v);
    else if (!strcmp(k, "sink_m")) P.sink_m = atof(v);
    else if (!strcmp(k, "edge_sink")) P.edge_sink = atof(v);
    else if (!strcmp(k, "src_lo")) P.src_lo = atof(v);
    else if (!strcmp(k, "src_hi")) P.src_hi = atof(v);
    else if (!strcmp(k, "absorb_rate")) P.absorb_rate = atof(v);
    else if (!strcmp(k, "wall_detune")) P.wall_detune = atof(v);
    else if (!strcmp(k, "which_path")) P.which_path = atoi(v);
    else if (!strcmp(k, "det_tap")) P.det_tap = atof(v);
    else if (!strcmp(k, "sigma_det")) P.sigma_det = atof(v);
    else if (!strcmp(k, "e_click")) P.e_click = atof(v);
    else if (!strcmp(k, "block_slit")) P.block_slit = atoi(v);
    else if (!strcmp(k, "t_expose")) P.t_expose = atof(v);
    else if (!strcmp(k, "field_J")) P.field_J = atof(v);
    else if (!strcmp(k, "n_quanta")) P.n_quanta = atoi(v);
    else if (!strcmp(k, "tag_rate")) P.tag_rate = atof(v);
    else if (!strcmp(k, "analyzer_deg")) P.analyzer_deg = atof(v);
    else if (!strcmp(k, "t_choice")) P.t_choice = atof(v);
    else if (!strcmp(k, "hom_c0")) P.hom_c0 = atof(v);
    else if (!strcmp(k, "hom_c1")) P.hom_c1 = atof(v);
    else if (!strcmp(k, "hom_barrier")) P.hom_barrier = atof(v);
    else if (!strcmp(k, "hom_dx")) P.hom_dx = atof(v);
    else if (!strcmp(k, "hom_solo")) P.hom_solo = atoi(v);
    else if (!strcmp(k, "hom_det")) P.hom_det = atof(v);
    else if (!strcmp(k, "quant_A0")) P.quant_A0 = atof(v);
    else if (!strcmp(k, "quant_mode")) P.quant_mode = atoi(v);
    else fprintf(stderr, "# WARN unknown key '%s'\n", k);
}

static void load_cfg(const char *path)
{
    FILE *f = fopen(path, "r");
    if (!f) { fprintf(stderr, "# ERROR cannot open cfg '%s'\n", path); exit(1); }
    char line[512];
    while (fgets(line, sizeof(line), f)) {
        char *h = strchr(line, '#');
        if (h) *h = 0;
        char *eq = strchr(line, '=');
        if (!eq) continue;
        *eq = 0;
        char *k = line, *v = eq + 1;
        while (*k == ' ' || *k == '\t') k++;
        char *e = k + strlen(k);
        while (e > k && (e[-1] == ' ' || e[-1] == '\t' || e[-1] == '\n' || e[-1] == '\r')) *--e = 0;
        while (*v == ' ' || *v == '\t') v++;
        e = v + strlen(v);
        while (e > v && (e[-1] == ' ' || e[-1] == '\t' || e[-1] == '\n' || e[-1] == '\r')) *--e = 0;
        if (*k) set_kv(k, v);
    }
    fclose(f);
}

static void print_cfg(void)
{
    const char *modes[] = { "field", "bell", "ladder" };
    const char *inits[] = { "vacuum", "noise", "pulse", "blob", "pairs", "slit", "hom", "ring", "shell" };
    printf("# cellfab — v89 cell-fabric kernel (no lattice, no prior code)\n");
    printf("# cfg seed=%lu mode=%s init=%s\n", P.seed, modes[P.mode], inits[P.init]);
    printf("# cfg L=%g dmin=%g r0=%g rjit=%g\n", P.L, P.dmin, P.r0, P.rjit);
    printf("# cfg C=%g w1=%g w2=%g q_detune=%g gamma_res=%g gamma_res_m=%g p_gate=%g\n",
           P.C, P.w1, P.w2, P.q_detune, P.gamma_res,
           P.gamma_res_m < 0 ? P.gamma_res : P.gamma_res_m, P.p_gate);
    printf("# cfg lock_floor=%g\n", P.lock_floor);
    printf("# cfg k_dep=%g k_dep_m=%g cap=%g e_s0=%g es_floor=%g\n",
           P.k_dep, P.k_dep_m, P.cap, P.e_s0, P.es_floor);
    printf("# cfg e_cond=%g f_conv=%g f_evap=%g s_pull=%g\n",
           P.e_cond, P.f_conv, P.f_evap, P.s_pull);
    printf("# cfg kappa_lock=%g kappa_align=%g kappa_freq=%g kappa_reac=%g sigma_tumble=%g\n",
           P.kappa_lock, P.kappa_align, P.kappa_freq, P.kappa_reac, P.sigma_tumble);
    printf("# cfg space: s_k=%g s_disp=%g es_gx=%g\n", P.s_k, P.s_disp, P.es_gx);
    printf("# cfg dt=%g T=%g diag_every=%d snap_every=%d\n",
           P.dt, P.T, P.diag_every, P.snap_every);
    printf("# cfg center=(%g,%g,%g) amp=%g sigma=%g k=(%g,%g,%g) prealign=%d noise_amp=%g\n",
           P.px, P.py, P.pz, P.amp, P.sigma, P.kx, P.ky, P.kz, P.prealign, P.noise_amp);
    if (P.aux_amp > 0)
        printf("# cfg aux pulse: amp=%g at (%g,%g,%g) k=(%g,%g,%g) sigma=%g\n",
               P.aux_amp, P.aux_x, P.aux_y, P.aux_z,
               P.aux_kx, P.aux_ky, P.aux_kz, P.aux_sigma);
    printf("# cfg music: comb_limit=%d rough_k=%g gamma_rough=%g mob_sym=%d mob_floor=%g\n",
           P.comb_limit, P.rough_k, P.gamma_rough, P.mob_sym, P.mob_floor);
    printf("# cfg field sector: two-component signed amplitude, unitary hops, field_J=%g\n",
           P.field_J);
    if (P.init == 4)
        printf("# cfg pairs: n=%d x0=%g x1=%g pq=%d:%d seedlock=%d d=[%g,%g] gap=%g\n",
               P.npairs, P.pair_x0, P.pair_x1, P.pair_p, P.pair_q,
               P.pair_seedlock, P.pair_dlo, P.pair_dhi, P.pair_gap);
    if (P.init == 5)
        printf("# cfg slit: wall_x=%g th=%g sep=%g hw=%g screen_x=%g src=[%g,%g] absorb=%g detune=%g which_path=%d tap=%g sigma_det=%g block=%d\n",
               P.wall_x, P.wall_th, P.slit_sep, P.slit_hw, P.screen_x,
               P.src_lo, P.src_hi, P.absorb_rate, P.wall_detune,
               P.which_path, P.det_tap, P.sigma_det, P.block_slit);
    if (P.n_quanta > 0 || P.tag_rate > 0)
        printf("# cfg quantum layer: n_quanta=%d tag_rate=%g analyzer_deg=%g t_choice=%g\n",
               P.n_quanta, P.tag_rate, P.analyzer_deg, P.t_choice);
    if (P.init == 6)
        printf("# cfg hom: coupler x=[%g,%g] barrier=%g dx=%g solo=%d det=%g\n",
               P.hom_c0, P.hom_c1, P.hom_barrier, P.hom_dx, P.hom_solo, P.hom_det);
}

/* C1: the partial comb — coprime ratios p:q up to complexity comb_limit,
 * ordered simplest first so equal resonance scores resolve to the simpler
 * interval (the ear's rule). */
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
/* rng + small math                                                    */
/* ------------------------------------------------------------------ */

static uint64_t rng_s;

static uint64_t xrand(void)
{
    uint64_t x = rng_s;
    x ^= x << 13; x ^= x >> 7; x ^= x << 17;
    return rng_s = x;
}

static double frand(void) { return (double)(xrand() >> 11) * (1.0 / 9007199254740992.0); }

static double grand(void)
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

static double gate_of(double dphi)
{
    double g = 0.5 * (1.0 + cos(dphi));
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
    do { a = grand(); b = grand(); c = grand(); n2 = a * a + b * b + c * c; } while (n2 < 1e-12);
    double inv = 1.0 / sqrt(n2);
    *x = a * inv; *y = b * inv; *z = c * inv;
}

/* ------------------------------------------------------------------ */
/* state                                                               */
/* ------------------------------------------------------------------ */

static int NC = 0, NL = 0;
/* cells */
static double *cx, *cy, *cz;         /* scaffold positions (init + diagnostics only) */
static double *cr0, *cr;             /* nominal / live area-of-effect radius */
static double *n1x, *n1y, *n1z, *n2x, *n2y, *n2z;
static double *th1, *th2, *cbeta;
static double *w1e, *w2e;
static double *Es, *Em, *Ee;         /* Ee is now a derived cache: |psi|^2 */
static double *fa1, *fa2;            /* the repair: field amplitude on the plane pair (H) */
static double *fb1, *fb2;            /* tier 2: second chirality/polarization component (V) */
static int pol_on = 0;               /* V component active (tag_rate > 0) */
static double *exA, *exB;            /* screen exposure, analyzer-basis resolved */
static double *req0, *req1, *scl0, *scl1;
/* links */
static int *li, *lj;
static double *ld, *lux, *luy, *luz;
static double *lA, *lA0;
static double *lem, *lph, *lwant;    /* [NL][chan][dir] flattened */
static double *lflux;                /* [NL][chan] deposits this step   */
static signed char *lp, *lq;         /* dense sector: locked-in partial ratio p:q */
static double rough_total = 0.0;     /* C2 ledger: dissonance radiated D->F */
static double cond_total = 0.0;      /* M0 ledger: condensed F->D (pass 6) */
static double evap_total = 0.0;      /* M0 ledger: over-full evaporation D->F */
static double backs_total = 0.0;     /* M0 ledger: space returned with D->F */
static double A0eff = 0.0;           /* action atom (0 = continuous limit) */
static double *qcnvD = NULL;         /* per-cell F->D conversion credit (variant 2) */
static double *qcnvF = NULL;         /* per-cell D->F conversion credit (variant 2) */
static long qfire_n = 0;             /* conversion fire counter (QATOM diag) */
/* instrument flags (init=slit): 0 vacuum, 1 wall, 2 screen, 3 sink, 4 recorder */
static signed char *cflag;
static double *fsum;                 /* per-cell total geometric joining weight */
static double *flload;               /* per-cell pitch share of in-flight dense energy */
static double *sprq, *swl;           /* space transport: per-cell wants, per-link flux */
static double *sscl;                 /* space transport: per-cell outflow scale */
static double *th1s, *th2s;          /* clock snapshot: sources read pre-pass values */
static double *nsnap;                /* normals snapshot for the parallel alignment */
static int want_th1 = 1;             /* refresh the diagnostic th1 cache this step */
static double *roughq;               /* per-cell rough fires this step (QATOM flush) */
static double *rngbuf;               /* serial gaussian draws for the parallel tumble */
static int *lcol, *colstart, *colidx, ncolors;  /* edge coloring (pass F hops) */
static double srad[8];               /* radial space flux accumulator (G4 shells) */
static double srad_t0 = 0;           /* window start for the flux rate */
static int *clickn;                  /* record grains already logged per cell */
static double sim_t = 0.0;
static double *expose_frz = NULL;    /* shutter: frozen screen record */

/* tier 1 — claim rule: the absorption process, recorded per screen cell
 * per interval. For non-interacting quanta in the linear regime, sampling
 * first-claims from this process is EXACT: each quantum follows the same
 * unitary evolution, the hazard is the absorb law h_i = rate*|psi_i|^2,
 * and the wave's own norm depletion implements the survival factor. */
static int nscr = 0;
static int *scr_id = NULL;
static double *scr_prev = NULL;
static float *qinc = NULL;
static double *qtms = NULL;
static int qints = 0;
#define QINT_MAX 600
/* cell -> incident links (CSR) */
static int *cls, *clidx;

#define SLOT(l, c, d) (4 * (l) + 2 * (c) + (d))

static double E0_total = 0.0;
static double cenx, ceny, cenz;      /* seed center for diagnostics */

/* diagnostic sample store */
static int nsamp = 0, nsamp_max = 0;
static double *ds_t, *ds_em, *ds_front, *ds_cmx, *ds_cmy, *ds_cmz, *ds_defA;

/* standing-pair registry (init=pairs): the consonance experiments */
#define NP_MAX 256
static int NP = 0;
static int ppi[NP_MAX], ppj[NP_MAX], ppl[NP_MAX];
static double pE0[NP_MAX];

/* ------------------------------------------------------------------ */
/* init: dart-throw cells, overlap channels, seed energy               */
/* ------------------------------------------------------------------ */

static void build_field(void)
{
    /* --- blue-noise cell scaffold: dart throwing, min separation dmin --- */
    double V = P.L * P.L * P.L;
    int nmax = (int)(0.8 * V / (P.dmin * P.dmin * P.dmin)) + 4096;
    cx = malloc(nmax * sizeof(double));
    cy = malloc(nmax * sizeof(double));
    cz = malloc(nmax * sizeof(double));
    int gn = (int)ceil(P.L / P.dmin);
    if (gn < 1) gn = 1;
    int nbins = gn * gn * gn;
    int *head = malloc(nbins * sizeof(int));
    int *nxt = malloc(nmax * sizeof(int));
    for (int b = 0; b < nbins; b++) head[b] = -1;

    int n = 0, fails = 0;
    double d2min = P.dmin * P.dmin;
    while (fails < 30000 && n < nmax) {
        double x = frand() * P.L, y = frand() * P.L, z = frand() * P.L;
        int bx = (int)(x / P.L * gn), by = (int)(y / P.L * gn), bz = (int)(z / P.L * gn);
        if (bx >= gn) bx = gn - 1;
        if (by >= gn) by = gn - 1;
        if (bz >= gn) bz = gn - 1;
        int ok = 1;
        for (int ax = bx - 1; ax <= bx + 1 && ok; ax++)
            for (int ay = by - 1; ay <= by + 1 && ok; ay++)
                for (int az = bz - 1; az <= bz + 1 && ok; az++) {
                    if (ax < 0 || ay < 0 || az < 0 || ax >= gn || ay >= gn || az >= gn) continue;
                    for (int q = head[(ax * gn + ay) * gn + az]; q >= 0; q = nxt[q]) {
                        double dx = cx[q] - x, dy = cy[q] - y, dz = cz[q] - z;
                        if (dx * dx + dy * dy + dz * dz < d2min) { ok = 0; break; }
                    }
                }
        if (!ok) { fails++; continue; }
        cx[n] = x; cy[n] = y; cz[n] = z;
        int b = (bx * gn + by) * gn + bz;
        nxt[n] = head[b]; head[b] = n;
        n++; fails = 0;
    }
    NC = n;
    free(head); free(nxt);

    /* --- cell state --- */
    cr0 = malloc(NC * sizeof(double));
    cr  = malloc(NC * sizeof(double));
    n1x = malloc(NC * sizeof(double)); n1y = malloc(NC * sizeof(double)); n1z = malloc(NC * sizeof(double));
    n2x = malloc(NC * sizeof(double)); n2y = malloc(NC * sizeof(double)); n2z = malloc(NC * sizeof(double));
    th1 = malloc(NC * sizeof(double)); th2 = malloc(NC * sizeof(double)); cbeta = malloc(NC * sizeof(double));
    w1e = malloc(NC * sizeof(double)); w2e = malloc(NC * sizeof(double));
    Es = malloc(NC * sizeof(double)); Em = calloc(NC, sizeof(double)); Ee = calloc(NC, sizeof(double));
    fa1 = calloc(NC, sizeof(double)); fa2 = calloc(NC, sizeof(double));
    fb1 = calloc(NC, sizeof(double)); fb2 = calloc(NC, sizeof(double));
    exA = calloc(NC, sizeof(double)); exB = calloc(NC, sizeof(double));
    pol_on = P.tag_rate > 0 || P.init == 6;
    req0 = malloc(NC * sizeof(double)); req1 = malloc(NC * sizeof(double));
    scl0 = malloc(NC * sizeof(double)); scl1 = malloc(NC * sizeof(double));

    for (int i = 0; i < NC; i++) {
        cr0[i] = P.r0 * (1.0 + P.rjit * (2.0 * frand() - 1.0));
        cr[i] = cr0[i];
        rand_unit(&n1x[i], &n1y[i], &n1z[i]);
        rand_unit(&n2x[i], &n2y[i], &n2z[i]);
        th1[i] = frand() * TWO_PI;
        th2[i] = frand() * TWO_PI;
        cbeta[i] = frand() * TWO_PI;
        Es[i] = P.e_s0;
        /* apparatus: imposed space gradient (G-battery initial condition) */
        if (P.es_gx != 0) {
            Es[i] = P.e_s0 * (1.0 + P.es_gx * (cx[i] - 0.5 * P.L) / (0.5 * P.L));
            double flo = P.es_floor + 0.02;
            if (Es[i] < flo) Es[i] = flo;
        }
    }

    /* --- candidate channels: areas of effect overlap (1.15 margin) --- */
    double rmax = 0;
    for (int i = 0; i < NC; i++) if (cr0[i] > rmax) rmax = cr0[i];
    double hcut = 1.15 * 2.0 * rmax;
    int g2 = (int)ceil(P.L / hcut);
    if (g2 < 1) g2 = 1;
    int nb2 = g2 * g2 * g2;
    int *head2 = malloc(nb2 * sizeof(int));
    int *nxt2 = malloc(NC * sizeof(int));
    for (int b = 0; b < nb2; b++) head2[b] = -1;
    for (int i = 0; i < NC; i++) {
        int bx = (int)(cx[i] / P.L * g2), by = (int)(cy[i] / P.L * g2), bz = (int)(cz[i] / P.L * g2);
        if (bx >= g2) bx = g2 - 1;
        if (by >= g2) by = g2 - 1;
        if (bz >= g2) bz = g2 - 1;
        int b = (bx * g2 + by) * g2 + bz;
        nxt2[i] = head2[b]; head2[b] = i;
    }

    for (int pass = 0; pass < 2; pass++) {
        NL = 0;
        for (int i = 0; i < NC; i++) {
            int bx = (int)(cx[i] / P.L * g2), by = (int)(cy[i] / P.L * g2), bz = (int)(cz[i] / P.L * g2);
            if (bx >= g2) bx = g2 - 1;
            if (by >= g2) by = g2 - 1;
            if (bz >= g2) bz = g2 - 1;
            for (int ax = bx - 1; ax <= bx + 1; ax++)
                for (int ay = by - 1; ay <= by + 1; ay++)
                    for (int az = bz - 1; az <= bz + 1; az++) {
                        if (ax < 0 || ay < 0 || az < 0 || ax >= g2 || ay >= g2 || az >= g2) continue;
                        for (int q = head2[(ax * g2 + ay) * g2 + az]; q >= 0; q = nxt2[q]) {
                            if (q <= i) continue;
                            double dx = cx[q] - cx[i], dy = cy[q] - cy[i], dz = cz[q] - cz[i];
                            double d2 = dx * dx + dy * dy + dz * dz;
                            double cut = 1.15 * (cr0[i] + cr0[q]);
                            if (d2 >= cut * cut) continue;
                            if (pass == 1) {
                                double d = sqrt(d2);
                                li[NL] = i; lj[NL] = q;
                                ld[NL] = d;
                                lux[NL] = dx / d; luy[NL] = dy / d; luz[NL] = dz / d;
                            }
                            NL++;
                        }
                    }
        }
        if (pass == 0) {
            li = malloc(NL * sizeof(int)); lj = malloc(NL * sizeof(int));
            ld = malloc(NL * sizeof(double));
            lux = malloc(NL * sizeof(double)); luy = malloc(NL * sizeof(double)); luz = malloc(NL * sizeof(double));
        }
    }
    free(head2); free(nxt2);

    lA = calloc(NL, sizeof(double));
    lA0 = calloc(NL, sizeof(double));
    lem = calloc(4 * NL, sizeof(double));
    lph = calloc(4 * NL, sizeof(double));
    lwant = calloc(4 * NL, sizeof(double));
    lflux = calloc(2 * NL, sizeof(double));
    lp = malloc(NL); lq = malloc(NL);
    memset(lp, 1, NL); memset(lq, 1, NL);
    cflag = calloc(NC, 1);
    clickn = calloc(NC, sizeof(int));
    fsum = calloc(NC, sizeof(double));
    flload = calloc(NC, sizeof(double));
    sprq = calloc(NC, sizeof(double));
    swl = calloc(NL, sizeof(double));
    sscl = calloc(NC, sizeof(double));
    th1s = calloc(NC, sizeof(double));
    th2s = calloc(NC, sizeof(double));
    roughq = calloc(NC, sizeof(double));
    rngbuf = calloc(6 * (size_t)NC, sizeof(double));
    nsnap = calloc(6 * (size_t)NC, sizeof(double));
    qcnvD = calloc(NC, sizeof(double));
    qcnvF = calloc(NC, sizeof(double));

    /* the action atom: A0 from the space-cell grain (no new constant) */
    double dsum = 0;
    for (int l = 0; l < NL; l++) dsum += ld[l];
    double dbar = NL ? dsum / NL : 1.0;
    A0eff = P.quant_A0 < 0 ? P.e_s0 * dbar / P.C : P.quant_A0;
    if (A0eff > 0)
        printf("# action atom: A0=%.4f (dbar=%.4f) mode=%d — eps(w)=A0*w/2pi; atoms at mode boundaries, transport continuous\n",
               A0eff, dbar, P.quant_mode);

    /* CSR incidence */
    cls = calloc(NC + 1, sizeof(int));
    clidx = malloc(2 * NL * sizeof(int));
    for (int l = 0; l < NL; l++) { cls[li[l] + 1]++; cls[lj[l] + 1]++; }
    for (int i = 0; i < NC; i++) cls[i + 1] += cls[i];
    int *fill = malloc(NC * sizeof(int));
    memcpy(fill, cls, NC * sizeof(int));
    for (int l = 0; l < NL; l++) {
        clidx[fill[li[l]]++] = l;
        clidx[fill[lj[l]]++] = l;
    }

    /* edge coloring for the parallel field hops: links in one color
     * share no cell, so a color batch applies its pair rotations
     * conflict-free; batches run in fixed color order — deterministic
     * for any thread count (the operator ordering differs from the old
     * sequential sweep at roundoff level; physics-gated). */
    lcol = malloc(NL * sizeof(int));
    ncolors = 0;
    {
        unsigned char *used = calloc(4096, 1);
        for (int l = 0; l < NL; l++) {
            int hi = 0;
            for (int q = cls[li[l]]; q < cls[li[l] + 1]; q++) {
                int o = clidx[q];
                if (o < l) { used[lcol[o]] = 1; if (lcol[o] > hi) hi = lcol[o]; }
            }
            for (int q = cls[lj[l]]; q < cls[lj[l] + 1]; q++) {
                int o = clidx[q];
                if (o < l) { used[lcol[o]] = 1; if (lcol[o] > hi) hi = lcol[o]; }
            }
            int c = 0;
            while (used[c]) c++;
            lcol[l] = c;
            if (c + 1 > ncolors) ncolors = c + 1;
            for (int k = 0; k <= hi; k++) used[k] = 0;
        }
        free(used);
    }
    colstart = calloc(ncolors + 1, sizeof(int));
    colidx = malloc(NL * sizeof(int));
    for (int l = 0; l < NL; l++) colstart[lcol[l] + 1]++;
    for (int c = 0; c < ncolors; c++) colstart[c + 1] += colstart[c];
    {
        int *cf = malloc((ncolors + 1) * sizeof(int));
        memcpy(cf, colstart, (ncolors + 1) * sizeof(int));
        for (int l = 0; l < NL; l++) colidx[cf[lcol[l]]++] = l;
        free(cf);
    }
    printf("# edge coloring: %d colors for %d links\n", ncolors, NL);
    free(fill);

    /* --- seed center --- */
    cenx = P.px < 0 ? 0.5 * P.L : P.px;
    ceny = P.py < 0 ? 0.5 * P.L : P.py;
    cenz = P.pz < 0 ? 0.5 * P.L : P.pz;

    /* --- prealign: transverse plane normals (init choice, not dynamics) --- */
    if (P.prealign) {
        /* two distinct transverse normals: the planes' intersection line
         * (the cells' joint axis) is the tilt direction */
        double kk = sqrt(P.kx * P.kx + P.ky * P.ky + P.kz * P.kz);
        double ux = 1, uy = 0, uz = 0;
        if (kk > 0) { ux = P.kx / kk; uy = P.ky / kk; uz = P.kz / kk; }
        double ax = 0, ay = 0, az = 1;
        if (fabs(uz) > 0.9) { ax = 1; az = 0; }
        double dp = ax * ux + ay * uy + az * uz;
        double t1x = ax - dp * ux, t1y = ay - dp * uy, t1z = az - dp * uz;
        double tn = sqrt(t1x * t1x + t1y * t1y + t1z * t1z);
        t1x /= tn; t1y /= tn; t1z /= tn;
        double t2x = uy * t1z - uz * t1y;
        double t2y = uz * t1x - ux * t1z;
        double t2z = ux * t1y - uy * t1x;
        for (int i = 0; i < NC; i++) {
            n1x[i] = t1x; n1y[i] = t1y; n1z[i] = t1z;
            n2x[i] = t2x; n2y[i] = t2y; n2z[i] = t2z;
        }
    }

    /* --- energy seeding (field seeds go into the amplitude pair) --- */
    if (P.noise_amp > 0) {
        /* ambient field churn (any init): the stressor / the room tone */
        for (int i = 0; i < NC; i++) {
            double e0 = frand() * P.noise_amp, ph = frand() * TWO_PI;
            fa1[i] += sqrt(e0) * cos(ph);
            fa2[i] += sqrt(e0) * sin(ph);
        }
    }
    if (P.init == 4) {
        /* standing pairs on existing links: two cells, each the other's
         * bow. Occupancy fixed (pair_x0) or on the tuning curve x*(d)
         * (pair_x0 < 0); phases at the m=1 lock or random. */
        int *used = calloc(NC, sizeof(int));
        for (int l = 0; l < NL && NP < P.npairs && NP < NP_MAX; l++) {
            double d = ld[l];
            if (d < P.pair_dlo || d > P.pair_dhi) continue;
            int i = li[l], j = lj[l];
            if (used[i] || used[j]) continue;
            int clash = 0;
            for (int p = 0; p < NP && !clash; p++) {
                int cand[2] = { i, j }, prev[2] = { ppi[p], ppj[p] };
                for (int a = 0; a < 2 && !clash; a++)
                    for (int b = 0; b < 2 && !clash; b++) {
                        double dx = cx[cand[a]] - cx[prev[b]];
                        double dy = cy[cand[a]] - cy[prev[b]];
                        double dz = cz[cand[a]] - cz[prev[b]];
                        if (dx * dx + dy * dy + dz * dz < P.pair_gap * P.pair_gap)
                            clash = 1;
                    }
            }
            if (clash) continue;
            double x0 = P.pair_x0;
            if (x0 < 0) {
                x0 = (P.w2 * d / (M_PI * P.C) - 1.0) / P.q_detune;  /* x*(d), m=1 */
                if (x0 < 0.02 || x0 > 0.9) continue;                /* no rung reachable */
            }
            double x1 = P.pair_x1 >= 0 ? P.pair_x1 : x0;
            double xa[2] = { x0, x1 };
            int cc[2] = { i, j };
            for (int a = 0; a < 2; a++) {
                int u = cc[a];
                /* target occupancy exactly: the space pull joins Em too */
                double add = xa[a] * P.cap / (1.0 + P.s_pull);
                Em[u] += add;
                double pull = P.s_pull * add;
                double avail = Es[u] - P.es_floor;
                if (pull > avail) pull = avail > 0 ? avail : 0;
                Es[u] -= pull;
                Em[u] += pull;
                /* both planes transverse to the link: axis = the channel */
                double ux = lux[l], uy = luy[l], uz = luz[l];
                double ax = 0, ay = 0, az = 1;
                if (fabs(uz) > 0.9) { ax = 1; az = 0; }
                double dp = ax * ux + ay * uy + az * uz;
                double t1x = ax - dp * ux, t1y = ay - dp * uy, t1z = az - dp * uz;
                double tn = sqrt(t1x * t1x + t1y * t1y + t1z * t1z);
                t1x /= tn; t1y /= tn; t1z /= tn;
                n1x[u] = t1x; n1y[u] = t1y; n1z[u] = t1z;
                n2x[u] = uy * t1z - uz * t1y;
                n2y[u] = uz * t1x - ux * t1z;
                n2z[u] = ux * t1y - uy * t1x;
            }
            double we = P.w2 / (1.0 + P.q_detune * (Em[i] / P.cap));
            th2[i] = frand() * TWO_PI;
            th2[j] = P.pair_seedlock
                     ? fmod((P.pair_q * (th2[i] - we * d / P.C)) / P.pair_p
                            + 8.0 * TWO_PI, TWO_PI)
                     : frand() * TWO_PI;
            lp[l] = (signed char)P.pair_p;
            lq[l] = (signed char)P.pair_q;
            ppi[NP] = i; ppj[NP] = j; ppl[NP] = l;
            pE0[NP] = Em[i] + Em[j];
            used[i] = used[j] = 1;
            NP++;
        }
        free(used);
        printf("# pairs: seeded %d of %d requested, d=[%g,%g], x0=%s%.3g, phases=%s\n",
               NP, P.npairs, P.pair_dlo, P.pair_dhi,
               P.pair_x0 < 0 ? "x*(d), example " : "",
               P.pair_x0 < 0 ? (P.w2 * 1.25 / (M_PI * P.C) - 1.0) / P.q_detune : P.pair_x0,
               P.pair_seedlock ? "locked" : "random");
    } else if (P.init == 7) {
        /* MASS R3 — the ring lock: N voices on a closed chord. Adjacent
         * pairs sit near the m=1 pair rung (omega*d/C = pi at 1:1), so
         * the loop closes when N is even: sum omega*d/C = N*pi. Loads
         * from the tuning curve at each voice's ACTUAL mean spacing;
         * phases seeded in the locked relation around the loop with
         * ring_wind extra turns distributed; axes transverse to the
         * local chord (kappa_align adapts them from there). */
        int n = P.ring_n > 0 ? P.ring_n : 6;
        int *pick = malloc(n * sizeof(int));
        double Rr = P.ring_d / (2.0 * sin(M_PI / n));
        for (int k = 0; k < n; k++) {
            double phk = TWO_PI * k / n;
            double tx = 0.5 * P.L + Rr * cos(phk);
            double ty = 0.5 * P.L + Rr * sin(phk);
            double tz = 0.5 * P.L;
            int best = -1; double bd = 1e30;
            for (int i = 0; i < NC; i++) {
                if (cflag[i]) continue;
                int used = 0;
                for (int q = 0; q < k; q++) if (pick[q] == i) used = 1;
                if (used) continue;
                double dx = cx[i] - tx, dy = cy[i] - ty, dz = cz[i] - tz;
                double dd = dx * dx + dy * dy + dz * dz;
                if (dd < bd) { bd = dd; best = i; }
            }
            pick[k] = best;
        }
        double Lring = 0;
        for (int k = 0; k < n; k++) {
            int u = pick[k], un = pick[(k + 1) % n];
            double ex = cx[un] - cx[u], ey = cy[un] - cy[u], ez = cz[un] - cz[u];
            Lring += sqrt(ex * ex + ey * ey + ez * ez);
        }
        double om_m = P.ring_m >= 0 ? TWO_PI * P.ring_m * P.C / Lring : -1;
        for (int k = 0; k < n; k++) {
            int u = pick[k], up = pick[(k + n - 1) % n], un = pick[(k + 1) % n];
            double dpx = cx[u] - cx[up], dpy = cy[u] - cy[up], dpz = cz[u] - cz[up];
            double dnx = cx[un] - cx[u], dny = cy[un] - cy[u], dnz = cz[un] - cz[u];
            double dp_ = sqrt(dpx * dpx + dpy * dpy + dpz * dpz);
            double dn_ = sqrt(dnx * dnx + dny * dny + dnz * dnz);
            double dbar = 0.5 * (dp_ + dn_);
            /* ring_m: uniform pitch from the ACTUAL loop length — closure
             * exact by construction; the closure integer is the seeded
             * topology class m = N/2 - w (B1 design note) */
            double xk = om_m > 0 ? (P.w2 / om_m - 1.0) / P.q_detune
                      : P.ring_x >= 0 ? P.ring_x
                       : (P.w2 * dbar / (M_PI * P.C) - 1.0) / P.q_detune;
            if (xk < 0.02) xk = 0.02;
            double add = xk * P.cap / (1.0 + P.s_pull);
            Em[u] += add;
            double pull = P.s_pull * add;
            double avail = Es[u] - P.es_floor;
            if (pull > avail) pull = avail > 0 ? avail : 0;
            Es[u] -= pull;
            Em[u] += pull;
            double ux = dnx / dn_, uy = dny / dn_, uz = dnz / dn_;
            double ax = 0, ay = 0, az = 1;
            if (fabs(uz) > 0.9) { ax = 1; az = 0; }
            double dp2 = ax * ux + ay * uy + az * uz;
            double t1x = ax - dp2 * ux, t1y = ay - dp2 * uy, t1z = az - dp2 * uz;
            double tn = sqrt(t1x * t1x + t1y * t1y + t1z * t1z);
            t1x /= tn; t1y /= tn; t1z /= tn;
            n1x[u] = t1x; n1y[u] = t1y; n1z[u] = t1z;
            n2x[u] = uy * t1z - uz * t1y;
            n2y[u] = uz * t1x - ux * t1z;
            n2z[u] = ux * t1y - uy * t1x;
        }
        double closure = 0;
        th2[pick[0]] = frand() * TWO_PI;
        for (int k = 0; k < n; k++) {
            int u = pick[k], un = pick[(k + 1) % n];
            double dnx = cx[un] - cx[u], dny = cy[un] - cy[u], dnz = cz[un] - cz[u];
            double dn_ = sqrt(dnx * dnx + dny * dny + dnz * dnz);
            double we = P.w2 / (1.0 + P.q_detune * (Em[u] / P.cap));
            closure += we * dn_ / P.C;
            if (k < n - 1)
                th2[un] = fmod(th2[u] - we * dn_ / P.C
                               + TWO_PI * (double)P.ring_wind / n
                               + 8.0 * TWO_PI, TWO_PI);
        }
        printf("# ring: seeded n=%d R=%.2f d_target=%.2f Lring=%.3f "
               "closure/2pi=%.4f wind=%d ring_m=%d\n",
               n, Rr, P.ring_d, Lring, closure / TWO_PI, P.ring_wind, P.ring_m);
        free(pick);
    } else if (P.init == 8) {
        /* H-series — the pressurized bubble (user proposition 2026-07-28:
         * a stable particle has strong surface tension and a density
         * hole, the "higgs hole" — dense-empty, space-rich, load-
         * bearing). Consonant closed surfaces must be BIPARTITE: the
         * reachable pair rung is the antiphase pi-rung, which cannot
         * close around odd cycles — triangulated shells are frustrated;
         * the CUBE is the minimal consonant shell (all 4-cycles; face
         * diagonals beyond the contact ceiling). Uniform pitch from the
         * ACTUAL mean edge (pi-rung), phases by spanning tree in the
         * bipartite classes; the core cells get a space overpressure
         * (an initial condition; E0 is computed after seeding). */
        double ccx = 0.5 * P.L, ccy = 0.5 * P.L, ccz = 0.5 * P.L;
        int nsh = P.shell_n == 8 ? 8 : 0;
        int pick[8];
        if (nsh == 8) {
            double h = 0.5 * P.shell_a;
            for (int k = 0; k < 8; k++) {
                double tx = ccx + ((k & 1) ? h : -h);
                double ty = ccy + ((k & 2) ? h : -h);
                double tz = ccz + ((k & 4) ? h : -h);
                int best = -1; double bd = 1e30;
                for (int i = 0; i < NC; i++) {
                    if (cflag[i]) continue;
                    /* the hole is load-bearing: vertex picks may not
                     * consume the core ball */
                    double hx = cx[i] - ccx, hy = cy[i] - ccy, hz = cz[i] - ccz;
                    if (P.r_core > 0
                        && hx * hx + hy * hy + hz * hz < P.r_core * P.r_core)
                        continue;
                    int used = 0;
                    for (int q = 0; q < k; q++) if (pick[q] == i) used = 1;
                    if (used) continue;
                    double dx = cx[i] - tx, dy = cy[i] - ty, dz = cz[i] - tz;
                    double dd = dx * dx + dy * dy + dz * dz;
                    if (dd < bd) { bd = dd; best = i; }
                }
                pick[k] = best;
            }
            /* refine picks for EDGE UNIFORMITY: for each vertex try its
             * runner-up candidates, keep swaps that reduce the spread
             * of the 12 edge lengths (foam jitter otherwise concentrates
             * cycle defects on co-tree edges — measured gates min 0) */
            for (int pass = 0; pass < 3; pass++) {
                for (int k = 0; k < 8; k++) {
                    double h = 0.5 * P.shell_a;
                    double tx = ccx + ((k & 1) ? h : -h);
                    double ty = ccy + ((k & 2) ? h : -h);
                    double tz = ccz + ((k & 4) ? h : -h);
                    double bestscore = 1e30; int bestc = pick[k];
                    for (int i = 0; i < NC; i++) {
                        if (cflag[i]) continue;
                        double hx = cx[i] - ccx, hy = cy[i] - ccy, hz = cz[i] - ccz;
                        if (P.r_core > 0
                            && hx * hx + hy * hy + hz * hz < P.r_core * P.r_core)
                            continue;
                        double dx = cx[i] - tx, dy = cy[i] - ty, dz = cz[i] - tz;
                        if (dx * dx + dy * dy + dz * dz > 1.44) continue;
                        int used = 0;
                        for (int q = 0; q < 8; q++)
                            if (q != k && pick[q] == i) used = 1;
                        if (used) continue;
                        double sc = 0;
                        for (int b = 0; b < 3; b++) {
                            int k2 = k ^ (1 << b);
                            double ex = cx[pick[k2]] - cx[i];
                            double ey = cy[pick[k2]] - cy[i];
                            double ez = cz[pick[k2]] - cz[i];
                            double dd = sqrt(ex * ex + ey * ey + ez * ez);
                            sc += (dd - P.shell_a) * (dd - P.shell_a);
                        }
                        if (sc < bestscore) { bestscore = sc; bestc = i; }
                    }
                    pick[k] = bestc;
                }
            }
            /* cube edges: vertices differing in exactly one axis bit */
            double abar = 0; int ne = 0;
            for (int k = 0; k < 8; k++)
                for (int b = 0; b < 3; b++) {
                    int k2 = k ^ (1 << b);
                    if (k2 < k) continue;
                    double ex = cx[pick[k2]] - cx[pick[k]];
                    double ey = cy[pick[k2]] - cy[pick[k]];
                    double ez = cz[pick[k2]] - cz[pick[k]];
                    abar += sqrt(ex * ex + ey * ey + ez * ez);
                    ne++;
                }
            abar /= ne;
            double om = M_PI * P.C / abar;
            double xs = P.shell_x >= 0 ? P.shell_x
                       : (P.w2 / om - 1.0) / P.q_detune;
            if (xs < 0.02) xs = 0.02;
            for (int k = 0; k < 8; k++) {
                int u = pick[k];
                double add = xs * P.cap / (1.0 + P.s_pull);
                Em[u] += add;
                double pull = P.s_pull * add;
                double avail = Es[u] - P.es_floor;
                if (pull > avail) pull = avail > 0 ? avail : 0;
                Es[u] -= pull;
                Em[u] += pull;
                /* planes transverse to the radial direction (the surface) */
                double ux = cx[u] - ccx, uy = cy[u] - ccy, uz = cz[u] - ccz;
                double un = sqrt(ux * ux + uy * uy + uz * uz);
                if (un > 1e-9) { ux /= un; uy /= un; uz /= un; }
                double ax = 0, ay = 0, az = 1;
                if (fabs(uz) > 0.9) { ax = 1; az = 0; }
                double dp2 = ax * ux + ay * uy + az * uz;
                double t1x = ax - dp2 * ux, t1y = ay - dp2 * uy, t1z = az - dp2 * uz;
                double tn = sqrt(t1x * t1x + t1y * t1y + t1z * t1z);
                t1x /= tn; t1y /= tn; t1z /= tn;
                n1x[u] = t1x; n1y[u] = t1y; n1z[u] = t1z;
                n2x[u] = uy * t1z - uz * t1y;
                n2y[u] = uz * t1x - ux * t1z;
                n2z[u] = ux * t1y - uy * t1x;
            }
            /* phases: BFS over cube edges from vertex 0, seedlock each */
            int seen[8] = {0}, qb[8], hq = 0, tq = 0;
            th2[pick[0]] = frand() * TWO_PI;
            seen[0] = 1; qb[tq++] = 0;
            while (hq < tq) {
                int k = qb[hq++];
                for (int b = 0; b < 3; b++) {
                    int k2 = k ^ (1 << b);
                    if (seen[k2]) continue;
                    double ex = cx[pick[k2]] - cx[pick[k]];
                    double ey = cy[pick[k2]] - cy[pick[k]];
                    double ez = cz[pick[k2]] - cz[pick[k]];
                    double dd = sqrt(ex * ex + ey * ey + ez * ez);
                    th2[pick[k2]] = fmod(th2[pick[k]] - om * dd / P.C
                                         + 8.0 * TWO_PI, TWO_PI);
                    seen[k2] = 1; qb[tq++] = k2;
                }
            }
            /* report seed gate quality over all 12 edges */
            double gmin = 1, gsum = 0;
            for (int k = 0; k < 8; k++)
                for (int b = 0; b < 3; b++) {
                    int k2 = k ^ (1 << b);
                    if (k2 < k) continue;
                    double ex = cx[pick[k2]] - cx[pick[k]];
                    double ey = cy[pick[k2]] - cy[pick[k]];
                    double ez = cz[pick[k2]] - cz[pick[k]];
                    double dd = sqrt(ex * ex + ey * ey + ez * ez);
                    double ps = wrap_pi(th2[pick[k]] - om * dd / P.C
                                        - th2[pick[k2]]);
                    double g = gate_of(ps);
                    if (g < gmin) gmin = g;
                    gsum += g;
                }
            printf("# shell: cube a_target=%.2f abar=%.3f omega=%.4f "
                   "x=%.4f gates min=%.3f mean=%.3f\n",
                   P.shell_a, abar, om, xs, gmin, gsum / 12.0);
        }
        /* the hole: core space overpressure (skip shell cells) */
        if (P.es_core > 0) {
            int ncore = 0;
            for (int i = 0; i < NC; i++) {
                if (cflag[i]) continue;
                int isv = 0;
                for (int q = 0; q < nsh; q++) if (pick[q] == i) isv = 1;
                if (isv) continue;
                double dx = cx[i] - ccx, dy = cy[i] - ccy, dz = cz[i] - ccz;
                if (dx * dx + dy * dy + dz * dz < P.r_core * P.r_core) {
                    Es[i] = P.es_core;
                    ncore++;
                }
            }
            while (ncore < 1) {
                /* never seed an empty hole: take the nearest free
                 * non-vertex cell */
                int best = -1; double bd = 1e30;
                for (int i = 0; i < NC; i++) {
                    if (cflag[i] || Es[i] == P.es_core) continue;
                    int isv = 0;
                    for (int q = 0; q < nsh; q++) if (pick[q] == i) isv = 1;
                    if (isv) continue;
                    double dx = cx[i] - ccx, dy = cy[i] - ccy, dz = cz[i] - ccz;
                    double dd = dx * dx + dy * dy + dz * dz;
                    if (dd < bd) { bd = dd; best = i; }
                }
                if (best < 0) break;
                Es[best] = P.es_core;
                ncore++;
            }
            printf("# shell core: es_core=%.2f r_core=%.2f cells=%d\n",
                   P.es_core, P.r_core, ncore);
        }
    } else if (P.init == 5) {
        /* double slit: wall (detuned medium) with two vacuum windows,
         * screen + edge sinks (record media), optional which-path
         * recorders inside the slits, plane-wave source slab. */
        double yA = 0.5 * P.L - 0.5 * P.slit_sep;
        double yB = 0.5 * P.L + 0.5 * P.slit_sep;
        int nwall = 0, nscr = 0, nsink = 0, ndet = 0, nsrc = 0;
        for (int i = 0; i < NC; i++) {
            double m = P.sink_m;
            if (cy[i] < m || cy[i] > P.L - m || cz[i] < m || cz[i] > P.L - m
                || cx[i] < 0.8) {
                cflag[i] = 3; nsink++;
            } else if (cx[i] > P.screen_x) {
                cflag[i] = 2; nscr++;
            } else if (fabs(cx[i] - P.wall_x) < 0.5 * P.wall_th) {
                int inA = fabs(cy[i] - yA) < P.slit_hw;
                int inB = fabs(cy[i] - yB) < P.slit_hw;
                if (P.block_slit == 1) inA = 0;
                if (P.block_slit == 2) inB = 0;
                if (inA || inB) {
                    int det = (P.which_path == 1) || (P.which_path == 2 && inA);
                    if (det) { cflag[i] = 4; ndet++; }
                } else {
                    cflag[i] = 1; nwall++;
                }
            }
            if (cflag[i] == 0 && cx[i] > P.src_lo && cx[i] < P.src_hi) {
                double ph = -P.kx * (cx[i] - P.src_lo);
                fa1[i] += sqrt(P.amp) * cos(ph);
                fa2[i] += sqrt(P.amp) * sin(ph);
                nsrc++;
            }
        }
        printf("# slit: wall=%d screen=%d sink=%d recorder=%d source_cells=%d slitA_y=%.2f slitB_y=%.2f\n",
               nwall, nscr, nsink, ndet, nsrc, yA, yB);
    } else if (P.init == 6) {
        /* tier 4 — HOM: two parallel waveguides along x (mute elsewhere),
         * a tunneling coupler where the gap detune drops to hom_barrier,
         * end sinks. Quantum A launched in wire 1 (component H), quantum
         * B in wire 2 (component V), offset by hom_dx (the delay). */
        double yw1 = 0.5 * P.L - 3.0, yw2 = 0.5 * P.L + 3.0, hw = 1.5;
        double zc = 0.5 * P.L;
        double x0 = P.px < 0 ? 7.0 : P.px;
        int nw = 0, ncp = 0, nsk = 0;
        for (int i = 0; i < NC; i++) {
            int in1 = fabs(cy[i] - yw1) < hw && fabs(cz[i] - zc) < hw;
            int in2 = fabs(cy[i] - yw2) < hw && fabs(cz[i] - zc) < hw;
            int ingap = cy[i] > yw1 + hw && cy[i] < yw2 - hw
                        && fabs(cz[i] - zc) < hw
                        && cx[i] > P.hom_c0 && cx[i] < P.hom_c1;
            if (in1 || in2) {
                if (cx[i] < 1.5 || cx[i] > P.L - 1.5) { cflag[i] = 3; nsk++; }
                else { cflag[i] = 0; nw++; }
            } else if (ingap) {
                cflag[i] = 5; ncp++;
            } else {
                cflag[i] = 1;
            }
        }
        double s2 = 2.0 * P.sigma * P.sigma;
        for (int i = 0; i < NC; i++) {
            if (cflag[i] != 0) continue;
            int in1 = fabs(cy[i] - yw1) < hw;
            double xs = in1 ? x0 : x0 + P.hom_dx;
            double g = exp(-(cx[i] - xs) * (cx[i] - xs) / s2);
            if (g < 1e-4) continue;
            double ph = -P.kx * (cx[i] - xs);
            if (in1) {
                fa1[i] += sqrt(P.amp * g) * cos(ph);
                fa2[i] += sqrt(P.amp * g) * sin(ph);
            } else if (!P.hom_solo) {
                fb1[i] += sqrt(P.amp * g) * cos(ph);
                fb2[i] += sqrt(P.amp * g) * sin(ph);
            }
        }
        printf("# hom: wire_cells=%d coupler_cells=%d sinks=%d wires_y=%.1f/%.1f coupler_x=[%g,%g]\n",
               nw, ncp, nsk, yw1, yw2, P.hom_c0, P.hom_c1);
    } else if (P.init == 2 || P.init == 3) {
        double s2 = 2.0 * P.sigma * P.sigma;
        for (int i = 0; i < NC; i++) {
            double dx = cx[i] - cenx, dy = cy[i] - ceny, dz = cz[i] - cenz;
            double rr2 = dx * dx + dy * dy + dz * dz;
            double g = exp(-rr2 / s2);
            if (g < 1e-3) continue;
            double tilt = -(P.kx * dx + P.ky * dy + P.kz * dz);
            if (P.init == 2) {
                fa1[i] += sqrt(P.amp * g) * cos(tilt);
                fa2[i] += sqrt(P.amp * g) * sin(tilt);
            } else {
                double add = P.amp * g;
                if (Em[i] + add > 0.95 * P.cap) add = 0.95 * P.cap - Em[i];
                if (add < 0) add = 0;
                Em[i] += add;
                /* the blob was made by converting local space: consistent start */
                double pull = P.s_pull * add;
                double avail = Es[i] - P.es_floor;
                if (pull > avail) pull = avail > 0 ? avail : 0;
                Es[i] -= pull;
                Em[i] += pull;
                th2[i] = tilt;
            }
        }
    }

    /* absorbing edges (apparatus): flag a margin of cells on every face
     * as sinks (record media), for open-space experiments under any
     * init — e.g. so a radiation-pressure strike is not followed by the
     * unabsorbed remainder reverberating in the closed box. */
    if (P.edge_sink > 0) {
        int nes = 0;
        double m = P.edge_sink;
        for (int i = 0; i < NC; i++) {
            if (cflag[i]) continue;
            if (cx[i] < m || cx[i] > P.L - m || cy[i] < m || cy[i] > P.L - m
                || cz[i] < m || cz[i] > P.L - m) { cflag[i] = 3; nes++; }
        }
        printf("# edge_sink: margin=%g cells=%d\n", m, nes);
    }

    /* auxiliary field pulse (apparatus): a second packet on top of any
     * init — e.g. aimed at a seeded blob for the radiation-pressure
     * experiment. Same construction as init=pulse. */
    if (P.aux_amp > 0) {
        double axc = P.aux_x < 0 ? 0.5 * P.L : P.aux_x;
        double ayc = P.aux_y < 0 ? 0.5 * P.L : P.aux_y;
        double azc = P.aux_z < 0 ? 0.5 * P.L : P.aux_z;
        double asg = P.aux_sigma > 0 ? P.aux_sigma : P.sigma;
        double s2 = 2.0 * asg * asg;
        for (int i = 0; i < NC; i++) {
            double dx = cx[i] - axc, dy = cy[i] - ayc, dz = cz[i] - azc;
            double rr2 = dx * dx + dy * dy + dz * dz;
            double g = exp(-rr2 / s2);
            if (g < 1e-3) continue;
            double tilt = -(P.aux_kx * dx + P.aux_ky * dy + P.aux_kz * dz);
            fa1[i] += sqrt(P.aux_amp * g) * cos(tilt);
            fa2[i] += sqrt(P.aux_amp * g) * sin(tilt);
        }
    }

    /* tier 1 — claim rule: index the screen for the absorption process */
    if (P.init == 5 && P.n_quanta > 0) {
        nscr = 0;
        for (int i = 0; i < NC; i++) if (cflag[i] == 2) nscr++;
        scr_id = malloc(nscr * sizeof(int));
        scr_prev = calloc(nscr, sizeof(double));
        qinc = calloc((size_t)QINT_MAX * (size_t)nscr, sizeof(float));
        qtms = calloc(QINT_MAX, sizeof(double));
        int q = 0;
        for (int i = 0; i < NC; i++) if (cflag[i] == 2) scr_id[q++] = i;
    }

    /* field caches from the amplitude pair */
    for (int i = 0; i < NC; i++) {
        Ee[i] = fa1[i] * fa1[i] + fa2[i] * fa2[i]
              + fb1[i] * fb1[i] + fb2[i] * fb2[i];
        if (Ee[i] > 1e-20) th1[i] = atan2(fa2[i], fa1[i]);
    }

    /* initial geometry -> A0 baseline; effective frequencies for t=0 diag */
    for (int i = 0; i < NC; i++) {
        double ratio = Es[i] / P.e_s0;
        cr[i] = cr0[i] * cbrt(ratio > 0 ? ratio : 0);
        double det = 1.0 + P.q_detune * Em[i] / P.cap;
        w1e[i] = P.w1 / det;
        w2e[i] = P.w2 / det;
    }
    for (int l = 0; l < NL; l++) {
        double d = ld[l], ri = cr[li[l]], rj = cr[lj[l]];
        double A = 0;
        if (d < ri + rj) {
            double t = d * d - rj * rj + ri * ri;
            double a2 = (4.0 * d * d * ri * ri - t * t) / (4.0 * d * d);
            if (a2 > 0) A = M_PI * a2;
            else if (d < fabs(ri - rj)) { double rm = ri < rj ? ri : rj; A = M_PI * rm * rm; }
        }
        lA[l] = lA0[l] = A;
    }

    double mdeg = NC ? 2.0 * (double)NL / NC : 0;
    printf("# fabric: cells=%d channels=%d mean_degree=%.2f box=%g dmin=%g\n",
           NC, NL, mdeg, P.L, P.dmin);
    if (P.snap_every > 0) mkdir(P.snap_dir, 0755);
}

/* ------------------------------------------------------------------ */
/* one step                                                            */
/* ------------------------------------------------------------------ */

/* Exact paired ledger moves between the dense/space ledgers and the
 * field amplitude: injection scales |psi| up by the exact energy (seeding
 * the phase from the dense clock when the field is silent [P]). */
static void field_inject(int i, double dE)
{
    if (dE <= 0) return;
    double e = fa1[i] * fa1[i] + fa2[i] * fa2[i];
    if (e > 1e-20) {
        double fac = sqrt((e + dE) / e);
        fa1[i] *= fac;
        fa2[i] *= fac;
    } else {
        fa1[i] = sqrt(dE) * cos(th2[i]);
        fa2[i] = sqrt(dE) * sin(th2[i]);
    }
    Ee[i] = e + dE;
}

/* S1/U2 — atoms at mode boundaries (ROADMAP §6.4), the uniform law:
 * transport WITHIN a mode is continuous; energy CONVERTING between modes
 * moves only in whole action atoms eps = A0*w/2pi, evaluated at the
 * converting channel's own completed cycles (beat cycles for the beat
 * channel, flight-cycle deliveries for roughness). quant_mode selects the
 * variant law — one value for a whole battery, never per experiment.
 * Two independent choices, four uniform variants:
 *   sizing: SRC (the sending mode's quantum) or DST (the receiving
 *           mode's quantum — a structural gap when eps_dst > eps_src);
 *   memory: FLOOR (sub-atom demand lapses with its cycle) or CREDIT
 *           (demand integrates across cycles, lapses at 2 atoms).
 *   1 = SRC+FLOOR   2 = SRC+CREDIT   3 = DST+FLOOR   4 = DST+CREDIT
 * A0eff = 0 is the continuous limit (all variants coincide, exactly). */
static double atoms_eps(double eps_src, double eps_dst)
{
    return P.quant_mode >= 3 ? eps_dst : eps_src;
}

/* the frequency that sizes the atom (for the QATOM linearity ledger) */
static double atoms_w(double w_src, double w_dst)
{
    return P.quant_mode >= 3 ? w_dst : w_src;
}

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

/* clamp a fired amount to an affordability ceiling while staying on the
 * atom grid; unfired atoms return to credit under the credit variant */
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
    if ((qfire_n++ % 50) == 0)
        printf("# QATOM t=%.2f dir=%s w=%.9g e=%.12g\n",
               sim_t, fd ? "FD" : "DF", w, e);
}

static void step_field(void)
{
    double dt = P.dt;
    double Aref = M_PI * P.r0 * P.r0;
    double dref = 2.0 * P.r0;
    double gm = P.gamma_res_m < 0 ? P.gamma_res : P.gamma_res_m;
    double G2c[2] = { P.gamma_res * P.gamma_res, gm * gm };
    double lockf = P.lock_floor * P.cap;

    /* pass 0: pitch share of in-flight dense energy. A cell's pitch load
     * counts what it is bound INTO: its store plus half of each incident
     * channel's in-flight dense energy — the channel is a joint process
     * of its two ends, so energy in transit between a pair's voices does
     * not lighten the pair (found via e7: the mutual-flight shelter was
     * detuning the very pair it shelters — every pair sharp of the rung). */
#pragma omp parallel for schedule(static)
    for (int i = 0; i < NC; i++) {
        double fl = 0;
        for (int q = cls[i]; q < cls[i + 1]; q++) {
            int l = clidx[q];
            fl += lem[SLOT(l, 1, 0)] + lem[SLOT(l, 1, 1)];
        }
        flload[i] = 0.5 * fl;
    }

    /* pass S: space-mode transport — PRESSURE PUSHES, nothing reaches out.
     * Space is a mode like the others, so it flows through the same
     * channels: each link carries space from its higher-pressure end
     * toward its lower-pressure end, at the channel's own conductance.
     * The flux ORIGINATES from the pressurized side's store (outflow-
     * limited above the conversion floor); an empty cell never draws —
     * there is no suction. Pressure carries a displacement term from
     * load: matter IS converted space, so a loaded cell pushes space out
     * (pi = Es + s_disp*(Em+Ee)), which is how a mass MAINTAINS its
     * depression instead of having it merely frozen in. Exactly paired
     * ledger moves; the depression a mass holds is the gravitational
     * footprint the G-battery measures. */
    if (P.s_k > 0) {
        /* three deterministic phases (Jacobi form): per-link wants from
         * the pre-step state; per-cell outflow scales; per-cell gather
         * apply — every link's resolved move computed identically from
         * both endpoints, so the paired books close with no scatter. */
#pragma omp parallel for schedule(static)
        for (int l = 0; l < NL; l++) {
            swl[l] = 0.0;
            if (lA[l] <= 0) continue;
            int i = li[l], j = lj[l];
            double pi_ = Es[i] + P.s_disp * (Em[i] + Ee[i]);
            double pj_ = Es[j] + P.s_disp * (Em[j] + Ee[j]);
            double dp = pi_ - pj_;
            if (dp == 0) continue;
            double w = (lA[l] / Aref) * (dref / ld[l]);
            swl[l] = P.s_k * dt * w * dp;
        }
#pragma omp parallel for schedule(static)
        for (int i = 0; i < NC; i++) {
            double rq = 0;
            for (int q = cls[i]; q < cls[i + 1]; q++) {
                int l = clidx[q];
                double f = swl[l];
                if ((f > 0 && li[l] == i) || (f < 0 && lj[l] == i))
                    rq += fabs(f);
            }
            sprq[i] = rq;
            double avail = Es[i] - P.es_floor;
            sscl[i] = (avail <= 0) ? 0.0
                     : (rq > avail ? avail / rq : 1.0);
        }
#pragma omp parallel for schedule(static)
        for (int i = 0; i < NC; i++) {
            double de = 0;
            for (int q = cls[i]; q < cls[i + 1]; q++) {
                int l = clidx[q];
                double f = swl[l];
                if (f == 0) continue;
                int src = f > 0 ? li[l] : lj[l];
                double mag = fabs(f) * sscl[src];
                de += (src == i) ? -mag : mag;
            }
            Es[i] += de;
        }
        /* G4 instrument (apparatus-gated, serial): radial projection of
         * the resolved moves, binned by shell (positive = outward) */
        if (P.rad_diag) {
            for (int l = 0; l < NL; l++) {
                double f = swl[l];
                if (f == 0) continue;
                int src = f > 0 ? li[l] : lj[l];
                double mag = fabs(f) * sscl[src];
                if (mag == 0) continue;
                double mx = 0.5 * (cx[li[l]] + cx[lj[l]]) - cenx;
                double my = 0.5 * (cy[li[l]] + cy[lj[l]]) - ceny;
                double mz = 0.5 * (cz[li[l]] + cz[lj[l]]) - cenz;
                double rr = sqrt(mx * mx + my * my + mz * mz);
                int k = (int)(rr / (0.5 * P.L / 8.0));
                if (rr > 1e-9 && k < 8) {
                    double sgn = src == li[l] ? 1.0 : -1.0;
                    srad[k] += mag * sgn * (lux[l] * mx + luy[l] * my
                                            + luz[l] * mz) / rr;
                }
            }
        }
    }

    /* pass 1: live radii, effective frequencies, clear requests */
#pragma omp parallel for schedule(static)
    for (int i = 0; i < NC; i++) {
        double ratio = Es[i] / P.e_s0;
        cr[i] = cr0[i] * cbrt(ratio > 0 ? ratio : 0);
        /* S1/U3 — pitch load is BOUND energy: the dense store (plus the
         * cell's share of in-flight dense processes) re-pitches the
         * cell's harmonics; passing field amplitude does not. One
         * definition everywhere: vacuum optics is linear automatically,
         * Kerr nonlinearity is a property of loaded matter. */
        double x = (Em[i] + flload[i]) / P.cap;
        if (cflag[i] >= 2) x = 0;   /* record media: deep and pitch-stable */
        double det = 1.0 + P.q_detune * x;
        double wb = P.w1;
        if (cflag[i] == 1) wb += P.wall_detune;   /* the wall: a mute — different pitch */
        else if (cflag[i] == 5) wb += P.hom_barrier;   /* coupler: a translucent mute */
        w1e[i] = wb / det;
        w2e[i] = P.w2 / det;
        req0[i] = req1[i] = 0.0;
    }

    /* pass 2: channel geometry + desired deposits (resonant joining) */
#pragma omp parallel for schedule(static)
    for (int l = 0; l < NL; l++) {
        int i = li[l], j = lj[l];
        double d = ld[l], ri = cr[i], rj = cr[j];
        double A = 0;
        if (d < ri + rj) {
            double t = d * d - rj * rj + ri * ri;
            double a2 = (4.0 * d * d * ri * ri - t * t) / (4.0 * d * d);
            if (a2 > 0) A = M_PI * a2;
            else if (d < fabs(ri - rj)) { double rm = ri < rj ? ri : rj; A = M_PI * rm * rm; }
        }
        lA[l] = A;
        lwant[SLOT(l, 0, 0)] = lwant[SLOT(l, 0, 1)] = 0.0;
        lwant[SLOT(l, 1, 0)] = lwant[SLOT(l, 1, 1)] = 0.0;
        lflux[2 * l] = lflux[2 * l + 1] = 0.0;
        if (A <= 0) continue;

        double geo = (A / Aref) * (dref / d);
        double occ_i = (Em[i] + Ee[i]) / P.cap;
        double occ_j = (Em[j] + Ee[j]) / P.cap;
        double head_i = 1.0 - occ_i; if (head_i < 0) head_i = 0; if (head_i > 1) head_i = 1;
        double head_j = 1.0 - occ_j; if (head_j < 0) head_j = 0; if (head_j > 1) head_j = 1;
        /* record media are deep: their grains do not throttle acceptance */
        if (cflag[i] >= 2 && cflag[i] <= 3) head_i = 1.0;
        if (cflag[j] >= 2 && cflag[j] <= 3) head_j = 1.0;

        /* the two planes of a cell intersect in a line — the cell's axis.
         * A channel carries only along directions lying in BOTH planes:
         * this collimates transfer to the axis the planes jointly define. */
        double d1i = n1x[i] * lux[l] + n1y[i] * luy[l] + n1z[i] * luz[l];
        double d2i = n2x[i] * lux[l] + n2y[i] * luy[l] + n2z[i] * luz[l];
        double d1j = n1x[j] * lux[l] + n1y[j] * luy[l] + n1z[j] * luz[l];
        double d2j = n2x[j] * lux[l] + n2y[j] * luy[l] + n2z[j] * luz[l];
        double axi = (1.0 - d1i * d1i) * (1.0 - d2i * d2i); if (axi < 0) axi = 0;
        double axj = (1.0 - d1j * d1j) * (1.0 - d2j * d2j); if (axj < 0) axj = 0;
        double axij = axi * axj;
        if (axij < 1e-10) continue;

        for (int c = 0; c < 2; c++) {
            /* the field sector no longer uses gated transport (pass F:
             * unitary amplitude dynamics); dense transport unchanged */
            if (c == 0) continue;
            /* locked grains do not migrate: no dense transfer at instruments */
            if (c == 1 && (cflag[i] || cflag[j])) continue;
            double mob_i = c == 0 ? Ee[i] : Em[i];
            double mob_j = c == 0 ? Ee[j] : Em[j];
            if (mob_i < 1e-15 && mob_j < 1e-15) continue;

            double nix = c == 0 ? n1x[i] : n2x[i], niy = c == 0 ? n1y[i] : n2y[i], niz = c == 0 ? n1z[i] : n2z[i];
            double njx = c == 0 ? n1x[j] : n2x[j], njy = c == 0 ? n1y[j] : n2y[j], njz = c == 0 ? n1z[j] : n2z[j];
            double dnn = nix * njx + niy * njy + niz * njz;
            double gpl = axij * dnn * dnn;
            if (gpl < 1e-8) continue;

            double wi = c == 0 ? w1e[i] : w2e[i];
            double wj = c == 0 ? w1e[j] : w2e[j];
            double thi = c == 0 ? th1[i] : th2[i];
            double thj = c == 0 ? th1[j] : th2[j];

            double res, g_ij, g_ji;
            double mi_eff = mob_i, mj_eff = mob_j;
            if (c == 1) {
                /* C1: sympathetic joining through coincident partials
                 * q*wi ~ p*wj; tongues narrow and weaken with complexity
                 * (Tenney height) — the interval vocabulary as dynamics */
                double best = -1.0;
                int bp = 1, bq = 1;
                for (int tt = 0; tt < ncomb; tt++) {
                    double pp = combp[tt], qq = combq[tt];
                    double gw = gm / (pp * qq);
                    double det = qq * wi - pp * wj;
                    double sc = (gw * gw / (gw * gw + det * det)) / (pp * qq);
                    if (sc > best) { best = sc; bp = combp[tt]; bq = combq[tt]; }
                }
                res = best;
                lp[l] = (signed char)bp;
                lq[l] = (signed char)bq;
                g_ij = gate_of(bq * thi - bq * wi * d / P.C - bp * thj);
                g_ji = gate_of(bp * thj - bp * wj * d / P.C - bq * thi);
                if (P.mob_sym) {
                    /* C3: harmony is mutual — exchange takes two; a silent
                     * partner responds at its sympathetic readiness floor */
                    double fl0 = P.mob_floor * P.cap;
                    double mjr = mob_j > fl0 ? mob_j : fl0;
                    double mir = mob_i > fl0 ? mob_i : fl0;
                    mi_eff = sqrt(mob_i * mjr);
                    mj_eff = sqrt(mob_j * mir);
                }
            } else {
                double dw = wi - wj;
                res = G2c[0] / (G2c[0] + dw * dw);
                g_ij = gate_of(thi - wi * d / P.C - thj);
                g_ji = gate_of(thj - wj * d / P.C - thi);
            }

            double kd = c == 0 ? P.k_dep : P.k_dep * P.k_dep_m;
            double base = kd * dt * geo * gpl * res;
            double w_ij = base * g_ij * head_j * mi_eff;
            double w_ji = base * g_ji * head_i * mj_eff;
            if (c == 1 && P.kappa_freq > 0) {
                /* A1 — the choir's correction: sympathetic exchange
                 * carries the DISPERSIVE partner of the comb resonance
                 * (the odd lineshape 2*det*G/(det^2+G^2)): net flow is
                 * biased toward the direction that pulls the coincident
                 * partials into coincidence. det>0 means voice i is
                 * sharp (too light) — feeding i flattens it. Vanishes
                 * at exact tune and far off resonance, so silent-room
                 * bleed is not amplified. Wants only; every deposit
                 * remains an exactly paired ledger move. */
                double detw = lq[l] * wi - lp[l] * wj;
                double gwb = gm / (lp[l] * lq[l]);
                double g2b = gwb * gwb;
                /* windowed inside the acceptance (x Lorentzian): the
                 * reactive pull belongs to the resonance it serves —
                 * far-off-resonance links (blob rim vs vacuum) carry no
                 * bias, so recruitment pressure cannot creep a sealed
                 * blob. Falls as 1/det^3 outside, ~linear inside. */
                double Db = 2.0 * detw * gwb * g2b
                            / ((detw * detw + g2b) * (detw * detw + g2b));
                /* ...and only for those singing together: the correction
                 * rides the mutual gate closure, so a locked pair is
                 * pulled into tune while an unlocked graded profile (a
                 * blob rim) feels no homogenizing pressure. */
                double kb = P.kappa_freq * Db * g_ij * g_ji;
                double bi = 1.0 - kb; if (bi < 0) bi = 0;
                double bj = 1.0 + kb; if (bj < 0) bj = 0;
                w_ij *= bi;
                w_ji *= bj;
            }
            if (c == 1 && P.kappa_reac > 0) {
                /* S2 — the choir's correction, DERIVED: sympathetic
                 * exchange carries its reactive (odd) component — the
                 * interference cross-flow of the coincident partials,
                 * which even-gate rate compression discards (with even
                 * gates on a rung, g_ij == g_ji for a locked pair, so
                 * no odd-in-det flow exists at rate level; the odd
                 * channel is interference energy). kappa is the
                 * channel's own conductance read as amplitude coupling:
                 * kappa_reac = 1 is the unitarity point, not a tuned
                 * strength (v89/s2/choir_pull.c: sign restoring, window
                 * ~2*Gamma_b inherent in res x lock loss, rim
                 * protection inherent in sin(psi) averaging to zero for
                 * unlocked pairs). */
                double ps_ij = wrap_pi(lq[l] * thi - lq[l] * wi * d / P.C
                                       - lp[l] * thj);
                double ps_ji = wrap_pi(lp[l] * thj - lp[l] * wj * d / P.C
                                       - lq[l] * thi);
                double Sm = sqrt(mi_eff * mj_eff);
                double hh = sqrt(head_i * head_j);
                double reac = P.kappa_reac * 0.5 * base * hh * Sm;
                w_ij -= reac * g_ij * sin(ps_ij);
                w_ji -= reac * g_ji * sin(ps_ji);
                if (w_ij < 0) w_ij = 0;
                if (w_ji < 0) w_ji = 0;
            }
            if (w_ij > 0) lwant[SLOT(l, c, 0)] = w_ij;
            if (w_ji > 0) lwant[SLOT(l, c, 1)] = w_ji;
        }
    }

    /* pass 3: outflow limiter (never overdraw a ledger) — per-cell
     * gather of this cell's outbound wants over its incident links */
#pragma omp parallel for schedule(static)
    for (int i = 0; i < NC; i++) {
        double r0 = 0, r1 = 0;
        for (int q = cls[i]; q < cls[i + 1]; q++) {
            int l = clidx[q];
            int dir = li[l] == i ? 0 : 1;
            r0 += lwant[SLOT(l, 0, dir)];
            r1 += lwant[SLOT(l, 1, dir)];
        }
        req0[i] = r0; req1[i] = r1;
        double a0 = 0.98 * Ee[i], a1 = 0.98 * Em[i];
        scl0[i] = (r0 > a0 && r0 > 0) ? a0 / r0 : 1.0;
        scl1[i] = (r1 > a1 && r1 > 0) ? a1 / r1 : 1.0;
    }

    /* pass 4: apply deposits into channel flight (paired ledger moves).
     * Resonant joining begins at the boundary: each deposit entrains the
     * receiver's clock toward the retarded tail, weighted by how much is
     * arriving against what the receiver already coherently holds.
     * Three deterministic phases: resolve deposits per link (each slot
     * owned by one link), debit sources per cell, entrain receivers per
     * cell — same books, no scatter. Foreign clocks are read from a
     * pre-pass snapshot so receiver updates cannot race source reads. */
#pragma omp parallel for schedule(static)
    for (int i = 0; i < NC; i++) { th1s[i] = th1[i]; th2s[i] = th2[i]; }
#pragma omp parallel for schedule(static)
    for (int l = 0; l < NL; l++) {
        if (lA[l] <= 0) continue;
        int i = li[l], j = lj[l];
        for (int c = 0; c < 2; c++) {
            for (int dir = 0; dir < 2; dir++) {
                int s = SLOT(l, c, dir);
                double f = lwant[s];
                if (f <= 0) continue;
                int src = dir == 0 ? i : j;
                f *= c == 0 ? scl0[src] : scl1[src];
                /* S1/U2: transport within a mode is continuous — the
                 * action atom lives at mode boundaries (conversions),
                 * not on same-mode deposits. (The earlier transport
                 * quantization froze few-quantum flow; see HBAR.md §6.) */
                lwant[s] = f;               /* resolved deposit */
                if (lem[s] <= 0) lph[s] = 0;
                lem[s] += f;
                lflux[2 * l + c] += f;
            }
        }
    }
#pragma omp parallel for schedule(static)
    for (int i = 0; i < NC; i++) {
        double d0 = 0, d1 = 0;
        for (int q = cls[i]; q < cls[i + 1]; q++) {
            int l = clidx[q];
            if (lA[l] <= 0) continue;
            int dir = li[l] == i ? 0 : 1;
            d0 += lwant[SLOT(l, 0, dir)] > 0 ? lwant[SLOT(l, 0, dir)] : 0;
            d1 += lwant[SLOT(l, 1, dir)] > 0 ? lwant[SLOT(l, 1, dir)] : 0;
        }
        Ee[i] -= d0;
        Em[i] -= d1;
    }
#pragma omp parallel for schedule(static)
    for (int i = 0; i < NC; i++) {
        for (int q = cls[i]; q < cls[i + 1]; q++) {
            int l = clidx[q];
            if (lA[l] <= 0) continue;
            int dir = li[l] == i ? 1 : 0;   /* the direction ARRIVING at i */
            int src = dir == 0 ? li[l] : lj[l];
            double d = ld[l];
            for (int c = 0; c < 2; c++) {
                double f = lwant[SLOT(l, c, dir)];
                if (f <= 0) continue;
                double mobr = c == 0 ? Ee[i] : Em[i];
                double wsrc = c == 0 ? w1e[src] : w2e[src];
                double thsrc = c == 0 ? th1s[src] : th2s[src];
                double *thr = c == 0 ? &th1[i] : &th2[i];
                double ms = 1.0, mr = 1.0;
                if (c == 1) {
                    ms = dir == 0 ? lq[l] : lp[l];
                    mr = dir == 0 ? lp[l] : lq[l];
                }
                double err = wrap_pi(ms * (thsrc - wsrc * d / P.C) - mr * *thr) / mr;
                double mix = f / (f + mobr + lockf);
                *thr += P.kappa_lock * mix * err;
            }
        }
    }

    /* pass 5: advance cycles; deliver only on completion (rate C, gate).
     * Gathered by receiver: each slot's deliveries are handled by the
     * cell that receives them (slot ownership is unique), foreign
     * clocks from the snapshot, rough fires deferred to the serial
     * QATOM flush in pass 6. */
#pragma omp parallel for schedule(static)
    for (int i2 = 0; i2 < NC; i2++) { th1s[i2] = th1[i2]; th2s[i2] = th2[i2]; }
#pragma omp parallel for schedule(static) reduction(+:rough_total,backs_total)
    for (int recv = 0; recv < NC; recv++) {
        for (int q = cls[recv]; q < cls[recv + 1]; q++) {
            int l = clidx[q];
            if (lA[l] <= 0) continue;   /* pinched channel: flight holds, conserved */
            int i = li[l], j = lj[l];
            int dir = li[l] == recv ? 1 : 0;   /* the direction ARRIVING here */
            int send = dir == 0 ? i : j;
            double adv = dt * P.C / ld[l];
            for (int c = 0; c < 2; c++) {
                int s = SLOT(l, c, dir);
                if (lem[s] <= 0) continue;
                lph[s] += adv;
                if (lph[s] < 1.0) continue;
                double freec = P.cap - (Em[recv] + Ee[recv]);
                double take = lem[s];
                if (take > freec) take = freec > 0 ? freec : 0;
                if (take > 0) {
                    double mobprev = c == 0 ? Ee[recv] : Em[recv];
                    lem[s] -= take;
                    if (c == 0) {
                        Ee[recv] += take;
                    } else {
                        /* C2: dissonance radiates (Tartini / Plomp-Levelt).
                         * The rough part of an off-comb delivery lands as
                         * field, not dense; a just interval is not rough.
                         * Energy leaving the dense mode un-converts its
                         * space share (matter is converted space — the
                         * pattern's space returns with it, as in evap).
                         * Paired ledger move: take = dense+field+space. */
                        double det = lq[l] * w2e[i] - lp[l] * w2e[j];
                        double R = 2.0 * fabs(det) * P.gamma_rough
                                   / (det * det + P.gamma_rough * P.gamma_rough);
                        double rough = take * P.rough_k * R;
                        /* U2: roughness is a D->F conversion — it fires in
                         * whole atoms at this channel's completed flight
                         * cycle, through the receiver's D->F account. */
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
                    }
                    /* completion reinforces the join: lock receiver's clock */
                    double wsend = c == 0 ? w1e[send] : w2e[send];
                    double thsend = c == 0 ? th1s[send] : th2s[send];
                    double *thr = c == 0 ? &th1[recv] : &th2[recv];
                    double ms = 1.0, mr = 1.0;
                    if (c == 1) {
                        ms = dir == 0 ? lq[l] : lp[l];
                        mr = dir == 0 ? lp[l] : lq[l];
                    }
                    double err = wrap_pi(ms * (thsend - wsend * ld[l] / P.C) - mr * *thr) / mr;
                    double mix = take / (take + mobprev + lockf);
                    *thr += P.kappa_lock * mix * err;
                }
                if (lem[s] <= 1e-17) { lem[s] = 0; lph[s] = 0; }
                else if (take <= 0) lph[s] = 0;      /* full receiver: wait a cycle */
                else lph[s] -= 1.0;                  /* residual rides again */
            }
        }
    }

    /* pass F: the repaired field sector — unitary amplitude dynamics.
     * Onsite rotation at w1e (load + wall detuning), then exact pairwise
     * hop rotations over live channels (each conserves the pair norm to
     * roundoff), then instrument drains and kicks. Superposition is
     * native: anti-phased arrivals cancel in components while the norm
     * redistributes to the constructive loci — interference as consonance
     * geometry, with energy never destroyed. */
#pragma omp parallel for schedule(static)
    for (int i = 0; i < NC; i++) {
        double ang = w1e[i] * dt;
        double cc, ss;
        sincos(ang, &ss, &cc);
        double a1 = fa1[i], a2 = fa2[i];
        fa1[i] = cc * a1 + ss * a2;
        fa2[i] = -ss * a1 + cc * a2;
        if (pol_on) {
            double b1 = fb1[i], b2 = fb2[i];
            fb1[i] = cc * b1 + ss * b2;
            fb2[i] = -ss * b1 + cc * b2;
        }
    }
    /* a cell's total joining bandwidth is a property of the cell, not of
     * its accidental contact geometry: symmetric normalization of the hop
     * weights removes the foam's diagonal disorder (the dominant wave
     * scatterer), leaving only the geometric (off-diagonal) texture */
#pragma omp parallel for schedule(static)
    for (int i = 0; i < NC; i++) {
        double s = 0;
        for (int q = cls[i]; q < cls[i + 1]; q++) {
            int l = clidx[q];
            if (lA[l] <= 0) continue;
            s += (lA[l] / Aref) * (dref / ld[l]);
        }
        fsum[i] = s;
    }
    /* hops in edge-color batches: within a color no two links share a
     * cell, so the pair rotations apply conflict-free; colors run in
     * fixed order (thread-count-independent result) */
    for (int cb = 0; cb < ncolors; cb++) {
#pragma omp parallel for schedule(static)
        for (int q = colstart[cb]; q < colstart[cb + 1]; q++) {
            int l = colidx[q];
            if (lA[l] <= 0) continue;
            int i = li[l], j = lj[l];
            double si = fsum[i], sj = fsum[j];
            if (si <= 1e-12 || sj <= 1e-12) continue;
            /* hop U = exp(-i tau X). With the seed convention theta = -k*x
             * (down-path clocks lag, kappa = -k), this pairing propagates the
             * packet along +k: v_g = +t*sum d*sin(k*d). Verified by centroid. */
            double w = (lA[l] / Aref) * (dref / ld[l]);
            double tau = P.field_J * w / sqrt(si * sj) * dt;
            double cc, ss;
            sincos(tau, &ss, &cc);
            double a1i = fa1[i], a2i = fa2[i], a1j = fa1[j], a2j = fa2[j];
            fa1[i] = cc * a1i + ss * a2j;
            fa2[i] = cc * a2i - ss * a1j;
            fa1[j] = cc * a1j + ss * a2i;
            fa2[j] = cc * a2j - ss * a1i;
            if (pol_on) {
                double b1i = fb1[i], b2i = fb2[i], b1j = fb1[j], b2j = fb2[j];
                fb1[i] = cc * b1i + ss * b2j;
                fb2[i] = cc * b2i - ss * b1j;
                fb1[j] = cc * b1j + ss * b2i;
                fb2[j] = cc * b2j - ss * b1i;
            }
        }
    }
    for (int i = 0; i < NC; i++) {
        if (!cflag[i]) continue;
        if (cflag[i] >= 2 && cflag[i] <= 3) {
            /* condensation medium: arriving field locks into record */
            double eH = fa1[i] * fa1[i] + fa2[i] * fa2[i];
            double eV = pol_on ? fb1[i] * fb1[i] + fb2[i] * fb2[i] : 0;
            if (eH + eV > 0) {
                double sf = exp(-P.absorb_rate * dt);
                double drained = (eH + eV) * (1.0 - sf * sf);
                fa1[i] *= sf; fa2[i] *= sf;
                if (pol_on) { fb1[i] *= sf; fb2[i] *= sf; }
                Em[i] += drained;   /* the conserved ledger: total energy */
                if (cflag[i] == 2 && (P.t_expose <= 0 || sim_t <= P.t_expose)) {
                    /* analyzer-basis decomposition of the record (tier 2/3):
                     * native H/V, or a rotated basis; t_choice flips the
                     * basis mid-run (delayed choice). Diagnostic split —
                     * exA + exB = drained exactly. */
                    double dA = drained, dB = 0;
                    if (pol_on) {
                        double adeg = P.analyzer_deg;
                        if (P.t_choice > 0 && sim_t < P.t_choice) adeg = -999;
                        if (adeg > -900) {
                            double ca = cos(adeg * M_PI / 180.0);
                            double sa = sin(adeg * M_PI / 180.0);
                            double c1 = ca * fa1[i] + sa * fb1[i];
                            double c2 = ca * fa2[i] + sa * fb2[i];
                            double ep = c1 * c1 + c2 * c2;
                            double frac = (eH + eV) > 0 ? ep / (eH + eV) * (1.0/(sf*sf)) : 0;
                            if (frac > 1) frac = 1;
                            dA = drained * frac;
                            dB = drained - dA;
                        } else {
                            dA = eH + eV > 0 ? drained * eH / (eH + eV) : 0;
                            dB = drained - dA;
                        }
                    }
                    exA[i] += dA;
                    exB[i] += dB;
                    /* U2: the click grain is the action atom of the
                     * absorbing mode, eps(w1e) — a knob (e_click) only in
                     * the continuous legacy limit. */
                    double grain = A0eff > 0 ? A0eff * w1e[i] / TWO_PI
                                             : P.e_click;
                    if (grain > 0) {
                        int n = (int)(Em[i] / grain);
                        for (int qn = clickn[i]; qn < n; qn++)
                            printf("# CLICK t=%.2f x=%.2f y=%.2f z=%.2f\n",
                                   sim_t, cx[i], cy[i], cz[i]);
                        if (n > clickn[i]) clickn[i] = n;
                    }
                }
            }
        } else if (cflag[i] == 4) {
            if (pol_on && P.tag_rate > 0) {
                /* tier 2 — unitary which-path tag: rotate H->V in the
                 * chirality pair. Lossless marking: no energy taken, no
                 * phase noise added; the path is written into the tag */
                double tr = P.tag_rate * dt;
                double ct = cos(tr), st = sin(tr);
                double h1 = fa1[i], h2 = fa2[i], v1 = fb1[i], v2 = fb2[i];
                fa1[i] = ct * h1 - st * v1;
                fa2[i] = ct * h2 - st * v2;
                fb1[i] = st * h1 + ct * v1;
                fb2[i] = st * h2 + ct * v2;
            }
            if (P.sigma_det > 0) {
                /* kick recorder (round-2 style back-action) */
                double kick = P.sigma_det * sqrt(dt) * grand();
                double cc = cos(kick), ss = sin(kick);
                double a1 = fa1[i], a2 = fa2[i];
                fa1[i] = cc * a1 + ss * a2;
                fa2[i] = -ss * a1 + cc * a2;
            }
            if (P.det_tap > 0) {
                double e = fa1[i] * fa1[i] + fa2[i] * fa2[i]
                         + (pol_on ? fb1[i] * fb1[i] + fb2[i] * fb2[i] : 0);
                double sf = exp(-P.det_tap * dt);
                fa1[i] *= sf; fa2[i] *= sf;
                if (pol_on) { fb1[i] *= sf; fb2[i] *= sf; }
                Em[i] += e * (1.0 - sf * sf);
            }
        }
    }
#pragma omp parallel for schedule(static)
    for (int i = 0; i < NC; i++) {
        Ee[i] = fa1[i] * fa1[i] + fa2[i] * fa2[i]
              + (pol_on ? fb1[i] * fb1[i] + fb2[i] * fb2[i] : 0);
        /* th1 is a diagnostic cache at runtime (the gated c==0 transport
         * is gone since the field repair): refresh only on diag steps */
        if (want_th1 && Ee[i] > 1e-20) th1[i] = atan2(fa2[i], fa1[i]);
    }

    /* pass 6: dense clock + beat-gated conversion (complete cycles only).
     * SERIAL by design: the QATOM diagnostic counter and print stream
     * stay deterministic. Deferred rough fires from pass 5 flush here
     * first, in cell order (aggregated per cell per step — sums of
     * whole atoms remain whole multiples on the eps(w) grid). */
    for (int i = 0; i < NC; i++) {
        if (roughq[i] > 0) {
            qatom_diag(0, atoms_w(w2e[i], w1e[i]), roughq[i]);
            roughq[i] = 0;
        }
    }
    for (int i = 0; i < NC; i++) {
        th2[i] = fmod(th2[i] + w2e[i] * dt, TWO_PI);
        if (cflag[i]) continue;   /* instruments do not beat-convert */
        cbeta[i] += (w1e[i] - w2e[i]) * dt;
        int beat_fire = 0;
        if (cbeta[i] >= TWO_PI) { cbeta[i] -= TWO_PI; beat_fire = 1; }
        else if (cbeta[i] <= -TWO_PI) { cbeta[i] += TWO_PI; beat_fire = 1; }
        if (beat_fire) {
            /* condensation F -> D at the cell's completed beat cycle,
             * consuming the cell's own space. U2: the conversion moves
             * whole action atoms only — the threshold structure
             * (photoelectric) IS the atom, not a tuned e_cond. */
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
                    fa1[i] *= fac;
                    fa2[i] *= fac;
                    if (pol_on) { fb1[i] *= fac; fb2[i] *= fac; }
                    Ee[i] -= d1;
                    Es[i] -= dsp;
                    Em[i] += d1 + dsp;
                }
            }
            /* over-full: harmonics refuse the load, D -> F,S — the same
             * uniform law through the cell's D->F account */
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

    /* pass 7: plane re-alignment (flux-weighted, headless) + tumble.
     * Two pulls per plane: toward the partners' plane orientation, and
     * toward containing the flux direction (n perpendicular to u), so the
     * planes' intersection line rotates onto the direction of transfer. */
    double sq = P.sigma_tumble * sqrt(dt);
    /* the gaussian stream is drawn SERIALLY in the exact legacy order
     * (cell-major, plane, xyz) so the tumble noise is bit-identical to
     * the single-thread kernel; the apply loop then runs parallel,
     * reading neighbor normals from a snapshot (no read/write race). */
    if (sq > 0)
        for (int i = 0; i < NC; i++)
            for (int k = 0; k < 6; k++) rngbuf[6 * (size_t)i + k] = grand();
#pragma omp parallel for schedule(static)
    for (int i = 0; i < NC; i++) {
        nsnap[6 * (size_t)i + 0] = n1x[i];
        nsnap[6 * (size_t)i + 1] = n1y[i];
        nsnap[6 * (size_t)i + 2] = n1z[i];
        nsnap[6 * (size_t)i + 3] = n2x[i];
        nsnap[6 * (size_t)i + 4] = n2y[i];
        nsnap[6 * (size_t)i + 5] = n2z[i];
    }
#pragma omp parallel for schedule(static)
    for (int i = 0; i < NC; i++) {
        for (int c = 0; c < 2; c++) {
            double mx = c == 0 ? n1x[i] : n2x[i];
            double my = c == 0 ? n1y[i] : n2y[i];
            double mz = c == 0 ? n1z[i] : n2z[i];
            double ax = 0, ay = 0, az = 0, fsum = 0;
            for (int q = cls[i]; q < cls[i + 1]; q++) {
                int l = clidx[q];
                double fl = lflux[2 * l + c];
                if (fl <= 0) continue;
                int o = li[l] == i ? lj[l] : li[l];
                const double *no = nsnap + 6 * (size_t)o + 3 * c;
                double ox = no[0];
                double oy = no[1];
                double oz = no[2];
                double sgn = (mx * ox + my * oy + mz * oz) >= 0 ? 1.0 : -1.0;
                double du = mx * lux[l] + my * luy[l] + mz * luz[l];
                ax += fl * (sgn * ox - du * lux[l]);
                ay += fl * (sgn * oy - du * luy[l]);
                az += fl * (sgn * oz - du * luz[l]);
                fsum += fl;
            }
            double vx = mx, vy = my, vz = mz;
            if (fsum > 0) {
                double w = P.kappa_align * dt / fsum;
                vx += w * ax; vy += w * ay; vz += w * az;
            }
            if (sq > 0) {
                const double *g = rngbuf + 6 * (size_t)i + 3 * c;
                vx += sq * g[0]; vy += sq * g[1]; vz += sq * g[2];
            }
            double nn = sqrt(vx * vx + vy * vy + vz * vz);
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
/* diagnostics                                                         */
/* ------------------------------------------------------------------ */

static double ksum(const double *a, int n)
{
    double s = 0, comp = 0;
    for (int i = 0; i < n; i++) {
        double y = a[i] - comp, t = s + y;
        comp = (t - s) - y;
        s = t;
    }
    return s;
}

static void snapshot(int step)
{
    char path[512];
    snprintf(path, sizeof(path), "%s/cells_%06d.tsv", P.snap_dir, step);
    FILE *f = fopen(path, "w");
    if (!f) { fprintf(stderr, "# WARN cannot write %s\n", path); return; }
    fprintf(f, "id\tx\ty\tz\tr\tEs\tEm\tEe\tth1\tth2\tw1e\tw2e\tn1x\tn1y\tn1z\tn2x\tn2y\tn2z\n");
    for (int i = 0; i < NC; i++)
        fprintf(f, "%d\t%.4f\t%.4f\t%.4f\t%.4f\t%.6g\t%.6g\t%.6g\t%.4f\t%.4f\t%.5f\t%.5f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n",
                i, cx[i], cy[i], cz[i], cr[i], Es[i], Em[i], Ee[i], th1[i], th2[i],
                w1e[i], w2e[i], n1x[i], n1y[i], n1z[i], n2x[i], n2y[i], n2z[i]);
    fclose(f);
}

/* MASS instrument — the lump census: connected components of free dense
 * energy above lump_thr, walked over live channels. The campaign's
 * primary eye: per-lump mass, centroid, rms radius (top 10 by mass). */
static int *lv_vis = NULL, *lv_q = NULL;
static void lumps_row(double t)
{
    if (!lv_vis) {
        lv_vis = malloc(NC * sizeof(int));
        lv_q = malloc(NC * sizeof(int));
    }
    for (int i = 0; i < NC; i++) lv_vis[i] = 0;
    double thr = P.lump_thr;
    double emfree = 0;
    for (int i = 0; i < NC; i++) if (!cflag[i]) emfree += Em[i];
    double lm[10], lx[10], ly[10], lz[10], lr[10];
    int nl = 0, ntot = 0;
    for (int i = 0; i < NC; i++) {
        if (lv_vis[i] || cflag[i] || Em[i] <= thr) continue;
        int h = 0, tq = 0;
        lv_q[tq++] = i; lv_vis[i] = 1;
        double m = 0, mx = 0, my = 0, mz = 0;
        while (h < tq) {
            int u = lv_q[h++];
            m += Em[u];
            mx += Em[u] * cx[u]; my += Em[u] * cy[u]; mz += Em[u] * cz[u];
            for (int q = cls[u]; q < cls[u + 1]; q++) {
                int l = clidx[q];
                if (lA[l] <= 0) continue;
                int o = li[l] == u ? lj[l] : li[l];
                if (!lv_vis[o] && !cflag[o] && Em[o] > thr) {
                    lv_vis[o] = 1; lv_q[tq++] = o;
                }
            }
        }
        ntot++;
        mx /= m; my /= m; mz /= m;
        double rg = 0;
        for (int q = 0; q < tq; q++) {
            int u = lv_q[q];
            double dx = cx[u] - mx, dy = cy[u] - my, dz = cz[u] - mz;
            rg += Em[u] * (dx * dx + dy * dy + dz * dz);
        }
        rg = sqrt(rg / m);
        int at = nl < 10 ? nl : -1;
        if (at < 0) {
            for (int k = 0; k < 10; k++)
                if (m > lm[k] && (at < 0 || lm[k] < lm[at])) at = k;
        }
        if (at >= 0) {
            if (nl < 10) nl++;
            lm[at] = m; lx[at] = mx; ly[at] = my; lz[at] = mz; lr[at] = rg;
        }
    }
    printf("# LUMP t=%.2f n=%d Emfree=%.6g", t, ntot, emfree);
    for (int k = 0; k < nl; k++) {
        int b = 0;
        for (int k2 = 1; k2 < nl; k2++) if (lm[k2] > lm[b]) b = k2;
        printf(" | m=%.6g x=%.2f y=%.2f z=%.2f rg=%.2f",
               lm[b], lx[b], ly[b], lz[b], lr[b]);
        lm[b] = -1;
    }
    printf("\n");
}

static void diag_row(double t)
{
    double tEs = ksum(Es, NC), tEm = ksum(Em, NC), tEe = ksum(Ee, NC);
    double tET = ksum(lem, 4 * NL);
    double tot = tEs + tEm + tEe + tET;
    double drift = E0_total != 0 ? (tot - E0_total) / E0_total : 0;
    if (!isfinite(tot)) { fprintf(stderr, "# FATAL non-finite energy at t=%g\n", t); exit(2); }

    /* dense centroid + containment over FREE cells only — instrument
     * cells (walls, screens, sinks) hold recorded ledger, not moving
     * matter (positions: diagnostic reconstruction only) */
    double cmx = 0, cmy = 0, cmz = 0, cin = 0, cmE = 0;
    double rcont = 2.5 * P.sigma;
    for (int i = 0; i < NC; i++) {
        if (cflag[i]) continue;
        cmE += Em[i];
        cmx += Em[i] * cx[i]; cmy += Em[i] * cy[i]; cmz += Em[i] * cz[i];
        double dx = cx[i] - cenx, dy = cy[i] - ceny, dz = cz[i] - cenz;
        if (dx * dx + dy * dy + dz * dz < rcont * rcont) cin += Em[i];
    }
    if (cmE > 1e-12) {
        cmx /= cmE; cmy /= cmE; cmz /= cmE; cin /= cmE;
    } else {
        cmx = cmy = cmz = cin = 0;
    }

    /* field front radius (threshold on seed amplitude) */
    double front = 0;
    double thr = 0.05 * (P.amp > 0 ? P.amp : 1.0);
    for (int i = 0; i < NC; i++) {
        if (Ee[i] > thr) {
            double dx = cx[i] - cenx, dy = cy[i] - ceny, dz = cz[i] - cenz;
            double rr = sqrt(dx * dx + dy * dy + dz * dz);
            if (rr > front) front = rr;
        }
    }

    /* contact-area defect + live channels */
    double defA = 0;
    int nact = 0;
    for (int l = 0; l < NL; l++) {
        defA += lA0[l] - lA[l];
        if (lA[l] > 0) nact++;
    }

    /* core vs far shells: radius + tick rate (curvature / time dilation) */
    double rin = 0, rout = 0, win = 0, wout = 0;
    int nin = 0, nout = 0;
    double s1 = P.sigma, s4 = 4.0 * P.sigma;
    for (int i = 0; i < NC; i++) {
        double dx = cx[i] - cenx, dy = cy[i] - ceny, dz = cz[i] - cenz;
        double rr = sqrt(dx * dx + dy * dy + dz * dz);
        if (rr < s1) { rin += cr[i]; win += w1e[i]; nin++; }
        else if (rr > s4) { rout += cr[i]; wout += w1e[i]; nout++; }
    }
    if (nin) { rin /= nin; win /= nin; }
    if (nout) { rout /= nout; wout /= nout; }

    /* momentum instruments (P-battery): field-energy centroid, and the
     * flux moment — the first moment of in-flight energy over link
     * directions. Positions/directions are diagnostic reconstruction
     * only; the laws never see them. */
    double cex = 0, cey = 0, cez = 0;
    if (tEe > 1e-12) {
        for (int i = 0; i < NC; i++) {
            cex += Ee[i] * cx[i]; cey += Ee[i] * cy[i]; cez += Ee[i] * cz[i];
        }
        cex /= tEe; cey /= tEe; cez /= tEe;
    }
    double pfx = 0, pfy = 0, pfz = 0;
    for (int l = 0; l < NL; l++) {
        double imb = (lem[SLOT(l, 0, 0)] + lem[SLOT(l, 1, 0)])
                   - (lem[SLOT(l, 0, 1)] + lem[SLOT(l, 1, 1)]);
        if (imb == 0) continue;
        pfx += imb * lux[l]; pfy += imb * luy[l]; pfz += imb * luz[l];
    }
    /* field-sector momentum current: the repaired field moves by unitary
     * hops, so its momentum lives in the amplitude phase gradient, not
     * in flight slots. Per link: J = 2*J_hop*Im[psi_i^* psi_j] (both
     * chirality components), summed over link directions. */
    double ffx = 0, ffy = 0, ffz = 0;
    {
        double Aref_d = M_PI * P.r0 * P.r0, dref_d = 2.0 * P.r0;
        for (int l = 0; l < NL; l++) {
            if (lA[l] <= 0) continue;
            int i = li[l], j = lj[l];
            double w = P.field_J * (lA[l] / Aref_d) * (dref_d / ld[l]);
            /* sign pairs with the hop generator U=exp(-i tau X) and the
             * seed convention theta=-k*x (DOUBLESLIT round 2): forward
             * current along +u is MINUS this Im-product. Calibrated
             * against the ballistic p1 packet. */
            double cur = fa1[i] * fa2[j] - fa2[i] * fa1[j]
                       + fb1[i] * fb2[j] - fb2[i] * fb1[j];
            double J = -2.0 * w * cur;
            ffx += J * lux[l]; ffy += J * luy[l]; ffz += J * luz[l];
        }
    }

    printf("%.3f\t%.9g\t%.9g\t%.9g\t%.9g\t%.12g\t%.3e\t%d\t%.3f\t%.3f\t%.3f\t%.3f\t%.4f\t%.4f\t%.5g\t%.4f\t%.4f\t%.5f\t%.5f\t%.4f\t%.4f\t%.4f\t%.5g\t%.5g\t%.5g\t%.5g\t%.5g\t%.5g\n",
           t, tEs, tEm, tEe, tET, tot, drift, nact,
           cmx, cmy, cmz, cin, front, rin, defA, rout, rin > 0 ? rin / (rout > 0 ? rout : 1) : 0,
           win, wout, cex, cey, cez, pfx, pfy, pfz, ffx, ffy, ffz);

    printf("# CONV t=%.2f cond=%.6g evap=%.6g rough=%.6g back_s=%.6g\n",
           t, cond_total, evap_total, rough_total, backs_total);
    if (P.lump_diag) lumps_row(t);

    /* G4 — radial throughput report: per shell, the accumulated space-
     * transport flux rate (positive = outward), the instantaneous field
     * radial current (the relativistic channel: evaporation leaves as
     * field at c), the space mean and the free dense sum. Separating
     * these is the point: a leaking blob's radiative wind must not be
     * mistaken for a static gravitational monopole. */
    if (P.rad_diag && P.s_k > 0) {
        double dr = 0.5 * P.L / 8.0;
        double fr[8], esh[8], emh[8];
        int nsh[8];
        for (int k = 0; k < 8; k++) { fr[k] = 0; esh[k] = 0; emh[k] = 0; nsh[k] = 0; }
        double Aref_d = M_PI * P.r0 * P.r0, dref_d = 2.0 * P.r0;
        for (int l = 0; l < NL; l++) {
            if (lA[l] <= 0) continue;
            int i = li[l], j = lj[l];
            double mx = 0.5 * (cx[i] + cx[j]) - cenx;
            double my = 0.5 * (cy[i] + cy[j]) - ceny;
            double mz = 0.5 * (cz[i] + cz[j]) - cenz;
            double rr = sqrt(mx * mx + my * my + mz * mz);
            int k = (int)(rr / dr);
            if (rr < 1e-9 || k >= 8) continue;
            double w = P.field_J * (lA[l] / Aref_d) * (dref_d / ld[l]);
            double cur = fa1[i] * fa2[j] - fa2[i] * fa1[j]
                       + fb1[i] * fb2[j] - fb2[i] * fb1[j];
            double J = -2.0 * w * cur;
            fr[k] += J * (lux[l] * mx + luy[l] * my + luz[l] * mz) / rr;
        }
        for (int i = 0; i < NC; i++) {
            if (cflag[i]) continue;
            double dx = cx[i] - cenx, dy = cy[i] - ceny, dz = cz[i] - cenz;
            int k = (int)(sqrt(dx * dx + dy * dy + dz * dz) / dr);
            if (k < 8) { esh[k] += Es[i]; emh[k] += Em[i]; nsh[k]++; }
        }
        double win = t - srad_t0;
        printf("# RAD t=%.2f", t);
        for (int k = 0; k < 8; k++)
            printf(" %.2f:%.4g,%.4g,%.4f,%.4g",
                   (k + 0.5) * dr, win > 0 ? srad[k] / win : 0.0, fr[k],
                   nsh[k] ? esh[k] / nsh[k] : 0, emh[k]);
        printf("\n");
        for (int k = 0; k < 8; k++) srad[k] = 0;
        srad_t0 = t;
    }

    if (nsamp < nsamp_max) {
        ds_t[nsamp] = t; ds_em[nsamp] = tEm; ds_front[nsamp] = front;
        ds_cmx[nsamp] = cmx; ds_cmy[nsamp] = cmy; ds_cmz[nsamp] = cmz;
        ds_defA[nsamp] = defA;
        nsamp++;
    }
}

/* gate/flux statistics over the energetic region per sector, against the
 * +x axis: where is the throttle — axis, resonance, tail gate, headroom? */
static void debug_row(double t)
{
    double gmw = P.gamma_res_m < 0 ? P.gamma_res : P.gamma_res_m;
    for (int c = 0; c < 2; c++) {
        double gp = 0, gn = 0, axs = 0, rs = 0, hd = 0, n = 0;
        double fp = 0, fn = 0, ep = 0, en = 0;
        double *EB = c == 0 ? Ee : Em;
        double *TH = c == 0 ? th1 : th2;
        double *WE = c == 0 ? w1e : w2e;
        double G2 = c == 0 ? P.gamma_res * P.gamma_res : gmw * gmw;
        for (int l = 0; l < NL; l++) {
            if (lA[l] <= 0) continue;
            int i = li[l], j = lj[l];
            double d = ld[l];
            int fwd_is_ij = lux[l] >= 0;    /* slot moving toward +x */
            double fij = lwant[SLOT(l, c, 0)] * (c == 0 ? scl0[i] : scl1[i]);
            double fji = lwant[SLOT(l, c, 1)] * (c == 0 ? scl0[j] : scl1[j]);
            fp += fwd_is_ij ? fij : fji;
            fn += fwd_is_ij ? fji : fij;
            ep += lem[SLOT(l, c, fwd_is_ij ? 0 : 1)];
            en += lem[SLOT(l, c, fwd_is_ij ? 1 : 0)];
            if (EB[i] + EB[j] < 0.02) continue;
            double d1i = n1x[i] * lux[l] + n1y[i] * luy[l] + n1z[i] * luz[l];
            double d2i = n2x[i] * lux[l] + n2y[i] * luy[l] + n2z[i] * luz[l];
            double d1j = n1x[j] * lux[l] + n1y[j] * luy[l] + n1z[j] * luz[l];
            double d2j = n2x[j] * lux[l] + n2y[j] * luy[l] + n2z[j] * luz[l];
            double axi = (1 - d1i * d1i) * (1 - d2i * d2i); if (axi < 0) axi = 0;
            double axj = (1 - d1j * d1j) * (1 - d2j * d2j); if (axj < 0) axj = 0;
            double dw = WE[i] - WE[j];
            double g_ij = gate_of(wrap_pi(TH[i] - WE[i] * d / P.C - TH[j]));
            double g_ji = gate_of(wrap_pi(TH[j] - WE[j] * d / P.C - TH[i]));
            gp += fwd_is_ij ? g_ij : g_ji;
            gn += fwd_is_ij ? g_ji : g_ij;
            axs += axi * axj;
            rs += G2 / (G2 + dw * dw);
            int r = fwd_is_ij ? j : i;
            double h = 1.0 - (Em[r] + Ee[r]) / P.cap;
            hd += h < 0 ? 0 : (h > 1 ? 1 : h);
            n += 1;
        }
        if (n < 1) n = 1;
        printf("# DBG%d t=%.2f links=%d gate+x=%.3f gate-x=%.3f axij=%.3f res=%.3f head+x=%.3f dep+x=%.4g dep-x=%.4g em+x=%.4g em-x=%.4g\n",
               c, t, (int)n, gp / n, gn / n, axs / n, rs / n, hd / n, fp, fn, ep, en);
    }
}

/* standing-pair diagnostics: per pair, the rung offset delta (the comma),
 * lock error (should temper to delta/2), two-way gate product (the tongue),
 * mutual-flight fraction (the shelter), retention, action per joint cycle. */
static void pair_report(double t, int final)
{
    if (NP <= 0) return;
    double s_gg[3] = { 0, 0, 0 }, s_ret[3] = { 0, 0, 0 }, s_fl[3] = { 0, 0, 0 };
    int s_n[3] = { 0, 0, 0 };
    for (int p = 0; p < NP; p++) {
        int i = ppi[p], j = ppj[p], l = ppl[p];
        double d = ld[l];
        double bp = lp[l], bq = lq[l];
        double wsum = bq * w2e[i] + bp * w2e[j];
        double delta = wrap_pi(wsum * d / P.C);
        double dth = wrap_pi(bq * th2[i] - bq * w2e[i] * d / P.C - bp * th2[j]);
        double gg = gate_of(dth)
                    * gate_of(wrap_pi(bp * th2[j] - bp * w2e[j] * d / P.C - bq * th2[i]));
        double fl = lem[SLOT(l, 1, 0)] + lem[SLOT(l, 1, 1)];
        double ret = pE0[p] > 0 ? (Em[i] + Em[j] + fl) / pE0[p] : 0;
        double flf = pE0[p] > 0 ? fl / pE0[p] : 0;
        double act = fl * 2.0 * d / P.C;
        double xm = 0.5 * (Em[i] + flload[i] + Em[j] + flload[j]) / P.cap;   /* pitch load (U3) */
        double ratio = w2e[j] > 1e-12 ? w2e[i] / w2e[j] : 0;
        double shed = pE0[p] - (Em[i] + Em[j] + fl);
        double flq = A0eff > 0
                     ? fl / (A0eff * 0.5 * (w2e[i] + w2e[j]) / TWO_PI)
                     : -1.0;
        printf("# PAIR t=%.1f p=%d d=%.3f pq=%d:%d delta=%+.3f x=%.3f ratio=%.4f dth=%+.3f gg=%.4f fl=%.3f ret=%.3f shed=%.3f act=%.4f flq=%.2f\n",
               t, p, d, (int)bp, (int)bq, delta, xm, ratio, dth, gg, flf, ret, shed, act, flq);
        int bin = fabs(delta) < 0.6 ? 0 : (fabs(delta) < 1.2 ? 1 : 2);
        s_gg[bin] += gg; s_ret[bin] += ret; s_fl[bin] += flf; s_n[bin]++;
    }
    if (final) {
        printf("# RESULT pair_tongue |delta|<0.6: n=%d gg=%.3f ret=%.3f fl=%.3f | 0.6-1.2: n=%d gg=%.3f ret=%.3f fl=%.3f | >1.2: n=%d gg=%.3f ret=%.3f fl=%.3f\n",
               s_n[0], s_n[0] ? s_gg[0] / s_n[0] : 0, s_n[0] ? s_ret[0] / s_n[0] : 0, s_n[0] ? s_fl[0] / s_n[0] : 0,
               s_n[1], s_n[1] ? s_gg[1] / s_n[1] : 0, s_n[1] ? s_ret[1] / s_n[1] : 0, s_n[1] ? s_fl[1] / s_n[1] : 0,
               s_n[2], s_n[2] ? s_gg[2] / s_n[2] : 0, s_n[2] ? s_ret[2] / s_n[2] : 0, s_n[2] ? s_fl[2] / s_n[2] : 0);
    }
}

/* mode=ladder: the harmonic limit, computed before it is measured.
 * Equal-pair rung m: w_eff*d/C = pi*m with w_eff = w2/(1+q*x)
 *   -> tuning curve x*(d,m) = (w2*d/(pi*m*C) - 1)/q
 * Tongue (entrained equilibrium splits the comma delta evenly):
 *   G^2(delta) = ((1+cos(delta/2))/2)^(2*p_gate)
 * Strangulation: channel dies at d > 2*r0*(1 - s_pull*x*cap/e_s0)^(1/3). */
static void run_ladder(void)
{
    double p2 = 2.0 * P.p_gate;
    double dhalf = 2.0 * acos(2.0 * pow(2.0, -1.0 / p2) - 1.0);
    printf("# ladder mode — computed harmonic limit (no dynamics run)\n");
    printf("# dense sector: w2=%g q=%g cap=%g p_gate=%g s_pull=%g C=%g\n",
           P.w2, P.q_detune, P.cap, P.p_gate, P.s_pull, P.C);
    printf("# tongue half-width (half-max of G^2): delta_half=%.3f rad\n", dhalf);
    for (int m = 1; m <= 4; m++) {
        double d_lo = M_PI * m * P.C / P.w2;                       /* x*=0 */
        double d_hi = M_PI * m * P.C * (1.0 + P.q_detune * 0.95) / P.w2; /* x*=0.95 */
        printf("# rung m=%d: reachable d=[%.3f, %.3f]  (x* from 0 to 0.95)\n", m, d_lo, d_hi);
    }
    double dcut = 2.0 * P.r0 * cbrt(1.0 - P.s_pull * 0.4 * P.cap / P.e_s0);
    printf("# strangulation at x=0.4: channel dies beyond d=%.3f\n", dcut);
    double gmr = P.gamma_res_m < 0 ? P.gamma_res : P.gamma_res_m;
    printf("# interval vocabulary at q_detune=%g (comb_limit=%d, reference x_i=0.10):\n",
           P.q_detune, P.comb_limit);
    for (int t = 0; t < ncomb; t++) {
        int pp = combp[t], qq = combq[t];
        if (pp < qq) continue;               /* each ratio >= 1 reported once */
        double wi = P.w2 / (1.0 + P.q_detune * 0.10);
        double wj = wi * (double)qq / (double)pp;
        double xj = (P.w2 / wj - 1.0) / P.q_detune;
        int ok = xj >= 0.0 && xj <= 0.95;
        printf("# %d:%d  partner_x=%.3f  %s  tongue_width=%.4f",
               pp, qq, xj, ok ? "REACHABLE" : "unreachable", gmr / (pp * qq));
        if (ok) {
            double gsum = qq * wi + pp * wj;
            for (int m = 1; m <= 3; m++) {
                double dr = TWO_PI * m / gsum;
                if (dr >= 0.9 && dr <= 1.8) printf("  d*(m=%d)=%.3f", m, dr);
            }
        }
        printf("\n");
    }
    printf("# d\tx*(d,m=1)\tdelta(d;x0=%.2f)\tG2_pred\n", P.pair_x0 > 0 ? P.pair_x0 : 0.4);
    double x0 = P.pair_x0 > 0 ? P.pair_x0 : 0.4;
    double we = P.w2 / (1.0 + P.q_detune * x0);
    for (double d = P.pair_dlo; d <= P.pair_dhi + 1e-9; d += 0.05) {
        double xs = (P.w2 * d / (M_PI * P.C) - 1.0) / P.q_detune;
        double delta = wrap_pi(2.0 * we * d / P.C);
        double g2 = pow(0.5 * (1.0 + cos(delta / 2.0)), p2);
        printf("%.3f\t%.3f\t%+.3f\t%.4f\n", d, xs, delta, g2);
    }
    printf("# RESULT ladder m_max=1 for nearest-neighbour links (2*pi*2/(w_sum_max*d_max) > 1)\n");
}

/* tier 4 — HOM coincidence registry. The two-quantum object is one
 * process with two ends: for non-interacting quanta the joint amplitude
 * stays the (anti)symmetrized product of the two single-quantum fields,
 * so coincidences are computed exactly from phiA, phiB at the detectors:
 *   boson:   |phiA(1)phiB(2) + phiB(1)phiA(2)|^2 / 2
 *   fermion: |phiA(1)phiB(2) - phiB(1)phiA(2)|^2 / 2
 *   distinguishable baseline: |phiA(1)|^2|phiB(2)|^2 + |phiB(1)|^2|phiA(2)|^2 */
static double hom_sym = 0, hom_asym = 0, hom_dist = 0;
static double hom_p1 = 0, hom_p2 = 0;   /* port populations (solo calibration) */
static int *homD1 = NULL, *homD2 = NULL;
static int nD1 = 0, nD2 = 0;

static void hom_accum(double dtint)
{
    if (P.init != 6) return;
    if (!homD1) {
        double yw1 = 0.5 * P.L - 3.0, yw2 = 0.5 * P.L + 3.0, hw = 1.5;
        double xlo = P.L - 1.5 - P.hom_det, xhi = P.L - 1.6;
        homD1 = malloc(NC * sizeof(int));
        homD2 = malloc(NC * sizeof(int));
        for (int i = 0; i < NC; i++) {
            if (cflag[i] != 0 || cx[i] < xlo || cx[i] > xhi) continue;
            if (fabs(cy[i] - yw1) < hw) homD1[nD1++] = i;
            else if (fabs(cy[i] - yw2) < hw) homD2[nD2++] = i;
        }
        printf("# hom detectors: D1=%d cells D2=%d cells x=[%.1f,%.1f]\n",
               nD1, nD2, xlo, xhi);
    }
    for (int a = 0; a < nD1; a++) {
        int i = homD1[a];
        double A1r = fa1[i], A1i = fa2[i], B1r = fb1[i], B1i = fb2[i];
        hom_p1 += dtint * (A1r * A1r + A1i * A1i);
        for (int b = 0; b < nD2; b++) {
            int j = homD2[b];
            double A2r = fa1[j], A2i = fa2[j], B2r = fb1[j], B2i = fb2[j];
            /* u = phiA(1)*phiB(2), v = phiB(1)*phiA(2), complex */
            double ur = A1r * B2r - A1i * B2i, ui = A1r * B2i + A1i * B2r;
            double vr = B1r * A2r - B1i * A2i, vi = B1r * A2i + B1i * A2r;
            double sr = ur + vr, si = ui + vi;
            double ar = ur - vr, ai = ui - vi;
            hom_sym += dtint * 0.5 * (sr * sr + si * si);
            hom_asym += dtint * 0.5 * (ar * ar + ai * ai);
            hom_dist += dtint * ((A1r * A1r + A1i * A1i) * (B2r * B2r + B2i * B2i)
                                 + (B1r * B1r + B1i * B1i) * (A2r * A2r + A2i * A2i));
        }
    }
    for (int b = 0; b < nD2; b++) {
        int j = homD2[b];
        hom_p2 += dtint * (fa1[j] * fa1[j] + fa2[j] * fa2[j]);
    }
}

static void hom_report(void)
{
    if (P.init != 6) return;
    double split = (hom_p1 + hom_p2) > 0 ? hom_p2 / (hom_p1 + hom_p2) : 0;
    printf("# RESULT hom coupler_split(A cross->D2)=%.3f\n", split);
    if (hom_dist > 0)
        printf("# RESULT hom g_boson=%.4f g_fermion=%.4f (1=distinguishable, 0=full dip, 2=full peak) dist=%.4g\n",
               hom_sym / hom_dist, hom_asym / hom_dist, hom_dist);
}

/* tier 1 — record one interval of the screen absorption process */
static void record_qinc(double t)
{
    if (qints >= QINT_MAX || nscr <= 0) return;
    for (int s = 0; s < nscr; s++) {
        int i = scr_id[s];
        double inc = Em[i] - scr_prev[s];
        if (inc < 0) inc = 0;
        qinc[(size_t)qints * nscr + s] = (float)inc;
        scr_prev[s] = Em[i];
    }
    qtms[qints++] = t;
}

/* tier 1 — sample n_quanta single-quantum claims from the recorded
 * absorption process. Exact for non-interacting quanta in the linear
 * regime: identical unitary evolution per quantum, hazard = absorb law,
 * survival = the wave's own norm depletion. One claim per quantum —
 * indivisibility is CONSTRUCTION R3 (delivery is atomic), not a fit. */
static void sample_quanta(void)
{
    if (P.n_quanta <= 0 || nscr <= 0 || qints <= 0) return;
    double tot = 0;
    for (size_t q = 0; q < (size_t)qints * nscr; q++) tot += qinc[q];
    if (tot <= 0) { printf("# QCLICK none: no screen absorption recorded\n"); return; }
    enum { NB = 36 };
    int qhist[NB];
    memset(qhist, 0, sizeof(qhist));
    double ylo = P.sink_m, yhi = P.L - P.sink_m, bw = (yhi - ylo) / NB;
    printf("# claim rule: %d quanta, each claimed whole by the first completed cycle\n",
           P.n_quanta);
    for (int n = 0; n < P.n_quanta; n++) {
        double r = frand() * tot, acc = 0;
        int hit_int = qints - 1, hit_s = nscr - 1;
        for (size_t q = 0; q < (size_t)qints * nscr; q++) {
            acc += qinc[q];
            if (acc >= r) { hit_int = (int)(q / nscr); hit_s = (int)(q % nscr); break; }
        }
        int i = scr_id[hit_s];
        printf("# QCLICK q=%d t=%.1f y=%.2f z=%.2f\n",
               n, qtms[hit_int], cy[i], cz[i]);
        int b = (int)((cy[i] - ylo) / bw);
        if (b >= 0 && b < NB) qhist[b]++;
    }
    double sb = 0, sd = 0;
    int nbr = 0, nda = 0;
    double D = P.screen_x - P.wall_x;
    double yA = 0.5 * P.L - 0.5 * P.slit_sep, yB = 0.5 * P.L + 0.5 * P.slit_sep;
    double kmag = sqrt(P.kx * P.kx + P.ky * P.ky + P.kz * P.kz);
    if (kmag < 1e-9) kmag = P.w1 / P.C;
    printf("# QHIST y_center clicks predicted_cos2\n");
    for (int b = 0; b < NB; b++) {
        double yc = ylo + (b + 0.5) * bw;
        double rA = sqrt(D * D + (yc - yA) * (yc - yA));
        double rB = sqrt(D * D + (yc - yB) * (yc - yB));
        double cph = cos(kmag * (rA - rB) / 2.0);
        double pred = cph * cph;
        printf("# QHIST %.2f\t%d\t%.3f\n", yc, qhist[b], pred);
        if (fabs(yc - 0.5 * P.L) < 5.5) {
            if (pred > 0.85) { sb += qhist[b]; nbr++; }
            else if (pred < 0.15) { sd += qhist[b]; nda++; }
        }
    }
    double mb = nbr ? sb / nbr : 0, md = nda ? sd / nda : 0;
    printf("# RESULT qclick_visibility_central=%.3f centre=%.1f minima=%.1f (discrete events)\n",
           (mb + md) > 0 ? (mb - md) / (mb + md) : 0, mb, md);
}

static void linfit(const double *x, const double *y, int n, double *m, double *b, double *r2)
{
    double sx = 0, sy = 0, sxx = 0, sxy = 0;
    for (int i = 0; i < n; i++) { sx += x[i]; sy += y[i]; sxx += x[i] * x[i]; sxy += x[i] * y[i]; }
    double den = n * sxx - sx * sx;
    *m = den != 0 ? (n * sxy - sx * sy) / den : 0;
    *b = n ? (sy - *m * sx) / n : 0;
    double ssr = 0, sst = 0, ym = n ? sy / n : 0;
    for (int i = 0; i < n; i++) {
        double e = y[i] - (*m * x[i] + *b);
        ssr += e * e;
        sst += (y[i] - ym) * (y[i] - ym);
    }
    *r2 = sst > 0 ? 1.0 - ssr / sst : 0;
}

static void final_report(void)
{
    double tEs = ksum(Es, NC), tEm = ksum(Em, NC), tEe = ksum(Ee, NC), tET = ksum(lem, 4 * NL);
    double tot = tEs + tEm + tEe + tET;
    /* exit face (G-battery): the +x absorbing face is a recorder — the
     * transverse centroid of what it recorded is the transmitted beam's
     * exit position (lensing: compare against the no-mass baseline) */
    if (P.edge_sink > 0) {
        double fm = P.edge_sink + 0.8;
        double fE = 0, fy = 0, fz = 0;
        for (int i = 0; i < NC; i++) {
            if (cflag[i] != 3) continue;
            if (cx[i] < P.L - fm) continue;
            if (cy[i] < fm || cy[i] > P.L - fm) continue;
            if (cz[i] < fm || cz[i] > P.L - fm) continue;
            fE += Em[i];
            fy += Em[i] * cy[i];
            fz += Em[i] * cz[i];
        }
        if (fE > 1e-12) { fy /= fE; fz /= fE; }
        printf("# RESULT exit_face E=%.5g y=%.4f z=%.4f\n", fE, fy, fz);
    }

    /* space profile (G-battery): shell means of Es around the seed
     * center — the gravitational footprint, if any law maintains one */
    {
        double shE[8]; int shN[8];
        double dr = 0.5 * P.L / 8.0;
        for (int k = 0; k < 8; k++) { shE[k] = 0; shN[k] = 0; }
        for (int i = 0; i < NC; i++) {
            if (cflag[i]) continue;
            double dx = cx[i] - cenx, dy = cy[i] - ceny, dz = cz[i] - cenz;
            int k = (int)(sqrt(dx * dx + dy * dy + dz * dz) / dr);
            if (k < 8) { shE[k] += Es[i]; shN[k]++; }
        }
        printf("# RESULT es_shells");
        for (int k = 0; k < 8; k++)
            printf(" %.2f:%.4f", (k + 0.5) * dr, shN[k] ? shE[k] / shN[k] : 0);
        printf("\n");
    }
    printf("# RESULT conservation E0=%.12g Efinal=%.12g rel_drift=%.3e\n",
           E0_total, tot, E0_total != 0 ? (tot - E0_total) / E0_total : 0);

    printf("# RESULT roughness_radiated=%.6g (dissonant dense deliveries converted to field, conserved)\n",
           rough_total);

    int h = nsamp / 2;
    double kk = sqrt(P.kx * P.kx + P.ky * P.ky + P.kz * P.kz);

    if (P.init == 2 && kk > 0 && nsamp - h >= 3) {
        double m, b, r2;
        linfit(ds_t + h, ds_front + h, nsamp - h, &m, &b, &r2);
        printf("# RESULT front_speed v=%.4f v_over_C=%.4f r2=%.4f (foam tortuosity keeps v below C)\n",
               m, m / P.C, r2);
    }
    if (P.init == 3 && nsamp - h >= 3) {
        double mx, my, mz, b, r2x, r2y, r2z;
        linfit(ds_t + h, ds_cmx + h, nsamp - h, &mx, &b, &r2x);
        linfit(ds_t + h, ds_cmy + h, nsamp - h, &my, &b, &r2y);
        linfit(ds_t + h, ds_cmz + h, nsamp - h, &mz, &b, &r2z);
        double sp = sqrt(mx * mx + my * my + mz * mz);
        double ck = 0;
        if (kk > 0 && sp > 0) ck = (mx * P.kx + my * P.ky + mz * P.kz) / (sp * kk);
        printf("# RESULT blob_drift v=(%.5f,%.5f,%.5f) speed=%.5f cos_to_kdir=%.4f\n",
               mx, my, mz, sp, ck);
    }
    if (nsamp >= 4) {
        double vm = 0, mm = 0;
        for (int i = 0; i < nsamp; i++) mm += ds_em[i];
        mm /= nsamp;
        for (int i = 0; i < nsamp; i++) vm += (ds_em[i] - mm) * (ds_em[i] - mm);
        if (vm / nsamp > 1e-8) {
            double m, b, r2;
            linfit(ds_em, ds_defA, nsamp, &m, &b, &r2);
            printf("# RESULT curvature_fit defA_per_Em=%.5f r2=%.4f (linear defect ~ converted energy)\n", m, r2);
        }
    }

    if (P.init == 5) {
        /* screen exposure vs y (z-integrated), against the parameter-free
         * two-source consonance map cos^2(w*(rA-rB)/2C) */
        enum { NB = 36 };
        double hist[NB];
        memset(hist, 0, sizeof(hist));
        double ylo = P.sink_m, yhi = P.L - P.sink_m, bw = (yhi - ylo) / NB;
        double e_scr = 0, e_recA = 0, e_recB = 0;
        for (int i = 0; i < NC; i++) {
            if (cflag[i] == 2) {
                double ex = expose_frz ? expose_frz[i] : Em[i];
                int b = (int)((cy[i] - ylo) / bw);
                if (b >= 0 && b < NB) hist[b] += ex;
                e_scr += ex;
            } else if (cflag[i] == 4) {
                if (cy[i] < 0.5 * P.L) e_recA += Em[i]; else e_recB += Em[i];
            }
        }
        double D = P.screen_x - P.wall_x;
        double yA = 0.5 * P.L - 0.5 * P.slit_sep;
        double yB = 0.5 * P.L + 0.5 * P.slit_sep;
        double sb = 0, sd = 0, sbc = 0, sdc = 0;
        int nbr = 0, nda = 0, nbc = 0, ndc = 0;
        /* fringe loci from the SEEDED wavevector: lambda = 2*pi/|k| */
        double kmag = sqrt(P.kx * P.kx + P.ky * P.ky + P.kz * P.kz);
        if (kmag < 1e-9) kmag = P.w1 / P.C;
        printf("# SCREEN y_center exposure predicted_cos2 (lambda=%.3f)\n",
               TWO_PI / kmag);
        for (int b = 0; b < NB; b++) {
            double yc = ylo + (b + 0.5) * bw;
            double rA = sqrt(D * D + (yc - yA) * (yc - yA));
            double rB = sqrt(D * D + (yc - yB) * (yc - yB));
            double cph = cos(kmag * (rA - rB) / 2.0);
            double pred = cph * cph;
            printf("# SCREEN %.2f\t%.4f\t%.3f\n", yc, hist[b], pred);
            if (fabs(yc - 0.5 * P.L) < 13.0) {
                if (pred > 0.85) { sb += hist[b]; nbr++; }
                else if (pred < 0.15) { sd += hist[b]; nda++; }
            }
            /* central template: centre maximum vs first minima — the
             * discriminator a single-slit hump cannot fake */
            if (fabs(yc - 0.5 * P.L) < 5.5) {
                if (pred > 0.85) { sbc += hist[b]; nbc++; }
                else if (pred < 0.15) { sdc += hist[b]; ndc++; }
            }
        }
        double mb = nbr ? sb / nbr : 0, md = nda ? sd / nda : 0;
        double vis = (mb + md) > 0 ? (mb - md) / (mb + md) : 0;
        double mbc = nbc ? sbc / nbc : 0, mdc = ndc ? sdc / ndc : 0;
        double visc = (mbc + mdc) > 0 ? (mbc - mdc) / (mbc + mdc) : 0;
        printf("# RESULT fringe_visibility_central=%.3f centre_max=%.3f (n=%d) first_minima=%.3f (n=%d)\n",
               visc, mbc, nbc, mdc, ndc);
        printf("# RESULT fringe_visibility_wide=%.3f bright_mean=%.3f (n=%d) dark_mean=%.3f (n=%d)\n",
               vis, mb, nbr, md, nda);
        printf("# RESULT screen_exposure=%.4f whichpath_record A=%.4f B=%.4f\n",
               e_scr, e_recA, e_recB);

        if (pol_on) {
            /* tier 2/3 — analyzer-basis histograms and per-ledger
             * visibility (eraser: fringes in one, anti-fringes in the
             * other, flat total) */
            double hA[NB], hB[NB];
            memset(hA, 0, sizeof(hA));
            memset(hB, 0, sizeof(hB));
            double tA = 0, tB = 0;
            for (int i = 0; i < NC; i++) {
                if (cflag[i] != 2) continue;
                int b = (int)((cy[i] - ylo) / bw);
                if (b >= 0 && b < NB) { hA[b] += exA[i]; hB[b] += exB[i]; }
                tA += exA[i]; tB += exB[i];
            }
            double sbA = 0, sdA = 0, sbB = 0, sdB = 0;
            int nb2 = 0, nd2 = 0;
            printf("# SCREENAB y_center exA exB pred\n");
            for (int b = 0; b < NB; b++) {
                double yc = ylo + (b + 0.5) * bw;
                double rA2 = sqrt(D * D + (yc - yA) * (yc - yA));
                double rB2 = sqrt(D * D + (yc - yB) * (yc - yB));
                double cph = cos(kmag * (rA2 - rB2) / 2.0);
                double pred = cph * cph;
                printf("# SCREENAB %.2f\t%.4f\t%.4f\t%.3f\n", yc, hA[b], hB[b], pred);
                if (fabs(yc - 0.5 * P.L) < 5.5) {
                    if (pred > 0.85) { sbA += hA[b]; sbB += hB[b]; nb2++; }
                    else if (pred < 0.15) { sdA += hA[b]; sdB += hB[b]; nd2++; }
                }
            }
            double vA = 0, vB = 0;
            if (nb2 && nd2) {
                double a1 = sbA / nb2, a0 = sdA / nd2, b1 = sbB / nb2, b0 = sdB / nd2;
                vA = (a1 + a0) > 0 ? (a1 - a0) / (a1 + a0) : 0;
                vB = (b1 + b0) > 0 ? (b1 - b0) / (b1 + b0) : 0;
            }
            printf("# RESULT analyzer ledgers: totalA=%.3f totalB=%.3f V_A=%.3f V_B=%.3f tagfrac=%.4f\n",
                   tA, tB, vA, vB, (tA + tB) > 0 ? tB / (tA + tB) : 0);
        }
        sample_quanta();
    }
}

/* ------------------------------------------------------------------ */
/* bell mode: joint unfinished harmonic vs LHV control                 */
/* ------------------------------------------------------------------ */

static void run_bell(void)
{
    double A[2] = { P.aA1 * M_PI / 180.0, P.aA2 * M_PI / 180.0 };
    double B[2] = { P.aB1 * M_PI / 180.0, P.aB2 * M_PI / 180.0 };
    double Ej[2][2], El[2][2];

    printf("# bell mode: analyzers pass with the kernel's plane-overlap law cos^2(delta)\n");
    printf("# joint: pair = ONE unfinished conversion; first completed cycle collapses the shared phase\n");
    printf("# lhv:   same cos^2 responder, but a private pre-assigned phase per packet\n");
    printf("# aA\taB\tE_joint_mc\tE_joint_th\tE_lhv_mc\tE_lhv_th\n");

    for (int ia = 0; ia < 2; ia++) {
        for (int ib = 0; ib < 2; ib++) {
            double a = A[ia], b = B[ib];
            long acc_j = 0, acc_l = 0;
            for (long n = 0; n < P.trials; n++) {
                /* joint harmonic: unpolarized until the first cycle completes */
                int Ap = frand() < 0.5;
                double lam = Ap ? a : a + M_PI / 2.0;
                double cb = cos(b - lam);
                int Bp = frand() < cb * cb;
                acc_j += (Ap == Bp) ? 1 : -1;
                /* LHV control: private phase fixed at the source */
                double lam0 = frand() * M_PI;
                double ca2 = cos(a - lam0), cb2 = cos(b - lam0);
                int Ap2 = frand() < ca2 * ca2;
                int Bp2 = frand() < cb2 * cb2;
                acc_l += (Ap2 == Bp2) ? 1 : -1;
            }
            Ej[ia][ib] = (double)acc_j / (double)P.trials;
            El[ia][ib] = (double)acc_l / (double)P.trials;
            printf("%.1f\t%.1f\t%+.4f\t%+.4f\t%+.4f\t%+.4f\n",
                   A[ia] * 180 / M_PI, B[ib] * 180 / M_PI,
                   Ej[ia][ib], cos(2.0 * (a - b)),
                   El[ia][ib], 0.5 * cos(2.0 * (a - b)));
        }
    }
    double Sj = Ej[0][0] + Ej[1][0] + Ej[1][1] - Ej[0][1];
    double Sl = El[0][0] + El[1][0] + El[1][1] - El[0][1];
    printf("# RESULT S_joint=%.4f S_lhv=%.4f classical_bound=2 tsirelson=%.4f trials=%ld\n",
           Sj, Sl, 2.0 * sqrt(2.0), P.trials);
    printf("# reading: the violation is bought by the joint harmonic (one process, two ends),\n");
    printf("# not by the response law — the identical cos^2 responder with private phases stays below 2.\n");
}

/* ------------------------------------------------------------------ */

int main(int argc, char **argv)
{
    cfg_defaults();
    for (int a = 1; a < argc; a++) {
        char *eq = strchr(argv[a], '=');
        if (eq) {
            *eq = 0;
            set_kv(argv[a], eq + 1);
        } else {
            load_cfg(argv[a]);
        }
    }
    rng_s = P.seed ? P.seed : 88172645463325252ULL;
    for (int q = 0; q < 8; q++) xrand();

    comb_build();
    print_cfg();
    if (P.mode == 1) { run_bell(); return 0; }
    if (P.mode == 2) { run_ladder(); return 0; }

    build_field();

    int NS = (int)(P.T / P.dt + 0.5);
    nsamp_max = NS / (P.diag_every > 0 ? P.diag_every : 1) + 8;
    ds_t = malloc(nsamp_max * sizeof(double));
    ds_em = malloc(nsamp_max * sizeof(double));
    ds_front = malloc(nsamp_max * sizeof(double));
    ds_cmx = malloc(nsamp_max * sizeof(double));
    ds_cmy = malloc(nsamp_max * sizeof(double));
    ds_cmz = malloc(nsamp_max * sizeof(double));
    ds_defA = malloc(nsamp_max * sizeof(double));

    E0_total = ksum(Es, NC) + ksum(Em, NC) + ksum(Ee, NC) + ksum(lem, 4 * NL);

    printf("# t\tE_space\tE_dense\tE_field\tE_transfer\tE_total\trel_drift\tlive_ch\tcm_x\tcm_y\tcm_z\tcontain\tfront_r\tr_core\tdefA\tr_far\tr_ratio\tw_core\tw_far\n");
    diag_row(0.0);
    pair_report(0.0, 0);
    if (P.snap_every > 0) snapshot(0);

    for (int s = 1; s <= NS; s++) {
        sim_t = s * P.dt;
        want_th1 = (P.diag_every > 0 && s % P.diag_every == 0)
                || (P.snap_every > 0 && s % P.snap_every == 0) || s == NS;
        step_field();
        if (P.init == 5 && P.t_expose > 0 && !expose_frz && sim_t >= P.t_expose) {
            /* the shutter closes: the record so far is the photograph */
            expose_frz = malloc(NC * sizeof(double));
            for (int i = 0; i < NC; i++) expose_frz[i] = cflag[i] == 2 ? Em[i] : 0;
            printf("# shutter closed at t=%.1f\n", sim_t);
        }
        if (P.diag_every > 0 && s % P.diag_every == 0) {
            diag_row(s * P.dt);
            if (P.debug) debug_row(s * P.dt);
            pair_report(s * P.dt, 0);
            if (P.init == 5 && P.n_quanta > 0
                && (P.t_expose <= 0 || sim_t <= P.t_expose))
                record_qinc(sim_t);
            if (P.init == 6) hom_accum(P.diag_every * P.dt);
        }
        if (P.snap_every > 0 && s % P.snap_every == 0) snapshot(s);
    }
    final_report();
    pair_report(P.T, 1);
    hom_report();
    return 0;
}
