/* freecell.c — v91 KERNEL OF RECORD (carried from v90/kernel/freecell.c,
 * itself carried from v89/freecell.c; both originals stay frozen as
 * evidence).
 *
 * v91 role: the production C kernel for the RADIANCE + COHERENT-CHANNEL
 * programme (v91/README.md is the charter; v91/THEORY.md the law spec).
 * ONE law-candidate term is added over v90, GUARDED so that the default
 * table is laws_V2g byte-exactly: the graded sub-cap RADIANCE channel
 * (k_rad=0 default = V2g; see "v91 RADIANCE" in pass 6). Everything
 * else — apparatus, meters, experiments — is the v90 surface verbatim.
 * The Go kernel (v91/fab) remains the A/B experiment.
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
#include <zstd.h>
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
    /* --- v91 LAW CANDIDATE (laws_V3r): graded sub-cap radiance.
     * FORGE E1 measured V2g's radiance as a cap-wall step (outtake
     * exactly 0 below cap) — no interior fixed point, no stable mass.
     * Candidate A: at each beat, demand k_rad*cap*x^p_rad (x = Em/cap),
     * times det (rad_clock=0: compensates the beat-period stretch so
     * the PER-TIME rate follows x^p_rad; rad_clock=1: rides the raw
     * slowing beat). Grain-quantized through the atoms machinery,
     * routed like evaporation. k_rad=0 => V2g byte-exactly. --- */
    double k_rad, p_rad;
    int    rad_clock;
    /* --- v91 LAW CANDIDATE B (CANTUS, coherent channel): the
     * superimposed harmonic-lock field (CANTUS.md; user-directed
     * 2026-08-04, "atoms are NOT cells"). Per-cell order parameter
     * ca (lives on the cells' own bonds, no background, no energy)
     * + holdings memory cxl. Pass H: (a) Kuramoto lock on the matter
     * CLOCKS th2 correcting the differential ladder-closure error
     * (never the wants — the anti-kappa_reac decision C-D1); (b)
     * within-mode retuning current on holdings driven by the
     * memory-deviation difference (pairwise-conserving; flavor-
     * preserving). k_cant=0 && k_tune=0 => byte-identical step. --- */
    double k_cant, k_tune, cant_tau;
    int    cant_seed;
    int    cant_grow;   /* 1 = self-growth (law candidate); 0 = seeded-only (INSTRUMENT: honest-medium probe, CANTUS.md §3.4) */
    /* --- v91 EXCHANGE REGISTRY (REGISTRY.md; user-opened 2026-08-05).
     * Identity-carrying transfers at slot grain: a per-slot ledger of
     * reciprocal DELIVERIES (the slot is the continuous identity pair
     * — it dies with either endpoint). reg_tau=0 => no stamps, no
     * prints, byte-identical step. reg_gate=1 gates the cantus gauge
     * growth by the registry match r (bond-vs-churn = identity, the
     * CANTUS-measured requirement). --- */
    double reg_tau;     /* delivery-ledger memory; 0 = registry OFF   */
    int    reg_gate;    /* cantus growth target *= r_s. 0 off;
                         * 1 = F-B  rho * gross/(gross+f0)  (f0=0 => F-A);
                         * 2 = F-D  s/(s+f0), s = 2*min = reciprocal flow */
    double reg_f0;      /* the flow half-saturation constant (units E/t) */
    /* --- v91 IDENTITY lane (IDENTITY.md; user-opened 2026-08-06).
     * Parcel-carried ONTOLOGICAL identity: an episode gid born with
     * the matter (x crosses par_hi), carried by its dense-flight
     * parcels (depositor label on the slot in-flight), registered at
     * arrival, dying with the matter (x below par_lo). par_tau=0 =>
     * no gid state, no stamps, no prints, byte-identical step.
     * par_gate=1 gates the cantus gauge growth target by the
     * maturity-clocked identity continuity r_id — stamp AGE is the
     * gate variable the lock cannot manufacture (it has no force
     * before arming). --- */
    double par_tau;     /* identity-ledger memory; 0 = lane OFF       */
    int    par_gate;    /* cantus growth target *= r_id; 0 off        */
    int    par_form;    /* 0 = I-A binary; 1 = I-B flow-shaped        */
    double par_lo;      /* episode retire threshold on x (hysteresis) */
    double par_hi;      /* episode mint threshold on x                */
    double par_mature;  /* stamp maturity time (anti-ignition clock)  */
    /* --- v91 WORKFN lane (WORKFN.md; user task list 2026-08-06).
     * Emergent conversion threshold: the condensation demand lives in
     * BOUND MATTER, not space (B7 fix). wf_on=0 => the econd line
     * reads exactly as today (byte-inert). wf_on=1: form W-A,
     * econd_i = Em[i] >= wf_floor ? 0 : wf_far — a PRESENCE
     * threshold (W-M1: foam Em == 0 exactly; matter > 0). --- */
    /* --- QUENCH-2 apparatus (QUENCH.md §6; seed-time only) ---
     * conf_r > 0 carves a circular pinned-fixture shell (cavity) of
     * radius conf_r about box centre with a leak gap; spin_m adds an
     * azimuthal phase winding to the seeded packet. Both inert at 0. */
    double conf_r, conf_gap, conf_th, conf_pinw;
    int    spin_m;
    double seed_mw;   /* v93 item 4: seed MATTER azimuthal winding m*phi on the blob (no field/door) */
    double imp_k;       /* QUENCH-2 fix: inward radial phase tilt on the
                         * seeded packet — implosion focusing = the
                         * LAWFUL whole-box compress-release (seeded
                         * energy only; conservative). 0 = off. */
    double qp_phase;    /* QUENCH-3 (QUENCH.md §7): the conversion door
                         * writes the wave's phase into the matter clock
                         * it feeds — th2 pulled toward arg(fa) by the
                         * field-derived fraction of the new holdings.
                         * No energy content; 0 = byte-inert. */
    /* --- HORIZON apparatus (HORIZON.md; assertion instrument): the
     * forced hole = the prover's NEW_LEDGER_X + ROUTE structure by
     * fiat. Cells inside bh_r of box centre are eaten each step into
     * the uncapped, pitchless, pi-invisible accumulator Eh (in the
     * conservation sum and nowhere else); their clocks freeze (no
     * beat, no door, no emission — I2 makes one-way structural).
     * bh_r=0 = byte-inert. --- */
    double bh_r, bh_k;
    double bh_sep;      /* FLOW: two hole centres at cx ± sep/2 (0 = one) */
    /* --- FLOW apparatus (FLOW.md; the ASYM.md bed-digging channel
     * law, v1 = space channel only): per-slot conductance weight
     * sbed grown by |lowpassed signed net flow| under a zero-sum
     * per-cell budget (anti-ignition structure); slot-borne, mortal,
     * clamped [0.2,5]; no energy content. bed_k=0 = byte-inert. --- */
    double bed_k, bed_tau;
    int    wf_on;       /* 0 = lane OFF (byte-inert)                  */
    double wf_floor;    /* presence floor (W-M1-selected 0.01)        */
    double wf_far;      /* empty-space demand (optics-grade 99)       */
    /* --- v91 AMPLITUDE lane Phase M (AMPLITUDE.md; user task #24).
     * SHADOW AMPLITUDE meter: a complex amplitude rides the existing
     * dense flows (deposit sqrt(f)e^{i m theta_src}, transit rotation
     * -m w dt, delivery composed into the dst chart frame) and moves
     * NOTHING. The meter is the coherence-deficit map rho_coh =
     * |sum dA| / sum |dA| per link/class. amp_tau=0 => byte-inert. */
    double amp_tau;     /* shadow window; 0 = lane OFF                */
    /* --- v92 AMPLITUDE Phase L lane L-1 (L0_DESIGN.md §3): the shadow
     * PROMOTED FROM METER TO DRIVER. Where the shadow amplitude composes
     * coherently on a chord (m>=amp_mmin) link, it biases the dense want
     * toward the coherent direction — "translation IS the current"
     * (momentum as the first moment of conversion). amp_drv=0 => the
     * shadow stays a Phase-M meter and the step is byte-inert. */
    double amp_drv;     /* L-1 amplitude-driven transport strength; 0 = OFF */
    double amp_mmin;    /* L-1 chart-order gate floor (2 = fifths+, not unison) */
    /* --- v93 THE UNITARY DENSE CHANNEL (v93/README.md PART II): the dense
     * sector gets the field sector's transport algebra. Within-mode dense
     * transport becomes a product of UNITARY PAIRWISE PLANE ROTATIONS (pass
     * F's cousin) on the dense amplitude psi_m = sqrt(Em) e^{i th2},
     * replacing the additive magnitude want. Each Givens hop conserves the
     * two-cell norm exactly (conservation is a theorem of the update, not a
     * patched ledger); the cross term 2 Im(psi_i* psi_j) that the additive
     * Em-ledger rejected IS the link current = where momentum lives.
     * amp_nat=0 => the additive want path runs unchanged (byte-inert vs the
     * 87-bar V3a surface); >0 => the unitary dense hop engages, the additive
     * want+inflight is bypassed. The door (pass 6) is NEVER unitarized. */
    double amp_nat;     /* unitary dense channel strength; 0 = additive want (byte-inert) */
    double amp_logate;  /* linearize probe: >0 drops the phase-dependent gate from tau_s */
    /* --- v93 arg(psi) door (v93/README.md §II.7): at a condensation event
     * the field click (sqrt(d1) at phase atan2(fa2,fa1)) composes COHERENTLY
     * with the existing matter amplitude psi_m, and the result is
     * renormalized to the conserved energy -- arg(psi_m_new) is the coherent-
     * sum direction (carries the field phase + the interference cross-term =
     * the current), |psi_m| fixed by Em. Replaces qp_phase's partial cell-
     * clock pull: the fired atom carries arg(psi_m), not m*th2 written
     * piecemeal (IV.6). amp_door=0 => byte-inert (qp_phase path unchanged). */
    double amp_door;    /* arg(psi) coherent-amplitude door; 0 = qp_phase path (byte-inert) */
    int    qatom_every;   /* apparatus (print-only): QATOM sampler period */
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
    int    tri2_kind;               /* second triangle's kind; -1 = follow tri_kind
                                     * (mixed-species pair: 0=UUD,1=UDD)  */
    double oct_x, oct_doff;
    double shear_eps, shear_t;      /* instantaneous deviatoric strain test */
    /* DS — double slit on the free substrate (v90 P1; DS.md). Wall is
     * CARVED VACUUM (no cells = no contact = no transport); 2D bath. */
    double slit_wallx, slit_th, slit_sep, slit_hw;
    double slit_screenx, slit_srcx, slit_sy;
    double slit_t0, slit_t1;        /* screen gate window */
    double slit_pinw;               /* fixture half-width about the wall plane */
    /* COMPOSITE (COMPOSITE.md): rings as many-celled quarks.
     * kind 0 = unison pair, 1 = UUD, 2 = UDD. Internal spacing per ring:
     * the m=1 unison rung if it fits inside contact (molecule), else
     * 0.98*contact (droplet) — the U/D structural asymmetry (T6). */
    int    rings_kind, rings_nv, rings_branch;
    double rings_xU, rings_xD, rings_gapoff;
    /* MOTION #28: two-blob collision (exp=blob2) */
    double blob2_sep;               /* COM separation along x */
    double blob2_kx;                /* opposing tilts +-kx toward each other */
    /* MOTION #29: embedded object in the slit beam (slit_obj=1: dense
     * blob at box centre; measure attenuation vs slit_obj=0) */
    int    slit_obj;
    double obj_amp, obj_sigma;
    /* DS tier 1: condensing screen (slit_clicks=1: screen-strip cells
     * become condensation-active; conversions there are CLICKS) */
    int    slit_clicks;
    /* MOTION: all-modes center-of-energy meter (p1_meter=1) — cumulative
     * first-moment FLOW per channel (space/flight/field/geometry).
     * Momentum is the first moment of conversion (v89 center-of-energy
     * theorem); flow bookkeeping is torus-safe: every term is a local
     * minimal-image displacement times the energy that crossed it.
     * Pure ledger — no feedback on the dynamics; default OFF so every
     * standing log stays byte-identical. */
    int    p1_meter;
    /* XSEC (MOTION.md): angular cross-section apparatus. sect_meter=1
     * enables an annular per-sector exposure meter (time-integrated
     * Ee*dt per angular bin), centred on (sect_x,sect_y) (<0 = box
     * centre), r in [sect_r0,sect_r1), gated [sect_t0,sect_t1]; sector
     * 0 is centred on +x (downbeam). obj_y (<0 = box centre) moves the
     * slit_obj occulter in y — the impact parameter. Pure ledger,
     * default OFF: every standing log stays byte-identical. */
    int    sect_meter, sect_n;
    double sect_r0, sect_r1, sect_x, sect_y, sect_t0, sect_t1;
    double obj_y;
    /* FORGE apparatus: convtag=1 forces the tag-split ledger print for
     * any experiment (slit_obj still implies it); grad_r0/grad_frac
     * thin the dart bath inside radius r0 of the box centre (a
     * matter-free foam-density structure); tag_r tags cells within
     * radius of centre at init WITHOUT loading them (region ledger). */
    int    convtag;
    double grad_r0, grad_frac, tag_r;
    double rings_xcomp;             /* CQ3b: whole-ring load reduction, U-flavor rings */
    double rings_xcompD;            /* CQ3b: whole-ring load shift, D-flavor rings */
    double rings_sep;               /* kind 3: nucleon COM separation */
    int    slit_mask;               /* 0 both, 1 A only, 2 B only, 3 no wall */
    char snap_dir[256];
    char snap_file[256];            /* FCS v3 single-stream output */
    int  snap_comp;                 /* 1 = compress CELL/LINK chunks (CMPD) */
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
    P.k_rad = 0; P.p_rad = 4; P.rad_clock = 0;
    P.k_cant = 0; P.k_tune = 0; P.cant_tau = 50; P.cant_seed = 0; P.cant_grow = 1;
    P.reg_tau = 0; P.reg_gate = 0; P.reg_f0 = 0;
    P.par_tau = 0; P.par_gate = 0; P.par_form = 0;
    P.par_lo = 0.002; P.par_hi = 0.02; P.par_mature = 400;
    P.wf_on = 0; P.wf_floor = 0.01; P.wf_far = 99;
    P.conf_r = 0; P.conf_gap = 0.3; P.conf_th = 1.6; P.conf_pinw = 3.0;
    P.spin_m = 0; P.imp_k = 0; P.qp_phase = 0;
    P.seed_mw = 0;
    P.bh_r = 0; P.bh_k = 1.0; P.bh_sep = 0;
    P.bed_k = 0; P.bed_tau = 30;
    P.amp_tau = 0;
    P.amp_drv = 0; P.amp_mmin = 2;
    P.amp_nat = 0;
    P.amp_logate = 0;
    P.amp_door = 0;
    P.qatom_every = 200;
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
    P.tri2_sep = 2.6; P.tri2_k2 = 0; P.tri2_kind = -1;
    P.oct_x = 0.325; P.oct_doff = 0.0;
    P.shear_eps = 0; P.shear_t = 0;
    P.slit_wallx = 16.0; P.slit_th = 1.6; P.slit_sep = 9.0; P.slit_hw = 2.0;
    P.slit_screenx = 28.0; P.slit_srcx = 8.0; P.slit_sy = 10.0;
    P.slit_t0 = 16.0; P.slit_t1 = 60.0; P.slit_mask = 0;
    P.slit_pinw = 3.0;
    P.rings_kind = 1; P.rings_nv = 6; P.rings_branch = 0;
    P.rings_xU = 0.28; P.rings_xD = -1; P.rings_gapoff = 0;
    P.rings_xcomp = 0; P.rings_xcompD = 0; P.rings_sep = 9.0;
    P.blob2_sep = 10.0; P.blob2_kx = 0;
    P.slit_obj = 0; P.obj_amp = 1.2; P.obj_sigma = 1.6;
    P.slit_clicks = 0;
    P.p1_meter = 0;
    P.sect_meter = 0; P.sect_n = 24;
    P.sect_r0 = 7.0; P.sect_r1 = 11.0; P.sect_x = -1; P.sect_y = -1;
    P.sect_t0 = 0; P.sect_t1 = 1e18;
    P.obj_y = -1;
    P.convtag = 0; P.grad_r0 = 0; P.grad_frac = 1.0; P.tag_r = 0;
    strcpy(P.snap_dir, "");
    strcpy(P.snap_file, "");
    P.snap_comp = 1;
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
    else if (!strcmp(k, "k_rad")) P.k_rad = atof(v);
    else if (!strcmp(k, "p_rad")) P.p_rad = atof(v);
    else if (!strcmp(k, "rad_clock")) P.rad_clock = atoi(v);
    else if (!strcmp(k, "k_cant")) P.k_cant = atof(v);
    else if (!strcmp(k, "k_tune")) P.k_tune = atof(v);
    else if (!strcmp(k, "cant_tau")) P.cant_tau = atof(v);
    else if (!strcmp(k, "cant_seed")) P.cant_seed = atoi(v);
    else if (!strcmp(k, "cant_grow")) P.cant_grow = atoi(v);
    else if (!strcmp(k, "reg_tau")) P.reg_tau = atof(v);
    else if (!strcmp(k, "reg_gate")) P.reg_gate = atoi(v);
    else if (!strcmp(k, "reg_f0")) P.reg_f0 = atof(v);
    else if (!strcmp(k, "par_tau")) P.par_tau = atof(v);
    else if (!strcmp(k, "par_gate")) P.par_gate = atoi(v);
    else if (!strcmp(k, "par_form")) P.par_form = atoi(v);
    else if (!strcmp(k, "par_lo")) P.par_lo = atof(v);
    else if (!strcmp(k, "par_hi")) P.par_hi = atof(v);
    else if (!strcmp(k, "par_mature")) P.par_mature = atof(v);
    else if (!strcmp(k, "wf_on")) P.wf_on = atoi(v);
    else if (!strcmp(k, "wf_floor")) P.wf_floor = atof(v);
    else if (!strcmp(k, "wf_far")) P.wf_far = atof(v);
    else if (!strcmp(k, "amp_tau")) P.amp_tau = atof(v);
    else if (!strcmp(k, "amp_drv")) P.amp_drv = atof(v);
    else if (!strcmp(k, "amp_mmin")) P.amp_mmin = atof(v);
    else if (!strcmp(k, "amp_nat")) P.amp_nat = atof(v);
    else if (!strcmp(k, "amp_logate")) P.amp_logate = atof(v);
    else if (!strcmp(k, "amp_door")) P.amp_door = atof(v);
    else if (!strcmp(k, "conf_r")) P.conf_r = atof(v);
    else if (!strcmp(k, "conf_gap")) P.conf_gap = atof(v);
    else if (!strcmp(k, "conf_th")) P.conf_th = atof(v);
    else if (!strcmp(k, "conf_pinw")) P.conf_pinw = atof(v);
    else if (!strcmp(k, "spin_m")) P.spin_m = atoi(v);
    else if (!strcmp(k, "seed_mw")) P.seed_mw = atof(v);
    else if (!strcmp(k, "imp_k")) P.imp_k = atof(v);
    else if (!strcmp(k, "qp_phase")) P.qp_phase = atof(v);
    else if (!strcmp(k, "bh_r")) P.bh_r = atof(v);
    else if (!strcmp(k, "bh_k")) P.bh_k = atof(v);
    else if (!strcmp(k, "bh_sep")) P.bh_sep = atof(v);
    else if (!strcmp(k, "bed_k")) P.bed_k = atof(v);
    else if (!strcmp(k, "bed_tau")) P.bed_tau = atof(v);
    else if (!strcmp(k, "qatom_every")) P.qatom_every = atoi(v);
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
    else if (!strcmp(k, "tri2_kind")) P.tri2_kind = atoi(v);
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
    else if (!strcmp(k, "rings_kind")) P.rings_kind = atoi(v);
    else if (!strcmp(k, "rings_nv")) P.rings_nv = atoi(v);
    else if (!strcmp(k, "rings_branch")) P.rings_branch = atoi(v);
    else if (!strcmp(k, "rings_xU")) P.rings_xU = atof(v);
    else if (!strcmp(k, "rings_xD")) P.rings_xD = atof(v);
    else if (!strcmp(k, "rings_gapoff")) P.rings_gapoff = atof(v);
    else if (!strcmp(k, "rings_xcomp")) P.rings_xcomp = atof(v);
    else if (!strcmp(k, "rings_xcompD")) P.rings_xcompD = atof(v);
    else if (!strcmp(k, "blob2_sep")) P.blob2_sep = atof(v);
    else if (!strcmp(k, "blob2_kx")) P.blob2_kx = atof(v);
    else if (!strcmp(k, "slit_obj")) P.slit_obj = atoi(v);
    else if (!strcmp(k, "obj_amp")) P.obj_amp = atof(v);
    else if (!strcmp(k, "obj_sigma")) P.obj_sigma = atof(v);
    else if (!strcmp(k, "slit_clicks")) P.slit_clicks = atoi(v);
    else if (!strcmp(k, "p1_meter")) P.p1_meter = atoi(v);
    else if (!strcmp(k, "sect_meter")) P.sect_meter = atoi(v);
    else if (!strcmp(k, "sect_n")) P.sect_n = atoi(v);
    else if (!strcmp(k, "sect_r0")) P.sect_r0 = atof(v);
    else if (!strcmp(k, "sect_r1")) P.sect_r1 = atof(v);
    else if (!strcmp(k, "sect_x")) P.sect_x = atof(v);
    else if (!strcmp(k, "sect_y")) P.sect_y = atof(v);
    else if (!strcmp(k, "sect_t0")) P.sect_t0 = atof(v);
    else if (!strcmp(k, "sect_t1")) P.sect_t1 = atof(v);
    else if (!strcmp(k, "obj_y")) P.obj_y = atof(v);
    else if (!strcmp(k, "convtag")) P.convtag = atoi(v);
    else if (!strcmp(k, "grad_r0")) P.grad_r0 = atof(v);
    else if (!strcmp(k, "grad_frac")) P.grad_frac = atof(v);
    else if (!strcmp(k, "tag_r")) P.tag_r = atof(v);
    else if (!strcmp(k, "rings_sep")) P.rings_sep = atof(v);
    else if (!strcmp(k, "snap_dir")) { strncpy(P.snap_dir, v, 255); P.snap_dir[255] = 0; }
    else if (!strcmp(k, "snap_file")) { strncpy(P.snap_file, v, 255); P.snap_file[255] = 0; }
    else if (!strcmp(k, "snap_comp")) P.snap_comp = atoi(v);
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
static double *ampre;   /* L-1 zero-sum renorm: pre-bias per-cell outflow total */
static double *sprq, *sscl;
static double *sbed, *bednet;   /* FLOW: per-slot bed weight + net memory */
static double *bedf_;           /* FLOW: per-cell renorm factor scratch */
static double *fsum_;
static double *dm1_, *dm2_;   /* v93 unitary dense channel: psi_m per cell scratch */
static unsigned char *tag;
static unsigned char *pin;   /* apparatus fixture: pass D skips pinned cells */
static unsigned char *scond; /* condensation-active override (DS tier 1 screen) */
static double *fxb, *fyb, *fzb;      /* geometric force gather buffers  */
static double *rngbuf, *nsnap, *th2s;
/* v91 cantus state + pass-H buffers (allocated always; zero when off).
 * v1.1 (CANTUS.md §3.3): the order parameter is LINK-BORNE — sgg_[s]
 * is the slot's own gauge memory (low-passed two-sided gate quality;
 * holds through lens blinks; dies with the slot). The per-cell
 * amplitude is a pure diagnostic (cant_amp_of = max incident sgg). */
static double *cxl_, *dthH;
static double *sgg_;                 /* per-slot cantus amplitude       */
/* v91 exchange registry (REGISTRY.md §1.1): per-slot ledger of
 * reciprocal deliveries. rfp_ = low-passed delivered rate i->j,
 * rfm_ = j->i, rdel_ = per-step delivered scratch [slot][dir].
 * Born 0 at slot_new, dies with the slot — no background, no energy. */
static double *rfp_, *rfm_, *rdel_;

/* v91 IDENTITY lane state (IDENTITY.md §1.1): episode gids on the
 * cells, parcel labels + bond identity stamps + identity-carried
 * delivery ledgers on the slots. Born at mint/slot_new, dies with the
 * episode/slot — no background, no energy moved, nothing anywhere is
 * indexed BY gid (gid_next is a label spring, not an index). */
static long long *cgid_;               /* [cell] episode gid, 0=none  */
static double *cbirth_;                /* [cell] episode birth time   */
static long long gid_next = 1;
static long long *slgid_;              /* [slot][dir] parcel label    */
static double *sborn_;                 /* [slot] slot birth time      */
static long long *pstampa_, *pstampb_; /* [slot] bond identity stamp  */
static double *pstampt_;               /* [slot] stamp time           */
static double *pdp_, *pdm_, *pdel_;    /* id-carried delivered rates  */
static double par_del_tot = 0, par_del_id = 0;
static long long par_mints = 0, par_retires = 0;
static double *par_tmp;                /* diag scratch (quartiles)    */
/* v91 AMPLITUDE Phase M state (AMPLITUDE.md §2): shadow in-flight
 * amplitude per slot-dir, windowed delivered vector+magnitude per
 * slot. Born 0 at slot_new, dies with the slot; moves nothing. */
static double *sre_, *sim_;            /* [slot][dir] in-flight shadow */
static double *amre_, *amim_, *amc_;   /* [slot] windowed delivery     */
static double *amdre_, *amdim_, *amdc_;/* [slot] per-step scratch      */
#define PAR_RING 4096
static double par_aged[2][PAR_RING];   /* slot age at death: 0=tag-pair, 1=bath */
static long long par_agedn[2];
static double *regq1_, *regq2_;      /* diag scratch (quantiles)        */

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
static double *shau_;                /* v93 [slot] unitary dense hop angle tau_s */
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
static double rad_total = 0;   /* v91 graded sub-cap radiance */
static double tune_total = 0;  /* v91 cantus within-mode retune |J| ledger */
static double A0eff = 0;
static long qfire_n = 0;
static double sim_t = 0;
static double E0_total = 0;
static double Eh_total = 0;   /* HORIZON: the hole ledger (pitchless) */
static double Eh_comp = 0;    /* Kahan compensation — the ledger must
                               * hold the drift floor while it grows */
static double bh_eat_f = 0, bh_eat_m = 0, bh_eat_s = 0;
static long bh_nin = 0;

static void eh_add(double v)
{
    double y = v - Eh_comp, t = Eh_total + y;
    Eh_comp = (t - Eh_total) - y;
    Eh_total = t;
}
static double cenx, ceny, cenz;
static long ncand_last = 0;
static double fs_t[512], fs_r[512], fs_y[512], fs_z[512];
static int nfsamp = 0;
#define DS_NBIN 96
static double ds_I[DS_NBIN];
static double ds_expo = 0;
static double *ds_cellI = NULL;   /* per-cell time-integrated Ee on the strip */
/* XSEC sector meter: per-sector time-integrated Ee*dt in the annulus,
 * plus time-integrated cell occupancy n*dt (exposure deficits must be
 * separable from foam-population shifts) */
#define SECT_MAX 96
static double sectE[SECT_MAX], sectN[SECT_MAX];
static double sect_cx = 0, sect_cy = 0;   /* resolved centre */
/* XSEC tag-split conversion ledger: the four global counters again,
 * accumulated where the event's cell is tagged (the occulter) */
static double ct_rough = 0, ct_cond = 0, ct_evap = 0, ct_backs = 0;
static double b2_ax[512], b2_bx[512], b2_tx[512];

/* p1_meter accumulators: cumulative first-moment flow, one 3-vector per
 * channel. Flight sits at the link midpoint (the convention cancels
 * over any completed transfer: deposit + arrival telescope to
 * dE * (x_dst - x_src) minimal-image). */
static double p1sp[3], p1fl[3], p1fd[3], p1gm[3];
#define CLICK_MAX 8192
static double click_t[CLICK_MAX], click_y[CLICK_MAX], click_e[CLICK_MAX];
static int nclick = 0;

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
static void qatom_diag(int fd, double w, double e, int ci, double em)
{
    int qe = P.qatom_every > 0 ? P.qatom_every : 1;
    if (A0eff <= 0 || e <= 0) return;
    if ((qfire_n++ % qe) == 0) {
        if (P.par_tau > 0)
            printf("# QATOM t=%.2f dir=%s w=%.9g e=%.12g i=%d Em=%.4f gid=%lld\n",
                   sim_t, fd ? "FD" : "DF", w, e, ci, em, cgid_[ci]);
        else
            printf("# QATOM t=%.2f dir=%s w=%.9g e=%.12g i=%d Em=%.4f\n",
                   sim_t, fd ? "FD" : "DF", w, e, ci, em);
    }
}

static void field_inject(int i, double dE)   /* cellfab.c:2711 verbatim */
{
    if (dE <= 0) return;
    double e = fa1[i]*fa1[i] + fa2[i]*fa2[i];
    if (P.amp_door > 0) {
        /* v93 symmetric reverse door (§II.7): the evaporated matter click
         * (sqrt(amp_door*dE) at phase th2) composes COHERENTLY with the
         * existing field amplitude; arg(fa_new) = the coherent-sum direction,
         * |fa_new| = sqrt(e+dE) (conserved). Reduces to the real-scale path
         * as amp_door->0. Fixes the one-way door: matter->field now carries
         * th2 (so winding imprinted at condensation is not erased on evap). */
        double argf = (e > 1e-20) ? atan2(fa2[i], fa1[i]) : th2[i];
        double ro = sqrt(e), rc = sqrt(P.amp_door * dE);
        double s1 = ro * lut_cos(argf) + rc * lut_cos(th2[i]);
        double s2 = ro * lut_sin(argf) + rc * lut_sin(th2[i]);
        double rm = sqrt(e + dE);
        double np = atan2(s2, s1);
        fa1[i] = rm * lut_cos(np);
        fa2[i] = rm * lut_sin(np);
    } else if (e > 1e-20) {
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
    ampre = malloc(nc * sizeof(double));
    sprq = malloc(nc * sizeof(double)); sscl = malloc(nc * sizeof(double));
    bedf_ = malloc(nc * sizeof(double));
    fsum_ = malloc(nc * sizeof(double));
    dm1_ = malloc(nc * sizeof(double)); dm2_ = malloc(nc * sizeof(double));
    tag = calloc(nc, 1);
    pin = calloc(nc, 1);
    scond = calloc(nc, 1);
    fxb = malloc(nc * sizeof(double)); fyb = malloc(nc * sizeof(double)); fzb = malloc(nc * sizeof(double));
    rngbuf = malloc(6 * (size_t)nc * sizeof(double));
    nsnap = malloc(6 * (size_t)nc * sizeof(double));
    th2s = malloc(nc * sizeof(double));
    cxl_ = calloc(nc, sizeof(double));
    dthH = calloc(nc, sizeof(double));
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
    shau_ = calloc(NLMAX, sizeof(double));
    sflux = calloc(NLMAX, sizeof(double));
    sldd = calloc(NLMAX, sizeof(double));
    swl = calloc(NLMAX, sizeof(double));
    sbed = malloc(NLMAX * sizeof(double));
    for (int s0 = 0; s0 < NLMAX; s0++) sbed[s0] = 1.0;
    bednet = calloc(NLMAX, sizeof(double));
    sgg_ = calloc(NLMAX, sizeof(double));
    rfp_ = calloc(NLMAX, sizeof(double));
    rfm_ = calloc(NLMAX, sizeof(double));
    rdel_ = calloc((size_t)2 * NLMAX, sizeof(double));
    cgid_ = calloc(nc, sizeof(long long));
    cbirth_ = calloc(nc, sizeof(double));
    slgid_ = calloc((size_t)2 * NLMAX, sizeof(long long));
    sborn_ = calloc(NLMAX, sizeof(double));
    pstampa_ = calloc(NLMAX, sizeof(long long));
    pstampb_ = calloc(NLMAX, sizeof(long long));
    pstampt_ = calloc(NLMAX, sizeof(double));
    pdp_ = calloc(NLMAX, sizeof(double));
    pdm_ = calloc(NLMAX, sizeof(double));
    pdel_ = calloc((size_t)2 * NLMAX, sizeof(double));
    par_tmp = malloc(nc > NLMAX ? nc * sizeof(double) : NLMAX * sizeof(double));
    sre_ = calloc((size_t)2 * NLMAX, sizeof(double));
    sim_ = calloc((size_t)2 * NLMAX, sizeof(double));
    amre_ = calloc(NLMAX, sizeof(double));
    amim_ = calloc(NLMAX, sizeof(double));
    amc_ = calloc(NLMAX, sizeof(double));
    amdre_ = calloc(NLMAX, sizeof(double));
    amdim_ = calloc(NLMAX, sizeof(double));
    amdc_ = calloc(NLMAX, sizeof(double));
    regq1_ = malloc(NLMAX * sizeof(double));
    regq2_ = malloc(NLMAX * sizeof(double));
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
    sbed[s] = 1.0; bednet[s] = 0;    /* FLOW: a reborn link has no bed */
    sgg_[s] = 0;                     /* cantus: a reborn bond starts mute */
    rfp_[s] = 0; rfm_[s] = 0;        /* registry: a reborn pair has no past */
    rdel_[2*s] = 0; rdel_[2*s+1] = 0;
    sborn_[s] = sim_t;               /* identity lane: fresh slot, no past */
    slgid_[2*s] = 0; slgid_[2*s+1] = 0;
    pstampa_[s] = 0; pstampb_[s] = 0; pstampt_[s] = 0;
    pdp_[s] = 0; pdm_[s] = 0; pdel_[2*s] = 0; pdel_[2*s+1] = 0;
    sre_[2*s] = 0; sre_[2*s+1] = 0; sim_[2*s] = 0; sim_[2*s+1] = 0;
    amre_[s] = 0; amim_[s] = 0; amc_[s] = 0;
    amdre_[s] = 0; amdim_[s] = 0; amdc_[s] = 0;
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
            if (P.par_tau > 0) {
                /* identity lane meter: slot age at death, by class
                 * (M-I2 feeds the par_mature selection I-G2b) */
                int pc = tag[sli[s]] && tag[slj[s]] ? 0
                       : (!tag[sli[s]] && !tag[slj[s]] ? 1 : -1);
                if (pc >= 0) {
                    par_aged[pc][par_agedn[pc] % PAR_RING] = sim_t - sborn_[s];
                    par_agedn[pc]++;
                }
            }
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
            if (P.bed_k > 0) swl[s] *= sbed[s];   /* FLOW: the bed */
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
        if (P.p1_meter) {
            /* P1 site 1: space transport src->dst = sgn(f) * d * u_hat */
            for (int s = 0; s < NSLOT; s++) {
                double f = swl[s];
                if (f == 0) continue;
                int src = f > 0 ? sli[s] : slj[s];
                double mag = fabs(f) * sscl[src];
                if (mag <= 0) continue;
                double m = (f > 0 ? 1.0 : -1.0) * mag * sd[s];
                p1sp[0] += m * sux[s]; p1sp[1] += m * suy[s]; p1sp[2] += m * suz[s];
            }
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
        /* FLOW (FLOW.md): the bed-digging channel law. Net memory
         * from the ACTUAL moved flow; growth on |net|; ZERO-SUM
         * per-cell renorm (the anti-ignition structure); clamp. */
        if (P.bed_k > 0) {
            double kb = P.bed_tau > dt ? dt / P.bed_tau : 1.0;
            for (int s = 0; s < NSLOT; s++) {
                if (sst[s] == S_FREE) continue;
                double f = swl[s], flow = 0;
                if (f != 0) {
                    int src = f > 0 ? sli[s] : slj[s];
                    flow = (f > 0 ? 1.0 : -1.0) * fabs(f) * sscl[src];
                }
                bednet[s] += kb * (flow / dt - bednet[s]);
                sbed[s] *= 1.0 + P.bed_k * dt * fabs(bednet[s]);
            }
            for (int i = 0; i < NC; i++) {
                double sum = 0; int n = 0;
                for (int q = cls_[i]; q < cls_[i + 1]; q++) {
                    int s = clidx[q];
                    if (sst[s] == S_FREE || sA[s] <= 0) continue;
                    sum += sbed[s]; n++;
                }
                bedf_[i] = (n > 0 && sum > 0) ? (double)n / sum : 1.0;
            }
            for (int s = 0; s < NSLOT; s++) {
                if (sst[s] == S_FREE || sA[s] <= 0) continue;
                sbed[s] *= sqrt(bedf_[sli[s]] * bedf_[slj[s]]);
                if (sbed[s] < 0.2) sbed[s] = 0.2;
                if (sbed[s] > 5.0) sbed[s] = 5.0;
            }
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
    if (P.amp_drv > 0 && P.amp_tau > 0)
        for (int i = 0; i < NC; i++) ampre[i] = 0;   /* L-1 zero-sum bookkeeping */
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
        shau_[s] = 0;
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

        /* v92 AMPLITUDE Phase L lane L-1 (L0_DESIGN.md §3, rev phase-current
         * per L1_FINDINGS.md finding 3): promote the Phase-M shadow from
         * meter to driver. The driver is the shadow's IN-PHASE current —
         * project each direction's shadow amplitude onto the RECEIVER's
         * transport frame (the closure phase bq*thi - bq*wi*d/C - bp*thj),
         * so the bias responds to a phase GRADIENT (the e3b tilt), which
         * the magnitude form was blind to. Signed, saturating, slot-borne.
         * Gate amp_mmin: 1 = include unison (e3b/p1 live here); 2 = chords
         * only (the fifth/nv=6 track). Byte-inert at amp_drv=0. */
        if (P.amp_drv > 0 && P.amp_tau > 0) {
            int ms = slp[s] > slq[s] ? slp[s] : slq[s];
            if (ms >= (int)P.amp_mmin) {
                /* record pre-bias outflow for the zero-sum renorm */
                ampre[i] += w_ij;
                ampre[j] += w_ji;
                /* J0 = Re(A_{i->j} e^{-i p th_j}) = |A0| cos(q th_i - p th_j)
                 * J1 = Re(A_{j->i} e^{-i q th_i});  net forward = J0 - J1   */
                double phr0 = slp[s] * th2[j];
                double J0 = sre_[2*s]*lut_cos(phr0)   + sim_[2*s]*lut_sin(phr0);
                double phr1 = slq[s] * th2[i];
                double J1 = sre_[2*s+1]*lut_cos(phr1) + sim_[2*s+1]*lut_sin(phr1);
                double Jnet = J0 - J1;
                double Jabs = fabs(J0) + fabs(J1) + 1e-15;
                double bias = P.amp_drv * Jnet / Jabs;   /* signed, |.|<=amp_drv */
                w_ij *= (1.0 + bias);
                w_ji *= (1.0 - bias);
            }
        }

        if (P.amp_nat > 0) {
            /* v93 UNITARY DENSE CHANNEL: fold the want ingredients into a
             * pairwise hop angle tau_s (the existing want-computation maps
             * onto the angle rather than a magnitude to debit). Mobility is
             * carried by the amplitude magnitude itself (as in pass F: a
             * cell with Em~0 has psi_m~0 and transports nothing). The door
             * (pass 6) enforces cap, so NO head factor here -- face 1 used
             * sqrt(head_i*head_j) but that FREEZES the dense core (cap-
             * saturated matter has head~0), leaving only the random-phase
             * surface to transport (e3b cos seed-variant, incoherent). The
             * closure gate (cos^p) survives as the angle envelope. tau_s=0
             * when closure vanishes. Byte-inert at amp_nat=0. */
            double gsym = sqrt(g_ij * g_ji);
            /* amp_logate (LINEARIZE probe, reviewer item 1): drop the phase-
             * dependent gate from tau_s -> a phase-independent linear hop
             * (tau = amp_nat*base). Tests whether the gate's phase->tau
             * feedback (the parametric drive) is the artifact source. The
             * gated form (cos^p closure envelope) is the default. */
            double taub = P.amp_nat * base * gsym;
            if (P.amp_logate > 0) taub = P.amp_nat * base;   /* LINEAR: no gate */
            if (taub > 0.5) taub = 0.5;   /* integrator guard: deep-overlap cap */
            shau_[s] = taub;
        } else {
            if (w_ij > 0) swant[2*s]   = w_ij;
            if (w_ji > 0) swant[2*s+1] = w_ji;
        }
    }

    /* L-1 zero-sum renorm (FLOW architecture, L1_FINDINGS/grok fix): hold
     * each cell's total outflow at its pre-bias level so the bias only
     * REDISTRIBUTES direction (anti-ignition + coherence-runaway bound),
     * never amplifying throughput. Energy conserved by construction. */
    if (P.amp_drv > 0 && P.amp_tau > 0) {
        for (int i = 0; i < NC; i++) {
            if (ampre[i] <= 0) continue;
            double post = 0;
            for (int q = cls_[i]; q < cls_[i+1]; q++) {
                int s = clidx[q];
                int dir = (sli[s] == i) ? 0 : 1;
                post += swant[2*s + dir];
            }
            if (post > 1e-15) {
                double fac = ampre[i] / post;
                for (int q = cls_[i]; q < cls_[i+1]; q++) {
                    int s = clidx[q];
                    int dir = (sli[s] == i) ? 0 : 1;
                    swant[2*s + dir] *= fac;
                }
            }
        }
    }

    /* pass U: v93 UNITARY DENSE CHANNEL — within-mode dense transport as a
     * product of UNITARY PAIRWISE PLANE ROTATIONS (pass F's cousin,
     * v93/README.md §II.3). Engages only when amp_nat>0; in that regime the
     * additive magnitude want above is bypassed (swant stays 0, so passes
     * 3-5 are no-ops). The dense amplitude psi_m = sqrt(Em) e^{i th2} hops
     * between link endpoints by a Givens rotation of angle shau_[s]:
     *   dm1[i]' = cc dm1[i] + ss dm2[j];  dm2[i]' = cc dm2[i] - ss dm1[j]
     *   dm1[j]' = cc dm1[j] + ss dm2[i];  dm2[j]' = cc dm2[j] - ss dm1[i]
     * Each hop conserves the two-cell norm EXACTLY (conservation is a
     * theorem of the update, not a patched ledger) -- the same mechanism
     * that makes the field unitary to roundoff. No sum-then-square ever
     * occurs (IV.4). The cross term 2 Im(psi_i* psi_j) that the additive
     * Em-ledger rejected IS the link current J_s = the dense momentum
     * (II.4/II.6). The door (pass 6) is never unitarized (IV.3). Hops are
     * applied in canonical link order (from the i side), like pass F.
     * Byte-inert at amp_nat=0 (this whole pass is skipped). */
    if (P.amp_nat > 0) {
        for (int i = 0; i < NC; i++) {
            double e = Em[i];
            if (e > 0) {
                double r = sqrt(e);
                dm1_[i] = r * lut_cos(th2[i]);
                dm2_[i] = r * lut_sin(th2[i]);
            } else {
                dm1_[i] = dm2_[i] = 0.0;
            }
        }
        /* Local clock precession IN-pass (face 3, mirroring pass F exactly):
         * rotate psi_m by +w2e*dt (== th2 += w2e dt) BEFORE the hops. The
         * out-of-pass advance (pass 6) is skipped when amp_nat>0 -- the
         * precess-then-hop order removes the first-order Trotter split that
         * face-2 showed as centroid wobble. Still byte-inert at amp_nat=0
         * (this whole block, and the pass-6 skip, are amp_nat-gated). */
        for (int i = 0; i < NC; i++) {
            double ang = w2e[i] * dt;
            double cc = lut_cos(ang), ss = lut_sin(ang);
            double a = dm1_[i], b = dm2_[i];
            dm1_[i] = cc * a - ss * b;
            dm2_[i] = ss * a + cc * b;
        }
        for (int i = 0; i < NC; i++) {
            for (int q = cls_[i]; q < cls_[i + 1]; q++) {
                int s = clidx[q];
                if (sli[s] != i) continue;       /* canonical: apply from the i side */
                if (sst[s] == S_FREE || sA[s] <= 0) continue;
                double tau = shau_[s];
                if (tau <= 0.0) continue;
                int j = slj[s];
                double cc, ss;
                sincos(tau, &ss, &cc);
                double m1i = dm1_[i], m2i = dm2_[i], m1j = dm1_[j], m2j = dm2_[j];
                dm1_[i] = cc * m1i + ss * m2j;
                dm2_[i] = cc * m2i - ss * m1j;
                dm1_[j] = cc * m1j + ss * m2i;
                dm2_[j] = cc * m2j - ss * m1i;
                if (P.p1_meter) {
                    /* P1 site (dense): energy the rotation moved into j,
                     * pairwise-conserved by construction; displacement i -> j */
                    double mj = (dm1_[j]*dm1_[j] + dm2_[j]*dm2_[j])
                              - (m1j*m1j + m2j*m2j);
                    double m = mj * sd[s];
                    p1fd[0] += m * sux[s]; p1fd[1] += m * suy[s]; p1fd[2] += m * suz[s];
                }
            }
        }
        for (int i = 0; i < NC; i++) {
            double a = dm1_[i], b = dm2_[i];
            double enew = a * a + b * b;
            Em[i] = enew;
            if (enew > 1e-12) {              /* preserve old th2 for empty cells
                * (atan2(0,0)=0 would pin every briefly-empty cell's clock to 0
                * every step = a phase-slip machine; reviewer fix) */
                double ph = atan2(b, a);
                if (ph < 0) ph += TWO_PI;    /* keep [0, 2pi) like pass 6 */
                th2[i] = ph;
            }
        }
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
            /* identity lane (IDENTITY.md §1.3.2): the departing parcel
             * carries the depositor's episode gid (I-D1 last-depositor
             * label; M-I4 measures the approximation honestly). */
            if (P.par_tau > 0) slgid_[2*s + dir] = cgid_[src];
            /* AMPLITUDE Phase M: the shadow deposit sqrt(f) e^{i m th} */
            if (P.amp_tau > 0) {
                double m_ = dir == 0 ? (double)slq[s] : (double)slp[s];
                double ph = m_ * th2[src];
                sre_[2*s + dir] += sqrt(f) * lut_cos(ph);
                sim_[2*s + dir] += sqrt(f) * lut_sin(ph);
            }
            sflux[s] += f;
            sfluxd[2*s + dir] += f;   /* directed ledger for circulation */
            if (P.p1_meter && f > 0) {
                /* P1 site 2: dense deposit src -> link midpoint */
                double m = (dir == 0 ? 0.5 : -0.5) * f * sd[s];
                p1fl[0] += m * sux[s]; p1fl[1] += m * suy[s]; p1fl[2] += m * suz[s];
            }
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
            /* AMPLITUDE Phase M: transit rotation -m w dt (the psi_e
             * flight phase, applied per step while in flight) */
            if (P.amp_tau > 0 && (sre_[sslot] != 0 || sim_[sslot] != 0)) {
                double m_ = dir == 0 ? (double)slq[s] : (double)slp[s];
                double dph = -m_ * w2e[send] * dt;
                double c_ = lut_cos(dph), s_ = lut_sin(dph);
                double re0 = sre_[sslot];
                sre_[sslot] = re0 * c_ - sim_[sslot] * s_;
                sim_[sslot] = re0 * s_ + sim_[sslot] * c_;
            }
            if (slph[sslot] < 1.0) continue;
            double freec = P.cap - (Em[recv] + Ee[recv]);
            double take = slem[sslot];
            if (take > freec) take = freec > 0 ? freec : 0;
            if (take > 0) {
                double mobprev = Em[recv];
                /* v91 registry stamp (REGISTRY.md §1.3): a delivered
                 * parcel between two continuous identities. Deliveries
                 * only — the rule-alpha flush below is a RETURN, not an
                 * exchange, and is not stamped. */
                if (P.reg_tau > 0) rdel_[sslot] += take;
                /* identity lane (IDENTITY.md §1.3.3): identity-carried
                 * delivery — counts only when the parcel's label is
                 * the source endpoint's LIVING episode (stale flight
                 * from a dead episode is foreign). Deliveries only;
                 * the rule-α flush below is a return, not stamped. */
                if (P.par_tau > 0) {
                    par_del_tot += take;
                    if (slgid_[sslot] != 0 && slgid_[sslot] == cgid_[send]) {
                        pdel_[sslot] += take;
                        par_del_id += take;
                    }
                }
                /* AMPLITUDE Phase M: deliver the amplitude fraction,
                 * composed into the dst chart frame (rotate by
                 * -n th_dst); magnitude and vector accumulate
                 * separately — their ratio is rho_coh. */
                if (P.amp_tau > 0) {
                    double Efl = slem[sslot];
                    double A2 = sre_[sslot]*sre_[sslot] + sim_[sslot]*sim_[sslot];
                    if (Efl > 1e-30 && A2 > 0) {
                        double frac = take / Efl;
                        if (frac > 1) frac = 1;
                        double sf = sqrt(frac);
                        double dre = sre_[sslot] * sf, dim = sim_[sslot] * sf;
                        double keep = sqrt(1.0 - frac);
                        sre_[sslot] *= keep; sim_[sslot] *= keep;
                        double n_ = dir == 0 ? (double)slp[s] : (double)slq[s];
                        double ph = -n_ * th2[recv];
                        double c_ = lut_cos(ph), s2_ = lut_sin(ph);
                        amdre_[s] += dre * c_ - dim * s2_;
                        amdim_[s] += dre * s2_ + dim * c_;
                        amdc_[s] += sqrt(dre*dre + dim*dim);
                    }
                }
                slem[sslot] -= take;
                if (P.p1_meter) {
                    /* P1 site 3: flight arrival, link midpoint -> recv */
                    double m = (dir == 0 ? 0.5 : -0.5) * take * sd[s];
                    p1fl[0] += m * sux[s]; p1fl[1] += m * suy[s]; p1fl[2] += m * suz[s];
                }
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
                if (tag[recv]) { ct_rough += rough; ct_backs += back_s; }
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
                if (P.p1_meter) {
                    /* P1 site 4: flush, link midpoint -> send (return) */
                    double m = (dir == 0 ? -0.5 : 0.5) * slem[sslot] * sd[s];
                    p1fl[0] += m * sux[s]; p1fl[1] += m * suy[s]; p1fl[2] += m * suz[s];
                }
                Em[send] += slem[sslot];
                beta_energy += slem[sslot];
                beta_returns++;
                slem[sslot] = 0; slph[sslot] = 0;
                sre_[sslot] = 0; sim_[sslot] = 0;  /* flush = return, not delivery */
            }
            if (slem[sslot] <= 1e-17) { slem[sslot] = 0; slph[sslot] = 0;
                                        sre_[sslot] = 0; sim_[sslot] = 0; }
            else if (take <= 0) slph[sslot] = 0;
            else slph[sslot] -= 1.0;
        }
    }

    /* pass H: v91 CANTUS (coherent-channel candidate B, CANTUS.md).
     * The superimposed harmonic-lock field, v1.1 (CANTUS.md §3.3):
     * sgg_[s] = the LINK-BORNE order parameter (low-passed two-sided
     * gate quality — one gauge per bond; grows only where locked
     * exchange PERSISTS on that bond, holds through lens blinks,
     * dies with the slot; no background, no energy); cxl = holdings
     * memory (each voice's remembered part in the chord — R-D1-
     * aligned: what radiance reads). v1's cell-borne amplitude
     * (max-support) ignited a bath-wide Kuramoto transition —
     * measured, rejected, recorded in CANTUS.md §4.1.
     * (a) the lock corrects the DIFFERENTIAL ladder-closure
     * error on the matter clocks th2 (the common error is pure
     * geometry — the bond walk's job; the phase sum contains no
     * phase). (b) the within-mode retune current moves holdings
     * along links on the memory-DEVIATION difference (fast piles
     * flatten; slow chord structure/flavor untouched), pairwise-
     * conserving, unquantized (within a mode — the standing law),
     * COE-metered. Slots visited once from the lower endpoint in
     * CSR-canonical order (the pass-4/5 convention, A/B-mirrored).
     * k_cant=0 && k_tune=0 => this pass does not execute. */

    /* pass H0: v91 EXCHANGE REGISTRY ledger (REGISTRY.md §1.3 item 2).
     * Low-pass the per-step delivered parcels into directed rates.
     * Runs over ALL non-free slots (including pinched sA=0 ones, so
     * the memory decays through lens blinks instead of freezing);
     * s ascending — A/B-canonical. reg_tau=0 => does not execute. */
    if (P.reg_tau > 0) {
        double kreg = P.reg_tau > dt ? dt / P.reg_tau : 1.0;
        for (int s = 0; s < NSLOT; s++) {
            if (sst[s] == S_FREE) continue;
            rfp_[s] += kreg * (rdel_[2*s]     / dt - rfp_[s]);
            rfm_[s] += kreg * (rdel_[2*s + 1] / dt - rfm_[s]);
            rdel_[2*s] = 0; rdel_[2*s + 1] = 0;
        }
    }

    /* pass H0b: v91 IDENTITY lane (IDENTITY.md §1.3 items 1+4).
     * First the cells: episode gids minted/retired by the hysteresis
     * pair on x = (Em+flload)/cap (flload = the start-of-step
     * snapshot — the same currency the diag reads). Then the slots:
     * low-pass the identity-carried deliveries, arm/clear the bond
     * identity stamps. Serial (local-clock kernel), s ascending —
     * A/B-canonical. par_tau=0 => does not execute. */
    if (P.par_tau > 0) {
        for (int i = 0; i < NC; i++) {
            double xep = (Em[i] + flload[i]) / P.cap;
            if (cgid_[i] == 0) {
                if (xep >= P.par_hi) {
                    cgid_[i] = gid_next++;
                    cbirth_[i] = sim_t;
                    par_mints++;
                }
            } else if (xep < P.par_lo) {
                cgid_[i] = 0;
                par_retires++;
            }
        }
        double kpar = P.par_tau > dt ? dt / P.par_tau : 1.0;
        for (int s = 0; s < NSLOT; s++) {
            if (sst[s] == S_FREE) continue;
            pdp_[s] += kpar * (pdel_[2*s]     / dt - pdp_[s]);
            pdm_[s] += kpar * (pdel_[2*s + 1] / dt - pdm_[s]);
            pdel_[2*s] = 0; pdel_[2*s + 1] = 0;
            long long ga = cgid_[sli[s]], gb = cgid_[slj[s]];
            if (pstampa_[s] == 0) {
                /* arm: first mutual identity-carried exchange between
                 * two living episodes stamps WHO is bonded */
                if (ga != 0 && gb != 0 && pdp_[s] > 0 && pdm_[s] > 0) {
                    pstampa_[s] = ga; pstampb_[s] = gb; pstampt_[s] = sim_t;
                }
            } else if (ga != pstampa_[s] || gb != pstampb_[s]) {
                /* either identity died or changed: the bond's history
                 * ends; a new pair must re-stamp and re-mature */
                pstampa_[s] = 0; pstampb_[s] = 0; pstampt_[s] = 0;
            }
        }
    }

    /* pass H0c: v91 AMPLITUDE Phase M (AMPLITUDE.md §1) — low-pass the
     * per-step delivered shadow into windowed vector + magnitude
     * rates. amp_tau=0 => does not execute. */
    if (P.amp_tau > 0) {
        double kam = P.amp_tau > dt ? dt / P.amp_tau : 1.0;
        for (int s = 0; s < NSLOT; s++) {
            if (sst[s] == S_FREE) continue;
            amre_[s] += kam * (amdre_[s] / dt - amre_[s]);
            amim_[s] += kam * (amdim_[s] / dt - amim_[s]);
            amc_[s]  += kam * (amdc_[s]  / dt - amc_[s]);
            amdre_[s] = 0; amdim_[s] = 0; amdc_[s] = 0;
        }
    }

    if (P.k_cant > 0 || P.k_tune > 0) {
        double ktau = P.cant_tau > dt ? dt / P.cant_tau : 1.0;
        for (int i = 0; i < NC; i++) {
            th2s[i] = th2[i];
            dthH[i] = 0;
        }
        for (int i = 0; i < NC; i++) {
            for (int q2 = cls_[i]; q2 < cls_[i + 1]; q2++) {
                int s = clidx[q2];
                if (sli[s] != i) continue;       /* visit once, from i */
                if (sst[s] == S_FREE || sA[s] <= 0) continue;
                int j = slj[s];
                if (Em[i] <= 1e-15 || Em[j] <= 1e-15) continue;
                double pp = slp[s], qq = slq[s];
                double d = sd[s];
                double ps_f = wrap_pi(qq*th2s[i] - qq*w2e[i]*d/P.C - pp*th2s[j]);
                double ps_b = wrap_pi(pp*th2s[j] - pp*w2e[j]*d/P.C - qq*th2s[i]);
                double gg = gate_of(ps_f) * gate_of(ps_b);
                /* v1.1: the LINK's own gauge memory (holds through
                 * lens blinks — non-eligible steps skip this update).
                 * reg_gate=1 (REGISTRY.md R-G3): the growth TARGET is
                 * gated by the registry match r = 2*min/(sum) of the
                 * directed delivery rates (form F-A) — coherence may
                 * only grow on identity-continuous reciprocal
                 * exchange; the bath stays dark by construction. */
                if (P.cant_grow || sgg_[s] > 0) {
                    double tgt = gg;
                    if (P.reg_gate) {
                        double gross = rfp_[s] + rfm_[s];
                        double mn = rfp_[s] < rfm_[s] ? rfp_[s] : rfm_[s];
                        double mult = 0.0;
                        if (P.reg_gate == 1) {
                            /* F-B: balance x flow saturation (f0=0 => F-A) */
                            if (gross > 0)
                                mult = (2.0 * mn / gross) * (gross / (gross + P.reg_f0));
                        } else {
                            /* F-D: reciprocal-flow saturation, s = 2*min */
                            double s2 = 2.0 * mn;
                            if (s2 > 0) mult = s2 / (s2 + P.reg_f0);
                        }
                        tgt *= mult;
                    }
                    if (P.par_gate) {
                        /* IDENTITY.md §1.3.5–6: growth only on
                         * maturity-clocked identity continuity. The
                         * gauge has no force before arming, so
                         * nothing it does can extend a bath slot's
                         * life to par_mature — stamp AGE is the gate
                         * variable the lock cannot manufacture. */
                        double rid = 0.0;
                        if (pstampa_[s] != 0
                            && cgid_[i] == pstampa_[s]
                            && cgid_[j] == pstampb_[s]
                            && pdp_[s] > 0 && pdm_[s] > 0) {
                            double page = sim_t - pstampt_[s];
                            double mat = (P.par_mature > 0 && page < P.par_mature)
                                         ? page / P.par_mature : 1.0;
                            if (P.par_form == 1) {
                                double gr2 = pdp_[s] + pdm_[s];
                                double mn2 = pdp_[s] < pdm_[s] ? pdp_[s] : pdm_[s];
                                rid = gr2 > 0 ? mat * 2.0 * mn2 / gr2 : 0.0;
                            } else rid = mat;
                        }
                        tgt *= rid;
                    }
                    sgg_[s] += ktau * (tgt - sgg_[s]);
                }
                double amp = sgg_[s];
                if (amp <= 0) continue;
                if (P.k_cant > 0) {
                    double e = 0.5 * wrap_pi(ps_f - ps_b);
                    double wl = P.k_cant * dt * amp;
                    double n2q = pp * pp + qq * qq;
                    dthH[i] -= wl * (qq / n2q) * e;
                    dthH[j] += wl * (pp / n2q) * e;
                }
                if (P.k_tune > 0) {
                    double ui = Em[i] / P.cap - cxl_[i];
                    double uj = Em[j] / P.cap - cxl_[j];
                    double J = P.k_tune * dt * amp * P.cap * (ui - uj);
                    int src = J > 0 ? i : j, dst = J > 0 ? j : i;
                    double mag = fabs(J);
                    double av = 0.98 * Em[src];
                    if (mag > av) mag = av;
                    double freec = P.cap - (Em[dst] + Ee[dst]);
                    if (freec < 0) freec = 0;
                    if (mag > freec) mag = freec;
                    if (mag > 0) {
                        Em[src] -= mag;
                        Em[dst] += mag;
                        tune_total += mag;
                        if (P.p1_meter) {
                            /* P1 site H: within-mode current src->dst
                             * over the full link (deposit+arrival
                             * telescoped) */
                            double m = (src == i ? 1.0 : -1.0) * mag * d;
                            p1fl[0] += m * sux[s];
                            p1fl[1] += m * suy[s];
                            p1fl[2] += m * suz[s];
                        }
                    }
                }
            }
        }
        for (int i = 0; i < NC; i++) {
            if (dthH[i] != 0) th2[i] += dthH[i];
            cxl_[i] += ktau * (Em[i] / P.cap - cxl_[i]);
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
        if (P.p1_meter) {
            /* P1 site 6: geometry — cells carry their energy when they
             * move; link flight (midpoint) moves with the mean of its
             * endpoints */
            for (int i = 0; i < NC; i++) {
                double ei = Es[i] + Em[i] + Ee[i];
                p1gm[0] += ei * fxb[i]; p1gm[1] += ei * fyb[i]; p1gm[2] += ei * fzb[i];
            }
            for (int s = 0; s < NSLOT; s++) {
                if (sst[s] == S_FREE) continue;
                double le = slem[2*s] + slem[2*s+1];
                if (le <= 0) continue;
                int i2 = sli[s], j2 = slj[s];
                p1gm[0] += le * 0.5 * (fxb[i2] + fxb[j2]);
                p1gm[1] += le * 0.5 * (fyb[i2] + fyb[j2]);
                p1gm[2] += le * 0.5 * (fzb[i2] + fzb[j2]);
            }
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
            if (P.p1_meter) {
                /* P1 site 5: field hop — energy the rotation moved into j
                 * (pairwise-conserved), displacement i -> j */
                double mj = (fa1[j]*fa1[j] + fa2[j]*fa2[j]) - (a1j*a1j + a2j*a2j);
                double m = mj * sd[s];
                p1fd[0] += m * sux[s]; p1fd[1] += m * suy[s]; p1fd[2] += m * suz[s];
            }
        }
    }
    for (int i = 0; i < NC; i++)
        Ee[i] = fa1[i]*fa1[i] + fa2[i]*fa2[i];

    /* pass 6: dense clock + beat-gated conversion (cellfab.c:3448),
     * serial by design */
    for (int i = 0; i < NC; i++) {
        if (roughq[i] > 0) {
            qatom_diag(0, atoms_w(w2e[i], w1e[i]), roughq[i], i, Em[i]);
            roughq[i] = 0;
        }
    }
    bh_nin = 0;
    for (int i = 0; i < NC; i++) {
        if (P.bh_r > 0) {
            double hdx = wr(px_[i] - 0.5 * P.L), hdy = wr(py_[i] - 0.5 * P.L);
            double hdz = wr(pz_[i] - 0.5 * P.L);
            if (P.bh_sep > 0) {
                /* FLOW: two centres at cx +- sep/2; use the nearer */
                double dxa = wr(px_[i] - (0.5 * P.L - 0.5 * P.bh_sep));
                double dxb = wr(px_[i] - (0.5 * P.L + 0.5 * P.bh_sep));
                hdx = fabs(dxa) < fabs(dxb) ? dxa : dxb;
            }
            if (hdx*hdx + hdy*hdy + hdz*hdz < P.bh_r * P.bh_r) {
                bh_nin++;
                if (Ee[i] > 0) {
                    eh_add(Ee[i]); bh_eat_f += Ee[i];
                    Ee[i] = 0; fa1[i] = 0; fa2[i] = 0;
                }
                if (Em[i] > 0) {
                    eh_add(Em[i]); bh_eat_m += Em[i];
                    Em[i] = 0;
                }
                double avail = Es[i] - P.es_floor;
                if (avail > 0) {
                    double d = P.bh_k * avail * dt;
                    if (d > avail) d = avail;
                    Es[i] -= d; eh_add(d); bh_eat_s += d;
                }
            continue;   /* pitchless: no clock, no beat, no door */
            }
        }
        /* v93 face 3: when amp_nat>0 the dense clock precesses IN pass U
         * (precess-then-hop, mirroring pass F); skip the out-of-pass advance
         * here. Byte-inert (amp_nat==0 takes the original branch). */
        if (P.amp_nat == 0)
            th2[i] = fmod(th2[i] + w2e[i] * dt, TWO_PI);
        cbeta[i] += (w1e[i] - w2e[i]) * dt;
        int beat_fire = 0;
        if (cbeta[i] >= TWO_PI) { cbeta[i] -= TWO_PI; beat_fire = 1; }
        else if (cbeta[i] <= -TWO_PI) { cbeta[i] += TWO_PI; beat_fire = 1; }
        if (beat_fire) {
            double econd_i = scond[i] ? 0.0
                : (P.wf_on ? (Em[i] >= P.wf_floor ? 0.0 : P.wf_far)
                           : P.e_cond);
            if (Ee[i] > econd_i) {
                double d1 = P.f_conv * (Ee[i] - econd_i);
                double eF = A0eff * w1e[i] / TWO_PI;
                double eD = A0eff * w2e[i] / TWO_PI;
                d1 = atoms_fire(d1, eF, eD, &qcnvD[i]);
                d1 = atoms_clamp(d1, 0.98 * Ee[i], eF, eD, &qcnvD[i]);
                if (d1 > 0) {
                    cond_total += d1;
                    if (tag[i]) ct_cond += d1;
                    if (scond[i] && nclick < CLICK_MAX) {
                        click_t[nclick] = sim_t;
                        click_y[nclick] = py_[i];
                        click_e[nclick] = d1;
                        nclick++;
                    }
                    qatom_diag(1, atoms_w(w1e[i], w2e[i]), d1, i, Em[i]);
                    double dsp = P.s_pull * d1;
                    double avail = Es[i] - P.es_floor;
                    if (avail < 0) avail = 0;
                    if (dsp > avail) dsp = avail;
                    double fac = Ee[i] > 0 ? sqrt((Ee[i] - d1) / Ee[i]) : 0;
                    fa1[i] *= fac; fa2[i] *= fac;
                    Ee[i] -= d1;
                    Es[i] -= dsp;
                    Em[i] += d1 + dsp;
                    if (P.amp_door > 0) {
                        /* v93 arg(psi) door (§II.7): the field click composes
                         * COHERENTLY with the existing matter amplitude;
                         * arg(psi_m_new) = the coherent-sum direction (carries
                         * the field phase + the interference cross-term),
                         * |psi_m| fixed by conserved Em. amp_door is a CONTINUOUS
                         * mix weight: the click amplitude is sqrt(amp_door*d1),
                         * so the door reduces to th2_old as amp_door->0 and is
                         * the full merge at amp_door=1. Fired atom carries
                         * arg(psi_m), not m*th2. */
                        double aph = atan2(fa2[i], fa1[i]);
                        double Em_old = Em[i] - (d1 + dsp);
                        if (Em_old < 0) Em_old = 0;
                        double ro = sqrt(Em_old), rc = sqrt(P.amp_door * d1);
                        double s1 = ro * lut_cos(th2[i]) + rc * lut_cos(aph);
                        double s2 = ro * lut_sin(th2[i]) + rc * lut_sin(aph);
                        th2[i] = fmod(atan2(s2, s1) + 8.0 * TWO_PI, TWO_PI);
                    } else if (P.qp_phase > 0) {
                        /* QUENCH-3: phase crosses the door (mix ≤ 1 by
                         * construction: d1 ≤ Em after the add) */
                        double aph = atan2(fa2[i], fa1[i]);
                        double mix = P.qp_phase * d1 / Em[i];
                        th2[i] = fmod(th2[i] + mix * wrap_pi(aph - th2[i])
                                      + 8.0 * TWO_PI, TWO_PI);
                    }
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
                    qatom_diag(0, atoms_w(w2e[i], w1e[i]), d2, i, Em[i]);
                    double bs = d2 * P.s_pull / (1.0 + P.s_pull);
                    backs_total += bs;
                    if (tag[i]) { ct_evap += d2; ct_backs += bs; }
                    Em[i] -= d2;
                    field_inject(i, d2 - bs);
                    Es[i] += bs;
                }
            }
            /* v91 RADIANCE (laws_V3r candidate A): graded SUB-CAP
             * emission — the steep-radiance law FORGE E1 measured
             * missing. x = holdings Em/cap (flight is on the links,
             * not radiatable from the cell — decision point R-D1).
             * det factor (rad_clock=0) cancels the beat slowdown so
             * the per-time rate is prop. to x^p_rad. Grains via the
             * atoms machinery (hbar-linearity preserved); routed as
             * evaporation (field inject + backsplash); folded into
             * the evap column of the tag ledger. k_rad=0 = V2g. */
            if (P.k_rad > 0 && Em[i] > 0) {
                double xr = Em[i] / P.cap;
                double comp = P.rad_clock ? 1.0 : (P.w2 / w2e[i]);
                double dr = P.k_rad * P.cap * pow(xr, P.p_rad) * comp;
                double eFr = A0eff * w2e[i] / TWO_PI;
                double eDr = A0eff * w1e[i] / TWO_PI;
                dr = atoms_fire(dr, eFr, eDr, &qcnvF[i]);
                dr = atoms_clamp(dr, Em[i], eFr, eDr, &qcnvF[i]);
                if (dr > 0) {
                    rad_total += dr;
                    qatom_diag(0, atoms_w(w2e[i], w1e[i]), dr, i, Em[i]);
                    double bsr = dr * P.s_pull / (1.0 + P.s_pull);
                    backs_total += bsr;
                    if (tag[i]) { ct_evap += dr; ct_backs += bsr; }
                    Em[i] -= dr;
                    field_inject(i, dr - bsr);
                    Es[i] += bsr;
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
                /* stable form: |ax/fs| <= 2 by construction (flux-weighted
                 * mean of unit-bounded vectors), but kappa_align*dt/fs
                 * overflows to inf when the only flux on a cell is
                 * subnormal dust (fs ~ 1e-311) — divide first */
                double w = P.kappa_align * dt;
                vx += w * (ax / fs); vy += w * (ay / fs); vz += w * (az / fs);
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

/* v91 cantus v1.1 diagnostic: a cell's amplitude = max incident
 * live-slot gauge memory (pure meter; the physical field is sgg_) */
static double cant_amp_of(int i)
{
    double a = 0;
    for (int q = cls_[i]; q < cls_[i + 1]; q++) {
        int s = clidx[q];
        if (sst[s] == S_FREE) continue;
        if (sgg_[s] > a) a = sgg_[s];
    }
    return a;
}

/* v91 registry diagnostics (REGISTRY.md §1.3 item 5, pure meter).
 * Match rho = 2*min/(sum) of the directed delivery rates. */
static double reg_rho(int s)
{
    double gross = rfp_[s] + rfm_[s];
    if (gross <= 0) return 0;
    double mn = rfp_[s] < rfm_[s] ? rfp_[s] : rfm_[s];
    return 2.0 * mn / gross;
}
static int dcmp_(const void *a, const void *b)
{
    double x = *(const double *)a, y = *(const double *)b;
    return x < y ? -1 : x > y ? 1 : 0;
}
/* class 0 = both endpoints tagged (the seeded body's own bonds);
 * class 1 = neither tagged (the bath). Live slots only. */
static int reg_stats(int cls, double *q25, double *q50, double *q75,
                     double *q90, double *grmed, double *flow)
{
    int n = 0, nf = 0;
    for (int s = 0; s < NSLOT; s++) {
        if (sst[s] != S_ALIVE) continue;
        int ti = tag[sli[s]] ? 1 : 0, tj = tag[slj[s]] ? 1 : 0;
        if (cls == 0 ? !(ti && tj) : (ti || tj)) continue;
        double gross = rfp_[s] + rfm_[s];
        regq1_[n] = reg_rho(s);
        regq2_[n] = gross;
        if (gross > 0) nf++;
        n++;
    }
    *q25 = *q50 = *q75 = *q90 = *grmed = *flow = 0;
    if (!n) return 0;
    qsort(regq1_, n, sizeof(double), dcmp_);
    qsort(regq2_, n, sizeof(double), dcmp_);
    *q25 = regq1_[(int)(0.25 * (n - 1))];
    *q50 = regq1_[(int)(0.50 * (n - 1))];
    *q75 = regq1_[(int)(0.75 * (n - 1))];
    *q90 = regq1_[(int)(0.90 * (n - 1))];
    *grmed = regq2_[(int)(0.50 * (n - 1))];
    *flow = (double)nf / n;
    return n;
}

/* v91 identity lane diag helpers (IDENTITY.md §1.3.7) */
static double par_q(double *a, int n, double f)
{ return n ? a[(int)(f * (n - 1))] : 0; }

static int par_ep_ages(double *q25, double *q50, double *q75)
{   /* ages of the LIVE episodes */
    int n = 0;
    for (int i = 0; i < NC; i++)
        if (cgid_[i] != 0) par_tmp[n++] = sim_t - cbirth_[i];
    *q25 = *q50 = *q75 = 0;
    if (!n) return 0;
    qsort(par_tmp, n, sizeof(double), dcmp_);
    *q25 = par_q(par_tmp, n, 0.25);
    *q50 = par_q(par_tmp, n, 0.50);
    *q75 = par_q(par_tmp, n, 0.75);
    return n;
}

static int par_death_ages(int cls, double *q25, double *q50, double *q90)
{   /* recent slot ages at death (ring window), by endpoint class */
    long long tot = par_agedn[cls];
    int n = tot > PAR_RING ? PAR_RING : (int)tot;
    for (int k = 0; k < n; k++) par_tmp[k] = par_aged[cls][k];
    *q25 = *q50 = *q90 = 0;
    if (!n) return 0;
    qsort(par_tmp, n, sizeof(double), dcmp_);
    *q25 = par_q(par_tmp, n, 0.25);
    *q50 = par_q(par_tmp, n, 0.50);
    *q90 = par_q(par_tmp, n, 0.90);
    return n;
}

static int par_stamp_ages(double *q50, double *mx, int *ntagpair)
{   /* ages of the LIVE bond identity stamps */
    int n = 0; *ntagpair = 0;
    for (int s = 0; s < NSLOT; s++) {
        if (sst[s] != S_ALIVE || pstampa_[s] == 0) continue;
        par_tmp[n++] = sim_t - pstampt_[s];
        if (tag[sli[s]] && tag[slj[s]]) (*ntagpair)++;
    }
    *q50 = *mx = 0;
    if (!n) return 0;
    qsort(par_tmp, n, sizeof(double), dcmp_);
    *q50 = par_q(par_tmp, n, 0.50);
    *mx = par_tmp[n - 1];
    return n;
}

/* v91 AMPLITUDE Phase M diag helper: per-class rho_coh quartiles
 * over slots with delivered magnitude (cls 0 = tag pair, 1 = bath) */
static int amp_stats(int cls, double *q25, double *q50, double *q75,
                     double *q90, double *magmed)
{
    int n = 0;
    for (int s = 0; s < NSLOT; s++) {
        if (sst[s] != S_ALIVE) continue;
        int ti = tag[sli[s]] ? 1 : 0, tj = tag[slj[s]] ? 1 : 0;
        if (cls == 0 ? !(ti && tj) : (ti || tj)) continue;
        if (amc_[s] <= 1e-12) continue;
        regq1_[n] = sqrt(amre_[s]*amre_[s] + amim_[s]*amim_[s]) / amc_[s];
        regq2_[n] = amc_[s];
        n++;
    }
    *q25 = *q50 = *q75 = *q90 = *magmed = 0;
    if (!n) return 0;
    qsort(regq1_, n, sizeof(double), dcmp_);
    qsort(regq2_, n, sizeof(double), dcmp_);
    *q25 = regq1_[(int)(0.25 * (n - 1))];
    *q50 = regq1_[(int)(0.50 * (n - 1))];
    *q75 = regq1_[(int)(0.75 * (n - 1))];
    *q90 = regq1_[(int)(0.90 * (n - 1))];
    *magmed = regq2_[(int)(0.50 * (n - 1))];
    return n;
}

static int bed_cmp(const void *a, const void *b)
{
    double d = *(const double *)a - *(const double *)b;
    return d < 0 ? -1 : d > 0 ? 1 : 0;
}

/* FLOW meters: live-slot sbed quartiles + max + grown count */
static void bed_stats(double *q25, double *q50, double *q75, double *mx,
                      int *ngrown, int *nlive)
{
    static double *buf = NULL;
    if (!buf) buf = malloc(NSLOT * sizeof(double));
    int n = 0, ng = 0;
    double m = 0;
    for (int s = 0; s < NSLOT; s++) {
        if (sst[s] == S_FREE || sA[s] <= 0) continue;
        buf[n++] = sbed[s];
        if (sbed[s] > m) m = sbed[s];
        if (fabs(sbed[s] - 1.0) > 0.05) ng++;
    }
    *nlive = n; *ngrown = ng; *mx = m;
    if (n == 0) { *q25 = *q50 = *q75 = 0; return; }
    qsort(buf, n, sizeof(double), bed_cmp);
    *q25 = buf[n / 4]; *q50 = buf[n / 2]; *q75 = buf[(3 * n) / 4];
}

static void bed_map_dump(void)
{
    int nout = 0;
    for (int s = 0; s < NSLOT && nout < 5000; s++) {
        if (sst[s] == S_FREE || sA[s] <= 0) continue;
        if (fabs(sbed[s] - 1.0) <= 0.05) continue;
        printf("# BEDMAP t=%.2f %d %d %.4f\n", sim_t, sli[s], slj[s], sbed[s]);
        nout++;
    }
}

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
    {
        double y = Eh_total - comp, t = s + y;
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

/* codec 1: columnar transpose (4-byte units, per record) + byte shuffle
 * + zstd. Per-chunk self-contained: random access and truncation safety
 * are preserved; measured ~0.46-0.55 ratio at ~175 MB/s (FCS.md). */
static unsigned char *fcs_s1 = NULL, *fcs_s2 = NULL, *fcs_s3 = NULL;
static size_t fcs_cap = 0;

static void fcs_ensure(size_t n)
{
    if (n <= fcs_cap) return;
    fcs_cap = n + n / 2 + 4096;
    fcs_s1 = realloc(fcs_s1, fcs_cap);
    fcs_s2 = realloc(fcs_s2, fcs_cap);
    fcs_s3 = realloc(fcs_s3, ZSTD_compressBound(fcs_cap));
}

static void fcs_colT(const unsigned char *p, unsigned char *out, size_t n,
                     int hdr, int stride)
{
    memcpy(out, p, hdr);
    size_t body = n - hdr;
    size_t nrec = body / (4 * (size_t)stride);
    const unsigned char *b = p + hdr;
    unsigned char *o = out + hdr;
    for (int c = 0; c < stride; c++)
        for (size_t i = 0; i < nrec; i++)
            memcpy(o + (c * nrec + i) * 4, b + (i * stride + c) * 4, 4);
}

static void fcs_shuffle4(const unsigned char *p, unsigned char *out, size_t n)
{
    size_t m = n / 4;
    for (size_t i = 0; i < m; i++) {
        out[i] = p[4*i];
        out[m+i] = p[4*i+1];
        out[2*m+i] = p[4*i+2];
        out[3*m+i] = p[4*i+3];
    }
    memcpy(out + 4*m, p + 4*m, n - 4*m);
}

/* write one chunk; CELL/LINK optionally wrapped in CMPD (codec 1) */
static void fcs_chunk_out(FILE *f, const char cc[4], const unsigned char *p,
                          size_t n, int hdr, int stride)
{
    if (P.snap_comp && stride > 0) {
        fcs_ensure(n);
        fcs_colT(p, fcs_s1, n, hdr, stride);
        fcs_shuffle4(fcs_s1, fcs_s2, n);
        size_t cn = ZSTD_compress(fcs_s3, ZSTD_compressBound(fcs_cap),
                                  fcs_s2, n, 1);
        if (!ZSTD_isError(cn) && cn + 13 < n) {
            uint64_t plen = 4 + 8 + 1 + cn;
            fwrite("CMPD", 1, 4, f); fcs_wu64(f, plen);
            fwrite(cc, 1, 4, f);
            fcs_wu64(f, (uint64_t)n);
            unsigned char codec = 1;
            fwrite(&codec, 1, 1, f);
            fwrite(fcs_s3, 1, cn, f);
            fflush(f);
            return;
        }
    }
    fwrite(cc, 1, 4, f);
    fcs_wu64(f, (uint64_t)n);
    fwrite(p, 1, n, f);
    fflush(f);
}

static void fcs_cell_frame(FILE *f)
{
    size_t clen = 8 + 8 + 4 + (size_t)NC * FCS_NCELLCOLS * 4;
    unsigned char *cp = malloc(clen);
    size_t o = 0;
    memcpy(cp + o, &sim_t, 8); o += 8;
    memcpy(cp + o, &P.L, 8); o += 8;
    uint32_t nc32 = (uint32_t)NC;
    memcpy(cp + o, &nc32, 4); o += 4;
    for (int i = 0; i < NC; i++) {
        float rec[FCS_NCELLCOLS];
        rec[0] = (float)px_[i]; rec[1] = (float)py_[i]; rec[2] = (float)pz_[i];
        rec[3] = (float)cr[i]; rec[4] = (float)Es[i]; rec[5] = (float)Em[i];
        rec[6] = (float)Ee[i]; rec[7] = (float)((Em[i] + flload[i]) / P.cap);
        rec[8] = (float)tag[i];
        rec[9] = (float)fa1[i]; rec[10] = (float)fa2[i]; rec[11] = (float)th2[i];
        memcpy(cp + o, rec, 4 * FCS_NCELLCOLS); o += 4 * FCS_NCELLCOLS;
    }
    fcs_chunk_out(f, "CELL", cp, clen, 20, FCS_NCELLCOLS);
    free(cp);

    uint32_t nl = 0;
    for (int s = 0; s < NSLOT; s++) if (sst[s] != S_FREE) nl++;
    size_t llen = 8 + 4 + (size_t)nl * (8 + FCS_NLINKCOLS * 4);
    unsigned char *lp = malloc(llen);
    o = 0;
    memcpy(lp + o, &sim_t, 8); o += 8;
    memcpy(lp + o, &nl, 4); o += 4;
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
        memcpy(lp + o, ij, 8); o += 8;
        memcpy(lp + o, lr, 4 * FCS_NLINKCOLS); o += 4 * FCS_NLINKCOLS;
    }
    fcs_chunk_out(f, "LINK", lp, llen, 12, 2 + FCS_NLINKCOLS);
    free(lp);
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
    if (P.sect_meter) {
        static const char *XC[] = {"th","E","n"};
        double rows[SECT_MAX * 3];
        for (int k = 0; k < P.sect_n; k++) {
            double thc = k * TWO_PI / P.sect_n;
            if (thc > M_PI) thc -= TWO_PI;
            rows[3*k] = thc;
            rows[3*k+1] = sectE[k];
            rows[3*k+2] = sectN[k];
        }
        fcs_anlz_table(f, "sector", XC, 3, rows, P.sect_n);
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
static int truss_e[128][2];
static double truss_dstar[128];
static double truss_shape0[3];
/* FCQ triangle registry: declared chart (p,q) per cycle edge, cycle
 * direction, and the second triangle for exp=tri2 */
static int tri_on = 0, tri_v[2][3];
/* COMPOSITE inter-ring boundary edges (directed, with charts) */
static int rint_n = 0;
static int rint_e[8][2];
static signed char rint_p[8], rint_q[8];
static double rint_dstar[8];
static int rings_nring = 0, rings_v0[6], rings_mode[6]; /* 1=molecule 0=droplet */
static double rings_w[6];
static int rings_grp[6];
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
    printf("# v91 radiance (laws_V3r candidate A): k_rad=%g p_rad=%g rad_clock=%d\n",
           P.k_rad, P.p_rad, P.rad_clock);
    printf("# v91 cantus (coherent-channel candidate B): k_cant=%g k_tune=%g cant_tau=%g cant_seed=%d cant_grow=%d\n",
           P.k_cant, P.k_tune, P.cant_tau, P.cant_seed, P.cant_grow);
    printf("# v91 registry (exchange-registry lane, REGISTRY.md): reg_tau=%g reg_gate=%d reg_f0=%g\n",
           P.reg_tau, P.reg_gate, P.reg_f0);
    printf("# v91 identity (parcel-identity lane, IDENTITY.md): par_tau=%g par_gate=%d par_form=%d par_lo=%g par_hi=%g par_mature=%g\n",
           P.par_tau, P.par_gate, P.par_form, P.par_lo, P.par_hi, P.par_mature);
    printf("# v91 workfn (emergent-threshold lane, WORKFN.md): wf_on=%d wf_floor=%g wf_far=%g\n",
           P.wf_on, P.wf_floor, P.wf_far);
    printf("# v91 amplitude (field-side-identity lane Phase M, AMPLITUDE.md): amp_tau=%g\n",
           P.amp_tau);
    printf("# v92 amplitude Phase L lane L-1 (L0_DESIGN.md): amp_drv=%g amp_mmin=%g\n",
           P.amp_drv, P.amp_mmin);
    printf("# v93 UNITARY DENSE CHANNEL (v93/README.md PART II): amp_nat=%g\n",
           P.amp_nat);
    printf("# v93 arg(psi) door (v93/README.md §II.7): amp_door=%g\n",
           P.amp_door);
    printf("# QUENCH-2 apparatus: conf_r=%g conf_gap=%g conf_th=%g conf_pinw=%g spin_m=%d imp_k=%g qp_phase=%g\n",
           P.conf_r, P.conf_gap, P.conf_th, P.conf_pinw, P.spin_m, P.imp_k, P.qp_phase);
    printf("# HORIZON apparatus (HORIZON.md): bh_r=%g bh_k=%g bh_sep=%g\n",
           P.bh_r, P.bh_k, P.bh_sep);
    printf("# FLOW apparatus (FLOW.md): bed_k=%g bed_tau=%g\n",
           P.bed_k, P.bed_tau);
    if (P.par_gate && P.par_tau <= 0)
        printf("# CONFIG ERROR: par_gate=1 with par_tau=0 — r_id == 0, gauge dark everywhere\n");
    if (P.par_gate && P.reg_gate) {
        printf("# CONFIG ERROR: par_gate and reg_gate both set — refusing to run law arms with two gates\n");
        exit(1);
    }
    printf("# GEOMETRY (apparatus): cfac=%g k_rep=%g mob_geo=%g kappa_bond=%g freeze_geo=%d\n",
           P.cfac, P.k_rep, P.mob_geo, P.kappa_bond, P.freeze_geo);
    printf("# bath=%d bath_frac=%g jam_sweeps=%d jam_k=%g L=%g dt=%g T=%g seed=%lu diag_every=%d\n",
           P.bath, P.bath_frac, P.jam_sweeps, P.jam_k, P.L, P.dt, P.T, P.seed,
           P.diag_every);

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
    if (P.grad_r0 > 0 && P.grad_frac < 1.0 && nb > 0) {
        /* FORGE: thin the dart bath inside r0 of the box centre — a
         * matter-free foam-density structure (dark-lensing apparatus) */
        double gc = 0.5 * P.L;
        int m = 0, nin = 0;
        for (int i = 0; i < nb; i++) {
            double dx = wr(bx[i]-gc), dy = wr(by[i]-gc), dz = wr(bz[i]-gc);
            int inside = dx*dx + dy*dy + dz*dz < P.grad_r0 * P.grad_r0;
            if (inside) nin++;
            if (!inside || frand() < P.grad_frac) {
                bx[m] = bx[i]; by[m] = by[i]; bz[m] = bz[i];
                m++;
            }
        }
        printf("# GRAD carved: removed %d of %d cells inside r=%g (frac=%g)\n",
               nb - m, nin, P.grad_r0, P.grad_frac);
        nb = m;
    }

    int extra = 0;
    if (!strcmp(P.exp, "pair") && !P.bath) extra = 2;
    if (!strcmp(P.exp, "ring")) extra = P.ring_n;
    if (!strcmp(P.exp, "oct")) extra = 6;
    if (!strcmp(P.exp, "tri")) extra = 3;
    if (!strcmp(P.exp, "tri2")) extra = 6;
    if (!strcmp(P.exp, "rings"))
        extra = (P.rings_kind == 0 ? 2 : P.rings_kind == 3 ? 6 : 3) * P.rings_nv;

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
        if (P.kx != 0 && !strcmp(P.exp, "ring")) {
            /* MOTION #31: e3b-style tilt ON a bound object (drive
             * detunes the bonds by design; pre-registered null) */
            for (int k = 0; k < nvo; k++) {
                int i = base_i + k;
                th2[i] = fmod(th2[i] - P.kx * wr(px_[i] - cx0) + 8.0 * TWO_PI, TWO_PI);
            }
            printf("# SEED ring tilt: kx=%g applied on top of the chain\n", P.kx);
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
            /* e3b tilt: phase ramp along x (cellfab.c:1750 convention);
             * OR seed_mw: a MATTER azimuthal winding m*phi about the blob
             * centre (v93 item 4 -- hand-seeded matter vortex, no field/door,
             * isolates hold-vs-imprint for retention). */
            if (P.seed_mw != 0)
                th2[i] = fmod(P.seed_mw * atan2(dy, dx) + 8.0 * TWO_PI, TWO_PI);
            else
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
            int kind_t = (t == 1 && P.tri2_kind >= 0) ? P.tri2_kind
                                                      : P.tri_kind;
            int v0 = base_i + 3 * t;
            for (int kk = 0; kk < 3; kk++) {
                int i = v0 + kk;
                seed_cell_defaults(i);
                n1x[i] = 0; n1y[i] = 0; n1z[i] = 1;  /* coplanar, gpl=1 */
                n2x[i] = 0; n2y[i] = 0; n2z[i] = 1;
                tag[i] = 1;
                tri_v[t][kk] = i;
            }
            if (kind_t == 0) {
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
        if (ntri == 2 && P.tri2_kind >= 0)
            printf("# SEED tri: MIXED pair T0=%s T1=%s\n",
                   P.tri_kind ? "UDD" : "UUD", P.tri2_kind ? "UDD" : "UUD");
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
            if (P.spin_m) tilt += P.spin_m * atan2(dy, dx);
            if (P.imp_k > 0) tilt += P.imp_k * sqrt(dx*dx + dy*dy);
            fa1[i] += sqrt(P.amp * g) * cos(tilt);
            fa2[i] += sqrt(P.amp * g) * sin(tilt);
            Ee[i] = fa1[i]*fa1[i] + fa2[i]*fa2[i];
        }
        for (int i = 0; i < NC; i++) normals_transverse(i, 1, 0, 0);
        printf("# SEED pulse: amp=%g sigma=%g kx=%g (prealigned transverse)\n",
               P.amp, P.sigma, P.kx);
    }
    else if (!strcmp(P.exp, "blob2")) {
        /* MOTION #28: two dense blobs with opposing tilts (approach),
         * per-blob tags for COM/exchange meters; total dense COM =
         * the Delta-p bookkeeping (centre-of-energy). */
        double cxA = fold(cx0 - 0.5 * P.blob2_sep), cxB = fold(cx0 + 0.5 * P.blob2_sep);
        int nA = 0, nB = 0;
        for (int i = 0; i < NC; i++) {
            double dxa = wr(px_[i]-cxA), dya = wr(py_[i]-cy0), dza = wr(pz_[i]-cz0);
            double dxb = wr(px_[i]-cxB), dyb = wr(py_[i]-cy0), dzb = wr(pz_[i]-cz0);
            double ga = P.amp * exp(-(dxa*dxa+dya*dya+dza*dza)/(2.0*P.sigma*P.sigma));
            double gb = P.amp * exp(-(dxb*dxb+dyb*dyb+dzb*dzb)/(2.0*P.sigma*P.sigma));
            double add = ga + gb;
            if (add < 1e-4) continue;
            if (Em[i] + add > 0.95 * P.cap) add = 0.95 * P.cap - Em[i];
            if (add <= 0) continue;
            Em[i] += add;
            double pull = P.s_pull * add;
            double avail = Es[i] - P.es_floor;
            if (pull > avail) pull = avail > 0 ? avail : 0;
            Es[i] -= pull;
            Em[i] += pull;
            /* tilt toward each other; phase ramp of the DOMINANT blob */
            if (P.blob2_kx != 0) {
                if (ga >= gb) th2[i] = fmod(-(P.blob2_kx * dxa) + 8.0*TWO_PI, TWO_PI);
                else          th2[i] = fmod(+(P.blob2_kx * dxb) + 8.0*TWO_PI, TWO_PI);
            }
            if (ga >= gb && ga > 0.05 * P.amp) { tag[i] = 1; nA++; }
            else if (gb > 0.05 * P.amp) { tag[i] = 2; nB++; }
        }
        if (P.blob2_kx != 0)
            for (int i = 0; i < NC; i++) normals_transverse(i, 1, 0, 0);
        printf("# SEED blob2: amp=%g sigma=%g sep=%g kx=%g tagged A=%d B=%d\n",
               P.amp, P.sigma, P.blob2_sep, P.blob2_kx, nA, nB);
    }
    else if (!strcmp(P.exp, "rings")) {
        /* COMPOSITE.md CQ0/CQ3b/CQ9: groups of rings as many-celled
         * quarks. kind 0 pair, 1 UUD, 2 UDD, 3 = TWO UUD nucleons at
         * rings_sep. Facing cells loaded at x - rings_xcomp (flight
         * compensation, CQ3b). Bath: spherical carve per ring. */
        int nv = P.rings_nv;
        int ngrp = P.rings_kind == 3 ? 2 : 1;
        int npg = P.rings_kind == 0 ? 2 : 3;
        int nring = ngrp * npg;
        double xU = P.rings_xU;
        double xD = P.rings_xD < 0 ? (0.5 + 1.8 * xU) / 1.2 : P.rings_xD;
        double wU = P.w2 / (1.0 + P.q_detune * xU);
        double wD = P.w2 / (1.0 + P.q_detune * xD);
        double dUD = 2.0 * TWO_PI / (2.0 * wU + 3.0 * wD);
        double xs[6], ws[6];
        for (int k = 0; k < nring; k++) {
            int lk = k % npg;
            double x = xU;
            if (P.rings_kind == 2 && lk >= 1) x = xD;
            if ((P.rings_kind == 1 || P.rings_kind == 3) && lk == 2) x = xD;
            xs[k] = x;
            ws[k] = P.w2 / (1.0 + P.q_detune * x);
        }
        double rload[6], dint[6], R[6];
        int mode[6];
        for (int k = 0; k < nring; k++) {
            double pull = P.s_pull * (xs[k] * P.cap / (1.0 + P.s_pull));
            double cap2 = P.e_s0 - P.es_floor;
            double es = P.e_s0 - (pull < cap2 ? pull : cap2);
            rload[k] = P.r0 * cbrt(es / P.e_s0);
            double rung = M_PI / ws[k];
            double contact = 2.0 * rload[k];
            if (rung < contact * P.cfac) { dint[k] = rung; mode[k] = 1; }
            else { dint[k] = 0.98 * contact; mode[k] = 0; }
            R[k] = dint[k] / (2.0 * sin(M_PI / nv));
        }
        /* group-local centers (SSS), then group offsets */
        double cxr[6], cyr[6];
        for (int g = 0; g < ngrp; g++) {
            int b0 = g * npg;
            double gx = ngrp == 2 ? (g == 0 ? -0.5 : 0.5) * P.rings_sep : 0.0;
            if (npg == 2) {
                double bg = M_PI / ws[b0] + P.rings_gapoff;
                double cdist = R[b0] + R[b0+1] + bg;
                cxr[b0] = gx - 0.5 * cdist; cyr[b0] = 0;
                cxr[b0+1] = gx + 0.5 * cdist; cyr[b0+1] = 0;
            } else {
                double bg01 = (fabs(ws[b0] - ws[b0+1]) < 1e-12 ? M_PI / ws[b0] : dUD) + P.rings_gapoff;
                double bg12 = (fabs(ws[b0+1] - ws[b0+2]) < 1e-12 ? M_PI / ws[b0+1] : dUD) + P.rings_gapoff;
                double bg02 = (fabs(ws[b0] - ws[b0+2]) < 1e-12 ? M_PI / ws[b0] : dUD) + P.rings_gapoff;
                double c01 = R[b0] + R[b0+1] + bg01;
                double c12 = R[b0+1] + R[b0+2] + bg12;
                double c02 = R[b0] + R[b0+2] + bg02;
                double x0 = 0, y0 = 0, x1 = c01, y1 = 0;
                double x2 = (c01 * c01 + c02 * c02 - c12 * c12) / (2.0 * c01);
                double y2s = c02 * c02 - x2 * x2;
                double y2 = y2s > 0 ? sqrt(y2s) : 0;
                double mx = (x0 + x1 + x2) / 3.0, my = (y0 + y1 + y2) / 3.0;
                cxr[b0] = gx + x0 - mx; cyr[b0] = y0 - my;
                cxr[b0+1] = gx + x1 - mx; cyr[b0+1] = y1 - my;
                cxr[b0+2] = gx + x2 - mx; cyr[b0+2] = y2 - my;
            }
        }
        /* bath: spherical carve per ring */
        if (P.bath && nb > 0) {
            int m = 0, rem = 0;
            for (int i = 0; i < nb; i++) {
                int keep = 1;
                for (int k = 0; k < nring && keep; k++) {
                    double dx = wr(px_[i] - (cx0 + cxr[k]));
                    double dy = wr(py_[i] - (cy0 + cyr[k]));
                    double dz = wr(pz_[i] - cz0);
                    double clear = R[k] + 0.8;
                    if (dx*dx + dy*dy + dz*dz < clear * clear) keep = 0;
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
                } else rem++;
            }
            printf("# CARVE rings: removed %d bath cells\n", rem);
            nb = m;
            NC = nb + extra;
        }
        rings_nring = nring;
        int base_i = nb;
        truss_n = 0;
        for (int k = 0; k < nring; k++) {
            rings_v0[k] = base_i + k * nv;
            rings_mode[k] = mode[k];
            rings_w[k] = ws[k];
            rings_grp[k] = k / npg;
            int g = k / npg, lk = k % npg;
            int nbr = g * npg + (lk + 1) % npg;
            double az = atan2(cyr[nbr] - cyr[k], cxr[nbr] - cxr[k]);
            for (int q = 0; q < nv; q++) {
                int i = rings_v0[k] + q;
                seed_cell_defaults(i);
                double a = az + TWO_PI * q / nv;
                px_[i] = fold(cx0 + cxr[k] + R[k] * cos(a));
                py_[i] = fold(cy0 + cyr[k] + R[k] * sin(a));
                pz_[i] = cz0;
                n1x[i] = 0; n1y[i] = 0; n1z[i] = 1;
                n2x[i] = 0; n2y[i] = 0; n2z[i] = 1;
                tag[i] = 1;
            }
            for (int q = 0; q < nv; q++) {
                truss_e[truss_n][0] = rings_v0[k] + q;
                truss_e[truss_n][1] = rings_v0[k] + (q + 1) % nv;
                truss_dstar[truss_n] = dint[k];
                truss_n++;
            }
        }
        /* facing vertices (within group) */
        int face[6][6];
        for (int a = 0; a < nring; a++)
            for (int b2 = 0; b2 < nring; b2++) {
                face[a][b2] = -1;
                if (a == b2 || rings_grp[a] != rings_grp[b2]) continue;
                double az = atan2(cyr[b2] - cyr[a], cxr[b2] - cxr[a]);
                int g = rings_grp[a], la = a % npg;
                int nbr0 = g * npg + (la + 1) % npg;
                double az0 = atan2(cyr[nbr0] - cyr[a], cxr[nbr0] - cxr[a]);
                double rel = az - az0;
                int q = (int)floor(rel / (TWO_PI / nv) + 0.5);
                q = ((q % nv) + nv) % nv;
                face[a][b2] = rings_v0[a] + q;
            }
        /* loads: per-flavor WHOLE-RING flight compensation (CQ3b v2 —
         * the facing-only variant measurably broke the fifth arithmetic;
         * the dynamic xbar excess is ring-wide internal flight, so the
         * compensation is ring-wide and per flavor) */
        for (int k = 0; k < nring; k++) {
            double xl = xs[k];
            if (fabs(ws[k] - wU) < 1e-12) xl -= P.rings_xcomp;
            else xl -= P.rings_xcompD;
            if (xl < 0) xl = 0;
            for (int q = 0; q < nv; q++)
                load_voice(rings_v0[k] + q, xl);
        }
        /* inter-ring boundary edges (within groups) */
        rint_n = 0;
        for (int g = 0; g < ngrp; g++) {
            int nedge = npg == 2 ? 1 : 3;
            for (int e = 0; e < nedge; e++) {
                int a = g * npg + e;
                int b2 = g * npg + (e + 1) % npg;
                int ia = face[a][b2], ib = face[b2][a];
                rint_e[rint_n][0] = ia; rint_e[rint_n][1] = ib;
                if (fabs(ws[a] - ws[b2]) < 1e-12) { rint_p[rint_n] = 1; rint_q[rint_n] = 1; rint_dstar[rint_n] = M_PI / ws[a]; }
                else if (ws[a] > ws[b2]) { rint_p[rint_n] = 3; rint_q[rint_n] = 2; rint_dstar[rint_n] = dUD; }
                else { rint_p[rint_n] = 2; rint_q[rint_n] = 3; rint_dstar[rint_n] = dUD; }
                rint_n++;
            }
        }
        /* phases: per-group chains; branch on the first U->D edge */
        if (P.seedlock) {
            for (int g = 0; g < ngrp; g++) {
                int b0 = g * npg;
                th2[rings_v0[b0]] = frand() * TWO_PI;
                int nedge = npg == 2 ? 1 : 3;
                for (int k = b0; k < b0 + npg; k++) {
                    if (k > b0) {
                        int e = g * nedge + (k - b0 - 1);
                        int ia = rint_e[e][0], ib = rint_e[e][1];
                        double dx = wr(px_[ib] - px_[ia]), dy = wr(py_[ib] - py_[ia]);
                        double d = sqrt(dx * dx + dy * dy);
                        if (rint_p[e] == 1 && rint_q[e] == 1)
                            th2[ib] = fmod(th2[ia] - ws[k-1] * d / P.C + 8.0 * TWO_PI, TWO_PI);
                        else if (rint_p[e] == 3)
                            th2[ib] = fmod((2.0 * (th2[ia] - ws[k-1] * d / P.C)
                                            + TWO_PI * P.rings_branch) / 3.0
                                           + 8.0 * TWO_PI, TWO_PI);
                        else
                            th2[ib] = fmod((3.0 * (th2[ia] - ws[k-1] * d / P.C)) / 2.0
                                           + 8.0 * TWO_PI, TWO_PI);
                    }
                    int start = k == b0 ? rings_v0[b0]
                              : rint_e[g * nedge + (k - b0 - 1)][1];
                    int sq = start - rings_v0[k];
                    for (int step = 1; step < nv; step++) {
                        int qa = (sq + step - 1) % nv, qb = (sq + step) % nv;
                        int ia = rings_v0[k] + qa, ib = rings_v0[k] + qb;
                        double dx = wr(px_[ib] - px_[ia]), dy = wr(py_[ib] - py_[ia]);
                        double d = sqrt(dx * dx + dy * dy);
                        th2[ib] = fmod(th2[ia] - ws[k] * d / P.C + 8.0 * TWO_PI, TWO_PI);
                    }
                }
            }
        }
        printf("# SEED rings: kind=%d ngrp=%d nv=%d xU=%.4f xD=%.4f wU=%.4f wD=%.4f xcomp=%.4f sep=%.2f\n",
               P.rings_kind, ngrp, nv, xU, xD, wU, wD, P.rings_xcomp, P.rings_sep);
        for (int k = 0; k < nring; k++)
            printf("# SEED ring%d: grp=%d mode=%s x=%.4f w=%.4f d_int=%.4f R=%.4f\n",
                   k, rings_grp[k], mode[k] ? "MOLECULE" : "DROPLET", xs[k], ws[k], dint[k], R[k]);
        for (int e = 0; e < rint_n; e++)
            printf("# SEED rint%d: %d-%d p:q=%d:%d d*=%.6f\n",
                   e, rint_e[e][0], rint_e[e][1], rint_p[e], rint_q[e], rint_dstar[e]);
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
            if (P.conf_r > 0) {
                /* QUENCH-2 cavity: carve the shell annulus except the
                 * leak gap (half-angle conf_gap about +x) */
                double rdx = wr(px_[i] - cx0), rdy = wr(py_[i] - cy0);
                double rr = sqrt(rdx*rdx + rdy*rdy);
                if (fabs(rr - P.conf_r) < 0.5 * P.conf_th
                    && fabs(atan2(rdy, rdx)) > P.conf_gap) keep = 0;
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
        if (P.conf_r > 0) {
            int ncp = 0;
            for (int i = 0; i < NC; i++) {
                double rdx = wr(px_[i] - cx0), rdy = wr(py_[i] - cy0);
                double rr = sqrt(rdx*rdx + rdy*rdy);
                if (fabs(rr - P.conf_r) < P.conf_pinw
                    && fabs(atan2(rdy, rdx)) > P.conf_gap) { pin[i] = 1; ncp++; }
            }
            printf("# PIN conf cavity: %d cells, r=%g th=%g gap=%g\n",
                   ncp, P.conf_r, P.conf_th, P.conf_gap);
        }
        /* quasi-plane packet: Gaussian sigma_x = P.sigma, sigma_y = slit_sy */
        for (int i = 0; i < NC; i++) {
            double dx = wr(px_[i] - P.slit_srcx), dy = wr(py_[i] - cy0);
            double g = exp(-dx*dx / (2.0 * P.sigma * P.sigma)
                           - dy*dy / (2.0 * P.slit_sy * P.slit_sy));
            if (g < 1e-8) continue;
            double tilt = -(P.kx * dx);
            if (P.spin_m) tilt += P.spin_m * atan2(dy, dx);
            if (P.imp_k > 0) tilt += P.imp_k * sqrt(dx*dx + dy*dy);
            fa1[i] += sqrt(P.amp * g) * cos(tilt);
            fa2[i] += sqrt(P.amp * g) * sin(tilt);
            Ee[i] = fa1[i]*fa1[i] + fa2[i]*fa2[i];
        }
        for (int i = 0; i < NC; i++) normals_transverse(i, 1, 0, 0);
        if (P.slit_obj) {
            /* MOTION #29/XSEC: dense object in the beam; obj_y (<0 = box
             * centre) sets the impact parameter b = obj_y - cy0 */
            int nobj = 0;
            double oy = P.obj_y < 0 ? cy0 : P.obj_y;
            for (int i = 0; i < NC; i++) {
                double dx = wr(px_[i]-cx0), dy = wr(py_[i]-oy);
                double add = P.obj_amp * exp(-(dx*dx+dy*dy)/(2.0*P.obj_sigma*P.obj_sigma));
                if (add < 1e-4) continue;
                if (Em[i] + add > 0.95 * P.cap) add = 0.95 * P.cap - Em[i];
                if (add <= 0) continue;
                Em[i] += add;
                double pull = P.s_pull * add;
                double avail = Es[i] - P.es_floor;
                if (pull > avail) pull = avail > 0 ? avail : 0;
                Es[i] -= pull;
                Em[i] += pull;
                if (add > 0.05 * P.obj_amp) { tag[i] = 1; nobj++; }
            }
            printf("# SEED slit_obj: amp=%g sigma=%g y=%.2f tagged=%d (dense occulter)\n",
                   P.obj_amp, P.obj_sigma, oy, nobj);
        }
        if (P.slit_clicks) {
            int nsc = 0;
            for (int i = 0; i < NC; i++)
                if (fabs(wr(px_[i] - P.slit_screenx)) <= 0.75) { scond[i] = 1; nsc++; }
            printf("# SEED slit_clicks: %d screen cells condensation-active (e_cond=0 there)\n", nsc);
        }
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

    if (P.tag_r > 0) {
        /* FORGE: region tag WITHOUT matter — the convtag ledger then
         * attributes conversions to the region (forging meter) */
        int nt = 0;
        for (int i = 0; i < NC; i++) {
            double dx = wr(px_[i]-cx0), dy = wr(py_[i]-cy0), dz = wr(pz_[i]-cz0);
            if (dx*dx + dy*dy + dz*dz < P.tag_r * P.tag_r) { tag[i] = 1; nt++; }
        }
        printf("# TAG region: %d cells within r=%g of centre (no load added)\n",
               nt, P.tag_r);
    }

    /* v91 cantus init: holdings memory starts AT the seeded part
     * (no startup retune transient) */
    for (int i = 0; i < NC; i++) cxl_[i] = Em[i] / P.cap;

    /* initial radii + topology + ledger */
    for (int i = 0; i < NC; i++) {
        double ratio = Es[i] / P.e_s0;
        cr[i] = cr0[i] * cbrt(ratio > 0 ? ratio : 0);
    }
    topo_refresh();
    /* v91 cantus v1.1: cant_seed=1 arms the BONDS between tagged
     * voices at full gauge memory (apparatus; default = self-growth) */
    if (P.cant_seed) {
        int nseed = 0;
        for (int s = 0; s < NSLOT; s++) {
            if (sst[s] == S_FREE) continue;
            if (tag[sli[s]] && tag[slj[s]]) { sgg_[s] = 1.0; nseed++; }
        }
        printf("# SEED cantus: sgg=1 on %d tagged-pair slots\n", nseed);
    }
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

    if (P.sect_meter) {
        sect_cx = P.sect_x < 0 ? cx0 : P.sect_x;
        sect_cy = P.sect_y < 0 ? cy0 : P.sect_y;
        if (P.sect_n > SECT_MAX) P.sect_n = SECT_MAX;
        if (P.sect_n < 1) P.sect_n = 1;
        printf("# SECT meter: centre=(%.2f,%.2f) r=[%g,%g) n=%d gate=[%g,%g]\n",
               sect_cx, sect_cy, P.sect_r0, P.sect_r1, P.sect_n, P.sect_t0, P.sect_t1);
    }

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
            if (!strcmp(P.exp, "blob2")) {
                double ax = 0, aw = 0, bx2 = 0, bw = 0, tx = 0, tw = 0;
                for (int i2 = 0; i2 < NC; i2++) {
                    double w = Em[i2] + flload[i2];
                    double xx = cenx + wr(px_[i2] - cenx);
                    if (w > 1e-12) { tx += w * xx; tw += w; }
                    if (tag[i2] == 1 && w > 1e-12) { ax += w * xx; aw += w; }
                    if (tag[i2] == 2 && w > 1e-12) { bx2 += w * xx; bw += w; }
                }
                double comA = aw > 0 ? ax / aw : 0, comB = bw > 0 ? bx2 / bw : 0;
                printf(" | A=%.4f B=%.4f sepx=%.4f tot=%.5f EmA=%.3f EmB=%.3f",
                       comA, comB, comB - comA, tw > 0 ? tx / tw : 0, aw, bw);
                if (nfsamp < 512) {
                    fs_t[nfsamp] = sim_t;
                    b2_ax[nfsamp] = comA; b2_bx[nfsamp] = comB;
                    b2_tx[nfsamp] = tw > 0 ? tx / tw : 0;
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
            if (rings_nring == 6) {
                int base = rings_v0[0];
                int ncross = 0; double dmin2 = 1e300;
                for (int s = 0; s < NSLOT; s++) {
                    if (sst[s] == S_FREE || sA[s] <= 0) continue;
                    int i = sli[s], j = slj[s];
                    if (i < base || j < base) continue;
                    int gi = rings_grp[(i - base) / P.rings_nv];
                    int gj = rings_grp[(j - base) / P.rings_nv];
                    if (gi != gj) { ncross++; if (sd[s] < dmin2) dmin2 = sd[s]; }
                }
                double g1x = 0, g1y = 0, g2x = 0, g2y = 0; int n1_ = 0, n2_ = 0;
                for (int k = 0; k < 6; k++)
                    for (int q = 0; q < P.rings_nv; q++) {
                        int i = rings_v0[k] + q;
                        double xx = cenx + wr(px_[i] - cenx), yy = ceny + wr(py_[i] - ceny);
                        if (rings_grp[k] == 0) { g1x += xx; g1y += yy; n1_++; }
                        else { g2x += xx; g2y += yy; n2_++; }
                    }
                double sepg = sqrt((g2x/n2_ - g1x/n1_)*(g2x/n2_ - g1x/n1_)
                                   + (g2y/n2_ - g1y/n1_)*(g2y/n2_ - g1y/n1_));
                printf(" | NG ncross=%d dmin=%.3f sep=%.4f", ncross,
                       ncross ? dmin2 : -1.0, sepg);
            }
            if (rint_n > 0) {
                printf(" | IR");
                for (int e = 0; e < rint_n; e++) {
                    int i = rint_e[e][0], j = rint_e[e][1];
                    int s = hfind(i < j ? i : j, i < j ? j : i);
                    if (s < 0 || sst[s] == S_FREE) { printf(" e%d=DEAD", e); continue; }
                    double p = rint_p[e], q = rint_q[e];
                    double psi = wrap_pi(q * th2[i] - q * w2e[i] * sd[s] / P.C - p * th2[j]);
                    double gf = gate_of(q * th2[i] - q * w2e[i] * sd[s] / P.C - p * th2[j]);
                    double gb = gate_of(p * th2[j] - p * w2e[j] * sd[s] / P.C - q * th2[i]);
                    printf(" e%d:psi=%+.3f gg=%.3f d=%.4f", e, psi, gf * gb, sd[s]);
                }
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
            if (P.p1_meter) {
                double tx = p1sp[0] + p1fl[0] + p1fd[0] + p1gm[0];
                double ty = p1sp[1] + p1fl[1] + p1fd[1] + p1gm[1];
                double tz = p1sp[2] + p1fl[2] + p1fd[2] + p1gm[2];
                printf("# P1 t=%.2f tot=(%+.6e,%+.6e,%+.6e) sp=%+.3e fl=%+.3e fd=%+.3e gm=%+.3e (x)\n",
                       sim_t, tx, ty, tz, p1sp[0], p1fl[0], p1fd[0], p1gm[0]);
            }
            if (P.slit_obj || P.convtag)
                printf("# CONVTAG t=%.2f rough=%.6f cond=%.6f evap=%.6f backs=%.6f\n",
                       sim_t, ct_rough, ct_cond, ct_evap, ct_backs);
            if (P.k_cant > 0 || P.k_tune > 0) {
                double at = 0, am = 0, xt = 0; int ntg = 0, nl = 0;
                for (int i2 = 0; i2 < NC; i2++) {
                    double av = cant_amp_of(i2);
                    if (av > am) am = av;
                    if (av > 0.5) nl++;
                    if (tag[i2]) { at += av; xt += cxl_[i2]; ntg++; }
                }
                printf("# CANT t=%.2f a_tag=%.4f a_max=%.4f xl_tag=%.4f nlock=%d tune=%.6f\n",
                       sim_t, ntg ? at / ntg : 0, am, ntg ? xt / ntg : 0,
                       nl, tune_total);
            }
            if (P.reg_tau > 0) {
                double a25, a50, a75, a90, agr, afl, b25, b50, b75, b90, bgr, bfl;
                int ntp = reg_stats(0, &a25, &a50, &a75, &a90, &agr, &afl);
                int nba = reg_stats(1, &b25, &b50, &b75, &b90, &bgr, &bfl);
                printf("# REG t=%.2f tp n=%d rho=[%.3f %.3f %.3f %.3f] gr=%.5f fl=%.3f | ba n=%d rho=[%.3f %.3f %.3f %.3f] gr=%.5f fl=%.3f\n",
                       sim_t, ntp, a25, a50, a75, a90, agr, afl,
                       nba, b25, b50, b75, b90, bgr, bfl);
            }
            if (P.par_tau > 0) {
                double ea25, ea50, ea75, d25t, d50t, d90t, d25b, d50b, d90b;
                double sa50, samx;
                int nst, nstp, ngid = 0;
                for (int i2 = 0; i2 < NC; i2++) if (cgid_[i2] != 0) ngid++;
                par_ep_ages(&ea25, &ea50, &ea75);
                int ndt = par_death_ages(0, &d25t, &d50t, &d90t);
                int ndb = par_death_ages(1, &d25b, &d50b, &d90b);
                nst = par_stamp_ages(&sa50, &samx, &nstp);
                printf("# PAR t=%.2f gids n=%d mint=%lld ret=%lld age=[%.0f %.0f %.0f] | sdeath tp n=%d [%.1f %.1f %.1f] ba n=%d [%.1f %.1f %.1f] | stamps n=%d tp=%d age50=%.0f max=%.0f | idfrac=%.4f\n",
                       sim_t, ngid, par_mints, par_retires, ea25, ea50, ea75,
                       ndt, d25t, d50t, d90t, ndb, d25b, d50b, d90b,
                       nst, nstp, sa50, samx,
                       par_del_tot > 0 ? par_del_id / par_del_tot : 0);
                if (tri_on)
                    for (int t2 = 0; t2 < ntri; t2++)
                        printf("# PAR tri T%d gid=(%lld,%lld,%lld) age=(%.0f,%.0f,%.0f)\n",
                               t2, cgid_[tri_v[t2][0]], cgid_[tri_v[t2][1]],
                               cgid_[tri_v[t2][2]],
                               cgid_[tri_v[t2][0]] ? sim_t - cbirth_[tri_v[t2][0]] : 0,
                               cgid_[tri_v[t2][1]] ? sim_t - cbirth_[tri_v[t2][1]] : 0,
                               cgid_[tri_v[t2][2]] ? sim_t - cbirth_[tri_v[t2][2]] : 0);
            }
            if (P.amp_tau > 0) {
                double a25, a50, a75, a90, amg, b25, b50, b75, b90, bmg;
                int ntp = amp_stats(0, &a25, &a50, &a75, &a90, &amg);
                int nba = amp_stats(1, &b25, &b50, &b75, &b90, &bmg);
                printf("# AMP t=%.2f tp n=%d rho=[%.3f %.3f %.3f %.3f] mag=%.5f | ba n=%d rho=[%.3f %.3f %.3f %.3f] mag=%.5f\n",
                       sim_t, ntp, a25, a50, a75, a90, amg,
                       nba, b25, b50, b75, b90, bmg);
            }
            if (P.bh_r > 0)
                printf("# HOLE t=%.2f Eh=%.6f nin=%ld eatF=%.4f eatM=%.4f eatS=%.4f\n",
                       sim_t, Eh_total, bh_nin, bh_eat_f, bh_eat_m, bh_eat_s);
            if (P.bed_k > 0) {
                double b25, b50, b75, bmx; int bng, bnl;
                bed_stats(&b25, &b50, &b75, &bmx, &bng, &bnl);
                printf("# BED t=%.2f n=%d q=[%.3f %.3f %.3f] max=%.3f grown=%d\n",
                       sim_t, bnl, b25, b50, b75, bmx, bng);
            }
        }
        if (P.snap_every > 0 && P.snap_dir[0] && st % P.snap_every == 0)
            write_fcs(st / P.snap_every);
        if (fcs_stream && P.snap_every > 0 && st % P.snap_every == 0) {
            fcs_cell_frame(fcs_stream);
            fcs_instrument(fcs_stream);
        }
        if (P.bed_k > 0 && P.snap_every > 0 && st % P.snap_every == 0)
            bed_map_dump();
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
        if (P.sect_meter) {
            double tn = (st + 1) * P.dt;
            if (tn >= P.sect_t0 && tn <= P.sect_t1) {
                int NS = P.sect_n;
                for (int i = 0; i < NC; i++) {
                    double dx = wr(px_[i] - sect_cx), dy = wr(py_[i] - sect_cy);
                    double r2 = dx*dx + dy*dy;
                    if (r2 < P.sect_r0 * P.sect_r0 || r2 >= P.sect_r1 * P.sect_r1) continue;
                    /* half-bin rotation: sector k is CENTRED at k*2pi/NS */
                    double u = (atan2(dy, dx) + M_PI / NS) / TWO_PI;
                    u -= floor(u);
                    int k = (int)(u * NS);
                    if (k >= NS) k = NS - 1;
                    sectE[k] += Ee[i] * P.dt;
                    sectN[k] += P.dt;
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
        printf("# RESULT conv rough=%.6f cond=%.6f evap=%.6f backs=%.6f rad=%.6f\n",
               rough_total, cond_total, evap_total, backs_total, rad_total);
        {   /* pad-31 apparatus: the standing sub-atom credit inventory */
            double cd = 0, cf = 0;
            for (int ii = 0; ii < NC; ii++) { cd += qcnvD[ii]; cf += qcnvF[ii]; }
            printf("# RESULT credit qcnvD=%.6f qcnvF=%.6f\n", cd, cf);
        }
        if (P.k_cant > 0 || P.k_tune > 0) {
            double at = 0, am = 0, xt = 0; int ntg = 0, nl = 0;
            for (int ii = 0; ii < NC; ii++) {
                double av = cant_amp_of(ii);
                if (av > am) am = av;
                if (av > 0.5) nl++;
                if (tag[ii]) { at += av; xt += cxl_[ii]; ntg++; }
            }
            printf("# RESULT cantus a_tag=%.4f a_max=%.4f xl_tag=%.4f nlock=%d tune_total=%.6f\n",
                   ntg ? at / ntg : 0, am, ntg ? xt / ntg : 0, nl, tune_total);
        }
        if (P.reg_tau > 0) {
            double a25, a50, a75, a90, agr, afl, b25, b50, b75, b90, bgr, bfl;
            int ntp = reg_stats(0, &a25, &a50, &a75, &a90, &agr, &afl);
            int nba = reg_stats(1, &b25, &b50, &b75, &b90, &bgr, &bfl);
            printf("# RESULT reg tp n=%d rho=[%.4f %.4f %.4f %.4f] gr=%.6f fl=%.4f | ba n=%d rho=[%.4f %.4f %.4f %.4f] gr=%.6f fl=%.4f\n",
                   ntp, a25, a50, a75, a90, agr, afl,
                   nba, b25, b50, b75, b90, bgr, bfl);
        }
        if (P.par_tau > 0) {
            int ngid = 0;
            for (int i2 = 0; i2 < NC; i2++) if (cgid_[i2] != 0) ngid++;
            printf("# RESULT par mints=%lld retires=%lld live=%d del_id=%.6f del_tot=%.6f idfrac=%.4f\n",
                   par_mints, par_retires, ngid, par_del_id, par_del_tot,
                   par_del_tot > 0 ? par_del_id / par_del_tot : 0);
        }
        if (P.bh_r > 0)
            printf("# RESULT hole Eh=%.6f nin=%ld eatF=%.6f eatM=%.6f eatS=%.6f\n",
                   Eh_total, bh_nin, bh_eat_f, bh_eat_m, bh_eat_s);
        if (P.bed_k > 0) {
            double b25, b50, b75, bmx; int bng, bnl;
            bed_stats(&b25, &b50, &b75, &bmx, &bng, &bnl);
            printf("# RESULT bed n=%d q=[%.4f %.4f %.4f] max=%.4f grown=%d\n",
                   bnl, b25, b50, b75, bmx, bng);
        }
        if (P.amp_tau > 0) {
            double a25, a50, a75, a90, amg, b25, b50, b75, b90, bmg;
            int ntp = amp_stats(0, &a25, &a50, &a75, &a90, &amg);
            int nba = amp_stats(1, &b25, &b50, &b75, &b90, &bmg);
            printf("# RESULT amp tp n=%d rho=[%.4f %.4f %.4f %.4f] mag=%.6f | ba n=%d rho=[%.4f %.4f %.4f %.4f] mag=%.6f\n",
                   ntp, a25, a50, a75, a90, amg,
                   nba, b25, b50, b75, b90, bmg);
        }
        if (P.slit_obj || P.convtag)
            /* net field capture at the occulter = cond - evap - rough + backs
             * (evap returns d2-bs to field; rough returns rough-back_s) */
            printf("# RESULT convtag rough=%.6f cond=%.6f evap=%.6f backs=%.6f net=%.6f\n",
                   ct_rough, ct_cond, ct_evap, ct_backs,
                   ct_cond - ct_evap - ct_rough + ct_backs);
        if (P.p1_meter) {
            /* cumulative first-moment flow (energy*length) per channel,
             * plus the rate: tot/T = the all-modes center-of-energy
             * velocity times total energy (the honest Delta-p number) */
            double tx = p1sp[0] + p1fl[0] + p1fd[0] + p1gm[0];
            double ty = p1sp[1] + p1fl[1] + p1fd[1] + p1gm[1];
            double tz = p1sp[2] + p1fl[2] + p1fd[2] + p1gm[2];
            printf("# RESULT p1 tot=(%+.6e,%+.6e,%+.6e) rate=(%+.6e,%+.6e,%+.6e)\n",
                   tx, ty, tz, tx / P.T, ty / P.T, tz / P.T);
            printf("# RESULT p1x sp=%+.6e fl=%+.6e fd=%+.6e gm=%+.6e\n",
                   p1sp[0], p1fl[0], p1fd[0], p1gm[0]);
        }
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
                /* directed ledger (PAULI-0 rate control): cumulative
                 * deposits each way + the standing (never-landed) flight */
                printf("# RESULT pairflux dep_ij=%.6e dep_ji=%.6e lem_ij=%.6e lem_ji=%.6e\n",
                       sfluxd[2*s], sfluxd[2*s+1], slem[2*s], slem[2*s+1]);
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
        if (rint_n > 0) {
            for (int e = 0; e < rint_n; e++) {
                int i = rint_e[e][0], j = rint_e[e][1];
                int s = hfind(i < j ? i : j, i < j ? j : i);
                if (s < 0 || sst[s] == S_FREE) {
                    printf("# RESULT rint%d %d-%d CHANNEL DEAD\n", e, i, j);
                    continue;
                }
                double p = rint_p[e], q = rint_q[e];
                double psi = wrap_pi(q * th2[i] - q * w2e[i] * sd[s] / P.C - p * th2[j]);
                double gf = gate_of(q * th2[i] - q * w2e[i] * sd[s] / P.C - p * th2[j]);
                double gb = gate_of(p * th2[j] - p * w2e[j] * sd[s] / P.C - q * th2[i]);
                printf("# RESULT rint%d %d-%d p:q=%d:%d d=%.6f dev=%+.6f psi=%+.4f gg=%.4f lem=%.6e\n",
                       e, i, j, rint_p[e], rint_q[e], sd[s], sd[s] - rint_dstar[e],
                       psi, gf * gb, slem[2*s] + slem[2*s+1]);
            }
            for (int k = 0; k < rings_nring; k++) {
                double ggs = 0, xm = 0, ems = 0; int ne = 0;
                for (int e = 0; e < truss_n; e++) {
                    int i = truss_e[e][0];
                    if (i < rings_v0[k] || i >= rings_v0[k] + P.rings_nv) continue;
                    int j = truss_e[e][1];
                    int s = hfind(i < j ? i : j, i < j ? j : i);
                    if (s >= 0 && sst[s] != S_FREE) {
                        double gf = gate_of(slq[s]*th2[i] - slq[s]*w2e[i]*sd[s]/P.C - slp[s]*th2[j]);
                        double gb = gate_of(slp[s]*th2[j] - slp[s]*w2e[j]*sd[s]/P.C - slq[s]*th2[i]);
                        ggs += gf * gb; ne++;
                    }
                }
                for (int q = 0; q < P.rings_nv; q++) {
                    int i = rings_v0[k] + q;
                    xm += (Em[i] + flload[i]) / P.cap;
                    ems += Em[i];
                }
                printf("# RESULT ringq%d mode=%s w=%.4f gg_int=%.4f xbar=%.4f Em=%.5f\n",
                       k, rings_mode[k] ? "MOLECULE" : "DROPLET", rings_w[k],
                       ne ? ggs / ne : 0, xm / P.rings_nv, ems);
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
        if (!strcmp(P.exp, "blob2") && nfsamp >= 6) {
            int h = nfsamp / 2;
            /* approach speeds from the FIRST half (before contact) */
            double sx=0,sxx=0,sya=0,syb=0,syt=0,sxa=0,sxb=0,sxt=0;
            for (int k2 = 0; k2 < h; k2++) {
                sx += fs_t[k2]; sxx += fs_t[k2]*fs_t[k2];
                sya += b2_ax[k2]; sxa += fs_t[k2]*b2_ax[k2];
                syb += b2_bx[k2]; sxb += fs_t[k2]*b2_bx[k2];
                syt += b2_tx[k2]; sxt += fs_t[k2]*b2_tx[k2];
            }
            double den = h*sxx - sx*sx;
            double vA = den != 0 ? (h*sxa - sx*sya)/den : 0;
            double vB = den != 0 ? (h*sxb - sx*syb)/den : 0;
            double vT = den != 0 ? (h*sxt - sx*syt)/den : 0;
            printf("# RESULT blob2 vA=%.6f vB=%.6f closing=%.6f vTotalCOM=%.2e sep_final=%.4f\n",
                   vA, vB, vA - vB, vT, b2_bx[nfsamp-1] - b2_ax[nfsamp-1]);
        }
        if (P.slit_clicks) {
            double esum = 0;
            for (int k2 = 0; k2 < nclick; k2++) esum += click_e[k2];
            printf("# RESULT clicks n=%d e_sum=%.4f\n", nclick, esum);
            for (int k2 = 0; k2 < nclick; k2++)
                printf("# CLICK t=%.2f y=%.4f e=%.5f\n",
                       click_t[k2], click_y[k2], click_e[k2]);
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
                "RESULT conv rough=%.6f cond=%.6f evap=%.6f backs=%.6f rad=%.6f\n",
                (Etot - E0_total) / (E0_total != 0 ? E0_total : 1),
                births, deaths, beta_returns,
                rough_total, cond_total, evap_total, backs_total, rad_total);
            fcs_anlz_text(fcs_stream, fin);
            fclose(fcs_stream);
            fcs_stream = NULL;
        }
        if (P.sect_meter) {
            double stot = 0;
            for (int k = 0; k < P.sect_n; k++) stot += sectE[k];
            printf("# RESULT sect Etot=%.6f n=%d r0=%g r1=%g centre=(%.2f,%.2f) gate=[%g,%g]\n",
                   stot, P.sect_n, P.sect_r0, P.sect_r1, sect_cx, sect_cy,
                   P.sect_t0, P.sect_t1);
            for (int k = 0; k < P.sect_n; k++) {
                double thc = k * TWO_PI / P.sect_n;
                if (thc > M_PI) thc -= TWO_PI;
                printf("# SECT k=%d th=%+.4f E=%.8f n=%.2f\n", k, thc, sectE[k], sectN[k]);
            }
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
