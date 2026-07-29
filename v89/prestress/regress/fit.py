#!/usr/bin/env python3
# v89 LEAK-LAW REGRESSION — sparse symbolic discovery over the existing run corpus.
#
# Re-runnable end to end:
#   python3 v89/prestress/regress/fit.py
# Reads:  the MASS/H-campaign logs in LOGDIR (scratchpad corpus)
# Writes: corpus.tsv (one row per run) next to this file, and prints every
#         table quoted in REGRESS.md / PREDICTIONS.md.
#
# Doctrine: v89 material only. No simulator runs. numpy/scipy only
# (PySR/Julia unavailable); STLSQ + exhaustive best-subset are hand-rolled.
#
# ---------------------------------------------------------------------------
# DEFINITIONS (documented here once; REGRESS.md quotes them)
#
# M_sum(t)    sum of the census lump masses listed on the '# LUMP' row
#             (top-10 shown; n<=10 for every run here except transient
#             blob debris, so the truncation is negligible). This is "the
#             object": census arcs of one ring are threshold artifacts
#             (MASS.md census caveat), so the sum, not the top lump, is
#             the object's mass.
# M_top(t)    largest lump — used ONLY to reproduce MASS.md's quoted rates
#             (the campaign log used the top lump) as a parser sanity gate.
# t_death     last t at which the 21-sample (+/-10 t.u.) rolling MEDIAN of
#             M_sum >= M_DEAD (0.15, just above the census threshold 0.1:
#             one sub-threshold crumb = dead). The median kills the
#             bath-flicker at the death tail (single samples popping over
#             threshold). Censored if that last t is within 30 t.u. of the
#             horizon T. Cross-checks kept: t_last_raw (last raw sample
#             >= M_DEAD) and t_n0 (first sustained census n==0).
#             NOTE vs MASS.md: the campaign log's death quotes are
#             approximate ("dissolved BY t~1900" — the census shows
#             sustained n=0 from t=1672); this rule is the reproducible
#             sharp version, and it reproduces the v3/a1c/h1-v3 quotes
#             to within ~1-2%.
# leak_early  100*(M_sum(300)-M_sum(50))/(250*M_sum(50)) [%/t.u.] over
#             [50,300] — the M1-comparable window every run has (for runs
#             shorter-lived than 300, window ends at min(300, t_death)).
# dMdt_late   linear-fit slope of M_sum over the last 20% of life
#             ([0.8*t_end, t_end], t_end = t_death or T). Units/t.u.
# rough_share (rough(t_end)-rough(50)) / (M_sum(50)-M_sum(t_end)); CONV
#             counters are cumulative. Share of the object's loss that
#             left as roughness radiation (D->F off-rung deliveries).
# x_voice(t)  M_sum(t)/(N*cap): mean per-voice load. N = seeded voice
#             count (rings: ring n; shells: 8 cube cells). Blobs have no
#             N -> no x_voice. cap = 2.5 (laws_V2g).
# t_skirt     first t with 0 < x_voice(t) <= X_SKIRT (0.0617, e7 boundary).
# x50         x_voice at t=50 (post-transient settled load; the seed's
#             t=0 surplus radiates in the first ~10-50 t.u.).
#
# PREDICTED SEED GATES (B1 algebra; laws_V2g p_gate=8):
#   g(psi) = ((1+cos psi)/2)^8
#   phi_bar = 2*pi*closure/N  (closure = printed closure/2pi; mean per-link
#             retardation omega*d/C actually seeded)
#   delta   = closure - round(closure)  (seam defect, dumped on the seam
#             link by the lock recursion; ring_m seeds have delta = 0)
#   regular links (lock recursion, wind=0): fwd psi=0 -> g=1
#                                           back psi = -2*phi_bar
#   seam link:  fwd psi = -2*pi*delta ; back psi = -(2*phi_bar + 2*pi*delta)
#   naive wound (ring_wind=w): every link fwd psi = -2*pi*w/N,
#       back psi = -(2*phi_bar - 2*pi*w/N); seam gets the extra -2*pi*delta.
#   ring_m seeds: uniform omega = 2*pi*m*C/Lring -> phi = 2*pi*m/N exactly,
#       fwd g=1 on all N links, back g(2*phi), no seam defect.
#   shells: seed line prints measured gates min/mean directly (cube: 12
#       edges, mutual pi-rung -> back=fwd). Directed count = 24.
#   Gate stats over the 2N directed seed links:
#       gmin, gmean, Dgate = sum(1-g), Bback = sum of back-direction gates.
# ---------------------------------------------------------------------------

import os, re, sys, math, itertools
import numpy as np

HERE   = os.path.dirname(os.path.abspath(__file__))
LOGDIR = "/tmp/claude-1000/-home-d-code-scp/63080320-b2b9-401b-b061-05d977d0dcbf/scratchpad"
CAP      = 2.5
X_SKIRT  = 0.0617
M_DEAD   = 0.15
P_GATE   = 8

def g_of(psi):
    return ((1.0 + math.cos(psi)) / 2.0) ** P_GATE

def wrap(psi):
    while psi >  math.pi: psi -= 2*math.pi
    while psi < -math.pi: psi += 2*math.pi
    return psi

# ---------------------------------------------------------------- parsing --
LUMP_RE  = re.compile(r'^# LUMP t=([\d.]+) n=(\d+) Emfree=([\d.eE+-]+)(.*)$')
MASS_RE  = re.compile(r'm=([\d.eE+-]+)')
CONV_RE  = re.compile(r'^# CONV t=([\d.]+) cond=([\d.eE+-]+) evap=([\d.eE+-]+) rough=([\d.eE+-]+) back_s=([\d.eE+-]+)')
RING_RE  = re.compile(r'^# ring: seeded n=(\d+) R=([\d.]+) d_target=([\d.]+)(?: Lring=([\d.]+))? closure/2pi=([\d.]+) wind=(\d+)(?: ring_m=(-?\d+))?')
SHELL_RE = re.compile(r'^# shell: cube a_target=([\d.]+) abar=([\d.]+) omega=([\d.]+) x=([\d.]+) gates min=([\d.]+) mean=([\d.]+)')
CORE_RE  = re.compile(r'^# shell core: es_core=([\d.]+) r_core=([\d.]+) cells=(\d+)')

def parse_log(path):
    d = dict(cfg={}, ring=None, shell=None, core=None, edge_sink=False,
             lump_t=[], lump_n=[], lump_sum=[], lump_top=[], lump_emfree=[],
             conv_t=[], conv=[])
    with open(path) as f:
        for line in f:
            if line.startswith('# cfg '):
                for k, v in re.findall(r'([A-Za-z_0-9]+)=([^\s]+)', line):
                    d['cfg'].setdefault(k, v)
                continue
            if line.startswith('# edge_sink:'):
                d['edge_sink'] = True; continue
            m = RING_RE.match(line)
            if m:
                d['ring'] = dict(n=int(m.group(1)), R=float(m.group(2)),
                                 d_target=float(m.group(3)),
                                 Lring=float(m.group(4)) if m.group(4) else None,
                                 closure=float(m.group(5)), wind=int(m.group(6)),
                                 ring_m=int(m.group(7)) if m.group(7) else None)
                continue
            m = SHELL_RE.match(line)
            if m:
                d['shell'] = dict(a_target=float(m.group(1)), abar=float(m.group(2)),
                                  omega=float(m.group(3)), x=float(m.group(4)),
                                  gmin=float(m.group(5)), gmean=float(m.group(6)))
                continue
            m = CORE_RE.match(line)
            if m:
                d['core'] = dict(es_core=float(m.group(1)), r_core=float(m.group(2)),
                                 cells=int(m.group(3)))
                continue
            m = LUMP_RE.match(line)
            if m:
                t = float(m.group(1)); n = int(m.group(2)); em = float(m.group(3))
                masses = [float(x) for x in MASS_RE.findall(m.group(4))]
                d['lump_t'].append(t); d['lump_n'].append(n)
                d['lump_sum'].append(sum(masses))
                d['lump_top'].append(masses[0] if masses else 0.0)
                d['lump_emfree'].append(em)
                continue
            m = CONV_RE.match(line)
            if m:
                d['conv_t'].append(float(m.group(1)))
                d['conv'].append(tuple(float(m.group(i)) for i in (2, 3, 4, 5)))
    for k in ('lump_t', 'lump_n', 'lump_sum', 'lump_top', 'lump_emfree', 'conv_t'):
        d[k] = np.array(d[k])
    d['conv'] = np.array(d['conv']) if len(d['conv']) else np.zeros((0, 4))
    return d

# ------------------------------------------------------------- gate algebra --
def ring_gates(ring):
    """Predicted directed seed gates from the B1 algebra. Returns dict."""
    N   = ring['n']
    clo = ring['closure']
    w   = ring['wind']
    phib = 2*math.pi*clo/N            # mean per-link retardation
    delta = clo - round(clo)          # seam defect (fraction of 2pi)
    exact = ring['ring_m'] is not None and ring['ring_m'] > 0
    fwd, back = [], []
    if exact:
        m = ring['ring_m']
        phi = 2*math.pi*m/N
        for _ in range(N):
            fwd.append(g_of(0.0))
            back.append(g_of(wrap(2*phi)))
    elif w == 0:
        for _ in range(N-1):
            fwd.append(g_of(0.0))
            back.append(g_of(wrap(2*phib)))
        fwd.append(g_of(wrap(2*math.pi*delta)))            # seam fwd
        back.append(g_of(wrap(2*phib + 2*math.pi*delta)))  # seam back
    else:  # naive winding: phase kick 2*pi*w/N on every link
        kick = 2*math.pi*w/N
        for _ in range(N-1):
            fwd.append(g_of(wrap(kick)))
            back.append(g_of(wrap(2*phib - kick)))
        fwd.append(g_of(wrap(kick + 2*math.pi*delta)))
        back.append(g_of(wrap(2*phib - kick + 2*math.pi*delta)))
    fwd, back = np.array(fwd), np.array(back)
    allg = np.concatenate([fwd, back])
    mm = round(clo) if not exact else ring['ring_m']
    return dict(gfwd_reg=fwd[0], gback_reg=back[0],
                gseam_fwd=fwd[-1], gseam_back=back[-1],
                gmin=allg.min(), gmean=allg.mean(),
                Dgate=float(np.sum(1.0 - allg)), Bback=float(back.sum()),
                m=mm, defect=abs(delta), links=2*N)

def shell_gates(shell):
    gmean, gmin = shell['gmean'], shell['gmin']
    links = 24                       # cube: 12 edges, mutual (back=fwd)
    return dict(gfwd_reg=gmean, gback_reg=gmean,
                gseam_fwd=float('nan'), gseam_back=float('nan'),
                gmin=gmin, gmean=gmean,
                Dgate=links*(1.0-gmean), Bback=12*gmean,
                m=float('nan'), defect=float('nan'), links=links)

# ------------------------------------------------------------ run registry --
# (family, class, N_voices, use_fit, note). Duplicates verified by LUMP-series
# diff (see REGRESS.md "corpus hygiene"): shorter runs are byte-identical
# physics prefixes of the longer run with the same seed (deterministic kernel).
RUNS = [
    # name                class        use   dup_of          note
    ("m0_blob",          "blob",       True,  None,  "M0 heavy blob, leak taxonomy run"),
    ("viz_blob",         "blob",       False, "m0_blob", "T=120 render duplicate"),
    ("m1_ctrl",          "blob",       True,  None,  "M1 light gaussian control (voice-scale)"),
    ("h1_solid",         "blob",       True,  None,  "H1 solid ball, shell-mass-matched"),
    ("m1_ring6",         "ring_auto",  False, "a1_ring6", "T=300 prefix"),
    ("m2_ring6",         "ring_auto",  False, "a1_ring6", "T=1000 prefix"),
    ("viz_ring6",        "ring_auto",  False, "a1_ring6", "T=120 render duplicate"),
    ("a1_ring6",         "ring_auto",  True,  None,  "ring6 m=3, open box, T=5000"),
    ("a1c_ring6_closed", "ring_auto",  True,  None,  "ring6 m=3, CLOSED box (no edge_sink)"),
    ("a2_s1",            "ring_auto",  True,  None,  "ring6 foam seed 20260721"),
    ("a2_s2",            "ring_auto",  True,  None,  "ring6 foam seed 20260722 (fragmented)"),
    ("a2_s3",            "ring_auto",  True,  None,  "ring6 foam seed 20260723"),
    ("a2_s4",            "ring_auto",  True,  None,  "ring6 foam seed 20260724"),
    ("a2_s5",            "ring_auto",  True,  None,  "ring6 foam seed 20260725 (died)"),
    ("m1_ring12",        "ring_auto",  False, "m2_ring12", "T=300 prefix"),
    ("m2_ring12",        "ring_auto",  True,  None,  "ring12 auto (closure 5.917: defect 0.083)"),
    ("m1_ring6w",        "ring_naive", True,  None,  "ring6 naive wind=1 (all gates ~0.100)"),
    ("b1_naive12",       "ring_naive", True,  None,  "ring12 naive wind=1"),
    ("b1_comp12",        "chain",      True,  None,  "fleet-v2 comp12: seam dead -> open chain"),
    ("b1_comp6",         "chain",      True,  None,  "fleet-v2 comp6: seam dead -> open chain"),
    ("b1_unwound_match", "chain",      True,  None,  "fleet-v2 unwound: seam dead -> open chain"),
    ("v3_comp12",        "ring_exact", False, "v3_comp12_L", "T=1000 prefix"),
    ("v3_comp6",         "ring_exact", False, "v3_comp6_L",  "T=1000 prefix"),
    ("v3_unwound12",     "ring_exact", False, "v3_unwound12_L", "T=1000 prefix"),
    ("v3_comp12_L",      "ring_exact", True,  None,  "exact m=5 wound mutual (longest-lived)"),
    ("v3_comp6_L",       "ring_exact", True,  None,  "exact m=2 one-way (back gate 1.5e-5)"),
    ("v3_unwound12_L",   "ring_exact", True,  None,  "exact m=6 unwound"),
    ("h1_shell0",        "shell",      True,  None,  "cube v1 small (empty-hole accident twin)"),
    ("h1_bubble",        "shell",      False, "h1_shell0", "LUMP-identical to shell0 (hole empty)"),
    ("h1_shell0_v2",     "shell",      True,  None,  "cube v2 tight/light (skirt-class)"),
    ("h1_shell0_v3",     "shell",      True,  None,  "cube v3 heavy, matched to bubble_v2"),
    ("h1_bubble_v2",     "shell",      True,  None,  "cube heavy + 2-cell pocket"),
    ("h1_pocket",        "pocket",     False, None,  "naked Es pocket, no dense object"),
    ("h1_pocket_v2",     "pocket",     False, None,  "naked Es pocket, no dense object"),
]

def analyze(name, cls):
    p = os.path.join(LOGDIR, name + ".log")
    d = parse_log(p)
    t, Ms, Mt = d['lump_t'], d['lump_sum'], d['lump_top']
    r = dict(run=name, cls=cls)
    r['seed'] = int(d['cfg'].get('seed', 0)); r['L'] = float(d['cfg'].get('L', 24))
    r['T'] = float(d['cfg'].get('T', t[-1] if len(t) else 0))
    r['closed_box'] = 0 if d['edge_sink'] else 1

    # ---- descriptors
    if d['ring']:
        rg = d['ring']; gates = ring_gates(rg)
        r.update(N=rg['n'], d_ab=rg['d_target'], closure=rg['closure'],
                 wind=rg['wind'], **gates)
    elif d['shell'] and d['lump_emfree'][0] > 0:
        sh = d['shell']; gates = shell_gates(sh)
        r.update(N=8, d_ab=sh['abar'], closure=float('nan'), wind=0, **gates)
    else:
        r.update(N=float('nan'), d_ab=float('nan'), closure=float('nan'), wind=0,
                 gfwd_reg=float('nan'), gback_reg=float('nan'),
                 gseam_fwd=float('nan'), gseam_back=float('nan'),
                 gmin=float('nan'), gmean=float('nan'), Dgate=float('nan'),
                 Bback=float('nan'), m=float('nan'), defect=float('nan'), links=0)
    r['es_core'] = d['core']['es_core'] if (d['core'] and d['core']['cells'] > 0) else 0.0

    if len(t) == 0:
        return r
    # ---- responses
    r['M0']  = Ms[0]
    i50  = int(np.searchsorted(t, 50.0));  i50 = min(i50, len(t)-1)
    r['M50'] = Ms[i50]
    Nv = r['N']
    r['x0']  = Ms[0]  / (Nv*CAP) if Nv == Nv else float('nan')
    r['x50'] = Ms[i50] / (Nv*CAP) if Nv == Nv else float('nan')

    # death time: last t where the rolling median (21 samples) of M_sum >= M_DEAD
    if len(Ms) >= 21:
        pad = 10
        Mmed = np.array([np.median(Ms[max(0, i-pad):min(len(Ms), i+pad+1)])
                         for i in range(len(Ms))])
    else:
        Mmed = Ms.copy()
    alive_raw = np.where(Ms >= M_DEAD)[0]
    r['t_last_raw'] = t[alive_raw[-1]] if len(alive_raw) else 0.0
    alive = np.where(Mmed >= M_DEAD)[0]
    if len(alive) == 0:
        r['t_death'] = 0.0; r['censored'] = 0
    else:
        tl = t[alive[-1]]
        cens = (t[-1] - tl) <= 30.0
        r['t_death'] = t[-1] if cens else tl
        r['censored'] = 1 if cens else 0
    # n==0 sustained cross-check
    z = np.where(d['lump_n'] == 0)[0]
    tn0 = float('nan')
    for i in z:
        j = np.searchsorted(t, t[i] + 50.0)
        if np.all(d['lump_n'][i:min(j, len(t))] == 0):
            tn0 = t[i]; break
    r['t_n0'] = tn0
    # structured-phase end: last t with rolling-median M_sum >= 1.0 (one
    # voice-mass). Robustness alternative to t_death: parked sub-voice
    # remnant crumbs (common on favorable foam cells) extend t_death but
    # not t_struct.
    aliv1 = np.where(Mmed >= 1.0)[0]
    if len(aliv1) == 0:
        r['t_struct'] = 0.0; r['cens_struct'] = 0
    else:
        tl1 = t[aliv1[-1]]
        c1 = (t[-1] - tl1) <= 30.0
        r['t_struct'] = t[-1] if c1 else tl1
        r['cens_struct'] = 1 if c1 else 0

    t_end = r['t_death']
    # early leak window [50, min(300, t_end)]
    e1 = min(300.0, t_end)
    j1 = int(np.searchsorted(t, e1)); j1 = min(j1, len(t)-1)
    if t[j1] > t[i50] and Ms[i50] > 0:
        r['leak_early'] = 100.0*(Ms[j1]-Ms[i50])/((t[j1]-t[i50])*Ms[i50])
        r['leak_early_top'] = 100.0*(Mt[j1]-Mt[i50])/((t[j1]-t[i50])*Mt[i50]) if Mt[i50] > 0 else float('nan')
    else:
        r['leak_early'] = float('nan'); r['leak_early_top'] = float('nan')
    # late window: last 20% of life
    w0 = 0.8*t_end
    sel = (t >= w0) & (t <= t_end)
    if sel.sum() >= 5:
        A = np.vstack([t[sel], np.ones(sel.sum())]).T
        slope = np.linalg.lstsq(A, Ms[sel], rcond=None)[0][0]
        r['dMdt_late'] = slope
    else:
        r['dMdt_late'] = float('nan')
    # mid-life exponential rate (for the skirt-integral death law):
    # fit ln M_sum on [50, 0.8*t_end] where M_sum > 0
    sel = (t >= 50.0) & (t <= 0.8*t_end) & (Ms > 0.05)
    if sel.sum() >= 10:
        A = np.vstack([t[sel], np.ones(sel.sum())]).T
        r['k_exp'] = -np.linalg.lstsq(A, np.log(Ms[sel]), rcond=None)[0][0]
    else:
        r['k_exp'] = float('nan')
    # roughness share over [50, t_end]
    if len(d['conv_t']):
        ct = d['conv_t']; rough = d['conv'][:, 2]
        c0 = int(np.searchsorted(ct, 50.0)); c1 = int(np.searchsorted(ct, t_end)); c1 = min(c1, len(ct)-1)
        dM = Ms[i50] - Ms[min(int(np.searchsorted(t, t_end)), len(t)-1)]
        r['rough_tot'] = rough[-1]
        r['rough_share'] = (rough[c1]-rough[c0])/dM if dM > 1e-9 else float('nan')
    else:
        r['rough_tot'] = float('nan'); r['rough_share'] = float('nan')
    # time to skirt
    if Nv == Nv:
        xv = Ms/(Nv*CAP)
        below = np.where((xv > 0) & (xv <= X_SKIRT))[0]
        r['t_skirt'] = t[below[0]] if len(below) else float('nan')
    else:
        r['t_skirt'] = float('nan')
    # effective per-voice leak current over the structured life:
    # c_eff = cap*(x50 - x_skirt)/(t_death - 50). For censored runs this is
    # an UPPER BOUND. NaN when the run starts at the skirt (margin < 0.02)
    # or has no voice count.
    marg = (r['x50'] - X_SKIRT) if Nv == Nv else float('nan')
    if marg == marg and marg > 0.02 and r['t_death'] > 100:
        r['c_eff'] = CAP*marg/(r['t_death'] - 50.0)
    else:
        r['c_eff'] = float('nan')
    return r

# ------------------------------------------------------------------- build --
rows = []
for name, cls, use, dup, note in RUNS:
    r = analyze(name, cls)
    r['use_fit'] = 1 if use else 0
    r['dup_of'] = dup or "-"
    r['note'] = note
    rows.append(r)

COLS = ["run","cls","use_fit","dup_of","seed","L","T","closed_box","N","d_ab",
        "closure","m","defect","wind","es_core","gfwd_reg","gback_reg",
        "gseam_fwd","gseam_back","gmin","gmean","Dgate","Bback","links",
        "M0","M50","x0","x50","t_death","censored","t_struct","cens_struct",
        "t_last_raw","t_n0","leak_early","leak_early_top","dMdt_late",
        "k_exp","c_eff","rough_tot","rough_share","t_skirt","note"]

def fmt(v):
    if isinstance(v, float):
        if v != v: return "nan"
        if v == 0: return "0"
        if abs(v) >= 1e4 or abs(v) < 1e-3: return f"{v:.4g}"
        return f"{v:.4f}"
    return str(v)

with open(os.path.join(HERE, "corpus.tsv"), "w") as f:
    f.write("\t".join(COLS) + "\n")
    for r in rows:
        f.write("\t".join(fmt(r.get(c, "nan")) for c in COLS) + "\n")
print(f"corpus.tsv written: {len(rows)} rows "
      f"({sum(r['use_fit'] for r in rows)} unique fit rows)")

# ------------------------------------------------- parser sanity vs MASS.md --
print("\n=== SANITY GATE — reproduce MASS.md quoted numbers ===")
BY = {r['run']: r for r in rows}
checks = [
    ("a1_ring6 death (MASS 'by ~1900'; census n=0 from 1672)",
                                    BY['a1_ring6']['t_death'],        1672),
    ("a1c closed death ~1600",      BY['a1c_ring6_closed']['t_death'],1600),
    ("v3_unwound12_L death ~2221",  BY['v3_unwound12_L']['t_death'],  2221),
    ("v3_comp6_L death ~3836",      BY['v3_comp6_L']['t_death'],      3836),
    ("v3_comp12_L censored @5000",  BY['v3_comp12_L']['t_death'],     5000),
    ("h1_shell0 death ~1222",       BY['h1_shell0']['t_death'],       1222),
    ("h1_solid death ~1158 (MASS; census lump to ~1270)",
                                    BY['h1_solid']['t_death'],        1270),
    ("h1_shell0_v2 death ~237 (MASS; flicker tail to ~310)",
                                    BY['h1_shell0_v2']['t_death'],    280),
    ("h1_shell0_v3 death ~1814",    BY['h1_shell0_v3']['t_death'],    1814),
    ("h1_bubble_v2 death ~1749",    BY['h1_bubble_v2']['t_death'],    1749),
    ("a1_ring6 leak_top ~-0.058",   BY['a1_ring6']['leak_early_top'], -0.058),
    ("m2_ring12 leak_top ~-0.051",  BY['m2_ring12']['leak_early_top'],-0.051),
]
for label, got, want in checks:
    ok = abs(got - want) <= max(0.12*abs(want), 0.012)
    print(f"  [{'ok' if ok else 'XX'}] {label}: got {got:.4g} (ref {want})")
# basis notes (MASS.md mixed measurement bases; reproduce each explicitly):
d_ctrl = parse_log(os.path.join(LOGDIR, "m1_ctrl.log"))
te, Em = d_ctrl['lump_t'], d_ctrl['lump_emfree']
i5, i3 = int(np.searchsorted(te, 50.0)), int(np.searchsorted(te, 300.0)); i3 = min(i3, len(te)-1)
ctrl_emrate = 100.0*(Em[i3]-Em[i5])/((te[i3]-te[i5])*Em[i5])
print(f"  [note] m1_ctrl MASS -0.232 is the Emfree basis: got {ctrl_emrate:.4g} "
      f"(census-sum basis gives {BY['m1_ctrl']['leak_early']:.4g})")
print(f"  [FLAG] m1_ring6w MASS -0.119 NOT reproducible from the census series "
      f"(top {BY['m1_ring6w']['leak_early_top']:.4g}, sum {BY['m1_ring6w']['leak_early']:.4g}); "
      f"the 'naive winding 2x leak' claim rests on an unrecoverable basis — "
      f"treated as unverified in the regression.")
for cname in ("b1_comp12", "b1_comp6", "b1_unwound_match"):
    print(f"  [note] {cname} (chain) leak sum-basis [50,300]: {BY[cname]['leak_early']:.4g} "
          f"(MASS quoted chains -0.092..-0.097)")

# ------------------------------------------------------------ fit matrices --
# Structured corpus (rings + chains + shells), unique physics only.
fitrows = [r for r in rows if r['use_fit'] and r['cls'] not in ('blob', 'pocket')]
print(f"\nstructured fit rows: {len(fitrows)}")

def build_matrix(rows_, features):
    X = np.array([[rw[f] for f in features] for rw in rows_], float)
    return X

FEATURES_ALL = ["x50", "margin", "N", "invN", "gmin", "gmean", "Dgate",
                "Bback", "gmeanN", "DgateN", "marg_gmean"]
FEATURES_RING = FEATURES_ALL + ["m", "defect", "wind"]

for r in fitrows:
    r['margin'] = r['x50'] - X_SKIRT
    r['invN'] = 1.0/r['N']
    r['gmeanN'] = r['gmean']*r['N']
    r['DgateN'] = r['Dgate']/r['N']
    r['marg_gmean'] = r['margin']*r['gmean']
    r['ln_xr'] = math.log(r['x50']/X_SKIRT) if r['x50'] > 0 else float('nan')

# --------------------------------------------------------------- regressors --
def ols(X, y):
    coef, *_ = np.linalg.lstsq(X, y, rcond=None)
    return coef

def r2(y, yhat):
    ss = np.sum((y - yhat)**2); tt = np.sum((y - np.mean(y))**2)
    return 1.0 - ss/tt if tt > 0 else float('nan')

def stlsq(X, y, lam=1e-6, thr=0.1, iters=20):
    """Standardized STLSQ: ridge + hard threshold on standardized coefs."""
    mu, sd = X.mean(0), X.std(0); sd[sd == 0] = 1
    Xs = (X - mu)/sd
    ymu, ysd = y.mean(), y.std() if y.std() > 0 else 1.0
    ys = (y - ymu)/ysd
    n, p = Xs.shape
    active = np.ones(p, bool)
    coef = np.zeros(p)
    for _ in range(iters):
        if not active.any(): break
        A = Xs[:, active]
        w = np.linalg.solve(A.T@A + lam*np.eye(A.shape[1]), A.T@ys)
        coef[:] = 0; coef[active] = w
        newactive = np.abs(coef) >= thr
        if (newactive == active).all(): break
        active = newactive
    # unstandardize
    beta = coef*ysd/sd
    b0 = ymu - np.dot(beta, mu)
    return b0, beta, active

def loto_cv(rows_, features, target, fitter):
    """Leave-one-topology-class-out CV. Returns pooled out-of-fold R2."""
    classes = sorted(set(rw['cls'] for rw in rows_))
    y_all, yh_all = [], []
    for c in classes:
        tr = [rw for rw in rows_ if rw['cls'] != c]
        te = [rw for rw in rows_ if rw['cls'] == c]
        if len(tr) < len(features) + 2:
            continue
        Xtr = build_matrix(tr, features); ytr = np.array([rw[target] for rw in tr])
        Xte = build_matrix(te, features); yte = np.array([rw[target] for rw in te])
        b0, beta = fitter(Xtr, ytr)
        yh = b0 + Xte@beta
        y_all += list(yte); yh_all += list(yh)
    y_all, yh_all = np.array(y_all), np.array(yh_all)
    return r2(y_all, yh_all), y_all, yh_all

def subset_search(rows_, features, target, kmax=2):
    """Exhaustive best-subset (k=1..kmax) with LOTO CV score. Honest for n~15."""
    y = np.array([rw[target] for rw in rows_])
    out = []
    for k in range(1, kmax+1):
        for combo in itertools.combinations(features, k):
            X = build_matrix(rows_, list(combo))
            if not np.all(np.isfinite(X)): continue
            Xd = np.column_stack([np.ones(len(X)), X])
            beta = ols(Xd, y)
            yh = Xd@beta
            fit_r2 = r2(y, yh)
            cv_r2, _, _ = loto_cv(rows_, list(combo), target,
                                  lambda Xt, yt: (ols(np.column_stack([np.ones(len(Xt)), Xt]), yt)[0],
                                                  ols(np.column_stack([np.ones(len(Xt)), Xt]), yt)[1:]))
            out.append((cv_r2, fit_r2, combo, beta))
    out.sort(key=lambda z: -(z[0] if z[0] == z[0] else -9))
    return out

# =========================================================== FIT 1: leak ==
print(f"\n=== FIT 1 — early leak rate  leak_early [pct/t.u.] (structured, n={len(fitrows)}) ===")
target = 'leak_early'
ok = [rw for rw in fitrows if np.isfinite(rw[target])]
res = subset_search(ok, FEATURES_ALL, target, kmax=2)
print("  top best-subset (CV = leave-one-class-out, pooled R2):")
for cv, fr, combo, beta in res[:6]:
    terms = " + ".join(f"{b:+.4g}*{f}" for b, f in zip(beta[1:], combo))
    print(f"    CV_R2={cv:+.3f}  fit_R2={fr:.3f}   leak = {beta[0]:+.4g} {terms}")

# rings/chains only, with closure/wind features
ringrows = [rw for rw in fitrows if rw['cls'] in ('ring_auto','ring_naive','chain','ring_exact')
            and np.isfinite(rw[target])]
print(f"\n  rings+chains only (n={len(ringrows)}), library + m/defect/wind:")
res2 = subset_search(ringrows, FEATURES_RING, target, kmax=2)
for cv, fr, combo, beta in res2[:6]:
    terms = " + ".join(f"{b:+.4g}*{f}" for b, f in zip(beta[1:], combo))
    print(f"    CV_R2={cv:+.3f}  fit_R2={fr:.3f}   leak = {beta[0]:+.4g} {terms}")

# STLSQ on the deduplicated library (margin/x50 and the two gate products are
# collinear pairs — one of each kept; lam raised so ridge actually damps)
FEATURES_STLSQ = ["x50", "N", "invN", "gmin", "gmean", "Dgate", "Bback"]
print("\n  STLSQ survivors (threshold sweep, standardized coefs, lam=1e-2):")
X = build_matrix(ok, FEATURES_STLSQ); y = np.array([rw[target] for rw in ok])
for thr in (0.1, 0.2, 0.3, 0.5, 0.8):
    b0, beta, act = stlsq(X, y, lam=1e-2, thr=thr)
    surv = [(FEATURES_STLSQ[i], beta[i]) for i in range(len(beta)) if act[i]]
    yh = b0 + X@beta
    cv_stlsq, _, _ = loto_cv(ok, [FEATURES_STLSQ[i] for i in range(len(act)) if act[i]] or FEATURES_STLSQ[:1],
                             target,
                             lambda Xt, yt: (ols(np.column_stack([np.ones(len(Xt)), Xt]), yt)[0],
                                             ols(np.column_stack([np.ones(len(Xt)), Xt]), yt)[1:]))
    print(f"    thr={thr}: fit_R2={r2(y,yh):.3f} CV_R2={cv_stlsq:+.3f}  " +
          (", ".join(f"{f}={b:+.4g}" for f, b in surv) if surv else "(none)"))

# ====================================================== FIT 2: death time ==
print("\n=== FIT 2 — death time (uncensored only) ===")
dead = [rw for rw in fitrows if rw['censored'] == 0 and rw['t_death'] > 0]
print(f"  uncensored structured deaths: n={len(dead)}  "
      f"({', '.join(rw['run'] for rw in dead)})")
for rw in dead:
    rw['ln_td'] = math.log(rw['t_death'])
res3 = subset_search(dead, ["x50", "margin", "ln_xr", "gmean", "gmin", "Dgate", "N"],
                     'ln_td', kmax=2)
print("  ln(t_death) best-subset:")
for cv, fr, combo, beta in res3[:6]:
    terms = " + ".join(f"{b:+.4g}*{f}" for b, f in zip(beta[1:], combo))
    print(f"    CV_R2={cv:+.3f}  fit_R2={fr:.3f}   ln_td = {beta[0]:+.4g} {terms}")

# the chosen (sparsest, best-CV) 1-term law, refit for the record:
X1 = np.column_stack([np.ones(len(dead)), [rw['ln_xr'] for rw in dead]])
y1 = np.array([rw['ln_td'] for rw in dead])
bb = ols(X1, y1)
resid = y1 - X1@bb; sd_ln = float(np.std(resid, ddof=2))
print(f"\n  CHOSEN LAW:  t_death = {math.exp(bb[0]):.0f} * (x50/{X_SKIRT})^{bb[1]:.3f}"
      f"   fit_R2={r2(y1, X1@bb):.3f}  sigma_ln={sd_ln:.3f} (x/ {math.exp(sd_ln):.2f})")

# censored consistency check for the chosen law
print("\n  censored rows vs the chosen law (consistent iff prediction >~ horizon;")
print("  a prediction far BELOW the horizon = the run beat the law = structure):")
for rw in fitrows:
    if rw['censored'] == 1 and np.isfinite(rw['ln_xr']):
        pred = math.exp(bb[0] + bb[1]*rw['ln_xr'])
        ratio = rw['T']/pred
        flag = "consistent" if pred >= 0.7*rw['T'] else f"BEAT THE LAW x{ratio:.1f}"
        print(f"    {rw['run']}: predicted {pred:.0f} vs alive-at {rw['T']:.0f}  [{flag}]")

# robustness: same law on t_struct (structured-phase end, kills the
# parked-crumb objection)
dead_s = [rw for rw in fitrows if rw['cens_struct'] == 0 and rw['t_struct'] > 100
          and np.isfinite(rw['ln_xr'])]
Xs_ = np.column_stack([np.ones(len(dead_s)), [rw['ln_xr'] for rw in dead_s]])
ys_ = np.log([rw['t_struct'] for rw in dead_s])
bbs = ols(Xs_, ys_)
print(f"\n  robustness on t_struct (last t with M_sum >= 1 voice-mass), n={len(dead_s)}:")
print(f"    t_struct = {math.exp(bbs[0]):.0f} * (x50/{X_SKIRT})^{bbs[1]:.3f}"
      f"   fit_R2={r2(ys_, Xs_@bbs):.3f}")
for rw in fitrows:
    if rw['cens_struct'] == 1 and np.isfinite(rw['ln_xr']):
        pred = math.exp(bbs[0] + bbs[1]*rw['ln_xr'])
        flag = "consistent" if pred >= 0.7*rw['T'] else f"BEAT THE LAW x{rw['T']/pred:.1f}"
        print(f"    censored: {rw['run']}: t_struct predicted {pred:.0f} vs alive-at "
              f"{rw['T']:.0f}  [{flag}]")

# ================================== FIT 3: the universal per-voice current ==
print("\n=== FIT 3 — the per-voice leak current c_eff = cap*(x50-x_skirt)/(t_death-50) ===")
print("  (the whole-life absolute leak per voice; censored rows give upper bounds)")
cvals = []
for rw in sorted(fitrows, key=lambda z: (z['censored'], z['cls'])):
    c = rw.get('c_eff', float('nan'))
    if not np.isfinite(c):
        continue
    if rw['censored'] == 0:
        print(f"    {rw['run']:<18} c_eff = {c:.3e}   (death {rw['t_death']:.0f})")
        cvals.append((rw, c))
    else:
        print(f"    {rw['run']:<18} c_eff <= {c:.3e}  (censored at {rw['T']:.0f})")
carr = np.array([c for _, c in cvals])
print(f"  uncensored: n={len(carr)}  median c0 = {np.median(carr):.3e}  "
      f"spread {carr.min():.2e}..{carr.max():.2e}  (MAD {np.median(np.abs(carr-np.median(carr))):.1e})")
# does ANY descriptor explain the residual spread of ln(c_eff)?
for rw, c in cvals:
    rw['ln_c'] = math.log(c)
res4 = subset_search([rw for rw, _ in cvals],
                     ["N", "invN", "gmin", "gmean", "Dgate", "Bback", "x50"],
                     'ln_c', kmax=1)
print("  ln(c_eff) one-term screens (CV = leave-one-class-out):")
for cv, fr, combo, beta in res4[:4]:
    print(f"    CV_R2={cv:+.3f} fit_R2={fr:.3f}  ln_c = {beta[0]:+.4g} "
          f"{beta[1]:+.4g}*{combo[0]}")
c0 = float(np.median(carr))

# ============================================== per-class reference table ==
print("\n=== per-class summary (fit rows) ===")
hdr = (f"{'run':<18}{'cls':<11}{'N':>4}{'x50':>8}{'gmin':>10}{'gmean':>7}"
       f"{'Dgate':>8}{'leak%':>9}{'c_eff':>10}{'t_death':>9}{'cens':>5}"
       f"{'t_struct':>9}{'t_skirt':>8}")
print(hdr)
def s(v, w, prec="{:.4g}"):
    if isinstance(v, float) and v != v: return "nan".rjust(w)
    if isinstance(v, float): return prec.format(v).rjust(w)
    return str(v).rjust(w)
for rw in sorted(fitrows + [r for r in rows if r['cls'] == 'blob' and r['use_fit']],
                 key=lambda z: (z['cls'], z['run'])):
    print(f"{rw['run']:<18}{rw['cls']:<11}" + s(rw['N'], 4, "{:.0f}") +
          s(rw.get('x50', float('nan')), 8) + s(rw['gmin'], 10, "{:.3g}") +
          s(rw['gmean'], 7, "{:.3f}") + s(rw['Dgate'], 8, "{:.3g}") +
          s(rw.get('leak_early', float('nan')), 9) +
          s(rw.get('c_eff', float('nan')), 10, "{:.2e}") +
          s(rw['t_death'], 9, "{:.0f}") + s(rw.get('censored', -1), 5) +
          s(rw.get('t_struct', float('nan')), 9, "{:.0f}") +
          s(rw.get('t_skirt', float('nan')), 8, "{:.0f}"))

# ===================================================== PREDICTIONS block ==
print("\n=== PHASE-2 PREDICTIONS (evaluated live from the fitted laws) ===")
# Law A (primary): t_death = A0 * (x50/x_skirt)^p        [FIT 2 1-term]
# Law B (current): t_death = 50 + cap*(x50-x_skirt)/c0   [FIT 3 universal]
# Uncertainty: fit sigma_ln (x/1.1) is NOT the honest band — foam chaos is
# +-30% (MASS.md ensemble) and only one a2 seed died. Quote x/1.5.
UNC = 1.5
print(f"  Law A: t_death = {math.exp(bb[0]):.0f}*(x50/{X_SKIRT})^{bb[1]:.3f}   "
      f"Law B: t_death = 50 + {CAP:.1f}*(x50-{X_SKIRT})/{c0:.2e}")
print(f"  band quoted: x/ {UNC} (fit sigma_ln {sd_ln:.2f} + foam chaos +-30%)")
def tdA(x):  return math.exp(bb[0] + bb[1]*math.log(x/X_SKIRT))
def tdB(x):  return 50 + CAP*(x - X_SKIRT)/c0
print("\n  scoring curve (Law A | Law B):")
for x in (0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.50, 0.60, 0.70):
    print(f"    x50={x:.2f}:  {tdA(x):6.0f} | {tdB(x):6.0f}")
cands = [
    ("(a) cube shell retuned, gates min>=0.95", 0.375,
     "same heavy-cube geometry (abar 1.586 -> x50~0.375 measured on the v3 twin); "
     "load-line null: gate quality bought nothing across shells v1/v2/v3"),
    ("(b) hex prism N=12 struts-only", 0.40,
     "assume ring_m-exact seeding at comp12-class load; bipartite prism, all links mutual"),
    ("(c) wound tube 2x ring12 co-rotating", 0.40,
     "comp12-class links + exact m=1 axial rungs; WOUND-MUTUAL class -> expect "
     "the comp12 exception: beats Law A by >=2.3x (comp12 bound)"),
    ("(d) ring12 + 2 consonant chords", 0.40,
     "comp12 + chords; extra mutual closure paths -> between Law A and the comp12 excess"),
    ("(e) torus net N~24", 0.35,
     "2D closed bipartite net; first N-scaling test of the c0*N current"),
]
for label, x, note in cands:
    a, b_ = tdA(x), tdB(x)
    print(f"\n  {label}\n    x50={x:.3f} -> Law A {a:.0f}, Law B {b_:.0f}  "
          f"(band x/ {UNC})\n    {note}")
print("\nDone. corpus.tsv + printed tables are the REGRESS.md source.")
