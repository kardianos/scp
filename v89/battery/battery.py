#!/usr/bin/env python3
# v89 unification battery — the ROADMAP §6.5 protocol, operational.
#
# One LAWS file, shared byte-identically by every experiment in a run;
# apparatus files may not touch law keys (purity-checked before anything
# executes); acceptance is physics, not bytes. Variants (laws_V*.cfg)
# compete as complete law tables — switching a law per experiment is
# structurally impossible here.
#
# usage: battery.py --laws laws_V1.cfg [--jobs 8] [--only e6_pairs ...]
#                   [--skip-run]   (re-parse existing logs only)

import argparse
import concurrent.futures as cf
import os
import re
import subprocess
import sys

ROOT = os.path.dirname(os.path.abspath(__file__))
V89 = os.path.dirname(ROOT)
# DEFAULT kernel is integer-ledger cellfab.c (mode 3). FP64 reference:
# cellfabf.c → build separately as cellfabf for A/B only.
SRC = os.path.join(V89, "cellfab.c")
BIN = os.path.join(V89, "cellfab")

LAW_KEYS = {
    "C", "dmin", "r0", "rjit", "w1", "w2", "q_detune", "gamma_res",
    "gamma_res_m", "p_gate", "lock_floor", "k_dep", "k_dep_m", "cap",
    "e_s0", "es_floor", "e_cond", "f_conv", "f_evap", "s_pull",
    "kappa_lock", "kappa_align", "kappa_freq", "kappa_reac", "sigma_tumble",
    "comb_limit", "rough_k",
    "gamma_rough", "mob_sym", "mob_floor", "field_J", "quant_A0",
    "quant_mode", "s_k", "s_disp",
}
BANNED_IN_APPARATUS = LAW_KEYS | {"e_click"}  # e_click: retired into eps(w)

EXPS = ["e1_conserve", "e2_pulse", "e3a_blob", "e3b_blob_tilt", "e4_curve",
        "e5_bell", "e6_pairs", "e7_tune", "e8_comma", "e9_fifth", "d1_slit",
        "t1_tonomura", "q2_eraser", "t4_hom", "qt_lo", "qt_hi",
        "p1_beam", "g1_footprint", "g3_shadow", "g4_throughput"]
# p2_press: recorded, NOT gated (ratchet rule: tests enter as they PASS).
# Measured 2026-07-28: absorption happens but the field's momentum does
# not survive the conversion — recoil deficit ~100x (predicted 0.13C,
# measured <=1e-3; the dense translation ceiling itself is ~5e-3).
# Third S2-full acceptance criterion (ROADMAP §7). Run via --only p2_press.

LAWS = {}

DRIFT_MAX = 5e-14


def cfg_keys(path):
    keys = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if "=" in line:
                keys.append(line.split("=", 1)[0].strip())
    return keys


def purity_check(laws_path):
    errs = []
    lk = cfg_keys(laws_path)
    if set(lk) != LAW_KEYS:
        missing = LAW_KEYS - set(lk)
        extra = set(lk) - LAW_KEYS
        if missing:
            errs.append(f"laws file missing law keys: {sorted(missing)}")
        if extra:
            errs.append(f"laws file has non-law keys: {sorted(extra)}")
    if len(lk) != len(set(lk)):
        errs.append("laws file has duplicate keys")
    for e in EXPS:
        ap = os.path.join(ROOT, "apparatus", e + ".cfg")
        bad = [k for k in cfg_keys(ap) if k in BANNED_IN_APPARATUS]
        if bad:
            errs.append(f"apparatus {e}: law keys present: {bad}")
    return errs


def build_kernel():
    # Serialize concurrent builds (e.g. several battery.py invocations in
    # parallel) via an flock on a lock file beside the binary. Without this
    # two builds race on the shared BIN path -> ETXTBSY "Text file busy".
    import fcntl
    lockpath = BIN + ".build.lock"
    with open(lockpath, "w") as lf:
        fcntl.flock(lf, fcntl.LOCK_EX)
        # Build to a temp path then atomically rename over BIN: a process
        # already executing the old binary keeps running unaffected, and any
        # concurrent exec attempt sees either the old or the new binary,
        # never a half-written one.
        tmp = BIN + ".build.tmp"
        cmd = ["gcc", "-O2", "-march=native", "-fopenmp", "-o", tmp, SRC, "-lm"]
        r = subprocess.run(cmd, capture_output=True, text=True)
        if r.returncode != 0:
            sys.exit(f"kernel build failed:\n{r.stderr}")
        os.replace(tmp, BIN)


RUN_THREADS = "1"   # per-run OpenMP threads; set from --jobs in main
EXTRAS = []         # apparatus k=v lines appended to every merged cfg
EXTRAS_FOR = {}     # per-experiment apparatus k=v lines (variant sets)


def run_one(laws_path, variant, name):
    vdir = os.path.join(ROOT, "runs", variant)
    os.makedirs(os.path.join(vdir, "cfg"), exist_ok=True)
    cfg = os.path.join(vdir, "cfg", name + ".cfg")
    log = os.path.join(vdir, name + ".log")
    with open(laws_path) as fh:
        laws = fh.read()
    with open(os.path.join(ROOT, "apparatus", name + ".cfg")) as fh:
        app = fh.read()
    with open(cfg, "w") as fh:
        fh.write(laws + "\n# --- apparatus ---\n" + app)
        if EXTRAS:
            fh.write("\n# --- extras ---\n" + "\n".join(EXTRAS) + "\n")
        if name in EXTRAS_FOR:
            fh.write("# --- extras (this experiment) ---\n"
                     + "\n".join(EXTRAS_FOR[name]) + "\n")
    env = dict(os.environ, OMP_NUM_THREADS=RUN_THREADS)
    with open(log, "w") as fh:
        r = subprocess.run([BIN, cfg], stdout=fh, stderr=subprocess.STDOUT,
                           env=env)
    return name, r.returncode


# ---------------------------------------------------------------- parsing

def grab(log, pat, idx=None):
    """last match of pat in log; floats of all groups (or group idx)."""
    m = None
    for m_ in re.finditer(pat, log):
        m = m_
    if m is None:
        return None
    g = [float(x) for x in m.groups()]
    return g if idx is None else g[idx]


def fmt(x, spec="%.2e"):
    """format a possibly-None numeric for result strings (None -> 'NA')."""
    return spec % x if x is not None else "NA"


def last_diag_em(log):
    em = None
    for line in log.splitlines():
        if line.startswith("#") or "\t" not in line:
            continue
        parts = line.split("\t")
        if len(parts) >= 15:
            try:
                em = float(parts[2])
            except ValueError:
                pass
    return em


def pair_series(log):
    rows = []
    for line in log.splitlines():
        if not line.startswith("# PAIR "):
            continue
        kv = dict(m.group(1, 2) for m in re.finditer(r"(\w+)=([^\s]+)", line))
        try:
            rows.append({
                "t": float(kv["t"]), "p": int(kv["p"]), "d": float(kv["d"]),
                "delta": float(kv["delta"]), "gg": float(kv["gg"]),
                "ratio": float(kv["ratio"]), "shed": float(kv["shed"]),
                "ret": float(kv["ret"]),
            })
        except (KeyError, ValueError):
            pass
    return rows


def qatom_points(log):
    return [(float(w), float(e)) for w, e in
            re.findall(r"# QATOM t=[-\d.e+]+ dir=\w+ w=([\d.e+-]+) e=([\d.e+-]+)", log)]


def conservation_ok(log):
    d = grab(log, r"# RESULT conservation .*rel_drift=([-\d.e+]+)", 0)
    # Return nan (not None) on missing lines so downstream f-strings that
    # format d with a spec (:.2e) never crash on incomplete logs.
    if d is None:
        return float("nan"), False
    return d, abs(d) < DRIFT_MAX


def diag_table(log):
    """diag rows as column lists (P-battery: needs the appended centroid
    and flux-moment columns, indices 19.. — older logs return fewer)."""
    rows = []
    for line in log.splitlines():
        if line.startswith("#") or "\t" not in line:
            continue
        parts = line.split("\t")
        try:
            rows.append([float(x) for x in parts])
        except ValueError:
            pass
    return rows


def linfit(ts, ys):
    n = len(ts)
    if n < 3:
        return 0.0
    mt = sum(ts) / n
    my = sum(ys) / n
    num = sum((t - mt) * (y - my) for t, y in zip(ts, ys))
    den = sum((t - mt) ** 2 for t in ts)
    return num / den if den > 0 else 0.0


# ------------------------------------------------------- acceptance checks

def chk_e1(log):
    # e1's claim: the ledger closes with every mechanism live (conversions
    # included when peaks clear the atom). Threshold physics is qt_lo/qt_hi.
    d, ok = conservation_ok(log)
    nq = len(qatom_points(log))
    em = last_diag_em(log)
    return ok, f"drift={d:.2e} qatoms={nq} Em_final={fmt(em, "%.3g")}"


def chk_qt_lo(log):
    d, ok = conservation_ok(log)
    nq = len(qatom_points(log))
    em = last_diag_em(log)
    sub = nq == 0 and (em is not None and em < 1e-9)
    return ok and sub, f"drift={d:.2e} qatoms={nq} Em_final={fmt(em, "%.3g")}"


def chk_e2(log):
    d, ok = conservation_ok(log)
    v = grab(log, r"# RESULT front_speed v=[\d.-]+ v_over_C=([\d.-]+)", 0)
    return ok and v is not None and v >= 0.3, f"drift={d:.2e} v/C={v}"


def chk_e3a(log):
    d, ok = conservation_ok(log)
    s = grab(log, r"# RESULT blob_drift .*speed=([\d.-]+)", 0)
    return ok and s is not None and s <= 0.002, f"drift={d:.2e} speed={s}"


def chk_e3b(log):
    d, ok = conservation_ok(log)
    g = grab(log, r"# RESULT blob_drift .*speed=([\d.-]+) cos_to_kdir=([\d.-]+)")
    if g is None:
        return False, f"drift={d:.2e} no blob_drift"
    s, c = g
    return ok and s >= 0.003 and c >= 0.8, f"drift={d:.2e} speed={s} cos={c}"


def chk_e4(log):
    d, ok = conservation_ok(log)
    g = grab(log, r"# RESULT curvature_fit defA_per_Em=([\d.-]+) r2=([\d.-]+)")
    if g is None:
        return False, f"drift={d:.2e} no curvature_fit"
    k, r2 = g
    return ok and k > 0 and r2 >= 0.95, f"drift={d:.2e} defA/Em={k} r2={r2}"


def chk_e5(log):
    g = grab(log, r"# RESULT S_joint=([\d.-]+) S_lhv=([\d.-]+)")
    if g is None:
        return False, "no S line"
    s, lhv = g
    return s >= 2.7 and lhv <= 2.1, f"S={s} S_lhv={lhv}"


def tongue(log):
    return grab(log, r"# RESULT pair_tongue \|delta\|<0\.6: n=(\d+) gg=([\d.]+) "
                     r"ret=([\d.]+) fl=([\d.]+) \| 0\.6-1\.2: n=(\d+) gg=([\d.]+) "
                     r"ret=([\d.]+) fl=([\d.]+) \| >1\.2: n=(\d+) gg=([\d.]+)")


def chk_e6(log):
    d, ok = conservation_ok(log)
    g = tongue(log)
    if g is None:
        return False, f"drift={d:.2e} no tongue"
    n_in, gg_in = g[0], g[1]
    n_far, gg_far = g[8], g[9]
    far_ok = n_far < 3 or gg_far <= 0.15
    return (ok and n_in >= 10 and gg_in >= 0.5 and far_ok), \
        f"drift={d:.2e} in: n={n_in:.0f} gg={gg_in} far: n={n_far:.0f} gg={gg_far}"


def chk_e7(log):
    # P4 (CONSONANCE IV): the computed tuning curve x*(d) puts pairs ON
    # the rung across all separations — scored on the pairs the law
    # table itself says can exist. Two claims, both scored:
    #  (1) living pairs (ret >= 0.3) end delta-pinned (frac >= 0.75);
    #  (2) the vacuum skirt PREDICTION: a voice within ~2*Gamma of the
    #      vacuum pitch dissolves into the room — every death must lie
    #      inside 1.5x the computed boundary x_skirt = 2G/(q*(w2-2G)).
    d, ok = conservation_ok(log)
    rows = pair_series(log)
    if not rows:
        return False, f"drift={d:.2e} no pairs"
    import math
    w2 = float(LAWS["w2"]); q = float(LAWS["q_detune"])
    gm = float(LAWS["gamma_res_m"])
    xskirt = 2 * gm / (q * (w2 - 2 * gm))
    tN = max(r["t"] for r in rows)
    seen = set(); fin = []
    for r in rows:
        if r["t"] == tN and r["p"] not in seen:
            seen.add(r["p"]); fin.append(r)
    for r in fin:
        r["xs"] = (w2 * r["d"] / math.pi - 1) / q
    alive = [r for r in fin if r["ret"] >= 0.3]
    dead = [r for r in fin if r["ret"] < 0.3]
    frac = (sum(1 for r in alive if abs(r["delta"]) < 0.15) / len(alive)) \
        if alive else 0
    gg = sum(r["gg"] for r in alive) / len(alive) if alive else 0
    dead_max = max((r["xs"] for r in dead), default=0.0)
    skirt_ok = dead_max <= 1.5 * xskirt
    return (ok and len(fin) >= 40 and len(alive) >= 30 and frac >= 0.75
            and skirt_ok), \
        f"drift={d:.2e} alive={len(alive)}/{len(fin)} frac|delta|<0.15={frac:.2f} " \
        f"gg={gg:.2f} dead_x*max={dead_max:.3f} skirt1.5={1.5*xskirt:.3f}"


def chk_e8(log):
    d, ok = conservation_ok(log)
    rows = pair_series(log)
    if not rows:
        return False, f"drift={d:.2e} no pairs"
    t0 = min(r["t"] for r in rows)
    tN = max(r["t"] for r in rows)
    d0 = {r["p"]: abs(r["delta"]) for r in rows if r["t"] == t0}
    sh = {r["p"]: r["shed"] for r in rows if r["t"] == tN}
    common = sorted(set(d0) & set(sh), key=lambda p: d0[p])
    if len(common) < 10:
        return False, f"drift={d:.2e} too few pairs"
    half = len(common) // 2
    lo = sum(sh[p] for p in common[:half]) / half
    hi = sum(sh[p] for p in common[half:]) / (len(common) - half)
    tot = sum(sh.values())
    return (ok and tot > 1.0 and hi > lo), \
        f"drift={d:.2e} shed_lo|delta0|={lo:.3f} shed_hi|delta0|={hi:.3f} total={tot:.1f}"


def chk_e9(log):
    d, ok = conservation_ok(log)
    rows = pair_series(log)
    locked = {}
    for r in rows:
        if r["gg"] > 0.3 and 1.45 < r["ratio"] < 1.55:
            locked[r["t"]] = locked.get(r["t"], 0) + 1
    n10 = max((n for t, n in locked.items() if 5 <= t <= 15), default=0)
    t_last = max((t for t, n in locked.items() if n >= 1), default=0)
    return (ok and n10 >= 5 and t_last >= 20), \
        f"drift={d:.2e} locked@t10={n10} t_last={t_last:.0f}"


def chk_d1(log):
    d, ok = conservation_ok(log)
    v = grab(log, r"# RESULT fringe_visibility_central=([\d.-]+)", 0)
    ex = grab(log, r"# RESULT screen_exposure=([\d.-]+)", 0)
    return (ok and v is not None and v >= 0.25 and ex is not None and ex >= 10), \
        f"drift={d:.2e} V={v} exposure={ex}"


def chk_t1(log):
    d, ok = conservation_ok(log)
    v = grab(log, r"# RESULT qclick_visibility_central=([\d.-]+)", 0)
    ncl = len(re.findall(r"^# QCLICK ", log, re.M))
    return (ok and v is not None and v >= 0.2 and ncl >= 100), \
        f"drift={d:.2e} V_click={v} clicks={ncl}"


def chk_q2(log):
    d, ok = conservation_ok(log)
    g = grab(log, r"# RESULT analyzer ledgers: totalA=[\d.-]+ totalB=[\d.-]+ "
                  r"V_A=([\d.-]+) V_B=([\d.-]+)")
    if g is None:
        return False, f"drift={d:.2e} no analyzer line"
    va, vb = g
    return ok and va >= 0.2 and vb <= -0.1, f"drift={d:.2e} V_A={va} V_B={vb}"


def chk_hom(log):
    d, ok = conservation_ok(log)
    sp = grab(log, r"# RESULT hom coupler_split\(A cross->D2\)=([\d.-]+)", 0)
    g = grab(log, r"# RESULT hom g_boson=([\d.-]+) g_fermion=([\d.-]+)")
    if g is None or sp is None:
        return False, f"drift={d:.2e} incomplete hom results"
    gb, gf = g
    return (ok and gb <= 0.47 and gf >= 0.53 and 0.35 <= sp <= 0.65), \
        f"drift={d:.2e} split={sp} g_b={gb} g_f={gf}"


def chk_qt_hi(log):
    d, ok = conservation_ok(log)
    nq = len(qatom_points(log))
    em = last_diag_em(log)
    return (ok and nq > 0 and em is not None and em > 0.5), \
        f"drift={d:.2e} qatoms={nq} Em_final={fmt(em, "%.3g")}"


def chk_p1(log):
    # momentum of light: ballistic direction persistence of the field
    # energy centroid + the flux moment forward-biased during transit
    d, ok = conservation_ok(log)
    rows = [r for r in diag_table(log) if len(r) >= 28 and 3.0 <= r[0] <= 18.0]
    if len(rows) < 8:
        return False, f"drift={d:.2e} too few diag rows"
    ts = [r[0] for r in rows]
    vx = linfit(ts, [r[19] for r in rows])
    vy = linfit(ts, [r[20] for r in rows])
    vz = linfit(ts, [r[21] for r in rows])
    import math
    sp = math.sqrt(vx * vx + vy * vy + vz * vz)
    cos = vx / sp if sp > 0 else 0
    C = float(LAWS["C"])
    ffx = sum(r[25] for r in rows) / len(rows)
    return (ok and cos >= 0.9 and 0.3 * C <= sp <= 1.0 * C and ffx > 0), \
        f"drift={d:.2e} v/C={sp / C:.3f} cos={cos:.3f} ffx={ffx:.3g}"


def chk_p2(log):
    # radiation pressure: still before arrival; absorbs; recoils ALONG
    # the packet direction (pushed, never pulled). Self-controlled.
    d, ok = conservation_ok(log)
    rows = [r for r in diag_table(log) if len(r) >= 25]
    if len(rows) < 60:
        return False, f"drift={d:.2e} too few diag rows"
    pre = [r for r in rows if 2.0 <= r[0] <= 8.0]
    post = [r for r in rows if 30.0 <= r[0] <= 115.0]
    v_pre = linfit([r[0] for r in pre], [r[8] for r in pre])
    v_post = linfit([r[0] for r in post], [r[8] for r in post])
    em0 = rows[0][2]
    em_pk = max(r[2] for r in rows)
    absorbed = em_pk - em0
    return (ok and abs(v_pre) <= 1e-3 and absorbed > 5.0 and v_post >= 1.5e-3), \
        f"drift={d:.2e} v_pre={v_pre:.2e} v_post={v_post:.2e} " \
        f"absorbed={absorbed:.1f}"


def chk_g1(log):
    # the gravitational footprint: a sealed mass maintains a graded space
    # depression against transport (displacement holds it open)
    d, ok = conservation_ok(log)
    m = re.search(r"# RESULT es_shells((?: [\d.]+:[\d.]+)+)", log)
    s = grab(log, r"# RESULT blob_drift .*speed=([\d.-]+)", 0)
    if m is None or s is None:
        return False, f"drift={d:.2e} incomplete"
    sh = [float(p.split(":")[1]) for p in m.group(1).split()]
    core = sh[0]
    far = sum(sh[-3:]) / 3
    dip = core / far if far > 0 else 1
    return (ok and dip <= 0.8 and s <= 0.002), \
        f"drift={d:.2e} core/far={dip:.3f} sealed_speed={s}"


def chk_g3(log):
    # occultation: the transmitted beam exits displaced AWAY from the
    # mass (the loaded core is a mirror — matter is emergently opaque)
    d, ok = conservation_ok(log)
    g = grab(log, r"# RESULT exit_face E=([\d.e+-]+) y=([\d.-]+) z=([\d.-]+)")
    if g is None:
        return False, f"drift={d:.2e} no exit_face"
    fE, fy, _ = g
    return (ok and fE >= 1.0 and fy >= 24.0), \
        f"drift={d:.2e} E_exit={fE:.2f} y_exit={fy:.2f} (launch 22.5, mass at 18)"


def rad_rows(log):
    rows = []
    for line in log.splitlines():
        if not line.startswith("# RAD"):
            continue
        t = float(re.search(r"t=([\d.]+)", line).group(1))
        sh = [tuple(float(x) for x in m.groups()) for m in
              re.finditer(r"([\d.]+):([-\d.e+]+),([-\d.e+]+),([-\d.e+]+),"
                          r"([-\d.e+]+)", line)]
        if len(sh) == 8:
            rows.append((t, sh))
    return rows


def chk_g4(log):
    # blob space throughput: NO steady monopole — far-shell space flux is
    # mass-rate-driven (equilibration wind / refill), decays with the
    # transients, and stays subdominant to the radiative field channel
    # while the blob leaks. A steady 1/r-sourcing accretion flux would
    # violate both bars. (The 1/r far field awaits a STABLE particle
    # whose internal cycle gives throughput at constant mass.)
    d, ok = conservation_ok(log)
    rows = rad_rows(log)
    if len(rows) < 100:
        return False, f"drift={d:.2e} too few RAD rows"
    def wavg(lo, hi, col, k):
        w = [sh[k][col] for t, sh in rows if lo <= t <= hi]
        return sum(w) / len(w) if w else 0.0
    sf_e = wavg(20, 60, 1, 5)
    sf_l = wavg(100, 160, 1, 5)
    ff_l = wavg(100, 160, 2, 5)
    def blobm(sh): return sum(s[4] for s in sh[:4])
    late = [(t, sh) for t, sh in rows if 100 <= t <= 160]
    mdot = (blobm(late[-1][1]) - blobm(late[0][1])) / (late[-1][0] - late[0][0])
    # RELAXED 2026-07-31 (user decision, until a stable mass exists):
    # decay factor 0.5 -> 0.75 (sparser substrates relax slower; the old
    # windows encoded the foam's timescale) and the radiative-subdominance
    # clause -> mass-rate BOOKKEEPING |sf_late| <= |dM/dt| (the actual
    # claim: far space flux is leak bookkeeping, no steady monopole; on
    # quiet substrates the leak legitimately stops radiating). The sharp
    # form returns with g5 once a C1 particle exists.
    decays = abs(sf_l) <= 0.75 * abs(sf_e) + 1e-6
    booked = abs(sf_l) <= abs(mdot) + 1e-6
    return (ok and decays and booked), \
        f"drift={d:.2e} sflux(r8.25) early={sf_e:+.3g} late={sf_l:+.3g} " \
        f"fflux_late={ff_l:+.3g} dM/dt={mdot:+.3f}"


CHECKS = {"e1_conserve": chk_e1, "e2_pulse": chk_e2, "e3a_blob": chk_e3a,
          "e3b_blob_tilt": chk_e3b, "e4_curve": chk_e4, "e5_bell": chk_e5,
          "e6_pairs": chk_e6, "e7_tune": chk_e7, "e8_comma": chk_e8,
          "e9_fifth": chk_e9, "d1_slit": chk_d1, "t1_tonomura": chk_t1,
          "q2_eraser": chk_q2, "t4_hom": chk_hom, "qt_lo": chk_qt_lo,
          "qt_hi": chk_qt_hi, "p1_beam": chk_p1, "p2_press": chk_p2,
          "g1_footprint": chk_g1, "g3_shadow": chk_g3,
          "g4_throughput": chk_g4}


def chk_linearity(all_logs, a0):
    """cross-battery: every fired grain is a whole atom of eps(w) = A0*w/2pi,
    across the spread of detuned frequencies the runs produced."""
    pts = []
    for log in all_logs.values():
        pts += qatom_points(log)
    if len(pts) < 10:
        return None, f"only {len(pts)} QATOM points (need 10)"
    ws = sorted(w for w, _ in pts)
    spread = (ws[-1] - ws[0]) / ws[-1] if ws[-1] > 0 else 0
    worst = 0.0
    for w, e in pts:
        eps = a0 * w / (2 * 3.141592653589793)
        k = e / eps
        dev = abs(k - round(k))
        worst = max(worst, dev)
    ok = worst < 1e-6 and spread >= 0.10
    return ok, f"n={len(pts)} w=[{ws[0]:.3f},{ws[-1]:.3f}] spread={spread:.2f} max_offgrid={worst:.2e}"


# ------------------------------------------------------------------- main

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--laws", required=True)
    ap.add_argument("--jobs", type=int, default=8)
    ap.add_argument("--only", nargs="*", default=None)
    ap.add_argument("--skip-run", action="store_true")
    ap.add_argument("--extra", nargs="*", default=[],
                    help="apparatus k=v lines appended to every merged cfg "
                         "(e.g. geom_relax=400); NOT law keys")
    ap.add_argument("--extra-for", nargs="+", action="append", default=[],
                    dest="extra_for",
                    help="EXP k=v [k=v ...]: apparatus lines appended only "
                         "to that experiment's merged cfg (the variant "
                         "apparatus set; bars never change here)")
    ap.add_argument("--tag", default="",
                    help="suffix for the variant run dir (keeps baselines)")
    args = ap.parse_args()

    laws_path = os.path.join(ROOT, args.laws)
    variant = re.sub(r"^laws_|\.cfg$", "", os.path.basename(args.laws))
    if args.tag:
        variant = variant + "_" + args.tag
    global EXTRAS, EXTRAS_FOR
    EXTRAS = args.extra
    EXTRAS_FOR = {ef[0]: ef[1:] for ef in args.extra_for}
    exps = args.only if args.only else EXPS

    errs = purity_check(laws_path)
    if errs:
        for e in errs:
            print("PURITY FAIL:", e)
        sys.exit(2)
    print(f"# purity OK: laws={sorted(LAW_KEYS).__len__()} keys, "
          f"apparatus clean; variant={variant}")

    global RUN_THREADS
    RUN_THREADS = str(max(1, (os.cpu_count() or 8) // args.jobs))

    if not args.skip_run:
        build_kernel()
        with cf.ThreadPoolExecutor(max_workers=args.jobs) as ex:
            futs = {ex.submit(run_one, laws_path, variant, e): e for e in exps}
            for fu in cf.as_completed(futs):
                name, rc = fu.result()
                print(f"# run done: {name} rc={rc}")

    laws = dict((k.strip(), v.strip()) for k, v in
                (l.split("=", 1) for l in open(laws_path)
                 if "=" in l and not l.strip().startswith("#")))
    a0 = float(laws["quant_A0"])
    global LAWS
    LAWS = laws

    logs = {}
    for e in exps:
        p = os.path.join(ROOT, "runs", variant, e + ".log")
        logs[e] = open(p).read() if os.path.exists(p) else ""

    print(f"\n== battery {variant} ==")
    fails = 0
    rows = []
    for e in exps:
        if not logs[e]:
            rows.append((e, False, "no log"))
            fails += 1
            continue
        ok, note = CHECKS[e](logs[e])
        rows.append((e, ok, note))
        fails += 0 if ok else 1
    lin_ok, lin_note = chk_linearity(logs, a0)
    if lin_ok is not None:
        rows.append(("LIN grain=eps(w)", lin_ok, lin_note))
        fails += 0 if lin_ok else 1
    else:
        rows.append(("LIN grain=eps(w)", True, "SKIP " + lin_note))

    w = max(len(r[0]) for r in rows)
    for name, ok, note in rows:
        print(f"{name:<{w}}  {'PASS' if ok else 'FAIL'}  {note}")
    print(f"== {variant}: {len(rows)-fails}/{len(rows)} pass ==")
    with open(os.path.join(ROOT, "runs", variant, "summary.tsv"), "w") as fh:
        for name, ok, note in rows:
            fh.write(f"{name}\t{'PASS' if ok else 'FAIL'}\t{note}\n")
    sys.exit(1 if fails else 0)


if __name__ == "__main__":
    main()
