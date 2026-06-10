#!/usr/bin/env python3
"""v67 Workpackage D — numerical side of DEBROGLIE.md.

1. Action-quantum table: E/(omega*Q) across the v66 branch (scan.tsv),
   h_eff = E/omega vs Q.
2. Verify the exact identity  E - omega*Q = int f'^2 dV  [thm, deBroglie.mac D8]
   on every radial profile in v66/results/ (PASS/FAIL).
3. Clock-gradient drift number + V51 order-of-magnitude comparison.
4. Orbit-quantization estimate table (Bohr-Sommerfeld phase closure).

Run: python3 deBroglie.py   (output: deBroglie_py.out)
"""
import glob, math, os, re, sys

ROOT = "/home/d/code/scp"
RES = os.path.join(ROOT, "v66", "results")
M2, MU, KAP = 2.25, -41.345, 50.0
M = math.sqrt(M2)
nfail = 0


def chk(name, ok, detail=""):
    global nfail
    print(("PASS " if ok else "FAIL ") + name + ("  " + detail if detail else ""))
    if not ok:
        nfail += 1


def simpson(y, h):
    n = len(y) - 1
    if n % 2 == 1:          # drop last interval into a trapezoid
        return simpson(y[:-1], h) + 0.5 * h * (y[-2] + y[-1])
    s = y[0] + y[-1] + 4.0 * sum(y[1:-1:2]) + 2.0 * sum(y[2:-1:2])
    return s * h / 3.0


# ----------------------------------------------------------------------
# 1. scan.tsv action-quantum table
# ----------------------------------------------------------------------
rows = []
with open(os.path.join(RES, "scan.tsv")) as fh:
    hdr = fh.readline().split()
    for line in fh:
        t = line.split()
        if len(t) < 6:
            continue
        rows.append(dict(zip(hdr, [float(x) for x in t])))

print("=" * 78)
print("1. ACTION-QUANTUM TABLE (v66/results/scan.tsv)")
print("   h_eff = E/omega (the ball's effective hbar);  ratio = E/(omega*Q)")
print("   [thm] E - omega*Q = int f'^2 dV  =>  ratio = 1 + G/(omega*Q) > 1")
print("=" * 78)
print(f"{'omega':>7} {'Q':>10} {'E':>10} {'h_eff=E/w':>10} {'E/(wQ)':>8} "
      f"{'E/(mQ)':>8} {'branch':>8}")
ratios_stable, ratios_all = [], []
for r in rows:
    w, Q, E = r["omega"], r["Q"], r["E"]
    ratio = E / (w * Q)
    stable = (r["dQ_domega_sign"] < 0) and (E < M * Q)   # fission + evaporation
    tag = "STABLE" if stable else ("thick" if r["dQ_domega_sign"] > 0 else "evap>")
    ratios_all.append((w, ratio))
    if stable:
        ratios_stable.append((w, ratio))
    print(f"{w:7.4f} {Q:10.1f} {E:10.1f} {E/w:10.1f} {ratio:8.4f} "
          f"{E/(M*Q):8.4f} {tag:>8}")

lo = min(x[1] for x in ratios_stable)
hi = max(x[1] for x in ratios_stable)
print(f"\nstable branch (dQ/dw<0 AND E<mQ): omega in "
      f"[{ratios_stable[0][0]:.4f},{ratios_stable[-1][0]:.4f}]")
print(f"E/(omega*Q) range on stable branch: [{lo:.4f}, {hi:.4f}]  "
      f"(deviation from 1: {100*(lo-1):.2f}% .. {100*(hi-1):.2f}%)")
chk("P1_E_eq_wQ_to_5pct_on_stable_branch", 1.0 < lo and hi < 1.05,
    f"range [{lo:.4f},{hi:.4f}]")
hi_all = max(x[1] for x in ratios_all)
chk("P1b_E_eq_wQ_to_5pct_whole_scan", hi_all < 1.05, f"max {hi_all:.4f}")

# ----------------------------------------------------------------------
# 2. profile check of E - wQ = int f'^2 dV  (4 pi int f'^2 r^2 dr)
# ----------------------------------------------------------------------
print()
print("=" * 78)
print("2. PROFILE VERIFICATION of E - omega*Q = int f'^2 dV  [thm -> verified]")
print("=" * 78)
print("   (omega <= 1.320 thin-wall edge rows excluded from the assert: shooting")
print("    at the long-double bisection floor, v66/THEORY.md SS4 footnote)")
print(f"{'omega':>7} {'E_prof':>11} {'wQ_prof':>11} {'E-wQ':>9} {'G=int fp^2':>10} "
      f"{'rel.err':>9} {'E vs scan':>9}")
worst = 0.0
worst_edge = 0.0
files = sorted(glob.glob(os.path.join(RES, "profile_omega*.txt")))
scanE = {round(r["omega"], 4): r["E"] for r in rows}
for fn in files:
    m = re.search(r"omega([0-9.]+)\.txt", fn)
    w = float(m.group(1))
    rr, ff = [], []
    with open(fn) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            a, b = line.split()
            rr.append(float(a)); ff.append(float(b))
    h = rr[1] - rr[0]
    n = len(rr)
    # central-difference f'
    fp = [0.0] * n
    for i in range(1, n - 1):
        fp[i] = (ff[i + 1] - ff[i - 1]) / (2 * h)
    fp[0] = 0.0
    fp[-1] = (ff[-1] - ff[-2]) / h
    def vt(s):
        return 0.5 * MU * s / (1 + KAP * s)
    eint = [( (1.5 * (w * w * f * f + p * p + M2 * f * f) + vt(f ** 6)) * r * r)
            for r, f, p in zip(rr, ff, fp)]
    qint = [3 * w * f * f * r * r for r, f in zip(rr, ff)]
    gint = [p * p * r * r for r, p in zip(rr, fp)]
    fourpi = 4 * math.pi
    E_p = fourpi * simpson(eint, h)
    wQ_p = w * fourpi * simpson(qint, h)
    G_p = fourpi * simpson(gint, h)
    rel = abs(E_p - wQ_p - G_p) / G_p
    edge = w <= 1.3201
    if edge:
        worst_edge = max(worst_edge, rel)
    else:
        worst = max(worst, rel)
    escan = scanE.get(round(w, 4), float("nan"))
    print(f"{w:7.4f} {E_p:11.1f} {wQ_p:11.1f} {E_p-wQ_p:9.2f} {G_p:10.2f} "
          f"{rel:9.2e} {abs(E_p-escan)/escan:9.2e}" + ("  (edge)" if edge else ""))
chk("P2_identity_E_minus_wQ_eq_G_profiles", worst < 1e-3,
    f"worst rel.err {worst:.2e} over {len(files)-2} profiles "
    f"(edge rows excluded: worst {worst_edge:.2e})")

# headline numbers at the reference ball omega=1.39
r139 = [r for r in rows if abs(r["omega"] - 1.39) < 1e-9][0]
w, Q, E = r139["omega"], r139["Q"], r139["E"]
print(f"\nreference ball omega={w}: E0={E}, Q={Q}, h_eff=E/w={E/w:.1f}, "
      f"E/(wQ)={E/(w*Q):.4f}")
print(f"de Broglie: k_dB = (w/E0) p = p/h_eff = p/(Q*{E/(w*Q):.4f})")
chk("P2b_heff_eq_Q_within_4pct_at_139", abs(E / (w * Q) - 1) < 0.04,
    f"h_eff/Q = {E/(w*Q):.4f}")

# ----------------------------------------------------------------------
# 3. clock-gradient force numbers
# ----------------------------------------------------------------------
print()
print("=" * 78)
print("3. CLOCK-GRADIENT FORCE  a = -c^2 d(dw/w)/dx   [thm kinematics, estimate dw(e)]")
print("=" * 78)
dwde = lambda e: 0.0119 - 2 * 0.0032 * e
dedx = 0.2 / 25.0
for e0 in (0.0, 0.1, 0.2):
    a = -dwde(e0) * dedx / 1.39
    print(f"  e0={e0:4.1f}:  a = {a:+.3e} code units (toward LOW e_bath)")
a0 = dwde(0.0) * dedx / 1.39
aV51 = 2 * 6.0 / 400 ** 2
print(f"  V51 (FUTURE.md F25): +6 units in T=400 => a_V51 ~ 2x/T^2 = {aV51:.2e}; "
      f"avg v = {6/400:.3f}/t.u.")
print(f"  ratio a_pred/a_V51 = {a0/aV51:.2f}  (same order; DIFFERENT backgrounds)")
chk("P3_same_order_as_V51", 0.1 < a0 / aV51 < 10, f"ratio {a0/aV51:.2f}")

# ----------------------------------------------------------------------
# 4. orbit quantization estimate
# ----------------------------------------------------------------------
print()
print("=" * 78)
print("4. ORBIT QUANTIZATION (Bohr-Sommerfeld phase closure)  [estimate]")
print("   closure: D*v_rel = 2n/omega; centripetal F_req = 2M/(w^2 D^3), M=E0=692")
print("   tail force scale: F_tail ~ mu_t*U0*exp(-mu_t(D-12))*(12/D), "
      "mu_t = sqrt(m^2-w^2)")
print("=" * 78)
mu_t = math.sqrt(M2 - 1.39 ** 2)
print(f"  tail decay mu_t = {mu_t:.4f};  U0 = |U(D=12)| in [0.2, 2] code units "
      "(tb2 ejection energetics, crude)")
print(f"{'D':>4} {'v_rel(n=1)':>10} {'F_req(n=1)':>10} {'F_tail lo':>10} "
      f"{'F_tail hi':>10} {'T_orb':>8}")
Mball = 692.0
wb = 1.39
for D in (8, 10, 12, 14, 16, 20):
    vrel = 2 / (wb * D)
    Freq = 2 * Mball / (wb ** 2 * D ** 3)
    Flo = mu_t * 0.2 * math.exp(-mu_t * (D - 12)) * 12 / D
    Fhi = mu_t * 2.0 * math.exp(-mu_t * (D - 12)) * 12 / D
    Torb = math.pi * wb * D ** 2
    print(f"{D:4d} {vrel:10.3f} {Freq:10.3f} {Flo:10.3f} {Fhi:10.3f} {Torb:8.0f}")
# n=1 orbit exists where F_req crosses the [Flo,Fhi] band
ok = False
for D in [8 + 0.1 * i for i in range(121)]:
    Freq = 2 * Mball / (wb ** 2 * D ** 3)
    Flo = mu_t * 0.2 * math.exp(-mu_t * (D - 12)) * 12 / D
    Fhi = mu_t * 2.0 * math.exp(-mu_t * (D - 12)) * 12 / D
    if Flo <= Freq <= Fhi:
        ok = True
print(f"\n  n=1 orbit radius crosses the tail-force band inside D in [8,20]: {ok}")
chk("P4_n1_orbit_in_accessible_range", ok)

print()
if nfail == 0:
    print("ALL NUMERIC CHECKS PASS")
else:
    print(f"{nfail} CHECKS FAILED")
    sys.exit(1)
