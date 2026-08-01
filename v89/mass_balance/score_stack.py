#!/usr/bin/env python3
"""M-B3 STACK test scorer: does contrast stack on a heavy wound ring?

Three arms, physics-identical to W1 c8_ring12 (same net, same V2g laws)
except edge_sink=0 (closed box) and the piston:
  bare     — closed box only
  es0p85   — closed box + piston pumping ambient Es toward 0.85 (deepen)
  es1p15   — closed box + piston pumping ambient Es toward 1.15 (fill)

Pre-registered (MASS.md M-B2 entry, stacking hypothesis): if the piston
contrast stacks on the ring's self-contrast, es0p85 outlives bare and
approaches the C1 range; es1p15 dies earlier (fill = destroy the pocket).
W1 absorbing-edge baseline: census death at t=449.
"""
import os, re, sys

HERE = os.path.dirname(os.path.abspath(__file__))
ARMS = ["bare", "es0p85", "es1p15"]
W1_LOG = os.path.join(HERE, "..", "prestress", "runs", "w1_c8_ring12.log")
W1_BASELINE = 449.0   # W1 scorer rule: first t with m1 < 0.3*m0 (COHESION)
LAWA = 2080.0         # from WAVE1_SCORED.tsv (x/voice 0.413)
C0 = 4.25e-4  # per-voice corpus bleed (V2g microscopic)
NVOICE = 12.0


def parse(log):
    lumps, rads, pis = [], [], []
    for line in open(log):
        if line.startswith("# LUMP"):
            t = float(re.search(r"t=([\d.]+)", line).group(1))
            n = int(re.search(r"n=(\d+)", line).group(1))
            ef = float(re.search(r"Emfree=([\d.e+-]+)", line).group(1))
            ms = [float(m) for m in re.findall(r"m=([\d.e+-]+)", line)]
            lumps.append((t, n, ef, sum(ms), max(ms) if ms else 0.0))
        elif line.startswith("# RAD"):
            t = float(re.search(r"t=([\d.]+)", line).group(1))
            sh = re.findall(r"([\d.]+):([\d.e+-]+),([\d.e+-]+),([\d.]+),([\d.e+-]+)", line)
            rads.append((t, [(float(r), float(es)) for r, _, _, es, _ in sh]))
        elif line.startswith("# PIS"):
            m = re.search(r"t=([\d.]+) es_tgt=([\d.]+) R=([\d.e+-]+)", line)
            if m:
                pis.append((float(m.group(1)), float(m.group(3))))
    return lumps, rads, pis


def linfit(pts):
    n = len(pts)
    if n < 3:
        return None
    mt = sum(t for t, _ in pts) / n
    mv = sum(v for _, v in pts) / n
    den = sum((t - mt) ** 2 for t, _ in pts)
    return sum((t - mt) * (v - mv) for t, v in pts) / den if den else None


print(f"{'arm':>12} {'cohesion':>9} {'deep0.1':>8} {'census':>7} {'/LawA':>6} "
      f"{'dM/dt[100,.8L]':>15} {'/voice/c0':>10} {'Es_min@500':>11} {'contrast':>9} {'pistR':>8}")
res = {}
for arm in ARMS + ["W1"]:
    log = W1_LOG if arm == "W1" else os.path.join(HERE, "runs", f"stack_{arm}.log")
    if not os.path.exists(log):
        continue
    lumps, rads, pis = parse(log)
    m0 = lumps[0][4] if lumps else 0.0
    t_coh = next((t for t, n, _, _, big in lumps if big < 0.3 * m0), None)
    t_deep = next((t for t, n, _, _, big in lumps if big < 0.1 * m0), None)
    life = max((t for t, n, *_ in lumps if n > 0), default=0.0)
    win = [(t, mtot) for t, n, _, mtot, _ in lumps if 100 <= t <= 0.8 * life]
    dmdt = linfit(win)
    pv = -dmdt / NVOICE / C0 if dmdt is not None else float("nan")
    esmin = contrast = float("nan")
    for t, sh in rads:
        if t >= 500:
            inner = min(es for r, es in sh if r < 6)
            outer = sum(es for r, es in sh if r > 8) / max(1, len([1 for r, _ in sh if r > 8]))
            esmin, contrast = inner, outer - inner
            break
    rend = pis[-1][1] if pis else 0.0
    res[arm] = (life, t_coh)
    print(f"{arm:>12} {t_coh or 0:>9.0f} {t_deep or 0:>8.0f} {life:>7.0f} {life/LAWA:>6.2f} "
          f"{dmdt if dmdt is not None else float('nan'):>15.3e} {pv:>10.2f} "
          f"{esmin:>11.4f} {contrast:>9.4f} {rend:>8.1f}")

print()
if "bare" in res and "es0p85" in res:
    b, p = res["bare"][0], res["es0p85"][0]
    print(f"# STACK verdict:   es0p85/bare census = {p/b:.2f} "
          f"({'STACKS' if p > 1.15 * b else 'NO STACK — self-contrast saturates'})")
if "bare" in res and "es1p15" in res:
    b, f = res["bare"][0], res["es1p15"][0]
    print(f"# OVERPRESSURE:    es1p15/bare census = {f/b:.2f} "
          f"({'LETHAL — raised ambient level doubles bleed' if f < 0.7 * b else 'no effect'})")
if "bare" in res and "W1" in res:
    print(f"# BATH (closed):   cohesion {res['W1'][1]:.0f}->{res['bare'][1]:.0f} "
          f"(+{100*(res['bare'][1]/res['W1'][1]-1):.0f}%), census {res['W1'][0]:.0f}->{res['bare'][0]:.0f} "
          f"({100*(res['bare'][0]/res['W1'][0]-1):+.0f}%) — cohesion only")
print("# NOTE: 'cohesion' = W1 scorer rule (first t with m1<0.3*m0); 'census' = last t with n>0.")
print("# The W1 449-vs-comp12-4879 anomaly is a COHESION split; census split is only 2972 vs 4879.")

# M-B3b: mass-vs-pocket calibration from the bare arm's own decay sweep.
# depth = outer-shell mean(Es) - inner-shell min(Es); fit vs M_tot for M>0.5.
bare = os.path.join(HERE, "runs", "stack_bare.log")
if os.path.exists(bare):
    lumps, rads, _ = parse(bare)
    mtot = {t: m for t, n, _, m, _ in lumps}
    pts = []
    first_depth = None
    for t, sh in rads:
        inner = min(es for r, es in sh if r < 6)
        outer = sum(es for r, es in sh if r > 8) / max(1, len([1 for r, _ in sh if r > 8]))
        if first_depth is None:
            first_depth = outer - inner
        if t in mtot and mtot[t] > 0.5:
            pts.append((mtot[t], outer - inner))
    n = len(pts)
    mx = sum(p[0] for p in pts) / n
    my = sum(p[1] for p in pts) / n
    sxx = sum((x - mx) ** 2 for x, _ in pts)
    sxy = sum((x - mx) * (y - my) for x, y in pts)
    syy = sum((y - my) ** 2 for _, y in pts)
    a, b = sxy / sxx, my - sxy / sxx * mx
    print(f"\n# M-B3b (bare arm, n={n}): depth = {a:.3e}*M + {b:.4f}  R2={sxy*sxy/(sxx*syy):.3f}")
    print(f"# intercept ~ 0 -> pocket exists only while the mass does (g4 no-steady-monopole);")
    print(f"# establishment: t=0 depth {first_depth:.3f} -> steady {max(d for _, d in pts):.3f} over ~40-100 t.u. (s_k scale)")
