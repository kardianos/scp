#!/usr/bin/env python3
"""M-B5 scorer: the W1 re-fight as desert riders on V2s.

Pre-registered bars (MASS.md M-B5 entry, written before launch):
  1. CURVE TRANSFER: in-situ per-voice drain -dEm/dt/N at the object's
     current x matches the M-B4 pinned-pair curve B_pv(x) within 2x.
  2. ADDRESS A/B: drain ratio m6/m5 ~ 2.3 on identical cells+topology.
  3. TONGUE TRANSIT: free1/m5/ring8 cross x=0.4167+-0.011 with a brief
     ~100x drain episode; m6 (starts below) shows none.
  4. CTRL: within 10% of m5 drain (large split = texture beyond rate).
  5. census (if reached) >= 0.5 * t_pred; t_pred: free1 33.6k,
     ring12_m5 32.2k, ring12_m6 7.7k, ring8_m3 45.5k.
"""
import math
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
CAP, X_SKIRT, X_FIFTH = 2.5, 0.0617, 0.4167
RUNS = [("free1", 1, 33610), ("ring12_m5", 12, 32180),
        ("ring12_m5_ctrl", 12, 32180), ("ring12_m6", 12, 7678),
        ("ring8_m3", 8, 45461)]

# M-B4 measured pair curve (x, B_pair) — verified 2026-08-01; spike point
# (0.4156, 4.009e-3) excluded from the smooth line (fast transit).
MB4 = [(0.2356, 1.852e-4), (0.2856, 1.437e-4), (0.3356, 9.664e-5),
       (0.3856, 7.404e-5), (0.4356, 5.252e-5), (0.4856, 4.860e-5),
       (0.4940, 4.507e-5), (0.5150, 4.002e-5), (0.5356, 3.522e-5),
       (0.6356, 4.026e-5), (0.7356, 1.113e-4), (0.7956, 3.866e-4),
       (0.8326, 3.730e-4)]
_XS = [p[0] for p in MB4]
_LB = [math.log(p[1] / 2.0) for p in MB4]


def B_pv(x):
    if x <= _XS[0]:
        s = (_LB[1] - _LB[0]) / (_XS[1] - _XS[0])
        return math.exp(_LB[0] + s * (x - _XS[0]))
    if x >= _XS[-1]:
        return math.exp(_LB[-1])
    i = max(0, min(len(_XS) - 2, next(k for k in range(len(_XS) - 1) if _XS[k + 1] > x)))
    f = (x - _XS[i]) / (_XS[i + 1] - _XS[i])
    return math.exp(_LB[i] + f * (_LB[i + 1] - _LB[i]))


def lump_series(log, nvoice):
    """(t, n_frag, M_obj, M_sys). M_obj = median-filtered sum of listed
    lump masses + E_transfer (the object's own holdings incl. its
    circulating flight; the M-B4 pin paid for GROSS shedding, so the
    free analog is the object book, not the box inventory). M_sys =
    E_dense + E_transfer (box dense inventory; declines only at the
    sink). Median filter (5) rides over cluster-listing dropouts —
    on V2s (degree 7.3) rings list as N singles with no haze bridges."""
    lumps = []
    for line in open(log):
        if line.startswith("# LUMP"):
            t = float(re.search(r"t=([\d.]+)", line).group(1))
            n = int(re.search(r"n=(\d+)", line).group(1))
            ms = [float(m) for m in re.findall(r"m=([\d.e+-]+)", line)]
            # per-voice load: total listed mass, rescaled for the top-K
            # listing cap when n exceeds it (near-uniform fragments), over
            # the seeded voice count. Robust to cluster merges (m6) and to
            # the cap (m5: n=12 lists 10).
            tot = sum(ms) * ((n / len(ms)) if ms and n > len(ms) else 1.0)
            lumps.append((t, n, (tot / nvoice) if ms else 0.0))
    diag = []
    for line in open(log):
        if line.startswith("#") or "\t" not in line:
            continue
        p = line.split("\t")
        try:
            diag.append((float(p[0]), float(p[2]), float(p[4])))
        except (ValueError, IndexError):
            pass
    out = []
    for i, (t, n, mpv) in enumerate(lumps):
        w = [lumps[j][2] for j in range(max(0, i - 2), min(len(lumps), i + 3))]
        mf = sorted(w)[len(w) // 2]
        k = min(range(len(diag)), key=lambda j: abs(diag[j][0] - t)) if diag else None
        msys = diag[k][1] + diag[k][2] if k is not None else float("nan")
        out.append((t, n, mf, msys))
    return out


def windowed_drain(series, nvoice, half=5):
    """[(t_mid, x_mid, drain_pv)] centered slopes over +-half samples."""
    out = []
    for i in range(half, len(series) - half):
        t0, _, m0, _ = series[i - half]
        t1, _, m1, _ = series[i + half]
        if t1 <= t0:
            continue
        drain = -(m1 - m0) / (t1 - t0)          # series is per-voice already
        tm, _, mm, _ = series[i]
        out.append((tm, mm / CAP, drain))
    return out


results = {}
for name, nv, tpred in RUNS:
    log = os.path.join(HERE, "runs", f"w1s_{name}.log")
    if not os.path.exists(log):
        continue
    s = lump_series(log, nv)
    if not s:
        print(f"{name}: no LUMP rows yet")
        continue
    census = max((t for t, n, _, _ in s if n > 0), default=0.0)
    alive = s[-1][1] > 0
    sysloss = s[0][3] - s[-1][3]
    wd = windowed_drain(s, nv)
    # curve-transfer ratio over the quiet segment: away from the tongue
    # coincidence (either side), above the skirt, early transient excluded
    quiet = [(t, x, d) for t, x, d in wd
             if t > 200 and abs(x - X_FIFTH) > 0.03 and x > X_SKIRT + 0.05 and d > 0]
    ratios = [d / B_pv(x) for _, x, d in quiet]
    med = sorted(ratios)[len(ratios) // 2] if ratios else float("nan")
    # tongue transit: max windowed drain inside the coincidence band vs
    # the quiet median drain
    band = [d for _, x, d in wd if abs(x - X_FIFTH) <= 0.011]
    qmed = sorted(d for _, _, d in quiet)[len(quiet) // 2] if quiet else float("nan")
    spike = (max(band) / qmed) if band and quiet and qmed > 0 else 0.0
    x_now = s[-1][2] / CAP if alive else 0.0
    results[name] = dict(census=census, alive=alive, med=med, spike=spike,
                         x_now=x_now, wd=wd, tpred=tpred)
    tail = f"ALIVE@{s[-1][0]:.0f} x={x_now:.3f}" if alive else f"census={census:.0f} ({census/tpred:.2f}*t_pred)"
    print(f"{name:16s} drain/B_pv median={med:5.2f}  tongue_spike={spike:6.1f}x  "
          f"sysloss={sysloss:+.4f}  {tail}")

print()
if "ring12_m5" in results and "ring12_m6" in results:
    # address A/B on the early window, both at their start addresses
    def early(name):
        w = [d for t, x, d in results[name]["wd"] if 500 <= t <= 2000 and d > 0]
        return sorted(w)[len(w) // 2] if w else float("nan")
    r = early("ring12_m6") / early("ring12_m5")
    print(f"# ADDRESS A/B: m6/m5 early drain = {r:.2f} (pre-reg ~2.3)")
if "ring12_m5" in results and "ring12_m5_ctrl" in results:
    a, b = results["ring12_m5"]["med"], results["ring12_m5_ctrl"]["med"]
    print(f"# CTRL: ctrl/m5 curve-ratio = {b/a:.3f} (pre-reg within 10%)")
for name in results:
    ok = 0.5 <= results[name]["med"] <= 2.0
    print(f"# BAR curve-transfer {name}: {'PASS' if ok else 'FAIL'} (median {results[name]['med']:.2f}, bar [0.5,2])")
