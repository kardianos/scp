#!/usr/bin/env python3
"""M-B1/M-B4 balance curves: net drain B(delta0) = -dR/dt of pinned pairs.

Series "mb1"  (M-B1, foam V2g):     pre-registered minimum at the rung
  (delta0=0); no zero crossing (microscopic no-particle); B(0) ~
  2*c0*cap ~ 8.5e-4 Em/t.u. MEASURED: B_min=8.85e-4 at delta0=0 (+4%),
  minimum rides HEAVY of the rung.
Series "mb1s" (M-B4, annealed V2s): same instrument on quiet geometry.
  MEASURED: floor is disorder (rung 9.66e-5 = 9.2x below foam); tongue
  drain spikes at the vocabulary coincidences x=0.4167 (3:2, 40x) and
  x=0.833 (2:1, 8x); tongue-free desert x~0.44-0.78 with min 3.5e-5 =
  0.041*c0/voice near x=0.536; noble-address (x=0.515, KAM) prediction
  REFUTED — monotone through the noble point (M-B4b).
"""
import glob, os, re, sys

HERE = os.path.dirname(os.path.abspath(__file__))
C0 = 4.25e-4          # per-voice corpus bleed (V2g microscopic)
X_FIFTH = 0.4167      # (3/2-1)/q_detune — the 3:2 vacuum coincidence
X_OCT = 0.8333        # (2-1)/q_detune — the 2:1 vacuum coincidence
X_NOBLE = 0.515       # 1/phi address (LIT.md transfer #3) — refuted


def pin_series(log):
    return [(float(t), float(r)) for t, r in
            re.findall(r"# PIN t=([\d.]+) R=([-\d.e+]+)", open(log).read())]


def em_series(log):
    rows = []
    for line in open(log):
        if line.startswith("#") or "\t" not in line:
            continue
        p = line.split("\t")
        try:
            rows.append((float(p[0]), float(p[2])))
        except (ValueError, IndexError):
            pass
    return rows


def net_x(log):
    """Load address x = first vertex load of the net the run was seeded from
    (net_file recorded in the run's sibling .cfg)."""
    for src in (log[:-4] + ".cfg", log):
        if not os.path.exists(src):
            continue
        m = re.search(r"net_file=(\S+)", open(src).read())
        if m and os.path.exists(m.group(1)):
            for line in open(m.group(1)):
                if line.startswith("V "):
                    return float(line.split()[4])
    m = re.search(r"x0=([\d.]+)", open(log).read())
    return float(m.group(1)) if m else float("nan")


def slope(series, tlo, thi):
    pts = [(t, v) for t, v in series if tlo <= t <= thi]
    if len(pts) < 3:
        return None
    n = len(pts)
    mt = sum(t for t, _ in pts) / n
    mv = sum(v for _, v in pts) / n
    num = sum((t - mt) * (v - mv) for t, v in pts)
    den = sum((t - mt) ** 2 for t, _ in pts)
    return num / den if den else None


for series, label in (("mb1", "M-B1 foam V2g"), ("mb1s", "M-B4 annealed V2s")):
    logs = sorted(glob.glob(os.path.join(HERE, "runs", f"{series}_d*_pin1.log")))
    if not logs:
        continue
    print(f"== {label} ==")
    print(f"{'delta0':>8} {'x':>7} {'B=-dR/dt [300,600]':>19} {'B/2/c0':>8}")
    rows = []
    for log in logs:
        m = re.search(rf"{series}_d([PM])(\d+)p(\d+)_", os.path.basename(log))
        d0 = (1 if m.group(1) == "P" else -1) * float(m.group(2) + "." + m.group(3))
        x = net_x(log)
        s = slope(pin_series(log), 300, 600)
        B = -s if s is not None else float("nan")
        rows.append((d0, x, B))
    rows.sort()
    for d0, x, B in rows:
        print(f"{d0:>+8.4f} {x:>7.4f} {B:>19.3e} {B/2/C0:>8.3f}")
    bmin = min(rows, key=lambda r: r[2])
    print(f"# minimum drain at delta0={bmin[0]:+.4f} x={bmin[1]:.4f} "
          f"B_min={bmin[2]:.3e} ({bmin[2]/2/C0:.3f}*c0/voice)")
    print(f"# zero crossing: "
          f"{'NO — all draining (no-particle statement holds)' if all(r[2] > 0 for r in rows) else 'YES?!'}")
    if series == "mb1s":
        spikes = [(x, B) for _, x, B in rows if abs(x - X_FIFTH) < 0.01 or abs(x - X_OCT) < 0.01]
        floor = min(B for _, x, B in rows if 0.44 <= x <= 0.78)
        for x, B in spikes:
            which = "3:2 fifth" if abs(x - X_FIFTH) < 0.01 else "2:1 octave"
            print(f"# tongue spike at x={x:.4f} ({which}): B={B:.3e} = {B/floor:.0f}x desert floor")
        nbd = {}
        for _, x, B in rows:
            if 0.49 <= x <= 0.55:
                nbd.setdefault(round(x, 3), []).append(B)
        nb = sorted((x, sum(v) / len(v)) for x, v in nbd.items())
        if len(nb) >= 3:
            mono = all(nb[i][1] > nb[i+1][1] for i in range(len(nb)-1))
            print(f"# noble probe (x={X_NOBLE}): "
                  + " ".join(f"{x:.3f}:{B:.2e}" for x, B in nb)
                  + f" -> {'monotone THROUGH noble point — KAM refuted (M-B4b)' if mono else 'non-monotone?'}")
    ctrl = os.path.join(HERE, "runs", f"{series}_dP0p00_pin0.log")
    if os.path.exists(ctrl):
        s = slope(em_series(ctrl), 300, 600)
        print(f"# free-pair control (unpinned, rung): dEm/dt = {s:.3e} (whole-box dense)")
    print()
