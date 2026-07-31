#!/usr/bin/env python3
"""M-B1 balance curve: net drain B(delta0) = -dR/dt of pinned pairs.

Pre-registered (before first scoring): minimum at the rung (delta0=0);
no zero crossing anywhere (the microscopic V2g no-particle measurement);
B(0) ~ 2*c0*cap ~ 8.5e-4 Em/t.u. (the corpus per-voice bleed, x2).
"""
import glob, os, re

HERE = os.path.dirname(os.path.abspath(__file__))

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

print(f"{'delta0':>8} {'x':>7} {'B=-dR/dt [300,600]':>19} {'conv rough':>11}")
rows = []
for log in sorted(glob.glob(os.path.join(HERE, "runs", "mb1_*_pin1.log"))):
    fn = os.path.basename(log)
    m = re.search(r"mb1_d([PM])(\d+)p(\d+)_", fn)
    d0 = (1 if m.group(1) == "P" else -1) * float(m.group(2) + "." + m.group(3))
    xm = re.search(r"x0=([\d.]+)", open(log).read())
    x = float(xm.group(1)) if xm else float("nan")
    s = slope(pin_series(log), 300, 600)
    cv = re.findall(r"# CONV t=600.* rough=([\d.e+-]+)", open(log).read())
    rows.append((d0, x, -s if s is not None else float("nan")))
    print(f"{d0:>+8.3f} {x:>7.4f} {-s:>19.3e} {cv[0] if cv else '-':>11}")
rows.sort()
bmin = min(rows, key=lambda r: r[2])
print(f"\n# minimum drain at delta0={bmin[0]:+.3f} (pre-reg: 0.000), "
      f"B_min={bmin[2]:.3e} (pre-reg ~8.5e-4)")
print(f"# zero crossing: {'NO — all draining (V2g no-particle, microscopic)' if all(r[2] > 0 for r in rows) else 'YES?!'}")

ctrl = os.path.join(HERE, "runs", "mb1_dP0p00_pin0.log")
if os.path.exists(ctrl):
    s = slope(em_series(ctrl), 300, 600)
    print(f"# free-pair control (unpinned, rung): dEm/dt = {s:.3e} (whole-box dense)")
