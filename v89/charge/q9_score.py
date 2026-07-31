#!/usr/bin/env python3
"""Q9 scorer: parse # RESULT slip rows -> I-V table + fits + plot.

Bars (pre-registered, CHARGE.md §7.5):
  B1 locked branch exists: nu ~ 0 below a critical detune dw_c
  B2 onset: (2*pi*nu)^2 linear in dw^2 (saddle-node / Adler form)
  B3 asymptote: nu -> dw_bare/2pi (the Josephson line, torque-free)
All axes use the MEASURED bare detune dw_bare (flight-load corrected),
never the naive seeded value.
"""
import glob, math, os, re

HERE = os.path.dirname(os.path.abspath(__file__))
RUNS = os.path.join(HERE, "runs")

rows = []
for log in sorted(glob.glob(os.path.join(RUNS, "q9_dw*.log"))):
    txt = open(log).read()
    m = re.search(r"# RESULT slip e=0 nu_slip=([+-][\d.]+) dw_bare=([+-][\d.]+) "
                  r"dw_over_2pi=([+-][\d.]+) slips=([+-][\d.]+) locked=(\d)", txt)
    if not m:
        print(f"!! no result in {log}")
        continue
    nu, dwb, dw2pi, slips, locked = (float(m.group(1)), float(m.group(2)),
                                     float(m.group(3)), float(m.group(4)),
                                     int(m.group(5)))
    seed_dw = float(re.search(r"dw=([\d.]+)", open(
        os.path.join(HERE, "nets", os.path.basename(log)[:-4] + ".net")).read()).group(1))
    rows.append((seed_dw, dwb, nu, dw2pi, slips, locked))

rows.sort()
print(f"{'dw_seed':>8} {'dw_meas':>9} {'nu_slip':>10} {'dw/2pi':>9} "
      f"{'nu/(dw/2pi)':>11} {'slips':>7} {'locked':>6}")
for sd, dwb, nu, dw2pi, slips, locked in rows:
    ratio = nu / dw2pi if abs(dw2pi) > 1e-9 else float('nan')
    print(f"{sd:>8.3f} {dwb:>9.4f} {nu:>10.6f} {dw2pi:>9.6f} "
          f"{ratio:>11.3f} {slips:>7.1f} {locked:>6d}")

# fits
run = [(dwb, nu) for _, dwb, nu, _, _, lk in rows if not lk and nu > 0]
lock = [dwb for _, dwb, nu, _, _, lk in rows if lk]
if run:
    # Adler form: (2*pi*nu)^2 = a*dw^2 + b  ->  slope a (Josephson: 1), dw_c = sqrt(-b/a)
    xs = [d * d for d, _ in run]
    ys = [(2 * math.pi * n) ** 2 for _, n in run]
    n = len(xs)
    sx, sy = sum(xs), sum(ys)
    sxx = sum(x * x for x in xs)
    sxy = sum(x * y for x, y in zip(xs, ys))
    a = (n * sxy - sx * sy) / (n * sxx - sx * sx) if n > 1 else float('nan')
    b = (sy - a * sx) / n if n > 1 else float('nan')
    dwc = math.sqrt(-b / a) if n > 1 and a > 0 and b < 0 else float('nan')
    ym = sy / n
    r2 = 1 - sum((y - (a * x + b)) ** 2 for x, y in zip(xs, ys)) / \
        max(1e-30, sum((y - ym) ** 2 for y in ys))
    print(f"\n# Adler fit (running branch, n={n}): (2pi*nu)^2 = a*dw^2 + b")
    print(f"#   a = {a:.4f}  (Josephson asymptote bar: a = 1)")
    print(f"#   dw_c = {dwc:.4f}  (locked branch extends to {max(lock) if lock else 0:.4f})")
    print(f"#   R^2 = {r2:.4f}")

# plot
try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    fig, ax = plt.subplots(figsize=(7, 5))
    d_run = [d for d, _ in run]
    n_run = [n_ for _, n_ in run]
    d_lk = lock
    ax.plot(d_run, n_run, "o", ms=7, label="measured (running)", zorder=5)
    ax.plot(d_lk, [0] * len(d_lk), "s", ms=7, label="measured (locked)", zorder=5)
    dd = np.linspace(0, max(d_run) * 1.05 if d_run else 0.6, 300)
    ax.plot(dd, dd / (2 * np.pi), "--", lw=1, label=r"Josephson line $\nu=\Delta\omega/2\pi$")
    if run and not math.isnan(dwc):
        da = dd[dd > dwc]
        ax.plot(da, np.sqrt(a * da**2 + b) / (2 * np.pi), "-", lw=1,
                label=fr"Adler fit  $\Delta\omega_c$={dwc:.3f}")
    ax.set_xlabel(r"measured bare detune  $\Delta\omega$  (pitch difference = 'voltage')")
    ax.set_ylabel(r"slip rate  $\nu$  (holonomy/t.u. = 'current')")
    ax.set_title("Q9 — the v89 weak link I–V  (laws_V2g, no new constants)")
    ax.legend()
    ax.grid(alpha=0.3)
    fig.tight_layout()
    out = os.path.join(RUNS, "q9_iv.png")
    fig.savefig(out, dpi=140)
    print(f"# plot: {out}")
except Exception as e:
    print(f"# plot skipped: {e}")
