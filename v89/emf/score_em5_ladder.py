#!/usr/bin/env python3
"""EM5 in-kernel dispersion gate: cone fit on the box ladder.

Reads em5_lad_n*.log (kernel runs with em5=1, em5_seed_n=n), collects
# RESULT em5_disp rows, fits omega(k) for k <= 1.05 with the
pre-registered bar |b| <= 0.1*m and R2 >= 0.99 (the maxfab bar).
Prototype reference (emf/runs/lg_ladder_foam.log): omega = 0.672*k + 0.033.
Also checks Hdrift (staggered-invariant book) and sum_err per run.
usage: score_em5_ladder.py <logdir>
"""
import glob
import os
import re
import sys

logdir = sys.argv[1] if len(sys.argv) > 1 else "."
rows = []
for log in sorted(glob.glob(os.path.join(logdir, "em5_lad_n*.log"))):
    txt = open(log).read()
    m = re.search(r"# RESULT em5_disp n=(\d+) kx=([\d.]+) omega_peak=([\d.]+)", txt)
    h = re.search(r"# RESULT em5 NF=\d+ Hdrift=([-\d.e+]+) sum_err=(\d+)", txt)
    if m:
        rows.append((float(m.group(2)), float(m.group(3)),
                     float(h.group(1)) if h else float("nan"),
                     int(h.group(2)) if h else -1))
        print(f"n={m.group(1):>2} kx={m.group(2)} omega={m.group(3)} "
              f"Hdrift={h.group(1) if h else '?'} sum_err={h.group(2) if h else '?'}")
if len(rows) < 3:
    sys.exit("need >=3 ladder points")
ks = [k for k, w, _, _ in rows if k <= 1.05]
ws = [w for k, w, _, _ in rows if k <= 1.05]
n = len(ks)
mk = sum(ks) / n
mw = sum(ws) / n
sxx = sum((k - mk) ** 2 for k in ks)
sxy = sum((k - mk) * (w - mw) for k, w in zip(ks, ws))
m_, b_ = sxy / sxx, mw - sxy / sxx * mk
sst = sum((w - mw) ** 2 for w in ws)
sse = sum((w - (m_ * k + b_)) ** 2 for k, w in zip(ks, ws))
r2 = 1 - sse / sst if sst > 0 else 1.0
ok = abs(b_) <= 0.1 * abs(m_) and r2 >= 0.99
print(f"# CONE fit (k<=1): omega = {m_:.3f}*k {b_:+.3f}, R2={r2:.4f} "
      f"-> {'CONE — EM1-flat gate instrument LIVE in-kernel' if ok else 'BAR FAILED'}")
print(f"# prototype reference: 0.672*k +0.033 R2=0.9997 (lg_ladder_foam)")
hmax = max(abs(h) for _, _, h, _ in rows)
smax = max(s for _, _, _, s in rows)
print(f"# books: max|Hdrift|={hmax:.2e} max sum_err={smax} "
      f"({'staggered-invariant class holds' if hmax < 1e-10 else 'DRIFT — investigate'})")
