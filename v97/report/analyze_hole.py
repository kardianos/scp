#!/usr/bin/env python3
"""HORIZON offline meters (HORIZON.md §4).
Input: `fcsdump -mode cells run.fcs` on stdin or file arg.
Per frame about the box centre (args 2,3; default 32,32):
  pi(r) radial profile (pi = es + s_disp*(em+ee), s_disp=0.3),
  live-cell and total-cell counts inside r<6 and r<bh_r (frame infall),
  Es(r) first bins (the well), field centroid + Ee-weighted radius
  (probe tracking / trapping), tagged-blob separation when tags present
  (HZ-M matter-infall arm)."""
import sys
import numpy as np

CX = float(sys.argv[2]) if len(sys.argv) > 2 else 32.0
CY = float(sys.argv[3]) if len(sys.argv) > 3 else 32.0
SD = 0.3

def main():
    f = open(sys.argv[1]) if len(sys.argv) > 1 else sys.stdin
    f.readline()
    frames = {}
    for line in f:
        p = line.split()
        if len(p) < 14:
            continue
        frames.setdefault(float(p[0]), []).append(
            (float(p[2]), float(p[3]), float(p[6]), float(p[7]),
             float(p[8]), float(p[10])))
    print("t | pi: r<2  2-4  4-6  6-8  8-12 12-16 16-24 | ncell<4 ncell<6 "
          "nlive<6 | probe: cx cy rEe | sepAB")
    for t in sorted(frames):
        a = np.array(frames[t])
        x, y, es, em, ee, tag = (a[:, i] for i in range(6))
        dx, dy = x - CX, y - CY
        r = np.hypot(dx, dy)
        pi = es + SD * (em + ee)
        bins = [(0, 2), (2, 4), (4, 6), (6, 8), (8, 12), (12, 16), (16, 24)]
        prof = []
        for r1, r2 in bins:
            m = (r >= r1) & (r < r2)
            prof.append(pi[m].mean() if m.sum() > 3 else float('nan'))
        n4 = int((r < 4).sum())
        n6 = int((r < 6).sum())
        nl6 = int(((r < 6) & (em >= 0.05)).sum())
        lit = ee > 1e-3
        if lit.sum() > 10:
            w = ee[lit]
            pcx = np.average(x[lit], weights=w)
            pcy = np.average(y[lit], weights=w)
            pr = np.hypot(pcx - CX, pcy - CY)
            probe = f"{pcx:6.2f} {pcy:6.2f} {pr:6.2f}"
        else:
            probe = "   -      -      - "
        sep = ""
        ta, tb = tag == 1, tag == 2
        if ta.sum() > 3 and tb.sum() > 3 and (em[ta].sum() > 0.1) \
                and (em[tb].sum() > 0.1):
            ax = np.average(x[ta], weights=em[ta] + 1e-9)
            bx = np.average(x[tb], weights=em[tb] + 1e-9)
            ay = np.average(y[ta], weights=em[ta] + 1e-9)
            by = np.average(y[tb], weights=em[tb] + 1e-9)
            sep = f"{np.hypot(ax-bx, ay-by):6.3f}"
        print(f"{t:7.1f} | " + " ".join(f"{v:5.3f}" for v in prof) +
              f" | {n4:4d} {n6:4d} {nl6:4d} | {probe} | {sep}")

if __name__ == "__main__":
    main()
