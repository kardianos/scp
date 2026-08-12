#!/usr/bin/env python3
"""v93 DNLS melt meters (RESUME §7 route B).
Input: `fcsdump -mode cells run.fcs` text on stdin or as file arg.
Per frame, over ALL cells (the packet leaks into untagged cells — the
tagged-centroid lesson):
  Em_tot   total dense energy (conservation cross-check)
  Em_max   peak cell occupancy — a bound packet holds its peak; a melting
           one decays it
  PR       participation number (sum Em)^2 / sum Em^2 — the number of cells
           effectively carrying the packet (self-trap holds PR small)
  rms      Em-weighted RMS radius about the Em-weighted centroid (minimum-
           image about the box centre; melt = rms growth)
  comx     Em-weighted centroid x (translation meter, all-cell — not tagged)
Usage: analyze_melt.py [dump.txt] [L]
"""
import sys, math
import numpy as np

L = float(sys.argv[2]) if len(sys.argv) > 2 else 16.0

def main():
    f = open(sys.argv[1]) if len(sys.argv) > 1 else sys.stdin
    f.readline()
    frames = {}
    for line in f:
        p = line.split()
        if len(p) < 14:
            continue
        # fcsdump cells: t i x y z r es em ee xload tag fa1 fa2 th2
        frames.setdefault(float(p[0]), []).append(
            (float(p[2]), float(p[3]), float(p[4]), float(p[7])))
    print("t | Em_tot Em_max PR rms comx")
    for t in sorted(frames):
        a = np.array(frames[t])
        x, y, z, em = a[:, 0], a[:, 1], a[:, 2], a[:, 3]
        tot = em.sum()
        if tot <= 0:
            print(f"{t:7.1f} | empty")
            continue
        pr = tot * tot / (em * em).sum()
        c = 0.5 * L
        # minimum-image displacements about the box centre, then Em-weighted
        # centroid and RMS about that centroid (packet stays near centre)
        dx = (x - c + 0.5 * L) % L - 0.5 * L
        dy = (y - c + 0.5 * L) % L - 0.5 * L
        dz = (z - c + 0.5 * L) % L - 0.5 * L
        mx, my, mz = (np.average(v, weights=em) for v in (dx, dy, dz))
        r2 = (dx - mx) ** 2 + (dy - my) ** 2 + (dz - mz) ** 2
        rms = math.sqrt(np.average(r2, weights=em))
        print(f"{t:7.1f} | {tot:9.4f} {em.max():.4f} {pr:7.1f} {rms:6.3f} {c+mx:7.3f}")

if __name__ == "__main__":
    main()
