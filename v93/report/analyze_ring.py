#!/usr/bin/env python3
"""v93 ring-retention meters (RESUME §7 route A).
Input: `fcsdump -mode cells run.fcs` text on stdin or as file arg.
The ring voices are the TAGGED cells (tag=1, set by the exp=ring seed).
Per frame, on the tagged cycle:
  W      discrete cycle winding of th2: sort tagged cells by azimuth about
         the box centre, sum wrapped consecutive phase increments / 2pi
         (includes the closing edge) — the topological charge.
  R_m    coherence |<e^{i(th2 - m*phi)}>| over tagged cells at the SEEDED m
         (the ring analog of item 4's R2d).
  pk     argmax_m R_m over m in [-6, 6] — where the spectrum actually peaks.
  Em_min/Em_mean   on the tagged cycle (the |psi|>0-on-cycle condition:
         winding is only defined/protected while Em_min stays > 0).
  keep   Em_tag / Em_tot — 1 minus the leak fraction into the bath.
Usage: analyze_ring.py [dump.txt] CX CY M
"""
import sys, math
import numpy as np

CX = float(sys.argv[2]) if len(sys.argv) > 2 else 8.0
CY = float(sys.argv[3]) if len(sys.argv) > 3 else 8.0
M  = int(float(sys.argv[4])) if len(sys.argv) > 4 else 2

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
            (float(p[2]), float(p[3]), float(p[7]), float(p[10]), float(p[13])))
    print("t | n_ring W R_m(seeded) pk | Em_min Em_mean keep")
    for t in sorted(frames):
        a = np.array(frames[t])
        x, y, em, tg, th = a[:, 0], a[:, 1], a[:, 2], a[:, 3], a[:, 4]
        ring = tg > 0.5
        n = int(ring.sum())
        if n < 3:
            print(f"{t:7.1f} | ring lost (n={n})")
            continue
        phi = np.arctan2(y[ring] - CY, x[ring] - CX)
        order = np.argsort(phi)
        tr = th[ring][order]
        d = np.diff(np.concatenate([tr, tr[:1]]))
        d = (d + np.pi) % (2 * np.pi) - np.pi
        W = d.sum() / (2 * np.pi)
        spec = {m: abs(np.exp(1j * (th[ring] - m * phi)).mean())
                for m in range(-6, 7)}
        Rm = spec.get(M, float('nan'))
        pk = max(spec, key=spec.get)
        emr = em[ring]
        keep = emr.sum() / em.sum() if em.sum() > 0 else float('nan')
        print(f"{t:7.1f} | {n:3d} {W:+6.3f} R{M}={Rm:.3f} pk={pk:+d} | "
              f"{emr.min():.4f} {emr.mean():.4f} {keep:.4f}")

if __name__ == "__main__":
    main()
