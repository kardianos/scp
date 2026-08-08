#!/usr/bin/env python3
"""v93 hoard spectrum (FACE C). Input: `fcsdump -mode cells run.fcs` final
frame on stdin/file. Prints active-cell count, Em percentiles, and a coarse
histogram of the upper tail (the condensed hoards) — is there a
characteristic hoard mass?"""
import sys, numpy as np
f = open(sys.argv[1]) if len(sys.argv) > 1 and sys.argv[1] != "-" else sys.stdin
f.readline()
frames = {}
for line in f:
    p = line.split()
    if len(p) < 14:
        continue
    frames.setdefault(float(p[0]), []).append(float(p[7]))
tmax = max(frames) if frames else None
em = np.array(frames[tmax]) if frames else np.array([])
print(f"(last frame t={tmax})")
act = em[em > 0.05]
print(f"n_active(>0.05)={len(act)} Em_tot={em.sum():.3f} Em_max={em.max():.4f}")
if len(act):
    for q in [50, 75, 90, 95, 99]:
        print(f"  p{q}={np.percentile(act,q):.4f}")
    hoards = act[act > 0.5]
    print(f"  n_hoards(>0.5)={len(hoards)} hoard_Em_tot={hoards.sum():.3f}")
    if len(hoards):
        hi = max(hoards.max(), 2.0)
        h, e = np.histogram(hoards, bins=12, range=(0.5, hi))
        for i, c in enumerate(h):
            if c:
                print(f"    [{e[i]:.2f},{e[i+1]:.2f}): {c}")
