#!/usr/bin/env python3
"""Compare braid separation D(t) across three simulations."""

import numpy as np

d3 = np.loadtxt("data/sep_3field.tsv", skiprows=1)
d6s = np.loadtxt("data/sep_6field_same.tsv", skiprows=1)
d6o = np.loadtxt("data/sep_6field_opp.tsv", skiprows=1)

t3, D3 = d3[:, 1], d3[:, 6]
t6s, D6s = d6s[:, 1], d6s[:, 6]
t6o, D6o = d6o[:, 1], d6o[:, 6]

# Ensure same number of frames
n = min(len(t3), len(t6s), len(t6o))
t = t3[:n]
D3 = D3[:n]
D6s = D6s[:n]
D6o = D6o[:n]

print("=" * 80)
print("Two-Braid Separation D(t): 3-field vs 6-field-same vs 6-field-opposite")
print("=" * 80)
print(f"{'t':>8s} {'D_3field':>10s} {'D_6f_same':>10s} {'D_6f_opp':>10s} "
      f"{'dD_same':>10s} {'dD_opp':>10s}")
print("-" * 58)
for i in range(n):
    dd_s = D6s[i] - D3[i]
    dd_o = D6o[i] - D3[i]
    print(f"{t[i]:8.2f} {D3[i]:10.3f} {D6s[i]:10.3f} {D6o[i]:10.3f} "
          f"{dd_s:10.3f} {dd_o:10.3f}")

# Summary statistics
print()
print("=" * 80)
print("SUMMARY")
print("=" * 80)

# Initial and final D
print(f"D(t=0):   3-field={D3[0]:.3f}  6f-same={D6s[0]:.3f}  6f-opp={D6o[0]:.3f}")
print(f"D(t=end): 3-field={D3[-1]:.3f}  6f-same={D6s[-1]:.3f}  6f-opp={D6o[-1]:.3f}")
print(f"Delta-D:  3-field={D3[-1]-D3[0]:.3f}  6f-same={D6s[-1]-D6s[0]:.3f}  6f-opp={D6o[-1]-D6o[0]:.3f}")

# Average D over second half of simulation
half = n // 2
mean_D3 = np.mean(D3[half:])
mean_D6s = np.mean(D6s[half:])
mean_D6o = np.mean(D6o[half:])
print(f"\nMean D (late time, t>{t[half]:.0f}):")
print(f"  3-field:     {mean_D3:.3f}")
print(f"  6f-same:     {mean_D6s:.3f}  (delta = {mean_D6s - mean_D3:+.3f})")
print(f"  6f-opposite: {mean_D6o:.3f}  (delta = {mean_D6o - mean_D3:+.3f})")

# Minimum D (closest approach)
min_D3 = np.min(D3)
min_D6s = np.min(D6s)
min_D6o = np.min(D6o)
print(f"\nMinimum D (closest approach):")
print(f"  3-field:     {min_D3:.3f} at t={t[np.argmin(D3)]:.1f}")
print(f"  6f-same:     {min_D6s:.3f} at t={t6s[np.argmin(D6s[:n])]:.1f}  "
      f"(delta = {min_D6s - min_D3:+.3f})")
print(f"  6f-opposite: {min_D6o:.3f} at t={t6o[np.argmin(D6o[:n])]:.1f}  "
      f"(delta = {min_D6o - min_D3:+.3f})")

# Infall rate: D(t=50) - D(t=0) / 50
mask_early = t <= 50
if np.sum(mask_early) > 2:
    v3 = np.polyfit(t[mask_early], D3[mask_early], 1)[0]
    v6s = np.polyfit(t[mask_early], D6s[mask_early], 1)[0]
    v6o = np.polyfit(t[mask_early], D6o[mask_early], 1)[0]
    print(f"\nEarly infall rate dD/dt (t<50):")
    print(f"  3-field:     {v3:.4f} /time")
    print(f"  6f-same:     {v6s:.4f} /time  ({100*(v6s-v3)/abs(v3):+.1f}%)")
    print(f"  6f-opposite: {v6o:.4f} /time  ({100*(v6o-v3)/abs(v3):+.1f}%)")

# Does theta add/change force?
print()
print("=" * 80)
print("CONCLUSIONS")
print("=" * 80)
if mean_D6s < mean_D3 - 0.1:
    print("6f-same < 3-field mean D: theta ADDS attraction (parallel currents attract)")
elif mean_D6s > mean_D3 + 0.1:
    print("6f-same > 3-field mean D: theta REDUCES attraction")
else:
    print("6f-same ~ 3-field mean D: theta has negligible effect on mean separation")

if mean_D6o > mean_D6s + 0.1:
    print("6f-opp > 6f-same mean D: opposite winding attracts LESS (antiparallel repel)")
elif mean_D6o < mean_D6s - 0.1:
    print("6f-opp < 6f-same mean D: opposite winding attracts MORE")
else:
    print("6f-opp ~ 6f-same mean D: winding direction has negligible effect")
