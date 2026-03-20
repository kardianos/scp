#!/usr/bin/env python3
"""Compare signed theta_phi between positive and negative winding simulations."""

import numpy as np

pos = np.loadtxt("data/theta_radial_pos_signed.tsv", skiprows=1)
neg = np.loadtxt("data/theta_radial_neg_signed.tsv", skiprows=1)

r_p, tp_p = pos[:, 0], pos[:, 1]  # theta_phi (signed)
r_n, tp_n = neg[:, 0], neg[:, 1]

# Match the radii (they should be the same)
n = min(len(r_p), len(r_n))
r_p, tp_p = r_p[:n], tp_p[:n]
r_n, tp_n = r_n[:n], tp_n[:n]

print("=" * 70)
print("Winding Reversal: Signed theta_phi comparison")
print("=" * 70)
print(f"{'r':>6s} {'W=+1':>12s} {'W=-1':>12s} {'sum':>12s} {'ratio':>8s}")
print("-" * 52)
for i in range(n):
    s = tp_p[i] + tp_n[i]
    ratio = tp_p[i] / tp_n[i] if abs(tp_n[i]) > 1e-20 else float('inf')
    print(f"{r_p[i]:6.2f} {tp_p[i]:12.4e} {tp_n[i]:12.4e} {s:12.4e} {ratio:8.3f}")

# Correlation
corr = np.corrcoef(tp_p, tp_n)[0, 1]
print(f"\nCorrelation coefficient: {corr:.4f}")
print(f"  (If exact sign flip: corr = -1.0)")

# Mean of |theta_phi(+) + theta_phi(-)| vs mean of |theta_phi|
mean_sum = np.mean(np.abs(tp_p + tp_n))
mean_abs = np.mean(np.abs(tp_p) + np.abs(tp_n)) / 2
print(f"\nMean |sum| = {mean_sum:.4e}")
print(f"Mean |individual| = {mean_abs:.4e}")
print(f"Cancellation ratio: {mean_sum / mean_abs:.3f}")
print(f"  (If perfect sign flip: ratio -> 0)")

# Also compare theta_phi^2 (energy) - should be similar magnitude
tp2_p = pos[:n, 2]  # theta_phi^2
tp2_n = neg[:n, 2]
ratio_energy = np.mean(tp2_p) / np.mean(tp2_n)
print(f"\nEnergy check: mean(theta_phi^2(+)) / mean(theta_phi^2(-)) = {ratio_energy:.3f}")
print(f"  (Should be ~1 if same dynamics)")

# Check at specific representative radii
print("\n--- Key radii ---")
for ri in [2.25, 4.25, 6.25, 8.25]:
    idx = np.argmin(np.abs(r_p - ri))
    print(f"r={r_p[idx]:.2f}: W=+1: {tp_p[idx]:+.4e}, W=-1: {tp_n[idx]:+.4e}, "
          f"ratio: {tp_p[idx]/tp_n[idx]:.3f}" if abs(tp_n[idx]) > 1e-20 else
          f"r={r_p[idx]:.2f}: W=+1: {tp_p[idx]:+.4e}, W=-1: {tp_n[idx]:+.4e}")
