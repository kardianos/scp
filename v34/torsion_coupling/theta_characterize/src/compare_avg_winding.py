#!/usr/bin/env python3
"""Compare time-averaged theta_phi between positive and negative winding."""

import numpy as np

pos = np.loadtxt("data/theta_avg_pos.tsv", skiprows=1)
neg = np.loadtxt("data/theta_avg_neg.tsv", skiprows=1)

r = pos[:, 0]
mean_p = pos[:, 1]
rms_p = pos[:, 2]
mean_n = neg[:, 1]
rms_n = neg[:, 2]

n = min(len(r), len(neg))
r = r[:n]
mean_p = mean_p[:n]
mean_n = mean_n[:n]

print("=" * 70)
print("TIME-AVERAGED theta_phi: W=+1 vs W=-1")
print("=" * 70)
print(f"{'r':>6s} {'<tphi>_pos':>12s} {'<tphi>_neg':>12s} {'sum':>12s} {'diff':>12s}")
print("-" * 56)
for i in range(n):
    print(f"{r[i]:6.2f} {mean_p[i]:12.4e} {mean_n[i]:12.4e} "
          f"{mean_p[i]+mean_n[i]:12.4e} {mean_p[i]-mean_n[i]:12.4e}")

# Test: if sign flips, mean_p + mean_n -> 0, mean_p - mean_n -> 2*mean_p
sum_signal = np.mean(np.abs(mean_p + mean_n))
diff_signal = np.mean(np.abs(mean_p - mean_n))
indiv_signal = np.mean(np.abs(mean_p) + np.abs(mean_n)) / 2

print(f"\nMean |sum(pos+neg)|  = {sum_signal:.4e}")
print(f"Mean |diff(pos-neg)| = {diff_signal:.4e}")
print(f"Mean |individual|    = {indiv_signal:.4e}")
print(f"Cancellation (sum/indiv): {sum_signal/indiv_signal:.3f}")
print(f"  If sign flip: cancellation -> 0, diff/indiv -> 2")
print(f"  If same sign: cancellation -> 2, diff/indiv -> 0")

# But the real question: is the DC part significant at all?
mean_ratio_pos = np.mean(pos[:n, 3])  # |mean|/rms
mean_ratio_neg = np.mean(neg[:n, 3])
print(f"\nMean |<tphi>|/rms_pos = {mean_ratio_pos:.5f}")
print(f"Mean |<tphi>|/rms_neg = {mean_ratio_neg:.5f}")
print(f"  (if >> 0.01: DC component significant; if << 0.01: pure oscillation)")

# Correlation of DC residuals
corr = np.corrcoef(mean_p, mean_n)[0, 1]
print(f"\nCorrelation of DC residuals: {corr:.4f}")
print(f"  -1.0 = perfect sign flip")
print(f"  +1.0 = same sign")
print(f"   0.0 = uncorrelated (noise)")

# RMS comparison (energy, should be equal)
rms_ratio = np.mean(rms_p[:n]) / np.mean(rms_n[:n])
print(f"\nRMS ratio (pos/neg): {rms_ratio:.3f}")
print(f"  (should be ~1 for same dynamics)")
