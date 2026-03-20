#!/usr/bin/env python3
"""Fit power-law to theta_phi^2(r) radial profile."""

import numpy as np

data = np.loadtxt("data/theta_radial.tsv", skiprows=1)
r = data[:, 0]
theta_phi2 = data[:, 1]
theta_r2 = data[:, 2]
theta_z2 = data[:, 3]
theta_tot2 = data[:, 4]
phi2 = data[:, 5]

# Fit theta_phi^2 ~ A / r^n in the range r=2 to r=12 (outside core, inside boundary)
mask = (r >= 2.0) & (r <= 12.0)
r_fit = r[mask]
tp2_fit = theta_phi2[mask]

# Log-log linear fit
log_r = np.log(r_fit)
log_tp2 = np.log(tp2_fit)
coeffs = np.polyfit(log_r, log_tp2, 1)
n_phi = -coeffs[0]  # theta_phi^2 ~ r^{-n}
A_phi = np.exp(coeffs[1])

# Also fit total theta^2
tt2_fit = theta_tot2[mask]
log_tt2 = np.log(tt2_fit)
coeffs_tot = np.polyfit(log_r, log_tt2, 1)
n_tot = -coeffs_tot[0]

# Also fit phi^2
p2_fit = phi2[mask]
log_p2 = np.log(p2_fit)
coeffs_phi = np.polyfit(log_r, log_p2, 1)
n_p = -coeffs_phi[0]

# Fraction that is azimuthal
frac_phi = theta_phi2 / theta_tot2

print("=" * 60)
print("Power-law fits: quantity ~ r^(-n)")
print("=" * 60)
print(f"theta_phi^2:  n = {n_phi:.3f}  (A = {A_phi:.4e})")
print(f"theta_total^2: n = {n_tot:.3f}")
print(f"phi^2:        n = {n_p:.3f}")
print()
print("For Biot-Savart (B ~ 1/r): expect B^2 ~ 1/r^2, so n=2")
print(f"theta_phi^2 exponent: {n_phi:.3f}  {'CONSISTENT' if abs(n_phi - 2) < 0.5 else 'NOT consistent'} with 1/r^2")
print()

# Component fractions at key radii
print("Azimuthal fraction theta_phi^2 / theta_total^2:")
for ri in [2.25, 4.25, 6.25, 8.25, 10.25]:
    idx = np.argmin(np.abs(r - ri))
    print(f"  r={r[idx]:.2f}: {frac_phi[idx]:.3f}")

# Dominance: compare theta_phi^2 vs theta_r^2
print()
print("Dominance ratio theta_phi^2 / theta_r^2:")
for ri in [2.25, 4.25, 6.25, 8.25, 10.25]:
    idx = np.argmin(np.abs(r - ri))
    ratio = theta_phi2[idx] / theta_r2[idx] if theta_r2[idx] > 0 else float('inf')
    print(f"  r={r[idx]:.2f}: {ratio:.2f}x")

# Also fit theta_phi^2 in tighter range r=3 to r=8 (cleaner)
mask2 = (r >= 3.0) & (r <= 8.0)
r_fit2 = r[mask2]
tp2_fit2 = theta_phi2[mask2]
log_r2 = np.log(r_fit2)
log_tp22 = np.log(tp2_fit2)
coeffs2 = np.polyfit(log_r2, log_tp22, 1)
n_phi2 = -coeffs2[0]
print(f"\nTighter fit (r=3-8): theta_phi^2 exponent = {n_phi2:.3f}")

# RMS residual of fit
predicted = np.exp(coeffs[1] + coeffs[0] * log_r)
resid = np.sqrt(np.mean((log_tp2 - (coeffs[1] + coeffs[0]*log_r))**2))
print(f"Log-log RMS residual (r=2-12): {resid:.4f}")
