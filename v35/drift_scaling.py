#!/usr/bin/env python3
"""V35: Energy drift scaling with epsilon."""
import numpy as np

L = 25.0
N = 80
dx = 2.0 * L / (N - 1)

epsilons = [0.0001, 0.001, 0.005, 0.01, 0.05, 0.1]
drifts_pct = []  # total drift at t=50

for eps in epsilons:
    tag = str(eps)
    path = f"data/quant_eps_{tag}_timeseries.tsv"
    data = np.loadtxt(path, skiprows=1)
    E = data[:, 9]
    drift = 100.0 * (E[-1] - E[0]) / abs(E[0])
    drifts_pct.append(drift)

# Fit log(drift) vs log(eps)
log_eps = np.log10(epsilons)
log_drift = np.log10(np.abs(drifts_pct))
p = np.polyfit(log_eps, log_drift, 1)

print("Energy drift scaling: |dE/E| ~ eps^alpha")
print(f"  alpha = {p[0]:.2f}")
print(f"  (best fit: log10|dE/E| = {p[0]:.2f} * log10(eps) + {p[1]:.2f})")
print()
print("Individual points:")
for i, eps in enumerate(epsilons):
    hbar = eps * dx
    print(f"  eps={eps:8.4f}  hbar={hbar:.4e}  drift={drifts_pct[i]:+8.1f}%")

# Also check: drift rate in energy units per time step
print()
print("Drift rate (code energy per time unit):")
for eps in epsilons:
    tag = str(eps)
    data = np.loadtxt(f"data/quant_eps_{tag}_timeseries.tsv", skiprows=1)
    t = data[:, 0]
    E = data[:, 9]
    # Rate from second half
    n2 = len(t) // 2
    rate = np.polyfit(t[n2:], E[n2:], 1)[0]
    # Energy injected per rounding step
    # Each step: N^3 * 6 fields * (eps/2)^2 / 12 average rounding error
    # Variance of rounding error = eps^2/12 per value
    # Number of values = N^3 * 6
    N3 = N**3
    dt_val = 0.1 * dx
    var_per_step = N3 * 6 * (eps**2 / 12.0)
    rate_expected = var_per_step / dt_val  # rough estimate
    print(f"  eps={eps:8.4f}  dE/dt={rate:+10.2f}  expected~{rate_expected:.2f}")
