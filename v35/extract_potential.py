#!/usr/bin/env python3
"""
V35 Phase 1: Extract effective potential V_eff(r) from V34 data.

Inputs:
  - v34/phonon_test/data/depletion_t0100.tsv  (delta_rho(r) profile)
  - v34/torsion_coupling/theta_characterize/data/theta_radial.tsv  (theta fields)

Output:
  - data/V_eff.tsv  (r, V_grav, V_em, V_total)
"""

import numpy as np
from scipy.optimize import curve_fit
import os

# ── Load V34 data ──────────────────────────────────────────────────────────

v34_base = os.path.join(os.path.dirname(__file__), '..', 'v34')

# Depletion profile (gravity)
dep_data = np.loadtxt(
    os.path.join(v34_base, 'phonon_test/data/depletion_t0100.tsv'),
    skiprows=4  # 3 comment lines + 1 header line
)
r_dep = dep_data[:, 0]
rho_dep = dep_data[:, 1]
drho_dep = dep_data[:, 2]

# Theta profile (EM)
theta_data = np.loadtxt(
    os.path.join(v34_base, 'torsion_coupling/theta_characterize/data/theta_radial.tsv'),
    skiprows=1
)
r_theta = theta_data[:, 0]
theta_phi2 = theta_data[:, 1]
theta_tot2 = theta_data[:, 4]
phi2 = theta_data[:, 5]

# ── Fit power laws ─────────────────────────────────────────────────────────

# 1. Gravity: delta_rho ~ B / r^n for r=5 to 25 (reliable range from V34 results)
mask_g = (r_dep >= 5.0) & (r_dep <= 25.0) & (drho_dep > 0)
r_gfit = r_dep[mask_g]
drho_gfit = drho_dep[mask_g]

log_r_g = np.log(r_gfit)
log_drho_g = np.log(drho_gfit)
coeffs_g = np.polyfit(log_r_g, log_drho_g, 1)
n_grav = -coeffs_g[0]
B_grav = np.exp(coeffs_g[1])

# R^2 for the gravity fit
predicted_g = np.exp(coeffs_g[1] + coeffs_g[0] * log_r_g)
ss_res_g = np.sum((drho_gfit - predicted_g)**2)
ss_tot_g = np.sum((drho_gfit - drho_gfit.mean())**2)
r2_grav = 1 - ss_res_g / ss_tot_g

print("=" * 65)
print("GRAVITY POTENTIAL: delta_rho(r) = B / r^n")
print("=" * 65)
print(f"  Fit range: r = {r_gfit[0]:.1f} to {r_gfit[-1]:.1f}")
print(f"  N points:  {len(r_gfit)}")
print(f"  B = {B_grav:.4f}")
print(f"  n = {n_grav:.4f}")
print(f"  R^2 = {r2_grav:.4f}")
print()

# The gravitational potential from delta_rho:
#   F_grav = -C * grad(rho), with C ~ 186 from V34
#   For spherical: F_grav(r) = -C * d(delta_rho)/dr
#   delta_rho ~ B/r^n  =>  F_grav ~ C*n*B / r^{n+1}
#   V_grav = -int F dr = -C*B / r^n  (up to constant)
#
# But simpler: the potential energy of a test particle in the
# density-gradient force field is proportional to delta_rho itself.
# V_grav(r) = -C_grav * delta_rho(r) = -C_grav * B_grav / r^n_grav
#
# From V34: the force between braids at D=15 gives dD/dt ~ 0.048.
# With M_braid ~ 5000 (energy units), F = M*a.
# The depletion at r=15 is delta_rho ~ 0.015.
# So C_grav = F / grad(delta_rho).
#
# For now: normalize V_grav so that V_grav(r=5) = -1.0 (at braid surface)
# This sets the energy scale. We'll scan hbar to find bound states.

V0_grav = B_grav / 5.0**n_grav  # delta_rho at r=5
print(f"  delta_rho(r=5) = {V0_grav:.4f}")
print(f"  V_grav will be normalized to V(5) = -1.0")

# 2. EM: theta force adds 27% to gravity for same-winding braids
#    From V34 theta_characterize:
#    - Same-winding: 23% faster infall => F_total = 1.23 * F_grav
#    - So F_em = 0.23 * F_grav (same spatial profile approximately)
#    Actually: the theta_phi^2(r) profile gives the spatial shape.
#    But theta is oscillatory (standing wave), so the DC force follows
#    the gravity profile shape.
#
# Conservative approach: V_em(r) = 0.27 * V_grav(r) for same winding
# (from the 27% force enhancement measured in V34)

em_fraction = 0.27  # from V34: same-winding adds 27% attraction

print()
print("=" * 65)
print("EM POTENTIAL: V_em(r) = em_fraction * V_grav(r)")
print("=" * 65)
print(f"  em_fraction = {em_fraction} (from V34 two-braid force measurement)")
print(f"  V_em has same radial shape as V_grav")
print()

# ── Construct V_eff on logarithmic grid ────────────────────────────────────

# Grid: r_min=1 to r_max=500000, N=100000 points, logarithmic
r_min = 1.0
r_max = 500000.0
N_r = 100000
r_grid = np.logspace(np.log10(r_min), np.log10(r_max), N_r)

# V_grav(r): use measured data where available, power-law extrapolation elsewhere
# For r < 5 (inside braid): repulsive core
# For r >= 5: V_grav = -B_grav / r^n_grav, normalized to V(5) = -1

def V_grav_func(r):
    """Gravitational potential from depletion profile."""
    V = np.zeros_like(r, dtype=float)

    # Outside braid (r >= 5): attractive power-law
    outer = r >= 5.0
    V[outer] = -B_grav / r[outer]**n_grav

    # Normalize so V(5) = -1.0
    V[outer] /= V0_grav

    # Inside braid (r < 5): transition to repulsive core
    # Use a smooth interpolation: V = -1 + repulsive at small r
    # Model: V(r) = -1 + C_core * (5/r - 1)^2 for r < 5
    # This gives V(5) = -1 (matches outer), V → +inf as r → 0
    inner = r < 5.0
    # Repulsive core strength: set so V(r=2) = 0 (zero crossing)
    # -1 + C*(5/2 - 1)^2 = 0 => C = 1/(1.5)^2 = 0.444
    C_core = 1.0 / (5.0/2.0 - 1.0)**2
    V[inner] = -1.0 + C_core * (5.0/r[inner] - 1.0)**2

    return V

def V_em_func(r):
    """EM potential: 27% of gravity, same spatial profile."""
    return em_fraction * V_grav_func(r)

def V_total_func(r):
    """Total potential: gravity + EM."""
    return (1.0 + em_fraction) * V_grav_func(r)

V_grav = V_grav_func(r_grid)
V_em = V_em_func(r_grid)
V_total = V_total_func(r_grid)

# ── Save to file ───────────────────────────────────────────────────────────

os.makedirs('data', exist_ok=True)
outfile = os.path.join(os.path.dirname(__file__), 'data', 'V_eff.tsv')
with open(outfile, 'w') as f:
    f.write("# V35 Phase 1: Effective potential from V34 data\n")
    f.write(f"# Gravity: delta_rho ~ {B_grav:.4f} / r^{n_grav:.4f} (R^2={r2_grav:.4f})\n")
    f.write(f"# EM fraction: {em_fraction} of gravity (same-winding)\n")
    f.write(f"# Normalized: V_grav(r=5) = -1.0\n")
    f.write(f"# Core: repulsive for r < 5 (zero crossing at r=2)\n")
    f.write("r\tV_grav\tV_em\tV_total\n")

    # Write a subset (every 100th point) for readability
    for i in range(0, N_r, 100):
        f.write(f"{r_grid[i]:.6e}\t{V_grav[i]:.6e}\t{V_em[i]:.6e}\t{V_total[i]:.6e}\n")

print(f"Saved {outfile} ({N_r//100} rows)")
print()

# ── Print key values ───────────────────────────────────────────────────────

print("=" * 65)
print("KEY POTENTIAL VALUES")
print("=" * 65)
print(f"  {'r':>10s}  {'V_grav':>12s}  {'V_em':>12s}  {'V_total':>12s}")
for rv in [1, 2, 3, 5, 10, 20, 50, 100, 500, 1000, 5000, 10000, 50000, 100000]:
    idx = np.argmin(np.abs(r_grid - rv))
    print(f"  {r_grid[idx]:10.1f}  {V_grav[idx]:12.6f}  {V_em[idx]:12.6f}  {V_total[idx]:12.6f}")

print()
print(f"Potential well depth: V_min = {V_total.min():.4f} at r = {r_grid[np.argmin(V_total)]:.2f}")
print(f"Zero crossing: r ~ 2.0 (by construction)")
print(f"Power-law exponent: n = {n_grav:.3f} (between Coulomb n=1 and Newton n=2)")
