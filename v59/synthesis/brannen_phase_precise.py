#!/usr/bin/env python3
"""
Test structural Brannen phase conjectures:
  φ_l = 2/9     (lepton, structural; was already established)
  φ_d ≈ 1/9     (CONJECTURE — half the lepton phase, found by Maxima at 1%)
  φ_u ≈ -1/14   (CONJECTURE — = -1/dim G_2, found at 1%)

Do a HIGH-PRECISION numerical fit of Brannen phase per sector with v59 t² fixed,
then compare to the conjectured structural values.
"""
import math
import numpy as np
from scipy.optimize import brentq, minimize_scalar

# Empirical masses (PDG 2024)
m_e   = 0.51099895069
m_mu  = 105.6583755
m_tau = 1776.86
m_u   = 2.16
m_d   = 4.67
m_s   = 93.4
m_c   = 1273.0
m_t   = 172570.0
m_b   = 4180.0

def brannen_masses(a, t, phi):
    return sorted([(a*(1+2*t*math.cos(phi+2*math.pi*k/3)))**2 for k in range(3)])

def fit_phi(masses, t_sq):
    """Find phi minimizing relative-mass deviation from Brannen formula."""
    t = math.sqrt(t_sq)
    a = sum(math.sqrt(m) for m in masses)/3.0
    obs = sorted(masses)
    def err(phi):
        pred = brannen_masses(a, t, phi)
        return sum(((p-o)/o)**2 for p, o in zip(pred, obs))
    # Scan many phi values; the minimum is unique up to Z_3
    best_phi = 0
    best_err = 1e20
    for trial in np.linspace(-math.pi, math.pi, 100000):
        e = err(trial)
        if e < best_err:
            best_err = e
            best_phi = trial
    # Refine
    from scipy.optimize import minimize
    res = minimize(lambda p: err(p[0]), [best_phi], method='Nelder-Mead',
                   options={'xatol':1e-10, 'fatol':1e-15})
    return res.x[0], res.fun

# Fit each sector
print("Sector  | empirical t²  | best-fit φ_X (rad)  | relative err²")
print("-"*72)
for sec, m_list, t_sq in [
    ('lepton', [m_e, m_mu, m_tau], 0.5),
    ('d-quark', [m_d, m_s, m_b], 0.6),
    ('u-quark', [m_u, m_c, m_t], 7/9),
]:
    phi_fit, err = fit_phi(m_list, t_sq)
    # Reduce mod Z_3
    phi_z3 = phi_fit % (2*math.pi/3)
    if phi_z3 > math.pi/3:
        phi_z3 -= 2*math.pi/3
    if phi_z3 < -math.pi/3:
        phi_z3 += 2*math.pi/3
    print(f"{sec:8s}|  {t_sq:11.5f}  |  {phi_fit:+.7f}  |  {err:.3e}")
    print(f"           Z_3-reduced phase = {phi_z3:+.7f}")

print()
print("=" * 72)
print("Comparison to v59 conjectures")
print("=" * 72)

structural_candidates = [
    ("φ_l = 2/9",       2/9),
    ("φ_l = Q_l/3",      (2/3)/3),
    ("φ_d = 1/9",        1/9),
    ("φ_d = φ_l/2",      (2/9)/2),
    ("φ_d = (1-Q_d)/(...)", 4/15/3),     # (1-Q_d)/3 = (4/15)/3
    ("φ_d = 11/100",     11/100),
    ("φ_d = (1-Q_d)/2.5", 4/15/2.5),
    ("φ_u = -1/14",      -1/14),
    ("φ_u = -1/dimG2",   -1/14),
    ("φ_u = -1/16",      -1/16),         # -1/Cl(3,1)
    ("φ_u = -2/27",      -2/27),         # -2/gen³
    ("φ_u = -π/(2·63)",  -math.pi/126),  # related to D_u
    ("φ_u = (1-Q_u)/(-4)", (1-23/27)/(-4)),
]

for name, val in structural_candidates:
    print(f"  {name:30s} = {val:+.7f}")

# Now fit each sector and compare
print()
print("Best-fit phases vs conjectured structural values:")
print("-" * 72)
phi_l_fit, _ = fit_phi([m_e, m_mu, m_tau], 0.5)
phi_d_fit, _ = fit_phi([m_d, m_s, m_b], 0.6)
phi_u_fit, _ = fit_phi([m_u, m_c, m_t], 7/9)

# Reduce mod Z_3 and pick the |φ| < π/3 representative
def reduce_z3(phi):
    while phi > math.pi/3:
        phi -= 2*math.pi/3
    while phi < -math.pi/3:
        phi += 2*math.pi/3
    return phi

phi_l_z3 = reduce_z3(phi_l_fit)
phi_d_z3 = reduce_z3(phi_d_fit)
phi_u_z3 = reduce_z3(phi_u_fit)

print(f"\nlepton best-fit (Z_3-reduced): {phi_l_z3:+.7f}")
print(f"  vs 2/9        = {2/9:+.7f}     gap {abs(phi_l_z3 - 2/9)/(2/9)*100:.4f}%")
print(f"  vs (2/9)·(-1) = {-2/9:+.7f}     gap {abs(phi_l_z3 + 2/9)/(2/9)*100:.4f}%")

print(f"\nd-quark best-fit (Z_3-reduced): {phi_d_z3:+.7f}")
print(f"  vs 1/9        = {1/9:+.7f}     gap {abs(phi_d_z3 - 1/9)/(1/9)*100:.4f}%")
print(f"  vs -1/9       = {-1/9:+.7f}     gap {abs(phi_d_z3 + 1/9)/(1/9)*100:.4f}%")
print(f"  vs (1-Q_d)/4  = {(4/15)/4:+.7f}     gap {abs(phi_d_z3 - 1/15)/((1/15))*100:.4f}%")

print(f"\nu-quark best-fit (Z_3-reduced): {phi_u_z3:+.7f}")
print(f"  vs -1/14         = {-1/14:+.7f}     gap {abs(phi_u_z3 + 1/14)/(1/14)*100:.4f}%")
print(f"  vs +1/14         = {1/14:+.7f}      gap {abs(phi_u_z3 - 1/14)/(1/14)*100:.4f}%")
print(f"  vs -2/27         = {-2/27:+.7f}     gap {abs(phi_u_z3 + 2/27)/(2/27)*100:.4f}%")
print(f"  vs -π/63·2       = {-math.pi/(63):+.7f}     gap {abs(phi_u_z3 + math.pi/63)/(math.pi/63)*100:.4f}%")

print()
print("=" * 72)
print("KEY OBSERVATIONS")
print("=" * 72)
print(f"""
Lepton (high precision): φ_l = {phi_l_z3:+.7f} rad
  Conjecture φ_l = 2/9 = 0.2222... is EXACT to 6+ significant digits.

d-quark: φ_d = {phi_d_z3:+.7f} rad
  Conjecture φ_d = 1/9 = 0.1111... :  match at {abs(phi_d_z3 - 1/9)/(1/9)*100:.3f}%

u-quark: φ_u = {phi_u_z3:+.7f} rad
  Conjecture φ_u = -1/14 = -0.0714... :  match at {abs(phi_u_z3 + 1/14)/(1/14)*100:.3f}%

If the pattern is:
  φ_l = 2/9
  φ_d = 1/9
  φ_u = -1/14

Then there's NO uniform scheme — the d-quark uses 9 (= gen²) but u-quark uses 14 (= dim G_2).

OR alternatively, look for a uniform pattern in terms of v59 integers:
  φ_l = +2/9   (Q_l num / 9 = 2/9; or dim G_2/(3·dim Spin(7)) = 14/63 = 2/9)
  φ_d = +1/9   (Q_l num / (2·9) = 1/9; or factor 1/2 of lepton)
  φ_u = -1/14  (= -1/dim G_2; matches at 1%)
""")

# Sum of phases
print(f"\nSum of fitted phases (mod 2π): {phi_l_z3 + phi_d_z3 + phi_u_z3:.6f}")
print(f"Conjecture-sum (2/9 + 1/9 - 1/14): {2/9 + 1/9 - 1/14:.6f}")
print(f"Cabibbo θ_C: {math.asin(0.2253):.6f}")
print(f"Sin θ_C: {0.2253:.6f}")
print()
print(f"Difference (φ_d - φ_u) under conjecture: 1/9 - (-1/14) = {1/9 + 1/14:.6f} = 23/126")
print(f"  Match with Cabibbo θ_C = 0.227? {abs(1/9 + 1/14 - 0.227)/0.227*100:.3f}%")
print(f"  23 = u-quark Koide numerator, 126 = 9·14 = gen²·dim G_2")
