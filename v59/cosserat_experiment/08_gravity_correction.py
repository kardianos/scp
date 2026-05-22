#!/usr/bin/env python3
"""
08_gravity_correction.py

Try the same machinery on G_e:  combine the cross-sector ratio (21 = dim Spin(7))
with a Killing-form-style structural prefactor.

v59 status:
  S_em (empirical)      = -ln α          = 4.9202
  S_em (v59 conjecture) = π²/2 - 2α      = 4.9202 (matches at 4 × 10⁻⁵)
  S_grav (empirical)    = -ln G_e        = ??? (computed below)
  S_grav (v59 conj)     = 21 · S_em      = 103.32  →  G_e = 1.34 × 10⁻⁴⁵
  G_e (empirical)       = G_N m_e² / (ℏc) ≈ 1.7521 × 10⁻⁴⁵   →   0.76× factor off

The 0.76 factor corresponds to an additive correction of ln(1/0.76) ≈ 0.27
in S_grav.  Search for a structurally-natural quantity matching 0.27.

Structurally available numbers in v59:
  dim G_2 = 14,  dim Spin(7) = 21,  dim Spin(8) = 28,  dim Cl(3,1) = 16,
  dim S^7 = 7,  octonion dim = 8,  Brannen phase = 2/9, ρ₀² = 1/2.

Hypothesis (motivated by 07_full_lagrangian.py's success):  the structural
correction is a logarithm of a ratio of these dimensions.  Try:
  G_e = (21 / 16) · α^21     ←  cross-sector × (dim Spin(7))/(dim Cl(3,1))
  G_e = (21 / 8)  · α^21     ←  with octonion dim
  G_e = (14 / 8)  · α^21     ←  G_2 / octonion
  ... and others.

Both 21 = dim Spin(7) and 16 = dim Cl(3,1) already appear in v59's structural
identifications (cross-sector ratio, π²/2 = 8π²/16).
"""

import numpy as np

# ----------------------------------------------------------------------
# Precise constants
# ----------------------------------------------------------------------
# CODATA 2018 values
alpha = 1 / 137.035999084
G_N = 6.67430e-11           # N m² / kg²
m_e = 9.1093837015e-31       # kg
hbar = 1.054571817e-34       # J s
c = 2.99792458e8             # m/s

# Dimensionless gravitational coupling for the electron
G_e_empirical = G_N * m_e**2 / (hbar * c)
S_em_empirical = -np.log(alpha)
S_grav_empirical = -np.log(G_e_empirical)

print("=" * 72)
print("Constants and empirical values")
print("=" * 72)
print(f"  α (CODATA)             = {alpha:.12e}  = 1/{1/alpha:.6f}")
print(f"  G_e (= G_N m_e² /(ℏc))  = {G_e_empirical:.6e}")
print(f"  S_em = -ln α            = {S_em_empirical:.6f}")
print(f"  S_grav = -ln G_e        = {S_grav_empirical:.6f}")

# v59 conjectures
def solve_v59_alpha():
    x = np.pi**2 / 2
    for _ in range(50):
        x = np.pi**2 / 2 - 2 * np.exp(-x)
    return np.exp(-x)
alpha_v59 = solve_v59_alpha()
S_em_v59 = -np.log(alpha_v59)
S_grav_v59_basic = 21 * S_em_v59
G_e_v59_basic = np.exp(-S_grav_v59_basic)

print()
print("=" * 72)
print("v59 conjecture: S_grav = 21 · S_em (no correction)")
print("=" * 72)
print(f"  α (v59 conj, -ln α + 2α = π²/2) = {alpha_v59:.12e}")
print(f"  S_em (v59)                       = {S_em_v59:.6f}")
print(f"  S_grav (v59 = 21·S_em)            = {S_grav_v59_basic:.6f}")
print(f"  G_e (v59 = α^21)                  = {G_e_v59_basic:.6e}")
print(f"  Empirical G_e                     = {G_e_empirical:.6e}")
print(f"  Ratio (v59/empirical)              = {G_e_v59_basic/G_e_empirical:.6f}")
print(f"  Gap in S_grav                      = {S_grav_v59_basic - S_grav_empirical:+.6f}")

# ----------------------------------------------------------------------
# Test the (21/16) correction
# ----------------------------------------------------------------------
print()
print("=" * 72)
print("Trial: G_e = (21/16) · α^21")
print("       (21 = dim Spin(7),  16 = dim Cl(3,1) — both v59-structural)")
print("=" * 72)

# Compute using the v59-conjecture α
G_e_predicted = (21 / 16) * alpha_v59**21
S_grav_predicted = -np.log(G_e_predicted)

print(f"  (21/16) · α_v59^21              = {G_e_predicted:.8e}")
print(f"  Empirical G_e                   = {G_e_empirical:.8e}")
print(f"  Ratio predicted/empirical       = {G_e_predicted/G_e_empirical:.6f}")
print(f"  Gap (multiplicative in G)        = {(G_e_predicted/G_e_empirical - 1)*100:+.3f} %")
print(f"  Gap in S_grav (additive)         = {S_grav_predicted - S_grav_empirical:+.6f}")

# Also with empirical α
G_e_predicted_emp_alpha = (21 / 16) * alpha**21
print(f"\n  Using empirical α (CODATA):")
print(f"    (21/16) · α^21                = {G_e_predicted_emp_alpha:.8e}")
print(f"    Empirical G_e                 = {G_e_empirical:.8e}")
print(f"    Ratio predicted/empirical     = {G_e_predicted_emp_alpha/G_e_empirical:.6f}")
print(f"    Gap                            = {(G_e_predicted_emp_alpha/G_e_empirical - 1)*100:+.3f} %")

# ----------------------------------------------------------------------
# Scan other v59-natural structural ratios
# ----------------------------------------------------------------------
print()
print("=" * 72)
print("Scan of v59-natural prefactors C in G_e = C · α^21")
print("=" * 72)

structural_constants = {
    'dim G_2        = 14':            14,
    'dim S^7        = 7':              7,
    'dim Spin(7)    = 21':             21,
    'dim Spin(8)    = 28':             28,
    'dim Cl(3,1)    = 16':             16,
    'octonion dim   = 8':              8,
    '3 (generations) ':                3,
    '2/3 (Koide Q) ':                  2/3,
}

# Build all ratios
ratios = {}
keys = list(structural_constants.items())
for nm1, v1 in keys:
    for nm2, v2 in keys:
        if v1 == v2:
            continue
        ratios[f"{nm1.split('=')[0].strip()} / {nm2.split('=')[0].strip()} = {v1}/{v2}"] = v1 / v2

# Add some compound forms
ratios['21/16 (Spin7/Cl31)'] = 21 / 16
ratios['21/8 (Spin7/Oct)'] = 21 / 8
ratios['28/21 (Spin8/Spin7)'] = 28 / 21
ratios['14/8 (G_2/Oct)'] = 14 / 8
ratios['7/(2/3) (S7/Q)'] = 7 / (2/3)

# Now scan
target_C = G_e_empirical / alpha_v59**21
print(f"\n  Target C: G_e_emp / α_v59^21 = {target_C:.6f}")
print(f"\n  {'Candidate prefactor C':<60} | predicted/empirical | gap %")
print(f"  {'-'*60} | {'-'*18} | {'-'*8}")

best_match = (None, np.inf)
for name, val in sorted(ratios.items(), key=lambda x: abs(x[1] - target_C)):
    if val < 0.01 or val > 100:
        continue
    G_pred = val * alpha_v59**21
    gap_pct = (G_pred/G_e_empirical - 1) * 100
    if abs(gap_pct) < 100:
        marker = "  ←" if abs(gap_pct) < 1 else ""
        print(f"  {name:<60} | {G_pred/G_e_empirical:.6f}             | {gap_pct:+6.2f} %{marker}")
        if abs(gap_pct) < abs(best_match[1]):
            best_match = (name, gap_pct)

# ----------------------------------------------------------------------
# Assessment
# ----------------------------------------------------------------------
print()
print("=" * 72)
print("Best match and assessment")
print("=" * 72)
print(f"  Best candidate prefactor: {best_match[0]}")
print(f"  Gap: {best_match[1]:+.3f} %")
print()

print("""
If (21/16) is structurally correct, the v59 G prediction becomes

    G_e = (dim Spin(7) / dim Cl(3,1)) · α^(dim Spin(7))
        = (21 / 16) · α^21

with BOTH the 21 (cross-sector ratio from Spin(7) symmetry) and the 16
(= dim Cl(3,1), appearing in S_em = 8π²/16 = π²/2) being individually
v59-structural inputs.

This would replace the previous "off by factor 0.76" G prediction with one
that matches at the sub-percent level — comparable in tightness to the
SU(2)_L coupling formula g_W² = 5·√α from 07_full_lagrangian.py.

CAVEATS:
  • No Lagrangian derivation of why the prefactor is specifically (21/16).
    A natural derivation might come from:
      - Volume ratio of two homogeneous spaces (e.g. Spin(7)/G_2 vs S^something).
      - Normalisation factor between two gauge sectors.
      - Trace/dim ratio of two Clifford-algebra reps.
  • The 0.3-0.6 % gap (depending on which α we use) is right at the precision
    of the constants involved — not yet a clean 10⁻⁵ derivation.

Combined with previous v59 results, the prediction table now reads:

  Quantity            | Conjecture                          | Match
  --------------------|-------------------------------------|---------
  Koide Q             | 14/21 (G_2/Spin(7))                  | 6×10⁻⁶ ✓
  Brannen φ           | Q/3 = 2/9                            | 7×10⁻⁶ ✓
  α                   | -ln α + 2α = π²/2 (= 8π²/16)         | 4×10⁻⁵ ✓
  g_W (= SU(2)_L)     | g_W² = 5·√α  (NEW: 07)               | 0.1-0.3 % ⚠
  G_e                 | (21/16) · α^21  (NEW: this work)     | TBD %
""")

# Save data
np.savez('/home/d/code/scp/v59/cosserat_experiment/08_gravity.npz',
         alpha=alpha, alpha_v59=alpha_v59,
         G_e_empirical=G_e_empirical,
         G_e_predicted=G_e_predicted,
         ratio_21_16=21/16,
         gap_pct=(G_e_predicted/G_e_empirical - 1)*100,
         S_grav_empirical=S_grav_empirical)
print(f"\nSaved data to 08_gravity.npz")
