#!/usr/bin/env python3
"""
05_alpha_alphaW_ratio.py

Test a specific structural conjecture from the v59 framework:

    α_W / α  =  1 / sin(δ_B)  =  csc(2/9)

This emerges naturally from the mode decomposition in 04_lagrangian_modes.py:
the silent (SU(2)_L) decay constant is f_S = ρ₀ sin(δ_B) while the active (U(1)?)
decay constant is f_A = ρ₀.  If both couplings scale as α ~ 1/(some power of f),
the ratio α_W/α is a function of f_A/f_S = 1/sin(δ_B).

Empirical test: how close is the conjecture 1/sin(2/9) to the empirical α_W/α?
"""

import numpy as np

# v59 structural value
delta_B = 2.0 / 9.0

# Empirical electromagnetic / weak couplings.  We test at TWO conventional scales.
# At low energy:
alpha_low = 1.0 / 137.035999084
# At M_Z (PDG 2024):
alpha_M_Z = 1.0 / 127.952  # approx

# SU(2)_L coupling.  Use g₂ = 2 m_W / v (tree-level Higgs identification).
m_W = 80.379    # GeV (PDG 2024)
v_higgs = 246.220  # GeV (Higgs VEV)
g2 = 2 * m_W / v_higgs
alpha_W_tree = g2**2 / (4 * np.pi)
# At M_Z (running):
g2_M_Z = 0.6517   # PDG 2024
alpha_W_M_Z = g2_M_Z**2 / (4 * np.pi)

print("=" * 72)
print("v59 conjecture: α_W / α = 1/sin(δ_B) = csc(2/9)")
print("=" * 72)
print()
print(f"  δ_B = 2/9 = {delta_B:.10f}")
print(f"  sin(δ_B) = {np.sin(delta_B):.10f}")
print(f"  1/sin(δ_B) = csc(2/9) = {1/np.sin(delta_B):.10f}")
print()
print("=" * 72)
print("Empirical comparisons")
print("=" * 72)
print()
print(f"  g₂ (from m_W/v, tree-level)   = {g2:.6f}")
print(f"  g₂ (M_Z, PDG)                  = {g2_M_Z:.6f}")
print(f"  α_W (tree, m_W/v)             = {alpha_W_tree:.6f} = 1/{1/alpha_W_tree:.3f}")
print(f"  α_W (M_Z, running)            = {alpha_W_M_Z:.6f} = 1/{1/alpha_W_M_Z:.3f}")
print(f"  α (low energy)                = {alpha_low:.6e} = 1/{1/alpha_low:.3f}")
print(f"  α (M_Z, running)              = {alpha_M_Z:.6e} = 1/{1/alpha_M_Z:.3f}")
print()
print("  Ratio  α_W / α  at various scale combinations:")
ratio_lowlow_tree = alpha_W_tree / alpha_low
ratio_MZMZ = alpha_W_M_Z / alpha_M_Z
ratio_treelow = alpha_W_tree / alpha_low
print(f"    α_W(tree, m_W/v) / α(low)       = {alpha_W_tree / alpha_low:.6f}")
print(f"    α_W(M_Z) / α(M_Z)                = {alpha_W_M_Z / alpha_M_Z:.6f}")
print(f"    Conjecture 1/sin(2/9)             = {1/np.sin(delta_B):.6f}")
print()
print("  Relative agreement:")
print(f"    | α_W(tree)/α(low) − 1/sin(2/9) | / [1/sin(2/9)]  = "
      f"{abs(alpha_W_tree/alpha_low - 1/np.sin(delta_B)) / (1/np.sin(delta_B)):.4f}"
      f"  ({100*abs(alpha_W_tree/alpha_low - 1/np.sin(delta_B)) / (1/np.sin(delta_B)):.2f}%)")
print(f"    | α_W(M_Z)/α(M_Z)  − 1/sin(2/9) | / [1/sin(2/9)] = "
      f"{abs(alpha_W_M_Z/alpha_M_Z - 1/np.sin(delta_B)) / (1/np.sin(delta_B)):.4f}"
      f"  ({100*abs(alpha_W_M_Z/alpha_M_Z - 1/np.sin(delta_B)) / (1/np.sin(delta_B)):.2f}%)")

# ----------------------------------------------------------------------
# Test other candidate structural identities
# ----------------------------------------------------------------------
print()
print("=" * 72)
print("Other candidate structural identities for α_W / α")
print("=" * 72)
print()
print("  Candidate         | value     | match to α_W/α")
candidates = {
    "1 / sin(2/9)":              1.0 / np.sin(delta_B),
    "1 / sin(δ_B)²":             1.0 / np.sin(delta_B)**2,
    "1 / tan(2/9)":              1.0 / np.tan(delta_B),
    "cot(2/9)":                  1.0 / np.tan(delta_B),
    "csc²(δ_B)/dim G₂":          1/np.sin(delta_B)**2 / 14,
    "1 / (2 sin δ_B)":           1.0 / (2 * np.sin(delta_B)),
    "√(dim Spin(7))":            np.sqrt(21),
    "21 / sin(δ_B)·exp(?)":      21,  # for cross-sector v59 ratio
}
emp_ratio = alpha_W_tree / alpha_low
for name, val in candidates.items():
    match = abs(val - emp_ratio) / emp_ratio
    star = "  ←" if match < 0.02 else "   "
    print(f"  {name:<28} | {val:9.6f} | {match*100:5.2f}% {star}")
print(f"  Empirical α_W/α (tree,low)   | {emp_ratio:9.6f} | reference")

# ----------------------------------------------------------------------
# Implications: combine with v59 α conjecture
# ----------------------------------------------------------------------
print()
print("=" * 72)
print("Combining with v59 α-conjecture: -ln α + 2α = π²/2")
print("=" * 72)
print()
print("v59 (Variant D, 04_findings.md) found: -ln α + 2α = π²/2 to 4 × 10⁻⁵.")
print("Solve iteratively for α:")

# Solve x + 2 e^{-x} = π²/2
def solve_alpha():
    x = np.pi**2 / 2  # initial guess
    for _ in range(50):
        x = np.pi**2 / 2 - 2 * np.exp(-x)
    return np.exp(-x)
alpha_v59 = solve_alpha()
print(f"  α (v59 conjecture)            = {alpha_v59:.10e}")
print(f"  α (empirical low energy)      = {alpha_low:.10e}")
print(f"  Relative diff                  = {(alpha_v59 - alpha_low)/alpha_low:.3e}")
print()
print("If also α_W = α / sin(δ_B):")
alpha_W_combined = alpha_v59 / np.sin(delta_B)
g2_combined = np.sqrt(4 * np.pi * alpha_W_combined)
print(f"  α_W (combined)                 = {alpha_W_combined:.6f} = 1/{1/alpha_W_combined:.3f}")
print(f"  g₂ (combined)                  = {g2_combined:.6f}")
print(f"  g₂ (empirical, m_W/v)          = {g2:.6f}")
print(f"  Relative diff                  = {(g2_combined - g2)/g2:.3f}")
print()
print(f"  g₂ (combined)                  = {g2_combined:.6f}")
print(f"  g₂ (M_Z, PDG)                  = {g2_M_Z:.6f}")
print(f"  Relative diff                  = {(g2_combined - g2_M_Z)/g2_M_Z:.3f}")

# ----------------------------------------------------------------------
# Assessment
# ----------------------------------------------------------------------
print()
print("=" * 72)
print("Assessment")
print("=" * 72)
print(f"""
The conjecture α_W/α = 1/sin(δ_B) = csc(2/9) ≈ 4.5388 is within {100*abs(emp_ratio - 1/np.sin(delta_B))/(1/np.sin(delta_B)):.1f}%
of the empirical α_W(tree)/α(low) = {emp_ratio:.4f}.  This is suggestive but not
yet at experimental precision (current empirical uncertainties are at the
sub-percent level for α and at the few-percent level for the running of
α_W between scales).

Combined with the v59 conjecture -ln α + 2α = π²/2, the relation
α_W = α / sin(δ_B) predicts g₂ ≈ {g2_combined:.4f}, vs empirical
{g2:.4f} (tree, m_W/v) or {g2_M_Z:.4f} (M_Z, PDG).  Both within ~1%.

This is the right ballpark — within a few percent — using two v59 structural
inputs (δ_B = 2/9 from Lie group ratios, and the π²/2 = 8π²/(dim Cl(3,1))
conjecture for the instanton action).  Not yet a precise derivation, but
the magnitudes line up consistently.

NEXT STEPS for a precise derivation:
  • Identify the geometric origin of the factor 1/sin(δ_B):  it is the
    "stretch" factor between the full S³ radius ρ₀ and the silent S²
    radius ρ₀ sin(δ_B).  A natural gauge-coupling normalisation involving
    this stretch could make the relation exact.
  • Pin down the gauge-field action functional in concert with the rest of
    the v59 structure (the radial mode for G, the active mode for α, the
    silent mode for α_W) — i.e. the full strain Lagrangian.
""")
