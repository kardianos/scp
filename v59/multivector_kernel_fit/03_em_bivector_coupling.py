#!/usr/bin/env python3
"""
v59/multivector_kernel_fit/03_em_bivector_coupling.py

Step 5: cross-sector test.  The same Cl(3,1) algebra that gave us the
charged lepton spectrum in steps 1-4 should, if v58's "unified
multivector force" claim is to be made quantitative, predict the
electromagnetic coupling alpha = 1/137.036 from its bivector-grade
content.

Honest framing: this is the white whale of GA-based physics.  Eddington,
Wyler, Atiyah and many others have attempted to derive alpha from
geometric volume ratios.  Wyler (1969) reached numerical agreement to
~5 digits but was rejected because the choice of homogeneous space was
arbitrary.  No clean derivation exists.

This script does NOT claim to solve the problem.  It does:
  1. Compute the natural bivector-vector coupling strength inside
     Cl(3,1) under standard normalization.
  2. Compare to the empirical alpha.
  3. Identify the gap and list mechanisms (Wyler-style volume ratios,
     larger algebras, topology) that could close it.

The goal is to know exactly what's missing, so that step 6 or a later
session can attack it with concrete machinery.
"""

import numpy as np

ALPHA_INV = 137.035999084
ALPHA = 1 / ALPHA_INV

print("="*72)
print("Step 5: Cross-sector test -- alpha from Cl(3,1) bivector grade")
print("="*72)

# =========================================================================
# Part A: Natural bivector-vector couplings in Cl(3,1)
# =========================================================================
print()
print("-"*72)
print("Part A: Natural bivector-vector coupling magnitudes in Cl(3,1)")
print("-"*72)

# In Cl(3,1) the bivector grade is 6-dimensional (e_12, e_13, e_14, e_23, e_24, e_34).
# A vector field psi (spatial part) is 3-dim.  The natural coupling is
# the geometric product, which for a bivector B and a vector v gives:
#     B v   =  (B . v)  +  (B ^ v)
# decomposed into vector and trivector parts.
#
# A Lagrangian of the form L_int = g * psi-bar B psi has coupling g.
# In natural geometric-algebra normalization (B normalized to unit norm,
# psi to unit norm), g is determined by the structure constants.
#
# For a unit bivector B with B^2 = -1, and unit vector psi,
#     |B psi|^2 = |B|^2 |psi|^2 - (B.psi)^2  -- but B.psi = 0 (bivector and
# vector are orthogonal grades), so |B psi|^2 = 1.
#
# The coupling magnitude in natural units is 1.

print("\nNatural normalization in Cl(3,1):")
print("  Unit bivector B, |B|^2 = 1 (B^2 = -1 for bivector)")
print("  Unit vector psi, |psi|^2 = 1")
print("  Geometric product B psi has |B psi| = 1 in this normalization.")
print()
print("So the 'natural' coupling magnitude is g_nat = 1.")
print()

# Compare to electromagnetic coupling
g_em = np.sqrt(4 * np.pi * ALPHA)
print(f"Empirical EM coupling: g_em = sqrt(4 pi alpha) = {g_em:.6f}")
print(f"Natural Cl(3,1) coupling: g_nat = 1.0")
print(f"Ratio g_nat / g_em = {1 / g_em:.6f}")
print()
print("So Cl(3,1) bivector coupling is naturally O(1), while EM is O(1/3).")
print("Order-of-magnitude is right, but the precise 1/137 is not predicted.")


# =========================================================================
# Part B: Wyler-style volume ratios
# =========================================================================
print()
print("-"*72)
print("Part B: Wyler-style volume ratios")
print("-"*72)
print()
print("Wyler (1969) conjectured:")
print("  alpha^-1 = (9 pi^4 / 16) * (2 pi^5 / 120)^(1/4)")
print("with the factors arising from volumes of specific homogeneous spaces.")
print()

# Wyler's formula
W = (9 * np.pi**4 / 16) * (2 * np.pi**5 / 120)**(1/4)
print(f"  Wyler's value:        {W:.6f}")
print(f"  alpha^-1 (empirical): {ALPHA_INV:.6f}")
print(f"  delta:                {W - ALPHA_INV:+.6f}")
print(f"  relative delta:       {(W - ALPHA_INV) / ALPHA_INV:+.3e}")
print()
print("Wyler's formula is famously close (~5 digits) but the choice of")
print("spaces is widely regarded as arbitrary.  No first-principles")
print("derivation is accepted.")
print()
print("Within Cl(3,1), the analogous volumes are:")
print(f"  Volume of S^3 (unit 3-sphere)    = 2 pi^2     = {2 * np.pi**2:.6f}")
print(f"  Volume of S^4 (unit 4-sphere)    = 8 pi^2 / 3 = {8 * np.pi**2 / 3:.6f}")
print(f"  Volume of S^5                    = pi^3       = {np.pi**3:.6f}")
print(f"  Volume of S^7                    = pi^4 / 3   = {np.pi**4 / 3:.6f}")
print(f"  Volume of SU(2) (= S^3)           = 2 pi^2     = {2 * np.pi**2:.6f}")
print(f"  Volume of SO(3,1) connected      = (infinite, non-compact)")
print(f"  Volume of Spin(4) / Spin(3,1)?    = need analytic continuation")


# =========================================================================
# Part C: Brannen phase as a possible source
# =========================================================================
print()
print("-"*72)
print("Part C: Could the Brannen phase phi = 2/9 rad be related to alpha?")
print("-"*72)
print()

phi_B = 2.0 / 9.0
print(f"Brannen phi = {phi_B:.10f} rad")
print(f"alpha       = {ALPHA:.10f}")
print()
print("Direct ratios:")
print(f"  phi / alpha           = {phi_B / ALPHA:.6f}")
print(f"  phi * alpha           = {phi_B * ALPHA:.6e}")
print(f"  phi - alpha           = {phi_B - ALPHA:.6f}")
print(f"  phi^2                 = {phi_B**2:.6f}")
print(f"  phi^2 / alpha         = {phi_B**2 / ALPHA:.6f}")
print(f"  alpha / phi^2         = {ALPHA / phi_B**2:.6f}")
print(f"  exp(-phi/alpha)       = {np.exp(-phi_B / ALPHA):.6e}")
print()
print(f"  1/alpha = {ALPHA_INV:.6f}")
print(f"  1/phi   = {1/phi_B:.6f}")
print(f"  1/(alpha * phi)       = {1/(ALPHA * phi_B):.6f}")
print()
print("No clean ratio jumps out.  alpha^-1 ~ 137 and phi ~ 0.22 rad have")
print("very different scales; there is no obvious dimensionless combination")
print("that lands on either from the other plus simple math.")


# =========================================================================
# Part D: Where would alpha come from in this framework?
# =========================================================================
print()
print("-"*72)
print("Part D: Mechanisms that could give alpha in this framework")
print("-"*72)
print()
print("Cl(3,1) by itself gives O(1) couplings.  To get alpha = 1/137 from")
print("the same algebra requires additional structure.  Candidate sources:")
print()
print("1. VACUUM MANIFOLD VOLUME RATIO.  If the gauge group orbit through")
print("   the vacuum has a specific volume V, the effective coupling is")
print("   suppressed by 1/V.  For SU(2)/U(1) = S^2 the volume is 4 pi;")
print("   for larger groups much bigger.  Wyler's value 137 comes from")
print("   a homogeneous space volume.")
print()
print("2. RENORMALIZATION FLOW.  alpha runs from a UV value alpha_UV to")
print("   the IR value 1/137.  The UV value could be O(1) and natural,")
print("   while 1/137 emerges from RG flow.  This is the standard story.")
print()
print("3. INSTANTON / TOPOLOGICAL FACTOR.  Exponentially small factors")
print("   exp(-S_inst) can give 1/137 if the instanton action S_inst ~ 5.")
print("   This is unrelated to the lepton sector at first order.")
print()
print("4. EMBEDDING IN A LARGER ALGEBRA.  Cl(0,7) (octonions), Cl(3,1)*Cl(0,3),")
print("   or exceptional algebras may have natural volume / structure-constant")
print("   ratios that give 137.  This is the Furey/Dixon line of inquiry.")
print()
print("5. WYLER-STYLE GEOMETRIC RATIO.  Volumes of homogeneous spaces of")
print("   the Lorentz / Poincare group, when combined with normalizations,")
print("   have given 1/137 numerically.  No first-principles derivation.")


# =========================================================================
# Part E: What this scan accomplishes
# =========================================================================
print()
print("="*72)
print("Part E: Honest assessment")
print("="*72)
print("""
What this step has established:

  - Cl(3,1) naturally produces O(1) bivector-vector couplings.  This is
    the correct order of magnitude for a UV gauge coupling, but
    qualitatively different from the IR alpha = 1/137 we observe.

  - No simple combination of the Brannen phase phi = 2/9 rad and the
    fine-structure constant alpha gives a clean dimensionless number.
    The lepton sector and the EM sector, in this framework, do not have
    an obvious cross-sector identity at the 'TYCHO ellipse' level.

  - Wyler-style volume ratios are the historical near-miss approach.
    Numerically close but mathematically arbitrary.

What this step did NOT do:

  - It did NOT predict alpha.  Doing so requires either:
      (a) writing the full Lagrangian including gauge fields, computing
          the RG flow from a natural UV value to the IR, OR
      (b) identifying a specific homogeneous-space volume ratio whose
          value is 137.036 by construction, OR
      (c) embedding Cl(3,1) in a larger algebra where 1/137 emerges
          as a natural structure constant.

  - None of these are one-session projects.

What this step DOES rule in/out:

  - Cl(3,1) Z_3 kernel + Z_3 invariant potential is structurally
    incomplete for predicting alpha.  Additional content (gauge group,
    larger algebra, or specific topology) is required.

  - The framework is NOT internally inconsistent at the level tested.
    The lepton-sector parameters do not constrain alpha within Cl(3,1)
    alone, so we can't yet falsify the unification claim.  We can also
    not confirm it.

CONSEQUENCE for the v58 'unified multivector force' claim:

  v58 claimed Newton + Maxwell emerge from one multivector equation.
  At the qualitative/formal level, yes -- different grade projections
  give different sector behaviors.  At the QUANTITATIVE level, the
  unification predicts no relation between the lepton spectrum and the
  EM coupling.  The two sectors are decoupled in Cl(3,1) without
  additional structure.

  This is the gap.  To close it requires either accepting that the
  unification is qualitative (in which case TYCHO-fitting in cross-
  sectors must be done sector-by-sector) or extending the framework
  beyond Cl(3,1).
""")

print("="*72)
print("Step 5 complete.")
print("="*72)
