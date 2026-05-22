#!/usr/bin/env python3
"""
09_unified_prefactors.py

Investigate the structural origin of the prefactors 5 (for g_W² = 5√α) and
21/16 (for G_e = (21/16)·α²¹).  Critical observation:

    5 = dim Spin(7) − dim Cl(3,1) = 21 − 16
    21/16 = dim Spin(7) / dim Cl(3,1)

Both prefactors use the SAME two structural numbers, but via different
operations (difference vs ratio).  This suggests a unified Lagrangian
mechanism in which these two ratios appear as different traces or
characters of the same underlying gauge structure.

We also look at the α prediction
    −ln α + 2α = π²/2 = 8π² / 16 = 8π² / dim Cl(3,1)
to see whether the "16 = dim Cl(3,1)" denominator there is the same as the
16 in the G prefactor.

The script also tests Lean's potential utility for each piece, since the
user asked whether Lean could help.
"""

import numpy as np
from fractions import Fraction

print("=" * 72)
print("Observation: both v59-tier prefactors use {dim Spin(7), dim Cl(3,1)}")
print("=" * 72)

dim_Spin7 = 21
dim_Cl31 = 16
dim_G2 = 14
dim_Spin8 = 28
dim_S7 = 7

print(f"\n  dim Spin(7)  = {dim_Spin7}    (= (7 choose 2), the rotation planes of ℝ⁷)")
print(f"  dim Cl(3,1)  = {dim_Cl31}    (= 2⁴, the spacetime Clifford algebra)")
print(f"  dim G_2       = {dim_G2}")
print(f"  dim Spin(8)   = {dim_Spin8}")
print(f"  dim S^7       = {dim_S7}    (= Spin(7)/G_2)")

print()
print(f"  difference:   dim Spin(7) − dim Cl(3,1) = {dim_Spin7 - dim_Cl31}   ← prefactor in g_W² = 5·√α")
print(f"  ratio:        dim Spin(7) / dim Cl(3,1) = {dim_Spin7/dim_Cl31:.6f} = {Fraction(dim_Spin7, dim_Cl31)}   ← prefactor in G_e = (21/16)·α²¹")
print()
print(f"  Note: in v59's α conjecture, π²/2 = 8π²/16 = 8π² / dim Cl(3,1).")
print(f"        So 16 is already structurally relevant in the EM-sector formula.")

# ----------------------------------------------------------------------
# Test 1: pattern of "natural" prefactor forms
# ----------------------------------------------------------------------
print()
print("=" * 72)
print("Pattern: a unified Lagrangian prefactor formula?")
print("=" * 72)

print("""
α prediction:   exp(−S_em)  =  α
                S_em        =  π²/2  =  8π² / dim Cl(3,1)        (no prefactor)

g_W prediction: g_W²        =  5 · √α  =  (dim Spin(7) − dim Cl(3,1)) · α^(1/2)
                α_W         =  (5/(4π)) · √α
                           =  ((dim Spin(7) − dim Cl(3,1)) / Vol(S²)) · α^(1/2)

G prediction:   G_e         =  (21/16) · α²¹
                            =  (dim Spin(7) / dim Cl(3,1)) · α^(dim Spin(7))

Trial unified form:
    (a) α       :  prefactor = ?  exponent = 1.
    (b) α_W     :  prefactor =  (dim Spin(7) − dim Cl(3,1)) / Vol(S²) =  5/(4π)
                  exponent = 1/2.
    (c) G       :  prefactor =  dim Spin(7) / dim Cl(3,1) =  21/16
                  exponent = dim Spin(7) = 21.

Conjectured uniform pattern: exponent equals the dim of the "natural
gauge object" for the sector.  For EM (a single U(1) generator): exponent 1.
For SU(2)_L on the silent S² with curvature: exponent 1/2 (from Pauli rep?).
For Spin(7) on the radial sector: exponent 21.

The prefactor structure is more elusive — it's a *difference* for g_W (with
S² denominator) and a *ratio* for G (with Cl(3,1) denominator).  We try
several functional forms to see if a single pattern fits.
""")

# ----------------------------------------------------------------------
# Test 2: Try unified algebraic forms
# ----------------------------------------------------------------------
print("=" * 72)
print("Test: can both prefactors arise from a single bilinear in (dim Spin(7), dim Cl(3,1))?")
print("=" * 72)

# Empirical values
alpha = 1 / 137.035999084
g2_emp = 0.6517
alpha_W_emp = g2_emp**2 / (4 * np.pi)
G_e_emp = 6.67430e-11 * (9.1093837015e-31)**2 / (1.054571817e-34 * 2.99792458e8)

# Trial: prefactor for X = A·dim_Spin(7) + B·dim_Cl(3,1) for some (A, B)
print("\n  For g_W: empirical α_W·4π/√α = {:.6f}".format(alpha_W_emp * 4 * np.pi / np.sqrt(alpha)))
print(f"  Target prefactor (g_W): 5 = (1)·21 + (-1)·16 ✓ ")
print(f"  Target prefactor (G):  21/16 = (1)·21 / (1)·16 ✓ ")

# A common form: numerator and denominator are linear combinations of (21, 16)?
print("\n  Try: G_e prefactor = num/denom with num, denom integer combinations of (21, 16)")
print()
print(f"  21/16 (Spin(7)/Cl(3,1))   = {21/16:.8f}  ← G match")
print(f"  (21-16)/16 = 5/16          = {5/16:.8f}")
print(f"  (21+16)/16 = 37/16         = {37/16:.8f}")
print(f"  21/(21-16) = 21/5           = {21/5:.6f}")
print(f"  (21-16)·4π / (21+16) ?     = {5 * 4 * np.pi / 37:.6f}")
print(f"  21·4π / 16² = 84π/256       = {84*np.pi/256:.6f}")

# Check the inverse: 1 + (Spin(7) - Cl(3,1))/Cl(3,1) = 1 + 5/16 = 21/16 ✓
print("\n  IDENTITY:  21/16 = 1 + (21−16)/16 = 1 + (dim Spin(7) − dim Cl(3,1))/dim Cl(3,1)")
print(f"             1 + 5/16 = {1 + 5/16}  (matches 21/16)")

# ----------------------------------------------------------------------
# Test 3: Maybe both prefactors come from a "stretch factor" / "(1 + δ)"
# ----------------------------------------------------------------------
print()
print("=" * 72)
print("Unified form candidate: (1 + δ_X) where δ_X depends on sector X")
print("=" * 72)
print(f"""
For G_e: prefactor = 1 + 5/16 = 1 + (dim Spin(7) − dim Cl(3,1))/dim Cl(3,1)

For g_W: prefactor = 5/(4π).  Is this also (1 + δ_W)?
  1 + δ_W = 5/(4π) = {5/(4*np.pi):.6f}  →  δ_W = {5/(4*np.pi) - 1:.6f}  (negative — wrong shape)

  Alternative: α_W (4π/5) = 1/√α
                4π/5 = 1/√α at v59-natural α implies √α = 5/(4π)?
                √α empirical (low-E) = {np.sqrt(alpha):.6f}, 5/(4π) = {5/(4*np.pi):.6f}.
                Off by factor {5/(4*np.pi) / np.sqrt(alpha):.4f}.

  Hmm — the "5/(4π)" prefactor in g_W is NOT of the form (1 + δ).  It's a
  rational number divided by 4π (the area of S²).

Possible interpretation:
  • G_e prefactor = "Clifford-algebra correction factor" 1 + (5/dim Cl(3,1))
  • g_W prefactor = "Killing-form / Vol(S²) factor" 5/(4π)
  • Both contain the SAME 5 = dim Spin(7) − dim Cl(3,1), but in different
    geometric normalisations (one wrt the Clifford algebra, the other wrt the
    silent target sphere).

This is the cleanest unified shape we can identify without a full
Lagrangian derivation.
""")

# ----------------------------------------------------------------------
# Test 4: α formula — what's the origin of π²/2 = 8π²/16?
# ----------------------------------------------------------------------
print("=" * 72)
print("α prediction: −ln α + 2α = π²/2 = 8π²/16")
print("=" * 72)

print(f"""
The v59 conjecture S_em = π²/2 is structurally interpreted as

    S_em = 8π² / dim Cl(3,1) = 8π² / 16

The 8π² is the Yang-Mills instanton coefficient — for an SU(N) BPST instanton
on ℝ⁴, the action is 8π² k / g², where k is the topological charge.

If we identify g² = 1 (some natural normalisation) and k = 1 (single instanton),
S_inst = 8π².  Then dividing by dim Cl(3,1) gives π²/2.

The 2α correction: empirically the formula −ln α + 2α = π²/2 matches at
4 × 10⁻⁵.  The 2α has the form of a one-loop running correction, but the
coefficient 2 doesn't obviously come from a standard β-function.

Candidate origins of the "2":
  • The number of polarisations of a massless photon (2).
  • Two copies of something in the algebra (e.g., even and odd Cl grades).
  • Or just an empirical correction matching to the precision we have.

Lean's role for this prediction:
  • Cannot derive π²/2 from a Lagrangian — that requires the full
    quantum-field-theory machinery for instantons, which is not in Mathlib.
  • CAN state the conjecture as a transcendental equation and verify
    numerical solutions to machine precision.
  • CAN encode the relation 8π² / dim Cl(3,1) = π²/2 once dim Cl(3,1) = 16
    is given as a structural input.
""")

# ----------------------------------------------------------------------
# Test 5: Where Lean would actually help
# ----------------------------------------------------------------------
print("=" * 72)
print("What Lean CAN do for these conjectures")
print("=" * 72)
print("""
For the EXPLORATORY phase (finding the right prefactors):
  • Lean is NOT the right tool.  We're searching for structural relations
    by symbolic and numerical exploration — Python is faster.
  • Once a candidate is found (like 21/16), Lean can VERIFY it as an
    arithmetic identity (e.g., `dimSO 7 / dimCl 31 = 21/16` is one line).

For the DERIVATION phase (deriving prefactors from Lagrangians):
  • Lean could in principle formalise the action functional and prove that
    its variational equations give the predicted couplings.
  • This requires substantial Mathlib formalisation of:
      - Yang-Mills theory on principal bundles
      - Sigma models with non-trivial target spaces
      - Instanton solutions in finite-dimensional cases
    None of which is currently in Mathlib at sufficient depth.
  • A reasonable intermediate goal: formalise the *algebraic* identities
    (like Killing-form embedding index calculations) that the predictions
    rely on, even if the full QFT derivation is out of reach.

For the v59 PREDICTION TABLE (consolidating results):
  • Lean is GREAT for this: encode each conjecture (with its structural
    inputs and the empirical match) as a theorem statement, even if the
    proof requires an empirical input.  This makes the prediction tier
    machine-checkable for consistency.
  • E.g.: theorem v59_G_prediction : (21/16 : ℝ) * alpha^21 = G_e ± ε
    can be stated, with ε quantifying the gap to empirical G_e.

CONCRETE NEXT STEP for Lean: encode the v59 prediction table as a Lean file
with each conjecture as a (currently unproved) theorem.  This is a few
hundred lines of work that gives us a machine-readable consolidated
prediction record.

For DERIVING the prefactors: substantial more Lie algebra in Lean would
help — specifically, formalising the Killing form on so(7), the embedding
index calculation, and the trace identities.  This is a real Mathlib lift.
""")

# Save
np.savez('/home/d/code/scp/v59/cosserat_experiment/09_unified.npz',
         dim_Spin7=21, dim_Cl31=16, dim_G2=14,
         diff_21_16=5, ratio_21_16=21/16,
         alpha=alpha, alpha_W_emp=alpha_W_emp, G_e_emp=G_e_emp)
print("Saved data to 09_unified.npz")
