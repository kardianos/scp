#!/usr/bin/env python3
"""
13_single_source.py

Test the single-source hypothesis: a SINGLE algebra projects to give all
three sector dimensions D_N = 28, 35, 63.

CANDIDATE: Cl(7)_even ≅ Cl(6) ≅ ℂ⊗𝕆 = Furey color algebra (dim 64).

Decomposition of Cl(7)_even by Clifford grade:
    Λ⁰ = 1   (identity)
    Λ² = 21  = dim Spin(7)        ← already v59-structural
    Λ⁴ = 35  = dim Λ⁴ℝ⁷ = (7 choose 4)
    Λ⁶ = 7   = dim S⁷
    Total: 64 = dim Cl(6) = dim Furey color algebra

The three sectors' D_N values arise as DIFFERENT subspaces:

    Lepton (N=0):  Λ² ⊕ Λ⁶            = 21 + 7  = 28  ← D_0
    d-quark (N=1):  Λ⁴                 =      35    = 35  ← D_1
    u-quark (N=2):  Λ² ⊕ Λ⁴ ⊕ Λ⁶      = 21+35+7 = 63  ← D_2

ALL THREE D_N values are now subspaces of the SAME parent: Cl(7)_even ≅ Cl(6).

This is the "single projected link" — one algebra, three projections (one
per Furey N-sector).

The G₂ orbit (14 dims) lives inside Cl(7)_even as a sub-piece (G₂ ⊂ Spin(7)
acts on Cl(7) preserving the associative 3-form).  Its 14-dim sub-orbit
restricts to each sub-grade structure of Cl(7)_even, giving the universal
14 in every numerator.

This script verifies the decomposition and proposes a rule for which grades
each sector picks.
"""

import math
from fractions import Fraction

# ---------------------------------------------------------------------
# Cl(7) grade decomposition
# ---------------------------------------------------------------------
Cl7_grades = {k: math.comb(7, k) for k in range(8)}
Cl7_even_grades = {k: v for k, v in Cl7_grades.items() if k % 2 == 0}

print("=" * 72)
print("Cl(7) grade structure")
print("=" * 72)
print(f"  Total dim Cl(7): {sum(Cl7_grades.values())} = 2^7 = 128\n")
print(f"  Grade-by-grade:")
for k, v in Cl7_grades.items():
    parity = "even" if k % 2 == 0 else "odd "
    star = " ← IN OUR ANALYSIS" if k in {2, 4, 6} else ""
    print(f"    Λ^{k}ℝ⁷ = {v:3} ({parity}){star}")
print()
print(f"  Cl(7)_even (sum of even grades): {sum(Cl7_even_grades.values())} = 64 = dim Cl(6)")
print(f"  Cl(7)_even ≅ Cl(6) ≅ ℂ⊗𝕆 (Furey color algebra)")

# ---------------------------------------------------------------------
# The single-source decomposition
# ---------------------------------------------------------------------
print()
print("=" * 72)
print("Sector-by-sector decomposition INSIDE Cl(7)_even (single source)")
print("=" * 72)

dim_G2 = 14
Lambda_2 = Cl7_grades[2]  # 21
Lambda_4 = Cl7_grades[4]  # 35
Lambda_6 = Cl7_grades[6]  # 7
Lambda_0 = Cl7_grades[0]  # 1

# The three sector decompositions
lepton_grades = [2, 6]
d_quark_grades = [4]
u_quark_grades = [2, 4, 6]

def sum_grades(grade_list):
    return sum(Cl7_grades[k] for k in grade_list)

print(f"\n  Lepton (N=0):   uses grades {lepton_grades} → dim {sum_grades(lepton_grades)}")
print(f"                  Λ² ⊕ Λ⁶ = {Lambda_2} + {Lambda_6} = {Lambda_2 + Lambda_6}  (= D_0 = 28 ✓)")
print()
print(f"  d-quark (N=1):  uses grades {d_quark_grades} → dim {sum_grades(d_quark_grades)}")
print(f"                  Λ⁴ = {Lambda_4}  (= D_1 = 35 ✓)")
print()
print(f"  u-quark (N=2):  uses grades {u_quark_grades} → dim {sum_grades(u_quark_grades)}")
print(f"                  Λ² ⊕ Λ⁴ ⊕ Λ⁶ = {Lambda_2} + {Lambda_4} + {Lambda_6} = {Lambda_2 + Lambda_4 + Lambda_6}  (= D_2 = 63 ✓)")

# ---------------------------------------------------------------------
# Verify the t² and Q predictions
# ---------------------------------------------------------------------
print()
print("=" * 72)
print("Brannen t² and Q from the single-source decomposition")
print("=" * 72)

for name, grades in [("lepton (N=0)", lepton_grades),
                     ("d-quark (N=1)", d_quark_grades),
                     ("u-quark (N=2)", u_quark_grades)]:
    D = sum_grades(grades)
    one_minus_t = Fraction(dim_G2, D)
    t_sq = 1 - one_minus_t
    Q = (1 + 2 * t_sq) / 3
    print(f"  {name}:")
    print(f"    D = {D}, 14/D = {one_minus_t}, t² = {t_sq}, Q = {Q}")

# ---------------------------------------------------------------------
# Selection-rule analysis
# ---------------------------------------------------------------------
print()
print("=" * 72)
print("Grade-selection rule by sector")
print("=" * 72)
print("""
Which Cl(7)_even grades does each sector select?

    Sector   N   Selects grades   Excludes from Cl(7)_even
    ──────  ──  ──────────────   ──────────────────────────
    e_R      0   {2, 6}           {0, 4}
    d_R      1   {4}              {0, 2, 6}
    u_R      2   {2, 4, 6}        {0}      (only excludes identity)

Pattern observations:
  • All sectors exclude Λ⁰ (identity, dim 1) — the "trivial" sub-grade.
  • d-quark selects ONLY Λ⁴ — the central even grade of Cl(7).
  • Lepton excludes Λ⁴ → selects {Λ², Λ⁶}.
  • u-quark selects everything non-identity: {Λ², Λ⁴, Λ⁶}.

Possible structural rule:
  • Λ⁴ is the "form grade" where the G₂ structure naturally lives
    (Λ⁴ is dual to Λ³ in Cl(7); Λ³ holds the associative 3-form).
  • d-quark (color triplet, N=1) couples to Λ⁴ specifically because
    its quantum number selects the trivector/4-vector content.
  • Lepton (color singlet) couples to Λ²+Λ⁶ — the "structural" grades
    (Λ² = so(7) = dim Spin(7), Λ⁶ = dim S⁷ via Hodge dual).
  • u-quark couples to everything non-identity (the full color sector).

The Furey N quantum number plausibly DRIVES this selection via the
Witt decomposition of Cl(6) — operators with different Fock-number act on
different Cl-grades.
""")

# ---------------------------------------------------------------------
# The deep correspondence: each sector vs Cl(7) grade
# ---------------------------------------------------------------------
print()
print("=" * 72)
print("Deeper interpretation: each grade has a v59-structural identity")
print("=" * 72)
print(f"""
The three non-identity even grades of Cl(7) each have a v59 identity:

  • Λ²ℝ⁷ = 21  =  dim so(7)  =  dim Spin(7)
       = v59 cross-sector ratio (S_grav / S_em = 21)
       = THE Spin(7) Lie algebra

  • Λ⁴ℝ⁷ = 35  =  (7 choose 4) = (7 choose 3) by Hodge duality
       = the 4-form space, dual to Λ³ℝ⁷ where the associative 3-form lives
       = appeared in dim G₂ = 49 − 35

  • Λ⁶ℝ⁷ = 7   =  dim S⁷ = Spin(7)/G₂ = Spin(8)/Spin(7)
       = the imaginary octonion dimension
       = appeared in dim G₂ = dim Spin(7) − dim S⁷

So ALL THREE non-identity even grades of Cl(7) are individually
v59-structural numbers.

EACH SECTOR DRAWS A DIFFERENT COMBINATION:
  • Lepton: the "structure-bearing" grades — Spin(7) and S⁷.
  • d-quark: ONLY the central form grade Λ⁴ ≅ Λ³.
  • u-quark: ALL THREE — the full structural content of Cl(7)_even minus
    identity.

THE UNIFIED SOURCE: Cl(7)_even (= Cl(6) = ℂ⊗𝕆 = Furey color algebra).
THE UNIVERSAL ORBIT: G₂ (dim 14) acts on each subspace with orbit dim 14.
THE EMERGENT DIMENSIONS: 28, 35, 63 from the sector-specific projections.
""")

# Save
import json
results = {
    "source_algebra": "Cl(7)_even = Cl(6) = ℂ⊗𝕆 = Furey color algebra",
    "source_dim": 64,
    "Cl7_grades": Cl7_grades,
    "sector_decomposition": {
        "lepton_N0": {"grades": lepton_grades, "D": 28, "decomp": "Λ² ⊕ Λ⁶ = 21 + 7"},
        "d_quark_N1": {"grades": d_quark_grades, "D": 35, "decomp": "Λ⁴ alone"},
        "u_quark_N2": {"grades": u_quark_grades, "D": 63, "decomp": "Λ² ⊕ Λ⁴ ⊕ Λ⁶ = 21 + 35 + 7"},
    },
    "common_orbit": "G₂ (dim 14)",
    "universal_invariant": "14 = dim G₂ adjoint orbit",
}
with open('/home/d/code/scp/v59/cosserat_experiment/13_single_source.json', 'w') as f:
    json.dump(results, f, indent=2)
print("Saved data to 13_single_source.json")
