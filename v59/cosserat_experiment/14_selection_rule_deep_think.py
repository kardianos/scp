#!/usr/bin/env python3
"""
14_selection_rule_deep_think.py

Deep think on the Furey-N → Cl(7)_even grade selection rule.

The user's skepticism is reasonable: there's no obvious mechanism that forces
sector N=0 to couple to {Λ², Λ⁶}, N=1 to {Λ⁴}, and N=2 to {Λ², Λ⁴, Λ⁶}.

The selection rule, if real, must come from one of:
  (a) A specific Lagrangian whose Yukawa structure mandates the rule.
  (b) A symmetry breaking pattern in the Furey ℂ⊗ℍ⊗𝕆 algebra that picks
      specific grades for each N.
  (c) A topological invariant (e.g. winding number) that varies with N.
  (d) A genuine empirical coincidence that we're over-interpreting.

This script explores several hypotheses without claiming derivation.
"""

import numpy as np
import math
from itertools import combinations

print("=" * 72)
print("Selection rule deep think")
print("=" * 72)

# Structural data
sectors = {
    "lepton_N0":   {"N": 0, "grades": {2, 6}, "D": 28, "physical": "color singlet"},
    "d_quark_N1":  {"N": 1, "grades": {4},    "D": 35, "physical": "color triplet, charge ±1/3"},
    "u_quark_N2":  {"N": 2, "grades": {2, 4, 6}, "D": 63, "physical": "color triplet, charge ∓1/3"},
}

print()
print("Observed selection patterns:")
for name, info in sectors.items():
    print(f"  {name}:  N = {info['N']}, grades = {info['grades']}, D = {info['D']}")
    print(f"          ({info['physical']})")

# ---------------------------------------------------------------------
# Hypothesis 1: Witt bidegree (number-preserving operators)
# ---------------------------------------------------------------------
print()
print("=" * 72)
print("Hypothesis 1: Witt bidegree counts")
print("=" * 72)
print("""
In Cl(6) ≅ ℂ⊗𝕆, the Witt basis (α_i, α̅_i) gives a bidegree (p, q) on each
operator: p = creation count, q = annihilation count.  Number-preserving
operators have p = q.

Number-preserving operators by Cl(6) grade:
  Grade 0:  (0,0) bidegree only — 1 element (identity)
  Grade 2:  (1,1) — 3 × 3 = 9 operators (α_i α̅_j)
  Grade 4:  (2,2) — (3 choose 2)² = 9 operators
  Grade 6:  (3,3) — 1 operator (full top form)

Total number-preserving: 1 + 9 + 9 + 1 = 20.

This 20 is NOT 28, 35, or 63.  And it's the SAME for all Furey N (since
it's just the subalgebra commuting with the number operator).

CONCLUSION: Hypothesis 1 doesn't give the D_N values.  The "ambient" D_N
is NOT the number-preserving subalgebra.
""")

# Verify counting
def count_bidegree(p, q):
    """Count operators with p creation and q annihilation (using 3 mode types)."""
    if p > 3 or q > 3:
        return 0
    return math.comb(3, p) * math.comb(3, q)

print("  Bidegree counts (3 α's, 3 α̅'s available):")
for k in range(7):
    pairs = [(p, k-p) for p in range(k+1) if k-p >= 0 and k-p <= 3 and p <= 3]
    elements = sum(count_bidegree(p, q) for p, q in pairs)
    np_elements = sum(count_bidegree(p, p) for p in [k//2] if k % 2 == 0)
    print(f"    Cl(6) grade {k}: total {elements}, number-preserving {np_elements}")

# ---------------------------------------------------------------------
# Hypothesis 2: Hodge dual / Poincaré duality structure
# ---------------------------------------------------------------------
print()
print("=" * 72)
print("Hypothesis 2: Hodge duality structure of Cl(7)")
print("=" * 72)
print("""
Cl(7) has Hodge dual (∗): Λ^k ↔ Λ^(7-k).
Within Cl(7)_even {Λ⁰, Λ², Λ⁴, Λ⁶}, Hodge duals leave even grades but go to
odd grades of Cl(7):  Λ⁰ ↔ Λ⁷, Λ² ↔ Λ⁵, Λ⁴ ↔ Λ³, Λ⁶ ↔ Λ¹.

So Hodge dual is NOT a symmetry of Cl(7)_even.

But within Cl(7)_even there's a different involution: complement-to-7-grade.
For Λ²k, the "co-grade" is 6−2k (within even grades):
  Λ⁰ ↔ Λ⁶  (k=0 ↔ k=3)
  Λ² ↔ Λ⁴  (k=1 ↔ k=2)

Under this Z₂ involution:
  Λ⁰ pair with Λ⁶
  Λ² pair with Λ⁴

The lepton selection {Λ², Λ⁶}: ONE element from each pair.
The d-quark selection {Λ⁴}: ONE element from one pair.
The u-quark selection {Λ², Λ⁴, Λ⁶}: three elements.

Doesn't immediately give a clean rule via this pairing.
""")

# ---------------------------------------------------------------------
# Hypothesis 3: Z_3 triality interaction with Furey N
# ---------------------------------------------------------------------
print()
print("=" * 72)
print("Hypothesis 3: Z_3 triality × Furey N")
print("=" * 72)
print("""
Spin(8) has Z_3 triality (cyclic on V, S+, S-).  Three generations come
from this Z_3.

Spin(7) ⊂ Spin(8) is the SUBGROUP that breaks triality (preserves one of
the three 8-dim reps).

The Cl(7) decomposition we're using is associated with Spin(7), so
triality is *already* broken in our decomposition.  The Furey N quantum
number is a SECOND grading on top of this.

Possible interpretation: under the triality × Furey N combined structure,
each fermion sector ends up in a specific (triality direction, N-grade)
combination that selects specific Cl(7) grades.

This is suggestive but underdetermined without a Lagrangian.
""")

# ---------------------------------------------------------------------
# Hypothesis 4: Topological winding number
# ---------------------------------------------------------------------
print()
print("=" * 72)
print("Hypothesis 4: Topological winding number selecting grades")
print("=" * 72)
print("""
A topological invariant W ∈ ℤ might select Cl(7)_even grades by
"W mod something".

If W = 2N (twice the Furey N, since α/α̅ is grade 1 in Cl(6) but grade 2
under Cl(7)_even ≅ Cl(6) re-grading):
  N=0: W=0 → ?
  N=1: W=2 → ?
  N=2: W=4 → ?

If W ↔ Cl(7) grade directly:
  N=0 (W=0): Λ⁰?  But lepton skips Λ⁰.
  N=1 (W=2): Λ²?  But d-quark uses Λ⁴, not Λ².
  N=2 (W=4): Λ⁴?  But u-quark uses all three.

Doesn't fit directly.

If W is "Pontryagin number" of a G_2-bundle, it's a different integer
invariant.  Without more concrete identification, this is speculation.
""")

# ---------------------------------------------------------------------
# Hypothesis 5: NOVEL — Selection by "complement of N-graded operator algebra"
# ---------------------------------------------------------------------
print()
print("=" * 72)
print("Hypothesis 5 (NOVEL): selection by complement of N-graded subalgebra")
print("=" * 72)
print("""
The Cl(6) algebra acts on the Fock space ℂ⁸.  Each Fock state |Ω_N⟩ is
created by a specific number of α's, hence has Furey N as a quantum number.

NEW CANDIDATE RULE:
  For fermion type N, the AMBIENT for its Brannen ξ is the subspace of
  Cl(7)_even *complementary* to the operators that ANNIHILATE its
  Witt-state |Ω_N⟩.

Equivalently: the ambient is "operators that have non-trivial action on |f⟩".

For |Ω_0⟩ = |0⟩ (lepton):
  Annihilated by anything containing α̅_i on the right (since α̅_i |0⟩ = 0).
  Cl(6) elements that annihilate |0⟩: ⟨α̅_i, ᾱ_iᾱ_j, ᾱ_1ᾱ_2ᾱ_3⟩ and
  combinations with extra α's on the LEFT (which doesn't help).

Actually: an operator A annihilates |0⟩ iff A|0⟩ = 0.  For A in Cl(6),
this happens iff A contains an α̅ as its rightmost-acting operator (modulo
re-ordering).

The "non-annihilating" operators for |0⟩: those that DON'T have an α̅ at
the right.  Equivalently: products of α_i (creation only) or identity.

These have Cl(6) bidegrees (0,0), (1,0), (2,0), (3,0).
Counts: 1 + 3 + 3 + 1 = 8.

So the "non-annihilating" subalgebra for |0⟩ has dim 8.  Not 28.

CONCLUSION: This naive rule doesn't give D_0 = 28 either.

Modified rule: perhaps the ambient is the subspace of operators that
COMMUTE WITH the projector onto |Ω_N⟩ (= the *centralizer* of the
projector).

Cl(6) acting on ℂ⁸ ≅ End(ℂ⁸).  Projector P_0 onto |0⟩ is rank 1.
Its centralizer in End(ℂ⁸) has dim 1 + 7² = 50 (block-diagonal in (1, 7)
decomposition).  Not 28.

For 3-dim N=1 block: P_1 (projector onto N=1 subspace) is rank 3.
Centralizer: 3² + 5² = 9 + 25 = 34.  Close to 35 but not exact.

For 1-dim N=3 block: 1² + 7² = 50.

So centralizer dims of Furey N-projectors: 50, 34, 34, 50.  Not 28, 35, 63.

Closer to D_1 = 35 (34 vs 35) but not exact.
""")

# Verify centralizer counts
print("  Centralizer dim of rank-r projector in End(ℂ⁸):")
for r in [1, 3, 3, 1]:
    cent_dim = r**2 + (8 - r)**2
    print(f"    rank {r}: centralizer dim = {r}² + {8-r}² = {cent_dim}")

# ---------------------------------------------------------------------
# Hypothesis 6: NOVEL — selection by orbit structure under SU(3)_c
# ---------------------------------------------------------------------
print()
print("=" * 72)
print("Hypothesis 6 (NOVEL): SU(3)_c orbit structure")
print("=" * 72)
print("""
The Furey N-states transform under SU(3)_c as:
  N=0:  singlet  (1-dim)
  N=1:  triplet  (3-dim) = color triplet
  N=2:  anti-triplet  (3-dim)
  N=3:  singlet  (1-dim)

Maybe the ambient D_N is connected to the SU(3)_c-INVARIANT subspace of a
specific representation.

For SU(3) acting on:
  • Cl(7)_even (64-dim total): we can decompose into SU(3) reps.
  • Λ⁰ ⊕ Λ² ⊕ Λ⁴ ⊕ Λ⁶ as Cl-grades but SU(3) acts on them too.

For SU(3) ⊂ G_2 ⊂ Spin(7) ⊂ Spin(8):
  • Λ²ℝ⁷ = so(7) decomposes under SU(3) as 1 + 8 + 6 + 6̄ or similar
    (depending on which SU(3) we choose).
  • Λ⁴ℝ⁷ decomposes similarly.
  • Λ⁶ℝ⁷ = 7 decomposes as 1 + 3 + 3̄ (the fundamental of SU(3))... maybe.

If lepton sector restricts to SU(3)_c-SINGLET part of {Λ², Λ⁶}: the
singlet dimensions might add to give a different number than 28.

If d-quark sector restricts to TRIPLET part of {Λ⁴}: dim might match
something related to 35 or its color-triplet sub-component.

This is the most physically motivated hypothesis: different SU(3)_c
content selects different Cl-grade content.

Concrete test would require computing the SU(3)_c decomposition of each
Cl(7) grade explicitly, then checking whether the sector-specific selection
matches a "color content" filter.

For now, this is a CANDIDATE MECHANISM that would need a detailed
Spin(7) ⊃ SU(3) branching calculation to verify.
""")

# ---------------------------------------------------------------------
# Assessment
# ---------------------------------------------------------------------
print()
print("=" * 72)
print("Honest assessment")
print("=" * 72)
print("""
After exploring six hypotheses for the Furey N → Cl-grade selection rule:

  • Hypothesis 1 (Witt bidegree): DOES NOT FIT
  • Hypothesis 2 (Hodge duality): SUGGESTIVE BUT INCOMPLETE
  • Hypothesis 3 (Z_3 triality): UNDERDETERMINED
  • Hypothesis 4 (winding number): NOT DERIVED
  • Hypothesis 5 (centralizer): CLOSE BUT WRONG (34 vs 35)
  • Hypothesis 6 (SU(3)_c branching): MOST PHYSICALLY MOTIVATED, REQUIRES CALCULATION

The user's skepticism is justified — the selection rule is NOT
derivable from any simple structural observation.  It's likely encoded in
a specific Lagrangian or symmetry-breaking pattern that requires more
concrete v59 theory than we have.

WHAT WE HAVE CONFIRMED:
  • Single source Cl(7)_even ≅ ℂ⊗𝕆 (the Furey color algebra) — REAL.
  • Decomposition Λ⁰⊕Λ²⊕Λ⁴⊕Λ⁶ = 1+21+35+7 = 64 — REAL.
  • Three sector ambients arise as graded sums:
      28 = Λ²+Λ⁶
      35 = Λ⁴
      63 = Λ²+Λ⁴+Λ⁶
  • All three D_N values are inside the same parent algebra.
  • Match to empirical Q values: 0.26 – 0.34 %.

WHAT IS NOT YET DERIVED:
  • WHY each Furey N picks its specific subset of Cl(7)_even grades.

WHAT WOULD BE NEEDED FOR DERIVATION:
  1. An explicit Lagrangian for the Brannen Yukawa coupling that includes
     the Furey ℂ⊗ℍ⊗𝕆 structure with sector-specific projectors.
  2. A symmetry-breaking pattern (from G_2 → SU(3)_c × U(1) or similar)
     that naturally selects Cl-grades by SU(3)_c content.
  3. A topological invariant (e.g. Chern class of a v59-natural bundle)
     that differs across sectors and selects grades accordingly.

The most physical next step would be a Hypothesis 6 calculation:
explicit SU(3)_c decomposition of Λ²ℝ⁷, Λ³ℝ⁷, Λ⁴ℝ⁷, Λ⁶ℝ⁷ under the natural
SU(3) ⊂ G_2 ⊂ Spin(7), looking for sector-specific filtration.

This is beyond a Python script — it requires explicit representation
theory computations or a Lie-algebra computer algebra system (LiE,
Sage, or similar).

The user's skepticism that I'd find this in a session is well-founded.
""")

# Save
import json
results = {
    "hypotheses_tested": 6,
    "matches_found": 0,
    "closest_match": "Hypothesis 5 (centralizer)",
    "closest_match_gap": "34 vs 35 (3% off)",
    "what_we_confirmed": "Single source Cl(7)_even with 28+35+7=63 decomposition",
    "what_we_didnt_find": "The Furey-N → Cl-grade selection mechanism",
    "next_concrete_step": "Hypothesis 6: SU(3)_c decomposition of Cl(7) grades",
}
with open('/home/d/code/scp/v59/cosserat_experiment/14_selection.json', 'w') as f:
    json.dump(results, f, indent=2)
print()
print("Saved data to 14_selection.json")
