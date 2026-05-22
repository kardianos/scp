#!/usr/bin/env python3
"""
15_su3_branching.py

Attempt the SU(3)_c branching calculation for Cl(7)_even grades.

If the selection rule comes from SU(3)_c-content filtering, we need to know
how each Λ^k ℝ⁷ decomposes under SU(3) ⊂ G₂ ⊂ Spin(7).

Standard G₂ ⊃ SU(3) branching (long-root SU(3)):
  1   → 1
  7   → 1 + 3 + 3̄              (the imaginary octonion)
  14  → 8 + 3 + 3̄              (the G₂ adjoint = octonion derivations)
  27  → 1 + 8 + (3+3̄) + (6+6̄)?  (the 2-symmetric tensor; needs verification)

For Cl(7), we use the G₂-decomposition of each Λ^kℝ⁷:

  Λ⁰ℝ⁷ = 1                        (G₂ scalar)
  Λ¹ℝ⁷ = 7                        (G₂ fundamental)
  Λ²ℝ⁷ = 14 + 7                   (G₂ adjoint + fundamental, since so(7) = g₂ ⊕ 7)
  Λ³ℝ⁷ = 1 + 7 + 27               (G₂: trivector contains assoc 3-form + ... )
  Λ⁴ℝ⁷ = 1 + 7 + 27               (Hodge dual of Λ³)
  Λ⁵ℝ⁷ = 14 + 7                   (Hodge dual of Λ²)
  Λ⁶ℝ⁷ = 7                        (Hodge dual of Λ¹)
  Λ⁷ℝ⁷ = 1                        (Hodge dual of Λ⁰)

Verifying dims:
  Λ² = 14 + 7 = 21 ✓
  Λ³ = 1 + 7 + 27 = 35 ✓
  Λ⁴ = 1 + 7 + 27 = 35 ✓
  Λ⁶ = 7 ✓
"""

# Branching for each Cl(7) grade under G₂:
# Each Λ^k decomposes into G₂ irreps:
G2_branching = {
    0: [(1, "trivial")],
    1: [(7, "fundamental")],
    2: [(14, "adjoint"), (7, "fundamental")],
    3: [(1, "trivial"), (7, "fundamental"), (27, "27-irrep")],
    4: [(1, "trivial"), (7, "fundamental"), (27, "27-irrep")],
    5: [(14, "adjoint"), (7, "fundamental")],
    6: [(7, "fundamental")],
    7: [(1, "trivial")],
}

print("G₂-decomposition of Λ^k ℝ⁷:")
print()
print(f"  {'Grade':<6} {'G₂ irreps':<35} {'Dim'}")
print(f"  {'-'*6} {'-'*35} {'-'*5}")
for k, irreps in G2_branching.items():
    total = sum(dim for dim, _ in irreps)
    irrep_str = " + ".join(f"{dim}" for dim, _ in irreps)
    print(f"  Λ^{k}    {irrep_str:<35} {total}")

# Now SU(3) ⊂ G₂ branching of G₂ irreps:
# (Long-root SU(3) embedding)
SU3_branching = {
    1:  [(1, "1")],
    7:  [(1, "1"), (3, "3"), (3, "3̄")],
    14: [(8, "8"), (3, "3"), (3, "3̄")],
    27: [(1, "1"), (8, "8"), (3, "3"), (3, "3̄"), (6, "6"), (6, "6̄")],  # needs verification
}

print()
print("SU(3) ⊂ G₂ branching of G₂ irreps:")
print()
print(f"  {'G₂ irrep':<10} {'SU(3) content':<35} {'Dim'}")
print(f"  {'-'*10} {'-'*35} {'-'*5}")
for g2_dim, su3_content in SU3_branching.items():
    total = sum(d for d, _ in su3_content)
    content_str = " + ".join(f"{d}" for d, _ in su3_content)
    print(f"  {g2_dim:<10} {content_str:<35} {total}")

# Combine: full SU(3) decomposition of each Cl(7) grade
print()
print("=" * 72)
print("Combined: SU(3)_c decomposition of each Λ^k ℝ⁷")
print("=" * 72)

def su3_content_total(grade):
    """For Cl(7)-grade k, return aggregated SU(3) reps."""
    content = {}
    for g2_dim, _ in G2_branching[grade]:
        for su3_dim, su3_name in SU3_branching[g2_dim]:
            content[su3_name] = content.get(su3_name, 0) + su3_dim
    return content

print()
print(f"  {'Grade':<6} {'SU(3)_c content':<40} {'Dim check'}")
print(f"  {'-'*6} {'-'*40} {'-'*10}")
for k in range(8):
    content = su3_content_total(k)
    total = sum(content.values())
    content_str = " + ".join(f"{v}({n})" for n, v in content.items())
    print(f"  Λ^{k}    {content_str:<40} {total}")

# ---------------------------------------------------------------------
# Now check: do v59 sector D_N values match specific SU(3) content filters?
# ---------------------------------------------------------------------
print()
print("=" * 72)
print("Sector D_N values as SU(3) content filters")
print("=" * 72)

# For each Furey N sector, the ambient D_N is:
# N=0 (lepton): D = 28 = Λ²+Λ⁶ = 21 + 7
# N=1 (d-quark): D = 35 = Λ⁴
# N=2 (u-quark): D = 63 = Λ²+Λ⁴+Λ⁶

# Try filter: keep only SU(3)-singlet content
print("\n  Filter: SU(3)_c singlet content only")
for k in range(2, 7, 2):
    content = su3_content_total(k)
    singlet = content.get("1", 0)
    print(f"    Λ^{k} singlet content: {singlet}")
print()

# Try filter: SU(3)-singlet + octet (color-neutral content)
print("  Filter: singlet + octet (= color-neutral)")
for k in range(2, 7, 2):
    content = su3_content_total(k)
    neutral = content.get("1", 0) + content.get("8", 0)
    print(f"    Λ^{k} color-neutral content: {neutral}")
print()

# Try filter: only color-rich (triplet/anti-triplet/6/6̄)
print("  Filter: color-charged content only (3, 3̄, 6, 6̄)")
for k in range(2, 7, 2):
    content = su3_content_total(k)
    charged = sum(v for n, v in content.items() if n not in ("1", "8"))
    print(f"    Λ^{k} color-charged content: {charged}")

print("""
Observations:
  • Λ²ℝ⁷ singlet content = 1 (one singlet from "7" branching)
  • Λ²ℝ⁷ color-neutral (1+8) content = 1 + 8 = 9
  • Λ²ℝ⁷ color-charged content = 21 - 9 = 12 (= 2·3 + 2·3̄)

  • Λ⁴ℝ⁷ singlet content = 2 (one from "1", one from "7")  (?)
  • Λ⁴ℝ⁷ color-neutral content = (1+1) + (8) = 10  (?)
  • Λ⁴ℝ⁷ color-charged content = 35 - 10 = 25  (?)

These don't immediately give 28, 35, 63 by any clean filter.

Need to check: my SU(3) branching for 27 of G₂ might be wrong (it's
"6+6̄" + others — I should verify with a reference).
""")

# ---------------------------------------------------------------------
# Different approach: count "G₂-invariants" per grade
# ---------------------------------------------------------------------
print()
print("=" * 72)
print("G₂-invariant content per grade (G₂ singlets in each Λ^k)")
print("=" * 72)
print()

for k, irreps in G2_branching.items():
    g2_singlet_dim = sum(dim for dim, name in irreps if dim == 1)
    print(f"  Λ^{k}: contains {g2_singlet_dim} G₂-singlet(s) (each 1-dim)")

print("""
G₂-invariants on ℝ⁷:
  • Λ⁰: 1 (trivial)
  • Λ³: 1 (the associative 3-form φ — THE G₂-defining structure)
  • Λ⁴: 1 (the coassociative 4-form *φ — Hodge dual of φ)
  • Λ⁷: 1 (the volume form)

So there are 4 G₂-invariant directions in Λ*ℝ⁷: scalar, φ, *φ, volume.

These are NOT what the selection rule picks — sectors use larger subspaces.

So if the selection rule existed, it must be something more subtle than
simple "G₂-invariant content".
""")

# ---------------------------------------------------------------------
# Honest conclusion
# ---------------------------------------------------------------------
print()
print("=" * 72)
print("Honest conclusion")
print("=" * 72)
print("""
The SU(3)_c branching approach (Hypothesis 6 from 14_selection_rule_deep_think.py)
does NOT cleanly give the sector D_N values via any simple "filter".

The exact dimensions 28, 35, 63 appear to come from FULL Λ^k subspaces (not
filtered by SU(3)_c content).  The selection rule operates at the
"which Cl-grade" level, not at the "which sub-representation within a
grade" level.

This means the selection rule is more like:
  Sector N → choose specific SET of grades from {Λ⁰, Λ², Λ⁴, Λ⁶}

with the choice:
  N=0: {Λ², Λ⁶}
  N=1: {Λ⁴}
  N=2: {Λ², Λ⁴, Λ⁶}

This RULE itself doesn't decompose into simpler statements about SU(3)
content or G₂ invariants.  It's a discrete choice that needs an
underlying mechanism.

CANDIDATES for the discrete choice rule (none yet proved):
  1. A Lagrangian whose Yukawa term has explicit sector-dependent
     Clifford-grade projectors.
  2. A symmetry-breaking ladder Spin(8) → Spin(7) → G₂ → SU(3) that
     produces sector-specific "left-over" structure.
  3. A topological invariant (Chern class or similar) that's discrete
     and varies with Furey N.
  4. An UNDISCOVERED structural fact about the Furey ℂ⊗ℍ⊗𝕆 algebra
     that gives a grade-selecting role to N.

I cannot derive the selection rule in this session.  The empirical
identification of Cl(7)_even = Cl(6) = ℂ⊗𝕆 as the SINGLE PARENT
ALGEBRA stands as a real structural finding — but the SELECTION
RULE within that parent is not yet pinned to a mechanism.

The user's skepticism was right.
""")
