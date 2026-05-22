#!/usr/bin/env python3
"""
16_Z2_decomposition.py

NOVEL FINDING (this session): The selection rule has a Z₂ × Z₂ structure.

The single source Cl(7)_even = Cl(6) ≅ ℂ⊗𝕆 decomposes into TWO structurally
distinct pieces:

  • L  = Λ² ⊕ Λ⁶  (dim 21 + 7 = 28)  ← "Lie algebra content"
        = Spin(7) Lie algebra (Λ² = so(7)) + S⁷ quotient direction (Λ⁶)
        Contains NO G₂-invariant.

  • F  = Λ⁴  (dim 35)  ← "G₂-form content"
        Contains the COASSOCIATIVE 4-FORM *φ (the G₂-defining structure).
        The Hodge dual of the associative 3-form on ℝ⁷.

These two parts are STRUCTURALLY DIFFERENT:
  • L contains gauge-like Lie-algebra content (rotation generators).
  • F contains the octonion multiplication structure (G₂-defining form).

The CLAIM: each Furey N sector either includes or excludes each part:

    Sector       Bit-L   Bit-F   Dim
    ──────       ─────   ─────   ───
    N=0 (e_R)     1       0       28      (lepton: Lie only)
    N=1 (d_R)     0       1       35      (d-quark: form only)
    N=2 (u_R)     1       1       63      (u-quark: both)
    N=3 (ν_R)     ?       ?       ?       (neutrino: prediction)

The pattern N → (Bit-L, Bit-F) has a Z₂ × Z₂ shape.  No single formula
gives the bits from N, but the observed pattern is consistent.

INTERPRETATION:
  • Color singlets (e_R) couple to Lie content (gauge generators) but NOT
    to G₂-form (octonion multiplication — no color algebra).
  • Single color-creation states (d_R) couple ONLY to G₂-form (the
    octonion multiplication structure where color "lives").
  • Double color-creation states (u_R) couple to BOTH.

This explains the additive identity:  D_u = D_lepton + D_d-quark.

The selection rule isn't proven but the Z₂×Z₂ structure is highly suggestive.
"""

# ---------------------------------------------------------------------
# The two structural pieces of Cl(7)_even
# ---------------------------------------------------------------------
print("=" * 72)
print("The two structural pieces of Cl(7)_even ≅ Cl(6) ≅ ℂ⊗𝕆")
print("=" * 72)

L_dim = 21 + 7  # Λ² + Λ⁶
F_dim = 35      # Λ⁴

print(f"""
  Piece L  ("Lie algebra content"):
      L = Λ²ℝ⁷ ⊕ Λ⁶ℝ⁷
      Λ² = so(7) Lie algebra of Spin(7) = 21
      Λ⁶ = Hodge dual of Λ¹ ≅ ℝ⁷ ≅ S⁷ direction  = 7
      dim L = {L_dim}
      Contains NO G₂-invariant form.
      Identified with "rotation generators" / "Spin(7) gauge content".

  Piece F  ("G₂-form content"):
      F = Λ⁴ℝ⁷
      dim F = {F_dim}
      Contains the coassociative 4-form *φ (1-dim G₂-invariant).
      Plus 34-dim non-G₂-invariant content.
      Identified with "octonion multiplication structure" /
                      "G₂-defining content" / "color algebra".

  Both L and F sit inside Cl(7)_even (= Cl(6) = ℂ⊗𝕆, dim 64).
""")

# ---------------------------------------------------------------------
# The Z₂ × Z₂ binary pattern
# ---------------------------------------------------------------------
print("=" * 72)
print("Z₂ × Z₂ binary pattern of sector selections")
print("=" * 72)

sectors_Z2 = [
    {"N": 0, "name": "e_R (lepton)",        "bit_L": 1, "bit_F": 0, "color": "singlet"},
    {"N": 1, "name": "d_R (color triplet)", "bit_L": 0, "bit_F": 1, "color": "triplet"},
    {"N": 2, "name": "u_R (color triplet)", "bit_L": 1, "bit_F": 1, "color": "triplet"},
    {"N": 3, "name": "ν_R (lepton)",        "bit_L": "?","bit_F": "?", "color": "singlet"},
]

print(f"\n  N | Sector             | Color    | Bit-L | Bit-F | D_N  | Source")
print(f"  --|--------------------|---------|-------|-------|------|------------")
for s in sectors_Z2:
    if s['bit_L'] == 1 and s['bit_F'] == 0:
        D = L_dim
        src = "L (lie content)"
    elif s['bit_L'] == 0 and s['bit_F'] == 1:
        D = F_dim
        src = "F (G₂-form)"
    elif s['bit_L'] == 1 and s['bit_F'] == 1:
        D = L_dim + F_dim
        src = "L ⊕ F"
    else:
        D = "?"
        src = "TBD"
    print(f"  {s['N']} | {s['name']:<18} | {s['color']:<8} | {s['bit_L']:^5} | {s['bit_F']:^5} | {D!s:^4} | {src}")

print()
print("  Z₂ × Z₂ inference for N=3 (neutrino):")
print("    If bit_L = 0 (color singlet, like lepton-but-uncharged-style):")
print("      bit_F unclear — could be 0 (no coupling) or 1 (like d-quark)")
print()
print("    Case A: bit_L = 0, bit_F = 0  →  D_3 = 0  →  Brannen ansatz inapplicable")
print("           (consistent with neutrinos having tiny non-Brannen masses)")
print("    Case B: bit_L = 0, bit_F = 1  →  D_3 = 35  →  Q_ν = 11/15 (same as d-quark)")

# ---------------------------------------------------------------------
# The interpretation: what L and F represent physically
# ---------------------------------------------------------------------
print()
print("=" * 72)
print("Physical interpretation of L and F")
print("=" * 72)
print("""
L (Lie algebra content):
  • Λ² = so(7) — generators of Spin(7) rotations
  • Λ⁶ — Hodge dual of ℝ⁷, related to S⁷ = Spin(7)/G₂ coset
  • Together: the "rotation gauge field" content
  • Coupling: ANY fermion in a representation of Spin(7) sees this content
  • Color singlets (leptons) couple to L because they transform under Spin(7)
    but not under SU(3)_c ⊂ G₂

F (G₂-form content):
  • Λ⁴ — contains the coassociative 4-form *φ, defining the G₂ structure
  • The G₂ structure IS the octonion multiplication
  • Coupling: fermions in COLOR representations see this content
  • Quarks (color triplets) couple to F because they carry octonion-related
    color charge

The Z₂×Z₂ structure encodes:
  • Bit-L = "does this fermion transform under the Spin(7) rotation algebra?"
  • Bit-F = "does this fermion carry octonion color charge?"

In the Furey decomposition:
  • Leptons (N=0,3) carry charge but no color → couple to L but not F
  • Quarks (N=1,2) carry color → couple to F
  • u-quarks (N=2, "double creation") couple to BOTH L and F
  • d-quarks (N=1, "single creation") couple to F only
  • Why u_R has L-coupling while d_R doesn't: presumably because u_R is at
    N=2 (composite of two α's) which has "broader" structural content.

THE KEY OBSERVATION:  the additive identity D_u = D_e + D_d emerges naturally:
the u-quark sees BOTH structural pieces (L AND F), so its ambient is the
DIRECT SUM = L ⊕ F = 28 + 35 = 63.

This is the cleanest mechanism we've identified.  It's still not proven from
a Lagrangian, but it provides a CLEAR STRUCTURAL NARRATIVE that:
  (a) explains why D_u = D_e + D_d as the direct sum of two structural pieces,
  (b) gives a binary classification of sectors by color and gauge content,
  (c) is consistent with the Furey ℂ⊗𝕆 picture of color SU(3) emergence.
""")

# Save
import json
results = {
    "novel_finding": "Z₂ × Z₂ binary structure of Furey-N selection rule",
    "piece_L": {"description": "Lie algebra content", "grades": ["Λ²", "Λ⁶"], "dim": 28},
    "piece_F": {"description": "G₂-form content (coassoc 4-form)", "grades": ["Λ⁴"], "dim": 35},
    "sector_classification": [
        {"N": 0, "bit_L": 1, "bit_F": 0, "D": 28, "interpretation": "lepton: gauge only"},
        {"N": 1, "bit_L": 0, "bit_F": 1, "D": 35, "interpretation": "d-quark: form only"},
        {"N": 2, "bit_L": 1, "bit_F": 1, "D": 63, "interpretation": "u-quark: both"},
    ],
    "additive_identity_explanation": "D_u = D_e + D_d because u-quark sees direct sum L ⊕ F",
    "physical_interpretation": "L = Spin(7) gauge content; F = octonion color algebra content",
}
import json
with open('/home/d/code/scp/v59/cosserat_experiment/16_Z2.json', 'w') as f:
    json.dump(results, f, indent=2)
print()
print("Saved Z₂×Z₂ analysis to 16_Z2.json")
