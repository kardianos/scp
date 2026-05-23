#!/usr/bin/env python3
"""
v59/algebra/brannen_kernel.py

Brannen kernel M_X(ξ_X) realized algebraically with ξ_X embedded into the
appropriate Cl(7)_even subspace for each sector.

ξ_X ∈ ℍ ≅ Cl(0,2) ⊂ Cl(7)_even where the embedding selects a 4-dim ℍ-slice:
  ℍ = ⟨1, i, j, k⟩  →  ⟨e_∅, e_a, e_b, e_c⟩
where (e_a, e_b, e_c) are 3 mutually-orthogonal bivectors (Λ² elements).

For each sector, the bivector choice (e_a, e_b, e_c) lives in the sector's
ambient subspace:
  Lepton (Bit-L=1):   bivectors from L = Λ²⊕Λ⁶  (specifically Λ² part)
  d-quark (Bit-F=1):  ξ embedded via Λ⁴ via Cl(3) ≅ ℍ identification
  u-quark (Bit-L=1, Bit-F=1):  bivectors from L⊕F

The Brannen mass kernel M_X = a_X(I + ξ_X S + ξ̄_X S²) acts on a 3-flavor
generation space (Z_3 triality).
"""
import numpy as np
import json
import os
import math
from itertools import combinations

# Load Cl(7)_even structure
struct_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           "cl7_even_structure.json")
with open(struct_path) as f:
    struct = json.load(f)

n_basis = struct["n_basis"]
GRADES = struct["grades"]
BASIS = [tuple(b) for b in struct["basis_tuples"]]
BASIS_IDX = {b: i for i, b in enumerate(BASIS)}
L_indices = struct["L_indices"]
F_indices = struct["F_indices"]
L2_indices = struct["L2_indices"]
L6_indices = struct["L6_indices"]
# Reconstruct multiplication table
mult_data = struct["mult_table"]
MULT_TABLE = {}
for idx, (sign, k) in enumerate(mult_data):
    i, j = idx // 64, idx % 64
    MULT_TABLE[(i, j)] = (sign, k)


# =============================================================================
# Part 1: Embed ℍ into Cl(7)_even via choice of 3 bivectors
# =============================================================================

def quaternion_embedding(bivec_a, bivec_b, bivec_c):
    """
    Embed ℍ ↪ Cl(7)_even by mapping
      1 → e_∅,
      i → e_{bivec_a},
      j → e_{bivec_b},
      k → e_{bivec_c}.

    For this to be a quaternion algebra, the three bivectors must
    satisfy:  e_a · e_b = e_c (cyclically) and e_a² = -1.

    Returns (e_idx_1, e_idx_i, e_idx_j, e_idx_k) — 4 indices in BASIS.
    """
    e_1 = BASIS_IDX[()]
    e_i = BASIS_IDX[tuple(sorted(bivec_a))]
    e_j = BASIS_IDX[tuple(sorted(bivec_b))]
    e_k = BASIS_IDX[tuple(sorted(bivec_c))]

    # Check quaternion relations
    sign_ij, prod_ij = MULT_TABLE[(e_i, e_j)]
    # i·j should = k or -k
    if prod_ij != e_k:
        return None  # not a quaternion triple
    return (e_1, e_i, e_j, e_k), sign_ij


def find_quaternion_triples_in_L():
    """Find all (e_a, e_b, e_c) quaternion triples where e_a, e_b, e_c are bivectors in L."""
    triples = []
    # Try all (a, b) pairs of 2-element subsets in L2 (bivectors)
    bivec_subsets = [BASIS[idx] for idx in L2_indices]
    for i, sa in enumerate(bivec_subsets):
        for j, sb in enumerate(bivec_subsets):
            if i >= j:
                continue
            # i·j gives a Cl product
            idx_a = BASIS_IDX[sa]
            idx_b = BASIS_IDX[sb]
            sign_ab, prod_ab = MULT_TABLE[(idx_a, idx_b)]
            if GRADES[prod_ab] != 2:
                continue  # not a bivector
            sc = BASIS[prod_ab]
            triple, sign = quaternion_embedding(sa, sb, sc)
            if triple is not None:
                # Check sign relation: i·j = +k (canonical) or -k
                triples.append({
                    'bivecs': (sa, sb, sc),
                    'indices': triple,
                    'sign_ij': sign_ab,
                })
    return triples


print("=" * 70)
print("Quaternion ℍ embeddings into L = Λ²⊕Λ⁶ (lepton ambient)")
print("=" * 70)

triples = find_quaternion_triples_in_L()
print(f"Found {len(triples)} candidate quaternion triples in Λ² ⊂ L")

# Show first few
for t in triples[:5]:
    sa, sb, sc = t['bivecs']
    print(f"  e_{sa} · e_{sb} = {t['sign_ij']:+d}·e_{sc}")


# =============================================================================
# Part 2: The CANONICAL lepton-sector ℍ-slice
# =============================================================================
# We pick a specific quaternion triple for the lepton sector.
# Standard choice: identify ℍ with (spacetime) bivector algebra,
# e.g., (e_01, e_02, e_12) → (i, j, k) using indices {0, 1, 2}.

# Check: e_01 · e_02 = -e_12 (or +e_12, depending on sign convention)
e01 = BASIS_IDX[(0, 1)]
e02 = BASIS_IDX[(0, 2)]
e12 = BASIS_IDX[(1, 2)]
s, k = MULT_TABLE[(e01, e02)]
print(f"\nCanonical lepton ℍ-slice: bivectors (e_01, e_02, e_12)")
print(f"  e_01 · e_02 = {s:+d}·e_{BASIS[k]} = {s:+d}·e_12")
if k == e12 and s == -1:
    print("  → e_01 · e_02 = -e_12 (this gives a quaternion algebra with reversed orientation)")
    # Equivalent: use (i, j, k) ↔ (e_01, e_02, -e_12)
    LEPTON_QUATERNION_TRIPLE = (e01, e02, e12)
    LEPTON_QUATERNION_SIGNS = (1, 1, -1)  # i, j, k correspond to e_01, e_02, -e_12
else:
    print("  Check sign convention — may need adjustment.")

# Build the ℍ-slice basis indices for the lepton sector
e_empty = BASIS_IDX[()]
LEPTON_SLICE_BASIS = [e_empty, e01, e02, e12]
print(f"Lepton ℍ-slice basis indices in Cl(7)_even: {LEPTON_SLICE_BASIS}")
print(f"  → {[BASIS[idx] for idx in LEPTON_SLICE_BASIS]}")


# =============================================================================
# Part 3: Construct ξ as a Cl(7)_even element
# =============================================================================
def xi_to_cl7_element(xi_components, slice_indices, slice_signs=(1, 1, 1, 1)):
    """
    Given quaternion components ξ = (ξ_0, ξ_1, ξ_2, ξ_3) and
    a 4-index slice (e_1_idx, e_i_idx, e_j_idx, e_k_idx) into Cl(7)_even,
    construct ξ as a 64-dim vector.
    """
    v = np.zeros(64)
    for c, idx, sgn in zip(xi_components, slice_indices, slice_signs):
        v[idx] += c * sgn
    return v


# Sector-specific Brannen quaternion ξ_X
# t_X² = 1 - 14/D_X gives |ξ_X|² = t_X²
# φ_X: Brannen phase per sector
# We use (ξ_0, ξ_1, 0, 0) form for the COMPLEX slice (i, j, k all zero except i)

t_l_sq = 1/2
phi_l = 2/9         # structural EXACT
xi_l_quat = (math.sqrt(t_l_sq) * math.cos(phi_l),
             math.sqrt(t_l_sq) * math.sin(phi_l),
             0.0,
             0.0)
xi_l_cl7 = xi_to_cl7_element(xi_l_quat, LEPTON_SLICE_BASIS,
                              (1, 1, 1, -1))  # i = e_01, j = e_02, k = -e_12

print(f"\nLepton ξ_l: |ξ_l|² = {sum(c**2 for c in xi_l_quat):.4f} (should be 1/2)")
print(f"  Cl(7)_even components: {xi_l_cl7[:5]} ... (nonzero at indices {[i for i, v in enumerate(xi_l_cl7) if abs(v) > 1e-10]})")


# =============================================================================
# Part 4: The Brannen mass kernel as a 3×3 matrix on flavor
# =============================================================================
# M(ξ) = a·(I + ξ·S + ξ̄·S²) on 3-flavor (generation) space.
# S is the cyclic shift, a 3×3 matrix:
S = np.array([[0, 1, 0],
              [0, 0, 1],
              [1, 0, 0]], dtype=complex)

def brannen_kernel_complex(a, t, phi):
    """3×3 Brannen kernel with ξ = t·e^(iφ) (complex slice)."""
    xi = t * np.exp(1j * phi)
    return a * (np.eye(3, dtype=complex) + xi * S + np.conj(xi) * S.conj().T)


# Verify: lepton kernel reproduces the lepton mass triplet
a_l = 17.7156  # √MeV
t_l = math.sqrt(1/2)
M_l = brannen_kernel_complex(a_l, t_l, phi_l)
M_l_H = (M_l + M_l.conj().T) / 2
eigs = np.linalg.eigvalsh(M_l_H)
masses = eigs**2
print()
print("=" * 70)
print("Sanity: Brannen lepton kernel reproduces (m_e, m_μ, m_τ)?")
print("=" * 70)
print(f"  Brannen amplitudes: {sorted(eigs)}")
print(f"  Predicted masses: {sorted(masses)}  (expected: 0.511, 105.66, 1776.86)")
Q = sum(masses) / sum(np.sqrt(masses))**2
print(f"  Koide Q = {Q:.6f}  (expected 2/3 = {2/3:.6f}, EXACT)")


# =============================================================================
# Part 5: SECTOR-SPECIFIC ℍ-slice embeddings
# =============================================================================
# Lepton ℍ-slice: in L (Λ²) — three bivectors picked from L₂.
# d-quark ℍ-slice: in F (Λ⁴) — need three "quaternion-like" elements in Λ⁴.
# u-quark: ℍ-slice in L⊕F (combination).

# For d-quark: Λ⁴ has dim 35 = (7 choose 4).  We need a 4-dim ℍ-slice in it.
# The simplest: use Λ⁴ elements that span a quaternionic substructure.
# Specifically, Λ⁴ ⊃ {e_∅·, e_a, e_b, e_c·} won't work since e_∅ is not in Λ⁴.
# Better: use 4-element subsets that close under multiplication modulo Λ⁴.
# Hmm — squares of 4-element subsets in Cl(7,0) are scalars.  Products of two
# 4-element subsets are EVEN grade — possibly in Λ⁰ or Λ⁴ or Λ⁶ depending on overlap.

# Let me find quaternion structures within Λ⁴.
def find_quaternion_in_F():
    """Find three Λ⁴ elements (e_a, e_b, e_c) satisfying e_a·e_b = ±e_c."""
    Lambda4 = [BASIS[idx] for idx in F_indices]
    triples = []
    for i in range(len(Lambda4)):
        sa = Lambda4[i]
        idx_a = BASIS_IDX[sa]
        # e_a · e_a = ±e_∅ (scalar).
        # We need an "imaginary" direction inside Λ⁴.
        # An element X ∈ Λ⁴ satisfies X² = ±1.  If X² = -1, X plays the role of i.
        s, k = MULT_TABLE[(idx_a, idx_a)]
        if k == BASIS_IDX[()] and s == -1:
            # X² = -1 — candidate for "i" in ℍ
            pass

    return triples

# A simpler statement: Λ⁴ doesn't have an internal ℍ algebra in the
# canonical sense (it's not closed under multiplication).  So we identify
# the d-quark ℍ-slice via a SECTOR-SPECIFIC subset of Λ⁴.

# For now, place d-quark ξ_d on a 4-dim subspace of Λ⁴ using a single
# representative quartic + its "rotated" counterparts.
# Standard convention: the F = Λ⁴ sector has the coassociative 4-form *φ
# as G_2-invariant element.  This is a specific 4-element subset (or
# linear combination).  Define the *φ-aligned ℍ-slice.

# For simplicity here: PLACEHOLDER — use 4 specific Λ⁴ elements as the
# d-quark ℍ-slice basis.  The PHYSICAL choice would be the G_2-equivariant
# embedding.
print()
print("=" * 70)
print("d-quark ℍ-slice (placeholder; needs G_2-equivariant refinement)")
print("=" * 70)
# Pick four Λ⁴ elements: e_1234, e_1256, e_3456, e_1357 (placeholder)
d_quark_slice = [
    BASIS_IDX[(0, 1, 2, 3)],
    BASIS_IDX[(0, 1, 4, 5)],
    BASIS_IDX[(2, 3, 4, 5)],
    BASIS_IDX[(0, 2, 4, 6)],
]
print(f"  d-quark slice basis: {[BASIS[idx] for idx in d_quark_slice]}")
print(f"  (PLACEHOLDER — replace with G_2-equivariant choice for production)")

# u-quark slice: spans L⊕F. Combination of lepton + d-quark slices.
u_quark_slice = LEPTON_SLICE_BASIS + d_quark_slice  # 8 elements
print(f"\nu-quark slice basis (L+F): {len(u_quark_slice)} elements")


# =============================================================================
# Part 6: Save sector embeddings
# =============================================================================
sector_data = {
    "lepton_slice": LEPTON_SLICE_BASIS,
    "lepton_slice_signs": [1, 1, 1, -1],
    "lepton_slice_basis_tuples": [list(BASIS[idx]) for idx in LEPTON_SLICE_BASIS],
    "lepton_brannen": {
        "a_l_sqrtMeV": a_l,
        "t_l_sq": 1/2,
        "phi_l": 2/9,
        "predicted_masses": sorted(masses.tolist()),
        "Q_l": Q,
    },
    "d_quark_slice_placeholder": d_quark_slice,
    "u_quark_slice_placeholder": u_quark_slice,
    "notes": "Sector ℍ-slice embeddings.  Lepton is canonical (e_01, e_02, e_12). "
             "d-quark and u-quark placeholders need G_2-equivariant refinement "
             "based on the coassociative 4-form *φ structure.",
}

out_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "brannen_embeddings.json")
with open(out_path, "w") as f:
    json.dump(sector_data, f, indent=2)
print(f"\nSector embedding data saved to {out_path}")


# =============================================================================
# Part 7: Hooks for the 1-loop V_eff calculation
# =============================================================================
print()
print("=" * 70)
print("Hooks for 1-loop V_eff (next: v_eff_loop.py)")
print("=" * 70)
print("""
For the 1-loop V_eff(φ_X) calculation:
  - We have the Brannen kernel M(ξ_X) per sector
  - We have the sector ℍ-slice in Cl(7)_even
  - The loop integral becomes a FINITE SUM over the sector subspace
    (dim 28 for lepton, dim 35 for d-quark, dim 63 for u-quark)
  - The coefficient of cos(3φ_X) is the structural integer N_X
  - The COUPLING is α at the sector's natural dim-density scale

Next file `v_eff_loop.py` will:
  1. Compute Tr[M(ξ)^n · log M(ξ)] expansions
  2. Extract the coefficient of cos(3φ) at each loop order
  3. Verify the predicted N_X = dim G_2 (d) or dim Spin(5) (u)
""")
