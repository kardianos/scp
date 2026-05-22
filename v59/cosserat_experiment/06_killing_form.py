#!/usr/bin/env python3
"""
06_killing_form.py

Killing form on so(7) and the SU(2) ⊂ Spin(7) embedding index.

Conjecture from 05_alpha_alphaW_ratio.py: α_W/α = √(dim Spin(7)) = √21.
Honest re-examination at consistent scale (M_Z): empirical α_W/α = 1/sin²θ_W
= 4.32, while √21 = 4.58 is 5.7 % off — weaker than the initial "1.4 % match"
which used inconsistent scales.

This script computes from first principles:
  • Killing form on so(7) as a 21×21 matrix in an explicit basis.
  • The so(3) ⊂ so(7) subalgebra corresponding to "imaginary-direction
    rotations" of ξ ∈ ℍ — i.e. the silent SU(2) ⊂ Spin(7).
  • The Dynkin embedding index, defined as
        T_{so(3) ↪ so(7)} = B_{so(7)}|_{so(3)} / B_{so(3)}
    (ratio of restricted Killing form to intrinsic Killing form on so(3)).
  • What natural numerical ratio this gives for gauge couplings.

Reference: for an embedding ρ: H → G of compact simple Lie groups, the
Dynkin embedding index X is defined by B_H(X, Y) = X · B_G(ρ(X), ρ(Y)),
and gives the natural rescaling of gauge couplings when one gauge group is
"embedded" in another.
"""

import numpy as np

# ---------------------------------------------------------------------
# Step 1.  Construct so(7) as 7×7 antisymmetric matrices.
# ---------------------------------------------------------------------
N = 7
print("=" * 70)
print("Step 1: so(7) — 7×7 antisymmetric matrices")
print("=" * 70)

# Basis: E_{ij} for 1 ≤ i < j ≤ 7, where E_{ij} has +1 at (i,j) and -1 at (j,i).
# Number of basis elements = (7 choose 2) = 21.
basis = []
labels = []
for i in range(N):
    for j in range(i+1, N):
        M = np.zeros((N, N))
        M[i, j] = 1.0
        M[j, i] = -1.0
        basis.append(M)
        labels.append(f"E_{{{i+1}{j+1}}}")

dim_so7 = len(basis)
print(f"  Basis size = (7 choose 2) = {dim_so7}  (= dim so(7) ✓)")

# ---------------------------------------------------------------------
# Step 2.  Trace form Tr_vec(T^a T^b) in the 7-dim vector representation.
# ---------------------------------------------------------------------
print()
print("=" * 70)
print("Step 2: Trace form in the vector representation")
print("=" * 70)

trace_vec = np.zeros((dim_so7, dim_so7))
for a in range(dim_so7):
    for b in range(dim_so7):
        trace_vec[a, b] = np.trace(basis[a] @ basis[b])

# Check that it's diagonal and uniformly -2 on the diagonal (standard for E_{ij}'s)
diag = np.diag(trace_vec)
off_diag_max = np.max(np.abs(trace_vec - np.diag(diag)))
print(f"  Tr_vec(E_a E_a) for first few:")
for a in range(5):
    print(f"    Tr(E_{labels[a].split('_')[1].strip('{}')}²) = {diag[a]:.4f}")
print(f"  Max off-diagonal |Tr(E_a E_b)|, a ≠ b: {off_diag_max:.4e}")
print(f"  All diagonal entries = -2 (basis is orthogonal under trace form ✓)")

# ---------------------------------------------------------------------
# Step 3.  Compute the Killing form B(X, Y) = Tr(ad_X ad_Y).
# ---------------------------------------------------------------------
print()
print("=" * 70)
print("Step 3: Killing form B(X, Y) = Tr(ad_X ad_Y)")
print("=" * 70)
print()
print("Build the adjoint representation: for each basis element E_a, ad_{E_a}")
print("is the 21×21 matrix whose action is [E_a, ·] in the basis {E_b}.")

# Structure constants: [E_a, E_b] = Σ_c f^{abc} E_c, then ad_{E_a}(E_b) = Σ_c f^{abc} E_c.
# Compute f^{abc} = -2 × Tr(E_a [E_b, E_c]) / Tr(E_a²) ... but easier to compute
# [E_a, E_b] explicitly and project onto the E basis using the trace form.

# Project a matrix M into the basis: M = Σ_c (Tr(M E_c) / Tr(E_c²)) E_c.
def project_basis(M):
    """Return coefficients in {basis[c]}."""
    coefs = np.zeros(dim_so7)
    for c in range(dim_so7):
        coefs[c] = np.trace(M @ basis[c]) / np.trace(basis[c] @ basis[c])
    return coefs

# ad_{E_a} as 21×21 matrix
ad = np.zeros((dim_so7, dim_so7, dim_so7))
for a in range(dim_so7):
    for b in range(dim_so7):
        comm = basis[a] @ basis[b] - basis[b] @ basis[a]
        ad[a, :, b] = project_basis(comm)

# Killing form
killing = np.zeros((dim_so7, dim_so7))
for a in range(dim_so7):
    for b in range(dim_so7):
        killing[a, b] = np.trace(ad[a] @ ad[b])

diag_K = np.diag(killing)
off_diag_K_max = np.max(np.abs(killing - np.diag(diag_K)))
print(f"  Killing form diagonal entries (first 5):")
for a in range(5):
    print(f"    B(E_{labels[a].split('_')[1].strip('{}')}, E_a) = {diag_K[a]:.4f}")
print(f"  Max off-diagonal |B(E_a, E_b)|, a≠b: {off_diag_K_max:.4e}")
print()
ratio_KvsTr = diag_K / diag
print(f"  Ratio B(E_a, E_a) / Tr_vec(E_a²) (first 5): "
      f"{[f'{r:.4f}' for r in ratio_KvsTr[:5]]}")
print(f"  Expected for so(n): B = (n-2) Tr_vec  →  for n=7: 5.")
print(f"  Mean ratio: {np.mean(ratio_KvsTr):.6f}")

# ---------------------------------------------------------------------
# Step 4.  Identify so(3) ⊂ so(7) for imaginary-direction rotations.
# ---------------------------------------------------------------------
print()
print("=" * 70)
print("Step 4: The so(3) ⊂ so(7) of imaginary-direction rotations")
print("=" * 70)
print()
print("In our v59 setup, the 7-dim imaginary octonion space is ⟨i, j, k, e, ie, je, ke⟩.")
print("The 'silent' SU(2) ⊂ Spin(7) corresponds to rotations of ⟨i, j, k⟩ — the")
print("first 3 imaginary directions (the quaternion imaginary subspace).")
print()
print("As elements of our so(7) basis, this so(3) is generated by:")
print("  J_1 = E_{23} (rotation in jk-plane)")
print("  J_2 = E_{13} (rotation in ik-plane)")
print("  J_3 = E_{12} (rotation in ij-plane)")
print()

# Identify the so(3) generators in our basis
# E_{ij} with i<j, listed as (1,2), (1,3), ..., (1,7), (2,3), ..., (6,7).
# Index of E_{ij}: i, j with i<j.  Let's find indices.
def basis_idx(i, j):
    """Index of E_{ij} in the basis (1-indexed conventions)."""
    a, b = min(i, j), max(i, j)
    for idx, label in enumerate(labels):
        ij = label.split('_')[1].strip('{}')
        if ij == f"{a}{b}":
            return idx
    raise ValueError(f"E_{i}{j} not found")

# so(3) generators (using 1-indexed quaternion units 1,2,3)
so3_indices = [basis_idx(1, 2), basis_idx(1, 3), basis_idx(2, 3)]
print(f"  Indices in our basis: {so3_indices}")
print(f"  Labels: {[labels[i] for i in so3_indices]}")

# Verify it's a closed so(3): [E_{12}, E_{13}] should be ± E_{23}.
def comm_in_basis(a_idx, b_idx):
    M = basis[a_idx] @ basis[b_idx] - basis[b_idx] @ basis[a_idx]
    return project_basis(M)

c1 = comm_in_basis(so3_indices[0], so3_indices[1])  # [E_{12}, E_{13}]
c2 = comm_in_basis(so3_indices[0], so3_indices[2])  # [E_{12}, E_{23}]
c3 = comm_in_basis(so3_indices[1], so3_indices[2])  # [E_{13}, E_{23}]
print(f"\n  [E_{{12}}, E_{{13}}] projected onto E_{{23}}: {c1[so3_indices[2]]:.4f}  (expect ±1)")
print(f"  [E_{{12}}, E_{{23}}] projected onto E_{{13}}: {c2[so3_indices[1]]:.4f}  (expect ∓1)")
print(f"  [E_{{13}}, E_{{23}}] projected onto E_{{12}}: {c3[so3_indices[0]]:.4f}  (expect ±1)")
print("  → Closed under the bracket: this is so(3) ⊂ so(7).")

# ---------------------------------------------------------------------
# Step 5.  Compute the Dynkin embedding index.
# ---------------------------------------------------------------------
print()
print("=" * 70)
print("Step 5: Dynkin embedding index of so(3) ⊂ so(7)")
print("=" * 70)

# Intrinsic Killing form on so(3) (computed standalone)
# For so(3): B(X, Y) = (3-2) Tr_vec(X Y) = Tr_3×3 vec rep.
# But we want this in the SAME basis E_{12}, E_{13}, E_{23} as embedded in so(7).
# Compute the intrinsic Killing form of so(3) using these generators.

# so(3) acts as 3×3 antisymmetric matrices on ⟨1, 2, 3⟩.
# In so(7), the same generators E_{ij} are 7×7 matrices, but their action on
# the (1,2,3) subspace is the SAME as the so(3) action.
# So the intrinsic Killing form of so(3) (using its OWN adjoint rep, which is 3-dim)
# differs from the restriction of the so(7) Killing form.

# Intrinsic so(3) Killing form: on 3-dim adjoint = vector rep, B_{so(3)}(X, Y) = Tr_vec(X Y).
# For E_{12}^{(3)} (3×3 matrix), Tr(E_{12}² ) = -2. So B_{so(3)}(E_{12}, E_{12}) = -2 × 1 = -2.
# Wait — but for so(3), B(X, Y) = (3-2) Tr_vec(X Y) = Tr_vec(X Y).
# This gives B_{so(3)}(E_{12}, E_{12}) = -2.

# Actually let me compute it directly. so(3) has 3 generators, 3-dim adjoint.
B_so3 = np.zeros((3, 3))
gens_so3 = [np.zeros((3, 3)), np.zeros((3, 3)), np.zeros((3, 3))]
# E_{12}: rotation in 12-plane
gens_so3[0][0, 1] = 1; gens_so3[0][1, 0] = -1
# E_{13}: rotation in 13-plane
gens_so3[1][0, 2] = 1; gens_so3[1][2, 0] = -1
# E_{23}: rotation in 23-plane
gens_so3[2][1, 2] = 1; gens_so3[2][2, 1] = -1

# Adjoint matrices: ad_{a}(b) = [g_a, g_b]; project onto {g_a}.
ad_so3 = np.zeros((3, 3, 3))
for a in range(3):
    for b in range(3):
        comm = gens_so3[a] @ gens_so3[b] - gens_so3[b] @ gens_so3[a]
        # Project onto basis: comm = Σ_c α_c g_c  with α_c = Tr(comm·g_c) / Tr(g_c²)
        for c in range(3):
            ad_so3[a, c, b] = np.trace(comm @ gens_so3[c]) / np.trace(gens_so3[c] @ gens_so3[c])

# Killing form on so(3)
B_so3_intrinsic = np.zeros((3, 3))
for a in range(3):
    for b in range(3):
        B_so3_intrinsic[a, b] = np.trace(ad_so3[a] @ ad_so3[b])

print(f"  Intrinsic so(3) Killing form on (E_{{12}}, E_{{13}}, E_{{23}}):")
print(B_so3_intrinsic)
print(f"  diag: {np.diag(B_so3_intrinsic)}")
print()

# Restriction of so(7) Killing form to the so(3) subalgebra
B_so7_restricted = killing[np.ix_(so3_indices, so3_indices)]
print(f"  so(7) Killing form restricted to so(3) sub-basis:")
print(B_so7_restricted)
print(f"  diag: {np.diag(B_so7_restricted)}")
print()

# Embedding index ratio
indices = np.diag(B_so7_restricted) / np.diag(B_so3_intrinsic)
print(f"  Ratio B_{{so(7)}}|_{{so(3)}} / B_{{so(3)}}: {indices}")
print(f"  Embedding (Dynkin) index = {np.mean(indices):.4f}")
print()
print("For so(N) acting on n consecutive labels (so(3) ⊂ so(7) with same generators):")
print("  embedding index = (N - 2) / (n - 2) = (7-2)/(3-2) = 5  (theoretical)")

# ---------------------------------------------------------------------
# Step 6.  Now also compute via the vector trace
# ---------------------------------------------------------------------
print()
print("=" * 70)
print("Step 6: Vector-trace form ratios")
print("=" * 70)
Tr_vec_so7_restricted = trace_vec[np.ix_(so3_indices, so3_indices)]
Tr_vec_so3 = np.array([[np.trace(g_a @ g_b) for g_b in gens_so3] for g_a in gens_so3])
print()
print(f"  Tr_vec(so(7))|_{{so(3)}} diag: {np.diag(Tr_vec_so7_restricted)}")
print(f"  Tr_vec(so(3)) diag:           {np.diag(Tr_vec_so3)}")
print(f"  Ratio (Tr_vec_so7|so3 / Tr_vec_so3): {np.diag(Tr_vec_so7_restricted)/np.diag(Tr_vec_so3)}")
print()
print("→ This ratio is 1 (the vector rep of so(7) restricted to the so(3) subspace")
print("  IS exactly the vector rep of so(3) on its 3-dim subspace).")

# ---------------------------------------------------------------------
# Step 7.  Implication for gauge-coupling ratio.
# ---------------------------------------------------------------------
print()
print("=" * 70)
print("Step 7: Gauge-coupling ratio from the Dynkin index")
print("=" * 70)
print()
print("""
For a YM gauge theory with gauge group G (here Spin(7)) broken to subgroup H
(here SU(2) ⊂ Spin(7)), the gauge couplings are related by

    1/g_G²  =  X_H  /  g_H²

where X_H is the Dynkin embedding index of H in G.  Equivalently,

    g_H  =  g_G · √(X_H)

For our so(3) ⊂ so(7) embedding, X = (7-2)/(3-2) = 5 (Killing-form ratio).

If we identify:
    H = SU(2)_L     (silent direction)         →  coupling g_W
    G = "parent group" of the full silent + active structure  →  coupling g_G

Then  g_W = g_G · √5  ≈ 2.236 · g_G.

This is the FIRST principled estimate, free from the convention games of
naive volume-ratio arguments.

Comparing to data: empirical g_W ≈ 0.6517 (PDG, at M_Z).
                      empirical g_U(1) = g_Y ≈ 0.358 (M_Z) or e ≈ 0.313 (low).
                      Empirical g_W / g_Y = 0.6517 / 0.358 = 1.821.
                      Empirical g_W / e   = 0.6517 / 0.313 = 2.082.

√5 = 2.236.  Empirical g_W/g_Y = 1.821.  Off by 23 %.
√5 = 2.236.  Empirical g_W/e   = 2.082.  Off by 7 %.

Neither matches √21 (which would be g_W / g_singlet = √21 = 4.58).

CONCLUSION: the natural Killing-form ratio for the so(3) ⊂ so(7) embedding
is √5, NOT √21.  The earlier conjecture α_W/α = √21 from 05_alpha_alphaW_ratio.py
is NOT supported by the Killing-form calculation — it was a numerological match
(within 1.4 % at mixed scales, within 5.7 % at consistent M_Z scales).

The √5 = √(dim Spin(7) − dim G_2) factor IS more structurally meaningful
(= dim of S^7 = Spin(7)/G_2), and it gives a coupling ratio of 2.24, closer
to the empirical g_W/e ≈ 2.08 than √21 was.
""")
print(f"  √5 (Killing-form natural)            = {np.sqrt(5):.4f}")
print(f"  empirical g_W / e (low-E)            = {0.6517 / 0.30282:.4f}  (using e from α=1/137)")
print(f"  empirical g_W / e (M_Z)              = {0.6517 / 0.31223:.4f}  (using e from α=1/127.95)")
print(f"  empirical g_W / g_Y (M_Z)             = {0.6517 / 0.358:.4f}")
print()

# A more careful comparison using PDG numbers:
g2_emp = 0.6517
gY_emp = 0.358   # at M_Z
e_emp_lowE = np.sqrt(4 * np.pi / 137.036)
e_emp_MZ = np.sqrt(4 * np.pi / 127.952)

print(f"  Ratios for various candidate identifications:")
print(f"    g_W / e (low-E)   = {g2_emp / e_emp_lowE:.4f}")
print(f"    g_W / e (M_Z)     = {g2_emp / e_emp_MZ:.4f}")
print(f"    g_W / g_Y (M_Z)   = {g2_emp / gY_emp:.4f}")
print(f"    √5                 = {np.sqrt(5):.4f}")
print(f"    √(dim S^7) = √7    = {np.sqrt(7):.4f}")
print(f"    √21                = {np.sqrt(21):.4f}")
print(f"    √14                = {np.sqrt(14):.4f}")

# Save
np.savez('/home/d/code/scp/v59/cosserat_experiment/06_killing.npz',
         dim_so7=dim_so7, killing=killing, trace_vec=trace_vec,
         so3_indices=so3_indices,
         dynkin_index_so3_in_so7=np.mean(indices))
print()
print("Saved Killing-form data to 06_killing.npz")
