#!/usr/bin/env python3
"""
v59/synthesis/octo_bridge_test.py

Operator-algebra checks in C⊗H⊗O, and what they say about the v_Higgs = 28²·a_l²
scale factor.

Two facts are established here numerically:

  (1) C⊗H⊗O is NOT alternative: (xx)y ≠ x(xy) in general.  So a Lagrangian written
      with literal algebraic squares like |Φ|² is ambiguous, and the natural
      well-defined replacement is the LEFT-REGULAR representation L_x(y) = x·y, whose
      operators compose associatively (composition of linear maps always does).

  (2) The associated trace form Tr(L_x L_x̃) = dim · |x|² is LINEAR in the dimension.
      A single operator trace therefore yields a factor of dim, i.e. √v ∝ √dim · a
      (the "equipartition" scaling) — NOT dim².

Conclusion: the operator-trace route does NOT explain the 28² factor.  The factor
28² = 784 is QUADRATIC in dim(L) and is a *component count* of the lepton mass
bilinear (a dim(L)×dim(L) Yukawa matrix), via its Frobenius² norm — the mechanism
pinned in integration_v58/03_higgs_bridge_result.md, where the 28 L-grade blades of
the real Cl(7) representation are shown to be Frobenius-orthogonal.  This script does
NOT reconstruct the true L = Λ²⊕Λ⁶ grade (the C⊗H⊗O basis here is not the Cl(7) blade
basis); for that, see bridge_vhiggs_cl7.py.

[History: an earlier version of this file claimed 28² emerged from a "bipartite vacuum
 polarization trace" Tr(P_L ⊗ P_L) of an arbitrary rank-28 projector.  That was removed:
 (a) the projector was the first 28 of 64 indices, not the L-grade; (b) Tr(P)² = dim² for
 ANY rank-28 projector and used none of the algebra; (c) a single vacuum-polarization loop
 is ONE trace (→ dim), so "bipartite loop = dim²" is internally inconsistent.  The honest
 result is the opposite: a single operator trace gives the LINEAR factor, below.]
"""

import numpy as np
import os

print("="*70)
print("Octospace operator-algebra test: squares and trace forms in C ⊗ H ⊗ O")
print("="*70)

# Load structure tensor
choh_path = os.path.join(os.path.dirname(__file__), "..", "furey_construction", "choh_structure.npz")
if not os.path.exists(choh_path):
    print(f"Error: Could not find {choh_path}")
    exit(1)

data = np.load(choh_path, allow_pickle=True)
C = data['struct_constants']  # shape (64, 64, 64): C[i, j, k] is k-th component of e_i * e_j
N_TOTAL = int(data['basis_decomp'].item()['N_TOTAL']) if 'basis_decomp' in data else 64


def CHO_mul(x, y):
    """Multiply x and y using the structure tensor: z_k = sum_ij x_i y_j C[i,j,k]."""
    return np.einsum('i,j,ijk->k', x, y, C)


def CHO_conj(x):
    """Tensor-product conjugation: conj(a⊗b⊗c) = ā⊗b̄⊗c̄. A basis element flips sign
    once per imaginary tensor factor (C: index a; H: index b; O: index c)."""
    result = np.zeros(N_TOTAL)
    for idx in range(N_TOTAL):
        a = idx // 32
        b = (idx % 32) // 8
        c = idx % 8
        sign_c = 1 if a == 0 else -1
        sign_h = 1 if b == 0 else -1
        sign_o = 1 if c == 0 else -1
        result[idx] = (sign_c * sign_h * sign_o) * x[idx]
    return result


def L_op(x):
    """64x64 matrix of left-multiplication L_x(y) = x·y, i.e. L_x[k,j] = sum_i x_i C[i,j,k]."""
    return np.einsum('i,ijk->kj', x, C)


print(f"\nLoaded C⊗H⊗O structure tensor: {C.shape}")

# ---------------------------------------------------------------------------
# Test 1: Alternative law  (xx)y = x(xy)
# ---------------------------------------------------------------------------
print("\n--- Test 1: Alternative (left) law  (xx)y = x(xy) ---")
np.random.seed(42)
violations = 0
max_diff = 0.0
for _ in range(10):
    x = np.random.randn(N_TOTAL)
    y = np.random.randn(N_TOTAL)
    lhs = CHO_mul(CHO_mul(x, x), y)   # (xx)y
    rhs = CHO_mul(x, CHO_mul(x, y))   # x(xy)
    diff = np.linalg.norm(lhs - rhs)
    max_diff = max(max_diff, diff)
    if diff > 1e-8:
        violations += 1

print(f"Violations of (xx)y = x(xy): {violations}/10  (max ||lhs-rhs|| = {max_diff:.3g})")
alternative = (violations == 0)
print(f"=> C⊗H⊗O is {'ALTERNATIVE' if alternative else 'NOT alternative'}.")

# ---------------------------------------------------------------------------
# Test 2: Left-regular operators and the trace form
# ---------------------------------------------------------------------------
print("\n--- Test 2: Left-regular operators L_x and the trace form Tr(L_x L_x̃) ---")

# Linear trace: Tr(L_x) = dim * (scalar part of x).
x = np.random.randn(N_TOTAL)
Lx = L_op(x)
trace_Lx = np.trace(Lx)
print(f"Tr(L_x) = {trace_Lx:.6f},   dim * x_0 = {N_TOTAL * x[0]:.6f}")
assert np.isclose(trace_Lx, N_TOTAL * x[0]), "Tr(L_x) should be dim * x_0"

# The norm trace form: Tr(L_x L_{conj(x)}) — the well-defined operator replacement for |x|².
ratios = []
for _ in range(5):
    z = np.random.randn(N_TOTAL)
    Lz = L_op(z)
    Lzbar = L_op(CHO_conj(z))
    tf = np.trace(Lz @ Lzbar)
    ratios.append(tf / (z @ z))      # tf / |z|^2
ratios = np.array(ratios)
print(f"Tr(L_x L_x̃) / |x|^2  over 5 random x: {ratios}")
print(f"=> trace form = ({ratios.mean():.3f}) * |x|^2   (= dim = {N_TOTAL})")
assert np.allclose(ratios, N_TOTAL, atol=1e-6), "trace form should equal dim * |x|^2"

# ---------------------------------------------------------------------------
# Test 3: What this says about the 28^2 scale factor
# ---------------------------------------------------------------------------
print("\n--- Test 3: Does the operator trace explain v_Higgs = 28^2 * a_l^2? ---")
dim_L = 28
print(f"A single operator trace is LINEAR in dimension:")
print(f"  restricting the trace form to a d={dim_L}-dim subspace gives a factor ~ d = {dim_L},")
print(f"  i.e. v ~ d * a^2  =>  sqrt(v) ~ sqrt(d) * a = sqrt({dim_L}) * a = {np.sqrt(dim_L):.3f} * a.")
print(f"  This is the EQUIPARTITION scaling, which 03_higgs_bridge_result.md EXCLUDES")
print(f"  (empirically sqrt(v)/a ~ 28, not sqrt(28) ~ 5.29).")
print()
print(f"The 28^2 = {dim_L**2} factor is QUADRATIC in d. It is the number of components of a")
print(f"  d x d lepton mass bilinear (Yukawa matrix Y), via its Frobenius^2 norm:")
print(f"  ||Y||_F^2 = (sum over d^2 entries) a^2 = d^2 * a^2  =>  sqrt(v) = d * a = {dim_L} * a.")
print(f"  Frobenius^2 of a matrix is NOT a single operator trace; the routes give different")
print(f"  powers of d ({dim_L} vs {dim_L**2}). The 28^2 belongs to the bilinear count, not the trace form.")

# ---------------------------------------------------------------------------
print("\n" + "="*70)
print("Summary (data-driven):")
print(f"1. C⊗H⊗O alternative law (xx)y=x(xy): {'holds' if alternative else 'FAILS'} "
      f"=> use left-regular operators L_x for a well-defined |Φ|².")
print(f"2. Operator trace form Tr(L_x L_x̃) = dim*|x|^2 = {N_TOTAL}*|x|^2  (LINEAR in dim).")
print(f"3. => a single operator trace gives sqrt(d) scaling (excluded). The 28^2=784 factor")
print(f"   is the component count of the lepton mass bilinear (Frobenius^2), per")
print(f"   integration_v58/03_higgs_bridge_result.md — NOT an operator trace.")
print("="*70)
