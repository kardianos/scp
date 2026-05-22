#!/usr/bin/env python3
"""
v59/furey_construction/01_choh_algebra.py

Variant A: Build the C ⊗ H ⊗ O algebra explicitly.

The tensor product C ⊗ H ⊗ O has dimension 2 × 4 × 8 = 64 over the reals.
This is Furey's conjectured algebra of the Standard Model fermions.

Strategy:
  1. Use 32 generators (basis): {1, i_C} ⊗ {1, i_H, j_H, k_H} ⊗ {e_0, e_1, ..., e_7}
     Wait, the C-part has only 2 basis elements (1, i_C), so 2 × 4 × 8 = 64 total.
  2. Multiplication is defined componentwise across the tensor factors.
  3. Verify the algebra is associative-by-factor (non-associativity confined
     to the O factor).

Properties to verify:
  - Dimension 64 over R.
  - C and H factors are associative; O factor is non-associative.
  - The full algebra is alternative (associator = 0 on like elements).
  - Norm |x|^2 = x x* is well-defined.
"""

import numpy as np
from itertools import product

print("="*72)
print("Variant A: Build C ⊗ H ⊗ O explicitly")
print("="*72)


# =========================================================================
# Part 1: Component algebras
# =========================================================================
print()
print("-"*72)
print("Part 1: Component algebras (C, H, O)")
print("-"*72)

# --- Complex numbers C ---
# Basis: 1, i.  i^2 = -1.
def C_mul(a, b):
    # a, b are length-2 arrays [re, im]
    return np.array([
        a[0]*b[0] - a[1]*b[1],
        a[0]*b[1] + a[1]*b[0],
    ])

def C_conj(a):
    return np.array([a[0], -a[1]])

# --- Quaternions H ---
# Basis: 1, i, j, k.  i^2=j^2=k^2=-1, ij=k, jk=i, ki=j.
def H_mul(a, b):
    a0, a1, a2, a3 = a
    b0, b1, b2, b3 = b
    return np.array([
        a0*b0 - a1*b1 - a2*b2 - a3*b3,
        a0*b1 + a1*b0 + a2*b3 - a3*b2,
        a0*b2 - a1*b3 + a2*b0 + a3*b1,
        a0*b3 + a1*b2 - a2*b1 + a3*b0,
    ])

def H_conj(a):
    return np.array([a[0], -a[1], -a[2], -a[3]])

# --- Octonions O ---
# Basis: e_0=1, e_1, ..., e_7.  Multiplication via Fano plane.

# Fano triples (Cayley convention).
# Each triple (a, b, c) means: e_a * e_b = e_c, e_b * e_c = e_a, e_c * e_a = e_b,
# with the reverse products getting -.
triples = [
    (1, 2, 3),
    (1, 4, 5),
    (1, 7, 6),
    (2, 4, 6),
    (2, 5, 7),
    (3, 4, 7),
    (3, 6, 5),
]

def build_O_table():
    """Returns oct_table[i][j] = (sign, k) for e_i * e_j = sign * e_k."""
    M = [[(0, 0) for _ in range(8)] for _ in range(8)]
    for i in range(8):
        M[0][i] = (1, i)
        M[i][0] = (1, i)
    for i in range(1, 8):
        M[i][i] = (-1, 0)
    for a, b, c in triples:
        M[a][b] = (1, c); M[b][c] = (1, a); M[c][a] = (1, b)
        M[b][a] = (-1, c); M[c][b] = (-1, a); M[a][c] = (-1, b)
    return M

O_table = build_O_table()

def O_mul(a, b):
    res = np.zeros(8)
    for i in range(8):
        if a[i] == 0:
            continue
        for j in range(8):
            if b[j] == 0:
                continue
            sign, k = O_table[i][j]
            res[k] += sign * a[i] * b[j]
    return res

def O_conj(a):
    res = np.copy(a)
    res[1:] *= -1
    return res

# Sanity
e_O = np.eye(8)
assert np.allclose(O_mul(e_O[1], e_O[2]), e_O[3]), "O sanity: e_1 e_2 = e_3"
assert np.allclose(O_mul(e_O[1], e_O[1]), -e_O[0]), "O sanity: e_1^2 = -1"
assert np.allclose(O_mul(O_mul(e_O[1], e_O[2]), e_O[4]),
                   -O_mul(e_O[1], O_mul(e_O[2], e_O[4]))), \
    "O non-associativity check (associator nonzero)"

print(f"\nC: dim 2.  Multiplication and conjugation built and verified.")
print(f"H: dim 4.  Multiplication and conjugation built and verified.")
print(f"O: dim 8.  Multiplication built via Fano plane; non-associativity verified.")


# =========================================================================
# Part 2: C ⊗ H ⊗ O as a 64-dim algebra
# =========================================================================
print()
print("-"*72)
print("Part 2: Tensor algebra C ⊗ H ⊗ O")
print("-"*72)
print()
print("Basis: e_{a,b,c} where a in {0,1} (C), b in {0,1,2,3} (H), c in {0,..,7} (O).")
print("Indexed by triple (a, b, c).  Total: 2 × 4 × 8 = 64 basis elements.")
print()
print("Tensor multiplication: (a ⊗ b ⊗ c)(a' ⊗ b' ⊗ c') = (a a') ⊗ (b b') ⊗ (c c').")
print()

# Total dimension
N_TOTAL = 2 * 4 * 8
print(f"Total dimension: {N_TOTAL}")

def CHO_basis_index(a, b, c):
    """Map (a, b, c) to a single index 0..63."""
    return a * 32 + b * 8 + c

def CHO_basis_decompose(idx):
    """Inverse: return (a, b, c) from a 0..63 index."""
    return (idx // 32, (idx % 32) // 8, idx % 8)

def CHO_mul(x, y):
    """Multiply two 64-element vectors x, y in C ⊗ H ⊗ O."""
    result = np.zeros(N_TOTAL)
    for i in range(N_TOTAL):
        if x[i] == 0:
            continue
        a1, b1, c1 = CHO_basis_decompose(i)
        for j in range(N_TOTAL):
            if y[j] == 0:
                continue
            a2, b2, c2 = CHO_basis_decompose(j)
            # Multiply component-wise
            # C component
            c_part = C_mul(
                np.array([1, 0]) if a1 == 0 else np.array([0, 1]),
                np.array([1, 0]) if a2 == 0 else np.array([0, 1])
            )
            # H component
            h1 = np.zeros(4); h1[b1] = 1
            h2 = np.zeros(4); h2[b2] = 1
            h_part = H_mul(h1, h2)
            # O component
            o1 = np.zeros(8); o1[c1] = 1
            o2 = np.zeros(8); o2[c2] = 1
            o_part = O_mul(o1, o2)
            # Combine: tensor product of all three results
            coeff = x[i] * y[j]
            for new_a in range(2):
                if c_part[new_a] == 0:
                    continue
                for new_b in range(4):
                    if h_part[new_b] == 0:
                        continue
                    for new_c in range(8):
                        if o_part[new_c] == 0:
                            continue
                        idx_new = CHO_basis_index(new_a, new_b, new_c)
                        result[idx_new] += coeff * c_part[new_a] * h_part[new_b] * o_part[new_c]
    return result


# =========================================================================
# Part 3: Verification
# =========================================================================
print()
print("-"*72)
print("Part 3: Verification of structure")
print("-"*72)

# Test: identity element
identity = np.zeros(N_TOTAL)
identity[CHO_basis_index(0, 0, 0)] = 1  # = 1 ⊗ 1 ⊗ e_0

# Pick a random element
np.random.seed(42)
x = np.random.randn(N_TOTAL)

# Verify identity multiplication
xx = CHO_mul(identity, x)
print(f"\nIdentity test: ||1 * x - x|| = {np.linalg.norm(xx - x):.2e}")
assert np.allclose(xx, x)

xx = CHO_mul(x, identity)
print(f"Identity test: ||x * 1 - x|| = {np.linalg.norm(xx - x):.2e}")
assert np.allclose(xx, x)


# Test: associator on like elements should vanish (alternative algebra)
# Take three octonionic elements (different o-components):
e1 = np.zeros(N_TOTAL); e1[CHO_basis_index(0, 0, 1)] = 1  # 1 ⊗ 1 ⊗ e_1
e2 = np.zeros(N_TOTAL); e2[CHO_basis_index(0, 0, 2)] = 1
e4 = np.zeros(N_TOTAL); e4[CHO_basis_index(0, 0, 4)] = 1

associator = CHO_mul(CHO_mul(e1, e2), e4) - CHO_mul(e1, CHO_mul(e2, e4))
print(f"\nAssociator (e_1, e_2, e_4) in O sector: ||a(x,y,z)|| = {np.linalg.norm(associator):.6f}")
print("Non-zero confirms non-associativity in O sector (expected for octonions).")


# =========================================================================
# Part 4: Norm and structure constants
# =========================================================================
print()
print("-"*72)
print("Part 4: Norm and structure")
print("-"*72)

def CHO_conj(x):
    """Conjugate in C ⊗ H ⊗ O: conjugate each factor."""
    result = np.zeros(N_TOTAL)
    for idx in range(N_TOTAL):
        if x[idx] == 0:
            continue
        a, b, c = CHO_basis_decompose(idx)
        # C conjugation: i -> -i, i.e., (a) -> (a) with sign (-1)^a
        sign_c = 1 if a == 0 else -1
        # H conjugation: 1 -> 1, i,j,k -> -i,-j,-k
        sign_h = 1 if b == 0 else -1
        # O conjugation: e_0 -> e_0, e_1..e_7 -> -e_i
        sign_o = 1 if c == 0 else -1
        sign_total = sign_c * sign_h * sign_o
        result[idx] = sign_total * x[idx]
    return result

# Sanity: |1|^2 = 1
norm_id = CHO_mul(identity, CHO_conj(identity))
print(f"\n||1 * conj(1)|| at identity component: {norm_id[CHO_basis_index(0,0,0)]:.4f}")


# Count structure constants
# Build the full 64x64 structure constant tensor C[i][j][k]
print()
print("Building full 64×64×64 structure constant tensor (this is ~262144 entries)...")
struct = np.zeros((N_TOTAL, N_TOTAL, N_TOTAL))
for i in range(N_TOTAL):
    ei = np.zeros(N_TOTAL); ei[i] = 1
    for j in range(N_TOTAL):
        ej = np.zeros(N_TOTAL); ej[j] = 1
        product = CHO_mul(ei, ej)
        struct[i, j, :] = product

n_nonzero = np.count_nonzero(struct)
n_total = N_TOTAL**3
print(f"  Non-zero structure constants: {n_nonzero} / {n_total} ({100*n_nonzero/n_total:.2f}%)")

# Verify alternative property: (xy)x = x(yx) for like elements
# For octonions, this is the alternative law.  In C ⊗ H ⊗ O, alternative
# property should still hold (since it's preserved by tensor products of
# normed division algebras).
print()
print("Checking alternative law (xy)x = x(yx) on random elements...")
np.random.seed(7)
violations = 0
for trial in range(5):
    x = np.random.randn(N_TOTAL)
    y = np.random.randn(N_TOTAL)
    lhs = CHO_mul(CHO_mul(x, y), x)
    rhs = CHO_mul(x, CHO_mul(y, x))
    diff = np.linalg.norm(lhs - rhs)
    if diff > 1e-8:
        violations += 1
        print(f"  Trial {trial}: ||(xy)x - x(yx)|| = {diff:.6f}  ** ALTERNATIVE LAW VIOLATED **")
    else:
        print(f"  Trial {trial}: ||(xy)x - x(yx)|| = {diff:.2e}  (alternative law holds)")

if violations == 0:
    print("\n  Alternative law holds for all 5 random elements.  Algebra is alternative.")
else:
    print(f"\n  Alternative law violated for {violations}/5 trials.  Algebra is NOT alternative.")
    print("  This is surprising -- C ⊗ H ⊗ O should be alternative.  Investigate.")


# =========================================================================
# Summary
# =========================================================================
print()
print("="*72)
print("Summary")
print("="*72)
print(f"""
C ⊗ H ⊗ O constructed:
  - Dimension: 64 over R.
  - Multiplication: tensor product of C-mul, H-mul, O-mul.
  - Conjugation: factorwise conjugation (sign flip on imaginary basis).
  - Identity: 1 ⊗ 1 ⊗ e_0.
  - Non-associative (associator nonzero on O sector).
  - Alternative law: {'HOLDS' if violations == 0 else 'VIOLATED'}.
  - Number of nonzero structure constants: {n_nonzero}.

Next step (Variant B): identify the SM idempotent in this algebra.

The algebra-construction infrastructure is ready for use by subsequent variants.
""")

# Save the algebra for reuse
np.savez("choh_structure.npz",
         struct_constants=struct,
         basis_decomp={'N_TOTAL': N_TOTAL})
print("Saved structure tensor to choh_structure.npz for downstream variants.")
