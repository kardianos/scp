#!/usr/bin/env python3
"""
gm4_operator_S2.py
Operator route to Tsirelson: for dichotomic observables A,A',B,B' with
A^2=A'^2=B^2=B'^2=I and [A,B]=0 etc. (tensor product Alice ⊗ Bob),

  S = A⊗B + A⊗B' + A'⊗B − A'⊗B'
  S^2 = 4 I − [A,A']⊗[B,B']   (sign fixed by the CHSH minus on A'B')

Hence ||S|| ≤ √(4 + ||[A,A']|| ||[B,B']||) ≤ √(4+4) = 2√2
since ||[X,Y]|| ≤ 2 for dichotomic observables.

This is the *same bound* as the vector parallelogram proof: both exploit a
2-norm / C*-norm quadratic estimate. Neither is the light cone.

Symbolic check with Pauli matrices (qubit singlet / CHSH optimal).

Tags: [D],[M]
"""
from __future__ import annotations
import math
import numpy as np
import sympy as sp


def pauli():
    I = sp.Matrix([[1, 0], [0, 1]])
    X = sp.Matrix([[0, 1], [1, 0]])
    Y = sp.Matrix([[0, -sp.I], [sp.I, 0]])
    Z = sp.Matrix([[1, 0], [0, -1]])
    return I, X, Y, Z


def kron(A, B):
    return sp.Matrix(sp.kronecker_product(A, B))


def operator_identity():
    I, X, Y, Z = pauli()
    # Standard CHSH operators: A=Z, A'=X, B=(Z+X)/√2, B'=(Z-X)/√2
    A = Z
    Ap = X
    B = (Z + X) / sp.sqrt(2)
    Bp = (Z - X) / sp.sqrt(2)
    # Check A^2 = I etc.
    assert sp.simplify(A * A - I) == sp.zeros(2)
    assert sp.simplify(Ap * Ap - I) == sp.zeros(2)
    # B^2: (Z+X)^2 / 2 = (2I + {Z,X})/2 = I since {Z,X}=0
    assert sp.simplify(sp.expand(B * B) - I) == sp.zeros(2)
    assert sp.simplify(sp.expand(Bp * Bp) - I) == sp.zeros(2)

    S = kron(A, B) + kron(A, Bp) + kron(Ap, B) - kron(Ap, Bp)
    S2 = sp.simplify(S * S)
    commA = sp.simplify(A * Ap - Ap * A)
    commB = sp.simplify(B * Bp - Bp * B)
    # Identity for this CHSH sign pattern: S² = 4I − [A,A']⊗[B,B']
    rhs = sp.simplify(4 * kron(I, I) - kron(commA, commB))
    residual = sp.simplify(S2 - rhs)
    # Spectral radius / eigenvalues of S
    ev = S.eigenvals()
    # Convert to floats for max |λ|
    norms = []
    for val, mult in ev.items():
        norms.append(abs(complex(val.evalf())))
    op_norm = max(norms)
    return {
        "S2_minus_rhs_is_zero": residual == sp.zeros(4),
        "op_norm_S": op_norm,
        "tsirelson": 2 * math.sqrt(2),
        "matches": abs(op_norm - 2 * math.sqrt(2)) < 1e-12,
        "commA": str(commA),
        "commB": str(commB),
    }


def vector_vs_operator_same_number():
    """Both routes give 2√2 — same algebraic object class (2-norm bound)."""
    return {
        "vector_route": "||b+b'|| + ||b-b'|| ≤ √2 √(||·||^2+||·||^2) = 2√2",
        "operator_route": "S²=4I−[A,A']⊗[B,B'] ⇒ ||S|| ≤ √(4+||[A,A']||||[B,B']||) ≤ 2√2",
        "light_cone_route": "no-signalling ⇒ |S|≤4 only (gm2)",
        "distinction": (
            "Spacetime light cone (Lorentzian quadratic form) enforces NS; "
            "Tsirelson uses a *positive-definite* 2-norm on outcome/Bloch space. "
            "Same *class* (quadratic form), different *instance*."
        ),
    }


if __name__ == "__main__":
    print("=== gm4 operator S^2 ===")
    r = operator_identity()
    print(r)
    print(vector_vs_operator_same_number())
    assert r["S2_minus_rhs_is_zero"]
    assert r["matches"]
    print("ALL CHECKS PASSED")
