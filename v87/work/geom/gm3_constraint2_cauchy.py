#!/usr/bin/env python3
"""
gm3_constraint2_cauchy.py
GEOM Constraint 2: field monism alone does NOT escape Bell.

Construction:
  λ = full field configuration on a Cauchy slice Σ in the common past.
  Dynamics local (finite light cone / stencil) ⇒ A = A(a, λ|_cone_A),
  B = B(b, λ|_cone_B). These are of the form A(a,λ), B(b,λ) with (R)+(L).
  Under (MI) ρ(λ|a,b)=ρ(λ), CHSH ≤ 2.

Holism ("the pair is one field structure") does not break this: the whole
slice is still a valid λ. What must break is free specification of slice data
OR measurement independence OR locality of responses.

Also re-prove the elementary |s(λ)|=2 fact.

Tags: [D]
"""
from __future__ import annotations
import itertools
import math


def chsh_s_lambda():
    """For every assignment of ±1 to A,A',B,B': s = A(B+B')+A'(B-B') has |s|=2."""
    results = []
    for A, Ap, B, Bp in itertools.product([-1, 1], repeat=4):
        s = A * (B + Bp) + Ap * (B - Bp)
        results.append(s)
    assert all(abs(s) == 2 for s in results)
    return {
        "n_assignments": len(results),
        "unique_s": sorted(set(results)),
        "all_abs_eq_2": all(abs(s) == 2 for s in results),
    }


def cauchy_slice_model_is_local():
    """
    Abstract model of Constraint 2:
    - λ is a vector in R^N representing slice data.
    - Alice's outcome depends only on coordinates in set IA and setting a.
    - Bob's on IB and b.
    - If IA ∩ IB may be nonempty (shared past region of the slice) — holism
      of the common past is allowed — but A still independent of b, B of a.

    Under MI, S ≤ 2. Numeric: sample random local responses on discrete λ.
    """
    # λ ∈ {0,1,...,K-1}; random deterministic A(a,λ), B(b,λ)
    rng_seed_vals = []
    import random

    random.seed(42)
    K = 32
    # For each of 4 settings pairs, same ρ under MI
    max_S = -1e9
    for trial in range(200):
        # responses
        A_tab = {
            0: [random.choice([-1, 1]) for _ in range(K)],
            1: [random.choice([-1, 1]) for _ in range(K)],
        }
        B_tab = {
            0: [random.choice([-1, 1]) for _ in range(K)],
            1: [random.choice([-1, 1]) for _ in range(K)],
        }
        # uniform ρ(λ)
        Es = {}
        for a in (0, 1):
            for b in (0, 1):
                e = sum(A_tab[a][k] * B_tab[b][k] for k in range(K)) / K
                Es[(a, b)] = e
        S = Es[(0, 0)] + Es[(0, 1)] + Es[(1, 0)] - Es[(1, 1)]
        if abs(S) > max_S:
            max_S = abs(S)
    return {
        "max_abs_S_over_200_random_local_models": max_S,
        "classical_bound": 2.0,
        "respects_CHSH": max_S <= 2.0 + 1e-12,
    }


def monism_does_not_break_factorisability():
    """
    Factorisability of the *measure* P(A,B|a,b) = ∫ ρ(λ) 1_{A=A(a,λ)} 1_{B=B(b,λ)}
    holds whenever (R)+(L)+(MI), whether or not the physical objects are
    'one field'. Holism of the *ontology* is orthogonal to factorisability
    of the *hidden-variable measure*.

    Howard 1985 separability vs locality is the conceptual anchor
    [UNVERIFIED primary this session; distinction used as definitional].
    """
    return {
        "claim": (
            "O2 holism does not imply non-factorisability; "
            "only MD, nonlocal responses, or non-Cauchy λ can."
        ),
        "mechanisms_that_work": [
            "Door A: ρ(λ|a,b) ≠ ρ(λ)",
            "Door B: A or B depends on distant setting",
            "Door C / global constraint: λ not free on a spacelike slice",
            "Topological holonomy: slice data over-determined (conjectural)",
        ],
    }


if __name__ == "__main__":
    print("=== gm3 Constraint 2 ===")
    r1 = chsh_s_lambda()
    print("s(λ) fact:", r1)
    r2 = cauchy_slice_model_is_local()
    print("random local MI models:", r2)
    r3 = monism_does_not_break_factorisability()
    print(r3)
    assert r1["all_abs_eq_2"]
    assert r2["respects_CHSH"]
    print("ALL CHECKS PASSED — Constraint 2 formalised")
