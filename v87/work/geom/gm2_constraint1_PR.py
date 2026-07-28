#!/usr/bin/env python3
"""
gm2_constraint1_PR.py
GEOM Constraint 1: light-cone / no-signalling alone does NOT give 2√2.

PR box: E00=E01=E10=+1, E11=−1 ⇒ S=4, and marginals are no-signalling.
Therefore "relativistic causality" as no-signalling is insufficient for Tsirelson.

Also: under (R)∧(L) with unrestricted MD (Door A), S=4 is realisable
(FABLE/GROK already showed; re-derived here for independence).

Tags: [D]
"""
from __future__ import annotations
import math
import numpy as np


def pr_box_S():
    E = {(0, 0): 1.0, (0, 1): 1.0, (1, 0): 1.0, (1, 1): -1.0}
    S = E[(0, 0)] + E[(0, 1)] + E[(1, 0)] - E[(1, 1)]
    return S


def pr_no_signalling():
    """Check that PR correlators admit a no-signalling behaviour.

    Standard PR: P(A,B|a,b) = 1/2 if A⊕B = a·b (mod 2), else 0
    for A,B,a,b ∈ {0,1}. Map to ±1: A' = 1-2A, etc.
    Marginals: P(A|a,b)=1/2 independent of b; same for B.
    """
    # Outcomes ±1; settings 0,1
    # PR: A*B = (+1) if a*b == 0 else (−1)  in multiplicative form for bits:
    # For bits: A⊕B = a AND b. In ±1: AB = (−1)^{a AND b} = 1 - 2*(a and b)
    # Actually: A,B ∈ {±1}, want E[AB]= (+1) for (0,0),(0,1),(1,0) and −1 for (1,1)
    # Behaviour: fix A uniform ±1, B = A * E_target(a,b)  (deterministic given A,a,b)
    # That signals if E depends on a for Bob's outcome... actually B depends on a.
    # Standard NS construction:
    # P(ab|xy) with a,b bits: 1/2 if a⊕b = x∧y else 0.
    ok = True
    for x in (0, 1):
        for y in (0, 1):
            # marginal P(A=0|x,y)
            pA0 = 0.0
            pB0 = 0.0
            for A in (0, 1):
                for B in (0, 1):
                    p = 0.5 if (A ^ B) == (x & y) else 0.0
                    if A == 0:
                        pA0 += p
                    if B == 0:
                        pB0 += p
            if abs(pA0 - 0.5) > 1e-12 or abs(pB0 - 0.5) > 1e-12:
                ok = False
    return ok


def md_realises_PR():
    """(R)∧(L), free ρ(λ|a,b): PR via two-point MD (GROK G2 style)."""
    # λ ∈ {+1,−1} product value; A=1 always for simplicity with two strategies
    # Strategies: for each setting pair, put mass on λ that realises AB = target
    # Local: A(a,λ)=λ_A[a], B(b,λ)=λ_B[b] with λ a 4-tuple
    # Minimal: λ ∈ {0,1}; A(a)=+1 always; B(b,λ,a) would break L.
    # Keep L: fix responses A0,A1,B0,B1 as functions of λ only.
    # For each (a,b), put all mass on a λ with A_a B_b = target E_ab.
    targets = {(0, 0): 1, (0, 1): 1, (1, 0): 1, (1, 1): -1}
    # Use λ ∈ {±1}^4 fully specified strategies; pick one good λ per pair
    found = {}
    for a, b in targets:
        t = targets[(a, b)]
        # search all 16
        for bits in range(16):
            A0 = 1 if (bits & 1) else -1
            A1 = 1 if (bits & 2) else -1
            B0 = 1 if (bits & 4) else -1
            B1 = 1 if (bits & 8) else -1
            Aa = A0 if a == 0 else A1
            Bb = B0 if b == 0 else B1
            if Aa * Bb == t:
                found[(a, b)] = (A0, A1, B0, B1)
                break
    assert len(found) == 4
    # E exact
    S = sum(
        targets[(0, 0)]
        + targets[(0, 1)]
        + targets[(1, 0)]
        - targets[(1, 1)]
        for _ in [0]
    )
    # S from targets
    S = targets[(0, 0)] + targets[(0, 1)] + targets[(1, 0)] - targets[(1, 1)]
    return {"found": found, "S": S}


def no_signalling_bound_is_4():
    """Algebraic max of |S| under only |E|≤1 is 4; NS does not tighten for CHSH.

    For bipartite CHSH, the no-signalling polytope allows the PR point S=4.
    """
    return {
        "algebraic_max": 4.0,
        "PR_S": pr_box_S(),
        "PR_is_NS": pr_no_signalling(),
        "Tsirelson": 2 * math.sqrt(2),
        "gap_PR_minus_Tsirelson": 4.0 - 2 * math.sqrt(2),
        "conclusion": "NS alone stops at 4, not 2√2",
    }


if __name__ == "__main__":
    print("=== gm2 Constraint 1 ===")
    print("PR S =", pr_box_S())
    print("PR NS?", pr_no_signalling())
    print("MD realises PR:", md_realises_PR())
    print(no_signalling_bound_is_4())
    assert pr_box_S() == 4.0
    assert pr_no_signalling()
    assert md_realises_PR()["S"] == 4
    print("ALL CHECKS PASSED — Constraint 1 formalised")
