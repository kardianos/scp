#!/usr/bin/env python3
"""
gm7_join_p1_geometry.py
Join FABLE's Door-A p-family with GEOM bound half.

FABLE F2: ρ_p(λ|b) ∝ |λ·b|^p, A=sgn(a·λ), B=−sgn(b·λ)
  p=0 → E=−1+2θ/π, |S|max=2
  p=1 → E=−cosθ = −a·b, |S|max=2√2   ← lands ON the geometric form
  p→∞ → E→−sgn(cosθ), |S|max→4

GEOM: once E=−a·b, |S|≤2√2 is theorem (gm1). So:
  - Reachability above 2: from MD tilt (Door A) — FABLE/GROK
  - Bound at 2√2: automatic IF the tilt lands on the inner-product form
  - But the p-family does not STOP at p=1: geometry of spacetime does not
    select p=1; the inner-product form is one point on a continuum.

This worksheet re-derives E_1(θ)=−cosθ (SymPy) and S_max(p) at key points,
confirming the join and the gap.

Tags: [D],[M]
"""
from __future__ import annotations
import math
import numpy as np
from scipy.special import beta as beta_fn


def E_p_closed(theta, p):
    """FABLE closed form:
    E_p(θ) = −1 + ((p+1)/π) * B((p+2)/2, 1/2) * ∫_0^θ sin^p(φ) dφ
    """
    from scipy.integrate import quad

    B = beta_fn((p + 2) / 2.0, 0.5)
    integ, _ = quad(lambda phi: math.sin(phi) ** p, 0.0, theta, epsabs=1e-14)
    return -1.0 + ((p + 1) / math.pi) * B * integ


def E_p1_is_minus_cos():
    """Symbolic / high-precision: p=1 ⇒ E=−cosθ."""
    # For p=1: B(3/2,1/2) = Γ(3/2)Γ(1/2)/Γ(2) = (√π/2)*√π / 1 = π/2
    # E = −1 + (2/π)*(π/2)*∫_0^θ sin φ dφ = −1 + (1 − cosθ) = −cosθ
    B = beta_fn(1.5, 0.5)
    assert abs(B - math.pi / 2) < 1e-12
    errs = []
    for th in np.linspace(0, math.pi, 25):
        e = E_p_closed(float(th), 1.0)
        errs.append(abs(e + math.cos(th)))
    return {"B_3_2_1_2": B, "max_err_vs_minus_cos": max(errs)}


def S_at_tsirelson_angles(p):
    """S = E(ab)+E(ab')+E(a'b)−E(a'b') at canonical angles.
    For rotationally symmetric E=E(θ):
    θ(ab)=π/4, θ(ab')=π/4, θ(a'b)=π/4, θ(a'b')=3π/4
    ⇒ S = 3 E(π/4) − E(3π/4).  (Typically negative; report value and |S|.)
    """
    e1 = E_p_closed(math.pi / 4, p)
    e3 = E_p_closed(3 * math.pi / 4, p)
    return 3 * e1 - e3


def scan_p():
    rows = []
    for p in [0, 0.5, 1.0, 1.5, 2.0, 3.0, 10.0]:
        S = S_at_tsirelson_angles(p)
        # Note: for p≠1 these angles may not maximise |S|; still shows overshoot
        rows.append((p, S, 2 * math.sqrt(2)))
    return rows


def geometric_join_statement():
    return {
        "at_p1": "E=−a·b exactly ⇒ |S|≤2√2 by geometry (gm1), equality attained",
        "reachability": "p>0 MD (Door A) — postulated, cost I≥0.046 bits (FABLE/GROK)",
        "bound_source_at_p1": "inner-product 2-norm, NOT light cone",
        "does_geometry_select_p1": False,
        "proof_geometry_does_not_select": "S(p) continuous and increasing through p=1",
        "hypothesis_status": (
            "BOUND half works CONDITIONAL on E being unit-vector bilinear; "
            "REACH half is Door A MD; GEOMETRY does not glue them — the "
            "glue (why p=1 / why E=−a·b) remains postulated."
        ),
    }


if __name__ == "__main__":
    print("=== gm7 join p=1 ↔ geometry ===")
    r = E_p1_is_minus_cos()
    print("p=1 identity:", r)
    assert r["max_err_vs_minus_cos"] < 1e-10
    print("S at Tsirelson angles vs p:")
    for p, S, tgt in scan_p():
        print(f"  p={p:4.1f}  S={S:.6f}  |S|={abs(S):.6f}  2√2={tgt:.6f}  |S|-2√2={abs(S)-tgt:+.6f}")
    S1 = S_at_tsirelson_angles(1.0)
    assert abs(abs(S1) - 2 * math.sqrt(2)) < 1e-8
    S0 = S_at_tsirelson_angles(0.0)
    assert abs(abs(S0) - 2.0) < 1e-6  # classical |S|=2 at these angles
    # overshoot for large p
    assert abs(S_at_tsirelson_angles(10.0)) > 2 * math.sqrt(2) + 0.5
    print("join:", geometric_join_statement())
    print("ALL CHECKS PASSED")
