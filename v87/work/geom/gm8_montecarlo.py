#!/usr/bin/env python3
"""
gm8_montecarlo.py
Monte Carlo of the geometric singlet model (FABLE p=1 / E=−a·b via MD).

Event source: explicit PRNG (numpy Generator PCG64).
Sampling: λ on S² with density ∝ |λ·b|  (p=1), responses A=sgn(a·λ), B=−sgn(b·λ).

Confirms S → 2√2 at canonical angles; S→2 at p=0; overshoot for large p.

Tags: [M]
"""
from __future__ import annotations
import math
import numpy as np


def sample_lambda_p(rng, b, p, n):
    """Sample n unit vectors with density ∝ |λ·b|^p on S².

    Method: in b-frame, z-component: pdf of u=λ·b on [-1,1] is
      ρ(u) ∝ |u|^p  (area element ⇒ factor flat in u after solid angle),
    actually surface measure: λ·b = cosθ, dA ∝ sinθ dθ dφ = −d(cosθ) dφ,
    so u ∈[-1,1] uniform under p=0; under density |u|^p, pdf ∝ |u|^p on [-1,1].
    Normalisation: ∫_{-1}^1 |u|^p du = 2/(p+1), pdf = ((p+1)/2) |u|^p.
    Sample |u| = U^{1/(p+1)}, sign fair; φ uniform; rotate so ê_z → b.
    """
    u_abs = rng.random(n) ** (1.0 / (p + 1.0))
    signs = rng.choice([-1.0, 1.0], size=n)
    u = u_abs * signs
    phi = rng.uniform(0, 2 * math.pi, size=n)
    # in b-frame: λ = (sqrt(1-u^2) cosφ, sqrt(1-u^2) sinφ, u)
    r_xy = np.sqrt(np.maximum(0.0, 1.0 - u * u))
    lx = r_xy * np.cos(phi)
    ly = r_xy * np.sin(phi)
    lz = u
    # rotate ê_z to b
    b = np.asarray(b, dtype=float)
    b = b / np.linalg.norm(b)
    ez = np.array([0.0, 0.0, 1.0])
    v = np.cross(ez, b)
    c = float(np.dot(ez, b))
    if np.linalg.norm(v) < 1e-12:
        # b parallel or anti to ez
        R = np.eye(3) if c > 0 else np.diag([1.0, -1.0, -1.0])
    else:
        vx = np.array(
            [
                [0, -v[2], v[1]],
                [v[2], 0, -v[0]],
                [-v[1], v[0], 0],
            ]
        )
        R = np.eye(3) + vx + vx @ vx * (1.0 / (1.0 + c))
    L = np.stack([lx, ly, lz], axis=1)  # (n,3)
    return L @ R.T


def correlator_mc(rng, a, b, p, n):
    a = np.asarray(a, dtype=float)
    a = a / np.linalg.norm(a)
    b = np.asarray(b, dtype=float)
    b = b / np.linalg.norm(b)
    L = sample_lambda_p(rng, b, p, n)
    A = np.sign(L @ a)
    B = -np.sign(L @ b)
    # sign(0)=0 rare
    prod = A * B
    # drop zeros
    prod = prod[prod != 0]
    E = float(np.mean(prod))
    se = float(np.std(prod, ddof=1) / math.sqrt(len(prod)))
    return E, se


def unit(th, ph=0.0):
    return np.array(
        [math.sin(th) * math.cos(ph), math.sin(th) * math.sin(ph), math.cos(th)]
    )


def run_suite():
    rng = np.random.default_rng(20260726)
    n = 2_000_000
    # Settings as vectors: a=ẑ, a'=x̂, b=(ẑ+x̂)/√2, b'=(ẑ−x̂)/√2
    a = np.array([0.0, 0.0, 1.0])
    ap = np.array([1.0, 0.0, 0.0])
    b = np.array([1.0, 0.0, 1.0]) / math.sqrt(2)
    bp = np.array([-1.0, 0.0, 1.0]) / math.sqrt(2)

    results = {}
    for p, label in [(1.0, "p1_quantum"), (0.0, "p0_classical"), (8.0, "p8_overshoot")]:
        Eab, se_ab = correlator_mc(rng, a, b, p, n)
        Eabp, se_abp = correlator_mc(rng, a, bp, p, n)
        Eapb, se_apb = correlator_mc(rng, ap, b, p, n)
        Eapbp, se_apbp = correlator_mc(rng, ap, bp, p, n)
        S = Eab + Eabp + Eapb - Eapbp
        se_S = math.sqrt(se_ab**2 + se_abp**2 + se_apb**2 + se_apbp**2)
        results[label] = {
            "p": p,
            "E": (Eab, Eabp, Eapb, Eapbp),
            "S": S,
            "se_S": se_S,
            "target_p1": 2 * math.sqrt(2),
        }
    return results


if __name__ == "__main__":
    print("=== gm8 Monte Carlo ===")
    res = run_suite()
    for k, v in res.items():
        print(k, v)
    # Checks
    p1 = res["p1_quantum"]
    assert abs(abs(p1["S"]) - 2 * math.sqrt(2)) < 5 * p1["se_S"]  # within ~5σ
    p0 = res["p0_classical"]
    assert abs(p0["S"]) <= 2.05  # classical ≤2 with MC noise
    p8 = res["p8_overshoot"]
    assert abs(p8["S"]) > 2 * math.sqrt(2) + 0.2  # clear overshoot
    print("ALL CHECKS PASSED")
    print(f"p=1: S={p1['S']:.6f} ± {p1['se_S']:.6f}  |S| vs 2√2={2*math.sqrt(2):.6f}")
    print(f"p=0: S={p0['S']:.6f} ± {p0['se_S']:.6f}")
    print(f"p=8: S={p8['S']:.6f} ± {p8['se_S']:.6f}")
