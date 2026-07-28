#!/usr/bin/env python3
"""
gm1_tsirelson_geometry.py
GEOM seat — geometric derivation of Tsirelson bound when E(a,b) = -a·b.

Hypothesis (bound half): once correlations are inner products of unit vectors
in a real inner-product space, |S| ≤ 2√2 follows from parallelogram + Cauchy–Schwarz.
This is elementary; the real question is WHY E is an inner product (gm5).

Tags: [D] algebra; [M] numeric max check.
"""
from __future__ import annotations
import math
import numpy as np
import sympy as sp

# ---------------------------------------------------------------------------
# Part A — symbolic: ||a+a'|| + ||a-a'|| ≤ 2√2 for unit a,a'
# ---------------------------------------------------------------------------
def symbolic_vector_bound():
    # Represent a, a' as unit vectors in R^2 (planar settings suffice for CHSH)
    # Use symbolic norms of sum/diff
    t = sp.symbols("t", real=True)  # angle between a and a'
    # ||a+a'||^2 = 2 + 2 cos t = 4 cos^2(t/2)
    # ||a-a'||^2 = 2 - 2 cos t = 4 sin^2(t/2)
    n_plus = 2 * sp.Abs(sp.cos(t / 2))
    n_minus = 2 * sp.Abs(sp.sin(t / 2))
    sum_nm = sp.simplify(n_plus + n_minus)
    # Max of 2(|cos(t/2)| + |sin(t/2)|) over t
    # On [0, π]: cos(t/2)≥0, sin(t/2)≥0 ⇒ f = 2(cos+sin)(t/2)
    # df=0: -sin+cos=0 ⇒ t/2=π/4 ⇒ t=π/2; f=2√2
    f = 2 * (sp.cos(t / 2) + sp.sin(t / 2))
    df = sp.diff(f, t)
    crit = sp.solve(df, t)
    crit_in_range = [c for c in crit if c.is_real and 0 <= float(c) <= math.pi + 1e-12]
    vals = [(c, float(f.subs(t, c))) for c in crit_in_range]
    # endpoints
    ends = [(0, float(f.subs(t, 0))), (sp.pi, float(f.subs(t, sp.pi)))]
    max_val = max([v for _, v in vals + ends])
    # Exact identity: max 2(cos+sin)(θ) for θ∈[0,π/2] = 2√2 at θ=π/4
    exact = sp.simplify(2 * (sp.cos(sp.pi / 4) + sp.sin(sp.pi / 4)))
    assert exact == 2 * sp.sqrt(2)
    # Parallelogram: ||u||^2 + ||v||^2 = 2(||a||^2+||a'||^2)=4 for u=a+a',v=a-a'
    # CS: ||u|| + ||v|| ≤ √2 √(||u||^2+||v||^2) = √8 = 2√2
    return {
        "sum_nm_expr": str(sum_nm),
        "critical_points": vals,
        "endpoints": ends,
        "max_numeric": max_val,
        "max_exact": str(exact),
        "equals_2sqrt2": bool(exact == 2 * sp.sqrt(2)),
    }


# ---------------------------------------------------------------------------
# Part B — CHSH when E(a,b) = -a·b
# S = E(a,b)+E(a,b')+E(a',b)-E(a',b')
#   = -(a·b + a·b' + a'·b - a'·b')
#   = -[(a+a')·b + (a-a')·b']   wait: a'·b - a'·b' = a'·(b-b')
# S = -a·b -a·b' -a'·b + a'·b' = -a·(b+b') - a'·(b - b')
# |S| ≤ ||a|| ||b+b'|| + ||a'|| ||b-b'|| = ||b+b'|| + ||b-b'|| ≤ 2√2
# ---------------------------------------------------------------------------
def chsh_bound_symbolic():
    # Unit b,b' in R^2: b=(1,0), b'=(cos φ, sin φ)
    phi = sp.symbols("phi", real=True)
    bp = sp.Matrix([1, 0])
    bm = sp.Matrix([sp.cos(phi), sp.sin(phi)])
    s_vec = bp + bm
    d_vec = bp - bm
    n_s = sp.simplify(sp.sqrt(s_vec.dot(s_vec)))
    n_d = sp.simplify(sp.sqrt(d_vec.dot(d_vec)))
    bound = sp.simplify(n_s + n_d)
    # At orthogonal settings φ=π/2: both norms = √2, sum = 2√2
    at_orth = bound.subs(phi, sp.pi / 2)
    return {
        "n_s": str(n_s),
        "n_d": str(n_d),
        "bound": str(bound),
        "at_phi_pi2": str(sp.simplify(at_orth)),
        "at_phi_pi2_num": float(at_orth),
        "is_2sqrt2": bool(sp.simplify(at_orth - 2 * sp.sqrt(2)) == 0),
    }


def chsh_at_canonical_angles():
    """Standard Tsirelson settings in plane: a,a',b,b' unit vectors."""
    # b = e_x rotated by π/4 etc. Use 2D angles from x-axis.
    # Common choice: a at 0, a' at π/2, b at π/4, b' at -π/4
    def u(th):
        return np.array([math.cos(th), math.sin(th)])

    a, ap = u(0.0), u(math.pi / 2)
    b, bp = u(math.pi / 4), u(-math.pi / 4)

    def E(x, y):
        return -float(np.dot(x, y))

    S = E(a, b) + E(a, bp) + E(ap, b) - E(ap, bp)
    return {
        "E_ab": E(a, b),
        "E_abp": E(a, bp),
        "E_apb": E(ap, b),
        "E_apbp": E(ap, bp),
        "S": S,
        "abs_S": abs(S),
        "S_exact": 2 * math.sqrt(2),
        "residual": abs(abs(S) - 2 * math.sqrt(2)),
    }


def numeric_global_max(n_starts=80, seed=20260726):
    """Maximise |S| over 4 angles on the circle for E=-cos(θ_i-θ_j)."""
    rng = np.random.default_rng(seed)

    def S_of(angs):
        aa, aap, bb, bbp = angs
        def ec(x, y):
            return -math.cos(x - y)
        return ec(aa, bb) + ec(aa, bbp) + ec(aap, bb) - ec(aap, bbp)

    best = -1e9
    best_ang = None
    for _ in range(n_starts):
        x0 = rng.uniform(0, 2 * math.pi, size=4)
        # simple coordinate ascent
        x = x0.copy()
        for _it in range(80):
            for k in range(4):
                grid = np.linspace(0, 2 * math.pi, 72, endpoint=False)
                vals = []
                for g in grid:
                    xx = x.copy()
                    xx[k] = g
                    vals.append(S_of(xx))
                x[k] = grid[int(np.argmax(vals))]
        s = S_of(x)
        if s > best:
            best = s
            best_ang = x.copy()
    return {
        "S_max_numeric": best,
        "S_target": 2 * math.sqrt(2),
        "abs_err": abs(best - 2 * math.sqrt(2)),
        "angles": best_ang.tolist() if best_ang is not None else None,
    }


def parallelogram_cs_proof_numeric():
    """Sample random unit a,a'; check ||a+a'||+||a-a'|| ≤ 2√2."""
    rng = np.random.default_rng(7)
    max_sum = 0.0
    n = 20000
    for _ in range(n):
        a = rng.normal(size=3)
        a /= np.linalg.norm(a)
        ap = rng.normal(size=3)
        ap /= np.linalg.norm(ap)
        s = np.linalg.norm(a + ap) + np.linalg.norm(a - ap)
        if s > max_sum:
            max_sum = s
    return {
        "max_over_samples": max_sum,
        "bound": 2 * math.sqrt(2),
        "all_below": max_sum <= 2 * math.sqrt(2) + 1e-12,
        "gap": 2 * math.sqrt(2) - max_sum,
    }


if __name__ == "__main__":
    print("=== gm1: geometric Tsirelson ===")
    r1 = symbolic_vector_bound()
    print("A symbolic vector bound:", r1)
    r2 = chsh_bound_symbolic()
    print("B CHSH vector bound:", r2)
    r3 = chsh_at_canonical_angles()
    print("C canonical angles:", r3)
    r4 = parallelogram_cs_proof_numeric()
    print("D parallelogram samples:", r4)
    r5 = numeric_global_max()
    print("E numeric S max:", r5)
    assert r1["equals_2sqrt2"]
    assert r2["is_2sqrt2"]
    assert r3["residual"] < 1e-12
    assert r4["all_below"]
    assert r5["abs_err"] < 1e-6
    print("ALL CHECKS PASSED")
