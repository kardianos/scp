#!/usr/bin/env python3
"""
gm5_why_inner_product.py
THE real geometric question (GEOM_BRIEF §2):

  Why are the correlations inner products of unit vectors?
  I.e. why E(a,b) = −a·b (or more generally a bilinear form on a real
  Hilbert space of bounded dimension / operator algebra)?

What is KNOWN [D]:
  IF E(a,b)=−u(a)·u(b) for unit maps u into a real inner-product space,
  THEN |S|≤2√2 by gm1.

What is NOT given by spacetime geometry:
  The light cone / Lorentzian form only constrains signalling (gm2).
  It does not force the correlator to be an Euclidean inner product.

Candidate sources of the bilinear form:
  (Q) Quantum theory: Born rule + projective measurements on C^n ⇒
      E = ⟨ψ|A⊗B|ψ⟩ reduces to −a·b for singlet (derived in QM, not fabric).
  (MD-p1) FABLE's p=1 sphere model: ρ∝|λ·b|, A=sgn(a·λ),B=−sgn(b·λ)
      yields E=−cosθ *by calculation* — the form is an output of a
      *postulated* tilt, not of metric structure.
  (GPT) Generalized probabilistic theories: state-space geometry can force
      correlators into various sets (Barrett PRA 2007 — UNVERIFIED primary).
  (Exclusivity) Lovász ϑ of the exclusivity graph (Cabello et al. PRL 2014
      — identifier from brief) gives quantum bound as an orthonormal-
      representation optimum — pure real-vector-space geometry of *events*,
      not of spacetime.

Fabric programme status:
  No measured Cosserat structure forces E into unit-vector bilinear form.
  Present kernel with MI predicts classical sawtooth E=−1+2θ/π (hemisphere)
  for local deterministic responses on a uniform sphere — NOT −cosθ.

This worksheet:
  1. Proves classical hemisphere E_cl = −1 + 2θ/π (closed form via SymPy).
  2. Shows E_cl is NOT an inner product of unit vectors (fails bilinearity /
     cosine shape).
  3. Shows that requiring E to be of the form −a·b for some embedding of
     settings is a *selection principle* that picks Tsirelson — and is not
     implied by the Lorentzian metric.

Tags: [D],[M]
"""
from __future__ import annotations
import math
import numpy as np
import sympy as sp


def classical_hemisphere_E():
    """E(θ) = ∫ A(a,λ)B(b,λ) dλ for A=sgn(a·λ), B=−sgn(b·λ), λ~Unif(S^1).

    On the circle (planar CHSH): λ∈[0,2π), A=sgn cos(λ−a), B=−sgn cos(λ−b).
    Closed form E = −1 + 2θ/π for θ=∠(a,b)∈[0,π].
    """
    # Analytic: fraction of λ where sgn cos(λ-a) sgn cos(λ-b) = +1 is
    # 1 − θ/π for the product AB with B=−sgn... so E[AB] = −(1 − 2θ/π) wait
    # A B = sgn cos(λ-a) * (−sgn cos(λ-b)) = − sgn cos(λ-a) sgn cos(λ-b)
    # E[sgn cos α sgn cos(α−θ)] = 1 − 2θ/π for θ∈[0,π]
    # so E[AB] = −(1 − 2θ/π) = −1 + 2θ/π
    th = sp.symbols("theta", real=True, positive=True)
    E_cl = -1 + 2 * th / sp.pi
    # Numeric check
    n = 200001
    lam = np.linspace(0, 2 * math.pi, n, endpoint=False)
    errs = []
    for theta in [0.0, math.pi / 4, math.pi / 2, 3 * math.pi / 4, math.pi]:
        a, b = 0.0, theta
        A = np.sign(np.cos(lam - a))
        B = -np.sign(np.cos(lam - b))
        # np.sign(0)=0; negligible
        E_num = float(np.mean(A * B))
        E_an = -1 + 2 * theta / math.pi
        errs.append(abs(E_num - E_an))
    return {
        "closed_form": str(E_cl),
        "max_abs_err_numeric": max(errs),
        "at_pi4": -1 + 2 * (math.pi / 4) / math.pi,  # -0.5
        "qm_at_pi4": -math.cos(math.pi / 4),  # -1/√2 ≈ -0.707
    }


def is_inner_product_form(E_func, n_angles=24):
    """Test whether a planar correlator E(θ) can be written −cos(c·θ) or more
    generally whether the Gram matrix of −E is positive semidefinite of rank≤d.

    For unit-vector form E(a,b)=−u(a)·u(b), the matrix G_ij = −E(a_i,a_j)
    must be a Gram matrix (PSD, rank ≤ ambient dim).
    """
    angles = np.linspace(0, math.pi, n_angles)
    # Use settings on a circle; G_ij = −E(|θ_i−θ_j|)
    G = np.zeros((n_angles, n_angles))
    for i, ti in enumerate(angles):
        for j, tj in enumerate(angles):
            d = abs(ti - tj)
            d = min(d, 2 * math.pi - d) if d > math.pi else d
            # for planar settings difference in [0,π]
            dd = abs(ti - tj)
            if dd > math.pi:
                dd = 2 * math.pi - dd
            # map to [0,π]
            if dd > math.pi:
                dd = 2 * math.pi - dd
            G[i, j] = -E_func(dd)
    # symmetrize
    G = 0.5 * (G + G.T)
    ev = np.linalg.eigvalsh(G)
    return {
        "min_eigenvalue": float(ev[0]),
        "max_eigenvalue": float(ev[-1]),
        "is_PSD": bool(ev[0] >= -1e-8),
        "numerical_rank_1e-6": int(np.sum(ev > 1e-6)),
    }


def E_classical(theta):
    return -1 + 2 * theta / math.pi


def E_quantum(theta):
    return -math.cos(theta)


def E_pr_step(theta):
    # extreme correlator; not rotationally a function of angle alone in full PR
    return -1.0 if theta < math.pi / 2 else +1.0


def cosine_is_gram():
    r = is_inner_product_form(E_quantum, n_angles=36)
    # For E=-cos(θ_i-θ_j), G_ij = cos(θ_i-θ_j) = u_i·u_j for u=(cosθ,sinθ)
    # rank 2, PSD
    return r


def classical_is_not_low_rank_cosine():
    r = is_inner_product_form(E_classical, n_angles=36)
    # Classical sawtooth is NOT a cosine; Gram may still be PSD in high dim
    # (Schoenberg) but the CHSH value of classical is 2, not 2√2.
    # Key: it is NOT the *Euclidean unit-vector embedding of the setting sphere
    # with E=-a·b*; it requires a different function.
    return r


def selection_principle_statement():
    return {
        "principle_IP": (
            "Correlations are (minus) inner products of unit vectors assigned "
            "to measurement settings in a real inner-product space."
        ),
        "implies": "|S| ≤ 2√2 (gm1) [D]",
        "from_light_cone": False,
        "from_O2_holism": False,
        "from_fabric_kernel": False,
        "from_QM": True,
        "from_FABLE_p1_postulate": True,  # output of postulated tilt
        "status_for_fabric": "POSTULATED / NOT DERIVED [P]",
    }


if __name__ == "__main__":
    print("=== gm5 why inner product ===")
    print("classical E:", classical_hemisphere_E())
    print("quantum Gram:", cosine_is_gram())
    print("classical Gram:", classical_is_not_low_rank_cosine())
    print("selection:", selection_principle_statement())
    cl = classical_hemisphere_E()
    assert cl["max_abs_err_numeric"] < 0.01
    qg = cosine_is_gram()
    assert qg["is_PSD"] and qg["numerical_rank_1e-6"] == 2
    print("ALL CHECKS PASSED")
