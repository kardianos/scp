#!/usr/bin/env python3
"""
01_toy_frequency_projection.py — Minimal formalization experiment for the
"octo-space frequency/phase projection" hypothesis.

Hypothesis (user-proposed, after v59-v61 dead-end assessment):
  The octonionic algebra (octo-space / Cl(7)_even ≅ ℂ⊗ℍ⊗𝕆) is the primary
  information substrate. What we perceive as physical quantities (Brannen
  masses, Koide Q, Brannen phase φ, structural integers, effective forces)
  are outputs of a frequency/spectral transform + phase projection applied
  to information in that algebra. The "before it comes real" layer is the
  algebraic frequency content; 3+1 observables are the projected result.

Goal of this toy:
  Test whether a *non-tautological* frequency model on the known order-3
  structures (sedenion S3 / triality Z3 acting on generation space) can
  produce the Brannen phase φ ≈ 2/9 and the equilibrium Koide Q = 2/3 as
  *outputs* of the projection, using only fixed algebra invariants
  (dim G2=14, dim Spin(7)=21, selection dims 28/35/63, etc.) as weights.
  No hand-insertion of 2/9 or the Brannen ratios into the projection rule.

Model (minimal, to avoid circularity):
  - Raw "octo-info" for the generation sector: a real 3-vector v on Z3
    (the sedenion/triality 3-cycle). Input v is *not* pre-loaded with
    φ=2/9 or Brannen eigenvalue ratios.
  - Frequency transform: DFT on Z3 (DC + two complex conjugate modes
    at frequencies ±1, with ω = exp(2π i /3)).
  - Projection rule (the thing being tested): a phase-extraction operator
    on the frequency components, weighted by algebra-derived factors
    (G2 content ratio 14/21 = Q, or 28/3 universal deviation constant,
    or selection-grade dims). The rule must not contain 2/9 or the
    target mass ratios as free parameters.
  - Output observables: extracted phase φ, effective "masses" m_k from
    the projected frequency components, computed Q from the m_k.
  - Success test: does there exist a natural (non-tuned) choice of the
    weighting such that the output φ satisfies 3φ ≈ 2/3 (the radian
    insert) *and* the effective Q is within 1% of 2/3, while the input
    v does not already encode these values?

This is deliberately a toy on the 3-generation Z3 factor. If it collapses
to tautology or requires the answer to be inserted, that is a clean
negative for this particular formalization. If it produces the phase
from frequency content + known dims, that is a positive signal worth
extending to the full algebra (L-grade complex structures, End(L), etc.).

References:
  - v59 G7 radian-insert gap (LeptonPhaseMagnitude.lean TODO for spectral
    character that produces cos(2/3) with argument = Q).
  - v60 gravity_recast/09 (Plebański parent is independent posit, not
    derived from OBE — motivation for changing ontology).
  - v60/gaps/rank_tension/01 (two-object resolution; 784 vs rank-3 is
    structural, not a single-Y problem).
  - sfa/analysis/freq_phase.c (existing phase extraction via atan2 on
    carrier; autocorrelation for frequencies — here re-interpreted as
    possible model of the projection step, not post-hoc diagnostics).
  - v2/field_knot and hopfion work (phase arg(Φ), winding, Hopf fibration
    as structure generators — precursors to projection ontology).

Run: python3 01_toy_frequency_projection.py
"""

import cmath
import math
import numpy as np
from typing import Tuple, List

# -----------------------------------------------------------------------------
# Fixed algebra invariants (from v59 theorem-grade results, no free params)
# -----------------------------------------------------------------------------
DIM_G2 = 14
DIM_SPIN7 = 21
Q_ALGEBRA = DIM_G2 / DIM_SPIN7          # 2/3 exactly
UNIVERSAL_DEVIATION = 28.0 / 3.0        # (1-Q)D for D=28,35,63 sectors
SELECTION_DIMS = (28, 35, 63)           # L, F, L⊕F

# ω for Z3 DFT
OMEGA = cmath.exp(2j * math.pi / 3)
OMEGA2 = OMEGA * OMEGA

def dft_z3(v: np.ndarray) -> np.ndarray:
    """
    DFT on Z3 (length-3 real or complex vector).
    F_k = (1/√3) Σ_j v_j ω^{k j}
    Returns complex array [F0 (DC), F1, F2].
    """
    N = 3
    F = np.zeros(3, dtype=complex)
    for k in range(N):
        s = 0.0 + 0.0j
        for j in range(N):
            s += v[j] * (OMEGA ** (k * j))
        F[k] = s / math.sqrt(N)
    return F

def idft_z3(F: np.ndarray) -> np.ndarray:
    """Inverse DFT (for round-trip sanity checks)."""
    N = 3
    v = np.zeros(3, dtype=complex)
    for j in range(N):
        s = 0.0 + 0.0j
        for k in range(N):
            s += F[k] * (OMEGA ** (-k * j))
        v[j] = s / math.sqrt(N)
    return v

# -----------------------------------------------------------------------------
# Brannen kernel reference form (for comparison only; not inserted into model)
# -----------------------------------------------------------------------------
def brannen_masses(a: float, t: float, phi: float) -> np.ndarray:
    """
    Brannen kernel (one conventional form used in v59-era work):
      m_k = a² |1 + t exp(i (φ + 2π k/3))|²   for k=0,1,2
    At equilibrium t²=1/2 and φ=2/9 one recovers Q=2/3 (phase-independent
    for the sum at this t; phase affects individual ratios and is pinned
    by precision data).
    This function is used *only* to generate target numbers for comparison.
    The projection model below must not use these formulas or the value 2/9.
    """
    ks = np.arange(3)
    z = 1.0 + t * np.exp(1j * (phi + 2 * math.pi * ks / 3))
    return (a * a) * np.abs(z) ** 2

def compute_q(m: np.ndarray) -> float:
    """Koide Q from three masses."""
    s1 = np.sum(np.sqrt(m))
    s2 = np.sum(m)
    return s2 / (s1 * s1) if s1 > 0 else 0.0

# -----------------------------------------------------------------------------
# The projection hypothesis (the thing under test)
# -----------------------------------------------------------------------------
def phase_projection_hypothesis(
    v: np.ndarray,
    weight_mode: str = "g2_content"
) -> Tuple[float, np.ndarray, float]:
    """
    Candidate "frequency + phase projection" operator.

    Inputs:
      v : length-3 real vector (raw octo-info on Z3 generation cycle).
          Must NOT already encode φ=2/9 or Brannen eigenvalue ratios.
      weight_mode : how to weight the frequency components using algebra
                    invariants. Options:
                      "g2_content"  — weight by Q_ALGEBRA = 14/21 (DC vs. ±1 modes)
                      "deviation"   — weight using UNIVERSAL_DEVIATION = 28/3
                      "selection"   — weight using selection dims (28,35,63)
                      "none"        — unweighted (sanity / circularity check)

    Projection rule (deliberately simple, to be non-circular):
      1. Compute DFT F = [F0, F1, F2].
      2. Apply algebra-derived weights w_k to the frequency modes.
      3. Form a complex "projected amplitude" A = Σ_k w_k * F_k * (phase factor).
      4. Extract effective phase φ_eff = arg(A)  (or a skewness measure on Re/Im).
      5. Form effective masses m_k from |F_k| or Re/Im projections after
         the weighting (low-frequency part after "filtering").
      6. Compute Q_eff from the m_k.

    The critical requirement: the weighting w_k and the phase-extraction
    formula must be built only from fixed algebra numbers (14,21,28,35,63,
    28/3, etc.). No free parameters that can be tuned to insert 2/9.

    Returns:
      phi_eff, m_eff, Q_eff
    """
    F = dft_z3(v)

    if weight_mode == "g2_content":
        # DC (k=0) gets weight 1; ±1 modes get weight proportional to G2 content
        # or to the ratio that already produces Q.
        w = np.array([1.0, Q_ALGEBRA, Q_ALGEBRA], dtype=float)
    elif weight_mode == "deviation":
        # Weight higher modes by the universal deviation constant
        w = np.array([1.0, UNIVERSAL_DEVIATION / 10.0, UNIVERSAL_DEVIATION / 10.0])
    elif weight_mode == "selection":
        # Weight by normalized selection dims (toy extension to "sectors")
        w = np.array([1.0, 28.0/63, 35.0/63])
    elif weight_mode == "none":
        w = np.ones(3)
    else:
        raise ValueError(f"unknown weight_mode: {weight_mode}")

    # Normalize weights to avoid trivial scaling
    w = w / np.sum(w)

    # Phase extraction: form a weighted complex sum from the frequency modes.
    # Use the imaginary part of the +1 mode (the "rotating" frequency) as the
    # source of phase, weighted by the algebra factor.
    # This is *not* atan2 on a pre-existing carrier; it is built from F.
    A = w[0] * F[0] + w[1] * F[1] + w[2] * F[2]

    # Effective phase: argument of the weighted rotating component,
    # folded into [0, 2π/3] to compare with the 3-cycle.
    phi_eff = math.atan2(A.imag, A.real) % (2 * math.pi)

    # Effective masses: |projection| of each frequency component after weighting.
    # (Interpretation: each frequency band contributes a mass-like density.)
    m_eff = np.abs(w * F) ** 2
    # Normalize so sum(m) is order-1 for comparison with Brannen a=1 cases
    if np.sum(m_eff) > 0:
        m_eff = m_eff / np.sum(m_eff)

    Q_eff = compute_q(m_eff)

    return phi_eff, m_eff, Q_eff

# -----------------------------------------------------------------------------
# Experiment harness
# -----------------------------------------------------------------------------
def run_experiment(v: np.ndarray, label: str) -> None:
    print(f"\n=== Experiment: {label} ===")
    print(f"Input v = {v}")
    print(f"  (sum v = {np.sum(v):.6f}, not pre-loaded with Brannen structure)")

    for mode in ["g2_content", "deviation", "selection", "none"]:
        phi, m, Q = phase_projection_hypothesis(v, weight_mode=mode)
        # Compare to target 2/9 ≈ 0.6981 rad (within one 3-cycle)
        target_phi = 2.0 / 9.0
        # The DFT phase lives in [0, 2π); fold target for comparison
        delta_phi = min(abs(phi - target_phi), abs(phi - (target_phi + 2*math.pi)))
        print(f"  [{mode:12s}] φ_eff = {phi:8.5f} rad  "
              f"(target 2/9≈0.69813, Δ={delta_phi:7.5f})  "
              f"m={m}  Q_eff={Q:.5f} (target 2/3≈0.66667)")

def main():
    print("=" * 78)
    print("01_toy_frequency_projection — octo-space frequency/phase projection toy")
    print("=" * 78)
    print()
    print("Fixed algebra invariants (v59 theorem-grade, no free parameters):")
    print(f"  dim G2 = {DIM_G2}, dim Spin(7) = {DIM_SPIN7}, Q = {Q_ALGEBRA:.10f}")
    print(f"  universal deviation (1-Q)D = {UNIVERSAL_DEVIATION:.10f}")
    print(f"  selection dims = {SELECTION_DIMS}")
    print()
    print("Reference Brannen (for comparison only; NOT inserted into model):")
    a, t, phi_ref = 1.0, math.sqrt(0.5), 2.0/9.0
    m_ref = brannen_masses(a, t, phi_ref)
    Q_ref = compute_q(m_ref)
    print(f"  a=1, t=√0.5, φ=2/9 → m_ref = {m_ref}  Q={Q_ref:.10f}")
    print()

    # Test vectors: simple, symmetric, or random — none encode 2/9 or Brannen ratios.
    # Case 1: uniform (DC only in frequency space)
    v_uniform = np.array([1.0, 1.0, 1.0])

    # Case 2: mildly asymmetric (realistic "info" variation across generations)
    v_asym = np.array([1.0, 0.9, 0.8])

    # Case 3: random (within reason)
    rng = np.random.default_rng(42)
    v_random = rng.uniform(0.5, 1.5, 3)

    # Case 4: one component dominant (extreme hierarchy test)
    v_dom = np.array([1.0, 0.1, 0.01])

    for v, label in [
        (v_uniform, "uniform (pure DC)"),
        (v_asym, "mildly asymmetric"),
        (v_random, "random (seed 42)"),
        (v_dom, "one dominant"),
    ]:
        run_experiment(v, label)

    print()
    print("=" * 78)
    print("Analysis notes (see 01_findings.md for full discussion):")
    print("  - If all modes produce φ_eff far from 2/9 *and* Q_eff far from 2/3,")
    print("    this formalization of the projection does not recover the targets.")
    print("  - If 'none' (unweighted) recovers them but weighted modes do not,")
    print("    the algebra invariants are not helping — possible circularity flag.")
    print("  - If a weighted mode recovers φ≈2/9 and Q≈2/3 from non-special v,")
    print("    that is a positive signal for this class of projection operator.")
    print("=" * 78)

if __name__ == "__main__":
    main()
