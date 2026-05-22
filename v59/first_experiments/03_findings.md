# Experiment 03 — Findings: Modular and Hyperbolic Scan for 2/9 rad

**Date**: 2026-05-22
**Script**: `03_modular_hyperbolic_scan.py`
**Target**: TYCHO_TABLE T1.2 (Brannen phase ≈ 2/9 rad) in less-standard mathematical structures, following the negative result from experiment 02.

---

## Summary

**Second consecutive null result for the 2/9-rad lead.** No natural value among Lobachevsky function ratios, small hyperbolic 3-manifold volumes (or arctan/log/ratios of them), Eisenstein / Dedekind η evaluations at CM points, Bloch-Wigner D(z) at natural arguments, or Rogers-Ramanujan q-series ratios matches 2/9 rad to within 10⁻⁴. A few values land within 10⁻² (close) but not within 10⁻⁴ (tight).

Combined with experiment 02, the 2/9 rad lead has now been tested against:
- Coxeter angles π/n, 2π/n
- Platonic and exceptional polytopes
- Lie algebra root systems (A_n n≤5, B_n n≤3, D_n n≤5, G₂, F₄, E₈)
- arctan/arcsin/arccos of all rationals p/q ≤ 49
- Simple algebraic combinations of √, ln, π, small integers
- Lobachevsky function values and ratios
- Hyperbolic 3-manifold volumes (Weeks, figure-8, Whitehead, Meyerhoff, 5_2, 6_1, Borromean, Bianchi)
- Dedekind η at CM points i, ω, i√2, (1+i√3)/2, (1+i√7)/2, i√3, (1+i√11)/2
- Bloch-Wigner D(z) at 10 natural complex arguments
- Rogers-Ramanujan q-series at multiple q values

**Nothing matches.** This is conclusive within the scope tested.

## Notable Near-Misses (>10⁻⁴ but <10⁻²)

| Expression                  | Value      | Δ from 2/9      |
|-----------------------------|------------|------------------|
| arctan(√(m_μ/m_τ))          | 0.23918    | +0.0170          |
| π/14                        | 0.22440    | +0.0022          |

These do not rise to tight-hit status.

## Notable Algebraic Observations (Side Findings)

- **arg(η((1+i√3)/2)) = π/24 = 0.13090** — exactly. This is the standard η transformation phase under T (τ → τ+1), giving e^(iπ/12) which has argument π/12 on η or π/24 on η/q^(1/24). Beautiful but unrelated to 2/9.
- **D(ω) = D(e^(2πi/3)) = 0.67663 = Cl₂(2π/3)** — Catalan-Clausen function value. Unrelated to 2/9.
- All tested D(z) values fall in the range [0, 1.015] = [0, V(ideal regular tetrahedron)], as expected since D is bounded by the maximal hyperbolic-tetrahedron volume.

These are sanity-check confirmations that the computation is correct; they are not 2/9 hits.

## What This Means for the 2/9 Hypothesis

Recall the two hypotheses from experiment 02:
- **H1 (coincidence)**: φ is close to but not exactly 2/9. Future m_τ improvements will reveal the deviation.
- **H2 (subtle origin)**: φ = 2/9 exactly, with origin in non-standard mathematics.

After three experiments and exhaustive scanning across standard and non-standard structures, **H2 is increasingly disfavored.** No natural mathematical object across modular, hyperbolic, Lie-algebraic, polytopal, or trigonometric structures yields 2/9 rad. If H2 is correct, the origin lies in something genuinely outside the standard mathematical toolkit (possibly: dynamics of a specific non-perturbative quantization, a special property of an as-yet-unidentified substrate, or higher arithmetic).

**Current best assessment**: φ = 2/9 is likely a coincidence at the current ~10⁻⁵ precision of m_τ. The next decimal place of m_τ — improvable by Belle II and tau-charm factories over the next decade — will likely show φ departing from exactly 2/9 at the seventh or eighth digit. Until that data arrives, the 2/9 lead is **exhausted** as an actionable target.

## Recommendation: Pivot to T1.5

The 2/9 investigation has produced strong negative information — pruning a wide swath of substrate candidates — but no positive Kepler ellipse. The cost-benefit balance favors moving to a different Tier-1 target while keeping 2/9 in reserve.

Next target: **T1.5 — muon-electron mass ratio = 206.7682830(46)**.

Reasons:
- Lower bar: only two natural invariants needed (no triple to align, no phase to extract).
- Different scale: 200× ratio is not a "small angle" problem; it's a ratio problem with much wider variation among candidates.
- Cross-checks: many ways to get an integer-like ratio (group orders, quantization indices, knot crossings, etc.).
- If a substrate passes T1.5, it can be extended to test T1.1 + T1.2 simultaneously.

## What Remains Worth Trying for 2/9 (Future)

Listed for the record in case the lead becomes actionable later:

1. **Atomic clocks or fundamental constants combinations**: e.g., 1/(2π × 9), small fractions of fundamental geometric phases in Mott scattering, etc.
2. **Specific RG flow trajectories**: scan QFT beta-function flows for ratios passing through 2/9 at some scale.
3. **Special angles of arithmetic 3-orbifolds**: angles in commensurability classes of Bianchi groups.
4. **Critical lines of L-functions**: imaginary parts of zeros of Riemann ζ, Dirichlet L-functions, modular L-functions — search for 2/9 rad among them.
5. **Specific cocycle values**: H³(SO(3); Z) cocycles evaluated on natural 3-cells.

Each is a deep rabbit hole. None should be pursued without independent motivation (e.g., a substrate emerging from another experiment that has 2/9 as a natural angle).

## Status After 3 Experiments

| Target | Status |
|--------|--------|
| T1.1 Koide Q = 2/3 | Confirmed empirically (~10⁻⁵). No substrate that gives |b|/a = 1/√2 from first principles. |
| T1.2 Brannen φ ≈ 2/9 rad | Confirmed empirically (consistent with exact 2/9 within m_τ uncertainty). **No structural origin identified across exhaustive scans.** Likely coincidence pending tighter m_τ. |
| T1.3 m_p/m_e | Not yet tested. |
| T1.4 α⁻¹ | Not yet tested. |
| T1.5 m_μ/m_e | Not yet tested. Next target. |

The numerological program has been rigorous about ruling things out. We have not yet found a Kepler ellipse, but we have decisively eliminated several wide classes of candidate substrates. That is the correct shape of a Kepler-stage search.
