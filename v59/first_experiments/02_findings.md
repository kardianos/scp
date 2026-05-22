# Experiment 02 — Findings: 2/9 rad Geometric Audit

**Date**: 2026-05-22
**Script**: `02_two_ninths_audit.py`
**Target**: TYCHO_TABLE T1.2 (Brannen phase ≈ 2/9 rad). Find a natural geometric or algebraic object whose canonical angle equals 2/9 rad to within 10⁻⁴.

---

## Summary

**Negative result with sharp consequence.** The Brannen phase φ = 2/9 rad is consistent with experiment (Part A), but it does not arise as the natural angle of any tested standard structure: no Coxeter angle π/n, no Platonic-solid vertex angle, no Lie-algebra root angle (A₂..A₅, B₂, B₃, D₄, D₅, G₂, F₄, E₈), no arctan/arcsin/arccos of any rational p/q with p,q ≤ 49.

This is informative: it rules out a wide class of substrate candidates and narrows the search to less-standard mathematical structures (modular forms, hyperbolic geometry, Berry phases on non-trivial bundles, RG flow trajectories) or to coincidence (φ converging to a value other than exactly 2/9 as m_τ measurement improves).

## Part A: φ = 2/9 vs. Experimental Uncertainty

Sensitivity of the Brannen-fit phase φ to the m_τ measurement:

| m_τ (MeV)   | n σ  | φ (rad)         | φ − 2/9     | |b|/a effective |
|-------------|------|-----------------|-------------|-----------------|
| 1776.50     | −3   | 0.2222547       | +3.25 × 10⁻⁵ | 0.9999603        |
| 1776.62     | −2   | 0.2222463       | +2.41 × 10⁻⁵ | 0.9999704        |
| 1776.74     | −1   | 0.2222380       | +1.58 × 10⁻⁵ | 0.9999806        |
| 1776.86     |  0   | 0.2222296       | +7.41 × 10⁻⁶ | 0.9999908        |
| 1776.98     | +1   | 0.2222213       | −0.94 × 10⁻⁶ | 1.0000009        |
| 1777.22     | +3   | 0.2222046       | −1.76 × 10⁻⁵ | 1.0000213        |

φ = 2/9 exactly corresponds to m_τ ≈ 1776.98 MeV — about +1 σ above the PDG central value. **This is well inside experimental tolerance.** The φ = 2/9 hypothesis is empirically tenable.

(One row at m_τ = 1777.10 produced a phi = 1.872 rad in the output; that is a permutation-labeling ambiguity in the least-squares fit, not a physical jump. The trend across the other rows is monotonic and clean.)

## Part B: Coxeter Angles π/n and 2π/n

Scanned n = 2..99 for π/n and 2π/n. Three entries land within 10⁻² of 2/9:

| Expression  | Value     | Δ from 2/9      |
|-------------|-----------|------------------|
| π/14        | 0.224399  | +2.18 × 10⁻³     |
| 2π/28       | 0.224399  | +2.18 × 10⁻³     |
| 2π/29       | 0.216662  | −5.56 × 10⁻³     |

**No tight hits.** π/14 is the closest at the 0.2% level — meaningfully different from 2/9 within experimental φ precision.

## Part C: arctan / arcsin / arccos of Rationals

Scanned all p/q with 0 ≤ p ≤ 49, 1 ≤ q ≤ 49 for arctan, arcsin, arccos. **No tight hit.**

Notable:
- arctan(2/9) = 0.2187 rad (Δ −3.6 × 10⁻³) — small-angle approximation off by 0.4%.
- arcsin(2/9) = 0.2241 rad (Δ +1.9 × 10⁻³) — likewise.

The number 2/9 is not the trigonometric image of any small-rational argument at tight precision.

## Part D: Platonic Solid Vertex Angles

Computed all distinct pairwise vertex angles for each Platonic solid:

| Solid          | Vertices | Distinct angles | Nearest to 2/9          |
|----------------|----------|-----------------|--------------------------|
| Tetrahedron    | 4        | 1               | 1.911 rad (Δ +1.69)      |
| Cube           | 8        | 3               | 1.231 rad (Δ +1.01)      |
| Octahedron     | 6        | 2               | 1.571 rad (Δ +1.35)      |
| Icosahedron    | 12       | 3               | 1.107 rad (Δ +0.88)      |
| Dodecahedron   | 20       | 5               | 0.730 rad (Δ +0.51)      |

**No hit anywhere.** Platonic solids have vertex angles in the range ~30° to ~110°, well above 12.7°. They do not produce small angles like 2/9 rad.

## Part E: Lie Algebra Root Angles

Computed all distinct pairwise root angles for A₂, A₃, A₄, A₅, B₂, B₃, D₄, D₅, G₂, F₄, E₈:

| Algebra | # roots | # distinct angles | Nearest to 2/9              |
|---------|---------|-------------------|------------------------------|
| A_n (n=2..5) | 6–30  | 3–4              | 60° = π/3 (Δ +47°)            |
| B₂, B₃, F₄  | 8–48 | 5–7              | 45° = π/4 (Δ +32°)            |
| D₄, D₅, E₈  | 24–240| 4                | 60° = π/3 (Δ +47°)            |
| G₂           | 12   | 6                | 30° = π/6 (Δ +17°)            |

**No hit.** All standard simply-laced Lie algebras have root angles that are multiples of π/3 (60°) and π/6 (30°). Non-simply-laced ones add π/4 (45°). None has a natural 12.7° angle.

This rules out a large class of "geometry of the gauge group" substrates.

## Part F: Algebraic Candidate Expressions

| Expression                          | Value           | Δ from 2/9       |
|-------------------------------------|-----------------|-------------------|
| 2/9 (target)                        | 0.2222222       | 0                 |
| arccos(1 − 2/81) ← Taylor 2/9       | 0.2226820       | +4.60 × 10⁻⁴      |
| ln(5/4)                             | 0.2231436       | +9.21 × 10⁻⁴      |
| arcsin(2/9)                         | 0.2240931       | +1.87 × 10⁻³      |
| arctan(2/9) = arctan(1/4.5)         | 0.2186689       | −3.55 × 10⁻³      |
| 2/(3 + √45)                         | 0.2060113       | −1.62 × 10⁻²      |

`arccos(1 − 2/81)` is the second-order Taylor expansion of cos(2/9) — it agrees with 2/9 only because of the small-angle Taylor series, not from any deeper relation. The miss at 4.6 × 10⁻⁴ confirms it is just the Taylor-series approximation, not an exact identity.

`ln(5/4) ≈ 0.2231` was checked because of the suggestive form (small log ratio); it misses by 10⁻³.

None of the algebraic expressions tested produces 2/9 exactly except by definition.

## Interpretation

The empirical Brannen phase φ for charged leptons is consistent with 2/9 rad to within current m_τ precision (~7 × 10⁻⁶). This is suggestive. But **2/9 rad is not a natural angle of any standard geometric or algebraic structure**:

- Not a Coxeter angle.
- Not an angle of any Platonic solid.
- Not a Lie-algebra root angle (rank ≤ 8 simply-laced, plus B, F, G).
- Not an arctan/arcsin/arccos of any low-denominator rational.
- Not a simple algebraic combination of √, π, ln of small integers.

## Two Hypotheses Going Forward

**H1 (Coincidence)**: φ is close to but not exactly 2/9. The next significant figure of m_τ — improvable by ongoing BES III, Belle II, and future tau-charm factories — will reveal whether the equality holds. If φ departs from 2/9 at the seventh or eighth digit, the apparent rational is a small-number coincidence, and we drop this lead.

**H2 (Subtle origin)**: φ = 2/9 exactly, and its origin is a structure not covered by standard polytopes, Lie algebras, or simple rationals. The candidate sources are:

1. **Modular / elliptic structures**: special values of modular forms, j-invariant, Eisenstein series, theta functions. φ could be a critical angle of a specific modular curve.
2. **Hyperbolic geometry**: closed-geodesic lengths or hyperbolic-triangle angles on a specific quotient surface (e.g., the modular surface H/SL(2,Z), or a 3-manifold like the Weeks manifold).
3. **Berry / Aharonov-Bohm phase**: a topologically nontrivial closed loop in some parameter space yields exactly 2π × (1/9π) = 2/9 rad of geometric phase. The 1/(9π) suggests a 9-fold cyclic structure with an irrational normalization.
4. **RG flow fixed-point ratio**: an RG trajectory in a coupling space passes through a fixed point where some ratio of couplings is exactly 2/9. This is reminiscent of how the Weinberg angle's GUT value 3/8 emerges.
5. **Non-perturbative quantization condition**: φ is the quantization of a topological action divided by ℏ in a substrate that produces ℏ as an emergent quantity.

## What the Audit Rules Out (Substrate-Side)

Substrates that rely solely on:
- Vertex / edge / face geometry of Platonic solids (or compounds).
- Root systems of low-rank simply-laced Lie algebras.
- Coxeter-group reflection angles.
- Simple rationals or small-radical combinations.

...cannot produce φ = 2/9 from their natural invariants alone. They could still be valid substrates with additional dynamical structure layered on top, but they will not predict 2/9 by structure alone.

This is a substantial pruning of the candidate-substrate space.

## Recommended Next Probes

In priority order, each ~one focused session:

1. **Modular-form audit**: scan special values (zeros, j-invariant residues, Eisenstein critical points) of modular forms over SL(2,Z) and small congruence subgroups for natural angles near 2/9 rad.

2. **Hyperbolic 3-manifold audit**: compute volumes, geodesic lengths, and cusp angles of small hyperbolic 3-manifolds (figure-8 complement, Weeks, etc.) and look for 2/9 rad in their invariants. Use SnapPy data (publicly available) if computational tools allow.

3. **Berry-phase scan**: for SU(2) coherent-state evolution along a closed loop on S² traced by a 3-fold-symmetric soliton, compute the Berry phase as a function of loop size. Find the loop that gives 2/9 rad.

4. **Switch to T1.5 (m_μ/m_e ratio)**: rather than chase 2/9 further, attack the cleaner two-mode test. The lower bar (only two natural invariants needed, no triple to align) may surface a substrate that the Koide triple search has missed.

5. **Refine m_τ**: track ongoing tau lepton mass measurements (BES III runs, Belle II) and re-run the Brannen fit when m_τ uncertainty improves by an order of magnitude. If φ stays at 2/9 to 7+ digits, H1 (coincidence) is ruled out and we commit to H2.

## Status of Tycho Hits So Far

After two experiments:

- **T1.1 (Koide Q = 2/3)**: numerically confirmed to 6 ppm in lepton data, naturally produced by any Hermitian Z₃-cyclic 3×3 matrix with |b|/a = 1/√2. **No substrate identified that gives |b|/a = 1/√2 from its own structure.**
- **T1.2 (Brannen φ ≈ 2/9 rad)**: empirically consistent with exactly 2/9. **No standard substrate gives this angle.** Either coincidence or a non-standard mathematical origin.

The numerological program is on solid footing (confirmed precision of the targets) but has not yet produced a Kepler ellipse — a substrate-derived match. The audit has, however, ruled out the most obvious substrate classes, which is exactly what a Kepler-stage program is supposed to do.
