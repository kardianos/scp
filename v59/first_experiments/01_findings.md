# Experiment 01 — Findings: First Koide Substrate Scan

**Date**: 2026-05-22
**Script**: `01_koide_substrate_scan.py`
**Target**: TYCHO_TABLE T1.1 (Koide) and T1.2 (Brannen phase φ).

---

## Summary

Numerical confirmation of the Brannen parametrization for charged leptons (Q = 2/3 and φ ≈ 2/9 rad both inside experimental uncertainty). No tested substrate naturally produces the required structure. One sharp residue identified: φ is within 7 × 10⁻⁶ of exactly 2/9 rad, demanding explanation.

## Experimental Results (Re-Derivation)

Using PDG 2024 charged lepton masses:

- m_e = 0.51099895069 MeV
- m_μ = 105.6583755 MeV
- m_τ = 1776.86 ± 0.12 MeV

Computed values:
- **Koide ratio Q** = 0.6666605115
- 2/3 = 0.6666666667
- Δ(Q − 2/3) = −6.155 × 10⁻⁶

The deviation from 2/3 is at the 6 ppm level, which is consistent with the 7 × 10⁻⁵ relative uncertainty in m_τ. The Koide identity Q = 2/3 is therefore **experimentally consistent with being exact** at current measurement precision.

**Brannen fit** to (√m_e, √m_μ, √m_τ) with the model √m_k = a(1 + √2 cos(2πk/3 + φ)):
- a = 17.7155617 √MeV
- φ = 0.22222963 rad = 12.7328°
- |b|/a effective = 0.99999077 (1.0 = perfect Brannen)
- fractional residual = 4.4 × 10⁻¹⁷ (machine precision)
- best permutation (k assignment to e, μ, τ) = (2, 0, 1)

The permutation means: at k = 0 the cosine is largest (predicts τ); at k = 1 most negative (predicts e); at k = 2 in between (predicts μ).

## The Sharp Residue: φ ≈ 2/9 rad

Tested simple rational and transcendental candidates against the fitted φ:

| Candidate         | Value          | Δ(φ_exp − cand) |
|-------------------|----------------|-----------------|
| 2/9 rad           | 0.22222222     | +7.41 × 10⁻⁶    |
| 1/4.5 rad         | 0.22222222     | +7.41 × 10⁻⁶    |
| arctan(2/9)       | 0.21867        | +3.56 × 10⁻³    |
| π/14              | 0.22440        | −2.17 × 10⁻³    |
| 1/(2π)            | 0.15915        | +6.31 × 10⁻²    |
| π/13              | 0.24166        | −1.94 × 10⁻²    |

**2/9 rad is closer than every other candidate by ~250×.** The discrepancy 7.4 × 10⁻⁶ is itself consistent with the m_τ measurement uncertainty. Within current experimental tolerance, **φ is consistent with being exactly 2/9 rad**.

### Why this is interesting

- The Koide identity needs Q = 2/3.
- The Brannen phase fits 2/9.
- 2/9 = (2/3)/3.

This is too tight to dismiss. If both are exact, the lepton sector encodes the rational pair (2/3, 2/9) which differ by a factor of 3 (the number of generations). Possible interpretations:

1. **Hidden algebraic relation**: Q and φ are not independent. Some deeper symmetry forces both.
2. **Coincidence at this precision**: φ might depart from 2/9 at the 7th digit as m_τ measurement tightens. Worth re-checking when m_τ is improved.
3. **Substrate signature**: a particular algebraic structure produces (Q, φ) = (2/3, 2/9) naturally.

This is a **Tycho-table-level lock**. It is now a target: a substrate that derives φ = 2/9 from its own structure (rather than fitting) would be a Kepler ellipse.

## Substrate Candidates Tested

### Hermitian Z₃-cyclic matrices (analytic baseline)

The 3-dimensional Hermitian matrix family M = aI + bS + b*S† (S = cyclic shift, b ∈ ℂ) has eigenvalues a + 2|b| cos(2πk/3 + arg b). Koide Q = 2/3 holds iff |b|/a = 1/√2.

**Substrate question reduced to**: what natural object predicts |b|/a = 1/√2 with arg b = 2/9 rad simultaneously?

### Candidates surveyed (all NOT passing without tuning)

| Candidate                          | Natural |b|/a | Matches? |
|------------------------------------|----------|----------|
| Cube edge / face-diag ratio        | 1/√2     | **suggestive** — but identification of "a" vs "|b|" ambiguous |
| Tetrahedron r_in/r_out             | 1/3      | no       |
| Cube r_in/r_out                    | 1/√3     | no       |
| Octahedron r_in/r_out              | 1/√3     | no       |
| Dodecahedron r_in/r_out            | ≈ 0.60   | no       |
| Cl(3,0) vector/bivector magnitudes | 1        | no       |
| Tight-binding hopping/on-site      | free     | requires tuning |

None of these naturally produces both the ratio AND the phase.

### Knot ropelengths (Pieranski/Cantarella data)

Tested all 56 triples among 8 prime knots up to crossing number 7. **Best Q = 0.354**, distance from 2/3 = 0.31 — not even close. All ropelengths sit in the narrow range 16–31, so triples cluster near Q = 1/3 (equal-magnitude limit).

**Conclusion**: ropelength is not the right knot invariant for sqrt(mass) under a linear identification. If knots are particles, the mass must depend non-linearly on knot topology (e.g., exponentially in some invariant, or polynomially in crossing number).

### Polyhedral vertex projections

Cube and icosahedron vertices projected on natural 3-fold axes give triples that are either fully degenerate (cube: three vertices at +1/3) or too few distinct values (icosahedron has only 4 distinct projection magnitudes). No Koide hits.

## What the Scan Did NOT Test (Yet)

- **Finite group representations**: A₄, S₄, A₅ Casimir spectra and irrep dimensions.
- **Octonions / Cl(0,7)**: Furey-program eigenvalue structures.
- **Exceptional Lie group root ratios**: G₂, F₄, E₆, E₇, E₈.
- **Non-linear functions of knot invariants** (exponential ropelength, polynomial of crossing number, Jones polynomial values at roots of unity).
- **Substrates with intrinsic SU(2) or SU(3) action** that could naturally produce Z₃-cyclic Hermitian couplings with a small phase.
- **Berry-phase / topological** sources of φ ≈ 2/9 rad on simple geometric objects.

## Recommended Next Probes

In priority order, each ~one Python session:

1. **2/9-rad geometric audit** — exhaustive search for natural occurrences of 2/9 rad (or 2/3 of 1/3 rad, etc.) in standard geometric objects: angles between specific Platonic solid vertices, geodesic lengths on quotient surfaces, Berry phases of small loops in S². If any natural object has exactly this angle, that's the lead.

2. **Finite-group Casimir scan** — eigenvalues of Casimir operators on irreducible representations of A₄, S₄, A₅, and the binary versions (2T, 2O, 2I). The icosahedral group A₅ has order 60 and irreps of dimensions 1, 3, 3, 4, 5 — interesting structure.

3. **Lie-algebra root ratios** — for G₂, F₄, E₆, E₇, E₈, compute the ratios of root lengths and weights. G₂ has rank 2 with a 6+6 root pattern that has natural 3-fold structure.

4. **Substrate with intrinsic chirality** — Hermitian Z₃-cyclic matrices arising from a system with chirality (Berry phase, magnetic flux through a triangular plaquette, SU(2) Wess-Zumino term). Find one whose flux is 2/9 rad and whose hopping/on-site ratio is 1/√2.

## What This Experiment Established

1. **Verified**: Koide Q = 2/3 holds for charged leptons to 6 ppm (experimental limit).
2. **Verified**: Brannen parametrization fits exactly (residual at machine precision).
3. **New numerological residue**: φ ≈ 2/9 rad to 7 × 10⁻⁶. **This is a sharp target.**
4. **Ruled out (without modification)**: knot ropelengths, polyhedral vertex projections, tight-binding without fine-tuning.
5. **Open**: cube edge/face-diagonal identification — suggestive (1/√2 ratio is there) but unclear which is "a" vs "|b|".
6. **Open**: finite-group spectra, octonions, exceptional Lie root systems.

The φ = 2/9 lock is the most actionable single residue. Next session should look exclusively for natural occurrences of 2/9 rad in well-defined geometric/algebraic settings.
