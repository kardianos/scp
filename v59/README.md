# v59 — Kepler Program: Tight Numerical Mysteries → Octonionic Standard Model

**Date**: 2026-05-22
**Status**: Active. Lepton-sector Kepler ellipse achieved (Brannen + Koide + cross-sector ratio structurally derived from $G_2 \subset \text{Spin}(7) \subset \text{Spin}(8)$). Furey-style full Standard Model construction now in progress.

---

## TL;DR

Three independent empirical residues in charged-lepton physics have been collapsed into one algebraic statement:

| Quantity | Structural identification | Empirical match |
|----------|---------------------------|-----------------|
| Koide ratio Q | $\dim G_2 / \dim \text{Spin}(7) = 14/21$ | $6.2 \times 10^{-6}$ |
| Brannen phase φ | $Q/3 = 2/9$ rad | $7.4 \times 10^{-6}$ |
| α–G hierarchy ratio | $\dim \text{Spin}(7) = 21$ | $2.6 \times 10^{-3}$ |

This is the first quantitative match between the SCP project's algebraic framework and empirical particle physics. **See [SUMMARY.md](SUMMARY.md) for the consolidated result.**

The Furey-style full Standard Model construction is being executed in [`furey_construction/`](furey_construction/) — see [`furey_construction/PLAN.md`](furey_construction/PLAN.md) for the multi-variant approach.

## Authoritative References

| File | Purpose |
|------|---------|
| [`SUMMARY.md`](SUMMARY.md) | Headline results across all steps (READ FIRST) |
| [`TYCHO_TABLE.md`](TYCHO_TABLE.md) | Curated experimental targets (PDG/CODATA values) |
| [`first_experiments/`](first_experiments/) | Steps 1–4: numerological scans (Koide, 2/9 rad, modular/hyperbolic, m_μ/m_e) |
| [`multivector_kernel_fit/`](multivector_kernel_fit/) | Steps 5–8: Cl(3,1) Z₃ kernel + constraint surface, α–G corroboration |
| [`octonionic_extension/`](octonionic_extension/) | Steps 9–11: dim Spin(7) = 21, Brannen from triality, instanton actions |
| [`furey_construction/`](furey_construction/) | Step 12+: full ℂ⊗ℍ⊗𝕆 construction (active) |

## Purpose

Identify, from established particle physics, the tightest dimensionless numerical relations that any viable particle-from-field substrate must reproduce. Then test candidate substrates by direct computation against those numbers. Substrates that cannot reproduce a single tight number with zero or one free parameter are ruled out. Substrates that reproduce one or more tight numbers become Kepler ellipses — empirical fits demanding explanation.

## Relation to Previous Work (v28–v58)

- [`../v58/LAYER_CRITIQUE.md`](../v58/LAYER_CRITIQUE.md) — fixed-grid soliton models cannot produce gravity as an emergent geometric effect.
- [`../v58/first_principles/EXPECTED_BEHAVIOR.md`](../v58/first_principles/EXPECTED_BEHAVIOR.md) — universal/non-dispersive propagation requirements.
- [`../v58/pregeometric/PARTICLES_AS_DENSITY_ACHIEVERS.md`](../v58/pregeometric/PARTICLES_AS_DENSITY_ACHIEVERS.md) — particles as density-elevation mechanisms.
- [`../v58/pregeometric/unified_multivector_force/`](../v58/pregeometric/unified_multivector_force/) — qualitative Newton+Maxwell recovery on 2D lattice; no quantitative TYCHO match.
- [`../v57/6FIELD_MECHANICAL_ASSESSMENT.md`](../v57/6FIELD_MECHANICAL_ASSESSMENT.md) — quantitative proton mechanical-structure failure of 6-field Cosserat.

## Why Kepler Before Newton

Kepler did not derive the ellipse from physics. He fit it to Tycho's observations and showed the fit was tight enough to demand an explanation. Newton's inverse-square law came afterward. Without Kepler's fit, Newton had no anchor.

The SCP project (v28–v58) accumulated rich phenomenology but no quantitative particle-physics match. The Kepler-stage approach: lock the targets first, then derive the structure that hits them.

## Methodology

A candidate substrate is tested by:

1. Enumerating its **natural invariants** — eigenvalues, topological invariants, mode spectra. Natural = not constructed to fit.
2. Selecting subsets in algebraically motivated ways.
3. Computing the **target relation** (Koide, Brannen, mass ratios) for that subset.
4. Comparing to the experimental value in `TYCHO_TABLE.md`.

A substrate passes a target if it reproduces the experimental value to within experimental precision (or 1% on first pass). It passes the whole program if it reproduces at least one Tier-1 target with ≤1 free parameter AND does not grossly violate any other.

## Falsifiability

**Ruled out** when no natural invariant subset satisfies any Tier-1 relation in `TYCHO_TABLE.md` better than chance.

**Promoted to Kepler ellipse** when it reproduces a Tier-1 target to experimental precision with ≤1 free parameter (preferably 0).

## Current Status (after 11 steps)

| Item | Status |
|------|--------|
| **Brannen form (eigenvalues at 120°)** | **Structural** (Cl(3,1) Z₃) |
| **Koide identity Q = 2/3** | **Structural** ($G_2/\text{Spin}(7)$) |
| **Brannen phase φ = 2/9 rad** | **Structural** ($Q/3$) |
| **Cross-sector ratio 21** | **Structural** ($\dim \text{Spin}(7)$) |
| **Three generations** | **Structural** (Z₃ ⊂ S₃ triality of Spin(8)) |
| Lepton mass scale $a$ | Empirical |
| α (absolute) | Empirical; $\pi^2/2$ conjecture (0.30% gap) |
| G (absolute) | Empirical; $21\pi^2/2$ conjecture (0.55% gap) |
| Full SM embedding | In progress ([`furey_construction/`](furey_construction/)) |

## What This Version Does Not Do

- No new fixed-grid simulations.
- No claims of substrate validity until a Tycho-table number is reproduced with ≤1 free parameter.
- No theory-building before fitting.

## Posture

A numerological program is honest only if it tracks failures as carefully as successes. Every substrate tested has its result documented — whether near-miss, clean rejection, or hit. A failed Kepler is more valuable than a successful epicycle.

The current state: **lepton-sector Kepler ellipse achieved**. Cross-sector unification (full SM + gravity) in progress.
