# Diagnostics and Constraint Layer (Placeholder)

**Date**: 2026-05-18  
**Relation to LAYER_CRITIQUE.md**: Supports Route 4 (reverse-engineering) and provides measurement infrastructure applicable to all routes, including pre-geometric and two-layer explorations.

## Purpose

This folder collects tools, thought experiments, and analysis approaches for:

- Extracting the effective metric or propagation characteristics that a candidate medium + excitation system actually produces.
- Testing universality, non-dispersion, equivalence, and long-range force laws against the quantitative requirements in `../first_principles/EXPECTED_BEHAVIOR.md`.
- Defining operational notions of "density," "local causal structure," and "tick rate" that can be computed from simulation output or algebraic data without assuming a background manifold.
- Cross-checking any emergent geometry against precision observables (GW speed coincidence, Shapiro delay, light deflection, von Laue stability, etc.).

## Status in v58

Light placeholder. No new analysis code or SFA readers yet. Existing project tools (sfa/ analysis suite, bin/ utilities, v39 BLV solver, v43/v51 depletion profilers, v57 mechanical analysis) remain available for targeted re-use on any new prototypes.

## Guiding Principle

Diagnostics must be designed *after* the physics requirements are stated, not before. The first task is to ensure that any measurement we perform is capable of falsifying a model against the list in EXPECTED_BEHAVIOR.md (non-dispersive to 10^{-15}, universal, 1/r² limit, realistic composite pressure profiles, etc.).

See `../DIRECTION.md` and the first-principles folder for the current target specification. Concrete diagnostic designs will be added only when a substrate proposal exists that can be run and measured.