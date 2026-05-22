# v58 — Post-Layer-Critique Exploration: First Principles and Pre-Geometric Directions

**Date**: 2026-05-18  
**Status**: Initial scaffolding. No kernels, no simulation code, no data products. Pure framing of the physics target and conceptual entry points.

---

## Purpose

This version directory opens the post-v57 phase following the structural diagnosis in [LAYER_CRITIQUE.md](./LAYER_CRITIQUE.md). The prior framework (6-field Cosserat on fixed grid, later GA multivector extensions, v28–v57) successfully explored soliton forces, depletion response, phase confinement, and binding analogs *inside* a passive fixed-background medium. It cannot produce gravity as an emergent geometric property of the medium's own state.

v58 deliberately shifts the layer:

- Lead with **first-principles derivation** of the *expected behavior from physics* under the hypothesis that (a) the fundamental entity is a dynamical field medium, (b) particles are its localized, stable excitations, and (c) gravity is a modulation of effective propagation / causal structure arising from the medium's density or state gradients (Route 4 reverse-engineering, per LAYER_CRITIQUE).
- Then explore whether a **pre-geometric** (combinatorial, algebraic, relational) substrate can naturally realize that behavior without presupposing the metric or causal structure it must produce (Route 3, the conceptually strongest long-term direction).

The split prioritizes physics requirements over implementation substrate. Any concrete model (pre-geometric, two-layer medium, or other) will later be judged by whether it meets the requirements articulated in `first_principles/EXPECTED_BEHAVIOR.md`.

## Directory Layout

- `LAYER_CRITIQUE.md` — hard-linked diagnostic (the "why" for the reset).
- `DIRECTION.md` — rationale for prioritization, folder structure, and decision to defer code.
- `first_principles/` — central requirements document and supporting notes. Leads the exploration.
- `pregeometric/` — framing and roadmap for algebraic/relational starting points.
- `two_layer/` — light placeholder for the pragmatic explicit-medium route (Route 1).
- `diagnostics/` — light placeholder for measurement, constraint, and validation considerations applicable across routes.

## Key References (Internal)

- LAYER_CRITIQUE.md (this dir and root) — structural incompatibilities (fixed background, passive medium, dispersive propagation, forces vs. geometry).
- v39/ (wavefront/BLV analysis) — demonstrated frequency-dependent propagation effects under amplitude/density variations.
- v43/, v51/ (depletion / gradient response) — genuine classical drift of composites toward self-induced low-density regions via asymmetric V(P); remains a force inside fixed space.
- v57/ (6FIELD_MECHANICAL_ASSESSMENT.md) — quantitative mismatch between 6-field braid mechanical profiles and lattice-QCD proton stress-energy (pressure signs, compactness).
- CONCEPT.md, DISCOVERIES.md (root) — current best understanding and historical record.
- CLAUDE.md — documentation standards (textbook style for theory docs; "it is NOT X, it IS Y"; quantitative targets; clear "what has not been established").

## Posture for Subsequent Work

- No new 6-field, GA-on-grid, or fixed-background kernels in this version or immediate follow-on without explicit justification against the requirements here.
- First-principles requirements are substrate-independent constraints. They must be satisfied before investment in discretization or dynamics.
- Pre-geometric exploration is encouraged even at small scale; a minimal demonstration of emergent causal structure + stable localized excitations would be high-value evidence.
- Two-layer and diagnostic work remain available as bridges or cross-checks but are secondary to the first-principles framing.

This directory is intentionally lightweight. Its value is in the clarity of the target it defines for the next modeling layer.