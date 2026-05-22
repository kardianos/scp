# v58 Direction — First-Principles Physics Expectations Before Substrate Choice

**Date**: 2026-05-18  
**Context**: Follows directly from LAYER_CRITIQUE.md (2026-05-18). Records the decision to lead the post-critique phase with explicit derivation of required physical behavior rather than immediate new code or literature-guided models.

---

## 1. The Core Decision

The modeling choice "write nonlinear wave equations for the things we think are inside the field and see if gravity falls out" is structurally incapable of producing gravity as an emergent property of the medium's own state. All accumulated analogs (depletion drift, charge-dependent forces, phase confinement, theta halos) remain forces acting on excitations *inside* a fixed causal structure whose node density, stiffness, or connectivity never dynamically defines the metric experienced by *all* signals.

**Therefore**:

- v58 and the immediate follow-on work begin with **first principles**: articulate, from established gravitational and particle phenomenology plus the hypothesis that particles are field excitations and gravity is density-gradient modulation of effective propagation, what *any* viable realization *must* deliver.
- Only after those requirements are stated rigorously do we evaluate or prototype candidate substrates.
- The first-principles folder leads; the pre-geometric folder follows as the direction that most cleanly avoids reintroducing the fixed-background problem.

This inverts the historical workflow (build a model, measure what forces or propagation effects appear) and replaces it with "state the target phenomenology first, then ask what minimal medium dynamics and excitations can produce it."

## 2. Mapping to LAYER_CRITIQUE Routes

| Description (LAYER_CRITIQUE) | LAYER § | Role in v58 | Priority |
|--------------------------------|----------------|-------------|----------|
| Reverse-engineering from phenomenology (medium properties required for GR + particles) | 2.4 | Implemented as `first_principles/` (EXPECTED_BEHAVIOR.md and follow-ons). Substrate-independent. | **Lead** |
| Pre-geometric / combinatorial / relational | 2.3 | `pregeometric/` (INITIAL_EXPLORATION.md, later roadmaps and small prototypes). | **Second** |
| Explicit two-layer medium + excitations | 2.1 | `two_layer/` (light placeholder). Pragmatic bridge for testing back-reaction and geometric propagation explicitly. | Tertiary |
| Condensed-matter explicit lattice | 2.2 | Not opened (disfavored: presupposes discrete mechanical substrate). | — |
| Literature-guided reset (hedgehog/Skyrme/etc. ignored by design) | 2.5 | Explicitly avoided. First-principles stance maintained. | — |

The synthesis in LAYER_CRITIQUE §3 is followed: use existing tools for targeted diagnostics only; do not refine the old layer as primary model; explore two-layer and Route 4 in parallel; give serious attention to pre-geometric.

## 3. Why First-Principles Leads (User-Explicit Rationale)

"First start with a first principles direction of expected behavior from physics, applied to field as particles, gravity as density gradient."

The user's orbital / tick-rate intuition is central:
- Local microscopic propagation speed (c between adjacent "nodes" or field elements) is constant.
- Effective macroscopic speed or clock rate varies because the *number of nodes* or the *correlation structure per coordinate distance* changes with local field density.
- In a denser region an orbiting clock (any internal oscillator built from field excitations) experiences more phase accumulation or more "hops" per revolution → runs slow relative to a distant observer using light signals.
- The same density field must affect light signals, gravitational waves, and test particles identically (equivalence, universality, non-dispersion).

This picture must be made quantitative against real experiments before any new dynamics are written.

## 4. What v58 Does Not Contain

- No C, CUDA, Go, or Lean kernels.
- No new SFA output or analysis binaries.
- No parameter scans or "try this Lagrangian" experiments.
- No direct porting of old v28–v57 mechanisms without re-examination against the requirements in `first_principles/`.

Code will appear only when a candidate substrate (pre-geometric or otherwise) has been shown, at least conceptually, to be capable of meeting the derived constraints.

## 5. Success Criteria for This Phase

A successful v58 outcome is a clear, citable statement of:

1. The minimal set of behaviors any "field-density → effective geometry" model must reproduce (with quantitative experimental anchors).
2. A well-scoped initial exploration plan for pre-geometric substrates that could satisfy those behaviors by construction.
3. Identification of the highest-leverage open questions (see LAYER_CRITIQUE §4) that the next modeling round must address.

Subsequent versions (v59+) will then prototype realizations measured against this target.

---

*This file is intentionally short. The physics content lives in the sub-folders.*