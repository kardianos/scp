# Phase 4 Completion Criteria — Stability and Self-Regulation

**Phase 4 is considered complete only when ALL of the following conditions are simultaneously satisfied:**

## 1. Coupled Density–Chirality Dynamics on Relational Graphs
- The Python relational framework has been extended (or new code written) to evolve graphs under the exact living candidate while tracking the coupled evolution of local density and protected chirality over time.
- Multiple long runs (on graphs with ≥ 500 nodes) have been performed to observe the long-term behavior of effective geometric quantities (local `d_eff`, isotropy, effective light-cone structure, protected density concentration).

## 2. Search for Self-Stabilizing Attractors
- Concrete evidence has been produced that certain regimes (within the safe band) act as attractors where the effective 3D geometry and local causal structure (`c`) are self-stabilizing.
- This includes demonstration that, after the system settles, small internal fluctuations do not cause runaway growth or collapse of the effective dimensionality or isotropy.

## 3. Perturbation Experiments + Invariance
- Small perturbations have been applied to evolved graphs while staying inside the safe band parameters (`|λ| ≤ 0.005`, `|μ| ≤ 0.001`).
- The system has been shown to return toward the same attractor (recovery of local `d_eff`, isotropy, and protected density statistics) after perturbation.
- Quantitative measures of stability (e.g., relaxation time, deviation bounds) have been recorded.

## 4. Lean Certification of Invariance Properties
- At least one non-trivial invariance or stability property has been machine-checked in Lean on real exported data from the Python runs.
- Examples of acceptable properties:
  - Under the living candidate + safe band, the effective dimensionality (or a proxy) remains within a stated interval after small perturbations.
  - The protected density concentration is non-decreasing or bounded below in the attractor regime on real data.
  - Some monotonicity or convergence property of the coupled density–chirality dynamics holds on exported trajectories.

## 5. Python ↔ Lean Alternation
- At least **two full alternation cycles** (Python dynamics/perturbation experiments + export → Lean proof of invariance/stability property → feedback) have been completed and logged.

## 6. Logging of Changes to Previous Phases
- Any modifications to code or data from Phases 1–3 must be logged in `PHASE4_CHANGES_TO_PREVIOUS_PHASES.md` with clear descriptions.

---

**Notes for the agent**
- The core scientific question of Phase 4 is whether the living candidate, through the interplay of density sourcing, chiral protection, and quadratic self-interaction inside the safe band, creates self-stabilizing regimes for effective 3D geometry.
- Reuse the existing relational graph machinery and living-candidate implementation as much as possible.
- All work must remain strictly background-free.
- The living candidate must be used with the exact locked parameters (no ad-hoc retuning unless explicitly justified and logged).

This file is the authoritative definition of when Phase 4 is complete. Do not declare completion until every condition above is met with evidence.