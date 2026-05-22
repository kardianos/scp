# Phase 5 — Lifting to the Living Candidate

This folder contains the work for **Phase 5** of the pre-geometric emergence program (as defined in `EMERGENCE_OF_EFFECTIVE_3D_GEOMETRY.md`, Section 6).

## Primary References
- `EMERGENCE_OF_EFFECTIVE_3D_GEOMETRY.md` (two levels up) — especially Section 6 (Phase 5 description) and the overall program goals.
- `PHASE5_COMPLETION_CRITERIA.md` (this folder) — authoritative target.

## Previous Phases
- **Phase 1**: `../phase1_minimal_relational_models/`
- **Phase 2**: `../phase2_quantitative_emergence_maps/`
- **Phase 3**: `../phase3_observer_centric_coarse_graining/`
- **Phase 4**: `../phase4_stability_and_self_regulation/`

You may read and extend artifacts from all previous phases. Any changes must be logged in `PHASE5_CHANGES_TO_PREVIOUS_PHASES.md`.

## Phase 5 Focus
Return to the full living candidate equation and demonstrate (through controlled ablation in Python + Lean certification on real data) that the quadratic self-interaction terms (`λ Ω² + μ ⟨Ω, Ω⟩`) and the ambient modulation `f_g(ρ)` are **necessary and sufficient** for the self-stabilizing attractor behavior observed in Phase 4.

This is the final “necessity” step that closes the loop with the original 2D results and the living candidate.

The agent must:
- Maintain strict Python ↔ Lean alternation.
- Keep all work background-free and use the exact living candidate (with clean ablations for the “without” cases).
- Log any modifications to prior phases.
- Continue autonomously until every condition in `PHASE5_COMPLETION_CRITERIA.md` is met.

## Starting Point
Begin by reading the completion criteria and the conceptual document. Review the Phase 4 stability harness (`stability_dynamics.py`) as the natural foundation for ablation experiments. The first practical task is typically to add clean ablation modes (disable quadratic terms, disable `f_g`, etc.) and run the same perturbation-recovery protocol from Phase 4 on the ablated versions for direct comparison.

This phase is intended to provide rigorous evidence that the specific features of the living candidate are not arbitrary but are required for the emergence and stability of effective 3D geometry in this pre-geometric setting.