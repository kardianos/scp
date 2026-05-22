# Phase 2 — Quantitative Emergence Maps

This folder contains the work for **Phase 2** of the pre-geometric emergence program (as defined in `EMERGENCE_OF_EFFECTIVE_3D_GEOMETRY.md`).

## Primary Reference Documents
- `EMERGENCE_OF_EFFECTIVE_3D_GEOMETRY.md` (parent directory) — conceptual foundations, criteria, and roadmap.
- `PHASE2_COMPLETION_CRITERIA.md` (this folder) — authoritative definition of when Phase 2 is complete.

## Phase 1 Artifacts (for reuse)
Phase 1 artifacts live in the sibling folder:
`../phase1_minimal_relational_models/`

Key reusable pieces:
- `minimal_graph_model.py` — background-free relational graph evolution under the living candidate + Ω-feedback.
- `phase1_snapshot_cycle1.json` (and `phase1_snapshot_example.json`) — rich, Lean-ready snapshots from 500+ node runs.
- `lean/UnifiedMultivector/Phase1Relational.lean` — certified extractors (`causalPastBall` family) and monotonicity theorems on real data.

**Rule**: You may read and import from Phase 1. If you need to modify anything in Phase 1, you **must** log every change in `PHASE2_CHANGES_TO_PHASE1.md` in this folder.

## Current Status (as of hand-off from Phase 1)
- Phase 1 is complete per its own criteria.
- We have working relational graphs (≥500 nodes), living-candidate evolution with self-regulation, rich export format, and Lean-certified core extractors on real data.
- The Phase 1 agent produced a concrete starting plan for Phase 2 (see its final log).

## Instructions for Agents Working Here
- Read `EMERGENCE_OF_EFFECTIVE_3D_GEOMETRY.md` (especially Section 6 on Phase 2 and Section 3.5 on dimensionality).
- Read `PHASE2_COMPLETION_CRITERIA.md` — this is your target.
- Work in strict Python ↔ Lean alternation loops.
- Produce concrete, auditable artifacts each cycle.
- Log any Phase 1 modifications.

Start by reading the completion criteria and the conceptual document, then begin the first alternation cycle toward the Phase 2 targets.