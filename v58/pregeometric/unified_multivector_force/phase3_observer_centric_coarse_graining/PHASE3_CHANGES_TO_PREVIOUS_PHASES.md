# Phase 3 Changes to Previous Phases Log

**Policy**: Any read-and-extend that becomes a modification (edit) to Phase 1 or Phase 2 code or data files must be recorded here with date, description, rationale, and exact diff summary. Pure imports and new files in the Phase 3 directory do **not** require an entry.

## Cycle 1 (2026-05-19)
- No modifications performed.
- `observer_model.py` used `import minimal_graph_model as pgm` (read-only execution of pure functions: create_seed_graph, add_nodes_..., compute_local_omega, causal_past_ball, RelationalNode, constants, etc.) and analogous re-implementation of small map helpers.
- All evolution, observer logic, exports, and numbers produced inside the Phase 3 folder.
- Lean module written as new file; self-contained to avoid triggering edits to Phase1Relational.lean.
- Result: 0 entries needed. Full fidelity to "exact living candidate" and prior artifacts without touching them.

## Cycle 2 (2026-05-19)
- Same as Cycle 1: no edits to `../phase1_minimal_relational_models/` or `../phase2_quantitative_emergence_maps/`.
- Minor internal refinement (weighted ruler) implemented locally in Phase 3 script.
- Exports and Lean updates remained local to Phase 3 + the shared lean/ tree (new Phase3 module only).

**Conclusion**: Throughout the autonomous Phase 3 work, the rule was followed with zero modifications to previous-phase source or data. All work reused prior logic via import or faithful local re-expression while keeping the living candidate byte-for-byte identical.

*This file is maintained for auditability even when empty.*