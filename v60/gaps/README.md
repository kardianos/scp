# v60 Gaps — Continuation and New Attack

This directory holds the updated gap attack for v60, focused on the two highest-priority items from the v59 closeout:

- G9 (induced-metric / emergent tensor gravity)
- G1 (rank tension / two-piece Y)

Subdirectories will be created per gap (e.g., `rank_tension/`, `induced_metric/`) or the work may live primarily under `../gravity_recast/` and `../synthesis/` with Lean modules under `../lean/`.

All work here follows the same standards as `v59/gaps/`:
- `README.md` (this file) + `ALTERNATIVES.md` + `FINDINGS.md` per sub-attack.
- Lean module with `sorry` only on the genuine open conjectures.
- Numerical scripts that run clean.
- Full documentation of both positive resolutions and clean falsifications.

See `../ROADMAP.md` and `../PLAN.md` for the current priority order and variant structure.

**Status (kickoff 2026-05-25)**: Scaffolding + first concrete attack launched.

- **G9 (top priority)**: Work has begun in `../gravity_recast/`. First artifact: `01_8space_to_spacetime_bindings.py` + `01_findings.md`. This directly addresses the 8-space (V^8 / Spin(8) / so(8) ≅ Λ²) binding problem and the required alteration (or re-interpretation via soldering) of the current fundamental OBE gravity equation `□ Ω_grav = f_g ρ_grav`. It enumerates the options from v59/gaps/gravity/ALTERNATIVES.md (pure Plebański A2 is the leading candidate), produces a viability table, and sets the immediate next steps (constraint implementation + DOF/helicity count reusing the v59 polarization machinery). See the script for the full restatement of the equation, the definition of V^8 consistent with the v59 algebra, and the binding ansätze (A1 KK-style split, A2 Plebański 2-form, B1 Sym² branching, C external fallback).

- **G1 (rank tension)**: Still scaffolding. Concrete work to begin in parallel under `rank_tension/` or `../synthesis/`.

The v59 `gaps/` tree (especially `gravity/ALTERNATIVES.md`, `FINDINGS.md`, `G8G9_Gravity.lean`, `g9_polarization_test.py`) remains the authoritative reference for all supporting theorems, helicity code, and the current equation. All v60 work must remain compatible with (or explicitly extend) the proved v59 skeleton (lepton=L forcing, Koide from G2/Spin(7), 784 End(L), EP-exact ρ_grav, etc.).

---

*This is a living directory. Update this README as the attack structure solidifies.*