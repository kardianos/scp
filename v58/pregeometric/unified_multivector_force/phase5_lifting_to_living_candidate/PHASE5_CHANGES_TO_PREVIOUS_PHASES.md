# Phase 5 Changes to Previous Phases Log

**Policy**: Any read-and-extend that becomes a modification (edit) to Phase 1–4 code or data files (or shared lean/ Phase*.lean files from prior) must be recorded here with date, description, rationale, and exact diff summary. Pure imports, new files in the Phase 5 directory, and new Lean certification modules (e.g. Phase5Lifting.lean) do **not** require an entry. The living candidate must remain byte-for-byte identical in its full ("unablated") case at all times. Ablation support is added as optional parameters with defaults that preserve 100% prior behavior and all previous numerical results.

## Rationale for minimal extension in Phase 5
Phase 5 requires *controlled clean ablations* of the quadratic self-interaction (λ Ω² + μ ⟨Ω,Ω⟩) and the ambient f_g(ρ) modulation while keeping the rest of the exact living candidate, graph evolution, protected rules, and all other machinery identical. 
- Quadratic ablation is already supported without edit (callers can pass λ=0, μ=0 to the existing update_all_omegas / compute_local_omega APIs).
- Ambient f_g ablation requires a small backward-compatible extension (optional `fg_override: Optional[float] = None` parameter) so that Phase 5 code can request a constant modulation (no ρ-dependent feedback) on otherwise identical code paths. This is the minimal change that allows rigorous "with vs without" comparison on the *same* implementation rather than a risky duplicate.

All prior-phase call sites, default behavior, and exported artifacts remain completely unchanged (verified by re-running prior demos if needed).

## Cycle 1 (2026-05-19)
- **Minimal API extension to Phase 1 `minimal_graph_model.py`** (justified and logged here):
  - Added optional parameter `fg_override: Optional[float] = None` to `compute_local_omega(...)` and forwarded through `update_all_omegas(...)`.
  - When `fg_override is None` (default): exact prior behavior — `fg = f_g_winning(rho_amb)`, `fem = 0.4 * fg`.
  - When `fg_override` provided (e.g. 0.75): use that constant value for fg/fem, removing density-dependent feedback while everything else (J_ρ/J_χ proxies, D mixing, protected suppression, quadratic iteration if λ/μ nonzero, growth bias, etc.) stays identical.
  - Updated the function docstring to document the Phase 5 use case for necessity/sufficiency experiments.
  - No other changes; no behavior change for any default call; all constants, f_g_winning, quad iteration, protected logic untouched.
  - Rationale: enables the core scientific deliverable of Phase 5 (demonstrating necessity of the two features) with maximal code reuse and fidelity. The alternative (local reimplementation of compute_local_omega in Phase 5) would duplicate ~100 lines and introduce divergence risk, violating "exact" and "clean ablation" requirements.
  - Diff summary (see search_replace performed): two function signatures + one internal conditional + docstring update. ~15 lines touched, all guarded by the new default=None path.
- No other files in phase1/, phase2/, phase3/, phase4/ were modified.
- New code (lifting_ablation.py, exports, Lean Phase5Lifting.lean, this log, PHASE5_LOG.md) created only inside phase5_lifting_to_living_candidate/ or as new shared Lean cert module.
- All full-living-candidate runs in Phase 5 use the default paths (fg_override=None, λ/μ=locked nonzero) → identical to Phase 4 data.

## Subsequent Cycles
- Entries added only on actual modifications (expected to be zero after the initial minimal extension).

**Conclusion**: The extension is the smallest possible that satisfies the Phase 5 requirement for rigorous, reproducible, side-by-side ablation experiments under a single canonical implementation of the living candidate. Full case is bit-for-bit the prior living candidate. All previous phases' results, exports, and Lean theorems remain valid and reproducible.

*This file is maintained for full auditability. The living candidate discipline is never compromised.*
