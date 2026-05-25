# Integration Plan: 7D Algebra into Stability Bounds, Fock-Space Forcing, and Option D

**Date**: 2026-05-24
**Status**: ⚠ LARGELY RETRACTED (2026-05-24) — see banner below.
**Owner**: Grok (executing now)

> **RETRACTION (2026-05-24).** The "Stability Bounds" half of this plan (Steps 1–4: wiring the
> octonion algebra into the Hessian/leakage calculation and *reproducing the Python
> f-amplitude crossover*) is **retracted**. The f-amplitude crossover is a **Cl(3,0)
> multivector-model** phenomenon (the Python sweeps use a Cl(3,0) GA library, not octonions);
> the octonion Lean only *appeared* to reproduce it because of an `octMultTable` sign bug
> (`e₁·e₆`, `e₁·e₇`). With the table corrected, the octonion algebra shows **no** crossover,
> and the certificate layer was deleted (`StabilityFromAlgebra.lean`, `PhaseC_Certificates.lean`).
> See `notes/2026-05-24-LivingCandidateCrossover-ReEvaluation.md`.
>
> **Still valid — the actual deliverable of this folder:** the structural Z₂×Z₂ **grade
> forcing** on genuine Cl(7) (`PhaseB_Theorems`, `Z2Z2Forcing`, `CliffordBladeGrade`) and the
> bug-immune **Shulga ratio** (`ShulgaParameters.geometric_ratio`). Live forward targets: the
> Fock/G₂-branching extension of the forcing, and the `21/16`,`5` coset derivation. Disregard
> the crossover-reproduction goal (Steps 1, 3, 4, and the crossover part of Step 2).

## Goal
Fully integrate the concrete 7D algebra realization (`SevenDAlgebra.lean` + supporting modules) into the existing stability bounds formalization, the Furey Fock-space Z₂×Z₂ derivation, and the Option D composite u-quark work.

The 7D Algebra agent has delivered the missing explicit matrices and Fock labeling. The previous 4-phase Lean effort and related high-effort agents now have the engine they need. This plan turns that foundation into connected, derivable results instead of parallel scaffolding.

## Four Concrete Steps (Execute in Order)

### Step 1: Immediate Bridge (Highest Leverage)
- Wire the explicit `sectorDiags`, `gamma` generators, `L_bivector_*`, and `F_fourform_*` from `SevenDAlgebra.lean` into the `density_algebra/lean/` modules.
- Update or extend `StabilityFromAlgebra.lean` (or create a clean bridge module) so that the Hessian and leakage calculations use the real matrices instead of schematic injection.
- Make the `fAmplitudeCrossoverDemo` and related stability checks compute from the actual algebra.

**Success criteria for Step 1**:
- `lake build` succeeds across both `furey/lean/` and `density_algebra/lean/`.
- At least one key claim (e.g., L-grade products leak into F-grade 4-forms) is now computed directly from the table rather than asserted.
- The Python `f_amplitude` crossover numbers can be reproduced (or bounded) using the Lean matrices.

### Step 2: Strengthen the Forcing Theorems
- Update the sketched theorem in `furey_construction/lean/Predictions.lean` (from Option D v2 work) to reference the concrete matrices on the |Ω_N⟩ blocks.
- Create or strengthen theorems in the stability work showing that only the observed L/F projectors yield positive-definite Hessians for the living-candidate potential when using the real multiplication table.
- Add at least one machine-checkable certificate or `rfl` / decidable example that ties a specific numerical regime (λ=0.005, μ=0.001, f_amplitude > 0.4) to the matrix results.

**Success criteria for Step 2**:
- The Z₂×Z₂ assignment and the stability bounds have at least partial theorems that depend on the explicit 7D matrices (not just axioms or simplified models).
- Clear traceability from the Fock labeling + generator action → the observed discrete structures.

### Step 3: Cross-Agent Synthesis Document
- Write a single, high-quality synthesis note that ties together the three recent high-effort threads:
  1. Furey Fock-space derivation of the Z₂×Z₂ assignment (`13_fock_mass_forcing_report.md` and updates)
  2. Option D structural/representation hypothesis (`FINDINGS_option_d_composite_v2.md`)
  3. 7D Algebra + 4-phase stability bounds formalization (this folder + `density_algebra/lean/` work)
- The document should show how the explicit matrices close the loop between representation theory, protection/density stacking, and machine-checkable stability conditions.

**Success criteria for Step 3**:
- A clear, citable document (probably in `v59/synthesis/`) that a future reader or agent can use to understand the current state of the core structural claims without re-reading five separate reports.

### Step 4: Next Technical Milestone — Real Hessian + Certificates Tied to Data
- Extend the current generators to a more complete set of L-grade and F-grade operators.
- Compute (or finite-difference) the actual Hessian of the quartic living-candidate potential restricted to different projector images, using the real matrices.
- Turn the Gershgorin / positive-definiteness certificates into something that can ingest (or be validated against) the concrete numbers coming from the Python `f_amplitude` sweeps.
- Produce at least one example where a specific Python run result is certified (or bounded) inside Lean.

**Success criteria for Step 4**:
- The stability bounds are no longer "modeled to match Python" but are computed from the algebra and then compared to Python.
- Clear path for future agents to replace finite differences with symbolic differentiation once Lean tooling improves.

## Execution Rules for This Plan
- Work step by step. Do not skip ahead.
- After each major step (or sub-step), verify with `lake build` and relevant `#eval` or Python runs.
- Write a short dated note in `notes/` after each step documenting what was done, what worked, and what the next concrete action is.
- Before declaring any step "done", explicitly check the work against the success criteria above.
- If a step is incomplete, continue iterating on it before moving to the next.
- Prioritize **derivation from the explicit matrices** over further modeling.

## Verification Cycle (Mandatory at the End)
Before considering the overall plan complete:
1. Re-read this `INTEGRATION_PLAN.md`.
2. Check every success criterion for Steps 1–4.
3. Confirm that the key claims (leakage, crossover, Z₂×Z₂ forcing, Option D compositeness) are now grounded in the real 7D multiplication table and explicit action on |Ω_N⟩ rather than schematic injection.
4. If any criterion is not met, continue working on the relevant step(s).

This plan turns the excellent foundation delivered by the 7D Algebra agent into connected, derivable results across the Furey, Option D, and stability-bounds threads.