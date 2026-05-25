# Continuation Plan: 4-Phase Lean Formalization of Stability Bounds (Post 7D Algebra Foundation)

**Date**: 2026-05-24  
**Status**: Hand-off document for the next high-effort agent.  
**Primary Location**: `v59/furey_construction/lean/7D_Algebra/`  
**Related Work**: `v59/density_algebra/lean/` and `v59/density_algebra/notes/`

---

## 1. Desired End Game (The Real Target)

The ultimate goal of the 4-phase Lean effort is **not** to build a nice model that reproduces observations we already like. The goal is to reach a state in which the key physical claims of the project become **derivable or machine-verifiable consequences** of the explicit algebra plus the explicit dynamics.

**Concrete success criteria**:

- A faithful formalization of the 8-dimensional algebra (with the project’s actual Fano multiplication table) inside Lean.
- The fermion states |Ω_N⟩ (N = 0,1,2,3) defined as explicit vectors in the 8-dimensional spinor representation, with correct labeling by form degree and color weight (lepton singlets on N=0/3, d-quarks on the 3, u-quarks on the 3̄).
- The living-candidate effective potential (or at minimum its second variation / Hessian) expressed in terms of the actual algebra: ρ_M + λ·scalar(M·M) + μ·norm term, plus ambient modulation.
- Protection modeled rigorously as projectors onto subsets of the 8 components (L, F, L+F, etc.).
- Theorems or strong machine-checkable certificates of the following form:

  > “For λ, μ inside the physically allowed small band, the quadratic form on the image of a projector P is positive-definite if and only if P corresponds to a closed subalgebra under the multiplication table, and the background amplitude lies at one of the discrete solutions (|ξ|² = 1/2, 3/5, or 7/9).”

- The observed discrete data (D = 28/35/63, the specific |ξ|² values, the f-amplitude crossover, the 21-factor, the Z₂×Z₂ bit assignment, etc.) should appear as **outputs** or **necessary conditions** rather than inputs.
- The Python numerical results (especially the protection differentiation and f_amplitude crossover) and the Maxima Hessian spectra should be recoverable as instances or consequences inside the same formal development.

In short: Lean should be able to **explain** why only certain discrete algebraic choices survive, using the actual multiplication table and the actual action on the fermion states.

---

## 2. Current State of the Work (What Has Already Been Delivered)

### 2.1 Work from the 4-Phase Lean Agent (019e575a-fb34-7572-b686-c905a424e870)

**Strengths**:
- Excellent process and documentation: 14+ high-quality dated notes in `v59/density_algebra/notes/` documenting 7+ angles per phase.
- Good scaffolding in `StabilityBounds.lean` and `DensityForcing.lean` (`StabilityCert`, `ProtTech`, `fAmplitudeCrossoverDemo`, `RadialZeroCondition`, etc.).
- Made simplified eigenvalue checks evaluable (`#eval` works on the diagonal model).
- Created `OctonionAlgebra.lean` that imports the Fano table and sketches the desired anisotropic story.
- `lake build` succeeds cleanly.

**Limitations / Gaps**:
- Core algebraic claims remain largely **modeling by construction** rather than derivation.
- `octMult` is a stub.
- `buildHessian8` starts from the simplified diagonal and manually injects negatives for non-closed masks.
- The f-amplitude crossover and “only closed subalgs work” statements are asserted via axioms or schematic injection, not computed from the multiplication table.
- No explicit 8×8 matrices for generator action on the labeled |Ω_N⟩ states (the exact piece repeatedly called out as missing by the Furey Fock-space agent and the 7D Algebra agent).

### 2.2 Work from the 7D Algebra Agent (019e5763-65f0-7443-849f-2b182960218a) — The Critical New Foundation

This agent created the dedicated folder you are now in and delivered the missing foundation:

- `SevenDAlgebra.lean` — faithful Fano multiplication table, correct Fock labeling of the |Ω_N⟩ states (N=0 lepton, N=1 d-quark, N=2 u-quark, N=3 lepton), `gamma(k)` generators, and the first explicit computable matrices for L-bivector and F-4-form action on the 8 states.
- `PLAN.md` and `ROADMAP.md` (living documents).
- Excellent notes documenting the work.

This is the exact missing piece. The previous 4-phase effort now has the concrete algebraic engine it was missing.

---

## 3. Rearticulated End Goal (After Seeing Both Efforts)

The target is now:

**A Lean development in which the second variation of the living-candidate effective potential can be computed (or rigorously bounded) from the multiplication table + the explicit form of the potential, for arbitrary protection projectors. From this one can prove (or obtain strong machine-checkable certificates for) the key forcing statements, and the Python/Maxima observations become recoverable consequences rather than inputs.**

The 7D Algebra work has delivered the necessary foundation (table + Fock labeling + generator matrices). The remaining work is to **bridge** from this foundation into the stability bounds and selection theorems.

---

## 4. Precise Gap Analysis

The gap is no longer “we don’t have the algebra.” The gap is now:

**Integration + Derivation**

1. The 7D Algebra module exists but is not yet wired into the stability bounds machinery (`StabilityBounds.lean` and `OctonionAlgebra.lean` in `density_algebra/lean/` still use a simplified/schematic model).
2. The explicit generator matrices (`gamma`, `L_bivector_01`, `F_fourform_0123`, etc.) exist in `SevenDAlgebra.lean`, but have not yet been used to compute or verify the Hessian of the living-candidate potential.
3. The “only closed subalgebras give stable wells” and “f-amplitude crossover” claims are still modeled rather than derived from the actual multiplication table + dynamics.
4. There is no yet a clean bridge (import, shared definitions, or certificate format) between the new 7D_Algebra work and the existing stability-bounds scaffolding.

---

## 5. Recommended Plan for the Next Agent (Prioritized)

**Primary Objective**: Use the new concrete 7D Algebra foundation to move the stability bounds from “modeled” to “derived/computed.”

### Phase A — Bridge the Two Bodies of Work (Highest Priority)

1. Create a clean integration module (suggested name: `StabilityFromAlgebra.lean` or `HessianFrom7D.lean`) inside `7D_Algebra/` or a new shared location.
2. Import or mirror the Fock labeling, `gamma` matrices, and `octMult` from `SevenDAlgebra.lean`.
3. Re-express (or compute) the living-candidate second variation using the actual multiplication table.
4. Compute (or symbolically sketch) the Hessian for the key projectors (L, F, LF) on mixed backgrounds.
5. Verify that the `fAmplitudeCrossoverDemo` behavior now emerges from the algebra rather than being injected.

### Phase B — Strengthen the Core Claims

- Turn `onlyClosedMasksStable` and the f-amplitude crossover into theorems or strong certified statements that depend on the actual table (not schematic injection).
- Add the first non-trivial proofs or machine-checkable certificates (e.g., “for λ=0.005, μ=0.001 and f_amplitude > 0.4, the L-only projector on a mixed source has a negative eigenvalue while LF does not”).
- Begin formalizing the link to the Z₂×Z₂ assignment using the explicit generator matrices on the |Ω_N⟩ states (tie directly to the work in `13_fock_mass_forcing_report.md`).

### Phase C — Computability and Certificates

- Make it possible to compute critical values inside Lean (min μ/λ for a mask, critical f_amplitude, etc.).
- Define a clean `StabilityCertificate` type that the Python sweeps can emit and Lean can verify.
- Write a small emitter script (Python) that can turn future sweep JSONs into Lean `example` or `theorem` statements.

### Phase D — Polish and Documentation

- Update `PLAN.md` and `ROADMAP.md` in this folder with the new status.
- Write a short integration note in `v59/density_algebra/notes/` explaining how the 7D_Algebra work now enables the next phase of the bounds formalization.
- Ensure everything builds cleanly from the parent lakefile when possible.

---

## 6. Suggested Working Style for the Next Agent

- Start by reading this file + `PLAN.md` + `ROADMAP.md` + the two most recent notes in this folder + the key notes from `v59/density_algebra/notes/` (especially the Lean roadmap completion summary).
- Use `todo_write` to track the four sub-phases above.
- After every significant piece of work (new module, first real Hessian from the table, first verified crossover, etc.), write a dated note.
- Prioritize **derivation over modeling**. Prefer computing the Hessian from the table (even if numerically or for specific generators) over asserting the desired behavior.
- Keep the code buildable and the documentation excellent — this work will be handed to future agents.

---

## 7. Immediate First Actions (Recommended)

1. Read the following files in this order:
   - This file (`CONTINUATION_4PHASE_LEAN_BOUNDS.md`)
   - `PLAN.md` and `ROADMAP.md` (this folder)
   - `SevenDAlgebra.lean` (focus on the gamma matrices, Fock labeling, and `printNonzero` examples)
   - `v59/density_algebra/lean/OctonionAlgebra.lean` and `StabilityBounds.lean`
   - The most recent notes in both `7D_Algebra/notes/` and `density_algebra/notes/`

2. Decide whether to keep all work inside `7D_Algebra/` or create a thin integration layer that both `7D_Algebra/` and `density_algebra/lean/` can depend on.

3. Begin Phase A (bridging) — this is the highest-leverage next step.

---

**This file is the hand-off.**  
Point the next agent at this document (and the folder) and tell it: “Read CONTINUATION_4PHASE_LEAN_BOUNDS.md first, then execute the plan.”

The previous 4-phase agent left excellent scaffolding and documentation.  
The 7D Algebra agent has now delivered the missing algebraic foundation.  
The next agent’s job is to connect the two and turn the modeled bounds into derived ones.

Good luck. This is one of the most important remaining pieces of formalization work in the project.