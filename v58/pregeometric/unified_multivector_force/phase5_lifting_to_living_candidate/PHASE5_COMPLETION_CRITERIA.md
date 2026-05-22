# Phase 5 Completion Criteria — Lifting to the Living Candidate

**Phase 5 is considered complete only when ALL of the following conditions are simultaneously satisfied:**

## 1. Controlled Ablation Experiments (Python)
- The Python relational framework has been extended to run **controlled ablation** experiments on graphs evolved under the exact living candidate.
- Specific terms have been disabled or varied while keeping everything else identical:
  - Quadratic self-interaction terms (`λ Ω² + μ ⟨Ω, Ω⟩`) turned off or set to zero.
  - Ambient modulation `f_g(ρ)` replaced by a constant (i.e., no density-dependent feedback).
- Multiple runs (on graphs with ≥ 500 nodes) have been performed with these ablations.
- Quantitative metrics show that the self-stabilizing attractor behavior observed in Phase 4 (bounded `d_eff`, isotropy, protected density concentration, recovery after perturbation) degrades significantly or disappears when these terms are removed.

## 2. Necessity and Sufficiency Demonstrated
- Clear, reproducible evidence (numbers + plots or tables) that the quadratic terms and the ambient modulation `f_g` are **necessary** for the stability properties found in Phase 4.
- Evidence that they are also **sufficient** (i.e., when the full living candidate with these terms inside the safe band is used, the attractor behavior reappears consistently).

## 3. Lean Certification on Real Exported Data
- At least one non-trivial property has been machine-checked in Lean demonstrating necessity or sufficiency on **real exported data** from the ablation experiments.
- Examples of acceptable properties:
  - Without the quadratic terms, a key invariance (e.g., bounded `d_eff` or protected density recovery) fails on the exported trajectories.
  - With the full living candidate (quadratic + `f_g`), a stability bound holds that does not hold in the ablated versions.
  - A theorem linking the presence of the quadratic + ambient terms to the existence of the attractor on concrete living-candidate-derived graphs.

## 4. Python ↔ Lean Alternation
- At least **two full alternation cycles** (Python ablation experiments + export → Lean proof on the exported ablation data → feedback into Python) have been completed and logged.

## 5. Logging of Changes to Previous Phases
- Any modifications to code or data from Phases 1–4 must be logged in `PHASE5_CHANGES_TO_PREVIOUS_PHASES.md` with clear descriptions.

---

**Notes for the agent**
- The scientific goal of Phase 5 is to close the loop by showing that the features that make the living candidate special (quadratic self-interaction + ambient density modulation) are not arbitrary but are required for the self-regulating effective geometry observed in earlier phases.
- Reuse the existing stability harness from Phase 4 (`stability_dynamics.py` or similar) as much as possible.
- All work must remain strictly background-free.
- The living candidate must be used with its exact locked parameters in the “full” case; ablations must be clean and clearly documented.
- When Phase 5 is complete, the program will have shown both sufficiency (Phases 1–4) and necessity (Phase 5) of the key terms in the living candidate for stable effective 3D geometry on relational graphs.

This file is the authoritative definition of when Phase 5 is complete. Do not declare completion until every condition above is met with evidence.