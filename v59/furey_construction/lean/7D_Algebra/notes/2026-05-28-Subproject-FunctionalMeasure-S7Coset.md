# Sub-Project Definition: Functional Measure / Jacobian of the Full S^7 Coset

**Date**: 2026-05-28  
**Owner / Lead**: TBD  
**Status**: Proposed  
**Related Documents**:
- `2026-05-28-GapRecheck.md` and `2026-05-28-GapRecheck-Resolution.md` (the absolute scale gap)
- `v59/shulga_integration/7D_PARAMETER_DERIVATION.md`
- `v59/shulga_integration/SHULGA_DYNAMICAL_MECHANISM.md`
- `v59/ROADMAP.md` (Q3-2 and related)
- Current `ShulgaParameters.lean` (the `shulga_abs_prefactor` gap)

---

## 1. Objective

Produce a well-defined, mathematically rigorous object — the **functional measure / Jacobian factor** arising from the path integral over the fast modes on the full internal geometry (S^7 or the Spin(8)/Spin(7) coset) — such that, when combined with the already-computed Green-function contribution (`mu_raw_exact`, `lambda_raw_exact`, and their ratio), it supplies a parameter-free (or dramatically less tuned) **overall absolute scale** for the effective potential that appears in the living-candidate dynamics.

In other words: derive, or at minimum isolate and make computable, the missing overall prefactor that currently sits as the hand-tuned `1/400` in `shulga_abs_prefactor`.

---

## 2. Background & Current State

The Lean module `ShulgaParameters.lean` currently computes the *ratio* λ/μ to high precision as the quotient of two regularized Green-function sums on S^7 (Gegenbauer / harmonic sums). This ratio is fully algebraic and derived from the internal geometry.

However, the **absolute overall scale** of the effective action
\[
S_\text{eff} = -\frac12 \langle J, G J \rangle
\]
still contains an external normalization constant. In the current implementation this appears as:
```lean
def shulga_abs_prefactor : Rat := 1 / 400
```
This single rational was tuned so that `|μ|` lands near the project's observed proton-scale value (~41). All prior gap analyses have identified replacing this tuned constant with a derived quantity from the functional integral as the central remaining research-scale task for the absolute energy scale.

The broader Shulga vision documents also identify functional determinants (Jacobians) over the coset as the source of the famous 21/16 and 5 prefactors. The absolute scale for the living-candidate potential is closely related but not identical.

---

## 3. Scope

### In Scope
- Precise definition of the functional measure on S^7 (or the relevant coset) induced by the fast-mode action.
- Identification and computation (exact or symbolic) of the one-loop / Gaussian Jacobian / functional determinant around a background.
- Isolation of the overall multiplicative prefactor that multiplies the entire Green-function contribution to S_eff.
- Clear separation between:
  - The already-computed harmonic sum (ratio)
  - The measure/Jacobian contribution (absolute scale)
  - The 21/16-type ratios that arise when further integrating out gauge degrees of freedom
- At least one concrete, computable expression or algorithm that can be turned into Lean code (even if only for the leading term or volume factor).
- Testable completeness indicators (see section 6).

### Out of Scope (for the initial definition of this sub-project)
- Full non-perturbative evaluation of the path integral.
- Derivation of the nonlinear saturation parameter κ (explicitly called out as a later frontier in 7D_PARAMETER_DERIVATION.md).
- Complete derivation of the 21/16 and 5 prefactors (related but separate deliverable).
- Direct numerical matching to the living-candidate Python reports (that comes after the mathematical object exists).

---

## 4. Key Mathematical Objects (to be made precise)

1. The fast-mode action on the coset (generalization of the 1D circle case in Shulga's original work).
2. The induced Riemannian volume form / measure on S^7 or Spin(8)/Spin(7).
3. The fluctuation operator (second variation of the fast action) around a background configuration sourced by a slow multivector M.
4. The functional determinant (or regularized log-det) of that operator.
5. The precise overall prefactor that appears in front of ⟨J, G J⟩ after Gaussian integration.

The sub-project must produce a document (and eventually Lean artifacts) that defines these objects rigorously in the language of the project.

---

## 5. Suggested Phased Approach

**Phase 0 – Mathematical Specification (required first)**
- Write a self-contained mathematical note defining the objects above in the project's notation.
- Clearly state all assumptions (e.g., one-loop / Gaussian approximation, choice of regularization, treatment of zero modes, etc.).

**Phase 1 – Leading / Volume Contribution**
- Compute or isolate the pure volume factor 1/Vol(S^7) and any simple overall constants coming from the normalization of the measure.
- Produce a concrete rational (or simple expression) that can be compared with the current 1/400.

**Phase 2 – One-Loop Determinant on S^7**
- Formalize the fluctuation operator for the simplest backgrounds.
- Compute (exactly or via zeta-function / heat-kernel methods) the leading contribution to the functional determinant.
- Produce a first derived expression for the absolute prefactor.

**Phase 3 – Lifting to the Coset + Gauge Integration**
- Extend the construction to the Spin(8)/Spin(7) coset.
- Isolate the additional Jacobian from integrating out the 21-dimensional Spin(7) gauge degrees of freedom.
- Reconcile with (or derive) the 21/16 factor discussed in the vision documents.

**Phase 4 – Lean Artifacts**
- Port the results (or at least the final expressions and the separation of ratio vs. scale) into a new Lean module (e.g., `S7FunctionalMeasure.lean` or `CosetJacobian.lean`).
- Wire the derived prefactor (even if only as an option) into the living-candidate potential and re-run the existing certificates.

---

## 6. Concrete Testable Completeness Indicators

These must be checkable by inspection, by running Lean code, or by producing a specific document. Each indicator should have a clear "done / not done" criterion.

### Indicator A – Mathematical Specification (Gate 1)
- [ ] A single, self-contained document exists that defines:
  - The fast action on the coset.
  - The precise functional measure.
  - The regularization scheme used for the determinant.
  - The exact expression (or algorithm) for the overall prefactor in front of S_eff.
- The document must explicitly state what part of the absolute scale this sub-project claims to derive versus what is still left as an external constant.

### Indicator B – Separation of Ratio vs. Scale (Gate 2)
- [ ] There exists a clean mathematical decomposition (in the spec document and ideally in Lean) of the form:
  \[
  \mu_\text{derived} = C_\text{measure/Jacobian} \times G_\text{self}( \theta_c )
  \]
  where \( G_\text{self} \) is the object already computed in `mu_raw_exact`, and \( C_\text{measure/Jacobian} \) is the new object supplied by this sub-project.
- The decomposition must be stated without circular reference to the final numerical value of |μ|.

### Indicator C – First Derived Expression (Gate 3)
- [ ] At least one non-trivial, explicit expression (or algorithm) for \( C_\text{measure/Jacobian} \) has been produced that is **not** a pure fit to the proton-scale depth.
- This expression must be written in a form that can be evaluated (even approximately) and compared with 1/400.

### Indicator D – Lean Portability (Gate 4)
- [ ] A new Lean file (or extension of `ShulgaParameters.lean`) exists that at minimum can:
  - Accept the mathematical expression for the Jacobian/measure factor as a parameter or constant.
  - Recompute the living-candidate certificates using this factor instead of the hardcoded `1/400`.
- The module must build cleanly and the existing `shulga_full_*` certificates must still type-check (even if they now require a different numerical threshold).

### Indicator E – Reconciliation with Vision Documents (Gate 5)
- [ ] The sub-project document explicitly relates its result to the 21/16 and 5 prefactors discussed in `SHULGA_DYNAMICAL_MECHANISM.md` and `CROSS_AGENT_SYNTHESIS.md`.
- It must state whether the absolute scale for the living-candidate potential is the same object as, or a different component of, those prefactors.

### Indicator F – Testable Numerical Anchor (optional but strongly recommended)
- [ ] When the derived \( C_\text{measure/Jacobian} \) is inserted, the resulting |μ| lies within a documented factor (e.g., within 2× or 5×) of the currently used phenomenological value (~41), **without having used that value as an input** to the derivation.
- This is a falsifiable prediction of the functional measure, not a fit.

---

## 7. Success Criteria for the Sub-Project (Overall)

The sub-project is considered complete when **Indicators A–E** are all satisfied, and the resulting object is referenced from the main living-candidate pipeline (even if only as an alternative path with a clear comment "derived from functional measure on S^7 coset").

A stronger success would also satisfy Indicator F with a reasonably tight numerical window.

---

## 8. Risks and Open Questions

- How much of the functional determinant can actually be computed exactly versus requiring approximation or zeta regularization that is hard to formalize in Lean?
- Is the absolute scale for the living-candidate potential best thought of as coming from the S^7 volume + one-loop det on S^7, or does it require the full coset + gauge integration from the beginning?
- What is the correct treatment of zero modes on the coset?
- How does this construction interact with (or subsume) the existing Green-function computation in `ShulgaParameters.lean`?

These should be addressed explicitly in the Phase 0 specification document.

---

## 9. Suggested First Deliverable (Minimum Viable)

A short but rigorous mathematical note (5–15 pages) titled something like:

> "Functional Measure and Jacobian Factor on S^7 / Spin(8)/Spin(7) for the Shulga Effective Potential"

This note must contain at least:
- A precise definition of the objects in section 4.
- A first candidate expression for the overall prefactor.
- Explicit statements of which of the completeness indicators (A–F) it claims to satisfy.

Once that note exists and has been reviewed, the sub-project can move to Lean formalization.

---

This sub-project definition is intentionally scoped to be research-scale but with clear, checkable exit criteria. It directly attacks the last major external constant remaining in the algebraic living-candidate potential after the ratio has been derived.