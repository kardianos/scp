# Phase 1 Completion Event — Minimal Relational Models

**Phase 1 is considered complete only when ALL of the following conditions are simultaneously satisfied:**

## 1. Scale & Background-Free Evolution
- At least one purely relational (background-free) graph with **≥ 500 nodes** has been successfully evolved under the **exact locked living candidate**:
  ```
  ⟨D Ω + λ Ω² + μ ⟨Ω, Ω⟩⟩_{0,2} = f_g(ρ_ambient) · J_ρ + f_em(ρ_ambient) · J_χ
  ```
  using the previously determined winning parameters (`f_g(ρ) = 1/(1 + ρ/ρ_crit)`, `ρ_crit ≈ 2.5`, `|λ| ≤ 0.005`, `|μ| ≤ 0.001`).

## 2. Python ↔ Lean Interoperability
- The Python simulator can produce serializable snapshots (JSON) that are directly ingestible by a Lean model in this folder.
- At least one non-trivial extraction algorithm or statistic (e.g., `causal_past_ball`, growth curve, or a simple d_eff proxy) has a corresponding Lean implementation.

## 3. Formal Certification (Lean)
- At least **one non-trivial property** of a Phase-1 statistic or extractor has been **machine-checked in Lean** on real exported snapshots from the Python simulator.
  Examples of acceptable certified properties:
  - Causal ball size is strictly monotonic with depth.
  - Growth curve satisfies a stated bound or inequality on the exported data.
  - Some algebraic identity used by an extractor holds on real snapshots.

## 4. Demonstrated Alternation
- The agent has completed **at least two full alternation cycles** of the form:
  Python simulation / evolution → export snapshot(s) → Lean work (ingestion + proof or certification) → insight/feedback back into the Python model or next simulation parameters.

---

**Notes for the agent:**
- The goal is not maximum scale, but a minimal yet meaningful demonstration that a background-free relational model under the living candidate can be evolved, measured, and partially certified in Lean.
- You may propose tighter or more insightful completion criteria later, but you must first reach (or exceed) the bar defined above before declaring Phase 1 complete.
- Keep all work strictly relational — no lattices or coordinate systems may be re-introduced.

This file is the authoritative definition of "Phase 1 done" for this sub-folder. Do not declare completion until every bullet is satisfied and documented.