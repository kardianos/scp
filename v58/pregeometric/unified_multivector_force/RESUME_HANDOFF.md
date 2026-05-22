# Resume / Handoff Guide — Unified Multivector Force Law Experiment

**Location**: `v58/pregeometric/unified_multivector_force/`  
**Date of this handoff**: 2026-05-19

## 1. Purpose of This Document

This file is written so that a future session (or new context) can resume the experiment with minimal re-reading. It points to everything important and states the exact next concrete tasks.

## 2. Key Documents (Read in This Order)

1. `EXPERIMENT_OUTLINE.md` — overall philosophy, dual-track structure, coordination protocol, success criteria.
2. `TERMINATION_CRITERIA_AND_CURRENT_STATUS.md` — precise definition of done + current concrete gap.
3. `BACKGROUND_AND_SPECULATIVE_EQUATIONS.md` — the current living set of candidate unified equations (A–D and any later additions).
4. `COORDINATION_LOG.md` — the running history of every round and the guidance issued between Python and Lean.
5. `PYTHON_FINDINGS.md` and `LEAN_FINDINGS.md` — the detailed outputs of the two tracks (most recent entries first).

## 3. Current State Summary (Concise)

- **Living candidate** (as of latest cycles): the projected form
  ```
  < D Ω + λ Ω² + μ ⟨Ω, Ω⟩ >_{0,2} = f_g(ρ_amb) · J_ρ(grade 1) + f_em(ρ_amb) · J_χ(grade 2)
  ```
  with the winning `f_g(ρ) = 1 / (1 + ρ_ambient / ρ_crit)` and conservative safe band |λ| ≤ 0.005, |μ| ≤ 0.001.

- **Python track**: Has ultra-dense (up to 500×500+) 2D retarded validation, thousands of real exported snapshots (ganja JSONs on disk), concrete numeric metrics, A-vs-B comparison, and protected-chirality origin observations with algebraic justification.

- **Lean track**: Has a working concrete Fin-8 model that ingests real exported snapshots, multiple machine-checked geometric-product identities on real exported data (including the density quadratic), and retarded implication theorems that now contain real (non-`sorry`) proof text. The last schematic piece is the concrete retarded operator realization on the largest snapshots.

- **Gap**: We do not yet have the first complete (fully non-schematic) machine-checked retarded implication theorems on real ultra-dense exported data.

## 4. Exact Next Concrete Tasks (to close the gap)

**Priority order**:

1. **Declare the current best form the official living candidate** (update `BACKGROUND_AND_SPECULATIVE_EQUATIONS.md` and `COORDINATION_LOG.md`).

2. **Python next actions** (choose one or more):
   - Generate the next richer snapshot batch (e.g., 200×200 or 300×300 with the full set of protected variants) and write the ganja JSON files to disk.
   - Produce full 40-step A-vs-B + protected trajectory tables for the new batch.
   - Export a few complete field snapshots in a format easy for the Lean model to ingest (full 8-component at multiple retarded times).

3. **Lean next actions**:
   - Ingest the newest snapshot batch into the concrete model.
   - Prove at least one more geometric-product identity on the new real snapshots.
   - Replace the remaining schematic piece (retarded operator realization) in at least one implication theorem (`candidateA_implies_newtonian_limit_retarded` or the B variant) and produce a complete (non-schematic) proof.
   - Update the trackers and post the new locked items.

4. **Coordination**:
   - After each agent cycle, append to the relevant findings file and to `COORDINATION_LOG.md`.
   - The coordinator (or the other agent) should continue the synthesis + guidance pattern already established.

## 5. How to Resume the Agents

- Python agent last subagent_id (as of this handoff): see the most recent `PYTHON_FINDINGS.md` entry.
- Lean agent last subagent_id: see the most recent `LEAN_FINDINGS.md` entry.

When resuming, give each agent the same instruction used in prior cycles:
- Read `EXPERIMENT_OUTLINE.md`, `BACKGROUND_AND_SPECULATIVE_EQUATIONS.md`, and the latest `COORDINATION_LOG.md`.
- Continue the tight feedback loop.
- Write a new structured entry in their findings file after the cycle.

## 6. Useful Background References (Outside This Folder)

- `v58/pregeometric/MULTIVECTOR_FORCE_LAW.md`
- `v58/pregeometric/PARTICLES_AS_DENSITY_ACHIEVERS.md`
- `v58/pregeometric/MEDIUM_DYNAMICS_TAILS_AND_WAKES.md`
- `v58/first_principles/EXPECTED_BEHAVIOR.md` (especially §3.6 and the new bullet about density as goal)

---

*If you are reading this later, the experiment is in a healthy, well-documented state with a clear gap and a clear path to the first defensible result. The loop is designed to be resumable with the documents above.*