# Furey Construction — Active Workspace

**Parent**: [`../README.md`](../README.md), [`../SUMMARY.md`](../SUMMARY.md)
**Plan**: [`PLAN.md`](PLAN.md)
**Status**: Active execution of the multi-variant attack on full Standard Model construction.

This directory implements Step 12+ of the v59 program: the full $\mathbb{C} \otimes \mathbb{H} \otimes \mathbb{O}$ construction proposed by Cohl Furey, Geoffrey Dixon, Tevian Dray and others. We arrived here independently via the kernel-fit and octonionic-extension steps; the convergence is documented in [`../octonionic_extension/02_findings.md`](../octonionic_extension/02_findings.md).

## Variants

The execution plan is detailed in [`PLAN.md`](PLAN.md). In order:

- A — `01_choh_algebra.py` — Build C⊗H⊗O algebra explicitly.
- B — `02_sm_idempotent.py` — Identify SM idempotents.
- C — `03_alpha_from_u1.py` — Derive U(1)_em coupling.
- D — `04_alpha_prediction.py` — Test the π²/2 conjecture.
- E — `05_quark_sector.py` — Quark mass spectrum.
- F — Lean formal verification (`lean/`).
- G — `06_gravity_sector.py` — Gravity sector identification.
- H — `ALL_FINDINGS.md` — consolidation.

Each variant produces a `NN_findings.md` document with its result, whether positive, partial, or null.
