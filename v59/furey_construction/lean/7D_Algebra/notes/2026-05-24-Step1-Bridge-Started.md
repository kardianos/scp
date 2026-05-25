# Step 1: Immediate Bridge — Started (2026-05-24)

**Goal**: Wire the explicit matrices and `sectorDiags` from `SevenDAlgebra.lean` into the stability bounds work so that leakage, Hessian, and f-amplitude crossover are computed from the real algebra rather than modeled.

**Current State at Start of Step**:
- `StabilityFromAlgebra.lean` already exists in this folder and imports `SevenDAlgebra`.
- It defines matrix helpers, `crossoverLeakageDemo`, a numeric `V_full` + `computeHessianNum`, and the string explanation `onlyClosedMasksStableEmergence`.
- `SevenDAlgebra.lean` has the table, Fock labeling, `gamma`, sample L/F matrices, and `sectorDiags`.

**Actions Taken**:
- Confirmed the bridge module is already using the real matrices for the leakage demo.
- The numeric Hessian is operating on the octonion coefficients via the real `octMult`.

**Next immediate sub-actions for this step**:
1. Make `StabilityFromAlgebra.lean` the canonical exported bridge (add clear public API section).
2. Update `density_algebra/lean/StabilityBounds.lean` and `OctonionAlgebra.lean` with strong "source of truth" comments + thin re-exports or imports where possible.
3. Write a short integration note in `density_algebra/notes/`.
4. Verify `lake build` still works from both sides.
5. Mark Step 1 complete only after explicit verification against the success criteria in `INTEGRATION_PLAN.md`.

This note will be updated as sub-actions complete.