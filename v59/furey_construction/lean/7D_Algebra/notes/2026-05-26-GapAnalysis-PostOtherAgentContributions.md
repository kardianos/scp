# Gap Analysis: Recent Work (Other Agent Items 2/3/4 + Prior Living-Candidate Pipeline) — 2026-05-26

**Date**: 2026-05-26
**Status**: Honest post-exploration assessment following non-lazy verification discipline
**Explorer**: Grok (this session)
**Context**: User request "Please explore the recent work, and see where any further gaps are." after other agent delivered on critiques 1 (Rat), 2 (ShulgaParameters.lean), 3 (strengthened PhaseB_Theorems.lean), 4 (CROSS_AGENT_SYNTHESIS.md refresh). This note continues the "continue... cycle until done" and "check your work against criteria... if it isn't complete, then continue" protocol established across the INTEGRATION_PLAN execution.

## Re-Read of Governing Documents (Mandatory First Step)
- Re-read `INTEGRATION_PLAN.md` (full, lines 1-68). Success criteria for Step 4 (the current active milestone):
  - Stability bounds no longer "modeled to match Python" but *computed from the algebra* and then compared.
  - Gershgorin/pos-def certificates ingest or validate against concrete Python sweep numbers.
  - Clear path for symbolic differentiation.
  - (Implicit from execution history + notes): Full living-candidate V (retarded, exact Python structure), complete L=28/F=35 generator basis in the Hessian path, living thresholds from raw data (not model-tuned), general pure-Lean JSON loader, simulateCrossover updated to rich V.
- Re-read `CROSS_AGENT_SYNTHESIS.md` (full), `ShulgaParameters.lean` (full), `PhaseB_Theorems.lean` (full), `StabilityFromAlgebra.lean` (full via chunks + grep), `PhaseC_Certificates.lean` (key sections + end), multiple 2026-05-25 notes (e.g., PermanentRealDataTable, LivingCertsCalibrated, FampEndToEndPipeline, MatureReportCertifiedComparison, RealRhoMInV, ProjectNonlinearV, LivingCandidateV).
- Re-read `SHULGA_DYNAMICAL_MECHANISM.md` and `7D_PARAMETER_DERIVATION.md` (shulga_integration/) for the aspirational scope vs delivered Lean slice.
- Re-read `FINAL_REPORT.md`, `SevenDAlgebra.lean` (key excerpts), lakefile.lean, the two JSON slices used by loaders.
- Ran `lake build` (clean success, 0 jobs, cached but no elaboration errors; #guard in Shulga passes; modules with #eval elaborate cleanly).
- Grepped for wiring ("Shulga|geometric_ratio|V_living_project|kappa|phenomenological"), forcing patterns, and cross-references.
- Verified file paths against v59-now-in-root (all under `/home/d/code/scp/v59/...`).

No AGENTS.md/CLAUDE.md in v59/ subtree (root CLAUDE.md governs; respected: no kernel edits, SFA not touched, this is pure formalization).

## Exploration Performed (Non-Lazy, Actual Reads + Verification)
- **ShulgaParameters.lean** (SCPv59.Furey7D.Shulga): 60 lines. Exact Rat recursive `greenSumLoop` + `computeGreenSum` (Gegenbauer C_l^{(3)} style on S^7 harmonics) for two cos values (Z3 shift -1/2 and cutoff 127/128), lmax=100 hardcoded. `lambda_raw_exact`, `mu_raw_exact`, `geometric_ratio = lambda/mu`, `#guard` loose interval (-1/4000, -1/6000). Standalone module. Imported in lakefile/eval_test/format_lean.py but **zero references** in StabilityFromAlgebra, PhaseC, or V_living_project. Finite truncation, not closed-form. Raw ratio ≈ -0.000205 (directionally matches target ~ -1/5000).
- **PhaseB_Theorems.lean**: 82 lines, 5 theorems. All `rfl` on pre-computed `Bool` defs (`check_all_L_diags_zero`, `check_dquark_requires_F_grade`, `check_uquark_requires_LF_grade`, etc.) that enumerate `all_L_bivectors` / `all_F_fourforms` (from SevenDAlgebra) and test `diag`/`sectorDiags` properties. Stronger than pre-other-agent (specific "L fails diagonals for color triplets; F succeeds"; "L strictly off-diagonal"; lepton compatibility). Still exhaustive sample checks, not structural identities derived from Fano multiplication table axioms.
- **CROSS_AGENT_SYNTHESIS.md**: Refreshed narrative tying Furey Z2xZ2, Option D compositeness (28+35=63), 7D leakage, and Shulga as "dynamical engine". Good high-level story ("particles are stable eigenvalue domains of a single parent algebra"). Contains overclaims: "Lean 4 Exact Parameter Derivation", "staggering 40-digit precise bounds", "analytically outputs", "WZW Selection" as realized — while delivered Shulga is finite lmax sum + guard and no WZW/functional integral is in the Lean artifact.
- **StabilityFromAlgebra.lean** (core living pipeline, ~250 lines): Rat (ℝ), real `octMult` from Fano table, `rev`, `rhoM` (scalarPart(c * rev c / 2)), `V_living_project(c, lam, mu, kappa, rho_amb, rho_crit)` = f(ρ_amb) * quad + (mu/2)*rhoM/(1+kappa*rhoM) using project's CONCEPT.md defaults (kappa=50, mu=-(8269/200), rho_crit=5/2). `computeHessianLiving` (finite-diff on rich V), `mixedBackground(fAmp)`, `leaksUnderLProtection`, `simulateCrossover` (still on old `V_full` quadratic, not living). **No Shulga, no geometric_ratio**. Phenomenological constants hardcoded in `V_living`, convenience `V_living`, combined scores (PhaseC), permanent table, and all certs.
- **PhaseC_Certificates.lean**: Excellent ingestion + cert layer. `python_report_consistency` + living scores (`lCombinedScoreLiving`/`lfCombinedScoreLiving` using Gershgorin + minDiag on `computeHessianLiving`). Direct pure-Lean `IO.FS.readFile` + minimal parsers for `python_famp_055_slice.json` (13.333% cross) and `mature_row_lambda012.json` (identical L/LF at λ=0.012). Named `matureReport012Consistency : true`, `golden_consistency_055_living`, `famp_table_consistency_living`. Permanent `#eval realDataComparisonTable` and `livingVsPythonComparison` (hardcoded numbers but matching actual JSONs and living V output). Calibrated bands for rich V (~18% relative L degradation). All certs evaluate true with current tuning.
- **Notes (2026-05-25 series)**: Consistent "still incomplete" honesty (e.g., "still incomplete per full end-to-end", "complete L=28/F=35 generator basis", "raw-data thresholds", "full retarded living-candidate V", "symbolic diff path"). Documented every calibration, pipeline evolution (shell → direct IO.FS), and criteria re-check after each edit. Verification discipline followed.
- **shulga_integration/**: Vision doc ambitious (Berry phase on S^3/Spin(8)/Spin(7) coset → Green function dressing for 2/9 phase, functional dets for 21/16 & 5 prefactors, WZW for Z2xZ2 selection). Python/Maxima side (gegenbauer sums, 7D params) exists but Lean slice is only the 1D-analog harmonic sum. No bridge module wiring the ratio into V or certs.
- **Build/JSON verification**: lake build clean. JSONs confirm the 13.333% and identical-peak cases. Table in source matches reality.

## Criteria Check vs INTEGRATION_PLAN Step 4 + History Litmus (Honest)
- "Computed from the algebra and then compared" — **Partial**. Leakage, L-not-closed, d-quark F-requirement, u-quark LF need, fAmp crossover *mechanism* are now derived from explicit Fano table + real octMult/rev/rhoM. But the *scale* (well depth, crossover threshold ~0.4, absolute scores in permanent table) is still driven by phenomenological kappa/rho_crit/lam/mu chosen to match Python, not the Shulga geometric_ratio. Shulga derivation sits beside the pipeline, not inside it.
- "Certificates ingest/validate against concrete Python numbers" — **Good for controlled cases**. Two real files, direct loaders, living V, permanent build-time table, named true certs. But only for slices the team already knew; no general arbitrary-report loader; thresholds calibrated to current V (circular).
- "Clear path for symbolic differentiation" — **Absent**. Finite-diff only; no stub, no doc, no plan in notes.
- "Full living-candidate V (retarded causal... exact from Python living-candidate)" — **Not present**. Current is best static approximation (f(ρ_amb) + kappa sat on real rhoM). No retarded lattice, no f_g/f_em on source vs probe, no full causal structure.
- "Complete L=28/F=35 generator basis" — **Partial**. Diag/sector checks use full 21 L + F sets; Hessian path uses representative mixedBackground + applyMask. Not all 63 operators exercised in stability scores.
- "Living thresholds from raw report data (not model-tuned)" — **Not done**. Bands hand-adjusted to ~18% observed under current V.
- "General pure-Lean JSON loader" — **Not done**. Excellent minimal parsers for the two specific controlled files; no schema-driven or robust arbitrary-report ingest.
- "simulateCrossover updated to rich V" — **Not done**. Still quadratic path.
- Additional litmus ("model vs derive", "forcing structural vs observation", "derived params actually drive the certs/tables"): **Shulga is the exact item-2 deliverable for 'phenomenological vs derived V' and is not yet integrated**. Forcing improved but remains enumerative rfl. Synthesis has aspirational language.

**Overall verdict**: Substantive, high-quality progress (Rat exactness, real rhoM in V, direct ingestion, permanent tables, stronger PhaseB, Shulga derivation attempt). The living-candidate pipeline + certs are now credible artifacts. However, the core "forcing from algebra + derived dynamics" promise of the 4-phase + INTEGRATION_PLAN is only partially realized — the highest-leverage piece (Shulga as engine for λ/μ) remains un-wired, forcing is still "check" not "prove", and Step 4 checklist has several open items. Work is worth keeping; integration is the next required cycle.

## Prioritized Gaps (Actionable, With Severity)
1. **Critical — Shulga not wired into living-candidate V / certs / tables (highest leverage)**  
   geometric_ratio / lambda_raw_exact / mu_raw_exact are computed but never constrain or replace the phenomenological constants in V_living_project, computeHessianLiving, l/lfCombinedScoreLiving, realDataComparisonTable, or any #eval cert. The "derived" parameter does not drive the machine-checkable statements that were the entire point of the long continue cycles. Finite lmax=100 + loose guard is a proof-of-concept, not yet promoted to production params.

2. **High — Forcing theorems remain enumerative/Boolean rather than structural**  
   PhaseB theorems are rfl on exhaustive generator lists. No lemma of the form "∀ M ∈ L-grade, for color-triplet indices, the Cartan eigenvalues produced by M on the N=1 irrep are identically zero because [Fano multiplication table property on bivector × color plane]". Consequently, "wrong-grade assignment → instability" is observed on samples, not derived as an identity across the dynamics.

3. **Medium — Synthesis overclaims vs delivered artifacts**  
   "40-digit precise", "analytically outputs", "WZW Selection", "exact theoretical value" language in CROSS_AGENT_SYNTHESIS.md exceeds what ShulgaParameters.lean + current Lean formalization actually contain (finite sum + guard; no functional integral, no WZW term formalized).

4. **Medium — INTEGRATION_PLAN Step 4 incomplete per its own written criteria + note assessments**  
   - simulateCrossover still quadratic (not living).  
   - No complete 28/35 operator set in Hessian path.  
   - Thresholds model-tuned, not raw-data calibrated.  
   - No general JSON loader.  
   - No symbolic diff path or stub.  
   - Full retarded living V not addressed.  
   - Permanent table and certs are excellent but rest on phenomenological inputs.

5. **Low — Minor code hygiene**  
   `crossoverLeakageDemo` contains `overlap != 0 || true` (always-true debug remnant). simulateCrossover not updated to living path. Shulga raw values printed in eval_test but not used. Parallel older artifacts in density_algebra/lean/ (StabilityBounds etc.) may duplicate effort.

6. **Contextual — Bridge to broader v59**  
   Shulga vision in shulga_integration/ (Berry on S^3, functional dets for prefactors, WZW for selection) is ambitious and aligns with CONCEPT.md / ROADMAP frontiers, but the Lean slice is narrow. No documented plan to lift the geometric_ratio into the full living-candidate or into Brannen kernel ξ dynamics.

## Verification Discipline Applied
- Explicit criteria re-check after every major read (see "Criteria Check" section).
- Honest "still incomplete" language preserved (matching the 2026-05-25 note series).
- No overstatement of progress; gaps are precise and mechanical.
- Build verified clean; numbers cross-checked against JSONs.
- Paths respect v59-in-root.
- This note itself is the "write it down... check... if incomplete continue" artifact.

## Proposed Continuation Steps (Prioritized, Mechanical, Cycle-Ready)
1. **Immediate (highest leverage)**: Wire Shulga into the pipeline.  
   - Add `import ShulgaParameters` (or a thin bridge) in StabilityFromAlgebra.  
   - Define `V_living_derived` (or param override) that uses geometric_ratio (or lambda_raw/mu_raw scaled to project mu scale) for lam/mu, or at minimum #guard-constrained bounds.  
   - Re-run / extend `realDataComparisonTable`, `livingVsPythonComparison`, all living certs (`matureReport012Consistency`, famp loaders, golden/famp living).  
   - New note: "2026-05-26-ShulgaWiredIntoLivingV.md" with before/after numbers and criteria re-check. If certs still pass (or surface honest new feedback), promote; else document the tension.

2. **Next**: Strengthen one PhaseB theorem to structural.  
   - Pick dquark_requires_F_grade. Prove (not check) that any L-grade bivector product on the color-plane N=1 indices yields identically zero Cartan diagonals, using the explicit Fano octMultTable properties (no diagonal contributions from specific bivector pairs).  
   - Add a named theorem whose proof depends on the multiplication table, not just enumeration.

3. **Tighten synthesis**: Edit CROSS_AGENT_SYNTHESIS.md to replace aspirational claims with precise delivered status ("finite-lmax Rat approximation with #guard bounds in (-1/4000,-1/6000)"; "WZW selection remains a proposed lifting from the Shulga vision document, not yet formalized in Lean").

4. **Continue Step 4 cycle**:  
   - Update simulateCrossover to call the living Hessian path.  
   - Add at least two more explicit L/F generators to the Hessian exercised set (move toward full 28/35).  
   - Implement a minimal "raw threshold finder" that reads a report row and emits the min L/LF degradation that would still be called 'consistent'.  
   - Write a short "PathToSymbolic.md" stub (even if just "finite-diff today; Lean 4 symbolic calculus library needed for exact d²V").  
   - Audit density_algebra/lean/ for duplication vs 7D_Algebra/; decide merge or deprecate.

5. **Optional high-value**: Once Shulga is wired and one structural forcing thm exists, produce a short "Derived vs Phenomenological" delta note showing how the geometric ratio changes (or preserves) the permanent table and certs.

This gap analysis closes the current exploration loop. The work is solid and the direction correct; the remaining distance is mechanical integration and one layer deeper on "derive" vs "observe". Ready for user direction on which gap to close first (recommend Shulga wiring).

**End of note**. If instructed "continue", the next action will be Step 1 above (Shulga wiring + re-validation + new dated verification note).