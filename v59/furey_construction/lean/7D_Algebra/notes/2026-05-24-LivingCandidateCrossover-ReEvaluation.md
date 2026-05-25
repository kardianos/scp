# Re-Evaluation: the Living-Candidate Crossover work, post octonion-table fix

**Date**: 2026-05-24
**Trigger**: the `octMultTable` sign-bug fix (`e₁·e₆`, `e₁·e₇`) made the algebra a genuine
Cl(7) and flipped 9 PhaseC crossover/protection certs `true→false`. This note re-evaluates the
*entire* living-candidate / crossover / stability-bound thread, retracts what was bug-driven,
keeps what survives, and sets a corrected plan. Supersedes the crossover claims in the
`2026-05-25-Step4-*` notes, the Step-4 sections of the GapAnalysis/GapAudit notes, and the
"re-run certs" gate of the S7-coset sub-project.

## 1. The decisive finding: two different algebras were conflated

The f-amplitude crossover ("L-only protection degrades as F-content rises; LF-stacking
preserves the well") originates in **Python** (`density_algebra/python/sweep_f_amplitude.py`,
`improved_density_protection_scan.py`). Those scripts compute with **`ga.MV` — a Cl(3,0)
geometric-algebra library** (basis `[s, e1,e2,e3, e12,e13,e23, e123]`), i.e. the **v58
unified-multivector model**. The crossover is therefore a **Cl(3,0) phenomenon**, computed
correctly there.

The Lean `7D_Algebra` PhaseC certificates tried to **reproduce that Cl(3,0) crossover inside
the octonion algebra** (`octMultTable`, `simulateCrossover`, the living-V Hessian). The
"match" worked only because the **buggy** octonion table produced a spurious non-uniform
diagonal that happened to mimic the Cl(3,0) signal. On the **corrected** octonion algebra the
match disappears (9 certs `→ false`, `crossoverLeakageDemo → false`).

**Conclusion:** the octonion algebra does **not** exhibit the f-amplitude crossover. The
crossover (if physically meaningful at all) belongs to the **Cl(3,0) layer**, not the
octonionic Furey layer. The certificate enterprise was bridging two distinct algebras and was
held together by the table bug.

## 2. Inventory & classification

### A. SURVIVES — genuine, correct-algebra, KEEP
- Genuine **Cl(7)** octonion algebra: fixed table + regression tests
  (`oct_anticommutative`, `oct_left_alternative`, `octMultTable_is_octonion`,
  `gamma_anticommute`, `gamma_sq_neg_id`).
- **Structural Z₂×Z₂ forcing** (now operator-rigorous on genuine Cl(7)):
  `L_grade_diagonal_free`, F color-Cartan rank ≥ 2, `color_splitting_requires_F`,
  `composite_color_requires_LF`, the `L_bivectors_product_eq_F` / `L_not_closed_reaches_F`
  facts, and the `CliffordBladeGrade` universal grade law + bridge.
- PhaseB grade-structure facts: `gamma_mul_diag_zero`, lepton/d/u sector diagonals,
  grade-intensity signs.
- **Shulga ratio** (`ShulgaParameters.geometric_ratio`): Gegenbauer S⁷ harmonic sums,
  independent of `octMultTable` — bug-immune.
- The Maxima density-well analyses (`octonion_sensitivity_analysis.mac`,
  `full_octonion_perturbation.mac`): use scalar/norm quantities → bug-immune (verified:
  identical output pre/post fix). NB: these are generic density-well analyses, not the
  f-amplitude crossover.

### B. DEAD — bug-driven octonion "reproduction" of a Cl(3,0) phenomenon, REMOVE/DEPRECATE
- `crossoverLeakageDemo`, `crossoverLeakageOverlap`, `simulateCrossoverLiving`,
  `mixedBackground`, the `computeHessianNum/Living` + Gershgorin "protection score" layer in
  `StabilityFromAlgebra.lean` (built to match the Cl(3,0) crossover via the octonion algebra).
- The PhaseC crossover/protection certificates: the 9 now-`false`
  (`python_crossover_cert_*`, `golden_consistency_055`, `famp_table_consistency`,
  `combined_cert_*`, `gershgorin_cert_*`) **and** the 7 living-V relative-difference certs —
  the latter "survive" only because their `|L−LF|/avg` bands are too loose to detect the
  vanished crossover; they are not evidence of anything and should not be presented as such.
- The Python-data ingestion machinery in PhaseC (the JSON loaders, golden tables) — it
  ingests Cl(3,0)-model numbers to "certify" the octonion algebra; the comparison is
  ill-founded.

### C. REAL BUT IN A DIFFERENT MODEL — clarify, don't conflate
- The Cl(3,0) Python crossover itself (`sweep_f_amplitude.py` etc.). Correctly computed in
  the v58 multivector model. If it is a result worth keeping, it lives in the **v58 / Cl(3,0)
  layer** and should be studied there. The octonion `7D_Algebra` layer should make **no claim**
  to reproduce or certify it.

## 3. Updated expectations (retractions)

- **RETRACT** the INTEGRATION_PLAN core promise that "the stability bounds are computed from
  the [octonion] algebra and then compared to Python" / "machine-checkable certificates that
  the explicit 7D algebra forces the f-amplitude crossover." The octonion algebra does not
  force the crossover; that promise was the table bug.
- **RETRACT** the protection-stacking / "L leaks, LF preserved" narrative as an *octonionic*
  result. (Its status in Cl(3,0) is a separate, v58-layer question.)
- The S7-coset sub-project's Indicator D ("wire derived prefactor in and re-run the
  crossover certs") is moot for the crossover certs (gone). Its absolute-scale half is gated
  behind there being a real correct-algebra stability structure to scale — there isn't one in
  the octonion layer. Its `21/16` & `5` half (Indicator E) remains valid and bug-immune.

## 4. Plan

### Remove / reframe (Lean `7D_Algebra`)
1. **PhaseC_Certificates.lean** — remove the crossover/stability certificate suite and its
   Python-ingestion machinery; keep nothing that "certifies" the octonion algebra against
   Cl(3,0) crossover data. (≈most of the file.)
2. **StabilityFromAlgebra.lean** — remove the living-V / Hessian / Gershgorin / simulateCrossover
   / leakage-demo layer. Keep only the genuine algebra helpers still used by the structural
   modules (`matMul`, `matAdd`, `octMult`, `diag`, etc. as needed by PhaseB/Z2Z2Forcing).
3. **PhaseB_Theorems.lean** — drop the leakage-demo theorems entirely (already reframed to
   `crossoverLeakageDemo_false_on_correct_table`; that can go too). Keep the grade-structure
   theorems and the new `gamma_anticommute`/`gamma_sq_neg_id` regression tests.
4. **ShulgaParameters.lean** — keep the ratio; drop `shulga_abs_*` (the tuned absolute scale)
   or clearly quarantine it, since it only fed the now-removed living-V.

### Reframe (docs)
5. Update `INTEGRATION_PLAN.md` Step-4 success criteria: retract the crossover-reproduction
   goal; state the genuine deliverables (structural forcing + Shulga ratio).
6. Mark the `2026-05-25-Step4-*` notes and the crossover sections of the GapAnalysis/GapAudit
   notes as **SUPERSEDED** (header pointer to this note) — keep as history, do not delete.
7. `CROSS_AGENT_SYNTHESIS.md` / synthesis docs: remove protection-stacking-as-octonion claims.

### Keep / forward
8. Keep all of §2.A. The 7D_Algebra package's real contribution is the **structural Z₂×Z₂
   forcing on genuine Cl(7)** + the `CliffordBladeGrade` bridge — that is solid and rigorous.
9. Forward targets (bug-immune, high-value):
   - extend the structural forcing toward the full Fock/Witt `N → SU(3)`-rep + G₂-branching
     (the remaining open part of task (1));
   - the Shulga ratio + the `21/16`,`5` coset-Jacobian *derivation* (S7-coset sub-project,
     Indicator E only);
   - if the Cl(3,0) crossover matters, pursue it in the v58 layer, not via octonions.

## 5. One-line summary
The living-candidate crossover was a **Cl(3,0)** phenomenon that the octonion Lean only
appeared to reproduce because of a multiplication-table sign bug. With the table corrected,
that whole sub-thread is retracted from the octonion layer; what remains (and is now more
rigorous) is the structural Z₂×Z₂ grade forcing and the bug-immune Shulga ratio.
