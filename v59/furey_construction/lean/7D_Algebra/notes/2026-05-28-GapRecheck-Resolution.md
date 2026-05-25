# Gap Recheck — Resolution pass (2026-05-28, follow-up to GapRecheck.md)

Addresses the concrete, code-correctness items raised in `2026-05-28-GapRecheck.md`.
Verified every claim against the live tree/build before acting; fixed the unambiguous
shipping bugs; reconciled the leakage narrative; tidied scratch. Larger research-scale
gaps are characterized honestly (NOT closed). Package builds clean.

## Fixed

1. **Prefactor desync (the "new gap" + gap 4) — FIXED.** Verified: `shulga_abs_prefactor = 1/400`
   in code, but several comments/strings said "1/16".
   - `ShulgaParameters.lean`: the two stale doc comments ("1/16 as a first structural guess",
     "1/16 is the first structural choice") rewritten to state `1/400` honestly as the single
     remaining TUNED external constant (ratio is algebraic; absolute scale is not yet derived).
   - `PhaseC_Certificates.lean`: the doc block (was "currently 1/16…") and the permanent
     `threeWayShulgaComparison` string (was `"(1/16) * G_S7(cutoff)"`) fixed. The string now
     prints `toString shulga_abs_prefactor` programmatically → **cannot desync again**. Build
     output confirms it now reads `1/400`.

2. **Leakage evidence "tension" (gap 2 / new gap 2) — RESOLVED (was not a contradiction).**
   Verified by eval: `crossoverLeakageDemo = true`, `crossoverLeakageOverlap = 0`,
   `diag(M²) = [-2,-2,-2,-2,0,0,0,0]`.
   The non-uniform diagonal IS the leakage: since `Z2Z2Forcing.L_grade_diagonal_free` proves the
   *entire* L-grade (Λ²⊕Λ⁶) is diagonal-free, and a scalar (Λ⁰) gives only a *uniform* diagonal,
   a non-uniform diagonal can come ONLY from the F-grade (Λ⁴). So squaring a mixed L source
   provably produces F content. The `crossoverLeakageOverlap = 0` is just Frobenius-orthogonality
   of M²'s *diagonal* F-content to the *off-diagonal* representative γ₀γ₁γ₂γ₃ — the wrong probe.
   - **Elevated the leakage from a narrow witness to a general structural forcing statement**
     (Z2Z2Forcing): the single `M²` computation was replaced by
     - `L_bivectors_product_eq_F` (no axioms): `(γ₀γ₁)(γ₂γ₃) = γ₀γ₁γ₂γ₃ = F_fourform_0123`;
     - `L_not_closed_reaches_F`: **L is not closed under multiplication — `L·L` reaches the
       F-grade** (the product of two disjoint bivectors is the 4-form on their union);
     - `composite_color_requires_LF`: the forcing combiner — (i) pure L has no color Cartan
       (`L_grade_diagonal_free`), (ii) products of L reach F (`L_not_closed_reaches_F`),
       (iii) F carries a rank-≥2 color Cartan (`F_color_cartan_rank_ge_two`); hence a
       color-charged *composite* (N=2 u-quark = product of L excitations) is pushed into F
       and must occupy L⊕F, not pure L.
   - Rewrote the notes on `crossoverLeakageOverlap` (StabilityFromAlgebra) and the PhaseB
     leakage theorem to point to this structural statement (the `M²` non-uniform diagonal is
     now just the special case `crossoverLeakageDemo`, whose cross term is `γ₀γ₁γ₂γ₃ ∈ F`).

3. **New top-level files (gap "new experiment files") — TRIAGED.**
   - `Z2Z2Forcing.lean`: NOT scratch — it is the task-(1) algebraic-forcing module (in the lakefile
     roots, builds, axiom-light). Documented in `../../notes/2026-05-24-phaseC-formalization.md`.
   - `export_shulga.lean`: legitimate — it is the `lean_exe shulga` root.
   - `test_forall.lean`, `eval_test.lean`: scratch leftovers (the former duplicated PhaseB's
     `L_bivectors_strictly_off_diagonal`; the latter bare `#eval`s). **Removed.**

## Still open (genuine, larger gaps — characterized, not closed)

- **Absolute Shulga scale not derived.** `1/400` is now clearly the single tuned external
  constant (no longer hidden by a desynced "1/16"). Deriving it needs the functional
  measure/Jacobian of the full S^7 coset — research-scale.
- **PhaseB forcing still partly enumerative — but materially advanced.** `gamma_mul_diag_zero`
  is universal over Fin 7; L-diag-free is now complete over all 28 generators; and the
  **compositeness direction is now an algebraic forcing argument** (`composite_color_requires_LF`):
  L not closed → `L·L` reaches F → F (not L) carries the color Cartan ⇒ a color-charged
  composite must be L⊕F. What still remains: the *full* "why each Fock sector takes its bits"
  (the Witt `N → SU(3)`-rep map + G₂ branching, e.g. why lepton is L and not F) — see task (1)
  in the 2026-05-24 note.
- **Step-4 items:** no complete L=28/F=35 basis wired into `computeHessianLiving` (uses
  representative masked sources); living thresholds are model-tuned, not raw-data-fitted;
  no general pure-Lean JSON loader (two controlled parsers only); no symbolic-diff
  implementation (finite-diff Hessian; `PathToSymbolic.md` is a stub); no full retarded
  causal living-candidate V.
- **Synthesis overclaims** (`CROSS_AGENT_SYNTHESIS.md`) and **density_algebra duplication** —
  docs/hygiene, not touched this pass.

## Follow-up: universal grade law + an octonionic-rep finding (same day)

Per the request to make the forcing "universal and non-enumerative, not 'we checked two'":

- **`lean/CliffordBladeGrade.lean`** (main package, axiom-clean): the universal,
  non-enumerative grade law. A basis blade = its index set; the blade product = symmetric
  difference; `grade = card`. Proved `grade_bladeMul` (`|S△T| = |S|+|T|−2|S∩T|`),
  `grade_bladeMul_disjoint`, and the headline `disjoint_bivectors_mul_isF` /
  `L_not_closed_on_disjoint_bivectors`: **for ALL disjoint bivectors, the product is
  grade-4 = F and is not L** — derived from `Finset.card_union_of_disjoint`, i.e. the grade
  law `2+2=4`, universally. This replaces the enumerative `Z2Z2Forcing.L_not_closed_reaches_F`
  as the structural statement.

- **~~IMPORTANT FINDING~~ → SUPERSEDED: it was a TABLE BUG (found + fixed 2026-05-28).**
  The "gammas don't Clifford-anticommute" observation, initially read as octonion
  non-associativity, was traced to **two sign errors in `octMultTable`** (the project's
  declared "single source of truth"):

  | entry | buggy table | correct (per its own Fano lines) |
  |---|---|---|
  | `e₁·e₆` | `+e₇` | `−e₇` |
  | `e₁·e₇` | `−e₆` | `+e₆` |

  Diagnosis: genuine octonions are **alternative**, and the linearized left-alternative law
  `a(bx)+b(ax)=(ab+ba)x` forces all imaginary left-mults to anticommute (a genuine Cl(7)).
  Only `e₁`'s pairs failed → not alternative → not octonions. The two entries even
  contradicted the table's own Fano line `(1,7,6)`. **Fixed both signs.** Now (verified):
  the table satisfies identity, squares=−1, anticommutativity, AND left-alternativity — a
  genuine octonion algebra — and the seven gammas form a **genuine Cl(7)** (all anticommute,
  square −I).

  **Regression tests added** (your suggestion): `SevenDAlgebra.oct_identity`,
  `oct_squares_neg_one`, `oct_anticommutative`, `oct_left_alternative`,
  `octMultTable_is_octonion`; matrix-level `PhaseB.gamma_anticommute`, `gamma_sq_neg_id`.
  These would have caught the bug and now guard against recurrence.

  **Consequences of the fix (every matrix changed):**
  - ✅ **Upside (big):** the gammas are now genuine Clifford generators, so the `Λ²/Λ⁴`
    ("L"/"F") labels ARE real Clifford grades — the operator-level grade bridge is now
    *rigorous* (the "needs the unformalized `Cl(7)_even ≅ ℂ⊗𝕆` iso" caveat was an artifact
    of the bug and is removed). `CliffordBladeGrade` applies directly.
  - ✅ **Core forcing survives:** `L_grade_diagonal_free` (all 28 L gens diagonal-free),
    F color-Cartan rank ≥ 2 (`F_A/F_B` diagonals `[1,1,-1]`,`[-1,1,-1]` unchanged),
    `color_splitting_requires_F`, `composite_color_requires_LF`, the `CliffordBladeGrade`
    bridge — all still hold.
  - ⚠ **Casualties (the crossover/leakage narrative was bug-driven):**
    `crossoverLeakageDemo` flipped `true→false` (genuine grade-4 element is diagonal-free);
    PhaseB's two leakage theorems were reframed (`crossoverLeakageDemo_false_on_correct_table`).
    In PhaseC, **9 crossover/protection-stacking certificates flipped `true→false`** — the
    f-amplitude crossover / "L-protection degrades, LF preserved" signal does NOT hold on the
    correct algebra; it was an artifact of the spurious non-uniform diagonal. Recorded
    honestly as `= false`; **the stability-bound narrative needs full re-derivation against
    Python re-run on the corrected algebra.**
  - 🔧 **`gamma_octonionic_not_anticommuting` removed** (it asserted `≠ 0`, now false); its
    docs in `Z2Z2Forcing`/`StabilityFromAlgebra`/`PhaseB` rewritten to the corrected picture.

  **Propagation of the fix (done 2026-05-24):** the same 2 sign errors were found and fixed
  in all copies — `density_algebra/lean/OctonionAlgebra.lean`, and the Maxima
  `octonion_sensitivity_analysis.mac`, `test_oct_mult.mac`, `full_octonion_perturbation.mac`
  (Maxima itself confirms the corrected table is anticommutative + left-alternative). No
  real Python carries the table: the density Python uses a Cl(3,0) GA library and
  `13_fock_mass_forcing.py` uses exterior/interior products; only throwaway `7D_Algebra/*.py`
  exploration scratch hold buggy copies.

  **Do the `.mac` numbers change? NO (verified by re-running both ways).** The two analysis
  `.mac` files (`octonion_sensitivity_analysis`, `full_octonion_perturbation`) give bit-for-bit
  identical output before/after the fix. Reason: they use only **scalar/norm** quantities
  (`scalar(M·M)=m₀²−Σmᵢ²`, `|M|²=Σmᵢ²`), and the scalar part of `eᵢeⱼ` is `−δᵢⱼ` regardless of
  the imaginary-part sign bug. So **no downstream Maxima-derived numbers need correction.**
  (Contrast: the Lean PhaseC certs DID flip because they use `|M·M|²` — the full product norm,
  including the imaginary part — and full-matrix diagonals, which are sign-dependent.)

- **`Z2Z2Forcing` now DEPENDS on the abstract law** (rigorous grade bridge). `Z2Z2Forcing`
  imports `CliffordBladeGrade` and derives the forcing grade fact from it:
  - `pairBlade`/`quadBlade` map the v59 generator labels to Clifford blades (index sets);
  - `pairBlade_grade` (a bivector label is grade 2), `pairBlade_bladeMul` (disjoint pair
    labels' blade product = the fourform label `{i,j,k,l}`);
  - `disjoint_pairBlades_are_F` — for ANY disjoint bivector labels the product is F-grade —
    proved as a *direct instance* of `CliffordBladeGrade.disjoint_bivectors_mul_isF` (no
    enumeration). So the grade/index identification between the two modules is now rigorous
    and universal.
  - **PACKAGING.** Because the bridge needs Mathlib + `CliffordBladeGrade`, `Z2Z2Forcing`
    now builds in the **parent** package (`lake build Furey7DAlgebra` from `lean/`), not the
    standalone Mathlib-free `7D_Algebra` package. Fixed the parent lakefile (added the
    missing `ShulgaParameters` root); removed `Z2Z2Forcing` from the standalone roots (which
    keeps the Mathlib-free octonion-matrix core). Also made one PhaseC file-IO `#eval`
    robust (`<|>` fallback) so the parent build — with a different CWD — doesn't fail on a
    relative JSON path.
  - **Still the irreducible gap:** the *operator-level* grade iso `Cl(7)_even ≅ ℂ⊗𝕆` (tying
    the octonion-left-mult matrix products to these blade labels) is NOT formalized; the
    rigorous bridge is at the grade/index-label level.

## Build / axioms
All `7D_Algebra` modules + the `shulga` exe build clean. New `mixed_L_source_leaks_to_F`
depends on no axioms. (PhaseC `native_decide` certs still additionally use `Lean.ofReduceBool`,
as before.)
