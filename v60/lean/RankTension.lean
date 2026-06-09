/-
Copyright (c) 2026  SCP Project
Released under Apache 2.0 license.

v60/lean/RankTension.lean

G1 — the rank tension of the EW bridge v_Higgs = dim(L)²·a² = 784 a².

Machine-checked arithmetic backbone of the analysis in
`v60/gaps/rank_tension/01_rank_tension.py`:
  * the dimension facts (784 = dim(L)², 25 = dim(L) − N_gen, max proper
    subalgebra of so(8) = 21 < 25);
  * the structural numerator 6 = N_gen² · Q (Q = 2/3 = dimG₂/dimSpin7);
  * the gravity↔EW ratio 9Q/784 = 6/784;
  * the DEFLATION as an exact algebraic identity: (Σm/v)/(9Q/784) = (784 a²)/v,
    i.e. the advertised "Σm = (6/784)v" match is *exactly* the EW bridge match
    v = 784 a² re-expressed (Σm = 9Qa² being definitional), carrying no
    independent information.

Builds against the v59 Mathlib:
  cd v59/furey_construction/lean && lake env lean ../../../v60/lean/RankTension.lean
-/

import Mathlib

namespace SCPv60.RankTension

/-! ## 1. Dimensions -/

def dimL : ℕ := 28          -- lepton ambient L = Λ²⊕Λ⁶ = dim so(8) = dim Spin(8)
def dimEndL : ℕ := dimL ^ 2 -- operator algebra End(L) = M₂₈(ℝ) (Burnside) → the 784
def Ngen : ℕ := 3
def dimG2 : ℕ := 14
def dimSpin7 : ℕ := 21      -- max proper subalgebra of so(8)

theorem endL_dim : dimEndL = 784 := by decide
theorem nine_is_Ngen_sq : Ngen ^ 2 = 9 := by decide

/-- **[thm] Koide numerator is structural:** `Q = 2/3 = dimG₂/dimSpin7`. -/
theorem koide_struct : (dimG2 : ℚ) / dimSpin7 = 2 / 3 := by
  unfold dimG2 dimSpin7; norm_num

/-! ## 2. The gravity↔EW ratio and its numerator -/

/-- `Σm = 9 Q a²` with `9 = N_gen²`; for `Q = 2/3` the numerator is `9·(2/3) = 6`. -/
theorem ratio_numerator : ((Ngen ^ 2 : ℚ)) * (2 / 3) = 6 := by
  unfold Ngen; norm_num

/-- **[thm] `Σm/v = 9Q/dim(L)² = 6/784`** (for `Q = 2/3`). -/
theorem grav_ew_ratio : ((Ngen ^ 2 : ℚ) * (2 / 3)) / (dimEndL : ℚ) = 6 / 784 := by
  unfold Ngen dimEndL dimL; norm_num

/-! ## 3. The deflation (exact algebraic identity)

`Σm = 9 Q a²` is DEFINITIONAL (`a = (Σ√m)/3`, `Q` = Koide).  Hence the ratio
`(Σm / v)` compared to the structural `9Q/784` is just the EW bridge ratio: -/

/-- **[thm] Deflation.** `(Σm/v) / (9Q/784) = (784 a²)/v` — the "gravity↔EW 6/784
    bonus" is exactly the EW bridge match `v = 784 a²` re-expressed, with no
    independent content.  (Cf. v59's `g_W² = 5√α = α(M_Z)`-in-disguise.) -/
theorem deflation (a v Q : ℝ) (hv : v ≠ 0) (hQ : Q ≠ 0) :
    (9 * Q * a ^ 2 / v) / (9 * Q / 784) = (784 * a ^ 2) / v := by
  field_simp

/-! ## 4. The rank tension and the 25-direction obstruction -/

/-- The "missing" directions: `dim(L) − N_gen = 25`. -/
def missing : ℕ := dimL - Ngen
theorem missing_eq : missing = 25 := by decide

/-- **[thm] Frobenius² ≠ vector reading:** `dim(L)² = 784 ≠ 28 = dim(L)`.  The EW
    scale uses the Frobenius² (784), not the singlet/vector reading (28 = √784·√28). -/
theorem readings_differ : dimEndL ≠ dimL := by decide

/-- **[thm] Subalgebra obstruction:** the max proper subalgebra of `so(8)` has
    dim 21 < 25, so a single-step `so(8) → H` cannot eat exactly 25 Goldstones
    (leave a dim-3 stabilizer); "3 light + 25 eaten" is not single-step. -/
theorem subalgebra_obstruction : dimSpin7 < missing := by decide

/-- The two Frobenius² readings live on spaces of DIFFERENT dimension: the EW
    bilinear on `End(L)` (784) vs the Brannen kernel on the `N_gen`-space (here its
    operator space `N_gen² = 9`).  `784 ≠ 9` — they cannot be one matrix. -/
theorem different_spaces : dimEndL ≠ Ngen ^ 2 := by decide

/-! ## 5. Reading

The arithmetic above is the skeleton of `01_rank_tension.py`'s verdict: the
gravity↔EW "bonus" is the bridge match in disguise (`deflation`); a single `Y`
cannot be both full-rank-democratic (784) and rank-3 (the spectrum); and the
25-direction problem has no single-step `so(8)` home (`subalgebra_obstruction`).
The only consistent reading is two objects on different spaces — which leaves
`v_Higgs` a conditional second scale (tied to `a_ℓ` only via the still-homeless
R1+R2 conjectures).  None of this is a *derivation* of `v` from `a_ℓ`; it is a
sharp characterization of why that derivation is still open. -/

end SCPv60.RankTension
