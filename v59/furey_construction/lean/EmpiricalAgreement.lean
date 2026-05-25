/-
  v59/furey_construction/lean/EmpiricalAgreement.lean

  Task (4d): turn the v59 "agrees with experiment to X%" claims — which the README
  and findings state only in prose — into machine-checked theorems.

  For each v59 structural prediction that is a *rational* number we encode the measured
  empirical central value (as a rational, citing the source) together with an explicit
  relative tolerance `τ`, and prove

      |structural − empirical| ≤ τ · empirical          (`relDevLE`)

  i.e. the structural value agrees with the empirical central value to relative
  precision `τ`.  This is a purely arithmetic (ℚ) statement, axiom-clean.

  HONESTY LABELS.  Two distinct kinds of agreement appear, and the theorems below are
  grouped accordingly:
    • WITHIN EXPERIMENTAL ERROR — the structural value sits inside the measurement
      uncertainty (lepton Koide Q, Brannen phase φ).
    • STRUCTURAL APPROXIMATION — the structural value matches the central value to
      O(0.1–0.6%), which is *larger* than the (small) experimental error on the EW
      observables; these are suggestive, not within error.  (For the quark Koide
      ratios the ~0.3% gap is plausibly inside the large quark-mass-scheme uncertainty.)

  Each structural literal below is the value proven elsewhere in this development; the
  citing theorem is named in the doc-string (e.g. `2/3 = Predictions.koide_Q_lepton`).
-/
import Mathlib.Tactic

namespace SCPv59.EmpiricalAgreement

/-- `relDevLE s e τ` : the value `s` agrees with the (positive) empirical central value
    `e` to relative tolerance `τ`, i.e. `|s − e| ≤ τ·e`. -/
def relDevLE (s e tol : ℚ) : Prop := |s - e| ≤ tol * e

/-- All agreement theorems reduce to this one decision procedure. -/
private lemma relDevLE_check {s e tol : ℚ}
    (h : -(tol * e) ≤ s - e ∧ s - e ≤ tol * e) : relDevLE s e tol := by
  unfold relDevLE; rw [abs_le]; exact h

/-! ## Within experimental error -/

/-- **Lepton Koide ratio.** Structural `Q = 2/3` (`= Predictions.koide_Q_lepton`,
    `= dim G₂ / dim Spin(7)`).  Empirical `Q = 0.6666605` from the charged-lepton masses
    (PDG); the `m_τ`-dominated uncertainty is ≈ 10⁻⁴, and `2/3` sits ≈ 6×10⁻⁶ from the
    central value — **within experimental error**. -/
theorem koide_lepton_agrees : relDevLE (2/3) (6666605/10000000) (1/10000) :=
  relDevLE_check (by constructor <;> norm_num)

/-- **Lepton Brannen phase.** Structural `φ = 2/9` (`= ScaleBridge`/`KoideAndBrannen`,
    `= Q/3`).  Empirical `0.2222296` — within experimental error (gap ≈ 7×10⁻⁶). -/
theorem brannen_phase_agrees : relDevLE (2/9) (2222296/10000000) (1/10000) :=
  relDevLE_check (by constructor <;> norm_num)

/-! ## Structural approximations (O(0.1–0.6%), not within EW experimental error) -/

/-- **Weak mixing angle (sin²).** Structural `sin²θ_W = 2/9`
    (`= ScaleBridge.sin_sq_thW`, the Brannen phase).  On-shell empirical `0.2231`;
    relative gap ≈ 0.39% (tolerance 0.5%).  NOT within PDG error. -/
theorem sin2_thetaW_agrees : relDevLE (2/9) (2231/10000) (5/1000) :=
  relDevLE_check (by constructor <;> norm_num)

/-- **Weak mixing angle (cos²).** Structural `cos²θ_W = 7/9 = t²_(u-quark)`
    (`= ScaleBridge.cos_sq_thW`).  On-shell empirical `0.7770`; gap ≈ 0.10% (tol 0.2%). -/
theorem cos2_thetaW_agrees : relDevLE (7/9) (7770/10000) (2/1000) :=
  relDevLE_check (by constructor <;> norm_num)

/-- **W/Z mass ratio.** Structural `(m_Z/m_W)² = 1/cos²θ_W = 9/7`
    (`= ScaleBridge.mZ_over_mW_sq`).  Empirical `(91.1876/80.369)² ≈ 1.28729` (PDG masses);
    gap ≈ 0.12% (tol 0.3%). -/
theorem mZ_over_mW_sq_agrees : relDevLE (9/7) (128729/100000) (3/1000) :=
  relDevLE_check (by constructor <;> norm_num)

/-- **Higgs mass.** Structural `m_H²/v² = 7/27`.  Empirical `(125.25/246.22)² ≈ 0.258773`
    (PDG `m_H`, electroweak `v`); gap ≈ 0.19% (tol 0.3%). -/
theorem mH_sq_over_v_sq_agrees : relDevLE (7/27) (258773/1000000) (3/1000) :=
  relDevLE_check (by constructor <;> norm_num)

/-- **Down-type quark Koide.** Structural `Q_d = 11/15` (`= Predictions.koide_Q_d_quark`,
    `t²_d = 3/5`).  Empirical `Q_down = 0.731` (from `m_d, m_s, m_b`); gap ≈ 0.32%
    (tol 0.5%) — plausibly within the large quark-mass-scheme uncertainty. -/
theorem koide_dquark_agrees : relDevLE (11/15) (731/1000) (5/1000) :=
  relDevLE_check (by constructor <;> norm_num)

/-- **Up-type quark Koide.** Structural `Q_u = 23/27` (`= Predictions.koide_Q_u_quark`,
    `t²_u = 7/9`).  Empirical `Q_up = 0.849`; gap ≈ 0.34% (tol 0.5%). -/
theorem koide_uquark_agrees : relDevLE (23/27) (849/1000) (5/1000) :=
  relDevLE_check (by constructor <;> norm_num)

/-- **Cabibbo angle.** Structural `sin²θ_C = 7·α(0) = 7000/137036` (`= ScaleBridge.sin_sq_cabibbo`
    at `α(0) = 1/137.036`, `7 = dim Im𝕆`).  Empirical `0.0508`; gap ≈ 0.55% (tol 1%). -/
theorem cabibbo_agrees : relDevLE (7000/137036) (508/10000) (1/100) :=
  relDevLE_check (by constructor <;> norm_num)

/-! ## Consolidated table -/

/-- All nine v59 rational predictions agree with their empirical central values to the
    stated tolerances (two within experimental error, seven to O(0.1–0.6%)). -/
theorem v59_empirical_agreement_table :
    relDevLE (2/3) (6666605/10000000) (1/10000)
    ∧ relDevLE (2/9) (2222296/10000000) (1/10000)
    ∧ relDevLE (2/9) (2231/10000) (5/1000)
    ∧ relDevLE (7/9) (7770/10000) (2/1000)
    ∧ relDevLE (9/7) (128729/100000) (3/1000)
    ∧ relDevLE (7/27) (258773/1000000) (3/1000)
    ∧ relDevLE (11/15) (731/1000) (5/1000)
    ∧ relDevLE (23/27) (849/1000) (5/1000)
    ∧ relDevLE (7000/137036) (508/10000) (1/100) :=
  ⟨koide_lepton_agrees, brannen_phase_agrees, sin2_thetaW_agrees, cos2_thetaW_agrees,
   mZ_over_mW_sq_agrees, mH_sq_over_v_sq_agrees, koide_dquark_agrees, koide_uquark_agrees,
   cabibbo_agrees⟩

end SCPv59.EmpiricalAgreement
