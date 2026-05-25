/-
  v59/furey_construction/lean/AlphaZero.lean

  Task (4a): the v59 low-energy fine-structure conjecture.

  v59 conjectures (see SUMMARY / FINDINGS_scale_bridge) the transcendental relation

      −ln α + 2α  =  π²/2  =  8π² / dim Cl(3,1)            (dim Cl(3,1) = 16)

  as the structural origin of α(0).  This file machine-checks the two checkable parts:

    * `rhs_eq_structural`        — the RHS identity `8π²/16 = π²/2` (so the "π²/2" is
                                   `8π²/dim Cl(3,1)`, the structural integer 16).
    * `alpha0_conjecture_bracket`/`alpha0_conjecture_holds` — at the *measured* value
                                   `α = 1/137.036` the relation holds to within `1/100`
                                   (the actual gap is ≈ 2.4×10⁻³).

  HONESTY.  This does NOT derive `α = 1/137.036` (that remains the empirical input), nor
  does it motivate the equation `−ln α + 2α = π²/2` itself (which is conjectural).  It
  certifies that the measured α is consistent with the conjectured relation to the stated
  precision — turning a prose "(4×10⁻⁵)/0.004% match" claim into a checked bound.  The
  `1/100` tolerance is loose because the `log 137.036` bound here is crude (`log 2` ±
  `log(1+x)` bracketing); it comfortably contains the true gap.

  CORRECTION (2026-05-24, `ALPHA_SCOPING.md`): the true gap of `−ln α + 2α` from `π²/2` is
  `≈3.6×10⁻⁵` (NOT `2.4×10⁻³` as previously stated here).  BUT this tight match relies on the
  `+2α` term, which `04_alpha_prediction.py` (Part 3) **reverse-engineered** as the correction
  that makes `π²/2` fit ("what correction makes π²/2 exact?" → scanned `1/70,1/64,1/8π²,α,2α`
  → `2α` fit).  The mechanism-backed *pure* instanton form `−ln α = π²/2 = 8π²/16` gives only
  `α⁻¹ = e^{π²/2} = 139.0` — `1.4%` off.  So the `0.004%` headline rides on a fitted `+2α`;
  the instanton (`g²=16=dim Cl(3,1)`) is itself ungrounded (`π₃(S⁷)=0` ⇒ no S⁷-based instanton).
-/
import Mathlib

namespace SCPv59.AlphaZero

open Real

/-- The measured low-energy fine-structure constant `α(0) = 1/137.036`. -/
noncomputable def alpha0 : ℝ := 1000 / 137036

/-- **Structural RHS.** `π²/2 = 8π² / 16`, and `16 = dim Cl(3,1)` (`Predictions.dimCl31`).
    So the conjecture's right-hand side is `8π² / dim Cl(3,1)`. -/
theorem rhs_eq_structural : 8 * Real.pi ^ 2 / 16 = Real.pi ^ 2 / 2 := by ring

/-- `−ln α + 2α` at `α = 1/137.036` equals `ln(137.036) + 2/137.036`. -/
theorem lhs_eq :
    -Real.log alpha0 + 2 * alpha0 = Real.log (137036 / 1000) + 2 / (137036 / 1000) := by
  unfold alpha0
  rw [show (1000 / 137036 : ℝ) = (137036 / 1000)⁻¹ by norm_num, Real.log_inv]
  norm_num

/-- **The numerical bracket.** With `α = 1/137.036`, `−ln α + 2α` is within `1/100` of
    `π²/2`.  Proof: `log 137.036 = 7·log 2 + log(34259/32000)`, bounded by Mathlib's
    `log 2` bounds and `x − x²/(…) ≤ log(1+x) ≤ x − 1` style bounds, against the π² bounds
    from `Real.pi_gt_d6` / `Real.pi_lt_d6`. -/
theorem alpha0_conjecture_bracket :
    |(Real.log (137036 / 1000) + 2 / (137036 / 1000)) - Real.pi ^ 2 / 2| ≤ 1 / 100 := by
  have hdec : Real.log (137036 / 1000 : ℝ) = 7 * Real.log 2 + Real.log (34259 / 32000) := by
    rw [show (137036 / 1000 : ℝ) = 2 ^ 7 * (34259 / 32000) by norm_num,
        Real.log_mul (by positivity) (by norm_num), Real.log_pow]
    push_cast; ring
  have hl2u := Real.log_two_lt_d9
  have hl2l := Real.log_two_gt_d9
  have hyu : Real.log (34259 / 32000 : ℝ) ≤ 34259 / 32000 - 1 :=
    Real.log_le_sub_one_of_pos (by norm_num)
  have hyl : (1 : ℝ) - 32000 / 34259 ≤ Real.log (34259 / 32000) := by
    have h := Real.log_le_sub_one_of_pos (show (0 : ℝ) < 32000 / 34259 by norm_num)
    rw [show (32000 / 34259 : ℝ) = (34259 / 32000)⁻¹ by norm_num, Real.log_inv] at h
    linarith
  have hpisq_l : (3.141592 : ℝ) ^ 2 ≤ Real.pi ^ 2 := by nlinarith [Real.pi_gt_d6, Real.pi_pos]
  have hpisq_u : Real.pi ^ 2 ≤ (3.141593 : ℝ) ^ 2 := by nlinarith [Real.pi_lt_d6, Real.pi_pos]
  rw [hdec, abs_le]
  constructor <;> nlinarith [hl2u, hl2l, hyu, hyl, hpisq_l, hpisq_u]

/-- **The v59 fine-structure conjecture holds at the measured α(0), to within 1/100.**
    `|(−ln α + 2α) − 8π²/dim Cl(3,1)| ≤ 1/100` at `α = 1/137.036`. -/
theorem alpha0_conjecture_holds :
    |(-Real.log alpha0 + 2 * alpha0) - 8 * Real.pi ^ 2 / 16| ≤ 1 / 100 := by
  rw [lhs_eq, rhs_eq_structural]
  exact alpha0_conjecture_bracket

end SCPv59.AlphaZero
