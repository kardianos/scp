/-
  v59/furey_construction/lean/PhaseAmbiguity.lean

  **The decisive negative result for `φ = 2/9`** (narrows the design hard).

  The charged-lepton masses do NOT determine the Brannen phase `φ` — only its `S₃`-orbit
  (the generation `Z₃` shift `φ → φ + 2π/3` and conjugation `φ → −φ`), equivalently only
  `cos 3φ`.  Empirical check: fitting the PDG masses returns `φ ≈ 2.3166 rad = 2/9 + 2π/3`
  — a full `Z₃` shift away from "2/9".  So **"φ = 2/9" is a choice of representative**, not a
  physical observable; the genuine invariant is `cos 3φ = cos(2/3) ≈ 0.7859`.

  Consequences for the program:
    * The contrast with the amplitude is sharp: `t² = 1/2` is invariant and solid (it comes
      from the `S₃`-invariant Koide `Q`); `φ = 2/9` is branch-fixed and convention-dependent.
    * A *structural* derivation must target `cos 3φ` (the not-nice `cos(2/3)`), not the nice
      "2/9" — which only looks nice in the principal branch and in radians.  This undercuts a
      naive "derive the rational 2/9" mechanism.
-/
import Mathlib
import BrannenKernel
import BrannenPhase

namespace SCPv59.PhaseAmbiguity

open Real SCPv59.BrannenKernel

/-! ## The generation Z₃ shift permutes the amplitudes -/

/-- **`φ → φ + 2π/3` cyclically permutes the Brannen amplitudes** `(s₀,s₁,s₂) → (s₁,s₂,s₀)`.
    Hence it leaves the *multiset* of masses — and every symmetric function of them — unchanged. -/
theorem s_shift (a t φ : ℝ) :
    s a t (φ + 2 * π / 3) 0 = s a t φ 1
    ∧ s a t (φ + 2 * π / 3) 1 = s a t φ 2
    ∧ s a t (φ + 2 * π / 3) 2 = s a t φ 0 := by
  refine ⟨by simp only [s], ?_, ?_⟩
  · simp only [s]; rw [show φ + 2 * π / 3 + 2 * π / 3 = φ + 4 * π / 3 by ring]
  · simp only [s]; rw [show φ + 2 * π / 3 + 4 * π / 3 = φ + 2 * π by ring, Real.cos_add_two_pi]

/-- **Conjugation `φ → −φ`** also permutes the amplitudes `(s₀,s₁,s₂) → (s₀,s₂,s₁)` (a
    reflection).  Together with `s_shift`, the full `S₃` acts on `φ` leaving the masses fixed. -/
theorem s_conj (a t φ : ℝ) :
    s a t (-φ) 0 = s a t φ 0
    ∧ s a t (-φ) 1 = s a t φ 2
    ∧ s a t (-φ) 2 = s a t φ 1 := by
  refine ⟨by simp only [s, Real.cos_neg], ?_, ?_⟩
  · simp only [s]
    rw [show -φ + 2 * π / 3 = -(φ + 4 * π / 3) + 2 * π by ring, Real.cos_add_two_pi, Real.cos_neg]
  · simp only [s]
    rw [show -φ + 4 * π / 3 = -(φ + 2 * π / 3) + 2 * π by ring, Real.cos_add_two_pi, Real.cos_neg]

/-! ## Hence the masses determine only `cos 3φ` -/

/-- The third power sum (the *only* phase-dependent moment, `BrannenPhase.sum_s_cube`) is
    invariant under the `Z₃` shift: the masses cannot tell `φ` from `φ + 2π/3`. -/
theorem moments_Z3_invariant (a t φ : ℝ) :
    (s a t (φ + 2 * π / 3) 0) ^ 3 + (s a t (φ + 2 * π / 3) 1) ^ 3 + (s a t (φ + 2 * π / 3) 2) ^ 3
      = (s a t φ 0) ^ 3 + (s a t φ 1) ^ 3 + (s a t φ 2) ^ 3 := by
  obtain ⟨h0, h1, h2⟩ := s_shift a t φ; rw [h0, h1, h2]; ring

/-- The genuine phase invariant `cos 3φ` is unchanged by the `Z₃` shift and by conjugation:
    `cos(3(φ+2π/3)) = cos 3φ` and `cos(3(−φ)) = cos 3φ`. -/
theorem cos3phi_invariant (φ : ℝ) :
    Real.cos (3 * (φ + 2 * π / 3)) = Real.cos (3 * φ)
    ∧ Real.cos (3 * (-φ)) = Real.cos (3 * φ) := by
  constructor
  · rw [show 3 * (φ + 2 * π / 3) = 3 * φ + 2 * π by ring, Real.cos_add_two_pi]
  · rw [show 3 * (-φ) = -(3 * φ) by ring, Real.cos_neg]

/-! ## The punchline: `φ = 2/9` is not unique -/

/-- **`φ = 2/9` and `φ = 2/9 + 2π/3` give identical charged-lepton masses.**  (The first two
    power sums are phase-independent (`sum_s`, `sum_s_sq`); the third agrees by
    `moments_Z3_invariant`.)  The naive PDG fit returns the *second* value (`≈ 2.3166`), so
    privileging "2/9" is a convention.  The physical content is `cos 3φ = cos(2/3)`. -/
theorem phase_2_9_not_unique (a t : ℝ) :
    (s a t (2 / 9 + 2 * π / 3) 0) ^ 3 + (s a t (2 / 9 + 2 * π / 3) 1) ^ 3
        + (s a t (2 / 9 + 2 * π / 3) 2) ^ 3
      = (s a t (2 / 9) 0) ^ 3 + (s a t (2 / 9) 1) ^ 3 + (s a t (2 / 9) 2) ^ 3 :=
  moments_Z3_invariant a t (2 / 9)

/-- The invariant value the masses actually fix is `cos(3·(2/9)) = cos(2/3)` — and the shifted
    representative gives the *same* `cos`.  So a structural derivation must explain `cos(2/3)`
    (about 0.7859, not a "nice" number), not the branch- and unit-dependent rational `2/9`. -/
theorem invariant_is_cos_two_thirds :
    Real.cos (3 * (2 / 9 : ℝ)) = Real.cos (2 / 3)
    ∧ Real.cos (3 * (2 / 9 + 2 * π / 3 : ℝ)) = Real.cos (2 / 3) := by
  refine ⟨by rw [show (3 : ℝ) * (2 / 9) = 2 / 3 by norm_num], ?_⟩
  rw [show 3 * (2 / 9 + 2 * π / 3 : ℝ) = 2 / 3 + 2 * π by ring, Real.cos_add_two_pi]

/-! ## A natural deep identity (Tier 3.1) is FALSIFIED

The most natural "mass phase = gauge angle" guess is `cos 3φ = cos²θ_W = 7/9`.  But the genuine
invariant `cos(2/3)` is strictly greater than `7/9`, so that identity is empirically wrong (the
gap ≈ 0.008 is ~100× the data precision). -/

/-- **`cos(2/3) > 7/9`** (strict), via `sin(1/3) < 1/3`:
    `cos(2/3) = 2cos²(1/3) − 1 = 1 − 2sin²(1/3) > 1 − 2/9 = 7/9`. -/
theorem cos_two_thirds_gt_seven_ninths : Real.cos (2 / 3) > 7 / 9 := by
  rw [show (2 : ℝ) / 3 = 2 * (1 / 3) by norm_num, Real.cos_two_mul]
  have hs : Real.sin (1 / 3) < 1 / 3 := Real.sin_lt (by norm_num)
  have hs0 : 0 < Real.sin (1 / 3) :=
    Real.sin_pos_of_pos_of_lt_pi (by norm_num) (by linarith [Real.pi_gt_three])
  have hpyth : Real.sin (1 / 3) ^ 2 + Real.cos (1 / 3) ^ 2 = 1 := Real.sin_sq_add_cos_sq (1 / 3)
  nlinarith [hs, hs0, hpyth]

/-- **The natural Tier-3.1 deep identity is FALSE.**  The phase invariant `cos 3φ = cos(2/3)`
    is NOT the Weinberg `cos²θ_W = 7/9` (`cos(2/3) > 7/9`).  So "the lepton mass-generation
    third-harmonic equals the gauge mixing angle" — the obvious way the two `2/9`'s could be
    one object — is ruled out.  The two `2/9`'s coincide only as the *values* `sin²θ_W` and
    `Q/3`, not through their cosines. -/
theorem phase_invariant_ne_cos_sq_thetaW : Real.cos (3 * (2 / 9 : ℝ)) ≠ 7 / 9 := by
  rw [show (3 : ℝ) * (2 / 9) = 2 / 3 by norm_num]
  have := cos_two_thirds_gt_seven_ninths; linarith

end SCPv59.PhaseAmbiguity
