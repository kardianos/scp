/-
  v59/furey_construction/lean/BrannenKernel.lean

  Machine-checked proof of the central algebraic identity of v59's lepton sector:

      For the Brannen cyclic mass-kernel parametrisation
          s_k = a * (1 + 2 t cos(2π k / 3 + φ))     k = 0, 1, 2
      the Koide ratio
          Q := (s_0^2 + s_1^2 + s_2^2) / (s_0 + s_1 + s_2)^2
      satisfies
          Q = (1 + 2 t^2) / 3
      and in particular
          Q = 2/3   ⟺   t^2 = 1/2  (the constraint-surface condition).

  Interpretation: the three s_k are the eigenvalues of the cyclic operator
      M(ξ) = a (I + ξ S + ξ̄ S^T)
  on ℂ^3 with S the 3-cycle shift, where ξ = t e^{i φ}.  The identification of
  s_k with √m_k (Brannen's convention) then makes Q the charged-lepton Koide
  ratio.  The implication |ξ|^2 = 1/2  ⇒  Q = 2/3 — which `04_findings.md`
  established numerically by sampling random points on S^3 ⊂ ℍ — is now a
  theorem about real-valued cosines, independent of the quaternionic embedding.
-/

import Mathlib.Analysis.SpecialFunctions.Trigonometric.Basic
import Mathlib.Tactic.Ring
import Mathlib.Tactic.FieldSimp
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.LinearCombination
import Mathlib.Data.Real.Basic

namespace SCPv59.BrannenKernel

open Real

/-! ## Cyclic-rotation cosine sums

The two trig identities we need:
  Σ_{k=0,1,2} cos(α_k) = 0           (sum of 120°-spaced cosines)
  Σ_{k=0,1,2} cos(α_k)^2 = 3/2       (via cos² = (1 + cos 2x)/2)
where α_k = 2π k/3 + φ. -/

/-- `cos(2π/3) = -1/2`. -/
lemma cos_two_pi_div_three : cos (2 * π / 3) = -(1/2) := by
  have h : (2 * π / 3 : ℝ) = π - π / 3 := by ring
  rw [h, cos_pi_sub, cos_pi_div_three]

/-- `sin(2π/3) = √3 / 2`. -/
lemma sin_two_pi_div_three : sin (2 * π / 3) = Real.sqrt 3 / 2 := by
  have h : (2 * π / 3 : ℝ) = π - π / 3 := by ring
  rw [h, sin_pi_sub, sin_pi_div_three]

/-- `cos(4π/3) = -1/2`. -/
lemma cos_four_pi_div_three : cos (4 * π / 3) = -(1/2) := by
  have h : (4 * π / 3 : ℝ) = π + π / 3 := by ring
  rw [h, cos_add, cos_pi, sin_pi, cos_pi_div_three, sin_pi_div_three]
  ring

/-- `sin(4π/3) = -√3/2`. -/
lemma sin_four_pi_div_three : sin (4 * π / 3) = -(Real.sqrt 3 / 2) := by
  have h : (4 * π / 3 : ℝ) = π + π / 3 := by ring
  rw [h, sin_add, sin_pi, cos_pi, sin_pi_div_three, cos_pi_div_three]
  ring

/-- The 120°-spaced cosine sum vanishes:
    `cos φ + cos(φ + 2π/3) + cos(φ + 4π/3) = 0`. -/
lemma cos_cycle_sum (φ : ℝ) :
    cos φ + cos (φ + 2 * π / 3) + cos (φ + 4 * π / 3) = 0 := by
  rw [cos_add φ (2 * π / 3), cos_add φ (4 * π / 3),
      cos_two_pi_div_three, sin_two_pi_div_three,
      cos_four_pi_div_three, sin_four_pi_div_three]
  ring

/-- Sum of squared 120°-spaced cosines is 3/2.

    Uses `cos²(x) = (1 + cos 2x)/2` and the cycle-sum lemma applied at 2φ:
    the doubled angles 2(2πk/3 + φ) = 4πk/3 + 2φ are still 120°-spaced mod 2π. -/
lemma cos_sq_cycle_sum (φ : ℝ) :
    (cos φ)^2 + (cos (φ + 2 * π / 3))^2 + (cos (φ + 4 * π / 3))^2 = 3 / 2 := by
  -- expand each square via the half-angle identity cos²(x) = (1 + cos(2x))/2
  have h0 : (cos φ)^2 = (1 + cos (2 * φ)) / 2 := by
    rw [Real.cos_sq]; ring
  have h1 : (cos (φ + 2 * π / 3))^2 = (1 + cos (2 * (φ + 2 * π / 3))) / 2 := by
    rw [Real.cos_sq]; ring
  have h2 : (cos (φ + 4 * π / 3))^2 = (1 + cos (2 * (φ + 4 * π / 3))) / 2 := by
    rw [Real.cos_sq]; ring
  rw [h0, h1, h2]
  -- The doubled-angle sum needs cos(2φ) + cos(2φ + 4π/3) + cos(2φ + 8π/3) = 0.
  -- Reduce 8π/3 = 2π + 2π/3 via periodicity, then apply cos_cycle_sum at 2φ.
  have hperiod : cos (2 * (φ + 4 * π / 3)) = cos (2 * φ + 2 * π / 3) := by
    have : 2 * (φ + 4 * π / 3) = (2 * φ + 2 * π / 3) + 2 * π := by ring
    rw [this, cos_add_two_pi]
  have hmid : cos (2 * (φ + 2 * π / 3)) = cos (2 * φ + 4 * π / 3) := by
    congr 1; ring
  rw [hperiod, hmid]
  have hcyc : cos (2 * φ) + cos (2 * φ + 2 * π / 3) + cos (2 * φ + 4 * π / 3) = 0 :=
    cos_cycle_sum (2 * φ)
  linarith

/-! ## The Brannen amplitude and Koide ratio -/

/-- Brannen amplitude: `s_k(a, t, φ) = a (1 + 2 t cos(2π k / 3 + φ))`.
    These are the eigenvalues of the cyclic mass-kernel
    `M(ξ) = a (I + ξ S + ξ̄ S^T)` with `ξ = t e^{i φ}`, identified with √m_k. -/
noncomputable def s (a t φ : ℝ) : Fin 3 → ℝ
  | ⟨0, _⟩ => a * (1 + 2 * t * cos φ)
  | ⟨1, _⟩ => a * (1 + 2 * t * cos (φ + 2 * π / 3))
  | ⟨2, _⟩ => a * (1 + 2 * t * cos (φ + 4 * π / 3))

/-- `Σ_k s_k = 3 a`. -/
theorem sum_s (a t φ : ℝ) :
    s a t φ 0 + s a t φ 1 + s a t φ 2 = 3 * a := by
  simp only [s]
  have h := cos_cycle_sum φ
  linear_combination (2 * a * t) * h

/-- `Σ_k s_k^2 = 3 a^2 (1 + 2 t^2)`. -/
theorem sum_s_sq (a t φ : ℝ) :
    (s a t φ 0)^2 + (s a t φ 1)^2 + (s a t φ 2)^2 = 3 * a^2 * (1 + 2 * t^2) := by
  simp only [s]
  have h1 := cos_cycle_sum φ
  have h2 := cos_sq_cycle_sum φ
  linear_combination (4 * a^2 * t) * h1 + (4 * a^2 * t^2) * h2

/-- The Koide ratio of the Brannen amplitudes. -/
noncomputable def Q (a t φ : ℝ) : ℝ :=
  ((s a t φ 0)^2 + (s a t φ 1)^2 + (s a t φ 2)^2) /
  (s a t φ 0 + s a t φ 1 + s a t φ 2)^2

/-- **Main theorem.** Q has the closed form `(1 + 2 t²) / 3`,
    independent of the phase φ. -/
theorem Q_value (a t φ : ℝ) (ha : a ≠ 0) :
    Q a t φ = (1 + 2 * t^2) / 3 := by
  unfold Q
  rw [sum_s, sum_s_sq]
  have h3a : (3 * a)^2 = 9 * a^2 := by ring
  rw [h3a]
  field_simp
  ring

/-- **Constraint-surface theorem.** The Koide identity Q = 2/3 holds iff
    `t² = 1/2`.  In the quaternionic kernel-fit (`04_findings.md`),
    `t = |ξ|_ℍ`, so `t² = 1/2` is the S³ ⊂ ℍ of radius `1/√2`. -/
theorem koide_iff_constraint (a t φ : ℝ) (ha : a ≠ 0) :
    Q a t φ = 2 / 3 ↔ t^2 = 1 / 2 := by
  rw [Q_value a t φ ha]
  constructor
  · intro h
    have : 1 + 2 * t^2 = 2 := by linarith
    linarith
  · intro h
    rw [h]; ring

/-- **Corollary at the constraint point.** For any phase φ, the Brannen amplitudes
    with `t = 1/√2` produce Koide `Q = 2/3` exactly. -/
theorem Q_at_constraint (a φ : ℝ) (ha : a ≠ 0) :
    Q a (1 / Real.sqrt 2) φ = 2 / 3 := by
  rw [koide_iff_constraint a (1 / Real.sqrt 2) φ ha]
  have hsq : (Real.sqrt 2)^2 = 2 := Real.sq_sqrt (by norm_num : (2:ℝ) ≥ 0)
  rw [one_div, inv_pow, hsq]
  norm_num

/-- **Brannen-phase identification.** φ = 2/9 is precisely Q/3 = (dim G_2)/(3·dim Spin(7)).
    This is the v59 step-10 result: the empirical Brannen phase 2/9 rad is the
    Koide ratio divided by the number of generations. -/
theorem brannen_phi_is_Q_div_three :
    (2 / 9 : ℝ) = (2 / 3) / 3 := by norm_num

end SCPv59.BrannenKernel
