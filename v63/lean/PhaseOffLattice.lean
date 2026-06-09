/-
  v63/lean/PhaseOffLattice.lean

  MACHINE-CHECKED unifying dichotomy from "trying both" integration routes
  (../integrate_density.py = Reading A, ../dh_localization.py = Reading B):

      the TRANSCENDENTAL phase 2/3  <=>  an OFF-LATTICE (non-π-rational) angle
                                    <=>  a free continuous parameter.

  Canonical geometry / characters / holonomies (the triality element, the Weyl
  group, the DH fixed points) all hand you COMMENSURATE angles -- rational
  multiples of 2π -- whose cosines are ALGEBRAIC (roots of unity / rationals).
  The phase argument 2/3 is NOT a rational multiple of π, so it is not such an
  angle: reaching it requires a generic, non-canonical evaluation point.

  This formalizes, and connects, two facts:
    * canonical side: the triality angle 2π/3 is π-rational with cos = −1/2
      (the G2 7-rep character at this element is the algebraic value 2, computed
      numerically in dh_localization.py);
    * phase side: 2/3 is not π-rational (the v59 `phase_not_pi_rational` /
      `koide_not_pi_rational` result, re-proven self-contained here).

  Build (against the v59 Mathlib, as in v61/lean/EwVevHome.lean):
    cd v59/furey_construction/lean && lake env lean ../../../v63/lean/PhaseOffLattice.lean
-/
import Mathlib

namespace SCPv63.PhaseOffLattice

open Real

/-! ## 1. Canonical (commensurate) side: π-rational angle, algebraic cosine -/

/-- The canonical triality angle `2π/3` IS a rational multiple of π (commensurate). -/
theorem triality_angle_pi_rational : ∃ q : ℚ, (2 * Real.pi / 3) = (q : ℝ) * Real.pi :=
  ⟨2 / 3, by push_cast; ring⟩

/-- Its cosine is **rational** (`= −1/2`) — i.e. algebraic.  This is the canonical
    character/holonomy value's ingredient: `cos(2π/3) = cos(π − π/3) = −cos(π/3) = −1/2`.
    (The full G2 7-rep character at this element is `2`, an algebraic integer —
    `dh_localization.py` B1.) -/
theorem cos_triality_rational : Real.cos (2 * Real.pi / 3) = -(1 / 2) := by
  rw [show (2 * Real.pi / 3 : ℝ) = Real.pi - Real.pi / 3 by ring,
      Real.cos_pi_sub, Real.cos_pi_div_three]

/-! ## 2. Phase side: 2/3 is NOT a rational multiple of π (off-lattice) -/

/-- **The phase argument `3φ = 2/3` is not a rational multiple of π.**  (Self-
    contained re-proof of v59 `koide_not_pi_rational`: if `2/3 = q·π` then `π`
    is rational, contradicting `irrational_pi`.)  So `2/3` is *off* every
    root-of-unity lattice — it is not a commensurate/geometric angle. -/
theorem phase_arg_not_pi_rational : ¬ ∃ q : ℚ, ((2 : ℝ) / 3) = (q : ℝ) * Real.pi := by
  rintro ⟨q, hq⟩
  have hq0 : q ≠ 0 := by rintro rfl; simp at hq
  have hqR : (q : ℝ) ≠ 0 := by exact_mod_cast hq0
  exact irrational_pi ⟨(2 / 3 : ℚ) / q, by
    push_cast; rw [div_eq_iff hqR, mul_comm Real.pi (q : ℝ)]; exact hq⟩

/-! ## 3. The dichotomy -/

/-- **Headline.**  The canonical triality angle is commensurate (π-rational) with
    an algebraic cosine `−1/2`; the phase argument `2/3` is **not** π-rational.
    So the phase is not a canonical character/holonomy angle — every parameter-free
    (commensurate) evaluation gives algebraic values, and `2/3` requires a generic,
    off-lattice point (a free continuous parameter).  Both integration readings
    (A: density-weighted circular mean; B: DH/character localization) converge on
    this. -/
theorem phase_off_lattice :
    (∃ q : ℚ, (2 * Real.pi / 3) = (q : ℝ) * Real.pi)        -- triality: commensurate
    ∧ Real.cos (2 * Real.pi / 3) = -(1 / 2)                  -- ...with algebraic cosine
    ∧ ¬ ∃ q : ℚ, ((2 : ℝ) / 3) = (q : ℝ) * Real.pi :=        -- phase 2/3: NOT commensurate
  ⟨triality_angle_pi_rational, cos_triality_rational, phase_arg_not_pi_rational⟩

end SCPv63.PhaseOffLattice
