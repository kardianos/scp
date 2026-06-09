/-
  v62/lean/PhaseNoGo.lean

  THE BRANNEN-PHASE NO-GO  (machine-checked core of v62).

  This file is the formal counterpart of ../no_go/NOGO.md and
  ../no_go/transcendence_nogo.py.  It rewrites v59's `LeptonPhaseEmpirical.lean`:
  that file ASSERTED, in prose, that the phase invariant `cos 3φ = cos(2/3)` is
  "transcendental, non-geometric".  Here that prose becomes a THEOREM, and we
  draw the consequence that closes the spectral-projection program for the phase.

  Number-type discriminant (../THESIS.md):
    * Every output of an ALGEBRAIC construction -- Casimir eigenvalues, group
      characters, dimensions, eigenvalues (hence eigenphase cosines) of a matrix
      with algebraic entries -- is an ALGEBRAIC number.
    * The Brannen phase invariant `cos(3φ) = cos(2/3)` and `cos φ = cos(2/9)` are
      TRANSCENDENTAL (Lindemann-Weierstrass: `cos q` is transcendental for every
      nonzero rational `q`).
    * Hence the lepton Brannen phase `φ = 2/9` is NOT producible by any algebraic
      (representation / character / spectral / J-eigenphase) map.  It is a
      dynamical (loop-ratio, P3) object or an input -- never an algebra output.

  The single `sorry` is the Lindemann-Weierstrass input.  Mathlib currently ships
  only its analytical part
  (`Mathlib/NumberTheory/Transcendental/Lindemann/AnalyticalPart.lean`), so the
  finished statement is cited classically, exactly as v59's Lean culture allows a
  cited `sorry` on a deep external theorem.  Everything downstream is proven.

  Build (against the v59 Mathlib, as in v61/lean/EwVevHome.lean):
    cd v59/furey_construction/lean && lake env lean ../../../v62/lean/PhaseNoGo.lean
-/
import Mathlib

namespace SCPv62.PhaseNoGo

open Real

/-! ## 1. The Lindemann-Weierstrass input (the only `sorry`) -/

/-- **Lindemann-Weierstrass (corollary).**  For every nonzero rational `q`, the
    real number `cos q` is transcendental over `ℚ`.

    Proof (classical): if `α ≠ 0` is algebraic then `e^{α}` is transcendental
    (Lindemann-Weierstrass).  Taking `α = i q` (algebraic, nonzero), `e^{iq}` is
    transcendental; were `cos q` algebraic, `e^{iq}` would be a root of
    `z² - 2(cos q) z + 1` with algebraic coefficients, hence algebraic -- a
    contradiction.  Mathlib has the analytical part of L-W
    (`…/Transcendental/Lindemann/AnalyticalPart`) but not yet the packaged
    theorem, so this corollary is taken as a cited input. -/
theorem cos_rat_transcendental (q : ℚ) (hq : q ≠ 0) :
    Transcendental ℚ (Real.cos (q : ℝ)) := by
  sorry

/-! ## 2. The two transcendental phase quantities -/

/-- The Brannen-phase **invariant** `cos(3φ) = cos(2/3)` is transcendental. -/
theorem cos_two_thirds_transcendental : Transcendental ℚ (Real.cos (2 / 3)) := by
  have h := cos_rat_transcendental (2 / 3) (by norm_num)
  norm_num at h
  exact h

/-- `cos φ = cos(2/9)` (cosine of the phase itself) is transcendental. -/
theorem cos_two_ninths_transcendental : Transcendental ℚ (Real.cos (2 / 9)) := by
  have h := cos_rat_transcendental (2 / 9) (by norm_num)
  norm_num at h
  exact h

/-! ## 3. The no-go consequences -/

/-- **No algebraic invariant equals the phase invariant.**  Any Casimir
    eigenvalue, dimension, character value, structure constant, or eigenphase
    cosine of a matrix with algebraic entries is an algebraic real; none of them
    can equal `cos(2/3)`.  (`Transcendental ℚ x` is by definition `¬ IsAlgebraic ℚ x`.) -/
theorem no_algebraic_eq_phase_invariant {x : ℝ} (hx : IsAlgebraic ℚ x) :
    x ≠ Real.cos (2 / 3) := by
  intro h
  exact absurd (h ▸ hx) cos_two_thirds_transcendental

/-- **No algebraic operator realizes the phase rotation.**  A real matrix with
    algebraic entries has an algebraic characteristic polynomial, hence algebraic
    eigenvalues; if `e^{iθ}` is such an eigenvalue then `cos θ = (e^{iθ}+e^{-iθ})/2`
    is algebraic.  Since `cos(2/9)` is transcendental, no such operator (no `J`,
    no Furey-pinned complex structure, no Shulga-modulated rotation) can have
    eigenphase `2/9`.  Formally: `cos(2/9)` is not algebraic. -/
theorem no_algebraic_eigenphase_two_ninths : ¬ IsAlgebraic ℚ (Real.cos (2 / 9)) :=
  cos_two_ninths_transcendental

/-- **Headline.**  Both phase quantities are transcendental, so the lepton
    Brannen phase `φ = 2/9` is closed to every algebraic construction. -/
theorem brannen_phase_no_go :
    Transcendental ℚ (Real.cos (2 / 3)) ∧ Transcendental ℚ (Real.cos (2 / 9)) :=
  ⟨cos_two_thirds_transcendental, cos_two_ninths_transcendental⟩

/-! ## 4. The POSITIVE side (carried from v59 LeptonPhaseEmpirical, self-contained):
       `φ = 2/9` is a tight, real empirical regularity = `Koide Q / 3`.
       The no-go says only that its MECHANISM is not algebraic, not that it is
       unreal. -/

/-- Measured charged-lepton Brannen phase (PDG-2024, principal `Z₃` branch). -/
def phi_measured : ℚ := 22222963 / 100000000

/-- Measured Koide ratio (PDG-2024). -/
def Q_measured : ℚ := 66666051 / 100000000

/-- `φ = 2/9` holds to `< 10⁻⁵` (Koide-class tightness). -/
theorem phase_agrees_2_9 : |(2 / 9 : ℚ) - phi_measured| < 1 / 100000 := by
  unfold phi_measured; rw [abs_sub_comm, abs_of_nonneg (by norm_num)]; norm_num

/-- `Q = 2/3` holds to the same `< 10⁻⁵`. -/
theorem koide_agrees_2_3 : |(2 / 3 : ℚ) - Q_measured| < 1 / 100000 := by
  unfold Q_measured; rw [abs_of_nonneg (by norm_num)]; norm_num

/-- The exact rational relation behind the empirics: `(2/3)/3 = 2/9`, i.e.
    `φ = Q / (#generations)`.  This is the relation whose value (`cos` of it) the
    no-go shows to be transcendental. -/
theorem phase_eq_Q_over_three : (2 / 3 : ℚ) / 3 = 2 / 9 := by norm_num

/-- The full phase story in one statement: the empirical regularity is tight and
    rational (`φ = Q/3 = 2/9`), yet its dynamical invariant `cos(3φ) = cos(2/3)`
    is transcendental -- so the regularity is real but its origin is not algebraic. -/
theorem phase_real_but_not_algebraic :
    |(2 / 9 : ℚ) - phi_measured| < 1 / 100000
    ∧ (2 / 3 : ℚ) / 3 = 2 / 9
    ∧ Transcendental ℚ (Real.cos (2 / 3)) :=
  ⟨phase_agrees_2_9, phase_eq_Q_over_three, cos_two_thirds_transcendental⟩

end SCPv62.PhaseNoGo
