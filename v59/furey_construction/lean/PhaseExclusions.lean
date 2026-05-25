/-
  v59/furey_construction/lean/PhaseExclusions.lean

  **Negative results for the φ = 2/9 program** — proving candidate mechanisms FALSE, to narrow
  the design space (the project's null-result methodology).  Companion to `BrannenPhase`
  (`z3_potential_does_not_select_2_9`).

  Established exclusions:

    (E1) φ = 2/9 is rational *in radians*, hence **not a rational multiple of π** — so it is NOT
         a geometric/holonomy/root angle (those are π-rational).  The phase must be a *ratio*
         (like Q/3) inserted as a radian value, not an angle read off a rotation.  This rules out
         the whole "tie a G₂⊂Spin(7) geometric angle to Q" class of mechanisms.

    (E2) Corollary: φ ≠ π/n for every n (no "nice fraction of π").

    (E3) The lowest Z₃-invariant harmonic potentials (`cos 3φ`, `cos 6φ`) have no critical point
         at 2/9 — naive vacuum alignment fails at both harmonics (extends `BrannenPhase`).

    (E4) Given cW = 5 (dual Coxeter), the gauge value `sin²θ_W = 2/9` **pins `cBL = 2` uniquely**
         (`cBL/(5+2cBL) = 2/9 ⟺ cBL = 2`).  So the open "2" of Tier 2.2 is not free — it is
         forced by 2/9; other small integers `1,3,4` give `1/7, 3/11, 4/13 ≠ 2/9`.
-/
import Mathlib

namespace SCPv59.PhaseExclusions

/-! ## (E1) The phase is rational in radians ⇒ not a rational multiple of π -/

/-- **φ = 2/9 is not a rational multiple of π.**  If `2/9 = q·π` with `q ∈ ℚ` then `π` would be
    rational, contradicting `irrational_pi`.  So the Brannen phase is *not* a geometric angle
    (holonomy / rotation / root angle), which are all π-rational. -/
theorem phase_not_pi_rational : ¬ ∃ q : ℚ, ((2 : ℝ) / 9) = (q : ℝ) * Real.pi := by
  rintro ⟨q, hq⟩
  have hq0 : q ≠ 0 := by rintro rfl; simp at hq
  have hqR : (q : ℝ) ≠ 0 := by exact_mod_cast hq0
  exact irrational_pi ⟨(2 / 9 : ℚ) / q, by
    push_cast
    rw [div_eq_iff hqR, mul_comm Real.pi (q : ℝ)]
    exact hq⟩

/-- **(E1′) The holonomy target `Q = 2/3` is also not a rational multiple of π.**  Same proof
    as `phase_not_pi_rational`.  Consequence: the closed generation-cycle flux `3φ = Q` cannot
    be a *structural/geometric* holonomy (discrete-group, root, or associator holonomies are
    rational multiples of π).  Together with the fact that the v58 connection `ω ∝ ∇log(MM̃)` is
    an exact 1-form (zero holonomy on any loop), this **rules out the geometric-phase mechanism**
    for `φ = Q/3`: the flux is a free Aharonov-Bohm phase, unforced by the algebra. -/
theorem koide_not_pi_rational : ¬ ∃ q : ℚ, ((2 : ℝ) / 3) = (q : ℝ) * Real.pi := by
  rintro ⟨q, hq⟩
  have hq0 : q ≠ 0 := by rintro rfl; simp at hq
  have hqR : (q : ℝ) ≠ 0 := by exact_mod_cast hq0
  exact irrational_pi ⟨(2 / 3 : ℚ) / q, by
    push_cast; rw [div_eq_iff hqR, mul_comm Real.pi (q : ℝ)]; exact hq⟩

/-- **(E2) Corollary**: `φ ≠ π/n` for every natural `n` — the phase is no "nice fraction of π". -/
theorem phase_ne_pi_div_nat (n : ℕ) : (Real.pi / n : ℝ) ≠ 2 / 9 := by
  rcases Nat.eq_zero_or_pos n with hn | hn
  · subst hn; rw [Nat.cast_zero, div_zero]; norm_num
  · intro h
    exact phase_not_pi_rational ⟨1 / n, by push_cast; rw [← h]; ring⟩

/-! ## (E3) The low Z₃-harmonic potentials do not select 2/9 -/

/-- **The `cos 6φ` potential also fails** (the second Z₃-invariant harmonic): its derivative
    `−6 sin(6φ)` is nonzero at `2/9` (since `0 < 4/3 < π`).  With
    `BrannenPhase.z3_potential_does_not_select_2_9` (the `cos 3φ` case), neither of the two lowest
    Z₃-invariant harmonics has an extremum at the Brannen phase. -/
theorem cos6_potential_does_not_select_2_9 :
    deriv (fun φ : ℝ => Real.cos (6 * φ)) (2 / 9) ≠ 0 := by
  have hf : HasDerivAt (fun φ : ℝ => 6 * φ) 6 (2 / 9) := by
    simpa using (hasDerivAt_id (2 / 9 : ℝ)).const_mul 6
  have hc : HasDerivAt (fun φ : ℝ => Real.cos (6 * φ)) (-Real.sin (6 * (2 / 9)) * 6) (2 / 9) :=
    hf.cos
  rw [hc.deriv]
  have hsin : Real.sin (6 * (2 / 9)) > 0 := by
    rw [show (6 : ℝ) * (2 / 9) = 4 / 3 by norm_num]
    exact Real.sin_pos_of_pos_of_lt_pi (by norm_num) (by linarith [Real.pi_gt_three])
  nlinarith [hsin]

/-! ## (E4) The gauge "2" is pinned by 2/9 -/

/-- The Pati-Salam mixing angle as a function of `cBL` (with `cW = 5`), simplified:
    `sin²θ_W = cBL/(5 + 2·cBL)`. -/
def sin2_of_cBL (cBL : ℚ) : ℚ := cBL / (5 + 2 * cBL)

/-- This matches `WeinbergPatiSalam` at `cBL = 2`: `2/(5+4) = 2/9`. -/
theorem sin2_of_cBL_two : sin2_of_cBL 2 = 2 / 9 := by unfold sin2_of_cBL; norm_num

/-- **(E4) The observed `2/9` forces `cBL = 2` uniquely** (given `cW = 5`):
    `cBL/(5+2cBL) = 2/9 ⟺ cBL = 2`.  So the open Tier-2.2 "2" is *not* a free parameter — it is
    determined by the value `2/9` and the dual-Coxeter `5`. -/
theorem gauge_cBL_pinned (cBL : ℚ) (h : 5 + 2 * cBL ≠ 0) :
    sin2_of_cBL cBL = 2 / 9 ↔ cBL = 2 := by
  unfold sin2_of_cBL
  rw [div_eq_div_iff h (by norm_num : (9 : ℚ) ≠ 0)]
  constructor
  · intro hcbl; linarith
  · intro hcbl; rw [hcbl]; ring

/-- The other small-integer candidates give the *wrong* angle (not `2/9`). -/
theorem gauge_small_ints_fail :
    sin2_of_cBL 1 = 1 / 7 ∧ sin2_of_cBL 3 = 3 / 11 ∧ sin2_of_cBL 4 = 4 / 13 := by
  refine ⟨?_, ?_, ?_⟩ <;> (unfold sin2_of_cBL; norm_num)

/-! ## (E5) Mass and gauge `2/9` come from different integers

The mass route reduces `14/63` by the common factor `7 = dimImO`; the gauge route reduces
`10/45` by `5 = h∨(Spin7)`.  The two `2/9`'s have *different* pre-reduction integers, so any
"single shared integer drives both" hypothesis is constrained (evidence the Tier-3.1 coincidence
is not a trivial identity). -/
theorem mass_gauge_distinct_reductions :
    (14 : ℚ) / 63 = 2 / 9 ∧ (10 : ℚ) / 45 = 2 / 9 ∧ (14 / 63 : ℚ) = 10 / 45 ∧ (14 ≠ 10) := by
  refine ⟨by norm_num, by norm_num, by norm_num, by norm_num⟩

end SCPv59.PhaseExclusions
