/-
  v59/furey_construction/lean/BrannenPhase.lean

  **Tier 0 of the φ = 2/9 program** (see `7D_Algebra/notes/2026-05-24-phase-2-9-task-scope.md`):
  de-circularize the Brannen phase and locate the gap precisely.

  In `LieDimensions.lean`, `brannen_phi := koide_Q / 3` is a *definition*, so the "φ = 2/9"
  theorems are arithmetic on an assumption.  Here we make the genuine content explicit:

    * `Q_phase_independent` — the Koide ratio is the SAME for every phase φ.  So Koide fixes
      the amplitude `t` (`t²=1/2 ⟺ Q=2/3`) but says **nothing** about φ: the relation φ=Q/3
      is independent content, not a consequence of the Koide identity.

    * The power-sum moments of the Brannen amplitudes:
          p₁ = Σ sₖ   = 3a                          (φ-independent)
          p₂ = Σ sₖ²  = 3a²(1+2t²)                   (φ-independent)
          p₃ = Σ sₖ³  = 3a³(1+6t²+2t³·cos 3φ)        (φ only via cos 3φ)
      The phase enters observables **only through `cos 3φ`** — the three generations (the Z₃
      120°-spacing) force the third harmonic.  So the natural phase variable is `3φ`, and the
      v59 conjecture `φ = Q/3` is exactly **`3φ = Q`**: the third-harmonic angle equals the
      Koide ratio.  *That* identity (the factor 3 = #generations) is the open Tier-1 target.
-/
import Mathlib
import BrannenKernel
import LieDimensions

namespace SCPv59.BrannenPhase

open Real SCPv59.BrannenKernel

/-! ## The phase is unconstrained by Koide -/

/-- **The Koide ratio is phase-independent.**  `Q(a,t,φ) = Q(a,t,φ')` for all phases — the
    Koide identity fixes only the amplitude `t`, never φ.  Hence `φ = Q/3` is genuine extra
    content (the Tier-1 target), not a corollary of `Q = 2/3`. -/
theorem Q_phase_independent (a t φ φ' : ℝ) (ha : a ≠ 0) : Q a t φ = Q a t φ' := by
  rw [Q_value a t φ ha, Q_value a t φ' ha]

/-! ## The third harmonic: φ enters observables only via `cos 3φ` -/

/-- Sum of cubed 120°-spaced cosines is `(3/4)·cos 3φ` (via `cos³x = (cos 3x + 3 cos x)/4`
    and the doubled/tripled-angle periodicity: `3(φ+2πk/3) = 3φ + 2πk`). -/
lemma cos_cube_cycle_sum (φ : ℝ) :
    (cos φ)^3 + (cos (φ + 2 * π / 3))^3 + (cos (φ + 4 * π / 3))^3 = 3 / 4 * cos (3 * φ) := by
  have e : ∀ x : ℝ, (cos x)^3 = (cos (3 * x) + 3 * cos x) / 4 := by
    intro x; rw [Real.cos_three_mul]; ring
  rw [e φ, e (φ + 2 * π / 3), e (φ + 4 * π / 3)]
  have hc := cos_cycle_sum φ
  have h1 : cos (3 * (φ + 2 * π / 3)) = cos (3 * φ) := by
    rw [show 3 * (φ + 2 * π / 3) = 3 * φ + 2 * π by ring, Real.cos_add_two_pi]
  have h2 : cos (3 * (φ + 4 * π / 3)) = cos (3 * φ) := by
    rw [show 3 * (φ + 4 * π / 3) = (3 * φ + 2 * π) + 2 * π by ring, Real.cos_add_two_pi,
        Real.cos_add_two_pi]
  rw [h1, h2]; linarith

/-- **Third power sum of the Brannen amplitudes.**
    `Σ sₖ³ = 3a³(1 + 6t² + 2t³·cos 3φ)` — the phase appears only through `cos 3φ`. -/
theorem sum_s_cube (a t φ : ℝ) :
    (s a t φ 0)^3 + (s a t φ 1)^3 + (s a t φ 2)^3
      = 3 * a^3 * (1 + 6 * t^2 + 2 * t^3 * cos (3 * φ)) := by
  simp only [s]
  have h1 := cos_cycle_sum φ
  have h2 := cos_sq_cycle_sum φ
  have h3 := cos_cube_cycle_sum φ
  linear_combination (6 * a^3 * t) * h1 + (12 * a^3 * t^2) * h2 + (8 * a^3 * t^3) * h3

/-- **The three moments and what fixes the phase.**  Bundling `sum_s`, `sum_s_sq`, `sum_s_cube`:
    the first two power sums are phase-independent; the third introduces the phase exactly via
    `cos 3φ`.  So data `(p₁,p₂,p₃)` ↔ `(a, t, cos 3φ)`: Koide `Q = p₂/p₁²` fixes `t`, and the
    third moment fixes `cos 3φ` — the phase observable is the third harmonic `3φ`. -/
theorem moments (a t φ : ℝ) :
    (s a t φ 0 + s a t φ 1 + s a t φ 2 = 3 * a)
    ∧ ((s a t φ 0)^2 + (s a t φ 1)^2 + (s a t φ 2)^2 = 3 * a^2 * (1 + 2 * t^2))
    ∧ ((s a t φ 0)^3 + (s a t φ 1)^3 + (s a t φ 2)^3
        = 3 * a^3 * (1 + 6 * t^2 + 2 * t^3 * cos (3 * φ))) :=
  ⟨sum_s a t φ, sum_s_sq a t φ, sum_s_cube a t φ⟩

/-! ## The Tier-1 target, stated cleanly: `3φ = Q`

The third harmonic `3φ` is the only way the phase enters; the v59 conjecture `φ = Q/3` is the
statement that this third-harmonic angle equals the Koide ratio. -/

/-- The structural phase target as `3φ = Q`: with `Q = dimG₂/dimSpin7 = 2/3`, this is `3φ = 2/3`,
    i.e. `φ = 2/9`.  (Pure restatement; the *derivation* is Tier 1.) -/
theorem target_three_phi_eq_Q (φ : ℝ) (h : φ = 2 / 9) : 3 * φ = (2 / 3 : ℝ) := by
  rw [h]; norm_num

/-- The structural value, matching `LieDimensions.brannen_phi`: `Q/3 = (dimG₂/dimSpin7)/3 = 2/9`. -/
theorem phase_target_eq : ((SCPv59.koide_Q : ℚ) / 3 : ℚ) = 2 / 9 := by
  rw [SCPv59.koide_Q_value]; norm_num

/-- **Lepton-specificity (Tier 4, numeric).**  The relation `φ = Q/3` holds for the charged
    leptons (`φ_l ≈ 0.2222 = 2/9 = Q_l/3`) but **NOT** for the quarks: the fitted Brannen
    phases (`05_quark_sector.py`: `φ_d ≈ 0.1102`, `φ_u ≈ −2.02`) are far from
    `Q_d/3 = 11/45 ≈ 0.244` and `Q_u/3 = 23/81 ≈ 0.284`.  So whatever forces `3φ = Q` is
    lepton / color-singlet specific (consistent with the `J_c`-complex-line living on the
    lepton singlet `{0,7}`), not a universal Brannen relation. -/
theorem quark_phase_not_Q_div_three :
    (11 / 45 : ℚ) - 110176 / 1000000 > 1 / 10        -- down-type: |gap| > 0.1
    ∧ (23 / 81 : ℚ) - (-2020000 / 1000000) > 2 := by  -- up-type:   |gap| > 2
  constructor <;> norm_num

/-- **Tier 1 attempt T1.2 — the simplest Z₃ phase potential fails.**  A `Z₃`-symmetric phase
    potential must be built from `cos 3φ` (the lowest `Z₃`-invariant harmonic; cf. `sum_s_cube`).
    Its critical points are `V'(φ) = −3 sin(3φ) = 0 ⟺ 3φ ∈ πℤ`, i.e. `φ ∈ {0, π/3, …}` — never
    `2/9`.  Concretely `2/9` is NOT a critical point: `V'(2/9) ≠ 0` (since `sin(2/3) > 0`).  So
    the naive vacuum-alignment mechanism does not select the Brannen phase; a working mechanism
    must couple the phase to more structure (the open `J∘Z₃` route). -/
theorem z3_potential_does_not_select_2_9 :
    deriv (fun φ : ℝ => Real.cos (3 * φ)) (2 / 9) ≠ 0 := by
  have hf : HasDerivAt (fun φ : ℝ => 3 * φ) 3 (2 / 9) := by
    simpa using (hasDerivAt_id (2 / 9 : ℝ)).const_mul 3
  have hc : HasDerivAt (fun φ : ℝ => Real.cos (3 * φ)) (-Real.sin (3 * (2 / 9)) * 3) (2 / 9) :=
    hf.cos
  rw [hc.deriv]
  have hsin : Real.sin (3 * (2 / 9)) > 0 := by
    rw [show (3 : ℝ) * (2 / 9) = 2 / 3 by norm_num]
    exact Real.sin_pos_of_pos_of_lt_pi (by norm_num) (by linarith [Real.pi_gt_three])
  nlinarith [hsin]

/-- **OPEN (Tier 1).**  The dynamical content is `3φ_phys = Q`, where `φ_phys` is the phase that
    reproduces the charged-lepton mass ratios (it controls `cos 3φ` via `sum_s_cube`, the only
    phase-dependent moment).  `Q_phase_independent` shows this is *not* implied by `Q = 2/3`,
    and `quark_phase_not_Q_div_three` shows it is lepton-specific.  No proof is asserted here;
    candidate mechanisms (J∘Z₃ alignment on the singlet, phase potential, …) are in the
    task-scope note and need the generation↔ideal (Witt) map. -/
theorem phase_is_open_content : True := trivial

end SCPv59.BrannenPhase
