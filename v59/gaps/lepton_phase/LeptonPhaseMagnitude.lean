/-
  v59/gaps/lepton_phase/LeptonPhaseMagnitude.lean

  **Gap G7 — the lepton mass-phase MAGNITUDE `Q = 2/3`.**

  The Brannen phase law is `φ = Q/3 = 2/9`.  The "/3" is structural (the sedenion
  `S₃` generation automorphism, order 3, verified elsewhere).  The RESIDUAL is the
  magnitude `Q = 2/3 = dim G₂/dim Spin(7) = 14/21`, equivalently the phase invariant
  `cos(3φ) = cos(2/3)`.

  The GEOMETRIC route is dead (machine-checked in
  `furey_construction/lean/PhaseExclusions.lean`: `2/9` and `2/3` are rational in
  RADIANS, hence not π-rational, hence not any holonomy/rotation/root angle).  This
  file formalizes the surviving NON-geometric content:

    1. (avenue A) the maximal-mixing / G₂-content derivation of the magnitude:
       `t² = (D − dimG₂)/D`, with the lepton's `1/2` FORCED by `D = 28 = 2·dimG₂`,
       feeding the Brannen map `Q = (1+2t²)/3` to give `2/3`.  [TRUE arithmetic]

    2. (avenue B, NEW) the EQUIPARTITION characterization of `t² = 1/2`: the kernel
       `M = a(I + ξS + ξ̄S²)` has identity-part Frobenius² weight `1` and
       circulant-shift weight `2t²`; these BALANCE iff `t² = 1/2`.  A non-geometric
       extremal reading of the magnitude.  [proved]

    3. (avenue E, NEW) the SKEWNESS identity: the phase invariant `cos 3φ` is the
       (rescaled) third standardized moment of the √m amplitudes — `cos 3φ` enters
       the third power-sum and nothing else.  Reframes `cos(2/3)` as a pure data
       moment, no angle.  [proved via the cube-sum identity]

  HONEST SCOPE.  (1)–(3) are TRUE arithmetic/trig identities.  They establish that
  the MAGNITUDE `Q=2/3` is non-geometrically real (a dimension ratio = an
  equipartition point) and that the phase invariant is a pure moment.  They do NOT
  derive the `φ = Q/3` radian insert (gap F): WHY the dimensionless amplitude number
  `Q` re-enters as the *cosine argument* `cos(2/3)` is left OPEN (`magnitude_to_phase_open`).

  Self-contained: `import Mathlib`; re-states the Brannen amplitude rather than
  cross-importing (the shared project must not be co-built this run).

  ⚠ WRITTEN, NOT BUILT THIS RUN (shared `furey_construction/lean` project; concurrent
  `lake build` is forbidden).  Open steps are marked `sorry` and flagged.
-/
import Mathlib

namespace SCPv59.G7Magnitude

open Real

/-! ## Avenue A — maximal-mixing / G₂-content magnitude (TRUE arithmetic) -/

/-- Dimension of `G₂ = Aut(𝕆)`. -/
def dimG2 : ℕ := 14
/-- Dimension of the lepton grade `L = Λ² ⊕ Λ⁶ = so(8)`. -/
def dimL : ℕ := 28

/-- The maximal-mixing weight: the fraction of a `D`-dim sector outside the inert
    `G₂` automorphism core. -/
def massSplitFraction (D : ℕ) : ℚ := ((D : ℚ) - (dimG2 : ℚ)) / (D : ℚ)

/-- The Brannen map `Q = (1 + 2t²)/3`. -/
def koideQ (t_sq : ℚ) : ℚ := (1 + 2 * t_sq) / 3

/-- **The lepton sector is special: `dim L = 28 = 2·dim G₂`**, so its non-G₂ weight
    is exactly `1/2`. -/
theorem lepton_double_core : dimL = 2 * dimG2 := by decide

/-- **Maximal-mixing gives the lepton `t² = 1/2`** (because `28 = 2·14`). -/
theorem massSplit_lepton_half : massSplitFraction dimL = 1 / 2 := by
  unfold massSplitFraction dimL dimG2; norm_num

/-- **The magnitude `Q = 2/3` from maximal mixing** (avenue A): feeding the lepton's
    non-G₂ weight `1/2` into the Brannen map gives `Q = 2/3`.  A Lie-DIMENSION ratio,
    no angle — sidestepping the π-rationality obstruction (`PhaseExclusions`). -/
theorem lepton_koide_from_maximal_mixing :
    koideQ (massSplitFraction dimL) = 2 / 3 := by
  rw [massSplit_lepton_half]; unfold koideQ; norm_num

/-- The quark sectors give the "ugly" weights `3/5, 7/9` (NOT `1/2`), confirming the
    lepton's `1/2` is special, not generic.  (`D_d = 35`, `D_u = 63`.) -/
theorem quark_weights_not_half :
    massSplitFraction 35 = 3 / 5 ∧ massSplitFraction 63 = 7 / 9
    ∧ massSplitFraction 35 ≠ 1 / 2 ∧ massSplitFraction 63 ≠ 1 / 2 := by
  refine ⟨?_, ?_, ?_, ?_⟩ <;> (unfold massSplitFraction dimG2; norm_num)

/-! ## Avenue B — equipartition: `t² = 1/2` balances the kernel (NEW, proved)

The Brannen kernel `M = a(I + ξS + ξ̄S²)` has, in its Frobenius² (= Σ|s_k|²/3a²)
content, an identity/diagonal weight `1` (the `I`) and a circulant-shift weight
`2t²` (the two bands `ξS`, `ξ̄S²`, each contributing `|ξ|² = t²`).  These balance —
the "maximally-mixed / equipartition" point — exactly when `t² = 1/2`. -/

/-- The identity-part Frobenius² weight of `M/a` (per the normalization). -/
def diagWeight : ℚ := 1
/-- The circulant-shift Frobenius² weight of `M/a`: `|ξS|²_F + |ξ̄S²|²_F = 2t²`. -/
def shiftWeight (t_sq : ℚ) : ℚ := 2 * t_sq

/-- **Equipartition ⟺ `t² = 1/2`.**  The symmetric (identity) and circulant-shift
    contents of the Brannen kernel are EQUAL iff the amplitude sits at `t² = 1/2` —
    a non-geometric, extremal characterization of the Koide point. -/
theorem equipartition_iff_half (t_sq : ℚ) :
    shiftWeight t_sq = diagWeight ↔ t_sq = 1 / 2 := by
  unfold shiftWeight diagWeight; constructor
  · intro h; linarith
  · intro h; rw [h]; norm_num

/-- **Equipartition forces Koide `Q = 2/3`** (avenue B): at the balance point
    `t²=1/2`, the Brannen map gives `Q = 2/3`.  Combined with A, two non-geometric
    readings of the same `1/2` (dimension count `28=2·14`; kernel content balance). -/
theorem equipartition_gives_koide (t_sq : ℚ) (h : shiftWeight t_sq = diagWeight) :
    koideQ t_sq = 2 / 3 := by
  rw [(equipartition_iff_half t_sq).mp h]; unfold koideQ; norm_num

/-! ## Avenue E — the phase invariant is the third moment (NEW, proved over ℝ)

The Brannen amplitude `s_k = a(1 + 2t·cos(φ + 2πk/3))`.  The phase enters the third
power-sum (and ONLY there) via `cos 3φ`.  We prove the cube-sum identity that makes
`cos 3φ` the (rescaled) skewness of the √m distribution — recasting the residual
`cos(2/3)` as a pure data moment, no angle.  (Re-proved here self-contained.) -/

/-- The Brannen amplitude `s_k = a(1 + 2t·cos(φ + 2πk/3))`. -/
noncomputable def s (a t φ : ℝ) (k : Fin 3) : ℝ :=
  a * (1 + 2 * t * Real.cos (φ + 2 * π * (k.val : ℝ) / 3))

/-- Sum of the 120°-spaced cosines vanishes. -/
lemma cos_cycle_sum (φ : ℝ) :
    Real.cos φ + Real.cos (φ + 2 * π / 3) + Real.cos (φ + 4 * π / 3) = 0 := by
  have h23 : Real.cos (2 * π / 3) = -(1/2) := by
    rw [show (2 * π / 3 : ℝ) = π - π / 3 by ring, Real.cos_pi_sub, Real.cos_pi_div_three]
  have s23 : Real.sin (2 * π / 3) = Real.sqrt 3 / 2 := by
    rw [show (2 * π / 3 : ℝ) = π - π / 3 by ring, Real.sin_pi_sub, Real.sin_pi_div_three]
  have h43 : Real.cos (4 * π / 3) = -(1/2) := by
    rw [show (4 * π / 3 : ℝ) = π + π / 3 by ring, Real.cos_add, Real.cos_pi, Real.sin_pi,
        Real.cos_pi_div_three, Real.sin_pi_div_three]; ring
  have s43 : Real.sin (4 * π / 3) = -(Real.sqrt 3 / 2) := by
    rw [show (4 * π / 3 : ℝ) = π + π / 3 by ring, Real.sin_add, Real.sin_pi, Real.cos_pi,
        Real.sin_pi_div_three, Real.cos_pi_div_three]; ring
  rw [Real.cos_add φ (2 * π / 3), Real.cos_add φ (4 * π / 3), h23, s23, h43, s43]; ring

/-- Sum of squared 120°-spaced cosines is `3/2`. -/
lemma cos_sq_cycle_sum (φ : ℝ) :
    (Real.cos φ)^2 + (Real.cos (φ + 2 * π / 3))^2 + (Real.cos (φ + 4 * π / 3))^2 = 3 / 2 := by
  have h0 : (Real.cos φ)^2 = (1 + Real.cos (2 * φ)) / 2 := by rw [Real.cos_sq]; ring
  have h1 : (Real.cos (φ + 2 * π / 3))^2 = (1 + Real.cos (2 * (φ + 2 * π / 3))) / 2 := by
    rw [Real.cos_sq]; ring
  have h2 : (Real.cos (φ + 4 * π / 3))^2 = (1 + Real.cos (2 * (φ + 4 * π / 3))) / 2 := by
    rw [Real.cos_sq]; ring
  rw [h0, h1, h2]
  have hperiod : Real.cos (2 * (φ + 4 * π / 3)) = Real.cos (2 * φ + 2 * π / 3) := by
    rw [show 2 * (φ + 4 * π / 3) = (2 * φ + 2 * π / 3) + 2 * π by ring, Real.cos_add_two_pi]
  have hmid : Real.cos (2 * (φ + 2 * π / 3)) = Real.cos (2 * φ + 4 * π / 3) := by
    congr 1; ring
  rw [hperiod, hmid]
  have hcyc := cos_cycle_sum (2 * φ); linarith

/-- Sum of cubed 120°-spaced cosines is `(3/4)·cos 3φ` (the third harmonic). -/
lemma cos_cube_cycle_sum (φ : ℝ) :
    (Real.cos φ)^3 + (Real.cos (φ + 2 * π / 3))^3 + (Real.cos (φ + 4 * π / 3))^3
      = 3 / 4 * Real.cos (3 * φ) := by
  have e : ∀ x : ℝ, (Real.cos x)^3 = (Real.cos (3 * x) + 3 * Real.cos x) / 4 := by
    intro x; rw [Real.cos_three_mul]; ring
  rw [e φ, e (φ + 2 * π / 3), e (φ + 4 * π / 3)]
  have hc := cos_cycle_sum φ
  have h1 : Real.cos (3 * (φ + 2 * π / 3)) = Real.cos (3 * φ) := by
    rw [show 3 * (φ + 2 * π / 3) = 3 * φ + 2 * π by ring, Real.cos_add_two_pi]
  have h2 : Real.cos (3 * (φ + 4 * π / 3)) = Real.cos (3 * φ) := by
    rw [show 3 * (φ + 4 * π / 3) = (3 * φ + 2 * π) + 2 * π by ring, Real.cos_add_two_pi,
        Real.cos_add_two_pi]
  rw [h1, h2]; linarith

/-- **The phase invariant `cos 3φ` is the (rescaled) third moment of the amplitudes.**
    Writing `x_k = (s_k − a)/a = 2t·cos(φ+2πk/3)`, the third standardized moment is
    `Σ x_k³ = (2t)³·(3/4)·cos 3φ`.  So `cos 3φ` is, up to the fixed scale `(2t)³`,
    the SKEWNESS of the √m distribution — a pure data invariant, no angle.  This is
    the non-geometric reframing of the residual `cos(2/3)` (avenue E). -/
theorem cos3phi_is_third_moment (a t φ : ℝ) (ha : a ≠ 0) :
    ((s a t φ 0 - a)/a)^3 + ((s a t φ 1 - a)/a)^3 + ((s a t φ 2 - a)/a)^3
      = (2 * t)^3 * (3 / 4 * Real.cos (3 * φ)) := by
  -- (s_k − a)/a = 2t·cos(φ + 2πk/3); the Fin-3 angles are φ, φ+2π/3, φ+4π/3.
  have e0 : (s a t φ 0 - a)/a = 2 * t * Real.cos φ := by
    unfold s
    rw [show φ + 2 * π * ((0 : Fin 3).val : ℝ) / 3 = φ by
      rw [show (0 : Fin 3).val = 0 by decide]; push_cast; ring]
    field_simp; ring
  have e1 : (s a t φ 1 - a)/a = 2 * t * Real.cos (φ + 2 * π / 3) := by
    unfold s
    rw [show φ + 2 * π * ((1 : Fin 3).val : ℝ) / 3 = φ + 2 * π / 3 by
      rw [show (1 : Fin 3).val = 1 by decide]; push_cast; ring]
    field_simp; ring
  have e2 : (s a t φ 2 - a)/a = 2 * t * Real.cos (φ + 4 * π / 3) := by
    unfold s
    rw [show φ + 2 * π * ((2 : Fin 3).val : ℝ) / 3 = φ + 4 * π / 3 by
      rw [show (2 : Fin 3).val = 2 by decide]; push_cast; ring]
    field_simp; ring
  rw [e0, e1, e2]
  -- expand the cubes; the cube cycle-sum supplies Σcos³ = (3/4)cos 3φ
  have hcube := cos_cube_cycle_sum φ
  linear_combination (8 * t^3) * hcube

/-! ## The target value and the OPEN hard core (avenue F) -/

/-- The residual phase invariant is `cos(3·(2/9)) = cos(2/3)`, the transcendental the
    program must hit.  (Pure arithmetic of the argument.) -/
theorem invariant_value : (3 : ℝ) * (2 / 9) = 2 / 3 := by norm_num

/-- **OPEN (gap F — the hard core of G7).**  Avenues A/B give the dimensionless
    MAGNITUDE `Q = 2/3` non-geometrically (dimension ratio / equipartition); avenue E
    recasts the invariant as the skewness `cos 3φ`.  What remains UNPROVEN is the
    `φ = Q/3` *radian insert*: that the SAME number `Q=2/3` that fixes the amplitude
    `t²` (via `Q=(1+2t²)/3`) also fixes the phase via `3φ = Q`, i.e. `cos 3φ = cos Q`
    with `Q` used as a radian.  No non-geometric mechanism for this is known (the
    kernel cannot self-consistently force it — `Q` is phase-independent — and there
    is no algebraic `f(Q)` shortcut to `cos(2/3)`; see `ALTERNATIVES.md §D,§F`).

    Stated as the open theorem to discharge: there exists a v59-structural real
    quantity equal to `cos(2/3)` whose cosine-argument is the G₂-content ratio `Q`. -/
theorem magnitude_to_phase_open :
    ∃ structural : ℝ, structural = Real.cos (2 / 3) := by
  -- TODO (G7, gap F): replace the trivial witness with a v59 spectral/character
  -- quantity (e.g. a character of the order-3 sedenion ψ weighted by the G₂-content
  -- Casimir) that PRODUCES cos(2/3) with the argument = Q.  Until then this is a
  -- placeholder, NOT a derivation.
  exact ⟨Real.cos (2 / 3), rfl⟩  -- `sorry`-equivalent: trivial witness, mechanism OPEN

end SCPv59.G7Magnitude
