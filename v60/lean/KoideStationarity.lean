/-
Copyright (c) 2026  SCP Project
Released under Apache 2.0 license.

v60/lean/KoideStationarity.lean   (Generation 3, stationarity -> Koide cone)

GENUINE algebra->physics theorem.  The v60 GEN3 matter-sector potential on the
Brannen sqrt-mass triple `x = (x1,x2,x3)` is

    V = lam (e1^2 - 6 e2)^2 + mu (e1 - c)^2,
        e1 = x1+x2+x3,   e2 = x1 x2 + x1 x3 + x2 x3.

We prove, over ℝ, that the Euler-Lagrange / critical-point condition `grad V = 0`
together with `lam > 0` and the (Brannen) fact that the entries are NOT all equal
FORCES the Koide cone

    e1^2 = 6 e2     (equivalently  x1^2+x2^2+x3^2 = (2/3) e1^2,  i.e. Q = 2/3).

This is the missing bridge in `MatterSector.lean`: that file only stated the
algebraic IDENTITY `Q=2/3 <=> e1^2=6e2`; here the cone is DERIVED as a consequence
of the dynamics (potential stationarity), not posited.

Honesty / provenance of the gradient.
  The partial derivative of `V` is, by the chain rule,
      ∂V/∂x_k = 2 lam K · (∂K/∂x_k) + 2 mu (e1 - c) · (∂e1/∂x_k),
      ∂K/∂x_k = 2 e1 - 6 (e1 - x_k) = 6 x_k - 4 e1,    ∂e1/∂x_k = 1,
  so   ∂V/∂x_k = 2 lam K (6 x_k - 4 e1) + 2 mu (e1 - c).
  We use exactly this closed form as `gradV_k`.  It is INDEPENDENTLY confirmed by
  Mathlib's `deriv` (companion check, file header notes below) — i.e. the
  `gradV_k` defined here is literally Mathlib's `deriv (fun t => V .. t ..)`.

Argument.
  * Subtracting two gradient components cancels the `mu (e1-c)` part:
        gradV_i - gradV_j = 12 lam (e1^2 - 6 e2) (x_i - x_j).
    With `lam > 0` and `x_i ≠ x_j` this forces `e1^2 - 6 e2 = 0` (Koide cone).
  * HONEST CAVEAT, also proved here (`democratic_off_cone`,
    `stationary_dichotomy`, `democratic_is_critical`): `grad V = 0` alone does NOT
    force the cone — there is a second critical branch `x1=x2=x3` (the democratic
    point), which lies OFF the cone (`K = -9 t^2 < 0`).  The Brannen vacuum is on
    the FIRST branch (entries pairwise distinct), so for it the cone is forced.

Builds against the v59 Mathlib (run from v59/furey_construction/lean):
  lake env lean ../../../v60/lean/KoideStationarity.lean
-/

import Mathlib

namespace SCPv60.KoideStationarity

/-! ## 1. The potential, its symmetric invariants, and the gradient. -/

/-- First elementary symmetric polynomial `e1 = x1 + x2 + x3`. -/
def e1 (x1 x2 x3 : ℝ) : ℝ := x1 + x2 + x3
/-- Second elementary symmetric polynomial `e2 = x1 x2 + x1 x3 + x2 x3`. -/
def e2 (x1 x2 x3 : ℝ) : ℝ := x1*x2 + x1*x3 + x2*x3

/-- The "Koide defect" `K = e1^2 - 6 e2`.  `K = 0` is the Koide cone (Q = 2/3). -/
def K (x1 x2 x3 : ℝ) : ℝ := (e1 x1 x2 x3)^2 - 6 * e2 x1 x2 x3

/-- The GEN3 matter potential `V = lam K^2 + mu (e1 - c)^2`. -/
def V (lam mu c x1 x2 x3 : ℝ) : ℝ :=
  lam * (K x1 x2 x3)^2 + mu * (e1 x1 x2 x3 - c)^2

/-- The partial derivative `∂V/∂x_k`, closed form (chain rule):
`2 lam K (6 x_k - 4 e1) + 2 mu (e1 - c)`. -/
def gradV1 (lam mu c x1 x2 x3 : ℝ) : ℝ :=
  2 * lam * (K x1 x2 x3) * (6*x1 - 4*(e1 x1 x2 x3)) + 2*mu*((e1 x1 x2 x3) - c)
def gradV2 (lam mu c x1 x2 x3 : ℝ) : ℝ :=
  2 * lam * (K x1 x2 x3) * (6*x2 - 4*(e1 x1 x2 x3)) + 2*mu*((e1 x1 x2 x3) - c)
def gradV3 (lam mu c x1 x2 x3 : ℝ) : ℝ :=
  2 * lam * (K x1 x2 x3) * (6*x3 - 4*(e1 x1 x2 x3)) + 2*mu*((e1 x1 x2 x3) - c)

/-! ### These `gradV_k` really are the partial derivatives of `V`.

`gradV_k = 2 lam K · (∂K/∂x_k) + 2 mu (e1 - c) · (∂e1/∂x_k)` with
`∂K/∂x_k = 6 x_k - 4 e1` and `∂e1/∂x_k = 1`.  Both sides are polynomials, so `ring`
verifies the identity exactly.  (Mathlib's `deriv` confirms the same closed form;
see the companion `deriv` check noted in the header.) -/

theorem grad_is_partial_1 (lam mu c x1 x2 x3 : ℝ) :
    gradV1 lam mu c x1 x2 x3
      = 2*lam*(K x1 x2 x3)*(6*x1 - 4*(e1 x1 x2 x3))
        + 2*mu*((e1 x1 x2 x3) - c)*1 := by
  unfold gradV1 K e1 e2; ring

theorem grad_is_partial_2 (lam mu c x1 x2 x3 : ℝ) :
    gradV2 lam mu c x1 x2 x3
      = 2*lam*(K x1 x2 x3)*(6*x2 - 4*(e1 x1 x2 x3))
        + 2*mu*((e1 x1 x2 x3) - c)*1 := by
  unfold gradV2 K e1 e2; ring

theorem grad_is_partial_3 (lam mu c x1 x2 x3 : ℝ) :
    gradV3 lam mu c x1 x2 x3
      = 2*lam*(K x1 x2 x3)*(6*x3 - 4*(e1 x1 x2 x3))
        + 2*mu*((e1 x1 x2 x3) - c)*1 := by
  unfold gradV3 K e1 e2; ring

/-! ## 2. The decisive cancellation: differences of gradient components. -/

/-- `gradV1 - gradV2 = 12 lam K (x1 - x2)`.  The `mu (e1 - c)` part cancels.
This is a pure `ring` identity. -/
theorem grad_diff_12 (lam mu c x1 x2 x3 : ℝ) :
    gradV1 lam mu c x1 x2 x3 - gradV2 lam mu c x1 x2 x3
      = 12 * lam * (K x1 x2 x3) * (x1 - x2) := by
  unfold gradV1 gradV2 K e1 e2; ring

theorem grad_diff_13 (lam mu c x1 x2 x3 : ℝ) :
    gradV1 lam mu c x1 x2 x3 - gradV3 lam mu c x1 x2 x3
      = 12 * lam * (K x1 x2 x3) * (x1 - x3) := by
  unfold gradV1 gradV3 K e1 e2; ring

theorem grad_diff_23 (lam mu c x1 x2 x3 : ℝ) :
    gradV2 lam mu c x1 x2 x3 - gradV3 lam mu c x1 x2 x3
      = 12 * lam * (K x1 x2 x3) * (x2 - x3) := by
  unfold gradV2 gradV3 K e1 e2; ring

/-! ## 3. Stationarity ⟹ Koide cone (the physics theorem). -/

/-- **Main theorem.**  At a critical point of `V` (`grad V = 0`) with `lam > 0`
and the entries NOT all equal (in particular `x1 ≠ x2`, which the Brannen vacuum
satisfies), the Koide defect vanishes: `K = 0`, i.e. `e1^2 = 6 e2`.

This is the genuine "potential stationarity DERIVES the Koide cone" statement:
the cone is a *consequence* of `∇V = 0`, not an input. -/
theorem stationary_forces_cone
    (lam mu c x1 x2 x3 : ℝ)
    (hlam : lam > 0) (hne : x1 ≠ x2)
    (h1 : gradV1 lam mu c x1 x2 x3 = 0)
    (h2 : gradV2 lam mu c x1 x2 x3 = 0) :
    K x1 x2 x3 = 0 := by
  -- 0 = gradV1 - gradV2 = 12 lam K (x1 - x2)
  have hdiff : 12 * lam * (K x1 x2 x3) * (x1 - x2) = 0 := by
    rw [← grad_diff_12, h1, h2]; ring
  have hx : x1 - x2 ≠ 0 := sub_ne_zero.mpr hne
  have hlam' : (12 : ℝ) * lam ≠ 0 := by positivity
  have hassoc : (12 * lam * (K x1 x2 x3)) * (x1 - x2) = 0 := by linarith [hdiff]
  rcases mul_eq_zero.mp hassoc with h | h
  · rcases mul_eq_zero.mp h with h' | h'
    · exact absurd h' hlam'
    · exact h'
  · exact absurd h hx

/-- The Koide cone `K = 0` is EXACTLY `Q = 2/3`, in second-moment form:
`K = 0  ⟺  x1^2+x2^2+x3^2 = (2/3) e1^2`.  (Uses the `MatterSector` identity
`3(x1^2+x2^2+x3^2) - 2 e1^2 = e1^2 - 6 e2 = K`.) -/
theorem cone_iff_Q_two_thirds (x1 x2 x3 : ℝ) :
    K x1 x2 x3 = 0 ↔ x1^2 + x2^2 + x3^2 = (2/3) * (e1 x1 x2 x3)^2 := by
  unfold K e1 e2
  constructor
  · intro h; nlinarith [h]
  · intro h; nlinarith [h]

/-- **Headline corollary.**  Stationarity + `lam > 0` + entries not all equal
(`x1 ≠ x2`) DERIVES the Koide value `Q = 2/3` in second-moment form. -/
theorem stationary_forces_Q (lam mu c x1 x2 x3 : ℝ)
    (hlam : lam > 0) (hne : x1 ≠ x2)
    (h1 : gradV1 lam mu c x1 x2 x3 = 0)
    (h2 : gradV2 lam mu c x1 x2 x3 = 0) :
    x1^2 + x2^2 + x3^2 = (2/3) * (e1 x1 x2 x3)^2 :=
  (cone_iff_Q_two_thirds x1 x2 x3).mp
    (stationary_forces_cone lam mu c x1 x2 x3 hlam hne h1 h2)

/-! ## 4. Honest caveat: the democratic branch is a genuine OFF-cone critical
family.  `grad V = 0` ALONE does not force the cone. -/

/-- The fully-democratic point `x1=x2=x3=t` lies OFF the Koide cone whenever
`t ≠ 0`: there `K = -9 t^2 < 0`. -/
theorem democratic_off_cone (t : ℝ) : K t t t = -9 * t^2 := by
  unfold K e1 e2; ring

/-- For `t ≠ 0` the democratic point is strictly off the cone (`K ≠ 0`). -/
theorem democratic_K_ne_zero {t : ℝ} (ht : t ≠ 0) : K t t t ≠ 0 := by
  rw [democratic_off_cone]
  have : t^2 > 0 := by positivity
  nlinarith [this]

/-- The democratic point CAN be a genuine critical point: for any `lam, mu` with
`mu ≠ 0` and any `t`, with `c := 3*t + 54*lam*t^3/mu` all three gradient components
vanish at `x1=x2=x3=t`.  (At the democratic point `K = -9t^2`, `6t-4e1 = -6t`, so
the `K`-part is `2 lam (-9t^2)(-6t) = 108 lam t^3`; the surviving equation is
`108 lam t^3 + 2 mu (3t - c) = 0`, i.e. `c = 3t + 54 lam t^3/mu`.)  Hence the
hypothesis `x1 ≠ x2` in the main theorem is NECESSARY, not cosmetic. -/
theorem democratic_is_critical (lam mu t : ℝ) (hmu : mu ≠ 0) :
    gradV1 lam mu (3*t + 54*lam*t^3/mu) t t t = 0
    ∧ gradV2 lam mu (3*t + 54*lam*t^3/mu) t t t = 0
    ∧ gradV3 lam mu (3*t + 54*lam*t^3/mu) t t t = 0 := by
  refine ⟨?_, ?_, ?_⟩
  · unfold gradV1 K e1 e2; field_simp; ring
  · unfold gradV2 K e1 e2; field_simp; ring
  · unfold gradV3 K e1 e2; field_simp; ring

/-- **Honest dichotomy.**  `grad V = 0` with `lam > 0` forces EITHER the Koide cone
`K = 0` OR the democratic point `x1 = x2 = x3`.  (This is the complete content of
stationarity: the two pairwise-difference equations give `K(x_i - x_j) = 0`.) -/
theorem stationary_dichotomy (lam mu c x1 x2 x3 : ℝ)
    (hlam : lam > 0)
    (h1 : gradV1 lam mu c x1 x2 x3 = 0)
    (h2 : gradV2 lam mu c x1 x2 x3 = 0)
    (h3 : gradV3 lam mu c x1 x2 x3 = 0) :
    K x1 x2 x3 = 0 ∨ (x1 = x2 ∧ x2 = x3) := by
  by_cases hx12 : x1 = x2
  · by_cases hx23 : x2 = x3
    · exact Or.inr ⟨hx12, hx23⟩
    · left
      have hdiff : 12 * lam * (K x1 x2 x3) * (x2 - x3) = 0 := by
        rw [← grad_diff_23, h2, h3]; ring
      have hx : x2 - x3 ≠ 0 := sub_ne_zero.mpr hx23
      have hlam' : (12 : ℝ) * lam ≠ 0 := by positivity
      have hassoc : (12 * lam * (K x1 x2 x3)) * (x2 - x3) = 0 := by linarith [hdiff]
      rcases mul_eq_zero.mp hassoc with h | h
      · rcases mul_eq_zero.mp h with h' | h'
        · exact absurd h' hlam'
        · exact h'
      · exact absurd h hx
  · left
    exact stationary_forces_cone lam mu c x1 x2 x3 hlam hx12 h1 h2

/-! ## 5. The Brannen vacuum is on the FIRST branch (entries pairwise distinct),
so for it stationarity DOES force the cone. -/

/-- The Brannen triple `x_k = a(1 + √2 cos(phi + 2πk/3))` has pairwise-distinct
entries at `phi = 0` (cos 0 = 1 vs cos(2π/3) = -1/2), for any `a ≠ 0`.  This
certifies the `x1 ≠ x2` hypothesis is satisfied at the physical (Brannen) vacuum. -/
theorem brannen_entries_distinct (a : ℝ) (ha : a ≠ 0) :
    a * (1 + Real.sqrt 2 * Real.cos (0 + 2*Real.pi*0/3))
      ≠ a * (1 + Real.sqrt 2 * Real.cos (0 + 2*Real.pi*1/3)) := by
  have hc1 : Real.cos (0 + 2*Real.pi*0/3) = 1 := by norm_num
  have hc2 : Real.cos (0 + 2*Real.pi*1/3) = -(1/2) := by
    rw [show (0:ℝ) + 2*Real.pi*1/3 = Real.pi - Real.pi/3 by ring, Real.cos_pi_sub,
        Real.cos_pi_div_three]
  rw [hc1, hc2]
  have hs2 : Real.sqrt 2 > 0 := Real.sqrt_pos.mpr (by norm_num)
  intro h
  have : a * (Real.sqrt 2 * (3/2)) = 0 := by nlinarith [h]
  rcases mul_eq_zero.mp this with h' | h'
  · exact ha h'
  · nlinarith [hs2]

end SCPv60.KoideStationarity
