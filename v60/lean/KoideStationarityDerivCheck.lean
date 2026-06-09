/-
Copyright (c) 2026  SCP Project
Released under Apache 2.0 license.

v60/lean/KoideStationarityDerivCheck.lean

Companion to KoideStationarity.lean.  A SELF-CONTAINED, GENUINE Mathlib `deriv`
verification that the closed-form gradient `gradV1` used there really IS the
analytic partial derivative `∂V/∂x1` of the GEN3 matter potential `V`.  We build
`HasDerivAt` from primitive Mathlib lemmas (`hasDerivAt_id`, `.pow`, `.const_mul`,
`.add/.sub`) and read off `deriv` via `HasDerivAt.deriv` — i.e. Lean's calculus,
not a hand-asserted formula.  This pins the correct coefficient `2 lam`
(chain rule: d(lam K^2) = 2 lam K K').

Run from v59/furey_construction/lean:
  lake env lean ../../../v60/lean/KoideStationarityDerivCheck.lean
-/
import Mathlib

noncomputable def Vf (lam mu c x1 x2 x3 : ℝ) : ℝ :=
  lam * ((x1 + x2 + x3)^2 - 6*(x1*x2 + x1*x3 + x2*x3))^2 + mu * ((x1 + x2 + x3) - c)^2

def gradV1f (lam mu c x1 x2 x3 : ℝ) : ℝ :=
  2 * lam * ((x1 + x2 + x3)^2 - 6*(x1*x2 + x1*x3 + x2*x3)) * (6*x1 - 4*(x1 + x2 + x3))
    + 2*mu*((x1 + x2 + x3) - c)

example (lam mu c x2 x3 x1 : ℝ) :
    deriv (fun t => Vf lam mu c t x2 x3) x1 = gradV1f lam mu c x1 x2 x3 := by
  unfold gradV1f
  have e : HasDerivAt (fun t : ℝ => t + x2 + x3) 1 x1 := by
    simpa using ((hasDerivAt_id x1).add_const x2).add_const x3
  have k : HasDerivAt
      (fun t : ℝ => (t + x2 + x3)^2 - 6*(t*x2 + t*x3 + x2*x3))
      (6*x1 - 4*(x1 + x2 + x3)) x1 := by
    have p2 : HasDerivAt (fun t : ℝ => (t + x2 + x3)^2) (2 * (x1 + x2 + x3) * 1) x1 := by
      simpa using e.pow 2
    have lin : HasDerivAt (fun t : ℝ => 6*(t*x2 + t*x3 + x2*x3)) (6*(x2 + x3)) x1 := by
      have a1 : HasDerivAt (fun t : ℝ => t*x2) x2 x1 := by simpa using (hasDerivAt_id x1).mul_const x2
      have a2 : HasDerivAt (fun t : ℝ => t*x3) x3 x1 := by simpa using (hasDerivAt_id x1).mul_const x3
      have s : HasDerivAt (fun t : ℝ => t*x2 + t*x3 + x2*x3) (x2 + x3) x1 := by
        simpa using (a1.add a2).add_const (x2*x3)
      simpa using s.const_mul 6
    have := p2.sub lin; convert this using 1; ring
  have v1 : HasDerivAt (fun t : ℝ => lam * ((t + x2 + x3)^2 - 6*(t*x2 + t*x3 + x2*x3))^2)
      (lam * (2 * ((x1 + x2 + x3)^2 - 6*(x1*x2 + x1*x3 + x2*x3)) * (6*x1 - 4*(x1 + x2 + x3)))) x1 := by
    have h := (k.pow 2).const_mul lam
    convert h using 1
    simp
  have v2 : HasDerivAt (fun t : ℝ => mu * ((t + x2 + x3) - c)^2)
      (mu * (2 * ((x1 + x2 + x3) - c) * 1)) x1 := by
    have h := ((e.sub_const c).pow 2).const_mul mu
    convert h using 1
    simp
  have hV : HasDerivAt (fun t : ℝ => Vf lam mu c t x2 x3)
      (lam * (2 * ((x1 + x2 + x3)^2 - 6*(x1*x2 + x1*x3 + x2*x3)) * (6*x1 - 4*(x1 + x2 + x3)))
        + mu * (2 * ((x1 + x2 + x3) - c) * 1)) x1 := by
    have := v1.add v2; unfold Vf; convert this using 1
  rw [hV.deriv]; ring
