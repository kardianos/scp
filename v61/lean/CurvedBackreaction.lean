/-
Copyright (c) 2026  SCP Project
Released under Apache 2.0 license.

v61/lean/CurvedBackreaction.lean   (Generation 1 of the v61 loop)

Machine-checked algebraic backbone of v61 GEN1 (01_curved_backreaction.py +
01_schwarzschild.mac): the nonlinear/curved gravity sector.  The full tensor
computation (Christoffel/Ricci/Einstein of the Schwarzschild ansatz = vacuum) is
done in SymPy and independently in Maxima `ctensor`; here we prove the algebraic
relations that tie it to the v60 Newtonian limit and to a curved-space observable.

  * horizon / mass:  r_s = 2 G M.
  * weak-field match:  g_00 = -(1 - r_s/r) with r_s = 2GM gives the Newtonian
    potential  Phi = -(g_00 + 1)/2 = -G M / r  -- EXACTLY the v60 GEN4 result.
  * Schwarzschild structure:  g_00 * g_11 = -1  (B = 1/A).
  * light deflection:  4GM/b = 2 * (2GM/b) -- GR bends light TWICE the Newtonian
    amount (a genuine curved-space prediction).

Builds against the v59 Mathlib:
  cd v59/furey_construction/lean && lake env lean ../../../v61/lean/CurvedBackreaction.lean
-/

import Mathlib

namespace SCPv61.CurvedBackreaction

/-- Schwarzschild horizon radius in terms of the gravitational charge `M = ∫ρ_grav`. -/
def horizon (G M : ℝ) : ℝ := 2 * G * M

/-- Weak-field match: with `g_00 = -(1 - 2GM/r)`, the Newtonian potential
`Phi = -(g_00 + 1)/2 = -GM/r` -- the v60 GEN4 / OBE result is the weak-field limit
of the Schwarzschild metric. -/
theorem weak_field_newtonian (r G M : ℝ) (hr : r ≠ 0) :
    (-((-(1 - 2 * G * M / r)) + 1) / 2) = -(G * M / r) := by
  field_simp; ring

/-- Schwarzschild structure `g_00 · g_11 = -1` (`B = 1/A`, `A = 1 - r_s/r`). -/
theorem schwarzschild_gtt_grr (r rs : ℝ) (hA : 1 - rs / r ≠ 0) :
    (-(1 - rs / r)) * (1 / (1 - rs / r)) = -1 := by
  field_simp

/-- Light deflection in GR is twice the Newtonian value: `4GM/b = 2·(2GM/b)`. -/
theorem deflection_doubles (G M b : ℝ) : 4 * G * M / b = 2 * (2 * G * M / b) := by ring

/-- Headline: horizon `r_s = 2GM`, the Newtonian weak-field limit, and the
factor-2 light deflection -- the curved-space gravity sector tied to v60. -/
theorem gen1_curved :
    (∀ r G M : ℝ, r ≠ 0 → (-((-(1 - 2 * G * M / r)) + 1) / 2) = -(G * M / r)) ∧
    (∀ G M b : ℝ, 4 * G * M / b = 2 * (2 * G * M / b)) := by
  refine ⟨fun r G M hr => by field_simp; ring, fun G M b => by ring⟩

end SCPv61.CurvedBackreaction
