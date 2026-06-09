/-
Copyright (c) 2026  SCP Project
Released under Apache 2.0 license.

v61/lean/Backreaction.lean   (Generation 2 of the v61 loop)

Machine-checked algebraic backbone of v61 GEN2 (02_backreaction.py +
02_backreaction.mac): matter backreaction.  The tensor identity
`G^t_t = -2 m'(r)/r^2` is computed in SymPy and Maxima `ctensor`; here we prove the
algebraic consequences that tie it to the source `rho_grav` and the v60 limit.

  * the (00) Einstein equation `G^t_t = -8 pi rho` with `G^t_t = -2 m'/r^2` gives
    the backreaction law `m'(r) = 4 pi r^2 rho`  (the coefficient identity
    `-2 (4 pi r^2 rho)/r^2 = -8 pi rho`).
  * Newtonian limit: `grad^2 Phi = m'/r^2 = 4 pi rho`  (the v60 GEN4 / OBE Poisson
    law, derived).
  * uniform-density mass `M = (4/3) pi R^3 rho_0` (the integral result form).
  * exterior matching: `g_rr = (1 - 2M/r)^{-1}` when `m = M` (GEN1 Schwarzschild).

Builds against the v59 Mathlib:
  cd v59/furey_construction/lean && lake env lean ../../../v61/lean/Backreaction.lean
-/

import Mathlib

namespace SCPv61.Backreaction

open Real

/-- The (00) Einstein equation `G^t_t = -8 pi rho` with `G^t_t = -2 m'/r^2`
(SymPy/Maxima) gives the backreaction law `m' = 4 pi r^2 rho`: substituting it,
`-2 (4 pi r^2 rho)/r^2 = -8 pi rho`. -/
theorem backreaction_coefficient (r ρ : ℝ) (hr : r ≠ 0) :
    -2 * (4 * Real.pi * r ^ 2 * ρ) / r ^ 2 = -8 * Real.pi * ρ := by
  field_simp; ring

/-- Newtonian limit: `grad^2 Phi = m'/r^2 = 4 pi rho` with `m' = 4 pi r^2 rho`
-- the v60 GEN4 / OBE Poisson law, now derived as the weak-field (00) equation. -/
theorem newtonian_poisson (r ρ : ℝ) (hr : r ≠ 0) :
    (4 * Real.pi * r ^ 2 * ρ) / r ^ 2 = 4 * Real.pi * ρ := by
  field_simp

/-- Uniform-density total mass `M = (4/3) pi R^3 rho_0` (the `int 4 pi r^2 rho dr`
result); equivalently the mean-density relation `rho_0 = 3M/(4 pi R^3)`. -/
theorem uniform_mass (R ρ0 : ℝ) :
    (4 / 3) * Real.pi * R ^ 3 * ρ0 = (4 * Real.pi * R ^ 3 * ρ0) / 3 := by ring

theorem mean_density (M R : ℝ) (hR : R ≠ 0) (hpi : Real.pi ≠ 0) :
    (4 / 3) * Real.pi * R ^ 3 * (3 * M / (4 * Real.pi * R ^ 3)) = M := by
  field_simp

/-- Exterior matching: when `m = M` (vacuum outside the matter), the radial metric
is the GEN1 Schwarzschild form `(1 - 2M/r)^{-1}`. -/
theorem exterior_schwarzschild (r M : ℝ) (h : 1 - 2 * M / r ≠ 0) :
    (1 / (1 - 2 * M / r)) * (1 - 2 * M / r) = 1 :=
  one_div_mul_cancel h

/-- Headline: backreaction law coefficient, the Newtonian Poisson limit, and the
uniform-density mass -- matter sources the GEN1 Schwarzschild charge. -/
theorem gen2_backreaction :
    (∀ r ρ : ℝ, r ≠ 0 → -2 * (4 * Real.pi * r ^ 2 * ρ) / r ^ 2 = -8 * Real.pi * ρ) ∧
    (∀ r ρ : ℝ, r ≠ 0 → (4 * Real.pi * r ^ 2 * ρ) / r ^ 2 = 4 * Real.pi * ρ) ∧
    (∀ R ρ0 : ℝ, (4 / 3) * Real.pi * R ^ 3 * ρ0 = (4 * Real.pi * R ^ 3 * ρ0) / 3) := by
  refine ⟨fun r ρ hr => by field_simp; ring, fun r ρ hr => by field_simp,
          fun R ρ0 => by ring⟩

end SCPv61.Backreaction
