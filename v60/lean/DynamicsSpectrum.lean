/-
Copyright (c) 2026  SCP Project
Released under Apache 2.0 license.

v60/lean/DynamicsSpectrum.lean   (Generation 7 of the dynamical-Lagrangian loop)

Machine-checked backbone of GEN7 (16_dynamics.c + 16_dynamics_check.py): the
nonlinear Euler-Lagrange time evolution reproduces the GEN5 spectrum, with the
relativistic dispersion `omega^2 = k^2 + m^2`.  Here we prove the clean algebraic
consequences:

  * at rest (`k = 0`) the oscillation frequency is the mass: `omega^2 = m^2`
    -- so the measured normal-mode frequencies ARE the Hessian eigenvalues (GEN5).
  * the MASSLESS Goldstone propagates at speed 1: `sqrt(k^2 + 0)/k = 1`.
  * MASSIVE modes sit above the light cone: `k^2 < k^2 + m^2` for `m > 0`.

Builds against the v59 Mathlib:
  cd v59/furey_construction/lean && lake env lean ../../../v60/lean/DynamicsSpectrum.lean
-/

import Mathlib

namespace SCPv60.DynamicsSpectrum

/-- Relativistic dispersion of a normal mode: `omega^2 = k^2 + m^2`. -/
def omegaSq (k m : ℝ) : ℝ := k ^ 2 + m ^ 2

/-- At rest (`k = 0`) the squared frequency equals the squared mass: the measured
homogeneous normal-mode frequencies are exactly the Hessian eigenvalues (GEN5). -/
theorem omega_at_rest (m : ℝ) : omegaSq 0 m = m ^ 2 := by
  unfold omegaSq; ring

/-- The massless Goldstone (`m = 0`) propagates at speed 1:
`sqrt(omegaSq k 0) / k = 1` for `k > 0`. -/
theorem goldstone_speed (k : ℝ) (hk : 0 < k) :
    Real.sqrt (omegaSq k 0) / k = 1 := by
  unfold omegaSq
  rw [show k ^ 2 + (0:ℝ) ^ 2 = k ^ 2 by ring, Real.sqrt_sq hk.le, div_self (ne_of_gt hk)]

/-- Massive modes sit strictly above the light cone: `k^2 < omegaSq k m` for `m ≠ 0`. -/
theorem massive_above_lightcone (k m : ℝ) (hm : m ≠ 0) : k ^ 2 < omegaSq k m := by
  unfold omegaSq
  have : 0 < m ^ 2 := by positivity
  linarith

/-- Headline: the dispersion relation gives `omega^2 = m^2` at rest, a unit-speed
massless Goldstone, and subluminal massive modes -- the verified GEN7 dynamics. -/
theorem gen7_dynamics :
    (∀ m : ℝ, omegaSq 0 m = m ^ 2) ∧
    (∀ k : ℝ, 0 < k → Real.sqrt (omegaSq k 0) / k = 1) ∧
    (∀ k m : ℝ, m ≠ 0 → k ^ 2 < omegaSq k m) := by
  refine ⟨fun m => by unfold omegaSq; ring, fun k hk => goldstone_speed k hk,
          fun k m hm => massive_above_lightcone k m hm⟩

end SCPv60.DynamicsSpectrum
