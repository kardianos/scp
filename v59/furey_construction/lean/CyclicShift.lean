/-
  v59/furey_construction/lean/CyclicShift.lean

  The Z_3 algebraic skeleton of the lepton kernel.

  The lepton mass operator of v59 is
      M(a, ξ) = a · (I_3 + ξ · S + ξ̄ · S^T)
  on ℂ^3, with S the 3-cycle (cyclic shift) matrix.  Its eigenvectors are the
  discrete-Fourier basis  v_k = (1, ω^k, ω^{2k})  and its eigenvalues are
      λ_k = a · (1 + ξ · ω^k + ξ̄ · ω^{-k}).
  Writing ξ = t · exp(i φ) and applying Euler's formula on the spectral pair,
      ξ ω^k + ξ̄ ω^{-k} = 2 t cos(2π k / 3 + φ),
  so λ_k coincides with the Brannen amplitude `s a t φ k` of `BrannenKernel.lean`.
  Combined with `BrannenKernel.Q_value`, this gives Koide  Q = (1 + 2 t²)/3,
  and `BrannenKernel.koide_iff_constraint` then gives  Q = 2/3 ⟺ t² = 1/2.

  This file establishes the Z_3 ring-theoretic facts the eigenvalue derivation
  rests on:
    (1) `ω` defined as `exp(2π i / 3)` is a primitive 3rd root of unity.
    (2) `ω^3 = 1`.
    (3) `ω ≠ 1`.
    (4) `1 + ω + ω^2 = 0`.
  The remaining matrix-diagonalisation and Euler-formula steps are standard
  linear algebra over ℂ and are not re-formalised here.
-/

import Mathlib.Analysis.SpecialFunctions.Complex.Analytic
import Mathlib.RingTheory.RootsOfUnity.Complex
import Mathlib.Tactic

namespace SCPv59.CyclicShift

open Complex Real

/-- The primitive 3rd root of unity `ω = exp(2π i / 3)`. -/
noncomputable def ω : ℂ := Complex.exp (2 * Real.pi * Complex.I / 3)

/-- `ω` is a primitive 3rd root of unity (Mathlib `Complex.isPrimitiveRoot_exp`). -/
lemma ω_isPrimitiveRoot : IsPrimitiveRoot ω 3 :=
  Complex.isPrimitiveRoot_exp 3 (by norm_num)

/-- `ω^3 = 1`. -/
theorem ω_pow_three : ω^3 = 1 := ω_isPrimitiveRoot.pow_eq_one

/-- `ω ≠ 1` (it has order 3, not 1). -/
theorem ω_ne_one : ω ≠ 1 := by
  intro h
  have hord : orderOf ω = 3 := ω_isPrimitiveRoot.eq_orderOf.symm
  have : orderOf ω = 1 := by rw [h]; exact orderOf_one
  omega

/-- **The Z_3 cyclic-shift identity.** `1 + ω + ω² = 0`.

    This is the algebraic root of v59's structural derivation: it is the reason
    the three Brannen amplitudes sum to `3 a` (the ω^k contributions cancel),
    and ultimately the reason the Koide ratio admits the closed form (1+2t²)/3. -/
theorem sum_one_omega_omega_sq : 1 + ω + ω^2 = 0 := by
  have h3 : ω^3 = 1 := ω_pow_three
  have hne : ω ≠ 1 := ω_ne_one
  have hsub : ω - 1 ≠ 0 := sub_ne_zero.mpr hne
  have hfact : (ω - 1) * (1 + ω + ω^2) = ω^3 - 1 := by ring
  rw [h3, sub_self] at hfact
  rcases mul_eq_zero.mp hfact with h | h
  · exact absurd h hsub
  · exact h

/-- `ω` is a cube root: `ω · ω² = 1` (= `ω^3`). -/
theorem ω_mul_ω_sq : ω * ω^2 = 1 := by
  have h3 : ω^3 = 1 := ω_pow_three
  linear_combination h3

end SCPv59.CyclicShift
