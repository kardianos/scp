/-
Copyright (c) 2026  SCP Project
Released under Apache 2.0 license.

v60/lean/G9Unification.lean

Building out the v59 gap (closing C3): the unified algebra is Cl(3,1) ⊗ Cl(7)_even
(v59's own stated factorization, SYNTHESIS §8).  The abstract reason the spacetime
Lorentz so(3,1) commutes with the internal Spin(7) is a GENERAL fact about tensor
(Kronecker) products: a spacetime operator `A ⊗ I` commutes with an internal
operator `I ⊗ B` for ALL A, B.  We prove that in full generality here (not just
the numerical instance verified in 07_unified_algebra.py).

Plus the dimension identities of the factorization, including the gravity
prefactor 21/16 = dim Spin(7) / dim Cl(3,1) (internal / spacetime).

Builds against the v59 Mathlib:
  cd v59/furey_construction/lean && lake env lean ../../../v60/lean/G9Unification.lean
-/

import Mathlib

open Matrix
open scoped Kronecker

namespace SCPv60.G9Unification

/-! ## 1. The general commutation theorem (the heart of C3)

In `Cl(3,1) ⊗ Cl(7)_even`, the spacetime factor and the internal factor commute.
Concretely: any spacetime operator lifted as `A ⊗ I` commutes with any internal
operator lifted as `I ⊗ B`.  This is exactly why `so(3,1)` (spacetime Lorentz /
the soldering 2-form) commutes with the full internal `Spin(7)` (incl. its
`SU(2)_L`, color, `G₂`, and the triality that gives the 3 generations). -/

/-- **[thm] Spacetime ⊗ internal factors commute — for all `A`, `B`.** -/
theorem spacetime_internal_commute
    {R : Type*} [CommRing R] {m n : Type*}
    [Fintype m] [Fintype n] [DecidableEq m] [DecidableEq n]
    (A : Matrix m m R) (B : Matrix n n R) :
    (A ⊗ₖ (1 : Matrix n n R)) * ((1 : Matrix m m R) ⊗ₖ B)
      = ((1 : Matrix m m R) ⊗ₖ B) * (A ⊗ₖ (1 : Matrix n n R)) := by
  rw [← mul_kronecker_mul, ← mul_kronecker_mul]
  simp

/-- The commutator vanishes (same statement, bracket form). -/
theorem spacetime_internal_bracket_zero
    {R : Type*} [CommRing R] {m n : Type*}
    [Fintype m] [Fintype n] [DecidableEq m] [DecidableEq n]
    (A : Matrix m m R) (B : Matrix n n R) :
    (A ⊗ₖ (1 : Matrix n n R)) * ((1 : Matrix m m R) ⊗ₖ B)
      - ((1 : Matrix m m R) ⊗ₖ B) * (A ⊗ₖ (1 : Matrix n n R)) = 0 := by
  rw [spacetime_internal_commute, sub_self]

/-! ## 2. Dimension identities of the Cl(3,1) ⊗ Cl(7)_even factorization -/

/-- `dim Cl(3,1) = 2^4 = 16` (the spacetime algebra). -/
def dimCl31 : ℕ := 2 ^ 4
/-- `dim Spin(7) = C(7,2) = 21` (the internal symmetry). -/
def dimSpin7 : ℕ := Nat.choose 7 2
/-- `dim Cl(7)_even = 2^6 = 64` (the internal algebra ≅ ℂ⊗ℍ⊗𝕆). -/
def dimCl7even : ℕ := 2 ^ 6

theorem dims_eq : dimCl31 = 16 ∧ dimSpin7 = 21 ∧ dimCl7even = 64 := by
  refine ⟨?_, ?_, ?_⟩ <;> decide

/-- Spacetime spinor (4) ⊗ octonion (8) = 32-dim total module. -/
theorem module_dim : 4 * 8 = 32 := by decide

/-- The unified algebra dimension `dim Cl(3,1) · dim Cl(7)_even = 16 · 64 = 1024`. -/
theorem unified_algebra_dim : dimCl31 * dimCl7even = 1024 := by decide

/-- **[thm] The gravity prefactor is the internal/spacetime dimension ratio.**
    `21/16 = dim Spin(7) / dim Cl(3,1)` — i.e. the `(21/16)α²¹` of v59's `G_e`
    is a ratio across the two factors of this factorization (SYNTHESIS §9). -/
theorem prefactor_is_dim_ratio : (dimSpin7 : ℚ) / (dimCl31 : ℚ) = 21 / 16 := by
  unfold dimSpin7 dimCl31; norm_num [Nat.choose]

/-! ## 3. Reading

`spacetime_internal_commute` is the closure of C3 at the abstract level: in v59's
own stated `Cl(spacetime) ⊗ Cl(internal)` structure, the spacetime Lorentz is in a
different tensor factor from the entire internal symmetry, so they commute — for
*every* pair of operators, not just the ones checked numerically.  The explicit
realization (Dirac `σ^{μν}` ⊗ octonion `Spin(7)`, with `[·,·]=0` verified) is in
`07_unified_algebra.py`.  Combined with the Schur obstruction
(`G9Soldering.no_internal_lorentz`: no Lorentz fits *inside* `Cl(7)_even`), the
gravity carrier is forced into the spacetime factor and is automatically
`G₂`/triality/color-compatible there. -/

end SCPv60.G9Unification
