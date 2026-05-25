/-
  v59/furey_construction/lean/KernelEigenvalues.lean

  Formalisation of item (3c) from the v59 Lean review: the **matrix
  diagonalisation** of the Brannen cyclic mass-kernel.

  `CyclicShift.lean` and `BrannenKernel.lean` derived the Koide ratio from the
  *assumed* fact that the cyclic operator

      M(a, ξ) = a · (I + ξ S + ξ̄ Sᵀ)   on ℂ³,   S the 3-cycle shift

  has the discrete-Fourier vectors  v_k = (1, ω^k, ω^{2k})  as eigenvectors with
  eigenvalues  λ_k = a (1 + ξ ω^k + ξ̄ ω^{2k}).  That diagonalisation was called
  "standard linear algebra ... not re-derived formally here" in both files.

  This module supplies the missing proof:

    * `Scyc_mulVec_v`   : `S  *ᵥ v_k = ω^k   • v_k`
    * `ScycT_mulVec_v`  : `Sᵀ *ᵥ v_k = ω^{2k} • v_k`   (ω^{2k} = ω^{-k})
    * `M_mulVec_eigen`  : `M(a,ξ) *ᵥ v_k = λ_k • v_k`   (**the diagonalisation**)
    * `lam_eq_brannen`  : with `ξ = t·e^{iφ}` (t,φ real), `λ_k = ↑(s a t φ k)`,
                          i.e. the complex eigenvalues are exactly the real
                          Brannen amplitudes of `BrannenKernel.s`.

  Composing `lam_eq_brannen` with `BrannenKernel.Q_value` / `koide_iff_constraint`
  makes the lepton-sector chain end-to-end: the *matrix* eigenvalues (not just an
  abstract `s` function) carry the Koide ratio `Q = (1 + 2t²)/3`.
-/
import Mathlib.Data.Matrix.Mul
import Mathlib.LinearAlgebra.Matrix.Notation
import Mathlib.Analysis.SpecialFunctions.Complex.Circle
import Mathlib.Tactic
import CyclicShift
import BrannenKernel

namespace SCPv59.KernelEigenvalues

open SCPv59.CyclicShift
open Matrix
open Complex

set_option maxRecDepth 4000
-- The eigenvector verifications are 3×3 case bashes; the `first | …` closers
-- legitimately leave some alternatives unused per goal.
set_option linter.unreachableTactic false
set_option linter.unusedTactic false

/-! ## The discrete-Fourier eigenvectors and the cyclic shift -/

/-- The k-th discrete-Fourier vector `v_k = (ω^{0·k}, ω^{1·k}, ω^{2·k})`. -/
noncomputable def v (k : Fin 3) : Fin 3 → ℂ := fun i => ω ^ (i.val * k.val)

/-- The 3-cycle (cyclic shift) matrix: `(S x)_i = x_{i+1}`. -/
def Scyc : Matrix (Fin 3) (Fin 3) ℂ := Matrix.of fun i j => if j = i + 1 then 1 else 0

/-! ## ω-power arithmetic mod 3 -/

/-- ω-powers depend only on the exponent mod 3 (since `ω³ = 1`). -/
lemma omega_pow_eq {a b : ℕ} (h : a % 3 = b % 3) : ω ^ a = ω ^ b := by
  conv_lhs => rw [← Nat.div_add_mod a 3]
  conv_rhs => rw [← Nat.div_add_mod b 3]
  rw [pow_add, pow_add, pow_mul, pow_mul, ω_pow_three, one_pow, one_pow, h]

/-- A multiple of 3 in the exponent gives `1`. -/
lemma omega_pow_eq_one {a : ℕ} (h : a % 3 = 0) : ω ^ a = 1 := by
  have h2 : ω ^ a = ω ^ 0 := omega_pow_eq (by simpa using h)
  simpa using h2

/-- Closes residual ω-power equalities `ω^a = ω^b` arising from the case bash,
    including the `ω^a = 1` / bare-`ω` (`= ω^1`) normal forms `ring_nf` produces. -/
macro "ω_close" : tactic =>
  `(tactic| first
     | (apply omega_pow_eq; decide)
     | (symm; apply omega_pow_eq_one; decide)
     | (apply omega_pow_eq_one; decide)
     | ((conv_lhs => rw [← pow_one ω]); apply omega_pow_eq; decide)
     | ((conv_rhs => rw [← pow_one ω]); apply omega_pow_eq; decide))

/-! ## The shift acts diagonally on the Fourier basis -/

/-- `S *ᵥ v_k = ω^k • v_k`. -/
lemma Scyc_mulVec_v (k : Fin 3) : Scyc.mulVec (v k) = ω ^ (k.val) • v k := by
  funext i
  fin_cases i <;> fin_cases k <;>
    (simp (config := {decide := true}) [Scyc, Matrix.mulVec, dotProduct,
        Fin.sum_univ_three, Matrix.of_apply, v, Pi.smul_apply, smul_eq_mul] <;>
     ring_nf <;> ω_close)

/-- `Sᵀ *ᵥ v_k = ω^{2k} • v_k` (and `ω^{2k} = ω^{-k}`). -/
lemma ScycT_mulVec_v (k : Fin 3) : Scycᵀ.mulVec (v k) = ω ^ (2 * k.val) • v k := by
  funext i
  fin_cases i <;> fin_cases k <;>
    (simp (config := {decide := true}) [Scyc, Matrix.transpose, Matrix.mulVec, dotProduct,
        Fin.sum_univ_three, Matrix.of_apply, v, Pi.smul_apply, smul_eq_mul] <;>
     ring_nf <;> ω_close)

/-! ## The Brannen mass-kernel and its eigenvalues -/

/-- The Brannen cyclic mass-kernel `M(a,ξ) = a (I + ξ S + ξ̄ Sᵀ)`. -/
noncomputable def M (a ξ : ℂ) : Matrix (Fin 3) (Fin 3) ℂ :=
  a • ((1 : Matrix (Fin 3) (Fin 3) ℂ) + ξ • Scyc + (starRingEnd ℂ ξ) • Scycᵀ)

/-- The closed-form eigenvalue `λ_k = a (1 + ξ ω^k + ξ̄ ω^{2k})`. -/
noncomputable def lam (a ξ : ℂ) (k : Fin 3) : ℂ :=
  a * (1 + ξ * ω ^ (k.val) + (starRingEnd ℂ ξ) * ω ^ (2 * k.val))

/-- **The diagonalisation.** Each Fourier vector `v_k` is an eigenvector of the
    Brannen kernel `M(a,ξ)` with eigenvalue `lam a ξ k`. -/
theorem M_mulVec_eigen (a ξ : ℂ) (k : Fin 3) :
    (M a ξ).mulVec (v k) = (lam a ξ k) • v k := by
  unfold M
  rw [smul_mulVec, add_mulVec, add_mulVec, smul_mulVec, smul_mulVec,
      one_mulVec, Scyc_mulVec_v, ScycT_mulVec_v]
  funext i
  simp only [lam, Pi.add_apply, Pi.smul_apply, smul_eq_mul]
  ring

/-! ## The eigenvalues are the real Brannen amplitudes

When the complex parameter is written in polar form `ξ = t·e^{iφ}` (t, φ real),
the spectral pair `ξ ω^k + ξ̄ ω^{2k}` collapses to `2 t cos(φ + 2π k/3)` (the
`ω^{2k} = ω^{-k}` identity is the `exp(2π k i) = 1` periodicity), so the complex
eigenvalue `lam` becomes the *real* Brannen amplitude `BrannenKernel.s`. -/

/-- The spectral pair in polar form: `ξ ω^k + ξ̄ ω^{2k} = 2 t cos(φ + 2π k/3)`. -/
lemma xi_sum (t φ : ℝ) (k : Fin 3) :
    (t : ℂ) * Complex.exp ((φ : ℂ) * Complex.I) * ω ^ (k.val)
    + (starRingEnd ℂ) ((t : ℂ) * Complex.exp ((φ : ℂ) * Complex.I)) * ω ^ (2 * k.val)
    = ((2 * t * Real.cos (φ + 2 * Real.pi * (k.val : ℝ) / 3) : ℝ) : ℂ) := by
  have hconj : (starRingEnd ℂ) ((t : ℂ) * Complex.exp ((φ : ℂ) * Complex.I))
      = (t : ℂ) * Complex.exp (-(φ : ℂ) * Complex.I) := by
    rw [map_mul, ← Complex.exp_conj]; simp [Complex.conj_ofReal]
  have homega : ∀ n : ℕ, ω ^ n = Complex.exp ((n : ℂ) * (2 * (Real.pi : ℂ) * Complex.I / 3)) := by
    intro n; rw [ω, ← Complex.exp_nat_mul]
  rw [hconj, homega k.val, homega (2 * k.val)]
  set θ : ℂ := (φ : ℂ) + 2 * (Real.pi : ℂ) * (k.val : ℂ) / 3 with hθ
  have e1 : (t : ℂ) * Complex.exp ((φ : ℂ) * Complex.I)
        * Complex.exp ((k.val : ℂ) * (2 * (Real.pi : ℂ) * Complex.I / 3))
      = (t : ℂ) * Complex.exp (θ * Complex.I) := by
    rw [mul_assoc, ← Complex.exp_add]; congr 2; rw [hθ]; ring
  have e2 : (t : ℂ) * Complex.exp (-(φ : ℂ) * Complex.I)
        * Complex.exp (((2 * k.val : ℕ) : ℂ) * (2 * (Real.pi : ℂ) * Complex.I / 3))
      = (t : ℂ) * Complex.exp (-(θ * Complex.I)) := by
    rw [mul_assoc, ← Complex.exp_add]
    rw [show -(φ : ℂ) * Complex.I + ((2 * k.val : ℕ) : ℂ) * (2 * (Real.pi : ℂ) * Complex.I / 3)
          = -(θ * Complex.I) + (k.val : ℂ) * (2 * (Real.pi : ℂ) * Complex.I) by
        push_cast; rw [hθ]; ring]
    rw [Complex.exp_add, Complex.exp_nat_mul, Complex.exp_two_pi_mul_I, one_pow, mul_one]
  rw [e1, e2, ← mul_add]
  rw [show Complex.exp (θ * Complex.I) + Complex.exp (-(θ * Complex.I)) = 2 * Complex.cos θ by
      rw [Complex.cos]; ring_nf]
  rw [hθ, show ((φ : ℂ) + 2 * (Real.pi : ℂ) * (k.val : ℂ) / 3)
        = ((φ + 2 * Real.pi * (k.val : ℝ) / 3 : ℝ) : ℂ) by push_cast; ring]
  rw [← Complex.ofReal_cos]; push_cast; ring

/-- **The eigenvalues are the Brannen amplitudes.** With `ξ = t·e^{iφ}`, the
    complex eigenvalue `lam (↑a) ξ k` equals the real Brannen amplitude
    `BrannenKernel.s a t φ k`.  Composed with `BrannenKernel.Q_value`, this shows
    the *matrix* `M` (not merely the abstract amplitude function) realises the
    Koide ratio `Q = (1 + 2 t²)/3`, and with `koide_iff_constraint` the value
    `2/3` at the constraint `t² = 1/2`. -/
theorem lam_eq_brannen (a t φ : ℝ) (k : Fin 3) :
    lam (a : ℂ) ((t : ℂ) * Complex.exp ((φ : ℂ) * Complex.I)) k
      = ((SCPv59.BrannenKernel.s a t φ k : ℝ) : ℂ) := by
  rw [lam, add_assoc, xi_sum]
  rw [show ((a : ℂ) * (1 + ↑(2 * t * Real.cos (φ + 2 * Real.pi * (k.val : ℝ) / 3))))
        = ((a * (1 + 2 * t * Real.cos (φ + 2 * Real.pi * (k.val : ℝ) / 3)) : ℝ) : ℂ) by
      push_cast; ring]
  norm_cast
  fin_cases k <;> (simp only [SCPv59.BrannenKernel.s]; ring_nf)

end SCPv59.KernelEigenvalues
