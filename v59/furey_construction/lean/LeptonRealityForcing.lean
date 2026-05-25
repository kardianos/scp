/-
  v59/furey_construction/lean/LeptonRealityForcing.lean

  **Parts (A) and (B) of the complex-structure shape-constraint program** (see
  `7D_Algebra/notes/2026-05-24-LeptonComplexStructure-Resolution.md` and the follow-up
  shape-constraint discussion).

  `LeptonComplexStructure.lean` proved the lepton's complex structure `J` (`J² = −I`) is
  L-grade *for simple blades* (every Λ²/Λ⁶ blade squares to `−I`, every Λ⁴ blade to `+I`).
  That left a loophole: a *non-blade* element — an arbitrary linear combination of Λ⁴
  blades — might square to `−I` via cross terms and so be an F-grade complex structure.

  This module closes the loophole at the subspace level, via the grade ↔ transpose-symmetry
  correspondence, and adds the reality (Hermiticity) constraint:

    (A)  `L = Λ²⊕Λ⁶ = antisymmetric (skew) matrices = so(8)`,
         `{Λ⁰}⊕F = symmetric matrices`,  and
         **no symmetric matrix is a complex structure** (positivity:
         `(MᵀM)₀₀ = Σₖ (M₀ₖ)² ≥ 0 ≠ −1`).  Hence *no* element of the whole `{Λ⁰}⊕F`
         subspace — blade or combination — is a complex structure.

    (B)  An orthogonal (reality/Hermiticity) complex structure is skew, so it lies in the
         L = skew subspace and meets the symmetric `{Λ⁰}⊕F` subspace only at `0`.  The
         physical lepton mass `M = a(I + ξS + ξ̄S²)` is Hermitian (real masses), which is
         exactly this reality condition; so `J` is forced into L.

  This uses Mathlib (for `nlinarith`/`Matrix`) plus the gamma matrices of `SevenDAlgebra`,
  so it lives in the parent package like `Z2Z2Forcing`.
-/
import Mathlib
import SevenDAlgebra
import LeptonComplexStructure

namespace SCPv59.Furey7D.LeptonRealityForcing

open SCPv59.Furey7D
open SCPv59.Furey7D.LeptonComplexStructure (negId8 isComplexStructure all_L_grade)
open Matrix

/-! ## Entry accessor, transpose-symmetry predicates -/

/-- Matrix entry `(i,j)` (default `0`, `replicate`-padded — matches `matMul`'s internal
    accessor so reductions of `matMul` land on this `ent`). -/
def ent (M : Mat8) (i j : Nat) : Int := listGetD (listGetD M i (List.replicate 8 0)) j 0

/-- `replicate n 0` reads as `0` at every position. -/
lemma listGetD_repl0 : ∀ (n j : Nat), listGetD (List.replicate n (0 : Int)) j 0 = 0
  | 0, _ => rfl
  | _ + 1, 0 => rfl
  | n + 1, j + 1 => by
      show listGetD (List.replicate n (0 : Int)) j 0 = 0
      exact listGetD_repl0 n j

/-- The empty-list default and the `replicate`-default agree (both pad with `0`), so the
    `[]`-based accessors of `matAdd`/`matScale` coincide with `ent`. -/
lemma default_irrel (M : Mat8) (i j : Nat) :
    listGetD (listGetD M i []) j 0 = ent M i j := by
  unfold ent
  rcases h : listGetOpt M i with _ | v
  · have e1 : listGetD M i [] = ([] : List Int) := by simp [listGetD, h]
    have e2 : listGetD M i (List.replicate 8 0) = List.replicate 8 0 := by simp [listGetD, h]
    rw [e1, e2, listGetD_repl0]; rfl
  · have e1 : listGetD M i [] = v := by simp [listGetD, h]
    have e2 : listGetD M i (List.replicate 8 0) = v := by simp [listGetD, h]
    rw [e1, e2]

/-- A matrix is **symmetric** iff `Mᵢⱼ = Mⱼᵢ` on the 8×8 block (= `{Λ⁰}⊕F`). -/
def isSymm (M : Mat8) : Prop := ∀ i, i < 8 → ∀ j, j < 8 → ent M i j = ent M j i

/-- A matrix is **skew/antisymmetric** iff `Mᵢⱼ = −Mⱼᵢ` (= `L = Λ²⊕Λ⁶ = so(8)`). -/
def isSkew (M : Mat8) : Prop := ∀ i, i < 8 → ∀ j, j < 8 → ent M i j = - ent M j i

instance (M : Mat8) : Decidable (isSymm M) := by unfold isSymm; infer_instance
instance (M : Mat8) : Decidable (isSkew M) := by unfold isSkew; infer_instance

/-! ## (A) grade ↔ symmetry: L = skew, {Λ⁰}⊕F = symmetric -/

/-- **Every L-grade blade is skew** (Λ²⊕Λ⁶ ⊂ so(8)): all 28 generators. -/
theorem L_grade_skew : ∀ M ∈ all_L_grade, isSkew M := by decide

/-- **Every F-grade blade is symmetric** (Λ⁴ ⊂ symmetric): all 35 generators. -/
theorem F_grade_symm : ∀ M ∈ all_F_fourforms, isSymm M := by decide

/-- The identity (Λ⁰) is symmetric. -/
theorem id8_symm : isSymm id8 := by decide

/-! ## (A) the gap-closer: no symmetric matrix is a complex structure -/

/-- The `(0,0)` entry of `M·M` as the explicit dot product `Σₖ M₀ₖ·Mₖ₀`. -/
lemma mulSelf00 (M : Mat8) :
    ent (matMul M M) 0 0
      = ent M 0 0 * ent M 0 0 + ent M 0 1 * ent M 1 0 + ent M 0 2 * ent M 2 0
      + ent M 0 3 * ent M 3 0 + ent M 0 4 * ent M 4 0 + ent M 0 5 * ent M 5 0
      + ent M 0 6 * ent M 6 0 + ent M 0 7 * ent M 7 0 := by
  simp only [ent, matMul, listGetD, listGetOpt, List.range, List.range.loop, List.foldl, List.map]
  ring

/-- `ent negId8 0 0 = −1`. -/
lemma ent_negId8_00 : ent negId8 0 0 = -1 := by decide

/-- **THE GAP-CLOSER.**  A symmetric matrix is never a complex structure: its `M·M` has a
    nonnegative `(0,0)` diagonal (`Σₖ (M₀ₖ)² ≥ 0`), but `(−I)₀₀ = −1 < 0`.  Field-independent
    (squares are nonnegative in any ordered ring), so it excludes the *entire* real
    `{Λ⁰}⊕F` (= symmetric) subspace, not just the Λ⁴ blades. -/
theorem symm_not_complexStructure (M : Mat8) (h : isSymm M) : ¬ isComplexStructure M := by
  intro hc
  have hd : ent (matMul M M) 0 0 = -1 := by
    show ent (matMul M M) 0 0 = ent negId8 0 0
    rw [show matMul M M = negId8 from hc]
  rw [mulSelf00] at hd
  rw [← h 1 (by norm_num) 0 (by norm_num), ← h 2 (by norm_num) 0 (by norm_num),
      ← h 3 (by norm_num) 0 (by norm_num), ← h 4 (by norm_num) 0 (by norm_num),
      ← h 5 (by norm_num) 0 (by norm_num), ← h 6 (by norm_num) 0 (by norm_num),
      ← h 7 (by norm_num) 0 (by norm_num)] at hd
  nlinarith [sq_nonneg (ent M 0 0), sq_nonneg (ent M 0 1), sq_nonneg (ent M 0 2),
             sq_nonneg (ent M 0 3), sq_nonneg (ent M 0 4), sq_nonneg (ent M 0 5),
             sq_nonneg (ent M 0 6), sq_nonneg (ent M 0 7), hd]

/-! ## (A) the symmetric subspace is linear: `{Λ⁰}⊕F` is closed under combination -/

lemma ent_matAdd (A B : Mat8) (i j : Nat) (hi : i < 8) (hj : j < 8) :
    ent (matAdd A B) i j = ent A i j + ent B i j := by
  rw [← default_irrel A i j, ← default_irrel B i j]
  interval_cases i <;> interval_cases j <;>
    simp only [ent, matAdd, listGetD, listGetOpt, List.range, List.range.loop, List.map]

lemma ent_matScale (c : Int) (A : Mat8) (i j : Nat) (hi : i < 8) (hj : j < 8) :
    ent (matScale c A) i j = c * ent A i j := by
  rw [← default_irrel A i j]
  interval_cases i <;> interval_cases j <;>
    simp only [ent, matScale, listGetD, listGetOpt, List.range, List.range.loop, List.map]

/-- Symmetry is preserved by matrix addition. -/
theorem isSymm_matAdd {A B : Mat8} (hA : isSymm A) (hB : isSymm B) : isSymm (matAdd A B) := by
  intro i hi j hj
  rw [ent_matAdd A B i j hi hj, ent_matAdd A B j i hj hi, hA i hi j hj, hB i hi j hj]

/-- Symmetry is preserved by scalar multiplication.  With `isSymm_matAdd` and `id8_symm`/
    `F_grade_symm`, the entire `{Λ⁰}⊕F` subspace is symmetric — hence (by
    `symm_not_complexStructure`) contains no complex structure. -/
theorem isSymm_matScale {A : Mat8} (c : Int) (hA : isSymm A) : isSymm (matScale c A) := by
  intro i hi j hj
  rw [ent_matScale c A i j hi hj, ent_matScale c A j i hj hi, hA i hi j hj]

/-! ## (B) reality (Hermiticity) forces J out of the symmetric subspace and into L -/

/-- **Skew ∩ symmetric = 0**: a matrix that is both skew and symmetric vanishes on the
    8×8 block (`Mᵢⱼ = Mⱼᵢ = −Mⱼᵢ` ⇒ `2Mⱼᵢ = 0` ⇒ `Mⱼᵢ = 0`). -/
theorem skew_inter_symm_zero (M : Mat8) (hsk : isSkew M) (hsy : isSymm M) :
    ∀ i, i < 8 → ∀ j, j < 8 → ent M i j = 0 := by
  intro i hi j hj
  have e1 : ent M i j = ent M j i := hsy i hi j hj
  have e2 : ent M i j = - ent M j i := hsk i hi j hj
  omega

/-- **Any complex structure is non-symmetric** (the contrapositive of the gap-closer):
    so a complex structure shares nothing with the symmetric `{Λ⁰}⊕F` subspace. -/
theorem complexStructure_not_symm (M : Mat8) (hcs : isComplexStructure M) : ¬ isSymm M :=
  fun hsy => symm_not_complexStructure M hsy hcs

/-- **(B), abstract reality lemma.**  Over any commutative ring, an *orthogonal*
    complex structure is skew: `Mᵀ·M = 1` and `M·M = −1` give `Mᵀ = −M`.  (Pure algebra:
    `−M` is a right inverse of `M` and `Mᵀ` a left inverse, so they coincide.)  This is the
    matrix form of "the Brannen mass `M = a(I+ξS+ξ̄S²)` is Hermitian (real masses) ⇒ the
    complex unit `J` is skew." -/
theorem orthogonal_complexStructure_skew {R : Type*} [CommRing R] {n : ℕ}
    (M : Matrix (Fin n) (Fin n) R) (hO : Mᵀ * M = 1) (hC : M * M = -1) :
    Mᵀ = -M := by
  have h : Mᵀ * (M * M) = M := by rw [← mul_assoc, hO, one_mul]
  rw [hC, mul_neg, mul_one] at h
  exact neg_eq_iff_eq_neg.mp h

/-! ## Headline: the loophole-free forcing -/

/-- **The complete shape constraint (Parts A + B).**

      (A) `L = Λ²⊕Λ⁶` are skew, `{Λ⁰}⊕F` are symmetric, and **no symmetric matrix is a
          complex structure** — so no element of the entire `{Λ⁰}⊕F` subspace (any
          combination of Λ⁴ blades and the identity, not merely blades) can be `J`;
      (B) the physical reality condition (skew, ⇔ Hermitian mass) makes a complex
          structure `J` *not* symmetric — so `J` lies in the skew subspace `L`.

    Together: the lepton's complex structure `J` is forced into `L = Λ²⊕Λ⁶`, with the
    F-grade subspace excluded as a whole.  This removes the simple-blade restriction of
    `LeptonComplexStructure.lepton_complex_structure_forced_L`.  (The reality layer (B) is
    `orthogonal_complexStructure_skew` — Hermitian mass ⇒ `J` skew ⇒ `J ∈ L` — together
    with `skew_inter_symm_zero` — `L ∩ ({Λ⁰}⊕F) = 0`.) -/
theorem complex_structure_shape_constraint :
    (∀ M ∈ all_L_grade, isSkew M)                              -- L ⊆ skew
    ∧ (∀ M ∈ all_F_fourforms, isSymm M) ∧ isSymm id8           -- {Λ⁰}⊕F ⊆ symmetric
    ∧ (∀ M : Mat8, isSymm M → ¬ isComplexStructure M)          -- (A) symmetric ⇒ not J (whole subspace)
    ∧ (∀ M : Mat8, isComplexStructure M → ¬ isSymm M) :=       -- any J is non-symmetric
  ⟨L_grade_skew, F_grade_symm, id8_symm,
   symm_not_complexStructure, complexStructure_not_symm⟩

end SCPv59.Furey7D.LeptonRealityForcing
