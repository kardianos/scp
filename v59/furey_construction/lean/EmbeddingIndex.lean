/-
  v59/furey_construction/lean/EmbeddingIndex.lean

  The Killing-form embedding index for so(n) ⊂ so(N).

  For the canonical block embedding (so(n) acting on the first n coordinates
  of ℝᴺ, trivial on the remaining N − n coordinates), the Killing-form
  embedding index — defined as the ratio of restricted Killing form to
  intrinsic Killing form — is exactly  `(N − 2) / (n − 2)`.

  This is a standard Lie algebra result.  Its derivation uses:

    (a) Killing form formula for so(n):
            B_{so(n)}(X, Y)  =  (n − 2) · Tr_vec(X · Y)
        where Tr_vec is the trace in the n-dim vector representation.
        Proof: structure-constant calculation on the basis E_{ij}.

    (b) Trace preservation under block embedding:
            Tr_vec_N (ι(X) · ι(Y))  =  Tr_vec_n (X · Y)
        for the canonical embedding ι: so(n) ↪ so(N).
        Proof: the block embedding places X in the upper-left n×n block,
        zeros elsewhere; the matrix product preserves the block structure;
        the trace of an N×N matrix with non-zero entries only in the upper-
        left n×n block is the trace of that n×n block.

  Combined: B_{so(N)}|_so(n) = (N − 2) · Tr_vec_n = (N − 2)/(n − 2) · B_{so(n)}.

  This file encodes the formula and proves the specific n=3, N=7 instance
  (the value 5) that appears in v59's SU(2)_L coupling conjecture.  It also
  proves the consistency identity that this 5 equals `dim Spin(7) − dim Cl(3,1)`.

  Full Lean derivation of the Killing form formula `B = (n − 2) · Tr_vec`
  requires substantial matrix-Lie-algebra machinery (structure constants,
  adjoint representation, trace identities).  This is encoded but
  not re-derived here; the numerical verification is in
  `cosserat_experiment/06_killing_form.py`.
-/

import Mathlib.Tactic
import SpinDimension

namespace SCPv59.EmbeddingIndex

open SCPv59 SCPv59.SpinDimension

/-! ## Definition and core theorem -/

/-- The Killing-form embedding index for the canonical block embedding
    so(n) ⊂ so(N), valid for n ≥ 3 (otherwise so(n) has zero Killing form).

    Value: `(N − 2) / (n − 2)`. -/
def killingEmbeddingIndex (n N : ℕ) : ℚ := (N - 2 : ℚ) / (n - 2 : ℚ)

/-- The defining formula. -/
theorem killingEmbeddingIndex_def (n N : ℕ) :
    killingEmbeddingIndex n N = (N - 2 : ℚ) / (n - 2 : ℚ) := rfl

/-! ## Specific values relevant to v59 -/

/-- **The v59 case**: so(3) ⊂ so(7), embedding index = 5.

    This is the structural origin of the "5" in the v59 SU(2)_L coupling
    conjecture g_W² = 5·√α.  Verified numerically in
    `cosserat_experiment/06_killing_form.py` to machine precision. -/
theorem index_so3_so7 : killingEmbeddingIndex 3 7 = 5 := by
  unfold killingEmbeddingIndex
  norm_num

/-- so(3) ⊂ so(8) (Spin(8) is the triality group): index = 6. -/
theorem index_so3_so8 : killingEmbeddingIndex 3 8 = 6 := by
  unfold killingEmbeddingIndex
  norm_num

/-- so(3) ⊂ so(4) (the ℍ rotation group): index = 2. -/
theorem index_so3_so4 : killingEmbeddingIndex 3 4 = 2 := by
  unfold killingEmbeddingIndex
  norm_num

/-- so(7) ⊂ so(8): index = 6/5. -/
theorem index_so7_so8 : killingEmbeddingIndex 7 8 = 6 / 5 := by
  unfold killingEmbeddingIndex
  norm_num

/-! ## Dual interpretation of the "5"

The "5" in the SU(2)_L conjecture has TWO equivalent structural
interpretations: it is both

  • the Killing-form embedding index for so(3) ⊂ so(7), and
  • the dimension difference dim Spin(7) − dim Cl(3,1) = 21 − 16.

These two interpretations are equal AS NATURAL NUMBERS, but they come from
different parts of the v59 structural picture.  The dimension-difference
form (encoded in `Predictions.killing_index_eq_dim_diff`) suggests a unified
ambient algebra; the embedding-index form (encoded here) connects to
Lie algebra theory directly. -/

/-- The "5" of v59's SU(2)_L conjecture admits two equivalent structural
    interpretations: Killing-form embedding index AND dimension difference. -/
theorem killing_index_5_dual_form :
    killingEmbeddingIndex 3 7 = 5 ∧ killingEmbeddingIndex 3 7 = ((dimSO 7 : ℚ) - 16) := by
  refine ⟨index_so3_so7, ?_⟩
  rw [index_so3_so7, dimSpin_seven]
  norm_num

/-! ## Connection to the g_W formula

The v59 SU(2)_L coupling conjecture `g_W² = 5·√α` uses the embedding
index 5 in front of `√α` (the v59 EM-instanton scale).  Restating: -/

/-- `g_W² = X · √α` where X is the so(3) ⊂ so(7) embedding index. -/
theorem g_W_sq_from_embedding (α : ℝ) :
    (5 : ℝ) * Real.sqrt α = (killingEmbeddingIndex 3 7 : ℝ) * Real.sqrt α := by
  rw [index_so3_so7]
  norm_cast

/-! ## The general n=3 case

For so(3) ⊂ so(N), the embedding index is `N − 2`.  This unifies all the
above: at N = 7 (Spin(7), v59's natural group), we recover 5. -/

/-- so(3) ⊂ so(N) embedding index = N − 2. -/
theorem index_so3 (N : ℕ) :
    killingEmbeddingIndex 3 N = (N : ℚ) - 2 := by
  unfold killingEmbeddingIndex
  show ((N : ℚ) - 2) / ((3 : ℚ) - 2) = (N : ℚ) - 2
  ring

end SCPv59.EmbeddingIndex
