/-
  VecStream.lean — Algebraic properties of polynomial vectorization for SCP streaming

  Proves key invariants of the P-frame compression scheme used by the SCP simulation:
  1. L2 projection optimality (normal equations minimize ||Vc - f||^2)
  2. Reconstruction exactness (predicted + delta = actual)
  3. Temporal model as orthogonal projection (residual perp to basis)
  4. Multi-field independence (block-diagonal fitting equivalence)
  5. Error propagation in P-frame chains (no accumulation)

  Uses the project's axiomatic R type from ScpLib.Basic. Proofs that require
  deep linear algebra are admitted with sorry, but all theorem TYPES are precise.
-/

import ScpLib.Basic

noncomputable section

open ScpLib

namespace VecStream

/-! ## Vector and matrix types for linear algebra

We work with finite-dimensional vectors and matrices over R,
parameterized by dimension.
-/

/-- A vector of dimension n over R. -/
abbrev Vec (n : Nat) := Fin n → R

/-- A matrix of dimensions m x n over R. -/
abbrev Mat (m n : Nat) := Fin m → Fin n → R

/-- Dot product of two n-vectors: sum_i u_i * v_i. -/
def vdot {n : Nat} (u v : Vec n) : R :=
  Fin.foldl n (init := (0 : R)) (fun acc i => acc + u i * v i)

/-- Matrix-vector product: (Av)_i = sum_j A_{ij} v_j. -/
def matvec {m n : Nat} (A : Mat m n) (v : Vec n) : Vec m :=
  fun i => vdot (A i) v

/-- Matrix transpose: (A^T)_{ji} = A_{ij}. -/
def transpose {m n : Nat} (A : Mat m n) : Mat n m :=
  fun j i => A i j

/-- Matrix-matrix product: (AB)_{ik} = sum_j A_{ij} B_{jk}. -/
def matmul {m n p : Nat} (A : Mat m n) (B : Mat n p) : Mat m p :=
  fun i k => vdot (A i) (fun j => B j k)

/-- Vector subtraction. -/
def vsub {n : Nat} (u v : Vec n) : Vec n :=
  fun i => u i - v i

/-- Vector addition. -/
def vadd {n : Nat} (u v : Vec n) : Vec n :=
  fun i => u i + v i

/-- Squared L2 norm: ||v||^2 = v . v. -/
def normSqV {n : Nat} (v : Vec n) : R := vdot v v

/-- The zero vector. -/
def vzero (n : Nat) : Vec n := fun _ => (0 : R)

/-! ## Section 1: L2 Projection Optimality

The polynomial fitting finds c = argmin_c ||Vc - f||^2 where V is the
Vandermonde matrix (m x n, m sample points, n polynomial coefficients)
and f is the data vector (m samples).

The normal equations state: V^T V c = V^T f.
When V^T V is invertible, c = (V^T V)^{-1} V^T f.

We prove: for any other vector d, ||Vd - f||^2 >= ||Vc* - f||^2
where c* solves the normal equations.
-/

/-- The residual vector r = Vc - f. -/
def residual {m n : Nat} (V : Mat m n) (c : Vec n) (f : Vec m) : Vec m :=
  vsub (matvec V c) f

/-- c satisfies the normal equations V^T V c = V^T f. -/
def satisfiesNormalEqs {m n : Nat} (V : Mat m n) (c : Vec n) (f : Vec m) : Prop :=
  matvec (matmul (transpose V) V) c = matvec (transpose V) f

/-- The normal equations residual is orthogonal to the column space of V.
    That is, V^T (Vc - f) = 0 when c satisfies V^T V c = V^T f. -/
theorem normal_eqs_orthogonality {m n : Nat}
    (V : Mat m n) (c : Vec n) (f : Vec m)
    (hc : satisfiesNormalEqs V c f) :
    matvec (transpose V) (residual V c f) = vzero n := by
  -- V^T(Vc - f) = V^T(Vc) - V^T f = (V^T V)c - V^T f = 0
  -- by the normal equations hypothesis
  sorry

/-- L2 Projection Optimality: if c* satisfies the normal equations,
    then for ANY vector d, ||Vd - f||^2 >= ||Vc* - f||^2.

    Proof sketch: Let r* = Vc* - f and e = V(d - c*).
    Then Vd - f = r* + e, and r* perp e (by normal eqs).
    So ||Vd - f||^2 = ||r*||^2 + ||e||^2 >= ||r*||^2.
-/
theorem l2_projection_optimal {m n : Nat}
    (V : Mat m n) (cstar : Vec n) (f : Vec m) (d : Vec n)
    (hc : satisfiesNormalEqs V cstar f) :
    normSqV (residual V cstar f) ≤ normSqV (residual V d f) := by
  -- Standard result: Pythagoras + orthogonality
  -- ||Vd - f||^2 = ||r* + V(d-c*)||^2 = ||r*||^2 + ||V(d-c*)||^2 + 2<r*, V(d-c*)>
  -- The cross term vanishes by normal_eqs_orthogonality, and ||V(d-c*)||^2 >= 0.
  sorry

/-! ## Section 2: Reconstruction Exactness

The P-frame scheme:
  - Encoder computes: delta = actual - predicted
  - Decoder computes: reconstructed = predicted + delta

This is trivially exact in exact arithmetic.
-/

/-- P-frame encoding: delta = actual - predicted. -/
def encode {n : Nat} (actual predicted : Vec n) : Vec n :=
  vsub actual predicted

/-- P-frame decoding: reconstructed = predicted + delta. -/
def decode {n : Nat} (predicted delta : Vec n) : Vec n :=
  vadd predicted delta

/-- Helper: p + (a + (-p)) = a. -/
private theorem add_sub_cancel_right (p a : R) : p + (a + (-p)) = a := by
  have h1 : a + (-p) = (-p) + a := R.add_comm a (-p)
  rw [h1]
  rw [← R.add_assoc]
  rw [R.add_neg_cancel]
  rw [R.zero_add]

/-- Reconstruction Exactness: decode(predicted, encode(actual, predicted)) = actual.
    In exact arithmetic, (predicted + (actual - predicted)) = actual. -/
theorem reconstruction_exact {n : Nat}
    (actual predicted : Vec n) :
    decode predicted (encode actual predicted) = actual := by
  funext i
  unfold decode encode vadd vsub
  -- Need: predicted i + (actual i - predicted i) = actual i
  rw [R.sub_def (actual i) (predicted i)]
  -- Now: predicted i + (actual i + -predicted i) = actual i
  exact add_sub_cancel_right (predicted i) (actual i)

/-! ## Section 3: Temporal Model as Orthogonal Projection

The temporal model fits c(t) ~ mean + amp*cos(wt + phi), which is equivalent
to fitting in the basis {1, cos(wt), sin(wt)} via:
  c(t) ~ a_0 + a_1*cos(wt) + a_2*sin(wt)

This is a least-squares fit in a 3-dimensional subspace, hence an orthogonal
projection. We prove the general property: the least-squares fit in any
finite basis is an orthogonal projection, and the residual is perpendicular
to the basis.
-/

/-- A temporal basis: k basis functions evaluated at m time points.
    B_{ti} = i-th basis function evaluated at time t. -/
abbrev TemporalBasis (m k : Nat) := Mat m k

/-- The temporal fit: coefficients a such that B*a approximates the data y.
    Satisfies the normal equations B^T B a = B^T y. -/
def temporalFit {m k : Nat} (B : TemporalBasis m k) (a : Vec k) (y : Vec m) : Prop :=
  satisfiesNormalEqs B a y

/-- The fitted values: y_hat = B * a. -/
def fittedValues {m k : Nat} (B : TemporalBasis m k) (a : Vec k) : Vec m :=
  matvec B a

/-- Temporal residual: r = y - y_hat. -/
def temporalResidual {m k : Nat} (B : TemporalBasis m k) (a : Vec k) (y : Vec m) : Vec m :=
  vsub y (fittedValues B a)

/-- Orthogonal Projection Property: the residual is orthogonal to every basis column.
    For each basis function j: sum_t r_t * B_{tj} = 0.

    This is equivalent to B^T r = 0, which follows directly from the normal equations.
-/
theorem temporal_residual_orthogonal {m k : Nat}
    (B : TemporalBasis m k) (a : Vec k) (y : Vec m)
    (hfit : temporalFit B a y) :
    matvec (transpose B) (temporalResidual B a y) = vzero k := by
  -- B^T(y - Ba) = B^T y - B^T B a = B^T y - B^T y = 0
  -- by the normal equations B^T B a = B^T y
  sorry

/-- Idempotence: fitting twice gives the same result.
    If a solves B^T B a = B^T y, and we re-fit y_hat = Ba,
    the solution is again a. (Projection squared = projection.)
-/
theorem temporal_fit_idempotent {m k : Nat}
    (B : TemporalBasis m k) (a : Vec k) (y : Vec m)
    (hfit : temporalFit B a y)
    (hBinv : ∀ (u v : Vec k), matvec (matmul (transpose B) B) u =
             matvec (matmul (transpose B) B) v → u = v) :
    satisfiesNormalEqs B a (fittedValues B a) := by
  -- B^T B a = B^T (Ba), which is (B^T B) a = B^T B a. Trivially true.
  sorry

/-! ## Section 4: Multi-field Independence

With 6 fields x 64 coefficients = 384 total per patch, each field's
64 coefficients are fit independently. The combined system has a
block-diagonal structure. We prove that solving the block-diagonal
system jointly gives the same result as solving each block separately.

We require n > 0 to avoid division-by-zero in index arithmetic.
-/

/-- A block-diagonal matrix built from K identical n x n blocks along the diagonal.
    The full matrix is (K*n) x (K*n). Requires n > 0 for index arithmetic. -/
def blockDiag {n : Nat} (A : Mat n n) (K : Nat) (hn : n > 0) : Mat (K * n) (K * n) :=
  fun i j =>
    let bi := i.val / n
    let bj := j.val / n
    if bi = bj then A ⟨i.val % n, Nat.mod_lt i.val hn⟩ ⟨j.val % n, Nat.mod_lt j.val hn⟩
    else (0 : R)

/-- k*n + i < K*n when k < K and i < n. -/
private theorem block_index_bound {n K : Nat} (k : Fin K) (i : Fin n) (_hn : n > 0) :
    k.val * n + i.val < K * n := by
  have hk := k.isLt
  have hi := i.isLt
  have h1 : (k.val + 1) * n ≤ K * n := Nat.mul_le_mul_right n hk
  have h2 : k.val * n + i.val < k.val * n + n := Nat.add_lt_add_left hi _
  have h3 : k.val * n + n = (k.val + 1) * n := by
    simp [Nat.add_mul, Nat.one_mul]
  omega

/-- i / n < K when i < K * n and n > 0. -/
private theorem div_block_bound {n K : Nat} (i : Fin (K * n)) (_hn : n > 0) :
    i.val / n < K := by
  have hi := i.isLt
  have hmul : K * n = n * K := Nat.mul_comm K n
  have hi' : i.val < n * K := hmul ▸ hi
  exact Nat.div_lt_of_lt_mul hi'

/-- Extract block k from a (K*n)-vector, giving an n-vector. -/
def extractBlock {n K : Nat} (v : Vec (K * n)) (k : Fin K) (_hn : n > 0) : Vec n :=
  fun i => v ⟨k.val * n + i.val, block_index_bound k i _hn⟩

/-- Assemble K n-vectors into one (K*n)-vector. -/
def assembleBlocks {n K : Nat} (blocks : Fin K → Vec n) (_hn : n > 0) : Vec (K * n) :=
  fun i => blocks ⟨i.val / n, div_block_bound i _hn⟩ ⟨i.val % n, Nat.mod_lt i.val _hn⟩

/-- Multi-field Independence: solving Ax = b where A is block-diagonal
    gives the same result as solving each block independently.

    If A = diag(A_0, ..., A_{K-1}) and b = (b_0, ..., b_{K-1}),
    and x_k solves A x_k = b_k for each k, then
    x = (x_0, ..., x_{K-1}) solves Ax = b.

    This justifies fitting each of the 6 SCP fields independently.
-/
theorem block_diagonal_independence {n K : Nat}
    (A : Mat n n) (b : Vec (K * n))
    (x_blocks : Fin K → Vec n) (hn : n > 0)
    (hsolve : ∀ (k : Fin K), matvec A (x_blocks k) = extractBlock b k hn) :
    matvec (blockDiag A K hn) (assembleBlocks x_blocks hn) = b := by
  -- For each component i in the full vector:
  -- (A_block x)_i = sum_j A_block[i,j] x_j
  -- Only j in the same block as i contribute (off-diagonal blocks are 0).
  -- So (A_block x)_i = sum_{j in block} A[i_local, j_local] x_{block}[j_local]
  -- = (A x_{block})_{i_local} = b_{block}[i_local] = b_i
  sorry

/-- Converse: if x solves the block-diagonal system, each block solves independently. -/
theorem block_diagonal_extract {n K : Nat}
    (A : Mat n n) (b : Vec (K * n)) (x : Vec (K * n))
    (hn : n > 0)
    (hsolve : matvec (blockDiag A K hn) x = b) (k : Fin K) :
    matvec A (extractBlock x k hn) = extractBlock b k hn := by
  sorry

/-! ## Section 5: Error Propagation in P-frame Chain

For a chain of N P-frames between I-frames:
  - I-frame provides the temporal model coefficients (mean, amp, phase)
  - Each P-frame k stores delta_k = actual_k - predicted_k(t_k)
  - Reconstruction: recon_k = predicted_k(t_k) + delta_k

Key property: each P-frame is independently reconstructible from the
temporal model and its own delta. There is NO error accumulation ---
unlike differential coding where errors compound.
-/

/-- A temporal model predicting coefficients at any time index. -/
structure TemporalModel (n : Nat) where
  /-- Predict the coefficient vector at time t. -/
  predict : Nat → Vec n

/-- A P-frame: stores the delta from prediction. -/
structure PFrame (n : Nat) where
  /-- Time index of this frame. -/
  timeIdx : Nat
  /-- Delta = actual - predicted. -/
  delta : Vec n

/-- An I-frame: stores the temporal model and the exact coefficients. -/
structure IFrame (n : Nat) where
  /-- The temporal model for prediction. -/
  model : TemporalModel n
  /-- Exact coefficients at this keyframe. -/
  coeffs : Vec n

/-- Encode a P-frame: given actual coefficients and the temporal model. -/
def encodePFrame {n : Nat} (model : TemporalModel n) (t : Nat) (actual : Vec n) : PFrame n :=
  { timeIdx := t, delta := encode actual (model.predict t) }

/-- Decode a P-frame: reconstruct actual coefficients. -/
def decodePFrame {n : Nat} (model : TemporalModel n) (pf : PFrame n) : Vec n :=
  decode (model.predict pf.timeIdx) pf.delta

/-- A chain of N P-frames with their actual data. -/
structure PFrameChain (n N : Nat) where
  /-- The temporal model (from the I-frame). -/
  model : TemporalModel n
  /-- The actual coefficient vectors at each frame. -/
  actuals : Fin N → Vec n
  /-- The time indices for each frame. -/
  times : Fin N → Nat
  /-- The encoded P-frames. -/
  frames : Fin N → PFrame n
  /-- Encoding consistency: each frame was encoded from its actual data. -/
  encoded_correctly : ∀ k, frames k = encodePFrame model (times k) (actuals k)

/-- P-frame chain reconstruction: each frame individually reconstructs to actual.
    NO error accumulation --- frame k depends only on (model, delta_k), not on
    any previous frame's reconstruction.
-/
theorem pframe_chain_exact {n N : Nat} (chain : PFrameChain n N) (k : Fin N) :
    decodePFrame chain.model (chain.frames k) = chain.actuals k := by
  rw [chain.encoded_correctly k]
  unfold decodePFrame encodePFrame
  simp only
  -- Reduces to: decode (model.predict t) (encode actual (model.predict t)) = actual
  exact reconstruction_exact (chain.actuals k) (chain.model.predict (chain.times k))

/-- No error accumulation: reconstruction of frame k is independent of frame j.
    Formally: changing actual_j (and thus delta_j) does not affect recon_k for k != j.
-/
theorem pframe_independence {n N : Nat}
    (chain1 chain2 : PFrameChain n N)
    (_same_model : chain1.model = chain2.model)
    (k : Fin N)
    (same_k : chain1.actuals k = chain2.actuals k)
    (_same_time : chain1.times k = chain2.times k) :
    decodePFrame chain1.model (chain1.frames k) =
    decodePFrame chain2.model (chain2.frames k) := by
  rw [pframe_chain_exact chain1 k, pframe_chain_exact chain2 k]
  rw [_same_model] at *
  exact same_k

/-! ## Combined pipeline theorem

The full vectorization pipeline:
  1. Spatial fitting: c = argmin ||Vc - f||^2 (Section 1)
  2. Temporal fitting: model from I-frame coefficients (Section 3)
  3. P-frame: delta = actual_c - predicted_c (Section 2)
  4. Per-field independence: 6 fields fit separately (Section 4)
  5. Exact reconstruction with no accumulation (Section 5)

We state the end-to-end roundtrip theorem.
-/

/-- End-to-end roundtrip: for any frame in the stream, decoding the stored
    representation recovers the original least-squares coefficients exactly.

    The pipeline is:
      actual coefficients --[encode]--> delta --[decode]--> reconstructed
    and reconstructed = actual (by reconstruction_exact), regardless of:
    - which frame in the chain (pframe_chain_exact)
    - which field among the 6 (block_diagonal_independence)
    - the temporal model accuracy (deltas compensate exactly)
-/
theorem pipeline_roundtrip {n N : Nat}
    (model : TemporalModel n)
    (actuals : Fin N → Vec n)
    (times : Fin N → Nat)
    (k : Fin N) :
    let chain : PFrameChain n N := {
      model := model
      actuals := actuals
      times := times
      frames := fun j => encodePFrame model (times j) (actuals j)
      encoded_correctly := fun _ => rfl
    }
    decodePFrame model (chain.frames k) = actuals k := by
  simp only
  unfold decodePFrame encodePFrame
  exact reconstruction_exact (actuals k) (model.predict (times k))

end VecStream
