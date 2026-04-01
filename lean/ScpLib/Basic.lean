/-
  ScpLib.Basic — Core definitions for the SCP field theory formalization

  Defines:
  - R: axiomatic real number type
  - FieldVec: a 3-component vector (Fin 3 → R)
  - Triple product P(φ) = φ₀ * φ₁ * φ₂
  - Partial derivatives of P
  - Key identities (cyclic symmetry, derivative formula)

  We axiomatize R with instances for Add, Mul, etc. and state
  algebraic identities as axioms. Without Mathlib, we cannot use
  `ring` or `linarith`, so algebraic steps use `sorry` — the
  STRUCTURE of the proofs is what matters here.
-/

noncomputable section

/-! ## Axiomatic real numbers

We define R as an opaque type with all the algebraic operations
and key properties axiomatized.
-/

axiom R : Type
axiom R.add : R → R → R
axiom R.mul : R → R → R
axiom R.neg : R → R
axiom R.sub : R → R → R
axiom R.div : R → R → R
axiom R.zero : R
axiom R.one : R
axiom R.le : R → R → Prop
axiom R.lt : R → R → Prop
axiom R.ofNat : Nat → R
axiom R.pow : R → Nat → R

instance : Add R := ⟨R.add⟩
instance : Mul R := ⟨R.mul⟩
instance : Neg R := ⟨R.neg⟩
instance : Sub R := ⟨R.sub⟩
instance : Div R := ⟨R.div⟩
instance : Zero R := ⟨R.zero⟩
instance : One R := ⟨R.one⟩
instance : LE R := ⟨R.le⟩
instance : LT R := ⟨R.lt⟩
instance (n : Nat) : OfNat R n := ⟨R.ofNat n⟩
instance : HPow R Nat R := ⟨R.pow⟩

-- Core algebraic axioms
axiom R.add_comm (a b : R) : a + b = b + a
axiom R.add_assoc (a b c : R) : (a + b) + c = a + (b + c)
axiom R.add_zero (a : R) : a + 0 = a
axiom R.zero_add (a : R) : 0 + a = a
axiom R.add_neg_cancel (a : R) : a + (-a) = 0
axiom R.mul_comm (a b : R) : a * b = b * a
axiom R.mul_assoc (a b c : R) : (a * b) * c = a * (b * c)
axiom R.mul_one (a : R) : a * 1 = a
axiom R.one_mul (a : R) : 1 * a = a
axiom R.mul_zero (a : R) : a * 0 = 0
axiom R.zero_mul (a : R) : 0 * a = 0
axiom R.left_distrib (a b c : R) : a * (b + c) = a * b + a * c
axiom R.right_distrib (a b c : R) : (a + b) * c = a * c + b * c
axiom R.sub_def (a b : R) : a - b = a + (-b)
axiom R.neg_mul (a b : R) : (-a) * b = -(a * b)
axiom R.mul_neg (a b : R) : a * (-b) = -(a * b)
axiom R.neg_neg (a : R) : -(-a) = a
axiom R.neg_zero : (-(0 : R)) = (0 : R)

-- Power axioms
axiom R.pow_zero (a : R) : a ^ (0 : Nat) = 1
axiom R.pow_one (a : R) : a ^ (1 : Nat) = a
axiom R.pow_two (a : R) : a ^ (2 : Nat) = a * a

-- Order axioms
axiom R.le_refl (a : R) : a ≤ a
axiom R.sq_nonneg (a : R) : (0 : R) ≤ a * a
axiom R.add_nonneg (a b : R) : (0 : R) ≤ a → (0 : R) ≤ b → (0 : R) ≤ a + b
axiom R.mul_nonneg (a b : R) : (0 : R) ≤ a → (0 : R) ≤ b → (0 : R) ≤ a * b
axiom R.mul_pos (a b : R) : (0 : R) < a → (0 : R) < b → (0 : R) < a * b
axiom R.neg_pos_of_neg (a : R) : a < 0 → 0 < (-a)
axiom R.neg_neg_of_pos (a : R) : (0 : R) < a → (-a) < 0
axiom R.lt_irrefl (a : R) : ¬(a < a)

-- OfNat compatibility
axiom R.ofNat_zero_eq : R.ofNat 0 = R.zero
axiom R.ofNat_one_eq : R.ofNat 1 = R.one
axiom R.two_pos : (0 : R) < (2 : R)
axiom R.two_ne_zero : (2 : R) ≠ (0 : R)
-- Numeral decomposition (needed to connect OfNat literals to algebra)
axiom R.two_eq : (2 : R) = (1 : R) + (1 : R)
axiom R.three_eq : (3 : R) = (2 : R) + (1 : R)

-- Cancellation
axiom R.mul_left_cancel (a b c : R) : a ≠ 0 → a * b = a * c → b = c
axiom R.add_left_cancel (a b c : R) : a + b = a + c → b = c

-- Division
axiom R.div_def (a b : R) : a / b = a * R.div 1 b  -- simplified
axiom R.mul_div_cancel (a : R) (b : R) : b ≠ 0 → a / b * b = a
axiom R.div_nonneg (a b : R) : (0 : R) ≤ a → (0 : R) < b → (0 : R) ≤ a / b
-- For our purposes: a/2 - a/2 = 0
axiom R.sub_self (a : R) : a - a = 0
axiom R.half_sub_half (a : R) : a / 2 - a / 2 = 0

-- Helper: a + a = 2 * a
theorem R.two_mul (a : R) : (2 : R) * a = a + a := by
  rw [R.two_eq, R.right_distrib, R.one_mul]

-- Helper: a + a + a = 3 * a
theorem R.three_mul (a : R) : (3 : R) * a = a + a + a := by
  rw [R.three_eq, R.right_distrib, R.two_mul, R.one_mul]

namespace ScpLib

/-! ## Field vector type -/

/-- A 3-component field vector. Represents (φ₀, φ₁, φ₂) or (θ₀, θ₁, θ₂) at a point. -/
abbrev FieldVec := Fin 3 → R

/-- Construct a FieldVec from three components. -/
def mkVec (x y z : R) : FieldVec
  | ⟨0, _⟩ => x
  | ⟨1, _⟩ => y
  | ⟨2, _⟩ => z

/-! ## Dot product -/

/-- Dot product of two field vectors: Σ_a u_a v_a -/
def dot (u v : FieldVec) : R :=
  u 0 * v 0 + u 1 * v 1 + u 2 * v 2

/-- Squared norm: |v|² = v · v -/
def normSq (v : FieldVec) : R := dot v v

/-- normSq expanded. -/
theorem normSq_expand (v : FieldVec) :
    normSq v = v 0 * v 0 + v 1 * v 1 + v 2 * v 2 := rfl

/-- normSq is non-negative (sum of squares). -/
theorem normSq_nonneg (v : FieldVec) : (0 : R) ≤ normSq v := by
  unfold normSq dot
  exact R.add_nonneg _ _
    (R.add_nonneg _ _ (R.sq_nonneg (v 0)) (R.sq_nonneg (v 1)))
    (R.sq_nonneg (v 2))

/-! ## Triple product -/

/-- The triple product P(φ) = φ₀ * φ₁ * φ₂. -/
def tripleProduct (φ : FieldVec) : R :=
  φ 0 * φ 1 * φ 2

/-- Partial derivative of P w.r.t. component a. -/
def dTripleProduct (φ : FieldVec) (a : Fin 3) : R :=
  match a with
  | ⟨0, _⟩ => φ 1 * φ 2
  | ⟨1, _⟩ => φ 0 * φ 2
  | ⟨2, _⟩ => φ 0 * φ 1

/-! ## Triple product identities -/

/-- P is symmetric under cyclic permutation: P(φ₀,φ₁,φ₂) = P(φ₁,φ₂,φ₀).
    Proof: φ₀*φ₁*φ₂ = φ₁*φ₂*φ₀ by commutativity and associativity. -/
theorem tripleProduct_cyclic (φ : FieldVec) :
    tripleProduct φ =
    tripleProduct (mkVec (φ 1) (φ 2) (φ 0)) := by
  unfold tripleProduct mkVec
  -- (φ0 * φ1) * φ2 = (φ1 * φ2) * φ0 -- by comm/assoc of R
  calc (φ 0 * φ 1) * φ 2
      = φ 0 * (φ 1 * φ 2) := R.mul_assoc _ _ _
    _ = (φ 1 * φ 2) * φ 0 := R.mul_comm _ _

/-- P is symmetric under the other cyclic permutation. -/
theorem tripleProduct_cyclic' (φ : FieldVec) :
    tripleProduct φ =
    tripleProduct (mkVec (φ 2) (φ 0) (φ 1)) := by
  unfold tripleProduct mkVec
  -- (φ0 * φ1) * φ2 = (φ2 * φ0) * φ1
  calc (φ 0 * φ 1) * φ 2
      = φ 2 * (φ 0 * φ 1) := R.mul_comm _ _
    _ = (φ 2 * φ 0) * φ 1 := (R.mul_assoc _ _ _).symm

/-- ∂P/∂φ₀ = φ₁ * φ₂ -/
theorem dTripleProduct_0 (φ : FieldVec) :
    dTripleProduct φ 0 = φ 1 * φ 2 := rfl

/-- ∂P/∂φ₁ = φ₀ * φ₂ -/
theorem dTripleProduct_1 (φ : FieldVec) :
    dTripleProduct φ 1 = φ 0 * φ 2 := rfl

/-- ∂P/∂φ₂ = φ₀ * φ₁ -/
theorem dTripleProduct_2 (φ : FieldVec) :
    dTripleProduct φ 2 = φ 0 * φ 1 := rfl

/-- Euler's theorem: Σ_a (∂P/∂φ_a) * φ_a = 3P.
    Each term is a rearrangement of φ₀*φ₁*φ₂, so the sum is 3P. -/
theorem euler_triple_product (φ : FieldVec) :
    dTripleProduct φ 0 * φ 0 +
    dTripleProduct φ 1 * φ 1 +
    dTripleProduct φ 2 * φ 2 =
    (3 : R) * tripleProduct φ := by
  unfold dTripleProduct tripleProduct
  -- LHS: (φ1*φ2)*φ0 + (φ0*φ2)*φ1 + (φ0*φ1)*φ2
  -- Each term equals (φ0*φ1)*φ2 by comm/assoc, so sum = 3*((φ0*φ1)*φ2)
  have h1 : (φ 1 * φ 2) * φ 0 = (φ 0 * φ 1) * φ 2 := by
    calc (φ 1 * φ 2) * φ 0
        = φ 0 * (φ 1 * φ 2) := R.mul_comm _ _
      _ = (φ 0 * φ 1) * φ 2 := (R.mul_assoc _ _ _).symm
  have h2 : (φ 0 * φ 2) * φ 1 = (φ 0 * φ 1) * φ 2 := by
    calc (φ 0 * φ 2) * φ 1
        = φ 0 * (φ 2 * φ 1) := R.mul_assoc _ _ _
      _ = φ 0 * (φ 1 * φ 2) := by rw [R.mul_comm (φ 2) (φ 1)]
      _ = (φ 0 * φ 1) * φ 2 := (R.mul_assoc _ _ _).symm
  rw [h1, h2]
  rw [R.three_mul]

/-- Triple product of zero field is zero. -/
theorem tripleProduct_zero :
    tripleProduct (fun (_ : Fin 3) => (0 : R)) = 0 := by
  simp [tripleProduct, R.mul_zero]

/-- dTripleProduct of zero field is zero for all components. -/
theorem dTripleProduct_zero (a : Fin 3) :
    dTripleProduct (fun (_ : Fin 3) => (0 : R)) a = 0 := by
  match a with
  | ⟨0, _⟩ => simp [dTripleProduct, R.zero_mul]
  | ⟨1, _⟩ => simp [dTripleProduct, R.zero_mul]
  | ⟨2, _⟩ => simp [dTripleProduct, R.zero_mul]

end ScpLib
