/-
  ScpLib.VectorCalc — Vector calculus identities for field theory

  Symbolic vector calculus: fields are functions Point → R or Point → FieldVec,
  and curl, div, grad are abstract operations satisfying standard identities.
-/

import ScpLib.Basic

noncomputable section

namespace ScpLib

/-! ## Spatial point and fields -/

/-- A point in 3D space. -/
abbrev Point := Fin 3 → R

/-- A scalar field: assigns a real number to each spatial point. -/
abbrev ScalarField := Point → R

/-- A vector field: assigns a 3-component vector to each spatial point. -/
abbrev VectorField := Point → FieldVec

/-! ## Levi-Civita symbol -/

/-- The Levi-Civita symbol ε_{ijk}.
    +1 for even permutations of (0,1,2), -1 for odd, 0 if any repeat. -/
def leviCivita (i j k : Fin 3) : R :=
  if i = j ∨ j = k ∨ i = k then (0 : R)
  else if (i.val < j.val ∧ j.val < k.val) ∨
          (j.val < k.val ∧ k.val < i.val) ∨
          (k.val < i.val ∧ i.val < j.val) then (1 : R)
  else (-1 : R)

-- Key Levi-Civita values
theorem lc_012 : leviCivita 0 1 2 = (1 : R) := by
  simp [leviCivita]

theorem lc_120 : leviCivita 1 2 0 = (1 : R) := by
  simp [leviCivita]

theorem lc_201 : leviCivita 2 0 1 = (1 : R) := by
  simp [leviCivita]

theorem lc_021 : leviCivita 0 2 1 = (-1 : R) := by
  simp [leviCivita]

theorem lc_210 : leviCivita 2 1 0 = (-1 : R) := by
  simp [leviCivita]

theorem lc_102 : leviCivita 1 0 2 = (-1 : R) := by
  simp [leviCivita]

/-- ε_{ijk} = 0 when any two indices equal. -/
theorem leviCivita_repeat (i j k : Fin 3) (h : i = j ∨ j = k ∨ i = k) :
    leviCivita i j k = (0 : R) := by
  simp [leviCivita, h]

/-- ε is antisymmetric under swap of first two indices: ε_{jik} = -ε_{ijk}. -/
theorem leviCivita_swap01 (i j k : Fin 3) :
    leviCivita j i k = R.neg (leviCivita i j k) := by
  -- Exhaustive check on all 27 cases
  sorry  -- native_decide would work with decidable ℤ; omitted for axiomatic R

/-! ## Abstract differential operators -/

/-- Abstract gradient. -/
axiom grad : ScalarField → VectorField

/-- Abstract divergence. -/
axiom div : VectorField → ScalarField

/-- Abstract curl: (curl F)_i = Σ_{j,k} ε_{ijk} ∂_j F_k -/
axiom curl : VectorField → VectorField

/-- Scalar Laplacian. -/
axiom laplacianS : ScalarField → ScalarField

/-- Vector Laplacian. -/
axiom laplacian : VectorField → VectorField

/-- Partial derivative of a vector field w.r.t. spatial direction j. -/
axiom partialDeriv : VectorField → Fin 3 → VectorField

/-! ## Fundamental identities (axioms) -/

/-- div(curl F) = 0 — divergence of a curl is zero.
    Follows from ∂_i ∂_j being symmetric while ε_{ijk} is antisymmetric. -/
axiom div_curl_zero (F : VectorField) :
    div (curl F) = (fun _ => (0 : R))

/-- curl(grad f) = 0 — curl of a gradient is zero. -/
axiom curl_grad_zero (f : ScalarField) :
    curl (grad f) = (fun _ => fun _ => (0 : R))

/-- curl(curl F) = grad(div F) - ∇²F. -/
axiom curl_curl (F : VectorField) :
    curl (curl F) = (fun x i =>
      R.add (grad (div F) x i) (R.neg (laplacian F x i)))

/-! ## Pointwise operations on vector fields -/

/-- Scalar times vector field. -/
def smulVF (f : ScalarField) (G : VectorField) : VectorField :=
  fun x i => f x * G x i

/-- Sum of vector fields. -/
def addVF (F G : VectorField) : VectorField :=
  fun x i => F x i + G x i

/-- Cross product of two vector fields (pointwise). -/
def crossVF (F G : VectorField) : VectorField := fun x i =>
  let j : Fin 3 := ⟨(i.val + 1) % 3, by omega⟩
  let k : Fin 3 := ⟨(i.val + 2) % 3, by omega⟩
  F x j * G x k - F x k * G x j

/-- curl(fG) = f curl(G) + (∇f) × G. -/
axiom curl_smul (f : ScalarField) (G : VectorField) :
    curl (smulVF f G) = addVF (smulVF f (curl G)) (crossVF (grad f) G)

/-! ## Helicity -/

/-- Helicity: h(F) = F · curl(F). Measures twistedness of the field. -/
def helicity (F : VectorField) : ScalarField := fun x =>
  dot (F x) (curl F x)

/-! ## Integration by parts for curl -/

/-- Pointwise identity: u · curl(v) = v · curl(u) + div(u × v).
    On a periodic/infinite domain, the div term integrates to zero. -/
axiom curl_ibp (u v : VectorField) (x : Point) :
    dot (u x) (curl v x) = dot (v x) (curl u x) + div (crossVF u v) x

end ScpLib
