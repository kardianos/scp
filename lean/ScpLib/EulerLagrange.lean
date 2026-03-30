/-
  ScpLib.EulerLagrange — Euler-Lagrange variation machinery

  For a Lagrangian density L(φ, ∂_j φ):
    F_a = ∂L/∂φ_a - ∂_j(∂L/∂(∂_j φ_a))

  Sets up the formal framework and defines standard force terms.
-/

import ScpLib.Basic
import ScpLib.VectorCalc

noncomputable section

namespace ScpLib

/-! ## Derivative tensor -/

/-- The spatial derivative tensor: (∂_j φ_a) is a 3x3 matrix at each point. -/
abbrev DerivTensor := Point → Fin 3 → Fin 3 → R

/-- Abstract extraction of derivatives from a field. -/
axiom fieldDeriv : VectorField → DerivTensor

/-! ## Standard force terms -/

/-- Mass force: F_mass_a = -m² φ_a.
    From L_mass = -(m²/2) |φ|². -/
def massForce (mSq : R) (φ : FieldVec) (a : Fin 3) : R :=
  -(mSq * φ a)

/-- Mass force on zero field is zero. -/
theorem massForce_zero (mSq : R) (a : Fin 3) :
    massForce mSq (fun _ => (0 : R)) a = 0 := by
  simp [massForce, R.mul_zero, R.neg_zero]

/-! ## Potential -/

/-- The SCP potential: V(P) = (μ/2) P² / (1 + κ P²).
    Saturating potential that prevents unbounded growth. -/
def potential (μ κ P : R) : R :=
  (μ / (2 : R)) * (P ^ (2 : Nat)) / ((1 : R) + κ * P ^ (2 : Nat))

/-- Derivative of potential: dV/dP = μ P / (1 + κ P²)². -/
def dPotential (μ κ P : R) : R :=
  μ * P / ((1 : R) + κ * P ^ (2 : Nat)) ^ (2 : Nat)

/-- dV/dP formula is a tautology (definition). -/
theorem dPotential_def (μ κ P : R) :
    dPotential μ κ P = μ * P / ((1 : R) + κ * P ^ (2 : Nat)) ^ (2 : Nat) := rfl

/-- Potential force: F_V_a = -(dV/dP)(dP/dφ_a). -/
def potentialForce (μ κ : R) (φ : FieldVec) (a : Fin 3) : R :=
  -(dPotential μ κ (tripleProduct φ) * dTripleProduct φ a)

/-- Potential force at φ = 0 is zero. -/
theorem potentialForce_zero (μ κ : R) (a : Fin 3) :
    potentialForce μ κ (fun _ => (0 : R)) a = 0 := by
  unfold potentialForce dPotential
  rw [tripleProduct_zero, dTripleProduct_zero]
  simp [R.mul_zero, R.neg_zero]

/-! ## Energy parameters -/

/-- Parameters for the SCP energy functional. -/
structure SCPParams where
  mSq : R        -- φ mass squared (default 2.25)
  mThetaSq : R   -- θ mass squared (default 0, massless)
  μ : R          -- potential coupling
  κ : R          -- potential saturation
  η : R          -- curl coupling strength

end ScpLib
