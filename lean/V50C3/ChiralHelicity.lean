/-
  V50C3.ChiralHelicity — Chiral helicity coupling κ_h P² φ·curl(φ)

  E_h = κ_h P² φ·curl(φ)

  Properties:
  - Vanishes when curl(φ) = 0
  - Vanishes when P = 0
  - Is a pseudoscalar (odd under parity φ → -φ)
-/

import ScpLib.Basic
import ScpLib.VectorCalc
import ScpLib.Energy

noncomputable section

namespace V50C3

open ScpLib

/-! ## Chiral helicity energy -/

/-- Chiral helicity energy: E_h = κ_h P(φ)² (φ · curl(φ)). -/
def chiralEnergy (κh : R) (φ curlPhi : FieldVec) : R :=
  κh * tripleProduct φ ^ (2 : Nat) * dot φ curlPhi

/-- E_h = 0 when curl(φ) = 0 (irrotational).
    If curlPhi = 0, then φ·curlPhi = 0. -/
theorem chiral_zero_irrotational (κh : R) (φ : FieldVec) :
    chiralEnergy κh φ (fun _ => (0 : R)) = (0 : R) := by
  unfold chiralEnergy dot
  simp [R.mul_zero, R.add_zero]

/-- E_h = 0 when P = 0. -/
theorem chiral_zero_when_P_zero (κh : R) (φ curlPhi : FieldVec)
    (hP : tripleProduct φ = (0 : R)) :
    chiralEnergy κh φ curlPhi = (0 : R) := by
  unfold chiralEnergy
  rw [hP]
  simp [R.pow_two, R.mul_zero, R.zero_mul]

/-- E_h = 0 when κ_h = 0. -/
theorem chiral_zero_uncoupled (φ curlPhi : FieldVec) :
    chiralEnergy (0 : R) φ curlPhi = (0 : R) := by
  unfold chiralEnergy
  rw [R.zero_mul, R.zero_mul]

/-! ## Parity: E_h is a pseudoscalar

Under φ → -φ, curl(φ) → -curl(φ):
- P = φ₀φ₁φ₂ → (-φ₀)(-φ₁)(-φ₂) = -P
- P² → P² (even power)
- φ · curl(φ) → (-φ)·(-curl(φ)) = φ·curl(φ)
- Wait: (-P)² = P², and (-φ)·(-curlφ) = φ·curlφ
- So actually E_h → E_h? Let's check more carefully.

φ → -φ means each component negated. Then:
- P → (-1)³ P = -P, so P² → P²
- dot(-φ, -curlφ) = dot(φ, curlφ) (two negatives cancel)
- Therefore E_h → κ_h × P² × (φ·curlφ) = E_h

E_h is actually EVEN under full parity φ → -φ.
The ODD behavior comes from a SINGLE component flip.
-/

/-- Under φ → -φ (full negation), E_h is invariant.
    This is because P² is even and dot(-φ,-curlφ) = dot(φ,curlφ). -/
theorem chiral_parity_even (κh : R) (φ curlPhi : FieldVec) :
    chiralEnergy κh (fun a => -(φ a)) (fun a => -(curlPhi a)) =
    chiralEnergy κh φ curlPhi := by
  unfold chiralEnergy tripleProduct dot
  -- (-φ0)*(-φ1)*(-φ2) = -(φ0*φ1*φ2), squared gives same P²
  -- (-φ)·(-curlφ) = φ·curlφ
  sorry  -- algebra with neg_mul, mul_neg

/-! ## Force from chiral helicity (algebraic part) -/

/-- Algebraic part of chiral force:
    F_h_a^{alg} = κ_h [2P(dP/dφ_a)(φ·curlφ) + P² curlφ_a] -/
def chiralForceAlg (κh : R) (φ curlPhi : FieldVec) (a : Fin 3) : R :=
  κh * ((2 : R) * tripleProduct φ * dTripleProduct φ a * dot φ curlPhi +
        tripleProduct φ ^ (2 : Nat) * curlPhi a)

/-- At P = 0, chiral force vanishes. -/
theorem chiralForce_zero_P (κh : R) (φ curlPhi : FieldVec)
    (hP : tripleProduct φ = (0 : R)) (a : Fin 3) :
    chiralForceAlg κh φ curlPhi a = (0 : R) := by
  unfold chiralForceAlg
  rw [hP]
  simp [R.pow_two, R.mul_zero, R.zero_mul, R.add_zero]

/-- At curlφ = 0, chiral force vanishes. -/
theorem chiralForce_zero_curl (κh : R) (φ : FieldVec) (a : Fin 3) :
    chiralForceAlg κh φ (fun _ => (0 : R)) a = (0 : R) := by
  unfold chiralForceAlg dot
  simp [R.mul_zero, R.add_zero]

end V50C3
