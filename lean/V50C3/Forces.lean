/-
  V50C3.Forces — Combined force expressions for full C3 equations

  Combines all force contributions and proves key structural results:
  - Complete φ and θ equations
  - Fixed-point theorem for θ (equilibrium at θ = curl(φ)/2)
  - Vacuum is stationary
-/

import ScpLib.Basic
import ScpLib.VectorCalc
import ScpLib.EulerLagrange
import ScpLib.Energy
import V44.Equations
import V50C3.Cosserat
import V50C3.ChiralHelicity

noncomputable section

namespace V50C3

open ScpLib

/-! ## C3 parameters -/

/-- Full C3 parameter set: V44 base + Cosserat strain + chiral. -/
structure C3Params extends SCPParams where
  α_strain : R    -- Cosserat strain coefficient
  κ_chiral : R    -- chiral helicity coupling

/-! ## Complete force expressions -/

/-- Complete φ force in C3:
    F_φ_a = ∇²φ - m²φ - (dV/dP)(dP/dφ) + η curl(θ) - α curl(M) + F_chiral -/
def phiForceC3 (p : C3Params) (φ _θ : FieldVec)
    (lapPhi curlPhi curlTheta curlM : FieldVec) (a : Fin 3) : R :=
  V44.phiForce p.toSCPParams φ lapPhi curlTheta a +
  (-(p.α_strain * curlM a)) +
  chiralForceAlg p.κ_chiral φ curlPhi a

/-- Complete θ force in C3:
    F_θ_a = ∇²θ - m_θ²θ + η curl(φ) + 2α M -/
def thetaForceC3 (p : C3Params) (θ : FieldVec)
    (lapTheta curlPhi : FieldVec) (a : Fin 3) : R :=
  V44.thetaForce p.toSCPParams θ lapTheta curlPhi a +
  thetaForceStrain p.α_strain curlPhi θ a

/-! ## Fixed-point theorem for theta

The Cosserat strain force 2α(curl(φ)/2 - θ) drives θ toward curl(φ)/2.
At θ = curl(φ)/2, the strain force vanishes: this is the fixed point.
-/

/-- The strain force has a fixed point at θ = curl(φ)/2. -/
theorem theta_fixed_point (p : C3Params) (curlPhi : FieldVec) (a : Fin 3) :
    thetaForceStrain p.α_strain curlPhi (fun b => curlPhi b / (2 : R)) a =
    (0 : R) :=
  thetaForce_zero_at_eq p.α_strain curlPhi a

/-- The strain force is linear in deviation from equilibrium:
    At θ = curl(φ)/2 + δθ: force = -2α δθ (restoring for α > 0). -/
theorem theta_linear_restoring (p : C3Params) (curlPhi δθ : FieldVec) (a : Fin 3) :
    thetaForceStrain p.α_strain curlPhi
      (fun b => curlPhi b / (2 : R) + δθ b) a =
    -((2 : R) * p.α_strain * δθ a) :=
  strain_force_linear p.α_strain curlPhi δθ a

/-- The restoring force opposes positive deviations when α > 0.
    If δθ_a > 0 and α > 0, then force_a < 0. -/
theorem theta_restoring_sign (α : R) (hα : (0 : R) < α) (δ : R) (hδ : (0 : R) < δ)
    (curlPhi : FieldVec) (a : Fin 3) :
    thetaForceStrain α curlPhi
      (fun b => curlPhi b / (2 : R) + (fun _ => δ) b) a < (0 : R) := by
  rw [strain_force_linear]
  -- Need to show -(2*α*δ) < 0, i.e., 2*α*δ > 0
  exact R.neg_neg_of_pos _ (R.mul_pos _ _ (R.mul_pos _ _ R.two_pos hα) hδ)

/-! ## Vacuum stationarity -/

/-- The trivial vacuum (φ = 0, θ = 0, all derivatives zero) is stationary. -/
theorem vacuum_phi (p : C3Params) (a : Fin 3) :
    phiForceC3 p (fun _ => (0 : R)) (fun _ => (0 : R))
      (fun _ => (0 : R)) (fun _ => (0 : R))
      (fun _ => (0 : R)) (fun _ => (0 : R)) a = (0 : R) := by
  unfold phiForceC3
  rw [V44.vacuum_stationary_phi, chiralForce_zero_curl]
  simp [R.mul_zero, R.neg_zero, R.add_zero]

theorem vacuum_theta (p : C3Params) (a : Fin 3) :
    thetaForceC3 p (fun _ => (0 : R)) (fun _ => (0 : R))
      (fun _ => (0 : R)) a = (0 : R) := by
  unfold thetaForceC3
  rw [V44.vacuum_stationary_theta]
  simp only [thetaForceStrain, mismatch]
  -- Need: 0 / 2 = 0
  have zero_div_two : (0 : R) / (2 : R) = 0 := by
    rw [R.div_def, R.zero_mul]
  simp [R.add_zero, R.sub_def, R.neg_zero, R.mul_zero, zero_div_two]

/-! ## Total C3 energy -/

/-- Total C3 static energy density. -/
def totalEnergy (p : C3Params) (φ θ : FieldVec)
    (gradPhiSq gradThetaSq : R)
    (curlPhi curlTheta : FieldVec) : R :=
  V44.staticEnergy p.toSCPParams φ θ gradPhiSq gradThetaSq curlTheta +
  energyStrain p.α_strain curlPhi θ +
  energyChiral p.κ_chiral (tripleProduct φ) (dot φ curlPhi)

/-! ## Summary of C3 force structure

The complete C3 equations are:

  φ̈_a = ∇²φ_a                    -- wave propagation
       - m²φ_a                    -- mass
       - (dV/dP)(dP/dφ_a)         -- nonlinear potential
       + η curl(θ)_a              -- V44 curl coupling
       - α curl(M)_a              -- Cosserat strain (KEY: minus sign!)
       + F_chiral_a               -- chiral helicity

  θ̈_a = ∇²θ_a                    -- wave propagation
       - m_θ²θ_a                  -- mass (usually 0)
       + η curl(φ)_a              -- V44 curl coupling
       + 2α M_a                   -- Cosserat strain restoring force

where M_a = curl(φ)_a/2 - θ_a is the Cosserat mismatch.

Key results proven in this formalization:
1. θ force from strain is 2α M_a (algebraic, definitional)
2. φ force from strain is -α curl(M)_a (EL + Levi-Civita identity)
3. Both forces vanish at θ = curl(φ)/2 (equilibrium)
4. E_strain ≥ 0 for α ≥ 0 (sum of squares)
5. Strain force is linear restoring in deviation from equilibrium
-/

end V50C3
