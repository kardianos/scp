/-
  V44.Equations — The V44 6-field Cosserat field equations

  The V44 system:
    φ̈_a = ∇²φ_a - m²φ_a - (dV/dP)(dP/dφ_a) + η curl(θ)_a
    θ̈_a = ∇²θ_a - m_θ²θ_a + η curl(φ)_a
-/

import ScpLib.Basic
import ScpLib.VectorCalc
import ScpLib.EulerLagrange
import ScpLib.Energy

noncomputable section

namespace V44

open ScpLib

/-! ## Force expressions -/

/-- Total force on φ_a in V44:
    F_φ_a = lap_φ_a - m²φ_a - (dV/dP)(dP/dφ_a) + η curl(θ)_a -/
def phiForce (p : SCPParams) (φ : FieldVec) (lapPhi curlTheta : FieldVec)
    (a : Fin 3) : R :=
  lapPhi a + massForce p.mSq φ a + potentialForce p.μ p.κ φ a + p.η * curlTheta a

/-- Total force on θ_a in V44:
    F_θ_a = lap_θ_a - m_θ²θ_a + η curl(φ)_a -/
def thetaForce (p : SCPParams) (θ : FieldVec) (lapTheta curlPhi : FieldVec)
    (a : Fin 3) : R :=
  lapTheta a + massForce p.mThetaSq θ a + p.η * curlPhi a

/-! ## Static energy density -/

/-- V44 static energy density:
    E = (1/2)|∇φ|² + (1/2)|∇θ|² + (m²/2)|φ|² + (m_θ²/2)|θ|² + V(P) - η φ·curl(θ) -/
def staticEnergy (p : SCPParams) (φ θ : FieldVec)
    (gradPhiSq gradThetaSq : R) (curlTheta : FieldVec) : R :=
  gradPhiSq / (2 : R) + gradThetaSq / (2 : R) +
  energyMassPhi p.mSq φ +
  energyMassTheta p.mThetaSq θ +
  energyPotential p.μ p.κ φ +
  energyCurlCoupling p.η φ curlTheta

/-! ## Properties -/

/-- At η = 0, φ force reduces to Klein-Gordon + potential. -/
theorem phiForce_zero_eta (p : SCPParams) (hη : p.η = (0 : R)) (φ : FieldVec)
    (lapPhi curlTheta : FieldVec) (a : Fin 3) :
    phiForce p φ lapPhi curlTheta a =
    lapPhi a + massForce p.mSq φ a + potentialForce p.μ p.κ φ a + (0 : R) := by
  unfold phiForce
  rw [hη]
  rw [R.zero_mul]

/-- At η = 0, θ force reduces to Klein-Gordon. -/
theorem thetaForce_zero_eta (p : SCPParams) (hη : p.η = (0 : R)) (θ : FieldVec)
    (lapTheta curlPhi : FieldVec) (a : Fin 3) :
    thetaForce p θ lapTheta curlPhi a =
    lapTheta a + massForce p.mThetaSq θ a + (0 : R) := by
  unfold thetaForce
  rw [hη, R.zero_mul]

/-- The trivial vacuum (φ = 0, θ = 0) is stationary. -/
theorem vacuum_stationary_phi (p : SCPParams) (a : Fin 3) :
    phiForce p (fun _ => (0 : R)) (fun _ => (0 : R)) (fun _ => (0 : R)) a = 0 := by
  unfold phiForce
  rw [massForce_zero, potentialForce_zero]
  simp [R.add_zero, R.mul_zero]

theorem vacuum_stationary_theta (p : SCPParams) (a : Fin 3) :
    thetaForce p (fun _ => (0 : R)) (fun _ => (0 : R)) (fun _ => (0 : R)) a = 0 := by
  unfold thetaForce
  rw [massForce_zero]
  simp [R.add_zero, R.mul_zero]

end V44
