/-
  ScpLib.Energy — Energy functional definitions for SCP field theory

  Defines all energy density contributions:
  - V_base(P): saturating potential
  - E_strain(φ,θ): Cosserat strain energy
  - E_chiral(φ): helicity-weighted coupling
  - E_curl: curl coupling energy
-/

import ScpLib.Basic
import ScpLib.VectorCalc
import ScpLib.EulerLagrange

noncomputable section

namespace ScpLib

/-! ## Energy density terms -/

/-- Potential energy density: V_base(P). -/
def energyPotential (μ κ : R) (φ : FieldVec) : R :=
  potential μ κ (tripleProduct φ)

/-- Mass energy density for φ: E_mass = (m²/2)|φ|². -/
def energyMassPhi (mSq : R) (φ : FieldVec) : R :=
  mSq / (2 : R) * normSq φ

/-- Mass energy density for θ: E_mass_θ = (m_θ²/2)|θ|². -/
def energyMassTheta (mThetaSq : R) (θ : FieldVec) : R :=
  mThetaSq / (2 : R) * normSq θ

/-! ## Cosserat strain energy -/

/-- The mismatch vector: M_a = curl(φ)_a / 2 - θ_a.
    Central quantity in Cosserat strain. When M = 0, the microspin θ
    tracks the macroscopic vorticity curl(φ)/2 perfectly. -/
def mismatch (curlPhi θ : FieldVec) : FieldVec := fun a =>
  curlPhi a / (2 : R) - θ a

/-- Cosserat strain energy density:
    E_strain = α Σ_a (curl(φ)_a / 2 - θ_a)² = α |M|² -/
def energyStrain (α : R) (curlPhi θ : FieldVec) : R :=
  α * normSq (mismatch curlPhi θ)

/-- Strain energy in terms of mismatch. -/
theorem energyStrain_eq (α : R) (curlPhi θ : FieldVec) :
    energyStrain α curlPhi θ = α * normSq (mismatch curlPhi θ) := rfl

/-- Strain energy is non-negative when α ≥ 0. -/
theorem energyStrain_nonneg (α : R) (hα : (0 : R) ≤ α) (curlPhi θ : FieldVec) :
    (0 : R) ≤ energyStrain α curlPhi θ := by
  unfold energyStrain
  exact R.mul_nonneg α (normSq (mismatch curlPhi θ)) hα (normSq_nonneg _)

/-! ## Chiral helicity energy -/

/-- Chiral helicity energy: E_h = κ_h P² φ·curl(φ). -/
def energyChiral (κh P φDotCurlPhi : R) : R :=
  κh * P ^ (2 : Nat) * φDotCurlPhi

/-! ## Curl coupling energy -/

/-- Curl coupling energy (V44): E_curl = -η φ·curl(θ). -/
def energyCurlCoupling (η : R) (φ : FieldVec) (curlTheta : FieldVec) : R :=
  -(η * dot φ curlTheta)

end ScpLib
