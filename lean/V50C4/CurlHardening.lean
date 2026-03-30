/-
  V50C4/CurlHardening.lean — Curl-squared hardening term

  E_harden = (β/2) |θ|² |∇×φ|²

  Theta force: -β|∇×φ|²θ_a  (theta heavy where twist is strong)
  Phi force:   -2 curl(Q)_a  where Q_a = (β/2)|θ|²curl(φ)_a
               (Maxima-verified: derive_harden.mac)
-/

import ScpLib.Basic
import ScpLib.VectorCalc

open ScpLib

-- The hardening energy density
noncomputable def energyHarden (β : R) (theta phi_curl : FieldVec) : R :=
  (β / 2) * normSq theta * normSq phi_curl

-- Theta force from hardening: -β|∇×φ|²θ_a
noncomputable def hardenThetaForce (β : R) (curl_phi theta : FieldVec) (a : Fin 3) : R :=
  -(β * normSq curl_phi * theta a)

-- Hardening intermediate vector: Q_a = (β/2)|θ|² curl(φ)_a
noncomputable def hardenQ (β : R) (theta phi_curl : FieldVec) : FieldVec :=
  fun a => (β / 2) * normSq theta * phi_curl a

-- KEY PROPERTY: E_harden = 0 in vacuum (curl(φ) = 0)
theorem harden_zero_irrotational (β : R) (theta : FieldVec) :
    energyHarden β theta (fun _ => (0 : R)) = (0 : R) := by
  unfold energyHarden normSq dot
  simp [R.mul_zero, R.zero_add]

-- KEY PROPERTY: E_harden = 0 when θ = 0
theorem harden_zero_no_theta (β : R) (phi_curl : FieldVec) :
    energyHarden β (fun _ => (0 : R)) phi_curl = (0 : R) := by
  unfold energyHarden normSq dot
  simp [R.mul_zero, R.zero_add, R.zero_mul]

-- KEY PROPERTY: E_harden ≥ 0 when β ≥ 0
theorem harden_nonneg (β : R) (hβ : (0 : R) ≤ β)
    (theta phi_curl : FieldVec) :
    (0 : R) ≤ energyHarden β theta phi_curl := by
  sorry -- requires ordered field multiplication

-- Theta force vanishes in vacuum
theorem hardenTheta_zero_in_vacuum (β : R) (theta : FieldVec) (a : Fin 3) :
    hardenThetaForce β (fun _ => (0 : R)) theta a = (0 : R) := by
  unfold hardenThetaForce normSq dot
  simp [R.mul_zero, R.zero_add, R.zero_mul]
  sorry -- neg_zero

-- Q vanishes in vacuum
theorem hardenQ_zero_in_vacuum (β : R) (theta : FieldVec) :
    hardenQ β theta (fun _ => (0 : R)) = fun _ => (0 : R) := by
  funext a
  unfold hardenQ
  simp [R.mul_zero]
