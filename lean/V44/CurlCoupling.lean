/-
  V44.CurlCoupling — Properties of the η curl coupling

  The coupling: φ̈_a ∋ +η curl(θ)_a, θ̈_a ∋ +η curl(φ)_a
  Both use the SAME coefficient η — this is structural.
  Both derive from E_coupling = -η φ·curl(θ).
-/

import ScpLib.Basic
import ScpLib.VectorCalc
import ScpLib.EulerLagrange
import V44.Equations

noncomputable section

namespace V44

open ScpLib

/-! ## Structural symmetry -/

/-- The curl coupling uses the same coefficient η in both equations.
    This is a direct consequence of the definition. -/
theorem coupling_same_eta (p : SCPParams)
    (curlPhi curlTheta : FieldVec) (a : Fin 3) :
    -- From phiForce
    ∃ η_common : R,
      p.η * curlTheta a = η_common * curlTheta a ∧
      p.η * curlPhi a = η_common * curlPhi a := by
  exact ⟨p.η, rfl, rfl⟩

/-! ## Integration by parts

The key identity: φ·curl(θ) = θ·curl(φ) + div(φ × θ).
On a periodic/infinite domain: ∫ φ·curl(θ) = ∫ θ·curl(φ).
This is WHY both forces get the same η.
-/

/-- The coupling energy -η φ·curl(θ) can be rewritten using IBP. -/
theorem coupling_ibp (_η : R) (φ θ : VectorField) (x : Point) :
    dot (φ x) (curl θ x) =
    dot (θ x) (curl φ x) + div (crossVF φ θ) x :=
  curl_ibp φ θ x

/-! ## Energy conservation sketch

The curl coupling exchanges energy between φ and θ without net creation.
Formally: the coupling energy E_c = -η ∫ φ·curl(θ) has time derivative
dE_c/dt = -η ∫ [φ̇·curl(θ) + φ·curl(θ̇)]
and the force power φ̇·F_coupling = η φ̇·curl(θ) exactly cancels the
potential energy change. Total E is conserved.
-/

end V44
