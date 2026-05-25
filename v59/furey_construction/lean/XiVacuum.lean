/-
  v59/furey_construction/lean/XiVacuum.lean

  Formalisation of task (4b) from the v59 Lean review: the **Goldstone structure of
  the Brannen-ξ vacuum** — the SM-Higgs-like mass spectrum of the dynamical lepton
  field ξ ∈ ℍ ≅ ℝ⁴.

  The v59 dynamical Higgs candidate (see `synthesis/DYNAMIC_XI.md`,
  `FINDINGS_dynamic_xi.md`) has the Mexican-hat potential

      V(ξ) = (λ/4) (|ξ|² − 1/2)²          on ℍ ≅ ℝ⁴,

  whose minima are the constraint surface |ξ|² = 1/2 — exactly the v59 Koide
  constraint `t² = 1/2` (`BrannenKernel.koide_iff_constraint`).  Numerically the
  Hessian at the vacuum was found to be `diag(λ, 0, 0, 0)`: one massive radial mode
  and three massless Goldstones.  This file proves that, rigorously:

    * `V_critical_point`   : the vacuum is a critical point (∂V vanishes there);
    * `V_hessian_directional` : the second directional derivative of V equals
                                `2λ (ξ·v)²` — the Hessian quadratic form;
    * `hessian_eq_massMatrix_quadForm` : that form is `vᵀ (massMatrix) v`, so the
                                explicit `massMatrix = 2λ ξ_a ξ_b` IS the Hessian;
    * `massMatrix_radial`  : the vacuum direction is an eigenvector, eigenvalue λ
                             (the massive radial / "Higgs" mode);
    * `massMatrix_goldstone` : every direction ⊥ to the vacuum is a null eigenvector
                             (massless — the Goldstones);
    * `massMatrix_trace`   : the trace (= Σ eigenvalues) is λ;
    * `xi_mass_spectrum`   : at the canonical vacuum, the three imaginary-quaternion
                             directions are massless Goldstones and the real direction
                             is the massive mode — spectrum `{λ, 0, 0, 0}`.

  Physically: the three Goldstones are the would-be longitudinal W±, Z (eaten when the
  silent SU(2)_L of `SilentDirection.lean` is gauged); the radial λ-mode is the Higgs.
-/
import Mathlib

namespace SCPv59.XiVacuum

open Matrix
open scoped BigOperators

noncomputable section

/-! ## The ξ field, its norm, and the Mexican-hat potential -/

/-- Euclidean inner product on ℍ ≅ ℝ⁴. -/
def dot (p v : Fin 4 → ℝ) : ℝ := ∑ a, p a * v a

/-- Squared norm `|ξ|²`. -/
def normSq (v : Fin 4 → ℝ) : ℝ := ∑ a, (v a) ^ 2

/-- The Brannen-ξ Mexican-hat potential `V(ξ) = (λ/4)(|ξ|² − 1/2)²`. -/
def V (lam : ℝ) (ξ : Fin 4 → ℝ) : ℝ := (lam / 4) * (normSq ξ - 1 / 2) ^ 2

/-- The vacuum manifold: `|ξ|² = 1/2` (the v59 Koide constraint `t² = 1/2`). -/
def IsVacuum (p : Fin 4 → ℝ) : Prop := normSq p = 1 / 2

/-- The straight ray `t ↦ p + t·v` used for directional derivatives. -/
def ray (p v : Fin 4 → ℝ) (t : ℝ) : Fin 4 → ℝ := fun a => p a + t * v a

/-! ## Calculus helpers: derivatives of a quadratic-leading quartic -/

/-- First derivative at 0 of `c₂t² + c₃t³ + c₄t⁴` vanishes. -/
lemma firstDeriv_poly (c2 c3 c4 : ℝ) :
    deriv (fun t : ℝ => c2 * t ^ 2 + (c3 * t ^ 3 + c4 * t ^ 4)) 0 = 0 := by
  apply HasDerivAt.deriv
  have e2 : HasDerivAt (fun t : ℝ => t ^ 2) (2 * (0 : ℝ)) 0 := by simpa using hasDerivAt_pow 2 (0 : ℝ)
  have e3 : HasDerivAt (fun t : ℝ => t ^ 3) (3 * (0 : ℝ) ^ 2) 0 := by simpa using hasDerivAt_pow 3 (0 : ℝ)
  have e4 : HasDerivAt (fun t : ℝ => t ^ 4) (4 * (0 : ℝ) ^ 3) 0 := by simpa using hasDerivAt_pow 4 (0 : ℝ)
  have := (e2.const_mul c2).add ((e3.const_mul c3).add (e4.const_mul c4))
  convert this using 1
  ring

/-- Second derivative at 0 of `c₂t² + c₃t³ + c₄t⁴` is `2c₂` (the Hessian coefficient). -/
lemma secondDeriv_poly (c2 c3 c4 : ℝ) :
    deriv (deriv (fun t : ℝ => c2 * t ^ 2 + (c3 * t ^ 3 + c4 * t ^ 4))) 0 = 2 * c2 := by
  have h1 : deriv (fun t : ℝ => c2 * t ^ 2 + (c3 * t ^ 3 + c4 * t ^ 4))
      = fun x => c2 * (2 * x) + (c3 * (3 * x ^ 2) + c4 * (4 * x ^ 3)) := by
    funext x
    apply HasDerivAt.deriv
    have e2 : HasDerivAt (fun t : ℝ => t ^ 2) (2 * x) x := by simpa using hasDerivAt_pow 2 x
    have e3 : HasDerivAt (fun t : ℝ => t ^ 3) (3 * x ^ 2) x := by simpa using hasDerivAt_pow 3 x
    have e4 : HasDerivAt (fun t : ℝ => t ^ 4) (4 * x ^ 3) x := by simpa using hasDerivAt_pow 4 x
    exact (e2.const_mul c2).add ((e3.const_mul c3).add (e4.const_mul c4))
  rw [h1]
  apply HasDerivAt.deriv
  have a1 : HasDerivAt (fun t : ℝ => 2 * t) 2 (0 : ℝ) := by simpa using (hasDerivAt_id (0 : ℝ)).const_mul 2
  have a2 : HasDerivAt (fun t : ℝ => 3 * t ^ 2) (3 * (2 * (0 : ℝ))) 0 := by
    have h := (show HasDerivAt (fun t : ℝ => t ^ 2) (2 * (0 : ℝ)) 0 by simpa using hasDerivAt_pow 2 (0 : ℝ))
    exact h.const_mul 3
  have a3 : HasDerivAt (fun t : ℝ => 4 * t ^ 3) (4 * (3 * (0 : ℝ) ^ 2)) 0 := by
    have h := (show HasDerivAt (fun t : ℝ => t ^ 3) (3 * (0 : ℝ) ^ 2) 0 by simpa using hasDerivAt_pow 3 (0 : ℝ))
    exact h.const_mul 4
  have := (a1.const_mul c2).add ((a2.const_mul c3).add (a3.const_mul c4))
  convert this using 1
  ring

/-! ## V along a ray through a vacuum point -/

/-- `|p + t v|² = |p|² + 2t (p·v) + t² |v|²`. -/
lemma normSq_ray (p v : Fin 4 → ℝ) (t : ℝ) :
    normSq (ray p v t) = normSq p + 2 * t * dot p v + t ^ 2 * normSq v := by
  simp only [normSq, dot, ray, Fin.sum_univ_four]; ring

/-- At a vacuum point (`|p|² = 1/2`), `V` along the ray is the quartic
    `λ(p·v)² t² + λ(p·v)|v|² t³ + (λ/4)|v|⁴ t⁴` — note the **vanishing constant and
    linear terms** (vacuum is a critical point) and the quadratic coefficient
    `λ(p·v)²` (half the Hessian form). -/
lemma V_along_ray (lam : ℝ) {p : Fin 4 → ℝ} (v : Fin 4 → ℝ) (t : ℝ) (hp : IsVacuum p) :
    V lam (ray p v t)
      = (lam * (dot p v) ^ 2) * t ^ 2
        + ((lam * (dot p v) * normSq v) * t ^ 3 + ((lam / 4) * (normSq v) ^ 2) * t ^ 4) := by
  have hp' : normSq p = 1 / 2 := hp
  unfold V
  rw [normSq_ray, hp']; ring

/-! ## The vacuum is a critical point with Hessian form `2λ (p·v)²` -/

/-- The first directional derivative of `V` at a vacuum point vanishes in every
    direction: the vacuum manifold consists of critical points. -/
theorem V_critical_point (lam : ℝ) {p : Fin 4 → ℝ} (v : Fin 4 → ℝ) (hp : IsVacuum p) :
    deriv (fun t => V lam (ray p v t)) 0 = 0 := by
  have hfun : (fun t => V lam (ray p v t))
      = fun t => (lam * (dot p v) ^ 2) * t ^ 2
          + ((lam * (dot p v) * normSq v) * t ^ 3 + ((lam / 4) * (normSq v) ^ 2) * t ^ 4) := by
    funext t; exact V_along_ray lam v t hp
  rw [hfun]; exact firstDeriv_poly _ _ _

/-- The Hessian quadratic form of `V` at a vacuum point: the second directional
    derivative along `v` is `2λ (p·v)²`.  It is positive along the radial direction
    (`v ∝ p`) and **zero in every direction orthogonal to `p`** — the Goldstones. -/
theorem V_hessian_directional (lam : ℝ) {p : Fin 4 → ℝ} (v : Fin 4 → ℝ) (hp : IsVacuum p) :
    deriv (deriv (fun t => V lam (ray p v t))) 0 = 2 * lam * (dot p v) ^ 2 := by
  have hfun : (fun t => V lam (ray p v t))
      = fun t => (lam * (dot p v) ^ 2) * t ^ 2
          + ((lam * (dot p v) * normSq v) * t ^ 3 + ((lam / 4) * (normSq v) ^ 2) * t ^ 4) := by
    funext t; exact V_along_ray lam v t hp
  rw [hfun, secondDeriv_poly]; ring

/-! ## The mass matrix (Hessian) and its eigenstructure -/

/-- The mass matrix `M_ab = ∂²V/∂ξ_a∂ξ_b` at the vacuum: the rank-one form `2λ p_a p_b`. -/
def massMatrix (lam : ℝ) (p : Fin 4 → ℝ) : Matrix (Fin 4) (Fin 4) ℝ :=
  Matrix.of fun a b => 2 * lam * (p a) * (p b)

/-- The quadratic form of `massMatrix` reproduces the Hessian form `2λ (p·v)²`. -/
theorem massMatrix_quadForm (lam : ℝ) (p v : Fin 4 → ℝ) :
    dotProduct v ((massMatrix lam p).mulVec v) = 2 * lam * (dot p v) ^ 2 := by
  simp only [massMatrix, Matrix.mulVec, dotProduct, Matrix.of_apply, dot, Fin.sum_univ_four]; ring

/-- **`massMatrix` is the Hessian of `V`**: its quadratic form equals the second
    directional derivative of `V` at the vacuum. -/
theorem hessian_eq_massMatrix_quadForm (lam : ℝ) {p : Fin 4 → ℝ} (v : Fin 4 → ℝ) (hp : IsVacuum p) :
    deriv (deriv (fun t => V lam (ray p v t))) 0 = dotProduct v ((massMatrix lam p).mulVec v) := by
  rw [V_hessian_directional lam v hp, massMatrix_quadForm]

/-- **Massive radial mode.** The vacuum direction is an eigenvector of the mass matrix
    with eigenvalue `λ` (the Higgs). -/
theorem massMatrix_radial (lam : ℝ) {p : Fin 4 → ℝ} (hp : IsVacuum p) :
    (massMatrix lam p).mulVec p = lam • p := by
  have hp' : (p 0) ^ 2 + (p 1) ^ 2 + (p 2) ^ 2 + (p 3) ^ 2 = 1 / 2 := by
    have := hp; simp only [IsVacuum, normSq, Fin.sum_univ_four] at this; linarith
  funext a
  simp only [massMatrix, Matrix.mulVec, dotProduct, Matrix.of_apply, Fin.sum_univ_four,
    Pi.smul_apply, smul_eq_mul]
  linear_combination (2 * lam * (p a)) * hp'

/-- **Massless Goldstone modes.** Every direction orthogonal to the vacuum is a null
    eigenvector of the mass matrix (eigenvalue 0). -/
theorem massMatrix_goldstone (lam : ℝ) {p w : Fin 4 → ℝ} (hw : dot p w = 0) :
    (massMatrix lam p).mulVec w = 0 := by
  have hw' : p 0 * w 0 + p 1 * w 1 + p 2 * w 2 + p 3 * w 3 = 0 := by
    simpa [dot, Fin.sum_univ_four] using hw
  funext a
  simp only [massMatrix, Matrix.mulVec, dotProduct, Matrix.of_apply, Fin.sum_univ_four, Pi.zero_apply]
  linear_combination (2 * lam * (p a)) * hw'

/-- The trace (= sum of eigenvalues) of the mass matrix at the vacuum is `λ`,
    consistent with the spectrum `{λ, 0, 0, 0}`. -/
theorem massMatrix_trace (lam : ℝ) {p : Fin 4 → ℝ} (hp : IsVacuum p) :
    (massMatrix lam p).trace = lam := by
  have hp' : (p 0) ^ 2 + (p 1) ^ 2 + (p 2) ^ 2 + (p 3) ^ 2 = 1 / 2 := by
    have := hp; simp only [IsVacuum, normSq, Fin.sum_univ_four] at this; linarith
  simp only [Matrix.trace, Matrix.diag, massMatrix, Matrix.of_apply, Fin.sum_univ_four]
  linear_combination (2 * lam) * hp'

/-! ## The canonical vacuum and the {λ, 0, 0, 0} spectrum -/

/-- Canonical vacuum: `(1/√2, 0, 0, 0)` — the "real" quaternion direction. -/
def vac : Fin 4 → ℝ := fun i => if i = 0 then Real.sqrt 2 / 2 else 0

/-- The three imaginary-quaternion directions `i, j, k` (the would-be SU(2)_L Goldstones). -/
def im1 : Fin 4 → ℝ := fun i => if i = 1 then 1 else 0
def im2 : Fin 4 → ℝ := fun i => if i = 2 then 1 else 0
def im3 : Fin 4 → ℝ := fun i => if i = 3 then 1 else 0

lemma vac_isVacuum : IsVacuum vac := by
  have h2 : Real.sqrt 2 ^ 2 = 2 := Real.sq_sqrt (by norm_num)
  simp only [IsVacuum, normSq, Fin.sum_univ_four]
  rw [show vac 0 = Real.sqrt 2 / 2 from rfl, show vac 1 = (0 : ℝ) from rfl,
      show vac 2 = (0 : ℝ) from rfl, show vac 3 = (0 : ℝ) from rfl]
  rw [div_pow, h2]; norm_num

/-- **The ξ mass spectrum is `{λ, 0, 0, 0}`** (one massive Higgs + three Goldstones).
    The real direction `vac` is the massive radial eigenvector (eigenvalue λ); the
    three independent imaginary directions `im1, im2, im3` are massless Goldstones. -/
theorem xi_mass_spectrum (lam : ℝ) :
    (massMatrix lam vac).mulVec vac = lam • vac
    ∧ (massMatrix lam vac).mulVec im1 = 0
    ∧ (massMatrix lam vac).mulVec im2 = 0
    ∧ (massMatrix lam vac).mulVec im3 = 0 := by
  refine ⟨massMatrix_radial lam vac_isVacuum, ?_, ?_, ?_⟩
  · exact massMatrix_goldstone lam (by simp [dot, vac, im1])
  · exact massMatrix_goldstone lam (by simp [dot, vac, im2])
  · exact massMatrix_goldstone lam (by simp [dot, vac, im3])

end

end SCPv59.XiVacuum
