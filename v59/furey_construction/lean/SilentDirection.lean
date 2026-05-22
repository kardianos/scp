/-
  v59/furey_construction/lean/SilentDirection.lean

  The Silent-Direction Theorem.

  Empirical finding (this session, `cosserat_experiment/03_scalar_dependence.py`):
  the Brannen mass-operator spectrum depends on ξ ∈ ℍ only via (Re ξ, |Im ξ|).
  Numerical verification: 4 × 10⁻¹⁵ deviation over 1000 random SO(3) rotations
  of Im ξ at fixed (Re ξ, |Im ξ|).

  Structural content:  for any unit quaternion q ∈ ℍ (i.e. |q|² = 1) and any
  ξ ∈ ℍ, the conjugation ξ ↦ q · ξ · q̄ preserves both
      Re(ξ)        (the scalar part)
      |Im(ξ)|²     (the squared imaginary magnitude, ≡ normSq − Re²)
  Hence the Brannen spectral pair (Re ξ, |Im ξ|) is invariant, and the spectrum
  of the cyclic mass operator
      M(ξ) = a · I + ξ · S + ξ̄ · Sᵀ      on ℍ³
  is unchanged.  The group of such rotations is the SU(2)/U(1) coset that
  Furey's ℂ ⊗ ℍ ⊗ 𝕆 picture identifies with the weak SU(2)_L gauge group.
  This file proves the two algebraic invariants.  The matrix-similarity
  consequence (eigenvalues unchanged) is recorded as a comment, as it is
  standard linear algebra (Diag(L_q) conjugates M(ξ) to M(qξq̄)).
-/

import Mathlib.Algebra.Quaternion
import Mathlib.Tactic

namespace SCPv59.SilentDirection

open Quaternion

variable (q ξ : ℍ[ℝ])

/-! ## Step 1 — Re is preserved up to a `normSq q` factor

For any q, ξ ∈ ℍ, the conjugation `q · ξ · star q` has real part `(normSq q) · (Re ξ)`.
This is a polynomial identity in the 8 real components and holds without any
unit constraint on `q`. -/

theorem re_conj : (q * ξ * star q).re = q.normSq * ξ.re := by
  rw [Quaternion.normSq_def']
  simp only [Quaternion.re_mul, Quaternion.imI_mul, Quaternion.imJ_mul,
             Quaternion.imK_mul,
             Quaternion.re_star, Quaternion.imI_star, Quaternion.imJ_star,
             Quaternion.imK_star]
  ring

/-! ## Step 2 — normSq is multiplicative

`normSq : ℍ →*₀ ℝ` is a `MonoidWithZeroHom`, so `normSq (a * b) = normSq a * normSq b`.
With `normSq (star q) = normSq q` this gives `normSq (q · ξ · star q) = (normSq q)² · normSq ξ`. -/

theorem normSq_conj : (q * ξ * star q).normSq = q.normSq^2 * ξ.normSq := by
  rw [map_mul, map_mul, normSq_star]
  ring

/-! ## Step 3 — Specialised to unit q -/

/-- **Silent direction, Re invariant.**  For unit `q` (|q|² = 1),
    conjugation `q · ξ · star q` preserves the real part of ξ. -/
theorem re_conj_unit (hq : q.normSq = 1) : (q * ξ * star q).re = ξ.re := by
  rw [re_conj, hq, one_mul]

/-- **Silent direction, normSq invariant.**  For unit `q`, conjugation preserves
    the full squared norm of ξ.  Combined with `re_conj_unit`, this gives
    |Im(q·ξ·star q)|² = normSq − Re² = |Im ξ|². -/
theorem normSq_conj_unit (hq : q.normSq = 1) : (q * ξ * star q).normSq = ξ.normSq := by
  rw [normSq_conj, hq, one_pow, one_mul]

/-- **The silent-direction invariant pair.**  For unit `q`, the conjugation
    `q · ξ · star q` preserves the *pair* `(Re ξ, normSq ξ)`, which is exactly the
    data on which the Brannen spectrum depends.

    Numerical evidence (`03_scalar_dependence.py`):  1000 random unit-q
    conjugations of a generic ξ leave each Brannen eigenvalue invariant to
    4 × 10⁻¹⁵ (machine precision).  This theorem is the algebraic explanation. -/
theorem silent_pair (hq : q.normSq = 1) :
    (q * ξ * star q).re = ξ.re ∧ (q * ξ * star q).normSq = ξ.normSq :=
  ⟨re_conj_unit q ξ hq, normSq_conj_unit q ξ hq⟩

/-- **Imaginary magnitude is preserved** (corollary of `silent_pair`).
    `|Im(q · ξ · star q)|² = |Im ξ|²` for unit q. -/
theorem im_normSq_conj_unit (hq : q.normSq = 1) :
    (q * ξ * star q).normSq - ((q * ξ * star q).re)^2 = ξ.normSq - ξ.re^2 := by
  rw [re_conj_unit q ξ hq, normSq_conj_unit q ξ hq]

/-! ## Step 4 — Geometric statement of the silent-direction theorem -/

/-- **Silent-direction theorem (geometric form).**  The conjugation action
    `ξ ↦ q · ξ · star q` of unit quaternions q on ℍ has the following property:
    for any ξ, the image lies in the 2-sphere
        S(r, β) = { η ∈ ℍ  :  η.re = r,  η.normSq − η.re² = β² }
    determined by ξ.  Hence the unit-quaternion conjugation acts within the
    2-sphere of imaginary directions at fixed (Re, |Im|), which is the silent
    submanifold of the Brannen kernel.

    In Furey's ℂ ⊗ ℍ ⊗ 𝕆 picture, the unit-quaternion group `{q : normSq q = 1}` ≅ SU(2)
    acts by conjugation on the imaginary directions exactly as SO(3) acts on ℝ³.
    The stabiliser of any imaginary direction is U(1), so the silent
    submanifold is SU(2)/U(1) ≅ S². -/
theorem silent_orbit
    (q ξ : ℍ[ℝ]) (hq : q.normSq = 1) :
    let η := q * ξ * star q
    η.re = ξ.re  ∧  η.normSq - η.re^2 = ξ.normSq - ξ.re^2 :=
  ⟨re_conj_unit q ξ hq, im_normSq_conj_unit q ξ hq⟩

/-! ## Step 5 — Consequence for the Brannen spectrum (documentation)

The cyclic Brannen mass operator on ℍ³ is
    M(ξ) = a · I + L_ξ · S + L_{star ξ} · Sᵀ
where `L_ξ : ℍ → ℍ` is left-multiplication.  Conjugation by the block matrix
    D_q = Diag(L_q, L_q, L_q)
yields
    D_q · M(ξ) · D_q⁻¹ = a · I + L_q L_ξ L_{q⁻¹} S + L_q L_{star ξ} L_{q⁻¹} Sᵀ
                       = a · I + L_{q ξ star q} · S + L_{q (star ξ) star q} · Sᵀ
                       = M(q · ξ · star q)
(using D_q ↔ S commutation, since S permutes flavors but does not act on each
ℍ-factor; and `star (q · ξ · star q) = q · star ξ · star q`).

For unit q, D_q is real-orthogonal, so M(ξ) and M(q · ξ · star q) are
orthogonally similar — they have the *same eigenvalues*.

Combined with the silent-pair theorem above, this gives the empirical result of
`03_scalar_dependence.py`: any two ξ with the same (Re, normSq) — equivalently
the same (Re, |Im|) — give the same Brannen spectrum.

The matrix-similarity step is standard linear algebra and is not re-derived
in Lean here; the algebraic invariance proved above is what makes that
similarity possible.
-/

end SCPv59.SilentDirection
