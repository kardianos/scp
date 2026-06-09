/-
  v59/gaps/gravity/G8G9_Gravity.lean

  Formalization of the STRUCTURAL identities behind the two gravity gaps:

    G8  gravity magnitude:   G_e = (21/16) · α^21,  21 = dim Spin(7),
                              16 = dim Cl(3,1),  5 = 21 − 16 (additive gauge),
                              and the "product over 21 generators" exponent law
                              f_g = α^(21/2),  G_e = (α^(1/2))^(2·21) = α^21.

    G9  scalar vs tensor:    the helicity content of candidate carriers.
                              A scalar carries h=0, a 4-vector {±1}, an
                              antisymmetric 2-form {±1,0}, a symmetric traceless
                              transverse 2-tensor {±2}.  Only the last is spin-2.
                              The propagating v59 gravity mode (scalar carrier)
                              therefore has NO h=±2 ⇒ fails LIGO.

  STATUS: written, NOT built this run (the shared furey_construction/lean
  project must not be `lake build`-ed concurrently). This file is self-contained
  (`import Mathlib`) so it can be checked independently. Open physical content
  is isolated behind `sorry`-free *encodings* where possible, and explicit
  `sorry` (flagged) where the step is a genuine open conjecture.

  Reading tags mirror UNIFIED_THEORY.md:  [thm] machine-checkable identity;
  [conj] structural ansatz, not derived;  [emp] data match (recorded as a real
  inequality with an explicit tolerance, NOT proved from a Lagrangian).
-/

import Mathlib

namespace SCPv59.Gravity

/-! ## Part 0 — the structural integers -/

/-- `dim Spin(7) = (7 choose 2) = 21`. -/
def dimSpin7 : ℕ := Nat.choose 7 2
/-- `dim Cl(3,1) = 2^4 = 16`. -/
def dimCl31 : ℕ := 2 ^ 4
/-- `dim G₂ = 14`. -/
def dimG2 : ℕ := 14
/-- `dim Spin(8) = (8 choose 2) = 28`. -/
def dimSpin8 : ℕ := Nat.choose 8 2

theorem dimSpin7_eq : dimSpin7 = 21 := by decide
theorem dimCl31_eq : dimCl31 = 16 := by decide
theorem dimSpin8_eq : dimSpin8 = 28 := by decide

/-! ## Part 1 — G8 structural identities -/

/-- **[thm] The dual-Coxeter / additive gauge index 5 = 21 − 16.**
    SU(2)_L coupling g_W² = 5√α is ADDITIVE over a 5-generator subset. -/
theorem gauge_index_additive : (5 : ℕ) = dimSpin7 - dimCl31 := by
  rw [dimSpin7_eq, dimCl31_eq]

/-- **[thm] The gravity prefactor 21/16 = dim Spin(7) / dim Cl(3,1).** -/
theorem gravity_prefactor : (dimSpin7 : ℚ) / (dimCl31 : ℚ) = 21 / 16 := by
  rw [dimSpin7_eq, dimCl31_eq]; norm_num

/-- **[thm] 21 is over-determined** (cautionary): the integer 21 equals four
    distinct v59 combinations, so "21 = dim Spin(7)" is not the *unique*
    structural reading of the exponent. -/
theorem twentyone_overdetermined :
    dimSpin7 = 21
    ∧ (28 : ℕ) - 7 = 21
    ∧ (14 : ℕ) + 7 = 21
    ∧ (35 : ℕ) - 14 = 21 := by
  refine ⟨dimSpin7_eq, ?_, ?_, ?_⟩ <;> decide

/-- **[conj→thm-form] The product/multiplicative exponent law.**
    If gravity couples as a coherent PRODUCT over all `dimSpin7` generators,
    each contributing an amplitude `√α`, then the squared coupling is
    `(√α)^(2·dimSpin7) = α^dimSpin7`.  We state this as the algebraic identity
    `(α^(1/2))^(2·21) = α^21` (rpow) for `α ≥ 0`. -/
theorem product_exponent_law (α : ℝ) (hα : 0 ≤ α) :
    (α ^ ((1:ℝ)/2)) ^ (2 * (dimSpin7 : ℝ)) = α ^ (dimSpin7 : ℝ) := by
  -- (α^(1/2))^(2n) = α^((1/2)·(2n)) = α^n, using Real.rpow_mul for α ≥ 0.
  rw [← Real.rpow_mul hα]
  congr 1
  ring

/-- **[thm] determinant realization of the product.**  The determinant of the
    `n × n` scalar matrix `c • I` is `c^n`.  With `c = √α` and `n = 21` this is
    `(√α)^21 = α^(21/2) = f_g`, the "top-form / volume over 21 generators".
    Here proved as the abstract det identity. -/
theorem det_scalar_matrix_pow (n : ℕ) (c : ℝ) :
    Matrix.det ((c : ℝ) • (1 : Matrix (Fin n) (Fin n) ℝ)) = c ^ n := by
  rw [Matrix.det_smul]
  simp

/-- The grav coupling² closed form: `G_e_form α = (21/16) α^21`. -/
noncomputable def G_e_form (α : ℝ) : ℝ :=
  ((dimSpin7 : ℝ) / (dimCl31 : ℝ)) * α ^ dimSpin7

theorem G_e_form_eq (α : ℝ) : G_e_form α = (21 / 16) * α ^ 21 := by
  unfold G_e_form
  rw [dimSpin7_eq, dimCl31_eq]
  norm_num

/-- **[emp] the magnitude match, stated as a tolerance inequality** (NOT
    derived).  With α = α(0) the IR value and `α_G(e) = (m_e/M_Pl)²`, the v59
    form matches at 0.25%.  We encode this as: there EXISTS an α in the
    physical CODATA band such that the relative error is below 3·10⁻³.
    (We give the witness numerically; this is an empirical statement.) -/
theorem G_e_empirical_match :
    ∃ α aGe : ℝ, 0 < α ∧ 0 < aGe ∧
      |G_e_form α - aGe| / aGe < (3 : ℝ) / 1000 := by
  -- witnesses: α = 1/137.035999084, aGe = (m_e/M_Pl)^2 ≈ 1.751810e-45
  refine ⟨1 / 137.035999084, 1.751810e-45, by norm_num, by norm_num, ?_⟩
  -- The arithmetic |(21/16)α^21 − aGe|/aGe ≈ 0.0025 < 0.003 is a finite decimal
  -- check.  norm_num can in principle discharge it but α^21 is a large power;
  -- we leave it as a flagged empirical step.
  sorry  -- [emp] flagged: 0.25% numeric match, verified in g8_exponent_test.py

/-! ## Part 2 — G9 helicity / polarization content

We model the on-shell little-group content of a massless mode propagating along
`z` by the SO(2)_z action on the relevant index space, and read the helicity
multiset.  Rather than build the full rep theory in Lean, we record the KNOWN
helicity multisets as definitions and prove the decisive *structural* facts:
  (i) the symmetric traceless transverse 2-tensor is the UNIQUE candidate
      containing ±2;
  (ii) a scalar / vector / antisymmetric-2-form do NOT contain ±2;
  (iii) the internal so(8) index does not add spacetime helicity (it is a
      spectator multiset of zeros under the SPACETIME little group). -/

/-- Helicity multisets (as `Multiset ℤ`) of the candidate massless carriers,
    propagating along z, physical (transverse) polarizations only. -/
def hel_scalar  : Multiset ℤ := {0}
def hel_vector  : Multiset ℤ := {-1, 1}
def hel_twoform : Multiset ℤ := {-1, 1}      -- antisymmetric rank-2: spin-1 content
def hel_sym2TT  : Multiset ℤ := {-2, 2}      -- symmetric traceless transverse: spin-2

/-- **[thm] only the symmetric-TT tensor carries h = ±2.** -/
theorem only_sym2_has_spin2 :
    (2 : ℤ) ∈ hel_sym2TT
    ∧ (2 : ℤ) ∉ hel_scalar
    ∧ (2 : ℤ) ∉ hel_vector
    ∧ (2 : ℤ) ∉ hel_twoform := by
  refine ⟨?_, ?_, ?_, ?_⟩ <;> decide

/-- **[thm] the scalar carrier (v59 gravity-as-density) is h = 0 only.**
    `box Ω_grav = f_g ρ_grav` with `ρ_grav = Tr(M†M)` a Lorentz SCALAR ⇒ the
    propagating spacetime carrier is a scalar ⇒ helicity {0}, no ±2. -/
theorem v59_gravity_scalar_no_spin2 :
    hel_scalar = {0} ∧ (2 : ℤ) ∉ hel_scalar ∧ (-2 : ℤ) ∉ hel_scalar := by
  refine ⟨rfl, by decide, by decide⟩

/-- **[thm] the internal index does not add spacetime helicity.**
    Tensoring the scalar carrier (helicity 0) with an `n`-dim INTERNAL (so(8))
    multiplet gives `n` copies of helicity 0 — still no ±2.  We model the
    n-fold internal multiplet as `Multiset.replicate n 0`; every member is 0,
    so ±2 are absent. -/
theorem internal_index_inert (n : ℕ) :
    (2 : ℤ) ∉ Multiset.replicate n 0 ∧ (-2 : ℤ) ∉ Multiset.replicate n 0 := by
  refine ⟨fun h => ?_, fun h => ?_⟩
  · have := Multiset.eq_of_mem_replicate h; norm_num at this
  · have := Multiset.eq_of_mem_replicate h; norm_num at this

/-- **[thm] a spacetime antisymmetric 2-form (bivector connection, EM-like)
    is spin-1, not spin-2.** Its physical helicities are {±1}; in particular no
    ±2.  So even promoting Ω to a SPACETIME bivector F_{μν} does not yield the
    LIGO tensor mode. -/
theorem spacetime_bivector_is_spin1 :
    hel_twoform = {-1, 1} ∧ (2 : ℤ) ∉ hel_twoform := by
  refine ⟨rfl, by decide⟩

/-- **[thm] the quadratic (stress-tensor) route DOES reach spin-2, but is
    quadratic.**  The symmetric bilinear `T_{μν} = Ω_μ^a Ω_ν^a` is symmetric,
    so its TT part is the spin-2 multiset {±2}.  We encode "symmetric ⇒ can be
    TT-projected to {±2}" by exhibiting the target multiset; the PHYSICS caveat
    (quadratic ⇒ O(α), short-range, a source not a fundamental graviton) is a
    Lagrangian statement, recorded as a comment, not provable here. -/
theorem bilinear_reaches_spin2 :
    hel_sym2TT = {-2, 2} ∧ (2 : ℤ) ∈ hel_sym2TT := by
  refine ⟨rfl, by decide⟩

/-! ## Part 3 — consolidated gravity-gap summary -/

/-- **The gravity-gap structural summary.**  Collects the machine-checkable
    parts of G8 (the integers and the product-exponent identity) and G9 (the
    helicity classification).  The genuinely OPEN parts — the empirical
    magnitude match (G8) and the existence of a fundamental symmetric-tensor
    massless mode (G9) — are NOT in this conjunction; they are flagged separately
    (`G_e_empirical_match` has a `sorry`; G9's LIGO failure is the *absence* of
    a spin-2 carrier in the present formulation). -/
theorem gravity_gap_summary :
    -- G8: prefactor and additive index
    ((dimSpin7 : ℚ) / (dimCl31 : ℚ) = 21 / 16)
    ∧ ((5 : ℕ) = dimSpin7 - dimCl31)
    -- G8: 21 over-determined (caution)
    ∧ (dimSpin7 = 21 ∧ (28:ℕ) - 7 = 21)
    -- G9: only sym-TT tensor has ±2; scalar/vector/2-form do not
    ∧ ((2 : ℤ) ∈ hel_sym2TT)
    ∧ ((2 : ℤ) ∉ hel_scalar)
    ∧ ((2 : ℤ) ∉ hel_twoform) := by
  refine ⟨gravity_prefactor, gauge_index_additive, ⟨dimSpin7_eq, by decide⟩, ?_, ?_, ?_⟩
  · decide
  · decide
  · decide

end SCPv59.Gravity
