/-
  v59/furey_construction/lean/HiggsVevReframe.lean

  **Integration target #1 (v_Higgs), step 1: the clean structural form.**

  The conjecture `v_Higgs = 28²·a_l²` (suspicious as a fitted "784") reframes into a relation
  between two small structural integers — the number of generations and dim(L):

      Σ√m_lepton / √v_Higgs  =  N_gen / dim(L)  =  3 / 28      (empirical: 0.03%).

  Equivalently `√v_Higgs = dim(L)·a_l` with `a_l = (Σ√m)/N_gen` the Brannen scale.  So the
  lepton √-mass quantum is `√v_Higgs / dim(L)` — the EW-vacuum scale per mass-bearing
  (`L = Λ²⊕Λ⁶`) direction.  This module proves the algebraic reframing and the integer value;
  the v58-dynamical derivation of *why* the quantum is `√v/dim(L)` is in
  `INTEGRATION_PLAN_v58.md` / the v_Higgs derivation note.
-/
import Mathlib

namespace SCPv59.HiggsVev

/-- `dim(L) = dim(Λ²⊕Λ⁶) = 28`, the mass-bearing (chiral/complex) grade of `Cl(7)_even`
    (proven to carry the complex structure / mass in `BladeSquareSign`, `LeptonRealityForcing`). -/
def dimL : ℕ := 28
/-- Number of generations (the sedenion `S₃`/`Z₃`; `sedenion_s3.py`). -/
def nGen : ℕ := 3

/-! ## The algebraic reframing -/

/-- **`v_Higgs = dim(L)²·a²  ⟺  √v_Higgs = dim(L)·a`** (for `v,a ≥ 0`).  Turns the fitted `28²`
    into the linear mode-count `√v = dim(L)·a`. -/
theorem vHiggs_sqrt_form (v a : ℝ) (hv : 0 ≤ v) (ha : 0 ≤ a) :
    v = ((dimL : ℝ)) ^ 2 * a ^ 2 ↔ Real.sqrt v = (dimL : ℝ) * a := by
  rw [show ((dimL : ℝ)) ^ 2 * a ^ 2 = ((dimL : ℝ) * a) ^ 2 by ring]
  constructor
  · intro h; rw [h, Real.sqrt_sq (by positivity)]
  · intro h; rw [← h, Real.sq_sqrt hv]

/-- **`√v = 28·a  ⟺  (Σ√m)/√v = 3/28`** when `a = (Σ√m)/3`, `v>0`  (here `28 = dim(L)`,
    `3 = N_gen`).  This is the cleanest form: the ratio of total √-lepton-mass to √-EW-scale is
    the integer ratio `N_gen/dim(L) = 3/28`. -/
theorem ratio_form (v S a : ℝ) (hv : 0 < v) (hgen : a = S / 3) :
    Real.sqrt v = 28 * a ↔ S / Real.sqrt v = 3 / 28 := by
  have hsv : Real.sqrt v ≠ 0 := ne_of_gt (Real.sqrt_pos.mpr hv)
  subst hgen
  rw [div_eq_div_iff hsv (by norm_num : (28 : ℝ) ≠ 0)]
  constructor
  · intro h; linarith [h]
  · intro h; linarith [h]

/-- The structural value of the ratio: `N_gen / dim(L) = 3/28`. -/
theorem ratio_value : (nGen : ℚ) / (dimL : ℚ) = 3 / 28 := by norm_num [nGen, dimL]

/-! ## (iii) backed into known physics: the lepton Yukawa sum rule

With the standard Higgs–Yukawa relation `m_f = y_f·(v/√2)`, the clean form `Σ√m/√v = 3/28`
becomes a **Yukawa sum rule** `Σ√y_lepton = (3/28)·2^{1/4}` — i.e. the lepton Yukawa is
*geometrically diluted* by the ambient dimension `dim(L)=28`.  This is the known-physics content
of piece (iii) of the v58 derivation (`integration_v58/01_higgs_vev.md`). -/

/-- `m_f = y_f·(v/√2)  ⇒  √m_f = √y_f · √v / 2^{1/4}`  (with `2^{1/4} = √√2`). -/
theorem sqrt_mass_from_yukawa (y v : ℝ) (hy : 0 ≤ y) (hv : 0 ≤ v) :
    Real.sqrt (y * (v / Real.sqrt 2)) = Real.sqrt y * (Real.sqrt v / Real.sqrt (Real.sqrt 2)) := by
  rw [Real.sqrt_mul hy, Real.sqrt_div hv]

/-- **The Yukawa-sum form.**  `Σ √m_f = (√v / 2^{1/4})·Σ √y_f`, so `Σ√m/√v = (Σ√y)/2^{1/4}`;
    hence the clean form `Σ√m/√v = 3/28` is the lepton Yukawa sum rule `Σ√y = (3/28)·2^{1/4}`. -/
theorem yukawa_sum_form (v y₁ y₂ y₃ : ℝ) (hv : 0 ≤ v) (h₁ : 0 ≤ y₁) (h₂ : 0 ≤ y₂) (h₃ : 0 ≤ y₃) :
    Real.sqrt (y₁ * (v / Real.sqrt 2)) + Real.sqrt (y₂ * (v / Real.sqrt 2))
        + Real.sqrt (y₃ * (v / Real.sqrt 2))
      = (Real.sqrt v / Real.sqrt (Real.sqrt 2)) * (Real.sqrt y₁ + Real.sqrt y₂ + Real.sqrt y₃) := by
  rw [sqrt_mass_from_yukawa y₁ v h₁ hv, sqrt_mass_from_yukawa y₂ v h₂ hv,
      sqrt_mass_from_yukawa y₃ v h₃ hv]; ring

/-! ## Empirical match (PDG) -/

/-- `Σ√m_lepton / √v_Higgs ≈ 0.1071061` (PDG masses, `v=246.22 GeV`) matches `3/28 = 0.1071428`
    to `~3×10⁻⁴`. -/
theorem ratio_agrees : |(3 / 28 : ℚ) - 1071062 / 10000000| < 1 / 2000 := by norm_num

/-! ## Vector vs. bilinear: WHY the factor is `dim(L)`, not `√dim(L)`

The Cl(7) bridge computation (`integration_v58/bridge_vhiggs_cl7.py`,
`03_higgs_bridge_result.md`) shows the 28 `L`-blades are Frobenius-**orthogonal**.  That
distinguishes two readings of "the vacuum lives on the `dim(L)`-dim grade with per-element
scale `a`", and only one reproduces the empirical `√v/a = 28`:

* **vector-in-`L` (the EXCLUDED "equipartition" reading):** a democratic vector over `dim(L)`
  *orthonormal* directions has squared energy norm `dim(L)·a²`, so `√v = √dim(L)·a` — for
  `dim(L)=28` this is `√28 ≈ 5.29`, not 28.
* **bilinear-on-`L` (the CORRECT reading):** a Yukawa is a *matrix* on `L`, with
  `dim(L)² = 784` components; its Frobenius² is `dim(L)²·a²`, so `√v = dim(L)·a = 28a`. -/

/-- **Bilinear (Frobenius²) reading.**  A democratic `dim(L)×dim(L)` matrix (every entry `a`)
    has Frobenius² `= Σᵢⱼ a² = dim(L)²·a²` — the `784 = dim(L)²` is the **number of components**
    of the mass bilinear on `L`. -/
theorem frobeniusSq_democratic (n : ℕ) (a : ℝ) :
    ∑ _i : Fin n, ∑ _j : Fin n, a ^ 2 = (n : ℝ) ^ 2 * a ^ 2 := by
  simp only [Finset.sum_const, Finset.card_univ, Fintype.card_fin, nsmul_eq_mul]
  ring

/-- **Bilinear reading gives `√v = dim(L)·a`** (the correct, natural-`L²` answer): the Frobenius²
    `dim(L)²·a²` has root `dim(L)·a`. -/
theorem bilinear_reading (a : ℝ) (ha : 0 ≤ a) :
    Real.sqrt (((dimL : ℝ)) ^ 2 * a ^ 2) = (dimL : ℝ) * a := by
  rw [show ((dimL : ℝ)) ^ 2 * a ^ 2 = ((dimL : ℝ) * a) ^ 2 by ring]
  exact Real.sqrt_sq (by positivity)

/-- **Vector reading gives `√v = √dim(L)·a`** (the excluded "equipartition" answer): a democratic
    vector over `dim(L)` orthonormal directions has squared norm `dim(L)·a²`. -/
theorem vector_reading (a : ℝ) (ha : 0 ≤ a) :
    Real.sqrt (((dimL : ℝ)) * a ^ 2) = Real.sqrt (dimL : ℝ) * a := by
  rw [Real.sqrt_mul (by positivity), Real.sqrt_sq ha]

/-- **The two readings genuinely differ**: `√dim(L) ≠ dim(L)` (`√28 ≈ 5.29 ≠ 28`), so the
    equipartition reading is excluded by the empirical `√v/a = 28` — only the bilinear
    (component-count) reading survives. -/
theorem readings_differ : Real.sqrt (dimL : ℝ) ≠ (dimL : ℝ) := by
  have hd : (dimL : ℝ) = 28 := by norm_num [dimL]
  rw [hd]; intro h
  have hsq : (Real.sqrt 28) ^ 2 = (28 : ℝ) ^ 2 := by rw [h]
  rw [Real.sq_sqrt (by norm_num : (0 : ℝ) ≤ 28)] at hsq
  norm_num at hsq

end SCPv59.HiggsVev
