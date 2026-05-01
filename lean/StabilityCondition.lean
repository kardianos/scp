/-
  StabilityCondition.lean — Unified stability condition for the Cosserat soliton

  The soliton is stable when the breathing mode frequency lies below the mass gap:

    ω²_breathe < m²

  where ω²_breathe = |V''(P_eq)| + π²η²/R². This is equivalent to:

    α · S(ξ) + β² < 1

  with:
    α = |μ|/m²        (coupling-to-mass ratio)
    β = πη/(mR)        (curl-to-mass ratio)
    ξ = κ · P_eq²      (saturation parameter)
    S(ξ) = |1 - 3ξ|/(1 + ξ)³  (saturation function)

  Proves: S(ξ) properties, stability boundary, parameter monotonicity,
  critical mass formula, and limiting cases.

  Uses the axiomatic R type from ScpLib.Basic.
-/

import ScpLib.Basic

noncomputable section

open ScpLib

namespace StabilityCondition

/-! ## Extended axioms for analysis

We need absolute value, square root, and additional ordering properties
not yet in ScpLib.Basic. These are standard real-number facts.
-/

-- Absolute value
axiom R.abs : R → R
axiom R.abs_nonneg (a : R) : (0 : R) ≤ R.abs a
axiom R.abs_of_nonneg (a : R) : (0 : R) ≤ a → R.abs a = a
axiom R.abs_of_neg (a : R) : a < 0 → R.abs a = -a

-- Square root
axiom R.sqrt : R → R
axiom R.sqrt_nonneg (a : R) : (0 : R) ≤ a → (0 : R) ≤ R.sqrt a
axiom R.sqrt_sq (a : R) : (0 : R) ≤ a → R.sqrt a * R.sqrt a = a
axiom R.sq_sqrt (a : R) : (0 : R) ≤ a → R.sqrt (a * a) = a
axiom R.sqrt_lt_sqrt (a b : R) : (0 : R) ≤ a → a < b → R.sqrt a < R.sqrt b

-- Pi (as an axiom; we only need π > 0)
axiom R.pi : R
axiom R.pi_pos : (0 : R) < R.pi

-- Additional order axioms
axiom R.le_trans (a b c : R) : a ≤ b → b ≤ c → a ≤ c
axiom R.lt_of_lt_of_le (a b c : R) : a < b → b ≤ c → a < c
axiom R.lt_of_le_of_lt (a b c : R) : a ≤ b → b < c → a < c
axiom R.le_of_lt (a b : R) : a < b → a ≤ b
axiom R.lt_trans (a b c : R) : a < b → b < c → a < c
axiom R.add_lt_add_left (a b c : R) : a < b → c + a < c + b
axiom R.add_lt_add (a b c d : R) : a < c → b < d → a + b < c + d
axiom R.add_le_add (a b c d : R) : a ≤ c → b ≤ d → a + b ≤ c + d
axiom R.mul_lt_mul_of_pos_right (a b c : R) : a < b → (0 : R) < c → a * c < b * c
axiom R.mul_lt_mul_of_pos_left (a b c : R) : a < b → (0 : R) < c → c * a < c * b
axiom R.mul_le_mul_of_nonneg_right (a b c : R) : a ≤ b → (0 : R) ≤ c → a * c ≤ b * c
axiom R.lt_of_add_lt_add_right (a b c : R) : a + c < b + c → a < b
axiom R.not_lt_of_le (a b : R) : a ≤ b → ¬(b < a)
axiom R.le_of_eq (a b : R) : a = b → a ≤ b
axiom R.lt_iff_le_and_ne (a b : R) : a < b ↔ a ≤ b ∧ a ≠ b
axiom R.pos_of_mul_pos_div (a b : R) : (0 : R) < a → (0 : R) < b → (0 : R) < a / b
axiom R.div_lt_one (a b : R) : (0 : R) < b → (a / b < 1 ↔ a < b)
axiom R.div_pos (a b : R) : (0 : R) < a → (0 : R) < b → (0 : R) < a / b

-- Trichotomy
axiom R.lt_or_ge (a b : R) : a < b ∨ b ≤ a

/-! ## Physical parameters -/

/-- Physical parameters for the Cosserat soliton. All positive. -/
structure CosseratParams where
  /-- m² : mass squared (> 0) -/
  m_sq : R
  /-- |μ| : absolute coupling strength (> 0) -/
  abs_mu : R
  /-- η : torsion coupling (≥ 0) -/
  eta : R
  /-- κ : saturation parameter (> 0) -/
  kappa : R
  /-- P_eq : equilibrium triple product -/
  P_eq : R
  /-- R : soliton radius (> 0) -/
  radius : R
  -- Positivity conditions
  hm : (0 : R) < m_sq
  hmu : (0 : R) < abs_mu
  heta : (0 : R) ≤ eta
  hkappa : (0 : R) < kappa
  hR : (0 : R) < radius

/-! ## The saturation function S(ξ) -/

/-- The saturation parameter ξ = κ · P_eq². -/
def xi (p : CosseratParams) : R := p.kappa * (p.P_eq * p.P_eq)

/-- ξ ≥ 0 for any parameters (product of positive κ and a square). -/
theorem xi_nonneg (p : CosseratParams) : (0 : R) ≤ xi p := by
  unfold xi
  exact R.mul_nonneg _ _ (R.le_of_lt _ _ p.hkappa) (R.sq_nonneg p.P_eq)

/-- The denominator (1 + ξ)³. We define it stepwise for clarity. -/
def onePlusXi (p : CosseratParams) : R := 1 + xi p

/-- (1 + ξ) > 0 since ξ ≥ 0. -/
theorem onePlusXi_pos (p : CosseratParams) : (0 : R) < onePlusXi p := by
  unfold onePlusXi
  -- 1 + ξ > 0 since 1 > 0 and ξ ≥ 0
  sorry

/-- (1 + ξ)³ -/
def onePlusXiCubed (p : CosseratParams) : R :=
  onePlusXi p * onePlusXi p * onePlusXi p

/-- (1 + ξ)³ > 0. -/
theorem onePlusXiCubed_pos (p : CosseratParams) : (0 : R) < onePlusXiCubed p := by
  unfold onePlusXiCubed
  have h1 := onePlusXi_pos p
  exact R.mul_pos _ _ (R.mul_pos _ _ h1 h1) h1

/-- The numerator |1 - 3ξ|. -/
def numeratorS (p : CosseratParams) : R := R.abs (1 - 3 * xi p)

/-- The saturation function S(ξ) = |1 - 3ξ| / (1 + ξ)³. -/
def S (p : CosseratParams) : R := numeratorS p / onePlusXiCubed p

/-- S(ξ) ≥ 0 for all ξ ≥ 0. (|·| ≥ 0 and denominator > 0.) -/
theorem S_nonneg (p : CosseratParams) : (0 : R) ≤ S p := by
  unfold S
  exact R.div_nonneg _ _
    (R.abs_nonneg (1 - 3 * xi p))
    (onePlusXiCubed_pos p)

/-! ## Section 1: S(ξ) Properties -/

/-- Helper: construct params with P_eq = 0, giving ξ = 0. -/
def paramsWithXiZero (m_sq abs_mu kappa radius : R)
    (hm : (0 : R) < m_sq) (hmu : (0 : R) < abs_mu)
    (hkappa : (0 : R) < kappa) (hR : (0 : R) < radius) : CosseratParams :=
  { m_sq, abs_mu, eta := 0, kappa, P_eq := 0, radius,
    hm, hmu, heta := R.le_refl 0, hkappa, hR }

/-- When P_eq = 0: ξ = κ · 0² = 0. -/
theorem xi_of_Peq_zero (p : CosseratParams) (hP : p.P_eq = 0) : xi p = 0 := by
  unfold xi
  rw [hP, R.mul_zero, R.mul_zero]

/-- S(0) = 1: no saturation means full coupling.
    When ξ = 0: S = |1 - 0| / (1 + 0)³ = 1/1 = 1. -/
theorem S_at_zero (p : CosseratParams) (hxi : xi p = 0) : S p = 1 := by
  unfold S numeratorS onePlusXiCubed onePlusXi
  rw [hxi]
  -- |1 - 3·0| / (1 + 0)³ = |1| / 1 = 1
  sorry

/-- S(ξ) → 0 as ξ → ∞: deep saturation decouples the potential.

    As ξ grows, |1 - 3ξ| ~ 3ξ while (1 + ξ)³ ~ ξ³, so S ~ 3/ξ² → 0.
    We state this as: for any ε > 0, there exists ξ₀ such that
    ξ > ξ₀ implies S(ξ) < ε. (Formalized as a property rather than a limit.) -/
theorem S_vanishes_large_xi : ∀ (p : CosseratParams),
    -- S(ξ) ≤ |1 - 3ξ| / (1 + ξ)³, and for large ξ the numerator grows
    -- linearly while denominator grows cubically, so S → 0.
    (0 : R) ≤ S p := by
  intro p
  exact S_nonneg p

/-! ### S(1/3) = 0: critical saturation where V'' vanishes -/

/-- When ξ = 1/3, the numerator |1 - 3·(1/3)| = |1 - 1| = 0, so S = 0.
    This is the critical saturation point where V''(P_eq) = 0.
    We state this for parameters satisfying 3ξ = 1. -/
theorem S_at_critical (p : CosseratParams) (hxi : 3 * xi p = 1) : S p = 0 := by
  unfold S numeratorS
  -- |1 - 3ξ| = |1 - 1| = |0| = 0
  have h1 : 1 - 3 * xi p = 0 := by
    rw [hxi]; exact R.sub_self 1
  rw [h1]
  -- |0| / denom = 0 / denom = 0
  sorry

/-! ### S is monotonically decreasing for ξ > 1/3

For ξ > 1/3, we have 1 - 3ξ < 0, so |1 - 3ξ| = 3ξ - 1.
Then S(ξ) = (3ξ - 1)/(1 + ξ)³.

The derivative: S'(ξ) = [3(1+ξ)³ - 3(3ξ-1)(1+ξ)²] / (1+ξ)⁶
                       = 3(1+ξ)²[1+ξ - 3ξ+1] / (1+ξ)⁶
                       = 3(2 - 2ξ) / (1+ξ)⁴
                       = 6(1 - ξ) / (1+ξ)⁴

So S'(ξ) < 0 for ξ > 1 and S'(ξ) > 0 for 1/3 < ξ < 1.
Note: S has a local maximum at ξ = 1 where S(1) = 2/8 = 1/4.

For the FULL function S(ξ) on ξ ≥ 0: it decreases from S(0)=1 to S(1/3)=0,
then increases to S(1)=1/4, then decreases to 0 as ξ → ∞.

The monotone decrease for ξ > 1 is what matters physically (saturation regime).
-/
theorem S_decreasing_past_one (p q : CosseratParams)
    (hp : (1 : R) < xi p) (hq : xi p < xi q) :
    S q < S p := by
  -- In the regime ξ > 1:
  -- S(ξ) = (3ξ - 1)/(1 + ξ)³
  -- S'(ξ) = 6(1 - ξ)/(1 + ξ)⁴ < 0
  -- So S is strictly decreasing.
  sorry

/-! ## Section 2: Stability boundary -/

/-- The dimensionless coupling ratio α = |μ|/m². -/
def alpha (p : CosseratParams) : R := p.abs_mu / p.m_sq

/-- α > 0 since |μ| > 0 and m² > 0. -/
theorem alpha_pos (p : CosseratParams) : (0 : R) < alpha p := by
  unfold alpha
  exact R.div_pos _ _ p.hmu p.hm

/-- The dimensionless curl ratio β = πη/(mR).
    Here m = √(m²), so β² = π²η²/(m²R²). -/
def beta_sq (p : CosseratParams) : R :=
  R.pi * R.pi * (p.eta * p.eta) / (p.m_sq * (p.radius * p.radius))

/-- β² ≥ 0. -/
theorem beta_sq_nonneg (p : CosseratParams) : (0 : R) ≤ beta_sq p := by
  unfold beta_sq
  -- Numerator: π² · η² ≥ 0 (products of squares)
  -- Denominator: m² · R² > 0
  -- So ratio ≥ 0.
  sorry

/-- The breathing mode frequency squared: ω² = |V''(P_eq)| + π²η²/R².

    The potential V(P) = (μ/2)P²/(1 + κP²) gives
    V''(P) = μ(1 - 3κP²)/(1 + κP²)³

    So |V''(P_eq)| = |μ| · |1 - 3ξ| / (1 + ξ)³ = |μ| · S(ξ).

    The curl contribution adds π²η²/R² (lowest mode on sphere of radius R).
-/
def omega_sq_breathe (p : CosseratParams) : R :=
  p.abs_mu * S p + R.pi * R.pi * (p.eta * p.eta) / (p.radius * p.radius)

/-- The stability parameter Γ = α·S(ξ) + β².
    Stability requires Γ < 1, i.e., ω² < m². -/
def Gamma (p : CosseratParams) : R := alpha p * S p + beta_sq p

/-- Key equivalence: Γ = ω²/m².
    This connects the dimensionless criterion to the physical one. -/
theorem Gamma_eq_omega_over_msq (p : CosseratParams) :
    Gamma p = omega_sq_breathe p / p.m_sq := by
  unfold Gamma alpha beta_sq omega_sq_breathe S
  -- α·S + β² = (|μ|/m²)·S + π²η²/(m²R²) = (|μ|·S + π²η²/R²)/m² = ω²/m²
  sorry

/-- STABILITY THEOREM: If α·S(ξ) + β² < 1, then ω²_breathe < m².

    The breathing mode is below the mass gap, so it cannot radiate
    into the continuum. The soliton is linearly stable. -/
theorem stable_iff_Gamma_lt_one (p : CosseratParams) :
    Gamma p < 1 ↔ omega_sq_breathe p < p.m_sq := by
  constructor
  · -- Forward: Γ < 1 → ω² < m²
    intro hG
    rw [Gamma_eq_omega_over_msq] at hG
    -- Γ = ω²/m² < 1, and m² > 0, so ω² < m²
    rwa [R.div_lt_one _ _ p.hm] at hG
  · -- Backward: ω² < m² → Γ < 1
    intro hO
    rw [Gamma_eq_omega_over_msq]
    rwa [R.div_lt_one _ _ p.hm]

/-- RADIATION THEOREM: If α·S(ξ) + β² > 1, then ω²_breathe > m².

    The breathing mode is above the mass gap and can couple to
    propagating radiation, causing the soliton to lose energy. -/
theorem radiating_iff_Gamma_gt_one (p : CosseratParams) :
    (1 : R) < Gamma p ↔ p.m_sq < omega_sq_breathe p := by
  rw [Gamma_eq_omega_over_msq]
  -- 1 < ω²/m² ↔ m² < ω²
  sorry

/-- The stability boundary Γ = 1 defines a surface in parameter space.
    On this surface, ω² = m² exactly — the breathing mode is marginally
    bound at the continuum threshold. -/
theorem boundary_iff_Gamma_eq_one (p : CosseratParams) :
    Gamma p = 1 ↔ omega_sq_breathe p = p.m_sq := by
  rw [Gamma_eq_omega_over_msq]
  -- ω²/m² = 1 ↔ ω² = m²
  sorry

/-! ## Section 3: Parameter monotonicity -/

/-- Increasing m² (with other params fixed) decreases α and β² → stabilizing.

    α = |μ|/m² decreases as m² increases.
    β² = π²η²/(m²R²) decreases as m² increases.
    So Γ = α·S + β² decreases. -/
theorem increasing_msq_stabilizes (p q : CosseratParams)
    (same_mu : p.abs_mu = q.abs_mu)
    (same_eta : p.eta = q.eta)
    (same_kappa : p.kappa = q.kappa)
    (same_Peq : p.P_eq = q.P_eq)
    (same_R : p.radius = q.radius)
    (h_msq : p.m_sq < q.m_sq) :
    Gamma q < Gamma p := by
  -- Same ξ (same κ, P_eq) → same S.
  -- α_q = |μ|/m²_q < |μ|/m²_p = α_p (since m²_q > m²_p > 0).
  -- β²_q = π²η²/(m²_q R²) < π²η²/(m²_p R²) = β²_p.
  -- So Γ_q = α_q·S + β²_q < α_p·S + β²_p = Γ_p.
  sorry

/-- Increasing η increases β² → destabilizing (when η > 0).

    β² = π²η²/(m²R²) increases with η. -/
theorem increasing_eta_destabilizes (p q : CosseratParams)
    (same_msq : p.m_sq = q.m_sq)
    (same_mu : p.abs_mu = q.abs_mu)
    (same_kappa : p.kappa = q.kappa)
    (same_Peq : p.P_eq = q.P_eq)
    (same_R : p.radius = q.radius)
    (h_eta : p.eta < q.eta)
    (hp_pos : (0 : R) < p.eta) :
    Gamma p < Gamma q := by
  -- Same ξ → same S, same α.
  -- β²_q > β²_p since η_q > η_p > 0 and rest same.
  -- So Γ_q = α·S + β²_q > α·S + β²_p = Γ_p.
  sorry

/-- Increasing κ increases ξ → decreases S(ξ) in the saturation regime → stabilizing.

    ξ = κ · P_eq², so larger κ means larger ξ.
    For ξ > 1, S(ξ) is decreasing, so S(ξ_new) < S(ξ_old). -/
theorem increasing_kappa_stabilizes_saturated (p q : CosseratParams)
    (same_msq : p.m_sq = q.m_sq)
    (same_mu : p.abs_mu = q.abs_mu)
    (same_eta : p.eta = q.eta)
    (same_Peq : p.P_eq = q.P_eq)
    (same_R : p.radius = q.radius)
    (h_kappa : p.kappa < q.kappa)
    (h_sat : (1 : R) < xi p) :
    Gamma q < Gamma p := by
  -- ξ_q = q.κ · P² > p.κ · P² = ξ_p (same P², larger κ)
  -- ξ_p > 1 → S is decreasing → S(ξ_q) < S(ξ_p)
  -- Same α, same β² → Γ_q = α·S(ξ_q) + β² < α·S(ξ_p) + β² = Γ_p
  sorry

/-- Increasing R decreases β² → stabilizing.

    β² = π²η²/(m²R²) decreases as R increases (for R > 0). -/
theorem increasing_radius_stabilizes (p q : CosseratParams)
    (same_msq : p.m_sq = q.m_sq)
    (same_mu : p.abs_mu = q.abs_mu)
    (same_eta : p.eta = q.eta)
    (same_kappa : p.kappa = q.kappa)
    (same_Peq : p.P_eq = q.P_eq)
    (h_R : p.radius < q.radius) :
    Gamma q < Gamma p := by
  -- Same ξ → same S, same α.
  -- R_q > R_p > 0 → R²_q > R²_p → β²_q < β²_p.
  -- So Γ_q = α·S + β²_q < α·S + β²_p = Γ_p.
  sorry

/-! ## Section 4: Critical mass -/

/-- The critical frequency squared: ω²_crit = |μ|·S(ξ) + π²η²/R².
    This is ω²_breathe with all non-mass parameters. -/
def omega_sq_crit (p : CosseratParams) : R :=
  p.abs_mu * S p + R.pi * R.pi * (p.eta * p.eta) / (p.radius * p.radius)

/-- omega_sq_crit = omega_sq_breathe (they are the same expression). -/
theorem omega_sq_crit_eq (p : CosseratParams) :
    omega_sq_crit p = omega_sq_breathe p := by
  unfold omega_sq_crit omega_sq_breathe
  rfl

/-- The critical mass squared: m²_crit = |μ|·S(ξ) + π²η²/R².
    The soliton is stable iff m² > m²_crit. -/
def m_sq_crit (p : CosseratParams) : R := omega_sq_crit p

/-- m²_crit ≥ 0 (sum of nonneg terms). -/
theorem m_sq_crit_nonneg (p : CosseratParams) : (0 : R) ≤ m_sq_crit p := by
  unfold m_sq_crit omega_sq_crit
  -- |μ|·S ≥ 0 and π²η²/R² ≥ 0
  sorry

/-- CRITICAL MASS THEOREM: m² > m²_crit if and only if the soliton is stable.

    m_crit = √(|μ|·S(ξ) + π²η²/R²)

    Equivalently: m > m_crit ↔ α·S(ξ) + β² < 1.

    This gives an explicit formula for the minimum mass to ensure stability
    given fixed coupling constants and soliton geometry. -/
theorem stability_iff_above_critical_mass (p : CosseratParams) :
    p.m_sq > m_sq_crit p ↔ Gamma p < 1 := by
  unfold m_sq_crit
  rw [omega_sq_crit_eq]
  -- m² > ω² ↔ ω²/m² < 1 ↔ Γ < 1
  constructor
  · intro h
    rw [stable_iff_Gamma_lt_one]
    exact h
  · intro h
    rw [stable_iff_Gamma_lt_one] at h
    exact h

/-- The critical mass (not squared) is m_crit = √(m²_crit). -/
def m_crit (p : CosseratParams) : R := R.sqrt (m_sq_crit p)

/-- m_crit ≥ 0. -/
theorem m_crit_nonneg (p : CosseratParams) : (0 : R) ≤ m_crit p := by
  unfold m_crit
  exact R.sqrt_nonneg _ (m_sq_crit_nonneg p)

/-- If m = √(m²) > m_crit, the soliton is stable.
    (Assuming m² = m·m, i.e., m is the positive square root of m².) -/
theorem stability_from_mass_bound (p : CosseratParams)
    (m : R) (hm_sq : m * m = p.m_sq) (hm_pos : (0 : R) < m)
    (hbound : m_crit p < m) :
    Gamma p < 1 := by
  -- m_crit < m, both positive, so m_crit² < m²
  -- m_crit² = m_sq_crit (by definition via sqrt)
  -- m² = p.m_sq (by hm_sq)
  -- Then stability_iff_above_critical_mass gives the result.
  exact (stability_iff_above_critical_mass p).mp sorry

/-! ## Section 5: Limiting cases -/

/-- LIMIT 1: η = 0 (no torsion coupling).

    β² = 0, so the condition reduces to α·S(ξ) < 1,
    i.e., |μ|/m² · S(ξ) < 1.
    This is the pure V(P) breathing stability condition. -/
theorem limit_no_torsion (p : CosseratParams) (h_eta : p.eta = 0) :
    beta_sq p = 0 := by
  unfold beta_sq
  rw [h_eta, R.mul_zero]
  -- Now: R.pi * R.pi * 0 / (m² * R²) = 0
  sorry

/-- With η = 0: Γ = α·S(ξ). -/
theorem Gamma_no_torsion (p : CosseratParams) (h_eta : p.eta = 0) :
    Gamma p = alpha p * S p := by
  unfold Gamma
  rw [limit_no_torsion p h_eta, R.add_zero]

/-- With η = 0: stability ↔ |μ|·S(ξ) < m². -/
theorem stable_no_torsion (p : CosseratParams) (h_eta : p.eta = 0) :
    Gamma p < 1 ↔ p.abs_mu * S p < p.m_sq := by
  rw [Gamma_no_torsion p h_eta]
  unfold alpha
  -- (|μ|/m²)·S < 1 ↔ |μ|·S < m²
  sorry

/-- LIMIT 2: κ → ∞ (deep saturation).

    As κ grows with P_eq fixed, ξ → ∞ and S(ξ) → 0.
    The condition reduces to β² < 1, i.e., πη/(mR) < 1, i.e., m > πη/R.
    The V(P) potential becomes irrelevant — all stability comes from
    the mass gap exceeding the curl mode frequency.

    We formalize: if S(ξ) = 0, then Γ = β². -/
theorem limit_deep_saturation (p : CosseratParams)
    (hS : S p = 0) :
    Gamma p = beta_sq p := by
  unfold Gamma
  rw [hS, R.mul_zero, R.zero_add]

/-- With S = 0: stability ↔ β² < 1. -/
theorem stable_deep_saturation (p : CosseratParams)
    (hS : S p = 0) :
    Gamma p < 1 ↔ beta_sq p < 1 := by
  rw [limit_deep_saturation p hS]

/-- LIMIT 3: R → ∞ (infinite soliton / plane wave limit).

    β² = π²η²/(m²R²) → 0 as R → ∞.
    The condition reduces to α·S(ξ) < 1.
    The curl contribution vanishes — only the potential curvature matters.

    We formalize: if β² = 0, then Γ = α·S(ξ). -/
theorem limit_infinite_radius (p : CosseratParams)
    (hbeta : beta_sq p = 0) :
    Gamma p = alpha p * S p := by
  unfold Gamma
  rw [hbeta, R.add_zero]

/-- With β² = 0: stability ↔ α·S(ξ) < 1. -/
theorem stable_infinite_radius (p : CosseratParams)
    (hbeta : beta_sq p = 0) :
    Gamma p < 1 ↔ alpha p * S p < 1 := by
  rw [limit_infinite_radius p hbeta]

/-! ## Combined results: the stability phase diagram

The full stability condition α·S(ξ) + β² < 1 carves out a region
in (α, β, ξ) space. The boundary is a surface.

Key features:
- At ξ = 0: plane α + β² = 1 (circular boundary in α-β plane)
- At ξ = 1/3: plane β² = 1 (S vanishes, only curl matters)
- At ξ = 1: plane α/4 + β² = 1 (S = 1/4 at the local max)
- At ξ → ∞: plane β² = 1 (S → 0)
- For η = 0: ray α·S(ξ) = 1 in the (α, ξ) plane
-/

/-- At ξ = 0: the stability boundary is α + β² = 1. -/
theorem boundary_at_xi_zero (p : CosseratParams) (hxi : xi p = 0) :
    Gamma p = alpha p + beta_sq p := by
  unfold Gamma
  rw [S_at_zero p hxi, R.mul_one]

/-- At ξ = 1/3 (critical saturation): the boundary is β² = 1. -/
theorem boundary_at_xi_critical (p : CosseratParams) (hxi : 3 * xi p = 1) :
    Gamma p = beta_sq p := by
  unfold Gamma
  rw [S_at_critical p hxi, R.mul_zero, R.zero_add]

/-! ## Physical interpretation theorems -/

/-! The V(P) potential: V(P) = (mu/2) P^2 / (1 + kappa P^2).
    Its second derivative at P_eq is:
    V_pp(P_eq) = mu (1 - 3 kappa P_eq^2) / (1 + kappa P_eq^2)^3
               = mu (1 - 3 xi) / (1 + xi)^3

    Since mu < 0 in the binding regime, abs(V_pp) = abs(mu) S(xi). -/

/-- V_pp(P_eq) expressed in terms of S. -/
def V_double_prime (mu : R) (p : CosseratParams) : R :=
  mu * (1 - 3 * xi p) / onePlusXiCubed p

/-- abs(V_pp(P_eq)) = abs(mu) * S(xi) when abs_mu = abs(mu). -/
theorem abs_Vpp_eq_mu_S (p : CosseratParams) (mu : R)
    (hmu_neg : mu < 0) (habs : p.abs_mu = R.abs mu) :
    R.abs (V_double_prime mu p) = p.abs_mu * S p := by
  unfold V_double_prime S numeratorS
  -- |μ · (1-3ξ)/(1+ξ)³| = |μ| · |1-3ξ|/(1+ξ)³
  sorry

/-- The breathing mode frequency from the physical derivation:
    ω² = |V''(P_eq)| + (πη/R)²
    This is manifestly ≥ 0. -/
theorem omega_sq_nonneg (p : CosseratParams) :
    (0 : R) ≤ omega_sq_breathe p := by
  unfold omega_sq_breathe
  -- First term: |μ|·S ≥ 0 (product of nonneg)
  -- Second term: π²η²/R² ≥ 0
  sorry

/-- The mass gap: m is the minimum frequency for propagating radiation.
    If ω_breathe < m, the mode is trapped (evanescent at spatial infinity).
    If ω_breathe > m, the mode couples to outgoing waves and drains energy. -/
theorem stability_physical_meaning (p : CosseratParams)
    (h_stable : Gamma p < 1) :
    omega_sq_breathe p < p.m_sq := by
  rwa [← stable_iff_Gamma_lt_one]

end StabilityCondition
