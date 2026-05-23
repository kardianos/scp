/-
  v59/furey_construction/lean/ScaleBridge.lean

  EW-sector closure: v_Higgs scale bridge, sin²θ_W, cos²θ_W, m_Z/m_W,
  and α(M_Z) — all from v59 structural inputs plus the empirical Brannen
  lepton scale a_l.

  Conjectures encoded (see `synthesis/FINDINGS_scale_bridge.md`):

    (S1) v_Higgs = D_lepton² · a_l²              where D_lepton = 28
    (S2) sin²θ_W = 2/9 = Brannen phase           (= Q_lepton / 3)
    (S3) cos²θ_W = 7/9 = t²_{u-quark}            (1 - dim G_2 / D_u-quark)
    (S4) m_Z / m_W = √(1 / cos²θ_W) = 3/√7        (custodial)
    (S5) α(M_Z) = 25 / (324 π²) = (5/(18π))²
         from g_W² = 5√α  ∧  sin²θ_W = 2/9  ∧  4πα = g_W²·sin²θ_W

  The arithmetic identities are proved; the empirical match (numerical
  agreement with PDG) is not — that is documented in the findings.
-/

import Mathlib.Tactic
import Predictions
import SpinDimension

namespace SCPv59.ScaleBridge

open SCPv59 SCPv59.SpinDimension SCPv59.Predictions

/-! ## Scale bridge: v_Higgs = D_lepton² · a_l² -/

/-- D_lepton = 28 = dim Spin(8) = the lepton-sector ambient dimension
    (Λ²ℝ⁷ ⊕ Λ⁶ℝ⁷, the "L" content of Cl(7)_even).
    Also the adjoint dimension of the algebra whose Z₃-triality generates
    the three generations. -/
def dimLepton : ℕ := 28

/-- Sanity: dim_lepton equals the L content of Cl(7)_even. -/
theorem dimLepton_eq_L : dimLepton = L_content := by
  unfold dimLepton; rw [L_dim]

/-- Sanity: dim_lepton equals (7 choose 2) + (7 choose 6) = dim Λ² + dim Λ⁶. -/
theorem dimLepton_eq_decomp :
    dimLepton = Nat.choose 7 2 + Nat.choose 7 6 := by decide

/-- v_Higgs (predicted from v59 structural inputs and the lepton Brannen
    scale a_l): `v = D_lepton² · a_l²`. -/
noncomputable def v_Higgs (a_l_sq : ℝ) : ℝ := (dimLepton : ℝ)^2 * a_l_sq

theorem v_Higgs_factor : (dimLepton : ℝ)^2 = 784 := by
  show ((28 : ℕ) : ℝ)^2 = 784
  norm_num

/-! ## Weak mixing angle: sin²θ_W = Brannen phase = 2/9 -/

/-- v59 conjecture: `sin²θ_W = 2/9`. -/
def sin_sq_thW : ℚ := 2/9

/-- sin²θ_W matches the Brannen phase from Tier 2 of the prediction table:
    `sin²θ_W = ((dim G_2 / dim Spin(7)) / 3) = (Q_lepton / 3) = 2/9`. -/
theorem sin_sq_thW_eq_brannen_phase :
    sin_sq_thW = ((dimG2 : ℚ) / (dimSO 7 : ℚ)) / 3 := by
  unfold sin_sq_thW
  rw [brannen_phase_structural]

/-- v59 conjecture: `cos²θ_W = 7/9`. -/
def cos_sq_thW : ℚ := 7/9

/-- sin² + cos² = 1, structurally. -/
theorem sin_cos_sq_thW_sum : sin_sq_thW + cos_sq_thW = 1 := by
  unfold sin_sq_thW cos_sq_thW; norm_num

/-- cos²θ_W matches the u-quark Brannen t²:
    `cos²θ_W = t²_{u-quark} = 1 − dim G_2 / dim U63 = 7/9`.

    This is the structural connection: the EW weak mixing angle is
    set by the u-quark Brannen parameter. -/
theorem cos_sq_thW_eq_t_sq_u_quark : cos_sq_thW = t_sq_u_quark := by
  unfold cos_sq_thW
  rw [t_sq_u_quark_eq]

/-! ## m_Z/m_W from custodial symmetry: 3/√7 -/

/-- Custodial-symmetry / on-shell relation: `m_Z² = m_W² / cos²θ_W`.
    With cos²θ_W = 7/9, this gives `m_Z² / m_W² = 9/7`, i.e. `m_Z/m_W = 3/√7`. -/
theorem mZ_over_mW_sq : (1 : ℚ) / cos_sq_thW = 9/7 := by
  unfold cos_sq_thW; norm_num

/-! ## α(M_Z) from EW-tree consistency

The SM tree-level relation `4πα = g²·g'²/(g² + g'²) = g²·sin²θ_W` combined
with the v59 conjectures `g_W² = 5√α` and `sin²θ_W = 2/9` forces a
specific value of `α` at the EW scale:

    4π·α = 5·√α · 2/9
    →  √α = 5·(2/9)/(4π) = 5/(18π)
    →  α  = 25 / (324 π²)
-/

/-- v59-derived α at the EW scale: `α(M_Z) = 25 / (324 π²)`.

    This is a CONSEQUENCE (not a free parameter) of two v59 conjectures
    plus the SM tree identity. -/
noncomputable def alpha_MZ : ℝ := 25 / (324 * Real.pi^2)

/-- √α(M_Z) = 5/(18π) — the cleanest structural form.
    Numerator: 5 = (dim Spin(7) − dim Cl(3,1)) (Killing-form embedding index).
    Denominator: 18·π = 2·9·π where 9 = denom of Brannen phase 2/9, 2 from
    rearranging 4πα = (g²)·(sin²θ_W). -/
theorem sqrt_alpha_MZ_form : Real.sqrt alpha_MZ = 5 / (18 * Real.pi) := by
  unfold alpha_MZ
  have hπ : (0 : ℝ) ≤ 5 / (18 * Real.pi) := by
    apply div_nonneg <;> [norm_num; positivity]
  rw [show (25 : ℝ) / (324 * Real.pi^2) = (5 / (18 * Real.pi))^2 by
    have : (18 : ℝ) * Real.pi ≠ 0 := by positivity
    field_simp
    ring]
  exact Real.sqrt_sq hπ

/-- The "factored" form of √α(M_Z): Killing-index × Brannen-phase / (4π). -/
theorem sqrt_alpha_MZ_factored :
    Real.sqrt alpha_MZ
      = ((dimSO 7 : ℝ) - (dimCl31 : ℝ)) * ((2 : ℝ) / 9) / (4 * Real.pi) := by
  rw [sqrt_alpha_MZ_form]
  rw [show ((dimSO 7 : ℝ) - (dimCl31 : ℝ)) = 5 by
    rw [dimSpin_seven]; norm_num [dimCl31]]
  have : (4 : ℝ) * Real.pi ≠ 0 := by positivity
  field_simp
  ring

/-- The α(M_Z) value (squared form). -/
theorem alpha_MZ_eq : alpha_MZ = 25 / (324 * Real.pi^2) := rfl

/-- EW-tree consistency: the conjecture pair (`g_W² = 5√α`, `sin²θ_W = 2/9`)
    together with `4πα = g_W² · sin²θ_W` forces `√α = 5/(18π)`.

    We state the consequence and prove the value follows. -/
theorem alpha_MZ_from_consistency
    (α : ℝ) (hα_pos : α > 0)
    (h_consistency : 5 * Real.sqrt α * (2/9 : ℝ) = 4 * Real.pi * α) :
    Real.sqrt α = 5 / (18 * Real.pi) := by
  -- 5·√α·(2/9) = 4πα  ⇒  10√α/9 = 4πα
  -- Using α = √α·√α:  10√α/9 = 4π·√α·√α  ⇒  (since √α > 0) 10/9 = 4π·√α
  -- ⇒ √α = 10/(36π) = 5/(18π).
  have hsa_pos : Real.sqrt α > 0 := Real.sqrt_pos.mpr hα_pos
  have hsa_ne : Real.sqrt α ≠ 0 := ne_of_gt hsa_pos
  have h_α_sq : α = Real.sqrt α * Real.sqrt α := (Real.mul_self_sqrt (le_of_lt hα_pos)).symm
  -- key intermediate: cancel √α on both sides of h_consistency
  -- h_consistency: 5·√α·(2/9) = 4π·α = 4π·(√α·√α)
  -- divide by √α: 5·(2/9) = 4π·√α
  have h1 : 5 * Real.sqrt α * (2/9 : ℝ) = 4 * Real.pi * (Real.sqrt α * Real.sqrt α) := by
    rw [← h_α_sq]; exact h_consistency
  have h2 : (5 * (2/9 : ℝ)) * Real.sqrt α = (4 * Real.pi * Real.sqrt α) * Real.sqrt α := by
    linarith [h1]
  have h3 : (5 * (2/9 : ℝ)) = 4 * Real.pi * Real.sqrt α :=
    mul_right_cancel₀ hsa_ne h2
  -- now solve: √α = (5·2/9)/(4π) = 10/(36π) = 5/(18π)
  have hπ_pos : (0 : ℝ) < 4 * Real.pi := by positivity
  have hπ_ne : (4 : ℝ) * Real.pi ≠ 0 := ne_of_gt hπ_pos
  have h_sqrt_val : Real.sqrt α = (5 * (2/9 : ℝ)) / (4 * Real.pi) := by
    field_simp
    linarith [h3]
  rw [h_sqrt_val]
  field_simp
  ring

/-! ## Pati-Salam Spin(7) decomposition (Step 5)

The Spin(7) symmetry of v59 decomposes (at the algebra level) as
  Spin(7) = G_2 × SU(2)_L × SU(2)_R × U(1)_{B-L}
with dim count 21 = 14 + 3 + 3 + 1.

Each factor plays a specific role:
  G_2          : octonion automorphisms, defines fermion content
  SU(2)_L      : silent direction (v59), gives W± / Z masses
  SU(2)_R      : right-handed SU(2), broken at high scale
  U(1)_{B-L}   : lepton-vs-quark distinguisher (matches μ-bisection of Step 4)

U(1)_Y emerges as the diagonal of SU(2)_R × U(1)_{B-L}:
  Y = 2·T_3^R + (B-L)

With Killing-form couplings (g_W² = 5√α, g_R² = 5√α, g_{B-L}² = 2√α),
the Pati-Salam reduction gives g'² = (10/7)·√α and sin²θ_W = 2/9 EXACTLY,
matching the v59 Brannen-phase identification of sin²θ_W. -/

/-- Spin(7) decomposition into G_2 × SU(2)_L × SU(2)_R × U(1)_{B-L}.
    Dimensions: 14 + 3 + 3 + 1 = 21 = dim Spin(7).  Axiom-free. -/
theorem spin7_pati_salam_decomp :
    dimG2 + 3 + 3 + 1 = dimSO 7 := by
  rw [dimSpin_seven, show dimG2 = 14 from dimG2_eq_14]

/-- Same identity stated more symmetrically as 14 + 7 = 21,
    with "7" = dim SU(2)_L + dim SU(2)_R + dim U(1)_{B-L}. -/
theorem spin7_g2_plus_seven :
    dimG2 + 7 = dimSO 7 := by
  rw [dimSpin_seven, show dimG2 = 14 from dimG2_eq_14]

/-- The "7" in spin7_g2_plus_seven matches dim Im𝕆 = top grade Λ⁶ of Cl(7)_even.
    Same number that appears in cos²θ_W = 7/9, sin²θ_C = 7·α, m_Z/m_W = 3/√7. -/
theorem seven_unifies_su2L_su2R_u1BL_eq_dimImO :
    (3 + 3 + 1 : ℕ) = cl7_grade_lambda6 := by decide

/-! ## Cabibbo angle from v59 structural integers

A new (2026-05-22 fermion-path session) finding: the Cabibbo angle of
the CKM matrix satisfies

    sin²θ_Cabibbo  ≈  (dim Im 𝕆) · α(0)  =  7 · α(0)

to 0.45% match (empirical sin θ_C = 0.22500, predicted = √(7·α(0)) = 0.22601).

This adds the Cabibbo angle to the v59 prediction tier, with 7 = dim Im𝕆
joining the structural-integer roster (= dim S⁷ = dim Λ⁶ℝ⁷ = top grade of
the v59 lepton ambient L). -/

/-- The "7" in the Cabibbo formula is dim Im𝕆 = dim S⁷ = dim Λ⁶ℝ⁷,
    the top grade of Cl(7)_even. -/
def dimImO : ℕ := 7

theorem dimImO_eq_choose_seven : dimImO = Nat.choose 7 6 := by decide

theorem dimImO_eq_lambda6 : dimImO = cl7_grade_lambda6 := by
  unfold dimImO; rfl

/-- Cabibbo squared angle as a function of α: sin²θ_C = 7·α. -/
noncomputable def sin_sq_cabibbo (α : ℝ) : ℝ := (dimImO : ℝ) * α

theorem sin_sq_cabibbo_value (α : ℝ) :
    sin_sq_cabibbo α = 7 * α := by
  unfold sin_sq_cabibbo dimImO
  norm_num

/-- The "7" in sin²θ_C is the SAME 7 as in cos²θ_W = 7/9.  Both come from
    the dim Im𝕆 = dim Λ⁶ℝ⁷ piece of Cl(7)_even. -/
theorem cabibbo_seven_eq_cos_thW_seven :
    (dimImO : ℚ) / 9 = cos_sq_thW := by
  unfold cos_sq_thW dimImO
  norm_num

/-! ## Universal Koide-deviation identity

A new (this-session) observation: across all three v59 fermion sectors,
the product `(1 - Q_N) · D_N` is INVARIANT and equals `D_lepton / 3 = 28/3`.

Lepton:   (1 - 2/3)  · 28 = (1/3)  · 28 = 28/3
d-quark:  (1 - 11/15)· 35 = (4/15) · 35 = 28/3
u-quark:  (1 - 23/27)· 63 = (4/27) · 63 = 28/3

Equivalently, `Q_N = 1 - (D_lepton/3)/D_N` for all N.

This is a CLEANER restatement of the Brannen pattern `t²_N = 1 − dim G₂/D_N`,
since Q = (1 + 2t²)/3 implies (1 - Q) = 2(1 - t²)/3 = 2·dim G_2 / (3·D_N) =
28/(3·D_N) when 2·dim G_2 = dim Spin(8) = D_lepton (a remarkable identity:
2·14 = 28).
-/

/-- For each v59 sector, (1 - Q_N) · D_N = 28/3.
    This is a uniform invariant across all three fermion sectors. -/
theorem koide_deviation_universal :
    (1 - koide_Q_from_t_sq t_sq_lepton) * 28 = (28 : ℚ)/3
    ∧ (1 - koide_Q_from_t_sq t_sq_d_quark) * dim3FormsR7 = (28 : ℚ)/3
    ∧ (1 - koide_Q_from_t_sq t_sq_u_quark) * dimU63 = (28 : ℚ)/3 := by
  refine ⟨?_, ?_, ?_⟩
  · rw [koide_Q_lepton]; norm_num
  · rw [koide_Q_d_quark]; unfold dim3FormsR7; norm_num
  · rw [koide_Q_u_quark]; unfold dimU63; norm_num

/-- The underlying numerical identity 2·dim G_2 = dim Spin(8) = D_lepton. -/
theorem two_dimG2_eq_dimSpin8 : 2 * (14 : ℕ) = 28 := by decide

/-- Underlying observation: 2·dim G_2 = dim_lepton (where dim_lepton = 28
    is the v59 lepton-sector ambient).  Equivalent to `dim G_2 = 14`. -/
theorem two_dimG2_eq_dimLepton : 2 * dimG2 = dimLepton := by
  unfold dimLepton
  rw [show dimG2 = 14 from dimG2_eq_14]

/-! ## Summary theorem -/

/-- The five new conjectures of the 2026-05-22 scale-bridge work,
    bundled into a single statement of their structural identities. -/
theorem scale_bridge_summary :
    -- (S1) D_lepton = 28 = L_content of Cl(7)_even
    dimLepton = L_content
    -- (S2) sin²θ_W = Brannen phase = 2/9
    ∧ sin_sq_thW = ((dimG2 : ℚ) / (dimSO 7 : ℚ)) / 3
    -- (S3) cos²θ_W = t²_{u-quark} = 7/9
    ∧ cos_sq_thW = t_sq_u_quark
    -- (S2+S3) sin² + cos² = 1
    ∧ sin_sq_thW + cos_sq_thW = 1
    -- (S4) m_Z²/m_W² = 1/cos²θ_W = 9/7
    ∧ (1 : ℚ) / cos_sq_thW = 9/7 := by
  refine ⟨?_, ?_, ?_, ?_, ?_⟩
  · exact dimLepton_eq_L
  · exact sin_sq_thW_eq_brannen_phase
  · exact cos_sq_thW_eq_t_sq_u_quark
  · exact sin_cos_sq_thW_sum
  · exact mZ_over_mW_sq

end SCPv59.ScaleBridge
