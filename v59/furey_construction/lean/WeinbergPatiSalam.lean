/-
  v59/furey_construction/lean/WeinbergPatiSalam.lean

  **Tier 2 of the φ = 2/9 program**: derive the gauge-sector `2/9` (the weak mixing angle)
  from the Pati-Salam coupling matching, rather than defining it.

  `ScaleBridge.sin_sq_thW := 2/9` is a *definition*.  Here we DERIVE it from the structural
  coupling-squared coefficients (in units of `√α`, which cancels in the ratio):

      g_W²     = cW·√α,   cW = 5  = h∨(Spin(7))   (dual Coxeter; `GaugePrefactorDualCoxeter`)
      g_{B-L}² = cBL·√α,  cBL = 2                  ((T2.2) the "2" — still structurally open)

  Pati-Salam matching of the hypercharge `U(1)_Y ⊂ SU(2)_R × U(1)_{B-L}` (with left-right
  symmetry `g_R = g_L = g_W`):

      1/g'² = 1/g_W² + 1/g_{B-L}²     ⇒     g'² = (10/7)·√α

  and then

      sin²θ_W = g'² / (g_W² + g'²) = (10/7)/(5 + 10/7) = 10/45 = **2/9**.

  The `√α` cancels, so this is exact rational arithmetic once `(cW, cBL) = (5, 2)`.
  The `5` is the dual Coxeter number (derived); the `2` is the remaining open input (T2.2).
-/
import Mathlib
import GaugePrefactorDualCoxeter

namespace SCPv59.WeinbergPatiSalam

/-! ## Coupling-squared coefficients (in units of `√α`) -/

/-- `g_W² = cW·√α` with `cW = 5 = h∨(Spin(7))` (dual Coxeter; `GaugePrefactorDualCoxeter`). -/
def cW : ℚ := 5
/-- `g_{B-L}² = cBL·√α` with `cBL = 2`.  **(T2.2) OPEN**: the structural origin of this `2`
    (B-L abelian factor / silent-direction count / Brannen numerator) is not yet derived. -/
def cBL : ℚ := 2

/-- `cW = 5` is the dual Coxeter number `h∨(Spin(7)) = 7 − 2`. -/
theorem cW_eq_dualCoxeter : cW = (SCPv59.GaugePrefactor.dualCoxeterSO 7 : ℚ) := by
  rw [SCPv59.GaugePrefactor.dualCoxeter_so7]; rfl

/-! ## Pati-Salam matching -/

/-- The Pati-Salam matching condition for the hypercharge coupling (coefficient form, `√α`
    cancelled): `1/c' = 1/cW + 1/cBL`. -/
def PatiSalam (cW cBL cP : ℚ) : Prop := 1 / cP = 1 / cW + 1 / cBL

/-- The hypercharge coefficient solving the matching: `c' = cW·cBL/(cW+cBL) = 10/7`. -/
def cPrime : ℚ := cW * cBL / (cW + cBL)

theorem cPrime_value : cPrime = 10 / 7 := by unfold cPrime cW cBL; norm_num

/-- `cPrime` genuinely satisfies the Pati-Salam matching. -/
theorem cPrime_satisfies_patiSalam : PatiSalam cW cBL cPrime := by
  unfold PatiSalam cPrime cW cBL; norm_num

/-! ## The weak mixing angle -/

/-- `sin²θ_W = g'²/(g_W² + g'²) = c'/(cW + c')` (the `√α` cancels). -/
def sin2_thetaW : ℚ := cPrime / (cW + cPrime)

/-- **Tier 2 result.**  `sin²θ_W = 2/9`, derived from `(cW, cBL) = (5, 2)` and the Pati-Salam
    matching — not assumed. -/
theorem sin2_thetaW_eq : sin2_thetaW = 2 / 9 := by
  unfold sin2_thetaW cPrime cW cBL; norm_num

/-- **General form.**  For any positive coefficients satisfying the matching, the mixing angle
    is `c'/(cW+c')`; specialised to `(5,2)` it is `2/9`.  This is the honest chain: the physics
    input is `(cW=5, cBL=2)` + Pati-Salam; the `2/9` is the consequence. -/
theorem sin2_from_coeffs (cP : ℚ) (hW : cW = 5) (hBL : cBL = 2)
    (hPS : PatiSalam cW cBL cP) : cP / (cW + cP) = 2 / 9 := by
  unfold PatiSalam at hPS
  simp only [hW, hBL] at hPS ⊢
  -- 1/cP = 1/5 + 1/2 = 7/10 ⇒ cP = 10/7
  have hcP : cP = 10 / 7 := by
    have hne : cP ≠ 0 := by
      rintro rfl; norm_num at hPS
    field_simp at hPS; linarith
  rw [hcP]; norm_num

/-! ## Connection to the existing (defined) value -/

/-- The Pati-Salam derivation reproduces the value `ScaleBridge.sin_sq_thW := 2/9` defines —
    i.e. the EW-sector `2/9` is now *derived* from `(5, 2)`, not posited. -/
theorem matches_scaleBridge : sin2_thetaW = (2 / 9 : ℚ) := sin2_thetaW_eq

end SCPv59.WeinbergPatiSalam
