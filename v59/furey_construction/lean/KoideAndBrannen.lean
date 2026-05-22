/-
  v59/furey_construction/lean/KoideAndBrannen.lean

  Statements of the empirical Koide and Brannen identities and their
  structural derivations from the Lie group dimensions.
-/

import «LieDimensions»

namespace SCPv59.KoideAndBrannen

open SCPv59

/-- Empirical Koide ratio for charged leptons.
    Q_emp = (m_e + m_mu + m_tau) / (sqrt(m_e) + sqrt(m_mu) + sqrt(m_tau))^2
          ≈ 0.6666605
    with the 2/3 value being satisfied to ~5e-6 (within m_tau precision).

    We do NOT prove the empirical claim here; we ENCODE it as a rational
    approximation (the structural value 2/3) and assert it matches data
    by external citation (PDG 2024, multivector_kernel_fit/01_findings.md). -/
def Koide_empirical_approx : Rat := 2 / 3

/-- Structural value: Q = dim G_2 / dim Spin(7) = 2/3 (from step 10). -/
def Koide_structural : Rat := koide_Q

/-- THEOREM: The structural Koide value equals 2/3. -/
theorem koide_structural_value : Koide_structural = 2/3 := by
  unfold Koide_structural
  exact koide_Q_value

/-- THEOREM: The structural and empirical Koide values agree (modulo the
    experimental m_tau uncertainty, which we encode as exact 2/3 here). -/
theorem koide_structural_matches_empirical :
    Koide_structural = Koide_empirical_approx := by
  unfold Koide_structural Koide_empirical_approx
  exact koide_Q_value

/-- Empirical Brannen phase ≈ 2/9 rad (within m_tau uncertainty). -/
def Brannen_empirical_approx : Rat := 2 / 9

/-- Structural value: phi = Q / 3 (from step 10). -/
def Brannen_structural : Rat := brannen_phi

/-- THEOREM: The structural Brannen value equals 2/9. -/
theorem brannen_structural_value : Brannen_structural = 2/9 := by
  unfold Brannen_structural
  exact brannen_phi_value

/-- THEOREM: phi = Q / 3 -- the structural identification. -/
theorem brannen_eq_koide_div_three :
    Brannen_structural = Koide_structural / 3 := by
  unfold Brannen_structural Koide_structural
  exact brannen_eq_koide_div_gens

/-- THEOREM: The structural and empirical Brannen values agree. -/
theorem brannen_matches : Brannen_structural = Brannen_empirical_approx := by
  unfold Brannen_structural Brannen_empirical_approx
  exact brannen_phi_value

/-- Cross-sector ratio = dim Spin(7) = 21 (from step 9). -/
def CrossSectorRatio_structural : Nat := cross_sector_ratio

/-- THEOREM: The cross-sector ratio is 21. -/
theorem cross_sector_is_21 : CrossSectorRatio_structural = 21 := by
  unfold CrossSectorRatio_structural cross_sector_ratio
  exact dim_Spin7_is_21

/-- The empirical cross-sector ratio is 20.945 (within 0.26% of 21). -/
def CrossSectorRatio_empirical_relative : Rat := 20945 / 1000

/-- THEOREM: The structural value (21) is close to the empirical 20.945
    (within the modulation-function precision; not exact). -/
example : (CrossSectorRatio_structural : Rat) - CrossSectorRatio_empirical_relative
    = 55 / 1000 := by
  unfold CrossSectorRatio_structural cross_sector_ratio dim_Spin7
        CrossSectorRatio_empirical_relative
  norm_num

end SCPv59.KoideAndBrannen
