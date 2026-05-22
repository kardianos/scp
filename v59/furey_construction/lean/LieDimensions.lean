/-
  v59/furey_construction/lean/LieDimensions.lean

  Formal statements of the exceptional Lie group dimensions and the
  v59 structural identifications:
    Q = dim G_2 / dim Spin(7) = 14/21 = 2/3
    phi = Q / 3 = 2/9
    cross-sector ratio = dim Spin(7) = 21
-/

import Mathlib.Tactic.NormNum
import Mathlib.Data.Rat.Defs

namespace SCPv59

/-- Dimension of the exceptional Lie group G_2 (automorphism group of octonions). -/
def dim_G2 : Nat := 14

/-- Dimension of Spin(7) (symmetry of the octonion imaginary sector). -/
def dim_Spin7 : Nat := 21

/-- Dimension of Spin(8) (the group with triality). -/
def dim_Spin8 : Nat := 28

/-- Number of fermion generations (Z_3 ⊂ S_3 triality of Spin(8)). -/
def n_generations : Nat := 3

/-- Cross-sector ratio appearing in v59 step 9 = dim Spin(7). -/
def cross_sector_ratio : Nat := dim_Spin7

/-- The Koide ratio Q = dim G_2 / dim Spin(7) as a rational number. -/
def koide_Q : Rat := (dim_G2 : Rat) / (dim_Spin7 : Rat)

/-- The Brannen phase phi = Q / 3 as a rational number. -/
def brannen_phi : Rat := koide_Q / (n_generations : Rat)

/-- THEOREM: dim Spin(7) factors as 3 * 7. -/
theorem dim_Spin7_factors : dim_Spin7 = 3 * 7 := by
  unfold dim_Spin7
  rfl

/-- THEOREM: dim Spin(7) is 21 (sanity check). -/
theorem dim_Spin7_is_21 : dim_Spin7 = 21 := by
  unfold dim_Spin7
  rfl

/-- THEOREM: Koide ratio Q = 14/21 = 2/3 as a rational. -/
theorem koide_Q_value : koide_Q = (2 : Rat) / 3 := by
  unfold koide_Q dim_G2 dim_Spin7
  norm_num

/-- THEOREM: Brannen phase phi = 2/9 as a rational. -/
theorem brannen_phi_value : brannen_phi = (2 : Rat) / 9 := by
  unfold brannen_phi koide_Q dim_G2 dim_Spin7 n_generations
  norm_num

/-- THEOREM: phi = Q / 3 (the structural identity of step 10). -/
theorem brannen_eq_koide_div_gens : brannen_phi = koide_Q / 3 := by
  unfold brannen_phi n_generations
  norm_num

/-- THEOREM: dim G_2 + dim Spin(7) = 35 (the larger Spin(8) chain). -/
theorem dim_chain_sum : dim_G2 + dim_Spin7 = 35 := by
  unfold dim_G2 dim_Spin7
  rfl

/-- THEOREM: dim Spin(8) - dim Spin(7) = 7 (the dim of the homogeneous space S^7). -/
theorem dim_S7 : dim_Spin8 - dim_Spin7 = 7 := by
  unfold dim_Spin8 dim_Spin7
  rfl

/-- THEOREM: dim Spin(7) - dim G_2 = 7 (the dim of S^7 = Spin(7)/G_2). -/
theorem dim_Spin7_minus_G2 : dim_Spin7 - dim_G2 = 7 := by
  unfold dim_Spin7 dim_G2
  rfl

end SCPv59
