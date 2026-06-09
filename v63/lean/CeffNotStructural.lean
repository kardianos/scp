/-
  v63/lean/CeffNotStructural.lean

  MACHINE-CHECKED FALSIFICATION of the "c_eff = v_min/c = 2/3" gating lead
  (see ../CEFF_TEST.md, ../ceff_test.py).

  The localization/holonomy idea would have the Brannen phase be exp(i ∮A) with
  ∮A = 2/3 supplied by a density-dependent meta-constant: 2/3 = c_eff·(weight).
  The proposed anchor was the BLV null-rotor core velocity ratio v_min/c ≈ 0.670.

  For c_eff to be a *structural, parameter-free* constant equal to 2/3, it must
  be the SAME for every soliton (it cannot depend on the pion mass).  The
  reproduced data (v3/bin/nullrotor_metric) says otherwise:
      massless soliton:  v_min/c = 0.67020151
      massive  soliton:  v_min/c = 0.62519338
  They differ, so v_min/c is density-dependent — it is not a structural constant,
  and in particular not a parameter-free 2/3.  This file proves that modus tollens.

  Build (against the v59 Mathlib, as in v61/lean/EwVevHome.lean):
    cd v59/furey_construction/lean && lake env lean ../../../v63/lean/CeffNotStructural.lean
-/
import Mathlib

namespace SCPv63.CeffNotStructural

/-- BLV core velocity ratio `v_min/c`, MASSLESS B=1 soliton
    (reproduced from `v3/bin/nullrotor_metric -profile profile_sigma_e1.dat`). -/
def vmin_massless : ℚ := 67020151 / 100000000

/-- BLV core velocity ratio `v_min/c`, MASSIVE soliton
    (`profile_massive_e1_mpi0.398.dat`). -/
def vmin_massive : ℚ := 62519338 / 100000000

/-! ## 1. The data: parameter dependence and the gap to 2/3 -/

/-- **`v_min/c` is density-dependent**: the massless and massive solitons do NOT
    share it (the pion mass moves it `0.670 → 0.625`). -/
theorem vmin_param_dependent : vmin_massless ≠ vmin_massive := by
  unfold vmin_massless vmin_massive; norm_num

/-- The massless value is **not** `2/3` — it misses by more than `3/1000`
    (`≈ 0.53 %`); it is a numerical profile minimum, not a structural `2/3`. -/
theorem massless_ne_two_thirds : |vmin_massless - (2 / 3 : ℚ)| > 3 / 1000 := by
  unfold vmin_massless
  rw [abs_of_nonneg (by norm_num)]
  norm_num

/-- The massive value misses `2/3` by more than `6 %`. -/
theorem massive_ne_two_thirds : |vmin_massive - (2 / 3 : ℚ)| > 4 / 100 := by
  unfold vmin_massive
  rw [abs_of_nonpos (by norm_num), neg_sub]
  norm_num

/-! ## 2. The falsification (modus tollens) -/

/-- **No structural `c_eff`.**  If some constant `κ` were the structural,
    density-INDEPENDENT `c_eff` realized by both solitons, the two solitons would
    share it.  They do not — so no such `κ` exists. -/
theorem no_structural_ceff (κ : ℚ)
    (h_massless : κ = vmin_massless) (h_massive : κ = vmin_massive) : False :=
  vmin_param_dependent (h_massless ▸ h_massive)

/-- **`2/3` in particular cannot be the structural `c_eff`.**  (Immediate from
    `no_structural_ceff`; the holonomy anchor `2/3 = c_eff` is excluded — and
    independently, neither value is even close, by §1.) -/
theorem two_thirds_not_structural_ceff :
    ¬ ((2 / 3 : ℚ) = vmin_massless ∧ (2 / 3 : ℚ) = vmin_massive) := by
  rintro ⟨h0, hM⟩
  exact no_structural_ceff (2 / 3) h0 hM

/-! ## 3. Contrast: the one EXACT structural constant in this sector is `2` -/

/-- The BLV `P/m` ratio for `K=0` modes is exactly `2` (machine precision in
    `gravity-null-rotor-metric.md`); the clean algebraic constant of the sector
    is `2`, not `2/3`. -/
def PoverM_K0 : ℚ := 2
theorem structural_constant_is_two : PoverM_K0 = 2 := rfl

/-- **Headline.**  `v_min/c` is density-dependent (`massless ≠ massive`), each
    value misses `2/3`, and no constant `κ` can be a structural parameter-free
    `c_eff`.  The holonomy lead `2/3 = c_eff` gets no parameter-free anchor from
    the BLV metric. -/
theorem ceff_lead_falsified :
    vmin_massless ≠ vmin_massive
    ∧ |vmin_massless - (2 / 3 : ℚ)| > 3 / 1000
    ∧ ¬ ((2 / 3 : ℚ) = vmin_massless ∧ (2 / 3 : ℚ) = vmin_massive) :=
  ⟨vmin_param_dependent, massless_ne_two_thirds, two_thirds_not_structural_ceff⟩

end SCPv63.CeffNotStructural
