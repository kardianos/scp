/-
  v59/furey_construction/lean/ChiralPhaseWindow.lean

  **M4 (hierarchy-from-phase / chiral protection) for the Brannen phase, formalized.**

  Pursuing the M4 lead (`koide_phase_law/physical_mechanisms.md`): the Brannen phase `φ` controls
  the *lightest* lepton mass.  In the Koide-normalised amplitude `√mₙ/a = 1 + √2·cos(φ + 2πn/3)`
  (`t²=½`, the proven amplitude), the electron generation is near-massless, and:

    * **chiral point:** the electron amplitude *vanishes* at `φ = π/12`
      (`1 + √2·cos(π/12 + 2π/3) = 0`, since the angle is `3π/4` and `cos(3π/4) = −√2/2`);
    * **physical phase is just inside the window:** `φ = 2/9 < π/12`, so the electron is **light
      but massive** (`1 + √2·cos(2/9 + 2π/3) > 0`).

  So the physical `φ = 2/9` sits a tiny gap `π/12 − 2/9 ≈ 0.04` below the chiral (electron-massless)
  point — "chiral protection keeps the electron light."

  HONEST SCOPE.  This is a **constraint / consistency relation**, NOT a derivation of `φ = 2/9`.
  Chiral protection explains why `φ` lies *near* `π/12` (small `m_e`), and pins the *window*; but the
  specific value `2/9 = Q/3` is the separate (unexplained) Brannen phase law, and the offset
  `π/12 − 2/9` is fixed by the measured `m_e` (an input), not derived.  `π/12` is π-rational
  (geometric); `2/9` is not (`PhaseExclusions.koide_not_pi_rational`) — so `2/9` is genuinely shifted
  off the chiral angle by a non-geometric amount.  What is formalised here is exactly the M4
  *window*: the chiral massless point and that the physical phase lies just inside it.
-/
import Mathlib

namespace SCPv59.ChiralPhaseWindow

open Real

/-- **The chiral (electron-massless) point.**  With the Koide amplitude `√2`, the electron Brannen
    amplitude `1 + √2·cos(φ + 2π/3)` vanishes at `φ = π/12` (the angle becomes `3π/4`,
    `cos(3π/4) = −√2/2`).  This is where chiral protection would make `m_e = 0`. -/
theorem chiral_massless_point :
    1 + Real.sqrt 2 * Real.cos (Real.pi / 12 + 2 * Real.pi / 3) = 0 := by
  have h2 : Real.sqrt 2 * Real.sqrt 2 = 2 := Real.mul_self_sqrt (by norm_num)
  rw [show Real.pi / 12 + 2 * Real.pi / 3 = Real.pi - Real.pi / 4 by ring,
      Real.cos_pi_sub, Real.cos_pi_div_four]
  have hkey : Real.sqrt 2 * -(Real.sqrt 2 / 2) = -1 := by
    rw [mul_neg, ← mul_div_assoc, h2]; norm_num
  rw [hkey]; norm_num

/-- **The physical phase is below the chiral point:** `2/9 < π/12`.  (Equivalent to `π > 8/3`.)
    So `φ = 2/9` lies *inside* the light-electron window — the electron is light but not massless. -/
theorem physical_phase_below_chiral : (2 : ℝ) / 9 < Real.pi / 12 := by
  linarith [Real.pi_gt_three]

/-- **The electron is light but massive at `φ = 2/9`:** the amplitude is strictly positive.
    Since `2/9 < π/12`, the electron angle `2/9 + 2π/3 < 3π/4`, and `cos` is strictly decreasing on
    `[0,π]`, so `cos(2/9+2π/3) > cos(3π/4) = −√2/2`, giving `1 + √2·cos(2/9+2π/3) > 0`. -/
theorem electron_light_but_massive :
    0 < 1 + Real.sqrt 2 * Real.cos (2 / 9 + 2 * Real.pi / 3) := by
  have hpi := Real.pi_gt_three
  have hab : (2 : ℝ) / 9 + 2 * Real.pi / 3 < 3 * Real.pi / 4 := by linarith
  have hbpi : (3 : ℝ) * Real.pi / 4 ≤ Real.pi := by linarith [Real.pi_pos]
  have ha_mem : (2 : ℝ) / 9 + 2 * Real.pi / 3 ∈ Set.Icc 0 Real.pi :=
    ⟨by positivity, le_of_lt (lt_of_lt_of_le hab hbpi)⟩
  have hb_mem : (3 : ℝ) * Real.pi / 4 ∈ Set.Icc 0 Real.pi :=
    ⟨by positivity, hbpi⟩
  have hcos : Real.cos (3 * Real.pi / 4) < Real.cos (2 / 9 + 2 * Real.pi / 3) :=
    Real.strictAntiOn_cos ha_mem hb_mem hab
  have hc34 : Real.cos (3 * Real.pi / 4) = -(Real.sqrt 2 / 2) := by
    rw [show (3 : ℝ) * Real.pi / 4 = Real.pi - Real.pi / 4 by ring,
        Real.cos_pi_sub, Real.cos_pi_div_four]
  rw [hc34] at hcos
  have h2 : Real.sqrt 2 * Real.sqrt 2 = 2 := Real.mul_self_sqrt (by norm_num)
  have hsqrt : (0 : ℝ) < Real.sqrt 2 := Real.sqrt_pos.mpr (by norm_num)
  nlinarith [hcos, h2, hsqrt]

/-- **The M4 window, bundled.**  Chiral massless point at `π/12`; physical phase `2/9` strictly
    below it; electron light but massive there.  (A *constraint*, not a derivation of `2/9` — see
    the module header.) -/
theorem m4_chiral_window :
    (1 + Real.sqrt 2 * Real.cos (Real.pi / 12 + 2 * Real.pi / 3) = 0)
    ∧ (2 : ℝ) / 9 < Real.pi / 12
    ∧ (0 < 1 + Real.sqrt 2 * Real.cos (2 / 9 + 2 * Real.pi / 3)) :=
  ⟨chiral_massless_point, physical_phase_below_chiral, electron_light_but_massive⟩

end SCPv59.ChiralPhaseWindow
