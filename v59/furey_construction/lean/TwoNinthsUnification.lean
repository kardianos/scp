/-
  v59/furey_construction/lean/TwoNinthsUnification.lean

  **Tier 3 of the φ = 2/9 program**: the value-level unification.

  The same number `2/9` is reached by two structurally independent routes:

    * **mass sector**  — the Brannen lepton phase `φ = Q/3`, `Q = dimG₂/dimSpin7 = 14/21`
      (`BrannenPhase`, `LieDimensions`):           `2/9 = 14/63 = dimG₂/(dimSpin7·3)`.
    * **gauge sector** — the weak mixing angle from Pati-Salam matching of `(cW,cBL)=(5,2)·√α`
      (`WeinbergPatiSalam`):                        `2/9 = (10/7)/(5+10/7) = 10/45`.
      (Complement: `cos²θ_W = 7/9 = dimImO/9`, `ScaleBridge`.)

  Here we prove the two routes give the *same* `2/9` (`mass_eq_gauge`), and record both
  structural decompositions.

  **(T3.1) OPEN**: whether this is a deep identity (the mass-generation phase *is* the EW
  mixing angle, sharing one structural origin) or a numerical coincidence at fit precision is
  unresolved — see the closing note.  What is proven is that they are numerically identical and
  each has its own structural derivation.
-/
import Mathlib
import LieDimensions
import WeinbergPatiSalam

namespace SCPv59.TwoNinths

open SCPv59

/-! ## The two routes -/

/-- **Mass-sector `2/9`**: the Brannen phase `Q/3`, with `Q = dimG₂/dimSpin7 = 2/3`. -/
theorem mass_sector : (koide_Q / 3 : ℚ) = 2 / 9 := by rw [koide_Q_value]; norm_num

/-- **Gauge-sector `2/9`**: `sin²θ_W`, derived from the Pati-Salam matching of `(cW,cBL)=(5,2)`. -/
theorem gauge_sector : WeinbergPatiSalam.sin2_thetaW = 2 / 9 := WeinbergPatiSalam.sin2_thetaW_eq

/-- **The unification (value level).**  Mass and gauge sectors yield the *same* `2/9`:
    `Q/3 = sin²θ_W`. -/
theorem mass_eq_gauge : (koide_Q / 3 : ℚ) = WeinbergPatiSalam.sin2_thetaW := by
  rw [mass_sector, gauge_sector]

/-! ## Structural decompositions of `2/9` -/

/-- Mass-sector decomposition: `2/9 = 14/63 = dimG₂ / (dimSpin7 · 3)` (Koide ratio over the
    three generations). -/
theorem mass_decomposition : (koide_Q / 3 : ℚ) = 14 / 63 := by rw [koide_Q_value]; norm_num

/-- Gauge-sector decomposition: `2/9 = (9 − 7)/9 = (9 − dimImO)/9` (complement of
    `cos²θ_W = dimImO/9 = 7/9`). -/
theorem gauge_decomposition : WeinbergPatiSalam.sin2_thetaW = ((9 : ℚ) - 7) / 9 := by
  rw [gauge_sector]; norm_num

/-- Both decompositions are the same rational, exhibited side by side: `dimG₂/(dimSpin7·3)`
    (mass) `= (9−dimImO)/9` (gauge) `= 2/9`. -/
theorem two_ninths_both_ways :
    (koide_Q / 3 : ℚ) = 14 / 63 ∧ WeinbergPatiSalam.sin2_thetaW = ((9 : ℚ) - 7) / 9
    ∧ (koide_Q / 3 : ℚ) = WeinbergPatiSalam.sin2_thetaW :=
  ⟨mass_decomposition, gauge_decomposition, mass_eq_gauge⟩

end SCPv59.TwoNinths
