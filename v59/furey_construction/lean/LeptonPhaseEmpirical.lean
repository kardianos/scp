/-
  v59/furey_construction/lean/LeptonPhaseEmpirical.lean

  **POSITIVE knowledge on the lepton phase `φ = 2/9`** — the empirical side, complementing the
  negative/structural results of `PhaseAmbiguity` and `PhaseExclusions`.

  The Brannen phase extracted from the PDG-2024 charged-lepton masses (reduced to its principal
  `Z₃` branch) is `φ = 0.22222963 ± 8.4e-6`.  Two facts that elevate `φ = 2/9` from "loose
  coincidence" to a Koide-class regularity:

    (P1) `φ = 2/9` holds to relative `≈ 3.3e-5` — the SAME `~10⁻⁵ / 0.9σ` level as the Koide
         relation `Q = 2/3` itself (both limited by the `m_τ` uncertainty).  So the phase
         relation is as tight as the celebrated Koide identity, not a vague `~0.22` near-miss.

    (P2) The "/3" is pinned to the **number of generations**: among `φ = Q/n`, only `n = 3`
         matches the data (`Q/2, Q/4, Q` miss by `> 1e-2`).  With `Q = dimG₂/dimSpin7 = 2/3`
         structural, this *selects* `φ = Q/(#generations) = 2/9`.

  Net (with the negatives): `φ = Q/3` is a precise, lepton-specific empirical regularity — a
  "second Koide relation" — that has NO identified structural mechanism (its invariant `cos 3φ
  = cos(2/3)` is transcendental, non-geometric, not the Weinberg angle).  Positive: the pattern
  is real and tight.  Negative: the mechanism is genuinely absent, not merely unfound-among-easy.
-/
import Mathlib
import LieDimensions

namespace SCPv59.LeptonPhaseEmpirical

open SCPv59

/-- Measured charged-lepton Brannen phase (PDG 2024, reduced to the principal `Z₃` branch):
    `φ = 0.22222963 ± 8.4e-6`.  (Central value; the `±` is dominated by `m_τ = 1776.86 ± 0.12`.) -/
def phi_measured : ℚ := 22222963 / 100000000

/-- Measured Koide ratio (PDG 2024): `Q = 0.66666051 ± 6.8e-6`. -/
def Q_measured : ℚ := 66666051 / 100000000

/-! ## (P1) `φ = 2/9` holds at Koide precision -/

/-- **The phase agrees with `2/9` to better than `10⁻⁵` (absolute).**  `|2/9 − φ_meas| ≈ 7.4e-6`. -/
theorem phase_agrees_2_9 : |(2 / 9 : ℚ) - phi_measured| < 1 / 100000 := by
  unfold phi_measured; rw [abs_sub_comm]; rw [abs_of_nonneg (by norm_num)]; norm_num

/-- **Koide agrees with `2/3` to the same `10⁻⁵`** (`|2/3 − Q_meas| ≈ 6.2e-6`) — exhibited
    side by side to show the phase relation is *as tight as* the Koide relation. -/
theorem koide_agrees_2_3 : |(2 / 3 : ℚ) - Q_measured| < 1 / 100000 := by
  unfold Q_measured; rw [abs_of_nonneg (by norm_num)]; norm_num

/-- **Both agreements are the same order** — the phase relation `φ = 2/9` and the Koide relation
    `Q = 2/3` deviate from data by comparable `~10⁻⁵` amounts (both `m_τ`-limited, `≈ 0.9σ`). -/
theorem phase_as_tight_as_koide :
    |(2 / 9 : ℚ) - phi_measured| < 1 / 100000 ∧ |(2 / 3 : ℚ) - Q_measured| < 1 / 100000 :=
  ⟨phase_agrees_2_9, koide_agrees_2_3⟩

/-! ## (P2) The "/3" is the number of generations -/

/-- **`φ = Q/n` selects `n = 3` uniquely.**  With structural `Q = 2/3`: `Q/3 = 2/9` matches the
    data to `< 10⁻⁵`, while `Q/2, Q/4, Q/1` miss by `> 10⁻²`.  So the divisor is pinned to the
    number of generations, giving `φ = Q/(#generations) = 2/9`. -/
theorem generation_count_pinned :
    |(koide_Q / 3 : ℚ) - phi_measured| < 1 / 100000        -- n = 3 : matches (= #generations)
    ∧ |(koide_Q / 2 : ℚ) - phi_measured| > 1 / 100         -- n = 2 : off
    ∧ |(koide_Q / 4 : ℚ) - phi_measured| > 1 / 100         -- n = 4 : off
    ∧ |(koide_Q / 1 : ℚ) - phi_measured| > 1 / 100 := by   -- n = 1 : off
  rw [koide_Q_value]
  refine ⟨?_, ?_, ?_, ?_⟩
  · unfold phi_measured; rw [abs_sub_comm, abs_of_nonneg (by norm_num)]; norm_num   -- n=3
  · unfold phi_measured; rw [abs_of_nonneg (by norm_num)]; norm_num                  -- n=2
  · unfold phi_measured; rw [abs_sub_comm, abs_of_nonneg (by norm_num)]; norm_num   -- n=4
  · unfold phi_measured; rw [abs_of_nonneg (by norm_num)]; norm_num                  -- n=1

/-- The matched divisor equals `n_generations` (`= 3`), and `koide_Q / n_generations = 2/9`. -/
theorem phase_eq_Q_over_generations : (koide_Q / (n_generations : ℚ) : ℚ) = 2 / 9 := by
  rw [koide_Q_value, show ((n_generations : ℚ)) = 3 by norm_num [n_generations]]; norm_num

end SCPv59.LeptonPhaseEmpirical
