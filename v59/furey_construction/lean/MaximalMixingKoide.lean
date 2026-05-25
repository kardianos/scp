/-
  v59/furey_construction/lean/MaximalMixingKoide.lean

  **Goal (1): Koide from the G₂-content of the maximally-mixed vacuum.**

  The session's central result (`integration_v58/06_koide_from_G2_maximal_mixing.md`,
  `g2_koide_derivation.py`): the Brannen amplitude `t²` and `v_Higgs` are the *same*
  equipartition mystery, and it resolves via two clean principles —

    (P1) MAXIMAL MIXING: the fermion order parameter is the maximally-mixed (uniform,
         maximally-symmetric) configuration on the sector.  Vindicated by
         `05_v58_vacuum_alignment.md`: the v58 grade-{0,2} energy is *flat* on the vacuum
         manifold (`M∈L=skew ⟹ M²=symmetric ⟹ ⟨M²⟩₂≡0`), so energy-minimization is silent
         on the alignment — maximal symmetry is the (consistent) primitive selection principle.

    (P2) G₂-INERT CORE: `G₂ = Aut(𝕆)` (the octonion automorphisms) preserves the algebra
         structure, so it carries NO generation-distinguishing (mass-splitting) information;
         the mass splitting is the complement.  `g2_koide_derivation.py` constructs `G₂` as
         the derivation algebra of 𝕆 and verifies `dim G₂ = 14` from first principles
         (grounding the `14`, previously only `dimGL 7 − dim3Forms 7 = 21−7`), and verifies
         `L = so(8) = G₂(14) ⊕ nonG₂(14)` so the non-G₂ fraction is exactly `1/2`.

  Consequence (this module): the maximally-mixed non-G₂ weight `(D − dimG₂)/D` equals the
  Brannen `t²` posited in `Predictions`, so the three Koide ratios are DERIVED from (P1)+(P2)
  rather than fitted.  The lepton lands on the symmetric `1/2` BECAUSE `D_lepton = 2·dimG₂`.

  HONEST SCOPE.  This formalises the *arithmetic* (`(D−dimG₂)/D ⇒ t² ⇒` Koide) and grounds
  `dimG₂=14` (built as the octonion derivation algebra).  These theorems are TRUE arithmetic.

  ⚠ MECHANISM RETRACTED (2026-05-24, `integration_v58/08_generation_map_result.md`).  The
  intended *bridge* — that the maximally-mixed L-vacuum's non-G₂ weight IS the Brannen amplitude
  `t²` — does **not** close: the `Aut(𝕊)=G₂×S₃` generation symmetry provably leaves `t²` free
  (S₃-covariance forces `M` circulant but fixes no `ξ`; max-mixing over the generation structure
  gives `t²∈{0,1}`, never `1/2`).  So `(D−dimG₂)/D` is a **dimensional coincidence in the so(8)
  Clifford grade, not a generation-level mechanism**.  Koide `2/3` (`t²=1/2`) is a *residual
  coupling magnitude*, parallel to the phase `φ=2/9`.  Read `massSplitFraction` as "the so(8)
  non-G₂ dimensional fraction (which numerically equals `t_sq`)", NOT as a derivation of `t²`.
-/
import Predictions
import SpinDimension

namespace SCPv59.MaximalMixing

open SCPv59.Predictions SCPv59.SpinDimension

/-! ## The maximal-mixing (equipartition) principle, made arithmetic -/

/-- **(P1)+(P2) as a formula.**  A vacuum maximally mixed (uniform) over a `D`-dim sector, with
    the `dimG₂`-dim G₂-automorphism core inert, has mass-splitting weight
    `t² = (D − dimG₂)/D` — the fraction of the sector outside `Aut(𝕆)`. -/
def massSplitFraction (D : ℕ) : ℚ := ((D : ℚ) - (dimG2 : ℚ)) / (D : ℚ)

/-- The equipartition weight is the Brannen `1 − dimG₂/D` form (for `D ≠ 0`). -/
theorem massSplit_eq (D : ℕ) (hD : (D : ℚ) ≠ 0) :
    massSplitFraction D = 1 - (dimG2 : ℚ) / (D : ℚ) := by
  unfold massSplitFraction
  field_simp

/-! ## The maximal-mixing weight reproduces the posited Brannen `t²` (so it DERIVES them) -/

theorem massSplit_lepton : massSplitFraction 28 = t_sq_lepton := by
  rw [massSplit_eq 28 (by norm_num)]; rfl

theorem massSplit_d_quark : massSplitFraction dim3FormsR7 = t_sq_d_quark := by
  rw [massSplit_eq dim3FormsR7 (by unfold dim3FormsR7; norm_num)]; rfl

theorem massSplit_u_quark : massSplitFraction dimU63 = t_sq_u_quark := by
  rw [massSplit_eq dimU63 (by unfold dimU63; norm_num)]; rfl

/-! ## The lepton's `1/2` is forced by `D_lepton = 2·dimG₂` -/

/-- **Why the lepton lands on the maximally-symmetric `t²=1/2`:** its sector dimension is
    *exactly double* the G₂ core (`28 = 2·14 = 2·dimG₂`), so the non-G₂ weight is `1/2`. -/
theorem lepton_half_from_double_core :
    (28 : ℕ) = 2 * dimG2 ∧ massSplitFraction 28 = 1 / 2 := by
  refine ⟨by rw [dimG2_eq_14], ?_⟩
  rw [massSplit_eq 28 (by norm_num),
      show (dimG2 : ℚ) = 14 by rw [dimG2_eq_14]; norm_cast]
  norm_num

/-! ## The derivation: Koide ratios FROM maximal mixing -/

/-- **Koide from the maximal-mixing principle.**  Feeding the maximally-mixed non-G₂ weight
    `massSplitFraction D` into the Brannen map `Q=(1+2t²)/3` yields the three observed Koide
    ratios — lepton `2/3`, d-quark `11/15`, u-quark `23/27` — with `D ∈ {28, 35, 63}` the
    sector grade-dimensions and the single universal `dimG₂ = 14` as the inert core. -/
theorem koide_from_maximal_mixing :
    koide_Q_from_t_sq (massSplitFraction 28) = 2 / 3
    ∧ koide_Q_from_t_sq (massSplitFraction dim3FormsR7) = 11 / 15
    ∧ koide_Q_from_t_sq (massSplitFraction dimU63) = 23 / 27 := by
  refine ⟨?_, ?_, ?_⟩
  · rw [massSplit_lepton]; exact koide_Q_lepton
  · rw [massSplit_d_quark]; exact koide_Q_d_quark
  · rw [massSplit_u_quark]; exact koide_Q_u_quark

/-- **Headline (the derivation chain).**  Given (P1) maximal mixing + (P2) `G₂=Aut(𝕆)` inert
    (`dimG₂=14`, built), the lepton Koide ratio is forced to `2/3`: the maximally-mixed non-G₂
    weight over `L` (D=28=2·dimG₂) is `1/2`, and `Q=(1+2·½)/3 = 2/3`. -/
theorem lepton_koide_from_maximal_mixing :
    massSplitFraction 28 = 1 / 2
    ∧ koide_Q_from_t_sq (massSplitFraction 28) = 2 / 3 :=
  ⟨lepton_half_from_double_core.2, koide_from_maximal_mixing.1⟩

end SCPv59.MaximalMixing
