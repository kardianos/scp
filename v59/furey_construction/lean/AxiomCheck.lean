/- AxiomCheck.lean — list axioms used by the headline v59 theorems.
   `#print axioms` writes its output to the compile log.  We expect only
   `propext`, `Classical.choice`, `Quot.sound` (the Mathlib trio). -/

import BrannenKernel
import CyclicShift
import KoideAndBrannen
import SpinDimension
import SilentDirection
import Predictions
import EmbeddingIndex
import ScaleBridge
import BrannenPhase
import WeinbergPatiSalam
import TwoNinthsUnification
import PhaseExclusions
import PhaseAmbiguity
import LeptonPhaseEmpirical
import HiggsVevReframe
import MaximalMixingKoide
import OctoHalf
import ChiralPhaseWindow
import KernelEigenvalues
import XiVacuum
import CliffordBladeGrade
import GaugePrefactorDualCoxeter
import LeptonGradeForcing
import LeptonComplexStructure
import BladeSquareSign
import LeptonRealityForcing
import ColorSU3

-- Kernel diagonalisation (the matrix M actually has the Brannen amplitudes as eigenvalues)
#print axioms SCPv59.KernelEigenvalues.M_mulVec_eigen
#print axioms SCPv59.KernelEigenvalues.lam_eq_brannen
-- Universal grade law (disjoint bivectors → F), the abstract basis of the Z₂×Z₂ forcing
#print axioms SCPv59.CliffordBladeGrade.disjoint_bivectors_mul_isF
#print axioms SCPv59.CliffordBladeGrade.grade_bladeMul
-- Gauge prefactor 5 = dual Coxeter number h∨(Spin(7))
#print axioms SCPv59.GaugePrefactor.gW_prefactor_is_dualCoxeter_spin7
-- lepton=L NOT forced by mass-channel availability (the honest non-forcing result)
#print axioms SCPv59.Furey7D.LeptonForcing.lepton_L_not_forced_by_availability
#print axioms SCPv59.Furey7D.LeptonForcing.F_chi_is_chirality
-- lepton=L FORCED via the complex structure: J=e₀₁∈Λ²⊂L; L=Λ²⊕Λ⁶ are exactly the
-- complex structures (square −I), F=Λ⁴ the real structures (square +I)
#print axioms SCPv59.Furey7D.LeptonComplexStructure.lepton_complex_structure_forced_L
#print axioms SCPv59.Furey7D.LeptonComplexStructure.L_grade_squares_neg
#print axioms SCPv59.Furey7D.LeptonComplexStructure.F_grade_squares_pos
-- the representation-independent reason: simple even k-blade² = (−1)^{k(k+1)/2}
#print axioms SCPv59.BladeSquareSign.prod_sq
#print axioms SCPv59.BladeSquareSign.even_grade_complex_structure_dichotomy
-- (A) gap-closer + (B) reality: no symmetric matrix is a complex structure (whole F⊕scalar
-- subspace, not just blades); L=skew, {Λ⁰}⊕F=symmetric; orthogonal CS ⇒ skew ⇒ J∈L
#print axioms SCPv59.Furey7D.LeptonRealityForcing.complex_structure_shape_constraint
#print axioms SCPv59.Furey7D.LeptonRealityForcing.symm_not_complexStructure
#print axioms SCPv59.Furey7D.LeptonRealityForcing.orthogonal_complexStructure_skew
-- (C1)+(C2): explicit color su(3) (closes, kills lepton singlet, commutes with J_c=γ₀γ₅);
-- color-invariance pins the lepton complex structure to ±canonical
#print axioms SCPv59.Furey7D.ColorSU3.colorSU3_closes
#print axioms SCPv59.Furey7D.ColorSU3.color_commutes_Jc
#print axioms SCPv59.Furey7D.ColorSU3.color_kills_lepton
#print axioms SCPv59.Furey7D.ColorSU3.colorInvariant_pins_lepton
#print axioms SCPv59.Furey7D.ColorSU3.Jc_pinned
-- quark single orbit + full classification of the color-invariant complex structures
#print axioms SCPv59.Furey7D.ColorSU3.quark_single_orbit
#print axioms SCPv59.Furey7D.ColorSU3.colorInvariant_classification
-- holonomy obstruction: the color Z₃ fixes the lepton singlet ⇒ φ=2/9 is not an intra-ideal holonomy
#print axioms SCPv59.Furey7D.ColorSU3.colorZ3_fixes_lepton
#print axioms SCPv59.Furey7D.ColorSU3.colorZ3_commutes_Jc
-- φ=2/9 program: Tier 0 (phase free w.r.t. Koide; enters only via cos 3φ ⇒ target 3φ=Q)
#print axioms SCPv59.BrannenPhase.Q_phase_independent
#print axioms SCPv59.BrannenPhase.sum_s_cube
-- Tier 2 (sin²θ_W = 2/9 derived from Pati-Salam (5,2)) and Tier 3 (mass = gauge unification)
#print axioms SCPv59.WeinbergPatiSalam.sin2_thetaW_eq
#print axioms SCPv59.WeinbergPatiSalam.sin2_from_coeffs
#print axioms SCPv59.TwoNinths.mass_eq_gauge
-- exclusions narrowing the design: φ not π-rational (no geometric angle); 2/9 pins gauge cBL=2
#print axioms SCPv59.PhaseExclusions.phase_not_pi_rational
#print axioms SCPv59.PhaseExclusions.koide_not_pi_rational
#print axioms SCPv59.PhaseExclusions.gauge_cBL_pinned
#print axioms SCPv59.PhaseExclusions.cos6_potential_does_not_select_2_9
-- DECISIVE: masses fix φ only mod the generation S₃ (only cos 3φ); φ=2/9 ≡ 2/9+2π/3 (the fit value)
#print axioms SCPv59.PhaseAmbiguity.phase_2_9_not_unique
#print axioms SCPv59.PhaseAmbiguity.s_shift
#print axioms SCPv59.PhaseAmbiguity.invariant_is_cos_two_thirds
-- the natural deep identity cos3φ = cos²θ_W = 7/9 is FALSIFIED (cos(2/3) > 7/9)
#print axioms SCPv59.PhaseAmbiguity.phase_invariant_ne_cos_sq_thetaW
-- POSITIVE: φ=2/9 holds at Koide precision (~10⁻⁵); the "/3" is uniquely the generation count
#print axioms SCPv59.LeptonPhaseEmpirical.phase_as_tight_as_koide
#print axioms SCPv59.LeptonPhaseEmpirical.generation_count_pinned
-- integration #1: v_Higgs=28²a² ⟺ √v=dim(L)·a ⟺ Σ√m/√v = N_gen/dim(L) = 3/28
#print axioms SCPv59.HiggsVev.vHiggs_sqrt_form
#print axioms SCPv59.HiggsVev.ratio_form
-- (iii) backed into known physics: the relation IS a lepton Yukawa sum rule Σ√y=(3/28)2^{1/4}
#print axioms SCPv59.HiggsVev.yukawa_sum_form
-- goal (1): Koide from the G₂-content of the maximally-mixed vacuum.  t²=(D−dimG₂)/D as the
-- max-mixing non-G₂ weight ⇒ Koide 2/3,11/15,23/27; lepton ½ forced by D_lepton=2·dimG₂.
#print axioms SCPv59.MaximalMixing.koide_from_maximal_mixing
#print axioms SCPv59.MaximalMixing.lepton_koide_from_maximal_mixing
#print axioms SCPv59.MaximalMixing.lepton_half_from_double_core
#print axioms SCPv59.MaximalMixing.massSplit_lepton
-- octomath: the Koide 1/2 = root-product of the L-grade complex-structure half-element (1+u)/2.
-- u²=−1 (L, complex) ⇒ P²=P−1/2 (root-product 1/2); u²=+1 (F, real) ⇒ idempotent (root-product 0).
-- The 1/2 is the L-grade signature, fixed by the grade square-sign (not arbitrary algebra).
#print axioms SCPv59.OctoHalf.half_element_law
#print axioms SCPv59.OctoHalf.complex_half_field
#print axioms SCPv59.OctoHalf.real_half_field
#print axioms SCPv59.OctoHalf.root_product_complex
-- M4 chiral-phase window (the phase φ=2/9 lead): electron-massless point at π/12; physical 2/9
-- below it (light-electron window); electron light-but-massive.  A CONSTRAINT, not a derivation of 2/9.
#print axioms SCPv59.ChiralPhaseWindow.chiral_massless_point
#print axioms SCPv59.ChiralPhaseWindow.physical_phase_below_chiral
#print axioms SCPv59.ChiralPhaseWindow.electron_light_but_massive
#print axioms SCPv59.ChiralPhaseWindow.m4_chiral_window
-- ξ vacuum Goldstone structure (mass spectrum {λ,0,0,0})
#print axioms SCPv59.XiVacuum.xi_mass_spectrum
#print axioms SCPv59.XiVacuum.V_hessian_directional
#print axioms SCPv59.XiVacuum.hessian_eq_massMatrix_quadForm

#print axioms SCPv59.BrannenKernel.Q_value
#print axioms SCPv59.BrannenKernel.koide_iff_constraint
#print axioms SCPv59.BrannenKernel.Q_at_constraint
#print axioms SCPv59.BrannenKernel.sum_s
#print axioms SCPv59.BrannenKernel.sum_s_sq
#print axioms SCPv59.CyclicShift.ω_pow_three
#print axioms SCPv59.CyclicShift.sum_one_omega_omega_sq
#print axioms SCPv59.CyclicShift.ω_mul_ω_sq
#print axioms SCPv59.KoideAndBrannen.koide_structural_value
#print axioms SCPv59.KoideAndBrannen.brannen_structural_value
#print axioms SCPv59.SpinDimension.dimSO_eq_choose
#print axioms SCPv59.SpinDimension.dimSpin_seven
#print axioms SCPv59.SpinDimension.twenty_one_threefold
#print axioms SCPv59.SpinDimension.dimG2_eq_14
#print axioms SCPv59.SpinDimension.dimG2_via_S7
#print axioms SCPv59.SpinDimension.g2_spin7_spin8_chain
#print axioms SCPv59.SpinDimension.koide_ratio_structural
#print axioms SCPv59.SpinDimension.brannen_phase_structural
#print axioms SCPv59.SilentDirection.re_conj
#print axioms SCPv59.SilentDirection.normSq_conj
#print axioms SCPv59.SilentDirection.re_conj_unit
#print axioms SCPv59.SilentDirection.normSq_conj_unit
#print axioms SCPv59.SilentDirection.silent_pair
#print axioms SCPv59.SilentDirection.im_normSq_conj_unit
#print axioms SCPv59.SilentDirection.silent_orbit
#print axioms SCPv59.Predictions.killing_index_eq_dim_diff
#print axioms SCPv59.Predictions.gravity_prefactor_value
#print axioms SCPv59.Predictions.G_e_conjecture_form
#print axioms SCPv59.Predictions.g_W_squared_form
#print axioms SCPv59.Predictions.v59_prediction_tier_summary
#print axioms SCPv59.Predictions.v59_three_structural_integers
#print axioms SCPv59.Predictions.alpha_conjecture_S_em
#print axioms SCPv59.EmbeddingIndex.index_so3_so7
#print axioms SCPv59.EmbeddingIndex.index_so3_so8
#print axioms SCPv59.EmbeddingIndex.index_so3_so4
#print axioms SCPv59.EmbeddingIndex.killing_index_5_dual_form
#print axioms SCPv59.EmbeddingIndex.index_so3
#print axioms SCPv59.Predictions.t_sq_lepton_eq
#print axioms SCPv59.Predictions.t_sq_d_quark_eq
#print axioms SCPv59.Predictions.t_sq_u_quark_eq
#print axioms SCPv59.Predictions.koide_Q_lepton
#print axioms SCPv59.Predictions.koide_Q_d_quark
#print axioms SCPv59.Predictions.koide_Q_u_quark
#print axioms SCPv59.Predictions.v59_D_sum_identity
#print axioms SCPv59.Predictions.brannen_pattern_uniform
#print axioms SCPv59.Predictions.dim_cl7_even
#print axioms SCPv59.Predictions.lepton_ambient_decomp
#print axioms SCPv59.Predictions.d_quark_ambient_decomp
#print axioms SCPv59.Predictions.u_quark_ambient_decomp
#print axioms SCPv59.Predictions.single_source_decomposition
#print axioms SCPv59.ScaleBridge.dimLepton_eq_L
#print axioms SCPv59.ScaleBridge.dimLepton_eq_decomp
#print axioms SCPv59.ScaleBridge.v_Higgs_factor
#print axioms SCPv59.ScaleBridge.sin_sq_thW_eq_brannen_phase
#print axioms SCPv59.ScaleBridge.cos_sq_thW_eq_t_sq_u_quark
#print axioms SCPv59.ScaleBridge.sin_cos_sq_thW_sum
#print axioms SCPv59.ScaleBridge.mZ_over_mW_sq
#print axioms SCPv59.ScaleBridge.sqrt_alpha_MZ_form
#print axioms SCPv59.ScaleBridge.sqrt_alpha_MZ_factored
#print axioms SCPv59.ScaleBridge.alpha_MZ_from_consistency
#print axioms SCPv59.ScaleBridge.scale_bridge_summary
#print axioms SCPv59.ScaleBridge.koide_deviation_universal
#print axioms SCPv59.ScaleBridge.two_dimG2_eq_dimSpin8
#print axioms SCPv59.ScaleBridge.two_dimG2_eq_dimLepton
#print axioms SCPv59.ScaleBridge.dimImO_eq_choose_seven
#print axioms SCPv59.ScaleBridge.dimImO_eq_lambda6
#print axioms SCPv59.ScaleBridge.sin_sq_cabibbo_value
#print axioms SCPv59.ScaleBridge.cabibbo_seven_eq_cos_thW_seven
#print axioms SCPv59.ScaleBridge.spin7_pati_salam_decomp
#print axioms SCPv59.ScaleBridge.spin7_g2_plus_seven
#print axioms SCPv59.ScaleBridge.seven_unifies_su2L_su2R_u1BL_eq_dimImO
