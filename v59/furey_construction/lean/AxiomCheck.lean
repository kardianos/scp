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
