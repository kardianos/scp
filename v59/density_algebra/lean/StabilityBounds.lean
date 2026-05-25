/-
  DensityAlgebra.StabilityBounds

  This module aims to state, in a machine-checkable way, the key mathematical
  bounds that force the observed discrete structures in the v59 algebraic
  patterns (specific |ξ|² values, L/F graded selections, protection projectors,
  the 21-factor, etc.).

  The bounds come primarily from the second-variation (Hessian) analysis of the
  effective density/ potential in the living candidate, when protection is
  modeled as a projector onto subsets of the 8-dimensional (octonion-like)
  space.

  Current philosophy (May 2026):
  - Full formalization of the 8D non-associative octonion multiplication table
    and its Hessian is heavy. We therefore proceed in layers:
      1. Precise axiomatic statements of the forcing relationships.
      2. Computable numerical checks for concrete (λ, μ, mask) tuples.
      3. Gradual internalization of the algebraic identities that the Maxima
         sensitivity analysis derived (anisotropic curvature, 7-fold degeneracy,
         radial vs transverse modes, etc.).

  This file is intended to grow into the place where we can eventually
  *compute* critical values (e.g. the f_amplitude crossover, the minimal μ/λ
  ratio for a given mask, the exact |ξ|² that makes the radial eigenvalue zero)
  inside Lean.
-/

-- Lightweight: no Mathlib for Phase 1. Projector and matrix ops deferred to Phase 2.
-- All core logic uses List Bool masks and explicit eigenvalue lists (diagonal model).

namespace DensityAlgebra

-- Lightweight numeric proxy (Phase 1): use Float for concrete checks and #eval.
-- Full mathematical ℝ + inequalities in Prop will use Mathlib later.
abbrev ℝ := Float

/-! ## 1. Basic types and protection masks -/

abbrev Dim8 := Fin 8

-- Index convention matching the Python/ga and Maxima work:
-- 0 : scalar
-- 1,2,3 : vectors (L-content)
-- 4,5,6 : bivectors (F-content)
-- 7 : trivector (F-content)

def L_mask : List Bool := [true, true, true, true, false, false, false, false]
def F_mask : List Bool := [true, false, false, false, true, true, true, true]
def LF_mask : List Bool := [true, true, true, true, true, true, true, true]
def Full_mask : List Bool := [true, true, true, true, true, true, true, true]

-- Projector representation (List-based for lightweight Phase 1; matrix form in Phase 2)
def mask_to_list (mask : List Bool) : List Bool := mask  -- identity for now; later can be more structured

/-! ## 2. The effective quadratic form (from Maxima sensitivity analysis)

From the octonion_sensitivity_analysis.mac and the living-candidate model,
the second variation of the effective potential around a protected background
takes (to first approximation) the anisotropic diagonal form:

  H_{ii} ≈ 2(μ + λ)   for the scalar/radial direction (i=0)
  H_{ii} ≈ 2(μ - λ)   for the "imaginary" directions that are *inside* the projector
  H_{ii} ≈ 0          (or very small) for directions killed by the projector

The 7-fold (or more generally the multiplicity) degeneracy of the (μ-λ) mode
is the algebraic origin of the "extra internal twists" counted by 7/21/etc.

A configuration is stable (positive-definite well) only when the eigenvalues
on the image of the projector are all positive for the chosen λ, μ.
-/

structure EffectiveHessian where
  lam : ℝ
  mu : ℝ
  scalar_eigenvalue : ℝ := 2 * (mu + lam)
  transverse_eigenvalue : ℝ := 2 * (mu - lam)
  -- (In the full 8D analysis the transverse space splits according to the Fano
  -- structure; we keep the simple model here as the leading term.)

/-- Returns the list of eigenvalues (with multiplicity) for a given mask.
    This is the computable heart of the bound. -/
def eigenvalues_for_mask (lam mu : ℝ) (mask : List Bool) : List ℝ :=
  let h : EffectiveHessian := { lam := lam, mu := mu }
  let n_transverse := (mask.filter (· = true)).length - 1  -- subtract scalar
  [h.scalar_eigenvalue] ++ List.replicate n_transverse h.transverse_eigenvalue

/-- The stability condition: all relevant eigenvalues must be positive. -/
def is_stable (lam mu : ℝ) (mask : List Bool) : Prop :=
  ∀ ev ∈ eigenvalues_for_mask lam mu mask, ev > 0

/-! ## 3. Concrete computable bounds (the ones we want to prove / compute)

These are the statements that directly correspond to the numerical observations
from the f_amplitude sweeps and the Maxima Hessians.

Example (to be strengthened with actual numbers from the sweeps):
For λ = 0.005 and μ = 0.001, only the LF and Full masks (and sufficiently
rich L masks when F-contamination is low) satisfy the stability condition,
while pure-F is marginal or unstable in the radial direction.
-/

def example_stability_bound : Prop :=
  let lam := (0.005 : ℝ)
  let mu := (0.001 : ℝ)
  is_stable lam mu LF_mask ∧
  is_stable lam mu Full_mask ∧
  ¬ is_stable lam mu F_mask   -- or only marginally stable depending on exact model

/-! ## 4. The "why only these discrete values" statements

The observed |ξ|² = 1/2, 3/5, 7/9 are the only radii at which the radial
eigenvalue crosses zero while the transverse eigenvalues (inside the projector)
remain positive for the living-candidate quadratic coefficients that also
keep the force channels cleanly separated (small commutation error).

This is the Lean counterpart of the Maxima statement:
"the radial mode crosses zero only at the discrete sector values for the
observed λ/μ ratios."
-/

structure RadialZeroCondition (lam mu r : ℝ) : Prop where
  radial_ev : 2 * (mu + lam * r) = 0   -- simplified leading term
  transverse_positive : 2 * (mu - lam) > 0
  -- plus the force-separation / commutation bound (still axiomatic for now)

/-- The claim that only the three observed r values satisfy the above for
    the physically allowed range of lam, mu. -/
axiom OnlyThreeRadialZeros :
  ∀ lam mu, (0 < lam ∧ lam < 0.01) → (0 < mu ∧ mu < 0.005) →
    let candidates := [(1/2 : ℝ), 3/5, 7/9]
    ∀ r, RadialZeroCondition lam mu r → r ∈ candidates

/-! ## 5. Roadmap toward computation (updated 2026-05-24)

The real 8×8 matrix foundation now lives in `v59/furey_construction/lean/7D_Algebra/`
(see `SevenDAlgebra.lean`, `StabilityFromAlgebra.lean`, and `From7DAlgebra.lean` in this folder for the thin bridge).

- Replace the simple `EffectiveHessian` with the actual matrices and multiplication
  table from the 7D Algebra work.
- Use the explicit generator action on the |Ω_N⟩ states to derive (not assert)
  the stability conditions for the different projectors.
- The f_amplitude crossover and "only closed subalgebras" claims are now being
  moved from schematic modeling to derivation (see `StabilityFromAlgebra.lean`
  and the notes in the 7D_Algebra folder).

See `INTEGRATION_PLAN.md` in `v59/furey_construction/lean/7D_Algebra/` for the current
integration plan and status.
-/

/-! ## Phase 1 Additions: Computable numerical checks (Float)

For concrete (lam, mu, mask) we can now #eval the eigenvalue lists and
a decidable stability predicate. This reproduces the spirit of the Maxima
spectra (anisotropic: scalar 2(μ+λ), transverse 2(μ-λ)) and allows checking
the Python crossover conditions once we adjust for the sign regime where
μ > λ (transverse positive).
-/

def is_stable_dec (lam mu : ℝ) (mask : List Bool) : Bool :=
  (eigenvalues_for_mask lam mu mask).all (fun ev => ev > 0.0)

-- Example concrete spectra (use μ > λ to have positive transverse for demo;
-- real data from sweeps use regimes where effective signs allow L/LF stability).
#eval eigenvalues_for_mask 0.001 0.005 L_mask
#eval eigenvalues_for_mask 0.001 0.005 LF_mask
#eval is_stable_dec 0.001 0.005 L_mask
#eval is_stable_dec 0.001 0.005 F_mask

-- Concrete numerical example (Phase 1): for μ > λ the diagonal model makes
-- all masks with >=1 active "stable" (positive evals). The distinction between
-- L / F / LF requires the full octonion couplings (Phase 2).
def example_stability_bound_v2 : Prop :=
  let lam := (0.001 : ℝ)
  let mu := (0.005 : ℝ)
  is_stable_dec lam mu LF_mask = true ∧
  is_stable_dec lam mu L_mask = true ∧
  is_stable_dec lam mu F_mask = true   -- all positive in simplified anisotropic model when μ>λ

/-! ## Additional Phase 1 Angles (7+ total)

Angle 4: Exact Rat arithmetic version (for norm_num / decide later).
-/

abbrev Q := Rat

def eigenvalues_for_mask_Q (lam mu : Q) (mask : List Bool) : List Q :=
  -- Phase 1 stub for exact arithmetic angle: return positive values for any mask
  -- (full reimplementation of 2*(mu ± lam) in Rat + norm_num proofs in follow-up work)
  let n_true := (mask.filter (· = true)).length
  let n_trans := n_true - 1
  let scalar_q : Q := 2 * (mu + lam)   -- Rat arith works directly
  let trans_q : Q := 2 * (mu - lam)
  [scalar_q] ++ List.replicate n_trans trans_q

def is_stable_dec_Q (lam mu : Q) (mask : List Bool) : Bool :=
  (eigenvalues_for_mask_Q lam mu mask).all (fun ev => ev > 0)

-- Angle 6: Inductive structure for protection masks / allowed grade sets
inductive ProtTech
  | L
  | F
  | LF
  | Custom (bits : List Bool)

def tech_to_mask : ProtTech → List Bool
  | .L => L_mask
  | .F => F_mask
  | .LF => LF_mask
  | .Custom b => b

def is_stable_for_tech (lam mu : ℝ) (t : ProtTech) : Bool :=
  is_stable_dec lam mu (tech_to_mask t)

-- Angle 7: Stability certificate (data that Python/Maxima can export for Lean to verify)
structure StabilityCert (lam mu : ℝ) (mask : List Bool) where
  evals : List ℝ
  matches_model : evals = eigenvalues_for_mask lam mu mask
  all_positive : evals.all (· > 0)
  -- In full: + proof term or hash of source data for audit

-- Example cert construction (verifiable from numeric output). Use sorry for Phase 1
-- because Float literals and computed 2*(..) not defeq in this toolchain for rfl/decide.
def example_L_cert : StabilityCert 0.001 0.005 L_mask := sorry

-- Angle 1/4 strengthened axiomatic: concrete bound on critical ratio
-- For the diagonal model, transverse >0 iff mu > lam (i.e. ratio λ/μ <1 )
def critical_mu_over_lam_for_transverse_positive : Prop :=
  ∀ (lam mu : ℝ), (0 < lam ∧ 0 < mu) →
    (is_stable_dec lam mu L_mask = true) ↔ (mu > lam)

end DensityAlgebra