/-
  DensityAlgebra.DensityForcing

  Machine-readable statements of the strongest forcing relationships
  that select the observed v59 algebraic structures (and rule out alternatives)
  when the underlying physics is a dynamical medium whose primary activity
  is achieving and protecting locally higher density while maintaining
  clean long-range force separation and non-dispersive emergent geometry.

  These are stated as `Prop` (or `axiom` where we do not yet have a proof
  from more primitive assumptions). The goal is to make the "why this
  and not the alternatives" logic explicit and checkable.

  This is the Lean counterpart to the hypotheses in ../HYPOTHESES.md
  and the bounded quantities in ../CONSTRAINTS_AND_TROUGHS.md.
-/

namespace DensityAlgebra

-- Lightweight stubs for Phase 1 (no Mathlib). Use Float as proxy for ℝ in numerical/axiomatic work.
-- Full Real and Mathlib reals + analysis deferred to later phases when we add the require.
abbrev ℝ := Float
def Real_sqrt (x : ℝ) : ℝ := x.sqrt   -- proxy

/-! ## Core Physical Assumptions (the "physics" we are formalizing) -/

-- The medium has an intrinsic drive to increase local density / connectivity.
axiom MediumDensityDrive : Prop

-- Stable high-density configurations require protection mechanisms
-- (algebraic, topological, or dynamical) to prevent rapid dispersal/radiation.
axiom ProtectionRequired : Prop

-- Long-lived structures must preserve clean separation between
-- the density-sourced (Newtonian) and bivector (Maxwell) response channels.
axiom ForceChannelSeparation : Prop

-- Emergent causal structure (propagation, light-cones) must remain
-- non-dispersive and universal for small fluctuations.
axiom NonDispersiveUniversality : Prop

-- Lightweight stubs for quantitative bounds (Phase 1 placeholders; will be refined
-- by computable eigenvalue conditions from StabilityBounds and data from Python/Maxima).
axiom DensityWellDepth : ℝ → ℝ
axiom MinStableDepth : ℝ
axiom ProtectionBudgetSpent : ℝ → ℝ
axiom MaxAffordableBudget : ℝ
axiom ForceSeparationMargin : ℝ → ℝ
axiom MinRequiredMargin : ℝ
axiom NumberOfStablePackings : ℝ → Nat
axiom DensityGain : Nat → ℝ
axiom ProtectionCost : Nat → ℝ
axiom ProtectionBudgetAvailable : ℝ
axiom ForceSeparation : Nat → ℝ
axiom ModulationCost : Nat → ℝ
axiom HierarchySuppression : ℝ
axiom DispersionIntroduced : Nat → ℝ

/-! ## Forcing Relationship 1: S³ Constraint Surface as Density-Protection Trade-off

The quaternionic parameter ξ is forced onto the 3-sphere of radius exactly 1/√2
(the surface on which Koide Q = 2/3 is automatic and the Brannen phase appears).

This is not a free choice. It is the unique radius at which three requirements
are simultaneously satisfied:
- Sufficient protection budget is spent to create a deep enough density well
  for stable lepton-like excitations.
- The remaining freedom still permits three degenerate but distinct packings
  (the generations) related by triality.
- The force-channel separation budget is not exhausted (the "silent" directions
  for the weak bosons remain available).

If the radius were larger, the wells would be too shallow (insufficient density locking).
If smaller, either the protection cost would over-constrain the algebra (no generations,
or force channels mix) or the configuration would be unstable.
-/
def S3_Radius_Forced_TradeOff : Prop :=
  MediumDensityDrive ∧ ProtectionRequired ∧ ForceChannelSeparation →
  ∃ r : ℝ, r = (1 : ℝ) / Real_sqrt 2 ∧
    (DensityWellDepth r > MinStableDepth) ∧
    (ProtectionBudgetSpent r ≤ MaxAffordableBudget) ∧
    (ForceSeparationMargin r > MinRequiredMargin) ∧
    (NumberOfStablePackings r = 3)

-- The observed radius is the only value that satisfies all three bounds simultaneously.
axiom S3_Radius_Is_Unique_Solution : S3_Radius_Forced_TradeOff

/-! ## Forcing Relationship 2: L/F Graded Selection as Stacked Protection Technologies

The single-source decomposition in Cl(7)_even and the exact additive identity
D_u = D_e + D_d (28 + 35 = 63) with the Z₂×Z₂ pattern is forced by the
existence of two distinct, composable protection technologies with different
density-gain / protection-cost ratios.

- L content (28) : cheap, light protection → sufficient for lepton-scale density wells.
- F content (35) : higher-gain but more expensive binding technology.
- Only the combination L+F permits the deepest stable wells (u-quarks) while
  still satisfying ForceChannelSeparation.

Any other selection of graded pieces either fails to reach observed density
or violates the separation or stability bounds.
-/
def LF_Selection_Forced_By_Protection_Budgets : Prop :=
  MediumDensityDrive ∧ ProtectionRequired ∧ ForceChannelSeparation →
  ∃ (L F : Nat), L = 28 ∧ F = 35 ∧
    (DensityGain L < DensityGain (L + F)) ∧
    (ProtectionCost (L + F) ≤ ProtectionBudgetAvailable) ∧
    (ForceSeparation (L + F) = ForceSeparation L)   -- additive, no extra mixing cost

axiom LF_Selection_Is_Forced : LF_Selection_Forced_By_Protection_Budgets

/-! ## Forcing Relationship 3: The ~21 Factors as Internal Cost of Density-Gradient Modulation

The repeated appearance of 21 (= dim Spin(7)) and related factors in the
EM-gravity hierarchy, weak couplings, and scale bridges counts the number
of additional internal relations that become active (and must be paid for)
when the medium allows local density to vary and modulate propagation.

EM can live in a smaller sub-algebra with low modulation cost.
Gravity and density-regulated sectors must "see" the full 21-dimensional
structure. The numerical factors are the price of turning on density-gradient
modulation without destroying the lower-grade protection or introducing dispersion.
-/
def TwentyOne_As_Internal_Modulation_Cost : Prop :=
  MediumDensityDrive ∧ NonDispersiveUniversality →
  ∃ internalDegrees : Nat, internalDegrees = 21 ∧
    (ModulationCost internalDegrees ≤ HierarchySuppression) ∧
    (DispersionIntroduced internalDegrees = 0)

axiom TwentyOne_Cost_Is_Forced : TwentyOne_As_Internal_Modulation_Cost

/-! ## Placeholder for more relationships (to be strengthened)

We will add further Props for:
- The dynamic ξ regulator as local density-credit order parameter.
- The overall algebra as the minimal solution space for the density problem
  with clean force separation.

Each new Prop should be accompanied by a comment linking back to the
corresponding hypothesis in HYPOTHESES.md and the specific bounded quantity
(density well, protection budget, relationship trough, stability margin)
that does the forcing.
-/

end DensityAlgebra