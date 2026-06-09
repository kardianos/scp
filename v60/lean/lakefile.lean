import Lake
open Lake DSL

package «scp_v60» where

/-
v60 Lean modules for the G9 (induced-metric / soldering) attack.

BUILD: Mathlib v4.29.0 is already compiled in the v59 tree, so the cheapest and
verified way to check these modules is against that environment (no fresh Mathlib
build needed):

  cd ../../v59/furey_construction/lean
  lake env lean ../../../v60/lean/G9Soldering.lean
  lake env lean ../../../v60/lean/G9InducedMetric.lean
  lake env lean ../../../v60/lean/G9ToyHelicity.lean
  lake env lean ../../../v60/lean/G9Unification.lean

All four compile clean (no errors, no warnings).  The headline theorems depend
only on the standard trio [propext, Classical.choice, Quot.sound] — no `sorryAx`
(several, e.g. no_internal_lorentz / graviton_dof_covariant, use no axioms at all).

A standalone `lake build` here would re-fetch and rebuild Mathlib from scratch;
prefer the `lake env lean` route above unless a fresh build is intended.
-/

@[default_target]
lean_lib G9Soldering where
  -- soldering = Minkowski sum of helicity charges; no-go vs resolution; charpoly X²+4

lean_lib G9InducedMetric where
  -- helicity classification (valid restatement of v59) + the single open claim

lean_lib G9ToyHelicity where
  -- corrected, non-circular nugget: TT generator charpoly = X² + 4

lean_lib G9Unification where
  -- Cl(3,1)⊗Cl(7)_even: general tensor-factor commutation + dim identities (closes C3)

lean_lib G9Plebanski where
  -- self-dual 't Hooft algebra (orthonormal triple) + 2-form dim counts + DOF=2

lean_lib G9ObeToPlebanski where
  -- NEGATIVE: Plebański not derivable from OBE (DOF/field-count obstruction)

lean_lib RankTension where
  -- G1: 784=28², the 6/784 deflation, readings differ, subalgebra obstruction

lean_lib ParentAction where
  -- GEN1: first-order parent dissolves the 09 obstruction (common-ancestor DOF lattice)

lean_lib CovariantFirstOrder where
  -- GEN2: covariant first-order Palatini -> 2 ghost-free TT DOF; B∧F = grade projection

lean_lib MatterSector where
  -- GEN3: potential vacuum derives Koide Q=2/3 (e1²=6e2); phase=Goldstone; 2nd moment=ρ_grav

lean_lib MatterGravityCoupling where
  -- GEN4: universal coupling (one f_g), EP exact (mass ratio 1), charge=9Qa², 1/r & 1/r² law

lean_lib SpectrumStability where
  -- GEN5: joint spectrum decoupled; no_tachyon (PSD via positivity); 5 modes, 0 ghost/tachyon

lean_lib SelectionRule where
  -- GEN6: selection rule (D_u=D_e+D_d); universal Koide deviation (1-Q)D=28/3; two-object G1

lean_lib DynamicsSpectrum where
  -- GEN7: nonlinear EL dynamics -> dispersion ω²=k²+m²; massless Goldstone (speed 1); ω²=m² at rest

require mathlib from git
  "https://github.com/leanprover-community/mathlib4.git" @ "v4.29.0"
