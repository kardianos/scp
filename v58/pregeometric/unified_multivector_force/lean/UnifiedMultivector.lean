/-
  UnifiedMultivector — Root module for the Lean 4 formalization of the
  unified multivector force law experiment (v58 pregeometric track).

  This package encodes:
  - A minimal fragment of 3D geometric algebra (Cl(3,0)) sufficient for
    grade projections and geometric products.
  - Speculative candidate unified equations (from BACKGROUND_AND_SPECULATIVE_EQUATIONS.md).
  - Theorems stating implications to Newtonian gravity (with density dependence)
    and static Maxwell / Coulomb limits, under explicitly documented assumptions.

  The style follows the project's existing conventions in root `lean/`:
  - Axiomatic reals (R) for algebraic reasoning without requiring Mathlib initially.
  - Clear namespace, structured comments, and `sorry` for unfinished steps.
  - Emphasis on what assumptions are required for proofs to succeed (these feed
    the Python discovery track).

  Future: Mathlib4 + Mathlib.Algebra.CliffordAlgebra can be added for richer
  abstract treatment once the concrete implications are clear.
-/

import UnifiedMultivector.Multivector
import UnifiedMultivector.Candidates
import UnifiedMultivector.KnownPhysics
import UnifiedMultivector.NewtonianLimit
import UnifiedMultivector.MaxwellLimit
import UnifiedMultivector.Assumptions
import UnifiedMultivector.Model   -- concrete Fin-8 realization (Cycle 2 start)
