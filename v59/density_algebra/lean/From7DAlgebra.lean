/-
  DensityAlgebra.From7DAlgebra

  Thin compatibility / re-export layer for the explicit 7D algebra realization.

  The authoritative implementation of the 8-dimensional algebra (Fano multiplication
  table, Fock |Ω_N⟩ labeling, gamma generators, L-grade and F-grade operator
  matrices, sector diagnostics) lives in:

      v59/furey_construction/lean/7D_Algebra/SevenDAlgebra.lean

  The canonical stability bridge that builds on those matrices lives in:

      v59/furey_construction/lean/7D_Algebra/StabilityFromAlgebra.lean

  This module exists so that code in `density_algebra/lean/` can have a single,
  well-documented place to reference the real algebraic foundation during the
  integration period.

  When the Lake dependency graph is cleaned up, this file can become a proper
  re-export or be removed in favor of direct imports.
-/ 

namespace DensityAlgebra

/-! ## Source of Truth Pointers

All new development of the real 7D algebra and its use for stability bounds,
Z₂×Z₂ forcing, and Option D compositeness should happen in the 7D_Algebra folder.

Key modules to import / study:
- `SevenDAlgebra`          — core algebra, Fock labeling, explicit 8×8 matrices
- `StabilityFromAlgebra`   — leakage demo, numeric Hessian, public stability helpers

See:
- `v59/furey_construction/lean/7D_Algebra/INTEGRATION_PLAN.md`
- `v59/furey_construction/lean/7D_Algebra/SevenDAlgebra.lean`
- `v59/furey_construction/lean/7D_Algebra/StabilityFromAlgebra.lean`

For the time being, the older schematic modules (`OctonionAlgebra.lean`,
`StabilityBounds.lean`) remain for compatibility and historical reference.
The real work is being migrated to the modules above.
-/

end DensityAlgebra