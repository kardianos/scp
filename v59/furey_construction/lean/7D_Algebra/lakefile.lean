import Lake
open Lake DSL

package Furey7DAlgebra where
  -- Dedicated package for the 7D Cl(7)_even / oct GA + 8D spinor action.
  -- Can be built standalone or later merged into parent furey/lean/ package.
  -- Mathlib intentionally omitted for this integration phase (all code is
  -- self-contained List-based; re-add when we need Matrix/Finset tactics).

lean_lib Furey7D where
  -- Explicit octonion Cl(7) algebra (genuine, post 2026-05-24 table fix) + grade-structure
  -- theorems + the bug-immune Shulga ratio. SevenDAlgebra is the source of truth.
  -- (The former StabilityFromAlgebra + PhaseC_Certificates crossover/stability layer was
  --  REMOVED 2026-05-24 — it modelled a Cl(3,0) phenomenon the octonion algebra does not
  --  exhibit; see notes/2026-05-24-LivingCandidateCrossover-ReEvaluation.md.)
  srcDir := "."
  -- NOTE: `Z2Z2Forcing` imports `CliffordBladeGrade` (Mathlib, in the parent package) for the
  -- rigorous grade bridge, so it builds in the PARENT package, not here.
  roots := #[`SevenDAlgebra, `PhaseB_Theorems, `ShulgaParameters, `LeptonGradeForcing,
             `LeptonComplexStructure]

lean_exe shulga where
  root := `export_shulga