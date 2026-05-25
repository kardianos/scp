import Lake
open Lake DSL

package density_algebra where
  -- Lightweight package for stating forcing relationships as Props
  -- No heavy mathlib dependency yet; we stay at the level of axioms and Prop

@[default_target]
lean_lib DensityAlgebra where
  roots := #[`DensityForcing, `StabilityBounds, `OctonionAlgebra, `From7DAlgebra]