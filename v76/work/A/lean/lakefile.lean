import Lake
open Lake DSL

package «v76A» where
  version := v!"0.1.0"
  leanOptions := #[⟨`autoImplicit, false⟩]

@[default_target]
lean_lib LocalityBudget where
  roots := #[`LocalityBudget]
  srcDir := "."
