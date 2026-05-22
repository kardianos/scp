import Lake
open Lake DSL

package «unified_multivector_force» where
  version := v!"0.1.0"
  leanOptions := #[⟨`autoImplicit, false⟩]

@[default_target]
lean_lib UnifiedMultivector where
  srcDir := "."
