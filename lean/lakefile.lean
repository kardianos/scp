import Lake
open Lake DSL

package ScpLean where
  version := v!"0.1.0"
  leanOptions := #[⟨`autoImplicit, false⟩]

@[default_target]
lean_lib ScpLib where
  srcDir := "."

@[default_target]
lean_lib V44 where
  srcDir := "."

@[default_target]
lean_lib V50C3 where
  srcDir := "."

@[default_target]
lean_lib V50C4 where
  srcDir := "."

@[default_target]
lean_lib Polariton where
  srcDir := "."

@[default_target]
lean_lib SinglePass where
  srcDir := "."
