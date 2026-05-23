import Lake
open Lake DSL

require "leanprover-community" / "mathlib" @ git "v4.29.0"

package SCPv59 where
  -- nothing fancy; minimal package

lean_lib «LieDimensions» where
  -- the abstract dim identities

lean_lib «Octonions» where
  -- the Fano plane structure

lean_lib «KoideAndBrannen» where
  -- the structural identifications

lean_lib «BrannenKernel» where
  -- the Brannen mass-operator eigenvalue identity (Q = (1+2t^2)/3)

lean_lib «CyclicShift» where
  -- ω = exp(2πi/3), the Z_3 structure underlying the lepton kernel

lean_lib «SpinDimension» where
  -- dim Spin(n) = n(n-1)/2 ⇒ dim Spin(7) = 21 (structural derivation)

lean_lib «SilentDirection» where
  -- silent-direction theorem: conjugation by unit q preserves (Re ξ, normSq ξ)

lean_lib «Predictions» where
  -- consolidated v59 prediction table: machine-checked structural identities

lean_lib «EmbeddingIndex» where
  -- Killing-form embedding index for so(n) ⊂ so(N)

lean_lib «AxiomCheck» where
  -- prints axioms used by the headline theorems

lean_lib «ScaleBridge» where
  -- EW-sector closure: v_Higgs = 28²·a², sin²θ_W = 2/9, α(M_Z) = 25/(324π²)
