/-
Copyright (c) 2026  SCP Project
Released under Apache 2.0 license.

v60/lean/G9ToyHelicity.lean   (CORRECTED & DOWNSCOPED — builds clean)

History note (honest): the original version of this file claimed a "full
derivation, no `sorry` on the spectrum, machine-checked".  It was NOT checked:
it imported a non-existent module path (`Mathlib.LinearAlgebra.Matrix.Charpoly`)
and tried to close `Multiset ℂ` goals with `decide` (complex equality is not
decidable) and to use `Matrix.eigenvalues` (a `IsHermitian`-only API) on a
non-Hermitian matrix.  It never compiled.

More importantly, its *physics* was circular: it built `Jz_TT ⊗ Id` from the
transverse-traceless (spin-2) generator inserted BY HAND, so finding helicity
±2 was tautological — the "constraint" did no work.  See the audit and
`v60/gravity_recast/04_soldering_helicity_honest.py`.

This rewrite keeps only the *correct, non-circular* nugget: the explicit
transverse-traceless rotation generator `!![0,-2; 2,0]` has characteristic
polynomial `X² + 4`, hence eigenvalues `±2i` (helicity ±2), DERIVED from the
matrix.  The real, non-circular G9 content (why ±2 needs soldering, and that the
soldered object reaches it) is in the companion `G9Soldering.lean`.

Builds against the v59 Mathlib:
  cd v59/furey_construction/lean && lake env lean ../../../v60/lean/G9ToyHelicity.lean
-/

import Mathlib

namespace SCPv60.G9ToyHelicity

/-- The transverse-traceless rotation generator on the graviton polarizations
    `(h_+, h_×)`, scaled so that `|eigenvalue| = 2` (spin-2 helicity). -/
def JzTT : Matrix (Fin 2) (Fin 2) ℂ := !![0, -2; 2, 0]

/-- **Charpoly `= X² + 4`** (trace 0, det 4), derived from the entries.
    Its roots are `±2i`, i.e. the helicities `±2`.  This is the corrected,
    non-fragile version of what the original file tried to do with a 4×4 matrix
    and `decide` over `Multiset ℂ`. -/
theorem JzTT_charpoly : JzTT.charpoly = Polynomial.X ^ 2 + Polynomial.C 4 := by
  rw [Matrix.charpoly_fin_two]
  simp only [JzTT, Matrix.trace_fin_two, Matrix.det_fin_two]
  ring_nf
  norm_num

/-- Honest caveat, recorded as a definition: the eigenvalue `±2` here is a
    property of a generator that already acts on a spacetime spin-2 object.  It
    becomes a *derivation* of LIGO content only once one explains why that
    object exists, i.e. soldering — see `G9Soldering.solder_reaches_two`. -/
def caveat_not_circular_only_with_soldering : Prop := True

end SCPv60.G9ToyHelicity
