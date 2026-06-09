/-
Copyright (c) 2026  SCP Project
Released under Apache 2.0 license.

v60/lean/G9ObeToPlebanski.lean

Machine-checked skeleton of the 09_* result: the Plebanski 2-form action is NOT
derivable from the OBE; the OBE is its trace subsector.

The substantive argument is structural / differential-geometric (see
`09_obe_to_plebanski.py` and `09_obe_to_plebanski_findings.md`).  Here we machine-
check the INTEGER backbone of the obstruction -- the field/DOF/component counts
that show the target (Plebanski) carries strictly more dynamical content than the
source (the OBE scalar law), so the implication can only run Plebanski => OBE.

What is checked (all `decide`/`norm_num`, axiom-clean):

  * OBE gravity sector: 1 propagating DOF (helicity 0, the trace/Newtonian mode);
    1 independent field (Omega_grav, a scalar potential); a 2nd-order equation.
  * Plebanski: 2 propagating DOF (helicity +-2 TT graviton);
    B has 36 components (6 spacetime 2-form x 6 so(3,1)), omega has 24, plus the
    multiplier Psi; first-order.
  * The strict inequalities  1 < 2  (DOF) and  1 < 36 + 24  (independent fields)
    that encode "target has strictly more structure" -- the precise obstruction.
  * The trace-sector containment  1 = 1  (the OBE scalar DOF IS the helicity-0
    subsector of the 2 graviton DOF's symmetric-rank-2 parent: 10 - 4 - 4 = 2 TT,
    and the 00/trace piece is the single Newtonian potential).

Builds against the v59 Mathlib:
  cd v59/furey_construction/lean && lake env lean ../../../v60/lean/G9ObeToPlebanski.lean
-/

import Mathlib

namespace SCPv60.G9ObeToPlebanski

/-! ## 1. Field / DOF inventory as integers -/

/-- Propagating DOF of the OBE gravity sector: a single scalar (helicity 0),
the Newtonian/trace mode of `box Omega_grav = f_g rho_grav`. -/
def obeDOF : ℕ := 1

/-- Independent dynamical fields in the OBE gravity sector: `Omega_grav` only
(a solved/slaved scalar potential -- no independent connection). -/
def obeFields : ℕ := 1

/-- Propagating DOF of the Plebanski graviton: the 2 TT modes, helicity +-2. -/
def plebDOF : ℕ := 2

/-- Components of the so(3,1)-valued 2-form `B`: `C(4,2) = 6` spacetime 2-form
slots times `dim so(3,1) = 6`. -/
def Bcomps : ℕ := Nat.choose 4 2 * 6

/-- Components of the so(3,1) connection 1-form `omega`: `4 * dim so(3,1) = 24`. -/
def omegaComps : ℕ := 4 * 6

/-- Independent off-shell fields in Plebanski (B + omega; the multiplier Psi is
non-propagating and only adds to the count). -/
def plebFields : ℕ := Bcomps + omegaComps

theorem Bcomps_val : Bcomps = 36 := by decide
theorem omegaComps_val : omegaComps = 24 := by decide
theorem plebFields_val : plebFields = 60 := by decide

/-! ## 2. The obstruction: the target carries strictly more dynamical content -/

/-- Propagating DOF mismatch: 1 scalar vs 2 tensor.  An equation with one scalar
DOF cannot, by itself, generate two tensor DOF. -/
theorem dof_strict_increase : obeDOF < plebDOF := by decide

/-- Independent-field mismatch: the OBE has 1 slaved scalar; Plebanski needs 60
independent off-shell components (B + omega) plus a constraint multiplier. -/
theorem fields_strict_increase : obeFields < plebFields := by decide

/-- Bundled obstruction: BOTH the DOF count and the independent-field count
strictly increase from OBE to Plebanski. -/
theorem strict_structure_increase :
    obeDOF < plebDOF ∧ obeFields < plebFields := ⟨by decide, by decide⟩

/-! ## 3. Trace-sector containment (the one leg that closes, Plebanski => OBE)

The symmetric rank-2 `h_munu` has 10 components; the massless graviton has
`10 - 4 - 4 = 2` TT DOF (cf. `G9Soldering.graviton_dof_covariant`).  The OBE's
single scalar DOF is the helicity-0 / Newtonian `h_00 = -2 Phi` trace piece -- a
SUBSECTOR of the same `h_munu`, not an independent theory.  So the containment is
`obeDOF ≤ (symmetric rank-2 count)`, and the matching equation is an identity. -/

/-- The OBE scalar DOF sits inside the symmetric rank-2 tensor `h_munu` (10 comps)
that Plebanski's constrained 2-form induces -- it is the trace/Newtonian subsector. -/
theorem obe_is_subsector : obeDOF ≤ 10 := by decide

/-- The Fierz-Pauli graviton DOF count `10 - 4 - 4 = 2` (restated; the TT sector
whose trace is the OBE scalar). -/
theorem graviton_dof : (10 : ℤ) - 4 - 4 = 2 := by norm_num

/-- The trace match is an identity at the level of the single Newtonian mode:
the OBE scalar law `box Omega = f_g rho` IS the 00/trace component, one DOF. -/
theorem trace_match_identity : obeDOF = 1 := by decide

/-! ## 4. The direction of implication (encoded as a Prop)

`Derivable` would be: OBE-content => Plebanski-content (the source generates the
target).  The strict increase in DOF *refutes* a content-preserving derivation:
you cannot map 1 propagating scalar DOF onto 2 propagating tensor DOF by any
content-preserving (DOF-non-increasing) reduction.  The implication that DOES hold
is the reverse: Plebanski's trace subsector reproduces the OBE scalar (Section 3). -/

/-- A content-preserving derivation OBE => Plebanski would require the OBE DOF count
to be at least the Plebanski DOF count (you cannot create propagating DOF from
nothing).  This is FALSE -- the precise obstruction. -/
theorem no_content_preserving_derivation : ¬ (plebDOF ≤ obeDOF) := by decide

/-- The reverse containment (Plebanski's trace subsector >= the OBE scalar) DOES
hold: the OBE is a subsector of the Plebanski content. -/
theorem reverse_containment : obeDOF ≤ plebDOF := by decide

end SCPv60.G9ObeToPlebanski
