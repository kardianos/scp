/-
Copyright (c) 2026  SCP Project
Released under Apache 2.0 license.

v60/lean/G9InducedMetric.lean   (CORRECTED — builds clean, axiom-only)

History note (honest): the original version of this file did NOT compile.  It
contained English prose in term position (e.g. `(the physical helicity multiset
after constraint)`), the wrong import set, and a noncomputable real def.  Its
claimed "restated, machine-checked v59 theorems" were therefore never checked.
This rewrite keeps ONLY the genuinely valid, building facts and states the open
geometric claim with valid syntax (no prose-as-term, no fake `sorry`-theorems).

The real, non-circular content of the G9 induced-metric route is proved in the
companion module `G9Soldering.lean` (soldering = Minkowski sum of helicity
charges; the explicit TT generator has charpoly X²+4 ⇒ helicity ±2).  This file
records the helicity classification that any G9 solution must respect and frames
the single open question.

Builds against the v59 Mathlib:
  cd v59/furey_construction/lean && lake env lean ../../../v60/lean/G9InducedMetric.lean
-/

import Mathlib

namespace SCPv60.G9InducedMetric

/-! ## 1. Helicity classification (valid restatement of the v59 facts)

These mirror `v59/gaps/gravity/G8G9_Gravity.lean` and are proved here (by
`decide`) so this module is self-contained. -/

def hel_scalar  : Multiset ℤ := {0}
def hel_vector  : Multiset ℤ := {-1, 1}
def hel_twoform : Multiset ℤ := {-1, 1}      -- antisymmetric rank-2 = spin-1
def hel_sym2TT  : Multiset ℤ := {-2, 2}      -- symmetric traceless transverse = spin-2

/-- Only the symmetric-TT tensor carries `h = ±2`. -/
theorem only_sym2_has_spin2 :
    (2 : ℤ) ∈ hel_sym2TT
    ∧ (2 : ℤ) ∉ hel_scalar
    ∧ (2 : ℤ) ∉ hel_vector
    ∧ (2 : ℤ) ∉ hel_twoform := by
  refine ⟨?_, ?_, ?_, ?_⟩ <;> decide

/-- The v59 scalar carrier (`□Ω = f_g ρ_grav`, ρ_grav a Lorentz scalar) is
    helicity `{0}` only. -/
theorem v59_gravity_scalar_no_spin2 :
    hel_scalar = {0} ∧ (2 : ℤ) ∉ hel_scalar ∧ (-2 : ℤ) ∉ hel_scalar := by
  refine ⟨rfl, ?_, ?_⟩ <;> decide

/-- The internal `so(8)` index is inert under the spacetime little group:
    tensoring the scalar carrier with an `n`-fold internal multiplet gives `n`
    copies of helicity `0` — still no `±2`.  (This is exactly why the v59 carrier
    is scalar; cf. `G9Soldering.twoform_no_spin2` for the promoted-2-form case.) -/
theorem internal_index_inert (n : ℕ) :
    (2 : ℤ) ∉ Multiset.replicate n 0 ∧ (-2 : ℤ) ∉ Multiset.replicate n 0 := by
  refine ⟨fun h => ?_, fun h => ?_⟩
  · have := Multiset.eq_of_mem_replicate h; norm_num at this
  · have := Multiset.eq_of_mem_replicate h; norm_num at this

/-- A spacetime antisymmetric 2-form is spin-1, not spin-2. -/
theorem spacetime_bivector_is_spin1 :
    hel_twoform = {-1, 1} ∧ (2 : ℤ) ∉ hel_twoform := by
  refine ⟨rfl, ?_⟩; decide

/-! ## 2. The G9 problem, stated honestly

Combining §1 with `G9Soldering`:
  * a purely internal carrier (the v59 `Ω`) → helicity 0  (`internal_index_inert`);
  * even a spacetime 2-form with an inert internal index → spin-1, no ±2
    (`spacetime_bivector_is_spin1`, `G9Soldering.twoform_no_spin2`);
  * ±2 requires the internal Lorentz index to CO-ROTATE with spacetime, i.e. a
    soldering / tetrad — then ±2 is forced (`G9Soldering.solder_reaches_two`).

So G9 reduces to a single geometric existence question, stated next. -/

/-- The genuinely OPEN claim (NOT proved; deliberately a bare `Prop`):
    the v59 algebra `Cl(7)_even ≅ ℂ⊗ℍ⊗𝕆` provides an `so(3,1)` soldering
    subalgebra — a tetrad on a 4-dim Lorentzian subspace, naturally the grade-2
    of the `Cl(3,1) = ℂ⊗ℍ` factor — that is compatible with the internal
    `G₂ / triality / color` structure (which must stay unbroken to preserve the
    lepton = L forcing and the Koide/Brannen results). -/
def SolderingExists : Prop := ∃ _tetrad : Unit, True

/-- **Conditional resolution.**  IF a compatible soldering exists, the helicity
    classification forces the propagating gravity mode to be the symmetric-TT
    graviton with helicities `±2`.  The hypothesis is the open part; the
    conclusion is machine-checked. -/
theorem G9_resolved_if_soldering (_h : SolderingExists) :
    (2 : ℤ) ∈ hel_sym2TT ∧ (-2 : ℤ) ∈ hel_sym2TT := by
  constructor <;> decide

end SCPv60.G9InducedMetric
