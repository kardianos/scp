/-
  v59/furey_construction/lean/CliffordBladeGrade.lean

  The **grade structure** behind the Z₂×Z₂ forcing, proved universally and
  non-enumeratively (no "we checked a few generators").

  Motivation.  `7D_Algebra/Z2Z2Forcing.lean` established the *representation*-level fact
  "the product of the two specific disjoint bivectors γ₀γ₁, γ₂γ₃ is an F-grade 4-form" by
  `decide`.  The structural content is a grade law that holds for **every** Clifford
  algebra and **all** disjoint bivectors, derived from the index-set combinatorics — not
  from checking matrices.

  Model.  In any Clifford algebra with anticommuting generators (`eᵢeⱼ = -eⱼeᵢ`, `i≠j`)
  and `eᵢ² = (scalar)`, a *basis blade* is the ordered product `e_S = ∏_{i∈S} eᵢ` indexed
  by a finite set `S`.  The geometric product of two basis blades is

        e_S · e_T  =  (±1) · e_{S △ T},

  because the shared indices `S ∩ T` square to scalars and cancel, leaving the symmetric
  difference `S △ T` (the sign is the parity of the reordering, irrelevant to the grade).
  Hence the **grade** (number of generators) of a product is `|S △ T|`, a function of the
  index sets alone.  We model a blade by its index set and `bladeMul` by `symmDiff`, and
  prove the grade law structurally.

  This is representation-independent: the explicit 8×8 Fano/octonion matrices of
  `SevenDAlgebra` are *one* realization (they satisfy exactly these generating relations),
  so the universal grade law below is the structural reason `L·L` (disjoint bivectors)
  lands in F — superseding the enumerative witness.
-/
import Mathlib

namespace SCPv59.CliffordBladeGrade

open Finset

/-- A Clifford basis blade over `n` generators, identified with its index set. -/
abbrev Blade (n : ℕ) := Finset (Fin n)

/-- The grade of a blade = number of generators it contains. -/
def grade {n} (S : Blade n) : ℕ := S.card

/-- The basis-blade (geometric) product, at the index-set level: `e_S · e_T = ±e_{S △ T}`,
    so the product blade's index set is the symmetric difference. -/
def bladeMul {n} (S T : Blade n) : Blade n := symmDiff S T

/-! ## The grade law -/

/-- **The grade multiplication rule.**  `grade(e_S · e_T) = |S| + |T| − 2|S∩T|`.
    Derived from `symmDiff = (S∪T) \ (S∩T)` and the cardinality identities — no
    enumeration over generators. -/
theorem grade_bladeMul {n} (S T : Blade n) :
    grade (bladeMul S T) = grade S + grade T - 2 * (S ∩ T).card := by
  have hsub : S ∩ T ⊆ S ∪ T := (Finset.inter_subset_left).trans Finset.subset_union_left
  have hu : (S ∪ T).card + (S ∩ T).card = S.card + T.card := Finset.card_union_add_card_inter S T
  unfold grade bladeMul
  rw [symmDiff_eq_sup_sdiff_inf, Finset.sup_eq_union, Finset.inf_eq_inter, Finset.card_sdiff,
      Finset.inter_eq_left.mpr hsub]
  omega

/-- For disjoint blades, grades simply add: `grade(e_S · e_T) = |S| + |T|`. -/
theorem grade_bladeMul_disjoint {n} (S T : Blade n) (h : Disjoint S T) :
    grade (bladeMul S T) = grade S + grade T := by
  unfold grade bladeMul
  rw [h.symmDiff_eq_sup]; exact Finset.card_union_of_disjoint h

/-! ## The Furey grades and the universal forcing -/

/-- F-grade content: the G₂-form grade Λ⁴. -/
def isF {n} (S : Blade n) : Prop := grade S = 4

/-- L-grade content: the "Lie-algebra" grades Λ²⊕Λ⁶ (`= dim Spin(8)` total). -/
def isL {n} (S : Blade n) : Prop := grade S = 2 ∨ grade S = 6

/-- **THE UNIVERSAL STRUCTURAL THEOREM.**  For *any* `n` and *any* two **disjoint
    bivectors** `S, T` (grade-2 L-blades), their product is a grade-4 **F-blade**.

    This is not "we checked two": it is the grade law `2 + 2 = 4` for disjoint index sets,
    holding for every disjoint pair of bivectors in every Clifford algebra. -/
theorem disjoint_bivectors_mul_isF {n} (S T : Blade n)
    (hS : grade S = 2) (hT : grade T = 2) (hd : Disjoint S T) :
    isF (bladeMul S T) := by
  unfold isF; rw [grade_bladeMul_disjoint S T hd, hS, hT]

/-- **The full bivector-product dichotomy.**  Two bivectors multiply to grade
    `4 − 2|S∩T|`: F-grade when disjoint (`|S∩T|=0`), an L-bivector when they share one
    index (`|S∩T|=1`), and a scalar when equal (`|S∩T|=2`).  The disjoint case is the only
    one that leaves L — and it is exactly the F-grade. -/
theorem bivectors_mul_grade {n} (S T : Blade n) (hS : grade S = 2) (hT : grade T = 2) :
    grade (bladeMul S T) = 4 - 2 * (S ∩ T).card := by
  rw [grade_bladeMul, hS, hT]

/-- **L is not closed under multiplication (universal form).**  The product of *any* two
    disjoint bivectors is F-grade and is **not** L-grade — so the bilinear product map
    sends every disjoint bivector pair out of L into F.  (This replaces the enumerative
    `Z2Z2Forcing.L_not_closed_reaches_F`, which exhibited a single such pair.) -/
theorem L_not_closed_on_disjoint_bivectors {n} (S T : Blade n)
    (hS : grade S = 2) (hT : grade T = 2) (hd : Disjoint S T) :
    isF (bladeMul S T) ∧ ¬ isL (bladeMul S T) := by
  refine ⟨disjoint_bivectors_mul_isF S T hS hT hd, ?_⟩
  unfold isL
  rw [grade_bladeMul_disjoint S T hd, hS, hT]
  omega

/-- Non-vacuity over the v59 algebra `Cl(7)`: disjoint bivectors exist (e.g. `{0,1}`,
    `{2,3}`), so the universal statement above is non-empty. -/
theorem disjoint_bivectors_exist :
    ∃ S T : Blade 7, grade S = 2 ∧ grade T = 2 ∧ Disjoint S T :=
  ⟨{0, 1}, {2, 3}, by decide, by decide, by decide⟩

/-- **The compositeness forcing, at the grade level.**  Combining the above: there exist
    disjoint bivectors (so L is genuinely not closed), and the product of *any* disjoint
    bivector pair lands in F (never back in L).  A state built by multiplying L-grade
    excitations — the N=2 u-quark in the Fock/Witt construction — is therefore pushed into
    the F-grade; with `Z2Z2Forcing.F_color_cartan_rank_ge_two` (only F carries the color
    Cartan) this forces a color-charged composite to occupy L⊕F rather than pure L. -/
theorem composite_forced_out_of_L :
    (∃ S T : Blade 7, grade S = 2 ∧ grade T = 2 ∧ Disjoint S T)
    ∧ (∀ S T : Blade 7, grade S = 2 → grade T = 2 → Disjoint S T →
          isF (bladeMul S T) ∧ ¬ isL (bladeMul S T)) :=
  ⟨disjoint_bivectors_exist, fun S T hS hT hd => L_not_closed_on_disjoint_bivectors S T hS hT hd⟩

end SCPv59.CliffordBladeGrade
