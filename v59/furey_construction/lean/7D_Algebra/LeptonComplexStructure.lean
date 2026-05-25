import SevenDAlgebra

namespace SCPv59.Furey7D.LeptonComplexStructure

open SCPv59.Furey7D

/-!
# The lepton complex structure lives in Λ² — this forces lepton = L

This module closes the open half of the Z₂×Z₂ forcing identified in
`LeptonGradeForcing.lean`.  That module proved lepton = L is *not* forced by
mass-channel **availability** (both grades furnish a lepton coupling, and F is if
anything richer), and isolated the precise missing ingredient:

> which grade carries the ℂ-imaginary unit `J` of Furey's minimal left ideal — if
> `J` lives in the Λ² bivector slice (as `brannen_kernel.py` embeds the lepton ξ),
> then lepton = L is forced.

Here that ingredient is supplied and the forcing is completed, by an exact fact
about the genuine (post sign-fix) Cl(7) generators.

## The mechanism

The Brannen lepton mass kernel `M = a(I + ξ S + ξ̄ S²)` (`KernelEigenvalues.lean`)
carries a **phase** `ξ = e^{iφ}`, `φ = 2/9`.  The nontrivial phase (`φ ≠ 0`, i.e.
`Q = 2/3 ≠ 1`) is what produces the Koide ratio; it requires a genuine **complex
structure** `J` (`J² = −1`) — the `i` in `ξ = cos φ + i sin φ`.

`brannen_kernel.py` (lines 9–16, 123–137) embeds the lepton's quaternion ℍ-slice
`⟨1, i, j, k⟩` as **three bivectors** `(e₀₁, e₀₂, e₁₂) ∈ Λ²`, so the complex unit
is `i = e₀₁`, a Λ²-blade.  In our Cl(7) realization `e₀₁ = γ₀γ₁`.

## Why it MUST be L, not F (the forcing)

A complex structure needs `J² = −I`.  For a simple even k-blade `B` (a product of
`k` distinct anticommuting generators with `γ² = −I`) one has

      B² = (−1)^{k(k+1)/2} · I,

so within the even subalgebra `Cl(7)_even = Λ⁰ ⊕ Λ² ⊕ Λ⁴ ⊕ Λ⁶`:

  * `Λ⁰` (k=0): `B² = +I`     — real (the identity)
  * `Λ²` (k=2): `B² = −I`     — **complex structure**   ⊂ L
  * `Λ⁴` (k=4): `B² = +I`     — real structure (involution)   = F
  * `Λ⁶` (k=6): `B² = −I`     — **complex structure**   ⊂ L

Hence **`L = Λ²⊕Λ⁶` is exactly the set of even grades whose blades are complex
structures, while `F = Λ⁴` furnishes only real structures.**  The lepton's phase
`ξ = e^{iφ}` therefore *cannot* be carried by an F-blade; its complex structure is
necessarily L-grade.  **lepton = L is forced.**

This is the complex-structure face of the dichotomy `LeptonGradeForcing` proved on
the singlet pair {0,7}: `J² = −I` ⇔ the antisymmetric ε/rotation block (L), while
`J² = +I` ⇔ the symmetric metric/involution block (F).

All statements below are `decide`-checked on the genuine Cl(7) matrices
(`propext`/`Quot.sound` only — no `native_decide`, no `sorry`, no extra axioms).
The square law is verified **exhaustively** over all 28 L-blades and all 35
F-blades of Cl(7), not on a sample.
-/

/-- `−I₈`, the matrix a complex structure squares to. -/
def negId8 : Mat8 := matScale (-1) id8

/-- The 7 six-form (`Λ⁶`) blades: each omits one generator (`γ^{(7)}` minus `γ_d`). -/
def hexads7 : List (List Nat) :=
  (List.range 7).map (fun d => (List.range 7).filter (fun k => k != d))

/-- All 7 `Λ⁶` operator matrices (the other half of `L = Λ²⊕Λ⁶`). -/
def all_L_sixforms : List Mat8 :=
  hexads7.map (fun ks => ks.foldl (fun acc k => matMul acc (gamma k)) id8)

/-- The full L-grade blade list `Λ² ⊕ Λ⁶` (21 + 7 = 28 generators). -/
def all_L_grade : List Mat8 := all_L_bivectors ++ all_L_sixforms

/-- A matrix `M` is a **complex structure** iff `M² = −I`.  (`abbrev` so `decide`
    sees through to the underlying decidable matrix equality.) -/
abbrev isComplexStructure (M : Mat8) : Prop := matMul M M = negId8

/-- A matrix `M` is a **real structure** (involution) iff `M² = +I`. -/
abbrev isRealStructure (M : Mat8) : Prop := matMul M M = id8

/-! ## Part 1 — the lepton ℍ-slice and its complex unit `i = e₀₁ = γ₀γ₁` -/

/-- The lepton quaternion-slice imaginary units, the three Λ²-bivectors that
    `brannen_kernel.py` uses for `(i, j, k)`: `(e₀₁, e₀₂, e₁₂)`. -/
def i_lep : Mat8 := matMul (gamma 0) (gamma 1)
def j_lep : Mat8 := matMul (gamma 0) (gamma 2)
def k_lep : Mat8 := matMul (gamma 1) (gamma 2)

/-- **The lepton complex unit `i = e₀₁` is a genuine complex structure.**
    `(γ₀γ₁)² = −I`: this is the `i` of the Brannen phase `ξ = e^{iφ}`. -/
theorem i_lep_isComplexStructure : isComplexStructure i_lep := by decide

/-- **The lepton ℍ-slice `(e₀₁, e₀₂, e₁₂)` is a genuine quaternion algebra ℍ.**
    Hamilton's relations: `i²=j²=k²=−I`, `ij=k`, `jk=i`, `ki=j`, and `ijk=−I`.
    So the lepton mass operator is built over `ℍ ⊂ {Λ⁰} ⊕ Λ²` — purely L-grade
    (plus the scalar). -/
theorem lepton_slice_is_quaternion :
    matMul i_lep i_lep = negId8 ∧ matMul j_lep j_lep = negId8 ∧
    matMul k_lep k_lep = negId8 ∧
    matMul i_lep j_lep = k_lep ∧ matMul j_lep k_lep = i_lep ∧
    matMul k_lep i_lep = j_lep ∧
    matMul (matMul i_lep j_lep) k_lep = negId8 := by decide

/-- **The lepton complex unit is an L-grade (Λ²) blade**, a member of the 21
    bivectors — not an F-grade element. -/
theorem i_lep_is_L_bivector : i_lep ∈ all_L_bivectors := by decide

/-! ## Part 2 — the exhaustive grade ⇒ square-sign law on Cl(7)

`L = Λ²⊕Λ⁶` is *exactly* the complex-structure grades; `F = Λ⁴` is real structures.
Checked over **every** blade of the genuine Cl(7) (28 L-blades, 35 F-blades). -/

/-- **Every L-grade blade is a complex structure** (`B² = −I`): all 21 bivectors
    and all 7 six-forms. -/
theorem L_grade_squares_neg : ∀ M ∈ all_L_grade, matMul M M = negId8 := by decide

/-- **Every F-grade blade is a real structure** (`B² = +I`): all 35 four-forms. -/
theorem F_grade_squares_pos : ∀ M ∈ all_F_fourforms, matMul M M = id8 := by decide

/-- `+I ≠ −I` in the 8×8 integer matrices (so the two structure types are genuinely
    disjoint, not a degenerate coincidence). -/
theorem id8_ne_negId8 : (id8 : Mat8) ≠ negId8 := by decide

/-- **No F-grade blade is a complex structure.**  Since every F-blade squares to
    `+I ≠ −I`, the F-grade cannot supply the `J` of `ξ = e^{iφ}`. -/
theorem F_no_complex_structure : ∀ M ∈ all_F_fourforms, ¬ isComplexStructure M := by
  intro M hM h
  exact id8_ne_negId8 (by rw [← F_grade_squares_pos M hM, h])

/-- **Every L-grade blade is a complex structure** (predicate form). -/
theorem L_all_complex_structures : ∀ M ∈ all_L_grade, isComplexStructure M :=
  L_grade_squares_neg

/-! ## Part 3 — the forcing: the lepton complex structure is L-grade -/

/-- **lepton = L is forced by the complex structure of the Brannen phase.**

    Bundling the facts:

      (1) the lepton ℍ-slice complex unit `i = e₀₁ = γ₀γ₁` is a genuine complex
          structure (`i² = −I`) that the Brannen phase `ξ = e^{iφ}` requires, and
          it is an L-grade (Λ²) blade (`i_lep_is_L_bivector`);
      (2) *every* L-grade blade is a complex structure (`B² = −I`), exhaustively
          over `Λ²⊕Λ⁶`;
      (3) *no* F-grade blade is a complex structure — every F-blade is a real
          structure (`B² = +I ≠ −I`).

    So the complex structure carrying the lepton's nontrivial Koide phase can only
    be L-grade, never F.  This is the missing ingredient flagged in
    `LeptonGradeForcing` (and `brannen_kernel.py`'s `Z₂×Z₂ undetermined`): the
    lepton's Hermitian/phase content is forced into L = Λ²⊕Λ⁶. -/
theorem lepton_complex_structure_forced_L :
    -- (1) the Brannen-phase complex unit is an L-grade complex structure
    (isComplexStructure i_lep ∧ i_lep ∈ all_L_bivectors)
    -- (2) L is exactly complex structures
    ∧ (∀ M ∈ all_L_grade, isComplexStructure M)
    -- (3) F carries no complex structure
    ∧ (∀ M ∈ all_F_fourforms, ¬ isComplexStructure M) :=
  ⟨⟨i_lep_isComplexStructure, i_lep_is_L_bivector⟩,
   L_all_complex_structures,
   F_no_complex_structure⟩

end SCPv59.Furey7D.LeptonComplexStructure
