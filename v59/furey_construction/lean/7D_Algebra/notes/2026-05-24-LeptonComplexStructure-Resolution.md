# lepton = L: the open half is closed (complex-structure forcing)

**Date**: 2026-05-24
**Modules**:
- `7D_Algebra/LeptonComplexStructure.lean` — concrete, exhaustive over genuine Cl(7).
- `BladeSquareSign.lean` (parent/Mathlib) — the representation-independent reason.
**Status**: builds in both packages; `decide`-checked / Mathlib-proved; 0 sorry; axiom-clean
(`propext`/`Quot.sound`, or none, for the concrete part; standard trio for the abstract).

## What the open half was
`LeptonGradeForcing.lean` proved lepton = L is *not* forced by mass-channel **availability**
(F is if anything richer), and isolated the missing ingredient:

> which grade carries the ℂ-imaginary unit `J` of Furey's minimal left ideal — if `J`
> lives in the Λ² bivector slice, lepton = L is forced.

## The resolution
The Brannen lepton mass kernel `M = a(I + ξS + ξ̄S²)` carries the phase `ξ = e^{iφ}`,
`φ = 2/9`.  A nontrivial phase (`Q = 2/3 ≠ 1`) **requires a genuine complex structure**
`J` (`J² = −1`) — the `i` of `cos φ + i sin φ`.  `brannen_kernel.py:9–16,123–137` embeds
the lepton ℍ-slice `(i,j,k)` as **three bivectors** `(e₀₁,e₀₂,e₁₂) ∈ Λ²`, so `J = e₀₁`.

### The forcing (the new content)
For a simple even k-blade, `B² = (−1)^{k(k+1)/2}·I`.  So in `Cl(7)_even = Λ⁰⊕Λ²⊕Λ⁴⊕Λ⁶`:

| grade | k(k+1)/2 | B² | type | sector |
|------|----------|----|------|--------|
| Λ⁰ (0) | 0  | +I | real (identity) | — |
| Λ² (2) | 3  | **−I** | **complex structure** | **L** |
| Λ⁴ (4) | 10 | +I | real (involution) | F |
| Λ⁶ (6) | 21 | **−I** | **complex structure** | **L** |

**`L = Λ²⊕Λ⁶` is exactly the even grades whose blades are complex structures; `F = Λ⁴`
furnishes only real structures.**  The lepton phase's `J` therefore *cannot* be F-grade.
**lepton = L is forced.**

This is the complex-structure face of `LeptonGradeForcing`'s singlet-pair dichotomy:
`J² = −I` ⇔ the antisymmetric ε/rotation block (L); `J² = +I` ⇔ the symmetric
metric/involution block (F).

## Lean content
Concrete (`LeptonComplexStructure.lean`, genuine Cl(7), `decide`):
- `i_lep_isComplexStructure` — `(γ₀γ₁)² = −I` (the Brannen-phase `i = e₀₁`).
- `lepton_slice_is_quaternion` — `(e₀₁,e₀₂,e₁₂)` is Hamilton's ℍ (`ij=k`, `jk=i`, `ki=j`, `ijk=−1`).
- `i_lep_is_L_bivector` — `e₀₁ ∈` the 21 Λ²-bivectors (so `J ∈ Λ² ⊂ L`).
- `L_grade_squares_neg` — **all 28** L-blades (Λ²⊕Λ⁶) square to `−I` (exhaustive).
- `F_grade_squares_pos` — **all 35** F-blades (Λ⁴) square to `+I` (exhaustive).
- `F_no_complex_structure` — no F-blade is a complex structure (since `+I ≠ −I`).
- `lepton_complex_structure_forced_L` — the bundle: `J` is an L-grade complex structure;
  L = complex structures; F has none. ⇒ lepton = L.

Universal (`BladeSquareSign.lean`, any ring, Mathlib):
- `prod_sq` — `(∏ l)² = (−1)^{tri |l|}` for pairwise-anticommuting square-(−1) elements
  (`tri n = n(n+1)/2`); the representation-independent `B² = (−1)^{k(k+1)/2}` law.
- `prod_sq_grade_two/four/six` + `even_grade_complex_structure_dichotomy` — grades 2,6 → `−1`
  (L), grade 4 → `+1` (F).  Not an 8×8-octonion accident: it is the parity of `k(k+1)/2`.

## Honest scope
- "Requires a complex structure" is the physical input (a nonzero Koide phase needs `J²=−1`);
  given that, the *grade* of `J` is forced by the algebra, which is what is proved.
- The concrete check is exhaustive over Cl(7)'s 63 even blades — complete for the actual
  algebra, not a sample.  The abstract law shows it holds in every Clifford algebra.
- Distinctness `−1 ≠ +1` is char-0 (holds for the ℝ/ℤ realization; would fail in char 2),
  so it is asserted in the concrete module (`id8_ne_negId8`), not the universal one.
