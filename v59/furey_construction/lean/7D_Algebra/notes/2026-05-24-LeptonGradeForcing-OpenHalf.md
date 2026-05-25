# Lepton = L vs F: the open half of the Z₂×Z₂ forcing (corrected Cl(7))

> **RESOLVED 2026-05-24** — the missing ingredient (§"The precise missing ingredient"
> below) is now supplied and the forcing completed in
> `LeptonComplexStructure.lean` + `BladeSquareSign.lean`.  See
> `2026-05-24-LeptonComplexStructure-Resolution.md`.  Summary: `L = Λ²⊕Λ⁶` is exactly
> the even grades whose blades square to `−I` (complex structures), `F = Λ⁴` squares to
> `+I` (real); the Brannen phase `ξ = e^{iφ}` needs a complex structure `J` (`J²=−1`),
> which is therefore L-grade — and `J = e₀₁ ∈ Λ²` concretely.  This document records the
> *open* state and remains accurate as the statement of the problem.

**Date**: 2026-05-24
**Module**: `LeptonGradeForcing.lean` (added to `Furey7D` lib; builds clean, 0 sorry/axiom).
**Scope**: the *singlet* half of the sector forcing — why leptons (Fock N=0,3 → spinor
indices 0,7) are assigned to L=Λ²⊕Λ⁶ and not F=Λ⁴. (The *composite/quark* half is already
forced by `Z2Z2Forcing.composite_color_requires_LF`; not touched.)

## The tension (recap)
G₂ branching: spinor 8 = 1⊕7, 7|_SU(3) = 1⊕3⊕3̄, so the two color singlets are the two
G₂-singlet directions. The unique even-grade G₂-invariant is the coassociative 4-form *φ ∈
**F**, NOT in L (which has no G₂-singlet). A "singlet couples to singlet" heuristic would
put leptons in **F** — opposite to the assignment.

## What was computed (all on the corrected, genuine Cl(7); `decide`-checked)
1. `L_lepton_diag_zero` — L is diagonal-free on indices 0 and 7 (L offers leptons NO
   diagonal/Majorana mass channel).
2. `F_lepton_diag_nonzero` — F is NOT diagonal-free on leptons (γ₃γ₄γ₅γ₆ gives a nonzero
   lepton diagonal). So "needs zero diagonal ⇒ L" is refuted; F is *richer* for leptons.
3. **Grade dichotomy on the {0,7} pair** (the clean exact fact):
   - `L_lepton_block_antisymmetric` — every L-generator: e₀₀=e₇₇=0 and e₀₇=−e₇₀ (ε / SO(2)
     rotation block — a Lie/connection coupling).
   - `F_lepton_block_symmetric` — every F-generator: e₀₇=e₇₀ (symmetric metric/mass block).
   - `lepton_block_dichotomy_strict` — both blocks are realized nonzero (witnesses γ₀γ₅ ∈ L,
     γ₀γ₁γ₂γ₆ ∈ F).
4. F's *unique* color-neutral diagonal is the global chirality operator:
   - `F_chi_is_chirality` — γ₃γ₄γ₅γ₆ = diag(1,1,1,1,−1,−1,−1,−1) (Cl(6) Weyl split, N≤1 vs N≥2).
   - `F_chi_acts_on_quarks` — it acts identically per Weyl half on the quark blocks, so it
     is NOT a lepton-isolating channel.
5. `lepton_L_not_forced_by_availability` — bundles (1)–(3): both grades furnish a lepton
   coupling; availability does not force lepton=L.

## Verdict
**lepton = L is NOT forced by the 8×8 mass-channel availability.** It is *consistent and
natural* but not forced. The earlier 7-angle report's "rep-theoretic forcing" was
interpretive prose on the buggy 3D model; on the genuine Cl(7) the diagonal-availability
argument actually points the other way (F has more lepton structure).

## The precise missing ingredient
The forcing reduces to a **complex-structure** question the real 8×8 Cl(7) does not fix:
which grade carries the ℂ-imaginary unit J of Furey's minimal left ideal in ℂ⊗ℍ⊗𝕆?
- The real L lepton-block is antisymmetric ([[0,1],[−1,0]] = J|_{0,7}); the real F
  lepton-block is symmetric. The physical Brannen lepton mass is *Hermitian*
  (`08_brannen_yukawa.py`: M=a(I+ξS+ξ̄S²), eigenvalues real).
- Hermitian = (real symmetric) directly, OR (real antisymmetric)·J after imposing J.
- So: **if J lives in the Λ² bivector slice** (as `brannen_kernel.py` *embeds* the lepton ξ),
  then L's antisymmetric block IS the genuine Hermitian lepton mass and lepton=L is forced.
  If J lives in F, the assignment flips.
- Concrete next target: derive that the complex structure J (Furey's i) sits in Λ² — a
  statement about the ℂ-structure of the ideal, outside the real 8×8 matrices.

Note: `08_brannen_yukawa.py:328` itself flags "Z₂×Z₂ pattern undetermined", consistent with
this verdict.
