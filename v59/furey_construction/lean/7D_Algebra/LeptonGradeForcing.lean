import SevenDAlgebra

namespace SCPv59.Furey7D.LeptonForcing

open SCPv59.Furey7D

/-!
# The open half of the Z₂×Z₂ forcing: is lepton = L (Λ²⊕Λ⁶) forced over F (Λ⁴)?

The *composite* direction is already settled (`Z2Z2Forcing`): a color triplet
needs the rank-≥2 color Cartan that only the F-grade supplies, so d/u-quarks
require F; pure L (diagonal-free) cannot split colors.  This module addresses the
**other half**, the genuinely open question: the leptons are color *singlets*
(Fock N=0, N=3 → spinor indices 0 and 7), and are assigned to L.  But the naive
G₂ representation theory points the *other way*:

  * G₂ branching: spinor `8 = 1 ⊕ 7`, and `7|_{SU(3)} = 1 ⊕ 3 ⊕ 3̄`, so the two
    color singlets are the two G₂-singlet directions (indices 0,7).
  * The unique G₂-invariant in the even grades is the coassociative 4-form `*φ`,
    which lives in **F = Λ⁴**, *not* in L = Λ²⊕Λ⁶ (which has no G₂-singlet).

So a "the color singlet couples to the G₂-singlet" heuristic would put leptons in
**F**, contradicting the assignment.  This module resolves the tension *by exact
computation on the corrected Cl(7)*, and the honest verdict is:

  **lepton = L is NOT forced by mass-channel availability.**  In fact F offers a
  color-singlet (= color-neutral) diagonal mass channel to the leptons that L
  cannot (L is diagonal-free).  What *is* an exact structural fact is a clean
  **grade dichotomy on the lepton singlet pair {0,7}**:

      every L-generator acts ANTISYMMETRICALLY on {0,7}  (a Lie/rotation ε-block),
      every F-generator acts   SYMMETRICALLY on {0,7}     (a metric/mass block),

  together with the identification of F's *unique* color-neutral diagonal as the
  global Cl(6) chirality operator (not a lepton-isolating channel).

These are all `decide`-checked on the genuine (post sign-fix) Cl(7) generators.
What remains to *force* lepton=L (rather than merely make it consistent) is a
choice of what counts as a fundamental mass term (symmetric/Hermitian metric pairing
vs. an antisymmetric connection) — an input the 8×8 algebra alone does not fix.  See
the closing comment for the precise missing ingredient.
-/

/-! ## Lepton-block accessors on the {0,7} singlet pair -/

/-- The `(0,7)` entry of an 8×8 matrix. -/
def e07 (M : Mat8) : Int := listGetD (listGetD M 0 []) 7 0
/-- The `(7,0)` entry of an 8×8 matrix. -/
def e70 (M : Mat8) : Int := listGetD (listGetD M 7 []) 0 0
/-- The `(0,0)` entry. -/
def e00 (M : Mat8) : Int := listGetD (listGetD M 0 []) 0 0
/-- The `(7,7)` entry. -/
def e77 (M : Mat8) : Int := listGetD (listGetD M 7 []) 7 0

/-- Rebuild the list of all 28 L-grade generators (Λ² ⊕ Λ⁶) here so this module is
    standalone (mirrors `Forcing.all_L_grade`). -/
def hexads7' : List (List Nat) :=
  (List.range 7).map (fun d => (List.range 7).filter (fun k => k != d))
def gammaProd' (ks : List Nat) : Mat8 := ks.foldl (fun acc k => matMul acc (gamma k)) id8
def all_L_sixforms' : List Mat8 := hexads7'.map gammaProd'
def all_L_grade' : List Mat8 := all_L_bivectors ++ all_L_sixforms'

/-! ## Fact 1 — L is diagonal-free on leptons, but F is NOT.

This is the fact that *refutes* the naive "leptons must be diagonal-free, hence L"
intuition: the F-grade supplies a genuine nonzero diagonal mass channel to the
lepton singlets.  So a lepton mass from F is not excluded by availability. -/

/-- Every L-grade generator has zero diagonal on both lepton indices (0 and 7). -/
theorem L_lepton_diag_zero :
    ∀ M ∈ all_L_grade', e00 M = 0 ∧ e77 M = 0 := by decide

/-- The F-grade is NOT diagonal-free on the leptons: some 4-form has a nonzero
    lepton diagonal — concretely `γ₃γ₄γ₅γ₆`.  Hence "needs zero diagonal ⇒ L" fails. -/
theorem F_lepton_diag_nonzero :
    ∃ M ∈ all_F_fourforms, e00 M ≠ 0 ∨ e77 M ≠ 0 :=
  ⟨matMul (matMul (gamma 3) (gamma 4)) (matMul (gamma 5) (gamma 6)),
   by decide, by decide⟩

/-! ## Fact 2 — the grade dichotomy on the lepton singlet pair {0,7}.

This is the clean, exact structural statement: the lepton-pair block is purely
antisymmetric for *every* L-generator and purely symmetric for *every* F-generator.
(`e₀₀=e₇₇=0`, `e₀₇=−e₇₀` for all of L; `e₀₇=e₇₀` for all of F.) -/

/-- **L acts antisymmetrically on the lepton pair {0,7}**: zero diagonal and
    `(0,7) = −(7,0)` — the ε / SO(2) rotation (Lie-algebra) block. -/
theorem L_lepton_block_antisymmetric :
    ∀ M ∈ all_L_grade', e00 M = 0 ∧ e77 M = 0 ∧ e07 M = - e70 M := by decide

/-- **F acts symmetrically on the lepton pair {0,7}**: `(0,7) = (7,0)` — the
    symmetric (metric / mass) block. -/
theorem F_lepton_block_symmetric :
    ∀ M ∈ all_F_fourforms, e07 M = e70 M := by decide

/-- The dichotomy is *strict*: an explicit L-generator has a nonzero antisymmetric
    lepton block (`(0,7)=1, (7,0)=−1`), while an explicit F-generator has a nonzero
    *symmetric* lepton block (`(0,7)=(7,0)=1`).  So the two grades furnish genuinely
    different (and both nonzero) lepton-pair couplings. -/
theorem lepton_block_dichotomy_strict :
    (∃ M ∈ all_L_grade', e07 M = 1 ∧ e70 M = -1) ∧
    (∃ M ∈ all_F_fourforms, e07 M = 1 ∧ e70 M = 1) :=
  ⟨⟨matMul (gamma 0) (gamma 5), by decide, by decide⟩,
   ⟨matMul (matMul (gamma 0) (gamma 1)) (matMul (gamma 2) (gamma 6)), by decide, by decide⟩⟩

/-! ## Fact 3 — F's unique color-neutral diagonal is the global Cl(6) chirality.

The element `γ₃γ₄γ₅γ₆` is the *only* F-generator giving a lepton diagonal that does
*not* split colors; but it equals `diag(1,1,1,1,−1,−1,−1,−1)`, the chirality operator
splitting the spinor 8 into its two Weyl halves (Fock N≤1 vs N≥2).  It acts identically
on the quark blocks, so it is a *global* grading, not a lepton-isolating mass.  Thus
even the one color-neutral F diagonal cannot be switched on for leptons alone. -/

/-- `F_chi := γ₃γ₄γ₅γ₆`, the color-neutral diagonal element of F. -/
def F_chi : Mat8 := matMul (matMul (gamma 3) (gamma 4)) (matMul (gamma 5) (gamma 6))

theorem F_chi_mem : F_chi ∈ all_F_fourforms := by decide

/-- `F_chi` is exactly the Cl(6) chirality operator `diag(1,1,1,1,−1,−1,−1,−1)`:
    it gives `+1` on the Fock-N≤1 states {0,1,2,3} and `−1` on the N≥2 states
    {4,5,6,7}.  In particular its lepton-block values are *opposite* (`+1` on N=0,
    `−1` on N=3) and tied to the *same* operator that grades the quarks. -/
theorem F_chi_is_chirality :
    diag F_chi = [1, 1, 1, 1, -1, -1, -1, -1] := by decide

/-- Consequently `F_chi` does NOT isolate the leptons: it acts (nontrivially and
    identically per Weyl half) on the quark sectors too.  Its d-quark block is the
    *constant* `[1,1,1]` (the N=1 states are all in the `+1` half) — color-neutral,
    but inseparable from the global grading. -/
theorem F_chi_acts_on_quarks :
    sectorDiags F_chi dQuarkIndices = [1, 1, 1] ∧
    sectorDiags F_chi uQuarkIndices = [-1, -1, -1] := by decide

/-! ## The honest verdict (machine-checked content) -/

/-- **lepton = L is NOT forced by mass-channel availability.**  Bundling the facts:

      (1) L is diagonal-free on the leptons (`L_lepton_diag_zero`) — L offers the
          leptons *no* diagonal (Majorana-type) mass channel at all;
      (2) F *does* offer the leptons a nonzero diagonal channel
          (`F_lepton_diag_nonzero`);
      (3) on the singlet pair {0,7}, L is purely antisymmetric and F purely
          symmetric (`L_lepton_block_antisymmetric`, `F_lepton_block_symmetric`).

    So both grades can furnish a lepton coupling, and if anything F is *richer*
    (it adds the diagonal that L lacks).  The assignment lepton=L is therefore a
    choice about *which kind of pairing counts as the fundamental mass* (the
    antisymmetric connection/ε-block of L vs. the symmetric metric block of F),
    not a consequence the 8×8 algebra forces.  This is the precise, machine-checked
    statement of the non-forcing. -/
theorem lepton_L_not_forced_by_availability :
    (∀ M ∈ all_L_grade', e00 M = 0 ∧ e77 M = 0) ∧
    (∃ M ∈ all_F_fourforms, e00 M ≠ 0 ∨ e77 M ≠ 0) ∧
    (∀ M ∈ all_L_grade', e07 M = - e70 M) ∧
    (∀ M ∈ all_F_fourforms, e07 M = e70 M) :=
  ⟨L_lepton_diag_zero,
   F_lepton_diag_nonzero,
   fun M hM => (L_lepton_block_antisymmetric M hM).2.2,
   F_lepton_block_symmetric⟩

/-! ## What WOULD force lepton = L — the precise missing ingredient

The computation pins the open question down to a single, well-posed input that the
8×8 representation does **not** supply:

  > Declare that the fundamental fermion **mass** term is the **symmetric/Hermitian
  > metric pairing** that, *after* the ℂ (or ℍ) structure of the Furey ideal is
  > imposed, has *real* eigenvalues — i.e. a Brannen kernel `M = a(I + ξS + ξ̄S²)`
  > that is Hermitian (`08_brannen_yukawa.py`).

Under the real 8×8 picture computed here, the symmetric lepton block is supplied by
**F** and the antisymmetric block by **L**.  So at face value the *real* metric mass
sits in F — the opposite of the assignment.  The assignment lepton=L can only be
forced once one fixes the complex structure `J` (the imaginary unit of ℂ in
`ℂ⊗ℍ⊗𝕆`) under which `J · (antisymmetric L-block) = Hermitian mass`.  Equivalently,
the *missing ingredient* is:

  **a derivation of which grade carries the ℂ-imaginary unit `J` of the Furey ideal.**

  • If `J ∈ L` (the bivector / `Λ²` slice, as `brannen_kernel.py` embeds the lepton
    ξ), then the L antisymmetric block `[[0,1],[−1,0]] = J|_{0,7}` is the genuine
    Hermitian lepton mass, and lepton=L *is* forced.
  • If `J ∈ F`, the assignment flips.

So the forcing reduces to: *show the complex structure `J` (Furey's `i`) lives in the
Λ² bivector grade.*  That is a statement about the ℂ-structure of the minimal left
ideal — outside the scope of the real 8×8 Cl(7) matrices — and is the concrete next
target.  Until it is established, lepton=L is **consistent and natural but not forced**
by the algebra computed here; the composite (quark) half of the Z₂×Z₂, by contrast,
*is* forced (`Z2Z2Forcing.composite_color_requires_LF`). -/

end SCPv59.Furey7D.LeptonForcing
