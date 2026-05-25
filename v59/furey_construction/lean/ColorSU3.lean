/-
  v59/furey_construction/lean/ColorSU3.lean

  **Steps (C1) and (C2)** of the complex-structure shape-constraint program
  (`7D_Algebra/notes/2026-05-24-StepC-and-phase-2-9.md`): pin the lepton's complex
  structure `J` to a unique element via the color SU(3)-commutant.

  Parts (A) and (B) (`LeptonRealityForcing.lean`) placed `J` in `L = Λ²⊕Λ⁶` (the skew /
  so(8) subspace).  Here we add the color symmetry and pin `J` on the lepton sector.

  ## The correct color complex structure (a correction)

  The *generation* ℍ-slice unit `γ₀γ₁` of `LeptonComplexStructure` does **not** preserve the
  color-singlet sector.  The **color** complex structure is `J_c = γ₀γ₅`: it is block-diagonal
  on the Fock/color sectors and pairs the states as `(0↔7), (1↔6), (2↔5), (3↔4)`, making
  `ℝ⁸ = ℂ⁴ = 1 ⊕ 3` — the lepton singlet `{0,7}` as the complex `1`, and the three quark
  modes `{(1,6),(2,5),(3,4)}` as the `3` (the Furey minimal-left-ideal split).  `J_c ∈ Λ² ⊂ L`.

  ## (C1) An explicit color `𝔰𝔲(3)`

  Eight generators (Cartan `H₁,H₂`; roots `A_ab, B_ab`) act on the three quark modes (as the
  standard `𝔲(3)→𝔰𝔲(3)`), are skew (`⊂ so(8)`), **annihilate the lepton singlet `{0,7}`**
  (so the lepton is a genuine color singlet), **commute with `J_c`** (so `J_c` is color-invariant),
  and **close** into `𝔰𝔲(3)` (the full A₂ structure-constant table, all coeffs `∈{±1,±2}`).

  ## (C2) Pinning the lepton complex structure

  The Cartan common kernel is *exactly* the lepton sector: `ker H₁ ∩ ker H₂ = span{e₀,e₇}`.
  Hence any operator commuting with `H₁,H₂` preserves `span{e₀,e₇}`; a *skew* (reality)
  complex structure then restricts there to a 2×2 skew complex structure, which is forced to
  `±[[0,1],[−1,0]]` (`skew_cs_2x2`).  So the lepton's complex structure is pinned, up to
  orientation, to the canonical form realized by `J_c` — `lepton_sector_pinned`.

  All gamma-model facts are `decide`-checked; `skew_cs_2x2` is a Mathlib `Fin 2` argument.
-/
import Mathlib
import SevenDAlgebra
import LeptonComplexStructure
import LeptonRealityForcing

namespace SCPv59.Furey7D.ColorSU3

open SCPv59.Furey7D
open SCPv59.Furey7D.LeptonComplexStructure (negId8 isComplexStructure)
open SCPv59.Furey7D.LeptonRealityForcing (ent isSkew mulSelf00 ent_matAdd ent_matScale)

/-! ## The color complex structure `J_c = γ₀γ₅` -/

def Jc : Mat8 := matMul (gamma 0) (gamma 5)

/-- `J_c² = −I`: a genuine complex structure. -/
theorem Jc_complexStructure : isComplexStructure Jc := by decide

/-- `J_c` is skew (in `L = so(8)`). -/
theorem Jc_isSkew : isSkew Jc := by decide

/-- `J_c` acts as the canonical complex structure `[[0,1],[−1,0]]` on the lepton sector:
    `e₀ ↦ −e₇`, `e₇ ↦ e₀`, with zero diagonal. -/
theorem Jc_lepton_block :
    ent Jc 0 0 = 0 ∧ ent Jc 7 7 = 0 ∧ ent Jc 7 0 = -1 ∧ ent Jc 0 7 = 1 := by decide

/-- `J_c` is block-diagonal on the color sectors: it preserves the lepton sector `{0,7}`
    (columns 0,7 supported on rows {0,7}) and the quark sector `{1..6}`. -/
theorem Jc_block_diagonal :
    (∀ r, r < 8 → r ≠ 0 → r ≠ 7 → ent Jc r 0 = 0 ∧ ent Jc r 7 = 0)
    ∧ (∀ c, c < 8 → c ≠ 0 → c ≠ 7 → ent Jc 0 c = 0 ∧ ent Jc 7 c = 0) := by decide

/-! ## (C1) The eight `𝔰𝔲(3)` generators -/

/-- Build an 8×8 matrix from a list of `(row, col, value)` entries. -/
def mkMat (es : List (Nat × Nat × Int)) : Mat8 :=
  (List.range 8).map (fun r => (List.range 8).map (fun c =>
    (es.filter (fun e => e.1 == r && e.2.1 == c)).foldl (fun a e => a + e.2.2) 0))

/-- Mode-mixing generator `A_ab` (real, `E_ab − E_ba` on the complex quark modes). -/
def Amat (ra ja rb jb : Nat) : Mat8 := mkMat [(rb,ra,1),(jb,ja,1),(ra,rb,-1),(ja,jb,-1)]
/-- Mode-mixing generator `B_ab` (`i(E_ab + E_ba)`). -/
def Bmat (ra ja rb jb : Nat) : Mat8 := mkMat [(jb,ra,1),(rb,ja,-1),(ja,rb,1),(ra,jb,-1)]
/-- Single-mode `U(1)` generator `T_a` (the `J_c` of mode `a`). -/
def Tmat (ra ja : Nat) : Mat8 := mkMat [(ja,ra,1),(ra,ja,-1)]

/-- The 3 quark modes under `J_c`: `m₁=(1,6)`, `m₂=(2,5)`, `m₃=(3,4)`. -/
def H1 : Mat8 := matAdd (Tmat 1 6) (matScale (-1) (Tmat 2 5))  -- T₁−T₂  (Cartan)
def H2 : Mat8 := matAdd (Tmat 2 5) (matScale (-1) (Tmat 3 4))  -- T₂−T₃  (Cartan)
def A12 : Mat8 := Amat 1 6 2 5
def A13 : Mat8 := Amat 1 6 3 4
def A23 : Mat8 := Amat 2 5 3 4
def B12 : Mat8 := Bmat 1 6 2 5
def B13 : Mat8 := Bmat 1 6 3 4
def B23 : Mat8 := Bmat 2 5 3 4

/-- The 8 color `𝔰𝔲(3)` generators. -/
def colorGens : List Mat8 := [H1, H2, A12, A13, A23, B12, B13, B23]

/-- Lie bracket (commutator) of two operators. -/
def commM (X Y : Mat8) : Mat8 := matAdd (matMul X Y) (matScale (-1) (matMul Y X))

/-- An operator **annihilates the lepton singlet** `{0,7}` iff its columns 0,7 vanish. -/
def killsLepton (M : Mat8) : Prop := ∀ r, r < 8 → ent M r 0 = 0 ∧ ent M r 7 = 0

instance (M : Mat8) : Decidable (killsLepton M) := by unfold killsLepton; infer_instance

/-- **Every color generator is skew** (`⊂ so(8) = L`). -/
theorem color_isSkew : ∀ G ∈ colorGens, isSkew G := by decide

/-- **Every color generator annihilates the lepton singlet `{0,7}`** — the lepton is a
    genuine color singlet (trivial `𝔰𝔲(3)`-rep). -/
theorem color_kills_lepton : ∀ G ∈ colorGens, killsLepton G := by decide

/-- **`J_c` is color-invariant**: it commutes with every `𝔰𝔲(3)` generator. -/
theorem color_commutes_Jc : ∀ G ∈ colorGens, commM G Jc = zero8 := by decide

/-! ### Closure: the full A₂ structure-constant table (so the 8 generators are `𝔰𝔲(3)`) -/

/-- **The `𝔰𝔲(3)` bracket relations** — all 28 upper-triangular commutators land back in the
    span of the 8 generators, with the standard A₂ structure constants (`∈ {±1,±2}`).  This
    certifies the 8 generators form a Lie subalgebra isomorphic to `𝔰𝔲(3)`. -/
theorem colorSU3_closes :
    commM H1 H2 = zero8
    ∧ commM H1 A12 = matScale (-2) B12 ∧ commM H1 A13 = matScale (-1) B13
    ∧ commM H1 A23 = B23 ∧ commM H1 B12 = matScale 2 A12 ∧ commM H1 B13 = A13
    ∧ commM H1 B23 = matScale (-1) A23
    ∧ commM H2 A12 = B12 ∧ commM H2 A13 = matScale (-1) B13 ∧ commM H2 A23 = matScale (-2) B23
    ∧ commM H2 B12 = matScale (-1) A12 ∧ commM H2 B13 = A13 ∧ commM H2 B23 = matScale 2 A23
    ∧ commM A12 A13 = A23 ∧ commM A12 A23 = matScale (-1) A13
    ∧ commM A12 B12 = matScale (-2) H1 ∧ commM A12 B13 = B23 ∧ commM A12 B23 = matScale (-1) B13
    ∧ commM A13 A23 = A12 ∧ commM A13 B12 = B23
    ∧ commM A13 B13 = matAdd (matScale (-2) H1) (matScale (-2) H2)
    ∧ commM A13 B23 = matScale (-1) B12
    ∧ commM A23 B12 = B13 ∧ commM A23 B13 = matScale (-1) B12 ∧ commM A23 B23 = matScale (-2) H2
    ∧ commM B12 B13 = A23 ∧ commM B12 B23 = A13 ∧ commM B13 B23 = A12 := by
  refine ⟨?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_,
         ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_⟩ <;> decide

/-! ## (C2) Pinning the lepton complex structure

The Cartan common kernel is exactly the lepton sector, so a color-invariant complex
structure preserves `{0,7}`; reality (skewness) then pins its restriction there. -/

/-- The lepton sector `{0,7}` is **preserved** by `M` iff `M` maps `e₀,e₇` back into
    `span{e₀,e₇}` (columns 0,7 supported on rows {0,7}). -/
def preservesLepton (M : Mat8) : Prop :=
  ∀ k, k < 8 → k ≠ 0 → k ≠ 7 → ent M k 0 = 0 ∧ ent M k 7 = 0

instance (M : Mat8) : Decidable (preservesLepton M) := by unfold preservesLepton; infer_instance

/-- The quark part of `M`'s column 0 vanishes (`M e₀ ∈ span{e₀,e₇}`). -/
def quarkCol0Zero (M : Mat8) : Prop := ∀ k, k < 8 → k ≠ 0 → k ≠ 7 → ent M k 0 = 0

instance (M : Mat8) : Decidable (quarkCol0Zero M) := by unfold quarkCol0Zero; infer_instance

/-- **The 2×2 pinning (the heart of C2).**  A skew complex structure whose column 0 is
    supported on the lepton sector restricts there to `±[[0,1],[−1,0]]`: `ent J 0 7 = ±1`.
    Proof: `(J·J)₀₀ = −1` (complex structure) collapses — column-0 support kills the quark
    cross terms, skewness kills the diagonal — to `−(J₀₇)² = −1`, so `J₀₇ = ±1`. -/
theorem lepton_block_canonical (J : Mat8)
    (hsk : isSkew J) (hcs : isComplexStructure J) (hpres : quarkCol0Zero J) :
    ent J 0 7 = 1 ∨ ent J 0 7 = -1 := by
  have h00 : ent (matMul J J) 0 0 = -1 := by
    show ent (matMul J J) 0 0 = ent negId8 0 0
    rw [show matMul J J = negId8 from hcs]
  rw [mulSelf00] at h00
  have e10 : ent J 1 0 = 0 := hpres 1 (by norm_num) (by norm_num) (by norm_num)
  have e20 : ent J 2 0 = 0 := hpres 2 (by norm_num) (by norm_num) (by norm_num)
  have e30 : ent J 3 0 = 0 := hpres 3 (by norm_num) (by norm_num) (by norm_num)
  have e40 : ent J 4 0 = 0 := hpres 4 (by norm_num) (by norm_num) (by norm_num)
  have e50 : ent J 5 0 = 0 := hpres 5 (by norm_num) (by norm_num) (by norm_num)
  have e60 : ent J 6 0 = 0 := hpres 6 (by norm_num) (by norm_num) (by norm_num)
  have d00 : ent J 0 0 = 0 := by have := hsk 0 (by norm_num) 0 (by norm_num); omega
  have a70 : ent J 7 0 = - ent J 0 7 := hsk 7 (by norm_num) 0 (by norm_num)
  rw [e10, e20, e30, e40, e50, e60, d00, a70] at h00
  have hsq : ent J 0 7 * ent J 0 7 = 1 := by nlinarith [h00]
  rcases mul_self_eq_one_iff.mp hsq with h | h
  · exact Or.inl h
  · exact Or.inr h

/-! ### The bridge: color-invariance ⇒ the lepton sector is preserved

The Cartan common kernel is exactly the lepton sector (`H₁` kills `{e₀,e₃,e₄,e₇}`, `H₂` kills
`{e₀,e₁,e₆,e₇}`, intersection `{e₀,e₇}`), so an operator commuting with `H₁,H₂` maps `e₀`
into that kernel `= span{e₀,e₇}`. -/

/-- `replicate` of a constant reads back as the constant at any index. -/
private lemma listGetD_repl_const {α} (a : α) : ∀ (n r : Nat),
    listGetD (List.replicate n a) r a = a
  | 0, _ => rfl
  | _ + 1, 0 => rfl
  | n + 1, r + 1 => by
      show listGetD (List.replicate n a) r a = a; exact listGetD_repl_const a n r

/-- Every entry of `zero8` is `0`. -/
private lemma ent_zero8 (r c : Nat) : ent zero8 r c = 0 := by
  show listGetD (listGetD zero8 r (List.replicate 8 0)) c 0 = 0
  rw [show zero8 = List.replicate 8 (List.replicate 8 (0 : Int)) from rfl,
      listGetD_repl_const, LeptonRealityForcing.listGetD_repl0]

/-- General `(r,c)` entry of a matrix product as the explicit length-8 dot product. -/
theorem ent_matMul (A B : Mat8) (r c : Nat) (hr : r < 8) (hc : c < 8) :
    ent (matMul A B) r c
      = ent A r 0 * ent B 0 c + ent A r 1 * ent B 1 c + ent A r 2 * ent B 2 c
      + ent A r 3 * ent B 3 c + ent A r 4 * ent B 4 c + ent A r 5 * ent B 5 c
      + ent A r 6 * ent B 6 c + ent A r 7 * ent B 7 c := by
  interval_cases r <;> interval_cases c <;>
    (simp only [ent, matMul, listGetD, listGetOpt, List.range, List.range.loop,
      List.foldl, List.map]; ring)

/-- From `[J, Hᵢ] = 0` with `Hᵢ` killing `e₀` (column 0 zero), the row-`r` dot product
    `Σₘ (Hᵢ)ᵣₘ Jₘ₀ = 0` — i.e. `Hᵢ (J e₀) = 0`, so `J e₀ ∈ ker Hᵢ`. -/
private lemma dotRow_zero (J Hi : Mat8) (hcomm : commM J Hi = zero8) (hk : killsLepton Hi)
    (r : Nat) (hr : r < 8) :
    ent Hi r 0 * ent J 0 0 + ent Hi r 1 * ent J 1 0 + ent Hi r 2 * ent J 2 0
    + ent Hi r 3 * ent J 3 0 + ent Hi r 4 * ent J 4 0 + ent Hi r 5 * ent J 5 0
    + ent Hi r 6 * ent J 6 0 + ent Hi r 7 * ent J 7 0 = 0 := by
  have hz : ent (commM J Hi) r 0 = 0 := by rw [hcomm]; exact ent_zero8 r 0
  rw [commM, ent_matAdd _ _ r 0 hr (by norm_num), ent_matScale _ _ r 0 hr (by norm_num),
      ent_matMul J Hi r 0 hr (by norm_num), ent_matMul Hi J r 0 hr (by norm_num)] at hz
  rw [(hk 0 (by norm_num)).1, (hk 1 (by norm_num)).1, (hk 2 (by norm_num)).1,
      (hk 3 (by norm_num)).1, (hk 4 (by norm_num)).1, (hk 5 (by norm_num)).1,
      (hk 6 (by norm_num)).1, (hk 7 (by norm_num)).1] at hz
  linarith [hz]

open LeptonRealityForcing in
/-- **The bridge.**  An operator commuting with the color Cartan `H₁,H₂` maps `e₀` into the
    Cartan common kernel `= span{e₀,e₇}`: its column-0 quark entries vanish. -/
theorem colorInvariant_quarkCol0Zero (J : Mat8)
    (h1 : commM J H1 = zero8) (h2 : commM J H2 = zero8) : quarkCol0Zero J := by
  have hk1 : killsLepton H1 := color_kills_lepton H1 (by simp [colorGens])
  have hk2 : killsLepton H2 := color_kills_lepton H2 (by simp [colorGens])
  -- H₁ rows isolate J₁₀,J₆₀,J₂₀,J₅₀ ;  H₂ rows isolate J₃₀,J₄₀
  have r6 := dotRow_zero J H1 h1 hk1 6 (by norm_num)
  have r1 := dotRow_zero J H1 h1 hk1 1 (by norm_num)
  have r5 := dotRow_zero J H1 h1 hk1 5 (by norm_num)
  have r2 := dotRow_zero J H1 h1 hk1 2 (by norm_num)
  have s4 := dotRow_zero J H2 h2 hk2 4 (by norm_num)
  have s3 := dotRow_zero J H2 h2 hk2 3 (by norm_num)
  simp only [show ent H1 6 0 = 0 from by decide, show ent H1 6 1 = 1 from by decide,
    show ent H1 6 2 = 0 from by decide, show ent H1 6 3 = 0 from by decide,
    show ent H1 6 4 = 0 from by decide, show ent H1 6 5 = 0 from by decide,
    show ent H1 6 6 = 0 from by decide, show ent H1 6 7 = 0 from by decide] at r6
  simp only [show ent H1 1 0 = 0 from by decide, show ent H1 1 1 = 0 from by decide,
    show ent H1 1 2 = 0 from by decide, show ent H1 1 3 = 0 from by decide,
    show ent H1 1 4 = 0 from by decide, show ent H1 1 5 = 0 from by decide,
    show ent H1 1 6 = -1 from by decide, show ent H1 1 7 = 0 from by decide] at r1
  simp only [show ent H1 5 0 = 0 from by decide, show ent H1 5 1 = 0 from by decide,
    show ent H1 5 2 = -1 from by decide, show ent H1 5 3 = 0 from by decide,
    show ent H1 5 4 = 0 from by decide, show ent H1 5 5 = 0 from by decide,
    show ent H1 5 6 = 0 from by decide, show ent H1 5 7 = 0 from by decide] at r5
  simp only [show ent H1 2 0 = 0 from by decide, show ent H1 2 1 = 0 from by decide,
    show ent H1 2 2 = 0 from by decide, show ent H1 2 3 = 0 from by decide,
    show ent H1 2 4 = 0 from by decide, show ent H1 2 5 = 1 from by decide,
    show ent H1 2 6 = 0 from by decide, show ent H1 2 7 = 0 from by decide] at r2
  simp only [show ent H2 4 0 = 0 from by decide, show ent H2 4 1 = 0 from by decide,
    show ent H2 4 2 = 0 from by decide, show ent H2 4 3 = -1 from by decide,
    show ent H2 4 4 = 0 from by decide, show ent H2 4 5 = 0 from by decide,
    show ent H2 4 6 = 0 from by decide, show ent H2 4 7 = 0 from by decide] at s4
  simp only [show ent H2 3 0 = 0 from by decide, show ent H2 3 1 = 0 from by decide,
    show ent H2 3 2 = 0 from by decide, show ent H2 3 3 = 0 from by decide,
    show ent H2 3 4 = 1 from by decide, show ent H2 3 5 = 0 from by decide,
    show ent H2 3 6 = 0 from by decide, show ent H2 3 7 = 0 from by decide] at s3
  -- r6: J₁₀=0, r1: −J₆₀=0, r5: −J₂₀=0, r2: J₅₀=0, s4: −J₃₀=0, s3: J₄₀=0
  intro k hk hk0 hk7
  interval_cases k <;> first | (exfalso; exact hk0 rfl) | (exfalso; exact hk7 rfl) | linarith

/-- **C2 — color-invariance pins the lepton complex structure.**  Any skew (reality) complex
    structure commuting with the color Cartan `H₁,H₂` has lepton block `±[[0,1],[−1,0]]`
    (`ent J 0 7 = ±1`).  So the lepton's complex structure is unique up to orientation — a
    single SU(3)-orbit — and `J_c = γ₀γ₅` realizes it. -/
theorem colorInvariant_pins_lepton (J : Mat8)
    (hsk : isSkew J) (hcs : isComplexStructure J)
    (h1 : commM J H1 = zero8) (h2 : commM J H2 = zero8) :
    ent J 0 7 = 1 ∨ ent J 0 7 = -1 :=
  lepton_block_canonical J hsk hcs (colorInvariant_quarkCol0Zero J h1 h2)

/-- `J_c = γ₀γ₅` preserves the lepton sector. -/
theorem Jc_preservesLepton : preservesLepton Jc := by decide

/-- **`J_c` realizes the canonical lepton complex structure**, as a color-invariant skew
    complex structure — the distinguished fixed point of the pinning. -/
theorem Jc_pinned : ent Jc 0 7 = 1 ∨ ent Jc 0 7 = -1 :=
  lepton_block_canonical Jc Jc_isSkew Jc_complexStructure (by decide)

/-- The lepton sector is in the Cartan kernel (`H₁,H₂` kill `e₀,e₇`). -/
theorem cartan_kills_lepton : killsLepton H1 ∧ killsLepton H2 :=
  ⟨color_kills_lepton H1 (by simp [colorGens]), color_kills_lepton H2 (by simp [colorGens])⟩

/-- No quark direction is in the Cartan common kernel: for each `k ∈ {1,…,6}` some `Hᵢ` has a
    nonzero column-`k` entry.  With `cartan_kills_lepton`, the common kernel is exactly
    `{e₀,e₇}`. -/
theorem cartan_quark_kernel_trivial :
    ∀ k, k < 8 → k ≠ 0 → k ≠ 7 →
      (∃ r, r < 8 ∧ ent H1 r k ≠ 0) ∨ (∃ r, r < 8 ∧ ent H2 r k ≠ 0) := by decide

/-! ## The quark single orbit

Schur on the (real-irreducible, complex-type) `3⊕3̄` quark block: a color-invariant skew
complex structure acts as `±J_c` on the quarks too.  Concretely, from `J e₁ = b·J_c e₁`
(seed), commuting with the 5 generators that send `e₁ ↦ e₂,e₃,e₄,e₅,e₆` propagates
`J = b·J_c` across all six quark columns; `J² = −I` forces `b = ±1`. -/

/-- **The symmetric commutator-entry identity.**  `[J,G] = 0` ⇒ for any `r,c < 8`,
    `Σₘ Jᵣₘ Gₘ_c = Σₘ Gᵣₘ Jₘ_c` (the `(r,c)` entries of `J·G` and `G·J` agree). -/
private lemma comm_entry (J G : Mat8) (c : Nat) (hc : c < 8) (hcomm : commM J G = zero8)
    (r : Nat) (hr : r < 8) :
    ent J r 0 * ent G 0 c + ent J r 1 * ent G 1 c + ent J r 2 * ent G 2 c
    + ent J r 3 * ent G 3 c + ent J r 4 * ent G 4 c + ent J r 5 * ent G 5 c
    + ent J r 6 * ent G 6 c + ent J r 7 * ent G 7 c
    = ent G r 0 * ent J 0 c + ent G r 1 * ent J 1 c + ent G r 2 * ent J 2 c
    + ent G r 3 * ent J 3 c + ent G r 4 * ent J 4 c + ent G r 5 * ent J 5 c
    + ent G r 6 * ent J 6 c + ent G r 7 * ent J 7 c := by
  have hz : ent (commM J G) r c = 0 := by rw [hcomm]; exact ent_zero8 r c
  rw [commM, ent_matAdd _ _ r c hr hc, ent_matScale _ _ r c hr hc,
      ent_matMul J G r c hr hc, ent_matMul G J r c hr hc] at hz
  linarith [hz]

/-- **Col-7 bridge** (mirror of `colorInvariant_quarkCol0Zero`): `J e₇ ∈ span{e₀,e₇}`. -/
theorem colorInvariant_quarkCol7Zero (J : Mat8)
    (h1 : commM J H1 = zero8) (h2 : commM J H2 = zero8) :
    ∀ k, k < 8 → k ≠ 0 → k ≠ 7 → ent J k 7 = 0 := by
  -- comm_entry at c=7: LHS uses ent Hᵢ m 7 = 0 (Hᵢ kills e₇), so RHS dot = 0; rows isolate
  have hk1 : ∀ m, m < 8 → ent H1 m 7 = 0 := fun m hm => (color_kills_lepton H1 (by simp [colorGens]) m hm).2
  have hk2 : ∀ m, m < 8 → ent H2 m 7 = 0 := fun m hm => (color_kills_lepton H2 (by simp [colorGens]) m hm).2
  have row : ∀ (Hi : Mat8), commM J Hi = zero8 → (∀ m, m < 8 → ent Hi m 7 = 0) → ∀ r, r < 8 →
      ent Hi r 0 * ent J 0 7 + ent Hi r 1 * ent J 1 7 + ent Hi r 2 * ent J 2 7
      + ent Hi r 3 * ent J 3 7 + ent Hi r 4 * ent J 4 7 + ent Hi r 5 * ent J 5 7
      + ent Hi r 6 * ent J 6 7 + ent Hi r 7 * ent J 7 7 = 0 := by
    intro Hi hc hcol r hr
    have he := comm_entry J Hi 7 (by norm_num) hc r hr
    rw [hcol 0 (by norm_num), hcol 1 (by norm_num), hcol 2 (by norm_num), hcol 3 (by norm_num),
        hcol 4 (by norm_num), hcol 5 (by norm_num), hcol 6 (by norm_num), hcol 7 (by norm_num)] at he
    linarith [he]
  have r6 := row H1 h1 hk1 6 (by norm_num)
  have r1 := row H1 h1 hk1 1 (by norm_num)
  have r5 := row H1 h1 hk1 5 (by norm_num)
  have r2 := row H1 h1 hk1 2 (by norm_num)
  have s4 := row H2 h2 hk2 4 (by norm_num)
  have s3 := row H2 h2 hk2 3 (by norm_num)
  simp only [show ent H1 6 0 = 0 from by decide, show ent H1 6 1 = 1 from by decide,
    show ent H1 6 2 = 0 from by decide, show ent H1 6 3 = 0 from by decide,
    show ent H1 6 4 = 0 from by decide, show ent H1 6 5 = 0 from by decide,
    show ent H1 6 6 = 0 from by decide, show ent H1 6 7 = 0 from by decide] at r6
  simp only [show ent H1 1 0 = 0 from by decide, show ent H1 1 1 = 0 from by decide,
    show ent H1 1 2 = 0 from by decide, show ent H1 1 3 = 0 from by decide,
    show ent H1 1 4 = 0 from by decide, show ent H1 1 5 = 0 from by decide,
    show ent H1 1 6 = -1 from by decide, show ent H1 1 7 = 0 from by decide] at r1
  simp only [show ent H1 5 0 = 0 from by decide, show ent H1 5 1 = 0 from by decide,
    show ent H1 5 2 = -1 from by decide, show ent H1 5 3 = 0 from by decide,
    show ent H1 5 4 = 0 from by decide, show ent H1 5 5 = 0 from by decide,
    show ent H1 5 6 = 0 from by decide, show ent H1 5 7 = 0 from by decide] at r5
  simp only [show ent H1 2 0 = 0 from by decide, show ent H1 2 1 = 0 from by decide,
    show ent H1 2 2 = 0 from by decide, show ent H1 2 3 = 0 from by decide,
    show ent H1 2 4 = 0 from by decide, show ent H1 2 5 = 1 from by decide,
    show ent H1 2 6 = 0 from by decide, show ent H1 2 7 = 0 from by decide] at r2
  simp only [show ent H2 4 0 = 0 from by decide, show ent H2 4 1 = 0 from by decide,
    show ent H2 4 2 = 0 from by decide, show ent H2 4 3 = -1 from by decide,
    show ent H2 4 4 = 0 from by decide, show ent H2 4 5 = 0 from by decide,
    show ent H2 4 6 = 0 from by decide, show ent H2 4 7 = 0 from by decide] at s4
  simp only [show ent H2 3 0 = 0 from by decide, show ent H2 3 1 = 0 from by decide,
    show ent H2 3 2 = 0 from by decide, show ent H2 3 3 = 0 from by decide,
    show ent H2 3 4 = 1 from by decide, show ent H2 3 5 = 0 from by decide,
    show ent H2 3 6 = 0 from by decide, show ent H2 3 7 = 0 from by decide] at s3
  intro k hk hk0 hk7
  interval_cases k <;> first | (exfalso; exact hk0 rfl) | (exfalso; exact hk7 rfl) | linarith

/-- **The column-1 seed.**  For a skew color-invariant `J`, the column `J e₁` lies in mode 1:
    `J e₁ = (J₆₁)·e₆`, i.e. `ent J m 1 = 0` for `m ∉ {6}`.  (`m=1` skew; `m∈{2,3,4,5}` from
    `J e₁ ∈ ker H₂`; `m=0,7` from the col-0/col-7 bridges via skew.) -/
private lemma Jcol1_seed (J : Mat8) (hsk : isSkew J)
    (h1 : commM J H1 = zero8) (h2 : commM J H2 = zero8) :
    ent J 0 1 = 0 ∧ ent J 1 1 = 0 ∧ ent J 2 1 = 0 ∧ ent J 3 1 = 0
    ∧ ent J 4 1 = 0 ∧ ent J 5 1 = 0 ∧ ent J 7 1 = 0 := by
  have c0 := colorInvariant_quarkCol0Zero J h1 h2  -- ent J k 0 = 0 (quark k)
  have c7 := colorInvariant_quarkCol7Zero J h1 h2  -- ent J k 7 = 0 (quark k)
  -- m∈{2,3,4,5}: J e₁ ∈ ker H₂ (H₂ kills e₁: column 1 of H₂ is zero)
  have hH2col1 : ∀ m, m < 8 → ent H2 m 1 = 0 := by decide
  have rowH2 : ∀ r, r < 8 →
      ent H2 r 0 * ent J 0 1 + ent H2 r 1 * ent J 1 1 + ent H2 r 2 * ent J 2 1
      + ent H2 r 3 * ent J 3 1 + ent H2 r 4 * ent J 4 1 + ent H2 r 5 * ent J 5 1
      + ent H2 r 6 * ent J 6 1 + ent H2 r 7 * ent J 7 1 = 0 := by
    intro r hr
    have he := comm_entry J H2 1 (by norm_num) h2 r hr
    rw [hH2col1 0 (by norm_num), hH2col1 1 (by norm_num), hH2col1 2 (by norm_num),
        hH2col1 3 (by norm_num), hH2col1 4 (by norm_num), hH2col1 5 (by norm_num),
        hH2col1 6 (by norm_num), hH2col1 7 (by norm_num)] at he
    linarith [he]
  have q5 := rowH2 5 (by norm_num)
  have q2 := rowH2 2 (by norm_num)
  have q4 := rowH2 4 (by norm_num)
  have q3 := rowH2 3 (by norm_num)
  simp only [show ent H2 5 0 = 0 from by decide, show ent H2 5 1 = 0 from by decide,
    show ent H2 5 2 = 1 from by decide, show ent H2 5 3 = 0 from by decide,
    show ent H2 5 4 = 0 from by decide, show ent H2 5 5 = 0 from by decide,
    show ent H2 5 6 = 0 from by decide, show ent H2 5 7 = 0 from by decide] at q5
  simp only [show ent H2 2 0 = 0 from by decide, show ent H2 2 1 = 0 from by decide,
    show ent H2 2 2 = 0 from by decide, show ent H2 2 3 = 0 from by decide,
    show ent H2 2 4 = 0 from by decide, show ent H2 2 5 = -1 from by decide,
    show ent H2 2 6 = 0 from by decide, show ent H2 2 7 = 0 from by decide] at q2
  simp only [show ent H2 4 0 = 0 from by decide, show ent H2 4 1 = 0 from by decide,
    show ent H2 4 2 = 0 from by decide, show ent H2 4 3 = -1 from by decide,
    show ent H2 4 4 = 0 from by decide, show ent H2 4 5 = 0 from by decide,
    show ent H2 4 6 = 0 from by decide, show ent H2 4 7 = 0 from by decide] at q4
  simp only [show ent H2 3 0 = 0 from by decide, show ent H2 3 1 = 0 from by decide,
    show ent H2 3 2 = 0 from by decide, show ent H2 3 3 = 0 from by decide,
    show ent H2 3 4 = 1 from by decide, show ent H2 3 5 = 0 from by decide,
    show ent H2 3 6 = 0 from by decide, show ent H2 3 7 = 0 from by decide] at q3
  -- assemble.  q5 isolates J₂₁, q4→J₃₁, q3→J₄₁, q2→J₅₁ (others via skew + col-0/7 bridges)
  refine ⟨?_, ?_, by linarith [q5], by linarith [q4], by linarith [q3], by linarith [q2], ?_⟩
  · -- ent J 0 1 = -ent J 1 0 (skew); ent J 1 0 = 0 (col-0)
    have ha := hsk 0 (by norm_num) 1 (by norm_num)
    have hb := c0 1 (by norm_num) (by norm_num) (by norm_num)
    omega
  · have := hsk 1 (by norm_num) 1 (by norm_num); omega   -- ent J 1 1 = -ent J 1 1
  · -- ent J 7 1 = -ent J 1 7 (skew); ent J 1 7 = 0 (col-7)
    have ha := hsk 7 (by norm_num) 1 (by norm_num)
    have hb := c7 1 (by norm_num) (by norm_num) (by norm_num)
    omega

/-- **C2 (quark) — the quark single orbit.**  A skew (reality) complex structure `J`
    commuting with the whole color `𝔰𝔲(3)` acts as `±J_c` on the quark block: there is
    `b = ±1` with `ent J r qc = b · ent J_c r qc` for every quark column `qc ∈ {1,…,6}`.
    Together with the lepton pinning (`colorInvariant_pins_lepton`) this completes the
    statement that the SU(3)-invariant orthogonal complex structures are `±J_c` up to the
    independent lepton/quark orientations — `J_c` representing the distinguished orbit. -/
theorem quark_single_orbit (J : Mat8) (hsk : isSkew J) (hcs : isComplexStructure J)
    (hcomm : ∀ G ∈ colorGens, commM J G = zero8) :
    ∃ b : Int, (b = 1 ∨ b = -1) ∧
      ∀ r qc, r < 8 → qc < 8 → qc ≠ 0 → qc ≠ 7 → ent J r qc = b * ent Jc r qc := by
  have h1 : commM J H1 = zero8 := hcomm H1 (by simp [colorGens])
  have h2 : commM J H2 = zero8 := hcomm H2 (by simp [colorGens])
  have hA13 : commM J A13 = zero8 := hcomm A13 (by simp [colorGens])
  have hB12 : commM J B12 = zero8 := hcomm B12 (by simp [colorGens])
  have hB13 : commM J B13 = zero8 := hcomm B13 (by simp [colorGens])
  have hA12 : commM J A12 = zero8 := hcomm A12 (by simp [colorGens])
  obtain ⟨s0, s1, s2, s3, s4, s5, s7⟩ := Jcol1_seed J hsk h1 h2
  set b := ent J 6 1 with hb
  -- b = ±1 from J² = −I at (1,1)
  have hbpm : b = 1 ∨ b = -1 := by
    have h11 : ent (matMul J J) 1 1 = -1 := by
      show ent (matMul J J) 1 1 = ent negId8 1 1; rw [show matMul J J = negId8 from hcs]
    rw [ent_matMul J J 1 1 (by norm_num) (by norm_num), s0, s1, s2, s3, s4, s5, s7] at h11
    have a16 : ent J 1 6 = - b := by have := hsk 1 (by norm_num) 6 (by norm_num); omega
    rw [a16] at h11
    have hbb : b * b = 1 := by nlinarith [h11]
    rcases mul_self_eq_one_iff.mp hbb with h | h
    · exact Or.inl h
    · exact Or.inr h
  refine ⟨b, hbpm, ?_⟩
  -- the RHS dot Σ Gᵣₘ Jₘ₁ collapses to Gᵣ₆·b (seed: only m=6 survives)
  have prop : ∀ (G : Mat8) (r : Nat),
      ent G r 0 * ent J 0 1 + ent G r 1 * ent J 1 1 + ent G r 2 * ent J 2 1
      + ent G r 3 * ent J 3 1 + ent G r 4 * ent J 4 1 + ent G r 5 * ent J 5 1
      + ent G r 6 * ent J 6 1 + ent G r 7 * ent J 7 1 = ent G r 6 * b := by
    intro G r; rw [s0, s1, s2, s3, s4, s5, s7]; ring
  intro r qc hr hqc hqc0 hqc7
  interval_cases qc
  · exact absurd rfl hqc0                                    -- qc = 0
  · -- qc = 1 : the seed column
    have hJ : ent J r 1 = (if r = 6 then b else 0) := by interval_cases r <;> simp_all
    have hJc : ent Jc r 1 = (if r = 6 then (1:Int) else 0) := by interval_cases r <;> decide
    rw [hJ, hJc]; split <;> ring
  · -- qc = 2 via A12 (A12 e₁ = e₂)
    have he := comm_entry J A12 1 (by norm_num) hA12 r hr
    rw [prop A12 r] at he
    simp only [show ent A12 0 1 = 0 from by decide, show ent A12 1 1 = 0 from by decide,
      show ent A12 2 1 = 1 from by decide, show ent A12 3 1 = 0 from by decide,
      show ent A12 4 1 = 0 from by decide, show ent A12 5 1 = 0 from by decide,
      show ent A12 6 1 = 0 from by decide, show ent A12 7 1 = 0 from by decide,
      mul_zero, mul_one, zero_add, add_zero] at he
    have hjc : ent Jc r 2 = ent A12 r 6 := by interval_cases r <;> decide
    rw [hjc]; linarith [he]
  · -- qc = 3 via A13
    have he := comm_entry J A13 1 (by norm_num) hA13 r hr
    rw [prop A13 r] at he
    simp only [show ent A13 0 1 = 0 from by decide, show ent A13 1 1 = 0 from by decide,
      show ent A13 2 1 = 0 from by decide, show ent A13 3 1 = 1 from by decide,
      show ent A13 4 1 = 0 from by decide, show ent A13 5 1 = 0 from by decide,
      show ent A13 6 1 = 0 from by decide, show ent A13 7 1 = 0 from by decide,
      mul_zero, mul_one, zero_add, add_zero] at he
    have hjc : ent Jc r 3 = ent A13 r 6 := by interval_cases r <;> decide
    rw [hjc]; linarith [he]
  · -- qc = 4 via B13
    have he := comm_entry J B13 1 (by norm_num) hB13 r hr
    rw [prop B13 r] at he
    simp only [show ent B13 0 1 = 0 from by decide, show ent B13 1 1 = 0 from by decide,
      show ent B13 2 1 = 0 from by decide, show ent B13 3 1 = 0 from by decide,
      show ent B13 4 1 = 1 from by decide, show ent B13 5 1 = 0 from by decide,
      show ent B13 6 1 = 0 from by decide, show ent B13 7 1 = 0 from by decide,
      mul_zero, mul_one, zero_add, add_zero] at he
    have hjc : ent Jc r 4 = ent B13 r 6 := by interval_cases r <;> decide
    rw [hjc]; linarith [he]
  · -- qc = 5 via B12
    have he := comm_entry J B12 1 (by norm_num) hB12 r hr
    rw [prop B12 r] at he
    simp only [show ent B12 0 1 = 0 from by decide, show ent B12 1 1 = 0 from by decide,
      show ent B12 2 1 = 0 from by decide, show ent B12 3 1 = 0 from by decide,
      show ent B12 4 1 = 0 from by decide, show ent B12 5 1 = 1 from by decide,
      show ent B12 6 1 = 0 from by decide, show ent B12 7 1 = 0 from by decide,
      mul_zero, mul_one, zero_add, add_zero] at he
    have hjc : ent Jc r 5 = ent B12 r 6 := by interval_cases r <;> decide
    rw [hjc]; linarith [he]
  · -- qc = 6 via H1 (H1 e₁ = e₆)
    have he := comm_entry J H1 1 (by norm_num) h1 r hr
    rw [prop H1 r] at he
    simp only [show ent H1 0 1 = 0 from by decide, show ent H1 1 1 = 0 from by decide,
      show ent H1 2 1 = 0 from by decide, show ent H1 3 1 = 0 from by decide,
      show ent H1 4 1 = 0 from by decide, show ent H1 5 1 = 0 from by decide,
      show ent H1 6 1 = 1 from by decide, show ent H1 7 1 = 0 from by decide,
      mul_zero, mul_one, zero_add, add_zero] at he
    have hjc : ent Jc r 6 = ent H1 r 6 := by interval_cases r <;> decide
    rw [hjc]; linarith [he]
  · exact absurd rfl hqc7                                    -- qc = 7

/-- **C2 (complete) — classification of the color-invariant complex structures.**  A skew
    (reality) complex structure commuting with the whole color `𝔰𝔲(3)` is determined by two
    signs: it acts as `±[[0,1],[−1,0]]` on the lepton singlet `{0,7}` and as `±J_c` on the
    quark block.  So the SU(3)-invariant orthogonal complex structures are exactly the four
    `(±_lepton, ±_quark)·J_c`, with `J_c = γ₀γ₅` the distinguished representative — the
    single Furey complex structure that makes `ℝ⁸ = ℂ⁴ = 1 ⊕ 3`.  This is the full C2
    pinning (lepton sector + quark single orbit). -/
theorem colorInvariant_classification (J : Mat8)
    (hsk : isSkew J) (hcs : isComplexStructure J) (hcomm : ∀ G ∈ colorGens, commM J G = zero8) :
    (ent J 0 7 = 1 ∨ ent J 0 7 = -1)
    ∧ (∃ b : Int, (b = 1 ∨ b = -1) ∧
        ∀ r qc, r < 8 → qc < 8 → qc ≠ 0 → qc ≠ 7 → ent J r qc = b * ent Jc r qc) :=
  ⟨colorInvariant_pins_lepton J hsk hcs
      (hcomm H1 (by simp [colorGens])) (hcomm H2 (by simp [colorGens])),
   quark_single_orbit J hsk hcs hcomm⟩

/-! ## The color Z₃ and the lepton-phase holonomy obstruction

The center/Weyl `Z₃` of the color `SU(3)` cyclically permutes the three quark modes and
**fixes the lepton singlet line `{0,7}` pointwise**.  Hence any color-`Z₃` holonomy on the
lepton is trivial: the Brannen lepton phase `φ = Q/3` cannot be a single-ideal color
holonomy — it must come from the (external) generation `Z₃`.  See
`koide_phase_law/witt_map.py`. -/

/-- The color `Z₃`: cycles the quark modes `(1,6)→(2,5)→(3,4)→(1,6)`, fixes the lepton line. -/
def colorZ3 : Mat8 :=
  mkMat [(0,0,1), (7,7,1), (2,1,1), (3,2,1), (1,3,1), (5,6,1), (4,5,1), (6,4,1)]

/-- `colorZ3` is order 3. -/
theorem colorZ3_cubed : matMul (matMul colorZ3 colorZ3) colorZ3 = id8 := by decide

/-- `colorZ3` commutes with the color complex structure `J_c` (it is a genuine color rotation). -/
theorem colorZ3_commutes_Jc : commM colorZ3 Jc = zero8 := by decide

/-- **The holonomy obstruction.**  `colorZ3` acts as the identity on the lepton singlet line
    `{0,7}` (column 0 is `e₀`, column 7 is `e₇`).  So transport of the lepton around the color
    `Z₃` cycle accumulates *no* phase — the lepton phase `φ=2/9` is not a color holonomy. -/
theorem colorZ3_fixes_lepton :
    (∀ r, r < 8 → ent colorZ3 r 0 = (if r = 0 then 1 else 0))
    ∧ (∀ r, r < 8 → ent colorZ3 r 7 = (if r = 7 then 1 else 0)) := by decide

end SCPv59.Furey7D.ColorSU3
