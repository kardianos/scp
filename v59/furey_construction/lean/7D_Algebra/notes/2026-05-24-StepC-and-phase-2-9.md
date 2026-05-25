# Next: (C) pinning the unique J, and the φ = 2/9 phase derivation

**Date**: 2026-05-24
**Status**: (C1) and (C2) DONE incl. the quark single-orbit (`ColorSU3.lean`); only the φ=2/9
phase remains.  Prerequisites (A) and (B) DONE (`LeptonRealityForcing.lean`): no symmetric
matrix is a complex structure (the whole `{Λ⁰}⊕F` subspace excluded, not just blades),
`L = skew = so(8)`, `{Λ⁰}⊕F = symmetric`, reality (Hermitian mass) ⇒ `J` skew ⇒ `J ∈ L`.

## DONE — (C1) + (C2)  (`ColorSU3.lean`, builds, decide / Mathlib, axiom-clean)

**Correction found:** the *color* complex structure is `J_c = γ₀γ₅`, NOT the generation
ℍ-slice `γ₀γ₁` (which does not preserve the color-singlet sector).  `J_c` is block-diagonal,
pairs `(0↔7),(1↔6),(2↔5),(3↔4)`, so `ℝ⁸ = ℂ⁴ = 1 ⊕ 3` — lepton singlet `{0,7}` ⊕ three quark
modes (the Furey minimal-left-ideal split); `J_c ∈ Λ² ⊂ L`.

- **(C1)** Explicit `𝔰𝔲(3)`: 8 generators (Cartan `H₁,H₂`; roots `A_ab,B_ab` on the 3 quark
  modes).  Proved (`decide`, 0 axioms): all skew (`color_isSkew`); all annihilate the lepton
  singlet `{0,7}` (`color_kills_lepton` — lepton IS a color singlet); all commute with `J_c`
  (`color_commutes_Jc` — `J_c` color-invariant); **full A₂ closure** — all 28 brackets, std
  structure constants `∈{±1,±2}` (`colorSU3_closes`).
- **(C2)** Pinning.  `lepton_block_canonical`: a skew complex structure with column-0
  supported on `{0,7}` has `J₀₇ = ±1` (via `(J·J)₀₀ = −1` ⇒ `−(J₀₇)² = −1`).  The bridge
  `colorInvariant_quarkCol0Zero`: commuting with `H₁,H₂` ⇒ `J e₀ ∈ ker H₁ ∩ ker H₂ =
  span{e₀,e₇}` (Cartan common kernel = lepton sector).  Hence `colorInvariant_pins_lepton`:
  **any skew (reality) complex structure commuting with the color Cartan has lepton block
  `±[[0,1],[−1,0]]`** — the lepton complex structure is unique up to orientation, realized by
  `J_c` (`Jc_pinned`).

### DONE — the quark single-orbit (`ColorSU3.quark_single_orbit`, `colorInvariant_classification`)
A skew color-invariant complex structure acts as `±J_c` on the quark block too:
`∃ b=±1, ∀ r, ∀ qc∈{1..6}, J r qc = b·(J_c) r qc`.  Schur made concrete (no abstract rep
theory): **seed** `J e₁ = b·J_c e₁` (`Jcol1_seed`: `m=1` skew; `m∈{2,3,4,5}` from
`J e₁ ∈ ker H₂`, since `H₂ e₁=0`; `m=0,7` from the col-0/col-7 bridges via skew) — needed the
**col-7 bridge** `colorInvariant_quarkCol7Zero`; **`b=±1`** from `J²=−I` at `(1,1)`;
**propagation** — the 5 generators send `e₁ ↦ e₂,e₃,e₄,e₅,e₆` (`A₁₂,A₁₃,B₁₃,B₁₂,H₁`), so
`[J,G]=0` carries `J=b·J_c` across all six quark columns (one `comm_entry` collapse each).
Capstone `colorInvariant_classification`: the SU(3)-invariant orthogonal CS's are the four
`(±_lepton,±_quark)·J_c`, with `J_c` the distinguished representative making `ℝ⁸=ℂ⁴=1⊕3`.
Builds, axiom-clean (decide / std trio), 0 sorry / 0 native_decide.  **(C) is now complete
except the φ phase below.**

So far the forcing is: **the lepton's complex structure `J` lies in `L = Λ²⊕Λ⁶`** (a
28-dim subspace).  Two things remain to make the picture complete.

## (C) Pin `J` to a unique element (the SU(3)-commutant)

`J ∈ L` is a 28-dim constraint; physically `J` is one specific operator (Furey's `i`,
realized as `e₀₁ = γ₀γ₁`).  The constraints that should cut `L` down to it:

1. **Color SU(3)-equivariance** `[J, 𝔰𝔲(3)] = 0`.  The color algebra `𝔰𝔲(3) ⊂ 𝔤₂ ⊂ 𝔰𝔬(7)`
   acts on the spinor `8 = 1 ⊕ 1 ⊕ 3 ⊕ 3̄` (lepton singlets {0,7}, d-triplet {1,2,3},
   u-antitriplet {4,5,6}).  Its real commutant is `ℂ ⊕ ℂ ⊕ M₂(ℝ)` (Schur: a complex
   structure on each of `3`, `3̄`, and a `2×2` real block on the singlet pair).  A complex
   structure in the commutant therefore acts as a scalar `i` on each color irrep and as
   `[[0,1],[−1,0]]` on {0,7} — exactly the L-block the agent found in `LeptonGradeForcing`.
2. **Charge / Z₃ normalization**: the U(1) charges `(0,+1,−1,0)` on `(1,3,3̄,1)` and the
   generation Z₃ phase fix `J` up to color conjugacy → the orbit of `e₀₁`.

**Formalization targets** (hardest piece — genuinely new representation theory):
- (C1) Construct the 8 generators of `𝔰𝔲(3) ⊂ 𝔤₂` explicitly in the gamma model (as
  bivectors / L-blades) and `decide` `[gen, color] = 0` for the candidate `J = e₀₁` and the
  per-irrep blocks.  Concretely: which Λ²-bivectors commute with the color Cartan, and show
  the color-neutral complex structures form (over the singlet pair + triplets) the expected
  `ℂ⊕ℂ⊕M₂(ℝ)`.
- (C2) Show the SU(3)-invariant *orthogonal* complex structures in `L` are a single SU(3)
  orbit, represented by `e₀₁`.  (Likely needs a small explicit basis-level argument rather
  than `decide`; the commutant is 8-dim, the invariant complex structures a low-dim subset.)
- Bridge note: `orthogonal_complexStructure_skew` (B) currently lives over Mathlib `Matrix`;
  (C) would benefit from a `toMatrix : Mat8 → Matrix (Fin 8) (Fin 8) ℚ` ring-hom bridge so
  the gamma blades, reality, and the commutant all live in one model.  That bridge
  (`toMatrix (matMul A B) = toMatrix A * toMatrix B`, etc.) is itself a worthwhile,
  reusable target.

## The φ = 2/9 phase derivation (the remaining lepton-sector input)

The Brannen mass kernel is `M = a(I + ξ S + ξ̄ S²)`, `ξ = e^{iφ}` (`KernelEigenvalues.lean`,
`M_mulVec_eigen`/`lam_eq_brannen`).  Two independent real parameters fit the charged
leptons: the **amplitude** `t = |ξ|` and the **phase** `φ`.

- **Amplitude** `t² = 1/2` ⇔ Koide `Q = 2/3` — already structural
  (`BrannenKernel.koide_iff_constraint`; `Q = dim G₂ / dim Spin(7) = 14/21`).
- **Phase** `φ = 2/9`.  What is *currently* proved is only the **arithmetic** identity
  `φ = Q/3 = (14/21)/3 = 2/9` (`LieDimensions` φ-rational, `SpinDimension.brannen_phase_structural`,
  `ScaleBridge.sin_sq_thW_eq_brannen_phase`).  This *assumes* the relation `φ = Q/3`; it does
  **not derive** why the mass-operator phase equals one third of the Koide ratio.  That
  `3·φ = Q` (the factor of 3 = the three generations / `Z₃`) is the open dynamical content.

**Why this connects to (C).**  The phase `φ` is the *argument* of `ξ`, and `ξ` rotates in the
complex direction set by `J`.  Once `J = e₀₁ ∈ Λ²` is pinned (C), the phase is a rotation
*in the `J`-plane*, and the `Z₃` generation shift `S` acts in the same minimal left ideal.
The conjecture to chase: **the compatibility of the single complex structure `J` with the
`Z₃` cyclic generation structure and the `Q = 2/3` amplitude forces `φ = Q/3`** — i.e. the
phase is not a second free input but is fixed by where `J` sits relative to `S`.  If so, the
last empirical lepton parameter (after the overall scale `a`) becomes structural.

**Formalization targets:**
- (P1) State the phase–amplitude relation `3·φ = Q` as a *theorem about the kernel*, not a
  definitional rational identity: e.g. that the unique Hermitian `M` with eigenvalues in
  Koide ratio `2/3` AND compatible with the `J`-plane / `Z₃` action has phase `φ = 2/9`.
- (P2) Tie `φ` to the `J`-rotation explicitly: express `ξ S` in terms of `J` and `S` and show
  the allowed phase is quantized by the `Z₃` (`ω = e^{2πi/3}`) already in `CyclicShift.lean`.
- (P3) If (P1)/(P2) resist, document the precise gap (as here): the value `2/9` is arithmetic
  given `φ=Q/3`; the *physical* `φ=Q/3` is the open relation.

## Suggested order
(C1) → (C2) → (P2) → (P1).  (C1) is `decide`-friendly (commutators of explicit blades);
(C2) and the phase derivation are the genuinely new content and the natural place for the
next deep push.  The `toMatrix` bridge is an enabling side-quest that pays off in both.
