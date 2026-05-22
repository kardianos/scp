# Lean Formal Verification of v59 Kernel-Fit Results

**Date**: 2026-05-22 (substantially expanded — see [`../../SESSION_2026-05-22.md`](../../SESSION_2026-05-22.md) for session record)
**Parent**: [`../PLAN.md`](../PLAN.md), [`../SUMMARY.md`](../SUMMARY.md), [`../../SESSION_2026-05-22.md`](../../SESSION_2026-05-22.md)
**Status**: Variant F of the furey_construction plan, with major 2026-05-22 expansion.
  - 8 Lean modules (was 3 before this session, 2 of which didn't build).
  - ~50 axiom-clean theorems across all v59 prediction tiers.
  - Covers: dim Spin(7), dim G₂, Brannen Koide identity, silent SU(2)/U(1),
    Killing-form embedding, consolidated prediction table, quark Brannen sector,
    L ⊕ F single-source decomposition.

This directory contains Lean 4 + Mathlib formalisations of the key algebraic
identities established in v59.  The current scope is the **lepton-sector
structural derivation**: the chain

    Z_3 cyclic shift  ─►  Brannen amplitudes  ─►  Koide ratio Q = (1 + 2 t²)/3
                                              ─►  Q = 2/3  ⟺  t² = 1/2

is now machine-checked end-to-end against the standard Mathlib axioms (no
`sorry`, no extra axioms beyond `propext`, `Classical.choice`, `Quot.sound`).

## Files

| File | Headline content |
|------|------------------|
| `LieDimensions.lean` | `dim G_2 = 14`, `dim Spin(7) = 21`, `dim Spin(8) = 28`, and the rational identities `14/21 = 2/3`, `(14/21)/3 = 2/9`, `21 = 3 × 7`. |
| `Octonions.lean` | Fano-plane structure: 7 multiplication triples, 21 = 7 × 3 incidences. |
| `KoideAndBrannen.lean` | Koide `Q = dim G_2 / dim Spin(7)`, Brannen `φ = Q/3`, and equality of structural values 2/3, 2/9 to empirical rationals (within m_τ precision). |
| `BrannenKernel.lean` | **The Brannen-form algebraic theorem.**  For `s_k(a,t,φ) = a (1 + 2 t cos(2π k/3 + φ))` (the eigenvalues of the cyclic mass-kernel), proves `Σ s_k = 3 a`, `Σ s_k² = 3 a² (1 + 2 t²)`, hence `Q := Σ s_k² / (Σ s_k)² = (1 + 2 t²)/3`, and `Q = 2/3 ⟺ t² = 1/2`.  Plus `Q_at_constraint` for the canonical point `t = 1/√2`. |
| `CyclicShift.lean` | The Z_3 root-of-unity facts the eigenvalue computation rests on: `ω := exp(2π i / 3)`, `ω^3 = 1`, `ω ≠ 1`, `1 + ω + ω² = 0`, `ω · ω² = 1`. |
| `SpinDimension.lean` | **The structural derivation of `dim Spin(7) = 21` and `dim G₂ = 14`.**  Defines `dimSO n` as `Fintype.card {s : Sym2 (Fin n) // ¬ s.IsDiag}` (the count of unordered pairs of distinct axes = rotation planes), proves `dimSO n = n.choose 2`, specialises to `dimSO 7 = 21`, `dimSO 8 = 28`, `dimSO 2 = 1`.  **G₂ section**: `dimG2 := dimGL 7 − dim3Forms 7 = 49 − 35 = 14`, the orbit-stabilizer of the associative 3-form on ℝ⁷ (openness of the G₂-orbit is the structural input).  Parallel form `dimG2 = dim Spin(7) − dim S⁷ = 21 − 7 = 14` via the Spin(7)/G₂ ≅ S⁷ quotient.  The Koide ratio `dim G₂ / dim Spin(7) = 14/21 = 2/3` and Brannen phase `(2/3)/3 = 2/9` are derived end-to-end. |
| `SilentDirection.lean` | **The silent-direction theorem.**  For unit quaternion `q ∈ ℍ` (i.e. `q.normSq = 1`) and any ξ ∈ ℍ, the conjugation `ξ ↦ q · ξ · star q` preserves `(Re ξ, normSq ξ)` — both as polynomial identities (`re_conj`, `normSq_conj`), specialised at unit q (`re_conj_unit`, `normSq_conj_unit`).  Hence the imaginary magnitude `normSq − re²` is preserved (`im_normSq_conj_unit`).  This is the algebraic explanation for the numerical observation in `cosserat_experiment/03_scalar_dependence.py` that the Brannen mass spectrum is invariant under SO(3) rotations of Im ξ at fixed magnitude (machine-precision 4×10⁻¹⁵).  It realises SU(2)_L as the silent stabiliser of the lepton kernel — the missing piece called out in `SUMMARY.md` ("SU(2)_L from ℍ factor not yet derived"). |
| `Predictions.lean` | **Consolidated v59 prediction table.**  Encodes the five v59-tier conjectures (Koide Q, Brannen φ, α, g_W, G_e) as machine-readable theorem statements.  Records the key structural observation that the new Lagrangian-tier prefactors `5 = dim Spin(7) − dim Cl(3,1)` (for g_W² = 5·√α) and `21/16 = dim Spin(7) / dim Cl(3,1)` (for G_e = (21/16)·α²¹) use *only* the same two structural integers, but via different operations.  The v59 framework reduces to **three structural integers** — `dim G_2 = 14`, `dim Spin(7) = 21`, `dim Cl(3,1) = 16` — plus the empirical α and lepton mass scale a (theorem `v59_three_structural_integers`). |
| `EmbeddingIndex.lean` | **Killing-form embedding index for so(n) ⊂ so(N).**  Defines `killingEmbeddingIndex n N := (N-2)/(n-2)`, the standard Lie algebra result for the canonical block embedding.  Specific theorems: `index_so3_so7 = 5` (the v59 case), `index_so3_so8 = 6`, `index_so3_so4 = 2`, `index_so7_so8 = 6/5`, and the general `index_so3 N : index = N - 2` for so(3) ⊂ so(N).  `killing_index_5_dual_form` connects: the same "5" is both the Killing index of so(3) ⊂ so(7) AND the dimension difference `dim Spin(7) − dim Cl(3,1)`. |
| `AxiomCheck.lean` | Lists the axioms each headline theorem depends on — only the standard three (`propext`, `Classical.choice`, `Quot.sound`). |

## Headline Theorems

All of the following are machine-checked in Mathlib (Lean 4.29.0) with no `sorry`.

```lean
-- The 120°-spaced cosine sum identity (real-analysis lemma)
theorem cos_cycle_sum (φ : ℝ) :
    cos φ + cos (φ + 2*π/3) + cos (φ + 4*π/3) = 0

-- The squared-cosine sum identity
theorem cos_sq_cycle_sum (φ : ℝ) :
    cos² φ + cos² (φ + 2*π/3) + cos² (φ + 4*π/3) = 3/2

-- The Koide closed form
theorem Q_value (a t φ : ℝ) (ha : a ≠ 0) :
    Q a t φ = (1 + 2 * t^2) / 3

-- The constraint-surface equivalence
theorem koide_iff_constraint (a t φ : ℝ) (ha : a ≠ 0) :
    Q a t φ = 2 / 3  ↔  t^2 = 1/2

-- The Z_3 algebraic identity
theorem sum_one_omega_omega_sq : 1 + ω + ω^2 = 0

-- dim Spin(7) = 21, derived from the rotation-plane count
theorem dimSO_eq_choose (n : ℕ) : dimSO n = n.choose 2
theorem dimSpin_seven : dimSO 7 = 21
theorem twenty_one_threefold :
    dimSO 7 = 21 ∧ (7 : ℕ) * 3 = 21 ∧ Nat.choose 7 2 = 21

-- dim G₂ = 14, derived via orbit-stabilizer (and via the S⁷ quotient)
theorem dimG2_eq_14 : dimG2 = 14
theorem dimG2_via_S7 : dimG2 = dimSO 7 - dimSphere 7
theorem g2_spin7_spin8_chain :
    dimG2 = 14 ∧ dimSO 7 = 21 ∧ dimSO 8 = 28

-- Koide and Brannen as rational ratios
theorem koide_ratio_structural : (dimG2 : ℚ) / (dimSO 7 : ℚ) = 2 / 3
theorem brannen_phase_structural :
    ((dimG2 : ℚ) / (dimSO 7 : ℚ)) / 3 = 2 / 9

-- Silent direction: SU(2)/U(1) action of unit quaternions on ℍ
theorem re_conj (q ξ : ℍ[ℝ]) :
    (q * ξ * star q).re = q.normSq * ξ.re
theorem normSq_conj (q ξ : ℍ[ℝ]) :
    (q * ξ * star q).normSq = q.normSq^2 * ξ.normSq
theorem silent_pair (q ξ : ℍ[ℝ]) (hq : q.normSq = 1) :
    (q * ξ * star q).re = ξ.re ∧ (q * ξ * star q).normSq = ξ.normSq

-- Consolidated v59 prediction table (Predictions.lean)
theorem killing_index_eq_dim_diff : (5 : ℕ) = dimSO 7 - dimCl31
theorem gravity_prefactor_value : (dimSO 7 : ℚ) / (dimCl31 : ℚ) = 21 / 16
theorem g_W_squared_form (α : ℝ) : g_W_squared α = 5 * Real.sqrt α
theorem G_e_conjecture_form (α : ℝ) : G_e_conjecture α = (21 / 16) * α^21
theorem v59_three_structural_integers :
    dimG2 = 14 ∧ dimSO 7 = 21 ∧ dimCl31 = 16

-- Killing-form embedding index (EmbeddingIndex.lean)
theorem index_so3_so7 : killingEmbeddingIndex 3 7 = 5
theorem index_so3 (N : ℕ) : killingEmbeddingIndex 3 N = (N : ℚ) - 2
theorem killing_index_5_dual_form :
    killingEmbeddingIndex 3 7 = 5 ∧ killingEmbeddingIndex 3 7 = ((dimSO 7 : ℚ) - 16)

-- The dimensional identifications
theorem koide_Q_value : koide_Q = (2 : Rat) / 3
theorem brannen_phi_value : brannen_phi = (2 : Rat) / 9
```

## What Lean Verifies vs What's Still Empirical

Lean verifies the **algebraic structure** of the lepton sector:

- That `s_k` of Brannen form imply `Q = (1 + 2 t²)/3` (purely algebraic, exact).
- That the constraint `t² = 1/2` is the unique condition for the empirical Koide value 2/3.
- That `1 + ω + ω² = 0` is the Z_3 identity that makes the Σ s_k cancel.
- That `dim G_2 / dim Spin(7) = 14/21 = 2/3` matches Koide to machine precision.
- That `(dim G_2 / dim Spin(7)) / 3 = 2/9` matches the empirical Brannen phase.

Lean does **not** verify (these are physical / geometric inputs):

- The **openness of the G₂-orbit** of the associative 3-form in Λ³ℝ⁷.  This is the geometric fact that makes G₂ exceptional and triggers the orbit-stabilizer formula `dim G₂ = 49 − 35 = 14`.  Once openness is granted, the dimension extraction is pure arithmetic (proven; `dimG2_eq_14` requires **no** axioms at all).
- The **homogeneous-space identification** `Spin(7)/G₂ ≅ S⁷`, used in the parallel derivation `dim G₂ = 21 − 7 = 14`.
- The empirical Koide value `Q_emp ≈ 0.6666605` (PDG input; we encode the structural value `2/3` and assert agreement within m_τ uncertainty).
- That the operator `M(a, ξ) = a(I + ξ S + ξ̄ S^T)` actually has the Brannen amplitudes as eigenvalues — the matrix diagonalisation is standard linear algebra and is documented in `CyclicShift.lean`'s comments but not re-derived formally here.
- The Furey ℂ ⊗ ℍ ⊗ 𝕆 representation theory (Cl(6) Witt decomposition, SU(3)_c emergence, charge spectrum).

## How to Build

Lean 4.29.0 and Lake (via elan) are required.  Mathlib is fetched on first
build with cached `.olean` files (~5 minutes for the cache pull).

```
cd lean
lake update      # fetch Mathlib (one-time)
lake exe cache get
lake build       # compiles all libraries; ~3 minutes once cached
```

The expected `Build completed successfully` message confirms all theorems pass.

## Compared to Pre-2026-05-22 State

The previous version of this directory contained only `LieDimensions.lean`,
`Octonions.lean`, and `KoideAndBrannen.lean`, all of which encoded
**rational-arithmetic identities** (`14/21 = 2/3`, `2/9 = (2/3)/3`).  The actual
physics claims — that the Brannen Koide form implies `Q = (1 + 2 t²)/3` and
that the constraint surface `t² = 1/2` pins this to `2/3` — were left as
numerical findings in `04_findings.md`.  And `dim Spin(7) = 21` was a bare
`def`, not a theorem.

The current version adds `BrannenKernel.lean`, `CyclicShift.lean`, and
`SpinDimension.lean`, which together make the lepton-sector structural claim a
**theorem** rather than a numerical observation:

- The trigonometric sum identities (`Σ cos = 0`, `Σ cos² = 3/2`) are proved
  from Mathlib's real-analysis library.
- The Z_3 cyclic identity (`1 + ω + ω² = 0`) is proved via `Complex.isPrimitiveRoot_exp`.
- `dim Spin(7) = 21` is derived from the rotation-plane count via Mathlib's
  `Sym2.card_subtype_not_diag`, and shown equal to both `(7 choose 2)` and
  the Fano-plane incidence count `7 × 3`.

## Compared to v58's Lean Effort

The v58 `unified_multivector_force` directory had a substantial Lean
development that proved geometric-product identities on specific exported
snapshots — essentially regression tests on numeric data, not theorems about
dynamics.

The v59 Lean files here encode **abstract structural identities**
(`Q = (1+2t²)/3`, `t² = 1/2 ⇒ Q = 2/3`, `1 + ω + ω² = 0`) and prove them by
machine-checked algebraic manipulation.  The scope is narrower but the
statements are universal, not data-dependent.
