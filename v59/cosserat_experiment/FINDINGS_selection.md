# Selection Rule — Z₂ × Z₂ Decomposition of Cl(7)_even

**Date**: 2026-05-22
**Parent**: `FINDINGS_single_source.md`
**Status**: Novel structural observation, NOT a Lagrangian derivation.

## Honest framing

The user's skepticism — "I'm skeptical you can find this selection rule" — was
correct.  I tried 6+ hypotheses and none yielded a clean derivation from first
principles.  But I did find a **non-trivial structural observation** in
trying.  This document records that observation honestly.

## The novel finding: L ⊕ F structural decomposition

The Cl(7)_even ≅ Cl(6) ≅ ℂ⊗𝕆 parent algebra decomposes into TWO
structurally distinct pieces (apart from the trivial identity Λ⁰):

```
  L  =  Λ²ℝ⁷ ⊕ Λ⁶ℝ⁷       dim 21 + 7 = 28      "Lie algebra content"
  F  =  Λ⁴ℝ⁷             dim 35              "G₂-form content"
```

These two pieces are **canonically distinguished** by the G₂ structure:

- **L contains NO G₂-invariant** — neither Λ² nor Λ⁶ has a G₂-singlet.
  Decomposes as `Λ² = 14 + 7` (= G₂ adjoint + fundamental) and `Λ⁶ = 7`
  (= G₂ fundamental).  All G₂-content is in NON-singlet reps.

- **F contains the G₂-defining structure** — the coassociative 4-form `*φ`
  lives in Λ⁴.  Decomposes as `Λ⁴ = 1 + 7 + 27` (G₂ singlet + fund + 27).
  The "1" is the G₂-defining tensor.

This binary distinction is the canonical G₂-structural decomposition of
Cl(7)_even — every Λ⁰⊕Λ²⊕Λ⁴⊕Λ⁶ element belongs to either L (no G₂-form)
or F (has G₂-form content).

## The observed Z₂ × Z₂ pattern

For the three Furey fermion sectors, the ambient subspace of Cl(7)_even
chosen by each sector corresponds to a 2-bit selection of (L, F):

| Sector | Furey N | Bit-L | Bit-F | Ambient | Dim |
|---|---|---|---|---|---|
| e_R (lepton) | 0 | 1 | 0 | L only | **28** |
| d_R (d-quark) | 1 | 0 | 1 | F only | **35** |
| u_R (u-quark) | 2 | 1 | 1 | L ⊕ F | **63** |
| ν_R (neutrino) | 3 | ? | ? | ? | ? |

For N = 0, 1, 2 the pattern is `(1,0), (0,1), (1,1)` — a Z₂ × Z₂ structure
with the trivial `(0,0)` slot empty (or possibly the neutrino, which has
non-Brannen mass spectrum).

The **additive identity D_u = D_e + D_d = 28 + 35 = 63** emerges naturally
because the u-quark ambient is the direct sum L ⊕ F.

## What L and F represent physically

**L (Λ² ⊕ Λ⁶, dim 28)** — "Spin(7) gauge / Lie algebra content":
- Λ²ℝ⁷ = so(7) Lie algebra of Spin(7) — the 21 rotation generators
- Λ⁶ℝ⁷ ≅ ℝ⁷ via Hodge dual — the 7-dim S⁷ = Spin(7)/G₂ direction
- These are the "structureless" (from G₂-perspective) parts of Cl(7)_even
- A fermion in a *Spin(7) representation* couples to L

**F (Λ⁴, dim 35)** — "G₂-form / octonion content":
- Λ⁴ℝ⁷ contains the coassociative 4-form *φ (1-dim G₂-invariant)
- Plus 34 dims of G₂-charged content (7 + 27 under G₂ branching)
- The G₂-form *φ ENCODES the octonion multiplication structure
- A fermion in a *G₂/color representation* couples to F

So the **Z₂ × Z₂ pattern reflects two independent gauge-content questions**:
- **Bit-L = "does this fermion couple to Spin(7) rotation generators?"**
- **Bit-F = "does this fermion couple to the octonion (color) algebra structure?"**

The pattern across sectors:
- **Leptons** (color singlets): YES gauge, NO color → L only
- **d-quarks** (single color creation): NO gauge?, YES color → F only
- **u-quarks** (double color creation): YES gauge, YES color → L ⊕ F

The "NO gauge" for d-quarks is the puzzling case — why don't they couple
to L?  This is the part the mechanism still doesn't explain cleanly.

## Lean encoding (planned)

The L ⊕ F decomposition is a real structural identity that can be added
to `Predictions.lean` as theorems:

```lean
def L_content : ℕ := cl7_grade_lambda2 + cl7_grade_lambda6    -- 28
def F_content : ℕ := cl7_grade_lambda4                         -- 35

theorem L_dim : L_content = 28 := by decide
theorem F_dim : F_content = 35 := by decide
theorem L_plus_F : L_content + F_content = dimU63 := by decide  -- 63

-- The Z₂ × Z₂ pattern (empirical observation)
-- (Couldn't derive, but encoded as a structural identity)
```

These would all be axiom-free (pure arithmetic).

## What this is and is not

**WHAT IT IS:**
- A CANONICAL STRUCTURAL DECOMPOSITION of Cl(7)_even based on G₂-invariant
  content.  L = "no G₂-invariant content", F = "G₂-invariant-bearing content".
- A clean explanation of the additive identity D_u = D_e + D_d.
- A Z₂ × Z₂ classification of fermion sectors by their L/F coupling.
- A physical interpretation: L = gauge content, F = color/octonion content.

**WHAT IT IS NOT:**
- A derivation of WHY each Furey N takes its specific (Bit-L, Bit-F) value.
  The pattern (1,0), (0,1), (1,1) for N=0, 1, 2 is observed but the
  underlying mechanism is not pinned to a Lagrangian.
- An explanation of why d-quark (N=1) DOESN'T couple to L.  Naively, color
  triplets should still feel Spin(7) rotations — so why is Bit-L = 0 for d?
- A prediction for the N=3 (neutrino) case.

## Novel relationships uncovered

1. **L vs F = no-G₂-invariant vs G₂-invariant-bearing**: this is a structural
   bisection of Cl(7)_even based on the G₂ defining tensor.  Not a Lie algebra
   decomposition (L is not closed under bracket) but a graded vector-space
   decomposition that maps onto the v59 sector ambient structure.

2. **The bits track the G₂ defining-tensor coupling**: F contains the G₂-form,
   L doesn't.  Fermions that "see" the G₂-form (carry octonion multiplication
   structure = color charge) couple to F.  This is the qualitative explanation
   for the empirical pattern.

3. **Additive identity is direct sum**: D_u = D_e + D_d isn't coincidental —
   the u-quark sector is the *direct sum* of the lepton and d-quark sectors
   in the L ⊕ F decomposition.

4. **The Brannen Koide formula `Q = (1 + 2t²)/3` with `t² = 1 − 14/D` becomes**
   particularly meaningful: the 14 = dim G₂ is the G₂ orbit dim, and the D_N
   is the dim of the L+F subspace each sector couples to.  Both numerator
   AND denominator are now graded by G₂ structure.

## The mechanism question (still open)

The cleanest open question:

> Given a fermion type X with Furey N-number N_X, what determines the
> bit-vector (B_L^X, B_F^X)?

For N = 0, 1, 2 we observe (1,0), (0,1), (1,1).  No formula `N → (B_L, B_F)`
is obvious.  Candidates:

(a) **Lagrangian-level**: a specific Yukawa Lagrangian with sector-dependent
    projection operators.  The projectors might come from natural
    Spin(8)/Spin(7)/G₂ symmetry breaking.

(b) **Fock-space level**: in the Cl(6) Witt construction, the operators
    that couple |Ω_N⟩ to itself within the 3-generation triality space
    might split into L-type and F-type contributions sector-by-sector.

(c) **Topological**: a Chern-class-like invariant of a v59-natural bundle
    might take different values in {L, F, L⊕F, ∅} for different sectors.

I haven't been able to derive any of these in this session.  The user's
skepticism was well-founded — this is a multi-week (or longer) project,
not a session-level derivation.

## What I'd want to do next (if continuing)

1. **Write the explicit Brannen Yukawa Lagrangian** for the Furey ℂ⊗ℍ⊗𝕆
   algebra with sector-dependent projectors.  See if the L ⊕ F structure
   emerges naturally from the Yukawa form.

2. **Compute Spin(8) → Spin(7) → G₂ symmetry breaking** explicitly.  At each
   step, identify the "broken generators" and see whether they map to L vs F.

3. **Look at the WZW / Wess-Zumino term** structure on the Spin(8) coset.
   Different sectors might have different WZW windings.

4. **Test against neutrinos**: even if Brannen ansatz doesn't apply, maybe
   neutrino oscillation mass-squared differences can distinguish the N=3
   (B_L, B_F) candidates.

But none of these is doable in a short session.  The user was right.

## Honest closing

What I CAN deliver in this session:
- ✓ The single source is Cl(7)_even ≅ Cl(6) (verified, Lean-encoded).
- ✓ The L ⊕ F decomposition is canonical (this finding, novel observation).
- ✓ The Z₂ × Z₂ pattern across sectors is documented.
- ✗ The mechanism that selects (Bit-L, Bit-F) for each Furey N is NOT derived.

This is honest progress, not a derivation.  The structural picture is
cleaner than before but the final mechanism question remains open.
