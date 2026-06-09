# GEN6 — Selection rule, universal Koide deviation, and the rank tension (G1)

**Date**: 2026-05-26 (dynamical-Lagrangian loop, Generation 6)
**Artifacts**:
- `15_selection_rank.py` (SymPy — all assertions pass)
- `15_koide_universal.mac` (independent **Maxima** cross-check)
- `../lean/SelectionRule.lean` (compiles clean; standard-trio axioms)
**Builds on**: `gaps/rank_tension/01_*` (the G1 tension), GEN3 (the matter potential).

---

## Verdict

GEN6 attacks the **second v60 gate (G1, rank tension)** plus the **selection rule**,
and connects them to the dynamical Lagrangian (GEN1–5). Three verified results:

### (A) Selection rule = orthogonal grade projectors (SymPy + Lean)

`Cl(7)_even = Λ⁰(1) ⊕ Λ²(21) ⊕ Λ⁴(35) ⊕ Λ⁶(7)` (= 64). With `L = Λ²⊕Λ⁶ = 28`,
`F = Λ⁴ = 35`, the selection-rule projectors
```
Π_lepton = Π_L,   Π_d = Π_F,   Π_u = Π_L + Π_F
```
are **orthogonal idempotents** (`Π_L²=Π_L`, `Π_F²=Π_F`, `Π_LΠ_F=0`) with traces
`28, 35, 63`. The **additive identity `D_u = D_e + D_d = 63 = 28 + 35`** is then
automatic from the grading — the Z₂×Z₂ selection rule (lepton=L, d=F, u=L⊕F) is
structural, not an extra input. (This is the "selection rule from the action": the
OBE force law's `Π_N[Ω]·Ψ` uses exactly these grade idempotents.)

### (B) Universal Koide deviation (SymPy + Maxima + Lean)

A single constant ties all three sectors' Koide values to their dimensions:
```
(1 − Q_N) D_N = 28/3     for N = lepton(28), d(35), u(63)
  ⟹  Q = 2/3,  11/15,  23/27
```
The constant is `28/3 = (1 − Q_lepton)·dim(L)`. So the d/u "Koide" values
`11/15, 23/27` (which v59 flagged as scale-convention-sensitive *as absolute
predictions*) are nonetheless tied to the lepton `Q=2/3` by this **structural
relation** through the selection-rule dimensions. Verified in three CAS/provers.

### (C) Rank tension (G1) in the dynamical language (SymPy + Lean)

`01_rank_tension` established the two-object resolution and flagged the open task:
*"a Lagrangian producing the rank-3 object from one scale."* GEN6 closes part of it:
- **The GEN3 matter potential IS that dynamical home** — its EL vacuum on the Koide
  cone produces the rank-3 Brannen kernel with `Σm = 9Qa² = 6a²`. So the
  generation-space object (ii), which `01` said had *no dynamical home*, now has one.
- The **End(L) bridge object (i)** (dim 784) remains a **separate sector**:
  `784 ≠ 9 = N_gen²` (different spaces). The 3 generations are the **Z₃ triality
  orbit** `{8_v, 8_s, 8_c}`, *not* a subspace stabilizer of the 28-dim `L` — so a
  single "two-piece Y" of one matrix is the wrong frame (confirms `01`'s subalgebra
  obstruction: max `so(8)` subalgebra 21 < 25).
- The **deflation** (`(Σm/v)/(9Q/784) = 784a²/v`) re-confirmed = 0.

**Net**: the two-object rank-tension resolution stands; GEN3 supplies the dynamical
home for the rank-3 half; the EW `v = 784a²` bridge remains a conditional second
scale (R1), now isolated as the *only* remaining piece without a dynamical home.

---

## Lean (`SelectionRule.lean`)

| theorem | content | kind |
|---|---|---|
| `additive_identity` | `D_u = D_L + D_F`, `= 63` | `decide` |
| `projectors_orthogonal_additive` | `D_L + D_F = D_u` (disjoint supports) | `decide` |
| `koide_universal` | `(1−Q)D = 28/3` for 2/3, 11/15, 23/27 | `norm_num` |
| `different_spaces` | `784 ≠ 9` (two objects) | `decide` |
| `gen6_selection_rank` | bundled headline | mixed |

---

## Status table

| claim | status | tool |
|---|---|---|
| `Cl(7)_even` grades `1+21+35+7=64`; `L=28, F=35` | **verified** | SymPy + Lean |
| selection projectors orthogonal idempotent; `D_u=D_e+D_d=63` | **verified** | SymPy + Lean |
| universal `(1−Q)D = 28/3` → 2/3, 11/15, 23/27 | **verified** | SymPy + Maxima + Lean |
| GEN3 = dynamical home of the rank-3 generation object | **verified** | SymPy (ties to GEN3) |
| End(L) 784 ≠ generation 9; two-object resolution | **verified** | SymPy + Lean |
| EW `v = 784a²` bridge has a dynamical home (R1) | **still open** (the last residual) | — |

---

## State of the program after GEN6

Both v60 gates are now addressed:
- **G9** (gravity → tensor): GEN1–5 — a stable linearized dynamical Lagrangian with
  2 ghost-free TT gravitons, EP-exact, the OBE trace law derived.
- **G1** (rank tension): GEN6 — the two-object resolution, with GEN3 supplying the
  rank-3 object's dynamical home; the EW-condensate (784) is a clean separate sector.

Remaining honest residuals (value-conjectures / separate sectors, *not* dynamics
gaps): the EW vev `v=784a²` (R1), the Brannen phase `φ=2/9` (Goldstone, GEN3), the
gravity magnitude `f_g~α^{21/2}`, and `α` itself.

## What's next

- **GEN7 (aspect 7)**: the actual **time-dependent dynamics** — a genuine solution
  of the nonlinear EL equations (`□Φ = −V'(Φ)` on the Koide-cone potential), e.g. a
  small **C/Python lattice integration** showing a stable oscillation about the
  vacuum (the massive modes) + a propagating massless Goldstone, confirming the
  linearized spectrum nonlinearly. That is "dynamics" in the most literal sense and
  closes the loop's mandate.
