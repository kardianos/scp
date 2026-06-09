# G1 — the rank tension: characterized, with a deflation

**Date**: 2026-05-25 (first concrete attack on G1, the second v59 blocker)
**Artifacts**: `01_rank_tension.py` (runs clean), `../../lean/RankTension.lean`
(compiles clean against v59 Mathlib, axiom-clean), this doc.
**Builds on**: v59 `gaps/ew_scale_bridge/` (README + FINDINGS — the tension was
stated there; here it is computed and sharpened).

## The problem (v59's framing)

The EW bridge `v_Higgs = dim(L)²·a_ℓ² = 784 a²` (0.068%) needs one bilinear
`Y ∈ End(L) = M₂₈(ℝ)` to be **full-rank democratic** (`‖Y‖²_F = 784a²`). The
physical lepton spectrum is **rank-3** (3 generations), with
`Σm = Tr(M†M) = 9Qa² = 6a²` (the gravity charge). One matrix cannot be both
(784 ≠ 6; full-rank ≠ rank-3), and no single-step `so(8)→H` leaves "3 light +
25 eaten." This is the rank tension.

## What this session establishes (all computed / machine-checked)

### 1. A deflation: the "6/784 gravity↔EW bonus" is not independent evidence
v59 advertises `Σm = (6/784)·v` at 0.07% as a relation tying the gravity charge to
the Higgs vev. But `Σm = 9Qa²` is **definitional** (`a = (Σ√m)/3`, `Q` = Koide), so

```
(Σm / v_obs) / (9Q/784)  =  784 a² / v_obs  =  the EW bridge match (0.068%).
```

The "bonus" carries **no information beyond `v = 784a²`** — it is the bridge match
re-expressed. (Directly analogous to v59's own `g_W² = 5√α = α(M_Z)`-in-disguise
deflation in `RIGOR_AUDIT`.) Machine-checked: `RankTension.deflation`. The numerator
`6 = N_gen²·Q` is structural (`ratio_numerator`, `grav_ew_ratio`).

### 2. The single-`Y` rank tension is real (not a bookkeeping artifact)
- A **democratic** vacuum (random ±a, "all 784 components comparable") has
  `‖Y‖²_F ≈ 784a²` ✓ but **28 comparable singular values** — full rank, no
  hierarchy. Wrong spectrum.
- The physical **Brannen kernel** is rank-3, hierarchical (`√m_k`), `‖M‖²_F = 6a²`.
- A **hand-built** `Y` (3 Brannen + 25 heavy at ≈5.6a) reproduces *both* norms
  (`‖Y‖²_F = 784a²` and a rank-3 light spectrum) — but it is **hierarchical, not
  the democratic vacuum** an `so(8)`-invariant potential selects, and its 25 heavy
  directions have **no symmetry home**. So such a `Y` is *posited*, not derived.

### 3. The 25-direction / subalgebra obstruction
`dim so(8) = 28`; "25 eaten Goldstones" ⇒ an unbroken `H` of `dim 3`. The maximal
proper subalgebra of `so(8)` is `dim Spin(7) = 21 < 25`, so this is **not a
single-step breaking** (`RankTension.subalgebra_obstruction`); and the 3 unbroken
directions of any such `H` are **not** the 3 generations — those come from the `Z₃`
**triality** on a separate generation space, not an `su(2)` stabilizer of `L`.

## Verdict

The only consistent reading is **two objects on different spaces**:
- an EW condensate/operator on `End(L) = M₂₈` (dim 784; sets `v` via `‖·‖²_F`), and
- the rank-3 Brannen kernel on the 3-dim **generation** space (from `Z₃` triality).

This **dissolves the contradiction** but does **not derive `v` from `a_ℓ`**: the two
share only the scale `a_ℓ` and the conjecture *"physical scale = Frobenius² of the
relevant mass bilinear"* (R1), which still has **no dynamical home** (v59's actual
Higgs `XiVacuum` is the 4-dim `ℍ` Mexican hat, not the 28-dim bilinear). So:

> **G1's "1 input `a_ℓ`" headline stays CONDITIONAL on R1+R2.** The rank tension is
> real for any single `Y`; the two-object reading is the only consistent resolution;
> and `v_Higgs` remains rigorously a **second scale** until a Lagrangian produces
> *both* objects from the one scale.

## Status table

| claim | status |
|---|---|
| bridge `v=784a²` (0.068%) | **[emp]** reproduced |
| `784 = dim End(L)` (Burnside) | **[thm]** (v59) |
| `Σm = 9Qa²`, `9 = N_gen²` | **[thm]** definitional (`nine_is_Ngen_sq`) |
| `6/784 = N_gen²Q/dim(L)²` | **[thm]** (`grav_ew_ratio`) |
| "6/784 gravity bonus" = bridge match re-expressed | **[thm] deflation** (`deflation`) |
| single `Y` can't be full-rank-784 **and** rank-3 | **[shown]** (§2) |
| no single-step `so(8)→H` eats 25 | **[thm-arith]** (`subalgebra_obstruction`) |
| R1 (`v=‖Y‖²_F`) + R2 (component `=a_ℓ`) | **[conj]** no dynamical home |
| two-object resolution | **consistent**, but leaves `v` a conditional 2nd scale |
| Lagrangian producing both objects from one scale | **[open]** — the real remaining task |

## What a resolution needs (carried from v59 FINDINGS, sharpened)
A dynamical principle that produces, from one scale `a_ℓ`: (i) the full-rank EW
condensate on `End(L)` with `‖·‖²_F = 784a²`, and (ii) the rank-3 Brannen kernel on
the triality generation space — with the 25 ambient directions made heavy. The
honest obstacle is the same one G9 hit on the dynamics side: there is no derived
Lagrangian yet, only the algebraic skeleton + posited scale conjecture. Building
that Lagrangian (and giving R1 a home, e.g. the deferred `ℝ⁴ ↪ End(L)` generation
embedding) is the next concrete step.
