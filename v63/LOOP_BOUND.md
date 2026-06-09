# v63 — What Bounds the Phase Parameter (loop calc)

**Date**: 2026-05-28
**Question**: in the surviving loop-ratio picture `cos(3φ) = −c₃/(4c₆)`, what
bounds the free parameter — and is it α-determined? (The one thread with upside:
if the phase collapses into `α`, the input count drops.)

---

## The free parameter and its two bounds

The phase is the vacuum alignment of `V_eff = c₃cos(3φ) + c₆cos(6φ)`; the free
parameter is `r = c₃/c₆`, with `cos(3φ) = −r/4`.

**(B1) Hard / kinematic bound.** A real in-range minimum needs `cos(3φ) ∈ [−1,1]`,
so

```
    r = c₃/c₆ ∈ [−4, 4]   ⇔   φ ∈ [0, π/3]   (the Z₃ fundamental wedge).
```

This is absolute — the only constraint the Z₃ potential structure imposes. All
three observed phases sit inside, near the commensurate `r=−4` (φ=0) end:

| sector | φ | `r = −4cos(3φ)` |
|---|---|---|
| lepton | +0.2222 | −3.144 |
| d-quark | +0.1086 | −3.790 |
| u-quark | −0.0725 | −3.906 |

**(B2) Loop / α bound.** If the phase is a 1-loop *gauge* correction to the
commensurate tree value `φ=0` (v59 `brannen_phase_alpha.py`), its size is bounded
by the loop parameter: `|φ| ≲ N·α`, with `N` a structural counting `≤ dim G₂ = 14`.
So the natural ceiling is `14·α(M_Z) ≈ 0.109`.

---

## Where the phases fall — the decisive split

| sector | \|φ\| | `N·α` fit | gap | ≤ 14·α(M_Z)? |
|---|---|---|---|---|
| d-quark | 0.1086 | `14·α(M_Z)` = 0.1094 | 0.76 % | **yes** |
| u-quark | 0.0725 | `10·α(0)` = 0.0730 | 0.64 % | **yes** |
| **lepton** | **0.2222** | needs `N ≈ 28.4` | — | **NO — exceeds by ×2.03** |

- **Quark phases** sit *inside* the loop ceiling and match `N·α` (N = 14 = dim G₂,
  10 = dim Spin(5)) to **<1 %**. So the quark phase parameter **is α-determined** —
  it collapses into `α` + a structural integer (α-downstream, not an independent
  input).
- **The lepton phase** `2/9 = 0.2222` **pierces the 1-loop ceiling 14·α ≈ 0.109 by
  2×**, and matching it as `N·α` needs `N ≈ 28–30` — no clean structural integer
  (cf. the unambiguous `14 = dim G₂`). So it is **not** a 1-loop α-correction.

---

## The dichotomy (and the honest answer)

- **lepton**: `φ_l = 2/9 = Q/3` *exactly* (rational, Koide-tied, ≈10⁻⁵). Not
  `N·α`, and not a tree Z₃ minimum (v59 `PhaseExclusions`: the tree potential gives
  commensurate `φ=0`). A clean **structural rational** with no derivation — the
  **irreducible free input** among the phases.
- **quarks**: `φ_d ≈ 14·α(M_Z)`, `φ_u ≈ −10·α(0)` — transcendental, 1-loop
  α-corrections off `φ=0`, **bounded by ≲14·α** and **α-determined**.

> **What the parameter is bounded by:**
> - absolutely — the **Z₃ wedge** `φ ∈ [0, π/3]` (`r ∈ [−4,4]`);
> - in the loop picture — **`|φ| ≲ (max structural N)·α ≈ 14α ≈ 0.11`**.
>
> The loop bound **contains and α-determines the quark phases** (they reduce to
> `α`), but the **lepton `2/9` pierces it 2×** — so it stays a genuine free
> (structural-rational) input. The hoped-for "phase collapses into α" is **real
> for the quarks, false for the lepton.**

---

## Verification

`loop_bound.py` (7/7 checks): all phases in the Z₃ wedge; `r ∈ (−4,4)`; quark
phases below `14·α(M_Z)` and matching `N·α` to <1%; lepton `> 14·α(M_Z)` needing
`N ≈ 28`; `φ_l = Q/3` exactly.

`lean/LoopBound.lean` (no `sorry`): `lepton_in_wedge` (`0 < 2/9 < π/3` via
`pi_gt_three`); `d_within_loop_bound` + `d_matches_14alpha` (<1%);
`lepton_exceeds_loop_bound` (`2/9 > 14·α(M_Z)`); `lepton_needs_large_N`
(`(2/9)/α(M_Z) > 20`); `phase_bound_dichotomy`.

---

## Net (the upside, quantified)

The loop calc delivers a **partial** win: the two **quark** Brannen phases are
α-bounded 1-loop corrections (`|φ|≲14α`) and collapse into the single input `α` +
structural integers — they are *not* independent conjectures. But the **lepton
phase `2/9 = Q/3`** is too large for a `N≤14` 1-loop correction, isn't a tree-level
minimum, and is a clean Koide-tied rational: it remains the **one irreducible free
input** in the phase sector. So the residual ledger tightens (quark phases →
α-downstream) without closing the lepton phase.

**Artifacts**: `loop_bound.py` (7/7), `lean/LoopBound.lean` (no `sorry`). Builds on
`brannen_phase_alpha.py`, `loop_ratio.py`, `g2_loop_precheck.py`,
`v59/.../PhaseExclusions.lean`.
