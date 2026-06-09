# v63 — Proceeding on the Loop-Ratio Route (first step + honest frontier)

**Date**: 2026-05-28
**Status**: the gating checks closed `c_eff`, holonomy, and canonical localization;
v62 left `cos(3φ) = −c₃/(4c₆)` (a dynamical loop-ratio) as the only number-type
that can equal the transcendental `cos(2/3)`. This is the first concrete step on
that route, plus an honest scope of what a real derivation requires.

---

## The target

Critical-point equation: `cos(3φ) = −c₃/(4c₆)`. For the observed `φ = 2/9`,
`cos(3·2/9) = cos(2/3) = 0.78589`. So a parameter-free loop must yield

```
    c₃/c₆ = −4 cos(2/3) = −3.1436.
```

## Step 1 — the parameter-free 1-loop ratio FAILS (`loop_ratio.py`)

Computed `c₃, c₆` from the 1-loop Coleman–Weinberg `V_eff(φ)` of the lepton
Brannen spectrum at the structural `t²=1/2` (the `m·log m → 0` integrand is finite
— my earlier "singular" worry was a `0·log0` numerical artifact):

| loop form | c₃ | c₆ | −c₃/(4c₆) | in [−1,1]? |
|---|---|---|---|---|
| `Σ m(log m−3/2)` | 1.980 | 0.0643 | **−7.70** | no |
| `Σ m²(log m²−3/2)` | 29.01 | −0.190 | **+38.07** | no |
| `Σ s⁴(log s⁴−3/2)` | 29.01 | −0.190 | **+38.07** | no |

Every form gives `−c₃/(4c₆)` **far outside `[−1,1]`** (target `0.786`) and even
**sign-unstable** across forms. So at 1-loop there is no in-range nontrivial
critical point — only the commensurate `sin(3φ)=0` minima (`φ=0, π/3, …`), which
are π-rational, hence `≠ 2/9`. (Already proven in v59 `PhaseExclusions`:
`z3_potential_does_not_select_2_9`, `cos6_potential_does_not_select_2_9`.)

## Step 2 — what an offset requires (and why it's a dial now)

To land on `cos(2/3)` the `cos(6φ)` coefficient must be boosted to
`|c₆| ≈ |c₃|/3.14` — comparable to `c₃`. But 1-loop **suppresses** the higher
harmonic (`c₆ ≪ c₃`). Boosting it needs a 2-loop / sub-leading term whose
coefficient is, right now, a **free fit** (the empirical "tilt" `ε` in
`v59/algebra/v_eff_loop.py`). A free coefficient = a dial, not a prediction.

## The make-or-break frontier (v59 Step-6)

v59's "loop" is **not** the continuum CW used above — it is a **finite sum over
the G₂-equivariant subspace** of `Cl(7)_even`. A finite, structured sum can have
different harmonic content (`c₆` need not be suppressed), so its ratio could be
fixed by the G₂ representation with no free loop scale. The single question that
decides everything:

> Does the **G₂-equivariant finite-sum loop** give `c₃/c₆ = −4cos(2/3)` with **no
> free coefficient** (only the structural counting `N_X`)?

- **Yes** → the phase is a genuine parameter-free prediction (the win).
- **Needs a fitted `ε` / loop order / `α`** → the phase is a dial, or α-downstream
  (like `f_g ~ α^{21/2}`), reducing to the irreducible input `α`.

**To compute it (unbuilt):** (1) build the `G₂ ⊂ Spin(7)` action on the sector
ambient (`Λ⁴(ℝ⁷)` for d-quark; the `L⊕F` structure for u-quark), identify the
`N_X`-dim invariant subspace; (2) evaluate the finite-sum `V_eff(φ)` on it;
(3) extract `c₃, c₆` and test the ratio with zero free coefficients.

**Honest risk**: the 1-loop form-sensitivity (`−7.70` vs `+38.1`) warns that the
ratio may stay scheme-dependent even in the finite sum. This is a genuine frontier
with uncertain payoff — not a sure thing.

---

## Verdict

**Yes, there is a way to proceed — and it is the *only* way left — but the simplest
version doesn't work, and the real calculation is hard with a real dial risk.**

- The phase definitively lives in the loop-ratio (every other route is closed).
- The parameter-free 1-loop ratio does **not** reproduce `2/9` (commensurate only).
- The decisive, well-posed next computation is the **G₂-equivariant finite-sum
  loop**: does `c₃/c₆ = −4cos(2/3)` come out structurally, or only by fitting?
  That single result is the make-or-break for the whole program's claim that the
  phase is derived rather than input.

**Artifacts**: `loop_ratio.py` (6/6 checks). Builds on `v62/residual_audit/
phase_dynamical_home.py`, `v59/algebra/v_eff_loop.py`, `v59/.../PhaseExclusions.lean`.
