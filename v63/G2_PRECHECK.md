# v63 — G₂-Loop Pre-Check (the build is NOT warranted)

**Date**: 2026-05-28
**Status**: cheap gating check before committing to the full G₂-equivariant
finite-sum loop (v59 Step-6) — the one route that could make the phase loop-ratio
`cos(3φ)=−c₃/(4c₆)` parameter-free. **Result: negative on both decisive questions.
Do not build it.**

---

## The two questions

The hope (`LOOP_RATIO.md`): the v59 "loop" is a **finite sum over a G₂-equivariant
subspace** of Cl(7), which might have un-suppressed `c₆` (unlike continuum CW),
making `c₃/c₆` structurally fixed and possibly `= −4cos(2/3)`.

- **(Q1)** Does the rep-theory premise hold — is there a clean `N_X`-dim
  (specifically the claimed `14`) G₂ subspace of `Λ⁴(ℝ⁷)`?
- **(Q2)** Even if so, can `N_X` change the ratio `c₃/c₆`?

## The computation (`g2_loop_precheck.py`, 6/6 checks)

Built `g₂ = Stab_{so(7)}(φ)` (the derivation algebra) from the repo's own Fano
table — `dim g₂ = 14` ✓ — and decomposed `Λᵖ(ℝ⁷)` by the quadratic Casimir
`Σₐ ρₚ(Tₐ)²`:

| space | Casimir spectrum (value : multiplicity) | irreps |
|---|---|---|
| `Λ¹` (7) | `{2.0 : 7}` | **7** |
| `Λ²` (21) | `{2.0 : 7, 4.0 : 14}` | **7 + 14** |
| `Λ⁴` (35) | `{0.0 : 1, 2.0 : 7, 4.667 : 27}` | **1 + 7 + 27** |

(`Casimir(7) = 2` matches the repo's noted `C₂(7)=2`; the `7`'s value is
consistent across `Λ¹,Λ²,Λ⁴` — the machinery is correct.)

## Findings — both negative

**(Q1) The premise is misattributed.** The adjoint **`14` lives in `Λ²`, not
`Λ⁴`**. `Λ⁴(ℝ⁷) = 1 + 7 + 27` (one singlet = the coassociative 4-form; no 14).
So there is no "14-dim G₂-equivariant subspace of `Λ⁴`" for the d-quark loop to
sum over — v59 `v_eff_loop.py`'s `N_X = 14 from Λ⁴` is a rep-theory error.

**(Q2) `N_X` cannot help anyway.** v59's loop is
`V_eff(φ) = −(N_X/64π²)·Tr[M(φ)²·log(…)]`: `N_X` is an **overall, φ-independent
prefactor**. Scaling `V` by a constant scales `c₃` and `c₆` *equally*, so it
**cancels** in `−c₃/(4c₆)` (verified: ratio `−7.70` for `N_X=1` and `N_X=14`,
identical). No G₂ counting — 14, 27, 7, anything — can un-suppress `c₆`.

**Why, structurally:** the `c₆` suppression is a **generation-Z₃ harmonic** effect
(`cos 6φ` enters at higher order in `t` than `cos 3φ`). The G₂/Λ⁴ structure is the
**internal/color** index — orthogonal to the generation Z₃. The two don't talk, so
the G₂ projection can't change the Z₃-harmonic content.

---

## Verdict

The finite-G₂-sum escape hatch is **illusory**: the claimed subspace isn't there,
and even a correct `N_X` cancels in the ratio. So the `c₆` boost required for
`φ=2/9` (`|c₆| ≈ |c₃|/3.14`) must come from a **genuinely free 2-loop coefficient**
(the empirical "tilt" `ε`) — a **dial**, not fixed by the algebra.

**Consequence for the program**: the phase-as-loop-ratio is, on current structure,
a **dial** — its value is set by a free sub-leading loop coefficient. It is not a
parameter-free prediction. The honest residual classification stands: the Brannen
phase is either (a) a free input, or (b) α-downstream *if* that tilt coefficient
turns out to be α-determined — which is itself a separate, unproven claim. The
large G₂-equivariant build would not change this, so it is **not warranted**.

This closes the last constructive route examined in the v62→v63 arc: every path to
a *parameter-free* `2/3` (algebraic no-go, c_eff, holonomy, canonical localization,
density-weighting, and now the G₂ loop counting) is closed. The phase's
transcendentality is inseparable from a free continuous parameter — consistent
throughout.

**Artifacts**: `g2_loop_precheck.py` (6/6 checks). Builds on `loop_ratio.py`,
`v59/algebra/v_eff_loop.py`, the repo Fano table (`SevenDAlgebra.lean`).
