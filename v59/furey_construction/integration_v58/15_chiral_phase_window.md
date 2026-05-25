# M4 chiral-protection for the phase φ=2/9 — pursued & Lean-formalized (a constraint, not a derivation)

*2026-05-24.  Per the steer "pursue the chiral protection / M4 lead for the phase, carry the
derivation into Lean."  Result: the M4 **light-electron window** is formalized and machine-checked
(`lean/ChiralPhaseWindow.lean`, axiom-clean) — the electron-massless point is `φ=π/12`, and the
physical `φ=2/9` sits just inside it (electron light-but-massive).  **But this is a *constraint*
(φ in the window), NOT a derivation of `2/9`** — confirming the note's honest verdict.*

## What M4 says, made precise

In the Koide-normalised Brannen amplitude `√mₙ/a = 1 + √2·cos(φ + 2πn/3)` (`t²=½`, the proven
amplitude), the phase `φ` controls the *lightest* mass:

- **Chiral (electron-massless) point:** the electron amplitude vanishes when its angle hits `3π/4`:
  `1 + √2·cos(φ+2π/3) = 0 ⟺ cos(φ+2π/3) = −1/√2 ⟺ φ = π/12`.  (`cos(3π/4)=−√2/2`.)
  So if an exact chiral symmetry forced `m_e=0`, it would sit `φ` at **π/12**.
- **Physical phase, just inside the window:** `φ = 2/9 < π/12` (⟺ `π > 8/3`), so the electron is
  **light but massive**: `1 + √2·cos(2/9+2π/3) > 0` (≈ 0.040 > 0).  The gap `π/12 − 2/9 ≈ 0.040`
  is tiny — "chiral protection keeps the electron light."

(Linear check: `√m_e/a ≈ π/12 − φ` near the chiral point — the electron √-mass *is* the angular
gap, `≈0.040`, matching PDG `√m_e/a = 0.0404`.)

## Lean (machine-checked, axiom-clean, in `AxiomCheck`)

`lean/ChiralPhaseWindow.lean`:
- `chiral_massless_point` : `1 + √2·cos(π/12 + 2π/3) = 0` — the electron-massless point is `π/12`.
- `physical_phase_below_chiral` : `2/9 < π/12` — physical phase inside the window.
- `electron_light_but_massive` : `0 < 1 + √2·cos(2/9 + 2π/3)` — light but nonzero (via strict
  monotonicity of `cos` on `[0,π]`).
- `m4_chiral_window` : the three bundled.

## Honest verdict — constraint, not derivation

This **pursues M4 to its honest limit and locks the geometry in Lean**, but it does **not** derive
`φ = 2/9`:

- Chiral protection explains why `φ` lies *near* `π/12` (small `m_e`) and pins the **window**; it
  does not select the point inside it.
- The specific value `φ = 2/9 = Q/3` is the **separate, unexplained Brannen phase law**; the offset
  `π/12 − 2/9` is fixed by the *measured* `m_e` (an **input**), not derived.
- `π/12` is **π-rational** (geometric); `2/9` is **not** (`PhaseExclusions.koide_not_pi_rational`).
  So `2/9` is genuinely shifted off the chiral angle by a non-geometric amount — exactly the residual.

So M4 is a **consistency relation**: it ties `φ` to whether `m_e` is small/zero (the hierarchy), and
confirms `2/9` lives in the light-electron window — but `m_e` (hence the offset, hence `2/9`) is
input, so the phase is not derived.  This matches the broader picture: the phase `φ=2/9` is a free
residual (not symmetry-fixed, not π-rational, not a holonomy), now also seen as *near the chiral
edge* — a place, not a derivation.

## Status of the lepton sector after this (honest roll-up)

- **Amplitude `t²=½`** (Koide `Q=2/3`): grounded as the **L-grade complex reversion-norm** (`OctoHalf`,
  `14_…md`), modulo the half-norm/grade-balance.  Proven-anchored to `mass∈L`.
- **Phase `φ=2/9`**: a **free residual** — near the chiral/electron-massless edge `π/12` (this note),
  but its value `2/9=Q/3` is unexplained, non-π-rational, input-dependent.  *Not derived.*

The chiral window is the most that M4 gives: a Lean-locked constraint, not the phase.
