# Constrained-OFE extremum: NEGATIVE — the lepton is ℍ-confined, so the F-channel is zero

*2026-05-24.  Set up the constrained-OFE extremum to test whether the dynamics fixes `t²=½`.
**Result: it does not** — and the careful setup *corrects* the optimism of `18_…md`.  The lepton
order parameter is the **color singlet**, confined to the 4-dim `ℍ`-slice, where the grade-4 (F)
channel — the only thing Cl(7) added over Cl(3,0) — is **identically zero**.  So the OFE on the
lepton is back to the `{0,2}=ℍ` dynamics (= v58's Cl(3,0)), whose alignment is flat; the `½` is the
Mexican-hat parameter (input), the phase a Goldstone (free).*

## The decisive fact: the lepton lives in `ℍ`, where F ≡ 0

The lepton order parameter (the Brannen quaternion `ξ`, the XiVacuum field) lives in the **`ℍ`-slice
`{1, e₀₁, e₀₂, e₁₂}`** — the color singlet.  This slice is **closed under multiplication (≅ℍ)**, so
it stays in grades `{0,2}`.  Verified: for any `ℍ`-element, `‖⟨ξ²⟩₄‖² = 0.000000` (the grade-4 / F
channel vanishes identically).  Contrast: a generic 28-dim `L`-config gives `‖⟨M²⟩₄‖² = 5.71` — but
**that is not the lepton.**

So the F-channel non-flatness of `18_…md` requires spreading over the **full** `L=so(8)` (28-dim);
the lepton uses only the **4-dim `ℍ`-slice**, where F is structurally zero.

> **⚠ This corrects `18_…md`.**  Lifting v58 to Cl(7) does **not** supply a working
> alignment-selecting channel *for the lepton*: on the `ℍ`-slice the even dynamics is grades
> `{0,2} = ℍ`, which is *exactly* v58's Cl(3,0) even subalgebra.  No new channel; same flatness.

## What the OFE on `ℍ` actually gives (= the XiVacuum)

On the `ℍ`-slice the OFE reduces to the v59 `XiVacuum` Mexican-hat dynamics (`XiVacuum.lean`):
- the **radial mode** fixes `|ξ|² = ½` — but the `½` is the potential parameter `μ²/λ`, an **input**
  (the hat is *built* with minimum at `½`; the value is not derived);
- the **three imaginary directions are Goldstones** (flat) — the **phase is free**.

So the dynamics yields exactly the *structure* (a fixed-norm vacuum manifold + flat phase) and
**not** the *value*: `t²=½` is the hat parameter (input), `φ` is a Goldstone (residual).

## Verdict on the constrained-OFE line

> **The constrained OFE does NOT fix `t²=½`.**  The lepton is `ℍ`-confined (color singlet); the
> grade-4 (F) channel — the sole Cl(7) novelty — is identically zero there; the residual `{0,2}=ℍ`
> dynamics is the XiVacuum, which delivers `|ξ|²=½` only as the Mexican-hat parameter (input) and the
> phase as a flat Goldstone (free).  The `18_…md` F-channel mechanism applies to colored/spread
> configs, **not** the color-singlet lepton.

This is the same boundary every route has hit, now closed for the dynamics too: **the octonionic
dynamics supplies vacuum STRUCTURE, not the coupling MAGNITUDE** (`17_…md`'s kinematics/dynamics
divide).  Even the F-channel hope is excluded for the lepton by its `ℍ`-confinement.

## The one technically-open crack (honest)

The F-channel *would* act if the lepton order parameter had support **outside `ℍ`** (in the colored
`L`-directions).  It doesn't, because it's the color singlet — *that* is why F vanishes.  So the only
way the OFE/F-channel could fix the lepton couplings is if there were a **dynamical reason the
singlet acquires colored (non-`ℍ`) components** — which would contradict it being the color singlet.
So this crack is, in effect, closed: **`ℍ`-confinement (color-singlet-ness) ⇒ F=0 ⇒ no OFE
constraint on `t²`.**

## Net

The constrained-OFE calculation returns a clean **negative**, and corrects `18`: the lepton's
color-singlet (`ℍ`) confinement makes the grade-4 channel vanish, so the Cl(7) lift adds nothing,
and the dynamics gives `|ξ|²=½` only as a potential parameter (input) with the phase free.  The
coupling magnitudes are not fixed by the dynamics — consistent with, and now a sharp instance of,
the kinematics/dynamics divide.  *The honest standing of the lepton sector is unchanged from the
roll-up in `15_…md`: amplitude `½` = the L-grade complex norm (structural, value = potential input);
phase `φ` = free residual near the chiral edge.*
