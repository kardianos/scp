# L-1 FINDINGS — first implementation + the design questions it surfaced (v92)

L-1 IMPLEMENTED in both kernels (C `freecell.c` + Go `fab/step.go`), inert at
defaults, per L0_DESIGN.md §3. Knobs `amp_drv`(0) / `amp_mmin`(2) /
`amp_tau`(0, existing). The Phase-M shadow is promoted from meter to driver:
where it composes coherently on a chord link, it biases the dense want toward
the coherent direction (signed, saturating, slot-borne).

## Verified (the engineering is sound)

- **Byte-inertness:** full battery at `amp_drv=0` is ALL GREEN 87 (the bias
  branch is skipped; only the purity header line added, now pinned).
- **C≡Go faithful with L-1 firing** (`amp_drv=0.5 amp_tau=30`): all 5 abx
  bars PASS (blob2/sl​it/ds1/rings/xsec identical up to drift col); Go drift
  floor 1.17e-13.
- **Energy-conserving:** every drift bar at the floor (conserve/blob/ds/ds1/
  p1). The deposit ledger balances; the outflow limiter caps. The zero-sum
  per-cell renorm (registered escalation) was NOT needed — conservation holds
  without it.
- **L1-C anti-ignition PASSES — and stronger than designed:** the warm bath
  is **byte-identical** at `amp_drv=0` vs `0.5` (T=40: births 4611; T=200:
  births 7780, deaths 7106, z_live 17.19 — every digit matches).

## Design finding 1 — the chart-order gate (D-2a) is redundant for anti-ignition

The bath is byte-identical at **both** `amp_mmin=2` AND `amp_mmin=1`. The
bath's transport is balanced (c0≈c1 per link ⇒ bias≈0), so the feedback is
structurally starved **by the balance, not by the chart gate**. D-2a (gate on
m≥2) is therefore redundant given D-2b (bounded bias on a balanced medium).
This is good news (anti-ignition is robust) but it changes the gate's role:
it is not the anti-ignition mechanism (the balance is).

## Design finding 2 — the gate makes unison transport unreachable

The e3b tilted blob and the p1_b2drv driven two-blob are **unison (m=1)**
matter (same-pitch cells, phase-tilt drive). At `amp_mmin=2` they are
byte-identical to `amp_drv=0` — the gate excludes them, so L1-A (e3b
translation) **cannot be reached** with the gate as specified. At `amp_mmin=1`
unison transport engages (p1_b2drv v_COE 5.16e-5→4.32e-5; the control
p1_b2ctl 3.71e-5→5.77e-5) — but the bath stays byte-identical (finding 1).

## Design finding 3 — the magnitude-bias is perturbative, not the coherent current

At `amp_mmin=1 amp_drv=0.5` the transport effect is small and ambiguous (the
driven signal decreased, the control noise rose) — not the robust coherent
translation L1-A wants. Mechanism read: the bias `(c0−c1)/(c0+c1)` uses the
shadow **magnitude** imbalance per direction. But the e3b tilt creates a
phase **gradient**, not a magnitude-coherence imbalance — forward deposits
from cells at different x carry different phases (kx·x), so they do NOT
sum to a large coherent magnitude. The bias is blind to a phase-gradient
current. D-1's "translation IS the current" likely needs the shadow
**phase** (the in-phase current Re(A) in the transport frame), not the
magnitude imbalance. This is the L0_DESIGN D-1 refinement the first run
surfaced.

## Status of L1-A (the transport acceptance)

Not yet measured cleanly. The e3b direct run is slow under radiance (the
glow era bloats the blob) and the short-T form gives `nfsamp<6` (tagged
centroid not captured before the blob churns). The p1 proxy (above) shows
the magnitude-bias is perturbative. **L1-A is not met** — the bias form
needs the phase-current refinement before the acceptance is testable.

## Proposed next step (the design refinement, for the user)

Revise the L-1 driver from **magnitude-imbalance** to **phase-current**:
bias the want by the shadow's in-phase projection along the link
(`bias ∝ Re(A_dir)` in the transport phase frame, or equivalently the
directional Poynting-like current), still saturating + slot-borne. This
responds to a phase gradient (the tilt drive) the magnitude-bias ignores.
Re-evaluate the gate: since anti-ignition holds at `amp_mmin=1` (finding 1),
the gate can be relaxed to include unison (where e3b/p1 live) — OR kept at
m≥2 and L1-A reframed as chord transport (merging L1 into L2). The phase-
current form + the gate decision are the two open design questions.

## What is committed

The L-1 apparatus (knobs + bias + C/Go mirror + purity pin), inert at
defaults — sound engineering, battery-green, C≡Go-faithful, anti-ignition-
safe. NOT adopted (amp_drv=0 default). The findings above are the design
feedback the first implementation owed the L-0 round.

---

## Addendum — the phase-current revision (user-directed: finding-3 fix + both gate tracks)

Driver revised magnitude-bias → **phase-current**: bias by the shadow's
in-phase projection onto the receiver transport frame
`J = Re(A·e^{-i·m_rcv·θ_rcv}) = |A|·cos(closure_phase)` (responds to the
transit-asymmetry term that a uniform tilt breaks). Both kernels mirrored.

**Verified:**
- **Byte-inert** at amp_drv=0 (battery ALL GREEN 87).
- **Anti-ignition ROBUST** across BOTH driver forms and BOTH gate settings:
  bath byte-identical (births 7780) at amp_drv=0 vs 0.1, mmin=1, phase-
  current. The bath's random closure phases average J→0; the feedback is
  starved by balance, independent of the gate or the driver form. **The
  Phase-L substrate is safe.**
- **Conservation** holds (drift floor).
- **The driver ENGAGES** (unlike the magnitude form): xsec headroom
  absorption rises monotonically with amp_drv (6.98 → 7.76 @0.1 → 8.12 @0.2
  → 10.45 @0.5). The phase-current moves the dynamics.

**Still open:**
- **C≡Go at strength:** physics stays byte-identical, but the feedback
  AMPLIFIES C/Go FP differences in meter/diagnostic lines (QATOM `e`, P1
  `tot`, credit) — 1-ulp, same envelope class as drift/p1. Masked at
  amp_drv≤0.1; degrades (more lines) at higher strength. A tolerance-based
  abx comparison may be needed for production-strength sweeps.
- **L1-A NOT yet met:** the phase-current engages but p1_b2drv v_COE
  *decreased* (5.16e-5→1.36e-5 @0.1) and xsec over-absorbs — coherent
  translation is not yet cleanly produced. The e3b direct measurement
  remains practically blocked (slow under radiance; nfsamp<6 at short T),
  so L1-A is not yet cleanly measured against either form.

**Reading.** The amplitude-completion transport (L1-A) is the programme's
central difficulty, now confirmed empirically: two driver forms (magnitude,
phase-current) are safe and engage the dynamics, but neither cleanly
delivers "translation IS the current." The safe substrate is built; the
transport PHYSICS needs deeper work — a clean e3b signal to tune against,
and likely a driver form that reads the transit-asymmetry current directly.
The chord track (mmin=2, L2's domain) is untested and may be where the
shadow's phase composition is actually load-bearing.
