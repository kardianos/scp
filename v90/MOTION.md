# MOTION — particle movement and collisions on the free substrate

Living campaign record (opened 2026-08-03, user directive: "encode
particle movement and collisions; quantify the current behavior").
Baseline first: what bound objects do **today**, on existing apparatus,
before any new mechanism. Kernel of record; standing law; no law change.

Standing context carried from v89: motion is regular conversion, not
displacement (PRINCIPLE); momentum is the first moment of conversion
(center-of-energy theorem); e3b translation is the one recorded-not-gated
experiment — the forward window exists at retuned tilt but seed variance
persists, and the mechanism is expected to sharpen only when "translation
IS the current" (amplitude completion). Collisions at rate level are
therefore expected to be **contact-mediated relaxation**, not ballistic
scattering.

## The task program

| # | task | state |
|---|---|---|
| 1 | baseline quantification (this file, §below) | **measured 2026-08-03** |
| 2 | two-blob collision apparatus (`exp=blob2`) + per-blob COM/exchange meters + total center-of-energy Δp bookkeeping | **executed (Round B)** |
| 3 | field–matter cross-section: aimed packet through an embedded object (occultation analog) | **executed (Round B): opacity is unfilled capacity** |
| 4 | capture/fusion: two rings swept through the bond capture range | **executed via COMPOSITE CQ1** (boundary lock, no fusion) |
| 5 | driven-object transport: e3b tilt on bound objects | **executed (Round B): null confirmed** |

## Pre-registered predictions (before the new runs)

* **P-MO1 (undriven wander):** bound objects show *bounded* COM wander —
  com_dev(t) saturates at a plateau ≪ object size; no ballistic
  component (the substrate has no velocity state to carry one).
* **P-MO2 (driven translation, reproduction):** the v90 kernel of record
  reproduces the committed e3b number exactly (speed 2.580e-3,
  cos 0.9508 at kx=3.2, L=24, amp=0.5, seed 20260802) — the slit/pin/FCS
  apparatus additions did not perturb any standing stream.
* **P-MO3 (two-object interaction):** the tri2 pair relaxes to the known
  fixed separation ≈ 3.009 from below (seeds 2.2–3.0 pushed OUT,
  overdamped, no oscillation); seeds beyond the bond capture range stay
  put — **no long-range attraction** (the space law pushes; nothing
  reaches out). Fixed point one-sided.

## Results (2026-08-03; logs `runs/motion/`; kernel of record)

**R1 — undriven wander (P-MO1 CONFIRMED).** Bound objects do not travel.
* ring6, vacuum, T=300: COM meanders *bounded* — com_dev 0.005 (t=50),
  0.013 (t=100), 0.024 (t=200), then **returns to 0.0015 (t=300)**;
  amplitude ≤ 0.025 = 1.6% of the object's rms radius. Zero ballistic
  component.
* UUD, vacuum, T=300: one plastic settle displaces the COM by 0.360 (the
  known slide-to-contact reconfiguration), then **exactly frozen**
  (0.3598 to 4 decimals for the remaining 250 t.u.).

**R2 — driven translation (P-MO2 CONFIRMED, with an instrument lesson).**
* At the committed cadence (diag_every=200) the v90 kernel reproduces the
  committed e3b number **to the last digit**: speed 0.002580,
  cos 0.9508, v=(2.45e-3, 3.24e-4, 7.30e-4) — and the v89 binary run
  today agrees byte-for-byte with the v90 kernel. The slit/pin/FCS/zstd
  apparatus additions perturbed nothing.
* Instrument lesson: the drift FIT is cadence-sensitive — cos 0.9508 /
  0.9403 / 0.9336 at diag_every 200/100/50 on the *same physics stream*
  (the x-component is robust at 2.45e-3 throughout; the transverse fit
  noise is sampling-cadence sensitivity, part of §9.7's recorded ±1e-3).
  diag_every is part of the fit apparatus and is now printed in every
  run header (both kernels).

**R3 — collisions (P-MO3 partially confirmed; sharpened).**
tri2 pair, vacuum, T=300, seed separations {2.2, 2.6, 3.0, 3.2, 3.6, 4.2}:

| seed sep | final sep (t=300) | behaviour |
|---|---|---|
| 2.2 | 2.8818 | pushed out, frozen by t≈30 |
| 2.6 | 3.0175 | pushed out, frozen by t≈30 |
| 3.0 | 3.1308 | pushed out, frozen by t≈30 |
| 3.2 | 3.2000 | **inert** (unmoved to 4 decimals, 300 t.u.) |
| 3.6 | 3.6000 | inert |
| 4.2 | 4.2000 | inert |

* **One-sided ✓ and exactly inert beyond contact**: no channels beyond
  the nearest-cell cutoff cfac·(rᵢ+rⱼ) ≈ 1.9 means no interaction of any
  kind — no long-range attraction, confirming the space-law thesis.
* **Overdamped ✓, fast**: relaxation constant τ ≈ 2 t.u. (deviation
  shrinks ~3× per 2 t.u.), frozen within ~30–50 t.u., no oscillation.
* **No universal fixed point ✗ (sharpens the earlier claim)**: finals
  are seed-dependent — the pair freezes into a *history-dependent
  plastic lock* in the range **[2.88, 3.13]**, not at one separation.
  The earlier committed 3.0093 (seed 2.6, Δk probe) is one member of
  this family. Any future gated bar must encode the range and the
  inertness boundary, not a single number.
* Collision taxonomy at rate level: **no bounce, no capture from beyond
  contact, no fusion — only contact-range plastic settle.** Ballistic
  scattering presupposes transport, which is the amplitude-completion
  frontier (task 5).

Battery after the header-print change: ALL GREEN (40 bars).

## Round B (2026-08-03): the queued tasks executed

**#31 driven-object transport — null confirmed, with a plastic nuance.**
Ring6 with the e3b tilt (kx=3.2): in vacuum, an EXACT null (com_dev
0.003 over 200 t.u. — the drive only decoheres bonds). In a bath: a
one-time plastic displacement (~0.7) during the tilt-induced
reconfiguration (the ring stays connected, conn = 1.000, but is crushed
— shape (1.13,1.13) → (0.43,0.66), retention 0.81), then bounded wander
0.6–0.9 with no trend. **Blobs drift continuously under tilt; bound
objects do not** — transport of structure remains an amplitude-completion
question, now measured from both sides.

**#29 field–matter cross-section — "opacity is unfilled capacity."**
Three stages, one law: (1) conversions off (the optics regime,
e_cond=99): a dense occulter is TRANSPARENT to the field (exposure
+5%/+3% vs control — the field couples to matter only through
conversion). (2) Conversions on, object at capacity (amp 1.2 ≈ 0.95cap):
still transparent — everything it converts it must re-shed (Δcond +0.5).
(3) Conversions on, object with headroom (amp 0.5): absorption is real
in the ledger — Δcond = **+4.9** (2.7% of the beam) into the object.
The screen-integral shadow is swamped by ±3–5% geometry effects at this
beam/object ratio (σy=14 beam vs σ=1.6 object) — angular-resolved
apparatus is the follow-up. Deep link: detection is conversion
(DOUBLESLIT §3), so **opacity is conversion capacity** — saturated
matter cannot absorb.

**#28 two-blob collision — approach is MERGING, not transport.**
`exp=blob2` (two e3b-class blobs, sep=7, per-blob tags, per-blob and
total dense-COM meters), driven (opposing tilts kx=3.2) vs undriven:

* Closing speed ≈ **0.0021–0.0024 in BOTH** — the undriven control
  closes as fast as the driven pair. Two adjacent blobs approach by
  *profile merging* (each spreads into the other's territory; the
  energy-weighted COMs drift together), not by ballistic motion. The
  tilt adds only ~1e-3 net closing (sep at t=300: 6.235 driven vs
  6.554 control, from 7.0−2σ overlap start ≈ 6.2 effective).
* No bounce, no repulsion — dense mounds interpenetrate freely (they
  are occupancy patterns, not bodies). Full merger extrapolates to
  ~2500 t.u. at this rate.
* Δp bookkeeping: total dense COM moves ≤ 7.8e-5/t.u. (driven) —
  ~30× below the per-blob speeds — but the control shows 5.3e-4:
  the dense-only COM is an INCOMPLETE center-of-energy (mode exchange
  with field/space shifts it without momentum meaning). The full
  all-modes center-of-energy meter is the queued refinement before any
  Δp bar is set. Conservation at the FP floor in both runs.

**Round B synthesis:** every mover at rate level is a *pattern* process
— blobs merge and drift as occupancy waves; bound structures do not
translate at all (one plastic kick, then bounded wander); the field
couples to matter only through conversion, and absorption requires
capacity headroom. Transport of STRUCTURE — the thing a collision
experiment ultimately wants — is on the amplitude completion, from
every direction we have measured.
