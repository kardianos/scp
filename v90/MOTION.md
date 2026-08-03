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
| 2 | two-blob collision apparatus (`exp=blob2`) + per-blob COM/exchange meters + total center-of-energy Δp bookkeeping | queued |
| 3 | field–matter cross-section: aimed packet through an embedded object (occultation analog) | queued |
| 4 | capture/fusion: two rings swept through the bond capture range (even+even parity predicted lockable) | queued |
| 5 | driven-object transport: e3b tilt on bound objects (pre-registered expectation: null at rate level; blobs translate, objects hold) | queued |

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
