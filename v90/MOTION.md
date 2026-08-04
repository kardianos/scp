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

## The all-modes center-of-energy meter (2026-08-03; battery `p1`, 8 bars)

The queued refinement, built. `p1_meter=1` (kernel apparatus key,
default OFF — every standing log stays byte-identical, verified by the
full battery). **Momentum is the first moment of conversion**, so the
meter is flow bookkeeping, not a global COM: every channel-borne
transfer contributes dE·Δx at its own link (minimal-image — torus-safe,
no unwrapping), split by channel:

| channel | sites |
|---|---|
| `sp` space | pass-S Es flows |
| `fl` flight | dense deposit (cell→link mid), arrival (mid→cell), rule-α flush (mid→source) — the midpoint convention telescopes to dE·(x_dst−x_src) over any completed transfer |
| `fd` field | pass-F pair rotations (energy moved is pairwise-conserved) |
| `gm` geometry | cells carry their energy when pass-D moves them; link flight rides the endpoint mean |

Conversions (cond/evap/rough/backsplash) change mode AT a cell — no
moment flow, correctly uninstrumented. Cumulative `# P1` per diag
cadence; `# RESULT p1` totals + rates. The all-modes center-of-energy
velocity is v_COE = P1_rate / E_total.

**Measured (L=16, T=80, blob2 Round-B geometry; battery-gated):**

| run | v_COE (x,y,z) | closing | \|v_COE_x\|/closing | dense-only vTotalCOM |
|---|---|---|---|---|
| null (2 cells beyond contact) | **(0, 0, 0) EXACTLY** | — | — | — |
| blob2 undriven | (−3.4e-6, +1.8e-5, +4.1e-6) | 3.5e-3 | **1.0e-3** | 1.1e-3 (**324× the honest number**) |
| blob2 driven kx=3.2 | (+1.7e-5, −2.8e-5, −5.6e-7) | 1.6e-2 | **1.1e-3** | −2.0e-4 |

Three claims, now gated (battery 56 bars GREEN):

1. **No channels = exactly zero moment flow.** Two cells beyond contact:
   every accumulator byte-zero over 60 t.u. Interaction is channels;
   nothing else moves energy, period.
2. **Blob approach is merging, not transport, in the honest ledger.**
   The all-modes center of energy moves ~1000× slower than the blobs
   close — the modes exchange moment locally and it cancels. The Δp of
   "two blobs colliding" is, to 1e-3 of the closing rate, zero.
3. **The Round-B dense-only COM number was mode-exchange artifact** —
   overstates the all-modes v_COE by 324× on the undriven control.
   The honest Δp bar is |v_COE| ≤ 5e-5 (measured ceiling 2.8e-5), and
   the flight+geometry channels dominate the (cancelling) internal
   flows; field ≈ 0, space ≈ 0 at blob scale.

Images: `runs/motion/viz/p1_b2_t0_em.png` / `p1_b2_t80_em.png` (the
mound pair, dense channel). Stream: `runs/streams/p1_b2.fcs`.

**Round-B replica at full scale (L=24, T=300, meter on; the dense
kinematics reproduce the committed Round-B numbers to the printed
digit — the meter is pure ledger):**

| arm | v_COE (x) | dense vTotalCOM | verdict on the dense-only number |
|---|---|---|---|
| undriven | −9.6e-6 | −5.27e-4 | **55× overstated** |
| driven kx=3.2 | +9.9e-6 | −7.82e-5 | **8× overstated, wrong sign** |

Both arms: \|v_COE_x\| ≈ 4e-3 of the closing speed. The honest
all-modes momentum of a blob collision is zero at the few-per-mille
level of the approach rate, at both box scales, driven or not; the
dense-only COM velocity is mode-exchange artifact in magnitude AND
sign. Channel split at T=300: flight dominates (\|fl\| ≈ 20), geometry
second (≈ 6), field ≈ 1 (driven) / 0 (undriven), space ≈ 0.02.

## XSEC — the angular cross-section apparatus (queue #5; opened 2026-08-04)

**Goal:** turn the #29 ledger result ("opacity is unfilled capacity";
Δcond = +4.9 was real but the screen-integral shadow drowned in ±3–5%
geometry effects) into a **shadow/differential measurement**: absorption
must have a *direction* (a forward deficit behind the absorber) and a
*profile* (it must fall off with impact parameter).

**Why #29 could not see a shadow (derived before the redesign):** the
beam was 9× wider than the object (σy=14 vs σ=1.6) AND the object was
sub-wavelength (2σ = 3.2 < λ = 6.98 — sub-wavelength obstacles cast no
shadow in any wave theory). The redesign must fix BOTH: a narrow beam
aimed at the object, and an object larger than the wavelength
(kx=2.0 → λ=3.14; obj_sigma=2.5 → 2σ=5 > λ).

**The threshold mechanism (derived from the kernel before first run):**
pitch detune divides BOTH pitches (freecell.c pass 1: w1e = w1/det,
w2e = w2/det, det = 1+q_detune·x), so the conversion grain
ε = A0eff·w1e/2π **shrinks with load** — a loaded cell fires grains at
a lower field threshold than a vacuum cell (Ee ≳ ε/f_conv scaled down
by det, softened by two-atom credit accumulation across beats). A
headroom object is therefore a *better absorber than the bath it sits
in* under laws_V2g verbatim — no regime override needed for the law
arms. This is the sharpest form of "detection is conversion".

**Apparatus (kernel additions, all default-off/neutral):**

* **Sector meter** (`sect_meter=1`, keys `sect_n/r0/r1/x/y/t0/t1`): an
  annular exposure meter centred on the object — per-sector
  time-integrated Ee·dt over cells in r ∈ [r0,r1), gated [t0,t1),
  sector 0 centred on +x (downbeam). The near-field angular
  distribution of the beam, before diffraction refills the deficit.
  `# RESULT sect` + per-bin `# SECT` rows; ANLZ table when streaming.
* **`obj_y`** — the occulter's y-centre (default box centre): the
  impact parameter b = obj_y − cy0 becomes a config scan.
* **Tag-split conversion ledger** (printed when `slit_obj=1`):
  cumulative rough/cond/evap/backs **at tagged (object) cells** —
  `# CONVTAG` timeline + `# RESULT convtag`. Net field capture at the
  object = cond − evap − rough + backs, per-run, no run-differencing.
* Narrow beam = existing `slit_sy` (no new code).

**Geometry (all arms):** 2D bath L=64, no wall (`slit_mask=3`), beam
src x=20 σx=4 **σy=3** kx=2.0, object at (32, obj_y) σ=2.5, annulus
r ∈ [7,11] on the object, screen x=40 (secondary), law table VERBATIM
(no optics overrides in the law arms), seed 20260802.

**Pre-registered predictions (before first run):**

* **P-XS1 (angular transparency floor, optics arm):** with conversions
  off (q_detune=0, e_cond=99), object-vs-control sector ratios are
  unity to within a small geometry floor: |r(θ)−1| ≤ 0.10 in both the
  forward and backward sector groups — #29's stage-1 transparency,
  made angular.
* **P-XS2 (the shadow exists and points forward):** law arms, headroom
  object (obj_amp=0.5): the downbeam group (|θ| ≤ 45°) shows a deficit
  r_fwd ≤ 0.90 and below the optics floor, while the upbeam group
  stays at unity within the floor — absorption has a direction.
* **P-XS3 (saturation is transparent, angularly):** object seeded at
  the cap (obj_amp=2.06 → peak Em = 0.95·cap·(1+s_pull)/… exactly at
  seed recipe): r_fwd within the optics floor of 1 and |net_tag| ≪
  the headroom arm's — opacity is unfilled capacity, now with a
  direction attached.
* **P-XS4 (differential profile):** net_tag(b) falls monotonically
  with impact parameter b ∈ {0, 2, 4, 8} (beam σy=3: b=8 ≈ 2.7σ is a
  miss): net(8) ≤ 0.1·net(0).
* **P-XS5 (ledger↔shadow consistency):** the forward missing exposure
  and net_tag agree in order of magnitude (ratio recorded; gated only
  if stable across seeds).
* **P-XS6:** conservation ≤ 1e-13 every arm; every standing battery
  log byte-identical with the new apparatus keys at defaults (full
  suite green before any XSEC run is interpreted).

Calibration (beam amp; gate timing at kx=2 group speed; annulus radii)
is apparatus tuning and will be recorded as instrument findings; the
bars gate the frozen design.

### XSEC results (2026-08-04; battery `xsec` 15 bars GREEN; logs `runs/motion/xsec/`)

Beam amp=4 (the ds1-calibrated linear ceiling); 32 runs total: 3-seed
law pairs, 3-seed optics pairs, 3-seed sat20, b-scan at two seeds with
matched controls, flat-top saturation, L=96 long-gate mechanism probe.

**The ledger claims (seed-robust; gated):**

* **A headroom object is a clean absorber.** net_tag = +7.274 / +2.968
  / +5.637 (seeds 20260802/111/314159) — always positive, and always
  PURE condensation: rough = evap = 0 exactly, every seed, every b.
* **Run-differencing understates absorption 4×** (the #29 method):
  global Δcond = +1.81 while the object took 7.27 — the object also
  SHADES the fog downstream, cancelling most of the global difference.
  The tag-split ledger was the missing instrument.
* **A truly saturated object is a net EMITTER.** The Gaussian near-cap
  seed (obj_amp=2.06) is *not* saturated — its shoulders have headroom
  and absorb (+6.95). The flat-top disc (obj_amp=20, clamped at
  0.95·cap, 1.09·cap after pull) inverts the sign: net = **−7.20 /
  −7.93 / −8.76** across seeds, evaporation-dominated (evap 11.8), and
  the shedding is beam-era (evap grows 2.1 → 11.8 over t 16 → 32), not
  seed-settle. **Opacity is unfilled capacity — with a sign: headroom
  absorbs, saturation emits.**
* **The differential cross-section profile.** net_tag falls
  monotonically with impact parameter b ∈ {0,2,4,8}: 7.27 > 4.21 >
  3.42 > 2.03 (seed 20260802) and 2.97 > 1.89 > 1.19 > 0.30 (seed
  111). The b=8 (2.7σ_beam) residual is real light: beam divergence at
  λ=3.14 (σ_eff ≈ 4.3 at the object) plus fog glow.
* **Absorption is a clock-rate effect, not a mass effect:** obj_amp
  0.8 absorbs *less* than 0.5 (6.58 vs 7.27) — load divides both
  pitches, slowing the beat clock: fewer conversion windows per unit
  time. More matter ≠ more opacity even below cap.

**The angular claims (same-seed differentials; s802 tripwires gated,
3-seed spread recorded):**

* **The absorber darkens its forward core**: rE(|θ|≤15°) = 0.500 vs
  control, with back sectors at 0.970 and sector OCCUPANCY at 0.986 —
  the shadow is light, not foam population (the n column separates
  them).
* **An inert object shades too — P-XS1 REFUTED as posed.** With
  conversions off (net = 0 exactly), the object still takes the core
  to rE = 0.786: a dense object is an **impedance defect** — load
  shrinks radii (s_pull), lens areas drop, and the beam partially
  REFLECTS/deflects. The L=96 double-gate probe rules out delay: the
  deficit does not refill (core 0.68 at 2× gate). #29's "transparent"
  was the integral number; angularly the inert object is a partial
  mirror.
* **At the gated seed, absorption deepens the core shadow 0.29 below
  the inert floor** (0.500 vs 0.786). BUT the 3-seed panels show
  single-seed angular ratios are **foam-speckle-dominated** (hdr core
  0.50 / 1.23 / 0.75; optics floor 0.79 / 1.14 / 0.75 — λ ≈ 3 light
  on a dmin=1 foam speckles, and the object's angular effect rides the
  speckle). The angular bars are therefore same-seed regression
  tripwires; the seed-robust statements are the ledger's. A many-seed
  angular harvest (ds1-style) is the upgrade path if an angular CLAIM
  is ever needed beyond the gated apparatus.
* **The emitter pours light sideways**: sat20 rE_side = 1.539 (its
  evaporated take re-radiates isotropically) while its core still
  shades (0.60, strongest lens). Absorber, inert lens, and emitter
  have distinct angular + ledger signatures — detection-is-conversion
  now has an angular face.

P-XS5 (ledger↔shadow consistency) recorded, not gated: screen-integral
deficit 2.1% vs absorbed 4.4% of beam; forward-annulus missing exposure
14.4 units vs net_tag 7.27 — same order, no clean transport factor on a
foam. Conservation ≤ 2.2e-15 across all 32 runs (P-XS6 ✓).

Images: `runs/motion/viz/xsec_hdr_avg_ee.png` (the beam with the
occulter — teal tags — and its darkened wake), `xsec_ctl_avg_ee.png`
(the unobstructed beam), `xsec_hdr_avg_em.png` (the condensation trail
— the fog is a cloud chamber: the beam leaves a matter track).
Stream: `runs/streams/xs_hdr.fcs` (battery-regenerated).

**Battery:** experiment `xsec`, 8 arms, **15 bars GREEN** (suite 79).
Kernel apparatus added this campaign (battery-verified byte-neutral at
defaults, twice): the annular sector meter (`sect_*`, exposure + occupancy
per sector), `obj_y` (impact parameter), and the tag-split conversion
ledger (`# CONVTAG` timeline + `# RESULT convtag`).
