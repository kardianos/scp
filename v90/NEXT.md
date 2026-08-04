# NEXT — standing state, interpretation, and the queue

Updated 2026-08-04. One page: what is HARD, what it MEANS, what to DO.
Detail lives in the campaign docs; this file is the index and the queue.

## Hard results (each gated or logged; doc → detail)

| result | where |
|---|---|
| Battery **85 bars ALL GREEN** through every kernel change (14 apparatus additions) | `runs/BATTERY.log`, `cmd/battery` |
| **Full A/B surface restored (queue #6)**: slit/rings/blob2 + p1/sect/convtag/clicks apparatus in BOTH kernels; battery `abx` 6 bars — blob2+p1 byte-identical INCLUDING drift; slit/ds1/xsec/rings identical in every physics digit (clicks, screen tables, convtag net 7.274023, sect Etot 961.008277), drift column at the FP floor is the only divergence | `VERIFY.md`, battery `abx` |
| **XSEC gated (15 bars)**: headroom object = clean absorber (net_tag +7.27/+2.97/+5.64 across seeds, rough=evap exactly 0; run-differencing understated 4×); true-saturated flat-top = net EMITTER (−7.20/−7.93/−8.76, evap-borne, side-glow rE 1.54) — **opacity is unfilled capacity, with a sign**; cross-section falls monotonically with impact parameter (7.27>4.21>3.42>2.03; second seed too); inert conversions-off object still shades the core (rE 0.79 — an impedance defect that REFLECTS, not delays: 2× gate probe); absorber core shadow 0.50 at the gated seed; per-seed angular ratios are foam-speckle-dominated — ledger claims are the robust ones; absorption is a CLOCK-RATE effect (heavier headroom object absorbs less) | `MOTION.md` XSEC |
| **All-modes COE meter** (`p1_meter`): no-channels null exactly zero; blobs merge with v_COE ~1000× below closing — approach carries no net momentum; dense-only COM was 55–324× inflated (wrong sign on the driven arm) | `MOTION.md` |
| **PAULI-0 gated**: saturation refusal is distinguishability-blind — at cap, admission exactly 0 for identical AND fifth-consonant senders; near cap both admitted, throttled (26× = gate quality, not door policy) | `README.md` stretch, battery `pauli0` |
| **CQ8-long: the bath's protection is transient** — the composite DISSOLVES into its protector (load leaks: Em −36–39% by t=480; bond dies of downstream detune at t≈150–190). Flux-machine interior forced from a fifth direction | `COMPOSITE.md` |
| Go kernel port faithful at the FP floor (bath/blob byte-identical; physics digits identical everywhere) | `VERIFY.md` |
| **DS tier 0 CONFIRMED on the live substrate** — fringes at parameter-free loci, V 0.652 vs frozen 0.667, 3-seed panel; Huygens wavelets visible in the phase view | `DS.md`, `runs/ds/viz/` |
| **DS tier 1 CLAIMED (battery `ds1`, 8 bars)** — 362 whole-grain clicks over 24 substrate seeds: fringe-phase score s̄ 0.683±0.017 (null 0.5) vs which-slit controls 0.380±0.036; minima 7.7× dark; every click k·ε, k∈{1,2}. Controls sit BELOW phase-blind (envelope peaks at the minima loci) | `DS.md` |
| Instrumented on the way: click threshold Ee ≥ ε/f_conv makes amp a *calibrated window* [3,4] — below it single slits go silent in-gate (first-pass control was reverb-lit); above it cap evap fires and minima fill. Plus a kernel overflow fixed (subnormal-dust flux in pass-7 alignment → stable form, both kernels, battery-gated) | `DS.md` findings 5–7 |
| FCS v3: chunked run streams, interleaved ANLZ instrumentation, lossless codec 1 (0.46–0.55 @174 MB/s), live `-follow` viewing | `FCS.md`, `cmd/volview` |
| Motion: bound objects do not travel (bounded wander; exact vacuum null under drive); blobs merge (undriven closes as fast as driven); collisions = contact-range plastic locks in [2.88,3.13], exactly inert beyond contact | `MOTION.md` |
| ~~Opacity is unfilled capacity — saturated matter transparent~~ **superseded by XSEC**: headroom absorbs (+), saturation EMITS (−), inert load reflects; the #29 integral numbers stand but hid all three signatures | `MOTION.md` XSEC |
| Composite: ~~macro-parity selection rule~~ **REFUTED by the 4-point sweep** (2026-08-03): the boundary edge meter weighs interior phases by the ℤ₃ ratio — hop parity is invisible to it; seed residual walks continuously; nv=6's π-pinning is real but unique; nv≥24 loses the boundary at T=200 | `COMPOSITE.md` |
| Composite: **T6 U/D asymmetry** (D-quarks are droplets — rung 2.17 > contact 1.53; U-quarks molecules), confirmed in-run | `COMPOSITE.md` |
| Composite: bond = standing flight (lem 0.09–0.34, E at FP floor); compensation cuts fifth drift 8×; **H-stiffness refuted as posed** (drift grows with N — transient flight chirps pitch) | `COMPOSITE.md` |
| **The bath protects composites** (boundary gg 0.9902 compressed onto the rung, interior 2–5× better than vacuum) | `COMPOSITE.md` |
| Nucleon–nucleon residual force = 1–2 contact channels below sep≈9; exactly inert at 11 | `COMPOSITE.md` |

## Current interpretation

**Scale-free structural laws** (won't change under any completion): the
channel-level ladder transfers to composites (boundary arithmetic = the
minimal triangle's numbers); π-parity constrains geometry at CELL scale
(odd rings dead — gated); the U/D molecule/droplet asymmetry; opacity =
conversion capacity; bonds carry their energy as in-flight inventory.
*(The composite-interior "even hop count" quantization was refuted by
the 2026-08-03 sweep — nv=6's π-pinning is real but does not
extrapolate; see COMPOSITE.md third execution.)*

**One frontier, measured from every side:** flight-load pitch drift vs
tongue width kills the fifth at cell scale, composite scale, and inside
bound composites; blobs move but structures don't; EMF emission is gated
on fifth transport; sharp charge needs the same completion. This is the
**S2 amplitude-completion acceptance surface**, now with quantified bars
waiting at each face: e3b translation coherent, the fifth holds, flavored
composites persist, composite EMF radiates, exclusion becomes testable.

**Second standing requirement (new this cycle):** a stationary composite
pitch requires **steady internal circulation** — the flux-machine
interior. Static seeds cannot hold flavor while their internal flight
decays. Matter that keeps its pitch is matter that keeps its books
flowing — the v89 particle thesis, forced at composite scale by
measurement.

**Where matter lives:** in the bath. Confinement pressure protects both
the boundary bond and the interior. Vacuum composites decohere; embedded
ones hold. Composite work should default to embedded apparatus.

**Matter–field coupling now has three measured faces (XSEC):** an
object below cap is an absorber (clean condensation, cross-section
falling with impact parameter); an object at cap is an emitter
(evaporation, side-glow); an object that cannot convert at all still
reflects (impedance defect — load shrinks lenses). And absorption runs
on the CLOCK, not the mass: load slows the beat, so heavier matter
absorbs more slowly. Opacity, emission, and reflection are all
statements about conversion capacity and conversion rate — never about
"solidity". The angular meter also taught a foam lesson: λ≈3 light on
a unit foam speckles, so single-seed angular ratios wobble ±0.4; the
ledger, not the far field, carries the seed-robust claims at this
scale.

## The queue (priority order; effort; blocking)

1. ~~**Tier-1 grain harvest**~~ **DONE 2026-08-03** — claimed, 8 bars
   gated (`ds1`), battery 48 GREEN. Tiers 2–3 analogs (eraser, delayed
   choice) are now unblocked (see Also-can-do).
2. ~~**All-modes center-of-energy meter**~~ **DONE 2026-08-03** —
   `p1_meter` flow bookkeeping (5 sites, 4 channels), battery `p1`
   8 bars (56 total GREEN): no-channels null EXACTLY zero; blob
   approach carries no net momentum (v_COE ~1000× below closing);
   dense-only COM was 324× inflated (MOTION.md).
3. ~~**Composite follow-ups**~~ **DONE 2026-08-03** (COMPOSITE.md third
   execution): the mod-12 sweep (nv 12/18/24/36 + branches) REFUTED
   the macro-parity rule as an extrapolating law — seed residual walks
   −1.11 rad/step, nv=18 not pinned, nv≥24 loses the boundary; no bar.
   CQ8 A/B: in the bath, vacuum-tuned compensation HURTS (bond 0.9980
   vs 0.9902, interiors 0.59/0.70 vs 0.42/0.47) — the bath retunes as
   it protects. CQ8-long T=480: protection is TRANSIENT — the
   composite dissolves into the bath (load leaks −36–39%; bond dies
   of downstream detune t≈150–190).
4. ~~**Pauli rate-level control**~~ **DONE 2026-08-03** — gated
   (battery `pauli0`, 8 bars): refusal at cap is exact and
   pitch-blind; README stretch bullet has the numbers. Full exclusion
   still blocked on S2.
5. ~~**Angular cross-section apparatus**~~ **DONE 2026-08-04** — gated
   (battery `xsec`, 15 bars; suite 79 GREEN): sector meter + obj_y +
   tag-split conversion ledger; headroom=absorber / flat-top-cap=
   emitter / inert=reflector; monotone impact-parameter profile;
   angular speckle lesson. MOTION.md XSEC has the record.
6. ~~**Go kernel ports** of `exp=slit`/`rings`/`blob2`~~ **DONE
   2026-08-04** — full apparatus surface in both kernels (incl.
   p1/sect/convtag/clicks); battery `abx` 6 bars gates the pairs
   (suite 85 GREEN). VERIFY.md queue-#6 section has the record.
7. ~~**P2 local-clock scheduler** (Go prototype)~~ **PROTOTYPE DONE
   2026-08-04** — `fab/localclock.go` on the real substrate; battery
   `p2lc` 8 bars: four conditions verified (R3 exactly 0; batch ==
   serial exactly 0 at 1–8 workers via the PENDING-min rule — a
   sharpening of v89's eligible-min, which was only exact under full
   eligibility); quiet-region economy measured (event ratio 2.37 →
   6.87 approaching the R=8 dilation bound as the box grows). P2.md
   has the record + the production-build checklist (event queues,
   full-law event decomposition, per-cell RNG streams, local topology).
   The C production engine remains open (P4-adjacent).
8. **S2 amplitude completion on the free substrate** — THE frontier.
   Acceptance surface as above; `v89/s2/` holds the derivation
   micro-model; kappa_reac=1 full pass remains the entry criterion.

## Also can-do (opportunistic)

* volview: ANLZ-table overlay in `-i` (plot ds_screen build-up live);
  turntable/GIF export for records; click-position overlay on the DS view.
* Battery: promote composite structural claims (parity rule, T6, bond
  attractor) to gated bars once the mod-12 sweep lands.
* DS tiers 2–3 analogs (eraser, delayed choice) after tier 1 is green.
* Streams: `runs/streams/` is regenerable by `./battery`; large `.fcs`
  archives to `/space/scp/` per the standing data policy.
