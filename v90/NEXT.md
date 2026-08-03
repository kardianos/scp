# NEXT — standing state, interpretation, and the queue

Updated 2026-08-03. One page: what is HARD, what it MEANS, what to DO.
Detail lives in the campaign docs; this file is the index and the queue.

## Hard results (each gated or logged; doc → detail)

| result | where |
|---|---|
| Battery **56 bars ALL GREEN** through every kernel change (9 apparatus additions this cycle) | `runs/BATTERY.log`, `cmd/battery` |
| **All-modes COE meter** (`p1_meter`): no-channels null exactly zero; blobs merge with v_COE ~1000× below closing — approach carries no net momentum; dense-only COM was 324× inflated | `MOTION.md` |
| Go kernel port faithful at the FP floor (bath/blob byte-identical; physics digits identical everywhere) | `VERIFY.md` |
| **DS tier 0 CONFIRMED on the live substrate** — fringes at parameter-free loci, V 0.652 vs frozen 0.667, 3-seed panel; Huygens wavelets visible in the phase view | `DS.md`, `runs/ds/viz/` |
| **DS tier 1 CLAIMED (battery `ds1`, 8 bars)** — 362 whole-grain clicks over 24 substrate seeds: fringe-phase score s̄ 0.683±0.017 (null 0.5) vs which-slit controls 0.380±0.036; minima 7.7× dark; every click k·ε, k∈{1,2}. Controls sit BELOW phase-blind (envelope peaks at the minima loci) | `DS.md` |
| Instrumented on the way: click threshold Ee ≥ ε/f_conv makes amp a *calibrated window* [3,4] — below it single slits go silent in-gate (first-pass control was reverb-lit); above it cap evap fires and minima fill. Plus a kernel overflow fixed (subnormal-dust flux in pass-7 alignment → stable form, both kernels, battery-gated) | `DS.md` findings 5–7 |
| FCS v3: chunked run streams, interleaved ANLZ instrumentation, lossless codec 1 (0.46–0.55 @174 MB/s), live `-follow` viewing | `FCS.md`, `cmd/volview` |
| Motion: bound objects do not travel (bounded wander; exact vacuum null under drive); blobs merge (undriven closes as fast as driven); collisions = contact-range plastic locks in [2.88,3.13], exactly inert beyond contact | `MOTION.md` |
| **Opacity is unfilled capacity** — field couples to matter only via conversion; saturated matter transparent; headroom object absorbs Δcond +4.9 | `MOTION.md` |
| Composite: **macro-parity selection rule** (interior rung hops contribute π each ⇒ nv ≡ 0 mod 12 for 60° triangles; nv=6 π-frustrated, nv=12 closes ×100) | `COMPOSITE.md` |
| Composite: **T6 U/D asymmetry** (D-quarks are droplets — rung 2.17 > contact 1.53; U-quarks molecules), confirmed in-run | `COMPOSITE.md` |
| Composite: bond = standing flight (lem 0.09–0.34, E at FP floor); compensation cuts fifth drift 8×; **H-stiffness refuted as posed** (drift grows with N — transient flight chirps pitch) | `COMPOSITE.md` |
| **The bath protects composites** (boundary gg 0.9902 compressed onto the rung, interior 2–5× better than vacuum) | `COMPOSITE.md` |
| Nucleon–nucleon residual force = 1–2 contact channels below sep≈9; exactly inert at 11 | `COMPOSITE.md` |

## Current interpretation

**Scale-free structural laws** (won't change under any completion): the
channel-level ladder transfers to composites (boundary arithmetic = the
minimal triangle's numbers); π-parity constrains geometry at BOTH scales
(odd rings dead; composite interiors quantized to even hop counts);
the U/D molecule/droplet asymmetry; opacity = conversion capacity;
bonds carry their energy as in-flight inventory.

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

## The queue (priority order; effort; blocking)

1. ~~**Tier-1 grain harvest**~~ **DONE 2026-08-03** — claimed, 8 bars
   gated (`ds1`), battery 48 GREEN. Tiers 2–3 analogs (eraser, delayed
   choice) are now unblocked (see Also-can-do).
2. ~~**All-modes center-of-energy meter**~~ **DONE 2026-08-03** —
   `p1_meter` flow bookkeeping (5 sites, 4 channels), battery `p1`
   8 bars (56 total GREEN): no-channels null EXACTLY zero; blob
   approach carries no net momentum (v_COE ~1000× below closing);
   dense-only COM was 324× inflated (MOTION.md).
3. **Composite follow-ups:** CQ8 long-run (does protection persist?);
   bath + compensation combined; the mod-12 parity sweep (nv 12/24/36)
   as a battery-bar candidate. *Cheap, unblocked.*
4. **Pauli rate-level control:** distinguishability-blind saturation
   refusal (the negative control that later separates exclusion from
   fullness). *Cheap, unblocked; full exclusion blocked on S2.*
5. **Angular cross-section apparatus** (narrow beam or sector meters)
   to turn ledger absorption into a shadow/differential measurement.
   *Moderate.*
6. **Go kernel ports** of `exp=slit`/`rings`/`blob2` → restore the full
   A/B surface over the new apparatus. *Mechanical, moderate.*
7. **P2 local-clock scheduler** (Go prototype first, per FREECELL §2's
   four conditions) → big boxes for big composites (nv=36+, 3D baths).
   *Large; the composite programme will start to need it at CQ7-scale.*
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
