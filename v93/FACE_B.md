# v93 — FACE B: lit-bath refill + the arg(psi) door (amp_door)

2026-08-08. RESUME §7-B / RING_DNLS §A.5. A.5 found the V3a dark bath gives
cond=0 (the workfn gate never opens), so amp_door was UNTESTED where it
matters. This face lights the medium (noise_amp seeds field -> door fires)
and asks: does amp_door=1 (coherent arg(psi) merge) retain the ring winding
under door traffic, where amp_door=0 (magnitude-only refill) scrambles it?
Prediction (RESUME §7-B): amp_door=0 scrambles; amp_door=1 the only hold.

## B.1 Lit bath (noise_amp=2, bath=1, V3a) — contact dephasing wins first

Door fires (cond=7-52 vs 0 dark). But the ring winding dies by t~8 in ALL
arms — the A.2 contact-dephasing timescale (~8 t.u.) beats the door's beat
cycle (~7 t.u.), so the door never fires enough ON THE RING before the
winding is gone. amp_door=0 vs 1 are indistinguishable (W=+1 R2=0.318 @t=4
both). The cond=41-52 is overwhelmingly BATH condensation (440:1 cells).

| arm | cond@80 | W@4 | R2@4 | W@8 | keep@80 |
|---|---|---|---|---|---|
| additive | 7.5 | +2.0 | 0.975 | -2.0 | 0.99 |
| uni d0 (amp_door=0) | 52.2 | +1.0 | 0.318 | 0.0 | 0.20 |
| uni d1 (amp_door=1) | 41.1 | +1.0 | 0.318 | 0.0 | 0.17 |
| dark control | **0.000** | — | — | — | — |

Dark control confirms A.5 (cond=0, door inert). The lit-bath test is
confounded: contact dephasing, not the door, kills the winding.

## B.2 Fair test: vacuum ring + Strang + field seed (door fires ON the ring)

Vacuum C6 (bath=0), Strang (hop_order=1), T=200. Face A showed Strang holds
the winding to t>=1000 here — so the door finally gets a fair test (the
beat cycle ~7 t.u. fires ~28 times over 200 t.u. on the ring itself).

| arm | cond@200 | W holds until | R2 trajectory |
|---|---|---|---|
| d0 na0 (no field, no door) | 0 | **t>=200** (R2 1.000->0.999) | baseline HELD |
| d1 na0 (amp_door=1, no field) | 0 | **t>=200** (R2 0.999) | **byte-inert** (== d0 na0) |
| d0 na2 (door fires, amp_door=0) | 5.44 | **~t=16** | -> wanders 0/pm1, R2 0.1-0.8 |
| d1 na2 (door fires, amp_door=1) | 5.44 | **~t=14** | -> wanders, same as d0 |

1. **Strang holds the vacuum winding to t=200** with no field (R2 0.999) —
   independent confirmation of Face A. amp_door=1 is byte-inert when the door
   doesn't fire (d1_na0 == d0_na0 exactly).
2. **Door traffic destroys the winding by t~14-16 regardless of amp_door.**
   The ~28 condensation events (cond=5.44) scramble W; amp_door=1 (coherent
   arg(psi) merge) is indistinguishable from amp_door=0 (magnitude-only).

## Why amp_door=1 does not help — the phase-random field

noise_amp seeds field at RANDOM phase. Each condensation injects matter at a
random field phase. amp_door=1 composes that phase coherently with psi_m --
but composing a RANDOM phase "coherently" is still a random phase injection,
identical in effect to amp_door=0's magnitude dilution. **The coherent door
only helps if the condensing field is phase-CORRELATED with the ring winding**
(a coherent driver), not for the thermal/random-phase field of a warm bath.

## Verdict

**The prediction is REFUTED.** amp_door=1 is NOT a retention rescue: under a
warm (phase-random) medium the door traffic destroys the winding by t~15-20
whether the door writes phase coherently or not. The arg(psi) door's phase
handling is irrelevant for a thermal refill. (It remains the right door for a
COHERENT driver — an untested, separately-seeded phase-matched field.)

The retention ceiling is NOT the door. It is (a) contact dephasing in the
bath (A.2, ~8 t.u.) and (b) the door's own random-phase injection when the
medium is thermal. Both are independent of amp_door.

## Files
- Script: `runs/faceB.sh` (lit-bath arms) + inline vacuum arms ->
  `runs/faceB/{lit_*,dark_*,vac_strang_*}`.
