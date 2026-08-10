# v94 M4 results (for claude face-3 second resume)

**Session ID to resume:** `aebbc7eb-4716-4e30-8c22-484e115428f2`
**Generated:** 2026-08-09.

## M4 — door phase retention on a lit-bath m=+2 ring

**Verdict: NOT MET (and a surprising negative result).** Target was R2 ≥ 0.5
at t=80 with door traffic firing. The door pair does NOT improve phase
retention on a winding ring — the stimulated outtake (face 2) actively
ERODES the winding.

### Setup

`exp=ring ring_n=6 bath=1 seed_mw=2 amp_nat=2 amp_logate=1 hop_order=1
noise_amp=1.0 k_rad=0.05 p_rad=4 wf_on=1 q_detune=1.2 L=16 T=80
seed=20260802`. Door_g swept over {0, 0.3, 0.5, 1.0}.

### Results (the key comparison: door_g=0 baseline vs door_g=0.5)

| t | door_g=0 (R2, W, Em_mean) | door_g=0.5 (R2, W, Em_mean) |
|---:|---|---|
| 0 | 1.000, +2, 0.81 | 1.000, +2, 0.81 |
| 2 | 0.851, +2, 0.81 | 0.851, +2, 0.81 |
| 4 | 0.314, +1, 0.81 | 0.314, +1, 0.81 |
| 8 | 0.142, −2, 0.80 | 0.151, −2, 0.80 |
| 60 | ~0.13, ±1, ~0.9 | ~0.17, ±1, ~0.1 |
| 78 | 0.226, +2, 1.08 | 0.055, +1, 0.07 |
| 80 | **0.620**, +1, 1.16 | **0.237**, −1, 0.07 |

### The mechanism (why the door pair HURTS here)

The m=+2 seed_mw winding creates a coherent phase pattern around the ring.
This phase pattern is precisely the **C > 0 geometry** that triggers the
stimulated outtake (face 2): matter is "ahead" of the register in the
complex plane. The outtake fires (`door_out = 8.56` over the run vs
`door_in = 1.83`), draining the ring's matter into the register and then
to the field. The ring's `Em_mean` drops from 0.81 to 0.07; phase-slip
sites open (`Em_min → 0.0001`); winding coherence collapses.

**This is a destructive feedback loop:** winding ⇒ outtake geometry ⇒
matter drain ⇒ support loss ⇒ winding death. Face 2 is a winding-erosion
machine on this configuration.

### The baseline surprise

The door_g=0 baseline on a LIT bath (`noise_amp=1.0`) actually achieves
R2 ≈ 0.62 at t=80 — far better than v93's cold-bath R2 ≈ 0.02. The hot
bath continuously feeds the ring via condensation; the ring's topology +
the unitary channel re-establish winding coherence after the initial
dephasing. The ring SELF-HEALS on a hot bath (without the door pair).

### Implications for face 3

1. **M4 does not grade face 3's (L1) "one algebra on every edge" claim
   favorably on winding configurations.** The door pair's outtake is a
   winding-erosion mechanism. The J-dependent link metric (RECKONING's
   next face) remains necessary for phase retention.

2. **The stimulated outtake needs a geometry gate** — it should NOT fire
   on configurations where it destroys the phase structure it's supposed
   to carry. Possible: gate on `Em_min > 0` (don't erode below support);
   or: gate on the local winding gradient (don't fire against the
   winding); or: the `chi_load` knob (face 3) self-detunes the register
   away from outtake resonance on bound structures.

3. **The face-2 outtake is geometry-asymmetric**: it helps on random-
   phase configurations (where it returns hot matter to the field — the
   M3b defense) but hurts on coherent-winding configurations (where it
   destroys the winding). The `chi_load` self-detuning (face 3 §III) may
   resolve this: a bound ring at high `|χ|²` would detune its own door,
   reducing the outtake rate exactly where it's destructive.

### What's still in-flight

- **M1** (metabolic propagation, v_M curve) — not yet run
- **M5** (antagonism release) — not yet run
- **M6** (metabolic fixed point) — not yet run

These will follow in a third resume if they land.
