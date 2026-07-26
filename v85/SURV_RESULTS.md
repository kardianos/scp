# SURV-1 / SURV-2 results — erosion-floor survival assay (Q-PIN test)
**Date:** 2026-07-25 · **Runs:** local CPU (topo1_local), N=64, L=15, g=0.05,
η=0, absorbing BC. Raw data: `/space/scp/v85/topo1/out/surv_Q*_diag.tsv`
(SURV-1, T=400, 8 balls) and `surv2_Q*_diag.tsv` (SURV-2, T=2000, 4 balls).
Seeds: `/space/scp/v85/topo1/seeds/surv_*.sfa`.

## Design (per TOE_v85.3 ladder rung 2, TOE_FEEDBACK R6)
Test the erosion-floor Q-PIN candidate: q₀ = Q_thr(g), the charge where
E/Q crosses m (157.4 at g=0.05) — below it, complete evaporation is
energetically open; the prediction was a **survivor floor** at Q_thr under
agitation. Agitation implementation: ×1.15 amplitude-scaled profiles
("self-agitation": ~30% excess energy+charge the ball must shed). Note the
scaling multiplies seeded Q by 1.3225: census points (actual seeded Q) =
{117, 136, 156 | 184, 212, 223, 265, 353}, bracketing Q_thr = 157.4.
(The θ-bath tool is inert at η=0 — matter-sector bath tool does not exist;
self-agitation was the config-pure variant.)

## SURV-1 (T=400, 8 balls)
| seed Q | final Q | Δ% |
|---|---|---|
| 117.1 | 114.5 | −2.2 |
| 135.6 | 134.2 | −1.0 |
| 156.1 | 154.9 | −0.8 |
| 184.5 | 183.0 | −0.8 |
| 212.4 | 210.5 | −0.9 |
| 223.3 | 221.2 | −0.9 |
| 264.9 | 262.0 | −1.1 |
| 353.6 | 348.7 | −1.4 |

No floor; mild rate gradient with the lightest eroding fastest (predicted
direction) — underpowered at this depth.

## SURV-2 (T=2000, sub-threshold trio + control)
| seed Q | final Q | Δ% | late-run rate |
|---|---|---|---|
| 117.1 | 97.3 | −16.9 | accelerating |
| 135.6 | 107.3 | −20.8 | accelerating |
| 156.1 | 109.4 | −29.9 | accelerating |
| 212.4 | 131.7 | −38.0 | decelerating late |

## Verdicts
1. **No erosion floor observed.** All agitated balls erode continuously
   through T=2000 with no parking at any special charge; the 212 control
   crossed Q_thr = 157.4 **smoothly, with no dynamical feature**. As a
   *dynamical selector under agitation*, the erosion-floor Q-PIN candidate
   is unsupported: NEGATIVE-leaning-inconclusive.
2. **Ordering inverted vs prediction at depth:** heavier agitated balls
   (carrying more absolute excess/ringing energy) shed hardest — erosion
   rate tracks stored agitation energy, not proximity to Q_thr.
3. **Caveat that keeps the floor alive as an end-state claim:** Q_thr
   governs *quiescent* energetics; these balls were still hot at T=2000
   (ringing energy funds emission regardless of Q). Erosion-to-exhaustion
   (T ≫ 2000, rates → 0) would be needed to test whether cooled survivors
   park preferentially above Q_thr. Not scheduled — low expected value
   after finding #2.
4. **Cross-experiment synthesis (with X11-PB, running):** the orbiting
   branch-bottom satellite in X11-PB — *continuously* agitated by tides
   and annihilation — showed accelerating runaway erosion at Q ≪ Q_thr,
   while SURV-2's isolated sub-threshold balls decelerated early. Sustained
   external agitation, not charge, is the operative erosion variable.
   Any surviving Q-PIN mechanism must therefore be sought in **agitated
   environments** (condensate fragmentation, collisions — Stage-5 epochs),
   not quiescent spectra ("selection is cosmogonic, not spectroscopic").

## Status of Q-PIN after this assay
The TOE_v85.3 Q-PIN lock stays open with its only coherent candidate
weakened: erosion-floor selection does not operate dynamically at the
tested agitation. Remaining directions: cosmogonic selection during
fragmentation (Stage-5 simulation), or the fiber/topological routes
(TOPO/G2 threads). The COUNT/HIER locks are unaffected.
