# P2 — the local-clock scheduler (queue #7): Go prototype record

Executed 2026-08-04. Spec: FREECELL §2's four measured conditions,
verbatim (README P2). Code: `fab/localclock.go` (`exp=p2lc`, Go-only).
Gate: battery `p2lc`, 8 bars. v89 evidence: `v89/localclock.c` proved
the execution half on a synthetic degree-7 graph; this prototype
re-proves it ON THE REAL SUBSTRATE and adds what v90 needs.

## What was built

* **The substrate is real:** a 2D dart-throw bath at box L, jam-settled,
  channels from the live contact rule (`topoRefresh`), pitches from the
  law (`w2e = w2/(1+q_detune·x)` under a seeded blob load — the active
  region), executed by the v89 reduced dynamics (phase advance +
  paired antisymmetric channel transfers). The reduced dynamics is the
  execution skeleton, NOT the 12-pass law — this is a scheduler
  prototype, exactly as the v89 suite was.
* **Activity-derived dilation:** active (loaded, detuned) cells step at
  dt; quiet bath cells at R·dt (R=8). Lookahead = 1.25·R·dt, in LOCAL
  TIME (condition 2).
* **Goroutine batch execution** with a deterministic local-min rule —
  the batch IS the serial schedule, executed wide.

## Measured (battery `p2lc`, seed 20260802, L=16 criteria / 16-32-64 scaling)

| criterion | result |
|---|---|
| R1 conservation | drift −5.6e-16 sync / −1.4e-16 async — FP floor, both engines |
| R2 convergence | err(dt)/err(dt/2) = **1.715** vs a dt/64 sync reference — first order, same solution (at Kc=0.2; see finding 2) |
| R3 determinism | rotated, reversed, and rotated+reversed scans: **0.000e+00 exactly** under the (t, index) total order; arrival-order control corrupts at 1.454 |
| condition 4 | max neighbour tick skew **219** while local times stayed inside the lookahead — the counter diverges by construction (a mod-M byte would need M > 438); ordering uses local TIME, never tau |
| batch | **bit-identical to the serial schedule at 1, 2, 4, 8 workers** (max\|dθ\| exactly 0), mean width 10.9, 6610 rounds |
| quiet-region economy | fixed active blob (~60 cells), growing bath: event ratio sync/async = **2.37 / 4.92 / 6.87** at L = 16/32/64 (N = 177/700/2739) — approaching the R=8 dilation bound as the quiet fraction grows. THE P2 claim, measured: cost concentrates on the active region, not the box |

## Findings on the way (both sharpen the v89 record)

1. **The batch conflict rule must take the minimum over PENDING events,
   not eligible ones.** v89's min-over-eligible rule measured 6.2 off
   the serial schedule here: an earlier-but-blocked neighbour must
   still hold you back, exactly as the serial (t, index) order would.
   v89 got exact results only because its lookahead made every event
   eligible at once (its R4). The pending-min rule is exact by the
   schedule-equivalence argument AND by measurement (0.000e+00), and
   degrades gracefully when eligibility is partial — the production
   scheduler must use it.
2. **R2 must be measured in a single-basin regime.** At Kc=0.8 the
   Kuramoto sector is multistable and the "error" saturates at the
   attractor scale (ratio 1.02 — integrator order invisible). At
   Kc=0.2 the ratio is 1.715. Convergence order is a property of the
   scheduler only where the dynamics doesn't amplify ulps into basin
   choices — the same lesson the ds A/B taught (linear optics preserves
   digit identity; chaos would not).

## What the production build still needs (the C engine, P4-adjacent)

* **Event queues.** The prototype scans all N+ne events per round —
  O(NT) per round is why its wall-ms lose to sync despite winning the
  event count 6.9×. The economy is real (events are the cost that
  scales); harvesting it needs a priority structure (bucket queue on
  local time) — mechanical, deferred to the C build.
* **The full law as events.** The 12-pass step must decompose into
  channel-owned transfers (pass S/4/5 already are transfers; pass F
  pair rotations pairwise; pass 6 per-cell beats have natural local
  times = the beat clock). The serial pass ORDER within a tick is a
  constraint the event decomposition must respect or re-derive.
* **Per-cell RNG streams.** The global serial xorshift stream (pass 7
  tumble) is incompatible with any work-skipping; the async engine
  needs per-cell counters (e.g. hash(seed, cell, tau)) — a determinism
  contract change that must be its own battery-gated step.
* **Live topology as events.** G0's global refresh must become local
  (channel birth/death on neighbourhood motion) — the contact rule is
  already local; the spatial hash update is the engineering.

Battery: 93 bars ALL GREEN with `p2lc` gated.
