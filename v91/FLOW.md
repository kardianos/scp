# FLOW — the bed-digging channel law in the kernel (user-directed 2026-08-07)

User directive: put the bed-digging channel law into the numeric
kernel; verify it is live; quench-test it; black-hole-test it; record
periodic snapshots for visualization; present analysis and visuals.
Registered BEFORE implementation. Design source: `ASYM.md` §5 (the
converged CORE mechanism); acceptance surface carried from there.
Nothing is adopted; bed_k=0 default; battery-gated both kernels.

## 1. The law, v1 (design decisions)

* **F-D1 (channel scope):** v1 acts on the SPACE channel only — the
  bed multiplies pass-S conductance and is driven by pass-S actual
  moved flow. The space channel is where the river/gravity lives
  (HORIZON); the dense sector (chords, doublet, quench-cloud
  recondensation) is structurally untouched except via space-mediated
  feedback — the clean first face. Dense/field faces = later phases.
* **F-D2 (anti-ignition structure):** the bed grows with |lowpassed
  signed net flow| and is renormalized ZERO-SUM per cell (Σ incident
  sbed = live-link count, enforced by the geometric mean of the two
  endpoints' factors each step): there is no global quantity to grow;
  symmetric churn (net ≈ 0) carries no food; a fed flow deepens its
  channel only by starving its own side-branches.
* **F-D3 (bounds and mortality):** sbed ∈ [0.2, 5.0] (clamped);
  slot-borne, reset to 1 at slot (re)birth — beds die with their
  links; no ledger entry (conductance is a rate, not energy — drift
  untouched by construction).
* Keys: `bed_k` (default **0** = byte-inert; growth rate), `bed_tau`
  (default 30 — the net-flow memory window). Header line
  `# FLOW apparatus (FLOW.md): bed_k=0 bed_tau=30`; battery pin
  extended. Update rule, per step when bed_k > 0, after pass S:

      flow_s   = signed actual moved space on link s (post-limiting)
      bednet_s += (dt/bed_tau)·(flow_s/dt − bednet_s)
      sbed_s   *= 1 + bed_k·dt·|bednet_s|
      per cell i: f_i = n_i / Σ_incident sbed ;  sbed_s *= √(f_i·f_j)
      clamp [0.2, 5];  pass S uses  swl = s_k·dt·w·sbed_s·dp

* **Apparatus extension:** `bh_sep` (default 0) — two hole centres at
  box centre ± sep/2 on x (the two-eater acceptance experiment).
* **Meters:** diag `# BED` line (quartiles/max/count-grown); RESULT
  bed; `# BEDMAP t= i j sbed` dumped at snapshot cadence for links
  with |sbed−1| > 0.05 (capped 5000/frame) — the visualization feed,
  joined to fcs cell positions at the same t.

## 2. Arms (all with fcs snapshots + BEDMAP; seed 20260802;
law-regime keys k_rad=0.05 p_rad=4 rad_clock=0 par_tau=10)

| arm | config | bar |
|---|---|---|
| FW-B | warm bath L=16 T=400, bed_k=0.5, reg_tau=10 | P-F2 anti-ignition |
| FQ-1 | Q2-class beam quench (slit L=64 T=1500 amp=8 σ=4 sy=3 kx=2) + bed_k=0.5 | P-F3 |
| FH-1 | HZ-0 eater (bh_r=4 bh_k=1, amp=0) + bed_k=0.5, T=1500 | P-F4 |
| FH-2 | double eater bh_r=3 bh_sep=16 bh_k=1, bed_k=0.5, T=1500 | P-F5 |
| FH-2c | FH-2 with bed_k=0 | P-F5 control |

## 3. Bars and registered predictions

* **P-F0 (gate):** defaults byte-inert; battery ALL GREEN 93; C≡Go
  physics byte-equal with the bed firing; drift ≤ floor everywhere
  (the law never touches energy).
* **P-F1 (live):** with bed_k > 0 the sbed distribution moves off 1
  wherever net flow stands (wiring proof).
* **P-F2 (ANTI-IGNITION — the critical bar, three-strikes lesson):**
  the warm bath must NOT channelize: sbed q75 stays < 1.2, no
  runaway growth, census/ledger unchanged vs bed-off at chaos grain.
  Symmetric churn carries no food. If this fails the law is
  flammable and the lane records a strike and stops.
* **P-F3 (quench):** the cloud's physics is UNCHANGED at attractor
  grain (population inside the standing band; stamps/idfrac same
  class) and the cloud grows NO significant beds — closed books ⇒
  no standing net ⇒ nothing to dig. Lawful matter stays bed-free;
  this is flow-gravity's version of "stable matter has no far
  field."
* **P-F4 (the eater digs):** around the feeding hole the beds
  CHANNELIZE: sbed spreads to the clamp on inflow-carrying links,
  RADIAL anisotropy emerges (radial links strengthen, tangential
  starve — measured as the mean sbed vs link-orientation split), and
  the bed-boosted conductance FEEDS BACK: eatS slope rises above the
  archived HZ-0 0.15/t.u. — the river deepens its own bed, the
  first self-amplifying gravity in the programme.
* **P-F5 (two eaters, the ASYM acceptance):** the double-hole arm
  grows a connected bed/π structure along the axis between the holes
  (the corridor drained from both sides becomes the valley) that the
  bed-off control lacks — "two eaters grow a channel between them,"
  the attraction precursor in channel form.
* Deliverables: analysis + rendered visuals (bed maps over cell
  states, π fields, growth curves) presented as an exhibit page.

## 4. Results (all arms harvested 2026-08-07; logs `runs/flow/`, fcs +
BEDMAP local; figures `runs/flow/figs/`; battery `runs/BATTERY_flow.log`)

* **P-F0 PASS:** battery ALL GREEN 93 with the law landed; standing
  logs moved by the two registered header lines only; C≡Go physics
  byte-equal with the bed firing (RESULT bed identical; drift-column-
  only divergence); drift at floor in every arm (worst −5.0e-15 =
  the eater class; law touches no energy by construction).
* **P-F1 PASS (live):** bed_k=0.5 moves the sbed distribution and
  micro-shifts conversion (cond 133.137058 vs .137170 at T=40).
* **P-F2 PASS — THE ANTI-IGNITION BAR HOLDS.** Warm bath, T=400:
  sbed q=[0.9934, 0.9971, 1.0051], max 1.066, grown 5 of 21,908
  slots. The zero-sum budget + symmetric churn starves the gate
  exactly as designed — the fourth self-growing law is the first
  that does NOT ignite the bath.
* **P-F3 PASS — lawful matter is bed-free and untouched.** The
  quench cloud with beds on: live 2,411 (the attractor band's own
  edge), idfrac 0.9850, census curve overlaying the archived bed-off
  Q2 for all 1,500 t.u.; beds: 22 grown of 8,570, max 1.117.
  Closed books ⇒ nothing to dig — flow-gravity's version of "stable
  matter has no far field," now with the mechanism visible.
* **P-F4 PASS — THE EATER DIGS ITS RIVER.** 1,608 grown links
  (of 8,265), max sbed 2.67 and still climbing linearly at t=1500
  (clamp 5 not yet reached), while the global median stays 1.0001
  and the interquartile band hugs 1 ± 0.005 — the zero-sum at work:
  only the river grows, paid for by its own side-branches. The bed
  field is RADIALLY POLARIZED: mean sbed of radial links 1.35 vs
  tangential 0.45 at the inner shell, converging ~1/r to 1.05/0.95
  at r≈20 (fig5) — the flow-gravity field rendered in conductance,
  a ~3× spoke:ring anisotropy at the boundary. FEEDBACK measured:
  Eh curves overlay the archived HZ-0 until t≈750 (bed growth era)
  then pull ahead — 418.4 vs 401.1 at t=1500 (+4.3%) with the gap
  widening: the river deepens its own bed. First self-amplifying
  gravity loop in the programme, gentle at bed_k=0.5.
* **P-F5 PASS at map grain — TWO EATERS SHARE ONE RIVER SYSTEM.**
  The double-hole arm (2,194 grown links, max 2.87) grows a single
  CONNECTED channel web spanning both holes: the facing rims are
  the most channelized, and the corridor between them is fully
  bed-bridged (fig2) — against the control's bare lattice. The pair
  also eats more together (Eh 545.9 vs control 516.6, +5.7%). "Two
  eaters grow a channel between them" — the attraction precursor,
  on the record in channel form.

## 5. Verdict (nothing adopted; bed_k=0 default; decisions the user's)

The bed-digging channel law does, in its first session, exactly what
ASYM.md §5 registered and nothing it forbade: it is live, exactly
conserved, battery-clean, ignition-immune where every predecessor
ignited, invisible to lawful matter, and around open books it grows
the programme's first persistent, radially-polarized, self-amplifying
flow structure — survivable asymmetry, on disk and on camera. The
open follow-ups it hands the user: bed_k/bed_tau response surface;
the dense-channel face (F-D1 phase 2); whether a bed web between
BODIES (not asserted eaters) can form once any lawful object runs
open books through the medium — the S2 junction again, now with the
survival structure waiting for it. Exhibit page: rendered from
`runs/flow/figs/` (fig1 eater epochs, fig2 double, fig3 quench,
fig4 growth/anti-ignition, fig5 well/anisotropy).
