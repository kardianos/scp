# DS — the double slit on the free-cell substrate (P1 first verification)

**Question:** does the free-cell fabric — cells whose geometry is state,
channels born and dying by the contact rule — carry the two-plane field
sector's interference exactly as the frozen graph did (v89
`DOUBLESLIT.md` §8, tier 0: fringes at parameter-free loci)?

**Kernel:** `kernel/freecell.c` (kernel of record), apparatus `exp=slit`.
**Logs:** `runs/ds_*.log`. **Analysis:** the `ds` experiment in
`cmd/battery` (bars below).

## Apparatus (declared before first run)

* **Medium:** a **2D** free-cell bath (dart throw in the z=L/2 plane, jam
  settled; coplanar dynamics is exact — all forces are in-plane). L=40,
  dmin=1 → ~1300 cells. 2D because the v89 screen convention is I(y) on a
  line, and because it makes the campaign cheap enough to run controls.
* **Optics regime, declared law-key overrides** (the v89 optics precedent,
  `cellfab_runs/d1_baseline.log`: q_detune=0 e_cond=99): **`q_detune=0`**
  (linear limit made exact) and **`e_cond=99`** (no condensation anywhere —
  tier 0 is the *wave* experiment; clicks are a later tier). Everything
  else laws_V2g verbatim. These are printed in every log header.
* **Wall = carved vacuum** (free-cell-native: transport and hops need
  contact; a slab with no cells is a wall — "matter is emergently opaque"
  in reverse). Thin (v89 lesson: thick walls cut off): wall_x=16,
  thickness 1.6. Slits: two kept bridges at y = 20 ± 4.5 (sep s=9), full
  width 4 (v89 lesson: ≥0.7λ; 4 = 0.70·5.71).
* **Source:** field packet (conservative — no drive), Gaussian
  σx=3, σy=10 (quasi-plane front), centred x=8, carrier tilt kx=1.1 →
  **λ = 2π/kx = 5.712 seeded** (the v89 DS wavelength).
* **Screen:** meter line at x=28 (D=12 behind the wall, the v89 distance):
  time-integrated Ee per y-bin (80 bins), gated t ∈ [16, 60] (packet
  arrival ~35 t.u. at the measured pulse speed 0.574·C; the periodic box
  wraps a re-arrival only after ~70 t.u. — the gate closes first; the box
  stays exactly conservative, no absorber needed at tier 0).
* **Controls:** slit A only (mask=1), slit B only (mask=2), no wall
  (mask=3, calibration: packet speed + on-medium wavelength), frozen
  substrate (freeze_geo=1) beside the live default.

## Parameter-free loci

Path difference δ(y) = √(D²+(y−y_A)²) − √(D²+(y−y_B)²) with y_A=15.5,
y_B=24.5, D=12, λ=5.712 seeded: **maxima at δ=mλ** → y = 20 (m=0),
y ≈ 12.5 / 27.5 (m=∓1, exact 2D formula); **minima at δ=±λ/2** →
y ≈ 16.3 / 23.7. Small-angle fringe spacing λD/s = 7.6.

## Pre-registered predictions (before first run)

* **P-DS1 (the claim):** I_AB(y) on the LIVE substrate shows the two-slit
  structure — a maximum within ±1.9 (Δy/4) of y=20, minima within ±1.9 of
  the predicted minima, raw central visibility V = (I_max−I_min)/(I_max+I_min)
  ≥ 0.15 (v89 frozen-graph V_norm was 0.32; the free foam is rougher).
* **P-DS2:** single-slit envelopes I_A, I_B are fringeless (no extremum
  alternation at the fringe period) and peak roughly behind their slits.
* **P-DS3 (the additivity breaker):** r(y) = I_AB/(I_A+I_B) < 0.8 at the
  predicted minima and > 1.2 at the central maximum — the round-1
  falsifier (I_AB = I_A+I_B) must NOT hold.
* **P-DS4:** conservation at the FP floor through the whole run (the
  periodic box pays for the no-absorber choice; nothing leaves the ledger).
* **P-DS5 (the free-cell claim):** the LIVE substrate reproduces the
  frozen-scaffold pattern — same maxima/minima loci within one bin (0.5),
  V_live within a factor 2 of V_frozen.

## Instrument findings on the way (kept for the record)

1. **A thin vacuum wall is transparent.** The contact rule links across any
   gap thinner than cfac·(rᵢ+rⱼ) ≈ 2.07 — the first wall (th=1.6) blocked
   nothing (all four masks: equal exposure). A vacuum wall must be thicker
   than the candidate cutoff: th=2.4.
2. **A vacuum gap in the LIVE foam heals in ~15 t.u.** Measured from `.fcs`
   snapshots: 16 slab cells at t=0 (the slit bridges), 68 by t=15, spanning
   every y — the jammed packing creeps into the free volume and the wave
   crosses at full strength. *The substrate defends no shape that no energy
   defends* — consistent with FREECELL's boundary thesis, and worth having
   measured. Consequence: the wall is a **pinned fixture** (`slit_pinw`:
   pass D skips cells within ±3 of the wall plane), the free-cell analog of
   the v89 imposed detuned-mute wall. The medium everywhere else is live.
3. **The periodic box has a back door.** The packet's backward/reflected
   components wrap in x and reach the screen from behind inside the gate —
   mask-independent exposure was the symptom (147–154 vs 180 no-wall).
   Cure: L=64 with the geometry placed so the back-door path arrives after
   the gate closes (no absorber needed; the box stays exactly conservative).
4. **2D group speed ≈ 1.05·C at kx=0.9** (the 3D bath gave 0.574 at kx=1.1)
   — the transit gate is timed from the measured 2D speed.

## Results (2026-08-03; logs `runs/ds_*.log`, panel `runs/ds/`)

Transmission sanity: no-wall 199.4; both slits 72.0; A-only 39.2; B-only
37.0 (I_A+I_B = 76.2 ≈ I_AB = 72.0 — total flux adds; the *pattern* is
what interference redistributes).

**Tier 0 confirmed on the live substrate — all five predictions:**

| run | y_peak (pred 32.0) | r(max0) | r(min−) | r(min+) | V_r |
|---|---|---|---|---|---|
| live, seed 20260802 | 33.10 | 1.67 | 0.37 | 0.34 | **0.652** |
| live, seed 111 | 32.30 | 1.59 | 0.47 | 0.29 | 0.614 |
| live, seed 314159 | 32.20 | 1.67 | 0.56 | 0.71 | 0.448 |
| **frozen control** | 32.60 | 1.51 | 0.30 | 0.30 | 0.667 |

* P-DS1 ✓ — central maximum at the parameter-free locus on every seed
  (|Δy| ≤ 1.1 < 1.9); V_r 0.45–0.65 ≥ 0.15. (v89 frozen-graph reference:
  V_norm 0.316.)
* P-DS2 ✓ — single-slit profiles are smooth two-hump envelopes.
* P-DS3 ✓ — additivity broken both ways: r ≤ 0.71 at every predicted
  minimum, r ≥ 1.51 at the central maximum, side maxima rise at the m=±1
  loci (r ≈ 1.35 at y≈18 and y≈46, predicted 17.9/46.1).
* P-DS4 ✓ — conservation −2.2e-15 through every run (the periodic box
  pays for the absence of an absorber).
* P-DS5 ✓ — **the live substrate carries interference as well as the
  frozen scaffold**: V 0.652 vs 0.667 (ratio 0.98), loci within 0.5.
  Substrate-freeing does not own the fringes — channel birth/death under
  the wave (live topology churn) does not decohere it.

**Gate:** battery experiment `ds`, 8 bars, GREEN (runs/BATTERY.log — 40
bars total). Field-transit images: `runs/ds/viz/` (volview, view=ee).
`ds_phase_t25.png` (FCS v2 phase view) shows the **mechanism**, not just
the screen result: two circular wavefront systems radiating from the slit
exits — Huygens wavelets at period λ — overlapping toward the screen.
`interactive_selfcheck.png` = the `volview -i` orbit view of the same
apparatus (wall slab + slit bridges + link graph).

## Tier 1 — first pass (2026-08-03): clicks avoid the minima

Apparatus: the condensing screen (`slit_clicks=1`): screen-strip cells
become condensation-active (per-cell e_cond override → 0) while the open
region keeps the optics regime; a **click** = one quantized conversion
event there (the carried atoms machinery fires in whole grains). Beam
amp=2.0, T=60; logs `runs/ds/clicks_two.log`, `clicks_A.log`.

* Two slits: **37 clicks**, and the predicted minima bands are EMPTY —
  0 clicks in y ∈ [26,30] and [36,40] — while clusters sit on the
  predicted maxima (5 near y=32, 6 near y=41.8).
* Single slit (control): **37 clicks**, and the same bands FILL
  (6 in [28,32], 2 in [36,40]) — no interference structure.

The tier-1 discriminator (clicks rebuild the wave's loci; which-slit
truncation removes them) is present at first-pass statistics. Honest
scope: ~37 grains per run is far below v89's 223; tier 1 is NOT yet
claimed — it needs a multi-seed grain harvest (≥200 clicks) and a
quantitative loci test before joining the battery.

## Tier 1 — instrument findings on the way to the harvest (2026-08-03)

5. **The first-pass single-slit control was reverb-lit.** Every one of
   its 37 clicks arrived at t ≥ 44.5 — after the exposure gate. At
   amp=2 a single slit never crosses the click threshold
   (0.25·Ee ≥ ε ⇒ Ee ≥ 1.21) during direct transit; the conservative
   periodic box then fills with wrapped field until the screen strip
   clicks on reverb. The first-pass "bands fill" observation stands
   only at reverb time; the harvest control must click IN-GATE —
   which sets the beam amplitude.
6. **The amp window is 3 ≤ amp ≤ 4 (calibrated, seed 20260802).**
   amp=2: two-slit clicks in-gate, single-slit zero. amp=3: controls
   click barely (2/5 per mask). amp=4: two-slit 20 in-gate with minima
   bands still EMPTY, controls 7/12 with bands FILLED — and evap = 0
   everywhere (optics linear). amp=5/6: cap evaporation fires on the
   hot central fringe (evap 0.53/1.06) and the two-slit minima begin
   to fill (threshold catches sub-fringe light): the discriminator
   erodes exactly where the declared linear regime ends. **amp=4 is
   the harvest point.**
7. **Hot beams exposed a kernel overflow (fixed, battery-gated).** At
   amp ≥ 4 the pass-7 alignment weight `kappa_align·dt/fs` overflows
   to inf when a cell's only flux is subnormal dust (fs ~ 1e-311),
   ending in NaN positions (the kernel's FATAL tripwire caught it).
   The bounded quantity is the flux-weighted mean direction |ax/fs| ≤ 2;
   the fix computes `kappa_align·dt·(ax/fs)` — same math, stable order,
   both kernels (C + Go step.go) identically. Full battery re-run:
   40 bars GREEN; dense-transport logs reseed at the ulp level
   (pair/ring/blob/uud), all measured values unchanged.

## Tier 1 — the grain harvest (pre-registered before the panel ran)

**Panel:** 12 two-slit seeds (20260802, 111, 314159, 271828, 141421,
173205, 577215, 662607, 299792, 137035, 161803, 618033) + the first 6
seeds × both single-slit masks (12 control runs). Tier-0 apparatus
verbatim except `slit_clicks=1 amp=4`. **In-gate clicks only**
(t ∈ [20,36], the tier-0 exposure gate — finding 5 makes the gate
load-bearing at tier 1).

**The quantitative loci test — the fringe-phase score.** Each click at
y scores s = cos²(π·δ(y)/λ), δ from the parameter-free loci (yA=27,
yB=37, D=14, λ=2π/0.9). A click distribution blind to fringe phase has
E[s] = 0.5 exactly; clicks sampling the two-slit intensity are
phase-locked (s̄ → high). Minima bands = predicted minima ±1.5
([25,28] ∪ [36,39] — where ideal fringe intensity < 15% of peak);
central max band [30,34].

**Pre-registered bars (battery experiment `ds1`, 8 bars):**
n_two ≥ 200; s̄_two ≥ 0.62 (≥ 7σ above the phase-blind null at n=200);
s̄_ctl ≤ 0.55 (controls ARE phase-blind); s̄_two − s̄_ctl ≥ 0.15;
minima dark two-slit: f_min_two ≤ 0.5·f_min_ctl; control fill:
nmin_ctl ≥ 8; R3 atomicity: every click e = k·ε(w₁), k ≥ 1 whole
(tolerance 2e-5 = the log's %.5f print floor); panel |drift| ≤ 1e-13.
Calibration expectation (12-seed yield check, two-slit only): n = 270,
s̄ = 0.683 ± 0.017, f_min = 0.048.

## Tier 1 — CLAIMED (2026-08-03; battery `ds1`, 8 bars GREEN)

First panel run, all eight pre-registered bars pass; floors then
sharpened to the measurement (ratchet: bars encode claims):

| statistic | measured | sharpened bar |
|---|---|---|
| grains, two-slit in-gate (12 seeds) | **270** | ≥ 200 |
| s̄_two (fringe-phase score; null 0.5) | **0.6827 ± 0.0168** (10.9σ) | ≥ 0.62 |
| s̄_ctl (12 single-slit runs, n=92) | **0.3803 ± 0.0358** | ≤ 0.50 |
| separation s̄_two − s̄_ctl | **0.3024** | ≥ 0.18 |
| minima-band fraction two-slit vs ctl | **0.048 vs 0.370** (7.7×) | ratio ≤ 0.33 |
| control minima-band clicks | **34** | ≥ 16 |
| atomicity | every click k·ε, k ∈ {1,2}, dev 3.5e-6 (print) | ≤ 2e-5 |
| worst panel drift | **7.1e-15** | ≤ 1e-13 |

Two readings beyond the pass: (1) the controls are not merely
phase-blind — s̄_ctl sits 3.3σ BELOW 0.5, because the single-slit
envelope peaks behind its slit, exactly where the two-slit minima lie;
which-slit truncation doesn't just remove the fringe preference, it
puts the clicks where the fringes forbid them. (2) Grains come in
whole atoms only: k=1 overwhelmingly, k=2 occasionally, nothing
fractional — the R3 claim rule (a click = a whole conversion grain)
holds at n=362 across 24 substrate realizations.

**The tier-1 statement:** on the live free-cell substrate, single
quantized conversion events at a condensing screen rebuild the wave's
parameter-free interference loci (fringe-phase-locked, minima dark),
and which-slit truncation removes the structure — at 362 grains,
gated. Images: `runs/ds/viz/ds1_t25_phase.png` (Huygens wavelets at
the harvest amp), `ds1_t25_ee.png`/`ds1_t30_ee.png` (transit),
`ds1_avg_ee.png` (time-average). Stream: `runs/streams/ds1_two.fcs`.

## What later tiers need

Eraser/delayed-choice analogs after tier 1 is green. The Go kernel does
not implement `exp=slit`/`exp=rings`/`exp=blob2` yet (C-only apparatus;
the A/B cross bars cover the shared experiments).
