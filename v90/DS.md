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

## What tier 1+ needs (not run yet)

Single-quantum clicks = a *condensing screen* (localized conversion: the
screen strip's cells given e_cond=0 while the open region keeps the optics
override) riding the carried atoms machinery; then the eraser/delayed-choice
analogs. The Go kernel does not implement `exp=slit` yet (C-only apparatus;
the A/B cross bars cover the shared experiments).
