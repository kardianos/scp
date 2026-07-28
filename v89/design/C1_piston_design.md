# C1' — Es reservoir / piston (design note)

*2026-07-28. Produced by a parallel design agent (v89-isolated,
read-only); archived verbatim-in-substance. Status: design accepted,
implementation pending (apparatus code behind piston_m>0, ratchet-gated).
MASS tech tree C1'.*

## Mechanism — cflag 6 "piston" skin + explicit reservoir account R

New instrument class cflag 6: every cell within `piston_m` of a box face.
The piston never touches an interior cell — no-suction is preserved
because the only exchange with the box is pass S's own pressure-pushed,
outflow-limited, exactly-paired link moves between skin cells and free
neighbours.

1. Init: Es[p] = pin(0) BEFORE E0_total is computed. Reservoirs
   R_s = R_f = 0.
2. Each step, immediately after pass S's gather-apply: for each piston
   cell, δ_p = Es[p] − pin(t); R_s += δ_p; Es[p] := pin(t). Pass S
   itself untouched. Raising pin = pressure source (outflow still
   limited by the skin's own avail); lowering pin = low-pressure
   destination. Interior outflow still originates from the interior
   side's own store — pressure pushes, nothing reaches out.
3. Field at the skin: record-media drain identity, booked to R_f (not
   local Em) — impedance-matched absorber at room pitch; keeps the
   skin's displacement pressure exactly π = pin.
4. cflag ≠ 0 already excludes piston cells from dense exchange, beat
   conversion, census, free sums.

## Bookkeeping — drift stays at the floor

Invariant every diag row: tEs + tEm + tEe + tET + R_s + R_f = E0_total.
The reset delta is exact by the Sterbenz lemma (per-step skin deviation
≲3e-3 against pin ∈ [0.8,1.2]); R_s/R_f are Kahan-compensated serial
accumulators. A 1.0→1.2 anneal at L=24 moves ~1.6e3 through R_s against
E0 ~ 1e4 — a 16% ledger line, which is why R is a COLUMN, not a
footnote. New diag line: `# PISTON t pin npist R_s R_f dRs dRf`.

## Params (apparatus only; laws untouched)

piston_m (0=off; ~1.2 = one-cell skin), piston_es0 (−1→e_s0),
piston_es1 (−1→hold), piston_t0/t1 (linear ramp window; t1≤t0 = step),
piston_absorb (−1→absorb_rate). Guards: pin ∈ [es_floor+0.02, 1.5·e_s0]
— the candidate-link graph is built to 1.15×radii ⇒ Es ≤ 1.15³ = 1.52 or
touching pairs silently lack channels. Warn on es_gx/edge_sink overlap.

## Geometry at ambient 1.2 / 0.8 (laws_V2g)

The load law counts BOUND energy only and d, C are Es-free ⇒ **rung loci
and the skirt locus are ambient-invariant by construction** (x*(d) and
x_skirt = 0.0617 unchanged). What moves: radii cbrt(Es) 1.063/0.928 →
contact ceiling d_cut 1.807/1.578, deepest rung load +17%/−20%, mean
conductance +28%/−30%, loaded death radius 1.736/1.484. High pressure
EXTENDS the species ladder; low pressure truncates it from the top (at
0.8, loaded pairs d ∈ (1.48,1.55] strangle). Ring6 at d=1.25 is
geometry-safe across the band. **Tripwire: any measured locus shift or
ring closure/2π ≠ 3.000 under pressure falsifies the bound-energy-only
pitch law at that point.**

## Risk table (condensed)

1. Skirt kinetics scale with conductance ±30% → run e7p (e7 + piston
   hold 1.2/0.8) as scored experiments before censuses; e7 unchanged.
2. pin ≤ 1.52 hard (candidate-graph ceiling); warn ≥1.4.
3. Ladder-end motion must not be misread as species physics — overlay
   computed d_cut/d_death on censuses; use pin 0.9 for rarefaction arms.
4. **Pre-existing artifact found during design: edge_sink cells
   accumulate recorded Em which enters pressure via s_disp — old sinks
   slowly push space inward on long runs.** Piston books to R so
   π_skin = pin exactly; prefer the piston for T > 500; any kernel-wide
   fix goes through the ratchet.
5. Quasi-static bound: box equilibration ~900 t.u. at L=24 — ramps
   t1−t0 ≥ 1000 or you launch pressure fronts; `# PISTON dRs → steady`
   is the settle proof.
6. G4 shells overlapping the skin read piston flux — interpret shells
   r < L/2 − piston_m − 1.
7. ε-grain discriminator for free: if the atom is secretly the local
   space grain e_s·d̄/C, pin 1.2 misprices conversions 20% — check
   QATOM grid residuals at pin ≠ 1 (HBAR candidate test).
8. Foam chaos ±30% → ≥5 shared seeds per arm; inert gate first
   (pin ≡ e_s0 arm matches no-piston within scatter; battery
   byte-identical at piston_m=0).

## First two experiments

EXP-1 (C2'/M3) blob census under scheduled ambient: heavy blob, T=3000,
arms A0 (pin 1.0) / A+ (ramp→1.2 over [200,1200]) / A− (→0.8) × 5
seeds. Success: A+ late |dM/dt| ≤ control/2 beyond scatter with the M0
roughness share FALLING; A− must accelerate (symmetric response = the
knob stirs, not pressurizes). Commissioning check: absolute core dip
Es_core ≈ Es_amb − s_disp·Em_core (0.72/0.32 predicted at 1.2/0.8).

EXP-2 (C3' = R1×R3) ring6 annealed: T=5000, arms B0/B+ (→1.2)/B−
(→0.9, above the strangulation confound), same 5 seeds as A2. Success:
B+ reaches the sharpened C1 plateau (|dM/dt| ≤ 1e-4) where B0 does not
and B− dies faster; closure/2π stays 3.000 (tripwire); QATOM residuals
checked at pin 1.2.
