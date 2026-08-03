# COMPOSITE — many-celled quarks: the primary focus

**User directive (2026-08-03):** many-celled quarks are the primary focus
of the program, with other many-cell experiments validating design and
simulation; **Charge and EMF run in parallel with Mass**, fed by shared
experiments. This document is the standing plan + running record.

**The ontology being built toward:** a cell is a grain of fabric, not a
particle. A **quark = a bound many-cell object whose collective mode
carries the flavor pitch** (U or D). A **nucleon = three such objects
bound by the incomplete harmonic** (the fifth triangle). Today's 3-cell
FCQ triangle is the minimal pattern-carrier; this program grows each
vertex into a composite and asks what survives, what sharpens, and what
breaks.

## Theory (stated before running)

**T1 — collective pitch.** A unison ring of N voices at per-cell load x
has every cell at w_e = w2/(1+q_detune·x); internally rung-locked, it is
one object with pitch w_e — *the same pitch a single cell at load x has*.
Composites change the substrate, not the arithmetic.

**T2 — the ladder is channel-level, so composite structure is boundary
arithmetic.** Gates live on channels between cells; the pair ladder
(q·wᵢ+p·wⱼ)d/C = 2πm constrains the **boundary-cell separations**, not
COM distances. Predicted inter-composite gaps are therefore the *same
numbers* as the minimal triangle (d*_UU(m=1) = π/w_U = d*_UD(m=2) at the
exact-fifth load; d*_DD = π/w_D beyond contact). A composite UUD is the
minimal triangle's arithmetic carried by facing boundary cells, with
bulk behind them.

**T3 — what many cells buy (the composite hypotheses).**
* **H-stiffness (the big one):** the measured killer of the fifth lock is
  flight-load pitch wander (the Γ/6 tongue, FREECELL §10.4). A collective
  pitch averages load fluctuation over N cells: σ(w̄) ~ σ_cell/√N. So
  composite flavor should be *stiffer*: the fifth-defect drift rate
  should FALL with composite size. If a critical N* exists where wander
  fits inside the tongue, **many-celled quarks would hold at rate level
  where single cells cannot** — the composite program's sharpest possible
  payoff, and testable by scaling N.
* **H-redundancy:** two composites can share several boundary channels;
  bond stiffness adds (~n_ch·K_b) and lock survival should improve.
* **H-mass-is-flight:** rate-level bonds are kinematic (P15 moves
  geometry with no energy entry), so naive binding energy = 0. The
  composite-specific energy component is the **standing flight on
  boundary channels** (in-flight inventory lem is real ledger energy).
  Predicted: E_total conserved exactly through binding; the "mass of the
  bond" is carried as cross-channel flight, measurable per channel.

**T4 — Charge lane.** The ℤ₃ branch structure lives on the interval
cycle; for composites the cycle runs over boundary channels. Branch
seeding, closure, and the antisymmetric fifth-defect precession
(minimal-triangle baseline: 0.016 rad/t.u.) should transfer — with
H-stiffness predicting *slower* precession at larger N.

**T5 — EMF lane.** Unison is silent (measured: exactly dark); detune
radiates (blob: rough 1.93/60 t.u.). A composite UUD's fifth boundary
channels detune dynamically → the composite nucleon-analog should have a
**measurable roughness emission rate and field halo** scaling with its
fifth-boundary count — the EMF observable of matter, produced by the same
runs that test Mass and Charge.

## The ladder (execute in order; all three lanes read the same runs)

| rung | experiment | lane(s) |
|---|---|---|
| CQ0 | `exp=rings` apparatus: R rings, per-ring loads, vertex-facing orientation, inter-ring gap offsets, boundary-cycle meters, per-ring coherence/energy/flight meters | — (battery-gated) |
| CQ1 | ring pair, unison, gap swept ±0.1 about π/w: boundary channel locks at the cell-pair rung? capture range? multi-channel redundancy? | MASS |
| CQ2 | composite UUD (ring_U, ring_U, ring_D; nv=6), branches k=0,1,2: ladder identity at boundaries, ℤ₃ closure, cycle circulation | CHARGE+MASS |
| CQ3 | precession stiffness: fifth-defect drift rate at nv=6 vs nv=12 vs the 3-cell baseline (0.016) — the H-stiffness test | CHARGE |
| CQ4 | composite UDD: D–D boundary beyond contact ⇒ open chain at composite scale | CHARGE |
| CQ5 | energy bookkeeping through binding: E_total at FP floor; cross-channel flight inventory = the bond's carried energy | MASS |
| CQ6 | emission: rough_total rate + Ee halo of composite UUD vs unison ring pair (silent control) | EMF |
| CQ7 | size scaling N = 6/12/24: pitch wander σ(w̄) vs 1/√N; extrapolate N* | all |
| CQ8 | composite in bath (live substrate): does bulk protect the boundary locks from bath scramble? | all |
| CQ9 | two composite nucleons: the inter-nucleon interaction (residual force analog) | MASS |

## Pre-registered predictions (before first run)

* **P-CQ1:** the unison ring pair's boundary channel locks on the
  cell-pair ladder (d* = π/w to the same precision as the 2-cell bond);
  beyond capture the pair is exactly inert (channels or nothing).
* **P-CQ2:** composite UUD closes: all boundary channels alive, ψ defects
  near 0 at seed for every ℤ₃ branch (closure is arithmetic, so exact at
  seed); UDD's D–D boundary forms no channel (open chain).
* **P-CQ3 (H-stiffness):** fifth-defect drift rate decreases from the
  3-cell baseline as N grows; direction is the claim, the exponent is
  the measurement.
* **P-CQ5:** E_total conserved at the FP floor through approach and
  locking; nonzero standing lem on locked boundary channels.
* **P-CQ6:** rough emission: UUD composite > 0; unison pair = exactly 0.

## Results (2026-08-03, first execution; logs `runs/composite/`)

**CQ0 ✓** `exp=rings` apparatus in the kernel of record (kinds pair/UUD/UDD,
per-ring molecule/droplet auto-selection, vertex-facing placement,
boundary-cycle meters). Battery after the kernel change: ALL GREEN (40).

**T6 confirmed in-run:** D-composites auto-select **DROPLET** (rung
2.171 > contact 1.529); U-composites are **MOLECULES** (rung 1.447 <
contact 1.647). The U/D structural asymmetry is real and printed in
every seed line.

**CQ1 ✓ (MASS) — the composite bond works; and it poisons its bulk.**
The unison boundary channel locks robustly: gg = 0.989, ψ = −0.053,
and all four gap offsets (±0.05, ±0.10) converge to the SAME attractor
(dev = +0.113 ± 0.001 above the rung — the composite pressure offset;
P-CQ1 confirmed, both sides). The bond carries **standing flight**
lem = 0.094 (pair) / 0.336 (UUD's UU edge) — H-mass-is-flight
confirmed: the bond's energy is in-transit inventory, and E_total holds
at the FP floor (1.7e-15) through binding (P-CQ5 ✓). **But the interior
pays:** every internal edge slid to exactly contact (dev +0.1997 =
contact − rung), lock surviving only adjacent to the facing vertices
(gg 0.38–0.40 there, ~0 elsewhere). Mechanism: the boundary channel's
standing flight load-detunes its carrier cells; the pitch discontinuity
breaks the interior unison chain. **Binding poisons bulk at rate
level** — the flight-load mechanism's third appearance.

**CQ2/CQ3 — THE MACRO-PARITY SELECTION RULE (new law, confirmed).**
Composite cycle closure includes interior path phases, and at the rung
each interior hop contributes exactly w·d = π. So ℤ₃ closure — automatic
for 3 cells — becomes a **parity constraint on composite geometry**:
the interior hop count between a composite's two facing vertices must be
EVEN. Measured: nv=6 (facing vertices adjacent → 1 hop, odd) is
π-frustrated — fifth boundaries at ψ ≈ ±2.5–2.9, gg = 0.0000,
identically for all three ℤ₃ branches. nv=12 (2 hops, even): ψ =
−0.05…−0.50, gg up ~100×. **A composite quark's internal geometry is
quantized** (for the 60° triangle arrangement: hops = nv/6, so
nv ≡ 0 mod 12). P-CQ2's "closure exact at seed" was wrong as stated —
and wrong in the most informative way possible.

**CQ3 (CHARGE) — composite precession exists; stiffness not yet.**
The nv=12 fifth defect precesses **coherently** through a full 2π in
~200 t.u. → 0.031 rad/t.u. (minimal-triangle baseline: 0.016). The
composite does NOT yet precess slower (P-CQ3 unconfirmed): the boundary
flight-load adds its own detune, dominating bulk averaging. Next rung
CQ3b: flight-compensated loads (seed x below target so pitch is exact
WITH the standing boundary flight) — the same correction FCQ needed.

**CQ4 ✓ (CHARGE)** — UDD at composite scale: the D–D boundary forms
**no channel** (CHANNEL DEAD; open chain transfers exactly), both
D-composites droplets.

**CQ6 ✗→sharpened (EMF)** — rough = 0.000000 for both UUD and the pair:
emission requires a LIVE fifth channel, and rate-level composite fifths
die before they transport. P-CQ6 refuted as stated; the sharpened claim:
**composite EMF emission is gated on fifth transport, i.e. on the
amplitude completion** — the EMF lane and the Charge lane share one
frontier, quantified.

**Standing synthesis after one execution:** the three lanes converge on
a single measured frontier — flight-load pitch drift vs tongue width —
now demonstrated at three scales (single-cell fifth, composite fifth,
composite interior). En route, the composite programme extracted three
scale-free structural laws that no amplitude completion will change:
the U/D molecule/droplet asymmetry (T6), the macro-parity quantization
of composite geometry (even interior hops; nv ≡ 0 mod 12 for 60°
triangles), and boundary-ladder transfer (the composite bond locks on
the cell ladder with a pressure offset, carrying its energy as flight).

## Ladder updates

* CQ0–CQ2, CQ4–CQ6 executed (above). CQ3 partially (precession measured;
  stiffness test blocked on flight compensation).
* **CQ3b (new, next):** flight-compensated composite loads, then re-test
  H-stiffness (drift rate vs N at nv = 12, 24).
* **CQ7 sharpened:** N-scaling must use nv ≡ 0 mod 12 (parity rule).
* CQ8 (bath), CQ9 (two nucleons) unchanged.
