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

## Second execution (2026-08-03, same day): CQ3b/CQ7/CQ8/CQ9

**CQ3b — compensation works; the facing-only variant refuted.**
Facing-cell-only compensation measurably BROKE the fifth (fifth-facing
cells carry no standing flight; shifting them un-tunes the arithmetic —
recorded, apparatus corrected). The dynamic diagnosis via per-ring xbar:
U-rings run +0.0093 above seed (internal flight prices into pitch —
flight-loads-pitch, as the law says); the D-droplet runs −0.0009. With
**per-flavor whole-ring compensation** (U seeds at x−0.0093), the nv=12
fifth-defect drift falls **8×**: −0.0646 → +0.0078 rad/t.u. (slightly
over-compensated — the residual is the transient, below).

**CQ7 — H-stiffness NOT confirmed; the real requirement identified.**
Drift rates: nv=12 −0.065, nv=24 −0.108 raw; compensated: +0.008 vs
−0.094. Drift GROWS with N, and static compensation fails at nv=24
because its flight excess is transient (final xbar ≈ nominal — the
interior decoheres and its flight dies, chirping the pitch on the way
down). **Measured conclusion: a stationary composite pitch requires
steady-state internal circulation — the flux-machine interior** (the
v89 particle thesis, now forced at composite scale by measurement).
Bulk alone does not average noise away; bulk with dying transport makes
it worse.

**CQ8 — the bath PROTECTS the composite.** nv=12 UUD embedded in the
live bath (L=20, carve 392): the UU boundary bond is the best measured
yet — gg = 0.9902 with the bath *compressing* it onto the rung
(dev = **−0.026** vs vacuum's +0.11 — confinement compression, same as
FREECELL §9.4), and ring-interior coherence is 2–5× better than vacuum
(gg_int 0.42 vs 0.09–0.23): cage pressure keeps internal channels
transporting. Fifths remain dead (as in vacuum). The composite's home
is the bath; the fifth's home is the amplitude completion.

**CQ9 — the residual force at nucleon scale.** Two UUD nucleons
(kind 3): inter-nucleon channels form only below COM sep ≈ 9–10
(ncross 1–2, dmin ≈ 1.5 = contact), producing a weak push-out
(7→7.13, 8→8.06, 9→9.02); exactly inert at 11 (ncross 0, sep frozen to
4 decimals). The rate-level residual interaction is 1–2 contact
channels' plastic settle — no capture, no fusion, one more instance of
the universal contact law.

## Third execution (2026-08-03): the parity sweep and the CQ8 A/B

**The mod-12 sweep REFUTES the macro-parity rule as an extrapolating
law** (nv = 12/18/24/36, vacuum, T=200, `runs/composite/cqp_*.log`;
nv=18 additionally at all three ℤ₃ branches). No battery bar — the
candidate failed its sweep. What the four points measure:

1. **The seed-time cycle residual walks continuously, not by π-parity.**
   The closing fifth's t=0 defect: −1.105 (nv=6), −2.213 (12), +2.959
   (18), +1.845 (24), −0.393 (36) — a steady ≈ −1.11 rad per nv-step
   of 6. The recorded rule conflated the interior CYCLE sum (where a
   rung hop contributes π) with the boundary EDGE meter, which weighs
   interior phases by the ℤ₃ ratio (q=3/p=2, mod 2π) — the hop parity
   is invisible to it. P-CQ2's arithmetic reading dies here.
2. **nv=18 (odd hops, predicted frustrated) is NOT nv=6.** Its fifths
   precess and settle near fit (−0.20/−0.97), not pinned at ±2.5–2.9;
   and its endings are BRANCH-DEPENDENT (k=1 keeps the UU bond at
   gg 0.986; k=0/2 lose it) — plastic-lock history, not arithmetic.
3. **nv ≥ 24 loses the composite boundary itself at T=200**: facing
   edges slide beyond contact, channels die (fifths CHANNEL DEAD, UU
   dead), interiors decohere — CQ7's transient-chirp mechanism eats
   the whole object. The committed cq7 runs show the same endings;
   the sweep just makes it size-systematic.
4. What survives, sharpened: **nv=6 is uniquely π-pinned** (fifths at
   max misfit ±2.5–2.9, gg ≡ 0, branch-identical — reproduced), and
   nv=12 is the largest size whose dead fifths rest near fit with the
   composite intact at T=200. The 6↔12 contrast is real; the CLAIMED
   quantization law behind it is not. The small-nv discriminator is
   an open question parked with the flux-machine frontier (a decaying
   interior has no stationary arithmetic to be quantized BY).

**CQ8 A/B — in the bath, compensation HURTS.** Same bath apparatus,
xcomp 0.0093 vs 0 (`cq8_bath.log` vs `cq8_noc.log`, T=120):

| meter | bath + comp | bath, no comp |
|---|---|---|
| UU bond gg / ψ | 0.9902 / +0.047 | **0.9980 / −0.017** |
| bond dev | −0.026 | +0.022 |
| bond standing flight lem | 5.0e-3 | **5.8e-2 (12×)** |
| interior gg_int (U,U) | 0.42 / 0.47 | **0.59 / 0.70** |
| D droplet gg_int | 0.0032 | 0.0116 |
| fifths | dead | dead |

The vacuum-derived offset (U rings −0.0093, tuned to vacuum standing
flight) OVER-corrects in the cage: the bath sets its own operating
point (bath xbar 0.25–0.26 vs vacuum ~0.27), and the mistuned interior
weakens every lock. **Compensation is a vacuum instrument; the bath
retunes as it protects.** Static pitch trims do not transfer between
environments — the flux-machine requirement again, from a fourth
direction.

**CQ8-long (T=480) — protection does NOT persist: the composite
DISSOLVES into its protector.** The UU bond stays locked through
t ≈ 144 (gg 0.991 at the t=144 sample) and is dead by t ≈ 192
(ψ jumps to +1.8, dev +0.165, lem → 1e-9). The cause is in the
ledger, not the lock: **the composite's dense load diffuses into the
bath** — U-ring xbar 0.251 → 0.16/0.18 and Em −36%, D droplet xbar
0.617 → 0.384 and Em −39% by t=480. As the carriers' load drains,
their pitch rises and the rung arithmetic moves out from under the
(plastically frozen) geometry — the bond dies of detune, downstream
of load leakage. Images: `cq8_long_t20_em.png` (three rings sharp
against a dark bath) vs `cq8_long_t480_em.png` (the bath glows with
the leaked load; the composite barely exceeds background);
`cq8_long_avg_ee.png` (time-averaged field — shell probe, diffuse).
**The bath protects the locks on the timescale it takes to rob the
load** (cage compression and load equalization are the same contact
physics). Static composites have no equilibrium in an open medium —
the flux-machine interior is forced from a FIFTH direction: only
steady throughput could hold xbar against the environment.

## Ladder state

* Executed: CQ0–CQ9 first passes; CQ3/CQ7 stiffness verdict; the
  mod-12 parity sweep (rule refuted, no bar); CQ8 A/B (bath beats
  bath+comp); CQ8-long (protection transient — dissolution measured).
  Battery through every kernel change: ALL GREEN (64).
* **Next frontier (all lanes):** the steady-interior composite — an
  internally circulating (flux-machine) composite whose pitch is
  stationary by throughput, not by seed. Rate-level candidates: driven
  boundary conditions are excluded (drives carry scales); the honest
  route is the S2 amplitude completion where transport IS the current.
  Five independent measurements now converge on it: H-stiffness
  transient (CQ7), compensation transience (CQ3b), bath retuning
  (CQ8 A/B), bath dissolution (CQ8-long), and the v89 particle thesis.
