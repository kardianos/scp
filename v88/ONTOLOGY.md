# v88 — Ontology and design direction: the cycle-trapped fabric

> **SUBORDINATE TO `v88/PRINCIPLE.md`.** That document is foundational:
> energy is never destroyed, space is a *mode* of energy, matter is
> converted space, and curvature follows from conservation. Anything here
> that assumes a background — a permanent index set with evolving contents
> — is superseded by it, including every instrument built so far.


**Date:** 2026-07-27 · Supersedes the Q-ball ontology (v66–v87) as the route to
emergent particles. Instruments: `v88/fabric_proc.c`, `fabric_mass.c`,
`fabric_trap.c`, `lattice_iso.c`, `fabric_pn.c`, `v86/n_battery/fabric_aniso.c`.

Tags: **[D]** derived · **[M]** measured · **[P]** postulated · **[C]** conjecture


---

## 0. CORRECTIONS LOG — read before anything else

This ontology has been revised four times by the user against measurement. Each
correction killed a concept that looked reasonable and was wrong. The superseded
concepts are named explicitly so they are not silently reintroduced.

### C1. Q-balls are not the emergent particle  *(supersedes v66–v87)*
**Old:** the stable gauged complex Q-ball is the nucleon template (Stage 1).
**Why it died:** zero bound internal modes for every ω ≥ 1.36; Q varies
continuously across the window with *every* sub-turn point stable, so the
stability criterion selects an **interval**; no species, no mass quantisation.
N8 then found winding rings at Q = 210 do not exist at all (n=1,2 drive ω above
m and disperse), removing the last discrete label. **[M]**
**Now:** a continuum PDE is intrinsically homogeneous and can only produce
continuous families. Discreteness must come from the fabric being discrete.

### C2. Discrete breathers are not the mechanism  *(supersedes `fabric_trap.c`)*
**Old:** model the cycle sector as a discrete nonlinear Schrödinger lattice;
localisation by amplitude pushing a mode outside the linear band.
**Why it died:** a breather is ONE mode with an amplitude. Its "discreteness" is
a threshold on a continuous dial — the Q-ball's failure in new clothes. It
*mimics* localisation without possessing structure that could label species.
Measured: it localises in **all** d = 1…7 with spread 0.25–0.31, never sharp.
**Now:** multi-dimensional harmonics. A cell internally multi-dimensional has a
harmonic *spectrum*; closure conditions between harmonics are **integer
relations**, not amplitude thresholds.

### C3. Energy does not hop; it is MEDIATED  *(supersedes the bilinear coupling
in `fabric_harmonic.c`)*
**Old:** cell↔cell bilinear hop `ε(z_i* z_j + c.c.)`, with detuning suppressing
transfer.
**Why it died — measured, `GROK_V88_SYMBOLIC.md`:** instantaneous transfer is
`A_sec = 2ε₁I = 0.075` **exactly, independent of detuning**, CV = 0.000, flat
fit. Detuning suppresses only the *time average*; on-resonance/far ratio 1.65 —
not a gate. With a direct cell–cell term the energy is already present the
instant the Hamiltonian is written; nothing can gate it afterwards. **[M]**
**Now:** cells couple **only to links**, never to each other. A link carries its
own harmonic; energy must be deposited into it, carried through its cycle, and
only then delivered. Transfer takes time and is gated **structurally** — a link
that cannot be driven resonantly never accepts the energy.
**Measured (`fabric_link.c`):** transfer falls **227×** across the detune window
(1.000 → 0.0044) against the bilinear hop's 1.65. **[M]**

### C4. c is emergent, not an input  **[D, given C3]**
The speed limit is *the rate at which cycles can couple*. Measured
`c ≈ 1.05·g·a`, constant over a 4× range in g; a detuned link drops it to
0.36–0.90 of that. Nothing puts a speed into the model. **[M]**

### C5. A link is not a scalar — it actuates in TWO PLANES along a chiral
harmonic  *(supersedes the scalar link of `fabric_link.c`)*
**Old:** one scalar oscillator per link.
**Why it died:** it cannot bind, and the reason is geometric — **in 1D a link
has no transverse plane**, so there is nowhere for a second actuation plane to
live, no handedness, and no handedness-dependent interaction. The linearity was
forced by the geometry the test was run in, not by the physics.
**Now:** a link along d̂ actuates in the two transverse directions with a ±90°
relative phase. `u± = (u₁ ∓ i u₂)/√2` are left/right circular and
`2 Im(u₁* u₂)` is a signed handedness carried by the geometry.
**Measured (`fabric_chiral.c`):** chirality 2.28×10⁻¹ at d=3 against
8.4×10⁻³ / 1.2×10⁻² at d=1,2 (20× suppression); same- vs opposite-handed pair
energies **split by 0.364**. A scalar link produces no such split, which is why
`fabric_harmonic` measured pure repulsion. **[M]**

### C6. "3–7" is NOT a range of workable dimensions — it is 3 EMERGENT
DIMENSIONS FROM 7 DEGREES OF FREEDOM
**My misreading (recorded because it wasted work):** I and the grok seat both
read "only workable in 3–7 dims" as *the mechanism functions for d = 3…7* and
tested that span. Both found no threshold — `fabric_trap` localises in all
d = 1…7. **That was answering a question the user never asked.** **[M, but
irrelevant]**
**Actual structure:** there are **7 degrees of freedom**, and **3 dimensions
emerge** from them. The connecting fact is that handedness requires a cross
product, and **the cross product exists only in 3 and 7 dimensions**. The 7
imaginary octonionic units carry the degrees of freedom; the 3 emergent
dimensions are what the chiral structure realises. This is the project's v59
octonionic thread meeting the fabric ontology. **[P — algebraic, not yet derived
from any simulation here]**
**Consequence:** the dimension count is not a tunable parameter to scan. Scans
over d are measuring the wrong thing.

### What is NOT established
* The chiral pair-energy split **saturates** (0.3311 at sep 6 = sep 8) instead
  of vanishing at large separation. A constant offset at large separation is a
  **self-energy difference between seed configurations, not an interaction**.
  Only the sep 2→6 variation (~0.033, 10% of the split) is separation-dependent.
  **Handedness changing the cost of a state is measured; handedness producing a
  binding force is NOT.**
* d = 1,2 chirality is 20× suppressed but **not identically zero**, so C5's
  geometric-necessity claim is demonstrated only as a strong suppression.
* C6 is algebraic and postulated. Nothing simulated here derives 3 from 7.

---

## 1. What exists

**The fabric is discrete.** Space is not a continuum sampled on a grid for
convenience; the cells are the substance. The lattice spacing is a physical
constant, not a numerical parameter. This is the single change from which
everything else follows.

**A cell has internal structure and its own cycle.** A cell is not a point
holding a number. It is multi-dimensional internally, and it carries a cyclic
internal process. The cycle is what maintains the cell's size.

**Cells have variable size.** The fabric's density is a field. In the mapping

```
    x(i) = a·i + ξ(i)                    ξ is the fabric's own degree of freedom
    θ(i) = ∇·ξ                           local density; θ < 0 = shrunk
```

the displacement ξ that defines the tessellation is the same object that carries
mass. There is no field living *on* the fabric — the fabric's geometry is the
field.

**Energy crosses between cells only in complete cycles, and it is MEDIATED.**
This is the governing rule and the source of every discrete quantity. Cells do
NOT exchange energy directly (see C3 — a bilinear cell–cell term makes the
transfer instantaneous and ungateable). Every link carries its own harmonic;
energy is deposited into the link, carried through its cycle, and only then
delivered. Transfer therefore takes **time**, and the rate at which cycles can
couple is what **c** is (C4).

**A link actuates in two transverse planes along a chiral harmonic** (C5). Its
handedness `2 Im(u₁* u₂)` is a signed label carried by the geometry. This
requires a cross product, which exists only in **3 and 7 dimensions** — the
structural origin of the dimension count (C6).

---

## 2. What a particle is

A particle is a **region whose cycle rate differs enough from its surroundings
that no complete cycle connects it to the exterior.** Energy inside cannot
leave: it is reflected internally at the boundary. The object is held together
by that mismatch, not by a balance of forces.

Three consequences, all structural rather than fitted:

1. **Mass is a cell count.** A particle occupies a definite number of cells.
   That integer is a genuine discrete label, which no continuum soliton
   possesses.
2. **Existence is a threshold, not a dial.** The trapping condition is that the
   region's rate lies *outside* the band its neighbours can support. A state
   either closes its cycle or does not exist. There is no continuous family.
3. **Species are configurations.** Different particles are different internal
   arrangements satisfying the closure condition — not different points on a
   branch.

**The interior configuration is what the "Higgs" is.** It is an arrangement of
degrees of freedom the fabric already has, holding the dense interior against
its own compression. It is emergent and intra-particle. **It is never a field in
the kernel.** A put-in-by-hand Higgs scalar with its own potential is an
imported species and is forbidden by the programme's own rule — the v71 sector
(`higgs_v`, `higgs_lam`, `higgs_kap`, described in the source as a
"self-compression bag") is therefore a design defect and must not be used.

**Cyclic energy tightens the cell.** Opposite to a fluid, where internal energy
expands. More cyclic energy → tighter cell, up to a point; past that point the
energy cannot go into further tightening and must be released into
configuration. **[M]** — measured on the existing branch: r_half falls
29.296 → 2.305 (12.7×) as ω rises 1.32 → 1.4955, and the "point" is the
Vakhitov–Kolokolov turn at ω ≈ 1.482, past which the object reorganises
(HC-6: r_core −10%, s_max +197%, binding +79%).

---

## 3. Why the Q-ball was the wrong object

Not a matter of tuning. Three independent measurements, all from this project:

| | measured |
|---|---|
| no internal spectrum | HC-1: **zero bound internal modes for every ω ≥ 1.36** — the whole production region |
| no mass quantisation | Q runs smoothly 144005 → 86.68 → 111 across the window; **every point below the turn is GSS-stable** |
| survivorship does not quantise | the stability criterion selects an **interval** (ω < 1.482), not a set of points |

A featureless droplet on a continuous family. Even granting Q ∈ ℤ from a later
quantisation, M(Q) is a smooth function of an integer — a dense tower, not
distinct species. And the last candidate for a discrete label inside the picture
is gone: **N8** found that winding rings at Q = 210 do not exist at all
(n=1 and n=2 both drive ω above m = 1.5 and disperse), so there is no E(n)
ladder.

Homogeneity is intrinsic to a continuum PDE, so continuous families are what it
must produce. The project ran at width/spacing ≈ 10–20 where the lattice is
invisible (D8b: pinning at R/E ≈ 1.2×10⁻⁸), which reconstituted that homogeneity
numerically. "Refine dx until converged" — the criterion gating most of v86 —
pushes *away* from the physics.

---

## 4. The fabric geometry

Measured constraints on what the discrete fabric may be:

* **Simple cubic is excluded.** Discreteness and crystallinity arrive in fixed
  proportion: the PN anisotropy is ~½–⅔ of the barrier at *every* width/spacing
  ratio. There is no regime where a cubic fabric is discrete but isotropic. **[M]**
* **Stencil tuning fixes the linear operator only.** BCC 8+6 with w₂/w₁ = 2/3
  (the Kelvin foam, whose truncated-octahedron cell has exactly those 14 faces)
  zeroes the 4th moment exactly: dispersion spread **1.49%** vs cubic's 22.69%.
  But nonlinear pinning anisotropy improves only ~1.6×, because pinning is set
  by the reciprocal lattice and stencil weights cannot move it. **[M]**
* **Procedural irregular geometry with regular computation works and is the
  chosen architecture.** `x(i) = a·i + A·δ(i mod M)` keeps a fixed logical
  stencil — memory layout and parallelism as cheap as a cubic grid — while the
  physical geometry is irregular. Per-site weights are *solved*
  (`Σⱼ wⱼ dⱼ⊗dⱼ = 2I` exactly), so consistency and 2nd-order isotropy hold at
  every cell however irregular. Every procedural fabric beats simple cubic on
  pinning anisotropy (best 1.6×), and M=13 wins both channels. **[M]**
* **δ is the replaceable core.** Frozen today. Everything downstream is derived
  from `x()` at build time, so making it dynamical — driven by local energy or
  regional strain — touches nothing else. That is the intended path to a fabric
  that responds to its own contents.

Still open: the Kelvin fabric has not been measured in the same statistic as the
procedural ones, so "Kelvin wins on isotropy" remains provisional; and a
**quasicrystalline (icosahedral)** fabric is the only construction that is both
deterministic and 5th-order isotropic (0.68% dispersion spread), untested.

---

## 5. Evidence for the mechanism

**Static density alone is not enough. [M]** `fabric_mass.c`: cell size dynamical,
double-well, with the geometric zero-sum rule `Σ ∇·ξ = 0`. It phase-separates
(40% shrunk) but **coarsens** — d=2 gives 78 lumps with sd/mean = 5.4, d=3 gives
2, d=4 gives a single blob of 26674 cells. No preferred size. The reason is that
it has no cycles: nothing sets a scale and nothing traps energy.

*(An earlier version of this file's potential had its minimum at θ = 0, so the
uniform state satisfied the zero-sum rule trivially and sat at the global
minimum — 0.00% shrunk in every dimension. The zero-sum rule forces lumpiness
only if the uniform state is disfavoured. Recorded because the error is easy to
repeat.)*

### 5a. BREATHERS ARE THE WRONG FRAMING — recorded as a correction

`fabric_trap.c` models the cycle sector as a discrete nonlinear Schrodinger
lattice, i.e. **discrete breathers**. That framing is rejected.

A breather is ONE mode with an amplitude. It is localised by nonlinearity and
quantised by **nothing** — its "discreteness" is a threshold on a continuous
dial, which is precisely the Q-ball's failure wearing different clothes. It
*mimics* the localisation without possessing the structure that would produce
species or mass ratios.

**The correct object is MULTI-DIMENSIONAL HARMONICS.** A cell that is internally
multi-dimensional (§1) has a harmonic *spectrum*, and closure/matching
conditions between harmonics are **integer relations**, not amplitude
thresholds. Integers enter by construction rather than by crossing a threshold.
That is where species, mass ratios and quantisation must come from. DNLS has no
internal harmonic structure at all.

The results below are retained as a record of what the *wrong* model produced,
not as evidence for the ontology.

**Cycles change the answer, in the wrong model. [M]** `fabric_trap.c`: one field ψ, resonant
neighbour transfer, cycle rate shifted by the cell's own amplitude
(ω = ω₀ − γ|ψ|²). From **random** initial conditions:

| d | sites | band top | peak \|ψ\|² | cells/lump | n_lumps | spread |
|---|---|---|---|---|---|---|
| 1 | 4096 | 4.00 | 1.5588 | 6.4 | 330 | **0.280** |
| 2 | 40000 | 8.00 | 1.5285 | 6.3 | 3257 | 0.293 |
| 3 | 39304 | 12.00 | 1.4251 | 6.1 | 3287 | 0.305 |
| 4 | 38416 | 16.00 | 1.7501 | 9.0 | 2204 | 0.262 |

Energy self-organises into thousands of localised objects of ~6–9 cells with a
**repeatable amplitude** (spread 26–31%), against the static model's 5.4.

**Dimension range: no threshold. [M]** Extended to d = 5,6,7 — cells/lump
8.7 / 9.3 / 11.3, spread 0.281 / 0.280 / 0.253, thousands of objects in each.
Combined with d = 1–4 the mechanism works in **every** dimension tested, with no
3–7 preference:

| d | 1 | 2 | 3 | 4 | 5 | 6 | 7 |
|---|---|---|---|---|---|---|---|
| cells/lump | 6.4 | 6.3 | 6.1 | 9.0 | 8.7 | 9.3 | 11.3 |
| spread | 0.280 | 0.293 | 0.305 | 0.262 | 0.281 | 0.280 | 0.253 |

This does not settle the 3–7 claim; it **locates** it. `fabric_trap` contains
only the cycle sector — no density/elastic coupling and no internal cell
dimensionality, which are the two ingredients the ontology says matter. The one
argument that gives a lower bound near 3 is elastic (a dilational inclusion's
strain falls as r^{1−d}, self-energy ∫r^{1−d}dr converging only for d > 2), and
that sector is absent here. **If the dimension restriction is real it lives in
the coupling between cycles and cell size, not in cycle-trapping alone.**
Note also that discrete NLS supports breathers in all d even though continuum
NLS collapses for d ≥ 2 — all-d behaviour is expected, not surprising.

**Three honest caveats on that table:**
1. This is the discrete nonlinear Schrödinger equation. Discrete breathers are
   established physics (MacKay–Aubry); nothing here is novel as mathematics.
   What is new is only its use as the fabric's particle mechanism.
2. Spread 0.26–0.31 is **narrow, not quantised**. A genuine quantum number would
   give a far sharper peak or discrete clusters of values.
3. `sd/mean` is the spread of **peak amplitude** at local maxima; `cells/lump`
   is an IPR-derived *average*, not a per-object count. The per-object size
   distribution has not been measured.
4. **Coarsening: checked, inconclusive.** Same initial condition to
   T = 200/400/800/1600 gives n_lumps = 3384 / 3257 / 2453 / 3592 — NON-monotonic.
   So it is not irreversible coarsening, but the ~30% wander means the peak
   threshold is too crude to establish a stable preferred size either. Moot in
   any case: see §5a.

---

## 6. Bell, and where 2√2 might come from

**[C], and explicitly speculative.** v87 B1 established that E(θ) = −cos θ ⟺ the
correlation is a *single angular harmonic* ⟺ S_max = 2√2, and that the √2 is an
ℓ¹/ℓ² ratio. Under this ontology the sinusoid is not a modelling choice: it is
the cell's internal cycle, and √2 is the peak-to-RMS ratio of that cycle. The
Tsirelson bound would then be a statement about how much correlation a cyclic
process can carry — a limit of cell size and cycle energy rather than an axiom.

This is **not** established and nothing in v87 tests it. What v87 *did* settle
stands unchanged and is the constraint any such account must respect: an offline
per-object dichotomic readout gives |S| ≤ 2 as an algebraic identity, so a real
test needs settings that are dynamical fabric degrees of freedom (B2).

---

## 7. Programme

### Immediate
0. **Separate self-energy from interaction in the chiral pair test.** The split
   saturates at large separation, so what is measured is a state-cost
   difference, not a force. Subtract the isolated-seed self-energy and keep only
   the separation-dependent part. **Until this is done there is no evidence of
   binding.** *Highest priority.*
0b. **Make the d=1,2 chirality identically zero** by projecting out degenerate
   directions rather than naming axes, so C5's geometric necessity is
   demonstrated rather than merely suggested.
0c. **Stop scanning over d.** Per C6 the dimension count is not a parameter.
1. **Formulate the multi-dimensional-harmonic cell** (§5a). What is a cell's
   state; what exactly is a "complete cycle" between two cells of different
   harmonic content; what is the closure condition; where do the integers enter.
   This replaces the breather line of work entirely. *Highest priority.*
2. **Does the 3–7 dimension restriction follow from the harmonic structure?**
   The cycle-only model shows no threshold at d = 1..7, so if the restriction is
   real it lives in the cycle↔cell-size coupling or the cell's internal
   dimensionality.
3. **Per-object size distribution**, not IPR averages — the actual test of
   whether mass is quantised, whatever the model.
4. **Kelvin fabric in the procedural statistic**, to settle §4's provisional
   ranking.

### Annealing — deriving larger particles of different types
The mechanism makes species = configurations, so the way to *find* species is to
anneal into them rather than to seed them by hand. Try, and compare what each
produces:

| schedule | intent |
|---|---|
| **slow thermal** — high effective noise, cooled over long time | ground-state configurations; the "lightest" species |
| **quench** — instantaneous cool from random | metastable / excited configurations, and the size distribution they freeze into |
| **cyclic (repeated heat/cool)** | shakes out defects; tests which configurations are genuinely stable vs merely reachable |
| **driven / pumped** — energy injected at a chosen rate | whether larger objects require sustained feeding, and where the growth threshold sits |
| **staged density** — anneal at successively higher fabric density | whether heavier species appear only above a density, i.e. a mass ladder from compression |
| **seeded coalescence** — anneal with existing lumps present | whether lumps bind into composite objects (the route to nucleon-like structures from cell-scale ones) |

Record for each: the spectrum of object sizes, whether it is discrete, and
whether the same species recur across seeds. **A species that recurs across
independent anneals is the first real evidence of a particle type.**

### Then
* **B2** — the dynamical-analyzer Bell test, now meaningful because the fabric
  has genuine internal cycles to act as choosers.
* Dynamical δ: let the fabric geometry respond to its own contents.
* Composite structure: whether trapped objects bind, and with what binding
  energy — the first contact with nuclei rather than with single particles.

---

## Appendix — closed questions from the superseded background models

Retained because they closed real questions, not because the models stand.
`fabric_trap`, `fabric_harmonic`, `fabric_mass`, `fabric_link`, `fabric_chiral`
and `fabric_proc` all carry a permanent index set and are therefore background
models, retired by `PRINCIPLE.md`.

**Lock-vector species were a truncation artifact. [M] — question closed.**
Rebuilding `fabric_harmonic` at m = 2, 3, 4 and running identical quenches:

| m | locks found | \|n·ω\| | θ |
|---|---|---|---|
| 2 | (2,−1) | 0.062 | −0.06…−0.08 |
| 3 | (2,−1,0), (−1,2,−1) | 0.063 | −0.10…−0.16 |
| 4 | (−2,0,2,−1) | 0.0000 | −0.65…−0.71 |

**No species recurs across m.** Each truncation yields its own set and no other,
and m = 4 is a different dynamical regime entirely (peak I of 22–52 against
0.2–0.7). The integer labels were counting the mode cutoff, not the fabric.

**Annealing found nothing to discover. [M]** Quench/thermal/detune × 3 seeds
produced no recurring species. Consistent with the symbolic back-tracing, which
had already killed annealing-as-discovery on this dynamics: with no preferred
size and no attractive channel there is nothing for an anneal to find.

**Why these failed, in hindsight:** all of them place structure *on* a stage and
evolve its values. Under the principle there is no stage — species are
isomorphism classes of closed integer locks in a combinatorial complex, not
occupancy patterns on immortal cells. The failures were not parameter choices.
