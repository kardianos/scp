# MASS — the stable particle (living lab notebook)

**Living document.** Need and design for obtaining mass in v89: a dense
pattern that persists at constant energy. Update with every campaign
entry; never rewrite history — append, date, and correct.

Subordinate to `PRINCIPLE.md`. Governed by the ratchet rule
(`battery/README.md`). Do not import pre-v89 constructions (shells,
Q-balls, profiles): candidates must be built from v89 laws alone.

---

## 1. Need — why mass is now THE bottleneck

Everything measured in the P/G campaigns terminates at the same missing
object:

* **The 1/r far field (G4, measured 2026-07-28).** A sealed leaking blob
  has NO steady space throughput — every far-shell space flux is
  d(mass)/dt bookkeeping (outward evaporation wind, inward refill), and
  extrapolating dM/dt → 0 all of it vanishes. A 1/r field requires
  steady space cycling AT CONSTANT MASS: internal conversion that takes
  up and lays back space continuously — "internal energy consuming c in
  E = mc²." That is a property our blobs do not have.
* **Radiation pressure (P2).** Absorbed light's momentum dies at the
  conversion door (~100× recoil deficit). Momentum-carrying matter needs
  wave-borne translation (S2-full), and testing it needs matter that
  survives being pushed.
* **Collisions are deferred** (user directive): two-body experiments
  need participants that outlive the encounter.
* **Equivalence-principle and free-fall tests** (g2, recorded): need
  masses that neither evaporate nor freeze during the fall.

**The measured leak (the enemy):** gaussian blobs lose mass at
dM/dt ≈ −0.13 (heavy, e3a class) to −0.54 (light) in early windows,
slowing but never stopping; a blob is a relaxing accident, not a
particle. Where the mass goes is itself unmeasured (M0 below).

## 2. What counts as a particle — certification bars

A candidate is a PARTICLE when, under the standing law table, unedited:

| bar | criterion |
|---|---|
| **C1 lifetime** | free dense mass constant: |dM/dt| ≤ 1e-3 sustained over T ≥ 1000 (vs ~0.1–0.5 today). *(Sharpened 2026-07-28 after ring6 met the letter while losing 38% of its mass: a PLATEAU also requires deceleration to |dM/dt| ≤ 1e-4 or extrapolated half-life ≥ 10⁴ — a slow steady death does not count.)* |
| **C2 self-force null** | centroid speed ≤ the e3a bar (no self-acceleration) |
| **C3 not a seed accident** | survives an ensemble of foam seeds (≥ 5), not one lucky foam |
| **C4 identity** | survives a perturbation (weak packet strike) and returns to its mass plateau |
| **C5 the prize** | measurable space throughput at constant mass (G4 instrument): the 1/r source. This is what NOTHING can currently attempt. |

When one passes C1–C4, it ratchets into the battery as the standing
particle; C5 becomes experiment g5.

## 3. Design — candidate routes

### R1. Pressure annealing (natural formation)

The space-transport law (laws_V2g: pressure pushes) gives a genuinely
new knob: ambient space pressure. A dense lump depresses its own Es
(displacement), shrinking its radii and channel areas — it
SELF-INSULATES (measured: e3a sealing improved 0.00174 → 0.00084 under
V2g). Higher ambient pressure deepens the contrast. Design:

* Apparatus: an Es reservoir — boundary cells with pinned space store
  (a "piston"), schedule the pinned value over the run (anneal:
  hot condensate → slow pressure rise → does the surviving census
  shift to longer-lived lumps?).
* Instruments: the lump census (below) + the g4 radial machinery.
* Risk to watch: pressure also compresses the vacuum skirt landscape —
  the e7 skirt boundary moves with ambient Es; keep the battery green.

### R2. Spin / winding (topological protection)

Leaking is untwisting: give the pattern something that must unwind
first. Two v89-native carriers:

* Field winding: the two-plane signed amplitude supports phase winding
  m·2π around a loop (a vortex of the field mode). Seeder: ring of
  cells, amplitude annulus, phase wound m times.
* Dense circulation: in-flight energy (lem) circulating around a closed
  cell loop — momentum-in-the-books around a ring; the flux-moment
  instrument sees it directly.

### R3. Ring locks (the fabric's chords — CONSONANCE open area A.3)

The standing pair is a two-voice lock; a RING is N voices in a closed
rung cycle: Σ ωᵢdᵢ/C = 2πm around the loop. Closure is the point: a
ring's conversion cycle has no free end to leak from, and its
circulating exchange is exactly the "internal energy consuming c" that
C5 needs. Narrower tongues (harder to form, per Tenney height) but the
payoff is a closed internal process. Design:

* Seeder: N cells at ladder-spaced separations d* (m = 1), pitches at
  the rung loads x*(d), phases seeded in the locked relation
  (generalize pair_seedlock around the loop), optional winding on top
  (R2 + R3 combined: circulating phase AND circulating flight).
* First scan: N ∈ {3, 4, 6}, m ∈ {1, 2}, d near the pair ladder rung;
  lifetime census vs the same-mass gaussian control.

### R4. What is EXCLUDED

Shell constructions, Q-ball profiles, imported potentials — pre-v89
material. If a route needs a new field or species, it is dead on
arrival (standing constraint 3).

## 4. Instruments needed (build before the first scan)

* **Lump census** (apparatus-gated, like rad_diag): connected components
  of free Em above a threshold; per-lump mass, centroid, radius;
  `# LUMP t id mass ...` rows + a final lifetime table. This is the
  campaign's primary eye.
* **Leak taxonomy** (M0): per-window ledger attribution of a blob's
  mass loss: dense spreading (Em leaving the lump but staying dense)
  vs evaporation (D→F conversions) vs space return (back_s). Mostly
  derivable from existing RAD columns + conversion counters; add what's
  missing.
* Ring seeder (R3) and winding seeder (R2) as new apparatus init modes.

## 5. Experiment ladder

| id | question | needs |
|---|---|---|
| M0 | where does a blob's mass GO (taxonomy of the leak) | leak instruments |
| M1 | ring locks: lifetime vs (N, m, d), vs gaussian control | ring seeder + census |
| M2 | winding: does m ≠ 0 slow the leak? vortex lifetime vs m | winding seeder + census |
| M3 | pressure annealing: census vs pressure schedule | Es reservoir + census |
| M4 | certification C1–C4 on the best candidate; ensemble | all |
| g5 | throughput at constant mass — the 1/r source | a particle |

Sequencing per user (2026-07-28): the speedup ladder first (real
measured speedup), then M0/M1. Ensembles are mandatory (foam chaos
±30%), which is exactly why the speedup matters.

## 6. Log

* **2026-07-28** — document opened. Standing table laws_V2g (21/21).
  Known: blobs leak (−0.13…−0.54), self-insulate under pressure,
  footprint maintained, no steady throughput (G4). Nothing yet meets
  C1. Speedup ladder in progress (bench/BENCH.md).

* **2026-07-28 (later) — M0 + first M1: THE LEAK IS ROUGHNESS, AND
  CLOSURE RETARDS IT.** Instruments landed (conversion ledgers # CONV,
  lump census # LUMP, init=ring seeder; 21/21 inert gate, f2c8d5d).
  M0 taxonomy, heavy blob [50,300]: lump loss −0.0879/t.u. of which
  roughness (off-rung deliveries radiating D→F) = −0.0455 (52%),
  bulk evaporation = 0.0000 EXACTLY (a settled blob never exceeds cap
  — pass 6's over-full door never opens), condensation reclaims +0.005,
  spread-to-crumbs −0.008; residual ≈ −0.039 consistent with in-flight
  dense inventory growth (tET check owed). So the stable-particle
  problem IS the consonance problem: stop the internal deliveries being
  rough and the drain stops — exactly what a ring of rung-locked pairs
  is designed to do (the R window vanishes on-rung).
  M1, equal-mass (~2–4 units): gaussian control −0.232 %/t.u.; ring6
  −0.058 (4.0× slower); ring12 −0.051 (4.5×); ring6 with naive winding
  m_w=1 −0.119 — the distributed twist DETUNES the locked gates
  (winding must be absorbed into the per-link rungs to help; redesign
  queued). Ring6 seeds close at exactly m=3 (closure/2π = 2.9988).
  T=1000: ring6 leak DECELERATES monotonically (−0.00118 → −0.00081 →
  −0.00063 across [50,300]/[300,600]/[600,1000]) and consolidates to
  n=1; M 2.01 → 1.25. Under the C1 rate bar for the whole horizon but
  38% mass lost → C1 sharpened (see §2). ring12: flat ~−0.0008, n=2.
  Census caveat: thr=0.1 sees voice arcs as separate lumps early
  (spacing > rim threshold); n>1 does not by itself mean the chord
  broke — per-lump gate diagnostics wanted.
  **Queued next:** M1b ring6 at T=5000 (plateau vs ever-slowing decay),
  C3 ensemble over foam seeds, winding-compensated ring seeds,
  pressure-annealed ring (R1×R3), the flight-inventory residual check.
