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

## 5b. The tech tree (2026-07-28, after M0/M1)

Standing practice on every node: **analytics AND rendered frames**
(`viz/render_slice.py`) — visuals have caught what numbers missed.

```
A — CERTIFICATION (the first particle)
├── A1  M1b: ring6 T=5000 — plateau vs ever-slowing decay   [RUNNING]
├── A2  C3 ensemble: ring6 × 5 foam seeds, T=1000           [RUNNING]
├── A3  C2 self-force null + C4 perturbation identity        [after A1,A2]
├── A4  CERTIFY → the particle experiment ratchets into
│        the battery                                          [gate A1–A3]
└── A5  g5: space throughput at CONSTANT mass — the 1/r
         source (closes G4's open question)                   [after A4]

B — SEED ENGINEERING (better closures)
├── B1  winding-compensated ring (twist absorbed into the
│        per-link rungs; naive winding measured 2× worse)     [design]
├── B2  species scan: N × d × interval × m — the first
│        species table (lifetime + mass spectrum)             [after A4]
├── B3  chords of chords: two rings in mutual consonance
│        (the composite / binding question)                   [after A4]
└── B4  collider probe: hard packet strike → debris census
         (the confinement test: never a lone free voice)      [after A4]

C — NATURAL FORMATION (annealing)
├── C1' Es-reservoir piston apparatus (scheduled boundary
│        pressure, LEDGERED as an instrument account)         [design]
├── C2' pressure-annealed condensate census                   [after C1']
└── C3' pressure-annealed ring (R1 × R3 combined)             [after C1']

D — INSTRUMENTS
├── D1  viz/render_slice.py — the fresh eye (v89-native;
│        SFA/volview are pre-v89 and excluded)                [BUILT]
├── D2  flight-inventory ledger check (closes M0's −0.039
│        residual attribution)                                [open]
└── D3  per-lump lock diagnostics (gg inside census lumps —
         is the chord still singing?)                         [open]

E — THEORY DEBTS FEEDING MASS
├── E1  S2-full amplitude completion (wave-borne translation;
│        kr=1 full pass; radiation pressure; moving clocks)   [large]
└── E2  species enumeration re-run with the REAL resonance
         rule (construct_species × the comb law)              [after B2]
```

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

* **2026-07-28 (Q&A, recorded at user request)** — four questions,
  answered from the measurements:
  1. *Do a cell's harmonics circulate outside its own cell in a mass?*
     Yes, two ways measured: in-flight circulation (M0's −0.039/t.u.
     residual = mass living in transit between cells; flight-loads-
     pitch already treats channel energy as bound content of both
     ends) and phase coherence (the lock — gate closure, the ring's
     m=3 — is a property of no single cell). The dense AMPLITUDE
     itself is still cell-local at rate level; extending it is
     S2-full, and the S2-full failure set (radiation pressure,
     wave-borne translation) is exactly the cost of that locality.
  2. *Different types of particles?* Yes by design: species = distinct
     closed structures, discrete because closure forces integers —
     knobs are N, interval content (p:q), ladder index m (ring6
     closes at exactly m=3), winding (once compensated), composites.
     Already measured: different chords are dynamically different
     (unison ≫ fifth > octave lifetimes). Masses from tuning-curve
     loads ⇒ a discrete mass spectrum if closure survives.
  3. *Protons and neutrons eventually?* The honest ladder: certify one
     particle → species table (B2) → composites (B3) → then ask which
     structures play nucleon ROLES. Missing structurally: charge (no
     candidate law yet — winding/chirality imbalance are the
     leads) and spin (the two-plane chirality pair is the natural
     double-valued structure). "Neutron→proton decay" is a shape the
     model can express (two chords differing by one interval, close
     masses, one unstable to the other); whether it does is far.
  4. *Composite, or blobs the colliders see?* The dichotomy dissolves:
     a ring IS composite (the census resolves its voices; a hard probe
     would scatter off them as point-like centers) but the voices are
     NOT particles — alone, a voice-scale blob dies at −0.232 %/t.u.
     (the measured control) while the same voices persist 4× longer
     inside closure. Sub-particles are roles in a process, not
     separable constituents; knocking one out should yield only new
     closed structures + radiation, never a lone free voice — the
     confinement phenomenology, testable as B4.

* **2026-07-28 (parallel design round, two agents)** — both design
  notes landed and are archived in `design/`:
  **B1 winding compensation** (`design/B1_winding_compensation.md`):
  the decisive reframe — *winding IS the choice of per-link
  retardation*; the naive seed's phases were right and its LOADS were
  wrong. Compensated seed: no phase kick (`ring_wind=0`), loads set so
  ω·d/C = π − 2πw/N per link; the closure integer m = N/2 − w becomes
  the topology/species label, verified by the seeder's own closure
  printout. Concrete: comp12 (N=12, d=1.25, ring_x=0.32054 → m=5) with
  a mass-AND-pitch-matched unwound control (d=1.5, same ring_x → m=6)
  — a clean topological A/B differing only in the closure integer.
  Protection argument: unwinding requires a link to cross the dead-gate
  desert (gate ~4e-10 at ψ≈−2.6). Risks: N=6 w=1 is fully one-way
  (back gates 1.5e-5); gg diagnostics read low BY DESIGN on wound
  rings (score by leak + circulation, not gg). Prediction worth
  chasing: a chirality-split mass spectrum (m = N/2−w heavy vs
  N/2+w light) — the first charge-like species axis.
  **C1' piston** (`design/C1_piston_design.md`): cflag-6 skin +
  explicit reservoir accounts R_s/R_f INSIDE the conservation sum
  (Sterbenz-exact resets; drift stays at the floor); apparatus params
  piston_m/es0/es1/t0/t1/absorb; rung LOCI proven ambient-invariant by
  the bound-energy-only pitch law — so any measured locus shift under
  pressure FALSIFIES that law (a free tripwire); quasi-static bound
  t1−t0 ≥ 1000 at s_k=0.06. Two experiments specified (blob census
  arms A0/A+/A−; ring anneal B0/B+/B−) with success criteria.
  **Instrument caveat found by the agent, recorded:** existing
  edge_sink cells accumulate recorded Em which enters pressure via
  s_disp — old sinks slowly push space inward on long runs; prefer the
  piston for T > 500 once built (kernel-wide fix only via the ratchet).
  comp12 probe fleet launched (zero code change).

* **2026-07-28 (B1 fleet v2 — the accidental open-chain control, and
  the ring_m fix).** The fixed-`ring_x` compensated seeds REJECTED
  themselves by their own closure printout: foam jitter inflated the
  actual loop length ~10%, the lock recursion dumped the whole defect
  (~π) on the seam link, seam gate ≈ 0 — comp12/comp6/unwound ran as
  OPEN CHAINS, not rings. The topological A/B did not happen. But the
  chains are a control we didn't plan: chains leaked −0.092…−0.097
  %/t.u. where the well-closed auto-seeded rings (closure 2.9988) ran
  −0.051…−0.058 — **closing the loop is worth ~1.7× by itself**,
  independent of winding. Fix per the design note: `ring_m` seeding
  (uniform ω = 2π·m·C/L_ring from the ACTUAL picked loop) — measured
  closures now EXACT (5.0000 / 6.0000 / 2.0000). Fleet v3 (m=5 wound
  vs m=6 unwound vs m=2 one-way) running; seeder change ratchet-gated
  in the same batch. Standing rule: a ring seed is accepted only if
  its closure printout is integer to <0.05.

* **2026-07-28 (A1 + A2: NO PLATEAU — THE SKIRT IS THE PARTICLE-KILLER).**
  A1 (ring6, T=5000, open box): the deceleration through t=1000
  (−0.00118 → −0.00063) was not a plateau approach — past t≈1000 the
  leak ACCELERATED (−0.00146) and the ring dissolved entirely by
  t≈1900 (M=0, n=0). At death the per-voice load was ≈0.05–0.06 — the
  vacuum-skirt boundary x_skirt = 0.0617. Mechanism: leak thins the
  voices → pitch approaches the room's → the skirt dissolves them
  (e7's law, now the mass campaign's central enemy). Sharpened C1:
  FAILED for ring6 — lifetime ≈ 1.5–2e3 t.u. A2 (5 foam seeds,
  T=1000): the slow-leak class is seed-robust (late rates −5e-5 …
  −6e-4; 3–4/5 surviving at t=1000, one dead, one fragmented) — not a
  seed accident, but all on the same road. VERDICT: closure retards
  (4×; closing the loop alone worth 1.7×) but does not stop; a
  particle requires leak → 0 at some x* > x_skirt, i.e. an
  EQUILIBRIUM, not a slow death.
  **PRE-REGISTERED PREDICTION (before fleet v3 lands):** if skirt-
  death is the mechanism, v3_comp12 (m=5, per-voice x ≈ 0.32 — 5×
  above the skirt) must vastly outlive the m=3 ring6 (x ≈ 0.128);
  lifetime ordering by seeded load, not by winding per se.
  Follow-ups queued: A1c closed-box (no edge_sink: does the ring's own
  shed bath recondense — self-atmosphere vs starvation); the piston B+
  arm (ambient pressure raising the dissolution barrier) becomes the
  strongest stabilization candidate.

* **2026-07-28 (fleet v3, exact closures — compensation WORKS; skirt
  endgame pending).** Seeder gate 21/21 (ring_m inert to the battery).
  With closures exact (5.0000/6.0000/2.0000): wound m=5 vs unwound m=6
  at matched geometry class: −0.054 vs −0.047 %/t.u. — within foam
  chaos; the naive-winding 2× penalty is GONE. Winding is now carried
  at no lifetime cost (topology for free) but buys no protection at
  this horizon. The one-way m=2 ring shows the best late rate
  (−0.0034). Heavy rings (x/voice 0.21–0.44) retain 44–53% at t=1000,
  no early collapse — consistent with skirt-death but not yet decisive:
  at −0.05 %/t.u. their voices reach the skirt only at t ≈ 3000+. The
  pre-registered prediction is therefore SCORED by the T=5000
  extensions (launched: comp12_L, unwound12_L, comp6_L) — skirt
  mechanism predicts death ordering by per-voice load: m1-class light
  ring (died 1.9e3) ≪ unwound12 ≲ comp12 < comp6. Chirality
  circulation signature not yet instrumented (open: per-ring flux
  loop sum). A1c closed-box control running (self-atmosphere vs
  starvation).

* **2026-07-28 (endgame: the skirt prediction scored; first structural
  lifetime effect; longest-lived object yet).** Death times, all
  T=5000 horizon: light m=3 ring (x/voice 0.128): open box 1900,
  closed box 1600 — **the closed-box control rules out starvation**
  (a ring in its own shed bath dies the same; environment is not the
  killer). Heavy rings: unwound12 (x≈0.21) died 2221; one-way comp6
  (x≈0.44) died 3836; **wound comp12 (x≈0.32, m=5, mutual back gates)
  ALIVE at t=5000** (M=0.59, n=1, still decaying ~−3e-4). VERDICT on
  the pre-registration: the load axis is CONFIRMED at class level —
  lifetime rises steeply with per-voice load (1.6–1.9k → 2.2k → 3.8k →
  >5k) exactly as skirt-death predicts. Within the heavy class it is
  NOT load-monotone: the wound mutual N=12 ring outlived the heavier
  one-way N=6 — the first evidence that STRUCTURE (winding/mutuality/
  N) buys lifetime beyond load. Attribution caveat: single seeds,
  comp12 vs comp6 differ in N, w, and mutuality simultaneously —
  ensemble + a one-axis A/B (wound vs unwound at matched N=12 ran:
  2221 vs >5000, differing ONLY in closure integer m=6 vs m=5 at
  slightly different loads 0.21/0.32 — the cleanest pair we have, and
  the wound heavier ring won by >2.3×). comp12 is the longest-lived
  object the program has produced. STILL NO EQUILIBRIUM: everything
  leaks toward the skirt; the stabilization routes standing are the
  piston (raise the ambient barrier) and roughness elimination.
  Next: piston build (C1'), then the B+ pressure arm on comp12.

### R5. The pressurized bubble — surface tension + the higgs hole
(user proposition, 2026-07-28)

The proposition: a stable particle has strong surface tension (a
consonant, self-insulating skin) and a DENSITY HOLE — the "higgs
hole", the reverse of matter: dense-empty, space-rich — which is
load-bearing in the circulation/tension pattern. v89 translation: an
enclosed Es overpressure pushes OUT (pressure law) while ambient
pushes IN — a Laplace bubble with an equilibrium radius: the restoring
term every ring lacked (a 1D loop encloses no volume; only a closed
2-surface can trap a pocket). NOT the excluded pre-v89 shell import:
this is built from the v89 space-transport law, which did not exist
before.

**Derived selection rule (new):** a consonant closed surface must be
BIPARTITE — the reachable pair rung is the antiphase π-rung, which
cannot close around odd cycles; triangulated shells (icosahedra) are
frustrated; the CUBE (all 4-cycles, face diagonals beyond the contact
ceiling) is the minimal consonant shell. If shells become species,
"why these shapes" starts here.

H-ladder: H1 the bubble (cube shell + pressurized core) vs three
controls (shell without pressure / solid ball at matched mass / naked
pocket); H2 tension spectroscopy (breathing mode + fracture threshold
— "strong glass" = stiff and brittle vs flowing); H3 the circulation
round (toroidal tube enclosing the pocket — winding + hole + surface
combined) if H1 shows the hole matters.

* **2026-07-28 (H1 pre-registration, before the fleet lands).** Seeder
  facts: foam discreteness inflates cell-scale cubes (nearest-vertex
  picks: a_target 1.25 → abar 1.43, x/voice 0.264); seed gates
  imperfect (min ≈ 0 on a co-tree edge — cycle defects concentrate);
  we let entrainment finish the seeding (e6: lock acquisition ~10
  t.u.) and MEASURE whether the shell self-locks. PREDICTIONS, per the
  proposition: lifetime ordering bubble > shell-without-pressure >
  solid ball (matched mass); the bubble's core Es excess persists
  (trapped by the skin) where the naked pocket's disperses at the bare
  pass-S rate; equilibrium signature = late |dM/dt| ≤ 1e-4 with core
  pressure maintained. If the hole is NOT load-bearing (bubble ≈
  shell0), the proposition's hole clause fails honestly while the
  shell clause may still stand.

* **2026-07-28 (H1 v1: the hole was empty — instruments caught it).**
  `# shell core: cells=0` — the tiny cube's vertex picks consumed the
  central cells, so the "bubble" arm ran with NO pocket: h1_bubble and
  h1_shell0 produced byte-identical logs (death t=1222 — incidental
  determinism proof, no hole test). Valid scraps: unpressurized
  mini-shell ≈ matched solid ball (1222 vs 1158) — no enclosure
  benefit at this size and seed quality (gates mean 0.64, min 0).
  Seeder fixed: vertex picks exclude the core ball; core seeding
  guarantees ≥1 cell (2 seeded); edge-uniformity refinement passes.
  The exclusion pushes the cube outward (abar 1.59, x/voice 0.39 —
  heavier, farther from the skirt; both arms share the seed so the
  A/B stays clean). H1 v2 running: bubble vs shell0 vs pocket,
  T=5000, ratchet gate chained. Pre-registration unchanged.

* **2026-07-28 (H1 v2: hole test confounded AGAIN — correction).** My
  v2 log claim "both arms share the seed" was WRONG: the core
  exclusion was conditional on es_core>0, so shell0 picked a tighter
  cube (abar 1.237, x/voice 0.118 — skirt-class light) while the
  bubble sat at abar 1.586, x 0.387. shell0's death at 237 is light-
  shell skirt physics, not a pressure result. What v2 does show:
  bubble (heavy shell + 2-cell pocket) died 1749 — below comp12's
  >5000: the bubble class does not yet beat rings at this seed
  quality; the naked pocket disperses by ~t=100 (the bare pass-S
  rate the skin must beat). Fix: core exclusion now unconditional on
  r_core (geometry identical across arms; only Es differs). h1_shell0
  v3 running with matched geometry + gate. The A/B that decides the
  hole clause is pressure on/off at IDENTICAL picks.

* **2026-07-28 (H1 v3: THE CLEAN A/B — the hole clause falsified at
  this scale, and the reason is a mechanism statement).** Matched
  geometry at last (identical seed lines, abar 1.586, x/voice 0.387):
  pressurized bubble died 1749; unpressurized twin died 1814 — no
  effect within foam chaos (±30%). The 2-cell ×1.3 pocket buys
  nothing, and the arithmetic says why it could not: (a) the pocket's
  excess (~0.6) is ~10× smaller than the shell's leak by t=1000; (b)
  more fundamentally, THE LEAK LEAVES AS FIELD (roughness radiation,
  D→F) — space pressure pushes the space mode and cannot dam a field
  channel. The hole cannot fix this killer at any scale UNLESS the
  skin first closes the field channel. VERDICT on the proposition:
  hole clause falsified as tested; the surface-tension clause is
  SHARPENED into the requirement — a particle's skin must be
  CONSONANT (on-rung or gate-closed everywhere: nothing rough to
  radiate). No jittered mini-cube achieves it (gates mean 0.6). The
  standing stabilization routes after H1: (1) consonant-skin
  engineering (larger/annealed shells where jitter averages out;
  entrainment-grown rather than seeded skins), (2) the piston
  (outside pressure — C1'), (3) S2-full (if translation is the
  current, rough deliveries may reorganize). H2 tension spectroscopy
  remains meaningful; H3 (circulation torus) deprioritized until a
  skin can hold.
