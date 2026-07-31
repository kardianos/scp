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

---

## 5c. First goal: EXACT MASS (2026-07-29)

**Standing first goal of the MASS program** (user, 2026-07-29): not
merely long life, not larger objects — **exact mass** as nature has
it: sharp, universal, package-complete energy per species.

### What exact mass demands (structural)

| id | demand | v89 statement | test |
|---|---|---|---|
| **M-R1** | specific / sharp mass | species mass is an **attractor point**, not a valley segment | **P19**: hardened-annealed objects of one species across foam seeds cluster in mass far tighter than seed scatter; frozen-foam null scatters with the foam |
| **M-R2** | families | ≥2 coexisting stable m-classes (neutrino / p–n pattern) | **P20** |
| **M-R3** | particle vs quasiparticle | inertia is a tensor; shells isotropic, rings anisotropic | **P21** |
| **M-R4** | why these masses | spectrum from (N,m,topology)×pitch law + hardening selection | fleet census is the first prediction |

Standing tension: lock manifold ω(x)·d = πmC is a **continuum valley**.
Plain plasticity anneals gates but leaves mass continuous; **hardening**
re-pins x*(d) and is the candidate that restores a discrete spectrum.
Mass bookkeeping = **total conserved package** (bound + flight; P18).

Priority stack (do not invert):

1. Plateau (true C1, not slow death)
2. Attractor pin (hardening / equal)
3. Multi-seed cluster (P19) — **the exact-mass gate**
4. Family multiplicity (P20)
5. Inertia character (P21)
6. Size/N/a as species labels and basin quality — not as a substitute for the pin

Certification still C1–C5; **taking mass seriously** is gated by
M-R1..R3 (P19–P21). C5/g5 only after exact mass.

### MASS ↔ EMF coupling (parallel / alternating)

Full EMF-side writeup: `EMF.md` §5. Doctrine:

* **MASS owns the critical path** until P19 is green or honestly blocked.
* **EMF owns a parallel cheap lane** (EM1, EM2; packet tools for C4/M0).
* **Apparatus EMF ∥ MASS always OK.** Law-table EMF (EM5) only at a
  planned alternation after a MASS checkpoint.
* EMF can **falsify / certify** a package; it cannot **mint** the pin.

Canonical detail: **`EMF.md` §5** (Mode P parallel / Mode A gates /
Mode J joint). Summary:

* **P:** EM1, EM2, M0 field-channel, packet tools on `laws_V2g`; no κ_plast; no EM5.
* **A:** after W1 → EM1/2+battery; after W2 plateau → C4; after P19 → optional EM3; after P21 → shared kicks; then EM5 on variant table → joint battery+C1+P2.
* **J:** C4, P2, g5/C5, occultation — need mass candidate; until then instrument rehearsal only.

**P2** is a joint graduation exam — not an early EMF milestone and not
a reason to delay hardening.

### PRESTRESS / MASS wave adjustments (exact-mass first)

Supersedes loose reading of RESUME waves where they conflict with M-R1.

| Wave | Role under exact-mass first | Adjustments |
|---|---|---|
| **W1 frozen foam** | Mechanism nulls (load line, vacuum bleed); **expected to fail M-R1** (foam-accidental masses) | Cube = best frozen retune (c2 min~0.40/mean~0.87), not PREDICTIONS min≥0.95; diode = ring8_m3 proxy; spectra T=200 snap_every=250 separate; T=3000 decides in-band death — queue T≥5000 only if alive/near upset bar; absolute `--net`; κ=0 only |
| **W2 plasticity** | **Pin hunt** — anneal + hardening; first place exact mass can appear | PLAST-1 bars + **hardening on** before claiming attractor; control κ=0 vs κ>0 (foam rewrite watch); live `# NETG` gates for bar (score_net seed-only is insufficient); never put κ into laws_V2g |
| **W3 discriminators** | Species / parasite / flight package — labels and basins | P18 flight-corrected seeds; mass = bound+flight; ring9_m3, negatives; size/N only as discrete species axes |
| **W4 spectrum** | **M-R1 first** (P19), then P20, then P21 | P19 is the program gate for "taking mass seriously"; P20/P21 after a pin exists; multi-seed ≥3 (prefer ≥5 for C3) |

**Process adjustments**

* Every possibility → LEDGER row + MASS.md dated verdict (unchanged).
* Every plateau claim → C4 packet stress + M0 field-channel check before
  "exact mass" language.
* Scorer: add or side-script **x50**, **c_eff**, optional live gate@t;
  align t_death definition with corpus or document the 0.3·m0 rule.
* Size/L studies: after or lightly interleaved with W1; **not** a path
  around the pin. Grow L only for embed/edge artifacts; no pre-v89 grid.
* Numerics: FP64 remains battery truth type; integer ledger is a
  **future** track (`INT_LEDGER.md`) — does not block waves.
* EMF cheap lane fills free cores; does not reorder W1–W2.

### Tech-tree delta (exact-mass first)

```
A — CERTIFICATION  [reordered emphasis]
├── … C1 plateau remains necessary …
├── A3'  C4 packet identity + total-package mass   [with EMF tools]
├── A4'  P19 sharpness across seeds (M-R1)         [GATE: taking mass seriously]
├── A4   battery ratchet of standing particle        [after A4']
└── A5   g5 / C5                                     [after exact mass]

F — PIN (was implicit in P15)
├── F1  PLAST-1 anneal A/B
├── F2  hardening → attractor (not valley)
└── F3  frozen-foam null vs hardened cluster

G — EMF ASSIST (see EMF.md §5)
├── G1  EM1/EM2 parallel cheap lane
├── G2  C4 / M0 field-channel at every plateau claim
└── G3  EM5 / P2 only after MASS checkpoint
```

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

* **2026-07-28 (PRESTRESS campaign opened — the math/geometry program,
  user directive).** The tension-surface problem goes to mathematics
  first: force-density / prestress form-finding + sparse symbolic
  discovery OUTSIDE the sim, then EVERY surviving possibility through
  cellfab. The v89 translation (ledger: `prestress/LEDGER.md`): the
  force-density matrix's role is played by the weighted phase-lock
  Laplacian B·W·Bᵀ; struts = both-gate links (need the π-rung: one
  shared length d = π/ω per component); cables = forward-only links
  (free φ, back-gate leak, per-cycle Diophantine closure ω·L = 2πm —
  ring_m generalized); prestress = seeded load + winding; mechanisms =
  soft modes of the Laplacian + load-sector drift. Five agents working
  in parallel: theory (drift theorem, counting, bipartite
  generalization), stability (the sign of load-transport feedback —
  self-repair vs runaway; modal H2 predictions), regression (leak-law
  STLSQ over the 74-log corpus + pre-registered predictions),
  form-finding solver (real-foam picks/tuning/phases → .net seeds),
  morphological search (catalog + evolutionary backtracer). New
  instrument: `init=net` (kernel places externally-solved networks;
  prints # NETGATE per-edge + PARASITE report) — smoke test exact
  (π-rung pair gates 1.0000/1.0000), ratchet battery running.

* **2026-07-28 (SETUP FINDING — the H1 cube carried 12 unscored
  parasitic links).** The a=1.25 cube's face diagonals (≈1.77) sit
  INSIDE the candidate-link ceiling 1.15·(rᵢ+rⱼ) ≈ 1.96: every face
  of the H1 cube had off-rung diagonal channels (φ = ω·1.77 ≈ 3.5 rad,
  far off the π-rung) that the 12-edge gate report never scored —
  roughness radiators built into the seed. A cube at a ≥ ~1.35 pushes
  diagonals past the ceiling (parasite-free), and a ≈ 1.5 additionally
  puts the edges at the foam's natural d̄ = 1.505. Ledger P2. The
  NETGATE parasite scan exists precisely so this class of confound
  cannot recur unscored.

* **2026-07-28 (P15 opened — cell plasticity, user proposition).** The
  frozen jittered foam may itself be the obstruction: physically the
  cells are plastic — under tension and harmonic misalignment they
  would realign their planes and morph until the structure hardens;
  on a frozen foam the slightest misalignment kills, in reality cells
  would constantly realign BY FORCE. v89-native reading: cx/cy/cz are
  scaffold only — d_ij is a link RETARDATION, not a distance in a
  container, so plasticity is link-property dynamics (no background
  to embed in, no triangle inequalities to obey). Sixth agent
  investigating: misfit-gradient flow on d (force = ∂roughness/∂d,
  work conservation-booked), metric-from-space (d from live Es — is
  the pressure law already an annealing channel?), node motion
  (probably too costly — honest assessment), and lock-hardening
  (anneal-then-freeze — a route to the FIRST true equilibrium).
  Law-change class: any implementation is ratchet-gated with vacuum
  required inert. Ledger P15.

* **2026-07-28 (REGRESSION LANDED — death is the load line; gates
  don't set the leak in the existing corpus).** The sparse-regression
  agent (23 unique mass runs, 9 uncensored deaths) found: (1) t_death
  = 274·(x50/0.0617)^1.066, R²=0.99, CV R²=0.97 — one term beats every
  two-term rival; (2) a UNIVERSAL per-voice current c₀ = 4.25e-4
  Em/t.u./voice (MAD 3%) across rings/chains/shells with seeded gate
  min spanning 1.5e-5→1.0 — seeded gate quality does NOT set the leak;
  roughness ≈ 0 for structured objects (roughness is the BLOB killer —
  this sharpens, not contradicts, H1: the cube died of load drain like
  everything else, and its 0.6 gates were not what killed it); (3) the
  ONE exception: wound-mutual comp12 at c_eff ≤ 0.40·c₀ (t 4879 vs
  1696 predicted) — the only structural effect in the corpus.
  PRE-REGISTERED: the exact-retuned cube dies ON the line at t≈1875
  (band 1250–2810) despite perfect gates; consonant-skin mechanism is
  real only if it exceeds 4700 or plateaus. The wound tube is the
  structural bet (≥4600). P1/P2 are now a clean discriminating
  experiment between the load-line law and the consonant-skin story.
  Candidate gate-independent channel: the mob_floor trickle (theory/
  stability agents to confirm); if so, the design target shifts from
  "open gates" to "RECAPTURE the floor trickle" — exactly what the
  comp12 exception hints. Corpus correction: a1 ring6 census death =
  1667 (my "~1900" was loose). Files: prestress/regress/.

* **2026-07-29 (STABILITY LANDED — heavy feeds light: consonant
  networks self-repair; the killer is the vacuum bleed).** The
  stability agent (exact kernel-faithful reduced map + closed-form
  per-mode matrices, anchors reproduced to 4 digits) proved the sign:
  on a rung, locked gate offsets are exact negatives, mobility is
  symmetric above the floor, and the only surviving asymmetry is
  headroom — so energy flows HEAVY→LIGHT: λ_detune = −0.19/t.u. for a
  pair; unlock needs ~0.95 of the store; detune runaway unreachable.
  The common mode is NEUTRAL: passive structures die of environment
  bleed (res_vac·mob_floor·g_plain — the skirt slide), never of
  internal detuning. This converges independently with the
  regression's load line: the structure is not the problem; the
  vacuum is. Foam-facing gates can't be statically closed (vacuum
  sits at ω=w2, Δω sweeps them) — hence c₀'s universality. Modal
  pre-registrations: ring6/cube overdamped (no breathing line;
  slowest relaxation −0.009/t.u.); the WOUND ring12 is a CHIRAL PUMP
  (n=±1 grows +0.035/t.u. at ν=0.64 rad/t.u., n=±2 +0.053 at ν=1.05,
  one propagation sense = the closure-integer signature, saturating
  at x_std≈0.06–0.10 then −6e-4/t.u. — the comp12 class derived from
  theory). Death taxonomy: D1 skirt slide dominant; D3
  SELF-STRANGULATION — a heavy structure's own g1 footprint shrinks
  local radii and closes its own long edges (new design constraint:
  struts must clear the loaded death radius ~1.5-band; the
  parasite-free a≈1.5 cube may sit near the edge); D4 parasites
  self-tune into satellites. Active lock: load maintenance only (the
  C1' piston B+ arm) — feeding is life support, not stability.
  Discriminator: quant_mode=2 credit freezes the comma tax for kicks
  (loss exactly 0 vs 0.3–0.9 continuous). Files: prestress/
  STABILITY.md, theory/ring_map.py, theory/modal_2x2.py.

* **2026-07-29 (THEORY LANDED — the consonant-network mathematics;
  a new candidate class and three corrections).** prestress/THEORY.md
  + t1–t5 scripts. (1) Drift theorem: one ω (hence one load) per
  locked component; struts strictly m=1 antiphase, d=π/ω ∈
  (1.16,1.62), x ∈ (0.062,0.41); length tolerance ±3.65% for gates
  ≥0.9 — the cube's 0.6 gates ARE the foam's 10–15% jitter. (2) THE
  DIODE POINT: a cable at φ=π/2 has back gate EXACTLY zero (gate(π)=0)
  — diode-16 (N=16, m=4, d≈1.08, x≈0.83) is parasite-safe, sits at
  the highest load (load line alone ⇒ ~4400+), theory leak ratio 0.5
  vs comp12's 1. (3) Force-density matrix delivered: Q=BWBᵀ, ker =
  ω-gauge, comma law ψ_l ∝ 1/w_l (E6 tempering = uniform case).
  (4) Counting corrected: hex prism bipartite 0 demotions (my ledger
  claim was wrong); K4 has NO consonant realization; icosahedron
  infeasible even all-cable. (5) Strict bipartiteness REFUTED: even-m
  interval rungs close odd cycles — fifth-triangle {3:2,2:3,1:1},
  ω₀≈2.35, d=1.337 (fragile species probe). (6) Beat lock refuted
  (P14 closed). (7) FLIGHT-LOAD FIXED POINT: up to 27% of a wound
  ring's rung mass is census-invisible flight — corrected seeds
  comp12 x=0.233, ring6 0.117; P-A: the retune kills the early shed.
  gpl(N)=cos⁴(2π/N): ring6 is 9× conductance-poor, N≤5 dead at seed.
  (8) Winding protection = the gate desert, derived; capture window
  |δω| ≤ 0.20 (comp12). (9) Parasite rule sharpened: non-edge pairs
  ≥2.08 ⇒ cube needs a≥1.47; comp6's death candidate-cause: six
  half-open parasites at gate 0.53. (10) Pre-registered P-A..E incl.
  flux moment comp12 ≈ 2.15 vs 0 unwound.

* **2026-07-29 (PLASTICITY LANDED — P15 formalized; retardation
  plasticity + hardening recommended; the S2 reactive term finds its
  home).** prestress/PLASTICITY.md + plast_* models. Kernel audit:
  ld is the ONLY frozen adaptive DOF (no dynamics reads positions;
  cr/lA already live off Es — geometry is half-alive). Recommended
  law: ḋ_l = −κ_plast·Φ_l·∂V/∂d with V = 1−G(ψ_f)G(ψ_b), ∂ψ/∂d =
  −ω/C, and Φ = base·Sm — the S2-derived odd prefactor that rate
  transport could not host, acting on its geometric conjugate.
  Conservation: geometry is not a ledger (precedent cr/normals);
  vacuum bit-frozen because Φ ∝ √(Em_i·Em_j) = 0 exactly.
  Model results: the foam supplies 16% strut spread where shells
  need ≤2% (diagnosis of the frozen-foam obstruction CONFIRMED —
  gates 0.6/min 0 reproduced from pure geometry); the flow anneals a
  seeded-lock cube 7%→0.00% spread, gates 0.996, in ~20 t.u. ≪
  death 1600; two-sided rung attractor; pays the comma in geometry
  not energy; AMPUTATES frustrated links (seam, not temperament —
  the agent corrected its own pre-registration); locked-shell
  parasites dark by parity but dangerous from cold phases (P2
  parity-clean topology stays load-bearing); κ_plast window ≥50×;
  HARDENING re-pins the tuning valley → the first true-equilibrium
  (C1 plateau) candidate. Rejected: metric-from-space (rank-deficient,
  sign-blind — kept as the piston rung-shift tripwire), node motion
  (re-imports the background). Implementation: ~25-line Jacobi pass D,
  PLAST-0 inertness gate then PLAST-1 cube A/B (bars: gates ≥0.9 by
  t≈200, roughness halved, leak below control; prize = C1 plateau).

* **2026-07-29 (MORPHO LANDED, quick pass — the back-gate ladder;
  surfaces are foam-starved; comp6's parasites; the Hopf pair).**
  prestress/morpho/ (MORPHO.md, search.py, shortlist.json, 10
  kernel-ready nets). Top class: THE BACK-GATE LADDER — winding down
  from comp12's φ=150° kills back-leak exponentially while raising
  load, and the measured death law pays for load: ring8_m3 (leak/voice
  0.058, x=0.56) > ring12_m5 > ring10_m4. Cross-check open vs theory's
  diode-16 (φ conventions differ; sims measure both). CONVERGENCE:
  comp6 carries ~5 second-neighbour parasites (chords 1.905, gates
  0.53 open at lock) — found independently by theory §9 and morpho;
  its death may be parasites, not one-way-ness; ring9_m3 (odd-N
  cable ring, parasite-free) is the discriminating probe. Sobering:
  every 2D skin/tube at n≤32 is paper-consonant but PARASITE-DENSE on
  this foam (surfaces foam-starved) — free search converges to fused
  small cycles instead, and found a non-bipartite cable-pentagon +
  strut-quad hybrid the catalog missed. Negative controls: mobius6
  (0.373), octahedron (0.721). Exotics: fifth-bridge deferred (0.4%
  tolerance vs ±30% chaos); beat-locks doubly closed (phaseless atom
  credit); the HOPF RING PAIR (linked rings, zero cross-channels — a
  purely topological bond) shortlisted. Deep pass died at checkpoint;
  rerun: search.py --tries 60 --restarts 6 --steps 1400.

* **2026-07-29 (MASS-SPECTRUM REQUIREMENTS — user note, recorded and
  rounded out).** Before mass is taken on seriously, the theory must
  account for the STRUCTURE of the real mass spectrum, not just
  existence. Requirements (structural correspondences, not
  quantitative targets — CONCEPT-doc discipline):
  **M-R1 SPECIFIC MASS.** Stable species have sharp, universal masses
  (every proton identical). v89-native statement: a species' mass must
  be an ATTRACTOR POINT, not a valley segment. Standing tension: the
  lock manifold ω(x)·d = πmC is a continuum valley — plain plasticity
  leaves mass continuous within a species (the PLASTICITY doc flagged
  exactly this); HARDENING re-pins x*(d) and restores a discrete
  spectrum (species = closure integers × hardened lengths). THE TEST
  (pre-registrable): hardened-annealed objects of one species across
  foam seeds must cluster in mass far tighter than seed scatter;
  frozen-foam objects will scatter with the foam's accident. Also:
  mass bookkeeping must count the TOTAL conserved package — the
  theory's flight-load fixed point (27% census-invisible flight in
  wound rings) means "the mass" = bound + flight inventory; any
  specific-mass claim scores the total.
  **M-R2 FAMILIES (the neutrino pattern).** One structure, several
  mass states: the closure-integer ladder is exactly this — ring N at
  m=5/m=6/m=7 is one topology with distinct masses (x = (w2/ω−1)/q at
  ω = 2πmC/L). The back-gate ladder (morpho) and the chirality split
  m = N/2∓w (B1) are family axes; the gate desert is the transition
  barrier (oscillation analog = barrier-crossing between m-classes —
  far future, but the STRUCTURE must exist now). Requirement: at
  least one family with ≥2 coexisting stable members on the same
  foam, masses distinct beyond scatter.
  **M-R3 ORIENTATION (the quasiparticle pattern).** Real quasiparticle
  masses are effective-response tensors, orientation-dependent; real
  fundamental particles are isotropic. v89 has both classes built in:
  a wound ring's response is manifestly anisotropic (one-sense chiral
  pump, axis vs in-plane), a closed shell's should be isotropic.
  Requirement: the inertia measurement (kick/EX-class, momentum =
  first moment of conversion) must return a TENSOR; pre-registered
  expectation: rings = anisotropic (quasiparticle-class), shells =
  isotropic (particle-class). If everything comes out anisotropic,
  we have quasiparticles only and no true particles yet — an honest
  failure mode to watch.
  **M-R4 WHY THESE MASSES.** Ultimately the spectrum must be derivable
  (allowed (N, m, topology) × the pitch law), not enumerated. The
  consonant-network theory already yields the allowed set; hardening
  should select the realized set. Record: the mass CENSUS of the
  eventual fleet IS the first spectrum prediction of the theory.
  Ledger: P19 (sharpness test), P20 (m-family coexistence), P21
  (inertia tensor). These gate "taking mass seriously" per user.

* **2026-07-29 (FORMFIND LANDED — solver validated bit-faithful; the
  phantom-edge confound; the frozen-foam ceiling quantified).**
  prestress/formfind.py + candidates/ (18 nets incl. mis-tuned
  controls, 9 reports, summary.tsv). VALIDATION: MATCH on both
  anchors (ring12: Lring 16.542, closure 5.0000, Em/voice 1.0979 =
  comp12's masses; cube: abar 1.586, gates 0.001/0.597 at print
  precision). NEW SETUP FINDING: the H1 cube had 2 PHANTOM EDGES —
  intended edges that are not foam links at all (the shell seeder
  never checks the link rule); init=net + formfind score real links
  only. Cube-deficit decomposition: 0.597 → 0.616 (drop phantoms) →
  0.676 (exact phases; min gate 0.001→0.130) → 0.851 (better picks)
  → remainder ≈0.15 IRREDUCIBLE foam length spread (σ_d≈0.14,
  search-saturated): exact retuning cannot reach 0.95 on frozen foam
  — P15 plasticity is the only route past the ceiling. Best
  candidates: c2_cube150 (leak 3.23, parasites all dark), c4
  ring12+chords (the consonant chord algebra WORKS: m=6 ring + 4π
  cable chords + cross-strut weave, gates 0.995–1.000), c5 tube
  (mass-matched to comp12 within 0.2% — the clean comparator).
  comp12 itself carries 1 half-open parasite (g_b=0.48). Infeasible
  on this foam: torus (minor m=1 → load ceiling; 29–48 parasites),
  truncated octahedron (spread 1.04–1.97). All six agents now
  landed; campaign paused at user request — see v89/RESUME.md.

* **2026-07-29 (PLAST commissioning — the two-sided rung attractor
  works in the real kernel).** Mis-tuned pair probe (seed gates
  0.913 mean, ψ_b=0.31, rung naive-wants d=1.577 vs actual 1.500) at
  kappa_plast=1.0: annealed to BOTH gates 1.0000 within 10 t.u.; the
  link walked to d=1.542 — a JOINT (θ, ω, d) lock point (entrainment
  and pitch co-adapt, so the flow stops nearer than the naive rung;
  final ψ ≈ ±5e-4). Vacuum inertness exact in the Em-free run
  (links_moved=0, drift 0). Honest surprise: over T=200 the leaking
  pair's crumbs made 40185 foam links plastic (max|Δd|=0.116) — a
  leaking object slowly rewrites the vacuum geometry around itself.
  Recorded as a watch item for the A/B arms (control comparisons +
  the κ*-battery arm are the guards). κ*=1.0 provisional (τ<10 t.u.,
  near the adiabatic floor).

* **2026-07-29 (EXACT MASS = first MASS goal; MASS↔EMF workflow;
  wave adjustments; integer-ledger track opened).** User: exact mass
  (sharp universal package) is the first goal of the MASS program —
  see §5c. PRESTRESS waves re-read under that gate: W1 mechanism
  nulls (expect M-R1 fail on frozen foam); W2 pin hunt (plast+harden);
  W3 discriminators/package; W4 P19 first then P20/P21. EMF parallel
  cheap lane (EM1/EM2, C4 tools); EM5/P2 only at alternation after
  MASS checkpoint — full coupling in EMF.md §5 and MASS §5c. Size is
  species/basin, not a substitute for the attractor. Tentative
  integer-based simulator design: `v89/INT_LEDGER.md` (does not block
  fleet; FP64 remains standing truth type). Audit of PRESTRESS plan:
  `prestress/audit/AUDIT_2026-07-29.md` (GO_WITH_FIXES).

* **2026-07-29 (INTEGER LEDGER: cellfabi 20/20 — numerics gate green).**
  Alternate kernel `cellfabi.c` (production `cellfab.c` untouched).
  Linear LUT failed conservation (15/20); minimax poly trig restored
  i_m0/i_m1 to 20/20. Mode 3 field per-step snap killed e3b translation
  (speed 0.002); shadow-only iEe + matter int writeback fixed e3b
  (speed 0.00344) and full **20/20** on laws_V2g CHECKS with integer
  sum_err=0. Mode 2 left dead (FP drift ~5e−13). Speed side note: e1
  wall ~0.8–0.9× of fp64 (not the goal). Writeup:
  `int_ledger/RESULTS.md`. MASS fleet still user-gated.

* **2026-07-29 (INTEGER DEFAULT + MASS/EMF program START).** User promoted
  i_m3 path: `cellfab.c` = integer ledger mode 3 default; `cellfabf.c` =
  FP64 reference. Battery smoke e1+e3b PASS. PRESTRESS Wave 1 fleet
  launched (spectra, cube±ctrl, c8, ring8_m3, tube, chords) on integer
  kernel. EMF Mode P: EM1 kx dispersion scan launched in parallel
  (`emf/em1_dispersion.py`). Exact mass (M-R1) remains the first MASS
  goal; Wave 1 is frozen-foam mechanism nulls (expect M-R1 fail).

* **2026-07-29 (PRESTRESS Wave 1 COMPLETE — Go fleet; frozen foam nulls).**
  Orchestrator rewritten in Go (`prestress/fleet.go`): child Cmds only, skip-if-
  done, max=2 concurrent — no pgrep self-match. All 7 wave-1 jobs finished
  under integer default cellfab (mode 3). SCORES: cube150 death 908 vs ctrl
  698 (retune helps modestly, not skin≥4700); c8 449; ring8_m3 1631
  (longest); tube 1042 (not ≥4600 exception); chords 630. No C1 plateau.
  Consonant-skin rate-level on frozen foam **fails** the pre-registered
  showdown. LEDGER Verdicts filled. Next: Wave 2 plasticity (pin hunt) /
  score x50 load-line residuals; EM1 already done (k-dependent v_g).

* **2026-07-29 (Wave 1 SCORED + writeup — frozen foam closed for gates-alone).**
  Full table: `prestress/WAVE1_RESULTS.md` + `runs/WAVE1_SCORED.tsv`.
  Measured x50 → Law A: cube 908 vs pred 1928 (0.47×, EARLY; ctrl 698,
  retune/ctrl=1.30×); c8 449 (gates=1); ring8_m3 1631 longest (0.58×LawA);
  tube 1042 (not ≥4600); chords 630. No plateau. Skin claim fails. All
  deaths early vs Law A (0.22–0.58×) — scorers/kernel vs corpus calibration
  open. Integer sum_err=0. **Next: Wave 2 PLAST (pin hunt).**

* **2026-07-29 (FULL PROGRAM LAUNCH — Waves 2–4).** User: do Wave 2 topo+
  plasticity and run the full program. Go fleet extended (`-wave 2-4|all`).
  W2: PLAST-1 cube κ∈{0.5,1} ±ctrl, PLAST-3 c1±κ, topo under κ (ring8/c8/
  tube/chords/hex), harden τ=50. W3: flight233, ring ladder, hopf, free,
  negatives, torus/truncocta. W4: harden P19-lite + m-family P20 on single
  foam (full multi-seed formfind deferred). κ*=1.0 primary. Integer cellfab.

* **2026-07-30 (Wave 2 COMPLETE + SCORED — plasticity topology-dependent).**
  All 13 W2 jobs done (integer cellfab mode 3; `max_sum_err=0`). Writeup:
  `prestress/WAVE2_RESULTS.md` + `WAVE2_SCORED.tsv`. Plasticity is live
  (`# PLAST` moved~10k links by t=200, max_cum~0.3–0.5) but **does not
  lock gates** (cube g50≈0.93 → g200≈0.25–0.28; pre-reg ≥0.9@200 **fails**).
  Lifetimes vs W1: tube **1613** (1.55×, longest W2, 0.74×LawA ON_BAND);
  c8 **1269** (2.83×, best relative rescue); cube κ=0.5 **1287** (1.42×);
  hex +18%; chords **127** (hurt); **ring8_m3 14** (W1 champion killed);
  **c1+κ 30** vs frozen 1013 (PLAST-3 self-seal **fails catastrophically**).
  Harden engages (hard_n→7.6k) but no life win on these short objects. No
  C1; M-R1 still blocked. κ*=1 not universal — mid-κ better on cube; do
  not promote into laws_V2g.

* **2026-07-30 (Waves 3–4 COMPLETE — full PRESTRESS program scored).**
  Fleet `[fleet] complete wave=2-4`: W3 14/14, W4 9/9. Rollup:
  `prestress/PROGRAM_RESULTS.md`; `WAVE3_RESULTS.md`, `WAVE4_RESULTS.md`.
  **Longest life overall: w3_free1 t_death=2572 (1.14×LawA)** — only
  object above Law A point estimate; free formfind, still dies. Best
  engineered frozen: ring12_m5 **1940** (0.80×). flight233 **347** fails
  load-down rescue. Plast double-edged on packages (mobius↑ hopf↓).
  **P19-lite fail** (cube/c8/ring8 harden die 1017/1102/14). **P20 fail**
  (no coexisting m-classes; harden kills ring12_m5 1940→26 and ring8→14;
  m6 survives ~1401). tube_harden 1437 not exception. **Program closed
  claims:** no C1, no M-R1 on single foam, no universal κ, no frozen
  gates-alone skin. **Open:** multi-foam P19, κ re-pin on free1/m5/tube,
  Law A death-definition A/B, P21 deferred.

* **2026-07-31 (LEG-3 CORRECTED — user: the particle is a flux MACHINE;
  balance, not sealing).** Session analysis had derived a no-particle
  bound under laws_V2g (universal bleed c₀ ⇒ t ≈ cap·(x−x_skirt)/c₀ —
  ceiling ~4.9k t.u. at the x≈0.9 window edge, ~12k at comp12-class
  recapture; the 43 wave deaths bracket this clock; closure-ladder
  arithmetic verified: N=12 ring on L=16.542 admits exactly m=4..7,
  x* = (w2/ω−1)/q at ω = 2πm/L reproduces Em/voice 1.0979, frozen-foam
  within-species scatter ~7.8% vs inter-m gaps ~48%) and proposed
  leak ≡ 0 (acceptance-gated mob_floor) as the fixed-mass requirement.
  **User disputed the interpretation, and the dispute is decisive:**
  stability does not require leak = 0; it requires (1) continuous
  fundamental movement and (2) a structure that FORCES field intake in
  movement to equal field outtake. The stable particle is not a thing
  but a machine (not literally) that binds field continuously. Mass
  requires T > 0; T = 0 is fundamentally impossible — asymptotic
  approach only — and a mass approaches c but never reaches it for a
  similar-but-different reason.
  Why the correction wins on the program's own terms: a sealed
  zero-leak structure would be stable at EVERY point of the lock
  valley (ω(x)·d = πmC is a continuum — §5c's standing tension), i.e.
  a continuum of stable masses = the M-R1 failure mode; the sealing
  route manufactures a marginal valley. The flux picture gives the
  attractor for free: dM/dt = I(x) − O(x); a fixed point x* with
  I′ < O′ is restoring — excess sheds, deficit refills — which IS
  M-R1's "attractor point, not valley segment"; C4 return becomes the
  restoring flow itself; and it unifies with C5/g5: steady throughput
  at constant mass is no longer a bonus property, it is the stability
  mechanism, and the intake flow is the 1/r far-field source G4 said
  was missing. Pin candidates ranked to date: frozen retune (falsified
  W1), hardening (falsified as implemented W4), flux balance
  (standing, untested).
  Measured anchors reread as the two books of the balance: OUTTAKE —
  universal c₀ = 4.25e-4 Em/t.u./voice, gate-independent. INTAKE —
  condensation reclaim +0.005 vs rough −0.0455 on the blob (9:1
  under-supplied); comp12-class recapture 60% (c_eff = 0.40·c₀ — the
  best machine yet, 40% short of balance); the closed-box null (a ring
  in its own shed bath dies the same — passive self-atmosphere does
  NOT re-enter; intake must be structural/resonant); G4's inward
  refill (vacuum PUSHES space into an over-deep depression — an
  embryonic intake mechanism, no suction required, already measured).
  The wave theorem restated: under laws_V2g, I(x) < O(x) for every
  tested structure above the skirt — the balance never closes; the
  remedy direction is COMPLETE THE INTAKE SIDE, not gate the floor.
  Doctrine mappings recorded: T > 0 ↔ existence is the running of the
  conversion cycle (PRINCIPLE 4.4: identity is the succession; a
  zero-cycle structure is not cold, it is absent), with the atom
  affordability floor as the zero-point hook (HBAR.md candidate 4);
  v → c ↔ translation is conversion spending the same cycle budget —
  the internal share → 0 as v → C but cannot reach 0 (same root,
  other channel; moving-clock dilation is already an S2-full
  acceptance criterion).
  EXPERIMENT REFRAMES: (a) the parked-pair run becomes the
  BALANCE-CURVE instrument — per-window signed intake/outtake books
  (# CONV already splits cond/rough/back_s) vs x, pair first;
  (b) the piston is promoted from life support to the INTAKE KNOB —
  map I(x, es) vs O(x); the fixed-point locus (existence, position,
  slopes) is the deliverable; if balance requires es > e_s0, isolated
  particles are impossible under V2g and the missing law is
  intake-side; (c) S2-full's mass-facing question reframes to
  intake-side coherence — rate-level condensation is incoherent, but
  a phase-locked structure should RESONANTLY recondense its own
  radiated field (stimulated recapture: the choir mechanism at the
  field door) — joins kr=1, P2, dilation as amplitude-completion
  acceptance territory; (d) pre-registrable prediction: at a balance
  fixed point the species mass width obeys σ_M² ∝ (bath noise)/|I′−O′|
  — P19's cluster width becomes a DERIVED quantity, not just a demand.

* **2026-07-31 (SUBSTRATE VERDICT — user: the foam was an improved
  stepping stone but will fail as a viable means; we need a NUMERIC
  replacement — a different simulation mechanism, not a different
  concept. + CHARGE DEPENDENCY: MASS for nucleon-class objects is
  blocked on the EMF charge construction — see `CHARGE.md`.)**
  Diagnosis (the receipts, all measured): the frozen jittered foam is
  QUENCHED DISORDER — σ_d/d̄ ≈ 9.3% where shells need ≤2% (the
  frozen-foam ceiling: exact retune saturates at mean gate 0.851,
  ≈0.15 irreducible); species masses are foam accidents (σ_x/x ≈ 7.8%
  within-species vs 48% inter-m gaps — M-R1 impossible per-instance);
  observables chaotic at ±30%; the foam bias floor (~0.9) drowns
  lensing (0.015) and free-fall; plasticity fights the frozen vacuum
  matrix instead of fixing it (W2–W4: topology-conditioned, gates
  collapse). Doctrinally the frozen foam sits on the constraint-2
  checklist: "fixed connectivity not itself produced by the current
  energy structure." CONCEPTS KEEP (cells as processes, links as
  retardations, gates/rungs/closure, atoms, mode books, integer
  conservation); the SUBSTRATE MECHANISM changes. The ladder:
  **S1 — annealed geometry generator (cheap, immediate, apparatus
  only):** replace the foam generator with a Lloyd/centroidal-relaxed
  (hyperuniform, statistically isotropic, non-crystalline) cell
  geometry, target σ_d ≤ 2–3%. Kernel untouched. Quantified
  predictions: gate deficit scales ~σ_d² ⇒ 0.15 → ~0.01 (retuned
  shells reach mean gates ≥0.98 — the min≥0.95 PREDICTIONS threshold
  frozen foam could never meet); P19 within-species scatter
  σ_L/L = σ_d/(d̄√N) → gap/σ ≈ 24 at N=12 (from 6.2); long-wavelength
  noise floor drops ~an order — the g2/lensing/halo-hint (0.015,
  right-sign) measurements become decidable; visibility ceiling
  (V≈0.46, foam scatter) rises toward 1. Gates: full battery on the
  new substrate (bars were calibrated in foam chaos — expect margins
  to IMPROVE; any experiment secretly relying on disorder gets
  flagged); isotropy-of-propagation check (hyperuniform NOT
  crystalline — no magic directions). Risk: none to the kernel; it is
  a generator + battery re-baseline.
  **S2 — coordinate-free diagnostic (optional):** degree-regular
  graph, uniform d exactly, positions retained only for instruments
  (P15 already established cx/cy/cz are scaffold). Tests whether ANY
  residual geometric noise matters; doctrinally cleanest frozen form.
  **S3 — the real replacement (kernel-2 design track, "livefab"):**
  the dynamic complex — link/cell creation-annihilation as conversion
  events (CONSTRUCTION.md rewrites), integer-ledger-booked per
  rewrite, deterministic via serial event flush (the QATOM pattern).
  The vacuum SELF-ANNEALS (σ_d → 0 dynamically — disorder removed
  without imposing a grid); structures own their local geometry (the
  M-R1 universality mechanism P15 wanted, without fighting a frozen
  matrix); space-mode books gain their missing face (cell creation =
  space conversion made literal). Prediction worth recording: a
  self-annealed vacuum sits at its own π-rung, d̄_vac = π/w2 = 1.083
  — close to the pinned quant_A0 = 1.15 — so on S3 the action atom
  A₀ = e_s0·d̄_vac/C becomes DERIVED from law constants instead of
  pinned per-box (resolves the recorded A₀ apparatus leak). S1
  results should carry over to S3 (S1 ≈ S3's fixed point, frozen).
  **Sequencing under both directives:** S1 build + battery → re-fight
  the W1 showdown and the flux-balance instrument on quiet geometry
  (MASS-simple: the pin at N=2 upward) → EM5 (law event, doubly
  motivated: Maxwell cone + polarization substrate for charge) →
  CHARGE.md ladder Q1–Q5 → nucleon-class MASS (3-voice machines,
  Q6–Q8, p/n split analog). The old "MASS owns the critical path until
  P19" is amended: for nucleon-class objects the critical path runs
  THROUGH the charge construction.
