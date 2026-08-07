# HORIZON — the forced black hole (user-directed 2026-08-07)

User directive (2026-08-07, condensed): create a modified simulation
that FORCES a black hole in an extreme situation; determine the missing
interaction pattern to sustain it and create the required event
horizon; use a symbolic prover with symbolic backtracing so the answer
is a required STRUCTURE, not tuned constants.

Registered BEFORE the apparatus lands. Companion: `BLACKHOLE.md` (the
five-wall no-go this lane inverts). Prover: `report/horizon_prover.py`
(output transcribed in §1 — every algebraic lemma sympy-checked on the
kernel's actual update forms; the structure search is exhaustive
backtracking over a move grammar from which constant-tuning is
excluded by rule). The forced hole is an ASSERTION INSTRUMENT in the
project's standing sense (frozen gauges, asserted identity): it
imposes the derived structure by fiat to measure what phenomenology
follows and what still refuses. Nothing is adopted; the apparatus is
default-inert and battery-gated; decisions are the user's.

## 1. The prover's verdict: the required structure

**Backtrace (requirement → blocking invariant of the standing law):**

| horizon requirement | blocked by | the sentence |
|---|---|---|
| H1 absorb-always | I1 refusal | intake multiplies (cap−Em−Ee); swallowing closes the door |
| H2 emit-never | I2 atoms | any pitched store fires atoms ε = A₀ω/2π; only ω=0 is silent |
| H2b no-return | I5 unitarity | a one-way boundary over a fixed interior contradicts pairwise rotations; the interior state space must GROW |
| H3 standing far field | I3 π-closure | closed books ⇒ π-flat identically (the Wall-3 theorem) |
| H4 clock stops | I6 clock floor | w/(1+q·x) at bounded x never reaches 0 |

**Minimal satisfying structures (backtracking search):**

* **size 1: ROSTER_DEATH** — cells can die; their capacity leaves the
  books with them. Uncapped (unbounded count), pitchless (no cells =
  no clocks = no atoms = structurally silent), π-invisible (not
  summed), books permanently open.
* **size 2: NEW_LEDGER_X + ROUTE_DOOR_TO_X** — an uncapped, pitchless,
  π-invisible storage mode X, with the door's headroom factor reading
  X's (infinite) headroom instead of Em+Ee.
* **Proven equivalent by property map:** X IS the ledger of consumed
  cells' capacity. **A black hole is the region where cells leave the
  roster** — THEORY.md §2.3's third lane (cell number as a dynamical
  variable), reached by backtracing from the horizon instead of from
  cosmology. Rejected structures, with their killing lemma: UNCAP_EM
  (pitched ⇒ shines, I2), STORE_IN_ES (π-visible ⇒ repels the medium,
  I4 — dπ/dEs = 1 > 0), either size-2 half alone (I1).
* **CLOCK_READS_SPACE (w_e = w·f(Es/e_s0)) is NOT required for
  trapping** (the store carries no clock at all) but IS required for
  the exterior GR analog — nearby matter redshifting with the well.
  Its absence under the standing law is a registered predicted NULL
  this campaign measures (P-HZ6). This is the derived missing
  interaction pattern beyond the store itself.
* **Quantitative (the trapped-surface condition):** inflow speed
  v(r) = Q̇/(2πr·σ_s); horizon where v = c_eff (measured 0.4C):
  r_h = Q̇/(2π·c_eff·σ_s). The pass-S wind ceiling at the standing
  s_k = 0.06 gives Q̇ ~ 2–4/t.u. ⇒ **wind-borne r_h ≈ 1–2 cells**
  (order-of-magnitude); r_h = 4 needs Q̇ ≈ 6.8/t.u. If a larger
  trapped surface appears it must ride FRAME MIGRATION (cell infall
  by contact mechanics, not bounded by s_k) — the probe ladder
  measures which channel carries the river.

## 2. The apparatus (default-inert, both kernels, battery-gated)

`bh_r` (default **0** = byte-inert; > 0 = hole radius about box
centre), `bh_k` (space-drain rate per t.u., default 1.0, read only
when bh_r > 0). One new global accumulator `Eh` — the by-fiat
realization of NEW_LEDGER_X (uncapped, pitchless, π-invisible: it
appears in the conservation sum and NOWHERE else). Per step, every
cell inside r < bh_r, in pass 6 BEFORE its clock would advance:

    Ee → Eh entirely (fa1=fa2=0);  Em → Eh entirely;
    Es → Eh at rate bh_k·(Es − es_floor)·dt  (the floor is respected:
        geometry survives — cells shrink toward r₀·∛es_floor ≈ 0.31);
    then CONTINUE — no clock advance, no beat, no door, no emission.

The eat is one-way by I2's own logic (Eh has no pitch to fire atoms
at); conservation is exact by routing (Eh joins the Kahan total-energy
sum; drift bar ≤ 2e-15 unchanged); π never reads Eh, so the drained
region is a STANDING π deficit the medium answers by its own law
(pass S pushes space in; contact mechanics migrates cells in). Flight
delivered into hole cells is eaten on the next pass-6 visit (≤ 1 step
residence — recorded). Hole cells remain full members of geometry,
field hops, and space transport: only their clock and door are gone —
pitchlessness IS the assertion. Header line
`# HORIZON apparatus (HORIZON.md): bh_r=0 bh_k=1` both kernels;
battery purity pin extended.

What the assertion is NOT: a law candidate. The lawful form of the
derived structure is ROSTER_DEATH (cells actually dying), which stays
design-only per THEORY.md §2.3 until the user opens that lane.

## 3. Arms (registered; L=64 2D law-regime foam, seed 20260802,
k_rad=0.05 p_rad=4 rad_clock=0 par_tau=10, unless noted)

| arm | config | question |
|---|---|---|
| HZ-0 | quiet foam + hole (amp=0, bh_r=4 bh_k=1, T=1500) | the pure space-eater: well depth, river rate, frame infall |
| HZ-X | EXTREME: centred implosion onto the hole (srcx=32, σ=16 sy=16 amp=2 imp_k=2 kx=0, bh_r=4 bh_k=1, T=1500) | maximum-drive feeding; Eh growth curve; census of survivors around a feeding hole |
| HZ-C | HZ-X with bh_k=0.01 | sub-threshold control: phenomena must scale with consumption, not with the apparatus's existence |
| HZ-P6/P10/P14 | hole on + outward probe packet at r∈{6,10,14} (σ=2 sy=2 amp=4 kx=2, T=250, snap 2.5 t.u.) | the trapped surface: escape vs turnaround; centroid-speed deficit vs paired bh_r=0 controls |
| HZ-M | 3D L=24 blob2 sep=12 kx=0 amp=1.6 σ=2.2 + hole bh_r=2 bh_k=1 T=1500 (+ bh_r=0 control) | does matter at range 6 fall in? |

## 4. Bars and registered predictions

* **P-HZ0 (gate):** defaults byte-inert (full battery ALL GREEN 93;
  standing logs move by the new header line only); C≡Go physics
  byte-equal with the hole firing; drift ≤ 2e-15 with Eh in the sum.
* **P-HZ1 (one-way):** Eh(t) strictly monotone; zero QATOM emission
  events from inside r < bh_r; eaten-channel split reported.
* **P-HZ2 (the well — the I3-converse made real):** a STANDING π
  deficit: interior π ≈ es_floor-class vs bulk ≈ 1.01–1.05
  (Δπ ≈ 0.9–1.0) sustained all run — the programme's first standing
  medium disturbance from a stationary object, exactly because its
  books never close; exterior π(r) monotone over several cells
  (contrast: ASTRO flat ±2×10⁻⁵–2×10⁻⁴).
* **P-HZ3 (the river, forked — both outcomes are results):** (a)
  sustained inward space wind: after the interior reservoir empties,
  eatS keeps growing at the wind-delivery rate (predicted Q̇ ~ 2–4
  /t.u. at r=4); (b) frame migration: live-cell count inside r < 6
  rises as cells fall in — OR contact repulsion jams the river (the
  substrate's degeneracy pressure — a Chandrasekhar-class choke,
  which would itself be a discovery).
* **P-HZ4 (the trapped surface):** probes at r ≥ 6 ESCAPE (wind-borne
  r_h ≈ 1–2 cells only) UNLESS frame migration extends the advection;
  any centroid-speed deficit vs the paired no-hole control = the
  first measured pull on light.
* **P-HZ5 (matter at range, forked):** (a) B2 holds — blobs stay
  parked despite the standing well ⇒ the first clean SPLIT of B2
  (no translation) from B6 (no field): a real force with no infall;
  (b) frame migration carries the blobs — the first gravity-like
  infall in the programme.
* **P-HZ6 (the redshift null — the derived missing interaction):**
  QATOM lines of matter near the hole track LOAD only (w_e reads x,
  never Es): no line shift vs range. The absence demonstrates
  CLOCK_READS_SPACE as the remaining structure for the exterior GR
  analog — registered as a future law-candidate shape, user-gated.
* Nothing adopted; bh_r=0 default; decisions the user's.

## 5. Results (all arms harvested 2026-08-07; logs `runs/horizon/`,
fcs local; instruments `report/horizon_prover.py`, `report/analyze_hole.py`)

* **P-HZ0 PASS:** battery ALL GREEN 93 (`runs/BATTERY_horizon.log`);
  standing logs move by the one new header line only; C≡Go physics
  byte-identical with the hole firing (RESULT hole lines byte-equal;
  drift-column-only divergence, and the Go drift is unchanged by the
  hole: 1.99e-14 defaults vs 1.97e-14 firing on the check config).
* **P-HZ1 PASS (one-way, conserved):** Eh strictly monotone in every
  arm — HZ-0 401.07 (pure space: eatF=eatM=0 exactly, the quiet foam
  is field/matter-dark as designed), HZ-X 1150.12 (492.8 field +
  284.2 matter + 373.1 space — the implosion fed it), HZ-M 1400.1.
  Drift with the Kahan-compensated ledger: HZ-X −1.9e-15, HZ-C
  −1.3e-15, both L=24 arms 0.000e+00 EXACTLY, probes ≤1.3e-15;
  HZ-0 −4.8e-15 (2.4× the standing 2e-15 wording — floor-class,
  recorded). No emission from the hole: structural (eaten cells
  never beat) and no in-hole QATOM events observed.
* **P-HZ2 LANDS — THE WELL EXISTS AND REACHES FAR.** π(r) about the
  hole, HZ-0 at t=1500: 0.050 / 0.051 / 0.314 / 0.552 / 0.734 /
  0.863 / 0.948 at bins r<2 … 16–24 — a MONOTONE STANDING gradient
  spanning ≥ 20 cells, deepening all run (bin 4–6: 0.602 → 0.314;
  bin 16–24: 1.000 → 0.948 — the front still propagating at
  t=1500), against the ±2×10⁻⁴ flat of every lawful object ASTRO
  measured. The I3-converse is real: permanently-open books source
  the programme's first standing far field. B6's phenomenology
  exists the moment the derived structure exists.
* **P-HZ3 — BOTH FORKS, in sequence: the river runs, then JAMS.**
  Frame infall: cells inside r<4: 36 → 187 by t≈250 → plateau ~200
  (5.6× densification, held by contact pressure — the registered
  Chandrasekhar-class choke: JAMMING IS THE SUBSTRATE'S DEGENERACY
  PRESSURE); the r<6 shell keeps accreting all run (79 → 318 HZ-X,
  → 279 HZ-0). Sustained consumption after the reservoir: ~0.15
  space/t.u. (Eh 148.6@t=98 → 400.9@t=1498, gently declining) —
  ~20–40× below the naive s_k wind ceiling: the jam and the
  shrunken-contact conductance throttle the river.
* **P-HZ4 — NO TRAPPED SURFACE; the eater is an absorber, not a
  horizon.** Law-regime probes are void: the B7 cloud-chamber foam
  eats the packet in ~20 cells in HOLE AND CONTROL arms alike (a
  measured opacity limit — horizon optics is unaskable in the law
  regime at this scale). The optics-regime ladder (e_cond=99, field
  flies free, hole still eats) answers cleanly: capture fraction
  falls GEOMETRICALLY with launch radius — 72% (r=3, launched
  INSIDE the eating radius) / 28% (r=6) / 13% (r=10) of packet
  energy eaten — while every surviving remnant ESCAPES at full
  speed, including from r=3. No turnaround, no drag resolvable
  above the selective-eating artifact (hole-arm lobes read FARTHER
  out because their inner tails are eaten). The trapped-surface
  condition (advection ≥ c_eff = 0.4C, prover: Q̇ ≥ 6.8/t.u. for
  r_h=4) is missed by the measured river ~40×.
* **P-HZ5 — MATTER DOES NOT FALL: the first clean B2/B6 split.**
  HZ-M blob separation over 1500 t.u.: 11.98 → 12.77 (hole) vs
  11.98 → 12.56 (control) — bodies at range 6 from a standing
  Δπ ≈ 0.9 well do not approach it (the settle-era dip 11.56@150
  recovers, exactly the COMBINE settle class). Yet the hole ate
  52.3 units of MATTER — the blobs' episode halos drained in
  (identity population 2000 vs control 3400): **accretion without
  infall — the hole eats the atmosphere, never moves the body.**
  A real force-field finally exists (B6 phenomenology) and B2
  (nothing translates) stands alone, measured in isolation for the
  first time.
* **P-HZ6 — redshift null, statistics-limited:** HZ-M D-line
  medians 1.582 (hole) vs 1.612 (control), n = 5/4 fires —
  no resolvable shift, consistent with the registered null; the
  structural fact stands (w_e reads x only; no Es term exists in
  the pitch law) — CLOCK_READS_SPACE remains the derived missing
  interaction for the exterior GR analog.

## 6. Verdict (nothing adopted; bh_r=0 default; decisions the user's)

1. **The forced hole produces the black-hole EXTERIOR and refuses
   the black-hole BOUNDARY.** With the prover's derived structure
   asserted (uncapped, pitchless, π-invisible store), the substrate
   immediately delivers: one-way absorption at exact conservation;
   the programme's first standing, monotone, ≥20-cell far field;
   frame accretion (36→200 cells; the packing front); accretion of
   matter's halos; a sustained river. What it will not deliver at
   the standing law is a TRAPPED SURFACE: light escapes from inside
   the eating radius at full speed, and matter does not fall at any
   radius.
2. **The horizon's missing structure is the same structure, needed
   harder.** The trapped surface requires the wave-carrying medium
   to advect inward at ≥ c_eff; the space channel is ~40× too slow
   at the measured river, and the frame channel CHOKES on jamming —
   immortal cells pile up at contact pressure (the substrate's
   degeneracy pressure) and stall the collapse. The prover's
   size-1 solution — ROSTER_DEATH — is therefore needed TWICE: as
   the lawful form of the store (sustainment) and as the jam-breaker
   (cells must LEAVE at the center for the river to run supersonic).
   **In this substrate the event horizon is the surface where the
   roster shrinks at wave speed** — cell death at rate
   Q̇ ≥ 2π·r_h·c_eff·σ_s. One structure closes both requirements;
   it is THEORY.md §2.3's third lane, now with its acceptance
   number.
3. **Independent walls isolated by the assertion:** B2 survives a
   real well (translation is missing physics of matter, not absence
   of force-fields — 'translation IS the current', S2); exterior
   redshift needs CLOCK_READS_SPACE (pitch must read local space —
   a law-candidate shape, user-gated); law-regime opacity (B7)
   forbids long-range horizon optics below the conversion scale —
   any horizon experiment must run its optics in the transparent
   regime or above the work-function threshold (WORKFN.md would
   change this arithmetic).
4. Instrument notes: Ee-weighted lobe tracking is biased outward by
   selective eating (recorded); box wrap contaminates probe tracks
   past t ≈ 90 (L=64, 0.4C); the sub-threshold control HZ-C at
   bh_k=0.01 still floors the interior over 1500 t.u. — a slower
   drain, not a qualitative control (its real contrast: nin 106 vs
   203, Eh 694 vs 1150).
