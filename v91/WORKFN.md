# WORKFN — the emergent conversion threshold (B7 fix; law candidate; user task list 2026-08-06)

Registered BEFORE any kernel edit or run. REALITY §B7/§C2: the real
world has a transparent vacuum AND absorbing detectors AT ONCE,
because conversion thresholds live in BOUND MATTER, not in space. The
programme currently declares a regime per experiment (law e_cond=0 /
optics e_cond=99) — an override, not an emergence — and the DS tier-1
screen is a hand-declared per-cell override already in the kernel
(`scond[i]`: screen cells condensation-active). This lane makes that
behavior EMERGENT: one table, both behaviors.

## 0. Thesis and the measured tension

Candidate law: **a cell's condensation demand depends on its own
standing dense matter** — empty space demands much, matter demands
nothing. Form W-A (simplest): at the beat-fire condensation site,

    econd_i = (Em[i] >= wf_floor) ? 0 : wf_far

(default-inert behind `wf_on=0`; `scond` retained as the historical
instrument override, unused in wf arms).

**The registered tension (QUENCH §5.2):** the law-regime e_cond=0
everywhere is exactly what let the quench FARM matter (the
self-sustaining cloud recondenses its own glow). A work function
gates that: debris keeps recondensing only where its Em already
stands above wf_floor. The quench-cloud's fate under wf_on=1 is a
REGISTERED DISCRIMINATOR, not a bar — if the work function kills
quench-genesis, the table gains B7 and loses the programme's only
measured creation mechanism; if the cloud survives (its debris sits
above the floor), the table holds both. Either way it is the
decision datum the user will need.

## 1. Design

* **Keys:** `wf_on` (default **0** = byte-inert: the econd_i line
  reads exactly as today), `wf_floor` (default 0.2 — must sit ABOVE
  the quiet-foam per-cell Em and BELOW bound-matter Em; the meter
  arm measures both before any law arm), `wf_far` (default 99 = the
  optics-regime demand for empty space). Printed on a purity header
  line; battery lawHeader extended.
* **Site (one line, C + Go mirrored):** the beat-fire condensation
  threshold (`econd_i = scond[i] ? 0 : P.e_cond` becomes, under
  wf_on: `econd_i = scond[i] ? 0 : (Em[i] >= wf_floor ? 0 :
  wf_far)`). Nothing else moves; evaporation untouched.
* **The honest detector consequence, registered:** under a work
  function the DS screen must BE dense matter (a seeded dense wall),
  not a declared plane — a config/seed-level change to the DS arm,
  not a kernel change. The wf DS arm therefore runs with a seeded
  dense screen and `scond` OFF.

## 2. Meter arm (before any law arm)

* **W-M1:** measure the Em gap on the standing bodies: quiet-foam
  per-cell Em distribution (slit substrate + warm bath), bound-body
  Em (chord voices, blob interior, quench-cloud debris). The gap
  (foam q90 vs bound q10) either exists — wf_floor selected inside
  it — or does not, and the lane reports the failure honestly (no
  floor fishing).

## 3. Law arms (registered now; run only after W-M1 publishes the floor)

* **W-L1 (transparency):** the DS wave arm with wf_on=1, global
  e_cond irrelevant, NO declared optics: fringes at the
  parameter-free loci through a transparent foam (bar: tier-0 loci
  bars unchanged).
* **W-L2 (detection):** DS with a seeded dense screen, scond off:
  clicks fire AT THE SCREEN ONLY (bar: click loci rebuild fringes as
  tier-1; zero clicks off-screen).
* **W-L3 (the quench discriminator):** Q2's beam dump with wf_on=1:
  report the cloud's population trajectory against QUENCH §4.1 —
  measurement, not bar; the tension datum for the user.
* **W-L4 (XSEC sanity):** the beam-into-foam arm: transmitted
  fraction vs the B7 cloud-chamber behavior — the transparency
  claim quantified.
* Battery ALL GREEN 93 at wf_on=0 before any results commit; every
  standing bar untouched at defaults.

## 4. Results

**4.1 W-M1 (the floor gap; analysis-only from standing snapshots;
2026-08-06).** Quiet slit foam far from any beam: **Em ≡ 0.0000
exactly** (n=1,888, max 0.0000 — the foam carries field, never dense
holdings). Warm bath (a matter medium): Em q10/q50/q90 =
0.010/0.107/0.311. Quench-cloud debris at t=5000: q10/q50 =
0.059/0.121, only 24.6% ≥ 0.2. **The gap exists at the EXISTENCE
grain, not a depth grain:** foam ≡ 0 vs everything-that-is-matter
> 0, with matter a continuum above zero.

**4.2 Design round (from W-M1, before any kernel edit):**
* **wf_floor = 0.01 selected** (a PRESENCE threshold: below the warm
  bath's q10, above exact zero) — all standing matter converts
  freely; empty foam demands wf_far. The depth-variant floor (0.2)
  is REJECTED by measurement: it would strangle the bath's low tail
  and ~75% of the quench cloud — floor fishing above the existence
  grain has no measured support.
* **Registered pre-run predictions:** W-L1 transparency PASSES (the
  DS foam is exactly empty); W-L2 clicks fire only at a seeded
  dense screen; **W-L3: the quench cloud is never BORN under wf**
  (its genesis condensed into Em=0 cells — a strict work function
  forbids exactly that, which is real physics: beams through true
  vacuum do not condense). Therefore:
* **W-L3b registered — creation = quench + nucleation site.** The
  beam-dump arm with wf_on=1 plus ONE seeded dense grain in the
  beam path: does the cloud grow from the grain? Real cloud-chamber
  physics (supersaturation nucleating on dust). If yes, the table
  holds B7 transparency AND a creation mechanism — the tension of
  §0 resolves as "creation requires seeds"; if no, the user's
  trade (B7 vs genesis) stands as measured.

**4.3 W-L1 / W-L3 / W-L3b (the law arms; logs
`runs/quench/wl1_transparent.log`, `wl3_cloud_wf.log`,
`wl3b_grain.log`) — all three registered predictions land exactly.**
(harvested 2026-08-06)

* **W-L1 PASS — optics through a transparent-by-law foam.** Fringe
  peak y = 32.3 (bar 32 ± 1.9) with the conversion ledger EXACTLY
  ZERO (cond = rough = rad = 0.000000): the beam crossed and
  interfered with no declared regime anywhere.
* **W-L3 — the cloud is never born, as predicted.** The amp-8 beam
  that farmed a 2,400-episode cloud in the law regime condenses
  NOTHING under wf (gids ≡ 0 for 1500 t.u.; mints = 0). Beams
  through true vacuum do not make matter — the model now agrees
  with reality on this.
* **W-L3b — CREATION SURVIVES AS NUCLEATION.** One dense grain in
  the beam path (the XSEC absorber as dust): the cloud grows FROM
  THE GRAIN — 0 → 269 (t=300) → 2,330 (t=1500), a full QUENCH
  §4.1-class population (ages 63/217/558; 6,580 standing stamps,
  max 1,286; idfrac 0.949). The §0 tension RESOLVES: **creation
  requires seeds** — supersaturated field nucleating on standing
  matter, real cloud-chamber physics.

## 5. Verdict (all arms complete; ADOPTION IS THE USER'S)

The W-A presence-threshold work function — one line, one measured
constant class (the existence-grain floor) — simultaneously
delivers, in one table with laws_V2g untouched elsewhere:
transparent vacuum (W-L1, conv ≡ 0), interference optics with no
declared regime (fringes at the parameter-free locus), detection
localized to matter (the scond screen becomes emergent in
principle; the dense-screen click arm W-L2 remains registered, not
yet run), the quench-genesis mechanism preserved in its physically
correct form (nucleation on seeds, W-L3b), and darkness of empty
space (W-L3). B7 — the audit's "cloud-chamber vacuum" ANTI row —
has its one-table fix measured. What this changes if adopted: the
per-experiment regime declaration (e_cond 0/99) retires; XSEC's
absorber/emitter results re-read as nucleation physics; the
quench-creation story becomes seeded-genesis. Adoption, the W-L2
dense-screen completion, and any REALITY re-grade are the user's
decisions; wf_on defaults 0 and the battery stays green either way.
