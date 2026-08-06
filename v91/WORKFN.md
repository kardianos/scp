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

(registered empty)

## 5. Verdict

(after all arms; adoption, as always, is the user's)
