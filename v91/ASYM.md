# ASYM — the asymmetry-driver pad campaign (user-directed 2026-08-07)

User directive (condensed): pursue flow-gravity; the largest change is
cell structure that makes asymmetry SURVIVABLE; the driving mechanism
is something already common to every cell but currently driven FROM
the simulation harness — it must be taken INTO the simulation so the
complete live state can drive it. Give 18 ideas including crazy ones;
then hit your head (metaphorically) against an RNG for more; test
them; refine each 3 times; then conclude.

Protocol (pad precedent, `pad/RESULTS.md`): each idea gets T0
(one-law/no-background coherence), T1 (the measured record), T2 (sim
probe where config-only reachable), and three refinement rounds — R1
steelman/fix, R2 confront the hardest standing fact, R3 reduce to the
minimal implementable form. Grades: DEAD / WEAK / FOLD(→n) / ALIVE /
STRONG / CORE. The RNG round is literal and recorded: python
`random.Random(20260807)` colliding {harness facilities} × {cell
attributes} × {verbs}; collisions verbatim in §2. Probe logs
`runs/horizon/pa*.log`. Nothing here is adopted; this campaign's
output is a converged design for user sign-off.

The inventory of things common to every cell but harness-driven
today: the global dt; the RNG stream (seed + in-run tumble draws);
the canonical link order; the pass splitting and its fixed order; the
Kahan compensation registers; the jam/contact solver; the candidate
(adjacency) rule; the periodic fold; the seed/init; the halting T;
the config table; thread domains; the float mantissa; the
snapshot/diag meters; the credit registers' SHAPE (shared scalars);
the registry/identity ledgers (meter-only). Eighteen ideas from that
inventory:

## 1. The eighteen

* **A1 scheduler-as-physics.** Update cadence is half-internal
  already (beats); internalize PRIORITY: contention resolves toward
  the hungrier cell — "who acts first" is flow direction.
* **A2 RNG-as-state.** The tumble noise is the only in-run RNG use;
  replace the external stream with a chaotic map of the cell's own
  complete live state — fluctuations become state-correlated, which
  is the one thing a ratchet needs and i.i.d. noise can never give.
* **A3 dt-as-proper-time.** Retire the global step; every cell
  advances on its own clock (the P2 scheduler promoted from
  infrastructure to law).
* **A4 order-as-state.** The canonical link order becomes
  state-ordered (largest-gradient-first).
* **A5 Kahan-vacuum.** The compensation registers become a physical
  sub-grain reservoir (the vacuum's loose change).
* **A6 stress-as-field.** The jam solver's contact forces become
  per-cell standing stress state (prestress) — stored directional
  force.
* **A7 mortal link roster.** Adjacency stops being a geometric
  recomputation: links are born, maintained, and die as state —
  and their formation follows delivered flow (rivers grow channels).
* **A8 physical boundary.** The fold becomes a lawful edge
  (inside/outside asymmetry source).
* **A9 continuous genesis.** The seed becomes an ongoing door
  (creation as physics, not init).
* **A10 meter back-action.** Observation disturbs (measurement as
  physics).
* **A11 vector credit.** The atoms credit registers — today shared
  per-cell SCALARS (R-D4 left their shape open) — become
  per-DIRECTION: a cell remembers which neighbor it owes. Standing
  credit dipoles = asymmetry stored BETWEEN atoms, below the
  quantum, without touching atom counting.
* **A12 state-ordered passes.** Operator order per cell driven by
  local state.
* **A13 constants-as-slow-fields.** s_k (etc.) locally dynamical —
  conductance responding to what flows.
* **A14 internal halting.** The run ends when the substrate says
  (heat death as law).
* **A15 per-cell seed.** Every cell carries its own hidden RNG state
  (character, not just name/gid).
* **A16 grain rounding with a sign.** Sub-ε rounding direction
  biased toward the flow (round-toward-the-debtor).
* **A17 thread-domains-as-causal-patches.** Parallel decomposition
  becomes physical patch structure.
* **A18 registry-as-law (flow hysteresis).** The exchange registry —
  today a physics-silent meter on every slot — becomes load-bearing:
  the standing NET on a link biases its future conductance. The
  river digs its own bed.

## 2. The RNG round (seed 20260807, collisions verbatim)

```
RNG-1: RNG-stream x links x mirrors      RNG-5: canonical-link-order x clock x rectifies
RNG-2: diag-meter x pitch x mirrors      RNG-6: float-mantissa x es-floor x orbits
RNG-3: candidate-rule x flight x delays  RNG-7: periodic-fold x beat-phase x orbits
RNG-4: credit-register x links x remembers  RNG-8: candidate-rule x tumble-planes x rectifies
```

* **R1 conservative noise.** Noise mirrored across links: every
  fluctuation is an EXCHANGE (one endpoint's gain is the other's
  loss) — the lawful form of A2: internal thermal micro-currents,
  exactly conserved, no injection.
* **R2 mutual meters.** Cells measure each other — and already do:
  the consonance gate IS a two-sided meter. (A10's lawful half
  exists; the collision closes it.)
* **R3 retarded adjacency.** The candidate rule gains flight delay:
  links form and die on c-lagged information — the contact rule made
  causal (a B8 repair).
* **R4 credit remembers links.** The RNG independently re-derived
  A11 verbatim — recorded as a two-hit coincidence and a vote.
* **R5 phase-ordered execution.** Cells act at their own clock-hour:
  the update sweep becomes a rotating wave — a CHIRAL substrate
  (parity violation for free; wound matter would finally have a
  preferred handedness).
* **R6 the breathing floor.** es_floor oscillates lawfully — a
  parametric pump on the space mode (Casimir-flavored vacuum drive).
* **R7 box-winding sectors.** The torus holds quantized global
  circulations of beat phase — cosmological angular momentum as
  topology.
* **R8 nematic tumble.** Tumble rectified by adjacency: planes align
  with neighbors — spontaneously ORIENTED cells, i.e. each cell
  grows the anisotropic tooth a ratchet needs.

## 3. Tests and three refinement rounds each

**A1 scheduler.** T0 clean. T1: P2 measured batch==serial EXACT at
1–8 workers — the current law is deliberately order-robust; beats
already internalize cadence. R1: the content left is PRIORITY under
contention. R2: hardest fact — nothing contends today (exactness
proves it); priority is inert without a contended resource. R3:
minimal form = pending-min key (time, want): fires only where wants
collide. Grade: ALIVE (needs A7/A18 to create contention).

**A2 internal noise.** T0: better than clean — the RNG stream is a
BACKGROUND the no-background rule never audited: an immortal,
extra-physical state influencing every cell. T1: in-run use is
tumble only; Probe A (§4) decides if it is load-bearing. R1:
steelman — i.i.d. external noise has ⟨asym⟩ ≡ 0 by symmetry; only
state-correlated fluctuation can rectify. R2: hardest fact —
Feynman: rectifying equilibrium noise breaks the second law; lawful
ONLY because the substrate is driven (radiance holds cells
off-equilibrium; the chord's standing circulation is measured lawful
rectification in curl form). R3: minimal form = R1's conservative
exchange-noise from a chaotic map of (th2, cbeta, fa1, fa2) — byte-
switchable, exactly conserved. Grade: STRONG (probe pending).

**A3 proper time.** R1: simultaneity artifacts vanish. R2: the
π-theorem is per-event, not per-time — retiming cannot open books;
dt is not where asymmetry dies. R3: fold. Grade: FOLD(→A1).

**A4 state order.** T1: C≡Go byte-identity REQUIRED canonical order;
P2 proved order-robustness — order is hygiene, not hidden physics.
R1–R3 collapse to contention priority. Grade: FOLD(→A1).

**A5 Kahan vacuum.** T0 VIOLATION: physics from representation
residue makes the law precision-dependent (the measured C↔Go 1-ulp
envelope would become physical divergence). R1: the honest content
is "a sub-grain reservoir" — which is A11's credit, explicit and
lawful. R2: floors sit 13 orders below atoms. R3: fold the metaphor,
kill the mechanism. Grade: DEAD (→A11 carries the idea).

**A6 stress field.** T1: v89 PRESTRESS — 43 jobs, all particle
claims null; plasticity topology-dependent; HORIZON measured
jamming as the river-choke. R1: stress is a SURVIVAL medium, not a
driver (stored force ≠ standing flow). R2: hardest — overdamped
geometry relaxes unless topologically locked. R3: minimal = the
existing bond d-walk given hysteresis (plastic rest lengths).
Grade: ALIVE (survival-half only).

**A7 mortal link roster + flux-grown channels.** T0 clean (slots
already die with endpoints — no immortal index). T1: the user's own
opening clause ("geometry should create density and flow") is this
idea; FORGE E3 measured forging-follows-density; slot machinery
exists. R1: a cell with more upstream than downstream links IS a
standing dipole — asymmetry survives in the GRAPH. R2: hardest —
three ignitions of self-growth gates; growth must not feed on
symmetric churn. R3: minimal = conductance redistribution
proportional to signed NET flow under a fixed per-cell budget
(zero-sum: growth here forces shrink there — ignition structurally
excluded because total conductance is conserved). Grade: CORE
(merged with A18 below).

**A8 boundary.** Cosmology-grade, not cell-grade asymmetry. Grade:
WEAK (E-G/cosmology lane note).

**A9 genesis.** WORKFN W-L3b already measured creation-as-
nucleation at the door — genesis IS door physics. Grade: FOLD(→
standing WORKFN lane).

**A10 meter back-action.** T0: violates the programme's central
hygiene (physics-silent meters — every campaign leans on it). R2:
RNG-2 shows the lawful half already exists (gates are mutual
meters). Grade: DEAD by discipline, closed honorably.

**A11 vector credit.** T0 clean (conserved bookkeeping). T1: R-D4
explicitly left register shape open; the README instrument lesson
("sub-grain demands accrue credit silently") shows the sector is
real; RNG-4 re-derived it independently. R1: the sub-atom ledger is
the one place a standing dipole can live WITHOUT touching atom
counting (ħ-linearity safe). R2: hardest — two-atom lapse bounds the
memory; fix: lapse per direction, not per cell. R3: minimal = split
qcnvD/qcnvF per incident link, per-direction lapse; acceptance =
every atoms/XSEC bar unchanged + standing credit dipoles measurable
at chords. Grade: STRONG.

**A12 state-ordered passes.** Same rock as A4; operator-order
robustness was engineered and verified. Grade: DEAD.

**A13 s_k as field.** T0: local law drift = background-flavored
unless the "constant" is lawful state — at which point it is a
conductance, i.e. A7/A18. The HORIZON threshold (s_k ≥ 0.119 for a
space-channel horizon) says WHY one would want it; the budgeted
channel law says HOW lawfully. Grade: FOLD(→A18).

**A14 internal halting.** No cell-grade asymmetry content. Grade:
DEAD (kept as a cosmology footnote).

**A15 per-cell seed.** A hidden immortal per-cell stream = a new
background; the lawful version (character = the cell's own live
state) is A2 verbatim. Grade: FOLD(→A2).

**A16 signed grain rounding.** Transport is never quantized (law);
rounding lives only at doors — where it is A11's directional credit
again. Grade: FOLD(→A11).

**A17 causal thread patches.** Deterministic-OpenMP was built
precisely to make domains invisible (byte-identical at any thread
count, measured). Physics from decomposition would un-measure that.
Grade: DEAD by measured discipline.

**A18 registry-as-law (with A7: THE BED-DIGGING CHANNEL LAW).** T0
clean (slot-borne, mortal). T1 support: per-link balance emerges
with the ledger window (ρ 0.49/0.68/0.79); the chord runs net =
0.75·gross lawfully; T1 danger: the ignition family (three strikes)
— any flow-begets-flow gate must not eat churn. R1: gate on SIGNED
NET, not closure/balance — bath churn is symmetric at long windows
(its net is the small 1−ρ tail), so net-gating starves churn
naturally. R2: hardest — REGISTRY measured bath standing net ≈ 0.21
at τ=100: nonzero food; amplification must be bounded. The fix is
structural, not a threshold: per-cell conductance BUDGET (Σw fixed,
redistribution zero-sum) — global ignition is impossible because
there is nothing global to grow. R3: minimal implementable form:
w_l ← w_l + k·net_l − (renormalize to budget), slot-borne, dying
with the link; acceptance surface: (i) battery green at defaults,
(ii) the warm bath does NOT channelize (zero-sum + symmetric churn),
(iii) a FED region channelizes (HZ-0-class eater grows persistent
inflow channels), (iv) two-body: a pair of eaters grows a channel
BETWEEN them = the first attraction. Grade: CORE.

**RNG-1 conservative noise.** Absorbed as A2's R3 form. CORE-half.
**RNG-2 mutual meters.** Closure note; no new law. FOLD.
**RNG-3 retarded adjacency.** Real B8 repair, orthogonal to
asymmetry; registered for the relativity lane. ALIVE (parked).
**RNG-5 phase-ordered execution.** The crazy keeper: a rotating
update wave = global chirality; would give wound matter a preferred
handedness (parity violation). Conflicts with P2 exactness — needs
the contention framework first. ALIVE (crazy, parked behind A1).
**RNG-6 breathing floor.** Parametric vacuum pump; conservation-safe
only as redistribution; adjacent to the E-F wall apparatus. WEAK.
**RNG-7 box winding sectors.** Topology footnote for cosmology.
WEAK.
**RNG-8 nematic tumble.** The ratchet's TOOTH: neighbor-aligned
planes = spontaneous orientation = per-cell anisotropy that A2's
drive can rectify against. Cheap form: tumble bias toward the local
delivered-flow direction. ALIVE (the third leg of the mechanism).

## 4. Probe results

**Probe A — is the harness RNG load-bearing for standing asymmetry?
NO — decisively.** (warm radiance bath L=16 T=400, reg_tau=10;
sigma_tumble ∈ {0, 0.01, 0.05}; logs `runs/horizon/pa*.log`.)
Per-link balance quartiles ρ[q25/q50/q75/q90] at T=400:

| arm | ρ quartiles | flowing | cond | drift |
|---|---|---|---|---|
| RNG OFF (fully deterministic) | 0.090 / 0.468 / 0.841 / 0.959 | 0.924 | 2333.8 | −1.3e-15 |
| table (0.01) | 0.094 / 0.473 / 0.840 / 0.959 | 0.922 | 2285.2 | −1.5e-15 |
| 5× (0.05) | 0.097 / 0.488 / 0.840 / 0.957 | 0.929 | 2064.8 | −1.1e-15 |

The standing exchange-asymmetry structure (1−ρ at every quartile,
the flowing fraction, the whole shape) is INVARIANT under removing
or quintupling the only in-run RNG use — differences sit at the
chaos band (~1–2%), while 5× noise does shave conversion totals
~10% (energy structure, not asymmetry structure). **The substrate
never needed dice: with sigma_tumble=0 the run is fully
deterministic and its churn still carries standing per-link
asymmetry 1−ρ50 ≈ 0.53 at τ=10.** The drive is already internal —
self-generated deterministic chaos. Internalizing the RNG (A2/RNG-1)
is therefore a pure BACKGROUND PURGE — the last unaudited immortal
state influencing every cell — free of physical cost, and it is NOT
the asymmetry driver.

## 5. Convergence and conclusion

The twenty-six ideas partition cleanly under test:

* **The noise family (A2, A15, RNG-1) is not the driver** — measured
  (Probe A). Keep the conservative internalization as hygiene.
* **The order/scheduler family (A1, A3, A4, A12, RNG-5) is not the
  driver** — measured (P2 batch==serial exactness was engineered and
  gated). Its content — contention priority, chirality — is real but
  parked until a contended resource exists.
* **The representation family (A5, A16, A17) is dead** by the
  programme's own measured disciplines.
* **The riddle's answer is the LEDGER: the mechanism common to every
  cell, already computed every step, and driven from the simulation
  is the exchange bookkeeping itself** — the registry's per-link
  net flow and the credit registers' shape. The simulation knows,
  at every instant, which way every link flows; the physics never
  reads it. Taking it into the simulation so the complete live
  state is used = **the bed-digging channel law** (A18+A7 merged,
  CORE): per-link conductance redistributed toward SIGNED NET flow
  under a fixed per-cell budget — slot-borne, mortal, zero-sum.
  Zero-sum is the anti-ignition structure (three strikes taught the
  constraint): there is no global quantity to grow, so symmetric
  churn cannot feed it, while a fed flow digs its channel deeper by
  starving its own side-channels. A11 (vector credit) is the same
  law at the sub-atom grain; RNG-8 (nematic tumble) is its
  cell-anisotropy face.
* **The ratchet triad, assembled:** DRIVE = the substrate's own
  deterministic churn (measured present, RNG-independent);
  TOOTH = per-cell anisotropy (RNG-8); MEMORY = budgeted signed
  channels (A18/A7/A11). Survivable asymmetry = channels that
  remember which way they flowed. The largest change, named
  honestly: **stationarity redefined from detailed balance to
  through-flow** — the cell's links stop being symmetric conduits
  and become budgeted, signed, mortal channels.

**Registered acceptance surface for the FLOW lane (awaits the
user):** (i) defaults inert, battery ALL GREEN 93, C≡Go; (ii) the
warm bath does NOT channelize (the anti-ignition bar — zero-sum +
symmetric churn must starve the gate); (iii) a fed eater (HZ-0
class) DOES channelize — persistent inflow channels; (iv) two
eaters grow a channel BETWEEN them — the programme's first
attraction; (v) chord books and the F2 doublet unchanged; (vi)
a_div(r) around a fed body becomes nonzero and ~1/r-class — the
first lawful flow-gravity field. Parked behind it: A1 contention
priority, RNG-5 chirality, RNG-3 retarded adjacency, A6 plastic
stress, RNG-6 breathing floor. Nothing adopted; the ledger and this
conclusion are the deliverable.
