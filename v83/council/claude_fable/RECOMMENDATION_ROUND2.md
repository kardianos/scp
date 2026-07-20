# RECOMMENDATION ROUND 2 — arc verdict, orbit-protocol defect, monist (M) audit

**From:** Claude Fable (independent design critic)
**To:** Program orchestrator
**Date:** 2026-07-20
**Re:** BACKGROUND_ROUND2 Q1–Q8 after hybrid force fix + E1b/E2 option matrix
**Prior memo:** `RECOMMENDATION.md` (Round 1). Stance: skeptical of the
orchestrator's framing and of monist over-reach; kill gates over construction.

---

## 0. Executive summary

1. **Verdict on the arc: progress, but less was measured than the brief claims.**
   The force-sign fix and E1a PASS genuinely unblock the coupling (my Round-1
   "E1a FAIL → stop everything" branch resolved the best possible way: it was a
   sign bug, not an interface gap). But the headline "orbit FAIL 0/31" is
   **over-determined by the protocol, not by physics**. Two independent defects
   confound it: (a) the revs≥0.75 gate was arithmetically out of reach for a
   *perfect* circular orbit at the scan durations used (a clean Kepler circle at
   r=10 accrues ~0.17–0.25 revs in T=280–400); (b) the measured E1a force law
   itself shows the far field is **steeper than 1/r³ beyond D≈8** (local exponent
   n≈4 between D=8 and 10), and by the standard stability criterion for central
   forces (F∝r⁻ⁿ stable circular orbits require n<3) **no seed velocity can
   produce a stable circular orbit there**. "Under-speed → core; over-speed →
   expand" is exactly the phenomenology of an n>3 field. Both defects are
   fixable without new physics. §2 has the arithmetic.
2. **Therefore: do not buy a new mechanism yet.** The brief's §5.6 suggestion
   that orbit FAIL signals an incomplete mass/inertia ledger is premature — the
   baseline E1b run also shed 20% of the ball's charge (Q 315→252) and 20% of
   total energy; no conclusion about missing physics survives that. The cheapest
   decisive experiment is one long, big-box run of the already-identified best
   cell (`B1_r10_vf1p0`-class) with the gate restated in orbital periods and the
   force exponent verified Coulomb-clean at the seed radius first. That is a
   ~3-day kill, and either outcome is decisive (§8, F1–F2).
3. **E2 / 2B: close static multi-center.** C1–C3 double-confirm v71: this theory
   fuses or separates; it does not statically stand off. One cheap kill remains
   that has never been run and is consistent with the only PASSes on the board:
   the **dynamical dinuclear rotor** — two co-phase balls at D≈10–12 (where E2
   itself measured net attraction via merge) seeded with *tangential* velocity,
   so the angular-momentum barrier replaces the missing standoff. If that fails,
   adopt C4 (liquid-drop only) as the standing 2B verdict and write the
   permanent negative (§5).
4. **Monist (M): sound as a research program, currently untested by these
   numerics — in either direction.** The locality-c reading of E=mc² is a real
   physical claim, not decoration, because both sides of the identity are
   independently measurable in the kernel. But the lock implementation does not
   instantiate it: a lock's mass is a config parameter, not a measured bound
   ledger — locks are point dualism with monist labels. So E1b's failure says
   nothing about (M), and E1a's success doesn't confirm it. The **strongest
   in-codebase falsifier** is the inertia identity on genuine fabric objects:
   push a parked Q-ball with a known external force and check
   m_inertial = ∫T₀₀/c² within tolerance, then check that the ~1.6% fusion mass
   defect shows up as an *inertial* defect of the droplet. Cheap, one week,
   bankable either way (§3).
5. **Stage 3 substrate revision (the biggest strategic recommendation):** the
   light orbiting object should be promoted from a registry lock to a **fabric
   anti-ball** — the conjugate-seed Q-ball with Q<0, which the complex kernel
   supports for free. It is ledger-complete (its inertia *is* ∫T₀₀), it tests
   (M) and orbit dynamics in one run, and it is exactly the sanctioned v75
   "Step 2 ± orbit/capture on E-lite" gate: if it fails by bag capture at
   contact, that is the documented evidence CLAUDE.md requires before opening
   the multi-fabric kernel. The lock sector remains a useful toy, not the
   Stage-3 acceptance substrate (§6).

---

## 1. Q1 — Verdict on the arc

**Progress, on three counts; dead end on one narrow thing only.**

- **Coupling is real.** F=+(gq)E with shared Gauss at the 1e-13 floor, attract/
  repel signs correct, Coulomb-class magnitude at D≤8. This was the load-bearing
  unknown from Round 1 (my gap #5) and it is now closed. That is the arc's
  durable asset.
- **The orbit FAIL is not yet a physics result.** See §2. A campaign whose gate
  cannot be passed by the very object it seeks, run in boxes where the far field
  is demonstrably non-Coulomb, has measured its own protocol. The honest
  scoreboard line is "hybrid orbit: **NOT YET TESTED under valid conditions**,"
  not FAIL. I am being deliberately hard on this because the brief's long-range
  section (§4.2) and theory section (§5.6) both build on "orbit FAIL" as an
  established physical fact, and downstream decisions (new co-fields, monist
  ledger doubts, program-freeze risk) are being priced off it.
- **The 2B NULL is real.** Unlike E1b, the E2 result matches an independent
  prior dataset (v71 interlock morals) in sign and morphology across gauged and
  ungauged sectors. Static multi-center is now a twice-confirmed negative. That
  is a legitimate close.
- **For monist Stage 3/4 specifically:** nothing in this arc bears on (M) yet,
  because no experiment in it involved an object whose mass is a fabric ledger
  (§3.3). The arc neither advanced nor damaged monism; it advanced plumbing.

## 2. The protocol defect (the evidence for §1)

Two independent problems, both derivable from the program's own tables.

**2.1 Gate arithmetic.** From the E1a force table (Q≈315, m_lock=2):

| r | F (measured) | v_circ=√(Fr/m) | T_orb=2πr/v | revs available in T=280 | in T=400 |
|---|---|---|---|---|---|
| 8 | 7.22e−4 | 0.054 | ~935 | 0.30 | 0.43 |
| 10 | 2.86e−4 | 0.038 | ~1660 | 0.17 | 0.24 |

Using continuum F instead of measured F shifts these by ≤40% and changes
nothing. The band gate demanded revs≥0.75 **and** late r∈[0.55,1.7]r₀. At
r≥8 and T=280–400, a *flawless* circular orbit scores 0.2–0.4 revs and FAILs.
The only trajectories that can reach 0.75 revs in that time are inspirals
(which spin up as r shrinks) — and the r-band correctly rejects those. So the
matrix's two gate clauses jointly exclude nearly every physical trajectory,
good or bad. **0/31 was close to guaranteed before the first run.**

Corroboration: the "flattest" cell, B1_r10_vf1p0, logged revs=0.34 with
r∈[10.0,10.5] over T=400 — i.e., **r held to ±2.5% while accruing roughly the
revolutions a perfect Kepler circle would accrue in that time.** That cell is
not the least-bad failure; it is plausibly a *passing orbit scored by a
mis-calibrated gate*. It must be re-run long before any new mechanism is
funded.

**2.2 Far-field screening / orbit stability.** Local force exponent n
(F∝r⁻ⁿ) from successive E1a pairs: n≈2.0 (D=5→6), n≈2.5 (6→8), **n≈4.2
(8→10)** — also visible as the measured/continuum ratio decaying 0.85 → 0.74 →
0.46. For central forces, circular orbits are dynamically stable only for
n<3. At D≥8 in these boxes (E1a at N=48 L=12; the B2–B5 scans at L=16, with
seeds at r up to 12), the field is past that bound — presumably box/sponge
screening of E — so **every** circular seed there sits on an instability:
under-speed falls in, over-speed escapes. Which is verbatim the observed
phenomenology. B4/B5 nulls are then expected: bag and soft core act at small
r, and the failure lives at large r. B3's "heavier → flatter" is likewise pure
kinematics, as the eval already suspected.

**Consequences:** (i) the 31-cell matrix cannot support "seed/box/mass/soft/
bag insufficient — needs new physics"; it supports "gate and box were
invalid"; (ii) the fix is one static force-exponent audit in a big box plus
one long run — days, not a mechanism program.

## 3. Q2 — The monist E=mc², taken seriously

**3.1 What (M) actually asserts, operationally.** Strip the metaphysics and
the locality-c reading makes two kernel-checkable claims:

1. **Inertia identity:** for any bound configuration, the resistance to
   acceleration equals the configuration's total field energy divided by the
   square of the free-signal bound: m_inertial = E_bound/c², with c measured
   from the free sector (lattice gauge/wave speed), not assumed.
2. **Warp corollary:** locks must rearrange the free medium such that free
   propagation near mass, described in a global chart at fixed local c, looks
   curved (path cost, Shapiro-like delay).

This is a genuine identification, not a unit convention, precisely because the
two sides are *independently measurable*: E from ∫T₀₀ (already computed), c
from free-wave dispersion on the lattice (already implicit in CFL), and
m_inertial from F=ma with a calibrated external force. Dualism does not force
these three numbers into one relation on a lattice; if they nonetheless agree
to a few percent — including for *composite* objects where binding energy has
been radiated away — that is nontrivial support for reading mass as a ledger
of the same continuum. I engage it as physics, and it deserves the test.

**3.2 What the kernel can and cannot test today.** Claim 1 is testable this
week (§3.4). Claim 2 is **not testable in this codebase**: the ψ path-cost
co-field is explicitly "optional bookkeeping — not in nuclear recipe"
(PHYSICS_RELATIONS R4), and scp_sim evolves on a fixed lattice with no
free-capacity depletion around locks. Be honest in all documents that the
gravity/warp half of (M) currently lives only in the v76/v77 sandboxes. Any
statement that kernel results support the *warp* reading is over-claim.

**3.3 Why this arc says nothing about (M) — in either direction.** The hybrid
light "particle" in E1b is a registry lock: its mass is an input parameter,
its E\* is imported, it has no bound-field configuration whose ledger could be
complete or incomplete. It is, ontologically, a dualist test particle wearing
a monist name tag. So:

- Orbit FAIL does **not** indicate "the mass/inertia ledger of the light lock
  is incomplete" (brief §5.6). There is no ledger to be incomplete. And per
  §2, the failure is anyway protocol-confounded.
- Symmetrically, E1a's force PASS does not confirm monism; a charged test
  particle in a field is the oldest dualist setup there is.

The B4 null (anti-lock bag inert on the hybrid light) is consistent with my
Round-1 position that hand-inserted, ledger-open forces are decoration — but
"consistent with" is all it is.

**3.4 Strongest falsifier of (M) in this codebase.** Ranked; all use existing
machinery, no kernel changes:

1. **Inertia identity on a parked Q-ball.** Apply a weak uniform external E
   (or a fixed distant source charge) to a parked gauged ball, measure a from
   the COM track, compute m_inertial = gQE_ext/a, compare to ∫T₀₀/c². (M)
   predicts agreement; a stable few-percent-calibrated disagreement that
   survives resolution refinement falsifies the identity as stated.
2. **Inertial mass defect.** Fusion defect is ~1.6% of rest energy (R2).
   (M) demands the *inertia* of the fused droplet be reduced by the same 1.6%
   relative to the summed constituents — energy that radiated away must have
   carried its inertia with it. Same push protocol on c6_light vs. free balls.
   This is the sharpest version, because it distinguishes "mass = ledger of
   current field content" from "mass = tag assigned at creation."
3. **Boost scaling.** v82's boosted-ball machinery: measure E(v) for a ball at
   v=0.1–0.4 and check γ-scaling against lattice c. Free locality as *the* c
   in m=E/c² predicts the relativistic form with the measured wave speed, not
   with any nominal one.

If 1–3 pass, (M)'s inertia leg is banked and the locality-c reading graduates
from slogan to measured relation in this codebase. If any robustly fails,
v76/v78 need revision before more monist language enters CONCEPT.md. Either
way, stop letting (M)'s status ride on lock-orbit experiments that cannot
touch it.

## 4. Q3 — Ranked next mechanisms for the "missing second scale"

Preliminary: **the evidence that a second scale is missing is currently weak**
(§2). Rank 0 is therefore not on the brief's list:

- **(0) No new mechanism: fix the far field and the gate, re-run one orbit.**
  Costs 3 days. Everything below is conditional on its failure.

Then, in order:

1. **(iii) Free ± pair with no nucleus — upgraded to the fabric anti-ball.**
   The complex kernel gives Q<0 balls by conjugate seeding at zero new code.
   A small anti-ball (Q≈−100) vs. a ball (Q≈+315) is a fully ledger-closed
   two-body Coulomb system where both masses are real ∫T₀₀ ledgers. It is
   simultaneously: the v75 Step-2 sanctioned experiment, a Stage-3 positronium
   precursor, and a live test of (M)'s dynamics. Its known risk — same-fabric
   bag attraction causing merge/annihilation at contact — is exactly the
   question CLAUDE.md says must be answered before the multi-fabric kernel is
   opened. Highest information per day on the board.
2. **(i) Free capacity co-field.** The theory's *declared* second-scale
   postulate (v82 F_core lineage), with a path to monist semantics
   (sequestration/depletion) rather than remaining a registry hack. Fund it
   **only if** rank-0 and rank-1 orbits fail under clean-Coulomb conditions —
   at which point it has a measured job (regularize close approach, absorb
   radiated energy) instead of being a patch for a screened box.
3. **(iv) Rigid composite C.** As *diagnostic scaffold only*: pinning internal
   DOF separates force-law failures from internal-response failures in one
   control run. Never a deliverable; label it as such in every plot.
4. **(ii) Matter-density surface wall at the nucleus.** Lowest rank of the
   constructive options: hand-localized, ledger-open, and aimed at small-r
   physics when the measured failure is at large r. B4/B5 nulls are directly
   discouraging. Revisit only if a *clean* orbit fails specifically by
   absorption at the ball surface — that failure mode, and only that one,
   gives a wall a motivated job (my Round-1 absorption branch).
5. **(v) Abandon locks / return to multi-fab.** Premature. Multi-fab is
   stopped on its own FAIL, and the sanctioned re-entry condition is precisely
   rank 1's failure-by-bag-capture. Don't reopen it by frustration; reopen it
   by that measurement.

## 5. Q4 — 2B call

**Pause static multi-center permanently; one cheap kill remains, and it is
dynamical.** C1–C3 plus v71 make the static verdict twice-measured: co-phase
merges (net attraction to contact), anti separates, large-D flatness is
kinematic stall. More phase/flavor/D scanning is confirmed low-yield.

The one untried cell consistent with everything measured: **C5, the dinuclear
rotor.** E2 itself shows net attraction at D≈10 (co merges from rest). Give
two co-phase balls tangential velocity from the same corrected circular-seed
protocol as v82-A (F_net=μv²/r, each ball carries v_rel/2), at D∈{10,12},
angular momentum supplying the barrier that statics lack.

- **PASS gate:** ≥2 relative revolutions with 2 charge clusters retained,
  per-cluster Q stable to <2%, Gauss floor. → 2B exists as a *rotational*
  multi-center state (physically respectable: real nuclear molecules are
  rotational resonances), and retained-A carbon becomes a rotor-chain
  question.
- **FAIL gate:** merge within 1 rev at both D, or separation. → Write the
  permanent negative: *retained-A multi-center carbon does not exist in this
  kernel, statically or dynamically, at g=0.05*; adopt C4 (liquid-drop only)
  as the standing 2B answer — which, per my Round 1, costs the atom program
  nothing, since no Stage-4 observable consuming retained centers has ever
  been named. (That challenge from Round 1 §Q7 remains open and standing.)

Budget: one GPU campaign, existing seeders, ≤3 days. After C5, 2B is closed
either way and should stop consuming council cycles.

## 6. Q5 — Stage 3 definition

**Reaffirm positronium-first, no nucleus — with one substrate amendment and
one gate amendment.**

- **Substrate:** primary acceptance object is the **fabric ± pair**
  (ball/anti-ball, §4 rank 1), because its masses are honest ledgers and its
  PASS/FAIL feeds the standing multi-fabric decision gate. The lock-pair
  (Round-1 Q8 protocol) is retained as the toy/development substrate; a
  lock-only PASS is a development result, not Stage-3 acceptance.
- **Gates (Round-1 Q8 list stands, plus):**
  6. **Time honesty:** every orbit gate is stated in measured orbital periods
     (T ≥ 3·T_orb from the *measured local* force), never in absolute revs
     over fixed T. This arc demonstrated why.
  7. **Field honesty:** before release, a static scan verifies the local force
     exponent n<2.5 across the intended radial band in the actual production
     box. No orbit claim, PASS or FAIL, from a band where n≥3.
- All other Round-1 clauses unchanged: ≥5 revs, ±10% seed-velocity
  robustness, closed ledger with radiated energy accounted, negative control,
  3D for acceptance. Annihilation still not required — with the fabric ±
  pair, contact annihilation is now a *possible* outcome; if observed, it is a
  bonus discovery, and the acceptance orbit must simply avoid contact.

## 7. Q6 — Long-range risk ("nuclear Q-ball museum + PIC force demo")

The freeze risk is real but its mechanism is specific: **the program keeps
converting protocol failures into physics conclusions, then pricing new
construction projects off them.** That is how a lab drifts into a museum — the
exhibits are fine, and every hallway ends in a mechanism proposal that exists
to explain an artifact. Countermeasures:

1. **Artifact-first triage rule:** no new mechanism may be proposed in
   response to a FAIL until the FAIL survives (a) a gate-arithmetic audit, (b)
   a box/BC audit, (c) a conservation audit (this arc's E1b baseline flunks
   (c) with −20% energy and −20% ball charge). Make this a standing checklist
   in every FINDINGS template.
2. **One program-level kill gate, dated:** if, by end of the next arc, no
   multi-rev orbit exists in *any* substrate under clean-Coulomb, calibrated
   conditions (locks 3D, fabric ± pair, hybrid), then conclude honestly that
   this kernel's dynamics do not support Coulomb orbits — likely a
   radiation/lattice-dissipation statement. That is itself a major structural
   finding; it would redirect the program toward either capacity-co-field-as-
   physics or a statics-only reframing of the atom goal. A museum is what you
   get by *never* letting the program-level question come to a verdict.
3. **Bank the monist result separately.** Run §3.4's inertia-identity tests
   regardless of orbit outcomes. If they pass, the program owns a real,
   citable claim — "mass measured as field ledger, including fusion defect
   inertia" — that does not depend on atoms ever working. That is the
   insurance against the museum scenario, and it is the *right* insurance
   because it is the monist thesis itself, tested where the codebase can
   actually test it.
4. **Keep the atom-critical path singular.** Everything funded next two weeks
   is on exactly one of: orbit validity (F1–F3), (M) inertia (F5), or a named
   closeout (C5). Nothing else.

## 8. Q7 — Two-week plan (PASS/FAIL gates only)

Ordering by information per day; F1 blocks F2/F3 seeding; F5 and C5 are
independent and parallel.

**F1 — Far-field audit (days 1–2, CPU/GPU).** Static pinned-lock force scan
vs. D around a parked ball in the intended production box (start L≥32, N
scaled), map local exponent n(D).
- **PASS:** a radial band spanning ≥1.5× in r with n∈[1.8,2.5]. → seed F2
  there.
- **FAIL:** no such band in feasible boxes. → the medium/BC cannot present a
  Coulomb exterior at orbit radii; escalate to box/sponge design (this becomes
  the blocking engineering problem — *named as engineering, not physics*).

**F2 — One honest hybrid orbit (days 2–6).** Single cell: lock seeded in the
F1-clean band, v_circ from measured local F, T ≥ 3·T_orb, box per F1.
- **PASS:** ≥2 revs, r band ±25%, ball Q drift <2%, Gauss floor, energy drift
  accounted as radiation. → Stage-4 pilot exists; H1-strong (droplet
  insufficient) is dead; proceed to multi-lock scale-up next arc.
- **FAIL — absorption:** with a clean far field and healthy ball, a
  nucleus-surface scale is *now* a measured need → fund §4 rank 4/2 with a
  defined job.
- **FAIL — instability despite n<2.5 and closed ledgers:** first genuine
  evidence for a dynamical/radiation gap → fund §4 rank 2 (capacity co-field).

**F3 — Fabric ± pair (days 3–8, GPU, parallel).** Conjugate-seed anti-ball
(Q≈−100) + ball (Q≈+315). (a) Static force sign/magnitude vs. D (does
opposite-Q attract through the shared gauge field, and where does the
charge-blind bag tail take over?). (b) Orbit attempt under the F2 protocol in
the F1-clean band.
- **PASS (multi-rev):** Stage-3 substrate settled as same-fabric; multi-fabric
  stays closed; positronium acceptance (§6) is next arc's centerpiece.
- **FAIL — bag capture/merge at distances where Coulomb should dominate:**
  this is the CLAUDE.md-sanctioned evidence that same-fabric bag physics
  blocks atoms → open the multi-fabric kernel *design discussion* (user
  authorization step), with data attached.

**F5 — (M) inertia identity (days 4–7, CPU, parallel).** §3.4 tests 1 and 2:
calibrated push on a parked ball (m_inertial vs. ∫T₀₀/c²), then on c6_light
(fusion-defect inertia).
- **PASS:** agreement ≤5% surviving one resolution step. → bank (M)'s inertia
  leg; monist language in CONCEPT.md may cite a measurement.
- **FAIL:** robust disagreement. → freeze monist claims in program documents;
  open a v76/v78 revision memo before further theory writing.

**C5 — Dinuclear rotor, the last 2B kill (days 8–11, GPU).** Per §5.
- Gates as stated there. Either outcome closes 2B.

**Day 12–14 — Decision memo.** Promote/kill on F2/F3/F5/C5 data; make the
multi-fabric call if F3 failed by capture; re-plan next arc against §7's
program-level kill gate.

Not in plan, deliberately: any new co-field, wall, or bag construction; any
E1b-style option matrix; any further static multi-center scan.

## 9. Q8 — Stop-doing list

1. **Stop seed-polish option matrices for orbits.** The 31-cell E1b matrix
   spent its budget measuring gate arithmetic and box screening. One validated
   cell run long beats 31 short cells, permanently.
2. **Stop stating orbit gates in absolute revolutions over fixed T.** All
   orbit gates are henceforth in measured orbital periods with T ≥ 3·T_orb
   (§6, gate 6).
3. **Stop claiming orbit results from radial bands where the measured force
   exponent is ≥3** — including retroactively: annotate the E1b FAIL in
   FINDINGS as protocol-confounded (§2), so future arcs don't cite it as a
   physics negative.
4. **Stop static multi-center scans** (phase, flavor, D, cold/long-T). Twice-
   confirmed negative. C5 is the sole remaining 2B expenditure.
5. **Stop bag/soft-core knob work on hybrid lights.** B4/B5 measured null;
   they act at small r, the problem was at large r.
6. **Stop drawing theory conclusions from runs with open ledgers.** −20%
   energy / −20% ball charge disqualifies a run from supporting any claim
   about missing physics (standing checklist, §7.1).
7. **Stop letting (M)'s status ride on lock experiments.** Locks cannot test
   it (§3.3); F5 can. Until F5 reports, monist wording in shared documents
   stays conditional.
8. **Don't reopen multi-fabric on frustration.** Reopen only on F3's
   capture-FAIL evidence, per the standing authorization gate.

## 10. Closing note

Round 1 said: prove the droplet insufficient before building its replacement.
This arc's true lesson is one level more basic: **prove the experiment valid
before believing its verdict.** The force fix is real progress — the coupling
that everything downstream needs is now measured and clean. But the arc's two
headline negatives are of different quality: the 2B NULL is a good
measurement, twice-confirmed, and should be closed with dignity (one rotor
kill, then done); the orbit FAIL is a protocol artifact wearing a physics
costume, and the program was one step from funding new mechanisms — and
doubting its own mass ledger — on its say-so. Run F1/F2 before any of that.
And run F5 regardless: the monist thesis deserves better than to be argued
about in memos while the one relation this codebase can actually measure —
that mass, including binding defect, is the inertia of a field ledger — sits
untested a week from an answer.

*— Fable*
