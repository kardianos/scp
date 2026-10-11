# v102 gather step — decision (how the background gathers)

**Written:** 2026-10-11, after the registered runs. Gates: `PREREG.md`.
A measurement in models; it rules nothing. No kernel, law or `v97/fab`
change.

- Tool: `cas/gx.c` (C, OpenMP, LAPACK/BLAS). Runs: `cas/run_t0.sh`,
  `cas/run_grid.sh`, `cas/run_grid_c32.sh` (the 3+2 c-on rerun, §7).
- Tables and bars: `cas/tab.go` (Go) → `go run tab.go` in `cas/`.
- The registered 3+2 supplement (n = 7): `cas/run_supp7.sh`, §8.
- Every candidate was run in every combination: none, a, b, c, ab, ac,
  bc, abc, in both fabrics. That is 16 grid processes plus 8 side-variant
  processes.

## Ending: **N**

| | |
|---|---|
| **G1 gathers** | fails for all 8 rule sets in both fabrics. The largest mean end w_lock in any X1 run is 0.0018 (3+2) and 0.0038 (W), against a bar of 0.5 |
| **G2–G4** | the four c-off sets (none, a, b, ab) pass G2, G3 and G4 in both fabrics. The four c-on sets (c, ac, bc, abc) fail G2 and G3 and pass G4 |
| **ending** | no set gathers (not Gc, not Gl). No c-on set holds (not H). So **N** |

**In one line:** no candidate, alone or combined, gathers the
background from field and space. Without a label the background does
not even hold: it evaporates from its open end as light in 13–22 t.u.
With a label it holds, and a refills its vacancies, but only from the
edge, and not into its ground state.

## 0. Tool gate T0: **PASS**

All 16 gate runs (8 rule sets × 2 fabrics; `cas/out/t0_*.out`) pass all
seven checks:

| gate | result |
|---|---|
| 1 Hermiticity | max asymmetry 0 in every sector |
| 2 Eigenpairs | residual ≤ 2.0e-14 |
| 3 LO against `vx.c` | lock energy less Σ E_f: −0.083316457937 (3+2, n = 4) and −0.100991827294 (W, n = 3), equal to vx.c in every printed digit |
| 4 Inertness of a and b | \|L0⟩ is an eigenstate under none, a, b and ab, with residual ≤ 7.6e-15 |
| 5 Bookkeeping | the quantum and energy ledgers have a maximum deviation of 0 in every grid run (each output's `CHECK` line) |
| 6 Completeness and clusters | ≤ 1.2e-14; cluster rates ≤ 7.0e-10 |
| 7 Broadened window (added by amendment, §7) | closed form against quadrature 1.13e-13; ⟨L0\|Γ\|L0⟩ = 0 under every c-off set |

Every grid output also ends in a `CHECK` line. They give ledger 0,
completeness ≤ 2.0e-14, and cluster rates ≤ 1.2e-9. The one exception
is the W ab side run (`gs_M4_ab`), at 3.8e-7. That check compares a
cluster's rate at each member's own energy with its rate at the
cluster's first energy, and the members differ by ≤ 1e-9. The rate the
run uses is the per-eigenstate sum, which is the first of the two, so
the larger value says only that one cluster sat on a steep part of the
window. No verdict is within orders of magnitude of being affected.

## 1. Bars

| rules | 3+2: G1 G2 G3 G4 | W: G1 G2 G3 G4 | G6 (c on) | gathers in both |
|---|---|---|---|---|
| none | ✘ ✔ ✔ ✔ | ✘ ✔ ✔ ✔ | ✘ | ✘ |
| a | ✘ ✔ ✔ ✔ | ✘ ✔ ✔ ✔ | ✘ | ✘ |
| b | ✘ ✔ ✔ ✔ | ✘ ✔ ✔ ✔ | ✘ | ✘ |
| c | ✘ ✘ ✘ ✔ | ✘ ✘ ✘ ✔ | ✔ | ✘ |
| ab | ✘ ✔ ✔ ✔ | ✘ ✔ ✔ ✔ | ✘ | ✘ |
| ac | ✘ ✘ ✘ ✔ | ✘ ✘ ✘ ✔ | ✔ | ✘ |
| bc | ✘ ✘ ✘ ✔ | ✘ ✘ ✘ ✔ | ✔ | ✘ |
| abc | ✘ ✘ ✘ ✔ | ✘ ✘ ✘ ✔ | ✔ | ✘ |

G3 under c fails for want of any run with w_lock ≥ 0.5 to score.

## 2. X1 — gathering from the empty chain

No run gathers. Each candidate stops at a different place.

| rule set | what space or light does to the empty chain | 3+2, dilute S | W, dilute S |
|---|---|---|---|
| none | nothing can arrive as space. Light passes through: 100 in, 100 out | — | — |
| **a** (c off) | **one** arrival, then nothing more. It lands on the boundary link (k−1, k): one quantum fills Λ's end cell (cap 1) and its partner leaves as light. The interior is unreachable: member links carry no event (R10), and the T-floor refuses to fill a lock cell from empty | 1.005 arrivals, then stuck (t ≈ 2.9e8 at ν = 1e-6) | 1.000 arrivals, then stuck |
| **b** (c off) | **a space-to-light converter.** Space arriving on the light link makes a pair. Outside the lock a pair's quanta hop at t_L, so both leave as light | 100 space in → 199.3 light out | 100 → 199.7 |
| **c** | nothing can arrive as space (an empty chain has no member, so no boundary link). Dense light saturates cells transiently (nsat 0.73 in 3+2, 1.41 in W), and they hop on at t_L | w_lock 0.0018 (light dense) | 0.0037 |
| **ab** | as a, then b converts every later arrival into light | 1 quantum held; 198 light out | 1 held; 198 light out |
| **ac** | **a converter at the wall.** Space arriving on link (0, 1) makes cell 0 a member in the target configuration, so the event is allowed. But cell 1 holds one quantum, below its cap of 2, so it hops at t_L, and once it leaves, cell 0 is no longer a member and follows it | 100 in → 199.4 light out | 100 → 199.1 |
| **bc**, **abc** | converters everywhere | 100 → 198.5 / 198.4 | 100 → 199.7 / 199.2 |

Nothing changes under bright light (SL dense). Every set's mean end
w_lock stays ≤ 0.002 (3+2) and ≤ 0.004 (W).

So no candidate makes field stay where space converted it. In the c-off
sets, a quantum is held only inside the placed Λ, and a cannot reach
Λ's interior. In the c-on sets, a quantum is held only in a saturated
cell, and the edge of anything that starts to gather is unsaturated.

## 3. X2 — holding the placed lock |L0⟩

| rule set | dark survival | S dense: w_lock, F/w | SL dense: w_lock |
|---|---|---|---|
| none, a | 1.000 (both fabrics); F/w 1.000 | 1.000, 1.000 | 3+2 0.0006 / 0.079; W 0.017 / 0.025 |
| b | 1.000; F/w 1.000 | 3+2 0.961, **0.103**; W 0.830, **0.311** | 3+2 0.0006; W 0.010 |
| ab | 1.000; F/w 1.000 | 3+2 0.991, **0.102**; W 0.968, **0.259** | 3+2 0.089; W 0.020 |
| c, ac, bc, abc | **0.000**: every quantum leaves as light (6.0 in 3+2, 4.0 in W) in a mean 22.4 t.u. (3+2) or 13.3 t.u. (W) | 0.000 | ≤ 0.004 |

- **Without a label the lock unzips.** Under c, the lock's light-facing
  end cell holds one quantum (below z = 2), so it is not a member.
  The link out of it then has no member end and hops at t_L. L0's
  residual is t_L = 0.9 exactly. When that quantum leaves, the next
  cell loses its member status in turn. The c-on spectrum has no sharp
  lock-like state: in the lock sector the best has F = 0.074 (3+2) or
  0.176 (W), and it leaks at 0.78 or 0.59. Adding a does not stop it.
  Space that arrives at the receding edge converts into light and
  leaves (ac's S-dense light out exceeds the lock's own count).
- **With a label the lock is dark and sharp.** L0 is an exact
  eigenstate with Γ = 0, so in the dark and in pure space (none, a)
  nothing ever happens to it.
- **b keeps the count in space but scrambles the state.** Under b,
  space converts into pairs on the light link. Their quanta exchange
  with the lock's end quantum across the boundary link, which keeps
  each N_i and changes modes, and then leave. The lock is left in an
  excited internal state: F/w falls to 0.10 (3+2) or 0.26–0.31 (W).
  (This mechanism is inferred, not traced: under b alone nothing can
  replace a lock quantum once it has gone, so the change of state must
  come through the exchange.)
  Formation (space out, always on) occasionally takes the end quantum
  too; ab's a refills it, which is why ab holds more weight than b.
- **Bright light defeats every lock** (SL dense, w_lock ≤ 0.09). That is
  already K3's finding for the kernel lock.

## 4. X3 — light is free

P(both light quanta leave as light) = 1.000 for all 16. No set binds
light into anything: not b's exchange on every link, not R12's core,
and not c.

## 5. X4 — the lock sector's spectrum

- **c off:** L0 is an exact eigenstate (residual 7.6e-15 in 3+2, 3.5e-15
  in W), with zero rate out. It is not the lowest state of its own
  sector. The lowest has w_lock = 0 and leaks at 0.82: it is a state
  with quanta in the light cells. The lock is a dark, metastable
  state, as in every v102 model.
- **c on:** as in §3, no lock-like eigenstate exists.

## 6. X5 — vacancy fate (reported, not scored)

| start | rule set | 3+2: P(refilled), F/w after | W: P(refilled), F/w after |
|---|---|---|---|
| end vacancy, 1-D | a | dilute 1.000, 0.31; dense 1.000, 0.32 | dilute 1.000, 0.10; dense 1.000, 0.00 |
| | ab | dilute 1.000, 0.01; dense 1.000, 0.07 | dilute **0.000** (b's conversion wins); dense 1.000, 0.02 |
| | none, b; every c-on set | 0 (c-on sets evaporate) | 0 |
| bulk vacancy, 1-D | every set | 0 | 0 |
| bulk vacancy, side variant | a | dilute 1.000 (t ≈ 7.2e7), dense 1.000 (71.8 t.u.), F/w 0.095 | dilute 1.000 (t ≈ 3.9e7), dense 1.000 (36.9 t.u.), F/w 0.000 |
| | ab | dilute 1.000, dense 0.995, F/w 0.02–0.04 | dilute 0.33, dense 1.000, F/w 0.00–0.03 |
| | none, b | 0 | 0 |
| | SL-dense supplement | ≤ 0.025 for every set | ≤ 0.083 for every set |

- **a gives the vacancy a finite life in a space environment.** Where a
  boundary link reaches the vacancy, a refills it every time: at the
  1-D end, and in the bulk once side links exist (the kernel's 3-D
  case). In 1-D the bulk vacancy is never refilled, because the T-floor
  may not empty the end cell to move the vacancy to the boundary.
- **The refilled lock is almost never L0** (F/w 0–0.32). Every event in
  the model adds or removes a quantum. In the dark, and in space under
  a, no event can act on a lock-pattern eigenstate: no light cell is
  occupied, and Λ's end cell is at its cap. Under b, space reaches the
  lock only through the light link, and that scrambles it rather than
  relaxing it (§3). So nothing carries an excited lock down its
  internal gap (1.40 J0 in 3+2 and 4.21 J0 in W, in L0's
  sector, from `vx.c`). This is K3's "a closed lock cannot relax", met
  again from the refill side.

## 7. What I expected (written before running), against what happened

| expectation (PREREG §6) | outcome |
|---|---|
| Under c off, G1 fails for every set, because Λ's interior is unreachable | **held**, by exactly that mechanism (§2) |
| b never starts a background from bare fiber, because R12's core lifts a (1,1) pair out of the space band | **outcome held, mechanism wrong.** Once broadened, b's arrivals fire freely (100 of 100). The pairs fail to stay because their quanta travel as light, not because R12 blocks them |
| Under c on, the open lock evaporates from its light end (G2 fails) unless a refills it | **held, and a does not rescue it** (§3) |

## 8. The 3+2 supplement (n = 7, Λ = 0 … 4)

Registered in `PREREG.md` §2, to run only once the main grid had
finished. Labelled, not scored.
- **Counts.** The supplement's counts were not registered. It used the
  main grid's 3+2 dilute count (200), and W's disclosed dense rule
  (12 trajectories at T = 1000), because its largest sector (6405) is
  at W's scale. Script: `cas/run_supp7.sh`. Tables: `go run tab.go s7`.
- **T0.** It passes for all eight rule sets. The lock energy less Σ E_f
  is −0.106222549798, equal to vx.c's LO open chain of 5 cells in
  every digit. Every output's `CHECK` line has ledger 0 and
  completeness ≤ 2.1e-14.

**It repeats the n = 6 picture exactly:**

| | n = 7 |
|---|---|
| bars | identical to n = 6: c-off sets ✘ ✔ ✔ ✔, c-on sets ✘ ✘ ✘ ✔ |
| X1 | nothing gathers: the largest mean end w_lock is 0.0005 |
| X2 under c | the lock unzips in the dark. All 8 quanta leave as light, in a mean 30.6 t.u. (c, ac) or 32.8 t.u. (bc, abc), against 22.4 t.u. at n = 6 |
| X2 under c off | dark survival 1.000. b's S-dense lock keeps its count (0.993), with F/w 0.51 |
| X5 | a refills the end vacancy every time (F/w 0.30 dilute, 0.17 dense). In the side variant, it refills the bulk vacancy every time (F/w 0.32 dilute, 0.25 dense; ab 0.02, 0.17) |
| X3 | light is free (1.000) |

The lock is one cell longer, and only its evaporation time grows. No
verdict moves with size.

## 9. What the gather step says

1. **Nothing gathers.** Field does not stay where space converts it.
   In every rule set, what holds a quantum is its hop amplitude: t_in
   = 0.1 J0 inside a lock, against t_L = 0.9 outside. That amplitude
   comes from a placed label (c off) or from saturation (c on), and
   every growing edge lacks saturation. The binding that might hold a
   quantum without either is the leading-order exchange, g²/Δ_type ≈
   0.005–0.010. That is about 1/400 of the light bandwidth 4t_L = 3.6,
   so nothing self-traps. (This is a statement about the registered
   constants, not a scan.)
2. **The label cannot simply be removed.** Membership by saturation
   (c) makes every open edge travel, so the background unzips into
   light. A literal "no label" reading of the background fails before
   gathering is even asked.
3. **a is repair, not growth.** It needs a member end. It refills what
   formation removes, at the edge in 1-D and in the bulk once side
   links exist. So a particle (a vacancy) has a finite life in a space
   environment. But a cannot grow an empty region, and its refills
   leave the lock in an excited state that nothing relaxes.
4. **b is conversion without retention.** Unlabelled fiber turns each
   arriving space quantum into two light quanta that leave
   (≈ 199 out per 100 in). It is the opposite of gathering: it makes
   light out of space.
5. **The background stays placed.** The open link of field → particle
   named in `../vacancy/DECIDE.md` remains open. This step narrows it:
   the missing piece is something that makes converted field stay put
   without a placed region and without saturation.

## 10. Process notes

1. **Model artifact found mid-run, and the deviation it forced.** The
   first launch used the sharp golden-rule window. The two-cell light
   region has a discrete spectrum (one quantum at E_f(1) ± t_L), so
   every event leaving light there landed at δ = 0.024 ± 0.9, outside
   the ±0.12 space band. That blocked every a-event from the empty Λ
   and made b look R12-blocked. The window was lifetime-broadened, and
   the change was disclosed in `PREREG.md` §6 with a new gate, T0-7,
   before any registered verdict was read. The sharp outputs are kept,
   unread for any verdict, in `cas/out/superseded_sharp/`.
2. **Tool faults caught by the gates and smoke tests,** all before any
   registered run:
   - round-off rates at exact band-edge coincidences (EDGE = 1e-9);
   - a rate cache that carried over between environments;
   - a sign error in the side-a kernel (side a never fired; now the T0
     check is Kf = Ka, which agree to 1.3e-3);
   - a fixed 64-coefficient state, overflowed by a cluster of 82.
3. **Cost.**
   - Broadened tails give every target state a rate. W dense runs
     therefore used 12 trajectories at T = 1000, as disclosed
     beforehand.
   - The 3+2 c-on processes thrashed against the rate-cache cap. The cap
     became a runtime key (`cache=`). A tiny and a large cache gave
     byte-identical output (3+2 ac, exp = hold, 3 trajectories). The
     four 3+2 c-on jobs were then stopped, by exact PID, and rerun
     with cache = 2e8 (`run_grid_c32.sh`). The cache changes speed
     only.
   - Wall time was about 4.1 h for the grid (00:32–04:40). The c-on
     rerun took 32 min (02:18–02:50, from the output files' mtimes).
   - `out/grid.log` has two faults. Its 02:17:59 "done" lines for the
     four 3+2 c-on jobs were written by xargs's wrappers when the
     original processes were killed. The rerun's own "done" lines were
     overwritten: the grid's stdout was redirected without append, so
     it wrote over them at its own offset. The outputs themselves are
     complete, and each ends in its `CHECK` line.
4. **A process-handling fault.** A waiter script used `pgrep -f` on a
   pattern that also matched its own command line. It was found and
   killed by exact PID, with no effect on any run. Since then,
   processes are identified by exact PID through `/proc/PID/cmdline`
   (`ps` shows gx's argv rewritten).
5. **Binaries.** Most of the grid ran on `gx` built at 00:20, before
   the `cache=` key was added. The c-on rerun used `gx_c`, the
   committed source. `gx` was rebuilt from the committed source after
   the grid, and its T0 output is identical to `gx_c`'s.

## 11. Limits

- **Small 1-D chains:** 3+2 n = 6 (lock of 4), W n = 5 (lock of 3), and
  the n = 7 supplement. Two light cells and a Markov lead stand in for
  the world. Every verdict is about these sizes.
- **Leading order only.** The exchange is −(g²/Δ)A†A at R12. Nothing
  here tests all-orders binding.
- **The broadened window is an approximation.** It gives the light
  cells their continuum widths, which the world gives them, but it
  remains a secular golden-rule picture. In particular, gathering that
  needs coherent multi-quantum processes could be missed. Nothing in
  the data hints at such a process.
- **W dense runs have 12 trajectories,** so their standard errors are
  large (W b S-dense w_lock 0.83 ± 0.10). No verdict sits near a bar.
- **c is tested in 1-D, its best case** (z = 2). In 3-D, c's literal
  reading needs six converted links per cell.

## Then

The background still has to be placed. This step names what would have
to change for it to gather: **converted field has to stay put by
itself**, without a placed region and without saturation. Candidates
for a next step, each to be tested in all combinations and each
needing a ruling before it could enter the law:

| | candidate | what it would do |
|---|---|---|
| d | **conversion marks the link** (withdrawn under R13: see the addendum) | a link made by a space conversion stays converted, so its quanta hop at t_in, until formation unconverts it. The label is then written by the event, not placed. This is v98 §1 ("locked field cannot travel as itself") read per link |
| e | **edge membership by contact** | a cell with N ≥ 1 belongs if a neighbour is a saturated member, so the open edge holds by its converted link. G4 (light stays free) is the test it must pass |
| f | **an internal relaxation channel for the lock** | a refilled or scrambled lock has no way down to L0 (§6). Something that keeps every N_i must carry its internal gap away, or the background's state is set by its history |

d with a is the combination that could, in principle, grow a background
link by link from a wall: each arrival converts one more link and sends
its partner out as light.

## Addendum, 2026-10-11: read under R13 (no label)

After this decision the user ruled: "Any encoding by label is an
instant fail intrinsically. If the setup of quanta cannot survive
without label, then it is dead." It is recorded as R13 in
`../README.md`. What it does to this step:

- **The c-off arm (none, a, b, ab) holds only because of a placed
  label, Λ.** Its G2–G4 passes and a's refills (§3, §6) are therefore
  void as evidence that anything survives. They describe a dead
  object.
- **The label-free arm (c, ac, bc, abc) is the whole verdict.** All
  four die: the lock unzips from its open edge into light in
  13–22 t.u., and nothing gathers. In these models, a finite piece of
  background cannot exist without a label.
- **The edge is what dies, not the interior.**
  - Under c the saturated interior keeps t_in. Only the unsaturated
    end cell travels.
  - The whole 3+2 lock's binding is 0.083. The light band's half-width
    is 1.8. So the lock sits inside the light continuum as a
    resonance, not below it as a bound state.
  - Nothing label-free in this model moves it out of the band, and the
    constants may not be scanned to do so.
- **Candidate d, as "Then" words it, is withdrawn.** It stores
  "converted" on the link, which is a label. Only a version that reads
  conversion from the field state itself would be admissible. e and f
  are rules on the field state, so R13 admits them.
- **The vacancy step's ring results** (`../vacancy/DECIDE.md`) may
  survive this. On a closed, saturated ring with isolated vacancies,
  c's membership is every saturated cell, which is the set the label
  named. So those results should carry over without a label, but that
  needs a check, not an assumption. Its open-chain results placed
  their edges, and they do not carry over as they stand.

**What remains possible:** a background with no edge. The medium
itself, converted, closed or unbounded, has nowhere to evaporate from.
No test yet asks whether such a medium survives coupling to light and
space, with a vacancy in it.
