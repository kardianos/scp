# v102 kernel step, K3 — decision (a lock in light under R12)

**Written:** 2026-10-09, after the registered runs (`cas/out/k3_*`).
Gates: `K3_PREREG.md`. No code change: the kernel is K2e's commit, under
`lockq_p11=2 lockq_proj=1`.
- Analysis: `cas/analyze_k3.py` (`cas/out/k3_analysis.out`).
- Runs: `cas/k3_run.sh` (six at a time), then `cas/k3_run2.sh` (the same
  commands for the runs not yet started, thirteen at a time). The
  binary is the K2e commit's kernel.
- Checks added after the runs (reported, not gated):
  - `cas/k2e_conv.py` on every run that reached the projection cap
    (`cas/out/k3_conv_*`);
  - `cas/k3_lastsite.py`, the window and R12 weight of every possible
    next jump in a final state (`cas/out/k3_last_*`, `k3s_last_*`);
  - `cas/k3_ker8.py`, the R12-allowed states with every bulk site holed
    (`cas/out/k3_ker8_*`).
- A seeding defect was found in the kernel's event draw (§5). Part A
  was run again as a labelled supplement with well-mixed seeds
  (`cas/out/k3s_*`, with the same checks).

## Ending: **Go**

- Every bar passes in all 26 registered runs.
- Q1–Q4 are answered as measured, below.
- **In one line:** light fills the lock until no jump fits the window
  (bright light) or as far as the light reaches (dim light). What it
  builds is an unbound block, whose size is set by the lock's length,
  the order in which holes formed, and the light, not by binding.
- A pair that is seeded bound stays bound and compact.
- A single hole walks in the 3+2 bulk and is pinned beside an end.

## 1. Bars

| bar | value | result |
|---|---|---|
| ledger | \|drift_rel\| ≤ 1e-12 | ≤ 3.5e-13 (part A); 0 (part B) |
| bookkeeping | `check_violations` = 0, `window_viol` = 0 | 0 and 0 in every run |
| replay | ⟨H⟩/J0 and every profile entry to 1e-8 at T = 2000 | ≤ 5.0e-11 |

**The projection's sweeps.**
- The most sweeps used:
  - 3+2 bright: 201 in every seed;
  - 3+2 dim: 119 (seed 1) and 201 (seeds 2–5);
  - W bright: 21–24;
  - W dim: 14;
  - part B: 1–2.
- 201 means that at least one candidate jump stopped at the cap.
  `cas/k2e_conv.py` replays each such run and, at every event, continues
  the capped projections to deep convergence:
  - the rate weights change by at most 1.0e-15;
  - the residual singlet weight is ≤ 3.5e-27;
  - the fired states move by at most 3.8e-14.
- In the dim runs the cap was reached between events, not at them:
  at the events the most sweeps were 158–161, with nothing capped.
- As in K2e §E4, the cap costs nothing measurable.

## 2. Q1: what light builds (part A, registered seeds)

| configuration | N_h per seed | blocks (nh > 0.5) | mean B (J0) | bound weight | class |
|---|---|---|---|---|---|
| 3+2 bright | 7, 7, 7, 7, 7 | [8], [7], [7], [8], [8] | +20.7 to +22.0 | 0.000 | block × 5 |
| 3+2 dim | 4, 6, 6, 6, 6 | [2, 2], [6], [6], [6], [6] | +10.1 to +17.2 | 0.000 | block × 5 |
| W bright | 4, 3, 3, 3, 3 | [3], [3], [3], [2], [3] | +5.5 to +16.6 | 0.000 | block × 5 |
| W dim | 1, 1, 1, 1, 1 | none | +4.30 | — | singles × 5 |

- **No nucleus anywhere.** The bound weight is 0.000 in every run with
  two or more holes.
- **Every lit configuration builds blocks** of adjacent holes, far from
  bound, except W dim, which makes one hole.
- **The box keeps relighting the lock.** Bright light makes most of its
  holes in the first burst (t < 3), but holes keep forming long after:
  - in 3+2 bright the seventh hole forms at t = 5.6 to 141;
  - in 3+2 dim holes form up to t = 1576;
  - in W dim the one hole forms at t = 895.6.
- In 3+2 bright the seven holes sit in eight bulk sites. The vacancy
  either stays beside an end cell (seeds 2, 3; block [7]) or spreads
  over the bulk (seeds 1, 4, 5). Then every site has nh ≈ 7/8 and the
  nh > 0.5 readout counts a block of 8.

**Why the lock stops** (`cas/k3_lastsite.py`, which reproduces the
kernel's printed δ of the last event to every digit in all four checked
runs):
- **3+2 bright is full at seven holes.**
  - R12 allows exactly one state with all eight bulk sites holed in
    each of the two charge sectors a jump can reach (`cas/k3_ker8.py`,
    seed 2's final state).
  - From the seven-hole state every next jump reaches it with weight
    ≤ 1e-21 (zero to roundoff).
  - The window refuses that jump as well: δ = −0.128 (seed 1) and
    −0.141 (seed 2), outside the ±0.12 band.
- **W bright at three holes ([2, 3, 4], seed 2) is full.** Every next
  jump keeps an R12 weight (up to about 1e-2), but its δ is −0.136 to
  −0.158, outside the band. Seed 1 formed its first holes beside the
  end cells and reached four.
- **3+2 dim is not full.** Seed 2's next jumps are inside the window
  (δ −0.094 to −0.110), with weights up to 2e-2. It is limited by the
  light.

## 3. Q2: does a nucleon walk? (part B1)

| run | start | D at T = 2000 | profile at T | verdict |
|---|---|---|---|---|
| 3+2 mid | site 4 | 0.627 | spread 0.97 around ⟨x⟩ = 3.99 | **walks** |
| 3+2 end | site 1 | 0.094 | 0.91 still on site 1 | **pinned** (just under 0.1) |
| W mid | site 3 | 0.090 | 0.91 still on site 3 | **pinned** |
| W end | site 1 | 0.134 | 0.87 still on site 1 | **slow** |

- 3+2 behaves as expected: the hole walks in the bulk and is pinned
  beside an end.
- **W did not.** Its mid hole is pinned. W n = 7 has five bulk sites,
  and its middle site is three cells from either end. The expectation
  (about 200 t.u. per site) assumed a bulk band, which a lock this
  short does not have.

## 4. Q3: does a pair hold together? (part B2)

| run | seeded | bound weight, t = 0 → T | ⟨d⟩ of the bound part | ⟨d⟩ of the unbound part | verdict |
|---|---|---|---|---|---|
| 3+2 | sites 3, 6 (d = 3), modes 0, 1 | 0.6419 → 0.6419 | 2.990 → 3.007 (+0.6%) | 2.982 → 3.106 | **holds** |
| W | sites 2, 4 (d = 2), modes 0, 3 | 0.0000 → 0.0000 | — | 2.000 → 1.928 | does not hold (registered rule) |

- **The 3+2 pair holds** for 2000 t.u.:
  - its bound weight, 0.64, is conserved (a closed lock), close to
    M3a's 0.66 on the ring;
  - the bound part stays compact at ⟨d⟩ = 3.0.
- **The W pair could not be tested.** In W n = 7 the lowest level of
  the two-hole (0, 0) sector has B = +0.175 J0, so no level is bound at
  all. The seeded state's bound weight is zero, and its holes stay
  where they were put (D = 0.075). The registered verdict is "does not
  hold", but what it shows is that the lock is too short for a W pair.

## 5. The seeding defect, and the supplement

**The defect.**
- The kernel's event draw is xorshift64, started from `lockq_seed`
  itself (`lockq.c`, `lq_rng = P.lockq_seed`), with no warm-up.
- From a small seed s, its first output is about 5.9e-11·s and its
  second is s/16. Draws from the third on are well mixed.
- So, with the registered seeds 1–5:
  - each run's first event fires at the first lock step with a positive
    rate;
  - its candidate is taken from the lowest 6–31% of the cumulative
    rate, in seed order.
- The bars are untouched, because the replay follows the kernel's own
  events.
- **The readouts are affected.**
  - W dim's five runs are one run: the same event at the same moment.
    Its "forces its build: True" says nothing.
  - In the other configurations, the first event of each seed is not an
    independent draw.
  - K2a–K2e also used lock seeds 1–3, so the same holds for their first
    events. Their gates were per event (ledger, window, sector, A/B) and
    do not depend on it. Their seed-to-seed readouts carry the same
    caveat. They are not re-run.
- **No code change is needed.** `lockq_seed` is read as a 64-bit
  integer, so a well-mixed seed avoids the defect.

**The supplement** (labelled, not pre-registered). Part A ran again with
the same commands and seeds splitmix64(1 … 5): 16294208416658607535,
7960286522194355700, 487617019471545679, 17909611376780542444,
1961750202426094747.

| configuration | N_h per seed | blocks (nh > 0.5) | mean B (J0) | bound weight |
|---|---|---|---|---|
| 3+2 bright | 7, 7, 7, 7, 7 | [3, 4], [6], [7], [7], [7] | +20.8 to +22.3 | 0.000 |
| 3+2 dim | 7, 6, 4, 5, 6 | [7], [6], [2, 2], [4], [6] | +7.1 to +21.0 | 0.000 |
| W bright | 4, 3, 3, 4, 3 | [4], [2], [2], [4], [3] | +5.5 to +16.2 | 0.000 |
| W dim | 1, 1, 1, 1, 1 | none | +4.30 | — |

- Every bar passes (`cas/out/k3s_analysis.out`): drift ≤ 3.4e-13,
  checks and window 0, replay ≤ 6.3e-11.
  The projection check on the eight capped runs (`cas/out/k3s_conv_*`):
  rate weights change by at most 6.3e-16, fired states move by at most
  5.8e-14.
- **The registered picture stands.** There is no nucleus, and every
  lit configuration builds unbound blocks except W dim.
- **The first events now differ.** Outside W dim they fire at t = 0.18
  to 0.58, at different sites and in different modes.
- **W dim is now a real result.** All five seeds make one hole, at site
  1 beside the end cell, at different moments (t = 981 to 1163) and in
  three different modes. The registered runs' single moment (t = 895.6)
  was the defect's.
- **3+2 bright again ends at N_h = 7** in every seed. Only where the one
  vacancy sits differs, which is what the block readout shows.
- W bright ends at 3 or 4 holes, as registered (registered 1 of 5 at
  4; supplement 2 of 5).

## 6. Q4: is a size forced?

**Under the registered rule** (all five seeds give the same N_h and
the same blocks), only W dim forces its build:
- in the registered runs this was the seeding defect;
- in the supplement it holds with independent seeds.

**Read across both sets (ten runs per configuration):**

| configuration | N_h | what sets it |
|---|---|---|
| 3+2 bright | **7 in all ten** | the lock is full. The eighth hole is refused by R12 (weight ≤ 1e-21) and by the window (§2; checked in registered seeds 1 and 2). |
| W bright | 3 (seven runs) or 4 (three runs) | the lock is full either way: every next jump is outside the window, at three holes ([2, 3, 4] and [1, 3, 4]: δ −0.131 to −0.161) and at four ([2, 3, 4, 5]: δ −0.132). Which one it ends at depends on the order the holes formed. |
| W dim | **1 in all ten**, at site 1 | the light. The lock is not full: every next jump is inside the window (δ −0.036 to −0.080), with R12 weight up to 0.5. The light reached the lock once in 2000 t.u., each time at site 1. |
| 3+2 dim | 4 to 7 | the light. Registered seed 2, at six holes, was still filling at T (§2); the supplement's seed 1 reached seven, where the bright runs stop. |

- **So a size is forced in two configurations, but not by binding.**
  - In 3+2 bright the lock is filled up to its length; that is
    jamming.
  - In W dim one hole forms beside the end, and no second one in
    2000 t.u.
  - The 3+2 block is far from bound (bound weight 0.000, mean B ≈
    +21 J0); the W dim hole is a single hole. Neither size belongs to
    an object: the first is the lock's bulk less one site, the second
    is where the light lands.
- **The size the energetics force (pairs, N\* = 2; `../bind/DECIDE.md`)
  never appears in light.** A pair forms only if two holes land at
  d ≈ 3, in complementary modes, with nothing else nearby. Bright light
  fills the gap first, and dim light makes one hole at a time, at
  random sites, with holes that walk.
- When a pair is seeded bound, it holds (§4).

## 7. Limits

- The exact sectors cap the lock at 3+2 n ≤ 12 and W n ≤ 7. In W n = 7
  a single hole has no bulk and a pair has no bound level, so W's Q2 and
  Q3 are about a short lock, not a long one.
- There is no reverse event: a hole cannot give its quantum back, so
  the lock only fills.
- There is no sink for binding inside the lock. A closed lock keeps its
  eigenstate weights, so a block formed unbound stays unbound.
- The box is periodic (L = 16), so a pulse relights the lock for the
  whole run. The late events depend on that.

## Then

K3 closes the kernel step's plan (`PLAN.md` §2), apart from K2c, the
MPS lock. What it shows:
- formation sees adjacency (about 3.6 J0 per N = 1 neighbour) but not
  binding (0.13–0.50 J0);
- a closed lock cannot relax;
- so light builds what lands, and fills until the window closes.

The energetic size cannot show up in light unless something lets a lock
reach its bound states. **Each candidate needs a ruling:**

| candidate | what it would do | needs |
|---|---|---|
| **a reverse boundary event** | a hole gives its quantum back through the boundary link, the time-reverse of formation. Unbound blocks could shed holes, and a lock in light would no longer only fill. | a rule for space coming in through the boundary link (R10 says only where space leaves). New kernel behaviour, off by default. |
| **core (A1) in place of R12** | no two N = 1 cells may share a link at all. Light could not build blocks. Holes would have to land at d ≥ 2, where a complementary pair at d = 3 is born bound (0.66 / 0.52). | replacing R12 with a stricter rule |
| **a sink inside the lock** | the binding goes into a lock quantum or a third hole | none is open (`../bind/` M3b, M3c) |

- **K2c (the MPS lock)** is authorized, but it does not address this.
  It would let W's Q2 and Q3 be asked in a lock long enough to have a
  bulk and a bound pair.
- **Seeds.** Every future lock run uses well-mixed seeds (splitmix64 of
  the run index, as in §5). A kernel fix that mixes the seed itself
  would be new behaviour; it is not made here.
