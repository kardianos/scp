# v102 kernel step, K3 — pre-registration (a lock in light under R12: what it builds, whether it walks, whether a pair holds)

**Written:** 2026-10-09, before any registered K3 run. The scouting runs
used seed 9 only, which is not a registered seed. They are disclosed in
§5.
**Authority:** user, 2026-10-09: "Go with singlet core, keep going"
(`../README.md` R12); `PLAN.md` §2 K3; `K2e_DECIDE.md` §Then.
**No code change.**
- K3 runs the kernel as committed in K2e, under R12: `lockq_p11=2
  lockq_proj=1`.
- The analysis tools gain readouts only: `cas/k2d_replay.py` applies
  seeded jumps, and it reports the profile after the last event and the
  pair's distances; `cas/analyze_k3.py` is new.
- **A replay fix, found while testing `analyze_k3.py` on scout outputs
  (seed 9 and K2e's E2 files, none registered).**
  - With `lockq_every=10` the replay missed the kernel by about 1e-3.
  - An event fires at the end of the kernel step in which its lock step
    runs (`lq_step` follows `step()`), so at the printed t + dt, not at
    t + dt·every.
  - All earlier runs used `lockq_every=1`, where the two agree, so no
    recorded number changes.
  - Fixed, the scouts replay to ≤ 4.4e-11. The same fix is in
    `cas/k2e_conv.py`.
- No law constant changes. No `v97/fab` change.

## 1. Questions

| | question | read from |
|---|---|---|
| Q1 | **What does light build in the lock** under R12, over long times? | part A: events, the final holes, blocks, mean B, bound weight |
| Q2 | **Does a nucleon walk?** | part B1: a single hole seeded mid-chain and beside an end, evolved to T = 2000 |
| Q3 | **Does a pair hold together?** | part B2: two holes seeded apart in complementary modes; the hole–hole distance of the bound part and of the unbound part, at t = 0 and at T |
| Q4 | **Is a size forced** in a given configuration? | part A: the final N_h and the blocks, seed by seed |

## 2. Configurations

**Shared by every run:**
- the V3a law (`k_rad=0.05 p_rad=4 rad_clock=0 wf_on=1 wf_floor=0.01
  wf_far=99`);
- `lockq=1 lockq_p11=2 lockq_proj=1`;
- T = 2000, `lockq_every=10` (dt_L = 0.2), `diag_every=5000` (a lock
  line every 100 t.u.);
- 3+2: `lockq_M=2 lockq_n=10 fmodes=2`; W: `lockq_M=4 lockq_n=7
  fmodes=4`.

**Part A, light** (`exp=pulse`, split light `fseed_mode=-1`):
- the chain placed so the pulse centre falls mid-chain: `lockq_dx=-7.5`
  (3+2) and `-5.0` (W);
- dim light (amp 1, σ 2.2) and bright light (amp 4, σ 3);
- `lockq_seed` = 1 … 5;
- 2 fabrics × 2 lights × 5 seeds = 20 runs.

**Part B, seeded** (`exp=bath`, which has no field; `lockq_rate=0`):

| run | 3+2 n = 10 | W n = 7 |
|---|---|---|
| B1-mid | one hole at site 4, mode 0 | one hole at site 3, mode 0 |
| B1-end | one hole at site 1 (beside the end cell), mode 0 | one hole at site 1, mode 0 |
| B2 | holes at sites 3 and 6 (d = 3), modes 0 and 1 | holes at sites 2 and 4 (d = 2), modes 0 and 3 |

- 6 runs. W n = 7 has five bulk sites, so its pair is seeded at d = 2.
  The W pair's size is about 3 (`../attract/`), so d = 2 overlaps it.

## 3. Bars (every run)

| bar | value |
|---|---|
| ledger | \|drift_rel\| ≤ 1e-12 |
| bookkeeping | `check_violations` = 0, `window_viol` = 0 |
| replay | `cas/k2d_replay.py … score proj` matches the kernel's last lock line (⟨H⟩/J0 and every profile entry) to 1e-8 at T = 2000 |

The projection's most-used sweep count is reported for every run.

## 4. Readouts and answers (registered, not gated)

**Part A, per run:**
- events: their times, sites, modes, δ, and the number of N = 1
  neighbours each had;
- the final N_h and hole sites (nh > 0.5);
- **blocks**: maximal runs of adjacent holes, with their lengths;
- adjacency, mean B, bound weight (the lock's weight on its sector's
  levels with B_k < −0.05 J0);
- the walk: D = ½ Σ_i |nh_i(T) − nh_i(after the last event)|.

**Q1 answer:** each run gets the first class that applies, and each
configuration (fabric × light) is summarized by its count of each class
over the five seeds:

| class | condition |
|---|---|
| **nucleus** | bound weight ≥ 0.10 |
| **block** | adjacency ≥ 0.5 and mean B > 0 |
| **singles** | N_h ≥ 1 |
| **empty** | N_h = 0 |

**Q4 answer:** a configuration "forces its build" if all five seeds end
with the same N_h and the same block lengths.

**Q2 (B1):**
- **walks:** D ≥ 0.5 by T;
- **pinned:** D < 0.1 by T;
- in between: "slow".
- The hole's ⟨x⟩ and the spread of its nh profile at T are reported.

**Q3 (B2):**
- the bound weight of the seeded state (conserved in a closed lock;
  the replay reports it at t = 0 and at T);
- ⟨d⟩ of the bound part and of the unbound part, at t = 0 and at T.
- **The pair holds** if the bound weight is ≥ 0.10 and the bound part's
  ⟨d⟩ changes by less than 15% from t = 0 to T.

## 5. Expectations, stated before the runs

**The scouts (seed 9, T = 20 unless stated):**
- At `lockq_dx=-7.5`:
  - 3+2 bright made 7 holes on sites 2–8, one block;
  - 3+2 dim made 4 holes on sites 4–7, one block;
  - a faint pulse (amp 0.5) made none.
- At `lockq_dx=-5.0`:
  - W bright made holes at sites 1, 3, 4;
  - W dim and faint made none.
- At the end placement (`lockq_dx=0`), T = 200, bright light:
  - 3+2 made 7 holes, the last at t = 6.8;
  - W made 3 holes by t = 2.8 and a fourth at t = 106. The box is
    periodic, so light still reaches the lock long after the pulse.
- Timing: T = 200 took 4.5 min (3+2) and 8 min (W).

**Expected answers** (from `K2e_DECIDE.md`, `K1_DECIDE.md` and the
scouts):
- **Q1:** blocks wherever light lights neighbouring cells. Dim W may be
  empty. No nucleus.
- **Q2:** the mid hole walks; the hole beside the end is pinned.
  - Its energy sits about 3.6 J0 from the bulk band, and the band is
    only 0.11 J0 wide (3+2) or 0.39 J0 (W).
  - At the bulk band's width, one site takes about 700 t.u. (3+2) or
    200 t.u. (W).
- **Q3:** a bound weight near M3a's (0.5–0.66 at d = 3 in 3+2). The
  bound part stays compact, and the unbound part spreads.
- **Q4:** bright light fills the lit stretch, so N_h is set by the light,
  not the lock. Whether it repeats across seeds is open.

## 6. Endings

| ending | condition |
|---|---|
| **Go** | every bar passes; Q1–Q4 are answered as measured |
| **Stop** | a bar fails; the cause is found and stated first |

## 7. Not in K3

- No reverse event (a hole giving its quantum back to the field). The
  lock only fills.
- No sink for binding inside the lock (K1-2).
- No MPS lock: exact sectors cap the lock at 3+2 n ≤ 12 and W n ≤ 7.
- No pricing other than R12.
