# v102 kernel step, K2d — decision (the (1,1) pricing in the kernel; a lock in light)

**Written:** 2026-10-09, after the runs in `cas/out/k2d_*`. Gates:
`K2d_PREREG.md`. Code:
- `v97/kernel/lockq.c`: `lockq_p11`, `lq_p11c`, the (1,1) block of
  `lq_two_build`;
- `v97/kernel/freecell.c`: the config field;
- the tools `cas/k2d_replay.py` and `cas/analyze_k2d.py`.

## Ending: **Go**

- All five pricings of `../bind/` now run in the kernel's lock as
  options.
- At the default the kernel is byte-identical.
- Under every pricing it matches `lo.py`, and every lock history
  replays exactly.

The lock-in-light runs give a clear physics split between the pricings
(§D3).

## D1. Inert at default: **PASS**

- **The full battery, both kernels** (`cas/out/k2d_battery.out`):
  **ALL GREEN (87 bars)**, bar table identical to K2b's.
  - 70 of 71 run logs are byte-identical.
  - `p2lc` (the Go kernel) differs only in its wall-clock timings
    (`ms_sync`, `ms_async`), as in K2a and K2b.
- **Kernel stdout against the K2b binary** (`cas/out/k2d_d1_stdout.out`):
  byte-identical for
  - `exp=bath`, `blob`, `slit`, `pulse` (T = 20);
  - K2b's two B4 lock configurations with `lockq=1`;
  - the same 3+2 run with `lockq_p11=0` set explicitly.
- The base binary, rebuilt from HEAD, is byte-identical to K2b's
  final binary.
- There are no new compiler warnings (9 before, 9 after).

## D2. A/B, every pricing: **PASS**

The setup:
- a seeded jump at bulk site 1, next to the end cell, so the lock has a
  (1,1) end link from the start;
- 2000 steps, tL = 40;
- against `cas/k2a_ab.py` with `P11` set
  (`cas/out/k2d_d2_*`).

| pricing | 3+2 n = 8: ⟨H⟩/J0 | max \|C − py\| | W n = 6: ⟨H⟩/J0 | max \|C − py\| |
|---|---:|---:|---:|---:|
| core | +89.0354999485 | 2.7e-11 | +89.0396547225 | 2.7e-11 |
| score | +62.0898841337 | 1.3e-11 | +49.3755523828 | 3.2e-11 |
| upper | −10.8946624999 | 4.5e-11 | −9.8166566115 | 3.0e-11 |
| lower | −21.5516126237 | 3.0e-12 | −22.3116533799 | 5.1e-11 |

The bar is 1e-6. Agreement is at the printed precision (10 digits).

## D3. A lock in light: **PASS** (bars), and the readouts

The setup:
- 60 runs: 3+2 n = 10 (eight bulk cells) and W n = 7 (five bulk cells);
- bright light (amp 4, σ 3) and dim light (amp 1, σ 2.2), split over
  the modes;
- three seeds × five pricings;
- outputs and analysis in `cas/out/k2d_d3_*`.

**The bars, every run:**
- |drift_rel| ≤ 3.6e-15;
- `check_violations` = 0 and `window_viol` = 0;
- the replay matches the kernel's final lock line to ≤ 5.0e-11 (bar
  1e-8).

**Readouts** (bright light; dim light made one hole in 3+2 and none in
W within T = 20, under every pricing):

| fabric | pricing | final N_h (seeds) | events beside a hole | adjacency at the end | bound weight |
|---|---|---|---:|---:|---|
| 3+2 | zero | 8, 7, 8 (of 8) | 14 | 6.3 | 0, 0, 0 |
| 3+2 | **core** | **3, 3, 3** | **0** | **0** | 0.27, 0, 0.27 |
| 3+2 | **score** | **3, 3, 3** | **0** | **0** | **0.51, 0, 0.51** |
| 3+2 | upper | 7, 6, 7 | 12 | 4.7 | 0, 0, 0 |
| 3+2 | lower | 8, 8, 8 | 16 | 7.0 | 0, 0, 0 |
| W | zero | 3, 3, 3 (of 5) | 6 | 2.0 | 0 |
| W | core | 1, 1, 1 | 0 | 0 | — |
| W | score | 1, 1, 1 | 0 | 0 | — |
| W | upper | 2, 2, 2 | 3 | 1.0 | 0 |
| W | lower | 3, 3, 3 | 6 | 2.0 | 0 |

**What this shows:**

1. **Zero, upper and lower let bright light fill the lock with
   side-by-side holes, and nothing binds.**
   - Upper prices only the singlet component of an adjacent pair. The
     other channels cost nothing, so formation next to a hole stays
     inside the window.
   - The end state is a lock of N = 1 cells, far above free holes:
     B = +21 to +45 J0 under zero and upper, +6 to +8 J0 under lower
     (3+2). It is an expensive configuration, not a nucleus.
2. **Core and score keep every hole apart.**
   - No event in any seed formed beside a hole or beside an end.
   - Every bright 3+2 lock stopped at exactly **three holes, at sites
     2, 4 and 7**. Every remaining site was next to a hole or an end,
     so every further event fell outside the window: the lock jams.
   - When one of the three holes had the complementary mode (seeds 1
     and 3: modes 0, 1, 0), the pair at sites 4 and 7 (d = 3) is
     **bound with weight 0.51** (score) or 0.27 (core).
   - Seed 2 made all three holes in mode 0, and no pair is possible.
3. **W n = 7 is too short to show a pair under core or score.** Sites 1
   and 5 are next to the ends. The only d = 3 placements, (1, 4) and
   (2, 5), each use one of them. Bright light made one hole in 20 time
   units.

## Limits stated

- The window uses the jump state's **mean** energy change (K2a's
  rule). Score blocks adjacent formation through that mean: the
  singlet part of the jump state costs about 104 J0. Upper does not:
  its singlet part costs only 5.9–8.8 J0.
  - An energy-resolved event rule would let the non-singlet part of an
    adjacent jump through under score too. Only core forbids adjacency
    outright.
- **Replay timing (a tool fix, made before scoring).** An event fires
  after its lock step's evolution, at t + dt_L, while its line prints
  the step's start t. The first version of `k2d_replay.py` applied
  jumps one step early. That gave 1e-5 disagreements in runs where the
  lock evolved after an event. Corrected, every run agrees to ≤ 5e-11.
- Exact sectors cap the lock at 3+2 n ≤ 12 and W n ≤ 7.
