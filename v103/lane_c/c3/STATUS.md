# c3 run status — two defects in the medium that holds three dimensions

**Read `../../FRAME.md` first** (R14: something with gaps, not
something in something). Then `PREREG.md`.

## State: T0 PASS, REGISTERED GRID RUNNING (started 2026-10-12 20:46:54)

**Grid processes (kill only by exact PID):** driver `run_c3.sh` PID
2075024, `xargs -P 12` PID 2075031. The c3pair children change as
processes finish; list them with `ps -o pid,etime,cmd -C c3pair`.
Progress is in `cas/out/run.log`, one line per finished process. The
grid is finished when `cas/out/run.done` exists.

Registered at commit 6f9b110.

The user, 2026-10-12: "Setup the next combintorial step, and commit (no
coauthor) and push, but do NOT run yet. I have another process running,
so you can read and write, but do not run anything, not even a smoke
test yet." Then, the same day: "Go."

| item | state |
|---|---|
| `PREREG.md` | written and committed. Registered; do not change it once any gate or run has executed, except by a dated amendment in §8 |
| `cas/c3pair.c` (C) | compiled 20:44 with `-Wall -Wextra`: no warnings, no change needed |
| `cas/tab.go` (Go) | `go vet` clean, no change needed. Not yet executed |
| `cas/run_t0.sh` | ran 20:44–20:46 (2.5 min) |
| `cas/run_c3.sh` | **running**, started 20:46:54 |
| T0 gates (PREREG §6) | **PASS** (`cas/out/t0.out`), every gate in every set, no tool change. Highlights: T0.3 restricted = full quench in 70/70 cases; T0.4 pairs at R = 16 are exactly their singles in all 8 sets × GG/XX/GX; T0.5 healing gives back the flat ⟨d⟩ exactly with no distance changed |
| registered grid (96 processes: 8 sets × L = 8, 12, 16 × 4 placements) | **running** since 20:46:54, 12 at a time, heaviest first |
| `DECIDE.md` | not written |

## Steps (the user said go 2026-10-12; steps 1 and 2 done, step 3 running)

1. **Compile first, without running anything else.**
   - `gcc -O3 -march=native -Wall -Wextra -o c3pair c3pair.c -lm` in
     `cas/`, then `go vet tab.go`.
   - The code has never been compiled. A compile fix that does not
     change behaviour needs no amendment. Any other change is recorded
     as a dated amendment in PREREG §8 before T0.
2. **From `cas/`, run `./run_t0.sh`.** It builds the tool and writes
   `out/t0.out`, whose last line must be `T0 PASS`.
   - **If a gate fails:** fix the tool, not the bars. Record the fix as
     a dated amendment in PREREG §8, then rerun T0.
   - My estimate, unmeasured: a few minutes. The largest gate is T0.4
     at L = 16.
3. **Run `./run_c3.sh` detached.** Use the CLAUDE.md pattern
   `setsid ... </dev/null >/dev/null 2>&1 &`, because the grid takes
   hours. It refuses to start unless T0 passed.
   - It runs 12 processes at a time, heaviest first (L = 16 with the
     link rule).
   - Progress: `out/run.log` gets one line per finished process. The
     marker `out/run.done` appears at the end. Each WINDOW line carries
     the CPU seconds of its medium.
   - My estimate, unmeasured: the heaviest process takes about 1.5 h,
     and the grid 4–6 h. An L = 16 process holds about 100 MB.
   - Kill only by exact PID, never by `pkill -f` on a pattern.
4. **Run `go run tab.go > out/tab.out`.**
5. **Write `DECIDE.md`** (ending U / L / Q / N / M, PREREG §7). Update
   this file and the status row in `../../README.md`. Commit and push
   with **no Co-Authored-By footer**.

**Update this file at every change of state:** compiled, T0 run, grid
started (with its PIDs), grid finished, DECIDE written.

## What this step decides

- **The sign test as an interaction, in three dimensions:** gap with
  gap, excess with excess, and gap with excess.
- **The medium** is g1's: bare quanta and contacts, no container, in all
  eight sets. Only the holding sets (a, ac, b, bc, ab, abc) are scored.
- **Readings:**
  - the pair's energy against separation (contact range is expected,
    because the energy is local);
  - the shared less-where (Δ₂ = d_pair − d_1 − d_2 + d_0);
  - "less space between" the two defects, counted through the medium;
  - whether a gap and an excess annihilate (heal the medium).
- **c3 cannot find a long-range force in a frozen crystal.** That
  question needs mobile defects (g2) or a coherent medium (a2).
