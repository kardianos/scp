# g1 run status — the 3-D geometry gate

**Read `../../FRAME.md` first** (R14: something with gaps, not
something in something). Then `PREREG.md`.

## State: REGISTERED, NOT RUN (2026-10-11)

The user, 2026-10-11: "Set it up, register it, note the current status
of the run so any additional agent will find the status, commit (no
coauthor) and push. Do not yet run."

**Do not start T0 or the grid without the user's explicit go-ahead.**

| item | state |
|---|---|
| `PREREG.md` | written and committed. Registered; do not change it once any gate or run has executed |
| `cas/geo3d.c` (C) | written. Compiles clean with `-Wall -Wextra`. **Never executed** |
| `cas/tab.go` (Go) | written. `go vet` clean. **Never executed** |
| `cas/run_t0.sh`, `cas/run_g1.sh` | written. **Never executed** |
| T0 gates (PREREG §6) | **not run** |
| registered grid (40 runs: 8 sets × L = 6, 8, 10, 12, 16) | **not run** |
| `DECIDE.md` | not written |

## When the user says go

1. From `cas/`, run `./run_t0.sh`. It builds the tool and writes
   `out/t0.out`, whose last line must be `T0 PASS`.
   - **If a gate fails:** fix the tool, not the bars. Record the fix as
     a dated amendment in PREREG §8 before any registered run. Then
     rerun T0.
   - My estimate, unmeasured: 10–30 min. T0.3b runs 8 × 10⁷ proposals
     on 9 quanta.
2. Run `./run_g1.sh` detached. Use the CLAUDE.md pattern
   `setsid ... </dev/null >/dev/null 2>&1 &`, because the grid takes
   hours. It refuses to start unless T0 passed.
   - It runs 12 processes at a time on the 16-core machine, heaviest
     first.
   - Progress: `out/run.log` gets one line per finished run. The marker
     `out/run.done` appears at the end.
   - My estimate, unmeasured: the L = 16 runs with the link rule (b)
     take about 1.5–2 h each, and the whole grid a few hours.
   - Kill only by exact PID, never by `pkill -f` on a pattern.
3. Run `go run tab.go > out/tab.out`.
4. Write `DECIDE.md` (ending G / H / F, PREREG §7). Update this file
   and the status row in `../../README.md`. Commit and push with **no
   Co-Authored-By footer**.

**Update this file at every change of state:** T0 run, grid started
(with its PIDs), grid finished, DECIDE written.

## What this gate decides

Whether a medium of bare quanta and contacts can hold three dimensions
by its own local rules, with no container. The candidates are flat
squares, an octahedral link and no triangles, in all eight
combinations. And, if it holds, what a gap does there: whether it is
readable and persistent, how far its defect reaches, its signed "less
where" and its field footprint, and the same for an excess. Every 3-D
step waits on this gate: a 3-D lane a, c3's 1/r test, and b2.
