# c3 run status — two defects in the medium that holds three dimensions

**Read `../../FRAME.md` first** (R14: something with gaps, not
something in something). Then `PREREG.md`.

## State: REGISTERED, NOT RUN, NOT COMPILED (2026-10-12)

Registered at commit 6f9b110.

The user, 2026-10-12: "Setup the next combintorial step, and commit (no
coauthor) and push, but do NOT run yet. I have another process running,
so you can read and write, but do not run anything, not even a smoke
test yet."

**Nothing has been executed for c3:** no compile, no `go vet`, no gate,
no shakedown, no registered run. Another process of the user's is
using the machine. **Do not start anything until the user says go.**

| item | state |
|---|---|
| `PREREG.md` | written and committed. Registered; do not change it once any gate or run has executed, except by a dated amendment in §8 |
| `cas/c3pair.c` (C) | written and read through. **Never compiled, never executed** |
| `cas/tab.go` (Go) | written. **Never vetted, never executed** |
| `cas/run_t0.sh`, `cas/run_c3.sh` | written. **Never executed** |
| T0 gates (PREREG §6) | **not run** |
| registered grid (96 processes: 8 sets × L = 8, 12, 16 × 4 placements) | **not run** |
| `DECIDE.md` | not written |

## When the user says go

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
