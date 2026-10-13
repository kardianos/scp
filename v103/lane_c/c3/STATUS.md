# c3 run status — two defects in the medium that holds three dimensions

**Read `../../FRAME.md` first** (R14: something with gaps, not
something in something). Then `PREREG.md`.

## State: DECIDED — ending U (2026-10-13). See `DECIDE.md`. Nothing of c3 is running

The registered grid ran 2026-10-12 20:46:54 to 2026-10-13 00:13:12
(3 h 26 min): 96/96 processes DONE, no FAIL, empty .err files, all 2552
CHECK lines passed. The driver (PID 2075024) and xargs (PID 2075031)
have exited.

**Ending U in all six scored sets** (a, b, ab, ac, bc, abc):
- GG and XX attract at contact;
- GX annihilates at contact (72/72, at the making; the flat medium and
  its where return exactly);
- the energy is exactly local (R_E 3–4; zero at R ≥ 6);
- every pair makes less where, sub-additively, with an overlap that
  reaches past the energy (post-hoc: about ½ per doubling of R).

Grid history:

- **21:18: lowered to nice 10** (renice by exact PID: the driver,
  xargs and the 12 running c3pair; later children inherit it). At
  21:16 another workload of the user's started (zig, testkit, corpus,
  acceptance processes), and the load reached 29 on 16 cores. The
  results do not depend on timing: each process is single-threaded
  and deterministic.
- **Measured pace:** 500 medium sweeps took 3 min (about 0.36 s per
  sweep at L = 16 with b). A window takes 140–160 CPU s. With the
  load, a heavy L = 16 process finished about 7 of its 32 media in
  30 min, so expect about 2–2.5 h each and a longer grid than
  registered (perhaps 7–9 h). A slower grid changes no result. In the
  event it took 3 h 26 min: the L = 12 and 8 processes are far cheaper.
- **tab.go, 21:25:** a tested-on-partial-output fix. §5's "at R = 8"
  column said "contact" when the R = 8 readings were still missing;
  it now says "(missing)". This changes no bar, and the whole-grid
  table is unaffected once every R = 8 reading exists.

Registered at commit 6f9b110.

The user, 2026-10-12: "Setup the next combintorial step, and commit (no
coauthor) and push, but do NOT run yet. I have another process running,
so you can read and write, but do not run anything, not even a smoke
test yet." Then, the same day: "Go."

| item | state |
|---|---|
| `PREREG.md` | written and committed. Registered; do not change it once any gate or run has executed, except by a dated amendment in §8 |
| `cas/c3pair.c` (C) | compiled 20:44 with `-Wall -Wextra`: no warnings, no change needed |
| `cas/tab.go` (Go) | `go vet` clean; one presentation fix during the grid (below); output `cas/out/tab.out` |
| `cas/posthoc.go` (Go) | post-hoc readings (persistence breakdown, relative where-interaction against R); output `cas/out/posthoc.out` |
| `cas/run_t0.sh` | ran 20:44–20:46 (2.5 min) |
| `cas/run_c3.sh` | ran 2026-10-12 20:46:54 to 2026-10-13 00:13:12 |
| T0 gates (PREREG §6) | **PASS** (`cas/out/t0.out`), every gate in every set, no tool change. Highlights: T0.3 restricted = full quench in 70/70 cases; T0.4 pairs at R = 16 are exactly their singles in all 8 sets × GG/XX/GX; T0.5 healing gives back the flat ⟨d⟩ exactly with no distance changed |
| registered grid (96 processes: 8 sets × L = 8, 12, 16 × 4 placements) | **done**, 96/96 |
| `DECIDE.md` | **written 2026-10-13: ending U** |

## Steps (the user said go 2026-10-12; all done 2026-10-13)

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
