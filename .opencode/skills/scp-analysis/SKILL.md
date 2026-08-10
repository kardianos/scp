---
name: scp-analysis
description: "Use when analyzing SCP fabrun/freecell simulation output — parsing diag rows, RESULT blocks, FCS spatial data, ring winding, correlations. Provides structured column-name-based parsing that avoids awk column-index bugs."
---

# scp-analysis

Use the `scpan` CLI for ALL structured analysis of SCP fabrun/freecell simulation
output. **NEVER use awk/sed/rg one-liners to extract diag or RESULT columns** —
positional `$N` indexing silently breaks when the log's pipe-delimited sections
change (a `$15`-vs-`$20` com_dev/maxstep column-swap bug went undetected for an
entire session). `scpan` parses the HEADER to build a column-name → value map,
so queries are by name and survive format drift.

## Build

```bash
cd /home/d/code/scp/v95 && go build -o scpan ./cmd/scpan
```

Binary: `/home/d/code/scp/v95/scpan`. Stdlib-only (no external deps).

## The two input formats

1. **Log files** (`*.log`, plain text): pipe-delimited diag rows (`t= ...`) plus
   `# RESULT ...` summary lines. Pass directly to `scpan diag` / `scpan result`.
2. **FCS cell dumps** (text output of `fcsdump -mode cells file.fcs`): the binary
   `.fcs` file needs zstd, so `scpan ring|cor|profile` consume the TEXT dump.
   Pipe it in (`-` = stdin) or save to a file:

```bash
./fcsdump -mode cells runs/j1_jm1.fcs | ./scpan ring --m 2 -
./fcsdump -mode cells runs/j1_jm1.fcs > /tmp/x.cells && ./scpan ring --m 2 /tmp/x.cells
```

## Subcommands

### diag — named columns over time

```bash
# specific columns (the com_dev/maxstep bug is impossible here — names, not indices)
./scpan diag --col com_dev,Em_tot --last 10 runs/j1_jm1.log

# time range
./scpan diag --col births,deaths --from 30 --to 50 runs/j1_jm1.log

# summary stats (min/max/mean/median) instead of time series
./scpan diag --stats --col com_dev,Em_tot runs/j1_jm1.log

# default (no --col): every column present in the file
./scpan diag runs/j1_jm1.log
```
Output: TSV with a header row, `t` first. Columns: `com_dev Em_tag ret rms conn
z_tag Em_tot births deaths drift_rel phi_nom ...` plus any `key=value` keys
actually present (e.g. `es_core/far`, `edge_dev_mean`). Flags: `--col --from
--to --last --stats`.

### result — RESULT key=value pairs

```bash
./scpan result runs/j1_jm1.log                 # all keys
./scpan result --key drift_rel,cond,rad runs/j1_jm1.log
```
Handles every RESULT shape: bare key/value (`drift_rel -5e-15`), alternating
(`births 18889 deaths 18550`), and `kind k=v k=v` (`conv rough=.. cond=..`).
Missing keys print a warning to stderr.

### ring — ring winding (port of report/analyze_ring.py)

```bash
./scpan ring --m 2 --cx 8 --cy 8 --frame all  /tmp/j1_jm1.cells
./scpan ring --m 2 --frame 80                  /tmp/j1_jm1.cells   # one frame
```
Per frame outputs `t W R_m Em_min Em_mean keep`, computed on the tagged
(`tag != 0`) cells: W = discrete th2 winding about the centre; R_m =
|⟨e^{i(th2 − m·φ)}⟩| at the seeded m; keep = Em_ring/Em_all. Values match the
Python analyzer to the printed precision. `--frame last` (default) | `all` | `<t>`.

### cor — spatial correlation between two FCS columns

```bash
./scpan cor --colA em --colB ee --x 8 --z 8 --dx 2 --dz 2 --axis y --bw 1 /tmp/x.cells
```
Selects the slab |x−cx|≤dx ∧ |z−cz|≤dz, bins along `--axis` (width `--bw`),
takes per-bin means of colA/colB, and prints the Pearson r across bins:
`t pearson_r n_bins mean_<colA> mean_<colB>`.

### profile — spatial profile of one column along an axis

```bash
./scpan profile --col em --axis y --x 8 --z 8 --dx 1 --dz 1 /tmp/x.cells
```
Same region selection as `cor`; per bin prints `t <axis> <col>_mean n`.

## When to use what

| Task | Command |
|------|---------|
| Read a diag column by name over time | `scpan diag --col NAME` |
| Get final/summary numbers from a run | `scpan result` |
| Min/max/mean/median of a column | `scpan diag --stats --col NAME` |
| Ring topological charge / coherence | `scpan ring --m M` (on fcsdump text) |
| Does column X correlate with Y in space? | `scpan cor --colA X --colB Y` |
| How does a field vary along an axis? | `scpan profile --col X --axis y` |

## Rule

**If you are about to write `awk '/^t=/{print $20}'` on a fabrun/freecell log —
stop. Use `scpan diag --col <name>` instead.** Column indices lie; names don't.
