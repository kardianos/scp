# Integer ledger kernel track

**DEFAULT development kernel (2026-07-29):** `v89/cellfab.c` with
`ledger_mode=3` (full integer matter ledger). Battery and prestress build
`cellfab` from that file.

**FP64 reference:** `v89/cellfabf.c` → binary `cellfabf` (A/B only).

| `ledger_mode` | Behaviour |
|---|---|
| **0** | FP dynamics + poly trig (no integer accounts) |
| **1** | **Shadow:** FP truth + parallel int64 accounts |
| **2** | dense+space+flight writeback (not battery-green as first written) |
| **3** | **DEFAULT** — Es/Em/lem integer truth; field unitary FP + iEe shadow |

```
# optional override (default is already 3 in cellfab.c)
ledger_mode=3
ledger_u=1e-12
```

Trig: `gen_trig_lut.py` → `trig_lut.inc` (minimax poly cos/sin).
Helpers: `ledger.h`. Design: `../INT_LEDGER.md`.

## Build

```bash
python3 v89/int_ledger/gen_trig_lut.py   # if trig_lut.inc missing
gcc -O2 -march=native -fopenmp -o v89/cellfab  v89/cellfab.c  -lm   # DEFAULT
gcc -O2 -march=native -fopenmp -o v89/cellfabf v89/cellfabf.c -lm  # FP64 ref
```

## Compare to FP64 battery

```bash
python3 v89/int_ledger/compare_battery.py \
  --laws v89/battery/laws_V2g.cfg \
  --modes 0,1,2,3 --jobs 2 --threads 4
```

Logs: `v89/int_ledger/runs/<fp64|i_m0|i_m1|i_m2|i_m3>/`
Summary: `v89/int_ledger/runs/SUMMARY.tsv`
Narrative: `RESULTS.md`

## Results

See **`RESULTS.md`** — post-fix battery: fp64 / i_m0 / i_m1 / i_m3 all
**20/20** physics. Mode 2 abandoned (stale). Speed: integer slightly
faster on e1 spot check (~0.8–0.9× wall of fp64).

## MASS waves

Were **blocked** during numerics (user directive 2026-07-29). Numerics
gate now green on mode 3; fleet launch is a separate user call.
