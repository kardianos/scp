# Integer ledger kernel track (`cellfabi`)

Production kernel **`v89/cellfab.c`** is unchanged (FP64).

Alternate binary **`v89/cellfabi`** (`cellfabi.c`):

| `ledger_mode` | Behaviour |
|---|---|
| **0** | FP64 dynamics + **LUT trig** (no integer accounts) |
| **1** | **Shadow:** FP64 truth + parallel int64 accounts (quantized Δ + residual) |
| **2** | Dense+space+flight truth: Es/Em/lem from int; field Ee still FP-derived |
| **3** | **Full** energy ledger: + Ee from int with amplitude rescale to match |

Config keys (cellfabi only):

```
ledger_mode=1
ledger_u=1e-12
```

Trig: `gen_trig_lut.py` → `trig_lut.inc` (quarter-wave cos, N=4096).
Helpers: `ledger.h`. Design: `../INT_LEDGER.md`.

## Build

```bash
python3 v89/int_ledger/gen_trig_lut.py
gcc -O2 -march=native -fopenmp -o v89/cellfabi v89/cellfabi.c -lm
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
