# Integer ledger campaign results (2026-07-29)

**Production kernel** `v89/cellfab.c` (FP64) is unchanged.  
**Alternate** `v89/cellfabi.c` → binary `v89/cellfabi`.

Laws: `battery/laws_V2g.cfg`. Acceptance: full physics CHECKS in
`battery/battery.py` (not process exit codes alone).  
Artifacts: `runs/<variant>/*.log`, `runs/SUMMARY.tsv`,
`runs/PHYSICS_SCORE.tsv`.

MASS waves remained **blocked** during this work (numerics first).

---

## 1. Modes

| mode | tag | behaviour |
|------|-----|-----------|
| — | **fp64** | production `cellfab` |
| 0 | **i_m0** | same dynamics; high-accuracy poly trig (`lut_cos`/`lut_sin`) |
| 1 | **i_m1** | shadow: FP64 truth + parallel int64 accounts (quantized Δ + residual) |
| 2 | **i_m2** | Es/Em/lem integer writeback (stale pre-fix suite; not re-run) |
| 3 | **i_m3** | matter+flight integer truth; field energy **shadowed** in int, unitary FP hops preserved |

Config (cellfabi only): `ledger_mode=N`, `ledger_u=1e-12` (default).

---

## 2. Fixes that landed mid-campaign

### Trig (was the i_m0/i_m1 failure mode)

- **Broken:** linear quarter-wave table N=4096 → phase error accumulated in
  field hops; `rel_drift` up to ~1e−8 (q2) and ~1e−13 on several blob runs;
  failed `DRIFT_MAX = 5e−14`.
- **Fixed:** Cody–Waite-style minimax poly on reduced angle
  (`int_ledger/gen_trig_lut.py` → `trig_lut.inc`). Self-test max\|err\| vs
  `math.cos` ~ **1.5e−15**. Optional coarse table path kept as `lut_cos_table`
  for future experiments only.

### e3b on mode 3

- **Broken:** per-step writeback `Ee ← iEe·u` + amplitude rescale destroyed
  coherent tilt transport: speed **0.00203** (need ≥0.003).
- **Fixed:** do **not** snap field amplitudes each step; `iEe` tracks FP
  field energy only. Es/Em/lem remain integer truth with residual + park.
  After fix: speed **0.00344**, cos **0.880** (both bars green).

---

## 3. Final physics matrix (post-fix full retest)

| variant | physics pass | max \|rel_drift\| (order) | integer sum_err |
|---------|--------------|---------------------------|-----------------|
| **fp64** | **20/20** | ~1e−15 | n/a |
| **i_m0** | **20/20** | ~1e−15 | n/a |
| **i_m1** | **20/20** | ~1e−15 (FP) | nonzero (shadow quantize noise; expected) |
| **i_m3** | **20/20** | ~1e−15 (FP) | **0** on scored runs |
| i_m2 (stale) | 2/20 | ~5e−13 | 0 int park, FP drift fails bar |

**e3b_blob_tilt (post-fix):**

| variant | speed (≥0.003) | cos (≥0.8) | drift |
|---------|----------------|------------|-------|
| fp64 / i_m0 / i_m1 | 0.00665 | 0.95 | ~1e−16 |
| i_m3 | 0.00344 | 0.88 | ~3e−15 |

i_m3 is **slower translation** than FP64 but **above the bar**.

---

## 4. What each mode taught us

1. **Trig accuracy is load-bearing for conservation.** Soft gates + unitary
   hops amplify small cos errors into ledger drift above 5e−14. A “fast
   table” is not free if it fails the ratchet.

2. **Shadow mode (1)** quantifies FP↔int disagreement without changing
   physics: FP bars match i_m0; `max_sum_err` grows (independent residual
   quantization of many accounts — not a claim of exact int conservation).

3. **Full matter ledger (3)** can pass the **entire battery** with
   `sum_err=0` on integer totals, **if** field unitarity is not smashed by
   per-step energy snaps. Exact mass bookkeeping (M-R1) cares about this
   split: pin matter packages in int; keep continuous amplitude where the
   model is continuous.

4. **Mode 2 as first written** (partial writeback without field policy)
   fails the drift bar systematically — recorded as a dead end until redesigned
   (finer `u`, paired residuals, or same field policy as mode 3).

---

## 5. Speed (side measurement, not a campaign goal)

`e1_conserve`, OMP=4, 3 repeats, median wall time (noisy; FP64 alone
ranged ~16–26 s):

| variant | median (s) | vs fp64 |
|---------|------------|---------|
| fp64 | 21.9 | 1.00× |
| i_m0 | 17.7 | ~0.81× |
| i_m1 | 18.1 | ~0.82× |
| i_m3 | 19.3 | ~0.88× |

Integer path **not slower** here; poly trig is likely the win, ledger tax
partially spends it. Treat as order-of-magnitude parity until a full
suite timing table exists.

---

## 6. Recommendation

| Use | Kernel |
|-----|--------|
| Standing science / MASS waves | **`cellfab` FP64** until cellfabi is promoted by explicit decision |
| Numerics / package honesty experiments | **`cellfabi` mode 3** (battery-green) |
| Trig-only A/B | mode 0 |
| Measure quantize noise | mode 1 shadow |
| Mode 2 | do not use until redesigned |

**Promotion criteria (suggested):** byte-stable OpenMP parity check;
mode 3 full battery green on a clean tree (done once); document e3b speed
margin; optional second foam seed on e3b/g1.

---

## 7. Files

| path | role |
|------|------|
| `v89/cellfab.c` | production FP64 (untouched) |
| `v89/cellfabi.c` | integer-ledger variant |
| `v89/int_ledger/gen_trig_lut.py` | regenerate `trig_lut.inc` |
| `v89/int_ledger/ledger.h` | quantize helpers |
| `v89/int_ledger/compare_battery.py` | multi-variant battery runner |
| `v89/INT_LEDGER.md` | design note |
| `v89/int_ledger/runs/` | logs + SUMMARY + PHYSICS_SCORE |

---

## 8. Log

* **2026-07-29** — first full matrix: i_m0/i_m1 15/20 (linear LUT drift);
  i_m3 19/20 (e3b only).  
* **2026-07-29 (later)** — poly trig + no field snap: **fp64 / i_m0 / i_m1 /
  i_m3 all 20/20**. i_m2 left stale. Speed spot-check: integer slightly
  faster on e1. RESULTS written; MASS still unblocked only by user call.
