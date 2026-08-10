# P2 dual-intake suite — FINAL (T=600)

**Dir:** `p2_T600_20260808_204313` · kernel fabrun + `rad_dual`  
**Apparatus:** quiet unitary blob + `k_rad=0.05` `p_rad=4` + vacuum `bath_xb` hold (`bath_em=0.05`)  
**Done:** 2026-08-08 ~21:43 CDT · **Not adopted into laws_V3a**

## Claim under test

With `rad_dual=1`, hoard load should track bath xb (slope ~1).  
Incumbent `rad_dual=0`: inventory independent of xb (slope ~0) under pure O tax.

## Results (last frame)

| arm | dual | xb | Em_end | max | nHot | class |
|---|---:|---:|---:|---:|---:|---|
| P2_ctrl_kr0 | 0 | 0.6* | 6258 | 7.76 | 2616 | box-fill (no rad tax) |
| P2_xb03_d0 | 0 | 0.3 | 4214 | 5.26 | 2538 | box-fill + O tax |
| P2_xb03_d1 | 1 | 0.3 | 4568 | 4.52 | 2595 | box-fill + dual |
| P2_xb06_d0 | 0 | 0.6 | 3884 | 4.96 | 2386 | box-fill + O |
| P2_xb06_d1 | 1 | 0.6 | 5361 | 7.14 | 2614 | box-fill + dual |
| P2_xb10_d0 | 0 | 1.0 | 3229 | 5.61 | 1965 | box-fill + O |
| P2_xb10_d1 | 1 | 1.0 | 5521 | 6.97 | 2622 | box-fill + dual |
| P2_xb15_d0 | 0 | 1.5 | 2229 | 5.80 | 1316 | box-fill + strong bath conflict |
| P2_xb15_d1 | 1 | 1.5 | 6098 | 8.67 | 2610 | box-fill + dual |

\*ctrl uses bath_xb=0.6, k_rad=0.

### Dual ratios (Em_d1 / Em_d0)

| xb | ratio | ΔEm |
|---:|---:|---:|
| 0.3 | **1.08** | +354 |
| 0.6 | **1.38** | +1477 |
| 1.0 | **1.71** | +2293 |
| 1.5 | **2.74** | +3869 |

**Direction correct:** dual intake grows inventory vs dual-off at same bath; effect strengthens with xb.

## Honest class judgment

1. **All arms are box-fill / bath-filled**, not compact HOLD (nHot≳1300–2600, Em≫148).  
2. **No quiet Em\* particle park** demonstrated — apparatus holds ambient field on vacuum cells and floods the box.  
3. Dual is a real **I channel** (ledger-side), not yet a monodispersity fix for the S00 clump.  
4. d0 Em **falls** as xb rises (stronger field bath changes partition / thruput); d1 Em **rises** with xb — qualitative dual signature.  
5. Local x\* = Em_cell/cap tracking of xb needs **per-cell / core** meters (not Em_tot) for the slope-1 claim; total Em is the wrong observable here.

## Battery gate

Prior battery RED 6/87: C≡Go abx failed largely because Go printed `| Em_tot=` and C did not.  
**Fix applied 22:00 CDT:** C `freecell.c` now mirrors Go Em_tot diag. Short abx_b2 masked body **PASS**. Full battery re-run still needed before any adoption.

## Recommendation

- Keep `rad_dual=0` default. **Do not adopt.**  
- Next P2 experiment: dual on **quiet no-bath** S00+rad (I balances O on the clump) — that is the Em\* park test mono failed.  
- Or measure core-mean x vs xb with dual, not box Em_tot.
