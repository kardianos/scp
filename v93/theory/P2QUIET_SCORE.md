# P2 quiet dual park (no bath flood) — FINAL T=2000

**Dir:** `p2quiet_T2000_20260808_230025`  
**Setup:** unitary blob, k_rad=0.05, p_rad=4, **bath_xb=0**, rad_dual∈{0,1}  
**Done:** 2026-08-09 ~00:30 CDT

## Full table

| arm | dual | Em path | Em_end | last-half slope | class |
|---|---:|---|---:|---:|---|
| Q_qd12_d0 | 0 | 149→131→121→116→**113** | 113.2 | −0.008 | **SLOW_DRAIN** |
| Q_qd12_d1 | 1 | 149→133→125→123→**134** | 133.6 | **+0.009** | **RECOVER** (not parked) |
| Q_qd36_d0 | 0 | 149→136→129→125→**122** | 122.5 | −0.007 | **SLOW_DRAIN** |
| Q_qd36_d1 | 1 | 149→137→131→128→**129** | 128.6 | −0.002 | **SLOW_DRAIN** slowed |

### Ratios d1/d0

| qd | ratio | ΔEm |
|---|---:|---:|
| 1.2 | **1.18** | +20.4 |
| 3.6 | **1.05** | +6.1 |

## Judgment

1. Quiet dual **does not box-fill** (no bath_xb) — correct apparatus for park test.  
2. Law trap + dual: late **inventory recovery** (slope flips +) — I overcompensates O at low Em, not a fixed point yet.  
3. Deep trap + dual: mild assist only (ratio 1.05); still draining.  
4. **Not adoption-ready** as Em* park; warrants longer T=5000 (launched).  
5. Dual remains default-off experimental.

*Follow-on: `p2quiet_T5000_*`.*
