# Monodispersity scan — INTERIM (4/8 done)

**Dir:** `mono_T2000_20260808_190046`  
**Setup:** quiet unitary blob, amp_nat=2, p_rad=4, T=2000, seed=20260802  
**Scored:** 2026-08-08 ~20:00 CDT

## Question

Does k_rad × q_detune produce a **monodisperse Em\* fixed point**, or only SLOW_DRAIN?

## Law trap q_detune=1.2 (DONE)

| arm | k_rad | Em | max | nHot | slope | class |
|---|---:|---:|---:|---:|---:|---|
| M_qd12_kr00 | 0 | **148.30** | 2.08 | 53 | ~0 | **HOLD** (S00 clone) |
| M_qd12_kr02 | 0.02 | 124.81 | 1.52 | 31 | −0.0069 | **SLOW_DRAIN** |
| M_qd12_kr05 | 0.05 | 113.15 | 1.31 | 26 | −0.0079 | **SLOW_DRAIN** (≡S05 class) |
| M_qd12_kr10 | 0.10 | 105.36 | 0.84 | 21 | −0.0071 | **SLOW_DRAIN** + peak weaken |

**No Em\* park.** Tax is monotone in k_rad: more k → lower Em, lower nHot, lower Em_max. No intermediate fixed point where O balances an internal I (there is no open I on quiet blob).

## Deep trap q_detune=3.6 (in flight ~t=900/2000)

| arm | Em@900 | note |
|---|---:|---|
| M_qd36_kr00 | 148.0 | HOLD-like |
| M_qd36_kr02 | 137.9 | drain slower than qd12@kr02 |
| M_qd36_kr05 | 129.8 | deep resists vs qd12@kr05 (113@2000) |
| M_qd36_kr10 | 124.0 | still draining |

Deep trap **slows** tax (matches S06 vs S05) but does not reverse it at mid-flight.

## Interim judgment

1. **Radiance alone does not monodispersify** the quiet unitary clump — it bleeds inventory.  
2. Monodispersity (if reachable) needs a **closed I=O cycle** or geometric selection, not open O on a closed cell.  
3. R1 `*_rad` at T=10000 still the long-horizon Em\* test (likely SLOW_DRAIN continues).  
4. Internal sustain path remains **unitary + qd, k_rad=0** (or very mild).

*Final mono table when batch-2 + `.done`.*
