# Monodispersity scan — FINAL (8/8)

**Dir:** `mono_T2000_20260808_190046`  
**Setup:** quiet unitary blob, amp_nat=2, p_rad=4, T=2000, seed=20260802  
**Done:** 2026-08-08 ~20:22 CDT

## Question

Does k_rad × q_detune produce a **monodisperse Em\* fixed point**, or only SLOW_DRAIN?

## Full table (structural classes — auto Em-HOLD retired)

| arm | qd | k_rad | Em | max | nHot | slope | class |
|---|---:|---:|---:|---:|---:|---:|---|
| M_qd12_kr00 | 1.2 | 0 | **148.30** | 2.08 | 53 | ~0 | **HOLD** |
| M_qd12_kr02 | 1.2 | 0.02 | 124.81 | 1.52 | 31 | −0.0068 | **SLOW_DRAIN** |
| M_qd12_kr05 | 1.2 | 0.05 | 113.15 | 1.31 | 26 | −0.0080 | **SLOW_DRAIN** |
| M_qd12_kr10 | 1.2 | 0.10 | 105.36 | 0.84 | 21 | −0.0073 | **SLOW_DRAIN** peaks weaken |
| M_qd36_kr00 | 3.6 | 0 | **147.66** | 2.51 | 59 | −0.0003 | **HOLD** deep |
| M_qd36_kr02 | 3.6 | 0.02 | 132.14 | 1.32 | 66 | −0.0051 | **SLOW_DRAIN** (slower) |
| M_qd36_kr05 | 3.6 | 0.05 | 122.49 | 1.27 | 45 | −0.0069 | **SLOW_DRAIN** |
| M_qd36_kr10 | 3.6 | 0.10 | 116.30 | 1.07 | 39 | −0.0064 | **SLOW_DRAIN** |

## Judgment

1. **No Em\* park** at any (qd, k_rad) on quiet closed blob — tax is monotone in k_rad.  
2. Deep trap **slows** drain (~+10–15 Em retained vs law trap at same k_rad) but does not reverse slope.  
3. Higher k_rad → fewer nHot and lower Em_max (mildly more mono peaks) **while inventory bleeds** — not a stable monodisperse particle.  
4. Monodispersity, if reachable, needs **closed I=O** (e.g. P2 dual intake under test) or geometry selection — not open O alone.  
5. Internal sustain path remains **unitary + qd, k_rad=0**.

*Auto probe SUSTAIN_REPORT Em-only HOLD labels are incorrect for kr>0 arms; this file is the score of record.*
