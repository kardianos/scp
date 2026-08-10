# R1 endurance scoreboard — FINAL (8/8 T=10000)

**Dir:** `R1_T10000_20260808_172409`  
**Done:** 2026-08-09 ~00:00 CDT  
**Pass bars:** last-half |slope|≲0.02; Em_end≳100; nHot≳20; Em_max≳1.0

## Full table

| arm | Em_end | max | nHot | last-half slope | verdict |
|---|---:|---:|---:|---:|---|
| **R1_S00** | 148.30 | 1.71 | 52 | ~0 | **PASS HOLD** |
| **R1_S04** | 147.29 | 1.99 | 55 | ~0 | **PASS HOLD** |
| **R1_S14** | 147.71 | 1.56 | 50 | ~0 | **PASS HOLD_ext** fd=+741 |
| **R1_S16** | 148.07 | 1.66 | 50 | ~0 | **PASS HOLD** |
| **R1_S00_s111** | 146.85 | 1.68 | 53 | ~0 | **PASS HOLD** seed-robust |
| **R1_S00_amp12** | 206.75 | 1.80 | 106 | −1.4e-4 | **PASS HOLD** higher mass |
| R1_S00_rad | 105.41 | 0.75 | 18 | −5.3e-4 | **SLOW_DRAIN** (no Em* park; peaks weaken) |
| R1_S04_rad | 113.90 | 0.66 | 22 | −6.7e-4 | **SLOW_DRAIN** (deep resists ~+8 Em) |

## Mission answers

| claim | status |
|---|---|
| Quiet uni+qd sustains to T=10k | **YES** |
| Deep trap endures | **YES** |
| Progressive kx endures | **YES** HOLD_ext |
| seed_mw endures | **YES** |
| Seed robustness | **YES** s111 |
| Mass dependence amp=1.2 | **YES** HOLD at Em~207 |
| Radiance Em* park at long T | **NO** — SLOW_DRAIN; Em_max drops below 1 |

## Headline

> **Internal compact HOLD is real and long-lived** under unitary + q_detune (and progressive / winding / seed / mass variants).  
> Radiance alone is a tax, not a monodisperse fixed point, out to T=10000.
