# Multi-Z: park-aware cost + parked nuclear templates

## Problem (F15)

Z2–Z6 fulls never scorecard-PASS because **c_Q vs seed Q** treats nuclear
liquid-drop parking (seed 908 → ~640 for c6_light) as FAIL_Q, even when massL
is stable.

## A — Park-aware cost (controller)

\[
Q_{\mathrm{park}} = Q(t_{\mathrm{mid}}),\quad
c_Q^{\mathrm{park}} = \mathrm{clip}\big((Q_{\mathrm{park}}-Q_{\mathrm{end}})/Q_{\mathrm{park}},0,1\big)
\]

Use mid-run (t ≈ T/2) as park reference; score only **post-park drift**.
Seed→mid drop is nuclear physics, not atom failure.

Also keep legacy c_Q_seed for comparison in ledger.

## B — Parked multi-Z seed

1. Generate c6_light (6×ω=1.46 octa D=10) as C fabric  
2. Optionally pre-relax T_park=300 (or use v74 parked snapshot if available)  
3. Add L6 shell R=22 rest (or soft L)  
4. Score with park-aware cost  

## Campaign

| ID | Setup | Score |
|----|--------|-------|
| MZ1 | Z4 R22 rest, park-aware on existing pattern | re-score if diag exists; else short run |
| MZ2 | Z6 ω=1.46 + L6 R22 rest, park-aware full T=400 | production |
| MZ3 | If MZ2 massL stable + c_Q^park ≤0.15 → PASS multi-Z packaging |

## Relation to B4

B4 proves **single-C + L cloud**. Multi-Z asks whether a **parked droplet**
+ L cloud also packages under the same private-bag B1 laws.
