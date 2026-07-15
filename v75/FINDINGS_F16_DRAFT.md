# F16 draft — B4 orbit + multi-Z park-aware (in progress)

**Date:** 2026-07-14/15  
**Depends on:** F15 B4_full PASS

## 1. Park-aware re-score of Stage-1 fulls (no new sims)

Using `v75/analysis/park_aware_score.py`:  
\(c_Q^{\mathrm{park}}=(Q_{\mathrm{mid}}-Q_{\mathrm{end}})/Q_{\mathrm{mid}}\).

| Trial | PASS_seed | PASS_park | cost_park | c_Q_park | c_L | Note |
|-------|-----------|-----------|-----------|---------|-----|------|
| B1a_full Z6 R22 | no | **no** | 0.215 | 0.224 | 0 | still post-park Q loss |
| **B2_full Z4 R18** | no | **YES** | 0.096 | 0.115 | 0.012 | multi-Z packaging under park score |
| **B3a_full Z2 R22** | no | **YES** | 0.057 | 0.123 | 0 | best multi-ball park score |
| B4_full Z1+L6 | YES | YES | 0.112 | 0 | 0 | F15 winner |

**Claim:** Soft-θ multi-Z (Z2–Z4) with stable L **already packages** under
park-aware cost; Stage-1 seed-based scorecard was too harsh.

## 2. B4 orbit campaign (running on v75st)

Frozen θ + tangential vt. Seeds: `gen_mf_shell_orbit`.

| ID | Setup | Status |
|----|--------|--------|
| b4o_pair_sub | C + 1 L @R=18, vt=0.7 v_c | running queue |
| b4o_pair_vc | vt=v_c≈0.071 | queued |
| b4o_pair_super | vt=1.3 v_c | queued |
| b4o_shell_vc | L6 co-rot vt=v_c | queued |

## 3. Multi-Z production (queued)

| ID | Setup |
|----|--------|
| mz2_z6_L6_R22 | c6 ω=1.46 octa + L6 R=22 rest; score park-aware |

## Tools added

- `sfa/seed/gen_mf_shell_orbit.c`
- `v75/analysis/park_aware_score.py`
- `v75/analysis/run_b4_orbit_mz.py`
- `v75/B4_ORBIT_PLAN.md`, `v75/MULTI_Z_PLAN.md`
