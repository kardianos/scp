# v75 Status

**Updated**: 2026-07-15 — F15 done; F16 orbit + multi-Z park-aware **running**

## Closed

| Item | Result |
|------|--------|
| B1 MF + atom ladder F11–F14 | isolation OK; Z6+L6 rest atom fail |
| **Self-tune Stage 1 F15** | **PASS** — B4 single heavy + L6 shell |
| Git | `29530c2` MF B1 + F11–F15 |

## Park-aware re-score (F16 partial)

| Trial | PASS_park |
|-------|-----------|
| B2_full Z4 | **YES** |
| B3a_full Z2 | **YES** |
| B1a_full Z6 | no |
| B4_full | YES (also seed PASS) |

## Active campaign (`v75st`)

`run_b4o_mz_campaign.sh`: b4o_pair_{sub,vc,super}, b4o_shell_vc, mz2_z6_L6_R22  
Log: remote `/root/b4o_mz_campaign.log` · seeds `/space/scp/v75/b4o_mz/`

## Next after F16 runs

1. Analyze D(t) for sub/vc/super; massL integrity  
2. Park-aware score mz2; finalize FINDINGS F16  
3. Teardown GPU when idle  
