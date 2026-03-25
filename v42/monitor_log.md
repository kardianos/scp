# Deuterium Simulation Monitor Log

## Instance Details
- Instance ID: 33489930
- GPU: Tesla V100-SXM2-32GB
- Simulation: deuterium N=512, T=500

## Status Checks

| Check | Status | Notes |
|-------|--------|-------|
| #1 | t=10/500 | 2% complete |
| #2 (+10m) | t=50/500 | 10% complete, ~4 steps/min |
| #3 (+20m) | t=85/500 | 17% |
| #4 (+30m) | t=115/500 | 23% |
| #5 (+40m) | t=155/500 | 31% |
| #6 (+50m) | t=190/500 | 38% |
| #7 (+60m) | t=220/500 | 44% |
| #8 (+70m) | t=260/500 | 52% |
| #9 (+80m) | t=295/500 | 59% |
| #10 (+90m) | t=325/500 | 65% |
| #11 (+100m) | t=365/500 | 73% |
| #12 (+110m) | t=395/500 | 79% |
| #13 (+120m) | t=430/500 | 86% |
| #14 (+130m) | t=470/500 | 94% |
| #15 (+135m) | t=490/500 | 98% |
| #16 (+137m) | t=495/500 | 99% |
| #17 (+139m) | t=500/500 | 100%, writing final output |
| #18 (+140m) | SIM_DONE | Simulation complete |

## Phase 2: Remote Analysis

Ran on remote after simulation completed:
- `analyze_sfa` -> deuterium_analysis.json (1198 bytes)
- `freq_phase` -> deuterium_freq.json (1374 bytes)
- `conv` (f32->f16) -> deuterium_f16.sfa (9,097,794,041 bytes)

Key results from analysis:
- S_mean = 0.6018, S_final = 0.5111
- E_pot(final) = -53.45
- P_int retention = 44.3%
- Cluster count: 9 -> 9 -> 7 -> 9 -> 8 -> 13 -> 13

## Phase 3/4: Download Verification

| File | Remote Size | Local Size | Status |
|------|-------------|------------|--------|
| deuterium_output.sfa (f32) | 29,569,722,810 | 29,569,722,810 | VERIFIED (via rsync) |
| deuterium_analysis.json | 1,198 | 1,198 | VERIFIED |
| deuterium_freq.json | 1,374 | 1,374 | VERIFIED |
| deuterium_diag.tsv | 17,975 | 17,975 | VERIFIED |
| deuterium_f16.sfa (re-conv) | 9,097,794,041 | 4,600,035,883 | INCOMPLETE (connection dropped) |
| deuterium_f16.sfa (orig rsync) | 4.3 GB | 4.3 GB | VERIFIED (from earlier rsync, frame-skipped) |

## Phase 5: Instance Status

Instance 33489930 was already terminated (no active instances found via vastai CLI).
No destroy needed -- instance gone before full f16 re-download completed.

## Summary

ALL critical data is safe locally:
- Full f32 SFA (29.5 GB) -- VERIFIED, authoritative source
- Analysis JSON -- VERIFIED
- Frequency/phase JSON -- VERIFIED
- Diagnostics TSV -- VERIFIED
- F16 SFA (4.3 GB, frame-skipped version from earlier rsync) -- valid SFA header confirmed

The incomplete 9.1 GB f16 (all-frames conversion) failed due to connection drop.
This is not a data loss -- the full f32 is available locally and can be re-converted.
