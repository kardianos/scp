# Deuterium Output Rsync Log

Remote: `root@ssh5.vast.ai:19930`
Started: 2026-03-24

## Sync History

| Cycle | Time | f32 Size | Diag Rows | Sim Time | Status |
|-------|------|----------|-----------|----------|--------|
| 1 | 20:47 | 31 MB | 4 | t=15 | OK - initial sync |
| 2 | 20:57 | 31 MB | 12 | t=55 | no new SFA frames yet (snap_dt=100) |
| 3 | 21:07 | 31 MB | 20 | t=95 | approaching first snap at t=100 |
| 4 | 21:17 | 2.7 GB (partial) | 29 | t=140 | broken pipe at 59%, remote=4.8GB |
| 4b | 21:28 | 4.5 GB | 32 | t=155 | resume OK, caught up to remote |
| 5 | 21:38 | 4.5 GB | 40 | t=195 | no change, t=200 snap imminent |
| 6 | 21:48 | 9.1 GB | 54 | t=265 | t=200 frame downloaded, caught up |
| 7 | 22:04 | 13.7 GB | 73 | t=360 | t=300 frame downloaded (4 of 6 frames) |
| 8 | 22:14 | 13.7 GB | 80 | t=395 | t=400 snap imminent, 75% sim complete |
| 9 | 22:30 | 18.3 GB | 91 | t=450 | t=400 frame downloaded (5 of 6 frames) |
| 10 | 22:40 | 18.3 GB | 99 | t=490 | sim nearly done, final snap at t=500 |
| 11 | 22:50 | 27.5 GB | 101 | t=500 | **SIM DONE**, f32 complete, analysis.json synced |
| 12 | 23:44 | 27.5 GB | 101 | t=500 | **ALL FILES SYNCED** f16 (4.3GB) + freq.json |

## Simulation Complete
- Wall time: 146.5 min (8789.6s), 172 ms/step
- Final file: 29,569,722,810 bytes (27.5 GB), 7 frames
- P_int retention: 44.3%, S_final: 0.511, aspect: 1.013
- Energy drift: -63.6% (E: 276,400 -> 100,610)

## Files Downloaded
| File | Size | Status |
|------|------|--------|
| deuterium_output.sfa (f32) | 27.5 GB | COMPLETE |
| deuterium_f16.sfa | 4.3 GB | COMPLETE |
| deuterium_diag.tsv | 18 KB | COMPLETE |
| deuterium_analysis.json | 1.2 KB | COMPLETE |
| deuterium_freq.json | 1.4 KB | COMPLETE |

## Config Summary (from deuterium_config.cfg)
- N=512, L=100, T=500, snap_dt=100, diag_dt=5, precision=f32
- Frame size: ~4.8 GB (compressed), Total frames: 6 (t=0,100,200,300,400,500)
- Expected total: ~28 GB
- Rate: ~40 sim-time / 10 min = ~2h remaining (t=265 of 500)
