# Campaign goal function — final scoreboard

**Updated:** 2026-07-19 after `QUEUE_DONE`

| Component | Weight | Value | Note |
|-----------|--------|------:|------|
| S_force | 0.25 | **0.55** | Jobs OK; charge sign OK; no \(a_{\mathrm{rel}}\) tracks |
| S_pair | 0.25 | **0.65** | Rest holds; head-on \(E_{\mathrm{em}}\) drop |
| S_orbit | 0.20 | **0.50** | 3/3 complete; vt=0.08 soft-fail energy |
| S_Lhold | 0.15 | **0.70** | S4 \(E_{\mathrm{em}}\) held (≠ v79 death) |
| S_ledger | 0.10 | **1.00** | Gauss floor all runs |
| S_morph | 0.05 | **0.25** | No renders |
| **G** | 1.00 | **0.62** | |
| hard_fail | | **false** | |

## Thresholds

- **G ≥ 0.55** — product path alive → **HIT**  
- **0.30 ≤ G < 0.55** — retune only  
- **G < 0.30** — soft-kill multi-fab atoms  

## Hard fail triggers

1. Gauss left 1e-13 class — **no**  
2. Missing diags — **no** (16 runlogs + diags local)  
3. Disk full — **no**  
