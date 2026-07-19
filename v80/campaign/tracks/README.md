# Tracks campaign — force (F,R,L1b) + low‑vt orbit (O)

**Version:** v80 product continuation (AND scorecard)  
**After:** S4 COM window PASS (`../S4_COM_WINDOW_ANALYSIS.md`)  
**Policy:** **KEEP SFAs** + run `mf_pair_track` / `sfa_qball_track` before any prune.

## Jobs

| Phase | Jobs | T | Purpose |
|-------|------|---|---------|
| Force F | mf_force D=8,12,16,20,24 | 150 | \(a_{\mathrm{rel}}(D)\) opposite C–L |
| Force R | elite_same D=12,16,20 | 100 | same-sign control tracks |
| Orbit O | vt=0.03, 0.05 R=16 | 400 | multi-rev / capture |

Grid: N=128 L=32 B1 multi-fab (force/orbit); E-lite same for R.

## Outputs

- Remote work: keep `*.sfa`, `*_diag.tsv`, `*.runlog`
- Local: `/space/scp/v80/tracks/` (SFAs) + `results/` (diags, tracks TSVs)
- Analysis: `analyze_tracks.py` → `FORCE_ORBIT_SCORE.md`
