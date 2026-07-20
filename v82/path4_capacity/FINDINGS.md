# path4_capacity FINDINGS

**Date:** 2026-07-20  
**Sandbox:** `bin/capacity_pic` (2D Yee + locks + `n_free`)

## Phase A — PASS

Pinned opposite pair, force probe with  
\(F = qE + k_{\mathrm{core}}(1-n_{\mathrm{free}}/n_{\mathrm{crit}})_+\nabla n_{\mathrm{free}}\).

- Zero-crossing of `F_along(D)` near **D\*≈5** for `k_core∈{2,8}`, `n_crit∈{0.3,0.6}`.  
- Non-monotonic force exists → **sign change at short range achieved**.  
- Artifact: at D=4 still attract (EM dominates / footprint scale); repulsive band at mid D then attract again at large D.

Artifact: `results/phaseA_scan.tsv`, `results/phaseA.log`.

## Phase B — FAIL (multi-rev)

Free pair D0=14, `k_core=2`, `n_crit=0.5`, vt scan, ~6000 steps:

| vt | sepmin–max | revs | note |
|----|------------|------|------|
| 0.08–0.15 | ~10–25 | ~0.18 | no collapse-through; slow expand |

Capacity **stops S1-style pass-through** (min sep≫0.7) but does not yet park multi-rev.

## Verdict

| | |
|--|--|
| P4 Phase A | **GREEN** |
| P4 Phase B | **not green** — continue tune / longer T / kernel 3D with capacity |
| Kernel port | **deferred** until orbit band + revs≥1 |
