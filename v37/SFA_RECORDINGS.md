# SFA Recordings — All Surviving Geometries

Generated 2026-03-20. All runs: N=128, L=15, T=300, snap=3 (102 frames), eta=0.5, mt=0.

## Files

| File | Size | Frames | Geometry | Final E_pot | Status |
|------|------|--------|----------|-------------|--------|
| `data/sfa_braid3z.sfa` | 7.9 GB | 102 | braid3(z) single-axis | -30.5 | SURVIVED |
| `data/sfa_UUU.sfa` | 8.0 GB | 102 | crossed UUU (all right-handed) | -159.0 | SURVIVED |
| `data/sfa_UUD.sfa` | 8.0 GB | 102 | crossed UUD (mixed) | -64.4 | SURVIVED |
| `data/sfa_UDD.sfa` | 8.0 GB | 102 | crossed UDD (mixed) | -128.3 | SURVIVED |
| `data/sfa_DDD.sfa` | 8.0 GB | 102 | crossed DDD (all left-handed) | -109.0 | SURVIVED |
| `data/truncated_sfa.sfa` | 7.7 GB | 102 | truncated helix (pre-existing) | -25.9 | SURVIVED |

**Total disk: 48 GB**

## Survival Ranking (by final |E_pot|)

1. **UUU** (crossed, all right-handed): E_pot = -159, Ep_avg = 151.6 — STRONGEST
2. **UDD** (crossed, mixed): E_pot = -128.3, Ep_avg = 85.7
3. **DDD** (crossed, all left-handed): E_pot = -109.0, Ep_avg = 141.2
4. **UUD** (crossed, mixed): E_pot = -64.4, Ep_avg = 91.3
5. **braid3(z)** (single-axis): E_pot = -30.5 — weakest but still alive at T=300
6. **truncated** (compact helix, T=200 run): E_pot = -25.9

## View Commands

```bash
# Interactive 3D viewer (drag=rotate, scroll=zoom, Left/Right=frame, Space=play)
/home/d/code/scp/sfa/viewer/volview data/sfa_braid3z.sfa
/home/d/code/scp/sfa/viewer/volview data/sfa_UUU.sfa
/home/d/code/scp/sfa/viewer/volview data/sfa_UUD.sfa
/home/d/code/scp/sfa/viewer/volview data/sfa_UDD.sfa
/home/d/code/scp/sfa/viewer/volview data/sfa_DDD.sfa
/home/d/code/scp/sfa/viewer/volview data/truncated_sfa.sfa
```

## SFA Columns (all files)

| Column | Type | Semantic | Component |
|--------|------|----------|-----------|
| phi_x | F64 | POSITION | 0 |
| phi_y | F64 | POSITION | 1 |
| phi_z | F64 | POSITION | 2 |
| theta_x | F64 | ANGLE | 0 |
| theta_y | F64 | ANGLE | 1 |
| theta_z | F64 | ANGLE | 2 |

## Notes

- Crossed chirality runs (UUU, UUD, UDD, DDD) all survive T=300. The crossed geometry
  is significantly more robust than single-axis braid3(z).
- DDD was NOT part of the original T=200 screening (only UUU/UUD/UDD were run by the
  other agent). We recorded DDD directly as an SFA run and it survived strongly.
- UUU and DDD (both homochiral) show the strongest survival. The mixed variants
  (UUD, UDD) are weaker but still survive.
- braid3(z) is the weakest survivor: E_pot decayed to 41% of initial by T=300.
- truncated_sfa.sfa was pre-existing from an earlier run (T=200 only).
- All files use zstd compression (linked with -lzstd at compile time).

## Build

```bash
# SFA version of compact/braid3 code (supports braid3 geometry)
gcc -O3 -march=native -fopenmp -o v37_sfa src/v37_sfa.c -lzstd -lm

# SFA version of crossed chirality code
gcc -O3 -march=native -fopenmp -o v37_crossed_sfa src/v37_crossed_sfa.c -lzstd -lm
```
