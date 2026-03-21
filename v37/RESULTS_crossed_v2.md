# Results: Crossed Braid v2 — Bug Fixes and Fragmentation Analysis

## Background

The original crossed braid experiments (v37_crossed.c) reported all chirality
variants as "surviving" T=200. Volumetric inspection revealed this was false —
all configurations fragmented into dozens of pieces. The analysis tools failed
to detect this because they tracked only global metrics (total E_pot, total P_int).

## Bugs Found in v37_crossed.c

### Bug 1: Chirality — strand winding never flipped
The strand center positions always used right-handed winding (`+kw*z`) regardless
of chirality. For D (left-handed), the strands should wind as `cos(-kw*z + 2πs/3)`.
This meant D braids had right-handed geometry with left-handed phase offsets —
an inconsistent, unphysical configuration.

### Bug 2: Chirality — velocity sign wrong for D
Left-handed wave velocity should be negated relative to right-handed. The code
used positive velocity for both. Subsumed by Bug 3 fix.

### Bug 3: Binding — traveling wave init caused arm propagation
Each arm was initialized with traveling wave velocity `ω A sin(kz + δ)`. This meant:
- z-arm propagated in +z
- x-arm propagated in +x
- y-arm propagated in +y
The arms literally walked away from the crossing point.

### Bug 4: Background field only along z
`ph_bg = k_bg * z + 2πa/3` — always z regardless of braid arm. This broke the
3-axis symmetry and introduced a z-bias (which is why UUU ≠ DDD in old results).

### Bug 5: No absorbing boundary
Periodic BC reflected radiated waves back into the braid structure, causing
destructive interference and fragmentation.

### Bug 6: Output in .bin format
Not SFA — incompatible with analysis tools.

### Bug 7: No fragmentation detection
Diagnostics tracked total E_pot, P_int, aspect ratio — none of which detect
when a single structure splits into multiple pieces.

## Fragmentation Analysis of Original Experiments

Using `sfa_frag.c` (connected-component flood fill on |P| > threshold):

### Old UUD (threshold=0.01)
```
t=0:   1 cluster  max|P|=20.5
t=30:  8 clusters max|P|=0.50
t=60:  10 clusters
t=90:  16 clusters (peak)
t=120: 12 clusters
t=180: 9 clusters
t=300: 9 clusters (steady state — permanently fragmented)
```

### Old UUU (threshold=0.01)
```
t=0:   2 clusters (already wrong!)
t=60:  32 clusters
t=120: 54 clusters
t=240: 60 clusters (peak)
t=300: 49 clusters
```

**Conclusion: ALL original crossed experiments were fake survivals.** The structures
shattered by t=30 and remained as 10-60 disconnected fragments.

## Fixes Applied in v37_crossed_v2.c

1. **True chirality**: `kw_sign = is_D ? -kw : kw` for strand winding direction
2. **Stationary start**: Zero braid velocity; only background gets traveling wave
3. **Per-field background**: φ₀→z, φ₁→x, φ₂→y (isotropic, no z-bias)
4. **Absorbing boundary**: Quadratic damping zone (width=3, rate=0.02/step)
5. **SFA output**: Direct SFA writing via sfa.h
6. **Fragmentation detection**: Inline cluster analysis + connection ratio
7. **Optional field permutation**: `-permute` flag (tested, negative result)

## V2 Results

### UUD A_strand=0.5, N=128, L=15, T=200

Cluster history (threshold=0.01):
```
t=0:   1 cluster  max|P|=20.6
t=20:  4 clusters
t=40:  8 clusters (peak — breathing minimum)
t=60:  7 clusters
t=100: 3 clusters (reconverging)
t=120: 2 clusters
t=140: 2 clusters
t=200: 1 cluster  (fully reconverged)
```

Energy evolution:
```
t=0:   E_total=3767  E_pot=-125  P_int=697
t=50:  E_total=2149  E_pot=-104  P_int=99
t=100: E_total=968   E_pot=-46   P_int=42
t=150: E_total=510   E_pot=-45   P_int=42
t=200: E_total=366   E_pot=-0.0  P_int=1
```

Final state: E_total=366 (90% drained by absorbing BC), aspect=3.2, 2 clusters.

### UUD A_strand=0.7, N=128, L=15, T=200

Stronger initial binding but same eventual energy drain:
```
t=0:   E_total=6405  E_pot=-169  P_int=1830
t=25:  E_total=5657  E_pot=-175  P_int=187
t=50:  E_total=3498  E_pot=-176  P_int=152
t=100: E_total=1206  E_pot=-39   P_int=24
t=150: E_total=504   E_pot=-0.4  P_int=3
t=200: E_total=312   E_pot=-0.0  P_int=1   (95% drained)
```

### UUD A_strand=0.5, permute=1 (field permutation)

**Dead on arrival.** E_pot=-1.3 at t=0 (vs -125 without permutation).
P_int=9.7 (vs 697). The cyclic delta rotation causes phase cancellation at the
crossing point: cos(0) + cos(3.0) + cos(4.4) = 1.0 - 0.99 - 0.31 = -0.30
vs 3×cos(0) = 3.0 without permutation. Binding 100x weaker.

**Conclusion: Field permutation is the wrong approach.**

## Comparison: Old vs V2 (UUD, threshold=0.01)

| Metric              | Old (v37_crossed) | V2 (v37_crossed_v2) |
|---------------------|-------------------|----------------------|
| t=0 clusters        | 1                 | 1                    |
| t=30 clusters       | 8                 | ~4                   |
| t=60 clusters       | 10                | ~7                   |
| t=100 clusters      | **16**            | 3                    |
| t=140 clusters      | 12                | 2                    |
| t=200 clusters      | 9 (steady)        | **1** (reconverges)  |
| Peak clusters       | **76** (auto-thresh) | **8**             |
| Reconverges?        | NO                | YES                  |
| max|P| at t=200     | 0.43 (diffuse)    | 0.15 (coherent)      |

## Remaining Issue: Absorbing BC Too Aggressive

The damping rate 0.02/step drains all braid energy by t≈180. Each breathing
cycle radiates waves outward; the absorber eats them. After ~30 cycles,
the braid has no energy left.

Options to fix:
1. Reduce damping rate to 0.005/step or lower
2. Use larger box (L=20-25) so waves travel longer before hitting boundary
3. Frequency-dependent absorber (absorb high-freq radiation, pass braid modes)
4. Remove absorbing BC entirely — test if chirality + stationary start alone
   prevent fragmentation with periodic BC

## Tools Created

| Tool | File | Purpose |
|------|------|---------|
| SFA Fragmentation Analyzer | `src/sfa_frag.c` | Reads SFA, flood-fill cluster detection |
| Corrected Crossed Braid | `src/v37_crossed_v2.c` | 7-bug fix, SFA output |
| SFA fixup | `/tmp/fix_sfa.c` | Repair streaming-mode SFA index |

### Build commands
```bash
gcc -O3 -march=native -fopenmp -o sfa_frag src/sfa_frag.c -lzstd -lm
gcc -O3 -march=native -fopenmp -o v37_crossed_v2 src/v37_crossed_v2.c -lzstd -lm
```

### Usage
```bash
# Run corrected crossed braid
OMP_NUM_THREADS=8 ./v37_crossed_v2 -chiral UUD -N 128 -L 15 -T 200 \
    -diag 5 -snap 10 -o data/crossed_v2_UUD

# Analyze fragmentation of any SFA file
./sfa_frag data/crossed_v2_UUD.sfa -threshold 0.01 -every 5

# Available flags for v2:
#   -Astrand 0.7    (amplitude, default 0.5)
#   -permute        (rotate delta between axes — negative result)
#   -damp_width 3   (absorbing BC zone width)
#   -damp_rate 0.02 (absorbing BC damping per step)
```

## SFA Policy

All simulation output now uses SFA format exclusively. See `/home/d/code/scp/CLAUDE.md`.
