# V47 — Topology-Dependent Coupling: Nuclear Binding Test

## Background

V46 analytical work established:
- η(P) = η₀ + η₁|P| with η₁ ≈ 115 (from per-voxel virial)
- η₁ ≈ √2|μ|/η₀ ≈ 117 (approximate relationship)
- η_eff at braid core ≈ 12 (above binding threshold η_crit ≈ 8)
- Predicted force curve: repulsive D<2, attractive D≈5-15, Coulomb D>20

## Experiment Design

Same setup as V45 but with η₁=115 added to the config:
- N=512, L=100, T=500, absorbing BC
- gen_deuterium seeds (same as V45 — allows direct comparison)
- COLZSTD output, auto-download

### Separation Matrix (reduced from V45 — focus on key distances)

| Run | D | Rationale |
|-----|---|-----------|
| D10 | 10 | Core boundary — where binding should appear |
| D15 | 15 | Interaction surface — V45 transition zone |
| D25 | 25 | Nuclear range — between V45 D20 and D40 |
| D80 | 80 | Non-interacting baseline |

### New config parameter
```
eta1 = 115.0
```
Added to all configs. All other parameters identical to V45.

## Execution Plan

1. Smoke test: CPU local (N=64, T=5) with η₁=115 — verify no crash
2. Smoke test: GPU remote (N=128, T=10) — verify CUDA kernel works
3. Verify frame count: confirm 6 frames (lround fix from V45)
4. Generate seeds (reuse V45 seeds at D=10,15,25,80)
5. Run 4 simulations sequentially on V100-32GB
6. Compare E_total(D) curve with V45 baseline
7. If binding detected at finite D: measure binding energy

## Success Criteria

- E_total(D=10 or D=15) < E_total(D=80) → **BINDING DETECTED**
- Force curve shows repulsion at D<5, attraction at D≈10-15 → **NUCLEAR FORCE**
- E_total(D=25) intermediate between D=15 and D=80 → **FORCE RANGE CORRECT**

## Cost Estimate

4 runs × 2.5 hr = 10 hr on V100-32GB ≈ $2
