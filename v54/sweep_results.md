# V54 Parameter Sweep Results

## Sweep Configuration

- Grid: N=96, L=15, periodic BC
- Seed: init=braid, A=0.4, R_tube=3.0, ellip=0.3325, delta=(0, 3.0005, 4.4325)
- Sweep time: T=50 per trial
- 750 parameter combinations (6×5×5×5)
- Analysis: flood-fill cluster detection at 3× mean |φ|² threshold
- Score = compactness × density_contrast × √P_max (rewards tight, dense, bound)

## Parameter Grid

| Parameter | Values |
|-----------|--------|
| μ | -20, -30, -41.345, -60, -80, -100 |
| η | 0.2, 0.35, 0.5, 0.7, 1.0 |
| m² | 1.5, 2.25, 4.0, 6.25, 9.0 |
| κ | 20, 35, 50, 75, 100 |

## Top 20 Results

| Rank | Score | μ | η | m² | κ | Clusters | Mass | rms_r | Compact | Contrast | P_max |
|------|-------|------|------|------|------|---------|------|-------|---------|----------|-------|
| 1 | 82.7 | -80 | 0.20 | 1.50 | 50 | 36 | 51.5 | 1.88 | 0.838 | 101.9 | 0.936 |
| 2 | 78.0 | -60 | 0.35 | 9.00 | 20 | 4 | 51.0 | 1.67 | 0.857 | 111.8 | 0.662 |
| 3 | 77.6 | -60 | 0.20 | 9.00 | 20 | 4 | 50.8 | 1.66 | 0.857 | 111.5 | 0.659 |
| 4 | 77.5 | -60 | 0.50 | 9.00 | 20 | 4 | 50.9 | 1.67 | 0.856 | 111.4 | 0.658 |
| 5 | 74.3 | -60 | 0.70 | 9.00 | 20 | 4 | 49.8 | 1.67 | 0.856 | 108.9 | 0.634 |
| 6 | 72.2 | -100 | 0.35 | 4.00 | 20 | 5 | 17.4 | 1.38 | 0.881 | 100.0 | 0.671 |
| 7 | 71.8 | -100 | 0.20 | 1.50 | 75 | 34 | 52.7 | 1.94 | 0.833 | 94.7 | 0.826 |
| 8 | 70.9 | -80 | 0.70 | 4.00 | 20 | 14 | 30.7 | 1.50 | 0.871 | 97.9 | 0.690 |
| 9 | 64.8 | -60 | 1.00 | 9.00 | 20 | 4 | 43.3 | 1.67 | 0.856 | 101.5 | 0.555 |
| 10 | 64.4 | -100 | 0.50 | 9.00 | 20 | 3 | 46.6 | 2.20 | 0.810 | 99.5 | 0.638 |
| 11 | 62.5 | -60 | 1.00 | 6.25 | 20 | 4 | 51.4 | 1.71 | 0.853 | 95.6 | 0.586 |
| 12 | 59.4 | -41.3 | 1.00 | 4.00 | 20 | 2 | 103.9 | 2.14 | 0.816 | 86.2 | 0.713 |
| 13 | 58.3 | -80 | 1.00 | 9.00 | 35 | 4 | 46.6 | 1.77 | 0.848 | 96.4 | 0.507 |
| 14 | 56.5 | -41.3 | 0.35 | 1.50 | 20 | 7 | 38.5 | 1.91 | 0.836 | 85.4 | 0.628 |
| 15 | 56.2 | -80 | 1.00 | 9.00 | 20 | 2 | 42.2 | 1.73 | 0.851 | 90.8 | 0.527 |
| 16 | 55.6 | -100 | 0.70 | 4.00 | 20 | 27 | 9.6 | 1.22 | 0.895 | 88.6 | 0.490 |
| 17 | 54.6 | -100 | 0.50 | 2.25 | 50 | 23 | 26.4 | 1.72 | 0.852 | 87.1 | 0.539 |
| 18 | 53.5 | -80 | 0.50 | 4.00 | 20 | 9 | 31.2 | 1.51 | 0.870 | 82.6 | 0.553 |
| 19 | 53.0 | -60 | 0.70 | 6.25 | 20 | 4 | 45.2 | 1.67 | 0.857 | 85.0 | 0.529 |
| 20 | 52.1 | -100 | 0.70 | 9.00 | 20 | 3 | 44.6 | 2.46 | 0.788 | 89.7 | 0.542 |

## Key Patterns

1. **Low η dominates**: η=0.2-0.5 produces better particles than η=0.7-1.0.
   Less curl coupling = less theta drain = longer-lived structures.

2. **Two parameter families**:
   - Family A: m²=1.5, μ=-80 to -100, κ=50-75 (low mass, strong binding, high saturation)
   - Family B: m²=9.0, μ=-60, κ=20 (high mass, moderate binding, low saturation)

3. **κ=20 appears in most top results**: Lower saturation allows deeper V(P) binding.

4. **Compactness 0.81-0.90**: All top particles have rms_r=1.2-2.5 (box L=15).
   These are genuinely localized, not box-filling.

5. **Density contrast 82-112**: 100× denser than background at peak.

6. **Cluster count varies**: Family A (low m²) produces many small clusters (27-36).
   Family B (high m²) produces few large clusters (2-4).

## Post-Sweep Long Run

Best parameters (μ=-80, η=0.2, m²=1.5, κ=50) were used for T=3000 with
tune_dt=50 auto-tuning. The particle formed at T=50 but dissolved by T=2400
under periodic BC — same dissolution pattern as all previous runs. The auto-tuner
detected the particle (compact=0.84, contrast=102) but couldn't prevent dissolution.

## Conclusion

The sweep identifies WHERE in parameter space particles form most readily.
The top parameters produce particles with 100× density contrast and 0.84+
compactness at T=50. But no parameter combination tested prevents the
eventual dissolution under periodic BC.

## Next Step

Test whether an opposite-chirality particle can be constructed. If both
chiralities exist, a bound pair might recycle theta between them and
achieve true stability.

## Data

- Sweep log: `/root/run_sweep1.log` on GPU (ssh4.vast.ai:11784)
- Post-sweep SFA: `/root/sweep1.sfa`
- Diag: `/root/sweep1_diag.tsv`
