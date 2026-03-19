# V32 SPH Two-Braid Gravity Test -- Phase 2 Results

## Configuration

- **Particles**: 242,326 total (240K uniform + 2,326 extra: 1,163 per braid)
- **Domain**: [-30, 30]^3 periodic in all directions
- **Braid 1**: centered at x = -10
- **Braid 2**: centered at x = +10 (initial separation D = 20)
- **Field**: 3-component scalar, m^2=2.25, mu=-41.3, kappa=50
- **SPH kernel**: Cubic spline, eta=1.2, adaptive smoothing length
- **Transport**: Energy-flux (Option 2), beta=0.01
- **Time**: T=300, CFL=0.3, dt adapted from 0.142 to ~0.14
- **Runtime**: ~4000 seconds (~67 min), ~2100 steps, 16 OMP threads
- **Separation tracking**: Energy-density-weighted centroid, threshold = 5x average E

## Key Result: D(t) DECREASES -- Weak Attractive Signal

The braid separation D(t) shows a slow but persistent decrease:

| Time window | D range | D mean | Trend |
|-------------|---------|--------|-------|
| T = 0-25    | 19.96-20.03 | 20.00 | Initial settling |
| T = 25-50   | 19.87-19.96 | 19.91 | Early decrease |
| T = 50-100  | 19.75-19.92 | 19.87 | Fluctuating low |
| T = 100-150 | 19.78-19.91 | 19.85 | Continuing decrease |
| T = 150-172 | 19.77-19.91 | 19.84 | Settled |

**Linear fit (T=0 to T=172):**
- D(t) = 19.978 - 0.00107 * t
- Rate: dD/dt = -0.00107 per time unit = **-0.107 per 100T**
- Total change: -0.18 over T=172 (from ~20.0 to ~19.8)

**Oscillation amplitude**: +/- 0.14 (peak-to-peak ~0.28)

**Signal-to-noise**: The linear trend (-0.18 over T=172) is comparable to the oscillation
amplitude (0.14). The signal is marginal -- roughly 1.3x the noise level.

## Detailed Separation Log

| T | D | x_left | x_right | n_left | n_right |
|---|---|--------|---------|--------|---------|
| 0 | 20.021 | -10.027 | 9.994 | 5041 | 3041 |
| 25 | 19.964 | -9.962 | 10.002 | 4466 | 1917 |
| 50 | 19.909 | -9.876 | 10.033 | 4627 | 1797 |
| 75 | 19.916 | -9.899 | 10.017 | 4280 | 1656 |
| 86 | 19.750 | -9.850 | 9.900 | 4209 | 1630 |
| 100 | 19.863 | -9.868 | 9.995 | 4207 | 1612 |
| 125 | 19.794 | -9.884 | 9.911 | 3951 | 1550 |
| 150 | 19.899 | -9.892 | 10.007 | 3799 | 1510 |
| 167 | 19.842 | -9.925 | 9.918 | 3749 | 1473 |
| 172 | 19.773 | -9.906 | 9.867 | 3729 | 1446 |

## Interpretation

### What the data shows

1. **D decreases by ~1%** over T=172 (20.02 to ~19.82), which is consistent with a weak
   attractive interaction between the two braids.

2. **Both centroids drift inward**, but not symmetrically:
   - Left braid: x moves from -10.03 to -9.91 (shift: +0.12, rightward)
   - Right braid: x moves from 9.99 to 9.87 (shift: -0.12, leftward)
   - The motion is roughly symmetric in the second half, though the left braid was initially
     more energetic (E_left ~ 5200 vs E_right ~ 1900).

3. **Particle count decays in both braids**: n_left drops from 5041 to 3729 (26% loss),
   n_right drops from 3041 to 1446 (52% loss). The right braid is dissipating faster,
   consistent with it being weaker initially. The total high-E particles (n_left + n_right)
   drops from 8082 to 5175 (36% loss).

### Why the signal is ambiguous

The 1% decrease in D could be:

1. **Gravitational attraction from metric contraction** (the hypothesis): Each braid creates
   a denser particle cloud, which provides better Laplacian resolution. When one braid's
   field extends into the other's denser cloud, the asymmetric neighbor count creates a
   net force toward the other braid.

2. **Centroid drift from braid dissipation**: As particles leak out of the high-E threshold
   zone, the centroid computation becomes noisier. If particles preferentially leak from
   the outward-facing side (which sees lower background density), the centroid could shift
   inward without any actual attraction.

3. **Asymmetric initialization artifact**: The left braid starts with 5041 high-E particles
   vs 3041 for the right braid. This 1.66x asymmetry (from the random seed) could bias
   the centroid measurements.

4. **Periodic boundary effect**: With domain L=30 and braids at +/-10, the periodic images
   are at distance 60-20=40 from each braid. The nearest image is 40 away (vs 20 between
   the actual braids), so PBC effects should be 4x weaker. But they could still contribute.

### Verdict: INCONCLUSIVE (leaning negative)

The signal (dD/dt = -0.001/T, total -1% over T=172) is:
- Comparable in magnitude to the oscillation noise (+/- 0.14)
- Confounded by braid dissipation (36% of high-E particles lost)
- Confounded by initial asymmetry (n_left >> n_right)

A clean gravitational signal would show:
- Steady, monotonic decrease in D (not oscillating)
- Acceleration (D'' < 0) as braids get closer
- No dependence on the energy threshold used for centroid computation

Instead we see noisy oscillations around a weakly declining mean. **The decline rate of
0.001/T is consistent with centroid bias from asymmetric particle loss, not gravitational
attraction.**

## Energy Conservation

Total energy grows from 11,611 to ~14,600 (26% increase over T=172). This is the same
Brookshaw Laplacian numerical heating seen in Phase 1, but milder because the larger domain
means lower average energy density (less numerical dissipation per particle).

## Technical Notes

- Source: `/home/d/code/scp/v32/src/v32_sph_two.c`
- Binary: `/home/d/code/scp/v32/v32_two`
- Data: `/home/d/code/scp/v32/data/sph2/`
- Separation log: `/home/d/code/scp/v32/data/sph2/separation.dat`
- Snapshots: `sph2_t{0000,0050,0100,0150,...}.bin` (21.3 MB each)

## Conclusion

The two-braid SPH simulation shows a **marginal decrease in braid separation** (~1% over
T=172), but the signal is comparable to centroid noise and confounded by particle
dissipation. The data does NOT provide convincing evidence for gravitational attraction from
geometric metric contraction. The SPH adaptive resolution mechanism creates particle
clustering (confirmed in Phase 1), but the inter-braid interaction from this clustering is
either absent or too weak to distinguish from numerical artifacts at the current resolution
and separation.

**To improve**: Run with symmetric initialization (equal particles in both braids), smaller
braid separation (D=10 instead of 20), and longer simulation time. If the signal is real,
halving the separation should quadruple the force (1/r^2 scaling), making it clearly visible
above the noise.
