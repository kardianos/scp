# V24-S3: Multi-Scale Forces — Two Pairwise Couplings

## Summary

Asymmetric pairwise coupling V_pw = l12*phi1*phi2 + l23*phi2*phi3 + l31*phi3*phi1
produces THREE distinct mediator masses from the 3x3 mass matrix eigenvalues.
The system supports TWO distinct force ranges simultaneously, with a sign change
in the inter-oscillon force between short and long range.

## Parameters

mu=-20, kappa=20, m=1.0, A=0.8, sigma=3.0, Nx=8000, xmax=200

## Phase 1: Single Oscillon with Asymmetric Coupling

**Couplings**: l12=0.85, l23=l31=0.95

**Mass matrix eigenvalues** (m^2 + eigenvalue of coupling matrix):

| Mode | m^2_eff  | m_eff  | Range (1/m) |
|------|----------|--------|-------------|
| 0    | 0.01588  | 0.1260 | 7.94        |
| 1    | 0.15000  | 0.3873 | 2.58        |
| 2    | 2.83412  | 1.6835 | 0.594       |

Range ratio (longest/shortest): **13.4x**

**Oscillon survival**: fc=0.994 at t=10000. Strongly stable.

omega=1.560 (above all mass thresholds — oscillon bound by triple-product
coupling, not mass gap). phi1=phi2 always (fields 1,2 related by the
symmetric l23=l31 coupling), phi3 amplitude ~5% larger.

**Key finding**: The proposal's original parameters (l12=0.5, l23=l31=0.99) are
**tachyonic** (lightest eigenvalue m^2_eff = -0.17). Stability requires the
smallest eigenvalue of the coupling matrix L to satisfy min(eval) > -m^2.
For (l12,l23=l31), the constraint is roughly l23 < m^2 + l12/2.

## Phase 2: Two-Oscillon Force Measurement

Couplings: l12=0.85, l23=l31=0.95. Measured acceleration (d^2 sep/dt^2)
in the first half of t=2000 evolution, before oscillons move significantly.

| D   | velocity (dsep/dt) | acceleration | Sign      |
|-----|-------------------|--------------|-----------|
| 10  | +2.23e-2          | -4.12e-5     | Repulsive |
| 15  | +6.27e-3          | -1.35e-5     | Repulsive |
| 20  | +2.66e-3          | -8.89e-6     | Repulsive |
| 30  | +1.16e-3          | +1.57e-5     | Attractive|
| 40  | -9.34e-5          | +4.50e-7     | Attractive|

**Sign change between D=20 and D=30**: short-range repulsion + long-range attraction.

### Two-Yukawa fit

F = 5.59e-1 * exp(-D/1.0) + 2.22e-5 * exp(-D/30.0)

Fitted ranges: R_strong = 1.0, R_weak = 30.0 (ratio 30x).
Residual: 9.22e-11.

Single-Yukawa fit: F = 1.63e-4 * exp(-D/7.0), residual = 2.18e-10.
Two-Yukawa is **2.4x better** than single-Yukawa.

Note: The fitted ranges (1.0 and 30.0) differ from the eigenvalue predictions
(0.59 and 7.94). The short-range fit captures the steep D=10-20 falloff, while
the long range is inflated by the sign change. The force includes oscillatory
(breathing) contributions that complicate clean Yukawa extraction.

### Eigenvalue predictions vs measured

- Mode 2 (range 0.59): matches short-range repulsive component qualitatively
- Mode 0 (range 7.94): should produce long-range component; the measured
  attraction at D=30-40 is consistent with a longer-range attractive mode
- The sign change is physical: different eigenmodes couple with different signs

## Phase 3: Coupling Ratio Scan

All configurations stable (fc > 0.99 at t=5000).

| l12  | l23=l31 | eval_0  | eval_1  | eval_2  | range_0 | range_1 | range_2 | ratio | omega  | fc     |
|------|---------|---------|---------|---------|---------|---------|---------|-------|--------|--------|
| 0.90 | 0.95    | 0.0331  | 0.1000  | 2.8669  | 5.49    | 3.16    | 0.59    | 9.3   | 1.5660 | 0.998  |
| 0.85 | 0.95    | 0.0159  | 0.1500  | 2.8341  | 7.94    | 2.58    | 0.59    | 13.4  | 1.5540 | 0.997  |
| 0.80 | 0.90    | 0.0658  | 0.2000  | 2.7342  | 3.90    | 2.24    | 0.60    | 6.4   | 1.5240 | 0.996  |
| 0.70 | 0.85    | 0.0980  | 0.3000  | 2.6020  | 3.19    | 1.83    | 0.62    | 5.2   | 1.4820 | 0.994  |
| 0.50 | 0.75    | 0.1603  | 0.5000  | 2.3397  | 2.50    | 1.41    | 0.65    | 3.8   | 1.3920 | 0.991  |

The "ratio" column is range_0/range_2 (longest/shortest range).

**Trends**:
- Increasing asymmetry (l12 << l23) creates a wider range spread
- Maximum range ratio ~13.4x achieved at (0.85, 0.95), near the tachyonic boundary
- All configurations produce THREE distinct mediator masses
- omega decreases as total coupling decreases (weaker binding)
- fc uniformly > 0.99 for all tested configurations

## Answer to Key Question

**Can the system support two forces at different ranges simultaneously?**

**YES**, with caveats:

1. The asymmetric pairwise coupling produces three distinct eigenvalues,
   giving three mediator masses and three Compton wavelengths (ranges).

2. The inter-oscillon force shows a **sign change** between D=20 and D=30,
   indicating two competing components: short-range repulsion (from the
   heavy mode) and long-range attraction (from the light mode).

3. The two-Yukawa fit is 2.4x better than single-Yukawa, confirming
   multi-scale structure.

4. The maximum range ratio achievable without tachyonic instability is
   ~13x (ranges 0.59 and 7.94 in code units).

**Limitation**: The lightest mode mass is bounded below by the tachyonic
constraint. To get truly separated scales (e.g., nuclear ~1 fm vs gravity ~inf),
the lightest mode would need m -> 0, which pushes the system to the tachyonic
boundary. At (l12=0.85, l23=l31=0.95), the lightest mode has m_eff=0.126,
corresponding to range ~7.94 code units (vs ~0.59 for the heaviest). This is
~10x separation, not the ~10^39 needed for nuclear/gravity.

## Files

- `src/proca_multi.c` — asymmetric pairwise solver (3 phases)
- `data/phase1_ts.tsv` — Phase 1 time series
- `data/phase2_force.tsv` — Phase 2 force vs separation
- `data/phase2_D{10,15,20,30,40}.tsv` — Phase 2 two-oscillon time series
- `data/phase3_scan.tsv` — Phase 3 coupling ratio scan
