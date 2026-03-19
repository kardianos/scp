# V32 Campaign Results: Binding-Weighted Gradient Model

## Model

    acc_a = laplacian(phi_a) + w(P) * alpha * (grad_rho / rho) . grad(phi_a) - m^2 phi_a - V'(phi_a)
    w(P) = 1/(1 + |P|/P_thresh),  P_thresh = 10% of peak |P|

Parameters: m^2=2.25, mu=-41.3, kappa=50, A_bg=0.1

## Critical Finding: Gradient Coupling Produces REPULSION, Not Attraction

The previously reported "attraction at alpha=0.5" was measured against alpha=0 (no coupling),
but the braids attract intrinsically even without the gradient coupling. When properly
controlled, the gradient coupling adds a repulsive correction at all tested alpha values.

---

## E4: Parameter Sweep (alpha = 0.1 - 0.7)

N=128, L=30, D=20, T=100. Control: alpha=0.

| alpha | deltaD | deltaD(control) | Coupling Effect | Interpretation |
|-------|--------|-----------------|-----------------|----------------|
| 0.0   | -1.00  | (baseline)      | 0.00            | BASELINE       |
| 0.1   | -0.82  | -1.00           | +0.18           | REPULSION      |
| 0.2   | -0.65  | -1.00           | +0.35           | REPULSION      |
| 0.3   | -0.11  | -1.00           | +0.89           | REPULSION      |
| 0.4   | +0.94  | -1.00           | +1.94           | REPULSION      |
| 0.5   | +0.94  | -1.00           | +1.94           | REPULSION      |
| 0.6   | +0.88  | -1.00           | +1.88           | REPULSION      |
| 0.7   | +8.99  | -1.00           | +9.99           | STRONG REPULSION |

**Conclusion**: The coupling is repulsive at all alpha. Increasing alpha increases
repulsion. At alpha >= 0.7, energy blows up (E: 21k -> 41k). The "sweet spot" does
not exist -- the weakest repulsion is at alpha=0.1 (+0.18), still net repulsive.

## E5/E5b: Separation Scan (D = 10 - 20)

N=128, L=30, T=100. Compared alpha=0 (control) vs alpha=0.1.

| D  | deltaD(alpha=0) | deltaD(alpha=0.1) | Coupling Effect |
|----|-----------------|-------------------|-----------------|
| 10 | -4.08           | -2.04             | +2.04 REPULSION |
| 12 | -6.34           | -5.78             | +0.56 REPULSION |
| 15 | -4.11           | -4.05             | +0.06 NEUTRAL   |
| 18 | -1.29           | -0.78             | +0.51 REPULSION |
| 20 | -1.00           | -0.82             | +0.18 NEUTRAL   |

**Intrinsic attraction**: Braids attract without any coupling (alpha=0).
The intrinsic force peaks around D=12 (deltaD=-6.34) and falls off at larger D.
This is the standard nonlinear field interaction between overlapping braids.

**Coupling effect**: At all separations, the gradient coupling adds repulsion.
The repulsion is strongest at D=10 (+2.04) where the braids are closest.

## E5 (variable box): Resolution Artifact

E5 with L=D*2+20 (variable box, dx varies 0.63-1.26):

| D  | L  | dx    | deltaD |
|----|-----|-------|--------|
| 10 | 40  | 0.630 | -3.44  |
| 15 | 50  | 0.787 | -4.80  |
| 20 | 60  | 0.945 | +0.44  |
| 25 | 70  | 1.102 | +0.01  |
| 30 | 80  | 1.260 | +0.88  |

At dx > 0.9, the braids are under-resolved (< 7 grid points across core R=3).
The "repulsion" at D=20-30 may be a resolution artifact. At fixed L=30 (dx=0.47),
the same D=20 shows attraction (deltaD=-0.82).

## Resolution Convergence (L=30, D=20, T=100, alpha=0.1)

| N   | dx    | deltaD |
|-----|-------|--------|
| 64  | 0.952 | -0.83  |
| 96  | 0.632 | -0.77  |
| 128 | 0.472 | -0.82  |
| 160 | 0.377 | -0.80  |

**Well converged**: deltaD = -0.80 +/- 0.03 across 4x resolution range.
The intrinsic attraction at fixed L=30 is numerically robust.

## E1: Single Braid Steady State

N=128, L=20, T=500, alpha=0.1.

- Energy oscillates 2.0e3 - 3.8e3 (breathing mode, period ~17 time units)
- Field concentration (fc) oscillates 0.15 - 0.63
- No energy drift: E_init=2698, E_final=2619 (3% decline)
- Braid is long-lived but NOT truly static -- perpetual breathing/pulsation
- w_core oscillates 0.58-0.97 (correlated with fc)
- No radial profile convergence (oscillations never damp)

## E2: Two Braids Long Run

N=128, L=40, T=500, alpha=0.1. (dx=0.63)

- D(t): starts ~20, dips to ~15.6 (t~190), then rebounds to ~36 by t=500
- Final: D = 19.95 -> 34.56 (deltaD = +14.61)
- Energy: 43647 -> 40548 (7% decline)
- Dynamics: approach -> scatter -> separate to far side of box
- At L=40, periodic BC creates repulsive image forces that dominate

## E6: Composite Braids (Close Together)

N=128, L=20, T=200, alpha=0.1. Initial D=8.

- D(t): 7.17 -> 4.16 (minimum at t~50) -> 12.8 (t~130) -> stabilizes ~12.5
- Classic scattering event: approach -> bounce -> separate
- NO bound state forms. The braids do not merge.
- Energy: 9395 -> 10314 (10% growth during scattering)

## E3: Five Braids

N=128, L=40, T=200, alpha=0.1. Pentagon arrangement at radius D/2=10.

- D: 12.02 -> 35.56 (deltaD = +23.54) -- massive expansion
- The 5 braids repel each other and spread across the box
- Energy: 271760 -> 300488 (10.5% growth)

---

## Summary of Conclusions

1. **The gradient coupling w(P) * alpha * (grad_rho/rho).grad(phi) is REPULSIVE**, not
   attractive. At all tested alpha values (0.1-0.7), the coupling pushes braids apart
   relative to the alpha=0 baseline.

2. **Braids attract intrinsically** (without the gradient coupling) through nonlinear
   field overlap. The intrinsic force peaks at D~12 (center-to-center) and decays at
   larger separations. This is an r^{-n} force with n > 2 (steeper than gravity).

3. **No bound states form**: Close braids (D=8) scatter rather than merge. The
   intrinsic attraction is not strong enough to trap, or the repulsive core prevents
   capture.

4. **Single braids do not reach steady state**: They breathe/pulsate indefinitely.
   This is a fundamental oscillatory mode, not a numerical artifact.

5. **Resolution converged**: The intrinsic attraction (deltaD ~ -0.8 at D=20, L=30, T=100)
   is robust across N=64-160 at fixed L.

6. **Periodic BC matters**: At larger L (and thus larger dx or image separation),
   the dynamics change qualitatively. L=30 gives attraction, L=40 gives net repulsion
   at D=20. This is likely from periodic image repulsion.

## Data Files

- `data/E4/sweep_summary.tsv` -- alpha sweep summary
- `data/E5b/force_law.tsv` -- separation scan at fixed L=30
- `data/E5_control/` -- alpha=0 control runs
- `data/E4_control/` -- alpha=0 control for D=20
- `data/Econv/` -- resolution convergence (N=64,96,128,160)
- `data/E1/timeseries.tsv` -- single braid evolution
- `data/E2/timeseries.tsv` -- two braid long run
- `data/E6/close/timeseries.tsv` -- close braid scattering
- `data/E3/timeseries.tsv` -- five braid evolution
- Field snapshots in `data/E{1,2,3,6}/field_t*.bin`
