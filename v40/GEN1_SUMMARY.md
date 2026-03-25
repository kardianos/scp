# Generation 1 Summary: Structural Variations

**Date**: 2026-03-23
**Grid**: N=128, L=15, T=20, dt_factor=0.025
**Physics**: m²=2.25, m_θ²=0, η=0.5, μ=-41.345, κ=50
**Baseline**: Single V34 braid (T=200) — S_final=0.81, S_mean=0.62, E_pot=-87

## Strategy

8 candidates testing single structural modifications to the baseline braid:
amplitude scaling, multi-braid configurations, chirality flips, geometry changes.

## Results

| Rank | ID | Description | S_final | S_mean | P_ret | E_pot | Verdict |
|------|-----|-------------|---------|--------|-------|-------|---------|
| 1 | 004 | Scale 1.5× | **1.475** | 1.067 | 63% | -233 | Winner — amplitude drives binding |
| 2 | 003 | Counter-braid (flipped chirality) | **1.240** | 0.919 | 76% | -179 | Winner — chirality pair binds |
| 3 | 002 | Perpendicular braids (90°) | **1.082** | 0.766 | 77% | -141 | Viable — orthogonal interference |
| 4 | 006 | Braid + oscillon | 0.758 | 0.578 | 73% | -86 | Marginal — oscillon weakly coupled |
| 5 | 001 | Two parallel braids (sep=8) | 0.732 | 1.090 | 54% | -82 | Peaks early, decays — leakage |
| 6 | 007 | High ellipticity (0.5) | 0.562 | 0.708 | 61% | -47 | Below baseline |
| 7 | 008 | Low background (A_bg=0.03) | 0.527 | 0.673 | 62% | -43 | Below baseline |
| 8 | 005 | Scale 0.7× | 0.143 | 0.235 | 22% | -3 | Dead — amplitude too low |

## Key Findings

### 1. Amplitude is the primary control knob

C4 (1.5×) achieved the highest score. The potential V(P) = μP²/(1+κP²) with μ=-41.345
rewards larger P. Since P = φ₀φ₁φ₂, higher amplitude means deeper potential wells.
Conversely, C5 (0.7×) collapsed — below a threshold amplitude (~0.6), the potential
well is too shallow to sustain binding against dispersive losses.

### 2. Counter-chirality pairing creates durable binding

C3 (braid + counter-braid) had the best P_int retention (76%) among the multi-structure
candidates. The counter-rotating traveling waves create constructive interference in
the triple product P, enhancing binding. Same-chirality parallel braids (C1) peaked
early but decayed — the same-direction waves create a node/antinode pattern that
eventually separates.

### 3. Non-collinear structures need more than superposition

C2 (perpendicular braids) showed moderate binding (S=1.08) but the structure was
caught at a breathing peak — the T=10 data had shown dispersal. The perpendicular
arrangement creates transient constructive interference but lacks a sustained binding
mechanism. This motivated the Gen 2 xyz 3-braid approach.

### 4. Geometry modifications alone don't help

C7 (high ellipticity) and C8 (low background) both underperformed the baseline.
Changing the braid's shape without increasing energy density or adding structural
complexity weakens the structure.

## Survivors for Gen 2

1. **004** — Scale 1.5× (highest S_final)
2. **003** — Counter-braid (best retention)
3. **007** — High ellipticity (user selected for exploration)
4. **002** — Perpendicular braids (inspired 3-braid xyz approach)

## Files

```
v40/gen_001/
  candidate_001/ — Two parallel braids (sep=8)
  candidate_002/ — Perpendicular braids (90°)
  candidate_003/ — Counter-braid (flipped chirality)
  candidate_004/ — Scale 1.5×
  candidate_005/ — Scale 0.7×
  candidate_006/ — Braid + oscillon
  candidate_007/ — High ellipticity (0.5)
  candidate_008/ — Low background (A_bg=0.03)
Each contains: config.cfg, seed.sfa, output.sfa, diag.tsv, analysis.json
```
