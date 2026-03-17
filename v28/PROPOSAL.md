# V28: Automated Search for Optimal Braided Soliton

## Goal
Find a field configuration that simultaneously satisfies:
1. **Spin-2 gravity**: l=2/l=0 > 20%
2. **EM charge**: torsion flux |Φ_T| > 0
3. **Topological confinement**: winding number conserved
4. **Localization**: field concentration fc > 0.3
5. **Stability**: no energy blowup over T=500

## Method
CMA-ES (Covariance Matrix Adaptation) + ASHA (Asynchronous Successive Halving)

### Parameter Space (16 dimensions)
**Initial condition (symmetry-breaking):**
- A1, A2, A3: per-field amplitudes (break 120° symmetry → nonzero Φ_T)
- delta2, delta3: phase offsets (includes 180°, 120°, arbitrary)
- R_tube: envelope width
- ellip, ellip_angle: elliptical cross-section (break azimuthal symmetry → l=2)
- k_factor: axial wavenumber (twist density)
- A_bg: background condensate amplitude
- A2_amp, A2_R: superimposed second braid layer

**Lagrangian:**
- mu, kappa: triple-product coupling
- mass: Klein-Gordon mass (0 = massless, bridges Regime 1↔2)
- lambda_pw: pairwise coupling

### Three-Tier ASHA
- **Tier 1 (screening)**: N=64, T=100, 1024 LHS candidates (~7 min)
- **Tier 2 (refinement)**: N=80, T=200, 4 CMA-ES populations × 5 gen × 16 (~10 min)
- **Tier 3 (validation)**: N=128, T=500, top 10 candidates (~7 min)

### LLM Training Analogies
- Large batch LHS = large batch SGD (explore broadly)
- ASHA tiers = learning rate scheduling (cheap→expensive)
- Multiple CMA-ES populations = population-based training (diversity)
- Post-hoc correlation analysis = feature importance (which params matter)

### Anchor Configurations (seeded into Tier 1)
0. V27 optimal: m=0, 120°, μ=-50, κ=50
1. V27 regime 1: m=1, 120°
2. 180° config: m=0, delta=(0,π)
3. 180° with mass: m=1
4. Asymmetric amplitudes: A=(1.2, 0.6, 0.3)
5. Elliptical braid: ellip=0.5
6. Intermediate mass with pairwise: m=0.5, λ_pw=0.5
7. Two-field dominant with background: A3=0.2, A_bg=0.3

## Output
- `data/tier1_results.tsv`: all 1024 screening results
- `data/tier2_results.tsv`: all CMA-ES results
- `data/tier3_results.tsv`: top 10 at full resolution
- `data/correlations.txt`: parameter-metric correlations
- Console: progress + top results + analysis
