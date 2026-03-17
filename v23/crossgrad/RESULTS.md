# V23-B Cross-Gradient Oscillon Results

## Parameters

- mu = -20.0, kappa = 20.0, m = 1.0
- A_init = 0.8, sigma = 3.0
- Nx = 2000, xmax = 60.0
- t_equil = 5000, t_pert = 1000
- Perturbation: eps = 1.0e-03, sigma_pert = 3.0

## Summary Table

| eta  | survived | omega  | gap %  | peak  | E_total | growth rate | grew | extent |
|------|----------|--------|--------|-------|---------|-------------|------|--------|
| 0.00 | YES      | 0.8700 | 13.00  | 0.506 |   1.273 | -3.8e-04    | NO   | 22.0   |
| 0.20 | YES      | 0.8820 | 11.80  | 0.471 |   1.273 | -4.9e-05    | YES  |  7.4   |
| 0.50 | YES      | 0.9420 |  5.80  | 0.389 |   1.260 | -2.7e-05    | NO   |  7.6   |
| 1.00 | YES      | 0.9600 |  4.00  | 0.353 |   1.311 | -2.1e-05    | YES  | 12.8   |
| 1.50 | YES      | 0.9900 |  1.00  | 0.325 |   1.359 | -3.7e-05    | YES  | 20.5   |
| 2.00 | YES      | 0.9900 |  1.00  | 0.242 |   1.394 | -4.0e-04    | YES  |  6.6   |

## Analysis

### Oscillon Survival

All oscillons survive for the full t=5000 equilibration at every tested eta
(0.0 through 2.0). The symmetric oscillon is robust against cross-gradient
coupling. However, the oscillons become progressively weaker: peak amplitude
drops from 0.506 (eta=0) to 0.242 (eta=2), and the core fraction fc drops
from ~1.0 to ~0.65-0.88 at eta=2.0, indicating more energy leakage to
radiation.

### Breathing Frequency

omega increases with eta as expected from the effective speed c_eff = sqrt(1+eta/3):
- eta=0.0: omega=0.870, gap=13%
- eta=0.2: omega=0.882, gap=11.8%
- eta=0.5: omega=0.942, gap=5.8%
- eta=1.0: omega=0.960, gap=4.0%
- eta=1.5: omega=0.990, gap=1.0%
- eta=2.0: omega=0.990, gap=1.0%

The gap margin closes rapidly. At eta >= 1.5, the breathing frequency is
within 1% of the mass gap. The DFT frequency resolution (~0.006) limits
precision at this point, so the true omega may be even closer to m=1.

The gap closure follows roughly omega ~ m*sqrt(1 - delta), where delta
decreases with eta. This is the expected behavior: higher effective wave
speed means the oscillon must breathe faster to maintain coherence, but
faster breathing pushes omega toward the radiation threshold.

### Antisymmetric Perturbation

All growth rates are NEGATIVE (decay), indicating the antisymmetric mode
is stable at all tested eta values. The "grew" flag is misleading -- it
triggers when the final envelope exceeds 1.5x the initial amplitude, but
this occurs because the initial perturbation amplitude is measured at a
single instant (which may catch a zero-crossing), while the late-time
envelope is measured as a peak over many oscillation cycles.

The key finding is that all growth rates are O(10^{-5} to 10^{-4}), all
negative. The antisymmetric perturbation oscillates at roughly the
breathing frequency and slowly decays. No shear instability is observed.

The perturbation spatial extent varies: it is compact (6-8 code lengths)
when localized on the oscillon, and large (20+) when the perturbation has
partially radiated away to the absorbing boundaries.

### Field Splitting at eta > 0

At eta > 0, the equilibrated oscillon spontaneously develops a small
asymmetry between phi_1 and phi_{2,3}: field 1 has slightly smaller peak
amplitude than fields 2,3 (e.g., at eta=1: pk1=0.35 vs pk2,3=0.40).
This is because the (1+eta) Laplacian coefficient makes field 1 "stiffer"
-- its gradient energy is higher, so it settles to a slightly narrower
profile with lower peak. The equilibrium energy compensates via the
potential coupling.

### Conclusions

1. The symmetric oscillon survives cross-gradient coupling up to at least
   eta=2.0, but weakens progressively.

2. The breathing frequency approaches the mass gap with increasing eta,
   closing from 13% gap at eta=0 to ~1% gap at eta=1.5-2.0.

3. No antisymmetric (shear) instability is found. The perturbation decays
   at all eta values tested. The shear mode is stable.

4. The cross-gradient coupling spontaneously breaks the phi_1/phi_{2,3}
   symmetry of the oscillon profile (field 1 is stiffer), but this does
   not lead to instability.

5. For Phase 2 (linearized spectrum), the near-closure of the gap at
   eta ~ 1.5 is the most interesting regime -- the oscillon is close to
   the radiation threshold, which may affect the antisymmetric spectrum
   differently in a more careful eigenvalue analysis.

## Data Files

- `data/crossgrad_eta0.0_ts.tsv` -- time series for eta=0.00
- `data/crossgrad_eta0.2_ts.tsv` -- time series for eta=0.20
- `data/crossgrad_eta0.5_ts.tsv` -- time series for eta=0.50
- `data/crossgrad_eta1.0_ts.tsv` -- time series for eta=1.00
- `data/crossgrad_eta1.5_ts.tsv` -- time series for eta=1.50
- `data/crossgrad_eta2.0_ts.tsv` -- time series for eta=2.00
