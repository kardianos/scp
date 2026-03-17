# V25 Phase 6: Deformed Oscillon — Tidal Spin-2 Radiation

## Parameters
- mu=-20, kappa=20, mass=1.0, A=0.8, sigma=3.0
- lambda_pw=0.5, eta=0.1, lambda_L=0.1, alpha_g=0.001
- Grid: N=96, L=15, dx=0.316, dt=0.063
- t_equil=500, t_tidal=500
- R_gw=12 (strain measurement shell)
- Threads: 16

## Equilibration
- Oscillon forms and stabilizes by t=500
- Final: E=183, fc=0.993, peak=0.156, R_osc=3.845
- Pre-tidal quadrupole: Q22=3.1e-12, Q20=-7.4e-13 (essentially zero, as expected for spherical oscillon)

## Phase 6a: Static Tidal Deformation

Static tidal potential V = -eps_T/2 * (x^2 - y^2) * Sum_a phi_a^2 applied with Hermite ramp over t=0..50.

| eps_T | Q22_max | Q20_max | h+_eq | h+_pole | l0% | l1% | l2% | k2 | E_final | Survived? |
|-------|---------|---------|-------|---------|-----|-----|-----|----|---------|-----------|
| 0.001 | 1.17e+03 | 4.71e+02 | 8.9e-06 | 5.8e-05 | 93.5 | 5.3 | 1.3 | 1390 | 96.2 | Partial |
| 0.01 | 3.01e+03 | 2.47e+03 | 3.9e-04 | 2.5e-04 | 74.1 | 1.6 | 24.3 | 358 | 0.00 | Destroyed |
| 0.1 | BLOWUP | -- | -- | -- | -- | -- | -- | -- | NaN | NaN |

**Key findings:**
- eps_T=0.001: Oscillon partially survives (E: 183->96), strain is 93.5% monopole with only 1.3% l=2. The tidal field creates a quadrupole but it is overwhelmed by the breathing l=0 mode.
- eps_T=0.01: The tidal field DESTROYS the oscillon (E->0 by t~200). However, during the destruction, significant l=2 content appears (24.3%). The quadrupole signal (Q22~3000) comes from the asymmetric dispersal.
- eps_T=0.1: Immediate instability and blowup. The (x^2-y^2) potential acts as an anti-confining force that overwhelms the self-interaction.

**Angular pattern (static):**
- At eps=0.001: h+_pole (5.8e-5) > h+_equator (8.9e-6) -- OPPOSITE of spin-2.
- At eps=0.01: h+_equator (3.9e-4) > h+_pole (2.5e-4) -- marginal equatorial dominance, BUT h+_45 (8.1e-4) is largest, not consistent with pure spin-2 pattern.

## Phase 6b: Oscillating Tidal Field (Omega_T=0.1)

Oscillating tidal potential V = -eps_T/2 * cos(Omega_T*t) * (x^2-y^2) * Sum phi_a^2.

| eps_T | Q22_max | h+_eq | h+_pole | l0% | l2% | k2 | E_final |
|-------|---------|-------|---------|-----|-----|----|---------|
| 0.001 | 1.68e+03 | 4.0e-05 | 1.3e-04 | 96.9 | 0.6 | 1995 | 62.6 |
| 0.01 | 3.92e+03 | 9.8e-04 | 6.4e-04 | 73.6 | 26.0 | 466 | 0.01 |
| 0.1 | BLOWUP | -- | -- | -- | -- | -- | NaN |

The oscillating tidal field produces similar results to static: the oscillon is destroyed at eps=0.01, producing transient l=2 content during dispersal. The oscillating field does not help sustain the oscillon.

## Phase 6c: Resonance Scan (eps_T=0.01, varying Omega_T)

| Omega_T | Q22_max | Q22_rms | h+_eq | h+_pole | l0% | l2% | k2 | E_final |
|---------|---------|---------|-------|---------|-----|-----|----|---------|
| 0.05 | 3.36e+03 | 9.02e+02 | 4.7e-04 | 6.9e-04 | 80.3 | 17.7 | 400 | 0.01 |
| 0.10 | 3.92e+03 | 9.73e+02 | 9.8e-04 | 6.4e-04 | 73.6 | 26.0 | 466 | 0.01 |
| 0.20 | 1.15e+03 | 4.18e+02 | 1.4e-04 | 1.6e-04 | 82.9 | 16.3 | 137 | 35.6 |
| 0.50 | 3.41e+02 | 6.44e+01 | 4.5e-05 | 8.6e-05 | 94.8 | 5.0 | 40.6 | 107.1 |

**Key findings:**
- **Omega_T=0.1 gives maximum l=2 content (26%) and largest Q22_max** -- this is the resonant frequency. However, it destroys the oscillon.
- **Omega_T=0.5**: oscillon survives (E=107), but l=2 content is only 5%. The fast oscillation averages out before producing significant deformation.
- **Omega_T=0.2**: intermediate -- oscillon partially survives (E=35.6), l=2=16.3%.
- The deformability k2 peaks at Omega_T=0.1 (466) and decreases with increasing frequency.
- **Trade-off**: strong l=2 signal <=> oscillon destruction. No regime gives BOTH large l=2 AND oscillon survival.

## Tidal Love Number

k2 = Q22_max / (eps_T * R_osc^5) where R_osc = 3.845.

| Regime | k2 |
|--------|-----|
| Static, eps=0.001 | 1390 |
| Static, eps=0.01 | 358 |
| Driven, Omega=0.1 | 466 |
| Driven, Omega=0.5 | 40.6 |

The Love number is very large (k2 >> 1), meaning the oscillon is extremely SOFT -- it deforms readily under tidal forces. For comparison, neutron stars have k2 = 0.05-0.15 and black holes have k2 = 0. Our oscillon's k2 ~ 100-1000 indicates it has essentially no rigidity against quadrupolar perturbations.

The k2 depends on eps_T (nonlinear response) and Omega_T (frequency-dependent), which is expected for a non-rigid, nonlinear object.

## Angular Pattern Analysis

For the best case (eps=0.01, Omega=0.1):
- h+_equator = 9.8e-04
- h+_pole = 6.4e-04
- h+_45deg = 1.3e-03 (largest!)
- Ratio h+_eq / h+_pole = 1.53

For pure spin-2 (l=2, m=2): h+_equator should be maximal and h+_pole should be zero. Instead, h+ is NON-ZERO at the poles and peaks at 45 degrees. This is NOT the spin-2 angular signature.

The mixed angular pattern indicates:
1. The dominant signal is still l=0 (breathing mode, 73.6%)
2. The l=2 component (26%) is present but coexists with the isotropic component
3. The strain measurement captures both the background breathing AND the tidal deformation

## Interpretation

The tidal deformation experiment reveals a fundamental problem: **the oscillon is too soft to sustain a quadrupolar deformation**. At eps_T large enough to produce significant l=2 content, the tidal field destroys the oscillon entirely. At small eps_T, the l=2 signal is tiny (1%) compared to the l=0 breathing mode.

This is consistent with the Phase 2 result (l=2 = 0.5% for unperturbed oscillon). The oscillon's internal binding is dominated by the scalar (l=0) channel. There is no elastic shear modulus that would resist quadrupolar deformation -- the fields simply disperse along the stretched axis.

**Comparison to GR**: In GR, the star's self-gravity provides the restoring force against tidal deformation, with k2 ~ 0.1. Our oscillon has k2 ~ 100-1000, meaning it is ~1000x softer than a neutron star. The lack of an elastic shear mode means the oscillon cannot sustain steady-state quadrupolar oscillations.

## Conclusion

**The tidal deformation approach does not produce clean spin-2 radiation.** While l=2 content reaches 26% during oscillon dispersal, this is a transient effect associated with destruction, not a steady-state quadrupolar oscillation. The oscillating oscillon is intrinsically a scalar (l=0) radiator.

The key missing ingredient is a **shear restoring force** -- the elastic couplings (eta, lambda_L) provide anisotropic wave speeds but not a potential energy barrier against quadrupolar deformation. A truly spin-2 emitter would need the quadrupole mode to be a normal mode of the oscillon, not a destruction mode.

## Data Files
- `data/tidal_static_eps{E}.tsv` -- time series for static tidal runs
- `data/tidal_driven_eps{E}_omega{O}.tsv` -- time series for driven tidal runs
- `data/tidal_summary.tsv` -- summary table (omega scan)
- `data/tidal_love.tsv` -- summary table (eps scan, overwritten)
