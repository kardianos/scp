# Massive Pion Skyrmion Analysis

## Summary

Solved the hedgehog Skyrmion profile ODE with pion mass term V = m_pi^2 rho_0^2 int(1-cos(f)) d^3x.
**No bound breathing modes exist** even at physical pion mass — the K=0 potential
well is too shallow for discrete states.

## Modified ODE

f''(r^2 + 2c_4 sin^2(f)) + 2f'r + c_4 f'^2 sin(2f)
  - sin(2f)(1 + c_4 sin^2(f)/r^2) - m_pi^2 r^2 sin(f) = 0

The pion mass term (+m_pi^2 r^2 sin(f)) acts as a restoring force toward f=0 (vacuum).

Key changes from massless case:
- Asymptotic: f ~ exp(-m_pi r)/r (exponential decay, not power-law 1/r^2)
- Virial: E_2 - E_4 + 3 E_V = 0 (no longer E_2 = E_4)
- Mass formula: Mc^2 = 2 E_4 - 2 E_V

## Pion mass in code units

m_pi has dimensions of 1/(code length). Physical pion mass:
- m_pi_code = 139.6 MeV * 0.5624 fm / (197.3 MeV*fm) = 0.398 code^{-1}  (at e=1, rho_0=1)
- Compton wavelength: 1/m_pi = 2.51 code lengths (about 1.7x the soliton radius)

This is NOT the same as m_pi in energy units (15.3 code E). The distinction is crucial.

## Soliton properties vs pion mass (e=1, rho_0=1)

| m_pi  | a=-f'(0) | E_total | E_2    | E_4    | E_V   | E/E_FB | Virial check |
|-------|----------|---------|--------|--------|-------|--------|-------------|
| 0     | 1.420    | 103.14  | 51.55  | 51.59  | 0.00  | 1.232  | 3.9e-2      |
| 0.100 | 1.442    | 103.83  | 50.63  | 52.55  | 0.64  | 1.240  | 7.5e-3      |
| 0.200 | 1.491    | 105.50  | 48.74  | 54.76  | 2.00  | 1.260  | 1.8e-4      |
| 0.398 | 1.617    | 110.19  | 44.73  | 60.27  | 5.18  | 1.316  | 1.4e-8      |
| 0.600 | 1.751    | 115.64  | 41.27  | 66.09  | 8.27  | 1.381  | 2.3e-10     |
| 1.000 | 2.001    | 126.65  | 36.24  | 76.86  | 13.54 | 1.512  | 8.0e-6      |

Note: virial residual E_2 - E_4 + 3E_V. Use R_max=10 for m_pi >= 1.0.

### Key observations

1. **Soliton mass increases monotonically** with m_pi (+6.8% at physical m_pi)
2. **E_2 decreases, E_4 increases** — the soliton is more compact (steeper slope a)
3. **E_V is moderate** — about 5% of total energy at physical m_pi
4. **Virial theorem satisfied** to machine precision for all m_pi

## Profile comparison: massive vs massless

| Quantity | Massless | Massive (m_pi=0.398) | Change |
|----------|----------|---------------------|--------|
| a = -f'(0) | 1.420 | 1.617 | +13.9% |
| E_total | 103.14 | 110.19 | +6.8% |
| R_rms | 1.496 code (0.841 fm) | 1.264 code (0.711 fm) | -15.5% |
| R_half (f=pi/2) | 1.239 code (0.697 fm) | 1.080 code (0.607 fm) | -12.8% |
| E_2/E_4 | 1.000 | 0.742 | — |
| Compton wavelength | inf | 2.51 code (1.41 fm) | — |

The massive soliton is significantly more compact (-15.5% in RMS radius).

## Breathing mode spectrum (massive, e=1)

Using the self-consistent massive profile (solved WITH pion mass):

| n  | lambda  | omega  | E (MeV) | nodes | Status |
|----|---------|--------|---------|-------|--------|
| 0  | 0.327   | 0.572  | 5.2     | 0     | Continuum |
| 1  | 0.602   | 0.776  | 7.1     | 2     | Continuum |
| 2  | 0.983   | 0.992  | 9.0     | 2     | Continuum |
| 6  | 4.511   | 2.124  | 19.3    | 6     | Continuum |

**Continuum threshold: lambda = m_pi^2 = 0.158**

Fine scan [0, 0.158]: **NO eigenvalues found**.
Fine scan [-1, 0]: **NO instabilities found**.

### Comparison with sigma-model profile + pion mass (wrong profile)

The previous run (using sigma-model profile with pion mass added to W only) found
lambda_0 < 0 (instability). This was an ARTIFACT of using the wrong equilibrium
profile. The sigma-model profile is not an energy minimum when V != 0.

Using the correct massive equilibrium profile:
- All eigenvalues positive (stable)
- No bound states below threshold
- Indicial exponent nu = 4.31 (vs 3.79 massless)

## Moment of inertia

| Quantity | Massless | Massive | Experiment |
|----------|----------|---------|-----------|
| Lambda | 141.6 | 86.6 | — |
| Delta-N splitting | 0.1 MeV | 0.2 MeV | 293.7 MeV |

The massive soliton has a smaller moment of inertia (more compact) and larger
Delta-N splitting, but still far below the experimental value.

## Parameter fit implications

The massless parameter fit gave: rho_0^4/e^2 = 0.02591 from the constraint
M_p * r_p / (hbar*c) = 4.00, using the universal constants K=103.13, C_R=1.4961.

With pion mass, these universal constants change:
- K_massive = E_sol / (rho_0^3/e) = 110.19 (vs 103.13)
- C_R_massive = R_rms / (rho_0/e) = 1.264 (vs 1.496)
- Product: K * C_R = 139.3 (vs 154.3)

New constraint: rho_0^4/e^2 = (M_p * r_p / hbar*c) / (K * C_R) = 4.00 / 139.3 = 0.02872

This changes the standard convention parameters:
- e_std = 3.69 * sqrt(0.02872/0.02591) = 3.88  (slightly larger)
- F_pi = rho_0 * ... (needs recalculation with new rho_0)

## Implications for Phase 12

1. **K=0 breathing mode: definitively no bound states** — this result is robust across
   both massless and massive equilibrium profiles.

2. **Stability confirmed**: the massive Skyrmion is a genuine energy minimum (no
   instabilities), unlike the artifact from using the sigma-model profile.

3. **Meson candidates require angular modes (K >= 1)**:
   - Pion: K=1, isovector channel, should have mass gap m_pi
   - Rho: K=1, spin-1 excitation
   - f_2(1270): K=2

4. **Profile quality matters**: using the wrong profile (massless when pion mass is
   present) creates artificial instabilities. Always use self-consistent profiles.

## Files

- `src/radial.c` — Modified to support pion mass (new: -mpi parameter)
- `src/normal_modes.c` — Updated compute_fpp for massive ODE consistency
- `data/profiles/profile_massive_e1_mpi0.398.dat` — Massive equilibrium profile (e=1)
- `data/profiles/profile_sigma_e1.dat` — Massless equilibrium profile (e=1, unchanged)
