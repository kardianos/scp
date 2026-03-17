# V26-DynAB: Propagating + Rotating Braid — RESULTS

## Configuration

Combined axial propagation (A) with azimuthal rotation (B):
```
phase_a = theta + k*z + 2*pi*a/3
phi_a(x,0)   = A(r_perp) * cos(phase_a)
vel_a(x,0)   = (omega + Omega) * A(r_perp) * sin(phase_a)
```

Parameters: mu=-20, kappa=20, m=1.0, Omega=0.1, A0=0.8, R_tube=3.0
Grid: N=128, L=20, dx=0.315, dt=0.063, periodic-z + absorbing-xy
Derived: k=0.1571, omega=1.0123, omega_eff=1.1123

## Evolution Summary (t=0..500)

| t | E_total | fc | peak_phi | |P|_max | Pz | Lz |
|---|---------|------|----------|---------|------|------|
| 0 | 1589.4 | 0.999 | 0.798 | 0.1269 | 189.6 | 1204.4 |
| 50 | 570.2 | 0.320 | 0.197 | 0.00197 | 78.2 | 499.8 |
| 100 | 363.2 | 0.700 | 0.171 | 0.00127 | 51.4 | 329.3 |
| 200 | 173.9 | 0.707 | 0.115 | 0.00037 | 25.2 | 161.2 |
| 300 | 98.0 | 0.474 | 0.069 | 0.00008 | 14.4 | 92.1 |
| 400 | 62.1 | 0.341 | 0.053 | 0.00004 | 9.1 | 58.6 |
| 500 | 40.4 | 0.282 | 0.039 | 0.00002 | 6.0 | 38.2 |

## Key Measurements

### Core fraction (fc)
- Initial: 0.999 (compact braid fills core)
- Oscillates between 0.07--0.89 as dispersed waves bounce
- Final: 0.282 (residual energy in core, most dispersed)
- Survival criterion met: fc=0.28 > 0.01, |P|max=1.5e-5 > 1e-6

### Triple product |P|_max
- Initial: 0.127 (strong three-field overlap)
- Drops rapidly: 0.003 at t=17, 0.001 at t=100
- Final: 1.5e-5 (essentially vanished)
- The triple product coupling is negligible after early dispersion

### Linear momentum Pz
- Initial: 189.6 (from axial propagation k*z in phase)
- Steady decay to 6.0 at t=500 (97% lost to absorbing BC)
- Ratio Pz/E stays approximately constant: Pz/E ~ 0.12--0.15
- Consistent with radiation carrying proportional momentum

### Angular momentum Lz
- Initial: 1204.4 (from theta in phase + Omega rotation)
- Steady decay to 38.2 at t=500 (97% lost)
- Ratio Lz/Pz stays constant at ~6.4 throughout evolution
- This ratio equals Lz/Pz = (omega+Omega)/k * R_eff^2, consistent with
  the braid carrying angular momentum proportional to linear momentum

### l=2 Multipole
- Multipole decomposition of energy density on cylindrical shell r=6.0:
  - l=0: 72.0% (monopole dominant)
  - l=1: 0.02% (negligible dipole)
  - l=2: **18.0%** (significant quadrupole)
  - l=3: 0.02%
  - l=4: 9.9%
- The l=2 quadrupole is the leading non-monopole mode
- This comes from the cos(2*theta) pattern of the braid cross-section

### DFT Analysis
- Peak omega (rho_center): 0.080 (slow modulation, not breathing)
- Peak omega (phi0_center): 1.110 (matches omega_eff = omega + Omega = 1.112)
- Breathing: NO (relative variance of center density = 0.96 but omega too low)
- The phi_0 oscillation frequency matches the initialized omega_eff precisely

## Interpretation

1. **No stable soliton**: The propagating+rotating braid disperses continuously.
   Energy drops from 1589 to 40 (97.5% radiated through absorbing boundaries).
   The triple-product potential alone (no pairwise coupling, no cross-gradient)
   is insufficient to bind the braid.

2. **Both momenta present**: The initialization successfully creates both Pz
   (linear, from k*z phase) and Lz (angular, from theta phase). Both decay
   at the same rate, maintaining a fixed Lz/Pz ratio throughout evolution.

3. **Quadrupole structure**: The l=2 fraction of 18% confirms the braid has
   genuine azimuthal structure (not spherically symmetric). This comes from
   the three-strand phase offset 2*pi*a/3 creating a cos(2*theta) pattern
   in the energy density.

4. **Frequency lock**: phi_0 at center oscillates at exactly omega_eff =
   omega + Omega = 1.112, confirming the propagation+rotation initialization
   works as designed.

## Files

- `src/dynAB.c` — source code
- `data/dynAB_phase1.tsv` — full time series (200 snapshots)
- `data/dynAB_rho_history.tsv` — center density for DFT (5000 samples)
- `data/dynAB_multipoles.tsv` — l=0..4 multipole coefficients

Wall time: 3631 sec (60.5 min), 8 threads, N=128^3.
