# Phonon Test: Radial Depletion Profile of a Single Braid

## Simulation Parameters

| Parameter | Value |
|-----------|-------|
| Grid | N=256, L=60 (box: [-60,60]^3) |
| Resolution | dx=0.4706, dt=0.05647 |
| Physics | m^2=2.25 (m=1.5), mu=-41.345, kappa=50 |
| Background | A_bg=0.1 (standing wave seed) |
| Runtime | T=200 (3541 steps, 23.9 min on 8 threads) |
| Energy conservation | drift = -0.061% over T=200 |
| Memory | 1.21 GB (single allocation) |
| Snapshots | t=0, t=100, t=200 |

Yukawa range for mass m=1.5: R_Yukawa = 1/m = 0.67 code units.

## Braid Center

| Time | x_c | y_c | z_c | Note |
|------|-----|-----|-----|------|
| t=100 | -0.000 | -0.000 | 10.631 | Braid oscillates in z |
| t=200 | -0.000 | 0.000 | 6.487 | Standing wave envelope |

The z-offset oscillation (period ~250 time units) indicates the braid is a
standing wave in z, not fully settled. The x-y center is fixed at origin.

## Radial Depletion Profile (Spherical Shells, dr=0.5)

Measurement: rho(r) = <sum_a phi_a^2> averaged in spherical shells centered on
the energy-weighted braid centroid.

### t=100 Profile

Background rho (r>40): 1.794e-02

| r | rho(r) | delta_rho | counts |
|-----|---------|-----------|--------|
| 0.25 | 7.670e-02 | +5.875e-02 | 4 |
| 1.25 | 7.629e-02 | +5.835e-02 | 100 |
| 2.25 | 7.826e-02 | +6.031e-02 | 308 |
| 3.25 | 8.423e-02 | +6.628e-02 | 624 |
| 4.25 | 8.310e-02 | +6.515e-02 | 1072 |
| 5.25 | 7.145e-02 | +5.351e-02 | 1628 |
| 6.25 | 5.860e-02 | +4.065e-02 | 2340 |
| 7.25 | 4.841e-02 | +3.047e-02 | 3160 |
| 8.25 | 4.590e-02 | +2.796e-02 | 4100 |
| 9.25 | 4.288e-02 | +2.493e-02 | 5108 |
| 10.25 | 3.931e-02 | +2.137e-02 | 6288 |
| 12.25 | 3.612e-02 | +1.818e-02 | 8984 |
| 14.25 | 3.401e-02 | +1.607e-02 | 12216 |
| 16.25 | 3.290e-02 | +1.496e-02 | 15968 |
| 18.25 | 3.331e-02 | +1.537e-02 | 19992 |
| 20.25 | 2.989e-02 | +1.195e-02 | 24664 |
| 22.25 | 2.621e-02 | +8.264e-03 | 29980 |
| 24.25 | 2.333e-02 | +5.384e-03 | 35320 |
| 26.25 | 2.158e-02 | +3.635e-03 | 41672 |
| 28.25 | 1.980e-02 | +1.861e-03 | 48152 |
| 30.25 | 2.033e-02 | +2.388e-03 | 55160 |
| 32.25 | 2.363e-02 | +5.683e-03 | 62712 |
| 34.25 | 2.412e-02 | +6.178e-03 | 70656 |
| 36.25 | 2.019e-02 | +2.243e-03 | 79160 |
| 38.25 | 1.733e-02 | -6.132e-04 | 88436 |
| 40.25 | 1.800e-02 | +5.845e-05 | 97760 |

### t=200 Profile

Background rho (r>40): 1.762e-02

| r | rho(r) | delta_rho | counts |
|-----|---------|-----------|--------|
| 0.25 | 2.733e-01 | +2.557e-01 | 8 |
| 1.25 | 2.970e-01 | +2.794e-01 | 108 |
| 2.25 | 2.292e-01 | +2.116e-01 | 312 |
| 3.25 | 1.087e-01 | +9.103e-02 | 656 |
| 4.25 | 4.895e-02 | +3.133e-02 | 1116 |
| 5.25 | 4.545e-02 | +2.783e-02 | 1668 |
| 7.25 | 6.488e-02 | +4.726e-02 | 3200 |
| 10.25 | 3.893e-02 | +2.130e-02 | 6320 |
| 12.25 | 1.913e-02 | +1.501e-03 | 9100 |
| 15.25 | 2.107e-02 | +3.445e-03 | 14128 |
| 20.25 | 3.218e-02 | +1.456e-02 | 24696 |
| 25.25 | 2.970e-02 | +1.207e-02 | 38456 |
| 30.25 | 2.064e-02 | +3.019e-03 | 55232 |
| 35.25 | 2.152e-02 | +3.892e-03 | 74896 |
| 40.25 | 2.024e-02 | +2.612e-03 | 97620 |

## Fit Results

### Fitting over r = 5 to 40

#### t=100 (cleaner, less boundary contamination)

| Model | R^2 | Parameters |
|-------|-----|------------|
| Yukawa (m=1.5 fixed) | **-0.44** | A=9.64e+02 |
| Yukawa (m free) | 0.952 | A=0.299, **m=0.025** (range=40) |
| Power law | 0.944 | B=0.451, **n=1.30** |
| Two-component (m=1.5 fixed) | 0.945 | A=30.9, B=0.428, n=1.29 |
| Two-component (m free) | **0.966** | A=-3.22, m=0.231, B=10.3, **n=2.27** |

#### t=200 (radiation shell contamination)

| Model | R^2 | Parameters |
|-------|-----|------------|
| Yukawa (m=1.5 fixed) | **-0.64** | A=6.45e+02 |
| Yukawa (m free) | 0.616 | A=0.265, **m=0.019** (range=53) |
| Power law | 0.617 | B=0.372, **n=1.24** |
| Two-component (m free) | **0.699** | A=-255, m=1.31, B=1.41, **n=1.77** |

### Fitting over restricted ranges (t=100)

| Range | Power law n | Power law R^2 | Yukawa(free) m | Yukawa(free) R^2 |
|-------|-------------|---------------|----------------|------------------|
| r=5-12 (core halo) | **1.37** | **0.986** | 0.046 | 0.979 |
| r=5-15 (near-field) | **1.27** | **0.982** | 0.029 | 0.974 |
| r=5-25 (inner half) | **1.16** | 0.963 | 0.013 | 0.959 |
| r=10-30 (mid-range) | 1.40 | 0.812 | 0.028 | 0.837 |
| r=20-40 (far-field) | 3.07 | 0.619 | 0.075 | 0.602 |

### Model-Independent Log-Log Slopes (t=100)

| r range | log-log slope | Implied n |
|---------|---------------|-----------|
| r=5-15 | -1.169 | **1.17** |
| r=10-25 | -1.346 | **1.35** |
| r=5-25 | -1.189 | **1.19** |

Point-by-point log-log slopes at t=100:

| r | slope | delta_rho |
|---|-------|-----------|
| 5.25 | -1.22 | 5.35e-02 |
| 6.25 | -1.74 | 4.07e-02 |
| 7.25 | -1.35 | 3.05e-02 |
| 8.25 | -0.82 | 2.80e-02 |
| 9.25 | -1.24 | 2.49e-02 |
| 10.25 | -1.21 | 2.14e-02 |
| 11.25 | -0.91 | 1.97e-02 |
| 12.25 | -0.70 | 1.82e-02 |
| 14.25 | -1.01 | 1.61e-02 |
| 16.25 | +0.19 | 1.50e-02 |
| 18.25 | -1.24 | 1.54e-02 |
| 20.25 | -2.99 | 1.19e-02 |
| 22.25 | -4.67 | 8.26e-03 |
| 24.25 | -4.58 | 5.38e-03 |
| 26.25 | -9.03 | 3.63e-03 |
| 28.25 | -3.55 | 1.86e-03 |

The slope steepens dramatically beyond r~20, from ~-1.2 to ~-5 and worse.
This steepening coincides with the radiation shell (visible as a bump at
r~30-35 at t=100). The slope changes sign at r>29, confirming the radiation shell.

## Analysis and Conclusions

### 1. Yukawa (m=1.5) is categorically excluded

Yukawa with the physical mass m=1.5 gives **negative R^2** in all fits (-0.44 to
-0.64), meaning it is worse than a flat line. At r=10, Yukawa m=1.5 predicts
delta_rho ~ e^{-15}/10 ~ 3e-8. The measured value is 2e-2, which is **500,000x
larger**. The field excess extends to r>40, whereas Yukawa m=1.5 would be
negligible beyond r~5.

### 2. The profile is NOT Yukawa at any mass

When the Yukawa mass is left free, the fit converges to m=0.019-0.025, giving a
Yukawa range of 40-53 code units. This is effectively m=0, meaning the fit is
really finding a power law (exp(-0.02r)/r ~ 1/r for r<40). The free-m Yukawa
and pure power law give nearly identical R^2.

### 3. Power-law exponent: n = 1.2 +/- 0.2

Across all fitting ranges and both time snapshots:
- Core halo (r=5-12): n = 1.37
- Near-field (r=5-15): n = 1.27
- Full range (r=5-40): n = 1.24-1.30
- Log-log fits: n = 1.17-1.35

**Best estimate: n = 1.2 +/- 0.2**

This is between 1/r (Coulomb/massless scalar potential, n=1) and 1/r^2
(gravitational force / radiation field, n=2).

### 4. Radiation shell complicates far-field measurement

The braid is not fully equilibrated even at T=200. A radiation shell propagates
outward:
- At t=100: bump at r=30-35, zero-crossing at r~38
- At t=200: broad excess r=20-40, shoulder structure at r=30-35

The steepening of log-log slopes beyond r~20 (from n~1.2 to n>5) reflects the
radiation shell's leading edge, not the equilibrium profile. The region r=5-20
at t=100 is the most reliable.

### 5. Cross-time evolution

Between t=100 and t=200:
- Core (r<4): phi^2 grows 2-6x (braid still oscillating, not settled)
- Mid-range (r=10-12): drops to 8% of t=100 value (node of standing wave passed)
- Far-field (r=20-30): grows 2-3x (radiation shell advancing)

This confirms the profile is time-dependent and not the equilibrium state.

### 6. Interpretation

**The field phi^2 excess around the braid decays as a POWER LAW with n~1.2,
not as a Yukawa exponential.** Despite m=1.5 being the fundamental field mass,
the braid generates a long-range field profile that is not exponentially
screened.

Possible explanations:
- **Collective/phonon mode**: The braid is a nonlinear bound state whose
  collective excitations may include a massless (or nearly massless) mode
  that propagates the phi^2 disturbance to long range.
- **Radiation halo**: The outgoing radiation from braid settling may establish
  a quasi-static halo with power-law character.
- **Nonlinear screening suppression**: The trilinear coupling V = mu*P^2/(1+kappa*P^2)
  may create an effective mass that vanishes at low amplitude, making the
  tail effectively massless.

The exponent n~1.2 does not match simple predictions:
- n=1 would be a 3D massless scalar (Coulomb-like)
- n=2 would be gravitational or dipole radiation
- n=4 would be the algebraic limit from MEMORY.md (1/r^6 for |omega|^2/r)
- n=1.2 may indicate a mix of monopole (n=1) and higher multipole components

### 7. What would be needed for a definitive answer

1. **Much longer runtime** (T=1000+) to let radiation escape and braid fully settle
2. **Damped evolution** (add friction term) to accelerate settling
3. **Larger box** (L=200) to separate radiation shell from equilibrium profile
4. **Multiple snapshots** (T=500,1000,2000) to verify the profile stabilizes

## Files

- `src/v33_G.c` — simulation source (copied from v34/G_metastability)
- `src/radial_depletion.c` — spherical shell analysis tool
- `fit_depletion.py` — fitting script (Yukawa, power law, two-component)
- `fit_detailed.py` — multi-range detailed fitting script
- `data/phonon/field_t{0000,0100,0200}.bin` — field snapshots (403 MB each)
- `data/phonon/timeseries.tsv` — energy time series
- `data/depletion_t0100.tsv` — radial profile at t=100
- `data/depletion_t0200.tsv` — radial profile at t=200
