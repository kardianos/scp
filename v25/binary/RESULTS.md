# V25 Phase 5: Binary Oscillon GW Emission — Results

## Setup

Two oscillons placed at x = +/-D/2 with tangential velocity v_orb in the y-direction.
Grid: N=128, L=30 (dx=0.47). Measurement sphere at R=20 (inside undamped zone R_abs=21).
All couplings from v25.c: eta=0.1, lambda_L=0.1, lambda_pw=0.5, mu=-20, kappa=20, alpha_g=0.001.

Equilibration: single oscillon at N=96, L=15 for t=200, then copied to two positions.

Two runs: v_orb = {0.01, 0.10}, both with D=12.

## Trajectory Results

Both runs show nearly identical dynamics — the orbital velocity is far too small to resist
the inter-oscillon attraction:

| Phase         | v=0.01          | v=0.10          |
|---------------|-----------------|-----------------|
| t=0 sep       | 11.81           | 11.81           |
| First bounce  | t~75, sep~3.6   | t~75, sep~3.6   |
| Max rebound   | t~110, sep~8.5  | t~110, sep~8.5  |
| Late sep      | ~4.0            | ~4.0            |
| E(t=0)        | 456.5           | 456.6           |
| E(t=500)      | ~299            | ~299            |

The oscillons immediately fall toward each other (attraction dominates over centrifugal
barrier), bounce at sep~3.6, rebound to sep~8.5, then settle into damped radial
oscillation around sep~4. They lose ~35% of their energy to radiation by t=500.

**No clean circular orbit forms.** The system is a radially oscillating binary,
not an orbiting one. The y-displacement is tiny (yA < 0.3 at all times).

## GW Strain Results

### v_orb = 0.10

Peak |h+| and |hx| at R=20 (6 angles):

| Direction              | |h+|_max    | |hx|_max    |
|------------------------|-------------|-------------|
| pole (th=0)            | 1.738e-3    | 1.897e-3    |
| equator x (th=pi/2,p=0) | 8.907e-4  | 2.487e-3    |
| equator y (th=pi/2,p=pi/2) | 1.497e-3 | 1.497e-3  |
| 45 deg (th=pi/4,p=0)  | 1.916e-3    | 3.134e-3    |
| south pole (th=pi)     | 1.600e-3    | 1.966e-3    |
| equator 45 (th=pi/2,p=pi/4) | 1.962e-3 | 3.311e-3 |

### v_orb = 0.01

Nearly identical to v=0.10 (confirming that orbital motion is negligible):

| Direction              | |h+|_max    | |hx|_max    |
|------------------------|-------------|-------------|
| pole (th=0)            | 1.725e-3    | 1.882e-3    |
| equator x (th=pi/2,p=0) | 9.047e-5  | 2.575e-3    |
| equator y (th=pi/2,p=pi/2) | 1.696e-3 | 1.787e-3  |
| 45 deg (th=pi/4,p=0)  | 2.002e-3    | 3.140e-3    |
| south pole (th=pi)     | 1.607e-3    | 1.897e-3    |
| equator 45 (th=pi/2,p=pi/4) | 2.013e-3 | 3.158e-3 |

## Angular Pattern Analysis

**SPIN-2 TEST: NEGATIVE.**

- h+ at pole / h+ at equator = 1.16 (v=0.1), 1.02 (v=0.01)
  Expected for spin-2: << 1 (suppressed at poles). Observed: ~1 (isotropic).
- h+ at 45 deg / h+ at equator = 1.28 (v=0.1), 1.18 (v=0.01)
  Expected for spin-2: ~0.5. Observed: >1.

The radiation pattern is NOT quadrupolar. It is approximately isotropic (monopolar),
consistent with spherically symmetric oscillon breathing dominating the signal.

## GW Frequency

- DFT peak: omega_GW = 0.193 (both runs)
- Expected 2*Omega: 0.033 (v=0.1) or 0.003 (v=0.01)
- Ratio: omega_GW / (2*Omega) = 5.8 (v=0.1) or 58.0 (v=0.01)

**The GW frequency does NOT match 2*Omega.** The measured frequency (0.193)
matches the oscillon's internal breathing frequency (omega/m ~ 0.193, consistent
with Phase 1 results). The strain signal is dominated by the oscillons' internal
dynamics, not their orbital motion.

## Multipole Decomposition

At R=20, decomposing |sigma^TT|^2 into Legendre polynomials:

| l | Fraction (v=0.1) | Fraction (v=0.01) |
|---|-------------------|---------------------|
| 0 | 93-97%            | 93-97%              |
| 1 | 0.5-2%            | 0.5-2%              |
| 2 | 2-7%              | 2-7%                |

**l=0 (monopole) dominates at 93-97%.** This is the same result as Phase 2 for
a single oscillon — the signal is overwhelmingly spin-0.

## Physical Interpretation

1. **Attraction too strong**: At D=12, the inter-oscillon attraction pulls the pair
   together in ~70 time units, regardless of v_orb. No stable orbit forms.

2. **Breathing dominates orbital**: The oscillon breathing frequency (omega~0.19)
   is 10-100x higher than the orbital frequency (Omega~0.01-0.002). The breathing
   mode radiates much more strongly than the orbital motion.

3. **Monopole radiation from breathing**: Each oscillon is nearly spherical and
   oscillates radially. This produces spin-0 (compression wave) radiation, not
   spin-2 (shear wave) radiation. The binary's orbital quadrupole moment is tiny
   compared to the breathing signal.

4. **Energy loss**: The system loses ~35% of energy by t=500 (456 -> 299), mainly
   from the initial collision/bounce and subsequent breathing radiation. This is
   NOT gravitational wave energy loss from orbital inspiral — it is acoustic
   radiation from the collision.

## What Would Be Needed for Spin-2 Detection

To observe quadrupolar GW at 2*Omega from a binary:

1. **Stable circular orbit**: Need much larger D (>30?) or a repulsive core to
   prevent merger. Alternatively, need to find the exact v_orb that balances
   the attraction (requires knowing F(D) precisely).

2. **Suppress breathing mode**: The breathing mode swamps the orbital signal.
   Would need oscillons that don't breathe (true solitons) or a way to subtract
   the monopole background.

3. **Longer evolution**: Even at v=0.1, D=12: T_orb=377. Would need t>1000 for
   multiple orbits, requiring a larger grid to avoid boundary effects.

4. **Stronger metric coupling**: alpha_g=0.001 is weak. The strain signal is
   ~1e-3, barely above the spherical oscillon background. Larger alpha_g might
   amplify the quadrupolar component but also destabilize the oscillons.

## Files

- `src/binary.c` — binary oscillon GW simulator
- `data/binary_v0.100_ts.tsv` — trajectory (positions, separation, energy) for v=0.1
- `data/binary_v0.100_gw.tsv` — TT strain at 6 angles vs time for v=0.1
- `data/binary_v0.100_multipoles.tsv` — multipole decomposition vs time for v=0.1
- `data/binary_v0.010_ts.tsv` — same for v=0.01
- `data/binary_v0.010_gw.tsv` — same for v=0.01
- `data/binary_v0.010_multipoles.tsv` — same for v=0.01

## Conclusion

**The binary oscillon system does NOT produce spin-2 gravitational waves at 2*Omega.**
The signal is dominated by monopolar (spin-0) radiation from the oscillons' internal
breathing mode at omega~0.19. No stable orbit forms at D=12 — the pair merges
within ~70 time units. The angular pattern of h+ is isotropic, not quadrupolar.

This is a physically correct null result: elastic oscillons that breathe spherically
will always emit predominantly monopole radiation. Quadrupolar emission requires
a sustained aspherical mass distribution that changes at the orbital frequency,
which requires either (a) a stable orbit or (b) non-breathing compact objects.
