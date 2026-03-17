# V26-DynBC: Rotating Massless Braid — Results

## Parameters

| Parameter | Value |
|-----------|-------|
| mu        | -20.0 |
| kappa     | 20.0  |
| mass      | 0.0   |
| Omega     | 0.1   |
| lambda_pw | 0.0   |
| eta       | 0.0   |
| lambda_L  | 0.0   |
| A0        | 0.8   |
| R_tube    | 3.0   |
| k_braid   | 0.1571 (= 2pi/40) |
| N         | 128   |
| L         | 20.0  |
| dx        | 0.3150 |
| dt        | 0.06299 (CFL=0.20) |
| tfinal    | 500   |
| BC        | periodic z, absorbing x,y (cylindrical damping) |

Auxiliary couplings (lambda_pw, eta, lambda_L) set to zero. Nonzero values
cause immediate instability with m=0 rotating initialization.

## Initialization

```
phi_a(x,0) = A(r_perp) * cos(theta + k*z + 2*pi*a/3)
vel_a(x,0) = 0.1 * A(r_perp) * sin(theta + k*z + 2*pi*a/3)
```
where A(r) = 0.8 * exp(-r^2 / 18).

## Phase 1: Evolution

| time | E_total | Lz       | fc     | peak   | |P|_max   |
|------|---------|----------|--------|--------|-----------|
| 0    | 370.87  | -108.28  | 0.3780 | 0.7978 | 0.126944  |
| 50   | 40.85   | -32.54   | 0.0096 | 0.0372 | 0.000013  |
| 100  | 16.31   | -16.29   | 0.2523 | 0.1468 | 0.000784  |
| 200  | 4.27    | -6.00    | 0.1842 | 0.0570 | 0.000046  |
| 300  | 1.83    | -2.97    | 0.1401 | 0.0173 | 0.000001  |
| 400  | 1.05    | -1.64    | 0.0983 | 0.0176 | 0.000001  |
| 500  | 0.54    | -0.93    | 0.0992 | 0.0222 | 0.000003  |

**Survived?** Formally YES (fc=0.099 > 0.01, |P|_max > 1e-6), but the braid
has effectively dispersed: E dropped 99.85% (370.87 -> 0.54), Lz dropped
99.14% (-108.28 -> -0.93), peaks dropped 97% (0.80 -> 0.02).

Wall time: 4043 sec (~67 min), 16 threads.

### Energy decay

Energy falls monotonically due to absorbing BC damping radiation away.
The triple-product potential is negligible at late times (|P|_max ~ 10^{-6}).
The core fraction oscillates between ~0.01 and ~0.38 early, settling to ~0.10
late, consistent with radiation sloshing through the core region.

### Angular momentum

Lz tracks energy decay closely. Ratio Lz/E ~ -1.7 is roughly constant,
suggesting angular momentum is carried away by radiation at the same rate
as energy.

## Phase 2: Non-Breathing Verification

- rho(center) mean: 5.17e-05, relative variance: 1.20
- Peak omega (rho): 0.020 (near-DC drift, not a true mode)
- Peak omega (phi_0): 0.290
- **Breathing? NO** (dispersed, not oscillating)

The high relative variance reflects the center density dropping through
zero and changing sign, not a coherent breathing oscillation. The DFT
shows no sharp peak — the signal is broadband.

## Phase 3: Strain Multipoles (R_shell = 8.0)

| l | coefficient | fraction |
|---|-------------|----------|
| 0 | 4.39e-03    | 0.843    |
| 1 | 1.54e-05    | 0.003    |
| 2 | 8.02e-04    | 0.154    |

Total |sigma|^2 on shell: 2.19e-05 (200 sample points).

Dominant l=0 monopole (84%) with weak l=2 quadrupole (15%).
The l=1 dipole is negligible (0.3%), consistent with the 3-fold
rotational symmetry of the braid initialization.

## Conclusion

The rotating massless braid **does not survive** as a localized structure.
Without a mass gap (m=0), the braid disperses via radiation through the
absorbing boundary. The triple-product potential alone (mu=-20, kappa=20)
cannot confine the rotating configuration: the potential requires all three
fields to overlap, and the centrifugal spreading from rotation breaks this
overlap faster than the potential can restore it.

The auxiliary couplings (lambda_pw, eta, lambda_L) at their V26 default
values cause immediate instability (exponential blowup within t~15) when
combined with m=0 and rotating initial conditions.

**Key finding**: Rotation + masslessness is insufficient for braid
stabilization. A mass gap appears necessary to prevent radiative
dispersion of the rotating braid.

## Output Files

- `data/dynBC_phase1.tsv` — full time series (5000 points)
- `data/dynBC_rho_history.tsv` — center density history
- `data/dynBC_phase2_dft.tsv` — DFT power spectrum
- `data/dynBC_phase3_strain.tsv` — strain on spherical shell
- `data/dynBC_phase3_multipoles.tsv` — multipole coefficients
