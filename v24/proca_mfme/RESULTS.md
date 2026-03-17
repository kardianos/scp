# V24-P5: Combined MF+ME Results

## Setup

Five fields: phi_1, phi_2, phi_3 (massive scalars, m=1), theta (massless Goldstone).
Parameters: mu=-20, kappa=20, A=0.8, sigma=3, Nx=8000, xmax=200, tfinal=10000.

EOM combines triple-product coupling, pairwise coupling lambda, and theta backreaction:
- phi_a: Lapl - m^2 phi_a - lambda(phi_b + phi_c) - dV_triple/dphi_a + 2g(dt_theta)P(dP/dphi_a)/(1+kP^2)^2
- theta: Lapl_theta - g * d_t[ P^2/(1+kP^2)^2 ]

Mass spectrum: m^2_A = m^2 - lambda (antisymmetric), m^2_S = m^2 + 2*lambda (symmetric).

## Phase 1: Baselines

| Run | g | lambda | omega | fc | E_phi | E_theta | survived |
|-----|---|--------|-------|----|-------|---------|----------|
| MF only | 1.0 | 0 | 0.888 | 0.999 | 1.222 | 5.9e-4 | YES |
| ME only | 0 | 0.99 | 1.602 | 0.999 | 11.718 | 0 | YES |

MF baseline: omega=0.888 < m=1.0 (below mass gap, true oscillon). Theta backreaction
drains negligible energy (E_theta ~ 6e-4 vs E_phi ~ 1.2).

ME baseline: omega=1.602 (this is the SYMMETRIC mode frequency, above gap).
m_A = sqrt(m^2 - lambda) = 0.100, predicted Proca range = 10.0.

## Phase 2: Combined (g=1, lambda=0.99)

| omega | fc | E_phi | E_theta | peak_phi | survived |
|-------|----|-------|---------|----------|----------|
| 1.602 | 0.998 | 12.088 | 2.2e-4 | 0.519 | YES |

Combined run is dominated by the symmetric mode (omega=1.602, same as ME-only).
The theta energy is small (2.2e-4) — backreaction is a perturbation, not dominant.
E_phi = 12.088 vs 11.718 (ME-only) — 3.2% higher with theta backreaction.

## Phase 3: Push lambda past m^2 = 1.0

### With backreaction (g=1)

| lambda | m^2_anti | m_A | omega | fc | E_phi | E_theta | survived |
|--------|----------|-----|-------|----|-------|---------|----------|
| 0.990 | +0.010 | 0.100 | 1.602 | 0.998 | 12.088 | 2.2e-4 | YES |
| 0.995 | +0.005 | 0.071 | 1.602 | 0.998 | 12.160 | 2.7e-4 | YES |
| 1.000 | 0.000 | 0.000 | 1.608 | 0.998 | 12.246 | 1.9e-4 | YES |
| 1.010 | -0.010 | tachy | 1.614 | 0.998 | 12.459 | 1.3e-4 | YES |
| 1.050 | -0.050 | tachy | 1.638 | 1.000 | 14.276 | 1.2e-3 | YES |

### Without backreaction (g=0, control)

| lambda | m^2_anti | omega | fc | E_phi | survived |
|--------|----------|-------|----|-------|----------|
| 0.990 | +0.010 | 1.602 | 0.999 | 11.718 | YES |
| 0.995 | +0.005 | 1.602 | 0.999 | 11.761 | YES |
| 1.000 | 0.000 | 1.608 | 0.999 | 11.806 | YES |
| 1.010 | -0.010 | 1.614 | 0.999 | 11.901 | YES |
| 1.050 | -0.050 | 1.638 | 0.998 | 12.478 | YES |

## Key Findings

### 1. Oscillon survives past the tachyonic boundary at ALL tested lambda values

Both g=1 (with theta) and g=0 (without theta) oscillons survive at lambda=1.01
and lambda=1.05, where m^2_anti = m^2 - lambda < 0 (tachyonic antisymmetric mode).

**Why?** The oscillon sits in the SYMMETRIC sector (all three phi equal). The
symmetric mode mass is m^2_S = m^2 + 2*lambda, which INCREASES with lambda.
The tachyonic instability is in the antisymmetric sector, which is not excited
by the symmetric initial condition. The oscillon naturally avoids the instability.

### 2. Theta backreaction effect is SMALL but detectable

E_phi is systematically higher with g=1 than g=0:
- lambda=0.99: 12.088 vs 11.718 (+3.2%)
- lambda=1.00: 12.246 vs 11.806 (+3.7%)
- lambda=1.05: 14.276 vs 12.478 (+14.4%)

The effect grows with lambda. At lambda=1.05, the theta backreaction injects
14% more energy into the phi sector. This is consistent with the positive
feedback loop hypothesized in the proposal.

### 3. At lambda=1.05 (tachyonic), E_phi is GROWING over time

With g=1 at lambda=1.05: E_phi goes from ~14.0 at t=0 to 14.276 at t=10000.
This is INCREASING — the oscillon is gaining energy, not losing it. The theta
field recycles energy back into the oscillon core through the backreaction term.
fc reaches 1.000 (perfect localization).

Without theta (g=0) at same lambda: E_phi goes from 13.6 to 12.478 (decreasing).

### 4. Omega stays above mass gap — NOT a true sub-gap oscillon

All runs show omega ~ 1.6 (above m=1.0 mass gap). This is the symmetric mode
frequency omega_S ~ sqrt(m^2 + 2*lambda) ~ 1.73 for lambda=1.0. The peak
frequency 1.6 is close to this.

The theta backreaction does NOT push omega below the gap (as originally hypothesized).
The dominant oscillation is the symmetric mode which cannot radiate efficiently
because it's coherent across all three fields.

### 5. The positive feedback loop EXISTS but is weak

At lambda=1.05 with g=1: E_theta = 1.2e-3 (highest of all runs), and E_phi
is increasing. The theta field is sourced by the oscillating P^2, and its
backreaction pumps energy back. But the effect is order ~10^-4 relative to
E_phi — the feedback loop operates but doesn't produce dramatic amplification.

## Conclusions

1. **The oscillon survives at lambda > m^2** — but this is NOT due to theta
   backreaction. It survives because the symmetric mode is inherently stable;
   only antisymmetric perturbations would trigger the tachyonic instability.

2. **Theta backreaction provides measurable energy injection** (3-14% excess
   E_phi at t=10000), growing with lambda. At lambda=1.05, E_phi is genuinely
   increasing over time when theta is active.

3. **The feedback loop is real but perturbative** — E_theta/E_phi ~ 10^-4.
   The Goldstone channel does not dramatically extend the stability range
   because the oscillon was already stable in the symmetric sector without it.

4. **The original hypothesis (theta pushes omega deeper below gap) is NOT
   confirmed.** Omega stays at ~1.6 > m for all runs. The oscillon at large
   lambda is a symmetric-mode configuration, not a sub-gap state.

5. **For a genuine test of theta-assisted stability**, one would need to
   excite an antisymmetric perturbation at lambda > m^2 and check whether
   theta backreaction suppresses the tachyonic growth. This was not tested here.
