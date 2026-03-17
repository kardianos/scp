# V26-DynAC: Propagating Massless Braid — Results

## Parameters

| Parameter | Value |
|-----------|-------|
| mu | -20.0 |
| kappa | 20.0 |
| mass | 0.0 (massless) |
| A0 | 0.8 |
| R_tube | 3.0 |
| k | 0.157080 (= 2pi/Lz) |
| N | 128 |
| L | 20.0 (Lz = 40.0) |
| dx | 0.3150 |
| dt | 0.06299 |
| tfinal | 500 |
| BC | periodic z, absorbing x,y |

## Result: NEGATIVE — Immediate Blowup

The propagating massless helical braid does **not** survive. The configuration
undergoes catastrophic instability within the first time step interval, with
energy crashing from +59.21 to -6096 by t=16.6.

### Timeline

| t | E_total | E_kin | E_grad | E_pot | fc | peak |phi| | |P|_max |
|---|---------|-------|--------|-------|----|-------------|---------|
| 0.0 | 59.21 | 13.50 | 73.79 | -28.08 | 0.993 | 0.798 | 0.127 |
| 16.6 | -6095.71 | 1514.02 | 5377.08 | -12986.81 | 0.000 | 1.934 | 2.106 |
| 33.3 | -7076.39 | 875.42 | 3891.92 | -11843.73 | 0.000 | 4.385 | 3.900 |
| 100 | -7623.80 | 607.94 | 3782.34 | -12014.08 | 0.000 | 1.898 | 1.871 |
| 250 | -7828.96 | 454.15 | 3257.35 | -11540.45 | 0.000 | 1.966 | 1.933 |
| 500 | -7915.41 | 410.76 | 3328.60 | -11654.77 | 0.000 | 1.540 | 1.742 |

### Failure Mechanism

1. **Instant collapse of core fraction**: fc drops from 0.993 to 0.000 by t=16.6.
   The transverse Gaussian envelope offers no confinement against the attractive
   triple-product potential (mu < 0).

2. **Energy runaway**: With m=0, there is no restoring force. The attractive
   potential V = mu P^2/(1 + kappa P^2) drives field amplitudes up from
   A0=0.8 to peak values of 1.5--4.4, feeding energy into the potential well.
   Total energy drops from +59 to -7900 as E_pot reaches -11,500.

3. **Turbulent saturation**: After the initial blowup (t < 50), the system
   reaches a turbulent quasi-steady state with energy slowly drifting downward.
   Peak amplitudes stabilize around 1.4--1.5 (late time) and energy density
   is spread uniformly across the grid (fc=0, Rp=0).

4. **No propagation signal**: The z-momentum Pz fluctuates around zero
   (+/-60) with no systematic drift, confirming that the helical propagation
   is immediately destroyed.

### Phase 2: DFT Analysis

- rho(center) mean: -0.412, variance: 8.3e-4, relative variance: 4.9e-3
- Peak omega (rho): 0.14 (low power)
- Peak omega (phi_0): 0.38 (moderate power)
- Breathing: NO (dispersed turbulent state, not breathing)

The center-point energy density is negative on average (attractive potential
dominates), and fluctuations are small relative to the mean — the system is
in a thermalized turbulent state, not a coherent oscillation.

## Interpretation

The massless propagating braid fails for the same fundamental reason as other
massless configurations: without a mass term (m=0), there is no mechanism to
prevent the attractive triple-product potential from driving runaway growth.
The kappa saturation term (kappa P^2 in the denominator) prevents literal
divergence, but the system still falls deep into the potential well and
thermalizes across the entire grid.

The key missing ingredient is a **confining term** — either:
- A mass term (m > 0) to provide a restoring force
- A repulsive pairwise coupling (lambda_pw) to resist collapse
- Sufficient kinetic pressure from high-k modes

At k = 2pi/40 = 0.157, the kinetic energy (Ek=13.5) is far too small compared
to the available potential energy depth (~12,000), so the wave cannot resist
collapse.

## Output Files

- `data/dynAC_phase1.tsv` — full time series (5000 snapshots)
- `data/dynAC_dft.tsv` — DFT power spectra
- `data/dynAC_rho_history.tsv` — center-point time history
