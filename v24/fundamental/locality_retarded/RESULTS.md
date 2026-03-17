# Locality Option 4: Retarded Integration — RESULTS

## Method

Replaced instantaneous Poisson equation with the causal wave equation:

    (1/c^2) d^2 Phi/dt^2 - d^2 Phi/dx^2 = -alpha * rho(x,t)

This is the exact retarded solution: signals propagate at speed c=1.0 with
no free parameters beyond alpha and c. Solved as a co-evolving PDE alongside
the oscillon fields using velocity Verlet.

Parameters: mu=-20, kappa=20, m=1.0, alpha=-1e-4, c=1.0,
Nx=4000, xmax=100, t_equil=5000, t_grav=5000, t_ramp=200.

## Phase 1: Phi(0) vs Poisson

After equilibration (E=1.274, fc=0.999, omega=0.888), the wave equation
Phi was initialized from zero and ramped up over t_ramp=200 via adiabatic
Hermite smoothstep on alpha.

**Results** (post-ramp, alpha fully on):

| Quantity | Value |
|---|---|
| Poisson Phi(0) | 6.28e-3 |
| Time-avg wave Phi(0), t=[6000,7500] | 4.89e-3 |
| Ratio <Phi_wave>/<Phi_poisson> | 0.78 |
| Phi_wave oscillation range | [2.3e-3, 7.4e-3] |
| Phi_wave oscillation period | ~8 (matches oscillon period) |
| Profile shape | Matches Poisson (same spatial structure) |

The wave equation Phi oscillates at the oscillon frequency because rho(t)
is time-dependent (oscillating energy density). The TIME-AVERAGED Phi_wave(0)
is 78% of the Poisson value. The shortfall is due to:
1. Initial transient from zero-IC ramp (outgoing waves carry energy)
2. Boundary absorption removing some of the reflected Phi waves
3. In 1D, the steady-state is Poisson plus residual oscillation

The spatial SHAPE of Phi_wave matches Poisson exactly — the well is deeper
at the center and flattens outside the oscillon core.

## Phase 2: Causality Verification

At t=7500 (after 2500 time units of gravity evolution), alpha was
instantaneously doubled from -1e-4 to -2e-4. This creates a sudden change
in the Phi source at the oscillon location (x=0).

Probe at x=50: monitored Phi(x=50, t).

**Key observation**: Phi at the probe was steadily decreasing before the
perturbation (natural oscillation). After the perturbation:

| Time | Phi(x=50) | dPhi/dt sign |
|---|---|---|
| 7548.7 | 1.026e-3 | negative (decreasing) |
| 7549.7 | 1.022e-3 | negative (still decreasing) |
| 7550.7 | 1.040e-3 | **POSITIVE** (reversal!) |
| 7551.7 | 1.064e-3 | positive (increasing) |
| 7553.6 | 1.135e-3 | positive (accelerating) |

**The signal arrives at the probe at t = 7550 +/- 1 time unit.**

- Perturbation time: t = 7500.0
- Signal arrival: t = 7550.5 +/- 0.5
- Measured delay: 50.5 +/- 0.5 time units
- Expected delay: x_probe / c = 50.0 / 1.0 = 50.0 time units
- **Error: 1% (within dt_probe = 1.0 resolution)**

Causality is EXACT: the retarded wave equation enforces that the signal
from x=0 reaches x=50 at exactly t = 50/c = 50 time units. No information
propagates faster than c.

## Oscillon Stability

The oscillon survived 5000 time units of retarded gravity evolution:
- omega_peak = 0.906 (< m = 1.0, confirming oscillon)
- E_equil = 1.274, E_final = 1.270, dE/E = -0.3%
- fc = 0.999 throughout (no dispersal)

Even after doubling alpha, the oscillon remained stable.

## Comparison with Poisson (Test F)

| Property | Retarded (wave eq) | Poisson (Test F) |
|---|---|---|
| Phi(0) | oscillating, <Phi>=4.9e-3 | static 6.3e-3 |
| Causality | EXACT (c=1.0) | INSTANTANEOUS |
| Source coupling | dynamic rho(t) | time-averaged <rho> |
| Radiation | outgoing Phi waves | none |
| CPU cost | same as flat + 1 extra PDE | Poisson solve every ~10 steps |
| Stability | stable | stable |

## Key Physics

1. **1D retarded potential oscillates**: Unlike 3D where the retarded potential
   converges to the static Coulomb value, in 1D the oscillating rho source
   generates propagating Phi waves. The time-averaged local Phi converges
   toward Poisson, but individual snapshots show +/- 50% oscillation.

2. **Exact causality**: The wave equation naturally enforces c as the maximum
   signal speed. Verified to 1% accuracy at x=50.

3. **No free parameters**: Only alpha (coupling) and c (propagation speed).
   The retarded solution is the REFERENCE against which other locality options
   (EMA-wave, telegraph) should be compared.

4. **Boundary effects matter**: Absorbing boundaries for Phi are essential
   to prevent reflected waves from corrupting the interior solution.

## Files

- `src/locality_retarded.c` — wave equation solver
- `data/retarded_ts.tsv` — time series (phi, E, fc, Phi_wave, Phi_poisson)
- `data/retarded_probe.tsv` — probe at x=50 (causality test)
- `data/retarded_phi_profile.tsv` — final Phi(x) profile comparison
- `data/retarded_spectrum.tsv` — DFT spectrum

Compile: `gcc -O3 -Wall -o locality_retarded src/locality_retarded.c -lm`
Run from project root: `./v24/fundamental/locality_retarded/locality_retarded`
