# V24-180A Results: 180-degree Anti-Phase Oscillon

## Executive Summary

The 180-degree anti-phase oscillon (phi1=phi2=+A*g, phi3=-A*g) is **exactly
degenerate** with the 0-degree state. Every measurable quantity -- frequency,
energy, radiation rate, profile, core fraction -- is bit-for-bit identical.
This is not approximate; it is an exact Z_2 symmetry of the Lagrangian.

**The 180-degree state is a genuine, stable oscillon, indistinguishable from
the 0-degree state in all dynamical properties.**

## Parameters

    mu = -20.0, kappa = 20.0, m = 1.0
    A = 0.8, sigma = 3.0
    Nx = 4000, xmax = 100, tfinal = 20000

## Phase 1+2: Side-by-Side Comparison

### Breathing Frequency

| State | omega_phi | omega/m | omega_P |
|-------|-----------|---------|---------|
| 0-deg | 0.8820    | 0.8820  | 0.8800  |
| 180-deg | 0.8820  | 0.8820  | 0.8800  |

Both states have omega < m (below mass gap) -- confirmed oscillons.
The P(t) spectrum has the same peak frequency as phi(t).

### Energy

| State | E_init | E(t=10000) | E(t=20000) | dE/dt (2nd half) |
|-------|--------|------------|------------|------------------|
| 0-deg | 3.3652 | 1.2638     | 1.2565     | -7.33e-07        |
| 180-deg | 3.3652 | 1.2638   | 1.2565     | -7.33e-07        |

Energy loss rate is identical: dE/dt = -7.3e-07 (extremely slow radiation).
Lifetime estimate: E/|dE/dt| ~ 1.7 million time units.

### Core Fraction and Peak Amplitude

| State | fc (avg, 2nd half) | peak (avg, 2nd half) |
|-------|--------------------|-----------------------|
| 0-deg | 0.9996            | 0.332                 |
| 180-deg | 0.9996          | 0.332                 |

### Harmonic Spectrum of P(t)

The phi_1(0,t) and P(0,t) spectra are **bit-for-bit identical** between 0-deg
and 180-deg (verified by diff). The only difference is sign:

    0-deg:   P(0,t) = +|f|^3 * oscillation
    180-deg: P(0,t) = -|f|^3 * oscillation

The sign of P does not affect any dynamical quantity because the potential
V = (mu/2) P^2 / (1 + kappa P^2) depends only on P^2.

### Profile Shape

Identical. The fields have the same spatial profile |phi_a(x)| at all times.
The only difference is the sign of phi_3.

## Phase 3: Stability Tests

### 3a. Flip Perturbation (push phi3 toward +)

    Init: phi3 = (-0.8 + 0.08) * g(x) = -0.72 * g(x)
    Result: omega = 0.8820, E_final = 1.2564, dE/dt = -7.14e-07, fc = 0.9996

The 180-degree state survives the flip perturbation completely. The oscillon
relaxes back to the 180-degree configuration. phi_3 remains negative at the
center while phi_1, phi_2 remain positive. The perturbed amplitude is shed as
radiation and the oscillon settles to the same attractor.

### 3b. Amplitude Perturbation (1.2x scaling)

    Init: all amplitudes scaled by 1.2 (0.96 instead of 0.8)
    Result: omega = 0.8910, E_final = 1.2505, dE/dt = -2.09e-06, fc = 0.990

The excess energy (E_init = 5.05 vs 3.37) is radiated away. The oscillon
settles to a slightly larger amplitude state (peak ~ 0.32) with slightly
higher frequency (0.891 vs 0.882). The radiation rate is ~3x higher than
the unperturbed case because the oscillon is still shedding excess energy
at t=20000 (it started with 50% more energy).

### 3c. Asymmetric Perturbation (phi1=0.9A, phi2=1.1A)

    Init: phi1 = 0.72*g, phi2 = 0.88*g, phi3 = -0.80*g
    Result: omega = 0.8820, E_final = 1.2565, dE/dt = -7.15e-07, fc = 0.9996

The asymmetry between phi1 and phi2 is quickly erased. By t~1000, the three
fields have equalized to |phi1| = |phi2| = |phi3|. The oscillon is a
symmetric attractor: any initial asymmetry is radiated away.

## Analysis: Why Exact Degeneracy?

The Z_2 symmetry phi_3 -> -phi_3 is an **exact symmetry** of the EOM:

    EOM for phi_a: d^2 phi_a/dt^2 = nabla^2 phi_a - m^2 phi_a + F_a
    F_a = -mu * P * (dP/dphi_a) / (1 + kappa*P^2)^2

Under phi_3 -> -phi_3:
  - P -> -P
  - dP/dphi_3 = phi_1*phi_2 unchanged
  - F_3 = -mu*(-P)*(phi_1*phi_2) / (1+kappa*P^2)^2 = +mu*P*(phi_1*phi_2)/... = -F_3_original

So -phi_3 satisfies the same EOM with the transformed force, which is just
the original force applied to -phi_3. The mapping (phi_1, phi_2, phi_3) ->
(phi_1, phi_2, -phi_3) maps solutions to solutions.

More generally, flipping any ODD number of fields is a symmetry (since V
depends on P^2, and P changes sign under an odd number of flips). The
six Z_2 orbits are:

    (+,+,+), (+,+,-), (+,-,+), (-,+,+)  [all related by permutation+flip]
    (+,-,-), (-,+,-), (-,-,+)            [two flipped -- same as one flip of complement]
    (-,-,-)                               [same as (+,+,+) since P -> -P]

All give identical dynamics. The "180-degree" label refers to the
relative phase, but the physics is exactly the same.

## Key Answers

1. **omega_180 = omega_0 = 0.882** (exactly). omega/m = 0.882.

2. **E_180 = E_0 = 1.257** at t=20000 (exactly). Same mass.

3. **dE/dt_180 = dE/dt_0 = -7.3e-07** (exactly). Same radiation rate.

4. **YES, the 180-degree state is truly stable.** It survives all three
   perturbation tests (flip, amplitude, asymmetric). It does NOT collapse
   to the 0-degree state -- it remains in the 180-degree configuration.
   The two states are distinct but degenerate minima.

5. **P(t) spectrum is identical** (same |P|^2 at all frequencies). The only
   difference is sign(P), which has no dynamical effect since V depends on P^2.

6. **Profile shape is identical**: |phi_a(x)| is the same for both states.

## Output Files

- `data/char180_0deg_ts.tsv` -- 0-degree time series (20000 records)
- `data/char180_180deg_ts.tsv` -- 180-degree time series
- `data/char180_stab_flip_ts.tsv` -- flip perturbation
- `data/char180_stab_amp_ts.tsv` -- amplitude perturbation
- `data/char180_stab_asym_ts.tsv` -- asymmetric perturbation
- `data/{label}_phi_spectrum.tsv` -- DFT of phi_1(0,t) for each run
- `data/{label}_P_spectrum.tsv` -- DFT of P(0,t) for each run
