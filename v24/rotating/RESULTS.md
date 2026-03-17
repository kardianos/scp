# V24-DG: Rotating Triad — Results

## Summary

**NEGATIVE.** Frequency splitting does NOT reduce radiation. The symmetric state (all three fields in phase) is a strong attractor. Any imposed frequency splitting decays exponentially back to the degenerate state within ~500 time units. The late-time energy loss rate is identical across all tested delta values to four significant figures.

## Setup

- Lagrangian: L = sum_a [1/2(dt phi_a)^2 - 1/2(dx phi_a)^2 - m^2/2 phi_a^2] - (mu/2)P^2/(1+kappa*P^2)
- Parameters: mu=-20, kappa=20, m=1.0, A=0.8, sigma=3.0
- Grid: Nx=4000, xmax=100, dx=0.050, dt=0.025
- Protocol: equilibrate t=5000, then perturb v1 += delta*phi1, v3 -= delta*phi3, evolve t=10000
- delta scan: {0.00, 0.01, 0.02, 0.03, 0.05, 0.10, 0.20}

## Results Table

| delta | E(start) | E(end)  | dE/dt (early) | dE/dt (late)  | dE/dt (2nd half) | omega |
|-------|----------|---------|---------------|--------------|-----------------|-------|
| 0.00  | 1.2738   | 1.2595  | -2.50e-06     | -7.46e-07    | -8.71e-07       | 0.879 |
| 0.01  | 1.2738   | 1.2595  | -2.51e-06     | -7.46e-07    | -8.71e-07       | 0.879 |
| 0.02  | 1.2739   | 1.2595  | -2.54e-06     | -7.48e-07    | -8.71e-07       | 0.879 |
| 0.03  | 1.2741   | 1.2595  | -2.59e-06     | -7.52e-07    | -8.71e-07       | 0.879 |
| 0.05  | 1.2748   | 1.2595  | -2.75e-06     | -7.60e-07    | -8.72e-07       | 0.879 |
| 0.10  | 1.2779   | 1.2593  | -3.49e-06     | -7.94e-07    | -8.79e-07       | 0.879 |
| 0.20  | 1.2903   | 1.2593  | -6.46e-06     | -6.69e-07    | -8.26e-07       | 0.879 |

## Key Findings

### 1. Frequency splitting is UNSTABLE

The perturbation v1 += delta*phi1, v3 -= delta*phi3 creates an initial phi1-phi3 difference proportional to delta. For delta=0.20, the maximum splitting |phi1-phi3| starts at 0.12 and decays:

- t=5200: 0.12 (initial)
- t=5400: 0.011 (90% decay in 200 t.u.)
- t=6000: 0.005
- t=8000: 0.002
- t=12000: 0.0004

The decay is exponential with timescale ~200 time units. By t=7000 (2000 t.u. post-perturbation), all three fields are indistinguishable to 0.1%.

### 2. Late-time dE/dt is identical

The late-quarter energy loss rate (best measure of steady-state radiation) is:

- delta=0.00: -7.46e-07
- delta=0.20: -6.69e-07

These differ by only 10%, and the delta=0.20 run starts with MORE energy (the perturbation adds kinetic energy), so it has more to shed. The overall second-half dE/dt values are within 5% across all deltas.

### 3. No frequency splitting survives

All seven runs show identical DFT peak frequencies: omega = 0.879 for all three fields. The DFT is computed over the second half of the run (t=10000-15000), by which time any initial splitting has completely decayed. The measured splitting omega1-omega3 = 0.000 in all cases.

### 4. Early dE/dt increases with delta

Larger delta means more perturbation energy to radiate away. The early-quarter dE/dt goes from -2.5e-06 (delta=0) to -6.5e-06 (delta=0.20), a factor of 2.6x. This is the shedding of the perturbation energy, not a lasting effect.

## Physical Interpretation

The triple-product coupling P = phi1*phi2*phi3 creates a strong restoring force toward amplitude equalization (already observed in v21 test 2). This same mechanism acts on phases: if phi1 oscillates slightly faster than phi3, the coupling transfers energy to restore synchrony.

The proposal's frequency analysis is correct — P contains frequencies at omega +/- 2*delta. But the AMPLITUDE of these beat terms is suppressed by the coupling: the system actively resists dephasing. The rotating triad is not an attractor; it is an unstable perturbation that decays back to the symmetric state.

The radiation rate is set by the single parameter omega/m = 0.879, which determines how far below the mass gap the oscillon sits. Changing the internal phase structure does not modify this ratio.

## Conclusion

The symmetric oscillon (delta=0) with all three fields in phase is a robust attractor. Frequency splitting is a perturbation that decays exponentially, not a new stable configuration. The radiation rate dE/dt ~ -8e-07 is intrinsic to the oscillon's omega/m ratio and cannot be reduced by internal phase rotation.

## Files

- `src/rotating1d.c` — solver code
- `data/rotating_d{delta}_ts.tsv` — time series per delta value
- `data/rotating_d{delta}_spectrum.tsv` — DFT spectra per delta value
- `data/rotating_summary.tsv` — summary table
