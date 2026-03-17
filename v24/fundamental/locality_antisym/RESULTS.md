# Locality Option 2: Antisymmetric Mode as Causal Gravity -- RESULTS

## Summary

**NEGATIVE**: The antisymmetric mode propagates causally and creates forces between
oscillons, but the mode is UNSTABLE (amplified ~350x from initial perturbation),
destroying the oscillon lattice on timescales ~1000. This rules out the antisymmetric
mode as a gravity mediator -- it is a runaway instability, not a gentle long-range force.

## Phase 1: Single Oscillon Equilibration

- Confining potential V = -sigma*sqrt(P^2+eps^2) + (kc/2)*P^4
- sigma=1.0, kc=0.01, lambda=0.5, m=1.0
- Oscillon mass: 10.605 (code units)
- Breathing frequency: omega = 1.158 (above mass gap m=1.0)
- Equilibrated at t=5000, 463 breathing maxima observed

## Phase 2+3: 8-Oscillon Chain

### A12 Propagation (CAUSAL)

The antisymmetric mode propagates causally through the lattice:

| Oscillon | Distance (lattice) | Distance (phys) | A12 arrival time | Signal speed |
|----------|-------------------|-----------------|------------------|--------------|
| 0        | 0                 | 0               | 0.0              | (source)     |
| 1        | 1                 | 16              | 10.0             | 1.60         |
| 7        | 1                 | 16              | 0.0*             | --           |
| 2        | 2                 | 32              | 29.9             | 1.07         |
| 6        | 2                 | 32              | 29.9             | 1.07         |
| 3        | 3                 | 48              | 49.9             | 0.96         |
| 5        | 3                 | 48              | 44.9             | 1.07         |
| 4        | 4                 | 64              | 54.9             | 1.17         |

*Osc 7 shows t=0 because it is the nearest neighbor and the initial Gaussian
tail may extend slightly.

The A12 group velocity is approximately c=1.0 (within 10-60% depending on distance),
consistent with causal propagation of the massive field.

### A12 Amplification (INSTABILITY)

The perturbation starts at eps=0.01 but max|A12| grows to ~3.5 at ALL oscillons
(amplification factor ~350x). This is the confining potential's instability:
the antisymmetric mode has near-zero gap, so small perturbations grow without bound.

Timeline:
- t=0: A12 = 0.015 at osc 0 only
- t=100: A12 ~ 0.005 at all oscillons (propagation complete)
- t=500: A12 ~ 0.005 (still small, oscillating)
- t=1000: A12 ~ 0.03 (growing)
- t=2500: A12 ~ 1-3 (symmetric structure destroyed)

### Force (Displacement) -- INCONCLUSIVE

The displacement signal is dominated by:
1. Pre-existing oscillon drift (oscillons already wobbling at ~0.9 code units from t=0)
2. The runaway A12 instability destroying the chain

All oscillons show displacement exceeding threshold at the first diagnostic step (t=5),
which is NOT a causal signal but pre-existing lattice dynamics.

Late-time average displacements are all negative (-3.9 to -9.0), indicating a collective
drift rather than a systematic attractive/repulsive pattern.

Energy conservation: dE/E = -12.7% over t=10000 (poor -- chain is unstable).

### Dispersion Relation

The antisymmetric mode spectrum is dominated by the breathing frequency omega=1.26:
- q=0: omega=1.26 (uniform breathing, very strong)
- q=1: omega=1.26 (same frequency, lattice phonon)
- q=2-4: omega~0 (DC offset from instability growth)

No clear gapless phonon dispersion was resolved -- the instability contaminates
the spectrum before the signal can be analyzed.

## Phase 4: Two Isolated Oscillons (D=40)

### Early-time behavior (t < 35)

Before the A12 signal arrives at oscillon 2:
- disp1 = +0.005 (rightward, TOWARD osc 2)
- disp2 = -0.005 (leftward, TOWARD osc 1)
- This is an ATTRACTIVE force at early times

However, this symmetric drift is tiny (0.005 code units in t=35) and likely
reflects the tail overlap between the two oscillons, not the antisymmetric mode.

### A12 arrival

- A12 arrives at oscillon 2 at t=34.9
- Signal speed: D/t = 40/34.9 = 1.15 (consistent with c=1, massive dispersion)

### Post-arrival behavior (t > 35)

After A12 reaches oscillon 2, the oscillons fly apart:
- Separation: 39.4 (initial) -> 204 (final, t=10000)
- Oscillon 2 displacement: +55 (late-time average, AWAY from osc 1)
- Direction: REPULSIVE

The "repulsion" is actually the A12 instability destroying both oscillons.
The antisymmetric mode grows, breaks the phi1=phi2=phi3 symmetry, and the
confined state collapses.

Note: absorbing boundaries also contribute to the outward drift.

## Conclusions

1. **Causal propagation: YES**. The antisymmetric mode propagates at approximately
   the speed of light through the oscillon medium.

2. **Force mediation: AMBIGUOUS at early times, DESTRUCTIVE at late times**.
   There is a hint of attractive force before A12 arrives (tail overlap),
   but the dominant effect is the A12 instability destroying the oscillons.

3. **Gravity candidate: NO**. The confining potential's antisymmetric mode is
   UNSTABLE -- it amplifies ~350x from a 1% perturbation. A gravitational
   mediator must produce a gentle, attractive, 1/r force. Instead we get a
   runaway instability that destroys the soliton structure.

4. **Root cause**: The near-zero gap of the antisymmetric mode (which seemed
   attractive for gravity -- massless mediator!) is actually the symptom of
   an instability. The confining potential creates a saddle point in the
   phi1=phi2=phi3 direction, and the antisymmetric mode rolls off the saddle.

## Parameters Used

```
sigma_c=1.0  kc=0.01  lambda=0.5  mass=1.0  eps_pert=0.01
N_osc=8  d_space=16  pert_osc=0
D_iso=40  t_equil=5000  t_settle=5000  t_evolve=10000
Nx_chain=2560  Nx_iso=8000
```

## Data Files

- `data/locality_profile.tsv` -- equilibrated single oscillon profile
- `data/chain_asym_ts.tsv` -- antisymmetric content vs time at each oscillon
- `data/chain_displacement_ts.tsv` -- center displacement vs time
- `data/chain_energy_ts.tsv` -- total energy vs time
- `data/chain_arrival.tsv` -- arrival times and force summary
- `data/chain_asym_spectrum.tsv` -- antisymmetric phonon dispersion
- `data/isolated_ts.tsv` -- two-oscillon time series
